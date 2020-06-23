!
! Program for 2D Poisson equation using linear triangular elements
!
! This implementation is for multiple processors
!
! This implementation use PetscSolver Module with a Class object
!
! In this implementation, several entries are inserted at a time,
! into PETSc matrices/vectors, as recommended by PETSc.
! This is efficient for large-scale models.
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 23-Oct-2017
! Place : Swansea, UK
!
!
!

      PROGRAM TriaMeshPoissonEquation

      ! IMPLICIT NONE

      USE ElementUtilitiesPoisson
      USE WriterVTK
      USE Module_SolverPetsc

#define PETSC_USE_FORTRAN_MODULES

#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include <petsc/finclude/petscmatdef.h>
#include <petsc/finclude/petsckspdef.h>
#include <petsc/finclude/petscpcdef.h>

#if defined(PETSC_USE_FORTRAN_MODULES)
      USE petscvec
      USE petscmat
      USE petscksp
      USE petscpc
#endif
      IMPLICIT NONE

#if !defined(PETSC_USE_FORTRAN_MODULES)
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
!#include "petsc/finclude/petscis.h"
!#include "petsc/finclude/petscis.h90"
#endif

! #include "metis.h"

      ! declare variables

      PetscErrorCode errpetsc;

      DOUBLE PRECISION :: tstart, tend

      LOGICAL :: FILEEXISTS, ISOPEN

      INTEGER :: ndim=2
      INTEGER :: ndof=1
      INTEGER :: npElem=3 ! number of nodes per element
      INTEGER :: io, nn
      DOUBLE PRECISION :: PI=ACOS(-1.0)
      DOUBLE PRECISION :: xc, yc, area, fact, val
      DOUBLE PRECISION :: xNode(3), yNode(3)
      DOUBLE PRECISION :: Klocal(3,3), Flocal(3)
      DOUBLE PRECISION :: valC(3), valDotC(3)
      DOUBLE PRECISION :: elemData(50), timeData(50)


      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coords
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DirichletBCs
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: solnApplied
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: solnVTK

      INTEGER, DIMENSION(:), ALLOCATABLE :: DirichletBC_row_nums
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: elemNodeConn
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NodeDofArrayOld
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NodeDofArrayNew
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ElemDofArray
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: forAssyMat
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NodeTypeOld
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NodeTypeNew
      INTEGER, DIMENSION(:), ALLOCATABLE :: assyForSoln
      INTEGER, DIMENSION(:), ALLOCATABLE :: diag_nnz, offdiag_nnz

      INTEGER, DIMENSION(:), ALLOCATABLE :: elem_procid
      INTEGER, DIMENSION(:), ALLOCATABLE :: node_procid
      INTEGER, DIMENSION(:), ALLOCATABLE :: locally_owned_nodes
      INTEGER, DIMENSION(:), ALLOCATABLE :: node_map_get_old
      INTEGER, DIMENSION(:), ALLOCATABLE :: node_map_get_new
      INTEGER, DIMENSION(:), ALLOCATABLE :: nNode_local_vector

      INTEGER, DIMENSION(:), ALLOCATABLE :: eptr, eind
      INTEGER, DIMENSION(:), ALLOCATABLE :: recvcounts, displs

      ! for METIS partitioning
      INTEGER, pointer     :: vwgt=>null(), vsize=>null(), options_metis=>null()
      DOUBLE PRECISION, pointer     :: tpwgts=>null()

      INTEGER :: forAssyVec(3)

      CHARACTER(len=32) :: arg

      TYPE(PetscSolver) :: solverpetsc

      Vec solnTemp


      PetscInt  nNode_global  ! number of nodes in the whole mesh
      PetscInt  nNode_local   ! number of nodes owned by the local processor
      PetscInt  nElem_global  ! number of global elements (in the whole model)
      PetscInt  nElem_local   ! number of local  elements (owned by the local processor)
      PetscInt  nDBC   ! number of Dirichlet BCs
      PetscInt  nFBC   ! number of force BCs - specified nodal forces

      PetscInt elem_start   ! starting element in the current processor
      PetscInt elem_end     ! end element in the current processor
      PetscInt size_global  ! total number of global DOF in the whole model
      PetscInt size_local   ! total number of local DOF owned by the current processor
      !
      ! Petsc stores matrices continuous row-wise owned by each processor
      PetscInt node_start   ! number of the first nodes owned by the current processor
      PetscInt node_end     ! number of the last nodes owned by the current processor
      PetscInt row_start    ! first row index of the matrix owned by the current processor
      PetscInt row_end      ! last row index of the matrix owned by the current processor


      PetscInt  ee, ii, jj, kk, ind
      PetscInt  count, row, col
      PetscInt  n1, n2, n3, nsize, e1, e2
      PetscInt  nodes_per_side
      PetscInt  nWeights
      PetscInt  objval, numflag
      ! PetscInt, pointer :: xadj, adjncy

      Vec            vec_SEQ;
      VecScatter     ctx;

      PetscScalar xx_v(1)
      ! PetscScalar, pointer :: xx_v(:)

      PetscOffset xx_i

      PetscInt  n_mpi_procs     ! total number of processors in the group
      PetscInt  this_mpi_proc   ! rank of the current processor

      PetscInt  its;


      CHARACTER (LEN=100) :: infileNodes
      CHARACTER (LEN=100) :: infileElems
      CHARACTER (LEN=100) :: infileDBCs
      CHARACTER (LEN=100) :: charTemp
      CHARACTER (LEN=100) :: outFileName

      !Set file names
      !The file names are specified as inputs from the command line
      IF( iargc() < 3 ) THEN
        WRITE(*,*) "Number of input files is not sufficient "
        WRITE(*,*) "You must enter names of THREE files"
        ! WRITE(*,*) "a.) Node file, b.) Element file, and c.) Dirichlet BC file"
        STOP "Aborting..."
      END IF

      ! intialise PETSc environment
      !
      call PetscInitialize(PETSC_NULL_CHARACTER, errpetsc)
      ! CHKERRQ(errpetsc)

      tstart = MPI_Wtime()

      call MPI_Comm_size(PETSC_COMM_WORLD, n_mpi_procs, errpetsc);
      call MPI_Comm_rank(PETSC_COMM_WORLD, this_mpi_proc, errpetsc);

      WRITE(charTemp,*) " this_mpi_proc = ", this_mpi_proc, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)


      CALL getarg(1, arg)
      infileNodes = arg
      !infileNodes = trim(infileNodes) // (this_mpi_proc+1)
      READ(arg,*) charTemp
      charTemp = trim(charTemp) // "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)


      CALL getarg(2, arg)
      infileElems = arg
      !infileElems = trim(infileElems) // (this_mpi_proc+1)
      READ(arg,*) charTemp
      charTemp = trim(charTemp) // "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)


      CALL getarg(3, arg)
      infileDBCs = arg
      READ(arg,*) charTemp
      charTemp = trim(charTemp) // "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)


      !
      ! Read nodal data files
      !

      ! check if the file exists

      INQUIRE(file=infileNodes, EXIST=FILEEXISTS)
      WRITE(charTemp,*) " FILEEXISTS = ", FILEEXISTS, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      IF(FILEEXISTS .NEQV. .TRUE.) THEN
        write(*,*) "File ... ", infileNodes, "does not exist"
        call EXIT(1)
      END IF

      ! Open the file and count number of nodes first
      nNode_local = 0
      OPEN(1, file=infileNodes,STATUS="OLD",ACTION="READ")
      DO
        READ(1,*, iostat=io)
        IF (io/=0) EXIT
        nNode_local = nNode_local + 1
      END DO 
      WRITE(charTemp,*) "Number of nodes = ", nNode_local, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      CLOSE(1)

      ALLOCATE(coords(nNode_local,ndim))

      ! Open the file and store nodal coordinates

      OPEN(1, file=infileNodes,STATUS="OLD",ACTION="READ")

      DO ii=1,nNode_local
        READ(1,*, iostat=io) nn, coords(ii,1), coords(ii,2)
        ! WRITE(*,*) nn, coords(ii,1), coords(ii,2)
      END DO 
      CLOSE(1)


      !
      ! Read element data files
      !

      ! check if the file exists

      INQUIRE(file=infileElems, EXIST=FILEEXISTS)
      WRITE(charTemp,*) " FILEEXISTS = ", FILEEXISTS, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      IF(FILEEXISTS .NEQV. .TRUE.) THEN
        write(*,*) "File ... ", infileElems, "does not exist"
        call EXIT(1)
      END IF

      ! Open the file and count number of nodes first
      nElem_local = 0
      OPEN(2, file=infileElems,STATUS="OLD",ACTION="READ")
      DO
        READ(2,*, iostat=io)
        IF (io/=0) EXIT
        nElem_local = nElem_local + 1
      END DO 
      WRITE(charTemp,*) "Number of Elements = ", nElem_local, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      CLOSE(2)

      ! This program is hardcoded for triangular elements
      ALLOCATE(elemNodeConn(nElem_local,npElem))

      ! Open the file and store nodal coordinates
       
      OPEN(2, file=infileElems,STATUS="OLD",ACTION="READ")
      DO ii=1,nElem_local
        READ(2,*, iostat=io) nn, n1, n2, n3
        elemNodeConn(ii,1) = n1
        elemNodeConn(ii,2) = n2
        elemNodeConn(ii,3) = n3
        ! write(*,*) nn, n1, n2, n3
      !   WRITE(*,*) nn, elemNodeConn(ii,1), elemNodeConn(ii,2), elemNodeConn(ii,3)
      END DO 
      CLOSE(2)


      !
      ! Read Dirichlet BC data
      !

      ! check if the file exists

      INQUIRE(file=infileDBCs, EXIST=FILEEXISTS)
      WRITE(charTemp,*) " FILEEXISTS = ", FILEEXISTS, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      IF(FILEEXISTS .NEQV. .TRUE.) THEN
        write(*,*) "File ... ", infileDBCs, "does not exist"
        call EXIT(1)
      END IF

      ! Open the file and count number of nodes first
      nDBC = 0
      OPEN(3, file=infileDBCs,STATUS="OLD",ACTION="READ")
      DO
        READ(3,*, iostat=io)
        IF (io/=0) EXIT
        nDBC = nDBC + 1
      END DO 
      WRITE(charTemp,*) "Number of Dirichlet BCs = ", nDBC, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      CLOSE(3)

      ! This program is hardcoded for triangular elements

      ALLOCATE(DirichletBCs(nDBC,3))

      size_global=0
      call MPI_Allreduce(size_local, size_global, 1, MPI_INT, MPI_SUM,
     1         PETSC_COMM_WORLD, errpetsc)

      WRITE(charTemp,*) "Number of Dirichlet BCs = ", nDBC, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)



      ! Initialize NodeType and NodeDofArray arrays
      !
      ALLOCATE(NodeTypeOld(nNode_global, ndof))
      ALLOCATE(NodeTypeNew(nNode_global, ndof))
      ALLOCATE(NodeDofArrayOld(nNode_global, ndof))
      ALLOCATE(NodeDofArrayNew(nNode_global, ndof))
      ind = nNode_global*ndof
      ALLOCATE(solnApplied(ind))

      kk=1
      DO ii=1,nNode_global
        DO jj=1,ndof
          NodeTypeOld(ii,jj)   = 0
          ! NodeTypeNew(ii,jj)   = 0
          NodeDofArrayOld(ii,jj) = 0
          NodeDofArrayNew(ii,jj) = 0
          solnApplied(kk) = 0.0
          kk = kk+1
        END DO
      END DO

      ! Open the file and process Dirichlet BCs
       
      OPEN(3, file=infileDBCs,STATUS="OLD",ACTION="READ")
      DO ii=1,nDBC

        READ(3,*, iostat=io) n1, n2, fact

        DirichletBCs(ii,1) = n1
        DirichletBCs(ii,2) = n2
        DirichletBCs(ii,3) = fact
        ! write(*,*) nn, xc, yc

        NodeTypeOld(n1, n2) = 1
        ind = (n1-1)*ndof + n2
        solnApplied(ind) = fact

      END DO 
      CLOSE(3)

      size_global = 1
      DO ii=1,nNode_global
        DO jj=1,ndof
          IF(NodeTypeOld(ii,jj) == 0) THEN
            NodeDofArrayOld(ii,jj) = size_global
            size_global = size_global + 1
          END IF
        END DO
      END DO

      size_global = size_global-1

      WRITE(charTemp,*) " Input files have been read successfully \n\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

      WRITE(charTemp,*) " Mesh statistics .....\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      WRITE(charTemp,*) " nElem_global   = ", nElem_global, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      WRITE(charTemp,*) " nNode_global   = ", nNode_global, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      WRITE(charTemp,*) " npElem         = ", npElem, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      WRITE(charTemp,*) " ndof           = ", ndof, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      WRITE(charTemp,*) " Total DOF      = ", size_global, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !! Partition the mesh. Here METIS is used.
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      WRITE(charTemp,*) " Partitioning the mesh \n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

      ALLOCATE(elem_procid(nElem_global))
      elem_procid = 0
      ALLOCATE(node_procid(nNode_global))
      node_procid = 0

      ALLOCATE( node_map_get_old(nNode_global) )
      ALLOCATE( node_map_get_new(nNode_global) )

      ! ind = nNode_global*ndof
      ! ALLOCATE( dof_map_get_old(ind) )
      ! ALLOCATE( dof_map_get_new(ind) )


      IF(n_mpi_procs == 1) THEN
        nElem_local = nElem_global
        size_local  = size_global

        node_start = 1
        node_end   = nNode_global

        row_start = 1
        row_end   = size_global

        ! WRITE(*,*) " elem_start = ", this_mpi_proc, elem_start, elem_end

        NodeTypeNew     = NodeTypeOld
        NodeDofArrayNew = NodeDofArrayOld

        kk=1
        DO ii=1, nNode_global
          node_map_get_old(ii) = ii
          node_map_get_new(ii) = ii
        END DO

      ELSE !

        ! array of size (nElem_global+1) which stores
        ! number of nodes per each element,
        ! with 0 as the first entry and
        ! nElem_global*npElem as the last entry
        ind = nElem_global+1
        ALLOCATE( eptr(ind) )

        ! array of size 'nElem_global*npElem'
        ! which stores eleme<->node connectivity information
        ind = nElem_global*npElem
        ALLOCATE( eind(ind) )

        ind = 0
        DO ee=1,nElem_global
          eptr(ee) = (ee-1)*npElem

          kk = (ee-1)*npElem;

          forAssyVec = elemNodeConn(ee, :)

          DO ii=1, npElem
            eind(kk+ii) = forAssyVec(ii) - 1
          END DO

          ind = ind+1
        END DO

        eptr(nElem_global+1) = nElem_global*npElem

        ! METIS options
        nodes_per_side = 2
        nWeights  = 1
        numflag=0;

        ! ALLOCATE(options_metis(METIS_NOPTIONS))

        ! call METIS_SetDefaultOptions(options_metis)

        ! !Specifies the partitioning method.
        ! options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;    ! Multilevel recursive bisectioning.
        ! options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;  ! Multilevel k-way partitioning.

        ! !options[METIS_OPTION_NSEPS] = 10;

        ! options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;  ! Edge-cut minimization
        ! !options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; ! Total communication volume minimization

        ! options[METIS_OPTION_NUMBERING] = 0;  ! C-style numbering is assumed that starts from 0.
        ! !options[METIS_OPTION_NUMBERING] = 1; ! Fortran-style numbering is assumed that starts from 1.

        ! METIS partition routine
        !int ret = METIS_PartMeshNodal(&nElem, &nNode, eptr, eind, NULL, NULL, &n_mpi_procs, NULL, options, &objval, elem_procid, node_procid);
        
        call PetscPrintf(PETSC_COMM_WORLD, " Before Metis \n", errpetsc)

        ind = n_mpi_procs
        ! ind = 20
        call METIS_PartMeshDual(
     1    nElem_global, nNode_global, eptr, eind, vwgt, vsize, 
     2    nodes_per_side, ind, tpwgts, options_metis,
     3    objval, elem_procid, node_procid)

        call PetscPrintf(PETSC_COMM_WORLD, " After Metis \n", errpetsc)

        ! IF(ind == METIS_OK) THEN
        !   WRITE(*,*) " METIS partition routine successful "
        ! ELSE
        !   STOP " METIS partition routine FAILED "
        ! END IF

        nElem_local=0
        DO ee=1,nElem_global
          IF( elem_procid(ee) == this_mpi_proc ) THEN
            nElem_local = nElem_local + 1
          END IF
        END DO

        WRITE(charTemp,*) " nElem_local = ",
     1    this_mpi_proc, nElem_local, "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

        nNode_local = 0
        DO ii=1, nNode_global
          IF( node_procid(ii) == this_mpi_proc ) THEN
            nNode_local = nNode_local + 1
          END IF
        END DO

        WRITE(charTemp,*) " nNode_local = ", 
     1   this_mpi_proc, nNode_local, "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

        ! create the vector (of size n_mpi_procs)
        ! consisting of nNode_local from all the processors
        ! in the communication
        ALLOCATE( nNode_local_vector(n_mpi_procs) )

        call MPI_Allgather(nNode_local, 1, MPI_INT, 
     1    nNode_local_vector, 1, MPI_INT, 
     2    PETSC_COMM_WORLD, errpetsc)

        WRITE(charTemp,*) "nNode_local_vector =",
     1    nNode_local_vector, "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

        ! compute the numbers of first and last nodes
        ! in the local processor
        node_start = 1
        node_end   = nNode_local_vector(1)

        Do ii=1, this_mpi_proc
          node_start = node_start + nNode_local_vector(ii)
          node_end   = node_end   + nNode_local_vector(ii+1)
        END DO

        WRITE(charTemp,*) " node_start = ", node_start, node_end, "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

        ! generate the list of locally owned nodes
        ALLOCATE( locally_owned_nodes(nNode_local) )

        kk=1
        DO ii=1, nNode_global
          IF( node_procid(ii) == this_mpi_proc ) THEN
            locally_owned_nodes(kk) = ii
            kk = kk+1
          END IF
        END DO

        ! write(*,*) "locall nodes", locally_owned_nodes

        ALLOCATE(recvcounts(n_mpi_procs))
        ALLOCATE(displs(n_mpi_procs))

        recvcounts = nNode_local_vector
        displs(1) = 0
        DO ii=1,n_mpi_procs
          displs(ii+1) = displs(ii) + nNode_local_vector(ii)
        END DO

        ! create a global list of locally_owned_nodes
        ! which will serve as a mapping from 
        ! NEW node numbers to OLD node numbers
        call MPI_Allgatherv(locally_owned_nodes, nNode_local, MPI_INT, 
     1    node_map_get_old, recvcounts, displs, MPI_INT, 
     2    PETSC_COMM_WORLD, errpetsc)

        ! IF(this_mpi_proc == 1) THEN
        ! WRITE(*,*) " node_map_get_old "
        ! DO ii=1,nNode_global
        !   WRITE(*,*) ii, node_map_get_old(ii)
        ! END DO

        ! WRITE(*,*) " "
        ! WRITE(*,*) " "
        ! WRITE(*,*) " "
        ! END IF


        ! create an array for mapping from
        ! OLD node numbers to NEW node numbers
        ! Also, generate NodeTypeNew array
        ! for computing the local and global DOF size
        ! as well as creating the element-wise array
        ! for element matrix/vector assembly
        DO ii=1,nNode_global
          n1 = node_map_get_old(ii)
          node_map_get_new(n1) = ii

          DO jj=1,ndof
            NodeTypeNew(ii,jj) = NodeTypeOld(n1,jj)
          END DO
        END DO

        ! compute the size (total DOF) of the local problem

        ! IF(this_mpi_proc == 1) THEN
        ! WRITE(*,*) " NodeTypeNew "
        ! DO ii=1,nNode_global
        !   WRITE(*,*) ii, NodeTypeNew(ii,1)
        ! END DO
        ! END IF

        ! WRITE(*,*) " "
        ! WRITE(*,*) " "
        ! WRITE(*,*) " "

        call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

        ! compute NodeDofArrayNew
        ind = 1
        ! DO ii=node_start, node_end
        DO ii=1, nNode_global
          DO jj=1,ndof
            IF(NodeTypeNew(ii,jj) == 0) THEN
              NodeDofArrayNew(ii,jj) = ind
              ind = ind + 1
            END IF
          END DO
        END DO

        ind = ind-1

        IF(ind /= size_global) THEN
          STOP "Something wrong with NodeDofArrayNew"
        END IF

        ! IF(this_mpi_proc == 1) THEN
        ! WRITE(*,*) " NodeTypeNew "
        ! DO ii=1,nNode_global
        !   WRITE(*,*) ii, NodeDofArrayNew(ii,1)
        ! END DO
        ! END IF

        ! compute first and last row indices
        ! of the rows owned by the local processor
        row_start  =  1e9
        row_end    = -1e9
        size_local = 1
        DO ii=node_start, node_end
          DO jj=1,ndof
            IF(NodeTypeNew(ii,jj) == 0) THEN
              ind = NodeDofArrayNew(ii,jj)
              row_start  = MIN(row_start, ind)
              row_end    = MAX(row_end,   ind)
              size_local = size_local + 1
            END IF
          END DO
        END DO

        size_local = size_local - 1

        WRITE(charTemp,*) "size_local = ",
     1    this_mpi_proc, size_local, size_global, "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

        WRITE(charTemp,*) "row_start = ",
     1    this_mpi_proc, row_start, row_end, "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

        ! check if the sum of local problem sizes 
        ! is equal to that of global problem size

        ind=0
        call MPI_Allreduce(size_local, ind, 1, MPI_INT, MPI_SUM,
     1         PETSC_COMM_WORLD, errpetsc)
        
        IF(ind /= size_global) THEN
          STOP "Sum of local problem sizes is not equal to global size"
        END IF

        ! update elem<->node connectivity 
        ! with new node numbers
        DO ee=1, nElem_global
          Do ii=1, npElem
            elemNodeConn(ee, ii) = 
     1      node_map_get_new( elemNodeConn(ee, ii) )
          END DO
        END DO

        ! update Dirichlet BC information
        ! with new node numbers
        Do ii=1, nDBC
          n1 = node_map_get_new( DirichletBCs(ii,1) )
          DirichletBCs(ii,1) = n1

          n2   = DirichletBCs(ii,2)
          fact = DirichletBCs(ii,3)

          ind = (n1-1)*ndof + n2
          solnApplied(ind) = fact
        END DO 

      END IF !!IF(n_mpi_procs == 1) THEN

      call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !prepare the global matrix pattern
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call PetscPrintf(PETSC_COMM_WORLD,
     1 " Preparing matrix pattern \n", errpetsc)

      ! write(*,*) " Creating arrays 1 "

      ! ElemDofArray is used for element matrix/vector assembly
      ! As the starting index in PETSc is ZERO (unlike ONE in Fortran),
      ! ONE is  substracted from each entry of ElemDofArray

      ind = npElem*ndof
      ALLOCATE(ElemDofArray(nElem_global, ind))

      DO ee=1,nElem_global
        DO ii=1,npElem
          n1 = ndof*(ii-1)
          n2 = elemNodeConn(ee,ii)

          DO jj=1,ndof
            ElemDofArray(ee, n1+jj) = NodeDofArrayNew(n2,jj) - 1
          END DO
        END DO
        ! IF(this_mpi_proc == 0 ) THEN
        !   write(*,*) ee, elemNodeConn(ee,:), ElemDofArray(ee,:)
        ! END IF
      END DO

      ! write(*,*) " ElemDofArray "
      ! nsize=npElem*ndof
      ! DO ee=1,nElem
      !   write(*,*)  ElemDofArray(ee,1), ElemDofArray(ee,2), ElemDofArray(ee,3)
      ! END DO

      ! write(*,*) " Creating arrays 2 "
      ALLOCATE(assyForSoln(size_global))
      ! write(*,*) " Creating arrays 2 "

      count = 1
      DO ii=1,nNode_global
        DO jj=1,ndof
          ! write(*,*) ii, jj, NodeDofArray(ii,jj)
          IF(NodeDofArrayNew(ii,jj) /= 0) THEN
            assyForSoln(count) = (ii-1)*ndof + jj;
            count = count + 1
          END IF
        END DO
      END DO

      call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

      call PetscPrintf(PETSC_COMM_WORLD,
     1 " Preparing matrix pattern DONE \n", errpetsc)


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! create PETSc variables
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      WRITE(charTemp,*) " size_local  = ",
     1  this_mpi_proc, size_local, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

      ! write(*,*) " Creating PETSc vectors 1"

      ! Here, the number of non-zeros per each row is set
      ! as the same for every row.
      ! This is not the best method but it avoids computing
      ! and storing temporary arrays for matrix pattern
      ! We might need to revisit this piece of code
      ! in the future

      ALLOCATE(diag_nnz(size_local))
      ALLOCATE(offdiag_nnz(size_local))

      n1 = 50; n2 = 25
      IF(size_local < 50) THEN
        n1 = size_local
        n2 = n1
      END IF

      ! write(*,*) "n1, n2 = ", this_mpi_proc, n1, n2

      DO ii=1,size_local
        diag_nnz(ii) = n1
        offdiag_nnz(ii) = n2
      END DO

      call PetscPrintf(PETSC_COMM_WORLD, 
     1 " Initialising petsc solver \n", errpetsc)

      !Initialize the petsc solver
      call solverpetsc%initialise(size_local, size_global,
     1  diag_nnz, offdiag_nnz)

      call PetscPrintf(PETSC_COMM_WORLD, 
     1 " Initialise the Matrix pattern \n", errpetsc)

      nsize = npElem*ndof
      Klocal = 0.0

      ! WRITE(*,*) Klocal

      ! Initialise the Matrix pattern
      LoopElem: DO ee=1, nElem_global
        ! write(*,*) "ee = ", ee, elem_procid(ee)
        IF(elem_procid(ee) == this_mpi_proc) THEN

          forAssyVec = ElemDofArray(ee,:)

          ! WRITE(*,*) forAssyVec
          call MatSetValues(solverpetsc%mtx, nsize, forAssyVec,
     1                          nsize, forAssyVec,
     2                          Klocal, INSERT_VALUES, errpetsc)
        END IF
      END DO LoopElem

      call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! compute element matrices and vectors and assemble them
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call PetscPrintf(PETSC_COMM_WORLD, 
     1  " Matrix pattern initialised. \n", errpetsc)
      call PetscPrintf(PETSC_COMM_WORLD, 
     1 " Computing element matrices now. \n", errpetsc)

      ! call solverpetsc%printInfo()

      call solverpetsc%setZero()

      call PetscPrintf(PETSC_COMM_WORLD, 
     1 " Generating element matrices and vectors \n", errpetsc)

      elemData(1) = 1.0;   elemData(2) = 1.0
      timeData(2) = 1.0;   timeData(3) = 0.0
      valC = 0.0;  valDotC = 0.0

      DO ee=1, nElem_global
        IF(elem_procid(ee) == this_mpi_proc) THEN
          ! write(*,*) " ee = ", ee

          n1 = node_map_get_old(elemNodeConn(ee,1))
          n2 = node_map_get_old(elemNodeConn(ee,2))
          n3 = node_map_get_old(elemNodeConn(ee,3))

          xNode(1) = coords(n1,1);  yNode(1) = coords(n1,2)
          xNode(2) = coords(n2,1);  yNode(2) = coords(n2,2)
          xNode(3) = coords(n3,1);  yNode(3) = coords(n3,2)

          ! Compute the element stiffness matrix and force vector
          call StiffnessResidualPoissonLinearTria(
     1      xNode, yNode, 
     2      elemData, timeData,  
     3      valC, valDotC, 
     4      Klocal, Flocal)

          ! Assemble the element matrix
          ! write(*,*) "assembling matrix for element ", ee
          forAssyVec = ElemDofArray(ee,:)

          call MatSetValues(solverpetsc%mtx, nsize, forAssyVec,
     1                            nsize, forAssyVec,
     2                            Klocal, ADD_VALUES, errpetsc)

          ! write(*,*) "applying BCs for element ", ee

          ! Assemble the element vector
          ! Also, apply Dirichlet BCs while assembling
          LoopI: DO ii=1,nsize
            row = forAssyVec(ii)
            IF( row == -1) THEN
              fact = solnApplied(elemNodeConn(ee,ii))
              LoopJ: DO jj=1,nsize
                col = forAssyVec(jj)
                IF( col /= -1) THEN
                  Flocal(jj) = Flocal(jj) - Klocal(jj,ii)*fact
                END IF
              END DO LoopJ
            END IF
          END DO LoopI

          ! IF(this_mpi_proc == 0) THEN
          ! write(*,*) "forAssyVec = ", forAssyVec
          ! write(*,*) " "
          ! write(*,*) Flocal
          ! write(*,*) " "
          ! END IF

          ! write(*,*) "assembling vector for element ", ee
          call VecSetValues(solverpetsc%rhsVec, nsize, forAssyVec,
     1                              Flocal, ADD_VALUES, errpetsc)

        END IF
      END DO

      call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

      call PetscPrintf(PETSC_COMM_WORLD,
     1 "Assembly done. Fianlising the solver. \n", errpetsc)

      ! call solverpetsc%printInfo()

      tstart = MPI_Wtime()

      call solverpetsc%factoriseAndSolve()

      tend = MPI_Wtime()

      WRITE(charTemp,*) "That took ", (tend-tstart), "seconds \n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call PetscPrintf(PETSC_COMM_WORLD, "Writing VTK file", errpetsc)

      ind = nNode_global*ndof
      ALLOCATE(solnVTK(ind))

      ! Applied BC values
      DO ii=1,nNode_global
        n1 = node_map_get_old(ii)
        solnVTK(n1) = solnApplied(ii)
      END DO

      call VecScatterCreateToAll(solverpetsc%solnVec, 
     1 ctx, vec_SEQ, errpetsc)
      call VecScatterBegin(ctx, solverpetsc%solnVec,
     1  vec_SEQ, INSERT_VALUES, SCATTER_FORWARD, errpetsc)
      call VecScatterEnd(ctx, solverpetsc%solnVec, 
     1 vec_SEQ, INSERT_VALUES, SCATTER_FORWARD, errpetsc)

      call VecScatterDestroy(ctx, errpetsc)

      ! Add solution for the free DOF
      call VecGetArray(vec_SEQ, xx_v, xx_i, errpetsc)

      OPEN(1, file="temp.dat", STATUS="UNKNOWN", ACTION="WRITE")

      DO ii=1,size_global
        ind = node_map_get_old(assyForSoln(ii))
        fact = xx_v(xx_i+ii)
        ! fact = xx_v(ii)
        ! x1 = coords(ind,1); y1 = coords(ind,2);
        ! val = (cosh(PI*y1)-sinh(PI*y1)/tanh(PI))*sin(PI*x1)
        write(1,*) ii, ind, fact
        solnVTK(ind) = fact
      END DO
      CLOSE(1)

      call VecRestoreArray(vec_SEQ, xx_v, xx_i, errpetsc)

      ! reset elem<->node connectivity 
      DO ee=1, nElem_global
        Do ii=1, npElem
          elemNodeConn(ee, ii) = 
     1      node_map_get_old( elemNodeConn(ee, ii) )
        END DO
      END DO

      WRITE(outFileName,'(A,I5.5,A)') "PoissonTria-soln.vtk"

      ! write the solution to the VTK file
      IF(this_mpi_proc == 0) THEN
        write(*,*) "Writing VTK file"
        call writeoutputvtk(
     1     ndim,
     2     nElem_global,
     3     nNode_global,
     4     npElem,
     5     ndof, 
     6     coords,
     7     elemNodeConn,
     8     elem_procid,
     9     solnVTK,
     +     outFileName)
      END IF

      call PetscPrintf(PETSC_COMM_WORLD, 
     1 "Deallocating the memory \n", errpetsc)


      IF( ALLOCATED(coords) )        DEALLOCATE(coords)
      IF( ALLOCATED(elemNodeConn) )  DEALLOCATE(elemNodeConn)
      IF( ALLOCATED(DirichletBCs) )  DEALLOCATE(DirichletBCs)
      IF( ALLOCATED(NodeTypeOld) )      DEALLOCATE(NodeTypeOld)
      IF( ALLOCATED(NodeTypeNew) )      DEALLOCATE(NodeTypeNew)
      IF( ALLOCATED(NodeDofArrayOld) )  DEALLOCATE(NodeDofArrayOld)
      IF( ALLOCATED(NodeDofArrayNew) )  DEALLOCATE(NodeDofArrayNew)
      IF( ALLOCATED(ElemDofArray) )  DEALLOCATE(ElemDofArray)
      IF( ALLOCATED(assyForSoln) )   DEALLOCATE(assyForSoln)
      IF( ALLOCATED(solnApplied) )   DEALLOCATE(solnApplied)
      IF( ALLOCATED(solnVTK) )       DEALLOCATE(solnVTK)
      IF( ALLOCATED(diag_nnz) )      DEALLOCATE(diag_nnz)
      IF( ALLOCATED(offdiag_nnz) )   DEALLOCATE(offdiag_nnz)


      IF( ALLOCATED(elem_procid) ) THEN
          DEALLOCATE(elem_procid)
      END IF

      IF( ALLOCATED(node_procid) ) THEN
          DEALLOCATE(node_procid)
      END IF

      IF( ALLOCATED(locally_owned_nodes) ) THEN
          DEALLOCATE(locally_owned_nodes)
      END IF

      IF( ALLOCATED(node_map_get_old) ) THEN
          DEALLOCATE( node_map_get_old )
      END IF
      IF( ALLOCATED(node_map_get_new) ) THEN
          DEALLOCATE( node_map_get_new )
      END IF
      ! IF( ALLOCATED(dof_map_get_old) ) THEN
      !     DEALLOCATE( dof_map_get_old )
      ! END IF
      ! IF( ALLOCATED(dof_map_get_new) ) THEN
      !     DEALLOCATE( dof_map_get_new )
      ! END IF

      IF( ALLOCATED(nNode_local_vector) ) THEN
          DEALLOCATE( nNode_local_vector )
      END IF

      IF( ALLOCATED(eptr) ) DEALLOCATE(eptr)
      IF( ALLOCATED(eind) ) DEALLOCATE(eind)

      IF( ALLOCATED(recvcounts) ) DEALLOCATE(recvcounts)
      IF( ALLOCATED(displs) )     DEALLOCATE(displs)


      call solverpetsc%free()

      ! tend = MPI_Wtime()

      ! write(*,*) "That took ", (tend-tstart), "seconds"

      call PetscPrintf(PETSC_COMM_WORLD,
     1  " Program is successful \n", errpetsc)

      call PetscFinalize(errpetsc)
      CHKERRQ(errpetsc)

      END PROGRAM TriaMeshPoissonEquation
