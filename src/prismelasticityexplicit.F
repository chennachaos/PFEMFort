!
! Program for 2D explicit dynamics using linear Quadrilateral elements
! Only linear elasticity is assumed in this implementation
!
! This implementation is for multiple processors
!
! This implementation uses Petsc Vectors
!
! In this implementation, several entries are inserted at a time,
! into PETSc matrices/vectors, as recommended by PETSc.
! This is efficient for large-scale models.
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 10-Jan-2018
! Place : Swansea, UK
!
!
!

      PROGRAM QuadElasticityExplicit

      USE ElementUtilitiesElasticity3D
      USE WriterVTK
      ! USE Module_SolverPetsc

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

      INTEGER :: ndim=3
      INTEGER :: ndof=3
      INTEGER :: npElem=6 ! number of nodes per element
      INTEGER :: io, nn
      DOUBLE PRECISION :: PI=ACOS(-1.0)
      DOUBLE PRECISION :: xc, yc, area, fact, val
      DOUBLE PRECISION :: xNode(6), yNode(6), zNode(6)
      DOUBLE PRECISION :: Klocal(18,18), Flocal(18)
      DOUBLE PRECISION :: dispElem(18), veloElem(18), acceElem(18)
      DOUBLE PRECISION :: elemData(50), timeData(50)


      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coords
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DirichletBCs
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ForceBCs
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

      INTEGER, DIMENSION(:), ALLOCATABLE :: elem_proc_id
      INTEGER, DIMENSION(:), ALLOCATABLE :: node_proc_id
      INTEGER, DIMENSION(:), ALLOCATABLE :: locally_owned_nodes
      INTEGER, DIMENSION(:), ALLOCATABLE :: node_map_get_old
      INTEGER, DIMENSION(:), ALLOCATABLE :: node_map_get_new
      INTEGER, DIMENSION(:), ALLOCATABLE :: nNode_local_vector

      INTEGER, DIMENSION(:), ALLOCATABLE :: eptr, eind
      INTEGER, DIMENSION(:), ALLOCATABLE :: displs

      ! for METIS partitioning
      INTEGER, pointer     :: vwgt=>null(), vsize=>null()
      DOUBLE PRECISION, pointer     :: tpwgts=>null()
      INTEGER :: nParts, options_metis(100)

      INTEGER :: forAssyVec(18), elemDofGlobal(18)

      CHARACTER(len=32) :: arg

      INTEGER :: stepsCompleted, stepsMax, fileCount
      DOUBLE PRECISION :: timeNow, timeFinal, dt
      DOUBLE PRECISION :: af, am, rhoInf, beta, gamm


      DOUBLE PRECISION :: DTT, IDTT
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  disp, dispPrev, dispPrev2, dispCur
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  velo, veloPrev, veloCur
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  acce, accePrev, acceCur
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  globalM, rhsVec

      ! Vec solnTemp
      ! Vec globalM, rhsVec


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
      PetscInt  n1, n2, n3, n4, n5, n6, nsize
      PetscInt  nodes_per_side, objval

      Vec            vec_SEQ;
      VecScatter     ctx;

      PetscScalar xx_v(1)
      PetscOffset xx_i

      PetscInt  n_mpi_procs     ! total number of processors in the group
      PetscInt  this_mpi_proc   ! rank of the current processor

      CHARACTER (LEN=100) :: infileNodes
      CHARACTER (LEN=100) :: infileElems
      CHARACTER (LEN=100) :: infileDBCs
      CHARACTER (LEN=100) :: infileFBCs
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
      ! call PetscInitialize("petsc_options.dat", errpetsc)
      CHKERRQ(errpetsc)

      ! tstart = MPI_Wtime()

      call MPI_Comm_size(PETSC_COMM_WORLD, n_mpi_procs, errpetsc);
      call MPI_Comm_rank(PETSC_COMM_WORLD, this_mpi_proc, errpetsc);

      WRITE(charTemp,*) " this_mpi_proc = ", this_mpi_proc, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)


      CALL getarg(1, arg)
      infileNodes = arg
      READ(arg,*) charTemp
      charTemp = trim(charTemp) // "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)


      CALL getarg(2, arg)
      infileElems = arg
      READ(arg,*) charTemp
      charTemp = trim(charTemp) // "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)


      CALL getarg(3, arg)
      infileDBCs = arg
      READ(arg,*) charTemp
      charTemp = trim(charTemp) // "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

      IF( iargc() == 4 ) THEN
        CALL getarg(4, arg)
        infileFBCs = arg
        READ(arg,*) charTemp
        charTemp = trim(charTemp) // "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      END IF


      call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)
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
      nNode_global = 0
      OPEN(1, file=infileNodes,STATUS="OLD",ACTION="READ")
      DO
        READ(1,*, iostat=io)
        IF (io/=0) EXIT
        nNode_global = nNode_global + 1
      END DO 
      WRITE(charTemp,*) "Number of nodes = ", nNode_global, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      CLOSE(1)

      ALLOCATE(coords(nNode_global,ndim))

      ! Open the file and store nodal coordinates
       
      OPEN(1, file=infileNodes,STATUS="OLD",ACTION="READ")

      DO ii=1,nNode_global
        READ(1,*, iostat=io) nn, coords(ii,1), coords(ii,2), coords(ii,3)
        ! WRITE(*,*) nn, coords(ii,1), coords(ii,2)
      END DO 
      CLOSE(1)

      call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)
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
      nElem_global = 0
      OPEN(2, file=infileElems,STATUS="OLD",ACTION="READ")
      DO
        READ(2,*, iostat=io)
        IF (io/=0) EXIT
        nElem_global = nElem_global + 1
      END DO 
      WRITE(charTemp,*) "Number of Elements = ", nElem_global, "\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
      CLOSE(2)

      call PetscPrintf(PETSC_COMM_WORLD, "reading elements \n", errpetsc)

      ! This program is hardcoded for triangular elements
      ALLOCATE(elemNodeConn(nElem_global,npElem))

      ! Open the file and store nodal coordinates

      OPEN(2, file=infileElems,STATUS="OLD",ACTION="READ")
      DO ii=1,nElem_global
        READ(2,*, iostat=io) nn, n1, n2, n3, n4, n5, n6
        elemNodeConn(ii,1) = n1
        elemNodeConn(ii,2) = n2
        elemNodeConn(ii,3) = n3
        elemNodeConn(ii,4) = n4
        elemNodeConn(ii,5) = n5
        elemNodeConn(ii,6) = n6
        ! write(*,*) nn, n1, n2, n3
      END DO 
      CLOSE(2)
      call PetscPrintf(PETSC_COMM_WORLD, "reading elements \n", errpetsc)

      call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)
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

        NodeTypeOld(n1, n2) = 1
        ind = (n1-1)*ndof + n2
        solnApplied(ind) = fact
      END DO 
      CLOSE(3)

      call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)
      !
      ! Read Force BC data
      ! Applied nodal forces only

      nFBC = 0
      IF( iargc() == 4 ) THEN
        ! check if the file exists

        INQUIRE(file=infileFBCs, EXIST=FILEEXISTS)
        WRITE(charTemp,*) " FILEEXISTS = ", FILEEXISTS, "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
        IF(FILEEXISTS .NEQV. .TRUE.) THEN
          write(*,*) "File ... ", infileFBCs, "does not exist"
          call EXIT(1)
        END IF

        ! Open the file and count number of nodes first
        nFBC = 0
        OPEN(4, file=infileFBCs,STATUS="OLD",ACTION="READ")
        DO
          READ(4,*, iostat=io)
          IF (io/=0) EXIT
          nFBC = nFBC + 1
        END DO 
        WRITE(charTemp,*) "Number of Force BCs = ", nFBC, "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)
        CLOSE(4)

        ! Open the file and process Dirichlet BCs
       
        ALLOCATE(ForceBCs(nFBC,3))

        OPEN(4, file=infileFBCs,STATUS="OLD",ACTION="READ")
        DO ii=1,nFBC
          READ(4,*, iostat=io) n1, n2, fact

          ForceBCs(ii,1) = n1
          ForceBCs(ii,2) = n2
          ForceBCs(ii,3) = fact
        END DO 
        CLOSE(4)
      END IF

      WRITE(charTemp,*) " Input files have been read successfully \n\n"
      call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

      call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !! Partition the mesh. Here METIS is used.
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ind = nNode_global*ndof
      ALLOCATE( disp(ind) )
      ALLOCATE( velo(ind) )
      ALLOCATE( acce(ind) )
      ALLOCATE( dispPrev(ind) )
      ALLOCATE( dispPrev2(ind) )
      ALLOCATE( veloPrev(ind) )
      ALLOCATE( accePrev(ind) )
      ALLOCATE( dispCur(ind) )
      ALLOCATE( veloCur(ind) )
      ALLOCATE( acceCur(ind) )


      ALLOCATE(elem_proc_id(nElem_global))
      elem_proc_id = 0
      ALLOCATE(node_proc_id(nNode_global))
      node_proc_id = 0

      ALLOCATE( node_map_get_old(nNode_global) )
      ALLOCATE( node_map_get_new(nNode_global) )

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
        WRITE(charTemp,*) " Partitioning the mesh \n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

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

          call PetscPrintf(PETSC_COMM_WORLD, " aaaaaaaaaaaa \n", errpetsc)
          call PetscPrintf(PETSC_COMM_WORLD, " Before Metis \n", errpetsc)

        IF(this_mpi_proc == 0) THEN
          DO ee=1,nElem_global
            eptr(ee) = (ee-1)*npElem

            kk = (ee-1)*npElem;

            forAssyVec = elemNodeConn(ee,:)

            DO ii=1, npElem
              eind(kk+ii) = forAssyVec(ii)-1
            END DO
          END DO

          eptr(nElem_global+1) = nElem_global*npElem
          call PetscPrintf(PETSC_COMM_WORLD, " aaaaaaaaaaaa \n", errpetsc)
          ! METIS options
          nodes_per_side = 3

          call METIS_SetDefaultOptions(options_metis)

          call PetscPrintf(PETSC_COMM_WORLD, " Before Metis \n", errpetsc)

          nParts = n_mpi_procs

          ! METIS partition routine
          call METIS_PartMeshNodal(
     1      nElem_global, nNode_global, eptr, eind, vwgt, vsize, 
     2      nParts, tpwgts, options_metis,
     3      objval, elem_proc_id, node_proc_id);

    !     call METIS_PartMeshDual(
    !  1    nElem_global, nNode_global, eptr, eind, vwgt, vsize, 
    !  2    nodes_per_side, nParts, tpwgts, options_metis,
    !  3    objval, elem_proc_id, node_proc_id)

        END IF

        call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

        call PetscPrintf(PETSC_COMM_WORLD, " After Metis \n", errpetsc)

        call MPI_Bcast(elem_proc_id, nElem_global, MPI_INT, 0, PETSC_COMM_WORLD, errpetsc);

        call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

        call MPI_Bcast(node_proc_id, nNode_global, MPI_INT, 0, PETSC_COMM_WORLD, errpetsc);

        call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

        nElem_local=0
        DO ee=1,nElem_global
          IF( elem_proc_id(ee) == this_mpi_proc ) THEN
            nElem_local = nElem_local + 1
          END IF
        END DO

        WRITE(charTemp,*) " nElem_local = ",
     1    this_mpi_proc, nElem_local, "\n"
        call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

        nNode_local = 0
        DO ii=1, nNode_global
          IF( node_proc_id(ii) == this_mpi_proc ) THEN
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
          IF( node_proc_id(ii) == this_mpi_proc ) THEN
            locally_owned_nodes(kk) = ii
            kk = kk+1
          END IF
        END DO

        ! write(*,*) "locall nodes", locally_owned_nodes

        ALLOCATE(displs(n_mpi_procs))

        displs(1) = 0
        DO ii=1,n_mpi_procs-1
          displs(ii+1) = displs(ii) + nNode_local_vector(ii)
        END DO

        ! create a global list of locally_owned_nodes
        ! which will serve as a mapping from 
        ! NEW node numbers to OLD node numbers
        call MPI_Allgatherv(locally_owned_nodes, nNode_local, MPI_INT, 
     1    node_map_get_old, nNode_local_vector, displs, MPI_INT, 
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

        WRITE(*,*) "row_start = ", this_mpi_proc, row_start, row_end
    !     WRITE(charTemp,*) "row_start = ",
    !  1    this_mpi_proc, row_start, row_end, "\n"
    !     call PetscPrintf(PETSC_COMM_WORLD, charTemp, errpetsc)

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
          ! write(*,*) ee, elemNodeConn(ee,:), ElemDofArray(ee,:)
        ! END IF
      END DO

      ! write(*,*) " ElemDofArray "
      nsize=npElem*ndof
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

      call PetscPrintf(PETSC_COMM_WORLD, 
     1 " Initialising petsc solver \n", errpetsc)

    !   call VecCreate(PETSC_COMM_WORLD, rhsVec, errpetsc)
    !   CHKERRQ(errpetsc)
    !   call VecSetSizes(rhsVec, size_local, size_global, errpetsc)
    !   CHKERRQ(errpetsc)
    !   call VecSetFromOptions(rhsVec, errpetsc)
    !   CHKERRQ(errpetsc)
    !   call VecDuplicate(rhsVec, globalM, errpetsc)
    !   CHKERRQ(errpetsc)

    !   call VecSetOption(rhsVec,
    !  1  VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, errpetsc)

    !   call VecSetOption(globalM,
    !  1  VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, errpetsc)


    !   call VecAssemblyBegin(globalM, errpetsc)
    !   CHKERRQ(errpetsc)
    !   call VecAssemblyEnd(globalM, errpetsc)
    !   CHKERRQ(errpetsc)
    !   call VecZeroEntries(globalM, errpetsc)
    !   CHKERRQ(errpetsc)

      ind = nNode_global*ndof
      ALLOCATE( globalM(ind) )
      ALLOCATE( rhsVec(ind) )

      globalM = 0.0

      ! Young's modulus and Poisson's ratio
      !elemData(1) = 500.0;   elemData(2) = 0.4
      !elemData(1) = 10.0**9;   elemData(2) = 0.0
      elemData(1) = 210.0*100000.0;   elemData(2) = 0.1

      ! Density
      !elemData(3) = 1000.0;
      elemData(3) = 7.6;
      ! Body force in X-, Y- and Z- direction
      elemData(4) = 100.0;   elemData(5) = 0.0; elemData(6) = 0.0

      ! Compute the global mass matrix (M)
      ! The Mass is assumed to be lumped so that
      ! the mass matrix is diagonal
      DO ee=1, nElem_global
        ! write(*,*) "ee = ", ee, elem_proc_id(ee)
        IF(elem_proc_id(ee) == this_mpi_proc) THEN

          forAssyVec = ElemDofArray(ee,:)

            kk=1
            DO ii=1, npElem
              ind = (elemNodeConn(ee,ii)-1)*ndof
              DO jj=1,ndof
                elemDofGlobal(kk) = ind+jj
                kk = kk+1
              END DO
            END DO

          ! WRITE(*,*) forAssyVec

          DO ii=1, npElem
            n1 = node_map_get_old(elemNodeConn(ee,ii))

            xNode(ii) = coords(n1,1)
            yNode(ii) = coords(n1,2)
            zNode(ii) = coords(n1,3)
          END DO

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! compute mass matrix and assemble it
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          call MassMatrixLinearPrism(xNode, yNode, zNode, elemData, Flocal)

          write(*,*) Flocal
          write(*,*) " "
          write(*,*) " "

          ! Assemble the element vector
          ! call VecSetValues(globalM, nsize, forAssyVec, Flocal, ADD_VALUES, errpetsc)

          DO ii=1,nsize
            ! globalM(forAssyVec(ii)) = globalM(forAssyVec(ii)) + Flocal(ii)
            globalM(elemDofGlobal(ii)) = globalM(elemDofGlobal(ii)) + Flocal(ii)
          END DO

        END IF
      END DO

      ! ind = nNode_global*ndof
      ! DO ii=1,ind
      !   write(*,*) ii, globalM(ii)
      ! END DO

      ! call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

      !DO ii=1,size_global
        !WRITE(*,*) ii, assyForSoln(ii)
      !END DO

      call PetscPrintf(PETSC_COMM_WORLD, "Computing the solution \n", errpetsc)

      ! time integration parameters
      timeData(2) = 1.0;   timeData(3) = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      disp = 0.0; dispPrev = 0.0; dispPrev2 = 0.0; dispCur = 0.0
      velo = 0.0; veloPrev = 0.0; veloCur = 0.0
      acce = 0.0; accePrev = 0.0; acceCur = 0.0


      af = 1.0;
      am = 1.0;
      stepsMax = 40000
      stepsCompleted=0
      dt = 0.00002
      timeNow=0.0
      timeFinal = stepsMax*dt

      DTT  = dt*dt
      IDTT = 1.0/DTT

      OPEN(1234, file="solnoutput.dat", STATUS="UNKNOWN", ACTION="WRITE")
      fileCount=1

      ! Time loop
      DO WHILE( (stepsCompleted < stepsMax ) .AND. (timeNow < timeFinal) )

        WRITE(*,*) " stepsCompleted = ", stepsCompleted, "timeNow = ", timeNow

        !elemData(4) = 0.0;
        !IF(timeNow <= 0.1) THEN
          !elemData(4) = 1.0
        !END IF

        ! compute time step
        ! dt = 0.01;
        ! dispCur = af*disp + (1.0-af)*dispPrev
        ! veloCur = af*velo + (1.0-af)*veloPrev
        ! acceCur = am*acce + (1.0-am)*accePrev

        ! Loop over elements and compute the RHS

        rhsVec = 0.0

        LoopElem: DO ee=1, nElem_global
          ! IF(elem_proc_id(ee) == this_mpi_proc) THEN
            ! write(*,*) " ee = ", ee

            dispElem = 0.0; veloElem = 0.0; acceElem = 0.0
            DO ii=1, npElem
              n1 = node_map_get_old(elemNodeConn(ee,ii))

              xNode(ii) = coords(n1,1)
              yNode(ii) = coords(n1,2)
              zNode(ii) = coords(n1,3)

              jj = (ii-1)*ndof
              kk = (n1-1)*ndof

              dispElem(jj+1) = disp(kk+1)
              dispElem(jj+2) = disp(kk+2)
              dispElem(jj+3) = disp(kk+3)

              veloElem(jj+1) = velo(kk+1)
              veloElem(jj+2) = velo(kk+2)
              veloElem(jj+3) = velo(kk+3)

              acceElem(jj+1) = acce(kk+1)
              acceElem(jj+2) = acce(kk+2)
              acceElem(jj+3) = acce(kk+3)
            END DO

            ! Compute the element force vector, including residual force
            call ResidualElasticityLinearPrism(
     1        xNode, yNode, zNode,
     2        elemData, timeData,  
     3        dispElem, veloElem, 
     4        Flocal)

            ! Assemble the element vector

            forAssyVec = ElemDofArray(ee,:)

            kk=1
            DO ii=1, npElem
              ind = (elemNodeConn(ee,ii)-1)*ndof
              DO jj=1,ndof
                elemDofGlobal(kk) = ind+jj
                kk = kk+1
              END DO
            END DO

            ! write(*,*) "assembling vector for element ", ee
            DO ii=1,nsize
              ! rhsVec(forAssyVec(ii)) = rhsVec(forAssyVec(ii)) + Flocal(ii)
              rhsVec(elemDofGlobal(ii)) = rhsVec(elemDofGlobal(ii)) + Flocal(ii)
            END DO

          ! END IF
        END DO LoopElem

        ! call MPI_Barrier(PETSC_COMM_WORLD, errpetsc)

        ! Apply Dirichlet BCs
        ! Dirichlet BCs are assumed to be all zeros. So, nothing to do here.


        call PetscPrintf(PETSC_COMM_WORLD,
     1   "Assembly done. Fianlising the solver. \n", errpetsc)

        ! Add specified nodal force 
        DO ii=1,nFBC
          n1 = ForceBCs(ii,1)
          n2 = ForceBCs(ii,2)

          row  = (n1-1)*ndof+n2
          fact = ForceBCs(ii,3)

          rhsVec(row) = rhsVec(row) + fact
        END DO 

        ! Add contribution from the Mass matrix terms, and
        ! solve for displacements
        ! dispTemp = 2.0*dispPrev - dispPrev2

        DO ii=1,size_global
          jj = assyForSoln(ii)

          rhsVec(jj) = rhsVec(jj) + IDTT*globalM(jj)*(2.0*dispPrev(jj)-dispPrev2(jj))

          disp(jj) = (DTT*rhsVec(jj))/globalM(jj);
          ! write(*,*) ii, jj, disp(jj)
        END DO

        WRITE(*,*) " Computing velocity and acceleration"

        ! compute velocity and acceleration
        velo = (disp - dispPrev2)/(2.0*dt);
        acce = (disp - 2.0*dispPrev + dispPrev2)/DTT;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !write(1234,*) timeNow, disp(305), disp(306), velo(305), velo(306)
        !write(1234,*) timeNow, disp(51), disp(52), velo(51), velo(52)
        !write(1234,*) timeNow, disp(142), disp(143), disp(144), velo(142), velo(143), velo(144)
        write(1234,*) timeNow, disp(3097), disp(3098), disp(3099), velo(3097), velo(3098), velo(3099)

        IF(MOD(stepsCompleted, 100) .EQ. 0) THEN
          call PetscPrintf(PETSC_COMM_WORLD, "Writing VTK file \n", errpetsc)

          WRITE(outFileName,'(A,I5.5,A)') "Elasticity-soln-", fileCount, ".vtk"

          ! write the solution to the VTK file
          call writeoutputvtk(
     1     ndim,
     2     nElem_global,
     3     nNode_global,
     4     npElem,
     5     ndof, 
     6     coords,
     7     elemNodeConn,
     8     elem_proc_id,
     9     disp,
     +     outFileName)

          fileCount = fileCount+1
        ENDIF

        ! store the variables
        dispPrev2 = dispPrev
        dispPrev  = disp
        veloPrev  = velo
        accePrev  = acce

        stepsCompleted = stepsCompleted + 1
        timeNow = timeNow + dt

      END DO ! Time loop

      CLOSE(1234)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! End of computations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call PetscPrintf(PETSC_COMM_WORLD, 
     1 "Deallocating the memory \n", errpetsc)


      IF( ALLOCATED(coords) )        DEALLOCATE(coords)
      IF( ALLOCATED(elemNodeConn) )  DEALLOCATE(elemNodeConn)
      IF( ALLOCATED(DirichletBCs) )  DEALLOCATE(DirichletBCs)
      IF( ALLOCATED(ForceBCs) )      DEALLOCATE(ForceBCs)
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


      IF( ALLOCATED(elem_proc_id) ) THEN
          DEALLOCATE(elem_proc_id)
      END IF

      IF( ALLOCATED(node_proc_id) ) THEN
          DEALLOCATE(node_proc_id)
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

      IF( ALLOCATED(nNode_local_vector) ) THEN
          DEALLOCATE( nNode_local_vector )
      END IF

      IF( ALLOCATED(eptr) ) DEALLOCATE(eptr)
      IF( ALLOCATED(eind) ) DEALLOCATE(eind)

      IF( ALLOCATED(displs) )     DEALLOCATE(displs)

      IF( ALLOCATED( disp ) ) DEALLOCATE(disp)
      IF( ALLOCATED( velo ) ) DEALLOCATE(velo)
      IF( ALLOCATED( acce ) ) DEALLOCATE(acce)
      IF( ALLOCATED( dispPrev ) ) DEALLOCATE(dispPrev)
      IF( ALLOCATED( dispPrev2 ) ) DEALLOCATE(dispPrev2)
      IF( ALLOCATED( veloPrev ) ) DEALLOCATE(veloPrev)
      IF( ALLOCATED( accePrev ) ) DEALLOCATE(accePrev)
      IF( ALLOCATED( dispCur ) ) DEALLOCATE(dispCur)
      IF( ALLOCATED( veloCur ) ) DEALLOCATE(veloCur)
      IF( ALLOCATED( acceCur ) ) DEALLOCATE(acceCur)
      IF( ALLOCATED( globalM ) ) DEALLOCATE(globalM)
      IF( ALLOCATED( rhsVec ) ) DEALLOCATE(rhsVec)

      ! call VecDestroy(globalM, errpetsc)
      ! call VecDestroy(rhsVec, errpetsc)


      ! tend = MPI_Wtime()

      ! write(*,*) "That took ", (tend-tstart), "seconds"

      call PetscPrintf(PETSC_COMM_WORLD,
     1  " Program is successful \n", errpetsc)

      call PetscFinalize(errpetsc)
      CHKERRQ(errpetsc)

      END PROGRAM QuadElasticityExplicit
