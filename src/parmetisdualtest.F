
      PROGRAM ParmetisDualTest

      USE, INTRINSIC :: ISO_C_BINDING

      IMPLICIT NONE

include "mpif.h"

      ! declare variables
      !BIND(C)

      DOUBLE PRECISION :: tstart, tend, fact
      CHARACTER (LEN=100) :: charTemp
      INTEGER :: ndim          ! number of dimensions
      INTEGER :: eType         ! element type
      INTEGER :: metisType     ! metis algorithm type - nodal/dual
      INTEGER :: nParts        ! number of partitions 
      INTEGER :: nNode_global  ! number of nodes in the whole mesh
      INTEGER :: nElem_global  ! number of global elements (in the whole model)
      INTEGER :: npElem        ! number of nodes per element
      INTEGER :: io            ! for checking input/output (from files)

      INTEGER, DIMENSION(:), ALLOCATABLE :: elmdist
      INTEGER, DIMENSION(:), ALLOCATABLE :: eptr
      INTEGER, DIMENSION(:), ALLOCATABLE :: eind
      INTEGER, DIMENSION(:), ALLOCATABLE :: xadj_expected
      INTEGER, DIMENSION(:), ALLOCATABLE :: adjncy_expected
      !INTEGER, DIMENSION(:,:), ALLOCATABLE :: xadj
      !INTEGER, DIMENSION(:,:), ALLOCATABLE :: adjncy
      !INTEGER, pointer     :: xadj(:), adjncy(:)
      INTEGER, target      :: xadj_tar(40), adjncy_tar(40)

      !TYPE(c_ptr)  :: xadj_ptr, adjncy_ptr
      INTEGER(KIND=C_INT), pointer :: xadj_ptr, adjncy_ptr

      ! for METIS partitioning
      INTEGER, pointer     :: vwgt=>null(), vsize=>null()
      DOUBLE PRECISION, pointer     :: tpwgts=>null()
      INTEGER :: options_metis(100), errpar, numflag, ncommonnodes

      INTEGER :: ee, ii, jj, kk, ind
      INTEGER :: n1, n2, n3, n4, n5, n6, n7, n8
      INTEGER :: ncommon_nodes, objval

      INTEGER ::  n_mpi_procs     ! total number of processors in the group
      INTEGER ::  this_mpi_proc   ! rank of the current processor


      call MPI_Init(errpar)
      call MPI_Comm_size(MPI_COMM_WORLD, n_mpi_procs, errpar);
      call MPI_Comm_rank(MPI_COMM_WORLD, this_mpi_proc, errpar);

      WRITE(*,*) " n_mpi_procs = ", n_mpi_procs
      WRITE(*,*) " this_mpi_proc = ", this_mpi_proc

      if(n_mpi_procs .EQ. 1) then
        WRITE(*,*) "This program should be called with at least two processors "
        STOP "Aborting now ..."
      end if


      allocate(elmdist(3))
      elmdist(1) = 0
      elmdist(2) = 4
      elmdist(3) = 8

      allocate(eptr(5))
      eptr(1) = 0
      eptr(2) = 3
      eptr(3) = 6
      eptr(4) = 9
      eptr(5) = 12

      allocate(eind(12))

      if(this_mpi_proc == 0) then
        eind(1)  = 1;      eind(2)  = 2;      eind(3)  = 5;
        eind(4)  = 1;      eind(5)  = 5;      eind(6)  = 4;
        eind(7)  = 2;      eind(8)  = 3;      eind(9)  = 6;
        eind(10) = 2;      eind(11) = 6;      eind(12) = 5;
      else
        eind(1)  = 4;      eind(2)  = 5;      eind(3)  = 8;
        eind(4)  = 4;      eind(5)  = 8;      eind(6)  = 7;
        eind(7)  = 5;      eind(8)  = 6;      eind(9)  = 9;
        eind(10) = 5;      eind(11) = 9;      eind(12) = 8;
      endif

      do ii=1,12
       eind(ii) = eind(ii)-1
      enddo

      WRITE(*,*) " Before Metis "
      call MPI_Barrier(MPI_COMM_WORLD, errpar)

      numflag=0
      ncommonnodes=2
      !allocate(xadj(40))
      !allocate(adjncy(40))

      !xadj=>xadj_tar
      !adjncy=>adjncy_tar

      call ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, numflag, 
     1   ncommonnodes, xadj_ptr, adjncy_ptr, MPI_COMM_WORLD)

      WRITE(*,*) " After Metis "

      ! IF(ind == METIS_OK) THEN
      !   WRITE(*,*) " METIS partition routine successful "
      ! ELSE
      !   STOP " METIS partition routine FAILED "
      ! END IF

      !DO ii=1,9
        WRITE(*,*) xadj_ptr
      !ENDDO
      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) " "
      call MPI_Barrier(MPI_COMM_WORLD, errpar)
      WRITE(*,*) adjncy_ptr
      !DO ii=1,15
        !WRITE(*,*) ii, adjncy_ptr(ii)
      !ENDDO
      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) " "

      ! these are the arrays for the global mesh
      ALLOCATE(xadj_expected(9))
      xadj_expected(1) =  0
      xadj_expected(2) =  2
      xadj_expected(3) =  4
      xadj_expected(4) =  5
      xadj_expected(5) =  7
      xadj_expected(6) = 10
      xadj_expected(7) = 11
      xadj_expected(8) = 13
      xadj_expected(9) = 15

      ALLOCATE(adjncy_expected(15))
      adjncy_expected(1)  =  2
      adjncy_expected(2)  =  4
      adjncy_expected(3)  =  1
      adjncy_expected(4)  =  5
      adjncy_expected(5)  =  4
      adjncy_expected(6)  =  1
      adjncy_expected(7)  =  7
      adjncy_expected(8)  =  2
      adjncy_expected(9)  =  6
      adjncy_expected(10) =  8
      adjncy_expected(11) =  5
      adjncy_expected(12) =  4
      adjncy_expected(13) =  8
      adjncy_expected(14) =  5
      adjncy_expected(15) =  7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      WRITE(*,*) "Deallocating the memory"

      IF( ALLOCATED(elmdist) )       DEALLOCATE(elmdist)
      IF( ALLOCATED(eptr) )          DEALLOCATE(eptr)
      IF( ALLOCATED(eind) )          DEALLOCATE(eind)
      !IF( ALLOCATED(xadj) )          DEALLOCATE(xadj)
      !IF( ALLOCATED(adjncy) )        DEALLOCATE(adjncy)
      IF( ALLOCATED(xadj_expected) )          DEALLOCATE(xadj_expected)
      IF( ALLOCATED(adjncy_expected) )        DEALLOCATE(adjncy_expected)

      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) "Program is successful"
      WRITE(*,*) " "
      WRITE(*,*) " "

      call MPI_Finalize(errpar)

      END PROGRAM ParmetisDualTest
