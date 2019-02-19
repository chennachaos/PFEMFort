!
! Program for paritioning the mesh
!
! This implementation reads the input mesh and
! partitions it into a specified number of domains.
! The assumption is that the mesh is composed of 
! same type of elements.
!
! Author: Dr. Chennakesava Kadapa
! Date  : 01-Nov-2017
! Place : Swansea, UK
!
!
! Inputs:
! ndim    --- number of dimensions
! eType   --- element type (numbers are as per VTK cell types)
!      5 - 3-noded Triangle
!      9 - 4-noded Quadrilateral
!     10 - 4-noded Tetrahedron
!     12 - 8-noded Hexahedron
! metisType    --- metis algorithm - nodal(=1)/dual
! nParts       --- number of partitions
! infileNodes  --- nodal coordinate file
! infileElems  --- element <-> node connectivity file
!
! Outpus:
! partition-mesh.vtk file which can be viewed using Paraview or Mayavi
!
! Note:
! Check the element<->node connectivity 
! and modify the code accordingly
!


      PROGRAM MeshPartition

      IMPLICIT NONE


      ! declare variables

      DOUBLE PRECISION :: tstart, tend, fact

      INTEGER :: ndim          ! number of dimensions
      INTEGER :: eType         ! element type
      INTEGER :: metisType     ! metis algorithm type - nodal/dual
      INTEGER :: nParts        ! number of partitions 
      INTEGER :: nNode_global  ! number of nodes in the whole mesh
      INTEGER :: nElem_global  ! number of global elements (in the whole model)
      INTEGER :: npElem        ! number of nodes per element
      INTEGER :: io            ! for checking input/output (from files)

      CHARACTER(len=32) :: arg

      CHARACTER (LEN=100) :: infileNodes
      CHARACTER (LEN=100) :: infileElems
      CHARACTER (LEN=100) :: outfileName
      LOGICAL :: FILEEXISTS, ISOPEN

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: coords
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: elemNodeConn

      INTEGER, DIMENSION(:), ALLOCATABLE :: elem_proc_id
      INTEGER, DIMENSION(:), ALLOCATABLE :: node_proc_id
      INTEGER, DIMENSION(:), ALLOCATABLE :: eptr
      INTEGER, DIMENSION(:), ALLOCATABLE :: eind
      INTEGER, DIMENSION(:), ALLOCATABLE :: intVecTemp

      ! for METIS partitioning
      INTEGER, pointer     :: vwgt=>null(), vsize=>null()
      DOUBLE PRECISION, pointer     :: tpwgts=>null()
      INTEGER :: options_metis(100)

      INTEGER :: ee, ii, jj, kk, ind
      INTEGER :: n1, n2, n3, n4, n5, n6, n7, n8
      INTEGER :: ncommon_nodes, objval


      call CPU_TIME(tstart)

      !Set file names
      !The file names are specified as inputs from the command line
      IF( iargc() < 4 ) THEN
        WRITE(*,*) "Number of input entries is not sufficient "
        STOP "Aborting..."
      END IF

      CALL getarg(1, arg)
      READ(arg,*) ndim

      CALL getarg(2, arg)
      READ(arg,*) eType

      CALL getarg(3, arg)
      READ(arg,*) metisType

      CALL getarg(4, arg)
      READ(arg,*) nParts

      CALL getarg(5, arg)
      infileNodes = arg

      CALL getarg(6, arg)
      infileElems = arg

      ! Check element type and set the parameters
      IF(eType == 5) THEN
        npElem = 3
        ncommon_nodes = 2
      ELSE IF(eType == 9) THEN
        npElem = 4
        ncommon_nodes = 2
      ELSE IF(eType == 10) THEN
        npElem = 4
        ncommon_nodes = 3
      ELSE IF(eType == 12) THEN
        npElem = 8
        ncommon_nodes = 4
      ELSE
        WRITE(*,*) " Wrong element type"
        STOP "Aborting..."
      END IF

      ALLOCATE( intVecTemp(npElem) )

      !
      ! Read nodal data files
      !

      ! check if the file exists
      INQUIRE(file=infileNodes, EXIST=FILEEXISTS)
      WRITE(*,*) " FILEEXISTS = ", FILEEXISTS
      IF(FILEEXISTS .NEQV. .TRUE.) THEN
        WRITE(*,*) "File ... ", infileNodes, "does not exist"
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
      CLOSE(1)

      ALLOCATE(coords(nNode_global,ndim))

      ! Open the file and store nodal coordinates
      OPEN(1, file=infileNodes,STATUS="OLD",ACTION="READ")

      IF(ndim == 2) THEN
        DO ii=1,nNode_global
          READ(1,*, iostat=io) ind,
     1      coords(ii,1), coords(ii,2)
        END DO 
      ELSE
        DO ii=1,nNode_global
          READ(1,*, iostat=io) ind, 
     1     coords(ii,1), coords(ii,2), coords(ii,3)
        END DO 
      END IF
      CLOSE(1)


      !
      ! Read element data files
      !

      ! check if the file exists
      INQUIRE(file=infileElems, EXIST=FILEEXISTS)
      WRITE(*,*) " FILEEXISTS = ", FILEEXISTS
      IF(FILEEXISTS .NEQV. .TRUE.) THEN
        WRITE(*,*) "File ... ", infileElems, "does not exist"
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
      CLOSE(2)

      ALLOCATE(elemNodeConn(nElem_global,npElem))

      ! Open the file and store nodal coordinates

      OPEN(2, file=infileElems,STATUS="OLD",ACTION="READ")
      IF(npElem == 3) THEN ! Triangular element
        DO ii=1,nElem_global
          READ(2,*, iostat=io) ind, n1, n2, n3
          elemNodeConn(ii,1) = n1
          elemNodeConn(ii,2) = n2
          elemNodeConn(ii,3) = n3
        END DO 
      ELSE IF(npElem == 4) THEN ! Quadrilateral or Tetra element
        DO ii=1,nElem_global
          READ(2,*, iostat=io) ind, n1, n2, n3, n4
          elemNodeConn(ii,1) = n1
          elemNodeConn(ii,2) = n2
          elemNodeConn(ii,3) = n3
          elemNodeConn(ii,4) = n4
        END DO 
      ELSE IF(npElem == 8) THEN ! Hexahedral element
        DO ii=1,nElem_global
          READ(2,*, iostat=io) ind, n1, n2, n3, n4, n5, n6, n7, n8
          elemNodeConn(ii,1) = n1
          elemNodeConn(ii,2) = n2
          elemNodeConn(ii,3) = n3
          elemNodeConn(ii,4) = n4
          elemNodeConn(ii,5) = n5
          elemNodeConn(ii,6) = n6
          elemNodeConn(ii,7) = n7
          elemNodeConn(ii,8) = n8
        END DO 
      ELSE
        STOP "The input file does not match with specified npElem"
      END IF
      CLOSE(2)

      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) " Input files have been read successfully "
      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) " Mesh statistics ..... "
      WRITE(*,*) " Number of elements  = ", nElem_global
      WRITE(*,*) " Number of nodes     = ", nNode_global
      WRITE(*,*) " Nodes per element   = ", npElem
      WRITE(*,*) " "
      WRITE(*,*) " "

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !! Partition the mesh. Here METIS is used.
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      WRITE(*,*) " Partitioning the mesh "

      ALLOCATE(elem_proc_id(nElem_global))
      elem_proc_id = 0
      ALLOCATE(node_proc_id(nNode_global))
      node_proc_id = 0

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

      DO ee=1,nElem_global
        eptr(ee) = (ee-1)*npElem

        kk = (ee-1)*npElem;

        intVecTemp = elemNodeConn(ee,:)

        DO ii=1, npElem
          eind(kk+ii) = intVecTemp(ii)-1
        END DO
      END DO

      eptr(nElem_global+1) = nElem_global*npElem


      ! Call the METIS's options and change them if necessary
      call METIS_SetDefaultOptions(options_metis)

      WRITE(*,*) " Before Metis "

      ! METIS partition routine
      IF(metisType == 1) THEN
        call METIS_PartMeshNodal(
     1      nElem_global, nNode_global, eptr, eind, vwgt, vsize, 
     2      nParts, tpwgts, options_metis,
     3      objval, elem_proc_id, node_proc_id);
      ELSE
        call METIS_PartMeshDual(
     1    nElem_global, nNode_global, eptr, eind, vwgt, vsize, 
     2    ncommon_nodes, nParts, tpwgts, options_metis,
     3    objval, elem_proc_id, node_proc_id)
      END IF

      WRITE(*,*) " After Metis "

      ! IF(ind == METIS_OK) THEN
      !   WRITE(*,*) " METIS partition routine successful "
      ! ELSE
      !   STOP " METIS partition routine FAILED "
      ! END IF

      DO ii=1,nElem_global
        WRITE(*,*) ii, elem_proc_id(ii)
      ENDDO
      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) " "
      DO ii=1,nNode_global
        WRITE(*,*) ii, node_proc_id(ii)
      ENDDO
      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      WRITE(*,*) "Writing VTK file"

      outfileName="partition-mesh.vtk"

      OPEN(11, file=outfileName, STATUS="UNKNOWN", ACTION="WRITE")

      ! Directives
      WRITE(11,'(A)') "# vtk DataFile Version 4.0"
      WRITE(11,'(A)') "Mesh parition"

      ! ASCII or Binary
      WRITE(11,*) "ASCII"

      ! Type of dataset : Structured/Unstructured/Rectilinear/Polydata...
      WRITE(11,'(A)') "DATASET UNSTRUCTURED_GRID"
      
      ! Coordinates of the points (nodes)
      WRITE(11,'(A,I8,A)') "POINTS ", nNode_global, " float"
      fact=0.0
      IF (ndim == 2) THEN
        DO ii=1,nNode_global
          WRITE(11,'(F12.6,F12.6,F12.6)')
     1      coords(ii,1), coords(ii,2), fact
        END DO
      ELSE
        DO ii=1,nNode_global
          WRITE(11,'(F12.6,F12.6,F12.6)')
     1     coords(ii,1), coords(ii,2), coords(ii,3)
        END DO
      END IF

      ! Element<->Nodes connectivity
      ! In VTK terminology, Cell<->Points
      ! <number of nodes> <n1> <n2> <n3> ...
      ! Starting index is 0
      ind = nElem_global*(npElem+1)
      WRITE(11,'(A,I8,I8)') "CELLS ", nElem_global, ind

      IF(ndim == 2) THEN
        IF(npElem == 3) THEN
          ! Triangular element
          n1 = 5
          DO ee=1,nElem_global
            WRITE(11,'(I10,I10,I10,I10)') npElem,
     1      elemNodeConn(ee,1)-1,
     2      elemNodeConn(ee,2)-1,
     3      elemNodeConn(ee,3)-1
          END DO
        ELSE
          ! Quadrilateral element
          n1 = 9
          DO ee=1,nElem_global
            WRITE(11,'(I10,I10,I10,I10,I10)') npElem,
     1      elemNodeConn(ee,1)-1,
     2      elemNodeConn(ee,2)-1,
     3      elemNodeConn(ee,3)-1,
     4      elemNodeConn(ee,4)-1
          END DO
        END IF
      ELSE
        IF(npElem == 4) THEN
          ! Tetrahedral element
          n1 = 10
          DO ee=1,nElem_global
            WRITE(11,'(I6,I12,I12,I12,I12)') npElem,
     1      elemNodeConn(ee,1)-1,
     2      elemNodeConn(ee,2)-1,
     3      elemNodeConn(ee,3)-1,
     4      elemNodeConn(ee,4)-1
          END DO
        ELSE
          ! Hexahedral element
          n1 = 12
          DO ee=1,nElem_global
            WRITE(11,'(I6,I12,I12,I12,I12,I12,I12,I12,I12)') npElem,
     1      elemNodeConn(ee,1)-1,
     2      elemNodeConn(ee,2)-1,
     3      elemNodeConn(ee,4)-1,
     4      elemNodeConn(ee,3)-1,
     5      elemNodeConn(ee,5)-1,
     6      elemNodeConn(ee,6)-1,
     7      elemNodeConn(ee,8)-1,
     8      elemNodeConn(ee,7)-1
          END DO
        END IF
      END IF

      ! Cell type, as per VTK
      WRITE(11,'(A,I8)') "CELL_TYPES", nElem_global
      DO ee=1,nElem_global
        WRITE(11,'(I3)') n1
      END DO

      ! Cell data. Processor ID
      WRITE(11,'(A,I8)') "CELL_DATA", nElem_global
      WRITE(11,'(A)') "SCALARS procid int 1"
      WRITE(11,'(A)') "LOOKUP_TABLE default"
      DO ee=1,nElem_global
        WRITE(11,'(I3)') elem_proc_id(ee)
      END DO

      ! close the file
      CLOSE(11)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      WRITE(*,*) "Deallocating the memory"

      IF( ALLOCATED(coords) )        DEALLOCATE(coords)
      IF( ALLOCATED(elemNodeConn) )  DEALLOCATE(elemNodeConn)
      IF( ALLOCATED(intVecTemp) )    DEALLOCATE(intVecTemp)
      IF( ALLOCATED(elem_proc_id) )  DEALLOCATE(elem_proc_id)
      IF( ALLOCATED(node_proc_id) )  DEALLOCATE(node_proc_id)
      IF( ALLOCATED(eptr) )          DEALLOCATE(eptr)
      IF( ALLOCATED(eind) )          DEALLOCATE(eind)


      call CPU_TIME(tend)
      WRITE(*,*) "That took ", (tend-tstart), "seconds"

      WRITE(*,*) " "
      WRITE(*,*) " "
      WRITE(*,*) "Program is successful"
      WRITE(*,*) " "
      WRITE(*,*) " "

      END PROGRAM MeshPartition
