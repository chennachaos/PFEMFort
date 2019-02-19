! Module for writing output as 
! legacy VTK format
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 25-Oct-2017
! Place : Swansea, UK
!
!
!

      MODULE WriterVTK
    !   IMPLICIT NONE

      CONTAINS

!***********************************************************************************************************************
!> \brief Write the mesh and the solution in VTK legacy format
!> \param  ndim         - dimension of the mesh
!> \param  nElem        - number of elements in the mesh
!> \param  nNode        - number of nodes in the mesh
!> \param  npElem       - nodes per element (assumed to be the same for every element)
!> \param  ndof         - number of dof per node
!> \param  coords       - coordinates of the nodes
!> \param  elemNodeConn - elemenet to node connectivity
!> \param  elem_procid  - the id of the processor to which the element belongs to
!> \param  soln         - solution at all the nodes

!> \author C.Kadapa
!> \date Oct 2017
!> \copyright Three Cliffs Geomechanical Analysis
!***********************************************************************************************************************
      SUBROUTINE writeoutputvtk(
     1     ndim,
     2     nElem,
     3     nNode,
     4     npElem,
     5     ndof, 
     6     coords,
     7     elemNodeConn,
     8     elem_procid,
     9     soln,
     +     fileName)

      INTEGER, INTENT(IN) :: ndim, nElem, nNode, npElem, ndof
      INTEGER, DIMENSION(:,:), INTENT(IN) :: elemNodeConn
      INTEGER, DIMENSION(:), INTENT(IN) :: elem_procid
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: coords
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: soln
      INTEGER :: ee, ii, jj, ind
      INTEGER :: n1, n2, n3
      DOUBLE PRECISION :: fact=0.0
      LOGICAL :: FILEEXISTS

      CHARACTER(LEN=*) :: fileName

      ! write(*,*) "iiiiiiiiiiiiiiii"
      ! INQUIRE(file=fileName, EXIST=FILEEXISTS)
      ! ! write(*,*) " FILEEXISTS = ", FILEEXISTS
      ! IF(FILEEXISTS .NEQV. .TRUE.) THEN
      !   write(*,*) "File ... ", infileNodes, "does not exist"
      !   call EXIT(1)
      ! END IF

      ! Open the file and write to it

      ! OPEN(1, file=fileName, STATUS="NEW", ACTION="WRITE")
      OPEN(1, file=fileName, STATUS="UNKNOWN", ACTION="WRITE")

      ! Directives
      write(1,'(A)') "# vtk DataFile Version 4.0"
      write(1,'(A)') "PoissonTwoD example"

      ! ASCII or Binary
      ! write(1,*) "ASCII"
      write(1,*) "ASCII"

      ! Type of dataset : Structured/Unstructured/Rectilinear/Polydata...
      write(1,'(A)') "DATASET UNSTRUCTURED_GRID"
      
      ! Coordinates of the points (nodes)
      write(1,'(A,I10,A)') "POINTS ", nNode, " float"
      IF(ndim == 2) THEN
        DO ii=1,nNode
          write(1,'(F12.6,F12.6,F12.6)') coords(ii,1), coords(ii,2), fact
        END DO
      ELSE
        DO ii=1,nNode
          write(1,'(F12.6,F12.6,F12.6)') coords(ii,1), coords(ii,2), coords(ii,3)
        END DO
      END IF

      ! Element<->Nodes connectivity
      ! In VTK terminology, Cell<->Points
      ! <number of nodes> <n1> <n2> <n3> ...
      ! Starting index is 0
      ind = nElem*(npElem+1)
      write(1,'(A,I10,I10)') "CELLS ", nElem, ind

      IF(ndim == 2) THEN
        IF(npElem == 3) THEN
          ! Linear Triangular element
          n1 = 5
          DO ee=1,nElem
            write(1,'(I10,I10,I10,I10)') npElem,
     1      elemNodeConn(ee,1)-1,
     2      elemNodeConn(ee,2)-1,
     3      elemNodeConn(ee,3)-1
          END DO
        ELSE IF(npElem == 6) THEN
          ! Quadratic Triangular element
          n1 = 22
          DO ee=1,nElem
            write(1,'(I10,I10,I10,I10,I10,I10,I10)') npElem,
     1      elemNodeConn(ee,1)-1,
     2      elemNodeConn(ee,2)-1,
     3      elemNodeConn(ee,3)-1,
     4      elemNodeConn(ee,4)-1,
     5      elemNodeConn(ee,5)-1,
     6      elemNodeConn(ee,6)-1
          END DO
        ELSE
          ! Quadrilateral element
          n1 = 9
          DO ee=1,nElem
            write(1,'(I10,I10,I10,I10,I10)') npElem,
     1      elemNodeConn(ee,1)-1,
     2      elemNodeConn(ee,2)-1,
     3      elemNodeConn(ee,3)-1,
     4      elemNodeConn(ee,4)-1
          END DO
        END IF
      ELSE ! ndim=3
        IF(npElem == 4) THEN
          ! Tetrahedral element
          n1 = 10
          DO ee=1,nElem
            write(1,'(I10,I10,I10,I10,I10)') npElem,
     1      elemNodeConn(ee,1)-1,
     2      elemNodeConn(ee,2)-1,
     3      elemNodeConn(ee,3)-1,
     4      elemNodeConn(ee,4)-1
          END DO
        ELSE IF(npElem == 6) THEN
          ! Prism/Wedge element
          n1 = 13
          DO ee=1,nElem
            write(1,'(I10,I10,I10,I10,I10)') npElem,
     1      elemNodeConn(ee,1)-1,
     2      elemNodeConn(ee,2)-1,
     3      elemNodeConn(ee,3)-1,
     4      elemNodeConn(ee,4)-1,
     5      elemNodeConn(ee,5)-1,
     6      elemNodeConn(ee,6)-1
          END DO
        ELSE
          ! Hexahedral element
          n1 = 12
        END IF
      END IF

      ! Cell type, as per VTK
      write(1,'(A,I10)') "CELL_TYPES", nElem
      DO ee=1,nElem
        write(1,'(I3)') n1
      END DO

      ! Cell data. Processor ID
      write(1,'(A,I10)') "CELL_DATA", nElem
      write(1,'(A)') "SCALARS procid int 1"
      write(1,'(A)') "LOOKUP_TABLE default"
      DO ee=1,nElem
        write(1,'(I3)') elem_procid(ee)
      END DO

      ! Point data
      write(1,'(A,I10)') "POINT_DATA", nNode
      IF(ndof == 1) THEN
        write(1,'(A)') "SCALARS solution float 1"
        write(1,'(A)') "LOOKUP_TABLE default"
        DO ii=1,nNode
          write(1,'(F12.6)') soln(ii)
        END DO
      ELSE IF(ndof == 2) THEN
        write(1,'(A)') "VECTORS solution float"
        DO ii=1,nNode
          ind = (ii-1)*ndof
          write(1,'(F12.6,F12.6,F12.6)') soln(ind+1),soln(ind+2),0.0
        END DO
      ELSE
        write(1,'(A)') "VECTORS solution float"
        DO ii=1,nNode
          ind = (ii-1)*ndof
          write(1,'(F12.6,F12.6,F12.6)') soln(ind+1),soln(ind+2),soln(ind+3)
        END DO
      END IF

      ! close the file
      CLOSE(1)

      END SUBROUTINE writeoutputvtk

      END MODULE WriterVTK
