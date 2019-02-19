!
! 3-noded Lagrange element for Poisson equation in 2D
!

MODULE LAGRANGEELEMENT_2DPOISSON_TRIA3NODE
IMPLICIT NONE
PUBLIC :: LagrangeElem2DPoissonTria3Node

TYPE LagrangeElem2DPoissonTria3Node
    PRIVATE

    !
    ! member variables
    !
    INTEGER :: elenum             ! element number
    INTEGER :: ndim=2             ! dimension
    INTEGER :: degree=1           ! degree of polynomial basis functions
    INTEGER :: npElem=3           ! nodes per element
    INTEGER :: ndof=1             ! number of DOF per node
    INTEGER :: nsize=3            ! number of DOF per element
    INTEGER :: elmType            ! element type
    INTEGER :: matType            ! material type
    INTEGER :: secType            ! section type
    INTEGER :: subdomId           ! subdomain/processor id for parallel code
    INTEGER :: elmTypeNameNum     ! element type name
    INTEGER :: matId              ! material id
    INTEGER :: finiteInt          ! flag for small/finite strain
    INTEGER :: sss                ! plane stress or plane strain or axisymmetric
    INTEGER :: nivGP=0            ! number of internal variables per Gauss point 
    INTEGER :: nGP=1              ! number of Gauss points

    ! internal variables, for plasicity problems
    DOUBLE PRECISION, DIMENSION(:), allocatable :: intVar1, intVar2

    INTEGER :: nodeNums(3), forAssyVec(3), globalDOFnums(3)

    !
    ! member functions
    !
    CONTAINS

        !PROCEDURE ::  Constructor
        PROCEDURE ::  prepareElemData
        PROCEDURE ::  calcStiffnessAndResidual
        PROCEDURE ::  assembleElementMatrixAndVector

END TYPE LagrangeElem2DPoissonTria3Node

CONTAINS 
    SUBROUTINE Constructor(this)
      class(LagrangeElem2DPoissonTria3Node) :: this
          INTEGER :: degree = 1
          INTEGER :: npElem = 3
          INTEGER :: nlbf   = 3
          INTEGER :: ndof   = 1
          INTEGER :: nsize  = 3
    END SUBROUTINE Constructor

    SUBROUTINE prepareElemData(this)
      class(LagrangeElem2DPoissonTria3Node) :: this
          INTEGER :: degree = 1
          INTEGER :: npElem = 3
          INTEGER :: nlbf   = 3
          INTEGER :: ndof   = 1
          INTEGER :: nsize  = 3
    END SUBROUTINE prepareElemData

    SUBROUTINE calcStiffnessAndResidual(this, vertices, Klocal, Flocal)
      IMPLICIT NONE
      class(LagrangeElem2DPoissonTria3Node) :: this

      DOUBLE PRECISION, DIMENSION(3,3) :: Klocal
      DOUBLE PRECISION, DIMENSION(3) :: Flocal
      DOUBLE PRECISION, DIMENSION(:,:) :: vertices
      DOUBLE PRECISION, DIMENSION(3,2) :: Bmat
      DOUBLE PRECISION, DIMENSION(2,3) :: BmatTrans
      DOUBLE PRECISION :: x1, x2, x3, y1, y2, y3
      DOUBLE PRECISION :: area
      DOUBLE PRECISION :: kx=1.0, ky=1.0

      x1 = vertices(this%nodeNums(1),1)
      y1 = vertices(this%nodeNums(1),2)
      x2 = vertices(this%nodeNums(2),1)
      y2 = vertices(this%nodeNums(2),2)
      x3 = vertices(this%nodeNums(3),1)
      y3 = vertices(this%nodeNums(3),2)

      area = 0.5*(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
      write(*,*) " area = ", area

      Bmat(1,1) = y2-y3; Bmat(2,1) = y3-y1; Bmat(3,1) = y1-y2
      Bmat(1,2) = x3-x2; Bmat(2,2) = x1-x3; Bmat(3,2) = x2-x1
      Bmat = Bmat/(2.0*area)

      Klocal = 0.0
      Flocal = 0.0

      BmatTrans = TRANSPOSE(Bmat)
  
      Klocal = MATMUL(Bmat, BmatTrans)
      Klocal = area*Klocal

      ! Flocal -= Klocal*dispC;

    END SUBROUTINE calcStiffnessAndResidual

    SUBROUTINE assembleElementMatrixAndVector(this, firstIter, flag, Klocal, Flocal, Mat, rhs, reac)
      IMPLICIT NONE
      class(LagrangeElem2DPoissonTria3Node) :: this
      LOGICAL :: firstIter, flag
      DOUBLE PRECISION :: fact=1.0
      DOUBLE PRECISION, dimension(:,:) :: Mat, Klocal
      DOUBLE PRECISION, dimension(:) :: Flocal, rhs, reac
      INTEGER :: i, j, r, c

      do i=1,this%npElem
        !add reaction forces
        !reac(globalDOFnums(ii)) = reac(globalDOFnums(ii)) + Flocal(ii);
        r = this%forAssyVec(i)

        if(r .ne. 0) then
          rhs(r) = rhs(r) + Flocal(i);

          do j=1,this%npElem
            c = this%forAssyVec(j)
            if(c .ne. 0) then
              Mat(r, c) = Mat(r,c) + Klocal(i,j)
            end if
          end do
        end if
      end do

      if(firstIter .eqv. .true.) then
        do i=1, this%npElem
          r = this%forAssyVec(i)

          if(r .ne. 0) then
            !fact = var1applied(globalDOFnums(ii))

            do j=1,this%npElem
              c = this%forAssyVec(j)
              if(c .ne. 0) then
                rhs(r) = rhs(r) - Klocal(j,i) * fact
              end if
            end do
          end if
        end do
      end if
    END SUBROUTINE assembleElementMatrixAndVector

END MODULE LAGRANGEELEMENT_2DPOISSON_TRIA3NODE


