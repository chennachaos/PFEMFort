! Subroutines for calculating the stiffness
! and residual (or force) for the simple
! Lagrange family of elements for
!
! 2D Linear Elasticity using displacement formulation
!
! Author: Dr. Chennakesava Kadapa
! Date  : 25-Oct-2017
! Place : Swansea, UK
!
!
!

      MODULE ElementUtilitiesElasticity2D

      USE QuadratureUtilities
      USE ElementUtilitiesBasisFuncs

      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE StiffnessResidualElasticityLinearTria(
     1   xNode, yNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   valC,            ! solution t_{n+af} for the 3 nodes
     5   valDotC,         ! d(solution)/dt at t_{n+am} for the 3 nodes
     6   Klocal, Flocal)  ! element stiffness and residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3)   :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:)   :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)   :: valC, valDotC
      DOUBLE PRECISION, DIMENSION(6,6) :: Klocal, Mlocal
      DOUBLE PRECISION, DIMENSION(3,3) :: Dmat
      DOUBLE PRECISION, DIMENSION(6)   :: Flocal

      DOUBLE PRECISION, DIMENSION(3,6) :: Bmat
      DOUBLE PRECISION, DIMENSION(6,3) :: BmatTrans
      DOUBLE PRECISION :: Jac, dvol, thick, b1, b2, b3, b4
      DOUBLE PRECISION :: grad(2,2), strain(3), stress(3)
      DOUBLE PRECISION :: area, timefact, af, E, nu, bforce(2)
      DOUBLE PRECISION :: N(3), dNdx(3), dNdy(3)

      INTEGER, PARAMETER :: nGP=1, nlbf=3
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TJ, TJp1

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(2)

      ! material parameters, 
      E  = elemData(1)
      nu = elemData(2)
      thick = elemData(3)
      bforce(1) = elemData(4)
      bforce(2) = elemData(5)

      b1 = E/(1.0-nu*nu)

      ! Plane stress
      Dmat(1,1) = b1;     Dmat(1,2) = b1*nu;   Dmat(1,3) = 0.0
      Dmat(2,1) = b1*nu;  Dmat(2,2) = b1   ;   Dmat(2,3) = 0.0
      Dmat(3,1) = 0.0;    Dmat(3,2) = 0.0;     Dmat(3,3) = b1*(1.0-nu)

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      ! call getGaussPointsTria(nGP, gpts1, gpts2, gwts)
      gpts1(1) = 1.0/3.0; gpts2(1) = 1.0/3.0;      gwts(1)  = 0.5


      ! loop over Gauss points and compute 
      ! element stiffness and force vectors
      Klocal = 0.0;      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)

        call computeBasisFunctions2D(.FALSE., 1, degree, param,
     1    xNode, yNode,
     2    N, dNdx, dNdy, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tria element in Elasticity" 
        END IF

        dvol = gwts(gp)*(Jac*thick)

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        ! compute the gradient of displacement first
        grad = 0.0
        DO ii=1, nlbf
          jj = (ii-1)*2
          b1 = valC(jj+1)
          b2 = valC(jj+2)

          grad(1,1) = grad(1,1) + b1*dNdx(ii)
          grad(1,2) = grad(1,2) + b1*dNdy(ii)
          grad(2,1) = grad(2,1) + b2*dNdx(ii)
          grad(2,2) = grad(2,2) + b2*dNdy(ii)

        END DO

        ! now compute the strain
        strain(1) = grad(1,1)
        strain(2) = grad(2,2)
        strain(3) = 0.5*(grad(1,2)+grad(2,1))

        stress = MATMUL(Dmat, strain)

        Bmat = 0.0
        Bmat(1,1) = dNdx(1); Bmat(1,3) = dNdx(2); Bmat(1,5) = dNdx(3)
        Bmat(2,2) = dNdy(1); Bmat(2,4) = dNdy(2); Bmat(2,6) = dNdy(3)
        Bmat(3,1) = dNdy(1); Bmat(3,3) = dNdy(2); Bmat(3,5) = dNdy(3)
        Bmat(3,2) = dNdx(1); Bmat(3,4) = dNdx(2); Bmat(3,6) = dNdx(3)

        Bmat = 0.0
        DO ii=1, nlbf
          TI   = (ii-1)*2+1
          TIp1 = TI+1

          Bmat(1,TI) = dNdx(ii);  Bmat(1,TIp1) = 0.0
          Bmat(2,TI) = 0.0;        Bmat(2,TIp1) = dNdy(ii)
          Bmat(3,TI) = dNdy(ii);  Bmat(3,TIp1) = dNdx(ii)
        END DO

        BmatTrans = TRANSPOSE(Bmat)
        Bmat = MATMUL(Dmat, Bmat)
  
        Klocal = MATMUL(BmatTrans, Bmat)
        Klocal = dvol*Klocal

        ! contributions to the RHS vector from the residuals
        DO ii=1, nlbf
          TI   = (ii-1)*2+1
          TIp1 = TI+1

          b4 = dvol*N(ii)

          Flocal(TI)   = Flocal(TI)   + b4*bforce(1)
          Flocal(TIp1) = Flocal(TIp1) + b4*bforce(2)
        END DO
      END DO

      END SUBROUTINE StiffnessResidualElasticityLinearTria

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ResidualElasticityLinearTria(
     1   xNode, yNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   dispC,           ! solution t_{n+af} for the 3 nodes
     5   veloC,           ! d(solution)/dt at t_{n+am} for the 3 nodes
     6   Flocal)          ! element residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3)   :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:)   :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)   :: dispC, veloC
      DOUBLE PRECISION, DIMENSION(3,3) :: Dmat
      DOUBLE PRECISION, DIMENSION(6)   :: Flocal

      DOUBLE PRECISION, DIMENSION(3,6) :: Bmat
      DOUBLE PRECISION, DIMENSION(6,3) :: BmatTrans
      DOUBLE PRECISION :: Jac, dvol, thick, b1, b2, b3, b4
      DOUBLE PRECISION :: grad(2,2), strain(3), stress(3)
      DOUBLE PRECISION :: area, timefact, af, E, nu, dens, bforce(2)
      DOUBLE PRECISION :: N(3), dNdx(3), dNdy(3)

      INTEGER, PARAMETER :: nGP=1, nlbf=3
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TJ, TJp1

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(2)

      ! material parameters, 
      E  = elemData(1)
      nu = elemData(2)
      dens = elemData(3)
      thick = 1.0;
      bforce(1) = elemData(4)
      bforce(2) = elemData(5)

      ! Plane stress
      ! b1 = E/(1.0-nu*nu)
      ! Dmat(1,1) = b1;     Dmat(1,2) = b1*nu;   Dmat(1,3) = 0.0
      ! Dmat(2,1) = b1*nu;  Dmat(2,2) = b1   ;   Dmat(2,3) = 0.0
      ! Dmat(3,1) = 0.0;    Dmat(3,2) = 0.0;     Dmat(3,3) = b1*(1.0-nu)

      ! Plane strain
      b1 = E/((1.0+nu)*(1.0-2.0*nu))
      Dmat(1,1) = b1*(1.0-nu);     Dmat(1,2) = b1*nu;            Dmat(1,3) = 0.0
      Dmat(2,1) = b1*nu;           Dmat(2,2) = b1*(1.0-nu)   ;   Dmat(2,3) = 0.0
      Dmat(3,1) = 0.0;             Dmat(3,2) = 0.0;              Dmat(3,3) = b1*(1.0-2.0*nu)*0.5

      ! WRITE(*,*) Dmat

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      ! call getGaussPointsTria(nGP, gpts1, gpts2, gwts)
      gpts1(1) = 1.0/3.0; gpts2(1) = 1.0/3.0;      gwts(1)  = 0.5


      ! loop over Gauss points and compute 
      ! element stiffness and force vectors
      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)

        call computeBasisFunctions2D(.FALSE., 1, degree, param,
     1    xNode, yNode,
     2    N, dNdx, dNdy, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tria element in Elasticity" 
        END IF

        dvol = gwts(gp)*(Jac*thick)

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        ! compute the gradient of displacement first
        grad = 0.0
        DO ii=1, nlbf
          jj = (ii-1)*2

          b1 = dispC(jj+1)
          b2 = dispC(jj+2)

          grad(1,1) = grad(1,1) + b1*dNdx(ii)
          grad(1,2) = grad(1,2) + b1*dNdy(ii)
          grad(2,1) = grad(2,1) + b2*dNdx(ii)
          grad(2,2) = grad(2,2) + b2*dNdy(ii)
        END DO

        ! now compute the strain
        strain(1) = grad(1,1)
        strain(2) = grad(2,2)
        strain(3) = 0.5*(grad(1,2)+grad(2,1))

        stress = MATMUL(Dmat, strain)

        ! contributions to the RHS vector from the residuals
        DO ii=1, nlbf
          TI   = (ii-1)*2+1
          TIp1 = TI+1

          b1 = dvol*dNdx(ii)
          b2 = dvol*dNdy(ii)
          b4 = (dens*dvol)*N(ii)

          Flocal(TI)   = Flocal(TI)   + b4*bforce(1) - b1*stress(1) - b2*stress(3)
          Flocal(TIp1) = Flocal(TIp1) + b4*bforce(2) - b1*stress(3) - b2*stress(2)
        END DO
      END DO

      END SUBROUTINE ResidualElasticityLinearTria

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Mass is assumed to be assumed.
      ! So, it is stored as a vector of diagonal vector

      SUBROUTINE MassMatrixLinearTria(
     1   xNode, yNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   Mlocal)          ! Mass matrix, stored as diagonal vector

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3)     :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData
      DOUBLE PRECISION, DIMENSION(6)    :: Mlocal
      DOUBLE PRECISION, DIMENSION(6,6) :: Klocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4, dens, fact
      DOUBLE PRECISION :: N(3), dNdx(3), dNdy(3)

      INTEGER, PARAMETER :: nGP=1, nlbf=3
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2, nsize=6

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(2)

      ! material parameters, 
      dens = elemData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      ! call getGaussPointsTria(nGP, gpts1, gpts2, gwts)
      gpts1(1) = 1.0/3.0; gpts2(1) = 1.0/3.0;      gwts(1)  = 0.5

      ! loop over Gauss points and compute 
      ! element mass matrix
      Klocal = 0.0
      Mlocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)

        call computeBasisFunctions2D(.FALSE., 1, degree, param,
     1    xNode, yNode,
     2    N, dNdx, dNdy, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*Jac

        ! WRITE(*,*) " VOLUME = ", Jac, dvol, dens

        DO ii=1, nlbf
          TI   = (ii-1)*2+1
          TIp1 = TI+1

          b4 = (dens*dvol)*N(ii)

          DO jj=1, nlbf
            TJ   = (jj-1)*2+1
            TJp1 = TJ+1

            fact = b4*N(jj)
            Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + fact
            Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + fact
          END DO
        END DO
      END DO

      ! WRITE(*,*) Klocal

      Mlocal = 0.0
      DO ii=1,nsize
        fact = 0.0
        DO jj=1,nsize
          fact = fact + Klocal(ii,jj)
        END DO
        Mlocal(ii) = fact
      END DO

      END SUBROUTINE MassMatrixLinearTria

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ResidualElasticityTria6(
     1   xNode, yNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   dispC,           ! solution t_{n+af} for the 3 nodes
     5   veloC,           ! d(solution)/dt at t_{n+am} for the 3 nodes
     6   Flocal)          ! element residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(6)   :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:)   :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)   :: dispC, veloC
      DOUBLE PRECISION, DIMENSION(3,3) :: Dmat
      DOUBLE PRECISION, DIMENSION(12)   :: Flocal

      DOUBLE PRECISION, DIMENSION(3,6) :: Bmat
      DOUBLE PRECISION, DIMENSION(6,3) :: BmatTrans
      DOUBLE PRECISION :: Jac, dvol, thick, b1, b2, b3, b4
      DOUBLE PRECISION :: grad(2,2), strain(3), stress(3)
      DOUBLE PRECISION :: area, timefact, af, E, nu, dens, bforce(2)
      DOUBLE PRECISION :: N(6), dNdx(6), dNdy(6)

      INTEGER, PARAMETER :: nGP=3, nlbf=6
      INTEGER :: degree=2, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TJ, TJp1

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(2)

      ! material parameters, 
      E  = elemData(1)
      nu = elemData(2)
      dens = elemData(3)
      thick = 1.0;
      bforce(1) = elemData(4)
      bforce(2) = elemData(5)

      ! Plane stress
      ! b1 = E/(1.0-nu*nu)
      ! Dmat(1,1) = b1;     Dmat(1,2) = b1*nu;   Dmat(1,3) = 0.0
      ! Dmat(2,1) = b1*nu;  Dmat(2,2) = b1   ;   Dmat(2,3) = 0.0
      ! Dmat(3,1) = 0.0;    Dmat(3,2) = 0.0;     Dmat(3,3) = b1*(1.0-nu)

      ! Plane strain
      b1 = E/((1.0+nu)*(1.0-2.0*nu))
      Dmat(1,1) = b1*(1.0-nu);     Dmat(1,2) = b1*nu;            Dmat(1,3) = 0.0
      Dmat(2,1) = b1*nu;           Dmat(2,2) = b1*(1.0-nu)   ;   Dmat(2,3) = 0.0
      Dmat(3,1) = 0.0;             Dmat(3,2) = 0.0;              Dmat(3,3) = b1*(1.0-2.0*nu)*0.5

      ! WRITE(*,*) Dmat

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      call getGaussPointsTria(nGP, gpts1, gpts2, gwts)

      ! loop over Gauss points and compute 
      ! element stiffness and force vectors
      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)

        call computeBasisFunctions2D(.FALSE., 1, degree, param,
     1    xNode, yNode,
     2    N, dNdx, dNdy, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tria element in Elasticity" 
        END IF

        dvol = gwts(gp)*(Jac*thick)

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        ! compute the gradient of displacement first
        grad = 0.0
        DO ii=1, nlbf
          jj = (ii-1)*2

          b1 = dispC(jj+1)
          b2 = dispC(jj+2)

          grad(1,1) = grad(1,1) + b1*dNdx(ii)
          grad(1,2) = grad(1,2) + b1*dNdy(ii)
          grad(2,1) = grad(2,1) + b2*dNdx(ii)
          grad(2,2) = grad(2,2) + b2*dNdy(ii)
        END DO

        ! now compute the strain
        strain(1) = grad(1,1)
        strain(2) = grad(2,2)
        strain(3) = 0.5*(grad(1,2)+grad(2,1))

        stress = MATMUL(Dmat, strain)

        ! contributions to the RHS vector from the residuals
        DO ii=1, nlbf
          TI   = (ii-1)*2+1
          TIp1 = TI+1

          b1 = dvol*dNdx(ii)
          b2 = dvol*dNdy(ii)
          b4 = (dens*dvol)*N(ii)

          Flocal(TI)   = Flocal(TI)   + b4*bforce(1) - b1*stress(1) - b2*stress(3)
          Flocal(TIp1) = Flocal(TIp1) + b4*bforce(2) - b1*stress(3) - b2*stress(2)
        END DO
      END DO

      END SUBROUTINE ResidualElasticityTria6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Mass is assumed to be lumped.
      ! So, it is stored as a vector of diagonal vector

      SUBROUTINE MassMatrixTria6(
     1   xNode, yNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   Mlocal)          ! Mass matrix, stored as diagonal vector

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(6)     :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData
      DOUBLE PRECISION, DIMENSION(12)    :: Mlocal
      DOUBLE PRECISION, DIMENSION(12,12) :: Klocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4, dens, fact
      DOUBLE PRECISION :: N(6), dNdx(6), dNdy(6)

      INTEGER, PARAMETER :: nGP=3, nlbf=6
      INTEGER :: degree=2, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2, nsize=12

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(2)

      ! material parameters, 
      dens = elemData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      call getGaussPointsTria(nGP, gpts1, gpts2, gwts)

      ! loop over Gauss points and compute 
      ! element mass matrix
      Klocal = 0.0
      Mlocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)

        call computeBasisFunctions2D(.FALSE., 1, degree, param,
     1    xNode, yNode,
     2    N, dNdx, dNdy, Jac)

        !DO ii=1,6
          !WRITE(*,*) ii, xNode(ii), yNode(ii)
        !END DO

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tria element in Elasticity" 
        END IF

        dvol = gwts(gp)*Jac

        ! WRITE(*,*) " VOLUME = ", Jac, dvol, dens

        DO ii=1, nlbf
          TI   = (ii-1)*2+1
          TIp1 = TI+1

          b4 = (dens*dvol)*N(ii)

          DO jj=1, nlbf
            TJ   = (jj-1)*2+1
            TJp1 = TJ+1

            fact = b4*N(jj)
            Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + fact
            Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + fact
          END DO
        END DO
      END DO

      ! WRITE(*,*) Klocal

      Mlocal = 0.0
      DO ii=1,nsize
        fact = 0.0
        DO jj=1,nsize
          fact = fact + Klocal(ii,jj)
        END DO
        Mlocal(ii) = fact
      END DO

      END SUBROUTINE MassMatrixTria6


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Subroutine using the 4-point Gauss integration scheme

      SUBROUTINE ResidualElasticityLinearQuad(
     1   xNode, yNode,    ! X and Y coordinates of the 4 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   dispC,           ! solution t_{n+af} for the 4 nodes
     5   veloC,           ! d(solution)/dt at t_{n+am} for the 4 nodes
     6   Flocal)          ! element residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(4)   :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:)   :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)   :: dispC, veloC
      DOUBLE PRECISION, DIMENSION(3,3) :: Dmat
      DOUBLE PRECISION, DIMENSION(8)   :: Flocal

      DOUBLE PRECISION, DIMENSION(3,8) :: Bmat
      DOUBLE PRECISION, DIMENSION(8,3) :: BmatTrans
      DOUBLE PRECISION :: Jac, dvol, thick, b1, b2, b3, b4
      DOUBLE PRECISION :: grad(2,2), strain(3), stress(3)
      DOUBLE PRECISION :: area, timefact, af, E, nu, dens, bforce(2)
      DOUBLE PRECISION :: N(4), dNdx(4), dNdy(4)

      INTEGER, PARAMETER :: nGP=4, nlbf=4
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TJ, TJp1

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(2)

      ! material parameters, 
      E  = elemData(1)
      nu = elemData(2)
      dens = elemData(3)
      thick = 1.0;
      bforce(1) = elemData(4)
      bforce(2) = elemData(5)

      ! Plane stress
      ! b1 = E/(1.0-nu*nu)
      ! Dmat(1,1) = b1;     Dmat(1,2) = b1*nu;   Dmat(1,3) = 0.0
      ! Dmat(2,1) = b1*nu;  Dmat(2,2) = b1   ;   Dmat(2,3) = 0.0
      ! Dmat(3,1) = 0.0;    Dmat(3,2) = 0.0;     Dmat(3,3) = b1*(1.0-nu)

      ! Plane strain
      b1 = E/((1.0+nu)*(1.0-2.0*nu))
      Dmat(1,1) = b1*(1.0-nu);     Dmat(1,2) = b1*nu;            Dmat(1,3) = 0.0
      Dmat(2,1) = b1*nu;           Dmat(2,2) = b1*(1.0-nu)   ;   Dmat(2,3) = 0.0
      Dmat(3,1) = 0.0;             Dmat(3,2) = 0.0;              Dmat(3,3) = b1*(1.0-2.0*nu)*0.5

      ! WRITE(*,*) Dmat

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      call getGaussPointsQuad(nGP, gpts1, gpts2, gwts)

      ! loop over Gauss points and compute element force vectors
      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)

        call computeBasisFunctions2D(.FALSE., 2, degree, param,
     1    xNode, yNode,
     2    N, dNdx, dNdy, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tria element in Elasticity" 
        END IF

        dvol = gwts(gp)*(Jac*thick)

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        ! compute the gradient of displacement first
        grad = 0.0
        DO ii=1, nlbf
          jj = (ii-1)*2

          b1 = dispC(jj+1)
          b2 = dispC(jj+2)

          grad(1,1) = grad(1,1) + b1*dNdx(ii)
          grad(1,2) = grad(1,2) + b1*dNdy(ii)
          grad(2,1) = grad(2,1) + b2*dNdx(ii)
          grad(2,2) = grad(2,2) + b2*dNdy(ii)
        END DO

        ! now compute the strain
        strain(1) = grad(1,1)
        strain(2) = grad(2,2)
        strain(3) = 0.5*(grad(1,2)+grad(2,1))

        stress = MATMUL(Dmat, strain)

        ! contributions to the RHS vector from the residuals
        DO ii=1, nlbf
          TI   = (ii-1)*2+1
          TIp1 = TI+1

          b1 = dvol*dNdx(ii)
          b2 = dvol*dNdy(ii)
          b4 = (dens*dvol)*N(ii)

          Flocal(TI)   = Flocal(TI)   + b4*bforce(1) - b1*stress(1) - b2*stress(3)
          Flocal(TIp1) = Flocal(TIp1) + b4*bforce(2) - b1*stress(3) - b2*stress(2)
        END DO
      END DO

      END SUBROUTINE ResidualElasticityLinearQuad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Subroutine using the 1-point Gauss integration scheme
      ! and with Hourglass stabilistion

      SUBROUTINE ResidualElasticityLinearQuadHG(
     1   xNode, yNode,    ! X and Y coordinates of the 4 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   dispC,           ! solution t_{n+af} for the 4 nodes
     5   veloC,           ! d(solution)/dt at t_{n+am} for the 4 nodes
     6   Flocal)          ! element residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(4)   :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:)   :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)   :: dispC, veloC
      DOUBLE PRECISION, DIMENSION(3,3) :: Dmat
      DOUBLE PRECISION, DIMENSION(8)   :: Flocal
      DOUBLE PRECISION, DIMENSION(4)   :: hgBase, hgVec

      DOUBLE PRECISION, DIMENSION(3,8) :: Bmat
      DOUBLE PRECISION, DIMENSION(8,3) :: BmatTrans
      DOUBLE PRECISION :: Jac, dvol, thick, b1, b2, b3, b4, fact
      DOUBLE PRECISION :: grad(2,2), strain(3), stress(3)
      DOUBLE PRECISION :: area, timefact, af, E, nu, dens, bforce(2)
      DOUBLE PRECISION :: N(4), dNdx(4), dNdy(4)
      DOUBLE PRECISION :: eps, epshat, alpha, ghat(2), qhat(2)
      DOUBLE PRECISION :: HGVISC, HG, HGV, H1, H2, m_WAVES

      INTEGER, PARAMETER :: nGP=1, nlbf=4
      INTEGER :: degree=1, type=2, gp, ii, jj
      INTEGER :: TI, TIp1

      ! material parameters, 
      E  = elemData(1)
      nu = elemData(2)
      dens = elemData(3)
      thick = 0.025
      bforce(1) = elemData(4)
      bforce(2) = elemData(5)

      ! Plane stress
      b1 = E/(1.0-nu*nu)
      Dmat(1,1) = b1;     Dmat(1,2) = b1*nu;   Dmat(1,3) = 0.0
      Dmat(2,1) = b1*nu;  Dmat(2,2) = b1   ;   Dmat(2,3) = 0.0
      Dmat(3,1) = 0.0;    Dmat(3,2) = 0.0;     Dmat(3,3) = b1*(1.0-nu)

      ! Plane strain
      !b1 = E/((1.0+nu)*(1.0-2.0*nu))
      !Dmat(1,1) = b1*(1.0-nu);     Dmat(1,2) = b1*nu;            Dmat(1,3) = 0.0
      !Dmat(2,1) = b1*nu;           Dmat(2,2) = b1*(1.0-nu)   ;   Dmat(2,3) = 0.0
      !Dmat(3,1) = 0.0;             Dmat(3,2) = 0.0;              Dmat(3,3) = b1*(1.0-2.0*nu)*0.5

      ! WRITE(*,*) Dmat

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! compute B matrix
      area = 0.5*((xNode(3)-xNode(1))*(yNode(4)-yNode(2)) + (xNode(2)-xNode(4))*(yNode(3)-yNode(1)))

      ! contribution to the force/residual from the standard part
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      fact = 2.0*area

      N(1) = 0.25;
      N(2) = 0.25;
      N(3) = 0.25;
      N(4) = 0.25;

      dNdx(1) = (yNode(2)-yNode(4))/fact
      dNdx(2) = (yNode(3)-yNode(1))/fact
      dNdx(3) = -dNdx(1)
      dNdx(4) = -dNdx(2)

      dNdy(1) = (xNode(4)-xNode(2))/fact
      dNdy(2) = (xNode(1)-xNode(3))/fact
      dNdy(3) = -dNdy(1)
      dNdy(4) = -dNdy(2)

      Bmat(1,1)  = dNdx(1)
      Bmat(1,2)  = dNdx(2)
      Bmat(1,3)  = dNdx(3)
      Bmat(1,4)  = dNdx(4)

      Bmat(2,1)  = dNdy(1)
      Bmat(2,2)  = dNdy(2)
      Bmat(2,3)  = dNdy(3)
      Bmat(2,4)  = dNdy(4)

      ! compute the gradient of displacement
      grad = 0.0
      DO ii=1, nlbf
        jj = (ii-1)*2

        b1 = dispC(jj+1)
        b2 = dispC(jj+2)

        grad(1,1) = grad(1,1) + b1*dNdx(ii)
        grad(1,2) = grad(1,2) + b1*dNdy(ii)
        grad(2,1) = grad(2,1) + b2*dNdx(ii)
        grad(2,2) = grad(2,2) + b2*dNdy(ii)
      END DO

      ! compute the strain
      strain(1) = grad(1,1)
      strain(2) = grad(2,2)
      strain(3) = 0.5*(grad(1,2)+grad(2,1))

      ! compute the stress
      stress = MATMUL(Dmat, strain)

      ! compute the RHS vector

      dvol = area*thick

      Flocal = 0.0
      DO ii=1, nlbf
        TI   = (ii-1)*2+1
        TIp1 = TI+1

        b1 = dvol*dNdx(ii)
        b2 = dvol*dNdy(ii)
        b4 = (dens*dvol)*N(ii)

        Flocal(TI)   = Flocal(TI)   + b4*bforce(1) - b1*stress(1) - b2*stress(3)
        Flocal(TIp1) = Flocal(TIp1) + b4*bforce(2) - b1*stress(3) - b2*stress(2)
      END DO

      ! contribution to the force/residual from the hourglass stabilisation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      hgBase(1) = 1.0;  hgBase(2) = -1.0;  hgBase(3) = 1.0;  hgBase(4) = -1.0;

!      b1 = xNode(1)-xNode(2)+xNode(3)-xNode(4)
!      b2 = yNode(1)-yNode(2)+yNode(3)-yNode(4)

!      ghat=0.0
!      fact=0.0
!      DO ii=1,4
!        hgVec(ii) = (hgBase(ii) - b1*dNdx(ii) - b2*dNdy(ii))/area

!        jj = (ii-1)*2

!        ghat(1) = ghat(1) + hgVec(ii)*dispC(jj+1)
!        ghat(2) = ghat(2) + hgVec(ii)*dispC(jj+2)

!        fact = fact + dNdx(ii)*dNdx(ii)+ dNdy(ii)*dNdy(ii)
!      ENDDO

!      ! add the contribution to the element vector
!      eps = 1.0
!      epshat = (area*area)*(fact*eps)/12.0

!      qhat = -epshat*ghat

!      !write(*,*) area, fact, ghat(1), ghat(2), qhat(1), qhat(2)

!      DO ii=1,4
!        TI   = (ii-1)*2+1
!        TIp1 = TI+1

!        Flocal(TI)    =  Flocal(TI)   + area*qhat(1)
!        Flocal(TIp1)  =  Flocal(TIp1) + area*qhat(2)
!      ENDDO


      ! Fhg = -0.25*Qhg*Rho*sqrt(A)*c
      HGVISC = 1.0      
      HGV=0.25*HGVISC
      m_WAVES=SQRT(E/dens)
      m_WAVES=50.0

      ! Compute hourglass factor
      HG=HGV*dens*m_WAVES*SQRT(area)

      ! Compute hourglass force Fi = HG*hi  where hi = ViJ*Gamma(J) where Gamma = [+1, -1, +1, -1 ]T
      H1=-(veloC(1)-veloC(3)+veloC(5)-veloC(7))*HG
      H2=-(veloC(2)-veloC(4)+veloC(6)-veloC(8))*HG

      !write(*,*) area, m_WAVES, HG, H1, H2

      Flocal(1) = Flocal(1) + H1
      Flocal(2) = Flocal(2) + H2
      Flocal(3) = Flocal(3) - H1
      Flocal(4) = Flocal(4) - H2
      Flocal(5) = Flocal(5) + H1
      Flocal(6) = Flocal(6) + H2
      Flocal(7) = Flocal(7) - H1
      Flocal(8) = Flocal(8) - H2

      END SUBROUTINE ResidualElasticityLinearQuadHG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Mass is assumed to be assumed.
      ! So, it is stored as a vector of diagonal vector

      SUBROUTINE MassMatrixLinearQuad(
     1   xNode, yNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   Mlocal)          ! Mass matrix, stored as diagonal vector

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(4)     :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData
      DOUBLE PRECISION, DIMENSION(8)    :: Mlocal
      DOUBLE PRECISION, DIMENSION(8,8) :: Klocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4, dens, fact
      DOUBLE PRECISION :: N(4), dNdx(4), dNdy(4)

      INTEGER, PARAMETER :: nGP=4, nlbf=4
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2, nsize=8

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(2)

      ! material parameters, 
      dens = elemData(3)

      ! Gauss point coordinates and weights

      CALL getGaussPointsQuad(nGP, gpts1, gpts2, gwts)

      ! loop over Gauss points and compute 
      ! element mass matrix
      Klocal = 0.0
      Mlocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)

        call computeBasisFunctions2D(.FALSE., 2, degree, param,
     1    xNode, yNode,
     2    N, dNdx, dNdy, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*Jac

        ! WRITE(*,*) " VOLUME = ", Jac, dvol, dens

        DO ii=1, nlbf
          TI   = (ii-1)*2+1
          TIp1 = TI+1

          b4 = (dens*dvol)*N(ii)

          DO jj=1, nlbf
            TJ   = (jj-1)*2+1
            TJp1 = TJ+1

            fact = b4*N(jj)
            Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + fact
            Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + fact
          END DO
        END DO
      END DO

      ! WRITE(*,*) Klocal

      Mlocal = 0.0
      DO ii=1,nsize
        fact = 0.0
        DO jj=1,nsize
          fact = fact + Klocal(ii,jj)
        END DO
        Mlocal(ii) = fact
      END DO

      END SUBROUTINE MassMatrixLinearQuad

      END MODULE ElementUtilitiesElasticity2D
