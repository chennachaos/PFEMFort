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

      MODULE ElementUtilitiesElasticity3D

      USE QuadratureUtilities
      USE ElementUtilitiesBasisFuncs

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ResidualElasticityLinearPrism(
     1   xNode, yNode, zNode,    ! X,Y and Z coordinates of the 6 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   dispC,           ! solution t_{n+af} for the 3 nodes
     5   veloC,           ! d(solution)/dt at t_{n+am} for the 6 nodes
     6   Flocal)          ! element residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(6)   :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)   :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)   :: dispC, veloC
      DOUBLE PRECISION, DIMENSION(6,6) :: Dmat
      DOUBLE PRECISION, DIMENSION(18)   :: Flocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4
      DOUBLE PRECISION :: grad(3,3), strain(6), stress(6)
      DOUBLE PRECISION :: area, timefact, af, E, nu, dens, bforce(3)
      DOUBLE PRECISION :: N(6), dNdx(6), dNdy(6), dNdz(6)

      INTEGER, PARAMETER :: nGP=1, nlbf=6
      INTEGER :: DEGREE=1, ETYPE=6, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      E    = elemData(1)
      nu   = elemData(2)
      dens = elemData(3)

      bforce(1) = elemData(4)
      bforce(2) = elemData(5)
      bforce(3) = elemData(6)


      b1 = E/((1.0+nu)*(1.0-2.0*nu))
      b2 = (1.0-2.0*nu)/2.0

      Dmat = 0.0
      Dmat(1,1) = b1*(1.0-nu); Dmat(1,2) = b1*nu;       Dmat(1,3) = b1*nu
      Dmat(2,1) = b1*nu;       Dmat(2,2) = Dmat(1,1);   Dmat(2,3) = b1*nu
      Dmat(3,1) = b1*nu;       Dmat(3,2) = b1*nu;       Dmat(3,3) = Dmat(1,1)
      Dmat(4,4) = b1*b2
      Dmat(5,5) = b1*b2
      Dmat(6,6) = b1*b2


      ! WRITE(*,*) Dmat

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      call getGaussPointsPrism(nGP, gpts1, gpts2, gpts3, gwts)


      ! loop over Gauss points and compute 
      ! element stiffness and force vectors
      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., ETYPE, DEGREE, param,
     1    xNode, yNode, zNode,
     2    N, dNdx, dNdy, dNdz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tria element in Elasticity" 
        END IF

        dvol = gwts(gp)*Jac

        !WRITE(*,*) " VOLUME = ", gp, Jac, dvol

        ! compute the gradient of displacement first
        grad = 0.0
        DO ii=1, nlbf
          jj = (ii-1)*3

          b1 = dispC(jj+1)
          b2 = dispC(jj+2)
          b3 = dispC(jj+3)

          grad(1,1) = grad(1,1) + b1*dNdx(ii)
          grad(1,2) = grad(1,2) + b1*dNdy(ii)
          grad(1,3) = grad(1,3) + b1*dNdz(ii)

          grad(2,1) = grad(2,1) + b2*dNdx(ii)
          grad(2,2) = grad(2,2) + b2*dNdy(ii)
          grad(2,3) = grad(2,3) + b2*dNdz(ii)

          grad(3,1) = grad(3,1) + b3*dNdx(ii)
          grad(3,2) = grad(3,2) + b3*dNdy(ii)
          grad(3,3) = grad(3,3) + b3*dNdz(ii)
        END DO

        ! now compute the strain
        strain(1) = grad(1,1) !sxx
        strain(2) = grad(2,2) !syy
        strain(3) = grad(3,3) !szz
        strain(4) = 0.5*(grad(1,2)+grad(2,1))  !sxy
        strain(5) = 0.5*(grad(2,3)+grad(3,2))  !syz
        strain(6) = 0.5*(grad(1,3)+grad(3,1))  !szx

        stress = MATMUL(Dmat, strain)

        ! contributions to the RHS vector from the residuals
        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b1 = dvol*dNdx(ii)
          b2 = dvol*dNdy(ii)
          b3 = dvol*dNdz(ii)
          b4 = (dens*dvol)*N(ii)

          Flocal(TI)   = Flocal(TI)   + b4*bforce(1) - b1*stress(1) - b2*stress(4) - b3*stress(6)
          Flocal(TIp1) = Flocal(TIp1) + b4*bforce(2) - b1*stress(4) - b2*stress(2) - b3*stress(5)
          Flocal(TIp2) = Flocal(TIp2) + b4*bforce(3) - b1*stress(6) - b2*stress(5) - b3*stress(3)
        END DO
      END DO

      END SUBROUTINE ResidualElasticityLinearPrism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Mass is assumed to be assumed.
      ! So, it is stored as a vector of diagonal vector

      SUBROUTINE MassMatrixLinearPrism(
     1   xNode, yNode, zNode,    ! X,Y and Z coordinates of the 6 nodes
     2   elemData,        ! element property data
     3   Mlocal)          ! Mass matrix, stored as diagonal vector

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(6)     :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData
      DOUBLE PRECISION, DIMENSION(18)    :: Mlocal
      DOUBLE PRECISION, DIMENSION(18,18) :: Klocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4, dens, fact
      DOUBLE PRECISION :: N(6), dNdx(6), dNdy(6), dNdz(6)

      INTEGER, PARAMETER :: nGP=2, nlbf=6
      INTEGER :: DEGREE=1, ETYPE=6, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2, nsize=18

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      dens = elemData(3)

      ! Gauss point coordinates and weights
      call getGaussPointsPrism(nGP, gpts1, gpts2, gpts3, gwts)

      ! loop over Gauss points and compute 
      ! element mass matrix
      Klocal = 0.0
      Mlocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., ETYPE, DEGREE, param,
     1    xNode, yNode, zNode,
     2    N, dNdx, dNdy, dNdz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*Jac

        ! WRITE(*,*) " VOLUME = ", Jac, dvol, dens

        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b4 = (dens*dvol)*N(ii)

          DO jj=1, nlbf
            TJ   = (jj-1)*3+1
            TJp1 = TJ+1
            TJp2 = TJ+2

            fact = b4*N(jj)
            Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + fact
            Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + fact
            Klocal(TIp2, TJp2) = Klocal(TIp2, TJp2) + fact
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

      END SUBROUTINE MassMatrixLinearPrism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE StiffnessResidualElasticityLinearTetra(
     1   xNode, yNode, zNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   valC,            ! solution t_{n+af} for the 3 nodes
     5   valDotC,         ! d(solution)/dt at t_{n+am} for the 3 nodes
     6   Klocal, Flocal)  ! element stiffness and residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(4)     :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)     :: valC, valDotC
      DOUBLE PRECISION, DIMENSION(12,12) :: Klocal
      DOUBLE PRECISION, DIMENSION(6,6)   :: Dmat
      DOUBLE PRECISION, DIMENSION(12)    :: Flocal

      DOUBLE PRECISION, DIMENSION(6,12) :: Bmat
      DOUBLE PRECISION, DIMENSION(12,6) :: BmatTrans
      DOUBLE PRECISION :: Jac, dvol, thick, b1, b2, b3, b4
      DOUBLE PRECISION :: grad(3,3), strain(6), stress(6)
      DOUBLE PRECISION :: timefact, af, E, nu, bforce(3)
      DOUBLE PRECISION :: N(4), dNdx(4), dNdy(4), dNdz(4)

      INTEGER, PARAMETER :: nGP=8, nlbf=4
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      E  = elemData(1)
      nu = elemData(2)

      bforce(1) = elemData(4)
      bforce(2) = elemData(5)
      bforce(3) = elemData(6)

      b1 = E/((1.0+nu)*(1.0-2.0*nu))
      b2 = (1.0-2.0*nu)/2.0

      Dmat = 0.0
      Dmat(1,1) = b1*(1.0-nu); Dmat(1,2) = b1*nu;       Dmat(1,3) = b1*nu
      Dmat(2,1) = b1*nu;       Dmat(2,2) = b1*(1.0-nu); Dmat(2,3) = b1*nu
      Dmat(3,1) = b1*nu;       Dmat(3,2) = b1*nu;       Dmat(3,3) = b1*(1.0-nu)
      Dmat(4,4) = b1*b2
      Dmat(5,5) = b1*b2
      Dmat(6,6) = b1*b2

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      ! call getGaussPointsTetra(nGP, gpts1, gpts2, gpts3, gwts)
      gpts1(1) = 0.25; gpts2(1) = 0.25;  gpts3(1) = 0.25;      gwts(1)  = 1.0/6.0

      ! loop over Gauss points and compute 
      ! element stiffness and force vectors
      Klocal = 0.0;      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., 1, degree, param,
     1    xNode, yNode, zNode,
     2    N, dNdx, dNdy, dNdz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*Jac

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        ! compute the gradient of displacement first
        grad = 0.0
        DO ii=1, nlbf
          jj = (ii-1)*3
          b1 = valC(jj+1)
          b2 = valC(jj+2)
          b3 = valC(jj+3)

          grad(1,1) = grad(1,1) + b1*dNdx(ii)
          grad(1,2) = grad(1,2) + b1*dNdy(ii)
          grad(1,3) = grad(1,3) + b1*dNdz(ii)

          grad(2,1) = grad(2,1) + b2*dNdx(ii)
          grad(2,2) = grad(2,2) + b2*dNdy(ii)
          grad(2,3) = grad(2,3) + b2*dNdz(ii)

          grad(3,1) = grad(3,1) + b3*dNdx(ii)
          grad(3,2) = grad(3,2) + b3*dNdy(ii)
          grad(3,3) = grad(3,3) + b3*dNdz(ii)
        END DO

        ! now compute the strain
        strain(1) = grad(1,1)
        strain(2) = grad(2,2)
        strain(3) = grad(3,3)
        strain(4) = 0.5*(grad(1,2)+grad(2,1))
        strain(5) = 0.5*(grad(2,3)+grad(3,2))
        strain(6) = 0.5*(grad(1,3)+grad(3,1))

        stress = MATMUL(Dmat, strain)

        Bmat = 0.0
        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2
        
          Bmat(1,TI) = dNdx(ii); Bmat(1,TIp1) = 0.0;       Bmat(1,TIp2) = 0.0;
          Bmat(2,TI) = 0.0;       Bmat(2,TIp1) = dNdy(ii); Bmat(2,TIp2) = 0.0;
          Bmat(3,TI) = 0.0;       Bmat(3,TIp1) = 0.0;       Bmat(3,TIp2) = dNdz(ii);
          Bmat(4,TI) = dNdy(ii); Bmat(4,TIp1) = dNdx(ii); Bmat(4,TIp2) = 0.0;
          Bmat(5,TI) = 0.0;       Bmat(5,TIp1) = dNdz(ii); Bmat(5,TIp2) = dNdy(ii)
          Bmat(6,TI) = dNdz(ii); Bmat(6,TIp1) = 0.0;       Bmat(6,TIp2) = dNdx(ii)
        END DO

        BmatTrans = TRANSPOSE(Bmat)
        Bmat = MATMUL(Dmat, Bmat)
  
        Klocal = MATMUL(BmatTrans, Bmat)
        Klocal = dvol*Klocal

        ! contributions to the RHS vector from the residuals
        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b4 = dvol*N(ii)

          Flocal(TI)   = Flocal(TI)   + b4*bforce(1)
          Flocal(TIp1) = Flocal(TIp1) + b4*bforce(2)
          Flocal(TIp2) = Flocal(TIp2) + b4*bforce(3)
        END DO
      END DO

      END SUBROUTINE StiffnessResidualElasticityLinearTetra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Mass is assumed to be assumed.
      ! So, it is stored as a vector of diagonal vector

      SUBROUTINE MassMatrixLinearTetra(
     1   xNode, yNode, zNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   Mlocal)          ! Mass matrix, stored as diagonal vector

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(4)     :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData
      DOUBLE PRECISION, DIMENSION(12)    :: Mlocal
      DOUBLE PRECISION, DIMENSION(12,12) :: Klocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4, dens, fact
      DOUBLE PRECISION :: N(4), dNdx(4), dNdy(4), dN_dz(4)

      INTEGER, PARAMETER :: nGP=8, nlbf=4
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2, nsize=12

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      dens = elemData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      ! call getGaussPointsTetra(nGP, gpts1, gpts2, gpts3, gwts)
      gpts1(1) = 0.25; gpts2(1) = 0.25;  gpts3(1) = 0.25;      gwts(1)  = 1.0/6.0

      ! loop over Gauss points and compute 
      ! element mass matrix
      Klocal = 0.0
      Mlocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., 1, degree, param,
     1    xNode, yNode, zNode,
     2    N, dNdx, dNdy, dN_dz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*(Jac*dens)

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b4 = dvol*N(ii)

          DO jj=1, nlbf
            TJ   = (jj-1)*3+1
            TJp1 = TJ+1
            TJp2 = TJ+2

            fact = b4*N(jj)
            Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + fact
            Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + fact
            Klocal(TIp2, TJp2) = Klocal(TIp2, TJp2) + fact
          END DO
        END DO
      END DO

      Mlocal = 0.0
      DO ii=1,nsize
        fact = 0.0
        DO jj=1,nsize
          fact = fact + Klocal(ii,jj)
        END DO
        Mlocal(ii) = fact
      END DO

      END SUBROUTINE MassMatrixLinearTetra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Mass is assumed to be assumed.
      ! So, it is stored as a vector of diagonal vector

      SUBROUTINE MassMatrixLinearHexa(
     1   xNode, yNode, zNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   Mlocal)          ! Mass matrix, stored as diagonal vector

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(8)     :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData
      DOUBLE PRECISION, DIMENSION(24)    :: Mlocal
      DOUBLE PRECISION, DIMENSION(24,24) :: Klocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4, dens, fact
      DOUBLE PRECISION :: N(8), dNdx(8), dNdy(8), dN_dz(8)

      INTEGER, PARAMETER :: nGP=8, nlbf=8
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2, nsize=24

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      dens = elemData(3)

      ! Gauss point coordinates and weights
      ! call getQuadraturePointsHex(nGP, gpts1, gpts2, gpts3, gwts)

      ! loop over Gauss points and compute 
      ! element mass matrix
      Klocal = 0.0
      Mlocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

    !     call computeBasisFunctions3D(.FALSE., 1, degree, param,
    !  1    xNode, yNode, zNode,
    !  2    N, dNdx, dNdy, dN_dz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*(Jac*dens)

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b4 = dvol*N(ii)

          DO jj=1, nlbf
            TJ   = (jj-1)*3+1
            TJp1 = TJ+1
            TJp2 = TJ+2

            fact = b4*N(jj)
            Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + fact
            Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + fact
            Klocal(TIp2, TJp2) = Klocal(TIp2, TJp2) + fact
          END DO
        END DO
      END DO

      Mlocal = 0.0
      DO ii=1,nsize
        fact = 0.0
        DO jj=1,nsize
          fact = fact + Klocal(ii,jj)
        END DO
        Mlocal(ii) = fact
      END DO

      END SUBROUTINE MassMatrixLinearHexa


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE ResidualElasticityLinearTetra(
     1   xNode, yNode, zNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   valC,            ! solution t_{n+af} for the 3 nodes
     5   valDotC,         ! d(solution)/dt at t_{n+am} for the 3 nodes
     6   Flocal)  ! element stiffness and residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(4)     :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)     :: valC, valDotC
      DOUBLE PRECISION, DIMENSION(6,6)   :: Dmat
      DOUBLE PRECISION, DIMENSION(12)    :: Flocal

      DOUBLE PRECISION, DIMENSION(6,12) :: Bmat
      DOUBLE PRECISION, DIMENSION(12,6) :: BmatTrans
      DOUBLE PRECISION :: Jac, dvol, thick, b1, b2, b3, b4
      DOUBLE PRECISION :: grad(3,3), strain(6), stress(6)
      DOUBLE PRECISION :: timefact, af, E, nu, bforce(3)
      DOUBLE PRECISION :: N(4), dNdx(4), dNdy(4), dN_dz(4)

      INTEGER, PARAMETER :: nGP=1, nlbf=4
      INTEGER :: degree=1, type=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(3)

      ! material parameters, 
      E  = elemData(1)
      nu = elemData(2)

      bforce(1) = elemData(4)
      bforce(2) = elemData(5)
      bforce(3) = elemData(6)

      b1 = E/((1.0+nu)*(1.0-2.0*nu))
      b2 = (1.0-2.0*nu)/2.0

      Dmat = 0.0
      Dmat(1,1) = b1*(1.0-nu); Dmat(1,2) = b1*nu;       Dmat(1,3) = b1*nu
      Dmat(2,1) = b1*nu;       Dmat(2,2) = b1*(1.0-nu); Dmat(2,3) = b1*nu
      Dmat(3,1) = b1*nu;       Dmat(3,2) = b1*nu;       Dmat(3,3) = b1*(1.0-nu)
      Dmat(4,4) = b1*b2
      Dmat(5,5) = b1*b2
      Dmat(6,6) = b1*b2

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      gpts1(1) = 0.25; gpts2(1) = 0.25;  gpts3(1) = 0.25
      gwts(1)  = 1.0/6.0

      ! loop over Gauss points and compute 
      ! element stiffness and force vectors
      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., 1, degree, param,
     1    xNode, yNode, zNode,
     2    N, dNdx, dNdy, dN_dz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Elasticity" 
        END IF

        dvol = gwts(gp)*Jac

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        ! compute the gradient of displacement first
        grad = 0.0
        DO ii=1, nlbf
          jj = (ii-1)*3
          b1 = valC(jj+1)
          b2 = valC(jj+2)
          b3 = valC(jj+3)

          grad(1,1) = grad(1,1) + b1*dNdx(ii)
          grad(1,2) = grad(1,2) + b1*dNdy(ii)
          grad(1,3) = grad(1,3) + b1*dN_dz(ii)

          grad(2,1) = grad(2,1) + b2*dNdx(ii)
          grad(2,2) = grad(2,2) + b2*dNdy(ii)
          grad(2,3) = grad(2,3) + b2*dN_dz(ii)

          grad(3,1) = grad(3,1) + b3*dNdx(ii)
          grad(3,2) = grad(3,2) + b3*dNdy(ii)
          grad(3,3) = grad(3,3) + b3*dN_dz(ii)
        END DO

        ! now compute the strain
        strain(1) = grad(1,1)
        strain(2) = grad(2,2)
        strain(3) = grad(3,3)
        strain(4) = 0.5*(grad(1,2)+grad(2,1))
        strain(5) = 0.5*(grad(2,3)+grad(3,2))
        strain(6) = 0.5*(grad(1,3)+grad(3,1))

        stress = MATMUL(Dmat, strain)

        Bmat = 0.0
        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2
        
          Bmat(1,TI) = dNdx(ii); Bmat(1,TIp1) = 0.0;       Bmat(1,TIp2) = 0.0;
          Bmat(2,TI) = 0.0;       Bmat(2,TIp1) = dNdy(ii); Bmat(2,TIp2) = 0.0;
          Bmat(3,TI) = 0.0;       Bmat(3,TIp1) = 0.0;       Bmat(3,TIp2) = dN_dz(ii);
          Bmat(4,TI) = dNdy(ii); Bmat(4,TIp1) = dNdx(ii); Bmat(4,TIp2) = 0.0;
          Bmat(5,TI) = 0.0;       Bmat(5,TIp1) = dN_dz(ii); Bmat(5,TIp2) = dNdy(ii)
          Bmat(6,TI) = dN_dz(ii); Bmat(6,TIp1) = 0.0;       Bmat(6,TIp2) = dNdx(ii)
        END DO

        BmatTrans = TRANSPOSE(Bmat)
  
        ! contributions to the RHS vector from the residuals
        DO ii=1, nlbf
          TI   = (ii-1)*3+1
          TIp1 = TI+1
          TIp2 = TI+2

          b1 = dvol*dNdx(ii)
          b2 = dvol*dNdy(ii)
          b3 = dvol*dN_dz(ii)
          b4 = dvol*N(ii)

          Flocal(TI)   = Flocal(TI)   + b4*bforce(1)
          Flocal(TIp1) = Flocal(TIp1) + b4*bforce(2)
          Flocal(TIp2) = Flocal(TIp2) + b4*bforce(3)

          Flocal(TI)   = Flocal(TI)   - ( b1*stress(1) + b2*stress(4) + b3*stress(6) )
          Flocal(TIp1) = Flocal(TIp1) - ( b1*stress(4) + b2*stress(2) + b3*stress(5) )
          Flocal(TIp2) = Flocal(TIp2) - ( b1*stress(6) + b2*stress(5) + b3*stress(3) )

        END DO
      END DO

      END SUBROUTINE ResidualElasticityLinearTetra

      END MODULE ElementUtilitiesElasticity3D
