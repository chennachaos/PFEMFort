! Subroutines for calculating the stiffness
! and residual (or force) for the simple
! Lagrange family of elements for
!
! Poisson equation
!
! Author: Dr. Chennakesava Kadapa
! Date  : 25-Oct-2017
! Place : Swansea, UK
!
!
!

      MODULE ElementUtilitiesPoisson

      USE ElementUtilitiesBasisFuncs

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE StiffnessResidualPoissonLinearTria(
     1   xNode, yNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   valC,            ! solution t_{n+af} for the 3 nodes
     5   valDotC,         ! d(solution)/dt at t_{n+am} for the 3 nodes
     6   Klocal, Flocal)  ! element stiffness and residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:) :: elemData, timeData
      DOUBLE PRECISION :: xNode(3), yNode(3)
      DOUBLE PRECISION :: valC(3), valDotC(3)
      DOUBLE PRECISION :: Klocal(3,3), Flocal(3)

      DOUBLE PRECISION :: kx, ky, b1, b2, b4, du(2)
      DOUBLE PRECISION :: area, timefact, af, Jac, dvol, force
      DOUBLE PRECISION :: N(3), dN_dx(3), dN_dy(3)

      INTEGER :: degree=1, nlbf=3, nGP=1, gp, ii, jj

      DOUBLE PRECISION :: gpts1(1), gpts2(1), gwts(1), param(2)

      ! material parameters, 
      ! for example, coefficients of heat conduction
      kx = elemData(1);       ky = elemData(2)

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      gpts1(1) = 1.0/3.0; gpts2(1) = 1.0/3.0
      gwts(1)  = 0.5

      ! loop over Gauss points and compute 
      ! element stiffness and force vectors
      Klocal = 0.0;      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)

        call computeBasisFunctions2D(.FALSE., 1, degree, param,
     1    xNode, yNode,
     2    N, dN_dx, dN_dy, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tria element in Poisson" 
        END IF

        dvol = gwts(gp) * Jac;
          
        du = 0.0
        DO ii=1, nlbf
          du(1) = du(1) + valC(ii)*dN_dx(ii)
          du(2) = du(2) + valC(ii)*dN_dy(ii)
        END DO

        force = 0.0

        Do ii=1, nlbf
          b1 = dN_dx(ii)*dvol;
          b2 = dN_dy(ii)*dvol;
          b4 = N(ii)*dvol;

          Flocal(ii) =  Flocal(ii) + b4*force - b1*du(1) - b2*du(2)

          DO jj=1, nlbf
            Klocal(ii, jj) = Klocal(ii, jj) +
     1        af*( b1*(kx*dN_dx(jj)) + b2*(ky*dN_dy(jj)) )
          END DO
        END DO
      END DO

    !   WRITE(*,*) Klocal

      END SUBROUTINE StiffnessResidualPoissonLinearTria


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE StiffnessResidualPoissonLinearTetra(
     1   xNode, yNode, zNode,    ! X, Y and Z coordinates of the 4 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   valC,            ! solution t_{n+af} for the 4 nodes
     5   valDotC,         ! d(solution)/dt at t_{n+am} for the 4 nodes
     6   Klocal, Flocal)  ! element stiffness and residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:) :: elemData, timeData
      DOUBLE PRECISION :: xNode(4), yNode(4), zNode(4)
      DOUBLE PRECISION :: valC(4), valDotC(4)
      DOUBLE PRECISION :: Klocal(4,4), Flocal(4)

      DOUBLE PRECISION :: kx, ky, kz, b1, b2, b3, b4, du(3)
      DOUBLE PRECISION :: area, timefact, af, Jac, dvol, force
      DOUBLE PRECISION :: N(4), dN_dx(4), dN_dy(4), dN_dz(4)

      INTEGER :: degree=1, type=4, nlbf=4, nGP=1, gp, ii, jj

      DOUBLE PRECISION :: gpts1(1), gpts2(1), gpts3(1)
      DOUBLE PRECISION :: gwts(1), param(3)

      ! material parameters, 
      ! for example, coefficients of heat conduction
      kx = elemData(1); ky = elemData(2); kz = elemData(3)

      ! time integration parameters
      af   = timeData(2)
      timefact = timeData(3)

      ! Gauss point coordinates and weights
      ! ONE Gauss point is assumed here
      gpts1(1) = 0.25; gpts2(1) = 0.25; gpts3(1) = 0.25
      gwts(1)  = 1.0/6.0

      ! loop over Gauss points and compute 
      ! element stiffness and force vectors
      Klocal = 0.0;      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., type, degree, param,
     1    xNode, yNode, zNode,
     2    N, dN_dx, dN_dy, dN_dz, Jac)

        IF(Jac < 0.0) THEN
          STOP " Negative Jacobian for the Tet element in Poisson" 
        END IF

        dvol = gwts(gp) * Jac

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        du = 0.0
        DO ii=1, nlbf
          du(1) = du(1) + valC(ii)*dN_dx(ii)
          du(2) = du(2) + valC(ii)*dN_dy(ii)
          du(3) = du(3) + valC(ii)*dN_dz(ii)
        END DO

        force = -6.0

        Do ii=1, nlbf
          b1 = dN_dx(ii)*dvol;
          b2 = dN_dy(ii)*dvol;
          b3 = dN_dz(ii)*dvol;
          b4 = N(ii)*dvol;

          Flocal(ii) =  Flocal(ii) + b4*force
          Flocal(ii) =  Flocal(ii) - b1*du(1) - b2*du(2) - b3*du(3)

          DO jj=1, nlbf
            Klocal(ii, jj) = Klocal(ii, jj) +
     1       af*( b1*(kx*dN_dx(jj)) +
     2       b2*(ky*dN_dy(jj))+ b3*(kz*dN_dz(jj)) )
          END DO
        END DO
      END DO

    !   WRITE(*,*) Klocal

      END SUBROUTINE StiffnessResidualPoissonLinearTetra

      END MODULE ElementUtilitiesPoisson
