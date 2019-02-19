! Subroutines for calculating the stiffness
! and residual (or force) for the simple
! Lagrange family of elements for
!
! Incompressible Navier-Stokes using
! the stabilised formulation
!
! Author: Dr. Chennakesava Kadapa
! Date  : 25-Oct-2017
! Place : Swansea, UK
!
!
!

      MODULE ElementUtilitiesCFD

      USE ElementUtilitiesBasisFuncs

      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE StiffnessResidualIncNavStokesLinearTriaFI(
     1   xNode, yNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   timeData,        ! time integration data
     4   valC,            ! velocity at t_{n+af} for the 3 nodes
     5   valDotC,         ! acceleration at t_{n+am} for the 3 nodes
     6   velP,            ! velocity at t_{n-1} for the 3 nodes
     7   presN,           ! pressure at t_{n+1} for the 3 nodes
     8   Klocal, Flocal)  ! element stiffness and residual

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3)   :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:)   :: elemData, timeData
      DOUBLE PRECISION, DIMENSION(:)   :: valC, valDotC
      DOUBLE PRECISION, DIMENSION(:)   :: velP, presN
      DOUBLE PRECISION, DIMENSION(6,6) :: Klocal
      DOUBLE PRECISION, DIMENSION(6)   :: Flocal

      DOUBLE PRECISION :: Jac, dvol, thick, b1, b2, b3, b4
      DOUBLE PRECISION :: vel(2), velPrev(2), Du(2), dp(2), pres=0.0
      DOUBLE PRECISION :: grad(2,2), gradN(2,2), stress(2,2)
      DOUBLE PRECISION :: Dj(2,3), tau(3)
      DOUBLE PRECISION :: area, timefact, af, rho, mu, bforce(2)
      DOUBLE PRECISION :: Da, Db, fact, fact2, urdr, urdr2
      DOUBLE PRECISION :: xx, yy
      DOUBLE PRECISION :: N(3), dN_dx(3), dN_dy(3)

      INTEGER :: degree=1, type=1, nlbf=3, nGP=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2

      DOUBLE PRECISION :: gpts1(1), gpts2(1)
      DOUBLE PRECISION :: gwts(1), param(2)

      ! material parameters, 
      rho  = elemData(1)
      mu = elemData(2)

      bforce(1) = elemData(4)
      bforce(2) = elemData(5)

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
          STOP " Negative Jacobian for the Tria element in Elasticity" 
        END IF

        dvol = gwts(gp)*(Jac*thick)

        ! WRITE(*,*) " VOLUME = ", Jac, dvol

        ! compute the gradient of displacement first
        xx = 0.0; yy = 0.0
        vel  = 0.0; velPrev = 0.0; grad = 0.0;  Du = 0.0
        velDot = 0.0
        pres = 0.0; dp = 0.0
        DO ii=1, nlbf
          xx = xx + xNode(ii)*N(ii)
          yy = yy + yNode(ii)*N(ii)

          TI = (ii-1)*2 + 1

          b1 = valC(TI)
          b2 = valC(TIp1)

          vel(0)    =  vel(0) + b1*N(ii)
          vel(1)    =  vel(1) + b2*N(ii)

          grad(1,1) = grad(1,1) + b1*dN_dx(ii)
          grad(1,2) = grad(1,2) + b1*dN_dy(ii)
          grad(2,1) = grad(2,1) + b2*dN_dx(ii)
          grad(2,2) = grad(2,2) + b2*dN_dy(ii)

          velPrev(0) = velPrev(0) + velP(TI)*N(ii)
          velPrev(1) = velPrev(1) + velP(TIp1)*N(ii)

          pres   = pres  + presN(ii)*N(ii)
          dp(0)  = dp(0) + presN(ii)*dN_dx(ii)
          dp(1)  = dp(1) + presN(ii)*dN_dy(ii)

          velDot(0) = velDot(0) + velDotC(TI)*N(ii)
          velDot(1) = velDot(1) + velDotC(TIp1)*N(ii)
        END DO

        stress = mu*grad; !this is pseudo-stress
        !stress = mu*(grad+grad.transpose());
        stress(1,1) = stress(1,1) - pres
        stress(2,2) = stress(2,2) - pres

          force = 0.0
          gradTvel = grad*vel ;

          res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;
          res2(2) = rho*(velDot(2) + gradTvel(2) - force(2)) ;

          rStab(1) = res2(1) - mu*Du(1) + dp(1) ;
          rStab(2) = res2(2) - mu*Du(2) + dp(2) ;

          IF(axsy .EQV. .TRUE.) THEN
            rad = xx;

            urdr  = vel(1)/rad;
            urdr2 = urdr/rad;
            dvol *= (2.0*PI*rad);

            rStab(1) -= mu*(grad(1,1)/rad - urdr2 );
            rStab(2) -= mu*(grad(2,1)/rad );
          END IF

          velTemp(1) = velPrev(1);
          velTemp(2) = velPrev(2);
          velTemp(3) = 0.0;

          !evaluateStabParams_algo1(&velTemp(0), h, rho, mu, dt,  beta, tau);
          !evaluateStabParams_algo2(&velTemp(0), h, rho, mu, dt,  beta, tau);
          !evaluateStabParams_algo3(velTemp, matG, dt, rho, mu, CI, tau);
          ! tau[0] *= elmDat[8];  // SUPG
          ! tau[1] *= elmDat[9];  // PSPG
          ! tau[2] *= elmDat[10]; // LSIC

          DO ii=1, nlbf
            TI   = (ii-1)*ndof+1;
            TIp1 = TI+1;
            TIp2 = TI+2;

            b1 = dN_dx(ii)*dvol;
            b2 = dN_dy(ii)*dvol;
            b4 = N(ii)*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;
            b8 = af*b4;

            Da = rho*(vel(1)*b1 + vel(2)*b2)*tau(1);

            DO jj=1, nlbf
              TJ   = (jj-1)*ndof+1;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              fact2 = rho*acceFact*N(jj);

              !time acceleration term
              fact = b4*fact2 ;

              !diffusion term
              fact = fact + ( b5*dN_dx(jj)+b6*dN_dy(jj) );

              Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + fact;
              Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + fact;

              Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + ( b5*dN_dx(jj) );
              Klocal(TI,   TJp1) = Klocal(TI,   TJp1) + ( b5*dN_dy(jj) );
              Klocal(TIp1, TJ)   = Klocal(TIp1, TJ)   + ( b6*dN_dx(jj) );
              Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + ( b6*dN_dy(jj) );

              ! convection term

              gradN = grad*(rho*N(jj));

              Db = rho*(vel(1)*dN_dx(jj) + vel(2)*dN_dy(jj));

              gradN(1,1) = gradN(1,1) + Db;
              gradN(2,2) = gradN(2,2) + Db;

              Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + b8*gradN(1,1)
              Klocal(TI,   TJp1) = Klocal(TI,   TJp1) + b8*gradN(1,2)
              Klocal(TIp1, TJ)   = Klocal(TIp1, TJ)   + b8*gradN(2,1)
              Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + b8*gradN(2,2)

              ! pressure term
              Klocal(TI,   TJp2) = Klocal(TI,   TJp2) - b1*N(jj)
              Klocal(TIp1, TJp2) = Klocal(TIp1, TJp2) - b2*N(jj)

              ! continuity equation
              Klocal(TIp2, TJ)   = Klocal(TIp2, TJ)   + b8*dN_dx(jj)
              Klocal(TIp2, TJp1) = Klocal(TIp2, TJp1) + b8*dN_dy(jj)

              ! SUPG and PSPG stabilisation terms

              gradN = gradN*af;

              Dj(1,1) = gradN(1,1) + fact2;
              Dj(1,2) = gradN(1,2);
              Dj(1,3) = dN_dx(jj);
              Dj(2,1) = gradN(2,1);
              Dj(2,2) = gradN(2,2) + fact2;
              Dj(2,3) = dN_dy(jj);

              IF(axsy .EQV. .TRUE.) THEN
                Dj(1,1) = Dj(1,1) - muTaf*(dN_dx(jj)/rad - N(jj)/rad/rad);
                Dj(2,2) = Dj(2,2) - muTaf*(dN_dx(jj)/rad);
              END IF

              ! SUPG
              Klocal(TI, TJ)     = Klocal(TI, TJ)     + Da*Dj(1,1);
              Klocal(TI, TJp1)   = Klocal(TI, TJp1)   + Da*Dj(1,2);
              Klocal(TI, TJp2)   = Klocal(TI, TJp2)   + Da*Dj(1,3);

              Klocal(TIp1, TJ)   = Klocal(TI, TJ)     + Da*Dj(2,1);
              Klocal(TIp1, TJp1) = Klocal(TI, TJp1)   + Da*Dj(2,2);
              Klocal(TIp1, TJp2) = Klocal(TI, TJp2)   + Da*Dj(2,3);

              Klocal(TI,   TJ)   = Klocal(TI, TJ)     + ( (tau(1)*af) * b1 * rStab(1) * N(jj) );
              Klocal(TI,   TJp1) = Klocal(TI, TJp1)   + ( (tau(1)*af) * b2 * rStab(1) * N(jj) );
              Klocal(TIp1, TJ)   = Klocal(TIp1, TJ)   + ( (tau(1)*af) * b1 * rStab(2) * N(jj) );
              Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + ( (tau(1)*af) * b2 * rStab(2) * N(jj) );

              ! PSPG stabilisation
              Klocal(TIp2, TJ)   = Klocal(TIp2, TJ)   + (b1*Dj(1,1) + b2*Dj(2,1))*tau(2);
              Klocal(TIp2, TJp1) = Klocal(TIp2, TJp1) + (b1*Dj(1,2) + b2*Dj(2,2))*tau(2);
              Klocal(TIp2, TJp2) = Klocal(TIp2, TJp2) + (b1*Dj(1,3) + b2*Dj(2,3))*tau(2);

              ! LSIC stabilisation

              fact = af*rho*tau(3);

              Klocal(TI,   TJ)   = Klocal(TI,   TJ)   + (b1*fact*dN_dx(jj));
              Klocal(TI,   TJp1) = Klocal(TI,   TJp1) + (b1*fact*dN_dy(jj));

              Klocal(TIp1, TJ)   = Klocal(TIp1, TJ)   + (b2*fact*dN_dx(jj));
              Klocal(TIp1, TJp1) = Klocal(TIp1, TJp1) + (b2*fact*dN_dy(jj));

              IF(axsy .EQV. .TRUE.) THEN
                  ! diffusion term
                  Klocal(TI, TJ)     = Klocal(TI,   TJ)   + (b4 * (mu/rad/rad) * (af*N(jj)) );
                  Klocal(TI, TJp2)   = Klocal(TI,   TJp2) - (b4 * N(jj)/rad);

                  ! continuity equation
                  Klocal(TIp2, TJ)   = Klocal(TIp2, TJ)   + (b4 * af*N(jj)/rad);
              END IF
            END DO

            Flocal(TI)   = Flocal(TI)   - (b4*res2(1) + b1*stress(1,1) + b2*stress(1,2) );
            Flocal(TIp1) = Flocal(TIp1) - (b4*res2(2) + b1*stress(2,1) + b2*stress(2,2) );
            Flocal(TIp2) = Flocal(TIp2) - (b4*TRACE(grad));

            ! SUPG stabilisation terms
            Flocal(TI)   = Flocal(TI)   - Da*rStab(1);
            Flocal(TIp1) = Flocal(TIp1) - Da*rStab(2);
             
            ! PSPG stabilisation terms
            Flocal(TIp2) = Flocal(TIp2) - (tau(2)*(b1*rStab(1)+b2*rStab(2)));

            ! LSIC stabilisation terms
            fact2 = tau(3)*rho*TRACE(grad);

            Flocal(TI)   = Flocal(TI)   - b1*fact2;
            Flocal(TIp1) = Flocal(TIp1) - b2*fact2;

            IF(axsy .EQV. .TRUE.) THEN
              Flocal(TI)   = Flocal(TI)   - (b4 * (mu/rad/rad) * vel(1) )
              Flocal(TI)   = Flocal(TI)   + (b4 * pres/rad)
              Flocal(TIp2) = Flocal(TIp2) - (b4 * vel(1)/rad)
            END IF
          END DO
      END DO

      END SUBROUTINE StiffnessResidualIncNavStokesLinearTria


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Mass is assumed to be lumped.
      ! So, it is stored as a vector of diagonal vector

      SUBROUTINE MassMatrixTria3Dof1(
     1   xNode, yNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   Mlocal)          ! Mass matrix, stored as diagonal vector

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3)     :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData
      DOUBLE PRECISION, DIMENSION(3)     :: Mlocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4, dens, fact

      INTEGER :: ii, jj

      ! compute element volume
      dvol = abs(0.5*(xNode(1)*(yNode(2)-yNode(3)) + xNode(2)*(yNode(3)-yNode(1)) + xNode(3)*(yNode(1)-yNode(2))));

      ! WRITE(*,*) " VOLUME = ", Jac, dvol, dens

      ! Compute element mass matrix

      ! material parameters
      dens = elemData(1)

      fact = dens*dvol/3.0;

      !Mlocal = 0.0
      DO ii=1,3
        Mlocal(ii) = fact
      END DO

      END SUBROUTINE MassMatrixTria3Dof1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Mass is assumed to be lumped.
      ! So, it is stored as a vector of diagonal vector

      SUBROUTINE MassMatrixTria6Dof1(
     1   xNode, yNode,    ! X and Y coordinates of the 3 nodes
     2   elemData,        ! element property data
     3   Mlocal)          ! Mass matrix, stored as diagonal vector

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(6)     :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:)     :: elemData
      DOUBLE PRECISION, DIMENSION(6)     :: Mlocal

      DOUBLE PRECISION :: Jac, dvol, b1, b2, b3, b4, dens, fact
      DOUBLE PRECISION :: N(6), dNdx(6), dNdy(6)

      INTEGER, PARAMETER :: nGP=3, nlbf=6
      INTEGER :: degree=2, type=1, gp, ii, jj

      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP)
      DOUBLE PRECISION :: gwts(nGP), param(2)

      ! Gauss point coordinates and weights
      call getGaussPointsTria(nGP, gpts1, gpts2, gwts)

      ! compute element volume

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
          STOP " Negative Jacobian for the Tria element in CFD" 
        END IF

        dvol = gwts(gp)*Jac

        ! WRITE(*,*) " VOLUME = ", Jac, dvol, dens
      END DO

      ! WRITE(*,*) Klocal

      ! Compute element mass matrix

      ! material parameters
      dens = elemData(1)

      fact = dens*dvol/6.0;

      !Mlocal = 0.0
      DO ii=1,6
        Mlocal(ii) = fact
      END DO

      END SUBROUTINE MassMatrixTria6Dof1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE StiffnessResidualIncNavStokesLinearTetra(
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
      DOUBLE PRECISION :: N(4), dN_dx(4), dN_dy(4), dN_dz(4)

      INTEGER :: degree=1, type=1, nlbf=4, nGP=1, gp, ii, jj
      INTEGER :: TI, TIp1, TIp2, TJ, TJp1, TJp2

      DOUBLE PRECISION :: gpts1(1), gpts2(1), gpts3(1)
      DOUBLE PRECISION :: gwts(1), param(3)

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
      Klocal = 0.0;      Flocal = 0.0

      DO gp=1,nGP
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

        call computeBasisFunctions3D(.FALSE., 1, degree, param,
     1    xNode, yNode, zNode,
     2    N, dN_dx, dN_dy, dN_dz, Jac)

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

          grad(1,1) = grad(1,1) + b1*dN_dx(ii)
          grad(1,2) = grad(1,2) + b1*dN_dy(ii)
          grad(1,3) = grad(1,3) + b1*dN_dz(ii)

          grad(2,1) = grad(2,1) + b2*dN_dx(ii)
          grad(2,2) = grad(2,2) + b2*dN_dy(ii)
          grad(2,3) = grad(2,3) + b2*dN_dz(ii)

          grad(3,1) = grad(3,1) + b3*dN_dx(ii)
          grad(3,2) = grad(3,2) + b3*dN_dy(ii)
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
        
          Bmat(1,TI) = dN_dx(ii); Bmat(1,TIp1) = 0.0;       Bmat(1,TIp2) = 0.0;
          Bmat(2,TI) = 0.0;       Bmat(2,TIp1) = dN_dy(ii); Bmat(2,TIp2) = 0.0;
          Bmat(3,TI) = 0.0;       Bmat(3,TIp1) = 0.0;       Bmat(3,TIp2) = dN_dz(ii);
          Bmat(4,TI) = dN_dy(ii); Bmat(4,TIp1) = dN_dx(ii); Bmat(4,TIp2) = 0.0;
          Bmat(5,TI) = 0.0;       Bmat(5,TIp1) = dN_dz(ii); Bmat(5,TIp2) = dN_dy(ii)
          Bmat(6,TI) = dN_dz(ii); Bmat(6,TIp1) = 0.0;       Bmat(6,TIp2) = dN_dx(ii)
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

      END SUBROUTINE StiffnessResidualIncNavStokesLinearTetra

      END MODULE ElementUtilitiesCFD
