!
! Utilities for Quadrature points
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 12-Nov-2017
! Place : Swansea, UK
!
!

      MODULE QuadratureUtilities

      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getGaussPoints1D(ngp, gausspoints, gaussweights)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ngp
      DOUBLE PRECISION, DIMENSION(:)  :: gausspoints, gaussweights

      SELECT CASE (ngp)
        CASE (1)  ! 1 Point quadrature rule

            gausspoints(1) = 0.0;        gaussweights(1) = 2.0;

        CASE (2)  ! 2 Point quadrature rule

            gausspoints(1) = -0.577350269189626;    gaussweights(1) = 1.0;
            gausspoints(2) =  0.577350269189626;    gaussweights(2) = 1.0;

        CASE (3)  ! 3 Point quadrature rule

            gausspoints(1) = -0.774596669241483;  gaussweights(1) = 0.555555555555556;
            gausspoints(2) =  0.0;                gaussweights(2) = 0.888888888888889;
            gausspoints(3) =  0.774596669241483;  gaussweights(3) = 0.555555555555556;

       CASE (4)  ! 4 Point quadrature rule

            gausspoints(1) = -0.861136311594953;    gaussweights(1) = 0.347854845137454;
            gausspoints(2) = -0.339981043584856;    gaussweights(2) = 0.652145154862546;
            gausspoints(3) =  0.339981043584856;    gaussweights(3) = 0.652145154862546;
            gausspoints(4) =  0.861136311594953;    gaussweights(4) = 0.347854845137454;

        CASE (5) ! 5 Point quadrature rule

            gausspoints(1) = -0.906179845938664;    gaussweights(1) = 0.236926885056189;
            gausspoints(2) = -0.538469310105683;    gaussweights(2) = 0.478628670499366;
            gausspoints(3) =  0.0;                  gaussweights(3) = 0.568888888888889;
            gausspoints(4) =  0.538469310105683;    gaussweights(4) = 0.478628670499366;
            gausspoints(5) =  0.906179845938664;    gaussweights(5) = 0.236926885056189;

        CASE (6) ! 6 Point quadrature rule

            gausspoints(1) = -0.932469514203152;     gaussweights(1) = 0.171324492379170;
            gausspoints(2) = -0.661209386466265;     gaussweights(2) = 0.360761573048139;
            gausspoints(3) = -0.238619186083197;     gaussweights(3) = 0.467913934572691;
            gausspoints(4) =  0.238619186083197;     gaussweights(4) = 0.467913934572691;
            gausspoints(5) =  0.661209386466265;     gaussweights(5) = 0.360761573048139;
            gausspoints(6) =  0.932469514203152;     gaussweights(6) = 0.171324492379170;

        CASE (7) ! 7 Point quadrature rule

            gausspoints(1) = -0.9491079123427585245261897;	gaussweights(1) = 0.1294849661688696932706114 ;
            gausspoints(2) = -0.7415311855993944398638648; 	gaussweights(2) = 0.2797053914892766679014678 ;
            gausspoints(3) = -0.4058451513773971669066064; 	gaussweights(3) = 0.3818300505051189449503698 ;
            gausspoints(4) =  0.0 ;			        gaussweights(4) = 0.4179591836734693877551020 ;
            gausspoints(5) =  0.4058451513773971669066064; 	gaussweights(5) = 0.3818300505051189449503698 ;
            gausspoints(6) =  0.7415311855993944398638648; 	gaussweights(6) = 0.2797053914892766679014678 ;
            gausspoints(7) =  0.9491079123427585245261897; 	gaussweights(7) = 0.1294849661688696932706114 ;

        CASE (8) ! 8 Point quadrature rule

            gausspoints(1) = -0.96028986; 	gaussweights(1) = 0.10122854 ;
            gausspoints(2) = -0.79666648; 	gaussweights(2) = 0.22238103 ;
            gausspoints(3) = -0.52553241; 	gaussweights(3) = 0.31370665 ;
            gausspoints(4) = -0.18343464; 	gaussweights(4) = 0.36268378 ;
            gausspoints(5) =  0.18343464; 	gaussweights(5) = 0.36268378 ;
            gausspoints(6) =  0.52553241; 	gaussweights(6) = 0.31370665 ;
            gausspoints(7) =  0.79666648; 	gaussweights(7) = 0.22238103 ;
            gausspoints(8) =  0.96028986; 	gaussweights(8) = 0.10122854 ;

        CASE DEFAULT
            STOP " invalid value of 'ngp' ! "

      END SELECT

      END SUBROUTINE getGaussPoints1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getGaussPointsQuad(ngp, gpts1, gpts2, gwts)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ngp
      DOUBLE PRECISION, DIMENSION(:) :: gpts1, gpts2, gwts

      INTEGER ::  nn, ii, jj, kk

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gpoints1, gweights1

      IF(ngp == 1) THEN
        nn = 1
      ELSE IF(ngp == 4) THEN
        nn = 2
      ELSE IF(ngp == 9) THEN
        nn = 3
      ELSE IF(ngp == 16) THEN
        nn = 4
      ELSE
        nn = 5
      ENDIF

      ALLOCATE(gpoints1(nn))
      ALLOCATE(gweights1(nn))

      call getGaussPoints1D(nn, gpoints1, gweights1)
  
      kk=1
      DO jj=1, nn
        DO ii=1, nn
          gpts1(kk) = gpoints1(ii)
          gpts2(kk) = gpoints1(jj)
          gwts(kk)  = gweights1(jj)*gweights1(ii);
        kk = kk+1
        END DO
      END DO

      DEALLOCATE(gpoints1)
      DEALLOCATE(gweights1)

      END SUBROUTINE getGaussPointsQuad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getGaussPointsHexa(ngp, gpts1, gpts2, gpts3, gwts)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ngp
      DOUBLE PRECISION, DIMENSION(:) :: gpts1, gpts2, gpts3, gwts

      INTEGER ::  nn, ii, jj, kk, ll

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gpoints1, gweights1

      IF(ngp == 1) THEN
        nn = 1
      ELSE IF(ngp == 8) THEN
        nn = 2
      ELSE IF(ngp == 27) THEN
        nn = 3
      ELSE IF(ngp == 64) THEN
        nn = 4
      ELSE
        nn = 5
      ENDIF

      ALLOCATE(gpoints1(nn))
      ALLOCATE(gweights1(nn))

      call getGaussPoints1D(nn, gpoints1, gweights1)

      ll=1
      DO kk=1, nn
        DO jj=1, nn
          DO ii=1, nn
            gpts1(ll) = gpoints1(ii)
            gpts2(ll) = gpoints1(jj)
            gpts3(ll) = gpoints1(kk)
            gwts(ll)  = gweights1(kk)*gweights1(jj)*gweights1(ii);
            ll = ll+1
          END DO
        END DO
      END DO

      DEALLOCATE(gpoints1)
      DEALLOCATE(gweights1)

      END SUBROUTINE getGaussPointsHexa
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getGaussPointsTria(ngp, gpts1, gpts2, gwts)
      ! weights are normalized to calculate the exact area of the triangle
      ! i.e. each weight is divided by 2.0

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ngp
      DOUBLE PRECISION, DIMENSION(:) :: gpts1, gpts2, gwts
      DOUBLE PRECISION :: r1d3 = 1.0/3.0, a1, fact=0.5;

      SELECT CASE (ngp)
        case (1)  ! 1 Point quadrature rule

            gpts1(1) = r1d3;        gpts2(1) = r1d3;        gwts(1) = fact*1.0;

        case (3)  ! 3 Point quadrature rule

            a1 = 1.0/3.0;

            !gpts1(1) = 0.5;        gpts2(1) = 0.0;        gwts(1) = fact*a1;
            !gpts1(2) = 0.5;        gpts2(2) = 0.5;        gwts(2) = fact*a1;
            !gpts1(3) = 0.0;        gpts2(3) = 0.5;        gwts(3) = fact*a1;

            gpts1(1) = 1.0/6.0;        gpts2(1) = 1.0/6.0;        gwts(1) = fact*a1;
            gpts1(2) = 1.0/6.0;        gpts2(2) = 4.0/6.0;        gwts(2) = fact*a1;
            gpts1(3) = 4.0/6.0;        gpts2(3) = 1.0/6.0;        gwts(3) = fact*a1;

        case (4)  ! 4 Point quadrature rule

            a1 = 25.0/48.0;

            gpts1(1) = r1d3;       gpts2(1) = r1d3;       gwts(1) = fact*(-27.0/48.0);
            gpts1(2) = 0.6;        gpts2(2) = 0.2;        gwts(2) = fact*a1;
            gpts1(3) = 0.2;        gpts2(3) = 0.6;        gwts(3) = fact*a1;
            gpts1(4) = 0.2;        gpts2(4) = 0.2;        gwts(4) = fact*a1;

        case (6)  ! 6 Point quadrature rule

            gpts1(1) = 0.10810301816807022736;    gpts2(1) = 0.44594849091596488632;    gwts(1) = fact*0.22338158967801146570;
            gpts1(2) = 0.44594849091596488632;    gpts2(2) = 0.10810301816807022736;    gwts(2) = fact*0.22338158967801146570;
            gpts1(3) = 0.44594849091596488632;    gpts2(3) = 0.44594849091596488632;    gwts(3) = fact*0.22338158967801146570;
            gpts1(4) = 0.81684757298045851308;    gpts2(4) = 0.09157621350977074346;    gwts(4) = fact*0.10995174365532186764;
            gpts1(5) = 0.09157621350977074346;    gpts2(5) = 0.81684757298045851308;    gwts(5) = fact*0.10995174365532186764;
            gpts1(6) = 0.09157621350977074346;    gpts2(6) = 0.09157621350977074346;    gwts(6) = fact*0.10995174365532186764;

        case (7)  ! 7 Point quadrature rule

            gpts1(1) = r1d3;                      gpts2(1) = r1d3;                      gwts(1) = fact*0.225;
            gpts1(2) = 0.79742698535308732240;    gpts2(2) = 0.10128650732345633880;    gwts(2) = fact*0.12593918054482715260;
            gpts1(3) = 0.10128650732345633880;    gpts2(3) = 0.79742698535308732240;    gwts(3) = fact*0.12593918054482715260;
            gpts1(4) = 0.10128650732345633880;    gpts2(4) = 0.10128650732345633880;    gwts(4) = fact*0.12593918054482715260;
            gpts1(5) = 0.05971587178976982045;    gpts2(5) = 0.47014206410511508977;    gwts(5) = fact*0.13239415278850618074;
            gpts1(6) = 0.47014206410511508977;    gpts2(6) = 0.05971587178976982045;    gwts(6) = fact*0.13239415278850618074;
            gpts1(7) = 0.47014206410511508977;    gpts2(7) = 0.47014206410511508977;    gwts(7) = fact*0.13239415278850618074;

        CASE DEFAULT
            STOP " invalid value of 'ngp' getGaussPointsTria ! "

      END SELECT

      END SUBROUTINE getGaussPointsTria
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getGaussPointsTetra(ngp, gpts1, gpts2, gpts3, gwts)
      ! weights are normalized to calculate the exact area of the tetrahedron
      ! i.e. each weight is divided by 6.0

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ngp
      DOUBLE PRECISION, DIMENSION(:)  :: gpts1, gpts2, gpts3, gwts
      DOUBLE PRECISION ::  fact=1.0/6.0

      SELECT CASE (ngp)

        case (1)

            gpts1(1) = 0.25;
            gpts2(1) = 0.25;
            gpts3(1) = 0.25;
            gwts(1) = fact;

        case (4)

            gpts1(1) = 0.5854101966249685;
            gpts1(2) = 0.1381966011250105;
            gpts1(3) = 0.1381966011250105;
            gpts1(4) = 0.1381966011250105;

            gpts2(1) = 0.1381966011250105;
            gpts2(2) = 0.1381966011250105;
            gpts2(3) = 0.1381966011250105;
            gpts2(4) = 0.5854101966249685;

            gpts3(1) = 0.1381966011250105;
            gpts3(2) = 0.1381966011250105;
            gpts3(3) = 0.5854101966249685;
            gpts3(4) = 0.1381966011250105;

            gwts(1)  = 0.2500000000000000*fact;
            gwts(2)  = 0.2500000000000000*fact;
            gwts(3)  = 0.2500000000000000*fact;
            gwts(4)  = 0.2500000000000000*fact;

        case (5)

            gpts1(1) = 0.2500000000000000;
            gpts1(2) = 0.5000000000000000;
            gpts1(3) = 0.1666666666666667;
            gpts1(4) = 0.1666666666666667;
            gpts1(5) = 0.1666666666666667;
 
            gpts2(1) = 0.2500000000000000;
            gpts2(2) = 0.1666666666666667;
            gpts2(3) = 0.1666666666666667;
            gpts2(4) = 0.1666666666666667;
            gpts2(5) = 0.5000000000000000;

            gpts3(1) = 0.2500000000000000;
            gpts3(2) = 0.1666666666666667;
            gpts3(3) = 0.1666666666666667;
            gpts3(4) = 0.5000000000000000;
            gpts3(5) = 0.1666666666666667;
 
            gwts(1)  = -0.8000000000000000*fact;
            gwts(2)  =  0.4500000000000000*fact;
            gwts(3)  =  0.4500000000000000*fact;
            gwts(4)  =  0.4500000000000000*fact;
            gwts(5)  =  0.4500000000000000*fact;

        CASE DEFAULT
            STOP " invalid value of 'ngp' in getGaussPointsTetra ! "

      END SELECT

      END SUBROUTINE getGaussPointsTetra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE getGaussPointsPrism(ngp, gpts1, gpts2, gpts3, gwts)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ngp
      DOUBLE PRECISION, DIMENSION(:)  :: gpts1, gpts2, gpts3, gwts
      DOUBLE PRECISION :: r1d3=1.0/3.0, r1d2=0.5, fact

      SELECT CASE (ngp)

        case (1)  ! 1 Point quadrature rule - 1 for triangle, 1 for the quad

            gpts1(1) = r1d3;
            gpts2(1) = r1d3;
            gpts3(1) = 0.0;
            gwts(1)  = r1d2*2.0; ! 2.0 for the weight in the normal direction

        case (2)  ! 2 Point quadrature rule - 1 for triangle, 2 for the quad

            gpts1(1) = r1d3;
            gpts2(1) = r1d3;
            gpts3(1) = -0.577350269189626;
            gwts(1)  = r1d2*1.0;

            gpts1(2) = r1d3;
            gpts2(2) = r1d3;
            gpts3(2) = 0.577350269189626;
            gwts(2)  = r1d2*1.0;

        case (3)  ! 3 Point quadrature rule - 3 for triangle, 1 for the quad

            fact = r1d2*r1d3*2.0

            gpts1(1) = 1.0/6.0;
            gpts2(1) = 1.0/6.0;
            gpts3(1) = 0.0;
            gwts(1)  = fact;

            gpts1(2) = 1.0/6.0;
            gpts2(2) = 4.0/6.0;
            gpts3(2) = 0.0;
            gwts(2)  = fact;

            gpts1(3) = 4.0/6.0;
            gpts2(3) = 1.0/6.0;
            gpts3(3) = 0.0;
            gwts(3)  = fact;

        case (6)  ! 6 Point quadrature rule - 3 for triangle, 2 for the quad

            fact = r1d2*r1d3*1.0;

            gpts1(1) = 1.0/6.0;
            gpts2(1) = 1.0/6.0;
            gpts3(1) = -0.577350269189626;
            gwts(1)  = fact;

            gpts1(2) = 1.0/6.0;
            gpts2(2) = 4.0/6.0;
            gpts3(2) = -0.577350269189626;
            gwts(2)  = fact;

            gpts1(3) = 4.0/6.0;
            gpts2(3) = 1.0/6.0;
            gpts3(3) = -0.577350269189626;
            gwts(3)  = fact;

            gpts1(4) = 1.0/6.0;
            gpts2(4) = 1.0/6.0;
            gpts3(4) = 0.577350269189626;
            gwts(4)  = fact;

            gpts1(5) = 1.0/6.0;
            gpts2(5) = 4.0/6.0;
            gpts3(5) = 0.577350269189626;
            gwts(5)  = fact;

            gpts1(6) = 4.0/6.0;
            gpts2(6) = 1.0/6.0;
            gpts3(6) = 0.577350269189626;
            gwts(6)  = fact;

        CASE DEFAULT
            STOP " invalid value of 'ngp' in getGaussPointsPrism ! "

      END SELECT

      END SUBROUTINE getGaussPointsPrism

      END MODULE QuadratureUtilities
