! Module for computing basis functions and
! their derivatives for 3D Lagrange elements
!
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 25-Oct-2017
! Place : Swansea, UK
!
!
!
      MODULE  ElementUtilitiesBasisFuncs

      CONTAINS

      SUBROUTINE  LagrangeBasisFunsTria(
     1  degree,           ! degree of the polynomial
     2  xi1, xi2,         ! parametric coordinates
     3  N,                ! values of the basis functions
     4  dN_dxi1, dN_dxi2)  ! derivatives of the basis functions

      IMPLICIT NONE

      INTEGER :: degree
      DOUBLE PRECISION :: xi1, xi2, xi3
      DOUBLE PRECISION :: fact1, fact2, val1, val2, val3, val4
      DOUBLE PRECISION, DIMENSION(:) :: N, dN_dxi1, dN_dxi2

      xi3 = 1.0-xi1-xi2;
  
      SELECT CASE (degree)
        CASE (0)

          N(1) = 1.0;

          dN_dxi1(1) = 0.0;
          dN_dxi2(1) = 0.0;

        CASE (1)

          N(1) = xi3;
          N(2) = xi1;
          N(3) = xi2;

          dN_dxi1(1) = -1.0;
          dN_dxi1(2) =  1.0;
          dN_dxi1(3) =  0.0;

          dN_dxi2(1) = -1.0;
          dN_dxi2(2) =  0.0;
          dN_dxi2(3) =  1.0;

        CASE (2)

!          val1 = xi*xi;

!          N(1) = -0.5*(xi-val1);
!          N(2) =  1.0-val1;
!          N(3) =  0.5*(xi+val1);
!          N(4) = -0.5*(xi-val1);
!          N(5) =  1.0-val1;
!          N(6) =  0.5*(xi+val1);

!          val1 = 2.0*xi;
!    
!          dN_dxi(1) = -0.5*(1.0-val1);
!          dN_dxi(2) = -val1;
!          dN_dxi(3) =  0.5*(1.0+val1);
!          dN_dxi(4) = -1.0;
!          dN_dxi(5) =  1.0;
!          dN_dxi(6) =  0.0;

!          dN_dzeta(1) = -1.0;
!          dN_dzeta(2) =  0.0;
!          dN_dzeta(3) =  1.0;
!          dN_dzeta(4) = -1.0;
!          dN_dzeta(5) =  0.0;
!          dN_dzeta(6) =  1.0;

          N(1) = xi3*xi3;
          N(2) = xi1*xi1;
          N(3) = xi2*xi2;
          N(4) = 2.0*xi1*xi3;
          N(5) = 2.0*xi1*xi2;
          N(6) = 2.0*xi2*xi3;

          dN_dxi1(1) = -2.0*xi3;
          dN_dxi1(2) =  2.0*xi1;
          dN_dxi1(3) =  0.0;
          dN_dxi1(4) =  2.0*(xi3 - xi1);
          dN_dxi1(5) =  2.0*xi2;
          dN_dxi1(6) = -2.0*xi2;

          dN_dxi2(1) = -2.0*xi3;
          dN_dxi2(2) =  0.0;
          dN_dxi2(3) =  2.0*xi2;
          dN_dxi2(4) = -2.0*xi1;
          dN_dxi2(5) =  2.0*xi1;
          dN_dxi2(6) =  2.0*(xi3 -xi2);

        CASE DEFAULT
          STOP "no basis functions for degree "

      END SELECT

      END SUBROUTINE  LagrangeBasisFunsTria

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE  LagrangeBasisFunsQuad(degree, xi, eta, N, dN_dxi, dN_deta)

      IMPLICIT NONE

      INTEGER :: degree
      DOUBLE PRECISION :: xi, eta
      DOUBLE PRECISION, DIMENSION(:) :: N, dN_dxi, dN_deta

      DOUBLE PRECISION :: fact1, fact2, v1, v2, v3, v4, v5, v6

      SELECT CASE (degree)
        case (0)

          N(1) = 1.0;

          dN_dxi(1) = 0.0;
          dN_deta(1) = 0.0;

        case (1)

          v1 = 1.0 - xi;
          v2 = 1.0 + xi;
          v3 = 1.0 - eta;
          v4 = 1.0 + eta;

          N(1) = 0.25*v3*v1;
          N(2) = 0.25*v3*v2;
          N(3) = 0.25*v4*v2;
          N(4) = 0.25*v4*v1;

          dN_dxi(1) = -0.25*v3;
          dN_dxi(2) =  0.25*v3;
          dN_dxi(3) =  0.25*v4;
          dN_dxi(4) = -0.25*v4;

          dN_deta(1) = -0.25*v1;
          dN_deta(2) = -0.25*v2;
          dN_deta(3) =  0.25*v2;
          dN_deta(4) =  0.25*v1;

        CASE DEFAULT

          STOP "no basis functions defined for this degree \n"

      END SELECT

      END SUBROUTINE LagrangeBasisFunsQuad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE  computeBasisFunctions2D(
     1  flag, ETYPE, degree, 
     2  param, 
     3  xNode, yNode, 
     4  N,
     5  dN_dx, dN_dy,
     6  Jac)

      LOGICAL :: flag
      INTEGER :: ETYPE, degree
      DOUBLE PRECISION :: param(2), Jac
      DOUBLE PRECISION, DIMENSION(:) :: xNode, yNode
      DOUBLE PRECISION, DIMENSION(:) :: N, dN_dx, dN_dy

      DOUBLE PRECISION :: dN_du1(6), dN_du2(6)
      DOUBLE PRECISION :: xx, yy, detinv, B(2,2), Binv(2,2)

      INTEGER :: ii, jj, count, nlbf

      IF(ETYPE == 1) THEN !triangular elements
        IF(degree == 1) THEN
          count = 3
          nlbf = count
        ELSE
          count = 6
          nlbf = count
        END IF
      ELSE
        count = degree + 1
        nlbf = count*count
      END IF


      IF(ETYPE == 1) THEN !// triangular elements
        call LagrangeBasisFunsTria(
     1    degree,
     2    param(1), param(2),
     3    N, dN_du1, dN_du2)
      ELSE  !// quad elements
        call LagrangeBasisFunsQuad(degree, param(1), param(2), N, dN_du1, dN_du2)
      END IF

      !Gradient of mapping from parameter space to physical space
      B = 0.0
      DO ii=1, nlbf
        xx = xNode(ii)
        yy = yNode(ii)

        B(1,1) = B(1,1) + (xx * dN_du1(ii)) 
        B(2,1) = B(2,1) + (xx * dN_du2(ii)) 
        B(1,2) = B(1,2) + (yy * dN_du1(ii))
        B(2,2) = B(2,2) + (yy * dN_du2(ii)) 
      END DO

      Jac  = B(1,1)*B(2,2) - B(1,2)*B(2,1)

      detinv = 1.0/Jac

      Binv(1,1) =  B(2,2) * detinv
      Binv(1,2) = -B(1,2) * detinv
      Binv(2,1) = -B(2,1) * detinv
      Binv(2,2) =  B(1,1) * detinv

      !Compute derivatives of basis functions w.r.t physical coordinates
      DO ii=1, nlbf
        dN_dx(ii) = dN_du1(ii) * Binv(1,1) + dN_du2(ii) * Binv(1,2)
        dN_dy(ii) = dN_du1(ii) * Binv(2,1) + dN_du2(ii) * Binv(2,2)
      END DO

      END SUBROUTINE  computeBasisFunctions2D


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE LagrangeBasisFunsTet(
     1  degree, 
     2  xi1, xi2, xi3,
     3  N, dN_dxi1, dN_dxi2, dN_dxi3)

      INTEGER :: degree
      DOUBLE PRECISION :: xi1, xi2, xi3
      DOUBLE PRECISION, DIMENSION(:) :: N, dN_dxi1, dN_dxi2, dN_dxi3
      DOUBLE PRECISION :: fact1, fact2, val1, val2, val3, val4
  
       SELECT CASE (degree)
        CASE (0)

          N(1) = 1.0;

          dN_dxi1(1) = 0.0;
          dN_dxi2(1) = 0.0;
          dN_dxi3(1) = 0.0;

        CASE (1)

          N(1) = xi1;
          N(2) = xi2;
          N(3) = 1.0 - xi1 - xi2 - xi3;
          N(4) = xi3;

          dN_dxi1(1) =  1.0;
          dN_dxi1(2) =  0.0;
          dN_dxi1(4) =  0.0;
          dN_dxi1(3) = -1.0;

          dN_dxi2(1) =  0.0;
          dN_dxi2(2) =  1.0;
          dN_dxi2(4) =  0.0;
          dN_dxi2(3) = -1.0;

          dN_dxi3(1) =  0.0;
          dN_dxi3(2) =  0.0;
          dN_dxi3(4) =  1.0;
          dN_dxi3(3) = -1.0;

      CASE DEFAULT

          WRITE(*,*) "no basis functions for degree", degree

      END SELECT

      END SUBROUTINE LagrangeBasisFunsTet

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE  LagrangeBasisFunsHex(degree, xi1, xi2, xi3, N, dN_dxi1, dN_dxi2, dN_dxi3)

      IMPLICIT NONE

      INTEGER :: degree
      DOUBLE PRECISION :: xi1, xi2, xi3
      DOUBLE PRECISION, DIMENSION(:) :: N, dN_dxi1, dN_dxi2, dN_dxi3

      DOUBLE PRECISION ::  v11, v12, v21, v22, v31, v32;

      SELECT CASE (degree)
        case (0)

          N(1) = 1.0;

          dN_dxi1(1) = 0.0;
          dN_dxi2(1) = 0.0;
          dN_dxi3(1) = 0.0;

        case (1)

          v11 = 1.0 - xi1;
          v12 = 1.0 + xi1;
          v21 = 1.0 - xi2;
          v22 = 1.0 + xi2;
          v31 = 1.0 - xi3;
          v32 = 1.0 + xi3;

          N(1) = 0.125*v11*v21*v31;
          N(2) = 0.125*v12*v21*v31;
          N(3) = 0.125*v11*v22*v31;
          N(4) = 0.125*v12*v22*v31;
          N(5) = 0.125*v11*v21*v32;
          N(6) = 0.125*v12*v21*v32;
          N(7) = 0.125*v11*v22*v32;
          N(8) = 0.125*v12*v22*v32;

          dN_dxi1(1) = -0.125*v21*v31;
          dN_dxi1(2) =  0.125*v21*v31;
          dN_dxi1(3) = -0.125*v22*v31;
          dN_dxi1(4) =  0.125*v22*v31;
          dN_dxi1(5) = -0.125*v21*v32;
          dN_dxi1(6) =  0.125*v21*v32;
          dN_dxi1(7) = -0.125*v22*v32;
          dN_dxi1(8) =  0.125*v22*v32;

          dN_dxi2(1) = -0.125*v11*v31;
          dN_dxi2(2) = -0.125*v12*v31;
          dN_dxi2(3) =  0.125*v11*v31;
          dN_dxi2(4) =  0.125*v12*v31;
          dN_dxi2(5) = -0.125*v11*v32;
          dN_dxi2(6) = -0.125*v12*v32;
          dN_dxi2(7) =  0.125*v11*v32;
          dN_dxi2(8) =  0.125*v12*v32;

          dN_dxi3(1) = -0.125*v11*v21;
          dN_dxi3(2) = -0.125*v12*v21;
          dN_dxi3(3) = -0.125*v11*v22;
          dN_dxi3(4) = -0.125*v12*v22;
          dN_dxi3(5) =  0.125*v11*v21;
          dN_dxi3(6) =  0.125*v12*v21;
          dN_dxi3(7) =  0.125*v11*v22;
          dN_dxi3(8) =  0.125*v12*v22;

        CASE DEFAULT

          STOP "no basis functions defined for this degree IN LagrangeBasisFunsHex"

      END SELECT

      END SUBROUTINE LagrangeBasisFunsHex


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE LagrangeBasisFunsPrism(degree, xi1, xi2, xi4, N, dN_dxi1, dN_dxi2, dN_dxi4)
      !xi1, xi2 and xi3 are the parametric coordinates for the triangle, and
      !xi4 is the parametric coordinate in the direction normal to the triangle

      IMPLICIT NONE

      INTEGER :: degree
      DOUBLE PRECISION, INTENT(IN) :: xi1, xi2, xi4
      DOUBLE PRECISION, DIMENSION(:) :: N, dN_dxi1, dN_dxi2, dN_dxi4
      DOUBLE PRECISION :: xi3

      xi3 = 1.0- xi1 - xi2;

      SELECT CASE (degree)
        case (1)

          N(1) = xi3*(0.5*(1.0-xi4));
          N(2) = xi1*(0.5*(1.0-xi4));
          N(3) = xi2*(0.5*(1.0-xi4));
          N(4) = xi3*(0.5*(1.0+xi4));
          N(5) = xi1*(0.5*(1.0+xi4));
          N(6) = xi2*(0.5*(1.0+xi4));


          dN_dxi1(1)  = (-1.0)*(0.5*(1.0-xi4));
          dN_dxi1(2)  = ( 1.0)*(0.5*(1.0-xi4));
          dN_dxi1(3)  = ( 0.0)*(0.5*(1.0-xi4));
          dN_dxi1(4)  = (-1.0)*(0.5*(1.0+xi4));
          dN_dxi1(5)  = ( 1.0)*(0.5*(1.0+xi4));
          dN_dxi1(6)  = ( 0.0)*(0.5*(1.0+xi4));

          dN_dxi2(1)  = (-1.0)*(0.5*(1.0-xi4));
          dN_dxi2(2)  = ( 0.0)*(0.5*(1.0-xi4));
          dN_dxi2(3)  = ( 1.0)*(0.5*(1.0-xi4));
          dN_dxi2(4)  = (-1.0)*(0.5*(1.0+xi4));
          dN_dxi2(5)  = ( 0.0)*(0.5*(1.0+xi4));
          dN_dxi2(6)  = ( 1.0)*(0.5*(1.0+xi4));

          dN_dxi4(1)  = xi3*(-0.5);
          dN_dxi4(2)  = xi1*(-0.5);
          dN_dxi4(3)  = xi2*(-0.5);
          dN_dxi4(4)  = xi3*( 0.5);
          dN_dxi4(5)  = xi1*( 0.5);
          dN_dxi4(6)  = xi2*( 0.5);

        CASE DEFAULT

          STOP "no basis functions defined for this degree IN LagrangeBasisFunsHex"

      END SELECT

      END SUBROUTINE LagrangeBasisFunsPrism


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE computeBasisFunctions3D(
     1  flag, ETYPE, degree, 
     2  param, 
     3  xNode, yNode, zNode,
     4  N,
     5  dN_dx, dN_dy, dN_dz,
     6  Jac)


      LOGICAL :: flag
      INTEGER :: ETYPE, degree
      DOUBLE PRECISION :: param(3), Jac
      DOUBLE PRECISION, DIMENSION(:) :: xNode, yNode, zNode
      DOUBLE PRECISION, DIMENSION(:) :: N, dN_dx, dN_dy, dN_dz

      INTEGER :: ii=0, jj=0, kk=0, count=0, nlbf=0

      DOUBLE PRECISION :: dN_du1(8), dN_du2(8), dN_du3(8)
      DOUBLE PRECISION :: xx=0.0, yy=0.0, zz=0.0, detinv
      DOUBLE PRECISION :: B(3,3), Binv(3,3)

      IF(ETYPE == 4) THEN !tetrahedral elements
        IF(degree == 1) THEN
          count = 4;
          nlbf = count;
        ELSE IF(degree == 2) THEN
          count = 10;
          nlbf = count;
        END IF
      ELSE IF(ETYPE == 5) THEN !pyramid elements
        count = 5;
        nlbf = count;
      ELSE IF(ETYPE == 6) THEN !penta elements
        count = 6;
        nlbf = count;
      ELSE IF(ETYPE == 8) THEN !hexa elements
        count = degree + 1;
        nlbf = count*count*count;
      ELSE
        STOP "Error in 'GeomDataLagrange::computeBasisFunctions3D' "
      END IF


      IF(ETYPE == 4) THEN !tetrahedral elements
        call LagrangeBasisFunsTet(
     1    degree, 
     2    param(1), param(2), param(3),
     3    N, dN_du1, dN_du2, dN_du3)
      ELSE IF(ETYPE == 6) THEN !penta elements
        call LagrangeBasisFunsPrism(
     1    degree,
     2    param(1), param(2), param(3),
     3    N, dN_du1, dN_du2, dN_du3)
      ELSE  !! hex elements
        call LagrangeBasisFunsHex(degree, param(1), param(2), param(3), N, dN_du1, dN_du2, dN_du3)
        !Lagrange_BasisFuns1D(degree, param(1), &N1(1), &dN1(1));
        !Lagrange_BasisFuns1D(degree, param(2), &N2(1), &dN2(1));
        !Lagrange_BasisFuns1D(degree, param(3), &N3(1), &dN3(1));
      END IF

      !Gradient of mapping from parameter space to physical space
      B = 0.0

      DO ii=1, nlbf
        xx = xNode(ii)
        yy = yNode(ii)
        zz = zNode(ii)

        B(1,1) = B(1,1) + (xx * dN_du1(ii))
        B(2,1) = B(2,1) + (xx * dN_du2(ii))
        B(3,1) = B(3,1) + (xx * dN_du3(ii))

        B(1,2) = B(1,2) + (yy * dN_du1(ii))
        B(2,2) = B(2,2) + (yy * dN_du2(ii))
        B(3,2) = B(3,2) + (yy * dN_du3(ii))

        B(1,3) = B(1,3) + (zz * dN_du1(ii))
        B(2,3) = B(2,3) + (zz * dN_du2(ii))
        B(3,3) = B(3,3) + (zz * dN_du3(ii))
      END DO

      ! determinant of the matrix
      Jac = B(1,1)*(B(2,2)*B(3,3) - B(2,3)*B(3,2))
      Jac = Jac + B(1,2)*(B(2,3)*B(3,1) - B(2,1)*B(3,3))
      Jac = Jac + B(1,3)*(B(2,1)*B(3,2) - B(2,2)*B(3,1))

      ! inverse determinant of the matrix
      detinv = 1.0/Jac

      ! Calculate the inverse of the matrix
      Binv(1,1) = +detinv * (B(2,2)*B(3,3) - B(2,3)*B(3,2))
      Binv(2,1) = -detinv * (B(2,1)*B(3,3) - B(2,3)*B(3,1))
      Binv(3,1) = +detinv * (B(2,1)*B(3,2) - B(2,2)*B(3,1))
      Binv(1,2) = -detinv * (B(1,2)*B(3,3) - B(1,3)*B(3,2))
      Binv(2,2) = +detinv * (B(1,1)*B(3,3) - B(1,3)*B(3,1))
      Binv(3,2) = -detinv * (B(1,1)*B(3,2) - B(1,2)*B(3,1))
      Binv(1,3) = +detinv * (B(1,2)*B(2,3) - B(1,3)*B(2,2))
      Binv(2,3) = -detinv * (B(1,1)*B(2,3) - B(1,3)*B(2,1))
      Binv(3,3) = +detinv * (B(1,1)*B(2,2) - B(1,2)*B(2,1))


      !!Compute derivatives of basis functions w.r.t physical coordinates
      DO ii=1, nlbf
        dN_dx(ii) = dN_du1(ii) * Binv(1,1) + dN_du2(ii) * Binv(1,2) + dN_du3(ii) * Binv(1,3)
        dN_dy(ii) = dN_du1(ii) * Binv(2,1) + dN_du2(ii) * Binv(2,2) + dN_du3(ii) * Binv(2,3)
        dN_dz(ii) = dN_du1(ii) * Binv(3,1) + dN_du2(ii) * Binv(3,2) + dN_du3(ii) * Binv(3,3)
      END DO

      END SUBROUTINE computeBasisFunctions3D

      END MODULE  ElementUtilitiesBasisFuncs
