! Base Class for Lagrange Elements used for
! Standard FEM

MODULE LAGRANGE_ELEMENT_BASE
IMPLICIT NONE
PRIVATE

PUBLIC :: LagrangeElement

TYPE LagrangeElement
    PRIVATE

    LOGICAL :: finite, axsy, followerLoadFlag, tracflag
    INTEGER :: elmType, matType, secType, subdomId, elmTypeNameNum
    INTEGER :: matId, finiteInt, sss, degree, npElem, ndim;
    INTEGER :: nlbf, ndof, nsize, nivGP, nGP, elenum;

!    DOUBLE PRECISION *intVar1, *intVar2, *elmDat, *matDat;

    DOUBLE PRECISION :: resi(3), primvar(3)
    INTEGER :: nodeNums(3), forAssyVec(3), globalDOFnums(3)

!    SolutionData  *SolnData;
!    GeomDataLagrange  *GeomData;

    CONTAINS

        PROCEDURE :: Constructor => Base_Constructor
    ! PROCEDURE  ::  getDimension()
    ! PROCEDURE  ::  getPolynomialDegree()
    ! PROCEDURE  ::  setSubdomainId(int sid)
    ! PROCEDURE  ::  getSubdomainId()
    ! PROCEDURE  ::  getNodesPerElement()
    ! PROCEDURE  ::  getNdofPerNode()
    ! PROCEDURE  ::  getNdofPerElement()
    ! PROCEDURE  ::  getNodeNumbers()
    ! PROCEDURE  ::  getVectorForAssembly()
    ! PROCEDURE  ::  getElmTypeNameNum()
    ! PROCEDURE  ::  printStiffnessMatrix()
    ! PROCEDURE  ::  printForceVector()
    ! PROCEDURE  ::  PrepareElemData()
    ! PROCEDURE  ::  prepareElemData2()
    ! PROCEDURE  ::  printPrimVariable()
    ! PROCEDURE  ::  getError()
    ! PROCEDURE  ::  initialiseDOFvalues()
    ! PROCEDURE  ::  initialiseKnotsAtGPs()
    ! PROCEDURE  ::  calcOutput(double u1, double v1)
    ! PROCEDURE  ::  initialiseIntVar()
    ! PROCEDURE  ::  createTractionDataVariable()
    ! PROCEDURE  ::  calcStiffnessAndResidual()
    ! PROCEDURE  ::  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
    ! PROCEDURE  ::  assembleElementMatrixAndVector(int, SparseMatrixXd&, double*)
    ! PROCEDURE  ::  assembleElementMatrix(int, Mat)
    PROCEDURE  ::  assembleElementVector
    ! PROCEDURE  ::  assembleElementVector(bool, bool, Vec, Vec, int start1=0, int start2=0);
    ! PROCEDURE  ::  volume(bool init = false)
    ! PROCEDURE  ::  calcError(int index)

END TYPE LagrangeElement

CONTAINS 
    SUBROUTINE Base_Constructor(this, finite)
      IMPLICIT NONE
      class(LagrangeElement) :: this
      INTEGER :: finite
    END SUBROUTINE Base_Constructor

    SUBROUTINE assembleElementVector(this, flag1, flag2, rhs, reac)
      IMPLICIT NONE
      class(LagrangeElement) :: this
      LOGICAL :: flag1, flag2
      DOUBLE PRECISION :: rhs, reac
    END SUBROUTINE assembleElementVector

END MODULE LAGRANGE_ELEMENT_BASE