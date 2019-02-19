MODULE LAGRANGE_ELEMENT_2D_NAVIERSTOKES_TRIA3NODE
IMPLICIT NONE
PUBLIC :: LagrangeElem2DNavierStokesTria3Node

TYPE LagrangeElem2DNavierStokesTria3Node 
{
  public:

    LagrangeElem2DNavierStokesTria3Node();

    virtual ~LagrangeElem2DNavierStokesTria3Node();

    virtual int getElmTypeNameNum()
    {  return 13; }

    virtual int calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal);

    virtual void prepareElemData();

    void prepareElemData2();


    virtual int calcInternalForces();

    virtual int calcLoadVector();

    virtual int calcOutput(double u1, double v1);
    
    virtual void discreteContourplot(int, int, int, int, double, double);

    virtual void projectToKnots(bool, int, int, int);

    virtual void projectStrain(int, int, double*);

    virtual void projectStress(int, double*);

    virtual void projectIntVar(int, double*);

    virtual void toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&);

END TYPE LagrangeElem2DNavierStokesTria3Node 



END MODULE LAGRANGE_ELEMENT_2DNAVIERSTOKES_TRIA3NODE