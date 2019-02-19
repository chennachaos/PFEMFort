/*=========================================================================

Create triangulation for the base grid and perform 
subtriangulation for the cut-cells


=========================================================================*/

#include "headersVTK.h"
#include "headersBasic.h"
#include <algorithm>


using namespace std;


//typedef float REAL;

typedef double REAL;



int main(int argc, char* argv[])
{
    //////////////////////////////////////////////
    //
    // declare vtk variables
    //
    //////////////////////////////////////////////


    time_t tstart, tend;

    //tstart = time(0);

        vtkSmartPointer<vtkDataSetMapper>        mapperVTK;
        vtkSmartPointer<vtkActor>                actorVTK;
        vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK;
        vtkSmartPointer<vtkPoints>               pointsVTK;
        vtkSmartPointer<vtkPoints>               pointsVTK2;
        vtkSmartPointer<vtkVertex>               vertexVTK;
        vtkSmartPointer<vtkLine>                 lineVTK;
        vtkSmartPointer<vtkQuad>                 quadVTK;
        vtkSmartPointer<vtkHexahedron>           hexVTK;

        vtkSmartPointer<vtkTriangle>             triaVTK;
        vtkSmartPointer<vtkPolygon>              polygonVTK;
        vtkSmartPointer<vtkTetra>                tetraVTK;
        vtkSmartPointer<vtkPyramid>              pyramidVTK;
        vtkSmartPointer<vtkWedge>                wedgeVTK;


        vtkSmartPointer<vtkIntArray>          nodeInOutVTK;
        vtkSmartPointer<vtkIntArray>          cutCellTypeVTK;
        vtkSmartPointer<vtkIntArray>          cellOrientationVTK;
        //vtkSmartPointer<vtkDoubleArray>          vectors, vectors2, scalars, scalars2;
        vtkSmartPointer<vtkFloatArray>          vecVTK, vecVTK2, distVTK, scaVTK2, cellDataVTK, cellDataVTK2;
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK;

        vtkSmartPointer<vtkExtractEdges>     extractEdgesVTK;

    mapperVTK    =  vtkSmartPointer<vtkDataSetMapper>::New();
    actorVTK     =  vtkSmartPointer<vtkActor>::New();
    uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New(); 
    pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    pointsVTK2   =  vtkSmartPointer<vtkPoints>::New();
    lineVTK      =  vtkSmartPointer<vtkLine>::New();
    quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    hexVTK       =  vtkSmartPointer<vtkHexahedron>::New();
    vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    polygonVTK   =  vtkSmartPointer<vtkPolygon>::New();
    tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    pyramidVTK   =  vtkSmartPointer<vtkPyramid>::New();
    wedgeVTK     =  vtkSmartPointer<vtkWedge>::New();

    vecVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    distVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    scaVTK2      =  vtkSmartPointer<vtkFloatArray>::New();
    vecVTK2      =  vtkSmartPointer<vtkFloatArray>::New();
    cellDataVTK  =  vtkSmartPointer<vtkFloatArray>::New();
    cellDataVTK2 =  vtkSmartPointer<vtkFloatArray>::New();
  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    nodeInOutVTK       =  vtkSmartPointer<vtkIntArray>::New();
    cutCellTypeVTK       =  vtkSmartPointer<vtkIntArray>::New();
  cellOrientationVTK    =  vtkSmartPointer<vtkIntArray>::New();

     extractEdgesVTK  =    vtkSmartPointer<vtkExtractEdges>::New();


  vtkSmartPointer<vtkUnstructuredGrid> uGrid =   vtkSmartPointer<vtkUnstructuredGrid>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writeruGrid   =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writeruGrid2  =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();


    //////////////////////////////////////////////
    //
    // generate base triangulation
    //
    //////////////////////////////////////////////

    REAL  x0, x1, y0, y1, z0, z1;

    int nEx, nEy, nEz, nNx, nNy, nNz, nNode, nElem;

    x0 = -1.6;    x1 =  1.6;    nEx = 5;
    y0 = -1.6;    y1 =  1.6;    nEy = 5;
    z0 = -1.6;    z1 =  1.6;    nEz = 5;

//
    if(argc == 0)
      cerr << " Input data " << endl;
    else
    {
       x0  = atof(argv[1]);
       x1  = atof(argv[2]);
       nEx = atoi(argv[3]);

       y0  = atof(argv[4]);
       y1  = atof(argv[5]);
       nEy = atoi(argv[6]);

       z0  = atof(argv[7]);
       z1  = atof(argv[8]);
       nEz = atoi(argv[9]);
    }
//


    nElem = nEx*nEy*nEz;

    nNx = nEx+1;
    nNy = nEy+1;
    nNz = nEz+1;

    nNode = nNx*nNy*nNz;


    uGridVTK->Reset();
    pointsVTK->Reset();
    cellDataVTK->Reset();
    cellDataVTK->Reset();
    cutCellTypeVTK->Reset();
    nodeInOutVTK->Reset();

    int  ii, jj, kk, ll, n1, n2, n3, n4, nn;
    int  ind, ind1, ind2, ind3, ind4, ind5, ind6;

    REAL dx = (x1-x0)/nEx;
    REAL dy = (y1-y0)/nEy;
    REAL dz = (z1-z0)/nEz;
    REAL  xx, yy, zz, fact;

    cout << x0 << '\t' << x1 << endl;
    cout << y0 << '\t' << y1 << endl;
    cout << z0 << '\t' << z1 << endl;


    vtkIdType pts[10];

    distVTK->SetName("dist");
    distVTK->SetNumberOfTuples(nNode);
    nodeInOutVTK->SetName("InOut");
    nodeInOutVTK->SetNumberOfTuples(nNode);
    cutCellTypeVTK->SetName("cellType");
    cutCellTypeVTK->SetNumberOfTuples(nElem);
    cellOrientationVTK->SetName("orient");
    cellOrientationVTK->SetNumberOfTuples(nElem);

    //////////////////////////////////////////////
    //
    // create nodes/points
    //
    //////////////////////////////////////////////

      ofstream fout("mesh-nodes.dat");

      if(fout.fail())
      {
         cout << " Could not open the Output file" << endl;
      exit(1);
      }

      fout.setf(ios::fixed);
      fout.setf(ios::showpoint);
      fout.precision(8);

    //fout << "\n\n\n" << endl;
    //fout << "nodes \n" << endl;

    ind = 1;
    zz = z0;
    for(kk=0; kk<nNz; kk++)
    {
      yy = y0;
      for(jj=0;jj<nNy;jj++)
      {
        xx = x0;
        for(ii=0;ii<nNx;ii++)
        {
          pts[0] = pointsVTK->InsertNextPoint(xx, yy, zz);

          distVTK->InsertTuple1(ind-1, 0.0);

          fout << ind << '\t' << xx << '\t' << yy << '\t' << zz << endl;

          xx += dx;
          ind++;
        }
        yy += dy;
      }
      zz += dz;
    }

    fout.close();

    //////////////////////////////////////////////
    //
    // create cells/elements
    //
    //////////////////////////////////////////////

    vtkIdType  cellId, ptId;

    cout << " AAAAAAAAAAAAA " << endl;
    cout << " pointsVTK->GetNumberOfPoints() = " <<  pointsVTK->GetNumberOfPoints() << endl;

    //fout << "\n\n\n" << endl;
    //fout << "elements \n" << endl;

      ofstream fout2("mesh-elems.dat");

      if(fout2.fail())
      {
         cout << " Could not open the Output file" << endl;
         exit(1);
      }

      fout2.setf(ios::fixed);
      fout2.setf(ios::showpoint);
      fout2.precision(8);


    nn = nNx*nNy;
    cellId = 1;
    ind=0;
    for(kk=0; kk<nEz; kk++)
    {
      ind5 = nn*kk;
      ind6 = nn*(kk+1);

      for(jj=0; jj<nEy; jj++)
      {
        ind1 = ind5 + nNx*jj;
        ind2 = ind5 + nNx*(jj+1);

        ind3 = ind6 + nNx*jj;
        ind4 = ind6 + nNx*(jj+1);

        for(ii=0;ii<nEx;ii++)
        {
          pts[0] = ind1+ii;          pts[4] = ind3+ii;
          pts[1] = pts[0]+1;         pts[5] = pts[4]+1;
          pts[2] = ind2+ii;          pts[6] = ind4+ii;
          pts[3] = pts[2]+1;         pts[7] = pts[6]+1;

          tetraVTK->GetPointIds()->SetId(0, pts[0]);
          tetraVTK->GetPointIds()->SetId(1, pts[1]);
          tetraVTK->GetPointIds()->SetId(2, pts[3]);
          tetraVTK->GetPointIds()->SetId(3, pts[5]);

          uGrid->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());

          tetraVTK->GetPointIds()->SetId(0, pts[0]);
          tetraVTK->GetPointIds()->SetId(1, pts[3]);
          tetraVTK->GetPointIds()->SetId(2, pts[2]);
          tetraVTK->GetPointIds()->SetId(3, pts[5]);

          uGrid->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());

          tetraVTK->GetPointIds()->SetId(0, pts[2]);
          tetraVTK->GetPointIds()->SetId(1, pts[3]);
          tetraVTK->GetPointIds()->SetId(2, pts[7]);
          tetraVTK->GetPointIds()->SetId(3, pts[5]);

          uGrid->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());

          tetraVTK->GetPointIds()->SetId(0, pts[4]);
          tetraVTK->GetPointIds()->SetId(1, pts[6]);
          tetraVTK->GetPointIds()->SetId(2, pts[7]);
          tetraVTK->GetPointIds()->SetId(3, pts[2]);

          uGrid->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());

          tetraVTK->GetPointIds()->SetId(0, pts[4]);
          tetraVTK->GetPointIds()->SetId(1, pts[7]);
          tetraVTK->GetPointIds()->SetId(2, pts[5]);
          tetraVTK->GetPointIds()->SetId(3, pts[2]);

          uGrid->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());

          tetraVTK->GetPointIds()->SetId(0, pts[0]);
          tetraVTK->GetPointIds()->SetId(1, pts[4]);
          tetraVTK->GetPointIds()->SetId(2, pts[5]);
          tetraVTK->GetPointIds()->SetId(3, pts[2]);

          uGrid->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());

          for(ll=0;ll<8;ll++)
          {
            pts[ll] += 1;
          }

          fout2 << cellId++ << '\t' << pts[0] << '\t' << pts[1] << '\t' << pts[3] << '\t' << pts[5] << endl;
          fout2 << cellId++ << '\t' << pts[0] << '\t' << pts[3] << '\t' << pts[2] << '\t' << pts[5] << endl;
          fout2 << cellId++ << '\t' << pts[2] << '\t' << pts[3] << '\t' << pts[7] << '\t' << pts[5] << endl;
          fout2 << cellId++ << '\t' << pts[4] << '\t' << pts[6] << '\t' << pts[7] << '\t' << pts[2] << endl;
          fout2 << cellId++ << '\t' << pts[4] << '\t' << pts[7] << '\t' << pts[5] << '\t' << pts[2] << endl;
          fout2 << cellId++ << '\t' << pts[0] << '\t' << pts[4] << '\t' << pts[5] << '\t' << pts[2] << endl;

          //fout2 << cellId++ << "\t 1 \t 1 \t 1 \t" << pts[0] << '\t' << pts[1] << '\t' << pts[3] << '\t' << pts[5] << endl;
          //fout2 << cellId++ << "\t 1 \t 1 \t 1 \t" << pts[0] << '\t' << pts[3] << '\t' << pts[2] << '\t' << pts[5] << endl;
          //fout2 << cellId++ << "\t 1 \t 1 \t 1 \t" << pts[2] << '\t' << pts[3] << '\t' << pts[7] << '\t' << pts[5] << endl;
          //fout2 << cellId++ << "\t 1 \t 1 \t 1 \t" << pts[4] << '\t' << pts[6] << '\t' << pts[7] << '\t' << pts[2] << endl;
          //fout2 << cellId++ << "\t 1 \t 1 \t 1 \t" << pts[4] << '\t' << pts[7] << '\t' << pts[5] << '\t' << pts[2] << endl;
          //fout2 << cellId++ << "\t 1 \t 1 \t 1 \t" << pts[0] << '\t' << pts[4] << '\t' << pts[5] << '\t' << pts[2] << endl;

          //cellId++;
        }
      }
    }

    fout2.close();

    ////////////////////////////////////////////////////////////////////
    //
    // boundary conditions
    //
    ////////////////////////////////////////////////////////////////////


    //fout << "\n\n\n" << endl;
    //fout << "prescribed boundary conditions \n" << endl;

    vector<int>  boundNodes;


    int nnXY = nNx*nNy;
    int nnYZ = nNy*nNz;
    int nnXZ = nNx*nNz;

    int side1=1, side2=1, side3=1, side4=1, side5=1, side6=1, ndof=1;

    /*
    if(side == 0) // normal to X, left side
    {
      ind=0;
      for(kk=0; kk<nNz; kk++)
      {
        ind5 = nn*kk;
        ind6 = nn*(kk+1);

        for(jj=0; jj<1; jj++)
        {
          ind1 = ind5 + nNx*jj;
          ind2 = ind5 + nNx*(jj+1);

          ind3 = ind6 + nNx*jj;
          ind4 = ind6 + nNx*(jj+1);

          for(ii=0;ii<nNx;ii++)
          {
            ptId = ind1+ii+1;

            for(ll=1;ll<=ndof;ll++)
            {
              fout << ptId << '\t' << ll << '\t' << 0.0 << endl;
            }

            cellId++;
          }
        }
      }
    }
    */


    if(side1 == 1) // normal to X, left side
    {
      for(kk=0; kk<nNz; kk++)
      {
        ind5 = nnXY*kk;

        for(jj=0; jj<nNy; jj++)
        {
          ind1 = ind5 + nNx*jj;

          for(ii=0; ii<1; ii++)
          {
            boundNodes.push_back(ind1+ii+1);
          }
        }
      }
    }
    if(side2 == 1) // normal to X, right side
    {
      for(kk=0; kk<nNz; kk++)
      {
        ind5 = nnXY*kk;

        for(jj=0; jj<nNx; jj++)
        {
          ind1 = ind5 + nNx*jj;

          for(ii=nNx-1; ii<nNx; ii++)
          {
            boundNodes.push_back(ind1+ii+1);
          }
        }
      }
    }
    if(side3 == 1) // normal to Y, bottom side
    {
      for(kk=0; kk<nNz; kk++)
      {
        ind5 = nnXY*kk;

        for(jj=0; jj<1; jj++)
        {
          ind1 = ind5 + nNx*jj;

          for(ii=0; ii<nNx; ii++)
          {
            boundNodes.push_back(ind1+ii+1);
          }
        }
      }
    }

    if(side4 == 1) // normal to Y, top side
    {
      ind=0;
      for(kk=0; kk<nNz; kk++)
      {
        ind5 = nnXY*kk;
        ind6 = nnXY*(kk+1);

        for(jj=nNy-1; jj<nNy; jj++)
        {
          ind1 = ind5 + nNx*jj;

          for(ii=0; ii<nNx; ii++)
          {
            boundNodes.push_back(ind1+ii+1);
          }
        }
      }
    }
    if(side5 == 1) // normal to Z, back side
    {
      for(kk=0; kk<1; kk++)
      {
        ind5 = nnXY*kk;

        for(jj=0; jj<nNy; jj++)
        {
          ind1 = ind5 + nNx*jj;

          for(ii=0;ii<nNx;ii++)
          {
            boundNodes.push_back(ind1+ii+1);
          }
        }
      }
    }
    if(side6 == 1) // normal to Z, front side
    {
      for(kk=nNz-1; kk<nNz; kk++)
      {
        ind5 = nnXY*kk;

        for(jj=0; jj<nNy; jj++)
        {
          ind1 = ind5 + nNx*jj;
          for(ii=0;ii<nNx;ii++)
          {
            boundNodes.push_back(ind1+ii+1);
          }
        }
      }
    }


      ofstream fout3("mesh-DirichBC.dat");

      if(fout3.fail())
      {
         cout << " Could not open the Output file" << endl;
         exit(1);
      }

      fout3.setf(ios::fixed);
      fout3.setf(ios::showpoint);
      fout3.precision(8);


    double  coord[3], val;

    sort(boundNodes.begin(), boundNodes.end());
    boundNodes.erase(unique(boundNodes.begin(), boundNodes.end()), boundNodes.end());

    for(ii=0; ii<boundNodes.size(); ii++)
    {
      n1 = boundNodes[ii]-1;
      pointsVTK->GetPoint(n1, coord);

      val = coord[0]*coord[0]+coord[1]*coord[1]+coord[2]*coord[2];

      distVTK->InsertTuple1(n1, val);

      fout3 << (n1+1) << '\t' << 1 << '\t' << val << endl;
    }

    fout3.close();

    //////////////////////////////////////////////
    //
    // setup and write polyData
    //
    //////////////////////////////////////////////

    char fname1[100];
    sprintf(fname1,"%s", "Tetmesh.vtu");

    uGrid->SetPoints(pointsVTK);

    uGrid->GetPointData()->AddArray(distVTK);

    writeruGrid->SetFileName(fname1);
    writeruGrid->SetInputData(uGrid);
    writeruGrid->Write();

    cout << " uGrid->GetNumberOfPoints() = " <<  uGrid->GetNumberOfPoints() << endl;
    cout << " uGrid->GetNumberOfCells()  = " <<  uGrid->GetNumberOfCells() << endl;


    return 0;
}










