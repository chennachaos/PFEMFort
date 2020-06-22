//
// Program for paritioning the mesh
//
// This implementation reads the input mesh and
// partitions into specified number of domains
//
// Author: Dr. Chennakesava Kadapa
// Date  : 19-Oct-2017
// Place : Swansea, UK
//
//
// Inputs:
// ndim    --- number of dimensions
// eType   --- element type (numbers are as per VTK cell types)
//      5 - 3-noded Triangle
//      9 - 4-noded Quadrilateral
//     10 - 4-noded Tetrahedron
//     12 - 8-noded Hexahedron
// metisType    --- metis algorithm - nodal/dual
// node_proc_ids       --- number of partitions
// infileNodes  --- nodal coordinate file
// infileElems  --- element <-> node connectivity file
//
// Outpus:
// partition-mesh.vtk file which can be viewed using Paraview or Mayavi
//
// Note:
// Check the element<->node connectivity 
// and modify the code accordingly
//


#include "headersVTK.h"
#include "headersBasic.h"
#include "metis.h"
#include "petscmat.h"
#include "petscksp.h"


using namespace std;




int main(int argc, char* argv[])
{
    PetscErrorCode errpetsc;
  
    // initialise PETSc. This will internally initialises MPI  environment
    errpetsc = PetscInitialize(&argc, &argv, (char *)0, NULL); CHKERRQ(errpetsc);
    //errpetsc = PetscInitialize(NULL, NULL, "petsc_options.dat", NULL); CHKERRQ(errpetsc);
  
    //////////////////////////////////////////////
    //
    // declare vtk variables
    //
    //////////////////////////////////////////////


    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkLine>                 lineVTK      =  vtkSmartPointer<vtkLine>::New();
    vtkSmartPointer<vtkQuad>                 quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkHexahedron>           hexVTK       =  vtkSmartPointer<vtkHexahedron>::New();
    vtkSmartPointer<vtkTriangle>             triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkTetra>                tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkPyramid>              pyramidVTK   =  vtkSmartPointer<vtkPyramid>::New();
    vtkSmartPointer<vtkWedge>                wedgeVTK     =  vtkSmartPointer<vtkWedge>::New();
    vtkSmartPointer<vtkIntArray>             nodeInOutVTK =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray>             procIdVTK    =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray>      cellOrientationVTK  =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkFloatArray>           vecVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           vecVTK2      =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           distVTK      =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTK  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTK2 =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    uGridVTK->Reset();
    pointsVTK->Reset();
    cellDataVTK->Reset();
    cellDataVTK->Reset();
    procIdVTK->Reset();
    nodeInOutVTK->Reset();

    int  ndim, eType, metisType, nparts, ncommon_nodes;
    int  nNode, nElem, npElem, count, ind;
    int  ee, ii, jj, kk, e1, e2, ind1, ind2;
    string nodeFileName, elemFileName;


    if(argc == 0)
      cerr << " Error in input data " << endl;
    else
    {
       ndim      = atoi(argv[1]);
       eType     = atoi(argv[2]);
       metisType = atoi(argv[3]);
       nparts    = atoi(argv[4]);

       nodeFileName = argv[5];
       elemFileName = argv[6];
    }


    //Check element type and set the parameters
    if(eType == 5)
    {
      npElem = 3;
      ncommon_nodes = 2;
    }
    else if(eType == 9)
    {
      npElem = 4;
      ncommon_nodes = 2;
    }
    else if(eType == 10)
    {
      npElem = 4;
      ncommon_nodes = 3;
    }
    else if(eType == 12)
    {
      npElem = 8;
      ncommon_nodes = 4;
    }
    else
    {
      cerr << " Wrong element type" << endl;
    }

    std::ifstream  infile_nodes(nodeFileName);
    std::ifstream  infile_elems(elemFileName);


    if(infile_nodes.fail())
    {
       cout << " Could not open the input nodes file " << endl;
       exit(1);
    }

    double  val[10];

    std::string line;
    
    // read nodal coordinates
    ////////////////////////////////////////////

    cout << " reading nodes " << endl;

    nNode = 0;
    while (std::getline(infile_nodes, line))
      ++nNode;

    vector<vector<double> >  node_coords(nNode, vector<double>(3));

    infile_nodes.clear();
    infile_nodes.seekg(0, infile_nodes.beg);

    if(ndim == 2)
    {
      ii=0;
      while(infile_nodes >> val[0] >> val[1] >> val[2] )
      {
        //printf("%12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2]);

        node_coords[ii][0] = val[1];
        node_coords[ii][1] = val[2];
        node_coords[ii][2] = 0.0;

        ii++;
      }
    }
    else
    {
      ii=0;
      while(infile_nodes >> val[0] >> val[1] >> val[2] >> val[3] )
      {
        //printf("%12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2]);

        node_coords[ii][0] = val[1];
        node_coords[ii][1] = val[2];
        node_coords[ii][2] = val[3];

        ii++;
      }
    }


    // read elements
    ////////////////////////////////////////////
    cout << " reading elements " << endl;

    if(infile_elems.fail())
    {
       cout << " Could not open the input elements file " << endl;
       exit(1);
    }


    nElem = 0;
    while (std::getline(infile_elems, line))
      ++nElem;

    cout << " nElem   " << nElem << endl;

    PetscInt  val2[10];

    infile_elems.clear();
    infile_elems.seekg(0, infile_elems.beg);
    
    vector<vector<int> >  elemNodeConn(nElem, vector<int>(npElem));

    if(npElem == 3)
    {
      ee=0;
      while(infile_elems >> val2[0] >> val2[1] >> val2[2] >> val2[3] )
      {
        //printf("%6d \t %6d \t %6d \t %6d \n", val2[4], val2[5], val2[6], val2[7]);

        for(ii=1; ii<=npElem; ii++)
          elemNodeConn[ee][ii-1] = val2[ii]-1;

        ee++;
      }
    }
    else if(npElem == 4)
    {
      ee=0;
      while(infile_elems >> val2[0] >> val2[1] >> val2[2] >> val2[3] >> val2[4] )
      {
        //printf("%6d \t %6d \t %6d \t %6d \n", val2[4], val2[5], val2[6], val2[7]);

        for(ii=1; ii<=npElem; ii++)
          elemNodeConn[ee][ii] = val2[ii]-1;

        ee++;
      }
    }
    else if(npElem == 8)
    {
      ee=0;
      while(infile_elems >> val2[0] >> val2[1] >> val2[2] >> val2[3] >> val2[4] >> val2[5] >> val2[6] >> val2[7] >> val2[8] )
      {
        //printf("%6d \t %6d \t %6d \t %6d \n", val2[4], val2[5], val2[6], val2[7]);

        for(ii=1; ii<=npElem; ii++)
          elemNodeConn[ee][ii] = val2[ii]-1;

        ee++;
      }
    }


    infile_nodes.close();
    infile_elems.close();

    cout << " nElem    " << nElem << endl;
    cout << " nNode    " << nNode << endl;
    cout << " npElem   " << npElem << endl;

    procIdVTK->SetName("procIde");
    procIdVTK->SetNumberOfTuples(nElem);


      cout << " partitioning the mesh " << endl;

      //
      /////////////////////////////////////////////////////////////////////////////
      //
      // Partition the mesh. This can be done using software libraries
      // Chaco, Jostle, METIS and Scotch.
      //
      // Here, METIS partitioning routines are used.
      // 
      /////////////////////////////////////////////////////////////////////////////

      idx_t nWeights  = 1;
      idx_t objval;
      idx_t *xadj, *adjncy;
      idx_t numflag=0;

      PetscInt  *eptr, *eind, *elem_proc_id, *node_proc_id;

      errpetsc  = PetscMalloc1(nElem+1,       &eptr);  CHKERRQ(errpetsc);
      errpetsc  = PetscMalloc1(nElem,         &elem_proc_id); CHKERRQ(errpetsc);
      errpetsc  = PetscMalloc1(nNode,         &node_proc_id); CHKERRQ(errpetsc);


      idx_t npElem_total = nElem*npElem;

      cout << " npElem_total = " << npElem_total << endl;

      errpetsc  = PetscMalloc1(npElem_total,  &eind); CHKERRQ(errpetsc);

      vector<int>  vecTemp2(npElem);

      eptr[0] = 0;

      kk=0;
      for(ee=0; ee<nElem; ee++)
      {
        eptr[ee+1] =  (ee+1)*npElem;

        vecTemp2 = elemNodeConn[ee] ;

        for(ii=0; ii<npElem; ii++)
          eind[kk+ii] = vecTemp2[ii] ;

        kk += npElem;
      }

      // for(ee=0; ee<nElem; ee++)
      // {
      //   ii = ee*npElem;
      //   cout << eind[ii] << '\t' << eind[ii+1] << '\t' << eind[ii+2] << endl;
      // }
      // cout << " \n\n\n " << endl;

      idx_t options[METIS_NOPTIONS];

      METIS_SetDefaultOptions(options);

      // Specifies the partitioning method.
      //options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;    // Multilevel recursive bisectioning.
      options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;  // Multilevel k-way partitioning.

      //options[METIS_OPTION_NSEPS] = 10;

      options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;   // Edge-cut minimization
      //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; // Total communication volume minimization

      options[METIS_OPTION_NUMBERING] = 0;  // C-style numbering is assumed that starts from 0.

      // METIS partition routine
      int ret;
      if(metisType == 1)
        ret = METIS_PartMeshNodal(&nElem, &nNode, eptr, eind, NULL, NULL, &nparts, NULL, options, &objval, elem_proc_id, node_proc_id);
      else
        ret = METIS_PartMeshDual(&nElem, &nNode, eptr, eind, NULL, NULL, &ncommon_nodes, &nparts, NULL, options, &objval, elem_proc_id, node_proc_id);

      if(ret == METIS_OK)
        cout << " METIS partition routine successful "  << endl;
      else
        cout << " METIS partition routine FAILED "  << endl;

      for(ee=0; ee<nNode; ee++)
        cout << ee << '\t' << node_proc_id[ee] << endl;

      for(ee=0; ee<nElem; ee++)
      {
        cout << ee << '\t' << elem_proc_id[ee] << endl;
        procIdVTK->InsertTuple1(ee, elem_proc_id[ee]);
      }

      errpetsc  = PetscFree(eptr);  CHKERRQ(errpetsc);
      errpetsc  = PetscFree(eind);  CHKERRQ(errpetsc);
      errpetsc  = PetscFree(elem_proc_id); CHKERRQ(errpetsc);
      errpetsc  = PetscFree(node_proc_id); CHKERRQ(errpetsc);

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    //
    // setup and write polyData
    //
    //////////////////////////////////////////////

    vtkIdType pt[10];

    cout << " Writing VTK file  " << endl;

    for(ii=0; ii<nNode; ii++)
      pt[0] = pointsVTK->InsertNextPoint(node_coords[ii][0], node_coords[ii][1], node_coords[ii][2]);

    if(ndim == 2)
    {
      if(npElem == 3)
      {
        for(ee=0; ee<nElem; ee++)
        {
          for(ii=0; ii<npElem; ii++)
            triaVTK->GetPointIds()->SetId(ii, elemNodeConn[ee][ii] );
  
          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
      }
      else if(npElem == 4)
      {
        for(ee=0; ee<nElem; ee++)
        {
          for(ii=0; ii<npElem; ii++)
            quadVTK->GetPointIds()->SetId(ii, elemNodeConn[ee][ii] );
  
          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
      }
    }
    else
    {
      if(npElem == 4)
      {
        for(ee=0; ee<nElem; ee++)
        {
          for(ii=0; ii<npElem; ii++)
          tetraVTK->GetPointIds()->SetId(ii, elemNodeConn[ee][ii] );
  
          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());
        }
      }
      else if(npElem == 8)
      {
        for(ee=0; ee<nElem; ee++)
        {
          for(ii=0; ii<npElem; ii++)
            hexVTK->GetPointIds()->SetId(ii, elemNodeConn[ee][ii] );
  
          uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
        }
      }
    }

      char fname1[100];

      //if(elType == QUAD)
        //sprintf(fname1,"%s%d%s%d%s%d%s", "mesh-quad-",nEx,"-",nEy,"-",n_mpi_procs,".vtp");
      //else
        //sprintf(fname1,"%s%d%s%d%s%d%s", "mesh-tria-",nEx,"-",nEy,"-",n_mpi_procs,".vtp");

      sprintf(fname1,"%s", "mesh-partition-cpp.vtu");
      //sprintf(fname1,"%s%d%s%d%s%d%s%d%s", "mesh-tria-",nEx,"-",nEy,"-",n_mpi_procs,"-",this_mpi_proc,".vtp");

      uGridVTK->SetPoints(pointsVTK);
      uGridVTK->GetCellData()->AddArray(procIdVTK);

      //Write the file.
      writerUGridVTK->SetFileName(fname1);
      writerUGridVTK->SetInputData(uGridVTK);
      writerUGridVTK->Write();

    // finalise  PETSc. This will internally finalises MPI  environment       //MPI_Finalize();
    errpetsc = PetscFinalize();CHKERRQ(errpetsc);

    return 0;
}












