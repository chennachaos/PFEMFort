/*=========================================================================

Create triangulation for the base grid and perform 
subtriangulation for the cut-cells


=========================================================================*/

#include <iostream>
#include <limits.h>
#include <float.h>
#include <vector>
#include <assert.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <set>
#include <fstream>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>

#include <numeric>
#include <iterator>
#include <iostream>
#include <functional>



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

    int  ii, jj, kk, ll, n1, n2, n3, n4, nn;
    int  ind, ind1, ind2, ind3, ind4, ind5, ind6;

    REAL dx = (x1-x0)/nEx;
    REAL dy = (y1-y0)/nEy;
    REAL dz = (z1-z0)/nEz;
    REAL  xx, yy, zz, fact;

    cout << x0 << '\t' << x1 << endl;
    cout << y0 << '\t' << y1 << endl;
    cout << z0 << '\t' << z1 << endl;

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

    cout << " AAAAAAAAAAAAA " << endl;

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
    int cellId = 1, pts[20];
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

          for(ll=0;ll<8;ll++)
           pts[ll] += 1;


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

    nn = nNx*nNy;

    int side1=1, side2=1, side3=1, side4=1, side5=1, side6=1, ndof=1;

    //int side1=0, side2=0, side3=1, side4=0, side5=0, side6=0, ndof=3;

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
        ind5 = nn*kk;

        for(jj=0; jj<nNx; jj++)
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
        ind5 = nn*kk;

        for(jj=0; jj<nNy; jj++)
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
        ind5 = nn*kk;

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
        ind5 = nn*kk;
        ind6 = nn*(kk+1);

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
        ind5 = nn*kk;

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
        ind5 = nn*kk;

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

      
      if(ndof == 1)
        val = coord[0]*coord[0]+coord[1]*coord[1]+coord[2]*coord[2];
      else
        val = 0.0;

      for(jj=1; jj<=ndof; jj++)
        fout3 << (n1+1) << '\t' << jj << '\t' << val << endl;
    }

    fout3.close();

    return 0;
}










