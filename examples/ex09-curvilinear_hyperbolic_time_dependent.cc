/*!
\file ex09-curvilinear_hyperbolic_time_dependent.cc

\brief Example of an hyperbolic equation on curvilinear grids.

This file performs the solution of an hyperbolic equation on curvilinear grids.

\date: Sunday, September 09, 2012
\version: 2012-09-03.

\author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
\author: Vincent Verardi:
*/
/*
Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
State University (SDSU).

http://www.csrc.sdsu.edu/mtk/

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification are permitted provided that the following conditions are met:

-# Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
-# Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
-# Neither the name of the CSRC, SDSU nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.
-# Modifications whether they are partial or complete; whether they are
additions or eliminations should be reported through an email at:
esanchez@sciences.sdsu.edu

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <locale>
#include <sstream>

#include "mtk.h"

using namespace std;

void CreateSnapshotInTime (int current_iteration,
                           double transient_problem,
                           MTKGrid1D *q,
                           MTKGrid1D *qinit) {

  char datfilename[10];
  char gplfilename[10];
  char pic[100];
  char tmpstr[100];
  char tmp2[100];
  FILE *datfile;
  int i;
  int jj;

  i = current_iteration;

  if ((i + 1) % 100 == 0) {

    // Create snapshot data file:
    sprintf(datfilename, "%d.dat", (i + 1));

    printf("Snapshot created: %s\n", datfilename);
    datfile = fopen(datfilename, "w");
    if (datfile == NULL) {
      fprintf(stderr, "Error creating file\n");
    }

    for (jj = 0; jj < q->nNodes() + 2*q->nDummy(); jj++) {
      if (jj < q->nDummy()) {
        fprintf(datfile, "%f\t%f\t%f\n", q->GetGridValue(jj), 0.0, 0.0);
      } else {
        if (jj > q->nNodes() - 1) {
          fprintf(datfile, "%f\t%f\t%f\n", q->GetGridValue(jj), 0.0, 0.0);
        } else {
          fprintf(datfile, "%f\t%f\t%f\n", q->GetGridValue(jj),
                  qinit->GetDataValue(jj + q->nDummy()),
                  q->GetDataValue(jj + q->nDummy()));
        }
      }
    }
    fclose(datfile);

    sprintf(gplfilename, "%d.gpl", (i + 1));
    datfile = fopen(gplfilename, "w");
    fprintf(datfile, "set terminal png\n");
    sprintf(pic, "set output 'snsp_%d.png'\n", i + 1);
    fprintf(datfile, "%s", pic);
    sprintf(pic, "set title 'Advection Equation with Cr = %f'\n",
            transient_problem);
    fprintf(datfile, "%s", pic);
    fprintf(datfile, "set grid\n");
    fprintf(datfile, "set yrange [0:1]\n");
    fprintf(datfile, "set key center outside bottom title 'Legend' box 3\n");
    sprintf(tmpstr,
      "plot '%d.dat' u 1:2 title 'Initial', '%d.dat' u 1:3 title 'Current'\n",
      i + 1,
      i + 1);
    fprintf(datfile, "%s", tmpstr);
    fclose(datfile);
    sprintf(tmp2, "gnuplot %d.gpl", i + 1);
    system(tmp2);
  }
}

void GenerateAndSeeAnimation () {

  system("./make_video.pl");
  system("gwenview anim-ex09.gif &");
}

int main() {

  //! Variables:
  double Cr;
  double dx;
  double dt;
  double u;
  double xmin;
  double xmax;
  // Number of iterations required to pass through the entire domain k times:
  double k = 2.0;
  double xc;
  double z;
  double Step;
  FILE *datfile;
  int i;
  int nDummy;
  int numItt;
  int cont;

  //! Begin execution:
  Cr = 1.0;
  dx = 0.002;
  dt = 0.02;
  u = Cr*dx/dt;
  xmin = 0.0;
  xmax = 1.0;
  nDummy = 2;
  numItt = (int)((xmax - xmin)*k/fabs(u)/dt);

  MTKGrid1D *q = new MTKGrid1D(nDummy, dx, xmin, xmax);

  q->SetGridValues();

  BindingType1D CB_NoBC = CenterBound_NoBC;

  q->ComputeMetric(CB_NoBC);

  BindingType1D btype = CenterBound;

  if (!q->BindDataonGrid1D(btype)) {
    fprintf(stderr, "ERROR: Cannot Bind Data!\n");
    return EXIT_FAILURE;
  }

  // Establish initial condition:
  for(i = 0; i < q->nNodes(); i++) {
    xc = 0.5*(q->GetGridValue(i + q->nDummy()) +
          q->GetGridValue(i + 1 +q->nDummy()));
    z = abs(xc - 0.5);
    q->SetDataValues(i + q->nDummy(), 1.0/(1.0 + exp(80.0*(z - 0.15))));
  }

  MTKGridManager *gm = new MTKGridManager();

  MTKGrid1D *qinit = gm->MultData1D (1, q);

  for (i = 0; i < numItt; i++) {

    MTKGrid1D *tmp = gm->MultData1D (1, q);

    Step = 3.0;
    for (cont = 3; cont >= 1; cont--) {

      MTKGrid1D *uq = gm->MultData1D (u, q);
      gm->Apply1DPeriodic(uq);
      MTKGrid1D *divuq = gm->CGM1DDiv(uq);
      q = gm->Set2Sum1D (tmp, divuq, (-dt/Step));
      Step--;

      delete uq;
      delete divuq;
    }

    CreateSnapshotInTime(i, Cr, q, qinit);

    delete tmp;
  }

  // Create data file with final solution:
  datfile = fopen("MainAdvec.dat", "w");
  if (datfile == NULL) {
    fprintf(stderr, "Error creating file\n");
  }

  for (i=0; i<q->nNodes()+2*q->nDummy(); i++) {
    if (i < q->nDummy()) {
      fprintf(datfile, "%f\t%f\t%f\n", q->GetGridValue(i), 0.0, 0.0);
    } else {
      if (i > q->nNodes() - 1) {
        fprintf(datfile, "%f\t%f\t%f\n", q->GetGridValue(i), 0.0, 0.0);
      } else {
        fprintf(datfile, "%f\t%f\t%f\n", q->GetGridValue(i),
                qinit->GetDataValue(i + q->nDummy()),
                q->GetDataValue(i + q->nDummy()));
      }
    }
  }
  fclose(datfile);

  GenerateAndSeeAnimation();

  return(0);
}
