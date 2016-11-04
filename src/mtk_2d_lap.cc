/*!
\file mtk_2d_lap.cc

\brief Constructs a default MTK_2DLap.

Constructs a default MTK_2DMLap. Note how the default order is
at least 2.

\date: Monday, September 03, 2012
\version: 2012-09-03.
\author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
*/
/* Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
State University (SDSU).

http:www.csrc.sdsu.edu/mtk/

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

#ifdef MACOS_SOLVERS_ON
  #include <Accelerate/Accelerate.h>
#else
extern "C" {
  #include <cblas.h>
}
#endif

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>

#include "mtk_roots.h"
#include "mtk_2d_lap.h"
#include "mtk_dense_matrix.h"

using namespace std;
using namespace mtk;


MTK_2DLap::MTK_2DLap() {
  kk_=mtk::MTK_DefaultOrder;
  mimetic_tol_= mtk::MTK_DefaultMimTol;
  order_=kk_;

}

MTK_2DLap::MTK_2DLap(int order) {
  kk_=order;
  mimetic_tol_= mtk::MTK_DefaultMimTol;


}


MTK_2DLap::MTK_2DLap(int order, mtk::MTK_Real mim_tol) {

  kk_=order;
  mimetic_tol_=mim_tol;

}

MTK_2DLap::MTK_2DLap(int order, mtk::MTK_Real mim_tol, int NumCellsX, mtk::MTK_Real West, mtk::MTK_Real East, int NumCellsY, mtk::MTK_Real South, mtk::MTK_Real North) {

  kk_=order;
  mimetic_tol_=mim_tol;
  number_of_cells_X_=NumCellsX;
  number_of_cells_Y_=NumCellsY;
  East_=East;
  West_= West;
  South_=South;
  North_=North;
}




MTK_2DLap::~MTK_2DLap() {


}

MTK_2DLap* MTK_2DLap::Construct2DLap() {

  return this;
}


MTK_DenseMatrix* MTK_2DLap::ReturnAsMatrix(int NumCellsX, MTK_Real West, MTK_Real East, int NumCellsY, MTK_Real South, MTK_Real North)
{

   int m=(NumCellsX+2)*(NumCellsY+2);
 //  int n=NumCellsY+1;
   MTK_DenseMatrix*    lap;
   MTK_DenseMatrix*    dd;
   MTK_DenseMatrix*    gg;

   mtk::MTK_2DDiv *div = new mtk::MTK_2DDiv(kk_,mimetic_tol_);
   mtk::MTK_2DGrad *grad = new mtk::MTK_2DGrad(kk_,mimetic_tol_);
   number_of_cells_X_=NumCellsX;
   number_of_cells_Y_=NumCellsY;
        East_ = East;
       West_=West;
      South_=South;
      North_=North;


   div = div->Construct2DDiv();
     if (div == nullptr) {
    return nullptr;
  }

  grad = grad->Construct2DGrad();
  if (grad == nullptr) {
    return nullptr;
  }


   gg=grad->ReturnAsMatrix(NumCellsX, West, East,NumCellsY, South, North);
   dd=div->ReturnAsMatrix(NumCellsX, West, East,NumCellsY, South, North);

   lap=new MTK_DenseMatrix(m,m);


  if (lap==nullptr)
   {
      cout << "Problem allocating memory for MTK_DenseMatrix"<<endl;
     return lap;
   }



   lap=lap->MatrixMultiplication(*dd,*gg);

  return lap;
}

MTK_Real* MTK_2DLap::ReturnAsReal(int NumCellsX, MTK_Real West, MTK_Real East, int NumCellsY, MTK_Real South, MTK_Real North)
{

  MTK_DenseMatrix *lap;
  MTK_Real *lapR;
  lap=ReturnAsMatrix(NumCellsX, West, East, NumCellsY, South, North);

    lapR=lap->DenseMatrix_To_MTK_Real(*lap);

  return lapR;
}

int MTK_2DLap::LapSize(void)
{
    return (number_of_cells_X_+2)*(number_of_cells_Y_+2);
}