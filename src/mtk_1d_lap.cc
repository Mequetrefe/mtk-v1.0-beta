/*!
\file mtk_1d_lap.cc

\brief Constructs a default MTK_1DCGMLaplacian.

Constructs a default MTK_1DCGMLaplacian. Note how the default order is
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
#include "mtk_1d_lap.h"

using namespace std;
using namespace mtk;


MTK_1DLap::MTK_1DLap() {
  kk_=mtk::MTK_DefaultOrder;
  mimetic_tol_= mtk::MTK_DefaultMimTol;
  order_=kk_;
  
}

MTK_1DLap::MTK_1DLap(int order) {
  kk_=order;
  mimetic_tol_= mtk::MTK_DefaultMimTol;

  
}


MTK_1DLap::MTK_1DLap(int order, mtk::MTK_Real mim_tol) {

  kk_=order;
  mimetic_tol_=mim_tol;
  
}


MTK_1DLap::~MTK_1DLap() {

  
}

MTK_1DLap* MTK_1DLap::Construct1DLap() {
  
  return this;
}



MTK_DenseMatrix* MTK_1DLap::ReturnAsMatrix(int NumCells, MTK_Real A, MTK_Real B)
{
  
   int m=NumCells+2;  //Laplacian dimension is m * m
   MTK_DenseMatrix*    lap;
   MTK_DenseMatrix*    dd;
   MTK_DenseMatrix*    gg;

   mtk::MTK_1DDiv *div = new mtk::MTK_1DDiv(kk_,mimetic_tol_);
   mtk::MTK_1DGrad *grad = new mtk::MTK_1DGrad(kk_,mimetic_tol_);
  
  lap=new MTK_DenseMatrix(m,m);

   div = div->Construct1DDiv();
     if (div == nullptr) {
    return nullptr;
  }

  grad = grad->Construct1DGrad();
  if (grad == nullptr) {
    return nullptr;
  }
  
  
  
  if (lap==nullptr)
   {
      cout << "Problem allocating memory for lap MTK_DenseMatrix"<<endl;
     return lap;
   }
   gg=grad->ReturnAsMatrix(NumCells, A, B); //Dim NumCells+1
   dd=div->ReturnAsMatrix(NumCells, A, B); // Dim NumCells+2;
    
// Matrix Multiplication routine
    lap=lap->MatrixMultiplication(*dd,*gg);
     return lap;
}



MTK_Real* MTK_1DLap::ReturnAsReal(int NumCellsX, MTK_Real West, MTK_Real East)
{

  MTK_DenseMatrix *DD;
  MTK_Real *DDR;
  DD=ReturnAsMatrix(NumCellsX, West, East);
  DDR=DD->DenseMatrix_To_MTK_Real(*DD);
  return DDR;
}

