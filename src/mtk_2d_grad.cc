// /*!
// \file mtk_2d_grad.cc
// 
// \brief Constructs a default MTK_2DGrad.
// 
// Constructs a default MTK_2DGrad. Note how the default order is
// at least 2.
// 
// \date: Thursday, September 12, 2013
// \version: 2013-09-12.
// \author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
// */
// /*
// Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
// State University (SDSU).
// 
// http://www.csrc.sdsu.edu/mtk/
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification are permitted provided that the following conditions are met:
// 
// -# Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// -# Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// -# Neither the name of the CSRC, SDSU nor the names of its contributors may
// be used to endorse or promote products derived from this software without
// specific prior written permission.
// -# Modifications whether they are partial or complete; whether they are
// additions or eliminations should be reported through an email at:
// esanchez@sciences.sdsu.edu
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// */
// 

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
#include "mtk_2d_grad.h"


using namespace std;
using namespace mtk;


MTK_2DGrad::MTK_2DGrad() {
  kk_=mtk::MTK_DefaultOrder;
  mimetic_tol_= mtk::MTK_DefaultMimTol;
  order_=kk_;
  
}

MTK_2DGrad::MTK_2DGrad(int order) {
  kk_=order;
  mimetic_tol_= mtk::MTK_DefaultMimTol;

  
}


MTK_2DGrad::MTK_2DGrad(int order, mtk::MTK_Real mim_tol) {

  kk_=order;
  mimetic_tol_=mim_tol;
  
}


MTK_2DGrad::~MTK_2DGrad() {

  
}

MTK_2DGrad* MTK_2DGrad::Construct2DGrad() {
  
  return this;
}



MTK_DenseMatrix* MTK_2DGrad::ReturnAsMatrix(int NumCellsX, MTK_Real West, MTK_Real East, int NumCellsY, MTK_Real South, MTK_Real North)
{
  int mx=NumCellsX+1;  //Gx vertical dimension 
  int nx=NumCellsX+2;  //Gx horizontal dimension
  int my=NumCellsY+1;  //Gy vertical dimension
  int ny=NumCellsY+2;  //Gy horizontal dimension

  mtk::MTK_DenseMatrix* Gx; //Dense Matrix to store the 1D Gradergence
  mtk::MTK_DenseMatrix* Gy; //Dense Matrix to store the 1D Gradergence
  mtk::MTK_DenseMatrix* Ix; // Extended Identity matrix rank = NumcellsX
  mtk::MTK_DenseMatrix* Iy; //Extenden Identity matrix rank=NumCellsY
  mtk::MTK_DenseMatrix* TIx; // Transposed Extended Identity matrix rank = NumcellsX
  mtk::MTK_DenseMatrix* TIy; // Transposed Extenden Identity matrix rank=NumCellsY
  mtk::MTK_DenseMatrix* Gxy; //GradX 2D
  mtk::MTK_DenseMatrix* Gyx; //GradY 2D
  mtk::MTK_DenseMatrix* G2D; //Grad 2D

  mtk::MTK_1DGrad *grad = new mtk::MTK_1DGrad(kk_,mimetic_tol_); 
  
  grad = grad->Construct1DGrad();
     if (grad == nullptr) {
    return nullptr;
  }


  Gx=grad->ReturnAsMatrix(NumCellsX, West, East); // Dx;
  Gy=grad->ReturnAsMatrix(NumCellsY, South, North); // Dy;
  Ix=new MTK_DenseMatrix(NumCellsX, true);
  Iy=new MTK_DenseMatrix(NumCellsY, true);
  
  TIx= new MTK_DenseMatrix(NumCellsX,NumCellsX+2);
  TIy= new MTK_DenseMatrix(NumCellsY,NumCellsY+2);
  
  TIx->Transpose(*Ix,*TIx);
  TIy->Transpose(*Iy,*TIy);

  Gxy = new MTK_DenseMatrix(mx*NumCellsY,nx*ny);
  Gyx = new MTK_DenseMatrix(my*NumCellsX,nx*ny);
  G2D = new MTK_DenseMatrix(mx*NumCellsY + my*NumCellsX,nx*ny );
#if MTK_DEBUG_LEVEL > 0  
  
  cout << "Gx :" << mx << "by " << nx << endl;

  cout << "Transpose Iy : " << NumCellsY<< " by " << ny  << endl;

  
  cout << "Gy :" << my << "by " << ny << endl;
  
  cout << "Transpose Ix : " << NumCellsX<< " by " << nx  << endl;
  
  
  cout << "Kronecker dimensions Grad 2D" <<mx*NumCellsY + my*NumCellsX << " by " <<  nx*ny <<endl;
#endif    
  
  Gxy->Kronecker(*TIy,*Gx, *Gxy);

  Gyx->Kronecker(*Gy,*TIx, *Gyx);

  
  for(auto ii=0;ii<nx*ny;ii++)
  {
    for(auto jj=0;jj<mx*NumCellsY;jj++)
    {
      G2D->SetValue(jj,ii,Gxy->GetValue(jj,ii));
        
    }
    for(auto kk=0; kk<my*NumCellsX; kk++)
    {
        G2D->SetValue(kk+mx*NumCellsY,ii,Gyx->GetValue(kk,ii) );
//        G2D->SetValue(kk+mx*NumCellsY,ii,314 );
    }
  }
  
  return G2D;
}



MTK_Real* MTK_2DGrad::ReturnAsReal(int NumCellsX, MTK_Real West, MTK_Real East, int NumCellsY, MTK_Real South, MTK_Real North)
{
   MTK_DenseMatrix *G2D;
   MTK_Real* G2DR;
   G2D=ReturnAsMatrix(NumCellsX, West, East, NumCellsY, South, North);
    G2DR=G2D->DenseMatrix_To_MTK_Real(*G2D);
  return G2DR;
}

