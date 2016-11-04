// /*!
// \file mtk_2d_div.cc
// 
// \brief Constructs a default MTK_2DDiv.
// 
// Constructs a default MTK_2DDiv. Note how the default order is
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
#include "mtk_2d_div.h"


using namespace std;
using namespace mtk;


MTK_2DDiv::MTK_2DDiv() {
  kk_=mtk::MTK_DefaultOrder;
  mimetic_tol_= mtk::MTK_DefaultMimTol;
  order_=kk_;
  
}

MTK_2DDiv::MTK_2DDiv(int order) {
  kk_=order;
  mimetic_tol_= mtk::MTK_DefaultMimTol;

  
}


MTK_2DDiv::MTK_2DDiv(int order, mtk::MTK_Real mim_tol) {

  kk_=order;
  mimetic_tol_=mim_tol;
  
}


MTK_2DDiv::~MTK_2DDiv() {

  
}

MTK_2DDiv* MTK_2DDiv::Construct2DDiv() {
  
  return this;
}



MTK_DenseMatrix* MTK_2DDiv::ReturnAsMatrix(int NumCellsX, MTK_Real West, MTK_Real East, int NumCellsY, MTK_Real South, MTK_Real North)
{
  int mx=NumCellsX+2;  //Dx vertical dimension 
  int nx=NumCellsX+1;  //Dx horizontal dimension
  int my=NumCellsY+2;  //Dy vertical dimension
  int ny=NumCellsY+1;  //Dy horizontal dimension

  mtk::MTK_DenseMatrix* Dx; //Dense Matrix to store the 1D Divergence
  mtk::MTK_DenseMatrix* Dy; //Dense Matrix to store the 1D Divergence
  mtk::MTK_DenseMatrix* Ix; // Extended Identity matrix rank = NumcellsX
  mtk::MTK_DenseMatrix* Iy; //Extenden Identity matrix rank=NumCellsY
  mtk::MTK_DenseMatrix* Dxy; //DivX 2D
  mtk::MTK_DenseMatrix* Dyx; //DivY 2D
  mtk::MTK_DenseMatrix* D2D; //Div 2D

  mtk::MTK_1DDiv *div = new mtk::MTK_1DDiv(kk_,mimetic_tol_); 
  
  div = div->Construct1DDiv();
     if (div == nullptr) {
    return nullptr;
  }


  Dx=div->ReturnAsMatrix(NumCellsX, West, East); // Dx;
  Dy=div->ReturnAsMatrix(NumCellsY, South, North); // Dy;
  Ix=new MTK_DenseMatrix(NumCellsX, true);
  Iy=new MTK_DenseMatrix(NumCellsY, true);

  Dxy = new MTK_DenseMatrix(mx*my,nx*NumCellsY);
  Dyx = new MTK_DenseMatrix(mx*my,ny*NumCellsX);
  D2D = new MTK_DenseMatrix(mx*my,nx*NumCellsY + ny*NumCellsX );
  
//  Dimensions of MTK_2DDiv would be mx*my rows and (nx*NumCellsY+by*NumCellsX) columns  
#if MTK_DEBUG_LEVEL > 0  
  
  cout << "Dx :" << mx << "by " << nx << endl;

  cout << "Iy : " << my<< " by " << NumCellsY  << endl;

  

  cout << "Dy :" << my << "by " << ny << endl;

  cout << "Iy : " << my<< " by " << NumCellsY  << endl;

  
  cout << "Kronecker dimensions Div 2D" << mx*my << " by " << nx*NumCellsY + ny*NumCellsX <<endl;
    
#endif
  
  Dxy->Kronecker(*Iy,*Dx, *Dxy);

  Dyx->Kronecker(*Dy,*Ix, *Dyx);

  
  for(auto ii=0;ii<mx*my;ii++)
  {
    for(auto jj=0;jj<nx*NumCellsY;jj++)
    {
      D2D->SetValue(ii,jj,Dxy->GetValue(ii,jj));
    }
    for(auto kk=0; kk<ny*NumCellsX; kk++)
    {
        D2D->SetValue(ii,kk+nx*NumCellsY,Dyx->GetValue(ii,kk) );
    }
  }
  
  return D2D;
}



MTK_Real* MTK_2DDiv::ReturnAsReal(int NumCellsX, MTK_Real West, MTK_Real East, int NumCellsY, MTK_Real South, MTK_Real North)
{

  MTK_DenseMatrix *D2D;
  MTK_Real *D2DR;
  D2D=ReturnAsMatrix(NumCellsX, West, East, NumCellsY, South, North);
  D2DR=D2D->DenseMatrix_To_MTK_Real(*D2D);
  return D2DR;
}

