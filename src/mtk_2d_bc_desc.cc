/*!
  \file mtk_2d_bc_desc.cc

  */
/* Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
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

#include <iostream>
#include <iomanip>

#include "mtk_2d_bc_desc.h"
#include "mtk_roots.h"

using namespace std;
using namespace mtk;

MTK_2DBCDesc::MTK_2DBCDesc(void):
  number_of_cells_X_(0),
  number_of_cells_Y_(0),
  East_(0),
  West_(0),
  South_(0),
  North_(0),
  lap_(nullptr)
  {

}

MTK_2DBCDesc::MTK_2DBCDesc(
                 MTK_2DLap *lap,
                 MTK_2DUniformStaggeredGrid *GRID,
                 MTK_Real (*South)(MTK_Real xx, MTK_Real yy),
                 MTK_Real (*North)(MTK_Real xx, MTK_Real yy),
                 MTK_Real (*West)(MTK_Real xx, MTK_Real yy),
                 MTK_Real (*East)(MTK_Real xx, MTK_Real yy) )
{
    // Extract the dimension from the laplacian matrix
  number_of_cells_X_ = lap->number_of_cells_X_;
  East_=lap->East_;
  West_=lap->West_;
   
  number_of_cells_Y_ = lap->number_of_cells_Y_ ;
  South_=lap->South_;
  North_=lap->North_;
  
  lap_=lap;

  fW_=West;
  fE_=East;
  fS_=South;
  fN_=North;
  GRID_=GRID;
  
  cout<<number_of_cells_X_<<endl;
  cout<<number_of_cells_Y_<<endl;

    
}

// MTK_2DBCDirichlet takes a 2D dense matrix Laplacian
// and fill the boundary conditions required to get Dirichlet conditions.

void MTK_2DBCDesc::MTK_2DBCDirichlet(MTK_DenseMatrix *LL)
{

    

    for(auto ii =0; ii < number_of_cells_X_+2;ii++)
    {
        // South
        LL->SetValue(ii,ii,1.0);
    }
    
    auto counter = number_of_cells_X_+1;
    
    for(auto ii =LL->matrix_data_->num_rows_-1; ii >= LL->matrix_data_->num_rows_-number_of_cells_X_-2;ii--)
    {
        //North
        LL->SetValue(ii,ii,1.0);
        counter--;
    }
    
    counter=0;
    for(auto ii =number_of_cells_X_+2; ii <LL->matrix_data_->num_rows_-number_of_cells_X_-2; ii=ii+number_of_cells_X_+2)
    {
        //West
        LL->SetValue(ii,ii,1.0);
        //East
        LL->SetValue(ii + number_of_cells_X_+1,ii+number_of_cells_X_+1,1.0);
        counter++;
    }

    
}

void MTK_2DBCDesc::MTK_2DBoundaryRHS(MTK_Real *RHS )
    {

        
      for (auto ii=0; ii < number_of_cells_X_+2; ii++)
      {
// South
	RHS[ii]=fS_(GRID_->valuesX(ii),South_);

// North
          RHS[(number_of_cells_X_+ 2)*(number_of_cells_Y_+1)+ii]=fN_(GRID_->valuesX(ii),North_);
	
      }

      for (auto ii=0; ii < number_of_cells_Y_; ii++)
      {
// West
	RHS[(ii+1)*(number_of_cells_X_+2)]=fW_(West_,GRID_->valuesY(ii));
//East
	RHS[(ii+1)*(number_of_cells_X_+2)+number_of_cells_X_+1]=fE_(East_,GRID_->valuesY(ii+1));
	
      }


      
    }
