/*!
\file mtk_2d_cgm_gradient.cc

\brief Constructs a default MTK_2DCGMGradient.

Constructs a default MTK_2DCGMGradient. Note how the default order is
at least 2.

\date: Thursday, September 12, 2013
\version: 2013-09-12.
\author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
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

#ifdef MACOS_SOLVERS_ON
  #include <Accelerate/Accelerate.h>
#else
extern "C" {
  #include <cblas.h>
}
#endif

#include <iostream>
#include <iomanip>

#include "mtk_roots.h"
#include "mtk_2d_cgm_gradient.h"

using namespace std;

/*! Constructs a default MTK_2DCGMGradient operator. */

/*! Constructs a default MTK_2DCGMGradient. Note how the default order is
 * at least 2.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 */
MTK_2DCGMGradient::MTK_2DCGMGradient(void):
  dense_values_(NULL),
  representation_(0),
  order_(2),
  size_(0) {

}

/*! Constructs a MTK_2DCGMGradient operator. */

/*! Constructs a default MTK_2DCGMGradient. Note how the default order is
 * at least 2.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 * \todo Saturday, April 21 2012 07:02 PM: Validate order/size relation.
 */
MTK_2DCGMGradient::MTK_2DCGMGradient(int order,
                                     int Nx, int Ny,
                                     double hx, double hy) {

  int ii;
  int nn;

  int size_D_in = (Nx + 1)*Ny + Nx*(Ny + 1);
  this->rows_ = (Nx + 1)*Ny + Nx*(Ny + 1);
  int size_G_in = (Nx + 2)*Ny + 2*Nx;
  this->cols_ = (Nx + 2)*Ny + 2*Nx;

  int row_idx;
  int col_idx;

  this->order_ = order;
  this->representation_ = 0;

  cout << "Creating CGM Gradient..." << endl;
  cout << "size_D_in = " << size_D_in << endl;
  cout << "size_G_in = " << size_G_in << endl;

  /* Allocate space for the values: */
  nn = size_D_in*size_G_in;
  this->size_ = nn;
  this->dense_values_ = NULL;
  this->dense_values_ = (MTK_Number *) calloc (nn,sizeof(MTK_Number));
  if (this->dense_values_ == (MTK_Number *) NULL) {
    fprintf(stderr, "Allocation error at line %d.\n", __LINE__ - 2);
  } else {
    fprintf(stdout, "Allocation successful at line %d.\n", __LINE__ - 4);
  }

  /* Vertical blocks: */
  for (int kk = 1; kk <= Ny; kk++) {
    row_idx = (kk-1)*(Nx+1) + 1;
    row_idx--;
    col_idx = (kk-1)*(Nx+2) + 1;
    col_idx--;

    /* First row: -8/3 3 -1/3 */
    dense_values_[row_idx*cols_ + col_idx + 0] = -(8.0/3.0)/hx;
    dense_values_[row_idx*cols_ + col_idx + 1] =  3.0/hx;
    dense_values_[row_idx*cols_ + col_idx + 2] = -(1.0/3.0)/ hx;

    for (ii = 1; ii <= (Nx - 0); ii++) {
      dense_values_[(row_idx+ii-1)*cols_ + (col_idx+ii-1) + 0] = -1.0/hx;
      dense_values_[(row_idx+ii-1)*cols_ + (col_idx+ii-1) + 1] = 1.0/hx;
    }
    // Last row: 1/3 -3 8/3
    int local_last_row = row_idx + Nx;
    int local_last_col = col_idx + Nx + 1;
    dense_values_[local_last_row*cols_ + local_last_col-2] = 1.0/3.0 / hx;
    dense_values_[local_last_row*cols_ + local_last_col-1] = -3.0 / hx;
    dense_values_[local_last_row*cols_ + local_last_col] = 8.0/3.0 / hx;
  }

  // Horizontal Blocks
  int sec_col_idx;
  int local_h_offset;
  int local_last_row;
  int offset_v = Ny*(Nx + 1);
  int offset_h = Ny*(Nx + 2);
  for (int kk = 1; kk <= Nx; kk++) {
    row_idx = (kk-1)*(Ny+1) + offset_v + 1;
    row_idx--;
    col_idx = (kk-1)*2 + offset_h + 1;
    col_idx--;
    sec_col_idx = (kk-1) + 2;
    sec_col_idx--;
    // First component: -8/3 3 -1/3
    dense_values_[row_idx*cols_ + col_idx] = -8.0/3.0 / hy;
    dense_values_[row_idx*cols_ + sec_col_idx] = 3.0 / hy;
    dense_values_[row_idx*cols_ + sec_col_idx + Nx + 2] = -1.0/3.0 / hy;
    local_h_offset = sec_col_idx;
    for (ii = 1; ii <= (Ny - 1); ii++) {
      dense_values_[(row_idx+ii)*cols_ + local_h_offset] = -1.0/hy;
      dense_values_[(row_idx+ii)*cols_ + local_h_offset+Nx+2] = 1.0/hy;
      local_h_offset = local_h_offset + Nx + 2;
    }
    // Last component: 1/3 -3 8/3
    local_last_row = row_idx + Ny;
    dense_values_[local_last_row*cols_ + local_h_offset-(Nx+2)] = 1.0/3.0 / hy;
    dense_values_[local_last_row*cols_ + local_h_offset] = -3.0 / hy;
    dense_values_[local_last_row*cols_ + col_idx + 1] = 8.0/3.0 / hy;
  }
}

/*! Prints a MTKIDCGGradient operator. */

/*! Prints a default MTKIDCGGradient.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 */
void MTK_2DCGMGradient::Print(int display) {

  int ii;
  int jj;
  int nn;

  switch (this->representation_) {
    case 0: {
      cout << "2D CG Gradient... dense representation:"  << endl;
      nn = this->size_;
      cout << "nn = " << nn << endl;
      for (ii = 0; ii < rows_; ii++) {
        for (jj = 0; jj < cols_; jj++) {
          cout << setw(5) << setprecision(4) <<
            this->dense_values_[ii*cols_ + jj];
        }
        cout << endl;
      }
      cout << setw(10) << this->order_ << " order 2D CG Gradient." << endl;
      cout << setw(10) << "\t Dimensions: " <<
        this->rows_ << " x " << this->cols_ << endl;
      break;
    }
    case 1: {
      if (display == 0) {
        cout << "1D CG Neumann... sparse rep; dense display:"  << endl;
      } else {
        cout << "1D CG Neumann... sparse rep: sparse display:"  << endl;
      }
      break;
    }
  }
}
