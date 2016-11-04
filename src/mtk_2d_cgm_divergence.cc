/*!
\file mtk_2d_cgm_divergence.cc

\brief Constructs a default MTK_2DCGMDivergence.

Constructs a default MTK_2DCGMDivergence. Note how the default order is
at least 2.

\date: Thursday, September 12, 2013
\version: 2013-09-12.
\author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
*/
/*
Copyright (C) 2013 Computational Science Research Center (CSRC) at San Diego
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
#include "mtk_2d_cgm_divergence.h"

using namespace std;

/*! Constructs a default MTK_2DCGMDivergence operator. */

/*! Constructs a default MTK_2DCGMDivergence. Note how the default order is
 * at least 2.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 */
MTK_2DCGMDivergence::MTK_2DCGMDivergence(void):
  dense_values_(NULL),
  representation_(0),
  order_(2),
  size_(0) {

}

/*! Constructs a MTK_2DCGMDivergence operator. */

/*! Constructs a default MTK_2DCGMDivergence. Note how the default order is
 * at least 2.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 * \todo Saturday, April 21 2012 07:02 PM: Validate order/size relation.
 */
MTK_2DCGMDivergence::MTK_2DCGMDivergence(int order,
                                     int Nx, int Ny,
                                     double hx, double hy) {

  int ii;
  int nn;

  int size_G_in = (Nx+2)*Ny+2*Nx;
  this->rows_ = (Nx+2)*Ny+2*Nx;
  int size_D_in = (Nx+1)*Ny+Nx*(Ny+1);
  this->cols_ = (Nx+1)*Ny+Nx*(Ny+1);

  int row_idx;
  int col_idx;

  this->order_ = order;
  this->representation_ = 0;

  /* Allocate space for the values: */
  nn = size_G_in*size_D_in;
  this->size_ = nn;
  this->dense_values_ = NULL;
  this->dense_values_ = (MTK_Number *) calloc (nn,sizeof(MTK_Number));
  if (this->dense_values_ == (MTK_Number *) NULL) {
    fprintf(stderr, "Allocation error at line %d.\n", __LINE__ - 2);
  } else {
    fprintf(stdout, "Allocation successful at line %d.\n", __LINE__ - 4);
  }

  /* Vertical blocks: */
  int offset_h = Ny*(Nx+1);
  int sec_col_idx;
  int local_col_idx;
  for (int kk = 1; kk <= Ny; kk++) {
    row_idx = (kk-1)*(Nx+2) + 1;
    row_idx--;
    col_idx = (kk-1)*(Nx+1);
    col_idx--;
    sec_col_idx = offset_h + kk;
    for (ii = 1; ii <= (Nx - 0); ii++) {
      local_col_idx = sec_col_idx + (ii-1)*(Ny+1) - 1;
      dense_values_[(row_idx+ii)*cols_ + (col_idx+ii)] = -1.0/hx;
      dense_values_[(row_idx+ii)*cols_ + (col_idx+ii+1)] = 1.0/hx;
      dense_values_[(row_idx+ii)*cols_ + (local_col_idx)] = -1.0/hy;
      dense_values_[(row_idx+ii)*cols_ + (local_col_idx+1)] = 1.0/hy;
    }
  }
}

/*! Prints a MTK2DCGDivergence operator. */

/*! Prints a default MTK2DCGDivergence.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 */
void MTK_2DCGMDivergence::Print(int display) {

  int ii;
  int jj;
  int nn;

  switch (this->representation_) {
    case 0: {
      cout << "2D CG Divergence... dense representation:"  << endl;
      nn = this->size_;
      cout << "nn = " << nn << endl;
      for (ii = 0; ii < rows_; ii++) {
        for (jj = 0; jj < cols_; jj++) {
          cout << setw(5) << setprecision(4) <<
            this->dense_values_[ii*cols_ + jj];
        }
        cout << endl;
      }
      cout << setw(10) << this->order_ << " order 2D CG Divergence." << endl;
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
