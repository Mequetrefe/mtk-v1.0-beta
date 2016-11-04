/*!
\file mtk_1d_cgm_laplacian.cc

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
#include "mtk_1d_cgm_laplacian.h"

using namespace std;

/*! Constructs a default MTK_1DCGMLaplacian. */

/*! Constructs a default MTK_1DCGMLaplacian. Note how the default order is
 * at least 2.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 */
MTK_1DCGMLaplacian::MTK_1DCGMLaplacian(void):
  dense_values_(NULL),
  representation_(0),
  order_(2),
  size_(0) {

}

/*! Constructs a default MTK_1DCGMLaplacian. */

/*! Constructs a default MTK_1DCGMLaplacian. Note how the default order is
 * at least 2.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 * \todo Saturday, April 21 2012 07:02 PM: Validate order/size relation.
 */
MTK_1DCGMLaplacian::MTK_1DCGMLaplacian(int order, int size, int representation) {

  int ii;
  int nn;

  this->order_ = order;
  this->size_ = size;
  this->representation_ = representation;

  /* Allocate space for the values: */
  nn = this->size_;
  this->dense_values_ = (MTK_Number *) calloc ((nn + 2)*(nn + 2),sizeof(MTK_Number));
  if (this->dense_values_ == (MTK_Number *) NULL) {
    fprintf(stderr, "Allocation error at line %d.\n", __LINE__ - 2);
  }
  /* Define the stencil: */
  this->dense_values_[nn + 2] = 8.0/3.0;
  this->dense_values_[nn + 3] = -4.0;
  this->dense_values_[nn + 4] = 4.0/3.0;
  for (ii = 2; ii < nn; ii++) {
      this->dense_values_[ii*(nn + 2) + ii - 1] = 1.0;
      this->dense_values_[ii*(nn + 2) + ii] = -2.0;
      this->dense_values_[ii*(nn + 2) + ii + 1] = 1.0;
  }
  this->dense_values_[(((nn + 2)*(nn + 2)) - 1) - (nn + 2) - 2] = 4.0/3.0;
  this->dense_values_[(((nn + 2)*(nn + 2)) - 1) - (nn + 2) - 1] = -4.0;
  this->dense_values_[(((nn + 2)*(nn + 2)) - 1) - (nn + 2)] = 8.0/3.0;

  /* Multiply by the number of cell centers: */
  /* Note that we multiply the Laplacian by N^2. */
  cblas_dscal((nn + 2)*(nn + 2), nn*nn, this->dense_values_, 1);
}

/*! Prints a MTKIDCGLaplacian operator. */

/*! Prints a default MTKIDCGLaplacian.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 */
void MTK_1DCGMLaplacian::Print(int display) {

  int ii;
  int jj;
  int nn;

  switch (this->representation_) {
    case 0: {
      cout << "1D CG Laplacian... dense representation:"  << endl;
      nn = this->size_ + 2;
      for (ii = 0; ii < nn; ii++) {
        for (jj = 0; jj < nn; jj++) {
          cout << setw(10) << this->dense_values_[ii*nn + jj];
        }
        cout << endl;
      }
      cout << setw(10) << this->order_ << " order CG Laplacian." << endl;
      cout << setw(10) << "\t Dimensions: (" <<
        this->size_ << " + 2) x (" << this->size_ << " + 2) " << endl;
      break;
    }
    case 1: {
      if (display == 0) {
        cout << "1D CG Laplacian... sparse rep; dense display:"  << endl;
      } else {
        cout << "1D CG Laplacian... sparse rep: sparse display:"  << endl;
      }
      break;
    }
  }
}
