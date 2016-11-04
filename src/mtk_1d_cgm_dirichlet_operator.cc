/*!
\file mtk_1d_cgm_dirichlet_operator.cc

\brief Constructs a default MTK_1DCGMDirichletOperator.

Constructs a default MTK_1DCGMDirichletOperator. Note how the default order is
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

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>

#include "mtk_roots.h"
#include "mtk_1d_cgm_dirichlet_operator.h"

using namespace std;


/*! Default constructor.  */

/*! This functions CONSTRUCTS a default Dirichlet operator.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 */
MTK_1DCGMDirichletOperator::MTK_1DCGMDirichletOperator(void):
  dense_values_(NULL),
  representation_(0),
  order_(2),
  size_(0) {

}

/*! Constructor #2. */
/*! Constructs an operator starting from the usual parameters.
    \todo Saturday, April 21 2012 07:02 PM: Validate order/size relation.
  */
MTK_1DCGMDirichletOperator::MTK_1DCGMDirichletOperator(int order, int size,
                                                   int representation) {

  int ii;
  int nn;

  this->order_ = order;
  this->size_ = size;
  this->representation_ = representation;

  /* Allocate space for the values: */
  nn = this->size_;
  this->dense_values_ = (MTK_Number *) malloc (sizeof(MTK_Number)*(nn + 2)*(nn + 2));
  if (this->dense_values_ == (MTK_Number *) NULL) {
    fprintf(stderr, "Allocation error at line %d.\n", __LINE__ - 2);
  }
  for (ii = 0; ii <= (((nn + 2)*(nn + 2)) - 1); ii++) {
    this->dense_values_[ii] = 0.0;
  }
  /* Define the stencil: */
  this->dense_values_[0] = 1.0;
  this->dense_values_[(((nn + 2)*(nn + 2)) - 1)] = 1.0;
}

/*! Prints an operator accordingly. */

/*! Prints an operator accordingly.
  \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 */
void MTK_1DCGMDirichletOperator::Print(int display) {

  int ii;
  int jj;
  int nn;

  switch (this->representation_) {
    case 0: {
      cout << "1D CG Dirichlet... dense representation:"  << endl;
      nn = this->size_ + 2;
      for (ii = 0; ii < nn; ii++) {
        for (jj = 0; jj < nn; jj++) {
          cout << setw(10) << this->dense_values_[ii*nn + jj];
        }
        cout << endl;
      }
      cout << setw(10) << this->order_ << " order CG Dirichlet." << endl;
      cout << setw(10) << "\t Dimensions: (" <<
        this->size_ << " + 2) x (" << this->size_ << " + 2) " << endl;
      break;
    }
    case 1: {
      if (display == 0) {
        cout << "1D CG Dirichlet... sparse rep; dense display:"  << endl;
      } else {
        cout << "1D CG Dirichlet... sparse rep: sparse display:"  << endl;
      }
      break;
    }
  }
}
