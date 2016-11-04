/*!
 * \file mtk_1d_array_dense_stencil_matrix.cc
 *
 * \brief Constructs a default MTK_1DArrayDenseStencilMatrix.
 *
 * Constructs a default MTK_1DArrayDenseStencilMatrix.
 *
 * \date: Monday, September 03, 2012
 * \version: 2012-09-03.
 * \author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
 */
 /* Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
 * State University (SDSU).
 *
 * http:www.csrc.sdsu.edu/mtk/
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification are permitted provided that the following conditions are met:
 *
 * -# Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * -# Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * -# Neither the name of the CSRC, SDSU nor the names of its contributors may
 * be used to endorse or promote products derived from this software without
 * specific prior written permission.
 * -# Modifications whether they are partial or complete; whether they are
 * additions or eliminations should be reported through an email at:
 * esanchez@sciences.sdsu.edu
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
#include "mtk_1d_cgm_dirichlet_operator.h"
#include "mtk_1d_cgm_neumann_operator.h"
#include "mtk_1d_cgm_gradient.h"
#include "mtk_1d_cgm_laplacian.h"
#include "mtk_1d_array_dense_stencil_matrix.h"

using namespace std;
using namespace mtk;
/*! Constructs a default MTK_1DArrayDenseStencilMatrix. */

/*! Constructs a default MTK_1DArrayDenseStencilMatrix. Note how the default order is
 * at least 2.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 */
MTK_1DArrayDenseStencilMatrix::MTK_1DArrayDenseStencilMatrix(void):
  dense_values_(NULL),
  representation_(0),
  order_(2),
  size_(0),
  num_rows_(0),
  num_cols_(0) {

}

/*! Constructs a default MTK_1DArrayDenseStencilMatrix. */

/*! Constructs a default MTK_1DArrayDenseStencilMatrix. Note how the default order is
 * at least 2. It assumes that \f$ M=\alpha A+\beta BG-L \f$.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 * \todo Saturday, April 21 2012 07:02 PM: Validate order/size relation.
 */
MTK_1DArrayDenseStencilMatrix::MTK_1DArrayDenseStencilMatrix (MTK_Real alpha,
                                       MTK_1DCGMDirichletOperator *AA,
                                       MTK_Real beta,
                                       MTK_1DCGMNeumannOperator *BB,
                                       MTK_1DCGMGradient *GG,
                                       MTK_1DCGMLaplacian *LL) {

  MTK_Real *identity;
  int ii;
  int nn;

  this->order_ = LL->order();
  this->size_ = LL->size();
  this->representation_ = 0;
  this->num_rows_ = size_ + 2;
  this->num_cols_ = size_ + 2;

  /* Allocate space for the values: */
  nn = this->size_;
  this->dense_values_ = (MTK_Real *) malloc (sizeof(MTK_Real)*(nn + 2)*(nn + 2));
  if (this->dense_values_ == (MTK_Real *) NULL) {
    fprintf(stderr, "Allocation error at line %d.\n", __LINE__ - 2);
  }
  for (ii = 0; ii <= (((nn + 2)*(nn + 2)) - 1); ii++) {
    this->dense_values_[ii] = 0.0;
  }
  /* Define the stencil: */
  this->dense_values_[0] = 1.0;
  this->dense_values_[(((nn + 2)*(nn + 2)) - 1)] = 1.0;

  /* Prepare the linear system: This mean, do: AA = alpha*AA + beta*BB*GG: */
  /* We will compute AA = beta*BB*GG + alpha*AA, using only one routine: */
  int ldBB = (nn + 1);
  int ldGG = (nn + 2);
  int ldAA = (nn + 2);
  int ldLL = (nn + 2);
  cblas_dgemm(  CblasRowMajor, CblasNoTrans, CblasNoTrans,
                (nn + 2), (nn + 2), (nn + 1), beta, BB->dense_values_, ldBB,
                GG->dense_values_, ldGG, alpha, AA->dense_values_, ldAA);

  /* Since the result gets stored in AA, we need to copy to this! */
  for (ii = 0; ii <= (((nn + 2)*(nn + 2)) - 1); ii++) {
    this->dense_values_[ii] = AA->dense_values_[ii];
  }

  /* Generate a nn*nn Identity matrix: */
  identity = (MTK_Real *) calloc ((nn + 2)*(nn + 2), sizeof(MTK_Real));
  if (identity == (MTK_Real *) NULL) {
    fprintf(stderr, "Allocation error at line %d.\n", __LINE__ - 2);
  }
  for (ii = 0; ii < (nn + 2); ii++) {
    identity[ii*(nn + 2) + ii] = 1.0;
  }

//   printf("Created identity:\n");
//   for (ii = 0; ii < (nn + 2); ii++) {
//     for (jj = 0; jj < (nn + 2); jj++) {
//       printf("%8.4f ", identity[ii*(nn + 2) + jj]);
//     }
//     putchar('\n');
//   }
//   putchar('\n');

  /* We now compute the subtraction of the Laplacian: M = M - L. */
  cblas_dgemm(  CblasRowMajor, CblasNoTrans, CblasNoTrans,
              (nn + 2), (nn + 2), (nn + 2), 1.0, identity, (nn + 2),
                this->dense_values_, (nn + 2), -1.0, LL->dense_values_, ldLL);

  /* Since the result gets stored in LL, we need to copy to this! */
  for (ii = 0; ii <= (((nn + 2)*(nn + 2)) - 1); ii++) {
    this->dense_values_[ii] = LL->dense_values_[ii];
  }
}

/*! Constructor #3.  */
/*! This functions CONSTRUCTS a stencil matrix based on common dims.
  \todo Make sure we validate things!
 */
MTK_1DArrayDenseStencilMatrix::MTK_1DArrayDenseStencilMatrix(int nrows, int ncols) {

  int ii;

  this->size_ = (nrows);
  num_rows_ = nrows;
  num_cols_ = ncols;
  order_ = 2;
  representation_ = 0;
  this->dense_values_ = (MTK_Real *) malloc (sizeof(MTK_Real)*(nrows + 2)*(ncols + 2));
  if (this->dense_values_ == (MTK_Real *) NULL) {
    fprintf(stderr, "Allocation error at line %d.\n", __LINE__ - 2);
  }
  for (ii = 0; ii <= (((nrows + 2)*(ncols + 2)) - 1); ii++) {
    this->dense_values_[ii] = 0.0;
  }
}

/*! Prints. */

/*! Prints.
 * \todo Friday, April 20 2012 11:27 PM: Change the representation to an enum.
 */
void MTK_1DArrayDenseStencilMatrix::Print(int display) {

  int ii;
  int jj;
  int nn;

  switch (this->representation_) {
    case 0: {
      cout << "1D Stencil Matrix... dense representation:"  << endl;
      nn = this->size_ + 2;
      for (ii = 0; ii < nn; ii++) {
        for (jj = 0; jj < nn; jj++) {
          cout << setw(10) << this->dense_values_[ii*nn + jj];
        }
        cout << endl;
      }
      cout << setw(10) << this->order_ << " order stencil matrix." << endl;
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
