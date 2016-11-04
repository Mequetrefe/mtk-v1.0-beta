/*!
 * \file mtk_dense_matrix.cc
 *
 * \brief Definition of the default constructor for the
 * MTKCCSSparseMatrix class.
 *
 * This file contains the definition of the default constructor for the
 * MTKCCSSparseMatrix.
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

#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "mtk_roots.h"
#include "mtk_2d_array_dense_matrix.h"

using namespace std;
using namespace mtk;

/*! Constructs a default MTK_2DArrayDenseMatrix instance. */

/*! Constructs a default MTK_2DArrayDenseMatrix instance. */
MTK_2DArrayDenseMatrix::MTK_2DArrayDenseMatrix(void):
  Values_(NULL),
  MemoryHolder_(NULL),
  nRow_(0),
  nCol_(0) {

}

/*! Constructs MTK_2DArrayDenseMatrix from a 1D array dense matrix. */

/*! Constructs MTKCCSSparseMatrix from a dense matrix which is implemented
 * as a 1D array.
 * \param *array INPUT: Pointer to the beginning of the 1D array implementing
 * the dense matrix.
 * \param num_rows INPUT: Number of rows of the dense matrix.
 * \param num_cols INPUT: Number of cols of the dense matrix.
 * \param unitary_offset INPUT: States if ...
 */
MTK_2DArrayDenseMatrix::MTK_2DArrayDenseMatrix(int num_rows,int num_cols) {

  int ii;

  /* Verify first that the requested sizes make sense: */
  if (num_rows <= 0) {
    fprintf(stderr, "Incorrect parameter at line %d.\n", __LINE__ - 1);
  }
  if (num_cols <= 0) {
    fprintf(stderr, "Incorrect parameter at line %d.\n", __LINE__ - 1);
  }

  nRow_ = num_rows;
  nCol_ = num_cols;
  number_values_ = nRow_*nCol_;

  MemoryHolder_ = (MTK_Real *) malloc(nRow_*nCol_*sizeof(MTK_Real));
  if (MemoryHolder_ == (MTK_Real *) NULL) {
    fprintf(stderr, "Allocation error at line %d.\n", __LINE__ - 2);
  }

  Values_ = (MTK_Real **) malloc(nRow_*sizeof(MTK_Real));
  if (Values_ == (MTK_Real **) NULL) {
    fprintf(stderr, "Allocation error at line %d.\n", __LINE__ - 2);
  }

  for (ii = 0; ii < nRow_; ii++) {
    Values_[ii] = MemoryHolder_ + (ii*nCol_);
  }
}

/*! Constructs MTK_2DArrayDenseMatrix from a 1D array dense matrix. */

/*! Constructs MTKCCSSparseMatrix from a dense matrix which is implemented
 * as a 1D array.
 * \param *array INPUT: Pointer to the beginning of the 1D array implementing
 * the dense matrix.
 * \param num_rows INPUT: Number of rows of the dense matrix.
 * \param num_cols INPUT: Number of cols of the dense matrix.
 * \param unitary_offset INPUT: States if ...
 */
void MTK_2DArrayDenseMatrix::PrintBlock() {

  int ii;
  int jj;

  fprintf(stdout, "Matrix... %d rows x %d cols.\n", nRow_, nCol_);
  fprintf(stdout, "Matrix... %d elements.\n", number_values_);
  for (ii = 0; ii < nRow_; ii++) {
    for (jj = 0; jj < nCol_; jj++) {
      fprintf(stdout, "%5.5g", Values_[ii][jj]);
    }
    putchar('\n');
  }
}

/*! Constructs MTK_2DArrayDenseMatrix from a 1D array dense matrix. */

/*! Constructs MTKCCSSparseMatrix from a dense matrix which is implemented
 * as a 1D array.
 * \param *array INPUT: Pointer to the beginning of the 1D array implementing
 * the dense matrix.
 * \param num_rows INPUT: Number of rows of the dense matrix.
 * \param num_cols INPUT: Number of cols of the dense matrix.
 * \param unitary_offset INPUT: States if ...
 */
MTK_Real MTK_2DArrayDenseMatrix::GetValue(int rr, int cc) {

  if (0 > rr || rr > number_values_ || 0 > cc || cc > number_values_) {
    fprintf(stderr, "Error!\n");
  }
  return Values_[rr][cc];
}

/*! Constructs MTK_2DArrayDenseMatrix from a 1D array dense matrix. */

/*! Constructs MTKCCSSparseMatrix from a dense matrix which is implemented
 * as a 1D array.
 * \param *array INPUT: Pointer to the beginning of the 1D array implementing
 * the dense matrix.
 * \param num_rows INPUT: Number of rows of the dense matrix.
 * \param num_cols INPUT: Number of cols of the dense matrix.
 * \param unitary_offset INPUT: States if ...
 */
MTK_2DArrayDenseMatrix::~MTK_2DArrayDenseMatrix() {

  free(MemoryHolder_);
  free(Values_);
}
