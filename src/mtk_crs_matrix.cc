/*!
 * \file mtk_crs_sparse_matrix.cc
 *
 * \brief Definition of the default constructor for the
 * MTK_CRSSparseMatrix class.
 *
 * This file contains the definition of the default constructor for the
 * MTK_CRSSparseMatrix.
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
#include <iomanip>

#include "mtk_tool_manager.h"
#include "mtk_crs_matrix.h"
#include "mtk_roots.h"

using namespace std;
using namespace mtk;	

/*! Constructs a default MTK_CRSSparseMatrix instance. */

/*! Constructs a default MTK_CRSSparseMatrix instance.
 */
MTK_CRSSparseMatrix::MTK_CRSSparseMatrix(void):
  values(NULL),
  col_ind(NULL),
  row_ptr(NULL),
  num_rows(0),
  num_cols(0),
  number_values(0),
  number_non_zero(0),
  number_zero(0),
  unitary_offset(0) {

}

/*! Constructs MTK_CRSSparseMatrix from a 1D array dense matrix. */

/*! Constructs MTK_CRSSparseMatrix from a dense matrix which is implemented
 * as a 1D array.
 * \param *array INPUT: Pointer to the beginning of the 1D array implementing
 * the dense matrix.
 * \param num_rows INPUT: Number of rows of the dense matrix.
 * \param num_cols INPUT: Number of cols of the dense matrix.
 * \param unitary_offset INPUT: States if the dense matrix should be treated as
 * if it were indexed in C-style or in Fortran style -that is in \f$[0,n-1]\f$
 * or in \f$[1,n]\f$- for constructing purposes. A unitary offset set to zero,
 * will consider the dense matrix as it were indexed from 0, therefore the
 * attained sparse matrix will have information in row 0 and column 0, whereas a
 * unitary offset of 1 will NOT include information from row 0 or column 0, but
 * will include information from row \f$m\f$ and column \f$n\f$, where \f$m\f$
 * and \f$n\f$ are the total number of rows or columns, respectively.
 * \todo Thursday, March 15 2012 11:07 PM: Verify if this actually works by
 * comparing it with a working example.
 */
MTK_CRSSparseMatrix::MTK_CRSSparseMatrix(
  mtk::MTK_Real *array,
  int num_rows,
  int num_cols,
  int unitary_offset) {

  int ii;                     // Iterator.
  int jj;                     // Iterator.
  int nnz_jj;                 // Iterator for the nnz values array.
  int rpt_jj;                 // Iterator for the row_ptr array.
  int offset;                 // Offset of the 1D array-based matrix.
  int new_row;                // Are we about to create a new row?
//  MTKToolManager *tool_mngr;  // Used to handle execution exceptions.

  // Initialize the tool manager:
  //tool_mngr = new MTKToolManager;

  // Verify that the matrix that we try to create actually exists: */
  if (array == (MTK_Real*) NULL) {
//     tool_mngr->MTKAbort("Incorrect parameter", MTK_THIS_LINE - 1, __FILE__);
  }
  // Verify first that the requested sizes make sense:
  if (num_rows <= 0) {
//     tool_mngr->MTKAbort("Incorrect parameter", MTK_THIS_LINE - 1, __FILE__);
  }
  if (num_cols <= 0) {
//     tool_mngr->MTKAbort("Incorrect parameter", MTK_THIS_LINE - 1, __FILE__);
  }
  this->number_non_zero = 0;
  this->number_zero = 0;
  this->num_rows = num_rows;
  this->num_cols = num_cols;
  this->number_values = num_rows*num_cols;
  this->unitary_offset = unitary_offset;
  // Traverse the given matrix to determine how many non zero values are:
  for (ii = 0; ii < this->number_values; ii++) {
    if (array[ii] != 0.0) {
      (this->number_non_zero)++;
    } else {
      (this->number_zero)++;
    }
  }

  // Allocate resources:
  this->values = (MTK_Real *) malloc (sizeof(MTK_Real)*this->number_non_zero);
  if (this->values == (MTK_Real *) NULL) {
//     tool_mngr->MTKAbort("Allocation error.", MTK_THIS_LINE - 2, __FILE__);
  }
  this->col_ind = (MTK_Real *) malloc (sizeof(MTK_Real)*this->number_non_zero);
  if (this->col_ind == (MTK_Real *) NULL) {
//     tool_mngr->MTKAbort("Allocation error.", MTK_THIS_LINE - 2, __FILE__);
  }
  this->row_ptr = (MTK_Real *) malloc (sizeof(MTK_Real)*this->num_rows);
  if (this->row_ptr == (MTK_Real *) NULL) {
//     tool_mngr->MTKAbort("Allocation error.", MTK_THIS_LINE - 2, __FILE__);
  }
  // Analyze and fill the new sparse matrix:
  /* 3. Analyze and fill the new sparse matrix: */
  nnz_jj = 0;
  rpt_jj = 0;
  /* For each row in the given matrix: */
  for (ii = 0; ii < (this->num_rows); ii++) {
    offset = ii*(this->num_cols);
    new_row = true;
    /* For each column in the given matrix: */
    for (jj = 0; jj < (this->num_cols); jj++) {
      /* Analyze current element: */
      if (array[offset + jj] != 0.0) {
        if (new_row) {
          new_row = 0;
          this->row_ptr[rpt_jj] = nnz_jj + this->unitary_offset;
          rpt_jj++;
        }
        this->values[nnz_jj] = array[offset + jj];
        this->col_ind[nnz_jj] = jj + this->unitary_offset;
        nnz_jj++;
      }
    }
  }
  /* 4. At the end, the final position in row_ptr stores the number of non-zero
   e lements:* */
  this->row_ptr[this->num_rows] = this->number_non_zero;
}

/*! Prints a MTK_CRSSparseMatrix as a collection of arrays. */

/*! Prints a MTK_CRSSparseMatrix as a collection of arrays.
 */
void MTK_CRSSparseMatrix::MTKArrayPrint(void) {

  int ii; /* Iterator. */

  fprintf(stdout, "Number of rows: %d\n", this->num_rows);
  fprintf(stdout, "Number of cols: %d\n", this->num_cols);
  fprintf(stdout, "Number of non-zero values: %d\n", this->number_non_zero);
  fprintf(stdout, "Non zero values:\n  ");
  for (ii = 0; ii < this->number_non_zero; ii++) {
    fprintf(stdout, "%.4f ", this->values[ii]);
  }
  fprintf(stdout, "\nColumn indexes:\n  ");
  for (ii = 0; ii < this->number_non_zero; ii++) {
    fprintf(stdout, "%.4f ", this->col_ind[ii]);
  }
  fprintf(stdout, "\nRow pointers:\n  ");
  for (ii = 0; ii < this->num_rows; ii++) {
    fprintf(stdout, "%.4f ", this->row_ptr[ii]);
  }
  fprintf(stdout, "\n");
}

/*! Prints a MTK_CRSSparseMatrix as a dense block matrix.
 * \param *array INPUT: Shows dots instead of the values so that the sparse
 * structure can be appreciated.
 * \todo Thursday, March 15 2012 11:07 PM: Implement converting algorithm just
 * to store the matrix in memory and let the user decide if he wants to store
 * the matrix.
 */
void MTK_CRSSparseMatrix::MTKBlockPrint(int dots) {

  int ii;         /* Iterator. */
  int jj;         /* Iterator. */
  int col_jj;     /* Iterator for the values array. */
  int val_jj;     /* Iterator for the columns indexes array. */

  val_jj = 0;
  col_jj = 0;
  for (ii = 0; ii < this->num_rows; ii++) {
    for (jj = 0; jj <  this->num_cols; jj++) {
      // If we've found a non-zero element:
      if (jj == (this->col_ind[col_jj] - this->unitary_offset)) {
        if (dots) {
          cout << "* ";
        } else {
          cout << setw(10) << this->values[val_jj];
        }
        val_jj++;
        col_jj++;
      } else {
        if (dots) {
          cout << " ";
        } else {
          cout << setw(10) << (MTK_Real) 0.0;
        }
      }
    }
    cout << endl;
  }
}
