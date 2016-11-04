/*!
 * \file mtk_ccs_sparse_matrix.cc
 *
 * \brief Definition of the default constructor for the
 * MTK_CCSSparseMatrix class.
 *
 * This file contains the definition of the default constructor for the
 * MTK_CCSSparseMatrix.
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

#include "mtk_ccs_matrix.h"

using namespace std;
using namespace mtk;

/*! Constructs a default MTK_CCSSparseMatrix instance. */

/*! Constructs a default MTK_CCSSparseMatrix instance.
 */
MTK_CCSSparseMatrix::MTK_CCSSparseMatrix(void):
  values(NULL),
  row_ind(NULL),
  col_ptr(NULL),
  num_rows(0),
  num_cols(0),
  number_values(0),
  number_non_zero(0),
  number_zero(0),
  unitary_offset(0) {

}

/*! Constructs MTK_CCSSparseMatrix from a 1D array dense matrix. */

/*! Constructs MTK_CCSSparseMatrix from a dense matrix which is implemented
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
MTK_CCSSparseMatrix::MTK_CCSSparseMatrix( MTK_Real *array,
                                                    int num_rows,
                                                    int num_cols,
                                                    int unitary_offset) {

  int ii;       /* Iterator. */
  int jj;       /* Iterator. */
  int nnz_jj;   /* Iterator for the nnz values array. */
  int rpt_jj;   /* Iterator for the row_ptr array. */
  int new_col;  /* Are we about to create a new col?. */

  /* Verify that the matrix that we try to create actually exists: */
  if (array == (MTK_Real*) NULL) {
    fprintf(stderr, "Incorrect parameter at line %d.\n", __LINE__ - 1);
  }
  /* Verify first that the requested sizes make sense: */
  if (num_rows <= 0) {
    fprintf(stderr, "Incorrect parameter at line %d.\n", __LINE__ - 1);
  }
  if (num_cols <= 0) {
    fprintf(stderr, "Incorrect parameter at line %d.\n", __LINE__ - 1);
  }

  /* Verify that the unitary offset equals 1 or 0: */
  if (unitary_offset == 1 || unitary_offset == 0) {
    this->unitary_offset = unitary_offset;
  } else {
    fprintf(stderr, "Incorrect parameter at line %d.\n", __LINE__ - 3);
  }

  this->number_non_zero = 0;
  this->num_rows = num_rows;
  this->num_cols = num_cols;
  this->number_values = num_rows*num_cols;
  /* 1. Traverse the given matrix to determine how many non-zero values are: */
  for (ii = 0; ii < (this->number_values); ii++) {
    if (array[ii] != 0.0) {
      (this->number_non_zero)++;
    } else {
      (this->number_zero)++;
    }
  }
  /* 2. Allocate resources: */
  /* 2.1. Allocate values array: */
  this->values = (MTK_Real *) malloc (sizeof(MTK_Real)*(this->number_non_zero));
  if (this->values == (MTK_Real *) NULL) {
    fprintf(stderr, "Incorrect parameter at line %d.\n", __LINE__ - 1);
  }
  /* 2.2. Allocate row indexes array: */
  this->row_ind = (int *) malloc (sizeof(int)*(this->number_non_zero));
  if (this->row_ind == (int *) NULL) {
    fprintf(stderr, "Incorrect parameter at line %d.\n", __LINE__ - 1);
  }
  /* 2.3. Allocate row pointers array: */
  /* We need an extra memory space to store the number of non-zero values: */
  this->col_ptr = (int *) malloc (sizeof(int)*(this->num_rows + 1));
  if (this->col_ptr == (int *) NULL) {
    fprintf(stderr, "Incorrect parameter at line %d.\n", __LINE__ - 1);
  }
  /* 3. Analyze and fill the new sparse matrix: */
  nnz_jj = 0;
  rpt_jj = 0;
  /* For each row in the given matrix: */
  for (ii = 0; ii < (this->num_rows); ii++) {
    new_col = 1;
    /* For each column in the given matrix: */
    for (jj = 0; jj < (this->num_cols); jj++) {
      /* Analyze current element: */
      if (array[jj*(this->num_cols) + ii] != 0.0) {
        if (new_col) {
          new_col = 0;
          this->col_ptr[rpt_jj] = nnz_jj + this->unitary_offset;
          rpt_jj++;
        }
        this->values[nnz_jj] = array[jj*(this->num_cols) + ii];
        this->row_ind[nnz_jj] = jj + this->unitary_offset;
        nnz_jj++;
      }
    }
  }
  /* 4. At the end, the final position in row_ptr stores the number of non-zero
   e lements:* */
  this->col_ptr[this->num_rows] = this->number_non_zero;
}

/*! Prints a MTK_CCSSparseMatrix as a collection of arrays. */

/*! Prints a MTK_CCSSparseMatrix as a collection of arrays.
 */
void MTK_CCSSparseMatrix::MTKArrayPrint(void) {

  int ii;         // Iterator.
  int groups_of;  // Number of values per output row.

  groups_of = 5;
  cout <<  "Number of rows: " << this->num_rows << endl;
  cout <<  "Number of cols: " << this->num_cols << endl;
  cout <<  "Number of non-zero values: " << this->number_non_zero << endl;
  cout <<  "Non zero values (1D array printed in groups of 5): " << endl << endl;
  for (ii = 0; ii < this->number_non_zero; ii++) {
    if ((ii != 0) && (ii % groups_of == 0)) {
      cout <<  endl;
    }
    cout <<  this->values[ii] << " ";
  }
  cout <<  endl << endl << "Row indexes:" << endl << " ";
  for (ii = 0; ii < this->number_non_zero; ii++) {
    cout <<  " " << this->row_ind[ii];
  }
  cout << endl << "Column pointers:" << endl << " ";
  for (ii = 0; ii < this->num_cols; ii++) {
    cout << this->col_ptr[ii] << " ";
  }
  cout <<  endl;
  cout <<  "Again, number of non-zero values:" << this->col_ptr[this->num_cols] << endl;
  cout <<  endl;
}

/*! Prints a MTK_CCSSparseMatrix as a dense block matrix. */

/*! Prints a MTK_CCSSparseMatrix as a dense block matrix.
 * \param *array INPUT: Shows dots instead of the values so that the sparse
 * structure can be appreciated.
 * \todo Thursday, March 15 2012 11:07 PM: Implement converting algorithm.
 */
void MTK_CCSSparseMatrix::MTKBlockPrint(int dots){

  cout << "This one is tricky!" << endl;
}
