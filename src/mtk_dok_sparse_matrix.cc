/*!
 * \file mtk_dok_sparse_matrix.cc
 *
 * \brief Definition of the default constructor for the
 * MTK_DOKSparseMatrix class.
 *
 * This file contains the definition of the default constructor for the
 * MTK_DOKSparseMatrix.
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
#include <cstring>

#include <iostream>

#include "mtk_dok_sparse_matrix.h"

using namespace std;
using namespace mtk;

/*! Constructs a default MTK_DOKSparseMatrix instance. */

/*! Constructs a default MTK_DOKSparseMatrix instance.
 */
MTK_DOKSparseMatrix::MTK_DOKSparseMatrix(void):
    values_(NULL),
    tol_(0.0),
    rows_(NULL),
    cols_(NULL),
    memory_chunk_(0),
    num_rows_(0),
    num_cols_(0),
    num_elements_(0),
    num_non_zero_(0),
    adapted_for_mumps_(0),
    allocated_memory_values_(0),
    allocated_memory_rows_cols_(0) {

}

/*! Adds.  */
/*! Adds a value.
  \param rr Row.
  \param cc Col.
  \param vv Value.
*/
int MTK_DOKSparseMatrix::add_value(int rr, int cc, double vv) {

  int pp_rr;  // Position in the rows array.
  int pp_cc;  // Position in the cols array.

  if (vv == 0.0) {
    cerr << "Makes no sense to store a 0 value.\n";
    return MTK_FALSE;
  }

  // Do we need to re-size the matrix?
  if (num_non_zero_ >= num_elements_) {
    // If we need to resize/initialize the matrix, we will reallocate the
    // descriptive arrays.
    double *tmp_values;
    int *tmp_rows;
    int *tmp_cols;
    int tmp_num_elements = (num_elements_) + (memory_chunk_);
    int tmp_requested_memory_values = tmp_num_elements*sizeof(double);
    int tmp_requested_memory_rows_cols = tmp_num_elements*sizeof(int);

    // If this were the first time, we would be allocating these for the first
    // time, since realloc works like a malloc for the first time.
    tmp_values = (double*) realloc(values_, tmp_requested_memory_values);
    if (tmp_values != NULL) {
      memset(((char *)(tmp_values)) + allocated_memory_values_,
             0.0,
             tmp_requested_memory_values - allocated_memory_values_);
    } else {
      return MTK_FALSE;
    }

    tmp_rows = (int*) realloc(rows_, tmp_requested_memory_rows_cols);
    if (tmp_rows != NULL) {
      memset(((char *)(tmp_rows)) + allocated_memory_rows_cols_,
              -1,
              tmp_requested_memory_rows_cols - allocated_memory_rows_cols_);
    } else {
      return MTK_FALSE;
    }

    tmp_cols = (int*) realloc(cols_, tmp_requested_memory_rows_cols);
    if (tmp_cols != NULL) {
      memset(((char *)(tmp_cols)) + allocated_memory_rows_cols_,
              -1,
              tmp_requested_memory_rows_cols - allocated_memory_rows_cols_);
    } else {
      return MTK_FALSE;
    }

    allocated_memory_values_ = tmp_requested_memory_values;
    allocated_memory_rows_cols_ = tmp_requested_memory_rows_cols;

    // Check if the allocation works. Since we are considering two cases at the
    // same time, we need to do the following:
    // 1. Is there any problem?:
    if (!tmp_values || !tmp_rows || !tmp_cols) {

      if (tmp_values)
        free(tmp_values);
      else if (values_)
        free(values_);

      if (tmp_rows)
        free(tmp_rows);
      else if (rows_)
        free(rows_);

      if (tmp_cols)
        free(tmp_cols);
      else if (cols_)
        free(cols_);

      values_ = NULL;
      rows_ = NULL;
      cols_ = NULL;
      num_rows_ = 0;
      num_cols_ = 0;
      num_elements_ = 0;
      num_non_zero_ = 0;
    }

    // If there weren't any problems with allocating the new extra space, update
    // the attributes with the new allocated space:
    values_ = tmp_values;
    rows_ = tmp_rows;
    cols_ = tmp_cols;
    num_elements_ = tmp_num_elements;
  }

  // Regardless of requiring a re-size OR NOT, we are now ready to add the new
  // value:
  // Is this a new row?
  if (rr >= num_rows_) {
    num_rows_++;
  }
  // Is this a new column?
  if (cc >= num_cols_) {
    num_cols_++;
  }

  // Locate position in the values_ array:
  pp_rr = 0;
  pp_cc = 0;
  // Search first for the position in the rows_ array, since our matrices will
  // be assumed to be traversed row-wise:
  int curr;
  int tail = 0;
  int head = num_elements_ - 1;
  while(head >= tail) {
    curr = midpoint(head, tail);
    // In this case, not to found an element means it belongs at the beginning:
    if (curr == head) {
      // But we must ensure that we are indeed using the far extreme that has a
      // zero value:

      while (values_[tail] != 0.0 && rr >= rows_[tail]) {
        tail++;
      }
      pp_rr = tail;
      break;
    }
    if (rr < rows_[curr] && values_[curr] != 0.0) {
      head = curr;  // Search lower array.
    } else {
      if (rr > rows_[curr] && values_[curr] != 0.0) {
        tail = curr;  // Search upper array.
      } else {
        /* Notice that the following conditional, REQUIRES for value to be
        initialized to 0's all along!!! God, I could use some sort of recalloc!
        */
        if (values_[curr] == 0.0) {
          head = curr;  // Search lower array.
        } else {
          pp_rr = curr;
          break;
        }
      }
    }
  }

  // Now search for the position in the cols_ array, but assuming it starts
  // in the same position we have just found for the rows_ array:
  tail = pp_rr;
  head = num_elements_ - 1;
  while(head >= tail) {
    curr = midpoint(head, tail);
    // In this case, not to found an element means it belongs at the beginning:
    if (curr == head) {
      // But we must ensure that we are indeed using the far extreme that has a
      // zero value:
      while (values_[tail] != 0.0 && cc >= cols_[tail]) {
        tail++;
      }
      pp_cc = tail;
      break;
    }
    if (cc < cols_[curr] && values_[curr] != 0.0) {
      head = curr;  // Search lower array.
    } else {
      if (cc > cols_[curr] && values_[curr] != 0.0) {
        tail = curr;  // Search upper array.
      } else {
        /* Notice that the following conditional, REQUIRES for value to be
        initialized to 0's all along!!! God, I could use some sort of recalloc!
        */
        if (values_[curr] == 0.0) {
          head = curr;  // Search lower array.
        } else {
          pp_cc = curr;
          break;
        }
      }
    }
  }
  // In this point, the value has to be stored in values_[pp_cc]!
  // But first, if pp_cc <= (num_non_zero_ - 1), then we need to move all the
  // elements, one position, to make room for the new element:
  if (pp_cc <= (num_non_zero_ - 1)) {
    cout << "Adapting memory for this input!\n";
  }

  // Insert the value:
  values_[pp_cc] = vv;
  rows_[pp_cc] = rr;
  cols_[pp_cc] = cc;
  // Increment the number of non-zero values:
  num_non_zero_++;

  return MTK_TRUE;
}

/*! Gets.  */
/*! Gets a value.
  \param rr Row.
  \param cc Col.
*/
double MTK_DOKSparseMatrix::get_value(int rr, int cc) {

  int pp_rr;
  int ii;

  if (!(adapted_for_mumps_)) {
    /* Search for the rows in the ORDERED arrays of rows: */
    /* Search for the rows in the ORDERED arrays of rows: */
    pp_rr = 0;
    while (rows_[pp_rr] != rr && pp_rr < num_non_zero_) {
      pp_rr++;
    }
    if (pp_rr == num_non_zero_) {
      return 0.0;
    }

    /* If we have found the row coordinate, look for the value of the column
    coordinate; however, cols_ may NOT be ordered, so perform a linear search:
    */
    ii = pp_rr;
    while (cols_[ii] != cc && ii < num_non_zero_) {
      ii++;
    }
    if (rows_[ii] == rr && cols_[ii] == cc) {
      return values_[ii];
    } else {
      return 0.0;
    }

  } else {

    rr++;
    cc++;

    /* Search for the rows in the ORDERED arrays of rows: */
    /* Search for the rows in the ORDERED arrays of rows: */
    pp_rr = 0;
    while (rows_[pp_rr] != rr && pp_rr < num_non_zero_) {
      pp_rr++;
    }
    if (pp_rr == num_non_zero_) {
      return 0.0;
    }

    /* If we have found the row coordinate, look for the value of the column
    coordinate; however, cols_ may NOT be ordered, so perform a linear search:
    */
    ii = pp_rr;
    while (cols_[ii] != cc && ii < num_non_zero_) {
      ii++;
    }
    if (rows_[ii] == rr && cols_[ii] == cc) {
      return values_[ii];
    } else {
      return 0.0;
    }

  }
}

/*! Prints.  */
/*! Prints.
*/
void MTK_DOKSparseMatrix::print_sparse(int groups_of) {

  cout << "Printing matrix in groups of" << groups_of << "...\n";
  if (!(values_)) {
    cerr << "The matrix has NOT been initialized.\n";

  } else if (!(num_elements_)) {
    cerr << "The matrix has NOT been filled yet.\n";

  } else {
    int ii; // Iterator.

    cout << "[" << num_rows_ << "x" << num_cols_ << "]\n";
    cout << "[" << num_non_zero_ << "nz out of" << num_elements_ << "]\n";
    if (adapted_for_mumps_)
      cout << "This matrix is adapted for MUMPS! Notice indexing!\n";
    cout << "Values:\n";
    for (ii = 0; ii < num_non_zero_; ii++) {
      cout << values_[ii] << ", ";
      if ((ii + 1) % groups_of == 0) {
        cout << endl;
      }
    }
    putchar('\n');
    cout << "Rows:\n";
    for (ii = 0; ii < num_non_zero_; ii++) {
      cout << rows_[ii] << " ";
      if ((ii + 1) % groups_of == 0) {
        cout << endl;
      }
    }
    cout << endl;
    cout << "Cols:\n";
    for (ii = 0; ii < num_non_zero_; ii++) {
      cout << cols_[ii] << ", ";
      if ((ii + 1) % groups_of == 0) {
        cout << endl;
      }
    }
    cout << endl;
  }
  cout << endl;
}

/*! Prints.  */
/*! Prints.
*/
void MTK_DOKSparseMatrix::print_dense(void) {

  cout << "Printing matrix...\n";
  if (!(values_)) {
    cerr << "The matrix has NOT been initialized.\n";

  } else if (!(num_elements_)) {
    cerr << "The matrix has NOT been filled yet.\n";

  } else {
    int ii; // Iterator.
    int jj; // Iterator.

    cout << "[" << num_rows_ << "x" << num_cols_ << "]\n";
    cout << "[" << num_non_zero_ << "nz out of" << num_elements_ << "]\n";
    if (adapted_for_mumps_) {
      cerr << "This matrix is adapted for MUMPS! Notice indexing!\n";
    }
    for (ii = 0; ii < num_rows_; ii++) {
      for (jj = 0; jj < num_cols_ ; jj++) {
        cout << this->get_value(ii, jj);
      }
      cout << endl;
    }
    cout << endl;
  }
}

/*! Adapts.  */
/*! Adapts.
*/
int MTK_DOKSparseMatrix::adapt_for_mumps(void) {

  if (adapted_for_mumps_) {
    cerr << "Can't adapt for MUMPS more than once!\n";
    return MTK_FALSE;
  } else {
    int ii; // Iterator.
    for (ii = 0; ii < num_non_zero_; ii++) {
      rows_[ii]++;
      cols_[ii]++;
    }
    adapted_for_mumps_ = MTK_TRUE;
    return MTK_TRUE;
  }
}

/*! Undo.  */
/*! Undo.
*/
int MTK_DOKSparseMatrix::undo_adapt_for_mumps(void) {

  if (!(adapted_for_mumps_)) {
    cerr << "Can't UNDO the adapt for MUMPS more than once!\n";
    return MTK_FALSE;
  } else {
    int ii; // Iterator.
    for (ii = 0; ii < num_non_zero_; ii++) {
      rows_[ii]--;
      cols_[ii]--;
    }
    adapted_for_mumps_ = MTK_FALSE;
    return MTK_TRUE;
  }
}

/*! Destroys.  */
/*! This procedure destroys.
*/
MTK_DOKSparseMatrix::~MTK_DOKSparseMatrix() {

  if (rows_)
    free(rows_);
  else
    cerr << "Can't destroy!\n";

  if (values_)
    free(values_);
  else
    cerr << "Can't destroy!\n";

  if (cols_)
    free(cols_);
  else
    cerr << "Can't destroy!\n";

  memory_chunk_ = 0;
  num_rows_ = 0;
  num_cols_ = 0;
  num_elements_ = 0;
  num_non_zero_ = 0;
  adapted_for_mumps_ = MTK_FALSE;
  rows_  = NULL;
  values_  = NULL;
  cols_ = NULL;
}
