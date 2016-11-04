/*!
\file ex13-dok_sparse_matrices_driver.cc

\brief Example on using the DOK sparse matrices data structures in MTK.

This file performs basic operations in DOK sparse matrices as defined in MTK.

\date: Monday, September 03, 2012
\version: 2011-09-03.

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

#include <iostream>
#include <iomanip>

#include "mtk.h"

using namespace std;

/*!
\brief MAIN module.
*/
int main () {

  double *sample_matrix;
  int rand_val;
  int n_rows;
  int n_cols;
  int size;
  int ii;
  int jj;
  int offset;

  srand (time(NULL));
  cout << "Testing sparse matrices in MTK." << endl;

  // Building the original matrix:
  n_rows = 4;
  n_cols = 4;
  size = n_rows*n_cols;
  sample_matrix = (double *) malloc(sizeof(double)*size);
  if (sample_matrix == (double *) NULL) {
    cerr << "ERROR: Allocation of the orig. matrix failed!\n" << endl;
  }

  // Initializing the original matrix:
  for (ii = 0; ii < size/2; ii++) {
    rand_val = rand() % 100;
    if (rand_val % 2 == 0) {
      sample_matrix[ii] = 0.0;
    } else {
      sample_matrix[ii] = rand_val;
    }
    sample_matrix[(size - 1) - ii] = sample_matrix[ii];
  }

  // Printing the original matrix:
  cout << "Original matrix:" << endl;
  for (ii = 0; ii < n_rows; ii++) {
    offset = ii*n_rows;
    for (jj = 0; jj < n_cols; jj++) {
      cout << setw(10) << sample_matrix[offset + jj] << " ";
    }
    cout << endl;
  }

  // Create DOK sparse matrix:
  MTK_DOKSparseMatrix *sample_sparse_matrix = new MTK_DOKSparseMatrix();

  for (ii = 0; ii < n_rows; ii++) {
    offset = ii*n_rows;
    for (jj = 0; jj < n_cols; jj++) {
      sample_sparse_matrix->add_value(ii, jj, sample_matrix[offset + jj]);
    }
    cout << endl;
  }

  sample_sparse_matrix->print_dense();
  sample_sparse_matrix->print_sparse(5);

  delete sample_sparse_matrix;

  free(sample_matrix);
  return EXIT_SUCCESS;
}
