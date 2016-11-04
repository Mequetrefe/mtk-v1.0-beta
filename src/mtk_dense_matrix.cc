/*!
\file mtk_dense_matrix.cc

\brief Defines a common dense matrix, using a 1-D array.

For developing purposes, it is better to have a not-so-intrincated data
structure implementing matrices. This is the purpose of this class: to be used
for prototypes of new code for small test cases. In every other instance, this
should be replaced by the most appropriate sparse matrix.

\date: Monday, September 03, 2012
\version: 2012-09-03.
\author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
 */
/*
Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
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

esanchez at sciences dot sdsu dot edu

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
#include <cstring>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "mtk_roots.h"
#include "mtk_dense_matrix.h"

using namespace std;
using namespace mtk;

/*! Constructs a default MTK_DenseMatrix instance. */

/*! Constructs a default MTK_DenseMatrix instance. */



MTK_DenseMatrix::MTK_DenseMatrix(void):
 Values_(nullptr),
  continuum_(nullptr),
 nRow_(0),
 nCol_(0) 
 {

}


/*! Constructs MTK_DenseMatrix from a 1D array dense matrix. */

/*! Constructs MTKCCSSparseMatrix from a dense matrix which is implemented
 * as a 1D array.
 * \param *array INPUT: Pointer to the beginning of the 1D array implementing
 * the dense matrix.
 * \param num_rows INPUT: Number of rows of the dense matrix.
 * \param num_cols INPUT: Number of cols of the dense matrix.
 * \param unitary_offset INPUT: States if ...
 */
MTK_DenseMatrix::MTK_DenseMatrix(const int num_rows, const int num_cols) {

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

/*! Constructs MTK_DenseMatrix from a 1D array dense matrix. */

/*! Constructs MTKCCSSparseMatrix from a dense matrix which is implemented
 * as a 1D array.
 * \param *array INPUT: Pointer to the beginning of the 1D array implementing
 * the dense matrix.
 * \param num_rows INPUT: Number of rows of the dense matrix.
 * \param num_cols INPUT: Number of cols of the dense matrix.
 * \param unitary_offset INPUT: States if ...
 */
void MTK_DenseMatrix::PrintBlock(void) {

  int ii;
  int jj;

  fprintf(stdout, "Matrix... %d rows x %d cols.\n", nRow_, nCol_);
  fprintf(stdout, "Matrix... %d elements.\n", number_values_);
  for (ii = 0; ii < nRow_; ii++) {
    for (jj = 0; jj < nCol_; jj++) {
      fprintf(stdout, "%5.5g ", Values_[ii][jj]);
    }
    putchar('\n');
  }
}

/*! Constructs MTK_DenseMatrix from a 1D array dense matrix. */

/*! Constructs MTKCCSSparseMatrix from a dense matrix which is implemented
 * as a 1D array.
 * \param *array INPUT: Pointer to the beginning of the 1D array implementing
 * the dense matrix.
 * \param num_rows INPUT: Number of rows of the dense matrix.
 * \param num_cols INPUT: Number of cols of the dense matrix.
 * \param unitary_offset INPUT: States if ...
 */
MTK_Real MTK_DenseMatrix::GetValue(int rr, int cc) {

  if (0 > rr || rr > number_values_ || 0 > cc || cc > number_values_) {
    fprintf(stderr, "Error!\n");
  }
  return Values_[rr][cc];
}


void MTK_DenseMatrix::SetValue(int rr,int cc,MTK_Real val)
{
	  if (0 > rr || rr > number_values_ || 0 > cc || cc > number_values_) {
	    fprintf(stderr, "Error!\n");
	  }
	Values_[rr][cc]=val;
}
/*! Constructs MTK_DenseMatrix from a 1D array dense matrix. */

/*! Constructs MTKCCSSparseMatrix from a dense matrix which is implemented
 * as a 1D array.
 * \param *array INPUT: Pointer to the beginning of the 1D array implementing
 * the dense matrix.
 * \param num_rows INPUT: Number of rows of the dense matrix.
 * \param num_cols INPUT: Number of cols of the dense matrix.
 * \param unitary_offset INPUT: States if ...
 */
MTK_DenseMatrix::~MTK_DenseMatrix() {

  free(MemoryHolder_);
  free(Values_);
}

// Used to compute a 1D index (idx) from a 2D set of indexes.
// WARNING: Not needed if we allocate the matrices so that they are contiguous
// in memory.
// TODO: Implement the MTK_DenseMatrix class with this constructor, so that we
// do not have to use this routine.



inline int idx(const int ii, const int offset, const int jj) {

  return ii*offset + jj;
}

// Transposes a matrix AA with rr rows and cc columns:

double* Transpose(double *AA, int rr, int cc) {

  double* AT {};  // Output transpose matrix.

  if (rr <= 0 || cc <= 0 || AA == nullptr) {
    return nullptr;
  }
  try {
    AT = new double[rr*cc];
    memset(AT, 0.0, sizeof(AT[0])*(rr*cc));
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  for (auto ii = 0; ii < rr; ii++) {
    for (auto jj = 0; jj < cc; jj++) {
      AT[idx(jj,rr,ii)] = AA[idx(ii,cc,jj)];
    }
  }

  return AT;
}

// DEF. In linear algebra, a VANDERMONDE MATRIX is a matrix with the terms of a
// geometric progression in each row. This progression uses the terms of a
// given GENERATOR VECTOR.

// Generates a Vandermonde matrix with a given generator vector gen:

double* Vandermonde(const double *gen,      // Given generator vector.
                    const int gen_length,   // Length of the generator vector.
                    const int pro_length,   // Length of the progression.
                    const bool transpose) { // Should create its transpose?

  double* TT {};  // Output Vandermonde matrix.

  // Check for the integrity of the arguments:
  if (gen == nullptr || gen_length < 1 || pro_length < 1) {
    return nullptr;
  }

  try {
    TT = new double[gen_length*pro_length];
    memset(TT, 0.0, sizeof(TT[0])*gen_length*pro_length);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  if (!transpose) {
    for (auto ii = 0; ii < gen_length; ii++) {
      for (auto jj = 0; jj < pro_length; jj++) {
        TT[idx(ii,pro_length,jj)] = pow(gen[ii], (double) jj);
      }
    }
  } else {
    for (auto ii = 0; ii < pro_length; ii++) {
      for (auto jj = 0; jj < gen_length; jj++) {
        TT[idx(ii,gen_length,jj)] = pow(gen[jj], (double) ii);
      }
    }
  }

  return TT;
}

// Prints a dense matrix:

bool MTK_DenseMatrix_Print(const double* aa,      // The matrix.
                           const int num_rows,    // Number of rows.
                           const int num_cols) {  // Number of columns.

  if (aa == nullptr || num_rows < 1 || num_cols < 1) {
    cerr << "ERROR printing matrix of line " << __LINE__ << endl;
    cerr << "Exiting..." << endl;
    return false;
  }

  for (auto ii = 0; ii < num_rows; ii++) {
    for (auto jj = 0; jj < num_cols; jj++) {
      cout << setw(12) << aa[idx(ii,num_cols,jj)];
    }
    cout << endl;
  }
  cout << endl;

  return true;
}

// In the MTK version, we should only receive two dense matrix objects, which
// already contain their information, regarding then number of rows and
// columns. Which is why, our constructor should only receive the two matrix
// operands.
double *kron(const double *aa,
             const int &mm,
             const int &nn,
             const double *bb,
             const int &pp,
             const int &qq) {

           double *kk{};            // Output matrix, as a 1D array.
  register double aa_factor{};      // Kept in register (KiR).
           int kk_num_rows{mm*pp};  // Number of rows of output matrix.
  register int kk_num_cols{nn*qq};  // Number of cols of output matrix (KiR).
  register int row_offset{};        // Row-offset to access output matrix (KiR).
  register int col_offset{};        // Col-offset to access output matrix (KiR).

  // We do not perform any validation on the sizes of the matrix, because it
  // is assumed that these will make sense since they will be provided within
  // the input matrix objects, which are assumed to always be in a consistent
  // state.

  kk = new double[kk_num_rows*kk_num_cols];

  for (auto ii = 0; ii < mm; ++ii) {
    row_offset = ii*pp;
    for (auto jj = 0; jj < nn; ++jj) {
      col_offset = jj*qq;
      aa_factor = aa[ii*nn + jj];
      for (auto ll = 0; ll < pp; ++ll) {
        for (auto oo = 0; oo < qq; ++oo) {
          kk[(ll + row_offset)*kk_num_cols + (oo + col_offset)] =
            aa_factor*bb[ll*qq + oo];
        }
      }
    }
  }

  // In the MTK, we will NOT return. This will be a constructor.
  return kk;
}
