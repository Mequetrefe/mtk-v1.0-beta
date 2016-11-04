/*!
  \file mtk_dok_sparse_matrix.h

  \brief Includes the definition of the class MTK_DOKSparseMatrix.

  This file contains the definition of the class MTK_DOKSparseMatrix.
  This class implements a sparse matrix using Dictionary of Keys (DOK) to
  ensure storing only the non-zero elements of a given dense matrix. As a
  consequence of the selected storage scheme, no assumptions are being made in
  terms of the actual structure of the matrix.

  \date: Tuesday, March 19, 2013 11:04 PM
  \version: 2013-03-19-01.
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
  esanchez@sciences.sdsu.edu

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#ifndef MTK_INCLUDE_MTK_DOK_SPARSE_MATRIX_H
#define MTK_INCLUDE_MTK_DOK_SPARSE_MATRIX_H

#include "mtk_enums.h"

class MTK_DOKSparseMatrix {

  public:
    /*! Constructs.  */
    /*! This procedure constructs.
    */
    MTK_DOKSparseMatrix();

    /*! Destroys.  */
    /*! This procedure destroys.
    */
    ~MTK_DOKSparseMatrix();

    /*! Adds.  */
    /*! Adds a value.
      \param rr Row.
      \param cc Col.
      \param vv Value.
    */
    int add_value(int rr, int cc, double vv);

    /*! Gets.  */
    /*! Gets a value.
      \param rr Row.
      \param cc Col.
    */
    double get_value(int rr, int cc);

    /*! Prints.  */
    /*! Prints.
      \param groups_of Groups.
    */
    void print_sparse(int groups_of);

    /*! Prints.  */
    /*! Prints.
    */
    void print_dense(void);

    /*! Adapts.  */
    /*! Adapts.
    */
    int adapt_for_mumps(void);

    /*! Undo.  */
    /*! Undo.
    */
    int undo_adapt_for_mumps(void);

  private:
    int midpoint (int imin, int imax) {

      // This implementation will return the ceil of the attained real number.
      return imin + (imax - imin)/2;
    }

    int binary_search(int A[], int key, int imin, int imax) {

      int key_not_found = -1;

      // Continue searching while [imin,imax] is not empty:
      while (imax >= imin) {
        // Calculate the midpoint for roughly equal partition:
        int imid = midpoint(imin, imax);
        // determine which subarray to search
        if (A[imid] <  key) {
          // Change min index to search upper subarray:
          imin = imid + 1;
        } else if (A[imid] > key ) {
          // Change max index to search lower subarray:
          imax = imid - 1;
        } else {
          // Key found at index imid:
          return imid;
        }
      }
      // Key not found:
      return key_not_found;
    }

    double *values_;                  /*!< Array of non-zero values. */
    double tol_;                      /*!< Array of non-zero values. */
    int *rows_;                       /*!< Array of non-zero values. */
    int *cols_;                       /*!< Array of non-zero values. */
    int memory_chunk_;                /*!< Array of non-zero values. */
    int num_rows_;                    /*!< Array of non-zero values. */
    int num_cols_;                    /*!< Array of non-zero values. */
    int num_elements_;                /*!< Array of non-zero values. */
    int num_non_zero_;                /*!< Array of non-zero values. */
    int adapted_for_mumps_;           /*!< Array of non-zero values. */
    int allocated_memory_values_;     /*!< Array of non-zero values. */
    int allocated_memory_rows_cols_;  /*!< Array of non-zero values. */
};

#endif
