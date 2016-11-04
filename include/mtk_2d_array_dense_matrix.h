/*!
  \file mtk_2d_array_dense_matrix.h

  \brief Includes the definition of the class MTK_2DArrayDenseMatrix.

  This file contains the definition of the class MTKCDenseMatrix.
  This class implements a sparse matrix using Compressed Column Scheme (CCS) to
  ensure storing only the non-zero elements of a given dense matrix. As a
  consequence of the selected storage scheme, no assumptions are being made in
  terms of the actual structure of the matrix.

  \date: Sunday, August 14 2011 11:42 PM
  \version: 2011-14-08-01.
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

#ifndef MTK_INCLUDE_MTK_DENSE_MATRIX_H
#define MTK_INCLUDE_MTK_DENSE_MATRIX_H

#include "mtk_roots.h"
namespace mtk{

/*! \brief Compatible class MTKCSSparseMatrixDouble.

  This class implements the Compressed Row Scheme (CCS) format for sparse
  matrices. This is a general data type that makes no assumptions about the
  sparsity of the original matrix. This is also compatible with well known
  sparse matrices implementations such as the one in
  <a href="http://crd.lbl.gov/~xiaoye/SuperLU/">SuperLU</a>.
 */
class MTK_2DArrayDenseMatrix {

public:
  /*! Default constructor.  */
  /*! This functions CONSTRUCTS a default CRS sparse matrix. */
  MTK_2DArrayDenseMatrix(void);

  /*! Constructs a CRS sparse matrix form a 1D array.  */
  /*! This functions CONSTRUCTS a CRS sparse matrix from an already predefined
   * matrix which is stored as an one dimensional array.
   * \param *array Input pointer to the array which hold the original matrix.
   * \param num_rows Number of rows of the original matrix.
   * \param num_cols Number of rows of the original matrix.
   * \return The status of the creation as a logical value.
   */
  MTK_2DArrayDenseMatrix(int num_rows,int num_cols);

  /*! Default destructor.  */
  /*! This functions DESTRUCTS an already created dense matrix. */
  bool SetValue (MTK_Real xx, int rr, int cc) {

    if (0 > rr || rr > number_values_ || 0 > cc || cc > number_values_) {
      return false;
    } else {
      Values_[rr][cc] = xx;
      return true;
    }
  }

  /*! Default destructor.  */
  /*! This functions DESTRUCTS an already created dense matrix. */
  int nRow() {

    return nRow_;
  }

  /*! Default destructor.  */
  /*! This functions DESTRUCTS an already created dense matrix. */
  int nCol() {

    return nCol_;
  }

  /*! Default destructor.  */
  /*! This functions DESTRUCTS an already created dense matrix. */
  MTK_Real GetValue(int rr, int cc);

  /*! Default destructor.  */
  /*! This functions DESTRUCTS an already created dense matrix. */
  void PrintBlock();

  /*! Default destructor.  */
  /*! This functions DESTRUCTS an already created dense matrix. */
  ~MTK_2DArrayDenseMatrix();

private:
  MTK_Real **Values_;      /*!< Array of non-zero values. */
  MTK_Real *MemoryHolder_; /*!< Array of column indexes. */
  int nRow_;                /*!< Number of rows. */
  int nCol_;                /*!< Number of columns. */
  int number_values_;
};

}
#endif
