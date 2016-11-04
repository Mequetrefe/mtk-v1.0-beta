/*!
  \file mtk_crs_sparse_matrix.h

  \brief Includes the definition of the class MTK_CRSSparseMatrix.

  This file contains the definition of the class MTK_CRSSparseMatrix.
  This class implements a sparse matrix using Compressed Row Scheme (CRS) to
  ensure storing only the non-zero elements of a given dense matrix. As a
  consequence of the selected storage scheme, no assumptions are being made in
  terms of the actual structure of the matrix.

  \date: Thursday, October 27 2011, 09:16 PM
  \version: 2011-10-27-01.

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

#ifndef MTK_INCLUDE_MTK_CRS_SPARSE_MATRIX_H
#define MTK_INCLUDE_MTK_CRS_SPARSE_MATRIX_H

/*! \brief Compatible class MTKCRSparseMatrix.

  This class implements the Compressed Row Scheme (CRS) for sparse matrices.
  This is a general data type that makes no assumptions about the sparsity of
  the original matrix. This is also compatible with well known sparse matrices
  implementations such as the one in
  <a href="http://crd.lbl.gov/~xiaoye/SuperLU/">SuperLU</a>.
*/
class MTK_CRSSparseMatrix {

  public:
    /*! Default constructor.  */
    /*! This functions CONSTRUCTS a default CRS sparse matrix. */
    MTK_CRSSparseMatrix(void);

    /*! Constructs a CRS sparse matrix form a 1D array.  */
    /*! This functions CONSTRUCTS a CRS sparse matrix from an already predefined
        matrix which is stored as an one dimensional array.
      \param *array Input pointer to the array which hold the original matrix.
      \param num_rows Number of rows of the original matrix.
      \param num_cols Number of rows of the original matrix.
      \return The status of the creation as a logical value.
    */
    MTK_CRSSparseMatrix(double *array,
                            int num_rows,
                            int num_cols,
                            int unitary_offset);

    /*! Default destructor.  */
    /*! This functions DESTRUCTS an already created CRS sparse matrix. */
    ~MTK_CRSSparseMatrix() {
      delete [] values;
      delete [] col_ind;
      delete [] row_ptr;
    }

    /*! Prints a CRS sparse matrix in standard matrix algebraic notation.  */
    /*! This procedure prints a CRS sparse matrix in standard algebraic (or
      block) notation.
      \param *in Input pointer to the matrix which one wants to print.
      \param dots Prints only the sparse structure.
    */
    void MTKBlockPrint(int dots);

    /*! Prints a CRS sparse matrix as it is stored: array by array.  */
    /*! This procedure prints the collection of arrays that conform the data
      structure used to represent a CRS sparse matrix.
      \param *in Input pointer to the matrix which one wants to print.
    */
    void MTKArrayPrint(void);

  private:
    double *values;       /*!< Array of non-zero values. */
    double *col_ind;      /*!< Array of column indexes. */
    double *row_ptr;      /*!< Array of values that start a row. */
    int num_rows;         /*!< Number of rows. */
    int num_cols;         /*!< Number of columns. */
    int number_values;    /*!< Amount of values in the matrix. */
    int number_non_zero;  /*!< Number of non-zero values. */
    int number_zero;      /*!< Number of non-zero values. */
    int unitary_offset;   /*!< Should position 0 be position 1? */
};

#endif
