/*!
\file mtk_dense_matrix.h

\brief Defines a common dense matrix, using a 1-D array.

For developing purposes, it is better to have a not-so-intrincated data
structure implementing matrices. This is the purpose of this class: to be used
for prototypes of new code for small test cases. In every other instance, this
should be replaced by the most appropriate sparse matrix.

\date: June 20, 2014, 12:00 PM

\version: 1.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at sciences dot sdsu dot edu

\bug No known bugs.
*/

/*
Copyright (C) 2015 Computational Science Research Center (CSRC) at San Diego
State University (SDSU).

Website for the project: http://www.csrc.sdsu.edu/mtk/

All rights reserved.

Redistribution and use in source and binary forms, with or without modification
are permitted provided that the following conditions are met:

-# Redistributions of source code must retain the above copyright notice, this
list of conditions and this disclaimer.
-# Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.
-# Neither the name of the CSRC, SDSU nor the names of its contributors may be
used to endorse or promote products derived from this software without specific
prior written permission.
-# Modifications whether they are partial or complete; whether they are
additions or eliminations should be reported through an email at:

esanchez at sciences dot sdsu dot edu

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef MTK_INCLUDE_MTK_DENSE_MATRIX_H
#define MTK_INCLUDE_MTK_DENSE_MATRIX_H

#include "mtk_roots.h"
#include "mtk_matrix.h"
#include "mtk_lapack_facade.h"


namespace mtk {

/*!
\class MTK_DenseMatrix

\ingroup c03-data_structures

\brief Defines a common dense matrix, using a 1-D array.

For developing purposes, it is better to have a not-so-intrincated data
structure implementing matrices. This is the purpose of this class: to be used
for prototypes of new code for small test cases. In every other instance, this
should be replaced by the most appropriate sparse matrix.
*/
class MTK_DenseMatrix {
  friend class MTK_LAPACKFacade;
 public:
  /*!
  \brief Default constructor.
  */
  MTK_DenseMatrix();

  /*!
  \brief Construct a dense matrix based on the given dimensions.
  \param[i] num_rows Number of rows of the required matrix.
  \param[i] num_cols Number of rows of the required matrix.
  \exception std::bad_alloc
  */
  MTK_DenseMatrix(const int num_rows, const int num_cols);

  /*!
  \brief Construct a zero-rows-padded identity matrix.
  \param[i] rank Rank or number of rows/cols in square matrix.
  \exception std::bad_alloc
  */
// MTK_DenseMatrix(const int &rank);

  /*!
  \brief Construct a dense matrix based on the given number of unknowns.
  \param[in] num_unknowns Number of unknowns to solve with the required matrix.
  \param[in] rank Rank or number of rows/cols in square matrix.
  \exception std::bad_alloc
  */
  //MTK_DenseMatrix(const int &num_unknowns, const int &rank);

  /*!
  \brief Construct a dense Vandermonde matrix.

  DEF. In linear algebra, a VANDERMONDE MATRIX is a matrix with terms of a
  geometric progression in each row. This progression uses the terms of a
  given GENERATOR VECTOR.

  This constructor generates a Vandermonde matrix, as defined above.

  \param[in] gen Given generator vector.
  \param[in] gen_length Length generator vector.
  \param[in] pro_length Length the progression.
  \param[in] transpose Should the transpose be created instead?

  \exception std::bad_alloc
  */
  MTK_DenseMatrix(const MTK_Real *gen,
                  const int gen_length,
                  const int pro_length,
                  const bool transpose);

  /*!
  \brief Construct a dense matrix based on the Kronecker product of arguments.
  \param[i] aa First matrix.
  \param[i] bb Second matrix.
  \exception std::bad_alloc
  */
  MTK_DenseMatrix(const MTK_DenseMatrix &aa, const MTK_DenseMatrix &bb);

  /*!
  \brief Destructor.
  */
  ~MTK_DenseMatrix();
  
  MTK_Real GetValue(int rr, int cc) ;

  /*!
  \brief Gets the number of rows.

  \return Number of rows of the matrix.
  */
  int num_rows() const;

  /*!
  \brief Gets the number of columns.

  \return Number of columns of the matrix.
  */
  int num_cols() const;

  /*!
  \brief Sets a value on the given coordinates.

  \param[in] row_coord Row coordinate.
  \param[in] col_coord Column coordinate.
  \param[in] val Row Actual value to be inserted.
  */
  void SetValue(int row_coord, int col_coord, MTK_Real val);

  /*!
  \brief Provides access to the matrix data.

  \return Pointer to a MTK_Matrix.
  */
  MTK_Matrix* matrix_data();

  /*!
  \brief Gets a value on the given coordinates.

  \param[in] row_coord Row coordinate.
  \param[in] col_coord Column coordinate.

  \return The required value at the specified coordinates.
  */
  MTK_Real GetValue(int row_coord, int col_coord) const;

  /*!
  \brief Prints the matrix as a block of number (standard way).
  */
  void PrintBlock(void);

  /*!
  \brief Transposes a given matrix.

  \param[in] m_in Input matrix.

  \return The transpose of m_in.

  \exception std::bad_alloc
  */
  static MTK_DenseMatrix* Transpose(const MTK_DenseMatrix* m_in);

 private:
  /*!
  \brief Used to compute a linear index.

  \todo Analyze how needed this function really is.
  */
  int idx(const int ii, const int offset, const int jj) const {

    return ii*offset + jj;
  }

  MTK_Matrix *matrix_data_; //!< Data related to the matrix nature.
  MTK_Real **Values_;       //!< Array of values.
  MTK_Real *continuum_;     //!< Array continuous memory positions.
  int nRow_;
  int nCol_;
  int number_values_;
  MTK_Real *MemoryHolder_;
};

}

#endif  // End of: MTK_INCLUDE_MTK_DENSE_MATRIX_H
