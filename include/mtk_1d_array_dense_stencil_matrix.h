/*!
  \file mtk_1d_array_dense_stencil_matrix.h

  \brief Includes the definition of the class MTK_1DArrayDenseStencilMatrix.

  This class implements a 1D Stencil matrix \f$ M \f$ which is then used to solve
  for the problem by means of solving \f$ M\mathbf{x}=\mathbf{s} \f$.

  \date: Monday, September 03, 2012
  \version: 2011-09-03.
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

#ifndef MTK_INCLUDE_MTK_1D_ARRAY_DENSE_STENCIL_MATRIX_H
#define MTK_INCLUDE_MTK_1D_ARRAY_DENSE_STENCIL_MATRIX_H

class MTK_1DCGMDirichletOperator;
class MTK_1DCGMNeumannOperator;
class MTK_1DCGMGradient;
class MTK_1DCGMLaplacian;

/*! \brief Includes the definition of the class MTK_1DArrayDenseStencilMatrix.
 *
 *   This class implements a 1D Stencil matrix \f$ M \f$ which is then used to
 * solve for the problem by means of solving \f$ M\mathbf{x}=\mathbf{s} \f$.
 * \todo Monday, September 03, 2012: Change the representation to an enum.
 * \todo Monday, September 03, 2012: Add the sparse connection.
 */
class MTK_1DArrayDenseStencilMatrix {

  public:
    /*! Default constructor.  */
    /*! This functions CONSTRUCTS a default stencil matrix. */
    MTK_1DArrayDenseStencilMatrix(void);

    /*! Constructor #2. */
    /*! Constructs a matrix using the fact: \f$ M=\alpha A+\beta BG-L \f$.*/
    MTK_1DArrayDenseStencilMatrix(double alpha, MTK_1DCGMDirichletOperator *AA,
                       double beta,  MTK_1DCGMNeumannOperator *BB,
                       MTK_1DCGMGradient *GG, MTK_1DCGMLaplacian *LL);

    /*! Constructor #3.  */
    /*! This functions CONSTRUCTS a stencil matrix based on common dims. */
    MTK_1DArrayDenseStencilMatrix(int nrows, int ncols);

    /*! Prints the matrix. */
    /*! Prints the matrix.
     */
    void Print(int display);

    /*! Default destructor.  */
    /*! This functions DESTRUCTS an already created stencil matrix. */
    ~MTK_1DArrayDenseStencilMatrix() {
      delete [] dense_values_;
    }

    /*! Sets a specific value.  */
    /*! Sets a specific value.  */
    void SetValue(int ii, int jj, double xx) {

      dense_values_[ii*size_ + jj] = xx;
    }

    /*! Gets the collection of non-zero values.  */
    /*! Gets the collection of non-zero values. */
    double* dense_values(void) { return dense_values_; }

    /*! Gets the type of representation.  */
    /*! Gets the type of representation. */
    int representation(void) { return representation_; }

    /*! Gets the amount of both 0s and n0s elements withing the matrix.  */
    /*! Gets the amount of both 0s and n0s elements withing the matrix. */
    int size(void) { return size_; }

    /*! Gets the number of rows of the matrix.  */
    /*! Gets the number of rows of the matrix. */
    int num_rows(void) { return num_rows_; }

    /*! Gets the number of columns of the matrix.  */
    /*! Gets the number of columns of the matrix. */
    int num_cols(void) { return num_cols_; }

    /*! Gets a value form within the matrix.  */
    /*! Gets a value form within the matrix. */
    double GetValue(int ii, int jj) {

      return dense_values_[ii*num_cols_ + jj];
    }

  private:
    double *dense_values_;  /*! < Collection of values for dense operators. */
    int representation_;    /*! < Selected storage option. */
    int order_;             /*! < Selected order of accuracy.*/
    int size_;              /*! < Number of nodes I will be acting upon. */
    int num_rows_;          /*! < Number of rows. */
    int num_cols_;          /*! < Number of columns. */
};

#endif
