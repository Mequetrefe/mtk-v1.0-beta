/*!
\file mtk_1d_div.h

\brief Includes the definition of the class MTK_1DDiv.

This class implements a 1D DIVERGENCE matrix operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm.

\date: Sunday, September 02, 2012

\version: 2012-09-02.

\author: Eduardo J. Sanchez: esanchez at mail dot sdsu dot edu
*/

/*
Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
State University (SDSU).

Website for the entire project: http://www.csrc.sdsu.edu/mtk/

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

esanchez at mail dot sdsu dot edu

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef MTK_INCLUDE_MTK_1D_DIV_H
#define MTK_INCLUDE_MTK_1D_DIV_H

#include <iostream>
#include <iomanip>

#include "glpk.h"

#include "mtk_roots.h"
#include "mtk_tool_manager.h"
#include "mtk_dense_matrix.h"
#include "mtk_blas_facade.h"
#include "mtk_lapack_facade.h"
#include "mtk_glpk_facade.h"

/****************************
 * These functions should be moved to an outside class
 * They are here temporarily
 * **********************/

inline int idx(const int ii, const int offset, const int jj);

mtk::MTK_Real* Vandermonde(const mtk::MTK_Real *gen,      // Given generator vector.
                    const int gen_length,   // Length of the generator vector.
                    const int pro_length,   // Length of the progression.
                    const bool transpose); // Should create its transpose?
// Prints a dense matrix:

bool DenseMatrix_Print(const mtk::MTK_Real* aa,      // The matrix.
                           const int num_rows,    // Number of rows.
                           const int num_cols);

mtk::MTK_Real CalculateNorm(mtk::MTK_Real *A,int n);


mtk::MTK_Real* Transpose(mtk::MTK_Real *AA, int rr, int cc);

namespace mtk {

/*!
\brief Definition of the MTK_1DDiv class.

This class implements a 1D DIVERGENCE matrix operator, constructed using the
Castillo-Blomgren-Sanchez (CBS) Algorithm.
*/
class MTK_1DDiv {
 public:
  friend MTK_DenseMatrix;
  /*!
  \brief
  */
  friend std::ostream& operator<< (std::ostream& stream,
                                   const MTK_1DDiv& matrix);
  /*!
  \brief
  */
  MTK_1DDiv(void);

  /*!
  \brief
  */
  MTK_1DDiv(int);

  /*!
  \brief
  */
  MTK_1DDiv(int, MTK_Real);
 
  /*!
  \brief
  */
  MTK_1DDiv* Construct1DDiv(void);


  /*!
  \brief
  */
  ~MTK_1DDiv();

  /*!
  \brief
  */
  MTK_DenseMatrix* ReturnAsMatrix(int, MTK_Real, MTK_Real);
  
  MTK_Real* ReturnWeights(void);
  
  MTK_Real* ReturnStencil(void);
  
  MTK_Real* ReturnExtraRows(void);
  MTK_DenseMatrix* ReturnMimeticCoefficients(void);

  
 private:
  /*!
  \brief
  */
  bool ComputeStencilInteriorGrid(void);

  /*!
  \brief
  */
  bool ComputeScaledNullSpace(void);

  /*!
  \brief
  */
  bool ComputePreliminaryApproximations(void);

  /*!
  \brief
  */
  bool ComputeWeights(void);

  /*!
  \brief
  */
  bool ComputeStencilBoundaryGrid(void);

  /*!
  \brief
  */
  bool AssembleOperator(void);

  /*!
  \brief Vandermonde matrix for interior stencil.
  */
 
  char side_;   // Used to recover Q from QR factorization. See Table 1.
  char trans_;  // Should matrix be transposed when solving?
  char transa_; // Is A the transposed in the matrix-matrix multiply?
  char transb_; // Is B the transposed in the matrix-matrix multiply?

  MTK_Real* DD_;  // Output array containing the operator and the weights.

  MTK_Real* dense_values_;  //< Collection of values for dense operators.
  MTK_Real* pp_;        // Spatial coordinates to create interior stencil.
  MTK_Real* oo_;        // Order-selector vector for interior stencil.
  MTK_Real* gg_;        // Spatial coordinates for the boundary points for AA.
  MTK_Real* AA_;        // Boundary Vandermonde to get the systems.
  MTK_Real* ob_;        // Order-selector vector for the boundary stencils.
  MTK_Real* work_;      // Work array for dgels_ solver and dgeqrf_ QR factor.
  MTK_Real* aa_;        // Matrix aa used to compute basis for null-space.
  MTK_Real* tau_;       // Elementary scalar Householder reflectors.
  MTK_Real* QQ_;        // Matrix Q containing QR factorization of matrix
  MTK_Real* KK_  ;      // Extracted basis for the null-space out of matrix Q.
  MTK_Real* SUBK_ ;     // Sub-matrix used to scale KK.
  MTK_Real* SUBKT_  ;   // Transposed sub-matrix used to scale KK.
  MTK_Real* TT_;
  MTK_Real* II_;        // Collection of right-hand sides used to scale KK.
  MTK_Real* KKT_;       // Transpose of the matrix containing KK.
  MTK_Real* NULLS_;     // Scaled basis for the null-space.
  MTK_Real* prem_apps_; // 2D array of boundary preliminary approximations.
  MTK_Real* ob_bottom_; // Last dim_null values of pre-scaled boundary sol.
  MTK_Real* PI_;        // Matrix used to compute the weights.
  MTK_Real* qq_;        // Array containing the weights.
  MTK_Real* qq_lp_;     // Array containing the weights.
  MTK_Real* PIT_;       // Transpose of the PI matrix to solve with LAPACK.
  MTK_Real* lambdas_;   // Collection scalars achieved along with weights.
  MTK_Real* alphas_;    // Collection of scalars to get final operator.
  MTK_Real* mim_bndy_;  // Array containing mimetic boundary approximations.

  MTK_Real mimetic_tol_;  //< Mimetic tolerance.
  MTK_Real mimetic_threshold_;  //< Mimetic tolerance.

  MTK_Real alpha_;    // Scalar in the matrix-matrix multiplication.
  MTK_Real beta_;     // Second scalar in the matrix-matrix multiplication.
  MTK_Real norm_;     //
  MTK_Real minnorm_;  //
  MTK_Real aux_;      //
  MTK_Real normerr_;  //

  int* ipiv_; // Contains information about pivoting.

  int order_; //< Selected order of accuracy.

  int kk_;                //< Order of numerical accuracy of the operator.
  int nrhs_;              // Number of right-hand sides to solve for.
  int TT_ld_;             // Leading dimension of the Vandermonde matrix.
  int ldoo_;              // Leading dimension of the order-selector vector.
  int dim_null_;          // Dimension null-space for boundary approximations.
  int num_bndy_approxs_;  // Required boundary points uniform order accuracy.
  int AA_num_rows_;       // Number rows for the boundary Vandermonde matrices.
  int AA_num_cols_;       // Number cols for the boundary Vandermonde matrices.
  int AA_ld_;             // Leading dimension of boundary Vandermonde matrix.
  int ob_ld_;             // Leading dimension order-selector vector boundary.
  int info_;              // Information regarding solution of the system.
  int lwork_;             // Length of work array for the dgels_ solver.
  int aa_num_rows_;       // Number of rows of matrix aa.
  int aa_num_cols_;       // Number of columns of matrix aa.
  int aaT_num_rows_;      // Number of rows of matrix aa transposed.
  int aaT_num_cols_;      // Number of columns of matrix aa transposed.
  int ltau_;              // Length of the tau array.
  int KK_num_rows_;       // Number of rows of the KK matrix.
  int KK_num_cols_;       // Number of columns of the KK matrix.
  int ldSUBKT_;           // Leading dimensions of the SUBKT matrix.
  int incx_;              // Specifies increment for the elements of X on DAXPY.
  int PI_num_cols_;       // Number of columns of the PI matrix.
  int lDD_;               // Length of the output array.
  int minrow_;            // Row from the optimizer with the minimum rel. nor.
  int row_;               // Row currently processed by the optimizer
};

}
#endif
