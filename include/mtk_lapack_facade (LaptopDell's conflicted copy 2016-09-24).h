/*!
\file mtk_lapack_facade.h

\brief Adapter class for the LAPACK API.

This class contains a collection of static classes, that posses direct access
to the underlying structure of the matrices, thus allowing programmers to
exploit some of the numerical methods implemented in the LAPACK.

The LAPACK is written in Fortran 90 and provides routines for solving systems of
simultaneous linear equations, least-squares solutions of linear systems of
equations, eigenvalue problems, and singular value problems.

\sa http://www.netlib.org/lapack/

\date: June 30, 2014, 11:33 AM

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

#ifndef MTK_INCLUDE_MTK_LAPACK_FACADE_H
#define	MTK_INCLUDE_MTK_LAPACK_FACADE_H

#include "mtk_roots.h"

namespace mtk {

class MTK_DenseMatrix;

extern "C" {
  

void dgels_(char* trans,  // 1. Am I giving the transpose of the matrix?
            int* m,       // 2. The number of rows of the matrix a.  m >= 0.
            int* n,       // 3. The number of columns of the matrix a.  n >= 0.
            int* nrhs,    // 4. The number of right-hand sides.
            MTK_Real* a,    // 5. On entry, the m-by-n matrix a.
            int* lda,     // 6. The leading dimension of a. lda >= max(1,m).
            MTK_Real* b,    // 7. On entry, matrix b of right-hand side vectors.
            int* ldb,     // 8. The leading dimension of b. ldb >= max(1,m,n).
            MTK_Real* work, // 9. On exit, if info = 0, work(1) is optimal lwork.
            int* lwork,   // 10. The dimension of the array work.
            int* info);   // 11. If info = 0, then successful exit.

// Double-Precision General QR Factorization to compute the null-space:
// Theory: http://www.netlib.org/lapack/lug/node69.html
// Documentation: http://www.math.utah.edu/software/lapack/lapack-d/dgeqrf.html
// Double-Precision Orthogonal Make Q from QR:
// dormqr_ overwrites the general real M-by-N matrix C with (Table 1):
//
//                 SIDE = 'L'     SIDE = 'R'
// TRANS = 'N':      Q * C          C * Q
// TRANS = 'T':      Q**T * C       C * Q**T
//
// where Q is a real orthogonal matrix defined as the product of k
// elementary reflectors
//
//       Q = H(1) H(2) . . . H(k)
//
// as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
// if SIDE = 'R'.
// Documentation: http://www.netlib.no/netlib/lapack/MTK_Real/dormqr.f



/*!
External from the LAPACK:

The LAPACK can be installed from either:
1. http://www.netlib.org/lapack/
2. https://software.intel.com/en-us/non-commercial-software-development
*/

/*!
\brief Double-Precision General Matrix Factorization & Multiple Right-Hand Side.

Solver for the Vandermonde system approximating at the interior of the grid.
http://www.netlib.org/lapack/explore-html/d8/d72/dgesv_8f.html

\param[in]   n     The order of the matrix A.
\param[in]   nrhs  The number of rows of the matrix a.  m >= 0.
\param[in,out]  a     On entry, the n-by-n matrix a.
\param[in]   lda   The leading dimension of a. lda >= max(1,m).
\param[in,out]  ipiv  The pivots that define the permutation matrix P.
\param[in,out]  b     On entry, matrix b of right-hand side vectors.
\param[in]   ldb   The leading dimension of b. ldb >= max(1,m,n).
\param[in,out]  info  If info = 0, then successful exit.
}
*/
void dgesv_(int* n,
            int* nrhs,
            MTK_Real* a,
            int* lda,
            int* ipiv,
            MTK_Real* b,
            int* ldb,
            int* info);

/*!
\brief Double-Precision General QR Factorization to compute the null-space.

Double-Precision General QR Factorization to compute the null-space:
Theory: http://www.netlib.org/lapack/lug/node69.html
Documentation: http://www.math.utah.edu/software/lapack/lapack-d/dgeqrf.html

\param[i]   m     The number of columns of the matrix a.  n >= 0.
\param[i]   n     The number of columns of the matrix a.  n >= 0.
\param[io]  a     On entry, the n-by-n matrix a.
\param[i]   lda   Leading dimension matrix.  LDA >= max(1,M).
\param[io]  tau   Scalars from elementary reflectors. min(M,N).
\param[io]  work  Workspace. info = 0, work(1) is optimal lwork.
\param[i]   lwork The dimension of work. lwork >= max(1,n).
\param[i]   info  info = 0: successful exit.
*/
void dgeqrf_(int *m,
             int *n,
             MTK_Real *a,
             int *lda,
             MTK_Real *tau,
             MTK_Real *work,
             int *lwork,
             int *info);

/*!
\brief Double-Precision Orthogonal Make Q from QR.

Double-Precision Orthogonal Make Q from QR:
dormqr_ overwrites the general real M-by-N matrix C with (Table 1):

                SIDE = 'L'     SIDE = 'R'
TRANS = 'N':      Q * C          C * Q
TRANS = 'T':      Q**T * C       C * Q**T

where Q is a real orthogonal matrix defined as the product of k
elementary reflectors

      Q = H(1) H(2) . . . H(k)

as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
if SIDE = 'R'.
Documentation: http://www.netlib.no/netlib/lapack/MTK_Real/dormqr.f

\param[i]   side  See Table 1 above.
\param[i]   trans See Table 1 above.
\param[i]   m     Number of rows of the C matrix.
\param[i]   n     Number of columns of the C matrix.
\param[i]   k     Number of reflectors.
\param[io]  a     The matrix containing the reflectors.
\param[i]   lda   The dimension of work. lwork >= max(1,n).
\param[i]   tau   Scalar factors of the elementary reflectors.
\param[i]   c     Output matrix.
\param[i]   ldc   Leading dimension of the output matrix.
\param[io]  work  Workspace. info = 0, work(1) optimal lwork.
\param[i]   lwork The dimension of work.
\param[io]  info  info = 0: successful exit.
*/
void dormqr_(char *side,
             char *trans,
             int *m,
             int *n,
             int *k,
             MTK_Real *a,
             int *lda,
             MTK_Real *tau,
             MTK_Real *c,
             int *ldc,
             MTK_Real *work,
             int *lwork,
             int *info);
}

/*!
\class MTK_LAPACKFacade

\ingroup c07-num_methods

\brief Adapter class for the LAPACK API

This class contains a collection of static classes, that posses direct access
to the underlying structure of the matrices, thus allowing programmers to
exploit the numerical methods implemented in the LAPACK.

The LAPACK is written in Fortran 90 and provides routines for solving systems
of simultaneous linear equations, least-squares solutions of linear systems of
equations, eigenvalue problems, and singular value problems.

\sa http://www.netlib.org/lapack/
*/
class MTK_LAPACKFacade {

  public:
    /*!
    \brief Solves a dense system of linear equations.

    Adapts the MTK to LAPACK's dgesv_ routine.

    \param[i] matrix Input matrix.
    \param[i] rhs Input right-hand side vector.
    \param[i] nrhs Number of input right-hand side vectors.

    \exception std::bad_alloc
    */
    static void SolveDenseSystem(MTK_DenseMatrix* matrix,
                                 MTK_Real *rhs,
                                 int nrhs);

    /*!
    \brief Performs a QR factorization on a dense matrix.

    Adapts the MTK to LAPACK's  routine.

    \param[i] matrix Input matrix.
    \param[i] rhs Input right-hand side vector.
    \param[i] nrhs Number of input right-hand side vectors.

    \return Success of the solution.

    \exception std::bad_alloc
    */
    static MTK_Real* QRFactorDenseMatrix(MTK_Real* aa,
                                         int aa_num_rows,
                                         int aa_num_cols,
                                         int aaT_num_rows,
                                         int aaT_num_cols,
                                         int num_bndy_approxs);
};
}

#endif	// End of: MTK_INCLUDE_MTK_LAPACK_ADAPTER_H
