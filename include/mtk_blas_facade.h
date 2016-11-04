/*!
\file mtk_blas_adapter.h

\brief Adapter class for the BLAS API.

This class contains a collection of static classes, that posses direct access
to the underlying structure of the matrices, thus allowing programmers to
exploit some of the numerical methods implemented in the BLAS.

The BLAS (Basic Linear Algebra Subprograms) are routines that provide standard
building blocks for performing basic vector and matrix operations. The Level 1
BLAS perform scalar, vector and vector-vector operations, the Level 2 BLAS
perform matrix-vector operations, and the Level 3 BLAS perform matrix-matrix
operations.

\sa http://www.netlib.org/blas/

\date: July 1, 2014, 5:00 PM

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

#ifndef MTK_INCLUDE_MTK_BLAS_ADAPTER_H
#define	MTK_INCLUDE_MTK_BLAS_ADAPTER_H

namespace mtk {
  
extern "C" {
    
// External from the BLAS:

// The BLAS can be installed from either:
// 1. http://www.netlib.org/blas/
// 2. https://software.intel.com/en-us/non-commercial-software-development

// Double precision General Matrix-Vector Multiplier:
// Performs one of the following matrix-vector operations:
// y := alpha*A*x + beta*y, or
// y := alpha*A'*x + beta*y.
// http://www.math.utah.edu/software/lapack/lapack-blas/dgemv.html

void dgemv_(char *trans,    // 1. Are the matrices transposed?
            int *m,         // 2.
            int *n,         // 3.
            double *alpha,  // 4.
            double *a,      // 5.
            int *lda,       // 6.
            double *x,      // 7.
            int *incx,      // 8.
            double *beta,   // 9.
            double *y,      // 10.
            int *incy);     // 11.

// Double-Precision General Matrix-Matrix multiplier:
// http://www.math.utah.edu/software/lapack/lapack-blas/dgemm.html

void dgemm_(char *transa,   // 1. Is the a matrix transposed?
            char* transb,   // 2. Is the b matrix transposed?
            int *m,         // 3.
            int *n,         // 4.
            int *k,         // 5.
            double *alpha,  // 6. Scalar modifying the first matrix.
            double *a,      // 7. First matrix.
            int *lda,       // 8. Leading dimension of the a matrix.
            double *b,      // 9. Second matrix.
            int *ldb,       // 10. Leading dimension of the b matrix.
            double *beta,   // 11. Scalar modifying the second matrix.
            double *c,      // 12. Third and output matrix.
            int *ldc);      // 13. Leading dimension of the c matrix.
}

/*!
\class MTK_BLASAdapter

\ingroup c07-num_methods

\brief Adapter class for the BLAS API.

This class contains a collection of static classes, that posses direct access
to the underlying structure of the matrices, thus allowing programmers to
exploit some of the numerical methods implemented in the BLAS.

The BLAS (Basic Linear Algebra Subprograms) are routines that provide standard
building blocks for performing basic vector and matrix operations. The Level 1
BLAS perform scalar, vector and vector-vector operations, the Level 2 BLAS
perform matrix-vector operations, and the Level 3 BLAS perform matrix-matrix
operations.

\sa http://www.netlib.org/blas/
*/
class MTK_BLASAdapter {

    /*! Default solver.  */
  /*! This functions solves the wave equation in 2D. */
  public:
    static void ScalarVectorProduct(double scalar, unsigned int nn,
                                    double vv[], double zz[]);
  
};
}

#endif	// End of: MTK_INCLUDE_MTK_BLAS_ADAPTER_H
