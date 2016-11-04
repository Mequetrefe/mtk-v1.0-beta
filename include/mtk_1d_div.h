/*!
\file mtk_1d_div.h

\brief Includes the definition of the class MTK_1DDiv.

This class implements a 1D DIVERGENCE matrix operator, based on the
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

#include "mtk_roots.h"

/*! \brief Definition of the MTK_1DDiv class.

This class implements a 1D DIVERGENCE matrix operator, based on the
Castillo-Blomgren-Sanchez (CBS) Algorithm.
*/
class MTK_1DDiv {

  public:

      friend std::ostream &operator<<(std::ostream &output,
                                  const MTK_1DDiv &in ) {

//     output << in.XXX;

    return output;
  }


  // Prints a mimetic operator as a collection of arrays:

void MTK_PrintOperatorArrays(void) {

  cout << "Constructed operator DD:" << endl;
  cout << "DD[0] = " << setw(9) << DD[0] << endl;
  cout << "DD[1:" << kk << "] = ";
  for (auto ii = 1; ii <= kk; ii++) {
    cout << setw(9) << DD[ii] << " ";
  }
  cout << endl;
  if (kk > 2) {
    cout << "DD[" << kk + 1 << ":" << 2*kk << "] = ";
    for (auto ii = kk + 1; ii <= 2*kk; ii++) {
      cout << setw(9) << DD[ii] << " ";
    }
    cout << endl;

    auto offset = (2*kk + 1);
    int mm {};
    for (auto ii = 0; ii < dim_null; ii++) {
      cout << "DD[" << offset + mm << ":" <<
        offset + mm + num_bndy_approxs - 1 << "] = ";
      for (auto jj = 0; jj < num_bndy_approxs; jj++) {
        auto value = DD[offset + (mm)];
        cout << setw(9) << ((abs(value - 0.0) < mtk_tol)? 0.0: value) << " ";
        mm++;
      }
      cout << endl;
    }
  }
  cout << endl;
}


    /*! Default constructor.  */
    /*! This functions CONSTRUCTS a default Divergence operator. */
    MTK_1DDiv(void);

    /*! Constructor #2. */
    /*! Constructs an operator starting from the usual parameters. */
    MTK_1DDiv(const MTK_1DDiv&);

    MTK_1DDiv(int);

    MTK_1DDiv(int,double);

    /*! Default destructor.  */
    /*! This functions DESTRUCTS an already created Divergence operator. */
    ~MTK_1DCGMDivergence() {
      delete [] dense_values_;
    }

  private:
    double
optimizer (double *A, int nrows, int ncols, int kk, double *hh, double *qq, int robjective, double mimetic_tolerance, int copy);

// Stage 1 of the CRS Algorithm:
// Compute the stencil approximating the interior of the staggered grid:

bool ComputeStencilInteriorGrid(void);

// Stage 2.1 of the CRS Algorithm:
// Compute a scaled basis of the null-space of the Vandermonde matrix
// approximating at the west boundary:

bool ComputeScaledNullSpace(void);

// Stage 2.2 of the CRS Algorithm:
// Compute the set of preliminary approximation on the boundary neighborhood:

bool ComputePreliminaryApproximations(void);

// Stage 2.3 of the CRS Algorithm:
// Assemble the PI matrix and compute the weights:

bool ComputeWeights(void);

// Stage 2.4 of the CRS Algorithm:
// Compute mimetic stencil approximating at boundary and assemble operator:

bool ComputeStencilBoundaryGrid(void);

// Stage 3: Final stage. Construct the output array with the operator and its
// weights:

bool AssembleOperator(void);





    double *dense_values_;  /*! < Collection of values for dense operators. */
    int representation_;    /*! < Selected storage option. */
    int order_;             /*! < Selected order of accuracy.*/
    int size_;              /*! < Number of nodes I will be acting upon. */


    // Input:
int kk {}; // Order of numerical accuracy of the operator.

// Output:
double* DD {};  // Output array containing the operator and the weights.

// Variables:
char side {'L'};    // Used to recover Q from QR factorization. See Table 1.
char trans {'N'};   // Should matrix be transposed when solving?
char transa {'N'};  // Is A the transposed in the matrix-matrix multiply?
char transb {'N'};  // Is B the transposed in the matrix-matrix multiply?

double* pp {};        // Spatial coordinates to create interior stencil.
double* TT {};        // Vandermon  de matrix for interior stencil.
double* oo {};        // Order-selector vector for interior stencil.
double* gg {};        // Spatial coordinates for the boundary points for AA.
double* AA {};        // Boundary Vandermonde to get the systems.
double* ob {};        // Order-selector vector for the boundary stencils.
double* work {};      // Work array for dgels_ solver and dgeqrf_ QR factor.
double* aa {};        // Matrix aa used to compute the basis for the null-space.
double* tau {};       // Elementary scalar Householder reflectors.
double* QQ {};        // Matrix Q containing the QR factorization of matrix aa.
double* KK {};        // Extracted basis for the null-space out of matrix Q.
double* SUBK {};      // Sub-matrix used to scale KK.
double* SUBKT {};     // Transposed sub-matrix used to scale KK.
double* II {};        // Collection of right-hand sides used to scale KK.
double* KKT {};       // Transpose of the matrix containing KK.
double* NULLS {};     // Scaled basis for the null-space.
double* prem_apps {}; // 2D array of boundary preliminary approximations.
double* ob_bottom {}; // Last dim_null values of the pre-scaled boundary sol.
double* PI {};        // Matrix used to compute the weights.
double* qq {};        // Array containing the weights.
double* qq_lp {};        // Array containing the weights.
double* PIT {};       // Transpose of the PI matrix to solve with LAPACK.
double* lambdas {};   // Collection of scalars achieved along with the weights.
double* alphas {};    // Collection of scalars to get final operator.
double* mim_bndy {};  // Array containing the mimetic boundary approximations.

double alpha {1.0}; // Scalar in the matrix-matrix multiplication.
double beta {0.0};  // Second scalar in the matrix-matrix multiplication.

int* ipiv {}; // Contains information about pivoting.

int nrhs {1};         // Number of right-hand sides to solve for.
int TT_ld {};       // Leading dimension of the Vandermonde matrix.
int ldoo {};        // Leading dimension of the order-selector vector.
int dim_null;         // Dimension of null-space for boundary approximations.
int num_bndy_approxs; // Required boundary points for uniform order accuracy.
int AA_num_rows;      // Number rows for the boundary Vandermonde matrices.
int AA_num_cols;      // Number cols for the boundary Vandermonde matrices.
int AA_ld;            // Leading dimension of the boundary Vandermonde matrix.
int ob_ld;            // Leading dimension order-selector vector for boundary.
int info {0};         // Information regarding solution of the system.
int lwork {-1};       // Length of work array for the dgels_ solver.
int aa_num_rows;      // Number of rows of matrix aa.
int aa_num_cols;      // Number of columns of matrix aa.
int aaT_num_rows;     // Number of rows of matrix aa transposed.
int aaT_num_cols;     // Number of columns of matrix aa transposed.
int ltau;             // Length of the tau array.
int KK_num_rows;      // Number of rows of the KK matrix.
int KK_num_cols;      // Number of columns of the KK matrix.
int ldSUBKT;          // Leading dimensions of the SUBKT matrix.
int incx {1};         // Specifies increment for the elements of X on DAXPY.
int PI_num_cols;      // Number of columns of the PI matrix.
int lDD;              // Length of the output array.

double norm;
double minnorm;
int minrow;
int row;
double aux;
double normerr;


};

#endif
