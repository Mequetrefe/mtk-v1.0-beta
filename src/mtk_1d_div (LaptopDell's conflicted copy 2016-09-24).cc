/*!
\file mtk_1d_div.cc

\brief Constructs a default MTK_1DDiv operator.

Constructs a default MTK_1DDiv. Note how the default order is at least 2.

\date: Thursday, June 03, 2012
\version: 2012-09-03.
\author: Raul Vargas Navarro: rvargasnavarro at nu.edu
*/

/*
Copyright (C) 2015 Computational Science Research Center (CSRC) at San Diego
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

esanchez at sciences.sdsu.edu

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>

#include "mtk_roots.h"

#include "mtk_1d_div.h"

using namespace std;
using namespace mtk;

/*! \brief Constructor #1. */

/*! Constructs an operator starting from default parameters.
*/
MTK_1DDiv::MTK_1DDiv() {

  side_ = 'L';
  trans_ = 'N';
  transa_ = 'N';
  transb_ = 'N';

  DD_ = nullptr;

  TT_=nullptr;
  pp_= nullptr;
  oo_ = nullptr;
  gg_ = nullptr;
  AA_ = nullptr;
  ob_ = nullptr;
  work_ = nullptr;
  aa_ = nullptr;
  tau_ = nullptr;
  QQ_ = nullptr;
  KK_ = nullptr;
  SUBK_ = nullptr;
  SUBKT_ = nullptr;
  II_ = nullptr;
  KKT_ = nullptr;
  NULLS_ = nullptr;
  prem_apps_ = nullptr;
  ob_bottom_ = nullptr;
  PI_ = nullptr;
  qq_ = nullptr;
  qq_lp_ = nullptr;
  PIT_ = nullptr;
  lambdas_ = nullptr;
  alphas_ = nullptr;
  mim_bndy_ = nullptr;

  mimetic_tol_ = mtk::MTK_DefaultMimTol;
  mimetic_threshold_ = mtk::MTK_DefaultMimeticThreshold;

  alpha_ = 1.0;
  beta_  = 0.0;
  norm_ = 0.0;
  minnorm_ = 0.0;
  aux_ = 0.0;
  normerr_ = 0.0;

  ipiv_ = nullptr;

  order_ = mtk::MTK_DefaultOrder;

  kk_ = mtk::MTK_DefaultOrder;
  nrhs_ = 1;
  TT_ld_ = 0;
  ldoo_ = 0;
  dim_null_ = 0;
  num_bndy_approxs_ = 0;
  AA_num_rows_ = 0;
  AA_num_cols_ = 0;
  AA_ld_ = 0;
  ob_ld_ = 0;
  info_ = 0;
  lwork_ = -1;
  aa_num_rows_ = 0;
  aa_num_cols_ = 0;
  aaT_num_rows_ = 0;
  aaT_num_cols_ = 0;
  ltau_ = 0;
  KK_num_rows_ = 0;
  KK_num_cols_ = 0;
  ldSUBKT_ = 0;
  incx_ = 1;
  PI_num_cols_ = 0;
  lDD_ = 0;
  minrow_ = 0;
  row_ = 0;
}

/*! \brief Constructor #2. */

/*! Constructs an operator starting from given order of accuracy.
    \param [in] order Order of numerical accuracy.
*/
MTK_1DDiv::MTK_1DDiv(int order) {

  side_ = 'L';
  trans_ = 'N';
  transa_ = 'N';
  transb_ = 'N';

  DD_ = nullptr;

  
  pp_= nullptr;
  TT_=nullptr;
  oo_ = nullptr;
  gg_ = nullptr;
  AA_ = nullptr;
  ob_ = nullptr;
  work_ = nullptr;
  aa_ = nullptr;
  tau_ = nullptr;
  QQ_ = nullptr;
  KK_ = nullptr;
  SUBK_ = nullptr;
  SUBKT_ = nullptr;
  II_ = nullptr;
  KKT_ = nullptr;
  NULLS_ = nullptr;
  prem_apps_ = nullptr;
  ob_bottom_ = nullptr;
  PI_ = nullptr;
  qq_ = nullptr;
  qq_lp_ = nullptr;
  PIT_ = nullptr;
  lambdas_ = nullptr;
  alphas_ = nullptr;
  mim_bndy_ = nullptr;

  mimetic_tol_ = mtk::MTK_DefaultMimTol;
  mimetic_threshold_ = mtk::MTK_DefaultMimeticThreshold;

  alpha_ = 1.0;
  beta_  = 0.0;
  norm_ = 0.0;
  minnorm_ = 0.0;
  aux_ = 0.0;
  normerr_ = 0.0;

  ipiv_ = nullptr;

  kk_ = order;
  order_ = kk_;

  nrhs_ = 1;
  TT_ld_ = 0;
  ldoo_ = 0;
  dim_null_ = 0;
  num_bndy_approxs_ = 0;
  AA_num_rows_ = 0;
  AA_num_cols_ = 0;
  AA_ld_ = 0;
  ob_ld_ = 0;
  info_ = 0;
  lwork_ = -1;
  aa_num_rows_ = 0;
  aa_num_cols_ = 0;
  aaT_num_rows_ = 0;
  aaT_num_cols_ = 0;
  ltau_ = 0;
  KK_num_rows_ = 0;
  KK_num_cols_ = 0;
  ldSUBKT_ = 0;
  incx_ = 1;
  PI_num_cols_ = 0;
  lDD_ = 0;
  minrow_ = 0;
  row_ = 0;
}

/*! \brief Constructor #3. */

/*! Constructs an operator starting from given order of accuracy and threshold.
    \param [in] order Order of numerical accuracy.
    \param [in] mimetic_threshold Given mimetic threshold.
*/
MTK_1DDiv::MTK_1DDiv(int order, mtk::MTK_Real mimetic_threshold) {

  side_ = 'L';
  trans_ = 'N';
  transa_ = 'N';
  transb_ = 'N';

  DD_ = nullptr;

  pp_= nullptr;
  TT_=nullptr;
  oo_ = nullptr;
  gg_ = nullptr;
  AA_ = nullptr;
  ob_ = nullptr;
  work_ = nullptr;
  aa_ = nullptr;
  tau_ = nullptr;
  QQ_ = nullptr;
  KK_ = nullptr;
  SUBK_ = nullptr;
  SUBKT_ = nullptr;
  II_ = nullptr;
  KKT_ = nullptr;
  NULLS_ = nullptr;
  prem_apps_ = nullptr;
  ob_bottom_ = nullptr;
  PI_ = nullptr;
  qq_ = nullptr;
  qq_lp_ = nullptr;
  PIT_ = nullptr;
  lambdas_ = nullptr;
  alphas_ = nullptr;
  mim_bndy_ = nullptr;

  mimetic_tol_ = mtk::MTK_DefaultMimTol;
  mimetic_threshold_ = mimetic_threshold;

  alpha_ = 1.0;
  beta_  = 0.0;
  norm_ = 0.0;
  minnorm_ = 0.0;
  aux_ = 0.0;
  normerr_ = 0.0;

  ipiv_ = nullptr;

  kk_ = order;
  order_ = kk_;

  nrhs_ = 1;
  TT_ld_ = 0;
  ldoo_ = 0;
  dim_null_ = 0;
  num_bndy_approxs_ = 0;
  AA_num_rows_ = 0;
  AA_num_cols_ = 0;
  AA_ld_ = 0;
  ob_ld_ = 0;
  info_ = 0;
  lwork_ = -1;
  aa_num_rows_ = 0;
  aa_num_cols_ = 0;
  aaT_num_rows_ = 0;
  aaT_num_cols_ = 0;
  ltau_ = 0;
  KK_num_rows_ = 0;
  KK_num_cols_ = 0;
  ldSUBKT_ = 0;
  incx_ = 1;
  PI_num_cols_ = 0;
  lDD_ = 0;
  minrow_ = 0;
  row_ = 0;
}

/*! \brief Destructor. */

/*! Destroys the built operator. */
MTK_1DDiv::~MTK_1DDiv() {
    
  if (oo_ != nullptr) {
    delete[] oo_;
    oo_ = nullptr;
  }
  if (qq_ != nullptr) {
    delete[] qq_;
    qq_ = nullptr;
  }
  if (qq_lp_ != nullptr) {
    delete[] qq_lp_;
    qq_lp_ = nullptr;
  }
  if (mim_bndy_ != nullptr) {
    delete[] mim_bndy_;
    mim_bndy_ = nullptr;
  }
  if (DD_ != nullptr) {
    delete[] DD_;
    DD_ = nullptr;
  }
  if (TT_ != nullptr) {
    delete[] TT_;
    TT_ = nullptr;
  }
  
  if (lambdas_ != nullptr) {
    delete[] lambdas_;
    lambdas_ = nullptr;
  }

  if (pp_ != nullptr) {
    delete[] pp_;
    pp_ = nullptr;
  }
  if (gg_ != nullptr) {
    delete[] gg_;
    gg_ = nullptr;
  }
  if (AA_ != nullptr) {
    delete[] AA_;
    AA_ = nullptr;
  }

  if (ob_ != nullptr) {
    delete[] ob_;
    ob_ = nullptr;
  }
  if (work_ != nullptr) {
    delete[] work_;
    work_ = nullptr;
  }
  if (tau_ != nullptr) {
    delete[] tau_;
    tau_ = nullptr;
  }
  if (QQ_ != nullptr) {
    delete[] QQ_;
    QQ_ = nullptr;
  }
  if (KK_ != nullptr) {
    delete[] KK_;
    KK_ = nullptr;
  }
  if (SUBK_ != nullptr) {
    delete[] SUBK_;
    SUBK_ = nullptr;
  }
  if (SUBKT_ != nullptr) {
    delete[] SUBKT_;
    SUBKT_ = nullptr;
  }
  if (II_ != nullptr) {
    delete[] II_;
    II_ = nullptr;
  }
  if (KKT_ != nullptr) {
    delete[] KKT_;
    KKT_ = nullptr;
  }
  if (NULLS_ != nullptr) {
    delete[] NULLS_;
    NULLS_ = nullptr;
  }
  if (prem_apps_ != nullptr) {
    delete[] prem_apps_;
    prem_apps_ = nullptr;
  }
  if (ob_bottom_ != nullptr) {
    delete[] ob_bottom_;
    ob_bottom_ = nullptr;
  }
  if (PI_ != nullptr) {
    delete[] PI_;
    PI_ = nullptr;
  }
  if (PIT_ != nullptr) {
    delete[] PIT_;
    PIT_ = nullptr;
  }
  if (lambdas_ != nullptr) {
    delete[] lambdas_;
    lambdas_ = nullptr;
  }
  if (alphas_ != nullptr) {
    delete[] alphas_;
    alphas_ = nullptr;
  }
  if (mim_bndy_ != nullptr) {
    delete[] mim_bndy_;
    mim_bndy_ = nullptr;
  }

  

  
}

/*! \brief Factory method to construct a divergence operator. */

/*! Constructs an operator starting from given order of accuracy and threshold.
    \param [in] kk_ Order of numerical accuracy.
    \param [in] mimetic_threshold Given mimetic threshold.
*/
MTK_1DDiv* MTK_1DDiv::Construct1DDiv() {

#if MTK_DEBUG_LEVEL > 0
  if (kk_ < 2) {
    cerr << "Order of accuracy should be greater or equal than 2." << endl;
    cerr << "Exiting..." << endl;
    return nullptr;
  }
  if ((kk_%2) != 0) {
    cerr << "Order of accuracy should be an even number." << endl;
    cerr << "Exiting..." << endl;
    return nullptr;
  }
  if (kk_ >= 8) {
    cout << "WARNING: Numerical accuracy is too high." << endl;
  }
#endif

  //! Step 1: Compute stencil for the interior cells.
  bool abort_construction = ComputeStencilInteriorGrid();


  if (!abort_construction) {
    cerr << "Could NOT complete stage 1." << endl;
    cerr << "Exiting..." << endl;
    return nullptr;
  }

     
  // At this point, we already have the values for the interior stencil stored
  // in the oo array.

  // It is noteworthy, that the 2nd order accurate divergence operator has NO
  // approximation at the boundary, thus it has no weights. For this case, the
  // dimension of the null-space of the Vandermonde matrices used to compute the
  // approximating coefficients at the boundary is 0. Ergo, we compute this
  // number first and then decide if we must compute anything at the boundary:

  dim_null_ = kk_/2 - 1;

  if (dim_null_ > 0) {

      num_bndy_approxs_ = (int) (3.0*((MTK_Real) kk_)/2.0);

      //! Step 2: Compute scaled null space for boundary approximation.

      // Compute the null-space from the first matrix transposed. For this we
      // will follow recommendations given in:
      //
      // http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=5&t=4506
      //
      // We will compute the QR Factorization of the transpose, as in the
      // following (MATLAB) pseudo-code:
      //
      // [Q,R] = qr(V'); % Full QR as defined in
      // % http://www.stanford.edu/class/ee263/notes/qr_matlab.pdf
      // null-space = Q(:, last (kk_/2 - 1) columns of Q );
      //
      // However, given the nature of the Vandermonde matrices we've just
      // computed, they all posses the same null-space. Therefore, we impose the
      // convention of computing the null-space of the first Vandermonde matrix
      // (west boundary).

      abort_construction = ComputeScaledNullSpace();

#if MTK_DEBUG_LEVEL > 0
      if (!abort_construction) {
        cerr << "Could NOT complete stage 2.1." << endl;
        cerr << "Exiting..." << endl;
        return nullptr;
      }
#endif

      //! Step 3: Compute boundary conditions scaled approximation.
      abort_construction = ComputePreliminaryApproximations();

#if MTK_DEBUG_LEVEL > 0
      if (!abort_construction) {
        cerr << "Could NOT complete stage 2.2." << endl;
        cerr << "Exiting..." << endl;
        return nullptr;
      }
#endif

      //! Step 4: Compute weights.

      // Assemble the PI_ matrix using:
      // 1. The collection of scaled preliminary approximations.
      // 2. The collection of coefficients approximating at the interior.
      // 3. The scaled basis for the null-space.
      abort_construction = ComputeWeights();

#if MTK_DEBUG_LEVEL > 0
      if (!abort_construction) {
        cerr << "Could NOT complete stage 2.3." << endl;
        cerr << "Exiting..." << endl;
        return nullptr;
      }
#endif

      //! Step 5: Compute stencil boundary grid.
      abort_construction = ComputeStencilBoundaryGrid();

      // Compute mimetic boundary approximations:
#if MTK_DEBUG_LEVEL > 0
      if (!abort_construction) {
        cerr << "Could NOT complete stage 2.4." << endl;
        cerr << "Exiting..." << endl;
        return nullptr;
      }
#endif

      if (NULLS_ != nullptr) {
        delete[] NULLS_;
        NULLS_ = nullptr;
      }
      if (prem_apps_ != nullptr) {
        delete[] prem_apps_;
        prem_apps_ = nullptr;
      }
    } // End of: if (dim_null_ > 0);

    // Once we have the following three collections of data:
    //   (a) the coefficients for the interior,
    //   (b) the coefficients for the boundary (if it applies),
    //   (c) and the weights (if it applies),
    // we will store everything in the output array:

    //! Step 6: Assemble operator.

    abort_construction = AssembleOperator();

#if MTK_DEBUG_LEVEL > 0
    if (!abort_construction) {
      cerr << "Could NOT complete stage 3." << endl;
      cerr << "Exiting..." << endl;
      return nullptr;
    }
#endif

    return this;
}

/*! \brief Stage 1 of the CRS Algorithm. */

/*! Compute the stencil approximating the interior of the staggered grid.
*/
bool MTK_1DDiv::ComputeStencilInteriorGrid() {

  // Interior spatial coordinates vector:
  try {
    pp_ = new MTK_Real[kk_];
    memset (pp_, 0.0, sizeof(pp_[0])*kk_);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  pp_[0] = 1.0/2.0 - ((MTK_Real) kk_)/2.0;
  for (auto ii = 1; ii < kk_; ii++) {
    pp_[ii] = pp_[ii - 1] + 1.0;
  }

#if MTK_DEBUG_LEVEL >0
  cout << "pp_ = " << endl;
  for (auto ii = 0; ii < kk_; ii++) {
    cout << pp_[ii] << " ";
  }
  cout << endl << endl;
#endif

  // Vandermonde matrix (using interior coordinates as generator):
  TT_= Vandermonde  (pp_, kk_, kk_, false);
  if (TT_ == nullptr) {
    cerr << "ERROR constructing Vandermonde matrix at line " <<
      __LINE__ - 2 << endl;
    cerr << "Exiting..." << endl;
    return false;
  }
#if MTK_DEBUG_LEVEL>0
  cout << "kk_ = "<< kk_ <<endl;
  cout << "Mimetic Threshold = " << mimetic_threshold_ << endl;
  cout << "TT = " << endl;
  if (!DenseMatrix_Print(TT_, kk_, kk_)) {
    return false;
  }
#endif

  // Order-selector vector:
  try {
    oo_ = new MTK_Real[kk_];
    memset (oo_, 0.0, sizeof(oo_[0])*kk_);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }
  oo_[1] = 1.0;

#if MTK_DEBUG_LEVEL>0
  cout << "oo =" << endl;
  for (auto ii = 0; ii < kk_; ii++) {
    cout << setw(12) << oo_[ii];
  }
  cout << endl << endl;
#endif

  // Solve systems of equations to obtain the stencil:
  try {
    ipiv_ = new int[kk_];
    memset(ipiv_, 0, sizeof(ipiv_[0])*kk_);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }
  ldoo_=kk_;
  TT_ld_=kk_;
  dgesv_(&kk_, &nrhs_, TT_, &TT_ld_, ipiv_, oo_, &ldoo_, &info_);
#if MTK_DEBUG_LEVEL >0
  if (!info_) {
    cout << "System successfully solved! Interior stencil attained!" << endl;
  }
  else {
    cerr << "Something went wrong solving system! info = " << info_ << endl;
    cerr << "Exiting..." << endl;
    return false;
  }
#endif

#if MTK_DEBUG_LEVEL>0
  cout << endl;
  cout << "ss =" << endl;
  for (auto ii = 0; ii < kk_; ii++) {
      cout << setw(12) << oo_[ii];
  }
  cout << endl << endl;
#endif
  free (ipiv_);
  free (pp_);
  free(TT_);
  return true;
}

// Stage 2.1 of the CRS Algorithm:
// Compute a scaled basis of the null-space of the Vandermonde matrix
// approximating at the west boundary:
bool MTK_1DDiv::ComputeScaledNullSpace(void) {

    // Generator vector for the first approximation:
    try {
        gg_ = new MTK_Real[num_bndy_approxs_];
        memset(gg_, 0.0, sizeof(gg_[0])*num_bndy_approxs_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    gg_[0] = -1.0/2.0;
    for (auto ii = 1; ii < num_bndy_approxs_; ii++) {
        gg_[ii] = gg_[ii - 1] + 1.0;
    }
#if MTK_DEBUG_LEVEL >0
    cout << "gg_west =" << endl;
    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        cout << setw(12) << gg_[ii];
    }
    cout << endl << endl;
#endif


// Construct the TRANSPOSE of the first Vandermonde matrix:
    aa_num_rows_ = kk_ + 1;
    aa_num_cols_ = num_bndy_approxs_;

    try {
        aa_ = new MTK_Real[aa_num_rows_*aa_num_cols_];
        memset(aa_, 0.0, sizeof(aa_[0])*(aa_num_rows_*aa_num_cols_));
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    for (auto ii = 0; ii < aa_num_rows_; ii++) {
        for (auto jj = 0; jj < aa_num_cols_; jj++) {
            aa_[ii*aa_num_cols_ + jj] = pow(gg_[jj], ii);
        }
    }

#if MTK_DEBUG_LEVEL>0
    cout << "aa_west = " << endl;
    for (auto ii = 0; ii < aa_num_rows_; ii++) {
        for (auto jj = 0; jj < aa_num_cols_; jj++) {
            cout << setw(12) << aa_[ii*aa_num_cols_ + jj];
        }
        cout << endl;
    }
    cout << endl;
#endif

    aaT_num_rows_ = aa_num_cols_;
    aaT_num_cols_ = aa_num_rows_;

// Factorize aa:

// Prepare to factorize: allocate and inquire for the value of lwork:
    try {
        work_ = new MTK_Real[1];
        memset(work_, 0.0, sizeof(aa_[0])*1);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    lwork_ = -1;
    info_ = 0;

    dgeqrf_(&aaT_num_rows_, &aaT_num_cols_, aa_, &aaT_num_rows_, tau_,
            work_, &lwork_, &info_);

    if (info_ == 0) {
        lwork_ = (int) work_[0];
    } else {
        cerr << "Could not get value for lwork on line " << __LINE__ - 5 << endl;
        cerr << "Exiting..." << endl;
        return EXIT_FAILURE;
    }

#if MTK_DEBUG_LEVEL>0
    cout << "lwork_ = " << endl << setw(12)<< lwork_ << endl << endl;
#endif

    if (work_ != nullptr) {
        delete[] work_;
        work_ = nullptr;
    }

// Once we know lwork, we can actually invoke the factorization:
    try {
        work_ = new MTK_Real [lwork_];
        memset(work_, 0.0, sizeof(work_[0])*lwork_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    ltau_ = min(aaT_num_rows_,aaT_num_cols_);

    try {
        tau_ = new MTK_Real [ltau_];
        memset(tau_, 0.0, sizeof(tau_[0])*ltau_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    dgeqrf_(&aaT_num_rows_, &aaT_num_cols_, aa_, &aaT_num_rows_,
            tau_, work_, &lwork_, &info_);

    if (!info_) {
#if MTK_DEBUG_LEVEL>0
        cout << "QR factorization successfully completed!" << endl << endl;
#endif
    } else {
        cerr << "Something went wrong solving system! info = " << info_ << endl;
        cerr << "Exiting..." << endl;
        return EXIT_FAILURE;
    }

#if MTK_DEBUG_LEVEL>0
    cout << "aa_west AFTER factorization:" << endl;
    for (auto ii = 0; ii < aa_num_rows_; ii++) {
        for (auto jj = 0; jj < aa_num_cols_; jj++) {
            cout << setw(13) << aa_[ii*aa_num_cols_ + jj];
        }
        cout << endl;
    }
    cout << endl;
#endif

// We now generate the real matrix Q with orthonormal columns. This has to
// be done separately since the actual output of dgeqrf_ (AA_) represents
// the orthogonal matrix Q as a product of min(aa_num_rows,aa_num_cols)
// elementary Householder reflectors. Notice that we must re-inquire the new
// value for lwork that is used.
    try {
        QQ_ = new MTK_Real[num_bndy_approxs_*num_bndy_approxs_];
        memset(QQ_, 0.0, sizeof(QQ_[0])*num_bndy_approxs_*num_bndy_approxs_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

// QQ_ must be initialized to the identity matrix:
    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        for (auto jj = 0; jj < num_bndy_approxs_; jj++) {
            QQ_[ii*num_bndy_approxs_ + jj] = (ii == jj)? 1.0: 0.0;
        }
    }

#if MTK_DEBUG_LEVEL >0
    cout << "Initialized QQ_T:" << endl;
    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        for (auto jj = 0; jj < num_bndy_approxs_; jj++) {
            cout << setw(13) << QQ_[ii*num_bndy_approxs_ + jj];
        }
        cout << endl;
    }
    cout << endl;
#endif

// Assemble the QQ_ matrix:
    lwork_ = -1;

    if (work_ != nullptr) {
        delete[] work_;
        work_ = nullptr;
    }

    try {
        work_ = new MTK_Real[1];
        memset(work_, 0.0, sizeof(work_[0])*1);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    dormqr_(&side_, &trans_, &aa_num_cols_, &aa_num_cols_, &ltau_, aa_,
            &aaT_num_rows_,
            tau_, QQ_, &num_bndy_approxs_, work_, &lwork_, &info_);

    if (info_ == 0) {
        lwork_ = (int) work_[0];
    } else {
        cerr << "Could not get value for lwork_ on line " << __LINE__ - 5 << endl;
        cerr << "Exiting..." << endl;
        return EXIT_FAILURE;
    }

#if MTK_DEBUG_LEVEL>0
    cout << "lwork = " << endl << setw(12)<< lwork_ << endl << endl;
#endif

    if (work_ != nullptr) {
        delete[] work_;
        work_ = nullptr;
    }

    try {
        work_ = new MTK_Real[lwork_];
        memset(work_, 0.0, sizeof(work_[0])*lwork_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    dormqr_(&side_, &trans_, &aa_num_cols_, &aa_num_cols_, &ltau_, aa_,
            &aaT_num_rows_,
            tau_, QQ_, &num_bndy_approxs_, work_, &lwork_, &info_);

    if (!info_) {
#if MTK_DEBUG_LEVEL>0
        cout << "Q matrix successfully assembled!" << endl << endl;
#endif
    } else {
        cerr << "Something went wrong solving system! info = " << info_ << endl;
        cerr << "Exiting..." << endl;
        return EXIT_FAILURE;
    }

#if MTK_DEBUG_LEVEL >0
    cout << "QQ_T =" << endl;
    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        for (auto jj = 0; jj < num_bndy_approxs_; jj++) {
            cout << setw(13) << QQ_[ii*num_bndy_approxs_ + jj];
        }
        cout << endl;
    }
    cout << endl;
#endif

// Extract the null-space from QQ, and save it in KK_:

    KK_num_rows_ = num_bndy_approxs_;
    KK_num_cols_ = dim_null_;

    try {
        KK_ = new MTK_Real[KK_num_rows_*KK_num_cols_];
        memset(KK_, 0.0, sizeof(KK_[0])*(KK_num_rows_*KK_num_cols_));
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    for (auto ii = (num_bndy_approxs_ - dim_null_); ii < num_bndy_approxs_; ii++) {
        for (auto jj = 0; jj < num_bndy_approxs_; jj++) {
            KK_[idx(jj,dim_null_,(ii - (num_bndy_approxs_ - dim_null_)))] =
                QQ_[idx(ii,num_bndy_approxs_,jj)];
        }
    }

#if MTK_DEBUG_LEVEL >0
    cout << "KK =" << endl;
    for (auto ii = 0; ii < KK_num_rows_; ii++) {
        for (auto jj = 0; jj < KK_num_cols_; jj++) {
            cout << setw(15) << KK_[idx(ii,KK_num_cols_,jj)];
        }
        cout << endl;
    }
    cout << endl;
#endif

// Scale KK:

// Scale thus requesting that the last entries of the attained basis for the
// null-space, adopt the pattern we require.
// Essentially we will implement the following MATLAB pseudo-code:
//  scalers = KK(num_bndy_approxs - (dim_null - 1):num_bndy_approxs,:)\B
//  SK = KK*scalers
// where SK is the scaled null-space.

// In this point, we almost have all the data we need correctly allocated
// in memory. We will create the matrix II_, and elements we wish to scale in
// the KK array. Using the concept of the leading dimension, we could just
// use KK, with the correct leading dimension and that is it. BUT I DO NOT
// GET how does it work. So I will just create a matrix with the content of
// this array that we need, solve for the scalers and then scale the
// whole KK:

// We will then create memory for that sub-matrix of KK (SUBK_):

    try {
        SUBK_ = new MTK_Real[dim_null_*dim_null_];
        memset(SUBK_, 0.0, sizeof(SUBK_[0])*dim_null_*dim_null_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

// Extracting information from KK:
    for (auto ii = num_bndy_approxs_ - dim_null_; ii < num_bndy_approxs_; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            SUBK_[idx(ii - (num_bndy_approxs_ - dim_null_),dim_null_,jj)] =
                KK_[idx(ii,dim_null_,jj)];
        }
    }

#if MTK_DEBUG_LEVEL >0
    cout << "Required SUBK_ for scaling:" << endl;
    for (auto ii = 0; ii < dim_null_; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            cout << setw(15) << SUBK_[idx(ii,dim_null_,jj)];
        }
        cout << endl;
    }
    cout << endl;
#endif

// Once we have the matrix for the scaling, create the collection of rhss,
// which looks like an identity matrix.
    try {
        II_ = new MTK_Real[dim_null_*dim_null_];
        memset(II_, 0.0, sizeof(II_[0])*dim_null_*dim_null_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    for (auto ii = 0; ii < dim_null_; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            II_[ii*dim_null_ + jj] = (ii == jj)? 1.0: 0.0;
        }
    }

#if MTK_DEBUG_LEVEL >0
    cout << "Required II_ for scaling:" << endl;
    for (auto ii = 0; ii < dim_null_; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            cout << setw(15) << II_[ii*dim_null_ + jj];
        }
        cout << endl;
    }
    cout << endl;
#endif

// Solve the system to compute the scalers.
// An example of the system to solve, for k = 8, is:
//
// SUBK*scalers = II_ or
//
// |  0.386018  -0.0339244   -0.129478 |           | 1 0 0 |
// | -0.119774   0.0199423   0.0558632 |*scalers = | 0 1 0 |
// | 0.0155708 -0.00349546 -0.00853182 |           | 0 0 1 |
//
// Notice this is a nrhs = 3 system.
// Noteworthy: we do NOT ACTUALLY ALLOCATE space for the scalers... they
// will be stored in the created identity matrix.
// Let us first transpose SUBK (because of LAPACK):

    SUBKT_ = Transpose(SUBK_, dim_null_, dim_null_);
    if (SUBKT_ == nullptr) {
        cerr << "ERROR constructing  matrix at line " << __LINE__ - 2 << endl;
        cerr << "Exiting..." << endl;
        return false;
    }

#if MTK_DEBUG_LEVEL >0
    cout << "Required SUBK transposed for scaling:" << endl;
    for (auto ii = 0; ii < dim_null_; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            cout << setw(15) << SUBK_[idx(ii,dim_null_,jj)];
        }
        cout << endl;
    }
    cout << endl;
#endif

// All right! Ready to solve using the LAPACK:
    nrhs_ = dim_null_;
    ldSUBKT_ = dim_null_;
    info_ = 0;

    try {
        ipiv_ = new int[ldSUBKT_];
        memset(ipiv_, 0, sizeof(ipiv_[0])*ldSUBKT_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    dgesv_(&ldSUBKT_, &nrhs_,SUBKT_, &ldSUBKT_, ipiv_, II_, &ldSUBKT_,
           &info_);
    if (!info_) {
#if MTK_DEBUG_LEVEL >0
        cout << "System successfully solved! Scalers attained!" << endl << endl;
#endif
    } else {
        cerr << "Something went wrong solving system! info = " << info_ << endl;
        return EXIT_FAILURE;
    }

#if MTK_DEBUG_LEVEL >0
    cout << "Computed scalers:" << endl;
    for (auto ii = 0; ii < dim_null_; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            cout << setw(15) << II_[idx(ii,dim_null_,jj)];
        }
        cout << endl;
    }
    cout << endl;
#endif

// We now multiply the two matrices to attain a scaled basis for null-space:

// Transpose the basis matrix KK:
    KKT_ = Transpose(KK_, num_bndy_approxs_, dim_null_);
    if (KKT_ == nullptr) {
        cerr << "ERROR constructing matrix at " << __LINE__ - 2 << endl;
        cerr << "Exiting..." << endl;
        return false;
    }

// Allocate space for the product matrix:
    try {
        NULLS_ = new MTK_Real[num_bndy_approxs_*dim_null_];
        memset(NULLS_, 0.0, sizeof(NULLS_[0])*num_bndy_approxs_*dim_null_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

// Perform matrix-matrix multiplication:
    dgemm_(&transa_, &transb_, &num_bndy_approxs_, &dim_null_, &dim_null_,
           &alpha_, KKT_, &num_bndy_approxs_, II_, &dim_null_, &beta_,
           NULLS_, &num_bndy_approxs_);
#if MTK_DEBUG_LEVEL >0
    cout << "Scaled basis for the null-space (tricky print):" << endl;
    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            auto value = NULLS_[idx(jj,num_bndy_approxs_,ii)];
            if (abs(value - 0.0) < mimetic_tol_) {
                cout << setw(12) << 0.0;
            } else {
                cout << setw(12) << value;
            }
        }
        cout << endl;
    }
    cout << endl;
#endif

// At this point, we have a scaled basis for the null-space, with the
// pattern we need!

// Local garbage collection :)
    if (gg_ != nullptr) {
        delete[] gg_;
        gg_ = nullptr;
    }
    if (aa_ != nullptr) {
        delete[] aa_;
        aa_ = nullptr;
    }
    if (work_ != nullptr) {
        delete[] work_;
        work_ = nullptr;
    }
    if (tau_ != nullptr) {
        delete[] tau_;
        tau_ = nullptr;
    }
    if (QQ_ != nullptr) {
        delete[] QQ_;
        QQ_ = nullptr;
    }
    if (KK_ != nullptr) {
        delete[] KK_;
        KK_ = nullptr;
    }
    if (SUBK_ != nullptr) {
        delete[] SUBK_;
        SUBK_ = nullptr;
    }
    if (II_ != nullptr) {
        delete[] II_;
        II_ = nullptr;
    }
    if (SUBKT_ != nullptr) {
        delete[] SUBKT_;
        SUBKT_ = nullptr;
    }
    if (ipiv_ != nullptr) {
        delete[] ipiv_;
        ipiv_ = nullptr;
    }
    if (KKT_ != nullptr) {
        delete[] KKT_;
        KKT_ = nullptr;
    }

    return true;
}

// Stage 2.2 of the CRS Algorithm:
// Compute the set of preliminary approximation on the boundary neighborhood:
bool MTK_1DDiv::ComputePreliminaryApproximations(void) {

    // Generator vector for the first approximation:
    try {
        gg_ = new MTK_Real[num_bndy_approxs_];
        memset(gg_, 0.0, sizeof(gg_[0])*num_bndy_approxs_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    gg_[0] = -1.0/2.0;
    for (auto ii = 1; ii < num_bndy_approxs_; ii++) {
        gg_[ii] = gg_[ii - 1] + 1.0;
    }

#if MTK_DEBUG_LEVEL>0
    cout << "gg_0 =" << endl;
    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        cout << setw(12) << gg_[ii];
    }
    cout << endl << endl;
#endif

    // Allocate 2D array to store the collection of preliminary approximations:
    try {
        prem_apps_ = new MTK_Real[num_bndy_approxs_*dim_null_];
        memset(prem_apps_, 0.0, sizeof(prem_apps_[0])*num_bndy_approxs_*dim_null_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    // Compute the dim_null near-the-boundary columns of the PI_ matrix:
    for (auto ll = 0; ll < dim_null_; ll++) {

        // Re-check new generator vector for every iteration except for the first:
        if (ll > 0) {
#if MTK_DEBUG_LEVEL>0
            cout << "gg_" << ll << " =" << endl;
            for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
                cout << setw(12) << gg_[ii];
            }
            cout << endl << endl;
#endif
        }

        // Create the Vandermonde matrix for this iteration:
        AA_ = Vandermonde(gg_, num_bndy_approxs_, kk_ + 1, false);
        if (AA_ == nullptr) {
            cerr << "ERROR constructing matrix at line " << __LINE__ - 2 << endl;
            cerr << "Exiting..." << endl;
            return EXIT_FAILURE;
        }

        AA_num_rows_ = kk_ + 1;
        AA_num_cols_ = num_bndy_approxs_;
        AA_ld_ = AA_num_rows_;

#if MTK_DEBUG_LEVEL >0
        // Although it is generated transposed, we print it non-transposed:
        cout << "AA^T_" << ll << " =" << endl;
        if (!DenseMatrix_Print(AA_, AA_num_rows_, AA_num_cols_)) {
            return EXIT_FAILURE;
        }
#endif

        // New order-selector vector (it gets re-written with LAPACK solutions):
        ob_ld_ = num_bndy_approxs_;

        try {
            ob_ = new MTK_Real[ob_ld_];
            memset(ob_, 0.0, sizeof(ob_[0])*ob_ld_);
        } catch (bad_alloc &memory_allocation_exception) {
            cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
            cerr << memory_allocation_exception.what() << endl;
        }

        ob_[1] = 1.0;

#if MTK_DEBUG_LEVEL >0
        cout << "ob = " << endl << endl;
        for (auto ii = 0; ii < ob_ld_; ii++) {
            cout << setw(12) << ob_[ii] << endl;
        }
        cout << endl;
#endif

        // Solving TT*rr = ob_ yields the columns rr of the KK matrix. However,
        // this is an under-determined system of equations. So we can not use the
        // same LAPACK routine (dgesv_). We will instead use dgels_.
        // We first invoke the solver to query for the value of lwork. For this,
        // we must at least allocate enough space to allow access to WORK(1), or
        // work[0]:
        // If LWORK = -1, then a workspace query is assumed; the routine only
        // calculates the optimal size of the WORK array, returns this value as
        // the first entry of the WORK array, and no error message related to
        // LWORK is issued by XERBLA.
        if (work_ != nullptr) {
            delete[] work_;
            work_ = nullptr;
        }
        try {
            work_ = new MTK_Real[1];
            memset(work_, 0.0, sizeof(work_[0])*1);
        } catch (bad_alloc &memory_allocation_exception) {
            cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
            cerr << memory_allocation_exception.what() << endl;
        }

        // IMPORTANT: Never forget to re-set the nrhs variable to the right
        // number:
        nrhs_ = 1;
        info_ = 0;
        lwork_ = -1;

        dgels_(&trans_, &AA_num_rows_, &AA_num_cols_, &nrhs_, AA_, &AA_ld_, ob_,
               &ob_ld_,
               work_, &lwork_, &info_);

        if (info_ == 0) {
            lwork_ = (int) work_[0];
        } else {
            cerr << "Could not get value for lwork on line " << __LINE__ - 2 << endl;
            cerr << "Exiting..." << endl;
            return EXIT_FAILURE;
        }

#if MTK_DEBUG_LEVEL >0
        cout << "lwork = " << endl << setw(12)<< lwork_ << endl << endl;
#endif

        // We then use lwork's new value to create the work array:
        if (work_ != nullptr) {
            delete[] work_;
            work_ = nullptr;
        }

        try {
            work_ = new MTK_Real[lwork_];
            memset(work_, 0.0, sizeof(work_[0])*lwork_);
        } catch (bad_alloc &memory_allocation_exception) {
            cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
            cerr << memory_allocation_exception.what() << endl;
        }

        // We now invoke the solver again:
        dgels_(&trans_, &AA_num_rows_, &AA_num_cols_, &nrhs_, AA_, &AA_ld_, ob_,
               &ob_ld_,
               work_, &lwork_, &info_);

        if (!info_) {
#if MTK_DEBUG_LEVEL >0
            cout << "System successfully solved!" << endl << endl;
#endif
        } else {
            cerr << "Something went wrong solving system! info = " << info_ << endl;
        }
#if MTK_DEBUG_LEVEL >0
        cout << "ob =" << endl;
        for (auto ii = 0; ii < ob_ld_; ii++) {
            cout << setw(12) << ob_[ii] << endl;
        }
        cout << endl;
#endif

        // Prior to save this solution, we scale it, using the scaled basis for
        // the null-space. This implies a DAXPY operation. However, we must
        // construct the arguments for this operation:
        // First, we extract the last dim_null values of the pre-scaled ob, and
        // save them into the ob_bottom_ array:

        try {
            ob_bottom_ = new MTK_Real[dim_null_];
            memset(ob_bottom_, 0.0, sizeof(ob_bottom_[0])*dim_null_);
        } catch (bad_alloc &memory_allocation_exception) {
            cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
            cerr << memory_allocation_exception.what() << endl;
        }

        for (auto ii = 0; ii < dim_null_; ii++) {
            ob_bottom_[(dim_null_ - 1) - ii] = ob_[num_bndy_approxs_ - ii - 1];
        }

#if MTK_DEBUG_LEVEL >0
        cout << "ob_bottom_ =" << endl;
        for (auto ii = 0; ii < dim_null_; ii++) {
            cout << setw(12) << ob_bottom_[ii] << endl;
        }
        cout << endl << endl;
#endif

        // Once we posses the bottom elements, we proceed with the scaling. We
        // must computed an scaled ob, sob, using the scaled null-space NULLS_.
        // Such operation is: sob = ob - NULLS_*ob_bottom_
        // or:                 ob = -1.0*NULLS_*ob_bottom_ + 1.0*ob
        // thus:                Y =    a*A    *x         +   b*Y (DAXPY).
        // Let us invoke this from the BLAS:
        trans_ = 'N';
        alpha_ = -1.0;
        beta_ = 1.0;
        dgemv_(&trans_, &num_bndy_approxs_, &dim_null_,
               &alpha_, NULLS_, &num_bndy_approxs_,
               ob_bottom_, &incx_, &beta_, ob_, &incx_);
#if MTK_DEBUG_LEVEL >0
        cout << "scaled ob_ (tricky print):" << endl;
        for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
            auto value = ob_[ii];
            if (abs(value - 0.0) < mimetic_tol_) {
                cout << setw(12) << 0.0 << endl;
            } else {
                cout << setw(12) << value << endl;
            }
        }
        cout << endl;
#endif

        // We save the recently scaled solution, into an array containing these.
        // We can NOT start building the PI_ matrix, simply because I want that part
        // to be separated since its construction depends on the algorithm we want
        // to implement.

        for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
            prem_apps_[ii*dim_null_ + ll] = ob_[ii];
        }

        // After the first iteration, simply shift the entries of the last
        // generator vector used:
        for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
            gg_[ii]--;
        }

        // Garbage collection for this loop:
        if (AA_ != nullptr) {
            delete[] AA_;
            AA_ = nullptr;
        }
        if (ob_ != nullptr) {
            delete[] ob_;
            ob_ = nullptr;
        }
        if (work_ != nullptr) {
            delete[] work_;
            work_ = nullptr;
        }
        if (ob_bottom_ != nullptr) {
            delete[] ob_bottom_;
            ob_bottom_ = nullptr;
        }
    } // End of: for (ll = 0; ll < dim_null; ll++);

#if MTK_DEBUG_LEVEL >0
    cout << "Matrix post-scaled preliminary apps (tricky print): " << endl;
    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            auto value = prem_apps_[idx(ii,dim_null_,jj)];
            cout << setw(12) << ((abs(value - 0.0) < mimetic_tol_)? 0.0: value);
        }
        cout << endl;
    }
    cout << endl;
#endif

    if (gg_ != nullptr) {
        delete[] gg_;
        gg_ = nullptr;
    }

    return true;
}

// Stage 2.3 of the CRS Algorithm:
// Assemble the PI_ matrix and compute the weights:
bool MTK_1DDiv::ComputeWeights(void) {


    // Collect the preliminary approximations on the PI_ matrix:
    try {
        PI_ = new MTK_Real[num_bndy_approxs_*(num_bndy_approxs_ - 1)];
        memset(PI_, 0.0, sizeof(PI_[0])*num_bndy_approxs_*(num_bndy_approxs_ - 1));
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }


    // We first process the queue of scaled near-the-boundary preliminary
    // solutions, which are queued in scaled_solutions. Each one of these,
    // will be a column of the PI_ matrix:
    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            PI_[ii*(2*dim_null_ + (kk_/2 + 1)) + jj] = prem_apps_[ii*dim_null_ + jj];
        }
    }
    // We now add the columns from the known stencil approximating at the
    // interior, however, these must be padded by zeros, according to their
    // position in the final PI_ matrix:
    auto mm = 0;
    for (auto jj = dim_null_; jj < kk_; jj++) {
        for (auto ii = 0; ii < kk_; ii++) {
            PI_[(ii + mm)*(2*dim_null_ + (kk_/2 + 1)) + jj] = oo_[ii];
        }
        mm++;
    }

    // We finally add the final set of columns: the scaled basis for the
    // null-space:
    for (auto jj = dim_null_ + (kk_/2 + 1); jj < num_bndy_approxs_ - 1; jj++) {
        for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
            PI_[ii*(2*dim_null_ + (kk_/2 + 1)) + jj] =
                NULLS_[(jj - (dim_null_ + (kk_/2 + 1)))*num_bndy_approxs_ + ii];
        }
    }
 //   cout.precision(4);

#if MTK_DEBUG_LEVEL >0
    cout << endl;
    cout << "ss =" << endl;
    for (auto ii = 0; ii < kk_; ii++) {
        cout << setw(11) << oo_[ii];
    }
    cout << endl << endl;
#endif

    PI_num_cols_ = num_bndy_approxs_ - 1;

#if MTK_DEBUG_LEVEL >0
    cout << "Constructed PI_ matrix with CRS Algorithm: " << endl;
    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        for (auto jj = 0; jj < PI_num_cols_; jj++) {
            if (abs(PI_[ii*PI_num_cols_ + jj] - 0.0) < mimetic_tol_) {
                cout << setw(11) << 0.0;
            } else {
                cout << setw(11) << PI_[ii*PI_num_cols_ + jj];
            }
        }
        cout << endl;
    }
    cout << endl;
#endif

    // Once we have constructed the matrix, we use the interior stencil to
    // build the proper RHS vector h that imposes the mimetic condition.
    // That is, we construct a system PI_*q = h, to solve for the weights:
    try {
        qq_ = new MTK_Real[num_bndy_approxs_];
        memset(qq_, 0.0, sizeof(qq_[0])*num_bndy_approxs_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }
    ///********** RAUL: Adding allocation for qq_lp_ (optimizer)
    try {
        qq_lp_ = new MTK_Real[num_bndy_approxs_];
        memset(qq_lp_, 0.0, sizeof(qq_lp_[0])*num_bndy_approxs_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }
    //*****************************

    qq_[0] = -1.0;
    for (auto ii = (kk_/2 + 2 - 1); ii < num_bndy_approxs_; ii++) {
        auto aux_xx = 0.0;
        for (auto jj = 0; jj < ((ii - (kk_/2 - 1)) - 1); jj++) {
            aux_xx += oo_[jj];
        }
        qq_[ii] = -1.0*aux_xx;
    }
    // RAUL: Copy qq_ to qq_lp_ for the optimizer
    for(auto ii =0; ii <num_bndy_approxs_; ii++)
    {
        qq_lp_[ii]=qq_[ii];
    }
    //**************************

#if MTK_DEBUG_LEVEL >0
    cout << "hh =" << endl;
    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        cout << setw(11) << qq_[ii] << endl;
    }
    cout << endl;
#endif

    // Since we intend to solve, we must transpose PI_.
    // TODO: We could create PI_, already transposed.
    PIT_ = Transpose(PI_, num_bndy_approxs_, PI_num_cols_);
    if (PIT_ == (MTK_Real*) nullptr) {
        cerr << "Could not create matrix at line " << __LINE__ - 2 << endl;
        cerr << "Exiting..." << endl;
        return EXIT_FAILURE;
    }

#if MTK_DEBUG_LEVEL >0
    cout << "Required PIT_ Transposed for solving for weights:" << endl;
    for (auto ii = 0; ii < PI_num_cols_; ii++) {
        for (auto jj = 0; jj < num_bndy_approxs_; jj++) {
            if (abs(PIT_[ii*num_bndy_approxs_ + jj] - 0.0) < mimetic_tol_) {
                cout << setw(12) << 0.0;
            } else {
                cout << setw(12) << PIT_[ii*num_bndy_approxs_ + jj];
            }
        }
        cout << endl;
    }
    cout << endl;
#endif

    // At this point, we are ready to solve for the weights. Once again we face
    // the challenge of solving with LAPACK. However, for the CRSA, this matrix
    // PI is over-determined, since it has more rows than unknowns. However,
    // according to the theory, the solution to this system is unique. We will use
    // dgels_.
    // Recall that we must first invoke the solver to query for the value of
    // lwork. For this, we must at least allocate enough space to allow access to
    // WORK(1), or work[0], because:

    // If LWORK = -1, then a workspace query is assumed; the routine only
    // calculates the optimal size of the WORK array, returns this value as the
    // first entry of the WORK array, and no error message related to LWORK is
    // issued by XERBLA.

    try {
        work_ = new MTK_Real[1];
        memset(work_, 0.0, sizeof(aa_[0])*1);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    nrhs_ = 1;
    lwork_ = -1;

    dgels_(&trans_, &num_bndy_approxs_, &PI_num_cols_, &nrhs_, PIT_,
           &num_bndy_approxs_,
           qq_, &num_bndy_approxs_, work_, &lwork_, &info_);

    if (info_ == 0) {
        lwork_ = (int) work_[0];
    } else {
        cerr << "Could not get value for lwork on line " << __LINE__ - 5 << endl;
        cerr << "Exiting..." << endl;
        return EXIT_FAILURE;
    }

#if MTK_DEBUG_LEVEL >0
    cout << "lwork = " << endl << setw(12)<< lwork_ << endl << endl;
#endif

    if (work_ != nullptr) {
        delete[] work_;
        work_ = nullptr;
    }

    try {
        work_ = new MTK_Real[lwork_];
        memset(work_, 0.0, sizeof(work_[0])*lwork_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    // We now invoke the solver again:
    trans_ = 'N';
    nrhs_ = 1;
    info_ = 0;
    dgels_(&trans_, &num_bndy_approxs_, &PI_num_cols_, &nrhs_, PIT_,
           &num_bndy_approxs_,
           qq_, &num_bndy_approxs_, work_, &lwork_, &info_);

    if (!info_) {
#if MTK_DEBUG_LEVEL >0
        cout << "System successfully solved!" << endl << endl;
#endif
    } else {
        cerr << "Something went wrong solving system! info = " << info_ <<
             endl;
    }
#if MTK_DEBUG_LEVEL >0
    cout << "qq_CRS =" << endl;
    for (auto ii = 0; ii < PI_num_cols_; ii++) {
        cout << setw(11) << qq_[ii] << endl;
    }
    cout << endl;
#endif

    // Then qq_lp contains hh and the last dim_null-1 elements are the lambdas_

       lambdas_=(mtk::MTK_Real *)malloc(dim_null_*sizeof(mtk::MTK_Real));
       for (auto ii=0; ii<dim_null_; ii++)
       {
        
          lambdas_[ii]=qq_[ii+kk_];

           
      }
    //  This part should be executed only if kk_ > 8 for div and 10 for grad

    if(kk_>=8)
    {
        mtk::MTK_Real temp;
        mtk::MTK_Real *lamed;
        norm_=CalculateNorm(qq_,kk_);
        minnorm_=1.0e6;
        minrow_=100;
    mtk::MTK_Real *PHI=nullptr;

    try {
        PHI = new MTK_Real[kk_*(kk_+ 1)];
        memset(PHI, 0.0, sizeof(PI_[0])*kk_*(kk_+1));
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }
#if MTK_DEBUG_LEVEL>1

    cout << "+++++++++++ PHI Inicia++++++++++++"<<endl;
    for (auto ii=0; ii<kk_+1;ii++)
    {
        for (auto jj=0;jj<kk_;jj++)
        {
            cout <<setw(12)<< PHI[ii*kk_+jj] << " " ;
        }
    cout << endl;
        
    }
    cout << endl;
    cout << "+++++++++++ PHI++++++++++++"<<endl;
    getchar();
#endif

    // PHI
    for (auto ii = 0; ii < kk_+1; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            PHI[ii*(kk_) + jj] = prem_apps_[ii*dim_null_ + jj];
        }
    }
    
#if MTK_DEBUG_LEVEL>1
    cout << "+++++++++++ PHI after prem apps++++++++++++"<<endl;
    for (auto ii=0; ii<kk_+1;ii++)
    {
        for (auto jj=0;jj<kk_;jj++)
        {
            cout <<setw(12)<< PHI[ii*kk_+jj] << " " ;
        }
    cout << endl;
        
    }
    cout << endl;
    cout << "+++++++++++ PHI++++++++++++"<<endl;
    getchar();
#endif
    auto mm = 0;
    for (auto jj = dim_null_; jj < dim_null_+2; jj++) {
        for (auto ii = 0; ii < kk_; ii++) {
            PHI[(ii + mm)*kk_ + jj] = oo_[ii];
        }
        mm++;
    }
#if MTK_DEBUG_LEVEL>1

    cout << "+++++++++++ PHI after placing oo++++++++++++"<<endl;
    for (auto ii=0; ii<kk_+1;ii++)
    {
        for (auto jj=0;jj<kk_;jj++)
        {
            cout <<setw(12)<< PHI[ii*kk_+jj] << " " ;
        }
    cout << endl;
        
    }
    cout << endl;
    cout << "+++++++++++ PHI++++++++++++"<<endl;

 
#endif 
    

    for(auto jj=kk_-1; jj >=kk_-dim_null_; jj--)
    {
        for(auto ii=0; ii<kk_+1; ii++)
        {
            PHI[ii*kk_+jj]=0.0;

        }
    }
#if MTK_DEBUG_LEVEL>1

    cout << "+++++++++++ PHI after setting zeros++++++++++++"<<endl;
    for (auto ii=0; ii<kk_+1;ii++)
    {
        for (auto jj=0;jj<kk_;jj++)
        {
            cout <<setw(12)<< PHI[ii*kk_+jj] << " " ;
        }
    cout << endl;
        
    }
    cout << endl;
    cout << "+++++++++++ PHI++++++++++++"<<endl;

  getchar();
#endif
 
    // PHI
        for (auto jj = 0; jj < kk_+1; jj++)
        {
            for (auto ii = 0; ii <dim_null_; ii++)
            {
#if MTK_DEBUG_LEVEL>1
                cout << "PHI Index = "<< (ii+kk_-dim_null_+jj*kk_) << "\t";
                cout << "PHI Value = "<< PHI[ii+kk_-dim_null_+jj*kk_ ]<<"\t";
                cout << "Prem_Apps Index = "<< (dim_null_-ii-1+jj*dim_null_) <<"\t";
                cout << "Prem_Apps Value = " << -prem_apps_[(dim_null_-ii-1+jj*dim_null_)]<<endl;
#endif
                PHI[(ii+kk_-dim_null_+jj*kk_)] = -prem_apps_[(dim_null_-ii-1+jj*dim_null_)];

                
            }
        }

#if MTK_DEBUG_LEVEL>1
    cout << "+++++++++++ PHI after negative prem apps++++++++++++"<<endl;
    for (auto ii=0; ii<kk_+1;ii++)
    {
        for (auto jj=0;jj<kk_;jj++)
        {
            cout <<setw(12)<< PHI[ii*kk_+jj] << " " ;
        }
    cout << endl;
        
    }
    cout << endl;
    cout << "+++++++++++ PHI++++++++++++"<<endl;
    getchar();
#endif
    for(auto ii=0; ii<kk_/2; ii++)
    {
        for (auto jj=dim_null_+2; jj<kk_; jj++)
        {
            auto swap=PHI[ii*kk_+jj];
            PHI[ii*kk_+jj]=PHI[(kk_-ii)*kk_+jj];
            PHI[(kk_-ii)*kk_+jj]=swap;
        }
    }

#if MTK_DEBUG_LEVEL>1
    cout << "+++++++++++ PHI after swapping and ready++++++++++++"<<endl;
    for (auto ii=0; ii<kk_+1;ii++)
    {
        for (auto jj=0;jj<kk_;jj++)
        {
            cout <<setw(12)<< PHI[ii*kk_+jj] << " " ;
        }
    cout << endl;
        
    }
    cout << endl;
    cout << "+++++++++++ PHI++++++++++++"<<endl;
    getchar();
#endif    

#if MTK_DEBUG_LEVEL >0
    cout << "Constructed PHI_ for the Optimizer: " << endl;
    for (auto ii = 0; ii < kk_+1; ii++) {
        for (auto jj = 0; jj < kk_; jj++) {
            if (abs(PHI[ii*(kk_) + jj] - 0.0) < mimetic_tol_) {
                cout << setw(11) << 0.0;
            } else {
                cout << setw(11) << PHI[ii*(kk_) + jj];
            }
        }
        cout << endl;
    }
    cout << endl;
#endif

        
      lamed=(mtk::MTK_Real *)malloc(dim_null_*sizeof(mtk::MTK_Real));

       for (auto ii=0; ii<dim_null_; ii++)
       {
        
          lamed[ii]=qq_lp_[ii+kk_+1];
        }


        
#if MTK_DEBUG_LEVEL>0
        cout<<"INSIDE OPTIMIZER"<<endl;

        for(auto ii=0; ii<num_bndy_approxs_; ii++)
        {
            cout << "qq"<<ii<<" = "<< qq_[ii]<<"\t"<<"qq_lp"<<ii<<" = " <<qq_lp_[ii]<<endl;

        }
        cout <<"NULLSPACE :"<<endl;
        for(auto ii=0; ii<dim_null_; ii++)
        {
            for (auto jj=0; jj<num_bndy_approxs_; jj++)
            {
                cout<<NULLS_[ii*num_bndy_approxs_+jj]<<" ";
            }
            cout<<endl;

        }
        cout <<"LAMBDAS:"<<endl;


        for (auto ii=0; ii<dim_null_; ii++)
        {
            cout<<lamed[ii]<<" "<<endl;
        }

        #endif
        for (auto ii=0; ii<num_bndy_approxs_; ii++)
        {
            temp=0;
            for(auto jj=0; jj<dim_null_; jj++)
            {
                temp=temp+lamed[jj]*NULLS_[jj*num_bndy_approxs_+ii];
                
            }
            qq_lp_[ii]=qq_lp_[ii]-temp;

        }
#if MTK_DEBUG_LEVEL>0
        for(auto ii=0; ii<kk_+1; ii++)
        {
            cout << "qq"<<ii<<" = "<< qq_[ii]<<"\t\t"<<"qq_lp"<<ii<<" = "
                 <<qq_lp_[ii]<<endl;

        }
#endif
        // First, it loops over all rows and record all the relative norms
        //  without changing anything yet.
        for(row_=0; row_<kk_+1; row_++)
        {
            normerr_=Optimizer(PHI,kk_+1,
                               kk_,kk_,	qq_lp_,qq_,row_,mimetic_threshold_,0);
            aux_=normerr_/norm_;
#if MTK_DEBUG_LEVEL>0
            cout <<"Relative norm:"<<aux_<<endl<<endl;
#endif
            if(aux_<minnorm_)
            {
                minnorm_ = aux_;
                minrow_=row_;
            }

            
        }
        // After we know which row yields the smallest relative norm
        // that row is chosen to be the objective function
        // and the result of the optimizer is chosen to be the new qq_

#if MTK_DEBUG_LEVEL>0
        cout<<"Minimum Relative Norm "<< minnorm_<<"found at row " << minrow_+1<<endl;
#endif
        normerr_=Optimizer(PHI,kk_+1,
                           kk_,kk_,qq_lp_,qq_,minrow_,mimetic_threshold_,1);
        aux_=normerr_/norm_;
#if MTK_DEBUG_LEVEL>0
            cout <<"Relative norm:"<<aux_<<endl;
#endif
        if (lamed != nullptr) {
        delete[] lamed;
        lamed = nullptr;
        }
    if (PHI != nullptr) {
        delete[] PHI;
        PHI = nullptr;
    }

    }

    if (PI_ != nullptr) {
        delete[] PI_;
        PI_ = nullptr;
    }

    if (PIT_ != nullptr) {
        delete[] PIT_;
        PIT_ = nullptr;
    }
    if (work_ != nullptr) {
        delete[] work_;
        work_ = nullptr;
    }
    return true;
}

// Stage 2.4 of the CRS Algorithm:	
// Compute mimetic stencil approximating at boundary and assemble operator:

bool MTK_1DDiv::ComputeStencilBoundaryGrid(void) {

    // Extract the scalars that have been computed within the weights vector:
    // WARNING: These are already store in memory. We could just adjust this
    // implementation so that it uses those that are in the qq_ array. But for now
    // we will stick to the original algorithm.

    
#if MTK_DEBUG_LEVEL >0
    cout << "lambdas_ =" << endl;
    for (auto ii = 0; ii < dim_null_; ii++) {
        cout << setw(11) << lambdas_[ii] << endl;
    }
    cout << endl;
#endif

    // Compute alpha_ values:
    try {
        alphas_ = new MTK_Real[dim_null_];
        memset(alphas_, 0.0, sizeof(alphas_[0])*dim_null_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    for (auto ii = 0; ii < dim_null_; ii++) {
        alphas_[ii] = lambdas_[ii]/qq_[ii] ;
    }

#if MTK_DEBUG_LEVEL >0
    cout << "alphas_ =" << endl;
    for (auto ii = 0; ii < dim_null_; ii++) {
        cout << setw(11) << alphas_[ii] << endl;
    }
    cout << endl;
#endif

    // Compute the mimetic boundary approximations:
    try {
        mim_bndy_ = new MTK_Real[num_bndy_approxs_*dim_null_];
        memset(mim_bndy_, 0.0, sizeof(mim_bndy_[0])*num_bndy_approxs_*dim_null_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            mim_bndy_[idx(ii,dim_null_,jj)] =
                prem_apps_[idx(ii,dim_null_,jj)] +
                alphas_[jj]*NULLS_[idx(jj,num_bndy_approxs_,ii)];
        }
    }

#if MTK_DEBUG_LEVEL >0
    cout << "Collection of mimetic approximations (Tricky print):" << endl;
    for (auto ii = 0; ii < num_bndy_approxs_; ii++) {
        for (auto jj = 0; jj < dim_null_; jj++) {
            auto value = mim_bndy_[idx(ii,dim_null_,jj)];
            if (abs(value - 0.0) < mimetic_tol_) {
                cout << setw(12) << 0.0;
            } else {
                cout << setw(12) << value;
            }
        }
        cout << endl;
    }
    cout << endl;
#endif

    // Collect garbage :)
//     if (lambdas_ != nullptr) {
//         delete[] lambdas_;
//         lambdas_ = nullptr;
//     }
    if (alphas_ != nullptr) {
        delete[] alphas_;
        alphas_ = nullptr;
    }
    return true;
}

// Stage 3: Final stage. Construct the output array with the operator and its
// weights:

bool MTK_1DDiv::AssembleOperator(void) {

    // The output array will have this form:
    // 1. The first entry of the array will contain the used order kk_.
    // 2. The second entry of the array will contain the collection of
    // approximating coefficients for the interior of the grid.
    // 3. IF kk_ > 2, then the third entry will contain a collection of weights.
    // 4. IF kk_ > 2, the next dim_null_ entries will contain the collections of
    // approximating coefficients for the west boundary of the grid.
    
  if (kk_ > 2) {
        lDD_ = 1 + kk_ + kk_ + dim_null_*num_bndy_approxs_;
    } else {
        lDD_ = 1 + kk_;
    }

    try {
        DD_ = new MTK_Real[lDD_];
        memset(DD_, 0.0, sizeof(DD_[0])*lDD_);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    // 1. The first entry of the array will contain the used order kk_.
    DD_[0] = kk_;

    // 2. The second entry of the array will contain the collection of
    // approximating coefficients for the interior of the grid.
    for (auto ii = 0; ii < kk_; ii++) {
        DD_[ii + 1] = oo_[ii];
    }

    // 3. IF kk_ > 2, then the third entry will contain a collection of weights.
    if (kk_ > 2) {
        for (auto ii = 0; ii < kk_; ii++) {
            DD_[(1 + kk_) + ii] = qq_[ii];
        }
    }

    // 4. IF kk_ > 2, the next dim_null_ entries will contain the collections of
    // approximating coefficients for the west boundary of the grid.

    if (kk_ > 2) {
        auto offset = (1+kk_+kk_);
        int mm =0;
        for (auto ii = 0; ii < dim_null_; ii++) {
            for (auto jj = 0; jj < num_bndy_approxs_; jj++) {
                DD_[offset + (mm)] = mim_bndy_[jj*dim_null_+ii];
                mm++;
            }
        }
    }


    return true;
}


MTK_DenseMatrix* MTK_1DDiv::ReturnAsMatrix(int num_cells, MTK_Real A, MTK_Real B)
{
   MTK_DenseMatrix *XX;

    MTK_Real *ee;
    MTK_Real *ss;
    MTK_Real InvDeltaX;

    int n=num_cells+2;
    int m=num_cells+1;

    int num_extra_rows=kk_/2-1;
    int elements_per_row=3*kk_/2;
    int cc=0;
    int ee_index=0;
    int stencil_size=kk_;

    if (num_cells <3*kk_-1)
    {
       cout<< "Number of cells is too small for the required order!"<<endl;
       return nullptr; 
    }


    XX=new MTK_DenseMatrix(n,m);
    ee= ReturnExtraRows();
    ss= ReturnStencil();
    InvDeltaX=(MTK_Real)(num_cells)/(B-A);

// Filling the upper extra rows:
#if MTK_DEBUG_LEVEL >1
//    for(cc=0; cc<num_extra_elements; cc++)
//        cout << ee[cc] << endl;
#endif
    // Padding first and last rows wirh zeros
    for(auto jj =0; jj<m; jj++)
    {
        XX->SetValue(0,jj,0.0);
        XX->SetValue((n-1),jj,0.0);
    }

    ee_index=0;
    for (auto ii = 1; ii < num_extra_rows +1; ii++)
    {
        cc=0;
        for(auto jj =0 ; jj<m ; jj++)
        {
            if(cc>=elements_per_row)
            {
                XX->SetValue(ii,jj,0.0);
            }
            else
            {
                XX->SetValue(ii,jj,ee[ee_index++]*InvDeltaX);
                cc++;
            }
        }
    }

// Filling the lower extra lows (flipping and negative)
    ee_index=0;
    for (auto ii=n-2; ii >= n - num_extra_rows-1; ii-- )
    {
        cc=0;
        for(auto jj =m-1 ; jj >=0 ; jj--)
        {

            if(cc>=elements_per_row)
            {
                XX->SetValue(ii,jj,0.0);
            }
            else
            {
                XX->SetValue(ii,jj,-ee[ee_index++]*InvDeltaX);
                cc++;
            }
        }
    }
// Filling the Stencil
    for (auto ii=num_extra_rows+1; ii< n-num_extra_rows-1; ii++)
    {
        for (auto jj=0; jj<m; jj++)
        {
            XX->SetValue(ii,jj,0.0);
        }
    }

    for (auto ii=num_extra_rows+1; ii< n-num_extra_rows-1; ii++)
    {
        auto jj = ii-num_extra_rows-1;
        for (cc=0; cc<stencil_size; cc++,jj++)
        {
            XX->SetValue(ii,jj,ss[cc]*InvDeltaX);
        }
    }

// cout<<"Returning DIV as a dense Matrix!"<<endl;
    free (ee);
    free (ss);
    return XX;
}




MTK_Real* MTK_1DDiv::ReturnWeights(void)
{
    MTK_Real *XX;
    XX= (MTK_Real *)malloc(kk_*sizeof(MTK_Real));
    if (XX==nullptr)
    {
        cout <<"Error allocating Memory to return weights"<<endl;
        return nullptr;
    }
    for (auto ii = 0; ii < kk_; ii++)
        XX[ii]=DD_[ii+kk_ +1];
    return XX;
}


MTK_Real* MTK_1DDiv::ReturnStencil(void)
{
    MTK_Real *XX;
    XX= (MTK_Real *)malloc(kk_*sizeof(MTK_Real));
    if (XX==nullptr)
    {
        cout <<"Error allocating Memory to return stencil"<<endl;
        return nullptr;
    }
    for (auto ii = 0; ii < kk_; ii++)
        XX[ii]=DD_[ii+1];
    return XX;
}

MTK_Real* MTK_1DDiv::ReturnExtraRows(void)
{
    MTK_Real *XX;
    int NumberOfValues = (kk_/2 -1)*(3*kk_/2);
    XX= (MTK_Real *)malloc(NumberOfValues*sizeof(MTK_Real));
    if (XX==nullptr)
    {
        cout <<"Error allocating Memory to return extra rows"<<endl;
        return nullptr;
    }
    for (auto ii = 0; ii < NumberOfValues; ii++)
    {
        XX[ii]=DD_[ii+kk_+kk_ +1];

        
    }
    return XX;

}

MTK_DenseMatrix* MTK_1DDiv::ReturnMimeticCoefficients(void)
{
    MTK_DenseMatrix *XX;
//    int NumberOfValues = (kk_/2 -1)*(3*kk_/2);
    XX= new MTK_DenseMatrix(kk_/2 -1, 3*kk_/2);
    if (XX==nullptr)
    {
        cout <<"Error allocating Memory to return extra rows"<<endl;
        return nullptr;
    }
    auto counter =0;
    for (auto ii = 0; ii < kk_/2-1; ii++)
    {
        for(auto jj =0;jj<3*kk_/2; jj++)
        {
               XX->SetValue(ii,jj, DD_[2*kk_+1+counter]);
               counter++;
        }
    }
    return XX;

}


// Transposes a matrix AA with rr rows and cc columns:

MTK_Real* Transpose(MTK_Real *AA, int rr, int cc) {

    MTK_Real* AT {};  // Output transpose matrix.

    if (rr <= 0 || cc <= 0 || AA == nullptr) {
        return nullptr;
    }
    try {
        AT = new MTK_Real[rr*cc];
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

MTK_Real* Vandermonde(const MTK_Real *gen,      // Given generator vector.
                      const int gen_length,   // Length of the generator vector.
                      const int pro_length,   // Length of the progression.
                      const bool transpose) { // Should create its transpose?

    MTK_Real* TT {};  // Output Vandermonde matrix.

    // Check for the integrity of the arguments:
    if (gen == nullptr || gen_length < 1 || pro_length < 1) {
        return nullptr;
    }

    try {
        TT = new MTK_Real[gen_length*pro_length];
        memset(TT, 0.0, sizeof(TT[0])*gen_length*pro_length);
    } catch (bad_alloc &memory_allocation_exception) {
        cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
        cerr << memory_allocation_exception.what() << endl;
    }

    if (!transpose) {
        for (auto ii = 0; ii < gen_length; ii++) {
            for (auto jj = 0; jj < pro_length; jj++) {
                TT[idx(ii,pro_length,jj)] = pow(gen[ii], (MTK_Real) jj);
            }
        }
    } else {
        for (auto ii = 0; ii < pro_length; ii++) {
            for (auto jj = 0; jj < gen_length; jj++) {
                TT[idx(ii,gen_length,jj)] = pow(gen[jj], (MTK_Real) ii);
            }
        }
    }

    return TT;
}

// Prints a dense matrix:

bool DenseMatrix_Print(const MTK_Real* aa,      // The matrix.
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

inline int idx(const int ii, const int offset, const int jj) {

    return ii*offset + jj;
}

MTK_Real CalculateNorm(MTK_Real *A,int n)
{
    int ii;
    MTK_Real sum_sq=0;
    for (ii=0; ii<n; ii++)
    {
        sum_sq=sum_sq+A[ii]*A[ii];
    }
    return sqrt(sum_sq);
}
