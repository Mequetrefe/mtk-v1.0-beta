#include <iostream> // Output.
#include <iomanip>  // MATLAB-like formatted text output.
#include <cmath>    // Use of pow.
#include <cstring>  // Use of memset.
#include <glpk.h>

#include "mtk_roots.h"
#include "mtk_1d_div.h"
#include "mtk_blas_facade.h"
#include "mtk_lapack_facade.h"
#include "mtk_glpk_facade.h"

#define VERBOSE 0 // Should I verbose the execution?
#define MTK_DEBUG_LEVEL 1




CONSTURCTOR {

    if(argc ==2)
  {
    kk=atoi(argv[1]);// Order as argument
  }
  else
  {
    kk=2; //Default
  }

  cout << "kk " << kk << endl;

  if (kk < 2) {
    cerr << "Order of accuracy should be greater or equal than 2." << endl;
    cerr << "Exiting..." << endl;
    return EXIT_FAILURE;
  }
  if ((kk%2) != 0) {
    cerr << "Order of accuracy should be an even number." << endl;
    cerr << "Exiting..." << endl;
    return EXIT_FAILURE;
  }
  if (kk >= 8) {
    cout << "WARNING: Numerical accuracy is too high." << endl;
  }

  // Compute stencil for the interior cells:

  if (!ComputeStencilInteriorGrid()) {
    cerr << "Could NOT complete stage 1." << endl;
    cerr << "Exiting..." << endl;
    return EXIT_FAILURE;
  }

  // At this point, we already have the values for the interior stencil stored
  // in the oo array.

  // It is noteworthy, that the 2nd order accurate divergence operator has NO
  // approximation at the boundary, thus it has no weights. For this case, the
  // dimension of the null-space of the Vandermonde matrices used to compute the
  // approximating coefficients at the boundary is 0. Ergo, we compute this
  // number first and then decide if we must compute anything at the boundary:

  dim_null = kk/2 - 1;

  if (dim_null > 0) {

    num_bndy_approxs = (int) (3.0*((double) kk)/2.0);

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
    // null-space = Q(:, last (kk/2 - 1) columns of Q );
    //
    // However, given the nature of the Vandermonde matrices we've just
    // computed, they all posses the same null-space. Therefore, we impose the
    // convention of computing the null-space of the first Vandermonde matrix
    // (west boundary).

    if (!ComputeScaledNullSpace()) {
      cerr << "Could NOT complete stage 2.1." << endl;
      cerr << "Exiting..." << endl;
      return EXIT_FAILURE;
    }

    // Compute the preliminary approximations for the boundary. These
    // preliminary approximations will then be used to compute the real
    // approximations by means of the weights.

    if (!ComputePreliminaryApproximations()) {
      cerr << "Could NOT complete stage 2.2." << endl;
      cerr << "Exiting..." << endl;
      return EXIT_FAILURE;
    }

    // Assemble the PI matrix using:
    // 1. The collection of scaled preliminary approximations.
    // 2. The collection of coefficients approximating at the interior.
    // 3. The scaled basis for the null-space.

    if (!ComputeWeights()) {
      cerr << "Could NOT complete stage 2.3." << endl;
      cerr << "Exiting..." << endl;
      return EXIT_FAILURE;
    }

    // Compute mimetic boundary approximations:
    if (!ComputeStencilBoundaryGrid()) {
      cerr << "Could NOT complete stage 2.4." << endl;
      cerr << "Exiting..." << endl;
      return EXIT_FAILURE;
    }

    if (NULLS != nullptr) {
      delete[] NULLS;
      NULLS = nullptr;
    }
    if (prem_apps != nullptr) {
      delete[] prem_apps;
      prem_apps = nullptr;
    }
  } // End of: if (dim_null > 0);

  // Once we have the following three collections of data:
  //   (a) the coefficients for the interior,
  //   (b) the coefficients for the boundary (if it applies),
  //   (c) and the weights (if it applies),
  // we will store everything in the output array:

  if (!AssembleOperator()) {
    cerr << "Could NOT complete stage 3." << endl;
    cerr << "Exiting..." << endl;
    return EXIT_FAILURE;
  }

  // Print resulting operator:
  MTK_PrintOperatorArrays();

  // Garbage collection :)
  if (oo != nullptr) {
    delete[] oo;
    oo = nullptr;
  }
  if (qq != nullptr) {
    delete[] qq;
    qq = nullptr;
  }
  if (mim_bndy != nullptr) {
    delete[] mim_bndy;
    mim_bndy = nullptr;
  }
  if (DD != nullptr) {
    delete[] DD;
    DD = nullptr;
  }
}

double
optimizer (double *A, int nrows, int ncols, int kk, double *hh, double *qq, int robjective, double mimetic_tolerance, int copy)
{
  glp_prob *lp;
  int problem_size, lp_nrows, lp_ncols; //Number of rows and columns
  int *ia, *ja;                 // of the LP problem`
  double *ar, obj_value, *objective, *rhs;
  int matsize;
  char rname[5], cname[5];
  char mps_file_name[18];
  int glp_index;
  int ii, jj;
  double x1;
  lp_nrows = kk;                // LP problem size is the same as the order.
  lp_ncols = kk;
  double *err;
  matsize = lp_nrows * lp_ncols;


//IMPORTANT: GLPK indexs from 1, so get the extra space needed.
// Memory allocation
  problem_size = lp_nrows * lp_ncols + 1;

  ia = (int *) malloc (sizeof (int) * (problem_size));

  if (ia == NULL)
  {
    printf ("Problem with malloc\n");
    return 2;
  }
  ja = (int *) malloc (sizeof (int) * (problem_size));
  if (ja == NULL)
  {
   printf ("Problem with malloc\n");
    return 3;
  }
  ar = (double *) malloc (sizeof (double) * (problem_size));
  if (ar == NULL)
  {
    printf ("Problem with malloc\n");
    return 4;
  }
  objective = (double *) malloc (sizeof (double) * (lp_ncols+1));
  if (objective == NULL)
  {
    printf ("Problem with malloc\n");
    return 5;
  }
  rhs = (double *) malloc (sizeof (double) * (lp_nrows + 1));
  if (rhs == NULL)
  {
    printf ("Problem with malloc\n");
    return 6;
  }
  err = (double *) malloc (sizeof (double) * (lp_nrows ));
  if (err == NULL)
  {
    printf ("Problem with malloc\n");
    return 7;
  }
if(MTK_DEBUG_LEVEL>0)
{
  printf ("Problem size = %d, rows = %d cols= %d\n", problem_size,
          lp_nrows, lp_ncols);
}
  lp = glp_create_prob ();
  glp_set_prob_name (lp, "MTK_optimizer");

  glp_set_obj_dir (lp, GLP_MIN);
  /* fill problem */
  glp_add_rows (lp, lp_nrows);
  for (ii = 1; ii <= lp_nrows; ii++)
  {
    sprintf(rname, "R%02d",ii);
    glp_set_row_name (lp, ii, rname);
  }
  glp_add_cols (lp, lp_ncols);

  for (ii = 1; ii <= lp_ncols; ii++)
  {
    sprintf(cname, "Q%02d",ii);
    glp_set_col_name (lp, ii, cname);
  }

      glp_index = 1;            //Index of the objective function`
/* Copy the row to the vector objective */
if(MTK_DEBUG_LEVEL>0)
{
    printf ("\nUsing row %d as objective.\n", robjective + 1);
}
    for (jj = 0; jj < kk; jj++)
    {
      objective[glp_index] = A[jj + robjective * ncols];
      glp_index++;
    }
    if(MTK_DEBUG_LEVEL>0)
    {
    printf ("\n");
    }

//  Forming the RHS
    glp_index = 1;
    rhs[0] = 0.0;
    for (ii = 0; ii <= lp_nrows; ii++)
    {
      if (ii != robjective)
      {
        rhs[glp_index] = hh[ii];
        glp_set_row_bnds (lp, glp_index, GLP_UP, 0.0, rhs[glp_index]);
        glp_index++;
      }
    }

    // Setting up the objective function
    for (ii = 1; ii <= lp_ncols; ii++)
    {
      glp_set_obj_coef (lp, ii, objective[ii]);
    }


// Constrains
    for (ii = 1; ii <= lp_ncols; ii++)
    {
      glp_set_col_bnds (lp, ii, GLP_LO, mimetic_tolerance, 0.0);
    }

    /* Copy the matrix minus the row objective to the glpk problem */
    glp_index = 1;
    for (ii = 0; ii <= kk; ii++)
    {
      for (jj = 0; jj < kk; jj++)
      {
        if (ii != robjective)
        {
          ar[glp_index] = A[jj + ii * ncols];
          glp_index++;
        }
      }
    }


    glp_index = 0;
    for (ii = 1; ii < problem_size; ii++)
    {
      if (((ii - 1) % lp_ncols) == 0)
        glp_index++;
      ia[ii] = glp_index;
      ja[ii] = (ii - 1) % lp_ncols + 1;
    }

    glp_load_matrix (lp, matsize, ia, ja, ar);
    /* solve problem */
    glp_simplex (lp, NULL);
    /* print the LP problem to a MPS file */
    if(MTK_DEBUG_LEVEL>0)
    {
    sprintf(mps_file_name,"LP_MPS_row_%02d.mps",robjective);
    glp_write_mps(lp, GLP_MPS_FILE,NULL,mps_file_name);
    }
    /* recover and display results */
    obj_value = glp_get_obj_val (lp);

  for(ii=1;ii<=lp_ncols;ii++)
  {

    err[ii-1] = qq[ii-1]- glp_get_col_prim(lp,ii);
  }
  if(MTK_DEBUG_LEVEL>0)
  {
    printf ("     CRS\tCBS\n");
    for (ii = 0; ii < lp_ncols; ii++)
    {
      printf ("q%d : %g\t%g\n", ii+1, glp_get_col_prim(lp,ii+1) ,qq[ii]);
    }
    printf ("Objective function value (row %d) = %g\n", robjective + 1,
              obj_value);
  }
 if(copy)
{
  for(ii=0;ii<lp_ncols;ii++)
  {

    qq[ii]= glp_get_col_prim(lp,ii+1);
  }
}

  x1=CalculateNorm(err,lp_ncols);
  /* housekeeping */
  glp_delete_prob (lp);
  glp_free_env ();
  free (ia);
  free (ja);
  free (ar);
  free (objective);
  free (rhs);
  return x1;
}


// Stage 1 of the CRS Algorithm:
// Compute the stencil approximating the interior of the staggered grid:

bool ComputeStencilInteriorGrid(void) {

  // Interior spatial coordinates vector:
  try {
    pp = new double[kk];
    memset (pp, 0.0, sizeof(pp[0])*kk);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  pp[0] = 1.0/2.0 - ((double) kk)/2.0;
  for (auto ii = 1; ii < kk; ii++) {
    pp[ii] = pp[ii - 1] + 1.0;
  }

  if (VERBOSE) {
    cout << "pp =" << endl;
    for (auto ii = 0; ii < kk; ii++) {
      cout << setw(12) << pp[ii];
    }
    cout << endl << endl;
  }

  // Vandermonde matrix (using interior coordinates as generator):
  TT = Vandermonde(pp, kk, kk, false);
  if (TT == nullptr) {
    cerr << "ERROR constructing matrix at line " << __LINE__ - 2 << endl;
    cerr << "Exiting..." << endl;
    return false;
  }

  if (VERBOSE) {
    cout << "TT =" << endl;
    if (!MTK_DenseMatrix_Print(TT, kk, kk)) {
      return false;
    }
  }

  // Order-selector vector:
  try {
    oo = new double[kk];
    memset (oo, 0.0, sizeof(oo[0])*kk);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  oo[1] = 1.0;

  if (VERBOSE) {
    cout << "oo =" << endl;
    for (auto ii = 0; ii < kk; ii++) {
      cout << setw(12) << oo[ii];
    }
    cout << endl << endl;
  }

  // Solve systems of equations to obtain the stencil:
  try {
    ipiv = new int[kk];
    memset(ipiv, 0, sizeof(ipiv[0])*kk);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }
  ldoo=kk;
  TT_ld=kk;
  dgesv_(&kk, &nrhs, TT, &TT_ld, ipiv, oo, &ldoo, &info);
  if (!info) {
    if (VERBOSE) {
      cout << "System successfully solved! Interior stencil attained!" << endl;
    }
  } else {
    cerr << "Something went wrong solving system! info = " << info << endl;
    cerr << "Exiting..." << endl;
    return false;
  }

  if (VERBOSE) {
    cout << endl;
    cout << "ss =" << endl;
    for (auto ii = 0; ii < kk; ii++) {
      cout << setw(12) << oo[ii];
    }
    cout << endl << endl;
  }

  // Garbage collection :)
  if (pp != nullptr) {
    delete[] pp;
    pp = nullptr;
  }
  if (TT != nullptr) {
    delete[] TT;
    TT = nullptr;
  }
  if (ipiv != nullptr) {
    delete[] ipiv;
    ipiv = nullptr;
  }

  return true;
}

// Stage 2.1 of the CRS Algorithm:
// Compute a scaled basis of the null-space of the Vandermonde matrix
// approximating at the west boundary:

bool ComputeScaledNullSpace(void) {

  // Generator vector for the first approximation:
  try {
    gg = new double[num_bndy_approxs];
    memset(gg, 0.0, sizeof(gg[0])*num_bndy_approxs);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  gg[0] = -1.0/2.0;
  for (auto ii = 1; ii < num_bndy_approxs; ii++) {
    gg[ii] = gg[ii - 1] + 1.0;
  }

  if (VERBOSE) {
    cout << "gg_west =" << endl;
    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      cout << setw(12) << gg[ii];
    }
    cout << endl << endl;
  }

  // Construct the TRANSPOSE of the first Vandermonde matrix:
  aa_num_rows = kk + 1;
  aa_num_cols = num_bndy_approxs;

  try {
    aa = new double[aa_num_rows*aa_num_cols];
    memset(aa, 0.0, sizeof(aa[0])*(aa_num_rows*aa_num_cols));
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  for (auto ii = 0; ii < aa_num_rows; ii++) {
    for (auto jj = 0; jj < aa_num_cols; jj++) {
      aa[ii*aa_num_cols + jj] = pow(gg[jj], ii);
    }
  }

  if (VERBOSE) {
    cout << "aa_west = " << endl;
    for (auto ii = 0; ii < aa_num_rows; ii++) {
      for (auto jj = 0; jj < aa_num_cols; jj++) {
        cout << setw(12) << aa[ii*aa_num_cols + jj];
      }
      cout << endl;
    }
    cout << endl;
  }

  aaT_num_rows = aa_num_cols;
  aaT_num_cols = aa_num_rows;

  // Factorize aa:

  // Prepare to factorize: allocate and inquire for the value of lwork:
  try {
    work = new double[1];
    memset(work, 0.0, sizeof(aa[0])*1);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  lwork = -1;
  info = 0;

  dgeqrf_(&aaT_num_rows, &aaT_num_cols, aa, &aaT_num_rows, tau,
          work, &lwork, &info);

  if (info == 0) {
    lwork = (int) work[0];
  } else {
    cerr << "Could not get value for lwork on line " << __LINE__ - 5 << endl;
    cerr << "Exiting..." << endl;
    return EXIT_FAILURE;
  }

  if (VERBOSE) {
    cout << "lwork = " << endl << setw(12)<< lwork << endl << endl;
  }

  if (work != nullptr) {
    delete[] work;
    work = nullptr;
  }

  // Once we know lwork, we can actually invoke the factorization:
  try {
    work = new double [lwork];
    memset(work, 0.0, sizeof(work[0])*lwork);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  ltau = min(aaT_num_rows,aaT_num_cols);

  try {
    tau = new double [ltau];
    memset(tau, 0.0, sizeof(0.0)*ltau);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  dgeqrf_(&aaT_num_rows, &aaT_num_cols, aa, &aaT_num_rows,
          tau, work, &lwork, &info);

  if (!info) {
    if (VERBOSE) {
      cout << "QR factorization successfully completed!" << endl << endl;
    }
  } else {
    cerr << "Something went wrong solving system! info = " << info << endl;
    cerr << "Exiting..." << endl;
    return EXIT_FAILURE;
  }

  if (VERBOSE) {
    cout << "aa_west AFTER factorization:" << endl;
    for (auto ii = 0; ii < aa_num_rows; ii++) {
      for (auto jj = 0; jj < aa_num_cols; jj++) {
        cout << setw(13) << aa[ii*aa_num_cols + jj];
      }
      cout << endl;
    }
    cout << endl;
  }

  // We now generate the real matrix Q with orthonormal columns. This has to
  // be done separately since the actual output of dgeqrf_ (AA) represents
  // the orthogonal matrix Q as a product of min(aa_num_rows,aa_num_cols)
  // elementary Householder reflectors. Notice that we must re-inquire the new
  // value for lwork that is used.
  try {
    QQ = new double[num_bndy_approxs*num_bndy_approxs];
    memset(QQ, 0.0, sizeof(QQ[0])*num_bndy_approxs*num_bndy_approxs);
    } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  // QQ must be initialized to the identity matrix:
  for (auto ii = 0; ii < num_bndy_approxs; ii++) {
    for (auto jj = 0; jj < num_bndy_approxs; jj++) {
      QQ[ii*num_bndy_approxs + jj] = (ii == jj)? 1.0: 0.0;
    }
  }

  if (VERBOSE) {
    cout << "Initialized QQ_T:" << endl;
    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      for (auto jj = 0; jj < num_bndy_approxs; jj++) {
        cout << setw(13) << QQ[ii*num_bndy_approxs + jj];
      }
      cout << endl;
    }
    cout << endl;
  }

  // Assemble the QQ matrix:
  lwork = -1;

  if (work != nullptr) {
    delete[] work;
    work = nullptr;
  }

  try {
    work = new double[1];
    memset(work, 0.0, sizeof(work[0])*1);
    } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  dormqr_(&side, &trans, &aa_num_cols, &aa_num_cols, &ltau, aa, &aaT_num_rows,
          tau, QQ, &num_bndy_approxs, work, &lwork, &info);

  if (info == 0) {
    lwork = (int) work[0];
  } else {
    cerr << "Could not get value for lwork on line " << __LINE__ - 5 << endl;
    cerr << "Exiting..." << endl;
    return EXIT_FAILURE;
  }

  if (VERBOSE) {
    cout << "lwork = " << endl << setw(12)<< lwork << endl << endl;
  }

  if (work != nullptr) {
    delete[] work;
    work = nullptr;
  }

  try {
    work = new double[lwork];
    memset(work, 0.0, sizeof(work[0])*lwork);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  dormqr_(&side, &trans, &aa_num_cols, &aa_num_cols, &ltau, aa, &aaT_num_rows,
          tau, QQ, &num_bndy_approxs, work, &lwork, &info);

  if (!info) {
    if (VERBOSE) {
      cout << "Q matrix successfully assembled!" << endl << endl;
    }
  } else {
    cerr << "Something went wrong solving system! info = " << info << endl;
    cerr << "Exiting..." << endl;
    return EXIT_FAILURE;
  }

  if (VERBOSE) {
    cout << "QQ_T =" << endl;
    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      for (auto jj = 0; jj < num_bndy_approxs; jj++) {
        cout << setw(13) << QQ[ii*num_bndy_approxs + jj];
      }
      cout << endl;
    }
    cout << endl;
  }

  // Extract the null-space from QQ, and save it in KK:
  KK_num_rows = num_bndy_approxs;
  KK_num_cols = dim_null;

  try {
    KK = new double[KK_num_rows*KK_num_cols];
    memset(KK, 0.0, sizeof(KK[0])*(KK_num_rows*KK_num_cols));
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  for (auto ii = (num_bndy_approxs - dim_null); ii < num_bndy_approxs; ii++) {
    for (auto jj = 0; jj < num_bndy_approxs; jj++) {
      KK[idx(jj,dim_null,(ii - (num_bndy_approxs - dim_null)))] =
        QQ[idx(ii,num_bndy_approxs,jj)];
    }
  }

  if (VERBOSE) {
    cout << "KK =" << endl;
    for (auto ii = 0; ii < KK_num_rows; ii++) {
      for (auto jj = 0; jj < KK_num_cols; jj++) {
        cout << setw(15) << KK[idx(ii,KK_num_cols,jj)];
      }
      cout << endl;
    }
    cout << endl;
  }

  // Scale KK:

  // Scale thus requesting that the last entries of the attained basis for the
  // null-space, adopt the pattern we require.
  // Essentially we will implement the following MATLAB pseudo-code:
  //  scalers = KK(num_bndy_approxs - (dim_null - 1):num_bndy_approxs,:)\B
  //  SK = KK*scalers
  // where SK is the scaled null-space.

  // In this point, we almost have all the data we need correctly allocated
  // in memory. We will create the matrix II, and elements we wish to scale in
  // the KK array. Using the concept of the leading dimension, we could just
  // use KK, with the correct leading dimension and that is it. BUT I DO NOT
  // GET how does it work. So I will just create a matrix with the content of
  // this array that we need, solve for the scalers and then scale the
  // whole KK:

  // We will then create memory for that sub-matrix of KK (SUBK):
  try {
    SUBK = new double[dim_null*dim_null];
    memset(SUBK, 0.0, sizeof(SUBK[0])*dim_null*dim_null);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  // Extracting information from KK:
  for (auto ii = num_bndy_approxs - dim_null; ii < num_bndy_approxs; ii++) {
    for (auto jj = 0; jj < dim_null; jj++) {
      SUBK[idx(ii - (num_bndy_approxs - dim_null),dim_null,jj)] =
        KK[idx(ii,dim_null,jj)];
    }
  }

  if (VERBOSE) {
    cout << "Required SUBK for scaling:" << endl;
    for (auto ii = 0; ii < dim_null; ii++) {
      for (auto jj = 0; jj < dim_null; jj++) {
        cout << setw(15) << SUBK[idx(ii,dim_null,jj)];
      }
      cout << endl;
    }
    cout << endl;
  }

  // Once we have the matrix for the scaling, create the collection of rhss,
  // which looks like an identity matrix.
  try {
    II = new double[dim_null*dim_null];
    memset(II, 0.0, sizeof(II[0])*dim_null*dim_null);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  for (auto ii = 0; ii < dim_null; ii++) {
    for (auto jj = 0; jj < dim_null; jj++) {
      II[ii*dim_null + jj] = (ii == jj)? 1.0: 0.0;
    }
  }

  if (VERBOSE) {
    cout << "Required II for scaling:" << endl;
    for (auto ii = 0; ii < dim_null; ii++) {
      for (auto jj = 0; jj < dim_null; jj++) {
        cout << setw(15) << II[ii*dim_null + jj];
      }
      cout << endl;
    }
    cout << endl;
  }

  // Solve the system to compute the scalers.
  // An example of the system to solve, for k = 8, is:
  //
  // SUBK*scalers = II or
  //
  // |  0.386018  -0.0339244   -0.129478 |           | 1 0 0 |
  // | -0.119774   0.0199423   0.0558632 |*scalers = | 0 1 0 |
  // | 0.0155708 -0.00349546 -0.00853182 |           | 0 0 1 |
  //
  // Notice this is a nrhs = 3 system.
  // Noteworthy: we do NOT ACTUALLY ALLOCATE space for the scalers... they
  // will be stored in the created identity matrix.
  // Let us first transpose SUBK (because of LAPACK):

  SUBKT = Transpose(SUBK, dim_null, dim_null);
  if (SUBKT == nullptr) {
    cerr << "ERROR constructing  matrix at line " << __LINE__ - 2 << endl;
    cerr << "Exiting..." << endl;
    return false;
  }

  if (VERBOSE) {
    cout << "Required SUBK transposed for scaling:" << endl;
    for (auto ii = 0; ii < dim_null; ii++) {
      for (auto jj = 0; jj < dim_null; jj++) {
        cout << setw(15) << SUBKT[idx(ii,dim_null,jj)];
      }
      cout << endl;
    }
    cout << endl;
  }

  // All right! Ready to solve using the LAPACK:
  nrhs = dim_null;
  ldSUBKT = dim_null;
  info = 0;

  try {
    ipiv = new int[ldSUBKT];
    memset(ipiv, 0, sizeof(ipiv[0])*ldSUBKT);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  dgesv_(&ldSUBKT, &nrhs, SUBKT, &ldSUBKT, ipiv, II, &ldSUBKT, &info);
  if (!info) {
    if (VERBOSE) {
      cout << "System successfully solved! Scalers attained!" << endl << endl;
    }
  } else {
    cerr << "Something went wrong solving system! info = " << info << endl;
    return EXIT_FAILURE;
  }

  if (VERBOSE) {
    cout << "Computed scalers:" << endl;
    for (auto ii = 0; ii < dim_null; ii++) {
      for (auto jj = 0; jj < dim_null; jj++) {
        cout << setw(15) << II[idx(ii,dim_null,jj)];
      }
      cout << endl;
    }
    cout << endl;
  }

  // We now multiply the two matrices to attain a scaled basis for null-space:

  // Transpose the basis matrix KK:
  KKT = Transpose(KK, num_bndy_approxs, dim_null);
  if (KKT == nullptr) {
    cerr << "ERROR constructing matrix at " << __LINE__ - 2 << endl;
    cerr << "Exiting..." << endl;
    return false;
  }

  // Allocate space for the product matrix:
  try {
    NULLS = new double[num_bndy_approxs*dim_null];
    memset(NULLS, 0.0, sizeof(NULLS[0])*num_bndy_approxs*dim_null);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  // Perform matrix-matrix multiplication:
  dgemm_(&transa, &transb, &num_bndy_approxs, &dim_null, &dim_null,
          &alpha, KKT, &num_bndy_approxs, II, &dim_null, &beta,
          NULLS, &num_bndy_approxs);
  if (VERBOSE) {
    cout << "Scaled basis for the null-space (tricky print):" << endl;
    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      for (auto jj = 0; jj < dim_null; jj++) {
        auto value = NULLS[idx(jj,num_bndy_approxs,ii)];
        if (abs(value - 0.0) < mtk_tol) {
          cout << setw(12) << 0.0;
        } else {
          cout << setw(12) << value;
        }
      }
      cout << endl;
    }
    cout << endl;
  }

  // At this point, we have a scaled basis for the null-space, with the
  // pattern we need!

  // Local garbage collection :)
  if (gg != nullptr) {
    delete[] gg;
    gg = nullptr;
  }
  if (aa != nullptr) {
    delete[] aa;
    aa = nullptr;
  }
  if (work != nullptr) {
    delete[] work;
    work = nullptr;
  }
  if (tau != nullptr) {
    delete[] tau;
    tau = nullptr;
  }
  if (QQ != nullptr) {
    delete[] QQ;
    QQ = nullptr;
  }
  if (KK != nullptr) {
    delete[] KK;
    KK = nullptr;
  }
  if (SUBK != nullptr) {
    delete[] SUBK;
    SUBK = nullptr;
  }
  if (II != nullptr) {
    delete[] II;
    II = nullptr;
  }
  if (SUBKT != nullptr) {
    delete[] SUBKT;
    SUBKT = nullptr;
  }
  if (ipiv != nullptr) {
    delete[] ipiv;
    ipiv = nullptr;
  }
  if (KKT != nullptr) {
    delete[] KKT;
    KKT = nullptr;
  }

  return true;
}

// Stage 2.2 of the CRS Algorithm:
// Compute the set of preliminary approximation on the boundary neighborhood:

bool ComputePreliminaryApproximations(void) {

  // Generator vector for the first approximation:
  try {
    gg = new double[num_bndy_approxs];
    memset(gg, 0.0, sizeof(gg[0])*num_bndy_approxs);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  gg[0] = -1.0/2.0;
  for (auto ii = 1; ii < num_bndy_approxs; ii++) {
    gg[ii] = gg[ii - 1] + 1.0;
  }

  if (VERBOSE) {
    cout << "gg_0 =" << endl;
    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      cout << setw(12) << gg[ii];
    }
    cout << endl << endl;
  }

  // Allocate 2D array to store the collection of preliminary approximations:
  try {
    prem_apps = new double[num_bndy_approxs*dim_null];
    memset(prem_apps, 0.0, sizeof(prem_apps[0])*num_bndy_approxs*dim_null);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  // Compute the dim_null near-the-boundary columns of the PI matrix:
  for (auto ll = 0; ll < dim_null; ll++) {

    // Re-check new generator vector for every iteration except for the first:
    if (ll > 0) {
      if (VERBOSE) {
        cout << "gg_" << ll << " =" << endl;
        for (auto ii = 0; ii < num_bndy_approxs; ii++) {
          cout << setw(12) << gg[ii];
        }
        cout << endl << endl;
      }
    }

    // Create the Vandermonde matrix for this iteration:
    AA = Vandermonde(gg, num_bndy_approxs, kk + 1, false);
    if (AA == nullptr) {
      cerr << "ERROR constructing matrix at line " << __LINE__ - 2 << endl;
      cerr << "Exiting..." << endl;
      return EXIT_FAILURE;
    }

    AA_num_rows = kk + 1;
    AA_num_cols = num_bndy_approxs;
    AA_ld = AA_num_rows;

    if (VERBOSE) {
      // Although it is generated transposed, we print it non-transposed:
      cout << "AA^T_" << ll << " =" << endl;
      if (!MTK_DenseMatrix_Print(AA, AA_num_rows, AA_num_cols)) {
        return EXIT_FAILURE;
      }
    }

    // New order-selector vector (it gets re-written with LAPACK solutions):
    ob_ld = num_bndy_approxs;

    try {
      ob = new double[ob_ld];
      memset(ob, 0.0, sizeof(ob[0])*ob_ld);
    } catch (bad_alloc &memory_allocation_exception) {
      cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
      cerr << memory_allocation_exception.what() << endl;
    }

    ob[1] = 1.0;

    if (VERBOSE) {
      cout << "ob = " << endl << endl;
      for (auto ii = 0; ii < ob_ld; ii++) {
        cout << setw(12) << ob[ii] << endl;
      }
      cout << endl;
    }

    // Solving TT*rr = ob yields the columns rr of the KK matrix. However,
    // this is an under-determined system of equations. So we can not use the
    // same LAPACK routine (dgesv_). We will instead use dgels_.
    // We first invoke the solver to query for the value of lwork. For this,
    // we must at least allocate enough space to allow access to WORK(1), or
    // work[0]:
    // If LWORK = -1, then a workspace query is assumed; the routine only
    // calculates the optimal size of the WORK array, returns this value as
    // the first entry of the WORK array, and no error message related to
    // LWORK is issued by XERBLA.
    if (work != nullptr) {
      delete[] work;
      work = nullptr;
    }
    try {
      work = new double[1];
      memset(work, 0.0, sizeof(work[0])*1);
    } catch (bad_alloc &memory_allocation_exception) {
      cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
      cerr << memory_allocation_exception.what() << endl;
    }

    // IMPORTANT: Never forget to re-set the nrhs variable to the right
    // number:
    nrhs = 1;
    info = 0;
    lwork = -1;

    dgels_(&trans, &AA_num_rows, &AA_num_cols, &nrhs, AA, &AA_ld, ob, &ob_ld,
            work, &lwork, &info);

    if (info == 0) {
      lwork = (int) work[0];
    } else {
      cerr << "Could not get value for lwork on line " << __LINE__ - 2 << endl;
      cerr << "Exiting..." << endl;
      return EXIT_FAILURE;
    }

    if (VERBOSE) {
      cout << "lwork = " << endl << setw(12)<< lwork << endl << endl;
    }

    // We then use lwork's new value to create the work array:
    if (work != nullptr) {
      delete[] work;
      work = nullptr;
    }

    try {
      work = new double[lwork];
      memset(work, 0.0, sizeof(work[0])*lwork);
    } catch (bad_alloc &memory_allocation_exception) {
      cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
      cerr << memory_allocation_exception.what() << endl;
    }

    // We now invoke the solver again:
    dgels_(&trans, &AA_num_rows, &AA_num_cols, &nrhs, AA, &AA_ld, ob, &ob_ld,
            work, &lwork, &info);

    if (!info) {
      if (VERBOSE) {
        cout << "System successfully solved!" << endl << endl;
      }
    } else {
      cerr << "Something went wrong solving system! info = " << info << endl;
    }
    if (VERBOSE) {
      cout << "ob =" << endl;
      for (auto ii = 0; ii < ob_ld; ii++) {
        cout << setw(12) << ob[ii] << endl;
      }
      cout << endl;
    }

    // Prior to save this solution, we scale it, using the scaled basis for
    // the null-space. This implies a DAXPY operation. However, we must
    // construct the arguments for this operation:
    // First, we extract the last dim_null values of the pre-scaled ob, and
    // save them into the ob_bottom array:

    try {
      ob_bottom = new double[dim_null];
      memset(ob_bottom, 0.0, sizeof(ob_bottom[0])*dim_null);
    } catch (bad_alloc &memory_allocation_exception) {
      cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
      cerr << memory_allocation_exception.what() << endl;
    }

    for (auto ii = 0; ii < dim_null; ii++) {
      ob_bottom[(dim_null - 1) - ii] = ob[num_bndy_approxs - ii - 1];
    }

    if (VERBOSE) {
      cout << "ob_bottom =" << endl;
      for (auto ii = 0; ii < dim_null; ii++) {
        cout << setw(12) << ob_bottom[ii] << endl;
      }
      cout << endl << endl;
    }

    // Once we posses the bottom elements, we proceed with the scaling. We
    // must computed an scaled ob, sob, using the scaled null-space NULLS.
    // Such operation is: sob = ob - NULLS*ob_bottom
    // or:                 ob = -1.0*NULLS*ob_bottom + 1.0*ob
    // thus:                Y =    a*A    *x         +   b*Y (DAXPY).
    // Let us invoke this from the BLAS:
    trans = 'N';
    alpha = -1.0;
    beta = 1.0;
    dgemv_(&trans, &num_bndy_approxs, &dim_null,
            &alpha, NULLS, &num_bndy_approxs,
            ob_bottom, &incx, &beta, ob, &incx);
    if (VERBOSE) {
      cout << "scaled ob (tricky print):" << endl;
      for (auto ii = 0; ii < num_bndy_approxs; ii++) {
        auto value = ob[ii];
        if (abs(value - 0.0) < mtk_tol) {
          cout << setw(12) << 0.0 << endl;
        } else {
          cout << setw(12) << value << endl;
        }
      }
      cout << endl;
    }

    // We save the recently scaled solution, into an array containing these.
    // We can NOT start building the PI matrix, simply because I want that part
    // to be separated since its construction depends on the algorithm we want
    // to implement.

    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      prem_apps[ii*dim_null + ll] = ob[ii];
    }

    // After the first iteration, simply shift the entries of the last
    // generator vector used:
    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      gg[ii]--;
    }

    // Garbage collection for this loop:
    if (AA != nullptr) {
      delete[] AA;
      AA = nullptr;
    }
    if (ob != nullptr) {
      delete[] ob;
      ob = nullptr;
    }
    if (work != nullptr) {
      delete[] work;
      work = nullptr;
    }
    if (ob_bottom != nullptr) {
      delete[] ob_bottom;
      ob_bottom = nullptr;
    }
  } // End of: for (ll = 0; ll < dim_null; ll++);

  if (VERBOSE) {
    cout << "Matrix post-scaled preliminary apps (tricky print): " << endl;
    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      for (auto jj = 0; jj < dim_null; jj++) {
        auto value = prem_apps[idx(ii,dim_null,jj)];
        cout << setw(12) << ((abs(value - 0.0) < mtk_tol)? 0.0: value);
      }
      cout << endl;
    }
    cout << endl;
  }

  if (gg != nullptr) {
    delete[] gg;
    gg = nullptr;
  }

  return true;
}

// Stage 2.3 of the CRS Algorithm:
// Assemble the PI matrix and compute the weights:

bool ComputeWeights(void) {

  // Collect the preliminary approximations on the PI matrix:
  try {
    PI = new double[num_bndy_approxs*(num_bndy_approxs - 1)];
    memset(PI, 0.0, sizeof(PI[0])*num_bndy_approxs*(num_bndy_approxs - 1));
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  // We start constructing the PI matrix...

  // We first process the queue of scaled near-the-boundary preliminary
  // solutions, which are queued in scaled_solutions. Each one of these,
  // will be a column of the PI matrix:
  for (auto ii = 0; ii < num_bndy_approxs; ii++) {
    for (auto jj = 0; jj < dim_null; jj++) {
      PI[ii*(2*dim_null + (kk/2 + 1)) + jj] = prem_apps[ii*dim_null + jj];
    }
  }

  // We now add the columns from the known stencil approximating at the
  // interior, however, these must be padded by zeros, according to their
  // position in the final PI matrix:
  auto mm = 0;
  for (auto jj = dim_null; jj < kk; jj++) {
    for (auto ii = 0; ii < kk; ii++) {
      PI[(ii + mm)*(2*dim_null + (kk/2 + 1)) + jj] = oo[ii];
    }
    mm++;
  }

  // We finally add the final set of columns: the scaled basis for the
  // null-space:
  for (auto jj = dim_null + (kk/2 + 1); jj < num_bndy_approxs - 1; jj++) {
    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      PI[ii*(2*dim_null + (kk/2 + 1)) + jj] =
        NULLS[(jj - (dim_null + (kk/2 + 1)))*num_bndy_approxs + ii];
    }
  }

  cout.precision(4);

  if (VERBOSE) {
    cout << endl;
    cout << "ss =" << endl;
    for (auto ii = 0; ii < kk; ii++) {
      cout << setw(11) << oo[ii];
    }
    cout << endl << endl;
  }

  PI_num_cols = num_bndy_approxs - 1;

  if (VERBOSE) {
    cout << "Constructed PI matrix with CRS Algorithm: " << endl;
    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      for (auto jj = 0; jj < PI_num_cols; jj++) {
        if (abs(PI[ii*PI_num_cols + jj] - 0.0) < mtk_tol) {
          cout << setw(11) << 0.0;
        } else {
          cout << setw(11) << PI[ii*PI_num_cols + jj];
        }
      }
      cout << endl;
    }
    cout << endl;
  }

  // Once we have constructed the matrix, we use the interior stencil to
  // build the proper RHS vector h that imposes the mimetic condition.
  // That is, we construct a system PI*q = h, to solve for the weights:
  try {
    qq = new double[num_bndy_approxs];
    memset(qq, 0.0, sizeof(qq[0])*num_bndy_approxs);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }
///********** RAUL: Adding allocation for qq_lp (optimizer)
  try {
    qq_lp = new double[num_bndy_approxs];
    memset(qq_lp, 0.0, sizeof(qq_lp[0])*num_bndy_approxs);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }
//*****************************

  qq[0] = -1.0;
  for (auto ii = (kk/2 + 2 - 1); ii < num_bndy_approxs; ii++) {
    auto aux_xx = 0.0;
    for (auto jj = 0; jj < ((ii - (kk/2 - 1)) - 1); jj++) {
      aux_xx += oo[jj];
    }
    qq[ii] = -1.0*aux_xx;
  }
// RAUL: Copy qq to qq_lp for the optimizer
  for(auto ii =0; ii <num_bndy_approxs;ii++)
  {
    qq_lp[ii]=qq[ii];
  }
//**************************

  if (VERBOSE) {
    cout << "hh =" << endl;
    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      cout << setw(11) << qq[ii] << endl;
    }
    cout << endl;
  }

  // Since we intend to solve, we must transpose PI.
  // TODO: We could create PI, already transposed.
  PIT = Transpose(PI, num_bndy_approxs, PI_num_cols);
  if (PIT == (double*) NULL) {
    cerr << "Could not create matrix at line " << __LINE__ - 2 << endl;
    cerr << "Exiting..." << endl;
    return EXIT_FAILURE;
  }

  if (VERBOSE) {
    cout << "Required PIT Transposed for solving for weights:" << endl;
    for (auto ii = 0; ii < PI_num_cols; ii++) {
      for (auto jj = 0; jj < num_bndy_approxs; jj++) {
        if (abs(PIT[ii*num_bndy_approxs + jj] - 0.0) < mtk_tol) {
          cout << setw(12) << 0.0;
        } else {
          cout << setw(12) << PIT[ii*num_bndy_approxs + jj];
        }
      }
      cout << endl;
      }
    cout << endl;
  }

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
    work = new double[1];
    memset(work, 0.0, sizeof(aa[0])*1);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  nrhs = 1;
  lwork = -1;

  dgels_(&trans, &num_bndy_approxs, &PI_num_cols, &nrhs, PIT, &num_bndy_approxs,
         qq, &num_bndy_approxs, work, &lwork, &info);

  if (info == 0) {
    lwork = (int) work[0];
  } else {
    cerr << "Could not get value for lwork on line " << __LINE__ - 5 << endl;
    cerr << "Exiting..." << endl;
    return EXIT_FAILURE;
  }

  if (VERBOSE) {
    cout << "lwork = " << endl << setw(12)<< lwork << endl << endl;
  }

  if (work != nullptr) {
    delete[] work;
    work = nullptr;
  }

  try {
    work = new double[lwork];
    memset(work, 0.0, sizeof(work[0])*lwork);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  // We now invoke the solver again:
  trans = 'N';
  nrhs = 1;
  info = 0;
  dgels_(&trans, &num_bndy_approxs, &PI_num_cols, &nrhs, PIT, &num_bndy_approxs,
         qq, &num_bndy_approxs, work, &lwork, &info);

  if (!info) {
    if (VERBOSE) {
      cout << "System successfully solved!" << endl << endl;
    }
  } else {
    cerr << "Something went wrong solving system! info = " << info <<
    endl;
  }
  if (VERBOSE) {
    cout << "qq_CRS =" << endl;
    for (auto ii = 0; ii < PI_num_cols; ii++) {
      cout << setw(11) << qq[ii] << endl;
    }
    cout << endl;
  }

//  This part should be executed only if kk > 8 for div and 10 for grad

if(kk>=8)
{  
norm=CalculateNorm(qq,kk);
minnorm=1.0e6;
minrow=100;

// First, it loops over all rows and record all the relative norms
//  without changing anything yet.

for(row=0;row<kk;row++)
{
  normerr=optimizer(PI,num_bndy_approxs, PI_num_cols,kk,qq_lp,qq,row,mimetic_tolerance,0);
  aux=normerr/norm;
  if(MTK_DEBUG_LEVEL >0)
  {
  printf("Relative norm: %g \n",aux);
  }
  if(aux<minnorm)
   {
     minnorm = aux;
     minrow=row;
   }
}
// After we know which row yields the smallest relative norm
// that row is chosen to be the objective function
// and the result of the optimizer is chosen to be the new qq

if(MTK_DEBUG_LEVEL >0)
{
printf("Minimum Relative Norm %g found at row %d:\n", minnorm,minrow+1);
}
normerr=optimizer(PI,num_bndy_approxs, PI_num_cols,kk,qq_lp,qq,minrow,mimetic_tolerance,1);
  aux=normerr/norm;
  if(MTK_DEBUG_LEVEL>0)
  {
  printf("Relative norm: %g \n",aux);
  }
}

  if (PI != nullptr) {
    delete[] PI;
    PI = nullptr;
  }
  if (PIT != nullptr) {
    delete[] PIT;
    PIT = nullptr;
  }
  if (work != nullptr) {
    delete[] work;
    work = nullptr;
  }
  return true;
}

// Stage 2.4 of the CRS Algorithm:
// Compute mimetic stencil approximating at boundary and assemble operator:

bool ComputeStencilBoundaryGrid(void) {

  // Extract the scalars that have been computed within the weights vector:
  // WARNING: These are already store in memory. We could just adjust this
  // implementation so that it uses those that are in the qq array. But for now
  // we will stick to the original algorithm.
  try {
    lambdas = new double[dim_null];
    memset(lambdas, 0.0, sizeof(lambdas[0])*dim_null);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  for (auto ii = 0; ii < dim_null; ii++) {
    lambdas[ii] = qq[kk + ii] ;
  }

  if (VERBOSE) {
    cout << "lambdas =" << endl;
    for (auto ii = 0; ii < dim_null; ii++) {
      cout << setw(11) << lambdas[ii] << endl;
    }
    cout << endl;
  }

  // Compute alpha values:
  try {
    alphas = new double[dim_null];
    memset(alphas, 0.0, sizeof(alphas[0])*dim_null);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  for (auto ii = 0; ii < dim_null; ii++) {
    alphas[ii] = lambdas[ii]/qq[ii] ;
  }

  if (VERBOSE) {
    cout << "alphas =" << endl;
    for (auto ii = 0; ii < dim_null; ii++) {
      cout << setw(11) << alphas[ii] << endl;
    }
    cout << endl;
  }

  // Compute the mimetic boundary approximations:
  try {
    mim_bndy = new double[num_bndy_approxs*dim_null];
    memset(mim_bndy, 0.0, sizeof(mim_bndy[0])*num_bndy_approxs*dim_null);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  for (auto ii = 0; ii < num_bndy_approxs; ii++) {
    for (auto jj = 0; jj < dim_null; jj++) {
      mim_bndy[idx(ii,dim_null,jj)] =
        prem_apps[idx(ii,dim_null,jj)] +
          alphas[jj]*NULLS[idx(jj,num_bndy_approxs,ii)];
    }
  }

  if (VERBOSE) {
    cout << "Collection of mimetic approximations (tricky print):" << endl;
    for (auto ii = 0; ii < num_bndy_approxs; ii++) {
      for (auto jj = 0; jj < dim_null; jj++) {
        auto value = mim_bndy[idx(ii,dim_null,jj)];
        if (abs(value - 0.0) < mtk_tol) {
          cout << setw(12) << 0.0;
        } else {
          cout << setw(12) << value;
        }
      }
      cout << endl;
    }
    cout << endl;
  }

  // Collect garbage :)
  if (lambdas != nullptr) {
    delete[] lambdas;
    lambdas = nullptr;
  }
  if (alphas != nullptr) {
    delete[] alphas;
    alphas = nullptr;
  }
  return true;
}

// Stage 3: Final stage. Construct the output array with the operator and its
// weights:

bool AssembleOperator(void) {

  // The output array will have this form:
  // 1. The first entry of the array will contain the used order kk.
  // 2. The second entry of the array will contain the collection of
  // approximating coefficients for the interior of the grid.
  // 3. IF kk > 2, then the third entry will contain a collection of weights.
  // 4. IF kk > 2, the next dim_null entries will contain the collections of
  // approximating coefficients for the west boundary of the grid.

  if (kk > 2) {
    lDD = 1 + kk + kk + dim_null*num_bndy_approxs;
  } else {
    lDD = 1 + kk;
  }

  try {
    DD = new double[lDD];
    memset(DD, 0.0, sizeof(DD[0])*lDD);
  } catch (bad_alloc &memory_allocation_exception) {
    cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
    cerr << memory_allocation_exception.what() << endl;
  }

  // 1. The first entry of the array will contain the used order kk.
  DD[0] = kk;

  // 2. The second entry of the array will contain the collection of
  // approximating coefficients for the interior of the grid.
  for (auto ii = 0; ii < kk; ii++) {
    DD[ii + 1] = oo[ii];
  }

  // 3. IF kk > 2, then the third entry will contain a collection of weights.
  if (kk > 2) {
    for (auto ii = 0; ii < kk; ii++) {
      DD[(1 + kk) + ii] = qq[ii];
    }
  }

  // 4. IF kk > 2, the next dim_null entries will contain the collections of
  // approximating coefficients for the west boundary of the grid.
  if (kk > 2) {
    auto offset = (2*kk + 1);
    int mm {};
    for (auto ii = 0; ii < dim_null; ii++) {
      for (auto jj = 0; jj < num_bndy_approxs; jj++) {
        DD[offset + (mm)] = mim_bndy[jj*dim_null + ii];
        mm++;
      }
    }
  }

  return true;
}



// Main module. Implementation of the CRS Algorithm to compute a k-th order
// accurate mimetic 1D divergence operator.
// Algorithm: http://www.csrc.sdsu.edu/research_reports/CSRCR2013-02.pdf
