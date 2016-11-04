#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>

#include "mtk_roots.h"
#include "mtk_1d_div.h"
#include "mtk_glpk_facade.h"


using namespace std;
using namespace mtk;

mtk::MTK_Real Optimizer (mtk::MTK_Real *A, int nrows, int ncols,\
		      int kk, mtk::MTK_Real *hh, mtk::MTK_Real *qq, \
		      int robjective, mtk::MTK_Real mimetic_threshold, int copy)
{
  glp_prob *lp;
  int problem_size, lp_nrows, lp_ncols; //Number of rows and columns
  int *ia, *ja;                 // of the LP problem`
  mtk::MTK_Real *ar, *objective, *rhs;
#if MTK_DEBUG_LEVEL>0
  mtk::MTK_Real  obj_value;
#endif
  int matsize;
  char rname[5], cname[5];
#if MTK_DEBUG_LEVEL >1
  char lp_file_name[18];
#endif
  int glp_index;
  int ii, jj;
  mtk::MTK_Real x1;
  lp_nrows = kk;                // LP problem size is the same as the order.
  lp_ncols = kk;
  mtk::MTK_Real *err;
  matsize = lp_nrows * lp_ncols;
  mtk::MTK_ToolManager calculate_norm;

//IMPORTANT: GLPK indexs from 1, so get the extra space needed.
// Memory allocation
  problem_size = lp_nrows * lp_ncols + 1;

  ia = (int *) malloc (sizeof (int) * (problem_size));

  if (ia == nullptr)
  {
    cout <<"Problem with malloc"<<endl;
    return 2;
  }
  ja = (int *) malloc (sizeof (int) * (problem_size));
  if (ja == nullptr)
  {
    cout <<"Problem with malloc"<<endl;
    return 3;
  }
  ar = (MTK_Real *) malloc (sizeof (MTK_Real) * (problem_size));
  if (ar == nullptr)
  {
    cout <<"Problem with malloc"<<endl;
    return 4;
  }
  objective = (MTK_Real *) malloc (sizeof (MTK_Real) * (lp_ncols+1));
  if (objective == nullptr)
  {
    cout <<"Problem with malloc"<<endl;
    return 5;
  }
  rhs = (MTK_Real *) malloc (sizeof (MTK_Real) * (lp_nrows + 1));
  if (rhs == nullptr)
  {
    cout <<"Problem with malloc"<<endl;
    return 6;
  }
  err = (MTK_Real *) malloc (sizeof (MTK_Real) * (lp_nrows ));
  if (err == nullptr)
  {
    cout <<"Problem with malloc"<<endl;
    return 7;
  }
#if MTK_DEBUG_LEVEL>0
  cout <<"Problem size = "<< problem_size<< " rows = "<<lp_nrows 
    << "cols = "<<lp_ncols << endl;
#endif
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
#if MTK_DEBUG_LEVEL>0
    cout <<"Using row "<< robjective+1<<" as objective."<< endl;
#endif
    for (jj = 0; jj < kk; jj++)
    {
      objective[glp_index] = A[jj + robjective * ncols];
      glp_index++;
    }
#if MTK_DEBUG_LEVEL >0
    cout << endl;
#endif
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
//    for(ii=0;ii<=lp_nrows;ii++)
//    {
//        cout<<rhs[ii]<<" "<<endl;
 //   }
 //   cout<<endl;
    
    // Setting up the objective function
    for (ii = 1; ii <= lp_ncols; ii++)
    {
      glp_set_obj_coef (lp, ii, objective[ii]);
    }


// Constrains
    for (ii = 1; ii <= lp_ncols; ii++)
    {
      glp_set_col_bnds (lp, ii, GLP_LO, mimetic_threshold, 0.0);
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
#if MTK_DEBUG_LEVEL >1
    sprintf(lp_file_name,"LP_row_%02d.lp",robjective+1);
    glp_write_lp(lp, nullptr,lp_file_name);
#endif
    /* solve problem */
    glp_simplex (lp, nullptr);
    
    /* Check status */
    
    if (glp_get_status(lp) == GLP_OPT)
    {    
    
    /* recover and display results */
    

      for(ii=1;ii<=lp_ncols;ii++)
      {
        err[ii-1] = qq[ii-1]- glp_get_col_prim(lp,ii);
      }
#if MTK_DEBUG_LEVEL>0
      obj_value = glp_get_obj_val (lp);

      cout <<"     CBS" << "\t"<<"CRS"<<endl;
      for (ii = 0; ii < lp_ncols; ii++)
      {
        cout <<"q"<<ii+1<<" : "<<glp_get_col_prim(lp,ii+1) <<"\t"<<qq[ii]<<endl;
      }
      cout <<"Objective function value (row"<<robjective+1<<") = " << obj_value << endl;
#endif
     if(copy)
     {
       for(ii=0;ii<lp_ncols;ii++)
       {
          qq[ii]= glp_get_col_prim(lp,ii+1);  // collect the solutions
   
       }
       for(ii=0;ii<kk/2-1;ii++)
       {
          qq[ii+kk+1]=hh[ii+kk+1];  // collect the lambdas
       }
     }
     x1=calculate_norm.CalculateNorm(err,lp_ncols);
  
    }
    else
    {
       x1=1e6; 
       
    }
    
    /* housekeeping */
    glp_delete_prob (lp);
    glp_free_env ();
    free (ia);
    free (ja);
    free (ar);
    free (objective);
    free (rhs);
    free (err);
    return x1;
}


