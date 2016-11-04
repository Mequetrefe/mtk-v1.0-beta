/*!
 * \file mtk_system.cc
 *
 * \brief Constructs a default MTK_System.
 *
 * Constructs a default MTK_System. Note how the default order is
 * at least 2.
 *
 * \date: Monday, September 03, 2012
 * \version: 2012-09-03.
 * \author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
 */
 /* Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
 * State University (SDSU).
 *
 * http:www.csrc.sdsu.edu/mtk/
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification are permitted provided that the following conditions are met:
 *
 * -# Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * -# Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * -# Neither the name of the CSRC, SDSU nor the names of its contributors may
 * be used to endorse or promote products derived from this software without
 * specific prior written permission.
 * -# Modifications whether they are partial or complete; whether they are
 * additions or eliminations should be reported through an email at:
 * esanchez@sciences.sdsu.edu
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <iostream>
#include <iomanip>

#include "mtk_1d_node.h"
#include "mtk_1d_array_dense_stencil_matrix.h"
#include "mtk_lr_1d_uniform_source_grid.h"
#include "mtk_system.h"

using namespace std;

#define ZERO 1.0E-20
#define true 1
#define false 0

/*! Constructs a system out of a stencil matrix and a right-hand side. */
/*! Constructs a system out of a stencil matrix and a right-hand side.
*/
MTK_System::MTK_System(MTK_1DArrayDenseStencilMatrix *AA, MTK_LR1DUniformSourceGrid *bb) {

  representation_ = AA->representation();
  AA_ = AA;
  bb_ = bb;
}

/*! Solves a system. */
/*! Solves a system.
 */
void MTK_System::Solve(void) {

  double *matrix;
  double *rhs;
  double **A;
  double *X;

  int matrix_num_rows;
  int matrix_num_cols;
  int size_rhs;
  int ii;
  int jj;

  matrix = AA_->dense_values();
  rhs = bb_->nodes();
  matrix_num_rows = AA_->num_rows();
  matrix_num_cols = matrix_num_rows;
  size_rhs = matrix_num_rows;

  // BEGIN GAUSS ELIMINATION!
  double C,XM,SUM;
  int N,M,I,J,ICHG,NN,IP,JJ,K,KK,OK;

  // 6. Generate augmented matrix:
  A = (double **) malloc(matrix_num_rows * sizeof(double*));
  for (ii = 0; ii < matrix_num_rows; ii++) {
    A[ii] = (double *) malloc((matrix_num_cols + 1)*sizeof(double));
  }

  // 7. Fill AA:
  for (ii = 0; ii < matrix_num_rows; ii++) {
    for (jj = 0; jj < matrix_num_cols; jj++) {
      A[ii][jj] = matrix[ii*matrix_num_cols + jj];
    }
    A[ii][jj] = rhs[ii];
  }

#if ( Debuglevel >= 1 )
  cout << "Augmented matrix of the system to solve:" << endl;
  // 8. Print augmented matrix:
  for (ii = 0; ii < matrix_num_rows; ii++) {
    for (jj = 0; jj < matrix_num_cols + 1; jj++) {
      cout << setw(10) << A[ii][jj] << " ";
    }
    cout << endl;
  }
#endif

  // 9. Allocate space for the solution:
  X = (double *) calloc(size_rhs, sizeof(double));

  // 10. Solve:
  OK = 1;
  N = matrix_num_rows;

  if (OK) {
    /* STEP 1 */
    NN = N - 1;
    M = N + 1;
    ICHG = 0;
    I = 1;
    while ((OK) && (I <= NN)) {
      /* STEP 2 */
      /* use IP instead of p */
      IP = I;
      while ((absval(A[IP-1][I-1]) <= ZERO) && (IP <= N))
        IP++;
      if (IP == M)
        OK = false;
      else {
        /* STEP 3 */
        if (IP != I) {
          for (JJ=1; JJ<=M; JJ++) {
            C = A[I-1][JJ-1];
            A[I-1][JJ-1] = A[IP-1][JJ-1];
            A[IP-1][JJ-1] = C;
          }
          ICHG = ICHG + 1;
        }
        /* STEP 4 */
        JJ = I + 1;
        for (J=JJ; J<=N; J++) {
          /* STEP 5 */
          /* use XM in place of m(J,I) */
          XM = A[J-1][I-1] / A[I-1][I-1];
          /* STEP 6 */
          for (K=JJ; K<=M; K++)
            A[J-1][K-1] = A[J-1][K-1] - XM * A[I-1][K-1];
          /*  Multiplier XM could be saved in A[J,I].  */
          A[J-1][I-1] = 0.0;
        }
      }
      I = I + 1;
    }
    if (OK) {
      /* STEP 7 */
      if (absval(A[N-1][N-1]) <= ZERO)
        OK = false;
      else {
        /* STEP 8 */
        /* start backward substitution */
        X[N-1] = A[N-1][M-1] / A[N-1][N-1];
        /* STEP 9 */
        for (K=1; K<=NN; K++) {
          I = NN - K + 1;
          JJ = I + 1;
          SUM = 0.0;
          for (KK=JJ ; KK<=N; KK++)
            SUM = SUM - A[I-1][KK-1] * X[KK-1];
          X[I-1] = (A[I-1][M-1]+SUM) / A[I-1][I-1];
        }
        /* STEP 10 */
        /* procedure completed successfully */
#if ( Debuglevel >= 1 )
        fprintf(stdout, "\n\nSolution vector:\n");
        for (I=1; I<=N; I++)
          fprintf(stdout, "  %12.8f", X[I-1]);
        fprintf (stdout, "\n\nwith %d row interchange(s)\n", ICHG);
#endif
        solution_ = X;
      }
    }
    if (!OK)
      printf("System has no unique solution\n");

    free(matrix);
    free(rhs);
    for (ii = 0; ii < matrix_num_cols; ii++) {
      free(A[ii]);
    }
    free(A);
  }
}

/*! Prints a system. */
/*! Prints a system.
 */
void MTK_System::Print(void) {

  int ii;
  int jj;
  int nn;

  AA_->Print(representation_);
  bb_->Print1DArray();
  cout << "Arising system of equation:" << endl;
  nn = AA_->size() + 2;
  cout.precision(2);
  for (ii = 0; ii < nn; ii++) {
    cout << "\t|";
    for (jj = 0; jj < nn; jj++) {
      cout << setw(6) << AA_->GetValue(ii, jj) << " ";
    }
    if (ii == floor((double) nn/2.0)) {
      cout << "\t| x =\t| " << fixed << bb_->GetValue(ii) << "\t|" << endl;
    } else {
      cout << "\t|\t| " << fixed << bb_->GetValue(ii) << "\t|" << endl;
    }
  }
  cout << endl;
}

/*! Get the \f$ i \f$-th element of the solution. */
/*! Get the \f$ i \f$-th element of the solution. */
double MTK_System::GetValueFromSolution(int ii) {

  if ((ii < 0) || (ii > (AA_->size() + 2 - 1))) {
    cout << "ERROR: Value out of range!" << endl;
    return -1000000000.0;
  } else {
    return solution_[ii];
  }
}

/*! Gets the attained solution.  */
/*! Gets the attained solution. */
MTK_1DNode** MTK_System::Solution(void) {

  int ii;
  int nn;
  MTK_1DNode **solution;

  nn = AA_->size();
  /* Allocate outgoing array of pointers to Nodes: */
  solution = (MTK_1DNode**) malloc((nn + 2)*sizeof(MTK_1DNode));
  if (solution == (MTK_1DNode**) NULL) {
    printf("ERROR 01: Nodal grid 1-D array could not be allocated!\n");
    printf("Exiting...\n");
  }
  /* Fill in the nodes values: */
  for (ii = 0; ii < (nn + 2); ii++) {
    solution[ii] = new MTK_1DNode(solution_[ii]);
  }
  return solution;
}
