/*!
\file mtk_dense_matrix.cc

\brief Implements a common dense matrix, using a 1-D array.

For developing purposes, it is better to have a not-so-intrincated data
structure implementing matrices. This is the purpose of this class: to be used
for prototypes of new code for small test cases. In every other instance, this
should be replaced by the most appropriate sparse matrix.

\date: Monday, September 03, 2012

\version: 2012-09-03.

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

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iostream>
#include <iomanip>

#include "mtk_roots.h"
#include "mtk_dense_matrix.h"

namespace mtk {

std::ostream& operator <<(std::ostream &stream, mtk::MTK_DenseMatrix &in) {

  for (auto ii = 0; ii < in.matrix_data_->num_rows(); ii++) {
    for (auto jj = 0; jj < in.matrix_data_->num_cols(); jj++) {
      stream << setw(12) << in.values_[ii][jj];
    }
  stream << endl;
}
stream << endl;

return stream;
}

}

mtk::MTK_DenseMatrix::MTK_DenseMatrix(): values_(nullptr), continuum_(nullptr) {

  matrix_data_ = new MTK_Matrix();
}

mtk::MTK_DenseMatrix::MTK_DenseMatrix(const  mtk::MTK_DenseMatrix &in) {}

mtk::MTK_DenseMatrix::MTK_DenseMatrix(
    const int &num_rows,
    const int &num_cols) {

  int ii;

  if (num_rows <= 0) {
    cerr << "Incorrect parameter at line " <<  __LINE__ - 1 << endl;
  }
  if (num_cols <= 0) {
    cerr << "Incorrect parameter at line " <<  __LINE__ - 1 << endl;
  }

  matrix_data_ = new MTK_Matrix();

  matrix_data_->num_rows_ = num_rows;
  matrix_data_->num_cols_ = num_cols;
  matrix_data_->num_values_ = num_rows*num_cols;

  continuum_ = (MTK_Real *) malloc(num_rows*num_cols*sizeof(MTK_Real));

  values_ = (MTK_Real **) malloc(num_rows*sizeof(MTK_Real));

  for (ii = 0; ii < num_rows; ii++) {
    values_[ii] = continuum_ + (ii*num_cols);
  }
}

mtk::MTK_DenseMatrix::MTK_DenseMatrix(const int &rank, const bool &padded) {

  int aux{};

  matrix_data_ = new MTK_Matrix();

  if (padded) {
    aux = 1;
  }

  matrix_data_->num_rows_ = aux + rank + aux;
  matrix_data_->num_cols_ = rank;
  matrix_data_->num_values_ = matrix_data_->num_rows_*matrix_data_->num_cols_;

  continuum_ = (MTK_Real *) malloc(matrix_data_->num_values_*sizeof(MTK_Real));

  values_ = (MTK_Real **) malloc(matrix_data_->num_rows_*sizeof(MTK_Real));

  for (auto ii = 0; ii < matrix_data_->num_rows_; ii++) {
    values_[ii] = continuum_ + (ii*matrix_data_->num_cols_);
  }
  for (auto ii =0 ; ii < matrix_data_->num_rows_; ++ii) {
    for (auto jj = 0; jj < matrix_data_->num_cols_; ++jj) {
      values_[ii][jj] = (ii == jj+aux)? 1.0: 0.0;
    }
  }
}

// DEF. In linear algebra, a VANDERMONDE MATRIX is a matrix with the terms of a
// geometric progression in each row. This progression uses the terms of a
// given GENERATOR VECTOR.

// Generates a Vandermonde matrix with a given generator vector gen:
mtk::MTK_DenseMatrix::MTK_DenseMatrix(const MTK_Real *gen,
                                      const int gen_length,
                                      const int pro_length,
                                      const bool transpose) {

  matrix_data_ = new MTK_Matrix();

  auto num_rows = gen_length;
  auto num_cols = pro_length;

  continuum_ = (MTK_Real *) malloc(num_rows*num_cols*sizeof(MTK_Real));

  values_ = (MTK_Real **) malloc(num_rows*sizeof(MTK_Real));

  for (auto ii = 0; ii < num_rows; ii++) {
    values_[ii] = continuum_ + (ii*num_cols);
  }

//   try {
//     TT = new MTK_Real[gen_length*pro_length];
//     memset(TT, 0.0, sizeof(TT[0])*gen_length*pro_length);
//   } catch (bad_alloc &memory_allocation_exception) {
//     cerr << "Memory allocation exception on line " << __LINE__ - 3 << endl;
//     cerr << memory_allocation_exception.what() << endl;
//   }

  if (!transpose) {
    for (auto ii = 0; ii < gen_length; ii++) {
      for (auto jj = 0; jj < pro_length; jj++) {
        (*values_)[idx(ii,pro_length,jj)] = pow(gen[ii], (MTK_Real) jj);
      }
    }
  } else {
    for (auto ii = 0; ii < pro_length; ii++) {
      for (auto jj = 0; jj < gen_length; jj++) {
        (*values_)[idx(ii,gen_length,jj)] = pow(gen[jj], (MTK_Real) ii);
      }
    }
  }

}

mtk::MTK_Real mtk::MTK_DenseMatrix::GetValue(
    const int &rr,
    const int &cc) const {

  return values_[rr][cc];
}

void  mtk::MTK_DenseMatrix::SetValue(
    const int &rr,
    const int &cc,
    const mtk::MTK_Real &val) {

  values_[rr][cc] = val;
}

mtk::MTK_DenseMatrix::~MTK_DenseMatrix() {

}

mtk::MTK_DenseMatrix mtk::MTK_DenseMatrix::Kron(
    const mtk::MTK_DenseMatrix &aa,
    const mtk::MTK_DenseMatrix &bb) {

  mtk::MTK_DenseMatrix output;
  register MTK_Real aa_factor{};  // Kept in register (KiR).
  register int row_offset{};    // Row-offset to access output matrix (KiR).
  register int col_offset{};    // Col-offset to access output matrix (KiR).

  output.matrix_data_->num_rows_ = aa.matrix_data_->num_rows_;
  output.matrix_data_->num_rows_ *= bb.matrix_data_->num_rows_;
  output.matrix_data_->num_cols_ = aa.matrix_data_->num_cols_;
  output.matrix_data_->num_cols_ *= bb.matrix_data_->num_cols_;

  register int kk_num_cols{output.matrix_data_->num_cols_};

  auto mm = aa.matrix_data_->num_rows_;
  auto nn = aa.matrix_data_->num_cols_;
  auto pp = bb.matrix_data_->num_rows_;
  auto qq = bb.matrix_data_->num_cols_;

  for (auto ii = 0; ii < mm; ++ii) {
    row_offset = ii*pp;
    for (auto jj = 0; jj < nn; ++jj) {
      col_offset = jj*qq;
      aa_factor = (*aa.values_)[ii*nn + jj];
      for (auto ll = 0; ll < pp; ++ll) {
        for (auto oo = 0; oo < qq; ++oo) {
          auto index = (ll + row_offset)*kk_num_cols + (oo + col_offset);
          (*output.values_)[index] = aa_factor*((*bb.values_)[ll*qq + oo]);
        }      }
    }
  }
  return output;
}

mtk::MTK_DenseMatrix mtk::MTK_DenseMatrix::Transpose() const {


//  mtk::MTK_DenseMatrix output;

  auto rr = matrix_data_->num_rows_;
  auto cc = matrix_data_->num_cols_;

  MTK_DenseMatrix output(rr,cc);

  for (auto ii = 0; ii < rr; ii++) {
    for (auto jj = 0; jj < cc; jj++) {
      (*output.values_)[idx(jj,rr,ii)] = (*values_)[idx(ii,cc,jj)];
    }
  }

  return output;
}

void mtk::MTK_DenseMatrix::Transpose( 
    const mtk::MTK_DenseMatrix &aa,
     mtk::MTK_DenseMatrix &bb) 
{


//  mtk::MTK_DenseMatrix output;

  auto rr = aa.matrix_data_->num_rows_;
  auto cc = aa.matrix_data_->num_cols_;

//  MTK_DenseMatrix output(rr,cc);

  for (auto ii = 0; ii < rr; ii++) {
    for (auto jj = 0; jj < cc; jj++) {
      bb.SetValue(jj,ii,aa.GetValue(ii,jj));
    }
  }

//  return output;
}



void mtk::MTK_DenseMatrix::Kronecker(
    const mtk::MTK_DenseMatrix &aa,
    const mtk::MTK_DenseMatrix &bb,
    mtk::MTK_DenseMatrix &cc) 
{

  auto rra = aa.matrix_data_->num_rows_;
  auto cca = aa.matrix_data_->num_cols_;
  auto rrb = bb.matrix_data_->num_rows_;
  auto ccb = bb.matrix_data_->num_cols_;

    
  for (auto ii=0;ii<rra;ii++)
  {
      for(auto jj=0;jj<cca;jj++)
      {
          for(auto kk=0;kk<rrb;kk++)
          {
              for(auto ll=0;ll<ccb;ll++)
              {
                  cc.SetValue(kk+ii*rrb,ll+jj*ccb, aa.GetValue(ii,jj)*bb.GetValue(kk,ll));
              }
          }
      }
  }
  
}




mtk::MTK_Real *mtk::MTK_DenseMatrix::DenseMatrix_To_MTK_Real (const mtk::MTK_DenseMatrix &aa)
{
 return(aa.continuum_);   
}



mtk::MTK_DenseMatrix* mtk::MTK_DenseMatrix::MatrixMultiplication(
    const mtk::MTK_DenseMatrix &aa,
    const mtk::MTK_DenseMatrix &bb)
{
    
  auto rr1 = aa.matrix_data_->num_rows_;
  auto cc1 = aa.matrix_data_->num_cols_;
  auto rr2 = bb.matrix_data_->num_rows_;
  auto cc2 = bb.matrix_data_->num_cols_;

  if(cc1 != rr2 )
  {
      cout << "Problem with Matrix dimensions..." <<endl;
 //     &cc= nullptr;
  }
  
  
  auto transa='N';
  auto alpha=1.0;
  auto beta=0.0;
  
  mtk::MTK_DenseMatrix *aaT, *bbT, *ccT, *cc;
  aaT=new mtk::MTK_DenseMatrix(cc1,rr1);
  bbT=new mtk::MTK_DenseMatrix(cc2,rr2);
  ccT=new mtk::MTK_DenseMatrix(cc2,rr1);
  cc=new mtk::MTK_DenseMatrix(rr1,cc2);
  
  cc->Transpose(aa,*aaT);
  cc->Transpose(bb,*bbT);
#if MTK_DEBUG_LEVEL >1
  std::cout << "mtk_dense_matrix.cc Matrix A Transposed "<<std::endl;
  
  std::cout <<*aaT << " " ;


  std::cout << "mtk_dense_matrix.cc Matrix B Transposed "<<std::endl;
  
  std::cout <<*bbT << " " ;
#endif

      mtk::dgemm_(&transa, &transa,&rr1,&cc2,&cc1,&alpha, aaT->DenseMatrix_To_MTK_Real(*aaT),&rr1, bbT->DenseMatrix_To_MTK_Real(*bbT),&rr2,&beta,ccT->DenseMatrix_To_MTK_Real(*ccT),&rr1);

  ccT->Transpose(*ccT,*cc);

  free (aaT);
  free (bbT);
  free (ccT);
  return cc;
    
}




