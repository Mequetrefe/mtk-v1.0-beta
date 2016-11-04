/*!
\file mtk_matrix.cc

\brief implementing the representation of a matrix in the MTK.

Definition of the representation for the matrices implemented in the MTK.

\date: June 20, 2014, 11:37 AM

\version: 2015-02-23

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at sciences dot sdsu dot edu

\bug No known bugs.
*/

/*
Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
State University (SDSU).

Website for the project: http://www.csrc.sdsu.edu/mtk/

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

esanchez at sciences dot sdsu dot edu

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "mtk_roots.h"
#include "mtk_matrix.h"

using namespace mtk;
using namespace std;

MTK_Matrix::MTK_Matrix():
  num_rows_(),
  num_cols_(),
  num_values_(),
  num_non_zero_(),
  num_zero_(),
  num_null_(),
  unitary_offset_(),
  rank_(),
  num_unknowns_(),
  bandwith_(),
  abs_density_(),
  rel_density_(),
  abs_sparsity_(),
  rel_sparsity_() {}

MTK_Matrix::MTK_Matrix(const MTK_Matrix &in):
  num_rows_(in.num_rows_),
  num_cols_(in.num_cols_),
  num_values_(in.num_values_),
  num_non_zero_(in.num_non_zero_),
  num_zero_(in.num_zero_),
  num_null_(in.num_null_),
  unitary_offset_(in.unitary_offset_),
  rank_(in.rank_),
  num_unknowns_(in.num_unknowns_),
  bandwith_(in.bandwith_),
  abs_density_(in.abs_density_),
  rel_density_(in.rel_density_),
  abs_sparsity_(in.abs_sparsity_),
  rel_sparsity_(in.rel_sparsity_) {}

MTK_Matrix::~MTK_Matrix() {}

int MTK_Matrix::num_rows() const {
  
  return num_rows_;
}

int MTK_Matrix::num_cols() const {
  
  return num_cols_;
}

int MTK_Matrix::num_values() const {
  
  return num_values_;
}

int MTK_Matrix::num_non_zero() const {
  
  return num_non_zero_;
}

int MTK_Matrix::num_zero() const {
  
  return num_zero_;
}

int MTK_Matrix::num_null() const {
  
  return num_null_;
}

int MTK_Matrix::unitary_offset() const {
  
  return unitary_offset_;
}

int MTK_Matrix::rank() const {
  
  return rank_;
}

int MTK_Matrix::num_unknowns() const {
  
  return num_unknowns_;
}

int MTK_Matrix::bandwith() const {
  
  return bandwith_;
}

MTK_Real MTK_Matrix::abs_density() const {
  
  return abs_density_;
}

MTK_Real MTK_Matrix::rel_density() const {
  
  return rel_density_;
}

MTK_Real MTK_Matrix::abs_sparsity() const {
  
  return abs_sparsity_;
}

MTK_Real MTK_Matrix::rel_sparsity() const {
  
  return rel_sparsity_;
}

void MTK_Matrix::set_unitary_offset(int unitary_offset) {
  
  unitary_offset_ = unitary_offset;
}

void MTK_Matrix::UpdateMatrixState() {
  
  
}
