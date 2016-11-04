/*!
\file mtk_glpk_facade.h

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

#ifndef MTK_INCLUDE_MTK_GLPK_FACADE_H
#define MTK_INCLUDE_MTK_GLPK_FACADE_H

#include <iostream>
#include <iomanip>

#include "glpk.h"

#include "mtk_roots.h"
#include "mtk_tool_manager.h"
#include "mtk_dense_matrix.h"
#include "mtk_blas_facade.h"
#include "mtk_lapack_facade.h"
#include "mtk_1d_div.h"

	mtk::MTK_Real Optimizer (mtk::MTK_Real *A,
			int nrows,
			int ncols,
			int kk,
			mtk::MTK_Real *hh,
			mtk::MTK_Real *qq,
			int robjective,
			mtk::MTK_Real mimetic_tol,
			int copy);
#endif