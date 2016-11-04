/*!
  \file mtk_2d_grad.h

  \brief Includes the definition of the class MTK_1DLap.

  This class implements a 1D Castillo-Grone (CG) LAPLACIAN matrix operator.

  \date: Sunday, September 02, 2012
  \version: 2012-09-02.
  \author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
 */
 /*
  Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
  State University (SDSU).

  http:www.csrc.sdsu.edu/mtk/

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
  esanchez@sciences.sdsu.edu

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
 */

#ifndef MTK_INCLUDE_MTK_2D_GRAD_H
#define MTK_INCLUDE_MTK_2D_GRAD_H


#include <iostream>
#include <iomanip>

#include "glpk.h"

#include "mtk_roots.h"
#include "mtk_tool_manager.h"
#include "mtk_dense_matrix.h"
#include "mtk_blas_facade.h"
#include "mtk_lapack_facade.h"
#include "mtk_glpk_facade.h"
#include "mtk_1d_grad.h"
#include "mtk_1d_div.h"

namespace mtk{
class MTK_2DGrad {
 public:
  friend MTK_DenseMatrix;
  /*!
  \brief
  */
  friend std::ostream& operator<< (std::ostream& stream,
                                   const MTK_2DGrad& matrix);
  /*!
  \brief
  */
  MTK_2DGrad(void);

  /*!
  \brief
  */
  MTK_2DGrad(int);

  /*!
  \brief
  */
  MTK_2DGrad(int, MTK_Real);
 
  /*!
  \brief
  */
  MTK_2DGrad* Construct2DGrad(int, MTK_Real);

  MTK_2DGrad* Construct2DGrad();


  /*!
  \brief
  */
  ~MTK_2DGrad();

  /*!
  \brief
  */
  MTK_DenseMatrix* ReturnAsMatrix(int , MTK_Real, MTK_Real, int , MTK_Real, MTK_Real);
  
  
  MTK_Real* ReturnAsReal(int, MTK_Real, MTK_Real, int, MTK_Real, MTK_Real);


  
 private:
 
  int order_; //< Selected order of accuracy.
  MTK_Real mimetic_tol_;  //< Mimetic tolerance.

  int kk_;                //< Order of numerical accuracy of the operator.
  
};
}
#endif
