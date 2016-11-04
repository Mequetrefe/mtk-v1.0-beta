/*!
  \file mtk_lr_1d_uniform_known_solution_holder_grid.h

  \brief Includes the definition of MTK_LR1DUniformKnownSolutionHolderGrid

  This file contains the definition of the class
  MTK_LR1DUniformKnownSolutionHolderGrid, which is a class intended for
  testing the numerical precision of the method of study as it is compared with
  a different method or known solution.

  \date: Monday, September 03, 2012
  \version: 2011-09-03.
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

#ifndef MTK_INCLUDE_MTK_LR_1D_UNIFORM_KNOWN_SOLUTION_HOLDER_GRID_DOUBLE_H
#define MTK_INCLUDE_MTK_LR_1D_UNIFORM_KNOWN_SOLUTION_HOLDER_GRID_DOUBLE_H

#include "mtk_enums.h"
#include "mtk_roots.h"
namespace mtk{

class MTK_1DNode;
class MTK_LR1DUniformStaggeredGrid;

/*! \brief Class to compare methods and attained numerical solution.

  This file contains the definition of the class
  MTK_LR1DUniformKnownSolutionHolderGrid, which is a class intended for
  testing the numerical precision of the method of study as it is compared with
  a different method or known solution.
  \todo Monday, September 03, 2012: Move SolveInitialValue to a solver class.
 */
}
#endif
