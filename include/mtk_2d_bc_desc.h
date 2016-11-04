/*!
  \file mtk_2d_bc_desc.h


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

#ifndef MTK_INCLUDE_MTK_2D_BC_DESC_H
#define MTK_INCLUDE_MTK_2D_BC_DESC_H
#include "mtk_2d_uniform_cartesian_staggered_grid.h"
#include "mtk_2d_lap.h"
#include "mtk_roots.h"
namespace mtk{
    class MTK_2DBCDesc {

    public:
    
    MTK_2DBCDesc(void);
    
    MTK_2DBCDesc(

                 MTK_2DLap* lap,
                 MTK_2DUniformStaggeredGrid *GRID,
                 MTK_Real (*South)(MTK_Real xx, MTK_Real yy),
                 MTK_Real (*North)(MTK_Real xx, MTK_Real yy),
                 MTK_Real (*West)(MTK_Real xx, MTK_Real yy),
                 MTK_Real (*East)(MTK_Real xx, MTK_Real yy) 
    );
    
    ~MTK_2DBCDesc()
    {
        
    }
    
// MTK_2DBCDirichlet takes a 2D dense matrix Laplacian
// and fill the boundary conditions required to get Dirichlet conditions.

    void MTK_2DBCDirichlet(MTK_DenseMatrix *);

    void MTK_2DBoundaryRHS(MTK_Real *RHS);

    private:
      int   number_of_cells_X_;
      int   number_of_cells_Y_;
      MTK_Real East_;
      MTK_Real West_;
      MTK_Real South_;
      MTK_Real North_;
      MTK_2DUniformStaggeredGrid *GRID_;
      
      MTK_2DLap* lap_;
      MTK_Real (*fS_)(MTK_Real xx, MTK_Real yy);
      MTK_Real (*fN_)(MTK_Real xx, MTK_Real yy);
      MTK_Real (*fW_)(MTK_Real xx, MTK_Real yy);
      MTK_Real (*fE_)(MTK_Real xx, MTK_Real yy); 

                 
    };
}
#endif

    
    
    