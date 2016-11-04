/*!
 * \file mtk_tool_manager.h
 *
 * \brief Class intended to perform several general routines.
 *
 * This file contains the definition of a class intended to perform some general
 * tasks along the entire algorithmical layout of the library.
 *
 * \date: Thursday, October 27 2011, 09:19 PM
 * \version: 2011-10-27-01.
 *
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

#ifndef MTK_INCLUDE_MTK_TOOL_MANAGER_H
#define MTK_INCLUDE_MTK_TOOL_MANAGER_H

#include "mtk_roots.h"

class MTK_1DNode;

namespace mtk{

/*! \brief Tool manager.
 *
 * This class implements all the required functionalities that help to ensure
 * correct execution of the includes routines. This class performs some general
 * tasks along the entire algorithmic layout of the library.
 */
class MTK_ToolManager {

  public:
    /*! Default constructor.  */
    /*! This functions CONSTRUCTS a default Tool Manager. */
    MTK_ToolManager();

    /*! Used to abort execution.  */
    /*! Used to abort execution. */
    void MTKAbort(char *error, int line, char *file);

    /*! Used to compute the 2-norm of the difference between two vectors. */
    /*! Used to compute the 2-norm of the difference between two vectors. */
    MTK_Real RelativeNorm2Difference(MTK_1DNode **xx, MTK_1DNode **yy, int NN);

    /*! Used to compute the 2-norm of the difference between two vectors. */
    /*! Used to compute the 2-norm of the difference between two vectors. */
    MTK_Real MagnitudeDifferenceWest(MTK_1DNode **xx, MTK_1DNode **yy, int NN);

    /*! Used to compute the 2-norm of the difference between two vectors. */
    /*! Used to compute the 2-norm of the difference between two vectors. */
    MTK_Real MagnitudeDifferenceEast(MTK_1DNode **xx, MTK_1DNode **yy, int NN);


    MTK_Real CalculateNorm(MTK_Real *A,int n);
  private:
    int tmp;  /*!< Temporal integer value??? */
};
}
#endif
