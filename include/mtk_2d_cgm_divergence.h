/*!
\file mtk_2d_cgm_divergence.h

\brief Includes the definition of the class MTK_2DCGMDivergence.

This class implements a 2D DIVERGENCE matrix operator, based on the
Castillo-Grone Method (CGM).

\date: Thursday, September 12
\version: 2013-09-12.
\author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
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
esanchez@sciences.sdsu.edu

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef MTK_INCLUDE_MTK_2D_CGM_DIVERGENCE_H
#define MTK_INCLUDE_MTK_2D_CGM_DIVERGENCE_H

/*! \brief Compatible class MTK_2DCGMDivergence.
 *
 * This class implements a 2D Castillo-Grone (CG) Divergence matrix operator.
 */
class MTK_2DCGMDivergence {

  public:
    /*! Default constructor.  */
    /*! This functions CONSTRUCTS a default Gradient operator. */
    MTK_2DCGMDivergence(void);

    /*! Constructor #2. */
    /*! Constructs an operator starting from the usual parameters. */
    MTK_2DCGMDivergence(int order, int Nx, int Ny, double hx, double hy);

    /*! Prints an operator accordingly. */
    /*! Prints an operator accordingly.
     */
    void Print(int display);

    int order() {return order_;};

    /*! Default destructor.  */
    /*! This functions DESTRUCTS an already created Gradient operator. */
    ~MTK_2DCGMDivergence() {
      //delete [] dense_values_;
    }

    double *dense_values_;  /*! < Collection of values for dense operators. */

  private:
    int representation_;    /*! < Selected storage option. */
    int order_;             /*! < Selected order of accuracy.*/
    int size_;              /*! < Number of nodes I will be acting upon. */
    int rows_;
    int cols_;
};

#endif
