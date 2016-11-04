/*!
  \file mtk_system.h

  \brief Includes the definition of the class MTK_System.

  This file contains the definition of a class implementing a system of
  equations, which quite often will yield the solution to a problem of interest.

  \date: Saturday, April 14 2012 08:06 PM
  \version: 2011-14-08-01.
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

#ifndef MTK_SYSTEM_H
#define MTK_SYSTEM_H

//class MTK_1DNode;
class MTK_1DArrayDenseStencilMatrix;
class MTK_LR1DUniformSourceGrid;
namespace mtk
{
/*! \brief Definition of the class MTK_System.

  This file contains the definition of a class implementing a system of
  equations, which quite often will yield the solution to a problem of interest.
  \todo Monday, September 03, 2012: Add the option on how to solve.
*/
class MTK_System {

  public:
    /*! Default constructor.  */
    /*! This procedure CONSTRUCTS a default instance. */
    MTK_System(void):
      AA_(0x0),
      bb_(0x0) {}

    /*! Constructs a system out of a stencil matrix and a right-hand side. */
    /*! Constructs a system out of a stencil matrix and a right-hand side.
    */
    MTK_System(MTK_1DArrayDenseStencilMatrix *AA, MTK_LR1DUniformSourceGrid *bb);

    /*! Default destructor.  */
    /*! This procedure DESTRUCTS an already created instance. */
    ~MTK_System() {

    }

    /*! Prints the numerical information about an instance.  */
    /*! Prints the numerical information about an instance. */
    void Print(void);

    /*! Solves the system.  */
    /*! Solves the system. */
    void Solve(void);

    /*! Gets the attained solution.  */
    /*! Gets the attained solution. */
    MTK_1DNode** Solution(void);

    /*! Get the \f$ i \f$-th element of the solution. */
    /*! Get the \f$ i \f$-th element of the solution. */
    double GetValueFromSolution(int ii);

    /*! Gets the RHS.  */
    /*! Gets the RHS. */
    MTK_LR1DUniformSourceGrid* rhs(void) {return bb_;}

  private:
    /*! Our own absolute value function.  */
    /*! Our own absolute value function. */
    double absval(double val) {
      if (val >= 0)
        return val;
      else
        return -val;
    }

    int representation_;            /*!< Is this a dense or a sparse system? */
    MTK_1DArrayDenseStencilMatrix *AA_;        /*!< Matrix of coefficients. */
    MTK_LR1DUniformSourceGrid *bb_;  /*!< Right-hand side. */
    double *solution_;              /*!< Attained solution. */
};
}
#endif
