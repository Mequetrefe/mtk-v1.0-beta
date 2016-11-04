/*!
\file mtk_steady_1d_problem.h

\brief Defines a steady 1D problem.

This file contains the definition of a class which allows for the
representation of a steady-state 1D problem in any domain of interest.

\date: Wednesday, April 18 2012 09:08 PM
\version: 2012-18-04-01.
\author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
*/
/*
Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
State University (SDSU).

http://www.csrc.sdsu.edu/mtk/

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

#ifndef MTK_INCLUDE_MTK_STEADY_1D_PROBLEM_H
#define MTK_INCLUDE_MTK_STEADY_1D_PROBLEM_H

/*! \brief Steady 1D problem.
 *
 * This class implements a steady-state problem. It contents the collection of
 * parameters that, as today, describe a problem. We are interested in modeling
 * a Robin's Boundary Value Problem since this represents the most general
 * scenario. The modeled Robin's problem has the general form:
 * \f[
 * \alpha u(a) + \beta u^\prime(a) = \textrm{\texttt{west\_bndy\_value}}
 * \f]
 * and
 * \f[
 * \gamma u(b) + \delta u^\prime(b) = \textrm{\texttt{east\_bndy\_value}}
 * \f]
 * where \f$ \Omega=[a,b]\f$ and the coefficients are given, along with the
 * values at the boundary by means of the constructor member function.
 */
class MTK_Steady1DProblem {

  public:
    /*! Default constructor.  */
    /*! This functions CONSTRUCTS a default problem. */
    MTK_Steady1DProblem(void);

    /*! Constructs a steady 1D problem from the collection of values. */
    /*! This procedure CONSTRUCTS an instance using the most commonly used
     * parameters given the context.
     */
    MTK_Steady1DProblem(double alpha, double beta, double west_bndy_value,
                       double gamma, double delta,double east_bndy_value);

    /*! Prints a steady 1D problem. */
    /*! Prints a steady 1D problem.
     */
    void MTKPrint(void);

  private:
    double alpha_;            /*!< Coefficient of the function at west. */
    double beta_;             /*!< Coefficient of the derivative at west. */
    double west_bndy_value_;  /*!< Value at the west boundary. */
    double gamma_;            /*!< Coefficient of the function at east. */
    double delta_;            /*!< Coefficient of the derivative at east. */
    double east_bndy_value_;  /*!< Value at the east boundary. */
};

#endif
