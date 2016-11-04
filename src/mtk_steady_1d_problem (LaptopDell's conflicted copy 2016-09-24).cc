/*!
\file mtk_steady_1d_problem.cc

\brief Creates a default instance of the MTK_Steady1DProblem class.

This file contains the definition of a constructor procedure for the
MTK_Steady1DProblem class. Specifically, it constructs the default instance.

\date: Monday, September 03, 2012
\version: 2012-09-03.
\author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
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

#include <iostream>

#include "mtk_steady_1d_problem.h"

using namespace std;
using namespace mtk;
/*! Constructs a default MTK_Steady1DProblem. */

/*! Constructs a default MTK_Steady1DProblem.
 */
MTK_Steady1DProblem::MTK_Steady1DProblem(void):
  alpha_(0.0),
  beta_(0.0),
  west_bndy_value_(0.0),
  gamma_(0.0),
  delta_(0.0),
  east_bndy_value_(0.0) {

}

/*! Constructs a MTK_Steady1DProblem object. */

/*! Constructs a MTK_Steady1DProblem object of the form
 * \f[
 * \alpha u(a) + \beta u^\prime(a) = \textrm{\texttt{west\_bndy\_value}}
 * \f]
 * \f[
 * \gamma u(b) + \delta u^\prime(b) = \textrm{\texttt{east\_bndy\_value}}
 * \f]
 * in \f$ \Omega=[a,b]\f$ from a collection of defined values.
 * \param alpha INPUT: \f$ \alpha \f$ coefficient for the function at west,
 * \param beta INPUT: \f$ \beta \f$ coefficient for the derivative at west,
 * \param west_bndy_value INPUT: Boundary condition at west.
 * \param gamma INPUT: \f$ \gamma \f$ coefficient for the function at east,
 * \param delta INPUT: \f$ \delta \f$ coefficient for the derivative at east,
 * \param east_bndy_value INPUT: Boundary condition at east.
 */
MTK_Steady1DProblem::MTK_Steady1DProblem(MTK_Real alpha, MTK_Real beta,
                                       MTK_Real west_bndy_value,
                                       MTK_Real gamma, MTK_Real delta,
                                       MTK_Real east_bndy_value) {
  this->alpha_ = alpha;
  this->beta_ = beta;
  this->west_bndy_value_ = west_bndy_value;
  this->gamma_ = gamma;
  this->delta_ = delta;
  this->east_bndy_value_ = east_bndy_value;
}

/*! Prints a steady 1D problem. */

/*! Prints a steady 1D problem.
 */
void MTK_Steady1DProblem::MTKPrint(void) {

  cout << "Problem of interest:" << endl;
  cout << "\t" << this->alpha_ << " u(a) - " << this->beta_ << " u'(a) = " <<
    this->west_bndy_value_ << endl;
  cout << "\t" << this->gamma_ << " u(b) + " << this->delta_ << " u'(b) = " <<
    this->east_bndy_value_ << endl;
}
