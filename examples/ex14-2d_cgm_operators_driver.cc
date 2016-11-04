/*!
\file ex14-2d_cgm_operators_driver.cc

\brief Example on using the 2D CGM Operators in the MTK.

Example on using the 2D CGM Operators in the MTK.

\date: Thursday, September 12, 2013
\version: 2013-09-12.

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

#include <cstdlib>

#include <iostream>
#include <iomanip>

#include "mtk.h"

using namespace std;

/*!
\brief MAIN module.
*/
int main () {

  cout << "Testing 2D CGM Operators." << endl;

  MTK_2DCGMGradient *gg;
  MTK_2DCGMDivergence *dd;

  gg = new MTK_2DCGMGradient(2, 3, 3, 0.6667, 0.6667);
  dd = new MTK_2DCGMDivergence(2, 3, 3, 0.6667, 0.6667);

  gg->Print(MTK_DENSE);
  dd->Print(MTK_DENSE);

  return EXIT_SUCCESS;
}
