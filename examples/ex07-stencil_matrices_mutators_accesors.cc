/*!
\file ex07-stencil_matrices_mutators_accesors.cc

\brief Example on accessing data within the stencil matrices.

This file performs accessing operations within the stencil matrices.

\date: Sunday, September 09, 2012
\version: 2012-09-03.

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

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "mtk.h"

using namespace std;

/*!
\brief MAIN module.
*/
int main () {

  //!Variables:
  int nrows;
  int ncols;
  MTK_1DArrayDenseStencilMatrix *MM;

  //! Begin execution:
  cout << "Example driver #7." << endl;
  cout << "Accessing data in stencil matrices." << endl;

  nrows = 5;
  ncols = 5;
  MM = new MTK_1DArrayDenseStencilMatrix(nrows, ncols);
  MM->Print(0);
  MM->SetValue(0,0,-1.0);
  MM->Print(0);

  return EXIT_SUCCESS;
}
