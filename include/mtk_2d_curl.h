/*!
\file mtk_2d_curl.h

\brief Includes the definition of the class MTK_2DCurl.

This class implements a 2D CULR matrix operator, based on the

\date: Monday, September 03, 2012

\version: 2012-09-03.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at sciences dot sdsu dot edu

\bug No known bugs.
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

#ifndef MTK_INCLUDE_MTK_2D_CURL_H
#define MTK_INCLUDE_MTK_2D_CURL_H

class MTK_2DCurl {
 public:
  friend class MTK_DenseMatrix;
  
  MTK_2DCurl();
  MTK_2DCurl(const MTK_2DCurl&);
  MTK_2DCurl(int num_x_cells, int num_y_cells);
  MTK_2DCurl(int num_x_cells, int num_y_cells, int order_accuracy);
  MTK_2DCurl(int num_x_cells,
            int num_y_cells,
            int order_accuracy,
            MTK_Real mimetic_tolerance);
  int order_accuracy();
  MTK_Real mimetic_tolerance();
  MTK_Real* GetMatrixForm();
   
 private:
  MTK_DenseMatrix curl_;
  int order_accuracy_;
  MTK_Real mimetic_tolerance_;
};

#endif
