/*!
  \file mtk_grid_manager.h

  \brief Includes the definition of the class MTKDenseMatrix.

  This file contains the definition of the class MTKDenseMatrix.
  This class implements a sparse matrix using Compressed Column Scheme (CCS) to
  ensure storing only the non-zero elements of a given dense matrix. As a
  consequence of the selected storage scheme, no assumptions are being made in
  terms of the actual structure of the matrix.

  \date: Sunday, August 14 2011 11:42 PM
  \version: 2011-14-08-01.
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

#ifndef MTK_INCLUDE_MTK_GRID_MANAGER_H
#define MTK_INCLUDE_MTK_GRID_MANAGER_H

typedef double MTKNumber;

class MTKGrid1D;

class MTKGridManager {

public:
  /*! Default constructor.  */
  /*! This functions CONSTRUCTS a default CRS sparse matrix. */
  MTKGridManager(){int dummy_ = 0; dummy_++;}

  MTKGrid1D* MultData1D (MTKNumber mult, MTKGrid1D* grid_in);

  MTKGrid1D* CGM1DDiv ( MTKGrid1D* const grid_in );

  MTKGrid1D* Set2Sum1D (MTKGrid1D* const grid_in1,MTKGrid1D* const grid_in2,
                        MTKNumber c);

  bool Apply1DPeriodic(MTKGrid1D *grid_in);

  ~MTKGridManager();

private:
  int dummy_;

};

#endif
