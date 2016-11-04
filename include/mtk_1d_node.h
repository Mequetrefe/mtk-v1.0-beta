/*!
  \file mtk_1d_node.h

  \brief Includes the definition of the class MTK_1DNode.

  This file contains the definition of the class MTK_1DNode.
  This class implements a 1D node, that is, a node as it is defined in 1D grid.
  It contains a value and an associated unit.

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

#ifndef MTK_INCLUDE_MTK_1D_NODE_H
#define MTK_INCLUDE_MTK_1D_NODE_H

/*! \brief Definition of a 1D node.

  Implements a node in 1D grid. A node in one dimension is nothing but a number
  which has an associated unit.
*/
class MTK_1DNode {
  public:
    /*! Default constructor.  */
    /*! This functions CONSTRUCTS a default CRS sparse matrix. */
    MTK_1DNode(void) {
      this->value_ = 0.0;
    }

    /*! Default constructor.  */
    /*! This functions CONSTRUCTS a default CRS sparse matrix. */
    MTK_1DNode(double value) {
      this->value_ = value;
    }

    /*! Default constructor.  */
    /*! This functions CONSTRUCTS a default CRS sparse matrix. */
    void set_value(double value) {
      this->value_ = value;
    }

    /*! Default constructor.  */
    /*! This functions CONSTRUCTS a default CRS sparse matrix. */
    double value(void) {
      return this->value_;
    }

  private:
    double value_;  /*!< Stored value. */
};

#endif
