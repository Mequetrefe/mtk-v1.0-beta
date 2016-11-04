/*!
  \file mtk_lr_1d_uniform_source_grid.h

  \brief Includes the definition of the class MTK_LR1DUniformSourceGrid.

  This file contains the definition of a class implementing a grid containing
  data about the source term.

  \date: Monday, September 03, 2012
  \version: 2011-09-03.
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

#ifndef MTK_INCLUDE_MTK_LR_1D_SOURCE_GRID_H
#define MTK_INCLUDE_MTK_LR_1D_SOURCE_GRID_H

class MTK1DNode;
class MTK_LR1DUniformStaggeredGrid;

/*! \brief Includes the definition of the class MTK_LR1DUniformSourceGrid.

  This class implements a one-dimensional uniform source grid of values with
  data type double.This grid contains information about the source term
  of a given PDE of interest.
  \todo Monday, September 03, 2012: Move SolveInitialValue to a solver class.
  \todo Monday, September 03, 2012: Change double to MTKNode in the boundaries.
*/
class MTK_LR1DUniformSourceGrid {

public:
  /*! Default constructor.  */
  /*! This procedure CONSTRUCTS a default instance. */
  MTK_LR1DUniformSourceGrid(void):
    nodes_(0x0),
    west_boundary_node_(0.0),
    east_boundary_node_(0.0),
    nodes_grid_step_size_(0.0),
    number_of_cells_(0){}

  /*! Constructs a source grid with start, final and step size values. */
  /*! This procedure CONSTRUCTS an instance using the most commonly used
    parameters given the context.
  */
  MTK_LR1DUniformSourceGrid(MTK_LR1DUniformStaggeredGrid *grid,
                           double west_bndy_value,
                           double east_bndy_value,
                           double (*source)(double));

  /*! Default destructor.  */
  /*! This procedure DESTRUCTS an already created instance. */
  ~MTK_LR1DUniformSourceGrid() {
    delete [] nodes_;
  }

  /*! Print formatted as a one-dimensional array.  */
  /*! This procedure prints the nodal grid in the format of a 1D array. */
  void Print1DArray(void);

  /*! Gets the nodes.  */
  /*! Gets the nodes. */
  double* nodes() { return nodes_; }

  /*! Gets the value form the \f$ i \f$-th node.  */
  /*! Gets the nodes. */
  double GetValue(int ii) { return nodes_[ii]; }

private:
  double *nodes_;               /*!< Collection of nodes. */
  double west_boundary_node_;   /*!< West boundary node. */
  double west_bndy_value_;      /*!< Value at west boundary node. */
  double east_boundary_node_;   /*!< East boundary node. */
  double east_bndy_value_;      /*!< Value at west boundary node. */
  double nodes_grid_step_size_; /*!< Step size defining the grid. */
  int number_of_cells_;         /*!< Number of cells. */
};

#endif
