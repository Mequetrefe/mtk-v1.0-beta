/*!
  \file mtk_lr_1d_uniform_staggered_grid.h

  \brief Includes the definition of the class MTK_LR1DUniformStaggeredGrid.

  This file contains the definition of a class implementing a staggered grid.

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

#ifndef MTK_INCLUDE_MTK_LR_1D_STAGGERED_GRID_H
#define MTK_INCLUDE_MTK_LR_1D_STAGGERED_GRID_H

class MTK_1DNode;

/*! \brief One-dimensional uniform staggered grid of double.

  This class implements a one-dimensional uniform staggered grid of values with
  data type double. Notice that a nodal grid differs of a staggered grid in the
  sense that no centers are defined.
  \todo Monday, September 03, 2012: Move SolveInitialValue to a solver class.
  \todo Monday, September 03, 2012: Change double to MTKNode in the boundaries.
*/
class MTK_LR1DUniformStaggeredGrid {

public:
  /*! Default constructor.  */
  /*! This procedure CONSTRUCTS a default instance. */
  MTK_LR1DUniformStaggeredGrid(void);

  /*! Constructs a staggered grid with start, final and step size values. */
  /*! This procedure CONSTRUCTS an instance using the most commonly used
    parameters given the context.
  */
  MTK_LR1DUniformStaggeredGrid(double west_boundary_node,
                          double east_boundary_node,
                          int number_of_cells);

  /*! Default destructor.  */
  /*! This procedure DESTRUCTS an already created instance. */
  ~MTK_LR1DUniformStaggeredGrid() {
    delete [] nodes_;
  }

  /*! Print formatted as a one-dimensional array.  */
  /*! This procedure prints the nodal grid in the format of a 1D array. */
  void Print1DArray(void);

  /*! Gets the nodes.  */
  /*! Gets the nodes. */
  MTK_1DNode** nodes(void) { return nodes_; }

  /*! Gets the node at the west boundary.  */
  /*! Gets the node at the west boundary. */
  double west_boundary_node(void) { return west_boundary_node_; }

  /*! Gets the node at the east boundary.  */
  /*! Gets the node at the east boundary. */
  double east_boundary_node(void) { return east_boundary_node_; }

  /*! Gets the nodes step size.  */
  /*! Gets the nodes step size. */
  double nodes_grid_step_size(void) { return nodes_grid_step_size_; }

  /*! Print formatted as a one-dimensional array.  */
  /*! This procedure prints the nodal grid in the format of a 1D array. */
  int number_of_cells(void) { return number_of_cells_; }

private:
  MTK_1DNode **nodes_;     /*!< Collection of nodes. */
  double west_boundary_node_;   /*!< West boundary node. */
  double east_boundary_node_;   /*!< East boundary node. */
  double nodes_grid_step_size_; /*!< Step size defining the grid. */
  int number_of_cells_;         /*!< Number of cells. */
};

#endif
