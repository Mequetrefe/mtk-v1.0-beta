/*!
  \file mtk_1d_uniform_staggered_grid.cc

  \brief Constructs a default MTK_LR1DUniformStaggeredGrid object.

  This file contains the definition of the default constructor for the class
  MTK_LR1DUniformStaggeredGrid.

  \date: Monday, September 03, 2012
  \version: 2012-09-03.
  \author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
*/
/* Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
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
#include <cstdio>

#include <iostream>
#include <iomanip>

#include "mtk_1d_grid_node.h"
#include "mtk_1d_uniform_cartesian_staggered_grid.h"
#include "mtk_roots.h"
using namespace std;
using namespace mtk;
/*! Constructs a default MTK_LR1DUniformStaggeredGrid. */

/*! Constructs a default MTK_LR1DUniformStaggeredGrid. */
MTK_1DUniformStaggeredGrid::MTK_1DUniformStaggeredGrid(void):
  nodes_(NULL),
  west_boundary_node_(0.0),
  east_boundary_node_(0.0),
  nodes_grid_step_size_(0.0),
  number_of_cells_(0){

}

/*! Constructs a MTK_LR1DUniformStaggeredGrid. */

/*! Constructs a MTK_LR1DUniformStaggeredGrid based on its given
 * boundaries.
 * This procedure has Debugging info!
 * \param west_boundary_node INPUT: Left boundary.
 * \param east_boundary_node INPUT: Right boundary.
 * \param number_of_cells INPUT: Number of cells.
 */
MTK_1DUniformStaggeredGrid::MTK_1DUniformStaggeredGrid(
  MTK_Real west_boundary_node,
  MTK_Real east_boundary_node,
  int number_of_cells) {

  int ii;

  this->nodes_ = NULL;
  this->west_boundary_node_ = west_boundary_node;
  this->east_boundary_node_ = east_boundary_node;
  this->number_of_cells_ = number_of_cells;
  this->nodes_grid_step_size_ =
  (this->east_boundary_node_ - this->west_boundary_node_)/
    this->number_of_cells_;

#if ( MTK_DEBUG_LEVEL > 0 )
printf(
"Create grid: from %lf to %lf with nodes_grid_step_size = %lf and %i nodes: \n",
this->west_boundary_node_,
this->east_boundary_node_,
this->nodes_grid_step_size_,
this->number_of_cells_);
#endif

  // We sum 2 because, for example, for 5 cells (5 cell centers), we would have
  // 5 cell centers PLUS 2 boundaries:
  this->nodes_ =
    (MTK_1DNode**)
      malloc((this->number_of_cells_ + 2)*sizeof(MTK_1DNode));
  if (this->nodes_ == (MTK_1DNode**) NULL) {
    printf("ERROR 01: Staggered grid 1-D array could not be allocated!\n");
    printf("Exiting...\n");
  }
  /* Fill in the nodes values: */
  /* The fact that we divide by 2.0 when storing the first cell center, is
   because of the staggering! Cells centers are shifted with respect the
   actual nodes (cell edges). */
  (this->nodes_)[0] = new MTK_1DNode(this->west_boundary_node_);
  (this->nodes_)[1] =
    new MTK_1DNode(
      this->west_boundary_node_ + this->nodes_grid_step_size_/2.0);
  for (ii = 2; ii <= this->number_of_cells_; ii++) {
    nodes_[ii] = new MTK_1DNode((this->nodes_[ii - 1])->value() +
      this->nodes_grid_step_size_);
  }
  (this->nodes_)[this->number_of_cells_ + 1] =
    new MTK_1DNode(this->east_boundary_node_);
#if ( MTK_DEBUG_LEVEL > 0 )
/* Print them: */
for (ii = 0; ii < this->number_of_cells_ + 2; ii++) {
  printf("%lf ", (nodes_[ii])->value());
}
printf("\n\n");
#endif

}

/*! Prints a MTK_LR1DUniformStaggeredGrid as a 1D array. */

/*! Prints a MTK_LR1DUniformStaggeredGrid formated as a one-dimensional
 * array.
 * \todo Monday, March 19 2012 10:50 PM: Add groups_of as a parameter.
 */
void MTK_1DUniformStaggeredGrid::Print1DArray(void) {

  int ii;
  int groups_of;

  groups_of = 10;
  cout << "[a,b] = " << endl;
  for (ii = 0; ii < this->number_of_cells_ + 2; ii++) {
    if ((ii != 0) && (ii % groups_of == 0)) {
      cout <<  endl;
    }
    cout << setw(10) << (nodes_[ii])->value() << " ";
  }
  cout << endl;
}
