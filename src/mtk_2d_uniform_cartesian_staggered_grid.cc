/*!
  \file mtk_2d_uniform_staggered_grid.cc

  \brief Constructs a default MTK_2DUniformStaggeredGrid object.

  This file contains the definition of the default constructor for the class
  MTK_2DUniformStaggeredGrid.

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
#include "mtk_2d_uniform_cartesian_staggered_grid.h"
#include "mtk_roots.h"
using namespace std;
using namespace mtk;
/*! Constructs a default MTK_LR1DUniformStaggeredGrid. */

/*! Constructs a default MTK_LR1DUniformStaggeredGrid. */
MTK_2DUniformStaggeredGrid::MTK_2DUniformStaggeredGrid(void):
  nodes_X_(NULL),

  nodes_Y_(NULL),
  
  west_boundary_node_(0.0),
  east_boundary_node_(0.0),
  nodes_grid_step_size_X_(0.0),
  number_of_cells_X_(0),

  south_boundary_node_(0.0),
  north_boundary_node_(0.0),
  nodes_grid_step_size_Y_(0.0),
  number_of_cells_Y_(0){

}

/*! Constructs a MTK_LR1DUniformStaggeredGrid. */

/*! Constructs a MTK_LR1DUniformStaggeredGrid based on its given
 * boundaries.
 * This procedure has Debugging info!
 * \param west_boundary_node INPUT: Left boundary.
 * \param east_boundary_node INPUT: Right boundary.
 * \param number_of_cells INPUT: Number of cells.
 */
MTK_2DUniformStaggeredGrid::MTK_2DUniformStaggeredGrid(
  MTK_Real west_boundary_node,
  MTK_Real east_boundary_node,
  int number_of_cells_X,
  MTK_Real south_boundary_node,
  MTK_Real north_boundary_node,
  int number_of_cells_Y
    
) {

  int ii;

  this->nodes_X_ = NULL;
  this->west_boundary_node_ = west_boundary_node;
  this->east_boundary_node_ = east_boundary_node;
  this->number_of_cells_X_ = number_of_cells_X;
  this->nodes_grid_step_size_X_ =
  (this->east_boundary_node_ - this->west_boundary_node_)/
    this->number_of_cells_X_;

  
  this->nodes_Y_ = NULL;
  this->south_boundary_node_ = south_boundary_node;
  this->north_boundary_node_ = north_boundary_node;
  this->number_of_cells_Y_ = number_of_cells_Y;
  this->nodes_grid_step_size_Y_ =
  (this->north_boundary_node_ - this->south_boundary_node_)/
    this->number_of_cells_Y_;

#if ( MTK_DEBUG_LEVEL > 0 )
printf(
"Create 2Dgrid: from X %lf to %lf with nodes_grid_step_size = %lf and %i nodes: \n",
this->west_boundary_node_,
this->east_boundary_node_,
this->nodes_grid_step_size_X_,
this->number_of_cells_X_);
printf(
    "\nand Y from %lf to %lf with nodes_grid_step_size = %lf and %i nodes: \n",
this->south_boundary_node_,
this->north_boundary_node_,
this->nodes_grid_step_size_Y_,
this->number_of_cells_Y_);
#endif

  // We sum 2 because, for example, for 5 cells (5 cell centers), we would have
  // 5 cell centers PLUS 2 boundaries:
  this->nodes_X_ =
    (MTK_1DNode**)
      malloc((this->number_of_cells_X_ + 2)*sizeof(MTK_1DNode));
  if (this->nodes_X_ == (MTK_1DNode**) NULL) {
    printf("ERROR 01: Horizontal Staggered grid 2-D array could not be allocated!\n");
    printf("Exiting...\n");
  }
  this->nodes_Y_ =
    (MTK_1DNode**)
      malloc((this->number_of_cells_Y_ + 2)*sizeof(MTK_1DNode));
  if (this->nodes_Y_ == (MTK_1DNode**) NULL) {
    printf("ERROR 01: Horizontal Staggered grid 2-D array could not be allocated!\n");
    printf("Exiting...\n");
  }
  /* Fill in the nodes values: */
  /* The fact that we divide by 2.0 when storing the first cell center, is
   because of the staggering! Cells centers are shifted with respect the
   actual nodes (cell edges). */
  (this->nodes_X_)[0] = new MTK_1DNode(this->west_boundary_node_);
  (this->nodes_X_)[1] =
    new MTK_1DNode(
      this->west_boundary_node_ + this->nodes_grid_step_size_X_/2.0);
  for (ii = 2; ii <= this->number_of_cells_X_; ii++) {
    nodes_X_[ii] = new MTK_1DNode((this->nodes_X_[ii - 1])->value() +
      this->nodes_grid_step_size_X_);
  }
  (this->nodes_X_)[this->number_of_cells_X_ + 1] =
    new MTK_1DNode(this->east_boundary_node_);
#if ( MTK_DEBUG_LEVEL > 0 )
/* Print them: */
for (ii = 0; ii < this->number_of_cells_X_ + 2; ii++) {
  printf("%lf ", (nodes_X_[ii])->value());
}
printf("\n\n");
#endif

(this->nodes_Y_)[0] = new MTK_1DNode(this->south_boundary_node_);
  (this->nodes_Y_)[1] =
    new MTK_1DNode(
      this->south_boundary_node_ + this->nodes_grid_step_size_Y_/2.0);
  for (ii = 2; ii <= this->number_of_cells_Y_; ii++) {
    nodes_Y_[ii] = new MTK_1DNode((this->nodes_Y_[ii - 1])->value() +
      this->nodes_grid_step_size_Y_);
  }
  (this->nodes_Y_)[this->number_of_cells_Y_ + 1] =
    new MTK_1DNode(this->north_boundary_node_);
#if ( MTK_DEBUG_LEVEL > 0 )
/* Print them: */
for (ii = 0; ii < this->number_of_cells_Y_ + 2; ii++) {
  printf("%lf ", (nodes_Y_[ii])->value());
}
printf("\n\n");
#endif

}

/*! Prints a MTK_LR1DUniformStaggeredGrid as a 1D array. */

/*! Prints a MTK_LR1DUniformStaggeredGrid formated as a one-dimensional
 * array.
 * \todo Monday, March 19 2012 10:50 PM: Add groups_of as a parameter.
 */
void MTK_2DUniformStaggeredGrid::Print2DArray(void) {

  int ii;
  int groups_of;

  groups_of = 10;
  cout << "[a,b] = " << endl;
  for (ii = 0; ii < this->number_of_cells_X_ + 2; ii++) {
    if ((ii != 0) && (ii % groups_of == 0)) {
      cout <<  endl;
    }
    cout << setw(10) << (nodes_X_[ii])->value() << " ";
  }
  cout << endl;
  cout << "[c,d] = " << endl;
  for (ii = 0; ii < this->number_of_cells_Y_ + 2; ii++) {
    if ((ii != 0) && (ii % groups_of == 0)) {
      cout <<  endl;
    }
    cout << setw(10) << (nodes_Y_[ii])->value() << " ";
  }
  cout << endl;
}

 MTK_Real MTK_2DUniformStaggeredGrid::valuesX(int ii)
{
    return this->nodes_X_[ii]->value();
}

 MTK_Real MTK_2DUniformStaggeredGrid::valuesY(int ii)
{
    return this->nodes_Y_[ii]->value();
}

