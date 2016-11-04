/*!
 * \file mtk_lr_1d_uniform_nodal_grid.cc
 *
 * \brief Constructs a default MTK_LR1DUniformNodalGrid object.
 *
 * This file contains the definition of the default constructor for the class
 * MTK_LR1DUniformNodalGrid.
 *
 * \date: Monday, September 03, 2012
 * \version: 2012-09-03.
 * \author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
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

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <iomanip>

#include "mtk.h"
#include "mtk_1d_node.h"
#include "mtk_lr_1d_uniform_nodal_grid.h"

using namespace std;

/*! Constructs a default MTK_LR1DUniformNodalGrid. */

/*! Constructs a default MTK_LR1DUniformNodalGrid. */
MTK_LR1DUniformNodalGrid::MTK_LR1DUniformNodalGrid(void):
  nodes_(NULL),
  start_boundary_node_(0.0),
  final_boundary_node_(0.0),
  nodes_grid_step_size_(0.0),
  number_of_nodes_(0){

}

/*! Constructs a MTK_LR1DUniformNodalGrid. */

/*! Constructs a MTK_LR1DUniformNodalGrid based on its given boundaries.
 * This procedure has Debugging info!
 * \param start_boundary_node INPUT: Left boundary.
 * \param final_boundary_node INPUT: Right boundary.
 * \param nodes_grid_step_size INPUT: Step size.
 */
MTK_LR1DUniformNodalGrid::MTK_LR1DUniformNodalGrid(
  double start_boundary_node,
  double final_boundary_node,
  double nodes_grid_step_size) {

  int ii;

  this->nodes_ = NULL;
  this->start_boundary_node_ = start_boundary_node;
  this->final_boundary_node_ = final_boundary_node;
  this->nodes_grid_step_size_ = nodes_grid_step_size;
  this->number_of_nodes_ = (int)
  (
    (this->final_boundary_node_ - this->start_boundary_node_)/
      (this->nodes_grid_step_size_ + 2.0)
  );

#if ( Debuglevel >= 1 )
printf(
"Create grid: from %lf to %lf with nodes_grid_step_size = %lf and %i nodes: \n",
this->start_boundary_node_,
this->final_boundary_node_,
this->nodes_grid_step_size_,
this->number_of_nodes_);
#endif

  this->nodes_ =
    (MTK_1DNode**) malloc(this->number_of_nodes_*sizeof(MTK_1DNode));
  if (this->nodes_ == (MTK_1DNode**) NULL) {
    printf("ERROR 01: Nodal grid 1-D array could not be allocated!\n");
    printf("Exiting...\n");
  }
  /* Fill in the nodes values: */
  (this->nodes_)[0] = new MTK_1DNode(this->start_boundary_node_);
  for (ii = 1; ii < this->number_of_nodes_; ii++) {
    nodes_[ii] = new MTK_1DNode(
      (this->nodes_[ii - 1])->value() + this->nodes_grid_step_size_);
  }

#if ( Debuglevel >= 1 )
/* Print them: */
for (ii = 0; ii < this->number_of_nodes_; ii++) {
  printf("%lf ", (nodes_[ii])->value());
}
printf("\n\n");
#endif

}

/*! Prints a MTK_LR1DUniformNodalGrid as a 1D array. */

/*! Prints a MTK_LR1DUniformNodalGrid formated as a one-dimensional
 * array.
 * \todo Monday, March 19 2012 10:50 PM: Add groups_of as a parameter.
 */
void MTK_LR1DUniformNodalGrid::Print1DArray(void) {

  int ii;
  int groups_of;

  groups_of = 10;
  for (ii = 0; ii < this->number_of_nodes_; ii++) {
    if ((ii != 0) && (ii % groups_of == 0)) {
      cout <<  endl;
    }
    cout << setw(10) << (nodes_[ii])->value() << " ";
  }
  cout << endl;
}
