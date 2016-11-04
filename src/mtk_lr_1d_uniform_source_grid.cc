/*!
 * \file mtk_lr_1d_uniform_source_grid.cc
 *
 * \brief Constructs a default MTK_LR1DUniformSourceGrid.
 *
 * Constructs a default MTK_LR1DUniformSourceGrid. Note how the default order is
 * at least 2.
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

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <iomanip>

#include "mtk_1d_node.h"
#include "mtk_lr_1d_uniform_staggered_grid.h"
#include "mtk_lr_1d_uniform_source_grid.h"

using namespace std;

/*! Default constructor.  */
/*! This procedure CONSTRUCTS a default instance. */
MTK_LR1DUniformSourceGrid::MTK_LR1DUniformSourceGrid(
    MTK_LR1DUniformStaggeredGrid *grid,
    double west_bndy_value,
    double east_bndy_value,
    double (*source)(double)) {

  int ii;
  int nn;

  nodes_ = NULL;
  west_boundary_node_ = grid->west_boundary_node();
  east_boundary_node_ = grid->east_boundary_node();
  number_of_cells_ = grid->number_of_cells();
  nodes_grid_step_size_ = grid->nodes_grid_step_size();
  west_bndy_value_ = west_bndy_value;
  east_bndy_value_ = east_bndy_value;

  /* Allocate space for the values: */
  nn = number_of_cells_;
  nodes_ = (double *) calloc ((nn + 2),sizeof(double));
  if (nodes_ == (double *) NULL) {
    fprintf(stderr, "Allocation error at line %d.\n", __LINE__ - 2);
  }
  /* Assign the values at the boundaries: */
  nodes_[0] = west_bndy_value_;
  nodes_[(nn + 2) - 1] = east_bndy_value_;
  /* Fill the interior nodes using the provided function: */
  for (ii = 1; ii < number_of_cells_ + 1; ii++) {
    nodes_[ii] = source(((grid->nodes())[ii])->value());
  }
}

/*! Prints a MTK_LR1DUniformSourceGrid as a 1D array. */

/*! Prints a MTK_LR1DUniformSourceGrid formated as a one-dimensional
 * array.
 * \todo Monday, March 19 2012 10:50 PM: Add groups_of as a parameter.
 */
void MTK_LR1DUniformSourceGrid::Print1DArray(void) {

  int ii;
  int groups_of;

  groups_of = 10;
  cout << "S([a,b]) = [";
  for (ii = 0; ii < (number_of_cells_ + 2); ii++) {
    if ((ii != 0) && (ii % groups_of == 0)) {
      cout <<  endl;
    }
    cout << setw(10) << nodes_[ii] << " ";
  }
  cout << " ]'" << endl;
  cout << endl;
}
