/*!
 * \file mtk_lr_1d_uniform_known_solution_holder_grid.cc
 *
 * \brief Definition of the default constructor for the
 * MTK_LR1DUniformKnownSolutionHolderGrid class.
 *
 * This file contains the definition of the default constructor for the
 * MTK_LR1DUniformKnownSolutionHolderGrid.
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
#include <cmath>

#include <iostream>
#include <iomanip>


#include "mtk_1d_node.h"
#include "mtk_lr_1d_uniform_known_solution_holder_grid.h"
#include "mtk_lr_1d_uniform_staggered_grid.h"

using namespace std;
using namespace mtk;
/*! Constructs a default MTK_LR1DUniformKnownSolutionHolderGrid instance.
 */

/*! Constructs a default MTK_LR1DUniformKnownSolutionHolderGrid instance.
 */
MTK_LR1DUniformKnownSolutionHolderGrid::
  MTK_LR1DUniformKnownSolutionHolderGrid(void):
  nodes_(NULL),
  solution_(NULL),
  start_boundary_node_(0.0),
  final_boundary_node_(0.0),
  nodes_grid_step_size_(0.0),
  number_of_nodes_(0){
}

/*! Constructs a MTK_LR1DUniformKnownSolutionHolderGrid instance.
 */

 /*! Constructs a MTK_LR1DUniformKnownSolutionHolderGrid instance from a
  * given discretized interval information of the form \f$[a,b]\;b>a \f$.
  * \param start_boundary_node INPUT: West boundary node, i.e. \f$ a \f$.
  * \param final_boundary_node INPUT: East boundary node, i.e. \f$ b \f$.
  * \param nodes_grid_step_size INPUT: Grid step size, i.e. \f$ \Delta k \f$,
  * where stands for the discretized physical quantity of interest,
 */
MTK_LR1DUniformKnownSolutionHolderGrid::
  MTK_LR1DUniformKnownSolutionHolderGrid(
    double start_boundary_node,
    double final_boundary_node,
    double nodes_grid_step_size) {

  int ii;

  this->nodes_ = NULL;
  this->start_boundary_node_ = start_boundary_node;
  this->final_boundary_node_ = final_boundary_node;
  this->nodes_grid_step_size_ = nodes_grid_step_size;
  this->number_of_nodes_ = (int)
    ((this->final_boundary_node_ - this->start_boundary_node_)/
    this->nodes_grid_step_size_ + 2.0);

#if ( Debuglevel >= 1 )
printf(
"Create grid: from %lf to %lf with nodes_grid_step_size = %lf and %i nodes: \n",
this->start_boundary_node_,
this->final_boundary_node_,
this->nodes_grid_step_size_,
this->number_of_nodes_);
if (this->uniform_)
  cout << "The grid is an UNIFORM grid." << endl;
else
  cout << "The grid is NOT an UNIFORM grid." << endl;
#endif

  this->nodes_ =
    (MTK_1DNode**) malloc(this->number_of_nodes_*sizeof(MTK_1DNode));
  if (this->nodes_ == (MTK_1DNode**) NULL) {
  printf("ERROR 01: Solution holder grid 1-D array could not be allocated!\n");
  printf("Exiting...\n");
  }
  this->solution_ =
  (MTK_1DNode**) malloc(this->number_of_nodes_*sizeof(MTK_1DNode));
  if (this->solution_ == (MTK_1DNode**) NULL) {
  printf("ERROR 01: Solution holder grid 1-D array could not be allocated!\n");
    printf("Exiting...\n");
  }
  /* Fill in the nodes values: */
  (this->nodes_)[0] = new MTK_1DNode(this->start_boundary_node_);
  solution_[0] = new MTK_1DNode(0.0);
  for (ii = 1; ii < this->number_of_nodes_; ii++) {
    nodes_[ii] =
      new MTK_1DNode(
        (this->nodes_[ii - 1])->value() + this->nodes_grid_step_size_);
    solution_[ii] = new MTK_1DNode(0.0);
  }

#if ( Debuglevel >= 1 )
/* Print them: */
for (ii = 0; ii < this->number_of_nodes_; ii++) {
  printf("%lf ", (nodes_[ii])->value());
}
printf("\n\n");
#endif

}

/*! Constructs a MTK_LR1DUniformKnownSolutionHolderGrid instance.
 */

 /*! Constructs a MTK_LR1DUniformKnownSolutionHolderGrid instance from a
  * given discretized interval information of the form \f$[a,b]\;b>a \f$.
  * \param start_boundary_node INPUT: West boundary node, i.e. \f$ a \f$.
  * \param final_boundary_node INPUT: East boundary node, i.e. \f$ b \f$.
  * \param nodes_grid_step_size INPUT: Grid step size, i.e. \f$ \Delta k \f$,
  * where stands for the discretized physical quantity of interest,
 */
MTK_LR1DUniformKnownSolutionHolderGrid::
  MTK_LR1DUniformKnownSolutionHolderGrid(
    double start_boundary_node,
    double final_boundary_node,
    double nodes_grid_step_size,
    double (*known_sol)(double)) {

  int ii;

  this->nodes_ = NULL;
  this->start_boundary_node_ = start_boundary_node;
  this->final_boundary_node_ = final_boundary_node;
  this->nodes_grid_step_size_ = nodes_grid_step_size;
  this->number_of_nodes_ = (int)
    ((this->final_boundary_node_ - this->start_boundary_node_)/
    this->nodes_grid_step_size_ + 2.0);

#if ( Debuglevel >= 1 )
printf(
"Create grid: from %lf to %lf with nodes_grid_step_size = %lf and %i nodes: \n",
this->start_boundary_node_,
this->final_boundary_node_,
this->nodes_grid_step_size_,
this->number_of_nodes_);
if (this->uniform_)
  cout << "The grid is an UNIFORM grid." << endl;
else
  cout << "The grid is NOT an UNIFORM grid." << endl;
#endif

  this->nodes_ =
    (MTK_1DNode**) malloc(this->number_of_nodes_*sizeof(MTK_1DNode));
  if (this->nodes_ == (MTK_1DNode**) NULL) {
  printf("ERROR 01: Solution holder grid 1-D array could not be allocated!\n");
    printf("Exiting...\n");
  }
  this->solution_ =
  (MTK_1DNode**) malloc(this->number_of_nodes_*sizeof(MTK_1DNode));
  if (this->solution_ == (MTK_1DNode**) NULL) {
  printf("ERROR 01: Solution holder grid 1-D array could not be allocated!\n");
    printf("Exiting...\n");
  }

  /* Fill in the nodes values: */
  nodes_[0] = new MTK_1DNode(this->start_boundary_node_);
  solution_[0] = new MTK_1DNode(known_sol(nodes_[0]->value()));
  for (ii = 1; ii < this->number_of_nodes_; ii++) {
    nodes_[ii] = new MTK_1DNode(
      (this->nodes_[ii - 1])->value() + this->nodes_grid_step_size_);
    solution_[ii] = new MTK_1DNode(known_sol(nodes_[ii]->value()));
  }

#if ( Debuglevel >= 1 )
/* Print them: */
for (ii = 0; ii < this->number_of_nodes_; ii++) {
  printf("%lf ", (nodes_[ii])->value());
}
printf("\n\n");
#endif

}

/*! Constructs a MTK_LR1DUniformKnownSolutionHolderGrid instance.
 */

 /*! Constructs a MTK_LR1DUniformKnownSolutionHolderGrid instance from a
  * given discretized interval information of the form \f$[a,b]\;b>a \f$.
  * \param start_boundary_node INPUT: West boundary node, i.e. \f$ a \f$.
  * \param final_boundary_node INPUT: East boundary node, i.e. \f$ b \f$.
  * \param nodes_grid_step_size INPUT: Grid step size, i.e. \f$ \Delta k \f$,
  * where stands for the discretized physical quantity of interest,
 */
MTK_LR1DUniformKnownSolutionHolderGrid::
  MTK_LR1DUniformKnownSolutionHolderGrid(
    double start_boundary_node,
    double final_boundary_node,
    int number_of_cells,
    double (*known_sol)(double)) {

  int ii;

  this->nodes_ = NULL;
  this->start_boundary_node_ = start_boundary_node;
  this->final_boundary_node_ = final_boundary_node;
  this->number_of_nodes_ = number_of_cells;
  this->nodes_grid_step_size_ =
  (final_boundary_node_ - start_boundary_node_)/number_of_nodes_;

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
  printf("ERROR 01: Solution holder grid 1-D array could not be allocated!\n");
    printf("Exiting...\n");
  }
  this->solution_ =
  (MTK_1DNode**) malloc(this->number_of_nodes_*sizeof(MTK_1DNode));
  if (this->solution_ == (MTK_1DNode**) NULL) {
  printf("ERROR 01: Solution holder grid 1-D array could not be allocated!\n");
    printf("Exiting...\n");
  }

  /* Fill in the nodes values: */
  nodes_[0] = new MTK_1DNode(this->start_boundary_node_);
  solution_[0] = new MTK_1DNode(known_sol(nodes_[0]->value()));
  for (ii = 1; ii < this->number_of_nodes_; ii++) {
    nodes_[ii] = new MTK_1DNode(
      (this->nodes_[ii - 1])->value() + this->nodes_grid_step_size_);
    solution_[ii] = new MTK_1DNode(known_sol(nodes_[ii]->value()));
  }

#if ( Debuglevel >= 1 )
  /* Print them: */

  for (ii = 0; ii < this->number_of_nodes_; ii++) {
    printf("%lf ", (nodes_[ii])->value());
  }
  printf("\n\n");
#endif
}

/*! Constructs a MTK_LR1DUniformKnownSolutionHolderGrid instance.
 */

 /*! Constructs a MTK_LR1DUniformKnownSolutionHolderGrid instance from a
  * given discretized interval information of the form \f$[a,b]\;b>a \f$.
  * \param start_boundary_node INPUT: West boundary node, i.e. \f$ a \f$.
  * \param final_boundary_node INPUT: East boundary node, i.e. \f$ b \f$.
  * \param nodes_grid_step_size INPUT: Grid step size, i.e. \f$ \Delta k \f$,
  * where stands for the discretized physical quantity of interest,
 */
MTK_LR1DUniformKnownSolutionHolderGrid::
  MTK_LR1DUniformKnownSolutionHolderGrid(
    MTK_LR1DUniformStaggeredGrid *grid,
    double (*known_sol)(double)) {

  int ii;

  this->nodes_ = NULL;
  this->start_boundary_node_ = grid->west_boundary_node();
  this->final_boundary_node_ = grid->east_boundary_node();
  this->number_of_nodes_ = grid->number_of_cells() + 2;
  this->nodes_grid_step_size_ = grid->nodes_grid_step_size();

#if ( DebugLevel >= 1 )
printf(
"Create grid: from %lf to %lf with nodes_grid_step_size = %lf and %i nodes: \n",
this->start_boundary_node_,
this->final_boundary_node_,
this->nodes_grid_step_size_,
this->number_of_nodes_);
#endif

  this->nodes_ = grid->nodes();

  this->solution_ =
  (MTK_1DNode**) malloc(this->number_of_nodes_*sizeof(MTK_1DNode));
  if (this->solution_ == (MTK_1DNode**) NULL) {
  printf("ERROR 01: Solution holder grid 1-D array could not be allocated!\n");
    printf("Exiting...\n");
  }

  /* Fill in the solution values: */
  solution_[0] = new MTK_1DNode(known_sol(nodes_[0]->value()));
  for (ii = 1; ii < this->number_of_nodes_; ii++) {
    solution_[ii] = new MTK_1DNode(known_sol(nodes_[ii]->value()));
  }

#if ( Debuglevel >= 1 )
  /* Print them: */
  for (ii = 0; ii < this->number_of_nodes_; ii++) {
    printf("%lf ", (nodes_[ii])->value());
  }
  printf("\n\n");
#endif
}

/*! Get the \f$ i \f$-th element of the solution. */
/*! Get the \f$ i \f$-th element of the solution. */
double MTK_LR1DUniformKnownSolutionHolderGrid::GetValue(int ii) {

  if (ii < 0 || ii > number_of_nodes_ + 2) {
    cout << "ERROR: Value out of range!" << endl;
    return -1000000000.0;
  } else {
    return (solution_[ii])->value();
  }
}

/*! Prints a MTK_LR1DUniformKnownSolutionHolderGrid as a 1D array. */

/*! Prints a MTK_LR1DUniformKnownSolutionHolderGrid formated as a
 * one-dimensional array.
 * \todo Wednesday, April 4 2012 11:51 PM: Add groups_of as a parameter.
 */
void MTK_LR1DUniformKnownSolutionHolderGrid::Print1DArray(void) {

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

  for (ii = 0; ii < this->number_of_nodes_; ii++) {
    if ((ii != 0) && (ii % groups_of == 0)) {
      cout <<  endl;
    }
    cout << setw(10) << (solution_[ii])->value() << " ";
  }
  cout << endl;
}

/*! Solves an Initial Value Problem on the grid object.
 */

 /*! Solves an Initial Value Problem on the grid object.
 * \param method INPUT: Selected method to solve for an IVP.
 * \param initial_value INPUT: Required initial value.
 * \todo Monday, September 03, 2012: Code other methods.
 */
bool MTK_LR1DUniformKnownSolutionHolderGrid::
  SolveInitialValue(
      MTK_ODEMethod method,
      double initial_value,
      double (*rhs)(double, double)) {

  int ii;
  double aux1;
  double yy;
  double hh;

  hh = this->nodes_grid_step_size_;

  switch (method) {
    case MTK_FORWARD_EULER: {
      break;
    }
    case MTK_BACKWARD_EULER: {
      break;
    }
    case MTK_TRAPEZOIDAL: {
      break;
    }
    case MTK_HEUN: {
      // Set initial value:
      (this->solution_[0])->set_value(initial_value);
      // Iterate:
      for (ii = 1; ii < this->number_of_nodes_; ii++) {
        // Predict using Forward Euler:
        yy = (this->solution_[ii-1])->value();
        aux1 = yy + hh*rhs(yy,yy);
        // Correct using Crank - Nicholson (Trapezoidal):
        (this->solution_[ii])->set_value(yy + (hh/2.0)*(rhs(yy,yy) + aux1));
      }
      break;
    }
  }
  return true;
}
