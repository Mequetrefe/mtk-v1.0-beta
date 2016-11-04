/*!
\file ex06-solution_holder_grid_access.cc

\brief Example on accessing data within the grids holding solution.

This file performs accessing operations on grids containing solutions.

\date: Sunday, September 09, 2012
\version: 2012-09-03.

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

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "mtk.h"

using namespace std;

/*!
\brief Creation of the function implementing the known analytical solution.
*/
double known_solution (double xx) {

  double lambda;

  lambda = -1.0;
  return (exp(lambda*xx) - 1.0)/(exp(lambda) - 1.0);
}

/*!
\brief MAIN module.
*/
int main () {

  //! Variables:
  bool save_plot;
  double west_bndy;
  double east_bndy;
  int number_of_cell_centers;
  int groups_of;
  int ii;

  MTK_LR1DUniformStaggeredGrid *space;
  MTK_LR1DUniformKnownSolutionHolderGrid *solution_grid;
  MTK_1DPlotProperties *props;
  MTK_1DPlotter *plotter;

  //! Begin execution:
  cout << "Example driver #6." << endl;
  cout << "Accessing data in solution holder grids." << endl;

  // Grids are collection of nodes; or (computationally speaking), collection
  // of pointers to nodes. Because of this, if we want to access the values
  // within a grid, we must access first each element by means of a
  // de-referentiation operation, and then use the accessor function of a node,
  // to get the actual value.

  // Let us first create a solution holder grid, holding the discrete values of
  // a function we claim to be the known analytical solution of a given PDE. For
  // this, we must:

  // 1. Discretize the space:
  west_bndy = 0.0;
  east_bndy = 10.0;
  number_of_cell_centers = 10;
  space = new MTK_LR1DUniformStaggeredGrid(west_bndy, east_bndy,
                                                number_of_cell_centers);
  cout << "Defined a LR 1D uniform staggered grid:" << endl << "\t";
  space->Print1DArray();

  // 2. Create the solution holder grid:
  solution_grid =
  new MTK_LR1DUniformKnownSolutionHolderGrid(space,
                                                  &known_solution);
  cout << "Defined a LR 1D sol. holder grid:" << endl << "\t";
  solution_grid->Print1DArray();

  // 3. Visualize, to see if the solution makes sense:
  plotter = new MTK_1DPlotter(solution_grid->independent(),
                                   solution_grid->dependent(),
                                   solution_grid->number_of_nodes());
  save_plot = false;
  props = new MTK_1DPlotProperties("x","u(x)","Known solution","lines", "red");
  plotter->set_plot_properties(props, save_plot, false, MTK_PNG);
  plotter->See();

  // Now, what if we want to access EACH value in the grids? Say you want to
  // print them in groups of 10 per row, since they can be a lot of them:
  groups_of = 10;
  cout << "[a,b] = " << endl;
  for (ii = 0; ii < solution_grid->number_of_nodes(); ii++) {
    if ((ii != 0) && (ii % groups_of == 0)) {
      cout <<  endl;
    }
    cout << setw(10) << ((solution_grid->dependent())[ii])->value() << " ";
  }
  cout << endl;

  // Notice how we have accessed the actual value:
  // 1. go to the solution grid, and get the collection of dependent values:
  // solution_grid->dependent()
  // 2. Index the specific value we want (in this case we used an iterator ii):
  // (solution_grid->dependent())[ii]
  // 3. Since this is pointing to a node, dereference it:
  // ((solution_grid->dependent())[ii])->
  // 4. AND ask for its value:
  // ((solution_grid->dependent())[ii])->value()

  // However, this is quite intricated, therefore, we will make a method
  // for the grids, to give away any indexed element! Using it we could do:
  cout << "Again, [a,b] = " << endl;
  for (ii = 0; ii < solution_grid->number_of_nodes(); ii++) {
    if ((ii != 0) && (ii % groups_of == 0)) {
      cout <<  endl;
    }
    cout << setw(10) << solution_grid->GetValue(ii) << " ";
  }
  cout << endl;

  return EXIT_SUCCESS;
}
