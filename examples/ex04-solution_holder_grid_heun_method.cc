/*!
\file ex04-solution_holder_grid_heun_method.cc

\brief Example solving an IVP and visualizing the solution.

Example solving an IVP and visualizing the solution.

\date: Monday, September 03, 2012
\version: 2011-09-03.

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

#include <cstdlib>
#include <cmath>

#include <iostream>

#include "mtk.h"

using namespace std;

/*!
\brief Creation of the function implementing the source term of the problem.
*/
double source_term(double tt, double xx) {

  return tt + sin(tt);
}

/*!
\brief Creation of the function implementing the known analytical solution.
*/
double known_solution(double tt) {

  return tt + sin(tt);
}

/*!
\brief MAIN module.
*/
int main () {

  //! Variables:
  bool save_plot;
  bool converge;
  double aa;
  double bb;
  double time_step;
  double initial_value;

  MTK_LR1DUniformNodalGrid *nodal_grid;
  MTK_LR1DUniformKnownSolutionHolderGrid *solution_grid;
  MTK_1DPlotProperties *props;
  MTK_1DPlotter *plotter;

  //! Begin execution:
  cout << "Example driver." << endl;
  cout << "Creating a solution holder grid using Heun's Method." << endl;

  aa = 0.0;
  bb = 5.0;
  time_step = 0.1;
  initial_value = 1.0;

  nodal_grid = new MTK_LR1DUniformNodalGrid(aa, bb, time_step);
  cout << "Nodal grid:" << endl;
  nodal_grid->Print1DArray();

  solution_grid = new MTK_LR1DUniformKnownSolutionHolderGrid(aa, bb, time_step);
  cout << "Solution holder grid:" << endl;

  solution_grid->Print1DArray();

  converge = solution_grid->SolveInitialValue(MTK_HEUN,
                                              initial_value,
                                              &source_term);

  if (converge) {;}

  plotter = new MTK_1DPlotter(solution_grid->independent(),
                              solution_grid->dependent(),
                              solution_grid->number_of_nodes());
  save_plot = false;
  props = new MTK_1DPlotProperties("t",
                                   "f(t)",
                                   "Attained solution",
                                   "linespoints", "red");
  plotter->set_plot_properties(props,save_plot, false, MTK_PNG);
  plotter->See();

  solution_grid = new MTK_LR1DUniformKnownSolutionHolderGrid(aa, bb,
                                                             time_step,
                                                             &known_solution);
  cout << "Solution holder grid:" << endl;
  solution_grid->Print1DArray();

  plotter = new MTK_1DPlotter(solution_grid->independent(),
                              solution_grid->dependent(),
                              solution_grid->number_of_nodes());
  save_plot = false;
  props = new MTK_1DPlotProperties("t",
                                   "F(t)",
                                   "Right-hand side",
                                   "linespoints",
                                   "green");
  plotter->set_plot_properties(props,save_plot,false,MTK_PNG);
  plotter->See();
}
