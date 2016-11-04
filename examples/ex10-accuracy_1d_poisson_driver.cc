/*!
\file ex10-accuracy_1d_poisson_driver.cc

\brief Accuracy when solving a one dimensional Poisson's equation.

This file contains an example consisting of a driver code which uses MTK to
solve the following one dimensional Poisson's equation over a
specified one-dimensional numerical domain. Particularly, the code solves:
\f[
-lap\;p(x)=F(x),
\f]
where
\f[
F(x)=-\frac{\lambda^2\exp(\lambda x)}{\exp(\lambda)-1}.
\f]
The problem uses specified Robin's Boundary Conditions of the form:
\f{eqnarray*}{
\alpha p(a)-\beta p^{\prime}(a)&=&\omega\\
\alpha p(b)+\beta p^{\prime}(b)&=&\epsilon,
\f}
where \f$ \omega =-1 \f$, \f$ \epsilon=0 \f$, \f$ \alpha=-\exp(\lambda)\f$
and \f$ \beta=(\exp(\lambda)-1)/\lambda \f$

If we consider \f$ \Omega=[0,1] \f$ and \f$ \lambda = -1 \f$, the problem has
known analytical solution given by:
\f[
  p(x)=\frac{e^{\lambda x}-1}{e^\lambda - 1}.
\f]

Particularly, in this file, we explore the attained order of accuracy when
solving this example.

\date: Sunday, November 04, 2012
\version: 2012-11-04.

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
#include <cstdlib>
#include <cmath>

#include "mtk.h"

using namespace std;

/*!
\brief Creation of the function implementing the source term of the problem.
*/
double source_term (double xx) {

  double lambda;

  lambda = -1.0;
  return -(lambda*lambda)*exp(lambda*xx)/(exp(lambda) - 1.0);
}

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
  bool loglog;
  double centers[] = {5, 10, 20, 50, 100, 200, 250, 500};
  double sizes[8];
  double errors_interior[8];
  double errors_west[8];
  double errors_east[8];
  double west_bndy;
  double east_bndy;
  double alpha;
  double beta;
  double west_bndy_value;
  double east_bndy_value;
  double lambda;
  double norm_diff_sol;
  int number_of_cell_centers;
  int desired_order;
  int ii;

  MTK_ToolManager *tools;
  MTK_LR1DUniformStaggeredGrid *space;
  MTK_LR1DUniformSourceGrid *source;
  MTK_1DCGMDirichletOperator *dir_comp;
  MTK_1DCGMNeumannOperator *neu_comp;
  MTK_1DCGMGradient *grad;
  MTK_1DCGMLaplacian *lap;
  MTK_1DArrayDenseStencilMatrix *stencil_matrix;
  MTK_LR1DUniformKnownSolutionHolderGrid *solution_grid;
  MTK_1DPlotProperties *props;
  MTK_1DPlotter *plotter;

  //! Begin execution:
  cout << "Example driver #10." << endl;
  cout << "Accuracy of a Poisson Driver." << endl;
  desired_order = 2;
  lambda = -1.0;
  west_bndy = 0.0;
  east_bndy = 1.0;
  alpha = -exp(lambda);
  beta = (exp(lambda) - 1.0)/(lambda);
  west_bndy_value = -1.0;
  east_bndy_value = 0.0;

  tools = new MTK_ToolManager();

  for (ii = 0; ii < 8; ii++) {
    number_of_cell_centers = centers[ii];
    sizes[ii] = (abs(west_bndy - east_bndy))/number_of_cell_centers;

    // Initialize mimetic operators:
    dir_comp = new MTK_1DCGMDirichletOperator(desired_order,
                                              number_of_cell_centers,
                                              MTK_DENSE);
    neu_comp = new MTK_1DCGMNeumannOperator(desired_order,
                                            number_of_cell_centers,
                                            MTK_DENSE);
    grad = new MTK_1DCGMGradient(desired_order,
                                 number_of_cell_centers,
                                 MTK_DENSE);
    lap = new MTK_1DCGMLaplacian(desired_order,
                                 number_of_cell_centers,
                                 MTK_DENSE);

    // Create the staggered grid:
    space = new MTK_LR1DUniformStaggeredGrid(west_bndy, east_bndy,
                                            number_of_cell_centers);

    // Create stencil matrix:
    stencil_matrix = new MTK_1DArrayDenseStencilMatrix(alpha, dir_comp,
                                                      beta, neu_comp,
                                                      grad, lap);

    // Discretize source term:
    source = new MTK_LR1DUniformSourceGrid(space, west_bndy_value,
                                          east_bndy_value, source_term);

    // Create & solve system of equations using stencil matrix and source term:
    MTK_System system(stencil_matrix, source);
    system.Solve();

    // Create a grid holding the known analytic solution:
    solution_grid =
    new MTK_LR1DUniformKnownSolutionHolderGrid(space, &known_solution);

    cout.precision(5);
    // Compute norm of the difference between known and computed solutions:
    norm_diff_sol =
      tools->RelativeNorm2Difference(solution_grid->dependent(),
                                    system.Solution(),
                                    solution_grid->number_of_nodes());
    errors_interior[ii] = norm_diff_sol;
    cout << "Relative Norm 2 of the difference = " << scientific <<
    norm_diff_sol << endl;

    errors_west[ii] =
      tools->MagnitudeDifferenceWest(solution_grid->dependent(),
                                    system.Solution(),
                                    solution_grid->number_of_nodes());
    cout << "Difference west= " << errors_west[ii] << endl;

    errors_east[ii] =
      tools->MagnitudeDifferenceEast(solution_grid->dependent(),
                                    system.Solution(),
                                    solution_grid->number_of_nodes());
    cout << "Difference east= " << scientific << errors_east[ii] << endl;
  }

  plotter = new MTK_1DPlotter(centers, errors_interior, 8);
  save_plot = false;
  loglog = true;
  props =
    new MTK_1DPlotProperties("log of the amount of cell centers",
                             "log of the relative norm-2 error",
                             "Behavior of the rel. error in the interior nodes",
                             "linespoints", "blue");
  plotter->set_plot_properties(props, save_plot, loglog, MTK_PNG);
  plotter->See();

  plotter = new MTK_1DPlotter(sizes, errors_interior, 8);
  save_plot = false;
  loglog = true;
  props =
    new MTK_1DPlotProperties("log of the grid step size",
                             "log of the relative norm-2 error",
                             "Behavior of the rel. error in the interior nodes",
                             "linespoints", "blue");
  plotter->set_plot_properties(props, save_plot, loglog, MTK_PNG);
  plotter->See();

  plotter = new MTK_1DPlotter(centers, errors_west, 8);
  save_plot = false;
  loglog = true;
  props =
    new MTK_1DPlotProperties("log of the amount of cell centers",
                             "log of the relative norm-2 error",
                             "Behavior of the rel. error at the west boundary",
                             "linespoints", "blue");
  plotter->set_plot_properties(props, save_plot, loglog, MTK_PNG);
  plotter->See();

  plotter = new MTK_1DPlotter(centers, errors_east, 8);
  save_plot = false;
  loglog = true;
  props =
    new MTK_1DPlotProperties("log of the amount of cell centers",
                             "log of the relative norm-2 error",
                             "Behavior of the rel. error at the east boundary",
                             "linespoints", "blue");
  plotter->set_plot_properties(props, save_plot, loglog, MTK_PNG);
  plotter->See();

  // Create a basic table depicting the behavior of the error:
  cout << "N h int o west o east o" << endl;

  for (ii = 0; ii < 8; ii++) {
    if ( ii == 0 )
      cout << (int) centers[ii] << " & " << sizes[ii] << " & " <<
      errors_interior[ii] << " & " << "-" << " & " << errors_west[ii] <<
      " & " << "-" << " & " << errors_east[ii] << " & " << "-"<< " \\\\" <<
      endl;
    else
      cout << (int) centers[ii] << " & " << sizes[ii] << " & " <<
      errors_interior[ii] << " & " <<
      (log(errors_interior[ii - 1]) - log(errors_interior[ii]))/
        (log(sizes[ii - 1]) - log(sizes[ii])) << " & " << errors_west[ii] <<
      " & " << (log(errors_west[ii - 1]) - log(errors_west[ii]))/
        (log(sizes[ii - 1]) - log(sizes[ii])) << " & " << errors_east[ii] <<
      " & " << (log(errors_east[ii - 1]) - log(errors_east[ii]))/
        (log(sizes[ii - 1]) - log(sizes[ii]))<< " \\\\" << endl;
  }
  return EXIT_SUCCESS;
}
