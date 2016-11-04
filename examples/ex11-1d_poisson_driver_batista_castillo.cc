/*!
  \file ex11-1d_poisson_driver_batista_castillo.cc

  \brief Driver code for solving a one dimensional Poisson's equation.

  This file contains an example consisting of a driver code which uses MTK to
  solve the following one dimensional Poisson's equation over a
  specified one-dimensional numerical domain. Particularly, the code solves:
  \f[
  -lap\;f(x)=F(x),
  \f]
  where
  \f[
  F(x)=\frac{2\times 10^{6}x}{\arctan(100)(1+1\times 10^{4}x^2)^2}.
  \f]
  The problem uses specified Robin's Boundary Conditions of the form:
  \f{eqnarray*}{
f(0)+ f^{\prime}(0)&=&-\frac{100}{\arctan(100)},\\
f(1)+ f^{\prime}(1)&=&1+\frac{100}{\arctan(100)(1+1\times 10^{4})}.
  \f}

  If we consider \f$ x\in\Omega=[0,1] \f$, the problem has
  known analytical solution given by:
  \f[
    f(x)=\frac{\arctan(100x)}{\arctan(100)}.
  \f]

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

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "mtk.h"

using namespace std;

/*!
\brief Creation of the function implementing the source term of the problem.
*/
double source_term (double xx) {

  return ((2.0e6)*xx)/
    (atan(100.0)*((1.0 + (1.0e4)*(xx*xx))*(1.0 + (1.0e4)*(xx*xx))));
}

/*!
\brief Creation of the function implementing the known analytical solution.
*/
double known_solution (double xx) {

  return atan(100.0*xx)/atan(100.0);
}

/*!
\brief MAIN module.
*/
int main () {

  //! Variables:
  bool save_plot;
  double west_bndy;
  double east_bndy;
  double alpha;
  double beta;
  double west_bndy_value;
  double gamma;
  double delta;
  double east_bndy_value;
  double norm_diff_sol;
  int number_of_cell_centers;
  int desired_order;

  MTK_ToolManager *tools;
  MTK_LR1DUniformStaggeredGrid *space;
  MTK_LR1DUniformSourceGrid *source;
  MTK_1DCGMDirichletOperator *dir_comp;
  MTK_1DCGMNeumannOperator *neu_comp;
  MTK_1DCGMGradient *grad;
  MTK_1DCGMLaplacian *lap;
  MTK_Steady1DProblem *problem_to_solve;
  MTK_1DArrayDenseStencilMatrix *stencil_matrix;
  MTK_LR1DUniformKnownSolutionHolderGrid *solution_grid;
  MTK_1DPlotProperties *props;
  MTK_1DPlotter *plotter;

  //! Begin execution:
  cout << "Example driver #11." << endl;
  cout << "A Diffusion-Reaction (Poisson) Driver." << endl;
  number_of_cell_centers = 1000;
  desired_order = 2;
  west_bndy = 0.0;
  east_bndy = 1.0;
  alpha = 1.0;
  beta = 1.0;
  west_bndy_value = -100.0/atan(100.0);
  gamma = alpha;
  delta = beta;
  east_bndy_value = 1 + ( 100/(atan(100)*(1 + 1e4)) );

  tools = new MTK_ToolManager();

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

  // Describe the problem:
  problem_to_solve = new MTK_Steady1DProblem(alpha, beta, west_bndy_value,
                                             gamma, delta, east_bndy_value);
  problem_to_solve->MTKPrint();

  // Create the staggered grid:
  space = new MTK_LR1DUniformStaggeredGrid(west_bndy, east_bndy,
                                           number_of_cell_centers);
  cout << "Defined a LR 1D uniform staggered grid:" << endl << "\t";
  space->Print1DArray();

  // Create stencil matrix:
  stencil_matrix = new MTK_1DArrayDenseStencilMatrix(alpha, dir_comp,
                                                     beta, neu_comp,
                                                     grad, lap);

  // Discretize source term:
  source = new MTK_LR1DUniformSourceGrid(space, west_bndy_value,
                                         east_bndy_value, source_term);

  // Create ans solve system of equations using stencil matrix and source term:
  MTK_System system(stencil_matrix, source);
  system.Solve();

  // Create a grid holding the known analytic solution:
  solution_grid = new MTK_LR1DUniformKnownSolutionHolderGrid(space,
                                                             &known_solution);

  // Visualize it and DO NOT save the generated plot:
  plotter = new MTK_1DPlotter(solution_grid->independent(),
                              solution_grid->dependent(),
                              solution_grid->number_of_nodes());
  save_plot = false;
  props = new MTK_1DPlotProperties("x",
                                  "f(x)",
                                  "Known analytic solution (1002 cells)",
                                  "lines",
                                   "red");
  plotter->set_plot_properties(props, save_plot, false, MTK_PNG);
  plotter->See();

  // Visualize computed solution and save the generated plot:
  plotter = new MTK_1DPlotter(solution_grid->independent(),
                              system.Solution(),
                              solution_grid->number_of_nodes());
  save_plot = false;
  props =
    new MTK_1DPlotProperties("x",
                             "f(x)",
    "Computed solution: 2nd order CGM operators (1002 cells)",
                             "linespoints", "blue");
  plotter->set_plot_properties(props,save_plot, false, MTK_PNG);
  plotter->See();

  // Compute norm of the difference between known and computed solutions:
  norm_diff_sol =
    tools->RelativeNorm2Difference(solution_grid->dependent(),
                                   system.Solution(),
                                   solution_grid->number_of_nodes());
  cout << "Relative Norm 2 of the difference = " << norm_diff_sol << endl;
  return EXIT_SUCCESS;
}
