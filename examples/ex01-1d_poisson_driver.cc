
#if __cplusplus == 201103L
#else
#include <iostream>
using namespace std;
int main () {
  cout << "This code HAS to be compiled to support C++11." << endl;
  cout << "On your shell, do: sudo apt-get install build-essentials" << endl;
  cout << "Exiting..." << endl;
}
#endif



/*!
\file ex01-1d_poisson_driver.cc

\brief Driver code for solving a one dimensional Poisson's equation.

This file contains an example consisting of a driver code which uses the MTK to
solve the following one dimensional Poisson's equation over a specified
one-dimensional numerical domain. Particularly, the code solves:
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

  //! Variables: 25.
  bool save_plot;             // Should we save the generated plot?
  double west_bndy;           // Spatial coordinate of the west boundary.
  double east_bndy;           // Spatial coordinate of the east boundary.
  double alpha;               // Dirichlet coefficient at the west boundary.
  double beta;                // Neumann coefficient at the west boundary.
  double west_bndy_value;     // Boundary value at the west boundary.
  double gamma;               // Dirichlet coefficient at the east boundary.
  double delta;               // Neumann coefficient at the east boundary.
  double east_bndy_value;     // Boundary value at the east boundary.
  double lambda;              // Lambda parameter.
  double norm_diff_sol;       // 2-norm of the difference for computing error.
  int number_of_cell_centers; // Number of cell centers for discretization.
  int desired_order;          // Desired order of accuracy for the solution.

  //! MTK classes: 12.
  MTK_ToolManager *tools;                 // Used for execution control.
  MTK_LR1DUniformStaggeredGrid *space;    // Staggered grid for discretization.
  MTK_LR1DUniformSourceGrid *source;      // Staggered grid for source term.
  MTK_1DCGMDirichletOperator *dir_comp;   // Dirichlet mimetic operator.
  MTK_1DCGMNeumannOperator *neu_comp;     // Neumann mimetic operator.
  MTK_1DCGMGradient *grad;                // Mimetic gradient operator.
  MTK_1DCGMLaplacian *lap;                // Mimetic Laplacian operator.
  MTK_Steady1DProblem *problem_to_solve;  // Actual problem to be solved.
  // Stencil matrix for system:
  MTK_1DArrayDenseStencilMatrix *stencil_matrix;
  // Grid that will contain the solution of the problem:
  MTK_LR1DUniformKnownSolutionHolderGrid *solution_grid;
  // Properties of the plots:
  MTK_1DPlotProperties *props;
  MTK_1DPlotter *plotter;                 // Creation of plots.

  //! 1. Begin execution (at line 134).
  cout << "Example driver number 1." << endl;
  cout << "A diffusion-reaction (Poisson) driver." << endl;

  number_of_cell_centers = 5;
  desired_order = 2;
  lambda = -1.0;
  west_bndy = 0.0;
  east_bndy = 1.0;
  alpha = -exp(lambda);
  beta = (exp(lambda) - 1.0)/(lambda);
  west_bndy_value = -1.0;
  gamma = alpha;
  delta = beta;
  east_bndy_value = 0.0;

  tools = new MTK_ToolManager();

  //! 1. Discretize: Create the staggered grid.
  space = new MTK_LR1DUniformStaggeredGrid(west_bndy, east_bndy,
                                           number_of_cell_centers);
  cout << "Defined a LR 1D uniform staggered grid:" << endl << "\t";
  space->Print1DArray();

  //! 2. Initialize mimetic operators.
  dir_comp = new MTK_1DCGMDirichletOperator(desired_order,
                                            number_of_cell_centers,
                                            MTK_DENSE);
  dir_comp->Print(MTK_DENSE);
  neu_comp = new MTK_1DCGMNeumannOperator(desired_order,
                                          number_of_cell_centers,
                                          MTK_DENSE);
  neu_comp->Print(MTK_DENSE);
  grad = new MTK_1DCGMGradient(desired_order,
                               number_of_cell_centers,
                               MTK_DENSE);
  grad->Print(MTK_DENSE);
  lap = new MTK_1DCGMLaplacian(desired_order,
                               number_of_cell_centers,
                               MTK_DENSE);
  lap->Print(MTK_DENSE);

  //! 3. Configure the problem.
  problem_to_solve = new MTK_Steady1DProblem(alpha, beta, west_bndy_value,
                                             gamma, delta, east_bndy_value);
  problem_to_solve->MTKPrint();

  //! 5. Create an stencil stencil matrix based on the operators.
  stencil_matrix = new MTK_1DArrayDenseStencilMatrix(alpha, dir_comp,
                                                     beta, neu_comp,
                                                     grad, lap);

  stencil_matrix->Print(MTK_DENSE);

  //! 6. Discretize the source term.
  source = new MTK_LR1DUniformSourceGrid(space, west_bndy_value,
                                         east_bndy_value, source_term);
  cout << "Defined a LR 1D uniform source grid:" << endl << "\t";
  source->Print1DArray();

  //! 7. Create and solve system of equations using stencil matrix and source.
  MTK_System system(stencil_matrix, source);
  system.Print();
  system.Solve();

  //! 8. Create a grid holding the known analytic solution.
  solution_grid = new MTK_LR1DUniformKnownSolutionHolderGrid(space,
                                                             &known_solution);
  cout << "Solution holder grid:" << endl;
  solution_grid->Print1DArray();

  //! 9. Visualize the known solution but DO NOT save the generated plot:
  plotter = new MTK_1DPlotter(solution_grid->independent(),
                              solution_grid->dependent(),
                              solution_grid->number_of_nodes());
  save_plot = false;
  props = new MTK_1DPlotProperties("x", "p(x)",
                                   "Known analytic solution (5 cells)",
                                   "lines", "red");
  plotter->set_plot_properties(props, save_plot, false, MTK_PNG);
  plotter->See();

  //! 10. Visualize computed solution and save the generated plot:
  plotter = new MTK_1DPlotter(solution_grid->independent(),
                              system.Solution(),
                              solution_grid->number_of_nodes());
  save_plot = false;
  props =
    new MTK_1DPlotProperties("x", "p(x)",
                             "Computed sol.: 2nd order CGM operators (5 cells)",
                             "linespoints", "blue");
  plotter->set_plot_properties(props, save_plot, false, MTK_PNG);
  plotter->See();

  //! 11. Compute norm of the difference between known and computed solutions.
  norm_diff_sol =
    tools->RelativeNorm2Difference(solution_grid->dependent(),
                                   system.Solution(),
                                   solution_grid->number_of_nodes());
  cout << "Relative 2-norm of the difference = " << norm_diff_sol << endl;
  return EXIT_SUCCESS;
}
