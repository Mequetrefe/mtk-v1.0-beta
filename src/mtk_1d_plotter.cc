/*!
 * \file mtk_1d_plotter.cc
 *
 * \brief Creates a default instance of the MTK_1DPlotter class.
 *
 * This file contains the definition of a constructor procedure for the
 * MTK_1DPlotter class. Specifically, it constructs the default instance.
 *
 * \date: Sunday, August 14 2011 11:42 PM
 * \version: 2011-14-08-01.
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

#include <iostream>
#include <fstream>

#include "mtk_1d_node.h"
#include "mtk_1d_plotter.h"
#include "mtk_1d_plot_properties.h"

using namespace std;

/*! Constructs a default MTK_1DPlotter. */

/*! Constructs a default MTK_1DPlotter.
 */
MTK_1DPlotter::MTK_1DPlotter(void):
  independent_variable_(NULL),
  dependent_variable_(NULL),
  number_of_nodes_(0),
  save_image_(false),
  plot_properties_(NULL),
  export_type_(MTK_PNG) {

}

/*! Constructs MTK_1DPlotter from a solution information. */

/*! Constructs MTK_1DPlotter from a solution information.
 * \param **independent_variable INPUT: Pointer to the beginning of the 1D array
 * containing the discretized independent variable of interest.
 * \param **dependent_variable INPUT: Pointer to the beginning of the 1D array
 * containing the attained solution.
 * \param number_of_nodes INPUT: Number of elements in both arrays.
 * \todo Thursday, March 15 2012 11:07 PM: Parametrize the data file creation.
 */
MTK_1DPlotter::MTK_1DPlotter( MTK_1DNode **independent_variable,
                                        MTK_1DNode **dependent_variable,
                                        int number_of_nodes) {

  int ii;
  ofstream file;

  // Part 1. Create the data file!
  this->independent_variable_ = independent_variable;
  this->dependent_variable_ = dependent_variable;
  cout << "Creating data file..." << endl;
  file.open ("solution.dat");
  file << "# Known solution.\n";
  file << "# x f(x)\n";
  for (ii = 0; ii < number_of_nodes; ii++) {
    file  << (this->independent_variable_[ii])->value() << " "
          << (this->dependent_variable_[ii])->value() << endl;
  }
  file << endl;
  file.close();
  cout << "Data file created: solution.dat" << endl;
}

/*! Constructs MTK_1DPlotter from a solution information. */

/*! Constructs MTK_1DPlotter from a solution information.
 * \param **independent_variable INPUT: Pointer to the beginning of the 1D array
 * containing the discretized independent variable of interest.
 * \param **dependent_variable INPUT: Pointer to the beginning of the 1D array
 * containing the attained solution.
 * \param number_of_nodes INPUT: Number of elements in both arrays.
 * \todo Thursday, March 15 2012 11:07 PM: Parametrize the data file creation.
 */
MTK_1DPlotter::MTK_1DPlotter( double *independent_variable,
                                        double *dependent_variable,
                                        int number_of_nodes) {

  int ii;
  ofstream file;

  // Part 1. Create the data file!
  cout << "Creating data file..." << endl;
  file.open ("solution.dat");
  file << "# Known solution.\n";
  file << "# x f(x)\n";
  for (ii = 0; ii < number_of_nodes; ii++) {
    file  << independent_variable[ii] << " "
          << dependent_variable[ii] << endl;
  }
  file << endl;
  file.close();
  cout << "Data file created: solution.dat" << endl;
}

/*! Sets the properties up for a MTK_1DPlotter. */

/*! Sets the properties up for a MTK_1DPlotter AND interfaces with the user
 * regarding the storage options.
 * \param plot_properties INPUT:Contains the actual properties.
 * \param save_image INPUT: Should the image be stored?
 * \param export_type INPUT: IF stored, which format should I use?
 * \todo Thursday, March 15 2012 11:07 PM: Make this actually use the properties
 * class.
 * \todo Monday, September 03, 2012: Select from the export type!
 */
void MTK_1DPlotter::set_plot_properties(MTK_1DPlotProperties *plot_properties,
                                        bool save_image,
                                        bool loglog,
                                        MTK_PlotExport export_type) {

  ofstream file;
  bool legend;

  plot_properties_ = plot_properties;
  legend = true;

  if (plot_properties_->default_ == false) {

    cout << "Creating gnuplot script file... data file created: plot.gpl" << endl;
    this->save_image_ = save_image;
    this->export_type_ = export_type;
    //! \todo: Select from the export type!
    file.open ("plot.gpl");

    if (!this->save_image_) {
      file << "set terminal png notransparent truecolor enhanced size 1024,600" << endl;
      file << "set output \'| display png:-\'" << endl;
    } else {
      file << "set terminal png transparent truecolor enhanced" << endl;
      file << "set output \'plot.png\'" << endl;
    }
    file << "set grid" << endl;
    if (loglog)
      file << "set log xy" << endl;
    file << "set xlabel \" " << plot_properties_->xlabel_ << "\"" << endl;
    file << "set ylabel \"" << plot_properties_->ylabel_ << "\"" << endl;
    if (legend)
      file << "set key center outside bottom title \'Legend\' box 3" << endl;
    file << "set title \"" << plot_properties_->title_ << "\"" << endl;
    file << "plot \'solution.dat\' lt rgb \"" << plot_properties_->color_ << "\" with " << plot_properties_->style_ << endl;
    file.close();

  } else {

    cout << "Creating gnuplot script file... data file created: plot.gpl" << endl;
    this->save_image_ = save_image;
    this->export_type_ = export_type;
    //! \todo: Select from the export type!
    file.open ("plot.gpl");
    if (!this->save_image_) {
      file << "set terminal png notransparent truecolor enhanced size 1024,600" << endl;
      file << "set output \'| display png:-\'" << endl;
    } else {
      file << "set terminal png transparent truecolor enhanced" << endl;
      file << "set output \'plot.png\'" << endl;
    }
    file << "set grid" << endl;
    if (loglog)
      file << "set log xy" << endl;
    file << "set xlabel \" " << plot_properties_->xlabel_ << "\"" << endl;
    file << "set ylabel \"" << plot_properties_->ylabel_ << "\"" << endl;
    file << "set key right outside Left title \'Legend\' box 3" << endl;
    file << "set title \"" << plot_properties_->title_ << "\"" << endl;
    file << "plot \'solution.dat\' with " << plot_properties_->style_ << endl;
    file.close();
  }
}

/*! System call to visualize the MTK_1DPlotter. */

/*! System call to visualize the MTK_1DPlotter.
 */
void MTK_1DPlotter::See(void) {

  system("gnuplot \'plot.gpl\'");
}
