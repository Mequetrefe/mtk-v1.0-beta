/*!
  \file mtk_1d_plotter.h

  \brief Includes the definition of the class MTK_1DPlotter.

  This file contains the definition of the class MTK_1DPlotter class. This
  is MTK's internal mechanism to help provide visual output tailored to user's
  need. For this, the class interfaces with <a href="http://www.gnuplot.info/">
  GNUPlot</a>, thus establishing the latter as an optional dependence to the
  toolkit. This class plots only graphics attained from the solution on 1D
  scenarios.

  \date: Sunday, September 09, 2012
  \version: 2011-09-09.
  \author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu
 */
 /*
  Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
  State University (SDSU).

  http:www.csrc.sdsu.edu/mtk/

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

#ifndef MTK_INCLUDE_MTK_1D_PLOTTER_H
#define MTK_INCLUDE_MTK_1D_PLOTTER_H

#include "mtk_constants.h"
#include "mtk_enums.h"

class MTK_1DNode;
class MTK_1DPlotProperties;

/*! \brief Plotting mechanism.

  This class interfaces with <a href="http://www.gnuplot.info/">
  GNUPlot</a> to help produce visual output tailored to user's need.
  \todo Monday, September 03, 2012: Is GNUPLOT really optional?
*/
class MTK_1DPlotter {

public:
  /*! Default constructor.  */
  /*! This procedure CONSTRUCTS a default MTKPlotter. */
  MTK_1DPlotter(void);

  /*! Constructor #2. */
  /*! This procedure CONSTRUCTS a MTKPlotter based on two collections of values;
    namely, a collection of independent values, a collection of dependent
    values and the information of the number of nodes.
  */
  MTK_1DPlotter( MTK_1DNode **independent_variable,
                MTK_1DNode **dependent_variable,
                int number_of_nodes);

  MTK_1DPlotter( double *independent_variable,
                double *dependent_variable,
                int number_of_nodes);

  /*! Default destructor.  */
  /*! This procedure DESTRUCTS an already created MTKPlotter. */
  ~MTK_1DPlotter() {
    delete [] independent_variable_;
    delete [] dependent_variable_;
  }

  /*! Set properties of the 1D plot to construct.  */
  /*! This procedure establishes the properties any plot can have. */
  void set_plot_properties(MTK_1DPlotProperties *plot_properties,
                                        bool save_image,
                                        bool loglog,
                                        MTK_PlotExport export_type);

  /*! See plot.  */
  /*! Once the plot has been successfully constructed, one most explicitly order
    its visualization.
  */
  void See(void);

private:
  MTK_1DNode **independent_variable_;  /*!< Array indep. values. */
  MTK_1DNode **dependent_variable_;    /*!< Array depen. values. */
  int number_of_nodes_;                     /*!< Number of nodes involved. */
  bool save_image_;                         /*!< Save attained plot? */
  MTK_1DPlotProperties *plot_properties_;    /*!< Properties of the plot. */
  MTK_PlotExport export_type_;               /*!< Exportation of the plot. */
};

#endif
