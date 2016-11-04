/*!
  \file mtk_1d_plot_properties.h

  \brief Includes the definition of the class MTK_1DPlotProperties.

  This file contains the definition of the class MTK_1DPlotProperties class. This
  help to contains the required properties to build a plot by the
  MTK1DPlotterDouble class.

  \date: Monday, September 03, 2012
  \version: 2011-09-03.
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

#ifndef MTK_INCLUDE_MTK_1D_PLOT_PROPERTIES_H
#define MTK_INCLUDE_MTK_1D_PLOT_PROPERTIES_H

#include <string>
using std::string;

#include "mtk_constants.h"
#include "mtk_enums.h"

class MTK_1DNode;

/*! \brief Plotting properties holding mechanism.

  This class interfaces with the MTK1DPlotterDouble class.
  \todo Monday, September 03, 2012: Member variables should really be private.
*/
class MTK_1DPlotProperties {

public:
  /*! Default constructor.  */
  /*! This procedure CONSTRUCTS a default MTKPlotProperties object. */
  MTK_1DPlotProperties(void):
    default_(true),
    xlabel_("x"),
    ylabel_("f(x)"),
    title_("Graph of f"),
    style_("lines"),
    color_("red"){}

  /*! Constructor #2. */
  /*! This procedure CONSTRUCTS a MTKPlotProperties based on the collection of
   * the provided properties.
  */
  MTK_1DPlotProperties(string xlabel,  string ylabel, string title,
                      string style,   string color) {

    xlabel_ = xlabel;
    ylabel_ = ylabel;
    title_ = title;
    style_ = style;
    color_ = color;
  }

  /*! Default destructor.  */
  /*! This procedure DESTRUCTS an already created MTKPlotProperties. */
  ~MTK_1DPlotProperties() {}

  /*! See plot.  */
  /*! Once the plot has been successfully constructed, one most explicitly order
    its visualization.
  */
  void See(void);

  bool default_;    /*!< Are these the default properties? */
  string xlabel_;   /*!< Label for the x-axis (independent variable). */
  string ylabel_;   /*!< Label for the y-axis (dependent variable). */
  string title_;    /*!< Title of the produced plot. */
  string style_;    /*!< Style for the line of the curve. */
  string color_;    /*!< Color of the line of the plot. */
};

#endif
