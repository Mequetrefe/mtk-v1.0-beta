/*!
  \file mtk_enums.h

  \brief Includes the definition of the enumerated data types.

  Is a convention within MTK that if the conceptualization of the ADT is not
  large enough in terms of its operational behavior, then it should be
  implemented as an enum. In this file, we've gather all the necessary enums.

  \date: Tuesday, March 13 2012 12:02 AM
  \version: 2011-14-08-01.
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

#ifndef MTK_INCLUDE_MTK_ENUMS_H
#define MTK_INCLUDE_MTK_ENUMS_H

/*! \brief Representations of the logic constants.
 *
 * A matrix operator can be represented as a dense or sparse matrix.
 */
enum MTK_Boolean {
  MTK_FALSE,
  MTK_TRUE
};

/*! \brief Representations of the operators.
 *
 * A matrix operator can be represented as a dense or sparse matrix.
 */
enum MTK_MatrixRepresentation {
  MTK_DENSE,
  MTK_SPARSE
};

/*! \brief Uniformity of the meshes.
 *
 * A mesh can be uniform or non-uniform as described in different works.
 */
enum MTK_Uniformity {
  MTK_NON_UNIFORM,
  MTK_UNIFORM
};

// Binding.

/*! \brief Methods to solve Ordinary Differential Equation's Init. Val. Probs.
 *
 * Collection of methods to solve Initial Value Problems in Ordinary
 * Differential Equations.
 */
enum MTK_ODEMethod {
  MTK_FORWARD_EULER,
  MTK_BACKWARD_EULER,
  MTK_TRAPEZOIDAL,
  MTK_HEUN
};

/*! \brief Plot properties.
 *
 * Collection of properties, the user can specify to tailor the resulting plot.
 */
enum MTK_PlotProperties {
  MTK_DEFAULT,
  MTK_CUSTOM
};

/*! \brief File formats to export generated plots.
 *
 * Collection of file formats (images) to which the plots can be exported.
 */
enum MTK_PlotExport {
  MTK_PNG
};

typedef enum  {
  NullBind,
  GradNaturalGrid,
  DivNaturalGrid,
  NodeBound,
  CenterBound,
  NodeBound_NoBC,
  CenterBound_NoBC

} BindingType1D;

#endif
