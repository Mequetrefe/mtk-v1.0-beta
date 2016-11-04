/*!
\file mtk_enums.h

\brief Considered enumeration types in the MTK.

Enumeration types are used throughout the MTK to differentiate instances of
derived classes, as well as for mnemonic purposes. In this file, the enumeration
types are listed alphabetically.

\date: June 18, 2014, 2:57 PM

\version: 1.

\author: Eduardo J. Sanchez (ejspeiro) - esanchez at sciences dot sdsu dot edu

\bug No known bugs.
*/

/*
Copyright (C) 2015 Computational Science Research Center (CSRC) at San Diego
State University (SDSU).

Website for the project: http://www.csrc.sdsu.edu/mtk/

All rights reserved.

Redistribution and use in source and binary forms, with or without modification
are permitted provided that the following conditions are met:

-# Redistributions of source code must retain the above copyright notice, this
list of conditions and this disclaimer.
-# Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.
-# Neither the name of the CSRC, SDSU nor the names of its contributors may be
used to endorse or promote products derived from this software without specific
prior written permission.
-# Modifications whether they are partial or complete; whether they are
additions or eliminations should be reported through an email at:

esanchez at sciences dot sdsu dot edu

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef MTK_ENUMS_H
#define MTK_ENUMS_H

namespace mtk {

/*!
\enum MTK_Boolean

\ingroup c02-enums

\brief Representations of the logic constants.
*/
enum MTK_Boolean {
  MTK_FALSE,
  MTK_TRUE
};

/*!
\enum MTK_BoundaryType

\ingroup c02-enums

\brief This describes a collection of supported boundary conditions (BCs).

PDEs can be solved considering one of the different types of BCs we hereby
include. The selection of the adequate type of BC affects the allocation of
memory within the data structures that are utilized to solve for the PDE(s) of
interest.
*/
typedef enum {
  MTK_DIRICHLET,
  MTK_NEUMANN,
  MTK_ROBIN,
  MTK_PERIODIC

} MTK_BoundaryType;

/*!
\enum MTK_MatrixStorageScheme

\ingroup c02-enums

\brief Considered matrix storage schemes to implement sparse matrices.

The considered sparse storage schemes are selected so that these are compatible
with some of the most used mathematical APIs, as follows: MTK_BANDED for
<a href="http://www.netlib.org/lapack/">LAPACK</a> and
<a href="http://www.netlib.org/scalapack/">ScaLAPACK</a>, MTK_{CRS,CCS} for
<a href="http://crd.lbl.gov/~xiaoye/SuperLU/">SuperLU</a>, and MTK_DOK for
<a href="http://graal.ens-lyon.fr/MUMPS/">MUMPS</a>. Finally, the MTK_DENSE
format is intended to be compatible with
<a href="http://www.netlib.org/blas/">CBLAS</a>.
*/
enum MTK_MatrixStorageScheme {
  MTK_NULL_MATRIX_STORAGE_SCHEME,
  MTK_DENSE,
  MTK_BANDED,
  MTK_CRS,
  MTK_CCS,
  MTK_DOK
};

/*!
\enum MTK_IntendedDimensionality

\ingroup c02-enums

\brief This data type describes if a 1-, 2- or 3-D context is being described.

Certain classes within the MTK, inherit from a base class, given the fact that
their computational behavior is different, since they are trying to model
either 1-, 2- or 3-D scenarios. This data type permits a distinction among
these.
*/
enum MTK_IntendedDimensionality {

  MTK_1D,
  MTK_2D,
  MTK_3D
};

typedef enum   {
  NullBind,
  GradNaturalGrid,
  DivNaturalGrid,
  NodeBound,
  CenterBound,
  NodeBound_NoBC,
  CenterBound_NoBC

} BindingType1D;

}

#endif  // End of: MTK_INCLUDE_MTK_ENUMS_H
