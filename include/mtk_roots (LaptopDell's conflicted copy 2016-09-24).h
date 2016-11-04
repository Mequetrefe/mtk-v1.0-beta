/*!
\file mtk_roots.h

\brief Fundamental definitions to be used across all classes of the MTK.

This file contains the fundamental definitions that classes of the MTK rely on
to be implemented. Examples of these definitions are the definition of
fundamental data types, parameters of functioning handling certain assumptions
dealing with memory allocation, and so forth.

\date: June 18, 2014, 2:29 PM

\version: 2015-02-23.

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

#ifndef MTK_ROOTS_H
#define MTK_ROOTS_H

/*!
\namespace mtk

\brief Mimetic Methods Toolkit namespace.
*/
namespace mtk {

/*!
\typedef MTK_Real

\ingroup c01-roots

\brief Users can simply change this to build a double or single precision MTK.
*/
#ifdef MTK_PRECISION_DOUBLE
typedef double MTK_Real;
#endif
#ifdef MTK_PRECISION_SINGLE
typedef float MTK_Real;
#endif

/*!
\def Max. length for a string.

\ingroup c01-roots
*/
const int MTK_MAX_CHARS =  100;

/*!
\def \def Max. length for a string holding a number:

\ingroup c01-roots
*/
const int MTK_MAX_DIGITS = 10;

/*!
\def Considered tolerance for comparisons all along the MTK (global):

\ingroup c01-roots
*/
const MTK_Real MTK_GlobTol = 1.0e-6;

/*!
\def Considered tolerance for comparisons the MTK matrices:

\ingroup c01-roots
*/
const MTK_Real MTK_MatrixTol = 1.0e-7;

/*!
\def Considered tolerance for higher-order mimetic operators:

\ingroup c01-roots
*/
const MTK_Real MTK_DefaultMimTol = 1.0e-6;

const MTK_Real MTK_DefaultMimeticThreshold = 1.0e-6;

/*!
\def Default order of accuracy for mimetic operators:

\ingroup c01-roots
*/
const int MTK_DefaultOrder = 2;
}

#endif  // End of: MTK_INCLUDE_MTK_ROOTS_H
