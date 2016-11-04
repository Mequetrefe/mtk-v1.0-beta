/*!
\file mtk.h

\brief Includes the entire API.

This file contains every required header file, thus containing the entire API.
In this way, client codes only have to instruct #include "mtk.h".

\warning IT IS EXTREMELY IMPORTANT THAT THE HEADERS ARE ADDED IN A SPECIFIC
ORDER; THAT IS, CONSIDERING THE DEPENDENCE BETWEEN THE CLASSES THESE CONTENT!

\date: Wednesday, August 17 2011 07:36 PM

\version: 2011-08-17-01.

\author: Eduardo J. Sanchez: esanchez at mail dot sdsu dot edu

Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
State University (SDSU).

Website for the entire project: http://www.csrc.sdsu.edu/mtk/

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

esanchez at mail dot sdsu dot edu

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*!
\mainpage Introduction

Computational Science has emerged as a cutting edge field when it comes to
advancing the knowledge of diverse sciences. Its interdisciplinary nature has
proven to be an efficient bridge towards the understanding of the diverse
physical phenomena that comprise nature. These physical phenomena are typically
modeled as a set of Partial Differential Equations, which usually correspond
to a particular conservation law. Therefore, numerical methods used to solve
these equations are of vital importance in the current paradigm of
interdisciplinary science.

The Mimetic Methods Toolkit (MTK) is an API which allows for an intuitive
implementation of CGM-based MDMs for the resolution of PDEs, yielding numerical
solutions that guarantee uniform order of accuracy, all along the modeled
physical domain. These numerical solutions ensure the satisfaction of
conservation laws, thus remaining faithful to the underlying physics of the
problem. The existing prototype of MTK has been developed in C++; therefore, it
exploits all the well-known advantages of both object-oriented application
models, and the extensive collection of data structure capabilities of this
language.

\section section_mtk_concerns MTK Concerns

Since collaborative development efforts are definitely important in achieving
the level of generality we intend the library to possess, we have divided the
API's source code according to the designated purpose the classes possess within
the API. These divisions (or concerns) are grouped by layers, and are
hierarchically related by the dependence they have among them.
One concern is said to depend on another one, if the classes the first concern
includes, rely on the classes the second concern includes. The following figure
depicts these concerns as well as the dependence among them:

\image html concerns.PNG "MTK Concerns" width=10cm
\image latex concerns.PNG "MTK Concerns" width=\textwidth
\image rtf concerns.PNG "MTK Concerns" width=\textwidth

\section section_flavors MTK Flavors

The MTK's numerical core is written in C++ but the toolkit is intended to keep
growing. Different and diverse computational needs have been taken into
account thus yielding the following different ''flavors'' or related
libraries to MTK:
-# CMTK: C wrappers collection for MTK; intended for sequential computations.
-# FMTK: Fortran wrappers collection for MTK; intended for sequential
computations.
-# MMTK: MATLAB wrappers collection for MTK; intended for sequential
computations.
-# RMTK: R software environment for statistical computing and graphics
wrappers collection for MTK; intended for sequential computations.
-# PyMTK: Python wrappers collection for MTK; intended for sequential
computations.

\section section_authors Authors and development team

The MTK is an effort mainly carried on by researchers at the Computational
Science Research Center (CSRC) at San Diego State University (SDSU).
Specifically, given the diversity of application examples and also of
theoretical contributors, the development in both, applied and theoretical
aspects are members of:
-# Mimetic Methods Research and Development Group.
-# Carbon Utilization and Sequestration Research and Development Group.

Currently the developers are:
-# Eduardo Sanchez, M.Sc.
-# Jose Castillo, Ph.D.
-# Guillermo Miranda, Ph.D.
-# Christopher Paolini, Ph.D.

\subsection subsection_acknowledgements Acknowledgements and contributions

The authors would like to acknowledge valuable advising, contributions and
feedback, from research personnel at the Computational Science Research Center
at San Diego State University, which were vital to the fruition of this work.
Specifically, our thanks go to (alphabetical order):
-# Mohammad Abouali
-# Dany De Cecchis
-# Julia Rossi
-# Raul Vargas--Navarro

\section section_prog_tools Programming tools

The development of MTK has been made possible through the use of the following
applications:
-# Editor: Kate - KDE Advanced Text Editor v-3.5.1. (c) 2000-2005 The Kate
Authors.
-# Compiler: gcc version 4.4.5 (Ubuntu/Linaro 4.4.4-14ubuntu5).

\section section_license_mod Licensing and modifications

Copyright (C) 2015 Computational Science Research Center (CSRC) at San Diego
State University (SDSU).

Website for the entire project: http://www.csrc.sdsu.edu/mtk/

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

\page page_readme Read me file and installation instructions

<pre>
Mimetic Methods Toolkit (MTK) Version 1.0 (beat) README File

Version: Version 1.0 (Beta).

Contributions and donations: http://www.csrc.sdsu.edu/mtk/

By: Eduardo J. Sanchez - esanchez at mail dot sdsu dot edu
    __________________________________________________________________

1. Description

Thank you for downloading the MTK. The Mimetic Methods Toolkit (MTK) is an API 
that allows for an intuitive implementation of Mimetic Discretization Methods 
for the resolution of Partial Differential Equations, yielding numerical 
solutions that guarantee uniform order of accuracy all along the modeled 
physical domain, while ensuring the satisfaction of conservation laws, thus 
remaining faithful to the underlying physics of the problem. The library is
fully developed in C++11, thus exploiting all the well-known advantages of both,
object-oriented programming and the extensive collection of data structures 
capabilities of this language.
    __________________________________________________________________

2. Dependencies

This README assumes all of these libraries are installed in the following
folder:

$(HOME)/Libraries/

In this version, the MTK optionally uses ATLAS-optimized BLAS routines for the 
internal computation on some of the layers. However, ATLAS requires both BLAS 
and LAPACK in order to create their optimized distributions. Therefore, the 
following dependencies tree arises:

1. LAPACK - Available from: http://www.netlib.org/lapack/
  1.1. BLAS - Available from: http://www.netlib.org/blas/

2. The GLPK - Available from: https://www.gnu.org/software/glpk/

3. (Optional) ATLAS - Available from: http://math-atlas.sourceforge.net/
  2.1. BLAS - Available from: http://www.netlib.org/blas/
  2.2. LAPACK - Available from: http://www.netlib.org/lapack/

4. (Optional) Doxygen - Available from http://www.stack.nl/~dimitri/doxygen/ -
This is not an explicit dependence in order to build the MTK, but it is desired
in case the user wants to generate his/her own documentation (for development
purposes).
    __________________________________________________________________

3. Installation

PART 1. CONFIGURATION OF THE MAKEFILE.

The following steps are required the build and test the MTK. Please use the
provided ''makefile_inc'' file, which should provide a solid template to
start with. The following command provides help on the options for make:

$ make help
-----
Makefile for MTK version 1.0.

Options are:
- make: builds only the library and the examples.
- all: builds the library, the examples and the documentation.
- mtklib: builds the library, i.e. generates the archive files.
- tests: generates the tests.
- examples: generates the examples.
- gendoc: generates the documentation for the library.
- checkheaders: checks syntax of the header files.

- clean: cleans ALL the generated files.
- cleanlib: cleans the generated archive and object files.
- cleantests: cleans the generated tests executables.
- cleanexamples: cleans the generated examples executables.
-----

PART 2. BUILD THE LIBRARY.

$ make

If successful you'll read (before building the examples):

----- Library created! Check:
/home/ejspeiro/Dropbox/Public/mtk/mtk-v1.0/lib/libmtk.a

Examples and tests will also be built.
    __________________________________________________________________

4. Frequently Asked Questions

Q: Why haven't you guys implemented GBS to build the library?
A: I'm on it as we speak! ;)

Q: When will the other flavors be ready?
A: Soon! I'm working on getting help on developing those.

Q: Is there any main reference when it comes to the theory on Mimetic Disc.
Methods?
A: Yes! Check: http://www.csrc.sdsu.edu/mtk/ and look for the bibliography
within the documentation for the API.

Q: Do I need to generate the documentation myself?
A: You can if you want to... but if you DO NOT want to, just go to our website.
    __________________________________________________________________

5. Contact and credits

MTK is a effort mainly carried on by researchers at the Computational Science
Research Center (CSRC) at San Diego State University (SDSU). Specifically, given
the diversity of application examples and also of theoretical contributors, the
development in both, applied and theoretical aspects are members of:

  1. Mimetic Methods Research and Development Group.
  2. Computational Geoscience Research and Development Group.
  3. Ocean Modeling Research and Development Group.

Currently the developers are:

  Jose E. Castillo, Ph.D. - jcastillo at mail dot sdsu dot edu
  Guillermo F. Miranda, Ph.D. - unigrav at hotmail dot com
  Christopher P. Paolini, Ph.D. - paolini at engineering dot sdsu dot edu

  Eduardo J. Sanchez, M.Sc. - esanchez at mail dot sdsu dot edu

Finally, FEEL FREE to contact me with suggestions or corrections:

Eduardo J. Sanchez - esanchez at mail dot sdsu dot edu

Thanks and happy coding!
</pre>

\page page_examples Examples

Examples are given in the <a href="files.html">files list</a> section. They are
provided in the /examples/ folder within the distributed software.

\page page_architectures Test architectures

In this page we intend to make a summary of all of the architectures in where
the MTK has been tested. The MTK is intended to be as portable as possible
throughout architectures. The following architectures have provided flawless
installations of the API and correct execution of the examples:

<pre>
1. Linux 3.2.0-23-generic-pae #36-Ubuntu SMP i386 GNU/Linux
   Intel(R) Pentium(R) M processor 1.73GHz 2048 KB of cache and stepping of 8
   gcc version 4.6.3 (Ubuntu/Linaro 4.6.3-1ubuntu5)
</pre>

Further architectures will be tested!

\page page_ref_the User manual, references and theory

The main source of references for this work can be found in:

http://www.csrc.sdsu.edu/mimetic-book/

However, a .PDF copy of this manual can be found
<a href="../latex/refman.pdf">here</a>.
*/

#ifndef MTK_INCLUDE_MTK_H
#define MTK_INCLUDE_MTK_H

//! CONCERN #1 - Roots. Fundamental constants AND functioning parameters.
#include "mtk_roots.h"

//! CONCERN #2 - Enumerations.
#include "mtk_enums.h"

//! CONCERN #3 - Data structures. Self-explanatory.
#include "mtk_matrix.h"
#include "mtk_dense_matrix.h"
#include "mtk_crs_matrix.h"
#include "mtk_ccs_matrix.h"
#include "mtk_dok_matrix.h"

//! CONCERN #3 - Tools. Execution assistance.
#include "mtk_tool_manager.h"

//! CONCERN #5 - Meshes and grids. Meshes control.
#include "mtk_1d_grid_node.h"
#include "mtk_1d_grid.h"
#include "mtk_grid_manager.h"
#include "mtk_1d_uniform_cartesian_staggered_grid.h"
#include "mtk_1d_nonuniform_curvilinear_staggered_grid.h"

//! CONCERN #6 - Mimetic operators.
#include "mtk_1d_dirichlet.h"
#include "mtk_1d_neumann.h"
#include "mtk_1d_grad.h"
#include "mtk_1d_div.h"
#include "mtk_1d_lap.h"
#include "mtk_2d_grad.h"
#include "mtk_2d_div.h"

//! CONCERN #7 - Numerical methods.
#include "mtk_blas_facade.h"
#include "mtk_lapack_facade.h"
#include "mtk_glpk_facade.h"

//! CONCERN #8 - Input.
#include "mtk_steady_1d_problem.h"

//! CONCERN #9 - Solvers.

//! CONCERN #10 - Output.

#endif
