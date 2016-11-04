/*!
  \file mtk_constants.h

  \brief Includes the definition of the used constant throughout the entire
  library.

  This file contains the definition of the most important constant used to
  define the static data structures along the whole library.

  \date: Sunday, August 14 2011 11:42 PM
  \version: 2011-14-08-01.
  \author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu

  \todo: Monday, March 12 2012 10:47 PM: Figure out a method to store the
  content of the __FILE__ macro.
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

#ifndef MTK_INCLUDE_MTK_CONSTANTS_H
#define MTK_INCLUDE_MTK_CONSTANTS_H

//! \def Max. length for a string:
const int MTK_MAX_CHARS =  100;
//! \def Max. length for a string holding a number:
const int MTK_MAX_DIGITS = 10;
//! \def Prints this line number:
const int MTK_THIS_LINE = __LINE__;
//! \def Prints the path to this file:
// const char *MTK_THIS_FILE = __FILE__;
//! \todo: Figure out a method to store the content of the __FILE__ macro.

#endif
