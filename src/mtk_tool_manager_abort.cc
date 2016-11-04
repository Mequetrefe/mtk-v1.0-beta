/*!
  \file mtk_tool_manager_abort.cc
  \date: Monday, August 15 2011 12:00 AM
  \version: 2011-08-14-01.
  \brief Aborts the execution.

  Aborts the execution stating the reason.

  \author: Eduardo J. Sanchez: esanchez@sciences.sdsu.edu

  Copyright (C) 2011 Computational Science Research Center (CSRC) at San Diego
  State University (SDSU).

  http:www.csrc.sdsu.edu/mmtk/

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

#include <cstdio>
#include <cstring>

#include "mtk_tool_manager.h"

/*! Prints a message explaining the abortion of the execution and aborts. */

/*! This procedure prints a message with is comprised of the line and the files
 * in which an error occurred thus leading to the abortion of the execution.
 * \param *error INPUT: Contains the original error message.
 * \param line INPUT: Contains the line where the error occurred.
 * \param *file INPUT: Contains the file which code led to the abortion.
 * \todo Tuesday, August 16 2011 01:50 PM: Analyze the use of signals in
 * signals.h.
 * \todo Tuesday, August 16 2011 01:50 PM: Fix name: CMTKToolManagerAbort (...)
 */
void MTKToolManager::MTKAbort(char *error, int line, char *file) {

  char explanation[MTK_MAX_CHARS]; /* Actual string to be printed. */
  char buffer[MTK_MAX_DIGITS];     /* Helps to make a string from the line. */

  /* Create the explanation: */
  strcpy(explanation, "ERROR: ");
  strcat(explanation, error);
  strcat(explanation, " in ");
  strcat(explanation, file);
  strcat(explanation, " at line ");
  snprintf(buffer, sizeof(buffer), "%d", line);
  strcat(explanation, buffer);
  strcat(explanation, ".\n\0");
  fprintf(stderr, "%s", explanation);
}
