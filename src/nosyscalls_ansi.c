/*                             nosyscalls_ansi.c
                              -------------------
*/

/* Date:   20-Mar-2010
   Author: Bourdin.KIS (Bourdin@KIS.Uni-Freiburg.de)
   Description:
 Dummy stubs for ANSI C function calls.
 Implement only functions here, which are absolutely necessary for linking.
*/

#include <stdio.h>

#include "headers_c.h"

/* ---------------------------------------------------------------------- */

void FTNIZE(file_size_c)
     (char *filename, FINT *bytes)
/* Determines the size of a file.
   Returns:
   * -1 if retrieving the file size failed
*/
{
  *bytes=-1;
  printf("NoSysCalls: file_size_c is not available\n");
}

/* ---------------------------------------------------------------------- */

