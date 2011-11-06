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
  *bytes = -1;
  printf ("NoSysCalls: file_size_c is not available\n");
}

/* ---------------------------------------------------------------------- */

void FTNIZE(get_pid_c)
     (FINT *pid)
/* Determines the PID of the current process.
   Returns:
   * integer containing the PID of the current process
   * -1 if retrieving the PID failed
*/
{
  *pid = -1;
  printf ("NoSysCalls: get_pid_c is not available\n");
}

/* ---------------------------------------------------------------------- */

void FTNIZE(get_env_var_c)
     (char *name, char *value)
/* Gets the content of an environment variable.
   Returns:
   * unchanged input string, if retrieving the environment variable failed
*/
{
  printf ("NoSysCalls: get_env_var_c is not available\n");
}

/* ---------------------------------------------------------------------- */

void FTNIZE(is_nan_c)
     (REAL *value, FINT *result)
/* Determine if value is not a number.
   Returns:
   * 1, if value is not a number
   * -1 on failure
*/
{
  *result = 1;
  printf ("NoSysCalls: is_nan_c is not available\n");
}

/* ---------------------------------------------------------------------- */
void FTNIZE(system_c) (char *command)
/* Date:   04-Nov-2011
   Author: MR (matthias.rheinhardt@helsinki.fi)
   Description: Dummy stub for ANSI C function system.
*/
{
  printf ("NoSysCalls: system_c is not available\n");
}

/* ---------------------------------------------------------------------- */

