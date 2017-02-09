/*                             gpu_astaroth_ansi.c
                              --------------------
*/

/* Date:   8-Feb-2017
   Author: M. Rheinhardt
   Description:
 ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
*/
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <dlfcn.h>

#include "headers_c.h"

/* ---------------------------------------------------------------------- */
void FTNIZE(initialize_gpu_c)
     (FINT *par)
/* Initializes GPU.
*/
{
}
/* ---------------------------------------------------------------------- */
void FTNIZE(finalize_gpu_c)
     (FINT *par)
/* Frees memory allocated on GPU.
*/
{
}
/* ---------------------------------------------------------------------- */
void FTNIZE(rhs_gpu_c)
     (REAL *uu, REAL *lnrho, REAL *duu, REAL *dlnrho)
/* Communication between CPU and GPU; calculation of rhs by GPU kernels.
*/
{
}
/* ---------------------------------------------------------------------- */
