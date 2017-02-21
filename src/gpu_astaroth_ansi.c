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

extern REAL cdata_mp_omega_;
extern REAL cdata_mp_theta_;
extern REAL viscosity_mp_nu_;
static FINT nx_, ny_, nz_;

extern REAL cdata_mp_dx_;
extern REAL cdata_mp_dy_;
extern REAL cdata_mp_dz_;
/* ---------------------------------------------------------------------- */
void FTNIZE(initialize_gpu_c)
     (FINT *nx, FINT *ny, FINT *nz, REAL *x, REAL *y, REAL *z )
/* Initializes GPU.
*/
{
  /*
  printf("omega = %e\n", cdata_mp_omega_);
  printf("nu = %e\n", viscosity_mp_nu_);
  printf("nx = %d\n", *nx);
  printf("ny = %d\n", *ny);
  printf("nz = %d\n", *nz);
  */
  nx_=*nx;
  ny_=*ny;
  nz_=*nz;

  printf("xmin = %e\n", x[4]);
  printf("xmax = %e\n", x[*nx-1+3]);
  printf("ymin = %e\n", y[4]);
  printf("ymax = %e\n", y[*ny-1+3]);
  printf("zmin = %e\n", z[4]);
  printf("zmax = %e\n", z[*nz-1+3]);
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
/* Communication between CPU and GPU; calculation of rhs of momentum eq., duu,
   and of continuity eq., dlnrho, by GPU kernels.
*/
{
  printf("nx_ = %d\n", nx_);
  printf("No GPU implementation yet");
// add calls to ASTAROTH routines here
}
/* ---------------------------------------------------------------------- */
