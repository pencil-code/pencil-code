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

void torch_train_c_api(REAL*); 
void torch_infer_c_api(int *);
void initGPU();
void registerGPU();
void initializeGPU(REAL*, FINT);
void finalizeGPU();
void getFArrayIn(REAL **);
void substepGPU(int );
void sourceFunctionAndOpacity(int);
void copyFarray(REAL*);
void loadFarray();
void reloadConfig();
void updateInConfigArr(int);
int  updateInConfigArrName(char *);
void updateInConfigScal(int,REAL);
int  updateInConfigScalName(char *, REAL);
void testRHS(REAL*,REAL*);
void gpuSetDt();
void random_initial_condition(void);

// for Gnu Compiler
extern char *__cparam_MOD_coornames;
extern REAL __cdata_MOD_y[14];
extern REAL __cdata_MOD_dx, __cdata_MOD_dy, __cdata_MOD_dz;

typedef struct real3{
  REAL x,y,z;
} real3;

typedef struct int3{
  int x,y,z;
} int3;

/* ---------------------------------------------------------------------- */
void FTNIZE(torchtrain_c)(REAL* loss_val)
{
	torch_train_c_api(loss_val);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(torchinfer_c)(int *flag)
{
	torch_infer_c_api(flag);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(initialize_gpu_c)(REAL* f, FINT* comm_fint)
{
// Initializes GPU.  
  /*
  printf("nx = %d\n", *nx);
  printf("ny = %d\n", *ny);
  printf("nz = %d\n", *nz);
  */
  //printf("coornames(1)= %s", __cparam_MOD_coornames[0]);

  //printf("dx = %f\n", __cdata_MOD_dx);
  //printf("dy = %f\n", __cdata_MOD_dy);
  //printf("dz = %f\n", __cdata_MOD_dz);

  initializeGPU(f, *comm_fint);
/*
  printf("xmin = %e\n", x[4]);
  printf("xmax = %e\n", x[nx-1+3]);
  printf("ymin = %e\n", y[4]);
  printf("ymax = %e\n", y[ny-1+3]);
  printf("zmin = %e\n", z[4]);
  printf("zmax = %e\n", z[nz-1+3]);
*/
}
/* ---------------------------------------------------------------------- */
void FTNIZE(init_gpu_c)()
{
// Initializes GPU use.
  initGPU();
}
/* ---------------------------------------------------------------------- */
void FTNIZE(register_gpu_c)()
{
// Allocates memory on GPU according to setup needs.

  registerGPU();
}
/* ---------------------------------------------------------------------- */
void FTNIZE(finalize_gpu_c)()
{
// Frees memory allocated on GPU.

  finalizeGPU();
}
/* ---------------------------------------------------------------------- */
void FTNIZE(get_farray_ptr_gpu_c)(REAL** p_f_in)
{
  getFArrayIn(p_f_in);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(rhs_gpu_c)(FINT *isubstep)

/* Communication between CPU and GPU: copy (outer) halos from CPU to GPU, 
   copy "inner halos" from GPU to CPU; calculation of rhss of momentum eq.
   and of continuity eq. by GPU kernels. Perform the Runge-Kutta substep 
   with number isubstep.
   Value at position ix=1,...,nx, iy=1,...,ny, iz=1,...,nz in the grid
   is found at the position ix-1+nghost + mx*(iy+nghost-1) + mx*my*(iz+nghost-1) in uu or lnrho.
   Here mx=nx+2*nghost etc.

   At beginning of substep: copy (outer) halos from host memory to device memory, that is (Fortran indexing)

   uu_x(1   :nghost,:,:), 
   uu_x(nx+1:mx,    :,:),

   uu_x(nghost+1:nghost+nx,1   :nghost,:), 
   uu_x(nghost+1:nghost+nx,ny+1:my,    :),

   uu_x(nghost+1:nghost+nx,nghost+1:nghost+ny,1   :nghost),
   uu_x(nghost+1:nghost+nx,nghost+1:nghost+ny,nz+1:mz    ).

   At end of substep: copy "inner halos" from device memory to host memory, that is

   uu_x(nghost+1:2*nghost ,nghost+1:nghost+ny,nghost+1:nghost+nz), 
   uu_x(nx+1    :nx+nghost,nghost+1:nghost+ny,nghost+1:nghost+nz),

   uu_x(2*nghost+1:nx,nghost+1:2*nghost ,nghost+1:nghost+nz), 
   uu_x(2*nghost+1:nx,ny+1    :ny+nghost,nghost+1:nghost+nz), 

   uu_x(2*nghost+1:nx,2*nghost+1:ny,nghost+1:2*nghost), 
   uu_x(2*nghost+1:nx,2*nghost+1:ny,nz+1    :nz+nghost) 

   If full=1, however, copy the full arrays.
*/
{
// Performs integration substep on GPU.

  substepGPU(*isubstep);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(copy_farray_c)(REAL* f)
{
// Copies vertex buffers from GPU into f-array on CPU.

  copyFarray(f);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(load_farray_c)()
{
// Copies f-array on CPU to vertex buffers on GPU.

  loadFarray();
}
/* ---------------------------------------------------------------------- */
void FTNIZE(reload_gpu_config_c)()
{
  reloadConfig();
}
/* ---------------------------------------------------------------------- */
void FTNIZE(update_on_gpu_scal_by_ind_c)(int *index, REAL* value)
{
  updateInConfigScal(*index,*value);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(update_on_gpu_arr_by_ind_c)(int *index)
{
  updateInConfigArr(*index);
}
/* ---------------------------------------------------------------------- */
int FTNIZE(update_on_gpu_scal_by_name_c)(char *varname, REAL* value)
{
  return updateInConfigScalName(varname,*value);
}
/* ---------------------------------------------------------------------- */
int FTNIZE(update_on_gpu_arr_by_name_c)(char *varname)
{
  return updateInConfigArrName(varname);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(test_rhs_c)(REAL* f_in, REAL* df_truth)
{
  testRHS(f_in,df_truth);
}
/* ---------------------------------------------------------------------- */
void FTNIZE(gpu_set_dt_c)()
{
  gpuSetDt();
}
/* ---------------------------------------------------------------------- */
void FTNIZE(calcq_gpu_c)(int *idir, int3 *dir, int3 *stop, real3 *unit_vec, int *lperiodic){
 // performs ray integration along direction dir for all possible starting points in subdomain,
 // communication and final correction of Q 

}
/* ---------------------------------------------------------------------- */
void FTNIZE(source_function_and_opacity_gpu_c)(int *inu){
	sourceFunctionAndOpacity(*inu);
}
/* ---------------------------------------------------------------------- */
