#pragma once
#include "common/config.h"


//Memory operations
void init_cuda_generic(CParamConfig* cparamconf, RunConfig* runconf);
void destroy_cuda_generic();
void load_grid_cuda_generic(real* lnrho, real* uu_x, real* uu_y, real* uu_z);
void store_grid_cuda_generic(real* lnrho, real* uu_x, real* uu_y, real* uu_z);

//Solver functions
void integrate_cuda_generic(real dt);
void integrate_step_cuda_generic(int isubstep, real dt);
real reduce_cuda_generic(ReductType t, ArrType a);

//Misc
void get_slice_cuda_generic(real* slice_lnrho, real* slice_uu, real* slice_uu_x, real* slice_uu_y, real* slice_uu_z);
