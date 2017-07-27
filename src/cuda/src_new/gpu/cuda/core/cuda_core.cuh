#pragma once
#include "common/config.h"

#ifdef INCLUDED_FROM_CUDA_CORE
    #define EXTERN
#else
    #define EXTERN extern 
#endif

EXTERN real *d_lnrho, *d_uu_x, *d_uu_y, *d_uu_z;
EXTERN real *d_lnrho_dst, *d_uu_x_dst, *d_uu_y_dst, *d_uu_z_dst;

void init_cuda_core(CParamConfig* cparams, RunConfig* runconf);
void destroy_cuda_core();

void load_grid_cuda_core(real* lnrho, real* uu_x, real* uu_y, real* uu_z);
void store_grid_cuda_core(real* lnrho, real* uu_x, real* uu_y, real* uu_z);
