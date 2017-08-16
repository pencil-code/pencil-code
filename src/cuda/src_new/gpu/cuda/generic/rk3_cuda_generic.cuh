#pragma once
#include "common/config.h"


void init_rk3_cuda_generic(CParamConfig* cparams);

void destroy_rk3_cuda_generic();

void rk3_cuda_generic(  real* d_lnrho,      real* d_uu_x,       real* d_uu_y,       real* d_uu_z, 
                        real* d_lnrho_dst,  real* d_uu_x_dst,   real* d_uu_y_dst,   real* d_uu_z_dst,
                        const int step_number, const real dt);

