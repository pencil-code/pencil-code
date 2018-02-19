#pragma once
#include "common/config.h"
#include "common/grid.h"


void rk3_cuda_generic(Grid* d_grid, Grid* d_grid_dst, 
                      const int step_number, const real dt, 
                      CParamConfig* cparams, cudaStream_t stream=0);

void rk3_inner_cuda_generic(Grid* d_grid, Grid* d_grid_dst, 
                      const int step_number, const real dt, 
                      CParamConfig* cparams, cudaStream_t hydro_stream=0, cudaStream_t induct_stream=0);

void rk3_outer_cuda_generic(Grid* d_grid, Grid* d_grid_dst, 
                      const int step_number, const real dt, 
                      CParamConfig* cparams, cudaStream_t stream=0);
