#pragma once
#include "common/config.h"
#include "common/grid.h"
#include "common/slice.h"
#include "common/forcing.h"
#include "gpu/cuda/cuda_generic.cuh"

//Init and destroy/////////////////////////////////////////////////////////////
void init_grid_cuda_core(Grid* d_grid, Grid* d_grid_dst, CParamConfig* cparams);
void destroy_grid_cuda_core(Grid* d_grid, Grid* d_grid_dst);
void init_halo_cuda_core(GPUContext & ctx, bool firstGPU, bool lastGPU);
void destroy_halo_cuda_core(GPUContext & ctx);

//Load and store///////////////////////////////////////////////////////////////
void load_grid_cuda_core(Grid* d_grid, CParamConfig* d_cparams, vec3i* h_start_idx, 
                         Grid* h_grid, CParamConfig* h_cparams);

void store_grid_cuda_core(Grid* h_grid, CParamConfig* h_cparams, 
                          Grid* d_grid, CParamConfig* d_cparams, vec3i* h_start_idx);

// Forcing
void load_forcing_dconsts_cuda_core(ForcingParams* forcing_params);
void update_forcing_coefs_cuda_PC(ForcingParams* fp, CParamConfig* cparams, int start_idx);

void store_slice_cuda_core(Slice* h_slice, CParamConfig* h_cparams, RunConfig* h_run_params, 
                           Slice* d_slice, CParamConfig* d_cparams, vec3i* h_start_idx);
#ifdef GPU_ASTAROTH 
void load_outer_halo_cuda_core(const GPUContext & ctx, Grid* h_grid, real* h_halobuffer, bool lfirstGPU, bool llastGPU);
void store_internal_halo_cuda_core(const GPUContext & ctx, Grid* h_grid, real* h_halobuffer, bool lfirstGPU, bool llastGPU);
#else
void load_outer_halo_cuda_core(Grid* d_grid, real* d_halobuffer, CParamConfig* d_cparams,
                               Grid* h_grid, real* h_halobuffer, CParamConfig* h_cparams, 
                               vec3i* h_start_idx);

void store_internal_halo_cuda_core(Grid* h_grid, real* h_halobuffer, CParamConfig* h_cparams, 
                                   vec3i* h_start_idx, 
                                   Grid* d_grid, real* d_halobuffer, CParamConfig* d_cparams);
#endif

//Device constant load and store///////////////////////////////////////////////
void load_hydro_dconsts_cuda_core(CParamConfig* cparams, RunConfig* run_params, 
                                  const vec3i start_idx);

//Utils////////////////////////////////////////////////////////////////////////
void print_gpu_config_cuda_core();
