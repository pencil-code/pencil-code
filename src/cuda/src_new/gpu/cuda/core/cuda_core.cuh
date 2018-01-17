#pragma once
#include "common/config.h"
#include "common/grid.h"
#include "common/slice.h"
#include "common/forcing.h"

//Init and destroy/////////////////////////////////////////////////////////////
void init_grid_cuda_core(Grid* d_grid, Grid* d_grid_dst, CParamConfig* cparams);
void destroy_grid_cuda_core(Grid* d_grid, Grid* d_grid_dst);

void init_halo_cuda_core(real* d_halo);
void destroy_halo_cuda_core(real* d_halo);

//Load and store///////////////////////////////////////////////////////////////
void load_grid_cuda_core(Grid* d_grid, CParamConfig* d_cparams, vec3i* h_start_idx, 
                         Grid* h_grid, CParamConfig* h_cparams);

void store_grid_cuda_core(Grid* h_grid, CParamConfig* h_cparams, 
                          Grid* d_grid, CParamConfig* d_cparams, vec3i* h_start_idx);

void store_slice_cuda_core(Slice* h_slice, CParamConfig* h_cparams, RunConfig* h_run_params, 
                           Slice* d_slice, CParamConfig* d_cparams, vec3i* h_start_idx);

void load_outer_halo_cuda_core(Grid* d_grid, real* d_halobuffer, CParamConfig* d_cparams,
                               Grid* h_grid, real* h_halobuffer, CParamConfig* h_cparams, 
                               vec3i* h_start_idx);

void store_internal_halo_cuda_core(Grid* h_grid, real* h_halobuffer, CParamConfig* h_cparams, 
                                   vec3i* h_start_idx, 
                                   Grid* d_grid, real* d_halobuffer, CParamConfig* d_cparams);


//Device constant load and store///////////////////////////////////////////////
void load_hydro_dconsts_cuda_core(CParamConfig* cparams, RunConfig* run_params, 
                                  const vec3i start_idx);
void load_forcing_dconsts_cuda_core(ForcingParams* forcing_params);

//Utils////////////////////////////////////////////////////////////////////////
void print_gpu_config_cuda_core();
