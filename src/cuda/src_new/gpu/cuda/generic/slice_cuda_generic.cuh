#pragma once
#include "common/config.h"
#include "common/grid.h"
#include "common/slice.h"

void init_slice_cuda_generic(Slice* d_slice, CParamConfig* cparamconf, RunConfig* run_params);

void destroy_slice_cuda_generic(Slice* d_slice);

void update_slice_cuda_generic(Slice* d_slice, Grid* d_grid, CParamConfig* cparams, RunConfig* run_params);
