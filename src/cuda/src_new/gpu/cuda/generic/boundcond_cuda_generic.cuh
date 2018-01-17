#pragma once
#include "common/config.h"
#include "common/grid.h"


void boundcond_cuda_generic(Grid* d_grid, CParamConfig* cparams, cudaStream_t stream=0);

void periodic_xy_boundconds_cuda_generic(Grid* d_grid, CParamConfig* cparams, cudaStream_t stream=0);
