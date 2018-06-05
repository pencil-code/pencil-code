#pragma once
#include "common/config.h"
#include "common/grid.h"

void rk3_entropy_step(const int step_number, const Grid* d_grid, Grid* d_grid_dst, const real dt,
                                   const CParamConfig* cparams,
                                   const cudaStream_t stream);
