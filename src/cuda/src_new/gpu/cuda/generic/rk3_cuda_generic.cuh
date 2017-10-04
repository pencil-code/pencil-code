#pragma once
#include "common/config.h"
#include "common/grid.h"


void rk3_cuda_generic(Grid* d_grid, Grid* d_grid_dst, const int step_number, const real dt, CParamConfig* cparams);

