#pragma once
#include "common/config.h"
#include "common/grid.h"


void boundcond_cuda_generic(Grid* d_grid, CParamConfig* cparams);

void periodic_xy_boundconds_cuda_generic(Grid* d_grid, CParamConfig* cparams);
