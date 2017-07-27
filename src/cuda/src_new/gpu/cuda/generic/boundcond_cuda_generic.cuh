#pragma once
#include "common/config.h"

void init_boundcond_cuda_generic(CParamConfig* cparamcfg);

void destroy_boundcond_cuda_generic();

void boundcond_cuda_generic(real* d_lnrho, real* d_uu_x, real* d_uu_y, real* d_uu_z);
