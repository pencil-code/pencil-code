#pragma once
#include "common/config.h"


void init_collectiveops_cuda_generic(CParamConfig* cparams);
void destroy_collectiveops_cuda_generic();

real get_reduction_cuda_generic(ReductType t, real* d_a, real* d_b = NULL, real* d_c = NULL);
