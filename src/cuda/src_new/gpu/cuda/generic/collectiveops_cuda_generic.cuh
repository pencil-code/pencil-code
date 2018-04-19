#pragma once
#include "common/config.h"

struct ReductionArray{
    real* d_vec_res;
    real* d_partial_result;
};

void init_reduction_array_cuda_generic(ReductionArray* reduct_arr, CParamConfig* cparams);
void destroy_reduction_array_cuda_generic(ReductionArray* reduct_arr);

real get_reduction_cuda_generic(ReductionArray* reduct_arr, ReductType t, CParamConfig* cparams, 
                                real* d_a, real* d_b = NULL, real* d_c = NULL);
