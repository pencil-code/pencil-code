#pragma once
#include "common/config.h"


void init_slice_cuda_generic(CParamConfig* cparams, RunConfig* runconf);

void destroy_slice_cuda_generic();

void update_slice_cuda_generic(real* d_lnrho, real* d_uu_x, real* d_uu_y, real* d_uu_z);

void store_slice_cuda_generic(real* slice_lnrho, real* slice_uu, real* slice_uu_x, real* slice_uu_y, real* slice_uu_z);
