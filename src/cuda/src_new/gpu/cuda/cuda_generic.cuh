#pragma once
#include "common/config.h"
#include "common/grid.h"
#include "common/slice.h"
#include "common/forcing.h"


//Memory operations
void init_cuda_generic(CParamConfig* cparamconf, RunConfig* runconf);
void destroy_cuda_generic();
void load_grid_cuda_generic(Grid* h_grid);
void store_grid_cuda_generic(Grid* h_grid);
void load_forcing_params_cuda_generic(ForcingParams* forcing_params);

//Solver functions
void integrate_cuda_generic(real dt);
void integrate_step_cuda_generic(int isubstep, real dt);
void boundcond_step_cuda_generic();
real reduce_cuda_generic(ReductType t, GridType a);

//Misc
void get_slice_cuda_generic(Slice* h_slice);
