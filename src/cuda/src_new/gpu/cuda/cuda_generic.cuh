#pragma once
#include <cuda_runtime.h>
#
#include "common/datatypes.h"
#include "common/config.h"
#include "common/grid.h"
#include "common/slice.h"
#include "common/forcing.h"

#include "core/concur_cuda_core.cuh"
#include "generic/collectiveops_cuda_generic.cuh"
 
enum streamType {FRONT=0, BACK, BOT, TOP, LEFTRIGHT, NUM_COPY_STREAMS};

/*
*       Struct for storing the necessary information for running any of the
*       single-GPU implementations on some specific GPU.
*
*       Contains f.ex. the "local" grid dimensions (d_cparams, dimensions after
*       decomposing the host grid) and the starting index (start_idx) in the host
*       grid for mapping from local grid to host grid coordinates.
*
*       These contexts are stored in a static global array, gpu_contexts.
*       If we wanted to integrate on, say, device 1 then we would set the current
*       device to (1) and pass the necessary information from gpu_contexts[1] to the
*       integration kernel.
*/

typedef struct{
    vec3i           start_idx;          //The starting coordinates in the host array (node-wide grid)
    CParamConfig    d_cparams;          //Local CParamConfig for the device (GPU-specific grid)
    Grid            d_grid, d_grid_dst; //Core arrays
    Slice           d_slice;            //Slice arrays
    ReductionArray  d_reduct_arr;       //Reduction arrays
    real*           d_halobuffer;       //Buffer used for multi-node halo transfers
    ConcurContext   concur_ctx;         //Device-specific streams and events
    int             d_halobuffer_size;
    int             d_halo_widths_x[3], d_halo_widths_y[3], d_halo_widths_z[3];
    cudaStream_t    d_copy_streams[NUM_COPY_STREAMS];  //Streams for halo copying.
} GPUContext;

//Memory operations
void init_cuda_generic(CParamConfig* cparamconf, RunConfig* runconf);
void initialize_cuda_generic(CParamConfig* cparamconf, RunConfig* runconf, const Grid & h_grid);
void destroy_cuda_generic();
void load_grid_cuda_generic(Grid* h_grid);
void store_grid_cuda_generic(Grid* h_grid);
void load_forcing_params_cuda_generic(ForcingParams* forcing_params);
void update_forcing_coefs_cuda_generic(ForcingParams* forcing_params);
void load_outer_halos_cuda_generic(Grid* g, real* halo);
void store_internal_halos_cuda_generic(Grid* g, real* halo);

//Solver functions
void integrate_cuda_generic(real dt);
void integrate_step_cuda_generic(int isubstep, real dt);
void boundcond_step_cuda_generic();
void exchange_halos_cuda_generic(bool circular);
#ifdef GPU_ASTAROTH
real reduce_cuda_PC(ReductType t, GridType a);
#else
real reduce_cuda_generic(ReductType t, GridType a);
#endif
#
//Misc
void get_slice_cuda_generic(Slice* h_slice);
