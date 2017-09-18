/*
*   Implementation for the generic cuda solution.
*   Manages multiple GPUs on a single node.
*/
#include "cuda_generic.cuh"

#include "core/cuda_core.cuh"
#include "core/dconsts_core.cuh"
#include "core/errorhandler_cuda.cuh"

#include "generic/collectiveops_cuda_generic.cuh"
#include "generic/rk3_cuda_generic.cuh"
#include "generic/boundcond_cuda_generic.cuh"
#include "generic/slice_cuda_generic.cuh"

//Host configs
static CParamConfig h_cparams;
static RunConfig h_run_params;
static bool is_initialized;

typedef struct {
    vec3i           start_idx;    //The starting coordinates in the host array (node-wide grid)
    CParamConfig    d_cparams;    //Local CParamConfig for the device (GPU-specific grid)
    Grid            d_grid, d_grid_dst; //Core arrays
    Slice           d_slice;            //Slice arrays
    ReductionArray  d_reduct_arr;       //Reduction arrays
} GPUContext;

static GPUContext gpu_contexts[NUM_DEVICES];


void init_cuda_generic(CParamConfig* cparamconf, RunConfig* runconf)
{   
    if (is_initialized) CRASH("cuda_generic already initialized!");
    is_initialized = true;

    //Copy the structs in case the caller deallocates them prematurely
    h_cparams = *cparamconf;
    h_run_params = *runconf;

    #pragma unroll
    for (int device_id=0; device_id < NUM_DEVICES; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];


        //Decompose the problem
        ctx->d_cparams = h_cparams;
        ctx->d_cparams.ny = h_cparams.ny / NUM_DEVICES; //Slice the y axis
        ctx->start_idx = (vec3i){0, device_id * ctx->d_cparams.ny, 0};
        printf("%d and start %d\n", ctx->d_cparams.ny, ctx->start_idx.y);
        ///////////TODO fix below        
        ctx->d_cparams.compute_missing_values(); //Purkka
        ctx->d_cparams.dsx = h_cparams.dsx;
        ctx->d_cparams.dsy = h_cparams.dsy;
        ctx->d_cparams.dsz = h_cparams.dsz;
        ctx->d_cparams.dsmin = h_cparams.dsmin;
        /////////


        load_hydro_dconsts_cuda_core(&ctx->d_cparams, &h_run_params, ctx->start_idx);
        init_grid_cuda_core(&ctx->d_grid, &ctx->d_grid_dst, &ctx->d_cparams);
        init_slice_cuda_generic(&ctx->d_slice, &ctx->d_cparams, &h_run_params);
        init_reduction_array_cuda_generic(&ctx->d_reduct_arr, &ctx->d_cparams);
    }
}


void destroy_cuda_generic()
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");
    is_initialized = false;

    #pragma unroll
    for (int device_id=0; device_id < NUM_DEVICES; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
    
        destroy_slice_cuda_generic(&ctx->d_slice);
        destroy_reduction_array_cuda_generic(&ctx->d_reduct_arr);
        destroy_grid_cuda_core(&ctx->d_grid, &ctx->d_grid_dst);
        cudaDeviceSynchronize(); cudaDeviceReset();
    }
}


void load_grid_cuda_generic(Grid* h_grid)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    //If we wanted to use another layout, we would do it here instead of using the core interface
    #pragma unroll
    for (int device_id=0; device_id < NUM_DEVICES; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
    
        load_grid_cuda_core(&ctx->d_grid, &ctx->d_cparams, &ctx->start_idx, h_grid, &h_cparams); 
    }
}


void store_grid_cuda_generic(Grid* h_grid)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    #pragma unroll
    for (int device_id=0; device_id < NUM_DEVICES; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];

        store_grid_cuda_core(h_grid, &h_cparams, &ctx->d_grid, &ctx->d_cparams, &ctx->start_idx); 
    }
}


inline void swap_ptrs(real** a, real** b)
{
	real* temp = *a;
	*a = *b;
	*b = temp;
}


static inline void swap_grid_ptrs(Grid* d_grid, Grid* d_grid_dst)
{
    for (int i=0; i < NUM_ARRS; ++i) 
        swap_ptrs(&(d_grid->arr[i]), &(d_grid_dst->arr[i]));
}

/*
//
// Boundary conditions used for standalone Astaroth
//
#include "cpu/model/model_boundcond.h"
//Bruteforce the bounds by transferring the whole grid to CPU, 
//solving the boundaries and then transferring the grid back to GPU
static void bruteforce_bounds()
{
    Grid tmp;
    grid_malloc(&tmp, &h_cparams);
    store_grid_cuda_generic(&tmp);

    boundcond(&tmp, &h_cparams);    

    load_grid_cuda_generic(&tmp);   
    grid_free(&tmp);
}*/


void boundcond_step_cuda_generic()
{
    printf("WARNING: boundcond_step_cuda_generic() not supported in this minimal Astaroth build");
    /*
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    if (NUM_DEVICES > 1)
        bruteforce_bounds();
    else
       boundcond_cuda_generic(&gpu_contexts[0].d_grid, &gpu_contexts[0].d_cparams);
    */     
}


void integrate_step_cuda_generic(int isubstep, real dt)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    //For all GPUs in the node
    #pragma unroll
    for (int device_id=0; device_id < NUM_DEVICES; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];

        //Integrate
        rk3_cuda_generic(&ctx->d_grid, &ctx->d_grid_dst, isubstep, dt, &ctx->d_cparams);

        //Swap src and dst device array pointers
        swap_grid_ptrs(&ctx->d_grid, &ctx->d_grid_dst);
    }
}

//Integration function used in standalone Astaroth
void integrate_cuda_generic(real dt)
{
    printf("WARNING: integrate_cuda_generic() not supported in this minimal Astaroth build");
    /*
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    cudaSetDevice(0);
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

    //Step 1
    boundcond_step_cuda_generic();
    integrate_step_cuda_generic(0, dt);

    //Step 2
    boundcond_step_cuda_generic();
    integrate_step_cuda_generic(1, dt);   

    //Step 3
    boundcond_step_cuda_generic();
    integrate_step_cuda_generic(2, dt);

	//TIME END
    #pragma unroll
    for (int device_id=0; device_id < NUM_DEVICES; ++device_id) {
        cudaSetDevice(device_id);
	    cudaDeviceSynchronize();
    }

    cudaSetDevice(0);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	printf("Full RK3 time elapsed: \t%f ms\n", time);
    */
}

#include "utils/utils.h" //For max/min/sum
real reduce_cuda_generic(ReductType t, GridType grid_type)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    real res[NUM_DEVICES];

    #pragma unroll
    for (int device_id=0; device_id < NUM_DEVICES; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];

        if (t == MAX_VEC_UU || t == MIN_VEC_UU || t == RMS_VEC_UU) {
            if (grid_type != NOT_APPLICABLE) {
                printf("Note: other than NOT_APPLICABLE passed to reduce_cuda_generic as ArrType." 
                       "This has no effect when a vector ReductType is selected\n");
            }
            res[device_id] = get_reduction_cuda_generic(&ctx->d_reduct_arr, t, &ctx->d_cparams, 
                                              ctx->d_grid.arr[UUX], ctx->d_grid.arr[UUY], ctx->d_grid.arr[UUZ]);
        } else {
            if (grid_type == NOT_APPLICABLE) CRASH("Invalid GridType in reduce_cuda_generic");
            res[device_id] = get_reduction_cuda_generic(&ctx->d_reduct_arr, t, &ctx->d_cparams, ctx->d_grid.arr[grid_type]);
        }
    }

    //Bruteforce: find max, min or rms from the gpu results
    for (int i=1; i < NUM_DEVICES; ++i) {
        if (t == MAX_VEC_UU || t == MAX_SCAL)
            res[0] = max(res[0], res[i]);
        else if (t == MIN_VEC_UU || t == MIN_SCAL)
            res[0] = min(res[0], res[i]);
        else if (t == RMS_VEC_UU || t == RMS_SCAL || t == RMS_EXP)
            res[0] = sum(res[0], res[i]);
        else
            CRASH("Unexpected ReductType in reduce_cuda_generic()");
    }


    return res[0];   
}


void get_slice_cuda_generic(Slice* h_slice)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    //For each GPU in the node
    #pragma unroll
    for (int device_id=0; device_id < NUM_DEVICES; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
        update_slice_cuda_generic(&ctx->d_slice, &ctx->d_grid, &ctx->d_cparams, &h_run_params);
    }

    //For each GPU in the node
    #pragma unroll
    for (int device_id=0; device_id < NUM_DEVICES; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
        store_slice_cuda_core(h_slice, &h_cparams, &h_run_params, &ctx->d_slice, &ctx->d_cparams, &ctx->start_idx);
    }
}


void load_forcing_params_cuda_generic(ForcingParams* forcing_params)
{
    #pragma unroll
    for (int device_id=0; device_id < NUM_DEVICES; ++device_id) {
        cudaSetDevice(device_id);
        load_forcing_dconsts_cuda_core(forcing_params);
    }
}


























































