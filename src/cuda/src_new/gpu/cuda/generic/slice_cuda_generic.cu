#include "slice_cuda_generic.cuh"
#include "gpu/cuda/core/dconsts_core.cuh"
#include "gpu/cuda/core/errorhandler_cuda.cuh"
#include "common/errorhandler.h"


static real *d_slice_lnrho, *d_slice_uu, *d_slice_uu_x, *d_slice_uu_y, *d_slice_uu_z;

static CParamConfig* cparams = NULL;
static RunConfig* run_params = NULL;
static bool is_initialized = false;


void init_slice_cuda_generic(CParamConfig* cparamconf, RunConfig* runconf)
{
    if (is_initialized) CRASH("slice_cuda_generic() already initialized!");
    is_initialized = true;

    //Store config
    cparams = cparamconf;
    run_params = runconf;

    if (run_params->slice_axis != 'z') CRASH("Slice axis other that z not yet supported!");
    const int slice_size = cparams->mx * cparams->my;

    //Allocate device memory
    CUDA_ERRCHK( cudaMalloc(&d_slice_lnrho, sizeof(real)*slice_size) );
    CUDA_ERRCHK( cudaMalloc(&d_slice_uu,    sizeof(real)*slice_size) );
    CUDA_ERRCHK( cudaMalloc(&d_slice_uu_x,  sizeof(real)*slice_size) );
    CUDA_ERRCHK( cudaMalloc(&d_slice_uu_y,  sizeof(real)*slice_size) );
    CUDA_ERRCHK( cudaMalloc(&d_slice_uu_z,  sizeof(real)*slice_size) );
}


void destroy_slice_cuda_generic()
{
    if (!is_initialized) CRASH("slice_cuda_generic() wasn't initialized!");
    is_initialized = false;

    //Get rid of the config ptr
    cparams = NULL;
    run_params = NULL;    

    //Deallocate device memory
    CUDA_ERRCHK( cudaFree(d_slice_lnrho) );
    CUDA_ERRCHK( cudaFree(d_slice_uu)    );
    CUDA_ERRCHK( cudaFree(d_slice_uu_x)  );
    CUDA_ERRCHK( cudaFree(d_slice_uu_y)  );
    CUDA_ERRCHK( cudaFree(d_slice_uu_z)  );
}


//Puts a debug slice in to d_slice_lnrho etc that contains
//the boundary zones (no paddings)
template <char slice_axis>
__global__ void slice_cuda_generic(real* d_lnrho, real* d_uu_x, real* d_uu_y, real* d_uu_z,
                                  real* d_slice_lnrho, real* d_slice_uu, 
                                  real* d_slice_uu_x, real* d_slice_uu_y, real* d_slice_uu_z)
{
    const int i = threadIdx.x + blockIdx.x*blockDim.x;
    const int j = threadIdx.y + blockIdx.y*blockDim.y;

    if (i >= d_mx || j >= d_my) //If out of bounds
        return;

    const int slice_idx = i + j*d_mx;
    const int grid_idx = slice_idx + (d_mz/2)*d_mx*d_my;//Take the from the middle

    //Load lnrho
    d_slice_lnrho[slice_idx] = d_lnrho[grid_idx];

    //Load uu, uu_x, uu_y, uu_z
    const real uu_x = d_uu_x[grid_idx];
    const real uu_y = d_uu_y[grid_idx];
    const real uu_z = d_uu_z[grid_idx];
    d_slice_uu[slice_idx] = uu_x*uu_x + uu_y*uu_y + uu_z*uu_z;
	d_slice_uu_x[slice_idx] = uu_x;
	d_slice_uu_y[slice_idx] = uu_y;
	d_slice_uu_z[slice_idx] = uu_z;
}


//Slices the assigned axis to d_slice_lnrho etc in device memory
void update_slice_cuda_generic(real* d_lnrho, real* d_uu_x, real* d_uu_y, real* d_uu_z)
{
    if (!is_initialized) CRASH("slice_cuda_generic() wasn't initialized!");

    //CUDA call
    const dim3 threads_per_block = {32, 32, 1};
    const dim3 blocks_per_grid   = {(unsigned int)ceil((double) cparams->mx / threads_per_block.x),
                                    (unsigned int)ceil((double) cparams->my / threads_per_block.y),
                                    1};

    slice_cuda_generic<'z'><<<blocks_per_grid, threads_per_block>>>( d_lnrho, d_uu_x, d_uu_y, d_uu_z,
                                                                    d_slice_lnrho, d_slice_uu,
                                                                    d_slice_uu_x, d_slice_uu_y, d_slice_uu_z);
    CUDA_ERRCHK_KERNEL();
}


void store_slice_cuda_generic(real* slice_lnrho, real* slice_uu, real* slice_uu_x, real* slice_uu_y, real* slice_uu_z)
{
    if (run_params->slice_axis != 'z') CRASH("Slice axis other that z not yet supported!");
    const int slice_size = cparams->mx * cparams->my;

    CUDA_ERRCHK( cudaMemcpy(slice_lnrho, d_slice_lnrho, sizeof(real)*slice_size, cudaMemcpyDeviceToHost) );
    CUDA_ERRCHK( cudaMemcpy(slice_uu,    d_slice_uu,    sizeof(real)*slice_size, cudaMemcpyDeviceToHost) );
    CUDA_ERRCHK( cudaMemcpy(slice_uu_x,  d_slice_uu_x,  sizeof(real)*slice_size, cudaMemcpyDeviceToHost) );
    CUDA_ERRCHK( cudaMemcpy(slice_uu_y,  d_slice_uu_y,  sizeof(real)*slice_size, cudaMemcpyDeviceToHost) );
    CUDA_ERRCHK( cudaMemcpy(slice_uu_z,  d_slice_uu_z,  sizeof(real)*slice_size, cudaMemcpyDeviceToHost) );
}




























































