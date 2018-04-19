#include "slice_cuda_generic.cuh"
#include "gpu/cuda/core/dconsts_core.cuh"
#include "gpu/cuda/core/errorhandler_cuda.cuh"
#ifdef GPU_ASTAROTH
#include "common/PC_moduleflags.h"
#endif 


void init_slice_cuda_generic(Slice* d_slice, CParamConfig* cparams, RunConfig* run_params)
{
    if (run_params->slice_axis != 'z') CRASH("Slice axis other that z not yet supported!");
    const int slice_size = sizeof(real) * cparams->mx * cparams->my;

    //Allocate device memory
    for (int i=0; i < NUM_SLICES; ++i)
        CUDA_ERRCHK( cudaMalloc(&d_slice->arr[i], slice_size) );
}


void destroy_slice_cuda_generic(Slice* d_slice)
{
    //Deallocate device memory
    for (int i=0; i < NUM_SLICES; ++i)
        CUDA_ERRCHK( cudaFree(d_slice->arr[i]) );
}


//Puts a debug slice in to d_slice_lnrho etc that contains
//the boundary zones (no paddings)
template <char slice_axis>
__global__ void slice_cuda_generic(Slice & slice, Grid & grid)
{
    const int i = threadIdx.x + blockIdx.x*blockDim.x;
    const int j = threadIdx.y + blockIdx.y*blockDim.y;

    if (i >= d_mx || j >= d_my) //If out of bounds
        return;

    const int slice_idx = i + j*d_mx;
    const int grid_idx = slice_idx + (d_mz/2)*d_mxy;//Take the from the middle

#ifdef HYDRO
    if (grid.LNRHO>=0){ 
    	real* d_lnrho = grid.arr[grid.LNRHO];
        real* d_slice_lnrho = slice.arr[SLICE_LNRHO];

    //Load lnrho
        d_slice_lnrho[slice_idx] = d_lnrho[grid_idx];
    }
    if (grid.UUX>=0){ 
      real* d_uux = grid.arr[grid.UUX];
      real* d_uuy = grid.arr[grid.UUY];
      real* d_uuz = grid.arr[grid.UUZ];
      real* d_slice_uu = slice.arr[SLICE_UU];
      real* d_slice_uux = slice.arr[SLICE_UUX];
      real* d_slice_uuy = slice.arr[SLICE_UUY];
      real* d_slice_uuz = slice.arr[SLICE_UUZ];


    //Load uu, uu_x, uu_y, uu_z
      const real uux = d_uux[grid_idx];
      const real uuy = d_uuy[grid_idx];
      const real uuz = d_uuz[grid_idx];
      d_slice_uu[slice_idx] = uux*uux + uuy*uuy + uuz*uuz;
      d_slice_uux[slice_idx] = uux;
      d_slice_uuy[slice_idx] = uuy;
      d_slice_uuz[slice_idx] = uuz;
    }
#endif

#ifdef MAGNETIC
    if (grid.AAX>=0){ 
        slice.arr[SLICE_AAX][slice_idx] = grid.arr[grid.AAX][grid_idx];
        slice.arr[SLICE_AAY][slice_idx] = grid.arr[grid.AAY][grid_idx];
        slice.arr[SLICE_AAZ][slice_idx] = grid.arr[grid.AAZ][grid_idx];
    }
#endif
}


//Slices the assigned axis to d_slice_lnrho etc in device memory
void update_slice_cuda_generic(Slice* d_slice, Grid* d_grid, CParamConfig* cparams, RunConfig* run_params)
{
    if (run_params->slice_axis != 'z') CRASH("Slice axis other that z not yet supported!");
    //CUDA call
    const dim3 threads_per_block(32, 32, 1);
    const dim3 blocks_per_grid((unsigned int)ceil((float) cparams->mx / threads_per_block.x),
                                    (unsigned int)ceil((float) cparams->my / threads_per_block.y),
                                    1);

    slice_cuda_generic<'z'><<<blocks_per_grid, threads_per_block>>>(*d_slice, *d_grid);
    CUDA_ERRCHK_KERNEL();
}
