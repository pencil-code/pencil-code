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
    real*           d_halobuffer;       //Buffer used for multi-node halo transfers
} GPUContext;

static GPUContext* gpu_contexts;
static int num_devices = -1;

void init_cuda_generic(CParamConfig* cparamconf, RunConfig* runconf)
{   
    if (is_initialized) CRASH("cuda_generic already initialized!");
    is_initialized = true;

    cudaGetDeviceCount(&num_devices);
    gpu_contexts = (GPUContext*) malloc(sizeof(GPUContext)*num_devices);
    printf("Using %d devices\n", num_devices);

    //Copy the structs in case the caller deallocates them prematurely
    h_cparams = *cparamconf;
    h_run_params = *runconf;

    //#pragma unroll
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];

        //Check p2p availability
        for (int peer=0; peer < num_devices; ++peer) {
            if (device_id == peer)
                continue;
            int can_access = 0;
            cudaDeviceCanAccessPeer(&can_access, device_id, peer);
            printf("%d can access peer %d? %d\n", device_id, peer, can_access);
        }
        const int peer_id = (device_id+1) % num_devices;
        if (device_id != peer_id) {
            printf("Enabling peer access between %d and %d\n", device_id, peer_id);
            cudaDeviceEnablePeerAccess(peer_id, 0);//Note: last parameter here is "flags, reserved for future use and must be set to 0"
        }

        //Decompose the problem
        ctx->d_cparams = h_cparams;
        ctx->d_cparams.nz = h_cparams.nz / num_devices; //Slice the z axis
        ctx->d_cparams.compute_missing_values(); //Purkka
        ctx->d_cparams.dsx = h_cparams.dsx;
        ctx->d_cparams.dsy = h_cparams.dsy;
        ctx->d_cparams.dsz = h_cparams.dsz;
        ctx->d_cparams.dsmin = h_cparams.dsmin;
        ctx->start_idx = (vec3i){0, 0, device_id * ctx->d_cparams.nz};
        printf("%d and start %d\n", ctx->d_cparams.nz, ctx->start_idx.z);

        //Allocate and init memory on the GPU
        load_hydro_dconsts_cuda_core(&ctx->d_cparams, &h_run_params, ctx->start_idx);
        init_grid_cuda_core(&ctx->d_grid, &ctx->d_grid_dst, &ctx->d_cparams);
        init_slice_cuda_generic(&ctx->d_slice, &ctx->d_cparams, &h_run_params);
        init_reduction_array_cuda_generic(&ctx->d_reduct_arr, &ctx->d_cparams);
        init_halo_cuda_core(ctx->d_halobuffer); //Note: Called even without multi-node
    }
}


void destroy_cuda_generic()
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");
    is_initialized = false;

    //#pragma unroll
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
    
        destroy_slice_cuda_generic(&ctx->d_slice);
        destroy_reduction_array_cuda_generic(&ctx->d_reduct_arr);
        destroy_grid_cuda_core(&ctx->d_grid, &ctx->d_grid_dst);
        destroy_halo_cuda_core(ctx->d_halobuffer);
        cudaDeviceSynchronize(); cudaDeviceReset();
    }
    free(gpu_contexts);
}


void load_grid_cuda_generic(Grid* h_grid)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    //If we wanted to use another layout, we would do it here instead of using the core interface
    //#pragma unroll
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
    
        load_grid_cuda_core(&ctx->d_grid, &ctx->d_cparams, &ctx->start_idx, h_grid, &h_cparams); 
    }
}


void store_grid_cuda_generic(Grid* h_grid)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    //#pragma unroll
    for (int device_id=0; device_id < num_devices; ++device_id) {
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
}

static const int top_device_id = num_devices-1;
static const int bottom_device_id = 0;
*/
/*
*   Terminology for multi-GPU halo transfers (subject to change):
*   (1) Shared outer halo: Sides, edges and corners on top and bottom of the whole grid 
*   (2) Local inner halo: Parts of the edges and sides, for which periodic boundary 
*        conditions can be applied without inter-device communication (i.e. the left and right 
*        edges/sides of the grid held by device)
*   (3) Shared inner halo: the overlapping area of the halos of two neighbouring devices.
*
*   Dependencies: (2) must be completed before starting (3), (1) can be done whenever. All (1), 
*   (2) and (3) must be completed before integration.
*
*
*   Cheatsheet for mapping the boundary zone indices:
    ////Shared outer halo////////////////////////////////////////////
    i.e side (x and z have 1-to-1 map):
    BOTTOM:    nymin -> nymax
    TOP:       ny -> 0

    edge (z has 1-to-1 map):
    BOTTOM:     nxmin, nymin -> nxmax, nymax
    BOTTOM:     nx,    nymin -> 0,     nymax    
    TOP:        nxmin, ny    -> nxmax, 0
    TOP:        nx,    ny    -> 0,     0

    corner
    BOTTOM:     nxmin, nymin, nzmin -> nxmax, nymax, nzmax
    BOTTOM:     nx,    nymin, nzmin -> 0,     nymax, nzmax
    BOTTOM:     nx,    nymin, nz    -> 0,     nymax, 0
    BOTTOM:     nxmin, nymin, nz    -> nxmax, nymax, 0

    BOTTOM:     nxmin, ny, nzmin -> nxmax, 0, nzmax
    BOTTOM:     nx,    ny, nzmin -> 0,     0, nzmax
    BOTTOM:     nx,    ny, nz    -> 0,     0, 0
    BOTTOM:     nxmin, ny, nz    -> nxmax, 0, 0


    ////Local outer halo////////////////////////////////////////////
    side (y and z have 1-to-1 map):
    nxmin -> nxmax
    nx    -> 0

    edge (y has 1-to-1 map):
    nxmin, nzmin -> nxmax, nzmax
    nx,    nzmin -> 0,     nzmax
    nxmin, nz    -> nxmax, 0
    nx,    nz    -> 0,     0


    ////Shared inner halo////////////////////////////////////////////
    halo (x and z have 1-to-1 map, mx*mz-wide transfers)"
    ny -> 0 (lower->upper)
    nymin -> nymax  (Upper->lower)

    ////////////////////////////////////////////////////////////////
*/

static void local_boundconds_cuda_generic()
{
    //#pragma unroll
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
        periodic_xy_boundconds_cuda_generic(&ctx->d_grid, &ctx->d_cparams);    
    }    
}

//Exchange halos between neighbouring blocks 
static void exchange_halos_cuda_generic()
{
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];

        const int peer_id = (device_id+1) % num_devices;

        const size_t slab_size = ctx->d_cparams.mx * ctx->d_cparams.my;
        const size_t transfer_size_bytes = BOUND_SIZE * slab_size * sizeof(real);

        const size_t z_src0 = ctx->d_cparams.nz * slab_size;
        const size_t z_dst0 = 0; 
        const size_t z_src1 = BOUND_SIZE * slab_size;
        const size_t z_dst1 = (ctx->d_cparams.nz + BOUND_SIZE) * slab_size;
        for (int w=0; w < NUM_ARRS; ++w) {
            CUDA_ERRCHK( cudaMemcpyPeerAsync(&gpu_contexts[peer_id].d_grid.arr[w][z_dst0], peer_id, &ctx->d_grid.arr[w][z_src0], device_id, transfer_size_bytes) );
            CUDA_ERRCHK( cudaMemcpyPeerAsync(&ctx->d_grid.arr[w][z_dst1], device_id, &gpu_contexts[peer_id].d_grid.arr[w][z_src1], peer_id, transfer_size_bytes) );
        }
    }
}

//Boundary conditions applied on the ghost zone surrounding the 
//computational domain
static void boundcond_multigpu_cuda_generic()
{

    //bruteforce_bounds();
    local_boundconds_cuda_generic();
    exchange_halos_cuda_generic();
}


void boundcond_step_cuda_generic()
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    if (num_devices > 1) {
        boundcond_multigpu_cuda_generic();
    } else {
       boundcond_cuda_generic(&gpu_contexts[0].d_grid, &gpu_contexts[0].d_cparams);    
    } 
}


void integrate_step_cuda_generic(int isubstep, real dt)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    //For all GPUs in the node
    //#pragma unroll
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];

        //Integrate
        rk3_cuda_generic(&ctx->d_grid, &ctx->d_grid_dst, isubstep, dt, &ctx->d_cparams);

        //Swap src and dst device array pointers
        swap_grid_ptrs(&ctx->d_grid, &ctx->d_grid_dst);
    }
}


void integrate_cuda_generic(real dt)
{
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
    //#pragma unroll
    for (int device_id=0; device_id < num_devices; ++device_id) {
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
}

#include "utils/utils.h" //For max/min/sum
real reduce_cuda_generic(ReductType t, GridType grid_type)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    real* res = (real*) malloc(sizeof(real)*num_devices);

    //#pragma unroll
    for (int device_id=0; device_id < num_devices; ++device_id) {
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
    for (int i=1; i < num_devices; ++i) {
        if (t == MAX_VEC_UU || t == MAX_SCAL)
            res[0] = max(res[0], res[i]);
        else if (t == MIN_VEC_UU || t == MIN_SCAL)
            res[0] = min(res[0], res[i]);
        else if (t == RMS_VEC_UU || t == RMS_SCAL || t == RMS_EXP)
            res[0] = sum(res[0], res[i]);
        else
            CRASH("Unexpected ReductType in reduce_cuda_generic()");
    }

    if (t == RMS_VEC_UU || t == RMS_SCAL || t == RMS_EXP)
        res[0] = sqrt(res[0] * 1.0 / (h_cparams.nx*h_cparams.ny*h_cparams.nz));//TODO note, not correct for non-equidistant grids

    const real retval = res[0];
    free(res);

    return retval;   
}


void get_slice_cuda_generic(Slice* h_slice)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    //For each GPU in the node
    //#pragma unroll
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
        update_slice_cuda_generic(&ctx->d_slice, &ctx->d_grid, &ctx->d_cparams, &h_run_params);
    }

    //For each GPU in the node
    //#pragma unroll
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
        store_slice_cuda_core(h_slice, &h_cparams, &h_run_params, &ctx->d_slice, &ctx->d_cparams, &ctx->start_idx);
    }
}


void load_forcing_params_cuda_generic(ForcingParams* forcing_params)
{
    //#pragma unroll
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        load_forcing_dconsts_cuda_core(forcing_params);
    }
}


void load_outer_halos_cuda_generic(Grid* h_grid, real* h_halobuffer)
{
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
        load_outer_halo_cuda_core(&ctx->d_grid, ctx->d_halobuffer, &ctx->d_cparams, 
                                  h_grid, h_halobuffer, &h_cparams, &ctx->start_idx);
    }
}


void store_internal_halos_cuda_generic(Grid* h_grid, real* h_halobuffer)
{
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
        store_internal_halo_cuda_core(h_grid, h_halobuffer, &h_cparams, &ctx->start_idx,
                                     &ctx->d_grid, ctx->d_halobuffer, &ctx->d_cparams);
    }
}

























































