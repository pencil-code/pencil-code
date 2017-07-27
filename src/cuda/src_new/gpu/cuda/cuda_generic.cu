/*
*   Implementation for the generic cuda solution.
*   Manages multiple GPUs on a single node.
*/
#include "cuda_generic.cuh"
#include "core/dconsts_core.cuh"
#include "core/cuda_core.cuh"
#include "generic/collectiveops_cuda_generic.cuh"
#include "generic/rk3_cuda_generic.cuh"
#include "generic/boundcond_cuda_generic.cuh"
#include "generic/slice_cuda_generic.cuh"

#include "common/errorhandler.h"

static CParamConfig cparams;
static RunConfig run_params;
static bool is_initialized;


void init_cuda_generic(CParamConfig* cparamconf, RunConfig* runconf)
{   
    if (is_initialized) CRASH("cuda_generic already initialized!");
    is_initialized = true;

    //Copy the structs in case the caller deallocates them prematurely
    cparams = *cparamconf;
    run_params = *runconf;

    init_cuda_core(&cparams, &run_params);
    init_collectiveops_cuda_generic(&cparams);
    init_rk3_cuda_generic(&cparams);
    init_boundcond_cuda_generic(&cparams);
    init_slice_cuda_generic(&cparams, &run_params);
}


void destroy_cuda_generic()
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");
    is_initialized = false;

    destroy_slice_cuda_generic();
    destroy_boundcond_cuda_generic();
    destroy_rk3_cuda_generic();
    destroy_collectiveops_cuda_generic();
    destroy_cuda_core();
    cudaDeviceSynchronize(); cudaDeviceReset();
}


void load_grid_cuda_generic(real* lnrho, real* uu_x, real* uu_y, real* uu_z)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    //If we wanted to use another layout, we would do it here instead of using the core interface
    load_grid_cuda_core(lnrho, uu_x, uu_y, uu_z); 
}


void store_grid_cuda_generic(real* lnrho, real* uu_x, real* uu_y, real* uu_z)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    store_grid_cuda_core(lnrho, uu_x, uu_y, uu_z);
}


inline void swap_ptrs(real** a, real** b)
{
	real* temp = *a;
	*a = *b;
	*b = temp;
}


inline void swap_grid_ptrs()
{
    swap_ptrs(&d_lnrho, &d_lnrho_dst);
    swap_ptrs(&d_uu_x,  &d_uu_x_dst);
    swap_ptrs(&d_uu_y,  &d_uu_y_dst);
    swap_ptrs(&d_uu_z,  &d_uu_z_dst);
}


void integrate_step_cuda_generic(int isubstep, real dt)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    //For all GPUs in the node

    //Compute local boundary conditions (share overlapping boundaries between GPUs)
    boundcond_cuda_generic(d_lnrho, d_uu_x, d_uu_y, d_uu_z);

    //Integrate
    rk3_cuda_generic(d_lnrho,     d_uu_x,     d_uu_y,   d_uu_z, 
                     d_lnrho_dst, d_uu_x_dst, d_uu_y_dst, d_uu_z_dst,
                     isubstep, dt);

    //Swap src and dst device array pointers
    swap_grid_ptrs();
}


void integrate_cuda_generic(real dt)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    //Step 1
    integrate_step_cuda_generic(0, dt);

    //Step 2
    integrate_step_cuda_generic(1, dt);   

    //Step 3
    integrate_step_cuda_generic(2, dt);
}


real reduce_cuda_generic(ReductType t, ArrType a)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    real res = NAN;
    if (t == MAX_VEC || t == MIN_VEC || t == RMS_VEC) {
        if (a != UU) {
            printf("Note: other than UU passed to reduce_cuda_generic as ArrType." 
                   "This has no effect when a vector ReductType is selected\n");
        }
        res = get_reduction_cuda_generic(t, d_uu_x, d_uu_y, d_uu_z);
    } else {
        switch (a) {
            case LNRHO:
                res = get_reduction_cuda_generic(t, d_lnrho);
                break;
            case UU_X:
                res = get_reduction_cuda_generic(t, d_uu_x);
                break;
            case UU_Y:
                res = get_reduction_cuda_generic(t, d_uu_y);
                break;
            case UU_Z:
                res = get_reduction_cuda_generic(t, d_uu_z);
                break;
            default:
                CRASH("Invalid ArrType"); 
        }
    }
    return res;   
}


void get_slice_cuda_generic(real* slice_lnrho, real* slice_uu, real* slice_uu_x, real* slice_uu_y, real* slice_uu_z)
{
    if (!is_initialized) CRASH("cuda_generic wasn't initialized!");

    //For each GPU in the node
        update_slice_cuda_generic(d_lnrho, d_uu_x, d_uu_y, d_uu_z);

    //For each GPU in the node
        store_slice_cuda_generic(slice_lnrho, slice_uu, slice_uu_x, slice_uu_y, slice_uu_z);
}


























































