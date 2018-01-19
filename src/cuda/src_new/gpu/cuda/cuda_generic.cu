/*
*   Implementation for the generic cuda solution.
*   Manages multiple GPUs on a single node using single-GPU implementations
*	defined in cuda subdirectories (cuda/core, cuda/generic etc)
*/
#include "cuda_generic.cuh"

#include "core/cuda_core.cuh"
#include "core/dconsts_core.cuh"
#include "core/errorhandler_cuda.cuh"
#include "core/concur_cuda_core.cuh"

#include "generic/collectiveops_cuda_generic.cuh"
#include "generic/rk3_cuda_generic.cuh"
#include "generic/boundcond_cuda_generic.cuh"
#include "generic/slice_cuda_generic.cuh"


/*
*	Host configs. 
*	These contain the information of the whole grid stored in this node.
* 	(f.ex. the grid dimensions before it has been decomposed for each GPU)
*/
static CParamConfig h_cparams;
static RunConfig h_run_params;
static bool is_initialized;

/*
*	Struct for storing the necessary information for running any of the
*	single-GPU implementations on some specific GPU. 
*
*	Contains f.ex. the "local" grid dimensions (d_cparams, dimensions after
*	decomposing the host grid) and the starting index (start_idx) in the host
*	grid for mapping from local grid to host grid coordinates.
*
*	These contexts are stored in a static global array, gpu_contexts.
*	If we wanted to integrate on, say, device 1 then we would set the current
* 	device to (1) and pass the necessary information from gpu_contexts[1] to the
*	integration kernel.
*
*/
typedef struct {
    vec3i           start_idx;    //The starting coordinates in the host array (node-wide grid)
    CParamConfig    d_cparams;    //Local CParamConfig for the device (GPU-specific grid)
    Grid            d_grid, d_grid_dst; //Core arrays
    Slice           d_slice;            //Slice arrays
    ReductionArray  d_reduct_arr;       //Reduction arrays
    real*           d_halobuffer;       //Buffer used for multi-node halo transfers
    ConcurContext   concur_ctx;         //Device-specific streams and events
} GPUContext;

static GPUContext* gpu_contexts;
static int num_devices = -1;


static inline void swap_ptrs(real** a, real** b)
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


static inline cudaStream_t get_stream(const int device_id, const StreamName str)
{
    return gpu_contexts[device_id].concur_ctx.streams[str];
}


static inline cudaEvent_t get_event(const int device_id, const EventName ev)
{
    return gpu_contexts[device_id].concur_ctx.events[ev];    
}


static void sync_devices() 
{
    int curr_device;
    cudaGetDevice(&curr_device);
    #pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
            cudaSetDevice(device_id);
	        cudaDeviceSynchronize();
    }
    cudaSetDevice(curr_device);
}


typedef enum {PEER_FRONT, PEER_BACK, NUM_PEERS} PeerType;
static int get_peer(PeerType pt, int device_id)
{
    switch (pt) {
        case PEER_FRONT:
            return (device_id+1) % num_devices;
        case PEER_BACK:
            return (num_devices+device_id-1) % num_devices;
        default:
            CRASH("Invalid PeerType");
    }
}


//TODO NOTE: peer access not supported between 4x p100, why?
//TEMP FIX: Commented out peer access enabling, now runs out-of-the-box
//on p100. Surprisingly peer access seems to work even without explicitly
//enabling it also on k80s.
static void set_peer_access(int device_id, bool enable_access)
{
    const int peer_front = get_peer(PEER_FRONT, device_id);
    const int peer_back  = get_peer(PEER_BACK, device_id);
    /*
    if (device_id != peer_front) {
        if (enable_access)
            cudaDeviceEnablePeerAccess(peer_front, 0);
        else
            cudaDeviceDisablePeerAccess(peer_front);
    }
    if (device_id != peer_back && peer_front != peer_back) {
        if (enable_access)
            cudaDeviceEnablePeerAccess(peer_back, 0);
        else
            cudaDeviceDisablePeerAccess(peer_back);    
    }*/
}


/*
*	Handles the allocation and initialization of the memories of all GPUs on
*	the node (incl. constant memory).
*/
__global__ void dummy_kernel() {}
void init_cuda_generic(CParamConfig* cparamconf, RunConfig* runconf)
{   
    if (is_initialized) { CRASH("cuda_generic already initialized!") }

    cudaGetDeviceCount(&num_devices);
    gpu_contexts = (GPUContext*) malloc(sizeof(GPUContext)*num_devices);
    printf("Using %d devices\n", num_devices);
    print_gpu_config_cuda_core();


    //Copy the structs in case the caller deallocates them prematurely
    h_cparams = *cparamconf;
    h_run_params = *runconf;

    //#pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];

        //printf("%d\n", __CUDA_ARCH__);
        printf("Trying to run a dummy kernel. If this fails, make sure that your\n"
                "device supports the CUDA architecture you are compiling for.\n"
                "Running dummy kernel... "); fflush(stdout);
        dummy_kernel<<<1, 1>>>();
        CUDA_ERRCHK_KERNEL_ALWAYS();
        printf("Success!\n");

        //Enable peer access
        set_peer_access(device_id, true);

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
        init_concur_ctx(&ctx->concur_ctx);
    }
    is_initialized = true;
}


/*
*	Deallocates all memory on the GPU
*/
void destroy_cuda_generic()
{
    if (!is_initialized) { CRASH("cuda_generic wasn't initialized!"); }
    
    //Sync all previous operations
    sync_devices();

    //Destroy everything
    #pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];    

        //Disable peer access
        set_peer_access(device_id, false);

        destroy_slice_cuda_generic(&ctx->d_slice);
        destroy_reduction_array_cuda_generic(&ctx->d_reduct_arr);
        destroy_grid_cuda_core(&ctx->d_grid, &ctx->d_grid_dst);
        destroy_halo_cuda_core(ctx->d_halobuffer);
        destroy_concur_ctx(&ctx->concur_ctx);
    }

    //Belt-and-suspenders-destroy-everything
    #pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        cudaDeviceReset();
    }

    free(gpu_contexts);

    is_initialized = false;
}


void load_grid_cuda_generic(Grid* h_grid)
{
    if (!is_initialized) { CRASH("cuda_generic wasn't initialized!") }

    //If we wanted to use another layout, we would do it here instead of using the core interface
    #pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
    
        load_grid_cuda_core(&ctx->d_grid, &ctx->d_cparams, &ctx->start_idx, h_grid, &h_cparams); 
    }
}


void store_grid_cuda_generic(Grid* h_grid)
{
    if (!is_initialized) { CRASH("cuda_generic wasn't initialized!") }

    #pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];

        store_grid_cuda_core(h_grid, &h_cparams, &ctx->d_grid, &ctx->d_cparams, &ctx->start_idx); 
    }
}


static void local_boundconds_cuda_generic()
{
    #pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];

        //Do local boundaries and signal when done
        periodic_xy_boundconds_cuda_generic(&ctx->d_grid, &ctx->d_cparams, 0);    
        const cudaEvent_t local_bc_done = get_event(device_id, EVENT_LOCAL_BC_DONE);
        cudaEventRecord(local_bc_done, 0);//Implicit synchronization with the default stream        
    }
}


static void fetch_halos_cuda_generic(GPUContext* ctx, const int device_id, cudaStream_t stream=0)
{
    const int front_id = get_peer(PEER_FRONT, device_id);
    const int back_id  = get_peer(PEER_BACK,  device_id);

    const size_t slab_size           = ctx->d_cparams.mx * ctx->d_cparams.my;
    const size_t transfer_size_bytes = BOUND_SIZE * slab_size * sizeof(real);

    const size_t z_src0 = ctx->d_cparams.nz * slab_size;
    const size_t z_dst0 = 0; 
    const size_t z_src1 = BOUND_SIZE * slab_size;
    const size_t z_dst1 = (ctx->d_cparams.nz + BOUND_SIZE) * slab_size;
    for (int w=0; w < NUM_ARRS; ++w) {
        CUDA_ERRCHK( cudaMemcpyPeerAsync(&ctx->d_grid.arr[w][z_dst0], device_id, 
                                         &gpu_contexts[back_id].d_grid.arr[w][z_src0], back_id,
                                         transfer_size_bytes, stream) ); //Back
        CUDA_ERRCHK( cudaMemcpyPeerAsync(&ctx->d_grid.arr[w][z_dst1], device_id, 
                                         &gpu_contexts[front_id].d_grid.arr[w][z_src1], front_id,
                                         transfer_size_bytes, stream) ); //Front
    }
}


static void exchange_halos_cuda_generic()
{
    #pragma omp parallel for num_threads(num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];

        //Wait until front and back neighbors are done with local boundary conditions
        const cudaStream_t global_stream = get_stream(device_id, STREAM_GLOBAL);
        const int peer_front = get_peer(PEER_FRONT, device_id);
        const int peer_back  = get_peer(PEER_BACK, device_id);
        cudaStreamWaitEvent(global_stream, get_event(peer_front, EVENT_LOCAL_BC_DONE), 0);
        cudaStreamWaitEvent(global_stream, get_event(peer_back, EVENT_LOCAL_BC_DONE), 0);

        //Get the updated halos from the front and back neighbor
        fetch_halos_cuda_generic(ctx, device_id, global_stream);
    }
}


void boundcond_step_cuda_generic()
{
    if (!is_initialized) { CRASH("cuda_generic wasn't initialized!") }

    local_boundconds_cuda_generic();
    exchange_halos_cuda_generic();
}


void integrate_step_cuda_generic(int isubstep, real dt)
{
    if (!is_initialized) { CRASH("cuda_generic wasn't initialized!") }

    //For all GPUs in the node in parallel
    #pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];

        //Integrate
        rk3_inner_cuda_generic(&ctx->d_grid, &ctx->d_grid_dst, isubstep, dt, 
                               &ctx->d_cparams, 
                               ctx->concur_ctx.streams[STREAM_LOCAL_HYDRO], 
                               ctx->concur_ctx.streams[STREAM_LOCAL_INDUCT]); 
        //WARNING: boundcond_step must have been called before rk3_outer.
        //If fetch_halos_cuda_generic() is not already be scheduled for execution
        //on the GPU, then the execution order will be wrong
        rk3_outer_cuda_generic(&ctx->d_grid, &ctx->d_grid_dst, isubstep, dt, 
                               &ctx->d_cparams, 
                                ctx->concur_ctx.streams[STREAM_GLOBAL]);

        //Swap src and dst device array pointers
        swap_grid_ptrs(&ctx->d_grid, &ctx->d_grid_dst);
    }

    //WARNING: this sync is not absolutely necessary but left here for safety:
    //without sync the host caller is able to execute other (potentially dangerous)
    //code in parallel with the GPU integration/memory transfers
    sync_devices(); //WARNING
}


void integrate_cuda_generic(real dt)
{
    if (!is_initialized) { CRASH("cuda_generic wasn't initialized!") }

    for (int isubstep=0; isubstep < 3; ++isubstep) {

        boundcond_step_cuda_generic();
        integrate_step_cuda_generic(isubstep, dt);

        //The original concurrency code, left here since it's easier to read
        //when boundary conditions and integration are not split up into separate
        //functions
        /*
        //Local boundaries and integration in the inner domain
        #pragma omp parallel for num_threads (num_devices)
        for (int device_id=0; device_id < num_devices; ++device_id) {
            cudaSetDevice(device_id);
            GPUContext* ctx = &gpu_contexts[device_id];

            //Do local boundaries and signal when done
            periodic_xy_boundconds_cuda_generic(&ctx->d_grid, &ctx->d_cparams, 0);    
            const cudaEvent_t local_bc_done = get_event(device_id, EVENT_LOCAL_BC_DONE);
            cudaEventRecord(local_bc_done, 0);//Implicit synchronization with the default stream

            //Start integrating in the inner computational domain
            rk3_inner_cuda_generic(&ctx->d_grid, &ctx->d_grid_dst, isubstep, dt, 
                                   &ctx->d_cparams, 
                                   ctx->concur_ctx.streams[STREAM_LOCAL_HYDRO], 
                                   ctx->concur_ctx.streams[STREAM_LOCAL_INDUCT]);            
        }

        //Communication of the outer halos among devices
        #pragma omp parallel for num_threads(num_devices)
        for (int device_id=0; device_id < num_devices; ++device_id) {
            cudaSetDevice(device_id);
            GPUContext* ctx = &gpu_contexts[device_id];

            //Wait until front and back neighbors are done with local boundary conditions
            const cudaStream_t global_stream = get_stream(device_id, STREAM_GLOBAL);
            const int peer_front = get_peer(PEER_FRONT, device_id);
            const int peer_back  = get_peer(PEER_BACK, device_id);
            cudaStreamWaitEvent(global_stream, get_event(peer_front, EVENT_LOCAL_BC_DONE), 0);
            cudaStreamWaitEvent(global_stream, get_event(peer_back, EVENT_LOCAL_BC_DONE), 0);

            //Get the updated halos from the front and back neighbor
            fetch_halos_cuda_generic(ctx, device_id, global_stream);
        }

        //Integrate in the outer computational domain
        #pragma omp parallel for num_threads(num_devices)
        for (int device_id=0; device_id < num_devices; ++device_id) {
            cudaSetDevice(device_id);
            GPUContext* ctx = &gpu_contexts[device_id];

            //Start integrating the outer domain after the updated halos
            //have arrived from neighbors
            rk3_outer_cuda_generic(&ctx->d_grid, &ctx->d_grid_dst, isubstep, dt, 
                                   &ctx->d_cparams, 
                                    ctx->concur_ctx.streams[STREAM_GLOBAL]);

            //We're done, swap src and dst device array pointers
            swap_grid_ptrs(&ctx->d_grid, &ctx->d_grid_dst);
        }*/
    }
}


#include "utils/utils.h" //For max/min/sum
real reduce_cuda_generic(ReductType t, GridType grid_type)
{
    if (!is_initialized) { CRASH("cuda_generic wasn't initialized!"); }

    real* res = (real*) malloc(sizeof(real)*num_devices);

    #pragma omp parallel for num_threads (num_devices)
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
            if (grid_type == NOT_APPLICABLE) { CRASH("Invalid GridType in reduce_cuda_generic"); }
            res[device_id] = get_reduction_cuda_generic(&ctx->d_reduct_arr, t, &ctx->d_cparams, ctx->d_grid.arr[grid_type]);
        }
    }

    //Bruteforce: find max, min or rms from the gpu results
    ////#pragma omp parallel  target teams distribute parallel for reduction(+:r)//TODO
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
    if (!is_initialized) { CRASH("cuda_generic wasn't initialized!"); }

    #pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
        update_slice_cuda_generic(&ctx->d_slice, &ctx->d_grid, &ctx->d_cparams, &h_run_params);
        cudaDeviceSynchronize();
        store_slice_cuda_core(h_slice, &h_cparams, &h_run_params, &ctx->d_slice, &ctx->d_cparams, &ctx->start_idx);
    }

//cd src/build/ && make -j && ac_srun_taito_multigpu 4 && cd ../../ && screen py_animate_data --nslices=100
}


void load_forcing_params_cuda_generic(ForcingParams* forcing_params)
{
    #pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        load_forcing_dconsts_cuda_core(forcing_params);
    }
}


void load_outer_halos_cuda_generic(Grid* h_grid, real* h_halobuffer)
{
    ////#pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
        load_outer_halo_cuda_core(&ctx->d_grid, ctx->d_halobuffer, &ctx->d_cparams, 
                                  h_grid, h_halobuffer, &h_cparams, &ctx->start_idx);
    }
}


void store_internal_halos_cuda_generic(Grid* h_grid, real* h_halobuffer)
{
    ////#pragma omp parallel for num_threads (num_devices)
    for (int device_id=0; device_id < num_devices; ++device_id) {
        cudaSetDevice(device_id);
        GPUContext* ctx = &gpu_contexts[device_id];
        store_internal_halo_cuda_core(h_grid, h_halobuffer, &h_cparams, &ctx->start_idx,
                                     &ctx->d_grid, ctx->d_halobuffer, &ctx->d_cparams);
    }
}

























































