#include "cuda_core.cuh"
#define INCLUDED_FROM_DCONST_DEFINER
#include "dconsts_core.cuh"
#include "errorhandler_cuda.cuh"
#include "common/config.h"


void print_gpu_config()
{
    int n_devices;
    if (cudaGetDeviceCount(&n_devices) != cudaSuccess) CRASH("No CUDA devices found!");
    printf("Num CUDA devices found: %u\n", n_devices);

    int initial_device;
    cudaGetDevice(&initial_device);
    for (int i = 0; i < n_devices; i++) {
        cudaSetDevice(i);

        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);
        printf("--------------------------------------------------\n");
        printf("Device Number: %d\n", i);
        printf("  Device name: %s\n", props.name);
        printf("  Compute capability: %d.%d\n", props.major, props.minor);
        printf("  Global memory\n");
        printf("    Memory Clock Rate (MHz): %d\n", props.memoryClockRate / (1000));
        printf("    Memory Bus Width (bits): %d\n", props.memoryBusWidth);
        printf("    Peak Memory Bandwidth (GiB/s): %f\n", 2.0*props.memoryClockRate*(props.memoryBusWidth/8)/(1024*1024));
        printf("    ECC enabled: %d\n", props.ECCEnabled);
        //Memory usage
        size_t free_bytes, total_bytes;
        CUDA_ERRCHK( cudaMemGetInfo(&free_bytes, &total_bytes) );
        const size_t used_bytes = total_bytes - free_bytes;     
        printf("    Total global mem: %.2f GiB\n", props.totalGlobalMem / (1024.0*1024*1024));
        printf("    Gmem used (GiB): %.2f\n", used_bytes / (1024.0*1024*1024));
        printf("    Gmem memory free (GiB): %.2f\n", free_bytes / (1024.0*1024*1024));
        printf("    Gmem memory total (GiB): %.2f\n", total_bytes / (1024.0*1024*1024));
        printf("  Caches\n");
        printf("    L2 size: %d KiB\n", props.l2CacheSize / (1024));
        printf("    Total const mem: %ld KiB\n", props.totalConstMem / (1024));
        printf("    Shared mem per block: %ld KiB\n", props.sharedMemPerBlock / (1024));
        printf("  Other\n");
        printf("    Warp size: %d\n", props.warpSize);
        //printf("    Single to double perf. ratio: %dx\n", props.singleToDoublePrecisionPerfRatio); //Not supported with older CUDA versions
        printf("    Stream priorities supported: %d\n", props.streamPrioritiesSupported);
        printf("--------------------------------------------------\n");
    }
    cudaSetDevice(initial_device);
}


void load_forcing_dconsts_cuda_core(ForcingParams* forcing_params)
{
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_FORCING_ENABLED, &forcing_params->forcing_enabled, sizeof(bool)) );
    //Copy forcing coefficients to the device's constant memory
    const size_t k_idx = forcing_params->k_idx;
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_KK_VEC_X, &forcing_params->kk_x[k_idx], sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_KK_VEC_Y, &forcing_params->kk_y[k_idx], sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_KK_VEC_Z, &forcing_params->kk_z[k_idx], sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_FORCING_KK_PART_X, &forcing_params->kk_part_x, sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_FORCING_KK_PART_Y, &forcing_params->kk_part_y, sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_FORCING_KK_PART_Z, &forcing_params->kk_part_z, sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_PHI, &forcing_params->phi, sizeof(real)) );
}


void load_hydro_dconsts_cuda_core(CParamConfig* cparams, RunConfig* run_params, const vec3i start_idx)
{ 
    //Grid dimensions
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_nx, &(cparams->nx), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_ny, &(cparams->ny), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_nz, &(cparams->nz), sizeof(int)) );

    CUDA_ERRCHK( cudaMemcpyToSymbol(d_mx, &(cparams->mx), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_my, &(cparams->my), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_mz, &(cparams->mz), sizeof(int)) );

    CUDA_ERRCHK( cudaMemcpyToSymbol(d_nx_min, &(cparams->nx_min), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_nx_max, &(cparams->nx_max), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_ny_min, &(cparams->ny_min), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_ny_max, &(cparams->ny_max), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_nz_min, &(cparams->nz_min), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_nz_max, &(cparams->nz_max), sizeof(int)) );

    CUDA_ERRCHK( cudaMemcpyToSymbol(d_DSX, &(cparams->dsx), sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_DSY, &(cparams->dsy), sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_DSZ, &(cparams->dsz), sizeof(real)) );

    const real dsx_offset = cparams->dsx*start_idx.x;
    const real dsy_offset = cparams->dsy*start_idx.y;
    const real dsz_offset = cparams->dsz*start_idx.z;
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_DSX_OFFSET, &dsx_offset, sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_DSY_OFFSET, &dsy_offset, sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_DSZ_OFFSET, &dsz_offset, sizeof(real)) ); 

    const real xorig = XORIG;
    const real yorig = YORIG;
    const real zorig = ZORIG;
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_XORIG, &xorig, sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_YORIG, &yorig, sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_ZORIG, &zorig, sizeof(real)) );    
    

    //Diff constants
    const real diff1_dx = 1.0/(60.0*cparams->dsx);
	const real diff1_dy = 1.0/(60.0*cparams->dsy);
	const real diff1_dz = 1.0/(60.0*cparams->dsz);

	const real diff2_dx = 1.0/(180.0*cparams->dsx*cparams->dsx);
	const real diff2_dy = 1.0/(180.0*cparams->dsy*cparams->dsy);
	const real diff2_dz = 1.0/(180.0*cparams->dsz*cparams->dsz);

	const real diffmn_dxdy = (1.0/720.0)*(1.0/cparams->dsx)*(1.0/cparams->dsy); 
	const real diffmn_dydz = (1.0/720.0)*(1.0/cparams->dsy)*(1.0/cparams->dsz);
	const real diffmn_dxdz = (1.0/720.0)*(1.0/cparams->dsz)*(1.0/cparams->dsx);

	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF1_DX_DIV, &diff1_dx, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF1_DY_DIV, &diff1_dy, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF1_DZ_DIV, &diff1_dz, sizeof(real)) );

	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF2_DX_DIV, &diff2_dx, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF2_DY_DIV, &diff2_dy, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF2_DZ_DIV, &diff2_dz, sizeof(real)) );

	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFFMN_DXDY_DIV, &diffmn_dxdy, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFFMN_DYDZ_DIV, &diffmn_dydz, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFFMN_DXDZ_DIV, &diffmn_dxdz, sizeof(real)) );


    //Viscosity
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_NU_VISC, &(run_params->nu_visc), sizeof(real)) );

    //Speed of sound
    const real cs2_sound = pow(run_params->cs_sound, 2.0);
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_CS2_SOUND, &cs2_sound, sizeof(real)) );	

    //Induction
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_ETA, &(run_params->eta), sizeof(real)) );
}


void init_grid_cuda_core(Grid* d_grid, Grid* d_grid_dst, CParamConfig* cparams)
{
    //Print the GPU configuration
    print_gpu_config();

    const size_t grid_size_bytes = sizeof(real) * cparams->mx * cparams->my * cparams->mz; 

    //Init device arrays
    for (int i=0; i < NUM_ARRS; ++i) {
        CUDA_ERRCHK( cudaMalloc(&(d_grid->arr[i]), grid_size_bytes) );
        CUDA_ERRCHK( cudaMemset(d_grid->arr[i], INT_MAX, grid_size_bytes) );
        CUDA_ERRCHK( cudaMalloc(&(d_grid_dst->arr[i]), grid_size_bytes) );
        CUDA_ERRCHK( cudaMemset(d_grid_dst->arr[i], INT_MAX, grid_size_bytes) );
    }
}


void destroy_grid_cuda_core(Grid* d_grid, Grid* d_grid_dst)
{
    for (int i=0; i < NUM_ARRS; ++i) {
        CUDA_ERRCHK( cudaFree(d_grid->arr[i]) );
        CUDA_ERRCHK( cudaFree(d_grid_dst->arr[i]) );
    }
}

void load_grid_cuda_core(Grid* d_grid, CParamConfig* d_cparams, vec3i* h_start_idx, Grid* h_grid, CParamConfig* h_cparams)
{
    const size_t grid_size_bytes = sizeof(real)* d_cparams->mx * d_cparams->my * d_cparams->mz;
    const size_t slice_size = h_cparams->mx*h_cparams->my;
    for (int w=0; w < NUM_ARRS; ++w) {
        CUDA_ERRCHK( cudaMemcpyAsync(&(d_grid->arr[w][0]), &(h_grid->arr[w][h_start_idx->z*slice_size]), grid_size_bytes, cudaMemcpyHostToDevice) );//NOTE: if stream not specified, uses the default stream->non-async behaviour
    }
}


void store_grid_cuda_core(Grid* h_grid, CParamConfig* h_cparams, Grid* d_grid, CParamConfig* d_cparams, vec3i* h_start_idx)
{
    const size_t grid_size_bytes = sizeof(real)* d_cparams->mx * d_cparams->my * d_cparams->nz;
    const size_t slice_size = h_cparams->mx * h_cparams->my;
    const size_t z_offset = BOUND_SIZE * slice_size;
    for (int w=0; w < NUM_ARRS; ++w) {
        CUDA_ERRCHK( cudaMemcpyAsync(&(h_grid->arr[w][z_offset + h_start_idx->z*slice_size]), &(d_grid->arr[w][z_offset]), grid_size_bytes, cudaMemcpyDeviceToHost) ); //NOTE: if stream not specified, uses the default stream->non-async behaviour
    }    
}


void store_slice_cuda_core(Slice* h_slice, CParamConfig* h_cparams, RunConfig* h_run_params, Slice* d_slice, CParamConfig* d_cparams, vec3i* h_start_idx)
{
    if (h_run_params->slice_axis != 'z') CRASH("Slice axis other that z not yet supported!");

    Slice buffer;
    slice_malloc(&buffer, d_cparams, h_run_params);

    cudaStream_t stream;
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    const size_t slice_size_bytes = sizeof(real) * d_cparams->mx * d_cparams->my;
    for (int w=0; w < NUM_SLICES; ++w)
        CUDA_ERRCHK( cudaMemcpyAsync(buffer.arr[w], d_slice->arr[w], slice_size_bytes, cudaMemcpyDeviceToHost, stream) );
    
    const size_t row_size_bytes = sizeof(real) * d_cparams->nx;
    for (int w=0; w < NUM_SLICES; ++w)
        for (int j=d_cparams->ny_min; j < d_cparams->ny_max; ++j)
            memcpy(&h_slice->arr[w][h_cparams->nx_min + (j+h_start_idx->y)*h_cparams->mx], 
                   &buffer.arr[w][d_cparams->nx_min + j*d_cparams->mx], row_size_bytes);
    
    cudaStreamDestroy(stream);
    slice_free(&buffer);
}































































































