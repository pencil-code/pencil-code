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

        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("Device Number: %d\n", i);
        printf("  Device name: %s\n", prop.name);
        printf("  Memory Clock Rate (KHz): %d\n", prop.memoryClockRate);
        printf("  Memory Bus Width (bits): %d\n", prop.memoryBusWidth);
        printf("  Peak Memory Bandwidth (GiB/s): %f\n", 2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/(1024*1024));

        //Memory usage
        size_t free_bytes, total_bytes;
        CUDA_ERRCHK( cudaMemGetInfo(&free_bytes, &total_bytes) );
        const size_t used_bytes = total_bytes - free_bytes;     
        printf("  GPU memory used (MiB): %f\n", (double) used_bytes / (1024*1024));
        printf("  GPU memory free (MiB): %f\n", (double) free_bytes / (1024*1024));
        printf("  GPU memory total (MiB): %f\n", (double) total_bytes / (1024*1024));
    }
    cudaSetDevice(initial_device);
    
    if (n_devices < NUM_DEVICES) CRASH("Invalid number of devices requested!");
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
    //Create a host buffer to minimize the number of device-host-device memcpys (very high latency)
    Grid buffer;
    grid_malloc(&buffer, d_cparams);

    const size_t slab_size_bytes = sizeof(real) * d_cparams->mx * d_cparams->my;
    for (int w=0; w < NUM_ARRS; ++w)
        for (int k=0; k < d_cparams->mz; ++k)
            memcpy(&buffer.arr[w][k*d_cparams->mx*d_cparams->my], &h_grid->arr[w][h_start_idx->y*h_cparams->mx + k*h_cparams->mx*h_cparams->my], slab_size_bytes);

    cudaStream_t stream;
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    const size_t grid_size_bytes = sizeof(real)* d_cparams->mx * d_cparams->my * d_cparams->mz; 
    for (int w=0; w < NUM_ARRS; ++w)
        CUDA_ERRCHK( cudaMemcpyAsync(d_grid->arr[w], buffer.arr[w], grid_size_bytes, cudaMemcpyHostToDevice, stream) );

    cudaStreamDestroy(stream);
    grid_free(&buffer);
}


void store_grid_cuda_core(Grid* h_grid, CParamConfig* h_cparams, Grid* d_grid, CParamConfig* d_cparams, vec3i* h_start_idx)
{
    Grid buffer;
    grid_malloc(&buffer, d_cparams);

    cudaStream_t stream;
    cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    const size_t grid_size_bytes = sizeof(real) * d_cparams->mx * d_cparams->my * d_cparams->mz;
    for (int w=0; w < NUM_ARRS; ++w)
        cudaMemcpyAsync(buffer.arr[w], d_grid->arr[w], grid_size_bytes, cudaMemcpyDeviceToHost, stream);

    const size_t row_size_bytes = sizeof(real) * d_cparams->nx; 
    for (int w=0; w < NUM_ARRS; ++w) {
        for (int k=d_cparams->nz_min; k < d_cparams->nz_max; ++k)
            for (int j=d_cparams->ny_min; j < d_cparams->ny_max; ++j)
                memcpy(&h_grid->arr[w][h_cparams->nx_min + (j+h_start_idx->y)*h_cparams->mx + k*h_cparams->mx*h_cparams->my], 
                       &buffer.arr[w][d_cparams->nx_min + j*d_cparams->mx + k*d_cparams->mx*d_cparams->my], row_size_bytes);
    } 

    cudaStreamDestroy(stream);
    grid_free(&buffer);
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































































































