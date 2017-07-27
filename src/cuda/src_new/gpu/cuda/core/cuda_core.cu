#define INCLUDED_FROM_CUDA_CORE
#include "cuda_core.cuh"
#include "errorhandler_cuda.cuh"
#include "common/config.h"
#include "common/errorhandler.h"

#define INCLUDED_FROM_DCONST_DEFINER
#include "dconsts_core.cuh"


static CParamConfig* cparams = NULL;
static RunConfig* run_params = NULL;
static bool is_initialized = false;


void print_gpu_config()
{
    int n_devices;
    if (cudaGetDeviceCount(&n_devices) != cudaSuccess) CRASH("No CUDA devices found!");

    printf("Num CUDA devices found: %u\n", n_devices);

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
}


void load_dconsts_core()
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

    const int bound_size = BOUND_SIZE;
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_bound_size, &bound_size, sizeof(real)) );


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
}


void init_cuda_core(CParamConfig* cparamconf, RunConfig* runconf)
{
    if (is_initialized) CRASH("init_cuda_core() already initialized!");
    is_initialized = true;

    //Print the GPU configuration
    //print_gpu_config();

    //Store config
    cparams = cparamconf;
    run_params = runconf;

    //Load core constants
    load_dconsts_core();

    const int grid_size = cparams->mx * cparams->my * cparams->mz; 

    //Init device arrays
	CUDA_ERRCHK( cudaMalloc(&d_lnrho, sizeof(real)*grid_size) );
	CUDA_ERRCHK( cudaMalloc(&d_uu_x,  sizeof(real)*grid_size) );
	CUDA_ERRCHK( cudaMalloc(&d_uu_y,  sizeof(real)*grid_size) );
	CUDA_ERRCHK( cudaMalloc(&d_uu_z,  sizeof(real)*grid_size) );

	CUDA_ERRCHK( cudaMalloc(&d_lnrho_dst, sizeof(real)*grid_size) );
	CUDA_ERRCHK( cudaMalloc(&d_uu_x_dst,  sizeof(real)*grid_size) );
	CUDA_ERRCHK( cudaMalloc(&d_uu_y_dst,  sizeof(real)*grid_size) );
	CUDA_ERRCHK( cudaMalloc(&d_uu_z_dst,  sizeof(real)*grid_size) );

    CUDA_ERRCHK( cudaMemset(d_lnrho, INT_MAX, sizeof(real)*grid_size) );
    CUDA_ERRCHK( cudaMemset(d_uu_x,  INT_MAX, sizeof(real)*grid_size) );
    CUDA_ERRCHK( cudaMemset(d_uu_y,  INT_MAX, sizeof(real)*grid_size) );
    CUDA_ERRCHK( cudaMemset(d_uu_z,  INT_MAX, sizeof(real)*grid_size) );

    CUDA_ERRCHK( cudaMemset(d_lnrho_dst, INT_MAX, sizeof(real)*grid_size) );
    CUDA_ERRCHK( cudaMemset(d_uu_x_dst,  INT_MAX, sizeof(real)*grid_size) );
    CUDA_ERRCHK( cudaMemset(d_uu_y_dst,  INT_MAX, sizeof(real)*grid_size) );
    CUDA_ERRCHK( cudaMemset(d_uu_z_dst,  INT_MAX, sizeof(real)*grid_size) );
}


void destroy_cuda_core()
{
    if (!is_initialized) CRASH("destroy_cuda_core() wasn't initialized!");
    is_initialized = false;

    //Get rid of the config ptr
    cparams = NULL;
    run_params = NULL;

    //Free device arrays
    CUDA_ERRCHK( cudaFree(d_lnrho) );
    CUDA_ERRCHK( cudaFree(d_uu_x)  );
    CUDA_ERRCHK( cudaFree(d_uu_y)  );
    CUDA_ERRCHK( cudaFree(d_uu_z)  );

    CUDA_ERRCHK( cudaFree(d_lnrho_dst) );
    CUDA_ERRCHK( cudaFree(d_uu_x_dst)  );
    CUDA_ERRCHK( cudaFree(d_uu_y_dst)  );
    CUDA_ERRCHK( cudaFree(d_uu_z_dst)  );
}


void load_grid_cuda_core(real* lnrho, real* uu_x, real* uu_y, real* uu_z)
{
    if (!is_initialized) CRASH("destroy_cuda_core() wasn't initialized!");

    const int grid_size = cparams->mx * cparams->my * cparams->mz; 
    CUDA_ERRCHK( cudaMemcpy(d_lnrho, lnrho, sizeof(real)*grid_size, cudaMemcpyHostToDevice) );
	CUDA_ERRCHK( cudaMemcpy(d_uu_x,  uu_x,  sizeof(real)*grid_size, cudaMemcpyHostToDevice) );
	CUDA_ERRCHK( cudaMemcpy(d_uu_y,  uu_y,  sizeof(real)*grid_size, cudaMemcpyHostToDevice) );
	CUDA_ERRCHK( cudaMemcpy(d_uu_z,  uu_z,  sizeof(real)*grid_size, cudaMemcpyHostToDevice) );
}


void store_grid_cuda_core(real* lnrho, real* uu_x, real* uu_y, real* uu_z)
{
    if (!is_initialized) CRASH("destroy_cuda_core() wasn't initialized!");

    const int grid_size = cparams->mx * cparams->my * cparams->mz; 
    CUDA_ERRCHK( cudaMemcpy(lnrho, d_lnrho, sizeof(real)*grid_size, cudaMemcpyDeviceToHost) );
	CUDA_ERRCHK( cudaMemcpy(uu_x,  d_uu_x,  sizeof(real)*grid_size, cudaMemcpyDeviceToHost) );
	CUDA_ERRCHK( cudaMemcpy(uu_y,  d_uu_y,  sizeof(real)*grid_size, cudaMemcpyDeviceToHost) );
	CUDA_ERRCHK( cudaMemcpy(uu_z,  d_uu_z,  sizeof(real)*grid_size, cudaMemcpyDeviceToHost) );
}






























































































