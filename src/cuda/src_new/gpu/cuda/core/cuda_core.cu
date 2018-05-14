#include "cuda_core.cuh"
#ifdef GPU_ASTAROTH
  #include "common/PC_moduleflags.h"
  #define EXTERN
  #include "common/PC_module_parfuncs.h"
#endif
#define INCLUDED_FROM_DCONST_DEFINER
#include "dconsts_core.cuh"
#include "errorhandler_cuda.cuh"
#include "copyHalosConcur.cuh"

#ifdef FORCING
    //Copy forcing coefficients to the device's constant memory
#ifdef GPU_ASTAROTH
void update_forcing_coefs_cuda_PC(ForcingParams* fp, CParamConfig* cparams, int start_idx)
{
    int comSize=3*sizeof(real);
    CUDA_ERRCHK( cudaMemcpyToSymbol( d_FORCING_COEF1, fp->coef1, comSize ));
    CUDA_ERRCHK( cudaMemcpyToSymbol( d_FORCING_COEF2, fp->coef2, comSize ));
    CUDA_ERRCHK( cudaMemcpyToSymbol( d_FORCING_COEF3, fp->coef3, comSize ));
    CUDA_ERRCHK( cudaMemcpyToSymbol( d_FORCING_FDA, fp->fda, comSize ));
    comSize=2*sizeof(real);  // fx, fy, fz are semantically complex !
    CUDA_ERRCHK( cudaMemcpyToSymbol( d_FORCING_FX, fp->fx, comSize*cparams->mx ));
    CUDA_ERRCHK( cudaMemcpyToSymbol( d_FORCING_FY, fp->fy, comSize*cparams->my ));
    CUDA_ERRCHK( cudaMemcpyToSymbol( d_FORCING_FZ, fp->fz+2*start_idx, comSize*cparams->mz ));
}
#else
void load_forcing_dconsts_cuda_core(ForcingParams* forcing_params)
{
    const size_t k_idx = forcing_params->k_idx;
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_KK_VEC_X, &forcing_params->kk_x[k_idx], sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_KK_VEC_Y, &forcing_params->kk_y[k_idx], sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_KK_VEC_Z, &forcing_params->kk_z[k_idx], sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_FORCING_KK_PART_X, &forcing_params->kk_part_x, sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_FORCING_KK_PART_Y, &forcing_params->kk_part_y, sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_FORCING_KK_PART_Z, &forcing_params->kk_part_z, sizeof(real)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_PHI, &forcing_params->phi, sizeof(real)) );
}
#endif
#endif

void load_hydro_dconsts_cuda_core(CParamConfig* cparams, RunConfig* run_params, const vec3i start_idx)
{ 
    //Grid dimensions
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_nx, &(cparams->nx), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_ny, &(cparams->ny), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_nz, &(cparams->nz), sizeof(int)) );

    CUDA_ERRCHK( cudaMemcpyToSymbol(d_mx, &(cparams->mx), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_my, &(cparams->my), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_mz, &(cparams->mz), sizeof(int)) );
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_mxy, &(cparams->mxy), sizeof(int)) );

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
//printf("diff1_dx, diff1_dy, diff1_dz= %f %f %f \n",diff1_dx, diff1_dy, diff1_dz);
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF1_DX_DIV, &diff1_dx, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF1_DY_DIV, &diff1_dy, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF1_DZ_DIV, &diff1_dz, sizeof(real)) );

	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF2_DX_DIV, &diff2_dx, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF2_DY_DIV, &diff2_dy, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFF2_DZ_DIV, &diff2_dz, sizeof(real)) );

	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFFMN_DXDY_DIV, &diffmn_dxdy, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFFMN_DYDZ_DIV, &diffmn_dydz, sizeof(real)) );
	CUDA_ERRCHK( cudaMemcpyToSymbol(d_DIFFMN_DXDZ_DIV, &diffmn_dxdz, sizeof(real)) );

#ifdef GPU_ASTAROTH
  #include "common/PC_modulepars.h"
printf("lpg= %d \n",lpressuregradient_gas);
//printf("lpg= %d \n",*((int *)p_par_hydro[0]));
printf("nu= %f\n",nu); 
//printf("p_par_viscosity[1-1]= %f\n",*(p_par_viscosity[1-1])); 
printf("cs20= %f\n",cs20); 
//printf("alpha_ts= %f %f %f \n",p_par_timestep[0][0],p_par_timestep[0][1],p_par_timestep[0][2]); 
//printf("beta_ts= %f %f %f \n",p_par_timestep[1][0],p_par_timestep[1][1],p_par_timestep[1][2]); 
printf("alpha_ts= %f %f %f \n",alpha_ts[0],alpha_ts[1],alpha_ts[2]); 
printf("beta_ts= %f %f %f \n",beta_ts[0],beta_ts[1],beta_ts[2]); 
#else
  #ifdef VISCOSITY
    //Viscosity
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_NU, &(run_params->nu_visc), sizeof(real)) );
  #endif

  #ifdef HYDRO
    //Speed of sound
    real cs2=pow(run_params->cs_sound);
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_CS20, &cs2, sizeof(real)) );	
  #endif

  #ifdef MAGNETIC
    //Diffusivity
    CUDA_ERRCHK( cudaMemcpyToSymbol(d_ETA, &(run_params->eta), sizeof(real)) );
  #endif
#endif
}

void init_grid_cuda_core(Grid* d_grid, Grid* d_grid_dst, CParamConfig* cparams)
{
    const size_t grid_size_bytes = sizeof(real) * cparams->mw;

    //Init device arrays
    for (int i=0; i < d_grid->NUM_ARRS; ++i) {
        CUDA_ERRCHK( cudaMalloc(&(d_grid->arr[i]), grid_size_bytes) );
        CUDA_ERRCHK( cudaMemset(d_grid->arr[i], INT_MAX, grid_size_bytes) );  //MR: What is INT_MAX?
        CUDA_ERRCHK( cudaMalloc(&(d_grid_dst->arr[i]), grid_size_bytes) );
        CUDA_ERRCHK( cudaMemset(d_grid_dst->arr[i], INT_MAX, grid_size_bytes) );
    }
}

void destroy_grid_cuda_core(Grid* d_grid, Grid* d_grid_dst)
{
    for (int i=0; i < d_grid->NUM_ARRS; ++i) {
        CUDA_ERRCHK( cudaFree(d_grid->arr[i]) );
        CUDA_ERRCHK( cudaFree(d_grid_dst->arr[i]) );
    }
}

void load_grid_cuda_core(Grid* d_grid, CParamConfig* d_cparams, vec3i* h_start_idx, Grid* h_grid, CParamConfig* h_cparams)
{
    // Loads full data cubes to device(s) with decomposition
    const size_t grid_size_bytes = sizeof(real)* d_cparams->mw;
    const size_t slice_size = h_cparams->mxy;
    for (int w=0; w < d_grid->NUM_ARRS; ++w) {
        CUDA_ERRCHK( cudaMemcpy(&(d_grid->arr[w][0]), &(h_grid->arr[w][h_start_idx->z*slice_size]), grid_size_bytes, cudaMemcpyHostToDevice) );
    }
}

void store_grid_cuda_core(Grid* h_grid, CParamConfig* h_cparams, Grid* d_grid, CParamConfig* d_cparams, vec3i* h_start_idx)
{
    // Stores full data cubes from device(s) with decomposition to host
    const size_t grid_size_bytes = sizeof(real)* d_cparams->mxy * d_cparams->nz;
    const size_t slice_size = h_cparams->mxy;
    const size_t z_offset = BOUND_SIZE * slice_size;
    for (int w=0; w < d_grid->NUM_ARRS; ++w) {
//printf("w,z_offset,h_start_idx->z*slice_size= %d %d %d \n",w,z_offset,h_start_idx->z*slice_size);
//printf("w,&(h_grid->arr[w][z_offset + h_start_idx->z*slice_size]), &(d_grid->arr[w][z_offset])= %d %d %d \n",w,&(h_grid->arr[w][z_offset + h_start_idx->z*slice_size]), &(d_grid->arr[w][z_offset]));
        CUDA_ERRCHK( cudaMemcpy(&(h_grid->arr[w][z_offset + h_start_idx->z*slice_size]), &(d_grid->arr[w][z_offset]), grid_size_bytes, cudaMemcpyDeviceToHost) );
        /*printf("w= %d \n",w);int indx=0;
        for (int iz=0;iz<38;iz++){
          for (int iy=0;iy<38;iy++){
            for (int ix=0;ix<38;ix++){
              printf("%f ",h_grid->arr[w][indx]); indx++;
            }
            printf("\n");
          }
        }
        printf("-----------------\n",w);*/
    }    
}

void store_slice_cuda_core(Slice* h_slice, CParamConfig* h_cparams, RunConfig* h_run_params, Slice* d_slice, CParamConfig* d_cparams, vec3i* h_start_idx)
{
    if (h_run_params->slice_axis != 'z') { CRASH("Slice axis other that z not yet supported!"); }

    Slice buffer;
    slice_malloc(&buffer, d_cparams, h_run_params);

    const size_t slice_size_bytes = sizeof(real) * d_cparams->mxy;
    for (int w=0; w < NUM_SLICES; ++w)
        CUDA_ERRCHK( cudaMemcpy(buffer.arr[w], d_slice->arr[w], slice_size_bytes, cudaMemcpyDeviceToHost) );
    
    const size_t row_size_bytes = sizeof(real) * d_cparams->nx;
    for (int w=0; w < NUM_SLICES; ++w)
        for (int j=d_cparams->ny_min; j < d_cparams->ny_max; ++j)
            memcpy(&h_slice->arr[w][h_cparams->nx_min + (j+h_start_idx->y)*h_cparams->mx], 
                   &buffer.arr[w][d_cparams->nx_min + j*d_cparams->mx], row_size_bytes);
    
    slice_free(&buffer);
}

void init_halo_cuda_core(GPUContext & ctx, bool lfirstGPU, bool llastGPU)
{
    //printf("init_halo_cuda_core\n");
// Allocate memory for d_halo here
    initHaloConcur(ctx, lfirstGPU, llastGPU);
}

void destroy_halo_cuda_core(GPUContext & ctx)
{
    //printf("destroy_halo_cuda_core\n");
    CUDA_ERRCHK( cudaFree(ctx.d_halobuffer) );

    for (int i=0; i<NUM_COPY_STREAMS; i++) {
      CUDA_ERRCHK(cudaStreamDestroy(ctx.d_copy_streams[i]));
    }
}

void load_outer_halo_cuda_core(const GPUContext & ctx, Grid* h_grid, real* h_halobuffer, bool lfirstGPU, bool llastGPU)
{
    //printf("load_outer_halo_cuda_core\n");

    // lock buffer for left & right plates
    lockHostMemyz();

    //Update the outer halos of the device  by copying from host.
    for (int w=0; w < h_grid->NUM_ARRS; ++w){ 
        copyOxyPlates(ctx,w,h_grid->arr[w],lfirstGPU,llastGPU);
        copyOxzPlates(ctx,w,h_grid->arr[w],lfirstGPU,llastGPU);
        copyOyzPlates(ctx,w,h_grid->arr[w],lfirstGPU,llastGPU);
    }
    //CUDA_ERRCHK_KERNEL();
    synchronizeStreams(ctx,lfirstGPU,llastGPU);
 
    // unlock buffer for left & right plates
    unlockHostMemyz();

    //!!!for (int w=0; w < h_grid->NUM_ARRS; ++w)   // time-critical!
        //!!!unlockHostMemOuter(ctx,h_grid->arr[w],lfirstGPU,llastGPU);
}

void store_internal_halo_cuda_core(const GPUContext & ctx, Grid* h_grid, real* h_halobuffer, bool lfirstGPU, bool llastGPU)
{
    //printf("store_internal_halo_cuda_core\n");
    //Store internal halos in d_halobuffer to host memory in h_halobuffer/h_grid

    // lock buffer for left & right plates
    lockHostMemyz();

    //Update the outer halos of the device  by copying from host.
    for (int w=0; w < h_grid->NUM_ARRS; ++w){
        copyIxyPlates(ctx,w,h_grid->arr[w],lfirstGPU,llastGPU);
        copyIxzPlates(ctx,w,h_grid->arr[w],lfirstGPU,llastGPU);
        copyIyzPlates(ctx,w,h_grid->arr[w],lfirstGPU,llastGPU);
    }
    synchronizeStreams(ctx,lfirstGPU,llastGPU);

    // unlock buffer for left & right plates
    unlockHostMemyz();

    // unlock all other locked memory
    //!!!for (int w=0; w < h_grid->NUM_ARRS; ++w)
        //!!!unlockHostMemInner(ctx,h_grid->arr[w],lfirstGPU,llastGPU);
}

void print_gpu_config_cuda_core()
{
    int n_devices;
    if (cudaGetDeviceCount(&n_devices) != cudaSuccess) { CRASH("No CUDA devices found!"); }
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

