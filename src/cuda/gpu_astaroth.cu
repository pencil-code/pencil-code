/*                             gpu_astaroth_ansi.cu
                              --------------------
*/

/* Date:   8-Feb-2017
   Author: M. Rheinhardt
   Description:
 ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
Comments: 
DATE 17-Feb-2017: Omer Anjum: Added description of the functions
*/

//C libraries
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

//CUDA libraries
//#include <cfloat>

//Headers
#include "dconsts.cuh"
#include "integrators.cuh"
#include "boundcond.cuh"
#include "timestep.cuh"
#include "defines.h"
#include "io.h"
#include "slice.cuh"
#include "smem.cuh"
#include "forcing.cuh"
#include "copyhalos.cuh"


//DEBUG
#include "initutils.h"
#include "diagnostics.cuh"
#define dbug 0

//#define _GNU_SOURCE
#include <string.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <dlfcn.h>

//#include "headers_c.h" 
int halo_size, d_lnrho_size;
float *d_halo;
float *halo; 
float *d_lnrho, *d_uu_x, *d_uu_y, *d_uu_z;
float *d_w_lnrho, *d_w_uu_x, *d_w_uu_y, *d_w_uu_z;
float *d_lnrho_dest, *d_uu_x_dest, *d_uu_y_dest, *d_uu_z_dest;
float *d_umax, *d_umin, *d_partial_result;//Device pointer for max vec and partial result for the reductions
float *d_urms; //Device pointer for urms 
float *d_uxrms, *d_uyrms, *d_uzrms; //Device pointer for uxrms, uyrms, uzrms 
float *d_rhorms, *d_rhomax, *d_rhomin; //Device pointer for rhorms, rhomax, rhomin
float *d_uxmax, *d_uymax, *d_uzmax; //Device pointer for uxmax, uymax, uzmax
float *d_uxmin, *d_uymin, *d_uzmin; //Device pointer for uxmin, uymin, uzmin


/*cudaError_t checkErr(cudaError_t result) {
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s \n", 
            cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
  return result;
}*/

//Swaps pointers a,b (do xor swap if too slow)
inline void swap_ptrs(float** a, float** b)
{
	float* temp = *a;
	*a = *b;
	*b = temp;
}

float nu, cs2;
//Loads constants into device memory
void load_dconsts(float nu_visc, float cs2_sound){
	//---------Grid dims-----------------------
	const int nx = NX, ny = NY, nz = NZ;
	const int pad_size = PAD_SIZE, bound_size = BOUND_SIZE;
	const int 	comp_domain_size_x = COMP_DOMAIN_SIZE_X, 
			comp_domain_size_y = COMP_DOMAIN_SIZE_Y, 
			comp_domain_size_z = COMP_DOMAIN_SIZE_Z;
	//Needed for calculating averages
	const float nelements_float = COMP_DOMAIN_SIZE_X*COMP_DOMAIN_SIZE_Y*COMP_DOMAIN_SIZE_Z;

	const float 	domain_size_x = DOMAIN_SIZE_X, 
			domain_size_y = DOMAIN_SIZE_Y, 
			domain_size_z = DOMAIN_SIZE_Z;

	checkErr( cudaMemcpyToSymbol(d_NX, &nx, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_NY, &ny, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_NZ, &nz, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_PAD_SIZE, &pad_size, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_BOUND_SIZE, &bound_size, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_COMP_DOMAIN_SIZE_X, &comp_domain_size_x, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_COMP_DOMAIN_SIZE_Y, &comp_domain_size_y, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_COMP_DOMAIN_SIZE_Z, &comp_domain_size_z, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_NELEMENTS_FLOAT, &nelements_float, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_DOMAIN_SIZE_X, &domain_size_x, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DOMAIN_SIZE_Y, &domain_size_y, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DOMAIN_SIZE_Z, &domain_size_z, sizeof(float)) );

	const int h_w_grid_y_offset = COMP_DOMAIN_SIZE_X;
	const int h_w_grid_z_offset = COMP_DOMAIN_SIZE_X*COMP_DOMAIN_SIZE_Y;
	const int h_grid_y_offset = NX;
	const int h_grid_z_offset = NX*NY;

	checkErr( cudaMemcpyToSymbol(d_W_GRID_Y_OFFSET, &h_w_grid_y_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_W_GRID_Z_OFFSET, &h_w_grid_z_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_GRID_Y_OFFSET, &h_grid_y_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_GRID_Z_OFFSET, &h_grid_z_offset, sizeof(int)) );
	//-------------------------------------------

	//------Computational domain's bottom and top indices---------
	const int cx_top = CX_TOP, cy_top = CY_TOP, cz_top = CZ_TOP;
	const int cx_bot = CX_BOT, cy_bot = CY_BOT, cz_bot = CY_BOT;

	checkErr( cudaMemcpyToSymbol(d_CX_TOP, &cx_top, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CY_TOP, &cy_top, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CZ_TOP, &cz_top, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_CX_BOT, &cx_bot, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CY_BOT, &cy_bot, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CZ_BOT, &cz_bot, sizeof(int)) );
	//-------------------------------------------

	//-------------Real distance between grid points---------------
	const float dx = DX, dy = DY, dz = DZ; 
	checkErr( cudaMemcpyToSymbol(d_DX, &dx, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DY, &dy, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DZ, &dz, sizeof(float)) );
	//-------------------------------------------

	//----------Location of the grid coordinate origin-------------
	const float xorig = XORIG, yorig = YORIG, zorig = ZORIG;
        checkErr( cudaMemcpyToSymbol(d_XORIG, &xorig, sizeof(float)) );
        checkErr( cudaMemcpyToSymbol(d_YORIG, &yorig, sizeof(float)) );
        checkErr( cudaMemcpyToSymbol(d_ZORIG, &zorig, sizeof(float)) );
	//-------------------------------------------

	//----------Shearing parameters---------------
	const float q_shear = Q_SHEAR;
	const float omega = OMEGA; 
	const int interp_order = INTERP_ORDER;

	printf("compute INTERP_ORDER = %i \n", INTERP_ORDER);
	printf("compute interp_order = %i \n", interp_order); 

        checkErr( cudaMemcpyToSymbol(d_INTERP_ORDER, &interp_order, sizeof(int)) );
        checkErr( cudaMemcpyToSymbol(d_Q_SHEAR, &q_shear, sizeof(float)) );
        checkErr( cudaMemcpyToSymbol(d_OMEGA, &omega, sizeof(float)) );
	//-------------------------------------------

	//------------------Optional physics switches------------------------
	const int lforcing = LFORCING;
	const int lshear = LSHEAR;
	const int lcoriolis = LCORIOLIS;
        checkErr( cudaMemcpyToSymbol(d_LFORCING, &lforcing, sizeof(int)) );
        checkErr( cudaMemcpyToSymbol(d_LSHEAR, &lshear, sizeof(int)) );
        checkErr( cudaMemcpyToSymbol(d_LCORIOLIS, &lcoriolis, sizeof(int)) );
	//-------------------------------------------------------------------

	//------------Random constants for computation-----------------
	const float h_ALPHA1 = 0.0; 
	const float h_ALPHA2 = -0.53125; 
	const float h_ALPHA3 = -1.1851851851851851;
	
	checkErr( cudaMemcpyToSymbol(d_ALPHA1, &h_ALPHA1, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_ALPHA2, &h_ALPHA2, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_ALPHA3, &h_ALPHA3, sizeof(float)) );

	const float h_BETA1 = 0.25; 
	const float h_BETA2 = 0.88888888888888884;
	const float h_BETA3 = 0.75;

	checkErr( cudaMemcpyToSymbol(d_BETA1, &h_BETA1, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_BETA2, &h_BETA2, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_BETA3, &h_BETA3, sizeof(float)) );

	//const float nu_visc = NU_VISC;
	//const float cs2_sound = pow(CS_SOUND, 2.0);
	checkErr( cudaMemcpyToSymbol(d_NU_VISC, &nu_visc, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_CS2_SOUND, &cs2_sound, sizeof(float)) );	
	//-------------------------------------------

	//------------Constants for derivatives-----------------

	const float flt_9 = 9.0;
	const float flt_45 = 45.0; 
	const float flt_60 = 60.0; 

	const float flt_2 = 2.0; 
	const float flt_27 = 27.0; 
	const float flt_270 = 270.0; 
	const float flt_490 = 490.0; 
	const float flt_180 = 180.0; 

	const float diff1_dx = 1.0/(60.0*DX);
	const float diff1_dy = 1.0/(60.0*DY);
	const float diff1_dz = 1.0/(60.0*DZ);

	const float diff2_dx = 1.0/(180.0*DX*DX);
	const float diff2_dy = 1.0/(180.0*DY*DY);
	const float diff2_dz = 1.0/(180.0*DZ*DZ);

	const float diffmn_dxdy = (1.0/720.0)*(1.0/DX)*(1.0/DY); 
	const float diffmn_dydz = (1.0/720.0)*(1.0/DY)*(1.0/DZ);
	const float diffmn_dxdz = (1.0/720.0)*(1.0/DZ)*(1.0/DX);
	
	checkErr( cudaMemcpyToSymbol(d_FLT_9, &flt_9, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_45, &flt_45, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_60, &flt_60, sizeof(float)) );

	checkErr( cudaMemcpyToSymbol(d_FLT_2, &flt_2, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_27, &flt_27, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_270, &flt_270, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_490, &flt_490, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_180, &flt_180, sizeof(float)) );

	checkErr( cudaMemcpyToSymbol(d_DIFF1_DX_DIV, &diff1_dx, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFF1_DY_DIV, &diff1_dy, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFF1_DZ_DIV, &diff1_dz, sizeof(float)) );

	checkErr( cudaMemcpyToSymbol(d_DIFF2_DX_DIV, &diff2_dx, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFF2_DY_DIV, &diff2_dy, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFF2_DZ_DIV, &diff2_dz, sizeof(float)) );

	checkErr( cudaMemcpyToSymbol(d_DIFFMN_DXDY_DIV, &diffmn_dxdy, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFFMN_DYDZ_DIV, &diffmn_dydz, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFFMN_DXDZ_DIV, &diffmn_dxdz, sizeof(float)) );
	//-------------------------------------------
}

extern "C"
{
/* ---------------------------------------------------------------------- */
bool finalizeGpu(float *uu_x, float *uu_y, float *uu_z, float *lnrho){
/* Frees memory allocated on GPU.
*/
	//----------------------------------------------------------
	// Load from device memory back into host for saving the final snapshot (TODO: Might need changes after writing async transfers )
	//----------------------------------------------------------
	checkErr( cudaMemcpy(lnrho, d_lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	checkErr( cudaMemcpy(uu_x,  d_uu_x,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	checkErr( cudaMemcpy(uu_y,  d_uu_y,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	checkErr( cudaMemcpy(uu_z,  d_uu_z,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//----------------------------------------------------------
	//----------------------------------------------------------
	// Save the final snapshot. 
	//----------------------------------------------------------
	/*save_grid_information(t); // Save grid information for the most current .dat file TODO: test. 
	save_grid_data(lnrho, DATA_LNRHO_PATH);
	save_grid_data(uu_x, DATA_UU_X_PATH);
	save_grid_data(uu_y, DATA_UU_Y_PATH);
	save_grid_data(uu_z, DATA_UU_Z_PATH);*/
	//----------------------------------------------------------

	//Destroy timers
	//cudaEventDestroy( start );
	//cudaEventDestroy( stop );

	//Free device memory
	checkErr( cudaFree(d_lnrho) );
	checkErr( cudaFree(d_uu_x) );
	checkErr( cudaFree(d_uu_y) );
	checkErr( cudaFree(d_uu_z) );

	//Free diagnostic helper variables/arrays
	checkErr( cudaFree(d_umax) ); checkErr( cudaFree(d_umin) );
	checkErr( cudaFree(d_urms) );
 	checkErr( cudaFree(d_uxrms) ); checkErr( cudaFree(d_uyrms) ); checkErr( cudaFree(d_uzrms) );
	checkErr( cudaFree(d_rhorms) );
	checkErr( cudaFree(d_rhomax) ); checkErr( cudaFree(d_rhomin) );
	checkErr( cudaFree(d_uxmax) ); checkErr( cudaFree(d_uymax) ); checkErr( cudaFree(d_uzmax) );
	checkErr( cudaFree(d_uxmin) ); checkErr( cudaFree(d_uymin) ); checkErr( cudaFree(d_uzmin) );
	checkErr( cudaFree(d_partial_result) );
	checkErr( cudaFree(d_halo) );
	/*//Free pinned memory
	checkErr( cudaFreeHost(slice_lnrho) );
	checkErr( cudaFreeHost(slice_uu) );
	checkErr( cudaFreeHost(slice_uu_x) );
	checkErr( cudaFreeHost(slice_uu_y) );
	checkErr( cudaFreeHost(slice_uu_z) );*/

	


	//Reset device
	cudaDeviceReset();


	return EXIT_SUCCESS;
}

void RKintegration(float *uu_x, float *uu_y, float *uu_z, float *lnrho, int mx, int my, int mz, int nghost, int isubstep){
	
	//need to make those calls asynchronize
	copyouterhalostodevice(lnrho, d_lnrho, halo, d_halo, mx, my, mz, nghost);
	copyouterhalostodevice(uu_x, d_uu_x, halo, d_halo, mx, my, mz, nghost);
	copyouterhalostodevice(uu_y, d_uu_y, halo, d_halo, mx, my, mz, nghost);
	copyouterhalostodevice(uu_z, d_uu_z, halo, d_halo, mx, my, mz, nghost);

	rungekutta2N_cuda(d_lnrho, d_uu_x, d_uu_y, d_uu_z, d_w_lnrho, d_w_uu_x, d_w_uu_y, d_w_uu_z, d_lnrho_dest, d_uu_x_dest, d_uu_y_dest, d_uu_z_dest);

	copyinternalhalostohost(lnrho, d_lnrho, halo, d_halo, mx, my, mz, nghost);
	copyinternalhalostohost(uu_x, d_uu_x, halo, d_halo, mx, my, mz, nghost);
	copyinternalhalostohost(uu_y, d_uu_y, halo, d_halo, mx, my, mz, nghost);
	copyinternalhalostohost(uu_z, d_uu_z, halo, d_halo, mx, my, mz, nghost);

	//Swap array pointers
	swap_ptrs(&d_lnrho, &d_lnrho_dest);
	swap_ptrs(&d_uu_x, &d_uu_x_dest);
	swap_ptrs(&d_uu_y, &d_uu_y_dest);
	swap_ptrs(&d_uu_z, &d_uu_z_dest);
	return;
}

//void intitializeGPU(float *uu_x, float *uu_y, float *uu_z, float *lnrho, int nx, int ny, int nz, int nghost, float *x, float *y, float *z, float NU_VISC, float cs2_sound){ 
void intitializeGPU(float *uu_x, float *uu_y, float *uu_z, float *lnrho, int nx, int ny, int nz, int nghost, float *x, float *y, float *z, float nu, float cs2){ 
		// nx = mx, ny = my, nz = mz halo_depth = nghost
		int device;
		int halo_depth = nghost;
		cudaGetDevice(&device);
		printf("Using device %d\n", device);
		//cudaSetDevice(device); //Not yet enabled

		//Ensure that we're using a clean device
		cudaDeviceReset();
	
		//----------------------------------------------------------
		// Initialize global host variables
		//----------------------------------------------------------
		print_init_config();
		print_run_config();
		print_additional_defines();
		//----------------------------------------------------------

		//----------------------------------------------------------
		// Allocate host memory
		//----------------------------------------------------------
 
		/*float *lnrho; //Log density
		float *uu_x, *uu_y, *uu_z; //velocities

		lnrho = (float*) malloc(sizeof(float)*GRID_SIZE);
		uu_x  = (float*) malloc(sizeof(float)*GRID_SIZE);
		uu_y  = (float*) malloc(sizeof(float)*GRID_SIZE);
		uu_z  = (float*) malloc(sizeof(float)*GRID_SIZE);*/

		//Format the grids into 0.0 values to avoid potential noise in memory
		/*set_grids_zero(lnrho, uu_x, uu_y, uu_z);
		printf("Initializing grid to zero successful!");*/

		//----------------------------------------------------------
		// Allocating arrays for halos
		//----------------------------------------------------------

		halo_size = (halo_depth*nx*2 + halo_depth*(ny-halo_depth*2)*2)*(nz-halo_depth*2) + nx*ny*(halo_depth*2);
		d_lnrho_size = nx*ny*nz;
		
		halo = (float*) malloc(sizeof(float)*halo_size);
		checkErr(cudaMalloc ((void **) &d_halo, sizeof(float)*halo_size));
		//note: int GRID_SIZE ?
		//note: W_GRID_SIZE = ?
		//----------------------------------------------------------
		// Allocate device memory
		//----------------------------------------------------------

		checkErr( cudaMalloc(&d_lnrho, sizeof(float)*GRID_SIZE) );
		checkErr( cudaMalloc(&d_uu_x, sizeof(float)*GRID_SIZE) );
		checkErr( cudaMalloc(&d_uu_y, sizeof(float)*GRID_SIZE) );
		checkErr( cudaMalloc(&d_uu_z, sizeof(float)*GRID_SIZE) );


		checkErr( cudaMalloc(&d_w_lnrho, sizeof(float)*W_GRID_SIZE) );
		checkErr( cudaMalloc(&d_w_uu_x, sizeof(float)*W_GRID_SIZE) );
		checkErr( cudaMalloc(&d_w_uu_y, sizeof(float)*W_GRID_SIZE) );
		checkErr( cudaMalloc(&d_w_uu_z, sizeof(float)*W_GRID_SIZE) );

		//Temporary arrays

		checkErr( cudaMalloc(&d_lnrho_dest, sizeof(float)*GRID_SIZE) );
		checkErr( cudaMalloc(&d_uu_x_dest, sizeof(float)*GRID_SIZE) );
		checkErr( cudaMalloc(&d_uu_y_dest, sizeof(float)*GRID_SIZE) );
		checkErr( cudaMalloc(&d_uu_z_dest, sizeof(float)*GRID_SIZE) );

		

		checkErr( cudaMalloc((float**) &d_umax, sizeof(float)) );   //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_umin, sizeof(float)) );   //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_urms, sizeof(float)) );   //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_uxrms, sizeof(float)) );  //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_uyrms, sizeof(float)) );  //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_uzrms, sizeof(float)) );  //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_rhorms, sizeof(float)) ); //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_rhomax, sizeof(float)) ); //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_rhomin, sizeof(float)) ); //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_uxmax, sizeof(float)) );   //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_uxmin, sizeof(float)) );   //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_uymax, sizeof(float)) );   //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_uymin, sizeof(float)) );   //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_uzmax, sizeof(float)) );   //TODO this somewhere else
		checkErr( cudaMalloc((float**) &d_uzmin, sizeof(float)) );   //TODO this somewhere else
		printf("Device mem allocated: %f MiB\n", (4*sizeof(float)*GRID_SIZE + 4*sizeof(float)*W_GRID_SIZE)/powf(2,20));
		printf("Main array (d_lnrho) dims: (%d,%d,%d)\ntemporary result array dims (d_w_lnrho etc)(%d,%d,%d)\n", NX,NY,NZ, 										COMP_DOMAIN_SIZE_X,COMP_DOMAIN_SIZE_Y,COMP_DOMAIN_SIZE_Z);
		
		//----------------------------------------------------------
		// Load data into device memory
		//----------------------------------------------------------
		checkErr( cudaMemcpy(d_lnrho, lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
		checkErr( cudaMemcpy(d_uu_x, uu_x, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
		checkErr( cudaMemcpy(d_uu_y, uu_y, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
		checkErr( cudaMemcpy(d_uu_z, uu_z, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
		//Init also the dest arrays to avoid roaming NaN values
		checkErr( cudaMemcpy(d_lnrho_dest, lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
		checkErr( cudaMemcpy(d_uu_x_dest, uu_x, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
		checkErr( cudaMemcpy(d_uu_y_dest, uu_y, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
		checkErr( cudaMemcpy(d_uu_z_dest, uu_z, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) ); 
		//----------------------------------------------------------
		//----------------------------------------------------------
		//Load constants into device memory
		//----------------------------------------------------------
		load_dconsts(nu, cs2);
		//----------------------------------------------------------
return;
}
}
/* ---------------------------------------------------------------------- */
