//                             gpu_astaroth.cu
//                             ---------------

/* Functions for initializing, finalizing and performing of one integration substep with ASTAROTH-nucleus,
   to be called from PencilCode.

   Comments: 
   DATE March 17, 2017: 
   Omer Anjum: Added description of the functions
*/

//C libraries
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

//Headers
#include "dconsts.cuh"
#include "integrators.cuh"
#include "timestep.cuh"
#include "../cparam_c.h"
#include "smem.cuh"
#include "../cdata_c.h"
#include "../eos_c.h"
#include "../hydro_c.h"
#include "../viscosity_c.h"
#include "../forcing_c.h"
#include "defines_PC.h"
//#include "copyhalos.cuh"
#include "copyHalosConcur.cuh"

//DEBUG
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

int halo_size;
float *halo; 
float *d_halo;

float *d_lnrho, *d_uu_x, *d_uu_y, *d_uu_z;
float *d_w_lnrho, *d_w_uu_x, *d_w_uu_y, *d_w_uu_z;
float *d_lnrho_dest, *d_uu_x_dest, *d_uu_y_dest, *d_uu_z_dest;

// Device pointer for diagnostic quantities

float *d_umax, *d_umin; 
float *d_urms; 
float *d_uxrms, *d_uyrms, *d_uzrms; 
float *d_rhorms, *d_rhomax, *d_rhomin; 
float *d_uxmax, *d_uymax, *d_uzmax; 
float *d_uxmin, *d_uymin, *d_uzmin; 
float *d_partial_result;                    //Device pointer for partial result for the reductions

// Parameter of ohysics modules.

float nu, cs2, force, tforce_stop;

const int idiag_urms=0,
          idiag_uxrms=1,
          idiag_uzrms=2,
          idiag_umax=3,
          idiag_uxmin=4,
          idiag_uymin=5,
          idiag_uzmin=6,
          idiag_uxmax=7,
          idiag_uymax=8,
          idiag_uzmax=9;

const int ndiags_hydro=10;
int *p_diags_hydro[ndiags_hydro];

const int npars_visc=1;
float *p_pars_visc[npars_visc];

const int npars_eos=1;
float *p_pars_eos[npars_eos];

const int npars_force=2;
float *p_pars_force[npars_force];

/***********************************************************************************************/
inline void swap_ptrs(float** a, float** b)
{
//  Swaps pointers a,b (do xor swap if too slow)

	float* temp = *a;
	*a = *b;
	*b = temp;
}
/***********************************************************************************************/
//using namespace PC;

void load_dconsts()
{
//  Loads constants into device memory

	checkErr( cudaMemcpyToSymbol(d_NX, &NX, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_NY, &NY, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_NZ, &NZ, sizeof(int)) );

        const int pad_size=PAD_SIZE;
	checkErr( cudaMemcpyToSymbol(d_PAD_SIZE, &pad_size, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_BOUND_SIZE, &BOUND_SIZE, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_COMP_DOMAIN_SIZE_X, &COMP_DOMAIN_SIZE_X, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_COMP_DOMAIN_SIZE_Y, &COMP_DOMAIN_SIZE_Y, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_COMP_DOMAIN_SIZE_Z, &COMP_DOMAIN_SIZE_Z, sizeof(int)) );

        const float nelements_float = W_GRID_SIZE;
	checkErr( cudaMemcpyToSymbol(d_NELEMENTS_FLOAT, &nelements_float, sizeof(float)) );

	checkErr( cudaMemcpyToSymbol(d_DOMAIN_SIZE_X, &DOMAIN_SIZE_X, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DOMAIN_SIZE_Y, &DOMAIN_SIZE_Y, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DOMAIN_SIZE_Z, &DOMAIN_SIZE_Z, sizeof(float)) );

	const int h_w_grid_y_offset = COMP_DOMAIN_SIZE_X;
	const int h_w_grid_z_offset = COMP_DOMAIN_SIZE_X*COMP_DOMAIN_SIZE_Y;
	const int h_grid_y_offset = NX;
	const int h_grid_z_offset = NX*NY;

	checkErr( cudaMemcpyToSymbol(d_W_GRID_Y_OFFSET, &h_w_grid_y_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_W_GRID_Z_OFFSET, &h_w_grid_z_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_GRID_Y_OFFSET, &h_grid_y_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_GRID_Z_OFFSET, &h_grid_z_offset, sizeof(int)) );

	//------Computational domain's bottom and top indices---------

	checkErr( cudaMemcpyToSymbol(d_CX_TOP, &CX_TOP, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CY_TOP, &CY_TOP, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CZ_TOP, &CZ_TOP, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_CX_BOT, &CX_BOT, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CY_BOT, &CY_BOT, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CZ_BOT, &CZ_BOT, sizeof(int)) );

	//-------------Real distance between grid points---------------
	checkErr( cudaMemcpyToSymbol(d_DX, &DX, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DY, &DY, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DZ, &DZ, sizeof(float)) );

	//----------Location of the grid coordinate origin-------------
        checkErr( cudaMemcpyToSymbol(d_XORIG, &XORIG, sizeof(float)) );
        checkErr( cudaMemcpyToSymbol(d_YORIG, &YORIG, sizeof(float)) );
        checkErr( cudaMemcpyToSymbol(d_ZORIG, &ZORIG, sizeof(float)) );

	//----------Shearing parameters---------------
	const int interp_order = INTERP_ORDER;
        checkErr( cudaMemcpyToSymbol(d_INTERP_ORDER, &interp_order, sizeof(int)) );

        checkErr( cudaMemcpyToSymbol(d_Q_SHEAR, &Q_SHEAR, sizeof(float)) );
        checkErr( cudaMemcpyToSymbol(d_OMEGA, &OMEGA, sizeof(float)) );

	//------------------Optional physics switches------------------------
        checkErr( cudaMemcpyToSymbol(d_LFORCING, &LFORCING, sizeof(int)) );
        checkErr( cudaMemcpyToSymbol(d_LSHEAR, &LSHEAR, sizeof(int)) );

        const int lcoriolis = LCORIOLIS;
        checkErr( cudaMemcpyToSymbol(d_LCORIOLIS, &lcoriolis, sizeof(int)) );

	//-----------Coefficients of Runge-Kutta method-------------
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

	checkErr( cudaMemcpyToSymbol(d_NU_VISC, &NU_VISC, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_CS2_SOUND, &CS2_SOUND, sizeof(float)) );	

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

	const float diffmn_dxdy = 1.0/(720.0*DX*DY); 
	const float diffmn_dydz = 1.0/(720.0*DY*DZ);
	const float diffmn_dxdz = 1.0/(720.0*DZ*DX);
	
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
}
/***********************************************************************************************/
extern "C" void initializeGPU(){ 

	int device;
	cudaGetDevice(&device);
	//cudaSetDevice(device); //Not yet enabled

	//Ensure that we're using a clean device
	cudaDeviceReset();

/*if (iproc==0){
printf("mx,my,mz,nx,ny,nz,nghost= %d %d %d %d %d %d %d \n", mx,my,mz,nx,ny,nz,nghost);
printf("nxgrid,nygrid,nzgrid= %d %d %d \n", nxgrid,nygrid,nzgrid);
printf("l1,l2,m1,m2,n1,n2= %d %d %d %d %d %d \n", l1,l2,m1,m2,n1,n2);
printf("xyz0, xyz1 %f %f %f %f %f %f \n", xyz0[0], xyz0[1], xyz0[2], xyz1[0], xyz1[1], xyz1[2]); 
printf("Lxyz %f %f %f \n", lxyz[0], lxyz[1], lxyz[2]); 
printf(lcartesian_coords ? "CARTESIAN \n" : "NONCARTESIAN");
}
printf("[xyz]minmax %f %f %f %f %f %f \n", x[l1-1], x[l2-1], y[m1-1], y[m2-1], z[n1-1], z[n2-1]); 
printf("[xyz]minmax_ghost %f %f %f %f %f %f \n", x[0], x[mx-1], y[0], y[my-1], z[0], z[mz-1]); */

	// Allocating arrays for halos

	halo_size = (nghost*nx*2 + nghost*(ny-nghost*2)*2)*(nz-nghost*2) + nx*ny*(nghost*2);
	halo = (float*) malloc(sizeof(float)*halo_size);
	checkErr(cudaMalloc((void **) &d_halo, sizeof(float)*halo_size));

	// Allocate device memory

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

        // Diagnostics quantities. TODO this somewhere else?

	checkErr( cudaMalloc((float**) &d_umax, sizeof(float)) );   
	checkErr( cudaMalloc((float**) &d_umin, sizeof(float)) );   
	checkErr( cudaMalloc((float**) &d_urms, sizeof(float)) );   
	checkErr( cudaMalloc((float**) &d_uxrms, sizeof(float)) );  
	checkErr( cudaMalloc((float**) &d_uyrms, sizeof(float)) );  
	checkErr( cudaMalloc((float**) &d_uzrms, sizeof(float)) );  
	checkErr( cudaMalloc((float**) &d_rhorms, sizeof(float)) ); 
	checkErr( cudaMalloc((float**) &d_rhomax, sizeof(float)) ); 
	checkErr( cudaMalloc((float**) &d_rhomin, sizeof(float)) ); 
	checkErr( cudaMalloc((float**) &d_uxmax, sizeof(float)) );   
	checkErr( cudaMalloc((float**) &d_uxmin, sizeof(float)) );   
	checkErr( cudaMalloc((float**) &d_uymax, sizeof(float)) );   
	checkErr( cudaMalloc((float**) &d_uymin, sizeof(float)) );   
	checkErr( cudaMalloc((float**) &d_uzmax, sizeof(float)) );   
	checkErr( cudaMalloc((float**) &d_uzmin, sizeof(float)) );   
	checkErr( cudaMalloc((float**) &d_partial_result, sizeof(float)) );   

        if (iproc==0)
	{
	  printf(" Device mem allocated: %f MiB\n", (4*sizeof(float)*GRID_SIZE + 4*sizeof(float)*W_GRID_SIZE)/powf(2,20));
		  //printf("Main array (d_lnrho etc) dims: (%d,%d,%d)\ntemporary result array dims (d_w_lnrho etc)(%d,%d,%d)\n",
                  //       NX,NY,NZ,COMP_DOMAIN_SIZE_X,COMP_DOMAIN_SIZE_Y,COMP_DOMAIN_SIZE_Z);
        }
        // Get private data from physics modules.
        
 	if (lhydro){
        	hydro_push2c(p_diags_hydro); 
	}
        if (lviscosity){
        	viscosity_push2c(p_pars_visc);
        	nu=*p_pars_visc[0];
	}
	if (lforcing){
        	forcing_push2c(p_pars_force);
        	force=*p_pars_force[0];
		tforce_stop=*p_pars_force[1];
	}

        eos_push2c(p_pars_eos);
        cs2=*p_pars_eos[0];

        /*if (iproc==0){	
	print_init_config();
	print_run_config();
	print_additional_defines();
        }

        if (iproc==0) {
        printf("nu %f \n", nu);
        printf("idiag_urms= %d \n", *p_diags_hydro[idiag_urms]);
        printf("idiag_uxrms= %d \n", *p_diags_hydro[idiag_uxrms]);
        printf("idiag_uzrms= %d \n", *p_diags_hydro[idiag_uzrms]);
        printf("idiag_umax= %d \n", *p_diags_hydro[idiag_umax]);
        printf("idiag_uxmin= %d \n", *p_diags_hydro[idiag_uxmin]);
        printf("idiag_uymin= %d \n", *p_diags_hydro[idiag_uymin]);
        printf("idiag_uzmin= %d \n", *p_diags_hydro[idiag_uzmin]);
        printf("idiag_uxmax= %d \n", *p_diags_hydro[idiag_uxmax]);
        printf("idiag_uymax= %d \n", *p_diags_hydro[idiag_uymax]);
        printf("idiag_uzmax= %d \n", *p_diags_hydro[idiag_uzmax]);
	}*/

	// Load constants into device memory.
	load_dconsts();

        initializeCopying();
}
/***********************************************************************************************/
extern "C" void substepGPU(float *uu_x, float *uu_y, float *uu_z, float *lnrho, int isubstep, bool full_inner=false, bool full=false){
	//need to make those calls asynchronize
	/*copyouterhalostodevice(lnrho, d_lnrho, halo, d_halo, mx, my, mz, nghost);
	copyouterhalostodevice(uu_x, d_uu_x, halo, d_halo, mx, my, mz, nghost);
	copyouterhalostodevice(uu_y, d_uu_y, halo, d_halo, mx, my, mz, nghost);
	copyouterhalostodevice(uu_z, d_uu_z, halo, d_halo, mx, my, mz, nghost);*/

printf(full ? "full\n" : "not full\n");
        if (full) 
	{
          copyAll(lnrho, d_lnrho);
          copyAll(uu_x, d_uu_x);
          copyAll(uu_y, d_uu_y);
          copyAll(uu_z, d_uu_z);
   	}
	else
	{
          copyOuterHalos(lnrho, d_lnrho);
          copyOuterHalos(uu_x, d_uu_x);
          copyOuterHalos(uu_y, d_uu_y);
          copyOuterHalos(uu_z, d_uu_z);
	}

	rungekutta2N_cuda(d_lnrho, d_uu_x, d_uu_y, d_uu_z, d_w_lnrho, d_w_uu_x, d_w_uu_y, d_w_uu_z,
                          d_lnrho_dest, d_uu_x_dest, d_uu_y_dest, d_uu_z_dest, isubstep);

	//Swap array pointers
	swap_ptrs(&d_lnrho, &d_lnrho_dest);
	swap_ptrs(&d_uu_x, &d_uu_x_dest);
	swap_ptrs(&d_uu_y, &d_uu_y_dest);
	swap_ptrs(&d_uu_z, &d_uu_z_dest);

	/*copyinternalhalostohost(lnrho, d_lnrho, halo, d_halo, mx, my, mz, nghost);
	copyinternalhalostohost(uu_x, d_uu_x, halo, d_halo, mx, my, mz, nghost);
	copyinternalhalostohost(uu_y, d_uu_y, halo, d_halo, mx, my, mz, nghost);
	copyinternalhalostohost(uu_z, d_uu_z, halo, d_halo, mx, my, mz, nghost);*/
printf(full_inner ? "full_inner\n" : "not full_inner\n");
        if (full_inner) 
	{
        	copyInnerAll(lnrho, d_lnrho);
        	copyInnerAll(uu_x, d_uu_x);
        	copyInnerAll(uu_y, d_uu_y);
        	copyInnerAll(uu_z, d_uu_z);
	}
	else
	{
        	copyInnerHalos(lnrho, d_lnrho);
        	copyInnerHalos(uu_x, d_uu_x);
        	copyInnerHalos(uu_y, d_uu_y);
        	copyInnerHalos(uu_z, d_uu_z);
	}
printf(ldiagnos ? "out\n" : "not out\n");
        if (ldiagnos) {
                timeseries_diagnostics_cuda(it, dt, t);
/*d_umax, d_umin, d_urms, d_uxrms, d_uyrms, d_uzrms, d_rhorms,
                                            d_rhomax, d_uxmax, d_uymax, d_uzmax,
                                            d_rhomin, d_uxmin, d_uymin, d_uzmin,*/
        }
        if (lfirst && ldt) {
		max_diffus(); max_advec();
	}
}
/***********************************************************************************************/
extern "C" void finalizeGPU()
{
// Frees memory allocated on GPU.

        //Destroy timers
        //cudaEventDestroy( start );
        //cudaEventDestroy( stop );

        //Free device memory of grids
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

        //Free pinned memory
        /*checkErr( cudaFreeHost(slice_lnrho) );
        checkErr( cudaFreeHost(slice_uu) );
        checkErr( cudaFreeHost(slice_uu_x) );
        checkErr( cudaFreeHost(slice_uu_y) );
        checkErr( cudaFreeHost(slice_uu_z) );*/

        finalizeCopying();
        cudaDeviceReset();
}
/***********************************************************************************************/
