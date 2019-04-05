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
#define EXTERN
#include "dconsts.cuh"
#include "integrators.cuh"
#include "timestep.cuh"
#include "../cparam_c.h"
#include "smem.cuh"
#include "../cdata_c.h"
#include "../density_c.h"
#include "../eos_c.h"
#include "../hydro_c.h"
#include "../viscosity_c.h"
#include "../forcing_c.h"
#include "../sub_c.h"
#include "defines_PC.h"
#include "copyhalos.cuh"
//#include "copyHalosConcur.cuh"

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
float *output;

#include "dfdf.cuh"

// Device pointer for diagnostic quantities

float *d_umax, *d_umin; 
float *d_urms, *d_rhorms; 
float *d_uxrms, *d_uyrms, *d_uzrms; 
float *d_partial_result, *d_scaldiag;                    //Device pointer for partial result for the reductions

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
        const int cx_top = CX_TOP;
        const int cy_top = CY_TOP;
        const int cz_top = CZ_TOP;

	checkErr( cudaMemcpyToSymbol(d_W_GRID_Y_OFFSET, &h_w_grid_y_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_W_GRID_Z_OFFSET, &h_w_grid_z_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_GRID_Y_OFFSET, &h_grid_y_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_GRID_Z_OFFSET, &h_grid_z_offset, sizeof(int)) );

	//------Computational domain's bottom and top indices---------

	checkErr( cudaMemcpyToSymbol(d_CX_TOP, &cx_top, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CY_TOP, &cy_top, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CZ_TOP, &cz_top, sizeof(int)) );

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
extern "C" void initializeGPU(float *uu_x, float *uu_y, float *uu_z, float *lnrho){ 

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

	checkErr( cudaMalloc( &d_umax, sizeof(float)) );   
	checkErr( cudaMalloc( &d_umin, sizeof(float)) );   
	checkErr( cudaMalloc( &d_urms, sizeof(float)) );   
	checkErr( cudaMalloc( &d_uxrms, sizeof(float)) );  
	checkErr( cudaMalloc( &d_uyrms, sizeof(float)) );  
	checkErr( cudaMalloc( &d_uzrms, sizeof(float)) );  
	checkErr( cudaMalloc( &d_rhorms, sizeof(float)) ); 
	checkErr( cudaMalloc( &d_partial_result, sizeof(float)) );   
	checkErr( cudaMalloc( &d_scaldiag, sizeof(float)) );   

        if (iproc==0)
	{
	  printf("In gpu_astaroth.cu in initializeGPU Device mem allocated: %f MiB\n", (4*sizeof(float)*GRID_SIZE + 4*sizeof(float)*W_GRID_SIZE)/powf(2,20));
		  //printf("Main array (d_lnrho etc) dims: (%d,%d,%d)\ntemporary result array dims (d_w_lnrho etc)(%d,%d,%d)\n",
                  //       NX,NY,NZ,COMP_DOMAIN_SIZE_X,COMP_DOMAIN_SIZE_Y,COMP_DOMAIN_SIZE_Z);
        }
        // Get private data from physics modules.
        printf("Calling density_push2c \n");
 	if (ldensity){
        	density_push2c(p_diags_density); 
	}
	printf("Calling hydro_push2c \n");
 	if (lhydro){
        	hydro_push2c(p_diags_hydro); 
	}
	printf("Calling viscosity_push2c \n");
        if (lviscosity){
        	viscosity_push2c(p_pars_visc);
	}
	printf("Calling forcing_push2c \n");
	if (lforcing){
        	forcing_push2c(p_pars_force);
	}
        eos_push2c(p_pars_eos);

        /*if (iproc==0){	
	print_init_config();
	print_run_config();
	print_additional_defines();
        }

        if (iproc==0) {
        printf("nu %f \n", nu);
        printf("idiag_urms= %d \n", idiag_urms);
        printf("idiag_uxrms= %d \n", idiag_uxrms);
        printf("idiag_uzrms= %d \n", idiag_uzrms);
        printf("idiag_umax= %d \n", idiag_umax);
        printf("idiag_uxmin= %d \n", idiag_uxmin);
        printf("idiag_uymin= %d \n", idiag_uymin);
        printf("idiag_uzmin= %d \n", idiag_uzmin);
        printf("idiag_uxmax= %d \n", idiag_uxmax);
        printf("idiag_uymax= %d \n", idiag_uymax);
        printf("idiag_uzmax= %d \n", idiag_uzmax);
	}*/

	// Load constants into device memory.
	load_dconsts();	
	
	// Allocating arrays for halos

	//halo_size = (nghost*nx*2 + nghost*(ny-nghost*2)*2)*(nz-nghost*2) + nx*ny*(nghost*2);
	printf("mx = %d, my = %d, mz = %d, nghost = %d", mx, my, mz, nghost);
	halo_size = (mx*my*mz) - (mx-2*nghost)*(my-2*nghost)*(mz-2*nghost);
	printf("in initializeGPU halo_size = %d\n",halo_size);
	halo = (float*) malloc(sizeof(float)*halo_size);
	checkErr(cudaMalloc((float**)&d_halo, sizeof(float)*halo_size));
	printf("Inside initializeGPU in gpu_astaroth_v2 &halo, &d_halo  %p %p pointing to  %p %p\n",&halo,&d_halo, halo, d_halo);
        //initializeCopying();
	printf("Stop: GPU initialized success inside gpu_astaroth_v2.cu\n");
}
/***********************************************************************************************/
float max_advec()
{
        float uxmax, uymax, uzmax=0, maxadvec_;
        get_maxscal_from_device(uxmax,d_uu_x);
        get_maxscal_from_device(uymax,d_uu_y);
        //get_maxscal_from_device(uzmax,d_uu_z);
//printf("UYMAX= %f \n", uymax);
//return;
        if (lmaximal_cdt) {
                maxadvec_=max(abs(uxmax)/dx,max(abs(uymax)/dy,abs(uzmax)/dz));
                /*advec_uu[ix]=max(abs(p%uu(:,1))*dline_1[0][ix],
                                         abs(p%uu(:,2))*dline_1[1][ix],
                                         abs(p%uu(:,3))*dline_1[2][ix]);*/
        }
        else
        {
                maxadvec_=(abs(uxmax)/dx+abs(uymax)/dy+abs(uzmax)/dz);
                /*advec_uu[ix]=abs(p%uu(:,1))*dline_1[0][ix]+
                         abs(p%uu(:,2))*dline_1[1][ix]+
                         abs(p%uu(:,3))*dline_1[2][ix]; */
        }
//printf("maxadvec_= %f \n", maxadvec_);
        return maxadvec_;
}
/***********************************************************************************************/
float max_diffus()
{
        float maxdiffus_;
        for (int i=0;i<nx;i++) maxdiffus_=max(maxdiffus_,nu*dxyz_2[i]);
        return maxdiffus_;
}
/***********************************************************************************************/
extern "C" void substepGPU(float *uu_x, float *uu_y, float *uu_z, float *lnrho, int isubstep, bool full_inner=false, bool full=false){
	
	//cudaSetDevice(0);
	printf("Stop: Now inside substepGPU\n");
	//need to make those calls asynchronize
	
	//printf(full ? "full\n" : "not full\n");
	//int offset=54273 + 134^2+134+60;
	//int offset = 134*134 +(134*3)+3;
	
	//printf("Inside initializeGPU in gpu_astaroth_v2 &halo, &d_halo  %p %p pointing to  %p %p\n",&halo,&d_halo, halo, d_halo);
	printf("Inside substepGPU in gpu_astaroth_v2 &halo, &d_halo  %p %p pointing to  %p %p\n",&halo,&d_halo, halo, d_halo);
	printf("Stop: Going to copy grid to GPU or copyOuterHalos inside gpu_astaroth.cu\n");
	/*if (iproc==0) {
		printf("uu_x= %f %f %f \n", *(uu_x+offset),*(uu_x+offset+1), *(uu_x+offset+2));
		printf("uu_y= %f %f %f \n", *(uu_y+offset),*(uu_y+offset+1), *(uu_y+offset+2));
		printf("uu_z= %f %f %f \n", *(uu_z+offset),*(uu_z+offset+1), *(uu_z+offset+2));
		printf("lnrho= %f %f %f \n", *(lnrho+offset),*(lnrho+offset+1), *(lnrho+offset+2));
	}*/
        if (full) 
	{
		 //----------------------------------------------------------
		// Load data into device memory
		//----------------------------------------------------------
		//halo_size = (nghost*nx*2 + nghost*(ny-nghost*2)*2)*(nz-nghost*2) + nx*ny*(nghost*2);
		//halo = (float*) malloc(sizeof(float)*halo_size);
		printf("Stop: Going to copy grid to GPU inside gpu_astaroth.cu\n");
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
	
   	}
	else
	{
		printf("Stop: Going to call copyouterhalostodevice inside gpu_astaroth.cu\n");
		copyouterhalostodevice(lnrho, d_lnrho, halo, d_halo, mx, my, mz, nghost);
		copyouterhalostodevice(uu_x, d_uu_x, halo, d_halo, mx, my, mz, nghost);
		copyouterhalostodevice(uu_y, d_uu_y, halo, d_halo, mx, my, mz, nghost);
		copyouterhalostodevice(uu_z, d_uu_z, halo, d_halo, mx, my, mz, nghost);
	}

	if (lfirst && ldt) {
		//float dt1_ = 0; // debug
                float dt1_advec  = max_advec()/cdt;
                float dt1_diffus = max_diffus()/cdtv;
                float dt1_=sqrt(pow(dt1_advec,2) + pow(dt1_diffus,2));
		//printf("stop1: dt1_ = %f, dt = %f\n", dt1_,dt);
                set_dt(dt1_);
		//printf("stop2: dt1_ = %f, dt = %f\n", dt1_,dt); 
		checkErr( cudaMemcpyToSymbol(d_DT, &dt, sizeof(float)));
        }
       	checkErr( cudaMemcpyToSymbol(d_DT, &dt, sizeof(float)));
	//printf("stop3: dt = %f in gpu_astaroth.cu\n", dt);
	printf("Calling rungekutta2N_cuda\n");
	//cudaSetDevice(0);

	rungekutta2N_cuda(d_lnrho, d_uu_x, d_uu_y, d_uu_z, d_w_lnrho, d_w_uu_x, d_w_uu_y, d_w_uu_z, d_lnrho_dest, d_uu_x_dest, d_uu_y_dest, d_uu_z_dest, isubstep);


	/*float tmp;
	const int mx = 128 + 6;
	const int my = mx;
	//const float mz = mx;
	const int bound_size = 3;
	//cudaDeviceSynchronize();
	cudaMemcpy(&tmp, &d_uu_x_dest[bound_size + bound_size*mx + bound_size*mx*my], sizeof(float), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();*/
	//printf("AAAAAAAAAAAAAA %e\n", tmp);
	
	//degbug starts
	/*if (1){
		checkErr(cudaMemcpy(output, d_uu_x_dest, sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost));
		//printf("output[54728] = %f\n", output[54728]);
		int idx;
		for(int j = 6; j < 9; j++){
			for(int i = 0; i < mx; i++){
				idx = (mx*mx)*3+mx*j+i;
				printf("uu_x[%d], output[idx] = %f %f\n",idx, uu_x[idx], output[idx]);
			}
		}
	}*/

		printf("Stop: Going to Swap Pointers\n");
		//Swap array pointers
		swap_ptrs(&d_lnrho, &d_lnrho_dest);
		swap_ptrs(&d_uu_x, &d_uu_x_dest);
		swap_ptrs(&d_uu_y, &d_uu_y_dest);
		swap_ptrs(&d_uu_z, &d_uu_z_dest);
		//printf("Stop: starting cudaMemcpyDeviceToHost in gpu_astaroth.cu \n");


	//printf("Stop: Starting Copying halos in gpu_astaroth.cu \n");

        if (full_inner) 
	{	
		printf("Stop: Inside full_inner to copy grid to host or inner halos to host\n");
		
        	checkErr( cudaMemcpy(lnrho, d_lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
		checkErr( cudaMemcpy(uu_x,  d_uu_x,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
		checkErr( cudaMemcpy(uu_y,  d_uu_y,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
		checkErr( cudaMemcpy(uu_z,  d_uu_z,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	}
	else
	{
		//copyinternalhalostohost(lnrho, d_lnrho, halo, d_halo, mx, my, mz, nghost);
		copyinternalhalostohost(uu_x, d_uu_x, halo, d_halo, mx, my, mz, nghost);
		//copyinternalhalostohost(uu_y, d_uu_y, halo, d_halo, mx, my, mz, nghost);
		//copyinternalhalostohost(uu_z, d_uu_z, halo, d_halo, mx, my, mz, nghost);
		printf("Stop: after copyinternalhalostohost\n");
	}

	//degbug starts
	
	/*if (1){
		checkErr(cudaMemcpy(output, d_uu_x, sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost));
		//printf("output[54728] = %f\n", output[54728]);
		int idx;
		for(int j = 3; j < 6; j++){
			for(int i = 3; i < 35; i++){
				idx = (mx*mx)*3+mx*j+i;
				printf("@idx=%d, uu_x[idx], output[idx] = %f %f\n",idx, uu_x[idx], output[idx]);
			}
		}
	}*/
	//printf("Stop: Finished Copying halos in gpu_astaroth.cu \n");
	//printf("Stop: Now going inside timeseries_diagnostics_cuda(it, dt, t) in gpu_astaroth.cu \n");
        if (ldiagnos) timeseries_diagnostics_cuda(it, dt, t);
	//printf("Stop: Finished executing timeseries_diagnostics_cuda(it, dt, t) in gpu_astaroth.cu \n");
}
/***********************************************************************************************/
extern "C" void finalizeGPU()
{
	
	//checkErr( cudaMemcpy(lnrho, d_lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//checkErr( cudaMemcpy(uu_x,  d_uu_x,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//checkErr( cudaMemcpy(uu_y,  d_uu_y,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//checkErr( cudaMemcpy(uu_z,  d_uu_z,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
        //Destroy timers
        //cudaEventDestroy( start );
        //cudaEventDestroy( stop );
        //Free device memory of grids
	printf("stop1: inside finalizeGPU in gpy_astaroth_v2.cu\n");
        checkErr( cudaFree(d_lnrho) );
	
//printf("vor helper, iproc= %d\n",iproc);
        checkErr( cudaFree(d_uu_x) );
        checkErr( cudaFree(d_uu_y) );
        checkErr( cudaFree(d_uu_z) );
	
        //Free diagnostic helper variables/arrays
        checkErr( cudaFree(d_umax) ); checkErr( cudaFree(d_umin) );
        checkErr( cudaFree(d_urms) );
        checkErr( cudaFree(d_uxrms) ); checkErr( cudaFree(d_uyrms) ); checkErr( cudaFree(d_uzrms) );
        checkErr( cudaFree(d_rhorms) );
        checkErr( cudaFree(d_partial_result) );
        checkErr( cudaFree(d_halo) );
	
        //Free pinned memory
        /*checkErr( cudaFreeHost(slice_lnrho) );
        checkErr( cudaFreeHost(slice_uu) );
        checkErr( cudaFreeHost(slice_uu_x) );
        checkErr( cudaFreeHost(slice_uu_y) );
        checkErr( cudaFreeHost(slice_uu_z) );*/

        //finalizeCopying();
	free(halo);
	printf("stop2: inside finalizeGPU in gpy_astaroth_v2.cu\n");
	cudaDeviceSynchronize();
        //checkErr(cudaDeviceReset());
	printf("GPU finalized %d", iproc);
}
/***********************************************************************************************/
