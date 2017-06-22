
#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "timestep.cuh"
#include "collectiveops.cuh"
#include "../diagnostics_c.h"
#include "ddiagsextern.cuh"
#include "dfdfextern.cuh"
#include "../hydro_c.h"
#include "../cparam_c.h"
#include "defines_dims_PC.h"
#include "../eos_c.h"
#include "../cdata_c.h"
#include "defines_PC.h"
#include "../diagnostics_c.h"

extern float nu, cs2;

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

extern int *p_diags_hydro[];
extern const int n_diags_hydro;

//using namespace PC;

//----------------------------------------------------------
//  Calculate Courant timestep for the system. 
//----------------------------------------------------------
// Takes the domain's device pointers as parameters and return a new timestep
// d_umax = a single float containing the max velocity of the computational domain (allocated in compute.cu)
// d_partial_result = temporary array containing the partial results of the reduction (allocated in compute.cu)

float timestep_cuda(float* d_umax, float* d_partial_result, float* d_uu_x, float* d_uu_y, float* d_uu_z)
{
        //MV: It is better to calculate dt within the CPU after we get umax from max_vec_cuda, 
        //MV: because we need the information in the CPU lever anyway         

	// Determine the correct time step for the system
	static float dt, umax, uu_dt, visc_dt; //Initialize only once (static var lifetime == entire program)

	//Get uu max to d_umax
	max_vec_cuda(d_umax, d_partial_result, d_uu_x, d_uu_y, d_uu_z);
	cudaDeviceSynchronize();
	cudaMemcpy(&umax, (float*)d_umax, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();

	printf("UMAX: %F\n", umax);

	//Get dsmin (cannot be defined in defines.h)
	static float dsmin = DX;
	if (dsmin > DY)
		dsmin = DY;
	if (dsmin > DZ)	
		dsmin = DZ;

	//Courant timesteps
	//MV: DSMIN not yet defined in the constants.
	uu_dt = CDT*(dsmin/(umax + CS_SOUND)); 
	//MV: This stays actually conctant now, but if we add phenomena 
	//MV: like hyperviscosty and shock viscosity this could change 
	visc_dt = CDTV*dsmin*dsmin / NU_VISC; //TODO inverse NU_VISC in defines.h

	if (uu_dt < visc_dt) {
		dt = uu_dt;
	} else {
		dt = visc_dt;
	}
	return dt;
}
void get_maxscal_from_device(float &maxscal,float *d_src,float *d_max)
{
        max_scal_cuda(&maxscal, d_partial_result, d_src);
        cudaDeviceSynchronize();
        cudaMemcpy(&maxscal, (float*)d_max, sizeof(float), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
}
void get_minscal_from_device(float &minscal,float *d_src,float *d_min)
{
        min_scal_cuda(&minscal, d_partial_result, d_src);
        cudaDeviceSynchronize();
        cudaMemcpy(&minscal, (float*)d_min, sizeof(float), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
}
void max_advec()
{
	float uxmax, uymax, uzmax, maxadvec_;
	get_maxscal_from_device(uxmax,d_uu_x,d_uxmax);
	get_maxscal_from_device(uymax,d_uu_y,d_uymax);
	get_maxscal_from_device(uzmax,d_uu_z,d_uzmax);
      
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
	for (int i=0;i<nx;i++) maxadvec[i]=max(maxadvec[i],maxadvec_);
printf("maxadvec_= %f", maxadvec_);
}
void max_diffus()
{       
	for (int i=0;i<nx;i++) maxdiffus[i]=max(maxdiffus[i],nu*dxyz_2[i]);
}

/*void timeseries_diagnostics_cuda(float* d_umax, float* d_umin, float* d_urms, float* d_uxrms, 
                                 float* d_uyrms, float* d_uzrms, float* d_rhorms, 
                                 float* d_rhomax, float* d_uxmax, float* d_uymax, float* d_uzmax, 
                                 float* d_rhomin, float* d_uxmin, float* d_uymin, float* d_uzmin,*/ 
void timeseries_diagnostics_cuda(int step, float dt, double t)
{
	//Calculate, print and save all of the diagnostic variables calculated within the CUDA devices. 

	static float urms, umax, umin; 
	static float uxmax, uymax, uzmax, uxmin, uymin, uzmin; 
	static float rhomax, rhomin, rhorms;
	static float uxrms, uyrms, uzrms;
	
	//Get uu max from d_umax
	max_vec_cuda(d_umax, d_partial_result, d_uu_x, d_uu_y, d_uu_z);
	cudaDeviceSynchronize();
	cudaMemcpy(&umax, (float*)d_umax, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();
printf("umax= %f\n", umax);
        __diagnostics_MOD_save_name(&umax,p_diags_hydro[idiag_umax]);

	//Get uu_x max from d_uxmax 
        get_maxscal_from_device(uxmax,d_uu_x,d_uxmax);
        __diagnostics_MOD_save_name(&uxmax,p_diags_hydro[idiag_uxmax]);

	//Get uu_y max from d_uymax
        get_maxscal_from_device(uymax,d_uu_y,d_uymax);

	//Get uu_z max from d_uzmax 
        get_maxscal_from_device(uzmax,d_uu_z,d_uzmax);

	//Get rho max from d_lnrho
        get_maxscal_from_device(rhomax,d_lnrho,d_rhomax);
	if (!ldensity_nolog) rhomax = exp(rhomax);       //Change away from the logarithmic form

	//Get uu min from d_umin 
	min_vec_cuda(d_umin, d_partial_result, d_uu_x, d_uu_y, d_uu_z);
	cudaDeviceSynchronize();
	cudaMemcpy(&umin, (float*)d_umin, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();

	//Get uu_x min from d_uxmin 
        get_minscal_from_device(uxmin,d_uu_x,d_uxmin);

	//Get uu_y min from d_uymin 
        get_minscal_from_device(uymin,d_uu_y,d_uymin);

	//Get uu_z min from d_uzmin 
        get_minscal_from_device(uzmin,d_uu_z,d_uzmin);

	//Get lnrho min from d_lnrho 
        get_minscal_from_device(rhomin,d_lnrho,d_rhomin);
	if (!ldensity_nolog) rhomin = exp(rhomin);       //Change away from the logarithmic form

	//Get uu rms from d_umax
	vec_rms_cuda(d_urms, d_partial_result, d_uu_x, d_uu_y, d_uu_z);
	cudaDeviceSynchronize();
	cudaMemcpy(&urms, (float*)d_urms, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();

	//Get uu_x rms from d_uxrms
	scal_rms_cuda(d_uxrms, d_partial_result, d_uu_x);
	cudaDeviceSynchronize();
	cudaMemcpy(&uxrms, (float*)d_uxrms, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();

	//Get uu_y rms from d_uyrms
	scal_rms_cuda(d_uyrms, d_partial_result, d_uu_y);
	cudaDeviceSynchronize();
	cudaMemcpy(&uyrms, (float*)d_uyrms, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();

	//Get uu_z rms from d_uzrms
	scal_rms_cuda(d_uzrms, d_partial_result, d_uu_z);
	cudaDeviceSynchronize();
	cudaMemcpy(&uzrms, (float*)d_uzrms, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();

	//Get rho rms from d_lnrhorms
	scal_exp_rms_cuda(d_rhorms, d_partial_result, d_lnrho);
	cudaDeviceSynchronize();
	cudaMemcpy(&rhorms, (float*)d_rhorms, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize(); 
        if (iproc==0) 
        printf(" step = %i; t = %e; dt = %e; umax = %e; umin = %e; urms = %e",step, t, dt, umax, umin, urms);
	//printf(" step = %i; t = %e; dt = %e; umax = %e; umin = %e; urms = %e; \n uxrms = %e; uyrms = %e; uzrms = %e; \n uxmax = %e; uymax = %e; uzmax = %e; \n uxmin = %e; uymin = %e; uzmin = %e; \n rhomax = %e; rhomin = %e; rhorms = %e \n", 
        //        step, t, dt, umax, umin, urms, uxrms, uyrms, uzrms, uxmax, uymax, uzmax, uxmin, uymin, uzmin, rhomax, rhomin, rhorms);
	
	//Save the step into a file 
	//save_ts(t, dt, step, urms, uxrms, uyrms, uzrms, uxmax, uymax, uzmax, rhorms, umax, rhomax, uxmin, uymin, uzmin, rhomin, umin);
}


