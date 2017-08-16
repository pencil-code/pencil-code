
#include <stdio.h>
#include <math.h>
#include <algorithm>

#define EXTERN extern
#include "timestep.cuh"
#include "collectiveops.cuh"
#include "../diagnostics_c.h"
#include "ddiagsextern.cuh"
#include "dfdf.cuh"
#include "../density_c.h"
#include "../hydro_c.h"
#include "../cparam_c.h"
#include "defines_dims_PC.h"
#include "../cdata_c.h"
#include "defines_PC.h"
#include "../viscosity_c.h"
#include "../eos_c.h"

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
void get_maxscal_from_device(float & maxscal,float *d_src)
{
	max_scal_cuda(d_scaldiag, d_partial_result, d_src);
       	cudaDeviceSynchronize();
       	cudaMemcpy(&maxscal, d_scaldiag, sizeof(float), cudaMemcpyDeviceToHost);
}
void get_minscal_from_device(float &minscal,float *d_src)
{
        min_scal_cuda(d_scaldiag, d_partial_result, d_src);
        cudaDeviceSynchronize();
        cudaMemcpy(&minscal, d_scaldiag, sizeof(float), cudaMemcpyDeviceToHost);
        //cudaDeviceSynchronize();
}

/*void timeseries_diagnostics_cuda(float* d_umax, float* d_umin, float* d_urms, float* d_uxrms, 
                                 float* d_uyrms, float* d_uzrms, float* d_rhorms, */
void timeseries_diagnostics_cuda(int step, float dt, double t)
{
	//Calculate, print and save all of the diagnostic variables calculated within the CUDA devices. 

	static float urms, umax, umin, diag; 
	static float rhorms;
	static float uxrms, uyrms, uzrms;
	
        if (idiag_uxmax>0) {
        	get_maxscal_from_device(diag,d_uu_x);
		save_name(diag,idiag_uxmax);
        }
        if (idiag_uymax>0) {
        	get_maxscal_from_device(diag,d_uu_y);
		save_name(diag,idiag_uymax);
        }
        if (idiag_uzmax>0) {
        	get_maxscal_from_device(diag,d_uu_z);
		save_name(diag,idiag_uzmax);
        }
        if (idiag_uxmin>0) {
        	get_minscal_from_device(diag,d_uu_x);
		save_name(diag,idiag_uxmin);
        }
        if (idiag_uymin>0) {
        	get_minscal_from_device(diag,d_uu_y);
		save_name(diag,idiag_uymin);
        }
        if (idiag_uzmin>0) {
        	get_minscal_from_device(diag,d_uu_z);
		save_name(diag,idiag_uzmin);
        }
        if (idiag_umax>0) {
		max_vec_cuda(d_umax, d_partial_result, d_uu_x, d_uu_y, d_uu_z);
		cudaDeviceSynchronize();
		cudaMemcpy(&umax, (float*)d_umax, sizeof(float), cudaMemcpyDeviceToHost); 
		cudaDeviceSynchronize();
printf("umax= %f\n", umax);
        	save_name(umax,idiag_umax);
        }
        if (idiag_rhomax>0) {
        	get_maxscal_from_device(diag,d_lnrho);
		if (!ldensity_nolog) diag = exp(diag);       //Change away from the logarithmic form
		save_name(diag,idiag_rhomax);
        }
        if (idiag_rhomin>0) {
        	get_minscal_from_device(diag,d_lnrho);
		if (!ldensity_nolog) diag = exp(diag);       //Change away from the logarithmic form
		save_name(diag,idiag_rhomin);
        }
        if (idiag_umin>0) {
		min_vec_cuda(d_umin, d_partial_result, d_uu_x, d_uu_y, d_uu_z);
		cudaDeviceSynchronize();
		cudaMemcpy(&umin, (float*)d_umin, sizeof(float), cudaMemcpyDeviceToHost); 
		cudaDeviceSynchronize();
	}
 	if (idiag_urms){
		vec_rms_cuda(d_urms, d_partial_result, d_uu_x, d_uu_y, d_uu_z);
		cudaDeviceSynchronize();
		cudaMemcpy(&urms, (float*)d_urms, sizeof(float), cudaMemcpyDeviceToHost); 
		cudaDeviceSynchronize();
	}
 	if (idiag_uxrms){
		scal_rms_cuda(d_uxrms, d_partial_result, d_uu_x);
		cudaDeviceSynchronize();
		cudaMemcpy(&uxrms, (float*)d_uxrms, sizeof(float), cudaMemcpyDeviceToHost); 
		cudaDeviceSynchronize();
	}
 	if (idiag_uyrms){
		scal_rms_cuda(d_uyrms, d_partial_result, d_uu_y);
		cudaDeviceSynchronize();
		cudaMemcpy(&uyrms, (float*)d_uyrms, sizeof(float), cudaMemcpyDeviceToHost); 
		cudaDeviceSynchronize();
	}
 	if (idiag_uzrms){
		scal_rms_cuda(d_uzrms, d_partial_result, d_uu_z);
		cudaDeviceSynchronize();
		cudaMemcpy(&uzrms, (float*)d_uzrms, sizeof(float), cudaMemcpyDeviceToHost); 
		cudaDeviceSynchronize();
	}
 	if (idiag_rhorms){
		scal_exp_rms_cuda(d_rhorms, d_partial_result, d_lnrho);
		cudaDeviceSynchronize();
		cudaMemcpy(&rhorms, (float*)d_rhorms, sizeof(float), cudaMemcpyDeviceToHost); 
		cudaDeviceSynchronize(); 
	}
        //if (iproc==0) 
        //printf(" step = %i; t = %e; dt = %e; umax = %e; umin = %e; urms = %e",step, t, dt, umax, umin, urms);
	//printf(" step = %i; t = %e; dt = %e; umax = %e; umin = %e; urms = %e; \n uxrms = %e; uyrms = %e; uzrms = %e; \n uxmax = %e; uymax = %e; uzmax = %e; \n uxmin = %e; uymin = %e; uzmin = %e; \n rhomax = %e; rhomin = %e; rhorms = %e \n", 
        //        step, t, dt, umax, umin, urms, uxrms, uyrms, uzrms, uxmax, uymax, uzmax, uxmin, uymin, uzmin, rhomax, rhomin, rhorms);
	
	//Save the step into a file 
	//save_ts(t, dt, step, urms, uxrms, uyrms, uzrms, uxmax, uymax, uzmax, rhorms, umax, rhomax, uxmin, uymin, uzmin, rhomin, umin);
}


