
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
}
void get_minvec_from_device(float &minvec,float *d_src_x,float *d_src_y,float *d_src_z)
{
        min_vec_cuda(d_scaldiag, d_partial_result, d_src_x, d_src_y, d_src_z);
        cudaDeviceSynchronize();
        cudaMemcpy(&minvec, d_scaldiag, sizeof(float), cudaMemcpyDeviceToHost);
}
void get_maxvec_from_device(float &maxvec,float *d_src_x,float *d_src_y,float *d_src_z)
{
        max_vec_cuda(d_scaldiag, d_partial_result, d_src_x, d_src_y, d_src_z);
        cudaDeviceSynchronize();
        cudaMemcpy(&maxvec, d_scaldiag, sizeof(float), cudaMemcpyDeviceToHost);
}
void get_rmsvec_from_device(float &rmsvec,float *d_src_x,float *d_src_y,float *d_src_z)
{
        vec_rms_cuda(d_scaldiag, d_partial_result, d_src_x, d_src_y, d_src_z);
        cudaDeviceSynchronize();
        cudaMemcpy(&rmsvec, d_scaldiag, sizeof(float), cudaMemcpyDeviceToHost);
}
void get_sumsqscal_from_device(float & sumsquare,float *d_src,bool exp=false)
{ 	
	if (exp) 
		scal_exp_rms_cuda(d_scaldiag, d_partial_result, d_src, true, false);
	else
		scal_rms_cuda(d_scaldiag, d_partial_result, d_src, true, false);

	cudaDeviceSynchronize();
        cudaMemcpy(&sumsquare, d_scaldiag, sizeof(float), cudaMemcpyDeviceToHost);
}
void get_sumscal_from_device(float & sum,float *d_src,bool exp=false)
{
	if (exp) 
		scal_exp_rms_cuda(d_scaldiag, d_partial_result, d_src, false, false);
	else
		scal_rms_cuda(d_scaldiag, d_partial_result, d_src, false, false);

	cudaDeviceSynchronize();
        cudaMemcpy(&sum, d_scaldiag, sizeof(float), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

void timeseries_diagnostics_cuda(int step, float dt, double t)
{
	//Calculate and save all of the diagnostic variables calculated within the CUDA devices. 

	static float diag; 
	
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
                diag=-diag;
		save_name(diag,idiag_uxmin);
        }
        if (idiag_uymin>0) {
        	get_minscal_from_device(diag,d_uu_y);
                diag=-diag;
		save_name(diag,idiag_uymin);
        }
        if (idiag_uzmin>0) {
        	get_minscal_from_device(diag,d_uu_z);
                diag=-diag;
		save_name(diag,idiag_uzmin);
        }
        if (idiag_umax>0) {
		get_maxvec_from_device(diag, d_uu_x, d_uu_y, d_uu_z);
                save_name(diag,idiag_umax);
        }
        if (idiag_rhomax>0) {
        	get_maxscal_from_device(diag,d_lnrho);
		if (!ldensity_nolog) diag = exp(diag);       //Change away from the logarithmic form
		save_name(diag,idiag_rhomax);
        }
        if (idiag_rhomin>0) {
        	get_minscal_from_device(diag,d_lnrho);
		if (!ldensity_nolog) diag = exp(diag);       //Change away from the logarithmic form
                diag=-diag;
		save_name(diag,idiag_rhomin);
        }
        if (idiag_umin>0) {
		get_minvec_from_device(diag, d_uu_x, d_uu_y, d_uu_z);
                diag=-diag;
                save_name(diag,idiag_umin);
	}
 	if (idiag_urms){
		get_rmsvec_from_device(diag, d_uu_x, d_uu_y, d_uu_z);
                save_name(diag,idiag_urms);
	}
 	if (idiag_uxrms){
		get_sumsqscal_from_device(diag,d_uu_x);
                save_name(diag,idiag_uxrms);
	}
 	if (idiag_uyrms){
		get_sumsqscal_from_device(diag,d_uu_y);
                save_name(diag,idiag_uyrms);
	}
 	if (idiag_uzrms){
		get_sumsqscal_from_device(diag,d_uu_z);
                save_name(diag,idiag_uzrms);
	}
 	if (idiag_rhom){
		get_sumscal_from_device(diag,d_lnrho,!ldensity_nolog);
                save_name(diag,idiag_rhom);
	}
 	if (idiag_rhorms){
		get_sumsqscal_from_device(diag,d_lnrho,!ldensity_nolog);
                save_name(diag,idiag_rhorms);
	}
        //if (iproc==0) 
        //printf(" step = %i; t = %e; dt = %e; umax = %e; umin = %e; urms = %e",step, t, dt, umax, umin, urms);
	//printf(" step = %i; t = %e; dt = %e; umax = %e; umin = %e; urms = %e; \n uxrms = %e; uyrms = %e; uzrms = %e; \n uxmax = %e; uymax = %e; uzmax = %e; \n uxmin = %e; uymin = %e; uzmin = %e; \n rhomax = %e; rhomin = %e; rhorms = %e \n", 
        //        step, t, dt, umax, umin, urms, uxrms, uyrms, uzrms, uxmax, uymax, uzmax, uxmin, uymin, uzmin, rhomax, rhomin, rhorms);
	
	//Save the step into a file 
	//save_ts(t, dt, step, urms, uxrms, uyrms, uzrms, uxmax, uymax, uzmax, rhorms, umax, rhomax, uxmin, uymin, uzmin, rhomin, umin);
}


