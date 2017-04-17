
#include <stdio.h>
#include <math.h>

#include "timestep.cuh"
#include "collectiveops.cuh"
#include "defines.h"
#include "io.h"

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

void timeseries_diagnostics_cuda(float* d_umax, float* d_umin, float* d_urms, float* d_uxrms, float* d_uyrms, float* d_uzrms, float* d_rhorms, 
                                 float* d_rhomax, float* d_uxmax, float* d_uymax, float* d_uzmax, 
                                 float* d_rhomin, float* d_uxmin, float* d_uymin, float* d_uzmin, 
                                 int step, float dt, float t, 
                                 float* d_partial_result, float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z)
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


	//Get uu_x max from d_uxmax 
	max_scal_cuda(d_uxmax, d_partial_result, d_uu_x);
	cudaDeviceSynchronize();
	cudaMemcpy(&uxmax, (float*)d_uxmax, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();

	//Get uu_y max from d_uymax
	max_scal_cuda(d_uymax, d_partial_result, d_uu_y);
	cudaDeviceSynchronize();
	cudaMemcpy(&uymax, (float*)d_uymax, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();

	//Get uu_z max from d_uzmax 
	max_scal_cuda(d_uzmax, d_partial_result, d_uu_z);
	cudaDeviceSynchronize();
	cudaMemcpy(&uzmax, (float*)d_uzmax, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();


	//Get rho max from d_lnrho
	max_scal_cuda(d_rhomax, d_partial_result, d_lnrho);
	cudaDeviceSynchronize();
	cudaMemcpy(&rhomax, (float*)d_rhomax, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();
	rhomax = exp(rhomax); //Change away from the logarithmic form


	//Get uu min from d_umin 
	min_vec_cuda(d_umin, d_partial_result, d_uu_x, d_uu_y, d_uu_z);
	cudaDeviceSynchronize();
	cudaMemcpy(&umin, (float*)d_umin, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();

	//Get uu_x min from d_uxmin 
	min_scal_cuda(d_uxmin, d_partial_result, d_uu_x);
	cudaDeviceSynchronize();
	cudaMemcpy(&uxmin, (float*)d_uxmin, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();

	//Get uu_y min from d_uymin 
	min_scal_cuda(d_uymin, d_partial_result, d_uu_y);
	cudaDeviceSynchronize();
	cudaMemcpy(&uymin, (float*)d_uymin, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();

	//Get uu_z min from d_uzmin 
	min_scal_cuda(d_uzmin, d_partial_result, d_uu_z);
	cudaDeviceSynchronize();
	cudaMemcpy(&uzmin, (float*)d_uzmin, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();


	//Get lnrho min from d_lnrho 
	max_scal_cuda(d_rhomin, d_partial_result, d_lnrho);
	cudaDeviceSynchronize();
	cudaMemcpy(&rhomin, (float*)d_rhomin, sizeof(float), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();
	rhomin = exp(rhomin); //Change away from the logarithmic form


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

	printf(" step = %i; t = %e; dt = %e; umax = %e; umin = %e; urms = %e; \n uxrms = %e; uyrms = %e; uzrms = %e; \n uxmax = %e; uymax = %e; uzmax = %e; \n uxmin = %e; uymin = %e; uzmin = %e; \n rhomax = %e; rhomin = %e; rhorms = %e \n", 
                step, t, dt, umax, umin, urms, uxrms, uyrms, uzrms, uxmax, uymax, uzmax, uxmin, uymin, uzmin, rhomax, rhomin, rhorms);
	
	//Save the step into a file 
	save_ts(t, dt, step, urms, uxrms, uyrms, uzrms, uxmax, uymax, uzmax, rhorms, umax, rhomax, uxmin, uymin, uzmin, rhomin, umin);

}


