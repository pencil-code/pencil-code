


#include <float.h>
#include <math.h>

#include "model_collectiveops.h"
#include "defines.h"
//#include "../defines.h"


float model_max_vec(float* uu_x, float* uu_y, float* uu_z) {

	float max_val = -FLT_MAX;
	for (int k=CZ_BOT; k < CZ_TOP; k++) {
		for (int j=CY_BOT; j < CY_TOP; j++) {
			for (int i=CX_BOT; i < CX_TOP; i++) {
				int idx = i + j*NX + k*NX*NY;
				float candidate = sqrtf(uu_x[idx]*uu_x[idx] +
							uu_y[idx]*uu_y[idx] +
							uu_z[idx]*uu_z[idx]);
				if(candidate > max_val) 
					max_val = candidate;
			}
		}
	}
	return max_val;
}

float model_min_vec(float* uu_x, float* uu_y, float* uu_z) {

	float min_val = FLT_MAX;
	for (int k=CZ_BOT; k < CZ_TOP; k++) {
		for (int j=CY_BOT; j < CY_TOP; j++) {
			for (int i=CX_BOT; i < CX_TOP; i++) {
				int idx = i + j*NX + k*NX*NY;
				float candidate = sqrtf(uu_x[idx]*uu_x[idx] +
							uu_y[idx]*uu_y[idx] +
							uu_z[idx]*uu_z[idx]);
				if(candidate < min_val) 
					min_val = candidate;
			}
		}
	}
	return min_val;
}


float model_vec_rms(float* uu_x, float* uu_y, float* uu_z) {
	return FLT_MAX;
}

float model_max_scal(float* scal) {

	float max_scal = -FLT_MAX;
	for (int k=CZ_BOT; k < CZ_TOP; k++) {
		for (int j=CY_BOT; j < CY_TOP; j++) {
			for (int i=CX_BOT; i < CX_TOP; i++) {
				int idx = i + j*NX + k*NX*NY;
				if(scal[idx] > max_scal) 
					max_scal = scal[idx];
			}
		}
	}
	return max_scal;
}

float model_min_scal(float* d_scal) {
	return FLT_MAX;
}

float model_scal_rms(float* scal) {
	return FLT_MAX;
}

float model_scal_exp_rms(float* scal) {
	return FLT_MAX;
}













