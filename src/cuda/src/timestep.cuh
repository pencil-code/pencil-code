#pragma once

float timestep_cuda(float* d_umax, float* d_partial_result, float* d_uu_x, float* d_uu_y, float* d_uu_z);

void timeseries_diagnostics_cuda(/*float* d_umax, float* d_umin, float* d_urms, float* d_uxrms, 
                                 float* d_uyrms, float* d_uzrms, float* d_rhorms,*/
                                 int step, float dt, double t);
void get_maxscal_from_device(float & maxscal,float *d_src);
void get_minscal_from_device(float & minscal,float *d_src);
