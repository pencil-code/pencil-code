#pragma once

float timestep_cuda(float* d_umax, float* d_partial_result, float* d_uu_x, float* d_uu_y, float* d_uu_z);

void timeseries_diagnostics_cuda(/*float* d_umax, float* d_umin, float* d_urms, float* d_uxrms, 
                                 float* d_uyrms, float* d_uzrms, float* d_rhorms,
                                 float* d_rhomax, float* d_uxmax, float* d_uymax, float* d_uzmax,
                                 float* d_rhomin, float* d_uxmin, float* d_uymin, float* d_uzmin,*/
                                 int step, float dt, double t);
void max_diffus();
void max_advec();

