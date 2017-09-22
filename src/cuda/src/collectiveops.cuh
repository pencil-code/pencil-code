
#pragma once

void max_vec_cuda(float* d_vec_max, float* d_partial_result, float* d_vec_x, float* d_vec_y, float* d_vec_z);

void min_vec_cuda(float* d_vec_min, float* d_partial_result, float* d_vec_x, float* d_vec_y, float* d_vec_z);

void max_scal_cuda(float* d_scal_max, float* d_partial_result, float* d_scal);

void min_scal_cuda(float* d_scal_min, float* d_partial_result, float* d_scal);

void vec_rms_cuda(float* d_vec_rms, float* d_partial_result, float* d_vec_x, float* d_vec_y, float* d_vec_z, bool root=true);

void scal_rms_cuda(float* d_scal_rms, float* d_partial_result, float* d_scal, bool sqr=true, bool root=true);

void scal_exp_rms_cuda(float* d_scal_rms, float* d_partial_result, float* d_scal, bool sqr=true, bool root=true);
