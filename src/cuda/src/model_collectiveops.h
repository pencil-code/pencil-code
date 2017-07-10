
#pragma once


float model_max_vec(float* uu_x, float* uu_y, float* uu_z);

float model_min_vec(float* uu_x, float* uu_y, float* uu_z);

float model_vec_rms(float* vec_x, float* vec_y, float* vec_z);

float model_max_scal(float* scal);

float model_min_scal(float* scal);

float model_scal_rms(float* scal);

float model_scal_exp_rms(float* scal);
