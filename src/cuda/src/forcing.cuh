
#pragma once

#include "smem.cuh"

__device__ void forcing_simple(float* fx, float* fy, float* fz, float grid_idx_x, float grid_idx_y, float grid_idx_z);

__device__ void forcing_nonhelical(float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                               	   float s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                              	   float s_uu_z[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
				   int sid_row, int sid_column, int sid_depth,    
				   float grid_idx_x, float grid_idx_y, float grid_idx_z);

int calc_k_vector_size();

void initialize_k_vectors(float* kk_x, float* kk_y, float* kk_z, float* kaver, int nk);

void choose_random_k(float* kk_vec_x, float* kk_vec_y, float* kk_vec_z, float* kk_x, float* kk_y, float* kk_z, int nk);

float choose_random_phi();

void get_forcing_unit_vector(float* ex, float* ey, float* ez);

void fkt_forcing_coefficient(float* forcing_kk_part_x, float* forcing_kk_part_y, float* forcing_kk_part_z, float kk_vec_x, float kk_vec_y, float kk_vec_z, float dt);
