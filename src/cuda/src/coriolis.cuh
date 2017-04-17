
#pragma once
#include "smem.cuh"

__device__ void coriolis(float* f_cor_x, float* f_cor_y, float u_shear_y, int sid_row, int sid_column, int sid_depth, 
                         float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                         float s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH]);
