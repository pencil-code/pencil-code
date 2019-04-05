
#pragma once
#include "smem.cuh"

__device__ float interp_shear(float* array, int ix, int iy, int iz, int sid_depth, int sid_y,
                             float s_coord_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH], float s_val_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH]);

__device__ void epicyclic(float* u_epic_y, int sid_row, int sid_column, int sid_depth, float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH]);

__device__ void shear_flow(float* u_shear_y, float grid_idx_x);
