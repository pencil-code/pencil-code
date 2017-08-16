
#pragma once

#include "smem.cuh"

//JP: This was originally an extern __device__ float function, 
//JP: removed extern because it shouldn't be needed.
__device__ float continuity(int sid_row, int sid_column, int sid_depth,
				   float u_shear_y,
                            	   float s_lnrho[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
				   float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                            	   float s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
				   float s_uu_z[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH]);


