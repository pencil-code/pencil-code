//
// Coriolis force is calculated in this module 
//

#include "coriolis.cuh"
#define EXTERN extern
#include "dconsts.cuh"

__device__ void coriolis(float* f_cor_x, float* f_cor_y, float u_shear_y, int sid_row, int sid_column, int sid_depth, 
                         float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                         float s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	/*
	Calculate the coriolis component into the velocity.  

	with \mathbf{\Omega} = \[0, 0, \Omega_0\]
	*/

	*f_cor_x = -2.0*(-d_OMEGA*(s_uu_y[sid_row][sid_column][sid_depth] + u_shear_y) );
	*f_cor_y = -2.0*( d_OMEGA*s_uu_x[sid_row][sid_column][sid_depth] );

	//printf("f_cor_x = %e, f_cor_y = %e", f_cor_x, f_cor_y);

}
