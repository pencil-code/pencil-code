
#include "continuity.cuh"
#include "diff.cuh"
#include "shear.cuh"

#include "diagnostics.cuh"

#include <cstdio>

/*
* Contains the continuity equation and related subroutines
* (Currently a dummy function)
*/
__device__ float continuity(	int sid_row, int sid_col, int sid_depth,
				float u_shear_y,
                            	float s_lnrho[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
				float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                            	float s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
				float s_uu_z[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// calculate the time derivative of the continuity equation
	// dlnrho / dt = - u dot nabla lnrho - nabla dot u
	//

	return NAN;
}

