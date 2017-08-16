

/*
* Contains the Navier-Stokes equation, linear viscosity and related subroutines
* (currently dummy functions)
*/

#include "hydro.cuh"
#include "diff.cuh"
#define EXTERN extern
#include "dconsts.cuh"
#include "coriolis.cuh"
#include "shear.cuh"

__device__ void S_grad_lnrho(	float* sgrhox, float* sgrhoy, float* sgrhoz, 
				int sid_row, int sid_column, int sid_depth, 
				float s_lnrho[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH],
                             	float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
				float s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                             	float s_uu_z[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Traceless rate of strain tensor
	//

	*sgrhox = *sgrhoy = *sgrhoz = NAN;
}

__device__ void nu_const(	float* fviscx, float* fviscy, float* fviscz, 
				int sid_row, int sid_column, int sid_depth, 
				float s_lnrho[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH],
                         	float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
				float s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                         	float s_uu_z[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Calculate the viscosity parameter 
	// fvisc = nu (nabla^2 u + 1/3 nabla nabla dot u + 2 S dot nabla lnrho) 
	//
	*fviscx = *fviscy = *fviscz = NAN;
}


__device__ void navier_stokes(	float* momx, float* momy, float* momz, 
				int sid_row, int sid_column, int sid_depth,    
				float u_shear_y,
				float s_lnrho[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH],
                              	float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                              	float s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                              	float s_uu_z[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Navier-Stokes equation
	// du / dt = - u dot nabla u - cs^2 nabla lnrho + f_visc + f_force
	// 
	*momx = *momy = *momz = NAN;
}


