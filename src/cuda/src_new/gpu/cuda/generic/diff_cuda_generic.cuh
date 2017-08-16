/*
*   Differential equations used for integration in rk3_cuda_generic.
*
*   The function definitions must be done in the header file, such that the
*   implementation is visible for the caller during compile time.
*   This way the compiler can optimize these functions instead of replacing
*   them with expensive ABI calls.
*/
#pragma once
#include "gpu/cuda/core/dconsts_core.cuh"
#include "common/datatypes.h"


/*
*   Defines specific to this smem layout
*/
#if USE_DOUBLE_PRECISION == 0
    #define RK_ELEMS_PER_THREAD 8
    #define RK_THREADS_X 8
    #define RK_THREADS_Y 8
    #define RK_THREADS_Z 8
#else
    #define RK_ELEMS_PER_THREAD 32 
    #define RK_THREADS_X 4 
    #define RK_THREADS_Y 4
    #define RK_THREADS_Z 4
#endif
#define RK_THREADBLOCK_SIZE (RK_THREADS_X * RK_THREADS_Y * RK_THREADS_Z)

//Shared memory sizes for results (d_lnrho etc, NOTE: smem uses row-major, so fastest varying dim
//should be columns to avoid bank conflicts)
#define SHARED_SIZE_ROW (RK_THREADS_Y + 2*BOUND_SIZE)//Rows
#define SHARED_SIZE_COL (RK_THREADS_X + 2*BOUND_SIZE)//Columns
#define SHARED_SIZE_DEPTH (RK_THREADS_Z + 2*BOUND_SIZE)

//Helper definition which is used to wrap the depth values
//around the shared memory block
#define WRAP_SMEM_DEPTH(x) ((x + 2*SHARED_SIZE_DEPTH) % SHARED_SIZE_DEPTH)

//Shared memory sizes for intermediate results (d_w_lnrho etc)
#define SHARED_SIZE_W_ROW RK_THREADS_Y //Rows
#define SHARED_SIZE_W_COL RK_THREADS_X //Columns
#define SHARED_SIZE_W_DEPTH RK_THREADS_Z//Depth


//
// Derivative operators, 1st order 
//

//sid_row represents the threadIdx.y and y axis in grid coords
//sid_column represents threadIdx.x and x axis in grid coords 
//sid_depth represents threadIdx.z and z asxis in grid coords
//
//and shared memory is arranged like: 
//					s_scal[sid_row][sid_column][sid_depth] or
// 					s_scal[y-axis in grid][x-axis in grid][z-axis in grid]
static __device__ real der_scaly(	int sid_row, int sid_column, int sid_depth, 
                            real s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
   	//
   	// Single derivative in y-direction
	//
	real res ;

	res = (
	-            s_scal[sid_row-3][sid_column][sid_depth] 
	+ (real) 9.0  * s_scal[sid_row-2][sid_column][sid_depth] 
	- (real) 45.0 * s_scal[sid_row-1][sid_column][sid_depth] 
	+ (real) 45.0 * s_scal[sid_row+1][sid_column][sid_depth] 
	- (real) 9.0  * s_scal[sid_row+2][sid_column][sid_depth] 
	+            s_scal[sid_row+3][sid_column][sid_depth] )
	* d_DIFF1_DY_DIV;
	// / ( (real) 60.0*d_DY ); //MV: Made these divisions to go away. -> need only be calculated once and used as a constant. 

   return res;
}


static __device__ real der_scalx(	int sid_row, int sid_column, int sid_depth, 
                            real s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Single derivative in x-direction
	//
	real res ;

	res = (
	-            s_scal[sid_row][sid_column-3][sid_depth] 
	+ (real) 9.0  * s_scal[sid_row][sid_column-2][sid_depth] 
	- (real) 45.0 * s_scal[sid_row][sid_column-1][sid_depth] 
	+ (real) 45.0 * s_scal[sid_row][sid_column+1][sid_depth] 
	- (real) 9.0  * s_scal[sid_row][sid_column+2][sid_depth] 
	+            s_scal[sid_row][sid_column+3][sid_depth] )
	* d_DIFF1_DX_DIV;
	// / ( (real) 60.0*d_DX ); 

	return res;
}


//Bank conflicts: ~2ms slower than other der_scals
static __device__ real der_scalz(	int sid_row, int sid_column, int sid_depth, 
                            real s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Single derivative in z-direction
	//
	real res ;

	res = (
	-            s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth-3)] 
	+ (real) 9.0  * s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth-2)] 
	- (real) 45.0 * s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth-1)] 
	+ (real) 45.0 * s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth+1)] 
	- (real) 9.0  * s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth+2)] 
	+            s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth+3)] )
	* d_DIFF1_DZ_DIV;
	// / ( (real) 60.0*d_DZ );

	return res;
}


//
// Derivative operators, 2nd order 
//
static __device__ real der2_scaly(	int sid_row, int sid_column, int sid_depth, 
                            real s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Double derivative in y-direction
	//
	real res;

	res = (
	  (real) 2.0   * s_scal[sid_row-3][sid_column][sid_depth] 
	- (real) 27.0  * s_scal[sid_row-2][sid_column][sid_depth] 
	+ (real) 270.0 * s_scal[sid_row-1][sid_column][sid_depth] 
	- (real) 490.0 * s_scal[sid_row  ][sid_column][sid_depth] 
	+ (real) 270.0 * s_scal[sid_row+1][sid_column][sid_depth] 
	- (real) 27.0  * s_scal[sid_row+2][sid_column][sid_depth] 
	+ (real) 2.0   * s_scal[sid_row+3][sid_column][sid_depth] )
	* d_DIFF2_DY_DIV;
	// / ( d_FLT_180*d_DY*d_DY );

	return res;
}


static __device__ real der2_scalx(	int sid_row, int sid_column, int sid_depth, 
                            real s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Double derivative in x-direction
	//
	real res;

	res = (
	  (real) 2.0   * s_scal[sid_row][sid_column-3][sid_depth]
	- (real) 27.0  * s_scal[sid_row][sid_column-2][sid_depth] 
	+ (real) 270.0 * s_scal[sid_row][sid_column-1][sid_depth] 
	- (real) 490.0 * s_scal[sid_row][sid_column  ][sid_depth]
	+ (real) 270.0 * s_scal[sid_row][sid_column+1][sid_depth]
	- (real) 27.0  * s_scal[sid_row][sid_column+2][sid_depth]
	+ (real) 2.0   * s_scal[sid_row][sid_column+3][sid_depth] )
	* d_DIFF2_DX_DIV;
	// / ( d_FLT_180*d_DX*d_DX );

	return res;
}


//Bank conflicts: ~0.5ms slower than other der2_scals
static __device__ real der2_scalz(	int sid_row, int sid_column, int sid_depth, 
                            real s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Double derivative in z-direction
	//
	real res;

	res = (
	  (real) 2.0   * s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth-3)] 
	- (real) 27.0  * s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth-2)] 
	+ (real) 270.0 * s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth-1)] 
	- (real) 490.0 * s_scal[sid_row][sid_column][sid_depth] 
	+ (real) 270.0 * s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth+1)] 
	- (real) 27.0  * s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth+2)] 
	+ (real) 2.0   * s_scal[sid_row][sid_column][WRAP_SMEM_DEPTH(sid_depth+3)] )
	* d_DIFF2_DZ_DIV;
	// / ( d_FLT_180*d_DZ*d_DZ );

	return res;
}


static __device__ real dernm_scal(	int sid_row, int sid_column, int sid_depth, int m, int n, 
                            real s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// (d/dm)(d scal/dn)
	// Adapted bidiagonal mixed-derivative operator from the Pencil code
	// x == 1  y == 2  z == 3 Numeric labels of the axis
	//
	real res;


	if ((m == 1 && n == 2) || (m == 2 && n == 1)) {  
		res = d_DIFFMN_DXDY_DIV*( 
		  (real) 270.0*( s_scal[sid_row+1][sid_column+1][sid_depth] 
	  	            - s_scal[sid_row-1][sid_column+1][sid_depth]
		            + s_scal[sid_row-1][sid_column-1][sid_depth]
		            - s_scal[sid_row+1][sid_column-1][sid_depth] )
		  -(real) 27.0*( s_scal[sid_row+2][sid_column+2][sid_depth] 
		            - s_scal[sid_row-2][sid_column+2][sid_depth] 
		            + s_scal[sid_row-2][sid_column-2][sid_depth] 
		            - s_scal[sid_row+2][sid_column-2][sid_depth] )
		  + (real) 2.0*( s_scal[sid_row+3][sid_column+3][sid_depth] 
		            - s_scal[sid_row-3][sid_column+3][sid_depth]
		            + s_scal[sid_row-3][sid_column-3][sid_depth]
		            - s_scal[sid_row+3][sid_column-3][sid_depth] )  
		          );
	}             
	else if ((m == 2 && n == 3) || (m == 3 && n == 2)) {
		res = d_DIFFMN_DYDZ_DIV*( 
		  (real) 270.0*( s_scal[sid_row+1][sid_column][WRAP_SMEM_DEPTH(sid_depth+1)] 
	 	            - s_scal[sid_row+1][sid_column][WRAP_SMEM_DEPTH(sid_depth-1)]
		            + s_scal[sid_row-1][sid_column][WRAP_SMEM_DEPTH(sid_depth-1)]
		            - s_scal[sid_row-1][sid_column][WRAP_SMEM_DEPTH(sid_depth+1)] )
		  -(real) 27.0*( s_scal[sid_row+2][sid_column][WRAP_SMEM_DEPTH(sid_depth+2)] 
		            - s_scal[sid_row+2][sid_column][WRAP_SMEM_DEPTH(sid_depth-2)] 
		            + s_scal[sid_row-2][sid_column][WRAP_SMEM_DEPTH(sid_depth-2)] 
		            - s_scal[sid_row-2][sid_column][WRAP_SMEM_DEPTH(sid_depth+2)] )
		  + (real) 2.0*( s_scal[sid_row+3][sid_column][WRAP_SMEM_DEPTH(sid_depth+3)] 
		            - s_scal[sid_row+3][sid_column][WRAP_SMEM_DEPTH(sid_depth-3)] 
		            + s_scal[sid_row-3][sid_column][WRAP_SMEM_DEPTH(sid_depth-3)] 
		            - s_scal[sid_row-3][sid_column][WRAP_SMEM_DEPTH(sid_depth+3)] )  
		          );
	} 
	else if ((m == 3 && n == 1) || (m == 1 && n == 3)) {
		res = d_DIFFMN_DXDZ_DIV*( 
		  (real) 270.0*( s_scal[sid_row][sid_column+1][WRAP_SMEM_DEPTH(sid_depth+1)] 
	 	            - s_scal[sid_row][sid_column-1][WRAP_SMEM_DEPTH(sid_depth+1)]
		            + s_scal[sid_row][sid_column-1][WRAP_SMEM_DEPTH(sid_depth-1)]
		            - s_scal[sid_row][sid_column+1][WRAP_SMEM_DEPTH(sid_depth-1)] )
		  -(real) 27.0*( s_scal[sid_row][sid_column+2][WRAP_SMEM_DEPTH(sid_depth+2)] 
		            - s_scal[sid_row][sid_column-2][WRAP_SMEM_DEPTH(sid_depth+2)] 
		            + s_scal[sid_row][sid_column-2][WRAP_SMEM_DEPTH(sid_depth-2)] 
		            - s_scal[sid_row][sid_column+2][WRAP_SMEM_DEPTH(sid_depth-2)] )
		  + (real) 2.0*( s_scal[sid_row][sid_column+3][WRAP_SMEM_DEPTH(sid_depth+3)] 
		            - s_scal[sid_row][sid_column-3][WRAP_SMEM_DEPTH(sid_depth+3)] 
		            + s_scal[sid_row][sid_column-3][WRAP_SMEM_DEPTH(sid_depth-3)] 
		            - s_scal[sid_row][sid_column+3][WRAP_SMEM_DEPTH(sid_depth-3)] )  
		          );

	} else {
        res = NAN;
    }
                        
	return res;

}


//
//   Partial derivative oparations, simple
//
static __device__ real div_vec(	int sid_row, int sid_column, int sid_depth, 
                            real s_vecx[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                            real s_vecy[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                            real s_vecz[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Divergence of a vector (\nabla \cdot \vec{u})
	// div vec
	// 
	real res;
	
	res = 
	  der_scalx(sid_row, sid_column, sid_depth, s_vecx)
	+ der_scaly(sid_row, sid_column, sid_depth, s_vecy)
	+ der_scalz(sid_row, sid_column, sid_depth, s_vecz);

	return res;
}


static __device__ void grad(	real* gx, real* gy, real* gz, 
                        int sid_row, int sid_column, int sid_depth, 
                        real s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Gradient of a scalar (\nabla u)
	// Nabla scal
	//
	*gx = der_scalx(sid_row, sid_column, sid_depth, s_scal);

	*gy = der_scaly(sid_row, sid_column, sid_depth, s_scal);
	
	*gz = der_scalz(sid_row, sid_column, sid_depth, s_scal);
}


//
//   Partial derivative oparations, complex
//
static __device__ void laplace_vec(real* lux, real* luy, real* luz, 
                            int sid_row, int sid_column, int sid_depth, 
                            real s_vecx[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                            real s_vecy[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                            real s_vecz[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Laplacian of a vector (\nabla^2 \vec{u})
	// Nabla^2 vec
	//
	*lux = der2_scalx(sid_row, sid_column, sid_depth, s_vecx)
	     + der2_scaly(sid_row, sid_column, sid_depth, s_vecx)
	     + der2_scalz(sid_row, sid_column, sid_depth, s_vecx);
	
	*luy = der2_scalx(sid_row, sid_column, sid_depth, s_vecy)
	     + der2_scaly(sid_row, sid_column, sid_depth, s_vecy)
	     + der2_scalz(sid_row, sid_column, sid_depth, s_vecy);
	
	*luz = der2_scalx(sid_row, sid_column, sid_depth, s_vecz)
	     + der2_scaly(sid_row, sid_column, sid_depth, s_vecz)
	     + der2_scalz(sid_row, sid_column, sid_depth, s_vecz);

}


static __device__ real vec_dot_nabla_scal(	int sid_row, int sid_column, int sid_depth, 
                                    real s_vecx[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                                    real s_vecy[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                                    real s_vecz[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                                    real s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	// Dot product of a scalar gradient (\vec{u} \cdot \nabla u)
	// vec dot nabla scal
	//
	real res;

	res = 
	  s_vecx[sid_row][sid_column][sid_depth] * der_scalx(sid_row, sid_column, sid_depth, s_scal)
	+ s_vecy[sid_row][sid_column][sid_depth] * der_scaly(sid_row, sid_column, sid_depth, s_scal)
	+ s_vecz[sid_row][sid_column][sid_depth] * der_scalz(sid_row, sid_column, sid_depth, s_scal);
	 
	return res;
	
}


static __device__ real vec_dot_nabla_scal_shear_y( int sid_row, int sid_column, int sid_depth, real u_shear_y,
                                            real s_vecx[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                                            real s_vecy[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                                            real s_vecz[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                                            real s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{	
	//
	// Dot product of a scalar gradient (\vec{u} \cdot \nabla u)
	// With shearing included in velocity components
	// vec dot nabla scal
	//
	real res;

	res = 
	  s_vecx[sid_row][sid_column][sid_depth]               * der_scalx(sid_row, sid_column, sid_depth, s_scal)
	+ (s_vecy[sid_row][sid_column][sid_depth] + u_shear_y) * der_scaly(sid_row, sid_column, sid_depth, s_scal)
	+ s_vecz[sid_row][sid_column][sid_depth]               * der_scalz(sid_row, sid_column, sid_depth, s_scal);
	 
	return res;
	
}


static __device__ void grad_div_vec(	real* gdvx, real* gdvy, real* gdvz, 
                                int sid_row, int sid_column, int sid_depth,
                                real s_vecx[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                                real s_vecy[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                                real s_vecz[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	//
	//The Gradient of a divergence ( \nabla (\nabla \cdot \vec{u}) )
	// Nabla (Nabla dot Vec)
	//
	*gdvx = der2_scalx(sid_row, sid_column, sid_depth, s_vecx) 
	      + dernm_scal(sid_row, sid_column, sid_depth, 1, 2, s_vecy) 
	      + dernm_scal(sid_row, sid_column, sid_depth, 1, 3, s_vecz);

	*gdvy = dernm_scal(sid_row, sid_column, sid_depth, 2, 1, s_vecx)  
	      + der2_scaly(sid_row, sid_column, sid_depth, s_vecy) 
	      + dernm_scal(sid_row, sid_column, sid_depth, 2, 3, s_vecz);

	*gdvz = dernm_scal(sid_row, sid_column, sid_depth, 3, 1, s_vecx) 
	      + dernm_scal(sid_row, sid_column, sid_depth, 3, 2, s_vecy) 
	      + der2_scalz(sid_row, sid_column, sid_depth, s_vecz);

}










































