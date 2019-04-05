

/*
* UNUSED, located atm in integrators.cu
*/



/*
#include "diff.cuh"
#define EXTERN extern
#include "dconsts.cuh"

#include <cstdio>

//sid_row represents the threadIdx.y and y axis in grid coords
//sid_column represents threadIdx.x and x axis in grid coords 
//
//and shared memory is arranged like: 
//					s_scal[sid_row][sid_column] or
// 					s_scal[y-axis in grid][x-axis in grid]

//------------------------------------------------------------------------------------------------------
//
// Derivative operators, 1st order 
//
__device__ float der_scalx(	int sid_row, int sid_column, 
				float s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL])
{
	//
	// Single derivative in x-direction
	//

	float res ;

	res = (
	-            s_scal[sid_row][sid_column-3] 
	+ d_FLT_9  * s_scal[sid_row][sid_column-2] 
	- d_FLT_45 * s_scal[sid_row][sid_column-1] 
	+ d_FLT_45 * s_scal[sid_row][sid_column+1] 
	- d_FLT_9  * s_scal[sid_row][sid_column+2] 
	+            s_scal[sid_row][sid_column+3] )
	* d_DIFF1_DX_DIV;
	// / ( d_FLT_60*d_DX ); 

	return res;
}


__device__ float der_scaly(	int sid_row, int sid_column, 
				float s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL])
{
   	//
   	// Single derivative in y-direction
	//

	float res ;

	res = (
	-            s_scal[sid_row-3][sid_column] 
	+ d_FLT_9  * s_scal[sid_row-2][sid_column] 
	- d_FLT_45 * s_scal[sid_row-1][sid_column] 
	+ d_FLT_45 * s_scal[sid_row+1][sid_column] 
	- d_FLT_9  * s_scal[sid_row+2][sid_column] 
	+            s_scal[sid_row+3][sid_column] )
	* d_DIFF1_DY_DIV;
	// / ( d_FLT_60*d_DY ); //MV: Made these divisions to go away. -> need only be calculated once and used as a constant. 

   return res;
}

__device__ float der_scalz(	float behind3, float behind2, float behind1,
				float infront1, float infront2, float infront3)
{
	//
	// Single derivative in z-direction
	//

	float res ;

	res = (
	-            behind3 
	+ d_FLT_9  * behind2
	- d_FLT_45 * behind1 
	+ d_FLT_45 * infront1 
	- d_FLT_9  * infront2 
	+            infront3 )
	* d_DIFF1_DZ_DIV;
	// / ( d_FLT_60*d_DZ );

	return res;
}
//------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------
//
// Derivative operators, 2nd order 
//
__device__ float der2_scalx(int sid_row, int sid_column, float s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL])
{
	//
	// Double derivative in x-direction
	//
	float res;

	res = (
	  d_FLT_2   * s_scal[sid_row][sid_column-3]
	- d_FLT_27  * s_scal[sid_row][sid_column-2] 
	+ d_FLT_270 * s_scal[sid_row][sid_column-1] 
	- d_FLT_490 * s_scal[sid_row][sid_column  ]
	+ d_FLT_270 * s_scal[sid_row][sid_column+1]
	- d_FLT_27  * s_scal[sid_row][sid_column+2]
	+ d_FLT_2   * s_scal[sid_row][sid_column+3] )
	* d_DIFF2_DX_DIV;
	// / ( d_FLT_180*d_DX*d_DX );

	return res;

}

__device__ float der2_scaly(int sid_row, int sid_column, float s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL])
{
	//
	// Double derivative in y-direction
	//
	float res;

	res = (
	  d_FLT_2   * s_scal[sid_row-3][sid_column] 
	- d_FLT_27  * s_scal[sid_row-2][sid_column] 
	+ d_FLT_270 * s_scal[sid_row-1][sid_column] 
	- d_FLT_490 * s_scal[sid_row  ][sid_column] 
	+ d_FLT_270 * s_scal[sid_row+1][sid_column] 
	- d_FLT_27  * s_scal[sid_row+2][sid_column] 
	+ d_FLT_2   * s_scal[sid_row+3][sid_column] )
	* d_DIFF2_DY_DIV;
	// / ( d_FLT_180*d_DY*d_DY );

	return res;
}

__device__ float der2_scalz(int sid_row, int sid_column, float s_scal[SHARED_SIZE_ROW][SHARED_SIZE_COL],
				float behind3, float behind2, float behind1,
				float infront1, float infront2, float infront3)
{
	//
	// Double derivative in z-direction
	//
	float res;

	res = (
	  d_FLT_2   * behind3 
	- d_FLT_27  * behind2 
	+ d_FLT_270 * behind1 
	- d_FLT_490 * s_scal[sid_row][sid_column] 
	+ d_FLT_270 * infront1 
	- d_FLT_27  * infront2 
	+ d_FLT_2   * infront3 )
	* d_DIFF2_DY_DIV;
	// / ( d_FLT_180*d_DY*d_DY );

	return res;

}
//------------------------------------------------------------------------------------------------------
*/

