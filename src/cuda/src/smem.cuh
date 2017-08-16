#pragma once
#include "defines_dims_PC.h"

//-----------------------------------------------
//-----Defines for rungekutta2N_cuda-------------
//-----------------------------------------------
#define RK_ELEMS_PER_THREAD_FIRST 128 //Number of elements evaluated by one thread in _first_half
#define RK_ELEMS_PER_THREAD_SECOND 1 //Number of elements evaluated by one thread in _second_half

//Dimensions of the thread block used in rungekutta2N_cuda
#define RK_THREADS_X 32 //Returns incorrect results with higher than 32, bugged somewhere
#define RK_THREADS_Y 4 //Returns incorrect results with lower than 4, bugged somewhere

#define RK_THREADS_PER_BLOCK (RK_THREADS_X*RK_THREADS_Y)

//Shared memory sizes for results (d_lnrho etc, NOTE: smem uses row-major, so fastest varying dim
//should be columns to avoid bank conflicts)
#define ytraversing 1

#define SHARED_SIZE_ROW (ytraversing*RK_THREADS_Y + 2*BOUND_SIZE)

#define SHARED_SIZE_COL (RK_THREADS_X + 2*BOUND_SIZE)//Columns

#define RK_THREADS_Z 1
#define SHARED_SIZE_DEPTH 1
//-----------------------------------


//-----------------------------------------------
//-----Defines for collective.cu-------------
//-----------------------------------------------
#define COL_THREADS_X 32
#define COL_THREADS_Y 8

#define COL_ELEMS_PER_THREAD 8
//-----------------------------------------------

//-----------------------------------------------
//-----Defines for shear.cu-------------
//-----------------------------------------------

//Assumed thus because in the boundcond.cu threadsPerBlock.x = 6; threadsPerBlock.y = 4; threadsPerBlock.z = 1;
//How many points are used for intepolation 
#define INTERP_ORDER 6
//How many y-points are interpolated in one block, as of threadsPerBlock.y = 4; 
#define INTERP_NPOINTS 4
//How many boundary layes (TOP and BOTTOM) are handled in one block, as of threadsPerBlock.x = 6; 
#define INTERP_DEPTH 6
 


