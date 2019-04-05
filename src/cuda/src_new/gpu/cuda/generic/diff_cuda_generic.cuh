/*
*   Difference formulae used for integration in rk3_cuda_generic.
*
*   The function definitions must be done in the header file, such that the
*   implementation is visible for the caller during compile time.
*   This way the compiler can optimize these functions instead of replacing
*   them with expensive ABI calls.
*/
#pragma once
#include "gpu/cuda/core/dconsts_core.cuh"
#include "common/datatypes.h"
#include "common/defines.h"

#define RK_THREADS_X (128)
#define RK_THREADS_Y (4)
#define RK_THREADS_Z (1)

#define RK_THREADS_PER_BLOCK (RK_THREADS_X*RK_THREADS_Y*RK_THREADS_Z)

#define RK_ELEMS_PER_THREAD (32)

//Shared memory
#define UNROLL_FACTOR (1)
#define SMEM_WIDTH (RK_THREADS_X + 2*BOUND_SIZE)
#define SMEM_HEIGHT (RK_THREADS_Y + 2*BOUND_SIZE)
#define SMEM_DEPTH (UNROLL_FACTOR)

//Support for variable boundary sizes in case we want to add padding later
#define XBOUND_SIZE BOUND_SIZE
#define YBOUND_SIZE BOUND_SIZE
#define ZBOUND_SIZE BOUND_SIZE


static __device__ __forceinline__ real der_scalx(const int smem_idx, const real* __restrict__ s_scal)
{
    const real res = ((s_scal[smem_idx + 3] - s_scal[smem_idx - 3]) 
	  + (real) 9.0  * (s_scal[smem_idx - 2] - s_scal[smem_idx + 2])
	  + (real) 45.0 * (s_scal[smem_idx + 1] - s_scal[smem_idx - 1])) 
          * d_DIFF1_DX_DIV;
    return res;
}


static __device__ __forceinline__ real der_scaly(const int smem_idx, const real* __restrict__ s_scal)
{
    const real res = ((s_scal[smem_idx + 3*SMEM_WIDTH] - s_scal[smem_idx - 3*SMEM_WIDTH]) 
	  + (real) 9.0  * (s_scal[smem_idx - 2*SMEM_WIDTH] - s_scal[smem_idx + 2*SMEM_WIDTH]) 
	  + (real) 45.0 * (s_scal[smem_idx + 1*SMEM_WIDTH] - s_scal[smem_idx - 1*SMEM_WIDTH]))
      * d_DIFF1_DY_DIV;
    return res;
}


static __device__ __forceinline__ real der_scalz(const real* __restrict__ r_scal)
{
    const real res = ((r_scal[6] - r_scal[0]) 
	  + (real) 9.0  * (r_scal[1] - r_scal[5]) 
	  + (real) 45.0 * (r_scal[4] - r_scal[2])) 
      * d_DIFF1_DZ_DIV;
    return res;
}


static __device__ __forceinline__ real der2_scalx(const int smem_idx, const real* __restrict__ s_scal)
{
	const real res = ((real) 2.0   * (s_scal[smem_idx - 3] + s_scal[smem_idx + 3])
	                - (real) 27.0  * (s_scal[smem_idx - 2] + s_scal[smem_idx + 2]) 
	                + (real) 270.0 * (s_scal[smem_idx - 1] + s_scal[smem_idx + 1]) 
	                - (real) 490.0 * s_scal[smem_idx])
	                * d_DIFF2_DX_DIV;
    return res;
}



static __device__ __forceinline__ real der2_scaly(const int smem_idx, const real* __restrict__ s_scal)
{
	const real res = ((real) 2.0   * (s_scal[smem_idx - 3*SMEM_WIDTH] + s_scal[smem_idx + 3*SMEM_WIDTH])
	                - (real) 27.0  * (s_scal[smem_idx - 2*SMEM_WIDTH] + s_scal[smem_idx + 2*SMEM_WIDTH]) 
	                + (real) 270.0 * (s_scal[smem_idx - 1*SMEM_WIDTH] + s_scal[smem_idx + 1*SMEM_WIDTH]) 
	                - (real) 490.0 * s_scal[smem_idx])
	                * d_DIFF2_DY_DIV;
    return res;
}


static __device__ __forceinline__ real der2_scalz(const real* __restrict__ r_scal)
{
	const real res = ((real) 2.0   * (r_scal[0] + r_scal[6])
	                - (real) 27.0  * (r_scal[1] + r_scal[5]) 
	                + (real) 270.0 * (r_scal[2] + r_scal[4]) 
	                - (real) 490.0 * r_scal[3])
	                * d_DIFF2_DZ_DIV;
    return res;
}


static __device__ __forceinline__ real der2_scalxy(const int smem_idx, const real* __restrict__ s_scal)
{
	const real res = (
	  (real) 2.0   * (   s_scal[smem_idx - 3 - 3*SMEM_WIDTH]
                        -s_scal[smem_idx + 3 - 3*SMEM_WIDTH]
                        +s_scal[smem_idx + 3 + 3*SMEM_WIDTH]
                        -s_scal[smem_idx - 3 + 3*SMEM_WIDTH])
	- (real) 27.0  * (   s_scal[smem_idx - 2 - 2*SMEM_WIDTH]
                        -s_scal[smem_idx + 2 - 2*SMEM_WIDTH]
                        +s_scal[smem_idx + 2 + 2*SMEM_WIDTH]
                        -s_scal[smem_idx - 2 + 2*SMEM_WIDTH])
	+ (real) 270.0 * (   s_scal[smem_idx - 1 - 1*SMEM_WIDTH]
                        -s_scal[smem_idx + 1 - 1*SMEM_WIDTH]
                        +s_scal[smem_idx + 1 + 1*SMEM_WIDTH]
                        -s_scal[smem_idx - 1 + 1*SMEM_WIDTH])
	)* d_DIFFMN_DXDY_DIV;
	return res;
}


template<int level>
static __device__ __forceinline__ real der2_scalxz(const int smem_idx, const real* __restrict__ s_scal)
{
    if (level == 1) {
        const real res = (real) 270.0 * (s_scal[smem_idx - 1] - s_scal[smem_idx + 1]) * d_DIFFMN_DXDZ_DIV;
        return res;
    } else if (level == 2) {
        const real res = (real) -27.0 * (s_scal[smem_idx - 2] - s_scal[smem_idx + 2]) * d_DIFFMN_DXDZ_DIV;
        return res;
    } else {
        const real res = (real)   2.0 * (s_scal[smem_idx - 3] - s_scal[smem_idx + 3]) * d_DIFFMN_DXDZ_DIV;
        return res;
    }
}


template<int level>
static __device__ __forceinline__ real der2_scalyz(const int smem_idx, const real* __restrict__ s_scal)
{
    if (level == 1) {
        const real res = (real) 270.0 * (s_scal[smem_idx - 1*SMEM_WIDTH] - s_scal[smem_idx + 1*SMEM_WIDTH]) * d_DIFFMN_DYDZ_DIV;
        return res;
    } else if (level == 2) {
        const real res = (real) -27.0 * (s_scal[smem_idx - 2*SMEM_WIDTH] - s_scal[smem_idx + 2*SMEM_WIDTH]) * d_DIFFMN_DYDZ_DIV;
        return res;
    } else {
        const real res = (real)   2.0 * (s_scal[smem_idx - 3*SMEM_WIDTH] - s_scal[smem_idx + 3*SMEM_WIDTH]) * d_DIFFMN_DYDZ_DIV;
        return res;
    }
}

