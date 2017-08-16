/*
*   Reduction operations
*   This all seems intimidating, but it's pretty directly adapted from
*   http://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/reduction/doc/reduction.pdf
*/
#include "collectiveops_cuda_generic.cuh"
#include "gpu/cuda/core/dconsts_core.cuh"
#include "gpu/cuda/core/errorhandler_cuda.cuh"
#include "common/errorhandler.h"
#include "utils/utils.h"//For templated max/min/sum

#define COL_THREADS_X (32) //TODO read from config
#define COL_THREADS_Y (8)
#define COL_ELEMS_PER_THREAD (8) //TODO ifndef, define here, else get from the compile flag

static real* d_vec_res          = NULL;
static real* d_partial_result   = NULL;
static bool is_initialized      = false;
static CParamConfig* cparams    = NULL;


//Comparison funcs (template magic :D)
template<class T>
__device__ T dmax(const T a, const T b) { return a > b ? a : b; }

template<class T>
__device__ T dmin(const T a, const T b) { return a < b ? a : b; }

template<class T>
__device__ T dsum(const T a, const T b) { return a + b; }


//The initial values when starting the reduction
template<class T>
__device__ T ddist(const T a, const T b, const T c) { return sqrt(a*a + b*b + c*c); }

template<class T>
__device__ T dsqrsum(const T a, const T b, const T c) { return a*a + b*b + c*c; }

template<class T>
__device__ T dexpsqrscal(const T a, const T b, const T c) { return exp(a)*exp(a); }

template<class T>
__device__ T dscal(const T a, const T b, const T c) { return a; }

//Function pointer definitions
typedef real (*ReduceFunc)(real, real);
typedef real (*ReduceInitFunc)(real, real, real);


template<ReduceFunc reduce_op>
__device__ void reduce_warp(volatile real* shared, int tid)
{
    shared[tid] = reduce_op(shared[tid], shared[tid+32]);__syncthreads();
    shared[tid] = reduce_op(shared[tid], shared[tid+16]);__syncthreads();
    shared[tid] = reduce_op(shared[tid], shared[tid+8]);__syncthreads();
    shared[tid] = reduce_op(shared[tid], shared[tid+4]);__syncthreads();
    shared[tid] = reduce_op(shared[tid], shared[tid+2]);__syncthreads();
    shared[tid] = reduce_op(shared[tid], shared[tid+1]);__syncthreads();
}


void init_collectiveops_cuda_generic(CParamConfig* conf)
{
    if (is_initialized) CRASH("collectiveops_cuda_generic() already initialized!");
    is_initialized = true;

    cparams = conf;

    // The device variable where the maximum found value found is written to
    CUDA_ERRCHK( cudaMalloc((real**) &d_vec_res, sizeof(real)) );
 
    //An intermediate array used in computing the reduction
	dim3 blocks_per_grid;
	blocks_per_grid.x = ceil(cparams->nx / (double)COL_THREADS_X);
	blocks_per_grid.y = ceil(cparams->ny / (double)COL_THREADS_Y);
	blocks_per_grid.z = ceil(cparams->nz / (double)COL_ELEMS_PER_THREAD);
	const int blocks_total = blocks_per_grid.x * blocks_per_grid.y * blocks_per_grid.z;
	CUDA_ERRCHK( cudaMalloc((real**) &d_partial_result, sizeof(real)*blocks_total) );
}


void destroy_collectiveops_cuda_generic()
{
    if (!is_initialized) CRASH("collectiveops_cuda_generic() wasn't initialized!");
    is_initialized = false;

    cparams = NULL;
    CUDA_ERRCHK( cudaFree(d_vec_res) ); d_vec_res = NULL;
    CUDA_ERRCHK( cudaFree(d_partial_result) ); d_partial_result = NULL;
}


template <unsigned int block_size, ReduceFunc reduce_op>
__global__ void reduce(real* dest, real* src, int problem_size)
{
	int tid = threadIdx.x;
	int i = tid + blockDim.x;//tid + offset
	extern __shared__ real shared[];
	shared[tid] = src[tid];
	
	//Add sequantially all blocks above block size
	while (i < problem_size) {
		shared[tid] = reduce_op(shared[tid], src[i]);
		i += blockDim.x;
	}
	__syncthreads();

	if (block_size >= 1024) {
		if (tid < 512) { 			
			shared[tid] = reduce_op(shared[tid], shared[tid+512]); 
			__syncthreads(); 
		}
	}
	if (block_size >= 512) {
		if (tid < 256) { 			
			shared[tid] = reduce_op(shared[tid], shared[tid+256]); 
			__syncthreads(); 
		}
	}
	if (block_size >= 256) {
		if (tid < 128) { 			
			shared[tid] = reduce_op(shared[tid], shared[tid+128]); 
			__syncthreads(); 
		}
	}
	if (block_size >= 128) {
		if (tid < 64) { 			
			shared[tid] = reduce_op(shared[tid], shared[tid+64]); 
			__syncthreads(); 
		}
	}
	if (tid < 32)
		reduce_warp<reduce_op>(shared, tid);		

	if (tid == 0) *dest = shared[0];
}


//Calculates the maximum vector magnitude found in the system
template <unsigned int block_size, ReduceFunc reduce_op, ReduceInitFunc reduce_init_op>//inline const T& min(const T& a, const T& b)
__global__ void reduce_initial(real* d_partial_result, real* d_vec_x, real* d_vec_y = NULL, real* d_vec_z = NULL)
{
	int tid, grid_idx;
	extern __shared__ real vec_shared[];

	tid = threadIdx.x + threadIdx.y*blockDim.x;//index inside the shared mem array
	grid_idx = 	(threadIdx.x + blockIdx.x*blockDim.x + d_nx_min) 	+	 
			(threadIdx.y + blockIdx.y*blockDim.y + d_ny_min)*d_mx +
			(blockIdx.z*COL_ELEMS_PER_THREAD 	     + d_nz_min)*d_mx*d_my;

	
	
	real vec;

    const bool REDUCE_VEC = (d_vec_y != NULL);
    if (REDUCE_VEC)
	    vec_shared[tid] = reduce_init_op(d_vec_x[grid_idx], 
                                         d_vec_y[grid_idx], 
                                         d_vec_z[grid_idx]); //init first value
    else
        vec_shared[tid] = reduce_init_op(d_vec_x[grid_idx], 0, 0);

	for (int i=1; i < COL_ELEMS_PER_THREAD; i++)
	{	
		grid_idx += d_mx*d_my;
        if (REDUCE_VEC)
	        vec = reduce_init_op(d_vec_x[grid_idx], 
                                             d_vec_y[grid_idx], 
                                             d_vec_z[grid_idx]); //init first value
        else
            vec = reduce_init_op(d_vec_x[grid_idx], 0, 0);
	
		vec_shared[tid] = reduce_op(vec_shared[tid], vec);
	}
	__syncthreads();

	if (block_size >= 1024) {
		if (tid < 512) { 			
			vec_shared[tid] = reduce_op(vec_shared[tid], vec_shared[tid+512]); 
			__syncthreads(); 
		}
	}

	if (block_size >= 512) {
		if (tid < 256) { 			
			vec_shared[tid] = reduce_op(vec_shared[tid], vec_shared[tid+256]); 
			__syncthreads(); 
		}
	}
	
	if (block_size >= 256) {
		if (tid < 128) { 			
			vec_shared[tid] = reduce_op(vec_shared[tid], vec_shared[tid+128]); 
			__syncthreads(); 
		}
	}

	if (tid < 64) 
	{ 			
		vec_shared[tid] = reduce_op(vec_shared[tid], vec_shared[tid+64]); 
		__syncthreads(); 
	}

	if (tid < 32)
		reduce_warp<reduce_op>(vec_shared, tid);

	if (tid == 0) d_partial_result[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = vec_shared[0];
}


/*
* Calculates the max vec found in the grid
* Puts the result in d_vec_res (in device memory)
*/
template<ReduceFunc reduce_op, ReduceInitFunc reduce_init_op>
void reduce_cuda_generic(real* d_vec_x, real* d_vec_y = NULL, real* d_vec_z = NULL)
{
    if (!is_initialized) CRASH("collectiveops_cuda_generic() wasn't initialized!");

	//-------------------------------------------------
	dim3 threads_per_block, blocks_per_grid;
	threads_per_block.x = COL_THREADS_X;
	threads_per_block.y = COL_THREADS_Y;
	threads_per_block.z = 1; // 2D blockdims only 
	const int SMEM_PER_BLOCK = threads_per_block.x * threads_per_block.y *threads_per_block.z * sizeof(real);

	blocks_per_grid.x = ceil((real) cparams->nx / (real)COL_THREADS_X);
	blocks_per_grid.y = ceil((real) cparams->ny / (real)COL_THREADS_Y);
	blocks_per_grid.z = ceil((real) cparams->nz / (real)COL_ELEMS_PER_THREAD);
	const int BLOCKS_TOTAL = blocks_per_grid.x * blocks_per_grid.y * blocks_per_grid.z;
	//------------------------------------------------
	
	switch (threads_per_block.x*threads_per_block.y*threads_per_block.z)
	{
		case 512:
			reduce_initial<512, reduce_op, reduce_init_op><<<blocks_per_grid, threads_per_block, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z); CUDA_ERRCHK_KERNEL(); break;
		case 256:
			reduce_initial<256, reduce_op, reduce_init_op><<<blocks_per_grid, threads_per_block, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z); CUDA_ERRCHK_KERNEL(); break;
		case 128:
			reduce_initial<128, reduce_op, reduce_init_op><<<blocks_per_grid, threads_per_block, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z); CUDA_ERRCHK_KERNEL(); break;
		default:
			printf("INCORRECT THREAD SIZE!\n");
			exit(EXIT_FAILURE);
	}
    if (BLOCKS_TOTAL >= 1024) {
        reduce<1024, reduce_op><<<1, 1024, 1024*sizeof(real)>>>(d_vec_res, d_partial_result, BLOCKS_TOTAL); CUDA_ERRCHK_KERNEL();
    } else if (BLOCKS_TOTAL >= 512) {
        reduce<512, reduce_op><<<1, 512, 512*sizeof(real)>>>(d_vec_res, d_partial_result, BLOCKS_TOTAL); CUDA_ERRCHK_KERNEL();
    } else if (BLOCKS_TOTAL >= 256) {
        reduce<256, reduce_op><<<1, 256, 256*sizeof(real)>>>(d_vec_res, d_partial_result, BLOCKS_TOTAL); CUDA_ERRCHK_KERNEL();
    } else if (BLOCKS_TOTAL >= 128) {
        reduce<128, reduce_op><<<1, 128, 128*sizeof(real)>>>(d_vec_res, d_partial_result, BLOCKS_TOTAL); CUDA_ERRCHK_KERNEL();
    } else if (BLOCKS_TOTAL >= 16) {
        reduce<16, reduce_op><<<1, 16, 16*sizeof(real)>>>(d_vec_res, d_partial_result, BLOCKS_TOTAL); CUDA_ERRCHK_KERNEL();
    } else {
        printf("INCORRECT BLOCKS_TOTAL (= %d) IN collectiveops.cu!\n", BLOCKS_TOTAL);
				exit(EXIT_FAILURE);
    }
    cudaDeviceSynchronize();
	//We're done
	return;
}


//template<ReductType t>
real get_reduction_cuda_generic(ReductType t, real* d_a, real* d_b, real* d_c)
{
    real res;
    switch (t) {
        case MAX_VEC:
            reduce_cuda_generic<dmax, ddist>(d_a, d_b, d_c);
            break;
        case MIN_VEC:
            reduce_cuda_generic<dmin, ddist>(d_a, d_b, d_c);
            break;
        case RMS_VEC:
            reduce_cuda_generic<dsum, dsqrsum>(d_a, d_b, d_c);
            break;
        case MAX_SCAL:
            reduce_cuda_generic<dmax, dscal>(d_a);
            break;
        case MIN_SCAL:
            reduce_cuda_generic<dmin, dscal>(d_a);
            break;
        case RMS_SCAL:
            reduce_cuda_generic<dsum, dsqrsum>(d_a);
            break;
        case RMS_EXP:
            reduce_cuda_generic<dsum, dexpsqrscal>(d_a);
            break;
        default:
            CRASH("Invalid type!");
    }
    cudaDeviceSynchronize();
    CUDA_ERRCHK( cudaMemcpy(&res, (real*)d_vec_res, sizeof(real), cudaMemcpyDeviceToHost) );
    cudaDeviceSynchronize();
    return res;
}











































































































