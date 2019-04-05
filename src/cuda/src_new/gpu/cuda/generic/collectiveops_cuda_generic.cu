/*
*   Reduction operations
*   This all seems intimidating, but it's pretty directly adapted from
*   http://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/reduction/doc/reduction.pdf
*/
#include "collectiveops_cuda_generic.cuh"
#include "gpu/cuda/core/dconsts_core.cuh"
#include "gpu/cuda/core/errorhandler_cuda.cuh"
#include "utils/utils.h"         //For templated max/min/sum

#define COL_THREADS_X (32)       //TODO read from config
#define COL_THREADS_Y (4)
#define COL_ELEMS_PER_THREAD (1) //TODO ifndef, define here, else get from the compile flag


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
__device__ T dexpscal(const T a, const T b, const T c) { return exp(a); }

template<class T>
__device__ T dsqrscal(const T a, const T b, const T c) { return a*a; }

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


void init_reduction_array_cuda_generic(ReductionArray* reduct_arr, CParamConfig* cparams)
{
    // The device variable where the maximum found value found is written to
    CUDA_ERRCHK( cudaMalloc((real**) &reduct_arr->d_vec_res, sizeof(real)) );
 
    //An intermediate array used in computing the reduction
	dim3 bpg;
	bpg.x = ceil(cparams->nx / (double)COL_THREADS_X);
	bpg.y = ceil(cparams->ny / (double)COL_THREADS_Y);
	bpg.z = ceil(cparams->nz / (double)COL_ELEMS_PER_THREAD);
	const int blocks_total = bpg.x * bpg.y * bpg.z;
	CUDA_ERRCHK( cudaMalloc((real**) &reduct_arr->d_partial_result, sizeof(real)*blocks_total) );
}


void destroy_reduction_array_cuda_generic(ReductionArray* reduct_arr)
{
    CUDA_ERRCHK( cudaFree(reduct_arr->d_vec_res) ); reduct_arr->d_vec_res = NULL;
    CUDA_ERRCHK( cudaFree(reduct_arr->d_partial_result) ); reduct_arr->d_partial_result = NULL;
}


template <unsigned int block_size, ReduceFunc reduce_op>
__global__ void reduce(real* dest, real* src, int problem_size)
{
	int tid = threadIdx.x;
	int i = tid + blockDim.x;//tid + offset
	extern __shared__ real shared[];
	shared[tid] = src[tid];
	
	//Add sequentially all blocks above block size
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
	extern __shared__ real vec_shared[];

	const int tid = threadIdx.x + threadIdx.y*blockDim.x;//index inside the shared mem array

    const int tx = threadIdx.x + blockIdx.x*blockDim.x + d_nx_min;
    const int ty = threadIdx.y + blockIdx.y*blockDim.y + d_ny_min;
    const int tz = blockIdx.z*COL_ELEMS_PER_THREAD     + d_nz_min;

	const int base_idx = tx + ty*d_mx + tz*d_mxy;

    assert(tx >= d_nx_min);
	assert(tx < d_nx_max);
    assert(ty >= d_ny_min);
    assert(ty < d_ny_max);
    assert(tz >= d_nz_min);
    assert(tz < d_nz_max);
	
	real vec;

    const bool REDUCE_VEC = (d_vec_y != NULL);
    if (REDUCE_VEC)
	    vec_shared[tid] = reduce_init_op(d_vec_x[base_idx], 
                                         d_vec_y[base_idx], 
                                         d_vec_z[base_idx]); //init first value
    else
        vec_shared[tid] = reduce_init_op(d_vec_x[base_idx], 0, 0);

	for (int i=1; i < COL_ELEMS_PER_THREAD && tz+i < d_nz_max; i++)
	{	
		const int grid_idx = base_idx + i*d_mxy;
        assert(tz+i < d_nz_max);
        
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
    if (block_size >= 128) {
	    if (tid < 64) { 			
		    vec_shared[tid] = reduce_op(vec_shared[tid], vec_shared[tid+64]); 
		    __syncthreads(); 
	    }
    }

	if (tid < 32)
		reduce_warp<reduce_op>(vec_shared, tid);

	if (tid == 0) d_partial_result[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = vec_shared[0];
}


/*
* Calculates the max vec found in the grid
* Puts the result in reduct_arr->d_vec_res (in device memory)
*/
template<ReduceFunc reduce_op, ReduceInitFunc reduce_init_op>
void reduce_cuda_generic(ReductionArray* reduct_arr, CParamConfig* cparams, real* d_vec_x, real* d_vec_y = NULL, real* d_vec_z = NULL)
{
    const dim3 tpb(COL_THREADS_X, COL_THREADS_Y, 1);
    const dim3 bpg((unsigned int) ceil((real) cparams->nx / (real)COL_THREADS_X),
                   (unsigned int) ceil((real) cparams->ny / (real)COL_THREADS_Y),
                   (unsigned int) ceil((real) cparams->nz / (real)COL_ELEMS_PER_THREAD));

    const size_t SMEM_PER_BLOCK = tpb.x * tpb.y *tpb.z * sizeof(real);
    const int BLOCKS_TOTAL = bpg.x * bpg.y * bpg.z;


    //Collectiveops works only when BLOCKS_TOTAL is divisible by the thread block size.
    //This is not good and collectiveops should be rewritten to support arbitrary grid dims
    if (BLOCKS_TOTAL % (tpb.x * tpb.y * tpb.z) != 0)
        CRASH("Incorrect BLOCKS_TOTAL in reduce_cuda_generic()")
	
	switch (tpb.x*tpb.y*tpb.z)
	{
		case 512:
			reduce_initial<512, reduce_op, reduce_init_op><<<bpg, tpb, SMEM_PER_BLOCK>>>(reduct_arr->d_partial_result, d_vec_x, d_vec_y, d_vec_z); CUDA_ERRCHK_KERNEL(); break;
		case 256:
			reduce_initial<256, reduce_op, reduce_init_op><<<bpg, tpb, SMEM_PER_BLOCK>>>(reduct_arr->d_partial_result, d_vec_x, d_vec_y, d_vec_z); CUDA_ERRCHK_KERNEL(); break;
		case 128:
			reduce_initial<128, reduce_op, reduce_init_op><<<bpg, tpb, SMEM_PER_BLOCK>>>(reduct_arr->d_partial_result, d_vec_x, d_vec_y, d_vec_z); CUDA_ERRCHK_KERNEL(); break;
		default:
			printf("INCORRECT THREAD SIZE!\n");
			exit(EXIT_FAILURE);
	}
    if (BLOCKS_TOTAL >= 1024) {
        reduce<1024, reduce_op><<<1, 1024, 1024*sizeof(real)>>>(reduct_arr->d_vec_res, reduct_arr->d_partial_result, BLOCKS_TOTAL); CUDA_ERRCHK_KERNEL();
    } else if (BLOCKS_TOTAL >= 512) {
        reduce<512, reduce_op><<<1, 512, 512*sizeof(real)>>>(reduct_arr->d_vec_res, reduct_arr->d_partial_result, BLOCKS_TOTAL); CUDA_ERRCHK_KERNEL();
    } else if (BLOCKS_TOTAL >= 256) {
        reduce<256, reduce_op><<<1, 256, 256*sizeof(real)>>>(reduct_arr->d_vec_res, reduct_arr->d_partial_result, BLOCKS_TOTAL); CUDA_ERRCHK_KERNEL();
    } else if (BLOCKS_TOTAL >= 128) {
        reduce<128, reduce_op><<<1, 128, 128*sizeof(real)>>>(reduct_arr->d_vec_res, reduct_arr->d_partial_result, BLOCKS_TOTAL); CUDA_ERRCHK_KERNEL();
    } else if (BLOCKS_TOTAL >= 16) {
        reduce<16, reduce_op><<<1, 16, 16*sizeof(real)>>>(reduct_arr->d_vec_res, reduct_arr->d_partial_result, BLOCKS_TOTAL); CUDA_ERRCHK_KERNEL();
    } else {
        printf("INCORRECT BLOCKS_TOTAL (= %d) IN collectiveops.cu!\n", BLOCKS_TOTAL);
				exit(EXIT_FAILURE);
    }

	return;
}

//template<ReductType t>
real get_reduction_cuda_generic(ReductionArray* reduct_arr, ReductType t, CParamConfig* cparams, real* d_a, real* d_b, real* d_c)
{
    real res;
    switch (t) {
        case MAX_VEC:
            reduce_cuda_generic<dmax, ddist>(reduct_arr, cparams, d_a, d_b, d_c);
            break;
        case MIN_VEC:
            reduce_cuda_generic<dmin, ddist>(reduct_arr, cparams, d_a, d_b, d_c);
            break;
        case RMS_VEC:
            reduce_cuda_generic<dsum, dsqrsum>(reduct_arr, cparams, d_a, d_b, d_c);
            break;
        case MAX_SCAL:
            reduce_cuda_generic<dmax, dscal>(reduct_arr, cparams, d_a);
            break;
        case MIN_SCAL:
            reduce_cuda_generic<dmin, dscal>(reduct_arr, cparams, d_a);
            break;
        case RMS_SCAL:
            reduce_cuda_generic<dsum, dsqrscal>(reduct_arr, cparams, d_a);
            break;
        case RMS_EXP:
            reduce_cuda_generic<dsum, dexpsqrscal>(reduct_arr, cparams, d_a);
            break;
        case SUM_SCAL:
            reduce_cuda_generic<dsum, dscal>(reduct_arr, cparams, d_a);
            break;
        case SUM_EXP:
            reduce_cuda_generic<dsum, dexpscal>(reduct_arr, cparams, d_a);
            break;
        default:
            CRASH("Invalid type!");
    }
    CUDA_ERRCHK( cudaMemcpy(&res, (real*)reduct_arr->d_vec_res, sizeof(real), cudaMemcpyDeviceToHost) );

    return res;
}
