
#include <stdio.h>

#define EXTERN extern
#include "dconsts.cuh"
#include "../cparam_c.h"
#include "smem.cuh"

//DEBUG
#include "diagnostics.cuh"

//This all seems intimidating, but it's pretty directly adapted from
//http://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/reduction/doc/reduction.pdf


__device__ void reduce_warp_max(volatile float* max_shared, int tid)
{
	//Find the maximum within a warp. 
	//MIIKKA's note: As all threads should be executed in the warp before proceeding, syncronization is not an issue?
	if (max_shared[tid] < max_shared[tid+32])
		max_shared[tid] = max_shared[tid+32]; 
	if (max_shared[tid] < max_shared[tid+16])
		max_shared[tid] = max_shared[tid+16]; 
	if (max_shared[tid] < max_shared[tid+8])
		max_shared[tid] = max_shared[tid+8]; 
	if (max_shared[tid] < max_shared[tid+4])
		max_shared[tid] = max_shared[tid+4]; 
	if (max_shared[tid] < max_shared[tid+2])
		max_shared[tid] = max_shared[tid+2]; 
	if (max_shared[tid] < max_shared[tid+1])
		max_shared[tid] = max_shared[tid+1]; 
}

__device__ void reduce_warp_min(volatile float* min_shared, int tid)
{
	//Find the minumum within a warp. 
	if (min_shared[tid] > min_shared[tid+32])
		min_shared[tid] = min_shared[tid+32]; 
	if (min_shared[tid] > min_shared[tid+16])
		min_shared[tid] = min_shared[tid+16]; 
	if (min_shared[tid] > min_shared[tid+8])
		min_shared[tid] = min_shared[tid+8]; 
	if (min_shared[tid] > min_shared[tid+4])
		min_shared[tid] = min_shared[tid+4]; 
	if (min_shared[tid] > min_shared[tid+2])
		min_shared[tid] = min_shared[tid+2]; 
	if (min_shared[tid] > min_shared[tid+1])
		min_shared[tid] = min_shared[tid+1]; 
}

__device__ void reduce_warp_sum(volatile float* sum_shared, int tid)
{
	//Sum all values within a warp. 
	//MIIKKA's note: As all threads should be executed in the warp before proceeding, syncronization is not an issue?
	sum_shared[tid] += sum_shared[tid+32]; 
	sum_shared[tid] += sum_shared[tid+16]; 
	sum_shared[tid] += sum_shared[tid+8]; 
	sum_shared[tid] += sum_shared[tid+4]; 
	sum_shared[tid] += sum_shared[tid+2]; 
	sum_shared[tid] += sum_shared[tid+1]; 
}


//assumes 512/1024 threads 
//TODO more flexible
template <unsigned int block_size>
__global__ void reduce_max(float* dest, float* src, int problem_size)
{
	int tid = threadIdx.x;
	int i = tid + blockDim.x;//tid + offset
	extern __shared__ float max_shared[];
	max_shared[tid] = src[tid];
	
	//Add sequantially all blocks above block size
	while (i < problem_size) {
		if (max_shared[tid] < src[i])
			max_shared[tid] = src[i];
		i += blockDim.x;
	}
	__syncthreads();

	//TODO make sure works with all kinds of nx,ny,nz combinations
	if (block_size >= 1024) {
		if (tid < 512) { 			
			if (max_shared[tid] < max_shared[tid+512])
				max_shared[tid] = max_shared[tid+512]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 512) {
		if (tid < 256) { 			
			if (max_shared[tid] < max_shared[tid+256])
				max_shared[tid] = max_shared[tid+256]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 256) {
		if (tid < 128) { 			
			if (max_shared[tid] < max_shared[tid+128])
				max_shared[tid] = max_shared[tid+128]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 128) {
		if (tid < 64) { 			
			if (max_shared[tid] < max_shared[tid+64])
				max_shared[tid] = max_shared[tid+64]; 
			__syncthreads(); 
		}
	}

	if (tid < 32)
		reduce_warp_max(max_shared, tid);
		
	if (tid == 0) *dest = max_shared[0];
}

//assumes 512/1024 threads 
//TODO more flexible
template <unsigned int block_size>
__global__ void reduce_min(float* dest, float* src, int problem_size)
{
	int tid = threadIdx.x;
	int i = tid + blockDim.x;//tid + offset
	extern __shared__ float min_shared[];
	min_shared[tid] = src[tid];
	
	//Add sequantially all blocks above block size
	while (i < problem_size) {
		if (min_shared[tid] > src[i])
			min_shared[tid] = src[i];
		i += blockDim.x;
	}
	__syncthreads();

	//TODO make sure works with all kinds of nx,ny,nz combinations
	if (block_size >= 1024) {
		if (tid < 512) { 			
			if (min_shared[tid] > min_shared[tid+512])
				min_shared[tid] = min_shared[tid+512]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 512) {
		if (tid < 256) { 			
			if (min_shared[tid] > min_shared[tid+256])
				min_shared[tid] = min_shared[tid+256]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 256) {
		if (tid < 128) { 			
			if (min_shared[tid] > min_shared[tid+128])
				min_shared[tid] = min_shared[tid+128]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 128) {
		if (tid < 64) { 			
			if (min_shared[tid] > min_shared[tid+64])
				min_shared[tid] = min_shared[tid+64]; 
			__syncthreads(); 
		}
	}

	if (tid < 32)
		reduce_warp_min(min_shared, tid);
		
	if (tid == 0) *dest = min_shared[0];
}

//assumes 512/1024 threads 
//TODO more flexible
template <unsigned int block_size>
__global__ void reduce_rms(float* dest, float* src, int problem_size, bool root=true)
{
	int tid = threadIdx.x;
	int i = tid + blockDim.x;//tid + offset
	extern __shared__ float sum_shared[];

	sum_shared[tid] = src[tid];
	
	//Add sequentially all blocks above block size
	while (i < problem_size) {
		sum_shared[tid] += src[i];
		i += blockDim.x;
	}
	__syncthreads();

	//TODO make sure works with all kinds of nx,ny,nz combinations
	if (block_size >= 1024) {
		if (tid < 512) { 			
			sum_shared[tid] += sum_shared[tid+512]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 512) {
		if (tid < 256) { 			
			sum_shared[tid] += sum_shared[tid+256]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 256) {
		if (tid < 128) { 			
			sum_shared[tid] += sum_shared[tid+128]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 128) {
		if (tid < 64) { 			
			sum_shared[tid] += sum_shared[tid+64]; 
			__syncthreads(); 
		}
	}

	if (tid < 32)
	{
		reduce_warp_sum(sum_shared, tid);
		__syncthreads(); 
	}

	//Finishes the mean of scal^2 and calculates the sqrt of it to get rms.
	if (tid == 0) 
	{
		if (root) 
        	  *dest = sqrtf( sum_shared[0]/d_NELEMENTS_FLOAT );
     		else
        	  *dest = sum_shared[0];
	}
}

//Calculates the maximum vector magnitude found in the system
template <unsigned int block_size>
__global__ void max_vec(float* d_partial_result, float* d_vec_x, float* d_vec_y, float* d_vec_z, const int elementsPerThread)
{
	int tid, grid_idx;
	extern __shared__ float max_vec_shared[];

	tid = threadIdx.x + threadIdx.y*blockDim.x;//index inside the shared mem array
	grid_idx = 	(threadIdx.x + blockIdx.x*blockDim.x + d_CX_BOT) 	+	 
			(threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT)*d_NX +
			(blockIdx.z*elementsPerThread 	     + d_CZ_BOT)*d_NX*d_NY;

	
	
	float vec;
	max_vec_shared[tid] = sqrtf(	d_vec_x[grid_idx]*d_vec_x[grid_idx]   + 
					d_vec_y[grid_idx]*d_vec_y[grid_idx] + 
					d_vec_z[grid_idx]*d_vec_z[grid_idx] ); //init first value
	for (int i=1; i < elementsPerThread; i++)
	{	
		grid_idx += d_NX*d_NY;
		vec = sqrtf(	d_vec_x[grid_idx]*d_vec_x[grid_idx]   + 
				d_vec_y[grid_idx]*d_vec_y[grid_idx] + 
				d_vec_z[grid_idx]*d_vec_z[grid_idx] );
	
		if (vec > max_vec_shared[tid])
			max_vec_shared[tid] = vec;
	}
	__syncthreads();

	if (block_size >= 1024) {
		if (tid < 512) { 			
			if (max_vec_shared[tid] < max_vec_shared[tid+512])
				max_vec_shared[tid] = max_vec_shared[tid+512]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 512) {
		if (tid < 256) { 			
			if (max_vec_shared[tid] < max_vec_shared[tid+256])
				max_vec_shared[tid] = max_vec_shared[tid+256]; 
			__syncthreads(); 
		}
	}
	
	if (block_size >= 256) {
		if (tid < 128) { 			
			if (max_vec_shared[tid] < max_vec_shared[tid+128])
				max_vec_shared[tid] = max_vec_shared[tid+128]; 
			__syncthreads(); 
		}
	}

	if (tid < 64) 
	{ 			
		if (max_vec_shared[tid] < max_vec_shared[tid+64])
			max_vec_shared[tid] = max_vec_shared[tid+64]; 
		__syncthreads(); 
	}

	if (tid < 32)
		reduce_warp_max(max_vec_shared, tid);

	if (tid == 0) d_partial_result[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = max_vec_shared[0];
}

//Calculates the minimum vector magnitude found in the system
template <unsigned int block_size>
__global__ void min_vec(float* d_partial_result, float* d_vec_x, float* d_vec_y, float* d_vec_z, const int elementsPerThread)
{
	int tid, grid_idx;
	extern __shared__ float min_vec_shared[];

	tid = threadIdx.x + threadIdx.y*blockDim.x;//index inside the shared mem array
	grid_idx = (threadIdx.x + blockIdx.x*blockDim.x + d_CX_BOT) 	 +	 
		   (threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT)*d_NX +
	       	   (blockIdx.z*elementsPerThread 	+ d_CZ_BOT)*d_NX*d_NY;

	
	
	float vec;
	min_vec_shared[tid] = sqrtf(d_vec_x[grid_idx]*d_vec_x[grid_idx] + 
				    d_vec_y[grid_idx]*d_vec_y[grid_idx] + 
				    d_vec_z[grid_idx]*d_vec_z[grid_idx] ); //init first value
	for (int i=1; i < elementsPerThread; i++)
	{	
		grid_idx += d_NX*d_NY;
		vec = sqrtf(d_vec_x[grid_idx]*d_vec_x[grid_idx] + 
			    d_vec_y[grid_idx]*d_vec_y[grid_idx] + 
		   	    d_vec_z[grid_idx]*d_vec_z[grid_idx] );
	
		if (vec < min_vec_shared[tid])
			min_vec_shared[tid] = vec;
	}
	__syncthreads();

	if (block_size >= 1024) {
		if (tid < 512) { 			
			if (min_vec_shared[tid] > min_vec_shared[tid+512])
				min_vec_shared[tid] = min_vec_shared[tid+512]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 512) {
		if (tid < 256) { 			
			if (min_vec_shared[tid] > min_vec_shared[tid+256])
				min_vec_shared[tid] = min_vec_shared[tid+256]; 
			__syncthreads(); 
		}
	}
	
	if (block_size >= 256) {
		if (tid < 128) { 			
			if (min_vec_shared[tid] > min_vec_shared[tid+128])
				min_vec_shared[tid] = min_vec_shared[tid+128]; 
			__syncthreads(); 
		}
	}

	if (tid < 64) 
	{ 			
		if (min_vec_shared[tid] > min_vec_shared[tid+64])
			min_vec_shared[tid] = min_vec_shared[tid+64]; 
		__syncthreads(); 
	}

	if (tid < 32)
		reduce_warp_min(min_vec_shared, tid);

	if (tid == 0) d_partial_result[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = min_vec_shared[0];
}


//Calculates the maximum value from a scalar array
template <unsigned int block_size>
__global__ void max_scal(float* d_partial_result, float* d_scal, const int elementsPerThread)
{
	int tid, grid_idx;
	extern __shared__ float max_scal_shared[];

	tid = threadIdx.x + threadIdx.y*blockDim.x;//index inside the shared mem array
	grid_idx = 	(threadIdx.x + blockIdx.x*blockDim.x + d_CX_BOT) 	+	 
			(threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT)*d_NX +
			(blockIdx.z*elementsPerThread 	     + d_CZ_BOT)*d_NX*d_NY;

	
	
	float scal;
	max_scal_shared[tid] = d_scal[grid_idx] ; //init first value
	for (int i=1; i < elementsPerThread; i++)
	{	
		grid_idx += d_NX*d_NY;
		scal = d_scal[grid_idx];
	
		if (scal > max_scal_shared[tid])
			max_scal_shared[tid] = scal;
	}
	__syncthreads();

	if (block_size >= 1024) {
		if (tid < 512) { 			
			if (max_scal_shared[tid] < max_scal_shared[tid+512])
				max_scal_shared[tid] = max_scal_shared[tid+512]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 512) {
		if (tid < 256) { 			
			if (max_scal_shared[tid] < max_scal_shared[tid+256])
				max_scal_shared[tid] = max_scal_shared[tid+256]; 
			__syncthreads(); 
		}
	}
	
	if (block_size >= 256) {
		if (tid < 128) { 			
			if (max_scal_shared[tid] < max_scal_shared[tid+128])
				max_scal_shared[tid] = max_scal_shared[tid+128]; 
			__syncthreads(); 
		}
	}

	if (tid < 64) 
	{ 			
		if (max_scal_shared[tid] < max_scal_shared[tid+64])
			max_scal_shared[tid] = max_scal_shared[tid+64]; 
		__syncthreads(); 
	}

	if (tid < 32)
		reduce_warp_max(max_scal_shared, tid);

	if (tid == 0) d_partial_result[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = max_scal_shared[0];
}

//Calculates the minumum value from a scalar array
template <unsigned int block_size>
__global__ void min_scal(float* d_partial_result, float* d_scal, const int elementsPerThread)
{
	int tid, grid_idx;
	extern __shared__ float min_scal_shared[];

	tid = threadIdx.x + threadIdx.y*blockDim.x;//index inside the shared mem array
	grid_idx = 	(threadIdx.x + blockIdx.x*blockDim.x + d_CX_BOT) 	+	 
			(threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT)*d_NX +
			(blockIdx.z*elementsPerThread 	     + d_CZ_BOT)*d_NX*d_NY;

	
	
	float scal;
	min_scal_shared[tid] = d_scal[grid_idx] ; //init first value
	for (int i=1; i < elementsPerThread; i++)
	{	
		grid_idx += d_NX*d_NY;
		scal = d_scal[grid_idx];
	
		if (scal < min_scal_shared[tid])
			min_scal_shared[tid] = scal;
	}
	__syncthreads();

	if (block_size >= 1024) {
		if (tid < 512) { 			
			if (min_scal_shared[tid] > min_scal_shared[tid+512])
				min_scal_shared[tid] = min_scal_shared[tid+512]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 512) {
		if (tid < 256) { 			
			if (min_scal_shared[tid] > min_scal_shared[tid+256])
				min_scal_shared[tid] = min_scal_shared[tid+256]; 
			__syncthreads(); 
		}
	}
	
	if (block_size >= 256) {
		if (tid < 128) { 			
			if (min_scal_shared[tid] > min_scal_shared[tid+128])
				min_scal_shared[tid] = min_scal_shared[tid+128]; 
			__syncthreads(); 
		}
	}

	if (tid < 64) 
	{ 			
		if (min_scal_shared[tid] > min_scal_shared[tid+64])
			min_scal_shared[tid] = min_scal_shared[tid+64]; 
		__syncthreads(); 
	}

	if (tid < 32)
		reduce_warp_min(min_scal_shared, tid);

	if (tid == 0) d_partial_result[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = min_scal_shared[0];
}

//Calculates the root mean square velocity found in the system
template <unsigned int block_size>
__global__ void vec2_sum(float* d_partial_result, float* d_vec_x, float* d_vec_y, float* d_vec_z, const int elementsPerThread)
{
	int tid, grid_idx;
	extern __shared__ float vec2_sum_shared[];

	tid = threadIdx.x + threadIdx.y*blockDim.x;//index inside the shared mem array
	grid_idx = 	(threadIdx.x + blockIdx.x*blockDim.x + d_CX_BOT) 	+	 
			(threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT)*d_NX +
			(blockIdx.z*elementsPerThread 	     + d_CZ_BOT)*d_NX*d_NY;

	
	//Initialize	
	vec2_sum_shared[tid] = d_vec_x[grid_idx]*d_vec_x[grid_idx] +
        		       d_vec_y[grid_idx]*d_vec_y[grid_idx] +
           		       d_vec_z[grid_idx]*d_vec_z[grid_idx];

	for (int i=1; i < elementsPerThread; i++)
	{	
		grid_idx += d_NX*d_NY;
		vec2_sum_shared[tid] += d_vec_x[grid_idx]*d_vec_x[grid_idx] + 
		       			d_vec_y[grid_idx]*d_vec_y[grid_idx] + 
		       			d_vec_z[grid_idx]*d_vec_z[grid_idx];
	}
	__syncthreads();

	if (block_size >= 1024) {
		if (tid < 512) { 			
			vec2_sum_shared[tid] += vec2_sum_shared[tid+512]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 512) {
		if (tid < 256) { 			
			vec2_sum_shared[tid] += vec2_sum_shared[tid+256]; 
			__syncthreads(); 
		}
	}
	
	if (block_size >= 256) {
		if (tid < 128) { 			
			vec2_sum_shared[tid] += vec2_sum_shared[tid+128]; 
			__syncthreads(); 
		}
	}

	if (tid < 64) 
	{ 			
		vec2_sum_shared[tid] += vec2_sum_shared[tid+64]; 
		__syncthreads(); 
	}

	if (tid < 32)
	{
		reduce_warp_sum(vec2_sum_shared, tid);
		__syncthreads(); 
	}

	if (tid == 0) d_partial_result[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = vec2_sum_shared[0];
}

//  Calculates the sum of a scalar or the sum of its squares.

template <unsigned int block_size>
__global__ void scal2_sum(float* d_partial_result, float* d_scal, const int elementsPerThread, bool sqr=true)
{
	int tid, grid_idx;
	extern __shared__ float scal2_sum_shared[];

	tid = threadIdx.x + threadIdx.y*blockDim.x;//index inside the shared mem array
	grid_idx = 	(threadIdx.x + blockIdx.x*blockDim.x + d_CX_BOT) 	+	 
			(threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT)*d_NX +
			(blockIdx.z*elementsPerThread 	     + d_CZ_BOT)*d_NX*d_NY;

        if (sqr) 
		scal2_sum_shared[tid] = d_scal[grid_idx]*d_scal[grid_idx]; //init first values
        else
		scal2_sum_shared[tid] = d_scal[grid_idx]; //init first value

	for (int i=1; i < elementsPerThread; i++)
	{	
		grid_idx += d_NX*d_NY;
	
        	if (sqr) 
			scal2_sum_shared[tid] += d_scal[grid_idx]*d_scal[grid_idx];
		else
			scal2_sum_shared[tid] += d_scal[grid_idx];
	}
	__syncthreads();

	if (block_size >= 1024) {
		if (tid < 512) { 			
			scal2_sum_shared[tid] += scal2_sum_shared[tid+512]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 512) {
		if (tid < 256) { 			
			scal2_sum_shared[tid] += scal2_sum_shared[tid+256]; 
			__syncthreads(); 
		}
	}
	
	if (block_size >= 256) {
		if (tid < 128) { 			
			scal2_sum_shared[tid] += scal2_sum_shared[tid+128]; 
			__syncthreads(); 
		}
	}

	if (tid < 64) 
	{ 			
		scal2_sum_shared[tid] += scal2_sum_shared[tid+64]; 
		__syncthreads(); 
	}

	if (tid < 32)
	{
		reduce_warp_sum(scal2_sum_shared, tid);
		__syncthreads(); 
	}

	if (tid == 0) d_partial_result[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = scal2_sum_shared[0];
}

//Calculates the rms of a exp(scalar)
template <unsigned int block_size>
__global__ void scal2_exp_sum(float* d_partial_result, float* d_scal, const int elementsPerThread, bool sqr=true)
{
	int tid, grid_idx;
	extern __shared__ float scal2_sum_shared[];

	tid = threadIdx.x + threadIdx.y*blockDim.x;//index inside the shared mem array
	grid_idx = 	(threadIdx.x + blockIdx.x*blockDim.x + d_CX_BOT) 	+	 
			(threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT)*d_NX +
			(blockIdx.z*elementsPerThread 	     + d_CZ_BOT)*d_NX*d_NY;
        float tmp;
        tmp=exp(d_scal[grid_idx]);
	if (sqr) 
		scal2_sum_shared[tid] = tmp*tmp; //Initialize
        else
		scal2_sum_shared[tid] = tmp; //Initialize

	for (int i=1; i < elementsPerThread; i++)
	{	
		grid_idx += d_NX*d_NY;
		tmp=exp(d_scal[grid_idx]);	
		if (sqr) 
			scal2_sum_shared[tid] += tmp*tmp;
		else
			scal2_sum_shared[tid] += tmp;
	}
	__syncthreads();

	if (block_size >= 1024) {
		if (tid < 512) { 			
			scal2_sum_shared[tid] += scal2_sum_shared[tid+512]; 
			__syncthreads(); 
		}
	}

	if (block_size >= 512) {
		if (tid < 256) { 			
			scal2_sum_shared[tid] += scal2_sum_shared[tid+256]; 
			__syncthreads(); 
		}
	}
	
	if (block_size >= 256) {
		if (tid < 128) { 			
			scal2_sum_shared[tid] += scal2_sum_shared[tid+128]; 
			__syncthreads(); 
		}
	}

	if (tid < 64) 
	{ 			
		scal2_sum_shared[tid] += scal2_sum_shared[tid+64]; 
		__syncthreads(); 
	}

	if (tid < 32)
	{
		reduce_warp_sum(scal2_sum_shared, tid);
		__syncthreads(); 
	}

	if (tid == 0) d_partial_result[blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*gridDim.x*gridDim.y] = scal2_sum_shared[0];
}


/*
* Calculates the max vec found in the grid
* Puts the result in d_vec_max
* NOTE! 
* 	-User is responsible for synchronizing the threads after calling if needed (cudaDeviceSynchronize())
* 	-Result is put in *device* memory, d_vec_max
*/
void max_vec_cuda(float* d_vec_max, float* d_partial_result, float* d_vec_x, float* d_vec_y, float* d_vec_z)
{
	//-------------------------------------------------
	static dim3 threadsPerBlock, blocksPerGrid;
	threadsPerBlock.x = COL_THREADS_X;
	threadsPerBlock.y = COL_THREADS_Y;
	threadsPerBlock.z = 1; // 2D blockdims only 
	static const int SMEM_PER_BLOCK = threadsPerBlock.x * threadsPerBlock.y *threadsPerBlock.z * sizeof(float);

	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)COL_THREADS_X);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)COL_THREADS_Y);
	blocksPerGrid.z = ceil((float) COMP_DOMAIN_SIZE_Z / (float)COL_ELEMS_PER_THREAD);
	static const int BLOCKS_TOTAL = blocksPerGrid.x * blocksPerGrid.y * blocksPerGrid.z;
	//------------------------------------------------
	
	switch (threadsPerBlock.x*threadsPerBlock.y*threadsPerBlock.z)
	{
		case 512:
			max_vec<512><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z, COL_ELEMS_PER_THREAD); checkKernelErr(); break;
		case 256:
			max_vec<256><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z, COL_ELEMS_PER_THREAD); checkKernelErr(); break;
		case 128:
			max_vec<128><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z, COL_ELEMS_PER_THREAD); checkKernelErr(); break;
		default:
			printf("INCORRECT THREAD SIZE!\n");
			exit(EXIT_FAILURE);
	}
	switch (BLOCKS_TOTAL){
		case 1024:
			reduce_max<1024><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_max, d_partial_result, BLOCKS_TOTAL);
			checkKernelErr();			
			break;
		case 512:
			reduce_max<512><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_max, d_partial_result, BLOCKS_TOTAL);
			checkKernelErr();
			break;
		case 256:
			reduce_max<256><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_max, d_partial_result, BLOCKS_TOTAL);
			checkKernelErr();
			break;
		case 128:
			reduce_max<128><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_max, d_partial_result, BLOCKS_TOTAL);
			checkKernelErr();
			break;
		case 16:
			reduce_max<16><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_max, d_partial_result, BLOCKS_TOTAL);
			checkKernelErr();
			break;
		default:
			//TODO:#4
			//Quick'n'dirty purkka; not very practical with grid sizes +256 because the reduction evaluates only the first 1024 values
			//parallerly. Every value in indices +1024 is evaluated sequentially. With 256^3 grid there will be a workload of 8 elements 
			//per thread (d_partial_result would be of size 8192 or something like that) that is okay, but will probably affect on
			//performance with 512^3 grids or more when one thread has too much work to do.
			if (BLOCKS_TOTAL > 1024) {
				reduce_max<1024><<<1, 1024, 1024*sizeof(float)>>>(d_vec_max, d_partial_result, BLOCKS_TOTAL);		
				checkKernelErr();
			}
			else {
				//Todo support for other sizes
				printf("INCORRECT BLOCKS_TOTAL (= %d) IN collectiveops.cu!\n", BLOCKS_TOTAL);
				exit(EXIT_FAILURE);
			}
	}
}

void min_vec_cuda(float* d_vec_min, float* d_partial_result, float* d_vec_x, float* d_vec_y, float* d_vec_z)
{
	//-------------------------------------------------
	static dim3 threadsPerBlock, blocksPerGrid;
	threadsPerBlock.x = COL_THREADS_X;
	threadsPerBlock.y = COL_THREADS_Y;
	threadsPerBlock.z = 1; // 2D blockdims only 
	static const int SMEM_PER_BLOCK = threadsPerBlock.x * threadsPerBlock.y *threadsPerBlock.z * sizeof(float);

	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)COL_THREADS_X);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)COL_THREADS_Y);
	blocksPerGrid.z = ceil((float) COMP_DOMAIN_SIZE_Z / (float)COL_ELEMS_PER_THREAD);
	static const int BLOCKS_TOTAL = blocksPerGrid.x * blocksPerGrid.y * blocksPerGrid.z;
	//------------------------------------------------
	
	switch (threadsPerBlock.x*threadsPerBlock.y*threadsPerBlock.z)
	{
		case 512:
			min_vec<512><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z, COL_ELEMS_PER_THREAD); break;
		case 256:
			min_vec<256><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z, COL_ELEMS_PER_THREAD); break;
		case 128:
			min_vec<128><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z, COL_ELEMS_PER_THREAD); break;
		default:
			printf("INCORRECT THREAD SIZE!\n");
			exit(EXIT_FAILURE);
	}
	switch (BLOCKS_TOTAL){
		case 1024:
			reduce_min<1024><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_min, d_partial_result, BLOCKS_TOTAL);
			break;
		case 512:
			reduce_min<512><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_min, d_partial_result, BLOCKS_TOTAL);
			break;
		case 256:
			reduce_min<256><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_min, d_partial_result, BLOCKS_TOTAL);
			break;
		case 128:
			reduce_min<128><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_min, d_partial_result, BLOCKS_TOTAL);
			break;
		case 16:
			reduce_min<16><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_min, d_partial_result, BLOCKS_TOTAL);
			break;
		default:
			//TODO:#4
			if (BLOCKS_TOTAL > 1024) {
				reduce_min<1024><<<1, 1024, 1024*sizeof(float)>>>(d_vec_min, d_partial_result, BLOCKS_TOTAL);		
			}
			else {
				//Todo support for other sizes
				printf("INCORRECT BLOCKS_TOTAL (= %d) IN collectiveops.cu!\n", BLOCKS_TOTAL);
				exit(EXIT_FAILURE);
			}
	}
}

//Calculate the rms of the velocity vector 
void vec_rms_cuda(float* d_vec_rms, float* d_partial_result, float* d_vec_x, float* d_vec_y, float* d_vec_z, bool root=true)
{
	//-------------------------------------------------
	static dim3 threadsPerBlock, blocksPerGrid;
	threadsPerBlock.x = COL_THREADS_X;
	threadsPerBlock.y = COL_THREADS_Y;
	threadsPerBlock.z = 1; // 2D blockdims only 
	static const int SMEM_PER_BLOCK = threadsPerBlock.x * threadsPerBlock.y *threadsPerBlock.z * sizeof(float);

	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)COL_THREADS_X);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)COL_THREADS_Y);
	blocksPerGrid.z = ceil((float) COMP_DOMAIN_SIZE_Z / (float)COL_ELEMS_PER_THREAD);
	static const int BLOCKS_TOTAL = blocksPerGrid.x * blocksPerGrid.y * blocksPerGrid.z;
	//------------------------------------------------
	
	switch (threadsPerBlock.x*threadsPerBlock.y*threadsPerBlock.z)
	{
		case 512:
			vec2_sum<512><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z, COL_ELEMS_PER_THREAD); break;
		case 256:
			vec2_sum<256><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z, COL_ELEMS_PER_THREAD); break;
		case 128:
			vec2_sum<128><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_vec_x, d_vec_y, d_vec_z, COL_ELEMS_PER_THREAD); break;
		default:
			printf("INCORRECT THREAD SIZE!\n");
			exit(EXIT_FAILURE);
	}
	switch (BLOCKS_TOTAL){
		case 1024:
			reduce_rms<1024><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 512:
			reduce_rms<512><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 256:
			reduce_rms<256><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 128:
			reduce_rms<128><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 16:
			reduce_rms<16><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_vec_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		default:
			//TODO:#4
			if (BLOCKS_TOTAL > 1024) {
				reduce_rms<1024><<<1, 1024, 1024*sizeof(float)>>>(d_vec_rms, d_partial_result, BLOCKS_TOTAL, root);		
			}
			else {
				//Todo support for other sizes
				printf("INCORRECT BLOCKS_TOTAL (= %d) IN collectiveops.cu!\n", BLOCKS_TOTAL);
				exit(EXIT_FAILURE);
			}
	}
}

void max_scal_cuda(float* d_scal_max, float* d_partial_result, float* d_scal)
{
	//-------------------------------------------------
	static dim3 threadsPerBlock, blocksPerGrid;
	threadsPerBlock.x = COL_THREADS_X;
	threadsPerBlock.y = COL_THREADS_Y;
	threadsPerBlock.z = 1; // 2D blockdims only 
	static const int SMEM_PER_BLOCK = threadsPerBlock.x * threadsPerBlock.y *threadsPerBlock.z * sizeof(float);

	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)COL_THREADS_X);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)COL_THREADS_Y);
	blocksPerGrid.z = ceil((float) COMP_DOMAIN_SIZE_Z / (float)COL_ELEMS_PER_THREAD);
	static const int BLOCKS_TOTAL = blocksPerGrid.x * blocksPerGrid.y * blocksPerGrid.z;
	//------------------------------------------------
	
	switch (threadsPerBlock.x*threadsPerBlock.y*threadsPerBlock.z)
	{
		case 512:
			max_scal<512><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD); break;
		case 256:
			max_scal<256><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD); break;
		case 128:
			max_scal<128><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD); break;
		default:
			printf("INCORRECT THREAD SIZE!\n");
			exit(EXIT_FAILURE);
	}
	switch (BLOCKS_TOTAL){
		case 1024:
			reduce_max<1024><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_max, d_partial_result, BLOCKS_TOTAL);
			break;
		case 512:
			reduce_max<512><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_max, d_partial_result, BLOCKS_TOTAL);
			break;
		case 256:
			reduce_max<256><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_max, d_partial_result, BLOCKS_TOTAL);
			break;
		case 128:
			reduce_max<128><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_max, d_partial_result, BLOCKS_TOTAL);
			break;
		case 16:
			reduce_max<16><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_max, d_partial_result, BLOCKS_TOTAL);
			break;
		default:
			//TODO:#4
			//Quick'n'dirty purkka; not very practical with grid sizes +256 because the reduction evaluates only the first 1024 values
			//parallerly. Every value in indices +1024 is evaluated sequentially. With 256^3 grid there will be a workload of 8 elements 
			//per thread (d_partial_result would be of size 8192 or something like that) that is okay, but will probably affect on
			//performance with 512^3 grids or more when one thread has too much work to do.
			if (BLOCKS_TOTAL > 1024) {
				reduce_max<1024><<<1, 1024, 1024*sizeof(float)>>>(d_scal_max, d_partial_result, BLOCKS_TOTAL);		
			}
			else {
				//Todo support for other sizes
				printf("INCORRECT BLOCKS_TOTAL (= %d) IN collectiveops.cu!\n", BLOCKS_TOTAL);
				exit(EXIT_FAILURE);
			}
	}
}

void min_scal_cuda(float* d_scal_min, float* d_partial_result, float* d_scal)
{
	//-------------------------------------------------
	static dim3 threadsPerBlock, blocksPerGrid;
	threadsPerBlock.x = COL_THREADS_X;
	threadsPerBlock.y = COL_THREADS_Y;
	threadsPerBlock.z = 1; // 2D blockdims only 
	static const int SMEM_PER_BLOCK = threadsPerBlock.x * threadsPerBlock.y *threadsPerBlock.z * sizeof(float);

	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)COL_THREADS_X);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)COL_THREADS_Y);
	blocksPerGrid.z = ceil((float) COMP_DOMAIN_SIZE_Z / (float)COL_ELEMS_PER_THREAD);
	static const int BLOCKS_TOTAL = blocksPerGrid.x * blocksPerGrid.y * blocksPerGrid.z;
	//------------------------------------------------
	
	switch (threadsPerBlock.x*threadsPerBlock.y*threadsPerBlock.z)
	{
		case 512:
			min_scal<512><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD); break;
		case 256:
			min_scal<256><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD); break;
		case 128:
			min_scal<128><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD); break;
		default:
			printf("INCORRECT THREAD SIZE!\n");
			exit(EXIT_FAILURE);
	}
	switch (BLOCKS_TOTAL){
		case 1024:
			reduce_min<1024><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_min, d_partial_result, BLOCKS_TOTAL);
			break;
		case 512:
			reduce_min<512><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_min, d_partial_result, BLOCKS_TOTAL);
			break;
		case 256:
			reduce_min<256><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_min, d_partial_result, BLOCKS_TOTAL);
			break;
		case 128:
			reduce_min<128><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_min, d_partial_result, BLOCKS_TOTAL);
			break;
		case 16:
			reduce_min<16><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_min, d_partial_result, BLOCKS_TOTAL);
			break;
		default:
			//TODO:#4
			if (BLOCKS_TOTAL > 1024) {
				reduce_min<1024><<<1, 1024, 1024*sizeof(float)>>>(d_scal_min, d_partial_result, BLOCKS_TOTAL);		
			}
			else {
				//Todo support for other sizes
				printf("INCORRECT BLOCKS_TOTAL (= %d) IN collectiveops.cu!\n", BLOCKS_TOTAL);
				exit(EXIT_FAILURE);
			}
	}
}

//Calculate the rms of a scalar 
void scal_rms_cuda(float* d_scal_rms, float* d_partial_result, float* d_scal, bool sqr=true, bool root=true)
{
	//-------------------------------------------------
	static dim3 threadsPerBlock, blocksPerGrid;
	threadsPerBlock.x = COL_THREADS_X;
	threadsPerBlock.y = COL_THREADS_Y;
	threadsPerBlock.z = 1; // 2D blockdims only 
	static const int SMEM_PER_BLOCK = threadsPerBlock.x * threadsPerBlock.y *threadsPerBlock.z * sizeof(float);

	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)COL_THREADS_X);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)COL_THREADS_Y);
	blocksPerGrid.z = ceil((float) COMP_DOMAIN_SIZE_Z / (float)COL_ELEMS_PER_THREAD);
	static const int BLOCKS_TOTAL = blocksPerGrid.x * blocksPerGrid.y * blocksPerGrid.z;
	//------------------------------------------------
	
	switch (threadsPerBlock.x*threadsPerBlock.y*threadsPerBlock.z)
	{
		case 512:
			scal2_sum<512><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD, sqr); break;
		case 256:
			scal2_sum<256><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD, sqr); break;
		case 128:
			scal2_sum<128><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD, sqr); break;
		default:
			printf("INCORRECT THREAD SIZE!\n");
			exit(EXIT_FAILURE);
	}
	switch (BLOCKS_TOTAL){
		case 1024:
			reduce_rms<1024><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 512:
			reduce_rms<512><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 256:
			reduce_rms<256><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 128:
			reduce_rms<128><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 16:
			reduce_rms<16><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		default:
			//TODO:#4
			if (BLOCKS_TOTAL > 1024) {
				reduce_rms<1024><<<1, 1024, 1024*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);		
			}
			else {
				//Todo support for other sizes
				printf("INCORRECT BLOCKS_TOTAL (= %d) IN collectiveops.cu!\n", BLOCKS_TOTAL);
				exit(EXIT_FAILURE);
			}
	}
}

//Calculate the rms of a exp(scalar) 
void scal_exp_rms_cuda(float* d_scal_rms, float* d_partial_result, float* d_scal, bool sqr=true, bool root=true)
{
	//-------------------------------------------------
	static dim3 threadsPerBlock, blocksPerGrid;
	threadsPerBlock.x = COL_THREADS_X;
	threadsPerBlock.y = COL_THREADS_Y;
	threadsPerBlock.z = 1; // 2D blockdims only 
	static const int SMEM_PER_BLOCK = threadsPerBlock.x * threadsPerBlock.y *threadsPerBlock.z * sizeof(float);

	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)COL_THREADS_X);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)COL_THREADS_Y);
	blocksPerGrid.z = ceil((float) COMP_DOMAIN_SIZE_Z / (float)COL_ELEMS_PER_THREAD);
	static const int BLOCKS_TOTAL = blocksPerGrid.x * blocksPerGrid.y * blocksPerGrid.z;
	//------------------------------------------------
	switch (threadsPerBlock.x*threadsPerBlock.y*threadsPerBlock.z)
	{
		case 512:
			scal2_exp_sum<512><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD, sqr); break;
		case 256:
			scal2_exp_sum<256><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD, sqr); break;
		case 128:
			scal2_exp_sum<128><<<blocksPerGrid, threadsPerBlock, SMEM_PER_BLOCK>>>(d_partial_result, d_scal, COL_ELEMS_PER_THREAD, sqr); break;
		default:
			printf("INCORRECT THREAD SIZE!\n");
			exit(EXIT_FAILURE);
	}
	switch (BLOCKS_TOTAL){
		case 1024:
			reduce_rms<1024><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 512:
			reduce_rms<512><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 256:
			reduce_rms<256><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 128:
			reduce_rms<128><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		case 16:
			reduce_rms<16><<<1, BLOCKS_TOTAL, BLOCKS_TOTAL*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);
			break;
		default:
			//TODO:#4
			if (BLOCKS_TOTAL > 1024) {
				reduce_rms<1024><<<1, 1024, 1024*sizeof(float)>>>(d_scal_rms, d_partial_result, BLOCKS_TOTAL, root);		
			}
			else {
				//Todo support for other sizes
				printf("INCORRECT BLOCKS_TOTAL (= %d) IN collectiveops.cu!\n", BLOCKS_TOTAL);
				exit(EXIT_FAILURE);
			}
	}
}
