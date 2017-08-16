

#include <stdio.h>

#include "slice.cuh"
#define EXTERN extern
#include "dconsts.cuh"
#include "defines.h"

/*
//Slices the only the computational domain (not used)
template <char slice_axis>
__global__ void animation_slice(float* d_slice_lnrho, float* d_slice_uu, 
				float* d_slice_uu_x, float* d_slice_uu_y, float* d_slice_uu_z, 
				float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z)
{
	int slice_idx_x = threadIdx.x + (blockIdx.x*blockDim.x);
	int slice_idx_y = threadIdx.y + (blockIdx.y*blockDim.y);
	
	int width = -1;
	int grid_idx = -1;
	switch (slice_axis)
	{
		case 'y':
			grid_idx = (slice_idx_x+d_CX_BOT) + (d_NY/2)*d_NX + (slice_idx_y+d_CZ_BOT)*d_NX*d_NY;
			width = d_COMP_DOMAIN_SIZE_X; 
			break;
		case 'x':
			grid_idx = (d_NX/2) + (slice_idx_y+d_CY_BOT)*d_NX + (slice_idx_x+d_CZ_BOT)*d_NX*d_NY;
			width = d_COMP_DOMAIN_SIZE_Z;
			break;
		case 'z':
			grid_idx = (slice_idx_x+d_CX_BOT) + (slice_idx_y+d_CY_BOT)*d_NX + (d_NZ/2)*d_NX*d_NY;
			width = d_COMP_DOMAIN_SIZE_X;
			break;
	}

	//Get slice index
	int slice_idx = slice_idx_x + slice_idx_y*width;
	
	//Save data to slice (hopefully raises an error if trying to access grid_idx = -1)
	d_slice_lnrho[slice_idx] = d_lnrho[grid_idx];

	float uu_x, uu_y, uu_z; //Fetch uu once from global memory to (hopefully) registers
	uu_x = d_uu_x[grid_idx]; uu_y = d_uu_y[grid_idx]; uu_z = d_uu_z[grid_idx];
	
	d_slice_uu[slice_idx] = uu_x*uu_x + uu_y*uu_y + uu_z*uu_z;
	d_slice_uu_x[slice_idx] = uu_x;
	d_slice_uu_y[slice_idx] = uu_y;
	d_slice_uu_z[slice_idx] = uu_z;
}*/


//Puts a debug slice in to d_slice_lnrho etc that contains
//the boundary zones (no paddings)
template <char slice_axis>
__global__ void animation_slice_debug(	float* d_slice_lnrho, float* d_slice_uu, 
					float* d_slice_uu_x, float* d_slice_uu_y, float* d_slice_uu_z, 
					float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z)
{
	int slice_idx_x = threadIdx.x + (blockIdx.x*blockDim.x);
	int slice_idx_y = threadIdx.y + (blockIdx.y*blockDim.y);
	
	int width = -1;
	int grid_idx = -1;
	switch (slice_axis)
	{
		case 'y':
			grid_idx = (slice_idx_x+d_CX_BOT-d_BOUND_SIZE) + (d_NY/2)*d_NX + (slice_idx_y)*d_NX*d_NY;
			width = d_COMP_DOMAIN_SIZE_X + 2*BOUND_SIZE; 
			//Check bounds, we might not get thread dims that are divisible by d_NX/NY/NZ
			if (slice_idx_x >= d_COMP_DOMAIN_SIZE_X + 2*BOUND_SIZE || slice_idx_y >= d_COMP_DOMAIN_SIZE_Z + 2*BOUND_SIZE) return;
			break;
		case 'x':
			grid_idx = (d_NX/2) + (slice_idx_y)*d_NX + (slice_idx_x)*d_NX*d_NY;
			width = d_COMP_DOMAIN_SIZE_Z + 2*BOUND_SIZE; 
			if (slice_idx_x >= d_COMP_DOMAIN_SIZE_Z + 2*BOUND_SIZE || slice_idx_y >= d_COMP_DOMAIN_SIZE_Y + 2*BOUND_SIZE) return;
			break;
		case 'z':
			grid_idx = (slice_idx_x+d_CX_BOT-d_BOUND_SIZE) + (slice_idx_y)*d_NX + (d_NZ/2)*d_NX*d_NY;
			width = d_COMP_DOMAIN_SIZE_X + 2*BOUND_SIZE; 
			if (slice_idx_x >= d_COMP_DOMAIN_SIZE_X + 2*BOUND_SIZE || slice_idx_y >= d_COMP_DOMAIN_SIZE_Y + 2*BOUND_SIZE) return;
			break;
	}

	//Get slice index
	int slice_idx = slice_idx_x + slice_idx_y*width;
	
	//Save data to slice (hopefully raises an error if trying to access grid_idx = -1)
	d_slice_lnrho[slice_idx] = d_lnrho[grid_idx];

	float uu_x, uu_y, uu_z; //Fetch uu once from global memory to (hopefully) registers
	uu_x = d_uu_x[grid_idx]; uu_y = d_uu_y[grid_idx]; uu_z = d_uu_z[grid_idx];
	
	d_slice_uu[slice_idx] = uu_x*uu_x + uu_y*uu_y + uu_z*uu_z;
	d_slice_uu_x[slice_idx] = uu_x;
	d_slice_uu_y[slice_idx] = uu_y;
	d_slice_uu_z[slice_idx] = uu_z;

}


//Slices the assigned axis to d_slice_lnrho etc in device memory
void get_slice_cuda( 	char slice_axis, float* d_slice_lnrho, float* d_slice_uu, 
		 	float* d_slice_uu_x, float* d_slice_uu_y, float* d_slice_uu_z, 
		 	float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z)
{
	int x_size, y_size;
	x_size = y_size = 2*BOUND_SIZE;
	static dim3 threadsPerBlock, blocksPerGrid;
	threadsPerBlock.x = 32;
	threadsPerBlock.y = 32;
	threadsPerBlock.z = 1;

	switch (slice_axis) {
		case 'x':
			x_size += COMP_DOMAIN_SIZE_Z;
			y_size += COMP_DOMAIN_SIZE_Y;
			blocksPerGrid.x = ceil((float) x_size / (float)threadsPerBlock.x);
			blocksPerGrid.y = ceil((float) y_size / (float)threadsPerBlock.y);
			blocksPerGrid.z = 1;
	
			//Slice
			animation_slice_debug<'x'><<<blocksPerGrid, threadsPerBlock>>>(d_slice_lnrho, d_slice_uu, d_slice_uu_x, d_slice_uu_y, d_slice_uu_z, d_lnrho, d_uu_x, d_uu_y, d_uu_z);
			break;
		case 'y':
			x_size += COMP_DOMAIN_SIZE_X;
			y_size += COMP_DOMAIN_SIZE_Z;
			blocksPerGrid.x = ceil((float) x_size / (float)threadsPerBlock.x);
			blocksPerGrid.y = ceil((float) y_size / (float)threadsPerBlock.y);
			blocksPerGrid.z = 1;
	
			//Slice
			animation_slice_debug<'y'><<<blocksPerGrid, threadsPerBlock>>>(d_slice_lnrho, d_slice_uu, d_slice_uu_x, d_slice_uu_y, d_slice_uu_z, d_lnrho, d_uu_x, d_uu_y, d_uu_z);
			break;
		case 'z':
			x_size += COMP_DOMAIN_SIZE_X;
			y_size += COMP_DOMAIN_SIZE_Y;
			blocksPerGrid.x = ceil((float) x_size / (float)threadsPerBlock.x);
			blocksPerGrid.y = ceil((float) y_size / (float)threadsPerBlock.y);
			blocksPerGrid.z = 1;
	
			//Slice
			animation_slice_debug<'z'><<<blocksPerGrid, threadsPerBlock>>>(d_slice_lnrho, d_slice_uu, d_slice_uu_x, d_slice_uu_y, d_slice_uu_z, d_lnrho, d_uu_x, d_uu_y, d_uu_z);
			break;
		default:
			printf("Invalid slice axis in slice.cu:save_slice_cuda()!\n");
			exit(EXIT_FAILURE);	
	}

}










