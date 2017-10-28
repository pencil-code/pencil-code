/* Date:   01-03-2017
   Author: Omer Anjum
   Description:
   Copying internal halos from GPU to host
Comments: 
Date: March 10, 2017
Omer Anjum
Very first version of code written. 
*/
#include <stdio.h>
#include "copyhalos.cuh"
/****************************************************************************************/
__global__ void copy_internal_rows(float* d_halo, float* d_grid, int nx, int ny, int nz, int halo_depth, dim3 blocksPerGrid)
{ 
	//int halo_size = (halo_depth*nx*2 + halo_depth*(ny-halo_depth*2)*2)*(nz-halo_depth*2) + nx*ny*(halo_depth*2);
	const int halo_idx_x = threadIdx.x + blockIdx.x*blockDim.x; 
	const int halo_idx_y = threadIdx.y + blockIdx.y*blockDim.y;
	const int halo_idx_z = threadIdx.z + blockIdx.z*blockDim.z;
	
	int halo_idx = 	(halo_idx_x) + (halo_idx_y)*(nx-2*halo_depth) + (halo_idx_z)*((nx-2*halo_depth)*(halo_depth*2)+(ny-(halo_depth*4))*(halo_depth*2));//last term 128*6+128*6

	int d_grid_idx = (halo_idx_x+halo_depth) + (halo_idx_y+halo_depth)*nx + (halo_idx_z+halo_depth)*nx*ny;
	if(halo_idx_x < nx-2*halo_depth){
		d_halo[halo_idx] = d_grid[d_grid_idx];
		d_halo[halo_idx+((nx-2*halo_depth)*halo_depth+(ny-(halo_depth*4))*(halo_depth*2))] = d_grid[d_grid_idx+(ny-3*halo_depth)*nx];
	}
	
}
/****************************************************************************************/
__global__ void copy_internal_cols(float* d_halo, float* d_grid, int nx, int ny, int nz, int halo_depth, dim3 blocksPerGrid)
{ 
	//int halo_size = (halo_depth*nx*2 + halo_depth*(ny-halo_depth*2)*2)*(nz-halo_depth*2) + nx*ny*(halo_depth*2);
	const int halo_idx_x = threadIdx.x + blockIdx.x*blockDim.x; 
	const int halo_idx_y = threadIdx.y + blockIdx.y*blockDim.y;
	const int halo_idx_z = threadIdx.z + blockIdx.z*blockDim.z;

	int halo_idx = halo_depth*(nx-2*halo_depth) +	(halo_idx_x) + (halo_idx_y)*2*halo_depth + (halo_idx_z)*((nx-2*halo_depth)*(halo_depth*2)+(ny-(halo_depth*4))*(halo_depth*2));//last term 134*6+128*6, first term taking threads to where columns data starts
	
	int d_grid_idx = (halo_idx_x+halo_depth) + (halo_idx_y+2*halo_depth)*nx + (halo_idx_z+halo_depth)*nx*ny;
	if(halo_idx_y < ny-4*halo_depth){ 
		d_halo[halo_idx] = d_grid[d_grid_idx];
		d_halo[halo_idx+halo_depth] = d_grid[d_grid_idx+(nx-3*halo_depth)];//---|idx|------|nx|---|nx+idx|
	}
}
/****************************************************************************************/
__global__ void copy_internal_frtbk(float* d_halo, float* d_grid, int nx, int ny, int nz, int halo_depth, dim3 blocksPerGrid)
{ 
	//int halo_size = (halo_depth*nx*2 + halo_depth*(ny-halo_depth*2)*2)*(nz-halo_depth*2) + nx*ny*(halo_depth*2);
	const int halo_idx_x = threadIdx.x + blockIdx.x*blockDim.x; 
	const int halo_idx_y = threadIdx.y + blockIdx.y*blockDim.y;
	const int halo_idx_z = threadIdx.z + blockIdx.z*blockDim.z;

	int halo_idx = (halo_depth*(nx-2*halo_depth)*2 +(ny-(halo_depth*4))*(halo_depth*2))*(nz-2*halo_depth) + (halo_idx_x) + (halo_idx_y)*(nx-2*halo_depth) + (halo_idx_z)*(nx-2*halo_depth)*(ny-2*halo_depth);//last term 134*6+128*6, first term taking threads to where columns data starts

	int d_grid_idx = (halo_idx_x+halo_depth) + (halo_idx_y+halo_depth)*nx + (halo_idx_z)*nx*ny;
	if(halo_idx_x < nx - 2*halo_depth && halo_idx_y < ny - 2*halo_depth && halo_idx_z < nz){
		d_halo[halo_idx] = d_grid[d_grid_idx];
		d_halo[halo_idx+(nx-2*halo_depth)*(ny-2*halo_depth)*halo_depth] = d_grid[d_grid_idx+nx*ny*(nz-halo_depth)];
	}
	/*__syncthreads();

	if(threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0 && blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0) {
		printf("Writing thread (%d,%d,%d) at block (%d,%d,%d) \n",threadIdx.x, threadIdx.y, threadIdx.z, 	 				blockIdx.x,blockIdx.y,blockIdx.z );
		printf("\n printing halo\n");
		for (int k=0; k < halo_size; k++) {
			printf("%d, ",d_halo[k]);
		}	

	}*/
}
/****************************************************************************************/
void fillhalosinhost(float* d_halo, float* d_grid, int nx, int ny, int nz, int halo_depth)
{
	//int ELEMS_PER_THREAD_in_z = nz-(2*halo_depth);
	//TODO: Adapt for shearing-periodic case
	static dim3 blocksPerGrid, threadsPerBlock;

	//Create streams for executing the boundary copy 
	//kernels concurrently.
	/*static cudaStream_t per_row_stream = NULL; 
	if (per_row_stream == NULL)
		cudaStreamCreate(&per_row_stream);
	static cudaStream_t per_col_stream = NULL; 
	if (per_col_stream == NULL)
		cudaStreamCreate(&per_col_stream);
	static cudaStream_t per_frtbk_stream = NULL; 
	if (per_frtbk_stream == NULL)
		cudaStreamCreate(&per_frtbk_stream);*/

	//Copy the top and bottom halos around the compute grid
	threadsPerBlock.x = 6;// increase to 32
	threadsPerBlock.y = halo_depth; // do not change
	threadsPerBlock.z = 1; // do not change
	blocksPerGrid.x = (int)ceil((double)nx-2*halo_depth / (double)threadsPerBlock.x);
	printf("\n %d, %d,",blocksPerGrid.x, threadsPerBlock.y);
	blocksPerGrid.y = 1;
	blocksPerGrid.z = nz-(2*halo_depth);
	//printf(" %d block in z= %d",threadsPerBlock.z, blocksPerGrid.z);
	//printf("\n----------------------\ngoing inside the kernel to copy rows\n-----------------------------\n");
	copy_internal_rows<<<blocksPerGrid, threadsPerBlock>>>(d_halo, d_grid, nx, ny, nz, halo_depth, blocksPerGrid);
	checkKernelErr();
	cudaDeviceSynchronize();
	

	//Copy the top and bottom halos around the compute grid
	threadsPerBlock.x = halo_depth; // do not change
	threadsPerBlock.y = 2; // increase to 32
	threadsPerBlock.z = 1; // do not change
	//printf("\n %d \n",threadsPerBlock.y);
	blocksPerGrid.x = 1;
	blocksPerGrid.y = (int)ceil((double)(ny-2*halo_depth) / (double)threadsPerBlock.y);
	//printf("%d blocksPerGrid.y \n", blocksPerGrid.y);
	blocksPerGrid.z = nz-(2*halo_depth);
	//printf(" %d block in z= %d",threadsPerBlock.z, blocksPerGrid.z);
	//printf("\n----------------------\ngoing inside the kernel to copy cols\n-----------------------------\n");
	copy_internal_cols<<<blocksPerGrid, threadsPerBlock>>>(d_halo, d_grid, nx, ny, nz, halo_depth, blocksPerGrid);
	checkKernelErr();
	cudaDeviceSynchronize();
	

	//Copy the front and back halos around the compute grid
	threadsPerBlock.x = 4;// increase to 32
	threadsPerBlock.y = 6;// increase to 32
	threadsPerBlock.z = 1; // do not change
	//printf("\n %d \n",threadsPerBlock.y);
	blocksPerGrid.x = (int)ceil((double)(nx-2*halo_depth) / (double)threadsPerBlock.x);
	blocksPerGrid.y = (int)ceil((double)(ny-2*halo_depth) / (double)threadsPerBlock.y);
	//printf("%d blocksPerGrid.y \n", blocksPerGrid.y);
	blocksPerGrid.z = halo_depth;
	//printf(" %d block in z= %d",threadsPerBlock.z, blocksPerGrid.z);
	//printf("\n----------------------\ngoing inside the kernel to copy frtbk\n-----------------------------\n");
	copy_internal_frtbk<<<blocksPerGrid, threadsPerBlock>>>(d_halo, d_grid, nx, ny, nz, halo_depth, blocksPerGrid);
	checkKernelErr();
	cudaDeviceSynchronize();
	
	//checkErr(cudaStreamDestroy(per_row_stream));
	//checkErr(cudaStreamDestroy(per_col_stream));
	//checkErr(cudaStreamDestroy(per_frtbk_stream));
	//printf("\n came back \n");

return;
	
}
/****************************************************************************************/


