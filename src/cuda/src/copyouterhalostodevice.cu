/* Date:   01-03-2017
   Author: Omer Anjum
   Description:
   Copying Outer halos from host to GPU
Comments: 
Date: March 14, 2017
Omer Anjum
Very first version of code written. 
*/
#include <stdio.h>


__global__ void copy_rows(float* d_halo, float* d_grid, int nx, int ny, int nz, int halo_depth, dim3 blocksPerGrid)
{ 
	
	const int halo_idx_x = threadIdx.x + blockIdx.x*blockDim.x; 
	const int halo_idx_y = threadIdx.y + blockIdx.y*blockDim.y;
	const int halo_idx_z = threadIdx.z + blockIdx.z*blockDim.z;

	int halo_idx = 	(halo_idx_x) + (halo_idx_y)*nx + (halo_idx_z)*(nx*(halo_depth*2)+(ny-(halo_depth*2))*(halo_depth*2));//last term 134*6+128*6

	int d_grid_idx = (halo_idx_x) + (halo_idx_y)*nx + (halo_idx_z+halo_depth)*nx*ny;
	if(halo_idx_x < nx && halo_idx_y < ny){
		d_grid[d_grid_idx] = d_halo[halo_idx];
		d_grid[d_grid_idx+(ny-halo_depth)*nx] = d_halo[halo_idx+(nx*halo_depth+(ny-(halo_depth*2))*(halo_depth*2))];
	}

}

__global__ void copy_cols(float* d_halo, float* d_grid, int nx, int ny, int nz, int halo_depth, dim3 blocksPerGrid)
{ 
	
	const int halo_idx_x = threadIdx.x + blockIdx.x*blockDim.x; 
	const int halo_idx_y = threadIdx.y + blockIdx.y*blockDim.y;
	const int halo_idx_z = threadIdx.z + blockIdx.z*blockDim.z;

	int halo_idx = halo_depth*nx +	(halo_idx_x) + (halo_idx_y)*2*halo_depth + (halo_idx_z)*(nx*(halo_depth*2)+(ny-(halo_depth*2))*(halo_depth*2));//last term 134*6+128*6, first term taking threads to where columns data starts

	int d_grid_idx = (halo_idx_x) + (halo_idx_y+halo_depth)*nx + (halo_idx_z+halo_depth)*nx*ny;
	if(halo_idx_x < nx && halo_idx_y < ny-halo_depth){
		//printf("d_halo[%d] = %d",halo_idx, d_halo[halo_idx]);
		//printf("\n%d %d\n",d_grid_idx, halo_idx);
		//printf("\n%d %d\n",d_grid_idx+(nx-halo_depth), halo_idx+halo_depth);
		d_grid[d_grid_idx] = d_halo[halo_idx];
		d_grid[d_grid_idx+(nx-halo_depth)] = d_halo[halo_idx+halo_depth];
	}

}

__global__ void copy_frtbk(float* d_halo, float* d_grid, int nx, int ny, int nz, int halo_depth, dim3 blocksPerGrid)
{ 
	
	const int halo_idx_x = threadIdx.x + blockIdx.x*blockDim.x; 
	const int halo_idx_y = threadIdx.y + blockIdx.y*blockDim.y;
	const int halo_idx_z = threadIdx.z + blockIdx.z*blockDim.z;

	int halo_idx = (halo_depth*nx*2 +(ny-(halo_depth*2))*(halo_depth*2))*(nz-2*halo_depth) + (halo_idx_x) + (halo_idx_y)*nx + (halo_idx_z)*nx*ny;//last term 134*6+128*6, first term taking threads to where columns data starts

	int d_grid_idx = (halo_idx_x) + (halo_idx_y)*nx + (halo_idx_z)*nx*ny;
	if(halo_idx_x < nx && halo_idx_y < ny && halo_idx_z < nz){
		//printf("d_halo[%d] = %d",halo_idx, d_halo[halo_idx]);
		//printf("\n%d %d\n",d_grid_idx, halo_idx);
		//printf("\n%d %d\n",d_grid_idx+(nx-halo_depth), halo_idx+halo_depth);
		d_grid[d_grid_idx] = d_halo[halo_idx];
		d_grid[d_grid_idx+nx*ny*(nz-halo_depth)] = d_halo[halo_idx+nx*ny*halo_depth];
	}
	/*__syncthreads();
	if(threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0 && blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0) {
		printf("Writing thread (%d,%d,%d) at block (%d,%d,%d) \n",threadIdx.x, threadIdx.y, threadIdx.z, 	 				blockIdx.x,blockIdx.y,blockIdx.z );
		for (int k=0; k < nz; k++) {
			printf("\n		-------------\n");
			for (int j=0; j < ny; j++) {
				for (int i=0; i < nx; i++) {
					int idx = i + j*nx + k*nx*ny;
					printf("	%d  ",d_grid[idx]);
				}
				printf("\n");
			}
		}	

	}*/

}


void fillhalosindevice(float* d_halo, float* d_grid, int nx, int ny, int nz, int halo_depth)
{
	//int ELEMS_PER_THREAD_in_z = nz-(2*halo_depth);
	//TODO: Adapt for shearing-periodic case
	static dim3 blocksPerGrid, threadsPerBlock;

	//Create streams for executing the boundary copy 
	//kernels concurrently.
	static cudaStream_t per_row_stream = NULL; 
	if (per_row_stream == NULL)
		cudaStreamCreate(&per_row_stream);
	static cudaStream_t per_col_stream = NULL; 
	if (per_col_stream == NULL)
		cudaStreamCreate(&per_col_stream);
	static cudaStream_t per_frtbk_stream = NULL; 
	if (per_frtbk_stream == NULL)
		cudaStreamCreate(&per_frtbk_stream);

	//Copy the top and bottom halos around the compute grid
	threadsPerBlock.x = 4; // increase to 32
	threadsPerBlock.y = halo_depth; //do not  change
	threadsPerBlock.z = 1; // do not change
	blocksPerGrid.x = (int)ceil((double)nx / (double)threadsPerBlock.x);
	printf("\n %d, %d,",blocksPerGrid.x, threadsPerBlock.y);
	blocksPerGrid.y = 1;
	blocksPerGrid.z = nz-(2*halo_depth);
	printf(" %d block in z= %d",threadsPerBlock.z, blocksPerGrid.z);
	printf("\n----------------------\ngoing inside the kernel to copy rows\n-----------------------------\n");
	copy_rows<<<blocksPerGrid, threadsPerBlock, 0, per_row_stream>>>(d_halo, d_grid, nx, ny, nz, halo_depth, blocksPerGrid);
	cudaThreadSynchronize(); // needs to be commented out at all places after first verification of code

	//Copy the top and bottom halos around the compute grid
	threadsPerBlock.x = halo_depth; // do not change
	threadsPerBlock.y = 1; // increase to 32
	threadsPerBlock.z = 1; //do not change
	//printf("\n %d \n",threadsPerBlock.y);
	blocksPerGrid.x = 1;
	blocksPerGrid.y = (int)ceil((double)(ny-2*halo_depth) / (double)threadsPerBlock.y);
	//printf("%d blocksPerGrid.y \n", blocksPerGrid.y);
	blocksPerGrid.z = nz-(2*halo_depth);
	//printf(" %d block in z= %d",threadsPerBlock.z, blocksPerGrid.z);
	printf("\n----------------------\ngoing inside the kernel to copy cols\n-----------------------------\n");
	copy_cols<<<blocksPerGrid, threadsPerBlock, 0, per_col_stream>>>(d_halo, d_grid, nx, ny, nz, halo_depth, blocksPerGrid);
	cudaThreadSynchronize();
	
	//Copy the front and back halos around the compute grid
	threadsPerBlock.x = 4; // increase to 32
	threadsPerBlock.y = 6; // increase to 32
	threadsPerBlock.z = 1; // do not change
	//printf("\n %d \n",threadsPerBlock.y);
	blocksPerGrid.x = (int)ceil((double)(nx) / (double)threadsPerBlock.x);
	blocksPerGrid.y = (int)ceil((double)(ny) / (double)threadsPerBlock.y);
	//printf("%d blocksPerGrid.y \n", blocksPerGrid.y);
	blocksPerGrid.z = halo_depth;
	//printf(" %d block in z= %d",threadsPerBlock.z, blocksPerGrid.z);
	printf("\n----------------------\ngoing inside the kernel to copy frtbk\n-----------------------------\n");
	copy_frtbk<<<blocksPerGrid, threadsPerBlock, 0, per_frtbk_stream>>>(d_halo, d_grid, nx, ny, nz, halo_depth, blocksPerGrid);
	cudaThreadSynchronize();
	printf("\n came back after filling outer halos to device\n");

	return;
	
}


