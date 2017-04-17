
#include <stdio.h>

#include "boundcond.cuh"
#include "dconstsextern.cuh"
#include "shear.cuh"

#include "diagnostics.cuh"

/*
*	This module is used to copy the cover of the computational domain to
*	the appropriate boundary zones. For reference;
*		-UP is used when moving to the positive direction along y-axis (+y)
*		-RIGHT is used when moving to the positive direction along x-axis (+x)
*		-BACK is used when moving to the positive direction along z-axis (+z)
*
*	Check astaroth-code/doc/boundcond.png for an illustrative picture (deprecated, only 
*	the sides are copied in this version using 19-point stencils)
*/

//Copies the left and right sides of the computational domain to an appropriate
//boundary zone (does not include the edges and corners of the boundary zone)
__global__ void per_x_sides(float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z)
{

        int ix, ix_bound, sid_depth;
        if (threadIdx.x < 3) { //Copy left of the computational domain to the boundary zone at the right
                ix = threadIdx.x + d_CX_BOT;
                ix_bound = threadIdx.x + d_CX_TOP;
        }
        else { //Copy right of the computational domain to the boundary zone at the left
                ix = (threadIdx.x-3) + (d_CX_TOP - 3);
                ix_bound = (d_CX_BOT-3) +(threadIdx.x-3);
        }
	sid_depth = threadIdx.x;

        int iy,iz, sid_y;
        iz = threadIdx.z + blockIdx.z*blockDim.z + d_CZ_BOT;//Don't add edges
        iy = threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT;
	sid_y = threadIdx.y;

	/*
	if (iy < d_CY_TOP && iz < d_CZ_TOP) {

		int grid_idx =  ix + iy*d_NX + iz*d_NX*d_NY;
                int bound_idx =  ix_bound + iy*d_NX + iz*d_NX*d_NY;

                d_lnrho[bound_idx] = d_lnrho[grid_idx];
                d_uu_x[bound_idx] = d_uu_x[grid_idx];
                d_uu_y[bound_idx] = d_uu_y[grid_idx];
                d_uu_z[bound_idx] = d_uu_z[grid_idx];

        }
	*/

	//Normal periodic boundary if shearing is not included
	if (d_LSHEAR == 0 || (d_DELTA_Y > -1.0e-30 && d_DELTA_Y < 1.0e-30) ) { //Meanign if d_DELTA_Y == 0, but floats are not that stable. 
		if (iy < d_CY_TOP && iz < d_CZ_TOP) {

			int grid_idx =  ix + iy*d_NX + iz*d_NX*d_NY;
			int bound_idx =  ix_bound + iy*d_NX + iz*d_NX*d_NY;

			d_lnrho[bound_idx] = d_lnrho[grid_idx];
			d_uu_x[bound_idx] = d_uu_x[grid_idx];
			d_uu_y[bound_idx] = d_uu_y[grid_idx];
			d_uu_z[bound_idx] = d_uu_z[grid_idx];
		}
	} else {
		if (iy < d_CY_TOP && iz < d_CZ_TOP) {
			int bound_idx =  ix_bound + iy*d_NX + iz*d_NX*d_NY;

			//Allocate shared memory for interpolation arrays 
			__shared__ float s_coord_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH];
			__shared__ float s_val_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH];

			//Perform the interpolation in assign the values to the boundaries
			d_lnrho[bound_idx] = interp_shear(d_lnrho, ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);
			d_uu_x[bound_idx]  = interp_shear(d_uu_x,  ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);
			d_uu_y[bound_idx]  = interp_shear(d_uu_y,  ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);
			d_uu_z[bound_idx]  = interp_shear(d_uu_z,  ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);

		}
	}

}

//Copies the top and bottom of the computational domain to an appropriate
//boundary zone (does not include the edges and corners of the boundary zone)
__global__ void per_y_sides(float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z)
{
	int iy, iy_bound;
	if (blockIdx.z < 3) { //Copy bottom of the computational domain to the boundary zone at the top
		iy = blockIdx.z + d_CY_BOT;
		iy_bound = blockIdx.z + d_CY_TOP;
	} 
	else { //Copy top of the computational domain to the boundary zone at the bottom
		iy = (blockIdx.z-3) + (d_CY_TOP - 3);
		iy_bound = (blockIdx.z-3);
	}

	int ix,iz;
	ix = threadIdx.x + blockIdx.x*blockDim.x + d_CX_BOT;
	iz = threadIdx.y + blockIdx.y*blockDim.y + d_CZ_BOT;

	if (ix < d_CX_TOP && iz < d_CZ_TOP) {

		int grid_idx =  ix + iy*d_NX + iz*d_NX*d_NY;
		int bound_idx =  ix + iy_bound*d_NX + iz*d_NX*d_NY;

		d_lnrho[bound_idx] = d_lnrho[grid_idx];
		d_uu_x[bound_idx] = d_uu_x[grid_idx];
		d_uu_y[bound_idx] = d_uu_y[grid_idx];
		d_uu_z[bound_idx] = d_uu_z[grid_idx];

	}
}

//Copies the front and back of the computational domain to an appropriate
//boundary zone (does not include the edges and corners of the boundary zone)
__global__ void per_z_sides(float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z)
{
	int iz, iz_bound;
	if (blockIdx.z < 3) { //Copy front of the computational domain to the boundary zone at the back
		iz = blockIdx.z + d_CZ_BOT;
		iz_bound = blockIdx.z + d_CZ_TOP;
	} 
	else { //Copy back of the computational domain to the boundary zone at the front
		iz = (blockIdx.z-3) + d_CZ_TOP - 3;
		iz_bound = (blockIdx.z-3);
	}

	int ix,iy;
	ix = threadIdx.x + blockIdx.x*blockDim.x + d_CX_BOT;
	iy = threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT;

	if (ix < d_CX_TOP && iy < d_CY_TOP) {

		int grid_idx =  ix + iy*d_NX + iz*d_NX*d_NY;
		int bound_idx =  ix + iy*d_NX + iz_bound*d_NX*d_NY;

		d_lnrho[bound_idx] = d_lnrho[grid_idx];
		d_uu_x[bound_idx] = d_uu_x[grid_idx];
		d_uu_y[bound_idx] = d_uu_y[grid_idx];
		d_uu_z[bound_idx] = d_uu_z[grid_idx];

	}
}


/*
__global__ void per_x_sides_shear(float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z)
{

	int ix, ix_bound, sid_depth;
	if (threadIdx.x < 3) { //Do left
		ix = threadIdx.x + d_CX_BOT;
		ix_bound = threadIdx.x + d_CX_TOP;
	} 
	else { //Do right
		ix = (threadIdx.x-3) + (d_CX_TOP - 3);
		ix_bound = (d_CX_BOT-3) +(threadIdx.x-3);
	}
	sid_depth = threadIdx.x;

	int iy,iz, sid_y; 
	iz = threadIdx.z + blockIdx.z*blockDim.z + d_CZ_BOT;//Don't add edges
	iy = threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT;
	sid_y = threadIdx.y;

	if (d_DELTA_Y > -1.0e-30 && d_DELTA_Y < 1.0e-30 ) { //Meanign if d_DELTA_Y == 0, but floats are not that stable. 
		if (iy < d_CY_TOP && iz < d_CZ_TOP) {

			int grid_idx =  ix + iy*d_NX + iz*d_NX*d_NY;
			int bound_idx =  ix_bound + iy*d_NX + iz*d_NX*d_NY;

			d_lnrho[bound_idx] = d_lnrho[grid_idx];
			d_uu_x[bound_idx] = d_uu_x[grid_idx];
			d_uu_y[bound_idx] = d_uu_y[grid_idx];
			d_uu_z[bound_idx] = d_uu_z[grid_idx];
		}
	} else {
		if (iy < d_CY_TOP && iz < d_CZ_TOP) {
			//TODO interpolation here
			int bound_idx =  ix_bound + iy*d_NX + iz*d_NX*d_NY;

			//Allocate shared memory for interpolation arrays 
			__shared__ float s_coord_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH];
			__shared__ float s_val_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH];

			//Perform the interpolation in assign the values to the boundaries
			d_lnrho[bound_idx] = interp_shear(d_lnrho, ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);
			d_uu_x[bound_idx]  = interp_shear(d_uu_x,  ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);
			d_uu_y[bound_idx]  = interp_shear(d_uu_y,  ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);
			d_uu_z[bound_idx]  = interp_shear(d_uu_z,  ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);

		}
	}
}
*/


void boundcond_cuda(float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z)
{
	//TODO: Adapt for shearing-periodic case

	static dim3 blocksPerGrid, threadsPerBlock;

	//Create streams for executing the boundary copy 
	//kernels concurrently.
	static cudaStream_t per_x_stream = NULL; 
	if (per_x_stream == NULL)
		cudaStreamCreate(&per_x_stream);
	static cudaStream_t per_y_stream = NULL; 
	if (per_y_stream == NULL)
		cudaStreamCreate(&per_y_stream);
	static cudaStream_t per_z_stream = NULL; 
	if (per_z_stream == NULL)
		cudaStreamCreate(&per_z_stream);

	//Copy the left and right sides to the appropriate boundary zone
	threadsPerBlock.x = 6;
	threadsPerBlock.y = 4;
	threadsPerBlock.z = 1;
	blocksPerGrid.x = 1;
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)threadsPerBlock.y);
	blocksPerGrid.z = ceil((float) COMP_DOMAIN_SIZE_Z / (float)threadsPerBlock.z);
	per_x_sides<<<blocksPerGrid, threadsPerBlock, 0, per_x_stream>>>(d_lnrho, d_uu_x, d_uu_y, d_uu_z);
	//checkKernelErr();


	//Copy the top and bottom sides to the appropriate boundary zone
	threadsPerBlock.x = 32;
	threadsPerBlock.y = 32;
	threadsPerBlock.z = 1;
	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)threadsPerBlock.x);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Z / (float)threadsPerBlock.y);
	blocksPerGrid.z = 6;
	per_y_sides<<<blocksPerGrid, threadsPerBlock, 0, per_y_stream>>>(d_lnrho, d_uu_x, d_uu_y, d_uu_z);
	//checkKernelErr();


	//Copy the front and back sides to the appropriate boundary zone
	threadsPerBlock.x = 32;
	threadsPerBlock.y = 32;
	threadsPerBlock.z = 1;
	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)threadsPerBlock.x);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)threadsPerBlock.y);
	blocksPerGrid.z = 6;
	per_z_sides<<<blocksPerGrid, threadsPerBlock, 0, per_z_stream>>>(d_lnrho, d_uu_x, d_uu_y, d_uu_z);
	//checkKernelErr();
	
}



//Copies the left and right sides of the computational domain to an appropriate
//boundary zone (does not include the edges and corners of the boundary zone)
__global__ void periodic_x_sides_scal(float* d_scal)
{

        int ix, ix_bound;
        if (threadIdx.x < 3) { //Copy left of the computational domain to the boundary zone at the right
                ix = threadIdx.x + d_CX_BOT;
                ix_bound = threadIdx.x + d_CX_TOP;
        }
        else { //Copy right of the computational domain to the boundary zone at the left
                ix = (threadIdx.x-3) + (d_CX_TOP - 3);
                ix_bound = (d_CX_BOT-3) +(threadIdx.x-3);
        }

        int iy,iz;
        iz = threadIdx.z + blockIdx.z*blockDim.z + d_CZ_BOT;//Don't add edges
        iy = threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT;

	if (iy < d_CY_TOP && iz < d_CZ_TOP) {

		int grid_idx =  ix + iy*d_NX + iz*d_NX*d_NY;
		int bound_idx =  ix_bound + iy*d_NX + iz*d_NX*d_NY;

		d_scal[bound_idx] = d_scal[grid_idx];
	}
}

//Copies the top and bottom of the computational domain to an appropriate
//boundary zone (does not include the edges and corners of the boundary zone)
__global__ void periodic_y_sides_scal(float* d_scal)
{
	int iy, iy_bound;
	if (blockIdx.z < 3) { //Copy bottom of the computational domain to the boundary zone at the top
		iy = blockIdx.z + d_CY_BOT;
		iy_bound = blockIdx.z + d_CY_TOP;
	} 
	else { //Copy top of the computational domain to the boundary zone at the bottom
		iy = (blockIdx.z-3) + (d_CY_TOP - 3);
		iy_bound = (blockIdx.z-3);
	}

	int ix,iz;
	ix = threadIdx.x + blockIdx.x*blockDim.x + d_CX_BOT;
	iz = threadIdx.y + blockIdx.y*blockDim.y + d_CZ_BOT;

	if (ix < d_CX_TOP && iz < d_CZ_TOP) {

		int grid_idx =  ix + iy*d_NX + iz*d_NX*d_NY;
		int bound_idx =  ix + iy_bound*d_NX + iz*d_NX*d_NY;

		d_scal[bound_idx] = d_scal[grid_idx];

	}
}

//Copies the front and back of the computational domain to an appropriate
//boundary zone (does not include the edges and corners of the boundary zone)
__global__ void periodic_z_sides_scal(float* d_scal)
{
	int iz, iz_bound;
	if (blockIdx.z < 3) { //Copy front of the computational domain to the boundary zone at the back
		iz = blockIdx.z + d_CZ_BOT;
		iz_bound = blockIdx.z + d_CZ_TOP;
	} 
	else { //Copy back of the computational domain to the boundary zone at the front
		iz = (blockIdx.z-3) + d_CZ_TOP - 3;
		iz_bound = (blockIdx.z-3);
	}

	int ix,iy;
	ix = threadIdx.x + blockIdx.x*blockDim.x + d_CX_BOT;
	iy = threadIdx.y + blockIdx.y*blockDim.y + d_CY_BOT;

	if (ix < d_CX_TOP && iy < d_CY_TOP) {

		int grid_idx =  ix + iy*d_NX + iz*d_NX*d_NY;
		int bound_idx =  ix + iy*d_NX + iz_bound*d_NX*d_NY;

		d_scal[bound_idx] = d_scal[grid_idx];
	}
}

//Copies the boundaries of an array in periodic fashion
void periodic_boundcond_scal_cuda(float* d_scal) {
	static dim3 blocksPerGrid, threadsPerBlock;

	//Copy the left and right sides to the appropriate boundary zone
	threadsPerBlock.x = 6;
	threadsPerBlock.y = 4;
	threadsPerBlock.z = 1;
	blocksPerGrid.x = 1;
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)threadsPerBlock.y);
	blocksPerGrid.z = ceil((float) COMP_DOMAIN_SIZE_Z / (float)threadsPerBlock.z);
	periodic_x_sides_scal<<<blocksPerGrid, threadsPerBlock>>>(d_scal);
	//checkKernelErr();


	//Copy the top and bottom sides to the appropriate boundary zone
	threadsPerBlock.x = 32;
	threadsPerBlock.y = 32;
	threadsPerBlock.z = 1;
	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)threadsPerBlock.x);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Z / (float)threadsPerBlock.y);
	blocksPerGrid.z = 6;
	periodic_y_sides_scal<<<blocksPerGrid, threadsPerBlock>>>(d_scal);
	//checkKernelErr();


	//Copy the front and back sides to the appropriate boundary zone
	threadsPerBlock.x = 32;
	threadsPerBlock.y = 32;
	threadsPerBlock.z = 1;
	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)threadsPerBlock.x);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)threadsPerBlock.y);
	blocksPerGrid.z = 6;
	periodic_z_sides_scal<<<blocksPerGrid, threadsPerBlock>>>(d_scal);
	//checkKernelErr();	
}








