//
// Shear flow effects are calculated in this module 
//

#include <cstdio>

#define EXTERN extern
#include "dconsts.cuh"
#include "shear.cuh"
#include "smem.cuh"

__device__ float get_nearest_interp_index(float y_source_coordinate, int sid_depth, int sid_y,
			     float s_coord_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH])
{
	//Get the nearest index to the interpolated coordinate
	int ind; 

	ind = 0;
	for (int i=0; i < d_INTERP_ORDER; i++) {
		if (s_coord_interp[i][sid_y][sid_depth] <= y_source_coordinate) {
			ind = i; //Update the nearest coordinate
		}
	}

	return ind;
}

__device__ float lin_interp(float y_source_coordinate, int sid_depth, int sid_y,
			     float s_coord_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH], 
                             float s_val_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH])
{
	//Perform  a simple linear interpolation
	int ind;
	float aa, bb, cc;
	float interpolated_value;

	//Get the lower index
	ind = get_nearest_interp_index(y_source_coordinate, sid_depth, sid_y, s_coord_interp); 

	//Perform the interpolation
	if (s_val_interp[ind][sid_y][sid_depth] == s_val_interp[ind+1][sid_y][sid_depth]) {
		interpolated_value = s_val_interp[ind][sid_y][sid_depth]; //In the case of a horizontal line
	} else {
		aa = y_source_coordinate - s_coord_interp[ind][sid_y][sid_depth];
		bb = s_coord_interp[ind+1][sid_y][sid_depth] - s_coord_interp[ind][sid_y][sid_depth];
		cc = s_val_interp[ind+1][sid_y][sid_depth] - s_val_interp[ind][sid_y][sid_depth];
		interpolated_value = s_val_interp[ind][sid_y][sid_depth] + (aa/bb)*cc;
	}

	//Temporary diagnostics here
        /*if (blockIdx.x == 0 && blockIdx.y == 10 && blockIdx.z == 10 && threadIdx.x == 1 && threadIdx.y == 1 && threadIdx.z == 0) {
                printf("ind = %i, d_INTERP_ORDER = %i, aa = %e, bb = %e, cc = %e, interpolated_value = %e \n",
                ind, d_INTERP_ORDER, aa, bb, cc, interpolated_value);
        }*/

	return interpolated_value;
}


/**__device__ float interpolate(float y_source_coordinate, int sid_depth, int sid_y,
			     float s_coord_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH], 
                             float s_val_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH])
{
	// This kernel calls the chosen interpolation method.
	//TODO: Make as template.  

	switch (d_INTERP_METHOD)
        {
                case d_SWITCH_INTERP_LINEAR:
                        interpolated_value = 1.0; 
                        break;
                case d_SWITCH_INTERP_CUBIC_SPLINE:
			interpolated_value = 2.0;
                        break;
        }


	return interpolated_value;
}**/


__device__ void init_interp_array(float* array, int ix, int iz, int sid_depth, int sid_y, int iy_shifted,
                             float s_coord_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH], float s_val_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH])
{
	// Copies values from a device array to interpolation array in shared memory. 
	// The arrays consists of array values relevant to the block in question.

	int n_index = 0; //The current index int the copy loop. Corresponds to the main axis of the interpolation array.
	int send_idx; //The closest source grid points in y-direction based on d_DELTA_Y. Corresponds to the main axis of the interpolation array.

	int order_half = d_INTERP_ORDER/2;
	int im_start = iy_shifted - order_half; 
	int im_end = im_start + d_INTERP_ORDER;

	for (int im=im_start; im < im_end; im++) {
		send_idx =  ix + im*d_NX + iz*d_NX*d_NY; 
		s_coord_interp[n_index][sid_y][sid_depth] = (im - d_CY_BOT)*d_DY - d_YORIG; // Coordinates of the grid points.  
		s_val_interp[n_index][sid_y][sid_depth] = array[send_idx];                  // The corresponding values of the grid points. 
		//NOTE: Temporary debug printout. (TODO: Remove when not needed) 
        	/*if (blockIdx.x == 0 && blockIdx.y == 10 && blockIdx.z == 10 && threadIdx.x == 1 && threadIdx.y == 1 && threadIdx.z == 0) {
			printf("n_index = %i, im = %i, d_INTERP_ORDER = %i, iy_shifted = %i, s_coord_interp[%i][%i][%i] = %e, s_val_interp[%i][%i][%i] = %e, array[%i] = %e \n",
			n_index, im, d_INTERP_ORDER, iy_shifted, n_index, sid_y, sid_depth, s_coord_interp[n_index][sid_y][sid_depth],  
			n_index, sid_y, sid_depth, s_val_interp[n_index][sid_y][sid_depth], send_idx, array[send_idx]);
		}*/
                n_index++;
	} 
 
}

__device__ float get_y_shift(int iy_shift, int iy)
{
	// y_shift for shearbound_coord
	float y_shift;
	y_shift = ((float) (iy + iy_shift - d_CY_BOT))*d_DY - d_YORIG;

	return y_shift; 
}

__device__ void shearbound_coord(int* iy_shifted, float* y_source_coordinate, int ix, int iy)
{
	// Calculates the coordinate shift in boundary conditions according to shear.
	// y_source_coordinate is the coordinate in the source boundary 
	// iy_shifted is the nearest index to this coordinate.
	float y_shift, y_source_orig; 
	int iy_shift;
	
	iy_shift = 0; 

	y_source_orig = (((float) (iy - d_CY_BOT))*d_DY - d_YORIG); 

	// Check which side of the shearing box we are looking for
	if (ix < d_NX/2) {
		// Physical coordinate of the source of for X BOTTOM boundary
		(*y_source_coordinate) = y_source_orig - d_DELTA_Y;
		// Check if the boundary is crossed. The shift is periodic in y direction.
		if ((*y_source_coordinate) < (0.0 - d_YORIG)) { 
			(*y_source_coordinate) = (*y_source_coordinate) + d_DOMAIN_SIZE_Y;
		}
		//Nearest index in y direction
		y_shift = get_y_shift(iy_shift, iy); 
		while ( y_shift > (*y_source_coordinate)) {
			iy_shift -= 1; 
                        y_shift = get_y_shift(iy_shift, iy);
		}
	}
	else if (ix > d_NX/2) {
		//Physical coordinate of the source of for X TOP boundary
		(*y_source_coordinate) = y_source_orig + d_DELTA_Y;
		// Check if the boundary is crossed. The shift is periodic in y direction.
		if ((*y_source_coordinate) > (d_DOMAIN_SIZE_Y - d_YORIG)) {
			(*y_source_coordinate) = (*y_source_coordinate) - d_DOMAIN_SIZE_Y;
		}
		//Nearest index in y direction
		y_shift = get_y_shift(iy_shift, iy);
		while ( y_shift < (*y_source_coordinate)) {
			iy_shift += 1; 
			y_shift = get_y_shift(iy_shift, iy); 
		}
	} else {
		printf("Error in shearbound in shear.cu!!!\n");
	}

	(*iy_shifted) = iy + iy_shift; //from gdb: seems to work ok 

	/*if (blockIdx.x == 0 && blockIdx.y == 10 && blockIdx.z == 10 && threadIdx.x == 1 && threadIdx.y == 1 && threadIdx.z == 0) {
		printf("iy = %i, iy_shift = %i, iy_shifted = %i, y_source_orig = %e, y_source_coordinate = %e \n", 
			iy, iy_shift, (*iy_shifted), y_source_orig, (*y_source_coordinate));
	}*/
}

__device__ float interp_shear(float* array, int ix, int iy, int iz, int sid_depth, int sid_y, 
		             float s_coord_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH], float s_val_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH])
{
	//TODO: Adapt linear and cubic spline interpollation for shear here. 

	//Step 1: Figure the correct point of interpolation
	int iy_shifted; //Nearest y-index to the interpollation poin
	float y_source_coordinate; //The physical coordinate of the interpollation point
	float interpolated_value; //Result of the interpolation for the thread

	shearbound_coord(&iy_shifted, &y_source_coordinate, ix, iy);

	__syncthreads(); // To avoid unsychronized reads and writes of the shared memory blocks.  

	//Step 2: Init s_coord_interp and s_val_interp for the current memory block
	init_interp_array(array, ix, iz, sid_depth, sid_y, iy_shifted, s_coord_interp, s_val_interp);
	
	__syncthreads();	

	//Step 3: Apply the interpolation method base 
	//TODO: Need a swich when we have different interpolators
	
	//Perform linear interpolation
	interpolated_value = lin_interp(y_source_coordinate, sid_depth, sid_y, s_coord_interp, s_val_interp);
	
	__syncthreads();

	//Step 4: Return the interpolated value connected to the thread
	
	return interpolated_value; 
} 

__device__ void epicyclic(float* u_epic_y, int sid_row, int sid_column, int sid_depth, float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH])
{
	/*
	Calculate the epicyclic motion component into the Navier-Stokes equation. 
	*/

	*u_epic_y = -d_Q_SHEAR*d_OMEGA*s_uu_x[sid_row][sid_column][sid_depth];
}

__device__ void shear_flow(float* u_shear_y, float grid_idx_x)
{
	/*
	Calculate the shear flow component into the velocity. Shear allowed only in y-direction. 
	*/

	float coord_x;

	coord_x = d_DX*grid_idx_x - d_XORIG;

	*u_shear_y = -d_Q_SHEAR*d_OMEGA*coord_x;

}
