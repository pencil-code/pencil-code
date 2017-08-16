//
// Shear flow effects are calculated in this module 
//

#include <cstdio>

#include "dconstsextern.cuh"
#include "shear.cuh"
#include "smem.cuh"

// TODO: Adapt? Pencil Code method for cubic spline interpolation. 
// TODO: Examine Pencil Code method for polynomial interpolation.

/*
__device__ float poly_interpolate(float s_coord_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH], float s_val_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH], 
			    float C_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH], float D_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH], 
			    int sid_depth, int sid_y, int im_start, int im_end, float y_source_coordinate, float *y_interp, float* dy_interp, int* ierr, float lprint)
{
	// TODO: Spline interpollation???
	//  Polynomial interpolation adapted from numerical recipes. (Assuming n:th order)
	//  shear.f90:  integer :: norder_poly = 3 (The Pencil Code assumes 3th order)
	//  morder = max(norder, 0)
	//  moh = morder / 2    
	//  ix1 = ix - moh
	//  ix2 = ix + moh
	//  call poly_interp_one(xa(ix1:ix2), ya(ix1:ix2), x(i), y(i), dy(i), istat, msg)
	//
	int closest_ind;
   	float source_diff, nearest_to_interp_coord, new_nearest_to_interp_coord;
	float source_diff_i, source_diff_im, weight;

	nearest_to_interp_coord = fabs(y_source_coordinate-s_coord_interp[0][sid_y][sid_depth]);

	closest_ind = 0;
	for (int i = 0; i <= d_INTERP_ORDER; i++) {
		// Find the index of the closest table entry 
		new_nearest_to_interp_coord = fabs(y_source_coordinate-s_coord_interp[i][sid_y][sid_depth]);
		if (new_nearest_to_interp_coord < nearest_to_interp_coord) {
			closest_ind=i;
			nearest_to_interp_coord = new_nearest_to_interp_coord;
		}
		// Initialize c:s and d:s 
		C_interp[i][sid_y][sid_depth] = s_val_interp[i][sid_y][sid_depth];
		D_interp[i][sid_y][sid_depth] = s_val_interp[i][sid_y][sid_depth];
	}

	// So that index will not go bellow zero
	int initial;
	if (closest_ind - 1 >= 0 ) { 
		initial = closest_ind - 1; 
	} else {
		initial = 0;
	} 
   
	*y_interp = s_val_interp[initial][sid_y][sid_depth]; // The initial approximation of y

	// For each column of the tableau, we loop over current c's and d's and update them
	im_start = im_start + 1 ; //TODO: Why this???????!!!
	im_end = im_end + 1  ;    //TODO: Why this???????!!!
	for (int m=im_start; m < im_end; m++) {
		for (int i=im_start; i <= im_end-m; i++) {
			source_diff_i = s_coord_interp[i-1][sid_y][sid_depth]   - y_source_coordinate;
			source_diff_im = s_coord_interp[i+m-1][sid_y][sid_depth] - y_source_coordinate;
			weight  =  C_interp[i][sid_y][sid_depth] - D_interp[i-1][sid_y][sid_depth];
			source_diff = source_diff_i-source_diff_im;
			if (source_diff == 0.0) {
				*ierr = 1;
			}
			if (lprint == 1) { //TODO: Restrict the thread!!!
				printf("Error in the routine poly_interpolate! \n");
				printf("m = %i, i = %i, d_INTERP_ORDER = %i, d_INTERP_ORDER-m = %i.\n", m, i, d_INTERP_ORDER, d_INTERP_ORDER-m);
				printf("s_coord_interp[i-1][sid_y][sid_depth] = %f, s_coord_interp[i+m-1][sid_y][sid_depth] = %f, y_source_coordinate = %f, C[i] = %f, D[i-1] = %f, source_diff_i = %f source_diff_im = %f \n",
					s_coord_interp[i-1][sid_y][sid_depth], s_coord_interp[i+m-1][sid_y][sid_depth], y_source_coordinate, C[i], D[i-1], source_diff_i,
					source_diff_im);
			}
			source_diff = source_diff/weight;
			C_interp[i-1][sid_y][sid_depth] = source_diff_im*source_diff;
			D_interp[i-1][sid_y][sid_depth] = source_diff_i*source_diff;
			if (lprint == 1) {
				printf("The new values:\n");
				printf("C_interp[i-1][sid_y][sid_depth] = %f, D_interp[i-1][sid_y][sid_depth] = %f, source_diff = %f \n\n",
					C_interp[i-1][sid_y][sid_depth],      D_interp[i-1][sid_y][sid_depth],      source_diff);
			}
		}
		// Choose C or D depending on which correction (up or low)
		// is better for the accumulated value.
		// the last dy added is the interpolation error estimate.
		if (2*closest_ind < d_INTERP_ORDER-m) {
			*dy_interp = C_interp[closest_ind][sid_y][sid_depth];
		} else {
			*dy_interp = D_interp[closest_ind-2][sid_y][sid_depth];
		}

		if (lprint == 1) {
			printf("The shift:\n");
			printf("dy_interp = %f \n\n", *dy_interp);
		}      

		*y_interp += *dy_interp;

		if (lprint == 1) {
			printf("The result in step m = %i: \n", m);
			printf("y_interp = %f \n\n", *y_interp);
		}            
	}
 
	*dy_interp = fabs(*dy_interp); // Take the absolute value of the interpolation error.

	if (*dy_interp == 0.0) {
		*ierr = 1;
	}

	free(C); free(D); 

}*/


__device__ void interp_limits(int* interp_start, int* interp_end, int* interp_order_local, int sid_depth, int sid_y, 
                              float s_val_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH])
{
	//   
	// Make interpolation arrays shorter if there are repeating values. The interpolation methon cannot handle such things,  
	// as they cause div by 0. 
	//

	//Look at bottom  s_val_interp[0][sid_y][sid_depth]
	int n;
	n = 1;
	while (n < d_INTERP_ORDER && *interp_order_local > 1) {
		if (s_val_interp[n][sid_y][sid_depth] == s_val_interp[n-1][sid_y][sid_depth]) {
			*interp_start++;
			*interp_order_local--;
			n++;
		} else {
			n = d_INTERP_ORDER; //Terminates the loop
		}
	} 

	// TODO: Cannot yet handle constant values that might appear in the middle!!! Unnecessary???
	//Look at top 
  	n = *interp_end-1;
	while (n >= 0 && *interp_order_local > 1 && *interp_end > *interp_start) {
		if (s_val_interp[n][sid_y][sid_depth] == s_val_interp[n+1][sid_y][sid_depth]) {
			*interp_end--;
			*interp_order_local--;
			n--;
		} else {
			n = -1; //Terminates the loop
		}
	} 

}

__device__ float get_y_shift(int iy_shift, int iy)
{
	//y_shift for shearbound_coord
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

	if (blockIdx.x == 0 && blockIdx.y == 10 && blockIdx.z == 10 && threadIdx.x == 1 && threadIdx.y == 1 && threadIdx.z == 0) {
		printf("iy = %i, iy_shift = %i, iy_shifted = %i, y_source_orig = %e, y_source_coordinate = %e \n", 
			iy, iy_shift, (*iy_shifted), y_source_orig, (*y_source_coordinate));
	}
}


__device__ void interp_shear(float* array, int grid_idx, int bound_idx, int ix, int iy, int iz, int ix_bound, int sid_depth, int sid_y, 
		             float s_coord_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH], float s_val_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH])
{
	//TODO: MODIFY FOR CUDA
	//   Handles the interpolated shearing boundary condition for one thread 
   
   	int send_idx, iy_shifted;
   	float y_source_coordinate, array_value_y_interp, dy_interp; 
   
	// Calculate the shift in coordinates and the nearest index according to delta_y
	shearbound_coord(&iy_shifted, &y_source_coordinate, ix, iy);

	// Initialize interpolation vectors (in y-direction) based on the new inds. indr should be unchanged
	// TODO: Into a function?
	int n, interp_order_local, n_shift;
	n=0;
	n_shift = 0;
	interp_order_local = d_INTERP_ORDER; //If we shorten the interpolation

	int order_half = d_INTERP_ORDER/2;
	int im_start = iy_shifted - order_half; 
	int im_end = im_start + d_INTERP_ORDER;

	for (int im=im_start; im < im_end; im++) {
		send_idx =  ix + im*d_NX + iz*d_NX*d_NY; //The closest source grid points in y-direction based on d_DELTA_Y
		s_coord_interp[n][sid_y][sid_depth] = (im - d_CY_BOT)*d_DY - d_YORIG; // Coordinates of the grid points. d_CY_BOT,d_DY, d_YORIG work ok 
		s_val_interp[n][sid_y][sid_depth] = array[send_idx];               // Values of the grid points. 
        	if (blockIdx.x == 0 && blockIdx.y == 10 && blockIdx.z == 10 && threadIdx.x == 1 && threadIdx.y == 1 && threadIdx.z == 0) {
			printf("n = %i, im = %i, d_INTERP_ORDER = %i, iy_shifted = %i, s_coord_interp[%i][%i][%i] = %e, s_val_interp[%i][%i][%i] = %e, array[%i] = %e \n", 
			n, im, d_INTERP_ORDER, iy_shifted, n, sid_y, sid_depth, s_coord_interp[n][sid_y][sid_depth],  n, sid_y, sid_depth, s_val_interp[n][sid_y][sid_depth],
			send_idx, array[send_idx]);
		}
                n++;
	} 

	//Diagnostics TODO: optional
	if (blockIdx.x == 0 && blockIdx.y == 10 && blockIdx.z == 10 && threadIdx.x == 1 && threadIdx.y == 1 && threadIdx.z == 0) {
		printf("s_coord_interp[0 : 5][%i][%i] = [%e, %e, %e, %e, %e, %e] \ns_val_interp[0 : 5][%i][%i] = [%e, %e, %e, %e, %e, %e] \ny_source_coordinate = %e, iy = %i, iy_shifted = %i \n", 
			sid_y, sid_depth, 
                        s_coord_interp[0][sid_y][sid_depth], s_coord_interp[1][sid_y][sid_depth], 
			s_coord_interp[2][sid_y][sid_depth], s_coord_interp[3][sid_y][sid_depth], 
			s_coord_interp[4][sid_y][sid_depth], s_coord_interp[5][sid_y][sid_depth], 
			sid_y, sid_depth, 
                        s_val_interp[0][sid_y][sid_depth], s_val_interp[1][sid_y][sid_depth], 
			s_val_interp[2][sid_y][sid_depth], s_val_interp[3][sid_y][sid_depth], 
			s_val_interp[4][sid_y][sid_depth], s_val_interp[5][sid_y][sid_depth], 
			y_source_coordinate, iy, iy_shifted);
	} 


	//Find limits for interpolation vector if s_val_interp has repeating values
	//The interpolation method cannot handle "constants"
	int interp_start, interp_end; //Movable starting and ending points
	//Set the defaults to whole array 
	interp_start = 0;
   	interp_end = d_INTERP_ORDER-1;
	//Diagnostics TODO: optional
	if (blockIdx.x == 0 && blockIdx.y == 10 && blockIdx.z == 10 && threadIdx.x == 1 && threadIdx.y == 1 && threadIdx.z == 0) {
		printf("Before checking limits: interp_start = %i,  interp_end = %i, interp_order_local = %i\n" , interp_start,  interp_end, interp_order_local);  
	}
	//Adjust limists
	interp_limits(&interp_start, &interp_end, &interp_order_local, sid_depth, sid_y, s_val_interp); 
	//Diagnostics TODO: optional
	if (blockIdx.x == 0 && blockIdx.y == 10 && blockIdx.z == 10 && threadIdx.x == 1 && threadIdx.y == 1 && threadIdx.z == 0) {
		printf("After checking limits: interp_start = %i,  interp_end = %i, interp_order_local = %i\n\n" , interp_start,  interp_end, interp_order_local);  
	}


	// Do the interpolation
	int less_that_precision, ierr;
	ierr = 0;

	//Try linear interpollation TODO: Make a __device__ kernel
	float val_interp;
	if (s_val_interp[order_half][sid_y][sid_depth] == s_val_interp[order_half+1][sid_y][sid_depth]) { 
		val_interp = s_val_interp[order_half][sid_y][sid_depth];
	} else {
		val_interp = s_val_interp[order_half][sid_y][sid_depth] + ((s_coord_interp[order_half][sid_y][sid_depth] - y_source_coordinate)
                                                                          /(s_coord_interp[order_half+1][sid_y][sid_depth] - s_coord_interp[order_half][sid_y][sid_depth]))
                                                                          *(s_val_interp[order_half+1][sid_y][sid_depth] - s_val_interp[order_half][sid_y][sid_depth]) ;
	}
 
	/*
	if (interp_order_local > 1) {
		if ( s_coord_interp[interp_start][sid_y][sid_depth] < y_source_coordinate && s_coord_interp[interp_end][sid_y][sid_depth] > y_source_coordinate ) {
		//TODO poly_interpolate(xa_interp, s_val_interp, interp_start, interp_end, d_INTERP_ORDER, y_source_coordinate, &array_value_y_interp, &dy_interp, &ierr, 0);
		// TODO: Do we need to check for dy_interp that our interpolation is precise enough?
		// TODO: abs(ya_interp[n] - ya_interp[n-1]) > dy_interp then should be treated as constant
		//TODO Check of variation with ya_interp is higher than dy_interp between interp_start, interp_end
		//TODO If not, asume array_value_y_interp = ya_interp[0], as bellow 
		//TODO less_that_precision = 0;
		//TODO: Figure out what in hell this actually does 
		//TODO for (int m=interp_start; m <= interp_end; m++) {
		//TODO	for (int l=interp_start; l < m; l++) {
		//TODO 		if (fabs(s_val_interp[m][sid_y][sid_depth] - s_val_interp[l][sid_y][sid_depth]) < dy_interp) {
		//TODO 			less_that_precision++;
		//TODO 		}
		//TODO 	}
		//TODO 	for (int l=m+1; l <= interp_end; l++) {
		//TODO 		if (fabs(s_val_interp[m][sid_y][sid_depth] - s_val_interp[l][sid_y][sid_depth]) < dy_interp) {
		//TODO 		less_that_precision++;
		//TODO 		}
		//TODO 	}
		//TODO }
		//TODO if (less_that_precision >= interp_order_local-1 || ierr == 1) {
		//TODO 	//Assume the values are practically constant. Choose the middle value
			//TODO 	array_value_y_interp = s_val_interp[(interp_start + interp_end)/2][sid_y][sid_depth]; 
			//TODO }
		} else {
			//TODO: What then???
		}
	} else {
		// If interpolation is skipped, i.e. s_val_interp is constant. Assing a value from one point in ya
		array_value_y_interp = s_val_interp[0][sid_y][sid_depth];
	} 
	*/

/*
	// Assign the interpolated values to bound_idx. 
	array[bound_idx] = array_value_y_interp;

	//TODO: Make into a function
	// Check for mysterious sight changes. Print diagnotics only then. 
	ierr = 0;
	for (n=0; n < interp_order_local; n++) {
		if (s_ya_interp[n] > 0.0 && array_value_y_interp < 0.0) {
			ierr++;
		} 
		else if (s_ya_interp[n] < 0.0 && array_value_y_interp > 0.0) {  
			ierr++;
		}
	}
      
	//Print error messages concerning the sign changes
	if (ierr >= d_INTERP_ORDER) {
		printf("Sign error spotted! \n");
		printf("ix = %i, iy = %i, iz = %i, \n", 
			ix     , iy     , iz);   
		for (n=0; n<interp_order_local; n++) {
			printf("s_xa_interp[%i] = %e  ",
			n, s_xa_interp[n]);
		}
		printf("\n");
		for (n=0; n<interp_order_local; n++) {
			printf("s_ya_interp[%i] = %e  ",
			n, s_ya_interp[n]);
		}
		printf("\n");
		printf("d_INTERP_ORDER/2 = %i, interp_order_local = %i, iy_shifted = %i, y_source_coordinate = %f array_value_y_interp = %e, dy_interp = %e, less_that_precision = %i \n\n", 
		d_INTERP_ORDER/2,      interp_order_local,      iy_shifted,      y_source_coordinate,     array_value_y_interp,      dy_interp,      less_that_precision);

		//TODOpoly_interpolate(s_xa_interp, s_ya_interp, interp_start, interp_end, d_INTERP_ORDER, y_source_coordinate, &array_value_y_interp, &dy_interp, &ierr, 1);
	}
*/
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
