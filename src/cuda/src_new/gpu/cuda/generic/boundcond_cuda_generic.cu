/*
*	This module is used to copy the cover of the computational domain to
*	the appropriate boundary zones. For reference;
*		-UP is used when moving to the positive direction along y-axis (+y)
*		-RIGHT is used when moving to the positive direction along x-axis (+x)
*		-BACK is used when moving to the positive direction along z-axis (+z)
*
*	Check astaroth-code/doc/boundcond.png for an illustrative picture
* 	(Additionally, the logic behind "per_xy_edges" etc. is such that if both
*	x and y are periodic, then the edges between CZ_BOT and CZ_TOP
*	can be copied. This is based on the assumption, that if z is not periodic,
*	then special rules are applied when copying data to the planes at z coordinates 
*	CZ_BOT-BOUND_SIZE ... CZ_BOT and CZ_TOP ... CZ_TOP+BOUND_SIZE)
*/
#include "boundcond_cuda_generic.cuh"
#include "gpu/cuda/core/dconsts_core.cuh"
#include "gpu/cuda/core/errorhandler_cuda.cuh"

//Copies the front and back of the computational domain to an appropriate
//boundary zone (does not include the edges and corners of the boundary zone)
__global__ void per_z_sides(Grid d_grid)
{
	int iz, iz_bound;
	if (blockIdx.z < 3) { //Copy front of the computational domain to the boundary zone at the back
		iz = blockIdx.z + d_nz_min;
		iz_bound = blockIdx.z + d_nz_max;
	} 
	else { //Copy back of the computational domain to the boundary zone at the front
		iz = (blockIdx.z-3) + d_nz_max - 3;
		iz_bound = (blockIdx.z-3);
	}

	int ix,iy;
	ix = threadIdx.x + blockIdx.x*blockDim.x + d_nx_min;
	iy = threadIdx.y + blockIdx.y*blockDim.y + d_ny_min;

	if (ix < d_nx_max && iy < d_ny_max) {

		int grid_idx =  ix + iy*d_mx + iz*d_mxy;
		int bound_idx =  ix + iy*d_mx + iz_bound*d_mxy;

        for (int i=0; i < d_grid.NUM_ARRS; ++i)
            d_grid.arr[i][bound_idx] = d_grid.arr[i][grid_idx]; 
	}
}

//Copies the top and bottom of the computational domain to an appropriate
//boundary zone (does not include the edges and corners of the boundary zone)
__global__ void per_y_sides(Grid d_grid)
{
	int iy, iy_bound;
	if (blockIdx.z < 3) { //Copy bottom of the computational domain to the boundary zone at the top
		iy = blockIdx.z + d_ny_min;
		iy_bound = blockIdx.z + d_ny_max;
	} 
	else { //Copy top of the computational domain to the boundary zone at the bottom
		iy = (blockIdx.z-3) + (d_ny_max - 3);
		iy_bound = (blockIdx.z-3);
	}

	int ix,iz;
	ix = threadIdx.x + blockIdx.x*blockDim.x + d_nx_min;
	iz = threadIdx.y + blockIdx.y*blockDim.y + d_nz_min;

	if (ix < d_nx_max && iz < d_nz_max) {

		int grid_idx =  ix + iy*d_mx + iz*d_mxy;
		int bound_idx =  ix + iy_bound*d_mx + iz*d_mxy;

        for (int i=0; i < d_grid.NUM_ARRS; ++i)
            d_grid.arr[i][bound_idx] = d_grid.arr[i][grid_idx]; 
	}
}


//Copies the left and right sides of the computational domain to an appropriate
//boundary zone (does not include the edges and corners of the boundary zone)
__global__ void per_x_sides(Grid d_grid)
{

        int ix, ix_bound;
        if (threadIdx.x < 3) { //Copy left of the computational domain to the boundary zone at the right
                ix = threadIdx.x + d_nx_min;
                ix_bound = threadIdx.x + d_nx_max;
        }
        else { //Copy right of the computational domain to the boundary zone at the left
                ix = (threadIdx.x-3) + (d_nx_max - 3);
                ix_bound = (d_nx_min-3) +(threadIdx.x-3);
        }

        int iy,iz;
        iz = threadIdx.z + blockIdx.z*blockDim.z + d_nz_min;//Don't add edges
        iy = threadIdx.y + blockIdx.y*blockDim.y + d_ny_min;

        int grid_idx =  ix + iy*d_mx + iz*d_mxy;
        int bound_idx =  ix_bound + iy*d_mx + iz*d_mxy;

        for (int i=0; i < d_grid.NUM_ARRS; ++i)
            d_grid.arr[i][bound_idx] = d_grid.arr[i][grid_idx]; 
    /*//Uncommented this TODO recheck if causes issues
	//Normal periodic boundary if shearing is not included
    int sid_depth = threadIdx.x;
    int sid_y = threadIdx.y;

	if (d_LSHEAR == 0 || (d_DELTA_Y > -1.0e-30 && d_DELTA_Y < 1.0e-30) ) { //Meanign if d_DELTA_Y == 0, but reals are not that stable. 
		if (iy < d_ny_max && iz < d_nz_max) {

			int grid_idx =  ix + iy*d_mx + iz*d_mxy;
			int bound_idx =  ix_bound + iy*d_mx + iz*d_mxy;

			d_lnrho[bound_idx] = d_lnrho[grid_idx];
			d_uu_x[bound_idx] = d_uu_x[grid_idx];
			d_uu_y[bound_idx] = d_uu_y[grid_idx];
			d_uu_z[bound_idx] = d_uu_z[grid_idx];
		}
	} else {
		if (iy < d_ny_max && iz < d_nz_max) {
			int bound_idx =  ix_bound + iy*d_mx + iz*d_mxy;

			//Allocate shared memory for interpolation arrays 
			__shared__ real s_coord_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH];
			__shared__ real s_val_interp[INTERP_ORDER][INTERP_NPOINTS][INTERP_DEPTH];

			//Perform the interpolation in assign the values to the boundaries
			d_lnrho[bound_idx] = interp_shear(d_lnrho, ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);
			d_uu_x[bound_idx]  = interp_shear(d_uu_x,  ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);
			d_uu_y[bound_idx]  = interp_shear(d_uu_y,  ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);
			d_uu_z[bound_idx]  = interp_shear(d_uu_z,  ix, iy, iz, sid_depth, sid_y, s_coord_interp, s_val_interp);

		}
	}
    */

}


//Copy the edges from upper front & back and bottom front & back to
//the appropriate boundary zones
//(Requires thread dims of (32, 3, 3) and blockDims of (ceil((real) d_nx / (real)tpb.x), 1, 4)
__global__ void per_yz_edges(Grid d_grid)
{
	int ix, iy, iz;
	int grid_idx, bound_idx;
	ix = threadIdx.x + blockIdx.x*blockDim.x + d_nx_min; //x index skips the boundary and starts from the computational domain
	iy = threadIdx.y + d_ny_min; 
	iz = threadIdx.z + d_nz_min;

	switch(blockIdx.z)
	{
		case 0: //Copy upper front edge of the computational domain to the boundary zone at bottom back
			grid_idx = ix + (iy + d_ny-3)*d_mx + iz*d_mxy;
			bound_idx = 	ix + 
					(iy-BOUND_SIZE)*d_mx + 
					(iz+d_nz)*d_mxy;
			break;

		case 1: //Copy bottom front edge of the computational domain to the boundary zone at upper back
			grid_idx = ix + iy*d_mx + iz*d_mxy;
			bound_idx = 	ix + 
					(iy+d_ny)*d_mx + 
					(iz+d_nz)*d_mxy;
			break;

		case 2: //Copy upper back edge of the computational domain to the boundary zone at bottom front
			grid_idx = ix + (iy + d_ny-3)*d_mx + (iz + d_nz-3)*d_mxy;
			bound_idx = 	ix + 
					(iy-BOUND_SIZE)*d_mx + 
					(iz-BOUND_SIZE)*d_mxy;
			break;

		case 3: //Copy bottom back edge of the computational domain to the boundary zone at upper front
			grid_idx = ix + iy*d_mx + (iz + d_nz-3)*d_mxy;
			bound_idx = 	ix + 
					(iy+d_ny)*d_mx + 
					(iz-BOUND_SIZE)*d_mxy;
			break;


	}
	
	if (ix < d_nx_max) {
        for (int i=0; i < d_grid.NUM_ARRS; ++i)
            d_grid.arr[i][bound_idx] = d_grid.arr[i][grid_idx]; 
	}
}

//Copy the edges from front left & right and back left & right to
//the appropriate boundary zones
//(Requires thread dims of (3, 32, 3) and blockDims of (1, ceil((real) d_ny / (real)tpb.y), 4))
__global__ void per_xz_edges(Grid d_grid)
{
	int ix, iy, iz;
	int grid_idx, bound_idx;
	ix = threadIdx.x + d_nx_min; 
	iy = threadIdx.y + blockIdx.y*blockDim.y + d_ny_min;
	iz = threadIdx.z + d_nz_min; 

	switch(blockIdx.z)
	{
		case 0: //Copy front left edge of the computational domain to the boundary zone at right back
			grid_idx = ix + iy*d_mx + iz*d_mxy;
			bound_idx = 	(ix + d_nx) + 
					iy*d_mx + 
					(iz+d_nz)*d_mxy;
			break;

		case 1: //Copy right front edge of the computational domain to the boundary zone at left back
			grid_idx = (ix + d_nx-3) + iy*d_mx + iz*d_mxy;
			bound_idx = 	(ix-BOUND_SIZE) + 
					iy*d_mx + 
					(iz+d_nz)*d_mxy;
			break;

		case 2: //Copy left back edge of the computational domain to the boundary zone at right front
			grid_idx = ix + iy*d_mx + (iz + d_nz-3)*d_mxy;
			bound_idx = 	(ix + d_nx) + 
					iy*d_mx + 
					(iz-BOUND_SIZE)*d_mxy;
			break;

		case 3: //Copy right back edge of the computational domain to the boundary zone at left front
			grid_idx = (ix + d_nx-3) + iy*d_mx + (iz + d_nz-3)*d_mxy;
			bound_idx = 	(ix-BOUND_SIZE) + 
					iy*d_mx + 
					(iz-BOUND_SIZE)*d_mxy;
			break;


	}
	
	if (iy < d_ny_max) {
        for (int i=0; i < d_grid.NUM_ARRS; ++i)
            d_grid.arr[i][bound_idx] = d_grid.arr[i][grid_idx]; 
	}
}


//Copy the edges from upper left & right and bottom left & right to
//the appropriate boundary zones
//(Requires thread dims of (3, 3, 32) and blockDims of (1, 4, ceil((real) d_nz / (real)tpb.z)))
__global__ void per_xy_edges(Grid d_grid)
{
	int ix, iy, iz;
	int grid_idx, bound_idx;
	ix = threadIdx.x + d_nx_min; 
	iy = threadIdx.y + d_ny_min; 
	iz = threadIdx.z + blockIdx.z*blockDim.z + d_nz_min;

	switch(blockIdx.y)
	{
		case 0: //Copy upper left edge of the computational domain to the boundary zone at bottom right
			grid_idx = ix + (iy + d_ny-3)*d_mx + iz*d_mxy;
			bound_idx = 	(ix + d_nx) + 
					(iy-BOUND_SIZE)*d_mx + 
					iz*d_mxy;
			break;

		case 1: //Copy upper right edge of the computational domain to the boundary zone at bottom left
			grid_idx = (ix + d_nx-3) + (iy + d_ny-3)*d_mx + iz*d_mxy;
			bound_idx = 	(ix-BOUND_SIZE) + 
					(iy-BOUND_SIZE)*d_mx + 
					iz*d_mxy;
			break;

		case 2: //Copy bottom left edge of the computational domain to the boundary zone at upper right
			grid_idx = ix + iy*d_mx + iz*d_mxy;
			bound_idx = 	(ix + d_nx) + 
					(iy + d_ny)*d_mx + 
					iz*d_mxy;
			break;

		case 3: //Copy bottom right edge of the computational domain to the boundary zone at upper left
			grid_idx = (ix + d_nx-3) + iy*d_mx + iz*d_mxy;
			bound_idx = 	(ix-BOUND_SIZE) + 
					(iy + d_ny)*d_mx + 
					iz*d_mxy;
			break;


	}
	
	if (iz < d_nz_max) {
        for (int i=0; i < d_grid.NUM_ARRS; ++i)
            d_grid.arr[i][bound_idx] = d_grid.arr[i][grid_idx]; 
	}
}


//Copies the corners of the computational domain to appropriate boundary areas
//Uses x,y,z to determine the index inside the block and blockIdx.z to determine 
//which one of the eight corners to copy.  
//(Requires thread dims of (3, 3, 3) and blockDims of (1, 1, 8))
__global__ void per_xyz_corners(Grid d_grid)
{
	int ix, iy, iz;
	int grid_idx, bound_idx;
	ix = threadIdx.x + d_nx_min; 
	iy = threadIdx.y + d_ny_min;
	iz = threadIdx.z + d_nz_min; 

	switch(blockIdx.z)
	{
		case 0: //Copy the bottom left front corner to boundary zone at upper right back (x=0, y=0, z=0)
			grid_idx = ix + iy*d_mx + iz*d_mxy;
			bound_idx = 	(ix + d_nx) + 
					(iy + d_ny)*d_mx + 
					(iz + d_nz)*d_mxy;
			break;

		case 1: //Copy the bottom left back corner to boundary zone at upper right front (x=0, y=0, z=1)
			grid_idx = ix + iy*d_mx + (iz+d_nz-3)*d_mxy;
			bound_idx = 	(ix + d_nx) + 
					(iy + d_ny)*d_mx + 
					(iz - BOUND_SIZE)*d_mxy;
			break;

		case 2: //Copy the upper left front corner to boundary zone at bottom right back (x=0, y=1, z=0)
			grid_idx = ix + (iy+d_ny-3)*d_mx + iz*d_mxy;
			bound_idx = 	(ix + d_nx) + 
					(iy - BOUND_SIZE)*d_mx + 
					(iz + d_nz)*d_mxy;
			break;

		case 3: //Copy the upper left back corner to boundary zone at bottom right front (x=0, y=1, z=1)
			grid_idx = ix + (iy+d_ny-3)*d_mx + (iz+d_nz-3)*d_mxy;
			bound_idx = 	(ix + d_nx) + 
					(iy - BOUND_SIZE)*d_mx + 
					(iz - BOUND_SIZE)*d_mxy;
			break;

		case 4: //Copy the bottom right front corner to boundary zone at upper left back (Do x=1, y=0, z=0)
			grid_idx = (ix+d_nx-3) + iy*d_mx + iz*d_mxy;
			bound_idx = 	(ix - BOUND_SIZE) + 
					(iy + d_ny)*d_mx + 
					(iz + d_nz)*d_mxy;
			break;

		case 5: //Copy the bottom right back corner to boundary zone at upper left front (x=1, y=0, z=1)
			grid_idx = (ix+d_nx-3) + iy*d_mx + (iz+d_nz-3)*d_mxy;
			bound_idx = 	(ix - BOUND_SIZE) + 
					(iy + d_ny)*d_mx + 
					(iz - BOUND_SIZE)*d_mxy;
			break;

		case 6: //Copy the upper right front corner to boundary zone at bottom left back (x=1, y=1, z=0)
			grid_idx = (ix+d_nx-3) + (iy+d_ny-3)*d_mx + iz*d_mxy;
			bound_idx = 	(ix - BOUND_SIZE) + 
					(iy - BOUND_SIZE)*d_mx + 
					(iz + d_nz)*d_mxy;
			break;

		case 7: //Copy the upper right back corner to boundary zone at bottom left front (x=1, y=1, z=1)
			grid_idx = 	(ix+d_nx-3) + 
					(iy+d_ny-3)*d_mx + 
					(iz+d_nz-3)*d_mxy;
			
			bound_idx = 	(ix - BOUND_SIZE) + 
					(iy - BOUND_SIZE)*d_mx + 
					(iz - BOUND_SIZE)*d_mxy;
			break;
	}

    for (int i=0; i < d_grid.NUM_ARRS; ++i)
        d_grid.arr[i][bound_idx] = d_grid.arr[i][grid_idx]; 
}


//Define boundcond types
#define PERIODIC_BOUNDCONDS 0
#define SHEARING_BOUNDCONDS 1
#define BOUNDCOND_TYPE_X PERIODIC_BOUNDCONDS
#define BOUNDCOND_TYPE_Y PERIODIC_BOUNDCONDS
#define BOUNDCOND_TYPE_Z PERIODIC_BOUNDCONDS
void boundcond_cuda_generic(Grid* d_grid, CParamConfig* cparams, cudaStream_t stream)
{
	//Quick summary:
	//The point in a 3D cuboid is copied to a location, where the location index is
	//offset in 1, 2 or 3 axes
	//f.ex.
	//	-Points that are copied by adding an offset in only one axis 
	//	(for example from the front of the computational domain to the boundary in the back)
	//	(Functions: per_z_sides, per_x_sides, per_y_sides)
	//	
	//	-Points that are offset in two axes, for example the top left edge (not including the corner)
	//	of the computational domain is copied to the boundary to the bottom right of the grid
	//	(Functions: per_xy_edges per_xz_edges per_yz_edges)
	//
	//	-Points that are offset in all three axes, e.g. the corners. For example the front top right
	//	3*3*3 cube of the computational domain is copied to the boundary zone in back bottom left in the grid. 
	//	(Function: per_xyz_corners)
	//
	// BOUNDCOND_TYPE_X, BOUNDCOND_TYPE_Y and BOUNDCOND_TYPE_Z are used to determine how 
	// the boundaries in their respective axis are supposed to be copied.

	//--------X BOUNDS---------------
	switch	(BOUNDCOND_TYPE_X) {
		case PERIODIC_BOUNDCONDS: {
			//Copy periodic x sides
            const dim3 tpb(6, 4, 1);
            const dim3 bpg(
                       1, 
                       (unsigned int)ceil((real) cparams->ny / (real)tpb.y), 
                       (unsigned int)ceil((real) cparams->nz / (real)tpb.z));
			per_x_sides<<<bpg, tpb, 0, stream>>>(*d_grid);
			CUDA_ERRCHK_KERNEL();

			//Copy periodic xy edges
			if (BOUNDCOND_TYPE_Y == PERIODIC_BOUNDCONDS) {
                const dim3 tpb(3, 3, 32);
                const dim3 bpg(
                           1, 
                           4, 
                           (unsigned int)ceil((float) cparams->nz / tpb.z));
				per_xy_edges<<<bpg, tpb, 0, stream>>>(*d_grid);
				CUDA_ERRCHK_KERNEL();
			}
			//Copy periodic xz edges
			if (BOUNDCOND_TYPE_Z == PERIODIC_BOUNDCONDS) {
                const dim3 tpb(3, 32, 3);
                const dim3 bpg(
                           1, 
                           (unsigned int)ceil((real) cparams->ny / (real)tpb.y), 
                           4);
				per_xz_edges<<<bpg, tpb, 0, stream>>>(*d_grid);
				CUDA_ERRCHK_KERNEL();
			}
			//If fully periodic, copy all corners
			if ((BOUNDCOND_TYPE_Y == PERIODIC_BOUNDCONDS) && (BOUNDCOND_TYPE_Z == PERIODIC_BOUNDCONDS)) {	
                const dim3 tpb(3, 3, 3);
                const dim3 bpg(1, 1, 8);
				per_xyz_corners<<<bpg, tpb, 0, stream>>>(*d_grid);
				CUDA_ERRCHK_KERNEL();
			}
			break;
        }
		default:
			printf("INVALID X TYPE IN BOUNDCOND_CUDA!\n");
			exit(EXIT_FAILURE);
	}
	//--------------------------------

	//--------Y BOUNDS--------------
	switch	(BOUNDCOND_TYPE_Y) {

		//Do periodic bounds for y sides
		case PERIODIC_BOUNDCONDS: {
            const dim3 tpb(32, 32, 1);
            const dim3 bpg(
                       (unsigned int)ceil((real) cparams->nx / (real)tpb.x),
                       (unsigned int)ceil((real) cparams->nz / (real)tpb.y),
                       6);
			per_y_sides<<<bpg, tpb, 0, stream>>>(*d_grid);
			CUDA_ERRCHK_KERNEL();
	
			//Copy periodic yz edges
			if (BOUNDCOND_TYPE_Z == PERIODIC_BOUNDCONDS) {
                const dim3 tpb(32, 3, 3);
                const dim3 bpg(
                           (unsigned int)ceil((real) cparams->nx / (real)tpb.x),
                           1,
                           4);
				per_yz_edges<<<bpg, tpb, 0, stream>>>(*d_grid);
				CUDA_ERRCHK_KERNEL();

			}
			break;
        }
		default:
			printf("INVALID Y TYPE IN BOUNDCOND_CUDA!\n");
			exit(EXIT_FAILURE);
	}
	//--------------------------------


	//---------Z BOUNDS----------------
	switch	(BOUNDCOND_TYPE_Z) {

		//Do periodic bounds for z sides
		case PERIODIC_BOUNDCONDS: {
            const dim3 tpb(32, 32, 1);
            const dim3 bpg((unsigned int)ceil((real) cparams->nx / (real)tpb.x),
                           (unsigned int)ceil((real) cparams->ny / (real)tpb.y),
                           6);
			per_z_sides<<<bpg, tpb, 0, stream>>>(*d_grid);
			CUDA_ERRCHK_KERNEL();
			break;
        }
		default:
			printf("INVALID Z TYPE IN BOUNDCOND_CUDA!\n");
			exit(EXIT_FAILURE);
	}	
	//--------------------------------
}


void periodic_xy_boundconds_cuda_generic(Grid* d_grid, CParamConfig* cparams, cudaStream_t stream)
{
    //Copy periodic x sides
    {
        const dim3 tpb(6, 4, 1);
        const dim3 bpg(1,
                       (unsigned int)ceil((real) cparams->ny / tpb.y),
                       (unsigned int)ceil((real) cparams->nz / tpb.z));
        per_x_sides<<<bpg, tpb, 0, stream>>>(*d_grid);
        CUDA_ERRCHK_KERNEL();
    }

    //Copy periodic xy edges
    {
        const dim3 tpb(3, 3, 32);
        const dim3 bpg(1,
                       4,
                       (unsigned int)ceil((real) cparams->nz / tpb.z));
	    per_xy_edges<<<bpg, tpb, 0, stream>>>(*d_grid);
	    CUDA_ERRCHK_KERNEL();
    }

    //Copy periodic y sides
    {
        const dim3 tpb(32, 32, 1);
        const dim3 bpg((unsigned int)ceil((real) cparams->nx / tpb.x),
                       (unsigned int)ceil((real) cparams->nz / tpb.y),
                       6);
        per_y_sides<<<bpg, tpb, 0, stream>>>(*d_grid);
        CUDA_ERRCHK_KERNEL();
    }
}
