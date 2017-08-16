#include "rk3_cuda_generic.cuh"
#include "diff_cuda_generic.cuh"

#include "gpu/cuda/core/errorhandler_cuda.cuh"
#include "common/errorhandler.h"

static real *d_w_lnrho, *d_w_uu_x, *d_w_uu_y, *d_w_uu_z;
static bool is_initialized      = false;
static CParamConfig* cparams    = NULL;


void init_rk3_cuda_generic(CParamConfig* conf)
{
    if (is_initialized) CRASH("rk3_cuda_generic() already initialized!");
    is_initialized = true;
    cparams = conf;

    const int w_grid_size = cparams->nx * cparams->ny * cparams->nz;

	CUDA_ERRCHK( cudaMalloc(&d_w_lnrho, sizeof(real)*w_grid_size) );
	CUDA_ERRCHK( cudaMalloc(&d_w_uu_x,  sizeof(real)*w_grid_size) );
	CUDA_ERRCHK( cudaMalloc(&d_w_uu_y,  sizeof(real)*w_grid_size) );
	CUDA_ERRCHK( cudaMalloc(&d_w_uu_z,  sizeof(real)*w_grid_size) );    
}


void destroy_rk3_cuda_generic()
{
    if (!is_initialized) CRASH("rk3_cuda_generic() wasn't initialized!");
    is_initialized = false;
    cparams = NULL;

    CUDA_ERRCHK( cudaFree(d_w_lnrho) );
    CUDA_ERRCHK( cudaFree(d_w_uu_x)  );
    CUDA_ERRCHK( cudaFree(d_w_uu_y)  );
    CUDA_ERRCHK( cudaFree(d_w_uu_z)  );
}


//template <int step_number>
//__launch_bounds__(RK_THREADBLOCK_SIZE, 1)
static __global__ void rk3_step(const real* __restrict__ d_lnrho,  //SOURCE
                        const real* __restrict__ d_uu_x, 
                        const real* __restrict__ d_uu_y, 
                        const real* __restrict__ d_uu_z, 
                        real* __restrict__ d_w_lnrho,       //INTERMEDIATE
                        real* __restrict__ d_w_uu_x, 
                        real* __restrict__ d_w_uu_y, 
                        real* __restrict__ d_w_uu_z,
                  		real* __restrict__ d_lnrho_dst,     //DESTINATION
                        real* __restrict__ d_uu_x_dst, 
                        real* __restrict__ d_uu_y_dst, 
                        real* __restrict__ d_uu_z_dst,
                        const int step_number, const real dt)
{
    const real alphas[] = {0.0, -0.53125, -1.1851851851851851};
    const real betas[]  = {0.25, 0.88888888888888884, 0.75};
    const real ALPHA = alphas[step_number];
    const real BETA  = betas[step_number];

	//Allocate shared memory
	__shared__ real s_lnrho[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH];
	__shared__ real s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH];
	__shared__ real s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH];
	__shared__ real s_uu_z[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH];

	real r_w_lnrho;
	real r_w_uu_x;
	real r_w_uu_y;
	real r_w_uu_z;

	
	//Define shared memory indices for the halo containing array:

	//Intermediate result array (wsid)
	int wsid_col = threadIdx.x;
	int wsid_row = threadIdx.y;
	int wsid_depth = threadIdx.z;

	//Shared memory array indices 
	//(these are offset to start from first non-boundary value, eg.
	//s_lnrho[sid_col][3][3] and r_w_lnrho[wsid_col][0][0] represent 
	//the same grid point's final and partial results.)
	int sid_col = wsid_col + BOUND_SIZE; 
	int sid_row = wsid_row + BOUND_SIZE; 
	int sid_depth = wsid_depth + BOUND_SIZE;

	//Define grid indices (in global memory) *without* pad/boundary offset
	int grid_idx_x = threadIdx.x + blockIdx.x*blockDim.x;
	int grid_idx_y = threadIdx.y + blockIdx.y*blockDim.y;
	int grid_idx_z = threadIdx.z + blockIdx.z*blockDim.z*RK_ELEMS_PER_THREAD;
	
	//Load constants from constant memory to registers
	int w_grid_y_offset = d_nx;
	int w_grid_z_offset = d_nx*d_ny;

	int grid_y_offset = d_mx;
	int grid_z_offset = d_mx*d_my;

	//Index in the partial result array (doesn't include boundary zones)
	int w_grid_idx = (grid_idx_x) +
			 (grid_idx_y)*w_grid_y_offset +
			 (grid_idx_z)*w_grid_z_offset;

	//Index in the final result array 
	//(is offset to start from first index of the computational domain)
	int grid_idx = 	(grid_idx_x + d_nx_min) +
			(grid_idx_y + d_ny_min)*grid_y_offset +
			(grid_idx_z + d_nz_min)*grid_z_offset;


///////
	//Initialization done, start loading stuff to faster memory
///////

	//Load intermediate results to registers
	switch (step_number) {
		case 0:
			r_w_lnrho = 0.0;
			r_w_uu_x  = 0.0;
			r_w_uu_y  = 0.0;
			r_w_uu_z  = 0.0;
			break;
		default:
			r_w_lnrho = d_w_lnrho[w_grid_idx];
			r_w_uu_x  = d_w_uu_x [w_grid_idx];
			r_w_uu_y  = d_w_uu_y [w_grid_idx];
			r_w_uu_z  = d_w_uu_z [w_grid_idx];
			break;
	}

	//Load the whole smem block from global memory during the first iteration
	int cpy_idx;
   	for (int depth=0; depth < SHARED_SIZE_DEPTH; depth += SHARED_SIZE_W_DEPTH) {
		if (wsid_depth+depth >= SHARED_SIZE_DEPTH) break; //overflow
      		for (int row=0; row < SHARED_SIZE_ROW; row += SHARED_SIZE_W_ROW) {
			if (wsid_row+row >= SHARED_SIZE_ROW) break; //overflow
        		for (int col=0; col < SHARED_SIZE_COL; col += SHARED_SIZE_W_COL) {	
				if (wsid_col+col >= SHARED_SIZE_COL) break; //overflow
            			//Shift index depending on the quater
				//MV: Seems to be formally ok. 
            			cpy_idx = 	(grid_idx_x + d_nx_min-BOUND_SIZE+col) + 
						(grid_idx_y + d_ny_min-BOUND_SIZE+row)*grid_y_offset + 
						(grid_idx_z + d_nz_min-BOUND_SIZE+depth)*grid_z_offset;
           			//Copy memory
        	    		s_lnrho[wsid_row+row][wsid_col+col][wsid_depth+depth] = d_lnrho[cpy_idx];
            			s_uu_x[wsid_row+row][wsid_col+col][wsid_depth+depth] = d_uu_x[cpy_idx];
           			s_uu_y[wsid_row+row][wsid_col+col][wsid_depth+depth] = d_uu_y[cpy_idx];
            			s_uu_z[wsid_row+row][wsid_col+col][wsid_depth+depth] = d_uu_z[cpy_idx];
			}
		}
	}
	__syncthreads();
	
	//Continuity result
	real cont_res;  
	//Momentum
	real mom_x, mom_y, mom_z;
	//Shear 
	//real u_shear_y = 0.0;
	
	int curr_iteration = 0;
	//Compute
	while (true) {

		//Get the shearing flow  
		//if (d_LSHEAR) {
			//shear_flow(&u_shear_y, grid_idx_x);
		//} 


		//Solve derivatives
		real ddx_lnrho = der_scalx(sid_row, sid_col, sid_depth, s_lnrho);
		real ddx_uu_x = der_scalx(sid_row, sid_col, sid_depth, s_uu_x);
		real ddx_uu_y = der_scalx(sid_row, sid_col, sid_depth, s_uu_y);
		real ddx_uu_z = der_scalx(sid_row, sid_col, sid_depth, s_uu_z);

		real ddy_lnrho = der_scaly(sid_row, sid_col, sid_depth, s_lnrho);
		real ddy_uu_x = der_scaly(sid_row, sid_col, sid_depth, s_uu_x);
		real ddy_uu_y = der_scaly(sid_row, sid_col, sid_depth, s_uu_y);
		real ddy_uu_z = der_scaly(sid_row, sid_col, sid_depth, s_uu_z);	

		real ddz_lnrho = der_scalz(sid_row, sid_col, sid_depth, s_lnrho);
		real ddz_uu_x = der_scalz(sid_row, sid_col, sid_depth, s_uu_x);
		real ddz_uu_y = der_scalz(sid_row, sid_col, sid_depth, s_uu_y);
		real ddz_uu_z = der_scalz(sid_row, sid_col, sid_depth, s_uu_z);	

		//Continuity equation
		//cont_res = continuity(sid_row, sid_col, sid_depth, u_shear_y, s_lnrho, s_uu_x, s_uu_y, s_uu_z);
		cont_res = -(s_uu_x[sid_row][sid_col][sid_depth] * ddx_lnrho +
				s_uu_y[sid_row][sid_col][sid_depth] * ddy_lnrho +
				s_uu_z[sid_row][sid_col][sid_depth] * ddz_lnrho) 
				-(ddx_uu_x + ddy_uu_y + ddz_uu_z);

		r_w_lnrho = ALPHA*r_w_lnrho + dt*cont_res;

		//Navier-Stokes
		//navier_stokes(&mom_x, &mom_y, &mom_z, sid_row, sid_col, sid_depth, u_shear_y, s_lnrho, s_uu_x, s_uu_y, s_uu_z);
		real d2dx2_uu_x = der2_scalx(sid_row, sid_col, sid_depth, s_uu_x);
		real d2dy2_uu_x = der2_scaly(sid_row, sid_col, sid_depth, s_uu_x);
		real d2dz2_uu_x = der2_scalz(sid_row, sid_col, sid_depth, s_uu_x);

		real d2dx2_uu_y = der2_scalx(sid_row, sid_col, sid_depth, s_uu_y);
		real d2dy2_uu_y = der2_scaly(sid_row, sid_col, sid_depth, s_uu_y);
		real d2dz2_uu_y = der2_scalz(sid_row, sid_col, sid_depth, s_uu_y);

		real d2dx2_uu_z = der2_scalx(sid_row, sid_col, sid_depth, s_uu_z);
		real d2dy2_uu_z = der2_scaly(sid_row, sid_col, sid_depth, s_uu_z);
		real d2dz2_uu_z = der2_scalz(sid_row, sid_col, sid_depth, s_uu_z);

		mom_x = - (s_uu_x[sid_row][sid_col][sid_depth] * ddx_uu_x + //vec_dot_nabla_scal
			s_uu_y[sid_row][sid_col][sid_depth] * ddy_uu_x +
			s_uu_z[sid_row][sid_col][sid_depth] * ddz_uu_x)
			- d_CS2_SOUND*ddx_lnrho //grad lnrho
			+ d_NU_VISC * ( d2dx2_uu_x + d2dy2_uu_x + d2dz2_uu_x); //nu_const 

		mom_y = - (s_uu_x[sid_row][sid_col][sid_depth] * ddx_uu_y +
			s_uu_y[sid_row][sid_col][sid_depth] * ddy_uu_y +
			s_uu_z[sid_row][sid_col][sid_depth] * ddz_uu_y)
			- d_CS2_SOUND*ddy_lnrho
			+ d_NU_VISC * ( d2dx2_uu_y + d2dy2_uu_y + d2dz2_uu_y );

		mom_z = - (s_uu_x[sid_row][sid_col][sid_depth] * ddx_uu_z +
			s_uu_y[sid_row][sid_col][sid_depth] * ddy_uu_z +
			s_uu_z[sid_row][sid_col][sid_depth] * ddz_uu_z)
			- d_CS2_SOUND*ddz_lnrho
			+ d_NU_VISC * ( d2dx2_uu_z + d2dy2_uu_z + d2dz2_uu_z );
		
		//grad_div_vec       
		mom_x += d_NU_VISC*(1.0/3.0)*((d2dx2_uu_x + dernm_scal(sid_row, sid_col, sid_depth, 1, 2, s_uu_y) + dernm_scal(sid_row, sid_col, sid_depth, 1, 3, s_uu_z)));
		mom_y += d_NU_VISC*(1.0/3.0)*((dernm_scal(sid_row, sid_col, sid_depth, 2, 1, s_uu_x) + d2dy2_uu_y + dernm_scal(sid_row, sid_col, sid_depth, 2, 3, s_uu_z)));
		mom_z += d_NU_VISC*(1.0/3.0)*((dernm_scal(sid_row, sid_col, sid_depth, 3, 1, s_uu_x) + dernm_scal(sid_row, sid_col, sid_depth, 3, 2, s_uu_y) + d2dz2_uu_z));
		
		//S_grad_lnrho 
		real Sxx, Sxy, Sxz, Syy, Syz, Szz;
		Sxx = (2.0/3.0)*ddx_uu_x - (1.0/3.0)*(ddy_uu_y + ddz_uu_z);
		Sxy = 0.5*(ddy_uu_x + ddx_uu_y);
		Sxz = 0.5*(ddz_uu_x + ddx_uu_z);

		Syy = (2.0/3.0)*ddy_uu_y - (1.0/3.0)*(ddx_uu_x + ddz_uu_z);
		Syz = 0.5*(ddz_uu_y + ddy_uu_z);
		Szz = (2.0/3.0)*ddz_uu_z - (1.0/3.0)*(ddx_uu_x + ddy_uu_y);//Note, old CPU version was bugged here (ddz_uu_x instead of dd_z_uu_z)

		mom_x += 2.0*d_NU_VISC*(Sxx*ddx_lnrho + Sxy*ddy_lnrho + Sxz*ddz_lnrho);
		mom_y += 2.0*d_NU_VISC*(Sxy*ddx_lnrho + Syy*ddy_lnrho + Syz*ddz_lnrho);
		mom_z += 2.0*d_NU_VISC*(Sxz*ddx_lnrho + Syz*ddy_lnrho + Szz*ddz_lnrho);


		r_w_uu_x = ALPHA*r_w_uu_x + dt*mom_x;
		r_w_uu_y = ALPHA*r_w_uu_y + dt*mom_y;
		r_w_uu_z = ALPHA*r_w_uu_z + dt*mom_z;
		__syncthreads();

    
        //Write the result to global memory
		d_lnrho_dst[grid_idx] = s_lnrho[sid_row][sid_col][sid_depth] + BETA*r_w_lnrho;
		d_uu_x_dst[grid_idx] = s_uu_x[sid_row][sid_col][sid_depth] + BETA*r_w_uu_x;
		d_uu_y_dst[grid_idx] = s_uu_y[sid_row][sid_col][sid_depth] + BETA*r_w_uu_y;
		d_uu_z_dst[grid_idx] = s_uu_z[sid_row][sid_col][sid_depth] + BETA*r_w_uu_z;	

		d_w_lnrho[w_grid_idx] = r_w_lnrho;
		d_w_uu_x [w_grid_idx] = r_w_uu_x;
		d_w_uu_y [w_grid_idx] = r_w_uu_y;
		d_w_uu_z [w_grid_idx] = r_w_uu_z;

		curr_iteration++;
		if (curr_iteration < RK_ELEMS_PER_THREAD) {

			//Advance grid index by a block depth to z direction
			grid_idx_z += RK_THREADS_Z;
			grid_idx += RK_THREADS_Z*grid_z_offset;
			w_grid_idx += RK_THREADS_Z*w_grid_z_offset;

			int prev_smem_start_idx = WRAP_SMEM_DEPTH((curr_iteration-1)*RK_THREADS_Z);
			int curr_smem_start_idx = WRAP_SMEM_DEPTH(curr_iteration*RK_THREADS_Z);
			
			//Copy new data from global memory (todo smem->smem)
   			for (int depth=0; depth < SHARED_SIZE_DEPTH; depth += SHARED_SIZE_W_DEPTH) {

			int actual_depth = WRAP_SMEM_DEPTH(wsid_depth+depth + prev_smem_start_idx);
			if (prev_smem_start_idx > curr_smem_start_idx) {
				if (actual_depth >= curr_smem_start_idx && actual_depth < prev_smem_start_idx) {
					break;
				}
			}
			else { 
				if (actual_depth < prev_smem_start_idx || actual_depth >= curr_smem_start_idx) {
					break;
				}
			}

			int grid_z_update_offset;
			int normalized_smem_offset;
			if (actual_depth >= prev_smem_start_idx) {
				normalized_smem_offset = actual_depth - prev_smem_start_idx;
			}
			else {
				normalized_smem_offset = SHARED_SIZE_DEPTH - prev_smem_start_idx + actual_depth;
			}			
			grid_z_update_offset = grid_idx_z - threadIdx.z + 2*BOUND_SIZE + normalized_smem_offset;
			
			if (wsid_depth+depth >= SHARED_SIZE_DEPTH) break; //overflow
      				for (int row=0; row < SHARED_SIZE_ROW; row += SHARED_SIZE_W_ROW) {
				if (wsid_row+row >= SHARED_SIZE_ROW) break; //overflow
        				for (int col=0; col < SHARED_SIZE_COL; col += SHARED_SIZE_W_COL) {	
					if (wsid_col+col >= SHARED_SIZE_COL) break; //overflow
            					//Shift index depending on the quater
						//MV: Seems to be formally ok. 
            					cpy_idx = 	(grid_idx_x + d_nx_min-BOUND_SIZE+col) + 
								(grid_idx_y + d_ny_min-BOUND_SIZE+row)*grid_y_offset +
								(d_nz_min-BOUND_SIZE+grid_z_update_offset)*grid_z_offset;

           					 //Copy memory
            					s_lnrho[wsid_row+row][wsid_col+col][actual_depth] = d_lnrho[cpy_idx];
            					s_uu_x[wsid_row+row][wsid_col+col][actual_depth] = d_uu_x[cpy_idx];
           					s_uu_y[wsid_row+row][wsid_col+col][actual_depth] = d_uu_y[cpy_idx];
            					s_uu_z[wsid_row+row][wsid_col+col][actual_depth] = d_uu_z[cpy_idx];
					
         				}
      				}
   			}
		
			//Load intermediate results to registers
			switch (step_number) {
				case 0:
					r_w_lnrho = 0.0;
					r_w_uu_x  = 0.0;
					r_w_uu_y  = 0.0;
					r_w_uu_z  = 0.0;
					break;
				default:
					r_w_lnrho = d_w_lnrho[w_grid_idx];
					r_w_uu_x  = d_w_uu_x [w_grid_idx];
					r_w_uu_y  = d_w_uu_y [w_grid_idx];
					r_w_uu_z  = d_w_uu_z [w_grid_idx];
					break;
			}
			__syncthreads();

			sid_depth = WRAP_SMEM_DEPTH(sid_depth + RK_THREADS_Z);  

		}
		else {
			//We're done, break
			break;
		}

	}    
}


void rk3_cuda_generic(  real* d_lnrho,      real* d_uu_x,       real* d_uu_y,       real* d_uu_z, 
                        real* d_lnrho_dst,  real* d_uu_x_dst,   real* d_uu_y_dst,   real* d_uu_z_dst,
                        const int step_number, const real dt)
{
    if (!is_initialized) CRASH("rk3_cuda_generic() wasn't initialized!");

    const dim3 threads_per_block = {RK_THREADS_X, RK_THREADS_Y, RK_THREADS_Z};
    const dim3 blocks_per_grid   = {(unsigned int) ceil((real) cparams->nx / (real)threads_per_block.x),
                                    (unsigned int) ceil((real) cparams->ny / (real)threads_per_block.y),
                                    (unsigned int) ceil((real) cparams->nz / (real)(threads_per_block.z*RK_ELEMS_PER_THREAD))};

    //TIME START
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

    //INTEGRATE
    rk3_step<<<blocks_per_grid, threads_per_block>>>(  d_lnrho, d_uu_x, d_uu_y, d_uu_z, 
                                                       d_w_lnrho, d_w_uu_x, d_w_uu_y, d_w_uu_z, 
                                                       d_lnrho_dst, d_uu_x_dst, d_uu_y_dst, d_uu_z_dst,  
                                                       step_number, dt); 
    CUDA_ERRCHK_KERNEL();

	//TIME END
	cudaDeviceSynchronize();
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	//printf("RK integration step time elapsed: \t%f ms\n", time);
}


//Explicit instantiation for the integration template
/*
template void rk3_cuda_generic<0>(real* d_lnrho,      real* d_uu_x,       real* d_uu_y,       real* d_uu_z, 
                                  real* d_lnrho_dst,  real* d_uu_x_dst,   real* d_uu_y_dst,   real* d_uu_z_dst,
                                  const real dt);
template void rk3_cuda_generic<1>(real* d_lnrho,      real* d_uu_x,       real* d_uu_y,       real* d_uu_z, 
                                  real* d_lnrho_dst,  real* d_uu_x_dst,   real* d_uu_y_dst,   real* d_uu_z_dst,
                                  const real dt);
template void rk3_cuda_generic<2>(real* d_lnrho,      real* d_uu_x,       real* d_uu_y,       real* d_uu_z, 
                                  real* d_lnrho_dst,  real* d_uu_x_dst,   real* d_uu_y_dst,   real* d_uu_z_dst,
                                  const real dt);
*/




































































