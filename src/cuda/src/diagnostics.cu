
#include <stdio.h>
#include <float.h>

#define EXTERN extern
#include "dconsts.cuh"
#include "../cparam_c.h"
#include "smem.cuh"

/****************************************************************************************************************/
__device__ void check_for_nan_inf_variable(int code, float var)
{
	// Check for the presence of inf or nan in the shared memory in a variable

	if (isnan(var)) {
		if (code == 1) {
			printf("\n   ALERT! found NaN in thread (%i, %i, %i) in block (%i, %i, %i) in u_dot_grad_lnrho! \n\n",
				threadIdx.x, threadIdx.y, threadIdx.z, blockIdx.x, blockIdx.y, blockIdx.z);
		} else if (code == 2) {
                        printf("\n   ALERT! found NaN in thread (%i, %i, %i) in block (%i, %i, %i) in div_u! \n\n",
                                threadIdx.x, threadIdx.y, threadIdx.z, blockIdx.x, blockIdx.y, blockIdx.z);
                } else if (code == 3) {
                        printf("\n   ALERT! found NaN in thread (%i, %i, %i) in block (%i, %i, %i) in per_x_sides()! \n\n",
                                threadIdx.x, threadIdx.y, threadIdx.z, blockIdx.x, blockIdx.y, blockIdx.z);
                } else if (code == 4) {
                        printf("\n   ALERT! found NaN in thread (%i, %i, %i) in block (%i, %i, %i) in per_y_sides()! \n\n",
                                threadIdx.x, threadIdx.y, threadIdx.z, blockIdx.x, blockIdx.y, blockIdx.z);
		} else {
                        printf("\n   ALERT! found NaN in thread (%i, %i, %i) in block (%i, %i, %i) in NOT DEFINED!!! \n\n",
                                threadIdx.x, threadIdx.y, threadIdx.z, blockIdx.x, blockIdx.y, blockIdx.z);
		}
	}
}
/****************************************************************************************************************/
__device__ void check_for_nan_inf(int step_number, int grid_idx_x, int grid_idx_y, int grid_idx_z,
                                  int sid_row, int sid_col, int sid_depth,
                                  float s_lnrho[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH],
                                  float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH],
                                  float s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH],
                                  float s_uu_z[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH] )
{
        // Check for the presence of inf or nan in the shared memory. This might be useful for both stability 
        // cheking and debugging purposes.
        int nan_count = 0;
        int inf_count = 0;

        //Check and count for NaN
        nan_count += isnan( s_lnrho[sid_row][sid_col][sid_depth] );
        nan_count += isnan( s_uu_x[sid_row][sid_col][sid_depth] );
        nan_count += isnan( s_uu_y[sid_row][sid_col][sid_depth] );
        nan_count += isnan( s_uu_z[sid_row][sid_col][sid_depth] );

        //Check and count for inf
        inf_count += isinf( s_lnrho[sid_row][sid_col][sid_depth] );
        inf_count += isinf( s_uu_x[sid_row][sid_col][sid_depth] );
        inf_count += isinf( s_uu_y[sid_row][sid_col][sid_depth] );
        inf_count += isinf( s_uu_z[sid_row][sid_col][sid_depth] );

        if (nan_count > 0) {
                printf("\n   ALERT! found %i NaN in thread (%i, %i, %i) in block (%i, %i, %i) in the coordinate (x,y,z) = (%i, %i, %i) RK_step: %i \n\n",
                        nan_count, threadIdx.x, threadIdx.y, threadIdx.z, blockIdx.x, blockIdx.y, blockIdx.z, grid_idx_x, grid_idx_y, grid_idx_z, step_number);
        }

        if (inf_count > 0) {
                printf("\n   ALERT! found %i Inf in thread (%i, %i, %i) in block (%i, %i, %i) in the coordinate (x,y,z) = (%i, %i, %i)\n\n",
                        inf_count, threadIdx.x, threadIdx.y, threadIdx.z, blockIdx.x, blockIdx.y, blockIdx.z, grid_idx_x, grid_idx_y, grid_idx_z);
        }

}
/****************************************************************************************************************/
__global__ void check_grid_for_nan_cuda(float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z, int* d_nan_count) 
{
	//Look into device memory and locate all NaN

	//Define the grid coordinates 
	int grid_idx_x = threadIdx.x + blockIdx.x*blockDim.x;
	int grid_idx_y = threadIdx.y + blockIdx.y*blockDim.y;
	int grid_idx_z = threadIdx.z + blockIdx.z*blockDim.z;

	//The actual index in the array 
	//%JP: Check only the computational domain starting from d_C?_BOT
	int grid_idx = (grid_idx_x + d_CX_BOT) +
		       (grid_idx_y + d_CY_BOT)*d_NX +
		       (grid_idx_z + d_CZ_BOT)*d_NX*d_NY;

	//Detect overflow and skip those (this could be more smooth)
	//%JP: Checked for overflow too late (after trying to load dirty data), fixed
	if ((grid_idx_x < d_CX_TOP) && (grid_idx_y < d_CY_TOP) && (grid_idx_z < d_CZ_TOP)) {

		//Check for NaN in the define index
		int nan_tot, nan_lnrho, nan_uu_x, nan_uu_y, nan_uu_z;
		nan_lnrho = isnan(d_lnrho[grid_idx]);
		nan_uu_x = isnan(d_uu_x[grid_idx]);
		nan_uu_y = isnan(d_uu_y[grid_idx]);
		nan_uu_z = isnan(d_uu_z[grid_idx]);
		nan_tot = nan_lnrho + nan_uu_x + nan_uu_y + nan_uu_x;

		//Print alert if NaN found
		if (nan_tot > 0) {
			//Add to nan count. NOTE: Precision is not that important here. Hope this is enough. 
			*d_nan_count += nan_tot;
			printf("check_grid_for_nan: WARNING! %i NaN found in: %i d_lnrho %i d_uu_x %i d_uu_y %i d_uu_z, at (%i, %i, %i) \n",
				nan_tot, nan_lnrho, nan_uu_x, nan_uu_y, nan_uu_z, grid_idx_x+d_CX_BOT, grid_idx_y+d_CY_BOT, grid_idx_z+d_CZ_BOT);
		}
	}
	
}
/****************************************************************************************************************/
/*
float check_grid_for_nan(float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z)
{
        //Look into device memory and locate all NaN. Set here the threads and blocks and call the kernel itself. 

	//Determine threadblock dims
	static dim3 blocksPerGrid, threadsPerBlock;
	threadsPerBlock.x = RK_THREADS_X;
	threadsPerBlock.y = RK_THREADS_Y;
	threadsPerBlock.z = RK_THREADS_Z;
	//%JP: Check only the computational domain. We have a padding in the x-axis
	//between x=0...CX_BOT-BOUND_SIZE, which may have whatever values (including NaNs)
	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)threadsPerBlock.x);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)threadsPerBlock.y);
	blocksPerGrid.z = ceil((float) COMP_DOMAIN_SIZE_Z / (float)threadsPerBlock.z);

	//TODO: Check possible overflow. (Also in the integrators) 
	printf("\nIn check_grid_for_nan \n NX = %i, NY = %i, NZ = %i \n threadsPerBlock = (%i, %i, %i) blocksPerGrid = (%i, %i, %i) ... ... ... ",
		NX, NY, NZ, threadsPerBlock.x, threadsPerBlock.y, threadsPerBlock.z, blocksPerGrid.x, blocksPerGrid.y, blocksPerGrid.z);

	int found_nan = 0; 
	int nan_count = 0;
	int* d_nan_count;

	//Initialize the diagnostic counter variable
	cudaMalloc((int**) &d_nan_count, sizeof(int)); 
	cudaMemcpy((int*) d_nan_count, &nan_count, sizeof(int), cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();

	//Call the kernel 
	check_grid_for_nan_cuda<<<blocksPerGrid, threadsPerBlock>>>(d_lnrho, d_uu_x, d_uu_y, d_uu_z, d_nan_count);

	//Syncronize before continuing
	cudaDeviceSynchronize();
	cudaMemcpy(&nan_count, (int*) d_nan_count, sizeof(int), cudaMemcpyDeviceToHost); 
	cudaDeviceSynchronize();
	if (nan_count > 0) {
		printf("%i NaN detected! Terminating soon...\n\n", nan_count);
		found_nan = 1; 
	} else {
		//printf("No NaN detected. \n\n", nan_count);
		printf("No NaN detected. \n\n");
	}  

	return found_nan;
}*/
/****************************************************************************************************************/

#include "../cdata_c.h"
#include "../viscosity_c.h"
#include "../eos_c.h"
#include "../forcing_c.h"
#include "defines_PC.h"

/****************************************************************************************************************/
void print_init_config()
{
   	printf("Comp domain sizes: %d, %d, %d\n", COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z);
	printf("Actual grid sizes (with pads+bounds): %d, %d, %d\n", NX, NY, NZ);
   	printf("Lbox_x = %f, Lbox_y = %f, Lbox_z = %f \n", DOMAIN_SIZE_X, DOMAIN_SIZE_Y, DOMAIN_SIZE_Z);
/*   	printf("iboundx = %d, iboundy = %d, iboundz = %d \n", BOUNDCOND_TYPE_X, BOUNDCOND_TYPE_Y, BOUNDCOND_TYPE_Z);
   	printf("unit_length = %f \n", UNIT_LENGTH);
*/
}
/****************************************************************************************************************/
void print_run_config()
{
	printf("NU_VISC %f\n", NU_VISC);
	printf("CS_SOUND %f\n", CS_SOUND);

	printf("LFORCING %d\n", LFORCING);
	printf("FORCING %f\n", FORCING);
	printf("KK1 %f\n", KK1);
	printf("KK2 %f\n", KK2);
	printf("KMAX %f\n", KMAX);
	printf("DKX %f\n", DKX);
	printf("DKY %f\n", DKY);
	printf("DKZ %f\n", DKZ);

	printf("Q_SHEAR %f\n", Q_SHEAR);
	printf("OMEGA %f\n", OMEGA);
}
/****************************************************************************************************************/
void print_additional_defines()
{
	printf("bound size %d\n", BOUND_SIZE);
	//printf(" %d", NDIM 3 //useless?

	printf("Comp domain top indices xyz:\n%d\n", CX_TOP);
	printf("%d\n", CY_TOP);
	printf("%d\n", CZ_TOP);

	printf("Comp domain bot indices xyz:\n%d\n", CX_BOT);
	printf("%d\n", CY_BOT);
	printf("%d\n", CZ_BOT);

	printf("Distances between gridpoints:\n%f, ", DX);
	printf("%f, ", DY);
	printf("%f\n", DZ);
}
/****************************************************************************************************************/
void print_grid_data(float* grid)
{
	for(int k=0; k < NZ; k++) {
		for(int j=0; j < NY; j++) {
			for (int i=0; i < NX; i++) {
				printf("%f ", grid[i + j*NX + k*NX*NY]);
			}
			printf("\n");
		}
		printf("\n\n\n");
	}
}
/****************************************************************************************************************/
float check_grids(float* CPU_lnrho, float* CPU_uu_x, float* CPU_uu_y, float* CPU_uu_z,
    		  float* GPU_lnrho, float* GPU_uu_x, float* GPU_uu_y, float* GPU_uu_z) 
{
	float error = 0.0f;

	float lnrho_error = 0.0f;
	float uu_x_error = 0.0f;
	float uu_y_error = 0.0f;
	float uu_z_error = 0.0f;

	for (int k=CZ_BOT; k < CZ_TOP; k++) {
		for (int j=CY_BOT; j < CY_TOP; j++) {
			for (int i=CX_BOT; i < CX_TOP; i++) {
				int idx = i + j*NX + k*NX*NY;
				if (	isnan(GPU_lnrho[idx]) || 
					isnan(GPU_uu_x[idx]) || 	
					isnan(GPU_uu_x[idx]) || 
					isnan(GPU_uu_x[idx])) 
				{
					error = FLT_MAX;
					printf("GPU result contains nan!\n");
					return error;
				}
				
				if (lnrho_error < abs(CPU_lnrho[idx]-GPU_lnrho[idx]))
					lnrho_error = abs(CPU_lnrho[idx]-GPU_lnrho[idx]);
				if (uu_x_error < abs(CPU_uu_x[idx]-GPU_uu_x[idx]))
					uu_x_error = abs(CPU_uu_x[idx]-GPU_uu_x[idx]);
				if (uu_y_error < abs(CPU_uu_y[idx]-GPU_uu_y[idx]))
					uu_y_error = abs(CPU_uu_y[idx]-GPU_uu_y[idx]);
				if (uu_z_error < abs(CPU_uu_z[idx]-GPU_uu_z[idx]))
					uu_z_error = abs(CPU_uu_z[idx]-GPU_uu_z[idx]);				

				error += abs(CPU_lnrho[idx]-GPU_lnrho[idx]);
				error += abs(CPU_uu_x[idx]-GPU_uu_x[idx]);
				error += abs(CPU_uu_y[idx]-GPU_uu_y[idx]);
				error += abs(CPU_uu_z[idx]-GPU_uu_z[idx]);
			}
		}
	}	
	
	if (lnrho_error > 0.0f)
		printf("\n\t\tHighest lnrho_error: %f", lnrho_error);
	if (uu_x_error > 0.0f)
		printf("\n\t\tHighest uu_x_error: %f", uu_x_error);
	if (uu_y_error > 0.0f)	
		printf("\n\t\tHighest uu_y_error: %f", uu_y_error);
	if (uu_z_error > 0.0f)	
		printf("\n\t\tHighest uu_z_error: %f", uu_z_error);
	printf("\n\t\tCPU / GPU: %f, %f\n", CPU_uu_y[CX_BOT + CY_BOT*NX + CZ_BOT*NX*NY], GPU_uu_y[CX_BOT + CY_BOT*NX + CZ_BOT*NX*NY]);
	return error;
}
/****************************************************************************************************************/

typedef float (*VecReductionFunctionHostPointer)(float* vec_x, float* vec_y, float* vec_z);
typedef void (*VecReductionFunctionDevicePointer)(float* d_vec_max, float* d_partial_result, 
						float* d_vec_x, float* d_vec_y, float* d_vec_z);

typedef float (*ScalReductionFunctionHostPointer)(float* scal);
typedef void (*ScalReductionFunctionDevicePointer)(float* d_vec_max, float* d_partial_result, 
						float* d_scal);
/*
* Runs diagnostics by checking that both CPU and GPU versions of 
* certain functions return the same result.
* Note: Requires that the GPU memory is already allocated, constants are
* initialized etc and the device is otherwise ready to start computing stuff.
*/
/*
void run_diagnostics(	float* lnrho, float* uu_x, float* uu_y, float* uu_z,
			float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z, 
                  	float* d_w_lnrho, float* d_w_uu_x, float* d_w_uu_y, float* d_w_uu_z,
			float* d_lnrho_dest, float* d_uu_x_dest, float* d_uu_y_dest, float* d_uu_z_dest,
			float* d_div_uu,
			float* d_umax, float* d_partial_result, int isubstep) 
{
	printf("Running diagnostics...\n\n");
	
	const int NUM_DEBUG_GRIDS = 9;

	/////////////////////////////////////////////////////////////////////////////////////
	// These arrays contain the pointers to the vector reduction functions             //
	////////////////Add model and GPU functions in these arrays//////////////////////////
	VecReductionFunctionHostPointer CPU_vec_reduction_functions[] = {&model_max_vec, &model_min_vec};
	VecReductionFunctionDevicePointer GPU_vec_reduction_functions[] = {&max_vec_cuda, &min_vec_cuda};
	/////////////////////////////////////////////////////////////////////////////////////

	const int NUM_VEC_REDUCTION_FUNCTIONS = sizeof(CPU_vec_reduction_functions) / sizeof(VecReductionFunctionHostPointer);

	int failures = 0;

	printf("Testing vec reduction functions...\n");
	for (int j=0; j < NUM_VEC_REDUCTION_FUNCTIONS; j++) {
		printf("Testing function #%d using debug grid...\n", j);
		for (int i=0; i < NUM_DEBUG_GRIDS; i++) {
			//Init debug grid
			printf("\t%d/%d... ", i, NUM_DEBUG_GRIDS-1);
			debug_grid_init(lnrho, uu_x, uu_y, uu_z, i);
	
			//Copy the grid to the device
			checkErr( cudaMemcpy(d_lnrho, lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
			checkErr( cudaMemcpy(d_uu_x, uu_x, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
			checkErr( cudaMemcpy(d_uu_y, uu_y, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
			checkErr( cudaMemcpy(d_uu_z, uu_z, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );	
			cudaDeviceSynchronize();

			//compute
			float CPU_result = CPU_vec_reduction_functions[j](uu_x, uu_y, uu_z);

			float GPU_result;
			GPU_vec_reduction_functions[j](d_umax, d_partial_result, d_uu_x, d_uu_y, d_uu_z);
			cudaDeviceSynchronize();
			cudaMemcpy(&GPU_result, (float*)d_umax, sizeof(float), cudaMemcpyDeviceToHost); 
			cudaDeviceSynchronize();
	

			float error = abs(CPU_result - GPU_result);
			const float epsilon = 0.00001f; //Because of the CPU-GPU floating-point differences
			if (error > epsilon) {
				printf("FAIL!\n");
				printf("\tCPU result: %f\n", CPU_result);
				printf("\tGPU result: %f\n", GPU_result);
				failures++;
			}
			else { printf("OK!\n"); }
		}
	}


	/////////////////////////////////////////////////////////////////////////////////////
	// These arrays contain the pointers to the scalar reduction functions             //
	////////////////Add model and GPU functions in these arrays//////////////////////////
	ScalReductionFunctionHostPointer CPU_scal_reduction_functions[] = {&model_max_scal};
	ScalReductionFunctionDevicePointer GPU_scal_reduction_functions[] = {&max_scal_cuda};
	/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////
	const int NUM_SCAL_REDUCTION_FUNCTIONS = sizeof(CPU_scal_reduction_functions) / sizeof(ScalReductionFunctionHostPointer);

	printf("Testing scal reduction functions...\n");
	for (int j=0; j < NUM_SCAL_REDUCTION_FUNCTIONS; j++) {
		printf("Testing function #%d using debug grid...\n", j);
		for (int i=0; i < NUM_DEBUG_GRIDS; i++) {
			//Init debug grid
			printf("\t%d/%d... ", i, NUM_DEBUG_GRIDS-1);
			debug_grid_init(lnrho, uu_x, uu_y, uu_z, i);
	
			//Copy the grid to the device
			checkErr( cudaMemcpy(d_lnrho, lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
			checkErr( cudaMemcpy(d_uu_x, uu_x, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
			checkErr( cudaMemcpy(d_uu_y, uu_y, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
			checkErr( cudaMemcpy(d_uu_z, uu_z, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );	
			cudaDeviceSynchronize();

			//compute
			float CPU_result = CPU_scal_reduction_functions[j](uu_x);

			float GPU_result;
			GPU_scal_reduction_functions[j](d_umax, d_partial_result, d_uu_x);
			cudaDeviceSynchronize();
			cudaMemcpy(&GPU_result, (float*)d_umax, sizeof(float), cudaMemcpyDeviceToHost); 
			cudaDeviceSynchronize();
	

			float error = abs(CPU_result - GPU_result);
			const float epsilon = 0.00001f; //Because of the CPU-GPU floating-point differences
			if (error > epsilon) {
				printf("FAIL!\n");
				printf("\tCPU result: %f\n", CPU_result);
				printf("\tGPU result: %f\n", GPU_result);
				failures++;
			}
			else { printf("OK!\n"); }
		}
	}


	float dt;
	cudaMemcpyFromSymbol(&dt, d_DT, sizeof(float));
	cudaDeviceSynchronize();	
	printf("Checking Rungekutta_2N_cuda with d_DT = %f...\n", dt);

	float *GPU_lnrho; //Log density
        float *GPU_uu_x, *GPU_uu_y, *GPU_uu_z; //velocities

        GPU_lnrho = (float*) malloc(sizeof(float)*GRID_SIZE);
        GPU_uu_x  = (float*) malloc(sizeof(float)*GRID_SIZE);
        GPU_uu_y  = (float*) malloc(sizeof(float)*GRID_SIZE);
        GPU_uu_z  = (float*) malloc(sizeof(float)*GRID_SIZE);

	for (int i=0; i < NUM_DEBUG_GRIDS; i++) {
		//Init debug grid
		printf("\t%d/%d... ", i, NUM_DEBUG_GRIDS-1);
		debug_grid_init(lnrho, uu_x, uu_y, uu_z, i);

		//Copy the grid to the device
		checkErr( cudaMemcpy(d_lnrho, lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
		checkErr( cudaMemcpy(d_uu_x, uu_x, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
		checkErr( cudaMemcpy(d_uu_y, uu_y, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
		checkErr( cudaMemcpy(d_uu_z, uu_z, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );	
		cudaDeviceSynchronize();
		
		//GPU	
		rungekutta2N_cuda(d_lnrho, d_uu_x, d_uu_y, d_uu_z, 
				  d_w_lnrho, d_w_uu_x, d_w_uu_y, d_w_uu_z,
				  d_lnrho_dest, d_uu_x_dest, d_uu_y_dest, d_uu_z_dest, isubstep);
		cudaDeviceSynchronize();

		checkErr( cudaMemcpy(GPU_lnrho, d_lnrho_dest, sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
		checkErr( cudaMemcpy(GPU_uu_x,  d_uu_x_dest,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
		checkErr( cudaMemcpy(GPU_uu_y,  d_uu_y_dest,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
		checkErr( cudaMemcpy(GPU_uu_z,  d_uu_z_dest,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
		cudaDeviceSynchronize();

		float error = check_grids(lnrho, uu_x, uu_y, uu_z,
					GPU_lnrho, GPU_uu_x, GPU_uu_y, GPU_uu_z);
		
		const float epsilon = 0.00001f; //Cutoff in CPU/GPU floating-point error
		printf("\n\t\tTotal error: %f\n", error);
		if (error > epsilon) {
			printf("\t\tFAIL!\n");
			failures++;
		}
		else { printf("\t\tOK!\n"); }
	}
	
	//
	//printf("Checking Rungekutta_2N_cuda with Courant timestep...\n");
	//for (int i=0; i < NUM_DEBUG_GRIDS; i++) {
	//	//Init debug grid
	//	printf("\t%d/%d with d_DT = ", i, NUM_DEBUG_GRIDS-1);
	//	debug_grid_init(lnrho, uu_x, uu_y, uu_z, i);

	//	//Copy the grid to the device
	//	checkErr( cudaMemcpy(d_lnrho, lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
	//	checkErr( cudaMemcpy(d_uu_x, uu_x, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
	//	checkErr( cudaMemcpy(d_uu_y, uu_y, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
	//	checkErr( cudaMemcpy(d_uu_z, uu_z, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );	
	//	cudaDeviceSynchronize();
	//	
	//	//Compute timestep (Tested max_vec already, so no need to validate this if the
	//	//host part of timestep_cuda is correct)
	//	dt = timestep_cuda(d_umax, d_partial_result, d_uu_x, d_uu_y, d_uu_z);
	//	checkErr( cudaMemcpyToSymbol(d_DT, &dt, sizeof(float)) );
	//	printf("%f... ", dt);

	//	//Latest GPU version	
	//	rungekutta2N_cuda(d_lnrho, d_uu_x, d_uu_y, d_uu_z, 
	//			  d_w_lnrho, d_w_uu_x, d_w_uu_y, d_w_uu_z,
	//			  d_lnrho_dest, d_uu_x_dest, d_uu_y_dest, d_uu_z_dest, 
	//			  d_w_lnrho_dest, d_w_uu_x_dest, d_w_uu_y_dest, d_w_uu_z_dest,
	//			  d_ddx_uu_x, d_ddy_uu_y, d_ddz_uu_z);
	//	cudaDeviceSynchronize();
	//	
	//	checkErr( cudaMemcpy(GPU_lnrho, d_lnrho_dest, sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//	checkErr( cudaMemcpy(GPU_uu_x,  d_uu_x_dest,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//	checkErr( cudaMemcpy(GPU_uu_y,  d_uu_y_dest,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//	checkErr( cudaMemcpy(GPU_uu_z,  d_uu_z_dest,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//	cudaDeviceSynchronize();

	//	checkErr( cudaMemcpy(d_lnrho, lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
	//	checkErr( cudaMemcpy(d_uu_x, uu_x, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
	//	checkErr( cudaMemcpy(d_uu_y, uu_y, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
	//	checkErr( cudaMemcpy(d_uu_z, uu_z, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );	
	//	cudaDeviceSynchronize();

	//	//Model GPU version	
	//	model_rungekutta2N_cuda(d_lnrho, d_uu_x, d_uu_y, d_uu_z, 
	//			  d_w_lnrho, d_w_uu_x, d_w_uu_y, d_w_uu_z,
	//			  d_lnrho_dest, d_uu_x_dest, d_uu_y_dest, d_uu_z_dest, 
	//			  d_w_lnrho_dest, d_w_uu_x_dest, d_w_uu_y_dest, d_w_uu_z_dest);
	//	cudaDeviceSynchronize();

	//	checkErr( cudaMemcpy(lnrho, d_lnrho_dest, sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//	checkErr( cudaMemcpy(uu_x,  d_uu_x_dest,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//	checkErr( cudaMemcpy(uu_y,  d_uu_y_dest,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//	checkErr( cudaMemcpy(uu_z,  d_uu_z_dest,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//	cudaDeviceSynchronize();

	//	float error = check_grids(lnrho, uu_x, uu_y, uu_z,
	//				GPU_lnrho, GPU_uu_x, GPU_uu_y, GPU_uu_z);
	//	
	//	//const float epsilon = 0.00001f; //Cutoff in CPU/GPU floating-point error
	//	const float epsilon = 0.01f; //Cutoff in CPU/GPU floating-point error
	//	printf("\n\t\tTotal error: %f\n", error);
	//	if (error > epsilon) {
	//		printf("\t\tFAIL!\n");
	//		failures++;
	//	}
	//	else { printf("\t\tOK!\n"); }
	//}
	//

	//TODO boundary conditions check

        free(GPU_lnrho);
        free(GPU_uu_x);
        free(GPU_uu_y);
        free(GPU_uu_z);



	printf("Diagnostics done. Failures found: %d.\n", failures);
}
*/
