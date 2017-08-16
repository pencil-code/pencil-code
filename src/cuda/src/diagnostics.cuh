#pragma once
#include "smem.cuh"

#define DEBUG
#include <stdio.h>
#include <assert.h>

//Checks CUDA errors; replaced with nop if DEBUG not defined
inline
cudaError_t checkErr(cudaError_t result) {
#ifdef DEBUG
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s \n", 
            cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
#endif
  return result;
}

inline
void checkKernelErr() {
#ifdef DEBUG
	checkErr( cudaPeekAtLastError() );
	checkErr( cudaDeviceSynchronize() );	
#endif
	return;
}

void print_run_config();

void print_init_config();

void print_additional_defines();

void print_grid_data(float* grid);

__device__ void check_for_nan_inf_variable(int code, float var);

__device__ void check_for_nan_inf(      int step_number, int grid_idx_x, int grid_idx_y, int grid_idx_z,
                                        int sid_row, int sid_col, int sid_depth,
                                        float s_lnrho[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH],
                                        float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH],
                                        float s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH],
                                        float s_uu_z[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH]       );

__global__ void check_grid_for_nan_cuda(float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z, int* d_nan_count);

float check_grid_for_nan(float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z);

void run_diagnostics(	float* lnrho, float* uu_x, float* uu_y, float* uu_z,
			float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z, 
                  	float* d_w_lnrho, float* d_w_uu_x, float* d_w_uu_y, float* d_w_uu_z,
			float* d_lnrho_dest, float* d_uu_x_dest, float* d_uu_y_dest, float* d_uu_z_dest,
			float* d_div_uu,
			float* d_umax, float* d_partial_result, int isubstep);


