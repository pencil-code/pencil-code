

#pragma once

void fillhalosindevice(float* d_halo, float* d_lnrho, int nx, int ny, int nz, int halo_depth);
void fillhalosinhost(float* d_halo, float* d_lnrho, int nx, int ny, int nz, int halo_depth);
void copyouterhalostodevice(float *grid, float *d_grid, float *halo, float *d_halo, int Nx, int Ny, int Nz, int halo_depth);
void copyinternalhalostohost(float *grid, float *d_grid, float *halo, float *d_halo, int Nx, int Ny, int Nz, int halo_depth);
void checkKernelErr();
cudaError_t checkErr(cudaError_t result);
