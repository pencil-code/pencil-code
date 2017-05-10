extern "C" bool finalizeGpu(float *uu_x, float *uu_y, float *uu_z, float *lnrho);
extern "C" void RKintegration(float *uu_x, float *uu_y, float *uu_z, float *lnrho, int mx, int my, int mz, int nghost, int isubstep);
extern "C" void intitializeGPU(float *uu_x, float *uu_y, float *uu_z, float *lnrho, int nx, int ny, int nz, int nghost, float *x, float *y, float *z, float nu, float cs2);
void rungekutta2N_cuda(	float* d_lnrho, float* d_uu_x, float* d_uu_y, float* d_uu_z, 
                  	float* d_w_lnrho, float* d_w_uu_x, float* d_w_uu_y, float* d_w_uu_z,
			float* d_lnrho_dest, float* d_uu_x_dest, float* d_uu_y_dest, float* d_uu_z_dest, int isubstep);



