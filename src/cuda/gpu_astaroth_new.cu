/*                             gpu_astaroth_ansi.cu
                              --------------------
*/

/* Date:   8-Feb-2017
   Author: M. Rheinhardt
   Description:
 ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
Comments: 
DATE 17-Feb-2017: Omer Anjum: Added description of the functions
DATE 27-Jul-2017: Johannes Pekkilae: Removed all deprecated stuff
*/

extern "C"
{
/* ---------------------------------------------------------------------- */
void RKintegration( float *uu_x, float *uu_y, float *uu_z, float *lnrho, 
                    int mx, int my, int mz, 
                    int nghost, int isubstep)
{
	return;
}

//void intitializeGPU(float *uu_x, float *uu_y, float *uu_z, float *lnrho, int nx, int ny, int nz, int nghost, float *x, float *y, float *z, float NU_VISC, float cs2_sound){ 
void intitializeGPU(float *uu_x, float *uu_y, float *uu_z, float *lnrho, 
                    int nx, int ny, int nz, 
                    int nghost, 
                    float *x, float *y, float *z, 
                    float nu, float cs2)
{ 
    
}


bool finalizeGpu(float *uu_x, float *uu_y, float *uu_z, float *lnrho)
{
    //Destroy  

	//Reset device
	cudaDeviceReset();

	return EXIT_SUCCESS;
}

/* ---------------------------------------------------------------------- */
}
