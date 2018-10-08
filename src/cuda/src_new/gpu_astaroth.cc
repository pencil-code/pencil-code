/*                             gpu_astaroth_ansi.cu
                               --------------------

   Date:   8-Feb-2017
   Author: M. Rheinhardt && J. Pekkilae
   Description:
 ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
Comments: 
DATE 17-Feb-2017: Omer Anjum: Added description of the functions
DATE 27-Jul-2017: Johannes Pekkilae: Removed all deprecated stuff and added
basic functions for accessing the GPU interface
DATE 18-Sep-2017: Johannes Pekkilae: Rewrote parts and updated with the new Grid, ForcingParams etc. structures
*/

//General
#include "common/errorhandler.h"
#include "common/config.h"
#include "common/forcing.h"
#include "utils/utils.h"

//GPU interface
#include "gpu/gpu.h"
#include "gpu/cuda/cuda_generic.cuh"

//PC interface
#include "common/PC_moduleflags.h"
#include "../../cparam_c.h"
#include "../../cdata_c.h"
#include "../../sub_c.h"
#include "diagnostics/diagnostics.h"
#define EXTERN extern
#include "common/PC_module_parfuncs.h"

static Grid h_grid;
static real* halo_buffer;
static CParamConfig cparams(mx, my, mz, mw, nx, ny, nz, nw, l1, l2, m1, m2, n1, n2, nghost);
static RunConfig run_params;
static ForcingParams forcing_params;

bool copyOmer=false;

/***********************************************************************************************/
real max_advec()
{
        real uxmax, uymax, uzmax, maxadvec_=0.;
        return maxadvec_;
}
/***********************************************************************************************/
real max_diffus()
{
        real maxdiffus_=0.;
#ifdef VISCOSITY
        maxdiffus_=nu*dxyz_2[nghost];
#endif
#ifdef MAGNETIC
        maxdiffus_=max(maxdiffus_,eta*dxyz_2[nghost]);}
#endif
        return maxdiffus_;
}
/***********************************************************************************************/
//Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries
//TODO: note, isubstep starts from 0 on the GPU side (i.e. rk3 does substeps 0, 1, 2)
//
extern "C" void substepGPU(int isubstep, bool full=false)
{
#ifdef FORCING
    //Update forcing params
    if (isubstep == itorder) { 
        forcing_params.Update();    	                // calculate on CPU
        GPUUpdateForcingCoefs(&forcing_params);		// load into GPU
    }
#endif  
    if (lfirst && ldt) {
         real dt1_advec  = max_advec()/cdt;
         real dt1_diffus = max_diffus()/cdtv;
         real dt1_=sqrt(pow(dt1_advec,2) + pow(dt1_diffus,2));
         set_dt(dt1_);
    }
    //Transfer the updated ghost zone to the device(s) in the node 

    if (full)
    {
        GPULoad(&h_grid); 
    }
    else
    {
        //!!!GPULoadOuterHalos(&h_grid,halo_buffer);
        GPULoad(&h_grid); 
        //exchange_halos_cuda_generic(false);      // no circular halo exchange
    }

    //Integrate on the GPUs in this node
    //NOTE: In astaroth code, isubstep may be {0,1,2}, in PC it is {1,2,3}

/*if (isubstep==1){
printf("vor substep 1 \n");
real maxux=h_grid.arr[h_grid.UUX][0];
for (int i=1;i<134*134*134;i++){if (h_grid.arr[h_grid.UUX][i]>maxux){ maxux=(h_grid.arr[h_grid.UUX])[i];}}
printf("maxux= %.14f\n", maxux);
}*/
    if (ldiagnos) timeseries_diagnostics(h_grid);
    GPUIntegrateStep(isubstep-1, dt);

    //!!!GPUStoreInternalHalos(&h_grid,halo_buffer);
    GPUStore(&h_grid);

//CRASH("GPUDiagnostics");
}
/* ---------------------------------------------------------------------- */
extern "C" void registerGPU(real* farray)
{
        h_grid.Setup(farray);
}
/* ---------------------------------------------------------------------- */
//Setup the GPUs in the node to be ready for computation
extern "C" void initGPU()
{
    //Initialize GPUs in the node

    GPUSelectImplementation(CUDA_FOR_PENCIL);
    GPUInit(&cparams, &run_params);         //Allocs memory on the GPU and loads device constants
}
extern "C" void initializeGPU()
{
    //Setup configurations used for initializing and running the GPU code
#ifdef DOUBLE_PRECISION
printf("DOUBLE PRECISION!\n");
#endif
    	cparams.Setup(dx, dy, dz);

        GPUInitialize(&cparams, &run_params, h_grid);         // loads device constants

    // initialize diagnostics
       init_diagnostics();

       //cparams.print();
}
/***********************************************************************************************/
extern "C" void copyFarray(real *uu_x, real *uu_y, real *uu_z, real *lnrho) // pars only for backwards compatibility
{
        if (copyOmer)
        {
                printf("Stop: Inside full_inner to copy grid to host or inner halos to host\n");
        }
        else
        {
                GPUStore(&h_grid);
        }
}
/***********************************************************************************************/
//Destroy the GPUs in the node (not literally hehe)
extern "C" void finalizeGPU()
{
    //Deallocate everything on the GPUs and reset
    GPUDestroy();
}
/***********************************************************************************************/
