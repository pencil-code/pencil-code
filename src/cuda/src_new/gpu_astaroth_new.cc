/*                             gpu_astaroth_ansi.cu
                              --------------------
*/

/* Date:   8-Feb-2017
   Author: M. Rheinhardt
   Description:
 ANSI C and standard library callable function wrappers for ASTAROTH-nucleus functions to be called from Fortran.
Comments: 
DATE 17-Feb-2017: Omer Anjum: Added description of the functions
DATE 27-Jul-2017: Johannes Pekkilae: Removed all deprecated stuff and added
basic functions for accessing the GPU interface
*/


//General
#include "common/datatypes.h"
#include "common/errorhandler.h"
#include "common/config.h"

//#include "common/defines_dims_PC.h"//Could be used to load some global vars from PC to Astaroth
//#include "common/defines_PC.h"     //during runtime (in setup_condif())

//GPU interface
#include "gpu/gpu.h"


extern "C"
{
/* ---------------------------------------------------------------------- */

//Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries
//TODO: note, isubstep starts from 0 on the GPU side (i.e. rk3 does substeps 0, 1, 2)
void RKintegration( real *uu_x, real *uu_y, real *uu_z, real *lnrho, 
                    int mx, int my, int mz, 
                    int nghost, int isubstep)
{
    const real dt = -1.0;//TODO: read from PC

	integrate_step_cuda_generic(isubstep, dt);

    //Transfer inner halos back to host (TODO: not yet implemented)
    //GPULoadOuterHalos(lnrho, uu_x, uu_y, uu_z) //copyouterhalostodevice
    
    //Share the local boundaries (within the node) and integrate
    GPUIntegrateStep(isubstep, dt);

    //Transfer inner halos back to host (TODO: not yet implemented)
    //GPUStoreInnerHalos(lnrho, uu_x, uu_y, uu_z); //copyinternalhalostohost 
}


//Setup configuration structs used for initializing and running the GPU code
static void setup_configs(CParamConfig* cparams, 
                          StartConfig*  start_params,
                          RunConfig*    run_params,
                          int nx, int ny, int nz,
                          int nghost,
                          real nu, real cs2)
{
    //We can either read the configs from the files in src/configs/ with
    //read_configs(cparams, start_params, run_params); 


    //Or we can fill the config structs manually, for example:
    cparams->nx = nx; // nx is the size of the computational domain
    cparams->ny = ny;
    cparams->nz = nz;

    //Fill rest of the cparamconfig (not needed once we get all of these values from PC)
    cparams->compute_missing_values();//Deduces f.ex. mx, l1 etc from nx, ny and nz.

    //cparams.nghost = nghost //Astaroth currently requires that nghost is known
    //at compile time (optimization reasons). See BOUND_SIZE in common/config.h.
    //TODO see also common/defines_dims_PC.h and common/defines_PC.h, these headers 
    //could potentially be included in this file and loaded into the configs
    //in this setup function, f.ex.
    //cparams->nx = PC_NX; and so on


    run_params->max_steps  = 101; //TODO read these from PC
    run_params->save_steps = 20;  //TODO read these from PC
    run_params->cdt        = 0.4; //TODO read these from PC
    run_params->cdtv       = 0.08;//TODO read these from PC
    run_params->nu_visc    = nu;
    run_params->cs_sound   = cs2;    

    //Print the configs(NOTE: not yet implemented)
    //cparams.print();
    //start_params.print();
    //run_params.print();
}


//Setup the GPUs in the node to be ready for computation
void intitializeGPU(real *uu_x, real *uu_y, real *uu_z, real *lnrho, 
                    int nx, int ny, int nz, 
                    int nghost, 
                    real *x, real *y, real *z, 
                    real nu, real cs2)
{
    //Setup configs
    CParamConfig cparams;
    StartConfig  start_params;
    RunConfig    run_params;
    setup_configs(&cparams, &start_params, &run_params, nx, ny, nz, nghost, nu, cs2);

    //Initialize GPUs in the node
    GPUSelectImplementation(CUDA_GENERIC);
    GPUInit(&cparams, &run_params); //Allocs memory on the GPU and loads device constants
    GPULoad(lnrho, uu_x, uu_y, uu_z); //Loads the whole grid from host to device

    //TODO: Any optional steps, for example store the first GPU slice to slice arrays on host:
    //GPUGetSlice(slice_lnrho, slice_uu, slice_uu_x, slice_uu_y, slice_uu_z);
}


//Destroy the GPUs in the node (not literally hehe)
bool finalizeGpu(real *uu_x, real *uu_y, real *uu_z, real *lnrho)
{
    //Deallocate everything on the GPUs and reset
    GPUDestroy();

    //TODO: Deallocate all host memory possibly allocated in initializeGPU() etc

	return EXIT_SUCCESS;
}

/* ---------------------------------------------------------------------- */
}





























