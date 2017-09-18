/*                             gpu_astaroth_ansi.cu
                              --------------------
*/

/* Date:   8-Feb-2017
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
#include "common/datatypes.h"
#include "common/errorhandler.h"
#include "common/config.h"

//#include "common/defines_dims_PC.h"//Could be used to load some global vars from PC to Astaroth
//#include "common/defines_PC.h"     //during runtime (with func setup_configs())

//GPU interface
#include "gpu/gpu.h"


extern "C"
{

static ForcingParams forcing_params;

/* ---------------------------------------------------------------------- */

static ForcingParams forcing_params_init()
{
    ForcingParams fp;
    fp.forcing_enabled = true;

	//Determine the k-vector arrays size 
	fp.nk = 1; //This can be one since PC handles the actual random k vector array
	//Allocate vector arrays
	fp.kk_x = (real*) malloc(sizeof(real)*fp.nk);
	fp.kk_y = (real*) malloc(sizeof(real)*fp.nk);
	fp.kk_z = (real*) malloc(sizeof(real)*fp.nk);

	//Get the forcing vector from PC.
    fp.kk_x[0] = -1; //TODO
    fp.kk_y[0] = -1; //TODO
    fp.kk_z[0] = -1; //TODO
    
    return fp;
}


static void forcing_params_destroy(ForcingParams* fp)
{
    free(fp->kk_x); fp->kk_x = NULL;
    free(fp->kk_y); fp->kk_y = NULL;
    free(fp->kk_z); fp->kk_z = NULL;
}


void forcing_params_update(ForcingParams* fp)
{
    //TODO fetch these from PC
    /*
    fp->k_idx = choose_random_k(fp->nk); 
    fp->phi = choose_random_phi();
    fkt_forcing_coefficient(&fp->kk_part_x, &fp->kk_part_y, &fp->kk_part_z, fp->kk_x[fp->k_idx], fp->kk_y[fp->k_idx], fp->kk_z[fp->k_idx], dt, run_config);
    */
}

//Do the 'isubstep'th integration step on all GPUs on the node and handle boundaries
//TODO: note, isubstep starts from 0 on the GPU side (i.e. rk3 does substeps 0, 1, 2)
void RKintegration( real *uu_x, real *uu_y, real *uu_z, real *lnrho, 
                    int mx, int my, int mz, 
                    int nghost, int isubstep)
{
    const real dt = -1.0;//TODO: read from PC

    //Transfer inner halos back to host (TODO: not yet implemented)
    //GPUStoreInnerHalos(host_grid) //copyinnerhalostohost

    //Compute boundary conditions with Pencil Code (TODO: not yet implemented)
    //pencil_code_solve_boundary_conditions_on_cpu()

    //Transfer the updated ghost zone to the devices in the node (TODO: not yet implemented)
    //GPULoadOuterHalos(host_grid) //copyouterhalostodevice

    //Update forcing params
    //NOTE: In astaroth code, the isubstep may be {0, 1, 2}. This may not be the case in PC
    if (LFORCING && isubstep == 2) { 
        forcing_params_update(&forcing_params);
        GPULoadForcingParams(&forcing_params);
    }
    
    //Integrate on the GPUs in this node
    GPUIntegrateStep(isubstep, dt);
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


static Grid setup_grid()
{
    Grid grid;
    grid.arr[LNRHO] = NULL;// should be something like grid.arr[LNRHO] = some_farr_in_pc; 
    grid.arr[UUX] = NULL;
    //... etc. See the definition of Grid and GridType in common/grid.h 
    return grid;
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

    //Setup various parameters
    forcing_params = forcing_params_init();
    GPULoadForcingParams(&forcing_params);

    //Initialize GPUs in the node
    GPUSelectImplementation(CUDA_GENERIC);
    GPUInit(&cparams, &run_params); //Allocs memory on the GPU and loads device constants
    
    //Load the whole grid from host to device
    Grid grid = setup_grid();
    GPULoad(&grid); //Loads the whole grid from host to device

    //TODO: Any optional steps, for example store the first GPU slice to slice arrays on host:
    //GPUGetSlice(slice_lnrho, slice_uu, slice_uu_x, slice_uu_y, slice_uu_z);
}


//Destroy the GPUs in the node (not literally hehe)
bool finalizeGPU(real *uu_x, real *uu_y, real *uu_z, real *lnrho)
{
    //Deallocate everything on the GPUs and reset
    GPUDestroy();

    //TODO: Deallocate all host memory possibly allocated in initializeGPU() etc
    forcing_params_destroy(&forcing_params);

	return EXIT_SUCCESS;
}



int main()
{
    printf("Compiled and built successfully! Exiting\n");

    //initializeGPU(TODO the appropriate params here);
    //RKIntegration(TODO the appropriate params here);
    //finalizeGPU(TODO the appropriate params here);

    return EXIT_SUCCESS;
}

/* ---------------------------------------------------------------------- */
}





























