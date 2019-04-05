/*
*   Sets the GPU interface up by pointing the interface functions to the functions
*	implementing the interface.
*
*	For example: only GPUInit() is visible to someone who includes "gpu/gpu.h". 
*	GPUInit() is set here to point to some actual implementation, for example
*	init_cuda_generic() defined in "gpu/cuda/cuda_generic.cu". 	
*
*	The primary benefit of this approach is that the host code calling these
*	interface functions do not have to care whether we use, say, cuda, opencl, 
*	multi-GPU or some alternative boundary condition function. Everything
*	happens under the hood and the implementation details are hidden to the caller.
*	Therefore even if we make drastic changes in the GPU code, we don't have to
*	rewrite the CPU code parts.
*	
*/
#include "gpu.h"
#include "common/errorhandler.h"
#include "cuda/cuda_generic.cuh"

//Memory management interface
GPUInitFunc               GPUInit               = &init_cuda_generic;
GPUInitializeFunc         GPUInitialize         = &initialize_cuda_generic;
GPUDestroyFunc            GPUDestroy            = &destroy_cuda_generic;
GPULoadFunc               GPULoad               = &load_grid_cuda_generic;  //Load from host to device
GPUStoreFunc              GPUStore              = &store_grid_cuda_generic; //Store from device to host
GPULoadOuterHalosFunc     GPULoadOuterHalos     = &load_outer_halos_cuda_generic;
GPUStoreInternalHalosFunc GPUStoreInternalHalos = &store_internal_halos_cuda_generic;
#ifdef FORCING
#ifdef GPU_ASTAROTH
GPUUpdateForcingCoefsFunc GPUUpdateForcingCoefs = &update_forcing_coefs_cuda_generic;
#else
GPULoadForcingParamsFunc  GPULoadForcingParams  = &load_forcing_params_cuda_generic;
#endif
#endif

//GPU solver interface
GPUIntegrateFunc     GPUIntegrate     = &integrate_cuda_generic;
GPUIntegrateStepFunc GPUIntegrateStep = &integrate_step_cuda_generic;
GPUBoundcondStepFunc GPUBoundcondStep = &boundcond_step_cuda_generic;

//Misc
GPUReduceFunc    GPUReduce   = NULL;
GPUGetSliceFunc  GPUGetSlice = &get_slice_cuda_generic;

//Names for grid types. TODO note. Possibly error prone, but how else to define these at compile time
//without additional init functions or malloc/frees?

const char* impl_type_names[] = {"CUDA_GENERIC",
                                 "CUDA_19P",
                                 "CUDA_55P",
                                 "CUDA_MAXWELL",
                                 "CUDA_FOR_PENCIL"};

//#define DO_MINIMAL_BUILD
#ifdef DO_MINIMAL_BUILD
    void GPUSelectImplementation(ImplType type) 
    { 
        if (type != CUDA_GENERIC) {
            printf("Warning, tried to select some other implementation than CUDA_GENERIC, this has no effect since DO_MINIMAL_BUILD is set.\n");
            CRASH("Crashing just to make a point"); 
        }
    } 
#else
    //#include "cuda/cuda_19p.cuh"
    //#include "cuda/cuda_55p.cuh"
    //#include "cuda/cuda_maxwell.cuh"

    //Select the GPU implementation (yes, this could be done much more nicely
    //with classes and inheritance, but doing this the "Pure C Way" makes it
    //probably easier to interface with PC)
    void GPUSelectImplementation(ImplType type) {
        if (type == CUDA_GENERIC) {
            //Memory management interface
            GPUInit                 = &init_cuda_generic;
            GPUDestroy              = &destroy_cuda_generic;
            GPULoad                 = &load_grid_cuda_generic;  //Load from host to device
            GPUStore                = &store_grid_cuda_generic; //Store from device to host
#ifndef GPU_ASTAROTH
            #ifdef FORCING
            GPULoadForcingParams    = &load_forcing_params_cuda_generic;
            #endif
#endif
            GPULoadOuterHalos       = &load_outer_halos_cuda_generic;
            GPUStoreInternalHalos   = &store_internal_halos_cuda_generic;

            //GPU solver interface
            GPUIntegrate     = &integrate_cuda_generic;
            GPUBoundcondStep = &boundcond_step_cuda_generic;
            GPUIntegrateStep = &integrate_step_cuda_generic;
#ifndef GPU_ASTAROTH
            GPUReduce        = &reduce_cuda_generic;
#endif
            //Misc
            GPUGetSlice = &get_slice_cuda_generic;
        }/*else if (type == CUDA_19P) {
            //Memory management interface
            GPUInit    = &init_cuda_19p;
            GPUDestroy = &destroy_cuda_19p;
            GPULoad    = &load_grid_cuda_19p;  //Load from host to device
            GPUStore   = &store_grid_cuda_19p; //Store from device to host
            GPULoadForcingParams = NULL;

            //GPU solver interface
            GPUIntegrate     = &integrate_cuda_19p;
            GPUIntegrateStep = NULL;
            GPUBoundcondStep = NULL;
            GPUReduce        = NULL;

            //Misc
            GPUGetSlice = &get_slice_cuda_19p;
        } else if (type == CUDA_55P) {
            //Memory management interface
            GPUInit    = &init_cuda_55p;
            GPUDestroy = &destroy_cuda_55p;
            GPULoad    = &load_grid_cuda_55p;  //Load from host to device
            GPUStore   = &store_grid_cuda_55p; //Store from device to host
            GPULoadForcingParams = NULL;

            //GPU solver interface
            GPUIntegrate     = &integrate_cuda_55p;
            GPUIntegrateStep = NULL;
            GPUBoundcondStep = NULL;
            GPUReduce        = NULL;

            //Misc
            GPUGetSlice = &get_slice_cuda_55p;
        } else if (type == CUDA_MAXWELL) {
            //Memory management interface
            GPUInit              = &init_cuda_maxwell;
            GPUDestroy           = &destroy_cuda_maxwell;
            GPULoad              = &load_grid_cuda_maxwell;  //Load from host to device
            GPUStore             = &store_grid_cuda_maxwell; //Store from device to host
            GPULoadForcingParams = &load_forcing_params_cuda_maxwell;

            //GPU solver interface
            GPUIntegrate     = &integrate_cuda_maxwell;
            GPUBoundcondStep = &boundcond_step_cuda_maxwell;
            GPUIntegrateStep = &integrate_step_cuda_maxwell;
            GPUReduce        = &reduce_cuda_maxwell;

            //Misc
            GPUGetSlice = &get_slice_cuda_maxwell;
        }*/ else if (type == CUDA_FOR_PENCIL) {
            //Memory management interface
            GPUInit              = &init_cuda_generic;
	    GPUInitialize        = &initialize_cuda_generic;
            GPUDestroy           = &destroy_cuda_generic;
            GPULoad              = &load_grid_cuda_generic;
            GPUStore             = &store_grid_cuda_generic;
#ifdef FORCING
#ifndef GPU_ASTAROTH
            GPULoadForcingParams = NULL;
#endif
            GPUUpdateForcingCoefs= &update_forcing_coefs_cuda_generic;
#endif
            GPULoadOuterHalos    = &load_outer_halos_cuda_generic;
            GPUStoreInternalHalos= &store_internal_halos_cuda_generic;

            //GPU solver interface
            GPUIntegrate     = NULL;                   // only substeps called from PC
            GPUBoundcondStep = NULL;     	       // boundary conditions set in PC
            GPUIntegrateStep = &integrate_step_cuda_generic;
#ifdef GPU_ASTAROTH
            GPUReduce        = &reduce_cuda_PC;
#endif
            //Misc
            GPUGetSlice = NULL;
        } else {
            CRASH("Invalid implementation type!");
        }
    }
#endif
