/*
*   Interface for accessing the GPU functions
*/
#pragma once
#include "cuda/cuda_generic.cuh"


//Memory management functions
typedef void (*GPUInitFunc)(CParamConfig* cparamconf, RunConfig* runconf);
typedef void (*GPUDestroyFunc)();
typedef void (*GPULoadFunc)(real* lnrho, real* uu_x, real* uu_y, real* uu_z);
typedef void (*GPUStoreFunc)(real* lnrho, real* uu_x, real* uu_y, real* uu_z);

//GPU solver functions
typedef void (*GPUIntegrateFunc)(real dt);
typedef void (*GPUIntegrateStepFunc)(int isubstep, real dt);
typedef real (*GPUReduceFunc)(ReductType t, ArrType a);

//Misc GPU functions
typedef void (*GPUGetSliceFunc)(real* slice_lnrho, real* slice_uu, real* slice_uu_x, real* slice_uu_y, real* slice_uu_z);


//Memory management interface
static GPUInitFunc    GPUInit    = &init_cuda_generic;
static GPUDestroyFunc GPUDestroy = &destroy_cuda_generic;
static GPULoadFunc    GPULoad    = &load_grid_cuda_generic;  //Load from host to device
static GPUStoreFunc   GPUStore   = &store_grid_cuda_generic; //Store from device to host

//GPU solver interface
static GPUIntegrateFunc     GPUIntegrate     = &integrate_cuda_generic;
static GPUIntegrateStepFunc GPUIntegrateStep = &integrate_step_cuda_generic;
static GPUReduceFunc        GPUReduce        = &reduce_cuda_generic;

//Misc
static GPUGetSliceFunc  GPUGetSlice = &get_slice_cuda_generic;


#define DO_MINIMAL_BUILD
#ifdef DO_MINIMAL_BUILD
typedef enum {CUDA_GENERIC, NUM_IMPLEMENTATIONS} ImplType;
static const char* impl_type_names[] = {"CUDA_GENERIC"};
static void GPUSelectImplementation(ImplType type) { } 
#else
///////////////////////////////////////////////////////////////////////////////
#include "cuda/cuda_19p.cuh"
#include "cuda/cuda_55p.cuh"
#include "cuda/cuda_maxwell.cuh"

typedef enum {
    CUDA_GENERIC,
    CUDA_19P,
    CUDA_55P,
    CUDA_MAXWELL,
    NUM_IMPLEMENTATIONS
} ImplType;

//Names for grid types. TODO note. Possibly error prone, but how else to define these at compile time
//without additional init functions or malloc/frees?
static const char* impl_type_names[] = {"CUDA_GENERIC",
                                        "CUDA_19P",
                                        "CUDA_55P",
                                        "CUDA_MAXWELL"};


//Select the GPU implementation (yes, this could be done much more nicely
//with classes and inheritance, but doing this the "Pure C Way" makes it
//probably easier to interface with PC)
static void GPUSelectImplementation(ImplType type) {
    if (type == CUDA_GENERIC) {
        //Memory management interface
        GPUInit    = &init_cuda_generic;
        GPUDestroy = &destroy_cuda_generic;
        GPULoad    = &load_grid_cuda_generic;  //Load from host to device
        GPUStore   = &store_grid_cuda_generic; //Store from device to host

        //GPU solver interface
        GPUIntegrate     = &integrate_cuda_generic;
        GPUIntegrateStep = &integrate_step_cuda_generic;
        GPUReduce        = &reduce_cuda_generic;

        //Misc
        GPUGetSlice = &get_slice_cuda_generic;
    } else if (type == CUDA_19P) {
        //Memory management interface
        GPUInit    = &init_cuda_19p;
        GPUDestroy = &destroy_cuda_19p;
        GPULoad    = &load_grid_cuda_19p;  //Load from host to device
        GPUStore   = &store_grid_cuda_19p; //Store from device to host

        //GPU solver interface
        GPUIntegrate     = &integrate_cuda_19p;
        GPUIntegrateStep = NULL;
        GPUReduce        = NULL;

        //Misc
        GPUGetSlice = &get_slice_cuda_19p;
    } else if (type == CUDA_55P) {
        //Memory management interface
        GPUInit    = &init_cuda_55p;
        GPUDestroy = &destroy_cuda_55p;
        GPULoad    = &load_grid_cuda_55p;  //Load from host to device
        GPUStore   = &store_grid_cuda_55p; //Store from device to host

        //GPU solver interface
        GPUIntegrate     = &integrate_cuda_55p;
        GPUIntegrateStep = NULL;
        GPUReduce        = NULL;

        //Misc
        GPUGetSlice = &get_slice_cuda_55p;
    } else if (type == CUDA_MAXWELL) {
        //Memory management interface
        GPUInit    = &init_cuda_maxwell;
        GPUDestroy = &destroy_cuda_maxwell;
        GPULoad    = &load_grid_cuda_maxwell;  //Load from host to device
        GPUStore   = &store_grid_cuda_maxwell; //Store from device to host

        //GPU solver interface
        GPUIntegrate     = &integrate_cuda_maxwell;
        GPUIntegrateStep = NULL;
        GPUReduce        = &reduce_cuda_maxwell;

        //Misc
        GPUGetSlice = &get_slice_cuda_maxwell;
    } else {
        CRASH("Invalid implementation type!");
    }
}
///////////////////////////////////////////////////////////////////////////////
#endif








































