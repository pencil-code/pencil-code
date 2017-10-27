/*
*   Interface for accessing the GPU functions
*/
#pragma once
#include "common/config.h"
#include "common/grid.h"
#include "common/slice.h"
#include "common/forcing.h"


typedef enum {
    CUDA_GENERIC = 0,
    CUDA_19P,
    CUDA_55P,
    CUDA_MAXWELL,
    NUM_IMPLEMENTATIONS
} ImplType;

extern const char* impl_type_names[]; 


//Memory management functions
//Single-node
typedef void (*GPUInitFunc)(CParamConfig* cparamconf, RunConfig* runconf);
typedef void (*GPUDestroyFunc)();
typedef void (*GPULoadFunc)(Grid* g);
typedef void (*GPUStoreFunc)(Grid* g);
typedef void (*GPULoadForcingParamsFunc)(ForcingParams* fp);
//Multi-node
typedef void (*GPULoadOuterHalosFunc)(Grid* g, real* halo);
typedef void (*GPUStoreInternalHalosFunc)(Grid* g, real* halo);


//GPU solver functions
typedef void (*GPUIntegrateFunc)(real dt);
typedef void (*GPUIntegrateStepFunc)(int isubstep, real dt);
typedef void (*GPUBoundcondStepFunc)();
typedef real (*GPUReduceFunc)(ReductType t, GridType a);

//Misc GPU functions
typedef void (*GPUGetSliceFunc)(Slice* s);


//Memory management interface
extern GPUInitFunc               GPUInit;
extern GPUDestroyFunc            GPUDestroy;
extern GPULoadFunc               GPULoad;  //Load from host to device
extern GPUStoreFunc              GPUStore; //Store from device to host
extern GPULoadForcingParamsFunc  GPULoadForcingParams;
extern GPULoadOuterHalosFunc     GPULoadOuterHalos;
extern GPUStoreInternalHalosFunc GPULoadInternalHalos;

//GPU solver interface
extern GPUIntegrateFunc     GPUIntegrate;
extern GPUIntegrateStepFunc GPUIntegrateStep;
extern GPUBoundcondStepFunc GPUBoundcondStep;
extern GPUReduceFunc        GPUReduce;

//Misc
extern GPUGetSliceFunc  GPUGetSlice;
void GPUSelectImplementation(ImplType type);
















































































