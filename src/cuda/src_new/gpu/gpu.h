/*
*   Interface for accessing the GPU functions
*/
#pragma once
#include "common/config.h"
#include "common/grid.h"
#include "common/slice.h"
#include "common/forcing.h"
#include "diagnostics/diagnostics.h"
#ifdef GPU_ASTAROTH
  #include "common/PC_moduleflags.h"
#endif


typedef enum {
    CUDA_GENERIC = 0,
    CUDA_19P,
    CUDA_55P,
    CUDA_MAXWELL,
    CUDA_FOR_PENCIL,
    NUM_IMPLEMENTATIONS
} ImplType;

extern const char* impl_type_names[]; 

/////Description of the multi-GPU interface///////////////////////////////////////////
/*
NOTE: All of the following functions operate on all GPUs on the node unless otherwise stated.

GPUInitFunc: 
	Starting point of all GPU computation. Handles the allocation and
        initialization of *all memory needed on all GPUs in the node*. In other words,
        setups everything GPU-side so that calling any other GPU interface function
        afterwards does not result in illegal memory accesses. 

GPUDestroyFunc:
	Opposite of GPUInitFunc. Frees all GPU allocations and resets all devices in
        the node.

GPULoadFunc:
	Decomposes and loads the whole host grid g (as defined in "common/grid.h")
        into GPU memories.

GPUStoreFunc:
	Combines and stores the grids in GPU memories into the host grid g.

GPULoadForcingParamsFunc:
	Takes a ForcingParams struct and copies its contents into the constant
        memories.

GPULoadOuterHalosFunc:
	Similar to GPULoadFunc, but loads only the ghost zone of the host grid into
        appropriate locations in GPU memory. TODO review this description.

GPUStoreInternalHalosFunc:
	Similar to GPUStoreFunc, but combines the ghost zones of the GPU grid and
        stores them back to the ghost zone of the host grid. TODO review this description.

GPUIntegrateFunc:
	Does the full RK3 integration (including boundary conditions)

GPUIntegrateStepFunc:
	Does a single RK3 step (without computing boundary conditions)

GPUBoundcondStepFunc:
	Applies the boundary conditions on the grids in GPU memory

GPUReduceFunc:
	Performs a reduction on the GPU grids. The possible ReductTypes and GridTypes
        can be found in "common/datatypes.h" (subject to change) and "common/grid.h".
        Usage f.ex. GPUReduce(MAX_SCAL_UU, LNRHO); finds the maximum in the lnrho arrays
        on the GPUs.

GPUGetSliceFunc:
	Similar to GPUStoreFunc, but instead stores a 2D slices to the Slice struct
        given as parameter
*/
////////////////////////////////////////////////////////////////////////////////
//Memory management functions
//Single-node
typedef void (*GPUInitFunc)(CParamConfig* cparamconf, RunConfig* runconf);
typedef void (*GPUInitializeFunc)(CParamConfig* cparamconf, RunConfig* runconf, const Grid & h_grid);
typedef void (*GPUDestroyFunc)();
typedef void (*GPULoadFunc)(Grid* g);
typedef void (*GPUStoreFunc)(Grid* g);
typedef void (*GPULoadForcingParamsFunc)(ForcingParams* fp);
typedef void (*GPUUpdateForcingCoefsFunc)(ForcingParams* fp);

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
extern GPUInitializeFunc         GPUInitialize;
extern GPUDestroyFunc            GPUDestroy;
extern GPULoadFunc               GPULoad;             //Load from host to device
extern GPUStoreFunc              GPUStore;            //Store from device to host
#ifndef GPU_ASTAROTH
extern GPULoadForcingParamsFunc  GPULoadForcingParams;
#endif
extern GPUUpdateForcingCoefsFunc GPUUpdateForcingCoefs;
extern GPULoadOuterHalosFunc     GPULoadOuterHalos;
extern GPUStoreInternalHalosFunc GPUStoreInternalHalos;

//GPU solver interface
extern GPUIntegrateFunc     GPUIntegrate;
extern GPUIntegrateStepFunc GPUIntegrateStep;
extern GPUBoundcondStepFunc GPUBoundcondStep;
extern GPUReduceFunc        GPUReduce;

//Misc
extern GPUGetSliceFunc  GPUGetSlice;
void GPUSelectImplementation(ImplType type);
