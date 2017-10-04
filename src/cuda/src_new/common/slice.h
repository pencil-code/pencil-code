#pragma once
#include "datatypes.h"
#include "config.h"

#if LINDUCTION == 1 //TODO Note: remember to change GridType in grid.h also when modifying this part
    typedef enum {SLICE_LNRHO = 0, SLICE_UUX, SLICE_UUY, SLICE_UUZ, SLICE_UU, SLICE_AX, SLICE_AY, SLICE_AZ, NUM_SLICES} SliceType;
#else
    typedef enum {SLICE_LNRHO = 0, SLICE_UUX, SLICE_UUY, SLICE_UUZ, SLICE_UU, NUM_SLICES} SliceType;    
#endif

#ifdef INCLUDED_FROM_SLICE_NAME_DEFINER
    //Usually these would be defined in slice.cc, however for readability
    //and making sure that these have the same ordering as SliceTypes it's
    //better that we use this ifdef trick and define slice_names them here.
    const char* slice_names[] = {"lnrho", "uu_x", "uu_y", "uu_z", "uu", "ax", "ay", "az"};
#else
    extern const char* slice_names[];
#endif

typedef struct {
    real* arr[NUM_SLICES];
} Slice;

void slice_malloc(Slice* s, CParamConfig *cparams, RunConfig *run_params);

void slice_free(Slice* s);
