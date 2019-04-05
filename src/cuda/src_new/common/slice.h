#pragma once
#include "datatypes.h"
#include "config.h"
#ifdef GPU_ASTAROTH
  #include "common/PC_moduleflags.h"
#endif

//Here we use preprocessor magic to define arrays held by the Slice struct in
//a modular fashion. IMPORTANT NOTE: NUM_SLICES must be defined after the 
//computational arrays such that we can use it elsewhere in the code to
//loop over all these grid with generic loops.
//
// SLICE_TYPES holds all the enumerator-string tuples we need. At compile-time,
// macro #define SLICE(type, name) type resolves SLICE_TYPES to a list of types,
// which is then used when creating the enumerations. 
//  F.ex.   #define SLICE_TYPE(type, name) type 
//          { SLICE_TYPES } => 
//          { SLICE_TYPES_HYDRO SLICE_TYPES_INDUCTION ... }
//          { SLICE_TYPE(SLICE_LNRHO, "lnrho"), SLICE_TYPE(SLICE_UX, "uu_x"), ... }
//          { SLICE_LNRHO, SLICE_UUX, ...}
//
//  A macro is used similarly to extract the strings which correspond to each SliceType
//  s.t. the definition in slice.cc looks like
//          #define SLICE_TYPE(type, name) name
//          { SLICE_TYPES } =>
//          { "lnrho", "uu_x", ... }
//  Sidenote: we could also do just "#define SLICE_TYPE(type) #type" to extract
//  a string from the type names, but this would break the existing python scripts
//  since they accept only specific file names.
//
//Note: in pure C99 we could do this cleanly with designated initializers 
//(slice_names[] = {[LNRHO]="lnrho"}) but sadly there's no similar construct
//in C++ and we have to resort to this somewhat dirty hack.
//
//Note 2: when adding new arrays, modify the enumerators SliceType (slice.h),
//GridType (grid.h) and also possibly InitType()
extern const char* slice_names[];

#ifdef HYDRO
    #define SLICE_TYPES_HYDRO \
            SLICE_TYPE(SLICE_LNRHO, "lnrho"), \
            SLICE_TYPE(SLICE_UUX,   "uu_x"), \
            SLICE_TYPE(SLICE_UUY,   "uu_y"), \
            SLICE_TYPE(SLICE_UUZ,   "uu_z"), \
            SLICE_TYPE(SLICE_UU,    "uu"),
#else
    #define SLICE_TYPES_HYDRO
#endif

#ifdef MAGNETIC
    #define SLICE_TYPES_INDUCTION \
            SLICE_TYPE(SLICE_AAX,    "aa_x"), \
            SLICE_TYPE(SLICE_AAY,    "aa_y"), \
            SLICE_TYPE(SLICE_AAZ,    "aa_z"),
#else
    #define SLICE_TYPES_INDUCTION
#endif

#define SLICE_TYPES \
        SLICE_TYPES_HYDRO \
        SLICE_TYPES_INDUCTION \
        SLICE_TYPE(NUM_SLICES,  "num_slices")


#define SLICE_TYPE(type, name) type
typedef enum { SLICE_TYPES } SliceType;
#undef SLICE_TYPE

typedef struct {
    real* arr[NUM_SLICES];
} Slice;

void slice_malloc(Slice* s, CParamConfig *cparams, RunConfig *run_params);

void slice_free(Slice* s);
