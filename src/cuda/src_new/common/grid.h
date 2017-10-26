#pragma once
#include "datatypes.h"
#include "config.h"


#if LINDUCTION == 1 //TODO Note: remember to change SliceType in slice.h also when modifying this part
    typedef enum {LNRHO = 0, UUX, UUY, UUZ, AX, AY, AZ, NUM_ARRS, NOT_APPLICABLE} GridType;
#else
    typedef enum {LNRHO = 0, UUX, UUY, UUZ, NUM_ARRS, NOT_APPLICABLE} GridType;
#endif

#ifdef INCLUDED_FROM_GRID_NAME_DEFINER
    //See slice.h for justification why these are defined here
    const char* grid_names[] = {"lnrho", "uu_x", "uu_y", "uu_z", "a_x", "a_y", "a_z"};
#else
    extern const char* grid_names[];
#endif

typedef struct {
    real* arr[NUM_ARRS];
} Grid;

typedef enum { 
    GRID_INCREASING = 0, 
    GRID_DECREASING,
    GRID_RANDOM, 
    GRID_COMP_DOMAIN_ONLY,
    GRID_BOUND_ZONES_ONLY,
    GRID_NO_FLOAT_ARITHMETIC_ERROR,
    GRID_WAVE,
    GRID_RAND_WITH_NEG,
    GRID_ALL_UNIQUE,
    GRID_ALL_ZERO,
    GRID_GAUSSIAN_RADIAL_EXPL,
    GRID_XWAVE,
    NUM_GRID_TYPES
} InitType;

#ifdef INCLUDED_FROM_INIT_TYPE_NAME_DEFINER
    //Names for grid types. TODO note. Possibly error prone, but how else to define these at compile time
    //without additional init functions or malloc/frees?
    //See slice.h for justification why these are defined here
    const char* init_type_names[] = {"GRID_INCREASING", 
                                    "GRID_DECREASING",
                                    "GRID_RANDOM", 
                                    "GRID_COMP_DOMAIN_ONLY",
                                    "GRID_BOUND_ZONES_ONLY",
                                    "GRID_NO_FLOAT_ARITHMETIC_ERROR",
                                    "GRID_WAVE",
                                    "GRID_RAND_WITH_NEG",
                                    "GRID_ALL_UNIQUE",
                                    "GRID_ALL_ZERO",
                                    "GRID_GAUSSIAN_RADIAL_EXPL",
                                    "GRID_XWAVE"};
#else
    extern const char* init_type_names[];
#endif

void grid_malloc(Grid* g, CParamConfig *cparams);

void grid_free(Grid* g);

void grid_clear(Grid* g, CParamConfig* cparams);

void grid_init(Grid* g, InitType type, CParamConfig* cparams, StartConfig* start_params);




















