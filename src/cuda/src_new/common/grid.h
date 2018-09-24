/*
* Declarations for a generic Grid structure. 
*   See slice.h for more information about the macros used here
*/
#pragma once

#include "datatypes.h"
#include "config.h"

/*
//GridType macros///////////////////////////////////////////////////////////////
extern const char* grid_names[];

#define GRID_TYPES_HYDRO \
        GRID_TYPE(LNRHO, "lnrho"), \
        GRID_TYPE(UUX,   "uu_x"), \
        GRID_TYPE(UUY,   "uu_y"), \
        GRID_TYPE(UUZ,   "uu_z"),

#if LINDUCTION == 1
    #define GRID_TYPES_INDUCTION \
            GRID_TYPE(AX,    "ax"), \
            GRID_TYPE(AY,    "ay"), \
            GRID_TYPE(AZ,    "az"),
#else
    #define GRID_TYPES_INDUCTION
#endif

#define GRID_TYPES_OTHER \
        GRID_TYPE(NUM_ARRS,  "num_arrs"), \
        GRID_TYPE(NOT_APPLICABLE,  "N/A")

#define GRID_TYPES \
        GRID_TYPES_HYDRO \
        GRID_TYPES_INDUCTION \
        GRID_TYPES_OTHER

#define GRID_TYPE(type, name) type
typedef enum { GRID_TYPES } GridType;
#undef GRID_TYPE
*/
typedef int GridType;

////////////////////////////////////////////////////////////////////////////////

struct Grid {
public:
    //real* arr[NUM_ARRS];
    real** arr=NULL;
    int UUX=-1, UUY=-1, UUZ=-1, RHO=-1, LNRHO=-1, AAX=-1, AAY=-1, AAZ=-1, SS=-1, NUM_ARRS=0;
    const char* names[8] = {"uu_x", "uu_y", "uu_z", "lnrho", "ss", "aa_x", "aa_y", "aa_z"};

    Grid();
    Grid & operator=(const Grid &);
    void Setup(real*);
};

/*
//InitType macros///////////////////////////////////////////////////////////////
extern const char* init_type_names[];

#define INIT_TYPES \
        INIT_TYPE(GRID_INCREASING, "GRID_INCREASING"), \
        INIT_TYPE(GRID_DECREASING, "GRID_DECREASING"), \
        INIT_TYPE(GRID_RANDOM, "GRID_RANDOM"), \
        INIT_TYPE(GRID_COMP_DOMAIN_ONLY, "GRID_COMP_DOMAIN_ONLY"), \
        INIT_TYPE(GRID_BOUND_ZONES_ONLY, "GRID_BOUND_ZONES_ONLY"), \
        INIT_TYPE(GRID_NO_FLOAT_ARITHMETIC_ERROR, "GRID_NO_FLOAT_ARITHMETIC_ERROR"), \
        INIT_TYPE(GRID_WAVE, "GRID_WAVE"), \
        INIT_TYPE(GRID_RAND_WITH_NEG, "GRID_RAND_WITH_NEG"), \
        INIT_TYPE(GRID_ALL_UNIQUE, "GRID_ALL_UNIQUE"), \
        INIT_TYPE(GRID_ALL_ZERO, "GRID_ALL_ZERO"), \
        INIT_TYPE(GRID_GAUSSIAN_RADIAL_EXPL, "GRID_GAUSSIAN_RADIAL_EXPL"), \
        INIT_TYPE(GRID_XWAVE, "GRID_XWAVE"), \
        INIT_TYPE(NUM_INIT_TYPES,  "num_init_types")

#define INIT_TYPE(type, name) type
typedef enum { INIT_TYPES } InitType;
#undef INIT_TYPE
////////////////////////////////////////////////////////////////////////////////

void grid_init(Grid* g, InitType type, CParamConfig* cparams, StartConfig* start_params);
*/

//Grid allocator functions//////////////////////////////////////////////////////
void grid_malloc(Grid* g, CParamConfig *cparams);

void grid_free(Grid* g);

void grid_clear(Grid* g, CParamConfig* cparams);
////////////////////////////////////////////////////////////////////////////////
