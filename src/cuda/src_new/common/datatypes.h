/*
* Definitions for the datatypes used throughout the program
* 
*/
#pragma once
#include <float.h> // For DBL_MAX and FLT_MAX
#include "../../headers_c.h"

// for standalone Astaroth, DOUBLE_PRECISION is defined in CMakeLists.txt
// (See "add_definitions(-DDOUBLE_PRECISION)")
// for PC-Astaroth it is set by the PC Makefile automatically in PC_modulesources

#ifdef DOUBLE_PRECISION
    //typedef double real;
    #ifndef  REAL
       #define REAL double
    #endif
    const double REAL_MAX = DBL_MAX;
    const double REAL_MIN = DBL_MIN;
    const double REAL_EPSILON = DBL_EPSILON;
    #define REAL_FORMAT_SPECIFIER "%lf"
#else
    //typedef float real;
    #ifndef  REAL
       #define REAL float
    #endif
    const float REAL_MAX = FLT_MAX;
    const double REAL_MIN = FLT_MIN;
    const double REAL_EPSILON = FLT_EPSILON;
    #define REAL_FORMAT_SPECIFIER "%f"
#endif
#define real REAL

#ifdef GPU_ASTAROTH
    typedef enum {MAX_VEC=0, MIN_VEC, RMS_VEC, MAX_SCAL, MIN_SCAL, RMS_SCAL, RMS_EXP, SUM_SCAL, SUM_EXP, MAX_VEC_UU, MIN_VEC_UU,  RMS_VEC_UU, NUM_REDUCT_TYPES} ReductType;
#else
    typedef enum {MAX_VEC_UU, MIN_VEC_UU, RMS_VEC_UU, MAX_SCAL, MIN_SCAL, RMS_SCAL, RMS_EXP, NUM_REDUCT_TYPES} ReductType;
#endif
//This should be defined somewhere else (like reduce.h, but would prefer not to make a new header file just for a single enumerator)

//Vector of 3 real type values (float or double depending on the definition of "real")
struct vec3r {
    real x, y, z;
};

//Vector of 3 integers
struct vec3i {
    int x, y, z;
};
