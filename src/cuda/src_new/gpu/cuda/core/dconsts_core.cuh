/*
*   GPU constants shared among all different implementations.
*   In order to avoid code duplication, EXTERN macro is used to define
*   the extern keyword for constants used in all other files but the one
*   which does the actual definition.
*/
#pragma once
#include "common/datatypes.h"

#ifdef INCLUDED_FROM_DCONST_DEFINER
    #define EXTERN
#else
    #define EXTERN extern 
#endif

EXTERN __constant__ int d_mx, d_my, d_mz; // Grid dims including the boundaries
EXTERN __constant__ int d_nx, d_ny, d_nz; // Dimensions of the computational domain
EXTERN __constant__ int d_nx_min, d_ny_min, d_nz_min; // Starting indices of the computational domain
EXTERN __constant__ int d_nx_max, d_ny_max, d_nz_max; // Ending indices (exclusive) of the computational domain
EXTERN __constant__ int d_bound_size; // Size of the boundary zone
EXTERN __constant__ real d_DIFF1_DX_DIV, d_DIFF1_DY_DIV, d_DIFF1_DZ_DIV; //Constants for diff operations
EXTERN __constant__ real d_DIFF2_DX_DIV, d_DIFF2_DY_DIV, d_DIFF2_DZ_DIV;
EXTERN __constant__ real d_DIFFMN_DXDY_DIV, d_DIFFMN_DYDZ_DIV, d_DIFFMN_DXDZ_DIV; 
EXTERN __constant__ real d_NU_VISC; //Viscosity
EXTERN __constant__ real d_CS2_SOUND; //Speed of sound
