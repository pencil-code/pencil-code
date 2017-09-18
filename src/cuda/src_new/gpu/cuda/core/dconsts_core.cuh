/*
*   GPU constants shared among all different implementations.
*   In order to avoid code duplication, DCONST_EXTERN macro is used to define
*   the DCONST_EXTERN keyword for constants used in all other files but the one
*   which does the actual definition.
*/
#pragma once
#include "common/datatypes.h"

#ifdef INCLUDED_FROM_DCONST_DEFINER
    #define DCONST_EXTERN
#else
    #define DCONST_EXTERN extern
#endif

//Dimensions
DCONST_EXTERN __constant__ int d_mx, d_my, d_mz; // Grid dims including the boundaries
DCONST_EXTERN __constant__ int d_nx, d_ny, d_nz; // Dimensions of the computational domain
DCONST_EXTERN __constant__ int d_nx_min, d_ny_min, d_nz_min; // Starting indices of the computational domain
DCONST_EXTERN __constant__ int d_nx_max, d_ny_max, d_nz_max; // Ending indices (exclusive) of the computational domain

//Distances (in real units) between gridpoints
DCONST_EXTERN __constant__ real d_DSX, d_DSY, d_DSZ;
//Offset for the distance when using multiple GPUs
DCONST_EXTERN __constant__ real d_DSX_OFFSET, d_DSY_OFFSET, d_DSZ_OFFSET;

//Coordinate axis origin locations
DCONST_EXTERN __constant__ real d_XORIG;
DCONST_EXTERN __constant__ real d_YORIG;
DCONST_EXTERN __constant__ real d_ZORIG;

//Constant divisions in differentials in coeff form
DCONST_EXTERN __constant__ real d_DIFF1_DX_DIV, d_DIFF1_DY_DIV, d_DIFF1_DZ_DIV; //Constants for diff operations
DCONST_EXTERN __constant__ real d_DIFF2_DX_DIV, d_DIFF2_DY_DIV, d_DIFF2_DZ_DIV;
DCONST_EXTERN __constant__ real d_DIFFMN_DXDY_DIV, d_DIFFMN_DYDZ_DIV, d_DIFFMN_DXDZ_DIV; 

//Viscosity and the speed of sound
DCONST_EXTERN __constant__ real d_NU_VISC; //Viscosity
DCONST_EXTERN __constant__ real d_CS2_SOUND; //Speed of sound


//Forcing related 
DCONST_EXTERN __constant__ int d_FORCING_ENABLED;
DCONST_EXTERN __constant__ real d_PHI;
DCONST_EXTERN __constant__ real d_KK_VEC_X;
DCONST_EXTERN __constant__ real d_KK_VEC_Y;
DCONST_EXTERN __constant__ real d_KK_VEC_Z;
DCONST_EXTERN __constant__ real d_FORCING_KK_PART_X;
DCONST_EXTERN __constant__ real d_FORCING_KK_PART_Y;
DCONST_EXTERN __constant__ real d_FORCING_KK_PART_Z;

//Induction
DCONST_EXTERN __constant__ real d_ETA; //Magnetic diffusivity







































