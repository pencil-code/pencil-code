
#pragma once

// Declare constants used by device.

//Dimensions
EXTERN __constant__ int d_NX, d_NY, d_NZ; //Grid dims

//top and bottom indices of the computational domain (bottom inclusive, top exclusive)
EXTERN __constant__ int d_CX_BOT, d_CY_BOT, d_CZ_BOT;
EXTERN __constant__ int d_CX_TOP, d_CY_TOP, d_CZ_TOP;

//Computational domain's dims
EXTERN __constant__ int d_COMP_DOMAIN_SIZE_X, d_COMP_DOMAIN_SIZE_Y, d_COMP_DOMAIN_SIZE_Z;
EXTERN __constant__ int d_BOUND_SIZE; //Boundary zone size
EXTERN __constant__ float d_NELEMENTS_FLOAT; //Number of grid point elements for averaging. 
EXTERN __constant__ float d_DOMAIN_SIZE_X, d_DOMAIN_SIZE_Y, d_DOMAIN_SIZE_Z;

//Grid offsets
//eg. coordinates in the intermediate array are represented as:
// i + j*d_W_GRID_Y_OFFSET + k*d_W_GRID_Z_OFFSET
//
//and coordinates in the other arrays (which have paddings and whatnot):
// i + j*d_GRID_Y_OFFSET + k*d_GRID_Z_OFFSET
EXTERN __constant__ int d_W_GRID_Y_OFFSET;
EXTERN __constant__ int d_W_GRID_Z_OFFSET;
EXTERN __constant__ int d_GRID_Y_OFFSET;
EXTERN __constant__ int d_GRID_Z_OFFSET;

//Distances (in real units) between gridpoints
EXTERN __constant__ float d_DX, d_DY, d_DZ;

//Coordinate axis origin locations
EXTERN __constant__ float d_XORIG;
EXTERN __constant__ float d_YORIG;
EXTERN __constant__ float d_ZORIG;

//Shearing parameters
EXTERN __constant__ float d_Q_SHEAR;
EXTERN __constant__ float d_OMEGA;
EXTERN __constant__ float d_DELTA_Y;
EXTERN __constant__ int d_INTERP_ORDER;

//Logic switches for optional features
EXTERN __constant__ int d_LFORCING;
EXTERN __constant__ int d_LSHEAR;
EXTERN __constant__ int d_LCORIOLIS;
EXTERN __constant__ int d_FTRIGGER;

//Timestep
EXTERN __constant__ float d_DT;

//Forcing related 
EXTERN __constant__ float d_PHI;
EXTERN __constant__ float d_KK_VEC_X;
EXTERN __constant__ float d_KK_VEC_Y;
EXTERN __constant__ float d_KK_VEC_Z;
EXTERN __constant__ float d_FORCING_KK_PART_X;
EXTERN __constant__ float d_FORCING_KK_PART_Y;
EXTERN __constant__ float d_FORCING_KK_PART_Z;

//Viscosity and soundspeed
EXTERN __constant__ float d_NU_VISC;
EXTERN __constant__ float d_CS2_SOUND;

//Constants for RK3
EXTERN __constant__ float d_ALPHA1, d_ALPHA2, d_ALPHA3;//0.0; -0.53125; -1.1851851851851851;
EXTERN __constant__ float d_BETA1, d_BETA2, d_BETA3;//0.25; 0.88888888888888884; 0.75;

//Constants for derivatives
EXTERN __constant__ float d_FLT_9;
EXTERN __constant__ float d_FLT_45; 
EXTERN __constant__ float d_FLT_60; 

EXTERN __constant__ float d_FLT_2; 
EXTERN __constant__ float d_FLT_27; 
EXTERN __constant__ float d_FLT_270; 
EXTERN __constant__ float d_FLT_490; 
EXTERN __constant__ float d_FLT_180; 

//Constant divisions in differentials in coeff form
EXTERN __constant__ float d_DIFF1_DX_DIV;
EXTERN __constant__ float d_DIFF1_DY_DIV;
EXTERN __constant__ float d_DIFF1_DZ_DIV;

EXTERN __constant__ float d_DIFF2_DX_DIV;
EXTERN __constant__ float d_DIFF2_DY_DIV;
EXTERN __constant__ float d_DIFF2_DZ_DIV;

EXTERN __constant__ float d_DIFFMN_DXDY_DIV; 
EXTERN __constant__ float d_DIFFMN_DYDZ_DIV;
EXTERN __constant__ float d_DIFFMN_DXDZ_DIV;

//Single precision bounds (also found in <cfloat>)
//EXTERN __constant__ float d_FLT_MAX;
//EXTERN __constant__ float d_FLT_MIN;

EXTERN __constant__ int d_halo_widths_x[3], d_halo_widths_y[3], d_halo_widths_z[3];
//btw <cmath> has isnan function etc
