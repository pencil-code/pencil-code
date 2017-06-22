
#pragma once
/*
* Declare constants used by device as extern.
*/

//Explicit definition, must declare with "extern" keyword in
//separately compiled modules

//Dimensions
extern __constant__ int d_NX, d_NY, d_NZ; //Grid dims

//top and bottom indices of the computational domain (bottom inclusive, top exclusive)
extern __constant__ int d_CX_BOT, d_CY_BOT, d_CZ_BOT;
extern __constant__ int d_CX_TOP, d_CY_TOP, d_CZ_TOP;

//Computational domain's dims
extern __constant__ int d_COMP_DOMAIN_SIZE_X, d_COMP_DOMAIN_SIZE_Y, d_COMP_DOMAIN_SIZE_Z;
extern __constant__ int d_PAD_SIZE; //Pad size on X-axis (includes boundary zone)
extern __constant__ int d_BOUND_SIZE; //Boundary zone size
extern __constant__ float d_NELEMENTS_FLOAT; //Number of grid point elements for averaging. 
extern __constant__ float d_DOMAIN_SIZE_X, d_DOMAIN_SIZE_Y, d_DOMAIN_SIZE_Z;

//Grid offsets
//eg. coordinates in the intermediate array are represented as:
// i + j*d_W_GRID_Y_OFFSET + k*d_W_GRID_Z_OFFSET
//
//and coordinates in the other arrays (which have paddings and whatnot):
// i + j*d_GRID_Y_OFFSET + k*d_GRID_Z_OFFSET
extern __constant__ int d_W_GRID_Y_OFFSET;
extern __constant__ int d_W_GRID_Z_OFFSET;
extern __constant__ int d_GRID_Y_OFFSET;
extern __constant__ int d_GRID_Z_OFFSET;

//Distances (in real units) between gridpoints
extern __constant__ float d_DX, d_DY, d_DZ;

//Coordinate axis origin locations
extern __constant__ float d_XORIG;
extern __constant__ float d_YORIG;
extern __constant__ float d_ZORIG;

//Shearing parameters
extern __constant__ float d_Q_SHEAR;
extern __constant__ float d_OMEGA;
extern __constant__ float d_DELTA_Y;
extern __constant__ int d_INTERP_ORDER;

//Logic switches for optional features
extern __constant__ int d_LFORCING;
extern __constant__ int d_LSHEAR;
extern __constant__ int d_LCORIOLIS;
extern __constant__ int d_FTRIGGER;

//Timestep
extern __constant__ float d_DT;

//Forcing related 
extern __constant__ float d_PHI;
extern __constant__ float d_KK_VEC_X;
extern __constant__ float d_KK_VEC_Y;
extern __constant__ float d_KK_VEC_Z;
extern __constant__ float d_FORCING_KK_PART_X;
extern __constant__ float d_FORCING_KK_PART_Y;
extern __constant__ float d_FORCING_KK_PART_Z;

//Viscosity and soundspeed
extern __constant__ float d_NU_VISC;
extern __constant__ float d_CS2_SOUND;

//Constants for RK3
extern __constant__ float d_ALPHA1, d_ALPHA2, d_ALPHA3;//0.0; -0.53125; -1.1851851851851851;
extern __constant__ float d_BETA1, d_BETA2, d_BETA3;//0.25; 0.88888888888888884; 0.75;

//Constants for derivatives
extern __constant__ float d_FLT_9;
extern __constant__ float d_FLT_45; 
extern __constant__ float d_FLT_60; 

extern __constant__ float d_FLT_2; 
extern __constant__ float d_FLT_27; 
extern __constant__ float d_FLT_270; 
extern __constant__ float d_FLT_490; 
extern __constant__ float d_FLT_180; 

//Constant divisions in differentials in coeff form
extern __constant__ float d_DIFF1_DX_DIV;
extern __constant__ float d_DIFF1_DY_DIV;
extern __constant__ float d_DIFF1_DZ_DIV;

extern __constant__ float d_DIFF2_DX_DIV;
extern __constant__ float d_DIFF2_DY_DIV;
extern __constant__ float d_DIFF2_DZ_DIV;

extern __constant__ float d_DIFFMN_DXDY_DIV; 
extern __constant__ float d_DIFFMN_DYDZ_DIV;
extern __constant__ float d_DIFFMN_DXDZ_DIV;

//Single precision bounds (also found in <cfloat>)
//extern __constant__ float d_FLT_MAX;
//extern __constant__ float d_FLT_MIN;


//btw <cmath> has isnan function etc
