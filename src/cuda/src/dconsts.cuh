
#pragma once
/*
* Define constants used by device
*/

//Explicit definition, must declare with "extern" keyword in
//separately compiled modules

//Dimensions
__constant__ int d_NX, d_NY, d_NZ; //Grid dims

//top and bottom indices of the computational domain (bottom inclusive, top exclusive)
__constant__ int d_CX_BOT, d_CY_BOT, d_CZ_BOT;
__constant__ int d_CX_TOP, d_CY_TOP, d_CZ_TOP;

//Computational domain's dims
__constant__ int d_COMP_DOMAIN_SIZE_X, d_COMP_DOMAIN_SIZE_Y, d_COMP_DOMAIN_SIZE_Z;
__constant__ int d_PAD_SIZE; //Pad size on X-axis (includes boundary zone)
__constant__ int d_BOUND_SIZE; //Boundary zone size
__constant__ float d_NELEMENTS_FLOAT; //Number of grid point elements for averaging.
__constant__ float d_DOMAIN_SIZE_X, d_DOMAIN_SIZE_Y, d_DOMAIN_SIZE_Z;

//Grid offsets
//eg. coordinates in the intermediate array are represented as:
// i + j*d_W_GRID_Y_OFFSET + k*d_W_GRID_Z_OFFSET
//
//and coordinates in the other arrays (which have paddings and whatnot):
// i + j*d_GRID_Y_OFFSET + k*d_GRID_Z_OFFSET
__constant__ int d_W_GRID_Y_OFFSET;
__constant__ int d_W_GRID_Z_OFFSET;
__constant__ int d_GRID_Y_OFFSET;
__constant__ int d_GRID_Z_OFFSET;

//Distances (in real units) between gridpoints
__constant__ float d_DX, d_DY, d_DZ;

//Coordinate axis origin locations
__constant__ float d_XORIG;
__constant__ float d_YORIG;
__constant__ float d_ZORIG;

//Shearing parameters
__constant__ float d_Q_SHEAR;
__constant__ float d_OMEGA;
__constant__ float d_DELTA_Y;
__constant__ int d_INTERP_ORDER;

//Logic switches for optional features
__constant__ int d_LFORCING;
__constant__ int d_LSHEAR;
__constant__ int d_LCORIOLIS;
__constant__ int d_FTRIGGER;

//Timestep
__constant__ float d_DT;

//Forcing related 
__constant__ float d_PHI;
__constant__ float d_KK_VEC_X;
__constant__ float d_KK_VEC_Y;
__constant__ float d_KK_VEC_Z;
__constant__ float d_FORCING_KK_PART_X;
__constant__ float d_FORCING_KK_PART_Y;
__constant__ float d_FORCING_KK_PART_Z;

//Viscosity and soundspeed 
__constant__ float d_NU_VISC;
__constant__ float d_CS2_SOUND;

//Constants for RK3
__constant__ float d_ALPHA1, d_ALPHA2, d_ALPHA3;//0.0; -0.53125; -1.1851851851851851;
__constant__ float d_BETA1, d_BETA2, d_BETA3;//0.25; 0.88888888888888884; 0.75;

//Constants for derivatives
__constant__ float d_FLT_9;
__constant__ float d_FLT_45; 
__constant__ float d_FLT_60; 

__constant__ float d_FLT_2; 
__constant__ float d_FLT_27; 
__constant__ float d_FLT_270; 
__constant__ float d_FLT_490; 
__constant__ float d_FLT_180; 

//Constant divisions in differentials in coeff form
__constant__ float d_DIFF1_DX_DIV;
__constant__ float d_DIFF1_DY_DIV;
__constant__ float d_DIFF1_DZ_DIV;

__constant__ float d_DIFF2_DX_DIV;
__constant__ float d_DIFF2_DY_DIV;
__constant__ float d_DIFF2_DZ_DIV;

__constant__ float d_DIFFMN_DXDY_DIV; 
__constant__ float d_DIFFMN_DYDZ_DIV;
__constant__ float d_DIFFMN_DXDZ_DIV;

//Single precision bounds (also found in <cfloat>)
//__constant__ float d_FLT_MAX;
//__constant__ float d_FLT_MIN;


//btw <cmath> has isnan function etc
