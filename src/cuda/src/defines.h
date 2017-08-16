#pragma once

/*
* Definitions in addition to the one's defined in .conf files
*/
#define BOUND_SIZE 3 //Boundary zone size in y&z axis
#define PAD_SIZE 32 //Pad size in x axis (total padding 2*PAD_SIZE), includes the boundary zone

#define NX (COMP_DOMAIN_SIZE_X + 2*PAD_SIZE)
#define NY (COMP_DOMAIN_SIZE_Y + 2*BOUND_SIZE)
#define NZ (COMP_DOMAIN_SIZE_Z + 2*BOUND_SIZE)

//Grid size of the temporary arrays that contain only the computational domain
#define W_GRID_SIZE COMP_DOMAIN_SIZE_X*COMP_DOMAIN_SIZE_Y*COMP_DOMAIN_SIZE_Z
#define GRID_SIZE NX*NY*NZ //Total grid size including the padding zones

//Top indices of the computational domain (!!!EXCLUSIVE!!!)
#define CX_TOP (NX - PAD_SIZE)
#define CY_TOP (NY - BOUND_SIZE)
#define CZ_TOP (NZ - BOUND_SIZE)

//Bottom indices of the computational domain (!!!INCLUSIVE!!!)
#define CX_BOT (PAD_SIZE)
#define CY_BOT (BOUND_SIZE)
#define CZ_BOT (BOUND_SIZE)

//Distances between gridpoints
#define DX ( DOMAIN_SIZE_X / (float)(CX_TOP-CX_BOT))
#define DY ( DOMAIN_SIZE_Y / (float)(CY_TOP-CY_BOT))
#define DZ ( DOMAIN_SIZE_Z / (float)(CZ_TOP-CZ_BOT))

//Minimum grid spacing (DSMIN used in timestep_cuda) cannot be defined
//here as it requires floating point comparison (DX < DY etc). 

//Data paths
/*
const char* DATA_LNRHO_PATH = "data/density.dat";
const char* DATA_UU_X_PATH = "data/velx.dat";
const char* DATA_UU_Y_PATH = "data/vely.dat";
const char* DATA_UU_Z_PATH = "data/velz.dat";

const char* DATA_LNRHO_PATH = "data/lnrho.dat";
const char* DATA_UU_X_PATH = "data/uu_x.dat";
const char* DATA_UU_Y_PATH = "data/uu_y.dat";
const char* DATA_UU_Z_PATH = "data/uu_z.dat";
*/
