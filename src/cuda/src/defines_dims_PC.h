#pragma once

//#define NPAD 32-nghost
#define NPAD 0
#define BOUND_SIZE nghost
#define PAD_SIZE nghost+NPAD

#define COMP_DOMAIN_SIZE_X nx
#define COMP_DOMAIN_SIZE_Y ny
#define COMP_DOMAIN_SIZE_Z nz

#define NX mx + 2*NPAD
#define NY my
#define NZ mz

#define W_GRID_SIZE nw
#define GRID_SIZE   mw

#define CX_BOT l1-1+NPAD
#define CY_BOT m1-1
#define CZ_BOT n1-1

#define CX_TOP l2-1+NPAD
#define CY_TOP m2-1
#define CZ_TOP n2-1

