#include <stdio.h>

#define PAD_SIZE 32
#define BOUND_SIZE 3 //Boundary zone size in y&z axis

#define COMP_DOMAIN_SIZE_X 128
#define COMP_DOMAIN_SIZE_Y 128
#define COMP_DOMAIN_SIZE_Z 128

#define NX (COMP_DOMAIN_SIZE_X + 2*PAD_SIZE)
#define NY (COMP_DOMAIN_SIZE_Y + 2*BOUND_SIZE)
#define NZ (COMP_DOMAIN_SIZE_Z + 2*BOUND_SIZE)

#define GRID_SIZE NX*NY*NZ //Total grid size including the padding zones

int main(int argc, char** argv)
{
	//----------------------------------------------------------
	// Allocate host memory
	//----------------------------------------------------------

	printf("Hello: GRID_SIZE = %i \n", GRID_SIZE);

	// No need to malloc as dims are known at compile time
	float lnrho[GRID_SIZE]; //Log density
	float uu_x[GRID_SIZE], uu_y[GRID_SIZE], uu_z[GRID_SIZE]; //velocities

	printf("Hello! \n");

	return EXIT_SUCCESS;
}	
