#include <stdio.h>
#include <stdlib.h>
#include "gpu_astaroth.cuh"
int main(){
	
	float uu_x[128];
	float uu_y[128];
	float uu_z[128];
	float lnrho[128];
	finalizeGpu(uu_x, uu_y, uu_z, lnrho);
return 0;
}
