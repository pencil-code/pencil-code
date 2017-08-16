//C libraries
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

//CUDA libraries
//#include <cfloat>

//Headers
#include "dconsts.cuh"
#include "integrators.cuh"
#include "boundcond.cuh"
#include "timestep.cuh"
#include "defines.h"
#include "io.h"
#include "slice.cuh"
#include "smem.cuh"
#include "forcing.cuh"


//DEBUG
#include "initutils.h"
#include "diagnostics.cuh"
#define dbug 0

//Swaps pointers a,b (do xor swap if too slow)
inline void swap_ptrs(float** a, float** b)
{
	float* temp = *a;
	*a = *b;
	*b = temp;
}


//Loads constants into device memory
//const int nx = NX, ny = NY, nz = NZ;
//checkErr( cudaMemcpyToSymbol(d_NX, &nx, sizeof(int)) );
void load_dconsts()
{
	//---------Grid dims-----------------------
	const int nx = NX, ny = NY, nz = NZ;
	const int pad_size = PAD_SIZE, bound_size = BOUND_SIZE;
	const int 	comp_domain_size_x = COMP_DOMAIN_SIZE_X, 
			comp_domain_size_y = COMP_DOMAIN_SIZE_Y, 
			comp_domain_size_z = COMP_DOMAIN_SIZE_Z;
	//Needed for calculating averages
	const float nelements_float = COMP_DOMAIN_SIZE_X*COMP_DOMAIN_SIZE_Y*COMP_DOMAIN_SIZE_Z;

	const float 	domain_size_x = DOMAIN_SIZE_X, 
			domain_size_y = DOMAIN_SIZE_Y, 
			domain_size_z = DOMAIN_SIZE_Z;

	checkErr( cudaMemcpyToSymbol(d_NX, &nx, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_NY, &ny, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_NZ, &nz, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_PAD_SIZE, &pad_size, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_BOUND_SIZE, &bound_size, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_COMP_DOMAIN_SIZE_X, &comp_domain_size_x, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_COMP_DOMAIN_SIZE_Y, &comp_domain_size_y, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_COMP_DOMAIN_SIZE_Z, &comp_domain_size_z, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_NELEMENTS_FLOAT, &nelements_float, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_DOMAIN_SIZE_X, &domain_size_x, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DOMAIN_SIZE_Y, &domain_size_y, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DOMAIN_SIZE_Z, &domain_size_z, sizeof(float)) );

	const int h_w_grid_y_offset = COMP_DOMAIN_SIZE_X;
	const int h_w_grid_z_offset = COMP_DOMAIN_SIZE_X*COMP_DOMAIN_SIZE_Y;
	const int h_grid_y_offset = NX;
	const int h_grid_z_offset = NX*NY;

	checkErr( cudaMemcpyToSymbol(d_W_GRID_Y_OFFSET, &h_w_grid_y_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_W_GRID_Z_OFFSET, &h_w_grid_z_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_GRID_Y_OFFSET, &h_grid_y_offset, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_GRID_Z_OFFSET, &h_grid_z_offset, sizeof(int)) );
	//-------------------------------------------

	//------Computational domain's bottom and top indices---------
	const int cx_top = CX_TOP, cy_top = CY_TOP, cz_top = CZ_TOP;
	const int cx_bot = CX_BOT, cy_bot = CY_BOT, cz_bot = CY_BOT;

	checkErr( cudaMemcpyToSymbol(d_CX_TOP, &cx_top, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CY_TOP, &cy_top, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CZ_TOP, &cz_top, sizeof(int)) );

	checkErr( cudaMemcpyToSymbol(d_CX_BOT, &cx_bot, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CY_BOT, &cy_bot, sizeof(int)) );
	checkErr( cudaMemcpyToSymbol(d_CZ_BOT, &cz_bot, sizeof(int)) );
	//-------------------------------------------

	//-------------Real distance between grid points---------------
	const float dx = DX, dy = DY, dz = DZ; //MV: Was int. Might have caused some problems?
	checkErr( cudaMemcpyToSymbol(d_DX, &dx, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DY, &dy, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DZ, &dz, sizeof(float)) );
	//-------------------------------------------

	//----------Location of the grid coordinate origin-------------
	const float xorig = XORIG, yorig = YORIG, zorig = ZORIG;
        checkErr( cudaMemcpyToSymbol(d_XORIG, &xorig, sizeof(float)) );
        checkErr( cudaMemcpyToSymbol(d_YORIG, &yorig, sizeof(float)) );
        checkErr( cudaMemcpyToSymbol(d_ZORIG, &zorig, sizeof(float)) );
	//-------------------------------------------

	//----------Shearing parameters---------------
	const float q_shear = Q_SHEAR;
	const float omega = OMEGA; 
	const int interp_order = INTERP_ORDER;

	printf("compute INTERP_ORDER = %i \n", INTERP_ORDER);
	printf("compute interp_order = %i \n", interp_order); 

        checkErr( cudaMemcpyToSymbol(d_INTERP_ORDER, &interp_order, sizeof(int)) );
        checkErr( cudaMemcpyToSymbol(d_Q_SHEAR, &q_shear, sizeof(float)) );
        checkErr( cudaMemcpyToSymbol(d_OMEGA, &omega, sizeof(float)) );
	//-------------------------------------------

	//------------------Optional physics switches------------------------
	const int lforcing = LFORCING;
	const int lshear = LSHEAR;
	const int lcoriolis = LCORIOLIS;
        checkErr( cudaMemcpyToSymbol(d_LFORCING, &lforcing, sizeof(int)) );
        checkErr( cudaMemcpyToSymbol(d_LSHEAR, &lshear, sizeof(int)) );
        checkErr( cudaMemcpyToSymbol(d_LCORIOLIS, &lcoriolis, sizeof(int)) );
	//-------------------------------------------------------------------

	//------------Random constants for computation-----------------
	const float h_ALPHA1 = 0.0; 
	const float h_ALPHA2 = -0.53125; 
	const float h_ALPHA3 = -1.1851851851851851;
	
	checkErr( cudaMemcpyToSymbol(d_ALPHA1, &h_ALPHA1, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_ALPHA2, &h_ALPHA2, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_ALPHA3, &h_ALPHA3, sizeof(float)) );

	const float h_BETA1 = 0.25; 
	const float h_BETA2 = 0.88888888888888884;
	const float h_BETA3 = 0.75;

	checkErr( cudaMemcpyToSymbol(d_BETA1, &h_BETA1, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_BETA2, &h_BETA2, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_BETA3, &h_BETA3, sizeof(float)) );

	const float nu_visc = NU_VISC;
	const float cs2_sound = pow(CS_SOUND, 2.0);
	checkErr( cudaMemcpyToSymbol(d_NU_VISC, &nu_visc, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_CS2_SOUND, &cs2_sound, sizeof(float)) );	
	//-------------------------------------------

	//------------Constants for derivatives-----------------
	//JP: For some reason i have a very weird feeling when typing these up,
	//assigning numbers... into variables... that are named like numbers...duude...
	//MV: You can learn all about being the Dude in http://dudeism.com/ ^_^
	//JP: Nice :D
	const float flt_9 = 9.0;
	const float flt_45 = 45.0; 
	const float flt_60 = 60.0; 

	const float flt_2 = 2.0; 
	const float flt_27 = 27.0; 
	const float flt_270 = 270.0; 
	const float flt_490 = 490.0; 
	const float flt_180 = 180.0; 

	//MV: Added these to avoid unneseccary division in differentials
	const float diff1_dx = 1.0/(60.0*DX);
	const float diff1_dy = 1.0/(60.0*DY);
	const float diff1_dz = 1.0/(60.0*DZ);

	const float diff2_dx = 1.0/(180.0*DX*DX);
	const float diff2_dy = 1.0/(180.0*DY*DY);
	const float diff2_dz = 1.0/(180.0*DZ*DZ);

	const float diffmn_dxdy = (1.0/720.0)*(1.0/DX)*(1.0/DY); 
	const float diffmn_dydz = (1.0/720.0)*(1.0/DY)*(1.0/DZ);
	const float diffmn_dxdz = (1.0/720.0)*(1.0/DZ)*(1.0/DX);
	
	checkErr( cudaMemcpyToSymbol(d_FLT_9, &flt_9, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_45, &flt_45, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_60, &flt_60, sizeof(float)) );

	checkErr( cudaMemcpyToSymbol(d_FLT_2, &flt_2, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_27, &flt_27, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_270, &flt_270, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_490, &flt_490, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_FLT_180, &flt_180, sizeof(float)) );

	checkErr( cudaMemcpyToSymbol(d_DIFF1_DX_DIV, &diff1_dx, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFF1_DY_DIV, &diff1_dy, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFF1_DZ_DIV, &diff1_dz, sizeof(float)) );

	checkErr( cudaMemcpyToSymbol(d_DIFF2_DX_DIV, &diff2_dx, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFF2_DY_DIV, &diff2_dy, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFF2_DZ_DIV, &diff2_dz, sizeof(float)) );

	checkErr( cudaMemcpyToSymbol(d_DIFFMN_DXDY_DIV, &diffmn_dxdy, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFFMN_DYDZ_DIV, &diffmn_dydz, sizeof(float)) );
	checkErr( cudaMemcpyToSymbol(d_DIFFMN_DXDZ_DIV, &diffmn_dxdz, sizeof(float)) );
	//-------------------------------------------
}


int main(int argc, char** argv)
{
	//http://devblogs.nvidia.com/parallelforall/cuda-pro-tip-always-set-current-device-avoid-multithreading-bugs/
	//Hmmmm....
	int device;
	cudaGetDevice(&device);
	printf("Using device %d\n", device);
	//cudaSetDevice(device); //Not yet enabled

	//Ensure that we're using a clean device
	cudaDeviceReset();


	//----------------------------------------------------------
	// Initialize global host variables
	//----------------------------------------------------------
	print_init_config();
	print_run_config();
	print_additional_defines();
	//----------------------------------------------------------

	//----------------------------------------------------------
	// Allocate host memory
	//----------------------------------------------------------

	// No need to malloc as dims are known at compile time
	// MV: Requires malloc because otherwise causes stack overflow on some other machines!
	//float lnrho[GRID_SIZE]; //Log density
	//float uu_x[GRID_SIZE], uu_y[GRID_SIZE], uu_z[GRID_SIZE]; //velocities
	// MV: Allocating solution. 
        float *lnrho; //Log density
        float *uu_x, *uu_y, *uu_z; //velocities

        lnrho = (float*) malloc(sizeof(float)*GRID_SIZE);
        uu_x  = (float*) malloc(sizeof(float)*GRID_SIZE);
        uu_y  = (float*) malloc(sizeof(float)*GRID_SIZE);
        uu_z  = (float*) malloc(sizeof(float)*GRID_SIZE);

	//Format the grids into 0.0 values to avoid potential noise in memory
	set_grids_zero(lnrho, uu_x, uu_y, uu_z);


	//----------------------------------------------------------


	//----------------------------------------------------------
	// Load existing data into host memory or initialize the grid
	//----------------------------------------------------------
	//MV: We might want to have the data not under the the src folder 
	//MV: but instead like: ../data/density.dat (meaning run_dir/data/density.dat)
	//MV: But this is not a pririty.
	//JP: Yeah, 3 possible solutions for this: symbolic link in ac_build script from src/data to /data, user defines the path in init.conf or then some wildcard solution
	//TODO.txt: #1
	const char* DATA_LNRHO_PATH = "data/lnrho.dat";
	const char* DATA_UU_X_PATH = "data/uu_x.dat";
	const char* DATA_UU_Y_PATH = "data/uu_y.dat";
	const char* DATA_UU_Z_PATH = "data/uu_z.dat";

	//If the first data file loads ok, load the rest (lol not really a foolproof solution)
	//TODO.txt: #2
	//----Note! Never enters this if!---
   	if (0 && load_grid_data(lnrho, DATA_LNRHO_PATH)) {
   		load_grid_data(uu_x, DATA_UU_X_PATH);
		load_grid_data(uu_y, DATA_UU_Y_PATH);
		load_grid_data(uu_z, DATA_UU_Z_PATH);
		printf("Data load successful.\n");
	}
	else // Initialize the grid 
	{
		lnrho_init(lnrho);
		hydro_init(uu_x, uu_y, uu_z);

		save_grid_information(0.0); // Save grid information for the most current .dat file TODO: test. 

		save_grid_data(lnrho, DATA_LNRHO_PATH);
		save_grid_data(uu_x, DATA_UU_X_PATH);
		save_grid_data(uu_y, DATA_UU_Y_PATH);
		save_grid_data(uu_z, DATA_UU_Z_PATH);

		save_grid_data(lnrho, (char*) "data/lnrho0.dat"); // The initial state for comparison
  	 	save_grid_data(uu_x, (char*) "data/uu_x0.dat");
   		save_grid_data(uu_y, (char*) "data/uu_y0.dat");
   		save_grid_data(uu_z, (char*) "data/uu_z0.dat");
		printf("Grid initialization successful.\n");
	}
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if(dbug){
		for (int k=0; k < NZ; k++) {
			for (int j=0; j < NY; j++) {
				for (int i=0; i < NX; i++) {
					//int idx = i + j*NX + k*NX*NY;
					//lnrho[idx] = idx;
					//uu_x[idx] = idx;
					//uu_y[idx] = idx;
					//uu_z[idx] = idx;
				}
			}
		}	
	
		FILE *fp1;
		fp1 = fopen("uux_input.txt", "w");
		for (int k= 0; k < 1; k++) {
			fprintf(fp1, "Zplane %d\n",k);
			for (int j=0; j < NY; j++) {
				fprintf(fp1, "Row %d\n",j);
				for (int i=0; i < NX; i++) {
					fprintf(fp1, "(%d, %f), idx = %d  ",i, uu_x[i + j*NX + k*NX*NY], i + j*NX + k*NX*NY);	
				}
			fprintf(fp1, "\n---------------\n");
			}
		}
	}

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//----------------------------------------------------------
	

	//DEBUG DEBUG DEBUG
	//TODO: #3
	//debug_clear_grid(lnrho, uu_x, uu_y, uu_z);
	//debug_grid_init(lnrho, uu_x, uu_y, uu_z, 5);

	//----------------------------------------------------------
	// Allocate pinned host memory and get its device pointers
	//----------------------------------------------------------
	//TODO Slice axis should be determined from run.conf
	const char slice_axis = 'z'; //Determine slice axis used for .ani slices
   	int slice_size;
   		switch(slice_axis) 
		{
			case 'x':
				slice_size = 	(COMP_DOMAIN_SIZE_Y+2*BOUND_SIZE) *
						(COMP_DOMAIN_SIZE_Z+2*BOUND_SIZE);
				break;
			case 'y':
				slice_size = 	(COMP_DOMAIN_SIZE_X+2*BOUND_SIZE) *
						(COMP_DOMAIN_SIZE_Z+2*BOUND_SIZE);
				break;
			case 'z':
				slice_size = 	(COMP_DOMAIN_SIZE_X+2*BOUND_SIZE) *
						(COMP_DOMAIN_SIZE_Y+2*BOUND_SIZE);
				break;
			default:
				printf("INVALID slice axis in compute.cu!\n");
				exit(EXIT_FAILURE);
		}
	
	float *slice_lnrho, *slice_uu, *slice_uu_x, *slice_uu_y, *slice_uu_z;
	checkErr( cudaHostAlloc((float**)&slice_lnrho, sizeof(float)*slice_size, cudaHostAllocMapped) );
	checkErr( cudaHostAlloc((float**)&slice_uu, sizeof(float)*slice_size, cudaHostAllocMapped) );
	checkErr( cudaHostAlloc((float**)&slice_uu_x, sizeof(float)*slice_size, cudaHostAllocMapped) );
	checkErr( cudaHostAlloc((float**)&slice_uu_y, sizeof(float)*slice_size, cudaHostAllocMapped) );
	checkErr( cudaHostAlloc((float**)&slice_uu_z, sizeof(float)*slice_size, cudaHostAllocMapped) );

	//Create device pointers for pinned host memory
	float *d_slice_lnrho, *d_slice_uu, *d_slice_uu_x, *d_slice_uu_y, *d_slice_uu_z;
	checkErr( cudaHostGetDevicePointer((float **)&d_slice_lnrho, (float *)slice_lnrho, 0) );
	checkErr( cudaHostGetDevicePointer((float **)&d_slice_uu, (float *)slice_uu, 0) );
	checkErr( cudaHostGetDevicePointer((float **)&d_slice_uu_x, (float *)slice_uu_x, 0) );
	checkErr( cudaHostGetDevicePointer((float **)&d_slice_uu_y, (float *)slice_uu_y, 0) );
	checkErr( cudaHostGetDevicePointer((float **)&d_slice_uu_z, (float *)slice_uu_z, 0) );
	//----------------------------------------------------------
	
	//----------------------------------------------------------
	// Allocate device memory
	//----------------------------------------------------------
	float *d_lnrho, *d_uu_x, *d_uu_y, *d_uu_z;
	checkErr( cudaMalloc(&d_lnrho, sizeof(float)*GRID_SIZE) );
	checkErr( cudaMalloc(&d_uu_x, sizeof(float)*GRID_SIZE) );
	checkErr( cudaMalloc(&d_uu_y, sizeof(float)*GRID_SIZE) );
	checkErr( cudaMalloc(&d_uu_z, sizeof(float)*GRID_SIZE) );

	float *d_w_lnrho, *d_w_uu_x, *d_w_uu_y, *d_w_uu_z;
	checkErr( cudaMalloc(&d_w_lnrho, sizeof(float)*W_GRID_SIZE) );
	checkErr( cudaMalloc(&d_w_uu_x, sizeof(float)*W_GRID_SIZE) );
	checkErr( cudaMalloc(&d_w_uu_y, sizeof(float)*W_GRID_SIZE) );
	checkErr( cudaMalloc(&d_w_uu_z, sizeof(float)*W_GRID_SIZE) );

	//Temporary arrays
	float *d_lnrho_dest, *d_uu_x_dest, *d_uu_y_dest, *d_uu_z_dest;
	checkErr( cudaMalloc(&d_lnrho_dest, sizeof(float)*GRID_SIZE) );
	checkErr( cudaMalloc(&d_uu_x_dest, sizeof(float)*GRID_SIZE) );
	checkErr( cudaMalloc(&d_uu_y_dest, sizeof(float)*GRID_SIZE) );
	checkErr( cudaMalloc(&d_uu_z_dest, sizeof(float)*GRID_SIZE) );

	float *d_div_uu; //Divergence field
	checkErr( cudaMalloc(&d_div_uu, sizeof(float)*GRID_SIZE) );	

	float* d_umax, *d_umin, *d_partial_result;//Device pointer for max vec and partial result for the reductions
        float* d_urms; //Device pointer for urms 
	float *d_uxrms, *d_uyrms, *d_uzrms; //Device pointer for uxrms, uyrms, uzrms 
	float *d_rhorms, *d_rhomax, *d_rhomin; //Device pointer for rhorms, rhomax, rhomin
	float *d_uxmax, *d_uymax, *d_uzmax; //Device pointer for uxmax, uymax, uzmax
	float *d_uxmin, *d_uymin, *d_uzmin; //Device pointer for uxmin, uymin, uzmin

	checkErr( cudaMalloc((float**) &d_umax, sizeof(float)) );   //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_umin, sizeof(float)) );   //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_urms, sizeof(float)) );   //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_uxrms, sizeof(float)) );  //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_uyrms, sizeof(float)) );  //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_uzrms, sizeof(float)) );  //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_rhorms, sizeof(float)) ); //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_rhomax, sizeof(float)) ); //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_rhomin, sizeof(float)) ); //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_uxmax, sizeof(float)) );   //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_uxmin, sizeof(float)) );   //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_uymax, sizeof(float)) );   //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_uymin, sizeof(float)) );   //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_uzmax, sizeof(float)) );   //TODO this somewhere else
        checkErr( cudaMalloc((float**) &d_uzmin, sizeof(float)) );   //TODO this somewhere else

	//Allocate memory for the partial result used for calculating max velocity in timestep_cuda
	dim3 blocksPerGrid;
	blocksPerGrid.x = ceil((float) COMP_DOMAIN_SIZE_X / (float)COL_THREADS_X);
	blocksPerGrid.y = ceil((float) COMP_DOMAIN_SIZE_Y / (float)COL_THREADS_Y);
	blocksPerGrid.z = ceil((float) COMP_DOMAIN_SIZE_Z / (float)COL_ELEMS_PER_THREAD);
	static const int BLOCKS_TOTAL = blocksPerGrid.x * blocksPerGrid.y * blocksPerGrid.z;
	checkErr( cudaMalloc((float**) &d_partial_result, sizeof(float)*BLOCKS_TOTAL) );
	//----------------------------------------------------------
	printf("Device mem allocated: %f MiB\n", (4*sizeof(float)*GRID_SIZE + 4*sizeof(float)*W_GRID_SIZE)/powf(2,20));
	printf("Main array (d_lnrho) dims: (%d,%d,%d)\ntemporary result array dims (d_w_lnrho etc)(%d,%d,%d)\n", NX,NY,NZ, COMP_DOMAIN_SIZE_X,COMP_DOMAIN_SIZE_Y,COMP_DOMAIN_SIZE_Z);

	//----------------------------------------------------------
	// Load data into device memory
	//----------------------------------------------------------
	checkErr( cudaMemcpy(d_lnrho, lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
	checkErr( cudaMemcpy(d_uu_x, uu_x, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
	checkErr( cudaMemcpy(d_uu_y, uu_y, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
	checkErr( cudaMemcpy(d_uu_z, uu_z, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
	//Init also the dest arrays to avoin roaming NaN values
        checkErr( cudaMemcpy(d_lnrho_dest, lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
        checkErr( cudaMemcpy(d_uu_x_dest, uu_x, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
        checkErr( cudaMemcpy(d_uu_y_dest, uu_y, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
        checkErr( cudaMemcpy(d_uu_z_dest, uu_z, sizeof(float)*GRID_SIZE, cudaMemcpyHostToDevice) );
	//----------------------------------------------------------
	/*FILE *fp;
	fp = fopen("refzplane0.txt", "w");
		for(int i = 0; i < NX; i++){
			//fprintf(fp, "NX = %d , NY = %d \n\n\n ", NX, NY );
			for(int j = 0; j < NY; j++){
				fprintf(fp, "%f | ", uu_x[(NX*i) + j]);
				//fprintf(fp, "%d | ", (NX*i) + j);
			}
				fprintf(fp, "\n\n");
		}
	fclose(fp);*/
	

	//----------------------------------------------------------
	//Load constants into device memory
	//----------------------------------------------------------
	load_dconsts();
	//----------------------------------------------------------

	

	//----------------------------------------------------------
	// Save the first slice
	//----------------------------------------------------------
	get_slice_cuda('z', d_slice_lnrho, d_slice_uu, d_slice_uu_x, d_slice_uu_y, d_slice_uu_z, d_lnrho, d_uu_x, d_uu_y, d_uu_z);
	//Sync
	cudaThreadSynchronize();
	//Save slice
	save_anim_slice(slice_axis, 0, slice_lnrho, slice_uu, slice_uu_x, slice_uu_y, slice_uu_z);

	//----------------------------------------------------------

	//----------------------------------------------------------
	// Init delta_y for shearing boundaries
	//----------------------------------------------------------

	float delta_y;
	delta_y = 0.0; //Assuming now that here t = 0 TODO: If t read from previous save, needs to be adapted
        checkErr( cudaMemcpyToSymbol(d_DELTA_Y, &delta_y, sizeof(float)) );


	//----------------------------------------------------------

	//----------------------------------------------------------
	// Compute boundary conditions and save to a slice for comparison
	//----------------------------------------------------------
	//Boundconds
	boundcond_cuda(d_lnrho, d_uu_x, d_uu_y, d_uu_z);
	boundcond_cuda(d_lnrho_dest, d_uu_x_dest, d_uu_y_dest, d_uu_z_dest);
	//Slice
	get_slice_cuda('z', d_slice_lnrho, d_slice_uu, d_slice_uu_x, d_slice_uu_y, d_slice_uu_z, d_lnrho, d_uu_x, d_uu_y, d_uu_z);
	//Sync
	cudaThreadSynchronize();
	//Save slice
	save_anim_slice(slice_axis, 1, slice_lnrho, slice_uu, slice_uu_x, slice_uu_y, slice_uu_z);

	//----------------------------------------------------------

	//Initialize the timestep file.
	init_ts_file();
	//init_used_k_file();

	//----------------------------------------------------------
	//Forcing stuff initialized here 
	//
	int nk;
	float *kk_x, *kk_y, *kk_z, kaver;
	float kk_vec_x, kk_vec_y, kk_vec_z, forcing_kk_part_x, forcing_kk_part_y, forcing_kk_part_z, phi;
	int ftrigger = 0; //A trigger which can turn forcing of at determined time. 
	checkErr( cudaMemcpyToSymbol(d_FTRIGGER, &ftrigger, sizeof(int)) );

	if (LFORCING) {
		//Determine the k-vector arrays size 
		nk = calc_k_vector_size();
		//Allocate vector arrays
		kk_x = (float*) malloc(sizeof(float)*nk);
		kk_y = (float*) malloc(sizeof(float)*nk);
		kk_z = (float*) malloc(sizeof(float)*nk);
		//Initialize the forcing k-value cloud.
		initialize_k_vectors(kk_x, kk_y, kk_z, &kaver, nk);
		init_used_k_file(kaver);
		//Set a random seed for forcing purposes
		srand(666);
	}

	//----------------------------------------------------------
	// Main loop
	//----------------------------------------------------------

	float dt; //Timestep
	float time_elapsed;//Real time
	float t = 0.0;	//Simulation time      TODO: Should be read from the most recent step, if exists
	int step = 0;	//Current step number  TODO: Saving slices should be based on a time interval. Should be read from the most recent step, if exists.
	int found_nan = 0; //Diagnostic variable for terminating the code if NaN has been found. 

	float t_anim; //Time when the next animation slice will be saved. 
	int t_anim_number; //The number of the animation slice file. 
	if (step > 0) { 
		//read_anim_state(&t_anim, &t_anim_number); //TODO: Include something like this into the code. 
	} else { 
		t_anim = ANIM_STEP; t_anim_number = 0; 
	}
	
	/*
	run_diagnostics(lnrho, uu_x, uu_y, uu_z,
			d_lnrho, d_uu_x, d_uu_y, d_uu_z, 
			d_w_lnrho, d_w_uu_x, d_w_uu_y, d_w_uu_z,
			d_lnrho_dest, d_uu_x_dest, d_uu_y_dest, d_uu_z_dest,
			d_div_uu,
			d_umax, d_partial_result);
	return;
	*/

	int steps = 0;
	float accumulated_time = 0.0f;

	//Create timer
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	//Looooooop
	//TODO: Max number of steps should be read from run.conf. Not like this. 
	dt = timestep_cuda(d_umax, d_partial_result, d_uu_x, d_uu_y, d_uu_z); //// omer remove this line after debug
	//while (t <= MAX_TIME && step < MAX_STEPS) {
	while (t <= 0*dt && step < MAX_STEPS) {
		//Check for NaN. TODO: make optional
		printf("Checking d_lnrho, d_uu_x, d_uu_y, d_uu_z...\n");
		found_nan = check_grid_for_nan(d_lnrho, d_uu_x, d_uu_y, d_uu_z);
                printf("Checking d_lnrho_dest, d_uu_x_dest, d_uu_y_dest, d_uu_z_dest...\n");
                found_nan = check_grid_for_nan(d_lnrho_dest, d_uu_x_dest, d_uu_y_dest, d_uu_z_dest);

		//Start timing
		cudaEventRecord( start, 0 );
	
		//Calculate timestep
		dt = timestep_cuda(d_umax, d_partial_result, d_uu_x, d_uu_y, d_uu_z);
		//Update dt in device's constant memory
		checkErr( cudaMemcpyToSymbol(d_DT, &dt, sizeof(float)) );
		printf("DT: %f\n", dt);
		//Update simulation time elapsed	
		t += dt;

		if (LSHEAR) {
			// Calculate the shearing shift delta_y for a shearing periodic x-boundary
			delta_y = Q_SHEAR*OMEGA*DOMAIN_SIZE_X*t; 
			delta_y -= floor(delta_y/DOMAIN_SIZE_Y)*DOMAIN_SIZE_Y;
			printf("delta_y = %f", delta_y);
			checkErr( cudaMemcpyToSymbol(d_DELTA_Y, &delta_y, sizeof(float)) );
		}
		
		if (t >= T_STOP_FORCING && ftrigger == 0) {
			//Stop forcing after declared time, byt turning on the trigger
			ftrigger = 1;
			checkErr( cudaMemcpyToSymbol(d_FTRIGGER, &ftrigger, sizeof(int)) );
		}

		if (LFORCING && ftrigger == 0) {
			//Randomly choose a forcing k-vector
			choose_random_k(&kk_vec_x, &kk_vec_y, &kk_vec_z, kk_x, kk_y, kk_z, nk); 
			phi = choose_random_phi();
			fkt_forcing_coefficient(&forcing_kk_part_x, &forcing_kk_part_y, &forcing_kk_part_z, kk_vec_x, kk_vec_y, kk_vec_z, dt);
			//Write the k-vertor stuff in a file for analysis
			save_used_k(kk_vec_x, kk_vec_y, kk_vec_z, phi, forcing_kk_part_x, forcing_kk_part_y, forcing_kk_part_z);
			//Copy forcing coefficients to the device's constant memory
			checkErr( cudaMemcpyToSymbol(d_KK_VEC_X, &kk_vec_x, sizeof(float)) );
			checkErr( cudaMemcpyToSymbol(d_KK_VEC_Y, &kk_vec_y, sizeof(float)) );
			checkErr( cudaMemcpyToSymbol(d_KK_VEC_Z, &kk_vec_z, sizeof(float)) );
			checkErr( cudaMemcpyToSymbol(d_FORCING_KK_PART_X, &forcing_kk_part_x, sizeof(float)) );
			checkErr( cudaMemcpyToSymbol(d_FORCING_KK_PART_Y, &forcing_kk_part_y, sizeof(float)) );
			checkErr( cudaMemcpyToSymbol(d_FORCING_KK_PART_Z, &forcing_kk_part_z, sizeof(float)) );
			checkErr( cudaMemcpyToSymbol(d_PHI, &phi, sizeof(float)) );
		}

		//Integrate
		rungekutta2N_cuda(d_lnrho, d_uu_x, d_uu_y, d_uu_z, 
				  d_w_lnrho, d_w_uu_x, d_w_uu_y, d_w_uu_z,
				  d_lnrho_dest, d_uu_x_dest, d_uu_y_dest, d_uu_z_dest,
				  d_div_uu);

		//Swap array pointers
		swap_ptrs(&d_lnrho, &d_lnrho_dest);
		swap_ptrs(&d_uu_x, &d_uu_x_dest);
		swap_ptrs(&d_uu_y, &d_uu_y_dest);
		swap_ptrs(&d_uu_z, &d_uu_z_dest);

		//TODO: Copy the device data to host, and write to a file, if at the chose step interval. 

		//Get time elapsed
		cudaEventRecord( stop, 0 );
		cudaEventSynchronize( stop );
		cudaEventElapsedTime( &time_elapsed, start, stop );
		printf("Time elapsed per cycle: \t%f ms\n", time_elapsed);
		printf("Simulation time elapsed: \t%f/%f s\n", t, MAX_TIME);

		steps++;
		accumulated_time += time_elapsed;
		
		//Save animation slice
		if (t >= t_anim) {
			//Transfer slice to slice arrays
			get_slice_cuda('z', d_slice_lnrho, d_slice_uu, d_slice_uu_x, d_slice_uu_y, d_slice_uu_z, d_lnrho, d_uu_x, d_uu_y, d_uu_z);
			//Sync
			cudaThreadSynchronize();
			//Save slice
			save_anim_slice(slice_axis, t_anim_number, slice_lnrho, slice_uu, slice_uu_x, slice_uu_y, slice_uu_z);
			//Update the next time when a frame is saved
			t_anim = t_anim + ANIM_STEP; 
			//Update the number for the next animation file 
			t_anim_number++;
			//save_anim_state(t_anim, t_anim_number); //TODO: Include something like this into the code. 
		}

		if (step % DIAG_STEPS == 0) {
			timeseries_diagnostics_cuda(d_umax, d_umin, d_urms, d_uxrms, d_uyrms, d_uzrms, d_rhorms, 
                                                    d_rhomax, d_uxmax, d_uymax, d_uzmax,
                                                    d_rhomin, d_uxmin, d_uymin, d_uzmin, 
                                                    step, dt, t, 
                                                    d_partial_result, d_lnrho, d_uu_x, d_uu_y, d_uu_z);
		}

		if (found_nan) {
			printf("Terminating computational loop due to nan.\n");
			step += MAX_STEPS; //End the loop if nan is found.
		}
		
		++step;
	}
	

/*	

		//--------------------------
		//Overview of the main loop
		//--------------------------
		//Concurrency depicted with 1,2,3 etc, 
		//=> same numbers execute concurrently:	

		//and perhaps this is more confusing than clarifying...
		
		//Point is, that kernels are queued on the GPU and CPU stops 
		//to only to execute its own instructions		

		//GPU BOTTLENECK VERSION:
		/*
		MV: The gpu bottleneck version seems somehow more intuitive for me. 
		MV: Do you have any idea which one would be theoretically more efficient?

		JP: This GPU bottleneck version should be more efficient with smaller grids
		as all GPU computation can focus solely on RK; with bigger grids CPU probably calculates the diagnostics
		slower than the GPU calculates the RK, so we would have to transfer some of the computation to
		the GPU => GPU will have to sacrifice RK time for calculating diagnostics. 

		MV: However the CPU bottleneck might not really be an cripling issue. This 
		MV: is because we will not need to calculate the diagnostic variables for 
		MV: every timestep so the possible slowdown might not bee too bad with 
		MV: large grids. The only problem here is actually the timestep. If we 
		MV: calculate the collective operation within the CPU, then we would still 
		MV: need to calculate the maximum velocity every timestep. This problem is 
		MV: of course solved by calculating the maximun velocity within the device, 
		MV: but the we need to still get into calculating collective operations 
		MV: with CUDA, which we wanted to avoid. Need to think about this... 
		
		JP: Well it's not really a problem to calculate the timestep because i
		already have an ultra optimized GPU code ready for it, and we need it
		for RK anyways, so the faster it's done, the better (GPU reductions are 
		insanely faster than anything on CPU anyways).

		* Device-Device cpy of comp. domain
		* 1. Async copy of of old comp. domain to host
		* 1. Device-Device cpy of comp. domain
		* 2. Call timestep kernel (GPU)
		* 3. Call RK3 kernel (GPU)
		* 2. Calculate slices & diagnostics from old data(CPU)
		* 3. Wait for RK3 to finish
		*/

		//CPU BOTTLENECK VERSION:
		/* 
		* 1. Call timestep kernel (GPU)
		* 2. Call RK3 kernel (GPU)
		* 3. Calculate slices & diagnostics (GPU)
		* 1. Wait for asyncmemcpy to finish and write old domain to file (CPU)
		* 4. Device-Device cpy of comp. domain)
		* 4,1-?. Async copy of of old comp. domain to host
		*/
	//}
	//----------------------------------------------------------

	printf("Avg time: %f\n", accumulated_time / steps);

	//----------------------------------------------------------
	// Load from device memory back into host for saving the final snapshot (TODO: Might need changes after writing async transfers )
	//----------------------------------------------------------
	checkErr( cudaMemcpy(lnrho, d_lnrho, sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	checkErr( cudaMemcpy(uu_x,  d_uu_x,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	checkErr( cudaMemcpy(uu_y,  d_uu_y,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	checkErr( cudaMemcpy(uu_z,  d_uu_z,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	//----------------------------------------------------------

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	FILE *fp;
	fp = fopen("finalresult_uu_x_zplane0v6.txt", "w");
		for(int i = 0; i < NX; i++){
			//fprintf(fp, "NX = %d , NY = %d \n\n\n ", NX, NY );
			for(int j = 0; j < NY; j++){
				fprintf(fp, "(%d, %d) %f \n",i, j, uu_x[(NX*i) + j]);
				//fprintf(fp, "%d | ", (NX*i) + j);
			}
				fprintf(fp, "------------------------------------------\n");
		}
	fclose(fp);

	/*checkErr( cudaMemcpy(uu_x,  d_w_uu_x,  sizeof(float)*GRID_SIZE, cudaMemcpyDeviceToHost) );
	FILE *fp2;
	fp2 = fopen("d_uu_x_dest_integrators_debug.txt", "w");
		for(int i = 0; i < NY; i++){
			//fprintf(fp, "NX = %d , NY = %d \n\n\n ", NX, NY );
			for(int j = 0; j < NX; j++){
				fprintf(fp2, "(%d, %d) %f, ",i, j, uu_x[(NX*i) + j]);
				//fprintf(fp2, "%d | ", (NX*i) + j);
			}
				fprintf(fp2, "\n\n------------------------------------------");
		}
	fclose(fp2);*/
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//----------------------------------------------------------
	// Save the final snapshot. 
	//----------------------------------------------------------
	save_grid_information(t); // Save grid information for the most current .dat file TODO: test. 
	save_grid_data(lnrho, DATA_LNRHO_PATH);
	save_grid_data(uu_x, DATA_UU_X_PATH);
	save_grid_data(uu_y, DATA_UU_Y_PATH);
	save_grid_data(uu_z, DATA_UU_Z_PATH);
	//----------------------------------------------------------

	//Destroy timers
	cudaEventDestroy( start );
	cudaEventDestroy( stop );

	//Free device memory
	checkErr( cudaFree(d_lnrho) );
	checkErr( cudaFree(d_uu_x) );
	checkErr( cudaFree(d_uu_y) );
	checkErr( cudaFree(d_uu_z) );

	//Free diagnostic helper variables/arrays
	checkErr( cudaFree(d_umax) ); checkErr( cudaFree(d_umin) );
	checkErr( cudaFree(d_urms) );
 	checkErr( cudaFree(d_uxrms) ); checkErr( cudaFree(d_uyrms) ); checkErr( cudaFree(d_uzrms) );
	checkErr( cudaFree(d_rhorms) );
	checkErr( cudaFree(d_rhomax) ); checkErr( cudaFree(d_rhomin) );
	checkErr( cudaFree(d_uxmax) ); checkErr( cudaFree(d_uymax) ); checkErr( cudaFree(d_uzmax) );
	checkErr( cudaFree(d_uxmin) ); checkErr( cudaFree(d_uymin) ); checkErr( cudaFree(d_uzmin) );
	checkErr( cudaFree(d_partial_result) );

	//Free pinned memory
	checkErr( cudaFreeHost(slice_lnrho) );
	checkErr( cudaFreeHost(slice_uu) );
	checkErr( cudaFreeHost(slice_uu_x) );
	checkErr( cudaFreeHost(slice_uu_y) );
	checkErr( cudaFreeHost(slice_uu_z) );

	//Reset device
	cudaDeviceReset();

	//Free host arrays
        free(lnrho);
        free(uu_x);
        free(uu_y);
        free(uu_z);

	if (LFORCING) {
		free(kk_x);
		free(kk_y);
		free(kk_z);
	}

	return EXIT_SUCCESS;
}	

