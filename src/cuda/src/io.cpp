
#include <stdio.h>
#include <stdlib.h>
//#include <string.h>

#include "io.h"
#include "defines.h"


//See if a file already exists. 
int exists(const char *fname)
{
   FILE *file;
   if (file = fopen(fname, "r"))
   {
      fclose(file);
      return 1;
   } else {
   return 0;
   }
}

//Loads data from given path to float* grid
//Returns 1 if successful, 0 otherwise
//!!Note!! animation_slice in compute.cu actually doesn't include the pads and bounds;
//todo to include pads&bounds
void save_grid_information(float time) 
{		

   //Save the basic information concernign the most current result grid that is needed e.g. for the purposes on data visualization 

   FILE* nnfile;

   // Format:
   //   NX   NY   NZ   GRID_SIZE   COMP_DOMAIN_SIZE_X   COMP_DOMAIN_SIZE_Y   COMP_DOMAIN_SIZE_Z   BOUND_SIZE   DX   DY   DZ   time
   
   nnfile = fopen("data/grid_info.ac", "w");
   fprintf(nnfile, "%i %i %i %i %i %i %i %i %i %f %f %f %f \n", NX, NY, NZ, GRID_SIZE, COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z, 
           BOUND_SIZE, DX, DY, DZ, time);
   fclose(nnfile);

}

int load_grid_data(float* grid, const char* path)
{
	FILE* binfile;

	//Read binfile
	binfile = fopen(path, "rb");
	if (binfile == NULL) { printf("Error opening %s\n", path); return 0; }
	fread(grid, sizeof(float), GRID_SIZE, binfile);
	fclose(binfile);
	//Done
	return 1;
}


//Save data from grid to path
int save_grid_data(float* grid, const char* path)
{
	FILE* binfile;

	//Write binfile
	binfile = fopen(path, "wb");
	if (binfile == NULL) { printf("Error opening %s\n", path); return 0; }	
	fwrite(grid, sizeof(float), GRID_SIZE, binfile);
	fclose(binfile);
	//Done
	return 1;
}


//Saves animation slice of the given x,y,z axis
//template <char slice_axis>
void save_anim_slice(char slice_axis, int slice_step, float* slice_lnrho, float* slice_uu, float* slice_uu_x, float* slice_uu_y, float* slice_uu_z)
{
   FILE* binfile_lnrho, *binfile_uu, *binfile_uu_x, *binfile_uu_y, *binfile_uu_z;

   const char* nn_path = "data/nn.ac";
   const char base1[] = "data/animation/lnrho";
   const char base2[] = "data/animation/uu";
   const char base3[] = "data/animation/uu_x";
   const char base4[] = "data/animation/uu_y";
   const char base5[] = "data/animation/uu_z";
   const char ending[] = ".ani";
   char filename1[64];
   char filename2[64];
   char filename3[64];
   char filename4[64];
   char filename5[64];

   //Define the filename 
   sprintf(filename1, "%s_%c_%d%s", base1, slice_axis, slice_step, ending);
   sprintf(filename2, "%s_%c_%d%s", base2, slice_axis, slice_step, ending);
   sprintf(filename3, "%s_%c_%d%s", base3, slice_axis, slice_step, ending);
   sprintf(filename4, "%s_%c_%d%s", base4, slice_axis, slice_step, ending);
   sprintf(filename5, "%s_%c_%d%s", base5, slice_axis, slice_step, ending);

   printf("Density slice: \"%s\"\n", filename1);
   printf("Velocity slice: \"%s\"\n", filename2);
   printf("Velocity_x slice: \"%s\"\n", filename3);
   printf("Velocity_y slice: \"%s\"\n", filename4);
   printf("Velocity_z slice: \"%s\"\n", filename5);

   binfile_lnrho = fopen(filename1, "wb");
   binfile_uu = fopen(filename2, "wb");
   binfile_uu_x = fopen(filename3, "wb");
   binfile_uu_y = fopen(filename4, "wb");
   binfile_uu_z = fopen(filename5, "wb");

   //-------------
   //Write grid dims (redundant, called multiple times during prog exec, can
   //be moved elsewhere if causes problems)
   FILE* nnfile = fopen(nn_path, "w");
   if (nnfile != NULL) 
	fprintf(nnfile, "%i %i %i \n", NX, NY, NZ);
   fclose(nnfile);
   //-------------

    if (binfile_lnrho 	== NULL 	|| 
	binfile_uu 	== NULL 	|| 
	binfile_uu_x 	== NULL 	|| 	
	binfile_uu_y 	== NULL 	|| 
	binfile_uu_z 	== NULL) 	{ printf("Error opening .ani slices!\n"); exit(-1); }
  
   //Write grid (Include pads&bounds)
   int slice_size;
   	switch(slice_axis) 
	{
		case 'x':
			slice_size = (COMP_DOMAIN_SIZE_Y + 2*BOUND_SIZE)*(COMP_DOMAIN_SIZE_Z + 2*BOUND_SIZE);
			break;
		case 'y':
			slice_size = (COMP_DOMAIN_SIZE_X + 2*BOUND_SIZE)*(COMP_DOMAIN_SIZE_Z + 2*BOUND_SIZE);
			break;
		case 'z':
			slice_size = (COMP_DOMAIN_SIZE_X + 2*BOUND_SIZE)*(COMP_DOMAIN_SIZE_Y + 2*BOUND_SIZE);
			break;
		default:
			printf("INVALID slice axis in save_anim_slice!\n");
			exit(-1);
	}
	//Write grid (DON'T include pads&bounds)
   /*
	int slice_size;
   	switch(slice_axis) 
	{
		case 'x':
			slice_size = COMP_DOMAIN_SIZE_Y*COMP_DOMAIN_SIZE_Z;
			break;
		case 'y':
			slice_size = COMP_DOMAIN_SIZE_X*COMP_DOMAIN_SIZE_Z;
			break;
		case 'z':
			slice_size = COMP_DOMAIN_SIZE_X*COMP_DOMAIN_SIZE_Y;
			break;
		default:
			printf("INVALID slice axis in save_anim_slice!\n");
			exit(-1);
	}*/
	

   fwrite(slice_lnrho, sizeof(float), slice_size, binfile_lnrho);
   fwrite(slice_uu, sizeof(float), slice_size, binfile_uu);
   fwrite(slice_uu_x, sizeof(float), slice_size, binfile_uu_x);
   fwrite(slice_uu_y, sizeof(float), slice_size, binfile_uu_y);
   fwrite(slice_uu_z, sizeof(float), slice_size, binfile_uu_z);
   //Close 

   fclose(binfile_lnrho);
   fclose(binfile_uu);
   fclose(binfile_uu_x);
   fclose(binfile_uu_y);
   fclose(binfile_uu_z);


}

void init_ts_file()
{ 
	// Intialize the first line of the diagnostic file. 
	FILE* diagfile;
	int existence;

	existence = exists("data/ts.ac");
	printf("Existence %i  \n", existence);

	if (existence == 0) {
		diagfile = fopen("data/ts.ac", "w");
		fprintf(diagfile, "step t dt urms uxrms uyrms uzrms uxmax uymax uzmax uxmin uymin uzmin rhorms umax rhomax rhomin umin \n");
		fclose(diagfile);
	}

}

void save_ts(float t, float dt, int step, float urms, float uxrms, float uyrms, float uzrms, float uxmax, float uymax, float uzmax, float rhorms, float umax, float rhomax,
             float uxmin, float uymin, float uzmin, float rhomin, float umin)
{ 
	// Saves diagnostic variables a file during computation
	FILE* diagfile;
	diagfile = fopen("data/ts.ac", "a");

	fprintf(diagfile, "%i %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n", 
                step, t, dt, urms, uxrms, uyrms, uzrms, uxmax, uymax, uzmax, uxmin, uymin, uzmin, rhorms, umax, rhomax, rhomin, umin);

	fclose(diagfile);
}

void init_used_k_file(float kaver)
{
        // Intialize the first line of the k-vector diagnostic file. 
        FILE* used_k_file;
        int existence;

        existence = exists("data/kk_used.ac");

        if (existence == 0) {
                used_k_file = fopen("data/kk_used.ac", "w");
		fprintf(used_k_file, "kaver = %f \n", kaver);
                fprintf(used_k_file, "kk_vec_x kk_vec_y kk_vec_z phi forcing_kk_part_x forcing_kk_part_y forcing_kk_part_z \n");
                fclose(used_k_file);
        }

}


void save_used_k(float kk_vec_x, float kk_vec_y, float kk_vec_z, float phi, float forcing_kk_part_x, float forcing_kk_part_y, float forcing_kk_part_z)
{
	// Saves values of used k-values and phase angle per step for diagnostic purposes
	FILE* used_k_file;
	used_k_file = fopen("data/kk_used.ac", "a");
	fprintf(used_k_file, "%e %e %e %e %e %e %e \n", kk_vec_x, kk_vec_y, kk_vec_z, phi, forcing_kk_part_x, forcing_kk_part_y, forcing_kk_part_z);
	fclose(used_k_file);
}

