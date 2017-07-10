#include <cmath>
#include <cstdio>

#include "forcing.cuh"
#define EXTERN extern
#include "dconsts.cuh"

__device__ void forcing_simple(float* fx, float* fy, float* fz, float grid_idx_x, float grid_idx_y, float grid_idx_z)
{
	// Calculate simple non-helical forcing for a grid point
	float waves; 
	float k_dot_x;

	float coord_x, coord_y, coord_z;

	// Coordinate vector
	coord_x = d_DX*grid_idx_x - d_XORIG;
	coord_y = d_DY*grid_idx_y - d_YORIG; 
	coord_z = d_DZ*grid_idx_z - d_ZORIG; 

   	// Wave number vector k 

   	k_dot_x = coord_x*d_KK_VEC_X + coord_y*d_KK_VEC_Y + coord_z*d_KK_VEC_Z;

	waves = cos(k_dot_x)*cos(d_PHI) - sin(k_dot_x)*sin(d_PHI);

	*fx = d_FORCING_KK_PART_X*waves;
	*fy = d_FORCING_KK_PART_Y*waves;
	*fz = d_FORCING_KK_PART_Z*waves;

	//if (*fx != *fx) {
	//&    printf("GOT NAN IN FORCING FUNCTION!!! \n waves = %e, phi = %e, k_dot_x = %e, d_FORCING_KK_PART_X = %e \n", waves, phi, k_dot_x, d_FORCING_KK_PART_X);
	//    exit(1); 
	//}

	//printf(" k_dot_x fk_x fk_y fk_z : %e  %e  %e  %e  %e  %e  %e  %e  \n",  k_dot_x, fk_x, fk_y, fk_z, N, dt, cs, f0);
	//printf("realfx realfy realfz waves FX FY FZ: %e  %e  %e  %e  %e  %e  %e \n", realfx, realfy, realfz, waves, *fx, *fy, *fz);
}

__device__ void forcing_nonhelical(float s_uu_x[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                               	   float s_uu_y[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
                              	   float s_uu_z[SHARED_SIZE_ROW][SHARED_SIZE_COL][SHARED_SIZE_DEPTH], 
				   int sid_row, int sid_column, int sid_depth,    
				   float grid_idx_x, float grid_idx_y, float grid_idx_z) 
{
	// Calculate simple non-helical forcing for a grid point
	float waves; 
	float k_dot_x;

	float coord_x, coord_y, coord_z;

	// Coordinate vector
	coord_x = d_DX*grid_idx_x - d_XORIG;
	coord_y = d_DY*grid_idx_y - d_YORIG; 
	coord_z = d_DZ*grid_idx_z - d_ZORIG; 

   	// Wave number vector k 

   	k_dot_x = coord_x*d_KK_VEC_X + coord_y*d_KK_VEC_Y + coord_z*d_KK_VEC_Z;

	waves = cos(k_dot_x)*cos(d_PHI) - sin(k_dot_x)*sin(d_PHI);

	s_uu_x[sid_row][sid_column][sid_depth] = s_uu_x[sid_row][sid_column][sid_depth] + d_FORCING_KK_PART_X*waves;
	s_uu_y[sid_row][sid_column][sid_depth] = s_uu_y[sid_row][sid_column][sid_depth] + d_FORCING_KK_PART_Y*waves;
	s_uu_z[sid_row][sid_column][sid_depth] = s_uu_z[sid_row][sid_column][sid_depth] + d_FORCING_KK_PART_Z*waves;
}

int calc_k_vector_size()
{
	// Calculate the size of the k-vector
	int nk;
	float kxstep, kystep, kzstep;
	float k, kref;
        float kmax_x, kmax_y, kmax_z;
        float dkstep_x, dkstep_y, dkstep_z;

        // Check if we leave out some chosen axis for simplicity (ie. testing)
        // Still a bit redneck solution for this. :P 
        if (DKX != 0.0) {kmax_x = KMAX; dkstep_x = DKX;} else {kmax_x = 0.0; dkstep_x = 1.0;}
        if (DKY != 0.0) {kmax_y = KMAX; dkstep_y = DKY;} else {kmax_y = 0.0; dkstep_y = 1.0;}
        if (DKZ != 0.0) {kmax_z = KMAX; dkstep_z = DKZ;} else {kmax_z = 0.0; dkstep_z = 1.0;}

	nk = 0;


        kxstep = -kmax_x; kystep = -kmax_y; kzstep = -kmax_z;
        while (kxstep <= kmax_x) {
                while (kystep <= kmax_y) {
                        while (kzstep <= kmax_z) {
				k=sqrt(kxstep*kxstep+kystep*kystep+kzstep*kzstep);
				kref = k;
				if (kref > KK1 && kref < KK2) {
					nk++;
				}
                                kzstep = kzstep + dkstep_z;
                        }
                        kzstep = -kmax_z;
                        kystep = kystep + dkstep_y;
                }
                kystep = -kmax_y;
                kxstep = kxstep + dkstep_x;
        }

	return nk;

}

void initialize_k_vectors(float* kk_x, float* kk_y, float* kk_z, float* kaver, int nk)
{
	// This function generates a chosen group of k vectors for the forcing function 
	int i;
	float kxstep, kystep, kzstep;
	float k, k_sum, kref;
	float kmax_x, kmax_y, kmax_z;
	float dkstep_x, dkstep_y, dkstep_z;

	// Check if we leave out some chosen axis for simplicity (ie. testing)
	// Still a bit redneck solution for this. :P 
	if (DKX != 0.0) {kmax_x = KMAX; dkstep_x = DKX;} else {kmax_x = 0.0; dkstep_x = 1.0;}
        if (DKY != 0.0) {kmax_y = KMAX; dkstep_y = DKY;} else {kmax_y = 0.0; dkstep_y = 1.0;}
        if (DKZ != 0.0) {kmax_z = KMAX; dkstep_z = DKZ;} else {kmax_z = 0.0; dkstep_z = 1.0;}

	k_sum = 0.0;
	i = 0 ;

	kxstep = -kmax_x; kystep = -kmax_y; kzstep = -kmax_z;
	while (kxstep <= kmax_x) {
		while (kystep <= kmax_y) {
			while (kzstep <= kmax_z) {
				k=sqrt(kxstep*kxstep+kystep*kystep+kzstep*kzstep);
				kref = k;
				//   kref=sqrt((kxstep/ex)^2+(kystep/ey)^2+(kzstep/ez)^2) what are ex, ey, ez ?? TODO: check this 
				if (kref > KK1 && kref < KK2) {
					k_sum=k_sum+kref;
					kk_x[i] = kxstep;
					kk_y[i] = kystep;
					kk_z[i] = kzstep;
					//  printf("%f %f %f ||| %f/%f/%f %f %f %f %i \n", kx[i], ky[i], kz[i], KK1, kref, KK2, kxstep, kystep, kzstep, nk);
					i++;
				}
				kzstep = kzstep + dkstep_z;
			}
			kzstep = -kmax_z;
			kystep = kystep + dkstep_y;
		} 
		kystep = -kmax_y; 
		kxstep = kxstep + dkstep_x;  
	}
  
	float n;

	n      = (float) nk;
	*kaver = k_sum/n;
	printf("Number of k-vectors nk = %i \n The average wavenumber kaver = %f with KK1 = %f KK2 = %f \n", nk, *kaver, KK1, KK2);

	printf("Listing all possible k-vectors:\n");
        printf("kk_x kk_y kk_z\n");
	for (i = 0; i < nk; i++) {
		printf("%.2f %.2f %.2f \n", kk_x[i], kk_y[i], kk_z[i]);
	}
}

void choose_random_k(float* kk_vec_x, float* kk_vec_y, float* kk_vec_z, float* kk_x, float* kk_y, float* kk_z, int nk)
{
	// Chooses a random k vector for forcing during one timestep
	int rind, ind;
	float find, fnk;
	float max = (float) RAND_MAX;

	// Randomized index
	//std::srand(std::time(0));
	rind = std::rand();
	find = (float) rind;
	fnk  = (float) nk;
	find = (find/max)*fnk;
	ind  = (int) find;

	// Select the k-vector for the timestep
	*kk_vec_x = kk_x[ind];
	*kk_vec_y = kk_y[ind];
	*kk_vec_z = kk_z[ind];

	//  printf("k_ind-arvot: %i, %f, %i, %f \n", ind, find, rind, fnk);
	//  printf("k11-arvot: %f, %f, %f \n", *kk_vec_x, *kk_vec_y, *kk_vec_z);

}

float choose_random_phi()
{
	//Fetch randomized phase angle phi
	int rphi;
	float phi;
	float max = (float) RAND_MAX;

	//Get and convert random numbers
	rphi = std::rand();
	phi = (float) rphi;

	//Get the angle
	//GENERATES WRONG DISTRIBUTION: from -pi/2 to pi/2
	//phi = (phi-(max/2.0))/max; // from -1 to 1 
	//NEW METHOD DOES THIS:
	//phi = phi/max; //from 0 to 1
	//phi = phi-0.5; //from -0.5 to 0.5 
	//phi = phi*2.0; //from -1 to 1
	phi = 2.0*(phi/max-0.5); //from -1 to 1
	phi = M_PI*phi;          // from -pi to pi
	//TODO: Test. This might be the source for centre acceleration.

	return phi;
}

void cross(float* rvecx, float* rvecy, float* rvecz, float avecx, float avecy, float avecz, float bvecx, float bvecy, float bvecz) {
	//Cross product, r = a x b

	*rvecx = avecy*bvecz - avecz*bvecy;
	*rvecy = avecz*bvecx - avecx*bvecz;
	*rvecz = avecx*bvecy - avecy*bvecx;

}

void dot(float* rscal, float avecx, float avecy, float avecz, float bvecx, float bvecy, float bvecz) {
	//dot product, r = a . b

	*rscal = avecx*bvecx + avecy*bvecy + avecz*bvecz;

}

void get_forcing_unit_vector(float* ex, float* ey, float* ez, float kk_vec_x, float kk_vec_y, float kk_vec_z)
{
	//TODO: Do it like the Pencil code   

	float ex0, ey0, ez0;
	float ex1, ey1, ez1;
	float ex2, ey2, ez2;
	float norm, sqrt_norm, phi;

	if ((kk_vec_y == 0) && (kk_vec_z == 0)) {
		ex0 = 0.0; 
		ey0 = 1.0; 
		ez0 = 0.0;
        } else {
		ex0 = 1.0; 
		ey0 = 0.0; 
		ez0 = 0.0;
	}

	//Above unite vector component are relaxed as follows 
	//(1) constructing two basis vectors for the plane perpendicular to the wave vector kk_vec
	//(2) choosing a random direction in that plane (angle phi)

	//Unit vector perp. to kk
	cross(&ex1, &ey1, &ez1, kk_vec_x, kk_vec_y, kk_vec_z, ex0, ey0, ez0);
	dot(&norm, ex1, ey1, ez1, ex1, ey1, ez1);
	sqrt_norm = sqrt(norm); 
	ex1 = ex1/sqrt_norm;  ey1 = ey1/sqrt_norm;  ez1 = ez1/sqrt_norm;

	//Unit vector perp. to kk, e1
	cross(&ex2, &ey2, &ez2, kk_vec_x, kk_vec_y, kk_vec_z, ex1, ey1, ez1);
	dot(&norm, ex2, ey2, ez2, ex2, ey2, ez2);
	sqrt_norm = sqrt(norm); 
	ex2 = ex2/sqrt_norm;  ey2 = ey2/sqrt_norm;  ez2 = ez2/sqrt_norm;

	//Randomized direction
	phi = choose_random_phi();
	*ex = cos(phi)*ex1 + sin(phi)*ex2;
	*ey = cos(phi)*ey1 + sin(phi)*ey2;
	*ez = cos(phi)*ez1 + sin(phi)*ez2;


}

void get_forcing_unit_vector_old(float* ex, float* ey, float* ez)
{
	//The old (nonfunctional) version for comparison 
	//For forcing
	//Calculate value for the unit vector e, which fulfilles the criteria
	// (k \times e) \cdot k = 0  and |e| = 1
	// We get the arbitary general form:
	// e = (A, \sqrt{1-2A^2}, A)
	// And by assigning an arbitrary numerical value A = 0.5 we get
	// e = (0.5, \sqrt{0.5}, 0.5). So: 

	*ex = 0.5;
	
	*ey = sqrt(0.5);

	*ez = 0.5; 
}

void fkt_forcing_coefficient(float* forcing_kk_part_x, float* forcing_kk_part_y, float* forcing_kk_part_z, float kk_vec_x, float kk_vec_y, float kk_vec_z, float dt)
{ 
	//Pre-calcute the parts where the forcing function which do not vary for different coordinates. 
	float N_coeff, fk_under, kk; 
	float ex, ey, ez, k_cross_ex, k_cross_ey, k_cross_ez;
	float k_dot_e;

	//The magnitude of k-vector
	kk = sqrt(kk_vec_x*kk_vec_x + kk_vec_y*kk_vec_y + kk_vec_z*kk_vec_z);

	// Unit vector perpendicular to k x e
	get_forcing_unit_vector(&ex, &ey, &ez, kk_vec_x, kk_vec_y, kk_vec_z);

	printf("ex = %e, ey = %e, ez = %e \n", ex, ey, ez);
	printf("kk_vec_x = %e, kk_vec_y = %e, kk_vec_z = %e \n", kk_vec_x, kk_vec_y, kk_vec_z);

	dot(&k_dot_e, kk_vec_x, kk_vec_y, kk_vec_z, ex, ey, ez);

	cross(&k_cross_ex, &k_cross_ey, &k_cross_ez, kk_vec_x, kk_vec_y, kk_vec_z, ex, ey, ez);

	//Magnitude coefficient
	N_coeff = FORCING*CS_SOUND*pow(kk*CS_SOUND/dt, 0.5);

	//The underside of nonhelical f_k (forcing coefficient)
	fk_under = sqrt(kk*kk - k_dot_e*k_dot_e);

	//The components of nonhelical f_k (forcing coefficient) 
	if (fk_under > 0.0) {
		*forcing_kk_part_x = N_coeff*(k_cross_ex / fk_under);
		*forcing_kk_part_y = N_coeff*(k_cross_ey / fk_under);
		*forcing_kk_part_z = N_coeff*(k_cross_ez / fk_under);
	} else {
		printf("Forcing skipped this step because of div by zero!");
		*forcing_kk_part_x = 0.0;
		*forcing_kk_part_y = 0.0;
		*forcing_kk_part_z = 0.0;
	}

	if ((*forcing_kk_part_x) != (*forcing_kk_part_x)) {
		printf("GOT NAN IN FORCING FUNCTION!!! \n forcing_kk_part_x = %e, k_cross_ex %e, fk_under %e, kk %e, k_dot_e %e \n", 
                       *forcing_kk_part_x, k_cross_ex, fk_under, kk, k_dot_e);
		exit(1); 
	}

}

