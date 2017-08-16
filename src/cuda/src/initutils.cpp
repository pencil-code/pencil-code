
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defines.h"

/*
* Contains the functions for initializing the grid with CPU
*/

//Initialize all grid values into 0.0 to avoid weird errors
void set_grids_zero(float* lnrho, float* uu_x, float* uu_y, float* uu_z)
{
	int idx;

	printf("Formating all grid values to 0.0 ...\n");
	
	//Note: initialises also the pad zone and boundaries
	for (int k=0; k < NZ; k++) {
		for (int j=0; j < NY; j++) {
			for (int i=0; i < NX; i++) {
				idx = i + j*NX + k*NX*NY;
				lnrho[idx] = 0.0;
				uu_x[idx] = 0.0;
				uu_y[idx] = 0.0;
				uu_z[idx] = 0.0;
			}
		}
	}
}

//Set all lnrho's items to to AMPL_LNRHO
void lnrho_const(float* lnrho)	
{
	int idx;
	
	//Note: initialises also the pad zone and boundaries
	for (int k=0; k < NZ; k++) {
		for (int j=0; j < NY; j++) {
			for (int i=0; i < NX; i++) {
				idx = i + j*NX + k*NX*NY;
				lnrho[idx] = AMPL_LNRHO;
			}
		}
	}
}

//Set all velocities of the chosen component to  AMPL_UU
void uu_const(float* uu)
{
	int idx;

	//Note: initialises also the pad zone and boundaries
	for (int k=0; k < NZ; k++) {
		for (int j=0; j < NY; j++) {
			for (int i=0; i < NX; i++) {
				idx = i + j*NX + k*NX*NY;
				uu[idx] = AMPL_UU;
			}
		}
	}
}


//Initialize velocity grid with uniformly random velocities in range (-AMPL_UU, AMPL_UU)
void uniform_random_velocity(float* uu_x, float* uu_y, float* uu_z)
{
   	int idx;
	float rnd; //random float between (-1.0, 1.0)

	srand(1000); //Init random with seed

	for (int k=CZ_BOT; k < CZ_TOP; k++) {
		for (int j=CY_BOT; j < CY_TOP; j++) {
			for (int i=CX_BOT; i < CX_TOP; i++) {
				idx = i + j*NX + k*NX*NY;
			
				//Shuffle random float in range (-1, 1)
				rnd = 2.0*((float)rand() / (float)RAND_MAX) - 1.0; 
				uu_x[idx] = AMPL_UU * rnd;

				rnd = 2.0*((float)rand() / (float)RAND_MAX) - 1.0; 
				uu_y[idx] = AMPL_UU * rnd;

				rnd = 2.0*((float)rand() / (float)RAND_MAX) - 1.0; 
				uu_z[idx] = AMPL_UU * rnd;
			}
		}
	}
}

void xgaussian_wave(float* uu)	
{
	// Gaussian velocity wave in yz-plane
   	int idx;
	float xx;

	for (int k=CZ_BOT; k < CZ_TOP; k++) {
		for (int j=CY_BOT; j < CY_TOP; j++) {
			for (int i=CX_BOT; i < CX_TOP; i++) {
				idx = i + j*NX + k*NX*NY;
				xx = DX*(i-CX_BOT)-XORIG;
				uu[idx] = AMPL_UU*exp( -pow((xx-INIT_LOC_UU_X),2.0) / (2.0*pow(WIDTH_UU, 2.0)) );
			}
		}
	}
}

void ygaussian_wave(float* uu)	
{
	// Gaussian velocity wave in xz-plane
   	int idx;
	float yy;

	for (int k=CZ_BOT; k < CZ_TOP; k++) {
		for (int j=CY_BOT; j < CY_TOP; j++) {
			for (int i=CX_BOT; i < CX_TOP; i++) {
				idx = i + j*NX + k*NX*NY;
				yy = DY*(j-CY_BOT)-YORIG;
				uu[idx] = AMPL_UU*exp( -pow((yy-INIT_LOC_UU_Y),2.0) / (2.0*pow(WIDTH_UU, 2.0)) );
			}
		}
	}
}

void zgaussian_wave(float* uu)	
{
	// Gaussian velocity wave in xy-plane
   	int idx;
	float zz; 

	for (int k=CZ_BOT; k < CZ_TOP; k++) {
		for (int j=CY_BOT; j < CY_TOP; j++) {
			for (int i=CX_BOT; i < CX_TOP; i++) {
				idx = i + j*NX + k*NX*NY;
				zz = DZ*(k-CZ_BOT)-ZORIG;
				uu[idx] = AMPL_UU*exp( -pow((zz-INIT_LOC_UU_Z),2.0) / (2.0*pow(WIDTH_UU, 2.0)) );
			}
		}
	}
}

void sin_kx_wave(float* uu)
{
        // Initial condition for the dissipation test
        // u_y = sin(k*x)
        int idx;
        float xx;

        for (int k=CZ_BOT; k < CZ_TOP; k++) {
                for (int j=CY_BOT; j < CY_TOP; j++) {
                        for (int i=CX_BOT; i < CX_TOP; i++) {
                                //TODO: This might cause a boundary effect if things do not match. Investigate!
                                idx = i + j*NX + k*NX*NY;
                                //xx = DX*(i-CX_BOT)-XORIG;
                                xx = DX*((float) i - (float) CX_BOT)- XORIG; 
                                uu[idx] = AMPL_UU*sin(INIT_K_WAVE*(xx-INIT_LOC_UU_X));
                        }
                }
        }
}

void sin_ky_wave(float* uu)
{
        // Initial condition for the dissipation test
        // u_z = sin(k*y)
        int idx;
        float yy;

        for (int k=CZ_BOT; k < CZ_TOP; k++) {
                for (int j=CY_BOT; j < CY_TOP; j++) {
                        for (int i=CX_BOT; i < CX_TOP; i++) {
                                //TODO: This might cause a boundary effect if things do not match. Investigate!
                                idx = i + j*NX + k*NX*NY;
                                yy = DY*((float) j - (float) CY_BOT)- YORIG; 
                                uu[idx] = AMPL_UU*sin(INIT_K_WAVE*(yy-INIT_LOC_UU_Y));
                        }
                }
        }
}

void sin_kz_wave(float* uu)
{
        // Initial condition for the dissipation test
        // u_x = sin(k*z)
        int idx;
        float zz;

        for (int k=CZ_BOT; k < CZ_TOP; k++) {
                for (int j=CY_BOT; j < CY_TOP; j++) {
                        for (int i=CX_BOT; i < CX_TOP; i++) {
                                //TODO: This might cause a boundary effect if things do not match. Investigate!
                                idx = i + j*NX + k*NX*NY;
                                zz = DZ*((float) k - (float) CZ_BOT)- ZORIG; 
                                uu[idx] = AMPL_UU*sin(INIT_K_WAVE*(zz-INIT_LOC_UU_Z));
                        }
                }
        }
}

void periodic_step_function_x(float* uu)
{
	// A periodix step function to test the effect of resolution into the system. 
        // u_x = tanh(cos x / delta x)
        
        int idx;
        float xx;

        for (int k=CZ_BOT; k < CZ_TOP; k++) {
                for (int j=CY_BOT; j < CY_TOP; j++) {
                        for (int i=CX_BOT; i < CX_TOP; i++) {
                                idx = i + j*NX + k*NX*NY;
                                xx = DX*((float) i - (float) CX_BOT)- XORIG; 
                                uu[idx] = AMPL_UU*tanh(cos(xx)/WIDTH_UU);
                        }
                }
        }

}

void periodic_step_function_y(float* uu)
{
	// A periodix step function to test the effect of resolution into the system. 
        // u_y = tanh(cos y / delta y)
        
        int idx;
        float yy;

        for (int k=CZ_BOT; k < CZ_TOP; k++) {
                for (int j=CY_BOT; j < CY_TOP; j++) {
                        for (int i=CX_BOT; i < CX_TOP; i++) {
                                idx = i + j*NX + k*NX*NY;
                                yy = DY*((float) j - (float) CY_BOT)- YORIG; 
                                uu[idx] = AMPL_UU*tanh(cos(yy)/WIDTH_UU);
                        }
                }
        }

}


void periodic_step_function_z(float* uu)
{
	// A periodix step function to test the effect of resolution into the system. 
        // u_z = tanh(cos z / delta z)
        
        int idx;
        float zz;

        for (int k=CZ_BOT; k < CZ_TOP; k++) {
                for (int j=CY_BOT; j < CY_TOP; j++) {
                        for (int i=CX_BOT; i < CX_TOP; i++) {
                                idx = i + j*NX + k*NX*NY;
                                zz = DZ*((float) k - (float) CZ_BOT)- ZORIG; 
                                uu[idx] = AMPL_UU*tanh(cos(zz)/WIDTH_UU);
                        }
                }
        }

}

void gaussian_3d_ball(float* uu)	
{
	// Gaussian velocity wave in 3d-ball 
   	int idx;
	float xx, yy, zz, rr2; 

	for (int k=CZ_BOT; k < CZ_TOP; k++) {
		for (int j=CY_BOT; j < CY_TOP; j++) {
			for (int i=CX_BOT; i < CX_TOP; i++) {
				idx = i + j*NX + k*NX*NY;
                                xx = DX*(i-CX_BOT)-XORIG;
                                yy = DY*(j-CY_BOT)-YORIG;
                                zz = DZ*(k-CZ_BOT)-ZORIG;
				rr2 = pow((xx-INIT_LOC_UU_X),2.0) + pow((yy-INIT_LOC_UU_Y),2.0) + pow((zz-INIT_LOC_UU_Z),2.0);
				uu[idx] = AMPL_UU*exp( -rr2 / (2.0*pow(WIDTH_UU, 2.0)) );
			}
		}
	}
}

void gaussian_radial_explosion(float* uu_x, float* uu_y, float* uu_z)	
{
	// Outward explosion with gaussian initial velocity profile. 
   	int idx;
	float xx, yy, zz, rr2, rr, theta, phi; 
	float uu_radial, theta_old;

	theta_old = 0;

	for (int k=CZ_BOT; k < CZ_TOP; k++) {
		for (int j=CY_BOT; j < CY_TOP; j++) {
			for (int i=CX_BOT; i < CX_TOP; i++) {
				//Calculate the value of velocity in a particular radius.
				idx = i + j*NX + k*NX*NY;
				//Determine the coordinates
                                xx = DX*(i-CX_BOT)-XORIG;
				xx = xx-INIT_LOC_UU_X;

                                yy = DY*(j-CY_BOT)-YORIG;
                                yy = yy-INIT_LOC_UU_Y;

                                zz = DZ*(k-CZ_BOT)-ZORIG;
                                zz = zz-INIT_LOC_UU_Z;

				rr2 = pow(xx,2.0) + pow(yy,2.0) + pow(zz,2.0);
				rr = sqrt(rr2);

				// Origin is different!
				float xx_abs, yy_abs, zz_abs;
				if (rr > 0.0) { 
					// theta range [0, PI]
					if (zz >= 0.0) {
						theta = acos(zz/rr);
                                		if (theta > M_PI/2.0 || theta < 0.0) {
                                        		printf("Explosion THETA WRONG: zz = %.3f, rr = %.3f, theta = %.3e/PI, M_PI = %.3e\n",
                                        		zz, rr, theta/M_PI, M_PI);
						}
					} else {
						zz_abs = -zz; // Needs a posite value for acos
						theta = M_PI - acos(zz_abs/rr); 
                            			if (theta < M_PI/2.0 || theta > 2*M_PI) {
  		                         		printf("Explosion THETA WRONG: zz = %.3f, rr = %.3f, theta = %.3e/PI, M_PI = %.3e\n", zz, rr, theta/M_PI, M_PI);
						}
					}

					// phi range [0, 2*PI]i
					if (xx != 0.0) {
						if (xx < 0.0 && yy >= 0.0) {
							//-+
							xx_abs = -xx; // Needs a posite value for atan
							phi = M_PI - atan(yy/xx_abs);
							if (phi < (M_PI/2.0) || phi > M_PI) {
								printf("Explosion PHI WRONG -+: xx = %.3f, yy = %.3f, phi = %.3e/PI, M_PI = %.3e\n", xx, yy, phi/M_PI, M_PI);
							}
						} else if (xx > 0.0 && yy < 0.0) {
							//+-
							yy_abs = -yy;
							phi = 2.0*M_PI - atan(yy_abs/xx);
							if (phi < (3.0*M_PI)/2.0 || phi > (2.0*M_PI + 1e-6)) {
								printf("Explosion PHI WRONG +-: xx = %.3f, yy = %.3f, phi = %.3e/PI, M_PI = %.3e\n", xx, yy, phi/M_PI, M_PI);
                                                        }
						} else if (xx < 0.0 && yy < 0.0) {
							//--
							yy_abs = -yy;
							xx_abs = -xx;
							phi = M_PI + atan(yy_abs/xx_abs);
                                                        if (phi < M_PI || phi > ((3.0*M_PI)/2.0 + 1e-6)) {
                                                                printf("Explosion PHI WRONG --: xx = %.3f, yy = %.3f, xx_abs = %.3f, yy_abs = %.3f, phi = %.3e, (3.0*M_PI)/2.0 = %.3e\n", 
								xx, yy, xx_abs, yy_abs, phi, (3.0*M_PI)/2.0);
                                                        }
						} else {
							//++
							phi = atan(yy/xx);
                                                        if (phi < 0 || phi > M_PI/2.0) {
                                                                printf("Explosion PHI WRONG --: xx = %.3f, yy = %.3f, phi = %.3e,\n", 
                                                                xx, yy, phi);
                                                        }
						}
					} else { //To avoid div by zero with atan
						if (yy > 0.0) {
							phi = M_PI / 2.0;
						} else if (yy < 0.0 ) {
							phi = (3.0*M_PI) / 2.0 ;
						} else {
							phi = 0.0;
						}
					}

					//Set zero for explicit safekeeping
					if (xx == 0.0 && yy == 0.0) {
						phi = 0.0; 
					}

					//Gaussian velocity 
					//uu_radial = AMPL_UU*exp( -rr2 / (2.0*pow(WIDTH_UU, 2.0)) );
					//New distribution, where that gaussion wave is not in the exact centre coordinates
					//uu_radial = AMPL_UU*exp( -pow((rr - 4.0*WIDTH_UU),2.0) / (2.0*pow(WIDTH_UU, 2.0)) ); //TODO: Parametrize the peak location.
                                        uu_radial = AMPL_UU*exp( -pow((rr - UU_SHELL_R),2.0) / (2.0*pow(WIDTH_UU, 2.0)) ); 

				} else { 
					uu_radial = 0.0; //TODO: There will be a discontinuity in the origin... Should the shape of the distribution be different?
				}

				//Determine the carthesian velocity components

				uu_x[idx] = uu_radial*sin(theta)*cos(phi);
				uu_y[idx] = uu_radial*sin(theta)*sin(phi);
				uu_z[idx] = uu_radial*cos(theta);

				//Temporary diagnosticv output (TODO: Remove after not needed)
				//if (theta > theta_old) {
				//if (theta > M_PI || theta < 0.0 || phi < 0.0 || phi > 2*M_PI) {
				/*	printf("Explosion: xx = %.3f, yy = %.3f, zz = %.3f, rr = %.3f, phi = %.3e/PI, theta = %.3e/PI\n, M_PI = %.3e",
						xx, yy, zz, rr, phi/M_PI, theta/M_PI, M_PI);
					printf("           uu_radial = %.3e, uu_x[%i] = %.3e, uu_y[%i] = %.3e, uu_z[%i] = %.3e \n", 
						uu_radial, idx, uu_x[idx], idx, uu_y[idx], idx, uu_z[idx]);
					theta_old = theta;
				*/

			}
		}
	}
}



void debug_clear_grid(float* lnrho, float* uu_x, float* uu_y, float* uu_z)
{
	for (int k=0; k < NZ; k++) {
		for (int j=0; j < NY; j++) {
			for (int i=0; i < NX; i++) {
				int idx = i + j*NX + k*NX*NY;
				lnrho[idx] = 0.0;
				uu_x[idx] = 0.0;
				uu_y[idx] = 0.0;
				uu_z[idx] = 0.0;
			}
		}
	}
}

#define GRID_NORMAL 0
#define GRID_RANDOM 1
#define GRID_COMP_DOMAIN_ONLY 2
#define GRID_BOUND_ZONES_ONLY 3
#define GRID_NO_FLOAT_ARITHMETIC_ERROR 4
#define GRID_WAVE 5
#define GRID_SMEM_TEST 6
#define GRID_BOUND_CPY_ZONE_ONLY 7
#define GRID_ALL_ZERO 8
void debug_grid_init(float* lnrho, float* uu_x, float* uu_y, float* uu_z, int type) 
{
	int idx;
	float rnd;
	int magic_number = 1000;
	srand(magic_number);
	debug_clear_grid(lnrho, uu_x, uu_y, uu_z);

	switch (type){
	case GRID_NORMAL:
		for (int k=0; k < NZ; k++) {
			for (int j=0; j < NY; j++) {
				for (int i=0; i < NX; i++) {
					int idx = i + j*NX + k*NX*NY;
					lnrho[idx] = i+j-k;
					uu_x[idx] = i;
					uu_y[idx] = j;
					uu_z[idx] = k;
				}
			}
		}
		break;

	case GRID_RANDOM:
		for (int k=0; k < NZ; k++) {
			for (int j=0; j < NY; j++) {
				for (int i=0; i < NX; i++) {
					int idx = i + j*NX + k*NX*NY;
					lnrho[idx] = rand() % magic_number;
					uu_x[idx] = rand() % magic_number;
					uu_y[idx] = rand() % magic_number;
					uu_z[idx] = rand() % magic_number;
				}
			}
		}
		break;

	case GRID_COMP_DOMAIN_ONLY:
		for (int k=0; k < NZ; k++) {
			for (int j=0; j < NY; j++) {
				for (int i=0; i < NX; i++) {
					int idx = i + j*NX + k*NX*NY;
					lnrho[idx] = 0;
					uu_x[idx] = 0;
					uu_y[idx] = 0;
					uu_z[idx] = 0;
				}
			}
		}
		for (int k=CZ_BOT; k < CZ_TOP; k++) {
			for (int j=CY_BOT; j < CY_TOP; j++) {
				for (int i=CX_BOT; i < CX_TOP; i++) {
					int idx = i + j*NX + k*NX*NY;
					lnrho[idx] = 5;//rand() % magic_number;
					uu_x[idx] = 6; //rand() % magic_number;
					uu_y[idx] = 7; //rand() % magic_number;
					uu_z[idx] = 8; //rand() % magic_number;
				}
			}
		}
		break;

	case GRID_BOUND_ZONES_ONLY:
		debug_grid_init(lnrho, uu_x, uu_y, uu_z, GRID_NORMAL);
		for (int k=CZ_BOT; k < CZ_TOP; k++) {
			for (int j=CY_BOT; j < CY_TOP; j++) {
				for (int i=CX_BOT; i < CX_TOP; i++) {
					int idx = i + j*NX + k*NX*NY;
					lnrho[idx] = 0.0;
					uu_x[idx] = 0.0;
					uu_y[idx] = 0.0;
					uu_z[idx] = 0.0;
				}
			}
		}
		break;

	case GRID_NO_FLOAT_ARITHMETIC_ERROR:
		for (int k=0; k < NZ; k++) {
			for (int j=0; j < NY; j++) {
				for (int i=0; i < NX; i++) {
					int idx = i + j*NX + k*NX*NY;
					lnrho[idx] = 8;
					uu_x[idx] = 8;
					uu_y[idx] = 0;
					uu_z[idx] = 6;//8 ja 6 menee tasan
				}
			}
		}
		break;
	case GRID_WAVE:
	for (int k=CZ_BOT; k < CZ_TOP; k++) {
		for (int j=CY_BOT; j < CY_TOP; j++) {
			for (int i=CX_BOT+30; i < CX_BOT+35; i++) {
				idx = i + j*NX + k*NX*NY;
				lnrho[idx] = 2.0;
				lnrho[idx] = AMPL_UU*powf((abs(CX_BOT+35-i) + abs(CX_BOT+30-i)),2.0)/(10.0);
			}
		}
	}
		break;
	case GRID_SMEM_TEST:
		for (int k=CZ_BOT; k < CZ_TOP; k++) {
			for (int j=CY_BOT; j < CY_TOP; j++) {
				for (int i=CX_BOT; i < CX_TOP; i++) {
					int idx = i + j*NX + k*NX*NY;
					lnrho[idx] = i+j-k;
					uu_x[idx] = i;
					uu_y[idx] = j;
					uu_z[idx] = k;
				}
			}
		}
		break;

	case GRID_BOUND_CPY_ZONE_ONLY:
		for (int k=CZ_TOP-3; k < CZ_TOP; k++) {
			for (int j=CY_BOT; j < CY_BOT+3; j++) {
				for (int i=CX_BOT; i < CX_TOP; i++) {
					int idx = i + j*NX + k*NX*NY;
					lnrho[idx] = 10;//rand() % magic_number;
					uu_x[idx] = i; //rand() % magic_number;
					uu_y[idx] = j; //rand() % magic_number;
					uu_z[idx] = k; //rand() % magic_number;
				}
			}
		}
		for (int k=CZ_TOP-3; k < CZ_TOP; k++) {
			for (int j=CY_TOP-3; j < CY_TOP; j++) {
				for (int i=CX_BOT; i < CX_TOP; i++) {
					int idx = i + j*NX + k*NX*NY;
					lnrho[idx] = 20;//rand() % magic_number;
					uu_x[idx] = i; //rand() % magic_number;
					uu_y[idx] = j; //rand() % magic_number;
					uu_z[idx] = k; //rand() % magic_number;
				}
			}
		}
		for (int k=CZ_TOP-3; k < CZ_TOP; k++) {
			for (int j=CY_BOT; j < CY_TOP; j++) {
				for (int i=CX_TOP-3; i < CX_TOP; i++) {
					int idx = i + j*NX + k*NX*NY;
					lnrho[idx] = 40;//rand() % magic_number;
					uu_x[idx] = i; //rand() % magic_number;
					uu_y[idx] = j; //rand() % magic_number;
					uu_z[idx] = k; //rand() % magic_number;
				}
			}
		}
		for (int k=CZ_TOP-3; k < CZ_TOP; k++) {
			for (int j=CY_BOT; j < CY_TOP; j++) {
				for (int i=CX_BOT; i < CX_BOT+3; i++) {
					int idx = i + j*NX + k*NX*NY;
					lnrho[idx] = 50;//rand() % magic_number;
					uu_x[idx] = i; //rand() % magic_number;
					uu_y[idx] = j; //rand() % magic_number;
					uu_z[idx] = k; //rand() % magic_number;
				}
			}
		}
		break;

	case GRID_ALL_ZERO:
		debug_clear_grid(lnrho, uu_x, uu_y, uu_z);
		break;

	default:
		printf("Invalid type in init_grid!\n");
		exit(-1);
	}
}

//Initialize lnrho
void lnrho_init(float* lnrho)
{
	//Init types defined in init.conf
	switch(INITCOND_LNRHO)
	{
	case INITCOND_LNRHO_CONST:
		lnrho_const(lnrho);
		break;
	
	default:
		printf("Unknown init type in lnrho_init!\n");
		exit(-1);	
	}
}

//Initialize velocities
void hydro_init(float* uu_x, float* uu_y, float* uu_z)
{
	//Init types defined in init.conf
	switch(INITCOND_UU)
	{
	case INITCOND_UU_CONST_X:
		uu_const(uu_x);
		break;

        case INITCOND_UU_CONST_Y:
                uu_const(uu_y);
                break;

        case INITCOND_UU_CONST_Z:
                uu_const(uu_z);
                break;

        case INITCOND_UU_CONST_XY:
                uu_const(uu_x);
                uu_const(uu_y);
                break;

        case INITCOND_UU_CONST_XZ:
                uu_const(uu_x);
                uu_const(uu_z);
                break;

        case INITCOND_UU_CONST_YZ:
                uu_const(uu_y);
                uu_const(uu_z);
        	break;

        case INITCOND_UU_CONST_XYZ:
                uu_const(uu_x);
                uu_const(uu_y);
                uu_const(uu_z);
                break;

	case INITCOND_UU_UNIFORM_RANDOM_VEL:
		uniform_random_velocity(uu_x, uu_y, uu_z);
		break;

	case INITCOND_UU_XGAUSSIAN_X:
		xgaussian_wave(uu_x);	
		break;

	case INITCOND_UU_YGAUSSIAN_Y:
		ygaussian_wave(uu_y);	
		break;

	case INITCOND_UU_ZGAUSSIAN_Z:
		zgaussian_wave(uu_z);	
		break;

	case INITCOND_UU_3DGAUSSIAN_X:
		gaussian_3d_ball(uu_x);
		break;

	case INITCOND_UU_3DGAUSSIAN_Y:
		gaussian_3d_ball(uu_y);
		break;

	case INITCOND_UU_3DGAUSSIAN_Z:
		gaussian_3d_ball(uu_z);
		break;

	case INITCOND_UU_KINETIC_EXPLOSION:
		gaussian_radial_explosion(uu_x, uu_y, uu_z);
		break;

	case INITCOND_UU_Y_SIN_KX:
        	sin_kx_wave(uu_y);
		break;

	case INITCOND_UU_Z_SIN_KY:
        	sin_ky_wave(uu_z);
		break;

	case INITCOND_UU_X_SIN_KZ:
        	sin_kz_wave(uu_x);
		break;

	case INITCOND_UU_STEP_PER_X:
                periodic_step_function_x(uu_x);
		break;

	case INITCOND_UU_STEP_PER_Y:
                periodic_step_function_y(uu_y);
		break;

	case INITCOND_UU_STEP_PER_Z:
                periodic_step_function_z(uu_z);
		break;

	default:
		printf("Unknown init type in hydro_init!\n");
		exit(-1);	
	}
}
