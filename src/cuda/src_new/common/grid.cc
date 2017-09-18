#include "errorhandler.h"
#include "math.h"
#define INCLUDED_FROM_GRID_NAME_DEFINER
#define INCLUDED_FROM_INIT_TYPE_NAME_DEFINER
#include "grid.h"


void grid_malloc(Grid* g, CParamConfig *cparams)
{
    const size_t grid_size = cparams->mx * cparams->my * cparams->mz;
    const size_t grid_size_bytes = sizeof (real) * grid_size;

    for (int i=0; i < NUM_ARRS; ++i) {
        g->arr[i] = (real*) malloc(grid_size_bytes); 
        if (g->arr[i] == NULL) CRASH("Malloc fail");
        for (int j=0; j < grid_size; ++j) 
            (g->arr[i])[j] = NAN;
    }
}


void grid_free(Grid* g)
{
    for (int i=0; i < NUM_ARRS; ++i) {
        free(g->arr[i]); g->arr[i] = NULL; 
    } 
}


//Initialize all grid values into 0.0 to avoid weird errors
void grid_clear(Grid* g, CParamConfig* cparams)
{
    const int mx = cparams->mx;
    const int my = cparams->my;
    const int mz = cparams->mz;
    
	//Note: initialises also the pad zone and boundaries
    for (int w=0; w < NUM_ARRS; ++w) {
	    for (int k=0; k < mz; k++) {
		    for (int j=0; j < my; j++) {
			    for (int i=0; i < mx; i++) {
				    int idx = i + j*mx + k*mx*my;
                    g->arr[w][idx] = 0;
			    }
		    }
	    }
    }
}


void gaussian_radial_explosion(Grid* grid, CParamConfig* cparams, StartConfig* start_params)	
{
    real* uu_x = grid->arr[UUX];
    real* uu_y = grid->arr[UUY];
    real* uu_z = grid->arr[UUZ];

    const int mx = cparams->mx;
    const int my = cparams->my;
    const int mz = cparams->mz; 

    const int nx = cparams->nx;
    const int ny = cparams->ny;
    const int nz = cparams->nz;     

    const int nx_min = cparams->nx_min;
    const int nx_max = cparams->nx_max;
    const int ny_min = cparams->ny_min;
    const int ny_max = cparams->ny_max;
    const int nz_min = cparams->nz_min;
    const int nz_max = cparams->nz_max;

    const real DX = cparams->dsx;
    const real DY = cparams->dsy;
    const real DZ = cparams->dsz;
    const real xorig = XORIG-0.000001;//TODO triggers the "PHI WRONG" etc errors if this is the exact orig for some reason
    const real yorig = YORIG-0.000001;
    const real zorig = ZORIG-0.000001;
    
    const real INIT_LOC_UU_X = 0.0;
    const real INIT_LOC_UU_Y = 0.0;
    const real INIT_LOC_UU_Z = 0.0;

    const real AMPL_UU = start_params->ampl_uu;
    const real UU_SHELL_R = 0.8;
    const real WIDTH_UU = 0.2;

	// Outward explosion with gaussian initial velocity profile. 
   	int idx;
	real xx, yy, zz, rr2, rr, theta, phi; 
	real uu_radial, theta_old;

	theta_old = 0;

	for (int k=nz_min; k < nz_max; k++) {
		for (int j=ny_min; j < ny_max; j++) {
			for (int i=nx_min; i < nx_max; i++) {
				//Calculate the value of velocity in a particular radius.
				idx = i + j*mx + k*mx*my;
				//Determine the coordinates
                                xx = DX*(i-nx_min)-xorig;
				xx = xx-INIT_LOC_UU_X;

                                yy = DY*(j-ny_min)-yorig;
                                yy = yy-INIT_LOC_UU_Y;

                                zz = DZ*(k-nz_min)-zorig;
                                zz = zz-INIT_LOC_UU_Z;

				rr2 = pow(xx,2.0) + pow(yy,2.0) + pow(zz,2.0);
				rr = sqrt(rr2);

				// Origin is different!
				real xx_abs, yy_abs, zz_abs;
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
                                                                printf("Explosion PHI WRONG --: xx = %.3f, yy = %.3f, phi = %.3e, (3.0*M_PI)/2.0 = %.3e\n", 
                                                                xx, yy, phi, (3.0*M_PI)/2.0);
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


void grid_init(Grid* grid, InitType init_type, CParamConfig* cparams, StartConfig* start_params) 
{
    real* lnrho = grid->arr[LNRHO];
    real* uu_x = grid->arr[UUX];//TODO init all grids
    real* uu_y = grid->arr[UUY];
    real* uu_z = grid->arr[UUZ];

    const int mx = cparams->mx;
    const int my = cparams->my;
    const int mz = cparams->mz;    

    const int nx_min = cparams->nx_min;
    const int nx_max = cparams->nx_max;
    const int ny_min = cparams->ny_min;
    const int ny_max = cparams->ny_max;
    const int nz_min = cparams->nz_min;
    const int nz_max = cparams->nz_max;

    const real ampl_uu = start_params->ampl_uu;


	int idx;
	real rnd;
	int magic_number = 1000;
	srand(magic_number);
	grid_clear(grid, cparams);

	switch (init_type){
	case GRID_INCREASING:
        for (int w=0; w < NUM_ARRS; ++w) {
		    for (int k=0; k < mz; k++) {
			    for (int j=0; j < my; j++) {
				    for (int i=0; i < mx; i++) {
					    int idx = i + j*mx + k*mx*my;
                        grid->arr[w][idx] = idx / (long double)(mx*my*mz) + w / (long double) NUM_ARRS;
				    }
			    }
		    }
        }
		break;
	case GRID_DECREASING:
        for (int w=0; w < NUM_ARRS; ++w) {
		    for (int k=0; k < mz; k++) {
			    for (int j=0; j < my; j++) {
				    for (int i=0; i < mx; i++) {
					    int idx = i + j*mx + k*mx*my;
                        grid->arr[w][idx] = 1.0l - idx / (long double)(mx*my*mz) + w / (long double) NUM_ARRS;
				    }
			    }
		    }
        }
		break;

	case GRID_RANDOM: //Limit the range to 0...1 so that we do not compute with wildly different values (so that we do not lose too much fp precision)
        for (int w=0; w < NUM_ARRS; ++w) {
		    for (int k=0; k < mz; k++) {
			    for (int j=0; j < my; j++) {
				    for (int i=0; i < mx; i++) {
					    int idx = i + j*mx + k*mx*my;
                        grid->arr[w][idx] =  (real)(rand() % magic_number) / magic_number;
				    }
			    }
		    }
        }
		break;

	case GRID_COMP_DOMAIN_ONLY:
        for (int w=0; w < NUM_ARRS; ++w) {
		    for (int k=0; k < mz; k++) {
			    for (int j=0; j < my; j++) {
				    for (int i=0; i < mx; i++) {
					    int idx = i + j*mx + k*mx*my;
					    grid->arr[w][idx] = NAN;
				    }
			    }
		    }
		    for (int k=nz_min; k < nz_max; k++) {
			    for (int j=ny_min; j < ny_max; j++) {
				    for (int i=nx_min; i < nx_max; i++) {
					    int idx = i + j*mx + k*mx*my;
					    grid->arr[w][idx] = (real)(rand() % magic_number) / magic_number;
				    }
			    }
		    }
        }
		break;

	case GRID_BOUND_ZONES_ONLY:
		grid_init(grid, GRID_INCREASING, cparams, start_params);
        for (int w=0; w < NUM_ARRS; ++w) {
		    for (int k=nz_min; k < nz_max; k++) {
			    for (int j=ny_min; j < ny_max; j++) {
				    for (int i=nx_min; i < nx_max; i++) {
					    int idx = i + j*mx + k*mx*my;
                        grid->arr[w][idx] = 0;
				    }
			    }
		    }
        }
		break;

	case GRID_NO_FLOAT_ARITHMETIC_ERROR:
        for (int w=0; w < NUM_ARRS; ++w) {
		    for (int k=0; k < mz; k++) {
			    for (int j=0; j < my; j++) {
				    for (int i=0; i < mx; i++) {
					    int idx = i + j*mx + k*mx*my;
                        grid->arr[w][idx] = w+1;
				    }
			    }
		    }
        }
		break;
	case GRID_WAVE:
	for (int k=nz_min; k < nz_max; k++) {
		for (int j=ny_min; j < ny_max; j++) {
			for (int i=nx_min+30; i < nx_min+35; i++) {
				idx = i + j*mx + k*mx*my;
				lnrho[idx] = ampl_uu*pow((fabs(nx_min+35-i) + fabs(nx_min+30-i)),2.0)/(10.0);
			}
		}
	}
		break;
	case GRID_RAND_WITH_NEG:
        for (int w=0; w < NUM_ARRS; ++w) {
		    for (int k=nz_min; k < nz_max; k++) {
			    for (int j=ny_min; j < ny_max; j++) {
				    for (int i=nx_min; i < nx_max; i++) {
					    int idx = i + j*mx + k*mx*my;
                        grid->arr[w][idx] = (real)(rand() % magic_number) / magic_number - 0.5;
				    }
			    }
		    }
        }
		break;

	case GRID_ALL_UNIQUE:
        for (int w=0; w < NUM_ARRS; ++w) {
		    for (int k=nz_min; k < nz_max; k++) {
			    for (int j=ny_min; j < ny_max; j++) {
				    for (int i=nx_min; i < nx_max; i++) {
					    int idx = i + j*mx + k*mx*my;
                        grid->arr[w][idx] = NUM_ARRS*idx / (long double) (nx_max*ny_max*nz_max);
				    }
			    }
		    }
        }
		break;

	case GRID_ALL_ZERO:
		grid_clear(grid, cparams);
		break;

	case GRID_GAUSSIAN_RADIAL_EXPL:
		gaussian_radial_explosion(grid, cparams, start_params);
		break;

	default:
		printf("Invalid init_type in init_grid!\n");
		exit(-1);
	}
}


