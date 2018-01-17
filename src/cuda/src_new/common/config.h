#pragma once
#include "datatypes.h"
#include "utils/utils.h"//For max()


//Compile-time PC defines: anyone who reads Configs should also have access to these
#define BOUND_SIZE (3)
#define LFORCING (1)
#define LINDUCTION (1) //What might L be? Just following the convention set by LFORCING...
#define DOMAIN_SIZE_X (6.28318530718)
#define DOMAIN_SIZE_Y (6.28318530718)
#define DOMAIN_SIZE_Z (6.28318530718)

#define XORIG (DOMAIN_SIZE_X / 2.0)
#define YORIG (DOMAIN_SIZE_Y / 2.0)
#define ZORIG (DOMAIN_SIZE_Z / 2.0)

#define T_STOP_FORCING (1.0)
#define FORCING (1e-4)
#define KK1 (4.5)
#define KK2 (5.5)
#define KMAX (10.0)
#define DKX (1.0)
#define DKY (1.0)
#define DKZ (1.0)   

//GPU defines
//#define NUM_DEVICES (1)//Deprecated, determined now at runtime

//Purely virtual parent for the config structs
typedef struct Config {
    virtual void parse(const char* keyword, const char* value)  = 0;
} Config;


typedef struct CParamConfig : public Config {
    int nx, ny, nz; //The grid dimensions of the computational domain as seen by the current processor (i.e. the block stored in SDRAM)
    int mx, my, mz; //The grid dimensions as seen by the current processor (i.e. the block stored in SDRAM)

    int nx_min, nx_max; //Smallest x index in the computational domain (l1 and l2 in PC)
    int ny_min, ny_max; //Smallest y index in the computational domain (m1 and m2 in PC)
    int nz_min, nz_max; //Smallest z index in the computational domain (n1 and n2 in PC)

    real dsx, dsy, dsz; //Spacing (Old Astaroth names: DX, DY and DZ)
    real dsmin;         //Min spacing
    
    //Implementation for virtual functions
    virtual void parse(const char* keyword, const char* value); //overrides (c++11-only feature)

    void compute_missing_values() {
        mx = nx + 2*BOUND_SIZE;
        my = ny + 2*BOUND_SIZE; 
        mz = nz + 2*BOUND_SIZE;
    
        //Bounds for the computational domain, i.e. nx_min <= i < nx_max
        nx_min = BOUND_SIZE;     //Inclusive
        nx_max = nx + BOUND_SIZE;//Exclusive
        ny_min = BOUND_SIZE;
        ny_max = ny + BOUND_SIZE;
        nz_min = BOUND_SIZE;
        nz_max = nz + BOUND_SIZE;

        //Spacing in the grid (TODO read from file)
        dsx = (real)(DOMAIN_SIZE_X / nx);
        dsy = (real)(DOMAIN_SIZE_Y / ny);
        dsz = (real)(DOMAIN_SIZE_Z / nz);
        dsmin = min(dsx, min(dsy, dsz));
    }
} CParamConfig;


typedef struct StartConfig : public Config {
    real ampl_uu;
    real ampl_lnrho;

    //Implementation for virtual functions
    virtual void parse(const char* keyword, const char* value); //overrides (c++11-only feature)
} StartConfig;


typedef struct RunConfig : public Config {
    int max_steps;      //Maximum number of steps the program runs
    int save_steps;     //The number of steps between outputting slides to a file
    //int diag_steps;     //Then number of steps between ???(TODO)

    char slice_axis;


    //real max_time;      //Maximum time simulated per run (see also max_steps)
    //real anim_step;     //

    real cdt;           //
    real cdtv;          //

    real nu_visc;       //
    real cs_sound;      //

    //Magnetic
    real eta;           //Magnetic diffusivity
    
    //Forcing
    real relhel;        //Helicity of forcing

    //Implementation for virtual functions
    virtual void parse(const char* keyword, const char* value); //overrides (c++11-only feature)
} RunConfig;









































