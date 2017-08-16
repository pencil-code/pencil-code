#pragma once
#include "datatypes.h"
#include "utils/utils.h"//For max()


//Compile-time PC defines: anyone who reads Configs should also have access to these
#define BOUND_SIZE (3)

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
    virtual void parse(const char* keyword, const char* value) override;

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
        const real domain_size = 6.28318530718; 
        dsx = domain_size / nx;
        dsy = domain_size / ny;
        dsz = domain_size / nz;
        dsmin = domain_size / max(nx, max(ny, nz));
    }
} CParamConfig;


typedef struct StartConfig : public Config {
    real ampl_uu;

    //Implementation for virtual functions
    virtual void parse(const char* keyword, const char* value) override;
} StartConfig;


typedef struct RunConfig : public Config {
    int max_steps;      //Maximum number of steps the program runs
    int save_steps;     //The number of steps between outputting slides to a file
    //int diag_steps;     //Then number of steps between ???(TODO)

    char slice_axis = 'z';//TODO read from file


    //real max_time;      //Maximum time simulated per run (see also max_steps)
    //real anim_step;     //

    real cdt;           //
    real cdtv;          //

    real nu_visc;       //
    real cs_sound;      //
    
    //Implementation for virtual functions
    virtual void parse(const char* keyword, const char* value) override;
} RunConfig;









































