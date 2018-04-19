#pragma once
#include "common/datatypes.h"
#include "common/defines.h"

//Purely virtual parent for the config structs
class Config {
    virtual void parse(const char* keyword, const char* value) = 0;
};

class CParamConfig : public Config {

    public:
    int mx; // mx is the size of the computational domain with ghost zones
    int my;
    int mz;
    int mxy, mw;

    int nx; // nx is the size of the computational domain
    int ny;
    int nz;
    int nw;

    int nghost;

    //Smallest and largest indices in the computational domain. 
    int nx_min;
    int nx_max;
    int ny_min;
    int ny_max;
    int nz_min;
    int nz_max;

    real dsx;
    real dsy;
    real dsz;
    real dsmin;         //Min spacing
    
    //Implementation for virtual functions

    //virtual void parse(const char* keyword, const char* value) override;
    virtual void parse(const char* keyword, const char* value); //overrides (c++11-only feature)

    CParamConfig();
    CParamConfig(int mx_, int my_, int mz_, int mw_, int nx_, int ny_, int nz_, int nw_,
                 int l1_, int l2_, int m1_, int m2_, int n1_, int n2_, 
                 int nghost_);
    void Setup(real dx_, real dy_, real dz_);

    void compute_missing_values();
};

class StartConfig : public Config {

    public: 
    real ampl_uu;
    real ampl_lnrho;

    //Implementation for virtual functions
    virtual void parse(const char* keyword, const char* value); //overrides (c++11-only feature)
};

class RunConfig : public Config {

    public: 
    int max_steps;      //Maximum number of steps the program runs
    int save_steps;     //The number of steps between outputting slides to a file
    int diag_steps;     //Then number of steps between ???(TODO)

    char slice_axis='z';

    //real max_time;      //Maximum time simulated per run (see also max_steps)
    //real anim_step;     //
    real cdt;           //
    real cdtv;          //

    //Hydro
    real nu_visc;       //Viscosity

    //Equation of state
    real cs_sound;      //Sound speed

    //Magnetic
    real eta;           //Magnetic diffusivity
  
    //Forcing
    real relhel;        //Helicity of forcing

    //Implementation for virtual functions
    virtual void parse(const char* keyword, const char* value); //overrides (c++11-only feature)
    void Setup(real cdt_, real cdtv_, real cs2_=0., real nu_=0., real eta_=0.);
};
