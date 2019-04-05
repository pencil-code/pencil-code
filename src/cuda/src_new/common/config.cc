#include <string.h>

#include "config.h"
#include "errorhandler.h"
#include "utils/utils.h"      // For max, min

#define RFS REAL_FORMAT_SPECIFIER //See datatypes.h for details, basically %f or %lf

void load_value(const char* value, vec3i* dst)
{
    if (sscanf(value, "(%d,%d,%d)", &(dst->x), &(dst->y), &(dst->z)) < 3) 
        CRASH("Malformed value");
}

void load_value(const char* value, vec3r* dst)
{
    if (sscanf(value, "(" RFS "," RFS "," RFS ")", &(dst->x), &(dst->y), &(dst->z)) < 3) 
        CRASH("Malformed value");   
}

void load_value(const char* value, int* dst)
{
    if (sscanf(value, "%d", dst) < 1) 
        CRASH("Malformed value");
}

void load_value(const char* value, real* dst)
{
    if (sscanf(value, RFS, dst) < 1) 
        CRASH("Malformed value");
}

void load_value(const char* value, char* dst)
{
    if (sscanf(value, "%c", dst) < 1) 
        CRASH("Malformed value");
}


void CParamConfig::parse(const char* keyword, const char* value)
{
    if (strcmp(keyword, "nx") == 0)
        load_value(value, &nx);
    else if (strcmp(keyword, "ny") == 0)
        load_value(value, &ny);
    else if (strcmp(keyword, "nz") == 0)
        load_value(value, &nz);
    else
        CRASH("Invalid keyword!");
}

CParamConfig::CParamConfig(){}
 
CParamConfig::CParamConfig(int mx_, int my_, int mz_, int mw_, int nx_, int ny_, int nz_, int nw_,
                           int l1, int l2, int m1, int m2, int n1, int n2, 
                           int nghost_){

      mx = mx_; // nx is the size of the computational domain
      my = my_;
      mz = mz_;

      mxy = mx*my;
      mw  = mw_;
 
      nx = nx_; // nx is the size of the computational domain
      ny = ny_;
      nz = nz_;
      nw = nw_;

      nx_min=l1-1;
      nx_max=l2-1;           // problematic as non-const in PC
      ny_min=m1-1;
      ny_max=m2-1;           // problematic as non-const in PC
      nz_min=n1-1;
      nz_max=n2-1;           // problematic as non-const in PC

      nghost = nghost_; 
}
void CParamConfig::Setup(real dx, real dy, real dz){

      dsx=dx;
      dsy=dy;
      dsz=dz;
      dsmin = min(dsx, min(dsy, dsz));
}

void CParamConfig::compute_missing_values(){

        mx = nx + 2*nghost;
        my = ny + 2*nghost; 
        mz = nz + 2*nghost;
        mw = mx*my*mz;
 
        //Bounds for the computational domain, i.e. nx_min <= i < nx_max
        nx_min = nghost;     //Inclusive
        nx_max = nx + nghost;//Exclusive
        ny_min = nghost;
        ny_max = ny + nghost;
        nz_min = nghost;
        nz_max = nz + nghost;

        //Spacing in the grid (TODO read from file)
        /*dsx = DOMAIN_SIZE_X / nx;
        dsy = DOMAIN_SIZE_Y / ny;
        dsz = DOMAIN_SIZE_Z / nz;*/
        dsmin = min(dsx, min(dsy, dsz));
};

void StartConfig::parse(const char* keyword, const char* value)
{
    if (strcmp(keyword, "ampl_uu") == 0)
        load_value(value, &ampl_uu);
    else if (strcmp(keyword, "ampl_lnrho") == 0)
        load_value(value, &ampl_lnrho);
    else
        CRASH("Invalid keyword!");
}

void RunConfig::Setup(real cdt_, real cdtv_, real cs_, real nu_, real eta_)
{
    cdt  = cdt_;
    cdtv = cdtv_;
    
    cs_sound = cs_;
    nu_visc  = nu_;
    eta = eta_;
}

void RunConfig::parse(const char* keyword, const char* value)
{
    if (strcmp(keyword, "max_steps") == 0)
        load_value(value, &max_steps);
    else if (strcmp(keyword, "save_steps") == 0)
        load_value(value, &save_steps);
    else if (strcmp(keyword, "cdt") == 0)
        load_value(value, &cdt);
    else if (strcmp(keyword, "cdtv") == 0)
        load_value(value, &cdtv);
    else if (strcmp(keyword, "nu_visc") == 0)
        load_value(value, &nu_visc);
    else if (strcmp(keyword, "cs_sound") == 0)
        load_value(value, &cs_sound);
    else if (strcmp(keyword, "eta") == 0)
        load_value(value, &eta);
    else if (strcmp(keyword, "relhel") == 0)
        load_value(value, &relhel);
    else if (strcmp(keyword, "slice_axis") == 0)
        load_value(value, &slice_axis);
    else
        CRASH("Invalid keyword!");    
}
