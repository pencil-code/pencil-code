#include <string.h>

#include "config.h"
#include "errorhandler.h"

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


void StartConfig::parse(const char* keyword, const char* value)
{
    if (strcmp(keyword, "ampl_uu") == 0)
        load_value(value, &ampl_uu);
    else if (strcmp(keyword, "ampl_lnrho") == 0)
        load_value(value, &ampl_lnrho);
    else
        CRASH("Invalid keyword!");
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
    else
        CRASH("Invalid keyword!");    
}




















