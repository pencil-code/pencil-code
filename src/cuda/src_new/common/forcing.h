#pragma once
#include "common/datatypes.h"

/*
In standalone Astaroth kk_x etc are arrays containing the random k vectors,
however in the PC version the size of these can be 1 (i.e. nk == 1  && k_idx == 0)
and just do kk_x[0] = the_rand_kk_x_component_from_pencil_code
after each full integration step
*/
typedef struct {
    bool forcing_enabled;
    int nk;//TODO size_t
    int k_idx;
    real* kk_x;
    real* kk_y;
    real* kk_z;
    real kk_part_x;
    real kk_part_y;
    real kk_part_z;
    real phi; 
    real kaver;
    //Could also do these with a single array and enum the members like in Grid or Slice,
    //however there's no particular reason to do so (we're not iterating over forcing params anywhere)
} ForcingParams;
