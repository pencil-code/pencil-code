//Compile-time PC defines: anyone who reads Configs should also have access to these

//Astaroth requires BOUND_SIZE to be know at compile-time:
//Otherwise 
//  a) It seriously limits the optimizations we could do
//  b) Calling device functions with dynamically allocated memory is wonky (CUDA 6)
//  c) The implementations will break with bound size other than 3 (2017-07-20) 

#pragma once

#ifdef GPU_ASTAROTH

  #include "utils/utils.h"
  #include "../../cparam_c.h"
  #include "../../cdata_c.h"

  #define BOUND_SIZE nghost

  #define DOMAIN_SIZE_X (Lxyz[0])
  #define DOMAIN_SIZE_Y (Lxyz[1])
  #define DOMAIN_SIZE_Z (Lxyz[2])

  #define XORIG (xyz0[0])
  #define YORIG (xyz0[1])
  #define ZORIG (xyz0[2])

#else

  #define BOUND_SIZE 3

  #define DOMAIN_SIZE_X (6.28318530718)
  #define DOMAIN_SIZE_Y (6.28318530718)
  #define DOMAIN_SIZE_Z (6.28318530718)

  #define XORIG (DOMAIN_SIZE_X / 2.0)
  #define YORIG (DOMAIN_SIZE_Y / 2.0)
  #define ZORIG (DOMAIN_SIZE_Z / 2.0)

  #define HYDRO
  #define VISCOSITY
  //#define FORCING
  //#define MAGNETIC

#endif
