#pragma once

#include "../forcing_c.h"

#define TOVEC3(type,arr) ((type){arr[0],arr[1],arr[2]})
#define TOACREAL3(arr) TOVEC3(AcReal3,arr)

typedef struct ForcingParams{

// contains parameters returned by PC routine forcing_pars_hel

    real coef1[3], coef2[3], coef3[3], fda[3], kk[3];
    real phase, fact;

    void Update(AcMeshInfo &config){

      forcing_pars_hel(coef1,coef2,coef3,fda,kk,&phase,&fact);
//printf("phase,fact,kk= %f %f %f %f %f\n", phase, fact, kk[0],kk[1],kk[2]);
      config.real_params[AC_phase] = phase;
      config.real_params[AC_fact] = fact;
      config.real3_params[AC_coef1] = TOACREAL3(coef1);
      config.real3_params[AC_coef2] = TOACREAL3(coef2);
      config.real3_params[AC_coef3] = TOACREAL3(coef3);
      config.real3_params[AC_fda] = TOACREAL3(fda);
      config.real3_params[AC_kk] = TOACREAL3(kk);
      // push to device !?
    }
};
static ForcingParams forcing_params;
