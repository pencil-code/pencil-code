#pragma once

#include "../forcing_c.h"
#include "submodule/include/astaroth.h"

#define TOVEC3(type,arr) ((type){arr[0],arr[1],arr[2]})
#define TOACREAL3(arr) TOVEC3(AcReal3,arr)

// contains parameters returned by PC routine forcing_pars_hel
typedef struct {

    real coef1[3], coef2[3], coef3[3], fda[3], kk[3];
    real phase, fact;

    void Update(){

      forcing_pars_hel(coef1,coef2,coef3,fda,kk,&phase,&fact);
//printf("phase,fact,kk= %f %f %f %f %f\n", phase, fact, kk[0],kk[1],kk[2]);
      acLoadDeviceConstant(AC_phase, phase);
      acLoadDeviceConstant(AC_fact, fact);
      acLoadVectorConstant(AC_coef1, TOACREAL3(coef1));
      acLoadVectorConstant(AC_coef2, TOACREAL3(coef2));
      acLoadVectorConstant(AC_coef3, TOACREAL3(coef3));
      acLoadVectorConstant(AC_fda, TOACREAL3(fda));
      acLoadVectorConstant(AC_kk, TOACREAL3(kk));
    }
} ForcingParams;
