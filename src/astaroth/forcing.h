#pragma once

#include "../forcing_c.h"
#include "submodule/include/astaroth.h"

#define TOVEC3(type,arr) ((type){arr[0],arr[1],arr[2]})
#define TOACREAL3(arr) TOVEC3(AcReal3,arr)

// contains parameters returned by PC routine forcing_const_coefs_hel_pars
typedef struct {

    real coef1[3], coef2[3], coef3[3], fda[3], kk[3];
    real phase, fact;

    void Update(){
      forcing_coefs_hel_pars(&coef1[0],&coef2[0],&coef3[0],&fda[0],&kk[0],&phase,&fact);
    };
    void Load(){
      acLoadDeviceConstant(AC_forcing_phase, phase);
      acLoadDeviceConstant(AC_fact, fact);
      acLoadVectorConstant(AC_coef1, TOACREAL3(coef1));
      acLoadVectorConstant(AC_coef2, TOACREAL3(coef2));
      acLoadVectorConstant(AC_coef3, TOACREAL3(coef3));
      acLoadVectorConstant(AC_fda, TOACREAL3(fda));
      acLoadVectorConstant(AC_kk, TOACREAL3(kk));
    }
} ForcingParams;
