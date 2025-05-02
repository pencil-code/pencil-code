#pragma once

#include "../forcing_c.h"

typedef struct ForcingParams{

// contains parameters returned by PC routine forcing_pars_hel

    real coef1[3], coef2[3], coef3[3], fda[3], kk[3];
    real phase, fact;

    void Update(){

      if(lforce_helical[0])
      {
      	forcing_pars_hel(coef1,coef2,coef3,fda,kk,&phase,&fact);
//printf("phase,fact,kk= %f %f %f %f %f\n", phase, fact, kk[0],kk[1],kk[2]);

      	Device device = acGridGetDevice();
      	acDeviceLoadScalarUniform(device, STREAM_DEFAULT, AC_phase, phase);
      	acDeviceLoadScalarUniform(device, STREAM_DEFAULT, AC_fact, fact);
      	acDeviceLoadVectorUniform(device, STREAM_DEFAULT, AC_coef1, TOACREAL3(coef1));
      	acDeviceLoadVectorUniform(device, STREAM_DEFAULT, AC_coef2, TOACREAL3(coef2));
      	acDeviceLoadVectorUniform(device, STREAM_DEFAULT, AC_coef3, TOACREAL3(coef3));
      	acDeviceLoadVectorUniform(device, STREAM_DEFAULT, AC_fda, TOACREAL3(fda));
      	acDeviceLoadVectorUniform(device, STREAM_DEFAULT, AC_kk, TOACREAL3(kk));
      }
    }
} ForcingParams ;
static ForcingParams forcing_params;
