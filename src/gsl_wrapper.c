/* Wrapper C functions for the gsl library */
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>

void  sp_besselj_l_(float* y, int*l, float* x){
   *y =  gsl_sf_bessel_jl(*l,*x);
}
/* ------------------------------------------ */
void  sp_bessely_l_(float *y, int*l, float* x){
  
   *y =  gsl_sf_bessel_yl(*l,*x);
}
/* ------------------------------------------ */
void  sp_harm_real_(float *y, int *l, int *m, float *theta, float *phi){
  float Plm;
  int ell = *l;
  int emm = *m;
  float fi = *phi;
  float x =  cos(*theta);
  if(emm<0){
    Plm = gsl_sf_legendre_sphPlm(ell,-emm,x);
    *y = pow(-1,emm)*Plm*cos(emm*fi);}
  else{
    Plm = gsl_sf_legendre_sphPlm(ell,emm,x);
    *y = (float)pow(-1,emm)*Plm*cos(emm*fi);}
}
/* -------------------------------------------- */
void  sp_harm_imag_(float *y, int *l, int *m, float *theta, float *phi){
  float Plm;
  int ell = *l;
  int emm = *m;
  float fi = *phi;
  float x =  cos(*theta);
  if(emm<0){
    Plm = gsl_sf_legendre_sphPlm(ell,-emm,x);
    *y = pow(-1,emm)*Plm*sin(emm*fi);}
  else{
    Plm = gsl_sf_legendre_sphPlm(ell,emm,x);
    *y = (float)pow(-1,emm)*Plm*sin(emm*fi);}
}


