/* Wrapper C functions for the gsl library (double precision) */
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>

void  sp_besselj_l_(double * y, int*l, double* x){
  
   *y =  gsl_sf_bessel_jl(*l,*x);
}

void  sp_bessely_l_(double *y, int*l, double* x){
  
   *y =  gsl_sf_bessel_yl(*l,*x);
}

void  sp_harm_real_(double *y, int *l, int *m, double *theta, double *phi){
  double Plm;
  int ell = *l;
  int emm = *m;
  double fi = *phi;
  double x =  cos(*theta);
  if(emm<0){
    Plm = gsl_sf_legendre_sphPlm(ell,-emm,x);
    *y = pow(-1,emm)*Plm*cos(emm*fi);}
  else{
    Plm = gsl_sf_legendre_sphPlm(ell,emm,x);
    *y = pow(-1,emm)*Plm*cos(emm*fi);}
}

/* --------------------*/
