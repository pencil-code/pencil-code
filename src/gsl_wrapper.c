/* Wrapper C functions for the gsl library */
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>

#include "headers_c.h"

void FTNIZE(sp_besselj_l)
     (REAL* y, FINT*l, REAL* x) {
   *y =  gsl_sf_bessel_jl(*l,*x);
}
/* ------------------------------------------ */
void FTNIZE(sp_bessely_l)
     (REAL *y, FINT*l, REAL* x) {
   *y =  gsl_sf_bessel_yl(*l,*x);
}
/* ------------------------------------------ */
void FTNIZE(sp_harm_real)
     (REAL *y, FINT *l, FINT *m, REAL *theta, REAL *phi) {
  REAL Plm;
  FINT ell = *l;
  FINT emm = *m;
  REAL fi = *phi;
  REAL x =  cos(*theta);
  if(emm<0){
    Plm = gsl_sf_legendre_sphPlm(ell,-emm,x);
    *y = pow(-1,emm)*Plm*cos(emm*fi);}
  else{
    Plm = gsl_sf_legendre_sphPlm(ell,emm,x);
    *y = (REAL)pow(-1,emm)*Plm*cos(emm*fi);}
}
/* -------------------------------------------- */
void FTNIZE(sp_harm_imag)
     (REAL *y, FINT *l, FINT *m, REAL *theta, REAL *phi) {
  REAL Plm;
  FINT ell = *l;
  FINT emm = *m;
  REAL fi = *phi;
  REAL x =  cos(*theta);
  if(emm<0){
    Plm = gsl_sf_legendre_sphPlm(ell,-emm,x);
    *y = pow(-1,emm)*Plm*sin(emm*fi);}
  else{
    Plm = gsl_sf_legendre_sphPlm(ell,emm,x);
    *y = (REAL)pow(-1,emm)*Plm*sin(emm*fi);}
}


