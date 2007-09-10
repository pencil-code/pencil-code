/* Wrapper C functions for the gsl library */
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <math.h>

/* Choose single or double precision here (typically done from the Makefile) */
#ifdef DOUBLE_PRECISION
#  define REAL double
#  define FINT int		/* should this be long int? */
#  define NBYTES 8
#else
#  define REAL float
#  define FINT int
#  define NBYTES 4
#endif

/* Pick correct number of underscores here (2 for g77 without
   `-fno-second-underscore', 1 for most other compilers).
   Use the `-DFUNDERSC=1' option in the Makefile to set this.
*/
#if (FUNDERSC == 0)
#  define FTNIZE(name) name
#elif (FUNDERSC == 1)
#  define FTNIZE(name) name##_
#else
#  define FTNIZE(name) name##__
#endif

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


