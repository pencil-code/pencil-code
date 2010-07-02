/* Dummy routines for gsl */

#include "headers_c.h"

// Define dummy routines with as many (or few) underscores as even the
// most perverse Fortran compilers might append
#define DUMMY_ROUTINE(name) \
  void name() { } \
  void name##_() { } \
  void name##__() { } \

DUMMY_ROUTINE(sp_besselj_l)
DUMMY_ROUTINE(sp_bessely_l)
DUMMY_ROUTINE(cyl_bessel_jnu)
DUMMY_ROUTINE(cyl_bessel_ynu)
DUMMY_ROUTINE(sp_harm_real)
DUMMY_ROUTINE(sp_harm_imag)
DUMMY_ROUTINE(legendre_pl)
