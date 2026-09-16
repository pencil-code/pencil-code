// checked 17.6.
rhs=0. 
glnrho = gradient(LNRHO)   // grad(rho) or grad(lnrho)!

#include "../density/diffusivity.h"

#if LHYDRO
  if (ldensity_nolog){
    rhs += - F_RHO*divergence(F_UVEC)
  }
  else{
    rhs += - divergence(F_UVEC)
  }
  if (lupw_lnrho){
    return rhs - ugrad_upw(F_RHO,F_UVEC)
  }
  else{
    return rhs - dot(F_UVEC, glnrho)
  }
#else
  return rhs
#endif
