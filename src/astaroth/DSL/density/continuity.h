// checked 17.6.
rhs=0. 
glnrho = gradient(LNRHO)   // grad(rho) or grad(lnrho)!

#include "../density/diffusivity.h"

#if LHYDRO
  if (ldensity_nolog){
    rhs += - RHO*divergence(UU)
  }
  else{
    rhs += - divergence(UU)
  }
  if (lupw_lnrho){
    return rhs - ugrad_upw(LNRHO,UU)
  }
  else{
    return rhs - dot(UU, glnrho)
  }
#else
  return rhs
#endif
