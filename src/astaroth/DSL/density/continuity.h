// checked 17.6.
rhs=0. 

#include "../density/diffusivity.h"

#if LHYDRO
  glnrho = gradient(LNRHO)   // grad(rho) or grad(lnrho)!
  if (ldensity_nolog){
    rhs += - value(RHO)*divergence(UU)
  }
  else{
    rhs += - divergence(UU)
  }
  return rhs - dot(vecvalue(UU), glnrho)
#else
  return rhs
#endif
