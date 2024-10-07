// checked 17.6.
rhs=0. 
glnrho = gradient(LNRHO)   // grad(rho) or grad(lnrho)!

#include "../density/diffusivity.h"

#if LHYDRO
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
