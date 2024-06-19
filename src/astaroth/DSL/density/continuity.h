// checked 17.6.
rhs=0. 
glnrho = gradient(LNRHO)   // grad(rho) or grad(lnrho)!

#include "../density/diffusivity.h"

if (ldensity_nolog){
  return rhs - dot(vecvalue(UU), glnrho) - value(RHO)*divergence(UU)
}
else{
  return rhs - dot(vecvalue(UU), glnrho) - divergence(UU)
}
