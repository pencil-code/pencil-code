rhs=0. 
glnrho = gradient(LNRHO)
#include "../density/diffusivity.h"
if (ldensity_nolog){
  return rhs - dot(vecvalue(UU), glnrho) - value(RHO)*divergence(UU)
}
else{
  return rhs - dot(vecvalue(UU), glnrho) - divergence(UU)
}

