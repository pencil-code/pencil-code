rhs=0. 
glnrho = real3(0.,0.,0.)
grho = real3(0.,0.,0.)
//if (ldensity_nolog){
  grho = gradient(RHO)
  glnrho = gradient(RHO)/value(RHO)
  #include "../density/diffusivity.h"
  return rhs - dot(vecvalue(UU), grho) - value(RHO)*divergence(UU)
//}
//else{
//  glnrho = real3(0.,0.,0.)
//  #include "../density/diffusivity.h"
//  return rhs - dot(vecvalue(UU), glnrho) - divergence(UU)
//}
//
