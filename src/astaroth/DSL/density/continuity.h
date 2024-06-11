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

//=======
//glnrho = gradient(LNRHO)   // grad(rho) or grad(lnrho)!
//
//#include "../density/diffusivity.h"
//
//if (ldensity_nolog){
//  return rhs - dot(vecvalue(UU), glnrho) - value(RHO)*divergence(UU)
//}
//else{
//  return rhs - dot(vecvalue(UU), glnrho) - divergence(UU)
//}
//>>>>>>> b63cfac73 (MR: minor)
