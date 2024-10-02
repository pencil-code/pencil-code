//checked 18.6.
  rhs = 0.
  rho1 = 0.
  lnrho = value(LNRHO)       // rho or lnrho
  grho = real3(0.,0.,0.)
  glnrho = gradient(LNRHO)   // grad(rho) or grad(lnrho) 
  
  if (ldensity_nolog){
    lnrho = log(value(RHO))
    rho1 =  1./value(RHO)
    grho = glnrho
    glnrho = grho*rho1
  }
  else
  {
    rho1 =  exp(-lnrho)
    grho = glnrho/rho1
  }
  lnTT = lnTT0+cv1*value(SS)+gamma_m1*(lnrho-lnrho0)
  TT = exp(lnTT)
#if LHYDRO
  rhs +=  2. * nu * contract(traceless_rateof_strain(UU))
        + zeta * rho1 * divergence(UU) * divergence(UU)   // precalculated?
#endif
#if LMAGNETIC
  j = (gradient_of_divergence(AA) - laplace(AA))/mu0
  rhs += eta * mu0 * dot(j,j)*rho1
#endif

#if LINTERSTELLAR
  #include "../entropy/heatcool.h"
#endif

  rhs /= TT
  #include "../entropy/heat_cond_hyper3.h"
  #include "../entropy/heat_cond_const_chi.h"

#if LHYDRO
  rhs += -dot(UU, gradient(SS))
#endif
  return rhs
