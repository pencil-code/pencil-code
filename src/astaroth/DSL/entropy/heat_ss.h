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
  if (lvisc_nu_const){rhs += 2. * nu * contract(traceless_rateof_strain(UU))}
  if (lvisc_rho_nu_const_bulk){rhs += zeta * rho1 * divergence(UU) * divergence(UU)}   // precalculated?
  if (lvisc_nu_shock && lshock_heat) {rhs += nu_shock * SHOCK * divergence(UU) * divergence(UU)}
#endif
#if LMAGNETIC
  j = (gradient_of_divergence(AA) - laplace(AA))/mu0
  rhs += eta * mu0 * dot(j,j)*rho1
#endif

#if LINTERSTELLAR
  #include "../entropy/heatcool.h"
#endif

  rhs /= TT
  gss = gradient(SS)
  del2ss = laplace(SS)
  maxchi=0.
  //!!!#include "../entropy/heat_cond_hyper3.h"
  //!!!#include "../entropy/heat_cond_const_chi.h"
  #include "../entropy/heat_cond_kramers.h"

#if LHYDRO
  if (lupw_ss)
  {
     rhs += - (dot(UU, gss) - dot(abs(UU),gradient_upwd(SS)))
     //rhs += -ugrad_upwd(SS, UU)
  }
  else
  {
     rhs += - (dot(UU, gss))
  }
#endif
  if (lcourant_dt && ldt && step_num == 0)
  {
  	reduce_max(maxchi, AC_maxchi)
  }
  return rhs
