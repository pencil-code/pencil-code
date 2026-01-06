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
  if (lvisc_nu_profy_bound) {rhs += 2.*nu_y[vertexIdx.y]*contract(traceless_rateof_strain(UU))}

#if LSHOCK
  if (lvisc_nu_shock && lshock_heat) {rhs += nu_shock * SHOCK * divergence(UU) * divergence(UU)}
#endif
#endif
#if LMAGNETIC
  j = (gradient_of_divergence(AA) - laplace(AA))/mu0
  eta_heat=0.
  if (lresi_eta_const){eta_heat+=eta}
  if (lresi_ydep) {eta_heat+=eta_y[vertexIdx.y]}
  rhs += eta_heat * mu0 * dot(j,j)*rho1
#endif

#if LINTERSTELLAR
  #include "../entropy/heatcool.h"
#endif

  rhs /= TT
  gss = gradient(SS)
  del2ss = laplace(SS)
  chitot=0.
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
  if (step_num == 0 && ldt && lcourant_dt)
  {
    dline_1, dxyz_2, dxyz_4, dxyz_6 = get_grid_mn()
    reduce_max(real(chitot*dxyz_2), AC_maxdiffchi)
  }
  return rhs
