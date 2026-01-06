//checked 18.6.
if (lheatc_hyper3ss){
      rhs += chi_hyper3 * del6(SS)
}

#if LSHOCK
if (lheatc_shock){

  del2lnrho = laplace(LNRHO)                // laplace(rho) or laplace(lnrho)
  glnTT = cv1*gradient(SS)+gamma_m1*glnrho

  if (ldensity_nolog){
    del2lnrho = (del2lnrho-norm2(grho)/value(RHO))/value(RHO)
  }
  chi = chi_shock * value(SHOCK) * cv1
  rhs += chi_shock * (value(SHOCK) * (cv1*laplace(SS) + gamma_m1*del2lnrho + dot(glnrho+glnTT,glnTT)) + dot(gradient(SHOCK),glnTT))
  chitot += chi
}
#endif
