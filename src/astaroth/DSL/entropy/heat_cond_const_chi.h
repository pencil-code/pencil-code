// heat conduction for constant diffusivity.

/// lnTT=lnTT0+cv1*ss+gamma_m1*(lnrho-lnrho0)
/// del2lnTT=gamma_m1*p%del2lnrho+cv1*p%del2ss
/// del2lnrho=del2rho/rho - grho**2/rho**2
//grad_lnrho = gradient(LNRHO)
if (lheatc_chiconst) {
  glnT = (gamma-1) * glnrho + cv1 * gradient(SS)
  del2lnT = 0.
  
  if (ldensity_nolog){
      del2lnT = (gamma-1)*rho1 * ( laplace(RHO) - norm2(grho)*rho1 ) + cv1 * laplace(SS)
  }
  else
  {
      del2lnT = (gamma-1) * laplace(LNRHO) + cv1 * laplace(SS)
  }
  rhs += cp * chi * ( dot(glnrho+glnT,glnT) + del2lnT )
  //    return cp * chi * ( dot(grad_lnrho+grad_lnT,grad_lnT) + del2_lnT )
  chitot += cp*chi
}
