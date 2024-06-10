if (lheatc_hyper3ss){
      rhs += chi_hyper3 * del6(SS)
}
#if LSHOCK
del2lnrho = 0.
glnTT = real3(0.,0.,0.)
if (lheatc_shock){
if (ldensity_nolog){
glnTT = cv1*gradient(SS)+gamma_m1*gradient(RHO)/value(RHO)
del2lnrho = laplace(RHO)/value(RHO)-norm2(gradient(RHO))/(value(RHO)*value(RHO))
}
else
{
glnTT = cv1*gradient(SS)+gamma_m1*gradient(LNRHO)
del2lnrho = laplace(LNRHO)
}
/*  if (pretend_lnTT) then
    thdiff=gamma*chi_shock*(p%shock*(p%del2lnrho+g2)+gshockglnTT)
  else
    if (lchi_shock_density_dep) then
      call dot(p%gshock,p%gss,gshockgss)
      call dot(twothird*p%glnrho+p%glnTT,p%gss,g2)
      thdiff=exp(-onethird*p%lnrho)*chi_shock* &
             (p%shock*(exp(twothird*p%lnrho)*p%del2ss+g2)+gshockgss)
    else*/
      rhs += chi_shock*value(SHOCK) * (cv1*laplace(SS) + gamma_m1*del2lnrho + dot(glnrho+glnTT,glnTT) + dot(gradient(SHOCK),glnTT))
}
#endif
