//checked 18.6.

rho1 = 0.
lnrho = 0.
glnrho = gradient(LNRHO)     // grad(rho) or grad(lnrho)!
grho = real3(0.,0.,0.)
	
if (ldensity_nolog){
    rho1 = 1./value(RHO)
    grho = glnrho
    glnrho = grho*rho1
    lnrho = log(value(RHO))
}
else{
    rho1 = exp(-value(LNRHO))
    grho = glnrho/rho1
    lnrho = value(LNRHO)
}   // v
#if LMAGNETIC
    jj = (gradient_of_divergence(AA) - laplace(AA))/mu0
    bb = curl(AA)
    advec2 = norm2(bb)*rho1/mu0
#else
    advec2 = 0.
#endif
#if LENTROPY
    cs2 = cs20 * exp(value(SS)/cv + gamma_m1 * (lnrho - lnrho0))  //v
    advec2 = advec2 + cs2
#endif
    reduce_max(step_num==0, abs(sum(value(UU)/AC_ds)) + sqrt(advec2)/AC_dsx, AC_maxadvec)
    rhs=real3(0.,0.,0.)
#if LVISCOSITY
    #include "../hydro/viscosity.h"
#endif
#if LGRAVITY
    if (lgravz_gas){
      rhs.z += gravz_zpencil[vertexIdx.z-NGHOST]
    }
#endif
    if (ladvection_velocity) {
      if (lupw_uu){
        rhs -= real3(ugrad_upw(UUX,UU), ugrad_upw(UUY,UU), ugrad_upw(UUZ,UU))
      }
      else{
        rhs -= gradient_tensor(UU) * UU // order?
      }
    }
    if (lpressuregradient_gas) {
#if LENTROPY
      rhs -= cs2 * (gradient(SS)/cp + glnrho)
#else
      if (leos_isothermal){
        rhs -= cs20 * glnrho
      }
    }
#endif
#if LMAGNETIC
    if (llorentzforce) {
      rhs += rho1 * cross(jj,bb)
    }
#endif
    return rhs 
