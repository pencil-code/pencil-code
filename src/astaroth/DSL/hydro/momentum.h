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
}
    rhs=real3(0.,0.,0.)
#if LVISCOSITY
    #include "../hydro/viscosity.h"
#endif
#if LGRAVITY
    if (lgravz_gas){
      rhs.z += gravz_zpencil[vertexIdx.z]
    }
#endif
    if (ladvection_velocity) {
      if (lupw_uu){
        rhs -= real3(ugrad_upw(UUX,UU), ugrad_upw(UUY,UU), ugrad_upw(UUZ,UU))
      }
      else{
        rhs -= gradient_tensor(UU) * UU
      }
    }
    if (lcoriolis_force) {
        rhs -= 2.*cross(Omegav,value(UU))
    }
    advec2 = 0.
    advec_cs2 = 0.
    if (lpressuregradient_gas) {
#if LENTROPY
      cs2 = cs20 * exp(value(SS)/cv + gamma_m1 * (lnrho - lnrho0))  //v
      rhs -= cs2 * (gradient(SS)/cp + glnrho)
      advec_cs2 = cs2
#else
      if (leos_isothermal){
        rhs -= cs20 * glnrho
        advec_cs2 = cs20
      }
#endif
    }
    //if (!lmaximal_cdtv) {advec_cs2 *= 3.}     // otherwise max!!
    //advec_cs2 *= 3./(AC_dsmin*AC_dsmin)
    advec_cs2 *= 1.*(AC_inv_ds_2.x+AC_inv_ds_2.y+AC_inv_ds_2.z)
    advec2 += advec_cs2
#if LMAGNETIC
    if (llorentzforce) {
      jj = (gradient_of_divergence(AA) - laplace(AA))/mu0
      bb = curl(AA)
      advec2 += norm2(bb/real3(dx,dy,dz))*rho1/mu0
      //advec2 += ((bb.x/dx)**2 + (bb.y/dy)**2 + (bb.z/dz)**2)*rho1/mu0
      rhs += rho1 * cross(jj,bb)
    }
#endif
    if (step_num == 0 && ldt && lcourant_dt)
    {
	max_advec = sum(abs(value(UU))/AC_ds) + sqrt(advec2)
#if LSHEAR
	if(lshear)
	{
		max_advec += abs(uy0[vertexIdx.x-NGHOST]*dy_1[vertexIdx.y])
	}
#endif
    	reduce_max(max_advec, AC_maxadvec)
    }
    return rhs 
