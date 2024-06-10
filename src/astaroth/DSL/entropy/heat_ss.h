rhs = 0.
rho1 = 0.
lnrho = 0.
glnrho = real3(0.,0.,0.)
if (ldensity_nolog){
    lnrho = log(value(RHO))
    rho1 =  1./value(RHO)
    glnrho = gradient(RHO)/value(RHO)
}
else
{
    lnrho = value(LNRHO)
    rho1 =  exp(-lnrho)
    glnrho = gradient(LNRHO)
}
    cv1 = 1./cv
    lnTT = lnTT0+cv1*value(SS)+gamma_m1*(lnrho-lnrho0)

    rhs +=  2. * nu * contract(stress_tensor(UU))
          + zeta * divergence(UU) * divergence(UU)

#if LMAGNETIC
    j = (gradient_of_divergence(AA) - veclaplace(AA))/mu0
    rhs += eta * mu0 * dot(j, j)*rho1
#endif

#if LINTERSTELLAR
    #include "../entropy/heatcool.h"
    rhs += heatcool
#endif

    rhs *= exp(-lnTT)
    #include "../entropy/heat_cond_hyper3.h"
    return -dot(vecvalue(UU), gradient(SS)) + rhs //+ heat_conduction(step_num)
