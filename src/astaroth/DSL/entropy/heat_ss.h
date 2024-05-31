    lnrho = value(LNRHO)
    rho1 = 1. / exp(-lnrho)
    cv1 = 1. /cv
    lnTT = lnTT0+cv1*value(SS)+gamma_m1*(lnrho-lnrho0)

    rhs =  2. * nu * contract(stress_tensor(UU))
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
    return -dot(vecvalue(UU), gradient(SS)) + rhs + heat_conduction(step_num)
