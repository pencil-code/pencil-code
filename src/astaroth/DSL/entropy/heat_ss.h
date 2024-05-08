    lnrho = value(LNRHO)
    rho1 = 1. / exp(-lnrho)
    cv1 = 1. /cv

    j = (gradient_of_divergence(AA) - veclaplace(AA))/mu0
    lnTT = lnTT0+cv1*value(SS)+gamma_m1*(lnrho-lnrho0)

    rhs = (  eta * mu0 * dot(j, j)*rho1
           + 2. * nu * contract(stress_tensor(UU)) 
           + zeta * divergence(UU) * divergence(UU))*exp(-lnTT)

    return -dot(vecvalue(UU), gradient(SS)) + rhs + heat_conduction(step_num)

