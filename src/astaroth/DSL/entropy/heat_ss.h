    rho = exp(value(LNRHO))
    inv_pT = 1. / rho // * exp(lnT()))
    j = (gradient_of_divergence(AA) - veclaplace(AA))/mu0
    rhs =    eta * mu0 * dot(j, j)
           + 2. * rho * nu * contract(stress_tensor(UU)) 
           + zeta * rho * divergence(UU) * divergence(UU)

    return -dot(vecvalue(UU), gradient(SS)) + inv_pT * rhs + heat_conduction()
