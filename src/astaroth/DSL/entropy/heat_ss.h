    S = stress_tensor(UU)
    inv_pT = 1. / (exp(value(LNRHO)) * exp(lnT()))
    j = (1. / AC_mu0) * (gradient_of_divergence(AA) - veclaplace(AA))
    RHS = (0) - (0) + AC_eta * AC_mu0 * dot(j, j) +
                       2. * exp(value(LNRHO)) * AC_nu * contract(S) +
                       AC_zeta * exp(value(LNRHO)) * divergence(UU) * divergence(UU)

    return -dot(vecvalue(UU), gradient(SS)) + inv_pT * RHS + heat_conduction()
