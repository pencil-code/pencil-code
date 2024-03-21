    j =  (gradient_of_divergence(AA) - veclaplace(AA))/mu0
    rho1 = exp(-value(LNRHO))

    cs2 = cs20 * exp(gamma * value(SS) / cp + gamma_m1 * (value(LNRHO) - lnrho0))

    return - gradients(UU) * vecvalue(UU)
           - cs2 * (gradient(SS)/cp + gradient(LNRHO))    //???
           + rho1 * cross(j, curl(AA))
           + nu * (veclaplace(UU) + (1./3.) * gradient_of_divergence(UU) + 2. * stress_tensor(UU) * gradient(LNRHO))
           + zeta * gradient_of_divergence(UU)                                                 // | overloaded?

