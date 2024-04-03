#if LMAGNETIC
    j =  (gradient_of_divergence(AA) - veclaplace(AA))/mu0
    rho1 = exp(-value(LNRHO))
#endif
#if LENTROPY
    cs2 = cs20 * exp(gamma * value(SS) / cp + gamma_m1 * (value(LNRHO) - lnrho0))
#endif

    return - gradients(UU) * vecvalue(UU)
#if LENTROPY
           - cs2 * (gradient(SS)/cp + gradient(LNRHO))
#else
	   - cs20 * gradient(LNRHO)
#endif
#if LMAGNETIC
           + rho1 * cross(j, curl(AA))
#endif
           + nu * (veclaplace(UU) + (1./3.) * gradient_of_divergence(UU) + 2. * stress_tensor(UU) * gradient(LNRHO))
           + zeta * gradient_of_divergence(UU)

