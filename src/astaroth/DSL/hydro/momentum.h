#if LMAGNETIC
    jj =  (gradient_of_divergence(AA) - veclaplace(AA))/mu0
    rho1 = exp(-value(LNRHO))
    bb = curl(AA)
    advec2 = dot(bb,bb)*rho1/mu0
#else
    advec2 = 0.
#endif
#if LENTROPY
    cs2 = cs20 * exp(gamma * value(SS)/cp + gamma_m1 * (value(LNRHO) - lnrho0))
    advec2 = advec2 + cs2
#endif
    uu=vecvalue(UU)
    reduce_max(step_num==0, abs(uu.x/AC_dsx+uu.y/AC_dsy+uu.z/AC_dsz) + sqrt(advec2)/AC_dsx, AC_maxadvec)

    return - gradients(UU) * uu
#if LENTROPY
           - cs2 * (gradient(SS)/cp + gradient(LNRHO))
#else
	   - cs20 * gradient(LNRHO)
#endif
#if LMAGNETIC
           + rho1 * cross(jj,bb)
#endif
           + nu * (veclaplace(UU) + (1./3.) * gradient_of_divergence(UU) + 2. * stress_tensor(UU) * gradient(LNRHO))
           + zeta * gradient_of_divergence(UU)

