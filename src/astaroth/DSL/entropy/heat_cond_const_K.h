// heat conduction for constant conductivity.

    cp1 = 1./cp
    grad_lnrho = gradient(LNRHO)
    grad_ss = gradient(SS)

    first_term = gamma *  cp1 * laplace(SS) + (gamma - 1.) * laplace(LNRHO)
    second_term = gamma * cp1 * grad_ss + (gamma - 1.) * grad_lnrho
    third_term = gamma * (cp1 * grad_ss) + (gamma - 1.)*grad_lnrho

    chi = hcond_Kconst/(exp(value(LNRHO)) * cp)

    if (lupdate_courant_dt)
    {
    	reduce_max(chi,AC_maxchi)
    }

    rhs += cp * chi * (first_term + dot(second_term, third_term))
    chitot += cp*chi
