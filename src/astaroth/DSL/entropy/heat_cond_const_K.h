heat_conduction() {

// heat conduction for constant conductivity.

    cp1 = 1./cp
    grad_ln_chi = -gradient(LNRHO)

    first_term = gamma *  cp1 * laplace(SS) + (gamma - 1.) * laplace(LNRHO)
    second_term = gamma * cp1 * gradient(SS) + (gamma - 1.) * gradient(LNRHO)
    third_term = gamma * (cp1 * gradient(SS) + gradient(LNRHO)) + grad_ln_chi

    chi = hcond_Kconst/(exp(value(LNRHO)) * cp)

    return cp * chi * (first_term + dot(second_term, third_term))
}
