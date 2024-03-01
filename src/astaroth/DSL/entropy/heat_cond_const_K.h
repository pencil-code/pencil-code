heat_conduction_const_K() {
    inv_AC_cp = 1. / AC_cp
    grad_ln_chi = -gradient(LNRHO)

    first_term = AC_gamma * inv_AC_cp * laplace(SS) + (AC_gamma - 1.) * laplace(LNRHO)
    second_term = AC_gamma * inv_AC_cp * gradient(SS) + (AC_gamma - 1.) * gradient(LNRHO)
    third_term = AC_gamma * (inv_AC_cp * gradient(SS) + gradient(LNRHO)) + grad_ln_chi

    chi = AC_hcond_Kconst / (exp(value(LNRHO)) * AC_cp)

    return AC_cp * chi * (first_term + dot(second_term, third_term))
}
