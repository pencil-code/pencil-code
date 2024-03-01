heat_conduction_const_chi() {
    
    // lnTT=lnTT0+cv1*ss+gamma_m1*(lnrho-lnrho0)

    inv_AC_cv = 1. / AC_cv
    grad_lnrho = gradient(VTXBUF_LNRHO)
    grad_lnT = (AC_gamma-1) * grad_lnrho            + inv_AC_cv * gradient(VTXBUF_ENTROPY)
    del2_lnT = (AC_gamma-1) * laplace(VTXBUF_LNRHO) + inv_AC_cv * laplace(VTXBUF_ENTROPY)

    return AC_cp * AC_chi * ( dot(grad_lnrho+grad_lnT,grad_lnT) + del2_lnT )
}
