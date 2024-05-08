heat_conduction(int step_num) {
    
// heat conduction for constant diffusivity.

// lnTT=lnTT0+cv1*ss+gamma_m1*(lnrho-lnrho0)

    cv1 = 1./cv
    grad_lnrho = gradient(LNRHO)
    grad_lnT = (gamma-1) * grad_lnrho     + cv1 * gradient(SS)
    del2_lnT = (gamma-1) * laplace(LNRHO) + cv1 * laplace(SS)

    return cp * chi * ( dot(grad_lnrho+grad_lnT,grad_lnT) + del2_lnT )
}
