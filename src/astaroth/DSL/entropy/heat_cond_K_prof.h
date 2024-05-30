heat_conduction(int step_num) {

// heat conduction for conductivity with profile.

    grad_lnTT = (gamma-1) * gradient(LNRHO)+ cv1 * gradient(SS)
    del2_lnTT = (gamma-1) * laplace(LNRHO) + cv1 * laplace(SS)
    ind_z = vertexIdx.z - NGHOST
    glnThcond = grad_lnTT + dlnhcond_prof[ind_z]                         ! grad ln(T*hcond)

    chi = exp(-value(LNRHO)) * hcond_prof[ind_z]

    reduce_max(step_num==0,chi,AC_maxchi)

    return chi * ( del2_lnTT + dot(grad_lnTT,glnThcond) )

}
