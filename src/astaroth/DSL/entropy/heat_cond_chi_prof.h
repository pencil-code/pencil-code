heat_conduction(int step_num, real chi0) { //, real *prof, real *dprof) {

// heat conduction for chi with profile.

    grad_lnrho = gradient(LNRHO)
    grad_lnTT = (gamma-1) * grad_lnrho + gradient(SS)/cv
    grad_ss   = gradient(SS)
    g2 = dot(grad_lnrho+grad_lnTT,grad_ss)

    ind_z = vertexIdx.z - NGHOST
    //return chi0 * ( prof[ind_z]*(laplace(SS)+g2) + dprof[ind_z]*grad_ss.z)
    return chi0 * ( chit_prof_stored[ind_z]*(laplace(SS)+g2) + dchit_prof_stored[ind_z]*grad_ss.z)

}
