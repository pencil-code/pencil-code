// heat conduction for conductivity with profile.

    grad_lnTT = (gamma-1) * gradient(LNRHO)+ cv1 * gradient(SS)
    del2_lnTT = (gamma-1) * laplace(LNRHO) + cv1 * laplace(SS)
    ind_z = vertexIdx.z - NGHOST
    glnThcond = grad_lnTT + dlnhcond_prof[ind_z]                         ! grad ln(T*hcond)

    chi = exp(-value(LNRHO)) * hcond_prof[ind_z]

    rhs += chi * ( del2_lnTT + dot(grad_lnTT,glnThcond) )

    chitot += chi
    if (lupdate_courant_dt)
    {
    	reduce_max(step_num==0,chi,AC_maxchi)
    }
