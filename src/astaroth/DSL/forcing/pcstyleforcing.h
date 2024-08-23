// PC-style helical forcing with profiles
//MR: checked 8/24

real3 AC_kk
real3 AC_coef1
real3 AC_coef2
real3 AC_coef3
real3 AC_fda

real AC_phase
real AC_fact

// PC-style helical forcing with support for profiles
forcing(){
    real3 pos = grid_position()

//MR: no discrete phases yet considerd!
    complex fx = AC_fact * exp(complex(0.0, AC_kk.x * k1_ff * pos.x + AC_phase))
    complex fy = exp(complex(0.0, AC_kk.y * k1_ff * pos.y))
    complex fz

    if (iforcing_zsym == 0) {
        fz = exp(complex(0.0, AC_kk.z * k1_ff * pos.z))
    }
    else if (iforcing_zsym == 1) {
        fz = complex(cos(AC_kk.z * k1_ff * pos.z), 0.0)
    }
    else if (iforcing_zsym == -1) {
        fz = complex(sin(AC_kk.z * k1_ff * pos.z), 0.0)
    }

    fxyz = fx * fy * fz

    force_ampl    = profx_ampl[vertexIdx.x - NGHOST] * profy_ampl[vertexIdx.y] * profz_ampl[vertexIdx.z]
    prof_hel_ampl = profx_hel [vertexIdx.x - NGHOST] * profy_hel [vertexIdx.y] * profz_hel [vertexIdx.z]

    ptx = complex(AC_coef1.x, prof_hel_ampl * AC_coef2.x) * fxyz
    pty = complex(AC_coef1.y, prof_hel_ampl * AC_coef2.y) * fxyz
    ptz = complex(AC_coef1.z, prof_hel_ampl * AC_coef2.z) * fxyz

    return real3(force_ampl * AC_fda.x * ptx.x, force_ampl * AC_fda.y * pty.x, force_ampl * AC_fda.z * ptz.x)
}
