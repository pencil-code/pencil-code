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
#if LFORCING
forcing(){
    real3 pos = grid_position()

//MR: no discrete phases yet considerd!
    complex fx = AC_fact * exp(complex(0.0, AC_kk.x * AC_k1_ff__mod__forcing * pos.x + AC_phase))
    complex fy = exp(complex(0.0, AC_kk.y * AC_k1_ff__mod__forcing * pos.y))
    complex fz

    if (AC_iforcing_zsym__mod__forcing == 0) {
        fz = exp(complex(0.0, AC_kk.z * AC_k1_ff__mod__forcing * pos.z))
    }
    else if (AC_iforcing_zsym__mod__forcing == 1) {
        fz = complex(cos(AC_kk.z * AC_k1_ff__mod__forcing * pos.z), 0.0)
    }
    else if (AC_iforcing_zsym__mod__forcing == -1) {
        fz = complex(sin(AC_kk.z * AC_k1_ff__mod__forcing * pos.z), 0.0)
    }

    fxyz = fx * fy * fz

    force_ampl    = AC_profx_ampl__mod__forcing[vertexIdx.x - NGHOST] * AC_profy_ampl__mod__forcing[vertexIdx.y] * AC_profz_ampl__mod__forcing[vertexIdx.z]
    AC_prof_hel_ampl__mod__forcing = AC_profx_hel__mod__forcing [vertexIdx.x - NGHOST] * AC_profy_hel__mod__forcing [vertexIdx.y] * AC_profz_hel__mod__forcing [vertexIdx.z]

    ptx = complex(AC_coef1.x, AC_prof_hel_ampl__mod__forcing * AC_coef2.x) * fxyz
    pty = complex(AC_coef1.y, AC_prof_hel_ampl__mod__forcing * AC_coef2.y) * fxyz
    ptz = complex(AC_coef1.z, AC_prof_hel_ampl__mod__forcing * AC_coef2.z) * fxyz

    return real3(force_ampl * AC_fda.x * ptx.x, force_ampl * AC_fda.y * pty.x, force_ampl * AC_fda.z * ptz.x)
}
#else
forcing(){return real3(0.0,0.0,0.0)}
#endif

