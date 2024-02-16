
#include <cuacComplex.h>

// PC-style helical forcing with profiles
// tbc: AcReal3 -> device?
AcReal3 AC_kk
AcReal3 AC_coef1
AcReal3 AC_coef2
AcReal3 AC_coef3
AcReal3 AC_fda

real AC_phase
real AC_fact
real AC_k1_ff

ScalarArray AC_profx_ampl
ScalarArray AC_profy_ampl
ScalarArray AC_profz_ampl
ScalarArray AC_profx_hel
ScalarArray AC_profy_hel
ScalarArray AC_profz_hel

int AC_iforcing_zsym

// PC-style helical forcing with support for profiles
AcReal3
pcforcing(ScalarArray profx_ampl, ScalarArray profy_ampl, ScalarArray profz_ampl, 
          ScalarArray profx_hel, ScalarArray profy_hel, ScalarArray profz_hel)
{
    AcReal3 pos = (AcReal3){(globalVertexIdx.x - AC_nx_min) * AC_dsx,
                            (globalVertexIdx.y - AC_ny_min) * AC_dsy,
                            (globalVertexIdx.z - AC_nz_min) * AC_dsz}

    acComplex fx = AC_fact * exp(acComplex(0.0, AC_kk.x * AC_k1_ff * pos.z + AC_phase))
    acComplex fy = exp(acComplex(0.0, AC_kk.y * AC_k1_ff * pos.y))
    acComplex fz

    if (AC_iforcing_zsym == 0) {
        fz = exp(acComplex(0.0, AC_kk.z * AC_k1_ff * pos.z))
    }
    else if (AC_iforcing_zsym == 1) {
        fz = acComplex(cos(AC_kk.z * AC_k1_ff * pos.z), 0.0)
    }
    else if (AC_iforcing_zsym == -1) {
        fz = acComplex(sin(AC_kk.z * AC_k1_ff * pos.z), 0.0)
    }
    else {
       // Failure
    }

    acComplex fxyz = fx * fy * fz
    // TODO recheck indices
    force_ampl    = profx_ampl[vertexIdx.x - NGHOST] * profy_ampl[vertexIdx.y] * profz_ampl[vertexIdx.z]
    prof_hel_ampl = profx_hel[vertexIdx.x - NGHOST] * profy_hel[vertexIdx.y] * profz_hel [vertexIdx.z]

    return force_ampl * AC_fda * (acComplex(AC_coef1.x, prof_hel_ampl * AC_coef2.x) * fxyz).x

/*
    Vector rhs;

    rhs.x = force_ampl * AC_fda.x * (Complex(AC_coef1.x, prof_hel_ampl * AC_coef2.x) * fxyz).x;
    rhs.y = force_ampl * AC_fda.y * (Complex(AC_coef1.y, prof_hel_ampl * AC_coef2.y) * fxyz).x;
    rhs.z = force_ampl * AC_fda.z * (Complex(AC_coef1.z, prof_hel_ampl * AC_coef2.z) * fxyz).x;
*/
}

