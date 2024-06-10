rho1 = 0.
lnrho = 0.
glnrho = real3(0.,0.,0.)
if (ldensity_nolog){
        rho1 = 1/value(RHO)
        glnrho = rho1*gradient(RHO)
        lnrho = log(value(RHO))
}
else{
        rho1 = exp(-value(LNRHO))
        glnrho = gradient(LNRHO)
        lnrho = value(LNRHO)
}
#if LMAGNETIC
    jj =  (gradient_of_divergence(AA) - veclaplace(AA))/mu0
    bb = curl(AA)
    advec2 = norm2(bb)*rho1/mu0
#else
    advec2 = 0.
#endif
#if LENTROPY
    cs2 = cs20 * exp(gamma * value(SS)/cp + gamma_m1 * (lnrho - lnrho0))
    advec2 = advec2 + cs2
#endif
    uu=vecvalue(UU)
    reduce_max(step_num==0, abs(uu.x/AC_dsx+uu.y/AC_dsy+uu.z/AC_dsz) + sqrt(advec2)/AC_dsx, AC_maxadvec)
    rhs=real3(0.,0.,0.)
#if LVISCOSITY
#include "../hydro/viscosity.h"
#endif
    return rhs 
           - gradients(UU) * uu
#if LENTROPY
           - cs2 * (gradient(SS)/cp + glnrho)
#else
	   - cs20 * glnrho
#endif
#if LMAGNETIC
           + rho1 * cross(jj,bb)
#endif
#if LGRAVITY
           + gravz_zpencil[vertexIdx.z-NGHOST]
#endif
