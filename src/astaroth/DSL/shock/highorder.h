#if LSHOCK
// checked 18.6.
#include "../stdlib/smooth_max.h"

// Get divergence of velocity.
divu_shock()
{
    // con_bias in (0,1) 0 discards values which do not contain negative divergence
    // 1 absolute divergence applies
    // <1 divergence used with reduced factor con_bias**shock_div_pow 
    divu = divergence(UU)
    tmp = max(-divu, con_bias*divu)
    if (shock_div_pow != 1.) {tmp = dt_div_pow * pow(tmp,shock_div_pow)}

    return tmp
}
#endif
