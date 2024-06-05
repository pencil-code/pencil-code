// Get divergence of velocity.
//#include "../stdlib/smooth_max.h"

divu_shock()
{
    // con_bias in (0,1) 0 discards values which do not contain negative divergence
    // 1 absolute divergence applies
    // <1 divergence used with reduced factor con_bias**shock_div_pow 
    dt_div_pow = pow(dtfactor,shock_div_pow-1)
    divu = divergence(UU)
    tmp = std::max(-divu, con_bias*divu)
    return dt_div_pow * pow(tmp,shock_div_pow)
}

// Calculate local maximum.
/*max5_shock()
{
    return max5(SHOCK)
}*/
// Apply gaussian smoothing to SHOCK
smooth_shock()
{
}
