#if LSHOCK
// checked 18.6.
#include "../stdlib/smooth_max.h"

// Get divergence of velocity.
divu_shock()
{
    // con_bias \elementof [0,1] 
    //          = 0 discards values which do not contain negative divergence
    //          = 1 absolute divergence applies
    //          < 1 divergence used with reduced factor con_bias**shock_div_pow 
    divu = 0.0
    if(AC_low_order_divu__mod__shock)
    {
	divu = divergence_2nd(UU)
    }
    else
    {
    	divu = divergence(UU)
    }
//if (blockIdx.x==0 && blockIdx.y==0 && threadIdx.x==8 && threadIdx.y==8 && threadIdx.z==8) {print("div= %e \n",divu)}

    tmp = 0.
    if (AC_lconvergence_only__mod__shock) {
      tmp = max(0.,-divu)
    } 
    else if (AC_lconvergence_bias__mod__shock){
        tmp = max(0.,-divu) + max(0.,AC_con_bias__mod__shock*divu)
    }
    else {
      tmp = abs(divu)
    }
    if (AC_shock_div_pow__mod__shock != 1.) {tmp = AC_dt_div_pow__mod__shock * pow(tmp,AC_shock_div_pow__mod__shock)}

    return tmp
}
#endif
