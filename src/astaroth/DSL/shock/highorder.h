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
    divu = divergence(UU)
//if (blockIdx.x==0 && blockIdx.y==0 && threadIdx.x==8 && threadIdx.y==8 && threadIdx.z==8) {print("div= %e \n",divu)}

    tmp = 0.
    if (lconvergence_only) {
      tmp = max(0.,-divu)
    } else {
      if (con_bias != 0.) {
        tmp = max(0.,-divu) + max(0.,con_bias*divu)
      } else {
        tmp = abs(divu)
      }
    }

    if (shock_div_pow != 1.) {tmp = dt_div_pow * pow(tmp,shock_div_pow)}

    return tmp
}
#endif
