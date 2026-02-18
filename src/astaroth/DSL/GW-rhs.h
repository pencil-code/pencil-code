Kernel GW_rhs(real AC_t__mod__cdata,real AC_dt__mod__cdata){
const int step_num = 0
const bool AC_lrmv__mod__cdata = false
  #include "rhs.h"

  #if LGRAVITATIONAL_WAVES_HTXK
  if(AC_lsplit_gw_rhs_from_rest_on_gpu__mod__gravitational_waves_htxk)
  {
  	write(F_STRESS_0,DF_STRESS_0)
  	write(F_STRESS_1,DF_STRESS_1)
  	write(F_STRESS_2,DF_STRESS_2)
  	write(F_STRESS_3,DF_STRESS_3)
  	write(F_STRESS_4,DF_STRESS_4)
  	write(F_STRESS_5,DF_STRESS_5)
  }
  #endif
}
