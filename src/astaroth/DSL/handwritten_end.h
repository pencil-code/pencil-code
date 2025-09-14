#include "../freeze_df.h"
 if(AC_lcourant_dt__mod__cdata)
 {
  	if (AC_iuu__mod__cdata != 0)    write( F_UU,  rk_intermediate(F_UU , DF_UU,  step_num, AC_dt__mod__cdata) )
        if (AC_iuun__mod__cdata != 0)   write( F_UUN,  rk_intermediate(F_UUN , DF_UUN,  step_num, AC_dt__mod__cdata) )
  	if ((AC_ilnrho__mod__cdata + AC_irho__mod__cdata) != 0)  write(F_RHO, rk_intermediate(F_RHO, DF_RHO, step_num, AC_dt__mod__cdata) )
        if ((AC_ilnrhon__mod__cdata + AC_irhon__mod__cdata) != 0)  write(F_RHON, rk_intermediate(F_RHON, DF_RHON, step_num, AC_dt__mod__cdata) )
  	if (AC_iss__mod__cdata != 0)  write( F_SS,  rk_intermediate(F_SS, DF_SS,  step_num, AC_dt__mod__cdata) )
  	if (AC_iaa__mod__cdata != 0) write( F_AA,  rk_intermediate(F_AA , DF_AA,  step_num, AC_dt__mod__cdata) )
#if LGRAVITATIONAL_WAVES_HTXK
	if (AC_lfirst__mod__cdata)
	{
		write(F_STRESS_0,DF_STRESS_0)
		write(F_STRESS_1,DF_STRESS_1)
		write(F_STRESS_2,DF_STRESS_2)
		write(F_STRESS_3,DF_STRESS_3)
		write(F_STRESS_4,DF_STRESS_4)
		write(F_STRESS_5,DF_STRESS_5)
	}
#endif

	if (AC_lfirst__mod__cdata)
	{
		reduce_max(dt1_max__mod__cdata,AC_dt1_max)
	}
#if LCHIRAL
  	if (AC_ixx_chiral__mod__chiral != 0) write(F_XX_CHIRAL,rk_intermediate(F_XX_CHIRAL, DF_XX_CHIRAL,  step_num, AC_dt__mod__cdata) )
  	if (AC_iyy_chiral__mod__chiral != 0) write(F_YY_CHIRAL,rk_intermediate(F_YY_CHIRAL, DF_YY_CHIRAL,  step_num, AC_dt__mod__cdata) )
  	if (AC_izz_chiral__mod__chiral != 0) write(F_ZZ_CHIRAL,rk_intermediate(F_ZZ_CHIRAL, DF_ZZ_CHIRAL,  step_num, AC_dt__mod__cdata) )
#endif
	if (AC_iecr__mod__cdata != 0) write(F_ECR,rk_intermediate(F_ECR,DF_ECR,step_num,AC_dt__mod__cdata))
	if(lchemistry)
	{
		for i in 0:nchemspec
		{
			write(F_CHEMISTRY_SPECIES[i],rk_intermediate(F_CHEMISTRY_SPECIES[i],DF_CHEMISTRY_SPECIES[i],step_num,AC_dt__mod__cdata))
		}
	}
 }
 else
 {
 	maximum_error = 0.0 
 	if (AC_iuu__mod__cdata != 0) maximum_error = rkf4_update(DF_UU,step_num,AC_dt__mod__cdata,ERROR_UU,BETA_UU,F_UU,uux_initial_max,uuy_initial_max,uuz_initial_max,maximum_error,AC_dt_ratio__mod__cdata,AC_dt_epsi__mod__cdata)
 	if (AC_iaa__mod__cdata != 0) maximum_error = rkf4_update(DF_AA,step_num,AC_dt__mod__cdata,ERROR_AA,BETA_AA,F_AA,aax_initial_max,aay_initial_max,aaz_initial_max,maximum_error,AC_dt_ratio__mod__cdata,AC_dt_epsi__mod__cdata)
 	if ((AC_ilnrho__mod__cdata + AC_irho__mod__cdata) != 0) maximum_error = rkf4_update(DF_RHO,step_num,AC_dt__mod__cdata,ERROR_RHO,BETA_RHO,F_RHO,rho_initial_max,maximum_error,AC_dt_ratio__mod__cdata,AC_dt_epsi__mod__cdata)
 	if (AC_iss__mod__cdata != 0) maximum_error = rkf4_update(DF_SS,step_num,AC_dt__mod__cdata,ERROR_SS,BETA_SS,F_SS,ss_initial_max,maximum_error,AC_dt_ratio__mod__cdata,AC_dt_epsi__mod__cdata)
 	if(step_num == 4) reduce_max(maximum_error,AC_maximum_error)
 }

