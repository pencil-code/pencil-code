#include "../freeze_df.h"
 if(lcourant_dt)
 {
  	if (AC_iuu__mod__cdata != 0)    
	{
		if(AC_lkinflow_as_aux__mod__cdata)
		{
			write( F_UU, ac_transformed_pencil_uu)
		}
		else
		{
			write( F_UU,  rk_intermediate(F_UU , DF_UU,  step_num, AC_dt__mod__cdata) )
		}
	}
        if (AC_iuun__mod__cdata != 0)   write( F_UUN,  rk_intermediate(F_UUN , DF_UUN,  step_num, AC_dt__mod__cdata) )
  	if ((AC_ilnrho__mod__cdata + AC_irho__mod__cdata) != 0)  write(F_RHO, rk_intermediate(F_RHO, DF_RHO, step_num, AC_dt__mod__cdata) )
        if ((AC_ilnrhon__mod__cdata + AC_irhon__mod__cdata) != 0)  write(F_RHON, rk_intermediate(F_RHON, DF_RHON, step_num, AC_dt__mod__cdata) )
  	if (AC_iss__mod__cdata != 0)  write( F_SS,  rk_intermediate(F_SS, DF_SS,  step_num, AC_dt__mod__cdata) )
  	if (AC_iaa__mod__cdata != 0) write( F_AA,  rk_intermediate(F_AA , DF_AA,  step_num, AC_dt__mod__cdata) )
	if(lbfield) write(F_BVEC, rk_intermediate(F_BVEC, DF_BVEC,step_num,AC_dt__mod__cdata))
#if LGRAVITATIONAL_WAVES_HTXK
	if (AC_lfirst__mod__cdata && !AC_lsplit_gw_rhs_from_rest_on_gpu__mod__gravitational_waves_htxk)
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
  	if (AC_itt__mod__cdata != 0 || AC_ilntt__mod__cdata != 0)  write( F_TT,  rk_intermediate(F_TT, DF_TT, step_num,AC_dt__mod__cdata))
	if(lchemistry)
	{
		for i in 0:nchemspec
		{
			write(F_CHEMISTRY_SPECIES[i],rk_intermediate(F_CHEMISTRY_SPECIES[i],DF_CHEMISTRY_SPECIES[i],step_num,AC_dt__mod__cdata))
		}
	}
	if(ldustvelocity)
	{
		for i in 0:ndustspec
		{
			write(F_DUST_VELOCITY[i],rk_intermediate(F_DUST_VELOCITY[i],DF_DUST_VELOCITY[i],step_num,AC_dt__mod__cdata))
		}
	}
	if(ldustdensity)
	{
		for i in 0:ndustspec
		{
			write(F_DUST_DENSITY[i],rk_intermediate(F_DUST_DENSITY[i],DF_DUST_DENSITY[i],step_num,AC_dt__mod__cdata))
		}
	}
#if LALPHADISK
	write(F_SIGMA,rk_intermediate(F_SIGMA,DF_SIGMA,step_num,AC_dt__mod__cdata))
#endif
#if LAXIONSU2BACK
	if(AC_iaxi_psi__mod__axionsu2back != 0)    write(F_AXI_PSI     ,rk_intermediate(F_AXI_PSI    ,DF_AXI_PSI   ,step_num,AC_dt__mod__cdata)) 
	if(AC_iaxi_psidot__mod__axionsu2back != 0) write(F_AXI_PSIDOT  ,rk_intermediate(F_AXI_PSIDOT ,DF_AXI_PSIDOT,step_num,AC_dt__mod__cdata)) 
	if(AC_iaxi_impsi__mod__axionsu2back != 0)    write(F_AXI_IMPSI   , rk_intermediate(F_AXI_IMPSI    ,DF_AXI_IMPSI    ,step_num,AC_dt__mod__cdata)) 
	if(AC_iaxi_impsidot__mod__axionsu2back != 0) write(F_AXI_IMPSIDOT, rk_intermediate(F_AXI_IMPSIDOT ,DF_AXI_IMPSIDOT ,step_num,AC_dt__mod__cdata))  

	if(AC_iaxi_tr__mod__axionsu2back != 0)     write(F_AXI_TR      ,rk_intermediate(F_AXI_TR     ,DF_AXI_TR    ,step_num,AC_dt__mod__cdata)) 
	if(AC_iaxi_trdot__mod__axionsu2back != 0)  write(F_AXI_TRDOT   ,rk_intermediate(F_AXI_TRDOT  ,DF_AXI_TRDOT ,step_num,AC_dt__mod__cdata)) 
	if(AC_iaxi_imtr__mod__axionsu2back != 0)     write(F_AXI_IMTR    , rk_intermediate(F_AXI_IMTR     ,DF_AXI_IMTR     ,step_num,AC_dt__mod__cdata)) 
	if(AC_iaxi_imtrdot__mod__axionsu2back != 0)  write(F_AXI_IMTRDOT , rk_intermediate(F_AXI_IMTRDOT  ,DF_AXI_IMTRDOT  ,step_num,AC_dt__mod__cdata)) 
	
	if(AC_iaxi_tl__mod__axionsu2back != 0)        write(F_AXI_TL        , rk_intermediate(F_AXI_TL       ,DF_AXI_TL       ,step_num,AC_dt__mod__cdata)) 
	if(AC_iaxi_tldot__mod__axionsu2back != 0)     write(F_AXI_TLDOT     , rk_intermediate(F_AXI_TLDOT    ,DF_AXI_TLDOT    ,step_num,AC_dt__mod__cdata))
	if(AC_iaxi_imtl__mod__axionsu2back != 0)      write(F_AXI_IMTL      , rk_intermediate(F_AXI_IMTL     ,DF_AXI_IMTL     ,step_num,AC_dt__mod__cdata))
	if(AC_iaxi_imtldot__mod__axionsu2back != 0)   write(F_AXI_IMTLDOT   , rk_intermediate(F_AXI_IMTLDOT  ,DF_AXI_IMTLDOT  ,step_num,AC_dt__mod__cdata))

	if(AC_iaxi_psil__mod__axionsu2back != 0)      write(F_AXI_PSIL      , rk_intermediate(F_AXI_PSIL     ,DF_AXI_PSIL     ,step_num,AC_dt__mod__cdata))      
	if(AC_iaxi_psildot__mod__axionsu2back != 0)   write(F_AXI_PSILDOT   , rk_intermediate(F_AXI_PSILDOT  ,DF_AXI_PSILDOT  ,step_num,AC_dt__mod__cdata))
	if(AC_iaxi_impsil__mod__axionsu2back != 0)    write(F_AXI_IMPSIL    , rk_intermediate(F_AXI_IMPSIL   ,DF_AXI_IMPSIL   ,step_num,AC_dt__mod__cdata))  
	if(AC_iaxi_impsildot__mod__axionsu2back != 0) write(F_AXI_IMPSILDOT , rk_intermediate(F_AXI_IMPSILDOT,DF_AXI_IMPSILDOT,step_num,AC_dt__mod__cdata))
#endif

#if LBACKREACT_INFL
	if(AC_iinfl_phi__mod__backreact_infl != 0)  write(F_INFL_PHI ,rk_intermediate(F_INFL_PHI ,DF_INFL_PHI ,step_num,AC_dt__mod__cdata))
	if(AC_iinfl_dphi__mod__backreact_infl != 0) write(F_INFL_DPHI,rk_intermediate(F_INFL_DPHI,DF_INFL_DPHI,step_num,AC_dt__mod__cdata))
#endif

#if LKLEIN_GORDON
	if(AC_iphi__mod__klein_gordon != 0)  write(F_PHI ,rk_intermediate(F_PHI ,DF_PHI ,step_num,AC_dt__mod__cdata))
	if(AC_idphi__mod__klein_gordon != 0) write(F_DPHI,rk_intermediate(F_DPHI,DF_DPHI,step_num,AC_dt__mod__cdata))

	if(AC_iphi_up_im__mod__klein_gordon != 0)   write(F_PHI_UP_IM,rk_intermediate(F_PHI_UP_IM,DF_PHI_UP_IM,step_num,AC_dt__mod__cdata))
	if(AC_iphi_down_re__mod__klein_gordon != 0) write(F_PHI_DOWN_RE,rk_intermediate(F_PHI_DOWN_RE,DF_PHI_DOWN_RE,step_num,AC_dt__mod__cdata))
	if(AC_iphi_down_im__mod__klein_gordon != 0) write(F_PHI_DOWN_IM,rk_intermediate(F_PHI_DOWN_IM,DF_PHI_DOWN_IM,step_num,AC_dt__mod__cdata))

	if(AC_idphi_up_im__mod__klein_gordon != 0)   write(F_DPHI_UP_IM,rk_intermediate(F_DPHI_UP_IM,    DF_DPHI_UP_IM,step_num,AC_dt__mod__cdata))
	if(AC_idphi_down_re__mod__klein_gordon != 0) write(F_DPHI_DOWN_RE,rk_intermediate(F_DPHI_DOWN_RE,DF_DPHI_DOWN_RE,step_num,AC_dt__mod__cdata))
	if(AC_idphi_down_im__mod__klein_gordon != 0) write(F_DPHI_DOWN_IM,rk_intermediate(F_DPHI_DOWN_IM,DF_DPHI_DOWN_IM,step_num,AC_dt__mod__cdata))

	if(AC_ipsi__mod__klein_gordon != 0)  write(F_PSI ,rk_intermediate(F_PSI ,DF_PSI ,step_num,AC_dt__mod__cdata))
	if(AC_idpsi__mod__klein_gordon != 0) write(F_DPSI,rk_intermediate(F_DPSI,DF_DPSI,step_num,AC_dt__mod__cdata))
#endif

#if LDISP_CURRENT
	if(AC_igamma__mod__disp_current != 0) write(F_GAMMA,rk_intermediate(F_GAMMA,DF_GAMMA,step_num,AC_dt__mod__cdata)) 
	if(AC_ia0__mod__disp_current != 0)    write(F_A0   ,rk_intermediate(F_A0   ,DF_A0   ,step_num,AC_dt__mod__cdata))
	if(AC_irhoe__mod__disp_current  != 0) write(F_RHOE ,rk_intermediate(F_RHOE ,DF_RHOE ,step_num,AC_dt__mod__cdata))
	if(AC_idiva_name__mod__disp_current  != 0) write(F_DIVA_NAME,rk_intermediate(F_DIVA_NAME,DF_DIVA_NAME,step_num,AC_dt__mod__cdata))
	if(AC_iex__mod__cdata != 0)    write(F_EVEC   ,rk_intermediate(F_EVEC   ,DF_EVEC   ,step_num,AC_dt__mod__cdata)) 
#endif
#if LCHIRAL_MHD
        if(AC_imu5__mod__chiral_mhd != 0) write(F_MU5, rk_intermediate(F_MU5, DF_MU5, step_num,AC_dt__mod__cdata))
        if(AC_imus__mod__chiral_mhd != 0) write(F_MUS, rk_intermediate(F_MUS, DF_MUS, step_num,AC_dt__mod__cdata))
#endif
	if(lpolymer)
	{
		for i in 0:6
		{
			write(F_POLY[i],rk_intermediate(F_POLY[i],DF_IPOLY__MOD__CDATA[i],step_num,AC_dt__mod__cdata))
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

