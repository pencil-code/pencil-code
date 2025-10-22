Field F_CHEMISTRY_SPECIES[nchemspec]
Field F_CHEMISTRY_REACTIONS[nchemspec]
Field AC_cv_r_spec_full__mod__chemistry[nchemspec]
Field rhs_y_full__mod__chemistry[nchemspec]
Field AC_xx_full__mod__chemistry[nchemspec]
Field AC_diff_full__mod__chemistry[nchemspec]
Field AC_diff_full_add__mod__chemistry[nchemspec]

field_order(AC_iux__mod__cdata-1) Field UUX
field_order(AC_iuy__mod__cdata-1) Field UUY
field_order(AC_iuz__mod__cdata-1) Field UUZ
#define F_UX UUX
#define F_UY UUY
#define F_UZ UUZ

#if LHYDRO
field_order(AC_itij__mod__hydro != 0 ? AC_itij__mod__hydro+0-1 : -1) Field F_TIJ_0
field_order(AC_itij__mod__hydro != 0 ? AC_itij__mod__hydro+1-1 : -1) Field F_TIJ_1
field_order(AC_itij__mod__hydro != 0 ? AC_itij__mod__hydro+2-1 : -1) Field F_TIJ_2
field_order(AC_itij__mod__hydro != 0 ? AC_itij__mod__hydro+3-1 : -1) Field F_TIJ_3
field_order(AC_itij__mod__hydro != 0 ? AC_itij__mod__hydro+4-1 : -1) Field F_TIJ_4
field_order(AC_itij__mod__hydro != 0 ? AC_itij__mod__hydro+5-1 : -1) Field F_TIJ_5
#endif

field_order(AC_iunx__mod__cdata-1) Field UUNX
field_order(AC_iuny__mod__cdata-1) Field UUNY
field_order(AC_iunz__mod__cdata-1) Field UUNZ
#define F_UNX UUNX
#define F_UNY UUNY
#define F_UNZ UUNZ

field_order(AC_iglobal_gg__mod__cdata != 0 ? AC_iglobal_gg__mod__cdata-1+0 : -1) Field F_GLOBAL_GX
field_order(AC_iglobal_gg__mod__cdata != 0 ? AC_iglobal_gg__mod__cdata-1+1 : -1) Field F_GLOBAL_GY
field_order(AC_iglobal_gg__mod__cdata != 0 ? AC_iglobal_gg__mod__cdata-1+2 : -1) Field F_GLOBAL_GZ

const Field3 F_GLOBAL_GVEC = {F_GLOBAL_GX,F_GLOBAL_GY,F_GLOBAL_GZ}

field_order(AC_ilnrho__mod__cdata-1) Field RHO
#define LNRHO RHO
#define F_RHO  RHO
#define F_LNRHO F_RHO

field_order(AC_ilnrhon__mod__cdata-1) Field RHON
#define LNRHON RHON
#define F_RHON  RHON
#define F_LNRHON F_RHON


field_order(AC_iax__mod__cdata-1) Field AAX
field_order(AC_iay__mod__cdata-1) Field AAY
field_order(AC_iaz__mod__cdata-1) Field AAZ

field_order(AC_iss__mod__cdata-1) Field F_SS

field_order(AC_iggt__mod__cdata -1) Field F_GGT
field_order(AC_iggx__mod__cdata -1) Field F_GGX
field_order(AC_iggtim__mod__cdata -1) Field F_GGTIM
field_order(AC_iggxim__mod__cdata -1) Field F_GGXIM

field_order(AC_ihht__mod__cdata -1) Field F_HHT
field_order(AC_ihhx__mod__cdata -1) Field F_HHX
field_order(AC_ihhtim__mod__cdata -1) Field F_HHTIM
field_order(AC_ihhxim__mod__cdata -1) Field F_HHXIM

field_order(AC_istress_ij__mod__cdata != 0 ? AC_istress_ij__mod__cdata+0-1 : -1) Field F_STRESS_0
field_order(AC_istress_ij__mod__cdata != 0 ? AC_istress_ij__mod__cdata+1-1 : -1) Field F_STRESS_1
field_order(AC_istress_ij__mod__cdata != 0 ? AC_istress_ij__mod__cdata+2-1 : -1) Field F_STRESS_2
field_order(AC_istress_ij__mod__cdata != 0 ? AC_istress_ij__mod__cdata+3-1 : -1) Field F_STRESS_3
field_order(AC_istress_ij__mod__cdata != 0 ? AC_istress_ij__mod__cdata+4-1 : -1) Field F_STRESS_4
field_order(AC_istress_ij__mod__cdata != 0 ? AC_istress_ij__mod__cdata+5-1 : -1) Field F_STRESS_5

field_order(AC_istresst__mod__cdata-1) Field F_STRESST
field_order(AC_istressx__mod__cdata-1) Field F_STRESSX
field_order(AC_istresstim__mod__cdata-1) Field F_STRESSTIM
field_order(AC_istressxim__mod__cdata-1) Field F_STRESSXIM


field_order(AC_ipoly__mod__cdata != 0 ? AC_ipoly__mod__cdata+0-1 : -1) Field F_POLY_0
field_order(AC_ipoly__mod__cdata != 0 ? AC_ipoly__mod__cdata+1-1 : -1) Field F_POLY_1
field_order(AC_ipoly__mod__cdata != 0 ? AC_ipoly__mod__cdata+2-1 : -1) Field F_POLY_2
field_order(AC_ipoly__mod__cdata != 0 ? AC_ipoly__mod__cdata+3-1 : -1) Field F_POLY_3
field_order(AC_ipoly__mod__cdata != 0 ? AC_ipoly__mod__cdata+4-1 : -1) Field F_POLY_4
field_order(AC_ipoly__mod__cdata != 0 ? AC_ipoly__mod__cdata+5-1 : -1) Field F_POLY_5

const Field F_POLY = [
			F_POLY_0,
			F_POLY_1,
			F_POLY_2,
			F_POLY_3,
			F_POLY_4,
			F_POLY_5
		     ]

field_order(AC_ipoly_fr__mod__cdata-1) Field F_POLY_FR

#define F_AX AAX
#define F_AY AAY
#define F_AZ AAZ

field_order(AC_ishock__mod__cdata-1) Field SHOCK // shock
	    //

#define SS F_SS

Field3 F_DUST_VELOCITY[ndustspec]
Field  F_DUST_DENSITY[ndustspec]
Field  F_DUST_MASS[ndustspec]
Field  F_DUST_ICE_MASS[ndustspec]

Field F_PHIUU
Field F_U0X, F_U0Y, F_U0Z
field_order(AC_ilorentz__mod__cdata-1) Field F_LORENTZ
Field F_GUIJ11,F_GUIJ21,F_GUIJ31,F_GUIJ12,F_GUIJ22,F_GUIJ32,F_GUIJ13,F_GUIJ23,F_GUIJ33
Field F_OX, F_OY, F_OZ
Field F_UU_SPHX, F_UU_SPHY, F_UU_SPHZ
Field F_BB_SPHX, F_BB_SPHY, F_BB_SPHZ
Field F_HLESS
Field F_EOSVAR2
Field F_GLOBAL_CS2
Field F_PP
Field F_SS_B
Field F_RHO_B
Field F_ETH
Field F_GLOBAL_LNRHO0
Field F_GLOBAL_SS0
Field F_HYPVIS
Field F_NUSMAG
Field F_BX, F_BY, F_BZ
Field F_GLOBAL_EXT_BX, F_GLOBAL_EXT_BY, F_GLOBAL_EXT_BZ
Field F_GLOBAL_EXT_AX, F_GLOBAL_EXT_AY, F_GLOBAL_EXT_AZ

Field3 F_GLOBAL_EEXTVEC
Field3 F_GLOBAL_JEXTVEC

field_order(AC_igamma__mod__disp_current-1) Field F_GAMMA
field_order(AC_irhoe__mod__disp_current-1 ) Field F_RHOE
field_order(AC_idiva_name__mod__disp_current-1) Field F_DIVA_NAME
field_order(AC_ia0__mod__disp_current-1 ) Field F_A0
field_order(AC_iex__mod__disp_current-1 ) Field F_EX
field_order(AC_iey__mod__disp_current-1 ) Field F_EY
field_order(AC_iez__mod__disp_current-1 ) Field F_EZ
field_order(AC_iedotx__mod__disp_current-1) Field F_EDOTX
field_order(AC_iedoty__mod__disp_current-1) Field F_EDOTY
field_order(AC_iedotz__mod__disp_current-1) Field F_EDOTZ
field_order(AC_idive__mod__disp_current-1) Field F_DIVE
field_order(AC_isige__mod__disp_current-1) Field F_SIGE
field_order(AC_isigb__mod__disp_current-1) Field F_SIGB

const Field3 F_EDOT = {F_EDOTX,F_EDOTY,F_EDOTZ}
const Field3 F_EVEC = {F_EX,F_EY,F_EZ}
const Field3 F_EDOTVEC         = {F_EDOTX,F_EDOTY,F_EDOTZ}

Field3 F_GLOBAL_AX_EXVEC


Field F_JX,F_JY,F_JZ
Field F_LAM
field_order(AC_itt__mod__cdata != 0 ? AC_itt__mod__cdata-1 : AC_ilntt__mod__cdata-1) Field F_TT
field_order(AC_iyh__mod__cdata-1) Field F_YH
#define TT F_TT
#define LNTT F_TT
#define F_LNTT F_TT
Field F_GLOBAL_HCOND
Field F_SS_RUN_AVER
Field F_ADV_DERX
Field F_ADV_DERY
Field F_ADV_DERZ
Field F_GLOBAL_GLHX,F_GLOBAL_GLHY,F_GLOBAL_GLHZ
Field F_HYPREX, F_HYPREY, F_HYPREZ
Field F_ETAT



field_order(AC_iglobal_glnTT__mod__cdata == 0 ? -1 : 0+AC_iglobal_glnTT__mod__cdata-1) Field F_GLOBAL_GLNTX
field_order(AC_iglobal_glnTT__mod__cdata == 0 ? -1 : 1+AC_iglobal_glnTT__mod__cdata-1) Field F_GLOBAL_GLNTY
field_order(AC_iglobal_glnTT__mod__cdata == 0 ? -1 : 2+AC_iglobal_glnTT__mod__cdata-1) Field F_GLOBAL_GLNTZ
const Field3 F_GLOBAL_GLNTVEC = {F_GLOBAL_GLNTX,F_GLOBAL_GLNTY,F_GLOBAL_GLNTZ}
const Field3 F_AVEC    = {F_AX, F_AY, F_AZ}
#define F_AA F_AVEC
#define AA F_AVEC


Field F_AXTEST
Field F_AYTEST
Field F_AZTEST

const Field3 F_AXTESTVEC = {F_AXTEST,F_AYTEST,F_AZTEST}
#define F_AXTESVEC F_AXTESTVEC

const Field3 F_UVEC    = {F_UX,F_UY,F_UZ}
const Field3 F_UU      = {F_UX,F_UY,F_UZ}
#define F_UU F_UVEC
#define UU F_UU

const Field3 F_UNVEC    = {F_UNX,F_UNY,F_UNZ}
const Field3 F_UUN      = {F_UNX,F_UNY,F_UNZ}
#define F_UUN F_UNVEC
#define UUN F_UUN

const Field3 F_U0VEC   = {F_U0X, F_U0Y, F_U0Z}
const Field3 F_OVEC    = {F_OX, F_OY, F_OZ}
const Field3 F_UU_SPH_VEC  = {F_UU_SPHX, F_UU_SPHY, F_UU_SPHZ}
const Field3 F_UU_SPHVEC   = {F_UU_SPHX, F_UU_SPHY, F_UU_SPHZ}
const Field3 F_BB_SPHVEC   = {F_BB_SPHX, F_BB_SPHY, F_BB_SPHZ}
const Field3 F_BVEC        = {F_BX,F_BY,F_BZ}
const Field3 F_GLOBAL_GLHVEC = {F_GLOBAL_GLHX,F_GLOBAL_GLHY,F_GLOBAL_GLHZ}


const Field3 F_JVEC            = {F_JX,F_JY,F_JZ}
const Field3 F__ADV_DERVEC     = {F_ADV_DERX,F_ADV_DERY,F_ADV_DERZ}
const Field3 F_HYPREVEC        = {F_HYPREX, F_HYPREY, F_HYPREZ}
const Field3 F_GLOBAL_EXT_AVEC = {F_GLOBAL_EXT_AX, F_GLOBAL_EXT_AY, F_GLOBAL_EXT_AZ}

Field3 F_VVEC
Field3 F_FVEC
field_order(AC_iecr__mod__cdata-1) Field F_ECR
#if LCHIRAL
field_order(AC_ixx_chiral__mod__chiral-1) Field F_XX_CHIRAL
field_order(AC_iyy_chiral__mod__chiral-1) Field F_YY_CHIRAL
field_order(AC_izz_chiral__mod__chiral-1) Field F_ZZ_CHIRAL
#endif


not_implemented(message)
{
    print("NOT IMPLEMENTED: %s\n",message)
}


tini_sqrt_div(real numerator, real a, real b)
{
  return 
     (abs(a) <= AC_REAL_MIN*5.0 || abs(b) <= AC_REAL_MIN*5.0)
     ? 0.0
     : numerator/(sqrt(a*b)) 
}
tini_sqrt_div_separate(real numerator, real a, real b)
{
  return 
     (abs(a) <= AC_REAL_MIN*5.0|| abs(b) <= AC_REAL_MIN*5.0)
     ? 0.0
     : numerator/(sqrt(a)*sqrt(b)) 
}


#define AC_maux__mod__cparam maux__mod__cparam
#define iphF_UU F_PHIUU

#define AC_n2__mod__cparam AC_nx+NGHOST_VAL+1 
#define AC_m2__mod__cparam AC_ny+NGHOST_VAL+1 
#define AC_l2__mod__cparam AC_nz+NGHOST_VAL+1 


#define AC_NGHOST_VAL__mod__cparam NGHOST_VAL

#define AC_mx AC_mlocal.x
#define AC_my AC_mlocal.y
#define AC_mz AC_mlocal.z

#define AC_nx AC_nlocal.x
#define AC_ny AC_nlocal.y
#define AC_nz AC_nlocal.z

#define AC_nxgrid AC_ngrid.x
#define AC_nygrid AC_ngrid.y
#define AC_nzgrid AC_ngrid.z

#define AC_dsx AC_ds.x
#define AC_dsy AC_ds.y
#define AC_dsz AC_ds.z

output real AC_df_rho_sum
Field TEST_1, TEST_2,TEST_3,TEST_4
output real AC_maxadvec
global output real AC_maximum_error
const real tini = AC_REAL_MIN*5

Field AC_cp_full__mod__chemistry
Field AC_lambda_full__mod__chemistry
Field AC_cv_full__mod__chemistry
Field AC_tpq_re__mod__gravitational_waves_htxk[6]
Field AC_tpq_im__mod__gravitational_waves_htxk[6]
Field AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[6]
Field AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[6]


field_order(AC_irho_b__mod__cdata-1) Field RHO_B
field_order(AC_iss_b__mod__cdata-1) Field SS_B

field_order(AC_isld_char__mod__cdata-1) Field SLD_CHAR_SPEED
Field SRAD
#define AC_srad__mod__radiation SRAD
field_order(AC_ikapparho__mod__radiation-1) Field F_KAPPARHO
field_order(AC_iqrad__mod__radiation-1) Field F_QRAD

Field AC_cp_full__mod__equationofstate

field_order(AC_imu5__mod__chiral_mhd-1) Field F_MU5
field_order(AC_imus__mod__chiral_mhd-1) Field F_MUS

Field3 BETA_UU
Field3 BETA_AA
Field  BETA_RHO
Field  BETA_SS

Field3 ERROR_UU
Field3 ERROR_AA
Field  ERROR_RHO
Field  ERROR_SS
#include "../stdlib/map.h"
