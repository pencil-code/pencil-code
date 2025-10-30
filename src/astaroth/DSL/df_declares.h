real3 DF_UVEC          = rk_intermediate_split_first(F_UVEC,step_num)
real DF_GUIJ11         = rk_intermediate_split_first(F_GUIJ11,step_num)
real DF_GUIJ12         = rk_intermediate_split_first(F_GUIJ12,step_num)
real DF_GUIJ13         = rk_intermediate_split_first(F_GUIJ13,step_num)
real DF_GUIJ21         = rk_intermediate_split_first(F_GUIJ21,step_num) 
real DF_GUIJ22         = rk_intermediate_split_first(F_GUIJ22,step_num) 
real DF_GUIJ23         = rk_intermediate_split_first(F_GUIJ23,step_num) 
real DF_GUIJ31         = rk_intermediate_split_first(F_GUIJ31,step_num) 
real DF_GUIJ32         = rk_intermediate_split_first(F_GUIJ32,step_num) 
real DF_GUIJ33         = rk_intermediate_split_first(F_GUIJ33,step_num) 
real DF_CS             = rk_intermediate_split_first(F_CS,step_num)
real DF_VISC_HEAT      = rk_intermediate_split_first(F_VISC_HEAT,step_num)
real3 DF_VISC_FORCVEC  = rk_intermediate_split_first(F_VISC_FORCVEC,step_num)
real3 DF_BVEC          = rk_intermediate_split_first(F_BVEC,step_num)
real3 DF_JVEC          = rk_intermediate_split_first(F_JVEC,step_num)
real3 DF_AVEC          = rk_intermediate_split_first(F_AVEC,step_num)
real3 DF_JXBVEC        = rk_intermediate_split_first(F_JXBVEC,step_num)
real3 DF__ADV_DERVEC   = rk_intermediate_split_first(F__ADV_DERVEC,step_num)
real  DF_UU_SPHR       = rk_intermediate_split_first(F_UU_SPHR,step_num)
real  DF_UU_SPHT       = rk_intermediate_split_first(F_U_SPHT,step_num)
real  DF_UU_SPHP       = rk_intermediate_split_first(F_UU_SPHP,step_num)
real  DF_NRHO          = rk_intermediate_split_first(F_NRHO,step_num)
real  DF_SS            = rk_intermediate_split_first(F_SS,step_num)
real  DF_TT            = rk_intermediate_split_first(F_TT,step_num)
real  DF_ETH           = rk_intermediate_split_first(F_ETH,step_num)
real  DF_RUN_AVER      = rk_intermediate_split_first(F_RUN_AVER,step_num)
real DF_PHIUU         = rk_intermediate_split_first(F_PHIUU,step_num)
real DF_SS_RUN_AVER    = rk_intermediate_split_first(F_SS_RUN_AVER,step_num)
real DF_RHO            = rk_intermediate_split_first(F_RHO,step_num)
real DF_RHON            = rk_intermediate_split_first(F_RHON,step_num)
real3 DF_DUST_VELOCITY[ndustspec]
real  DF_DUST_DENSITY[ndustspec]
real  DF_DUST_MASS[ndustspec]
real  DF_DUST_ICE_MASS[ndustspec]
if(ldustvelocity)
{
	for dustspec in 0:ndustspec
	{
		DF_DUST_VELOCITY[dustspec] = rk_intermediate_split_first(F_DUST_VELOCITY[dustspec],step_num)
	}
}
if(ldustdensity)
{
	for dustspec in 0:ndustspec
	{
		DF_DUST_DENSITY[dustspec]  = rk_intermediate_split_first(F_DUST_DENSITY[dustspec],step_num) 
		DF_DUST_MASS[dustspec]     = rk_intermediate_split_first(F_DUST_MASS[dustspec],step_num) 
		DF_DUST_ICE_MASS[dustspec] = rk_intermediate_split_first(F_DUST_ICE_MASS[dustspec],step_num) 
	}
}
real3 DF_UNVEC = rk_intermediate_split_first(F_UNVEC,step_num)
real DF_CHEMISTRY_SPECIES[nchemspec]
if(lchemistry)
{
	for i in 0:nchemspec
	{
		DF_CHEMISTRY_SPECIES[i] = rk_intermediate_split_first(F_CHEMISTRY_SPECIES[i],step_num)
	}
}
real DF_CHEMISTRY_REACTIONS[nchemspec]
real DF_ECR = rk_intermediate_split_first(F_ECR,step_num)
real DF_XX_CHIRAL = rk_intermediate_split_first(F_XX_CHIRAL,step_num)
real DF_YY_CHIRAL = rk_intermediate_split_first(F_YY_CHIRAL,step_num)
real DF_ZZ_CHIRAL = rk_intermediate_split_first(F_ZZ_CHIRAL,step_num)

DF_SIGMA = rk_intermediate_split_first(F_SIGMA,step_num)
#define DF_LNTT DF_TT
#define DF_LNRHO DF_RHO
#define DF_LNRHON DF_RHON
#define DF_UU DF_UVEC
#define DF_UUN DF_UNVEC


#define DF_UX DF_UU.x
#define DF_UY DF_UU.y
#define DF_UZ DF_UU.z

#define DF_UNX DF_UUN.x
#define DF_UNY DF_UUN.y
#define DF_UNZ DF_UUN.z

#define DF_AX DF_AA.x
#define DF_AY DF_AA.y
#define DF_AZ DF_AA.z
#define DF_AA DF_AVEC
#define AC_0 (0)
real DF_STRESS_0 = 0.0
real DF_STRESS_1 = 0.0
real DF_STRESS_2 = 0.0
real DF_STRESS_3 = 0.0
real DF_STRESS_4 = 0.0
real DF_STRESS_5 = 0.0


#if LAXIONSU2BACK
real DF_AXI_PSI    = rk_intermediate_split_first(F_AXI_PSI,step_num)
real DF_AXI_PSIDOT = rk_intermediate_split_first(F_AXI_PSIDOT,step_num)
real DF_AXI_IMPSI  = rk_intermediate_split_first(F_AXI_IMPSI,step_num)
real DF_AXI_IMPSIDOT  = rk_intermediate_split_first(F_AXI_IMPSIDOT,step_num)

real DF_AXI_TR    = rk_intermediate_split_first(F_AXI_TR,step_num)
real DF_AXI_TRDOT = rk_intermediate_split_first(F_AXI_TRDOT,step_num)
real DF_AXI_IMTR    = rk_intermediate_split_first(F_AXI_IMTR,step_num)
real DF_AXI_IMTRDOT = rk_intermediate_split_first(F_AXI_IMTRDOT,step_num)

real DF_AXI_TL    = rk_intermediate_split_first(F_AXI_TL,step_num)
real DF_AXI_TLDOT = rk_intermediate_split_first(F_AXI_TLDOT,step_num)
real DF_AXI_IMTL    = rk_intermediate_split_first(F_AXI_IMTL,step_num)
real DF_AXI_IMTLDOT = rk_intermediate_split_first(F_AXI_IMTLDOT,step_num)

real DF_AXI_PSIL    = rk_intermediate_split_first(F_AXI_PSIL,step_num)
real DF_AXI_PSILDOT = rk_intermediate_split_first(F_AXI_PSILDOT,step_num)
real DF_AXI_IMPSIL    = rk_intermediate_split_first(F_AXI_IMPSIL,step_num)
real DF_AXI_IMPSILDOT = rk_intermediate_split_first(F_AXI_IMPSILDOT,step_num)
#endif


real DF_GAMMA = rk_intermediate_split_first(F_GAMMA,step_num)
real DF_A0    = rk_intermediate_split_first(F_A0,step_num)
real DF_RHOE  = rk_intermediate_split_first(F_RHOE,step_num)
real3 DF_EVEC          = rk_intermediate_split_first(F_EVEC,step_num)
real  DF_DIVA_NAME     = rk_intermediate_split_first(F_DIVA_NAME,step_num)
#define DF_EX DF_EVEC.x
#define DF_EY DF_EVEC.y
#define DF_EZ DF_EVEC.z

#if LBACKREACT_INFL
real DF_INFL_PHI  = rk_intermediate_split_first(F_INFL_PHI,step_num)
real DF_INFL_DPHI = rk_intermediate_split_first(F_INFL_DPHI,step_num)
#endif

#if LKLEIN_GORDON
//For klein gordon
real DF_PHI  = rk_intermediate_split_first(F_PHI,step_num)
real DF_DPHI = rk_intermediate_split_first(F_DPHI,step_num)
real DF_PSI  = rk_intermediate_split_first(F_PSI,step_num)
real DF_DPSI  = rk_intermediate_split_first(F_DPSI,step_num)

#define DF_PHI_UP_RE DF_PHI
real DF_PHI_UP_IM = rk_intermediate_split_first(F_PHI_UP_IM,step_num)
real DF_PHI_DOWN_RE = rk_intermediate_split_first(F_PHI_DOWN_RE,step_num)
real DF_PHI_DOWN_IM = rk_intermediate_split_first(F_PHI_DOWN_IM,step_num)

#define DF_DPHI_UP_RE DF_DPHI
real DF_DPHI_UP_IM = rk_intermediate_split_first(F_DPHI_UP_IM,step_num)
real DF_DPHI_DOWN_RE = rk_intermediate_split_first(F_DPHI_DOWN_RE,step_num)
real DF_DPHI_DOWN_IM = rk_intermediate_split_first(F_DPHI_DOWN_IM,step_num)
#endif

real3 DF_AXTESTVEC = rk_intermediate_split_first(F_AXTESTVEC,step_num)

real DF_MU5 = rk_intermediate_split_first(F_MU5,step_num)
real DF_MUS = rk_intermediate_split_first(F_MUS,step_num)

real3 DF_UXVEC = rk_intermediate_split_first(F_UXVEC,step_num)
real3 DF_UXSVEC = rk_intermediate_split_first(F_UXSVEC,step_num)
real3 DF_OXVEC = rk_intermediate_split_first(F_OXVEC,step_num)
real3 DF_OXSVEC = rk_intermediate_split_first(F_OXSVEC,step_num)

real3 DF_BXVEC = rk_intermediate_split_first(F_BXVEC,step_num)
real3 DF_JXVEC = rk_intermediate_split_first(F_JXVEC,step_num)

#define DF_AXTESVEC DF_AXTESTVEC

real DF_IPOLY__MOD__CDATA[6]
if(lpolymer)
{
	for i in 0:6 {DF_IPOLY__MOD__CDATA[i] = rk_intermediate_split_first(F_POLY[i],step_num)}
}

dt1_max__mod__cdata = 0.0

const int AC_itsub__mod__cdata  = step_num+1
#if LDUSTDENSITY
reac_dust__mod__cdata = 0.0
#endif
