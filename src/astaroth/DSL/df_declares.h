real3 DF_UVEC          = read_df(F_UVEC,step_num,AC_lsplit_update__mod__cdata)
real DF_GUIJ11         = read_df(F_GUIJ11,step_num,AC_lsplit_update__mod__cdata)
real DF_GUIJ12         = read_df(F_GUIJ12,step_num,AC_lsplit_update__mod__cdata)
real DF_GUIJ13         = read_df(F_GUIJ13,step_num,AC_lsplit_update__mod__cdata)
real DF_GUIJ21         = read_df(F_GUIJ21,step_num,AC_lsplit_update__mod__cdata) 
real DF_GUIJ22         = read_df(F_GUIJ22,step_num,AC_lsplit_update__mod__cdata) 
real DF_GUIJ23         = read_df(F_GUIJ23,step_num,AC_lsplit_update__mod__cdata) 
real DF_GUIJ31         = read_df(F_GUIJ31,step_num,AC_lsplit_update__mod__cdata) 
real DF_GUIJ32         = read_df(F_GUIJ32,step_num,AC_lsplit_update__mod__cdata) 
real DF_GUIJ33         = read_df(F_GUIJ33,step_num,AC_lsplit_update__mod__cdata) 
real DF_CS             = read_df(F_CS,step_num,AC_lsplit_update__mod__cdata)
real DF_VISC_HEAT      = read_df(F_VISC_HEAT,step_num,AC_lsplit_update__mod__cdata)
real3 DF_VISC_FORCVEC  = read_df(F_VISC_FORCVEC,step_num,AC_lsplit_update__mod__cdata)
real3 DF_BVEC          = read_df(F_BVEC,step_num,AC_lsplit_update__mod__cdata)
real3 DF_JVEC          = read_df(F_JVEC,step_num,AC_lsplit_update__mod__cdata)
real3 DF_AVEC          = read_df(F_AVEC,step_num,AC_lsplit_update__mod__cdata)
real3 DF_JXBVEC        = read_df(F_JXBVEC,step_num,AC_lsplit_update__mod__cdata)
real3 DF__ADV_DERVEC   = read_df(F__ADV_DERVEC,step_num,AC_lsplit_update__mod__cdata)
real  DF_UU_SPHR       = read_df(F_UU_SPHR,step_num,AC_lsplit_update__mod__cdata)
real  DF_UU_SPHT       = read_df(F_U_SPHT,step_num,AC_lsplit_update__mod__cdata)
real  DF_UU_SPHP       = read_df(F_UU_SPHP,step_num,AC_lsplit_update__mod__cdata)
real  DF_NRHO          = read_df(F_NRHO,step_num,AC_lsplit_update__mod__cdata)
real  DF_SS            = read_df(F_SS,step_num,AC_lsplit_update__mod__cdata)
real  DF_TT            = read_df(F_TT,step_num,AC_lsplit_update__mod__cdata)
real  DF_ETH           = read_df(F_ETH,step_num,AC_lsplit_update__mod__cdata)
real  DF_RUN_AVER      = read_df(F_RUN_AVER,step_num,AC_lsplit_update__mod__cdata)
real DF_PHIUU          = read_df(F_PHIUU,step_num,AC_lsplit_update__mod__cdata)
real DF_SS_RUN_AVER    = read_df(F_SS_RUN_AVER,step_num,AC_lsplit_update__mod__cdata)
real DF_RHO            = read_df(F_RHO,step_num,AC_lsplit_update__mod__cdata)
real DF_RHON           = read_df(F_RHON,step_num,AC_lsplit_update__mod__cdata)
real3 DF_DUST_VELOCITY[ndustspec]
real  DF_DUST_DENSITY[ndustspec]
real  DF_DUST_MASS[ndustspec]
real  DF_DUST_ICE_MASS[ndustspec]

real DF_LAMRA = read_df(F_LAMRA,step_num,AC_lsplit_update__mod__cdata)
if(ldustvelocity)
{
	for dustspec in 0:ndustspec
	{
		DF_DUST_VELOCITY[dustspec] = read_df(F_DUST_VELOCITY[dustspec],step_num,AC_lsplit_update__mod__cdata)
	}
}
real DF_ICC__MOD__CDATA[npscalar]
if(lpscalar)
{
	for pscalar in 0:npscalar
	{
		DF_ICC__MOD__CDATA[pscalar] = read_df(F_CVEC[pscalar],step_num,AC_lsplit_update__mod__cdata)
	}
}
if(ldustdensity)
{
	for dustspec in 0:ndustspec
	{
		DF_DUST_DENSITY[dustspec]  = read_df(F_DUST_DENSITY[dustspec],step_num,AC_lsplit_update__mod__cdata) 
		DF_DUST_MASS[dustspec]     = read_df(F_DUST_MASS[dustspec],step_num,AC_lsplit_update__mod__cdata) 
		DF_DUST_ICE_MASS[dustspec] = read_df(F_DUST_ICE_MASS[dustspec],step_num,AC_lsplit_update__mod__cdata) 
	}
}
real3 DF_UNVEC = read_df(F_UNVEC,step_num,AC_lsplit_update__mod__cdata)
real DF_CHEMISTRY_SPECIES[nchemspec]
if(lchemistry)
{
	for i in 0:nchemspec
	{
		DF_CHEMISTRY_SPECIES[i] = read_df(F_CHEMISTRY_SPECIES[i],step_num,AC_lsplit_update__mod__cdata)
	}
}
real DF_CHEMISTRY_REACTIONS[nchemspec]
real DF_ECR = read_df(F_ECR,step_num,AC_lsplit_update__mod__cdata)
real DF_XX_CHIRAL = read_df(F_XX_CHIRAL,step_num,AC_lsplit_update__mod__cdata)
real DF_YY_CHIRAL = read_df(F_YY_CHIRAL,step_num,AC_lsplit_update__mod__cdata)
real DF_ZZ_CHIRAL = read_df(F_ZZ_CHIRAL,step_num,AC_lsplit_update__mod__cdata)

DF_SIGMA = read_df(F_SIGMA,step_num,AC_lsplit_update__mod__cdata)
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
real DF_AXI_PSI    = read_df(F_AXI_PSI,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_PSIDOT = read_df(F_AXI_PSIDOT,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMPSI  = read_df(F_AXI_IMPSI,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMPSIDOT  = read_df(F_AXI_IMPSIDOT,step_num,AC_lsplit_update__mod__cdata)

real DF_AXI_TR    = read_df(F_AXI_TR,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_TRDOT = read_df(F_AXI_TRDOT,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMTR    = read_df(F_AXI_IMTR,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMTRDOT = read_df(F_AXI_IMTRDOT,step_num,AC_lsplit_update__mod__cdata)

real DF_AXI_TL    = read_df(F_AXI_TL,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_TLDOT = read_df(F_AXI_TLDOT,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMTL    = read_df(F_AXI_IMTL,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMTLDOT = read_df(F_AXI_IMTLDOT,step_num,AC_lsplit_update__mod__cdata)

real DF_AXI_PSIL      = read_df(F_AXI_PSIL,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_PSILDOT   = read_df(F_AXI_PSILDOT,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMPSIL    = read_df(F_AXI_IMPSIL,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMPSILDOT = read_df(F_AXI_IMPSILDOT,step_num,AC_lsplit_update__mod__cdata)

real DF_AXI_UR      = read_df(F_AXI_UR,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_URDOT   = read_df(F_AXI_URDOT,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMUR    = read_df(F_AXI_IMUR,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMURDOT = read_df(F_AXI_IMURDOT,step_num,AC_lsplit_update__mod__cdata)

real DF_AXI_UL      = read_df(F_AXI_UL,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_ULDOT   = read_df(F_AXI_ULDOT,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMUL    = read_df(F_AXI_IMUL,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMULDOT = read_df(F_AXI_IMULDOT,step_num,AC_lsplit_update__mod__cdata)
#endif

#if LAXIONU1BACK 
real DF_AXI_AR      = read_df(F_AXI_AR,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_ARDOT   = read_df(F_AXI_ARDOT,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMAR    = read_df(F_AXI_IMAR,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMARDOT = read_df(F_AXI_IMARDOT,step_num,AC_lsplit_update__mod__cdata)

real DF_AXI_AL      = read_df(F_AXI_AL,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_ALDOT   = read_df(F_AXI_ALDOT,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMAL    = read_df(F_AXI_IMAL,step_num,AC_lsplit_update__mod__cdata)
real DF_AXI_IMALDOT = read_df(F_AXI_IMALDOT,step_num,AC_lsplit_update__mod__cdata)
#endif



real DF_GAMMA = read_df(F_GAMMA,step_num,AC_lsplit_update__mod__cdata)
real DF_A0    = read_df(F_A0,step_num,AC_lsplit_update__mod__cdata)
real DF_RHOE  = read_df(F_RHOE,step_num,AC_lsplit_update__mod__cdata)
real3 DF_EVEC          = read_df(F_EVEC,step_num,AC_lsplit_update__mod__cdata)
real  DF_DIVA_NAME     = read_df(F_DIVA_NAME,step_num,AC_lsplit_update__mod__cdata)
#define DF_EX DF_EVEC.x
#define DF_EY DF_EVEC.y
#define DF_EZ DF_EVEC.z

#if LBACKREACT_INFL
real DF_INFL_PHI  = read_df(F_INFL_PHI,step_num,AC_lsplit_update__mod__cdata)
real DF_INFL_DPHI = read_df(F_INFL_DPHI,step_num,AC_lsplit_update__mod__cdata)
#endif

#if LKLEIN_GORDON
//For klein gordon
real DF_PHI  = read_df(F_PHI,step_num,AC_lsplit_update__mod__cdata)
real DF_DPHI = read_df(F_DPHI,step_num,AC_lsplit_update__mod__cdata)
real DF_PSI  = read_df(F_PSI,step_num,AC_lsplit_update__mod__cdata)
real DF_DPSI  = read_df(F_DPSI,step_num,AC_lsplit_update__mod__cdata)

#define DF_PHI_UP_RE DF_PHI
real DF_PHI_UP_IM = read_df(F_PHI_UP_IM,step_num,AC_lsplit_update__mod__cdata)
real DF_PHI_DOWN_RE = read_df(F_PHI_DOWN_RE,step_num,AC_lsplit_update__mod__cdata)
real DF_PHI_DOWN_IM = read_df(F_PHI_DOWN_IM,step_num,AC_lsplit_update__mod__cdata)

#define DF_DPHI_UP_RE DF_DPHI
real DF_DPHI_UP_IM = read_df(F_DPHI_UP_IM,step_num,AC_lsplit_update__mod__cdata)
real DF_DPHI_DOWN_RE = read_df(F_DPHI_DOWN_RE,step_num,AC_lsplit_update__mod__cdata)
real DF_DPHI_DOWN_IM = read_df(F_DPHI_DOWN_IM,step_num,AC_lsplit_update__mod__cdata)
#endif

real3 DF_AXTESTVEC = read_df(F_AXTESTVEC,step_num,AC_lsplit_update__mod__cdata)

real DF_MU5 = read_df(F_MU5,step_num,AC_lsplit_update__mod__cdata)
real DF_MUS = read_df(F_MUS,step_num,AC_lsplit_update__mod__cdata)

real3 DF_UXVEC = read_df(F_UXVEC,step_num,AC_lsplit_update__mod__cdata)
real3 DF_UXSVEC = read_df(F_UXSVEC,step_num,AC_lsplit_update__mod__cdata)
real3 DF_OXVEC = read_df(F_OXVEC,step_num,AC_lsplit_update__mod__cdata)
real3 DF_OXSVEC = read_df(F_OXSVEC,step_num,AC_lsplit_update__mod__cdata)

real3 DF_BXVEC = read_df(F_BXVEC,step_num,AC_lsplit_update__mod__cdata)
real3 DF_JXVEC = read_df(F_JXVEC,step_num,AC_lsplit_update__mod__cdata)

#define DF_AXTESVEC DF_AXTESTVEC

real DF_IPOLY__MOD__CDATA[6]
if(lpolymer)
{
	for i in 0:6 {DF_IPOLY__MOD__CDATA[i] = read_df(F_POLY[i],step_num,AC_lsplit_update__mod__cdata)}
}

dt1_max__mod__cdata = 0.0

const int AC_itsub__mod__cdata  = step_num+1
#if LDUSTDENSITY
reac_dust__mod__cdata = 0.0
#endif

const bool AC_lsubstepping_in_time__mod__cdata = false
