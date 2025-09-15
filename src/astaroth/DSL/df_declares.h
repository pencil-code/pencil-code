real3 DF_UVEC          = real3(0.0,0.0,0.0)
real DF_GUIJ11         = 0.0 
real DF_GUIJ12         = 0.0 
real DF_GUIJ13         = 0.0 
real DF_GUIJ21         = 0.0 
real DF_GUIJ22         = 0.0 
real DF_GUIJ23         = 0.0 
real DF_GUIJ31         = 0.0 
real DF_GUIJ32         = 0.0 
real DF_GUIJ33         = 0.0 
real DF_CS             = 0.0
real DF_VISC_HEAT      = 0.0
real3 DF_VISC_FORCVEC  = real3(0.0,0.0,0.0)
real3 DF_DIVA          = real3(0.0,0.0,0.0)
real3 DF_BVEC          = real3(0.0,0.0,0.0)
real3 DF_JVEC          = real3(0.0,0.0,0.0)
real3 DF_AVEC          = real3(0.0,0.0,0.0)
real3 DF_EVEC          = real3(0.0,0.0,0.0)
real3 DF_JXBVEC        = real3(0.0,0.0,0.0)
real3 DF_JXBVEC        = real3(0.0,0.0,0.0)
real3 DF__ADV_DERVEC   = real3(0.0,0.0,0.0)
real  DF_UU_SPHR       = 0.0
real  DF_UU_SPHT       = 0.0
real  DF_UU_SPHP       = 0.0
real  DF_NRHO          = 0.0
real  DF_SS            = 0.0
real  DF_TT            = 0.0
real  DF_ETH           = 0.0
real  DF_RUN_AVER      = 0.0
real DF_PHIUU         = 0.0
real DF_SS_RUN_AVER    = 0.0
real DF_RHO            = 0.0
real DF_RHON            = 0.0
real3 DF_DUST_VELOCITY[AC_ndustspec__mod__cparam]
real  DF_DUST_DENSITY[AC_ndustspec__mod__cparam]
real  DF_DUST_MASS[AC_ndustspec__mod__cparam]
real  DF_DUST_ICE_MASS[AC_ndustspec__mod__cparam]
real3 DF_UNVEC = real3(0.0,0.0,0.0)
real DF_CHEMISTRY_SPECIES[nchemspec]
if(lchemistry)
{
	for i in 0:nchemspec
	{
		DF_CHEMISTRY_SPECIES[i] = 0.0
	}
}
real DF_CHEMISTRY_REACTIONS[nchemspec]
real DF_ECR = 0.0
real DF_XX_CHIRAL = 0.0
real DF_YY_CHIRAL = 0.0
real DF_ZZ_CHIRAL = 0.0
if(ldustvelocity)
{
	for dustspec in 0:AC_ndustspec__mod__cparam
	{
		DF_DUST_VELOCITY[dustspec] = real3(0.0,0.0,0.0)
	}
}
if(ldustdensity)
{
	for dustspec in 0:AC_ndustspec__mod__cparam
	{
		DF_DUST_DENSITY[dustspec]  = 0.0
		DF_DUST_MASS[dustspec]     = 0.0
		DF_DUST_ICE_MASS[dustspec] = 0.0
	}
}

DF_SIGMA = 0.0
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
dt1_max__mod__cdata = 0.0
