Field UUX, UUY, UUZ
#define F_UX UUX
#define F_UY UUY
#define F_UZ UUZ

Field RHO
#define F_RHO  RHO
#define F_LNRHO F_RHO

Field AAX, AAY, AAZ
#define F_AX AAX
#define F_AY AAY
#define F_AZ AAZ

Field3 F_DUST_VELOCITY[AC_ndustspec__mod__cparam]
Field  F_DUST_DENSITY[AC_ndustspec__mod__cparam]
Field  F_DUST_MASS[AC_ndustspec__mod__cparam]
Field  F_DUST_ICE_MASS[AC_ndustspec__mod__cparam]

Field F_PHIUU
Field F_U0X, F_U0Y, F_U0Z
Field F_LORENTZ
Field F_GUIJ11,F_GUIJ21,F_GUIJ31,F_GUIJ12,F_GUIJ22,F_GUIJ32,F_GUIJ13,F_GUIJ23,F_GUIJ33
Field F_OX, F_OY, F_OZ
Field F_UU_SPHX, F_UU_SPHY, F_UU_SPHZ
Field F_BB_SPHX, F_BB_SPHY, F_BB_SPHZ
Field F_HLESS
Field F_EOSVAR2
Field F_GLOBAL_CS2
Field F_PP
Field F_SS
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

Field F_JX,F_JY,F_JZ
Field F_EDOTX,F_EDOTY,F_EDOTZ
Field F_LAM
Field F_TT
Field F_GLOBAL_HCOND
Field F_SS_RUN_AVER
Field F_ADV_DERX
Field F_ADV_DERY
Field F_ADV_DERZ
Field F_GLOBAL_GLNTX,F_GLOBAL_GLNTY,F_GLOBAL_GLNTZ
Field F_GLOBAL_GLHX,F_GLOBAL_GLHY,F_GLOBAL_GLHZ
Field F_HYPREX, F_HYPREY, F_HYPREZ
Field F_ETAT



const Field3 F_GLOBAL_GLNTVEC = {F_GLOBAL_GLNTX,F_GLOBAL_GLNTY,F_GLOBAL_GLNTZ}
const Field3 F_AVEC    = {F_AX, F_AY, F_AZ}
#define F_AA F_AVEC
const Field3 F_UVEC    = {F_UX,F_UY,F_UZ}
const Field3 F_UU      = {F_UX,F_UY,F_UZ}
#define F_UU F_UVEC
const Field3 F_U0VEC   = {F_U0X, F_U0Y, F_U0Z}
const Field3 F_OVEC    = {F_OX, F_OY, F_OZ}
const Field3 F_UU_SPH_VEC  = {F_UU_SPHX, F_UU_SPHY, F_UU_SPHZ}
const Field3 F_UU_SPHVEC   = {F_UU_SPHX, F_UU_SPHY, F_UU_SPHZ}
const Field3 F_BB_SPHVEC   = {F_BB_SPHX, F_BB_SPHY, F_BB_SPHZ}
const Field3 F_BVEC        = {F_BX,F_BY,F_BZ}
const Field3 F_GLOBAL_GLHVEC = {F_GLOBAL_GLHX,F_GLOBAL_GLHY,F_GLOBAL_GLHZ}


const Field3 F_JVEC            = {F_JX,F_JY,F_JZ}
const Field3 F_EDOTVEC         = {F_EDOTX,F_EDOTY,F_EDOTZ}
const Field3 F__ADV_DERVEC     = {F_ADV_DERX,F_ADV_DERY,F_ADV_DERZ}
const Field3 F_HYPREVEC        = {F_HYPREX, F_HYPREY, F_HYPREZ}
const Field3 F_GLOBAL_EXT_AVEC = {F_GLOBAL_EXT_AX, F_GLOBAL_EXT_AY, F_GLOBAL_EXT_AZ}


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


run_const int AC_nxyz_max


gmem real AC_hcond_prof__mod__energy[AC_nxyz_max]
gmem real AC_dlnhcond_prof__mod__energy[AC_nxyz_max]
gmem real AC_chit_prof_stored__mod__energy[AC_nxyz_max]
gmem real AC_dchit_prof_stored__mod__energy[AC_nxyz_max]
gmem real AC_chit_prof_fluct_stored__mod__energy[AC_nxyz_max]
gmem real AC_dchit_prof_fluct_stored__mod__energy[AC_nxyz_max]

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

#define n__mod__cdata (vertexIdx.z+1)
#define m__mod__cdata (vertexIdx.y+1)
#define AC_n__mod__cdata (vertexIdx.z+1)
#define AC_m__mod__cdata (vertexIdx.y+1)

output real AC_df_rho_sum
Field TEST_1, TEST_2,TEST_3,TEST_4
output real AC_maxadvec
output real AC_maximum_error
const real tini = AC_REAL_MIN*5

#include "../stdlib/map.h"
