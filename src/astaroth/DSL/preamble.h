#include "../../PC_moduleflags.h"

#if LENTROPY
#define LENERGY 1
#endif
const bool AC_lpencil_check__mod__cdata = false
const bool AC_lpencil_check_at_work__mod__cdata = false

#include "../typedefs.h"

#define AC_n1 n1
#define AC_m1 m1
#define AC_l1 l1

#define nhless AC_nhless__mod__hydro
#define n_odevars AC_n_odevars__mod__cdata

#include "PC_nghost.h"

#define NGHOST_VAL NGHOST
const int prof_nz = 150
#define mreactions AC_mreactions__mod__chemistry
#define nreactions AC_nreactions__mod__chemistry

#define LNRHO RHO
#define ks_modes AC_ks_modes__mod__forcing

#define AC_pretend_lnTT AC_pretend_lntt__mod__cdata

#define hcond_prof_size AC_hcond_prof_size__mod__energy
#define chit_prof_size AC_chit_prof_size__mod__energy
#define chit_prof_fluct_stored_size AC_chit_prof_fluct_stored_size__mod__energy

#define AC_l2 AC_l2__mod__cdata
#define AC_n2 AC_n2__mod__cdata
#define AC_m2 AC_m2__mod__cdata

#define AC_dx2_bound AC_dx2_bound__mod__cdata
#define AC_dy2_bound AC_dy2_bound__mod__cdata
#define AC_dz2_bound AC_dz2_bound__mod__cdata


#define AC_maux__mod__cparam maux__mod__cparam
#define iphF_UU F_PHIUU

#define AC_n2__mod__cparam AC_nx+NGHOST_VAL+1 
#define AC_m2__mod__cparam AC_ny+NGHOST_VAL+1 
#define AC_l2__mod__cparam AC_nz+NGHOST_VAL+1 

#define AC_NGHOST_VAL__mod__cparam NGHOST_VAL

#define AC_mx mx 
#define AC_my my 
#define AC_mz mz 

#define AC_nx nx 
#define AC_ny ny
#define AC_nz nz 

#define AC_nxgrid nxgrid 
#define AC_nygrid nygrid 
#define AC_nzgrid nzgrid 

#define AC_dsx AC_ds.x
#define AC_dsy AC_ds.y
#define AC_dsz AC_ds.z

#define n__mod__cdata ((NGHOST-AC_nmin.z)+vertexIdx.z+1)
#define m__mod__cdata ((NGHOST-AC_nmin.y)+vertexIdx.y+1)
#define AC_n__mod__cdata n__mod__cdata
#define AC_m__mod__cdata m__mod__cdata

#include "../stdlib/math"
#include "../stdlib/cg.h"
#include "../stdlib/general_grid"
#include "../stdlib/general_derivs.h"
//#include "../stdlib/pc_derivs.h"
#include "../stdlib/general_operators.h"
#define AC_NGHOST__mod__cparam nghost
#define REAL_MAX AC_REAL_MAX
//TP: nphis1 and nphis2 don't actually work. simply declared to compile the code
//
#define IN_DSL 1
#include "cparam.h"
#undef AC_bot
#undef AC_top
#include "../../../cparam_c.h"
const int AC_xtop__mod__equationofstate=nx
#include "../../../cparam_pencils.inc_c.h"
#include "PC_modulepardecs.h"
#include "../stdlib/optimized_integrators.h"
#include "../stdlib/slope_limited_diffusion.h"

#include "../fieldecs.h"

output real AC_maxchi
#if Ltimestep_rkf_lowsto_MODULE
enum PC_SUB_STEP_NUMBER
{
	PC_FIRST_SUB_STEP,
	PC_SECOND_SUB_STEP,
	PC_THIRD_SUB_STEP,
	PC_FOURTH_SUB_STEP,
	PC_FIFTH_SUB_STEP
}
#else
enum PC_SUB_STEP_NUMBER
{
	PC_FIRST_SUB_STEP,
	PC_SECOND_SUB_STEP,
	PC_THIRD_SUB_STEP
}
#endif

#define AC_lfirst__mod__cdata (step_num == 0)
#define AC_llast__mod__cdata (step_num == AC_num_substeps__mod__cdata-1)
maxval(real x) {return x}
maxval(real[] arr) 
{
	real res = -AC_REAL_MAX
	for i in 0:size(arr)
	{
		res = max(res,arr[i])
	}
	return res
}
minval(real x) {return x}

#define notanumber(x) (false)

output real AC_dt1_max
output real AC_dt1_advec
output real AC_dt1_diffus

#if Leos_idealgas_MODULE
#define AC_lntt0__mod__equationofstate AC_lnTT0__mod__equationofstate
#endif

#define AC_gamma1__mod__energy  AC_gamma1__mod__equationofstate 
#define AC_gamma1__mod__magnetc AC_gamma1__mod__equationofstate 
#define AC_lupdate_courant_dt__mod__cdata (AC_lfirst__mod__cdata && lcourant_dt)
#define lproc_print__mod__cdata (false)
#define AC_allproc_print__mod__cdata (false)
#define lfirstpoint__mod__cdata (false)
#define AC_headt__mod__cdata (false)


output real uux_initial_max
output real uuy_initial_max
output real uuz_initial_max
       
output real aax_initial_max
output real aay_initial_max
output real aaz_initial_max

output real rho_initial_max
output real ss_initial_max

output real AC_maxnu

ac_unused_real_array_1d(index) {suppress_unused_warning(index); return 0.0}
ac_unused_real_array_2d(index1,index2) {suppress_unused_warning(index1); suppress_unused_warning(index2); return 0.0}
ac_unused_real_array_3d(index1,index2,index3) {suppress_unused_warning(index1); suppress_unused_warning(index2); suppress_unused_warning(index3); return 0.0}
const real ac_unused_real_scalar = 0.0
const real ac_real_unused_scalar = 0.0
#define epsi AC_REAL_EPSILON
#define AC_nxgrid__mod__cparam nxgrid
#define AC_nygrid__mod__cparam nxgrid
#define AC_nzgrid__mod__cparam nxgrid
#define maux__mod__cparam maux

global output real AC_Arms
#include "../forcing"
#define ikx (vertexIdx.x-NGHOST+1)
#define iky (vertexIdx.y-NGHOST+1)
#define ikz (vertexIdx.z-NGHOST+1)
#define AC_ipx__mod__cdata (AC_domain_coordinates.x)
#define AC_ipy__mod__cdata (AC_domain_coordinates.y)
#define AC_ipz__mod__cdata (AC_domain_coordinates.z)
const int AC_0 = 0
#define AC_impossible__mod__cparam impossible
#define AC_tini__mod__cparam tini

const int AC_ij_table__mod__gravitational_waves_htxk = [
					[1,4,6],
					[4,2,5],
					[6,5,3]
				       ]
gmem real AC_nphis1__mod__cdata[AC_mlocal.y]
const real AC_arms__mod__magnetic = 0.0

#if LDUSTVELOCITY
#else
gmem real AC_reac_dust__mod__cdata[1]
#endif

#if LCHEMISTRY
#else
gmem real AC_reac_chem__mod__cdata[1]
#endif

const real AC_ascale__mod__cdata = 0.0

#define AC_ftopktop__mod__energy AC_FtopKtop__mod__energy
#define AC_fbotkbot__mod__energy AC_FbotKbot__mod__energy
#define AC_fbot__mod__energy AC_Fbot__mod__energy 
#define AC_ftop__mod__energy AC_Ftop__mod__energy 
#define FbotKbot AC_FbotKbot__mod__energy
#define Fbot AC_Fbot__mod__energy
#define fbot Fbot
#define ftop Ftop
#define FtopKtop AC_FtopKtop__mod__energy
#define Ftop AC_Ftop__mod__energy
#define fbcx AC_fbcx__mod__cdata
#define fbcy AC_fbcy__mod__cdata
#define fbcz AC_fbcz__mod__cdata
#define AC_ldensity_nolog AC_ldensity_nolog__mod__cdata
#define AC_cs20 AC_cs20__mod__equationofstate
#define AC_gamma_m1 AC_gamma_m1__mod__equationofstate
#define AC_gamma AC_gamma__mod__equationofstate
#define AC_lnrho0 AC_lnrho0__mod__equationofstate
#define AC_cv1 AC_cv1__mod__equationofstate
#define AC_lheatc_chiconst AC_lheatc_chiconst__mod__energy
#define AC_chi AC_chi__mod__energy
#define AC_lheatc_kramers AC_lheatc_kramers__mod__energy
#define nkramers AC_nkramers__mod__energy
#define AC_cp AC_cp__mod__equationofstate
#define AC_cv AC_cv__mod__equationofstate
#define hcond0_kramers  AC_hcond0_kramers__mod__energy
#define AC_pretend_lntt  AC_pretend_lntt__mod__cdata
#define AC_cs2bot AC_cs2bot__mod__equationofstate
#define cs2top AC_cs2top__mod__equationofstate
#define AC_lreference_state AC_lreference_state__mod__cdata
#define ltemperature_nolog AC_ltemperature_nolog__mod__cdata
#define AC_lread_oldsnap AC_lread_oldsnap__mod__cdata

#include "../stdlib/bc.h"
#include "../bcs/funcs_full.h"
#include "../bcs/funcs_overload.h"
#include "../hydro/before_boundary.h"

#ifndef LDISP_CURRENT
const bool AC_lphi_hom__mod__disp_current = false
const bool AC_lphi_linear_regime__mod__disp_current = false
const bool AC_lnoncollinear_eb__mod__disp_current = false
const bool AC_lnoncollinear_eb_aver__mod__disp_current = false
const bool AC_lcollinear_eb__mod__disp_current = false
const bool AC_lcollinear_eb_aver__mod__disp_current = false
const bool AC_lallow_bprime_zero__mod__disp_current = false
#endif

#define AC_iproc_world__mod__cdata (0)
#include "../axionSU2back.h"


#include "../backreact_infl.h"

#if LPOLYMER
#else
const real AC_trelax_poly__mod__cdata = 0.0
#endif
#if Leos_ionization_MODULE
const real yhmin = AC_REAL_MIN
const real yhmax = 1-AC_REAL_EPSILON
#endif

//Crucially lmultithread has to be false from the point of view GPU that we don't do some things twice
const bool AC_lmultithread__mod__cdata = false
#if Lentropy_MODULE
Profile<Z> AC_ssmz__mod__energy
Profile<Z> AC_del2ssmz__mod__energy
#endif

#if LDUSTVELOCITY
const real ueta = 0.0
const real teta = 0.0
const real ul0  = 0.0
const real tl0  = 0.0
const real tl01  = tl0/(tl0+tini)
const real teta1 = teta/(teta+tini)
#endif
