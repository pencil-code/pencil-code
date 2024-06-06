// automatically generated; do not edit!

int AC_lpressuregradient_gas // from hydro
#define lpressuregradient_gas AC_lpressuregradient_gas

int AC_ldiff_shock // from density
#define ldiff_shock AC_ldiff_shock
real AC_diffrho_shock // from density
#define diffrho_shock AC_diffrho_shock

real AC_eta // from magnetic
#define eta AC_eta
real AC_eta_hyper3 // from magnetic
#define eta_hyper3 AC_eta_hyper3
int AC_lresi_eta_const // from magnetic
#define lresi_eta_const AC_lresi_eta_const
int AC_lresi_hyper3 // from magnetic
#define lresi_hyper3 AC_lresi_hyper3

real AC_chi // from energy
#define chi AC_chi
real AC_nkramers // from energy
#define nkramers AC_nkramers
real AC_hcond0_kramers // from energy
#define hcond0_kramers AC_hcond0_kramers
real AC_hcond_Kconst // from energy
#define hcond_Kconst AC_hcond_Kconst
real AC_chi_hyper3 // from energy
#define chi_hyper3 AC_chi_hyper3
real AC_chi_t0 // from energy
#define chi_t0 AC_chi_t0
gmem real AC_hcond_prof[AC_nz] // from energy
#define hcond_prof AC_hcond_prof
gmem real AC_dlnhcond_prof[AC_nz] // from energy
#define dlnhcond_prof AC_dlnhcond_prof
gmem real AC_chit_prof_stored[AC_nz] // from energy
#define chit_prof_stored AC_chit_prof_stored
gmem real AC_dchit_prof_stored[AC_nz] // from energy
#define dchit_prof_stored AC_dchit_prof_stored
int AC_lheatc_hyper3ss // from energy
#define lheatc_hyper3ss AC_lheatc_hyper3ss
int AC_lheatc_shock // from energy
#define lheatc_shock AC_lheatc_shock

real AC_cs20 // from equationofstate
#define cs20 AC_cs20
real AC_gamma // from equationofstate
#define gamma AC_gamma
real AC_cv // from equationofstate
#define cv AC_cv
real AC_cp // from equationofstate
#define cp AC_cp
real AC_lnrho0 // from equationofstate
#define lnrho0 AC_lnrho0
real AC_lnTT0 // from equationofstate
#define lnTT0 AC_lnTT0
real AC_gamma_m1 // from equationofstate
#define gamma_m1 AC_gamma_m1
real AC_gamma1 // from equationofstate
#define gamma1 AC_gamma1

int AC_ishock_max // from shock
#define ishock_max AC_ishock_max
real AC_div_threshold // from shock
#define div_threshold AC_div_threshold
real AC_shock_linear // from shock
#define shock_linear AC_shock_linear
real AC_shock_div_pow // from shock
#define shock_div_pow AC_shock_div_pow
real AC_dt_div_pow // from shock
#define dt_div_pow AC_dt_div_pow
real AC_dtfactor // from shock
#define dtfactor AC_dtfactor
real AC_con_bias // from shock
#define con_bias AC_con_bias

real AC_GammaUV // from interstellar
#define GammaUV AC_GammaUV
real AC_cUV // from interstellar
#define cUV AC_cUV
real AC_T0UV // from interstellar
#define T0UV AC_T0UV
int AC_ncool // from interstellar
#define ncool AC_ncool
real AC_heatingfunction_scale[2] // from interstellar
#define heatingfunction_scale AC_heatingfunction_scale
gmem real AC_heat_SNI_profile[AC_nz] // from interstellar
#define heat_SNI_profile AC_heat_SNI_profile
gmem real AC_heat_SNII_profile[AC_nz] // from interstellar
#define heat_SNII_profile AC_heat_SNII_profile
real AC_lncoolT[11] // from interstellar
#define lncoolT AC_lncoolT
real AC_lncoolH[11] // from interstellar
#define lncoolH AC_lncoolH
real AC_coolB[11] // from interstellar
#define coolB AC_coolB
real AC_heating_rate_code // from interstellar
#define heating_rate_code AC_heating_rate_code

real AC_nu // from viscosity
#define nu AC_nu
real AC_zeta // from viscosity
#define zeta AC_zeta
real AC_nu_hyper3 // from viscosity
#define nu_hyper3 AC_nu_hyper3
real AC_nu_shock // from viscosity
#define nu_shock AC_nu_shock
int AC_lvisc_nu_const // from viscosity
#define lvisc_nu_const AC_lvisc_nu_const
int AC_lvisc_hyper3_nu_const // from viscosity
#define lvisc_hyper3_nu_const AC_lvisc_hyper3_nu_const
int AC_lvisc_nu_shock // from viscosity
#define lvisc_nu_shock AC_lvisc_nu_shock
