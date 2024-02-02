// automatically generated; do not edit!
const int n_pars_forcing=9;
 EXTERN void *p_par_forcing[n_pars_forcing];
 extern "C" void *pushpars2c$forcing_(void *(&p_par)[n_pars_forcing]);
#define k1_ff (*((real *) p_par_forcing[1-1]))
#define tforce_stop (*((real *) p_par_forcing[2-1]))
#define iforcing_zsym (*((real *) p_par_forcing[3-1]))
#define profx_ampl ((real *) p_par_par_forcing[4-1]) // (nx)
#define profy_ampl ((real *) p_par_par_forcing[5-1]) // (my)
#define profz_ampl ((real *) p_par_par_forcing[6-1]) // (mz)
#define profx_hel ((real *) p_par_par_forcing[7-1]) // (nx)
#define profy_hel ((real *) p_par_par_forcing[8-1]) // (my)
#define profz_hel ((real *) p_par_par_forcing[9-1]) // (mz)
const int n_pars_hydro=1;
 EXTERN void *p_par_hydro[n_pars_hydro];
 extern "C" void *__hydro_MOD_pushpars2c(void *(&p_par)[n_pars_hydro]);
#define lpressuregradient_gas (*((int *) p_par_hydro[1-1]))
const int n_pars_magnetic=1;
 EXTERN void *p_par_magnetic[n_pars_magnetic];
 extern "C" void *__magnetic_MOD_pushpars2c(void *(&p_par)[n_pars_magnetic]);
#define eta (*((real *) p_par_magnetic[1-1]))
const int n_pars_energy=3;
 EXTERN void *p_par_energy[n_pars_energy];
 extern "C" void *__energy_MOD_pushpars2c(void *(&p_par)[n_pars_energy]);
#define chi (*((real *) p_par_energy[1-1]))
#define dlnhcond_prof ((real *) p_par_par_energy[2-1]) // (nz)
#define hcond_prof ((real *) p_par_par_energy[3-1]) // (nz)
const int n_pars_equationofstate=6;
 EXTERN void *p_par_equationofstate[n_pars_equationofstate];
 extern "C" void *__equationofstate_MOD_pushpars2c(void *(&p_par)[n_pars_equationofstate]);
#define cs20 (*((real *) p_par_equationofstate[1-1]))
#define gamma (*((real *) p_par_equationofstate[2-1]))
#define cv (*((real *) p_par_equationofstate[3-1]))
#define cp (*((real *) p_par_equationofstate[4-1]))
#define lnrho0 (*((real *) p_par_equationofstate[5-1]))
#define lnTT0 (*((real *) p_par_equationofstate[6-1]))
const int n_pars_viscosity=2;
 EXTERN void *p_par_viscosity[n_pars_viscosity];
 extern "C" void *__viscosity_MOD_pushpars2c(void *(&p_par)[n_pars_viscosity]);
#define nu (*((real *) p_par_viscosity[1-1]))
#define zeta (*((real *) p_par_viscosity[2-1]))
#define alpha_ts ((real *) p_par_par_timestep[1-1]) // (3)
#define beta_ts ((real *) p_par_par_timestep[2-1]) // (3)
