!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)

  private
                                                                                                       
  public :: eoscalc,pressure_gradient,temperature_gradient,get_cp1
  public :: temperature_laplacian
  public :: ilnrho_ss, ilnrho_lnTT, ilnrho_pp, ilnrho_ee, ilnrho_TT
  public :: ipp_ss,ipp_cs2
  public :: eosperturb
  public :: get_soundspeed
  public :: getmu
  public :: getdensity, gettemperature, getpressure
  public :: get_average_pressure
                                                                                                       
  public :: register_eos
  public :: initialize_eos, units_eos
  public :: rprint_eos, get_slices_eos
  public :: read_eos_init_pars, write_eos_init_pars
  public :: read_eos_run_pars,  write_eos_run_pars

  public :: select_eos_variable
                                                                                                       
  public :: pencil_criteria_eos, pencil_interdep_eos
  public :: calc_pencils_eos

  public :: ioncalc, ioninit
  public :: temperature_hessian

!ajwm SHOULDN'T BE PUBLIC
  public :: cs0,cs20,lnrho0,rho0
!Shouldn't be public, certainly means don't add anymore!!
!,mu,Rgas    BREAKS THE AUTO-TEST

  public :: gamma,gamma_m1,gamma_inv,cs2top,cs2bot,cs2top_ini,dcs2top_ini
  public :: beta_glnrho_global, beta_glnrho_scaled
  public :: cs2cool
  public :: mpoly, mpoly0, mpoly1, mpoly2
  public :: isothtop
  public :: ieos_profile,profz_eos,dprofz_eos

