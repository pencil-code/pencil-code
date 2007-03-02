!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)

  private
                                                                                                       
  public :: eoscalc,pressure_gradient,temperature_gradient,get_cp1
  public :: temperature_laplacian
  public :: ilnrho_ss, ilnrho_lnTT, ilnrho_pp, ilnrho_ee
  public :: eosperturb
  public :: get_soundspeed
  public :: getmu
  public :: getdensity
                                                                                                       
  public :: register_eos
  public :: initialize_eos, units_eos
  public :: rprint_eos
  public :: read_eos_init_pars, write_eos_init_pars
  public :: read_eos_run_pars,  write_eos_run_pars

  public :: select_eos_variable
                                                                                                       
  public :: pencil_criteria_eos, pencil_interdep_eos
  public :: calc_pencils_eos

  public :: ioncalc, ioninit
  public :: temperature_hessian

! Boundary conditions
  public :: bc_ss_flux,bc_ss_temp_old,bc_ss_energy
  public :: bc_ss_temp_x,bc_ss_temp_y,bc_ss_temp_z,bc_ss_temp2_z
  public :: bc_ss_stemp_x,bc_ss_stemp_y,bc_ss_stemp_z
  public :: bc_lnrho_temp_z,bc_lnrho_pressure_z
  public :: bc_stellar_surface,bc_stellar_surface_2
  public :: bc_lnrho_hydrostatic_z
! Initial conditions
  public :: isothermal_entropy,isothermal_lnrho_ss

!ajwm SHOULDN'T BE PUBLIC
  public :: cs0,cs20,lnrho0,rho0
!Shouldn't be public, certainly means don't add anymore!!
!,mu,Rgas    BREAKS THE AUTO-TEST

  public :: gamma,gamma1,gamma11,cs2top,cs2bot
  public :: beta_glnrho_global, beta_glnrho_scaled
  public :: cs2cool
  public :: mpoly, mpoly0, mpoly1, mpoly2
  public :: isothtop
