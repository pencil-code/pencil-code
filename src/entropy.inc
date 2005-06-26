!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_entropy, initialize_entropy
  public :: read_entropy_init_pars, write_entropy_init_pars
  public :: read_entropy_run_pars, write_entropy_run_pars
  public :: rprint_entropy
  public :: init_ss, dss_dt
  public :: pencil_criteria_entropy, pencil_interdep_entropy
  public :: calc_pencils_entropy

!ajwm SHOULDN'T BE SHARED
  public :: hcond0,hcond1,Fbot,FbotKbot,Ftop,Kbot,FtopKtop,chi, lmultilayer
!!  public :: nu_turb
  public :: lheatc_chiconst

