!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_forcing, initialize_forcing
  public :: read_forcing_run_pars,  write_forcing_run_pars
  public :: output_persistent_forcing, input_persistent_forcing
  public :: rprint_forcing
  public :: addforce
  public :: forcing_continuous, forcing_cont_after_boundary
  public :: lhydro_forcing, ltestflow_forcing
  
  public :: pencil_criteria_forcing, pencil_interdep_forcing
  public :: calc_pencils_forcing
  public :: forcing_clean_up
  public :: forcing_cont
  public :: forcing_pars_hel

  public :: pushpars2c, pushdiags2c

  public :: n_forcing_cont  ! should be protected
