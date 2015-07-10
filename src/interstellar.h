!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_interstellar, initialize_interstellar
  public :: rprint_interstellar
  public :: init_interstellar
  public :: get_slices_interstellar
  public :: read_interstellar_init_pars, write_interstellar_init_pars
  public :: read_interstellar_run_pars,  write_interstellar_run_pars
  public :: pencil_criteria_interstellar
!  public :: interstellar_before_boundary
  public :: calc_heat_cool_interstellar, check_SN
!  public :: calc_snr_damping
!  public :: calc_snr_damp_int
!  public :: calc_snr_unshock
  public :: input_persistent_interstellar, output_persistent_interstellar
  public :: addmassflux

