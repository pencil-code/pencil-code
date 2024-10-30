!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
!
! functions
!
    public :: initialize_training
    public :: register_training
    public :: training_before_boundary
    public :: finalize_training
    public :: read_training_run_pars
    public :: write_training_run_pars
    public :: calc_diagnostics_training, rprint_training
    public :: div_reynolds_stress
    public :: get_slices_training
