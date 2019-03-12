!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: register_particles_lyapunov 
  public :: initialize_particles_lyapunov,init_particles_lyapunov
  public :: read_plyapunov_init_pars, write_plyapunov_init_pars
  public :: read_plyapunov_run_pars, write_plyapunov_run_pars
  public :: particles_stochastic_lyapunov
  public :: rprint_particles_lyapunov
  public :: dlyapunov_dt,dlyapunov_dt_pencil
  public :: calc_pencils_par_lyapunov
