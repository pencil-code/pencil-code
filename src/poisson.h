!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: initialize_poisson
  public :: inverse_laplacian
  public :: inverse_laplacian_fft_z
  public :: inverse_laplacian_z_2nd_neumann
  public :: inverse_laplacian_semispectral
  public :: read_poisson_init_pars, write_poisson_init_pars
  public :: read_poisson_run_pars, write_poisson_run_pars
  public :: get_acceleration
