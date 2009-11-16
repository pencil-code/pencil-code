!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private
 
  public :: interpolate_linear, interpolate_quadratic
  public :: interpolate_quadratic_spline
  public :: map_nearest_grid, map_xxp_grid, map_vvp_grid
  public :: sort_particles_imn, particle_pencil_index
  public :: shepherd_neighbour
  public :: interpolation_consistency_check
  public :: interpolate_quantities, cleanup_interpolated_quantities
