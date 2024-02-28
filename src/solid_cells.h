!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
     
private

public :: initialize_solid_cells
public :: init_solid_cells
public :: freeze_solid_cells
public :: update_solid_cells
public :: update_solid_cells_pencil
public :: read_solid_cells_init_pars
public :: read_solid_cells_run_pars
public :: write_solid_cells_init_pars
public :: write_solid_cells_run_pars
public :: in_solid_cell
public :: dsolid_dt, calc_diagnostics_solid
public :: dsolid_dt_integrate
public :: rprint_solid_cells
public :: pencil_criteria_solid_cells
public :: close_interpolation
public :: solid_cells_clean_up 
public :: register_solid_cells
public :: solid_cells_timestep_first
public :: solid_cells_timestep_second
public :: time_step_ogrid
public :: wsnap_ogrid

public :: r_ogrid, r_int_outer
public :: xorigo_ogrid
public :: map_nearest_grid_ogrid
public :: interpolate_particles_ogrid
public :: sc_diags_reductions
public :: sc_init_reduc_pointers
