!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)

  private

  public :: lun_input, lun_output, lcollective_IO, IO_strategy
  public :: register_io, finalize_io
  public :: output_snap, output_snap_finalize, output_part_snap, output_pointmass
  public :: output_stalker_init, output_stalker, output_part_finalize
  public :: input_snap, input_snap_finalize, input_part_snap, input_pointmass
  public :: output_globals, input_globals, output_slice_position, output_slice
  public :: init_write_persist, write_persist, write_persist_id
  public :: init_read_persist, read_persist, read_persist_id
  public :: output_timeseries

  public :: wgrid, rgrid, wdim, output_profile, input_profile
  public :: output_average, trim_average
  public :: wproc_bounds, rproc_bounds
  public :: directory_names, log_filename_to_file

  interface wdim
    module procedure wdim_default_grid
    module procedure wdim_default
    module procedure wdim
  endinterface

  interface output_average
    module procedure output_average_1D
    module procedure output_average_2D
    module procedure output_average_1D_chunked
    module procedure output_average_phi
  endinterface

