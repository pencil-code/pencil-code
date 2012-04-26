!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)

  private

  public :: lun_input, lun_output, lcollective_IO, IO_strategy
  public :: register_io, fseek_pos, output_form
  public :: output_snap, output_snap_finalize
  public :: input_snap, input_snap_finalize
  public :: output_globals, input_globals
  public :: init_write_persist, write_persist, write_persist_id
  public :: init_read_persist, read_persist, read_persist_id

  public :: wgrid, rgrid
  public :: wproc_bounds, rproc_bounds
  public :: directory_names, log_filename_to_file

