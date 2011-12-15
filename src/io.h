!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)

  private

  public :: lun_input, lun_output, lcollective_IO
  public :: register_io
  public :: output_snap, output_snap_finalize
  public :: input_snap, input_snap_finalize
  public :: output_globals, input_globals
  public :: init_write_persist, write_persist, write_persist_id
  public :: read_persist

  public :: wgrid, rgrid, wtime, rtime
  public :: wproc_bounds, rproc_bounds
  public :: directory_names, log_filename_to_file

