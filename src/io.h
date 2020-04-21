!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)

  private

  public :: lun_input, lun_output, lcollective_IO, IO_strategy
  public :: register_io, finalize_io
  public :: output_snap, output_snap_finalize, output_part_snap, output_pointmass
  public :: output_stalker_init, output_stalker, output_part_finalize
  public :: input_snap, input_snap_finalize, input_part_snap, input_pointmass
  public :: output_globals, input_globals
  public :: init_write_persist, write_persist, write_persist_id, &
            write_persist_logical_0D, write_persist_logical_1D, write_persist_int_0D, &
            write_persist_int_1D, write_persist_real_0D, write_persist_real_1D
  public :: init_read_persist, read_persist, read_persist_id, persist_exists
  public :: read_persist_logical_0D, read_persist_logical_1D, read_persist_int_0D, &
            read_persist_int_1D, read_persist_real_0D, read_persist_real_1D

  public :: wgrid, rgrid
  public :: wproc_bounds, rproc_bounds
  public :: directory_names, log_filename_to_file
!
  ! define unique logical unit number for input and output calls
  integer, parameter :: lun_input=88
  integer, parameter :: lun_output=91

  logical :: persist_initialized = .false.

  interface write_persist
    module procedure write_persist_logical_0D
    module procedure write_persist_logical_1D
    module procedure write_persist_int_0D                                                                                       
    module procedure write_persist_int_1D
    module procedure write_persist_real_0D
    module procedure write_persist_real_1D
  endinterface
!        
  interface read_persist
    module procedure read_persist_logical_0D
    module procedure read_persist_logical_1D
    module procedure read_persist_int_0D
    module procedure read_persist_int_1D
    module procedure read_persist_real_0D
    module procedure read_persist_real_1D
  endinterface
