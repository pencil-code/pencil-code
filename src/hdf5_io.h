  public :: initialize_hdf5, finalize_hdf5, file_open_hdf5, file_close_hdf5, create_group_hdf5
  public :: exists_in_hdf5, input_hdf5, output_hdf5, output_hdf5_double, wdim, input_dim
  public :: index_append, particle_index_append, pointmass_index_append, index_get, index_reset
  public :: input_profile, output_profile
  public :: input_slice, output_slice, output_slice_position
  public :: output_timeseries
  public :: output_average, trim_average

  interface output_average
    module procedure output_average_1D
    module procedure output_average_2D
    module procedure output_average_1D_chunked
    module procedure output_average_phi
  endinterface

  interface wdim
    module procedure wdim_default_grid
    module procedure wdim_default
    module procedure wdim
  endinterface

  ! file location settings
  character(len=*), parameter :: index_pro = 'index.pro'
  character(len=*), parameter :: particle_index_pro = 'particle_index.pro'
  character(len=*), parameter :: pointmass_index_pro = 'pointmass_index.pro'

