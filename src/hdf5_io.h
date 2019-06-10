  public :: initialize_hdf5, finalize_hdf5, file_open_hdf5, file_close_hdf5, create_group_hdf5
  public :: exists_in_hdf5, input_hdf5, output_hdf5, output_hdf5_double, output_dim
  public :: index_append, particle_index_append, pointmass_index_append, index_get, index_reset

  ! file location settings
  character(len=*), parameter :: index_pro = 'index.pro'
  character(len=*), parameter :: particle_index_pro = 'particle_index.pro'
  character(len=*), parameter :: pointmass_index_pro = 'pointmass_index.pro'

