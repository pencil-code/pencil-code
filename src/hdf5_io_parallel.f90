! $Id$
!
!  This module takes care of massive parallel HDF5 file Input/Output.
!  We use here only F95 features for HPC-friendly behaviour.
!
module HDF5_IO
!
  use Cdata
  use HDF5
  use Messages, only: fatal_error
!
  implicit none
!
  interface output_hdf5
    module procedure output_hdf5_string
    module procedure output_hdf5_int
    module procedure output_hdf5_0D
    module procedure output_hdf5_1D
    module procedure output_hdf5_3D
    module procedure output_hdf5_4D
  endinterface
!
  interface input_hdf5
    module procedure input_hdf5_0D
    module procedure input_hdf5_1D
    module procedure input_hdf5_3D
    module procedure input_hdf5_4D
  endinterface
!
  include 'hdf5_io.h'
  include 'mpif.h'
!
  private
!
  integer :: h5_err
  integer(HID_T) :: h5_file, h5_dset, h5_plist, h5_fspace, h5_mspace, h5_dspace, h5_ntype, h5_group
  integer, parameter :: n_dims = 3
  integer(kind=8), dimension(n_dims+1) :: local_size, local_subsize, local_start
  integer(kind=8), dimension(n_dims+1) :: global_size, global_subsize, global_start
!
  contains
!***********************************************************************
    subroutine initialize_hdf5
!
!  Initialize the HDF IO.
!
!  28-Oct-2016/PABoudin: coded
!
      use Mpicomm, only: mpi_precision
!
      ! dimensions for local data portion with ghost layers
      local_size(1) = mx
      local_size(2) = my
      local_size(3) = mz
      local_size(4:n_dims+1) = 1
!
      ! dimensions for local data portion without ghost layers
      local_subsize(1) = nx
      local_subsize(2) = ny
      local_subsize(3) = nz
      local_subsize(4:n_dims+1) = 1
!
      ! include the ghost layers only on the outer box boundaries
      if (lfirst_proc_x) local_subsize(1) = local_subsize(1) + nghost
      if (lfirst_proc_y) local_subsize(2) = local_subsize(2) + nghost
      if (lfirst_proc_z) local_subsize(3) = local_subsize(3) + nghost
      if (llast_proc_x)  local_subsize(1) = local_subsize(1) + nghost
      if (llast_proc_y)  local_subsize(2) = local_subsize(2) + nghost
      if (llast_proc_z)  local_subsize(3) = local_subsize(3) + nghost
!
      ! displacements in HDF5 use C-like format, ie. they start from zero
      local_start(1) = l1 - 1
      local_start(2) = m1 - 1
      local_start(3) = n1 - 1
      local_start(4:n_dims+1) = 0

      ! include lower ghost cells on the lower edge
      ! (upper ghost cells are taken care of by the increased 'local_subsize')
      if (lfirst_proc_x) local_start(1) = local_start(1) - nghost
      if (lfirst_proc_y) local_start(2) = local_start(2) - nghost
      if (lfirst_proc_z) local_start(3) = local_start(3) - nghost
!
      ! size of the data in the global file
      global_size(1) = mxgrid
      global_size(2) = mygrid
      global_size(3) = mzgrid
      global_size(4:n_dims+1) = 1
!
      ! starting position of this processor's data portion in the global file
      global_start(1) = nghost + ipx*nx
      global_start(2) = nghost + ipy*ny
      global_start(3) = nghost + ipz*nz
      global_start(4:n_dims+1) = 0
!
      ! include lower ghost layers on the lower edge
      ! (upper ghost cells are taken care of by the increased 'local_subsize')
      if (lfirst_proc_x) global_start(1) = global_start(1) - nghost
      if (lfirst_proc_y) global_start(2) = global_start(2) - nghost
      if (lfirst_proc_z) global_start(3) = global_start(3) - nghost
!
      ! initialize parallel HDF5 Fortran libaray
      call h5open_f (h5_err)
      if (h5_err /= 0) call fatal_error ('initialize_hdf5', 'initialize parallel HDF5 library', .true.)
      if (mpi_precision == MPI_REAL) then
        h5_ntype = H5T_NATIVE_REAL
      else
        h5_ntype = H5T_NATIVE_DOUBLE
      endif
!
    endsubroutine initialize_hdf5
!***********************************************************************
    subroutine finalize_hdf5
!
      ! close the HDF5 library
      call h5close_f (h5_err)
      if (h5_err /= 0) call fatal_error ('finalize_hdf5', 'close parallel HDF5 library', .true.)
!
    endsubroutine finalize_hdf5
!***********************************************************************
    subroutine file_open_hdf5(file, truncate, global, read_only)
!
      character (len=*), intent(in) :: file
      logical, optional, intent(in) :: truncate
      logical, optional, intent(in) :: global
      logical, optional, intent(in) :: read_only
!
      logical :: ltrunc, lglobal, lread_only
      integer :: h5_read_mode
!
      lread_only = .false.
      if (present (read_only)) lread_only = read_only
      h5_read_mode = H5F_ACC_RDWR_F
      if (lread_only) h5_read_mode = H5F_ACC_RDONLY_F
!
      ltrunc = .true.
      if (present (truncate)) ltrunc = truncate
      if (lread_only) ltrunc = .false.
!
      lglobal = .true.
      if (present (global)) lglobal = global
!
      if (lglobal) then
        ! setup file access property list
        call h5pcreate_f (H5P_FILE_ACCESS_F, h5_plist, h5_err)
        if (h5_err /= 0) call fatal_error ('file_open_hdf5', 'create global file access property list', .true.)
        call h5pset_fapl_mpio_f (h5_plist, MPI_COMM_WORLD, MPI_INFO_NULL, h5_err)
        if (h5_err /= 0) call fatal_error ('file_open_hdf5', 'modify global file access property list', .true.)
!
        if (ltrunc) then
          ! create empty (or truncated) HDF5 file
          call h5fcreate_f (trim (file), H5F_ACC_TRUNC_F, h5_file, h5_err, access_prp=h5_plist)
          if (h5_err /= 0) call fatal_error ('file_open_hdf5', 'create global file "'//trim (file)//'"', .true.)
        else
          ! open existing HDF5 file
          call h5fopen_f (trim (file), h5_read_mode, h5_file, h5_err, access_prp=h5_plist)
          if (h5_err /= 0) call fatal_error ('file_open_hdf5', 'open global file "'//trim (file)//'"', .true.)
        endif
!
        call h5pclose_f (h5_plist, h5_err)
        if (h5_err /= 0) call fatal_error ('file_open_hdf5', 'close global file access property list', .true.)
      else
        call h5fopen_f (trim (file), h5_read_mode, h5_file, h5_err)
        if (h5_err /= 0) call fatal_error ('file_open_hdf5', 'open local file "'//trim (file)//'"', .true.)
      endif
!
    endsubroutine file_open_hdf5
!***********************************************************************
    subroutine file_close_hdf5
!
      call h5fclose_f (h5_file, h5_err)
      if (h5_err /= 0) call fatal_error ('file_close_hdf5', 'close file', .true.)
!
    endsubroutine file_close_hdf5
!***********************************************************************
    subroutine create_group_hdf5(name)
!
      character (len=*), intent(in) :: name
!
      call h5gcreate_f (h5_file, trim (name), h5_group, h5_err)
      if (h5_err /= 0) call fatal_error ('create_group_hdf5', 'create group "'//trim (name)//'"', .true.)
      call h5gclose_f (h5_group, h5_err)
      if (h5_err /= 0) call fatal_error ('create_group_hdf5', 'close group "'//trim (name)//'"', .true.)
!
    endsubroutine create_group_hdf5
!***********************************************************************
    subroutine input_hdf5_0D(name, data)
!
      character (len=*), intent(in) :: name
      real, intent(out) :: data
!
      real, dimension(1) :: read
!
      call input_hdf5_1D (name, read, 1)
      data = read(1)
!
    endsubroutine input_hdf5_0D
!***********************************************************************
    subroutine input_hdf5_1D(name, data, nv)
!
!  Read HDF5 dataset as scalar or array.
!
!  26-Oct-2016/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (nv), intent(out) :: data
!
      integer(HSIZE_T), dimension(1) :: size
!
      size = (/ nv /)
!
      ! open dataset
      call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'open dataset "'//trim (name)//'"', .true.)
      ! read dataset
      call h5dread_f (h5_dset, h5_ntype, data, size, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'read data "'//trim (name)//'"', .true.)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'close dataset "'//trim (name)//'"', .true.)
!
    endsubroutine input_hdf5_1D
!***********************************************************************
    subroutine input_hdf5_3D(name, data)
!
      character (len=*), intent(in) :: name
      real, dimension (mx,my,mz), intent(out) :: data
!
      call input_hdf5_4D (name, data, 1)
!
    endsubroutine input_hdf5_3D
!***********************************************************************
    subroutine input_hdf5_4D(name, data, nv)
!
!  Read HDF5 dataset from a distributed 3D or 4D array.
!
!  26-Oct-2016/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(out) :: data
!
      integer(kind=8), dimension (n_dims+1) :: h5_stride, h5_count
!
      global_size(n_dims+1) = nv
      local_size(n_dims+1) = nv
      local_subsize(n_dims+1) = nv
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (n_dims+1, local_size, h5_mspace, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'create local memory space "'//trim (name)//'"', .true.)
!
      ! open the dataset
      call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'open dataset "'//trim (name)//'"', .true.)
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'get dataset for file space "'//trim (name)//'"', .true.)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start, h5_count, h5_err, h5_stride, local_subsize)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'select hyperslab within file "'//trim (name)//'"', .true.)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start, h5_count, h5_err, h5_stride, local_subsize)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'select hyperslab within file "'//trim (name)//'"', .true.)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'set data transfer properties "'//trim (name)//'"', .true.)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'select collective IO "'//trim (name)//'"', .true.)
!
      ! collectively read the data
      call h5dread_f (h5_dset, h5_ntype, data, &
          global_size, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'read dataset "'//trim (name)//'"', .true.)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'close file space "'//trim (name)//'"', .true.)
      call h5sclose_f (h5_mspace, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'close memory space "'//trim (name)//'"', .true.)
      call h5dclose_f (h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'close dataset "'//trim (name)//'"', .true.)
      call h5pclose_f (h5_plist, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'close parameter list "'//trim (name)//'"', .true.)
!
    endsubroutine input_hdf5_4D
!***********************************************************************
    subroutine output_hdf5_string(name, data)
!
      character (len=*), intent(in) :: name
      character (len=*), intent(in) :: data
!
      integer(HID_T) :: h5_strtype
      integer(HSIZE_T), dimension(2) :: size
      character (len=len(data)), dimension(1) :: str_data
      integer(SIZE_T), dimension(1) :: str_len
!
      str_len(1) = len_trim (data)
      size(1) = str_len(1)
      size(2) = 1
      str_data(1) = data
!
      ! create data space
      call H5Tcopy_f (H5T_STRING, h5_strtype, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'copy string data space type "'//trim (name)//'"', .true.)
      call H5Tset_strpad_f (h5_strtype, H5T_STR_NULLPAD_F, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'modify string data space type "'//trim (name)//'"', .true.)
      call h5screate_simple_f (1, size(1), h5_dspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create string data space "'//trim (name)//'"', .true.)
      ! create dataset
      call h5dcreate_f (h5_file, trim (name), h5_strtype, h5_dspace, h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create string dataset "'//trim (name)//'"', .true.)
      ! write dataset
      call h5dwrite_vl_f (h5_dset, h5_strtype, str_data, size, str_len, h5_err, h5_dspace)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'write string data "'//trim (name)//'"', .true.)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close string dataset "'//trim (name)//'"', .true.)
      call h5sclose_f (h5_dspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close string data space "'//trim (name)//'"', .true.)
!
    endsubroutine output_hdf5_string
!***********************************************************************
    subroutine output_hdf5_int(name, data)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: data
!
      integer(HID_T) :: h5_inttype
      integer(kind=8), dimension(1) :: size = (/ 1 /)
!
      ! create data space
      call h5screate_simple_f (1, size(1), h5_dspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create integer data space "'//trim (name)//'"', .true.)
      ! create dataset
      call h5dcreate_f (h5_file, trim (name), H5T_NATIVE_INTEGER, h5_dspace, h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create integer dataset "'//trim (name)//'"', .true.)
      ! write dataset
      call h5dwrite_f (h5_dset, H5T_NATIVE_INTEGER, data, size, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'write integer data "'//trim (name)//'"', .true.)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close integer dataset "'//trim (name)//'"', .true.)
      call h5sclose_f (h5_dspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close integer data space "'//trim (name)//'"', .true.)
!
    endsubroutine output_hdf5_int
!***********************************************************************
    subroutine output_hdf5_0D(name, data)
!
      character (len=*), intent(in) :: name
      real, intent(in) :: data
!
      call output_hdf5_1D (name, (/ data /), 1)
!
    endsubroutine output_hdf5_0D
!***********************************************************************
    subroutine output_hdf5_1D(name, data, nv)
!
!  Write HDF5 dataset as scalar or array.
!
!  24-Oct-2016/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (nv), intent(in) :: data
!
      integer(kind=8), dimension(1) :: size
!
      size = (/ nv /)
!
      ! create data space
      if (nv <= 1) then
        call h5screate_f (H5S_SCALAR_F, h5_dspace, h5_err)
        if (h5_err /= 0) call fatal_error ('output_hdf5', 'create scalar data space "'//trim (name)//'"', .true.)
        call h5sset_extent_simple_f (h5_dspace, 0, size(1), size(1), h5_err) 
      else
        call h5screate_f (H5S_SIMPLE_F, h5_dspace, h5_err)
        if (h5_err /= 0) call fatal_error ('output_hdf5', 'create simple data space "'//trim (name)//'"', .true.)
        call h5sset_extent_simple_f (h5_dspace, 1, size, size, h5_err) 
      endif
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'set data space extent "'//trim (name)//'"', .true.)
      ! create dataset
      call h5dcreate_f (h5_file, trim (name), h5_ntype, h5_dspace, h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create dataset "'//trim (name)//'"', .true.)
      ! write dataset
      if (nv <= 1) then
        call h5dwrite_f (h5_dset, h5_ntype, data(1), size, h5_err)
      else
        call h5dwrite_f (h5_dset, h5_ntype, data, size, h5_err)
      endif
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'write data "'//trim (name)//'"', .true.)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close dataset "'//trim (name)//'"', .true.)
      call h5sclose_f (h5_dspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close data space "'//trim (name)//'"', .true.)
!
    endsubroutine output_hdf5_1D
!***********************************************************************
    subroutine output_hdf5_3D(name, data)
!
      character (len=*), intent(in) :: name
      real, dimension (mx,my,mz), intent(in) :: data
!
      call output_hdf5_4D (name, data, 1)
!
    endsubroutine output_hdf5_3D
!***********************************************************************
    subroutine output_hdf5_4D(name, data, nv)
!
!  Write HDF5 dataset from a distributed 3D or 4D array.
!
!  26-Oct-2016/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: data
!
      integer(kind=8), dimension (n_dims+1) :: h5_stride, h5_count
!
      global_size(n_dims+1) = nv
      local_size(n_dims+1) = nv
      local_subsize(n_dims+1) = nv
!
      ! define 'file-space' to indicate the data portion in the global file
      call h5screate_simple_f (n_dims+1, global_size, h5_fspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create global file space "'//trim (name)//'"', .true.)
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (n_dims+1, local_size, h5_mspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create local memory space "'//trim (name)//'"', .true.)
!
      ! create the dataset
      call h5dcreate_f (h5_file, trim (name), h5_ntype, h5_fspace, h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create dataset "'//trim (name)//'"', .true.)
      call h5sclose_f (h5_fspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close global file space "'//trim (name)//'"', .true.)
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'get dataset for file space "'//trim (name)//'"', .true.)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start, h5_count, h5_err, h5_stride, local_subsize)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select hyperslab within file "'//trim (name)//'"', .true.)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start, h5_count, h5_err, h5_stride, local_subsize)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select hyperslab within file "'//trim (name)//'"', .true.)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'set data transfer properties "'//trim (name)//'"', .true.)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select collective IO "'//trim (name)//'"', .true.)
!
      ! collectively write the data
      call h5dwrite_f (h5_dset, h5_ntype, data, &
          global_size, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'write dataset "'//trim (name)//'"', .true.)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close file space "'//trim (name)//'"', .true.)
      call h5sclose_f (h5_mspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close memory space "'//trim (name)//'"', .true.)
      call h5dclose_f (h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close dataset "'//trim (name)//'"', .true.)
      call h5pclose_f (h5_plist, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close parameter list "'//trim (name)//'"', .true.)
!
    endsubroutine output_hdf5_4D
!***********************************************************************
endmodule HDF5_IO
