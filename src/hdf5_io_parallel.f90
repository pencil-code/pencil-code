! $Id$
!
!  This module takes care of massive parallel HDF5 file Input/Output.
!  We use here only F95 features for HPC-friendly behaviour.
!
module HDF5_IO
!
  use Cdata
  use Cparam, only: mvar, maux, labellen
  use General, only: loptest, itoa
  use HDF5
  use Messages, only: fatal_error
  use Mpicomm, only: lroot, mpi_precision, mpiscan_int, mpibcast_int
!
  implicit none
!
  interface output_hdf5
    module procedure output_hdf5_string
    module procedure output_hdf5_int_0D
    module procedure output_hdf5_int_1D
    module procedure output_hdf5_0D
    module procedure output_hdf5_1D
    module procedure output_hdf5_3D
    module procedure output_hdf5_4D
  endinterface
!
  interface input_hdf5
    module procedure input_hdf5_int_0D
    module procedure input_hdf5_int_1D
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
  integer(kind=8), dimension(n_dims+1) :: global_size, global_start
  logical :: lcollective = .false., lwrite = .false.
!
  type element
    character(len=labellen) :: label
    integer :: component
    type (element), pointer :: previous
  endtype element
  type (element), pointer :: last => null(), last_particle => null()
!
  contains
!***********************************************************************
    subroutine initialize_hdf5
!
!  Initialize the HDF IO.
!
!  28-Oct-2016/PABoudin: coded
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
    subroutine file_open_hdf5(file, truncate, global, read_only, write)
!
      character (len=*), intent(in) :: file
      logical, optional, intent(in) :: truncate
      logical, optional, intent(in) :: global
      logical, optional, intent(in) :: read_only
      logical, optional, intent(in) :: write
!
      logical :: ltrunc, lread_only
      integer :: h5_read_mode
!
      if (lcollective .or. lwrite) call file_close_hdf5 ()
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
      lcollective = .true.
      if (present (global)) lcollective = global
      lwrite = lroot
      if (present (write)) lwrite = write
!
      if (lcollective) then
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
      elseif (lwrite) then
        call h5fopen_f (trim (file), h5_read_mode, h5_file, h5_err)
        if (h5_err /= 0) call fatal_error ('file_open_hdf5', 'open local file "'//trim (file)//'"', .true.)
      endif
!
    endsubroutine file_open_hdf5
!***********************************************************************
    subroutine file_close_hdf5
!
      if (.not. (lcollective .or. lwrite)) return
!
      call h5fclose_f (h5_file, h5_err)
      if (h5_err /= 0) call fatal_error ('file_close_hdf5', 'close file', .true.)
!
      lcollective = .false.
      lwrite = .false.
!
    endsubroutine file_close_hdf5
!***********************************************************************
    subroutine create_group_hdf5(name)
!
      character (len=*), intent(in) :: name
!
      if (.not. (lcollective .or. lwrite)) return
!
      call h5gcreate_f (h5_file, trim (name), h5_group, h5_err)
      if (h5_err /= 0) call fatal_error ('create_group_hdf5', 'create group "'//trim (name)//'"', .true.)
      call h5gclose_f (h5_group, h5_err)
      if (h5_err /= 0) call fatal_error ('create_group_hdf5', 'close group "'//trim (name)//'"', .true.)
!
    endsubroutine create_group_hdf5
!***********************************************************************
    subroutine input_hdf5_int_0D(name, data)
!
      character (len=*), intent(in) :: name
      integer, intent(out) :: data
!
      integer, dimension(1) :: read
!
      call input_hdf5_int_1D (name, read, 1)
      data = read(1)
!
    endsubroutine input_hdf5_int_0D
!***********************************************************************
    subroutine input_hdf5_int_1D(name, data, nv)
!
!  Read HDF5 dataset as scalar or array.
!
!  05-Jun-2017/Fred: coded based on input_hdf5_1D
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      integer, dimension (nv), intent(out) :: data
!
      integer(HSIZE_T), dimension(1) :: size
!
      size = (/ nv /)
!
      ! open dataset
      call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'open dataset "'//trim (name)//'"', .true.)
      ! read dataset
      call h5dread_f (h5_dset, H5T_NATIVE_INTEGER, data, size, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'read data "'//trim (name)//'"', .true.)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'close dataset "'//trim (name)//'"', .true.)
!
    endsubroutine input_hdf5_int_1D
!***********************************************************************
    subroutine input_hdf5_0D(name, data)
!
      character (len=*), intent(in) :: name
      real, intent(out) :: data
!
      real, dimension(1) :: input
!
      call input_hdf5_1D (name, input, 1)
      data = input(1)
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
!  Read HDF5 dataset from a distributed 3D array.
!
!  26-Oct-2016/PABourdin: coded
!
      character (len=*), intent(in) :: name
      real, dimension (mx,my,mz), intent(out) :: data
!
      integer(kind=8), dimension (n_dims) :: h5_stride, h5_count
      integer, parameter :: n = n_dims
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (n, local_size(1:n), h5_mspace, h5_err)
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
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start(1:n), h5_count, h5_err, h5_stride, local_subsize(1:n))
      if (h5_err /= 0) call fatal_error ('input_hdf5', 'select hyperslab within file "'//trim (name)//'"', .true.)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start(1:n), h5_count, h5_err, h5_stride, local_subsize(1:n))
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
    endsubroutine input_hdf5_3D
!***********************************************************************
    subroutine input_hdf5_4D(name, data, nv)
!
!  Read HDF5 dataset from a distributed 4D array.
!
!  26-Oct-2016/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(out) :: data
!
      integer(kind=8), dimension (n_dims+1) :: h5_stride, h5_count
      integer :: pos
!
      if (name == 'f') then
        ! write components of f-array
        do pos=1, nv
          call input_hdf5_3D ('data/'//index_get(pos), data(:,:,:,pos))
        enddo
        return
      endif
!
      ! read other 4D array
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
      if (lcollective) call fatal_error ('output_hdf5', 'string output requires local file "'//trim (name)//'"', .true.)
      if (.not. lwrite) return
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
    subroutine output_local_hdf5_int_1D(name, data, nv)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      integer, dimension(nv), intent(in) :: data
!
      integer(kind=8), dimension(1) :: size
!
      if (.not. lwrite) return
!
      size = (/ nv /)
!
      ! create data space
      call h5screate_simple_f (1, size, h5_dspace, h5_err)
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
    endsubroutine output_local_hdf5_int_1D
!***********************************************************************
    subroutine output_hdf5_int_0D(name, data)
!
!  Write HDF5 dataset as scalar from one or all processor.
!
!  22-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: data
!
      integer, dimension(1) :: output = (/ 1 /)
!
      output = data
      call output_hdf5_int_1D(name, output, 1, .true.)
!
    endsubroutine output_hdf5_int_0D
!***********************************************************************
    subroutine output_hdf5_int_1D(name, data, nv, same_size)
!
!  Write HDF5 dataset as array from one or all processors.
!
!  24-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      integer, dimension (nv), intent(in) :: data
      logical, optional, intent(in) :: same_size
!
      logical :: lsame_size
      integer :: total, offset, last
      integer(kind=8), dimension (1) :: local_size_1D, local_subsize_1D, local_start_1D
      integer(kind=8), dimension (1) :: global_size_1D, global_start_1D
      integer(kind=8), dimension (1) :: h5_stride, h5_count
!
      if (.not. lcollective) then
        call output_local_hdf5_int_1D(name, data, nv)
        return
      endif
!
      lsame_size = .false.
      if (present (same_size)) lsame_size = same_size
      if (lsame_size) then
        last = nv * (iproc + 1) - 1
        total = nv * ncpus
        offset = nv * iproc
      else
        call mpiscan_int(nv, offset)
        last = offset - 1
        total = offset
        offset = offset - nv
        call mpibcast_int(total, ncpus-1)
      endif
      local_start_1D = 0
      local_size_1D = nv
      local_subsize_1D = nv
      global_size_1D = total
      global_start_1D = offset
!
      ! define 'file-space' to indicate the data portion in the global file
      call h5screate_simple_f (1, global_size_1D, h5_fspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create global file space "'//trim (name)//'"', .true.)
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (1, local_size_1D, h5_mspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create local memory space "'//trim (name)//'"', .true.)
!
      ! create the dataset
      call h5dcreate_f (h5_file, trim (name), H5T_NATIVE_INTEGER, h5_fspace, h5_dset, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create dataset "'//trim (name)//'"', .true.)
      call h5sclose_f (h5_fspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'close global file space "'//trim (name)//'"', .true.)
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'get dataset for file space "'//trim (name)//'"', .true.)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select hyperslab within file "'//trim (name)//'"', .true.)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select hyperslab within memory "'//trim (name)//'"', .true.)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'set data transfer properties "'//trim (name)//'"', .true.)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select collective IO "'//trim (name)//'"', .true.)
!
      ! collectively write the data
      call h5dwrite_f (h5_dset, H5T_NATIVE_INTEGER, data, &
          global_size_1D, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
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
    endsubroutine output_hdf5_int_1D
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
    subroutine output_local_hdf5_1D(name, data, nv)
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
      if (.not. (lcollective .or. lwrite)) return
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
    endsubroutine output_local_hdf5_1D
!***********************************************************************
    subroutine output_hdf5_1D(name, data, nv, same_size)
!
!  Write HDF5 dataset as scalar or array.
!
!  24-Oct-2016/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (nv), intent(in) :: data
      logical, optional, intent(in) :: same_size
!
      logical :: lsame_size
      integer :: total, offset, last
      integer(kind=8), dimension (1) :: local_size_1D, local_subsize_1D, local_start_1D
      integer(kind=8), dimension (1) :: global_size_1D, global_start_1D
      integer(kind=8), dimension (1) :: h5_stride, h5_count
!
      if (.not. lcollective) then
        call output_local_hdf5_1D(name, data, nv)
        return
      endif
!
      lsame_size = .false.
      if (present (same_size)) lsame_size = same_size
      if (lsame_size) then
        last = nv * (iproc + 1) - 1
        total = nv * ncpus
        offset = nv * iproc
      else
        call mpiscan_int(nv, offset)
        last = offset - 1
        total = offset
        offset = offset - nv
        call mpibcast_int(total, ncpus-1)
      endif
      local_start_1D = 0
      local_size_1D = nv
      local_subsize_1D = nv
      global_size_1D = total
      global_start_1D = offset
!
      ! define 'file-space' to indicate the data portion in the global file
      call h5screate_simple_f (1, global_size_1D, h5_fspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create global file space "'//trim (name)//'"', .true.)
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (1, local_size_1D, h5_mspace, h5_err)
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
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select hyperslab within file "'//trim (name)//'"', .true.)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select hyperslab within memory "'//trim (name)//'"', .true.)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'set data transfer properties "'//trim (name)//'"', .true.)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select collective IO "'//trim (name)//'"', .true.)
!
      ! collectively write the data
      call h5dwrite_f (h5_dset, h5_ntype, data, &
          global_size_1D, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
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
    endsubroutine output_hdf5_1D
!***********************************************************************
    subroutine output_hdf5_3D(name, data)
!
!  Write HDF5 dataset from a distributed 3D array.
!
!  17-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: name
      real, dimension (mx,my,mz), intent(in) :: data
!
      integer(kind=8), dimension (n_dims) :: h5_stride, h5_count
      integer, parameter :: n = n_dims
!
      if (.not. lcollective) call fatal_error ('output_hdf5', '3D array output requires global file "'//trim (name)//'"', .true.)
!
      ! define 'file-space' to indicate the data portion in the global file
      call h5screate_simple_f (n, global_size(1:n), h5_fspace, h5_err)
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'create global file space "'//trim (name)//'"', .true.)
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (n, local_size(1:n), h5_mspace, h5_err)
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
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start(1:n), h5_count, h5_err, h5_stride, local_subsize(1:n))
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select hyperslab within file "'//trim (name)//'"', .true.)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start(1:n), h5_count, h5_err, h5_stride, local_subsize(1:n))
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select hyperslab within memory "'//trim (name)//'"', .true.)
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
    endsubroutine output_hdf5_3D
!***********************************************************************
    subroutine output_hdf5_4D(name, data, nv)
!
!  Write HDF5 dataset from a distributed 4D array.
!
!  26-Oct-2016/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: data
!
      integer(kind=8), dimension (n_dims+1) :: h5_stride, h5_count
      integer :: pos
!
      if (.not. lcollective) call fatal_error ('output_hdf5', '4D array output requires global file "'//trim (name)//'"', .true.)
!
      if (name == 'f') then
        ! write components of f-array
        call create_group_hdf5 ('data')
        do pos=1, nv
          call output_hdf5_3D ('data/'//index_get(pos), data(:,:,:,pos))
        enddo
        return
      endif
!
      ! write other 4D array
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
      if (h5_err /= 0) call fatal_error ('output_hdf5', 'select hyperslab within memory "'//trim (name)//'"', .true.)
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
    subroutine set_part_dist_hdf5(name, nv)
!
      character (len=*), intent(in) :: name
      integer, intent(out) :: nv
!
      call fatal_error ('set_part_dist_hdf5', 'Not yet implemented.')
!
    endsubroutine set_part_dist_hdf5
!***********************************************************************
    subroutine get_part_dist_hdf5(name, nv)
!
      character (len=*), intent(in) :: name
      integer, intent(out) :: nv
!
      call fatal_error ('get_part_dist_hdf5', 'Not yet implemented.')
!
    endsubroutine get_part_dist_hdf5
!***********************************************************************
    subroutine index_append(varname,ivar,vector,array)
!
! 14-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: varname
      integer, intent(in) :: ivar
      integer, intent(in), optional :: vector
      integer, intent(in), optional :: array
!
      integer :: pos
!
      ! omit all unused variables
      if (ivar <= 0) return
!
      ! ignore vectors because they get expanded in 'farray_index_append'
      if (present (vector) .and. .not. present (array)) return
!
      if (lroot) open(3,file=trim(datadir)//'/'//trim(index_pro), POSITION='append')
      if (present (array)) then
        ! backwards-compatibile expansion: iuud => indgen(vector)
        if (lroot) write(3,*) trim(varname)//'=indgen('//trim(itoa(array))//')*'//trim(itoa(vector))//'+'//trim(itoa(ivar))
        ! expand array: iuud => iuud#=(#-1)*vector+ivar
        do pos=1, array
          if (lroot) write(3,*) trim(varname)//trim(itoa(pos))//'='//trim(itoa((pos-1)*vector+ivar))
          call index_register (trim(varname)//trim(itoa(pos)), (pos-1)*vector+ivar)
        enddo
      else
        if (lroot) write(3,*) trim(varname)//'='//trim(itoa(ivar))
        call index_register (trim(varname), ivar)
      endif
      if (lroot) close(3)
!
    endsubroutine index_append
!***********************************************************************
    subroutine particle_index_append(label,ilabel)
!
! 22-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: ilabel
!
      if (lroot) then
        open(3,file=trim(datadir)//'/'//trim(particle_index_pro), POSITION='append')
        write(3,*) trim(label)//'='//trim(itoa(ilabel))
        close(3)
      endif
      call index_register (trim(label), ilabel, particle=.true.)
!
    endsubroutine particle_index_append
!***********************************************************************
    function index_get(ivar,particle)
!
! 17-Oct-2018/PABourdin: coded
!
      character (len=labellen) :: index_get
      integer, intent(in) :: ivar
      logical, optional, intent(in) :: particle
!
      type (element), pointer, save :: current => null()
!
      index_get = ''
      current => last
      if (loptest (particle)) current => last_particle
      do while (associated (current))
        if (current%component == ivar) then
          index_get = current%label(2:len(current%label))
          exit
        endif
        current => current%previous
      enddo
!
    endfunction index_get
!***********************************************************************
    subroutine index_register(varname,ivar,particle)
!
! 17-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: varname
      integer, intent(in) :: ivar
      logical, optional, intent(in) :: particle
!
      type (element), pointer, save :: new => null()
!
      if (.not. loptest (particle)) then
        ! ignore variables that are not written
        if ((ivar < 1) .or. (ivar > mvar+maux)) return
      endif
!
      ! ignore non-index variables
      if ((varname(1:1) /= 'i') .or. (varname(2:2) == '_')) return
!
      ! append this entry to an internal list of written HDF5 variables
      allocate (new)
      nullify (new%previous)
      if (loptest (particle)) then
        if (associated (last_particle)) new%previous => last_particle
      else
        if (associated (last)) new%previous => last
      endif
      new%label = trim(varname)
      new%component = ivar
      if (loptest (particle)) then
        last_particle => new
      else
        last => new
      endif
!
    endsubroutine index_register
!***********************************************************************
    subroutine index_reset()
!
! 14-Oct-2018/PABourdin: coded
!
      type (element), pointer, save :: current => null()
!
      if (lroot) then
        open(3,file=trim(datadir)//'/'//trim(index_pro),status='replace')
        close(3)
        open(3,file=trim(datadir)//'/'//trim(particle_index_pro),status='replace')
        close(3)
      endif
!
      do while (associated (last))
        current => last
        last => last%previous
        deallocate (current)
        nullify (current)
      enddo
!
      do while (associated (last_particle))
        current => last_particle
        last_particle => last%previous
        deallocate (current)
        nullify (current)
      enddo
!
    endsubroutine index_reset
!***********************************************************************
endmodule HDF5_IO
