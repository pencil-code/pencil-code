! $Id$
!
!  This module takes care of massive parallel HDF5 file Input/Output.
!  We use here only F95 and MPI features for HPC-friendly behaviour.
!
module HDF5_IO
!
  use Cparam
  use Cdata
  use General, only: loptest, itoa, numeric_precision, keep_compiler_quiet
  use HDF5
  use Messages, only: fatal_error, warning
  use Mpicomm, only: lroot, mpiscan_int, mpibcast_int
!
  implicit none
!
  interface input_hdf5
    module procedure input_hdf5_string
    module procedure input_hdf5_int_0D
    module procedure input_hdf5_int_1D
    module procedure input_hdf5_0D
    module procedure input_hdf5_1D
    module procedure input_hdf5_part_2D
    module procedure input_hdf5_profile_1D
    module procedure input_hdf5_2D
    module procedure input_hdf5_3D
    module procedure input_hdf5_4D
  endinterface
!
  interface output_hdf5
    module procedure output_hdf5_string
    module procedure output_hdf5_int_0D
    module procedure output_hdf5_int_1D
    module procedure output_hdf5_0D
    module procedure output_hdf5_1D
    module procedure output_hdf5_pencil_1D
    module procedure output_hdf5_part_2D
    module procedure output_hdf5_profile_1D
    module procedure output_local_hdf5_2D
    module procedure output_hdf5_slice_2D
    module procedure output_local_hdf5_3D
    module procedure output_hdf5_3D
    module procedure output_hdf5_4D
    module procedure output_hdf5_torus_rect
  endinterface
!
  interface hdf5_input_slice
    module procedure input_slice_arr
    module procedure input_slice_real_arr
    module procedure input_slice_scat
  endinterface
!
  interface output_hdf5_double
    module procedure output_hdf5_double_0D
    module procedure output_hdf5_double_1D
  endinterface
!
  include 'hdf5_io.h'
!
  private
!
  integer :: h5_err
  integer(HID_T) :: h5_file, h5_dset, h5_plist, h5_fspace, h5_mspace, h5_dspace, h5_ntype, h5_dptype, h5_group
  integer, parameter :: n_dims = 3
  integer(kind=8), dimension(n_dims+1) :: local_size, local_subsize, local_start
  integer(kind=8), dimension(n_dims+1) :: global_size, global_start
  logical :: lcollective = .false., lwrite = .false.
  character (len=fnlen) :: current
  integer :: mvar_out, maux_out
!
! Number of open/close retries if file is locked, waiting time in seconds.
!
  integer, parameter :: ntries=10, nsleep=60
!
! The name of the calling subroutine.
!
  character(LEN=2*labellen) :: scaller=''
!
  type element
    character(len=labellen) :: label
    integer :: component
    type (element), pointer :: previous
  endtype element
  type (element), pointer :: last => null(), last_particle => null(), last_pointmass => null()
!
  contains
!***********************************************************************
    subroutine initialize_hdf5(nxyz,ngrid,mvar_,maux_)
!
!  Initialize the HDF IO.
!
!  28-Oct-2016/PABoudin: coded
!
      use Cdata
      use General, only: ioptest

      integer, dimension(3), optional :: nxyz,ngrid
      integer, optional :: mvar_,maux_

      ! dimensions for local data portion without ghost layers
      if (present(nxyz)) then
        local_subsize(1:3) = nxyz
      else
        local_subsize(1:3) = (/nx,ny,nz/)
      endif
      local_subsize(4:n_dims+1) = 1
!
      ! dimensions for local data portion with ghost layers
      local_size(1:3) = local_subsize(1:3)+2*nghost
      local_size(4:n_dims+1) = 1
!
      ! starting position of this processor's data portion in the global file
      global_start(1:3) = nghost + (/ipx,ipy,ipz/)*local_subsize(1:3)
      global_start(4:n_dims+1) = 0
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
      if (present(ngrid)) then
        global_size(1:3) = ngrid+2*nghost
      else
        global_size(1:3) = (/mxgrid,mygrid,mzgrid/)
      endif
      global_size(4:n_dims+1) = 1
!
      ! include lower ghost layers on the lower edge
      ! (upper ghost cells are taken care of by the increased 'local_subsize')
      if (lfirst_proc_x) global_start(1) = global_start(1) - nghost
      if (lfirst_proc_y) global_start(2) = global_start(2) - nghost
      if (lfirst_proc_z) global_start(3) = global_start(3) - nghost
!
      mvar_out = ioptest(mvar_,mvar)
      if (present(maux_)) then
        maux_out=maux_
      elseif (lwrite_aux) then
        maux_out = maux
      else
        maux_out = 0
      endif

    endsubroutine initialize_hdf5
!***********************************************************************
    subroutine init_hdf5
!
      use Syscalls, only: sizeof_real

      ! initialize parallel HDF5 Fortran libaray
      call h5open_f (h5_err)
      call check_error (h5_err, 'initialize parallel HDF5 library', caller='init_hdf5')
      h5_dptype = H5T_NATIVE_DOUBLE
      if (sizeof_real() < 8) then
        h5_ntype = H5T_NATIVE_REAL
      else
        h5_ntype = h5_dptype
      endif
!
    endsubroutine init_hdf5
!***********************************************************************
    subroutine finalize_hdf5
!
      ! close the HDF5 library
      call h5close_f (h5_err)
      call check_error (h5_err, 'close parallel HDF5 library', caller='finalize_hdf5')
!
    endsubroutine finalize_hdf5
!***********************************************************************
    subroutine check_error(code, message, dataset, caller, lerrcont)
!
!   7-May-2019/MR: made caller optional, added module variable scaller,
!                  so caller needs to be set only once in a subroutine
!
      use General, only: loptest

      integer, intent(in) :: code
      character (len=*), intent(in) :: message
      character (len=*), optional, intent(in) :: dataset, caller
      logical, optional, intent(INOUT) :: lerrcont
!
      if (present(caller)) scaller=caller

      ! check for an HDF5 error
      if (code /= 0) then
        if (loptest(lerrcont)) then
          if (present (dataset)) then
            call warning (scaller, message//' '//"'"//trim (dataset)//"'"//' in "'//trim (current)//'"')
          else
            call warning (scaller, message)
          endif
        else
          if (present (dataset)) then
            call fatal_error (scaller, message//' '//"'"//trim (dataset)//"'"//' in "'//trim (current)//'"', .true.)
          else
            call fatal_error (scaller, message, .true.)
          endif
        endif
      else
        if (loptest(lerrcont)) lerrcont=.false.
      endif
!
    endsubroutine check_error
!***********************************************************************
    subroutine file_open_hdf5(file, truncate, global, read_only, write, comm)
!
!   7-May-2019/MR: added optional par comm for use in h5pset_fapl_mpio_f (default: MPI_COMM_WORLD)
!
      use General, only: loptest, ioptest
      use Mpicomm, only: MPI_COMM_WORLD, MPI_INFO_NULL
!
      character (len=*), intent(inout) :: file
      logical, optional, intent(in) :: truncate
      logical, optional, intent(in) :: global
      logical, optional, intent(in) :: read_only
      logical, optional, intent(in) :: write
      integer, optional, intent(in) :: comm
!
      logical :: ltrunc, lread_only
      integer :: h5_read_mode, pos
      integer :: i
!
      if (lcollective .or. lwrite) call file_close_hdf5
!
      lread_only = loptest (read_only)
      h5_read_mode = H5F_ACC_RDWR_F
      if (lread_only) h5_read_mode = H5F_ACC_RDONLY_F
!
      ltrunc = loptest (truncate, .true.)
      if (lread_only) ltrunc = .false.
!
      lcollective = loptest (global, .true.)
      lwrite = loptest (write, lroot)
!
      pos = index (file, '.dat.h5')
      if (pos > 1) file = file(1:pos-1)//'.h5'
      current = trim (file)
!
      if (lcollective) then
        ! setup file access property list
        call h5pcreate_f (H5P_FILE_ACCESS_F, h5_plist, h5_err)
        call check_error (h5_err, 'create global file access property list', caller='file_open_hdf5')
        call h5pset_fapl_mpio_f (h5_plist, ioptest(comm,MPI_COMM_WORLD), MPI_INFO_NULL, h5_err)
        call check_error (h5_err, 'modify global file access property list')

        if (ltrunc) then
          ! create empty (or truncated) HDF5 file
          call h5fcreate_f (trim (file), H5F_ACC_TRUNC_F, h5_file, h5_err, access_prp=h5_plist)
          call check_error (h5_err, 'create global file "'//trim (file)//'"')
        else
          ! open existing HDF5 file
          i=0; h5_err=1
          do while (h5_err/=0.and.i<ntries)
            call h5fopen_f (trim (file), h5_read_mode, h5_file, h5_err, access_prp=h5_plist)
            i=i+1
            if (h5_err/=0) call sleep(nsleep)
          enddo
          call check_error (h5_err, 'open global file "'//trim (file)//'"')
        endif
!
        call h5pclose_f (h5_plist, h5_err)
        call check_error (h5_err, 'close global file access property list')
      elseif (lwrite) then
        if (ltrunc) then
          call h5fcreate_f (trim (file), H5F_ACC_TRUNC_F, h5_file, h5_err)
          call check_error (h5_err, 'create local file "'//trim (file)//'"', caller='file_open_hdf5')
        else

          i=0; h5_err=1
          do while (h5_err/=0.and.i<ntries)
            call h5fopen_f (trim (file), h5_read_mode, h5_file, h5_err)
            i=i+1
            if (h5_err/=0) call sleep(nsleep)
          enddo
          call check_error (h5_err, 'open local file "'//trim (file)//'"', caller='file_open_hdf5')
        endif
      endif
!
    endsubroutine file_open_hdf5
!***********************************************************************
    subroutine file_close_hdf5

      integer :: i
!
      if (.not. (lcollective .or. lwrite)) return

      i=0; h5_err=1

      do while (h5_err/=0.and.i<ntries)
        call h5fclose_f (h5_file, h5_err)
        i=i+1
        if (h5_err/=0) call sleep(nsleep)
      enddo
!
      call check_error (h5_err, 'close file "'//trim (current)//'"',caller='file_close_hdf5')
!
      current = repeat (' ', fnlen)
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
      if (exists_in_hdf5 (trim (name))) return
!
      call h5gcreate_f (h5_file, trim (name), h5_group, h5_err)
      call check_error (h5_err, 'create group', name, caller='create_group_hdf5')
      call h5gclose_f (h5_group, h5_err)
      call check_error (h5_err, 'close group', name)
!
    endsubroutine create_group_hdf5
!***********************************************************************
    logical function exists_in_hdf5(name)
!
      character (len=*), intent(in) :: name
!
      exists_in_hdf5 = .false.
      if (.not. (lcollective .or. lwrite)) return
!
      call h5lexists_f(h5_file, trim (name), exists_in_hdf5, h5_err)
      if (h5_err /= 0) exists_in_hdf5 = .false.
!
    endfunction exists_in_hdf5
!***********************************************************************
    subroutine input_hdf5_string(name, data)
!
      character (len=*), intent(in) :: name
      character (len=*), intent(out) :: data

      call fatal_error('input_hdf5_string','not yet implemented')
      call keep_compiler_quiet(name)
      data=''

    endsubroutine input_hdf5_string
!***********************************************************************
    subroutine input_hdf5_int_0D(name, data, lerrcont)
!
      character (len=*), intent(in) :: name
      integer, intent(out) :: data
      logical, optional, intent(INOUT) :: lerrcont
!
      integer, dimension(1) :: read
!
      call input_hdf5_int_1D (name, read, 1, lerrcont)
      data = read(1)
!
    endsubroutine input_hdf5_int_0D
!***********************************************************************
    subroutine input_local_hdf5_int_1D(name, data, nv, lerrcont)
!
!  Read HDF5 dataset as scalar or array.
!
!  05-Jun-2017/Fred: coded based on input_hdf5_1D
!
      use General, only: loptest

      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      integer, dimension (nv), intent(out) :: data
      logical, optional, intent(INOUT) :: lerrcont
!
      integer(HSIZE_T), dimension(1) :: size
!
      if (lcollective) call check_error (1, 'local input requires local file', caller='input_local_hdf5_int_1D')
      if (.not. lwrite) return
!
      size = (/ nv /)
!
      ! open dataset
      call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
      call check_error (h5_err, 'open dataset', name, caller='input_local_hdf5_int_1D', lerrcont=lerrcont)
      if (loptest(lerrcont)) return

      ! read dataset
      call h5dread_f (h5_dset, H5T_NATIVE_INTEGER, data, size, h5_err)
      call check_error (h5_err, 'read data', name)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
!
    endsubroutine input_local_hdf5_int_1D
!***********************************************************************
    subroutine input_hdf5_int_1D(name, data, nv, same_size, lerrcont)
!
!  Read HDF5 dataset as scalar or array.
!
!  24-Oct-2018/PABourdin: coded
!
      use General, only: loptest

      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      integer, dimension (nv), intent(out) :: data
      logical, optional, intent(in) :: same_size
      logical, optional, intent(INOUT) :: lerrcont
!
      integer :: total, offset, last
      integer(kind=8), dimension (1) :: local_size_1D, local_subsize_1D, local_start_1D
      integer(kind=8), dimension (1) :: global_size_1D, global_start_1D
      integer(kind=8), dimension (1) :: h5_stride, h5_count
!
      if (.not. lcollective) then
        call input_local_hdf5_int_1D(name, data, nv,lerrcont)
        return
      endif
!
      if (loptest(same_size)) then
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
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (1, local_size_1D, h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name, caller='input_hdf5_int_1D')
!
      ! open the dataset
      call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
      call check_error (h5_err, 'open dataset', name, lerrcont=lerrcont)
      if (loptest(lerrcont)) return
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      call check_error (h5_err, 'select hyperslab within memory', name)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively read the data
      call h5dread_f (h5_dset, H5T_NATIVE_INTEGER, data, &
          global_size_1D, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'read dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
!
    endsubroutine input_hdf5_int_1D
!***********************************************************************
    subroutine input_hdf5_0D(name, data, lerrcont)
!
      character (len=*), intent(in) :: name
      real, intent(out) :: data
      logical, optional, intent(INOUT) :: lerrcont
!
      real, dimension(1) :: input
!
      call input_hdf5_1D (name, input, 1, lerrcont=lerrcont)
      data = input(1)
!
    endsubroutine input_hdf5_0D
!***********************************************************************
    subroutine input_local_hdf5_1D(name, data, nv, lerrcont)
!
!  Read HDF5 dataset as scalar or array.
!
!  26-Oct-2016/PABourdin: coded
!
      use General, only: loptest

      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (nv), intent(out) :: data
      logical, optional, intent(INOUT) :: lerrcont
!
      integer(HSIZE_T), dimension(1) :: size
!
      if (lcollective) call check_error (1, 'local input requires local file')
      if (.not. lwrite) return
!
      size = (/ nv /)
!
      ! open dataset
      call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
      call check_error (h5_err, 'open dataset', name, caller='input_local_hdf5_1D',lerrcont=lerrcont)
      if (loptest(lerrcont)) return

      ! read dataset
      call h5dread_f (h5_dset, h5_ntype, data, size, h5_err)
      call check_error (h5_err, 'read data', name)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
!
    endsubroutine input_local_hdf5_1D
!***********************************************************************
    subroutine input_hdf5_1D(name, data, nv, same_size, lerrcont)
!
!  Read HDF5 dataset as scalar or array.
!
!  24-Oct-2016/PABourdin: coded
!
      use General, only: loptest

      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (nv), intent(out) :: data
      logical, optional, intent(in) :: same_size
      logical, optional, intent(INOUT) :: lerrcont
!
      integer :: total, offset, last
      integer(kind=8), dimension (1) :: local_size_1D, local_subsize_1D, local_start_1D
      integer(kind=8), dimension (1) :: global_size_1D, global_start_1D
      integer(kind=8), dimension (1) :: h5_stride, h5_count
!
      if (.not. lcollective) then
        call input_local_hdf5_1D(name, data, nv, lerrcont)
        return
      endif
!
      if (loptest(same_size)) then
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
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (1, local_size_1D, h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name, caller='input_hdf5_1D')
!
      ! open the dataset
      call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
      call check_error (h5_err, 'open dataset', name, lerrcont=lerrcont)
      if (loptest(lerrcont)) return
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      call check_error (h5_err, 'select hyperslab within memory', name)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively read the data
      call h5dread_f (h5_dset, h5_ntype, data, &
          global_size_1D, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'read dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
!
    endsubroutine input_hdf5_1D
!***********************************************************************
    subroutine input_hdf5_part_2D(name, data, mv, nc, nv, lerrcont)
!
!  Read HDF5 particle dataset into a distributed array.
!
!  24-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: mv, nc
      real, dimension (mv,nc), intent(out) :: data
      integer, intent(in) :: nv
      logical, optional, intent(INOUT) :: lerrcont
!
      integer :: pos
      character (len=labellen) :: label
!
      if (.not. lcollective) &
        call check_error (1, 'particle input requires a global file', name, caller='input_hdf5_part_2D')
!
      ! read components into particle data array
      do pos=1, nc
        if (name == 'fp') then
          label = 'part/'//trim(index_get(pos, particle=.true.))
        else
          label = trim(name)
          if (nc >= 2) label = trim(label)//'_'//trim(itoa(pos))
        endif
        call input_hdf5_1D (label, data(1:nv,pos), nv, lerrcont=lerrcont)
      enddo
!
    endsubroutine input_hdf5_part_2D
!***********************************************************************
    subroutine input_hdf5_profile_1D(name, data, ldim, gdim, np1, np2)
!
!  Write HDF5 dataset from a 1D profile.
!
!  08-Nov-2018/PABourdin: adapted from output_hdf5_slice_2D
!
      character (len=*), intent(in) :: name
      real, dimension (:) :: data
      integer, intent(in) :: ldim, gdim, np1, np2
!
      integer(kind=8), dimension (1) :: h5_stride, h5_count, loc_dim, glob_dim, loc_start, glob_start, loc_subdim
!
      if (.not. lcollective) &
        call check_error (1, '1D profile input requires global file', name, caller='input_hdf5_part_2D')
!
      loc_dim(1) = ldim
      glob_dim(1) = gdim
      loc_start(1) = 0
      glob_start(1) = np1 - 1
      loc_subdim(1) = np2 - np1 + 1
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (1, loc_dim, h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name, caller='input_hdf5_part_2D')
!
      ! open dataset
      call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
      call check_error (h5_err, 'open dataset', name)
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, glob_start, h5_count, h5_err, h5_stride, loc_subdim)
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, loc_start, h5_count, h5_err, h5_stride, loc_subdim)
      call check_error (h5_err, 'select hyperslab within memory', name)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively read the data
      call h5dread_f (h5_dset, h5_ntype, data, &
          glob_dim, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'read dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
!
    endsubroutine input_hdf5_profile_1D
!***********************************************************************
    subroutine input_hdf5_2D(name, gdims, iprocs, data, lerrcont)
!
!  Read HDF5 dataset from a distributed 2D array.
!
!  26-Oct-2016/MR: coded
!
      use General, only: loptest

      character (len=*),     intent(in) :: name
      integer, dimension(2), intent(in) :: gdims, iprocs 
      real, dimension (:,:), intent(out):: data
      logical, optional, intent(INOUT) :: lerrcont
!
      integer(kind=8), dimension(2), parameter :: h5_stride=1, h5_count=1
      integer(kind=8), dimension(2) :: ldims, i8dum
!
      ! define 'memory-space' to indicate the local data portion in memory
      ldims=(/size(data,1),size(data,2)/)
      call h5screate_simple_f (n, ldims, h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name, caller='input_hdf5_2D')
!
      ! open the dataset
      call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
      call check_error (h5_err, 'open dataset', name, lerrcont=lerrcont)
      if (loptest(lerrcont)) return
!
      ! define local 'hyper-slab' in the global file
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, iprocs*ldims, h5_count, h5_err, h5_stride, ldims)
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! define local 'hyper-slab' portion in memory
      i8dum=(/0,0/)
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, i8dum, h5_count, h5_err, h5_stride, ldims)
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively read the data
      i8dum=gdims
      call h5dread_f (h5_dset, h5_ntype, data, &
          i8dum, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'read dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
!
    endsubroutine input_hdf5_2D
!***********************************************************************
    subroutine input_hdf5_3D(name, data, lerrcont)
!
!  Read HDF5 dataset from a distributed 3D array.
!
!  26-Oct-2016/PABourdin: coded
!
      use General, only: loptest

      character (len=*), intent(in) :: name
      real, dimension (:,:,:), intent(out) :: data
      logical, optional, intent(INOUT) :: lerrcont
!
      integer(kind=8), dimension (n_dims) :: h5_stride, h5_count
      integer, parameter :: n = n_dims
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (n, local_size(1:n), h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name, caller='input_hdf5_3D')
!
      ! open the dataset
      call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
      call check_error (h5_err, 'open dataset', name, lerrcont=lerrcont)
      if (loptest(lerrcont)) return
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start(1:n), h5_count, h5_err, h5_stride, local_subsize(1:n))
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start(1:n), h5_count, h5_err, h5_stride, local_subsize(1:n))
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively read the data
      call h5dread_f (h5_dset, h5_ntype, data, &
          global_size, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'read dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
!
    endsubroutine input_hdf5_3D
!***********************************************************************
    subroutine input_hdf5_4D(name, data, nv, lerrcont)
!
!  Read HDF5 dataset from a distributed 4D array.
!
!  26-Oct-2016/PABourdin: coded
!
      use General, only: loptest

      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (:,:,:,:), intent(out) :: data
      logical, optional, intent(INOUT) :: lerrcont
!
      integer(kind=8), dimension (n_dims+1) :: h5_stride, h5_count
!
      ! read other 4D array
      global_size(n_dims+1) = nv
      local_size(n_dims+1) = nv
      local_subsize(n_dims+1) = nv
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (n_dims+1, local_size, h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name, caller='input_hdf5_4D')
!
      ! open the dataset
      call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
      call check_error (h5_err, 'open dataset', name, lerrcont=lerrcont)
      if (loptest(lerrcont)) return
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start, h5_count, h5_err, h5_stride, local_subsize)
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start, h5_count, h5_err, h5_stride, local_subsize)
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively read the data
      call h5dread_f (h5_dset, h5_ntype, data, &
          global_size, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'read dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
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
      character (len=len(data)+1), dimension(1) :: str_data
      integer(SIZE_T), dimension(1) :: str_len
!
      if (lcollective) call check_error (1, 'string output requires local file', name, caller='output_hdf5_string')
      if (.not. lwrite) return
!
      str_len(1) = len_trim (data)
      size(1) = str_len(1)
      size(2) = 1
      str_data(1) = data//char(0)
!
      ! create data space
      call H5Tcopy_f (H5T_STRING, h5_strtype, h5_err)
      call check_error (h5_err, 'copy string data space type', name, caller='output_hdf5_string')
      call H5Tset_strpad_f (h5_strtype, H5T_STR_NULLPAD_F, h5_err)
      call check_error (h5_err, 'modify string data space type', name)
      call h5screate_simple_f (1, size(2), h5_dspace, h5_err)
      call check_error (h5_err, 'create string data space', name)
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open string dataset', name)
      else
        ! create dataset
        call h5dcreate_f (h5_file, trim (name), h5_strtype, h5_dspace, h5_dset, h5_err)
        call check_error (h5_err, 'create string dataset', name)
      endif
      ! write dataset
      call h5dwrite_vl_f (h5_dset, h5_strtype, str_data, size, str_len, h5_err, h5_dspace)
      call check_error (h5_err, 'write string data', name)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close string dataset', name)
      call h5sclose_f (h5_dspace, h5_err)
      call check_error (h5_err, 'close string data space', name)
!
    endsubroutine output_hdf5_string
!***********************************************************************
    subroutine output_hdf5_torus_rect(name, data)
!
!  Outputs (potentially) varying parameters of rectangular toroid (persistent data).
!
!  16-May-2020/MR: coded
!
      use Geometrical_types
      use Iso_c_binding

      character (len=*), intent(in) :: name
      type(torus_rect), intent(in) :: data
      integer(KIND=ikind8) :: ptr
      integer(SIZE_T) :: offset

!
      integer(HID_T) :: h5_torustype, h5_vec3type
      integer(HSIZE_T), dimension(1) :: size

      real :: dummy

      return !  because of problem with compound data type
      ! create data type
      !call h5tcreate_f(H5T_COMPOUND_F, 8*sizeof(dummy), h5_torustype, h5_err)
      call check_error (h5_err, 'create torus data type', name)

      size=(/3/)
      !call h5tarray_create_f(h5_ntype, 1, size, h5_vec3type, h5_err)
      offset=OFFSETOF(data,data%center(1))
      !call h5tinsert_f(h5_torustype,"center",offset,h5_vec3type,h5_err)
      offset=OFFSETOF(data,data%th)
      !call h5tinsert_f(h5_torustype,"th",offset,h5_ntype,h5_err)
      offset=OFFSETOF(data,data%ph)
      !call h5tinsert_f(h5_torustype,"ph",offset,h5_ntype,h5_err)
      offset=OFFSETOF(data,data%r_in)
      !call h5tinsert_f(h5_torustype,"r_in",offset,h5_ntype,h5_err)
      offset=OFFSETOF(data,data%thick)
      !call h5tinsert_f(h5_torustype,"thick",offset,h5_ntype,h5_err)
      offset=OFFSETOF(data,data%height)
      !call h5tinsert_f(h5_torustype,"height",offset,h5_ntype,h5_err)
      call check_error (h5_err, 'populate torus data type', name)

      ! create data space
      size = (/ 1 /)
      !call h5screate_simple_f (1, size, h5_dspace, h5_err)
      call check_error (h5_err, 'create torus data space', name)

      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open torus dataset', name)
      else
        ! create dataset
        !call h5dcreate_f (h5_file, trim (name), h5_torustype, h5_dspace, h5_dset, h5_err)
        call check_error (h5_err, 'create torus dataset', name)
      endif

      ! [PAB] a pointer is not a valid argument for 'h5dwrite_f':
      !ptr = C_LOC(data)
      !call h5dwrite_f(h5_dset, h5_torustype, data(1), size, h5_err) ! was 'ptr' before
      !call check_error (h5_err, 'write torus dataset', name)
      ! [PAB] This will not work like this, sorry, we should do this in an easier way.
      ! Deactivating offending 'h5dwrite_f' line now, so that the autotests work again.

      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close torus dataset', name)
      call h5sclose_f (h5_dspace, h5_err)
      call check_error (h5_err, 'close torus data space', name)
      call h5tclose_f (h5_torustype, h5_err)
      call check_error (h5_err, 'close torus data type', name)

      contains
!----------------------------------------------------------------------
    function offsetof(base,comp) result(offset)

    use Geometrical_types, only: torus_rect

    integer(HSIZE_T) :: offset
    type(torus_rect) :: base
    real :: comp

    offset = loc(comp)-loc(base)

    endfunction offsetof

    endsubroutine output_hdf5_torus_rect
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
    subroutine output_local_hdf5_int_1D(name, data, nv)
!
!  Write HDF5 dataset as scalar from one or all processor.
!
!  23-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      integer, dimension(nv), intent(in) :: data
!
      integer(kind=8), dimension(1) :: size
!
      if (lcollective) call check_error (1, 'local output requires local file', caller='output_local_hdf5_int_1D')
      if (.not. lwrite) return
!
      size = (/ nv /)
!
      ! create data space
      call h5screate_simple_f (1, size, h5_dspace, h5_err)
      call check_error (h5_err, 'create integer data space', name, caller='output_local_hdf5_int_1D')
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open integer dataset', name)
      else
        ! create dataset
        call h5dcreate_f (h5_file, trim (name), H5T_NATIVE_INTEGER, h5_dspace, h5_dset, h5_err)
        call check_error (h5_err, 'create integer dataset', name)
      endif
      ! write dataset
      call h5dwrite_f (h5_dset, H5T_NATIVE_INTEGER, data, size, h5_err)
      call check_error (h5_err, 'write integer data', name)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close integer dataset', name)
      call h5sclose_f (h5_dspace, h5_err)
      call check_error (h5_err, 'close integer data space', name)
!
    endsubroutine output_local_hdf5_int_1D
!***********************************************************************
    subroutine output_hdf5_int_1D(name, data, nv, same_size)
!
!  Write HDF5 dataset as array from one or all processors.
!
!  24-Oct-2018/PABourdin: coded
!
      use General, only: loptest

      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      integer, dimension (nv), intent(in) :: data
      logical, optional, intent(in) :: same_size
!
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
      if (loptest(same_size)) then
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
      call check_error (h5_err, 'create global file space', name, caller='output_hdf5_int_1D')
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (1, local_size_1D, h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name)
!
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open integer dataset', name)
      else
        ! create the dataset
        call h5dcreate_f (h5_file, trim (name), H5T_NATIVE_INTEGER, h5_fspace, h5_dset, h5_err)
        call check_error (h5_err, 'create dataset', name)
      endif
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close global file space', name)
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      call check_error (h5_err, 'select hyperslab within memory', name)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively write the data
      call h5dwrite_f (h5_dset, H5T_NATIVE_INTEGER, data, &
          global_size_1D, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'write dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
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
      if (lcollective) call check_error (1, 'local output requires local file', caller='output_local_hdf5_1D')
      if (.not. lwrite) return
!
      size = (/ nv /)
!
      ! create data space
      if (nv <= 1) then
        call h5screate_f (H5S_SCALAR_F, h5_dspace, h5_err)
        call check_error (h5_err, 'create scalar data space', name, caller='output_local_hdf5_1D')
        call h5sset_extent_simple_f (h5_dspace, 0, size(1), size(1), h5_err)
      else
        call h5screate_f (H5S_SIMPLE_F, h5_dspace, h5_err)
        call check_error (h5_err, 'create simple data space', name, caller='output_local_hdf5_1D')
        call h5sset_extent_simple_f (h5_dspace, 1, size, size, h5_err)
      endif
      call check_error (h5_err, 'set data space extent', name)
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open dataset', name)
      else
        ! create dataset
        call h5dcreate_f (h5_file, trim (name), h5_ntype, h5_dspace, h5_dset, h5_err)
        call check_error (h5_err, 'create dataset', name)
      endif
      ! write dataset
      if (nv <= 1) then
        call h5dwrite_f (h5_dset, h5_ntype, data(1), size, h5_err)
      else
        call h5dwrite_f (h5_dset, h5_ntype, data, size, h5_err)
      endif
      call check_error (h5_err, 'write data', name)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5sclose_f (h5_dspace, h5_err)
      call check_error (h5_err, 'close data space', name)
!
    endsubroutine output_local_hdf5_1D
!***********************************************************************
    subroutine output_hdf5_1D(name, data, nv, same_size)
!
!  Write HDF5 dataset as scalar or array.
!
!  24-Oct-2016/PABourdin: coded
!
      use General, only: loptest

      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (nv), intent(in) :: data
      logical, optional, intent(in) :: same_size
!
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
      if (loptest(same_size)) then
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
      call check_error (h5_err, 'create global file space', name, caller='output_hdf5_1D')
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (1, local_size_1D, h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name)
!
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open dataset', name)
      else
        ! create the dataset
        call h5dcreate_f (h5_file, trim (name), h5_ntype, h5_fspace, h5_dset, h5_err)
        call check_error (h5_err, 'create dataset', name)
      endif
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close global file space', name)
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start_1D, h5_count, h5_err, h5_stride, local_subsize_1D)
      call check_error (h5_err, 'select hyperslab within memory', name)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively write the data
      call h5dwrite_f (h5_dset, h5_ntype, data, &
          global_size_1D, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'write dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
!
    endsubroutine output_hdf5_1D
!***********************************************************************
    subroutine output_hdf5_pencil_1D(name, data, py, pz)
!
!  Write HDF5 dataset from penciled data.
!
!  24-Jun-2019/PABourdin: coded
!
      character (len=*), intent(in) :: name
      real, dimension (nx), intent(in) :: data
      integer, intent(in) :: py, pz
!
      integer, parameter :: n = 3
      integer(kind=8), dimension (n) :: h5_stride, h5_count, loc_dim, glob_dim, loc_start, glob_start
!
      if (.not. lcollective) &
        call check_error (1, '1D pencil output requires global file', name, caller='output_hdf5_pencil_1D')
!
      loc_dim = (/nx,1,1/)
      glob_dim = global_size(1:3) - 2*nghost
      loc_start = 0
      glob_start = global_start(1:3) - nghost + (/0,py,pz/)
!
      ! define 'file-space' to indicate the data portion in the global file
      call h5screate_simple_f (n, glob_dim, h5_fspace, h5_err)
      call check_error (h5_err, 'create global file space', name, caller='output_hdf5_pencil_1D')
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (n, loc_dim, h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name)
!
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open dataset', name)
      else
        ! create the dataset
        call h5dcreate_f (h5_file, trim (name), h5_ntype, h5_fspace, h5_dset, h5_err)
        call check_error (h5_err, 'create dataset', name)
      endif
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close global file space', name)
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, glob_start, h5_count, h5_err, h5_stride, loc_dim)
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, loc_start, h5_count, h5_err, h5_stride, loc_dim)
      call check_error (h5_err, 'select hyperslab within memory', name)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively write the data
      call h5dwrite_f (h5_dset, h5_ntype, data, &
          glob_dim, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'write dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
!
    endsubroutine output_hdf5_pencil_1D
!***********************************************************************
    subroutine output_hdf5_part_2D(name, data, mv, nc, nv)
!
!  Write HDF5 dataset from a distributed particle array.
!
!  22-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: mv, nc
      real, dimension (mv,nc), intent(in) :: data
      integer, intent(in) :: nv
!
      integer :: pos
      character (len=labellen) :: label
!
      if (.not. lcollective) &
        call check_error (1, 'particle output requires a global file', name, caller='output_hdf5_part_2D')
!
      ! write components of particle data array
      do pos=1, nc
        if (name == 'fp') then
          label = 'part/'//trim(index_get(pos, particle=.true.))
        else
          label = trim(name)
          if (nc >= 2) label = trim(label)//'_'//trim(itoa(pos))
        endif
        call output_hdf5_1D (label, data(1:nv,pos), nv)
      enddo
!
    endsubroutine output_hdf5_part_2D
!***********************************************************************
    subroutine output_hdf5_profile_1D(name, data, ldim, gdim, ip, np1, np2, ng, lhas_data)
!
!  Write HDF5 dataset from a 1D profile.
!
!  08-Nov-2018/PABourdin: adapted from output_hdf5_slice_2D
!
      character (len=*), intent(in) :: name
      real, dimension (:), intent(in) :: data
      integer, intent(in) :: ldim, gdim, ip, np1, np2, ng
      logical, intent(in) :: lhas_data
!
      integer(kind=8), dimension (1) :: h5_stride, h5_count, loc_dim, glob_dim, loc_start, glob_start, loc_subdim
!
      if (.not. lcollective) &
        call check_error (1, '1D profile output requires global file', name, caller='output_hdf5_profile_1D')
!
      loc_dim(1) = ldim
      glob_dim(1) = gdim
      loc_start(1) = np1 - 1
      glob_start(1) = ip * (ldim - 2*ng) + loc_start(1)
      loc_subdim(1) = np2 - np1 + 1
!
      ! define 'file-space' to indicate the data portion in the global file
      call h5screate_simple_f (1, glob_dim, h5_fspace, h5_err)
      call check_error (h5_err, 'create global file space', name, caller='output_hdf5_profile_1D')
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (1, loc_dim, h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name)
!
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open dataset', name)
      else
        ! create the dataset
        call h5dcreate_f (h5_file, trim (name), h5_ntype, h5_fspace, h5_dset, h5_err)
        call check_error (h5_err, 'create dataset', name)
      endif
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close global file space', name)
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, glob_start, h5_count, h5_err, h5_stride, loc_subdim)
      call check_error (h5_err, 'select hyperslab within file', name)
      if (.not. lhas_data) then
        call h5sselect_none_f (h5_fspace, h5_err)
        call check_error (h5_err, 'set empty hyperslab within file', name)
      endif
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, loc_start, h5_count, h5_err, h5_stride, loc_subdim)
      call check_error (h5_err, 'select hyperslab within memory', name)
      if (.not. lhas_data) then
        call h5sselect_none_f (h5_mspace, h5_err)
        call check_error (h5_err, 'set empty hyperslab within memory', name)
      endif
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively write the data
      call h5dwrite_f (h5_dset, h5_ntype, data, &
          glob_dim, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'write dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
!
    endsubroutine output_hdf5_profile_1D
!***********************************************************************
    subroutine output_local_hdf5_2D(name, data, dim1, dim2)
!
!  Write HDF5 dataset from a local 2D array.
!
!  14-Nov-2018/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: dim1, dim2
      real, dimension (dim1,dim2), intent(in) :: data
!
      integer(kind=8), dimension(2) :: size
!
      if (lcollective) call check_error (1, 'local 2D output requires local file', caller='output_local_hdf5_2D')
      if (.not. lwrite) return
!
      size = (/ dim1, dim2 /)
!
      ! create data space
      call h5screate_f (H5S_SIMPLE_F, h5_dspace, h5_err)
      call check_error (h5_err, 'create simple data space', name, caller='output_local_hdf5_2D')
      call h5sset_extent_simple_f (h5_dspace, 2, size, size, h5_err)
      call check_error (h5_err, 'set data space extent', name)
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open dataset', name)
      else
        ! create dataset
        call h5dcreate_f (h5_file, trim (name), h5_ntype, h5_dspace, h5_dset, h5_err)
        call check_error (h5_err, 'create dataset', name)
      endif
      ! write dataset
      call h5dwrite_f (h5_dset, h5_ntype, data, size, h5_err)
      call check_error (h5_err, 'write data', name)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5sclose_f (h5_dspace, h5_err)
      call check_error (h5_err, 'close data space', name)
!
    endsubroutine output_local_hdf5_2D
!***********************************************************************
    subroutine output_hdf5_slice_2D(name, data, ldim1, ldim2, gdim1, gdim2, ip1, ip2, has_data)
!
!  Write HDF5 dataset from a 2D slice.
!
!  29-Oct-2018/PABourdin: coded
!   7-May-2019/MR: made has_data optional (default: .true.)
!
      use General, only: loptest

      character (len=*), intent(in) :: name
      real, dimension (:,:), pointer :: data
      integer, intent(in) :: ldim1, ldim2, gdim1, gdim2, ip1, ip2
      logical, optional, intent(in) :: has_data
!
      integer(kind=8), dimension (2) :: h5_stride, h5_count, loc_dim, glob_dim, loc_start, glob_start
      logical :: lhas_data
!
      if (.not. lcollective) &
        call check_error (1, '2D slice output requires global file', name, caller='output_hdf5_slice_2D')
!
      lhas_data=loptest(has_data,.true.)
      loc_dim(1) = ldim1
      loc_dim(2) = ldim2
      glob_dim(1) = gdim1
      glob_dim(2) = gdim2
      loc_start(1) = 0
      loc_start(2) = 0
      glob_start(1) = ip1 * ldim1
      glob_start(2) = ip2 * ldim2
!
      ! define 'file-space' to indicate the data portion in the global file
      call h5screate_simple_f (2, glob_dim, h5_fspace, h5_err)
      call check_error (h5_err, 'create global file space', name, caller='output_hdf5_slice_2D')
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (2, loc_dim, h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name)
!
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open dataset', name)
      else
        ! create the dataset
        call h5dcreate_f (h5_file, trim (name), h5_ntype, h5_fspace, h5_dset, h5_err)
        call check_error (h5_err, 'create dataset', name)
      endif
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close global file space', name)
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, glob_start, h5_count, h5_err, h5_stride, loc_dim)
      call check_error (h5_err, 'select hyperslab within file', name)
      if (.not. lhas_data) then
        call h5sselect_none_f (h5_fspace, h5_err)
        call check_error (h5_err, 'output_hdf5_slice_2D', 'set empty hyperslab within file', name)
      endif
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, loc_start, h5_count, h5_err, h5_stride, loc_dim)
      call check_error (h5_err, 'select hyperslab within memory', name)
      if (.not. lhas_data) then
        call h5sselect_none_f (h5_mspace, h5_err)
        call check_error (h5_err, 'output_hdf5_slice_2D', 'set empty hyperslab within memory', name)
      endif
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively write the data
      if (lhas_data) then
        call h5dwrite_f (h5_dset, h5_ntype, data, &
             glob_dim, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      else
        call h5dwrite_f (h5_dset, h5_ntype, 0, &
            glob_dim, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      endif
      call check_error (h5_err, 'write dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
!
    endsubroutine output_hdf5_slice_2D
!***********************************************************************
    subroutine output_local_hdf5_3D(name, data, dim1, dim2, dim3)
!
!  Write HDF5 dataset from a local 3D array.
!
!  26-Nov-2018/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: dim1, dim2, dim3
      real, dimension (dim1,dim2,dim3), intent(in) :: data
!
      integer(kind=8), dimension(3) :: size
!
      if (lcollective) call check_error (1, 'local 3D output requires local file', caller='output_local_hdf5_3D')                       
      if (.not. lwrite) return
!
      size = (/ dim1, dim2, dim3 /)
!
      ! create data space
      call h5screate_f (H5S_SIMPLE_F, h5_dspace, h5_err)
      call check_error (h5_err, 'create simple data space', name, caller='output_local_hdf5_3D')
      call h5sset_extent_simple_f (h5_dspace, 3, size, size, h5_err)
      call check_error (h5_err, 'set data space extent', name)
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open dataset', name)
      else
        ! create dataset
        call h5dcreate_f (h5_file, trim (name), h5_ntype, h5_dspace, h5_dset, h5_err)
        call check_error (h5_err, 'create dataset', name)
      endif
      ! write dataset
      call h5dwrite_f (h5_dset, h5_ntype, data, size, h5_err)
      call check_error (h5_err, 'write data', name)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5sclose_f (h5_dspace, h5_err)
      call check_error (h5_err, 'close data space', name)
!
    endsubroutine output_local_hdf5_3D
!***********************************************************************
    subroutine output_hdf5_3D(name, data)
!
!  Write HDF5 dataset from a distributed 3D array.
!
!  17-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: name
      real, dimension (:,:,:), intent(in) :: data
!
      integer(kind=8), dimension (n_dims) :: h5_stride, h5_count
      integer, parameter :: n = n_dims
!
      if (.not. lcollective) &
        call check_error (1, '3D array output requires global file', name, caller='output_hdf5_3D')
!
      ! define 'file-space' to indicate the data portion in the global file
      call h5screate_simple_f (n, global_size(1:n), h5_fspace, h5_err)
      call check_error (h5_err, 'create global file space', name, caller='output_hdf5_3D')
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (n, local_size(1:n), h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name)
!
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open dataset', name)
      else
        ! create the dataset
        call h5dcreate_f (h5_file, trim (name), h5_ntype, h5_fspace, h5_dset, h5_err)
        call check_error (h5_err, 'create dataset', name)
      endif
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close global file space', name)
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start(1:n), h5_count, h5_err, h5_stride, local_subsize(1:n))
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start(1:n), h5_count, h5_err, h5_stride, local_subsize(1:n))
      call check_error (h5_err, 'select hyperslab within memory', name)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively write the data
      call h5dwrite_f (h5_dset, h5_ntype, data, &
          global_size, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'write dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
!
    endsubroutine output_hdf5_3D
!***********************************************************************
    subroutine output_hdf5_4D(name, data, nv, compress)
!
!  Write HDF5 dataset from a distributed 4D array.
!
!  26-Oct-2016/PABourdin: coded
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (:,:,:,:), intent(in) :: data
      logical, optional, intent(in) :: compress
!
      integer(kind=8), dimension (n_dims+1) :: h5_stride, h5_count
      integer(kind=8), dimension(4), parameter :: chunk_dims=(/128,128,128,128/)
!
      if (.not. lcollective) &
        call check_error (1, '4D array output requires global file', name, caller='output_hdf5_4D')
!
      ! write other 4D array
      global_size(n_dims+1) = nv
      local_size(n_dims+1) = nv
      local_subsize(n_dims+1) = nv
!
      ! define 'file-space' to indicate the data portion in the global file
      call h5screate_simple_f (n_dims+1, global_size, h5_fspace, h5_err)
      call check_error (h5_err, 'create global file space', name, caller='output_hdf5_4D')
!
      ! define 'memory-space' to indicate the local data portion in memory
      call h5screate_simple_f (n_dims+1, local_size, h5_mspace, h5_err)
      call check_error (h5_err, 'create local memory space', name)

      if (loptest(compress)) then     ! not yet tested
        call h5pset_chunk_f(h5_plist, 4, chunk_dims, h5_err)
        call h5pset_deflate_f(h5_plist, 6, h5_err)
      endif
!
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open dataset', name)
      else
        ! create the dataset
        call h5dcreate_f (h5_file, trim (name), h5_ntype, h5_fspace, h5_dset, h5_err)
        call check_error (h5_err, 'create dataset', name)
      endif
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close global file space', name)
!
      ! define local 'hyper-slab' in the global file
      h5_stride(:) = 1
      h5_count(:) = 1
      call h5dget_space_f (h5_dset, h5_fspace, h5_err)
      call check_error (h5_err, 'get dataset for file space', name)
      call h5sselect_hyperslab_f (h5_fspace, H5S_SELECT_SET_F, global_start, h5_count, h5_err, h5_stride, local_subsize)
      call check_error (h5_err, 'select hyperslab within file', name)
!
      ! define local 'hyper-slab' portion in memory
      call h5sselect_hyperslab_f (h5_mspace, H5S_SELECT_SET_F, local_start, h5_count, h5_err, h5_stride, local_subsize)
      call check_error (h5_err, 'select hyperslab within memory', name)
!
      ! prepare data transfer
      call h5pcreate_f (H5P_DATASET_XFER_F, h5_plist, h5_err)
      call check_error (h5_err, 'set data transfer properties', name)
      call h5pset_dxpl_mpio_f (h5_plist, H5FD_MPIO_COLLECTIVE_F, h5_err)
      call check_error (h5_err, 'select collective IO', name)
!
      ! collectively write the data
      call h5dwrite_f (h5_dset, h5_ntype, data, &
          global_size, h5_err, file_space_id=h5_fspace, mem_space_id=h5_mspace, xfer_prp=h5_plist)
      call check_error (h5_err, 'write dataset', name)
!
      ! close data spaces, dataset, and the property list
      call h5sclose_f (h5_fspace, h5_err)
      call check_error (h5_err, 'close file space', name)
      call h5sclose_f (h5_mspace, h5_err)
      call check_error (h5_err, 'close memory space', name)
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5pclose_f (h5_plist, h5_err)
      call check_error (h5_err, 'close parameter list', name)
!
    endsubroutine output_hdf5_4D
!***********************************************************************
    subroutine output_hdf5_double_0D(name, data)
!
      character (len=*), intent(in) :: name
      real(KIND=rkind8), intent(in) :: data
!
      call output_hdf5_double_1D (name, (/ data /), 1)
!
    endsubroutine output_hdf5_double_0D
!***********************************************************************
    subroutine output_hdf5_double_1D(name, data, nv)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real(KIND=rkind8), dimension (nv), intent(in) :: data
!
      call output_local_hdf5_double_1D (name, data, nv)
!
    endsubroutine output_hdf5_double_1D
!***********************************************************************
    subroutine output_local_hdf5_double_1D(name, data, nv)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real(KIND=rkind8), dimension (nv), intent(in) :: data
!
      integer(KIND=ikind8), dimension(1) :: size
!
      if (lcollective) call check_error (1, 'local output requires local file', caller='output_local_hdf5_double_1D')
      if (.not. lwrite) return
!
      size = (/ 1 /)
!
      ! create data space
      if (nv <= 1) then
        call h5screate_f (H5S_SCALAR_F, h5_dspace, h5_err)
        call check_error (h5_err, 'create scalar data space', name, caller='output_local_hdf5_double_1D')
        call h5sset_extent_simple_f (h5_dspace, 0, size(1), size(1), h5_err)
      else
        call h5screate_f (H5S_SIMPLE_F, h5_dspace, h5_err)
        call check_error (h5_err, 'create simple data space', name, caller='output_local_hdf5_double_1D')
        call h5sset_extent_simple_f (h5_dspace, 1, size, size, h5_err)
      endif
      call check_error (h5_err, 'set data space extent', name)
      if (exists_in_hdf5 (name)) then
        ! open dataset
        call h5dopen_f (h5_file, trim (name), h5_dset, h5_err)
        call check_error (h5_err, 'open dataset', name)
      else
        ! create dataset
        call h5dcreate_f (h5_file, trim (name), h5_dptype, h5_dspace, h5_dset, h5_err)
        call check_error (h5_err, 'create dataset', name)
      endif
      ! write dataset
      if (nv <= 1) then
        call h5dwrite_f (h5_dset, h5_dptype, data(1), size, h5_err)
      else
        call h5dwrite_f (h5_dset, h5_dptype, data, size, h5_err)
      endif
      call check_error (h5_err, 'write data', name)
      ! close dataset and data space
      call h5dclose_f (h5_dset, h5_err)
      call check_error (h5_err, 'close dataset', name)
      call h5sclose_f (h5_dspace, h5_err)
      call check_error (h5_err, 'close data space', name)
!
    endsubroutine output_local_hdf5_double_1D
!***********************************************************************
    subroutine output_dim(file, mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar_out, maux_out, mglobal)
!
!  Write dimension to file.
!
!  02-Nov-2018/PABourdin: coded
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar_out, maux_out, mglobal
!
      character (len=fnlen) :: filename
!
      filename = trim(datadir)//'/'//trim(file)//'.h5'
      call file_open_hdf5 (filename, global=.false., truncate=.true.)
      call output_hdf5 ('nx', mx_out - 2*nghost)
      call output_hdf5 ('ny', my_out - 2*nghost)
      call output_hdf5 ('nz', mz_out - 2*nghost)
      call output_hdf5 ('mx', mx_out)
      call output_hdf5 ('my', my_out)
      call output_hdf5 ('mz', mz_out)
      call output_hdf5 ('nxgrid', mxgrid_out - 2*nghost)
      call output_hdf5 ('nygrid', mygrid_out - 2*nghost)
      call output_hdf5 ('nzgrid', mzgrid_out - 2*nghost)
      call output_hdf5 ('mxgrid', mxgrid_out)
      call output_hdf5 ('mygrid', mygrid_out)
      call output_hdf5 ('mzgrid', mzgrid_out)
      call output_hdf5 ('mvar', mvar_out)
      call output_hdf5 ('maux', maux_out)
      call output_hdf5 ('mglobal', mglobal)
      call output_hdf5 ('precision', numeric_precision())
      call output_hdf5 ('nghost', nghost)
      if (lprocz_slowest) then
        call output_hdf5 ('procz_slowest', 1)
      else
        call output_hdf5 ('procz_slowest', 0)
      endif
      call output_hdf5 ('nprocx', nprocx)
      call output_hdf5 ('nprocy', nprocy)
      call output_hdf5 ('nprocz', nprocz)
      call output_hdf5 ('ncpus', ncpus)
      call file_close_hdf5
!
    endsubroutine output_dim
 !***********************************************************************
    subroutine input_dim(file, mx_in, my_in, mz_in, mvar_in, maux_in, mglobal_in, &
                         prec_in, nghost_in, nprocx_in, nprocy_in, nprocz_in, local)
!
!  Read dimensions from dim.dat (local or global).
!
      use General, only: loptest

      character (len=*), intent(in) :: file
      integer, intent(out) :: mx_in, my_in, mz_in, mvar_in, maux_in, mglobal_in
      integer, intent(out) :: nprocx_in, nprocy_in, nprocz_in, nghost_in
      character, intent(out) :: prec_in
      logical, optional :: local
!
      character (len=fnlen) :: filename
      integer :: iproc_slowest
      integer :: mxgrid_in, mygrid_in, mzgrid_in, ncpus_in
!
      filename = trim(datadir)//'/'//trim(file)//'.h5'
!
      call file_open_hdf5 (filename, global=.false., truncate=.true.)
      call input_hdf5 ('mx', mx_in)
      call input_hdf5 ('my', my_in)
      call input_hdf5 ('mz', mz_in)
      call input_hdf5 ('mxgrid', mxgrid_in)
      call input_hdf5 ('mygrid', mygrid_in)
      call input_hdf5 ('mzgrid', mzgrid_in)
      call input_hdf5 ('mvar', mvar_in)
      call input_hdf5 ('maux', maux_in)
      call input_hdf5 ('mglobal', mglobal_in)
      call input_hdf5 ('precision', prec_in)
      call input_hdf5 ('nghost', nghost_in)
      call input_hdf5 ('procz_slowest', iproc_slowest)
      call input_hdf5 ('nprocx', nprocx_in)
      call input_hdf5 ('nprocy', nprocy_in)
      call input_hdf5 ('nprocz', nprocz_in)
      call input_hdf5 ('ncpus', ncpus_in)
      call file_close_hdf5

    endsubroutine input_dim
!***********************************************************************
    subroutine wdim_default_grid(file)
!
!  Write dimension to file.
!
!  02-Nov-2018/PABourdin: redesigned
!
      character (len=*), intent(in) :: file
!
      if (file == 'dim.dat') return
      call output_dim (file, mx, my, mz, mxgrid, mygrid, mzgrid, mvar, maux, mglobal)
!
    endsubroutine wdim_default_grid
!***********************************************************************
    subroutine wdim_default(file, mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out)
!
!  Write dimension to file.
!
!  02-Nov-2018/PABourdin: redesigned
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out
!
      call output_dim (file, mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar, maux, mglobal)
!
    endsubroutine wdim_default
!***********************************************************************
    subroutine wdim(file, mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar_out, maux_out)
!
!  Write dimension to file.
!
!   8-sep-01/axel: adapted to take my_out,mz_out
!   4-oct-16/MR: added optional parameters mvar_out,maux_out
!  02-Nov-2018/PABourdin: redesigned, moved to IO modules
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar_out, maux_out
!
      call output_dim (file, mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar_out, maux_out, mglobal)
!
    endsubroutine wdim
!***********************************************************************
    subroutine input_aver_2D(filename, time, variables, data)
!       
!  Read an 2D-average file at a given time.
!
!  01-april-21/MR: coded
!
      use File_IO, only: file_exists
      use Mpicomm, only: mpibarrier

      character (len=*),                    intent(inout):: filename
      character (len=*), dimension(:),      intent(in)   :: variables
      real,                                 intent(in)   :: time
      real,              dimension(:,:,:,:),intent(out)  :: data
!            
      integer :: it, nt, comm, slice_root
      real :: tt
      character(LEN=fnlen) :: group
!
      if (lroot) then
        if (file_exists (filename)) then
          ! find last written average
          call file_open_hdf5(filename,global=.false.,read_only=.true.,write=.false.)
          if (exists_in_hdf5 ('last')) then
            call input_hdf5 ('last', nt)
          else
            call fatal_error('input_aver_2D','no "last" group in HDF5 file '//trim(filename))
          endif
          call file_close_hdf5
        else
          call fatal_error('input_aver_2D','no HDF5 file '//trim(filename))
        endif
      endif
!
      call mpibarrier(comm)
      call file_open_hdf5(filename,truncate=.false.,read_only=.true.,write=.false.,comm=comm)

      do it=1,nt
        group=itoa(it)
        call input_hdf5 (trim(group)//'/time', tt)

!          call input_hdf5 (trim(group)//'data', (/nxgrid, nygrid/), (/ipx,
!          ipy/), data(:,:,it))
      enddo
      data=0.
      call file_close_hdf5
!
    endsubroutine input_aver_2D
!***********************************************************************
    subroutine output_average_1D(path, label, nc, name, data, time, lbinary, lwrite, header)
!
!   Output 1D average to a file.
!
!   26-Nov-2018/PABourdin: coded
!
      use File_io, only: file_exists
      use General, only: itoa
!
      character (len=*), intent(in) :: path, label
      integer, intent(in) :: nc
      character (len=fmtlen), dimension(nc), intent(in) :: name
      real, dimension(:,:), intent(in) :: data
      real, intent(in) :: time
      logical, intent(in) :: lbinary, lwrite
      real, dimension(:), optional, intent(in) :: header
!
      character (len=fnlen) :: filename
      character (len=intlen) :: group
      integer :: last, ia
      logical :: lexists
!
      if (.not. lwrite .or. (nc <= 0)) return
!
      filename = trim(datadir)//'/averages/'//trim(label)//'.h5'
      lexists = file_exists (filename)
      call file_open_hdf5 (filename, global=.false., truncate=(.not. lexists))
      if (exists_in_hdf5 ('last')) then
        call input_hdf5 ('last', last)
        last = last + 1
      else
        if (present (header)) call output_hdf5 ('r', header, size(header))
        last = 0
      endif
      group = trim (itoa (last))
      call create_group_hdf5 (group)
      do ia = 1, nc
        call output_hdf5 (trim(group)//'/'//trim(name(ia)), data(:,ia), size(data,1))
      enddo
      call output_hdf5 (trim(group)//'/time', time)
      call output_hdf5 ('last', last)
      call file_close_hdf5
!
    endsubroutine output_average_1D
!***********************************************************************
    subroutine output_average_1D_chunked(path, label, nc, name, data, full, time, lbinary, lwrite, header)
!
!   Output 1D chunked average to a file.
!
!   26-Nov-2018/PABourdin: coded
!
      use File_io, only: file_exists
      use General, only: itoa
!
      character (len=*), intent(in) :: path, label
      integer, intent(in) :: nc
      character (len=fmtlen), dimension(nc), intent(in) :: name
      real, dimension(:,:,:), intent(in) :: data
      integer, intent(in) :: full
      real, intent(in) :: time
      logical, intent(in) :: lbinary, lwrite
      real, dimension(:), optional, intent(in) :: header
!
      character (len=fnlen) :: filename
      character (len=intlen) :: group
      integer :: last, ia
      logical :: lexists
!
      if (.not. lwrite .or. (nc <= 0)) return
!
      filename = trim(datadir)//'/averages/'//trim(label)//'.h5'
      lexists = file_exists (filename)
      call file_open_hdf5 (filename, global=.false., truncate=(.not. lexists))
      if (exists_in_hdf5 ('last')) then
        call input_hdf5 ('last', last)
        last = last + 1
      else
        if (present (header)) call output_hdf5 ('r', header, size(header))
        last = 0
      endif
      group = trim (itoa (last))
      call create_group_hdf5 (group)
!
      do ia = 1, nc
        call output_hdf5 (trim(group)//'/'//trim(name(ia)), reshape (data(:,:,ia), (/ full /)), full)
      enddo
      call output_hdf5 (trim(group)//'/time', time)
      call output_hdf5 ('last', last)
      call file_close_hdf5
!
    endsubroutine output_average_1D_chunked
!***********************************************************************
    subroutine output_average_phi(path, number, nr, nc, name, data, time, r, dr)
!       
!   Output phi average to a file with these records:
!   1) nr_phiavg, nz_phiavg, nvars, nprocz
!   2) t, r_phiavg, z_phiavg, dr, dz
!   3) data
!   4) len(labels),labels
!
!   27-Nov-2014/PABourdin: cleaned up code from write_phiaverages
!   25-Nov-2018/PABourdin: coded
!         
      use File_io, only: file_exists
      use General, only: itoa
!
      character (len=*), intent(in) :: path, number
      integer, intent(in) :: nr, nc
      character (len=fmtlen), dimension(nc), intent(in) :: name
      real, dimension(:,:,:,:), intent(in) :: data 
      real, intent(in) :: time
      real, dimension(nr), intent(in) :: r
      real, intent(in) :: dr 
!   
      character (len=fnlen) :: filename
      character (len=intlen) :: group
      integer :: last, ia
      logical :: lexists
      real, dimension (nr,nzgrid) :: component
!
      if (.not. lroot .or. (nc <= 0)) return
!
      filename = trim(datadir)//'/averages/phi.h5'
      lexists = file_exists (filename)
      call file_open_hdf5 (filename, global=.false., truncate=(.not. lexists))
      if (exists_in_hdf5 ('last')) then
        call input_hdf5 ('last', last)
        last = last + 1
      else
        call output_hdf5 ('r', r, nr)
        call output_hdf5 ('dr', dr)
        last = 0
      endif
      group = trim (itoa (last))
      call create_group_hdf5 (group)
      do ia = 1, nc
        component = reshape (data(:,:,:,ia), (/ nr, nzgrid /))
        call output_hdf5 (trim(group)//'/'//trim(name(ia)), component, size(data,1), nzgrid)
      enddo
      call output_hdf5 (trim(group)//'/time', time)
      call output_hdf5 ('last', last)
      call file_close_hdf5
!
    endsubroutine output_average_phi
!***********************************************************************
    subroutine trim_average(path, plane, ngrid, nname)
!       
!  Trim a 1D-average file for times past the current time.
!         
!  25-apr-16/ccyang: coded
!  23-Nov-2018/PABourdin: moved to IO module
!       
      use File_io, only: file_exists, delete_file
      use General, only: itoa
!         
      character (len=*), intent(in) :: path, plane
      integer, intent(in) :: ngrid, nname
!         
      character(len=fnlen) :: filename
      real :: time_file, t_sp
      integer :: last, pos
!
      if (.not. lroot) return 
      if ((ngrid <= 0) .or. (nname <= 0)) return 
!        
      filename = trim(datadir)//'/averages/'//trim(plane)//'.h5'
      if (.not. file_exists (filename)) return
!       
      t_sp = real (t)
      call file_open_hdf5 (filename, global=.false., truncate=.false.)
      if (exists_in_hdf5 ('last')) then
        call input_hdf5 ('last', last)
        call input_hdf5 (trim(itoa(last))//'/time', time_file)
        if (time_file > t_sp) then
          do pos = last, 0, -1 
            if (pos < last) call input_hdf5 (trim(itoa(pos))//'/time', time_file)
            if (time_file < t_sp) then
              if (pos /= last) call output_hdf5 ('last', pos)
              call file_close_hdf5
              return
            endif
          enddo
        endif
      endif
      call file_close_hdf5
      if (t_sp == 0.0) call delete_file (filename)
!
    endsubroutine trim_average
!***********************************************************************
    subroutine output_timeseries(data, data_im)
!
!  Append diagnostic data to 'time_series.h5' file.
!
!  01-Apr-2019/PABourdin: coded
!
      use File_io, only: file_exists
      use General, only: itoa
!
      real, dimension(2*nname), intent(in) :: data, data_im
!
      integer :: pos, iteration
      character (len=fmtlen) label
      character (len=fnlen) :: filename
      logical :: lexists
      integer, save :: offset = -1
!
      filename = trim(datadir)//'/time_series.h5'
      lexists = file_exists (filename)
      call file_open_hdf5 (filename, global=.false., truncate=.not. lexists)
      if (offset == -1) then
        offset = 0
        if (lexists .and. exists_in_hdf5 ('last')) call input_hdf5 ('last', offset)
      endif
      iteration = offset + it-1
      do pos = 1, nname
        label = cname(pos)
        label = label(1:min(index(label,' '), index(label,'('))-1)
        if (label == 'it') cycle
        call create_group_hdf5 (label)
        call output_hdf5 (trim (label)//'/'//itoa (iteration), data(pos))
        if ((itype_name(pos) >= ilabel_complex) .and. (cform(pos) /= '')) then
          label = trim(label)//'/imaginary_part'
          call create_group_hdf5 (label)
          call output_hdf5 (trim (label)//'/'//itoa (iteration), data_im(pos))
        endif
      enddo
      call output_hdf5 ('last', iteration)
      call output_hdf5 ('step', it1)
      call file_close_hdf5
!
    endsubroutine output_timeseries
!***********************************************************************
    subroutine hdf5_output_slice_position
!
!  27-Oct-2018/PABourdin: no action required for HDF5 output
!
    endsubroutine hdf5_output_slice_position
!***********************************************************************
    subroutine input_slice_real_arr(file, time, pos, data)
!
!  dummy
!
!  24-may-19/MR: coded
!
      use File_io, only: file_exists

      character (len=*),     intent(in) :: file
      real,                  intent(out):: time
      real,                  intent(out):: pos
      real, dimension(:,:,:),intent(out):: data
!
      call fatal_error('input_slice_real_arr', 'not implemented for HDF5')

      call keep_compiler_quiet(file)
      call keep_compiler_quiet(pos)
      call keep_compiler_quiet(time)
      call keep_compiler_quiet(data)

    endsubroutine input_slice_real_arr
!***********************************************************************
    subroutine input_slice_scat(file,pos,data,ind,nt)
!       
!  dummy
!
!  24-may-19/MR: coded
!
      use General, only: scattered_array

      character (len=*),   intent(in) :: file
      real,                intent(out):: pos
      type(scattered_array),intent(out):: data
      integer :: ind,nt

      call fatal_error('input_slice_scat', 'Not implemented for HDF5')

      call keep_compiler_quiet(file)
      call keep_compiler_quiet(pos)
      call keep_compiler_quiet(ind,nt)

    endsubroutine input_slice_scat
!***********************************************************************
    subroutine input_slice_arr(datadir, time, label, suffix, pos, data)
!       
!  read a slice file
!
!  24-may-19/MR: coded
!
      use File_IO, only: file_exists
      use Mpicomm, only: MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE,mpibarrier
      use General, only: find_proc

      character (len=*),     intent(in) :: datadir,label,suffix
      real, dimension(:),    intent(out):: time
      real,                  intent(out):: pos
      real, dimension(:,:,:),intent(out):: data
!            
      integer :: it, nt, comm, slice_root
      character(LEN=fnlen) :: filename,group
!
      select case (suffix(1:2))
        case ('xy'); comm=MPI_COMM_XYPLANE; slice_root=find_proc(0,0,ipz)
        case ('xz'); comm=MPI_COMM_XZPLANE; slice_root=find_proc(0,ipy,0)
        case ('yz'); comm=MPI_COMM_YZPLANE; slice_root=find_proc(ipx,0,0)
      end select
    
      filename = trim(datadir)//'/slices/'//trim(label)//'_'//trim(suffix)//'.h5'
      if (iproc==slice_root) then
        if (file_exists (filename)) then
          ! find last written slice
          call file_open_hdf5 (filename, global=.false., read_only=.true.,write=.false.)
          if (exists_in_hdf5 ('last')) then
            call input_hdf5 ('last', nt)
          else
            call fatal_error('input_slice_arr','no "last" group in HDF5 file '//trim(filename))
          endif
        else
          call fatal_error('input_slice_arr','no HDF5 file '//trim(filename))
        endif
        call file_close_hdf5
      endif
!
      call mpibarrier(comm)
      call file_open_hdf5 (filename,truncate=.false.,read_only=.true.,write=.false.,comm=comm)

      do it=1,nt
        group=itoa(it)
        call input_hdf5 (trim(group)//'/time', time(it))
        call input_hdf5 (trim(group)//'position', pos)

      ! collect data along 'xy', 'xz', or 'yz'
        select case (suffix(1:2))
        case ('xy')
          call input_hdf5 (trim(group)//'data', (/nxgrid, nygrid/), (/ipx, ipy/), data(:,:,it))
        case ('xz')
          call input_hdf5 (trim(group)//'data', (/nxgrid, nzgrid/), (/ipx, ipz/), data(:,:,it))
        case ('yz')
          call input_hdf5 (trim(group)//'data', (/nygrid, nzgrid/), (/ipy, ipz/), data(:,:,it))
        case default
          call fatal_error ('input_slice', 'unknown 2D slice "'//trim (suffix)//'"', .true.)
        endselect
      enddo

      call file_close_hdf5
!
    endsubroutine input_slice_arr
!***********************************************************************
    subroutine hdf5_output_slice(lwrite, time, label, suffix, pos, grid_pos, data)
!
!  append to a slice file
!
!  30-Oct-2018/PABourdin: coded
!
      use File_io, only: file_exists
      use General, only: itoa, find_proc
      use Mpicomm, only: mpibcast_int, mpibarrier, &
                         MPI_COMM_XYPLANE, MPI_COMM_XZPLANE, MPI_COMM_YZPLANE
!
      logical,                        intent(in) :: lwrite
      real,                           intent(in) :: time, pos
      character (len=*),              intent(in) :: label, suffix
      integer,                        intent(in) :: grid_pos
      real, dimension (:,:), pointer             :: data
!
      character (len=fnlen) :: filename, group
      integer :: last, slice_root, comm, ndim1, ndim2
      real :: time_last
!
      if (.not.(lwrite.and.associated(data))) return

      select case (suffix(1:2))
        case ('xy'); comm=MPI_COMM_XYPLANE; slice_root=find_proc(0,0,ipz)
        case ('xz'); comm=MPI_COMM_XZPLANE; slice_root=find_proc(0,ipy,0)
        case ('yz'); comm=MPI_COMM_YZPLANE; slice_root=find_proc(ipx,0,0)
      end select

      filename = trim(datadir)//'/slices/'//trim(label)//'_'//trim(suffix)//'.h5'
      if (iproc==slice_root) then
        if (file_exists (filename)) then
          ! find last written slice
          call file_open_hdf5 (filename, global=.false., read_only=.true.,write=.true.)
          if (exists_in_hdf5 ('last')) then
            call input_hdf5 ('last', last)
            do last=last,1,-1
              call input_hdf5 (trim(itoa(last))//'/time', time_last)
              if (time > time_last) exit
            enddo
            last=last+1
          else
          endif
        else
          ! create empty file
          call file_open_hdf5 (filename, global=.false.,truncate=.true.,write=.true.)
          last=1
        endif
        call file_close_hdf5
      endif

      call mpibcast_int (last,proc=0,comm=comm)    ! proc=0 as the procressor rank w.r.t. communicator comm is needed
      group = trim(itoa(last))//'/'
  !
      if (iproc==slice_root) then
        call file_open_hdf5 (filename, global=.false., truncate=.false., write=.true.)
        call output_hdf5 ('last', last)
        call create_group_hdf5 (group)
        call output_hdf5 (trim(group)//'time', time)
        call output_hdf5 (trim(group)//'position', pos)
        call output_hdf5 (trim(group)//'coordinate', grid_pos)
        call file_close_hdf5
      endif
!
      call mpibarrier(comm)
      call file_open_hdf5 (filename, truncate=.false., write=.true.,comm=comm)
      ! collect data along 'xy', 'xz', or 'yz'
      ndim1=max(1,size(data,1)); ndim2=max(1,size(data,2))
      select case (suffix(1:2))
      case ('xy')
        call output_hdf5 (trim(group)//'data', data, ndim1, ndim2, nxgrid, nygrid, ipx, ipy)
      case ('xz')
        call output_hdf5 (trim(group)//'data', data, ndim1, ndim2, nxgrid, nzgrid, ipx, ipz)
      case ('yz')
        call output_hdf5 (trim(group)//'data', data, ndim1, ndim2, nygrid, nzgrid, ipy, ipz)
      case default
        call fatal_error ('output_slice', 'unknown 2D slice "'//trim (suffix)//'"', .true.)
      endselect

      call file_close_hdf5
!
    endsubroutine hdf5_output_slice
!***********************************************************************
    subroutine output_profile(name, coord, a, type, lsave_name, lhas_ghost)
!
!  Writes a profile along any direction.
!
!  07-Nov-2018/PABourdin: coded
!
      use General, only: loptest
      use File_io, only: parallel_file_exists
!
      real, dimension(:) :: coord, a
      character (len=*) :: name
      character :: type
      logical, optional :: lsave_name, lhas_ghost
!
      character (len=fnlen) :: filename
      integer :: np, ng, ip, np_global, np1, np2
      logical :: lexists, lwrite, lp1, lp2
!
      if (.not. lwrite_prof) return
!
      np = size(coord)
      ng = 0
      if (loptest (lhas_ghost)) ng = 3
      select case (type)
      case ('x')
        np_global = (np - 2*ng) * nprocx + 2*ng
        ip = ipx
        lp1 = lfirst_proc_x
        lp2 = llast_proc_x
        lwrite = lfirst_proc_yz
      case ('y')
        np_global = (np - 2*ng) * nprocy + 2*ng
        ip = ipy
        lp1 = lfirst_proc_y
        lp2 = llast_proc_y
        lwrite = lfirst_proc_xz
      case ('z')
        np_global = (np - 2*ng) * nprocz + 2*ng
        ip = ipz
        lp1 = lfirst_proc_z
        lp2 = llast_proc_z
        lwrite = lfirst_proc_xy
      case default
        call fatal_error ('output_profile', 'unknown direction "'//type//'"')
      endselect
      np1 = 1
      np2 = np
      if (.not. lp1) np1 = np1 + ng
      if (.not. lp2) np2 = np2 - ng
!
      ! write profile
      filename = trim(directory_snap)//'/'//'profile.h5'
      lexists = parallel_file_exists (filename)
      call file_open_hdf5 (filename, truncate=(.not. lexists))
      call create_group_hdf5 (type)
      call create_group_hdf5 (type//'/'//trim(name))
      call output_hdf5 (type//'/'//trim(name)//'/data', a, np, np_global, ip, np1, np2, ng, lwrite)
      call output_hdf5 (type//'/'//trim(name)//'/position', coord, np, np_global, ip, np1, np2, ng, lwrite)
      call file_close_hdf5
!
    endsubroutine output_profile
!***********************************************************************
    subroutine input_profile(name, type, a, np, lhas_ghost)
!
!  Reads a profile from a file along any direction.
!
!  07-Nov-2018/PABourdin: coded
!
      use General, only: loptest
      use Mpicomm, only: mpibcast_real, MPI_COMM_WORLD
!
      character (len=*), intent(in) :: name
      character, intent(in) :: type
      integer, intent(in) :: np
      real, dimension(np), intent(out) :: a
      logical, optional :: lhas_ghost
!
      character (len=fnlen) :: filename
      integer :: np_global, np1, np2, ng
!
      ng = 0
      if (loptest (lhas_ghost)) ng = 3
      select case (type)
      case ('x')
        np_global = (np - 2*ng) * nprocx + 2*ng
        np1 = 1 + ipx * (np - 2*ng)
      case ('y')
        np_global = (np - 2*ng) * nprocy + 2*ng
        np1 = 1 + ipy * (np - 2*ng)
      case ('z')
        np_global = (np - 2*ng) * nprocz + 2*ng
        np1 = 1 + ipz * (np - 2*ng)
      case default
        call fatal_error ('input_profile', 'unknown direction "'//type//'"')
      endselect
      np2 = np1 + np - 1
!
      ! read profile
      filename = trim(directory_snap)//'/'//'profile.h5'
      call file_open_hdf5 (filename, read_only=.true.)
      call input_hdf5 (type//'/'//trim(name)//'/data', a, np, np_global, np1, np2)
      call file_close_hdf5
!
!  Should we check that coord == z for type == 'z' etc.?
!
    endsubroutine input_profile
!***********************************************************************
    subroutine index_append(varname,ivar,vector,array)
!
! 14-Oct-2018/PABourdin: coded
! 09-Jul-2020/PAB: reworked
!
      character (len=*), intent(in) :: varname
      integer, intent(in) :: ivar
      integer, intent(in) :: vector
      integer, intent(in) :: array
!
      integer, parameter :: lun_output = 92
      character (len=len(varname)) :: quantity
      character (len=2), dimension (9) :: components
      integer :: pos, vec, arr, l
!
      ! omit all unused variables
      if (ivar <= 0) return
!
      if (lroot) open (lun_output, file=trim(datadir)//'/'//trim(index_pro), POSITION='append')
      if ((vector > 0) .and. (abs(array) > 0)) then
        ! expand array first: iuud => [iuud1,iuud2,iuud3,...]
        ! expand vector then: iuud# => [iuud#x,iuud#y,iuud#z] => ivar+(#-1)*vector+[0,1,2]
        do arr = 1, abs(array)
          pos = ivar + (arr-1) * vector
          if ('i'//trim(index_get (pos, quiet=.true.)) /= trim(varname)//trim(itoa(arr))//'x') then
            if (lroot) write (lun_output,*) trim(varname)//trim(itoa(arr))//'x='//trim(itoa(pos))
            call index_register (trim(varname)//trim(itoa(arr))//'x', pos)
          endif
          if ('i'//trim(index_get (pos+1, quiet=.true.)) /= trim(varname)//trim(itoa(arr))//'y') then
            if (lroot) write (lun_output,*) trim(varname)//trim(itoa(arr))//'y='//trim(itoa(pos+1))
            call index_register (trim(varname)//trim(itoa(arr))//'y', pos+1)
          endif
          if ('i'//trim(index_get (pos+2, quiet=.true.)) /= trim(varname)//trim(itoa(arr))//'z') then
            if (lroot) write (lun_output,*) trim(varname)//trim(itoa(arr))//'z='//trim(itoa(pos+2))
            call index_register (trim(varname)//trim(itoa(arr))//'z', pos+2)
          endif
        enddo
      elseif (array > 0) then
        ! backwards compatibility: ind => indgen(N) + ivar
        !!! if (lroot) write (lun_output,*) trim(varname)//'=indgen('//trim(itoa(array))//')+'//trim(itoa(ivar))
        ! expand array: ind => ind[1,...,N] = ivar + [0,...,N-1]
        do pos = 1, array
          if ('i'//trim(index_get (ivar+pos-1, quiet=.true.)) == trim(varname)//trim(itoa(pos))) cycle
          if (lroot) write (lun_output,*) trim(varname)//trim(itoa(pos))//'='//trim(itoa(ivar+pos-1))
          call index_register (trim(varname)//trim(itoa(pos)), ivar+pos-1)
        enddo
      elseif (vector > 0) then
        ! backwards compatibility: iuu => ivar
        !!! if (lroot) write (lun_output,*) trim(varname)//'='//trim(itoa(ivar))
        ! expand vectors
        if (vector == 3) then
          quantity = trim (varname)
          l = len (trim (quantity))
          if (l == 3) then
            ! shortcuts for two-letter names: iuu, iaa, etc.
            if (quantity(2:2) == quantity(3:3)) quantity = trim (quantity(1:2))
          endif
          ! expand vectors: iuu => [iux,iuy,iuz], iaa => [iax,iay,iaz], etc.
          if ('i'//trim(index_get (ivar, quiet=.true.)) /= trim(quantity)//'x') then
            if (lroot) write (lun_output,*) trim(quantity)//'x'//'='//trim(itoa(ivar))
            call index_register (trim(quantity)//'x', ivar)
          endif
          if ('i'//trim(index_get (ivar+1, quiet=.true.)) /= trim(quantity)//'y') then
            if (lroot) write (lun_output,*) trim(quantity)//'y'//'='//trim(itoa(ivar+1))
            call index_register (trim(quantity)//'y', ivar+1)
          endif
          if ('i'//trim(index_get (ivar+2, quiet=.true.)) /= trim(quantity)//'z') then
            if (lroot) write (lun_output,*) trim(quantity)//'z'//'='//trim(itoa(ivar+2))
            call index_register (trim(quantity)//'z', ivar+2)
          endif
        else
          ! expand tensors
          if (vector == 6) then
            ! expand symmetric 3x3 tensor (6 different components)
            components = (/ 'xx', 'xy', 'xz', 'yy', 'yz', 'zz', '  ', '  ', '  ' /)
          else
            ! expand asymmetric 3x3 tensor (9 different components)
            components = (/ 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz' /)
          endif
          do pos = 1, vector
            if ('i'//trim(index_get (ivar+pos-1, quiet=.true.)) == trim(varname)//'_'//trim(components(pos))) cycle
            if (lroot) write (lun_output,*) trim(varname)//'_'//trim(components(pos))//'='//trim(itoa(ivar+pos-1))
            call index_register (trim(varname)//'_'//trim(components(pos)), ivar+pos-1)
          enddo
        endif
      else
        ! scalar: ilnrho => ivar
        if ('i'//trim(index_get (ivar, quiet=.true.)) /= trim(varname)) then
          if (lroot) write (lun_output,*) trim(varname)//'='//trim(itoa(ivar))
          call index_register (trim(varname), ivar)
        endif
      endif
      if (lroot) close(lun_output)
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
      integer, parameter :: lun_output = 92
!
      if ('i'//index_get (ilabel, particle=.true., quiet=.true.) == label) return
      if (lroot) then
        open(lun_output,file=trim(datadir)//'/'//trim(particle_index_pro), POSITION='append')
        write(lun_output,*) trim(label)//'='//trim(itoa(ilabel))
        close(lun_output)
      endif
      call index_register (trim(label), ilabel, particle=.true.)
!
    endsubroutine particle_index_append
!***********************************************************************
    subroutine pointmass_index_append(label,ilabel)
!
! 13-Apr-2019/PABourdin: copied from 'particle_index_append'
!
      character (len=*), intent(in) :: label
      integer, intent(in) :: ilabel
!
      integer, parameter :: lun_output = 92
!
      if ('i'//index_get (ilabel, pointmass=.true., quiet=.true.) == label) return
      if (lroot) then
        open(lun_output,file=trim(datadir)//'/'//trim(pointmass_index_pro), POSITION='append')
        write(lun_output,*) trim(label)//'='//trim(itoa(ilabel))
        close(lun_output)
      endif
      call index_register (trim(label), ilabel, pointmass=.true.)
!
    endsubroutine pointmass_index_append
!***********************************************************************
    function index_get(ivar,particle,pointmass,quiet)
!
! 17-Oct-2018/PABourdin: coded
!
      character (len=labellen) :: index_get
      integer, intent(in) :: ivar
      logical, optional, intent(in) :: particle, pointmass, quiet
!
      type (element), pointer, save :: current => null()
      integer, save :: max_reported = -1
!
      index_get = ''
      current => last
      if (loptest (particle)) current => last_particle
      if (loptest (pointmass)) current => last_pointmass
      do while (associated (current))
        if (current%component == ivar) then
          index_get = current%label(2:len(current%label))
          exit
        endif
        current => current%previous
      enddo
!
      if (lroot .and. .not. loptest (quiet) .and. (index_get == '') .and. (max_reported < ivar)) then
        call warning ('index_get', 'f-array index #'//trim (itoa (ivar))//' not found!')
        if (max_reported == -1) then
          call warning ('index_get', &
              'This likely indicates a mismatch in the mvar/maux contributions of the modules that are active in this setup.')
          call warning ('index_get', &
              'Alternatively, some variables may not have been initialized correctly. Both is an error and should be fixed!')
        endif
        max_reported = ivar
      endif
!
    endfunction index_get
!***********************************************************************
    subroutine index_register(varname,ivar,particle,pointmass)
!
! 17-Oct-2018/PABourdin: coded
!
      character (len=*), intent(in) :: varname
      integer, intent(in) :: ivar
      logical, optional, intent(in) :: particle, pointmass
!
      type (element), pointer, save :: new => null()
!
      if (.not. loptest (particle) .and. .not. loptest (pointmass)) then
        ! ignore variables that are not written
        if ((ivar < 1) .or. (ivar > mfarray)) return
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
      elseif (loptest (pointmass)) then
        if (associated (last_pointmass)) new%previous => last_pointmass
      else
        if (associated (last)) new%previous => last
      endif
      new%label = trim(varname)
      new%component = ivar
      if (loptest (particle)) then
        last_particle => new
      elseif (loptest (pointmass)) then
        last_pointmass => new
      else
        last => new
      endif
!
    endsubroutine index_register
!***********************************************************************
    subroutine index_reset
!
! 14-Oct-2018/PABourdin: coded
!
      type (element), pointer, save :: current => null()
      integer, parameter :: lun_output = 92
!
      if (lroot) then
        open(lun_output,file=trim(datadir)//'/'//trim(index_pro),status='replace')
        close(lun_output)
        open(lun_output,file=trim(datadir)//'/'//trim(particle_index_pro),status='replace')
        close(lun_output)
        open(lun_output,file=trim(datadir)//'/'//trim(pointmass_index_pro),status='replace')
        close(lun_output)
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
      do while (associated (last_pointmass))
        current => last_pointmass
        last_pointmass => last%previous
        deallocate (current)
        nullify (current)
      enddo
!
    endsubroutine index_reset
!***********************************************************************
    subroutine output_settings(time, time_only)
!
!  Write additional settings and grid.
!
!  13-Nov-2018/PABourdin: moved from other functions
!
      use General, only: loptest
      use Mpicomm, only: collect_grid
      use Syscalls, only: sizeof_real
!
      real, optional, intent(in) :: time
      logical, optional, intent(in) :: time_only
!
      real, dimension(:), allocatable :: gx, gy, gz
      integer :: alloc_err, mxgrid_, mygrid_, mzgrid_
!
      if (lroot.and.present(time)) call output_hdf5 ('time', time)
      if (loptest(time_only)) return
!
      mxgrid_=global_size(1); mygrid_=global_size(2); mzgrid_=global_size(3);
!
      if (lroot) then
        allocate (gx(mxgrid_), gy(mygrid_), gz(mzgrid_), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('output_settings', 'allocate memory for gx,gy,gz', .true.)
      endif
!
      call collect_grid (x(:local_size(1)), y(:local_size(2)), z(:local_size(3)), gx, gy, gz)
      if (lroot) then
        call create_group_hdf5 ('grid')
        call output_hdf5 ('grid/x', gx, mxgrid_)
        call output_hdf5 ('grid/y', gy, mygrid_)
        call output_hdf5 ('grid/z', gz, mzgrid_)
        call output_hdf5 ('grid/dx', dx)
        call output_hdf5 ('grid/dy', dy)
        call output_hdf5 ('grid/dz', dz)
        call output_hdf5 ('grid/Lx', Lx)
        call output_hdf5 ('grid/Ly', Ly)
        call output_hdf5 ('grid/Lz', Lz)
        call output_hdf5 ('grid/Ox', x0)
        call output_hdf5 ('grid/Oy', y0)
        call output_hdf5 ('grid/Oz', z0)
      endif
      call collect_grid (dx_1(:local_size(1)), dy_1(:local_size(2)), dz_1(:local_size(3)), gx, gy, gz)
      if (lroot) then
        call output_hdf5 ('grid/dx_1', gx, mxgrid_)
        call output_hdf5 ('grid/dy_1', gy, mygrid_)
        call output_hdf5 ('grid/dz_1', gz, mzgrid_)
      endif
      call collect_grid (dx_tilde(:local_size(1)), dy_tilde(:local_size(2)), dz_tilde(:local_size(3)), gx, gy, gz)
      if (lroot) then
        call output_hdf5 ('grid/dx_tilde', gx, mxgrid_)
        call output_hdf5 ('grid/dy_tilde', gy, mygrid_)
        call output_hdf5 ('grid/dz_tilde', gz, mzgrid_)
        call create_group_hdf5 ('unit')
        call output_hdf5 ('unit/system', unit_system)
        call output_hdf5_double ('unit/density', unit_density)
        call output_hdf5_double ('unit/length', unit_length)
        call output_hdf5_double ('unit/velocity', unit_velocity)
        call output_hdf5_double ('unit/magnetic', unit_magnetic)
        call output_hdf5_double ('unit/temperature', unit_temperature)
        call output_hdf5_double ('unit/mass', unit_mass)
        call output_hdf5_double ('unit/energy', unit_energy)
        call output_hdf5_double ('unit/time', unit_time)
        call output_hdf5_double ('unit/flux', unit_flux)
        call create_group_hdf5 ('settings')
        call output_hdf5 ('settings/mx', mxgrid_)
        call output_hdf5 ('settings/my', mygrid_)
        call output_hdf5 ('settings/mz', mzgrid_)
        call output_hdf5 ('settings/nx', mxgrid_-2*nghost)
        call output_hdf5 ('settings/ny', mygrid_-2*nghost)
        call output_hdf5 ('settings/nz', mzgrid_-2*nghost)
        call output_hdf5 ('settings/l1', nghost)
        call output_hdf5 ('settings/m1', nghost)
        call output_hdf5 ('settings/n1', nghost)
        call output_hdf5 ('settings/l2', mxgrid_-nghost-1)
        call output_hdf5 ('settings/m2', mygrid_-nghost-1)
        call output_hdf5 ('settings/n2', mzgrid_-nghost-1)
        call output_hdf5 ('settings/nghost', nghost)
        call output_hdf5 ('settings/mvar', mvar_out)
        call output_hdf5 ('settings/maux', maux_out)
        call output_hdf5 ('settings/mglobal', mglobal)
        call output_hdf5 ('settings/nprocx', nprocx)
        call output_hdf5 ('settings/nprocy', nprocy)
        call output_hdf5 ('settings/nprocz', nprocz)
        if (sizeof_real() < 8) then
          call output_hdf5 ('settings/precision', 'S')
        else
          call output_hdf5 ('settings/precision', 'D')
        endif
        ! versions represent only non-compatible file formats
        ! 0 : experimental
        ! 1 : first public release
        call output_hdf5 ('settings/version', 0)
      endif
!
    endsubroutine output_settings
!***********************************************************************
endmodule HDF5_IO
