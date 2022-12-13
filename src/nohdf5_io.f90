! $Id$
!
!  This module is a dummy module for non-HDF5 IO.
!  28-Oct-2016/PABoudin: coded
!
module HDF5_IO
!
  use Cdata
  use General, only: keep_compiler_quiet, itoa, numeric_precision
  use Messages, only: fatal_error
  use Mpicomm, only: lroot
!
  implicit none
!
  interface hdf5_input_slice
    module procedure input_slice_real_arr
    module procedure input_slice_scat_arr
  endinterface
!
  interface input_hdf5
    module procedure input_hdf5_int_0D
    module procedure input_hdf5_int_1D
    module procedure input_hdf5_0D
    module procedure input_hdf5_1D
    module procedure input_hdf5_part_2D
    module procedure input_hdf5_profile_1D
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
    module procedure output_hdf5_part_2D
    module procedure output_hdf5_profile_1D
    module procedure output_local_hdf5_2D
    module procedure output_hdf5_slice_2D
    module procedure output_local_hdf5_3D
    module procedure output_hdf5_3D
    module procedure output_hdf5_4D
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
  integer, parameter :: lun_input = 89, lun_output = 92
!
  contains
!***********************************************************************
    subroutine initialize_hdf5(nxyz,ngrid)
!
      integer, dimension(3), optional :: nxyz,ngrid
      ! nothing to do
      call keep_compiler_quiet(nxyz)
      call keep_compiler_quiet(ngrid)
!
    endsubroutine initialize_hdf5
!***********************************************************************
    subroutine init_hdf5
!
      ! nothing to do
!
    endsubroutine init_hdf5
!***********************************************************************
    subroutine finalize_hdf5
!
      call fatal_error ('finalize_hdf5', 'You can not use HDF5 without setting an HDF5_IO module.')
!
    endsubroutine finalize_hdf5
!***********************************************************************
    subroutine file_open_hdf5(file, truncate, global, read_only, write, comm)
!
      character (len=*), intent(in) :: file
      logical, optional, intent(in) :: truncate
      logical, optional, intent(in) :: global
      logical, optional, intent(in) :: read_only
      logical, optional, intent(in) :: write
      integer, optional, intent(in) :: comm
!
      call keep_compiler_quiet(file)
      call keep_compiler_quiet(.true.,truncate, global, read_only)
      call keep_compiler_quiet(.true.,write)
      call keep_compiler_quiet(comm)

      call fatal_error ('file_open_hdf5', 'You can not use HDF5 without setting an HDF5_IO module.')
!
    endsubroutine file_open_hdf5
!***********************************************************************
    subroutine file_close_hdf5
!
      call fatal_error ('file_close_hdf5', 'You can not use HDF5 without setting an HDF5_IO module.')
!
    endsubroutine file_close_hdf5
!***********************************************************************
    subroutine create_group_hdf5(name)
!
      character (len=*), intent(in) :: name
!
      call fatal_error ('create_group_hdf5', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
!
    endsubroutine create_group_hdf5
!***********************************************************************
    logical function exists_in_hdf5(name)
!
      character (len=*), intent(in) :: name
!
      call fatal_error ('exists_in_hdf5', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      exists_in_hdf5 = .false.
!
    endfunction exists_in_hdf5
!***********************************************************************
    subroutine input_hdf5_int_0D(name, data)
!
      character (len=*), intent(in) :: name
      integer, intent(out) :: data
!
      call fatal_error ('input_hdf5_int_0D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      data = -1
!
    endsubroutine input_hdf5_int_0D
!***********************************************************************
    subroutine input_hdf5_int_1D(name, data, nv, same_size)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      integer, dimension (nv), intent(out) :: data
      logical, optional, intent(in) :: same_size
!
      call fatal_error ('input_hdf5_int_1D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      data(:) = -1
      call keep_compiler_quiet(.true., same_size)
!
    endsubroutine input_hdf5_int_1D
!***********************************************************************
    subroutine input_hdf5_0D(name, data)
!
      character (len=*), intent(in) :: name
      real, intent(out) :: data
!
      call fatal_error ('input_hdf5_0D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      data = -1.0
!
    endsubroutine input_hdf5_0D
!***********************************************************************
    subroutine input_hdf5_1D(name, data, nv, same_size)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (nv), intent(out) :: data
      logical, optional, intent(in) :: same_size
!
      call fatal_error ('input_hdf5_1D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      data(:) = -1.0
      call keep_compiler_quiet(.true., same_size)
!
    endsubroutine input_hdf5_1D
!***********************************************************************
    subroutine input_hdf5_part_2D(name, data, mv, nc, nv)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: mv, nc
      real, dimension (mv,mparray), intent(out) :: data
      integer, intent(out) :: nv
!
      call fatal_error ('input_hdf5_part_2D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(mv)
      call keep_compiler_quiet(nc)
      data(:,:) = -1.0
      nv = -1
!
    endsubroutine input_hdf5_part_2D
!***********************************************************************
    subroutine input_hdf5_profile_1D(name, data, ldim, gdim, np1, np2)
!
      character (len=*), intent(in) :: name
      real, dimension (:) :: data
      integer, intent(in) :: ldim, gdim, np1, np2
!
      call fatal_error ('input_hdf5_profile_1D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(ldim)
      call keep_compiler_quiet(gdim)
      call keep_compiler_quiet(np1)
      call keep_compiler_quiet(np2)
!
    endsubroutine input_hdf5_profile_1D
!***********************************************************************
    subroutine input_hdf5_3D(name, data)
!
      character (len=*), intent(in) :: name
      real, dimension (mx,my,mz), intent(out) :: data
!
      call fatal_error ('input_hdf5_3D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      data(:,:,:) = -1.0
!
    endsubroutine input_hdf5_3D
!***********************************************************************
    subroutine input_hdf5_4D(name, data, nv)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(out) :: data
!
      call fatal_error ('input_hdf5_4D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      data(:,:,:,:) = -1.0
!
    endsubroutine input_hdf5_4D
!***********************************************************************
    subroutine output_hdf5_string(name, data)
!
      character (len=*), intent(in) :: name
      character (len=*), intent(in) :: data
!
      call fatal_error ('output_hdf5_string', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name,data)
!
    endsubroutine output_hdf5_string
!***********************************************************************
    subroutine output_hdf5_torus_rect(name, data)
!
      use Geometrical_types, only: torus_rect

      character (len=*), intent(in) :: name
      type(torus_rect), intent(in) :: data

      call fatal_error ('output_hdf5_torus_rect', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data%height)
!
    endsubroutine output_hdf5_torus_rect
!***********************************************************************
    subroutine output_hdf5_int_0D(name, data)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: data
!
      call fatal_error ('output_hdf5_int_0D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
!
    endsubroutine output_hdf5_int_0D
!***********************************************************************
    subroutine output_hdf5_int_1D(name, data, nv, same_size)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      integer, dimension(nv), intent(in) :: data
      logical, optional, intent(in) :: same_size
!
      call fatal_error ('output_hdf5_int_1D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
      call keep_compiler_quiet(.true., same_size)
!
    endsubroutine output_hdf5_int_1D
!***********************************************************************
    subroutine output_hdf5_0D(name, data)
!
      character (len=*), intent(in) :: name
      real, intent(in) :: data
!
      call fatal_error ('output_hdf5_0D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
!
    endsubroutine output_hdf5_0D
!***********************************************************************
    subroutine output_hdf5_1D(name, data, nv, same_size)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (nv), intent(in) :: data
      logical, optional, intent(in) :: same_size
!
      call fatal_error ('output_hdf5_1D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
      call keep_compiler_quiet(.true., same_size)
!
    endsubroutine output_hdf5_1D
!***********************************************************************
    subroutine output_hdf5_part_2D(name, data, mv, nc, nv)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: mv, nc
      real, dimension (mv,mparray), intent(in) :: data
      integer, intent(in) :: nv
!
      call fatal_error ('output_hdf5_part_2D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
      call keep_compiler_quiet(mv)
      call keep_compiler_quiet(nc)
      call keep_compiler_quiet(nv)
!
    endsubroutine output_hdf5_part_2D
!***********************************************************************
    subroutine output_hdf5_profile_1D(name, data, ldim, gdim, ip, np1, np2, ng, lhas_data)
!
      character (len=*), intent(in) :: name
      real, dimension (:), intent(in) :: data
      integer, intent(in) :: ldim, gdim, ip, np1, np2, ng
      logical, intent(in) :: lhas_data
!
      call fatal_error ('output_hdf5_profile_1D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
      call keep_compiler_quiet(ldim)
      call keep_compiler_quiet(gdim)
      call keep_compiler_quiet(ip)
      call keep_compiler_quiet(np1)
      call keep_compiler_quiet(np2)
      call keep_compiler_quiet(ng)
      call keep_compiler_quiet(lhas_data)
!
    endsubroutine output_hdf5_profile_1D
!***********************************************************************
    subroutine output_local_hdf5_2D(name, data, dim1, dim2)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: dim1, dim2
      real, dimension (dim1,dim2), intent(in) :: data
!
      call fatal_error ('output_local_hdf5_2D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
      call keep_compiler_quiet(dim1)
      call keep_compiler_quiet(dim2)
!
    endsubroutine output_local_hdf5_2D
!***********************************************************************
    subroutine output_hdf5_slice_2D(name, data, ldim1, ldim2, gdim1, gdim2, ip1, ip2, lhas_data)
!
      character (len=*), intent(in) :: name
      real, dimension (:,:), pointer :: data
      integer, intent(in) :: ldim1, ldim2, gdim1, gdim2, ip1, ip2
      logical, intent(in) :: lhas_data
!
      call fatal_error ('output_hdf5_slice_2D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
      call keep_compiler_quiet(ldim1)
      call keep_compiler_quiet(ldim2)
      call keep_compiler_quiet(gdim1)
      call keep_compiler_quiet(gdim2)
      call keep_compiler_quiet(ip1)
      call keep_compiler_quiet(ip2)
      call keep_compiler_quiet(lhas_data)
!
    endsubroutine output_hdf5_slice_2D
!***********************************************************************
    subroutine output_local_hdf5_3D(name, data, dim1, dim2, dim3)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: dim1, dim2, dim3
      real, dimension (dim1,dim2,dim3), intent(in) :: data
!
      call fatal_error ('output_local_hdf5_3D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
      call keep_compiler_quiet(dim1)
      call keep_compiler_quiet(dim2)
      call keep_compiler_quiet(dim3)
!
    endsubroutine output_local_hdf5_3D
!***********************************************************************
    subroutine output_hdf5_3D(name, data)
!
      character (len=*), intent(in) :: name
      real, dimension (mx,my,mz), intent(in) :: data
!
      call fatal_error ('output_hdf5_3D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
!
    endsubroutine output_hdf5_3D
!***********************************************************************
    subroutine output_hdf5_4D(name, data, nv)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: data
!
      call fatal_error ('output_hdf5_4D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
!
    endsubroutine output_hdf5_4D
!***********************************************************************
    subroutine input_dim(wrkdir, mx_in, my_in, mz_in, mvar_in, maux_in, mglobal_in, &
                         prec_in, nghost_in, nprocx_in, nprocy_in, nprocz_in, local)
!
!  Read dimensions from dim.dat (local or global).
!
      use General, only: loptest

      character (len=*), intent(in) :: wrkdir
      integer, intent(out) :: mx_in, my_in, mz_in, mvar_in, maux_in, mglobal_in
      integer, intent(out) :: nprocx_in, nprocy_in, nprocz_in, nghost_in
      logical, optional :: local
!
      character (len=fnlen) :: filename
      character :: prec_in
      integer :: iprocz_slowest, ipx_in, ipy_in, ipz_in
!
      if (lroot) then
        ! local or global dimension file
        if (loptest(local)) then
          filename = trim(wrkdir)//'/data/proc0/dim.dat'
        else
          filename = trim(wrkdir)//'/data/dim.dat'
        endif
        open(lun_input,file=filename)
        read(lun_input,'(3i7,3i7)') mx_in, my_in, mz_in, mvar_in, maux_in, mglobal_in
        read(lun_input,'(a)') prec_in
        read(lun_input,'(3i5)') nghost_in
        if (loptest(local)) then
          read(lun_input,'(4i5)') ipx_in, ipy_in, ipz_in
          nprocx_in=ipx_in; nprocy_in=ipy_in; nprocz_in=ipz_in
        else
          read(lun_input,'(4i5)') nprocx_in, nprocy_in, nprocz_in, iprocz_slowest
        endif
        close(lun_input)
      endif

    endsubroutine input_dim
!***********************************************************************
    subroutine output_hdf5_double_0D(name, data)
!
      character (len=*), intent(in) :: name
      double precision, intent(in) :: data
!
      call fatal_error ('output_hdf5_double_0D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
!
    endsubroutine output_hdf5_double_0D
!***********************************************************************
    subroutine output_hdf5_double_1D(name, data, nv)
!
      character (len=*), intent(in) :: name
      integer, intent(in) :: nv
      double precision, dimension (nv), intent(in) :: data
!
      call fatal_error ('output_hdf5_double_1D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(nv)
!
    endsubroutine output_hdf5_double_1D
!***********************************************************************
    subroutine output_dim(file, mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar_out, maux_out, mglobal)
!
!  Write dimension to file.
!
!  02-Nov-2018/PABourdin: moved IO parts from wdim to low-level IO modules
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: mx_out, my_out, mz_out, mxgrid_out, mygrid_out, mzgrid_out, mvar_out, maux_out, mglobal
!
      character (len=fnlen) :: filename
      integer :: iprocz_slowest
!
      if (lroot) then
        ! global dimension file
        filename = trim(datadir)//'/'//trim(file)
        open(lun_output,file=filename)
        write(lun_output,'(3i7,3i7)') mxgrid_out, mygrid_out, mzgrid_out, mvar_out, maux_out, mglobal
        write(lun_output,'(a)') numeric_precision()
        write(lun_output,'(3i5)') nghost, nghost, nghost ! why do we write three always identical numbers?
        iprocz_slowest = 0
        if (lprocz_slowest) iprocz_slowest=1
        write(lun_output,'(4i5)') nprocx, nprocy, nprocz, iprocz_slowest
        close(lun_output)
      endif
!
      if (.not. lmonolithic_io) then
        ! write processor dimension files to "proc#/" for distributed IO cases
        filename = trim(directory)//'/'//trim(file)
        open(lun_output,file=filename)
        write(lun_output,'(3i7,3i7)') mx_out, my_out, mz_out, mvar_out, maux_out, mglobal
        write(lun_output,'(a)') numeric_precision()
        write(lun_output,'(3i5)') nghost, nghost, nghost ! why do we write three always identical numbers?
        write(lun_output,'(3i5)') ipx, ipy, ipz
        close(lun_output)
      endif
!
    endsubroutine output_dim
!***********************************************************************
    subroutine wdim_default_grid(file)
!
!  Write dimension to file.
!
!  02-Nov-2018/PABourdin: redesigned
!
      character (len=*), intent(in) :: file
!
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
    subroutine output_timeseries(data, data_im)
!   
!  Append diagnostic data to a binary file.
!   
!  01-Apr-2019/PABourdin: coded
! 
      real, dimension(2*nname), intent(in) :: data, data_im
! 
      ! dummy routine   
!                       
    endsubroutine output_timeseries
!***********************************************************************
    subroutine input_slice_position(directory,ix_bc_,iy_bc_,iy2_bc_,iz_bc_,iz2_bc_,iz3_bc_,iz4_bc_)
!
!  'data/procN/slice_position.dat' is distributed, but may not be synchronized
!  on I/O error (-> dist=0) as this would make it disfunctional; correct a posteriori if necessary.
!
!  24-May-2019/MR: cloned from output_slice_position
!
      character(LEN=*) :: directory
      integer, intent(out) :: ix_bc_,iy_bc_,iy2_bc_,iz_bc_,iz2_bc_,iz3_bc_,iz4_bc_
      logical :: lexist_slice_xy, lexist_slice_xy2, lexist_slice_xy3, lexist_slice_xy4, &
                 lexist_slice_xz, lexist_slice_xz2, lexist_slice_yz

      open (lun_input, file=trim(directory)//'/data/slice_position.dat', STATUS='unknown')
      read (lun_input, '(l5,i5)') lexist_slice_xy, iz_bc_
      read (lun_input, '(l5,i5)') lexist_slice_xy2, iz2_bc_
      read (lun_input, '(l5,i5)') lexist_slice_xy3, iz3_bc_
      read (lun_input, '(l5,i5)') lexist_slice_xy4, iz4_bc_
      read (lun_input, '(l5,i5)') lexist_slice_xz, iy_bc_
      read (lun_input, '(l5,i5)') lexist_slice_xz2, iy2_bc_
      read (lun_input, '(l5,i5)') lexist_slice_yz, ix_bc_
      close(lun_input)
!
    endsubroutine input_slice_position
!***********************************************************************
    subroutine input_slice_real_arr(file, time, pos, data)
!
!  read a slice file
!
!  24-may-19/MR: coded
!
      use File_io, only: file_exists

      character (len=*),   intent(in) :: file
      real,                intent(out):: time
      real,                intent(out):: pos
      real, dimension(:,:,:),intent(out):: data
!
      integer :: nt, ios
!
      if (.not.file_exists(file)) &
        call fatal_error('input_slice', 'no slices file '//trim(file))

      open(lun_input, file=file, form='unformatted')
      nt=0; ios=0
      do while(ios==0)
        read(lun_input,iostat=ios) data(:,:,nt+1), time, pos
        if (ios/=0) exit
        nt=nt+1
      enddo
      close(lun_input)
!
    endsubroutine input_slice_real_arr
!***********************************************************************
    subroutine input_slice_scat_arr(file, pos, data, ivar, nt)
!
!  read a slice file
!
!  24-may-19/MR: coded
!
      use General, only: scattered_array, store_scattered_array
      use File_io, only: file_exists

      character (len=*),   intent(in) :: file
      real,                intent(out):: pos
      type(scattered_array), pointer  :: data   !intent(inout)
      integer,             intent(in) :: ivar
      integer,             intent(out):: nt
!
      integer :: ios
      real :: time
      real,dimension(data%dim1,data%dim2) :: slc
!
      if (.not.file_exists(file)) &
        call fatal_error('input_slice', 'no slices file '//trim(file))

      open(lun_input, file=file, form='unformatted')
      nt=0; ios=0
      do while(ios==0)
        read(lun_input,iostat=ios) slc, time, pos
        if (ios/=0) exit
        nt=nt+1
!if (ivar==1.and.iproc==120) print*, 'nt=', nt
        call store_scattered_array(ivar,nt,slc,data,time)
      enddo
      close(lun_input)
!
    endsubroutine input_slice_scat_arr
!***********************************************************************
    subroutine hdf5_output_slice_position
!
!  'data/procN/slice_position.dat' is distributed, but may not be synchronized
!  on I/O error (-> dist=0) as this would make it disfunctional; correct a posteriori if necessary.
!
!  27-Oct-2018/PABourdin: cleaned up
!
      use Slices_methods, only: write_rslice_position

      open (lun_output, file=trim(directory)//'/slice_position.dat', STATUS='unknown')
      write (lun_output, '(l5,i5," XY")' ) lwrite_slice_xy, iz_loc
      write (lun_output, '(l5,i5," XY2")') lwrite_slice_xy2, iz2_loc
      write (lun_output, '(l5,i5," XY3")') lwrite_slice_xy3, iz3_loc
      write (lun_output, '(l5,i5," XY4")') lwrite_slice_xy4, iz4_loc
      write (lun_output, '(l5,i5," XZ")' ) lwrite_slice_xz, iy_loc
      write (lun_output, '(l5,i5," XZ2")') lwrite_slice_xz2, iy2_loc
      write (lun_output, '(l5,i5," YZ")' ) lwrite_slice_yz, ix_loc
      call write_rslice_position(lun_output)

      close (lun_output)
!
    endsubroutine hdf5_output_slice_position
!***********************************************************************
    subroutine hdf5_output_slice(lwrite, time, label, suffix, pos, grid_pos, data)
!
!  append to a slice file
!
!  12-nov-02/axel: coded
!  26-jun-06/anders: moved to slices
!  22-sep-07/axel: changed Xy to xy2, to be compatible with Mac
!  28-Oct-2018/PABourdin: cleaned up, moved to IO module
!
      logical, intent(in) :: lwrite
      real, intent(in) :: time
      character (len=*), intent(in) :: label, suffix
      real, intent(in) :: pos
      integer, intent(in) :: grid_pos
      real, dimension (:,:), pointer :: data
!
      if (.not. lwrite .or. .not. associated(data)) return
!
!  files data/procN/slice*.* are distributed and will be synchronized a-posteriori on I/O error
!
      open (lun_output, file=trim(directory)//'/slice_'//trim(label)//'.'//trim(suffix), &
            form='unformatted', position='append')
      write (lun_output) data, time, pos
      close (lun_output)
!
    endsubroutine hdf5_output_slice
!***********************************************************************
    subroutine index_append(varname,ivar,vector,array)
!
! 14-oct-18/PAB: coded
! 09-Jul-2020/PAB: reworked
!
      character (len=*), intent(in) :: varname
      integer, intent(in) :: ivar
      integer, intent(in) :: vector
      integer, intent(in) :: array
!
      character (len=len(varname)) :: component
      integer :: pos, l
!
      if (lroot) then
        open (lun_output, file=trim (datadir)//'/'//trim (index_pro), POSITION='append')
        if ((vector > 0) .and. (array > 0)) then
          ! expand array first: iuud => [iuud1,iuud2,iuud3,...]
          ! expand vector then: iuud# => [iuud#x,iuud#y,iuud#z] => ivar+(#-1)*vector+[0,1,2]
          do pos = 1, array
            write (lun_output,*) trim(varname)//trim(itoa(pos))//'x='//trim(itoa(ivar+(pos-1)*vector))
            write (lun_output,*) trim(varname)//trim(itoa(pos))//'y='//trim(itoa(ivar+(pos-1)*vector+1))
            write (lun_output,*) trim(varname)//trim(itoa(pos))//'z='//trim(itoa(ivar+(pos-1)*vector+2))
          enddo
        elseif ((vector > 0) .and. (array < 0)) then
          write (lun_output,*) trim(varname)//'=indgen('//trim(itoa(-array))//')*'//trim(itoa(vector))//'+'//trim(itoa(ivar))
        elseif (array /= 0) then    ! i.e. vector<0
          ! backwards compatibility: ind => indgen(array)+ivar
          write (lun_output,*) trim (varname)//'=indgen('//trim(itoa(abs(array)))//')+'//trim (itoa (ivar))
          ! expand array: ind => [ind1,ind2,...]
          do pos = 1, array
            write (lun_output,*) trim(varname)//trim(itoa(pos))//'='//trim(itoa(ivar+(pos-1)))
          enddo
        elseif (vector > 0) then    ! i.e. array=0
          ! expand vectors: iuu => [iux,iuy,iuz], iaa => [iax,iay,iaz], etc.
          if (vector == 3) then
            ! backwards compatibility: write original vector
            write (lun_output,*) trim(varname)//'='//trim (itoa(ivar))
            ! apply shortcuts
            component = trim (varname)
            l = len (trim (component))
            if (l == 3) then
              ! double endings: iuu, iaa, etc.
              if (component(2:2) == component(3:3)) l = 2
            endif
            write (lun_output,*) trim(component(1:l))//'x='//trim (itoa(ivar))
            write (lun_output,*) trim(component(1:l))//'y='//trim (itoa(ivar+1))
            write (lun_output,*) trim(component(1:l))//'z='//trim (itoa(ivar+2))
          elseif (vector == 6) then
            ! expand symmetric 3x3 tensor (6 different components)
            write (lun_output,*) trim(varname)//'_xx='//trim (itoa(ivar))
            write (lun_output,*) trim(varname)//'_xy='//trim (itoa(ivar+1))
            write (lun_output,*) trim(varname)//'_xz='//trim (itoa(ivar+2))
            write (lun_output,*) trim(varname)//'_yy='//trim (itoa(ivar+3))
            write (lun_output,*) trim(varname)//'_yz='//trim (itoa(ivar+4))
            write (lun_output,*) trim(varname)//'_zz='//trim (itoa(ivar+5))
          elseif (vector==9) then
            ! expand asymmetric 3x3 tensor (9 different components)
            write (lun_output,*) trim(varname)//'_xx='//trim (itoa(ivar))
            write (lun_output,*) trim(varname)//'_xy='//trim (itoa(ivar+1))
            write (lun_output,*) trim(varname)//'_xz='//trim (itoa(ivar+2))
            write (lun_output,*) trim(varname)//'_yx='//trim (itoa(ivar+3))
            write (lun_output,*) trim(varname)//'_yy='//trim (itoa(ivar+4))
            write (lun_output,*) trim(varname)//'_yz='//trim (itoa(ivar+5))
            write (lun_output,*) trim(varname)//'_zx='//trim (itoa(ivar+6))
            write (lun_output,*) trim(varname)//'_zy='//trim (itoa(ivar+7))
            write (lun_output,*) trim(varname)//'_zz='//trim (itoa(ivar+8))
          else
            call fatal_error('index_append','unknown tensor type')
          endif
        else
          ! scalar: ilnrho => ivar
          write (lun_output,*) trim(varname)//'='//trim(itoa(ivar))
        endif
        close (lun_output)
      endif
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
        open(lun_output,file=trim(datadir)//'/'//trim(particle_index_pro), POSITION='append')
        write(lun_output,*) trim(label)//'='//trim(itoa(ilabel))
        close(lun_output)
      endif
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
      if (lroot) then
        open(lun_output,file=trim(datadir)//'/'//trim(pointmass_index_pro), POSITION='append')
        write(lun_output,*) trim(label)//'='//trim(itoa(ilabel))
        close(lun_output)
      endif
!
    endsubroutine pointmass_index_append
!***********************************************************************
    function index_get(ivar,particle,pointmass,quiet)
!
! 13-Nov-2018/PABourdin: coded
!
      character (len=labellen) :: index_get
      integer, intent(in) :: ivar
      logical, optional, intent(in) :: particle, pointmass, quiet
!
      call fatal_error ('index_get', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(ivar)
      call keep_compiler_quiet(particle)
      call keep_compiler_quiet(pointmass)
      call keep_compiler_quiet(quiet)
!
      index_get=' '
!
    endfunction index_get
!***********************************************************************
    subroutine index_reset
!
! 14-Oct-2018/PABourdin: coded
!
      if (lroot) then
! Axel, please uncomment for debugging on Beskow:
! write (123,*) '============> replace "index.pro" with empty file'
        open(lun_output,file=trim(datadir)//'/'//trim(index_pro),status='replace')
        close(lun_output)
        open(lun_output,file=trim(datadir)//'/'//trim(particle_index_pro),status='replace')
        close(lun_output)
        open(lun_output,file=trim(datadir)//'/'//trim(pointmass_index_pro),status='replace')
        close(lun_output)
      endif
!
    endsubroutine index_reset
!***********************************************************************
    subroutine output_profile(fname, coord, a, type, lsave_name, lhas_ghost)
!
!  Writes a profile to a file.
!
!  10-jul-05/axel: coded
!  07-Nov-2018/PABourdin: moved to IO modules
!
      use General, only: loptest
!
      real, dimension(:) :: coord, a
      character (len=*) :: fname
      character :: type
      logical, optional :: lsave_name, lhas_ghost
!
      integer :: pos, np
!
!  If within a loop, do this only for the first step (indicated by lwrite_prof).
!
      if (lwrite_prof) then
        ! Write profile.
        open(lun_output,file=trim(directory)//'/'//type//'prof_'//trim(fname)//'.dat',position='append')
        np = size(coord)
        do pos=1, np
          write(lun_output,*) coord(pos), a(pos)
        enddo
        close(lun_output)
!
        ! Add file name to list of profiles if lsave_name is true.
        if (loptest(lsave_name)) then
          open(lun_output,file=trim(directory)//'/'//type//'prof_list.dat',position='append')
          write(lun_output,*) fname
          close(lun_output)
        endif
      endif
!
    endsubroutine output_profile
!***********************************************************************
    subroutine input_profile(fname, type, a, np, lhas_ghost)
!
!  Read a profile from a file (taken from stratification, for example).
!
!  15-may-15/axel: coded
!  05-Nov-2018/PABourdin: generalized to any direction and moved to IO module
!
      character (len=*), intent(in) :: fname
      character, intent(in) :: type
      integer, intent(in) :: np
      real, dimension(np), intent(out) :: a
      logical, optional :: lhas_ghost
!
      integer :: pos
      real, dimension(np) :: coord
!
      ! Read profile.
      open(lun_input,file=trim(directory)//type//'/prof_'//trim(fname)//'.dat')
      do pos=1, np
        read(lun_input,*) coord(pos), a(pos)
      enddo
      close(lun_input)
!
!  Should we check that coord == z for type == 'z'?
!
    endsubroutine input_profile
!***********************************************************************
    subroutine output_average_1D(path, label, nc, name, data, time, lbinary, lwrite, header)
!
!   Output 1D average to a file.
!
!   16-Nov-2018/PABourdin: coded
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
!
      if (.not. lwrite .or. (nc <= 0)) return
!
      filename = trim(path) // '/' // trim(label) // 'averages.dat'
      if (trim (label) == 'phi_z') filename = trim(path) // '/phizaverages.dat'
      if (lbinary) then
        open(lun_output, file=filename, form='unformatted', position='append')
        if (present (header)) write(lun_output) header
        write(lun_output) time
        write(lun_output) data(:,1:nc)
        close(lun_output)
      else
        open(lun_output, file=filename, position='append')
        if (present (header)) write(lun_output,'(1p,8'//trim(fmt_avgs)//')') header
        write(lun_output,'(1pe12.5)') time
        write(lun_output,'(1p,8'//trim(fmt_avgs)//')') data(:,1:nc)
        close(lun_output)
      endif
!
    endsubroutine output_average_1D
!***********************************************************************
    subroutine output_average_1D_chunked(path, label, nc, name, data, full, time, lbinary, lwrite, header)
!
!   Output 1D chunked average to a file.
!
!   16-Nov-2018/PABourdin: coded
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
      if (present (header)) then
        call output_average_2D(path, label, nc, name, data, time, lbinary, lwrite, header)
      else
        call output_average_2D(path, label, nc, name, data, time, lbinary, lwrite)
      endif
!
    endsubroutine output_average_1D_chunked
!***********************************************************************
    subroutine output_average_2D(path, label, nc, name, data, time, lbinary, lwrite, header)
!
!   Output average to a file.
!
!   16-Nov-2018/PABourdin: coded
!
      character (len=*), intent(in) :: path, label
      integer, intent(in) :: nc
      character (len=fmtlen), dimension(nc), intent(in) :: name
      real, dimension(:,:,:), intent(in) :: data
      real, intent(in) :: time
      logical, intent(in) :: lbinary, lwrite
      real, dimension(:), optional, intent(in) :: header
!
      character (len=fnlen) :: filename
      integer :: ia
!
      if (.not. lwrite .or. (nc <= 0)) return
!
      filename = trim(path) // '/' // trim(label) // 'averages.dat'
      if (lbinary) then
        open(lun_output, file=filename, form='unformatted', position='append')
        if (present (header)) write(lun_output) header
        write(lun_output) time
        if (label(1:1) == 'z') then
          write(lun_output) ( data(ia,:,:), ia=1, nc )
        else
          write(lun_output) data(:,:,1:nc)
        endif
        close(lun_output)
      else
        open(lun_output, file=filename, position='append')
        if (present (header)) write(lun_output,'(1p,8'//trim(fmt_avgs)//')') header
        write(lun_output,'(1pe12.5)') time
        if (label(1:1) == 'z') then
          write(lun_output,'(1p,8'//trim(fmt_avgs)//')') ( data(ia,:,:), ia=1, nc )
        else
          write(lun_output,'(1p,8'//trim(fmt_avgs)//')') data(:,:,1:nc)
        endif
        close(lun_output)
      endif
!
    endsubroutine output_average_2D
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
      use General, only: safe_character_append
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
      character (len=1024) :: labels
      integer :: pos
!     
      if (.not. lroot .or. (nc <= 0)) return
!
      filename = 'PHIAVG' // trim(number)
      open(lun_output, file=trim(path)//'/averages/'//trim(filename), form='unformatted', position='append')
      write(lun_output) nr, nzgrid, nc, nprocz
      write(lun_output) time, r, z(n1)+(/(pos*dz, pos=0, nzgrid-1)/), dr, dz
      ! note: due to passing data as implicit-size array,
      ! the indices (0:nz) are shifted to (1:nz+1),
      ! so that we have to write only the portion (2:nz+1).
      ! data has to be repacked to avoid writing an array temporary
      ! ngrs: writing an array temporary outputs corrupted data on copson
      write(lun_output) pack(data(:,2:nz+1,:,1:nc), .true.)
      labels = trim(name(1))
      if (nc >= 2) then
        do pos = 2, nc
          call safe_character_append (labels, ",", trim(name(pos)))
        enddo
      endif
      write(lun_output) len(labels), labels
      close(lun_output)
!
      ! append filename to list
      open (lun_output, FILE=trim(datadir)//'/averages/phiavg.files', POSITION='append')
      write (lun_output,'(A)') trim(filename)
      close (lun_output)
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
      use File_io, only: file_exists

      character (len=*), intent(in) :: path, plane
      integer, intent(in) :: ngrid, nname
!
      real, dimension(:), allocatable :: tmp
      character(len=fnlen) :: filename
      integer :: pos, num_rec, ioerr, alloc_err
      integer :: lun_input = 84
      real :: time
!
      if (.not. lroot) return
      if ((ngrid <= 0) .or. (nname <= 0)) return
      filename = trim(path) // '/' // trim(plane) // 'averages.dat'
      if (.not. file_exists (filename)) return
!
      allocate (tmp(ngrid * nname), stat = alloc_err)
      if (alloc_err > 0) call fatal_error('trim_average', 'Could not allocate memory for averages')
!
      if (lwrite_avg1d_binary) then
        open(lun_input, file=filename, form='unformatted', action='readwrite')
      else
        open(lun_input, file=filename, action='readwrite')
      endif
!
      ! count number of records
      num_rec = 0
      ioerr = 0
      time = 0.0
      do while ((t > time) .and. (ioerr == 0))
        num_rec = num_rec + 1
        if (lwrite_avg1d_binary) then
          read(lun_input, iostat=ioerr) time
          read(lun_input, iostat=ioerr) tmp
        else
          read(lun_input, *, iostat=ioerr) time
          read(lun_input, *, iostat=ioerr) tmp
        endif
      enddo
      if (time == t) num_rec = num_rec - 1
!
      if ((time >= t) .and. (ioerr == 0)) then
        ! trim excess data at the end of the average file
        rewind(lun_input)
        if (num_rec > 0) then
          do pos = 1, num_rec
            if (lwrite_avg1d_binary) then
              read(lun_input, iostat=ioerr) time
              read(lun_input, iostat=ioerr) tmp
            else
              read(lun_input, *, iostat=ioerr) time
              read(lun_input, *, iostat=ioerr) tmp
            endif
          enddo
        endif
        endfile (lun_input)
        if (ip <= 10) print *, 'trim_average: trimmed '//trim(plane)//'-averages for t >= ', t
      endif
!
      close (lun_input)
      deallocate (tmp)
!
    endsubroutine trim_average
!***********************************************************************
endmodule HDF5_IO
