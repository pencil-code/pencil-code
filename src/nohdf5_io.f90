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
    subroutine initialize_hdf5
!
      ! nothing to do
!
    endsubroutine initialize_hdf5
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
    subroutine output_hdf5_double_0D(name, data)
!
      character (len=*), intent(in) :: name
      double precision, intent(in) :: data
!
      call fatal_error ('output_hdf5_double_0D', 'You can not use HDF5 without setting an HDF5_IO module.')
      call keep_compiler_quiet(name)
      call keep_compiler_quiet(data)
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
      call keep_compiler_quiet(data)
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
    subroutine index_append(varname,ivar,vector,array)
!
! 14-oct-18/PAB: coded
!
      character (len=*), intent(in) :: varname
      integer, intent(in) :: ivar
      integer, intent(in), optional :: vector
      integer, intent(in), optional :: array
!
      if (lroot) then
        open(lun_output,file=trim(datadir)//'/'//trim(index_pro), POSITION='append')
! Axel, please uncomment for debugging on Beskow:
! write (123,*) '============> append to "index.pro" ivar,varname=', ivar, varname
        if (present (vector) .and. present (array)) then
          ! expand array: iuud => indgen(vector)
          write(lun_output,*) trim(varname)//'=indgen('//trim(itoa(array))//')*'//trim(itoa(vector))//'+'//trim(itoa(ivar))
        else
          write(lun_output,*) trim(varname)//'='//trim(itoa(ivar))
        endif
        close(lun_output)
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
endmodule HDF5_IO
