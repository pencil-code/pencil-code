! $Id$
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   debug_io_hdf5.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  HDF5 debug-IO module
!
!  24-Jun-2019/PABourdin: coded HDF5 variant
!
module Debug_IO
!
  use Cdata
  use Cparam
  use HDF5_IO
  use File_io, only: parallel_file_exists
!
  implicit none
!
  public :: lun_input, lun_output
  public :: output, output_pencil
!
  interface output              ! Overload the `output' function
    module procedure output_vect
    module procedure output_scal
  endinterface
!
  interface output_pencil        ! Overload the `output_pencil' function
    module procedure output_pencil_vect
    module procedure output_pencil_scal
  endinterface
!
  ! define unique logical unit number for input and output calls
  integer :: lun_input=89
  integer :: lun_output=92
!
contains
!***********************************************************************
    subroutine output_vect(file,a,nv)
!
!  write debug snapshot file from a vector quantity
!
      use General, only: itoa
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv
      real, dimension (mx,my,mz,nv), intent(in) :: a
!
      character (len=fnlen) :: filename, dataset
      integer :: pos
!
      filename = trim (datadir_snap)//'/'//trim (file)//'.h5'
      call file_open_hdf5 (filename, truncate=(imn == 1))
      if (nv == 1) then
        call output_hdf5 ('data', a(:,:,:,1))
      elseif (nv == 3) then
        call output_hdf5 ('data_x', a(:,:,:,1))
        call output_hdf5 ('data_y', a(:,:,:,2))
        call output_hdf5 ('data_z', a(:,:,:,3))
      else
        do pos = 1, nv
          dataset = 'data_'//trim (itoa (pos))
          call output_hdf5 (dataset, a(:,:,:,pos))
        enddo
      endif
      call file_close_hdf5
!
      if (lroot .and. (imn == 1)) then
        call file_open_hdf5 (filename, truncate=.false., global=.false.)
        call output_hdf5 ('time', real(t))
        call output_hdf5 ('num_components', nv)
        call create_group_hdf5 ('settings')
        call output_hdf5 ('settings/nprocx', nprocx)
        call output_hdf5 ('settings/nprocy', nprocy)
        call output_hdf5 ('settings/nprocz', nprocz)
        call file_close_hdf5
      endif
!
    endsubroutine output_vect
!***********************************************************************
    subroutine output_scal(file,a,nv)
!
!  write debug snapshot file from a scalar quantity
!
!
      integer, intent(in), optional :: nv
      real, dimension (mx,my,mz), intent(in) :: a
      character (len=*), intent(in) :: file
!
      call output_vect (file, (/ a /), 1)
!
    endsubroutine output_scal
!***********************************************************************
    subroutine output_pencil_vect(file,a,nv)
!
!  write debug snapshot file of penciled vector data
!
      use General, only: itoa
!
      integer, intent(in) :: nv
      real, dimension (nx,nv), intent(in) :: a
      character (len=*), intent(in) :: file
!
      character (len=fnlen) :: filename, dataset
      integer :: pos
!
      if (lroot .and. headt .and. (imn == 1)) &
          print *, 'output_pencil: Writing to ' // trim(file) // ' for debugging -- this may slow things down'
!
      filename = trim (datadir_snap)//'/'//trim (file)//'.h5'
      call file_open_hdf5 (filename, truncate=(imn == 1))
      if (nv == 1) then
        call output_hdf5 ('data', a(:,1), mm(imn)-m1, nn(imn)-n1)
      elseif (nv == 3) then
        call output_hdf5 ('data_x', a(:,1), mm(imn)-m1, nn(imn)-n1)
        call output_hdf5 ('data_y', a(:,2), mm(imn)-m1, nn(imn)-n1)
        call output_hdf5 ('data_z', a(:,3), mm(imn)-m1, nn(imn)-n1)
      else
        do pos = 1, nv
          dataset = 'data_'//trim (itoa (pos))
          call output_hdf5 (dataset, a(:,pos), mm(imn)-m1, nn(imn)-n1)
        enddo
      endif
      call file_close_hdf5
!
      if (lroot .and. (imn == 1)) then
        call file_open_hdf5 (filename, truncate=.false., global=.false.)
        call output_hdf5 ('time', real(t))
        call output_hdf5 ('num_components', nv)
        call create_group_hdf5 ('settings')
        call output_hdf5 ('settings/nprocx', nprocx)
        call output_hdf5 ('settings/nprocy', nprocy)
        call output_hdf5 ('settings/nprocz', nprocz)
        call file_close_hdf5
      endif
!
    endsubroutine output_pencil_vect
!***********************************************************************
    subroutine output_pencil_scal(file,a,nv)
!
!  write debug snapshot file of penciled scalar data
!
      integer, intent(in), optional :: nv
      real, dimension (nx), intent(in) :: a
      character (len=*), intent(in) :: file
!
      call output_pencil_vect (file, (/ a /), 1)
!
    endsubroutine output_pencil_scal
!***********************************************************************
endmodule Debug_IO
