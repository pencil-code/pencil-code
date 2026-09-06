! $Id$
!
!  Module for boundary conditions. Extracted from (no)mpicomm, since
!  all non-periodic (external) boundary conditions require the same
!  code for serial and parallel runs.
!
module Boundcond
!
  use Cdata
  use Mpicomm
!
  implicit none
!
  private
!
  public :: driving
!
  real, dimension (mcom) :: tau_inv=0.0
  integer, dimension (mcom) :: target_proc_x=-1, target_proc_y=-1, target_proc_z=-1
!
  type :: data_array
    real, pointer, dimension(:,:) :: ptr => null()
  end type data_array
  type(data_array), dimension(mcom) :: data_slots_xy => null(), data_slots_xz => null(), data_slots_yz => null()
!
!  Run parameters.
!
  character (len=fnlen), dimension(mcom) :: driver_xy="", driver_xz="", driver_yz=""
  logical dimension(mcom) :: ldrive_xy=.false., ldrive_xz=.false., ldrive_yz=.false.
  integer, dimension (mcom) :: driver_pos_x=0, driver_pos_y=0, driver_pos_z=0
  real, dimension (mcom) :: decay_time=0.0
!
  namelist /driver_run_pars/ &
      driver_xy, driver_xz, driver_yz, &
      driver_pos_x, driver_pos_y, driver_pos_z, &
      decay_time
!
  contains
!***********************************************************************
    subroutine initialize_driver
!
!  Initialize the driver.
!
! 05-Sep-2026/PABourdin: coded
!
      use Messages, only: svn_id
!
      integer :: f_index
!
!  Identify version number (generated automatically by SVN).
!
      call svn_id( &
           "$Id$")
!
      tau_inv(:) = decay_time(:)
!
      target_proc_x = (driver_pos_x-1) / nx
      target_proc_y = (driver_pos_y-1) / ny
      target_proc_z = (driver_pos_z-1) / nz
!
      ldrive_xy = (driver_xy(f_index) /= "") .and. (target_proc_z == ipz)
      ldrive_xz = (driver_xz(f_index) /= "") .and. (target_proc_y == ipy)
      ldrive_yz = (driver_yz(f_index) /= "") .and. (target_proc_x == ipx)
!
      do f_index = 1, mcom
!
        if (ldrive_xy) then
          if (not associated (data_slots_xy(f_index)%ptr)) then
            allocate (data_slots_xy(f_index)%ptr(nx,ny))
          endif
        endif
!
        if (ldrive_xz) then
          if (not associated (data_slots_xz(f_index)%ptr)) then
            allocate (data_slots_xz(f_index)%ptr(nx,nz))
          endif
        endif
!
        if (ldrive_yz) then
          if (not associated (data_slots_yz(f_index)%ptr)) then
            allocate (data_slots_yz(f_index)%ptr(ny,nz))
          endif
        endif
!
      enddo
!
    endsubroutine initialize_driver
!***********************************************************************
    subroutine finalize_driver
!
!  Initialize the driver.
!
! 05-Sep-2026/PABourdin: coded
!
      integer :: f_index
!
      do f_index = 1, mcom
        if (associated (data_slots_xy(f_index)%ptr)) deallocate (data_slots_xy(f_index)%ptr)
        if (associated (data_slots_xz(f_index)%ptr)) deallocate (data_slots_xz(f_index)%ptr)
        if (associated (data_slots_yz(f_index)%ptr)) deallocate (data_slots_yz(f_index)%ptr)
      enddo
!
    endsubroutine finalize_driver
!***********************************************************************
    subroutine read_driver_run_pars(iomsg)
!
      use File_io, only: parallel_unit
!
      character(LEN=*), intent(out) :: iomsg
      integer :: iostat
!
      read (parallel_unit, NML=driver_run_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat == 0) iomsg = ""
!
    endsubroutine read_driver_run_pars
!***********************************************************************
    subroutine write_driver_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write (unit, NML=driver_run_pars)
!
    endsubroutine write_driver_run_pars
!***********************************************************************
    subroutine find_frame (time, filename, frame_type, frame_pos, frame_time, plane, lreader)
!
!  Finds the position of the frame before/at (l) or after (r) the given time.
!  If a frame matches 'time', this frame is considered to be (l).
!  The time table is expected to be in SI units, not using F77 record markers.
!
!  'time' is the reference time that is being searched for.
!  'filename' is the time table file name.
!  'frame_type' indicates the desired frame type (l) or (r).
!  'frame_pos' is set to the position (record number) of the desired frame.
!  'frame_time' is set to the time of the corresponding frame.
!
!  06-Sep-2026/PABourdin: adapted from the solar_corona module
!
      use File_io, only: file_exists
      use Mpicomm, only: distribute_xy, distribute_xz, distribute_yz
!
      real, intent(in) :: time
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: frame_type
      integer, intent(out) :: frame_pos
      real, intent(out) :: frame_time
      character(len=2), intent(in) :: plane
      logical, intent(in) :: lreader
!
      integer :: rec_len, io_error
      real :: time_l, delta_t
      integer, parameter :: unit=17
!
      if ((.not. ldrive_xy) .and. (.not. ldrive_xz) .and. (.not. ldrive_yz)) return
!
      if (lreader) then
        ! read time and frame number from "_times.dat" file, if it exists
        if (.not. file_exists (filename)) then
          ! No time series => use only first frame, forever
          call warning ('driver', '"'//trim (filename)//'" not found, using only first frame, forever!')
          frame_pos = 1
          frame_time = huge (0.0)
        else
          ! Read the time table from file
          inquire (iolength=rec_len) time
          open (unit, file=filename, form='unformatted', recl=rec_len, access='direct')
!
          io_error = 0
          delta_t = 0.0
          frame_pos = 1
          read (unit, rec=frame_pos) time_l
          time_l = time_l / unit_time
          if (time_l < 0.0) call fatal_error ('find_frame', trim (filename)//' first frame must be >= 0.', .true.)
          if (time < time_l) then
            ! 'time' is still before the first frame
            ! => set following (r) frame to point to the first frame
            io_error = -1
            frame_time = time_l
            time_l = 0.0
          endif
!
          do while (io_error == 0)
            frame_pos = frame_pos + 1
            read (unit, rec=frame_pos, iostat=io_error) frame_time
            if (io_error == 0) then
              frame_time = frame_time / unit_time + delta_t
              ! Test if correct time step has been reached
              if ((time >= time_l) .and. (time < frame_time)) exit
              ! If not, continue searching...
              time_l = frame_time
            else
              ! There was an error while reading, check why
              if (frame_pos <= 2) call fatal_error ('find_frame', &
                  trim (filename)//' contains less than two frames.', .true.)
              if (time_l <= 0.0) call fatal_error ('find_frame', &
                  trim (filename)//' last frame must have time > 0.', .true.)
              ! EOF reached => read from beginning
              delta_t = delta_t + time_l
              frame_pos = 0
              io_error = 0
            endif
          enddo
          close (unit)
!
          if (frame_type == 'l') then
            frame_pos = frame_pos - 1
            frame_time = time_l
          endif
        endif
      endif
!
      if (plane == "xy") then
        ! Distribute results in the xy-plane
        if (lreader) then
          call distribute_xy (frame_pos, frame_pos)
          call distribute_xy (frame_time, frame_time)
        else
          call distribute_xy (frame_pos)
          call distribute_xy (frame_time)
        endif
      elseif (plane == "xz") then
        ! Distribute results in the xz-plane
        if (lreader) then
          call distribute_xz (frame_pos, frame_pos)
          call distribute_xz (frame_time, frame_time)
        else
          call distribute_xz (frame_pos)
          call distribute_xz (frame_time)
        endif
      elseif (plane == "yz") then
        ! Distribute results in the yz-plane
        if (lreader) then
          call distribute_yz (frame_pos, frame_pos)
          call distribute_yz (frame_time, frame_time)
        else
          call distribute_yz (frame_pos)
          call distribute_yz (frame_time)
        endif
      endif
!
    endsubroutine find_frame
!***********************************************************************
    subroutine driver_update(f, df, f_index)
!
!  Update the driving data, if needed, including time-interpolation.
!
! 05-Sep-2026/PABourdin: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      integer, intent(in) :: f_index, grid_pos
!
      if (ldrive_xy) then
        ! check if driver data needs to be updated from file
        
        ! interpolate driver data in time
        
      endif
!
      if (ldrive_xz) then
        ! check if driver data needs to be updated from file
        
        ! interpolate driver data in time
        
      endif
!
      if (ldrive_yz) then
        ! check if driver data needs to be updated from file
        
        ! interpolate driver data in time
        
      endif
!
    endsubroutine driver_update
!***********************************************************************
    subroutine driver_apply(f, df, f_index, grid_pos)
!
!  Apply the driving in the specified f-array component at specified positions.
!
! 05-Sep-2026/PABourdin: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      integer, intent(in) :: f_index, grid_pos
!
      integer, save :: px
!
      if (ldrive_xy) then
        ! apply driving in xy-plane at desired z position
        df(l1:l2,m,n,f_index) = df(l1:l2,m,n,f_index) - tau_inv * (f(l1:l2,m,n,f_index) - data_xy(:,m-nghost))
      endif
!
      if (ldrive_xz) then
        ! apply driving in xz-plane at desired y position
        df(l1:l2,m,n,f_index) = df(l1:l2,m,n,f_index) - tau_inv * (f(l1:l2,m,n,f_index) - data_xz(:,n-nghost))
      endif
!
      if (ldrive_yz) then
        ! apply driving in yz-plane at desired x position
        pos_l = driver_pos_x(f_index) + nghost
        df(pos_l,m,n,f_index) = df(pos_l,m,n,f_index) - tau_inv * (f(pos_l,m,n,f_index) - data_yz(m-nghost,n-nghost))
      endif
!
    endsubroutine driver_apply
!***********************************************************************

