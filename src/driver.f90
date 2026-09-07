! $Id$
!
!  Module for boundary driving from external files for any f-array component(s).
!
module Driver
!
  use Cdata
  use Mpicomm
!
  implicit none
!
  private
!
  public :: driver_apply
!
  real, dimension (mcom) :: tau_inv=0.0
  integer, dimension (mcom) :: target_proc_x=-1, target_proc_y=-1, target_proc_z=-1
!
  type :: data_array
    real, pointer, dimension(:,:) :: frame => null()
    real, pointer, dimension(:,:) :: frame_l => null()
    real, pointer, dimension(:,:) :: frame_r => null()
    real :: time_l, time_r
  end type data_array
  type(data_array), dimension(mcom) :: data_xy, data_xz, data_yz
!
!  Run parameters.
!
  character (len=fnlen), dimension(mcom) :: driver_xy="", driver_xz="", driver_yz=""
  logical, dimension(mcom) :: ldrive_xy=.false., ldrive_xz=.false., ldrive_yz=.false.
  integer, dimension (mcom) :: driver_pos_x=0, driver_pos_y=0, driver_pos_z=0
  real, dimension (mcom) :: data_unit=0.0, decay_time=0.0, time_offset=0.0
!
  namelist /driver_run_pars/ &
      driver_xy, driver_xz, driver_yz, &
      driver_pos_x, driver_pos_y, driver_pos_z, &
      data_unit, decay_time, time_offset
!
  contains
!***********************************************************************
    subroutine initialize_driver
!
!  Initialize the driver.
!
! 06-Sep-2026/PABourdin: coded
!
      use Messages, only: svn_id, fatal_error
!
      integer :: f_index, alloc_err
!
!  Identify version number (generated automatically by SVN).
!
      call svn_id( &
           "$Id$")
!
      tau_inv(:) = 1.0 / decay_time(:)
!
      target_proc_x(:) = (driver_pos_x(:)-1) / nx
      target_proc_y(:) = (driver_pos_y(:)-1) / ny
      target_proc_z(:) = (driver_pos_z(:)-1) / nz
!
      do f_index = 1, mcom
        ldrive_xy(f_index) = (driver_xy(f_index) /= "") .and. (target_proc_z(f_index) == ipz)
        ldrive_xz(f_index) = (driver_xz(f_index) /= "") .and. (target_proc_y(f_index) == ipy)
        ldrive_yz(f_index) = (driver_yz(f_index) /= "") .and. (target_proc_x(f_index) == ipx)
!
        if (ldrive_xy(f_index) .and. (data_unit(f_index) /= 0.0)) &
            call fatal_error ('initialize_driver', "Trying to use driving without setting the corresponding 'data_unit'.", .true.)
!
        if (ldrive_xy(f_index)) then
          if (.not. associated (data_xy(f_index)%frame)) then
            allocate (data_xy(f_index)%frame(nx,ny), data_xy(f_index)%frame_l(nx,ny), &
                data_xy(f_index)%frame_r(nx,ny), stat=alloc_err)
            if (alloc_err > 0) call fatal_error ('initialize_driver', 'Could not allocate "data_xy".', .true.)
          endif
        endif
!
        if (ldrive_xz(f_index)) then
          if (.not. associated (data_xz(f_index)%frame)) then
            allocate (data_xz(f_index)%frame(nx,nz), data_xz(f_index)%frame_l(nx,nz), &
                data_xz(f_index)%frame_r(nx,nz), stat=alloc_err)
            if (alloc_err > 0) call fatal_error ('initialize_driver', 'Could not allocate "data_xz".', .true.)
          endif
        endif
!
        if (ldrive_yz(f_index)) then
          if (.not. associated (data_yz(f_index)%frame)) then
            allocate (data_yz(f_index)%frame(ny,nz), data_yz(f_index)%frame_l(ny,nz), &
                data_yz(f_index)%frame_r(ny,nz), stat=alloc_err)
            if (alloc_err > 0) call fatal_error ('initialize_driver', 'Could not allocate "data_yz".', .true.)
          endif
        endif
      enddo
!
    endsubroutine initialize_driver
!***********************************************************************
    subroutine finalize_driver
!
!  Finalize the driver.
!
! 06-Sep-2026/PABourdin: coded
!
      integer :: f_index
!
      do f_index = 1, mcom
        if (associated (data_xy(f_index)%frame)) &
            deallocate (data_xy(f_index)%frame, data_xy(f_index)%frame_l, data_xy(f_index)%frame_r)
        if (associated (data_xz(f_index)%frame)) &
            deallocate (data_xz(f_index)%frame, data_xz(f_index)%frame_l, data_xz(f_index)%frame_r)
        if (associated (data_yz(f_index)%frame)) &
            deallocate (data_yz(f_index)%frame, data_yz(f_index)%frame_l, data_yz(f_index)%frame_r)
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
    subroutine interpolate_time_2D (time, time_l, time_r, data_l, data_r, data)
!
!  Interpolate 2D data frame in time.
!
!  06-Sep-2026/PABourdin: adapted from the "solar_corona" module
!
      real, intent(in) :: time, time_l, time_r
      real, dimension(:,:), intent(in) :: data_l, data_r
      real, dimension(:,:), intent(inout) :: data
!
      real :: factor
!
      if (time <= time_l) then
        data = data_l
      elseif (time >= time_r) then
        data = data_r
      else
        ! Interpolate data
        factor = (time - time_l) / (time_r - time_l)
        data = data_l * (1.0 - factor) + data_r * factor
      endif
!
    endsubroutine interpolate_time_2D
!***********************************************************************
    subroutine update_frame (time, time_l, time_r, f_index, plane, dim_1, dim_2, frame_l, frame_r, data_local)
!
!  Check if an update of the data frame is needed and load frame from file.
!  An interpolated data frame will be added to the given local frame.
!  The previous frame (l) and following frame (r) are updated.
!
!  06-Sep-2026/PABourdin: adapted from the "solar_corona" module
!
      real, intent(in) :: time
      real, intent(inout) :: time_l, time_r
      integer, intent(in) :: f_index
      character(len=2), intent(in) :: plane
      integer, intent(in) :: dim_1, dim_2
      real, dimension(dim_1,dim_2), intent(inout) :: frame_l, frame_r, data_local
!
      character (len=fnlen), save :: filename=""
      integer :: pos_l, pos_r
      logical :: lreader
      logical, save :: lfirst_call=.true.
!
      if (plane == "xy") then
        lreader = ldrive_xy(f_index) .and. lfirst_proc_xy
      elseif (plane == "xz") then
        lreader = ldrive_xz(f_index) .and. lfirst_proc_xz
      elseif (plane == "yz") then
        lreader = ldrive_yz(f_index) .and. lfirst_proc_yz
      endif
!
      if (lfirst_call) then
        ! Load previous (l) frame and store it in (r), will be shifted later
        filename = trim (driver_xy(f_index))//"_"//trim (plane)//"_times.dat"
        call find_frame (time, filename, 'l', pos_l, time_l, f_index, plane, lreader)
        if (pos_l == 0) then
          ! The simulation started before the first frame of the time series
          ! start from zero velocities
          frame_r = 0.0
          time_l = -time_offset(f_index)
        else
          filename = trim (driver_xy(f_index))//"_"//trim (plane)//".dat"
          call read_frame (pos_l, filename, f_index, plane, dim_1, dim_2, frame_r, lreader)
        endif
        ! Make sure that the following (r) frame will get loaded:
        time_r = time_l
        lfirst_call = .false.
      endif
!
      if (time >= time_r) then
        ! Shift data from following (r) to previous (l) frame
        frame_l = frame_r
        time_l = time_r
        ! Read new following (r) frame
        filename = trim (driver_xy(f_index))//"_"//trim (plane)//"_times.dat"
        call find_frame (time, filename, 'r', pos_r, time_r, f_index, plane, lreader)
        filename = trim (driver_xy(f_index))//"_"//trim (plane)//".dat"
        call read_frame (pos_r, filename, f_index, plane, dim_1, dim_2, frame_r, lreader)
      endif
!
    endsubroutine update_frame
!***********************************************************************
    subroutine read_frame (frame, filename, f_index, plane, dim_1, dim_2, data, lreader)
!
!  Reads one data frame from a given file at a given frame position
!  and distributes the results in the respective plane.
!  The data is expected to be in SI units, not using F77 record markers.
!
!  06-Sep-2026/PABourdin: adapted from the "solar_corona" module
!
      use Messages, only: fatal_error
      use Mpicomm, only: distribute_xy, distribute_xz, distribute_yz
!
      integer, intent(in) :: frame
      character(len=*), intent(in) :: filename
      integer, intent(in) :: f_index
      character(len=2), intent(in) :: plane
      integer, intent(in) :: dim_1, dim_2
      real, dimension(dim_1,dim_2), intent(out) :: data
      logical, intent(in) :: lreader
!
      integer, parameter :: unit=12
      real, dimension(:,:), allocatable :: buffer
      integer :: buffer_dim_1, buffer_dim_2
      integer :: rec_len, alloc_err
!
      if (lreader) then
        if (plane == "xy") then
          buffer_dim_1 = dim_1 * nprocx
          buffer_dim_2 = dim_2 * nprocy
        elseif (plane == "xz") then
          buffer_dim_1 = dim_1 * nprocx
          buffer_dim_2 = dim_2 * nprocz
        elseif (plane == "yz") then
          buffer_dim_1 = dim_1 * nprocy
          buffer_dim_2 = dim_2 * nprocz
        endif
        allocate (buffer(buffer_dim_1,buffer_dim_2), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_frame', 'Could not allocate buffer.', .true.)
!
        ! read data frame from file
        inquire (IOLENGTH=rec_len) 1.0d0
        rec_len = rec_len * dim_1 * dim_2
        open (unit, file=filename, form='unformatted', recl=rec_len, access='direct')
        read (unit, rec=frame) buffer
!
        ! distribute data along plane
        if (plane == "xy") then
          call distribute_xy (data, buffer)
        elseif (plane == "xz") then
          call distribute_xz (data, buffer)
        elseif (plane == "yz") then
          call distribute_yz (data, buffer)
        endif
        deallocate (buffer)
      else
        ! receive local portion of data
        if (plane == "xy") then
          call distribute_xy (data)
        elseif (plane == "xz") then
          call distribute_xz (data)
        elseif (plane == "yz") then
          call distribute_yz (data)
        endif
      endif
!
      ! convert SI to PC units
      data = data / data_unit(f_index)
!
    endsubroutine read_frame
!***********************************************************************
    subroutine find_frame (time, filename, frame_type, frame_pos, frame_time, f_index, plane, lreader)
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
!  06-Sep-2026/PABourdin: adapted from the "solar_corona" module
!
      use File_io, only: file_exists
      use Messages, only: warning, fatal_error
      use Mpicomm, only: distribute_xy, distribute_xz, distribute_yz
!
      real, intent(in) :: time
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: frame_type
      integer, intent(out) :: frame_pos
      real, intent(out) :: frame_time
      integer, intent(in) :: f_index
      character(len=2), intent(in) :: plane
      logical, intent(in) :: lreader
!
      integer :: out_int, rec_len, io_error
      real :: out_real, time_l, delta_t
      integer, parameter :: unit=17
!
      if ((.not. ldrive_xy(f_index)) .and. (plane == "xy")) return
      if ((.not. ldrive_xz(f_index)) .and. (plane == "xz")) return
      if ((.not. ldrive_yz(f_index)) .and. (plane == "yz")) return
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
          call distribute_xy (out_int, frame_pos)
          call distribute_xy (out_real, frame_time)
        else
          call distribute_xy (frame_pos)
          call distribute_xy (frame_time)
        endif
      elseif (plane == "xz") then
        ! Distribute results in the xz-plane
        if (lreader) then
          call distribute_xz (out_int, frame_pos)
          call distribute_xz (out_real, frame_time)
        else
          call distribute_xz (frame_pos)
          call distribute_xz (frame_time)
        endif
      elseif (plane == "yz") then
        ! Distribute results in the yz-plane
        if (lreader) then
          call distribute_yz (out_int, frame_pos)
          call distribute_yz (out_real, frame_time)
        else
          call distribute_yz (frame_pos)
          call distribute_yz (frame_time)
        endif
      endif
!
    endsubroutine find_frame
!***********************************************************************
    subroutine driver_update(f_index)
!
!  Update the driving data, if needed, including time-interpolation.
!
! 05-Sep-2026/PABourdin: coded
!
      integer, intent(in) :: f_index
!
      real :: time
      real, dimension(:,:), pointer :: frame_l, frame_r, data_local
!
      time = t - time_offset(f_index)
!
      if (ldrive_xy(f_index)) then
        frame_l = data_xy(f_index)%frame_l
        frame_r = data_xy(f_index)%frame_r
        ! check if driver data needs to be updated from file
        call update_frame (time, data_xy(f_index)%time_l, data_xy(f_index)%time_r, f_index, &
            "xy", nx, ny, frame_l, frame_r, data_xy(f_index)%frame)
        ! interpolate driver data in time
        call interpolate_time_2D (time, data_xy(f_index)%time_l, data_xy(f_index)%time_r, frame_l, frame_r, data_xy(f_index)%frame)
      endif
!
      if (ldrive_xz(f_index)) then
        frame_l = data_xz(f_index)%frame_l
        frame_r = data_xz(f_index)%frame_r
        ! check if driver data needs to be updated from file
        call update_frame (time, data_xz(f_index)%time_l, data_xz(f_index)%time_r, f_index, &
            "xz", nx, nz, frame_l, frame_r, data_xz(f_index)%frame)
        ! interpolate driver data in time
        call interpolate_time_2D (time, data_yz(f_index)%time_l, data_yz(f_index)%time_r, frame_l, frame_r, data_xz(f_index)%frame)
      endif
!
      if (ldrive_yz(f_index)) then
        frame_l = data_yz(f_index)%frame_l
        frame_r = data_yz(f_index)%frame_r
        ! check if driver data needs to be updated from file
        call update_frame (time, data_yz(f_index)%time_l, data_yz(f_index)%time_r, f_index, &
            "yz", ny, nz, frame_l, frame_r, data_yz(f_index)%frame)
        ! interpolate driver data in time
        call interpolate_time_2D (time, data_yz(f_index)%time_l, data_yz(f_index)%time_r, frame_l, frame_r, data_yz(f_index)%frame)
      endif
!
    endsubroutine driver_update
!***********************************************************************
    subroutine driver_apply(f, df)
!
!  Apply the driving in the specified f-array component at specified positions.
!
! 05-Sep-2026/PABourdin: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
!
      integer :: l, m, n
      integer :: f_index
      real :: tau_invers
!
      do f_index = 1, mcom
        call driver_update(f_index)
        tau_invers = tau_inv(f_index)
!
        if (ldrive_xy(f_index)) then
          ! apply driving in xy-plane at desired z position
          n = driver_pos_z(f_index) + nghost
          df(l1:l2,m1:m2,n,f_index) = df(l1:l2,m1:m2,n,f_index) - tau_invers * (f(l1:l2,m1:m2,n,f_index) - data_xy(f_index)%frame)
        endif
!
        if (ldrive_xz(f_index)) then
          ! apply driving in xz-plane at desired y position
          m = driver_pos_y(f_index) + nghost
          df(l1:l2,m,n1:n2,f_index) = df(l1:l2,m,n1:n2,f_index) - tau_invers * (f(l1:l2,m,n1:n2,f_index) - data_xz(f_index)%frame)
        endif
!
        if (ldrive_yz(f_index)) then
          ! apply driving in yz-plane at desired x position
          l = driver_pos_x(f_index) + nghost
          df(l,m1:m2,n1:n2,f_index) = df(l,m1:m2,n1:n2,f_index) - tau_invers * (f(l,m1:m2,n1:n2,f_index) - data_yz(f_index)%frame)
        endif
      enddo
!
    endsubroutine driver_apply
!***********************************************************************
endmodule Driver
