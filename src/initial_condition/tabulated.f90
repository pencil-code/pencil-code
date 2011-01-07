! $Id$
!
!  This module loads an existing snapshot and interpolates it to the current
!  grid as the initial condition.  The snapshot should be saved as a global
!  snapshot in a single file (default name: tabulated.dat) under the run
!  directory.  The user should specify its dimensions using nxtab, nytab, and nztab.
!  The number of variables in the snapshot must be consistent with the current
!  mfarray.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../initial_condition.h'
!
  character(len=256) :: file = 'tabulated.dat'
  integer :: nxtab = nxgrid, nytab = 1, nztab = 1
!
  namelist /initial_condition_pars/ nxtab, nytab, nztab, file
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f)
!
!  Initializes all the f arrays in one call.
!
!  21-dec-10/ccyang: coded
!
      use General, only: spline
      use Mpicomm
!
      real, dimension (mx,my,mz,mfarray), intent(out) :: f
!
      real, dimension(:,:,:,:), allocatable :: ftab
      real, dimension(:), allocatable :: xtab, ytab, ztab
      real, dimension(:), allocatable :: fint
      integer :: mxtab, mytab, mztab
      real :: tt, dxt, dyt, dzt, deltayt
!
      integer, parameter :: lun = 1
      integer :: status
      integer :: i, m
!
!  Allocate memory for read.
!
      mxtab = 2 * nghost + nxtab
      mytab = 2 * nghost + nytab
      mztab = 2 * nghost + nztab
      allocate (ftab(mxtab,mytab,mztab,mfarray), xtab(mxtab), ytab(mytab), ztab(mztab), fint(mx), stat=status)
      if (status /= 0) call fatal_error('initialize_initial_condition', 'failed to allocate memory')
!
      if (lroot) then
!
!  Read the snapshot.
!
        open (unit=lun, file=trim(file), form='unformatted', status='old', action='read', iostat=status)
        if (status /= 0) call fatal_error('initialize_initial_condition', 'failed to open the snapshot file')
        read(lun, iostat=status) ftab
        if (status /= 0) call fatal_error('initialize_initial_condition', 'failed to read the f arrays')
        if (lshear) then
          read(lun, iostat=status) tt, xtab, ytab, ztab, dxt, dyt, dzt, deltayt
        else
          read(lun, iostat=status) tt, xtab, ytab, ztab, dxt, dyt, dzt
        endif
        if (status /= 0) call fatal_error('initialize_initial_condition', 'failed to read the time and grid information')
        close (unit=lun)
      endif
!
!  Broadcast the data.
!
      call mpibcast_real(ftab, (/ mxtab, mytab, mztab, mfarray /))
      call mpibcast_real(xtab, mxtab)
      call mpibcast_real(ytab, mytab)
      call mpibcast_real(ztab, mztab)
!
!  Initialize the f arrays by interpolation.
!
      if (nytab /= 1 .or. nztab /= 1) call fatal_error('initialize_initial_condition', 'only 1D table in x is implemented.')
      f(:,:,:,irho) = 0.
      m = nghost + 1
      do i = 1, mfarray
        call spline(xtab, ftab(:,m,m,i), x, fint, mxtab, mx)
        f(:,:,:,i) = f(:,:,:,i) + spread(spread(fint,2,my),3,mz)
      enddo
!
!  Deallocate the memory.
!
      deallocate (ftab, xtab, ytab, ztab, fint, stat=status)
      if (status /= 0) call warning('initialize_initial_condition','failed to deallocate memory')
!
    endsubroutine initial_condition_all
!***********************************************************************
    subroutine read_initial_condition_pars(unit,iostat)
!
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=initial_condition_pars,IOSTAT=iostat)
      else
        read(unit,NML=initial_condition_pars)
      endif
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
!
      write(unit,NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
