! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
!   n_pivot = number of pivot points
!   udrive = strength of the forcing velocity
!   xp = pivot point
!   yp = pivot point
!   rad = distance of the footpoint to the pivot point
!   width = width of the footpoint
!   lam_twist = Gaussian RMS width
!   lam_u = coefficient for the exponential saturation for the velocity
!   r_profile = velocity profile in the radial direction
!   z_profile = velocity profile along the z-direction
!   lam_z = coefficient for the z-profiles
!   blinking = flag for the 'vortex blinking'
!   delay_blink = waiting time for the blinking process to commence
!   t_blink_up = time period for the driver to be positive
!   t_blink_0u = time period for the driver to be zero after positive driving
!   t_blink_down = time period for the driver to be negative
!   t_blink_0d= time period for the driver to be zero after negative driving
!
  integer :: n_pivot = 1
  real, dimension(6) :: udrive = (/0.1,0.,0.,0.,0.,0./)
  real, dimension(6) :: xp = (/0.,0.,0.,0.,0.,0./)
  real, dimension(6) :: yp = (/0.,0.,0.,0.,0.,0./)
  real, dimension(6) :: rad = (/2.,0.,0.,0.,0.,0./)
  real, dimension(6) :: width = (/1.,0.,0.,0.,0.,0./)
  real, dimension(6) :: lam_twist = (/0.17,0.,0.,0.,0.,0./)
  real :: lam_u = 1
  character (len=labellen) :: r_profile = 'gaussian'
  character (len=labellen) :: z_profile = 'exp'
  real :: lam_z = 1
  logical, save :: lblink = .False.
  real, dimension(6) :: delay_blink = (/0.,0.,0.,0.,0.,0./)
  real, dimension(6) :: t_blink_up = (/1.,0.,0.,0.,0.,0./)
  real, dimension(6) :: t_blink_0u = (/0.,0.,0.,0.,0.,0./)
  real, dimension(6) :: t_blink_down = (/0.,0.,0.,0.,0.,0./)
  real, dimension(6) :: t_blink_0d = (/0.,0.,0.,0.,0.,0./)
!    
! input parameters
  namelist /special_init_pars/ n_pivot
!
! run parameters
  namelist /special_run_pars/ &
    n_pivot, udrive, xp, yp, rad, width, lam_twist, lam_u, r_profile, z_profile, lam_z, &
    lblink, delay_blink, t_blink_up, t_blink_0u, t_blink_down, t_blink_0d
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
! Called with lreloading indicating a RELOAD
!
!  07-may-2015/iomsn (Simon Candelaresi): coded
!
      use Mpicomm, only: parallel_file_exists
!
      real, dimension(mx,my,mz,mfarray) :: f
!
! Consistency checks:
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  Initialize special condition; called by start.f90.
!
!  07-may-2015/iomsn (Simon Candelaresi): coded
!      
      real, dimension(mx,my,mz,mfarray), intent(out) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine finalize_special(f)
!
!  Called right before exiting.
!
!  07-may-2015/iomsn (Simon Candelaresi): coded
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine finalize_special
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      include '../unit.h'
      integer, intent(inout), optional :: iostat
!
      read (unit, NML=special_init_pars)
      if (present (iostat)) iostat = 0
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write (unit,NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
      include '../unit.h'
      integer, intent(inout), optional :: iostat
!
      read (unit, NML=special_run_pars)
      if (present (iostat)) iostat = 0
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write (unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!  07-may-2015/iomsn (Simon Candelaresi): coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      ! Apply driving of velocity field.
      if (.not. lpencil_check_at_work) &
          call vel_driver (f, df)
!
    endsubroutine special_calc_hydro
!***********************************************************************
  subroutine vel_driver (f, df)
!
! Drive bottom boundary horizontal velocities with given velocity field.
!
!  07-may-2015/iomsn (Simon Candelaresi): coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      integer :: ip
      real :: offset                    ! offset for u = 0 at the edge of the foopoint
      real, dimension(mx) :: dist       ! distance of point to pivot point
      real, dimension(mx) :: ux, uy     ! velocity in x and y direction      
      real :: z_factor                  ! multiplication factor for the z-dependence
      integer :: blink                  ! either -1, 0 or 1, depending on the blinking stage
      real :: t_blink_tot               ! total time for a blink switching on, off and negative
!
      do ip = 1, n_pivot
        offset = exp(-width(ip)**2/8./lam_twist(ip)**2)
        
        dist = sqrt((x(l1:l2)-xp(ip))**2 + (y(m)-yp(ip))**2)
        select case (r_profile)
        case ('gaussian')
            ux = (exp(-(dist-rad(ip))**2/(2*lam_twist(ip)**2))-offset)*udrive(ip)/(1-offset)*(-y(m)+yp(ip))/dist
            uy = (exp(-(dist-rad(ip))**2/(2*lam_twist(ip)**2))-offset)*udrive(ip)/(1-offset)*(x(l1:l2)-xp(ip))/dist
        case ('linear_exp')
            ux = exp(-abs(dist-rad(ip))/lam_twist(ip))*udrive(ip)*(-y(m)+yp(ip))
            uy = exp(-abs(dist-rad(ip))/lam_twist(ip))*udrive(ip)*(x(l1:l2)-xp(ip))
        case default
            ux = (exp(-(dist-rad(ip))**2/(2*lam_twist(ip)**2))-offset)*udrive(ip)/(1-offset)*(-y(m)+yp(ip))/dist
            uy = (exp(-(dist-rad(ip))**2/(2*lam_twist(ip)**2))-offset)*udrive(ip)/(1-offset)*(x(l1:l2)-xp(ip))/dist
        end select
        
        ! add z-dependence
        select case (z_profile)
        case ('sharp')
            z_factor = 0
            if (z(n) == xyz0(nghost)) z_factor = 1
        case ('exp')
            z_factor = exp(-(z(n)-xyz0(nghost))*lam_z)
            if (z(n) <= xyz0(nghost)) z_factor = 1
        case ('erf')
            z_factor = 1-erf((z(n)-xyz0(nghost))*lam_z)
            if (z(n) <= xyz0(nghost)) z_factor = 1
        case default
            z_factor = exp(-(z(n)-xyz0(nghost))*lam_z)
            if (z(n) <= xyz0(nghost)) z_factor = 1
        end select
        ux = z_factor*ux
        uy = z_factor*uy

        ! add 'vortex blinking' by switching the sign in time
        blink = 1
        t_blink_tot = t_blink_up(ip) + t_blink_0u(ip) + t_blink_down(ip) + t_blink_0d(ip)
        if (lblink) blink = 0
        if ((lblink) .and. (t >= delay_blink(ip))) then
            if (((t-delay_blink(ip))/t_blink_tot-floor((t-delay_blink(ip))/t_blink_tot))*t_blink_tot < t_blink_up(ip)) then
                blink = 1
            elseif ((((t-delay_blink(ip))/t_blink_tot-floor((t-delay_blink(ip))/t_blink_tot))*t_blink_tot >= t_blink_up(ip)) .and. &
                ((t-delay_blink(ip))/t_blink_tot-floor((t-delay_blink(ip))/t_blink_tot))*t_blink_tot < &
                (t_blink_up(ip)+t_blink_0u(ip))) then
                blink = 0
            elseif ((((t-delay_blink(ip))/t_blink_tot-floor((t-delay_blink(ip))/t_blink_tot))*t_blink_tot >= &
                (t_blink_up(ip)+t_blink_0u(ip))) .and. &
                ((t-delay_blink(ip))/t_blink_tot-floor((t-delay_blink(ip))/t_blink_tot))*t_blink_tot < &
                (t_blink_up(ip)+t_blink_0u(ip)+t_blink_down(ip))) then
                blink = -1
             else
                blink = 0
            endif
        endif
        
        df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) - (f(l1:l2,m,n,iux) - blink*ux)/lam_u*z_factor
        df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) - (f(l1:l2,m,n,iuy) - blink*uy)/lam_u*z_factor
      enddo
!
    endsubroutine vel_driver
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
