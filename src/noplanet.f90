! $Id: noplanet.f90,v 1.18 2006-04-05 14:35:18 wlyra Exp $
!
!  Dummy module
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lplanet = .false.
!
! PENCILS PROVIDED uu_kep
!
!***************************************************************

module Planet

  use Cparam
  use Cdata
  use Messages

  implicit none

  include 'planet.h'

  real :: gc=0.          !location and mass
  real :: b=0.           !peak radius for potential
  integer :: nc=2        !exponent of smoothed potential 
  integer :: n_periods=5. !periods for ramping
  logical :: lramp=.false.
  logical :: lwavedamp=.false.,llocal_iso=.false.

  !namelist /planet_init_pars/ dummy
  !namelist /planet_run_pars/ dummy
 
  contains

!***********************************************************************
    subroutine register_planet()
!
!  07-nov-05/wlad: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_noplanet called twice')
      first = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: noplanet.f90,v 1.18 2006-04-05 14:35:18 wlyra Exp $")
!
!      if (nvar > mvar) then
!        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
!        call stop_it('Register_nointerstellar: nvar > mvar')
!      endif
!
    endsubroutine register_planet
!***********************************************************************
    subroutine initialize_planet(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  05-nov-05/wlad: dummy
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      logical :: lstarting
!
      if (NO_WARN) print*, lstarting
!
    endsubroutine initialize_planet
!***********************************************************************
    subroutine rprint_planet(lreset,lwrite)
!
!  reads and registers print parameters relevant to entropy
!
!   7-nov-05/wlad: adapted from interstellar
!
      use Cdata
      use Sub
!
!      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  write column where which totalmass is stored
!
      if (lwr) then
        write(3,*) 'i_totalmass=0'
      endif
!
      if (NO_WARN) print*,lreset
!
    endsubroutine rprint_planet
!!***************************************************    
    subroutine read_planet_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit
    endsubroutine read_planet_init_pars
!***********************************************************************
    subroutine write_planet_init_pars(unit)
      integer, intent(in) :: unit
      if (NO_WARN) print*,unit
    endsubroutine write_planet_init_pars
!***********************************************************************
    subroutine read_planet_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit
    endsubroutine read_planet_run_pars
!***********************************************************************
    subroutine write_planet_run_pars(unit)
      integer, intent(in) :: unit
      if (NO_WARN) print*,unit
    endsubroutine write_planet_run_pars
!***********************************************************************
    subroutine pencil_criteria_planet()
      !
      !
    endsubroutine pencil_criteria_planet
!***********************************************************************
    subroutine calc_pencils_planet(f,p,g0,r0_pot)
!
! calculate keplerian velocity as a pencil, so I don't
! have to recalculate it on other routines from f,g0 and r0_pot
!
! 28-mar-06/wlad : dummy
!
     use Cdata
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
      real :: g0,r0_pot
!
      intent(in) :: f
      intent(inout) :: p
!
      if (lpencil(i_uu_kep)) p%uu_kep=0.
!   
    endsubroutine calc_pencils_planet
!***********************************************************************
    subroutine gravity_companion(f,df,fp,dfp,gs,r0_pot,n_pot,p)
!      
!8-nov-05/wlad : dummy      
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (mpar_loc,mpvar) :: fp,dfp
      type (pencil_case) :: p
      real :: gs,r0_pot
      integer :: n_pot
!
      call stop_it("noplanet.f90 - gravity_companion")
!      
    endsubroutine gravity_companion
!***********************************************************************
    subroutine calc_torque(unit)
      integer, intent(in) :: unit
      if (NO_WARN) print*,unit
    endsubroutine calc_torque
!**********************************************************************
    subroutine local_isothermal(cs20,corr)
!
!8-nov-05/wlad : dummy  
!
      use Mpicomm, only: stop_it
      use Cdata
!
      real, intent(in)  :: cs20
      real, dimension (nx) :: corr
!
      call stop_it("noplanet.f90 - local_isothermal")
!
    endsubroutine local_isothermal
!**********************************************************************
   subroutine gravity_star(gs,r0_pot,n_pot,g_r,xstar,ystar)
!
!8-nov-05/wlad : dummy  
!
     use Cdata
     use Mpicomm, only: stop_it
!
     real, dimension (nx), intent(out) :: g_r
     real, optional :: xstar,ystar !initial position of star
     real :: gs,r0_pot
     integer :: n_pot
!
     g_r=0.
     call stop_it("noplanet.f90 - gravity_star")
!
   endsubroutine gravity_star
!***************************************************************
    subroutine wave_damping(f,df)
!
!20-nov-05/wlad : dummy  
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(mx,my,mz,mvar) :: df
!
      call stop_it("noplanet.f90 - wave_damping")
!    
    endsubroutine wave_damping
!***************************************************************
    subroutine time_average(p)
!
! 02-03-06/wlad : dummy
!
      type (pencil_case) :: p
!
    endsubroutine time_average
!***************************************************************
    subroutine get_ramped_mass(gp,gs,g0,mdot,m2dot)
!
! 03-mar-06/wlad :: dummy
!
      use Mpicomm, only: stop_it
!     
      real :: gs,gp,g0,mdot,m2dot
!
      intent(out) :: gs,gp,mdot,m2dot
!
      gs=1. ; gp=0. ; mdot=0. ; m2dot = 0.
!
      call stop_it("noplanet.f90 - get_ramped_mass")
!
    endsubroutine get_ramped_mass
!***************************************************************
  endmodule Planet
