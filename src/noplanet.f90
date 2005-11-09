! $Id: noplanet.f90,v 1.3 2005-11-09 00:33:54 wlyra Exp $
!
!  Dummy module
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lplanet = .false.
!
!***************************************************************

module Planet

  use Cparam
  use Cdata
  use Messages

  implicit none

  include 'planet.h'

  real :: gc=0.

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
           "$Id: noplanet.f90,v 1.3 2005-11-09 00:33:54 wlyra Exp $")
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
    subroutine gravity_companion(f,df,p)
      
!8-nov-05/wlad : dummy      

      use Cdata
      use Mpicomm, only: stop_it

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      call stop_it("gravity_companion: noplanet is used")
      
    endsubroutine gravity_companion
!***********************************************************************
    subroutine calc_torque(unit)
      integer, intent(in) :: unit
      if (NO_WARN) print*,unit
    endsubroutine calc_torque
!**********************************************************************
!    subroutine local_isothermal(cs20,corr)
!      use Mpicomm, only: stop_it
!      use Cdata
!
!      real, intent(in)  :: cs20
!      real, dimension (nx) :: corr
!
!      call stop_it("local_isothermal: llocal_iso true but noplanet is used")
!
!    endsubroutine local_isothermal
!**********************************************************************
  endmodule Planet
