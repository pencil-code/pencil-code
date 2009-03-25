! $Id$

!  This modules implements viscous heating and diffusion terms
!  here for cases 1) nu constant, 2) mu = rho.nu 3) constant and

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lviscosity = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED gshock(3); shock; visc_heat
!
!***************************************************************
module Viscosity
!
  use Cdata
  use Messages
!
  implicit none
!
  include 'viscosity.h'
!
  logical :: lvisc_first=.false.
!
  contains
!***********************************************************************
    subroutine register_viscosity()
!
!  19-nov-02/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
!  Identify version number.
!
      if (lroot) call cvs_id( &
           "$Id$")
!
    endsubroutine register_viscosity
!***********************************************************************
    subroutine initialize_viscosity(lstarting)
!
!  02-apr-02/tony: coded
      logical, intent(in) :: lstarting

      if (NO_WARN) print*,lstarting
    endsubroutine initialize_viscosity
!***********************************************************************
    subroutine read_viscosity_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_viscosity_init_pars
!***********************************************************************
    subroutine write_viscosity_init_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit

    endsubroutine write_viscosity_init_pars
!***********************************************************************
    subroutine read_viscosity_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_viscosity_run_pars
!***********************************************************************
    subroutine write_viscosity_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit

    endsubroutine write_viscosity_run_pars
!*******************************************************************
    subroutine rprint_viscosity(lreset,lwrite)
!
!  Writes ishock to index.pro file
!
!  02-apr-03/tony: adapted from visc_const
!
      logical :: lreset
      logical, optional :: lwrite
!
      if (NO_WARN) print*,lreset,present(lwrite)  !(to keep compiler quiet)
    endsubroutine rprint_viscosity
!***********************************************************************
    subroutine pencil_criteria_viscosity()
!
!  All pencils that the Viscosity module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_viscosity
!***********************************************************************
    subroutine pencil_interdep_viscosity(lpencil_in)
!
!  Interdependency among pencils from the Viscosity module is specified here.
!
!  20-11-04/anders: coded
!
      use Cdata
!
      logical, dimension (npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in !(keep compiler quiet)
!
    endsubroutine pencil_interdep_viscosity
!***********************************************************************
    subroutine calc_pencils_viscosity(f,p)
!
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
!  define shock, gshock, visc_heat
!
      if (lpencil(i_shock)) p%shock=0.
      if (lpencil(i_gshock)) p%gshock=0.
      if (lpencil(i_visc_heat)) p%visc_heat=0.
!
      if (NO_WARN) print*, f
!
    endsubroutine calc_pencils_viscosity
!!***********************************************************************
    subroutine calc_viscosity(f)
!
!
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*,f  !(to keep compiler quiet)
!
    endsubroutine calc_viscosity
!!***********************************************************************
    subroutine calc_viscous_heat(f,df,p,Hmax)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  20-nov-02/tony: coded
!
      use Cdata

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar)    :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax
!
      intent(in) :: df,p,Hmax
!
      if (NO_WARN) print*,df,p,Hmax  !(keep compiler quiet)
!
    endsubroutine calc_viscous_heat
!***********************************************************************
    subroutine calc_viscous_force(df,p)
!
!  calculate viscous heating term for right hand side of entropy equation
!
!  20-nov-02/tony: coded
!
      use Cdata

      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent (in) :: df,p
!
      if (NO_WARN) print*,df,p  !(keep compiler quiet)
!
    endsubroutine calc_viscous_force
!***********************************************************************
    subroutine calc_visc_heat_ppd(f,df,p)
!    
!  Dummy routine
! 
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      intent (in) :: f,df,p
!
      if (NO_WARN) print*,f,df,p  !(keep compiler quiet)
!
    endsubroutine calc_visc_heat_ppd
!***********************************************************************
    subroutine getnu(nu_)
!
!  Dummy routine
!
      real,intent(out) :: nu_
!
      nu_=0.0
!
      if (NO_WARN) print*,nu_  !(keep compiler quiet)
!
    endsubroutine getnu
!***********************************************************************

endmodule Viscosity
