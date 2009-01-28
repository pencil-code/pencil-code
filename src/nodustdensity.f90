! $Id$

!  This module is used both for the initial condition and during run time.
!  It contains dlnrhod_dt and init_lnrhod, among other auxiliary routines.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Dustdensity

  use Cparam
  use Cdata
  use Messages

  implicit none

  include 'dustdensity.h'

  !integer :: dummy              ! We cannot define empty namelists
  logical :: ldustnulling=.false.

  !namelist /dustdensity_init_pars/ dummy
  !namelist /dustdensity_run_pars/  dummy

  ! diagnostic variables (needs to be consistent with reset list below)
  integer :: idiag_rhodm=0

  contains

!***********************************************************************
    subroutine register_dustdensity()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrhod; increase nvar accordingly.
!
!  18-mar-03/axel: adapted from dustdensity
!
      use Mpicomm, only: stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_dustdensity called twice')
      first = .false.
!
      ldustdensity = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id$")
!
    endsubroutine register_dustdensity
!***********************************************************************
    subroutine initialize_dustdensity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-mar-03/axel: adapted from dustdensity
!
!  do nothing
!
    endsubroutine initialize_dustdensity
!***********************************************************************
    subroutine init_nd(f)
!
!  initialise lnrhod; called from start.f90
!
!  18-mar-03/axel: adapted from dustdensity
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*,f ! keep compiler quiet
    endsubroutine init_nd
!***********************************************************************
    subroutine pencil_criteria_dustdensity()
!
!  All pencils that the Dustdensity module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_dustdensity
!***********************************************************************
    subroutine pencil_interdep_dustdensity(lpencil_in)
!
!  Interdependency among pencils provided by the Dustdensity module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in  !(keep compiler quiet)
!
    endsubroutine pencil_interdep_dustdensity
!***********************************************************************
    subroutine calc_pencils_dustdensity(f,p)
!
!  Calculate Dustdensity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  13-nov-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f, p
!
      if (NO_WARN) print*, f, p
!
    endsubroutine calc_pencils_dustdensity
!***********************************************************************
    subroutine dndmd_dt(f,df,p)
!
!  continuity equation
!  calculate dlnrhod/dt = - u.gradlnrhod - divud
!
!  18-mar-03/axel: adapted from dustdensity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      if (NO_WARN) print*,f,df,p !(keep compiler quiet)
!
    endsubroutine dndmd_dt
!***********************************************************************
    subroutine read_dustdensity_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_dustdensity_init_pars
!***********************************************************************
    subroutine write_dustdensity_init_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit

    endsubroutine write_dustdensity_init_pars
!***********************************************************************
    subroutine read_dustdensity_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_dustdensity_run_pars
!***********************************************************************
    subroutine write_dustdensity_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit
    endsubroutine write_dustdensity_run_pars
!***********************************************************************
    subroutine redist_mdbins(f)
!
!  Redistribute dust number density and dust density in mass bins
!
!  4-may-2004/wolf: Adapted from dustdensity.f90
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*,f
    endsubroutine redist_mdbins
!***********************************************************************
    subroutine null_dust_vars(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*,f
    endsubroutine null_dust_vars
!***********************************************************************
    subroutine reinit_criteria_dust
!
!  Force reiniting of dust variables if certain criteria are fulfilled
!
!  Dummy subroutine
!
    endsubroutine reinit_criteria_dust
!***********************************************************************
    subroutine rprint_dustdensity(lreset,lwrite)
!
!  reads and registers print parameters relevant for compressible part
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Sub
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  write column where which dust density variable is stored
!
      if (lwr) then
        write(3,*) 'i_rhodm=',idiag_rhodm
        write(3,*) 'nname=',nname
        write(3,*) 'ind=',ind
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_dustdensity
!***********************************************************************

endmodule Dustdensity
