! $Id$

!  This module is used both for the initial condition and during run time.
!  It contains dlnrhon_dt and init_lnrhon, among other auxiliary routines.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Neutraldensity
 
  use Cparam
  use Cdata
  use Messages

  implicit none

  include 'neutraldensity.h'

  !integer :: dummy              ! We cannot define empty namelists
  logical :: lneutraldensity_nolog=.false.

  !namelist /neutraldensity_init_pars/ dummy
  !namelist /neutraldensity_run_pars/  dummy

  ! diagnostic variables (needs to be consistent with reset list below)
  integer :: idiag_rhonm=0

  contains

!***********************************************************************
    subroutine register_neutraldensity()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrhon; increase nvar accordingly.
!
!  18-mar-03/axel: adapted from neutraldensity
!
      use Mpicomm, only: stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_neutraldensity called twice')
      first = .false.
!
      lneutraldensity = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id$")
!
    endsubroutine register_neutraldensity
!***********************************************************************
    subroutine initialize_neutraldensity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-mar-03/axel: adapted from neutraldensity
!
!  do nothing
!
    endsubroutine initialize_neutraldensity
!***********************************************************************
    subroutine init_lnrhon(f,xx,yy,zz)
!
!  initialise lnrhon; called from start.f90
!
!  18-mar-03/axel: adapted from neutraldensity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if (NO_WARN) print*,f,xx,yy,zz ! keep compiler quiet
    endsubroutine init_lnrhon
!***********************************************************************
    subroutine pencil_criteria_neutraldensity()
!
!  All pencils that the Neutraldensity module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_neutraldensity
!***********************************************************************
    subroutine pencil_interdep_neutraldensity(lpencil_in)
!
!  Interdependency among pencils provided by the Neutraldensity module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in  !(keep compiler quiet)
!
    endsubroutine pencil_interdep_neutraldensity
!***********************************************************************
    subroutine calc_pencils_neutraldensity(f,p)
!
!  Calculate Neutraldensity pencils.
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
    endsubroutine calc_pencils_neutraldensity
!***********************************************************************
    subroutine dlnrhon_dt(f,df,p)
!
!  continuity equation
!  calculate dlnrhon/dt = - u.gradlnrhon - divud
!
!  18-mar-03/axel: adapted from neutraldensity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      if (NO_WARN) print*,f,df,p !(keep compiler quiet)
!
    endsubroutine dlnrhon_dt
!***********************************************************************
    subroutine read_neutraldensity_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_neutraldensity_init_pars
!***********************************************************************
    subroutine write_neutraldensity_init_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit

    endsubroutine write_neutraldensity_init_pars
!***********************************************************************
    subroutine read_neutraldensity_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_neutraldensity_run_pars
!***********************************************************************
    subroutine write_neutraldensity_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit
    endsubroutine write_neutraldensity_run_pars
!***********************************************************************
    subroutine rprint_neutraldensity(lreset,lwrite)
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
!  write column where which neutral density variable is stored
!
      if (lwr) then
        write(3,*) 'i_rhonm=',idiag_rhonm
        write(3,*) 'nname=',nname
        write(3,*) 'ilnrhon=',ilnrhon
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_neutraldensity
!***********************************************************************

endmodule Neutraldensity
