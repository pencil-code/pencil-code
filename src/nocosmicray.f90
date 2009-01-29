! $Id$

!  This modules solves the passive scalar advection equation

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lcosmicray = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Cosmicray

  use Cparam
  use Cdata
  use Messages

  implicit none

  include 'cosmicray.h'

  !namelist /cosmicray_init_pars/ dummy
  !namelist /cosmicray_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
!  integer :: idiag_rhocrm=0,idiag_crmax=0,idiag_lncrm=0,idiag_lncrmz=0

  contains

!***********************************************************************
    subroutine register_cosmicray()
!
!  Initialise variables which should know that we solve for active
!  scalar: iecr - the cosmic ray energy density; increase nvar accordingly
!
!  09-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_cosmicray called twice')
      first = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id$")
!
    endsubroutine register_cosmicray
!***********************************************************************
    subroutine initialize_cosmicray(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  09-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if (NO_WARN) print*,'f=',f
!
    endsubroutine initialize_cosmicray
!***************************************:********************************
    subroutine read_cosmicray_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_cosmicray_init_pars
!***********************************************************************
    subroutine write_cosmicray_init_pars(unit)
!
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit
!
    endsubroutine write_cosmicray_init_pars
!***********************************************************************
    subroutine read_cosmicray_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat) .and. (NO_WARN)) print*,iostat
!
      if (NO_WARN) print*,unit
!
    endsubroutine read_cosmicray_run_pars
!***********************************************************************
    subroutine write_cosmicray_run_pars(unit)
!
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit
!
    endsubroutine write_cosmicray_run_pars
!***********************************************************************
    subroutine init_ecr(f)
!
!  initialise magnetic field; called from start.f90
!  AB: maybe we should here call different routines (such as rings)
!  AB: and others, instead of accummulating all this in a huge routine.
!  We have an init parameter (initlncc) to stear magnetic i.c. independently.
!
!   6-jul-02/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*,f !(prevent compiler warnings)
!
    endsubroutine init_ecr
!***********************************************************************
    subroutine pencil_criteria_cosmicray()
!
!  All pencils that the Cosmicray module depends on are specified here.
!
!  20-11-04/anders: coded
!
!
    endsubroutine pencil_criteria_cosmicray
!***********************************************************************
    subroutine pencil_interdep_cosmicray(lpencil_in)
!
!  Interdependency among pencils provided by the Cosmicray module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension (npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in !(keep compiler quiet)
!
    endsubroutine pencil_interdep_cosmicray
!***********************************************************************
    subroutine calc_pencils_cosmicray(f,p)
!
!  Calculate Cosmicray pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      if (NO_WARN) print*, f, p !(keep compiler quiet)
!
    endsubroutine calc_pencils_cosmicray
!***********************************************************************
    subroutine decr_dt(f,df,p)
!
!  cosmic ray density evolution
!
!   09-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,df,p
!
      if (NO_WARN) print*,f,df,p
!
    endsubroutine decr_dt
!***********************************************************************
    subroutine rprint_cosmicray(lreset,lwrite)
!
!  reads and registers print parameters relevant for cosmic rays
!
!   09-oct-03/tony: coded
!
      use Sub
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
!        idiag_rhoccm=0; idiag_ccmax=0; idiag_lnccm=0; idiag_lnccmz=0
      endif
!
!  write column where which cosmic ray variable is stored
!
      if (lwr) then
!        write(3,*) 'i_lnccmz=',idiag_lnccmz
        write(3,*) 'iecr=',iecr
      endif
!
    endsubroutine rprint_cosmicray
!***********************************************************************

endmodule Cosmicray
