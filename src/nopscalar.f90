! $Id$

!  This modules solves the passive scalar advection equation

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
! PENCILS PROVIDED cc; cc1
!
!***************************************************************

module Pscalar

  use Cparam
  use Cdata
  use Messages

  implicit none

  include 'pscalar.h'

!ajwm needed by prints.f90 for some reason.
  real :: rhoccm=0., cc2m=0., gcc2m=0.
  integer :: idiag_gcc2m=0, idiag_cc2m=0, idiag_rhoccm=0
  !namelist /pscalar_init_pars/ dummy
  !namelist /pscalar_run_pars/  dummy

  contains

!***********************************************************************
    subroutine register_pscalar()
!
!  Initialise variables which should know that we solve for passive
!  scalar: ilncc; increase nvar accordingly
!
!  6-jul-02/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_lncc called twice')
      first = .false.
!
      lpscalar = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id$")
!
    endsubroutine register_pscalar
!***********************************************************************
    subroutine initialize_pscalar(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  24-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if (NO_WARN) print*,'f=',f
    endsubroutine initialize_pscalar
!***********************************************************************
    subroutine init_lncc(f)
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
    endsubroutine init_lncc
!***********************************************************************
    subroutine pencil_criteria_pscalar()
!
!  All pencils that the Pscalar module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_pscalar
!***********************************************************************
    subroutine pencil_interdep_pscalar(lpencil_in)
!
!  Interdependency among pencils provided by the Pscalar module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in !(keep compiler quiet)
!
    endsubroutine pencil_interdep_pscalar
!***********************************************************************
    subroutine calc_pencils_pscalar(f,p)
!
!  Calculate Pscalar pencils.
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
! cc
      if (lpencil(i_cc)) p%cc=1.
! cc1
      if (lpencil(i_cc1)) p%cc1=1.
!
      if (NO_WARN) print*, f
!
    endsubroutine calc_pencils_pscalar
!***********************************************************************
    subroutine dlncc_dt(f,df,p)
!
!  magnetic field evolution
!
!  calculate dc/dt=-uu.gradlncc
!
!   6-jul-02/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)  :: f,df,p
!
      if (NO_WARN) print*,f,df,p
!
    endsubroutine dlncc_dt
!***********************************************************************
    subroutine read_pscalar_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_pscalar_init_pars
!***********************************************************************
    subroutine write_pscalar_init_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit

    endsubroutine write_pscalar_init_pars
!***********************************************************************
    subroutine read_pscalar_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_pscalar_run_pars
!***********************************************************************
    subroutine write_pscalar_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit
    endsubroutine write_pscalar_run_pars
!***********************************************************************
    subroutine rprint_pscalar(lreset,lwrite)
!
!  reads and registers print parameters relevant for passive scalar
!
!   6-jul-02/axel: coded
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
      endif
!
!  write column where which passive scalar variable is stored
!
      if (lwr) then
        write(3,*) 'i_rhoccm=0'
        write(3,*) 'i_ccmax=0'
        write(3,*) 'i_lnccm=0'
        write(3,*) 'i_lnccmz=0'
        write(3,*) 'ilncc=0'
        write(3,*) 'icc=0'
      endif
!
    endsubroutine rprint_pscalar
!***********************************************************************
    subroutine calc_mpscalar
    endsubroutine calc_mpscalar
!***********************************************************************

endmodule Pscalar
