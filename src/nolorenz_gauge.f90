! $Id$
!
!  Lorenz gauge, dphi/dt = -cphi2*divA
!
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: llorenz_gauge = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Lorenz_gauge

  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages

  implicit none

  include 'lorenz_gauge.h'
!
  contains
!***********************************************************************
    subroutine register_lorenz_gauge
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_lorenz_gauge
!***********************************************************************
    subroutine initialize_lorenz_gauge(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)
!
    endsubroutine initialize_lorenz_gauge
!***********************************************************************
    subroutine init_lorenz_gauge(f)
!
!  initialise lorenz_gauge condition; called from start.f90
!  06-oct-2003/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_lorenz_gauge
!***********************************************************************
    subroutine pencil_criteria_lorenz_gauge
!
!  All pencils that this lorenz_gauge module depends on are specified here.
!
!  25-feb-07/axel: adapted

    endsubroutine pencil_criteria_lorenz_gauge
!***********************************************************************
    subroutine pencil_interdep_lorenz_gauge(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_lorenz_gauge
!***********************************************************************
    subroutine calc_pencils_lorenz_gauge(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_lorenz_gauge
!***********************************************************************
    subroutine dlorenz_gauge_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dlorenz_gauge_dt
!***********************************************************************
    subroutine read_lorenz_gauge_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_lorenz_gauge_init_pars
!***********************************************************************
    subroutine write_lorenz_gauge_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_lorenz_gauge_init_pars
!***********************************************************************
    subroutine read_lorenz_gauge_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_lorenz_gauge_run_pars
!***********************************************************************
    subroutine write_lorenz_gauge_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_lorenz_gauge_run_pars
!***********************************************************************
    subroutine rprint_lorenz_gauge(lreset,lwrite)
!
!  reads and registers print parameters relevant to lorenz_gauge
!
!   06-oct-03/tony: coded
!
!   define counters
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_lorenz_gauge
!***********************************************************************
    subroutine get_slices_lorenz_gauge(f,slices)
!
!  Write slices for animation of electric potential
!
!  26-feb-07/axel: adapted from gross_pitaevskii
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_lorenz_gauge
!***********************************************************************
endmodule Lorenz_gauge

