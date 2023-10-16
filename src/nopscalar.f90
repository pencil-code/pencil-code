! $Id$
!
!  This modules solves the passive scalar advection equation.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpscalar = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
! PENCILS PROVIDED cc; cc1; gcc(3,0)
!
!***************************************************************
module Pscalar
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'pscalar.h'
!
  real :: rhoccm=0.0, cc2m=0.0, gcc2m=0.0
  integer :: idiag_gcc2m=0, idiag_cc2m=0, idiag_rhoccm=0
!
  contains
!***********************************************************************
    subroutine register_pscalar
!
!  Initialise variables which should know that we solve for passive
!  scalar: ilncc; increase nvar accordingly.
!
!  6-jul-02/axel: coded
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_pscalar
!***********************************************************************
    subroutine initialize_pscalar(f)
!
!  Perform any necessary post-parameter read initialization.
!
!  24-nov-02/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_pscalar
!***********************************************************************
    subroutine init_lncc(f)
!
!  Initialise passive scalar field.
!
!   6-jul-02/axel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_lncc
!***********************************************************************
    subroutine pencil_criteria_pscalar
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
      call keep_compiler_quiet(lpencil_in)
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
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
! cc
      if (lpencil(i_cc)) p%cc=1.0
! cc1
      if (lpencil(i_cc1)) p%cc1=1.0
! gcc
      if (lpencil(i_gcc)) p%gcc=0.0
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_pscalar
!***********************************************************************
    subroutine calc_diagnostics_pscalar(p)

      type (pencil_case) :: p

      call keep_compiler_quiet(p)
 
    endsubroutine calc_diagnostics_pscalar
!***********************************************************************
    subroutine dlncc_dt(f,df,p)
!
!  Passive scalar evolution.
!
!   6-jul-02/axel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)  :: f,df,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine dlncc_dt
!***********************************************************************
    subroutine read_pscalar_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_pscalar_init_pars
!***********************************************************************
    subroutine write_pscalar_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_pscalar_init_pars
!***********************************************************************
    subroutine read_pscalar_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_pscalar_run_pars
!***********************************************************************
    subroutine write_pscalar_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_pscalar_run_pars
!***********************************************************************
    subroutine rprint_pscalar(lreset,lwrite)
!
!  Reads and registers print parameters relevant for passive scalar.
!
!   6-jul-02/axel: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_pscalar
!***********************************************************************
    subroutine get_slices_pscalar(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_pscalar
!***********************************************************************
    subroutine pscalar_before_boundary(f)

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f

      call keep_compiler_quiet(f)

    endsubroutine pscalar_before_boundary
!***********************************************************************
    subroutine calc_mpscalar
!
    endsubroutine calc_mpscalar
!***********************************************************************
endmodule Pscalar
