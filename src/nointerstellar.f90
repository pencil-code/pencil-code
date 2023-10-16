! $Id$
!
!  Dummy module.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linterstellar = .false.
!
! PENCILS PROVIDED heat; cool; heatcool
!
!***************************************************************
module Interstellar
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'interstellar.h'
!
  contains
!***********************************************************************
    subroutine register_interstellar
!
!  19-nov-02/tony: coded
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_interstellar
!***********************************************************************
    subroutine initialize_interstellar(f)
!
!  Perform any post-parameter-read initialization eg. set derived
!  parameters
!
!  24-nov-02/tony: coded - dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_interstellar
!***********************************************************************
    subroutine input_persistent_interstellar(id,done)
!
!  Read in the stored time of the next SNI
!
      integer, optional :: id
      logical, optional :: done
!
      if (present (id)) call keep_compiler_quiet(id)
      if (present (done)) call keep_compiler_quiet(done)
!
    endsubroutine input_persistent_interstellar
!***********************************************************************
    logical function output_persistent_interstellar()
!
!  Writes out the time of the next SNI
!
      output_persistent_interstellar = .false.
!
    endfunction output_persistent_interstellar
!***********************************************************************
    subroutine rprint_interstellar(lreset,lwrite)
!
!  reads and registers print parameters relevant to interstellar
!
!   1-jun-02/axel: adapted from magnetic fields
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_interstellar
!***********************************************************************
    subroutine get_slices_interstellar(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_interstellar
!***********************************************************************
    subroutine read_interstellar_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_interstellar_init_pars
!***********************************************************************
    subroutine write_interstellar_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_interstellar_init_pars
!***********************************************************************
    subroutine read_interstellar_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_interstellar_run_pars
!***********************************************************************
    subroutine write_interstellar_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_interstellar_run_pars
!***********************************************************************
    subroutine init_interstellar(f)
!
!  initialise magnetic field; called from start.f90
!  30-jul-2006/tony: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_interstellar
!***********************************************************************
    subroutine pencil_criteria_interstellar
!
!  All pencils that the Interstellar module depends on are specified here.
!
!  26-03-05/tony: coded
!
    endsubroutine pencil_criteria_interstellar
!***********************************************************************
    subroutine pencil_interdep_interstellar(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_cool)) then
        lpencil_in(i_lnTT)=.true.
        lpencil_in(i_lnrho)=.true.
      endif
      if (lpencil_in(i_heat)) then
        lpencil_in(i_lnTT)=.true.
      endif
      if (lpencil_in(i_heatcool)) then
        lpencil_in(i_cool)=.true.
        lpencil_in(i_heat)=.true.
      endif
!
    endsubroutine pencil_interdep_interstellar
!***********************************************************************
    subroutine calc_pencils_interstellar(f,p)
!
!  Calculate Interstellar pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!
      real, dimension(mx,my,mz,mfarray), intent(IN)   :: f
      type(pencil_case),                 intent(INOUT):: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!      if (lpencil(i_cool)) call calc_cool_func(p%cool,p%lnTT,p%lnrho)
!!
!      if (lpencil(i_heat)) call calc_heat(p%heat,p%lnTT)
!!
!      if (lpencil(i_heatcool)) p%heatcool=p%heat-p%cool
!
    endsubroutine calc_pencils_interstellar
!***********************************************************************
    subroutine calc_diagnostics_interstellar(p)

      type(pencil_case) :: p

      call keep_compiler_quiet(p)

    endsubroutine calc_diagnostics_interstellar
!***********************************************************************
    subroutine interstellar_after_boundary(f)
!
!  This routine calculates and applies the optically thin cooling function
!  together with UV heating.
!
!  01-aug-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine interstellar_after_boundary
!***********************************************************************
    subroutine interstellar_before_boundary(f)
!
!  This routine calculates and applies the optically thin cooling function
!  together with UV heating.
!
!  01-aug-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine interstellar_before_boundary
!***********************************************************************
    subroutine check_SN(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine check_SN
!***********************************************************************
    subroutine calc_heat_cool_interstellar(f,df,p,Hmax)
!
!  adapted from calc_heat_cool
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(in) :: df
      type (pencil_case), intent(in) :: p
      real, dimension(nx), intent(in) :: Hmax
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(Hmax)
!
    endsubroutine calc_heat_cool_interstellar
!***********************************************************************
    subroutine addmassflux(f)
!
      real, intent(inout), dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine addmassflux
!***********************************************************************
endmodule interstellar
