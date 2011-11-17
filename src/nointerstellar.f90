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
!***************************************************************
module Interstellar
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'interstellar.h'
!
  contains
!***********************************************************************
    subroutine register_interstellar()
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
    subroutine initialize_interstellar(f,lstarting)
!
!  Perform any post-parameter-read initialization eg. set derived
!  parameters
!
!  24-nov-02/tony: coded - dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_interstellar
!***********************************************************************
    subroutine input_persistent_interstellar(id,lun,done)
!
!  Read in the stored time of the next SNI
!
      integer :: id,lun
      logical :: done
!
      call keep_compiler_quiet(id)
      call keep_compiler_quiet(lun)
      call keep_compiler_quiet(done)
!
    endsubroutine input_persistent_interstellar
!***********************************************************************
    logical function output_persistent_interstellar(lun)
!
!  Writes out the time of the next SNI
!
      integer :: lun
!
      call keep_compiler_quiet(lun)

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
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
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
    subroutine read_interstellar_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
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
    subroutine read_interstellar_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
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
    subroutine pencil_criteria_interstellar()
!
!  All pencils that the Interstellar module depends on are specified here.
!
!  26-03-05/tony: coded
!
    endsubroutine pencil_criteria_interstellar
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
    subroutine check_SN(f)
!
!  dummy routine for checking for SNe (interstellar)
!
    real, dimension(mx,my,mz,mfarray) :: f
!
    call keep_compiler_quiet(f)
!
    endsubroutine check_SN
!***********************************************************************
    subroutine calc_snr_unshock(penc)
!
      real, dimension(mx), intent(inout) :: penc
!
      call keep_compiler_quiet(penc)
!
    endsubroutine calc_snr_unshock
!***********************************************************************
    subroutine calc_snr_damping(p)
!
      type (pencil_case) :: p
!
      call keep_compiler_quiet(p)
!
    endsubroutine calc_snr_damping
!***********************************************************************
    subroutine calc_snr_damp_int(int_dt)
!
      real :: int_dt
!
      call keep_compiler_quiet(int_dt)
!
    endsubroutine calc_snr_damp_int
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
