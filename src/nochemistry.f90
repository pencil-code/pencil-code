! $Id$
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lchemistry = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Chemistry
!
  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'chemistry.h'
!
  contains
!***********************************************************************
    subroutine register_chemistry()
!
    endsubroutine register_chemistry
!***********************************************************************
    subroutine initialize_chemistry(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_chemistry
!***********************************************************************
    subroutine init_chemistry(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_chemistry
!***********************************************************************
    subroutine pencil_criteria_chemistry()
!
    endsubroutine pencil_criteria_chemistry
!***********************************************************************
    subroutine pencil_interdep_chemistry(lpencil_in)
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_chemistry
!***********************************************************************
    subroutine calc_pencils_chemistry(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_chemistry
!***********************************************************************
    subroutine calc_for_chem_mixture(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_for_chem_mixture
!**********************************************************************
    subroutine dchemistry_dt(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dchemistry_dt
!***********************************************************************
    subroutine read_chemistry_init_pars(unit,iostat)
!
      integer :: unit
      integer, optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_chemistry_init_pars
!***********************************************************************
   subroutine write_chemistry_init_pars(unit)
!
      integer :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_chemistry_init_pars
!***********************************************************************
    subroutine read_chemistry_run_pars(unit,iostat)
!
      integer :: unit
      integer, optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_chemistry_run_pars
!***********************************************************************
    subroutine write_chemistry_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_chemistry_run_pars
!***********************************************************************
    subroutine rprint_chemistry(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_chemistry
!***********************************************************************
    subroutine get_slices_chemistry(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices)
!
    endsubroutine get_slices_chemistry
!********************************************************************
    subroutine bc_nscbc_subin_x(f,df,topbot,val)
!
!   nscbc case 
!   subsonic inflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    ux,uy,T are fixed, drho  is calculated 
!
!   16-nov-08/natalia: coded.
!
      use MpiComm, only: stop_it
      use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice, der_pencil
!      use Chemistry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot
      real, dimension (mcom), optional :: val
!
      intent(inout) :: f
      intent(out) :: df
!

     call keep_compiler_quiet(f)
     call keep_compiler_quiet(df)

    endsubroutine bc_nscbc_subin_x
!********************************************************************
  subroutine bc_nscbc_nref_subout_x(f,df,topbot)
!
!   nscbc case 
!   subsonic non-reflecting outflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    drho. dT, dux, duy  are calculated, p_inf can be 
!   fixed (if nscbc_sigma <>0)
!
!   16-nov-08/natalia: coded.
!
      use MpiComm, only: stop_it
      !use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot

     
     call keep_compiler_quiet(f)
     call keep_compiler_quiet(df)
!
    endsubroutine bc_nscbc_nref_subout_x
!***********************************************************************
 subroutine bc_nscbc_nref_subout_y(f,df,topbot)
!
!   nscbc case 
!   subsonic non-reflecting outflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    drho. dT, dux, duy  are calculated, p_inf can be 
!   fixed (if nscbc_sigma <>0)
!
!   16-nov-08/natalia: coded.
!
      use MpiComm, only: stop_it
      !use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot

     
     call keep_compiler_quiet(f)
     call keep_compiler_quiet(df)
!
    endsubroutine bc_nscbc_nref_subout_y
!***********************************************************************
 subroutine bc_nscbc_nref_subout_z(f,df,topbot)
!
!   nscbc case 
!   subsonic non-reflecting outflow boundary conditions
!   now it is 2D case (should be corrected for 3D case)
!    drho. dT, dux, duy  are calculated, p_inf can be 
!   fixed (if nscbc_sigma <>0)
!
!   16-nov-08/natalia: coded.
!
      use MpiComm, only: stop_it
      !use EquationOfState, only: cs0, cs20
      use Deriv, only: der_onesided_4_slice,der_pencil
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      character (len=3) :: topbot

     
     call keep_compiler_quiet(f)
     call keep_compiler_quiet(df)
!
    endsubroutine bc_nscbc_nref_subout_z
!***********************************************************************
  subroutine chemistry_clean_up()
!
  endsubroutine chemistry_clean_up
!********************************************************************
endmodule Chemistry
