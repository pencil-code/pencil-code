! $Id$
!
!  This module contains routines both for delta-correlated
!  and continuous forcing. The fcont pencil is only provided
!  for continuous forcing.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lforcing = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED fcont(3,n_forcing_cont_max)
!
!***************************************************************
module Forcing
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  logical :: lhydro_forcing=.false.,ltestflow_forcing=.false.

  include 'forcing.h'
!
  integer :: n_forcing_cont=n_forcing_cont_max
  integer :: pushpars2c, pushdiags2c  ! should be procedure pointer (F2003)
!
  contains
!***********************************************************************
    subroutine register_forcing
!
!  add forcing in timestep()
!  11-may-2002/wolf: coded
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine initialize_forcing
!
!  initialize random number generator in processor-dependent fashion
!  see comments in start.f90 for details
!
    endsubroutine initialize_forcing
!***********************************************************************
    subroutine addforce(f)
!
!  add forcing in timestep()
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine addforce
!***********************************************************************
    subroutine forcing_cont_after_boundary(f)
!
!  precalculate parameters that are new at each timestep,
!  but the same for all pencils
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine forcing_cont_after_boundary
!***********************************************************************
    subroutine pencil_criteria_forcing
!
!  Dummy routine
!
    endsubroutine pencil_criteria_forcing
!***********************************************************************
    subroutine pencil_interdep_forcing(lpencil_in)
!
!  Dummy routine
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_forcing
!***********************************************************************
    subroutine calc_pencils_forcing(f,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (lpencil(i_fcont)) p%fcont=0.
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_forcing
!***********************************************************************
    subroutine forcing_continuous(df,p)
!
!  dummy routine
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine forcing_continuous
!***********************************************************************
    subroutine forcing_cont(force)
!
      real, dimension (nx,3), intent(out) :: force
!
      call keep_compiler_quiet(force)
!
    endsubroutine forcing_cont
!***********************************************************************
    subroutine read_forcing_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_forcing_run_pars
!***********************************************************************
    subroutine write_forcing_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_forcing_run_pars
!***********************************************************************
    subroutine input_persistent_forcing(id,done)
!
!  Read in the stored time of the next SNI
!
      integer :: id
      logical :: done
!
      call keep_compiler_quiet(id)
      call keep_compiler_quiet(done)
!
    endsubroutine input_persistent_forcing
!***********************************************************************
    logical function output_persistent_forcing()
!
!  Writes out the time of the next SNI
!
!  16-nov-11/MR: changed into logical function
!
      output_persistent_forcing = .false.

    endfunction output_persistent_forcing
!***********************************************************************
    subroutine rprint_forcing(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!  26-jan-04/axel: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)

    endsubroutine rprint_forcing
!***********************************************************************
    subroutine fconst_coefs_hel(force_fact,kkx,kky,kkz,nk,kav,coef1,coef2,coef3,kk,phase,fact,fda)
!
      use General, only: keep_compiler_quiet
!
      real,                   intent(in ) :: force_fact,kav
      integer,                intent(in ) :: nk
      real,    dimension (nk),intent(in ) :: kkx,kky,kkz
      real,    dimension (3), intent(out) :: coef1,coef2,coef3,kk,fda
      real,                   intent(out) :: phase,fact

      call keep_compiler_quiet(force_fact,kav,phase,fact)
      call keep_compiler_quiet(kkx,kky,kkz,fda)
      call keep_compiler_quiet(coef1,coef2,coef3,kk)
      call keep_compiler_quiet(nk)

    endsubroutine fconst_coefs_hel
!***********************************************************************
    subroutine forcing_clean_up
!
!   12-aug-09/dhruba: coded
!   dummy routine.
!
    endsubroutine forcing_clean_up
!***********************************************************************
endmodule Forcing
