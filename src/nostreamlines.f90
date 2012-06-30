! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************
module Streamlines
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  public :: tracers_prepare, wtracers, trace_streamlines
  public :: read_streamlines_init_pars, write_streamlines_init_pars
  public :: read_streamlines_run_pars, write_streamlines_run_pars
!
  integer, public :: ntracers
!
  contains
!***********************************************************************
    subroutine trace_streamlines(f,tracers,n_tracers,h_max,h_min,l_max,tol,vv)
!
!  trace stream lines of the vetor field stored in f(:,:,:,iaa)
!
!   13-feb-12/simon: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: n_tracers
      real :: h_max, h_min, l_max, tol
      real, dimension (mx,my,7) :: tracers
      real, pointer, dimension (:,:,:,:) :: vv
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(vv)
      call keep_compiler_quiet(n_tracers)
      call keep_compiler_quiet(h_max,h_min,l_max,tol)
      call keep_compiler_quiet(tracers)
!
    endsubroutine trace_streamlines
!***********************************************************************
    subroutine tracers_prepare()
!
!  Dummy routine
!
    endsubroutine tracers_prepare
!***********************************************************************
  subroutine wtracers(f,path)
!
    real, dimension (mx,my,mz,mfarray) :: f
    character(len=*) :: path
!
    call keep_compiler_quiet(path)
    call keep_compiler_quiet(f)
!
  endsubroutine wtracers
!***********************************************************************
  subroutine read_streamlines_init_pars(unit,iostat)
!
    integer, intent(in) :: unit
    integer, intent(inout), optional :: iostat
!
    call keep_compiler_quiet(unit)
    if (present(iostat)) call keep_compiler_quiet(iostat)
!
  endsubroutine read_streamlines_init_pars
!***********************************************************************
  subroutine write_streamlines_init_pars(unit)
!
    integer, intent(in) :: unit
!
    call keep_compiler_quiet(unit)
!
  endsubroutine write_streamlines_init_pars
!***********************************************************************
  subroutine read_streamlines_run_pars(unit,iostat)
!
    integer, intent(in) :: unit
    integer, intent(inout), optional :: iostat
!
    call keep_compiler_quiet(unit)
    if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_streamlines_run_pars
!***********************************************************************
  subroutine write_streamlines_run_pars(unit)
!
    integer, intent(in) :: unit
!
    call keep_compiler_quiet(unit)
!
  endsubroutine write_streamlines_run_pars
!***********************************************************************
endmodule Streamlines
