! $Id: streamlines.f90 17621 2011-09-05 07:05:20Z iomsn $
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************
module Streamlines
!
  use Cdata
  use Cparam
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  private 
!
  public :: tracers_prepare, wtracers, trace_streamlines
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
  endsubroutine wtracers
!***********************************************************************
endmodule Streamlines
