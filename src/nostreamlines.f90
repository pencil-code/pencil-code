! $Id: streamlines.f90 17621 2011-09-05 07:05:20Z iomsn $
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM integer, parameter :: nghost = 3
!
!***************************************************************
module Streamlines
!
  use Cdata
  use Cparam
  use Messages
!
  implicit none
  
  contains
!*********************************************************************** 
    subroutine trace_streamlines(f,iaa,tracers)
!
!  trace stream lines of the vetor field stored in f(:,:,:,iaa)
!
!   13-feb-12/simon: coded

      real, dimension (mx,my,mz,mfarray) :: f
      integer :: iaa, l, m
      real, dimension (mx,my,7) :: tracers
!
      intent(in)  :: f,iaa,tracers
      
      call keep_compiler_quiet(f)
!
    endsubroutine trace_streamlines
!***********************************************************************
endmodule Streamlines
