! $Id$
!
!  This module take care of WENO (weighted essentially non oscillatory)
!  transport.
!
module WENO_transport
!
  implicit none
!
  private
!
  public :: weno_transp
!
  contains
!***********************************************************************
    subroutine weno_transp(fq,m,n,iq,iq1,iux,iuy,iuz,dq,dx_1,dy_1,dz_1,ref,ref1)
!
!  Solve the equation dq/dt + div(u*q) = 0 using the WENO method.
!
!  29-dec-09/evghenii+anders: dummy
!
      use Cparam, only: impossible
      use General, only: keep_compiler_quiet
!
      real, dimension(:,:,:,:), intent(in ) :: fq
      integer, intent(in) :: m, n
      integer, intent(in) :: iq, iq1, iux, iuy, iuz
      real, dimension(:), intent(out) :: dq
      real, dimension(:), intent(in)  :: dx_1, dy_1, dz_1
      real, dimension(:), optional, intent(in) :: ref,ref1
!
      dq = impossible
!
      call keep_compiler_quiet(fq)
      call keep_compiler_quiet(m,n,iq,iq1)
      call keep_compiler_quiet(iux,iuy,iuz)
      call keep_compiler_quiet(dx_1,dy_1,dz_1)
      call keep_compiler_quiet(dq,ref,ref1)
!
    endsubroutine weno_transp
!***********************************************************************
endmodule WENO_transport
