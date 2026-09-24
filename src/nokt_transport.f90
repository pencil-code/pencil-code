! $Id$
!
!  Dummy counterpart of kt_transport.f90 (Kurganov-Tadmor flux-limited
!  transport for the conservative relativistic Higgsless hydro).  Selected
!  by default (KT_TRANSPORT = nokt_transport); see kt_transport.f90.
!
!  Kept intentionally minimal (Cparam + General only, no Messages) so it is
!  trivially compilable regardless of module build order, matching the
!  noweno_transport dummy contract.  If lkt_transport=T is set without
!  compiling the real module, kt_div returns `impossible`, which surfaces
!  immediately as a non-finite state -- the same fail-loud behaviour as the
!  other noXXX transport dummies.
!
!  02-sep-2026/Isak Stomberg: created
!
module KT_transport
!
  use Cparam, only: mx, my, mz, mfarray, impossible, ikind8
  use Quiet, only: keep_compiler_quiet
!
  implicit none
!
  include "kt_transport.h"
!
  contains
!***********************************************************************
    subroutine kt_init(irho_in,iux_in,ihless_in,eps_in,width_abs_in,theta_in)
!
!  02-sep-2026/Isak Stomberg: dummy
!
      integer, intent(in) :: irho_in, iux_in, ihless_in
      real, intent(in) :: eps_in, width_abs_in, theta_in
!
      call keep_compiler_quiet(irho_in,iux_in,ihless_in)
      call keep_compiler_quiet(eps_in,width_abs_in,theta_in)
!
    endsubroutine kt_init
!***********************************************************************
    subroutine kt_div(f,divergence)
!
!  02-sep-2026/Isak Stomberg: dummy
!  16-sep-2026/Isak Stomberg: renamed from kt_transp
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(:), intent(out) :: divergence
!
      divergence=impossible
      call keep_compiler_quiet(f)
!
    endsubroutine kt_div
!***********************************************************************
    subroutine kt_div_tensor(f,divergence)
!  09-sep-2026/Isak Stomberg: dummy
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(:,:), intent(out) :: divergence
      divergence=impossible
      call keep_compiler_quiet(f)
    endsubroutine kt_div_tensor
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=100
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    endsubroutine pushpars2c
!***********************************************************************
endmodule KT_transport
