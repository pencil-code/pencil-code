! $Id$
!
!  This module is the dummy for the opacity module(s).
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lopacity = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***********************************************************************
module Opacity
!
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'opacity.h'

  contains
!***********************************************************************
      subroutine initialize_opacity
!
      endsubroutine initialize_opacity
!***********************************************************************
      subroutine read_opacity_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      call keep_compiler_quiet(iostat)

      endsubroutine read_opacity_run_pars
!***********************************************************************
      subroutine write_opacity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)

      endsubroutine write_opacity_run_pars
!***********************************************************************
      subroutine get_opacity(tt,rho,kappa,dkap_dtt,dkap_drho)

      real, dimension(:), intent(in) :: tt, rho
      real, dimension(:), intent(out):: kappa,dkap_dtt,dkap_drho

      call keep_compiler_quiet(tt)
      call keep_compiler_quiet(rho)
      call keep_compiler_quiet(kappa)
      call keep_compiler_quiet(dkap_dtt)
      call keep_compiler_quiet(dkap_drho)

      endsubroutine get_opacity
!***********************************************************************
end module 
