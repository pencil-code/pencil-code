! $Id$
!
! MODULE_DOC: This module contains Yin-Yang related dummy types and functions.
!
!**************************************************************************
!
module Yinyang
!
  use Cdata, only: lroot
  use General, only: keep_compiler_quiet

  implicit none

  include 'yinyang.h'

  type ind_coeffs
    integer, dimension(1,1,1) :: inds
    real, dimension(1,1,1) :: coeffs
  end type

contains
!
!**************************************************************************
    subroutine bilin_interp(indcoeffs, ith, iph, f, buffer, i2buf, i3buf)
!
!  Dummy routine.
!
!  20-dec-15/MR: coded
! 
      type(ind_coeffs),         intent(IN) :: indcoeffs
      integer,                  intent(IN) :: ith, iph, i2buf, i3buf
      real, dimension(:,:,:,:), intent(IN) :: f
      real, dimension(:,:,:,:), intent(OUT):: buffer

      if (lroot) &
        print*, 'bilin_interp: not implemented in Fortran 95'
      stop

      buffer=0.

    endsubroutine bilin_interp
!**************************************************************************
    subroutine biquad_interp(indcoeffs, ith, iph, f, buffer, i2buf, i3buf)
!
!  Dummy routine.
!
!  20-dec-15/MR: coded
! 
      type(ind_coeffs),         intent(IN) :: indcoeffs
      integer,                  intent(IN) :: ith, iph, i2buf, i3buf
      real, dimension(:,:,:,:), intent(IN) :: f
      real, dimension(:,:,:,:), intent(OUT):: buffer

      if (lroot) &
        print*, 'bilin_interp: not implemented in Fortran 95'
      stop

      buffer=0.

    endsubroutine biquad_interp
!***********************************************************************
    function prep_interp(thphprime,indcoeffs,itype,n_onegap,n_twogap,th_range) result (nok)
!
!  Dummy routine.
!
!  20-dec-15/MR: coded
!
      real, dimension(:,:,:),          intent(IN) :: thphprime
      type(ind_coeffs),                intent(OUT):: indcoeffs
      integer,                         intent(IN) :: itype
      integer,               optional, intent(IN) :: n_onegap,n_twogap
      integer, dimension(2), optional, intent(OUT):: th_range

      integer :: nok

      nok=0
      if (lroot) &
        print*, 'prep_interp: not implemented in Fortran 95'
      stop

      indcoeffs%inds=0
      indcoeffs%coeffs=0.
      if (present(th_range)) th_range=0
!
    endfunction prep_interp
!**************************************************************************
    subroutine coeffs_to_weights(intcoeffs,indweights)
!
!  Dummy routine.
!
!   20-mar-16/MR: coded
!
      type(ind_coeffs), intent(IN) :: intcoeffs
      type(ind_coeffs), intent(OUT):: indweights

      if (lroot) &
        print*, 'coeffs_to_weights: not implemented in Fortran 95'
      stop

      indweights%inds=0
      indweights%coeffs=0.

    endsubroutine coeffs_to_weights
!*******************************************************************
    function in_overlap_mask(indth,indph) result (ok)

      integer, intent(IN) :: indth,indph
      logical :: ok

      ok=.true.

    endfunction in_overlap_mask
!*******************************************************************
endmodule Yinyang
