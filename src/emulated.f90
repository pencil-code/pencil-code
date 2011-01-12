!
!  Module containing emulating not universally available intrinsics
!
module Emulated 
!
  use Cparam
!
  implicit none
!
  contains
!***********************************************************************
  elemental function isnan(v)

! 23-Dec-10: MR coded
! emulates not yet standard isnan

  logical          :: isnan
  real, intent(in) :: v
  real             :: vl
  integer          :: iv
  equivalence(vl,iv)

  vl=v

  isnan = iv==inan

  endfunction isnan
!***********************************************************************
endmodule emulated
