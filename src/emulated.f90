!
!  Module emulating not universally available intrinsics
!
module Emulated 
!
  implicit none
!
  integer (KIND=4) :: inan
!  integer (KIND=8) :: inan 
  real (KIND=4), parameter :: rtype=0.
!  real (KIND=8), parameter :: rtype=0.d0
!
  data inan /Z'7FBFFFFF'/     !(KIND=4)
!  data inan /Z'7FF7FFFFFFFFFFFF'/      !(KIND=8)
! 
  real :: rnan
! 
 contains
!***********************************************************************
  subroutine initialize_nan()
!
    rnan=transfer(inan,rtype)
!
  endsubroutine initialize_nan
!**********************************************************************
  elemental function isnan(v)
!
! 23-Dec-10/MR: coded
! 15-Jan-11/MR: modified for portability; added alternative code for double precision
!
! emulates not yet standard function isnan

  logical          :: isnan
  real, intent(in) :: v
!
  integer (KIND=4), parameter :: itype=0
!  integer (KIND=8), parameter :: itype=0
  integer (KIND=4) :: iv
!  integer (KIND=8) :: iv
!
  iv=transfer(v,itype)
!
  isnan = iv==inan
!
  endfunction isnan
!***********************************************************************
endmodule emulated
