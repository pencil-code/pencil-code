! $Id: nopower_spectrum.f90,v 1.5 2003-06-16 04:41:11 brandenb Exp $
!
module  power_spectrum
  !
  use Cdata
  !
  implicit none
  !
  contains

!***********************************************************************
    subroutine power(f,sp)
!
  real, dimension (mx,my,mz,mvar+maux) :: f
  character (len=1) :: sp
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(sp=='') print*,f(1,1,1,1)  !(to keep compiler happy)
  endsubroutine power
!***********************************************************************
  subroutine powerhel(f,sp)
!
  real, dimension (mx,my,mz,mvar+maux) :: f
  character (len=3) :: sp
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(sp=='') print*,f(1,1,1,1)  !(to keep compiler happy)
  endsubroutine powerhel
!***********************************************************************
    subroutine power_1d(f,sp,ivec)
!
  real, dimension (mx,my,mz,mvar+maux) :: f
  character (len=1) :: sp
  integer :: ivec
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(sp=='') print*,f(1,1,1,1),ivec  !(to keep compiler happy)
  endsubroutine power_1d
!***********************************************************************

end module power_spectrum
