! $Id: nopower_spectrum.f90,v 1.1 2002-11-19 09:10:00 brandenb Exp $
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
  real, dimension (mx,my,mz,mvar) :: f
  character (len=1) :: sp
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(sp=='') print*,f(1,1,1,1)  !(to keep compiler happy)
  endsubroutine power
!***********************************************************************
  subroutine powerhel(f,sp)
!
  real, dimension (mx,my,mz,mvar) :: f
  character (len=3) :: sp
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(sp=='') print*,f(1,1,1,1)  !(to keep compiler happy)
  endsubroutine powerhel
!***********************************************************************

end module power_spectrum
