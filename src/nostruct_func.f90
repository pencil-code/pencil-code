! $Id: nostruct_func.f90,v 1.1 2002-12-26 16:47:46 nilshau Exp $
!
module  struct_func
  !
  use Cdata
  !
  implicit none
  !
  contains

!***********************************************************************
    subroutine structure(f)
!
  real, dimension (mx,my,mz,mvar) :: f
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(ip>25)  print*,f(1,1,1,1)  !(to keep compiler happy)
end subroutine structure
!***********************************************************************

end module struct_func
