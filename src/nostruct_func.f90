! $Id: nostruct_func.f90,v 1.2 2003-01-06 20:23:23 nilshau Exp $
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
  if(ip<=15) print*,'Use STRUCT_FUNC  = struct_func in Makefile.local'
  if(ip>25)  print*,f(1,1,1,1)  !(to keep compiler happy)
end subroutine structure
!***********************************************************************

end module struct_func
