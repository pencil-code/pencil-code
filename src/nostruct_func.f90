! $Id: nostruct_func.f90,v 1.3 2003-01-10 14:00:27 nilshau Exp $
!
module  struct_func
  !
  use Cdata
  !
  implicit none
  !
  contains

!***********************************************************************
    subroutine structure(f,variabl)
!
  real, dimension (mx,my,mz,mvar) :: f
  character (len=4) :: variabl
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(ip<=15) print*,'Use STRUCT_FUNC  = struct_func in Makefile.local'
  if(ip>25)  print*,f(1,1,1,1)  !(to keep compiler happy)
end subroutine structure
!***********************************************************************

end module struct_func
