! $Id: nostruct_func.f90,v 1.4 2003-01-13 19:00:57 mee Exp $
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
  character (len=*) :: variabl
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(ip<=15) print*,'Use STRUCT_FUNC  = struct_func in Makefile.local'
  if(ip>25)  print*,f(1,1,1,1)  !(to keep compiler happy)
end subroutine structure
!***********************************************************************

end module struct_func
