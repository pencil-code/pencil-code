! $Id: nostruct_func.f90,v 1.5 2003-01-20 16:56:53 nilshau Exp $
!
module  struct_func
  !
  use Cdata
  !
  implicit none
  !
  contains

!***********************************************************************
    subroutine structure(f,ivec,b_vec,variabl)
!
  real, dimension (mx,my,mz,mvar) :: f
  real, dimension (nx,ny,nz) :: b_vec
  integer :: ivec
  character (len=*) :: variabl
  !
  if(ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if(ip<=15) print*,'Use STRUCT_FUNC  = struct_func in Makefile.local'
  if(ip>25)  print*,f(1,1,1,1)  !(to keep compiler happy)
end subroutine structure
!***********************************************************************

end module struct_func
