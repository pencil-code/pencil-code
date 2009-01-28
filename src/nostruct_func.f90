! $Id$
!
module  struct_func
  !
  use Cdata
  !
  implicit none

  include 'struct_func.h'
  !
  contains

!***********************************************************************
    subroutine structure(f,ivec,b_vec,variabl)
!
  real, dimension (mx,my,mz,mfarray) :: f
  real, dimension (nx,ny,nz) :: b_vec
  integer :: ivec
  character (len=*) :: variabl
  !
  if (ip<=15) print*,'Use POWER=power_spectrum in Makefile.local'
  if (ip<=15) print*,'Use STRUCT_FUNC  = struct_func in Makefile.local'
  if (NO_WARN)  print*,f,ivec,b_vec,variabl  !(to keep compiler happy)
end subroutine structure
!***********************************************************************

endmodule struct_func
