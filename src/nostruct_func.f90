! $Id$
!
module struct_func
!
  use Cparam
  use General, only: keep_compiler_quiet
!
  implicit none
!
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
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ivec)
      call keep_compiler_quiet(b_vec)
      call keep_compiler_quiet(variabl)
!
    endsubroutine structure
!***********************************************************************
endmodule struct_func
