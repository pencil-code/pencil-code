! $Id$
!
!  Dummy
!
!***************************************************************
module ImplicitDiffusion
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: fatal_error
!
  implicit none
!
  include 'implicit_diffusion.h'
!
  contains
!***********************************************************************
    subroutine integrate_diffusion(get_diffus_coeff, f, ivar1, ivar2)
!
!  Dummy
!
      external :: get_diffus_coeff
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar1
      integer, intent(in), optional :: ivar2
!
      call fatal_error('integrate_diffusion', 'ImplicitDiffusion module is not plugged in. ')
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ivar1)
      if (present(ivar2)) call keep_compiler_quiet(ivar2)
!
    endsubroutine integrate_diffusion
!***********************************************************************
endmodule ImplicitDiffusion
