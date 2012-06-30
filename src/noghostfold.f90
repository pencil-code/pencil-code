! $Id$
!
module GhostFold
!
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  private
!
  public :: fold_df, fold_f
!
  contains
!***********************************************************************
    subroutine fold_df(df,ivar1,ivar2)
!
!  06-May-2009/anders: dummy
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: ivar1,ivar2
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(ivar1,ivar2)
!
    endsubroutine fold_df
!***********************************************************************
    subroutine fold_f(f,ivar1,ivar2)
!
!  06-May-2009/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ivar1,ivar2)
!
    endsubroutine fold_f
!***********************************************************************
endmodule GhostFold
