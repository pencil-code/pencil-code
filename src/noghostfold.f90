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
  public :: fold_df, fold_f, fold_df_3points
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
    subroutine fold_df_3points(df,ivar1,ivar2)
!
!  12-Jan-2015/jonas: dummy
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: ivar1,ivar2
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(ivar1,ivar2)
!
    endsubroutine fold_df_3points
!***********************************************************************
endmodule GhostFold
