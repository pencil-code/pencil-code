! $Id$
!
!  Module for boundary conditions. Extracted from (no)mpicomm, since
!  all non-periodic (external) boundary conditions require the same
!  code for serial and parallel runs.
!
module Boundcond
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Mpicomm
!
  implicit none
!
  private
!
  public :: update_ghosts
  public :: boundconds, boundconds_x, boundconds_y, boundconds_z
  public :: bc_per_x, bc_per_y, bc_per_z
!
  contains
!***********************************************************************
    subroutine update_ghosts(a)
!
!  23-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: a
!
      call keep_compiler_quiet(a)
!
    endsubroutine update_ghosts
!***********************************************************************
    subroutine boundconds(f,ivar1_opt,ivar2_opt)
!
!  23-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      call keep_compiler_quiet(f)
      if (present(ivar1_opt)) call keep_compiler_quiet(ivar1_opt)
      if (present(ivar2_opt)) call keep_compiler_quiet(ivar2_opt)
!
    endsubroutine boundconds
!***********************************************************************
    subroutine boundconds_x(f,ivar1_opt,ivar2_opt)
!
!  23-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      call keep_compiler_quiet(f)
      if (present(ivar1_opt)) call keep_compiler_quiet(ivar1_opt)
      if (present(ivar2_opt)) call keep_compiler_quiet(ivar2_opt)
!
    endsubroutine boundconds_x
!***********************************************************************
    subroutine boundconds_y(f,ivar1_opt,ivar2_opt)
!
!  23-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      call keep_compiler_quiet(f)
      if (present(ivar1_opt)) call keep_compiler_quiet(ivar1_opt)
      if (present(ivar2_opt)) call keep_compiler_quiet(ivar2_opt)
!
    endsubroutine boundconds_y
!***********************************************************************
    subroutine boundconds_z(f,ivar1_opt,ivar2_opt)
!
!  23-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      call keep_compiler_quiet(f)
      if (present(ivar1_opt)) call keep_compiler_quiet(ivar1_opt)
      if (present(ivar2_opt)) call keep_compiler_quiet(ivar2_opt)
!
    endsubroutine boundconds_z
!***********************************************************************
    subroutine bc_per_x(f,topbot,j)
!
!  23-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      character (len=3) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(j)
!
    endsubroutine bc_per_x
!***********************************************************************
    subroutine bc_per_y(f,topbot,j)
!
!  23-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      character (len=3) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(j)
!
    endsubroutine bc_per_y
!***********************************************************************
    subroutine bc_per_z(f,topbot,j)
!
!  23-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      character (len=3) :: topbot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(j)
!
    endsubroutine bc_per_z
!***********************************************************************
endmodule Boundcond
