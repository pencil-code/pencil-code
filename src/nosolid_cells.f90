! $Id: timestep.f90 9840 2008-09-05 07:29:37Z ajohan $
!
!  This module add solid (as in no-fluid) cells in the domain.
!  This can be used e.g. in order to simulate a cylinder in a cross flow.
!
module Solid_Cells

  use Cparam
  use Cdata
  use Sub, only: keep_compiler_quiet

  implicit none

  include 'solid_cells.h'

  contains
!***********************************************************************
    subroutine initialize_solid_cells
!
!  19-nov-08/nils: dummy
!
      lsolid_cells=.false.
!
    endsubroutine initialize_solid_cells
!***********************************************************************  
    subroutine update_solid_cells(f)
!
!  Set the boundary values of the solid area such that we get a 
!  correct fluid-solid interface.
!
!  19-nov-08/nils: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,k,idir
!
      call keep_compiler_quiet(f)
!
    endsubroutine update_solid_cells
!***********************************************************************  
    subroutine freeze_solid_cells(df)
!
!  If we are in a solid cell set df=0 for all variables
!
!  19-nov-08/nils: dummy
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: i,j,k
! 
      call keep_compiler_quiet(df)
!
    endsubroutine freeze_solid_cells
!***********************************************************************  
  endmodule Solid_Cells
