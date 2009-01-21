! $Id$
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
    subroutine init_solid_cells(f,xx,yy,zz)
!
!  Initial conditions for cases where we have solid structures in the domain.
!  Typically the flow field is set such that we have no-slip conditions
!  at the solid structure surface.
! 
!  28-nov-2008/nils: coded
!
      use Cdata
      use Sub

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(xx,yy,zz)
!
    end subroutine init_solid_cells
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
    function in_solid_cell(part_pos )
!
!  Check if the position px,py,pz is within a colid cell
!
!  02-dec-2008/nils: coded
!
      logical :: in_solid_cell
      real, dimension(3) :: part_pos
!
      in_solid_cell=.false.
!
      call keep_compiler_quiet(part_pos )
!
    end function in_solid_cell
!***********************************************************************
    subroutine read_solid_cells_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) call keep_compiler_quiet(iostat)
      call keep_compiler_quiet(unit)
!
    endsubroutine read_solid_cells_init_pars
!***********************************************************************
    subroutine read_solid_cells_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) call keep_compiler_quiet(iostat)
      call keep_compiler_quiet(unit)
!
    endsubroutine read_solid_cells_run_pars
!***********************************************************************
    subroutine write_solid_cells_init_pars(unit)
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_solid_cells_init_pars
!***********************************************************************
    subroutine write_solid_cells_run_pars(unit)
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_solid_cells_run_pars
!***********************************************************************  
  endmodule Solid_Cells
