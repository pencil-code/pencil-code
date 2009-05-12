! $Id$
!
!  This module add solid (as in no-fluid) cells in the domain.
!  This can be used e.g. in order to simulate a cylinder in a cross flow.
!
module Solid_Cells
!
  use Cdata
  use Cparam
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'solid_cells.h'
!
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
    subroutine init_solid_cells(f)
!
!  Initial conditions for cases where we have solid structures in the domain.
!  Typically the flow field is set such that we have no-slip conditions
!  at the solid structure surface.
! 
!  28-nov-2008/nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_solid_cells
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
    function in_solid_cell(part_pos,part_rad)
!
!  Check if the position px,py,pz is within a colid cell
!
!  02-dec-2008/nils: coded
!
      logical :: in_solid_cell
      real, dimension(3) :: part_pos
      real :: part_rad
!
      in_solid_cell=.false.
!
      call keep_compiler_quiet(part_rad)
      call keep_compiler_quiet(part_pos)
!
    endfunction in_solid_cell
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
    subroutine close_interpolation(f,ix0,iy0,iz0,icyl,ivar1,xxp,gpp,&
        fluid_point)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: ix0,iy0,iz0,icyl,ivar1
      real, intent(inout) :: gpp
      real, dimension(3), intent(in) :: xxp
      logical, intent(in) :: fluid_point
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ix0)
      call keep_compiler_quiet(iy0)
      call keep_compiler_quiet(iz0)
      call keep_compiler_quiet(icyl)
      call keep_compiler_quiet(ivar1)
      call keep_compiler_quiet(xxp)
      call keep_compiler_quiet(gpp)
      call keep_compiler_quiet(fluid_point)
!
    endsubroutine close_interpolation
!***********************************************************************  
  endmodule Solid_Cells
