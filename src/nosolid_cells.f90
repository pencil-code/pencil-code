! $Id$
!
!  This module add solid (as in no-fluid) cells in the domain.
!  This can be used e.g. in order to simulate a cylinder in a cross flow.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lsolid_cells = .false.
!
!***************************************************************
module Solid_Cells
!
  use Cdata
  use Cparam
  use Sub, only: keep_compiler_quiet
  use messages, only: warning
!
  implicit none
!
  include 'solid_cells.h'
!
  integer :: idiag_c_dragx=0       ! DIAG_DOC:
  integer :: idiag_c_dragy=0       ! DIAG_DOC:
  integer :: idiag_c_dragz=0       ! DIAG_DOC:
!
  contains
!***********************************************************************
    subroutine initialize_solid_cells(f)
!
!  19-nov-08/nils: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
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
      real    :: part_rad
!
      in_solid_cell=.false.
!
      call keep_compiler_quiet(part_pos)
      call keep_compiler_quiet(part_rad)
!
    endfunction in_solid_cell
!***********************************************************************
    subroutine read_solid_cells_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_solid_cells_init_pars
!***********************************************************************
    subroutine read_solid_cells_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_solid_cells_run_pars
!***********************************************************************
    subroutine write_solid_cells_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_solid_cells_init_pars
!***********************************************************************
    subroutine write_solid_cells_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_solid_cells_run_pars
!***********************************************************************
    subroutine dsolid_dt(f,df,p)
!
      real, dimension (mx,my,mz,mfarray):: f
      real, dimension (mx,my,mz,mvar):: df
      type (pencil_case):: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dsolid_dt
!***********************************************************************
    subroutine dsolid_dt_integrate
!
! dummy routine
!
    endsubroutine dsolid_dt_integrate
!***********************************************************************
    subroutine rprint_solid_cells(lreset,lwrite)
!
!  reads and registers print parameters relevant for solid cells
!  dummy routine
!
!   mar-2009/kragset: coded
!
      use Diagnostics, only: parse_name
!
      logical :: lreset,lwr
      logical, optional :: lwrite
      integer :: iname
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_c_dragx=0
        idiag_c_dragy=0
        idiag_c_dragz=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'c_dragx',idiag_c_dragx)
        call parse_name(iname,cname(iname),cform(iname),'c_dragy',idiag_c_dragy)
        call parse_name(iname,cname(iname),cform(iname),'c_dragz',idiag_c_dragz)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
!
      endif
!
    endsubroutine rprint_solid_cells
!***********************************************************************
    subroutine pencil_criteria_solid_cells()
!
!  All pencils that the Solid_Cells module depends on are specified here.
!
!  mar-2009/kragset: dummy
!
    endsubroutine pencil_criteria_solid_cells
!***********************************************************************
    subroutine close_interpolation(f,ix0,iy0,iz0,icyl,xxp,gpp,&
        fluid_point,lnew_interpolation_method)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: ix0,iy0,iz0,icyl
      real, dimension(mvar), intent(inout) :: f_tmp
      real, dimension(3), intent(in) :: xxp
      logical, intent(in) :: fluid_point,lnew_interpolation_method
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ix0)
      call keep_compiler_quiet(iy0)
      call keep_compiler_quiet(iz0)
      call keep_compiler_quiet(icyl)
      call keep_compiler_quiet(xxp)
      call keep_compiler_quiet(f_tmp)
      call keep_compiler_quiet(fluid_point)
      call keep_compiler_quiet(lnew_interpolation_method)
!
    endsubroutine close_interpolation
!***********************************************************************
    subroutine solid_cells_clean_up
!
!  This is a dummy routine.
!
!  7-oct-2010/dhruba: coded
!
      call warning('solid_cells_clean_up','dont call in nosolid_cells')
!
    endsubroutine solid_cells_clean_up
!***********************************************************************
endmodule Solid_Cells
