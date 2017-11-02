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
! CPARAM logical, parameter :: lsolid_ogrid = .false.
!
!***************************************************************
module Solid_Cells
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: warning
!
  implicit none
!
  include 'solid_cells.h'
!
  real :: r_ogrid
  real :: r_int_outer
  real, dimension(3) :: xorigo_ogrid
!
  contains
!***********************************************************************
    subroutine register_solid_cells
!
!  Dummy routine
!
    end subroutine register_solid_cells
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
    subroutine update_solid_cells_pencil(f)
!
!  Set the boundary values of the solid area such that we get a
!  correct fluid-solid interface.
!
!  30-mar-15/Jørgen+nils: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine update_solid_cells_pencil
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
    subroutine read_solid_cells_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_solid_cells_init_pars
!***********************************************************************
    subroutine write_solid_cells_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_solid_cells_init_pars
!***********************************************************************
    subroutine read_solid_cells_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_solid_cells_run_pars
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
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_solid_cells
!***********************************************************************
    subroutine pencil_criteria_solid_cells
!
!  All pencils that the Solid_Cells module depends on are specified here.
!
!  mar-2009/kragset: dummy
!
    endsubroutine pencil_criteria_solid_cells
!***********************************************************************
    subroutine close_interpolation(f,ix0,iy0,iz0,icyl,xxp,f_tmp,&
        fluid_point)
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: ix0,iy0,iz0,icyl
      real, dimension(mvar), intent(inout) :: f_tmp
      real, dimension(3), intent(in) :: xxp
      logical, intent(in) :: fluid_point
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ix0)
      call keep_compiler_quiet(iy0)
      call keep_compiler_quiet(iz0)
      call keep_compiler_quiet(icyl)
      call keep_compiler_quiet(xxp)
      call keep_compiler_quiet(f_tmp)
      call keep_compiler_quiet(fluid_point)
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
    subroutine solid_cells_timestep_first(f)
!
!  Setup dfs in the beginning of each itsub.
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    end subroutine solid_cells_timestep_first
!***********************************************************************
    subroutine solid_cells_timestep_second(f,int_dt,int_ds)
!
!  Time evolution of solid_cells variables.
!
      real, dimension(mx,my,mz,mfarray) :: f
      real :: int_dt, int_ds
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(int_dt)
      call keep_compiler_quiet(int_ds)
!
    endsubroutine solid_cells_timestep_second
!***********************************************************************
    subroutine time_step_ogrid(f)
!
!  Dummy routine
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    end subroutine time_step_ogrid
!***********************************************************************
    subroutine wsnap_ogrid(chsnap,enum,flist)
!
!  Dummy routine
!
      character(len=*), intent(in) :: chsnap
      character(len=*), intent(in), optional :: flist
      logical, intent(in), optional :: enum
!
      if (ALWAYS_FALSE) print*, chsnap
      if (ALWAYS_FALSE) then 
        if (present(enum)) print*, enum
        if (present(enum)) print*, flist
      endif
    endsubroutine wsnap_ogrid
!***********************************************************************
    subroutine map_nearest_grid_ogrid(xxp,ineargrid_ogrid,rthz)
!
!  Dummy routine
!
      real, dimension (3) :: xxp,rthz
      integer, dimension (4) :: ineargrid_ogrid
!
      intent(in)  :: xxp,rthz
      intent(out) :: ineargrid_ogrid
!      
      if(ALWAYS_FALSE) print*, xxp,rthz
      ineargrid_ogrid=0
!      
    endsubroutine map_nearest_grid_ogrid
!***********************************************************************
    subroutine interpolate_particles_ogrid(ivar1,ivar2,xxp,gp,inear_glob)
!
!  Dummy routine
!
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (4) :: inear_glob
!
      intent(in)  :: ivar1, ivar2, xxp, inear_glob
      intent(out) :: gp
!
      if (ALWAYS_FALSE) then
        print*, ivar1,ivar2,xxp
        print*, inear_glob
      endif
      gp=0.
!
    endsubroutine interpolate_particles_ogrid
!***********************************************************************
endmodule Solid_Cells
