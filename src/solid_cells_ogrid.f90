! $Id$
!
!  This module add solid (as in no-fluid) cells in the domain.
!  This can be used e.g. in order to simulate a cylinder in a cross flow.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lsolid_cells = .true.
!
!***************************************************************
!
!
!   New solid cells module
!   O-grid around cylinder or sphere, coupled to cartesian
!   grid outside the o-grid by interpolation.
!
!   Very fine resolution of boundary layer
!
!
module Solid_Cells

  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'solid_cells.h'
!
  integer, parameter            :: max_items=5
  integer                       :: ncylinders=0, dummy
  integer                       :: nobjects
  character(len=1)              :: flow_direction
  real, dimension(max_items)    :: cylinder_radius=0.
  real, dimension(max_items)    :: cylinder_temp=703.0
  real, dimension(max_items)    :: cylinder_xpos=0., cylinder_ypos=0.
  real, dimension(max_items)    :: cylinder_zpos=0.
  real :: skin_depth=0.

  type :: ogrid
    real, dimension(3) :: x
    real :: r_comp, r_inter
    real :: dr, dtheta
  endtype ogrid

  type solid_object
    character(len=10) :: form='cylinder'
    real :: r,T
    real, dimension(3) :: x
  endtype solid_object

  type (solid_object), dimension(max_items), save :: objects

  namelist /solid_cells_init_pars/ &
      cylinder_temp, ncylinders, cylinder_radius, cylinder_xpos, &
      cylinder_ypos, cylinder_zpos, flow_direction, skin_depth

  namelist /solid_cells_run_pars/  &
      dummy
  
  contains 
 !***********************************************************************
  subroutine register_solid_cells()
! 
!  Dummy routine
!
    end subroutine register_solid_cells
!***********************************************************************
subroutine initialize_grid()
end subroutine initialize_grid
!***********************************************************************
subroutine initialize_solid_cells(f)
!
  real, dimension(mx,my,mz,mfarray), intent(inout) :: f
  integer :: icyl, isph,i
  !
  !  Loop over all cylinders
  !
  do icyl = 1,ncylinders
    if (cylinder_radius(icyl) > 0) then
      objects(icyl)%r    = cylinder_radius(icyl)
      objects(icyl)%x(1) = cylinder_xpos(icyl)
      objects(icyl)%x(2) = cylinder_ypos(icyl)
      objects(icyl)%x(3) = cylinder_zpos(icyl)
      objects(icyl)%T    = cylinder_temp(icyl)
      objects(icyl)%form = 'cylinder'
    else
      call fatal_error('initialize_solid_cells', &
          'All cylinders must have non-zero radii!')
    endif
  enddo
  !
  ! Give all points inside the object a value of zero for all variables
  !
  !TODO
  !
end subroutine initialize_solid_cells
!***********************************************************************
subroutine init_solid_cells(f)
!
  real, dimension(mx,my,mz,mfarray), intent(in):: f
  call keep_compiler_quiet(f)
!
end subroutine init_solid_cells
!***********************************************************************
subroutine dsolid_dt(f,df,p)
  !
  real, dimension(mx,my,mz,mfarray), intent(in):: f
  real, dimension(mx,my,mz,mvar), intent(in)   :: df
  type (pencil_case), intent(in)                :: p
  !
  call keep_compiler_quiet(df,f)
  !
end subroutine dsolid_dt
!***********************************************************************
subroutine dsolid_dt_integrate
end subroutine dsolid_dt_integrate
!***********************************************************************
subroutine rprint_solid_cells(lreset,lwrite)
!
  logical :: lreset
  logical, optional :: lwrite
!
end subroutine rprint_solid_cells
!***********************************************************************
subroutine update_solid_cells(f)
!
!  Set the boundary values of the solid area such that we get a
!  correct fluid-solid interface.
!
!  19-nov-2008/nils: coded
!
  real, dimension(mx,my,mz,mfarray) :: f
!
end subroutine update_solid_cells
!***********************************************************************
subroutine update_solid_cells_pencil(f)
!
!  Set the boundary values of the solid area such that we get a
!  correct fluid-solid interface.
!
!  30-mar-15/Jørgen+nils: coded
!
  real, dimension(mx,my,mz,mfarray) :: f
!
end subroutine update_solid_cells_pencil
!***********************************************************************
    subroutine freeze_solid_cells(df)
!
!  If we are in a solid cell (or in a cell where the value of the variables are
!  found from interpolation) set df=0 for all variables
!
!  19-nov-2008/nils: coded
!
      real, dimension(mx,my,mz,mvar) :: df
      integer :: i
!
!      do i = l1,l2
!df(i,m,n,:) = 0
!      enddo
!
    endsubroutine freeze_solid_cells
!***********************************************************************
    function in_solid_cell(part_pos,part_rad)
!
!  Check if the position px,py,pz is within a solid cell
!
!  02-dec-2008/nils: coded
!TODO: This shouls be made more efficient (including the function call!)
!
      logical :: in_solid_cell
      real, dimension(3) :: obj_pos, part_pos
      real :: obj_rad, distance2, part_rad, rad_part
      integer :: iobj, i, ndims
!
      in_solid_cell = .false.
!
      do iobj = 1,nobjects
        obj_rad = objects(iobj)%r
        obj_pos = objects(iobj)%x(1:3)
        distance2 = 0
!
!  Loop only over the number of dimensions required
!
        ndims = 2
        if (objects(iobj)%form == 'sphere') ndims = 3
        do i = 1,ndims
          distance2 = distance2+(obj_pos(i)-part_pos(i))**2
        enddo
!
!  The object_skin is the closest a particle can get to the solid
!  cell before it is captured (this variable is normally zero).
!
        if (sqrt(distance2) < obj_rad+rad_part) then
          in_solid_cell = .true.
        endif
      enddo
!
    endfunction in_solid_cell
!***********************************************************************
    subroutine pencil_criteria_solid_cells()
!
!  All pencils that the Solid_Cells module depends on are specified here.
!
!  mar-2009/kragset: coded
!
!  Request p and sij-pencils here
!  Request rho-pencil
      lpenc_requested(i_pp) = .true.
      lpenc_requested(i_sij) = .true.
      lpenc_requested(i_rho) = .true.
!      if (idiag_Nusselt /= 0) lpenc_requested(i_gTT) = .true.
!
    endsubroutine pencil_criteria_solid_cells
!***********************************************************************
    subroutine solid_cells_clean_up()
!
!  Deallocate the variables allocated in solid_cells
!
!  7-oct-2010/dhruba: adeped from hydro_kinematic
!  21-jul-2011/bing: fixed, only deallocate variable if allocted
!
      print*, 'Deallocating some solid_cells variables ...'
!      if (allocated(fpnearestgrid)) deallocate(fpnearestgrid)
!      if (allocated(c_dragx)) deallocate(c_dragx)
!      if (allocated(c_dragy)) deallocate(c_dragy)
!      if (allocated(c_dragz)) deallocate(c_dragz)
!      if (allocated(c_dragx_p)) deallocate(c_dragx_p)
!      if (allocated(c_dragy_p)) deallocate(c_dragy_p)
!      if (allocated(c_dragz_p)) deallocate(c_dragz_p)
!      if (allocated(Nusselt)) deallocate(Nusselt)
      print*, '..Done.'
!
    endsubroutine solid_cells_clean_up
!***********************************************************************
subroutine recive_flow_info()
end subroutine recive_flow_info
!***********************************************************************
subroutine interpolate_cartesian_to_ogrid()
end subroutine interpolate_cartesian_to_ogrid
!***********************************************************************
subroutine timestep_solid()
end subroutine timestep_solid
!***********************************************************************
subroutine resolve_boundary()
end subroutine resolve_boundary
!***********************************************************************
subroutine resolve_flow()
end subroutine resolve_flow
!***********************************************************************
subroutine interpolate_ogrid_to_cartesian()
end subroutine interpolate_ogrid_to_cartesian
!***********************************************************************
subroutine send_flow_info()
end subroutine send_flow_info
!***********************************************************************
subroutine compute_draglift()
end subroutine compute_draglift
!***********************************************************************
subroutine print_solid()
end subroutine print_solid
!***********************************************************************
logical function within_ogrid_comp(xxp)
!
! Particles could check first if thery are within the o-grids 
! computational domain (closer than interpolation domain)
! when the particle positions are interpolated. If this is the 
! case they should use the o-grid for interpolation,
! but if they are outside, they should use the cartesian grid
! Need a variable R_grid for this, that is, the radius of the o-grid. 
! Mabye a different name than R_grid.
!
! 26-jan-2017/jorgen: coded
!
! Check if current possition is within the o-grids computational domain
!
  use Cdata ! or something, need to get solid and ogrid from somewhere
  real, dimension(3), intent(in) :: xxp
  
!  within_ogrid_comp = .false.
!
!  if(sum((xxp(1:2)-solid%x(1:2))**2) < ogrid%r_comp) then
!    within_ogrid_comp = .true.
!  endif
!  
endfunction within_ogrid_comp
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
    subroutine read_solid_cells_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      !read(parallel_unit, NML=solid_cells_init_pars, IOSTAT=iostat)
      iostat = 0
      read(parallel_unit, NML=solid_cells_init_pars)
!
    endsubroutine read_solid_cells_init_pars
!***********************************************************************
    subroutine write_solid_cells_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=solid_cells_init_pars)
!
    endsubroutine write_solid_cells_init_pars
!***********************************************************************
    subroutine read_solid_cells_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=solid_cells_run_pars, IOSTAT=iostat)
!
    endsubroutine read_solid_cells_run_pars
!***********************************************************************
    subroutine write_solid_cells_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=solid_cells_run_pars)
!
    endsubroutine write_solid_cells_run_pars
!***********************************************************************
    subroutine close_interpolation(f,ix0_,iy0_,iz0_,iobj,xxp,f_tmp, fluid_point)
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer, intent(in) :: ix0_, iy0_, iz0_, iobj
      real, dimension(mvar), intent(inout) :: f_tmp
      real, dimension(3), intent(in) :: xxp
      logical, intent(in) :: fluid_point
!
    end subroutine close_interpolation
!***********************************************************************
  end module Solid_Cells
