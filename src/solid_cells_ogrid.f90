! d$
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

! Remove:
! TODO TODO
! THE LINES BELOW ARE TAKEN FROM GRID.F90
! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! PENCILS PROVIDED x_mn; y_mn; z_mn; r_mn; r_mn1
! PENCILS PROVIDED phix; phiy
! PENCILS PROVIDED pomx; pomy
! PENCILS PROVIDED rcyl_mn; rcyl_mn1; phi_mn
! PENCILS PROVIDED evr(3); rr(3); evth(3)
!
!***************************************************************
! TODO TODO
!***********************************************************************
! /Remove:
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

! PARAMETERS NECESSARY FOR GRID CONSTRUCTION 
!Todo
  integer, parameter :: nxgrid_ogrid=128,nygrid_ogrid=128,nzgrid_ogrid=128 ! start.in
  !integer, parameter :: nxgrid_ogrid,nygrid_ogrid,nzgrid_ogrid ! start.in
  real, parameter :: r_cylinder=0.1                            ! start.in 
  real, parameter :: r_ogrid=2*r_cylinder                      ! start.in
  character (len=labellen), dimension(3) :: grid_func_ogrid='linear' ! start.in


  integer, parameter :: nx_ogrid=nxgrid_ogrid/nprocx,ny_ogrid=nygrid_ogrid/nprocy,nz_ogrid=nzgrid_ogrid/nprocz
  integer, parameter :: mx_ogrid=nx_ogrid+2*nghost,l1_ogrid=1+nghost,l2_ogrid=mx_ogrid-nghost
  integer, parameter :: my_ogrid=ny_ogrid+2*nghost,m1_ogrid=1+nghost,m2_ogrid=my_ogrid-nghost
  integer, parameter :: mz_ogrid=nz_ogrid+2*nghost,n1_ogrid=1+nghost,n2_ogrid=mz_ogrid-nghost
  integer, parameter :: mxgrid_ogrid=nxgrid_ogrid+2*nghost
  integer, parameter :: mygrid_ogrid=nygrid_ogrid+2*nghost
  integer, parameter :: mzgrid_ogrid=nzgrid_ogrid+2*nghost
  real, dimension (mx_ogrid) :: x_ogrid,dx_1_ogrid,dx2_ogrid,dx_tilde_ogrid,xprim_ogrid
  real, dimension (my_ogrid) :: y_ogrid,dy_1_ogrid,dy2_ogrid,dy_tilde_ogrid,yprim_ogrid
  real, dimension (mz_ogrid) :: z_ogrid,dz_1_ogrid,dz2_ogrid,dz_tilde_ogrid,zprim_ogrid
  real, dimension (mx_ogrid) :: dVol_x_ogrid,dVol1_x_ogrid
  real, dimension (my_ogrid) :: dVol_y_ogrid,dVol1_y_ogrid
  real, dimension (mz_ogrid) :: dVol_z_ogrid,dVol1_z_ogrid
  real, dimension (nz_ogrid) :: dxyz_2_ogrid,dxyz_4_ogrid,dxyz_6_ogrid,dVol_ogrid
  real :: dx_ogrid,dy_ogrid,dz_ogrid
  real :: dxmin_ogrid,dxmax_ogrid
  
  logical, dimension(3) :: lequidist_ogrid
  real, dimension(3) :: xyz_star_ogrid
  real, dimension(-nghost:nghost,2) :: coeffs_1_x_ogrid
  real, dimension(-nghost:nghost,2) :: coeffs_1_y_ogrid
  real, dimension(-nghost:nghost,2) :: coeffs_1_z_ogrid
  real, dimension(0:nprocx) :: procx_bounds_ogrid
  real, dimension(0:nprocy) :: procy_bounds_ogrid
  real, dimension(0:nprocz) :: procz_bounds_ogrid
  real :: box_volume_ogrid=1.0 
  real :: r_int_ogrid=0.,r_ext_ogrid=impossible
  real, dimension(3) :: coeff_grid_o=1.0

  real, dimension (nx_ogrid) :: rcyl_mn_ogrid,rcyl_mn1_ogrid,rcyl_mn2_ogrid,rcyl_weight_ogrid
  integer :: lpoint_ogrid=(mx_ogrid+1)/2, mpoint_ogrid=(my_ogrid+1)/2, npoint_ogrid=(mz_ogrid+1)/2
  integer :: lpoint2_ogrid=(mx_ogrid+1)/4,mpoint2_ogrid=(my_ogrid+1)/4,npoint2_ogrid=(mz_ogrid+1)/4
  real, dimension(nxgrid_ogrid) :: xgrid_ogrid, dx1grid_ogrid, dxtgrid_ogrid
  real, dimension(nygrid_ogrid) :: ygrid_ogrid, dy1grid_ogrid, dytgrid_ogrid
  real, dimension(nzgrid_ogrid) :: zgrid_ogrid, dz1grid_ogrid, dztgrid_ogrid
  real, dimension(mxgrid_ogrid) :: xglobal_ogrid
  real, dimension(mygrid_ogrid) :: yglobal_ogrid
  real, dimension(mzgrid_ogrid) :: zglobal_ogrid
  real, dimension (-nghost:nghost) :: dx2_bound_ogrid=0., dy2_bound_ogrid=0., dz2_bound_ogrid=0.
  real, dimension (nx_ogrid) :: dxmax_pencil_ogrid,dxmin_pencil_ogrid
  real, dimension (nx_ogrid,3) :: dline_1_ogrid
  
  ! For mn-loop over ogrid pencils
  integer :: n_ogrid, m_ogrid
  ! For periodic boundaries
  integer, parameter :: l1i_ogrid=l1_ogrid+nghost-1
  integer :: l2i_ogrid=mx_ogrid-2*nghost+1
  integer, parameter :: m1i_ogrid=m1_ogrid+nghost-1
  integer :: m2i_ogrid=my_ogrid-2*nghost+1
  integer, parameter :: n1i_ogrid=n1_ogrid+nghost-1
  integer :: n2i_ogrid=mz_ogrid-2*nghost+1

  contains 
!***********************************************************************
  subroutine register_solid_cells()
! 
!  Dummy routine
!
  end subroutine register_solid_cells
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
! ROUTINE NOT NECESSARY IF ROUTINES FROM GRID.F90 AND HYDRO.f90 ARE USED
!    subroutine metrics(xi,eta,zeta)
!      ! Compute metrics on the selected grid used around the solid 
!      ! Metrics are the computational coordinates differentiated wrt physical space
!      ! Ref Sec.5.6.2 in Tannehill,Anderson,Pletcher(1997),2nd Ed.
!      real,dimension(mxi),  intent(out) :: xi
!      real,dimension(meta), intent(out) :: eta
!      real,dimension(mzeta),intent(out) :: zeta
!      !real,dimension(:), allocatable :: dxdxi,dxdeta,dxdzeta,dydxi,dydeta,dydzeta,dzdxi,dzdeta,dzdzeta
!      real, parameter :: twopi=2*pi
!      real, parameter :: dr_grid=r_grid-r_solid
!      integer :: i
!
!      ! Set up uniformly spaced grid in computational space
!      do i=1,mxi
!        xi=i
!      enddo
!      do i=1,meta
!        eta=i
!      enddo
!      do i=1,mzeta
!        mzeta=i
!      enddo
!
!      !!TODO: Compine this with cylindrical grid by updating metrics
!      !if(streched_from_center) then
!      !  eta=dr_grid*((beta+1)-(beta-1)*((beta+1)/(beta-1))**(1-eta))//(((beta+1)/(beta-1))**(1-eta)+1)
!      !endif
!
!
!      if(cylinder_grid) then
!        ! _J indicates multiplied by Jacobi determinant
!        jacobi=1./(twopi*dr_grid*(r_solid+dr_grid*eta))
!        dxidx_J=dr_grid*sin(twopi*xi)
!        dxidy_J=dr_grid*cos(twopi*xi)
!        dxidz=0.
!        detadx=cos(twopi*xi)/dr_grid
!        detady=sin(twopi*xi)/dr_grid
!        detadz=0.
!        dzetadx=0.
!        dzetady=0.
!        dzetadz=1.
!      endif
!
!    end subroutine metrics
!
    subroutine construct_grid_ogrid
!
!  Constructs a non-equidistant cylindrical grid x_ogrid that surrounds the solid cylinder 
!  and partially overlaps with the cartesian grid x.
!  The grid x_ogrid(xi) is constructed from an equidistant grid xi with grid spacing dx=1.
!  grid spacing dxi=1. For grid_func_ogrid='linear' this is equivalent to an
!  equidistant grid. 
!
!  Grid is defined as (x_ogrid) = (r,theta,z)
!  Must be periodic and equidistant in theta-direction
!  Cannot be periodic in r-direction
!
!  dx_1_o and dx_tilde_o are the coefficients that enter the formulae for the
!  1st and 2nd derivative:
!
!    ``df/dx_o'' = ``df/dxi'' * dx_1_o
!    ``d2f/dx2_o'' = ``df2/dxi2'' * dx_1_o**2 + dx_tilde_o * ``df/dxi'' * dx_1_o
!
!  These coefficients are also very useful when adapting grid dependend stuff
!  such as the timestep. A simple substitution
!    1./dx_o -> dx_1_o
!  should suffice in most cases.
!
!  31-jan-17/jorgen: adapted from grid-module
!
      !TODO: Do we need both grid_profile_0D and grid_profile_1D here?
      use grid, only: grid_profile,find_star,calc_bound_coeffs

      real :: x00,y00,z00
      real :: xi1lo,xi1up,g1lo,g1up
      real :: xi2lo,xi2up,g2lo,g2up
      real :: xi3lo,xi3up,g3lo,g3up
      real, dimension(3,2) :: xi_step
      real, dimension(3,3) :: dxyz_step
      real :: xi1star,xi2star,xi3star,bound_prim1,bound_prim2
!
      real, dimension(mx_ogrid) :: g1,g1der1,g1der2,xi1,xprim2_ogrid
      real, dimension(my_ogrid) :: g2,g2der1,g2der2,xi2,yprim2_ogrid
      real, dimension(mz_ogrid) :: g3,g3der1,g3der2,xi3,zprim2_ogrid
!
      real, dimension(0:2*nprocx+1) :: xi1proc,g1proc
      real, dimension(0:2*nprocy+1) :: xi2proc,g2proc
      real, dimension(0:2*nprocz+1) :: xi3proc,g3proc
!
      real :: a,b,c
      integer :: i
!
      lequidist_ogrid=(grid_func_ogrid=='linear')
!
!  Abbreviations
!
      x0 = r_cylinder
      y0 = -pi
      z0 = xyz0(3)
      Lx = r_ogrid-r_cylinder
      Ly = 2*pi
      Lz = Lxyz(3)
!
!  Set the lower boundary and the grid size.
!
      x00 = x0
      y00 = y0
      z00 = z0
!
      dx_ogrid = Lx / merge(nxgrid_ogrid, max(nxgrid_ogrid-1,1), .false.)
      dy_ogrid = Ly / merge(nygrid_ogrid, max(nygrid_ogrid-1,1), .true.)
      dz_ogrid = Lz / merge(nzgrid_ogrid, max(nzgrid_ogrid-1,1), lperi(3))
!
!  REMOVED OPTION
!  Shift the lower boundary if requested, but only for periodic directions.
!  Shift the lower boundary if requested, but only for periodic directions.
!  Contrary to the upper case (lshift_origin)
!  REMOVED OPTION
!
!  Produce index arrays xi1, xi2, and xi3:
!    xi = 0, 1, 2, ..., N-1 for non-periodic grid
!    xi = 0.5, 1.5, 2.5, ..., N-0.5 for periodic grid
!
      do i=1,mx_ogrid; xi1(i)=i-nghost-1+ipx*nx_ogrid; enddo
      do i=1,my_ogrid; xi2(i)=i-nghost-1+ipy*ny_ogrid; enddo
      do i=1,mz_ogrid; xi3(i)=i-nghost-1+ipz*nz_ogrid; enddo
!
      xi2 = xi2 + 0.5
      if (lperi(3)) xi3 = xi3 + 0.5
!
!  Produce index arrays for processor boundaries, which are needed for
!  particle migration (see redist_particles_bounds). The select cases
!  should use these arrays to set g{2,3}proc using the grid function.
!
      do i=0,nprocx
        xi1proc(2*i)  =i*nx_ogrid-1
        xi1proc(2*i+1)=i*nx_ogrid
      enddo
      do i=0,nprocy
        xi2proc(2*i)  =i*ny_ogrid-1
        xi2proc(2*i+1)=i*ny_ogrid
      enddo
      do i=0,nprocz
        xi3proc(2*i)  =i*nz_ogrid-1
        xi3proc(2*i+1)=i*nz_ogrid
      enddo
      xi2proc = xi2proc + 0.5
      if (lperi(3)) xi3proc = xi3proc + 0.5
!
!  The following is correct for periodic and non-periodic case
!    Periodic: x(xi=0) = x0 and x(xi=N) = x1
!    Non-periodic: x(xi=0) = x0 and x(xi=N-1) = x1
!
      xi1lo=0.; xi1up=nxgrid_ogrid-merge(0.,1.,.false.)
      xi2lo=0.; xi2up=nygrid_ogrid-merge(0.,1.,.true.)
      xi3lo=0.; xi3up=nzgrid_ogrid-merge(0.,1.,lperi(3))
!
!  Construct nonequidistant grid
!
!  x coordinate
!
      if (nxgrid_ogrid==1) then
        x_ogrid = x00 + 0.5 * dx_ogrid
        ! hopefully, we will only ever multiply by the following quantities:
        xprim_ogrid = 0.
        xprim2_ogrid = 0.
        dx_1_ogrid = 0.
        dx_tilde_ogrid = 0.
        g1proc=x00
      else
!
        select case (grid_func_ogrid(1))
!
        case ('linear','sinh')
          a=coeff_grid_o(1)*dx_ogrid
          xi1star=find_star(a*xi1lo,a*xi1up,x00,x00+Lx,xyz_star_ogrid(1),grid_func_ogrid(1))/a
          call grid_profile(a*(xi1  -xi1star),grid_func_ogrid(1),g1,g1der1,g1der2)
          call grid_profile(a*(xi1lo-xi1star),grid_func_ogrid(1),g1lo)
          call grid_profile(a*(xi1up-xi1star),grid_func_ogrid(1),g1up)
!
          x_ogrid     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim_ogrid =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2_ogrid=    Lx*(g1der2*a**2)/(g1up-g1lo)
!
          ! Since lsolid_cells=True
          call grid_profile(a*(xi1proc-xi1star),grid_func_ogrid(1),g1proc)
          g1proc=x00+Lx*(g1proc  -  g1lo)/(g1up-g1lo)
!
        case default
          call fatal_error('construct_grid', &
                           'No such x grid function - '//grid_func_ogrid(1))
        endselect
!
        dx2_ogrid=xprim_ogrid**2
        dx_1_ogrid=1./xprim_ogrid
        dx_tilde_ogrid=-xprim2_ogrid/dx2_ogrid
!
        if (lfirst_proc_x) &
          dx2_bound_ogrid(-1:-nghost:-1)= 2.*(x_ogrid(l1_ogrid+1:l1_ogrid+nghost)-x_ogrid(l1_ogrid))
        if (llast_proc_x) &
          dx2_bound_ogrid(nghost:1:-1)  = 2.*(x_ogrid(l2_ogrid)-x_ogrid(l2-nghost:l2_ogrid-1))
!
        call calc_bound_coeffs(x_ogrid,coeffs_1_x_ogrid)

      endif
!
!  y coordinate
!
      if (nygrid_ogrid==1) then
        y_ogrid = y00 + 0.5 * dy_ogrid
        ! hopefully, we will only ever multiply by the following quantities:
        yprim_ogrid = 0.
        yprim2_ogrid = 0.
        dy_1_ogrid = 0.
        dy_tilde_ogrid = 0.
        g2proc=y00
      else
!
        select case (grid_func_ogrid(2))
!
        case ('linear')
!
          a=coeff_grid_o(2)*dy_ogrid
          xi2star=find_star(a*xi2lo,a*xi2up,y00,y00+Ly,xyz_star_ogrid(2),grid_func_ogrid(2))/a
          call grid_profile(a*(xi2  -xi2star),grid_func_ogrid(2),g2,g2der1,g2der2)
          call grid_profile(a*(xi2lo-xi2star),grid_func_ogrid(2),g2lo)
          call grid_profile(a*(xi2up-xi2star),grid_func_ogrid(2),g2up)
!
          y_ogrid     =y00+Ly*(g2  -  g2lo)/(g2up-g2lo)
          yprim_ogrid =    Ly*(g2der1*a   )/(g2up-g2lo)
          yprim2_ogrid=    Ly*(g2der2*a**2)/(g2up-g2lo)
!
          ! Since lsolid_cells=True
            call grid_profile(a*(xi2proc-xi2star),grid_func_ogrid(2),g2proc)
            g2proc=y00+Ly*(g2proc  -  g2lo)/(g2up-g2lo)
!
        case default
          call fatal_error('construct_grid', &
                           'No such y grid function - '//grid_func_ogrid(2))
        endselect
!
! Added parts for spherical coordinates and cylindrical coordinates.      ! JORGEN:Removed spherical part
! From now on dy = d\theta but dy_1 = 1/rd\theta and similarly for \phi.
! corresponding r and rsin\theta factors for equ.f90 (where CFL timesteps
! are estimated) are removed.
!
        dy2_ogrid=yprim_ogrid**2
        dy_1_ogrid=1./yprim_ogrid
        dy_tilde_ogrid=-yprim2_ogrid/dy2_ogrid
!
        if (lfirst_proc_y) &
          dy2_bound_ogrid(-1:-nghost:-1)= 2.*(y_ogrid(m1_ogrid+1:m1_ogrid+nghost)-y_ogrid(m1_ogrid))
        if (llast_proc_y) &
          dy2_bound_ogrid(nghost:1:-1)  = 2.*(y_ogrid(m2_ogrid)-y_ogrid(m2_ogrid-nghost:m2_ogrid-1))
!
        call calc_bound_coeffs(y_ogrid,coeffs_1_y_ogrid)

      endif
!
!  z coordinate
!
      if (nzgrid_ogrid==1) then
        z_ogrid = z00 + 0.5 * dz_ogrid
        ! hopefully, we will only ever multiply by the following quantities:
        zprim_ogrid = 0.
        zprim2_ogrid = 0.
        dz_1_ogrid = 0.
        dz_tilde_ogrid = 0.
        g3proc=z00
      else
!
        select case (grid_func_ogrid(3))
!
        case ('linear','sinh')
!
          a=coeff_grid_o(3)*dz_ogrid
          xi3star=find_star(a*xi3lo,a*xi3up,z00,z00+Lz,xyz_star_ogrid(3),grid_func_ogrid(3))/a
          call grid_profile(a*(xi3  -xi3star),grid_func_ogrid(3),g3,g3der1,g3der2)
          call grid_profile(a*(xi3lo-xi3star),grid_func_ogrid(3),g3lo)
          call grid_profile(a*(xi3up-xi3star),grid_func_ogrid(3),g3up)
!
          z_ogrid     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
          zprim_ogrid =    Lz*(g3der1*a   )/(g3up-g3lo)
          zprim2_ogrid=    Lz*(g3der2*a**2)/(g3up-g3lo)
!
          ! Since lsolid_cells is True
            call grid_profile(a*(xi3proc-xi3star),grid_func_ogrid(3),g3proc)
            g3proc=z00+Lz*(g3proc-g3lo)/(g3up-g3lo)
!
        case default
          call fatal_error('construct_grid', &
                           'No such z grid function - '//grid_func_ogrid(3))
        endselect
!
        dz2_ogrid=zprim_ogrid**2
        dz_1_ogrid=1./zprim_ogrid
        dz_tilde_ogrid=-zprim2_ogrid/dz2_ogrid
!
        if (lfirst_proc_z) &
          dz2_bound_ogrid(-1:-nghost:-1)= 2.*(z_ogrid(n1_ogrid+1:n1_ogrid+nghost)-z_ogrid(n1_ogrid))
        if (llast_proc_z) &
          dz2_bound_ogrid(nghost:1:-1)  = 2.*(z_ogrid(n2_ogrid)-z_ogrid(n2_ogrid-nghost:n2_ogrid-1))
!
        call calc_bound_coeffs(z_ogrid,coeffs_1_z_ogrid)

      endif
!
!  Compute averages across processor boundaries to calculate the physical
!  boundaries
!
      do i=0,nprocx
        procx_bounds_ogrid(i)=(g1proc(2*i)+g1proc(2*i+1))*0.5
      enddo
      do i=0,nprocy
        procy_bounds_ogrid(i)=(g2proc(2*i)+g2proc(2*i+1))*0.5
      enddo
      do i=0,nprocz
        procz_bounds_ogrid(i)=(g3proc(2*i)+g3proc(2*i+1))*0.5
      enddo
!
    endsubroutine construct_grid_ogrid
!***********************************************************************
    subroutine initialize_grid_ogrid
!
!  Coordinate-related issues: nonuniform meshes, different coordinate systems
!
!  31-jan-17/Jorgen: adapted from initialize_grid subroutine in grid.f90
!
      use Sub, only: remove_zprof
      use Mpicomm
      use IO, only: lcollective_IO
!
      real :: fact, dxmin_x, dxmin_y, dxmin_z, dxmax_x, dxmax_y, dxmax_z
      integer :: xj,yj,zj,itheta
      real :: Area_xy_ogrid, Area_yz_ogrid, Area_xz_ogrid
!
! CALL TO COORDS_AUX REMOVED, ONLY CYLINDRICAL NEEDED
!
      rcyl_mn_ogrid=x(l1_ogrid:l2_ogrid)
      if (x(l1_ogrid)==0.) then
        rcyl_mn1_ogrid(2:)=1./x(l1_ogrid+1:l2_ogrid)
        rcyl_mn1_ogrid(1)=0.
      else
        rcyl_mn1_ogrid=1./x(l1_ogrid:l2_ogrid)
      endif
      rcyl_mn2_ogrid=rcyl_mn1_ogrid**2
!
!  determine global minimum and maximum of grid spacing in any direction
!
      if (lequidist_ogrid(1) .or. nxgrid <= 1) then
        dxmin_x = dx_ogrid
        dxmax_x = dx_ogrid
      else
        dxmin_x = minval(xprim_ogrid(l1_ogrid:l2_ogrid))
        dxmax_x = maxval(xprim_ogrid(l1_ogrid:l2_ogrid))
      endif
!
      dxmin_y = dy_ogrid*minval(x_ogrid(l1_ogrid:l2_ogrid))
      dxmax_y = dy_ogrid*maxval(x_ogrid(l1_ogrid:l2_ogrid))
!
      if (lequidist_ogrid(3) .or. nzgrid <= 1) then
        dxmin_z = dz
        dxmax_z = dz
      else
        dxmin_z = minval(zprim_ogrid(n1_ogrid:n2_ogrid))
        dxmax_z = maxval(zprim_ogrid(n1_ogrid:n2_ogrid))
      endif
!
!  Find minimum/maximum grid spacing. Note that
!  minval( (/dxmin_x,dxmin_y,dxmin_z/), MASK=((/nxgrid,nygrid,nzgrid/) > 1) )
!  will be undefined if all n[xyz]grid==1, so we have to add the fourth
!  component with a test that is always true
!
      dxmin_ogrid = minval( (/dxmin_x, dxmin_y, dxmin_z, huge(dx_ogrid)/), &
                MASK=((/nxgrid_ogrid, nygrid_ogrid, nzgrid_ogrid, 2/) > 1) )

      call mpiallreduce_min(dxmin_ogrid,dxmin_x)
      dxmin_ogrid=dxmin_x
!
      if (dxmin_ogrid == 0) &
        call fatal_error ("initialize_grid", "check Lx,Ly,Lz: is one of them 0?", .true.)
!
      dxmax_ogrid = maxval( (/dxmax_x, dxmax_y, dxmax_z, epsilon(dx_ogrid)/), &
                MASK=((/nxgrid_ogrid, nygrid_ogrid, nzgrid_ogrid, 2/) > 1) )

      call mpiallreduce_max(dxmax,dxmax_x)
      dxmax_ogrid=dxmax_x
!       
!  Grid spacing. For non-equidistant grid or non-Cartesian coordinates
!  the grid spacing is calculated in the (m,n) loop.
!
!  HANDLING OF CARTESIAN REMOVED
!
! Box volume, cylinder symmetrical
!
! REMOVED CARTESIAN AND SPHERICAL
!
      box_volume_ogrid=1.
      if (nxgrid_ogrid/=1) then
          box_volume_ogrid = box_volume_ogrid*.5*(r_ogrid**2-r_cylinder**2)
      endif
      box_volume_ogrid = box_volume_ogrid*2.*pi
      if (nzgrid_ogrid/=1) box_volume_ogrid = box_volume_ogrid*Lxyz(3)
!
!  Volume element and area of coordinate surfaces.
!  Note that in the area a factor depending only on the coordinate x_i which defines the surface by x_i=const. is dropped.
!
      Area_xy_ogrid=1.; Area_yz_ogrid=1.; Area_xz_ogrid=1.
!
! CARTESIAN AND SPHERICAL REMOVED
!
!  Volume element.
!
      if (nxgrid_ogrid/=1) then
        dVol_x_ogrid=x_ogrid*xprim_ogrid
        Area_xy_ogrid=Area_xy_ogrid*1./2.*(r_ogrid**2-r_cylinder**2)
        Area_xz_ogrid=Area_xz_ogrid*(r_ogrid**2-r_cylinder**2)
      else
        dVol_x_ogrid=1./2.*(r_ogrid**2-r_cylinder**2)
        Area_xy_ogrid=Area_xy_ogrid*dVol_x_ogrid(1)
        Area_xz_ogrid=Area_xz_ogrid*Lxyz(1)
      endif
!
!  theta extent (cylindrically symmetric)
!
      if (nygrid/=1) then
        dVol_y_ogrid=yprim_ogrid
      else
        dVol_y_ogrid=2.*pi
      endif
      Area_xy_ogrid=Area_xy_ogrid*2*pi
      Area_yz_ogrid=Area_yz_ogrid*2*pi
!
!  z extent (vertically extended)
!
      if (nzgrid/=1) then
        dVol_z_ogrid=zprim_ogrid
        Area_xz_ogrid=Area_xz_ogrid*Lxyz(3)
        Area_yz_ogrid=Area_yz_ogrid*Lxyz(3)
      else
        dVol_z_ogrid=1.
      endif
!
!  Trapezoidal rule
!
      rcyl_weight_ogrid=rcyl_mn_ogrid
      if (lfirst_proc_x) rcyl_weight_ogrid( 1)=.5*rcyl_weight_ogrid( 1)
      if (llast_proc_x ) rcyl_weight_ogrid(nx_ogrid)=.5*rcyl_weight_ogrid(nx_ogrid)
!
!  Stop if no existing coordinate system is specified
!
!  Inverse volume elements
!
      dVol1_x_ogrid = 1./dVol_x_ogrid
      dVol1_y_ogrid = 1./dVol_y_ogrid
      dVol1_z_ogrid = 1./dVol_z_ogrid
!
!  Define inner and outer radii for non-cartesian coords.
!  If the user did not specify them yet (in start.in),
!  these are the first point of the first x-processor,
!  and the last point of the last x-processor.
!
      if (nprocx/=1) then
!
!  The root (iproc=0) has by default the first value of x
!
        if (lroot) then
          if (r_int_ogrid==0) r_int_ogrid=x_ogrid(l1_ogrid)
!
!  The root should also receive the value of r_ext from
!  from the last x-processor (which is simply nprocx-1
!  for iprocy=0 and iprocz=0) for broadcasting.
!
          if (r_ext_ogrid==impossible) &
            call mpirecv_real(r_ext_ogrid,nprocx-1,111)
          endif
!
!  The last x-processor knows the value of r_ext, and sends
!  it to root, for broadcasting.
!
          if ((r_ext_ogrid==impossible).and.&
              (llast_proc_x.and.lfirst_proc_y.and.lfirst_proc_z)) then
            r_ext_ogrid=x_ogrid(l2_ogrid)
            call mpisend_real(r_ext_ogrid,0,111)
          endif
!
!  Broadcast the values of r_int and r_ext
!
        call mpibcast_real(r_int_ogrid,comm=MPI_COMM_WORLD)
        call mpibcast_real(r_ext_ogrid,comm=MPI_COMM_WORLD)
      else
!
!  Serial-x. Just get the local grid values.
!
        if (r_int_ogrid == 0)         r_int_ogrid=x_ogrid(l1_ogrid)
        if (r_ext_ogrid ==impossible) r_ext_ogrid=x_ogrid(l2_ogrid)
      endif
      if (lroot) print*,'initialize_grid, r_int,r_ext=',r_int_ogrid,r_ext_ogrid
!
!  For a non-periodic mesh, multiply boundary points by 1/2.
!  Do it for each direction in turn.
!  If a direction has no extent, it is automatically periodic
!  and the corresponding step is therefore not called.
!
      if (lfirst_proc_x) dVol_x_ogrid(l1_ogrid)=.5*dVol_x_ogrid(l1_ogrid)
      if (llast_proc_x ) dVol_x_ogrid(l2_ogrid)=.5*dVol_x_ogrid(l2_ogrid)
!
      if (.not.lperi(3)) then
        if (lfirst_proc_z) dVol_z_ogrid(n1_ogrid)=.5*dVol_z_ogrid(n1_ogrid)
        if (llast_proc_z ) dVol_z_ogrid(n2_ogrid)=.5*dVol_z_ogrid(n2_ogrid)
      endif
!
!  Print the value for which output is being produced.
!  (Have so far only bothered about single processor output.)
!
      if (lroot) then
        lpoint_ogrid=min(max(l1_ogrid,lpoint_ogrid),l2_ogrid)
        mpoint_ogrid=min(max(m1_ogrid,mpoint_ogrid),m2_ogrid)
        npoint_ogrid=min(max(n1_ogrid,npoint_ogrid),n2_ogrid)
        lpoint2_ogrid=min(max(l1_ogrid,lpoint2_ogrid),l2_ogrid)
        mpoint2_ogrid=min(max(m1_ogrid,mpoint2_ogrid),m2_ogrid)
        npoint2_ogrid=min(max(n1_ogrid,npoint2_ogrid),n2_ogrid)
        print*,'(x,y,z)_ogrid(point)=',x_ogrid(lpoint_ogrid),y_ogrid(mpoint_ogrid),z_ogrid(npoint_ogrid)
        print*,'(x,y,z)_ogrid(point2)=',x_ogrid(lpoint2_ogrid),y_ogrid(mpoint2_ogrid),z_ogrid(npoint2_ogrid)
      endif
!
!  Clean up profile files.
!
!TODO: DO I NEED TO DO ANYTHING WITH THIS?
      if (lroot.or..not.lcollective_IO) call remove_zprof
      lwrite_prof=.true.
!
!  Set the the serial grid arrays, that contain the coordinate values
!  from all processors.
!
      call construct_serial_arrays
!
    endsubroutine initialize_grid_ogrid
!!***********************************************************************
! TODO: NOT A PRIORITY AT THE MOMENT
!
!    subroutine save_grid(lrestore)
!!
!!  Saves grid into local statics (needed for downsampled output)
!!
!!  6-mar-14/MR: coded
!!
!      use General, only : loptest
!!
!      logical, optional :: lrestore
!!
!      real, dimension(mx), save :: xs,dx_1s,dx_tildes
!      real, dimension(my), save :: ys,dy_1s,dy_tildes
!      real, dimension(mz), save :: zs,dz_1s,dz_tildes
!      real, dimension(nx), save :: r_mns,r1_mns,r2_mns
!      real, dimension(my), save :: sinths,sin1ths,sin2ths,cosths, &
!                                   cotths,cos1ths,tanths
!      real, dimension(mz), save :: sinphs, cosphs
!      real, dimension(nx), save :: rcyl_mns,rcyl_mn1s,rcyl_mn2s
!
!      real, save :: dxs, dys, dzs
!
!      logical, save :: lfirst=.true.
!
!      if (loptest(lrestore)) then
!        if (lfirst) then
!          call fatal_error('save_grid','first call must have lrestore=F')
!        else
!          dx=dxs; dy=dys; dz=dzs
!          x=xs; y=ys; z=zs
!          dx_1_ogrid=dx_1s; dy_1=dy_1s; dz_1=dz_1s
!          dx_tilde_ogrid=dx_tildes; dy_tilde=dy_tildes; dz_tilde=dz_tildes
!
!          if (lspherical_coords) then
!            r_mn=r_mns; r1_mn=r1_mns; r2_mn=r2_mns
!            sinth=sinths; sin1th=sin1ths; sin2th=sin2ths; costh=cosths
!            cotth=cotths; cos1th=cos1ths; tanth=tanths
!            sinphs=sinph; cosphs=cosph
!          elseif (lcylindrical_coords) then
!            rcyl_mn=rcyl_mns; rcyl_mn1=rcyl_mn1s; rcyl_mn2=rcyl_mn2s
!          endif
!        endif
!      elseif (lfirst) then
!        lfirst=.false.
!        dxs=dx; dys=dy; dzs=dz
!        xs=x; ys=y; zs=z
!        dx_1s=dx_1; dy_1s=dy_1; dz_1s=dz_1
!        dx_tildes=dx_tilde; dy_tildes=dy_tilde; dz_tildes=dz_tilde
!
!        if (lspherical_coords) then
!          r_mns=r_mn; r1_mns=r1_mn; r2_mns=r2_mn
!          sinths=sinth; sin1ths=sin1th; sin2ths=sin2th; cosths=costh
!          cotths=cotth; cos1ths=cos1th; tanths=tanth
!          sinph=sinphs; cosph=cosphs
!        elseif (lcylindrical_coords) then
!          rcyl_mns=rcyl_mn; rcyl_mn1s=rcyl_mn1; rcyl_mn2s=rcyl_mn2
!        endif
!      endif
!
!    endsubroutine save_grid
!!***********************************************************************
    subroutine pencil_criteria_grid
!
!  All pencils that this special module depends on are specified here.
!
!  01-feb-17/Jorgen: Adapted from grid.f90
!
      if (any(lfreeze_varext).or.any(lfreeze_varint)) then
        lpenc_requested(i_rcyl_mn)=.true.
      endif
!
    endsubroutine pencil_criteria_grid
!***********************************************************************
    subroutine pencil_interdep_grid(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  01-feb-17/Jorgen: Adapted from grid.f90
!
! TODO: where are the index variables i_***  initialized and set?
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_rcyl_mn1)) lpencil_in(i_rcyl_mn)=.true.
      if (lpencil_in(i_r_mn1)) lpencil_in(i_r_mn)=.true.
      if (lpencil_in(i_evr)) lpencil_in(i_r_mn)=.true.
      if (lpencil_in(i_evth).or.lpencil_in(i_evr)) then
         lpencil_in(i_pomx)=.true.
         lpencil_in(i_pomy)=.true.
         lpencil_in(i_rcyl_mn)=.true.
         lpencil_in(i_r_mn1)=.true.
      endif
!
      if (lpencil_in(i_rr)) then
        lpencil_in(i_x_mn)=.true.
        lpencil_in(i_y_mn)=.true.
        lpencil_in(i_z_mn)=.true.
      endif
!
    endsubroutine pencil_interdep_grid
!***********************************************************************
!TODO: Do I need to change p%rcyl_mn to p%rcyl_mn_ogrid, etc.
    subroutine calc_pencils_grid_ogrid(f,p,m_ogrid,n_ogrid)
!
!  Calculate Grid/geometry related pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   31-jan-17/Jorgen: Adapted from calc_pencils_grid in grid.f90
!                     Only cylindrical coodrinats included
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      type (pencil_case) :: p
      integer, intent(in) :: m_ogrid,n_ogrid
!
      intent(in) :: f
      intent(inout) :: p
!
      if (lpencil(i_x_mn))     p%x_mn    = x_ogrid(l1_ogrid:l2_ogrid)*cos(x_ogrid(m_ogrid))
      if (lpencil(i_y_mn))     p%y_mn    = x_ogrid(l1_ogrid:l2_ogrid)*sin(y_ogrid(m_ogrid))
      if (lpencil(i_z_mn))     p%z_mn    = spread(z_ogrid(n_ogrid),1,nx_ogrid)
      if (lpencil(i_r_mn))     p%r_mn    = sqrt(x_ogrid(l1_ogrid:l2_ogrid)**2+z_ogrid(n_ogrid)**2)
      if (lpencil(i_rcyl_mn))  p%rcyl_mn = x_ogrid(l1_ogrid:l2_ogrid)
      if (lpencil(i_phi_mn))   p%phi_mn  = spread(y_ogrid(m_ogrid),1,nx_ogrid)
      if (lpencil(i_rcyl_mn1)) p%rcyl_mn1=1./max(p%rcyl_mn,tini)
      if (lpencil(i_r_mn1))    p%r_mn1   =1./max(p%r_mn,tini)
      if (lpencil(i_pomx))     p%pomx    = 1.
      if (lpencil(i_pomy))     p%pomy    = 0.
      if (lpencil(i_phix))     p%phix    = 0.
      if (lpencil(i_phiy))     p%phiy    = 1.
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_grid_ogrid
!***********************************************************************
! TODO: COULD IT BE A PROBLEM HERE, CALLING A VARIABLE X WHEN X IS DEFINED IN CDATA.F90?
    subroutine real_to_index(n, x, xi)
!
!  Transforms coordinates in real space to those in index space.
!
!  10-sep-15/ccyang: coded.
!
      integer, intent(in) :: n
      real, dimension(n,3), intent(in) :: x
      real, dimension(n,3), intent(out) :: xi
!
      real, parameter :: ngp1 = nghost + 1
      integer :: i
!
!  Work on each direction.
!
      nonzero: if (n > 0) then
        dir: do i = 1, 3
          if (lactive_dimension(i)) then
            call inverse_grid(i, x(:,i), xi(:,i), local=.true.)
          else
            xi(:,i) = ngp1
          endif
        enddo dir
      endif nonzero
!
    endsubroutine real_to_index
!***********************************************************************
    subroutine inverse_grid(dir, x, xi, local)
!
!  Transform the x coordinates in real space to the xi coordinates in
!  index space in dir direction, where dir = 1, 2, or, 3.
!  If local is present and .true., the index space is with respect to
!  the local grid.
!
!  31-jan-17/Jorgen: adapted from inverse_grid in grid.f90
!
      use General, only: arcsinh
!
      integer, intent(in) :: dir
      real, dimension(:), intent(in) :: x
      real, dimension(:), intent(out) :: xi
      logical, intent(in), optional :: local
!
      character(len=linelen) :: msg
      logical :: loc
      integer :: shift
      real :: h, a, b, c
      ! TODO: Move these to top?
      real, dimension(3) :: xyz0_ogrid,Lxyz_ogrid
      ! Shorthand notation
      xyz0_ogrid(1)=r_cylinder
      xyz0_ogrid(2)=-pi
      xyz0_ogrid(3)=xyz0(3)
      Lxyz_ogrid(1)=r_ogrid-r_cylinder
      Lxyz_ogrid(2)=2*pi
      Lxyz_ogrid(2)=Lxyz(3)
      
!
!  Sanity check.
!
      if (any(lshift_origin) .or. any(lshift_origin_lower)) &
          call fatal_error('inverse_grid', 'lshift_origin and lshift_origin_lower are not supported. ')
!
!  Global or local index space?
!
      loc = .false.
      if (present(local)) loc = local
!
!  Check the direction.
!
      ckdir: select case (dir)
      case (1) ckdir
        h = dx_ogrid
        if (loc) shift = nx_ogrid * ipx
      case (2) ckdir
        h = dy_ogrid
        if (loc) shift = ny_ogrid * ipy
      case (3) ckdir
        h = dz_ogrid
        if (loc) shift = nz_ogrid * ipz
      case default ckdir
        write(msg,*) 'unknown direction dir = ', dir
        call fatal_error('inverse_grid', trim(msg))
      endselect ckdir
!
!  Make the inversion according to the grid function.
!
      func: select case (grid_func_ogrid(dir))
!
      case ('linear') func
        xi = (x - xyz0_ogrid(dir)) / h
!
      case ('sinh') func
        a = coeff_grid_o(dir) * Lxyz_ogrid(dir)
        b = sinh(a)
        c = cosh(a) - 1.0
        a = (xyz_star_ogrid(dir) - xyz0_ogrid(dir)) / Lxyz_ogrid(dir)
        a = a * b / sqrt(1.0 + 2.0 * a * (1.0 - a) * c)
        b = (sqrt(1.0 + a * a) * b - a * c) / Lxyz_ogrid(dir)
        xi = (arcsinh(a) + arcsinh(b * (x - xyz0_ogrid(dir)) - a)) / (coeff_grid_o(dir) * h)
!
      case default func
        call fatal_error('inverse_grid', 'unknown grid function ' // trim(grid_func_ogrid(dir)))
!
      endselect func
!
!  Shift to match the global index space.
!
      if (dir==2 .or. (dir==3 .and. lperi(dir))) then
        xi = xi + real(nghost) + 0.5
      else
        xi = xi + real(nghost + 1)
      endif
!
!  Convert to the local index space if requested.
!
      if (loc) xi = xi - real(shift)
!
    endsubroutine inverse_grid
!***********************************************************************
    subroutine construct_serial_arrays
!
!  The arrays xyz are local only, yet sometimes the serial array is
!  needed. Construct here the serial arrays out of the local ones,
!  but gathering them processor-wise and broadcasting the constructed
!  array. This is only done in start time, so legibility (3 near-copies
!  of the same code) is preferred over code-reusability (one general
!  piece of code called three times).
!
!  31-jan-17/Jorgen: Adapted from construct_serial_arrays in grid.f90
!
      use Mpicomm, only: mpisend_real,mpirecv_real,mpibcast_real, mpiallreduce_sum_int, MPI_COMM_WORLD
!
      real, dimension(nx_ogrid) :: xrecv, x1recv, x2recv
      real, dimension(ny_ogrid) :: yrecv, y1recv, y2recv
      real, dimension(nz_ogrid) :: zrecv, z1recv, z2recv
      integer :: jx,jy,jz,iup,ido,iproc_recv
      integer :: iproc_first, iproc_last
!
      xrecv=0.; yrecv=0.; zrecv=0.
      x1recv=0.; y1recv=0.; z1recv=0.
      x2recv=0.; y2recv=0.; z2recv=0.
!
!  Serial x array
!
      if (iproc/=root) then
!
!  All processors of the same row (ipx,ipy or ipz)
!  send their array values to the root.
!
        if ((ipy==0).and.(ipz==0)) then
          call mpisend_real(x_ogrid(l1_ogrid:l2_ogrid),nx_ogrid,root,111)
          call mpisend_real(dx_1_ogrid(l1_ogrid:l2_ogrid),nx_ogrid,root,112)
          call mpisend_real(dx_tilde_ogrid(l1_ogrid:l2_ogrid),nx_ogrid,root,113)
        endif
      else
!
!  The root processor, in turn, receives the data from the others
!
        do jx=0,nprocx-1
          !avoid send-to-self
          if (jx/=root) then
!
!  Formula of the serial processor number:
!  iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
!  Since for the x-row ipy=ipz=0, this reduces
!  to iproc_recv=jx.
!
            iproc_recv=jx
            call mpirecv_real(xrecv,nx_ogrid,iproc_recv,111)
            call mpirecv_real(x1recv,nx_ogrid,iproc_recv,112)
            call mpirecv_real(x2recv,nx_ogrid,iproc_recv,113)
!
            ido=jx    *nx_ogrid + 1
            iup=(jx+1)*nx_ogrid
            xgrid_ogrid(ido:iup)=xrecv
            dx1grid_ogrid(ido:iup)=x1recv
            dxtgrid_ogrid(ido:iup)=x2recv
          else
            !the root just copies its value to the serial array
            xgrid_ogrid(1:nx_ogrid)=x_ogrid(l1_ogrid:l2_ogrid)
            dx1grid_ogrid(1:nx_ogrid)=dx_1_ogrid(l1_ogrid:l2_ogrid)
            dxtgrid_ogrid(1:nx_ogrid)=dx_tilde_ogrid(l1_ogrid:l2_ogrid)
          endif
        enddo
      endif
!
!  Serial array constructed. Broadcast the result. Repeat the
!  procedure for y and z arrays.
!
      call mpibcast_real(xgrid_ogrid,nxgrid_ogrid,comm=MPI_COMM_WORLD)
      call mpibcast_real(dx1grid_ogrid,nxgrid_ogrid,comm=MPI_COMM_WORLD)
      call mpibcast_real(dxtgrid_ogrid,nxgrid_ogrid,comm=MPI_COMM_WORLD)
!
!  Serial y-array
!
      if (iproc/=root) then
        if (ipx==0.and.ipz==0) then
          call mpisend_real(y_ogrid(m1_ogrid:m2_ogrid),ny_ogrid,root,221)
          call mpisend_real(dy_1_ogrid(m1_ogrid:m2_ogrid),ny_ogrid,root,222)
          call mpisend_real(dy_tilde_ogrid(m1_ogrid:m2_ogrid),ny_ogrid,root,223)
        endif
      else
        do jy=0,nprocy-1
          if (jy/=root) then
            iproc_recv=nprocx*jy
            call mpirecv_real(yrecv,ny_ogrid,iproc_recv,221)
            call mpirecv_real(y1recv,ny_ogrid,iproc_recv,222)
            call mpirecv_real(y2recv,ny_ogrid,iproc_recv,223)
            ido=jy    *ny_ogrid + 1
            iup=(jy+1)*ny_ogrid
            ygrid_ogrid(ido:iup)=yrecv
            dy1grid_ogrid(ido:iup)=y1recv
            dytgrid_ogrid(ido:iup)=y2recv
          else
            ygrid_ogrid(1:ny_ogrid)=y_ogrid(m1_ogrid:m2_ogrid)
            dy1grid_ogrid(1:ny_ogrid)=dy_1_ogrid(m1_ogrid:m2_ogrid)
            dytgrid_ogrid(1:ny_ogrid)=dy_tilde_ogrid(m1_ogrid:m2_ogrid)
          endif
        enddo
      endif
      call mpibcast_real(ygrid_ogrid,nygrid_ogrid)
      call mpibcast_real(dy1grid_ogrid,nygrid_ogrid)
      call mpibcast_real(dytgrid_ogrid,nygrid_ogrid)
!
!  Serial z-array
!
      if (iproc/=root) then
        if (ipx==0.and.ipy==0) then
          call mpisend_real(z_ogrid(n1_ogrid:n2_ogrid),nz_ogrid,root,331)
          call mpisend_real(dz_1_ogrid(n1_ogrid:n2_ogrid),nz_ogrid,root,332)
          call mpisend_real(dz_tilde_ogrid(n1_ogrid:n2_ogrid),nz_ogrid,root,333)
        endif
      else
        do jz=0,nprocz-1
          if (jz/=root) then
            iproc_recv=nprocx*nprocy*jz
            call mpirecv_real(zrecv,nz_ogrid,iproc_recv,331)
            call mpirecv_real(z1recv,nz_ogrid,iproc_recv,332)
            call mpirecv_real(z2recv,nz_ogrid,iproc_recv,333)
            ido=jz    *nz_ogrid + 1
            iup=(jz+1)*nz_ogrid
            zgrid_ogrid(ido:iup)=zrecv
            dz1grid_ogrid(ido:iup)=z1recv
            dztgrid_ogrid(ido:iup)=z2recv
          else
            zgrid_ogrid(1:nz_ogrid)=z_ogrid(n1_ogrid:n2_ogrid)
            dz1grid_ogrid(1:nz_ogrid)=dz_1_ogrid(n1_ogrid:n2_ogrid)
            dztgrid_ogrid(1:nz_ogrid)=dz_tilde_ogrid(n1_ogrid:n2_ogrid)
          endif
        enddo
      endif
      call mpibcast_real(zgrid_ogrid,nzgrid_ogrid)
      call mpibcast_real(dz1grid_ogrid,nzgrid_ogrid)
      call mpibcast_real(dztgrid_ogrid,nzgrid_ogrid)
!
!  Check the first and last processors.
!
      iup = 0
      if (lfirst_proc_xyz) iup = iproc
      ido = 0
      if (llast_proc_xyz) ido = iproc
      call mpiallreduce_sum_int(iup, iproc_first)
      call mpiallreduce_sum_int(ido, iproc_last)
!
!  Communicate the ghost cells.
!
      xglobal_ogrid(nghost+1:mxgrid_ogrid-nghost) = xgrid_ogrid
      yglobal_ogrid(nghost+1:mygrid_ogrid-nghost) = ygrid_ogrid
      zglobal_ogrid(nghost+1:mzgrid_ogrid-nghost) = zgrid_ogrid
!
      xglobal_ogrid(1:nghost) = x_ogrid(1:nghost)
      yglobal_ogrid(1:nghost) = y_ogrid(1:nghost)
      zglobal_ogrid(1:nghost) = z_ogrid(1:nghost)
!
      xglobal_ogrid(mxgrid_ogrid-nghost+1:mxgrid_ogrid) = x_ogrid(mx_ogrid-nghost+1:mx_ogrid)
      yglobal_ogrid(mygrid_ogrid-nghost+1:mygrid_ogrid) = y_ogrid(my_ogrid-nghost+1:my_ogrid)
      zglobal_ogrid(mzgrid_ogrid-nghost+1:mzgrid_ogrid) = z_ogrid(mz_ogrid-nghost+1:mz_ogrid)
!
      call mpibcast_real(xglobal_ogrid(1:nghost), nghost, iproc_first)
      call mpibcast_real(yglobal_ogrid(1:nghost), nghost, iproc_first)
      call mpibcast_real(zglobal_ogrid(1:nghost), nghost, iproc_first)
!
      call mpibcast_real(xglobal_ogrid(mxgrid_ogrid-nghost+1:mxgrid_ogrid), nghost, iproc_last)
      call mpibcast_real(yglobal_ogrid(mygrid_ogrid-nghost+1:mygrid_ogrid), nghost, iproc_last)
      call mpibcast_real(zglobal_ogrid(mzgrid_ogrid-nghost+1:mzgrid_ogrid), nghost, iproc_last)
!
    endsubroutine construct_serial_arrays
!***********************************************************************
    subroutine get_grid_mn_ogrid(m_ogrid,n_ogrid)
!
    integer, intent(in) :: m_ogrid,n_ogrid
!  Gets the geometry of the pencil at each (m,n) in the mn-loop.
!
!  31-jan-17/Jorgen Adapted from get_grid_mn in grid.f90

      dline_1_ogrid(:,1) = dx_1_ogrid(l1_ogrid:l2_ogrid)
      dline_1_ogrid(:,2) = rcyl_mn1_ogrid * dy_1_ogrid(m_ogrid)
      dline_1_ogrid(:,3) = dz_1_ogrid(n_ogrid)
!
      dxmax_pencil_ogrid = 0.
      dxmin_pencil = 0.
      if (nxgrid_ogrid /= 1) then 
        dxmax_pencil_ogrid =     1.0 / dline_1_ogrid(:,1)
        dxmin_pencil_ogrid =     1.0 / dline_1_ogrid(:,1)
      endif
      if (nygrid_ogrid /= 1) then 
        dxmax_pencil_ogrid = max(1.0 / dline_1_ogrid(:,2), dxmax_pencil_ogrid)
        dxmin_pencil_ogrid = min(1.0 / dline_1_ogrid(:,2), dxmin_pencil_ogrid)
      endif
      if (nzgrid_ogrid /= 1) then 
        dxmax_pencil_ogrid = max(1.0 / dline_1_ogrid(:,3), dxmax_pencil_ogrid)
        dxmin_pencil_ogrid = min(1.0 / dline_1_ogrid(:,3), dxmin_pencil_ogrid)
      endif
!
      if (lmaximal_cdtv) then
        dxyz_2_ogrid = max(dline_1_ogrid(:,1)**2, dline_1_ogrid(:,2)**2, dline_1_ogrid(:,3)**2)
        dxyz_4_ogrid = max(dline_1_ogrid(:,1)**4, dline_1_ogrid(:,2)**4, dline_1_ogrid(:,3)**4)
        dxyz_6_ogrid = max(dline_1_ogrid(:,1)**6, dline_1_ogrid(:,2)**6, dline_1_ogrid(:,3)**6)
      else
        dxyz_2_ogrid = dline_1_ogrid(:,1)**2 + dline_1_ogrid(:,2)**2 + dline_1_ogrid(:,3)**2
        dxyz_4_ogrid = dline_1_ogrid(:,1)**4 + dline_1_ogrid(:,2)**4 + dline_1_ogrid(:,3)**4
        dxyz_6_ogrid = dline_1_ogrid(:,1)**6 + dline_1_ogrid(:,2)**6 + dline_1_ogrid(:,3)**6
      endif
!
      dVol_ogrid = dVol_x_ogrid(l1_ogrid:l2_ogrid)*dVol_y_ogrid(m_ogrid)*dVol_z_ogrid(n_ogrid)
!
    endsubroutine get_grid_mn_ogrid
!***********************************************************************
    subroutine grid_bound_data

!  31-jan-17/Jorgen: adapted from version in grid.f90

      use grid, only: calc_bound_coeffs

      if (nxgrid_ogrid>1) then
        if (lfirst_proc_x) &
          dx2_bound_ogrid(-1:-nghost:-1)= 2.*(x_ogrid(l1_ogrid+1:l1_ogrid+nghost)-x(l1_ogrid))
        if (llast_proc_x) &
          dx2_bound_ogrid(nghost:1:-1)  = 2.*(x_ogrid(l2_ogrid)-x(l2-nghost:l2_ogrid-1))
!
        call calc_bound_coeffs(x_ogrid,coeffs_1_x_ogrid)
      endif

      if (nygrid_ogrid>1) then
        if (lfirst_proc_y) &
          dy2_bound_ogrid(-1:-nghost:-1)= 2.*(y_ogrid(m1_ogrid+1:m1_ogrid+nghost)-y_ogrid(m1_ogrid))
        if (llast_proc_y) &
          dy2_bound_ogrid(nghost:1:-1)  = 2.*(y_ogrid(m2)-y_ogrid(m2-nghost:m2-1))
!
        call calc_bound_coeffs(y_ogrid,coeffs_1_y_ogrid)
      endif

      if (nzgrid_ogrid>1) then
        if (lfirst_proc_z) &
          dz2_bound_ogrid(-1:-nghost:-1)= 2.*(z_ogrid(n1_ogrid+1:n1_ogrid+nghost)-z_ogrid(n1_ogrid))
        if (llast_proc_z) &
          dz2_bound_ogrid(nghost:1:-1)  = 2.*(z_ogrid(n2_ogrid)-z_ogrid(n2_ogrid-nghost:n2_ogrid-1))
!
        call calc_bound_coeffs(z_ogrid,coeffs_1_z_ogrid)
      endif

    endsubroutine grid_bound_data
!***********************************************************************
    subroutine time_step_ogrid(f_ogrid,df_ogrid,p_ogrid)
!
      use Mpicomm, only: mpifinalize, mpiallreduce_max, MPI_COMM_WORLD
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f_ogrid
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar_ogrid) :: df_ogrid
      type (pencil_case) :: p_ogrid
      real :: ds, dtsub
      real :: dt1, dt1_local, dt1_last=0.0
      integer :: j
!
!  Coefficients for up to order 3.
!
      if (itorder==1) then
        alpha_ts=(/ 0.0, 0.0, 0.0 /)
        beta_ts =(/ 1.0, 0.0, 0.0 /)
      elseif (itorder==2) then
        alpha_ts=(/   0.0, -1/2.0, 0.0 /)
        beta_ts =(/ 1/2.0,    1.0, 0.0 /)
      elseif (itorder==3) then
        !alpha_ts=(/   0.0, -2/3.0, -1.0   /)
        !beta_ts =(/ 1/3.0,    1.0,  1/2.0 /)
        !  use coefficients of Williamson (1980)
        alpha_ts=(/   0.0, -5/9.0 , -153/128.0 /)
        beta_ts =(/ 1/3.0, 15/16.0,    8/15.0  /)
      else
        if (lroot) print*,'Not implemented: itorder=',itorder
        call mpifinalize
      endif
!
!  dt_beta_ts may be needed in other modules (like Dustdensity) for fixed dt.
!
      if (.not. ldt) dt_beta_ts=dt*beta_ts
!
!  Set up df and ds for each time sub.
!
      do itsub=1,itorder
        lfirst=(itsub==1)
        llast=(itsub==itorder)
        if (lfirst) then
          df=0.0
          ds=0.0
        else
          df=alpha_ts(itsub)*df !(could be subsumed into pde, but is dangerous!)
          ds=alpha_ts(itsub)*ds
        endif
!
!  Change df according to the chosen physics modules.
!
        call pde_ogrid(f_ogrid,df_ogrid,p_ogrid)
!
        ds=ds+1.0
!
!  If we are in the first time substep we need to calculate timestep dt.
!  Only do it on the root processor, then broadcast dt to all others.
!
        if (lfirst.and.ldt) then
          dt1_local=maxval(dt1_max(1:nx))
          ! Timestep growth limiter
          if (real(ddt) > 0.) dt1_local=max(dt1_local,dt1_last)
          call mpiallreduce_max(dt1_local,dt1,MPI_COMM_WORLD)
          dt=1.0/dt1
          if (loutput_varn_at_exact_tsnap) call shift_dt(dt)
          ! Timestep growth limiter
          if (ddt/=0.) dt1_last=dt1_local/ddt
        endif
!
!  Calculate dt_beta_ts.
!
        if (ldt) dt_beta_ts=dt*beta_ts
        if (ip<=6) print*, 'time_step: iproc, dt=', iproc_world, dt 
        dtsub = ds * dt_beta_ts(itsub)
!
!  Time evolution of grid variables.
!  (do this loop in pencils, for cache efficiency)
!
        do j=1,mvar 
          do n_ogrid=n1_ogrid,n2_ogrid
            do m_ogrid=m1_ogrid,m2_ogrid
              f_ogrid(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,j)=f_ogrid(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,j) &
                +dt_beta_ts(itsub)*df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,j)
            enddo
          enddo
        enddo
!
!  Increase time.
!
        t = t + dtsub
!
      enddo
!
    endsubroutine time_step_ogrid
!***********************************************************************
    subroutine pde_ogrid(f,df,p)
!
!  Call the different evolution equations.
!
      use General, only: notanumber
      use Mpicomm
      use Sub
!
      logical :: early_finalize
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar_ogrid) :: df
      real, dimension (nx_ogrid,3) :: df_iuu_pencil
      type (pencil_case) :: p
      real, dimension(1)   :: mass_per_proc
      integer :: iv,nyz
!
      intent(inout)  :: f       ! inout due to  lshift_datacube_x,
                                ! density floor, or velocity ceiling
      intent(out)    :: df, p
!
!  Print statements when they are first executed.
!
      headtt = headt .and. lfirst .and. lroot
!
      if (headtt.or.ldebug) print*,'pde_ogrid: ENTER'
      if (headtt) call svn_id( &
           "$Id$")
!
!  Prepare x-ghost zones; required before f-array communication
!  AND shock calculation
!
      call boundconds_x_ogrid(f)
!
!  Initiate communication and do boundary conditions.
!  Required order:
!  1. x-boundaries (x-ghost zones will be communicated) - done above
!  2. communication
!  3. y- and z-boundaries
!
      if (ldebug) print*,'pde: bef. initiate_isendrcv_bdry'
      call initiate_isendrcv_bdry_ogrid(f)
      call finalize_isendrcv_bdry_ogrid(f)
      call boundconds_y_ogrid(f)
      call boundconds_z_ogrid(f)
!
!------------------------------------------------------------------------------
!  Do loop over m and n.
!
      nyz=ny_ogrid*nz_ogrid
      mn_loop: do imn=1,nyz
        n_ogrid=nn(imn)
        m_ogrid=mm(imn)
        lfirstpoint=(imn==1)      ! true for very first m-n loop
        llastpoint=(imn==(nyz)) ! true for very last m-n loop
!
!  Grid spacing. In case of equidistant grid and cartesian coordinates
!  this is calculated before the (m,n) loop.
!
        call get_grid_mn_ogrid
!
!  Calculate grid/geometry related pencils.
!
        call calc_pencils_grid_ogrid(f,p)
!
!  Calculate pencils for the pencil_case.
!
        call calc_pencils_hydro_ogrid(f,p)
        call calc_pencils_density_ogrid(f,p)
        call calc_pencils_eos_ogrid(f,p)
        call calc_pencils_viscosity_ogrid(f,p)
        call calc_pencils_energy_ogrid(f,p)
!
!  --------------------------------------------------------
!  NO CALLS MODIFYING PENCIL_CASE PENCILS BEYOND THIS POINT
!  --------------------------------------------------------
!
!  hydro, density, and entropy evolution
!  Note that pressure gradient is added in denergy_dt of noentropy to momentum,
!  even if lentropy=.false.
!
        call duu_dt_ogrid(f,df,p)
        call dlnrho_dt_ogrid(f,df,p)
        call denergy_dt_ogrid(f,df,p)
!
!  End of loops over m and n.
!
        headtt_ogrid=.false.
      enddo mn_loop
!
    endsubroutine pde_ogrid
!***********************************************************************
    subroutine duu_dt_ogrid(f,df,p)
!
!  velocity evolution
!  calculate du/dt = - u.gradu - 2Omega x u + grav + Fvisc
!  pressure gradient force added in density and entropy modules.
!
      use Diagnostics
      use Sub, only: vecout, dot, dot2, identify_bcs, cross, multsv_mn_add
      use General, only: transform_thph_yy
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mfarray_ogrid) :: f
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mvar_ogrid) :: df
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(inout) :: f,df
!
!  Advection term.
!
      df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)-p%ugu
!
! calculate viscous force
!
      call calc_viscous_force_ogrid(df,p)
!
    endsubroutine duu_dt_ogrid
!***********************************************************************
    subroutine dlnrho_dt_ogrid(f,df,p)
!
!  Continuity equation.
!  Calculate dlnrho/dt = - u.gradlnrho - divu
!
      use Diagnostics
      use Sub
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mfarray_ogrid) :: f
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mvar_ogrid) :: df
      type (pencil_case) :: p
!
      intent(in)  :: p
      intent(inout) :: df,f
!
      real, dimension (nx_ogrid) :: density_rhs 
      integer :: j
!
!  Continuity equation.
!      
      density_rhs= - p%ugrho   - p%rho*p%divu      
!
!  Add the continuity equation terms to the RHS of the density df.
!
      df(l1:l2,m,n,irho) = df(l1:l2,m,n,irho) + density_rhs
!
    endsubroutine dlnrho_dt_ogrid
!***********************************************************************
    subroutine denergy_dt_ogrid(f,df,p)
!
!  Calculate pressure gradient term for isothermal/polytropic equation
!  of state.
!
      use EquationOfState, only: beta_glnrho_global, beta_glnrho_scaled
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension(nx) :: ufpres, uduu
      integer :: j,ju
      integer :: i
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Add isothermal/polytropic pressure term in momentum equation.
!
      if (lhydro.and.lpressuregradient_gas) then
        do j=1,3
          ju=j+iuu-1
          df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)+p%fpres(:,j)
        enddo
!
    endsubroutine denergy_dt_ogrid
!***********************************************************************
    subroutine calc_pencils_eos_std_ogrid(f,p)
!
! Envelope adjusting calc_pencils_eos_pencpar to the standard use with
! lpenc_loc=lpencil
!
      real, dimension (mx,my,mz,mfarray),intent(INOUT):: f
      type (pencil_case),                intent(OUT)  :: p
!
      call calc_pencils_eos_pencpar(f,p,lpencil)
!
    endsubroutine calc_pencils_eos_std_ogrid
!***********************************************************************
    subroutine calc_pencils_eos_pencpar_ogrid(f,p,lpenc_loc)
!
!  Calculate EquationOfState pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      use Sub
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mfarray_ogrid) :: f
      type (pencil_case),                intent(OUT)  :: p
      logical, dimension(:),             intent(IN)   :: lpenc_loc
!
      real, dimension(nx_ogrid) :: tmp
      integer :: i,j
      real, dimension(:,:), pointer :: reference_state
!
!  Inverse cv and cp values.
!
      if (lpenc_loc(i_cv1)) p%cv1=cv1
      if (lpenc_loc(i_cp1)) p%cp1=cp1
      if (lpenc_loc(i_cv))  p%cv=1/cv1
      if (lpenc_loc(i_cp))  p%cp=1/cp1
      if (lpenc_loc(i_cp1tilde)) p%cp1tilde=cp1
!
      if (lpenc_loc(i_glnmumol)) p%glnmumol(:,:)=0.0
!
      select case (ieosvars)
!
!  Work out thermodynamic quantities for given lnrho or rho and TT.
!
      case (ilnrho_TT,irho_TT)
        if (lpenc_loc(i_TT))   p%TT=f(l1:l2,m,n,ieosvar2)
        if (lpenc_loc(i_TT1).or.lpenc_loc(i_hlnTT))  p%TT1=1/f(l1:l2,m,n,ieosvar2)
        if (lpenc_loc(i_lnTT).or.lpenc_loc(i_ss).or.lpenc_loc(i_ee)) &
            p%lnTT=log(f(l1:l2,m,n,ieosvar2))
        if (lpenc_loc(i_cs2))  p%cs2=cp*gamma_m1*f(l1:l2,m,n,ieosvar2)
        if (lpenc_loc(i_gTT))  call grad(f,ieosvar2,p%gTT)
        if (lpenc_loc(i_glnTT).or.lpenc_loc(i_hlnTT)) then
          do i=1,3; p%glnTT(:,i)=p%gTT(:,i)*p%TT1; enddo
        endif
        if (lpenc_loc(i_del2TT).or.lpenc_loc(i_del2lnTT)) &
            call del2(f,ieosvar2,p%del2TT)
        if (lpenc_loc(i_del2lnTT)) then
          tmp=0.0
          do i=1,3
            tmp=tmp+p%glnTT(:,i)**2
          enddo
          p%del2lnTT=p%del2TT*p%TT1-tmp
        endif
        if (lpenc_loc(i_hlnTT)) then
          call g2ij(f,iTT,p%hlnTT)
          do i=1,3; do j=1,3
            p%hlnTT(:,i,j)=p%hlnTT(:,i,j)*p%TT1-p%glnTT(:,i)*p%glnTT(:,j)
          enddo; enddo
        endif
        if (lpenc_loc(i_del6TT)) call del6(f,ieosvar2,p%del6TT)
        if (lpenc_loc(i_ss)) p%ss=cv*(p%lnTT-lnTT0-gamma_m1*(p%lnrho-lnrho0))
        if (lpenc_loc(i_pp)) p%pp=cv*gamma_m1*p%rho*p%TT
        if (lpenc_loc(i_ee)) p%ee=cv*exp(p%lnTT)
!
!  Work out thermodynamic quantities for given lnrho or rho and cs2.
!
      case (ilnrho_cs2,irho_cs2)
        if (leos_isentropic) then
          call fatal_error('calc_pencils_eos', &
              'leos_isentropic not implemented for ilnrho_cs2, try ilnrho_ss')
        elseif (leos_isothermal) then
          if (lpenc_loc(i_cs2)) p%cs2=cs20
          if (lpenc_loc(i_lnTT)) p%lnTT=lnTT0
          if (lpenc_loc(i_glnTT)) p%glnTT=0.0
          if (lpenc_loc(i_hlnTT)) p%hlnTT=0.0
          if (lpenc_loc(i_del2lnTT)) p%del2lnTT=0.0
          if (lpenc_loc(i_ss)) p%ss=-(cp-cv)*(p%lnrho-lnrho0)
          if (lpenc_loc(i_del2ss)) p%del2ss=-(cp-cv)*p%del2lnrho
          if (lpenc_loc(i_gss)) p%gss=-(cp-cv)*p%glnrho
          if (lpenc_loc(i_hss)) p%hss=-(cp-cv)*p%hlnrho
          if (lpenc_loc(i_pp)) p%pp=gamma1*p%rho*cs20
        elseif (leos_localisothermal) then
          if (lpenc_loc(i_cs2)) p%cs2=f(l1:l2,m,n,iglobal_cs2)
          if (lpenc_loc(i_lnTT)) call fatal_error('calc_pencils_eos', &
              'temperature not needed for localisothermal')
          if (lpenc_loc(i_glnTT)) &
              p%glnTT=f(l1:l2,m,n,iglobal_glnTT:iglobal_glnTT+2)
          if (lpenc_loc(i_hlnTT)) call fatal_error('calc_pencils_eos', &
              'no gradients yet for localisothermal')
          if (lpenc_loc(i_del2lnTT)) call fatal_error('calc_pencils_eos', &
              'no gradients yet for localisothermal')
          if (lpenc_loc(i_ss)) call fatal_error('calc_pencils_eos', &
              'entropy not needed for localisothermal')
          if (lpenc_loc(i_del2ss)) call fatal_error('calc_pencils_eos', &
              'no gradients yet for localisothermal')
          if (lpenc_loc(i_gss)) call fatal_error('calc_pencils_eos', &
              'entropy gradient not needed for localisothermal')
          if (lpenc_loc(i_hss)) call fatal_error('calc_pencils_eos', &
              'no gradients yet for localisothermal')
          if (lpenc_loc(i_pp)) p%pp=p%rho*p%cs2
        else
          call fatal_error('calc_pencils_eos', &
              'Full equation of state not implemented for ilnrho_cs2')
        endif
!
      case default
        call fatal_error('calc_pencils_eos','case not implemented yet')
      endselect
!
    endsubroutine calc_pencils_eos_pencpar_ogrid
!***********************************************************************
    subroutine calc_pencils_hydro_ogrid(f,p)
!
! Envelope adjusting calc_pencils_hydro_pencpar to the standard use with
! lpenc_loc=lpencil
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mfarray_ogrid),intent(in) :: f
      type (pencil_case),intent(OUT):: p
!
      call calc_pencils_hydro_nonlinear_ogrid(f,p,lpencil)
!
    endsubroutine calc_pencils_hydro_ogrid
!***********************************************************************
    subroutine calc_pencils_hydro_nonlinear_ogrid(f,p,lpenc_loc)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      use Deriv
      use Sub
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mfarray_ogrid),intent(in) :: f
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      real, dimension (nx_ogrid) :: tmp, tmp2
      integer :: i, j, ju, ij,jj,kk,jk
!
      intent(in) :: lpenc_loc
      intent(out):: p
! uu
      if (lpenc_loc(i_uu)) p%uu=f(l1:l2,m,n,iux:iuz)
! u2
      if (lpenc_loc(i_u2)) call dot2_mn(p%uu,p%u2)
! uij
      if (lpenc_loc(i_uij)) then 
        call gij(f,iuu,p%uij,1)
        if (lparticles_lyapunov) then
          jk=0
          do jj=1,3; do kk=1,3
            jk=jk+1
            f(l1:l2,m,n,iguij+jk-1) = p%uij(:,jj,kk)
          enddo;enddo
        endif
      endif
! divu
      if (lpenc_loc(i_divu)) call div_mn(p%uij,p%divu,p%uu)
! sij
      if (lpenc_loc(i_sij)) call traceless_strain(p%uij,p%divu,p%sij,p%uu,lshear_rateofstrain)
! sij2
      if (lpenc_loc(i_sij2)) call multm2_sym_mn(p%sij,p%sij2)
! oo (=curlu)
      oo: if (lpenc_loc(i_oo)) then
        if (ioo /= 0) then
          p%oo = f(l1:l2,m,n,iox:ioz)
        else
          call curl_mn(p%uij,p%oo,p%uu)
        endif
      endif oo
! o2
      if (lpenc_loc(i_o2)) call dot2_mn(p%oo,p%o2)
! ou
      if (lpenc_loc(i_ou)) call dot_mn(p%oo,p%uu,p%ou)
! ugu
      if (lpenc_loc(i_ugu)) then
        if (headtt.and.lupw_uu) print *,'calc_pencils_hydro: upwinding advection term'
        call u_dot_grad(f,iuu,p%uij,p%uu,p%ugu,UPWIND=lupw_uu)
      endif
! ugu2
      if (lpenc_loc(i_ugu2)) call dot2_mn(p%ugu,p%ugu2)
! ogu ... ogu2
      if (lpenc_loc(i_ogu)) call u_dot_grad(f,iuu,p%uij,p%oo,p%ogu,UPWIND=lupw_uu)
! u3u21, u1u32, u2u13, u2u31, u3u12, u1u23
      if (lpenc_loc(i_u3u21)) p%u3u21=p%uu(:,3)*p%uij(:,2,1)
      if (lpenc_loc(i_u1u32)) p%u1u32=p%uu(:,1)*p%uij(:,3,2)
      if (lpenc_loc(i_u2u13)) p%u2u13=p%uu(:,2)*p%uij(:,1,3)
      if (lpenc_loc(i_u2u31)) p%u2u31=p%uu(:,2)*p%uij(:,3,1)
      if (lpenc_loc(i_u3u12)) p%u3u12=p%uu(:,3)*p%uij(:,1,2)
      if (lpenc_loc(i_u1u23)) p%u1u23=p%uu(:,1)*p%uij(:,2,3)
!
! del2u, graddivu
!
      if (lpenc_loc(i_graddivu)) then
        if (headtt.or.ldebug) print*,'calc_pencils_hydro: call gij_etc'
        call gij_etc(f,iuu,p%uu,p%uij,p%oij,GRADDIV=p%graddivu)
      endif
      if (lpenc_loc(i_del2u)) then
        call curl_mn(p%oij,p%curlo,p%oo)
        p%del2u=p%graddivu-p%curlo
      endif
!
!del2uj, d^u/dx^2 etc
!
      if (lpenc_loc(i_d2uidxj)) &
        call d2fi_dxj(f,iuu,p%d2uidxj)
!
    endsubroutine calc_pencils_hydro_nonlinear_ogrid
!***********************************************************************
    subroutine calc_pencils_density_ogrid(f,p)
!
!  Calculate Density pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mfarray_ogrid),intent(in) :: f
      type (pencil_case) :: p
      intent(inout) :: f,p
!
!  Differentiate between log density and linear density.
!
      if (ldensity_nolog) then
        call calc_pencils_linear_density(f,p)
      else
        call fatal_error('calc_pencils_density_ogrid',&
            'Must use linear density for solid_cells_ogrid')
      endif
! ekin
      if (lpencil(i_ekin)) p%ekin=0.5*p%rho*p%u2
!
    endsubroutine calc_pencils_density_ogrid
!***********************************************************************
    subroutine calc_pencils_linear_density_ogrid(f,p)
!
!  Calculate Density pencils for linear density.
!  Most basic pencils should come first, as others may depend on them.
!
      use Sub, only: grad,dot,dot2,u_dot_grad,del2,del6,multmv,g2ij, dot_mn
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mfarray_ogrid),intent(in) :: f
      type (pencil_case) :: p
      intent(inout) :: f,p
!
      integer :: i
!
! rho
      p%rho=f(l1:l2,m,n,irho)
! rho1
      if (lpencil(i_rho1)) p%rho1=1.0/p%rho
! lnrho
      if (lpencil(i_lnrho))p%lnrho=log(p%rho)
! glnrho and grho
      if (lpencil(i_glnrho).or.lpencil(i_grho)) then
!
        call grad(f,irho,p%grho)
! 
        if (lpencil(i_glnrho)) then
          do i=1,3
            p%glnrho(:,i)=p%rho1*p%grho(:,i)
          enddo
        endif
      endif
! ugrho
      if (lpencil(i_ugrho)) &
        call u_dot_grad(f,ilnrho,p%grho,p%uu,p%ugrho,UPWIND=lupw_rho)
! glnrho2
      if (lpencil(i_glnrho2)) call dot2(p%glnrho,p%glnrho2)
! del2rho
      if (lpencil(i_del2rho)) then
        call del2(f,irho,p%del2rho)
        if (lreference_state) p%del2rho=p%del2rho+reference_state(:,iref_d2rho)
      endif
! del2lnrho
      if (lpencil(i_del2lnrho)) p%del2lnrho=p%rho1*p%del2rho-p%glnrho2
! sglnrho
      if (lpencil(i_sglnrho)) call multmv(p%sij,p%glnrho,p%sglnrho)
!
    endsubroutine calc_pencils_linear_density_ogrid
!***********************************************************************
    subroutine calc_pencils_viscosity_ogrid(f,p)
!
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-nov-04/anders: coded
!  18-may-12/MR: calculation of viscous heat for boussinesq added
!  14-oct-15/MR: corrected viscous force for slope-limited flux
!
      use Sub
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mfarray_ogrid),intent(in) :: f
      type (pencil_case) :: p
      intent(inout) :: f,p
!
      real, dimension (nx_ogrid,3) :: tmp,tmp2,tmp5,gradnu,sgradnu,gradnu_shock
      real, dimension (nx_ogrid) :: murho1,zetarho1,muTT,tmp3,tmp4,pnu,pnu_shock
      real, dimension (nx_ogrid) :: lambda_phi,prof,prof2,derprof,derprof2,qfvisc
      real, dimension (nx_ogrid) :: gradnu_effective,fac
      real, dimension (nx_ogrid,3) :: deljskl2,fvisc_nnewton2
!
      integer :: i,j,ju,ii,jj,kk,ll
!
!  Viscous force and viscous heating are calculated here (for later use).
!
      p%fvisc=0.0                               !!! not needed
!
!  viscous force: nu*del2v
!  -- not physically correct (no momentum conservation), but
!  numerically easy and in most cases qualitatively OK,
!  for boussinesq (divu=0) yet exact
!
      p%fvisc=p%fvisc+nu*p%del2u
!
    end subroutine calc_pencils_viscosity_ogrid
!***********************************************************************
    subroutine calc_viscous_force_ogrid(df,p)
!
!  Calculate viscous force term for right hand side of equation of motion.
!
      use Sub, only: cross, dot2
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mfarray_ogrid),intent(in) :: f
      type (pencil_case) :: p
      integer :: i
!
      intent (in) :: p
      intent (inout) :: df
!
!  Add viscosity to equation of motion
!
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%fvisc
!
    endsubroutine calc_viscous_force_ogrid
!***********************************************************************
    subroutine calc_pencils_energy_ogrid(f,p)
!
!  Calculate Energy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mfarray_ogrid),intent(in) :: f
      type (pencil_case) :: p
!
      integer :: j
!
      intent(in) :: f
      intent(inout) :: p
! Ma2
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
!
!  fpres (=pressure gradient force)
!
      if (lpencil(i_fpres)) then
        do j=1,3
          if (llocal_iso) then
            p%fpres(:,j)=-p%cs2*(p%glnrho(:,j)+p%glnTT(:,j))
          else
            p%fpres(:,j)=-p%cs2*p%glnrho(:,j)
          endif
        enddo
      endif
!
! tcond (dummy)
!
      if (lpencil(i_tcond)) p%tcond=0.
! sglnTT (dummy)
      if (lpencil(i_sglnTT)) p%sglnTT=0.
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_energy_ogrid
!***********************************************************************
    subroutine boundconds_x_ogrid(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in x, except for periodic part handled by communication.
!  Remark: boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-'corners').
!
!  For ogrids: Only boundary conditions at cylinder surface is set. The BC on
!  the 'top' is set by interpolation from cartesian grid, outside the timestep.
!
!  06-feb-17/Jorgen: Adapted from boundcond.f90 to be used for ogrids
!
      use EquationOfState
!
      real, dimension (:,:,:,:) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j, k
      logical :: ip_ok
      type (boundary_condition) :: bc
      integer :: tester
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  Boundary conditions in x.
!
      topbot='bot'; ip_ok=lfirst_proc_x

      tester=1
!
      if (ip_ok) then
        do j=ivar1,ivar2
          if(tester)
            bc%bcname=bcx12(j,k)
            bc%ivar=j
            bc%location=(((k-1)*2)-1)   ! -1/1 for x bot/top
            bc%value1=fbcx(j,k)
            bc%value2=fbcx(j,k)
            bc%done=.false.
!
            call special_boundconds(f,bc)
!
            if (.not.bc%done) then
              call fatal_error_local("boundconds_x",'no such boundary condition')
            endif
          else
            ! TODO
            ! Do something smart to set body confined boundary condition properly
          endif
        enddo
      endif
!
    endsubroutine boundconds_x_ogrid
!***********************************************************************
    subroutine boundconds_y_ogrid(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in y, except for periodic part handled by communication.
!  Remark: boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-'corners').
!
!  06-feb-17/Jorgen: Adapted from boundcond.f90 to be used for ogrids where the
!                    y-dir is always periodic
!
      real, dimension (:,:,:,:) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j, k
      logical :: ip_ok
      character (len=bclen) :: topbot
      type (boundary_condition) :: bc
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  Boundary conditions in y
!  Periodic, with y being the theta direction for the cylindrical grid
!
      if (nprocy==1) then
        do j=ivar1,ivar2
!  Bottom boundary
          f(:,1:m1_ogrid-1,:,j) = f(:,m2i_ogrid:m2_ogrid,:,j)
!  Top boundary
          f(:,m2_ogrid+1:,:,j) = f(:,m1_ogrid:m1i_ogrid,:,j)
        enddo
      endif
!
    endsubroutine boundconds_y_ogrid
!***********************************************************************
    subroutine boundconds_z_ogrid(f,ivar1_opt,ivar2_opt)
!
!  Boundary conditions in z, except for periodic part handled by communication.
!  Remark: boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-'corners').
!
!  06-feb-17/Jorgen: Adapted from boundcond.f90 to be used for ogrids where the
!                    y-dir is always periodic as long as nzgrid=/0
!
      real, dimension (:,:,:,:) :: f
      integer, optional :: ivar1_opt, ivar2_opt
      real, dimension (size(f,4)) :: fbcz_zero
      integer :: ivar1, ivar2, j, k
      logical :: ip_ok
      character (len=bclen) :: topbot
      type (boundary_condition) :: bc
!
      if (ldebug) print*,'boundconds_z: ENTER: boundconds_z'
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      select case (nzgrid)
!
      case (1)
        if (ldebug) print*,'boundconds_z: no z-boundary'
!
!  Boundary conditions in z
!
      case default
        if (nprocz==1) then
          do j=ivar1,ivar2
!  Bottom boundary
            f(:,:,1:n1_ogrid-1,j) = f(:,:,n2i_ogrid:n2_ogrid,j)
!  Top boundary
            f(:,:,n2_ogrid+1:,j) = f(:,:,n1_ogrid:n1i_ogrid,j)
          enddo
        endif
      endselect
!
    endsubroutine boundconds_z_ogrid
!***********************************************************************
    subroutine initiate_isendrcv_bdry_ogrid(f,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values. Called in the beginning of pde.
!  Does not wait for the receives to finish (done in finalize_isendrcv_bdry)
!  leftneigh and rightneigh are initialized by mpicomm_init.
!
!  07-feb-17/Jorgen: adapted from mpicomm.f90
!
      real, dimension(:,:,:,:), intent(inout):: f
      integer, optional,        intent(in)   :: ivar1_opt, ivar2_opt
!
      integer :: tolowyr,touppyr,tolowys,touppys,tolowzr,touppzr,tolowzs,touppzs ! msg. tags placeholders
      integer :: TOllr,TOulr,TOuur,TOlur,TOlls,TOuls,TOuus,TOlus                 ! placeholder tags
      integer :: ivar1, ivar2, nbufy, nbufz, nbufyz, mxl, comm, bufact, dir
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
      if (ivar2==0) return
!
      mxl=size(f,1)
!
!  Set communication across x-planes.
!
      if (nprocx>1) call isendrcv_bdry_x(f,ivar1_opt,ivar2_opt)
!
!  Set communication across y-planes.
!  Standard message tags from defaults for surfaces and yz-corners.
!
      tolowyr=tolowy; tolowzr=tolowz; TOlls=TOll; TOllr=TOll
      tolowys=tolowy; tolowzs=tolowz; TOlus=TOlu; TOlur=TOlu
      touppyr=touppy; touppzr=touppz; TOuus=TOuu; TOuur=TOuu
      touppys=touppy; touppzs=touppz; TOuls=TOul; TOulr=TOul
!
!  Allocate and send/receive buffers across y-planes
!
      bufact=mxl*nghost*(ivar2-ivar1+1)

      if (nprocy>1) then
!
!  Internal, periodic exchange y-plane buffers.
!
        lbufyo(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m1i_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2)         !(lower y-zone)
        comm=MPI_COMM_GRID

        nbufy=bufact*bufsizes_yz(INYL,IRCV)
        call MPI_IRECV(lbufyi(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
                       ylneigh,touppyr,comm,irecv_rq_fromlowy,mpierr)

        nbufy=bufact*bufsizes_yz(INYL,ISND)
        call MPI_ISEND(lbufyo(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
                       ylneigh,tolowys,comm,isend_rq_tolowy,mpierr)

        ubufyo(:,:,:,ivar1:ivar2)=f(:,m2i_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(upper y-zone)
        comm=MPI_COMM_GRID
!
        nbufy=bufact*bufsizes_yz(INYU,IRCV)
        call MPI_IRECV(ubufyi(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
                       yuneigh,tolowyr,comm,irecv_rq_fromuppy,mpierr)

        nbufy=bufact*bufsizes_yz(INYU,ISND)
        call MPI_ISEND(ubufyo(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
                       yuneigh,touppys,comm,isend_rq_touppy,mpierr)
      endif
!
!  Set communication across z-planes.
!
      if (nprocz>1) then
!        
        lbufzo(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m2_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2) !lower z-planes
        comm=MPI_COMM_GRID

        nbufz=bufact*bufsizes_yz(INZL,IRCV)
        call MPI_IRECV(lbufzi(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
                       zlneigh,touppzr,comm,irecv_rq_fromlowz,mpierr)

        nbufz=bufact*bufsizes_yz(INZL,ISND)
        call MPI_ISEND(lbufzo(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
                       zlneigh,tolowz,comm,isend_rq_tolowz,mpierr)
!
        ubufzo(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m2_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2) !upper z-planes
        comm=MPI_COMM_GRID

        nbufz=bufact*bufsizes_yz(INZU,IRCV)
        call MPI_IRECV(ubufzi(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
                       zuneigh,tolowzr,comm,irecv_rq_fromuppz,mpierr)

        nbufz=bufact*bufsizes_yz(INZU,ISND)
        call MPI_ISEND(ubufzo(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
                       zuneigh,touppz,comm,isend_rq_touppz,mpierr)

      endif
!
!  The four corners (in counter-clockwise order).
!  (NOTE: this should work even for nprocx>1)
!
      if ((nprocz>1).and.(nprocy>1)) then
!
!  Internal and periodic yz-corner buffers.
!
        bufact=mxl*(ivar2-ivar1+1)
!
!  Lower y, lower z.
!
        llbufo(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m1i_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2)
        comm=MPI_COMM_GRID
!
        nbufyz=bufact*product(bufsizes_yz_corn(:,INLL,IRCV))
        if (llcornr>=0) call MPI_IRECV(llbufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       llcornr,TOuur,comm,irecv_rq_FRll,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INLL,ISND))
        if (llcorns>=0) call MPI_ISEND(llbufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       llcorns,TOlls,comm,isend_rq_TOll,mpierr)
!
!  Upper y, lower z.
!
        ulbufo(:,:,:,ivar1:ivar2)=f(:,m2i_ogrid:m2_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2)
        comm=MPI_COMM_GRID

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUL,IRCV))
        if (ulcornr>=0) call MPI_IRECV(ulbufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       ulcornr,TOlur,comm,irecv_rq_FRul,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUL,ISND))
        if (ulcorns>=0) call MPI_ISEND(ulbufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       ulcorns,TOuls,comm,isend_rq_TOul,mpierr)
!
!  Upper y, upper z.
!
        uubufo(:,:,:,ivar1:ivar2)=f(:,m2i_ogrid:m2_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2)
        comm=MPI_COMM_GRID

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUU,IRCV))
        if (uucornr>=0) call MPI_IRECV(uubufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       uucornr,TOllr,comm,irecv_rq_FRuu,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUU,ISND))
        if (uucorns>=0) call MPI_ISEND(uubufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       uucorns,TOuus,comm,isend_rq_TOuu,mpierr)
!
!  Lower y, upper z.
!
        lubufo(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m1i_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2)
        comm=MPI_COMM_GRID

        nbufyz=bufact*product(bufsizes_yz_corn(:,INLU,IRCV))
        if (lucornr>=0) call MPI_IRECV(lubufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       lucornr,TOulr,comm,irecv_rq_FRlu,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INLU,ISND))
        if (lucorns>=0) call MPI_ISEND(lubufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       lucorns,TOlus,comm,isend_rq_TOlu,mpierr)
!
      endif
!
    endsubroutine initiate_isendrcv_bdry_ogrid
!***********************************************************************
    subroutine finalize_isendrcv_bdry_ogrid(f,ivar1_opt,ivar2_opt)
!
!  Make sure the communications initiated with initiate_isendrcv_bdry are
!  finished and insert the just received boundary values.
!  Receive requests do not need to (and on OSF1 cannot) be explicitly
!  freed, since MPI_Wait takes care of this.
!
!  07-feb-17/Jorgen: Adapted from mpicomm.f90
!
      use General, only: transpose_mn, notanumber

      real, dimension(:,:,:,:), intent(inout):: f
      integer, optional,        intent(in)   :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
      if (ivar2==0) return
!
!  1. wait until data received
!  2. set ghost zones
!  3. wait until send completed, will be overwritten in next time step
!
!  Communication across y-planes (includes periodic bc)
!
      if (nprocy>1) then
        call MPI_WAIT(irecv_rq_fromuppy,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowy,irecv_stat_fl,mpierr)
        do j=ivar1,ivar2
          if (.not. lfirst_proc_y .or. bcy12(j,1)=='p') then  
            f(:,1:m1_ogrid-1,n1_ogrid:n2_ogrid,j)=lbufyi(:,:,:,j)       ! set lower buffer
          endif

          if (.not. llast_proc_y .or. bcy12(j,2)=='p') then
            f(:,m2_ogrid+1:,n1_ogrid:n2_ogrid,j)=ubufyi(:,:,:,j)        ! set upper buffer
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowy,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppy,isend_stat_tu,mpierr)
      endif
!
!  Communication in z (includes periodic bc)
!
      if (nprocz>1) then
        call MPI_WAIT(irecv_rq_fromuppz,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowz,irecv_stat_fl,mpierr)
        do j=ivar1,ivar2
          if (.not. lfirst_proc_z .or. bcz12(j,1)=='p') then
            f(:,m1_ogrid:m2_ogrid,1:n1_ogrid-1,j)=lbufzi(:,:,:,j)       ! set lower buffer
          endif

          if (.not. llast_proc_z .or. bcz12(j,2)=='p') then 
            f(:,m1_ogrid:m2_ogrid,n2_ogrid+1:,j)=ubufzi(:,:,:,j)        ! set upper buffer
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowz,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppz,isend_stat_tu,mpierr)
      endif
!
!  The four yz-corners (in counter-clockwise order)
!
       if (nprocz>1.and.nprocy>1) then

        if (uucornr>=0) call MPI_WAIT(irecv_rq_FRuu,irecv_stat_Fuu,mpierr)
        if (lucornr>=0) call MPI_WAIT(irecv_rq_FRlu,irecv_stat_Flu,mpierr)
        if (llcornr>=0) call MPI_WAIT(irecv_rq_FRll,irecv_stat_Fll,mpierr)
        if (ulcornr>=0) call MPI_WAIT(irecv_rq_FRul,irecv_stat_Ful,mpierr)

        do j=ivar1,ivar2
!
!  Set ll corner
!
          if (.not. lfirst_proc_z .or. bcz12(j,1)=='p') then
            if  (.not.lfirst_proc_y.or.bcy12(j,1)=='p') then    ! inner or periodic proc boundaries
              f(:,1:m1_ogrid-1,1:n1_ogrid-1,j)=llbufi(:,:,:,j)               ! fill lower left corner
            endif
!
!  Set ul corner
!
            if (.not.llast_proc_y .or.bcy12(j,2)=='p') then    ! inner or periodic proc boundaries
                f(:,m2_ogrid+1:,1:n1_ogrid-1,j)=ulbufi(:,:,:,j)                ! fill lower right corner
            endif
          endif
!
!  Set uu corner
!
          if (.not. llast_proc_z .or. bcz12(j,2)=='p') then
            if (.not.llast_proc_y.or.bcy12(j,2)=='p') then    ! inner or periodic proc boundaries
              f(:,m2+1:,n2+1:,j)=uubufi(:,:,:,j)                ! fill upper right corner
            endif
!
!  Set lu corner
!
            if (.not. lfirst_proc_y .or. bcy12(j,1)=='p') then
              f(:,1:m1-1,n2+1:,j)=lubufi(:,:,:,j)                ! fill upper left corner
            endif
          endif
        enddo
        
        if (llcorns>=0) call MPI_WAIT(isend_rq_TOll,isend_stat_Tll,mpierr)
        if (ulcorns>=0) call MPI_WAIT(isend_rq_TOul,isend_stat_Tul,mpierr)
        if (uucorns>=0) call MPI_WAIT(isend_rq_TOuu,isend_stat_Tuu,mpierr)
        if (lucorns>=0) call MPI_WAIT(isend_rq_TOlu,isend_stat_Tlu,mpierr)

      endif
!
!  make sure the other processors don't carry on sending new data
!  which could be mistaken for an earlier time
!
     call mpibarrier
!call mpifinalize
!stop
!
    endsubroutine finalize_isendrcv_bdry_ogrid
!***********************************************************************
    subroutine isendrcv_bdry_x_ogrid(f,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values for x-direction. Sends and receives
!  before continuing to y and z boundaries, as this allows the edges
!  of the grid to be set properly.
!
!   2-may-09/anders: coded
!
      real, dimension(:,:,:,:), intent(inout) :: f
      integer, intent(in), optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbufx, j
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      if (nprocx>1) then

        lbufxo(:,:,:,ivar1:ivar2)=f(l1_ogrid:l1i_ogrid,m_ogrid1:m2,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(lower x-zone)
        ubufxo(:,:,:,ivar1:ivar2)=f(l2i_ogrid:l2_ogrid,m_ogrid1:m2,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(upper x-zone)
        nbufx=ny_ogrid*nz_ogrid*nghost*(ivar2-ivar1+1)

        call MPI_IRECV(ubufxi(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xuneigh,tolowx,MPI_COMM_GRID,irecv_rq_fromuppx,mpierr)
        call MPI_IRECV(lbufxi(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xlneigh,touppx,MPI_COMM_GRID,irecv_rq_fromlowx,mpierr)
        call MPI_ISEND(lbufxo(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xlneigh,tolowx,MPI_COMM_GRID,isend_rq_tolowx,mpierr)
        call MPI_ISEND(ubufxo(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xuneigh,touppx,MPI_COMM_GRID,isend_rq_touppx,mpierr)
        call MPI_WAIT(irecv_rq_fromuppx,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowx,irecv_stat_fl,mpierr)
!
!  Inner communication or (shear-)periodic boundary conditions in x
!  MR: Communication should only happen under these conditions.
!
        do j=ivar1,ivar2
          if (.not. lfirst_proc_x .or. bcx12(j,1)=='p' .or. &
              (bcx12(j,1)=='she'.and.nygrid==1)) then
            f( 1:l1_ogrid-1,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j)=lbufxi(:,:,:,j)  !!(set lower buffer)
          endif
          if (.not. llast_proc_x .or. bcx12(j,2)=='p' .or. &
              (bcx12(j,2)=='she'.and.nygrid==1)) then
            f(l2_ogrid+1:,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j)=ubufxi(:,:,:,j)  !!(set upper buffer)
          endif
        enddo

        call MPI_WAIT(isend_rq_tolowx,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppx,isend_stat_tu,mpierr)

      endif
!
    endsubroutine isendrcv_bdry_x_ogrid
!***********************************************************************
  end module Solid_Cells

