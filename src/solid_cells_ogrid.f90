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
!
!   New solid cells module
!   Overlapping cylindrical grid (O-grid) around cylinder or sphere, 
!   coupled to cartesian grid outside the o-grid by interpolation.
!
!   Very fine resolution of boundary layer
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

  interface dot2_ogrid
    module procedure dot2_mn_ogrid
    module procedure dot2_0_ogrid
  endinterface
  interface u_dot_grad_ogrid
    module procedure u_dot_grad_scl_ogrid
    module procedure u_dot_grad_vec_ogrid
  endinterface

!  Cylinder parameters
  character(len=1) :: flow_direction
  real, parameter :: cylinder_radius=0.                        ! Set in start.in
  real, parameter :: cylinder_temp=703.0                       ! Set in start.in
  real, parameter :: cylinder_xpos=0., cylinder_ypos=0.        ! Set in start.in
  real, parameter :: cylinder_zpos=0.                          ! Set in start.in
  real, parameter :: skin_depth=0.                             ! Set in start.in
  real, init_uu=0., ampl_noise=0.                              ! Set in start.in
  character(len=labellen) :: initsolid_cells='cylinderstream_x'! Set in start.in

  real, dimension(3) :: xyz0_ogrid, Lxyz_ogrid, xorigo_ogrid

!  Fundamental grid parameters
  real :: r_ogrid=2*cylinder_radius                                  ! Set in start.in?
  character (len=labellen), dimension(3) :: grid_func_ogrid='linear' ! Set in start.in


!***************************************************
! Pencil case ogrid
  integer, parameter :: npencils_ogrid=37
  type pencil_case_ogrid
    real, dimension (nx_ogrid)     :: x_mn    
    real, dimension (nx_ogrid)     :: y_mn    
    real, dimension (nx_ogrid)     :: z_mn    
    real, dimension (nx_ogrid)     :: rcyl_mn 
    real, dimension (nx_ogrid)     :: phi_mn  
    real, dimension (nx_ogrid)     :: rcyl_mn1
    real, dimension (nx_ogrid,3)   :: fpres
    real, dimension (nx_ogrid,3)   :: fvisc
    real, dimension (nx_ogrid)     :: ugu
    real, dimension (nx_ogrid)     :: rho 
    real, dimension (nx_ogrid)     :: rho1 
    real, dimension (nx_ogrid)     :: lnrho
    real, dimension (nx_ogrid,3)   :: grho
    real, dimension (nx_ogrid,3)   :: glnrho
    real, dimension (nx_ogrid)     :: ugrho   
    real, dimension (nx_ogrid,3)   :: sglnrho 
    real, dimension (nx_ogrid,3)   :: uu
    real, dimension (nx_ogrid)     :: u2
    real, dimension (nx_ogrid,3,3) :: uij
    real, dimension (nx_ogrid)     :: divu
    real, dimension (nx_ogrid,3,3) :: sij
    real, dimension (nx_ogrid)     :: sij2
    real, dimension (nx_ogrid,3)   :: ugu
    real, dimension (nx_ogrid)     :: ugu2
    real, dimension (nx_ogrid,3)   :: del2u
    real, dimension (nx_ogrid)     :: cv1
    real, dimension (nx_ogrid)     :: cp1
    real, dimension (nx_ogrid)     :: cv
    real, dimension (nx_ogrid)     :: cp
    real, dimension (nx_ogrid)     :: divu      
    real, dimension (nx_ogrid)     :: TT
    real, dimension (nx_ogrid)     :: TT1
    real, dimension (nx_ogrid)     :: cs2
    real, dimension (nx_ogrid)     :: gTT
    real, dimension (nx_ogrid)     :: del2TT
    real, dimension (nx_ogrid)     :: pp
    real, dimension (nx_ogrid)     :: ee
  endtype pencil_case_ogrid
  
  integer :: i_og_x_mn    =1
  integer :: i_og_y_mn    =2
  integer :: i_og_z_mn    =3
  integer :: i_og_rcyl_mn =4
  integer :: i_og_phi_mn  =5
  integer :: i_og_rcyl_mn1=6
  integer :: i_og_fpres   =7
  integer :: i_og_fvisc   =8
  integer :: i_og_ugu     =9
  integer :: i_og_rho     =10
  integer :: i_og_rho1    =11
  integer :: i_og_lnrho   =12
  integer :: i_og_grho    =13
  integer :: i_og_glnrho  =14
  integer :: i_og_ugrho   =15
  integer :: i_og_sglnrho =16
  integer :: i_og_uu      =17
  integer :: i_og_u2      =18
  integer :: i_og_uij     =19
  integer :: i_og_divu    =20
  integer :: i_og_sij     =21
  integer :: i_og_sij2    =22
  integer :: i_og_ugu     =23
  integer :: i_og_ugu2    =24
  integer :: i_og_del2u   =25
  integer :: i_og_cv1     =26
  integer :: i_og_cp1     =27
  integer :: i_og_cv      =28
  integer :: i_og_cp      =29
  integer :: i_og_divu    =30
  integer :: i_og_TT      =31
  integer :: i_og_TT1     =32
  integer :: i_og_cs2     =33
  integer :: i_og_gTT     =34
  integer :: i_og_del2TT  =35
  integer :: i_og_pp      =36
  integer :: i_og_ee      =37

  character (len=15), parameter, dimension(npencils_ogrid) :: pencil_names_ogrid = &
    (/ 'x_mn          ', 'y_mn          ', 'z_mn          ', 'rcyl_mn       '  &
     , 'phi_mn        ', 'rcyl_mn1      ', 'fpres         ', 'fvisc         '  &
     , 'ugu           ', 'rho           '  &
     , 'rho1          ', 'lnrho         ', 'grho          ', 'glnrho        '  &
     , 'ugrho         ', 'sglnrho       '  &
     , 'uu            ', 'u2            ', 'uij           ', 'divu          '  &
     , 'sij           ', 'sij2          ', 'ugu           ', 'ugu2          '  &
     , 'del2u         ', 'divu          '  &
     , 'cv1           ', 'cp1           ', 'cp            ', 'cv            '  &
     , 'TT            ', 'TT1           '  &
     , 'cs2           ', 'gTT           ', 'del2TT        '  &
     , 'pp            ', 'ee            ' /)
  logical,dimension(npencils_ogrid):: lpencil = .true.
!***************************************************

!  Pencils to be used for curvelinear grid computations
  type(pencil_case_ogrid) p_ogrid 

! PARAMETERS NECESSARY FOR GRID CONSTRUCTION 
! TODO: Clean up
!  Global ogrid (TODO: NEEDED?)
  integer, parameter :: mxgrid_ogrid=nxgrid_ogrid+2*nghost
  integer, parameter :: mygrid_ogrid=nygrid_ogrid+2*nghost
  integer, parameter :: mzgrid_ogrid=nzgrid_ogrid+2*nghost
  real, dimension(nxgrid_ogrid) :: xgrid_ogrid, dx1grid_ogrid, dxtgrid_ogrid
  real, dimension(nygrid_ogrid) :: ygrid_ogrid, dy1grid_ogrid, dytgrid_ogrid
  real, dimension(nzgrid_ogrid) :: zgrid_ogrid, dz1grid_ogrid, dztgrid_ogrid
  real, dimension(mxgrid_ogrid) :: xglobal_ogrid
  real, dimension(mygrid_ogrid) :: yglobal_ogrid
  real, dimension(mzgrid_ogrid) :: zglobal_ogrid
!  Local ogrid and derivatives
  integer, parameter :: nx_ogrid=nxgrid_ogrid/nprocx,ny_ogrid=nygrid_ogrid/nprocy,nz_ogrid=nzgrid_ogrid/nprocz
  integer, parameter :: mx_ogrid=nx_ogrid+2*nghost,l1_ogrid=1+nghost,l2_ogrid=mx_ogrid-nghost
  integer, parameter :: my_ogrid=ny_ogrid+2*nghost,m1_ogrid=1+nghost,m2_ogrid=my_ogrid-nghost
  integer, parameter :: mz_ogrid=nz_ogrid+2*nghost,n1_ogrid=1+nghost,n2_ogrid=mz_ogrid-nghost
  real, dimension (mx_ogrid) :: x_ogrid,dx_1_ogrid,dx2_ogrid,dx_tilde_ogrid,xprim_ogrid
  real, dimension (my_ogrid) :: y_ogrid,dy_1_ogrid,dy2_ogrid,dy_tilde_ogrid,yprim_ogrid
  real, dimension (mz_ogrid) :: z_ogrid,dz_1_ogrid,dz2_ogrid,dz_tilde_ogrid,zprim_ogrid

!  Grid properties computed in grid routines and used in external routines
  real :: dxmin_ogrid,dxmax_ogrid
  logical, dimension(3) :: lequidist_ogrid
  real, dimension (nx_ogrid) :: rcyl_mn_ogrid,rcyl_mn1_ogrid,rcyl_mn2_ogrid,rcyl_weight_ogrid
!  Grid properties computed in grid routines and used other grid routines by call from external routines
  real, dimension (nx_ogrid,3) :: dline_1_ogrid
  real :: dx_ogrid,dy_ogrid,dz_ogrid
  real, dimension(3) :: xyz_star_ogrid
  real, dimension(-nghost:nghost,2) :: coeffs_1_x_ogrid
  real, dimension(-nghost:nghost,2) :: coeffs_1_y_ogrid
  real, dimension(-nghost:nghost,2) :: coeffs_1_z_ogrid
  real, dimension(3) :: coeff_grid_o=1.0
  real, dimension (-nghost:nghost) :: dx2_bound_ogrid=0., dy2_bound_ogrid=0., dz2_bound_ogrid=0.
!  Grid properties computed in grid routines, but not used in current implementation
  real, dimension (mx_ogrid) :: dVol_x_ogrid,dVol1_x_ogrid
  real, dimension (my_ogrid) :: dVol_y_ogrid,dVol1_y_ogrid
  real, dimension (mz_ogrid) :: dVol_z_ogrid,dVol1_z_ogrid
  real, dimension (nz_ogrid) :: dxyz_2_ogrid,dxyz_4_ogrid,dxyz_6_ogrid,dVol_ogrid
  real, dimension(0:nprocx) :: procx_bounds_ogrid
  real, dimension(0:nprocy) :: procy_bounds_ogrid
  real, dimension(0:nprocz) :: procz_bounds_ogrid
  real :: box_volume_ogrid=1.0 
  integer :: lpoint_ogrid=(mx_ogrid+1)/2, mpoint_ogrid=(my_ogrid+1)/2, npoint_ogrid=(mz_ogrid+1)/2
  integer :: lpoint2_ogrid=(mx_ogrid+1)/4,mpoint2_ogrid=(my_ogrid+1)/4,npoint2_ogrid=(mz_ogrid+1)/4

  real, dimension (nx_ogrid) :: dxmax_pencil_ogrid,dxmin_pencil_ogrid  ! Possible remove if no heat flux or viscosity module
  
  ! For mn-loop over ogrid pencils
  integer :: n_ogrid, m_ogrid
  ! For time-step
  real :: t_ogrid
  integer :: lfirst_ogrid, llast_ogrid

  real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray) ::  f_ogrid
  real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mvar)    :: df_ogrid
  

!  Read start.in file
  namelist /solid_cells_init_pars/ &
      cylinder_temp, ncylinders, cylinder_radius, cylinder_xpos, &
      cylinder_ypos, cylinder_zpos, flow_direction, skin_depth &
      initsolid_cells, init_uu

!  Read run.in file
  namelist /solid_cells_run_pars/  &
    
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
!  Define the geometry of the solids and construct the grid around the solid
!  object.
!  Currently only a single solid object with overlapping grid is implemented.
!  Use solid_cells.f90 without overlapping grids if several solids inside the
!  flow domain is desired.
!
!  XX-feb-17/Jorgen: Coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer :: icyl, isph,i

      real, dimension(3) :: xyz0_ogrid,Lxyz_ogrid

      if (cylinder_radius(icyl) <= 0) then
        call fatal_error('initialize_solid_cells_ogrid', &
            'All cylinders must have non-zero radii!')
      endif
!  Interlap parameters of the ogrid
      xyz0_ogrid(1)=cylinder_radius
      xyz0_ogrid(2)=-pi
      xyz0_ogrid(3)=xyz0(3)
      Lxyz_ogrid(1)=r_ogrid-xyz0_ogrid(1)
      Lxyz_ogrid(2)=2*pi
      Lxyz_ogrid(2)=Lxyz(3)
!  Position related to the flow domain and cartesian grid
      xorigo_ogrid(1) = cylinder_xpos(icyl)
      xorigo_ogrid(2) = cylinder_ypos(icyl)
      xorigo_ogrid(3) = cylinder_zpos(icyl)
!
!  Give all points inside the object a value of zero for all variables in the
!  cartesian array of unknowns
!
      do i = l1,l2
        do k = m1,m2
          if (sqrt((x(i)-xorigo_ogrid(1)**2+(y(i)-xorigo_ogrid(2)**2) < xyz0_ogrid(1)) then
            f(i,j,:,:)=0.
          endif
        enddo
      enddo
!
! Try to find flow direction
!
      flow_dir = 0
      if (fbcx(1,1) > 0) then; flow_dir = 1
      elseif (fbcx(1,2) < 0) then; flow_dir = -1
      elseif (fbcy(2,1) > 0) then; flow_dir = 2
      elseif (fbcy(2,2) < 0) then; flow_dir = -2
      elseif (fbcz(3,1) > 0) then; flow_dir = 3
      elseif (fbcz(3,2) < 0) then; flow_dir = -3
      endif
      if (flow_dir /= 0) then
        if (lroot) then
          print*,'By using fbc[x,y,z] I found the flow direction to be in the ', &
              flow_dir,' direction.'
        endif
      else
        do i = 1,3
          if (lperi(i)) then
            if (.not. lperi(mod(i,3)+1) .and. .not. lperi(mod(i+1,3)+1)) then
              flow_dir = i
              if (lroot) then
                print*,'By using lperi I found the flow direction to be in the ', &
                    flow_dir,' direction.'
              endif
            endif
          endif
        enddo
        if (lset_flow_dir) flow_dir = flow_dir_set
        if (flow_dir == 0) then
          call fatal_error('initialize_solid_cells', &
              'I was not able to determine the flow direction!')
        endif
      endif
!
! Find inlet temperature
!
      if (ilnTT /= 0) then
        if (flow_dir == 1) T0 = fbcx(ilnTT,1)
        if (flow_dir == -1) T0 = fbcx(ilnTT,2)
        if (flow_dir == 2) T0 = fbcy(ilnTT,1)
        if (flow_dir == -2) T0 = fbcy(ilnTT,2)
        if (flow_dir == 3) T0 = fbcz(ilnTT,1)
        if (flow_dir == -3) T0 = fbcz(ilnTT,2)
        if (.not. ltemperature_nolog) T0 = exp(T0)
      endif
!
!  Initialize the pencils overlapping grid
!
      call initialize_pencils_ogrid(p,0.0)
!
!  Construct overlapping grid
!
      call construct_grid_ogrid
!
!  Initialize overlapping grid
!
      call initialize_grid_ogrid
!
    end subroutine initialize_solid_cells
!***********************************************************************
    subroutine init_solid_cells(f)
!
!  Initial conditions for cases where we have solid structures in the domain.
!  Typically the flow field is set such that we have no-slip conditions
!  at the solid structure surface.
!
!  This routine is adapted to resolve the solid structures with overlapping 
!  curvelinear grids. Hence, the f-array, which is the array of flow variables
!  on the carteisian grid, must use potential flow to set the ghost zones that
!  will be set by intepolation for t>0.
! 
!  Only implemented for a single solid at present!
!
!  XX-feb-17/Jorgen: Coded
!
      use Initcond, only: gaunoise
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real :: a2, rr2, rr2_low, rr2_high
      real :: wall_smoothing,wall_smoothing_temp
      real :: Lorth,flowx,flowy,shift_flow,shift_orth,flow_r,orth_r
      integer i,j,cyl,iflow,iorth
!
!  Set up variables in the appropriate direction of the flow
!
      if(initsolid_cells='cylinderstream_x') then
        iflow= iux
        iorth= iuy
        Lorth= Lxyz(2)
        flowx=1.
        flowy=0.
        if (xorigo_ogrid(2) /= 0) then
          print*,'When using cylinderstream_x all cylinders must have'
          print*,'zero offset in y-direction!'
          call fatal_error('init_solid_cells:','')
        endif
      elseif(initsolid_cells='cylinderstream_y') then
        iflow= iuy
        iorth= iux
        Lorth= Lxyz(1)
        flowx=0.
        flowy=1.
        if (xorigo_ogrid(1) /= 0) then
          print*,'When using cylinderstream_y all cylinders must have'
          print*,'zero offset in x-direction!'
          call fatal_error('init_solid_cells:','')
        endif
      else
        if (lroot) print*,'No such value for init_solid_cells:', &
            initsolid_cells)
          call fatal_error('init_solid_cells','initialization of cylinderstream')
        endif
      endif
!
!  Stream functions for flow around a cylinder as initial condition.
!
      call gaunoise(ampl_noise,f,iux,iuz)
      f(:,:,:,iflow) = f(:,:,:,iflow)+init_uu
      shift_flow = 0
      do i = l1,l2
        do j = m1,m2
! Choose correct points depending on flow direction
          flow_r=(x(i)-xorigo_ogrid(1))*flowx+(y(j)-xorigo_ogrid(2))*flowx
          orth_r=(x(i)-xorigo_ogrid(1))*flowy+(y(j)-xorigo_ogrid(2))*flowx
          rr2 = flow_r**2+orth_r**2
          if (rr2 > a2) then
            do cyl = 0,100
              if (cyl == 0) then
                wall_smoothing = 1-exp(-(rr2-a2)/skin_depth**2)
                f(i,j,:,iorth) = f(i,j,:,iorth)-init_uu* &
                  2*flow_r*orth_r*a2/rr2**2*wall_smoothing
                f(i,j,:,iflow) = f(i,j,:,iflow)+init_uu* &
                  (0. - a2/rr2 + 2*orth_r**2*a2/rr2**2)*wall_smoothing
                if (ilnTT /= 0) then
                  wall_smoothing_temp = 1-exp(-(rr2-a2)/(sqrt(a2))**2)
                  f(i,j,:,ilnTT) = wall_smoothing_temp*f(i,j,:,ilnTT) &
                    +cylinder_temp*(1-wall_smoothing_temp)
                  f(i,j,:,ilnrho) = f(l2,m2,n2,ilnrho) &
                    *f(l2,m2,n2,ilnTT)/f(i,j,:,ilnTT)
                endif
              else
                shift_orth = cyl*Lorth
                rr2_low = (flow_r+shift_flow)**2+(orth_r+shift_orth)**2
                rr2_high = (flow_r-shift_flow)**2+(orth_r-shift_orth)**2
                f(i,j,:,iflow) = f(i,j,:,iflow)+init_uu*( &
                    +2*(orth_r-shift_orth)**2*a2/rr2_high**2-a2/rr2_high &
                    +2*(orth_r+shift_orth)**2*a2/rr2_low**2 -a2/rr2_low)
                f(i,j,:,iorth) = f(i,j,:,iorth)-init_uu*( &
                    +2*(flow_r-shift_flow)*(orth_r-shift_orth) &
                    *a2/rr2_high**2 &
                    +2*(flow_r+shift_flow)*(orth_r+shift_orth) &
                    *a2/rr2_low**2)
              endif
            enddo
          else
            if (ilnTT /= 0) then
              f(i,j,:,ilnTT) = cylinder_temp
              f(i,j,:,ilnrho) = f(l2,m2,n2,ilnrho) &
                *f(l2,m2,n2,ilnTT)/cylinder_temp
            endif
          endif
        enddo
      enddo
!
    endsubroutine init_solid_cells
!***********************************************************************
    
    subroutine initialize_pencils_ogrid(p,penc0)
!
!  Initialize all pencils that are necessary on ogrid
!
!  10-feb-17/Jorgen: Coded
!
      type(pencil_case_ogrid) :: p_ogrid
      real :: penc0
!
! TODO: SET ALL PENCILS
      p%XXXXX = penc0
      lpencil_ogrid(i_og_XXXXX) = .true.
!
    endsubroutine initialize_pencils_ogrid
!***********************************************************************
    subroutine dsolid_dt(f,df,p)
!
!  Dummy routine
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
!
!  Dummy routine
!
    end subroutine dsolid_dt_integrate
!***********************************************************************
    subroutine rprint_solid_cells(lreset,lwrite)
!
! TODO: Write routine
!
      logical :: lreset
      logical, optional :: lwrite
!
    end subroutine rprint_solid_cells
!***********************************************************************
    subroutine update_solid_cells(f)
!
!  Dummy routine
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    end subroutine update_solid_cells
!***********************************************************************
    subroutine update_solid_cells_pencil(f)
!
!  Dummy routine
!
      real, dimension(mx,my,mz,mfarray) :: f
!
    end subroutine update_solid_cells_pencil
!***********************************************************************
    subroutine freeze_solid_cells(df)
!
!  Dummy routine
!
      real, dimension(mx,my,mz,mvar) :: df
!
      call keep_compiler_quiet(f)
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
    subroutine flow_cartesian_to_curvlinear(f_cartesian,f_ogrid)
!
!   Interpolate all flow variables from cartesian to curvelinear grid
!
      real, dimension (mx,my,mz,mfarray) :: f_cartesian
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f_ogrid
    endsubroutine flow_cartesian_to_curvlinear(f_cartesian,f_ogrid)
!***********************************************************************
    subroutine flow_curvlinear_to_cartesian(f_cartesian,f_ogrid)
!
!   Interpolate all flow variables from curvelinear to cartesian grid
!
      real, dimension (mx,my,mz,mfarray) :: f_cartesian
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f_ogrid
    endsubroutine flow_curvlinear_to_cartesian(f_cartesian,f_ogrid)
!***********************************************************************
subroutine compute_draglift()
end subroutine compute_draglift
!***********************************************************************
subroutine print_solid()
end subroutine print_solid
!***********************************************************************
    subroutine solid_cells_timestep_first(f)
!
!  08-feb-17/Jorgen: Coded
!
!  No need to setup dfs in the beginning of each itsub.
!  Only save time, which is the time before timesteps on cartesian grid
!  is performed. Will need this to set timestep of ogrid iterations.
!
      real, dimension(mx,my,mz,mfarray) :: f

      t_ogrid = t
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
!TODO: Can probably skip this
    subroutine find_ogrid_boundaries(f)
!
!  Find the boundaries of the overlapping grid such that we can set the
!  zone of cartesian grid points inside the curvelinear grid that should
!  be set by interpolation.
!
!  Store data in arrays grid_overlap_cartesian and grid_overlap_curvlinear, 
!  for cartesian and curvlinear grid points, respectively.
!  For both arrays, the following is valid:
!  If grid_ovlap(ip,jp,kp,1)= 0  we are in a cell outside the curvelinear grid
!  If grid_ovlap(ip,jp,kp,1)= 1  we are in a cell outside the curvelinear grid
!                                that is used as a donor point for interpolation
!  If grid_ovlap(ip,jp,kp,1)= 10 we are in a cell inside the curvelinear grid
!                                that is computed in the same way as the cells
!                                outside the overlapping area
!  If grid_ovlap(ip,jp,kp,1)= 11 we are in a cell inside the curvelinear grid
!                                that is computed in the same way as the cells
!                                outside the overlapping area, and that is used
!                                as a donor point for interpolation
!  If grid_ovlap(ip,jp,kp,1)= 20 we are in a cell that should be computed by interpolation
!  If grid_ovlap(ip,jp,kp,1)= 21 we are in a cell that should no be computed at all
! 
!  09-feb-17/Jorgen: Adapted from find_solid_cell_boundaries in solid_cells.f90
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer :: i, j, k, iobj, cw
      real :: x2, y2, z2, xval_p, xval_m, yval_p, yval_m, zval_p, zval_m
      real :: dr, r_point, x_obj, y_obj, z_obj, r_obj
      character(len=10) :: form
!
!  Initialize overlapping grid array
!
      grid_ovlap_cartesian=0
      grid_ovlap_curvlinear=0

!  Set interpolated points in curvlinear grid, only boundary points interpolated
      if(x_ogrid(mx+1)>r_ogrid) then
        grid_ovlap_curvlinear(mx+1:nx)=20
        grid_ovlap_curvlinear(mx-int_stencil_len,mx)=11
      enddo
!  Set interpolation points on cartesian grid
! TODO
      r_grid2 = r_grid**2
      r_inter_outer2 = r_inter_outer**2
      r_inter_inner2 = r_inter_inner**2
      do j=m1,m2
        do i=l1,l2
          r_point2 = (x(i)-xorigo_ogrid(1))**2+(y(j)-xorigo_ogrid(2))**2
!
          if (r_point2<rgrid2) then
            grid_ovlap(i,j,:)=1
          elseif (r_point2>r_inter2)
          endif
        enddo
      enddo
  endsubroutine find_ogrid_boundaries
!***********************************************************************
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
      x0 = xyz0_ogrid(1)
      y0 = xyz0_ogrid(2)
      z0 = xyz0_ogrid(3)
      Lx = Lxyz_ogrid(1)
      Ly = Lxyz_ogrid(2)
      Lz = Lxyz_ogrid(3)
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
      use Mpicomm, only: mpiallreduce_min,mpiallreduce_max,mpibcast_real,mpirecv_real,mpisend_real
!
      real :: dxmin_x, dxmin_y, dxmin_z, dxmax_x, dxmax_y, dxmax_z
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
! Box volume, cylinder symmetrical
!
      box_volume_ogrid=1.
      if (nxgrid_ogrid/=1) then
          box_volume_ogrid = box_volume_ogrid*.5*(r_ogrid**2-xyz0_ogrid(1)**2)
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
        Area_xy_ogrid=Area_xy_ogrid*1./2.*(r_ogrid**2-xyz0_ogrid(1)**2)
        Area_xz_ogrid=Area_xz_ogrid*(r_ogrid**2-xyz0_ogrid(1)**2)
      else
        dVol_x_ogrid=1./2.*(r_ogrid**2-xyz0_ogrid(1)**2)
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
!  Set the the serial grid arrays, that contain the coordinate values
!  from all processors.
!
      call construct_serial_arrays
!
!  Construct interpolation arrays, to be used for interpolation between
!  the overlapping grids
!
      call construct_interpolation_arrays
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
!***********************************************************************
    subroutine calc_pencils_grid_ogrid(f,p)
!
!  Calculate Grid/geometry related pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   31-jan-17/Jorgen: Adapted from calc_pencils_grid in grid.f90
!                     Only cylindrical coodrinats included
!
      if (lpencil(i_og_x_mn))     p_ogrid%x_mn    = x_ogrid(l1_ogrid:l2_ogrid)*cos(x_ogrid(m_ogrid))
      if (lpencil(i_og_y_mn))     p_ogrid%y_mn    = x_ogrid(l1_ogrid:l2_ogrid)*sin(y_ogrid(m_ogrid))
      if (lpencil(i_og_z_mn))     p_ogrid%z_mn    = spread(z_ogrid(n_ogrid),1,nx_ogrid)
      if (lpencil(i_og_rcyl_mn))  p_ogrid%rcyl_mn = x_ogrid(l1_ogrid:l2_ogrid)
      if (lpencil(i_og_phi_mn))   p_ogrid%phi_mn  = spread(y_ogrid(m_ogrid),1,nx_ogrid)
      if (lpencil(i_og_rcyl_mn1)) p_ogrid%rcyl_mn1= 1./max(p_ogrid%rcyl_mn,tini)
!
    endsubroutine calc_pencils_grid_ogrid
!***********************************************************************
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
    subroutine get_grid_mn_ogrid()
!
!  Gets the geometry of the pencil at each (m,n) in the mn-loop.
!
!  31-jan-17/Jorgen: Adapted from get_grid_mn in grid.f90

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
    subroutine time_step_ogrid(f_cartesian)
!
!  Perform time steps on the curvelinear grid, including interpolation of 
!  flow variables back and forth between the overlapping grids.
!  The time iterations should equal to one time step on the cartesian grid
!
!  07-feb-17/Jorgen+Nils: Adapded from timestep.f90
!
      use Mpicomm, only: mpifinalize, mpiallreduce_max, MPI_COMM_WORLD
!
      real, dimension (mx,my,mz,mfarray) :: f_cartesian
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df_ogrid
      real :: ds, dtsub_ogrid, dt_ogrid, dt_cartesian
      integer :: timestep_factor, tstep_ogrid
      integer :: j
!!  ! TODO TODO TODO
!!  !  FROM MAIN TIME STEPPING LOOP IN RUN.F90
!!  !
!!  !  Find out which pencils to calculate at current time-step.
!!  !
!!      lpencil = lpenc_requested
!!      if (lout)   lpencil=lpencil .or. lpenc_diagnos
!!      if (l2davg) lpencil=lpencil .or. lpenc_diagnos2d
!!      if (lvideo) lpencil=lpencil .or. lpenc_video

!
!  Update boundary of the overlapping grid by interpolating
!  from cartesian grid to curvlinear grid
!
! TODO
      call flow_cartesian_to_curvlinear(f_cartesian,f_ogrid)
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
        !  use coefficients of Williamson (1980)
        alpha_ts=(/   0.0, -5/9.0 , -153/128.0 /)
        beta_ts =(/ 1/3.0, 15/16.0,    8/15.0  /)
      else
        if (lroot) print*,'Not implemented: itorder=',itorder
        call mpifinalize
      endif

!
!  Time step for ogrid
!
      timestep_factor = ceiling(dxmin/dxmin_ogrid)   ! number of timesteps on cylinder grid for one timestep on cartesian grid
      if(timestep_factor < 1)  then
        timestep_factor = 1
      endif
      dt_cartesian = t-t_ogrid
      dt_ogrid = dt_cartesian/timestep_factor
      dt_beta_ts=dt_ogrid*beta_ts
!
!  Perform a number of timestep equal to timestep_factor, such that the
!  endtime t_ogrid is equal to t after the timesteps
!
      do tstep_ogrid=1,timestep_factor

!
!  Set up df and ds for each time sub.
!

        do itsub=1,itorder
          lfirst_ogrid=(itsub==1)
          llast_ogrid=(itsub==itorder)
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
!  Calculate substep 
!
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
          t_ogrid = t_ogrid + dtsub_ogrid
!
        enddo
      enddo
!     TODO: Could perhaps need a test that checks that t_ogrid == t
!
!  Update boundary of the overlapping grid by interpolating
!  from curvlinear grid to cartesian grid
!
! TODO
      call flow_curvlinear_to_cartesian(f_cartesian,f_ogrid)
!
    endsubroutine time_step_ogrid
!***********************************************************************
    subroutine pde_ogrid(df)
!
!  06-feb-17/Jorgen+Nils: Adapded from equ.f90
!

      use Solid_cells_Mpicomm
!
!  Call the different evolution equations.
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df
      integer :: nyz
!
      intent(out)    :: df
!
!  Print statements when they are first executed.
!
      headtt = headt .and. lfirst_ogrid .and. lroot
!
      if (headtt.or.ldebug) print*,'pde_ogrid: ENTER'
      if (headtt) call svn_id("$Id$")
!
!  Prepare x-ghost zones; required before f-array communication
!  AND shock calculation
!
      call boundconds_x_ogrid(f_ogrid)
!
!  Initiate communication and do boundary conditions.
!  Required order:
!  1. x-boundaries (x-ghost zones will be communicated) - done above
!  2. communication
!  3. y- and z-boundaries
!
      if (ldebug) print*,'pde: bef. initiate_isendrcv_bdry'
      call initiate_isendrcv_bdry_ogrid(f_ogrid)
      call finalize_isendrcv_bdry_ogrid(f_ogrid)
!  Since only periodic implementation of boundconds in y- and z-dir, call only
!  if single processor in said direction. Otherwise, handled in MPI-communication.
      if (nprocy==1)                  call boundconds_y_ogrid(f_ogrid)
      if ((nprocz==1).and.(nzgrid>1)) call boundconds_z_ogrid(f_ogrid)
!
!------------------------------------------------------------------------------
!  Do loop over m and n.
!
      nyz=ny_ogrid*nz_ogrid
      mn_loop: do imn=1,nyz
        n_ogrid=nn(imn)
        m_ogrid=mm(imn)
!
!  Grid spacing. In case of equidistant grid and cartesian coordinates
!  this is calculated before the (m,n) loop.
!
        call get_grid_mn_ogrid
!
!  Calculate grid/geometry related pencils.
!
        call calc_pencils_grid_ogrid
!
!  Calculate pencils for the pencil_case.
!
        call calc_pencils_hydro_ogrid
        call calc_pencils_density_ogrid
        call calc_pencils_eos_ogrid
        call calc_pencils_viscosity_ogrid
        call calc_pencils_energy_ogrid
!
!  --------------------------------------------------------
!  NO CALLS MODIFYING PENCIL_CASE PENCILS BEYOND THIS POINT
!  --------------------------------------------------------
!
!  hydro, density, and entropy evolution
!  Note that pressure gradient is added in denergy_dt of noentropy to momentum,
!  even if lentropy=.false.
!
        call duu_dt_ogrid(df)
        call dlnrho_dt_ogrid(df)
        call denergy_dt_ogrid(df)
!
!  End of loops over m and n.
!
        headtt=.false.
      enddo mn_loop
!
    endsubroutine pde_ogrid
!***********************************************************************
    subroutine duu_dt_ogrid(df)
!
!  velocity evolution
!  calculate du/dt = - u.gradu - 2Omega x u + grav + Fvisc
!  pressure gradient force added in density and entropy modules.
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mvar) :: df
      intent(inout) :: df
!
!  Advection term.
!
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz)-p_ogrid%ugu
!
!  Viscous term
!
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p_ogrid%fvisc
!
    endsubroutine duu_dt_ogrid
!***********************************************************************
    subroutine dlnrho_dt_ogrid(df)
!
!  Continuity equation.
!  Calculate dlnrho/dt = - u.gradlnrho - divu
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mvar) :: df
      intent(inout) :: df
!
      real, dimension (nx_ogrid) :: density_rhs 
      integer :: j
!
!  Continuity equation.
!      
      density_rhs= - p_ogrid%ugrho   - p_ogrid%rho*p_ogrid%divu      
!
!  Add the continuity equation terms to the RHS of the density df.
!
      df(l1:l2,m,n,irho) = df(l1:l2,m,n,irho) + density_rhs
!
    endsubroutine dlnrho_dt_ogrid
!***********************************************************************
    subroutine denergy_dt_ogrid(df)
!
!  Calculate pressure gradient term for isothermal/polytropic equation
!  of state.
!
      real, dimension (mx_ogrid,my_ogrid,m_ogridz,mvar) :: df
      integer :: j,ju
      intent(inout) :: df
!
!  Add isothermal/polytropic pressure term in momentum equation.
!
      do j=1,3
        ju=j+iuu-1
        df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ju)=df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ju)+p_ogrid%fpres(:,j)
      enddo
!
    endsubroutine denergy_dt_ogrid
!***********************************************************************
    subroutine calc_pencils_eos_ogrid()
!
!  Calculate EquationOfState pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  10-feb-17/Jorgen+Nils: Adapted from calc_pencils_eos_pencpar in eos_idealgas.f90
!
      use eos, only: get_cp1, get_cv1
!
      real, dimension(nx_ogrid) :: tmp
      integer :: i,j
      real :: cp1, cv1
!
!  Inverse cv and cp values.
!
      call get_cp1(cp1)
      call get_cv1(cv1)
      if (lpenc_loc(i_cv1)) p_ogrid%cv1=cv1
      if (lpenc_loc(i_cp1)) p_ogrid%cp1=cp1
      if (lpenc_loc(i_cv))  p_ogrid%cv=1/cv1
      if (lpenc_loc(i_cp))  p_ogrid%cp=1/cp1
!
!  Work out thermodynamic quantities for given lnrho or rho and TT.
!
      if (iTT .gt. 0) then
        if (lpenc_loc(i_TT))   p_ogrid%TT=f_ogrid(l1:l2,m,n,iTT)
        if (lpenc_loc(i_TT1).or.lpenc_loc(i_hlnTT))  p_ogrid%TT1=1/f_ogrid(l1:l2,m,n,iTT)
        if (lpenc_loc(i_cs2))  p_ogrid%cs2=cp*gamma_m1*f_ogrid(l1:l2,m,n,iTT)
        if (lpenc_loc(i_gTT))  call grad_ogrid(f_ogrid,iTT,p_ogrid%gTT)
        if (lpenc_loc(i_del2TT).or.lpenc_loc(i_del2lnTT)) &
            call del2_ogrid(f_ogrid,iTT,p_ogrid%del2TT)
        if (lpenc_loc(i_pp)) p_ogrid%pp=cv*gamma_m1*p_ogrid%rho*p_ogrid%TT
        if (lpenc_loc(i_ee)) p_ogrid%ee=cv*exp(p_ogrid%lnTT)
!
!  Work out thermodynamic quantities for given lnrho or rho and cs2.
!
      else
        if (leos_isentropic) then
          call fatal_error('calc_pencils_eos', &
              'leos_isentropic not implemented for ilnrho_cs2, try ilnrho_ss')
        elseif (leos_isothermal) then
          if (lpenc_loc(i_cs2)) p_ogrid%cs2=cs20
          if (lpenc_loc(i_pp)) p_ogrid%pp=gamma1*p_ogrid%rho*cs20
        else
          call fatal_error('calc_pencils_eos', &
              'Full equation of state not implemented for ilnrho_cs2')
        endif
      endif
!
    endsubroutine calc_pencils_eos_ogrid
!***********************************************************************
    subroutine calc_pencils_hydro_ogrid()
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      real, dimension (nx_ogrid) :: tmp, tmp2
      integer :: i, j, ju, ij,jj,kk,jk
!
! Pencils: uu, u2, uij, divu, sij, sij2, ugu, ugu2, del2u
!
      if (lpencil(i_uu)) p_ogrid%uu=f_ogrid(l1:l2,m,n,iux:iuz)
      if (lpencil(i_u2)) call dot2_mn_ogrid(p_ogrid%uu,p_ogrid%u2)
      if (lpencil(i_uij)) call gij_ogrid(f_ogrid,iuu,p_ogrid%uij,1)
      if (lpencil(i_divu)) call div_mn_ogrid(p_ogrid%uij,p_ogrid%divu,p_ogrid%uu)
      if (lpencil(i_sij)) call traceless_strain_ogrid(p_ogrid%uij,p_ogrid%divu,p_ogrid%sij,p_ogrid%uu)
      if (lpencil(i_sij2)) call multm2_sym_mn_ogrid(p_ogrid%sij,p_ogrid%sij2)
      if (lpencil(i_ugu)) call u_dot_grad_ogrid(f_ogrid,iuu,p_ogrid%uij,p_ogrid%uu,p_ogrid%ugu)
      if (lpencil(i_ugu2)) call dot2_mn_ogrid(p_ogrid%ugu,p_ogrid%ugu2)
      if (lpencil(i_del2u)) call del2v_ogrid(f_ogrid,iux,p_ogrid%del2u)
!
    endsubroutine calc_pencils_hydro_ogrid
!***********************************************************************
    subroutine calc_pencils_density_ogrid()
!
!  Calculate Density pencils for linear density.
!  Most basic pencils should come first, as others may depend on them.
!
!  10-feb-17/Jorgen+Nils: Adapted from calc_pencils_linear_density in density.f90
!
      integer :: i
!
      if (.not. ldensity_nolog) then
        call fatal_error('calc_pencils_density_ogrid','Must use linear density for solid_cells_ogrid')
      endif
!
! Pencils: rho, rho1, lnrho, glnrho, grho, ugrho, sglnrho
!
      p_ogrid%rho=f_ogrid(l1:l2,m,n,irho)
      if (lpencil(i_rho1)) p_ogrid%rho1=1.0/p_ogrid%rho
      if (lpencil(i_lnrho))p_ogrid%lnrho=log(p_ogrid%rho)
      if (lpencil(i_glnrho).or.lpencil(i_grho)) then
        call grad_ogrid(f_ogrid,irho,p_ogrid%grho)
        if (lpencil(i_glnrho)) then
          do i=1,3
            p_ogrid%glnrho(:,i)=p_ogrid%rho1*p_ogrid%grho(:,i)
          enddo
        endif
      endif
      if (lpencil(i_ugrho)) call u_dot_grad_ogrid(f_ogrid,ilnrho,p_ogrid%grho,p_ogrid%uu,p_ogrid%ugrho)
      if (lpencil(i_sglnrho)) call multmv_mn_ogrid(p_ogrid%sij,p_ogrid%glnrho,p_ogrid%sglnrho)
!
    endsubroutine calc_pencils_density_ogrid
!***********************************************************************
    subroutine calc_pencils_viscosity_ogrid()
!
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  10-feb-17/Jorgen+Nils: Adapted from rountine in viscosity.f90
!
!  Viscous force and viscous heating are calculated here (for later use).
!
      p_ogrid%fvisc=0.0                              
!
!  viscous force: nu*del2v
!  -- not physically correct (no momentum conservation), but
!  numerically easy and in most cases qualitatively OK
!
      p_ogrid%fvisc=p_ogrid%fvisc+nu*p_ogrid%del2u
!
    end subroutine calc_pencils_viscosity_ogrid
!***********************************************************************
    subroutine calc_pencils_energy_ogrid()
!
!  Calculate Energy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      integer :: j
!
!  Pencils: fpres (=pressure gradient force)
!
      do j=1,3
        p_ogrid%fpres(:,j)=-p_ogrid%cs2*p_ogrid%glnrho(:,j)
      enddo
!
    endsubroutine calc_pencils_energy_ogrid
!***********************************************************************
!  FROM BOUNDCOND.F90
!***********************************************************************
!  ROUTINES
!    boundconds_x_ogrid
!    boundconds_y_ogrid
!    boundconds_z_ogrid
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
            !call special_boundconds(f,bc)
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
    subroutine boundconds_y_ogrid(f)
!
!  Periodic boundary condition for runs with a single processor in y-direction 
!
!  06-feb-17/Jorgen: Adapted from boundcond.f90 to be used for ogrids where the
!                    y-dir is always periodic
!
      real, dimension (:,:,:,:) :: f
!
      integer :: ivar1, ivar2, j
      integer :: m1i_ogrid=m1_ogrid+nghost-1
      integer :: m2i_ogrid=my_ogrid-2*nghost+1
!
      ivar1=1; ivar2=min(mcom,size(f,4))
!
!  Boundary conditions in y
!  Periodic, with y being the theta direction for the cylindrical grid
!
      do j=ivar1,ivar2
!  Bottom boundary
        f(:,1:m1_ogrid-1,:,j) = f(:,m2i_ogrid:m2_ogrid,:,j)
!  Top boundary
        f(:,m2_ogrid+1:,:,j) = f(:,m1_ogrid:m1i_ogrid,:,j)
      enddo
!
    endsubroutine boundconds_y_ogrid
!***********************************************************************
    subroutine boundconds_z_ogrid(f)
!
!  Periodic boundary condition for 3D-runs with a single processor in z-direction 
!
!  06-feb-17/Jorgen: Adapted from boundcond.f90 to be used for ogrids where the
!                    z-dir is always periodic as long as nzgrid=/1
!
      real, dimension (:,:,:,:) :: f
      integer :: ivar1, ivar2, j

      integer :: n1i_ogrid=n1_ogrid+nghost-1
      integer :: n2i_ogrid=mz_ogrid-2*nghost+1
!
      if (ldebug) print*,'boundconds_z: ENTER: boundconds_z'
!
      ivar1=1; ivar2=min(mcom,size(f,4))
!
!  Boundary conditions in z
!
      do j=ivar1,ivar2
!  Bottom boundary
        f(:,:,1:n1_ogrid-1,j) = f(:,:,n2i_ogrid:n2_ogrid,j)
!  Top boundary
        f(:,:,n2_ogrid+1:,j) = f(:,:,n1_ogrid:n1i_ogrid,j)
      enddo
!
    endsubroutine boundconds_z_ogrid
!***********************************************************************
! FROM SUB.f90
!***********************************************************************
! ROUTINES
!   grad_ogrid
!   del2_ogrid
!   g2ij
!   dot2_mn_ogrid
!   gij_ogrid
!   div_mn_ogrid
!   traceless_strain_ogrid
!   multm2_sym_mn_ogrid
!   curl_mn_ogrid
!   dot_mn_ogrid
!   dot2_0_ogrid
!   del2v_ogrid
!   u_dot_grad_vec_ogrid
!   u_dot_grad_scl_ogrid
!   multmv_mn_ogrid
!***********************************************************************
    subroutine grad_ogrid(f,k,g)
!
!  Calculate gradient of a scalar, get vector.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid,3) :: g
      integer :: k
!
      intent(in) :: f,k
      intent(out) :: g
!
      call der_ogrid(f,k,g(:,1),1)
      call der_ogrid(f,k,g(:,2),2)
      call der_ogrid(f,k,g(:,3),3)
!
    endsubroutine grad_ogrid
!***********************************************************************
    subroutine del2_ogrid(f,k,del2f)
!
!  Calculate del2 of a scalar, get scalar.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      intent(in) :: f,k
      intent(out) :: del2f
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid) :: del2f,d2fdx,d2fdy,d2fdz,tmp
      integer :: k
!
      call der2_ogrid(f,k,d2fdx,1)
      call der2_ogrid(f,k,d2fdy,2)
      call der2_ogrid(f,k,d2fdz,3)
      del2f=d2fdx+d2fdy+d2fdz
!
!  Since we have cylindrical coordinates
      call der_ogrid(f,k,tmp,1)
      del2f=del2f+tmp*rcyl_mn1_ogrid
!
    endsubroutine del2_ogrid
!***********************************************************************
    subroutine g2ij(f,k,g)
!
!  Calculates the Hessian, i.e. all second derivatives of a scalar.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid,3,3) :: g
      real, dimension (nx_ogrid) :: tmp
      integer :: i,j,k
!
      intent(in) :: f,k
      intent(out) :: g
!
!  Run though all 9 possibilities, treat diagonals separately.
!
      do j=1,3
        call der2_ogrid(f,k,tmp,j)
        g(:,j,j)=tmp
        do i=j+1,3
          call derij(f,k,tmp,i,j)
          g(:,i,j)=tmp
          g(:,j,i)=tmp
        enddo
      enddo
!
    endsubroutine g2ij
!***********************************************************************
    subroutine dot2_mn_ogrid(a,b)
!
!  Dot product with itself, to calculate max and rms values of a vector.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (nx_ogrid,3) :: a
      real, dimension (nx_ogrid) :: b
!
      intent(in) :: a
      intent(out) :: b

      b=a(:,1)**2+a(:,2)**2+a(:,3)**2
!
    endsubroutine dot2_mn_ogrid
!***********************************************************************
    subroutine gij_ogrid(f,k,g,nder)
!
!  Calculate gradient of a vector, return matrix.
!
!   07-feb-17/Jorgen: Adapted from sub.f90
!
      use Deriv, only: der,der2,der3,der4,der5,der6
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid,3,3) :: g
      real, dimension (nx_ogrid) :: tmp
      integer :: i,j,k,k1,nder
!
      intent(in) :: f,k,nder
      intent(out) :: g
!
      k1=k-1
      do i=1,3 
        do j=1,3
          if (nder == 1) then
            call der_ogrid(f,k1+i,tmp,j)
          elseif (nder == 2) then
            call der2_ogrid(f,k1+i,tmp,j)
          endif
          g(:,i,j)=tmp
        enddo 
      enddo
!
    endsubroutine gij_ogrid
!***********************************************************************
    subroutine div_mn_ogrid(aij,b,a)
!
!  Calculate divergence from derivative matrix.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (nx_ogrid,3,3) :: aij
      real, dimension (nx_ogrid,3) :: a
      real, dimension (nx_ogrid) :: b
!
      intent(in) :: aij,a
      intent(out) :: b
!
      b=aij(:,1,1)+aij(:,2,2)+aij(:,3,3)
!
!  Adjustments for other cylindrical coordinate system
!
      b=b+rcyl_mn1_ogrid*a(:,1)
!
    endsubroutine div_mn_ogrid
!***********************************************************************
    subroutine traceless_strain_ogrid(uij,divu,sij,uu)
!
!  Calculates traceless rate-of-strain tensor sij from derivative tensor uij
!  and divergence divu within each pencil;
!  curvilinear co-ordinates require velocity argument uu, so this is not optional here
!  In-place operation is possible, i.e. uij and sij may refer to the same array.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
    use General, only: loptest

    real, dimension (nx_ogrid,3,3)   :: uij, sij
    real, dimension (nx_ogrid)       :: divu
    real, dimension (nx_ogrid,3)     :: uu
!
    integer :: i,j
!
    intent(in)  :: uij, divu, lss
    intent(out) :: sij
!
    do j=1,3
      sij(:,j,j)=uij(:,j,j)-(1./3.)*divu
      do i=j+1,3
        sij(:,i,j)=.5*(uij(:,i,j)+uij(:,j,i))
        sij(:,j,i)=sij(:,i,j)
      enddo
    enddo
    sij(:,1,2)=sij(:,1,2)-.5*rcyl_mn1_ogrid*uu(:,2)
    sij(:,2,2)=sij(:,2,2)+.5*rcyl_mn1_ogrid*uu(:,1)
    sij(:,2,1)=sij(:,1,2)
!
    endsubroutine traceless_strain_ogrid
!***********************************************************************
    subroutine multm2_sym_mn_ogrid(a,b)
!
!  Symmetric matrix squared, gives scalar.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (nx_ogrid,3,3), intent(in) :: a
      real, dimension (nx_ogrid), intent(out) :: b
!
      integer :: i, j
!
      b = a(:,1,1)**2
      do i = 2, 3
        b = b + a(:,i,i)**2
        do j = 1, i-1
          b = b + 2 * a(:,i,j)**2
        enddo
      enddo
!
    endsubroutine multm2_sym_mn_ogrid
!***********************************************************************
    subroutine curl_mn_ogrid(aij,b,a)
!
!  Calculate curl from derivative matrix.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (nx,3,3), intent (in) :: aij
      real, dimension (nx,3), intent (in), optional :: a
      real, dimension (nx,3), intent (out) :: b
      integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7
!
      b(:,1)=aij(:,3,2)-aij(:,2,3)
      b(:,2)=aij(:,1,3)-aij(:,3,1)
      b(:,3)=aij(:,2,1)-aij(:,1,2)
!
!  Adjustments for cylindrical coordinate system.
!  If we go all the way to the center, we need to put a regularity condition.
!  We do this here currently up to second order, and only for curl_mn.
!
      b(:,3)=b(:,3)+a(:,2)*rcyl_mn1_ogrid
      if (rcyl_mn_ogrid(1)==0.) b(i1,3)=(360.*b(i2,3)-450.*b(i3,3)+400.*b(i4,3) &
                                  -225.*b(i5,3)+72.*b(i6,3)-10.*b(i7,3))/147.
!
    endsubroutine curl_mn_ogrid
!***********************************************************************
    subroutine dot_mn_ogrid(a,b,c)
!
!  Dot product, c=a.b, on pencil arrays
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (:,:) :: a,b
      real, dimension (:) :: c
!
      intent(in) :: a,b
      intent(inout) :: c
!
      integer :: i
      logical :: l0
!
      c=0.
      do i=1,size(a,2)
        c=c+a(:,i)*b(:,i)
      enddo
!
    endsubroutine dot_mn_ogrid
!***********************************************************************
    subroutine dot2_0_ogrid(a,b)
!
!  Dot product, c=a.b, of two simple 3-d arrays.
!
!  07-fab-15/Jorgen: Copied from sub.f90
!
      real, dimension (:) :: a
      real :: b
!
      intent(in) :: a
      intent(out) :: b
!
      b = dot_product(a,a)
!
    endsubroutine dot2_0_ogrid
!***********************************************************************
    subroutine del2v_ogrid(f,k,del2f)
!
!  Calculate del2 of a vector, get vector.
!
!  28-oct-97/axel: coded
!  15-mar-07/wlad: added cylindrical coordinates
!
      use Deriv, only: der
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid,3) :: del2f
      real, dimension (nx_ogrid) :: tmp
      integer :: i,k,k1
!
      intent(in) :: f,k
      intent(out) :: del2f
!
!  do the del2 diffusion operator
!
      k1=k-1
      do i=1,3
        call del2(f,k1+i,tmp)
        del2f(:,i)=tmp
      enddo
!
      !del2 already contains the extra term 1/r*d(uk)/dt
      call der(f,k1+2,tmp,2)
      del2f(:,1)=del2f(:,1) -(2*tmp+f(l1:l2,m,n,k1+1))*rcyl_mn2_ogrid
      call der(f,k1+1,tmp,2)
      del2f(:,2)=del2f(:,2) +(2*tmp-f(l1:l2,m,n,k1+2))*rcyl_mn2_ogrid
!
    endsubroutine del2v_ogrid
!***********************************************************************
    subroutine u_dot_grad_vec_ogrid(f,k,gradf,uu,ugradf)
!
!  u.gradu
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      intent(in) :: f,k,gradf,uu
      intent(out) :: ugradf
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid,3,3) :: gradf
      real, dimension (nx_ogrid,3) :: uu,ff,ugradf,grad_f_tmp
      real, dimension (nx_ogrid) :: tmp
      integer :: j,k
!
      if (k<1 .or. k>mfarray) then
        call fatal_error('u_dot_grad_vec_ogrid','variable index is out of bounds')
        return
      endif
!
      do j=1,3
        grad_f_tmp = gradf(:,j,:)
        call u_dot_grad_scl(f,k+j-1,grad_f_tmp,uu,tmp)
        ugradf(:,j)=tmp
      enddo
!
!  Adjust to cylindrical coordinate system
!  The following now works for general u.gradA.
!
      ff=f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,k:k+2)
      ugradf(:,1)=ugradf(:,1)-rcyl_mn1_ogrid*(uu(:,2)*ff(:,2))
      ugradf(:,2)=ugradf(:,2)+rcyl_mn1_ogrid*(uu(:,1)*ff(:,2))
!
    endsubroutine u_dot_grad_vec_ogrid
!***********************************************************************
    subroutine u_dot_grad_scl_ogrid(f,k,gradf,uu,ugradf)
!
!  Do advection-type term u.grad f_k.
!  Assumes gradf to be known, but takes f and k as arguments to be able
!  to calculate upwind correction.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      use General, only: loptest
!
      intent(in) :: f,k,gradf,uu
      intent(out) :: ugradf
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid,3) :: uu,gradf
      real, dimension (nx_ogrid) :: ugradf
      integer :: k
!
      if (k<1 .or. k>mfarray) then
        call fatal_error('u_dot_grad_scl_ogrid','variable index is out of bounds')
        return
      endif
!
      call dot_mn_ogrid(uu,gradf,ugradf)
!
    endsubroutine u_dot_grad_scl_ogrid
!***********************************************************************
    subroutine multmv_mn_ogrid(a,b,c)
!
!  Matrix multiplied with vector, gives vector.
!
!  C_i = A_{i,j} B_j
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (nx_ogrid,3,3) :: a
      real, dimension (nx_ogrid,3) :: b,c
      real, dimension (nx_ogrid) :: tmp
      integer :: i,j
      logical, optional :: ladd
!
      intent(in) :: a,b,ladd
      intent(out) :: c
!
      do i=1,3
        j=1
        tmp=a(:,i,j)*b(:,j)
        do j=2,3
          tmp=tmp+a(:,i,j)*b(:,j)
        enddo
        c(:,i)=tmp
      enddo
!
    endsubroutine multmv_mn_ogrid
!***********************************************************************
! FROM DERIV.F90
!***********************************************************************
! ROUTINES
!   der_ogrid
!   der2_ogrid
!   derij_ogrid
!***********************************************************************
    subroutine der_ogrid(f, k, df, j, ignoredx)
!
!  calculate derivative df_k/dx_j 
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!  07-feb-17/Jorgen: Adapted from deriv.f90
!
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray), intent(in) :: f
      real, dimension(nx_ogrid), intent(out) :: df
      integer, intent(in) :: j, k
      logical, intent(in), optional :: ignoredx
!
      real, parameter :: a = 1.0 / 60.0
      real, dimension(nx_ogrid) :: fac
      logical :: withdx
!
      if (present(ignoredx)) then
        withdx = .not. ignoredx
        if (ignoredx) fac = a
      else
        withdx = .true.
      endif
!
      if (j==1) then
        if (nxgrid_ogrid/=1) then
          if (withdx) fac = a * dx_1_ogrid(l1_ogrid:l2_ogrid)
          df=fac*(+ 45.0*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid,k)) &
                  -  9.0*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid,k)) &
                  +      (f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid,k)))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (withdx) fac = a * dy_1_ogrid(m_ogrid)
          df=fac*(+ 45.0*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid,k)) &
                  -  9.0*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid,k)) &
                  +      (f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid,k)))
          ! Since we have cylindrical coordinates
          if (withdx) df = df * rcyl_mn1_ogrid
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (withdx) fac = a * dz_1_ogrid(n_ogrid)
          df=fac*(+ 45.0*(f(l1:l2_ogrid,m_ogrid,n_ogrid+1,k)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-1,k)) &
                  -  9.0*(f(l1:l2_ogrid,m_ogrid,n_ogrid+2,k)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-2,k)) &
                  +      (f(l1:l2_ogrid,m_ogrid,n_ogrid+3,k)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-3,k)))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_ogrid
!***********************************************************************
    subroutine der2_ogrid(f,k,df2,j)
!
!  calculate 2nd derivative d^2f_k/dx_j^2
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!  07-feb-17/Jorgen: Adapted from deriv.f90
!
      use General, only: loptest

      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid) :: df2,fac,df
      integer :: j,k
!
      real :: der2_coef0, der2_coef1, der2_coef2, der2_coef3

      intent(in)  :: f,k,j
      intent(out) :: df2
!
      der2_coef0=-490./180.; der2_coef1=270./180.
      der2_coef2=-27./180.; der2_coef3=2./180.

      if (j==1) then
        if (nxgrid_ogrid/=1) then
          fac=dx_1_ogrid(l1_ogrid:l2_ogrid)**2
          df2=fac*(der2_coef0* f(l1_ogrid  :l2_ogrid  ,m_ogrid,n_ogrid,k) &
                  +der2_coef1*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid,k)+f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid,k)) &
                  +der2_coef2*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid,k)+f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid,k)) &
                  +der2_coef3*(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid,k)+f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid,k)))
          if (.not.lequidist_ogrid(j)) then
            call der_ogrid(f,k,df,j)
            df2=df2+dx_tilde_ogrid(l1_ogrid:l2_ogrid)*df
          endif
        else
          df2=0.
        endif
      elseif (j==2) then
        if (nygrid_ogrid/=1) then
          fac=dy_1_ogrid(m_ogrid)**2
          df2=fac*(der2_coef0* f(l1_ogrid:l2_ogrid,m_ogrid  ,n_ogrid,k) &
                  +der2_coef1*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid,k)) &
                  +der2_coef2*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid,k)) &
                  +der2_coef3*(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid,k)))
          df2=df2*rcyl_mn2_ogrid
          if (.not.lequidist_ogrid(j)) then
            call der(f,k,df,j)
            df2=df2+dy_tilde_ogrid(m_ogrid)*df
          endif
        else
          df2=0.
        endif
      elseif (j==3) then
        if (nzgrid_ogrid/=1) then
          fac=dz_1_ogrid(n_ogrid)**2
          df2=fac*( der2_coef0* f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid    ,k) &
                   +der2_coef1*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+1,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-1,k)) &
                   +der2_coef2*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+2,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-2,k)) &
                   +der2_coef3*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+3,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-3,k)))
          if (.not.lequidist_ogrid(j)) then
            call der(f,k,df,j)
            df2=df2+dz_tilde_ogrid(n_ogrid)*df
          endif
        else
          df2=0.
        endif
      endif
!
    endsubroutine der2_ogrid
!***********************************************************************
    subroutine derij_ogrid(f,k,df,i,j)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!
!  17-feb-17/Jorgen: Adapted from deriv.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid) :: df,fac
      integer :: i,j,k
!
      intent(in) :: f,k,i,j
      intent(out) :: df
!
      if (lbidiagonal_derij) then
        !
        ! Use bidiagonal mixed-derivative operator, i.e.
        ! employ only the three neighbouring points on each of the four
        ! half-diagonals. This gives 6th-order mixed derivatives as the
        ! version below, but involves just 12 points instead of 36.
        !
        if ((i==1.and.j==2).or.(i==2.and.j==1)) then
          if (nxgrid_ogrid/=1.and.nygrid_ogrid/=1) then
            fac=(1./720.)*dx_1_ogrid(l1_ogrid:l2_ogrid)*dy_1_ogrid(m_ogrid)
            df=fac*( &
                        270.*( f(l1_ogrid+1:l2_ogrid+1,m_ogrid+1,n_ogrid,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid+1,n_ogrid,k)  &
                              +f(l1_ogrid-1:l2_ogrid-1,m_ogrid-1,n_ogrid,k)-f(l1_ogrid+1:l2_ogrid+1,m_ogrid-1,n_ogrid,k)) &
                       - 27.*( f(l1_ogrid+2:l2_ogrid+2,m_ogrid+2,n_ogrid,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid+2,n_ogrid,k)  &
                              +f(l1_ogrid-2:l2_ogrid-2,m_ogrid-2,n_ogrid,k)-f(l1_ogrid+2:l2_ogrid+2,m_ogrid-2,n_ogrid,k)) &
                       +  2.*( f(l1_ogrid+3:l2_ogrid+3,m_ogrid+3,n_ogrid,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid+3,n_ogrid,k)  &
                              +f(l1_ogrid-3:l2_ogrid-3,m_ogrid-3,n_ogrid,k)-f(l1_ogrid+3:l2_ogrid+3,m_ogrid-3,n_ogrid,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid_ogrid/=1.and.nzgrid_ogrid/=1) then
            fac=(1./720.)*dy_1_ogrid(m_ogrid)*dz_1_ogrid(n_ogrid)
            df=fac*( &
                        270.*( f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid+1,k)-f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid-1,k)  &
                              +f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid-1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid+1,k)) &
                       - 27.*( f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid+2,k)-f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid-2,k)  &
                              +f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid-2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid+2,k)) &
                       +  2.*( f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid+3,k)-f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid-3,k)  &
                              +f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid-3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid+3,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid_ogrid/=1.and.nxgrid_ogrid/=1) then
            fac=(1./720.)*dz_1_ogrid(n_ogrid)*dx_1_ogrid(l1_ogrid:l2_ogrid)
            df=fac*( &
                        270.*( f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid+1,k)  &
                              +f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid-1,k)-f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid-1,k)) &
                       - 27.*( f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid+2,k)  &
                              +f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid-2,k)-f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid-2,k)) &
                       +  2.*( f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid+3,k)  &
                              +f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid-3,k)-f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid-3,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or z-direction'
          endif
        endif
!
      else                      ! not using bidiagonal mixed derivatives
        !
        ! This is the old, straight-forward scheme
        !
        if ((i==1.and.j==2).or.(i==2.and.j==1)) then
          if (nxgrid_ogrid/=1.and.nygrid_ogrid/=1) then
            fac=(1./60.**2)*dx_1_ogrid(l1_ogrid:l2_ogrid)*dy_1_ogrid(m_ogrid)
            df=fac*( &
              45.*((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid+1,n_ogrid,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid+1,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid+1,n_ogrid,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid+1,n_ogrid,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid+1,n_ogrid,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid+1,n_ogrid,k))) &
                  -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid-1,n_ogrid,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid-1,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid-1,n_ogrid,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid-1,n_ogrid,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid-1,n_ogrid,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid-1,n_ogrid,k))))&
              -9.*((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid+2,n_ogrid,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid+2,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid+2,n_ogrid,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid+2,n_ogrid,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid+2,n_ogrid,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid+2,n_ogrid,k))) &
                  -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid-2,n_ogrid,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid-2,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid-2,n_ogrid,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid-2,n_ogrid,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid-2,n_ogrid,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid-2,n_ogrid,k))))&
                 +((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid+3,n_ogrid,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid+3,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid+3,n_ogrid,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid+3,n_ogrid,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid+3,n_ogrid,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid+3,n_ogrid,k))) &
                  -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid-3,n_ogrid,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid-3,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid-3,n_ogrid,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid-3,n_ogrid,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid-3,n_ogrid,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid-3,n_ogrid,k))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid_ogrid/=1.and.nzgrid_ogrid/=1) then
            fac=(1./60.**2)*dy_1_ogrid(m_ogrid)*dz_1_ogrid(n_ogrid)
            df=fac*( &
              45.*((45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid+1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid+1,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid+1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid+1,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid+1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid+1,k))) &
                  -(45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid-1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid-1,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid-1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid-1,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid-1,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid-1,k))))&
              -9.*((45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid+2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid+2,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid+2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid+2,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid+2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid+2,k))) &
                  -(45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid-2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid-2,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid-2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid-2,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid-2,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid-2,k))))&
                 +((45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid+3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid+3,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid+3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid+3,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid+3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid+3,k))) &
                  -(45.*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid-3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid-3,k))  &
                    -9.*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid-3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid-3,k))  &
                       +(f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid-3,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid-3,k))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid_ogrid/=1.and.nxgrid_ogrid/=1) then
            fac=(1./60.**2)*dz_1_ogrid(n_ogrid)*dx_1_ogrid(l1_ogrid:l2_ogrid)
            df=fac*( &
              45.*((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid+1,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid+1,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid+1,k))) &
                  -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid-1,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid-1,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid-1,k))))&
              -9.*((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid+2,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid+2,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid+2,k))) &
                  -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid-2,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid-2,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid-2,k))))&
                 +((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid+3,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid+3,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid+3,k))) &
                  -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid-3,k))  &
                    -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid-3,k))  &
                       +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid-3,k))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or z-direction'
          endif
        endif
!
      endif                     ! bidiagonal derij

!  Since we have cylindrical coordinates
      if ( i+j==3 .or. i+j==5 ) df=df*rcyl_mn1_ogrid
!
    endsubroutine derij_ogrid
!***********************************************************************
end module Solid_Cells
