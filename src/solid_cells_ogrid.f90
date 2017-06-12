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

!  Cylinder parameters
  logical :: lset_flow_dir=.false.
  real :: cylinder_radius=0.                        ! Set in start.in
  real :: cylinder_temp=703.0                       ! Set in start.in
  real :: cylinder_xpos=0., cylinder_ypos=0.        ! Set in start.in
  real :: cylinder_zpos=0.                          ! Set in start.in
  real :: skin_depth=0.                             ! Set in start.in
  real :: init_uu=0., ampl_noise=0.
  character(len=labellen) :: initsolid_cells='cylinderstream_x'! Set in start.in
  real :: T0 ! Inlet temperature
  integer :: ncylinders=1,flow_dir=0, flow_dir_set
  real, dimension(3) :: xyz0_ogrid, Lxyz_ogrid, xorigo_ogrid
!  Boundary condition
  logical :: SBP=.true.
!  Fundamental grid parameters
  real :: r_ogrid=0.                                                 ! Set in start.in?
  character (len=labellen), dimension(3) :: grid_func_ogrid='linear' ! Set in start.in
  integer :: inter_stencil_len = 1                                ! set in start.in?
  integer, parameter :: nx_ogrid=nxgrid_ogrid/nprocx,ny_ogrid=nygrid_ogrid/nprocy,nz_ogrid=nzgrid_ogrid/nprocz
!***************************************************
! Pencil case ogrid
  integer, parameter :: npencils_ogrid=29
  type pencil_case_ogrid
    real, dimension (nx_ogrid)     :: x_mn    
    real, dimension (nx_ogrid)     :: y_mn    
    real, dimension (nx_ogrid)     :: z_mn    
    real, dimension (nx_ogrid)     :: rcyl_mn 
    real, dimension (nx_ogrid)     :: phi_mn  
    real, dimension (nx_ogrid)     :: rcyl_mn1
    real, dimension (nx_ogrid,3)   :: fpres
    real, dimension (nx_ogrid,3)   :: fvisc
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
    real, dimension (nx_ogrid,3)   :: graddivu
    real, dimension (nx_ogrid)     :: lnTT
    real, dimension (nx_ogrid)     :: cs2
    real, dimension (nx_ogrid)     :: pp
    real, dimension (nx_ogrid)     :: ss
  endtype pencil_case_ogrid
!  
  integer :: i_og_x_mn    =1
  integer :: i_og_y_mn    =2
  integer :: i_og_z_mn    =3
  integer :: i_og_rcyl_mn =4
  integer :: i_og_phi_mn  =5
  integer :: i_og_rcyl_mn1=6
  integer :: i_og_fpres   =7
  integer :: i_og_fvisc   =8
  integer :: i_og_rho     =9
  integer :: i_og_rho1    =10
  integer :: i_og_lnrho   =11
  integer :: i_og_grho    =12
  integer :: i_og_glnrho  =13
  integer :: i_og_ugrho   =14
  integer :: i_og_sglnrho =15
  integer :: i_og_uu      =16
  integer :: i_og_u2      =17
  integer :: i_og_uij     =18
  integer :: i_og_divu    =19
  integer :: i_og_sij     =20
  integer :: i_og_sij2    =21
  integer :: i_og_ugu     =22
  integer :: i_og_ugu2    =23
  integer :: i_og_del2u   =24
  integer :: i_og_graddivu=25
  integer :: i_og_lnTT    =26
  integer :: i_og_cs2     =27
  integer :: i_og_pp      =28
  integer :: i_og_ss      =29
!
  character (len=15), parameter, dimension(npencils_ogrid) :: pencil_names_ogrid = &
    (/ 'x_mn          ', 'y_mn          ', 'z_mn          ', 'rcyl_mn       '  &
     , 'phi_mn        ', 'rcyl_mn1      ', 'fpres         ', 'fvisc         '  &
     , 'rho           '  &
     , 'rho1          ', 'lnrho         ', 'grho          ', 'glnrho        '  &
     , 'ugrho         ', 'sglnrho       '  &
     , 'uu            ', 'u2            ', 'uij           ', 'divu          '  &
     , 'sij           ', 'sij2          ', 'ugu           ', 'ugu2          '  &
     , 'del2u         ', 'graddivu      '  &
     , 'lnTT          '  &
     , 'cs2           '  &
     , 'pp            ', 'ss            ' /)
  logical,dimension(npencils_ogrid):: lpencil_ogrid
!***************************************************
! DATA TYPE FOR INTERPOLATION AND INTERPOLATION STENCILS
!***************************************************
  type interpol_comm_metadata
    integer :: send_to_proc
    integer :: ip_id
    integer, dimension(3) :: i_low_corner
  endtype interpol_comm_metadata
!***************************************************
  type interpol_grid_data
    integer, dimension(3) :: i_xyz                    ! Indices on the grid we interpolate TO
    integer, dimension(3) :: ind_global_neighbour     ! Indices on the grid we interpolate FROM (global)
    integer, dimension(3) :: ind_local_neighbour      ! Indices on the grid we interpolate FROM (local)
    real, dimension(3)    :: xyz                      ! Coordinates on the grid we interpolate FROM
    integer :: from_proc                              ! Get data from other processor
  endtype interpol_grid_data
!***************************************************
  type(interpol_grid_data), dimension(:), allocatable :: curvilinear_to_cartesian
  type(interpol_grid_data), dimension(:), allocatable :: cartesian_to_curvilinear
  type(interpol_comm_metadata), dimension(:), allocatable :: send_curvilinear_to_cartesian
  type(interpol_comm_metadata), dimension(:), allocatable :: send_cartesian_to_curvilinear
!***************************************************
  character(len=9) :: grid_interpolation='trilinear'
  real, dimension(ncpus,3) :: xyz0_loc_all            ! Grid start, cartesian grid, all processors
  real, dimension(ncpus,3) :: xyz1_loc_all            ! Grid stop, cartesian grid, all processors
  real, dimension(3) :: xyz0_loc_ogrid                ! Grid start, curviliner grid, local
  real, dimension(3) :: xyz1_loc_ogrid                ! Grid stop, curviliner grid, local
  real, dimension(ncpus,3) :: xyz0_loc_all_ogrid      ! Grid start, curviliner grid, all processors
  real, dimension(ncpus,3) :: xyz1_loc_all_ogrid      ! Grid stop, curviliner grid, all processors
  integer :: n_ip_cart_to_curv=0                      ! Number of interpolations performed by this proc, from curvilinear to cartesian
  integer :: n_ip_curv_to_cart=0                      ! Number of interpolations performed by this proc, from curvilinear to cartesian
  logical :: lcheck_interpolation=.true.              ! Should we check that all interpolated points are withing data points?
  logical :: lcheck_init_interpolation=.true.         ! Should we check that pre-processed intepolation stencils are set up correctly?
  integer :: n_procs_send_curv_to_cart                ! Number of processors to send ip to
  integer :: n_procs_send_cart_to_curv                ! Number of processors to send ip to
  integer, dimension(:), allocatable :: n_ip_to_proc_curv_to_cart   ! Number of ip to send to each of these processors
  integer, dimension(:), allocatable :: n_ip_to_proc_cart_to_curv   ! Number of ip to send to each of these processors
  integer :: n_procs_recv_curv_to_cart                ! Number of processors to recieve ip from
  integer :: n_procs_recv_cart_to_curv                ! Number of processors to recieve ip from
  integer, dimension(:), allocatable :: n_ip_recv_proc_curv_to_cart ! Number of ip to recieve from each of these processors
  integer, dimension(:), allocatable :: n_ip_recv_proc_cart_to_curv ! Number of ip to recieve from each of these processors
  integer, dimension(:), allocatable :: procs_recv_curv_to_cart     ! Ordered array of processors to recieve ip from 
  integer, dimension(:), allocatable :: procs_recv_cart_to_curv     ! Ordered array of processors to recieve ip from 
  integer, dimension(:), allocatable :: procs_send_curv_to_cart     ! Ordered array of processors to send ip to
  integer, dimension(:), allocatable :: procs_send_cart_to_curv     ! Ordered array of processors to send ip to
  !!
  integer :: max_send_ip_curv_to_cart                 ! Maximum number of ip sent
  integer :: max_recv_ip_curv_to_cart                 ! Maximum number of ip recv
  integer :: max_send_ip_cart_to_curv                 ! Maximum number of ip sent
  integer :: max_recv_ip_cart_to_curv                 ! Maximum number of ip recv
  !
  real :: r_int_outer, r_int_inner
!***************************************************
! PARAMETERS NECESSARY FOR GRID CONSTRUCTION 
!  Global ogrid
!  A few of these are necessary for the interpolation between grids
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
  real, dimension(3) :: xyz_star_ogrid=0.0
  real, dimension(-nghost:nghost,2) :: coeffs_1_x_ogrid
  real, dimension(-nghost:nghost,2) :: coeffs_1_y_ogrid
  real, dimension(-nghost:nghost,2) :: coeffs_1_z_ogrid
  real, dimension(3) :: coeff_grid_o=1.0
  real, dimension (-nghost:nghost) :: dx2_bound_ogrid=0., dy2_bound_ogrid=0., dz2_bound_ogrid=0.
!  Grid properties computed in grid routines, but not used in current implementation
  real, dimension (mx_ogrid) :: dVol_x_ogrid,dVol1_x_ogrid
  real, dimension (my_ogrid) :: dVol_y_ogrid,dVol1_y_ogrid
  real, dimension (mz_ogrid) :: dVol_z_ogrid,dVol1_z_ogrid
  real, dimension (nx_ogrid) :: dxyz_2_ogrid,dxyz_4_ogrid,dxyz_6_ogrid,dVol_ogrid
  real, dimension(0:nprocx) :: procx_bounds_ogrid
  real, dimension(0:nprocy) :: procy_bounds_ogrid
  real, dimension(0:nprocz) :: procz_bounds_ogrid
  real :: box_volume_ogrid=1.0 
  integer :: lpoint_ogrid=(mx_ogrid+1)/2, mpoint_ogrid=(my_ogrid+1)/2, npoint_ogrid=(mz_ogrid+1)/2
  integer :: lpoint2_ogrid=(mx_ogrid+1)/4,mpoint2_ogrid=(my_ogrid+1)/4,npoint2_ogrid=(mz_ogrid+1)/4

  real, dimension (nx_ogrid) :: dxmax_pencil_ogrid,dxmin_pencil_ogrid  ! Possible remove if no heat flux or viscosity module
  
  ! For mn-loop over ogrid pencils
  integer :: imn_ogrid
  integer, target :: m_ogrid,n_ogrid
  integer, dimension (ny_ogrid*nz_ogrid) :: mm_ogrid,nn_ogrid
  logical, dimension (ny_ogrid*nz_ogrid) :: necessary_ogrid=.false.
  integer, dimension (my_ogrid,mz_ogrid) :: imn_array_ogrid
  ! For time-step
  logical :: llast_ogrid

!  Pencils and f-array to be used for curvilinear grid computations
  type(pencil_case_ogrid) p_ogrid 
  save p_ogrid
  real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray), save ::  f_ogrid=0.

!  Summation by parts arrays
  real, dimension(6,9) :: D1_SBP, D2_SBP

!  EOS parameters
  real :: rho0, lnrho0, lnTT0

!  Parameters related to interpolation and communication of interpolation points
  real, dimension(:,:), allocatable :: xyz_comm
  real, dimension(:,:), allocatable :: rthz_comm
!  Diagnostics for output
  integer :: idiag_c_dragx=0
  integer :: idiag_c_dragy=0
!
!  Read start.in file
  namelist /solid_cells_init_pars/ &
      cylinder_temp, cylinder_radius, cylinder_xpos, ncylinders, &
      cylinder_ypos, cylinder_zpos, flow_dir_set, skin_depth, &
      initsolid_cells, init_uu, r_ogrid, lset_flow_dir,ampl_noise, &
      grid_func_ogrid, coeff_grid_o, xyz_star_ogrid, grid_interpolation, &
      lcheck_interpolation, lcheck_init_interpolation, SBP, inter_stencil_len

!  Read run.in file
  namelist /solid_cells_run_pars/ &
      flow_dir_set, lset_flow_dir
    
  interface dot2_ogrid
    module procedure dot2_mn_ogrid
    module procedure dot2_0_ogrid
  endinterface
  interface u_dot_grad_ogrid
    module procedure u_dot_grad_scl_ogrid
    module procedure u_dot_grad_vec_ogrid
  endinterface
  interface read_snap_ogrid
    module procedure read_snap_double_ogrid
    module procedure read_snap_single_ogrid
  endinterface
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
!  Define the geometry of the solids and construct the grid around the solid
!  object.
!  Currently only a single solid object with overlapping grid is implemented.
!  Use solid_cells.f90 without overlapping grids if several solids inside the
!  flow domain is desired.
!
!  feb--apr-17/Jorgen: Coded
!
      use Solid_Cells_Mpicomm, only: initialize_mpicomm_ogrid
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer :: i
!
      if (cylinder_radius <= 0) then
        call fatal_error('initialize_solid_cells_ogrid', &
            'All cylinders must have non-zero radii!')
      endif
      if(r_ogrid <= 0) r_ogrid=3.*cylinder_radius

      if (lroot) then
        print*, 'nxgrid_ogrid, nygrid_ogrid, nzgrid_ogrid=', nxgrid_ogrid, nygrid_ogrid, nzgrid_ogrid
        print*, 'Cylidner radius=',cylinder_radius
        print*, 'Cylindrical grid radius=',r_ogrid
      endif
!  Interlap parameters of the ogrid
      xyz0_ogrid(1)=cylinder_radius
      xyz0_ogrid(2)=-pi
      xyz0_ogrid(3)=xyz0(3)
      Lxyz_ogrid(1)=r_ogrid-xyz0_ogrid(1)
      Lxyz_ogrid(2)=2*pi
      Lxyz_ogrid(3)=Lxyz(3)
!  Position related to the flow domain and cartesian grid
      xorigo_ogrid(1) = cylinder_xpos
      xorigo_ogrid(2) = cylinder_ypos
      xorigo_ogrid(3) = cylinder_zpos
!
! Try to find flow direction
!
      flow_dir = 0
      if (fbcx(1,1) > 0)     then; flow_dir = 1
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
!  If velocity for initial potential flow is not set, use boundary condition to set this
!
      if(init_uu==0.) then
        if(flow_dir==1)     then; init_uu=fbcx(1,1)
        elseif(flow_dir==2) then; init_uu=fbcy(2,1)
        elseif(flow_dir==3) then; init_uu=fbcz(3,1)
        endif
        print*, 'By using fbc[x,y,z] I set the initial velocity to ',init_uu
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
!  Initialize the pencils overlapping grid and construct arrays
!  for pencil indices
!
      call initialize_pencils_ogrid(0.0)
      call setup_mm_nn_ogrid
!
!  Construct overlapping grid
!
      call construct_grid_ogrid
!
!  Initialize overlapping grid
!
      call initialize_grid_ogrid
! 
!  Set interpolation zone for curvilinear to cartesian
!
      r_int_outer=r_ogrid-inter_stencil_len*dxmax_ogrid*sqrt(2.)
      r_int_inner=xyz0_ogrid(1)
      if(lroot) then
        if(.not.lequidist_ogrid(1)) then
          print*, 'Non-linear grid in radial direction - dx_rcyl, dx_rogrid:', &
              0.5*(xglobal_ogrid(nghost+1)-xglobal_ogrid(nghost-1)), &
              0.5*(xglobal_ogrid(mxgrid_ogrid-nghost+1)-xglobal_ogrid(mxgrid_ogrid-nghost-1))
          print*, 'Theta grid spacing - r_cyl*dy_ogrid,r_int_outer*dy_ogrid,r_ogrid*dy_ogrid',&
              xyz0_ogrid(1)*(yglobal_ogrid(2)-yglobal_ogrid(1)), &
              r_int_outer*(yglobal_ogrid(2)-yglobal_ogrid(1)), &
              r_ogrid*(yglobal_ogrid(2)-yglobal_ogrid(1))
          print*, 'Cartesian grid spacing - dx, dy, dz:', dx,dy,dz
          print*, 'Timestep factor:', ceiling(dxmin/dxmin_ogrid)
        else
          print*, 'Radial grid spacing - dx_ogrid:', &
              0.5*(xglobal_ogrid(nghost+1)-xglobal_ogrid(nghost-1))
          print*, 'Theta grid spacing - r_cyl*dy_ogrid,r_int_outer*dy_ogrid,r_ogrid*dy_ogrid',&
              xyz0_ogrid(1)*(yglobal_ogrid(2)-yglobal_ogrid(1)), &
              r_int_outer*(yglobal_ogrid(2)-yglobal_ogrid(1)), &
              r_ogrid*(yglobal_ogrid(2)-yglobal_ogrid(1))
          print*, 'Cartesian grid spacing - dx, dy, dz:', dx,dy,dz
          print*, 'Timestep factor:', ceiling(dxmin/dxmin_ogrid)
        endif
        if(.not.lequidist_ogrid(2)) then
          call fatal_error('initialize_solid_cells','non-linear grid in theta direction not allowed')
        endif
        if(.not.lequidist_ogrid(3))  print*, 'Non-linear grid in z-direction'
      endif

!
!  Read data.
!  Snapshot data are saved in the data subdirectory.
!  This directory must exist, but may be linked to another disk.
!
      if(lrun) call rsnap_ogrid('ogvar.dat',lread_nogrid)
!
!  The following is here to avoid division in sub.f90 for diagnostic
!  outputs of integrated values in the non equidistant case.
!  Do this even for uniform meshes, in which case xprim=dx, etc.
!  Remember that dx_1=0 for runs without extent in that direction.
!
      if (nxgrid==1) then; xprim=1.0; else; xprim=1./dx_1; endif
      if (nygrid==1) then; yprim=1.0; else; yprim=1./dy_1; endif
      if (nzgrid==1) then; zprim=1.0; else; zprim=1./dz_1; endif

!
!  Construct summation by parts-stencils, if SBP is on
!
      if(SBP) then
        D1_SBP(1,:)=(/ -21600./13649.  , 104009./54596.  , 30443./81894.   , & 
                       -33311./27298.  , 16863./27298.   , -15025./163788. , &
                       0.              , 0.              , 0.              /)
        D1_SBP(2,:)=(/ -104009./240260., 0.              , -311./72078.    , & 
                       20229./24026.   , -24337./48052.  , 36661./360390.  , &
                       0.              , 0.              , 0.              /)
        D1_SBP(3,:)=(/ -30443./162660.  , 311./32532.     , 0.              , & 
                       -11155./16266.  , 41287./32532.   , -21999./54220.  , &
                       0.              , 0.              , 0.              /)
        D1_SBP(4,:)=(/ 33311./107180.  , -20229./21436.  , 485./1398.      , & 
                       0.              , 4147./21436.    , 25427./321540.  , &
                       72./5359.       , 0.              , 0.              /)
        D1_SBP(5,:)=(/ -16863./78770.  , 24337./31508.   , -41287./47262.  , & 
                       -4147./15754.   , 0.              , 342523./472620. , &
                       -1296./7877.    , 144./7877.      , 0.              /)
        D1_SBP(6,:)=(/ 15025./525612.  , -36661./262806. , 21999./87602.   , & 
                       -25427./262806. , -342523./525612., 0.              , &
                       32400./43801.   , -6480./43801.   , 720./43801.     /)
        D2_SBP(1,:)=(/ 114170./40947.  , -438107./54596. ,  336409./40947. , & 
                       -276997./81894. ,  3747./13649.   , 21035./163788.  , &
                       0.              , 0.              , 0.              /)
        D2_SBP(2,:)=(/ 6173./5860.     , -2066./879.     ,  3283./1758.    , & 
                       -303./293.      ,  2111./3516.    , -601./4395.     , &
                       0.              , 0.              , 0.              /)
        D2_SBP(3,:)=(/ -52391./81330.  ,  134603./32532. ,  -21982./2711.  , & 
                       112915./16266.  , -46969./16266.  , 30409./54220.   , &
                       0.              , 0.              , 0.              /)
        D2_SBP(4,:)=(/ 68603./321540.  , -12423./10718.  ,  112915./32154. , & 
                       -75934./16077.  ,  53369./21436.  , -54899./160770. , &
                       48./5359.       , 0.              , 0.              /)
        D2_SBP(5,:)=(/ -7053./39385.   ,  86551./94524.  ,  -46969./23631. , & 
                       53369./15754.   , -87904./23631.  , 820271./472620. , &
                       -1296./7877.    , 96./7877.       , 0.              /)
        D2_SBP(6,:)=(/ 21035./525612.  , -24641./131403. ,  30409./87602.  , & 
                       -54899./131403. ,  820271./525612., -117600./43801. , &
                       64800./43801.   , -6480./43801.   , 480./43801.     /)
      endif
!
!  Set up necessary units for equation of state
!
      call initialize_eos
!
!  Check if it will be necessary to use communication between processors 
!  to perform the interpolation between the overlapping grids.
!  Build serial arrays and data structures necessary to performe communication
!  and interpolation efficiently.
!
      call construct_serial_bdry_cartesian
      call construct_serial_bdry_curv
      call initialize_interpolate_points
      !call initialize_send_ip_points
      call initialize_send_ip_points_alt
! 
!  Allocate arrays used for communications across processors internally on the ogrid
!
      call initialize_mpicomm_ogrid

      if(lroot) then
        print*, 'Interpolation zone: r_ogrid, r_int_outer, r_int_inner',r_ogrid,r_int_outer,r_int_inner
        print*, 'Timestep factor',  max(1,ceiling(dxmin/dxmin_ogrid))
      endif
    end subroutine initialize_solid_cells
!***********************************************************************
    subroutine init_solid_cells(f)
!
!  Initial conditions for cases where we have solid structures in the domain.
!  Typically the flow field is set such that we have no-slip conditions
!  at the solid structure surface.
!
!  This routine is adapted to resolve the solid structures with overlapping 
!  curvilinear grids. Hence, the f-array, which is the array of flow variables
!  on the carteisian grid, must use potential flow to set the ghost zones that
!  will be set by intepolation for t>0.
! 
!  Only implemented for a single solid at present!
!
!  XX-feb-17/Jorgen: Coded
!
      use Initcond, only: gaunoise
      use Sub,      only: wdim,control_file_exists
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real :: a2, rr2, rr2_low, rr2_high, theta_low, theta_high
      real :: wall_smoothing,wall_smoothing_temp
      real :: Lorth,flowx,flowy,shift_flow,shift_orth,flow_r,orth_r
      integer i,j,cyl,iflow,iorth
      real :: shift_top,shift_bot,r_k_top,r_k_bot,theta_k_top,theta_k_bot
      real :: ur_k_top ,ur_k_bot ,uth_k_top,uth_k_bot
      logical :: lnoerase=.false.
!
!  Cartesian array f:
!
!  Set up variables in the appropriate direction of the flow
      if(abs(flow_dir)==1.or.initsolid_cells=='cylinderstream_x') then
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
      elseif(abs(flow_dir)==2.or.initsolid_cells=='cylinderstream_y') then
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
        if (lroot) print*,'No such value for init_solid_cells:', initsolid_cells
        call fatal_error('init_solid_cells','initialization of cylinderstream')
      endif
!
!  Stream functions for flow around a cylinder as initial condition.
!
      call gaunoise(ampl_noise,f,iux,iuz)
      f(:,:,:,iflow) = f(:,:,:,iflow)+init_uu
      shift_flow = 0
      a2=xyz0_ogrid(1)**2
      do i = l1,l2
        do j = m1,m2
! Choose correct points depending on flow direction
          flow_r=(x(i)-xorigo_ogrid(1))*flowx+(y(j)-xorigo_ogrid(2))*flowy
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
!  Velocities inside the solid objects are set to zero 
            f(i,j,:,iux:iuz)=0.
            if (ilnTT /= 0) then
              f(i,j,:,ilnTT) = cylinder_temp
              f(i,j,:,ilnrho) = f(l2,m2,n2,ilnrho) &
                *f(l2,m2,n2,ilnTT)/cylinder_temp
            endif
          endif
        enddo
      enddo
!
!  Cylindrical array f_ogrid:
!
!  Stream functions for flow around a cylinder as initial condition.
!  Note that here ux and uy are the velocities in r and theta directions, respectively
!
!  Rotate system if the flow is in y-direction
      flowy=-flowy*pi*0.5
      f_ogrid(:,:,:,iux:iuz)=0.
      if(ldensity_nolog) then
        f_ogrid(:,:,:,irho)=1.
      else
        f_ogrid(:,:,:,ilnrho)=0.
      endif
      call gaunoise_ogrid(ampl_noise,iux,iuz)
      do i=l1_ogrid,l2_ogrid+nghost
        rr2=x_ogrid(i)**2
        wall_smoothing = 1-exp(-(rr2-a2)/skin_depth**2)
        do j=m1_ogrid,m2_ogrid
!  Compute potential flow past single cylinder
          f_ogrid(i,j,:,iux) = +init_uu*(1-a2/rr2)*cos(y_ogrid(j)+flowy)
          f_ogrid(i,j,:,iuy) = -init_uu*(1+a2/rr2)*sin(y_ogrid(j)+flowy)
          if (ilnTT /= 0) then
            wall_smoothing_temp = 1-exp(-(rr2-a2)/(sqrt(a2))**2)
            f_ogrid(i,j,:,ilnTT) = wall_smoothing_temp*f_ogrid(i,j,:,ilnTT) &
              +cylinder_temp*(1-wall_smoothing_temp)
            f_ogrid(i,j,:,ilnrho) = f_ogrid(l2_ogrid,m2_ogrid,n2_ogrid,ilnrho) &
              *f_ogrid(l2_ogrid,m2_ogrid,n2_ogrid,ilnTT)/f_ogrid(i,j,:,ilnTT)
          endif
!  Compute contribution to flow from cylinders above and below, due to periodic boundary conditions
          do cyl = 1,100
            shift_top = cyl*Lorth
            shift_bot = cyl*Lorth
            r_k_top=sqrt(x_ogrid(i)**2-2*x_ogrid(i)*shift_top*sin(y_ogrid(j)+flowy)+shift_top**2)
            r_k_bot=sqrt(x_ogrid(i)**2-2*x_ogrid(i)*shift_bot*sin(y_ogrid(j)+flowy)+shift_bot**2)
            theta_k_top=atan2(x_ogrid(i)*sin(y_ogrid(j)+flowy)-shift_top,(x_ogrid(i)*cos(y_ogrid(j)+flowy)))
            theta_k_bot=atan2(x_ogrid(i)*sin(y_ogrid(j)+flowy)-shift_bot,(x_ogrid(i)*cos(y_ogrid(j)+flowy)))
            ur_k_top =init_uu*a2*( x_ogrid(i)*cos(theta_k_top)-shift_top*sin(theta_k_top-(y_ogrid(j)+flowy)))/(r_k_top**3)
            ur_k_bot =init_uu*a2*( x_ogrid(i)*cos(theta_k_bot)-shift_bot*sin(theta_k_bot-(y_ogrid(j)+flowy)))/(r_k_bot**3)
            uth_k_top=init_uu*a2*(-x_ogrid(i)*sin(theta_k_top)+shift_top*cos(theta_k_top-(y_ogrid(j)+flowy)))/(r_k_top**3)
            uth_k_bot=init_uu*a2*(-x_ogrid(i)*sin(theta_k_bot)+shift_bot*cos(theta_k_bot-(y_ogrid(j)+flowy)))/(r_k_bot**3)
            f_ogrid(i,j,:,iux) = f_ogrid(i,j,:,iux)+(ur_k_top+ur_k_bot)
            f_ogrid(i,j,:,iuy) = f_ogrid(i,j,:,iuy)+(uth_k_top+uth_k_bot)
          enddo
        enddo
!  Force no-slip condition on the cylinder surface
        f_ogrid(i,:,:,iux:iuy)=f_ogrid(i,:,:,iux:iuy)*wall_smoothing
      enddo
!
!  Write initial condition to disk.
!
      if (lwrite_ic) then
        call wsnap_ogrid('OGVAR0',ENUM=.false.,FLIST='ogvarN.list')
      endif
!
!  The option lnowrite writes everything except the actual var.dat file.
!  This can be useful if auxiliary files are outdated, and don't want
!  to overwrite an existing var.dat.
!
      lnoerase = control_file_exists("NOERASE")
      if (.not.lnowrite .and. .not.lnoerase) then
        call wsnap_ogrid('ogvar.dat',ENUM=.false.)
      endif
!
!  Write ogdim.dat files, local and global
!
      call wdim(trim(directory)//'/ogdim.dat',mx_ogrid,my_ogrid,mz_ogrid)
      if (lroot) then
        call wdim(trim(datadir)//'/ogdim.dat', &
            mxgrid_ogrid,mygrid_ogrid,mzgrid_ogrid,lglobal=.true.)
      endif
    endsubroutine init_solid_cells
!***********************************************************************
    subroutine initialize_pencils_ogrid(penc0)
!
!  Initialize all pencils that are necessary on ogrid
!
!  14-feb-17/Jorgen: Coded
!
      real :: penc0
!
      p_ogrid%x_mn=penc0
      p_ogrid%y_mn=penc0
      p_ogrid%z_mn=penc0
      p_ogrid%rcyl_mn=penc0
      p_ogrid%phi_mn=penc0
      p_ogrid%rcyl_mn1=penc0
      p_ogrid%fpres=penc0
      p_ogrid%fvisc=penc0
      p_ogrid%rho=penc0
      p_ogrid%rho1=penc0
      p_ogrid%lnrho=penc0
      p_ogrid%grho=penc0
      p_ogrid%glnrho=penc0
      p_ogrid%ugrho=penc0
      p_ogrid%sglnrho=penc0
      p_ogrid%uu=penc0
      p_ogrid%u2=penc0
      p_ogrid%uij=penc0
      p_ogrid%sij=penc0
      p_ogrid%sij2=penc0
      p_ogrid%ugu=penc0
      p_ogrid%ugu2=penc0
      p_ogrid%del2u=penc0
      p_ogrid%divu=penc0  
      p_ogrid%lnTT=penc0
      p_ogrid%cs2=penc0
      p_ogrid%pp=penc0
      p_ogrid%ss=penc0
    
      lpencil_ogrid=.true.
!
    endsubroutine initialize_pencils_ogrid
!***********************************************************************
    subroutine initialize_eos
!  
!  Set up parameters necessary to compute the energy and pressure using
!  the ideal gas eos.
!  This is done in the unit_eos routine in eos_idealgas.f90 for the 
!  cartesian solver
!
!  4-apr-17/Jorgen: Coded
!
      use EquationOfState, only: get_cv1,get_cp1,cs20,gamma_m1
      real :: cp1, cp
!
!  Inverse cv and cp values.
!
      call get_cp1(cp1)
      cp=1./cp1
!
      rho0=1.0
      lnrho0=log(rho0)
      if (gamma_m1/=0.0) then
        lnTT0=log(cs20/(cp*gamma_m1))  !(general case)
      else
        lnTT0=log(cs20/cp)  !(isothermal/polytropic cases: check!)
      endif
    endsubroutine initialize_eos
!***********************************************************************
    subroutine initialize_interpolate_points
!
! Build arrays of interpolation data on processors that perform interpolation.
! Necessary to perform communications in an efficient manner.
!
! apr-17/Jorgen: Coded
!
      real, dimension(3) :: xyz,rthz
      real :: xr,yr
      integer :: i,j,k,ii
!
!  Only implemented for trilinear interpolation between grids
!
      if(.not.grid_interpolation=='trilinear') then
        call fatal_error('initialize_interpolate_points','only implemented for trilinear interpolation')
      endif
!
!  Set up interpolation stencil and data for interpolation from cartesian
!  to curvilinear grid
!  Only done for processors at containing the end points of the radial values.
!
!  If interpolation point requires data from outside this processors domain,
!  set find the grid points and processor id for communication.
!
      if(llast_proc_x) then
        n_ip_cart_to_curv=ny_ogrid*nz_ogrid*nghost
        allocate(cartesian_to_curvilinear(n_ip_cart_to_curv))
!
        ii=0
        do k=n1_ogrid,n2_ogrid
          do j=m1_ogrid,m2_ogrid
            do i=l2_ogrid+1,l2_ogrid+nghost
              ii=ii+1
              xyz=(/ x_ogrid(i)*cos(y_ogrid(j))+xorigo_ogrid(1), &
                      x_ogrid(i)*sin(y_ogrid(j))+xorigo_ogrid(2), &
                      z_ogrid(k) /)
              cartesian_to_curvilinear(ii)%i_xyz = (/ i,j,k /)
              cartesian_to_curvilinear(ii)%xyz = xyz
              call find_near_cart_ind_global(cartesian_to_curvilinear(ii)%ind_global_neighbour,xyz)
              if(.not. this_proc_cartesian(xyz)) then
                call find_proc_cartesian(xyz,cartesian_to_curvilinear(ii)%from_proc)
                if(cartesian_to_curvilinear(ii)%from_proc==iproc) then
                  ! Some ghost points might have this_proc!=iproc, but still return iproc from find_proc_cartesan
                  call ind_global_to_local_cart(cartesian_to_curvilinear(ii)%ind_global_neighbour, &
                        cartesian_to_curvilinear(ii)%ind_local_neighbour,lcheck_init_interpolation)
                endif
              else
                cartesian_to_curvilinear(ii)%from_proc=iproc
                call ind_global_to_local_cart(cartesian_to_curvilinear(ii)%ind_global_neighbour, &
                      cartesian_to_curvilinear(ii)%ind_local_neighbour,lcheck_init_interpolation)
              endif
            enddo
          enddo
        enddo
      endif
!
!  Set up interpolation stencil and data for interpolation from curvilinear
!  to cartesian grid
!  Here we do not know beforehand what points are needed, so we must iterate 
!  through all points. Do this twice, once to set the size of the arrays of 
!  interpolation data, and once to set the data values (after allocating arrays).
!
!  If interpolation point requires data from outside this processors domain,
!  find the grid points and processor id for communication.
!
      do k=n1,n2
        do j=m1,m2
          do i=l1,l2
            xr=x(i)-xorigo_ogrid(1)
            yr=y(j)-xorigo_ogrid(2)
            rthz(1)=sqrt(xr**2+yr**2)
            if((rthz(1)<=r_int_outer) .and. rthz(1)>=r_int_inner) then 
              n_ip_curv_to_cart=n_ip_curv_to_cart+1
            endif
          enddo
        enddo
      enddo
      allocate(curvilinear_to_cartesian(n_ip_curv_to_cart))
!
      ii=0
      do k=n1,n2
        do j=m1,m2
          do i=l1,l2
            xr=x(i)-xorigo_ogrid(1)
            yr=y(j)-xorigo_ogrid(2)
            rthz=(/ sqrt(xr**2+yr**2),atan2(yr,xr),z(k) /)
            if((rthz(1)<=r_int_outer) .and. rthz(1)>=r_int_inner) then 
              ii=ii+1
              curvilinear_to_cartesian(ii)%i_xyz = (/ i,j,k /)
              curvilinear_to_cartesian(ii)%xyz = rthz
              call find_near_curv_ind_global(curvilinear_to_cartesian(ii)%ind_global_neighbour,rthz)
              if(.not. this_proc_curvilinear(rthz,lcheck_init_interpolation)) then
                call find_proc_curvilinear(rthz,curvilinear_to_cartesian(ii)%from_proc)
              else
                curvilinear_to_cartesian(ii)%from_proc=iproc
                call ind_global_to_local_curv(curvilinear_to_cartesian(ii)%ind_global_neighbour, &
                      curvilinear_to_cartesian(ii)%ind_local_neighbour,lcheck_init_interpolation)
              endif
            endif
          enddo
        enddo
      enddo
!
    endsubroutine initialize_interpolate_points
!***********************************************************************
    subroutine initialize_send_ip_points
!
! Build arrays of interpolation data on processors that contain data 
! necessary for interpolation on other processors. 
!
! apr-17/Jorgen: Coded
!
      use Mpicomm, only: mpirecv_int, mpisend_nonblock_int, mpibarrier, mpiwait
      use Solid_Cells_Mpicomm, only: finalize_isend_init_interpol
      integer :: i,iip,npoint
      integer, dimension(ncpus) :: from_proc_curv_to_cart=0
      integer, dimension(ncpus) :: from_proc_cart_to_curv=0
      integer, dimension(:,:,:), allocatable :: ind_from_proc_curv
      integer, dimension(:,:,:), allocatable :: ind_from_proc_cart
      integer, dimension(:,:), allocatable :: ip_id_curv_to_cart
      integer, dimension(:,:), allocatable :: ip_id_cart_to_curv
      integer, dimension(:), allocatable   :: send_to_curv_to_cart
      integer, dimension(:), allocatable   :: send_to_cart_to_curv
      integer, dimension(:,:), allocatable :: send_data_curv_to_cart
      integer, dimension(:,:), allocatable :: send_data_cart_to_curv
      integer, dimension(:), allocatable   :: send_id_curv_to_cart
      integer, dimension(:), allocatable   :: send_id_cart_to_curv
      integer :: max_from_proc, from_proc
      integer, dimension(2) :: nelements
      integer :: size_arr, npoints_requested
      integer, dimension(:), allocatable   :: tmp_arr1D
      integer, dimension(:,:), allocatable :: tmp_arr2D
      integer :: nreq0D,nreq1D,nreq2D
      integer, dimension(ncpus-1) :: ireq0D,ireq1D
      integer, dimension(3*(ncpus-1)) :: ireq2D
      integer :: iter1,iter2
! TODO: COULD THIS BE MOVED INTO SOLID_CELLS_OGRID_MPICOMM?
      if(n_ip_curv_to_cart>0) then
        do i=1,n_ip_curv_to_cart
          from_proc=curvilinear_to_cartesian(i)%from_proc
          if(from_proc/=iproc) then
! Must access from_proc+1 instead of from_proc, to avoid accessing element 0
            from_proc_curv_to_cart(from_proc+1)=from_proc_curv_to_cart(from_proc+1)+1
          endif
        enddo
      endif
!
      max_from_proc=maxval(from_proc_curv_to_cart)
      if(max_from_proc>0) then
        allocate(ind_from_proc_curv(ncpus,max_from_proc,3))
        allocate(ip_id_curv_to_cart(ncpus,max_from_proc))
        do iip=0,ncpus-1
          if(from_proc_curv_to_cart(iip+1)>0) then
            npoint=0
            do i=1,n_ip_curv_to_cart
              if(curvilinear_to_cartesian(i)%from_proc==iip) then
                npoint=npoint+1
! Must access iip+1 instead of iip, to avoid accessing element 0
                ind_from_proc_curv(iip+1,npoint,:)=curvilinear_to_cartesian(i)%ind_global_neighbour
                ip_id_curv_to_cart(iip+1,npoint)=i
              endif
            enddo
          endif
        enddo
      endif
!
      if(n_ip_cart_to_curv>0) then
        do i=1,n_ip_cart_to_curv
          from_proc=cartesian_to_curvilinear(i)%from_proc
          if(from_proc/=iproc) then
! Must access from_proc+1 instead of from_proc, to avoid accessing element 0
            from_proc_cart_to_curv(from_proc+1)=from_proc_cart_to_curv(from_proc+1)+1
          endif
        enddo
      endif
!
      max_from_proc=maxval(from_proc_cart_to_curv)
      if(max_from_proc>0) then
        allocate(ind_from_proc_cart(ncpus,max_from_proc,3))
        allocate(ip_id_cart_to_curv(ncpus,max_from_proc))
        do iip=0,ncpus-1
         if(from_proc_cart_to_curv(iip+1)>0) then
            npoint=0
            do i=1,n_ip_cart_to_curv
              if(cartesian_to_curvilinear(i)%from_proc==iip) then
                npoint=npoint+1
! Must access iip+1 instead of iip, to avoid accessing element 0
                ind_from_proc_cart(iip+1,npoint,:)=cartesian_to_curvilinear(i)%ind_global_neighbour
                ip_id_cart_to_curv(iip+1,npoint)=i
              endif
            enddo
          endif
        enddo
      endif
! 
!  Arrays containing information about which points should be sent by what processor to this
!  processor has now been created. Now, there should be some communication to let all processors
!  know which grid points they should SEND and who should RECIEVE them.
!
!  Note: Code is repeated twice in stead of being programmed as a function, since some compilers do
!  not support allocatable arrays as in/out from subroutines/functions
!  Use som variant of processor number as unique MPI tag (iip,iip+ncpus,etc.) in communication.
!
!  Curvilinear to Cartesian
!
      nreq0D=0
      nreq1D=0
      nreq2D=0
      do iip=0,ncpus-1
!  Send number of points requested from each processors, and send what points are requested
!  if the number of points is larger than zero.
!  Avoid sending to oneself
        if(iip/=iproc) then
          nreq0D=nreq0D+1
          call mpisend_nonblock_int(from_proc_curv_to_cart(iip+1),iip,iip,ireq0D(nreq0D))
          if(from_proc_curv_to_cart(iip+1)>0) then
            nelements=(/ from_proc_curv_to_cart(iip+1),3 /)
            do i=1,3
              nreq2D=nreq2D+1
              call mpisend_nonblock_int(ind_from_proc_curv(iip+1,1:nelements(1),i),nelements(1),iip,200+i,ireq2D(nreq2D))
            enddo
            nreq1D=nreq1D+1
            call mpisend_nonblock_int(ip_id_curv_to_cart(iip+1,1:nelements(1)),nelements(1),iip,iip+2*ncpus,ireq1D(nreq1D))
          endif
        endif
      enddo
      allocate(send_to_curv_to_cart(0))
      allocate(send_data_curv_to_cart(0,3))
      allocate(send_id_curv_to_cart(0))
      do iip=0,ncpus-1
!  Recieve data from all processors. If any points are requested, create array of request.
!  Avoid recieving from oneself
        if(iip/=iproc) then
          call mpirecv_int(npoints_requested,iip,iproc)
!  Allocation/deallocation in a very inefficient manner, but this is only done during pre-processing
!  so memory effieient code is a priority.
          if(npoints_requested>0) then
!  Expand array
            size_arr=size(send_to_curv_to_cart)
            allocate(tmp_arr1D(size_arr))
            tmp_arr1D = send_to_curv_to_cart
            deallocate(send_to_curv_to_cart)
            allocate(send_to_curv_to_cart(size_arr+npoints_requested))
            send_to_curv_to_cart(1:size_arr)=tmp_arr1D
            deallocate(tmp_arr1D)
            !
            send_to_curv_to_cart(size_arr+1:size_arr+npoints_requested)=iip
            nelements=(/ npoints_requested,3 /)
!  Expand array
            allocate(tmp_arr2D(size_arr,3))
            tmp_arr2D = send_data_curv_to_cart
            deallocate(send_data_curv_to_cart)
            allocate(send_data_curv_to_cart(size_arr+npoints_requested,3))
          
            send_data_curv_to_cart(1:size_arr,:)=tmp_arr2D
            deallocate(tmp_arr2D)
            do i=1,3
              call mpirecv_int(send_data_curv_to_cart(size_arr+1:size_arr+npoints_requested,i),nelements(1),iip,200+i)
            enddo
!  Expand array
            allocate(tmp_arr1D(size_arr))
            tmp_arr1D=send_id_curv_to_cart
            deallocate(send_id_curv_to_cart)
            allocate(send_id_curv_to_cart(size_arr+npoints_requested))
            send_id_curv_to_cart(1:size_arr)=tmp_arr1D
            deallocate(tmp_arr1D)
            call mpirecv_int(send_id_curv_to_cart(size_arr+1:size_arr+npoints_requested),npoints_requested,iip,iproc+2*ncpus)
          endif
        endif
      enddo
      do i=1,nreq0D
        call mpiwait(ireq0D(i))
      enddo
      do i=1,nreq1D
        call mpiwait(ireq1D(i))
      enddo
      do i=1,nreq2D
        call mpiwait(ireq2D(i))
      enddo
      call mpibarrier
      !call finalize_isend_init_interpol(ireq1D,ireq2D,nreq1D,nreq2D)
!
!  Cartesian to curvilinear
!
      nreq1D=0
      nreq2D=0
      do iip=0,ncpus-1
!  Send number of points requested from each processors, and send what points are requested
!  if the number of points is larger than zero.
!  Avoid sending to oneself
        if(iip/=iproc) then
          nreq1D=nreq1D+1
          call mpisend_nonblock_int(from_proc_cart_to_curv(iip+1),iip,iip+3*ncpus,ireq1D(nreq1D))
          if(from_proc_cart_to_curv(iip+1)>0) then
            nelements=(/ from_proc_cart_to_curv(iip+1),3 /)
            nreq2D=nreq2D+2
            call mpisend_nonblock_int(ind_from_proc_cart(iip+1,1:nelements(1),:),nelements,iip,iip+4*ncpus,ireq2D(nreq2D-1))
            call mpisend_nonblock_int(ip_id_cart_to_curv(iip+1,1:nelements(1)),nelements(1),iip,iip+5*ncpus,ireq2D(nreq2D))
          endif
        endif
      enddo
      allocate(send_to_cart_to_curv(0))
      allocate(send_data_cart_to_curv(0,3))
      allocate(send_id_cart_to_curv(0))
      do iip=0,ncpus-1
!  Recieve data from all processors. If any points are requested, create array of request.
!  Avoid recieving from oneself
        if(iip/=iproc) then
          call mpirecv_int(npoints_requested,iip,iproc+3*ncpus)
!  Allocation/deallocation in a very inefficient manner, but this is only done during pre-processing
!  so memory effieient code is a priority.
          if(npoints_requested>0) then
!  Expand array
            size_arr=size(send_to_cart_to_curv)
            allocate(tmp_arr1D(size_arr))
            tmp_arr1D = send_to_cart_to_curv
            deallocate(send_to_cart_to_curv)
            allocate(send_to_cart_to_curv(size_arr+npoints_requested))
            send_to_cart_to_curv(1:size_arr)=tmp_arr1D
            deallocate(tmp_arr1D)
            !
            send_to_cart_to_curv(size_arr+1:size_arr+npoints_requested)=iip
            nelements=(/ npoints_requested,3 /)
!  Expand array
            allocate(tmp_arr2D(size_arr,3))
            tmp_arr2D = send_data_cart_to_curv
            deallocate(send_data_cart_to_curv)
            allocate(send_data_cart_to_curv(size_arr+npoints_requested,3))
            send_data_cart_to_curv(1:size_arr,:)=tmp_arr2D
            deallocate(tmp_arr2D)
            call mpirecv_int(send_data_cart_to_curv(size_arr+1:size_arr+npoints_requested,:),nelements,iip,iproc+4*ncpus)
!  Expand array
            allocate(tmp_arr1D(size_arr))
            tmp_arr1D=send_id_cart_to_curv
            deallocate(send_id_cart_to_curv)
            allocate(send_id_cart_to_curv(size_arr+npoints_requested))
            send_id_cart_to_curv(1:size_arr)=tmp_arr1D
            deallocate(tmp_arr1D)
            call mpirecv_int(send_id_cart_to_curv(size_arr+1:size_arr+npoints_requested),npoints_requested,iip,iproc+5*ncpus)
          endif
        endif
      enddo
      call finalize_isend_init_interpol(ireq1D,ireq2D,nreq1D,nreq2D)
!
!  Deallocate arrays not not needed later
!
      if(allocated(ind_from_proc_curv))  deallocate(ind_from_proc_curv)
      if(allocated(ind_from_proc_cart))  deallocate(ind_from_proc_cart)
      if(allocated(ip_id_curv_to_cart))  deallocate(ip_id_curv_to_cart)
      if(allocated(ip_id_cart_to_curv))  deallocate(ip_id_cart_to_curv)
!
!  Translate recieved global indices to local indices and save the module variables for communication 
!
      size_arr=size(send_data_curv_to_cart(:,1))
      allocate(send_curvilinear_to_cartesian(size_arr))
      do i=1,size_arr
        send_curvilinear_to_cartesian(i)%send_to_proc=send_to_curv_to_cart(i)
        send_curvilinear_to_cartesian(i)%ip_id=send_id_curv_to_cart(i)
        call ind_global_to_local_curv(send_data_curv_to_cart(i,:), &
            send_curvilinear_to_cartesian(i)%i_low_corner,lcheck_init_interpolation)
      enddo
      size_arr=size(send_data_cart_to_curv(:,1))
      allocate(send_cartesian_to_curvilinear(size_arr))
      do i=1,size_arr
        send_cartesian_to_curvilinear(i)%send_to_proc=send_to_cart_to_curv(i)
        send_cartesian_to_curvilinear(i)%ip_id=send_id_cart_to_curv(i)
        call ind_global_to_local_cart(send_data_cart_to_curv(i,:), &
            send_cartesian_to_curvilinear(i)%i_low_corner,lcheck_init_interpolation)
      enddo
!
!  Set some auxiliary parameters to help with the interpolation communication
!  Global to module, not to proc
!
      size_arr=size(send_data_curv_to_cart(:,1))
      n_procs_send_curv_to_cart=0
      if(size_arr>0) then
        n_procs_send_curv_to_cart=1
        do i=2,size_arr
          if(send_curvilinear_to_cartesian(i)%send_to_proc /= &
              send_curvilinear_to_cartesian(i-1)%send_to_proc) then
            n_procs_send_curv_to_cart=n_procs_send_curv_to_cart+1
          endif
        enddo
      endif
      allocate(n_ip_to_proc_curv_to_cart(n_procs_send_curv_to_cart))
      n_ip_to_proc_curv_to_cart=1
      do i=2,size_arr
        if(send_curvilinear_to_cartesian(i)%send_to_proc == &
            send_curvilinear_to_cartesian(i-1)%send_to_proc) then
          n_ip_to_proc_curv_to_cart=n_ip_to_proc_curv_to_cart+1
        endif
      enddo
      size_arr=size(send_data_cart_to_curv(:,1))
      n_procs_send_cart_to_curv=0
      if(size_arr>0) then
        n_procs_send_cart_to_curv=1
        do i=2,size_arr
          if(send_cartesian_to_curvilinear(i)%send_to_proc /= &
              send_cartesian_to_curvilinear(i-1)%send_to_proc) then
            n_procs_send_cart_to_curv=n_procs_send_cart_to_curv+1
          endif
        enddo
      endif
      allocate(n_ip_to_proc_cart_to_curv(n_procs_send_cart_to_curv))
      n_ip_to_proc_cart_to_curv=1
      do i=2,size_arr
        if(send_cartesian_to_curvilinear(i)%send_to_proc == &
            send_cartesian_to_curvilinear(i-1)%send_to_proc) then
          n_ip_to_proc_cart_to_curv=n_ip_to_proc_cart_to_curv+1
        endif
      enddo
!
      n_procs_recv_curv_to_cart=count(from_proc_curv_to_cart.gt.0)
      n_procs_recv_cart_to_curv=count(from_proc_cart_to_curv.gt.0)
      allocate(n_ip_recv_proc_curv_to_cart(n_procs_recv_curv_to_cart))
      allocate(n_ip_recv_proc_cart_to_curv(n_procs_recv_cart_to_curv))
      n_ip_recv_proc_curv_to_cart=pack(from_proc_curv_to_cart,from_proc_curv_to_cart.gt.0)
      n_ip_recv_proc_cart_to_curv=pack(from_proc_cart_to_curv,from_proc_cart_to_curv.gt.0)
      allocate(procs_recv_curv_to_cart(n_procs_recv_curv_to_cart))
      allocate(procs_recv_cart_to_curv(n_procs_recv_cart_to_curv))
      iter1=1
      iter2=1
      do iip=0,ncpus-1
        if(from_proc_cart_to_curv(iip+1)>0) then 
          procs_recv_cart_to_curv(iter1)=iip
          iter1=iter1+1
        endif
        if(from_proc_curv_to_cart(iip+1)>0) then
          procs_recv_curv_to_cart(iter2)=iip
          iter2=iter2+1
        endif
      enddo
      max_send_ip_curv_to_cart=maxval(n_ip_to_proc_curv_to_cart)
      max_recv_ip_curv_to_cart=maxval(n_ip_recv_proc_curv_to_cart)
      max_send_ip_cart_to_curv=maxval(n_ip_to_proc_cart_to_curv)
      max_recv_ip_cart_to_curv=maxval(n_ip_recv_proc_cart_to_curv)
!
!  Deallocate arrays
!
      deallocate(send_to_curv_to_cart)
      deallocate(send_to_cart_to_curv)
      deallocate(send_data_curv_to_cart)
      deallocate(send_data_cart_to_curv)
      deallocate(send_id_curv_to_cart)
      deallocate(send_id_cart_to_curv)
!
!  Make sure that all processors complete this initialization before continuing
!
      call mpibarrier
    endsubroutine initialize_send_ip_points
!***********************************************************************
    subroutine initialize_send_ip_points_alt
!
! Build arrays of interpolation data on processors that contain data 
! necessary for interpolation on other processors. 
!
! apr-17/Jorgen: Coded
!
      use Mpicomm, only: mpirecv_int, mpisend_nonblock_int, mpibarrier, &
                         mpirecv_nonblock_int, mpisend_int, mpiwait, mpibcast_int

      use Solid_Cells_Mpicomm, only: finalize_isend_init_interpol
      integer :: i,iip,npoint
      integer, dimension(ncpus) :: from_proc_curv_to_cart=0
      integer, dimension(ncpus) :: from_proc_cart_to_curv=0
      integer, dimension(:,:,:), allocatable :: ind_from_proc_curv
      integer, dimension(:,:,:), allocatable :: ind_from_proc_cart
      integer, dimension(:,:), allocatable :: ip_id_curv_to_cart
      integer, dimension(:,:), allocatable :: ip_id_cart_to_curv
      integer, dimension(:), allocatable   :: send_to_curv_to_cart
      integer, dimension(:), allocatable   :: send_to_cart_to_curv
      integer, dimension(:,:), allocatable :: send_data_curv_to_cart
      integer, dimension(:,:), allocatable :: send_data_cart_to_curv
      integer, dimension(:), allocatable   :: send_id_curv_to_cart
      integer, dimension(:), allocatable   :: send_id_cart_to_curv
      integer :: max_from_proc, from_proc
      integer, dimension(2) :: nelements
      integer :: size_arr, npoints_requested
      integer, dimension(:), allocatable   :: tmp_arr1D
      integer, dimension(:,:), allocatable :: tmp_arr2D
      integer :: nreq1D,nreq2D
      integer, dimension(ncpus-1) :: ireq1D
      integer, dimension(ncpus-1,3) :: ireq2D
      !
      integer, dimension(ncpus,ncpus) :: from_proc_curv_to_cart_glob=0
      integer, dimension(ncpus,ncpus) :: from_proc_cart_to_curv_glob=0
      integer :: iter, ind_start, ind_stop, ip_recv_tot, ip_send_tot, n_ip_proc
      integer, dimension(2) :: buf_size
      integer, dimension(:), allocatable :: id_bufi, id_bufo
      integer, dimension(:,:), allocatable :: ijk_bufi, ijk_bufo
      ! TODO TODO TODO : Ensure that every ID in cartesian_to_curvilinear array (etc.) are unique!!!
! TODO: COULD THIS BE MOVED INTO SOLID_CELLS_OGRID_MPICOMM?
      if(n_ip_curv_to_cart>0) then
        do i=1,n_ip_curv_to_cart
          from_proc=curvilinear_to_cartesian(i)%from_proc
          if(from_proc/=iproc) then
! Must access from_proc+1 instead of from_proc, to avoid accessing element 0
            from_proc_curv_to_cart(from_proc+1)=from_proc_curv_to_cart(from_proc+1)+1
          endif
        enddo
      endif
!
      max_from_proc=maxval(from_proc_curv_to_cart)

      if(max_from_proc>0) then
        allocate(ind_from_proc_curv(ncpus,max_from_proc,3))
        allocate(ip_id_curv_to_cart(ncpus,max_from_proc))
        do iip=0,ncpus-1
          if(from_proc_curv_to_cart(iip+1)>0) then
            npoint=0
            do i=1,n_ip_curv_to_cart
              if(curvilinear_to_cartesian(i)%from_proc==iip) then
                npoint=npoint+1
! Must access iip+1 instead of iip, to avoid accessing element 0
                ind_from_proc_curv(iip+1,npoint,:)=curvilinear_to_cartesian(i)%ind_global_neighbour
                ip_id_curv_to_cart(iip+1,npoint)=i
              endif
            enddo
          endif
        enddo
      endif
!
      if(n_ip_cart_to_curv>0) then
        do i=1,n_ip_cart_to_curv
          from_proc=cartesian_to_curvilinear(i)%from_proc
          if(from_proc/=iproc) then
! Must access from_proc+1 instead of from_proc, to avoid accessing element 0
            from_proc_cart_to_curv(from_proc+1)=from_proc_cart_to_curv(from_proc+1)+1
          endif
        enddo
      endif
!
      max_from_proc=maxval(from_proc_cart_to_curv)
      if(max_from_proc>0) then
        allocate(ind_from_proc_cart(ncpus,max_from_proc,3))
        allocate(ip_id_cart_to_curv(ncpus,max_from_proc))
        do iip=0,ncpus-1
         if(from_proc_cart_to_curv(iip+1)>0) then
            npoint=0
            do i=1,n_ip_cart_to_curv
              if(cartesian_to_curvilinear(i)%from_proc==iip) then
                npoint=npoint+1
! Must access iip+1 instead of iip, to avoid accessing element 0
                ind_from_proc_cart(iip+1,npoint,:)=cartesian_to_curvilinear(i)%ind_global_neighbour
                ip_id_cart_to_curv(iip+1,npoint)=i
              endif
            enddo
          endif
        enddo
      endif
!
!  Arrays containing information about which points should be sent by what processor to this
!  processor has now been created. Now, there should be some communication to let all processors
!  know which grid points they should SEND and who should RECIEVE them.
!
!  Note: Code is repeated twice in stead of being programmed as a function, since some compilers do
!  not support allocatable arrays as in/out from subroutines/functions
!  Use som variant of processor number as unique MPI tag (iip,iip+ncpus,etc.) in communication.
!
!  Curvilinear to Cartesian
!
      if(lroot) then
        do iip=0,ncpus-1
          if(iip/=root) then
            call mpirecv_int(from_proc_curv_to_cart_glob(iip+1,:),ncpus,iip,110)
            call mpirecv_int(from_proc_cart_to_curv_glob(iip+1,:),ncpus,iip,115)
          else
            from_proc_curv_to_cart_glob(root+1,:)=from_proc_curv_to_cart
            from_proc_cart_to_curv_glob(root+1,:)=from_proc_cart_to_curv
          endif
        enddo
      else
        call mpisend_int(from_proc_curv_to_cart,ncpus,root,110)
        call mpisend_int(from_proc_cart_to_curv,ncpus,root,115)
      endif

      call mpibcast_int(from_proc_curv_to_cart_glob,(/ncpus,ncpus/))
      call mpibcast_int(from_proc_cart_to_curv_glob,(/ncpus,ncpus/))
!
!  Set some auxiliary parameters to help with the interpolation communication
!  Global to module, not to proc
!
      n_procs_send_curv_to_cart=count(from_proc_curv_to_cart_glob(:,iproc+1).gt.0)
      n_procs_send_cart_to_curv=count(from_proc_cart_to_curv_glob(:,iproc+1).gt.0)
      allocate(n_ip_to_proc_curv_to_cart(n_procs_send_curv_to_cart))
      allocate(n_ip_to_proc_cart_to_curv(n_procs_send_cart_to_curv))
      n_ip_to_proc_curv_to_cart=pack(from_proc_curv_to_cart_glob(:,iproc+1),&
        from_proc_curv_to_cart_glob(:,iproc+1).gt.0)
      n_ip_to_proc_cart_to_curv=pack(from_proc_cart_to_curv_glob(:,iproc+1),&
        from_proc_cart_to_curv_glob(:,iproc+1).gt.0)
      allocate(procs_send_curv_to_cart(n_procs_send_curv_to_cart))
      allocate(procs_send_cart_to_curv(n_procs_send_cart_to_curv))
      iter=1
      do iip=0,ncpus-1
        if(from_proc_curv_to_cart_glob(iip+1,iproc+1)>0) then
          procs_send_curv_to_cart(iter)=iip
          iter=iter+1
        endif
      enddo
      iter=1
      do iip=0,ncpus-1
        if(from_proc_cart_to_curv_glob(iip+1,iproc+1)>0) then
          procs_send_cart_to_curv(iter)=iip
          iter=iter+1
        endif
      enddo
      max_send_ip_curv_to_cart=maxval(n_ip_to_proc_curv_to_cart)
      max_send_ip_cart_to_curv=maxval(n_ip_to_proc_cart_to_curv)
!
      n_procs_recv_curv_to_cart=count(from_proc_curv_to_cart.gt.0)
      n_procs_recv_cart_to_curv=count(from_proc_cart_to_curv.gt.0)
      allocate(n_ip_recv_proc_curv_to_cart(n_procs_recv_curv_to_cart))
      allocate(n_ip_recv_proc_cart_to_curv(n_procs_recv_cart_to_curv))
      n_ip_recv_proc_curv_to_cart=pack(from_proc_curv_to_cart,from_proc_curv_to_cart.gt.0)
      n_ip_recv_proc_cart_to_curv=pack(from_proc_cart_to_curv,from_proc_cart_to_curv.gt.0)
      allocate(procs_recv_curv_to_cart(n_procs_recv_curv_to_cart))
      allocate(procs_recv_cart_to_curv(n_procs_recv_cart_to_curv))
      iter=1
      do iip=0,ncpus-1
        if(from_proc_curv_to_cart(iip+1)>0) then
          procs_recv_curv_to_cart(iter)=iip
          iter=iter+1
        endif
      enddo
      iter=1
      do iip=0,ncpus-1
        if(from_proc_cart_to_curv(iip+1)>0) then
          procs_recv_cart_to_curv(iter)=iip
          iter=iter+1
        endif
      enddo
      max_recv_ip_curv_to_cart=maxval(n_ip_recv_proc_curv_to_cart)
      max_recv_ip_cart_to_curv=maxval(n_ip_recv_proc_cart_to_curv)
!
!  SEND/RECV CURVILINEAR TO CARTESIAN
!
!  Build recv buffers for ip_send data
!  That is, bufferes that RECIEVE data from processors REQUESTING interpolation points.
!  Size of buffers set by n_procs_send** etc., since they will contain the send data SENT
!  by this processor during interpolation.
!
      ip_send_tot=sum(from_proc_curv_to_cart_glob(:,iproc+1))
      allocate(ijk_bufi(ip_send_tot,3))
      allocate(id_bufi(ip_send_tot))
!
!  Post non-blocking recieves
!
      ind_start=1
      do iter=1,n_procs_send_curv_to_cart
        ind_stop=ind_start+n_ip_to_proc_curv_to_cart(iter)-1
        iip=procs_send_curv_to_cart(iter)
        buf_size=(/ind_stop-ind_start+1,3/)
        do i=1,3
          call mpirecv_nonblock_int(ijk_bufi(ind_start:ind_stop,i),buf_size(1),iip,200+i,ireq2D(iter,i))
        enddo
        call mpirecv_nonblock_int(id_bufi(ind_start:ind_stop),buf_size(1),iip,210,ireq1D(iter))
        ind_start=ind_stop+1
      enddo
!
!  Build send buffers for ip_send data and post blocking sends.
!  That is, bufferes that SEND data from processors REQUESTING interpolation points.
!  Size of buffers set by n_procs_recv** etc., since they will contain the data RECIEVED 
!  by this processor during interpolation.
!
      ip_recv_tot=sum(from_proc_curv_to_cart)
      allocate(ijk_bufo(ip_recv_tot,3))
      allocate(id_bufo(ip_recv_tot))
      ind_start=1
      do iter=1,n_procs_recv_curv_to_cart
        n_ip_proc=n_ip_recv_proc_curv_to_cart(iter)
        ind_stop=ind_start+n_ip_proc-1
        iip=procs_recv_curv_to_cart(iter)
        print*, 'iproc,iip,n_procs_recv', iproc,iip,n_procs_recv_curv_to_cart
        ijk_bufo(ind_start:ind_stop,:)=ind_from_proc_curv(iip+1,1:n_ip_proc,:)
        id_bufo(ind_start:ind_stop)=ip_id_curv_to_cart(iip+1,1:n_ip_proc)
        buf_size=(/ind_stop-ind_start+1,3/)
        do i=1,3
          call mpisend_int(ijk_bufo(ind_start:ind_stop,i),buf_size(1),iip,200+i)
        enddo
        call mpisend_int(id_bufo(ind_start:ind_stop),buf_size(1),iip,210)
        ind_start=ind_stop+1
      enddo
!
!  Wait for recieved data and build send_curvilinear_to_cartesian data container
!
      allocate(send_curvilinear_to_cartesian(ip_send_tot))
      ind_start=1
      do iter=1,n_procs_send_curv_to_cart
        ind_stop=ind_start+n_ip_to_proc_curv_to_cart(iter)-1
        iip=procs_send_curv_to_cart(iter)
        buf_size=(/ind_stop-ind_start+1,3/)
        send_curvilinear_to_cartesian(ind_start:ind_stop)%send_to_proc=iip
        do i=1,3
          call mpiwait(ireq2D(iter,i))
        enddo
        do i=ind_start,ind_stop
          call ind_global_to_local_curv(ijk_bufi(i,:),&
              send_curvilinear_to_cartesian(i)%i_low_corner,lcheck_init_interpolation)
        enddo
        call mpiwait(ireq1D(iter))
        send_curvilinear_to_cartesian(ind_start:ind_stop)%ip_id=id_bufi(ind_start:ind_stop)
        ind_start=ind_stop+1
      enddo
      deallocate(ijk_bufi)
      deallocate(id_bufi)
      deallocate(ijk_bufo)
      deallocate(id_bufo)
!
!  SEND/RECV CARTESIAN TO CURVILINEAR
!
!  Build recv buffers for ip_send data
!  That is, bufferes that RECIEVE data from processors REQUESTING interpolation points.
!  Size of buffers set by n_procs_send** etc., since they will contain the send data SENT
!  by this processor during interpolation.
!
      ip_send_tot=sum(from_proc_cart_to_curv_glob(:,iproc+1))
      allocate(ijk_bufi(ip_send_tot,3))
      allocate(id_bufi(ip_send_tot))
!
!  Post non-blocking recieves
!
      ind_start=1
      do iter=1,n_procs_send_cart_to_curv
        ind_stop=ind_start+n_ip_to_proc_cart_to_curv(iter)-1
        iip=procs_send_cart_to_curv(iter)
        buf_size=(/ind_stop-ind_start+1,3/)
        do i=1,3
          call mpirecv_nonblock_int(ijk_bufi(ind_start:ind_stop,i),buf_size(1),iip,1200+i,ireq2D(iter,i))
        enddo
        call mpirecv_nonblock_int(id_bufi(ind_start:ind_stop),buf_size(1),iip,1210,ireq1D(iter))
        ind_start=ind_stop+1
      enddo
!
!  Build send buffers for ip_send data and post blocking sends.
!  That is, bufferes that SEND data from processors REQUESTING interpolation points.
!  Size of buffers set by n_procs_recv** etc., since they will contain the data RECIEVED 
!  by this processor during interpolation.
!
      ip_recv_tot=sum(from_proc_cart_to_curv)
      allocate(ijk_bufo(ip_recv_tot,3))
      allocate(id_bufo(ip_recv_tot))
      ind_start=1
      do iter=1,n_procs_recv_cart_to_curv
        n_ip_proc=n_ip_recv_proc_cart_to_curv(iter)
        ind_stop=ind_start+n_ip_proc-1
        iip=procs_recv_cart_to_curv(iter)
        ijk_bufo(ind_start:ind_stop,:)=ind_from_proc_cart(iip+1,1:n_ip_proc,:)
        id_bufo(ind_start:ind_stop)=ip_id_cart_to_curv(iip+1,1:n_ip_proc)
        buf_size=(/ind_stop-ind_start+1,3/)
        do i=1,3
          call mpisend_int(ijk_bufo(ind_start:ind_stop,i),buf_size(1),iip,1200+i)
        enddo
        call mpisend_int(id_bufo(ind_start:ind_stop),buf_size(1),iip,1210)
        ind_start=ind_stop+1
      enddo
!
!  Wait for recieved data and build send_curvilinear_to_cartesian data container
!
      allocate(send_cartesian_to_curvilinear(ip_send_tot))
      ind_start=1
      do iter=1,n_procs_send_cart_to_curv
        ind_stop=ind_start+n_ip_to_proc_cart_to_curv(iter)-1
        iip=procs_send_cart_to_curv(iter)
        buf_size=(/ind_stop-ind_start+1,3/)
        send_cartesian_to_curvilinear(ind_start:ind_stop)%send_to_proc=iip
        do i=1,3
          call mpiwait(ireq2D(iter,i))
        enddo
        do i=ind_start,ind_stop
          call ind_global_to_local_cart(ijk_bufi(i,:),&
              send_cartesian_to_curvilinear(i)%i_low_corner,lcheck_init_interpolation)
        enddo
        call mpiwait(ireq1D(iter))
        send_cartesian_to_curvilinear(ind_start:ind_stop)%ip_id=id_bufi(ind_start:ind_stop)
        ind_start=ind_stop+1
      enddo
      deallocate(ijk_bufi)
      deallocate(id_bufi)
      deallocate(ijk_bufo)
      deallocate(id_bufo)

      call mpibarrier
    endsubroutine initialize_send_ip_points_alt
!***********************************************************************
    subroutine ind_global_to_local_curv(i_rthz_global,i_rthz_local,lcheck)
!
!  Translate global indices to local indices on the curvilinear grid
!
!  18-apr-17/Jorgen: Coded
!
      integer, dimension(3), intent(in) :: i_rthz_global
      integer, dimension(3), intent(out) :: i_rthz_local
      logical, intent(in) :: lcheck
      real, dimension(3) :: rthz_global
      integer :: ii,jj,kk
!
      rthz_global=(/ xglobal_ogrid(i_rthz_global(1)), yglobal_ogrid(i_rthz_global(2)), &
                     zglobal_ogrid(i_rthz_global(3)) /)
      do ii=1,mxgrid_ogrid
        if(abs(rthz_global(1)-x_ogrid(ii))<1e-12) then
          i_rthz_local(1)=ii
          exit
        endif
      enddo
      do jj=1,mygrid_ogrid
        if(rthz_global(2)==y_ogrid(jj)) then
          i_rthz_local(2)=jj
          exit
        endif
      enddo
!
      if(nzgrid_ogrid==1) then
        i_rthz_local(3)=n1
      else
        do kk=1,mzgrid_ogrid
          if(rthz_global(3)==z_ogrid(kk)) then
            i_rthz_local(3)=kk
            exit
          endif
        enddo
      endif
!
      if(i_rthz_local(1)==0) print*, 'ZERO INDEX IN R-DIRECTION!'
      if(lcheck) then
        if((rthz_global(1)-x_ogrid(i_rthz_local(1))>1e-12) .or. &
           (rthz_global(2)-y_ogrid(i_rthz_local(2))>1e-12) .or. &
           (rthz_global(3)-z_ogrid(i_rthz_local(3))>1e-12)) then
          print*, ''
          print*, 'iproc', iproc
          print*, 'Correct global to local conversion not performed'
          print*, 'rthz_global',rthz_global
          print*, 'rthz_local',x_ogrid(i_rthz_local(1)),y_ogrid(i_rthz_local(2)),z_ogrid(i_rthz_local(3))
          print*, 'i_rthz_global',i_rthz_global
          print*, 'i_rthz_local',i_rthz_local
          print*, 'xyz0_loc_ogrid[x,y,z]', xyz0_loc_ogrid
          print*, 'xyz0_loc_ogrid[x,y,z]', xyz1_loc_ogrid
          print*, ''
          call fatal_error('ind_global_to_local_curv','correct local point not found!')
        endif
      endif
    endsubroutine ind_global_to_local_curv
!***********************************************************************
    subroutine ind_global_to_local_cart(i_xyz_global,i_xyz_local,lcheck)
!
!  Translate global indices to local indices on the cartesian grid
!
!  18-apr-17/Jorgen: Coded
!
      integer, dimension(3), intent(in) :: i_xyz_global
      integer, dimension(3), intent(out) :: i_xyz_local
      logical, intent(in) :: lcheck
      real, dimension(3) :: xyz_global
      integer :: ii,jj,kk
!
      xyz_global=(/ xglobal(i_xyz_global(1)), yglobal(i_xyz_global(2)), &
                    zglobal(i_xyz_global(3)) /)
      do ii=1,mxgrid
        if(xyz_global(1)==x(ii)) then
          i_xyz_local(1)=ii
          exit
        endif
      enddo
      do jj=1,mygrid
        if(xyz_global(2)==y(jj)) then
          i_xyz_local(2)=jj
          exit
        endif
      enddo
      if(nzgrid==1) then
        i_xyz_local(3)=n1
      else
        do kk=1,mzgrid
          if(xyz_global(3)==z(kk)) then
            i_xyz_local(3)=kk
            exit
          endif
        enddo
      endif
!
      if(lcheck) then
        if((xyz_global(1)-x(i_xyz_local(1))/=0.) .or. &
           (xyz_global(2)-y(i_xyz_local(2))/=0.) .or. &
           (xyz_global(3)-z(i_xyz_local(3))/=0.)) then
          print*, ''
          print*, 'Correct global to local conversion not performed'
          print*, 'iproc', iproc
          print*, 'xyz_global',xyz_global
          print*, 'xyz_local',x(i_xyz_local(1)),y(i_xyz_local(2)),z(i_xyz_local(3))
          print*, 'i_xyz_global',i_xyz_global
          print*, 'i_xyz_local',i_xyz_local
          print*, 'xyz0_loc[x,y,z]', xyz0_loc
          print*, 'xyz0_loc[x,y,z]', xyz1_loc
          print*, ''
          call fatal_error('ind_global_to_local_cart','correct local point not found!')
        endif
      endif
    endsubroutine ind_global_to_local_cart
!***********************************************************************
!    subroutine expand_int_array(arr,nexpand,size_arr)
!!  Expand dynamic input array
!      integer :: size_arr
!      integer, dimension(size_arr), allocatable, intent(inout) :: arr
!      integer, intent(in) :: nexpand
!      integer, dimension(size_arr) :: tmp_arr
!!
!      tmp_arr = arr
!      deallocate(arr)
!      allocate(arr(size_arr+nexpand))
!      arr(1:size_arr)=tmp_arr
!    endsubroutine expand_int_array
!!***********************************************************************
!    subroutine expand_int_array2(arr,nexpand,size_arr)
!!  Expand dynamic input array
!      integer :: size_arr
!      integer, dimension(size_arr), allocatable, intent(inout) :: arr
!      integer, intent(in) :: nexpand
!      integer, dimension(size_arr) :: tmp_arr
!!
!      tmp_arr = arr
!      deallocate(arr)
!      allocate(arr(size_arr(1)+nexpand,size_arr(2)))
!      arr(1:size_arr(1),:)=tmp_arr
!    endsubroutine expand_int_array2
!!***********************************************************************
    subroutine print_grids_only(num)
!  Print to file
    character(len=16) :: xofile,yofile,xcfile,ycfile
    integer :: i,num

    print*,'iproc,x_og(l1_ogrid-1)-r_cyl,lfirst_proc_x,llast_proc_x',iproc,x_ogrid(l1_ogrid-1)-cylinder_radius,&
      lfirst_proc_x,llast_proc_x

    
    xofile='x_ogrid'
    yofile='y_ogrid'
    xcfile='x_cgrid'
    ycfile='y_cgrid'
    write(xofile,"(A7,I1)") trim(xofile),num
    write(yofile,"(A7,I1)") trim(yofile),num
    write(xcfile,"(A7,I1)") trim(xcfile),num
    write(ycfile,"(A7,I1)") trim(ycfile),num

    xofile=trim(xofile)//'.dat'
    yofile=trim(yofile)//'.dat'
    xcfile=trim(xcfile)//'.dat'
    ycfile=trim(ycfile)//'.dat'
    open(unit=1,file=trim(xofile))
    open(unit=2,file=trim(yofile))
    open(unit=3,file=trim(xcfile))
    open(unit=4,file=trim(ycfile))
    do i=l1_ogrid-nghost,l2_ogrid+nghost
      write(1,*) x_ogrid(i)*cos(y_ogrid(m1_ogrid:m2_ogrid))
      write(2,*) x_ogrid(i)*sin(y_ogrid(m1_ogrid:m2_ogrid))
    enddo
    do i=l1-nghost,l2+nghost
      write(3,*) x(i)
    enddo
    write(4,*) y(m1-nghost:m2+nghost)

    close(1)
    close(2)
    close(3)
    close(4)

    endsubroutine print_grids_only
!***********************************************************************
    logical function this_proc_cartesian(xyz)
      !TODO: Should perhaps use xyz0_loc and xyz1_loc for this
      !      to make it consistent with find_proc_cartesian
!
!  Check if the grid points needed for interpolation between from cartesian
!  to curvilinear grid are prestemt on this processor.
!  At present only valid for trilinear interpolation, where points the 
!  four points [(x_i,y_j,z_k) for i=ii:ii+1 etc] are needed
!
!  06-apr-17/Jorgen: Coded
!
      real, dimension(3), intent(in) :: xyz
!
      this_proc_cartesian=.true.
      if(xyz(1)<x(l1).or.xyz(1)>x(l2)) this_proc_cartesian=.false.
      if(xyz(2)<y(m1).or.xyz(2)>y(m2)) this_proc_cartesian=.false.
      if(xyz(3)<z(n1).or.xyz(3)>z(n2)) this_proc_cartesian=.false.
!
    end function this_proc_cartesian
!***********************************************************************
    logical function this_proc_curvilinear(rthz,lcheck)
      !TODO: Should perhaps use xyz0_loc_ogrid and xyz1_loc_ogrid for this
      !      to make it consistent with find_proc_curvilinear
!
!  Check if the grid points needed for interpolation between from curvilinear
!  to cartesian grid are prestemt on this processor.
!  At present only valid for trilinear interpolation, where points the 
!  four points [(x_i,y_j,z_k) for i=ii:ii+1 etc] are needed
!
!  06-apr-17/Jorgen: Coded
!
      real, dimension(3), intent(in) :: rthz
      logical, intent(in) :: lcheck
!
      this_proc_curvilinear=.true.
      if(rthz(1)<x_ogrid(1).or.rthz(1)>x_ogrid(mx_ogrid)) this_proc_curvilinear=.false.
      if(rthz(2)<y_ogrid(1).or.rthz(2)>y_ogrid(my_ogrid)) this_proc_curvilinear=.false.
      if(rthz(3)<z_ogrid(1).or.rthz(3)>z_ogrid(mz_ogrid)) this_proc_curvilinear=.false.
!
      if(lcheck) then
        if(rthz(1)<xyz0_ogrid(1)) then
          call fatal_error('this_proc_curvilinear','interpolation point is INSIDE the solid cylinder!')
        elseif(rthz(1)>r_ogrid) then
          call fatal_error('this_proc_curvilinear','interpolation point is OUTSIDE the curvilinear grid!')
        endif
      endif
    end function this_proc_curvilinear
!***********************************************************************
    subroutine drag_force_pencils(c_dragx,c_dragy)
!
!  Compute the total fluid force upon the cylinder 
!
!  \vec{F}=\vec{F_p}+\vec{F_s}
!
!  \vec{F_p}=\int{-p_{r=0} d\vec{A}}\limit_A 
!    dA=R*H*\Delta\theta
!    d\vec{A}=dA\hat{r}
!  \vec{F_s}=\int\vec{\tau}dA
!    \tau=\nu\rho_{r=0}(\frac{du_{\theta}}{dr})_{r=0}
!    \vec{\tau}=\tau\hat{\theta}
!   
!  F_x=\vec{F}.\hat{x}=F_p\cos{\theta}-F_s\sin{theta}
!  F_y=\vec{F}.\hat{y}=F_p\sin{\theta}+F_s\cos{theta}
!
!  Normalization computed in find_drag_coeff routine.
!
!  10-apr-17/Jorgen: Coded
!
      use Viscosity, only: getnu
!
      real, intent(inout) :: c_dragx,c_dragy
      real :: F_press,F_shear
      real :: nu
!
      call getnu(nu_input=nu)
!
      F_press=-p_ogrid%pp(1)
      F_shear=nu*p_ogrid%rho(1)*p_ogrid%uij(1,2,1)
      c_dragx=c_dragx+(F_press*cos(y_ogrid(m_ogrid))-F_shear*sin(y_ogrid(m_ogrid)))
      c_dragy=c_dragy+(F_press*sin(y_ogrid(m_ogrid))+F_shear*cos(y_ogrid(m_ogrid)))
!
    endsubroutine drag_force_pencils
!***********************************************************************
    subroutine drag_coeffs(c_dragx,c_dragy)
!
!  Sum up the computed drag on root processor.
!  Normalization done in the end of the computation.
!  Use initial velocity in flow direction and initial density to compute
!  drag coefficients. These should be equal to the inlet velocity. 
!
!  11-apr-17/Jorgen: Coded
!
      use Mpicomm, only: mpireduce_sum
      real, intent(inout) :: c_dragx,c_dragy
      real :: c_dragx_all,c_dragy_all
      real :: norm
!
      norm=dy_ogrid/(nzgrid_ogrid*rho0*init_uu**2)
!
      call mpireduce_sum(c_dragx,c_dragx_all)
      call mpireduce_sum(c_dragy,c_dragy_all)
!
      if(lroot) then
        c_dragx=c_dragx_all*norm
        c_dragy=c_dragy_all*norm
        if (idiag_c_dragx /= 0) fname(idiag_c_dragx)=c_dragx
        if (idiag_c_dragy /= 0) fname(idiag_c_dragy)=c_dragy
      endif
!
    endsubroutine drag_coeffs
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
!  Reads and registers print parameters relevant for solid cells
!
!   mar-2009/kragset: coded
!   nov-2010/kragset: generalized to include drag in z-direction
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset, lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr = lwrite
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_c_dragx = 0
        idiag_c_dragy = 0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname = 1,nname
        call parse_name(iname,cname(iname),cform(iname),'c_dragx',idiag_c_dragx)
        call parse_name(iname,cname(iname),cform(iname),'c_dragy',idiag_c_dragy)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
! TODO: Needed?
      if (lwr) then
      endif
!
    endsubroutine rprint_solid_cells
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
!
!  If we are inside the region of interpolation from curvilinear to cartesian
!  grid, the comutations of the f-array are frozen, by setting df=0 for all
!  variables
!
!  22-feb-17/Jorgen: Coded
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real :: r_int_outer2
      integer :: i
!
      r_int_outer2=r_int_outer**2
!
      do i=l1,l2
        if((x(i)-xorigo_ogrid(1))**2+(y(m)-xorigo_ogrid(2))**2 < r_int_outer2) then
          df(i,m,n,:)=0.
        endif
      enddo
!
    endsubroutine freeze_solid_cells
!***********************************************************************
  function in_solid_cell(part_pos,part_rad)
!
!  Check if the position px,py,pz is within a solid cylinder
!
!  14-feb-17/Jorgen: coded
!
!TODO: The call to this function, in particles_dust, should be made more efficient!
    logical :: in_solid_cell
    real, dimension(3), intent(in) :: part_pos
    real, intent(in) :: part_rad
    real :: r_solid_par
!
    in_solid_cell = .false.
!
    r_solid_par = sqrt(sum((xorigo_ogrid(1:2)-part_pos(1:2))**2))
    if(r_solid_par<=(xyz0_ogrid(1)+part_rad)) then
      in_solid_cell=.true.
    endif
!
  endfunction in_solid_cell
!***********************************************************************
  subroutine pencil_criteria_solid_cells()
!
!  Dummy routine
!
  endsubroutine pencil_criteria_solid_cells
!***********************************************************************
  subroutine solid_cells_clean_up
!
!  Dummy routine
!
  endsubroutine solid_cells_clean_up
!***********************************************************************
  subroutine communicate_ip_cart_to_curv(f_cartesian,ivar1,ivar2)
!
!  Send and recieve necessary information to perform interpolation from 
!  the cartesian to the curvilinear grid.
!  
!  apr-17/Jorgen: Coded
!  
    use Mpicomm, only: mpisend_nonblock_int,mpisend_nonblock_real,mpirecv_int,mpirecv_real,mpiwait,mpibarrier
    real, dimension(mx,my,mz,mfarray), intent(in) :: f_cartesian
    integer, intent(in) :: ivar1,ivar2
    integer, dimension(n_procs_send_cart_to_curv) :: ireq1D, ireq5D
    integer, dimension(5) :: nbuf_farr
    integer, dimension(max_recv_ip_cart_to_curv) :: id_bufi
    real, dimension(max_send_ip_cart_to_curv,2,2,2,ivar2-ivar1+1) :: f_bufo
    real, dimension(max_recv_ip_cart_to_curv,2,2,2,ivar2-ivar1+1) :: f_bufi
    real, dimension(2,2,2,ivar2-ivar1+1) :: farr
    integer :: i,j,k,id,ipp
    integer :: iter, send_to, recv_from
    integer, dimension(3) :: inear_loc
    integer :: ind_send_first, ind_send_last, ind_recv_first, ind_recv_last
!
    nbuf_farr(2:4)=2
    nbuf_farr(5)=ivar2-ivar1+1
!
    ind_send_first=1
    do iter=1,n_procs_send_cart_to_curv
      ind_send_last=n_ip_to_proc_cart_to_curv(iter)+ind_send_first-1
      send_to=send_cartesian_to_curvilinear(ind_send_last)%send_to_proc
      nbuf_farr(1)=ind_send_last-ind_send_first+1
      do ipp=1,nbuf_farr(1)
        i=send_cartesian_to_curvilinear(ind_send_first+ipp-1)%i_low_corner(1)
        j=send_cartesian_to_curvilinear(ind_send_first+ipp-1)%i_low_corner(2)
        k=send_cartesian_to_curvilinear(ind_send_first+ipp-1)%i_low_corner(3)
        f_bufo(ipp,:,:,:,:)=f_cartesian(i:i+1,j:j+1,k:k+1,ivar1:ivar2)
      enddo
      !print*, 'iproc: send id info', iproc,send_cartesian_to_curvilinear(ind_send_first:ind_send_last)%ip_id
      call mpisend_nonblock_int(send_cartesian_to_curvilinear(ind_send_first:ind_send_last)%ip_id, &
        nbuf_farr(1),send_to,send_to,ireq1D(iter))
      call mpisend_nonblock_real(f_bufo(1:nbuf_farr(1),:,:,:,:),nbuf_farr,send_to,send_to+ncpus,ireq5D(iter))
      ind_send_first=ind_send_last+1
    enddo
    ind_recv_first=1
    do iter=1,n_procs_recv_cart_to_curv
      ind_recv_last=n_ip_recv_proc_cart_to_curv(iter)
      recv_from=procs_recv_cart_to_curv(iter)
      nbuf_farr(1)=ind_recv_last-ind_recv_first+1
      call mpirecv_int(id_bufi(1:nbuf_farr(1)),nbuf_farr(1),recv_from,iproc)
      call mpirecv_real(f_bufi(1:nbuf_farr(1),:,:,:,:),nbuf_farr,recv_from,iproc+ncpus)
      do ipp=1,nbuf_farr(1)
        call interpolate_point_cart_to_curv(id_bufi(ipp),ivar1,ivar2,f_bufi(ipp,:,:,:,:),f_cartesian)
      enddo
    enddo
!
!  Interpolate remaining points 
!
    do id=1,n_ip_cart_to_curv
    ! TODO: Make this more efficient by only looping over id's not used above and eliminate if-statement 
      if(cartesian_to_curvilinear(id)%from_proc==iproc) then
        inear_loc=cartesian_to_curvilinear(id)%ind_local_neighbour
        farr(:,:,:,ivar1:ivar2)=f_cartesian(inear_loc(1):inear_loc(1)+1,inear_loc(2):inear_loc(2)+1, &
          inear_loc(3):inear_loc(3)+1,ivar1:ivar2)
        call interpolate_point_cart_to_curv(id,ivar1,ivar2,farr,f_cartesian)
      endif
    enddo
!
!  Finalize nonblocking sends
!
    do iter=1,n_procs_send_cart_to_curv
      call mpiwait(ireq1D(iter))
      call mpiwait(ireq5D(iter))
    enddo
    call mpibarrier
!
  endsubroutine communicate_ip_cart_to_curv
!***********************************************************************
  subroutine communicate_ip_curv_to_cart(f_cartesian,ivar1,ivar2)
!
!  Send and recieve necessary information to perform interpolation from 
!  the curvilinear to the cartesian grid.
!  
!  apr-17/Jorgen: Coded
!  
    use Mpicomm, only: mpisend_nonblock_int,mpisend_nonblock_real,mpirecv_int,mpirecv_real,mpiwait,mpibarrier
    real, dimension(mx,my,mz,mfarray), intent(inout) :: f_cartesian
    integer, intent(in) :: ivar1,ivar2
    integer, dimension(n_procs_send_curv_to_cart) :: ireq1D, ireq5D
    integer, dimension(5) :: nbuf_farr
    integer, dimension(max_recv_ip_curv_to_cart) :: id_bufi
    real, dimension(max_send_ip_curv_to_cart,2,2,2,ivar2-ivar1+1) :: f_bufo
    real, dimension(max_recv_ip_curv_to_cart,2,2,2,ivar2-ivar1+1) :: f_bufi
    real, dimension(2,2,2,ivar2-ivar1+1) :: farr
    integer :: i,j,k,id,ipp
    integer :: iter, send_to, recv_from
    integer, dimension(3) :: inear_loc
    integer :: ind_send_first, ind_send_last, ind_recv_first, ind_recv_last
!
    nbuf_farr(2:4)=2
    nbuf_farr(5)=ivar2-ivar1+1
!
    ind_send_first=1
    do iter=1,n_procs_send_curv_to_cart
      ind_send_last=n_ip_to_proc_curv_to_cart(iter)+ind_send_first-1
      send_to=send_curvilinear_to_cartesian(ind_send_last)%send_to_proc
      nbuf_farr(1)=ind_send_last-ind_send_first+1
      do ipp=1,nbuf_farr(1)
        i=send_curvilinear_to_cartesian(ind_send_first+ipp-1)%i_low_corner(1)
        j=send_curvilinear_to_cartesian(ind_send_first+ipp-1)%i_low_corner(2)
        k=send_curvilinear_to_cartesian(ind_send_first+ipp-1)%i_low_corner(3)
        f_bufo(ipp,:,:,:,:)=f_ogrid(i:i+1,j:j+1,k:k+1,ivar1:ivar2)
      enddo
      call mpisend_nonblock_int(send_curvilinear_to_cartesian(ind_send_first:ind_send_last)%ip_id, &
        nbuf_farr(1),send_to,send_to,ireq1D(iter))
      call mpisend_nonblock_real(f_bufo(1:nbuf_farr(1),:,:,:,:),nbuf_farr,send_to,send_to+ncpus,ireq5D(iter))
      ind_send_first=ind_send_last+1
    enddo
    ind_recv_first=1
    do iter=1,n_procs_recv_curv_to_cart
      ind_recv_last=n_ip_recv_proc_curv_to_cart(iter)
      recv_from=procs_recv_curv_to_cart(iter)
      nbuf_farr(1)=ind_recv_last-ind_recv_first+1
      call mpirecv_int(id_bufi(1:nbuf_farr(1)),nbuf_farr(1),recv_from,iproc)
      call mpirecv_real(f_bufi(1:nbuf_farr(1),:,:,:,:),nbuf_farr,recv_from,iproc+ncpus)
      do ipp=1,nbuf_farr(1)
        call interpolate_point_curv_to_cart(f_cartesian,id_bufi(ipp),ivar1,ivar2,f_bufi(ipp,:,:,:,:))
      enddo
    enddo
!
!  Interpolate remaining points 
!
    do id=1,n_ip_curv_to_cart
    ! TODO: Make more efficient
      if(curvilinear_to_cartesian(id)%from_proc==iproc) then
        inear_loc=curvilinear_to_cartesian(id)%ind_local_neighbour
        farr(:,:,:,ivar1:ivar2)=f_ogrid(inear_loc(1):inear_loc(1)+1,inear_loc(2):inear_loc(2)+1, &
          inear_loc(3):inear_loc(3)+1,ivar1:ivar2)
        call interpolate_point_curv_to_cart(f_cartesian,id,ivar1,ivar2,farr)
      endif
    enddo
!
!  Finalize nonblocking sends
!
    do iter=1,n_procs_send_curv_to_cart
      call mpiwait(ireq1D(iter))
      call mpiwait(ireq5D(iter))
    enddo
    call mpibarrier
!
  endsubroutine communicate_ip_curv_to_cart
!***********************************************************************
  subroutine interpolate_point_cart_to_curv(id,ivar1,ivar2,farr,f_cartesian)
!
!  Use linear interpolation routine to interpolate the values on the cartesian 
!  grid to the interpolation point on the curvilinear grid
!
    real, dimension(mx,my,mz,mfarray), intent(in) :: f_cartesian
    integer, intent(in) :: id,ivar1,ivar2
    real, dimension(2,2,2,ivar2-ivar1+1), intent(in) :: farr
    integer :: i,j,k
    real, dimension(3) :: xyz_ip
    integer, dimension(3) :: inear_glob
    real, dimension(ivar2-ivar1+1) :: f_ip
    !TODO
    integer, parameter :: order=6
    integer :: k_start,k_stop
    real, dimension(order,order,2,ivar2-ivar1+1) :: farr_large
    logical :: highorder
    integer, dimension(3) :: inear_loc
!
    xyz_ip=cartesian_to_curvilinear(id)%xyz
    inear_glob=cartesian_to_curvilinear(id)%ind_global_neighbour
! 
!  Perform interpolation on cartesian grid
!
!TODO TODO
!    k_start=-floor(order/2.)+1
!    k_stop=floor(order/2.)
!    highorder=.true.
!    inear_loc=curvilinear_to_cartesian(id)%ind_local_neighbour
!    farr_large(:,:,:,ivar1:ivar2)=f_cartesian(inear_loc(1)+k_start:inear_loc(1)+k_stop, &
!        inear_loc(2)+k_start:inear_loc(2)+k_stop,inear_loc(3):inear_loc(3)+1,ivar1:ivar2)
!!
!    if(highorder) then
!      if(.not. linear_interpolate_cart_HO(farr_large,ivar1,ivar2,xyz_ip,inear_glob,f_ip,lcheck_interpolation,order)) then
!        call fatal_error('linear_interpolate_curvilinear highorder','interpolation from curvilinear to cartesian')
!      endif
!    else
      if(.not. linear_interpolate_cartesian(farr,ivar1,ivar2,xyz_ip,inear_glob,f_ip,lcheck_interpolation)) then
        call fatal_error('linear_interpolate_cartesian','interpolation from cartesian to curvilinear')
      endif
!    endif
!
!  Update curvilinear grid with the new data values
!
    i=cartesian_to_curvilinear(id)%i_xyz(1)
    j=cartesian_to_curvilinear(id)%i_xyz(2)
    k=cartesian_to_curvilinear(id)%i_xyz(3)
    f_ogrid(i,j,k,iux) = f_ip(iux)*cos(y_ogrid(j)) + f_ip(iuy)*sin(y_ogrid(j))
    f_ogrid(i,j,k,iuy) = -f_ip(iux)*sin(y_ogrid(j)) + f_ip(iuy)*cos(y_ogrid(j))
    f_ogrid(i,j,k,iuz:ivar2) = f_ip(iuz:ivar2)
!
  endsubroutine interpolate_point_cart_to_curv
!***********************************************************************
  subroutine interpolate_point_curv_to_cart(f_cartesian,id,ivar1,ivar2,farr)
!
!  Use linear interpolation routine to interpolate the values on the cartesian 
!  grid to the interpolation point on the curvilinear grid
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f_cartesian
    integer, intent(in) :: id,ivar1,ivar2
    real, dimension(2,2,2,ivar2-ivar1+1), intent(in) :: farr
    integer :: i,j,k
    real, dimension(3) :: xyz_ip
    integer, dimension(3) :: inear_glob
    real, dimension(ivar2-ivar1+1) :: f_ip
    !TODO
    integer, parameter :: order=2
    integer :: k_start,k_stop
    real, dimension(order,order,2,ivar2-ivar1+1) :: farr_large
    logical :: highorder, hermite
    integer, dimension(3) :: inear_loc
!
!
    xyz_ip=curvilinear_to_cartesian(id)%xyz
    inear_glob=curvilinear_to_cartesian(id)%ind_global_neighbour
! 
!  Perform interpolation on cartesian grid
!
!TODO TODO
    !k_start=-floor(order/2.)+1
    !k_stop=floor(order/2.)
    !if(xyz_ip(1)>2*xyz0_ogrid(1)) then
    !  hermite=.false.
    !  !highorder=.false.
    !  highorder=.true.
    !  inear_loc=curvilinear_to_cartesian(id)%ind_local_neighbour
    !  farr_large(:,:,:,ivar1:ivar2)=f_ogrid(inear_loc(1)+k_start:inear_loc(1)+k_stop, &
    !      inear_loc(2)+k_start:inear_loc(2)+k_stop,inear_loc(3):inear_loc(3)+1,ivar1:ivar2)
    !else
    !  highorder=.false.
    !endif
!
    !if(hermite) then
    !  if(.not. hermite_interpolate_curv(farr,ivar1,ivar2,xyz_ip,inear_glob,inear_loc,f_ip,lcheck_interpolation)) then
    !    call fatal_error('hermite_interpolation','interpolation from curvilinear to cartesian')
    !  endif
    !else
    !  if(highorder) then
    !    if(.not. linear_interpolate_curv_HO(farr_large,ivar1,ivar2,xyz_ip,inear_glob,f_ip,lcheck_interpolation,order)) then
    !      call fatal_error('linear_interpolate_curvilinear highorder','interpolation from curvilinear to cartesian')
    !    endif
    !  else
        if(.not. linear_interpolate_curvilinear(farr,ivar1,ivar2,xyz_ip,inear_glob,f_ip,lcheck_interpolation)) then
          call fatal_error('linear_interpolate_curvilinear','interpolation from curvilinear to cartesian')
        endif
    !  endif
    !endif
!
!  Update curvilinear grid with the new data values
!
    i=curvilinear_to_cartesian(id)%i_xyz(1)
    j=curvilinear_to_cartesian(id)%i_xyz(2)
    k=curvilinear_to_cartesian(id)%i_xyz(3)
    f_cartesian(i,j,k,iux)=f_ip(iux)*cos(xyz_ip(2))-f_ip(iuy)*sin(xyz_ip(2))
    f_cartesian(i,j,k,iuy)=f_ip(iux)*sin(xyz_ip(2))+f_ip(iuy)*cos(xyz_ip(2))
    f_cartesian(i,j,k,iuz:ivar2)=f_ip(iuz:ivar2)
!
  endsubroutine interpolate_point_curv_to_cart
!***********************************************************************
  logical function HO_interp_curv_loc(farr,ivar1,ivar2,xxp,inear_glob,inear_loc,fp,lcheck,order)
!
!  Interpolate the value of f to arbitrary (xp, yp) CURVILINEAR coordinate
!  using the high-order lagrangian interpolation.
! 
!  The coefficients are determined by the 2xN grid points surrounding the
!  interpolation point.
! 
!  TEST VERSION: ONLY APPLICABLE FOR SERIAL RUNS
!
!  21-may-17/Jorgen: Adapted from linear_interpolate_curvilinear_HO
!
      integer :: ivar1, ivar2
      integer, intent(in) :: order
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: fp
      real, dimension (order,order,ivar2-ivar1+1) :: farr
      integer, dimension (3) :: inear_glob, inear_loc
      logical :: lcheck
!
      intent(in)  :: farr, ivar1, ivar2, xxp, inear_glob, inear_loc, lcheck
      intent(out) :: fp
!
      integer :: i,ix0,iy0,iz0
      real :: x0,y0,xN,yN,dx,dy,dxN,dyN,xp,yp

      real, dimension(order,ivar2-ivar1+1) :: gp
      integer :: j,k,l
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
      HO_interp_curv_loc=.true.
!
      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3)
!
!  Check if the grid point interval is really correct.
!
      if ((xglobal_ogrid(ix0)<=xxp(1) .and. xglobal_ogrid(ix0+1)>=xxp(1) .or. nxgrid==1) .and. &
          (yglobal_ogrid(iy0)<=xxp(2) .and. yglobal_ogrid(iy0+1)>=xxp(2) .or. nygrid==1) .and. &
          (zglobal_ogrid(iz0)<=xxp(3) .and. zglobal_ogrid(iz0+1)>=xxp(3) .or. nzgrid==1)) then
        ! Everything okay
      else
        print*, 'HO_interp_curv_loc: Global interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc_world
        print*, 'mxgrid_ogrid, xglobal_ogrid(1), xglobal_ogrid(mxgrid_ogrid) = ', & 
            mxgrid_ogrid, xglobal_ogrid(1), xglobal_ogrid(mxgrid_ogrid)
        print*, 'mygrid_ogrid, yglobal_ogrid(1), yglobal_ogrid(mygrid_ogrid) = ', &
            mygrid_ogrid, yglobal_ogrid(1), yglobal_ogrid(mygrid_ogrid)
        print*, 'mzgrid_ogrid, zglobal_ogrid(1), zglobal_ogrid(mzgrid_ogrid) = ', & 
            mzgrid_ogrid, zglobal_ogrid(1), zglobal_ogrid(mzgrid_ogrid)
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
        HO_interp_curv_loc=.false.
        return
      endif
      if ((xglobal_ogrid(ix0)-x_ogrid(inear_loc(1))<10.e-10 ).or. &
          (yglobal_ogrid(iy0)-y_ogrid(inear_loc(2))<10.e-10 )) then
        ! Everything okay
      else
        print*, 'HO_interp_curv_loc: Global and local interpolation point values do not match' 
        print*, 'GLOBAL:'
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
        print*, 'LOCAL:'
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), x_ogrid(inear_loc(1)), x_ogrid(inear_loc(1)+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), y_ogrid(inear_loc(2)), y_ogrid(inear_loc(2)+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), z_ogrid(inear_loc(3)), z_ogrid(inear_loc(3)+1)
        print*, 'DIFF:'
        print*, 'xglobal_ogrid(ix0)-x_ogrid(inear_loc(1)) =', xglobal_ogrid(ix0)-x_ogrid(inear_loc(1))
        print*, 'yglobal_ogrid(iy0)-y_ogrid(inear_loc(2)) =', yglobal_ogrid(iy0)-y_ogrid(inear_loc(2))
        HO_interp_curv_loc=.false.
        return
      endif
!
!  Mapping to quadratic area with corners at (0,0) -- (1,1)
!
    x0 = xglobal_ogrid(ix0-floor(0.5*order)+1)
    y0 = yglobal_ogrid(iy0-floor(0.5*order)+1)
    xN = xglobal_ogrid(ix0+floor(0.5*order))
    yN = yglobal_ogrid(iy0+floor(0.5*order))
!
    dx = xglobal_ogrid(ix0+1)-xglobal_ogrid(ix0)
    dy = yglobal_ogrid(iy0+1)-yglobal_ogrid(iy0)
!
!  Check that gridspacing is correct
!
    if((xN-x0 - (order-1)*dx)<10.e-10 .and. (yN-y0 - (order-1)*dy)<10.e-10) then
      !Do nothing
    else
      print*, 'HO_interp_curv_loc: Grid spacing error'
      print*, 'x0, x1, xN = ', x0,xglobal_ogrid(ix0-floor(0.5*order)+2),xN
      print*, 'dx, N*dx,xN-x0 = ', dx,(order-1)*dx,xN-x0
      print*, 'y0, y1, yN = ', y0,yglobal_ogrid(iy0-floor(0.5*order)+2),yN
      print*, 'dy, N*dy, yN-y0 = ', dy,(order-1)*dy, yN-y0
    endif
!
    dxN = dx*(order-1)
    dyN = dy*(order-1)
    xp = (xxp(1)-x0)/dxN
    yp = (xxp(2)-y0)/dyN
    x0 = 0.
    y0 = 0.
    xN = 1.
    yN = 1.
    dx = 1./(order-1)
    dy = 1./(order-1)

    do i=1,order
      call lagrange1D(farr(:,i,:),ivar1,ivar2,yp,y0,yN,dy,order,gp(i,:))
    enddo
    call lagrange1D(gp,ivar1,ivar2,xp,x0,xN,dx,order,fp) 

!
!  Do a reality check on the interpolation scheme.
!
    if (lcheck) then
      do i=1,2!ivar2-ivar1+1
        if (((fp(i)>maxval(farr(:,:,i)).and.i/=3) .or. (fp(i)<minval(farr(:,:,i)).and.i/=3)).and.ix0>floor(nx_ogrid*0.5)) then
          print*, 'linear_interpolate_curvilinear_HO: interpolated value is smaller or larger than'
          print*, 'linear_interpolate_curvilinear_HO: a values at the corner points, even after linearization!'
          print*, 'linear_interpolate_curvilinear: xxp=', xxp
          print*, 'linear_interpolate_curvilinear: ix0, iy0, iz0=', ix0,iy0,iz0
          print*, 'linear_interpolate_curvilinear: x0, y0, z0=', &
              xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
          print*, 'linear_interpolate_curvilinear: i, fp(i)=', i, fp(i)
          print*, 'linear_interpolate_curvilinear: farr=', farr(:,:,i)
          print*, '------------------'
          print*, 'DETAILS:'
          do j=1,order
            print*, 'j,x(j)',ix0-floor(0.5*order)+1+j-1,xglobal_ogrid(ix0-floor(0.5*order)+1+j-1)
            print*, 'j,y(j)',iy0-floor(0.5*order)+1+j-1,yglobal_ogrid(iy0-floor(0.5*order)+1+j-1)
            print*, 'farr(j,:,1)',farr(j,:,1)
            print*, 'farr(j,:,2)',farr(j,:,2)
          enddo
          print*, '------------------'
          HO_interp_curv_loc=.false.
        endif
        if (fp(i)/=fp(i)) then
          print*, 'linear_interpolate_curvilinear: interpolated value is NaN'
          print*, 'linear_interpolate_curvilinear: xxp=', xxp
          print*, 'linear_interpolate_curvilinear: x0, y0, z0=', &
              xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
          print*, 'linear_interpolate_curvilinear: i, fp(i)=', i, fp(i)
          print*, 'linear_interpolate_curvilinear: farr=', farr(:,:,i)
          print*, '------------------'
          HO_interp_curv_loc=.false.
        endif
      enddo
    endif
!
  endfunction HO_interp_curv_loc
!
!***********************************************************************
  subroutine lagrange1D(farr,ivar1,ivar2,xp,x0,xN,dx,order,fp)
!
      integer, intent(in) :: ivar1, ivar2
      integer, intent(in) :: order
      real, dimension(order,ivar2-ivar1+1), intent(in)  :: farr
      real, intent(in) :: xp, x0, xN, dx
!
      real, dimension (ivar2-ivar1+1), intent(out) :: fp
!
      real, dimension(order) :: l
      real :: xi
      integer :: i,j,ivar
      l = 1.
      do i=1,order
        xi = x0+dx*(i-1)
        do j=1,order
          if(i/=j) then
            l(i) = l(i)*(xp-(x0+dx*(j-1)))/(xi-(x0+dx*(j-1)))
          endif
        enddo
      enddo
      fp=0.
      do ivar=ivar1,ivar2
        do i=1,order
          fp(ivar)=fp(ivar)+l(i)*farr(i,ivar)
        enddo
      enddo
!
  endsubroutine lagrange1D
!***********************************************************************
  subroutine flow_cartesian_to_curvilinear(f_cartesian,f_og)

    use General, only: linear_interpolate
!
!  Interpolate all flow variables from cartesian to curvilinear grid
!  Only need to do this for the radial direction
!
!  Find position in (x,y,z)-coordinates from (r,theta,z)-system
!  Use this to interpolate (linearly) from nearest neighbours
!  Only works for iux:iuz and scalar values (rho,T,etc.) at present.
!
!  NOTE: Does not work for parallell runs
!
!  16-feb-17/Jorgen: Coded
!
    real, dimension (mx,my,mz,mfarray) :: f_cartesian
    real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f_og
    real, dimension (3) :: xyz
    integer :: lower_i,upper_i,lower_j,upper_j,lower_k,upper_k
    integer :: ivar1,ivar2
    integer, dimension (3) :: inear
    real, dimension (ilnrho-iux+1) :: gp
    integer :: i,j,k
!
    ivar1=iux
    ivar2=ilnrho
    do k=n1_ogrid,n2_ogrid
      do j=m1_ogrid,m2_ogrid
        do i=l2_ogrid+1,l2_ogrid+nghost
          xyz=(/ x_ogrid(i)*cos(y_ogrid(j))+xorigo_ogrid(1), &
                  x_ogrid(i)*sin(y_ogrid(j))+xorigo_ogrid(2), &
                  z_ogrid(k) /)
          call find_near_cartesian_indices(lower_i,upper_i,lower_j,upper_j, &
                lower_k,upper_k,xyz)
          inear=(/ lower_i, lower_j, lower_k /)

          if ( .not. linear_interpolate(f_cartesian,ivar1,ivar2,xyz,gp,inear,lcheck_interpolation) ) then
            call fatal_error('linear_interpolate','interpolation from cartesian to curvilinear')
          endif
          f_og(i,j,k,iux)=gp(iux)*cos(y_ogrid(j))+gp(iuy)*sin(y_ogrid(j))
          f_og(i,j,k,iuy)=-gp(iux)*sin(y_ogrid(j))+gp(iuy)*cos(y_ogrid(j))
          f_og(i,j,k,iuz)=gp(iuz)
          f_og(i,j,k,irho)=gp(irho)
        enddo
      enddo
    enddo

  endsubroutine flow_cartesian_to_curvilinear
!***********************************************************************
  subroutine flow_curvilinear_to_cartesian(f_cartesian)
!
!  Interpolate all flow variables from curvilinear to cartesian grid
!
!  NOTE: Does not work for parallell runs
!
    real, dimension (mx,my,mz,mfarray) :: f_cartesian
    real, dimension (3) :: rthz
    integer :: i,j,k
    real :: xr,yr
    integer :: lower_i,upper_i,lower_j,upper_j,lower_k,upper_k
    integer, dimension (3) :: inear
    real, dimension (irho-iux+1) :: gp
    ! TODO
    integer, parameter :: ivar1=1,ivar2=4
    integer, dimension(3) :: inear_glob
    integer, parameter :: order=4
    real, dimension(order,order,ivar2-ivar1+1) :: farr
  
    do k=n1,n2
      do j=m1,m2
        do i=l1,l2
          xr=x(i)-xorigo_ogrid(1)
          yr=y(j)-xorigo_ogrid(2)
          rthz=(/ sqrt(xr**2+yr**2),atan2(yr,xr),z(k) /)
          if((rthz(1)<=r_int_outer) .and.(rthz(1)>=r_int_inner)) then  
            call find_near_curvilinear_indices(lower_i,upper_i,lower_j,upper_j, &
                  lower_k,upper_k,rthz)
            inear=(/ lower_i, lower_j, lower_k /)
            if ( .not. linear_interpolate_ogrid(ivar1,ivar2,rthz,gp,inear,lcheck_interpolation) ) then
              call fatal_error('linear_interpolate_ogrid','interpolation from curvilinear to cartesian')
            endif
            !! ! TODO: Interface for high order interpolation
            !! call find_near_curv_ind_global(inear_glob,rthz)
            !! farr = f_ogrid(inear(1)-floor(order*0.5)+1:inear(1)+floor(order*0.5),&
            !!                inear(2)-floor(order*0.5)+1:inear(2)+floor(order*0.5),1,ivar1:ivar2)
            !! if( .not. HO_interp_curv_loc(farr,ivar1,ivar2,rthz,inear_glob,inear,gp,lcheck_interpolation,order)) then
            !!   call fatal_error('HO_interp_curv_loc','interpolation from curvilinear to cartesian')
            !! endif
            !! ! TODO
            f_cartesian(i,j,k,iux)=gp(iux)*cos(rthz(2))-gp(iuy)*sin(rthz(2))
            f_cartesian(i,j,k,iuy)=gp(iux)*sin(rthz(2))+gp(iuy)*cos(rthz(2))
            f_cartesian(i,j,k,iuz:ivar2)=gp(iuz:ivar2)
          endif
        enddo
      enddo
    enddo
!
  endsubroutine flow_curvilinear_to_cartesian
!***********************************************************************
  logical function linear_interpolate_ogrid(ivar1,ivar2,rthz,gp,inear,lcheck)
!
!  Interpolate the value of g to arbitrary (xp, yp, zp) coordinate on the ogrid
!  using the linear interpolation formula
!
!    g(x,y,z) = A*x*y*z + B*x*y + C*x*z + D*y*z + E*x + F*y + G*z + H .
!
!  The coefficients are determined by the 8 grid points surrounding the
!  interpolation point.
!
!  21-feb-17/Jorgen: Adapted from general.f90
!
    use Cdata
!
    integer :: ivar1, ivar2
    real, dimension (3) :: rthz
    real, dimension (ivar2-ivar1+1) :: gp
    integer, dimension (3) :: inear
!
    real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, g5, g6, g7, g8
    real :: rp0, thp0, zp0
    real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
    integer :: i, ix0, iy0, iz0
    logical :: lfirstcall=.true.,lcheck
!
    intent(in)  :: rthz, ivar1, ivar2, lcheck
    intent(out) :: gp
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
    linear_interpolate_ogrid = .true.
!
    ix0=inear(1); iy0=inear(2); iz0=inear(3)
    if ( (x_ogrid(ix0)>rthz(1)) .and. nxgrid_ogrid/=1) ix0=ix0-1
    if ( (y_ogrid(iy0)>rthz(2)) .and. nygrid_ogrid/=1) iy0=iy0-1
    if ( (z_ogrid(iz0)>rthz(3)) .and. nzgrid_ogrid/=1) iz0=iz0-1
!
!  Check if the grid point interval is really correct.
!
    if ((x_ogrid(ix0)<=rthz(1) .and. x_ogrid(ix0+1)>=rthz(1) .or. nxgrid_ogrid==1) .and. &
        (y_ogrid(iy0)<=rthz(2) .and. y_ogrid(iy0+1)>=rthz(2) .or. nygrid_ogrid==1) .and. &
        (z_ogrid(iz0)<=rthz(3) .and. z_ogrid(iz0+1)>=rthz(3) .or. nzgrid_ogrid==1)) then
      ! Everything okay
    else
      print*, 'linear_interpolate_ogrid: Interpolation point does not ' // &
          'lie within the calculated grid point interval.'
      print*, 'iproc = ', iproc_world
      print*, 'mx_ogrid, x_ogrid(1), x_ogrid(mx_ogrid) = ', mx_ogrid, x_ogrid(1), x_ogrid(mx_ogrid)
      print*, 'my_ogrid, y_ogrid(1), y_ogrid(my_ogrid) = ', my_ogrid, y_ogrid(1), y_ogrid(my_ogrid)
      print*, 'mz_ogrid, z_ogrid(1), z_ogrid(mz_ogrid) = ', mz_ogrid, z_ogrid(1), z_ogrid(mz_ogrid)
      print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
      print*, 'rp, rp0, rp1 = ', rthz(1), x_ogrid(ix0), x_ogrid(ix0+1)
      print*, 'thp, thp0, thp1 = ', rthz(2), y_ogrid(iy0), y_ogrid(iy0+1)
      print*, 'zp, zp0, zp1 = ', rthz(3), z_ogrid(iz0), z_ogrid(iz0+1)
      linear_interpolate_ogrid = .false.
      return
    endif
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!  Set it equal to 0 for dimensions having 1 grid points; this will make sure
!  that the interpolation is bilinear for 2D grids.
!
    rp0=0; thp0=0; zp0=0
    if (nxgrid/=1) rp0=rthz(1)-x_ogrid(ix0)
    if (nygrid/=1) thp0=rthz(2)-y_ogrid(iy0)
    if (nzgrid/=1) zp0=rthz(3)-z_ogrid(iz0)
!
!  Calculate derived grid spacing parameters needed for interpolation.
!  For an equidistant grid we only need to do this at the first call.
!
    if (lequidist_ogrid(1)) then
      if (lfirstcall) dx1=dx_1_ogrid(ix0) !1/dx
    else
      dx1=1/(x_ogrid(ix0+1)-x_ogrid(ix0))
    endif
!
    if (lequidist_ogrid(2)) then
      if (lfirstcall) dy1=dy_1_ogrid(iy0)
    else
      dy1=1/(y_ogrid(iy0+1)-y_ogrid(iy0))
    endif
!
    if (lequidist_ogrid(3)) then
      if (lfirstcall) dz1=dz_1_ogrid(iz0)
    else
      dz1=1/(z_ogrid(iz0+1)-z_ogrid(iz0))
    endif
!
    if ( (.not. all(lequidist_ogrid)) .or. lfirstcall) then
      dxdy1=dx1*dy1; dxdz1=dx1*dz1; dydz1=dy1*dz1
      dxdydz1=dx1*dy1*dz1
    endif
!
!  Function values at all corners.
!
    g1=f_ogrid(ix0  ,iy0  ,iz0  ,ivar1:ivar2)
    g2=f_ogrid(ix0+1,iy0  ,iz0  ,ivar1:ivar2)
    g3=f_ogrid(ix0  ,iy0+1,iz0  ,ivar1:ivar2)
    g4=f_ogrid(ix0+1,iy0+1,iz0  ,ivar1:ivar2)
    g5=f_ogrid(ix0  ,iy0  ,iz0+1,ivar1:ivar2)
    g6=f_ogrid(ix0+1,iy0  ,iz0+1,ivar1:ivar2)
    g7=f_ogrid(ix0  ,iy0+1,iz0+1,ivar1:ivar2)
    g8=f_ogrid(ix0+1,iy0+1,iz0+1,ivar1:ivar2)
!
!  Interpolation formula.
!
    gp = g1 + rp0*dx1*(-g1+g2) + thp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) + &
        rp0*thp0*dxdy1*(g1-g2-g3+g4) + rp0*zp0*dxdz1*(g1-g2-g5+g6) + &
        thp0*zp0*dydz1*(g1-g3-g5+g7) + &
        rp0*thp0*zp0*dxdydz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
!
!  Do a reality check on the interpolation scheme.
!
    if (lcheck) then
      do i=1,ivar2-ivar1+1
        if (gp(i)>max(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
          print*, 'linear_interpolate_ogrid: interpolated value is LARGER than'
          print*, 'linear_interpolate_ogrid: a values at the corner points!'
          print*, 'linear_interpolate_ogrid: xxp=', rthz
          print*, 'linear_interpolate_ogrid: r0, th0, z0=', &
              x_ogrid(ix0), y_ogrid(iy0), z_ogrid(iz0)
          print*, 'linear_interpolate_ogrid: r1, th1, z1=', &
              x_ogrid(ix0+1), y_ogrid(iy0+1), z_ogrid(iz0+1)
          print*, 'linear_interpolate_ogrid: i, gp(i)=', i, gp(i)
          print*, 'linear_interpolate_ogrid: g1...g8=', &
              g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
          print*, '------------------'
        endif
        if (gp(i)<min(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
          print*, 'linear_interpolate_ogrid: interpolated value is smaller than'
          print*, 'linear_interpolate_ogrid: a values at the corner points!'
          print*, 'linear_interpolate_ogrid: xxp=', rthz
          print*, 'linear_interpolate_ogrid: x0, y0, z0=', &
              x_ogrid(ix0), y_ogrid(iy0), z_ogrid(iz0)
          print*, 'linear_interpolate_ogrid: r1, th1, z1=', &
              x_ogrid(ix0+1), y_ogrid(iy0+1), z_ogrid(iz0+1)
          print*, 'linear_interpolate_ogrid: i, gp(i)=', i, gp(i)
          print*, 'linear_interpolate_ogrid: g1...g8=', &
              g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
          print*, '------------------'
        endif
      enddo
    endif
!
    if (lfirstcall) lfirstcall=.false.
!
  endfunction linear_interpolate_ogrid
!***********************************************************************
  logical function linear_interpolate_cartesian(farr,ivar1,ivar2,xxp,inear_glob,fp,lcheck)
!
!  Interpolate the value of f to arbitrary (xp, yp, zp) coordinate
!  using the linear interpolation formula
!
!    g(x,y,z) = A*x*y*z + B*x*y + C*x*z + D*y*z + E*x + F*y + G*z + H .
!
!  The coefficients are determined by the 8 grid points surrounding the
!  interpolation point.
!  Global coordinates are used for the interpolation, to allow interpolation of 
!  values outside this processors domain.
!
!  21-apr-17/Jorgen: Adapted from linear_interpolate in general.f90
!
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (2,2,2,ivar2-ivar1+1) :: farr
      real, dimension (ivar2-ivar1+1) :: fp
      integer, dimension (3) :: inear_glob
      logical :: lcheck
!
      real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, g5, g6, g7, g8
      real :: xp0, yp0, zp0
      real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
      logical :: lfirstcall=.true.
      integer :: ix0, iy0, iz0, i
!
      intent(in)  :: farr, ivar1, ivar2, xxp, inear_glob, lcheck
      intent(out) :: fp
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
      linear_interpolate_cartesian = .true.
!
      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3)
!
!  Check if the grid point interval is really correct.
!
      if ((xglobal(ix0)<=xxp(1) .and. xglobal(ix0+1)>=xxp(1) .or. nxgrid==1) .and. &
          (yglobal(iy0)<=xxp(2) .and. yglobal(iy0+1)>=xxp(2) .or. nygrid==1) .and. &
          (zglobal(iz0)<=xxp(3) .and. zglobal(iz0+1)>=xxp(3) .or. nzgrid==1)) then
        ! Everything okay
      else
        print*, 'linear_interpolate_cartesian: Global interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc_world
        print*, 'mxgrid, xglobal(1), xglobal(mx) = ', mxgrid, xglobal(1), xglobal(mxgrid)
        print*, 'mygrid, yglobal(1), yglobal(my) = ', mygrid, yglobal(1), yglobal(mygrid)
        print*, 'mzgrid, zglobal(1), zglobal(mz) = ', mzgrid, zglobal(1), zglobal(mzgrid)
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal(ix0), xglobal(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal(iy0), yglobal(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal(iz0), zglobal(iz0+1)
        linear_interpolate_cartesian = .false.
        return
      endif
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!  Set it equal to 0 for dimensions having 1 grid points; this will make sure
!  that the interpolation is bilinear for 2D grids.
!
      xp0=0; yp0=0; zp0=0
      if (nxgrid/=1) xp0=xxp(1)-xglobal(ix0)
      if (nygrid/=1) yp0=xxp(2)-yglobal(iy0)
      if (nzgrid/=1) zp0=xxp(3)-zglobal(iz0)
!
!  Calculate derived grid spacing parameters needed for interpolation.
!
      dx1=1/(xglobal(ix0+1)-xglobal(ix0))
      dy1=1/(yglobal(iy0+1)-yglobal(iy0))
      if(nzgrid/=1) then
        dz1=1/(zglobal(iz0+1)-zglobal(iz0))
      else 
        dz1=1
      endif
!
      dxdy1=dx1*dy1; dxdz1=dx1*dz1; dydz1=dy1*dz1
      dxdydz1=dx1*dy1*dz1
!
!  Function values at all corners.
!
      g1=farr(1,1,1,ivar1:ivar2)
      g2=farr(2,1,1,ivar1:ivar2)
      g3=farr(1,2,1,ivar1:ivar2)
      g4=farr(2,2,1,ivar1:ivar2)
      g5=farr(1,1,2,ivar1:ivar2)
      g6=farr(2,1,2,ivar1:ivar2)
      g7=farr(1,2,2,ivar1:ivar2)
      g8=farr(2,2,2,ivar1:ivar2)
!
!  Interpolation formula.
!
      fp = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) + &
          xp0*yp0*dxdy1*(g1-g2-g3+g4) + xp0*zp0*dxdz1*(g1-g2-g5+g6) + &
          yp0*zp0*dydz1*(g1-g3-g5+g7) + &
          xp0*yp0*zp0*dxdydz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
!
!  Do a reality check on the interpolation scheme.
!
      if (lcheck) then
        do i=1,ivar2-ivar1+1
          if (fp(i)>max(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
            print*, 'linear_interpolate_cartesian: interpolated value is LARGER than'
            print*, 'linear_interpolate_cartesian: a values at the corner points!'
            print*, 'linear_interpolate_cartesian: x0, y0, z0=', &
                xglobal(ix0), yglobal(iy0), zglobal(iz0)
            print*, 'linear_interpolate_cartesian: i, fp(i)=', i, fp(i)
            print*, 'linear_interpolate_cartesian: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
          endif
          if (fp(i)<min(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
            print*, 'linear_interpolate_cartesian: interpolated value is smaller than'
            print*, 'linear_interpolate_cartesian: a values at the corner points!'
            print*, 'linear_interpolate_cartesian: xxp=', xxp
            print*, 'linear_interpolate_cartesian: x0, y0, z0=', &
                xglobal(ix0), yglobal(iy0), zglobal(iz0)
            print*, 'linear_interpolate_cartesian: i, fp(i)=', i, fp(i)
            print*, 'linear_interpolate_cartesian: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
          endif
          if (fp(i)/=fp(i)) then
            print*, 'linear_interpolate_cartesian: interpolated value is NaN'
            print*, 'linear_interpolate_cartesian: xxp=', xxp
            print*, 'linear_interpolate_cartesian: x0, y0, z0=', &
                xglobal(ix0), yglobal(iy0), zglobal(iz0)
            print*, 'linear_interpolate_cartesian: i, fp(i)=', i, fp(i)
            print*, 'linear_interpolate_cartesian: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
            linear_interpolate_cartesian=.false.
          endif
        enddo
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
  endfunction linear_interpolate_cartesian
!***********************************************************************
  logical function linear_interpolate_curv_HO(farr,ivar1,ivar2,xxp,inear_glob,fp,lcheck,order)
!
!  Interpolate the value of f to arbitrary (xp, yp) CURVILINEAR coordinate
!  using the high-order lagrangian interpolation.
! 
!  TODO: Extend to 3D
!
!  The coefficients are determined by the 2xN grid points surrounding the
!  interpolation point.
!  Global coordinates are used for the interpolation, to allow interpolation of 
!  values outside this processors domain.
!
!  26-apr-17/Jorgen: Adapted from linear_interpolate_curvilinear
!
      integer :: ivar1, ivar2
      integer, intent(in) :: order
      real, dimension (3) :: xxp
      real, dimension (order,order,2,ivar2-ivar1+1) :: farr
      real, dimension (ivar2-ivar1+1) :: fp
      integer, dimension (3) :: inear_glob
      logical :: lcheck
!
      real :: xp0, yp0, zp0
!
      intent(in)  :: farr, ivar1, ivar2, xxp, inear_glob, lcheck
      intent(out) :: fp
      integer :: i,ix0,iy0,iz0
      real, dimension(ivar2-ivar1+1) :: g0   ,g1   ,g2   ,g3   
      real :: dx01 ,dx02 ,dx03 ,dx10 ,dx12 ,dx13 ,dx20 ,dx21 ,dx23 ,dx30 ,dx31 ,dx32 ,xp , &
              l0,l1,l2,l3,l4,l5,l6,l7 ,dy01 ,dy02 ,dy03 ,dy10 , &
              dy12 ,dy13 ,dy20 ,dy21 ,dy23 ,dy30 ,dy31 ,dy32,yp
      real, dimension(order,ivar2-ivar1+1) :: g_interp
      real, dimension(order) :: lagrange
      real, dimension(order,order) :: dx1,dy1
      integer :: j,k,l
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
      linear_interpolate_curv_HO= .true.
!
      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3)
!
!  Check if the grid point interval is really correct.
!
      if ((xglobal_ogrid(ix0)<=xxp(1) .and. xglobal_ogrid(ix0+1)>=xxp(1) .or. nxgrid==1) .and. &
          (yglobal_ogrid(iy0)<=xxp(2) .and. yglobal_ogrid(iy0+1)>=xxp(2) .or. nygrid==1) .and. &
          (zglobal_ogrid(iz0)<=xxp(3) .and. zglobal_ogrid(iz0+1)>=xxp(3) .or. nzgrid==1)) then
        ! Everything okay
      else
        print*, 'linear_interpolate_curvilinear_highorder: Global interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc_world
        print*, 'mxgrid_ogrid, xglobal_ogrid(1), xglobal_ogrid(mxgrid_ogrid) = ', & 
            mxgrid_ogrid, xglobal_ogrid(1), xglobal_ogrid(mxgrid_ogrid)
        print*, 'mygrid_ogrid, yglobal_ogrid(1), yglobal_ogrid(mygrid_ogrid) = ', &
            mygrid_ogrid, yglobal_ogrid(1), yglobal_ogrid(mygrid_ogrid)
        print*, 'mzgrid_ogrid, zglobal_ogrid(1), zglobal_ogrid(mzgrid_ogrid) = ', & 
            mzgrid_ogrid, zglobal_ogrid(1), zglobal_ogrid(mzgrid_ogrid)
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
        linear_interpolate_curv_HO= .false.
        return
      endif
!
!  Set up 1D Lagrange basis polynomials in x-direction
! 
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        l=-floor(order/2.)
        do j=1,order
          l=l+1
          dx1(i,j)=xglobal_ogrid(ix0+k)-xglobal_ogrid(ix0+l)
        enddo
        dx1(i,i)=1 ! To avoid division by zero
      enddo
      dx1=1./dx1
      xp=xxp(1)
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        lagrange(i)=1./(xp-xglobal_ogrid(ix0+k))
      enddo
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        lagrange(:)=lagrange(:)*(xp-xglobal_ogrid(ix0+k))*dx1(:,i)
      enddo
      g_interp=0
      do i=1,order
        g_interp(:,ivar1:ivar2)=g_interp(:,ivar1:ivar2)+farr(i,:,1,ivar1:ivar2)*lagrange(i)
      enddo
      ! y-dir
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        l=-floor(order/2.)
        do j=1,order
          l=l+1
          dy1(i,j)=yglobal_ogrid(iy0+k)-yglobal_ogrid(iy0+l)
        enddo
        dy1(i,i)=1 ! To avoid division by zero
      enddo
      dy1=1./dy1
      yp=xxp(2)
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        lagrange(i)=1./(yp-yglobal_ogrid(iy0+k))
      enddo
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        lagrange(:)=lagrange(:)*(yp-yglobal_ogrid(iy0+k))*dy1(:,i)
      enddo
      fp=0
      do i=1,order
        fp(ivar1:ivar2)=fp(ivar1:ivar2)+g_interp(i,ivar1:ivar2)*lagrange(i)
      enddo

!!       dx01=1/(xglobal_ogrid(ix0-1)-xglobal_ogrid(ix0))
!!       dx02=1/(xglobal_ogrid(ix0-1)-xglobal_ogrid(ix0+1))
!!       dx03=1/(xglobal_ogrid(ix0-1)-xglobal_ogrid(ix0+2))
!!       dx10=-dx01
!!       dx12=1/(xglobal_ogrid(ix0)-xglobal_ogrid(ix0+1))
!!       dx13=1/(xglobal_ogrid(ix0)-xglobal_ogrid(ix0+2))
!!       dx20=-dx02
!!       dx21=-dx12
!!       dx23=1/(xglobal_ogrid(ix0+1)-xglobal_ogrid(ix0+2))
!!       dx30=-dx03
!!       dx31=-dx12
!!       dx32=-dx23
!! !
!!       xp=xxp(1)
!! !
        
!!       l0=(xp-xglobal_ogrid(ix0  ))*dx01*(xp-xglobal_ogrid(ix0+1))*dx02!*(xp-xglobal_ogrid(ix0+2))*dx03 
!!       l1=(xp-xglobal_ogrid(ix0-1))*dx10*(xp-xglobal_ogrid(ix0+1))*dx12!*(xp-xglobal_ogrid(ix0+2))*dx03 
!!       l2=(xp-xglobal_ogrid(ix0-1))*dx20*(xp-xglobal_ogrid(ix0  ))*dx21!*(xp-xglobal_ogrid(ix0+2))*dx23 
!!       l3=(xp-xglobal_ogrid(ix0-1))*dx30*(xp-xglobal_ogrid(ix0  ))*dx31*(xp-xglobal_ogrid(ix0+1))*dx32 
!! !
!! !  Compute interpolation in x-direction (radial direction)
!! !
!!       g0=farr(1,1,1,ivar1:ivar2)*l0+farr(2,1,1,ivar1:ivar2)*l1+farr(3,1,1,ivar1:ivar2)*l2!+farr(4,1,1,ivar1:ivar2)*l3
!!       g1=farr(1,2,1,ivar1:ivar2)*l0+farr(2,2,1,ivar1:ivar2)*l1+farr(3,2,1,ivar1:ivar2)*l2!+farr(4,2,1,ivar1:ivar2)*l3
!!       g2=farr(1,3,1,ivar1:ivar2)*l0+farr(2,3,1,ivar1:ivar2)*l1+farr(3,3,1,ivar1:ivar2)*l2!+farr(4,3,1,ivar1:ivar2)*l3
!!       g3=farr(1,4,1,ivar1:ivar2)*l0+farr(2,4,1,ivar1:ivar2)*l1+farr(3,4,1,ivar1:ivar2)*l2+farr(4,4,1,ivar1:ivar2)*l3
!! !
!! !  Set up 1D Lagrange basis polynomials in y-direction
!! ! 
!!       dy01=1/(yglobal_ogrid(iy0-1)-yglobal_ogrid(iy0))
!!       dy02=1/(yglobal_ogrid(iy0-1)-yglobal_ogrid(iy0+1))
!!       dy03=1/(yglobal_ogrid(iy0-1)-yglobal_ogrid(iy0+2))
!!       dy10=-dy01
!!       dy12=1/(yglobal_ogrid(iy0)-yglobal_ogrid(iy0+1))
!!       dy13=1/(yglobal_ogrid(iy0)-yglobal_ogrid(iy0+2))
!!       dy20=-dy02
!!       dy21=-dy12
!!       dy23=1/(yglobal_ogrid(iy0+1)-yglobal_ogrid(iy0+2))
!!       dy30=-dy03
!!       dy31=-dy12
!!       dy32=-dy23
!! !
!!       yp=xxp(2)
!! !
!!       l4=(yp-yglobal_ogrid(iy0  ))*dy01*(yp-yglobal_ogrid(iy0+1))*dy02!*(yp-yglobal_ogrid(iy0+2))*dy03 
!!       l5=(yp-yglobal_ogrid(iy0-1))*dy10*(yp-yglobal_ogrid(iy0+1))*dy12!*(yp-yglobal_ogrid(iy0+2))*dy03 
!!       l6=(yp-yglobal_ogrid(iy0-1))*dy20*(yp-yglobal_ogrid(iy0  ))*dy21!*(yp-yglobal_ogrid(iy0+2))*dy23 
!!       l7=(yp-yglobal_ogrid(iy0-1))*dy30*(yp-yglobal_ogrid(iy0  ))*dy31*(yp-yglobal_ogrid(iy0+1))*dy32 
!! !
!! !  Compute interpolation in y-direction (radial direction)
!! !
!!       fp=g0*l4+g1*l5+g2*l6!+g3*l7
!! !
!! !  Do a reality check on the interpolation scheme.
!
      if (lcheck) then
        do i=1,ivar2-ivar1+1
          !if (fp(i)>maxval(farr(:,:,1,i)).and.i/=3) then
            !print*, 'linear_interpolate_curvilinear: interpolated value is LARGER than'
            !print*, 'linear_interpolate_curvilinear: a values at the corner points!'
            !print*, 'linear_interpolate_curvilinear: xp, yp =',xxp(1:2)
            !print*, 'linear_interpolate_curvilinear: x0, y0 =', &
            !    xglobal_ogrid(ix0), yglobal_ogrid(iy0)
            !print*, 'linear_interpolate_curvilinear: x0+1, y0+1=', &
            !    xglobal_ogrid(ix0+1), yglobal_ogrid(iy0+1)
            !print*, 'linear_interpolate_curvilinear: i, fp(i)=', i, fp(i)
            !print*, 'linear_interpolate_curvilinear: farr(1:order,1)=', farr(1:order,1,1,i)
            !print*, 'linear_interpolate_curvilinear: farr(1:order,2)=', farr(1:order,2,1,i)
            !print*, 'linear_interpolate_curvilinear: farr(1:order,3)=', farr(1:order,3,1,i)
            !print*, 'linear_interpolate_curvilinear: farr(1:order,4)=', farr(1:order,4,1,i)
            !print*, 'linear_interpolate_curvilinear: g0,g1,g2,g3=', g0(i),g1(i),g2(i),g3(i)
            !print*, 'linear_interpolate_curvilinear: l0,l1,l2,l3=', l0,l1,l2,l3
            !print*, '------------------'
          !endif
          !if (fp(i)<minval(farr(:,:,1,i)).and.i/=3) then
            !print*, 'linear_interpolate_curvilinear: interpolated value is smaller than'
            !print*, 'linear_interpolate_curvilinear: a values at the corner points!'
            !print*, 'linear_interpolate_curvilinear: xp, yp =',xxp(1:2)
            !print*, 'linear_interpolate_curvilinear: x0, y0 =', &
            !    xglobal_ogrid(ix0), yglobal_ogrid(iy0)
            !print*, 'linear_interpolate_curvilinear: x0+1, y0+1=', &
            !    xglobal_ogrid(ix0+1), yglobal_ogrid(iy0+1)
            !print*, 'linear_interpolate_curvilinear: i, fp(i)=', i, fp(i)
            !print*, 'linear_interpolate_curvilinear: farr(1:order,1)=', farr(1:order,1,1,i)
            !print*, 'linear_interpolate_curvilinear: farr(1:order,2)=', farr(1:order,2,1,i)
            !print*, 'linear_interpolate_curvilinear: farr(1:order,3)=', farr(1:order,3,1,i)
            !print*, 'linear_interpolate_curvilinear: farr(1:order,4)=', farr(1:order,4,1,i)
            !print*, 'linear_interpolate_curvilinear: g0,g1,g2,g3=', g0(i),g1(i),g2(i),g3(i)
            !print*, 'linear_interpolate_curvilinear: l0,l1,l2,l3=', l0,l1,l2,l3
            !print*, '------------------'
          !endif
          if (fp(i)>maxval(farr(:,:,1,i)).and.i/=3) then 
           l1=(xp-xglobal_ogrid(ix0+1))/(xglobal_ogrid(ix0)-xglobal_ogrid(ix0+1))
           l2=(xp-xglobal_ogrid(ix0  ))/(xglobal_ogrid(ix0+1)-xglobal_ogrid(ix0))
           g1=farr(k,k,1,ivar1:ivar2)*l1+farr(k+1,k,1,ivar1:ivar2)*l2
           g2=farr(k,k+1,1,ivar1:ivar2)*l1+farr(k,k+1,1,ivar1:ivar2)*l2
           l3=(yp-yglobal_ogrid(iy0+1))/(yglobal_ogrid(iy0)-yglobal_ogrid(iy0+1))
           l4=(yp-yglobal_ogrid(iy0  ))/(yglobal_ogrid(iy0+1)-yglobal_ogrid(iy0))
           fp=g1*l3+g2*l4
          elseif (fp(i)<minval(farr(:,:,1,i)).and.i/=3) then
           l1=(xp-xglobal_ogrid(ix0+1))/(xglobal_ogrid(ix0)-xglobal_ogrid(ix0+1))
           l2=(xp-xglobal_ogrid(ix0  ))/(xglobal_ogrid(ix0+1)-xglobal_ogrid(ix0))
           g1=farr(k,k,1,ivar1:ivar2)*l1+farr(k+1,k,1,ivar1:ivar2)*l2
           g2=farr(k,k+1,1,ivar1:ivar2)*l1+farr(k,k+1,1,ivar1:ivar2)*l2
           l3=(yp-yglobal_ogrid(iy0+1))/(yglobal_ogrid(iy0)-yglobal_ogrid(iy0+1))
           l4=(yp-yglobal_ogrid(iy0  ))/(yglobal_ogrid(iy0+1)-yglobal_ogrid(iy0))
           fp=g1*l3+g2*l4
          endif
          if ((fp(i)>maxval(farr(:,:,1,i)).and.i/=3) .or. (fp(i)<minval(farr(:,:,1,i)).and.i/=3)) then
            print*, 'linear_interpolate_curvilinear_HO: interpolated value is smaller or larger than'
            print*, 'linear_interpolate_curvilinear_HO: a values at the corner points, even after linearization!'
            print*, '------------------'
            linear_interpolate_curv_HO=.false.
          endif
          if (fp(i)/=fp(i)) then
            print*, 'linear_interpolate_curvilinear: interpolated value is NaN'
            print*, 'linear_interpolate_curvilinear: xxp=', xxp
            print*, 'linear_interpolate_curvilinear: x0, y0, z0=', &
                xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
            print*, 'linear_interpolate_curvilinear: i, fp(i)=', i, fp(i)
            print*, 'linear_interpolate_curvilinear: farr=', farr(:,:,1,i)
            print*, '------------------'
            linear_interpolate_curv_HO=.false.
          endif
        enddo
      endif
!
  endfunction linear_interpolate_curv_HO
!***********************************************************************
  logical function linear_interpolate_cart_HO(farr,ivar1,ivar2,xxp,inear_glob,fp,lcheck,order)
!
!  Interpolate the value of f to arbitrary (xp, yp) CARTESIAN coordinate
!  using the high-order lagrangian interpolation.
! 
!  TODO: Extend to 3D
!  TODO: Extend to arbitrary order
!
!  The coefficients are determined by the 2xN grid points surrounding the
!  interpolation point.
!  Global coordinates are used for the interpolation, to allow interpolation of 
!  values outside this processors domain.
!
!  26-apr-17/Jorgen: Adapted from linear_interpolate_curv_HO
!
      integer :: ivar1, ivar2
      integer, intent(in) :: order
      real, dimension (3) :: xxp
      real, dimension (order,order,2,ivar2-ivar1+1) :: farr
      real, dimension (ivar2-ivar1+1) :: fp
      integer, dimension (3) :: inear_glob
      logical :: lcheck
!
      real :: xp0, yp0, zp0
!
      intent(in)  :: farr, ivar1, ivar2, xxp, inear_glob, lcheck
      intent(out) :: fp
      integer :: i,ix0,iy0,iz0
      real, dimension(ivar2-ivar1+1) :: g0   ,g1   ,g2   ,g3   
      real :: dx01 ,dx02 ,dx03 ,dx10 ,dx12 ,dx13 ,dx20 ,dx21 ,dx23 ,dx30 ,dx31 ,dx32 ,xp , &
              l0,l1,l2,l3,l4,l5,l6,l7 ,dy01 ,dy02 ,dy03 ,dy10 , &
              dy12 ,dy13 ,dy20 ,dy21 ,dy23 ,dy30 ,dy31 ,dy32,yp
      real, dimension(order,ivar2-ivar1+1) :: g_interp
      real, dimension(order) :: lagrange
      real, dimension(order,order) :: dx1,dy1
      integer :: j,k,l
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
      linear_interpolate_cart_HO= .true.
!
      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3)
!
!  Check if the grid point interval is really correct.
!
      if ((xglobal(ix0)<=xxp(1) .and. xglobal(ix0+1)>=xxp(1) .or. nxgrid==1) .and. &
          (yglobal(iy0)<=xxp(2) .and. yglobal(iy0+1)>=xxp(2) .or. nygrid==1) .and. &
          (zglobal(iz0)<=xxp(3) .and. zglobal(iz0+1)>=xxp(3) .or. nzgrid==1)) then
        ! Everything okay
      else
        print*, 'linear_interpolate_cartesian: Global interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc_world
        print*, 'mxgrid, xglobal(1), xglobal(mx) = ', mxgrid, xglobal(1), xglobal(mxgrid)
        print*, 'mygrid, yglobal(1), yglobal(my) = ', mygrid, yglobal(1), yglobal(mygrid)
        print*, 'mzgrid, zglobal(1), zglobal(mz) = ', mzgrid, zglobal(1), zglobal(mzgrid)
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal(ix0), xglobal(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal(iy0), yglobal(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal(iz0), zglobal(iz0+1)
        linear_interpolate_cart_HO = .false.
        return
      endif
!
!  Set up 1D Lagrange basis polynomials in x-direction
! 
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        l=-floor(order/2.)
        do j=1,order
          l=l+1
          dx1(i,j)=xglobal(ix0+k)-xglobal(ix0+l)
        enddo
        dx1(i,i)=1 ! To avoid division by zero
      enddo
      dx1=1./dx1
      xp=xxp(1)
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        lagrange(i)=1./(xp-xglobal(ix0+k))
      enddo
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        lagrange(:)=lagrange(:)*(xp-xglobal(ix0+k))*dx1(:,i)
      enddo
      g_interp=0
      do i=1,order
        g_interp(:,ivar1:ivar2)=g_interp(:,ivar1:ivar2)+farr(i,:,1,ivar1:ivar2)*lagrange(i)
      enddo
      ! y-dir
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        l=-floor(order/2.)
        do j=1,order
          l=l+1
          dy1(i,j)=yglobal(iy0+k)-yglobal(iy0+l)
        enddo
        dy1(i,i)=1 ! To avoid division by zero
      enddo
      dy1=1./dy1
      yp=xxp(2)
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        lagrange(i)=1./(yp-yglobal(iy0+k))
      enddo
      k=-floor(order/2.)
      do i=1,order
        k=k+1
        lagrange(:)=lagrange(:)*(yp-yglobal(iy0+k))*dy1(:,i)
      enddo
      fp=0
      do i=1,order
        fp(ivar1:ivar2)=fp(ivar1:ivar2)+g_interp(i,ivar1:ivar2)*lagrange(i)
      enddo
!
!  Do a reality check on the interpolation scheme.
!
      if (lcheck) then
        do i=1,ivar2-ivar1+1
          !if (fp(i)>maxval(farr(:,:,1,i)).and.i/=3) then
          !  print*, 'linear_interpolate_cart_HO: interpolated value is LARGER than'
          !  print*, 'linear_interpolate_cart_HO: a values at the corner points!'
          !  print*, 'linear_interpolate_cart_HO: xp, yp =',xxp(1:2)
          !  print*, 'linear_interpolate_cart_HO: x0, y0 =', &
          !      xglobal(ix0), yglobal(iy0)
          !  print*, 'linear_interpolate_cart_HO: x0+1, y0+1=', &
          !      xglobal(ix0+1), yglobal(iy0+1)
          !  print*, 'linear_interpolate_cart_HO: i, fp(i)=', i, fp(i)
          !  print*, 'linear_interpolate_cart_HO: farr(1:order,1)=', farr(1:order,1,1,i)
          !  print*, 'linear_interpolate_cart_HO: farr(1:order,2)=', farr(1:order,2,1,i)
          !  print*, 'linear_interpolate_cart_HO: farr(1:order,3)=', farr(1:order,3,1,i)
          !  print*, 'linear_interpolate_cart_HO: farr(1:order,4)=', farr(1:order,4,1,i)
          !  print*, '------------------'
          !endif
          !if (fp(i)<minval(farr(:,:,1,i)).and.i/=3) then
          !  print*, 'linear_interpolate_cart_HO: interpolated value is smaller than'
          !  print*, 'linear_interpolate_cart_HO: a values at the corner points!'
          !  print*, 'linear_interpolate_cart_HO: xp, yp =',xxp(1:2)
          !  print*, 'linear_interpolate_cart_HO: x0, y0 =', &
          !      xglobal(ix0), yglobal(iy0)
          !  print*, 'linear_interpolate_cart_HO: x0+1, y0+1=', &
          !      xglobal(ix0+1), yglobal(iy0+1)
          !  print*, 'linear_interpolate_cart_HO: i, fp(i)=', i, fp(i)
          !  print*, 'linear_interpolate_cart_HO: farr(1:order,1)=', farr(1:order,1,1,i)
          !  print*, 'linear_interpolate_cart_HO: farr(1:order,2)=', farr(1:order,2,1,i)
          !  print*, 'linear_interpolate_cart_HO: farr(1:order,3)=', farr(1:order,3,1,i)
          !  print*, 'linear_interpolate_cart_HO: farr(1:order,4)=', farr(1:order,4,1,i)
          !  print*, '------------------'
          !endif
          if (fp(i)>maxval(farr(:,:,1,i)).and.i/=3) then 
           l1=(xp-xglobal(ix0+1))/(xglobal(ix0)-xglobal(ix0+1))
           l2=(xp-xglobal(ix0  ))/(xglobal(ix0+1)-xglobal(ix0))
           g1=farr(k,k,1,ivar1:ivar2)*l1+farr(k+1,k,1,ivar1:ivar2)*l2
           g2=farr(k,k+1,1,ivar1:ivar2)*l1+farr(k,k+1,1,ivar1:ivar2)*l2
           l3=(yp-yglobal(iy0+1))/(yglobal(iy0)-yglobal(iy0+1))
           l4=(yp-yglobal(iy0  ))/(yglobal(iy0+1)-yglobal(iy0))
           fp=g1*l3+g2*l4
          elseif (fp(i)<minval(farr(:,:,1,i)).and.i/=3) then
           l1=(xp-xglobal(ix0+1))/(xglobal(ix0)-xglobal(ix0+1))
           l2=(xp-xglobal(ix0  ))/(xglobal(ix0+1)-xglobal(ix0))
           g1=farr(k,k,1,ivar1:ivar2)*l1+farr(k+1,k,1,ivar1:ivar2)*l2
           g2=farr(k,k+1,1,ivar1:ivar2)*l1+farr(k,k+1,1,ivar1:ivar2)*l2
           l3=(yp-yglobal(iy0+1))/(yglobal(iy0)-yglobal(iy0+1))
           l4=(yp-yglobal(iy0  ))/(yglobal(iy0+1)-yglobal(iy0))
           fp=g1*l3+g2*l4
          endif
          if ((fp(i)>maxval(farr(:,:,1,i)).and.i/=3) .or. (fp(i)<minval(farr(:,:,1,i)).and.i/=3)) then
            print*, 'linear_interpolate_cart_HO: interpolated value is smaller or larger than'
            print*, 'linear_interpolate_cart_HO: a values at the corner points, even after linearization!'
            print*, '------------------'
            linear_interpolate_cart_HO=.false.
          endif
          if (fp(i)/=fp(i)) then
            print*, 'linear_interpolate_cart_HO: interpolated value is NaN'
            print*, 'linear_interpolate_cart_HO: xxp=', xxp
            print*, 'linear_interpolate_cart_HO: x0, y0, z0=', &
                xglobal(ix0), yglobal(iy0), zglobal(iz0)
            print*, 'linear_interpolate_cart_HO: i, fp(i)=', i, fp(i)
            print*, 'linear_interpolate_cart_HO: farr=', farr(:,:,1,i)
            print*, '------------------'
            linear_interpolate_cart_HO=.false.
          endif
        enddo
      endif
!
  endfunction linear_interpolate_cart_HO
!***********************************************************************
  logical function hermite_interpolate_curv(farr,ivar1,ivar2,xxp,inear_glob,inear_loc,fp,lcheck) 
!
!  Interpolate between 2D-grids using Hermite surface interpolation
!
    integer :: ivar1, ivar2
    real, dimension(3) :: xxp
    real, dimension(2,2,2,ivar2-ivar1+1) :: farr
    real, dimension(ivar2-ivar1+1) :: fp
    integer, dimension(3) :: inear_glob, inear_loc
    logical :: lcheck
!    
    intent(in)  :: farr, ivar1, ivar2, xxp, inear_glob, lcheck
    intent(out) :: fp
!
    real, dimension(4) :: Fu, Fv
    real, dimension(4,4,ivar2-ivar1+1) :: B
    real :: u,v
!
    integer :: ix0,iy0,iz0,ix0_loc,iy0_loc,iz0_loc,i
    real, parameter :: a = 1.0 / 60.0
    real :: fac
!
    ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3)
    ix0_loc=inear_loc(1); iy0_loc=inear_loc(2); iz0_loc=inear_loc(3)
    u = (xxp(1)-xglobal_ogrid(ix0))/(xglobal_ogrid(ix0+1)-xglobal_ogrid(ix0))
    v = (xxp(2)-yglobal_ogrid(iy0))/(yglobal_ogrid(iy0+1)-yglobal_ogrid(iy0))
    Fu(1) =  2*u**3 - 3*u*u +1
    Fu(2) = -2*u**3 + 3*u*u 
    Fu(3) =    u**3 - 2*u*u +u
    Fu(4) =    u**3 - u**2
    Fv(1) =  2*v**3 - 3*v*v +1
    Fv(2) = -2*v**3 + 3*v*v 
    Fv(3) =    v**3 - 2*v*v +v
    Fv(4) =    v**3 - v**2
    B(1,1,:) = farr(1,1,1,:)
    B(2,1,:) = farr(2,1,1,:)
    B(1,2,:) = farr(1,2,1,:)
    B(2,2,:) = farr(2,2,1,:)
!
    fac = a
    B(3,1,:) = fac*(+ 45.0*(f_ogrid(ix0_loc+1  ,iy0_loc  ,iz0_loc,:)-f_ogrid(ix0_loc-1  ,iy0_loc  ,iz0_loc,:)) &
                    -  9.0*(f_ogrid(ix0_loc+2  ,iy0_loc  ,iz0_loc,:)-f_ogrid(ix0_loc-2  ,iy0_loc  ,iz0_loc,:)) &
                    +      (f_ogrid(ix0_loc+3  ,iy0_loc  ,iz0_loc,:)-f_ogrid(ix0_loc-3  ,iy0_loc  ,iz0_loc,:)))
    B(4,1,:) = fac*(+ 45.0*(f_ogrid(ix0_loc+1+1,iy0_loc  ,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc  ,iz0_loc,:)) &
                    -  9.0*(f_ogrid(ix0_loc+2+1,iy0_loc  ,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc  ,iz0_loc,:)) &
                    +      (f_ogrid(ix0_loc+3+1,iy0_loc  ,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc  ,iz0_loc,:)))
    B(3,2,:) = fac*(+ 45.0*(f_ogrid(ix0_loc+1  ,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-1  ,iy0_loc+1,iz0_loc,:)) &
                    -  9.0*(f_ogrid(ix0_loc+2  ,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-2  ,iy0_loc+1,iz0_loc,:)) &
                    +      (f_ogrid(ix0_loc+3  ,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-3  ,iy0_loc+1,iz0_loc,:)))
    B(4,2,:) = fac*(+ 45.0*(f_ogrid(ix0_loc+1+1,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc+1,iz0_loc,:)) &
                    -  9.0*(f_ogrid(ix0_loc+2+1,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc+1,iz0_loc,:)) &
                    +      (f_ogrid(ix0_loc+3+1,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc+1,iz0_loc,:)))
!
    fac = a*rcyl_mn1_ogrid(ix0_loc)
    B(1,3,:) = fac*(+ 45.0*(f_ogrid(ix0_loc  ,iy0_loc+1  ,iz0_loc,:)-f_ogrid(ix0_loc  ,iy0_loc-1  ,iz0_loc,:)) &
                    -  9.0*(f_ogrid(ix0_loc  ,iy0_loc+2  ,iz0_loc,:)-f_ogrid(ix0_loc  ,iy0_loc-2  ,iz0_loc,:)) &
                    +      (f_ogrid(ix0_loc  ,iy0_loc+3  ,iz0_loc,:)-f_ogrid(ix0_loc  ,iy0_loc-3  ,iz0_loc,:)))
    B(1,4,:) = fac*(+ 45.0*(f_ogrid(ix0_loc  ,iy0_loc+1+1,iz0_loc,:)-f_ogrid(ix0_loc  ,iy0_loc-1+1,iz0_loc,:)) &
                    -  9.0*(f_ogrid(ix0_loc  ,iy0_loc+2+1,iz0_loc,:)-f_ogrid(ix0_loc  ,iy0_loc-2+1,iz0_loc,:)) &
                    +      (f_ogrid(ix0_loc  ,iy0_loc+3+1,iz0_loc,:)-f_ogrid(ix0_loc  ,iy0_loc-3+1,iz0_loc,:)))
    B(2,3,:) = fac*(+ 45.0*(f_ogrid(ix0_loc+1,iy0_loc+1  ,iz0_loc,:)-f_ogrid(ix0_loc+1,iy0_loc-1  ,iz0_loc,:)) &
                    -  9.0*(f_ogrid(ix0_loc+1,iy0_loc+2  ,iz0_loc,:)-f_ogrid(ix0_loc+1,iy0_loc-2  ,iz0_loc,:)) &
                    +      (f_ogrid(ix0_loc+1,iy0_loc+3  ,iz0_loc,:)-f_ogrid(ix0_loc+1,iy0_loc-3  ,iz0_loc,:)))
    B(2,4,:) = fac*(+ 45.0*(f_ogrid(ix0_loc+1,iy0_loc+1+1,iz0_loc,:)-f_ogrid(ix0_loc+1,iy0_loc-1+1,iz0_loc,:)) &
                    -  9.0*(f_ogrid(ix0_loc+1,iy0_loc+2+1,iz0_loc,:)-f_ogrid(ix0_loc+1,iy0_loc-2+1,iz0_loc,:)) &
                    +      (f_ogrid(ix0_loc+1,iy0_loc+3+1,iz0_loc,:)-f_ogrid(ix0_loc+1,iy0_loc-3+1,iz0_loc,:)))
!  
    fac=(1./60.**2)*rcyl_mn1_ogrid(ix0_loc)
    B(3,3,:)=fac*( &
              45.*((45.*(f_ogrid(ix0_loc+1,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc+1,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1,iy0_loc-1,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc-1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc-1,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc-1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc-1,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc-1,iz0_loc,:))))&
              -9.*((45.*(f_ogrid(ix0_loc+1,iy0_loc+2,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc+2,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc+2,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc+2,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc+2,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc+2,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1,iy0_loc-2,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc-2,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc-2,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc-2,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc-2,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc-2,iz0_loc,:))))&
                 +((45.*(f_ogrid(ix0_loc+1,iy0_loc+3,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc+3,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc+3,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc+3,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc+3,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc+3,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1,iy0_loc-3,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc-3,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc-3,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc-3,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc-3,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc-3,iz0_loc,:))))&
                   )
    B(4,3,:)=fac*( &
              45.*((45.*(f_ogrid(ix0_loc+1+1,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc+1,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc+1,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1+1,iy0_loc-1,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc-1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc-1,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc-1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc-1,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc-1,iz0_loc,:))))&
              -9.*((45.*(f_ogrid(ix0_loc+1+1,iy0_loc+2,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc+2,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc+2,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc+2,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc+2,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc+2,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1+1,iy0_loc-2,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc-2,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc-2,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc-2,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc-2,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc-2,iz0_loc,:))))&
                 +((45.*(f_ogrid(ix0_loc+1+1,iy0_loc+3,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc+3,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc+3,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc+3,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc+3,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc+3,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1+1,iy0_loc-3,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc-3,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc-3,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc-3,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc-3,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc-3,iz0_loc,:))))&
                   )
    B(3,4,:)=fac*( &
              45.*((45.*(f_ogrid(ix0_loc+1,iy0_loc+1+1,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc+1+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc+1+1,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc+1+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc+1+1,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc+1+1,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1,iy0_loc-1+1,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc-1+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc-1+1,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc-1+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc-1+1,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc-1+1,iz0_loc,:))))&
              -9.*((45.*(f_ogrid(ix0_loc+1,iy0_loc+2+1,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc+2+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc+2+1,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc+2+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc+2+1,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc+2+1,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1,iy0_loc-2+1,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc-2+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc-2+1,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc-2+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc-2+1,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc-2+1,iz0_loc,:))))&
                 +((45.*(f_ogrid(ix0_loc+1,iy0_loc+3+1,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc+3+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc+3+1,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc+3+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc+3+1,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc+3+1,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1,iy0_loc-3+1,iz0_loc,:)-f_ogrid(ix0_loc-1,iy0_loc-3+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2,iy0_loc-3+1,iz0_loc,:)-f_ogrid(ix0_loc-2,iy0_loc-3+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3,iy0_loc-3+1,iz0_loc,:)-f_ogrid(ix0_loc-3,iy0_loc-3+1,iz0_loc,:))))&
                   )
    B(4,4,:)=fac*( &
              45.*((45.*(f_ogrid(ix0_loc+1+1,iy0_loc+1+1,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc+1+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc+1+1,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc+1+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc+1+1,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc+1+1,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1+1,iy0_loc-1+1,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc-1+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc-1+1,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc-1+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc-1+1,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc-1+1,iz0_loc,:))))&
              -9.*((45.*(f_ogrid(ix0_loc+1+1,iy0_loc+2+1,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc+2+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc+2+1,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc+2+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc+2+1,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc+2+1,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1+1,iy0_loc-2+1,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc-2+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc-2+1,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc-2+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc-2+1,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc-2+1,iz0_loc,:))))&
                 +((45.*(f_ogrid(ix0_loc+1+1,iy0_loc+3+1,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc+3+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc+3+1,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc+3+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc+3+1,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc+3+1,iz0_loc,:))) &
                  -(45.*(f_ogrid(ix0_loc+1+1,iy0_loc-3+1,iz0_loc,:)-f_ogrid(ix0_loc-1+1,iy0_loc-3+1,iz0_loc,:))  &
                    -9.*(f_ogrid(ix0_loc+2+1,iy0_loc-3+1,iz0_loc,:)-f_ogrid(ix0_loc-2+1,iy0_loc-3+1,iz0_loc,:))  &
                       +(f_ogrid(ix0_loc+3+1,iy0_loc-3+1,iz0_loc,:)-f_ogrid(ix0_loc-3+1,iy0_loc-3+1,iz0_loc,:))))&
                   )
!
    do i=ivar1,ivar2
      fp(i)=dot_product(matmul(Fu,B(:,:,i)),Fv)
    enddo

    hermite_interpolate_curv=.true.
!
!  Do a reality check on the interpolation scheme.
!
    if (lcheck) then
      do i=1,ivar2-ivar1+1
        if (fp(i)>maxval(farr)) then
          print*, 'hermite_interpolate_curvilinear: interpolated value is LARGER than'
          print*, 'hermite_interpolate_curvilinear: a values at the corner points!'
          print*, 'hermite_interpolate_curvilinear: xxp=', xxp
          print*, 'hermite_interpolate_curvilinear: x0, y0, z0=', &
              xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
          print*, 'hermite_interpolate_curvilinear: x1, y1, z1=', &
              xglobal_ogrid(ix0+1), yglobal_ogrid(iy0+1), zglobal_ogrid(iz0+1)
          print*, 'hermite_interpolate_curvilinear: u, v =',u,v
          print*, 'hermite_interpolate_curvilinear: i, fp(i)=', i, fp(i)
          print*, 'hermite_interpolate_curvilinear: farr=', farr(:,:,1,i)
          print*, 'B-matrix', B(:,:,i)
          print*, '------------------'
        endif
        if (fp(i)<minval(farr)) then 
          print*, 'hermite_interpolate_curvilinear: interpolated value is smaller than'
          print*, 'hermite_interpolate_curvilinear: a values at the corner points!'
          print*, 'hermite_interpolate_curvilinear: xxp=', xxp
          print*, 'hermite_interpolate_curvilinear: x0, y0, z0=', &
              xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
          print*, 'hermite_interpolate_curvilinear: x1, y1, z1=', &
              xglobal_ogrid(ix0+1), yglobal_ogrid(iy0+1), zglobal_ogrid(iz0+1)
          print*, 'hermite_interpolate_curvilinear: u, v =',u,v
          print*, 'hermite_interpolate_curvilinear: i, fp(i)=', i, fp(i)
          print*, 'hermite_interpolate_curvilinear: farr=', farr(:,:,1,i)
          print*, 'B-matrix', B(:,:,i)
          print*, '------------------'
        endif
        if (fp(i)/=fp(i)) then
          print*, 'linear_interpolate_curvilinear: interpolated value is NaN'
          print*, 'linear_interpolate_curvilinear: xxp=', xxp
          print*, 'linear_interpolate_curvilinear: x0, y0, z0=', &
              xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
          print*, 'linear_interpolate_curvilinear: i, fp(i)=', i, fp(i)
          print*, 'hermite_interpolate_curvilinear: farr=', farr(:,:,1,i)
          print*, '------------------'
          hermite_interpolate_curv=.false.
        endif
      enddo
    endif
  
  
  endfunction hermite_interpolate_curv
!***********************************************************************
  logical function linear_interpolate_curvilinear(farr,ivar1,ivar2,xxp,inear_glob,fp,lcheck)
!
!  Interpolate the value of f to arbitrary (xp, yp, zp) CURVILINEAR coordinate
!  using the linear interpolation formula.
!
!    g(x,y,z) = A*x*y*z + B*x*y + C*x*z + D*y*z + E*x + F*y + G*z + H .
!
!  The coefficients are determined by the 8 grid points surrounding the
!  interpolation point.
!  Global coordinates are used for the interpolation, to allow interpolation of 
!  values outside this processors domain.
!
!  26-apr-17/Jorgen: Adapted from linear_interpolate_cartesian
!
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (2,2,2,ivar2-ivar1+1) :: farr
      real, dimension (ivar2-ivar1+1) :: fp
      integer, dimension (3) :: inear_glob
      logical :: lcheck
!
      real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, g5, g6, g7, g8
      real :: xp0, yp0, zp0
      real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
      logical :: lfirstcall=.true.
      integer :: ix0, iy0, iz0, i
!
      intent(in)  :: farr, ivar1, ivar2, xxp, inear_glob, lcheck
      intent(out) :: fp
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
      linear_interpolate_curvilinear= .true.
!
      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3)
!
!  Check if the grid point interval is really correct.
!
      if ((xglobal_ogrid(ix0)<=xxp(1) .and. xglobal_ogrid(ix0+1)>=xxp(1) .or. nxgrid==1) .and. &
          (yglobal_ogrid(iy0)<=xxp(2) .and. yglobal_ogrid(iy0+1)>=xxp(2) .or. nygrid==1) .and. &
          (zglobal_ogrid(iz0)<=xxp(3) .and. zglobal_ogrid(iz0+1)>=xxp(3) .or. nzgrid==1)) then
        ! Everything okay
      else
        print*, 'linear_interpolate_curvilinear: Global interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc_world
        print*, 'mxgrid_ogrid, xglobal_ogrid(1), xglobal_ogrid(mxgrid_ogrid) = ', & 
            mxgrid_ogrid, xglobal_ogrid(1), xglobal_ogrid(mxgrid_ogrid)
        print*, 'mygrid_ogrid, yglobal_ogrid(1), yglobal_ogrid(mygrid_ogrid) = ', &
            mygrid_ogrid, yglobal_ogrid(1), yglobal_ogrid(mygrid_ogrid)
        print*, 'mzgrid_ogrid, zglobal_ogrid(1), zglobal_ogrid(mzgrid_ogrid) = ', & 
            mzgrid_ogrid, zglobal_ogrid(1), zglobal_ogrid(mzgrid_ogrid)
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
        linear_interpolate_curvilinear= .false.
        return
      endif
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!  Set it equal to 0 for dimensions having 1 grid points; this will make sure
!  that the interpolation is bilinear for 2D grids.
!
      xp0=0; yp0=0; zp0=0
      if (nxgrid_ogrid/=1) xp0=xxp(1)-xglobal_ogrid(ix0)
      if (nygrid_ogrid/=1) yp0=xxp(2)-yglobal_ogrid(iy0)
      if (nzgrid_ogrid/=1) zp0=xxp(3)-zglobal_ogrid(iz0)
!
!  Calculate derived grid spacing parameters needed for interpolation.
!
      dx1=1/(xglobal_ogrid(ix0+1)-xglobal_ogrid(ix0))
      dy1=1/(yglobal_ogrid(iy0+1)-yglobal_ogrid(iy0))
      if(nzgrid/=1) then
        dz1=1/(zglobal_ogrid(iz0+1)-zglobal_ogrid(iz0))
      else 
        dz1=1
      endif
!
      dxdy1=dx1*dy1; dxdz1=dx1*dz1; dydz1=dy1*dz1
      dxdydz1=dx1*dy1*dz1
!
!  Function values at all corners.
!
      g1=farr(1,1,1,ivar1:ivar2)
      g2=farr(2,1,1,ivar1:ivar2)
      g3=farr(1,2,1,ivar1:ivar2)
      g4=farr(2,2,1,ivar1:ivar2)
      g5=farr(1,1,2,ivar1:ivar2)
      g6=farr(2,1,2,ivar1:ivar2)
      g7=farr(1,2,2,ivar1:ivar2)
      g8=farr(2,2,2,ivar1:ivar2)
!
!  Interpolation formula.
!
      fp = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) + &
          xp0*yp0*dxdy1*(g1-g2-g3+g4) + xp0*zp0*dxdz1*(g1-g2-g5+g6) + &
          yp0*zp0*dydz1*(g1-g3-g5+g7) + &
          xp0*yp0*zp0*dxdydz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
!
!  Do a reality check on the interpolation scheme.
!
      if (lcheck) then
        do i=1,ivar2-ivar1+1
          if (fp(i)>max(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
            print*, 'linear_interpolate_curvilinear: interpolated value is LARGER than'
            print*, 'linear_interpolate_curvilinear: a values at the corner points!'
            print*, 'linear_interpolate_curvilinear: x0, y0, z0=', &
                xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
            print*, 'linear_interpolate_curvilinear: i, fp(i)=', i, fp(i)
            print*, 'linear_interpolate_curvilinear: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
          endif
          if (fp(i)<min(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
            print*, 'linear_interpolate_curvilinear: interpolated value is smaller than'
            print*, 'linear_interpolate_curvilinear: a values at the corner points!'
            print*, 'linear_interpolate_curvilinear: xxp=', xxp
            print*, 'linear_interpolate_curvilinear: x0, y0, z0=', &
                xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
            print*, 'linear_interpolate_curvilinear: i, fp(i)=', i, fp(i)
            print*, 'linear_interpolate_curvilinear: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
          endif
          if (fp(i)/=fp(i)) then
            print*, 'linear_interpolate_curvilinear: interpolated value is NaN'
            print*, 'linear_interpolate_curvilinear: xxp=', xxp
            print*, 'linear_interpolate_curvilinear: x0, y0, z0=', &
                xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
            print*, 'linear_interpolate_curvilinear: i, fp(i)=', i, fp(i)
            print*, 'linear_interpolate_curvilinear: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
            linear_interpolate_curvilinear=.false.
          endif
        enddo
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
  endfunction linear_interpolate_curvilinear
!***********************************************************************
  subroutine solid_cells_timestep_first(f)
!
!  08-feb-17/Jorgen: Coded
!
!  Only save time, which is the time before timesteps on cartesian grid
!  is performed. Will need this to set timestep of ogrid iterations.
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

    real :: Lx_og,Ly_og,Lz_og
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
    Lx_og = Lxyz_ogrid(1)
    Ly_og = Lxyz_ogrid(2)
    Lz_og = Lxyz_ogrid(3)
!
!  Set the lower boundary and the grid size.
!
    x00 = xyz0_ogrid(1)
    y00 = xyz0_ogrid(2)
    z00 = xyz0_ogrid(3)
!
    dx_ogrid = Lx_og / merge(nxgrid_ogrid, max(nxgrid_ogrid-1,1), .false.)
    dy_ogrid = Ly_og / merge(nygrid_ogrid, max(nygrid_ogrid-1,1), .true.)
    dz_ogrid = Lz_og / merge(nzgrid_ogrid, max(nzgrid_ogrid-1,1), lperi(3))
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
        xi1star=find_star(a*xi1lo,a*xi1up,x00,x00+Lx_og,xyz_star_ogrid(1),grid_func_ogrid(1))/a
        call grid_profile(a*(xi1  -xi1star),grid_func_ogrid(1),g1,g1der1,g1der2)
        call grid_profile(a*(xi1lo-xi1star),grid_func_ogrid(1),g1lo)
        call grid_profile(a*(xi1up-xi1star),grid_func_ogrid(1),g1up)
!
        x_ogrid     =x00+Lx_og*(g1  -  g1lo)/(g1up-g1lo)
        xprim_ogrid =    Lx_og*(g1der1*a   )/(g1up-g1lo)
        xprim2_ogrid=    Lx_og*(g1der2*a**2)/(g1up-g1lo)
!
        ! Since lsolid_cells=True
        call grid_profile(a*(xi1proc-xi1star),grid_func_ogrid(1),g1proc)
        g1proc=x00+Lx_og*(g1proc  -  g1lo)/(g1up-g1lo)
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
        dx2_bound_ogrid(nghost:1:-1)  = 2.*(x_ogrid(l2_ogrid)-x_ogrid(l2_ogrid-nghost:l2_ogrid-1))
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
        a=coeff_grid_o(2)*dy_ogrid
        xi2star=find_star(a*xi2lo,a*xi2up,y00,y00+Ly_og,xyz_star_ogrid(2),grid_func_ogrid(2))/a
        call grid_profile(a*(xi2  -xi2star),grid_func_ogrid(2),g2,g2der1,g2der2)
        call grid_profile(a*(xi2lo-xi2star),grid_func_ogrid(2),g2lo)
        call grid_profile(a*(xi2up-xi2star),grid_func_ogrid(2),g2up)
!
        y_ogrid     =y00+Ly_og*(g2  -  g2lo)/(g2up-g2lo)
        yprim_ogrid =    Ly_og*(g2der1*a   )/(g2up-g2lo)
        yprim2_ogrid=    Ly_og*(g2der2*a**2)/(g2up-g2lo)
!
        ! Since lsolid_cells=True
          call grid_profile(a*(xi2proc-xi2star),grid_func_ogrid(2),g2proc)
          g2proc=y00+Ly_og*(g2proc  -  g2lo)/(g2up-g2lo)
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
        a=coeff_grid_o(3)*dz_ogrid
        xi3star=find_star(a*xi3lo,a*xi3up,z00,z00+Lz_og,xyz_star_ogrid(3),grid_func_ogrid(3))/a
        call grid_profile(a*(xi3  -xi3star),grid_func_ogrid(3),g3,g3der1,g3der2)
        call grid_profile(a*(xi3lo-xi3star),grid_func_ogrid(3),g3lo)
        call grid_profile(a*(xi3up-xi3star),grid_func_ogrid(3),g3up)
!
        z_ogrid     =z00+Lz_og*(g3  -  g3lo)/(g3up-g3lo)
        zprim_ogrid =    Lz_og*(g3der1*a   )/(g3up-g3lo)
        zprim2_ogrid=    Lz_og*(g3der2*a**2)/(g3up-g3lo)
!
        ! Since lsolid_cells is True
          call grid_profile(a*(xi3proc-xi3star),grid_func_ogrid(3),g3proc)
          g3proc=z00+Lz_og*(g3proc-g3lo)/(g3up-g3lo)
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
    rcyl_mn_ogrid=x_ogrid(l1_ogrid:l2_ogrid)
    if (x_ogrid(l1_ogrid)==0.) then
      rcyl_mn1_ogrid(2:)=1./x_ogrid(l1_ogrid+1:l2_ogrid)
      rcyl_mn1_ogrid(1)=0.
    else
      rcyl_mn1_ogrid=1./x_ogrid(l1_ogrid:l2_ogrid)
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
      call fatal_error ("initialize_grid", "check Lx_og,Ly_og,Lz_og: is one of them 0?", .true.)
!
    dxmax_ogrid = maxval( (/dxmax_x, dxmax_y, dxmax_z, epsilon(dx_ogrid)/), &
              MASK=((/nxgrid_ogrid, nygrid_ogrid, nzgrid_ogrid, 2/) > 1) )

    call mpiallreduce_max(dxmax,dxmax_x)
    dxmax_ogrid=dxmax_x
!
! Box volume, cylinder symmetrical
!
!TODO: Is this calculation correct? Is it needed?
    box_volume_ogrid=1.
    if (nxgrid_ogrid/=1) then
        box_volume_ogrid = box_volume_ogrid*.5*(r_ogrid**2-xyz0_ogrid(1)**2)
    endif
    box_volume_ogrid = box_volume_ogrid*2.*pi
    if (nzgrid_ogrid/=1) box_volume_ogrid = box_volume_ogrid*Lxyz_ogrid(3)
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
      Area_xz_ogrid=Area_xz_ogrid*Lxyz_ogrid(1)
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
      Area_xz_ogrid=Area_xz_ogrid*Lxyz_ogrid(3)
      Area_yz_ogrid=Area_yz_ogrid*Lxyz_ogrid(3)
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
  endsubroutine initialize_grid_ogrid
!!***********************************************************************
  subroutine calc_pencils_grid_ogrid
!
!  Calculate Grid/geometry related pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   31-jan-17/Jorgen: Adapted from calc_pencils_grid in grid.f90
!                     Only cylindrical coodrinats included
!
    if (lpencil_ogrid(i_og_x_mn))     p_ogrid%x_mn    = x_ogrid(l1_ogrid:l2_ogrid)*cos(y_ogrid(m_ogrid))
    if (lpencil_ogrid(i_og_y_mn))     p_ogrid%y_mn    = x_ogrid(l1_ogrid:l2_ogrid)*sin(y_ogrid(m_ogrid))
    if (lpencil_ogrid(i_og_z_mn))     p_ogrid%z_mn    = spread(z_ogrid(n_ogrid),1,nx_ogrid)
    if (lpencil_ogrid(i_og_rcyl_mn))  p_ogrid%rcyl_mn = x_ogrid(l1_ogrid:l2_ogrid)
    if (lpencil_ogrid(i_og_phi_mn))   p_ogrid%phi_mn  = spread(y_ogrid(m_ogrid),1,nx_ogrid)
    if (lpencil_ogrid(i_og_rcyl_mn1)) p_ogrid%rcyl_mn1= 1./max(p_ogrid%rcyl_mn,tini)
!
  endsubroutine calc_pencils_grid_ogrid
!***********************************************************************
  subroutine real_to_index_ogrid(n, x, xi)
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
          call inverse_grid_ogrid(i, x(:,i), xi(:,i), local=.true.)
        else
          xi(:,i) = ngp1
        endif
      enddo dir
    endif nonzero
!
  endsubroutine real_to_index_ogrid
!***********************************************************************
  subroutine inverse_grid_ogrid(dir, x, xi, local)
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
        call fatal_error('inverse_grid_ogrid', 'lshift_origin and lshift_origin_lower are not supported. ')
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
      call fatal_error('inverse_grid_ogrid', trim(msg))
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
      call fatal_error('inverse_grid_ogrid', 'unknown grid function ' // trim(grid_func_ogrid(dir)))
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
  endsubroutine inverse_grid_ogrid
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
    use Mpicomm, only: mpisend_real,mpirecv_real,mpibcast_real, mpiallreduce_sum_int
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
    call mpibcast_real(xgrid_ogrid,nxgrid_ogrid)
    call mpibcast_real(dx1grid_ogrid,nxgrid_ogrid)
    call mpibcast_real(dxtgrid_ogrid,nxgrid_ogrid)
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
  subroutine time_step_ogrid(f_cartesian)
!
!  Perform time steps on the curvilinear grid, including interpolation of 
!  flow variables back and forth between the overlapping grids.
!  The time iterations should equal to one time step on the cartesian grid
!
!  07-feb-17/Jorgen+Nils: Adapded from timestep.f90
!
    use Mpicomm, only: mpifinalize, mpiallreduce_max
!
    real, dimension (mx,my,mz,mfarray) :: f_cartesian
    real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df_ogrid
    real :: ds, dtsub_ogrid, dt_ogrid, dt_cartesian
    integer :: timestep_factor, tstep_ogrid
    integer :: j
    integer, save :: iterator=0
    real, dimension(3) :: alpha_ts_ogrid=0.,beta_ts_ogrid=0.,dt_beta_ts_ogrid=0.
!
!  Coefficients for up to order 3.
!
    if (itorder==1) then
      alpha_ts_ogrid=(/ 0.0, 0.0, 0.0 /)
      beta_ts_ogrid=(/ 1.0, 0.0, 0.0 /)
    elseif (itorder==2) then
      alpha_ts_ogrid=(/   0.0, -1/2.0, 0.0 /)
      beta_ts_ogrid=(/ 1/2.0,    1.0, 0.0 /)
    elseif (itorder==3) then
      !  use coefficients of Williamson (1980)
      alpha_ts_ogrid=(/   0.0, -5/9.0 , -153/128.0 /)
      beta_ts_ogrid=(/ 1/3.0, 15/16.0,    8/15.0  /)
    else
      if (lroot) print*,'Not implemented: itorder=',itorder
      call mpifinalize
    endif
!
!  Interpolate data from cartesian to curvilinear grid.
!  Before interpolating, necessary points outside this processors domain are
!  recieved from appropriate processor
!
    call communicate_ip_cart_to_curv(f_cartesian,iux,irho)
!
!  Time step for ogrid
!
    timestep_factor = ceiling(dxmin/dxmin_ogrid)   ! number of timesteps on cylinder grid for one timestep on cartesian grid
    if(timestep_factor < 1)  then
      timestep_factor = 1
    endif
    !timestep_factor=1
!
!  Uncomment for forced timestep factor = 1
!  Should be used with care, typically when dt is set in run.in
!
    !timestep_factor=1
!
    dt_ogrid = dt/timestep_factor
    dt_beta_ts_ogrid=dt_ogrid*beta_ts_ogrid
!
!  Perform a number of timestep equal to timestep_factor, such that the
!  endtime t_ogrid is equal to t after the timesteps
!
    do tstep_ogrid=1,timestep_factor
!
!  Set up df and ds for each time sub.
!
      do itsub=1,itorder
        llast_ogrid=(tstep_ogrid==timestep_factor).and.(itsub==itorder)
        if (itsub==1) then
          df_ogrid=0.0
          ds=0.0
        else
          df_ogrid=alpha_ts_ogrid(itsub)*df_ogrid !(could be subsumed into pde, but is dangerous!)
          ds=alpha_ts_ogrid(itsub)*ds
        endif
!
!  Change df according to the chosen physics modules.
!
        call pde_ogrid(df_ogrid)
!
        ds=ds+1.0
!
!  Calculate substep 
!
        dtsub_ogrid = ds * dt_beta_ts_ogrid(itsub)
!
!  Time evolution of grid variables.
!  (do this loop in pencils, for cache efficiency)
!
        do j=1,mvar 
          do n_ogrid=n1_ogrid,n2_ogrid
            do m_ogrid=m1_ogrid,m2_ogrid
              f_ogrid(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,j)=f_ogrid(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,j) &
                +dt_beta_ts_ogrid(itsub)*df_ogrid(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,j)
            enddo
          enddo
        enddo
      enddo
    enddo
!
!  Interpolate data from curvilinear to cartesian grid.
!  Before interpolating, necessary points outside this processors domain are
!  recieved from appropriate processor
!
    call update_ghosts_ogrid
    call communicate_ip_curv_to_cart(f_cartesian,iux,irho)
!
! !!TODO:Printing
! if(ncpus==1 .or. iproc==1) then
!     iterator = iterator+1
!   if(mod(iterator,25)==0) then
!     call print_ogrid(int(iterator/25))
!     call print_cgrid(int(iterator/25),f_cartesian)
!   endif
! endif

    call wsnap_ogrid('OGVAR',ENUM=.true.,FLIST='ogvarN.list')
!
  endsubroutine time_step_ogrid
!***********************************************************************
    subroutine print_ogrid(num)
!  Print to file
    character(len=16) :: xfile,yfile,fxfile,fyfile,rhfile
    integer :: i,num

    
    if(num<10) then
      xfile='x_ogrid0'
      yfile='y_ogrid0'
      fxfile='fx_ogrid0'
      fyfile='fy_ogrid0'
      rhfile='rh_ogrid0'
      write(fxfile,"(A9,I1)") trim(fxfile),num
      write(fyfile,"(A9,I1)") trim(fyfile),num
      write(rhfile,"(A9,I1)") trim(rhfile),num
    else
      xfile='x_ogrid'
      yfile='y_ogrid'
      fxfile='fx_ogrid'
      fyfile='fy_ogrid'
      rhfile='rh_ogrid'
      write(fxfile,"(A8,I2)") trim(fxfile),num
      write(fyfile,"(A8,I2)") trim(fyfile),num
      write(rhfile,"(A8,I2)") trim(rhfile),num
    endif

    xfile=trim(xfile)//'.dat'
    yfile=trim(yfile)//'.dat'
    fxfile=trim(fxfile)//'.dat'
    fyfile=trim(fyfile)//'.dat'
    rhfile=trim(rhfile)//'.dat'
    open(unit=1,file=trim(xfile))
    open(unit=2,file=trim(yfile))
    open(unit=11,file=trim(fxfile))
    open(unit=12,file=trim(fyfile))
    open(unit=13,file=trim(rhfile))
    !TODO: nghost
    do i=l1_ogrid-nghost,l2_ogrid+nghost
      write(1,*) x_ogrid(i)*cos(y_ogrid(m1_ogrid:m2_ogrid))
      write(2,*) x_ogrid(i)*sin(y_ogrid(m1_ogrid:m2_ogrid))
      write(11,*) f_ogrid(i,m1_ogrid:m2_ogrid,4,iux)*cos(y_ogrid(m1_ogrid:m2_ogrid)) &
            -f_ogrid(i,m1_ogrid:m2_ogrid,4,iuy)*sin(y_ogrid(m1_ogrid:m2_ogrid))
      write(12,*) f_ogrid(i,m1_ogrid:m2_ogrid,4,iux)*sin(y_ogrid(m1_ogrid:m2_ogrid)) &
            +f_ogrid(i,m1_ogrid:m2_ogrid,4,iuy)*cos(y_ogrid(m1_ogrid:m2_ogrid))
      write(13,*) f_ogrid(i,m1_ogrid:m2_ogrid,4,irho)
    enddo

    close(1)
    close(2)
    close(11)
    close(12)
    close(13)
    !write(*,*) 'Press any key to continue'
    !read(*,*) 

    endsubroutine print_ogrid
!***********************************************************************
    subroutine print_cgrid(num,f_cartesian)
!  Print to file
    character(len=16) :: xfile,yfile,fxfile,fyfile,rhfile
    real, dimension (mx,my,mz,mfarray) :: f_cartesian
    integer :: i,num

    if(num<10) then
      xfile='x_cgrid0'
      yfile='y_cgrid0'
      fxfile='fx_cgrid0'
      fyfile='fy_cgrid0'
      rhfile='rh_cgrid0'
      write(fxfile,"(A9,I1)") trim(fxfile),num
      write(fyfile,"(A9,I1)") trim(fyfile),num
      write(rhfile,"(A9,I1)") trim(rhfile),num
    else
      xfile='x_cgrid'
      yfile='y_cgrid'
      fxfile='fx_cgrid'
      fyfile='fy_cgrid'
      rhfile='rh_cgrid'
      write(fxfile,"(A8,I2)") trim(fxfile),num
      write(fyfile,"(A8,I2)") trim(fyfile),num
      write(rhfile,"(A8,I2)") trim(rhfile),num
    endif

    xfile=trim(xfile)//'.dat'
    yfile=trim(yfile)//'.dat'
    fxfile=trim(fxfile)//'.dat'
    fyfile=trim(fyfile)//'.dat'
    rhfile=trim(rhfile)//'.dat'
    open(unit=1,file=trim(xfile))
    open(unit=2,file=trim(yfile))
    open(unit=11,file=trim(fxfile))
    open(unit=12,file=trim(fyfile))
    open(unit=13,file=trim(rhfile))
    do i=l1,l2
      write(1,*) x(i)
      write(11,*) f_cartesian(i,m1:m2,4,iux)
      write(12,*) f_cartesian(i,m1:m2,4,iuy)
      write(13,*) f_cartesian(i,m1:m2,4,irho)
    enddo
    write(2,*) y(m1:m2)

    close(1)
    close(2)
    close(11)
    close(12)
    close(13)
    !write(*,*) 'Press any key to continue'
    !read(*,*) 

    endsubroutine print_cgrid
!***********************************************************************
    subroutine pde_ogrid(df)
!
!  06-feb-17/Jorgen+Nils: Adapded from equ.f90
!
!  Call the different evolution equations.
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df
      integer :: nyz_ogrid
!
      intent(out)    :: df
      real :: c_dragx,c_dragy
      c_dragx=0.
      c_dragy=0.
!
!  Initiate communication and do boundary conditions.
!
      call update_ghosts_ogrid
!
!------------------------------------------------------------------------------
!  Do loop over m and n.
!
      nyz_ogrid=ny_ogrid*nz_ogrid
      mn_loop: do imn_ogrid=1,nyz_ogrid
        n_ogrid=nn_ogrid(imn_ogrid)
        m_ogrid=mm_ogrid(imn_ogrid)
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
        call calc_pencils_energy_ogrid
        call calc_pencils_viscosity_ogrid
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
!  Compute drag and lift coefficient, if this is the last sub-timestep
!

        if(llast_ogrid.and.lfirst_proc_x.and.((idiag_c_dragx/=0).or.(idiag_c_dragy/=0))) then
          call drag_force_pencils(c_dragx,c_dragy)
        endif
!
!  End of loops over m and n.
!
      enddo mn_loop
!
      if(llast_ogrid.and.((idiag_c_dragx/=0).or.(idiag_c_dragy/=0))) call drag_coeffs(c_dragx,c_dragy)
!
!  -------------------------------------------------------------
!  NO CALLS MODIFYING DF BEYOND THIS POINT (APART FROM FREEZING)
!  -------------------------------------------------------------
!  Frerzing must be done after the full (m,n) loop, as df may be modified
!  outside of the considered pencil.
!
!  Freezing boundary conditions in x (radial direction)
!
      do imn_ogrid=1,nyz_ogrid
        n_ogrid=nn_ogrid(imn_ogrid)
        m_ogrid=mm_ogrid(imn_ogrid)
        df(l1_ogrid,m_ogrid,n_ogrid,iux:irho) = 0.
      enddo
!
    endsubroutine pde_ogrid
!***********************************************************************
    subroutine duu_dt_ogrid(df)
!
!  velocity evolution
!  calculate du/dt = - u.gradu - 2Omega x u + grav + Fvisc
!  pressure gradient force added in density and entropy modules.
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df
      intent(inout) :: df
!
!  Advection term.
!
      df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iux:iuz) = df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iux:iuz)-p_ogrid%ugu
!
!  Viscous term
!
      df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iux:iuz) = df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iux:iuz) + p_ogrid%fvisc
!
    endsubroutine duu_dt_ogrid
!***********************************************************************
    subroutine dlnrho_dt_ogrid(df)
!
!  Continuity equation.
!  Calculate dlnrho/dt = - u.gradlnrho - divu
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df
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
      df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,irho) = df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,irho) + density_rhs
!
    endsubroutine dlnrho_dt_ogrid
!***********************************************************************
    subroutine denergy_dt_ogrid(df)
!
!  Calculate pressure gradient term for isothermal/polytropic equation
!  of state.
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df
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
    subroutine calc_pencils_hydro_ogrid()
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
! Pencils: uu, u2, uij, divu, sij, sij2, ugu, ugu2, del2u
!
      if (lpencil_ogrid(i_og_uu)) p_ogrid%uu=f_ogrid(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iux:iuz)
      if (lpencil_ogrid(i_og_u2)) call dot2_mn_ogrid(p_ogrid%uu,p_ogrid%u2)
      if (lpencil_ogrid(i_og_uij)) call gij_ogrid(f_ogrid,iuu,p_ogrid%uij,1)
      if (lpencil_ogrid(i_og_divu)) call div_mn_ogrid(p_ogrid%uij,p_ogrid%divu,p_ogrid%uu)
      if (lpencil_ogrid(i_og_sij)) call traceless_strain_ogrid(p_ogrid%uij,p_ogrid%divu,p_ogrid%sij,p_ogrid%uu)
      if (lpencil_ogrid(i_og_sij2)) call multm2_sym_mn_ogrid(p_ogrid%sij,p_ogrid%sij2)
      if (lpencil_ogrid(i_og_ugu)) call u_dot_grad_ogrid(f_ogrid,iuu,p_ogrid%uij,p_ogrid%uu,p_ogrid%ugu)
      if (lpencil_ogrid(i_og_ugu2)) call dot2_mn_ogrid(p_ogrid%ugu,p_ogrid%ugu2)
      if (lpencil_ogrid(i_og_graddivu).and.lpencil_ogrid(i_og_del2u)) then
        call gij_etc_ogrid(f_ogrid,iuu,p_ogrid%uu,p_ogrid%uij,DEL2=p_ogrid%del2u,GRADDIV=p_ogrid%graddivu)
      elseif(lpencil_ogrid(i_og_graddivu)) then
        call gij_etc_ogrid(f_ogrid,iuu,p_ogrid%uu,p_ogrid%uij,GRADDIV=p_ogrid%graddivu)
      elseif(lpencil_ogrid(i_og_del2u)) then
        call gij_etc_ogrid(f_ogrid,iuu,p_ogrid%uu,p_ogrid%uij,DEL2=p_ogrid%del2u)
      endif
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
      use density, only:lupw_lnrho
!
      integer :: i
      logical :: lupw_rho=.false.
!
      if (.not. ldensity_nolog) then
        call fatal_error('calc_pencils_density_ogrid','Must use linear density for solid_cells_ogrid')
      endif
      if(lupw_lnrho) then
        lupw_rho=.true.
      endif
!
! Pencils: rho, rho1, lnrho, glnrho, grho, ugrho, sglnrho
!
      p_ogrid%rho=f_ogrid(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,irho)
      if (lpencil_ogrid(i_og_rho1)) p_ogrid%rho1=1.0/p_ogrid%rho
      if (lpencil_ogrid(i_og_lnrho)) p_ogrid%lnrho=log(p_ogrid%rho)
      if (lpencil_ogrid(i_og_glnrho)) then
        call grad_ogrid(f_ogrid,irho,p_ogrid%grho)
        do i=1,3
          p_ogrid%glnrho(:,i)=p_ogrid%rho1*p_ogrid%grho(:,i)
        enddo
      endif
      if (lpencil_ogrid(i_og_ugrho)) call u_dot_grad_ogrid(f_ogrid,ilnrho,p_ogrid%grho,p_ogrid%uu,p_ogrid%ugrho,UPWIND=lupw_rho)
      if (lpencil_ogrid(i_og_sglnrho)) call multmv_mn_ogrid(p_ogrid%sij,p_ogrid%glnrho,p_ogrid%sglnrho)
!
! rho
! rho1
! lnrho
! glnrho and grho
! uglnrho
! ugrho
! sglnrho
! uij5glnrho
      !if (lpencil(i_uij5glnrho)) call multmv(p%uij5,p%glnrho,p%uij5glnrho)
!
      !if (lpencil(i_uuadvec_grho)) call h_dot_grad(p%uu_advec,p%grho,p%uuadvec_grho)
!
    endsubroutine calc_pencils_density_ogrid
!***********************************************************************
    subroutine calc_pencils_eos_ogrid()
!
!  Calculate EquationOfState pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  10-feb-17/Jorgen+Nils: Adapted from calc_pencils_eos_pencpar in eos_idealgas.f90
!
      use EquationOfState, only: get_cv1,get_cp1,cs20,gamma_m1,gamma1
!
      real, dimension(nx_ogrid) :: tmp
      integer :: i,j
      real :: cp1, cv1, cp, cv
      logical :: leos_isentropic=.true.
      logical :: leos_isothermal=.false.
!
!  Inverse cv and cp values.
!
      call get_cp1(cp1)
      call get_cv1(cv1)
      cp=1./cp1
      cv=1./cv1
!!  !
!!  !  Work out thermodynamic quantities for given lnrho or rho and TT.
!!  !
!!        if (iTT .gt. 0) then
!!          if (lpencil_ogrid(i_og_TT))   p_ogrid%TT=f_ogrid(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT)
!!          if (lpencil_ogrid(i_og_TT1))  p_ogrid%TT1=1/f_ogrid(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT)
!!          if (lpencil_ogrid(i_og_cs2))  p_ogrid%cs2=cp*gamma_m1*f_ogrid(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT)
!!          if (lpencil_ogrid(i_og_gTT))  call grad_ogrid(f_ogrid,iTT,p_ogrid%gTT)
!!          if (lpencil_ogrid(i_og_del2TT)) &
!!              call del2_ogrid(f_ogrid,iTT,p_ogrid%del2TT)
!!          if (lpencil_ogrid(i_og_pp)) p_ogrid%pp=cv*gamma_m1*p_ogrid%rho*p_ogrid%TT
!!          if (lpencil_ogrid(i_og_ee)) p_ogrid%ee=cv*p_ogrid%TT
!!  !
!!  !  Work out thermodynamic quantities for given lnrho or rho and cs2.
!!  !
!!        else
!!          if (leos_isentropic) then
!!            call fatal_error('calc_pencils_eos', &
!!                'leos_isentropic not implemented for ilnrho_cs2, try ilnrho_ss')
!!          elseif (leos_isothermal) then
!!            if (lpencil_ogrid(i_og_cs2)) p_ogrid%cs2=cs20
!!            if (lpencil_ogrid(i_og_pp)) p_ogrid%pp=gamma1*p_ogrid%rho*cs20
!!          else
!!            call fatal_error('calc_pencils_eos', &
!!                'Full equation of state not implemented for ilnrho_cs2')
!!          endif
!!        endif
!
!  Work out thermodynamic quantities for given lnrho or rho and ss.
!
      if (leos_isentropic) then
        if (lpencil_ogrid(i_og_ss)) p_ogrid%ss=0.0
        if (lpencil_ogrid(i_og_cs2)) p_ogrid%cs2=cs20*exp(gamma_m1*(p_ogrid%lnrho-lnrho0))
      elseif (leos_isothermal) then
        if (lpencil_ogrid(i_og_ss)) p_ogrid%ss=-(cp-cv)*(p_ogrid%lnrho-lnrho0)
        if (lpencil_ogrid(i_og_cs2)) p_ogrid%cs2=cs20
      endif
!
      if (lpencil_ogrid(i_og_lnTT)) p_ogrid%lnTT=lnTT0+cv1*p_ogrid%ss+gamma_m1*(p_ogrid%lnrho-lnrho0)
      if (lpencil_ogrid(i_og_pp)) p_ogrid%pp=(cp-cv)*exp(p_ogrid%lnTT+p_ogrid%lnrho)
    endsubroutine calc_pencils_eos_ogrid
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
    subroutine calc_pencils_viscosity_ogrid()
!
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  10-feb-17/Jorgen+Nils: Adapted from rountine in viscosity.f90
!
      use viscosity, only:getnu
      real :: nu
      integer :: j
!
      call getnu(nu_input=nu)
!      
!  Viscous force and viscous heating are calculated here (for later use).
!
      p_ogrid%fvisc=0.0                              
!
!  viscous force: nu*del2v
!  -- not physically correct (no momentum conservation), but
!  numerically easy and in most cases qualitatively OK
!
      !p_ogrid%fvisc=p_ogrid%fvisc+nu*p_ogrid%del2u
!
      do j=1,3
        p_ogrid%fvisc(:,j) = p_ogrid%fvisc(:,j) + nu*(2*p_ogrid%sglnrho(:,j)+p_ogrid%del2u(:,j) + 1./3.*p_ogrid%graddivu(:,j))
      enddo
    end subroutine calc_pencils_viscosity_ogrid
!***********************************************************************
!  FROM BOUNDCOND.F90
!***********************************************************************
!  ROUTINES
!    boundconds_x_ogrid
!    boundconds_y_ogrid
!    boundconds_z_ogrid
!***********************************************************************
    subroutine boundconds_x_ogrid
!
!  Boundary conditions at the cylinder surface in the radial direction (xdir).
!  For ogrids, only boundary conditions at cylinder surface is set. The BC on
!  the 'top' is set by interpolation from cartesian grid, outside the timestep.
!  Only need to conpute boundary value for the density, using stencil that
!  satisfies the SBP energy conservation. No-slip on the surface is respected
!  automatically since we set df(l1_ogrid,:,:,:)=0 after the mn-loop (freeze).
!  If SBP is not used, the grid ponts inside the surface are computed using 
!  one-sided differences. Note that these do not guarantee stability!
!  
!  Remark: boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-'corners').
!
!  06-feb-17/Jorgen: Adapted from boundcond.f90 to be used for ogrids
!  06-apr-17/Jorgen: Cleanup and working with SBP property
!
!  Only set cyinder boundary here, not processor boundaries
!
      if(lfirst_proc_x) then
        if(SBP) then
          call bval_from_neumann_arr_ogrid_alt
        else
          call set_ghosts_onesided_ogrid(iux)
          call set_ghosts_onesided_ogrid(iuy)
          call set_ghosts_onesided_ogrid(iuz)
          call bval_from_neumann_arr_ogrid
          call set_ghosts_onesided_ogrid(irho)
        endif
      endif
!
    endsubroutine boundconds_x_ogrid
!***********************************************************************
    subroutine boundconds_y_ogrid
!
!  Periodic boundary condition for runs with a single processor in y-direction 
!
!  06-feb-17/Jorgen: Adapted from boundcond.f90 to be used for ogrids where the
!                    y-dir is always periodic
!
      integer :: ivar1, ivar2, j
      integer :: m1i_ogrid=m1_ogrid+nghost-1
      integer :: m2i_ogrid=my_ogrid-2*nghost+1
!
      ivar1=1; ivar2=min(mcom,size(f_ogrid,4))
!
!  Boundary conditions in y
!  Periodic, with y being the theta direction for the cylindrical grid
!
      do j=ivar1,ivar2
!  Bottom boundary
        f_ogrid(:,1:m1_ogrid-1,:,j) = f_ogrid(:,m2i_ogrid:m2_ogrid,:,j)
!  Top boundary
        f_ogrid(:,m2_ogrid+1:,:,j) = f_ogrid(:,m1_ogrid:m1i_ogrid,:,j)
      enddo
!
    endsubroutine boundconds_y_ogrid
!***********************************************************************
    subroutine boundconds_z_ogrid
!
!  Periodic boundary condition for 3D-runs with a single processor in z-direction 
!
!  06-feb-17/Jorgen: Adapted from boundcond.f90 to be used for ogrids where the
!                    z-dir is always periodic as long as nzgrid=/1
!
      integer :: ivar1, ivar2, j
      integer :: n1i_ogrid=n1_ogrid+nghost-1
      integer :: n2i_ogrid=mz_ogrid-2*nghost+1
!
      ivar1=1; ivar2=min(mcom,size(f_ogrid,4))
!
!  Boundary conditions in z
!
      do j=ivar1,ivar2
!  Bottom boundary
        f_ogrid(:,:,1:n1_ogrid-1,j) = f_ogrid(:,:,n2i_ogrid:n2_ogrid,j)
!  Top boundary
        f_ogrid(:,:,n2_ogrid+1:,j) = f_ogrid(:,:,n1_ogrid:n1i_ogrid,j)
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
!   graddiv_ogrid         ! Obsolete
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
    subroutine g2ij_ogrid(f,k,g)
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
          call derij_ogrid(f,k,tmp,i,j)
          g(:,i,j)=tmp
          g(:,j,i)=tmp
        enddo
      enddo
      if(SBP) then
        call fatal_error('solid_cells: g2ij_ogrid',&
          'not implemented with summation by parts property')
      endif
!
    endsubroutine g2ij_ogrid
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
    intent(in)  :: uij, divu
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
      real, dimension (nx_ogrid,3,3), intent (in) :: aij
      real, dimension (nx_ogrid,3),   intent (in), optional :: a
      real, dimension (nx_ogrid,3),   intent (out) :: b
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
!  14-feb-17/Jorgen: Adapted from del2v_ogrid in der.f90
!  15-mar-07/wlad: added cylindrical coordinates
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
        call del2_ogrid(f,k1+i,tmp)
        del2f(:,i)=tmp
      enddo
!
      !del2 already contains the extra term 1/r*d(uk)/dt
      call der_ogrid(f,k1+2,tmp,2)
      del2f(:,1)=del2f(:,1) -(2*tmp+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,k1+1))*rcyl_mn2_ogrid
      call der_ogrid(f,k1+1,tmp,2)
      del2f(:,2)=del2f(:,2) +(2*tmp-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,k1+2))*rcyl_mn2_ogrid
!
    endsubroutine del2v_ogrid
!***********************************************************************
    subroutine u_dot_grad_vec_ogrid(f,k,gradf,uu,ugradf,upwind)
!
!  u.gradu
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      use General, only: loptest
!
      intent(in) :: f,k,gradf,uu,upwind
      intent(out) :: ugradf
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid,3,3) :: gradf
      real, dimension (nx_ogrid,3) :: uu,ff,ugradf,grad_f_tmp
      real, dimension (nx_ogrid) :: tmp
      integer :: j,k
      logical, optional :: upwind
!
      if (k<1 .or. k>mfarray) then
        call fatal_error('u_dot_grad_vec_ogrid','variable index is out of bounds')
        return
      endif
!
      do j=1,3
        grad_f_tmp = gradf(:,j,:)
        call u_dot_grad_scl_ogrid(f,k+j-1,grad_f_tmp,uu,tmp,UPWIND=loptest(upwind))
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
    subroutine u_dot_grad_scl_ogrid(f,k,gradf,uu,ugradf,upwind)
!
!  Do advection-type term u.grad f_k.
!  Assumes gradf to be known, but takes f and k as arguments to be able
!  to calculate upwind correction.
!
!  07-feb-17/Jorgen: Adapted from sub.f90
!
      use General, only: loptest
!
      intent(in) :: f,k,gradf,uu,upwind
      intent(out) :: ugradf
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid,3) :: uu,gradf
      real, dimension (nx_ogrid) :: ugradf
      integer :: k
      logical, optional :: upwind
!
      if (k<1 .or. k>mfarray) then
        call fatal_error('u_dot_grad_scl_ogrid','variable index is out of bounds')
        return
      endif
!
      call dot_mn_ogrid(uu,gradf,ugradf)
!
!  Upwind correction
!
      if (present(upwind)) then
        if (upwind) call doupwind_ogrid(f,k,uu,ugradf)
      endif
!
    endsubroutine u_dot_grad_scl_ogrid
!***********************************************************************
    subroutine doupwind_ogrid(f,k,uu,ugradf)
!
!  Calculates upwind correction, works incrementally on ugradf
!
!  27-feb-17/Jorgen: Adapted from sub.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray), intent(IN)    :: f
      integer :: k
      real, dimension (nx_ogrid,3),                         intent(IN)    :: uu
      real, dimension (nx_ogrid),                           intent(INOUT) :: ugradf
!
      real, dimension (nx_ogrid,3) :: del6f
      integer                      :: ii,msk
      integer, dimension(nx_ogrid) :: indxs
!
      do ii=1,3
!
        if ( lequidist_ogrid(ii) ) then
          call der6_ogrid(f,k,del6f(1,ii),ii,UPWIND=.true.)
        else
          where( uu(:,ii)>=0 )
            indxs = 7
          elsewhere
            indxs = 8
          endwhere
          call deri_3d_inds_ogrid(f(1,1,1,k),del6f(1,ii),indxs,ii,lnometric=.true.)
        endif
!
        del6f(:,ii) = abs(uu(:,ii))*del6f(:,ii)
!
      enddo
!
      del6f(:,2) = rcyl_mn1_ogrid*del6f(:,2)
!
!
      ugradf = ugradf-sum(del6f,2)
!
    endsubroutine doupwind_ogrid
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
!
      intent(in) :: a,b
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
!!     subroutine graddiv_ogrid(f,k,graddiv)
!! !
!! !  Calculates a number of second derivative expressions of a vector
!! !  outputs a number of different vector fields.
!! !  gradcurl is not the vector gradient.
!! !  Surprisingly, calling derij only if graddiv or curlcurl are present
!! !  does not speed up the code on Mephisto @ 32x32x64.
!! !
!! !  05-apr-17/Jorgen - Adapted from del2v_etc in sub.f90
!! !
!!       real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
!!       real, dimension (nx_ogrid,3,3) :: fjji,fijj
!!       real, dimension (nx_ogrid,3) :: graddiv
!!       real, dimension (nx_ogrid) :: tmp
!!       integer :: i,j,k,k1
!! !
!!       intent(in) :: f,k
!!       intent(out) :: del2,graddiv,curlcurl,gradcurl
!! !
!! !  calculate f_{i,jj} and f_{j,ji}
!! !  AJ: graddiv needs diagonal elements from the first tmp (derij only sets
!! !      off-diagonal elements)
!! !
!!       k1=k-1
!!       do i=1,3
!!       do j=1,3
!!         call der2_ogrid(f,k1+i,tmp,  j) 
!!         fijj(:,i,j)=tmp  ! f_{i,jj}
!!         call derij_ogrid(f,k1+j,tmp,j,i) 
!!         fjji(:,i,j)=tmp  ! f_{j,ji}
!!         endif
!!       enddo
!!       enddo
!! !
!! !  the diagonal terms have not been set in derij; do this now
!! !  ** They are automatically set above, because derij   **
!! !  ** doesn't overwrite the value of tmp for i=j!       **
!! !
!! !     do j=1,3
!! !       fjji(:,j,j)=fijj(:,j,j)
!! !     enddo
!! !
!! !  calculate f_{i,jk} for i /= j /= k
!! !
!! !
!!       do i=1,3
!!         graddiv(:,i)=fjji(:,i,1)+fjji(:,i,2)+fjji(:,i,3)
!!       enddo
!!       ! Since we have cylindrical coordinates
!!       call der(f,k1+1,tmp,1)
!!       graddiv(:,1)=graddiv(:,1)+tmp*rcyl_mn1_ogrid-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,k1+1)*rcyl_mn2_ogrid
!!       call der(f,k1+1,tmp,2)
!!       graddiv(:,2)=graddiv(:,2)+tmp*rcyl_mn1_ogrid
!!       call der(f,k1+1,tmp,3)
!!       graddiv(:,3)=graddiv(:,3)+tmp*rcyl_mn1_ogrid
!! !
!!     endsubroutine graddiv_ogrid
!***********************************************************************
    subroutine gij_etc_ogrid(f,iref,aa,aij,bij,del2,graddiv,lcovariant_derivative)
!
!  Calculate B_i,j = eps_ikl A_l,jk and A_l,kk.
!
!  05-apr-17/Jorgen - Adapted from gij_etc in sub.f90
!
      use Deriv, only: der2,derij
      use General, only: loptest
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray), intent (in) :: f
      integer, intent (in) :: iref
      logical, intent (in), optional :: lcovariant_derivative
      real, dimension (nx_ogrid,3), intent (in)   :: aa
      real, dimension (nx_ogrid,3,3), intent (in) :: aij
      real, dimension (nx_ogrid,3,3), intent (out), optional :: bij
      real, dimension (nx_ogrid,3), intent (out), optional :: del2,graddiv
!
!  Locally used variables.
!
      real, dimension (nx_ogrid,3,3,3) :: d2A
      integer :: iref1,i,j
!
!  Reference point of argument.
!
      iref1=iref-1
!
!  Calculate all mixed and non-mixed second derivatives
!  of the vector potential (A_k,ij).
!
!  Do not calculate both d^2 A_k/(dx dy) and d^2 A_k/(dy dx).
!  (This wasn't spotted by me but by a guy from SGI...)
!  Note: for non-cartesian coordinates there are different correction terms,
!  see below.
!
      do i=1,3
        do j=1,3
          call der2_ogrid(f,iref1+i,d2A(:,j,j,i),j)
        enddo
        call derij_ogrid(f,iref1+i,d2A(:,2,3,i),2,3); d2A(:,3,2,i)=d2A(:,2,3,i)
        call derij_ogrid(f,iref1+i,d2A(:,3,1,i),3,1); d2A(:,1,3,i)=d2A(:,3,1,i)
        call derij_ogrid(f,iref1+i,d2A(:,1,2,i),1,2); d2A(:,2,1,i)=d2A(:,1,2,i)
      enddo
!
!  Since we have cylindrical coordinates
!  Psi_{,phi^ pom^} = Psi_{,pom^ phi^} - Psi_{,\phi^}/pom .
!
        do i=1,3
          d2A(:,2,1,i)=d2A(:,2,1,i)-aij(:,i,2)*rcyl_mn1_ogrid
        enddo
!
!  Calculate optionally b_i,j = eps_ikl A_l,kj,
!  del2_i = A_i,jj and graddiv_i = A_j,ji .
!
      if (present(bij)) then
!
        bij(:,1,:)=d2A(:,2,:,3)-d2A(:,3,:,2)
        bij(:,2,:)=d2A(:,3,:,1)-d2A(:,1,:,3)
        bij(:,3,:)=d2A(:,1,:,2)-d2A(:,2,:,1)
!  Corrections for cylindrical coordinates.
          bij(:,3,2)=bij(:,3,2)+ aij(:,2,2)*rcyl_mn1_ogrid
!          !bij(:,3,1)=bij(:,3,1)+(aij(:,2,1)+aij(:,1,2))*rcyl_mn1-aa(:,2)*rcyl_mn2
!  FAG:Partial correction to -d2A(:,2,1,1) already applied above +aij(:,i,2)*rcyl_mn1
!
          bij(:,3,1)=bij(:,3,1)+aij(:,2,1)*rcyl_mn1_ogrid-aa(:,2)*rcyl_mn2_ogrid
          if (loptest(lcovariant_derivative)) then
            !bij(:,1,1)=bij(:,1,1)
            bij(:,1,2)=bij(:,1,2)+(aij(:,3,1)-aij(:,1,3))*rcyl_mn1_ogrid
            !bij(:,1,3)=bij(:,1,3)
            !bij(:,2,1)=bij(:,2,1)
            bij(:,2,2)=bij(:,2,2)+(aij(:,3,2)-aij(:,2,3))*rcyl_mn1_ogrid
            !bij(:,2,3)=bij(:,2,3)
            !bij(:,3,1)=bij(:,3,1)
            !bij(:,3,2)=bij(:,3,2)
            bij(:,3,3)=bij(:,3,3)+aij(:,2,3)*rcyl_mn1_ogrid
          endif
      endif
!
!  Calculate del2 and graddiv, if requested.
!
! HEREHEREHERE: THIS CANNOT BE CORRECT FOR CYLINDER COORDINATES
!               MUST FIX THIS
      if (present(graddiv)) then
        graddiv(:,:)=d2A(:,1,:,1)+d2A(:,2,:,2)+d2A(:,3,:,3)
!!!! ONLY FOR SHPERICAL COORDINATES?
!!!! !  Since we have cylindrical coordinates
!!!!         graddiv(:,1)=graddiv(:,1)+aij(:,1,1)*rcyl_mn1_ogrid*2+ &
!!!!            aij(:,2,1)*rcyl_mn1_ogrid*cotth(m_ogrid) &
!!!!            -aa(:,2)*rcyl_mn2_ogrid*cotth(m_ogrid)-aa(:,1)*rcyl_mn2_ogrid*2
!!!!         graddiv(:,2)=graddiv(:,2)+aij(:,1,2)*rcyl_mn1_ogrid*2+ &
!!!!            aij(:,2,2)*rcyl_mn1_ogrid*cotth(m_ogrid) &
!!!!            -aa(:,2)*rcyl_mn2_ogrid*sin2th(m_ogrid)
!!!!         graddiv(:,3)=graddiv(:,3)+aij(:,1,3)*rcyl_mn1_ogrid*2+ &
!!!!            aij(:,2,3)*rcyl_mn1_ogrid*cotth(m_ogrid)
      endif
!
      if (present(del2)) then
        del2(:,:)=d2A(:,1,1,:)+d2A(:,2,2,:)+d2A(:,3,3,:)
!  Since we have cylindrical coordinates
        del2(:,1)= del2(:,1)+&
          rcyl_mn1_ogrid*(aij(:,1,1)-2*aij(:,2,2))-rcyl_mn2_ogrid*aa(:,1)
        del2(:,2)=del2(:,2)+&               
          rcyl_mn1_ogrid*(aij(:,2,1)-2*aij(:,1,2))-rcyl_mn2_ogrid*aa(:,2)
      endif
!
    endsubroutine gij_etc_ogrid
!***********************************************************************
!***********************************************************************
!***********************************************************************
! FROM DERIV.F90
!***********************************************************************
! ROUTINES
!   der_ogrid
!   der2_ogrid
!   derij_ogrid
!   der6_ogrid
!   deri_3d_inds
!***********************************************************************
    subroutine der_ogrid(f, k, df, j)
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
!
      integer :: i
      real, parameter :: a = 1.0 / 60.0
      real, dimension(nx_ogrid) :: fac
!
      if (j==1) then
        if (nxgrid_ogrid/=1) then
          if(SBP.and.lfirst_proc_x) then
            call der_ogrid_SBP(f(1:l1_ogrid+8,:,:,:),k,df(1:6))
            i=6
          else
            i=0
          endif
          fac = a * dx_1_ogrid(l1_ogrid:l2_ogrid)
          df(1+i:nx_ogrid)=fac(1+i:nx_ogrid) * &
                          (+ 45.0*(f(l1_ogrid+1+i:l2_ogrid+1,m_ogrid,n_ogrid,k)-f(l1_ogrid-1+i:l2_ogrid-1,m_ogrid,n_ogrid,k)) &
                           -  9.0*(f(l1_ogrid+2+i:l2_ogrid+2,m_ogrid,n_ogrid,k)-f(l1_ogrid-2+i:l2_ogrid-2,m_ogrid,n_ogrid,k)) &
                           +      (f(l1_ogrid+3+i:l2_ogrid+3,m_ogrid,n_ogrid,k)-f(l1_ogrid-3+i:l2_ogrid-3,m_ogrid,n_ogrid,k)))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid_ogrid/=1) then
          fac = a * dy_1_ogrid(m_ogrid)
          df=fac*(+ 45.0*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid,k)-f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid,k)) &
                  -  9.0*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid,k)-f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid,k)) &
                  +      (f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid,k)-f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid,k)))
          ! Since we have cylindrical coordinates
          df = df * rcyl_mn1_ogrid
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid_ogrid/=1) then
          fac = a * dz_1_ogrid(n_ogrid)
          df=fac*(+ 45.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+1,k)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-1,k)) &
                  -  9.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+2,k)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-2,k)) &
                  +      (f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+3,k)-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-3,k)))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in z-direction'
        endif
      endif
      !if(m_ogrid == 16 .and. j==1) then
      !  print*, ''
      !  print*, 'df,k'
      !  print*, df,k
      !endif
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
      integer :: j,k,i
!
      real :: der2_coef0, der2_coef1, der2_coef2, der2_coef3

      intent(in)  :: f,k,j
      intent(out) :: df2
!
      der2_coef0=-490./180.; der2_coef1=270./180.
      der2_coef2=-27./180.; der2_coef3=2./180.

      if (j==1) then
        if (nxgrid_ogrid/=1) then
          if(SBP.and.lfirst_proc_x) then
            call der2_ogrid_SBP(f(1:l1_ogrid+8,:,:,:),k,df2(1:6))
            i=6
          else
            i=0
          endif
          fac=dx_1_ogrid(l1_ogrid:l2_ogrid)**2
          df2(1+i:nx_ogrid)=fac(1+i:nx_ogrid) * &
                  (der2_coef0* f(l1_ogrid  +i:l2_ogrid  ,m_ogrid,n_ogrid,k) &
                  +der2_coef1*(f(l1_ogrid+1+i:l2_ogrid+1,m_ogrid,n_ogrid,k)+f(l1_ogrid-1+i:l2_ogrid-1,m_ogrid,n_ogrid,k)) &
                  +der2_coef2*(f(l1_ogrid+2+i:l2_ogrid+2,m_ogrid,n_ogrid,k)+f(l1_ogrid-2+i:l2_ogrid-2,m_ogrid,n_ogrid,k)) &
                  +der2_coef3*(f(l1_ogrid+3+i:l2_ogrid+3,m_ogrid,n_ogrid,k)+f(l1_ogrid-3+i:l2_ogrid-3,m_ogrid,n_ogrid,k)))
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
            call der_ogrid(f,k,df,j)
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
            call der_ogrid(f,k,df,j)
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
!  05-apr-17/Jorgen: Added summation by parts operator near cylinder surface
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid) :: df,fac
      real, dimension (6) :: facSBP
      integer :: i,j,k,ii,kk
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
        ! Note that the summation by parts operators are only implemented for this scheme
        !
        if ((i==1.and.j==2).or.(i==2.and.j==1)) then
          if (nxgrid_ogrid/=1.and.nygrid_ogrid/=1) then
            if(SBP.and.lfirst_proc_x) then
              ii=6
              do kk=1,6
                facSBP=(1./60.)*dx_1_ogrid(l1_ogrid:l1_ogrid+5)*dy_1_ogrid(m_ogrid)
                df(kk)=facSBP(kk)*( &
                  45.*(D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid+1,n_ogrid,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid+1,n_ogrid,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid+1,n_ogrid,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid+1,n_ogrid,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid+1,n_ogrid,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid+1,n_ogrid,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid+1,n_ogrid,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid+1,n_ogrid,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid+1,n_ogrid,k) &
                      -D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid-1,n_ogrid,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid-1,n_ogrid,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid-1,n_ogrid,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid-1,n_ogrid,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid-1,n_ogrid,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid-1,n_ogrid,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid-1,n_ogrid,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid-1,n_ogrid,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid-1,n_ogrid,k)) &
                  -9.*(D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid+2,n_ogrid,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid+2,n_ogrid,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid+2,n_ogrid,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid+2,n_ogrid,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid+2,n_ogrid,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid+2,n_ogrid,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid+2,n_ogrid,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid+2,n_ogrid,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid+2,n_ogrid,k) &
                      -D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid-2,n_ogrid,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid-2,n_ogrid,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid-2,n_ogrid,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid-2,n_ogrid,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid-2,n_ogrid,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid-2,n_ogrid,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid-2,n_ogrid,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid-2,n_ogrid,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid-2,n_ogrid,k)) &
                     +(D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid+3,n_ogrid,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid+3,n_ogrid,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid+3,n_ogrid,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid+3,n_ogrid,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid+3,n_ogrid,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid+3,n_ogrid,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid+3,n_ogrid,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid+3,n_ogrid,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid+3,n_ogrid,k) &
                      -D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid-3,n_ogrid,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid-3,n_ogrid,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid-3,n_ogrid,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid-3,n_ogrid,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid-3,n_ogrid,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid-3,n_ogrid,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid-3,n_ogrid,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid-3,n_ogrid,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid-3,n_ogrid,k)))
              enddo
            else
              ii=0
            endif
            fac=(1./60.**2)*dx_1_ogrid(l1_ogrid:l2_ogrid)*dy_1_ogrid(m_ogrid)
            df(1+ii:nx_ogrid)=fac(1+ii:nx_ogrid)*( &
              45.*((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid+1,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid+1,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid+1,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid+1,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid+1,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid+1,n_ogrid,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid-1,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid-1,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid-1,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid-1,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid-1,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid-1,n_ogrid,k))))&
              -9.*((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid+2,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid+2,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid+2,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid+2,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid+2,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid+2,n_ogrid,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid-2,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid-2,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid-2,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid-2,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid-2,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid-2,n_ogrid,k))))&
                 +((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid+3,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid+3,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid+3,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid+3,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid+3,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid+3,n_ogrid,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid-3,n_ogrid,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid-3,n_ogrid,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid-3,n_ogrid,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid-3,n_ogrid,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid-3,n_ogrid,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid-3,n_ogrid,k))))&
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
            if(SBP.and.lfirst_proc_x) then
              ii=6
              do kk=1,6
                facSBP=(1./60.)*dx_1_ogrid(l1_ogrid:l1_ogrid+5)*dz_1_ogrid(n_ogrid)
                df(kk)=facSBP(kk)*( &
                  45.*(D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid,n_ogrid+1,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid,n_ogrid+1,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid,n_ogrid+1,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid,n_ogrid+1,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid,n_ogrid+1,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid,n_ogrid+1,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid,n_ogrid+1,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid,n_ogrid+1,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid,n_ogrid+1,k)                                           &
                      -D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid,n_ogrid-1,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid,n_ogrid-1,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid,n_ogrid-1,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid,n_ogrid-1,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid,n_ogrid-1,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid,n_ogrid-1,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid,n_ogrid-1,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid,n_ogrid-1,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid,n_ogrid-1,k))                                          &
                  -9.*(D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid,n_ogrid+2,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid,n_ogrid+2,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid,n_ogrid+2,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid,n_ogrid+2,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid,n_ogrid+2,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid,n_ogrid+2,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid,n_ogrid+2,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid,n_ogrid+2,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid,n_ogrid+2,k)                                           &
                      -D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid,n_ogrid-2,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid,n_ogrid-2,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid,n_ogrid-2,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid,n_ogrid-2,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid,n_ogrid-2,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid,n_ogrid-2,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid,n_ogrid-2,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid,n_ogrid-2,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid,n_ogrid-2,k))                                          &
                     +(D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid,n_ogrid+3,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid,n_ogrid+3,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid,n_ogrid+3,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid,n_ogrid+3,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid,n_ogrid+3,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid,n_ogrid+3,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid,n_ogrid+3,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid,n_ogrid+3,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid,n_ogrid+3,k)                                           &
                      -D1_SBP(kk,1)*f(l1_ogrid  ,m_ogrid,n_ogrid-3,k)+D1_SBP(kk,2)*f(l1_ogrid+1,m_ogrid,n_ogrid-3,k) + &
                       D1_SBP(kk,3)*f(l1_ogrid+2,m_ogrid,n_ogrid-3,k)+D1_SBP(kk,4)*f(l1_ogrid+3,m_ogrid,n_ogrid-3,k) + &
                       D1_SBP(kk,5)*f(l1_ogrid+4,m_ogrid,n_ogrid-3,k)+D1_SBP(kk,6)*f(l1_ogrid+5,m_ogrid,n_ogrid-3,k) + &
                       D1_SBP(kk,7)*f(l1_ogrid+6,m_ogrid,n_ogrid-3,k)+D1_SBP(kk,8)*f(l1_ogrid+7,m_ogrid,n_ogrid-3,k) + &
                       D1_SBP(kk,9)*f(l1_ogrid+8,m_ogrid,n_ogrid-3,k)))
              enddo
            else
              ii=0
            endif
            fac=(1./60.**2)*dz_1_ogrid(n_ogrid)*dx_1_ogrid(l1_ogrid:l2_ogrid)
            df(1+ii:nx_ogrid)=fac(1+ii:nx_ogrid)*( &
              45.*((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid+1,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid+1,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid+1,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid-1,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid-1,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid-1,k))))&
              -9.*((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid+2,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid+2,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid+2,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid-2,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid-2,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid-2,k))))&
                 +((45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid+3,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid+3,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid+3,k))) &
                  -(45.*(f(l1_ogrid+1+ii:l2_ogrid+1,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-1+ii:l2_ogrid-1,m_ogrid,n_ogrid-3,k))  &
                    -9.*(f(l1_ogrid+2+ii:l2_ogrid+2,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-2+ii:l2_ogrid-2,m_ogrid,n_ogrid-3,k))  &
                       +(f(l1_ogrid+3+ii:l2_ogrid+3,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-3+ii:l2_ogrid-3,m_ogrid,n_ogrid-3,k))))&
                   )
            !fac=(1./60.**2)*dz_1_ogrid(n_ogrid)*dx_1_ogrid(l1_ogrid:l2_ogrid)
            !df=fac*( &
            !  45.*((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid+1,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid+1,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid+1,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid+1,k))) &
            !      -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid-1,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid-1,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid-1,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid-1,k))))&
            !  -9.*((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid+2,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid+2,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid+2,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid+2,k))) &
            !      -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid-2,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid-2,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid-2,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid-2,k))))&
            !     +((45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid+3,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid+3,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid+3,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid+3,k))) &
            !      -(45.*(f(l1_ogrid+1:l2_ogrid+1,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid-3,k))  &
            !        -9.*(f(l1_ogrid+2:l2_ogrid+2,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid-3,k))  &
            !           +(f(l1_ogrid+3:l2_ogrid+3,m_ogrid,n_ogrid-3,k)-f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid-3,k))))&
            !       )
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
    subroutine der6_ogrid(f,k,df,j,ignoredx,upwind)
!
!  Calculate 6th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!    The optional flag UPWIND is a variant thereof, which calculates
!  D^(6)*dx^5/60, which is the upwind correction of centered derivatives.
!
!   27-feb-17/Jorgen: Adapted from deriv.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid) :: df,fac
      integer :: j,k,i=0
      logical, optional :: ignoredx,upwind
      logical :: igndx,upwnd
!
      intent(in)  :: f,k,j,ignoredx,upwind
      intent(out) :: df
!
      if (present(ignoredx)) then
        igndx = ignoredx
      else
        if (.not. lequidist_ogrid(j)) then
          call fatal_error('der6','for non-equidistant grid only '//&
              'if dx is ignored.')
          igndx = .true.
        endif
        igndx = .false.
      endif
!
      if (present(upwind)) then
        if (.not. lequidist_ogrid(j)) then
          call fatal_error('der6','upwind cannot be used with '//&
              'non-equidistant grid.')
        endif
        upwnd = upwind
      else
        upwnd = .false.
!        if ((.not.lcartesian_coords).and.(.not.igndx)) then
!DM: non cartesian grids should not necessarily use upwinding. Wlad do you disagree ?
!         if (.not.igndx) then
!          call fatal_error('der6','in non-cartesian coordinates '//&
!              'just works if upwinding is used')
!        endif
      endif
!
      if(SBP.and.lfirst_proc_x) i=6
        
      if (j==1) then
        if (nxgrid_ogrid/=1) then
          if (igndx) then
            fac=1.
          else if (upwnd) then
            fac=(1.0/60)*dx_1_ogrid(l1_ogrid:l2_ogrid)
          else
            fac=dx_1_ogrid(l1_ogrid:l2_ogrid)**6
          endif
          df(1:i+1)=0
          df(1+i:nxgrid)=fac(1+i:nxgrid)*(- 20.0* f(l1_ogrid+i:l2_ogrid,m_ogrid,n_ogrid,k) &
                        + 15.0*(f(l1_ogrid+1+i:l2_ogrid+1,m_ogrid,n_ogrid,k)+f(l1_ogrid-1:l2_ogrid-1,m_ogrid,n_ogrid,k)) &
                        -  6.0*(f(l1_ogrid+2+i:l2_ogrid+2,m_ogrid,n_ogrid,k)+f(l1_ogrid-2:l2_ogrid-2,m_ogrid,n_ogrid,k)) &
                        +      (f(l1_ogrid+3+i:l2_ogrid+3,m_ogrid,n_ogrid,k)+f(l1_ogrid-3:l2_ogrid-3,m_ogrid,n_ogrid,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid_ogrid/=1) then
          if (igndx) then
            fac=1.
          else if (upwnd) then
            fac=(1.0/60)*dy_1_ogrid(m_ogrid)
          else
            fac=dy_1_ogrid(m_ogrid)**6
          endif
          df=fac*(- 20.0* f(l1_ogrid:l2_ogrid,m_ogrid  ,n_ogrid,k) &
                  + 15.0*(f(l1_ogrid:l2_ogrid,m_ogrid+1,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-1,n_ogrid,k)) &
                  -  6.0*(f(l1_ogrid:l2_ogrid,m_ogrid+2,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-2,n_ogrid,k)) &
                  +      (f(l1_ogrid:l2_ogrid,m_ogrid+3,n_ogrid,k)+f(l1_ogrid:l2_ogrid,m_ogrid-3,n_ogrid,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid_ogrid/=1) then
          if (igndx) then
            fac=1.
          else if (upwnd) then
            fac=(1.0/60)*dz_1_ogrid(n_ogrid)
          else
            fac=dz_1_ogrid(n_ogrid)**6
          endif
          df=fac*(- 20.0* f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid  ,k) &
                  + 15.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+1,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-1,k)) &
                  -  6.0*(f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+2,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-2,k)) &
                  +      (f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid+3,k)+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid-3,k)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der6_ogrid
!***********************************************************************
    subroutine deri_3d_inds_ogrid(f,df,inds,j,lignored,lnometric)
!
!  dummy routine for compatibility
!
!  27-feb-17/Jorgen: Adapted from deriv.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
      real, dimension (nx_ogrid)                           :: df
      integer                             :: j
      logical,                   optional :: lignored, lnometric
      integer, dimension(nx_ogrid)              :: inds
!
      intent(in)  :: f,j,inds,lignored,lnometric
      intent(out) :: df
!
      call fatal_error('deri_3d_inds_ogrid','Upwinding not implemented for nonuniform grids')
!
! dummy computation to avoid compiler warnings of unused variables
      if (present(lignored).and.present(lnometric)) &
          df  = inds + f(l1_ogrid:l2_ogrid,1,1,1) + j
!
    endsubroutine deri_3d_inds_ogrid
!************************************************************************
    subroutine set_ghosts_onesided_ogrid(ivar)
!
!   Set ghost points for onesided boundary conditions with Dirichlet BC
!   on the cylidner surface.
!   Only works for the radial direction.
!
!   16-feb-17/Jorgen: Adapted from deriv.f90.
!
      integer :: k,ivar,i
!
      do i=1,nghost
        k=l1_ogrid-i
        f_ogrid(k,:,:,ivar)=7*f_ogrid(k+1,:,:,ivar) &
                          -21*f_ogrid(k+2,:,:,ivar) &
                          +35*f_ogrid(k+3,:,:,ivar) &
                          -35*f_ogrid(k+4,:,:,ivar) &
                          +21*f_ogrid(k+5,:,:,ivar) &
                           -7*f_ogrid(k+6,:,:,ivar) &
                             +f_ogrid(k+7,:,:,ivar)
      enddo

    endsubroutine set_ghosts_onesided_ogrid
!***********************************************************************
    subroutine bval_from_neumann_arr_ogrid
!
!  Calculates the boundary value from the Neumann BC d f/d x_i = val employing
!  one-sided difference formulae. val depends on x,y.
!
!  16-feb-17/Jorgen: Adapted from deriv.f90
!
      real :: val=0.
      integer :: k

      k=l1_ogrid
      f_ogrid(k,:,:,ilnrho) = (-val*60.*dx_ogrid + 360.*f_ogrid(k+1,:,:,ilnrho) &
                                                 - 450.*f_ogrid(k+2,:,:,ilnrho) &
                                                 + 400.*f_ogrid(k+3,:,:,ilnrho) &
                                                 - 225.*f_ogrid(k+4,:,:,ilnrho) &
                                                 +  72.*f_ogrid(k+5,:,:,ilnrho) &
                                                 -  10.*f_ogrid(k+6,:,:,ilnrho) )/147.

    endsubroutine bval_from_neumann_arr_ogrid
!***********************************************************************
    subroutine bval_from_neumann_arr_ogrid_alt
!
!  Calculates the boundary value from the Neumann BC d f/d x_i = val employing
!  one-sided difference formulae. val depends on x,y.
!
!  16-feb-17/Jorgen: Adapted from deriv.f90
!
      real :: val=0.
      integer :: k

      k=l1_ogrid
      f_ogrid(k,:,:,irho) = -(D1_SBP(1,2)*f_ogrid(k+1,:,:,irho) + &
                              D1_SBP(1,3)*f_ogrid(k+2,:,:,irho) + &
                              D1_SBP(1,4)*f_ogrid(k+3,:,:,irho) + &
                              D1_SBP(1,5)*f_ogrid(k+4,:,:,irho) + &
                              D1_SBP(1,6)*f_ogrid(k+5,:,:,irho) + &
                              D1_SBP(1,7)*f_ogrid(k+6,:,:,irho) + &
                              D1_SBP(1,8)*f_ogrid(k+7,:,:,irho) + &
                              D1_SBP(1,9)*f_ogrid(k+8,:,:,irho) )/D1_SBP(1,1)

    endsubroutine bval_from_neumann_arr_ogrid_alt
!***********************************************************************
    subroutine gaunoise_ogrid(ampl,i1,i2)
!
!  Add Gaussian noise (= normally distributed) white noise for variables i1:i2
!
!  13-feb-17/Jorgen: Adapted from gaunoise_vect in initcond.f90
!
      use General, only: random_number_wrapper
      real :: ampl
      integer :: i1,i2
!
      real, dimension (mx_ogrid) :: r,p,tmp
      integer :: i
!
      intent(in)    :: ampl,i1,i2
!
!  set gaussian random noise vector
!
      if (ampl==0) then
        if (lroot) print*,'gaunoise_ogrid: ampl=0 for i1,i2=',i1,i2
      else
        if ((ip<=8).and.lroot) print*,'gaunoise_ogrid: i1,i2=',i1,i2
        do n=1,mz_ogrid; do m=1,my_ogrid
          do i=i1,i2
            if (lroot.and.m==1.and.n==1) print*,'gaunoise_ogrid: variable i=',i
            if (modulo(i-i1,2)==0) then
              call random_number_wrapper(r)
              call random_number_wrapper(p)
              tmp=sqrt(-2*log(r))*sin(2*pi*p)
            else
              tmp=sqrt(-2*log(r))*cos(2*pi*p)
            endif
            f_ogrid(:,m_ogrid,n_ogrid,i)=f_ogrid(:,m_ogrid,n_ogrid,i)+ampl*tmp
          enddo
        enddo; enddo
      endif
!
    endsubroutine gaunoise_ogrid
!***********************************************************************
    subroutine find_proc_cartesian(xyz,from_proc)
!
!  Find the processor that stores the grid points, and return the processor
!  id and the grid index of the bottom neighbouring point on a global grid.
!  Necessary for interpolation between grids on parallel systems
!
!  13-apr-17/Jorgen: Coded
!
      real, dimension(3), intent(in) :: xyz
      integer, intent(out) :: from_proc
      integer :: low_i_global,low_j_global,low_k_global
      integer :: i
      logical :: found_proc=.false.
!
      do i=1,ncpus
        if( ((xyz(1)>=xyz0_loc_all(i,1)).and.(xyz(1)<=xyz1_loc_all(i,1))) .and. &
            ((xyz(2)>=xyz0_loc_all(i,2)).and.(xyz(2)<=xyz1_loc_all(i,2))) .and. &
            ((xyz(3)>=xyz0_loc_all(i,3)).and.(xyz(3)<=xyz1_loc_all(i,3))) )  then
            from_proc=i-1
            found_proc=.true.
            exit
        endif
      enddo
      if(.not.found_proc) then
        print*, 'find_proc_cartesian: error when searching for interpolation point'
        print*, 'find_proc_cartesian: x,y,z',xyz
        print*, 'find_proc_cartesian: x0_loc_all',xyz0_loc_all(:,1)
        print*, 'find_proc_cartesian: x1_loc_all',xyz1_loc_all(:,1)
        print*, 'find_proc_cartesian: y0_loc_all',xyz0_loc_all(:,2)
        print*, 'find_proc_cartesian: y1_loc_all',xyz1_loc_all(:,2)
        print*, 'find_proc_cartesian: z0_loc_all',xyz0_loc_all(:,3)
        print*, 'find_proc_cartesian: z1_loc_all',xyz1_loc_all(:,3)
        call fatal_error('find_proc_cartesian', &
          'could not locate interpolation point on any processor!')
      endif
!
    endsubroutine find_proc_cartesian
!***********************************************************************
    subroutine find_proc_curvilinear(rthz,from_proc)
!
!  Find the processor that stores the grid points, and return the processor
!  id and the grid index of the bottom neighbouring point on a global grid.
!  Necessary for interpolation between grids on parallel systems
!
!  13-apr-17/Jorgen: Coded
!
      real, dimension(3), intent(in) :: rthz
      integer, intent(out) :: from_proc
      integer :: low_i_global,low_j_global,low_k_global
      integer :: i
      logical :: found_proc=.false.
!
      do i=1,ncpus
        if( ((rthz(1)>=xyz0_loc_all_ogrid(i,1)).and.(rthz(1)<=xyz1_loc_all_ogrid(i,1))) .and. &
            ((rthz(2)>=xyz0_loc_all_ogrid(i,2)).and.(rthz(2)<=xyz1_loc_all_ogrid(i,2))) .and. &
            ((rthz(3)>=xyz0_loc_all_ogrid(i,3)).and.(rthz(3)<=xyz1_loc_all_ogrid(i,3))) )  then
            from_proc=i-1
            found_proc=.true.
            exit
        endif
      enddo
      if(.not.found_proc) call fatal_error('find_proc_curvilinear', &
          'could not locate interpolation point on any processor!')
!
    endsubroutine find_proc_curvilinear
!***********************************************************************
    subroutine construct_serial_bdry_cartesian
!
!  Build arrays containing cartesian corner values of all processors
!  The arrays xyz0_loc_all and xyz1_loc_all are accessed by
!  (iproc+1,[1,2,3]), where [1,2,3] is [x,y,z] corner.
!  Need to use iproc+1 instead of iproc to avoid accessing zeroth element.
!
!  13-apr-17/Jorgen: Coded
!
      use Mpicomm, only: mpisend_int, mpisend_real, mpirecv_int, &
                         mpirecv_real, mpibcast_real
      real, dimension(3) :: xyz0_loc_recv
      real, dimension(3) :: xyz1_loc_recv
      integer, dimension(2) :: nbcast=(/ ncpus,3 /)
      integer :: iproc_recv,j
!
      if (iproc/=root) then
!
!  All processors send their array values to the root.
!
        call mpisend_int(iproc,root,990)
        call mpisend_real(xyz0_loc,3,root,991)
        call mpisend_real(xyz1_loc,3,root,992)
      else
!
!  The root processor, in turn, receives the data from the others
!
        do j=0,ncpus-1
        !avoid send-to-self
          if (j/=root) then
            call mpirecv_int(iproc_recv,j,990)
            call mpirecv_real(xyz0_loc_recv,3,iproc_recv,991)
            call mpirecv_real(xyz1_loc_recv,3,iproc_recv,992)
!
            xyz0_loc_all(iproc_recv+1,:)=xyz0_loc_recv
            xyz1_loc_all(iproc_recv+1,:)=xyz1_loc_recv
          else
!  The root just copies its value to the serial array
            xyz0_loc_all(root+1,:)=xyz0_loc
            xyz1_loc_all(root+1,:)=xyz1_loc
          endif
        enddo
      endif
!
!  Serial array constructed. Broadcast the result. 
!
      call mpibcast_real(xyz0_loc_all,nbcast)
      call mpibcast_real(xyz1_loc_all,nbcast)
    endsubroutine construct_serial_bdry_cartesian
!***********************************************************************
    subroutine construct_serial_bdry_curv
!
!  Build arrays containing curvilinear corner values of all processors
!  The arrays xyz0_loc_all and xyz1_loc_all are accessed by
!  (iproc+1,[1,2,3]), where [1,2,3] is [x,y,z] corner.
!  Need to use iproc+1 instead of iproc to avoid accessing zeroth element.
!
!  Unlike the cartesian version of this, we need to first construct the
!  local arrays xyz0_loc_ogrid and xyz1_loc_ogrid. This is done in
!  start.in/run.in for the cartesian grid.
!
!  13-apr-17/Jorgen: Coded
!
      use Mpicomm, only: mpisend_int, mpisend_real, mpirecv_int, &
                         mpirecv_real, mpibcast_real
      real, dimension(3) :: xyz0_loc_recv_ogrid
      real, dimension(3) :: xyz1_loc_recv_ogrid
      real, dimension(3) :: Lxyz_loc_ogrid
      integer, dimension(2) :: nbcast=(/ ncpus,3 /)
      integer :: iproc_recv,j
!
!  Constructing local arrays, with code copied from run.f90
!  Size of box at local processor. The if-statement is for
!  backward compatibility.
!
      if (lequidist_ogrid(1)) then
        Lxyz_loc_ogrid(1) = Lxyz_ogrid(1)/nprocx
        xyz0_loc_ogrid(1) = xyz0_ogrid(1)+ipx*Lxyz_loc_ogrid(1)
        xyz1_loc_ogrid(1) = xyz0_loc_ogrid(1)+Lxyz_loc_ogrid(1)
      else
!
!  In the equidistant grid, the processor boundaries (xyz[01]_loc) do NOT
!  coincide with the l[mn]1[2] points. Also, xyz0_loc[ipx+1]=xyz1_loc[ipx], i.e.,
!  the inner boundary of one is exactly the outer boundary of the other. Reproduce
!  this behavior also for non-equidistant grids.
!
        if (ipx==0) then
          xyz0_loc_ogrid(1) = x_ogrid(l1_ogrid)
        else
          xyz0_loc_ogrid(1) = x_ogrid(l1_ogrid) - .5/dx_1_ogrid(l1_ogrid)
        endif
        if (ipx==nprocx-1) then
          xyz1_loc_ogrid(1) = x_ogrid(l2_ogrid)
        else
          xyz1_loc_ogrid(1) = x_ogrid(l2_ogrid+1) - .5/dx_1_ogrid(l2_ogrid+1)
        endif
        Lxyz_loc_ogrid(1) = xyz1_loc_ogrid(1) - xyz0_loc_ogrid(1)
      endif
!
      if (lequidist_ogrid(2)) then
        Lxyz_loc_ogrid(2) = Lxyz_ogrid(2)/nprocy
        xyz0_loc_ogrid(2) = xyz0_ogrid(2)+ipy*Lxyz_loc_ogrid(2)
        xyz1_loc_ogrid(2) = xyz0_loc_ogrid(2)+Lxyz_loc_ogrid(2)
      else
        if (ipy==0) then
          xyz0_loc_ogrid(2) = y_ogrid(m1_ogrid)
        else
          xyz0_loc_ogrid(2) = y_ogrid(m1_ogrid) - .5/dy_1_ogrid(m1_ogrid)
        endif
        if (ipy==nprocy-1) then
          xyz1_loc_ogrid(2) = y_ogrid(m2_ogrid)
        else
          xyz1_loc_ogrid(2) = y_ogrid(m2_ogrid+1) - .5/dy_1_ogrid(m2_ogrid+1)
        endif
        Lxyz_loc_ogrid(2) = xyz1_loc_ogrid(2) - xyz0_loc_ogrid(2)
      endif
!
      if (lequidist_ogrid(3)) then 
        Lxyz_loc_ogrid(3) = Lxyz_ogrid(3)/nprocz
        xyz0_loc_ogrid(3) = xyz0_ogrid(3)+ipz*Lxyz_loc_ogrid(3)
        xyz1_loc_ogrid(3) = xyz0_loc_ogrid(3)+Lxyz_loc_ogrid(3)
      else
        if (ipz==0) then
          xyz0_loc_ogrid(3) = z_ogrid(n1_ogrid) 
        else
          xyz0_loc_ogrid(3) = z_ogrid(n1_ogrid) - .5/dz_1_ogrid(n1_ogrid)
        endif
        if (ipz==nprocz-1) then
          xyz1_loc_ogrid(3) = z_ogrid(n2_ogrid)
        else
          xyz1_loc_ogrid(3) = z_ogrid(n2_ogrid+1) - .5/dz_1_ogrid(n2_ogrid+1)
        endif
        Lxyz_loc_ogrid(3) = xyz1_loc_ogrid(3) - xyz0_loc_ogrid(3)
      endif
!
!  Communicate arrays and generate global arrays
!
      if (iproc/=root) then
!
!  All processors send their array values to the root.
!
        call mpisend_int(iproc,root,880)
        call mpisend_real(xyz0_loc_ogrid,3,root,881)
        call mpisend_real(xyz1_loc_ogrid,3,root,882)
      else
!
!  The root processor, in turn, receives the data from the others
!
        do j=0,ncpus-1
        !avoid send-to-self
          if (j/=root) then
            call mpirecv_int(iproc_recv,j,880)
            call mpirecv_real(xyz0_loc_recv_ogrid,3,iproc_recv,881)
            call mpirecv_real(xyz1_loc_recv_ogrid,3,iproc_recv,882)
!
            xyz0_loc_all_ogrid(iproc_recv+1,:)=xyz0_loc_recv_ogrid
            xyz1_loc_all_ogrid(iproc_recv+1,:)=xyz1_loc_recv_ogrid
          else
!  The root just copies its value to the serial array
            xyz0_loc_all_ogrid(root+1,:)=xyz0_loc_ogrid
            xyz1_loc_all_ogrid(root+1,:)=xyz1_loc_ogrid
          endif
        enddo
      endif
!
!  Serial array constructed. Broadcast the result. 
!
      call mpibcast_real(xyz0_loc_all_ogrid,nbcast)
      call mpibcast_real(xyz1_loc_all_ogrid,nbcast)
    endsubroutine construct_serial_bdry_curv
!***********************************************************************
    subroutine find_near_cartesian(xyz,xyz_neighbours)
!
!  Find the grid point values of all neighbouring points on cartesian mesh
!  Return array containing a cube with element 1 in the bottom left corner
!  (low z-plane) and element 8 in the top right corner (top z-plane)
!
!  12-apr-17/Jorgen: Coded
!
      real, dimension(3), intent(in)    :: xyz
      real, dimension(8,3), intent(out)  :: xyz_neighbours
      integer :: lower_i, upper_i, lower_j, upper_j, lower_k, upper_k
      integer :: i,j,k,ii=0
!
      call find_near_cartesian_indices(lower_i,upper_i,lower_j,upper_j, &
        lower_k,upper_k,xyz)
!
      do i=lower_i,upper_i
        do j=lower_j,upper_j
          do k=lower_k,upper_k
            ii=ii+1
            xyz_neighbours(ii,:)= (/ x(i),y(j),z(k) /)
          enddo
        enddo
      enddo
    endsubroutine find_near_cartesian
!***********************************************************************
    subroutine find_near_curvilinear(rthz,rthz_neighbours)
!
!  Find the grid point values of all neighbouring points on curvilinear mesh.
!  Return array containing a (curvilinear) cube with element 1 in the bottom left corner
!  (low z-plane) and element 8 in the top right corner (top z-plane)
!
!  12-apr-17/Jorgen: Coded
!
      real, dimension(3), intent(in)    :: rthz
      real, dimension(8,3), intent(out)  :: rthz_neighbours
      integer :: lower_i, upper_i, lower_j, upper_j, lower_k, upper_k
      integer :: i,j,k,ii=0
!
      call find_near_curvilinear_indices(lower_i,upper_i,lower_j,upper_j, &
        lower_k,upper_k,rthz)
!
      do i=lower_i,upper_i
        do j=lower_j,upper_j
          do k=lower_k,upper_k
            ii=ii+1
            rthz_neighbours(ii,:)= (/ x_ogrid(i),y_ogrid(j),z_ogrid(k) /)
          enddo
        enddo
      enddo
    endsubroutine find_near_curvilinear
!***********************************************************************
    subroutine find_near_cartesian_indices(lower_i,upper_i,lower_j,upper_j, &
        lower_k,upper_k,xyz)
!
!  Find i, j and k indices for all neighbouring grid points
!
      integer :: ii, jj, kk
      integer, intent(out) :: lower_i, upper_i, lower_j, upper_j, lower_k, upper_k
      real, intent(in), dimension(3)  :: xyz
!
      lower_i = 0
      upper_i = 0
      do ii = 2,mx
        if (x(ii) > xyz(1)) then
          lower_i = ii-1
          upper_i = ii
          exit
        endif
      enddo
!
      lower_j = 0
      upper_j = 0
      do jj = 2,my
        if (y(jj) > xyz(2)) then
          lower_j = jj-1
          upper_j = jj
          exit
        endif
      enddo
!
      if (nzgrid == 1) then
        lower_k = n1
        upper_k = n1
      else
        lower_k = 0
        upper_k = 0
        do kk = 2,mz
          if (z(kk) > xyz(3)) then
            lower_k = kk-1
            upper_k = kk
            exit
          endif
        enddo
      endif
!
    endsubroutine find_near_cartesian_indices
!***********************************************************************
    subroutine find_near_cart_ind_global(lower_indices,xyz)
!
!  Find i, j and k indices for lower neighbouring grid point on global grid
!
!  13-apr-17/Jorgen: Adapted from find_near_cartesian_indices
!
      integer :: ii, jj, kk
      integer, dimension(3), intent(out) :: lower_indices 
      real, dimension(3), intent(in)  :: xyz
!
!
      lower_indices = 0
      do ii = 1+nghost,mxgrid-nghost
        if (xglobal(ii) > xyz(1)) then
          lower_indices(1) = ii-1
          exit
        endif
      enddo
!
      do jj = 1+nghost,mygrid-nghost
        if (yglobal(jj) > xyz(2)) then
          lower_indices(2) = jj-1
          exit
        endif
      enddo
!
      if (nzgrid == 1) then
        lower_indices(3) = n1
      else
        do kk = 1+nghost,mzgrid-nghost
          if (zglobal(kk) > xyz(3)) then
            lower_indices(3) = kk-1
            exit
          endif
        enddo
      endif
!
    endsubroutine find_near_cart_ind_global
!***********************************************************************
    subroutine find_near_curvilinear_indices(lower_i,upper_i,lower_j,upper_j, &
        lower_k,upper_k,rthz)
!
!  Find i, j and k indices for all neighbouring grid points
!
      integer :: ii, jj, kk
      integer, intent(out) :: lower_i, upper_i, lower_j, upper_j, lower_k, upper_k
      real, intent(in), dimension(3)  :: rthz
!
      lower_i = 0
      upper_i = 0
      do ii = 2,mx_ogrid
        if (x_ogrid(ii) > rthz(1)) then
          lower_i = ii-1
          upper_i = ii
          exit
        endif
      enddo
!
      lower_j = 0
      upper_j = 0
      do jj = 2,my_ogrid
        if (y_ogrid(jj) > rthz(2)) then
          lower_j = jj-1
          upper_j = jj
          exit
        endif
      enddo
!
      if (nzgrid_ogrid == 1) then
        lower_k = n1
        upper_k = n1
      else
        lower_k = 0
        upper_k = 0
        do kk = 2,mz_ogrid
          if (z_ogrid(kk) > rthz(3)) then
            lower_k = kk-1
            upper_k = kk
            exit
          endif
        enddo
      endif
!
    endsubroutine find_near_curvilinear_indices
!***********************************************************************
    subroutine find_near_curv_ind_global(lower_indices,rthz)
!
!  Find i, j and k indices for lower neighbouring grid point on global grid
!
!  13-apr-17/Jorgen: Adapted from find_near_curvilinear_indices
!
      integer :: ii, jj, kk
      integer, dimension(3), intent(out) :: lower_indices
      real, dimension(3), intent(in) :: rthz
!
      lower_indices = 0
      do ii = 1,mxgrid_ogrid
        if (xglobal_ogrid(ii) > rthz(1)) then
          lower_indices(1) = ii-1
          exit
        endif
      enddo
!
      do jj = 1,mygrid_ogrid
        if (yglobal_ogrid(jj) > rthz(2)) then
          lower_indices(2) = jj-1
          exit
        endif
      enddo
!
      if (nzgrid == 1) then
        lower_indices(3) = n1
      else
        do kk = 1,mzgrid_ogrid
          if (zglobal_ogrid(kk) > rthz(3)) then
            lower_indices(3) = kk-1
            exit
          endif
        enddo
      endif
!
    endsubroutine find_near_curv_ind_global
!***********************************************************************
    subroutine wsnap_ogrid(chsnap,enum,flist)
!
!  Write snapshot file of overlapping grid, labelled consecutively if enum==.true.
!  Otherwise just write a snapshot without label (used for var.dat).
!
!  21-feb-17/Jorgen: Adapted from snapsjots.f90
!
      use General, only: safe_character_assign
      use IO, only: log_filename_to_file
      use Sub, only: read_snaptime, update_snaptime
!
      character(len=*), intent(in) :: chsnap
      character(len=*), intent(in), optional :: flist
      logical, intent(in), optional :: enum
!
      real, save :: tsnap
      integer, save :: nsnap
      logical, save :: lfirst_call=.true.
      logical :: enum_, lsnap
      character (len=fnlen) :: file
      character (len=intlen) :: ch
!
      if (present(enum)) then
        enum_=enum
      else
        enum_=.false.
      endif
!
!  Output snapshot with label in 'tsnap' time intervals.
!  File keeps the information about number and time of last snapshot.
!
      if (enum_) then
        call safe_character_assign(file,trim(datadir)//'/ogtsnap.dat')
!
!  At first call, need to initialize tsnap.
!  tsnap calculated in read_snaptime, but only available to root processor.
!
        if (lfirst_call) then
          call read_snaptime(file,tsnap,nsnap,dsnap,t)
          lfirst_call=.false.
        endif
!
!  Check whether we want to output snapshot. If so, then
!  update ghost zones for var.dat (cheap, since done infrequently).
!
        call update_snaptime(file,tsnap,nsnap,dsnap,t,lsnap,ch)
        if (lsnap) then
          call update_ghosts_ogrid
          call safe_character_assign(file,trim(chsnap)//ch)
          call output_snap_ogrid(f_ogrid,file=file)
          if (ip<=10.and.lroot) print*,'wsnap: written snapshot ',file
          if (present(flist)) call log_filename_to_file(file,flist)
        endif
!
      else
!
!  Write snapshot without label (typically, var.dat).
!
        call update_ghosts_ogrid
        call safe_character_assign(file,trim(chsnap))
        call output_snap_ogrid(f_ogrid,file=file)
        if (present(flist)) call log_filename_to_file(file,flist)
      endif
!
      if (lformat) call output_snap_form_ogrid (file)
!
    endsubroutine wsnap_ogrid
!***********************************************************************
    subroutine update_ghosts_ogrid
!
!  Update all ghost zones of f_ogrid.
!  Initiate communication and do boundary conditions.
!  Required order:
!  1. x-boundaries (x-ghost zones will be communicated) 
!  2. communication
!  3. y- and z-boundaries
!
!  21-feb-17/Jorgen: Adapted from boundcond.f90
!
      use Solid_Cells_Mpicomm, only: initiate_isendrcv_bdry_ogrid, finalize_isendrcv_bdry_ogrid
!
      call boundconds_x_ogrid
      call initiate_isendrcv_bdry_ogrid(f_ogrid)
      call finalize_isendrcv_bdry_ogrid(f_ogrid)
!  Since only periodic implementation of boundconds in y- and z-dir, call only
!  if single processor in said direction. Otherwise, handled in MPI-communication.
      if (nprocy==1)                  call boundconds_y_ogrid
      if ((nprocz==1).and.(nzgrid>1)) call boundconds_z_ogrid

!
    endsubroutine update_ghosts_ogrid
!***********************************************************************
    subroutine output_snap_ogrid(a,file)
!
!  Write snapshot file, always write time and mesh, could add other things.
!
!  21-feb-17/Jorgen: Adapted from io_dist.f90
!
      use Mpicomm, only: start_serialize, end_serialize
      use IO, only: lun_output
      use File_io, only: delete_file
!
      real, dimension (:,:,:,:),  intent(IN) :: a
      character (len=*), optional,intent(IN) :: file
!
      real :: t_sp   ! t in single precision for backwards compatibility
!
      t_sp = t
      if (lserial_io) call start_serialize
      if (present(file)) then
        call delete_file(trim(directory_snap)//'/'//file)
        open (lun_output, FILE=trim(directory_snap)//'/'//file, FORM='unformatted', status='new')
      endif
!
      if (lwrite_2d) then
        if (nz_ogrid == 1) then
          write (lun_output) a(:,:,n1_ogrid,:)
        else
          call fatal_error('output_snap_ogrid','lwrite_2d used for simulation with nz_ogri/=1!')
        endif
      else
        write (lun_output) a
      endif
      write (lun_output) t_sp, x_ogrid(1:size(a,1)), y_ogrid(1:size(a,2)), z_ogrid(1:size(a,3)), dx_ogrid, dy_ogrid, dz_ogrid
!
      close (lun_output)
      if (lserial_io) call end_serialize

    endsubroutine output_snap_ogrid
!***********************************************************************
    subroutine output_snap_form_ogrid(file)
!
!  Write FORMATTED snapshot file
!
!  21/feb-17/Jorgen: Adapted from snapshot.f90
!
      use IO, only: lun_output
!
      character (len=*), intent(in) :: file
      integer :: i, j, k
!
      open(lun_output,FILE=trim(directory_dist)//trim(file)//'.form')
!
      if (lwrite_2d) then
        if (nz_ogrid==1) then
          do i = l1_ogrid, l2_ogrid
            do j = m1_ogrid, m2_ogrid
              write(lun_output,'(40(f12.5))') x_ogrid(i),y_ogrid(j),z_ogrid(n1), &
                    dx_ogrid,dy_ogrid,dz_ogrid,f_ogrid(i,j,n1_ogrid,:)
            enddo
          enddo
        else
          call fatal_error('output_snap_form_ogrid','lwrite_2d used for simulation with nz_ogri/=1!')
        endif
!
      else
        do i = l1_ogrid, l2_ogrid
          do j = m1_ogrid, m2_ogrid
            do k = n1_ogrid, n2_ogrid
              write(lun_output,'(40(f12.5))') x_ogrid(i),y_ogrid(j),z_ogrid(k), &
                    dx_ogrid,dy_ogrid,dz_ogrid,f_ogrid(i,j,k,:)
            enddo
          enddo
        enddo
!
      endif
!
      close(lun_output)
!
    endsubroutine output_snap_form_ogrid
!***********************************************************************
    subroutine rsnap_ogrid(chsnap,lread_nogrid)
!
!  Read snapshot file.
!
!  21-feb-17/Jorgen: Adapted from snapshot.f90
!
      use IO, only: lun_input
      use Mpicomm, only: end_serialize
!
      logical :: lread_nogrid
      integer :: mode
      character (len=*) :: chsnap
!
!  possibility of not reading the mesh data nor the time
!  of the snapshot. The mesh information is then taken from
!  proc*/mesh.dat
!
      if (lread_nogrid) then
        mode=0
      else
        mode=1
      endif
!
      call input_snap_ogrid(chsnap,f_ogrid,mfarray,mode)
      close (lun_input)
      if (lserial_io) call end_serialize
!
!  Read data using lnrho, and now convert to rho.
!  This assumes that one is now using ldensity_nolog=T.
!
      if (lread_oldsnap_lnrho2rho) then
        print*,'convert lnrho -> rho',ilnrho,irho
        if (irho>0) &
          f_ogrid(:,:,:,irho)=exp(f_ogrid(:,:,:,ilnrho))
      endif
!
    endsubroutine rsnap_ogrid
!***********************************************************************
    subroutine input_snap_ogrid(file,a,nv,mode)
!
!  manages reading of snapshot from different precision
!
!  21-feb-17/Jorgen: Adapted from io_dist.f90
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv, mode
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,nv), intent(out) :: a

      real(KIND=rkind8), dimension(:,:,:,:), allocatable :: adb
      real(KIND=rkind4), dimension(:,:,:,:), allocatable :: asg

      real(KIND=rkind8), dimension(:), allocatable :: xdb,ydb,zdb
      real(KIND=rkind4), dimension(:), allocatable :: xsg,ysg,zsg

      real(KIND=rkind8) :: dxdb,dydb,dzdb,deltaydb
      real(KIND=rkind4) :: dxsg,dysg,dzsg,deltaysg
      real :: deltay_ogrid

      if (lread_from_other_prec) then
        if (kind(a)==rkind4) then
          allocate(adb(mx_ogrid,my_ogrid,mz_ogrid,nv),xdb(mx_ogrid),ydb(my_ogrid),zdb(mz_ogrid))
          call read_snap_ogrid(file,adb,xdb,ydb,zdb,dxdb,dydb,dzdb,deltaydb,nv,mode)
          a=adb; x_ogrid=xdb; y_ogrid=ydb; z_ogrid=zdb; dx_ogrid=dxdb; dy_ogrid=dydb; dz_ogrid=dzdb; deltay_ogrid=deltaydb
        elseif (kind(a)==rkind8) then
          allocate(asg(mx_ogrid,my_ogrid,mz_ogrid,nv),xsg(mx_ogrid),ysg(my_ogrid),zsg(mz_ogrid))
          call read_snap_ogrid(file,asg,xsg,ysg,zsg,dxsg,dysg,dzsg,deltaysg,nv,mode)
          a=asg; x_ogrid=xsg; y_ogrid=ysg; z_ogrid=zsg; dx_ogrid=dxsg; dy_ogrid=dysg; dz_ogrid=dzsg; deltay_ogrid=deltaysg
        endif
      else
        call read_snap_ogrid(file,a,x_ogrid,y_ogrid,z_ogrid,dx_ogrid,dy_ogrid,dz_ogrid,deltay_ogrid,nv,mode)
      endif

    endsubroutine input_snap_ogrid
!***********************************************************************
    subroutine read_snap_single_ogrid(file,a,x,y,z,dx,dy,dz,deltay,nv,mode)
!
!  Read snapshot file in single precision, possibly with mesh and time (if mode=1).
!
!  21-feb-17/Jorgen: Adapted from io_dist.f90
!
      use Mpicomm, only: start_serialize, mpibcast_real, mpiallreduce_or, &
                         stop_it, mpiallreduce_min, mpiallreduce_max, MPI_COMM_WORLD
      use IO, only: lun_input
!
      character (len=*), intent(in) :: file
      real :: deltay_ogrid=0.0
      integer, intent(in) :: nv, mode
      real(KIND=rkind4), dimension (mx_ogrid,my_ogrid,mz_ogrid,nv), intent(out) :: a
!
      real(KIND=rkind4) :: t_sp, t_sgl

      real(KIND=rkind4),                       intent(out) :: dx, dy, dz, deltay
      real(KIND=rkind4), dimension (mx_ogrid), intent(out) :: x
      real(KIND=rkind4), dimension (my_ogrid), intent(out) :: y
      real(KIND=rkind4), dimension (mz_ogrid), intent(out) :: z

      real :: t_test   ! t in single precision for backwards compatibility

      logical :: ltest
      ! set ireset_tstart to 1 or 2 to coordinate divergent timestamp
      integer :: MINT=1
      integer :: MAXT=2
!
      if (lserial_io) call start_serialize
      open (lun_input, FILE=trim(directory_snap)//'/'//file, FORM='unformatted', status='old')
      if (lwrite_2d) then
        if (nz == 1) then
          read (lun_input) a(:,:,4,:)
        else
          call fatal_error ('read_snap_single_ogrid','lwrite_2d used for simulation with nz_ogri/=1!')
        endif
      else
!
!  Possibility of reading data with different numbers of ghost zones.
!  In that case, one must regenerate the mesh with luse_oldgrid=T.
!
        if (nghost_read_fewer==0) then
          read (lun_input) a
        elseif (nghost_read_fewer>0) then
          read (lun_input) &
              a(1+nghost_read_fewer:mx_ogrid-nghost_read_fewer, &
                1+nghost_read_fewer:my_ogrid-nghost_read_fewer, &
                1+nghost_read_fewer:mz_ogrid-nghost_read_fewer,:)
!
!  The following 3 possibilities allow us to replicate 1-D data input
!  in x (nghost_read_fewer=-1), y (-2), or z (-3) correspondingly.
!
        elseif (nghost_read_fewer==-1) then
          read (lun_input) a(:,1:1+nghost*2,1:1+nghost*2,:)
          a=spread(spread(a(:,m1_ogrid,n1_ogrid,:),2,my_ogrid),3,mz_ogrid)
        elseif (nghost_read_fewer==-2) then
          read (lun_input) a(1:1+nghost*2,:,1:1+nghost*2,:)
          a=spread(spread(a(l1_ogrid,:,n1_ogrid,:),1,mx_ogrid),3,mz_ogrid)
        elseif (nghost_read_fewer==-3) then
          read (lun_input) a(1:1+nghost*2,1:1+nghost*2,:,:)
          a=spread(spread(a(l1_ogrid,m1_ogrid,:,:),1,mx_ogrid),2,my_ogrid)
        else
          call fatal_error('read_snap_single_ogrid','nghost_read_fewer must be >=0')
        endif
      endif

      if (mode == 1) then
!
!  Check whether we want to read deltay from snapshot.
!
        if (lshear) then
          read (lun_input) t_sp, x_ogrid, y_ogrid, z_ogrid, dx_ogrid, dy_ogrid, dz_ogrid, deltay_ogrid
        else
          if (nghost_read_fewer==0) then
            read (lun_input) t_sp, x_ogrid, y_ogrid, z_ogrid, dx_ogrid, dy_ogrid, dz_ogrid
          elseif (nghost_read_fewer>0) then
            read (lun_input) t_sp
          endif
        endif
!
!  Verify consistency of the snapshots regarding their timestamp,
!  unless ireset_tstart=T, in which case we reset all times to tstart.
!
        if ((ireset_tstart == 0) .or. (tstart == impossible)) then
!
          t_test = t_sp
          call mpibcast_real(t_test,comm=MPI_COMM_WORLD)
          call mpiallreduce_or((t_test /= t_sp) .and. .not. lread_from_other_prec &
                                .or. (abs(t_test-t_sp) > 1.e-6),ltest,MPI_COMM_WORLD)
!
!  If timestamps deviate at any processor
!
          if (ltest) then
            if (ireset_tstart > 0) then
!
!  If reset of tstart enabled and tstart unspecified, use minimum of all t_sp
!
              if (ireset_tstart == MINT) then
                call mpiallreduce_min(t_sp,t_sgl,MPI_COMM_WORLD)
                if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT.',&
                                       ' Using (min) t=', t_sgl,'with ireset_tstart=', MINT,'.'
              elseif (ireset_tstart >= MAXT) then
                call mpiallreduce_max(t_sp,t_sgl,MPI_COMM_WORLD)
                if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT.',&
                                       ' Using (max) t=', t_sgl,'with ireset_tstart=', MAXT,'.'
              endif
              tstart=t_sgl
            else
              write (*,*) 'ERROR: '//trim(directory_snap)//'/'//trim(file)// &
                          ' IS INCONSISTENT: t=', t_sp
              call stop_it('read_snap_single')
            endif
          else
            tstart=t_sp
          endif
!
!  Setting time is done in main snap reading rountine
!  Check that time read from overlapping grids match
!  
          if(t_sp/=t) then
            call fatal_error ('read_snap_single_ogrid', 'time differs for cylindrical and cartesian snapshot')
          endif
        endif
      endif
    endsubroutine read_snap_single_ogrid
!***********************************************************************
    subroutine read_snap_double_ogrid(file,a,x,y,z,dx,dy,dz,deltay,nv,mode)
!
!  Read snapshot file in double precision, possibly with mesh and time (if mode=1).
!
!  21-feb-17/Jorgen: Adapted from io_dist.f90
!                             
      use Mpicomm, only: start_serialize, mpibcast_real, mpiallreduce_or, &
                         stop_it, mpiallreduce_min, mpiallreduce_max, MPI_COMM_WORLD
      use IO, only: lun_input
!
      character (len=*), intent(in) :: file
      integer, intent(in) :: nv, mode
      real :: deltay_ogrid=0.0
      real(KIND=rkind8), dimension (mx_ogrid,my_ogrid,mz_ogrid,nv), intent(out) :: a
!
      real(KIND=rkind8) :: t_sp, t_dbl

      real(KIND=rkind8), intent(out) :: dx, dy, dz, deltay
      real(KIND=rkind8), dimension (mx_ogrid), intent(out) :: x
      real(KIND=rkind8), dimension (my_ogrid), intent(out) :: y
      real(KIND=rkind8), dimension (mz_ogrid), intent(out) :: z

      real :: t_test   ! t in single precision for backwards compatibility
      logical :: ltest
      ! set ireset_tstart to 1 or 2 to coordinate divergent timestamp
      integer :: MINT=1
      integer :: MAXT=2
!
      if (lserial_io) call start_serialize
      open (lun_input, FILE=trim(directory_snap)//'/'//file, FORM='unformatted', status='old')
      if (lwrite_2d) then
        if (nz_ogrid == 1) then
          read (lun_input) a(:,:,4,:)
        else
          call fatal_error ('read_snap_double_ogrid','lwrite_2d used for simulation with nz_ogri/=1!')
        endif
      else
!
!  Possibility of reading data with different numbers of ghost zones.
!  In that case, one must regenerate the mesh with luse_oldgrid=T.
!
        if (nghost_read_fewer==0) then
          read (lun_input) a
        elseif (nghost_read_fewer>0) then
          read (lun_input) &
              a(1+nghost_read_fewer:mx_ogrid-nghost_read_fewer, &
                1+nghost_read_fewer:my_ogrid-nghost_read_fewer, &
                1+nghost_read_fewer:mz_ogrid-nghost_read_fewer,:)
!
!  The following 3 possibilities allow us to replicate 1-D data input
!  in x (nghost_read_fewer=-1), y (-2), or z (-3) correspondingly.
!
        elseif (nghost_read_fewer==-1) then
          read (lun_input) a(:,1:1+nghost*2,1:1+nghost*2,:)
          a=spread(spread(a(:,m1_ogrid,n1_ogrid,:),2,my_ogrid),3,mz_ogrid)
        elseif (nghost_read_fewer==-2) then
          read (lun_input) a(1:1+nghost*2,:,1:1+nghost*2,:)
          a=spread(spread(a(l1_ogrid,:,n1_ogrid,:),1,mx_ogrid),3,mz_ogrid)
        elseif (nghost_read_fewer==-3) then
          read (lun_input) a(1:1+nghost*2,1:1+nghost*2,:,:)
          a=spread(spread(a(l1_ogrid,m1_ogrid,:,:),1,mx_ogrid),2,my_ogrid)
        else
          call fatal_error('read_snap_double','nghost_read_fewer must be >=0')
        endif
      endif

      if (ip <= 8) print *, 'read_snap: read ', file
      if (mode == 1) then
!
!  Check whether we want to read deltay from snapshot.
!
        if (lshear) then
          read (lun_input) t_sp, x, y, z, dx, dy, dz, deltay
        else
          if (nghost_read_fewer==0) then
            read (lun_input) t_sp, x, y, z, dx, dy, dz
          elseif (nghost_read_fewer>0) then
            read (lun_input) t_sp
          endif
        endif
!
!  Verify consistency of the snapshots regarding their timestamp,
!  unless ireset_tstart=T, in which case we reset all times to tstart.
!
        if ((ireset_tstart == 0) .or. (tstart == impossible)) then
!
          t_test = t_sp
          call mpibcast_real(t_test,comm=MPI_COMM_WORLD)
          call mpiallreduce_or((t_test /= t_sp) .and. .not. lread_from_other_prec &
                               .or. (abs(t_test-t_sp) > 1.e-6),ltest, MPI_COMM_WORLD)
!
!  If timestamp deviates at any processor
!
          if (ltest) then
            if (ireset_tstart > 0) then
!
!  If reset of tstart enabled and tstart unspecified, use minimum of all t_sp
!
              if (ireset_tstart == MINT) then
                call mpiallreduce_min(t_sp,t_dbl,MPI_COMM_WORLD)
                if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT.',&
                                       ' Using (min) t=', t_dbl,'with ireset_tstart=', MINT,'.'
              elseif (ireset_tstart >= MAXT) then
                call mpiallreduce_max(t_sp,t_dbl,MPI_COMM_WORLD)
                if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT.',&
                                       ' Using (max) t=', t_dbl,'with ireset_tstart=', MAXT,'.'
              endif
              tstart=t_dbl
              if (lroot) write (*,*) 'Timestamps in snapshot INCONSISTENT. Using t=', tstart, '.'
            else
              write (*,*) 'ERROR: '//trim(directory_snap)//'/'//trim(file)// &
                          ' IS INCONSISTENT: t=', t_sp
              call stop_it('read_snap_double')
            endif
          else
            tstart=t_sp
          endif
!
!  Setting time is done in main snap reading rountine
!  Check that time read from overlapping grids match
!  
          if(t_sp/=t) then
            call fatal_error ('read_snap_double_ogrid', 'time differs for cylindrical and cartesian snapshot')
          endif
        endif
      endif
!
    endsubroutine read_snap_double_ogrid
!***********************************************************************
    subroutine setup_mm_nn_ogrid
!
!  Produce index-array for the sequence of points to be worked through:
!  Before the communication has been completed, the nghost=3 layers next
!  to the processor boundary (m1, m2, n1, or n2) cannot be used yet.
!  In the mean time we can calculate the interior points sufficiently far
!  away from the boundary points. Here we calculate the order in which
!  m and n are executed. At one point, necessary(imn)=.true., which is
!  the moment when all communication must be completed.
!
!  24-feb-17/Jorgen: Adapted from general.f90
!
      integer :: min_m1i_m2,max_m2i_m1
      integer :: n1i_ogrid=n1_ogrid+nghost-1
      integer :: n2i_ogrid=mz_ogrid-2*nghost+1
      integer :: m1i_ogrid=m1_ogrid+nghost-1
      integer :: m2i_ogrid=my_ogrid-2*nghost+1
!
!  For non-parallel runs simply go through m and n from bottom left to to right.
!
      imn_array_ogrid=0
      if (ncpus==1) then
        imn_ogrid=1
        necessary_ogrid(1)=.true.
        do n_ogrid=n1_ogrid,n2_ogrid
          do m_ogrid=m1_ogrid,m2_ogrid
            mm_ogrid(imn_ogrid)=m_ogrid
            nn_ogrid(imn_ogrid)=n_ogrid
            imn_array_ogrid(m_ogrid,n_ogrid)=imn_ogrid
            imn_ogrid=imn_ogrid+1
          enddo
        enddo
      else
        imn_ogrid=1
        do n_ogrid=n1i_ogrid+2,n2i_ogrid-2
          do m_ogrid=m1i_ogrid+2,m2i_ogrid-2
            if (imn_array_ogrid(m_ogrid,n_ogrid) == 0) then
              mm_ogrid(imn_ogrid)=m_ogrid
              nn_ogrid(imn_ogrid)=n_ogrid
              imn_array_ogrid(m_ogrid,n_ogrid)=imn_ogrid
              imn_ogrid=imn_ogrid+1
            endif
          enddo
        enddo
        necessary_ogrid(imn_ogrid)=.true.
!
!  Do the upper stripe in the n-direction.
!
        do n_ogrid=max(n2i_ogrid-1,n1_ogrid+1),n2_ogrid
          do m_ogrid=m1i_ogrid+2,m2i_ogrid-2
            if (imn_array_ogrid(m_ogrid,n_ogrid) == 0) then
              mm_ogrid(imn_ogrid)=m_ogrid
              nn_ogrid(imn_ogrid)=n_ogrid
              imn_array_ogrid(m_ogrid,n_ogrid)=imn_ogrid
              imn_ogrid=imn_ogrid+1
            endif
          enddo
        enddo
!
!  Do the lower stripe in the n-direction.
!
        do n_ogrid=n1_ogrid,min(n1i_ogrid+1,n2_ogrid)
          do m_ogrid=m1i_ogrid+2,m2i_ogrid-2
            if (imn_array_ogrid(m_ogrid,n_ogrid) == 0) then
              mm_ogrid(imn_ogrid)=m_ogrid
              nn_ogrid(imn_ogrid)=n_ogrid
              imn_array_ogrid(m_ogrid,n_ogrid)=imn_ogrid
              imn_ogrid=imn_ogrid+1
            endif
          enddo
        enddo
!
!  Left and right hand boxes.
!  NOTE: need to have min(m1i,m2) instead of just m1i, and max(m2i,m1)
!  instead of just m2i, to make sure the case ny=1 works ok, and
!  also that the same m is not set in both loops.
!  ALSO: need to make sure the second loop starts not before the
!  first one ends; therefore max_m2i_m1+1=max(m2i,min_m1i_m2+1).
!
        min_m1i_m2=min(m1i_ogrid+1,m2_ogrid)
        max_m2i_m1=max(m2i_ogrid-1,min_m1i_m2+1)
!
        do n_ogrid=n1_ogrid,n2_ogrid
          do m_ogrid=m1_ogrid,min_m1i_m2
            if (imn_array_ogrid(m_ogrid,n_ogrid) == 0) then
              mm_ogrid(imn_ogrid)=m_ogrid  
              nn_ogrid(imn_ogrid)=n_ogrid
              imn_array_ogrid(m_ogrid,n_ogrid)=imn_ogrid
              imn_ogrid=imn_ogrid+1
            endif
          enddo
          do m_ogrid=max_m2i_m1,m2_ogrid
            if (imn_array_ogrid(m_ogrid,n_ogrid) == 0) then
              mm_ogrid(imn_ogrid)=m_ogrid
              nn_ogrid(imn_ogrid)=n_ogrid
              imn_array_ogrid(m_ogrid,n_ogrid)=imn_ogrid
              imn_ogrid=imn_ogrid+1
            endif
          enddo
        enddo
      endif
!
    endsubroutine setup_mm_nn_ogrid
!***********************************************************************
    subroutine der_ogrid_SBP(f,k,df)
! 
!  Summation by parts boundary condition for first derivative.
!  Only implemented in radial direction.
!
!  21-mar-17/Jorgen: Coded
      real, dimension(l1_ogrid+8,my_ogrid,mz_ogrid,mfarray), intent(in) :: f
      real, dimension(6), intent(out) :: df
      integer, intent(in) :: k
      integer :: i

      do i=1,6
        df(i)=dx_1_ogrid(i)*(D1_SBP(i,1)*f(l1_ogrid  ,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,2)*f(l1_ogrid+1,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,3)*f(l1_ogrid+2,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,4)*f(l1_ogrid+3,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,5)*f(l1_ogrid+4,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,6)*f(l1_ogrid+5,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,7)*f(l1_ogrid+6,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,8)*f(l1_ogrid+7,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,9)*f(l1_ogrid+8,m_ogrid,n_ogrid,k) )
      enddo

    endsubroutine der_ogrid_SBP
!***********************************************************************
    subroutine der2_ogrid_SBP(f,k,df2)
! 
!  Summation by parts boundary condition for second derivative.
!  Only implemented in radial direction.
!
!  21-mar-17/Jorgen: Coded
      real, dimension(l1_ogrid+8,my_ogrid,mz_ogrid,mfarray), intent(in) :: f
      real, dimension(6), intent(out) :: df2
      integer, intent(in) :: k
      integer :: i

      do i=1,6
        df2(i)=(dx_1_ogrid(i)**2)*(D2_SBP(i,1)*f(l1_ogrid  ,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,2)*f(l1_ogrid+1,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,3)*f(l1_ogrid+2,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,4)*f(l1_ogrid+3,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,5)*f(l1_ogrid+4,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,6)*f(l1_ogrid+5,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,7)*f(l1_ogrid+6,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,8)*f(l1_ogrid+7,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,9)*f(l1_ogrid+8,m_ogrid,n_ogrid,k) )
      enddo

    endsubroutine der2_ogrid_SBP
!*********************************************************************** 
    subroutine der_ogrid_SBP_experimental(f,k,df)
! 
!  Summation by parts boundary condition for first derivative.
!  Only implemented in radial direction.
!  This experimental routine is for the outer boundary
!
!  12-may-17/Jorgen: Coded
!
      real, dimension(9,my_ogrid,mz_ogrid,mfarray), intent(in) :: f
      real, dimension(6), intent(out) :: df
      integer, intent(in) :: k
      real :: fac
      integer :: i

      do i=1,6
        fac = dx_1_ogrid(mx_ogrid-i+1)
        df(7-i)=fac*        (D1_SBP(i,1)*f(9,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,2)*f(8,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,3)*f(7,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,4)*f(6,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,5)*f(5,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,6)*f(4,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,7)*f(3,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,8)*f(2,m_ogrid,n_ogrid,k) + &
                             D1_SBP(i,9)*f(1,m_ogrid,n_ogrid,k) )
      enddo

    endsubroutine der_ogrid_SBP_experimental
!***********************************************************************
    subroutine der2_ogrid_SBP_experimental(f,k,df2)
! 
!  Summation by parts boundary condition for second derivative.
!  Only implemented in radial direction.
!  This experimental routine is for the outer boundary
!
!  12-may-17/Jorgen: Coded
      real, dimension(9,my_ogrid,mz_ogrid,mfarray), intent(in) :: f
      real, dimension(6), intent(out) :: df2
      integer, intent(in) :: k
      real :: fac
      integer :: i

      do i=1,6
        fac = dx_1_ogrid(mx_ogrid-i+1)
        df2(7-i)=(fac**2)       *(D2_SBP(i,1)*f(9,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,2)*f(8,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,3)*f(7,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,4)*f(6,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,5)*f(5,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,6)*f(4,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,7)*f(3,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,8)*f(2,m_ogrid,n_ogrid,k) + &
                                  D2_SBP(i,9)*f(1,m_ogrid,n_ogrid,k) )
      enddo

    endsubroutine der2_ogrid_SBP_experimental
!*********************************************************************** 
end module Solid_Cells
