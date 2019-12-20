module solid_cells_ogrid_cdata
!
  use Cparam
!
  implicit none
!
  public
!
!  Force same timestep on all grids
  logical :: lock_dt=.false.
!  Cylinder parameters
  logical :: lset_flow_dir=.false.
  real :: cylinder_radius=0.                        ! Set in start.in
  real :: cylinder_temp=703.0                       ! Set in start.in
  real :: cylinder_xpos=0., cylinder_ypos=0.        ! Set in start.in
  real :: cylinder_zpos=0.                          ! Set in start.in
  real :: skin_depth_solid=0.                       ! Set in start.in
  real :: init_uu=0., ampl_noise=0.
  character(len=labellen) :: initsolid_cells='cylinderstream_x'! Set in start.in
  real :: T0 ! Inlet temperature
  integer :: ncylinders=1,flow_dir=0, flow_dir_set
  real, dimension(3) :: xyz0_ogrid, Lxyz_ogrid, xorigo_ogrid
!  Shift periodic grid, useful for accuracy assessment
  logical, dimension(3) :: lshift_origin_ogrid=.false. ! Set in start.in
  logical, dimension(3) :: lshift_origin_lower_ogrid=.false.   ! Set in start.in
!  Boundary condition
  logical :: SBP=.true.
  logical :: BDRY5=.false.
  logical :: SBP_optimized=.false.
!  Filtering
  logical :: lfilter_solution=.false.
  logical :: lfilter_rhoonly=.false.
  logical :: lfilter_TT=.false.
  real, dimension(:,:,:,:), allocatable ::  f_filterH_lowerx,f_filterH_upperx
  real, dimension(:,:,:,:), allocatable ::  f_filterH_lowery,f_filterH_uppery
  integer :: filter_frequency = 1
!  Free paramter in filter coefficients:
  real :: af=0.1
  integer, parameter :: filter_Hsize = 10/2-nghost
!  Particle interpolation scheme
  integer :: particle_interpolate=1     ! 1: linear, 2: pseudo_quadratic, 3: quadratic
  logical :: lparticle_uradonly=.false.  ! Turn on to have linear inteprolation of all but ur
  logical :: lspecial_rad_int=.false.  ! Turn on to have special quadratic interpolation for particles near the surface
  logical :: lspecial_rad_int_mom=.false.  ! Turn on to have special quadratic interpolation for particles near within momentum thickness
  real :: delta_momentum    ! Momentum thickness of boundary layer
!  Time discretization
  logical :: lrk_tvd=.false.
!  Interpolation method
  integer :: interpolation_method=1  ! Set in start.in
  integer :: inter_len=2             ! Length of interpolation stencil
  integer :: interp_shift=0       ! Distance between interpolation r_int_outer and r_ogrid
  integer :: interpol_filter=0       ! Shift by x cells in filtering
  integer :: interpol_order_poly=5   ! Interpolation order for polynomial interpolation
!  Fundamental grid parameters
  real :: r_ogrid=0.                                                 ! Set in start.in?
  character (len=labellen), dimension(3) :: grid_func_ogrid='linear' ! Set in start.in
  integer, parameter :: nx_ogrid=nxgrid_ogrid/nprocx,ny_ogrid=nygrid_ogrid/nprocy,nz_ogrid=nzgrid_ogrid/nprocz
!***************************************************
! Pencil case ogrid
  integer, parameter :: npencils_ogrid=55
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
    real, dimension (nx_ogrid)     :: ugTT
    real, dimension (nx_ogrid)     :: TT
    real, dimension (nx_ogrid,3)   :: gTT
    real, dimension (nx_ogrid)     :: del2TT
    real, dimension (nx_ogrid)     :: lambda
    real, dimension (nx_ogrid,3)   :: glambda
    real, dimension (nx_ogrid,nchemspec)     :: Diff_penc_add
    real, dimension (nx_ogrid,nchemspec)     :: DYDt_diff
    real, dimension (nx_ogrid,nchemspec)     :: DYDt_reac
    real, dimension (nx_ogrid,nchemspec)     :: H0_RT
    real, dimension (nx_ogrid,nchemspec)     :: hhk_full
    real, dimension (nx_ogrid,3,nchemspec)   :: ghhk
    real, dimension (nx_ogrid,nchemspec)     :: S0_R
    real, dimension (nx_ogrid,3)   :: glncp
    real, dimension (nx_ogrid)     :: cp
    real, dimension (nx_ogrid)     :: cv
    real, dimension (nx_ogrid)     :: nu
    real, dimension (nx_ogrid,3)   :: gradnu
    real, dimension (nx_ogrid)      :: cv1          
    real, dimension (nx_ogrid)      :: cp1           
    real, dimension (nx_ogrid,3)    :: glnTT         
    real, dimension (nx_ogrid)      :: del2lnTT      
    real, dimension (nx_ogrid,3)    :: rho1gpp    
    real, dimension (nx_ogrid)      :: TT1
    real, dimension (nx_ogrid)      :: RR
    real, dimension (nx_ogrid,3)    :: glnRR  
  endtype pencil_case_ogrid
!  
  integer :: i_og_x_mn    =1
  integer :: i_og_y_mn    =2
  integer :: i_og_z_mn    =3
  integer :: i_og_rcyl_mn =4
  integer :: i_og_phi_mn  =5
  integer :: i_og_rcyl_mn1=6
  integer :: i_og_fpres   =7  
  !integer :: i_og_fvisc   =8 ! Unused variables
  !integer :: i_og_rho     =9
  integer :: i_og_rho1    =10
  integer :: i_og_lnrho   =11
  !integer :: i_og_grho    =12
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
  integer :: i_og_ugTT    =30
  integer :: i_og_TT      =31
  integer :: i_og_gTT     =32
  integer :: i_og_del2TT  =33
  integer :: i_og_lambda  =34
  integer :: i_og_Diff_penc_add =35
  integer :: i_og_DYDt_diff     =36
  integer :: i_og_DYDt_reac     =37
  integer :: i_og_H0_RT         =38
  integer :: i_og_hhk_full      =39
  integer :: i_og_ghhk          =40
  integer :: i_og_glncp         =41
  integer :: i_og_cv            =42
  integer :: i_og_cp            =43
  integer :: i_og_nu            =44
  integer :: i_og_gradnu        =45
  integer :: i_og_glambda       =46
  integer :: i_og_cv1           =47
  integer :: i_og_cp1           =48
  integer :: i_og_glnTT         =49
  integer :: i_og_del2lnTT      =50
  integer :: i_og_glnRR         =51
  integer :: i_og_rho1gpp       =52
  integer :: i_og_RR            =53
  integer :: i_og_TT1           =54
  integer :: i_og_S0_R          =55
!
  ! Unused variables
  !character (len=15), parameter, dimension(npencils_ogrid) :: pencil_names_ogrid = &
  !  (/ 'x_mn          ', 'y_mn          ', 'z_mn          ', 'rcyl_mn       '  &
  !   , 'phi_mn        ', 'rcyl_mn1      ', 'fpres         ', 'fvisc         '  &
  !   , 'rho           '  &
  !   , 'rho1          ', 'lnrho         ', 'grho          ', 'glnrho        '  &
  !   , 'ugrho         ', 'sglnrho       '  &
  !   , 'uu            ', 'u2            ', 'uij           ', 'divu          '  &
  !   , 'sij           ', 'sij2          ', 'ugu           ', 'ugu2          '  &
  !   , 'del2u         ', 'graddivu      '  &
  !   , 'lnTT          '  &
  !   , 'cs2           '  &
  !   , 'pp            ', 'ss            ' /)
  logical,dimension(npencils_ogrid):: lpencil_ogrid
!***************************************************
! DATA TYPE FOR INTERPOLATION AND INTERPOLATION STENCILS
!***************************************************
  type interpol_comm_metadata
    integer :: send_to_proc
    integer :: ip_id
    integer, dimension(3) :: i_near_neighbour         ! Lower left corner for interpolation_method=1, closest point for > 1
    integer, dimension(3) :: i_global_neighbour       
    real, dimension(3) :: xyz
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
  real :: r_int_outer, r_int_inner, r_int_inner_poly, r_int_inner_vid=0.
  logical :: lbidiagonal_derij_ogrid=.false.
  integer :: interpol_max                             ! # of interpolated points when not writing video
  logical, pointer :: linterp_pressure                 ! interpolate pressure, recover TT from eos
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
  real, dimension(mxgrid_ogrid) :: xglobal_ogrid, dx1global_ogrid
  real, dimension(mygrid_ogrid) :: yglobal_ogrid, dy1global_ogrid
  real, dimension(mzgrid_ogrid) :: zglobal_ogrid, dz1global_ogrid
!  Necessary for particle runs
  real, dimension(:,:,:,:,:), allocatable ::  f_ogrid_procs
  integer, dimension(:,:), allocatable :: ip_proc
  integer, dimension(:), allocatable :: ip_proc_pointer
  integer, dimension(:), allocatable :: recv_part_data_from
  integer, dimension(:), allocatable :: send_part_data_to
  integer :: n_procs_recv_part_data 
  integer :: n_procs_send_part_data 
  integer :: ivar1_part=1
  integer :: ivar2_part=mvar
  !public :: xgrid_ogrid, ygrid_ogrid, zgrid_ogrid, xglobal_ogrid, yglobal_ogrid, zglobal_ogrid
  !public :: f_ogrid_procs
  !public :: r_ogrid, xorigo_ogrid
!  Local ogrid and derivatives
  integer, parameter :: mx_ogrid=nx_ogrid+2*nghost,l1_ogrid=1+nghost,l2_ogrid=mx_ogrid-nghost
  integer, parameter :: my_ogrid=ny_ogrid+2*nghost,m1_ogrid=1+nghost,m2_ogrid=my_ogrid-nghost
  integer, parameter :: mz_ogrid=nz_ogrid+2*nghost,n1_ogrid=1+nghost,n2_ogrid=mz_ogrid-nghost
  real, dimension (mx_ogrid) :: x_ogrid,dx_1_ogrid,dx2_ogrid,dx_tilde_ogrid,xprim_ogrid
  real, dimension (my_ogrid) :: y_ogrid,dy_1_ogrid,dy2_ogrid,dy_tilde_ogrid,yprim_ogrid
  real, dimension (mz_ogrid) :: z_ogrid,dz_1_ogrid,dz2_ogrid,dz_tilde_ogrid,zprim_ogrid
!  Grid properties computed in grid routines and used in external routines
  integer :: timestep_factor
  real :: dxmin_ogrid,dxmax_ogrid,drmax_ogrid
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
  integer, parameter :: mfarray_ogrid=mvar+maux
  real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), save ::  f_ogrid=0.
  real, dimension(:,:,:,:), allocatable ::  f_tmp ! Array allocated if lrk_tvd=.true.

!  Summation by parts arrays
  real, dimension(6,9) :: D1_SBP, D2_SBP

!  EOS parameters
  real ::  lnTT0!,rho0, lnrho0

!  Energy parameters
  real, pointer :: chi
  logical, pointer :: ladvection_temperature, lheatc_chiconst, lupw_lnTT
  logical :: TT_square_fit = .false.
  real :: Tgrad_stretch=1.

! Eos_chemistry + chemistry parameters
  integer :: ieosvars=-1, ieosvar1=-1, ieosvar2=-1, ieosvar_count=0
  integer :: ll1_ogrid,ll2_ogrid,mm1_ogrid,mm2_ogrid,nn1_ogrid,nn2_ogrid
  logical, pointer :: lheatc_chemistry, lflame_front_2D, lreac_as_aux
  real, dimension(nchemspec) :: chemspec0
  logical :: lreac_heter=.false.
! Reaction rate array needed for BC
  real, dimension (:,:,:), allocatable :: heter_reaction_rate
  real :: Pr_number1
  real, dimension(:), pointer :: Lewis_coef1
  real :: solid_reactions_intro_time=0.0
  real, pointer :: p_init
  logical :: ldist_CO2=.false., ldist_CO=.false.
  integer, pointer :: ireac

!  Diagnostics for output
  integer :: idiag_c_dragx=0
  integer :: idiag_c_dragy=0
  integer :: idiag_Nusselt=0
  integer :: idiag_mdot_C =0

  ! Index for auxiliary gradient of temperature on ogrid, as
  ! well as additional variables for thermophoresis cases
  real :: init_rho_cyl = 1.0
  logical, pointer :: lthermophoretic_forces
  real, dimension(:,:), allocatable ::  curv_cart_transform
  
!  Zero-gradiant boundary condition for rho
  logical :: lexpl_rho = .true.
!  For equation of state
  logical :: leos_isothermal
  logical :: leos_isentropic
!  Read start.in file
!***********************************************************************
endmodule solid_cells_ogrid_cdata
