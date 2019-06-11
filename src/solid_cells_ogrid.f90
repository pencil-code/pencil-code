!
!  This module add solid (as in no-fluid) cells in the domain.
!  This can be used e.g. in order to simulate a cylinder in a cross flow.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lsolid_cells = .true.
! CPARAM logical, parameter :: lsolid_ogrid = .true.
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
  use Solid_cells_ogrid_sub
  use Solid_cells_ogrid_cdata
  use Solid_cells_ogrid_chemistry
!
  implicit none
!
  include 'solid_cells.h'
!
  interface get_polar_coords 
    module procedure get_polar_coords_2D
    module procedure get_polar_coords_3D
    module procedure get_polar_coords_3D_alt
  endinterface
  
  namelist /solid_cells_init_pars/ &
      cylinder_temp, cylinder_radius, cylinder_xpos, ncylinders, &
      cylinder_ypos, cylinder_zpos, flow_dir_set, skin_depth_solid, &
      initsolid_cells, init_uu, r_ogrid, lset_flow_dir,ampl_noise, &
      grid_func_ogrid, coeff_grid_o, xyz_star_ogrid, &
      lcheck_interpolation, lcheck_init_interpolation, SBP, BDRY5, &
      interpolation_method, lock_dt, lexpl_rho, &
      lshift_origin_ogrid,lshift_origin_lower_ogrid, interpol_filter, &
      lrk_tvd, SBP_optimized, interp_shift, &
      particle_interpolate, lparticle_uradonly, &
      interpol_order_poly, lfilter_solution, af, lspecial_rad_int, &
      lfilter_rhoonly, lspecial_rad_int_mom, ivar1_part,ivar2_part, &
      lstore_ogTT, init_rho_cyl, lfilter_TT, r_int_inner_vid, &
      TT_square_fit, Tgrad_stretch, filter_frequency, lreac_heter, solid_reactions_intro_time

!  Read run.in file
  namelist /solid_cells_run_pars/ &
      flow_dir_set, lset_flow_dir, interpolation_method, lcheck_interpolation, &
      SBP, BDRY5, lrk_tvd, SBP_optimized, lexpl_rho, &
      particle_interpolate, lparticle_uradonly, lfilter_solution, lock_dt, af, &
      lspecial_rad_int, lfilter_rhoonly, lspecial_rad_int_mom, &
      ivar1_part,ivar2_part, solid_reactions_intro_time, &
      lstore_ogTT, filter_frequency
    
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
!---------------------------------
!
!  The two lines below can be included, together with the subroutine send_rcv_all_data
!  for MPI-testing purposes. Exchanges all data between processors. Should not be used
!  for other than testing, due to inefficiency.
!!  real, dimension (mxgrid_ogrid, mygrid_ogrid, mzgrid_ogrid,mfarray_ogrid), save ::  fgrid_ogrid=0.
!!  real, dimension (mxgrid, mygrid, mzgrid,mfarray), save ::  fgrid_cartesian=0.
!
!---------------------------------
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
      use SharedVariables, only: get_shared_variable
      use EquationOfState, only: lpres_grad, Pr_number
  !    use Energy, only: lpres_grad
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer :: i, ndims, k
!
      call keep_compiler_quiet(f)
      if (cylinder_radius <= 0) then
        call fatal_error('initialize_solid_cells_ogrid', &
            'All cylinders must have non-zero radii!')
      endif
      if(r_ogrid <= 0) r_ogrid=3.*cylinder_radius
!
      if (r_int_inner_vid <= cylinder_radius) &
          r_int_inner_vid = cylinder_radius+0.05*cylinder_radius
!
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

      call check_cyl_pos(xorigo_ogrid,xyz0,xyz1)
!
!  Currently only implemented for single cylinder 
!
      if(ALWAYS_FALSE) print*, ncylinders
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
! Initial chemical composition
!
! Chemistry
      if (lchemistry) then
        do k = 1,nchemspec
          if (flow_dir == 1) chemspec0(k) = fbcx(ichemspec(k),1)
          if (flow_dir == -1) chemspec0(k) = fbcx(ichemspec(k),2)
          if (flow_dir == 2) chemspec0(k) = fbcy(ichemspec(k),1)
          if (flow_dir == -2) chemspec0(k) = fbcy(ichemspec(k),2)
          if (flow_dir == 3) chemspec0(k) = fbcz(ichemspec(k),1)
          if (flow_dir == -3) chemspec0(k) = fbcz(ichemspec(k),2)
        enddo
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
!  Factor difference between timestep on background grid and overset grid 
!
      if(lock_dt) then
        timestep_factor = 1
      else
!Jorgen:  Removed +1 from timestep_factor, still not safe if ogrid limited by
!         diffusive timestep
        timestep_factor = ceiling(dxmin/dxmin_ogrid)
        if(timestep_factor < 1)  then
          timestep_factor = 1
        endif
        !timestep_factor = timestep_factor*timestep_factor
      endif
!
!  Set interpolation limits
!
      call set_interpolation_limits
!
!  Inform user
!
      if(lroot) then
        if(.not.lequidist_ogrid(1)) then
          print*, ''
          print*, 'Non-linear grid in radial direction - dx_rcyl, dx_rogrid:', &
              0.5*(xglobal_ogrid(nghost+1)-xglobal_ogrid(nghost-1)), &
              0.5*(xglobal_ogrid(mxgrid_ogrid-nghost+1)-xglobal_ogrid(mxgrid_ogrid-nghost-1))
          print*, 'Theta grid spacing - r_cyl*dy_ogrid,r_int_outer*dy_ogrid,r_ogrid*dy_ogrid',&
              xyz0_ogrid(1)*(yglobal_ogrid(2)-yglobal_ogrid(1)), &
              r_int_outer*(yglobal_ogrid(2)-yglobal_ogrid(1)), &
              r_ogrid*(yglobal_ogrid(2)-yglobal_ogrid(1))
          print*, ''
          print*, 'dtheta/dr_surf', &
                  xyz0_ogrid(1)*(yglobal_ogrid(2)-yglobal_ogrid(1)) /&
                  (0.5*(xglobal_ogrid(nghost+1)-xglobal_ogrid(nghost-1)))
          print*, 'dtheta/dr_rogrid',&
                  r_ogrid*(yglobal_ogrid(2)-yglobal_ogrid(1)) / &
                  (0.5*(xglobal_ogrid(mxgrid_ogrid-nghost+1)-xglobal_ogrid(mxgrid_ogrid-nghost-1)) )
          print*, 'dx/dr_rogrid', dx/ &
                  (0.5*(xglobal_ogrid(mxgrid_ogrid-nghost+1)-xglobal_ogrid(mxgrid_ogrid-nghost-1)) )
          print*, 'dx/dtheta_rogrid', dx/ &
                  (r_ogrid*(yglobal_ogrid(2)-yglobal_ogrid(1))) 
          print*, ''
          print*, 'Cartesian grid spacing - dx, dy, dz:', dx,dy,dz
        else
          print*, 'Radial grid spacing - dx_ogrid:', &
              0.5*(xglobal_ogrid(nghost+1)-xglobal_ogrid(nghost-1))
          print*, 'Theta grid spacing - r_cyl*dy_ogrid,r_int_outer*dy_ogrid,r_ogrid*dy_ogrid',&
              xyz0_ogrid(1)*(yglobal_ogrid(2)-yglobal_ogrid(1)), &
              r_int_outer*(yglobal_ogrid(2)-yglobal_ogrid(1)), &
              r_ogrid*(yglobal_ogrid(2)-yglobal_ogrid(1))
          print*, 'Cartesian grid spacing - dx, dy, dz:', dx,dy,dz
        endif
        print*, 'Timestep factor:', timestep_factor
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
      if(BDRY5) then 
        if(lroot) print*, 'Cylinder boundary condition: Fifth order boundary closures'
        SBP=.false.
      elseif(SBP) then
        if(lroot) print*, 'Cylinder boundary condition: Third order SBP boundary closures'
        if(.not. SBP_optimized) then
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
        else
          D1_SBP(1,1) = -21600./13649.; D1_SBP(1,2) = 81763./40947.;  D1_SBP(1,3) = 131./27298. 
          D1_SBP(1,4) = -9143./13649.;  D1_SBP(1,5) = 20539./81894.;  D1_SBP(1,6) = 0. 
          D1_SBP(1,7) = 0.;             D1_SBP(1,8) = 0.;             D1_SBP(1,9) = 0.
          D1_SBP(2,1) = -81763./180195.;D1_SBP(2,2) = 0.;             D1_SBP(2,3) = 7357./36039.
          D1_SBP(2,4) = 30637./72078.;  D1_SBP(2,5) = -2328./12013.;  D1_SBP(2,6) = 6611./360390.
          D1_SBP(2,7) = 0.;             D1_SBP(2,8) = 0.;             D1_SBP(2,9) = 0.
          D1_SBP(3,1) = -131./54220.;   D1_SBP(3,2) = -7357./16266.;  D1_SBP(3,3) = 0.
          D1_SBP(3,4) = 645./2711.;     D1_SBP(3,5) = 11237./32532.;  D1_SBP(3,6) = -3487./27110.
          D1_SBP(3,7) = 0.;             D1_SBP(3,8) = 0.;             D1_SBP(3,9) = 0.
          D1_SBP(4,1) = 9143./53590.;   D1_SBP(4,2) = -30637./64308.; D1_SBP(4,3) = -645./5359.
          D1_SBP(4,4) = 0.;             D1_SBP(4,5) = 13733./32154.;  D1_SBP(4,6) = -67./4660.
          D1_SBP(4,7) = 72./5359.;      D1_SBP(4,8) = 0.;             D1_SBP(4,9) = 0.
          D1_SBP(5,1) = -20539./236310.;D1_SBP(5,2) = 2328./7877.;    D1_SBP(5,3) =-11237./47262.
          D1_SBP(5,4) = -13733./23631.; D1_SBP(5,5) = 0.;             D1_SBP(5,6) = 89387./118155.
          D1_SBP(5,7) = -1296./7877.;   D1_SBP(5,8) = 144./7877.;     D1_SBP(5,9) = 0.
          D1_SBP(6,1) = 0.;             D1_SBP(6,2) = -6611./262806.; D1_SBP(6,3) = 3487./43801.
          D1_SBP(6,4) = 1541./87602.;   D1_SBP(6,5) = -89387./131403.;D1_SBP(6,6) = 0.
          D1_SBP(6,7) = 32400./43801.;  D1_SBP(6,8) = -6480./43801.;  D1_SBP(6,9) = 720./43801.
          D2_SBP(1,1) =0.3548420602490798e1   ; D2_SBP(3,1) =-0.5393903966319141e-1 ; D2_SBP(5,1) = 0.1623318041994786e-1
          D2_SBP(1,2) =-0.1162385694827807e2  ; D2_SBP(3,2) =0.1153943542621719e1   ; D2_SBP(5,2) =-0.8794616833597996e-1
          D2_SBP(1,3) =0.1480964237069501e2   ; D2_SBP(3,3) =-0.2040716873611299e1  ; D2_SBP(5,3) = 0.103577624811612e0
          D2_SBP(1,4) =-0.8968412049815223e1  ; D2_SBP(3,4) =0.698739734417074e0    ; D2_SBP(5,4) = 0.114967901600216e1
          D2_SBP(1,5) =0.2059642370694317e1   ; D2_SBP(3,5) =0.421429883414006e0    ; D2_SBP(5,5) =-0.2443599523155367e1
          D2_SBP(1,6) =0.3761430517226221e0   ; D2_SBP(3,6) =-0.2262171762222378e0  ; D2_SBP(5,6) = 0.1375113224609842e1
          D2_SBP(1,7) =-0.2015793975095019e0  ; D2_SBP(3,7) =0.5090670369467911e-1  ; D2_SBP(5,7) =-0.1218565837960692e0
          D2_SBP(1,8) =0.5117538641997827e-13 ; D2_SBP(3,8) =-0.4371323842747547e-2 ; D2_SBP(5,8) = 0.8668492495883396e-2
          D2_SBP(1,9) =-0.3386357570016522e-15; D2_SBP(3,9) =0.2245491919975288e-3  ; D2_SBP(5,9) = 0.1307369479706344e-3
          D2_SBP(2,1) =0.857883182233682e0    ; D2_SBP(4,1) =-0.2032638843942139e-1 ; D2_SBP(6,1) =-0.3185308684167192e-2
          D2_SBP(2,2) =-0.1397247220064007e1  ; D2_SBP(4,2) =0.4181668262047738e-1  ; D2_SBP(6,2) = 0.1943844988205038e-1
          D2_SBP(2,3) =0.3461647289468133e-1  ; D2_SBP(4,3) =0.1009041221554696e1   ; D2_SBP(6,3) =-0.3865422059089032e-1
          D2_SBP(2,4) =0.6763679122231971e0   ; D2_SBP(4,4) =-0.2044119911750601e1  ; D2_SBP(6,4) =-0.8123817099768654e-1
          D2_SBP(2,5) =-0.1325900419870384e0  ; D2_SBP(4,5) =0.9609112011420257e0   ; D2_SBP(6,5) = 0.1445296692538394e1
          D2_SBP(2,6) =-0.6345391502339508e-1 ; D2_SBP(4,6) =0.9142374273488277e-1  ; D2_SBP(6,6) =-0.2697689107917306e1
          D2_SBP(2,7) =0.244383001412735e-1   ; D2_SBP(4,7) =-0.4316909959745465e-1 ; D2_SBP(6,7) = 0.1494463382995396e1
          D2_SBP(2,8) =-0.2800316968929196e-4 ; D2_SBP(4,8) =0.4668725019017949e-2  ; D2_SBP(6,8) =-0.1495167135596915e0
          D2_SBP(2,9) =0.1331275129575954e-4  ; D2_SBP(4,9) =-0.2461732836225921e-3 ; D2_SBP(6,9) = 0.110849963339009e-1
        endif
      else
        lbidiagonal_derij_ogrid=.true.
        if(lroot) print*, 'WARNING: No cylinder boundary condition set'
      endif
!
!  Set up necessary units for equation of state
!
      if (.not. lchemistry) then
        call initialize_eos_ogr
      else
        call initialize_eos_chemistry
        Pr_number1 = 1./Pr_number
      endif       
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
!  If filterint is used, initialize additional boundary zones (halos) needed for high order filter
!  Must be done before initialization of mpicomm
!
      if(lfilter_solution) call initialize_pade_filter(f_ogrid)
! 
!  Allocate arrays used for communications across processors internally on the curvilinear grid
!
      call initialize_mpicomm_ogrid(lfilter_solution)
!
!  If particles are used in the simulation, set up 'local' arrays of f_ogrid variables to use
!  when computing particle velocities (etc.) for particles inside r_ogrid
!
      if(lparticles) then
        call initialize_particles_ogrid(ivar1_part,ivar2_part)
      endif

      if(lroot) then
        print*, 'Interpolation zone: r_ogrid, r_int_outer, r_int_inner',r_ogrid,r_int_outer,r_int_inner
      endif
! Check interpolation method
      if(lroot) then
        if(interpolation_method==1) then
          print*, 'interpolation_method==1: Linear interpolation used'
        elseif(interpolation_method==2) then
          print*, 'interpolation_method==2: Quadratic spline interpolation used'
        elseif(interpolation_method==3) then
          print*, 'interpolation_method==3: Quadratic spline interpolation for velocities, T and species'
          print*, '                       : Linear interpolation for density'
        elseif(interpolation_method==5) then
          print*, 'interpolation_method==5: Polynomial interpolation'
          print*, 'WARNING                : ONLY SERIAL AT THE MOMENT!!'
        elseif(mod(interpolation_method,2)==0) then
          print*, 'interpolation_method==',interpolation_method,': ',interpolation_method,&
                  'th order Lagrangian interpolation used'
        endif
        if(lparticles) then
            print*, ''
          if(particle_interpolate==1) then
            print*, 'particle_interpolate==1: Linear particle interpolation'
          elseif(particle_interpolate==4) then
            print*, 'particle_interpolate==4: Qubic particle interpolation'
            print*, 'WARNING                : NOT PROPERLY TESTED YET!'
          else
            if(particle_interpolate==2) then
              print*, 'particle_interpolate==2: Pseudo-quadratic particle interpolation'
            elseif(particle_interpolate==3) then
              print*, 'particle_interpolate==3: Quadratic particle interpolation'
            else
              call fatal_error('initialize_solid_cells','particle interpolation does not exist') 
            endif
            if(lparticle_uradonly) print*, '                       : Only for radial direction'
          endif
          if(lspecial_rad_int_mom) then
            print*, 'Particle interpolation with special handling near surface,'
            print*, 'withing (momentum thickness + radius):', delta_momentum
          endif
        endif
      endif
!
!  Get thermal diffusivity from energy module
!
      if (iTT .ne. 0) then
        if (.not. lchemistry) call get_shared_variable('chi',chi)
        call get_shared_variable('lheatc_chiconst',lheatc_chiconst)
        call get_shared_variable('ladvection_temperature',ladvection_temperature)
        call get_shared_variable('lupw_lnTT',lupw_lnTT)
      else
       !   call fatal_error('initialize_solid_cells',&
       !       'Must use linear temperature for solid_cells_ogrid') 
      endif
!
      if (lchemistry) then
        call get_shared_variable('lheatc_chemistry',lheatc_chemistry)
        call get_shared_variable('lflame_front_2D',lflame_front_2D)
        call get_shared_variable('linterp_pressure',linterp_pressure)
        call initialize_chemistry_og(f_ogrid)
        call get_shared_variable('Lewis_coef1',Lewis_coef1)
        call get_shared_variable('p_init',p_init)
      endif
!
!  If TVD Runge-Kutta method is used, temoporary array is needed for storage
!
      if(lrk_tvd) allocate(f_tmp(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid))
!
!  For particles with thermophoretic effects, the gradient of the temperature has
!  to be saved in as an auxiliary so that it can then be interpolated over
!  at the particles position      
!
      if (lstore_ogTT) then
         iogTTx = iTT + 1
         iogTTy = iogTTx + 1
         iogTTz = iogTTy + 1
! NILS: This test should be made more generic!!!!!!
         if (mogaux .ne. 3 .and. (.not. lpres_grad)) then
            call fatal_error('solid_cells_ogrid','mogaux .ne. ndims. Mogaux increased?')
         endif
         allocate(curv_cart_transform(my_ogrid,2))
         call create_curv_cart_transform(curv_cart_transform)
      endif   
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
!  curvilinear grids. Hence, the f-array, which is the array of flow variables
!  on the carteisian grid, must use potential flow to set the ghost zones that
!  will be set by intepolation for t>0.
! 
!  Only implemented for a single solid at present!
!
!  XX-feb-17/Jorgen: Coded
!
      use Initcond, only: gaunoise
      use IO,       only: wdim
      use Sub,      only: control_file_exists
      use EquationOfState, only: imass, species_constants
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz) :: mu1_full=0
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid) :: mu1_full_ogr=0
      real :: a2, rr2, rr2_low, rr2_high 
      real :: wall_smoothing,wall_smoothing_temp
      real :: Lorth,flowx,flowy,shift_flow,shift_orth,flow_r,orth_r
      integer i,j,cyl,iflow,iorth, k, j2, j3
      real :: shift_top,shift_bot,r_k_top,r_k_bot,theta_k_top,theta_k_bot
      real :: ur_k_top ,ur_k_bot ,uth_k_top,uth_k_bot
      logical :: lnoerase=.false.
      real :: coef1, coef2, coef0, r_gradT
      real :: lambda_Suth = 1.5e-5, Suth_const = 200.
      real :: Rgas_unit_sys, Rgas
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
      if (TT_square_fit) then
         r_gradT = Tgrad_stretch*r_ogrid
         coef2 = (cylinder_temp-f(l2,m2,n2,ilnTT))/ &
                 (a2-2*r_gradT*xyz0_ogrid(1)+r_gradT**2)
         coef1 = -2*coef2*r_gradT
         coef0 = f(l2,m2,n2,ilnTT) + coef2*r_gradT**2
      endif
      if (unit_system == 'cgs' .and. lchemistry) then
        Rgas_unit_sys = k_B_cgs/m_u_cgs
        Rgas = Rgas_unit_sys/unit_energy
      endif
      do i = l1,l2
        do j = m1,m2
! Choose correct points depending on flow direction
          flow_r=(x(i)-xorigo_ogrid(1))*flowx+(y(j)-xorigo_ogrid(2))*flowy
          orth_r=(x(i)-xorigo_ogrid(1))*flowy+(y(j)-xorigo_ogrid(2))*flowx
          rr2 = flow_r**2+orth_r**2
          if (rr2 > a2) then
            do cyl = 0,100
              if (cyl == 0) then
                wall_smoothing = 1-exp(-(rr2-a2)/skin_depth_solid**2)
                f(i,j,:,iorth) = f(i,j,:,iorth)-init_uu* &
                  2*flow_r*orth_r*a2/rr2**2*wall_smoothing
                f(i,j,:,iflow) = f(i,j,:,iflow)+init_uu* &
                  (0. - a2/rr2 + 2*orth_r**2*a2/rr2**2)*wall_smoothing
                if ((ilnTT /= 0 .and. .not. lchemistry) .or. (lchemistry .and. .not. lflame_front_2D)) then
                  if (TT_square_fit .and. sqrt(rr2) .le. r_gradT) then
                    f(i,j,:,ilnTT) = coef2*rr2 + coef1*sqrt(rr2) + coef0
                  else
                    wall_smoothing_temp = 1-exp(-(rr2-a2)/(sqrt(a2)*Tgrad_stretch)**2)
                    f(i,j,:,ilnTT) = wall_smoothing_temp*f(i,j,:,ilnTT) &
                      +cylinder_temp*(1-wall_smoothing_temp)
                  endif
                  if (.not. lchemistry) then
                    f(i,j,:,irho) = f(l2,m2,n2,irho) &
                                  * f(l2,m2,n2,ilnTT)/f(i,j,:,ilnTT)
                  endif
                endif
                if (lchemistry .and. .not. lflame_front_2D) then
                    do k = 1,nchemspec
                      f(i,j,:,ichemspec(k)) = chemspec0(k)
                      mu1_full(i,j,:)=mu1_full(i,j,:)+f(i,j,:,ichemspec(k))/species_constants(k,imass)
                    enddo
                    f(i,j,:,iRR) = mu1_full(i,j,:)*Rgas
                    f(i,j,:,irho) = p_init/f(i,j,:,iRR)/f(i,j,:,ilnTT)
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
            if ((ilnTT /= 0 .and. .not. lchemistry) .or. (lchemistry .and. .not. lflame_front_2D)) then
              f(i,j,:,ilnTT) = cylinder_temp
              f(i,j,:,irho) = f(l2,m2,n2,irho) &
                *f(l2,m2,n2,ilnTT)/cylinder_temp
            endif
            if (lchemistry .and. .not. lflame_front_2D) then
              do k = 1,nchemspec
                f(i,j,:,ichemspec(k)) = chemspec0(k)
              enddo
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
      if (iTT .ne. 0) then
        f_ogrid(:,:,:,iTT)=f(l1,m1,n1,iTT)
      endif
      if (lchemistry) then  
        do k = 1,nchemspec
          f_ogrid(:,:,:,ichemspec(k)) = f(l1,m1,n1,ichemspec(k))    
        enddo
      endif
      if(ldensity_nolog) then
        f_ogrid(:,:,:,irho)=init_rho_cyl
      else
        call fatal_error('init_solid_cells','Must use linear density for solid_cells_ogrid')
      endif
      call gaunoise_ogrid(ampl_noise,iux,iuz)
      do i=l1_ogrid,l2_ogrid+nghost
        rr2=x_ogrid(i)**2
        wall_smoothing = 1-exp(-(rr2-a2)/skin_depth_solid**2)
        do j=m1_ogrid,m2_ogrid
!  Compute potential flow past single cylinder
          f_ogrid(i,j,:,iux) = +init_uu*(1-a2/rr2)*cos(y_ogrid(j)+flowy)
          f_ogrid(i,j,:,iuy) = -init_uu*(1+a2/rr2)*sin(y_ogrid(j)+flowy)
          if (ilnTT /= 0) then
            if (TT_square_fit) then
              f_ogrid(i,j,:,ilnTT) = coef2*rr2 + coef1*sqrt(rr2) + coef0
            else
              wall_smoothing_temp = 1-exp(-(rr2-a2)/(sqrt(a2)*Tgrad_stretch)**2)
              f_ogrid(i,j,:,ilnTT) = wall_smoothing_temp*f_ogrid(i,j,:,ilnTT) &
                                   + cylinder_temp*(1-wall_smoothing_temp)
            endif
           ! wall_smoothing_temp = 1-exp(-(rr2-a2)/(sqrt(a2))**2)
            f_ogrid(i,j,:,irho) = f_ogrid(l2_ogrid,m2_ogrid,n2_ogrid,irho) &
              *f(l2,m2,n2,ilnTT)/f_ogrid(i,j,:,ilnTT)
! EWA: I changed f_ogrid to f_cartesian here (see line above and below
!      this comment) as it causes my simulations to crash
!      *f_ogrid(l2_ogrid,m2_ogrid,n2_ogrid,ilnTT)/f_ogrid(i,j,:,ilnTT)
!
! TODO: Set initial conditions for chemistry on the ogrid
            if (lchemistry) then  
              do k = 1,nchemspec
      !          if (ichemspec(k) == 7) then
      !            f_ogrid(i,j,:,ichemspec(k)) = (1-wall_smoothing)
      !          else
      !            f_ogrid(i,j,:,ichemspec(k)) = (1-f_ogrid(i,j,:,7))*chemspec0(k)
      !          endif
                f_ogrid(i,j,:,ichemspec(k)) = chemspec0(k)
              enddo
              if (lreac_heter) then
                ! Heterogeneous reactions only work with simplified mechanism
                f_ogrid(i,j,:,ichemspec(4)) = f_ogrid(i,j,:,ichemspec(4))*wall_smoothing
                f_ogrid(i,j,:,ichemspec(2)) = 1.-f_ogrid(i,j,:,ichemspec(5))-f_ogrid(i,j,:,ichemspec(4))&
                                                -f_ogrid(i,j,:,ichemspec(1))-f_ogrid(i,j,:,ichemspec(3))
                ! Set dp/dr = 0 at the cylinder surface when heter. reactions
                lexpl_rho = .true.
              endif
            endif
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
      if (lchemistry .and. lflame_front_2D) then
!
        f_ogrid(:,:,:,irho) = f(l2,m2,n2,irho)
        f_ogrid(:,:,:,iTT) = f(l2,m2,n2,iTT)
        do k = 1,nchemspec
          f_ogrid(:,:,:,ichemspec(k)) = f(l2,m2,n2,ichemspec(k))
        enddo
!
      elseif (lchemistry) then
!
        do k = 1,nchemspec
          mu1_full_ogr(:,:,:) = mu1_full_ogr(:,:,:)+f_ogrid(:,:,:,ichemspec(k))/species_constants(k,imass)
        enddo
!
        f_ogrid(:,:,:,iRR)  = mu1_full_ogr(:,:,:)*Rgas
        f_ogrid(:,:,:,irho) = p_init/f_ogrid(:,:,:,iRR)/f_ogrid(:,:,:,ilnTT)
!
        do j3 = 1,mz_ogrid
          do j2 = 1,my_ogrid
              f_ogrid(:,j2,j3,iviscosity) = lambda_Suth*f_ogrid(:,j2,j3,iTT)**(3./2.)&
                                          /(Suth_const+f_ogrid(:,j2,j3,iTT))/f_ogrid(:,j2,j3,irho)
          enddo
        enddo
      endif
!
!  Write initial condition to disk.
!
      if (lwrite_ic.and.lstart) then
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
!  Write ogdim.dat
!
      call wdim('ogdim.dat', mx_ogrid, my_ogrid, mz_ogrid, mxgrid_ogrid, mygrid_ogrid, mzgrid_ogrid)
!
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
      p_ogrid%ugTT=penc0
      p_ogrid%TT=penc0
      p_ogrid%gTT=penc0
      p_ogrid%del2TT=penc0
      p_ogrid%lambda=penc0
      p_ogrid%glambda=penc0
      p_ogrid%Diff_penc_add=penc0
      p_ogrid%DYDt_diff=penc0
      p_ogrid%DYDt_reac=penc0
      p_ogrid%H0_RT=penc0
      p_ogrid%hhk_full=penc0
      p_ogrid%ghhk=penc0
      p_ogrid%glncp=penc0
      p_ogrid%cv=penc0
      p_ogrid%cp=penc0
      p_ogrid%nu=penc0
      p_ogrid%gradnu=penc0
      p_ogrid%cv1=penc0
      p_ogrid%cp1=penc0
      p_ogrid%glnTT=penc0
      p_ogrid%del2lnTT=penc0
      p_ogrid%glnRR=penc0
      p_ogrid%RR=penc0
      p_ogrid%rho1gpp=penc0
      p_ogrid%TT1=penc0
      p_ogrid%S0_R=penc0
!    
!  Defined which pencils to solve
!
      lpencil_ogrid=.true.

      if (iTT == 0) then
        lpencil_ogrid(i_og_ugTT)=.false.
        lpencil_ogrid(i_og_TT)=.false.
        lpencil_ogrid(i_og_gTT)=.false.
        lpencil_ogrid(i_og_del2TT)=.false.
      endif

      if (nchemspec == 0) then
        lpencil_ogrid(i_og_lambda)=.false.
        lpencil_ogrid(i_og_glambda)=.false.
        lpencil_ogrid(i_og_Diff_penc_add)=.false.
        lpencil_ogrid(i_og_DYDt_diff)=.false.
        lpencil_ogrid(i_og_DYDt_reac)=.false.
        lpencil_ogrid(i_og_H0_RT)=.false.
        lpencil_ogrid(i_og_hhk_full)=.false.
        lpencil_ogrid(i_og_ghhk)=.false.
        lpencil_ogrid(i_og_glncp)=.false.
        lpencil_ogrid(i_og_cv)=.false.
        lpencil_ogrid(i_og_cp)=.false.
        lpencil_ogrid(i_og_nu)=.false.
        lpencil_ogrid(i_og_gradnu)=.false.
        lpencil_ogrid(i_og_cv1)=.false.
        lpencil_ogrid(i_og_cp1)=.false.
        lpencil_ogrid(i_og_glnTT)=.false.
        lpencil_ogrid(i_og_del2lnTT)=.false.
        lpencil_ogrid(i_og_glnRR)=.false.
        lpencil_ogrid(i_og_RR)=.false.
        lpencil_ogrid(i_og_rho1gpp)=.false.
        lpencil_ogrid(i_og_TT1)=.false.
        lpencil_ogrid(i_og_S0_R)=.false.
      endif
!
    endsubroutine initialize_pencils_ogrid
!***********************************************************************
    subroutine initialize_eos_ogr
!  
!  Set up parameters necessary to compute the energy and pressure using
!  the ideal gas eos.
!  This is done in the unit_eos routine in eos_idealgas.f90 for the 
!  cartesian solver
!
!  4-apr-17/Jorgen: Coded
!
      use EquationOfState, only: get_cv1,get_cp1,cs20,gamma_m1,rho0,lnrho0
      real :: cp1, cp
!
!  Inverse cv and cp values.
!
        call get_cp1(cp1)
        cp=1./cp1
!
!        rho0=1.0
!        lnrho0=log(rho0)
        if (gamma_m1/=0.0) then
          lnTT0=log(cs20/(cp*gamma_m1))  !(general case)
          leos_isentropic=.true.
        else
          lnTT0=log(cs20/cp)  !(isothermal/polytropic cases: check!)
          leos_isothermal=.true.
        endif
!
    endsubroutine initialize_eos_ogr
!***********************************************************************
    subroutine initialize_interpolate_points
!
! Build arrays of interpolation data on processors that perform interpolation.
! Necessary to perform communications in an efficient manner.
!
! apr-17/Jorgen: Coded
!
      real, dimension(3) :: xyz,rthz
      integer :: i,j,k,ii
      real, dimension(3) :: xyz_neigh, rthz_neigh
!
!  Set the length of the interpolation stencil based on choice of method 
!
      if(interpolation_method==1) then
        inter_len=2    
      elseif(interpolation_method==2 .or. interpolation_method==3) then
        inter_len=3
   !   elseif(interpolation_method==4) then
   !     inter_len=5
      elseif(interpolation_method==5) then
        inter_len=interpol_order_poly
      elseif(mod(interpolation_method,2)==0) then
        inter_len=interpolation_method+1
      else
        call fatal_error('initialize_interpolate_points','selected interpolation method does not exist!')
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
              call find_near_ind_global_cart(cartesian_to_curvilinear(ii)%ind_global_neighbour, &
                  xyz,lcheck_init_interpolation)
!
!  If higher order interpolation is used, adjust nearest index to be the index to point actually nearest
!  the interpolation point, NOT the index of bottom left point in cell containing interpolation point
!
              if(interpolation_method>1) then
                call adjust_inear_cart_glob(cartesian_to_curvilinear(ii)%ind_global_neighbour,xyz)
              endif
!
              xyz_neigh= (/ xglobal(cartesian_to_curvilinear(ii)%ind_global_neighbour(1)), &
                            yglobal(cartesian_to_curvilinear(ii)%ind_global_neighbour(2)), &
                            zglobal(cartesian_to_curvilinear(ii)%ind_global_neighbour(3)) /)
              call find_proc_cartesian(xyz_neigh,cartesian_to_curvilinear(ii)%from_proc)
              if (cartesian_to_curvilinear(ii)%from_proc == iproc) then
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
            rthz(1)=radius_ogrid(x(i),y(j))
            if((rthz(1)<=r_int_outer) .and. rthz(1)>=r_int_inner_vid) then
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
            call get_polar_coords(x(i),y(j),z(k),rthz)
            if((rthz(1)<=r_int_outer) .and. rthz(1)>=r_int_inner) then
              ii=ii+1
              curvilinear_to_cartesian(ii)%i_xyz = (/ i,j,k /)
              curvilinear_to_cartesian(ii)%xyz = rthz
              call find_near_ind_global_curv(curvilinear_to_cartesian(ii)%ind_global_neighbour, &
                  rthz,lcheck_init_interpolation)
!
!  If higher order interpolation is used, adjust nearest index to be the index to point actually nearest
!  the interpolation point, NOT the index of bottom left point in cell containing interpolation point
!  Also, make sure that no interpolation points try to use points inside the cylinder
!
              if(interpolation_method>1) then
                call adjust_inear_curv_glob(curvilinear_to_cartesian(ii)%ind_global_neighbour,rthz)
              endif
!
              rthz_neigh= (/ xglobal_ogrid(curvilinear_to_cartesian(ii)%ind_global_neighbour(1)), &
                             yglobal_ogrid(curvilinear_to_cartesian(ii)%ind_global_neighbour(2)), &
                             zglobal_ogrid(curvilinear_to_cartesian(ii)%ind_global_neighbour(3)) /)
              call find_proc_curvilinear(rthz_neigh,curvilinear_to_cartesian(ii)%from_proc)
              if(curvilinear_to_cartesian(ii)%from_proc==iproc) then
                call ind_global_to_local_curv(curvilinear_to_cartesian(ii)%ind_global_neighbour, &
                      curvilinear_to_cartesian(ii)%ind_local_neighbour,lcheck_init_interpolation)
              endif
            endif
          enddo
        enddo
      enddo
      interpol_max = ii
      do k=n1,n2
        do j=m1,m2
          do i=l1,l2
            call get_polar_coords(x(i),y(j),z(k),rthz)
            if((rthz(1)<r_int_inner) .and. rthz(1)>=r_int_inner_vid) then
              ii=ii+1
              curvilinear_to_cartesian(ii)%i_xyz = (/ i,j,k /)
              curvilinear_to_cartesian(ii)%xyz = rthz
              call find_near_ind_global_curv(curvilinear_to_cartesian(ii)%ind_global_neighbour, &
                  rthz,lcheck_init_interpolation)
!
!  If higher order interpolation is used, adjust nearest index to be the index to point actually nearest
!  the interpolation point, NOT the index of bottom left point in cell containing interpolation point
!  Also, make sure that no interpolation points try to use points inside the cylinder
!
              if(interpolation_method>1) then
                call adjust_inear_curv_glob(curvilinear_to_cartesian(ii)%ind_global_neighbour,rthz)
              endif
!
              rthz_neigh= (/ xglobal_ogrid(curvilinear_to_cartesian(ii)%ind_global_neighbour(1)), &
                             yglobal_ogrid(curvilinear_to_cartesian(ii)%ind_global_neighbour(2)), &
                             zglobal_ogrid(curvilinear_to_cartesian(ii)%ind_global_neighbour(3)) /)
              call find_proc_curvilinear(rthz_neigh,curvilinear_to_cartesian(ii)%from_proc)
              if(curvilinear_to_cartesian(ii)%from_proc==iproc) then
                call ind_global_to_local_curv(curvilinear_to_cartesian(ii)%ind_global_neighbour, &
                      curvilinear_to_cartesian(ii)%ind_local_neighbour,lcheck_init_interpolation)
              endif
            endif
          enddo
        enddo
      enddo
!print*, 'read ii'
!read(*,*) ii
!
    endsubroutine initialize_interpolate_points
!!***********************************************************************
    subroutine initialize_send_ip_points_alt
!
! Build arrays of interpolation data on processors that contain data 
! necessary for interpolation on other processors. 
!
! apr-17/Jorgen: Coded
!
      use Mpicomm, only: mpirecv_int, mpisend_nonblock_int, mpibarrier, &
                         mpirecv_nonblock_int, mpisend_int, mpiwait, mpibcast_int, &
                         mpirecv_nonblock_real, mpisend_real

      use Solid_Cells_Mpicomm, only: finalize_isend_init_interpol
      integer :: i,iip,npoint
      integer, dimension(ncpus) :: from_proc_curv_to_cart=0
      integer, dimension(ncpus) :: from_proc_cart_to_curv=0
      integer, dimension(:,:,:), allocatable :: ind_from_proc_curv
      integer, dimension(:,:,:), allocatable :: ind_from_proc_cart
      integer, dimension(:,:), allocatable :: ip_id_curv_to_cart
      integer, dimension(:,:), allocatable :: ip_id_cart_to_curv
      integer :: max_from_proc, from_proc
      integer, dimension(ncpus-1) :: ireq1D
      integer, dimension(ncpus-1,3) :: ireq2D
      integer, dimension(ncpus-1,3) :: ireq2D_xyz
      !
      integer, dimension(ncpus,ncpus) :: from_proc_curv_to_cart_glob=0
      integer, dimension(ncpus,ncpus) :: from_proc_cart_to_curv_glob=0
      integer, dimension(ncpus) :: from_proc_bufi =0
      integer :: iter, ind_start, ind_stop, ip_recv_tot, ip_send_tot, n_ip_proc
      integer, dimension(2) :: buf_size
      integer, dimension(:), allocatable :: id_bufi, id_bufo
      integer, dimension(:,:), allocatable :: ijk_bufi, ijk_bufo
      real, dimension(:,:), allocatable :: xyz_bufi, xyz_bufo
      integer, dimension(3) :: indices_global
      !TODO TEMP BELOW
      real, dimension(:,:,:), allocatable :: xyz_from_curv_to_cart


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
        allocate(xyz_from_curv_to_cart(ncpus,max_from_proc,3))
        do iip=0,ncpus-1
          if(from_proc_curv_to_cart(iip+1)>0) then
            npoint=0
            do i=1,n_ip_curv_to_cart
              if(curvilinear_to_cartesian(i)%from_proc==iip) then
                npoint=npoint+1
! Must access iip+1 instead of iip, to avoid accessing element 0
                ind_from_proc_curv(iip+1,npoint,:)=curvilinear_to_cartesian(i)%ind_global_neighbour
                ip_id_curv_to_cart(iip+1,npoint)=i
                xyz_from_curv_to_cart(iip+1,npoint,:)=curvilinear_to_cartesian(i)%xyz
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
            call mpirecv_int(from_proc_bufi,ncpus,iip,110)
            from_proc_curv_to_cart_glob(iip+1,:)=from_proc_bufi
            call mpirecv_int(from_proc_bufi,ncpus,iip,115)
            from_proc_cart_to_curv_glob(iip+1,:)=from_proc_bufi
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
      !max_send_ip_curv_to_cart=maxval(n_ip_to_proc_curv_to_cart)
      !max_send_ip_cart_to_curv=maxval(n_ip_to_proc_cart_to_curv)
      if(n_procs_send_curv_to_cart>0) then
        max_send_ip_curv_to_cart=maxval(n_ip_to_proc_curv_to_cart)
      else
        max_send_ip_curv_to_cart=0
      endif
      if(n_procs_send_cart_to_curv>0) then
        max_send_ip_cart_to_curv=maxval(n_ip_to_proc_cart_to_curv)
      else
        max_send_ip_cart_to_curv=0
      endif
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
      !max_recv_ip_curv_to_cart=maxval(n_ip_recv_proc_curv_to_cart)
      !max_recv_ip_cart_to_curv=maxval(n_ip_recv_proc_cart_to_curv)
      if(n_procs_recv_curv_to_cart>0) then
        max_recv_ip_curv_to_cart=maxval(n_ip_recv_proc_curv_to_cart)
      else
        max_recv_ip_curv_to_cart=0
      endif
      if(n_procs_recv_cart_to_curv>0) then
        max_recv_ip_cart_to_curv=maxval(n_ip_recv_proc_cart_to_curv)
      else
        max_recv_ip_cart_to_curv=0
      endif
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
      allocate(xyz_bufi(ip_send_tot,3))
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
          call mpirecv_nonblock_real(xyz_bufi(ind_start:ind_stop,i),buf_size(1),iip,190+i,ireq2D_xyz(iter,i))
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
      allocate(xyz_bufo(ip_recv_tot,3))
      allocate(id_bufo(ip_recv_tot))
      ind_start=1
      do iter=1,n_procs_recv_curv_to_cart
        n_ip_proc=n_ip_recv_proc_curv_to_cart(iter)
        ind_stop=ind_start+n_ip_proc-1
        iip=procs_recv_curv_to_cart(iter)
        ijk_bufo(ind_start:ind_stop,:)=ind_from_proc_curv(iip+1,1:n_ip_proc,:)
        !TODO TEMP
        xyz_bufo(ind_start:ind_stop,:)=xyz_from_curv_to_cart(iip+1,1:n_ip_proc,:)
        !
        id_bufo(ind_start:ind_stop)=ip_id_curv_to_cart(iip+1,1:n_ip_proc)
        buf_size=(/ind_stop-ind_start+1,3/)
        do i=1,3
          call mpisend_int(ijk_bufo(ind_start:ind_stop,i),buf_size(1),iip,200+i)
          call mpisend_real(xyz_bufo(ind_start:ind_stop,i),buf_size(1),iip,190+i)
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
          call mpiwait(ireq2D_xyz(iter,i))
        enddo
        do i=ind_start,ind_stop
          indices_global = ijk_bufi(i,:)
          send_curvilinear_to_cartesian(i)%i_global_neighbour=indices_global
          send_curvilinear_to_cartesian(i)%xyz=xyz_bufi(i,:)
          call ind_global_to_local_curv(indices_global, &
              send_curvilinear_to_cartesian(i)%i_near_neighbour,lcheck_init_interpolation)
        enddo
        call mpiwait(ireq1D(iter))
        send_curvilinear_to_cartesian(ind_start:ind_stop)%ip_id=id_bufi(ind_start:ind_stop)
        ind_start=ind_stop+1
      enddo
      deallocate(ijk_bufi)
      deallocate(xyz_bufi)
      deallocate(id_bufi)
      deallocate(ijk_bufo)
      deallocate(xyz_bufo)
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
          indices_global=ijk_bufi(i,:)
          call ind_global_to_local_cart(indices_global, &
              send_cartesian_to_curvilinear(i)%i_near_neighbour,lcheck_init_interpolation)
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
!
      i_rthz_local(1) = i_rthz_global(1) - nx_ogrid*ipx
      i_rthz_local(2) = i_rthz_global(2) - ny_ogrid*ipy
      i_rthz_local(3) = i_rthz_global(3) - nz_ogrid*ipz
!
      if(lcheck) then
        if(abs(x_ogrid(i_rthz_local(1))-xglobal_ogrid(i_rthz_global(1)))>1.e-12) &
          print*, 'ERROR: incorrect transformation of global to local coordinates in r-direction'
        if(abs(y_ogrid(i_rthz_local(2))-yglobal_ogrid(i_rthz_global(2)))>1.e-12) &
          print*, 'ERROR: incorrect transformation of global to local coordinates in th-direction'
        if(abs(z_ogrid(i_rthz_local(3))-zglobal_ogrid(i_rthz_global(3)))>1.e-12) &
          print*, 'ERROR: incorrect transformation of global to local coordinates in z-direction'
      endif
!
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
!
      i_xyz_local(1) = i_xyz_global(1) - nx*ipx
      i_xyz_local(2) = i_xyz_global(2) - ny*ipy
      i_xyz_local(3) = i_xyz_global(3) - nz*ipz

      if(lcheck) then
        if(abs(x(i_xyz_local(1))-xglobal(i_xyz_global(1)))>1.e-12) &
          print*, 'ERROR: incorrect transformation of global to local coordinates in r-direction'
        if(abs(y(i_xyz_local(2))-yglobal(i_xyz_global(2)))>1.e-12) &
          print*, 'ERROR: incorrect transformation of global to local coordinates in th-direction'
        if(abs(z(i_xyz_local(3))-zglobal(i_xyz_global(3)))>1.e-12) &
          print*, 'ERROR: incorrect transformation of global to local coordinates in z-direction'
      endif
!
    endsubroutine ind_global_to_local_cart
!!***********************************************************************
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
      use EquationOfState, only: rho0
!
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
    subroutine Nusselt_pencils(Nusselt)
!
!  Compute the Nusselt number of the cylinder 
!
!  23-aug-17/Ewa+Nils: Coded
!
      real, intent(inout) :: Nusselt
      real :: gradT
!
      gradT=p_ogrid%gTT(1,1)
      Nusselt = Nusselt - gradT  
!
    endsubroutine Nusselt_pencils
!***********************************************************************
    subroutine Nusselt_coeffs(Nusselt)
!
!  Sum up the computed Nusselt number on root processor.
!  Normalization done in the end of the computation.
!
!  23-aug-17/Ewa+Nils: Coded
!
      use Mpicomm, only: mpireduce_sum
      real, intent(inout) :: Nusselt
      real :: Nusselt_all
      real :: norm
!
      norm = 2.*cylinder_radius/(cylinder_temp-T0)/(ny_ogrid*nz_ogrid)
!
      call mpireduce_sum(Nusselt,Nusselt_all)
!
      if(lroot) then
        Nusselt=Nusselt_all*norm
        if (idiag_Nusselt /= 0) fname(idiag_Nusselt)=Nusselt
      endif
!
    endsubroutine Nusselt_coeffs
!***********************************************************************
    subroutine dsolid_dt(f,df,p)
!
!  Dummy routine
!
      real, dimension(mx,my,mz,mfarray), intent(in):: f
      real, dimension(mx,my,mz,mvar), intent(in)   :: df
      type (pencil_case), intent(in)               :: p
!
      call keep_compiler_quiet(df,f)
      call keep_compiler_quiet(p)
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
      use General, only: loptest
!
      integer :: iname
      logical :: lreset
      logical, optional :: lwrite
!
!  Dummy variable
!
      if(loptest(lwrite)) print*, lwrite
!
!  Reset everything in case of reset
!
      if (lreset) then
        idiag_c_dragx = 0
        idiag_c_dragy = 0
        idiag_Nusselt = 0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname = 1,nname
        call parse_name(iname,cname(iname),cform(iname),'c_dragx',idiag_c_dragx)
        call parse_name(iname,cname(iname),cform(iname),'c_dragy',idiag_c_dragy)
        call parse_name(iname,cname(iname),cform(iname),'Nusselt',idiag_Nusselt)
      enddo
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
      call keep_compiler_quiet(f)
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
      integer :: i
!
!
      do i=l1,l2
        if(radius_ogrid(x(i),y(m)) <= r_int_outer) then
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
    logical :: in_solid_cell
    real, dimension(:), intent(in) :: part_pos
    real, intent(in) :: part_rad
    real :: r_solid_par
!
    in_solid_cell = .false.
!
    r_solid_par = radius_ogrid(part_pos(1),part_pos(2)) 
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
    use Mpicomm, only: mpisend_int,mpisend_real,mpirecv_int,mpirecv_real
    real, dimension(mx,my,mz,mfarray), intent(in) :: f_cartesian
    integer, intent(in) :: ivar1,ivar2
    integer, dimension(5) :: nbuf_farr
    integer, dimension(max_recv_ip_cart_to_curv) :: id_bufi
    real, dimension(inter_len,inter_len,inter_len,ivar2-ivar1+1,max_send_ip_cart_to_curv) :: f_bufo
    real, dimension(inter_len,inter_len,inter_len,ivar2-ivar1+1,max_recv_ip_cart_to_curv) :: f_bufi
    real, dimension(inter_len,inter_len,inter_len,ivar2-ivar1+1) :: farr
    integer :: i,j,k,id,ipq
    integer :: ii1,ii2,jj1,jj2,kk1,kk2
    integer :: iter, send_to, recv_from
    integer, dimension(3) :: inear_loc
    integer :: ind_send_first, ind_send_last
    integer, dimension(max_send_ip_cart_to_curv) :: ip_bufo
    integer :: ind, ipoly
    integer :: ipp_int=0, iRR_int=0
!
    if(interpolation_method==1) then
      nbuf_farr(1:3)=2
      ii1=0; ii2=1; jj1=0; jj2=1; kk1=0; kk2=1
    elseif(interpolation_method==2 .or. interpolation_method==3) then
      nbuf_farr(1:3)=3
      ii1=1; ii2=1; jj1=1; jj2=1; kk1=1; kk2=1
    elseif(interpolation_method==5) then
      nbuf_farr(1:3)=interpol_order_poly
      ipoly=floor((interpol_order_poly)*0.5)
      ii1=ipoly; ii2=ipoly; jj1=ipoly; jj2=ipoly; kk1=ipoly; kk2=ipoly
    elseif(mod(interpolation_method,2)==0) then
      nbuf_farr(1:3)=interpolation_method+1
      ii1=interpolation_method/2; ii2=ii1; jj1=ii1; jj2=ii1; kk1=ii1; kk2=ii1
    endif
    nbuf_farr(4)=ivar2-ivar1+1
!
    if (lchemistry .and. linterp_pressure) then
      ipp_int  = mvar-ivar1+2
      iRR_int  = mvar-ivar1+3
    endif
!
!  Send to processors with proc > iproc
!
    ind_send_first=1
    do iter=1,n_procs_send_cart_to_curv
      ind_send_last=n_ip_to_proc_cart_to_curv(iter)+ind_send_first-1
      send_to=send_cartesian_to_curvilinear(ind_send_last)%send_to_proc
      if(send_to>iproc) then
        nbuf_farr(5)=ind_send_last-ind_send_first+1
        do ipq=1,nbuf_farr(5)
          ind=ind_send_first+ipq-1
          i=send_cartesian_to_curvilinear(ind)%i_near_neighbour(1)
          j=send_cartesian_to_curvilinear(ind)%i_near_neighbour(2)
          k=send_cartesian_to_curvilinear(ind)%i_near_neighbour(3)
          if (lchemistry .and. linterp_pressure) then
            f_bufo(:,:,:,ipp_int,ipq)=f_cartesian(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,ipp)
            f_bufo(:,:,:,iRR_int,ipq)=f_cartesian(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,iRR)
          endif
          f_bufo(:,:,:,ivar1:mvar,ipq)=f_cartesian(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,ivar1:mvar)
        enddo
        ip_bufo(1:nbuf_farr(5)) = send_cartesian_to_curvilinear(ind_send_first:ind_send_last)%ip_id
        call mpisend_int(ip_bufo(1:nbuf_farr(5)),nbuf_farr(5),send_to,send_to)
        call mpisend_real(f_bufo(:,:,:,:,1:nbuf_farr(5)),nbuf_farr,send_to,send_to+ncpus)
      endif
      ind_send_first=ind_send_last+1
    enddo
!
!  Recieve from processors with proc < iproc
!
    do iter=n_procs_recv_cart_to_curv,1,-1
      recv_from=procs_recv_cart_to_curv(iter)
      if(recv_from<iproc) then
        nbuf_farr(5)=n_ip_recv_proc_cart_to_curv(iter)
        call mpirecv_int(id_bufi(1:nbuf_farr(5)),nbuf_farr(5),recv_from,iproc)
        call mpirecv_real(f_bufi(:,:,:,:,1:nbuf_farr(5)),nbuf_farr,recv_from,iproc+ncpus)
        do ipq=1,nbuf_farr(5)
          farr = f_bufi(:,:,:,:,ipq)
           call interpolate_point_cart_to_curv(id_bufi(ipq),ivar1,ivar2,farr,f_cartesian,ipp_int,iRR_int)
        enddo
      endif
    enddo
!
!  Send to processors with proc < iproc
!
    ind_send_first=1
    do iter=1,n_procs_send_cart_to_curv
      ind_send_last=n_ip_to_proc_cart_to_curv(iter)+ind_send_first-1
      send_to=send_cartesian_to_curvilinear(ind_send_last)%send_to_proc
      if(send_to<iproc) then
        nbuf_farr(5)=ind_send_last-ind_send_first+1
        do ipq=1,nbuf_farr(5)
          ind=ind_send_first+ipq-1
          i=send_cartesian_to_curvilinear(ind)%i_near_neighbour(1)
          j=send_cartesian_to_curvilinear(ind)%i_near_neighbour(2)
          k=send_cartesian_to_curvilinear(ind)%i_near_neighbour(3)
          if (lchemistry .and. linterp_pressure) then
            f_bufo(:,:,:,ipp_int,ipq)=f_cartesian(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,ipp)
            f_bufo(:,:,:,iRR_int,ipq)=f_cartesian(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,iRR)
          endif
          f_bufo(:,:,:,ivar1:mvar,ipq)=f_cartesian(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,ivar1:mvar)
        enddo
        ip_bufo(1:nbuf_farr(5)) = send_cartesian_to_curvilinear(ind_send_first:ind_send_last)%ip_id
        call mpisend_int(ip_bufo(1:nbuf_farr(5)),nbuf_farr(5),send_to,send_to)
        call mpisend_real(f_bufo(:,:,:,:,1:nbuf_farr(5)),nbuf_farr,send_to,send_to+ncpus)
      endif
      ind_send_first=ind_send_last+1
    enddo
!
!  Recieve from processors with proc > iproc
!
    do iter=n_procs_recv_cart_to_curv,1,-1
      recv_from=procs_recv_cart_to_curv(iter)
      if(recv_from>iproc) then
        nbuf_farr(5)=n_ip_recv_proc_cart_to_curv(iter)
        call mpirecv_int(id_bufi(1:nbuf_farr(5)),nbuf_farr(5),recv_from,iproc)
        call mpirecv_real(f_bufi(:,:,:,:,1:nbuf_farr(5)),nbuf_farr,recv_from,iproc+ncpus)
        do ipq=1,nbuf_farr(5)
          farr = f_bufi(:,:,:,:,ipq)
          call interpolate_point_cart_to_curv(id_bufi(ipq),ivar1,ivar2,farr,f_cartesian,ipp_int,iRR_int)
        enddo
      endif
    enddo
!
!  Interpolate remaining points 
!
    do id=1,n_ip_cart_to_curv
      if(cartesian_to_curvilinear(id)%from_proc==iproc) then
        inear_loc=cartesian_to_curvilinear(id)%ind_local_neighbour
          if (lchemistry .and. linterp_pressure) then
            farr(:,:,:,ipp_int)=f_cartesian(inear_loc(1)-ii1:inear_loc(1)+ii2, &
          inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,ipp)
            farr(:,:,:,iRR_int)=f_cartesian(inear_loc(1)-ii1:inear_loc(1)+ii2, &
          inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,iRR)
          endif
          farr(:,:,:,ivar1:mvar)=f_cartesian(inear_loc(1)-ii1:inear_loc(1)+ii2, &
          inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,ivar1:mvar)
        call interpolate_point_cart_to_curv(id,ivar1,ivar2,farr,f_cartesian,ipp_int,iRR_int)
      endif
    enddo
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
    use Mpicomm, only: mpisend_int,mpisend_real,mpirecv_int,mpirecv_real
    use EquationOfState, only: lpres_grad
    real, dimension(mx,my,mz,mfarray), intent(inout) :: f_cartesian
    integer, intent(in) :: ivar1,ivar2
    integer, dimension(5) :: nbuf_farr
    integer, dimension(max_recv_ip_curv_to_cart) :: id_bufi
    real, dimension(inter_len,inter_len,inter_len,ivar2-ivar1+1,max_send_ip_curv_to_cart) :: f_bufo
    real, dimension(inter_len,inter_len,inter_len,ivar2-ivar1+1,max_recv_ip_curv_to_cart) :: f_bufi
    real, dimension(inter_len,inter_len,inter_len,ivar2-ivar1+1) :: farr
    integer :: i,j,k,id,ipq
    integer :: ii1,ii2,jj1,jj2,kk1,kk2
    integer :: iter, send_to, recv_from
    integer, dimension(3) :: inear_loc
    integer :: ind_send_first, ind_send_last 
    integer, dimension(max_send_ip_curv_to_cart) :: ip_bufo
    integer :: ind, ipoly
    integer :: igpx_int=0, igpy_int=0, ipp_int=0, iRR_int=0
!
    if(interpolation_method==1) then
      nbuf_farr(1:3)=2
      ii1=0; ii2=1; jj1=0; jj2=1; kk1=0; kk2=1
    elseif(interpolation_method==2 .or. interpolation_method==3) then
      nbuf_farr(1:3)=3
      ii1=1; ii2=1; jj1=1; jj2=1; kk1=1; kk2=1
    elseif(interpolation_method==5) then
      nbuf_farr(1:3)=interpol_order_poly
      ipoly=floor((interpol_order_poly)*0.5)
      ii1=ipoly; ii2=ipoly; jj1=ipoly; jj2=ipoly; kk1=ipoly; kk2=ipoly
    elseif(mod(interpolation_method,2)==0) then
      nbuf_farr(1:3)=interpolation_method+1
      ii1=interpolation_method/2; ii2=ii1; jj1=ii1; jj2=ii1; kk1=ii1; kk2=ii1
    endif
    nbuf_farr(4)=ivar2-ivar1+1
!
    if (lpres_grad .and. (lchemistry .and. linterp_pressure)) then
      igpx_int = mvar-ivar1+2
      igpy_int = mvar-ivar1+3
      ipp_int  = mvar-ivar1+4
      iRR_int  = mvar-ivar1+5
    elseif (lpres_grad) then
      igpx_int = mvar-ivar1+2
      igpy_int = mvar-ivar1+3
    elseif (lchemistry .and. linterp_pressure) then
      ipp_int  = mvar-ivar1+2
      iRR_int  = mvar-ivar1+3
    endif
!
! Send to proc > iproc
!
    ind_send_first=1
    do iter=1,n_procs_send_curv_to_cart
      ind_send_last=n_ip_to_proc_curv_to_cart(iter)+ind_send_first-1
      send_to=send_curvilinear_to_cartesian(ind_send_last)%send_to_proc
      if(send_to>iproc) then
        nbuf_farr(5)=ind_send_last-ind_send_first+1
        do ipq=1,nbuf_farr(5)
          ind=ind_send_first+ipq-1
          i=send_curvilinear_to_cartesian(ind)%i_near_neighbour(1)
          j=send_curvilinear_to_cartesian(ind)%i_near_neighbour(2)
          k=send_curvilinear_to_cartesian(ind)%i_near_neighbour(3)
          if (lpres_grad) then
            f_bufo(:,:,:,igpx_int,ipq)=f_ogrid(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,igpx)
            f_bufo(:,:,:,igpy_int,ipq)=f_ogrid(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,igpy)
          endif
          if (lchemistry .and. linterp_pressure) then
            f_bufo(:,:,:,ipp_int,ipq)=f_ogrid(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,ipp)
            f_bufo(:,:,:,iRR_int,ipq)=f_ogrid(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,iRR)
          endif
          f_bufo(:,:,:,ivar1:mvar,ipq)=f_ogrid(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,ivar1:mvar)
        enddo
        ip_bufo(1:nbuf_farr(5)) = send_curvilinear_to_cartesian(ind_send_first:ind_send_last)%ip_id
        call mpisend_int(ip_bufo(1:nbuf_farr(5)),nbuf_farr(5),send_to,send_to)
        call mpisend_real(f_bufo(:,:,:,:,1:nbuf_farr(5)),nbuf_farr,send_to,send_to+ncpus)
      endif
      ind_send_first=ind_send_last+1
    enddo
!
! Recieve from proc < iproc
!
    do iter=n_procs_recv_curv_to_cart,1,-1
      recv_from=procs_recv_curv_to_cart(iter)
      if(recv_from<iproc) then
        nbuf_farr(5)=n_ip_recv_proc_curv_to_cart(iter)
        call mpirecv_int(id_bufi(1:nbuf_farr(5)),nbuf_farr(5),recv_from,iproc)
        call mpirecv_real(f_bufi(:,:,:,:,1:nbuf_farr(5)),nbuf_farr,recv_from,iproc+ncpus)
        do ipq=1,nbuf_farr(5)
          farr = f_bufi(:,:,:,:,ipq)
          call interpolate_point_curv_to_cart(f_cartesian,id_bufi(ipq),ivar1,ivar2,farr,&
                                                       igpx_int,igpy_int,ipp_int,iRR_int)
        enddo
      endif
    enddo
!
! Send to proc < iproc
!
    ind_send_first=1
    do iter=1,n_procs_send_curv_to_cart
      ind_send_last=n_ip_to_proc_curv_to_cart(iter)+ind_send_first-1
      send_to=send_curvilinear_to_cartesian(ind_send_last)%send_to_proc
      if(send_to<iproc) then
        nbuf_farr(5)=ind_send_last-ind_send_first+1
        do ipq=1,nbuf_farr(5)
          ind=ind_send_first+ipq-1
          i=send_curvilinear_to_cartesian(ind)%i_near_neighbour(1)
          j=send_curvilinear_to_cartesian(ind)%i_near_neighbour(2)
          k=send_curvilinear_to_cartesian(ind)%i_near_neighbour(3)
          if (lpres_grad) then
            f_bufo(:,:,:,igpx_int,ipq)=f_ogrid(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,igpx)
            f_bufo(:,:,:,igpy_int,ipq)=f_ogrid(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,igpy)
          endif
          if (lchemistry .and. linterp_pressure) then
            f_bufo(:,:,:,ipp_int,ipq)=f_ogrid(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,ipp)
            f_bufo(:,:,:,iRR_int,ipq)=f_ogrid(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,iRR)
          endif
          f_bufo(:,:,:,ivar1:mvar,ipq)=f_ogrid(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,ivar1:mvar)
        enddo
        ip_bufo(1:nbuf_farr(5)) = send_curvilinear_to_cartesian(ind_send_first:ind_send_last)%ip_id
        call mpisend_int(ip_bufo(1:nbuf_farr(5)),nbuf_farr(5),send_to,send_to)
        call mpisend_real(f_bufo(:,:,:,:,1:nbuf_farr(5)),nbuf_farr,send_to,send_to+ncpus)
      endif
      ind_send_first=ind_send_last+1
    enddo
!
! Recieve from proc > iproc
!
    do iter=n_procs_recv_curv_to_cart,1,-1
      recv_from=procs_recv_curv_to_cart(iter)
      if(recv_from>iproc) then
        nbuf_farr(5)=n_ip_recv_proc_curv_to_cart(iter)
        call mpirecv_int(id_bufi(1:nbuf_farr(5)),nbuf_farr(5),recv_from,iproc)
        call mpirecv_real(f_bufi(:,:,:,:,1:nbuf_farr(5)),nbuf_farr,recv_from,iproc+ncpus)
        do ipq=1,nbuf_farr(5)
          farr = f_bufi(:,:,:,:,ipq)
          call interpolate_point_curv_to_cart(f_cartesian,id_bufi(ipq),ivar1,ivar2,farr,&
                                                       igpx_int,igpy_int,ipp_int,iRR_int)
        enddo
      endif
    enddo
!
!  Interpolate remaining points 
!
    if (lvideo .and. lwrite_slices) then
      do id=1,n_ip_curv_to_cart
      ! TODO: Make more efficient
        if(curvilinear_to_cartesian(id)%from_proc==iproc) then
          inear_loc=curvilinear_to_cartesian(id)%ind_local_neighbour
          if (lpres_grad) then
            farr(:,:,:,igpx_int)=f_ogrid(inear_loc(1)-ii1:inear_loc(1)+ii2, &
              inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,igpx)
            farr(:,:,:,igpy_int)=f_ogrid(inear_loc(1)-ii1:inear_loc(1)+ii2, &
              inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,igpy)
          endif
          if (lchemistry .and. linterp_pressure) then
            farr(:,:,:,ipp_int)=f_ogrid(inear_loc(1)-ii1:inear_loc(1)+ii2, &
              inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,ipp)
            farr(:,:,:,iRR_int)=f_ogrid(inear_loc(1)-ii1:inear_loc(1)+ii2, &
              inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,iRR)
          endif
          farr(:,:,:,ivar1:mvar)=f_ogrid(inear_loc(1)-ii1:inear_loc(1)+ii2, &
              inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,ivar1:mvar)
          call interpolate_point_curv_to_cart(f_cartesian,id,ivar1,ivar2,farr,&
                                                       igpx_int,igpy_int,ipp_int,iRR_int)
        endif
      enddo
    else
      do id=1,interpol_max
        if(curvilinear_to_cartesian(id)%from_proc==iproc) then
          inear_loc=curvilinear_to_cartesian(id)%ind_local_neighbour
          if (lpres_grad) then
            farr(:,:,:,igpx_int)=f_ogrid(inear_loc(1)-ii1:inear_loc(1)+ii2, &
              inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,igpx)
            farr(:,:,:,igpy_int)=f_ogrid(inear_loc(1)-ii1:inear_loc(1)+ii2, &
              inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,igpy)
          endif
          if (lchemistry .and. linterp_pressure) then
            farr(:,:,:,ipp_int)=f_ogrid(inear_loc(1)-ii1:inear_loc(1)+ii2, &
              inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,ipp)
            farr(:,:,:,iRR_int)=f_ogrid(inear_loc(1)-ii1:inear_loc(1)+ii2, &
              inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,iRR)
          endif
          farr(:,:,:,ivar1:mvar)=f_ogrid(inear_loc(1)-ii1:inear_loc(1)+ii2, &
              inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,ivar1:mvar)
          call interpolate_point_curv_to_cart(f_cartesian,id,ivar1,ivar2,farr,&
                                                       igpx_int,igpy_int,ipp_int,iRR_int)
        endif
      enddo
    endif
!
  endsubroutine communicate_ip_curv_to_cart
!***********************************************************************
!***********************************************************************
  subroutine interpolate_point_cart_to_curv(id,ivar1,ivar2,farr,f_cartesian,ipp_int,iRR_int)
!
    use EquationOfState, only: gamma_m1!,get_cv1
    real, dimension(mx,my,mz,mfarray), intent(in) :: f_cartesian
!
!  Use linear interpolation routine to interpolate the values on the cartesian 
!  grid to the interpolation point on the curvilinear grid
!
    integer, intent(in) :: id,ivar1,ivar2
    real, dimension(inter_len,inter_len,inter_len,ivar2-ivar1+1), intent(in) :: farr
    integer :: i,j,k
    real, dimension(3) :: xyz_ip
    integer, dimension(3) :: inear_glob
    real, dimension(ivar2-ivar1+1) :: f_ip
    integer :: ii,jj,kk
 !   real :: cv1
    integer, intent(in) :: ipp_int,iRR_int
!
    xyz_ip=cartesian_to_curvilinear(id)%xyz
    inear_glob=cartesian_to_curvilinear(id)%ind_global_neighbour
! 
!  Perform interpolation on cartesian grid
!
    if(interpolation_method==1) then
      if(.not. linear_interpolate_cartesian(farr,ivar1,ivar2,xyz_ip,inear_glob,f_ip,lcheck_interpolation)) then
        call fatal_error('linear_interpolate_cartesian','interpolation from cartesian to curvilinear')
      endif
    elseif(interpolation_method==2) then
      call interpolate_quadratic_spline(farr,ivar1,ivar2,xyz_ip,f_ip,inear_glob)
    elseif(mod(interpolation_method,2)==0) then
      if(.not. interp_lagrange(farr,ivar1,ivar2,xyz_ip,inear_glob,f_ip,.true.,.false.,lcheck_interpolation)) then
        call fatal_error('interp_lagrange','interpolation from cartesian to curvilinear')
      endif
    elseif(interpolation_method==3) then
      call interpolate_quadratic_spline(farr,ivar1,ivar2,xyz_ip,f_ip,inear_glob)
      !call interpolate_quadratic_spline(farr(:,:,:,iux:iuz),iux,iuz,xyz_ip,f_ip(iux:iuz),inear_glob)
      ! Adjust coordinates, if necessary
      if (xglobal(inear_glob(1))>xyz_ip(1)) then
        inear_glob(1) = inear_glob(1)-1
        ii=1
      else
        ii=2
      endif
      if (yglobal(inear_glob(2))>xyz_ip(2)) then
        inear_glob(2) = inear_glob(2)-1
        jj=1
      else
        jj=2
      endif
      if (zglobal(inear_glob(3))>xyz_ip(3)) then
        kk=1
        inear_glob(3) = inear_glob(3)-1
      else
        kk=2
      endif
      if(.not. linear_interpolate_cartesian(farr(ii:ii+1,jj:jj+1,2:3,4),4,4, &
              xyz_ip,inear_glob,f_ip(irho),lcheck_interpolation)) then
        call fatal_error('linear_interpolate_cartesian','interpolation from cartesian to curvilinear')
      endif
    elseif(interpolation_method==5) then
      call poly_interp_cart(ivar1,ivar2,xyz_ip,f_ip,id,f_cartesian,interpol_order_poly)
    endif
!
!  Update curvilinear grid with the new data values
!
    i=cartesian_to_curvilinear(id)%i_xyz(1)
    j=cartesian_to_curvilinear(id)%i_xyz(2)
    k=cartesian_to_curvilinear(id)%i_xyz(3)
    f_ogrid(i,j,k,iux) = f_ip(iux)*cos(y_ogrid(j)) + f_ip(iuy)*sin(y_ogrid(j))
    f_ogrid(i,j,k,iuy) = -f_ip(iux)*sin(y_ogrid(j)) + f_ip(iuy)*cos(y_ogrid(j))
    f_ogrid(i,j,k,iuz:ivar2) = f_ip(iuz:mvar)
    if (lchemistry .and. linterp_pressure) then
      f_ogrid(i,j,k,ipp)=f_ip(ipp_int)
      f_ogrid(i,j,k,iRR)=f_ip(iRR_int)
    endif
!
!   If pressure is interpolated, recover temperature field from eos
!
    if (lchemistry .and. linterp_pressure) then
        f_ogrid(i,j,k,iTT)=f_ip(ipp_int)/(f_ip(irho)*f_ip(iRR_int))
  ! does not work when no chemistry
  !      call get_cv1(cv1)
  !      f_ogrid(i,j,k,iTT)=f_ip(ipp_int)*cv1/(f_ip(irho)*gamma_m1)
    endif
!
  endsubroutine interpolate_point_cart_to_curv
!!***********************************************************************
  subroutine interpolate_point_curv_to_cart(f_cartesian,id,ivar1,ivar2,farr,igpx_int,igpy_int,ipp_int,iRR_int)
!
!  Use linear interpolation routine to interpolate the values on the cartesian 
!  grid to the interpolation point on the curvilinear grid
!
    use EquationOfState, only: lpres_grad,gamma_m1!,get_cv1
  !  use Energy, only: lpres_grad
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f_cartesian
    integer, intent(in) :: id,ivar1,ivar2
    real, dimension(inter_len,inter_len,inter_len,ivar2-ivar1+1), intent(in) :: farr
    integer :: i,j,k
    real, dimension(3) :: xyz_ip
    integer, dimension(3) :: inear_glob
    real, dimension(ivar2-ivar1+1) :: f_ip
    integer :: ii,jj,kk
!    real :: cv1
    integer :: igpx_int,igpy_int,ipp_int,iRR_int
!
    xyz_ip=curvilinear_to_cartesian(id)%xyz
    inear_glob=curvilinear_to_cartesian(id)%ind_global_neighbour
!
    if(interpolation_method==1) then
      if(.not. linear_interpolate_curvilinear(farr,ivar1,ivar2,xyz_ip,inear_glob,f_ip,lcheck_interpolation)) then
        call fatal_error('linear_interpolate_curvilinear','interpolation from curvilinear to cartesian')
      endif
    elseif(interpolation_method==2) then
      call interpolate_quadratic_sp_og(farr,ivar1,ivar2,xyz_ip,f_ip,inear_glob)
    elseif(mod(interpolation_method,2)==0) then
      if(.not. interp_lagrange(farr,ivar1,ivar2,xyz_ip,inear_glob,f_ip,.false.,.true.,lcheck_interpolation)) then
        call fatal_error('interp_lagrange','interpolation from curvilinear to cartesian')
      endif
    elseif(interpolation_method==3) then
      call interpolate_quadratic_sp_og(farr,ivar1,ivar2,xyz_ip,f_ip,inear_glob)
      !call interpolate_quadratic_sp_og(farr(:,:,:,iux:iuz),iux,iuz,xyz_ip,f_ip(iux:iuz),inear_glob)
      ! Adjust coordinates, if necessary
      if (xglobal_ogrid(inear_glob(1))>xyz_ip(1)) then
        inear_glob(1) = inear_glob(1)-1
        ii=1
      else
        ii=2
      endif
      if (yglobal_ogrid(inear_glob(2))>xyz_ip(2)) then
        inear_glob(2) = inear_glob(2)-1
        jj=1
      else
        jj=2
      endif
      if (zglobal_ogrid(inear_glob(3))>xyz_ip(3)) then 
        inear_glob(3) = inear_glob(3)-1
        kk=1
      else
        kk=2
      endif
      if(.not. linear_interpolate_curvilinear(farr(ii:ii+1,jj:jj+1,kk:kk+1,irho),irho,irho,&
            xyz_ip,inear_glob,f_ip(irho),lcheck_interpolation)) then
        call fatal_error('linear_interpolate_curvilinear','interpolation from curvilinear to cartesian ')
      endif
    elseif(interpolation_method==5) then
      if(xyz_ip(1)>=r_int_inner_poly) then 
        call poly_interp_curv(ivar1,ivar2,xyz_ip,f_ip,id,interpol_order_poly)
      else
!
!  Use linear inteprolation if too near surface to use high order polynomial
!  Preliminary solution to get inear_glob to point to bottom left corner of cell
!
        ii=floor(interpol_order_poly*0.5)
        jj=floor(interpol_order_poly*0.5)
        if(xyz_ip(1)<xglobal_ogrid(inear_glob(1))) then
          if(xyz_ip(2)<yglobal_ogrid(inear_glob(2))) then
            if(.not. linear_interpolate_curvilinear(farr(ii-1:ii,ii-1:ii,ii:ii+1,:),iux,irho,&
                  xyz_ip,(/inear_glob(1)-1,inear_glob(2)-1,inear_glob(3)/),f_ip,lcheck_interpolation)) then
              call fatal_error('linear_interpolate_curvilinear','interpolation from curvilinear to cartesian ')
            endif
          else
            if(.not. linear_interpolate_curvilinear(farr(ii-1:ii,ii:ii+1,ii:ii+1,:),iux,irho,&
                  xyz_ip,(/inear_glob(1)-1,inear_glob(2),inear_glob(3)/),f_ip,lcheck_interpolation)) then
              call fatal_error('linear_interpolate_curvilinear','interpolation from curvilinear to cartesian ')
            endif
          endif
        else
          if(xyz_ip(2)<yglobal_ogrid(inear_glob(2))) then
            if(.not. linear_interpolate_curvilinear(farr(ii:ii+1,ii-1:ii,ii:ii+1,:),iux,irho,&
                  xyz_ip,(/inear_glob(1),inear_glob(2)-1,inear_glob(3)/),f_ip,lcheck_interpolation)) then
              call fatal_error('linear_interpolate_curvilinear','interpolation from curvilinear to cartesian ')
            endif
          else
            if(.not. linear_interpolate_curvilinear(farr(ii:ii+1,ii:ii+1,ii:ii+1,:),iux,irho,&
                  xyz_ip,(/inear_glob(1),inear_glob(2),inear_glob(3)/),f_ip,lcheck_interpolation)) then
              call fatal_error('linear_interpolate_curvilinear','interpolation from curvilinear to cartesian ')
            endif
          endif
        endif
      endif
    endif
!
!  Update cartesian grid with the new data values
!
    i=curvilinear_to_cartesian(id)%i_xyz(1)
    j=curvilinear_to_cartesian(id)%i_xyz(2)
    k=curvilinear_to_cartesian(id)%i_xyz(3)
    f_cartesian(i,j,k,iux)=f_ip(iux)*cos(xyz_ip(2))-f_ip(iuy)*sin(xyz_ip(2))
    f_cartesian(i,j,k,iuy)=f_ip(iux)*sin(xyz_ip(2))+f_ip(iuy)*cos(xyz_ip(2))
    if (lchemistry .and. linterp_pressure) then
      f_cartesian(i,j,k,ipp)=f_ip(ipp_int)
      f_cartesian(i,j,k,iRR)=f_ip(iRR_int)
    endif
    if (lpres_grad) then
      f_cartesian(i,j,k,igpx)=f_ip(igpx_int)*cos(xyz_ip(2))-f_ip(igpy_int)*sin(xyz_ip(2))
      f_cartesian(i,j,k,igpy)=f_ip(igpx_int)*sin(xyz_ip(2))+f_ip(igpy_int)*cos(xyz_ip(2))
    endif
    f_cartesian(i,j,k,iuz:mvar)=f_ip(iuz:mvar)
!
!   If pressure is interpolated, recover temperature field from eos
!
    if (lchemistry .and. linterp_pressure) then
        f_cartesian(i,j,k,iTT)=f_ip(ipp_int)/(f_ip(irho)*f_ip(iRR_int))
  ! does not work when no chemistry
  !      call get_cv1(cv1)
  !      f_cartesian(i,j,k,iTT)=f_ip(ipp_int)*cv1/(f_ip(irho)*gamma_m1)

    endif
!
  endsubroutine interpolate_point_curv_to_cart
!***********************************************************************
  subroutine interp_point_curv_to_cart_alt(xyz_ip,inear_glob,ivar1,ivar2,farr,f_ip)
!
!  Use linear interpolation routine to interpolate the values on the curvilinear 
!  grid to the interpolation point on the cartesian grid
!
    real, dimension(3), intent(in) :: xyz_ip
    integer, dimension(3), intent(in) :: inear_glob
    integer, intent(in) :: ivar1,ivar2
    real, dimension(inter_len,inter_len,inter_len,ivar2-ivar1+1), intent(in) :: farr
    real, dimension(ivar2-ivar1+1), intent(out) :: f_ip
!
    if(interpolation_method==1) then
      if(.not. linear_interpolate_curvilinear(farr,ivar1,ivar2,xyz_ip,inear_glob,f_ip,lcheck_interpolation)) then
        call fatal_error('linear_interpolate_curvilinear','interpolation from curvilinear to cartesian')
      endif
    elseif(interpolation_method==2) then
      call interpolate_quadratic_sp_og(farr,ivar1,ivar2,xyz_ip,f_ip,inear_glob)
    elseif(mod(interpolation_method,2)==0) then
      if(.not. interp_lagrange(farr,ivar1,ivar2,xyz_ip,inear_glob,f_ip,.false.,.true.,lcheck_interpolation)) then
        call fatal_error('interp_lagrange','interpolation from curvilinear to cartesian')
      endif
    endif
!
  endsubroutine interp_point_curv_to_cart_alt
!***********************************************************************
  subroutine transform_curv_to_cart(f_ip,f_cartesian,id,ivar1,ivar2)
!
!  Update curvilinear grid with the new data values
!
!  02-okt-17/Jorgen: Coded
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f_cartesian
    integer, intent(in) :: id,ivar1,ivar2
    real, dimension(ivar2-ivar1+1), intent(in) :: f_ip
    integer :: i,j,k
    real, dimension(3) :: xyz_ip

    xyz_ip=curvilinear_to_cartesian(id)%xyz
    i=curvilinear_to_cartesian(id)%i_xyz(1)
    j=curvilinear_to_cartesian(id)%i_xyz(2)
    k=curvilinear_to_cartesian(id)%i_xyz(3)
    f_cartesian(i,j,k,iux)=f_ip(iux)*cos(xyz_ip(2))-f_ip(iuy)*sin(xyz_ip(2))
    f_cartesian(i,j,k,iuy)=f_ip(iux)*sin(xyz_ip(2))+f_ip(iuy)*cos(xyz_ip(2))
    f_cartesian(i,j,k,iuz:ivar2)=f_ip(iuz:ivar2)
!
  endsubroutine transform_curv_to_cart
!!***********************************************************************
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
      integer :: vari2
!
      intent(in)  :: farr, ivar1, ivar2, xxp, inear_glob, lcheck
      intent(out) :: fp

!  Set vari2 /= ivar2 to allow for farr with ivar1>1
      vari2=ivar2-ivar1+1

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
      g1=farr(1,1,1,1:vari2)
      g2=farr(2,1,1,1:vari2)
      g3=farr(1,2,1,1:vari2)
      g4=farr(2,2,1,1:vari2)
      g5=farr(1,1,2,1:vari2)
      g6=farr(2,1,2,1:vari2)
      g7=farr(1,2,2,1:vari2)
      g8=farr(2,2,2,1:vari2)
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
  logical function interp_lagrange(farr_in,ivar1,ivar2,xxp,inear_glob,fp,lcart_to_curv,lcurv_to_cart,lcheck)
!
!  Interpolate the value of f to (xp, yp) CURVILINEAR coordinate
!  using the Nth-order lagrangian interpolation.
! 
!  TODO: Extend to 3D
!  TODO: Adjust nearest point to be ACTUALLY NEAREST, not bottom left corner
!        Needed due to asymetric stencil
!
!  The coefficients are determined by the (N+1)x(N+1) grid points surrounding the
!  interpolation point.
!  Global coordinates are used for the interpolation, to allow interpolation of 
!  values outside this processors domain.
!
!  1-mar-19/Eva: Coded
!
      integer :: ivar1, ivar2, half_order
      real, dimension (3) :: xxp
      real, dimension (interpolation_method+1,interpolation_method+1,interpolation_method+1,ivar2-ivar1+1) :: farr_in
      !TODO
      real, dimension (-(interpolation_method/2):(interpolation_method/2),&
                       -(interpolation_method/2):(interpolation_method/2),ivar2-ivar1+1) :: farr
      real, dimension (ivar2-ivar1+1) :: fp
      integer, dimension (3) :: inear_glob
      logical :: lcart_to_curv, lcurv_to_cart, lcheck 
!
      intent(in)  :: farr_in, ivar1, ivar2, xxp, inear_glob, lcheck, lcart_to_curv, lcurv_to_cart
      intent(out) :: fp

      real, dimension(-(interpolation_method/2):(interpolation_method/2)) :: xglob, yglob, deltax, deltay, x_i, y_i
      real, dimension(-(interpolation_method/2):(interpolation_method/2),&
                      -(interpolation_method/2):(interpolation_method/2)) :: dx_ij, dy_ij
      real, dimension(-(interpolation_method/2):(interpolation_method/2)) :: lag
      real, dimension(-(interpolation_method/2):(interpolation_method/2),ivar2-ivar1+1) :: gp
      integer :: i,j,ix0,ix1,iy0,iy1
!
      interp_lagrange= .true.
      farr(:,:,:) = farr_in(:,:,(interpolation_method/2)+1,:)
      half_order = interpolation_method/2
!
!  Get grid points
!
      if(lcart_to_curv) then
        xglob(-half_order:half_order) = xglobal(inear_glob(1)-interpolation_method/2:inear_glob(1)+interpolation_method/2)
        yglob(-half_order:half_order) = yglobal(inear_glob(2)-interpolation_method/2:inear_glob(2)+interpolation_method/2)
      elseif(lcurv_to_cart) then
        xglob(-half_order:half_order) = xglobal_ogrid(inear_glob(1)-interpolation_method/2:inear_glob(1)+interpolation_method/2)
        yglob(-half_order:half_order) = yglobal_ogrid(inear_glob(2)-interpolation_method/2:inear_glob(2)+interpolation_method/2)
      else
        print*,'interp_lagrange: Not interpolated to any specific grid!'
        interp_lagrange= .false.
        return
      endif
!
!  Compute distance from xxp to surrounding grid points
!  Needed for checking that inear_glob is correct, and for lagrange polynomials
!
      x_i(:) = xxp(1)
      y_i(:) = xxp(2)
      deltax = x_i - xglob
      deltay = y_i - yglob
!
!  Check that inear_glob actually points to the grid point closest to xxp
!
      if(lcheck) then
        if ((any(abs(deltax)<abs(deltax(0)))) .or. (any(abs(deltay)<abs(deltay(0))))) then 
          print*, 'interp_lagrange: Interpolation point does not lie closest to center grid point.' 
          print*, 'ix0, iy0, iz0 = ', inear_glob(1:3)
          print*, 'xp, xglob(-3:3) = ', xxp(1), xglob
          print*, 'yp, yglob(-3:3) = ', xxp(2), yglob
          interp_lagrange = .false.
          return
        endif
        if (lcurv_to_cart .and. any(abs(xglob) > r_ogrid)) then
          print*,'Ghost point was used as interpolation point on curv. grid!'
          print*,'Shift the interpolation zone inside.'
          interp_lagrange = .false.
          print*,'xglob',xglob
          return
        endif
      endif
!
!  Interpolate in x-direction
! 
!  Compute distances
!
      do i = -half_order,half_order
        do j = i+1,half_order
            dx_ij(i,j) = xglob(i)-xglob(j)
            dx_ij(j,i) = -dx_ij(i,j)
            dy_ij(i,j) = yglob(i)-yglob(j)
            dy_ij(j,i) = -dy_ij(i,j)
        enddo
        dx_ij(i,i) = 0
        dy_ij(i,i) = 0
      enddo
!
!  Compute products of x-x_k/(x_i-x_k) for (i!=k)
!
      lag(:)=1
      do i = -half_order,half_order
        do j = -half_order,half_order
          if (i .ne. j) then
            lag(i) = lag(i)*deltax(j)/dx_ij(i,j)
          endif
        enddo
      enddo
!
!  Interpolate points in x-direction
!
      gp(:,:) = 0
      do i=ivar1,ivar2
        do j = -half_order,half_order
          gp(:,i) = gp(:,i)+lag(j)*farr(j,:,i)  
        enddo
      enddo
!
!  Interpolate in y-direction
! 
!  Compute distances
!
      lag(:)=1
      do i = -half_order,half_order
        do j = -half_order,half_order
          if (i .ne. j) then
            lag(i) = lag(i)*deltay(j)/dy_ij(i,j)
          endif
        enddo
      enddo
!
!  Interpolate points in y-direction
!
      do i=ivar1,ivar2
        do j=-half_order,half_order
          fp(i) = sum(lag(:)*gp(:,i))
        enddo
      enddo

      if (lcheck) then
        do i=1,ivar2-ivar1+1
          if ((fp(i)>maxval(farr(:,:,i)).and.i/=half_order) .or. (fp(i)<minval(farr(:,:,i)).and.i/=half_order)) then
!
!  Compensate for overshoots by linear interpolation
!
            ix0=0; ix1=1; iy0=0; iy1=1
            if(xglob(0)>xxp(1)) then
              ix0=half_order; ix1=half_order+1
            else
              ix0=half_order+1; ix1=half_order+2
            endif
            if(yglob(0)>xxp(2)) then
              iy0=half_order; iy1=half_order+1
            else
              iy0=half_order+1; iy1=half_order+2
            endif
            if(lcart_to_curv) then
              interp_lagrange= linear_interpolate_cartesian(farr_in(ix0:ix1,iy0:iy1,half_order+1:half_order+2,i),i,i,xxp, &
                                           (/inear_glob(1)+ix0-(half_order+1),inear_glob(2)+iy0-(half_order+1),inear_glob(3)/),&
                                           gp(0,i),lcheck_interpolation)
            else
              interp_lagrange= linear_interpolate_curvilinear(farr_in(ix0:ix1,iy0:iy1,half_order+1:half_order+2,i),i,i,xxp, &
                                           (/inear_glob(1)+ix0-(half_order+1),inear_glob(2)+iy0-(half_order+1),inear_glob(3)/),&
                                           gp(0,i),lcheck_interpolation)
            endif
!
            fp(i)=gp(0,i)
          endif
          if (fp(i)/=fp(i)) then
            print*, 'interp_interpolate: interpolated value is NaN'
            print*, 'interp_interpolate: xxp=', xxp
            print*, 'interp_interpolate: i, fp(i)=', i, fp(i)
            print*, '------------------'
            interp_lagrange=.false.
          endif
        enddo
      endif
!
  endfunction interp_lagrange
!!***********************************************************************
!!***********************************************************************
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
      integer :: vari2
!
      intent(in)  :: farr, ivar1, ivar2, xxp, inear_glob, lcheck
      intent(out) :: fp
!
!  Set vari2 /= ivar2 to allow for farr with ivar1>1
      vari2 = ivar2-ivar1+1
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
      g1=farr(1,1,1,1:vari2)
      g2=farr(2,1,1,1:vari2)
      g3=farr(1,2,1,1:vari2)
      g4=farr(2,2,1,1:vari2)
      g5=farr(1,1,2,1:vari2)
      g6=farr(2,1,2,1:vari2)
      g7=farr(1,2,2,1:vari2)
      g8=farr(2,2,2,1:vari2)
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
            print*, 'linear_interpolate_curvilinear: ix0, iy0, iz0=', ix0,iy0,iz0
            print*, 'linear_interpolate_curvilinear: x0, y0, z0=', &
                xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
            print*, 'linear_interpolate_curvilinear: x0+1, y0+1, z0+1=', &
                xglobal_ogrid(ix0+1), yglobal_ogrid(iy0+1), zglobal_ogrid(iz0+1)
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
  subroutine close_interpolation(f,ix0_,iy0_,iz0_,iobj,xxp,f_tmp,fluid_point)
!
    real, dimension(:,:,:,:) :: f
    integer, intent(in) :: ix0_, iy0_, iz0_, iobj
    real, dimension(:), intent(inout) :: f_tmp
    real, dimension(3), intent(in) :: xxp
    logical, intent(in) :: fluid_point
!
    call keep_compiler_quiet(f)
    call keep_compiler_quiet(f_tmp)
    call keep_compiler_quiet(ix0_,iy0_,iz0_)
    call keep_compiler_quiet(fluid_point)
    call keep_compiler_quiet(xxp)
    call keep_compiler_quiet(iobj)
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
    real :: x0,y0,z0
    real :: xi1lo,xi1up,g1lo,g1up
    real :: xi2lo,xi2up,g2lo,g2up
    real :: xi3lo,xi3up,g3lo,g3up
    real :: xi1star,xi2star,xi3star
!
    real, dimension(mx_ogrid) :: g1,g1der1,g1der2,xi1,xprim2_ogrid
    real, dimension(my_ogrid) :: g2,g2der1,g2der2,xi2,yprim2_ogrid
    real, dimension(mz_ogrid) :: g3,g3der1,g3der2,xi3,zprim2_ogrid
!
    real, dimension(0:2*nprocx+1) :: xi1proc,g1proc
    real, dimension(0:2*nprocy+1) :: xi2proc,g2proc
    real, dimension(0:2*nprocz+1) :: xi3proc,g3proc
!
    real :: a
    integer :: i
!
    lequidist_ogrid=(grid_func_ogrid=='linear')
!
!  Abbreviations
!
    x0 = xyz0_ogrid(1)
    y0 = xyz0_ogrid(2)
    z0 = xyz0_ogrid(3)
    Lx_og = Lxyz_ogrid(1)
    Ly_og = Lxyz_ogrid(2)
    Lz_og = Lxyz_ogrid(3)
!
!  Set the lower boundary and the grid size.
!
    x00 = x0
    y00 = y0
    z00 = z0
!
    dx_ogrid = Lx_og / merge(nxgrid_ogrid, max(nxgrid_ogrid-1,1), .false.)
    dy_ogrid = Ly_og / merge(nygrid_ogrid, max(nygrid_ogrid-1,1), .true.)
    dz_ogrid = Lz_og / merge(nzgrid_ogrid, max(nzgrid_ogrid-1,1), lperi(3))
!
!  Shift the lower boundary if requested, but only for periodic directions.
!
      if (lshift_origin_ogrid(1)) x00 = x0 + 0.5 * dx_ogrid
      if (lshift_origin_ogrid(2)) y00 = y0 + 0.5 * dy_ogrid
      if (lshift_origin_ogrid(3)) z00 = z0 + 0.5 * dz_ogrid
!
!  Shift the lower boundary if requested, but only for periodic directions.
!  Contrary to the upper case (lshift_origin)
!
      if (lshift_origin_lower_ogrid(1)) x00 = x0 - 0.5 * dx_ogrid
      if (lshift_origin_lower_ogrid(2)) y00 = y0 - 0.5 * dy_ogrid
      if (lshift_origin_lower_ogrid(3)) z00 = z0 - 0.5 * dz_ogrid
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
    real :: dxmin_x, dxmin_y, dxmin_z, dxmax_x, dxmax_y, dxmax_z, dxmax_r
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
    if (lequidist_ogrid(1) .or. nxgrid_ogrid <= 1) then
      dxmin_x = dx_ogrid
      dxmax_x = dx_ogrid
    else
      dxmin_x = minval(xprim_ogrid(l1_ogrid:l2_ogrid))
      dxmax_x = maxval(xprim_ogrid(l1_ogrid:l2_ogrid))
    endif
!
    if (lequidist_ogrid(2) .or. nygrid_ogrid <= 1) then
      dxmin_y = dy_ogrid*minval(x_ogrid(l1_ogrid:l2_ogrid))
      dxmax_y = dy_ogrid*maxval(x_ogrid(l1_ogrid:l2_ogrid))
    else
      dxmin_y = minval(yprim_ogrid(m1_ogrid:m2_ogrid))
      dxmax_y = maxval(yprim_ogrid(m1_ogrid:m2_ogrid))
    endif
!
    if (lequidist_ogrid(3) .or. nzgrid_ogrid <= 1) then
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
!  Maximum radial grid spacing needed for interpolation region
!
    drmax_ogrid=dxmax_x
    call mpiallreduce_max(drmax_ogrid,dxmax_r)
    drmax_ogrid=dxmax_r
!
    dxmax_ogrid = maxval( (/dxmax_x, dxmax_y, dxmax_z, epsilon(dx_ogrid)/), &
              MASK=((/nxgrid_ogrid, nygrid_ogrid, nzgrid_ogrid, 2/) > 1) )
!
    call mpiallreduce_max(dxmax_ogrid,dxmax_x)
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
!                     Only cylindrical coordinates included
!
    if (lpencil_ogrid(i_og_x_mn))     &
        p_ogrid%x_mn    = x_ogrid(l1_ogrid:l2_ogrid)*cos(y_ogrid(m_ogrid))
    if (lpencil_ogrid(i_og_y_mn))     &
        p_ogrid%y_mn    = x_ogrid(l1_ogrid:l2_ogrid)*sin(y_ogrid(m_ogrid))
    if (lpencil_ogrid(i_og_z_mn))     &
        p_ogrid%z_mn    = spread(z_ogrid(n_ogrid),1,nx_ogrid)
    if (lpencil_ogrid(i_og_rcyl_mn))  &
        p_ogrid%rcyl_mn = x_ogrid(l1_ogrid:l2_ogrid)
    if (lpencil_ogrid(i_og_phi_mn))   &
        p_ogrid%phi_mn  = spread(y_ogrid(m_ogrid),1,nx_ogrid)
    if (lpencil_ogrid(i_og_rcyl_mn1)) &
        p_ogrid%rcyl_mn1= 1./max(p_ogrid%rcyl_mn,tini)
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
!  Do the same for the 1/dx arrays
!
    dx1global_ogrid(nghost+1:mxgrid_ogrid-nghost) = dx1grid_ogrid
    dy1global_ogrid(nghost+1:mygrid_ogrid-nghost) = dy1grid_ogrid
    dz1global_ogrid(nghost+1:mzgrid_ogrid-nghost) = dz1grid_ogrid
!
    dx1global_ogrid(1:nghost) = dx_1_ogrid(1:nghost)
    dy1global_ogrid(1:nghost) = dy_1_ogrid(1:nghost)
    dz1global_ogrid(1:nghost) = dz_1_ogrid(1:nghost)
!
    dx1global_ogrid(mxgrid_ogrid-nghost+1:mxgrid_ogrid) = dx_1_ogrid(mx_ogrid-nghost+1:mx_ogrid)
    dy1global_ogrid(mygrid_ogrid-nghost+1:mygrid_ogrid) = dy_1_ogrid(my_ogrid-nghost+1:my_ogrid)
    dz1global_ogrid(mzgrid_ogrid-nghost+1:mzgrid_ogrid) = dz_1_ogrid(mz_ogrid-nghost+1:mz_ogrid)
!
    call mpibcast_real(dx1global_ogrid(1:nghost), nghost, iproc_first)
    call mpibcast_real(dy1global_ogrid(1:nghost), nghost, iproc_first)
    call mpibcast_real(dz1global_ogrid(1:nghost), nghost, iproc_first)
!
    call mpibcast_real(dx1global_ogrid(mxgrid_ogrid-nghost+1:mxgrid_ogrid), nghost, iproc_last)
    call mpibcast_real(dy1global_ogrid(mygrid_ogrid-nghost+1:mygrid_ogrid), nghost, iproc_last)
    call mpibcast_real(dz1global_ogrid(mzgrid_ogrid-nghost+1:mzgrid_ogrid), nghost, iproc_last)
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
    use Boundcond, only: update_ghosts
    use EquationOfState, only: lpres_grad
 !   use Energy, only: lpres_grad
!
    real, dimension (mx,my,mz,mfarray) :: f_cartesian
    real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df_ogrid
    real :: dt_ogrid
    integer :: tstep_ogrid
    integer :: j,jj=0, tss=0
    real, dimension(3) :: alpha_ts_ogrid=0.,beta_ts_ogrid=0.,dt_beta_ts_ogrid=0.

  !call  run_tests_ogrid
    if(.not.lrk_tvd) then
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
    endif
!
!  Interpolate data from cartesian to curvilinear grid.
!  Before interpolating, necessary points outside this processors domain are
!  recieved from appropriate processor
!
    call update_ghosts(f_cartesian,1,mvar)
    if (lchemistry .and. linterp_pressure) then
      call communicate_ip_cart_to_curv(f_cartesian,1,mvar+2)
    else
      call communicate_ip_cart_to_curv(f_cartesian,1,mvar)
    endif
!
    dt_ogrid = dt/timestep_factor
    !print*, 'dt_ogrid', dt_ogrid
    !print*, 'convective timestep', dxmin_ogrid/maxval(f_ogrid(:,:,4,iux:iuy))
    !print*, 'viscous timestep', dxmin_ogrid**2/1.e-3
    if(.not. lrk_tvd) dt_beta_ts_ogrid=dt_ogrid*beta_ts_ogrid
!
!  Perform a number of timestep equal to timestep_factor, such that the
!  endtime t_ogrid is equal to t after the timesteps
!
    tss = tss + 1 
    do tstep_ogrid=1,timestep_factor
!
!  Set up df for each time sub.
!
      df_ogrid=0.0
      if(lrk_tvd) then
!  First subtimestep
        call pde_ogrid(f_ogrid,df_ogrid,dt_ogrid)
        do j=1,mvar 
          f_tmp(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j) = &
              f_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j) &
            + dt_ogrid*df_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j)
        enddo
        f_tmp(l2_ogrid+1:mx_ogrid,:,:,:)=f_ogrid(l2_ogrid+1:mx_ogrid,:,:,:)
!  Second subtimestep
        df_ogrid=0.0
        call pde_ogrid(f_tmp,df_ogrid,dt_ogrid)
        do j=1,mvar 
          f_tmp(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j) = &
              (3./4.)*f_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j) &
            +(1./4.)*(f_tmp(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j) &
                      +dt_ogrid*df_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j))
        enddo
!  Third subtimestep
        llast_ogrid=(tstep_ogrid==timestep_factor)
        df_ogrid=0.0
        call pde_ogrid(f_tmp,df_ogrid,dt_ogrid)
        do j=1,mvar 
          f_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j) = &
              (1./3.)*f_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j) &
            +(2./3.)*(f_tmp(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j) &
                      +dt_ogrid*df_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j))
        enddo
      else
        do itsub=1,itorder
          df_ogrid=alpha_ts_ogrid(itsub)*df_ogrid
          llast_ogrid=(tstep_ogrid==timestep_factor).and.(itsub==itorder)
!
!  Change df according to the chosen physics modules.
!
          call pde_ogrid(f_ogrid,df_ogrid,dt_ogrid)
!
!  Time evolution of grid variables.
!
          do j=1,mvar 
            f_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j) = &
                f_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j) &
              + dt_beta_ts_ogrid(itsub)*df_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j)
          enddo
        enddo
      endif
  !    if(lfilter_solution) then
  !      call communicate_filter_zones(f_ogrid,f_filterH_lowerx,f_filterH_upperx,f_filterH_lowery,f_filterH_uppery)
  !      call pade_filter(f_ogrid)
  !      call update_ghosts_ogrid(f_ogrid)
  !    endif
    enddo
!
!  Interpolate data from curvilinear to cartesian grid.
!  Before interpolating, necessary points outside this processors domain are
!  recieved from appropriate processor
!
    call update_ghosts_ogrid(f_ogrid)
!
!  Filter solution if this option is set
!
    if(mod(tss,filter_frequency)==0 .and. lfilter_solution) then
      call communicate_filter_zones(f_ogrid,f_filterH_lowerx,f_filterH_upperx,f_filterH_lowery,f_filterH_uppery)
      call pade_filter(f_ogrid)
      call update_ghosts_ogrid(f_ogrid)
    endif

    if (lpres_grad .and. (lchemistry .and. linterp_pressure)) then
      call communicate_ip_curv_to_cart(f_cartesian,1,mvar+4)
    elseif (lpres_grad .or. (lchemistry .and. linterp_pressure)) then
      call communicate_ip_curv_to_cart(f_cartesian,1,mvar+2)
    else
      call communicate_ip_curv_to_cart(f_cartesian,1,mvar)
    endif
! 
!     !TODO: Should use the particle flow info in the interpolation point
!     !      computation above
     if(lparticles)  call update_ogrid_flow_info(ivar1_part,ivar2_part)
!
    call wsnap_ogrid('OGVAR',ENUM=.true.,FLIST='ogvarN.list')
!
!  Set silly values for f_cartesian inside r_int_inner
!
! for debugging
! do ii=l1,l2
! do jj=m1,m2
! do kk=n1,n2
!       if(radius_ogrid(x(ii),y(jj))<(r_int_inner)) then
!         f_cartesian(ii,jj,:,:)=huge_real
!       endif
! enddo
! enddo
! enddo

  endsubroutine time_step_ogrid
!***********************************************************************
    subroutine pde_ogrid(f_og,df,dt_ogrid)
!
!  06-feb-17/Jorgen+Nils: Adapded from equ.f90
!
!  Call the different evolution equations.
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df
      integer :: nyz_ogrid
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
      intent(out)    :: df
      real :: c_dragx,c_dragy, Nusselt, dt_ogrid
      c_dragx=0.
      c_dragy=0.
      Nusselt=0.
!
!  Initiate communication and do boundary conditions.
!
      call boundconds_x_ogrid(f_og)
      call update_ghosts_ogrid(f_og)
      if (lchemistry) call chemspec_normalization_N2_og(f_og)
!
      if (lchemistry) call calc_for_chem_mixture_ogrid(f_og)
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
        call calc_pencils_hydro_ogrid(f_og)
        call calc_pencils_density_ogrid(f_og)
        if(.not. lchemistry) then
          call calc_pencils_eos_ogrid(f_og)
        else
          call calc_pencils_eos_ogrid_chem(f_og)
        endif
        call calc_pencils_energy_ogrid(f_og)
        call calc_pencils_viscosity_ogrid
        if (lchemistry) then
          call calc_pencils_chemistry_ogrid(f_og)
          if (lreac_heter) call calc_heter_reaction_term(f_og)
        endif
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
        call denergy_dt_ogrid(df,f_og)
        if (lchemistry) call dYk_dt_ogrid(f_og,df,dt_ogrid)
!
!  Compute drag and lift coefficient, if this is the last sub-timestep
!
        if(llast_ogrid.and.lfirst_proc_x) then
          if ((idiag_c_dragx/=0).or.(idiag_c_dragy/=0)) then
            call drag_force_pencils(c_dragx,c_dragy)
          endif
          if (idiag_Nusselt/=0) then
            call Nusselt_pencils(Nusselt)
          endif
        endif
!
!  End of loops over m and n.
!
      enddo mn_loop
!
      if(llast_ogrid) then
        if ((idiag_c_dragx/=0).or.(idiag_c_dragy/=0)) then
          call drag_coeffs(c_dragx,c_dragy)
        endif
        if ((idiag_Nusselt/=0)) then
          call Nusselt_coeffs(Nusselt)
        endif
      endif

!
!  -------------------------------------------------------------
!  NO CALLS MODIFYING DF BEYOND THIS POINT (APART FROM FREEZING)
!  -------------------------------------------------------------
!  Frerzing must be done after the full (m,n) loop, as df may be modified
!  outside of the considered pencil.
!
!  Freezing boundary conditions in x (radial direction), only on points
!  at the surface
!
      if(lfirst_proc_x) then
        do imn_ogrid=1,nyz_ogrid
          n_ogrid=nn_ogrid(imn_ogrid)
          m_ogrid=mm_ogrid(imn_ogrid)
          df(l1_ogrid,m_ogrid,n_ogrid,iux:iuz) = 0.
          if (lexpl_rho) df(l1_ogrid,m_ogrid,n_ogrid,irho) = 0.
          if (lchemistry .and. lreac_heter) &
             df(l1_ogrid,m_ogrid,n_ogrid,ichemspec(1):ichemspec(nchemspec)) = 0.
        enddo
        if (iTT .gt. 0) df(l1_ogrid,:,:,iTT) = 0.
      endif
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
    subroutine denergy_dt_ogrid(df,f_og)
!
!  Calculate pressure gradient term for isothermal/polytropic equation
!  of state.
!
      use EquationOfState, only: gamma_m1
!
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid), intent(in) ::  f_og
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df
      integer :: j
      intent(inout) :: df
!
!  Add isothermal/polytropic pressure term in momentum equation.
!
   !   if (.not. lchemistry) then
          df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iux:iuz)= &
              df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iux:iuz)+p_ogrid%fpres
    !  else
! TODO: pressure gradient term when chemistry
     !     df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iux:iuz)= &
     !         df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iux:iuz)-p_ogrid%rho1gpp
     ! endif
!
!  Solve Energy equation (in case of non-isothermal equation of state)
!
      if (iTT .ne. 0) then
!
!  Advection term and PdV-work.
!
        if (ladvection_temperature) then
          df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT) = &
              df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT) - p_ogrid%ugTT
        endif
!
!  Add divu term.
!  If lchemistry divu is added in dYk_dt
!
        if (ldensity .and. (.not. lchemistry)) then
          df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT) = &
              df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT) - &
              gamma_m1*p_ogrid%TT*p_ogrid%divu
        endif
!
!  Calculate viscous contribution to temperature.
!
!      if (lviscosity.and.lviscosity_heat) call calc_viscous_heat(df,p,Hmax)
!
!  Thermal conduction
!
        if (lheatc_chiconst) then
          call calc_heatcond_constchi_ogrid(df,p_ogrid)
        elseif (lheatc_chemistry) then
          call calc_heatcond_chemistry_ogrid(f_og,df)
        else
          call fatal_error('denergy_dt_ogrid','Must use lheatc_chicons=T or lheatc_chemistry=T')
        endif
!
      endif
!
    endsubroutine denergy_dt_ogrid
!***********************************************************************
    subroutine calc_heatcond_constchi_ogrid(df,p_ogrid)
!
!  Calculate the radiative diffusion term for constant chi:
!  lnTT version: cp*chi*Div(rho*TT*glnTT)/(rho*cv*TT)
!           = gamma*chi*(g2.glnTT+g2lnTT) where g2=glnrho+glnTT
!    TT version: cp*chi*Div(rho*gTT)/(rho*cv)
!           = gamma*chi*(g2.gTT+g2TT) where g2=glnrho
!
!  17-aug-17/ewa+nils: adapted from temperature_idealgas
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: df
      type (pencil_case_ogrid) :: p_ogrid
      real, dimension (nx_ogrid) :: g2
!
      intent(in) :: p_ogrid
      intent(inout) :: df
!
      call dot(p_ogrid%glnrho,p_ogrid%gTT,g2)
      g2=g2+p_ogrid%del2TT
!
!  Add heat conduction to RHS of temperature equation.
!
      df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT) = &
          df(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT) + gamma*chi*g2
!
    endsubroutine calc_heatcond_constchi_ogrid
!***********************************************************************
    subroutine calc_pencils_hydro_ogrid(f_og)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
! Pencils: uu, u2, uij, divu, sij, sij2, ugu, ugu2, del2u
!
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og

      if (lpencil_ogrid(i_og_uu)) p_ogrid%uu=f_og   (l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iux:iuz)
      if (lpencil_ogrid(i_og_u2)) call dot2_mn_ogrid(p_ogrid%uu,p_ogrid%u2)
      if (lpencil_ogrid(i_og_uij)) call gij_ogrid(f_og   ,iuu,p_ogrid%uij)!,1)
      if (lpencil_ogrid(i_og_divu)) call div_mn_ogrid(p_ogrid%uij,p_ogrid%divu,p_ogrid%uu)
      if (lpencil_ogrid(i_og_sij)) call traceless_strain_ogrid(p_ogrid%uij,p_ogrid%divu,p_ogrid%sij,p_ogrid%uu)
      if (lpencil_ogrid(i_og_sij2)) call multm2_sym_mn_ogrid(p_ogrid%sij,p_ogrid%sij2)
      if (lpencil_ogrid(i_og_ugu)) call u_dot_grad_ogrid(f_og   ,iuu,p_ogrid%uij,p_ogrid%uu,p_ogrid%ugu)
      if (lpencil_ogrid(i_og_ugu2)) call dot2_mn_ogrid(p_ogrid%ugu,p_ogrid%ugu2)
      if (lpencil_ogrid(i_og_graddivu).and.lpencil_ogrid(i_og_del2u)) then
        call gij_etc_ogrid(f_og   ,iuu,p_ogrid%uu,p_ogrid%uij,DEL2=p_ogrid%del2u,GRADDIV=p_ogrid%graddivu)
      elseif(lpencil_ogrid(i_og_graddivu)) then
        call gij_etc_ogrid(f_og   ,iuu,p_ogrid%uu,p_ogrid%uij,GRADDIV=p_ogrid%graddivu)
      elseif(lpencil_ogrid(i_og_del2u)) then
        call gij_etc_ogrid(f_og   ,iuu,p_ogrid%uu,p_ogrid%uij,DEL2=p_ogrid%del2u)
      endif
!
    endsubroutine calc_pencils_hydro_ogrid
!***********************************************************************
    subroutine calc_pencils_density_ogrid(f_og)
!
!  Calculate Density pencils for linear density.
!  Most basic pencils should come first, as others may depend on them.
!
!  10-feb-17/Jorgen+Nils: Adapted from calc_pencils_linear_density in density.f90
!
      use density, only:lupw_lnrho
!
      integer :: i
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
!
! Pencils: rho, rho1, lnrho, glnrho, grho, ugrho, sglnrho
!
      p_ogrid%rho=f_og   (l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,irho)
      if (lpencil_ogrid(i_og_rho1)) p_ogrid%rho1=1.0/p_ogrid%rho
      if (lpencil_ogrid(i_og_lnrho)) p_ogrid%lnrho=log(p_ogrid%rho)
      if (lpencil_ogrid(i_og_glnrho)) then
        call grad_ogrid(f_og   ,irho,p_ogrid%grho)
        do i=1,3
          p_ogrid%glnrho(:,i)=p_ogrid%rho1*p_ogrid%grho(:,i)
        enddo
      endif
      if (lpencil_ogrid(i_og_ugrho)) call u_dot_grad_ogrid(f_og   ,irho,p_ogrid%grho,p_ogrid%uu,p_ogrid%ugrho,UPWIND=lupw_lnrho)
      if (lpencil_ogrid(i_og_sglnrho)) call multmv_mn_ogrid(p_ogrid%sij,p_ogrid%glnrho,p_ogrid%sglnrho)
!
    endsubroutine calc_pencils_density_ogrid
!***********************************************************************
    subroutine calc_pencils_eos_ogrid(f_og)
!
!  Calculate EquationOfState pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  10-feb-17/Jorgen+Nils: Adapted from calc_pencils_eos_pencpar in eos_idealgas.f90
!
      use EquationOfState, only: get_cv1,get_cp1,cs20,gamma_m1,lnrho0
!
      real :: cp1, cv1, cp, cv
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
!
!  Inverse cv and cp values.
!
      call get_cp1(cp1)
      call get_cv1(cv1)
      cp=1./cp1
      cv=1./cv1
!
      if (iTT .ne. 0) then
        if (lpencil_ogrid(i_og_TT)) &
            p_ogrid%TT=f_og   (l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT)
        if (lpencil_ogrid(i_og_lnTT)) &
            p_ogrid%lnTT=log(f_og   (l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT))
        if (lpencil_ogrid(i_og_gTT)) call grad_ogrid(f_og   ,iTT,p_ogrid%gTT)
        if (lpencil_ogrid(i_og_cs2))  &
            p_ogrid%cs2=cp*gamma_m1*f_og   (l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,iTT)
        if (lpencil_ogrid(i_og_pp)) p_ogrid%pp=cv*gamma_m1*p_ogrid%rho*p_ogrid%TT
        if (lpencil_ogrid(i_og_del2TT)) call del2_ogrid(f_og   ,iTT,p_ogrid%del2TT)
      else
        if (leos_isentropic) then
          if (lpencil_ogrid(i_og_ss)) p_ogrid%ss=0.0
          if (lpencil_ogrid(i_og_cs2)) &
              p_ogrid%cs2=cs20*exp(gamma_m1*(p_ogrid%lnrho-lnrho0))
          if (lpencil_ogrid(i_og_lnTT)) &
              p_ogrid%lnTT=lnTT0+cv1*p_ogrid%ss+gamma_m1*(p_ogrid%lnrho-lnrho0)
          if (lpencil_ogrid(i_og_pp)) &
              p_ogrid%pp=(cp-cv)*exp(p_ogrid%lnTT+p_ogrid%lnrho)
        elseif (leos_isothermal) then
          if (lpencil_ogrid(i_og_ss)) p_ogrid%ss=-(cp-cv)*(p_ogrid%lnrho-lnrho0)
          if (lpencil_ogrid(i_og_cs2)) p_ogrid%cs2=cs20
          if (lpencil_ogrid(i_og_pp)) p_ogrid%pp=p_ogrid%cs2*p_ogrid%rho
        endif
      endif

   !   if (linterp_pressure) then
   !     f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,ipp) = p_ogrid%pp
   !   endif

    endsubroutine calc_pencils_eos_ogrid
!***********************************************************************
    subroutine calc_pencils_energy_ogrid(f_og)
!
!  Calculate Energy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      use EquationOfState, only: gamma1, lpres_grad
!      use Energy, only: lpres_grad
!
      integer :: j
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
!
!  Pencils: fpres (=pressure gradient force)
!
      if (iTT .ne. 0) then
        if (lpencil_ogrid(i_og_ugTT)) then
          call u_dot_grad_ogrid(f_og   ,iTT,p_ogrid%gTT,p_ogrid%uu,&
              p_ogrid%ugTT,UPWIND=lupw_lnTT)
        endif
        if ((lpencil_ogrid(i_og_fpres)) .and. (.not. lchemistry)) then
          do j=1,3
            p_ogrid%fpres(:,j)=-gamma1*p_ogrid%cs2*&
                (p_ogrid%glnrho(:,j)+p_ogrid%gTT(:,j)/p_ogrid%TT)
          enddo
        endif
      else
         if (lpencil_ogrid(i_og_fpres)) then
            do j=1,3
               p_ogrid%fpres(:,j)=-p_ogrid%cs2*p_ogrid%glnrho(:,j)
            enddo
         endif
      endif
!
! Store pressure gradient as auxillary if requsted (this is done elswhere
! for cases with chemistry).
!
      if (lpres_grad .and. (.not. lchemistry)) then
        f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,igpx) = &
            -p_ogrid%fpres(:,1)*p_ogrid%rho
        f_og(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,igpy) = &
            -p_ogrid%fpres(:,2)*p_ogrid%rho
      endif
!
    endsubroutine calc_pencils_energy_ogrid
!***********************************************************************
    subroutine calc_pencils_viscosity_ogrid
!
!  Calculate Viscosity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  10-feb-17/Jorgen+Nils: Adapted from rountine in viscosity.f90
!
      use viscosity, only:getnu
      real :: nu
      integer :: j
      real, dimension (nx_ogrid,3) :: sgradnu
!      
!  Viscous force and viscous heating are calculated here (for later use).
!
      p_ogrid%fvisc=0.0                              
!
      if (.not. lchemistry) then
        call getnu(nu_input=nu)
        do j=1,3
          p_ogrid%fvisc(:,j) = p_ogrid%fvisc(:,j) + nu*(2*p_ogrid%sglnrho(:,j)+p_ogrid%del2u(:,j) &
                             + 1./3.*p_ogrid%graddivu(:,j))
          !p_ogrid%fvisc(:,j) = p_ogrid%fvisc(:,j) + nu*(p_ogrid%del2u(:,j))
        enddo
      else
        call multmv_mn_ogrid(p_ogrid%sij,p_ogrid%gradnu,sgradnu)
        do j=1,3
          p_ogrid%fvisc(:,j)=p_ogrid%nu*(p_ogrid%del2u(:,j)+1./3.*p_ogrid%graddivu(:,j)) + 2.*sgradnu(:,j)
        enddo
        if (ldensity) then
          do j=1,3
            p_ogrid%fvisc(:,j)=p_ogrid%fvisc(:,j) + 2.*p_ogrid%nu*p_ogrid%sglnrho(:,j)
          enddo
        endif
!
!  Viscous heating and time step.
!
!        if (lpencil(i_visc_heat)) p%visc_heat=p%visc_heat+2*p%nu*p%sij2
!        if (lfirst .and. ldt) p%diffus_total=p%diffus_total+p%nu
!
      endif
!
    end subroutine calc_pencils_viscosity_ogrid
!***********************************************************************
!  FROM BOUNDCOND.F90
!***********************************************************************
!  ROUTINES
!    boundconds_x_ogrid
!    boundconds_y_ogrid
!    boundconds_z_ogrid
!***********************************************************************
    subroutine boundconds_x_ogrid(f_og)
!
!  Boundary conditions at the cylinder surface in the radial direction (xdir).
!  For ogrids, only boundary conditions at cylinder surface is set. The BC on
!  the 'top' is set by interpolation from cartesian grid, outside the timestep.
!  Only need to compute boundary value for the density, using stencil that
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
      use density, only:lupw_lnrho
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
      integer :: k
!
      if(lfirst_proc_x) then
        if(SBP) then
          ! chemistry BCs are taken care of here
          call bval_from_neumann_SBP(f_og)
        elseif(BDRY5) then
          if(lexpl_rho) call bval_from_neumann_bdry5(f_og)
          if (lchemistry) call fatal_error('boundconds_x_ogrid', &
          'chemistry BCs set correctly only when SBP=T')
        else
          !TODO set f_og as input
          call set_ghosts_onesided_ogrid(iux)
          call set_ghosts_onesided_ogrid(iuy)
          call set_ghosts_onesided_ogrid(iuz)
          if (lfilter_TT) then
             call set_ghosts_onesided_ogrid(iTT)
          endif
          call bval_from_neumann_arr_ogrid
          call set_ghosts_onesided_ogrid(irho)
          if (lchemistry) then
            do k = 1,nchemspec
              call set_ghosts_onesided_ogrid(ichemspec(k))
            enddo
          endif
          if (lchemistry) call fatal_error('boundconds_x_ogrid', &
          'chemistry BCs set correctly only when SBP=T')
        endif
        !if(lupw_lnrho) then
        !  if(lexpl_rho) call bval_from_neumann_upw_ogrid
        !  call set_ghosts_onesided_upw_ogrid(irho)
        !endif
      endif
!
    endsubroutine boundconds_x_ogrid
!***********************************************************************
    subroutine boundconds_y_ogrid(f_og)
!
!  Periodic boundary condition for runs with a single processor in y-direction 
!
!  06-feb-17/Jorgen: Adapted from boundcond.f90 to be used for ogrids where the
!                    y-dir is always periodic
!
      integer :: ivar1, ivar2, j
      integer :: m1i_ogrid=m1_ogrid+nghost-1
      integer :: m2i_ogrid=my_ogrid-2*nghost+1
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
      !
      ! JONAS i did this change
      !
      ivar1=1; ivar2=min(mcom,size(f_og,4))
!
!  Boundary conditions in y
!  Periodic, with y being the theta direction for the cylindrical grid
!
      do j=ivar1,ivar2
!  Bottom boundary
        f_og   (:,1:m1_ogrid-1,:,j) = f_og   (:,m2i_ogrid:m2_ogrid,:,j)
!  Top boundary
        f_og   (:,m2_ogrid+1:,:,j) = f_og   (:,m1_ogrid:m1i_ogrid,:,j)
      enddo
!
      if (lchemistry .and. linterp_pressure) then
        f_og   (:,1:m1_ogrid-1,:,ipp) = f_og   (:,m2i_ogrid:m2_ogrid,:,ipp)
        f_og   (:,m2_ogrid+1:,:,ipp) = f_og   (:,m1_ogrid:m1i_ogrid,:,ipp)
        if (lchemistry) then
          f_og   (:,1:m1_ogrid-1,:,iRR) = f_og   (:,m2i_ogrid:m2_ogrid,:,iRR)
          f_og   (:,m2_ogrid+1:,:,iRR) = f_og   (:,m1_ogrid:m1i_ogrid,:,iRR)
        endif
      endif

    endsubroutine boundconds_y_ogrid
!***********************************************************************
    subroutine boundconds_z_ogrid(f_og)
!
!  Periodic boundary condition for 3D-runs with a single processor in z-direction 
!
!  06-feb-17/Jorgen: Adapted from boundcond.f90 to be used for ogrids where the
!                    z-dir is always periodic as long as nzgrid=/1
!
      integer :: ivar1, ivar2, j
      integer :: n1i_ogrid=n1_ogrid+nghost-1
      integer :: n2i_ogrid=mz_ogrid-2*nghost+1
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
!
      ivar1=1; ivar2=min(mcom,size(f_og,4))
!
!  Boundary conditions in z
!
      do j=ivar1,ivar2
!  Bottom boundary
        f_og   (:,:,1:n1_ogrid-1,j) = f_og   (:,:,n2i_ogrid:n2_ogrid,j)
!  Top boundary
        f_og   (:,:,n2_ogrid+1:,j) = f_og   (:,:,n1_ogrid:n1i_ogrid,j)
      enddo
!
      if (lchemistry .and. linterp_pressure) then
        f_og   (:,:,1:n1_ogrid-1,ipp) = f_og   (:,:,n2i_ogrid:n2_ogrid,ipp)
        f_og   (:,:,n2_ogrid+1:,ipp) = f_og   (:,:,n1_ogrid:n1i_ogrid,ipp)
        if (lchemistry) then
          f_og   (:,:,1:n1_ogrid-1,iRR) = f_og   (:,:,n2i_ogrid:n2_ogrid,iRR)
          f_og   (:,:,n2_ogrid+1:,iRR) = f_og   (:,:,n1_ogrid:n1i_ogrid,iRR)
        endif  
      endif
!
    endsubroutine boundconds_z_ogrid
!***********************************************************************
    subroutine boundconds_y_filter(f_og,f_Hloy,f_Hupy,Hsize)
!
!  Periodic boundary condition for runs with a single processor in y-direction 
!  Extended ghosts zones used if filtering is on
!
!  29-feb-17/Jorgen: Coded
!
      integer :: ivar1, ivar2, j
      integer :: m1i_ogrid=m1_ogrid+nghost-1
      integer :: m2i_ogrid=my_ogrid-2*nghost+1
      integer, intent(in) :: Hsize
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid) :: f_og
      real, dimension (mx_ogrid,filter_Hsize,nz_ogrid,mfarray_ogrid) :: f_Hloy,f_Hupy
      
      intent(in) :: f_og
      intent(inout) :: f_Hloy,f_Hupy
!
      ivar1=1; ivar2=min(mcom,size(f_og,4))
!
!  Boundary conditions in y
!  Periodic, with y being the theta direction for the cylindrical grid
!
      do j=ivar1,ivar2
!  Bottom boundary
        f_Hloy(l1_ogrid:l2_ogrid,:,:,j) = f_og(l1_ogrid:l2_ogrid,m2i_ogrid-Hsize:m2i_ogrid-1,n1_ogrid:n2_ogrid,j)
!  Top boundary
        f_Hupy(l1_ogrid:l2_ogrid,:,:,j) = f_og(l1_ogrid:l2_ogrid,m1i_ogrid+1:m1i_ogrid+Hsize,n1_ogrid:n2_ogrid,j)
      enddo

      if (lchemistry .and. linterp_pressure) then
        f_Hloy(l1_ogrid:l2_ogrid,:,:,ipp) = f_og(l1_ogrid:l2_ogrid,m2i_ogrid-Hsize:m2i_ogrid-1,n1_ogrid:n2_ogrid,ipp)
        f_Hupy(l1_ogrid:l2_ogrid,:,:,ipp) = f_og(l1_ogrid:l2_ogrid,m1i_ogrid+1:m1i_ogrid+Hsize,n1_ogrid:n2_ogrid,ipp)
        if (lchemistry) then
          f_Hloy(l1_ogrid:l2_ogrid,:,:,iRR) = f_og(l1_ogrid:l2_ogrid,m2i_ogrid-Hsize:m2i_ogrid-1,n1_ogrid:n2_ogrid,iRR)
          f_Hupy(l1_ogrid:l2_ogrid,:,:,iRR) = f_og(l1_ogrid:l2_ogrid,m1i_ogrid+1:m1i_ogrid+Hsize,n1_ogrid:n2_ogrid,iRR)
        endif
      endif

    endsubroutine boundconds_y_filter
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
    subroutine find_proc_cartesian(xyz,from_proc,incl_gp)
!
!  Find the processor that stores the grid points, and return the processor
!  id and the grid index of the bottom neighbouring point on a global grid.
!  Necessary for interpolation between grids on parallel systems
!
!  13-apr-17/Jorgen: Coded
!  19-apr-18/Jorgen: Added possibility to use ghost points
!
      real, dimension(3), intent(in) :: xyz
      integer, intent(out) :: from_proc
      integer :: i,ishift
      logical :: found_proc
      logical, optional :: incl_gp
      real, parameter :: fDP=1.e-15
!
      found_proc=.false. 
!
!  If ghost points are included, check first current processor for the point
!
      if (present(incl_gp)) then
        if(interpolation_method==1) then
          ishift=0
        else
          ishift=1
        endif
        if(incl_gp) then
          if( ((xyz(1)>(x(1+ishift))).and.(xyz(1)<(x(mx-ishift)))) .and. &
              ((xyz(2)>(y(1+ishift))).and.(xyz(2)<(y(my-ishift)))) .and. &
              ((xyz(3)>(z(1+ishift))).and.(xyz(3)<(z(mz-ishift)))) )  then
              from_proc=iproc
              found_proc=.true.
              return
          endif
        endif
      endif
      do i=1,ncpus
        if( ((xyz(1)>=(xyz0_loc_all(i,1)-fDP)).and.(xyz(1)<=(xyz1_loc_all(i,1)+fDP))) .and. &
            ((xyz(2)>=(xyz0_loc_all(i,2)-fDP)).and.(xyz(2)<=(xyz1_loc_all(i,2)+fDP))) .and. &
            ((xyz(3)>=(xyz0_loc_all(i,3)-fDP)).and.(xyz(3)<=(xyz1_loc_all(i,3)+fDP))) )  then
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
!  and for particle computations
!
!  13-apr-17/Jorgen: Coded
!
      real, dimension(3), intent(in) :: rthz
      integer, intent(out) :: from_proc
      integer :: i
      logical :: found_proc
      real, parameter :: fDP=1.e-15
!
      found_proc=.false.
      do i=1,ncpus
        if( ((rthz(1)>=(xyz0_loc_all_ogrid(i,1)-fDP)).and.(rthz(1)<=(xyz1_loc_all_ogrid(i,1)+fDP))) .and. &
            ((rthz(2)>=(xyz0_loc_all_ogrid(i,2)-fDP)).and.(rthz(2)<=(xyz1_loc_all_ogrid(i,2)+fDP))) .and. &
            ((rthz(3)>=(xyz0_loc_all_ogrid(i,3)-fDP)).and.(rthz(3)<=(xyz1_loc_all_ogrid(i,3)+fDP))) )  then
            from_proc=i-1
            found_proc=.true.
            exit
        endif
      enddo
      if(.not.found_proc) then
        call fatal_error('find_proc_curvilinear', &
          'could not locate interpolation point on any processor!')
      endif
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
    subroutine find_near_ind_local_cart(inear_loc,xxp,lcheck)
!
!  Find nearest local indices of point xxp
!  Return only indices correponding to low corner of cube containing the point xxp,
!  i.e, inear_glob = (/ ix0, iy0, iz0 /), where x(ix0) <= xxp(1) <= x(ix0+1), etc. 
!
!  01-aug-17/Jorgen: Coded
!
      integer, dimension(3), intent(out) :: inear_loc
      real, dimension(3), intent(in) :: xxp
      logical, intent(in) :: lcheck
!   
      integer, dimension(3) :: inear_glob
!
      call find_near_ind_global_cart(inear_glob,xxp,lcheck)
      call ind_global_to_local_cart(inear_glob,inear_loc,lcheck)
!
    endsubroutine find_near_ind_local_cart
!***********************************************************************
    subroutine find_near_ind_global_cart(ineargrid,xxp,lcheck)
!
!  Find nearest global indices of point xxp
!  Return only indices correponding to low corner of cube containing the point xxp,
!  i.e, inear_glob = (/ ix0, iy0, iz0 /), where xglobal(ix0) <= xxp(1) <= xglobal(ix0+1), etc. 
!
!  Find nearest grid point by bisection if the grid is not equidistant.
!
!  01-aug-17/Jorgen: Coded
!
      integer, dimension(3), intent(out) :: ineargrid
      real, dimension(3), intent(in) :: xxp
      logical, intent(in) :: lcheck

      integer :: ix0, iy0, iz0
      real, save :: dx1, dy1, dz1
      logical, save :: lfirstcall=.true.

!
!  Default values in case of missing directions.
!
      ix0=nghost+1; iy0=nghost+1; iz0=nghost+1
      if(lfirstcall) then
        dx1=dx_1(l1) 
        dy1=dy_1(m1) 
        dz1=dz_1(n1)
        lfirstcall=.false.
      endif
!
!  Find nearest grid point in x-direction.
!
      if (nxgrid/=1) then
        if (lequidist(1)) then
          ix0 = floor((xxp(1)-xglobal(1))*dx1) + 1
        else
          call find_low_gp_index_bisection(xxp(1),xglobal,ix0)
        endif
      endif
!
!  Find nearest grid point in y-direction.
!
      if (nygrid/=1) then
        if (lequidist(2)) then
          iy0 = floor((xxp(2)-yglobal(1))*dy1) + 1
        else
          call find_low_gp_index_bisection(xxp(2),yglobal,iy0)
        endif
      endif
!
!  Find nearest grid point in z-direction.
!
      if (nzgrid/=1) then
        if (lequidist(3)) then
          iz0 = floor((xxp(3)-zglobal(1))*dz1) + 1
        else
          call find_low_gp_index_bisection(xxp(3),zglobal,iz0)
        endif
      endif
!
      ineargrid=(/ ix0,iy0,iz0 /)
!
!  If requested, check if the correct grid points are found
!
      !if(lcheck) then
      !  if ((xglobal(ineargrid(1))-xxp(1)  )>1.e-14 .or. &
      !      (xxp(1)-xglobal(ineargrid(1)+1))>1.e-14 .or. &
      !      (yglobal(ineargrid(2))-xxp(2)  )>1.e-14 .or. & 
      !      (xxp(2)-yglobal(ineargrid(2)+1))>1.e-14 .or. & 
      !      (zglobal(ineargrid(3))-xxp(3)  )>1.e-14 .or. & 
      !      (xxp(3)-zglobal(ineargrid(3)+1))>1.e-14) then
      if(lcheck) then
        if ((xglobal(ineargrid(1))-xxp(1)  )>0. .or. &
            (xxp(1)-xglobal(ineargrid(1)+1))>0. .or. &
            (yglobal(ineargrid(2))-xxp(2)  )>0. .or. &
            (xxp(2)-yglobal(ineargrid(2)+1))>0. .or. &
            (zglobal(ineargrid(3))-xxp(3)  )>0. .or. &
            (xxp(3)-zglobal(ineargrid(3)+1))>0.) then 
!
!  Try adjusting. Might be needed if there is a very close overlap between grid points
!
          if ((xglobal(ineargrid(1))-xxp(1)  )>0. .and. & 
              (xglobal(ineargrid(1))-xxp(1)  )<1.e-12) then
            ineargrid(1) = ineargrid(1) - 1
          elseif ((xxp(1)-xglobal(ineargrid(1)+1))>0. .and. &
              (xxp(1)-xglobal(ineargrid(1)+1))<1.e-12) then
            ineargrid(1) = ineargrid(1) + 1
          endif
          if ((yglobal(ineargrid(2))-xxp(2)  )>0. .and. & 
              (yglobal(ineargrid(2))-xxp(2)  )<1.e-12) then
            ineargrid(2) = ineargrid(2) - 1
          elseif ((xxp(2)-yglobal(ineargrid(2)+1))>0. .and. &
              (xxp(2)-yglobal(ineargrid(2)+1))<1.e-12) then
            ineargrid(2) = ineargrid(2) + 1
          endif
          if ((zglobal(ineargrid(3))-xxp(3)  )>0. .and. & 
              (zglobal(ineargrid(3))-xxp(3)  )<1.e-12) then
            ineargrid(3) = ineargrid(3) - 1
          elseif ((xxp(3)-zglobal(ineargrid(3)+1))>0. .and. &
              (xxp(3)-zglobal(ineargrid(3)+1))<1.e-12) then
            ineargrid(3) = ineargrid(3) + 1
          endif
!             
!  If there is still a problem, return an error message
!
          if ((xglobal(ineargrid(1))-xxp(1)  )>0. .or. & 
              (xxp(1)-xglobal(ineargrid(1)+1))>0. .or. & 
              (yglobal(ineargrid(2))-xxp(2)  )>0. .or. & 
              (xxp(2)-yglobal(ineargrid(2)+1))>0. .or. & 
              (zglobal(ineargrid(3))-xxp(3)  )>0. .or. & 
              (xxp(3)-zglobal(ineargrid(3)+1))>0.) then  
            print*, 'Information about what went wrong:'
            print*, '----------------------------------'
            print*, 'ERROR: find nearest grid point, cartesian'
            print*, 'Information about what went wrong:'
            print*, '----------------------------------'
            print*, 'it, itsub, t=', it, itsub, t
            print*, 'iproc  =', iproc
            print*, 'xxp    =', xxp
            print*, 'ineargrid   =', ineargrid(:)
            print*, 'xglobal(ix0),xglobal(ix0+1)',xglobal(ineargrid(1)),xglobal(ineargrid(1)+1)
            print*, 'yglobal(iy0),yglobal(iy0+1)',yglobal(ineargrid(2)),yglobal(ineargrid(2)+1)
            print*, 'zglobal(iz0),zglobal(iz0+1)',zglobal(ineargrid(3)),zglobal(ineargrid(3)+1)
            call fatal_error_local('find_near_grid_point_cartesian','')
          endif
        endif
      endif
    endsubroutine find_near_ind_global_cart
!***********************************************************************
    subroutine find_near_ind_local_curv(inear_loc,xxp,lcheck)
!
!  Find nearest local indices of point xxp on curvilinear grid
!  Return only indices correponding to low corner of cube containing the point xxp,
!  i.e, inear_glob = (/ ix0, iy0, iz0 /), where x_ogrid(ix0) <= xxp(1) <= x_ogrid(ix0+1), etc. 
!
!  01-aug-17/Jorgen: Coded
!
      integer, dimension(3), intent(out) :: inear_loc
      real, dimension(3), intent(in) :: xxp
      logical, intent(in) :: lcheck
!   
      integer, dimension(3) :: inear_glob
!
      call find_near_ind_global_curv(inear_glob,xxp,lcheck)
      call ind_global_to_local_curv(inear_glob,inear_loc,lcheck)
!
    endsubroutine find_near_ind_local_curv
!***********************************************************************
    subroutine find_near_ind_global_curv(ineargrid,xxp,lcheck)
!
!  Find nearest global indices of point xxp on the curvilinear grid
!  Return only indices correponding to low corner of cube containing the point xxp,
!  i.e, inear_glob = (/ ix0, iy0, iz0 /), where xglobal_ogrid(ix0) <= xxp(1) <= xglobal_ogrid(ix0+1), etc. 
!
!  Find nearest grid point by bisection if the grid is not equidistant.
!
!  01-aug-17/Jorgen: Coded
!
      integer, dimension(3), intent(out) :: ineargrid
      real, dimension(3), intent(in) :: xxp
      logical, intent(in) :: lcheck

      integer :: ix0, iy0, iz0
      real, save :: dx1_ogrid, dy1_ogrid, dz1_ogrid
      logical, save :: lfirstcall=.true.
!
!  Default values in case of missing directions.
!
      ix0=nghost+1; iy0=nghost+1; iz0=nghost+1
!
      if(lfirstcall) then
        dx1_ogrid=dx_1_ogrid(l1_ogrid) 
        dy1_ogrid=dy_1_ogrid(m1_ogrid) 
        dz1_ogrid=dz_1_ogrid(n1_ogrid)
        lfirstcall=.false.
      endif
!
!  Find nearest grid point in x-direction.
!
      if (nxgrid_ogrid/=1) then
        if (lequidist_ogrid(1)) then
          ix0 = floor((xxp(1)-xglobal_ogrid(1))*dx1_ogrid) + 1
        else
          call find_low_gp_index_bisection(xxp(1),xglobal_ogrid,ix0)
        endif
      endif
!
!  Find nearest grid point in y-direction.
!
      if (nygrid_ogrid/=1) then
        if (lequidist_ogrid(2)) then
          iy0 = floor((xxp(2)-yglobal_ogrid(1))*dy1_ogrid) + 1
        else
          call find_low_gp_index_bisection(xxp(2),yglobal_ogrid,iy0)
        endif
      endif
!
!  Find nearest grid point in z-direction.
!
      if (nzgrid_ogrid/=1) then
        if (lequidist_ogrid(3)) then
          iz0 = floor((xxp(3)-zglobal_ogrid(1))*dz1_ogrid) + 1
        else
          call find_low_gp_index_bisection(xxp(3),zglobal_ogrid,iz0)
        endif
      endif
!
      ineargrid=(/ ix0,iy0,iz0 /)
!
!  If requested, check if the correct grid points are found
!
      if(lcheck) then
        if ((xglobal_ogrid(ineargrid(1))-xxp(1)  )>0. .or. & !1.e-14 .or. &  
            (xxp(1)-xglobal_ogrid(ineargrid(1)+1))>0. .or. & !1.e-14 .or. &  
            (yglobal_ogrid(ineargrid(2))-xxp(2)  )>0. .or. & !1.e-14 .or. &  
            (xxp(2)-yglobal_ogrid(ineargrid(2)+1))>0. .or. & !1.e-14 .or. &  
            (zglobal_ogrid(ineargrid(3))-xxp(3)  )>0. .or. & !1.e-14 .or. &  
            (xxp(3)-zglobal_ogrid(ineargrid(3)+1))>0.) then  !1.e-14) then
!
!  Try adjusting. Might be needed if there is a very close overlap between grid points
!
          if ((xglobal_ogrid(ineargrid(1))-xxp(1)  )>0. .and. & 
              (xglobal_ogrid(ineargrid(1))-xxp(1)  )<1.e-12) then
            ineargrid(1) = ineargrid(1) - 1
          elseif ((xxp(1)-xglobal_ogrid(ineargrid(1)+1))>0. .and. &
              (xxp(1)-xglobal_ogrid(ineargrid(1)+1))<1.e-12) then
            ineargrid(1) = ineargrid(1) + 1
          endif
          if ((yglobal_ogrid(ineargrid(2))-xxp(2)  )>0. .and. & 
              (yglobal_ogrid(ineargrid(2))-xxp(2)  )<1.e-12) then
            ineargrid(2) = ineargrid(2) - 1
          elseif ((xxp(2)-yglobal_ogrid(ineargrid(2)+1))>0. .and. &
              (xxp(2)-yglobal_ogrid(ineargrid(2)+1))<1.e-12) then
            ineargrid(2) = ineargrid(2) + 1
          endif
          if ((zglobal_ogrid(ineargrid(3))-xxp(3)  )>0. .and. & 
              (zglobal_ogrid(ineargrid(3))-xxp(3)  )<1.e-12) then
            ineargrid(3) = ineargrid(3) - 1
          elseif ((xxp(3)-zglobal_ogrid(ineargrid(3)+1))>0. .and. &
              (xxp(3)-zglobal_ogrid(ineargrid(3)+1))<1.e-12) then
            ineargrid(3) = ineargrid(3) + 1
          endif
!             
!  If there is still a problem, return an error message
!
          if ((xglobal_ogrid(ineargrid(1))-xxp(1)  )>0. .or. & 
              (xxp(1)-xglobal_ogrid(ineargrid(1)+1))>0. .or. & 
              (yglobal_ogrid(ineargrid(2))-xxp(2)  )>0. .or. & 
              (xxp(2)-yglobal_ogrid(ineargrid(2)+1))>0. .or. & 
              (zglobal_ogrid(ineargrid(3))-xxp(3)  )>0. .or. & 
              (xxp(3)-zglobal_ogrid(ineargrid(3)+1))>0.) then  
!
            print*, 'Information about what went wrong:'
            print*, '----------------------------------'
            print*, 'ERROR: find nearest grid point, curvilinear'
            print*, 'Information about what went wrong:'
            print*, '----------------------------------'
            print*, 'it, itsub, t=', it, itsub, t
            print*, 'iproc  =', iproc
            print*, 'xxp    =', xxp
            print*, 'ineargrid   =', ineargrid(:)
            print*, 'xglobal(ix0),xglobal(ix0+1)',xglobal_ogrid(ineargrid(1)),xglobal_ogrid(ineargrid(1)+1)
            print*, 'yglobal(iy0),yglobal(iy0+1)',yglobal_ogrid(ineargrid(2)),yglobal_ogrid(ineargrid(2)+1)
            print*, 'zglobal(iz0),zglobal(iz0+1)',zglobal_ogrid(ineargrid(3)),zglobal_ogrid(ineargrid(3)+1)
            call fatal_error_local('find_near_ind_global_curvilinear','')
            call fatal_error('find_near_ind_global_curv','nearest grid point error')
          endif
        endif
      endif
!
    endsubroutine find_near_ind_global_curv
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
          call update_ghosts_ogrid(f_ogrid)
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
        call update_ghosts_ogrid(f_ogrid)
        call safe_character_assign(file,trim(chsnap))
        call output_snap_ogrid(f_ogrid,file=file)
        if (present(flist)) call log_filename_to_file(file,flist)
      endif
!
      if (lformat) call output_snap_form_ogrid (file)
!
    endsubroutine wsnap_ogrid
!***********************************************************************
    subroutine update_ghosts_ogrid(f_og)
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
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
!
      !call boundconds_x_ogrid
      call initiate_isendrcv_bdry_ogrid(f_og)
      call finalize_isendrcv_bdry_ogrid(f_og)
      if (nprocy==1)                  call boundconds_y_ogrid(f_og)
      if ((nprocz==1).and.(nzgrid>1)) call boundconds_z_ogrid(f_og)
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
      if (lstore_ogTT) then
         call input_snap_ogrid(chsnap,f_ogrid(:,:,:,1:mfarray_ogrid-mogaux),mfarray_ogrid-mogaux,mode)
      else
         call input_snap_ogrid(chsnap,f_ogrid,mfarray_ogrid,mode)
      endif
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
    subroutine read_snap_single_ogrid(file,a,x_ogrid,y_ogrid,z_ogrid, &
          dx_ogrid,dy_ogrid,dz_ogrid,deltay_ogrid,nv,mode)
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
      integer, intent(in) :: nv, mode
      real(KIND=rkind4), dimension (mx_ogrid,my_ogrid,mz_ogrid,nv), intent(out) :: a
!
      real(KIND=rkind4) :: t_sp, t_sgl

      real(KIND=rkind4),                       intent(out) :: dx_ogrid, dy_ogrid, dz_ogrid, deltay_ogrid
      real(KIND=rkind4), dimension (mx_ogrid), intent(out) :: x_ogrid
      real(KIND=rkind4), dimension (my_ogrid), intent(out) :: y_ogrid
      real(KIND=rkind4), dimension (mz_ogrid), intent(out) :: z_ogrid

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
    subroutine read_snap_double_ogrid(file,a,x_ogrid,y_ogrid,z_ogrid, &
          dx_ogrid,dy_ogrid,dz_ogrid,deltay_ogrid,nv,mode)
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
      real(KIND=rkind8), dimension (mx_ogrid,my_ogrid,mz_ogrid,nv), intent(out) :: a
!
      real(KIND=rkind8) :: t_sp, t_dbl

      real(KIND=rkind8), intent(out) :: dx_ogrid, dy_ogrid, dz_ogrid, deltay_ogrid
      real(KIND=rkind8), dimension (mx_ogrid), intent(out) :: x_ogrid
      real(KIND=rkind8), dimension (my_ogrid), intent(out) :: y_ogrid
      real(KIND=rkind8), dimension (mz_ogrid), intent(out) :: z_ogrid

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
    subroutine map_nearest_grid_ogrid(xxp,ineargridproc_ogrid,rthz)
!
!  Find index (ix0, iy0, iz0) of nearest grid point of particle in global
!  coordinates.
!
!  jul/aug-17/Jorgen: Adapted from map_nearest_grid in particles_sub
!                     to work for overlapping grid
!
      real, dimension (3) :: xxp, rthz
      integer, dimension (4) :: ineargridproc_ogrid
!
      intent(in)  :: xxp
      intent(out) :: ineargridproc_ogrid,rthz
!
      call get_polar_coords(xxp,rthz)
      call find_near_ind_global_curv(ineargridproc_ogrid(1:3),rthz,lcheck_interpolation)
      call find_proc_curvilinear(rthz,ineargridproc_ogrid(4))
!
    endsubroutine map_nearest_grid_ogrid
!***********************************************************************
    subroutine find_low_gp_index_bisection(qp,q,iq0)
!
!  Find nearest local indices of point qp on grid q by bisecting the interval.
!  Its main use is for non-equidistant grids. 
!  Return only indices correponding lower point on line  containing qp
!  i.e, q(iq0) <= qp <= q(iq0+1), etc. 
!
!  01-aug-17/Jorgen: Adapted from particles_mpicomm
!
      real, dimension (:) :: q
      real :: qp
      integer :: iq0,jl,ju,jm
!
      intent (in) :: qp,q
      intent (out) :: iq0
!
      jl=1+nghost
      ju=size(q)-nghost
!
      do while((ju-jl)>1)
        jm=(ju+jl)/2
        if (qp > q(jm)) then
          jl=jm
        else
          ju=jm
        endif
      enddo
      iq0=jl
!
    endsubroutine find_low_gp_index_bisection
!***********************************************************************
    subroutine initialize_particles_ogrid(ivar1,ivar2)
!
!  Set up f_ogrid_procs, that sends information about ogrid to the appropriate
!  processor that needs this for computation of particle properties (velocity,
!  temperature, etc.)
!  
!  07-jul-17/Jorgen: Coded
!
      use Mpicomm, only: mpirecv_logical, mpisend_logical, mpibcast_logical
      use Viscosity, only: getnu
!
      integer, intent(in) :: ivar1,ivar2
      integer :: i,j,k,iip
      real :: rr
      logical, dimension(ncpus) :: linside_proc = .false.
      real, dimension(3) :: rthz
      integer :: from_proc
      integer :: procs_needed 
      logical, dimension(ncpus,ncpus) :: part_data_comm_glob
      real :: nu
! 
!  Check if any points on the cartesian grid on this processor is inside the part of 
!  the ogrid that is used in computing the flow 
!
      do i=l1,l2
        do j=m1,m2
          rr = radius_ogrid(x(i),y(j))
          if(rr<=r_int_outer .and. rr>=cylinder_radius) then
            do k=n1,n2
              call get_polar_coords(x(i),y(j),z(k),rthz)
              call find_proc_curvilinear(rthz,from_proc)
              linside_proc(from_proc+1) = .true.
            enddo
          endif
        enddo
      enddo
! 
!  Initialize arrays of points needed for particle properties 
!  The ip_proc is needed to transform global coordinates to coordinates local to 
!  the processor considered in the f_ogrid_procs_arrary
!  Procs_pointer is used to point to correct place in ip_proc and f_ogrid_procs array
!  when using point from specific processor
!
!  JONAS: The mogaux in the allocation of f_ogrid_procs is cruuuude, need to implement
!  it more elegantly      
!
      procs_needed = count(linside_proc)
      if(procs_needed>0) then
        allocate(f_ogrid_procs(procs_needed,mx_ogrid,my_ogrid,mz_ogrid,ivar2-ivar1+1+mogaux))
        allocate(ip_proc(procs_needed,3))
        allocate(recv_part_data_from(procs_needed))
        allocate(ip_proc_pointer(ncpus))
        ip_proc_pointer = -1
        k=1
        do iip=0,ncpus-1
          if(linside_proc(iip+1)) then
            recv_part_data_from(k)=iip
            if (lprocz_slowest) then
              ip_proc(k,1) = modulo(iip, nprocx)
              ip_proc(k,2) = modulo(iip/nprocx, nprocy)
              ip_proc(k,3) = iip/nprocxy
            else
              ip_proc(k,1) = modulo(iip, nprocx)
              ip_proc(k,2) = iip/nprocxz
              ip_proc(k,3) = modulo(iip/nprocx, nprocz)
            endif
            ip_proc_pointer(iip+1) = k
            k=k+1
          endif
        enddo
      endif
!
!  Communicate this to the other processors
!
      if(lroot) then
        do iip=0,ncpus-1
          if(iip/=root) then
            call mpirecv_logical(part_data_comm_glob(iip+1,:),ncpus,iip,899)
          else
            part_data_comm_glob(root+1,:)=linside_proc
          endif
        enddo
      else
        call mpisend_logical(linside_proc,ncpus,root,899)
      endif
      call mpibcast_logical(part_data_comm_glob,(/ncpus,ncpus/))
!
!  Determine if this processor must send flow data to any others
!  Never send data to onself
!
      n_procs_send_part_data = count(part_data_comm_glob(:,iproc+1))
      if(part_data_comm_glob(iproc+1,iproc+1)) then
        n_procs_send_part_data = n_procs_send_part_data-1
      endif
      if(n_procs_send_part_data>0) then
        allocate(send_part_data_to(n_procs_send_part_data))
        k=1
        do iip=0,ncpus-1
          if(part_data_comm_glob(iip+1,iproc+1)) then
            if(iip/=iproc) then
              send_part_data_to(k)=iip
              k=k+1
            endif
          endif
        enddo
      endif
      n_procs_recv_part_data = procs_needed
!
!  Send/recv information of this timestep (needed if continuing run from previous simulation)
!
      call update_ogrid_flow_info(ivar1,ivar2)
! 
!  Set momentum thickness, may be needed for particle interpolation
!  Thickness from Weber et al. 2013
!
      call getnu(nu_input=nu)
      delta_momentum = ((0.20669/sqrt(2.))/sqrt(2.*cylinder_radius*init_uu/nu))*(2*cylinder_radius)+cylinder_radius
!  Use only one type of special handling near the surface
      if(lspecial_rad_int_mom) lspecial_rad_int=.false.
!
    endsubroutine initialize_particles_ogrid
!***********************************************************************
    subroutine update_ogrid_flow_info(ivar1,ivar2)
!
!  Communicate f_ogrid to the processors that need it for computation of 
!  particle properties. 
!  
!  07-jul-17/Jorgen: Coded
!
      use Mpicomm, only: mpirecv_real, mpisend_nonblock_real, mpibarrier, mpiwait
!
      integer, intent(in) :: ivar1, ivar2
      integer :: iter, ivar
      integer :: recv_from, send_to
      integer, dimension(3) :: flow_buf_size = (/mx_ogrid,my_ogrid,mz_ogrid/)
      integer, dimension(n_procs_send_part_data,ivar2-ivar1+1) :: ireq2D
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,ivar2-ivar1+1) :: fbufi
!
!  Post non-blocking sends to send flow data to prosessors that need it
!  Note that send_to is never iproc, since this item is taken out of send_par_data_to
!  during initialization
!
      do iter=1,n_procs_send_part_data
        send_to = send_part_data_to(iter)
        do ivar=ivar1,ivar2
          call mpisend_nonblock_real(f_ogrid(:,:,:,ivar),flow_buf_size,send_to,800+ivar,ireq2D(iter,ivar))
        enddo
      enddo
!
!  Post blocking recieves, if this processor recieves any data from other than self
!  Use buffer array to avoid creating temporary array at runtime
!
      do iter=1,n_procs_recv_part_data
        recv_from = recv_part_data_from(iter)
        if(recv_from /= iproc) then
          do ivar=ivar1,ivar2
            call mpirecv_real(fbufi(:,:,:,ivar),flow_buf_size,recv_from,800+ivar)
          enddo
          f_ogrid_procs(iter,:,:,:,ivar1:ivar2) = fbufi
        else
          f_ogrid_procs(iter,:,:,:,ivar1:ivar2) = f_ogrid(:,:,:,ivar1:ivar2)
        endif
      enddo
!
      do iter=1,n_procs_send_part_data
        do ivar=ivar1,ivar2
          call mpiwait(ireq2D(iter,ivar))
        enddo
      enddo
      call mpibarrier
!
    endsubroutine update_ogrid_flow_info
!***********************************************************************
    subroutine interpolate_particles_ogrid(ivar1,ivar2,xxp,gp,inear_glob)
!
!  Use information from the overlapping curvilinear grid to interpolate the velocity
!  of the particle. Must transform to cartesian velocity after interpolation.
!
!  Use linear interpolation for all flow variables, except the 
!  radial component of the velocity. Quadratic interpolation is used for this
!  components, as this is consistent with how the radial velocity behaves near
!  the cylinder.
!
!  31-okt-17/Jorgen: Coded
!
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (4) :: inear_glob
      real :: tmp
!
      intent(in)  :: ivar1, ivar2, xxp, inear_glob
      intent(out) :: gp
!
      if(ivar1==iux) then
         if(ivar2/=iuz) call fatal_error('interpolate_particels_ogrid','should not interpolate rho here!')
         if(particle_interpolate==1) then
            call interpolate_linear_ogrid(ivar1,ivar1+1,xxp,gp(:),inear_glob)
         elseif(particle_interpolate==2) then
            call interpolate_pseudoquad(ivar1,xxp,gp(iux),inear_glob)
            if(lparticle_uradonly) then
               call interpolate_linear_ogrid(ivar1+1,ivar1+1,xxp,gp(iux+1),inear_glob)
            else
               call interpolate_pseudoquad(ivar1+1,xxp,gp(iux+1),inear_glob)
            endif
         elseif(particle_interpolate==3) then
            call interpolate_quad_ogrid(ivar1,xxp,gp(iux),inear_glob)
            if(lparticle_uradonly) then
               call interpolate_linear_ogrid(ivar1+1,ivar1+1,xxp,gp(iux+1),inear_glob)
            else
               call interpolate_quad_ogrid(ivar1+1,xxp,gp(iux+1),inear_glob)
            endif
         elseif(particle_interpolate==4) then
            call interpolate_pseudocubic(ivar1,xxp,gp(iux),inear_glob)
            call interpolate_pseudocubic(ivar1+1,xxp,gp(iux+1),inear_glob)
         endif
!
! Only update z-velocity if 3D run
!
         if(nzgrid_ogrid/=1) then
            call interpolate_linear_ogrid(ivar2,ivar2,xxp,gp(iuz),inear_glob)
         else
            gp(iuz)=0.
         endif
!
! Override interpolation scheme if special handling for particles very close to the surface
! is activated
!
         if(lspecial_rad_int) then
            if((xglobal_ogrid(inear_glob(1))==xyz0_ogrid(1))) then
               call interpolate_ogrid_near(iux,iux,xxp,gp(iux),inear_glob)
            endif
         elseif(lspecial_rad_int_mom) then
            if((xglobal_ogrid(inear_glob(1))<delta_momentum)) then
               call interpolate_ogrid_near_mom(iux,iux,xxp,gp(iux),inear_glob)
            endif
         endif
         tmp=gp(iux)
         gp(iux)=tmp*cos(xxp(2))-gp(iuy)*sin(xxp(2))
         gp(iuy)=tmp*sin(xxp(2))+gp(iuy)*cos(xxp(2))
      else
         call interpolate_linear_ogrid(ivar1,ivar2,xxp,gp,inear_glob)
         if(ivar1<irho) then
            call fatal_error('interpolate_particels_ogrid','should not interpolate anything but rho here!')
         endif
      endif
!   
    endsubroutine interpolate_particles_ogrid
!***********************************************************************
  subroutine interpolate_linear_ogrid(ivar1,ivar2,xxp,gp,inear_glob)
!
!  Use information from the overlapping curvilinear grid to interpolate the velocity
!  of the particle. Must transform to cartesian velocity after interpolation.
!  Interpolate the value of f to arbitrary (xp, yp, zp) curvilinear coordinate
!  using the linear interpolation formula.
!
!    g(x,y,z) = A*x*y*z + B*x*y + C*x*z + D*y*z + E*x + F*y + G*z + H .
!
!  The coefficients are determined by the 8 grid points surrounding the
!  interpolation point.
!  Global coordinates are used for the interpolation, such that only velocity
!  information needs to be communicated between the processors before interpolation.
!
!  06-jul-17/Jorgen: Adapted from linear_interpolate_curvilinear in solid_cells_ogrid.f90
!
    integer :: ivar1, ivar2
    real, dimension (3) :: xxp
    real, dimension (ivar2-ivar1+1) :: gp
    integer, dimension (4) :: inear_glob
!
    real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, g5, g6, g7, g8
    real, dimension (ivar2-ivar1+1) :: f0,f1,gp2
    real :: xp0, yp0, zp0
    real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
    integer :: ix0, iy0, iz0, i 
    integer :: ix0_proc, iy0_proc, iz0_proc, proc, ind_proc
!
    intent(in)  :: ivar1, ivar2, xxp, inear_glob
    intent(out) :: gp
!
    real :: dxx0,dxx1,dyy0,dyy1
    ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3); proc=inear_glob(4)
    ind_proc = ip_proc_pointer(proc+1)
    if(ind_proc<1) then
       print*, 'ERROR: Pointing to f_array that does not exist'
       print*, 'This can be due to too many processors in parralelization'
    endif
!
!  Check if the grid point interval is really correct.
!
    if ((xglobal_ogrid(ix0)<=xxp(1) .and. xglobal_ogrid(ix0+1)>=xxp(1) .or. nxgrid_ogrid==1) .and. &
        (yglobal_ogrid(iy0)<=xxp(2) .and. yglobal_ogrid(iy0+1)>=xxp(2) .or. nygrid_ogrid==1) .and. &
        (zglobal_ogrid(iz0)<=xxp(3) .and. zglobal_ogrid(iz0+1)>=xxp(3) .or. nzgrid_ogrid==1)) then
      ! Everything okay
    else
      print*, 'interpolate_linear_ogrid: Global interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
        call fatal_error('interpolate_linear_ogrid','point outside of interval for particle interpolation')
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
      dx1=1./(xglobal_ogrid(ix0+1)-xglobal_ogrid(ix0))
      dy1=1./(yglobal_ogrid(iy0+1)-yglobal_ogrid(iy0))
      if(nzgrid_ogrid/=1) then
        dz1=1./(zglobal_ogrid(iz0+1)-zglobal_ogrid(iz0))
      else 
        dz1=1.
      endif
!
      dxdy1=dx1*dy1; dxdz1=dx1*dz1; dydz1=dy1*dz1
!!      dxdydz1=dx1*dy1*dz1
      dxx0=xxp(1)-xglobal_ogrid(ix0)
      dxx1=xxp(1)-xglobal_ogrid(ix0+1)
      dyy0=xxp(2)-yglobal_ogrid(iy0)
      dyy1=xxp(2)-yglobal_ogrid(iy0+1)
!
!  Transform global coordinates to coordinates local to the f_ogrid_procs array
!
      ix0_proc=ix0-nx_ogrid*ip_proc(ind_proc,1)
      iy0_proc=iy0-ny_ogrid*ip_proc(ind_proc,2)
      iz0_proc=iz0-nz_ogrid*ip_proc(ind_proc,3)
      !TODO TODO : REMOVE THIS
      !if(ivar2>irho) then
        !print*, 'Debug', shape(f_ogrid_procs), 'ivars',ivar1,ivar2,irho,iTT
        !print*, 'ERROR: Variable not existing on f_ogrid requested'
      !endif
!
!  Function values at all corners.
!
      g1=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc  ,ivar1:ivar2)
      g2=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc  ,ivar1:ivar2)
      g3=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc  ,ivar1:ivar2)
      g4=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc  ,ivar1:ivar2)
      !print*, 'ogriddebug', g1,g2,g3,g4

      f0=g1*dxx1*(-dx1)+g2*dxx0*dx1
      f1=g3*dxx1*(-dx1)+g4*dxx0*dx1
      gp=f0*dyy1*(-dy1)+f1*dyy0*dy1
!!       gp2 = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + xp0*yp0*dxdy1*(g1-g2-g3+g4)
!!       if(gp(1)/=gp2(1)) then
!! print*, 'ERROR IN LINEAR INTERPOLATION'
!! print*, 'gp:',gp
!! print*, 'gp2:',gp2
!! call fatal_error('interpolate_linear_ogrid','particle velocity interpolation error')
!! endif
!!      g5=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc+1,ivar1:ivar2)
!!      g6=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc+1,ivar1:ivar2)
!!      g7=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc+1,ivar1:ivar2)
!!      g8=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc+1,ivar1:ivar2)
!
!  Interpolation formula.
!
!!      gp = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) + &
!!          xp0*yp0*dxdy1*(g1-g2-g3+g4) + xp0*zp0*dxdz1*(g1-g2-g5+g6) + &
!!          yp0*zp0*dydz1*(g1-g3-g5+g7) + &
!!          xp0*yp0*zp0*dxdydz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
!
!
!  Simplify if only a 2D-run
!
!
!  Do a reality check on the interpolation scheme.
!
      if (lcheck_interpolation) then
        do i=1,ivar2-ivar1+1
          if ((gp(i)>max(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) .or. &
            (gp(i)<min(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) .or. &
            (gp(i)/=gp(i))) then
            if (gp(i)>max(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
              print*, 'interpolate_linear_ogrid: interpolated value is LARGER than'
              print*, 'interpolate_linear_ogrid: a values at the corner points!'
            elseif (gp(i)<min(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
              print*, 'interpolate_linear_ogrid: interpolated value is smaller than'
              print*, 'interpolate_linear_ogrid: a values at the corner points!'
            elseif (gp(i)/=gp(i)) then
              print*, 'interpolate_linear_ogrid: interpolated value is NaN'
            endif
            print*, 'iproc = ', iproc
            !print*, 'ipar = ', ipar
            print*, 'interpolate_linear_ogrid: xxp=', xxp
            print*, 'interpolate_linear_ogrid: x0, y0, z0=', &
                xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
            print*, 'interpolate_linear_ogrid: i, gp(i)=', i, gp(i)
            print*, 'interpolate_linear_ogrid: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
            call fatal_error('interpolate_linear_ogrid','particle velocity interpolation error')
          endif
        enddo
      endif
!
    endsubroutine interpolate_linear_ogrid
!***********************************************************************
    subroutine interpolate_ogrid_near(ivar1,ivar2,xxp,gp,inear_glob)
!
!  Use information from the overlappint curvilinear grid to interpolate the velocity
!  of the particle. 
!  Special handling for particles very close to the surface. 
!
!  26-feb-18/Jorgen: Coded
!
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (4) :: inear_glob
!
      real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, f0, f1!g5, g6, g7, g8
      real :: xp0, yp0, zp0
      real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
      integer :: ix0, iy0, iz0, i 
      integer :: ix0_proc, iy0_proc, iz0_proc, proc, ind_proc
!
      intent(in)  :: ivar1, ivar2, xxp, inear_glob
      intent(out) :: gp
!
      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3); proc=inear_glob(4)
      ind_proc = ip_proc_pointer(proc+1)
      if(ind_proc<1) then
         print*, 'ERROR: Pointing to f_array that does not exist'
         print*, 'This may be due to too many processors in parallel'
      endif
!
!  Check if the grid point interval is really correct.
!
      if ((xglobal_ogrid(ix0)<=xxp(1) .and. xglobal_ogrid(ix0+1)>=xxp(1) .or. nxgrid_ogrid==1) .and. &
          (yglobal_ogrid(iy0)<=xxp(2) .and. yglobal_ogrid(iy0+1)>=xxp(2) .or. nygrid_ogrid==1) .and. &
          (zglobal_ogrid(iz0)<=xxp(3) .and. zglobal_ogrid(iz0+1)>=xxp(3) .or. nzgrid_ogrid==1)) then
        ! Everything okay
      else
        print*, 'interpolate_linear_ogrid: Global interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
        call fatal_error('interpolate_linear_ogrid','point outside of interval for particle interpolation')
        return
      endif
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!  Set it equal to 0 for dimensions having 1 grid points; this will make sure
!  that the interpolation is bilinear for 2D grids.
!
      xp0=0; yp0=0; zp0=0
      !if (nxgrid_ogrid/=1) xp0=xxp(1)-xyz0_ogrid(1)
      if (nxgrid_ogrid/=1) xp0=xxp(1)-xglobal_ogrid(ix0)
      if (nygrid_ogrid/=1) yp0=xxp(2)-yglobal_ogrid(iy0)
      if (nzgrid_ogrid/=1) zp0=xxp(3)-zglobal_ogrid(iz0)
!
!  Calculate derived grid spacing parameters needed for interpolation.
!
      !dx1=1./(xglobal_ogrid(ix0+1)-xyz0_ogrid(1))
      dx1=1./(xglobal_ogrid(ix0+1)-xglobal_ogrid(ix0))
      dy1=1./(yglobal_ogrid(iy0+1)-yglobal_ogrid(iy0))
      if(nzgrid_ogrid/=1) then
        dz1=1./(zglobal_ogrid(iz0+1)-zglobal_ogrid(iz0))
      else 
        dz1=1.
      endif
!
      dxdy1=dx1*dy1; dxdz1=dx1*dz1; dydz1=dy1*dz1
      dxdydz1=dx1*dy1*dz1
!
!  Transform global coordinates to coordinates local to the f_ogrid_procs array
!
      ix0_proc=ix0-nx_ogrid*ip_proc(ind_proc,1)
      iy0_proc=iy0-ny_ogrid*ip_proc(ind_proc,2)
      iz0_proc=iz0-nz_ogrid*ip_proc(ind_proc,3)
!
!  Function values at all corners.
!
      g1=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc  ,ivar1:ivar2)
      g2=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc  ,ivar1:ivar2)
      g3=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc  ,ivar1:ivar2)
      g4=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc  ,ivar1:ivar2)
      ! g5=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc+1,ivar1:ivar2)
      ! g6=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc+1,ivar1:ivar2)
      ! g7=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc+1,ivar1:ivar2)
      ! g8=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc+1,ivar1:ivar2)
!
!  Interpolation formula.
!  Linear interpolation along theta(y) and z-direction
!
      if(nzgrid_ogrid/=1) then
        call fatal_error('interpolate_ogrid_near','not implemented in 3D')
      endif
      
      f0 = g1 + yp0*dy1*(-g1+g3)
      f1 = g2 + yp0*dy1*(-g2+g4)
      
      if(any(f0/=0.)) then
        call fatal_error('interpolate_ogrid_near','interpolated value should be zero at the surface')
      endif
      
      gp = f1*xp0*xp0*dx1*dx1

      !gp = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) + &
      !    xp0*yp0*dxdy1*(g1-g2-g3+g4) + xp0*zp0*dxdz1*(g1-g2-g5+g6) + &
      !    yp0*zp0*dydz1*(g1-g3-g5+g7) + &
      !    xp0*yp0*zp0*dxdydz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
!
!  Do a reality check on the interpolation scheme.
!
      if (lcheck_interpolation) then
        do i=1,ivar2-ivar1+1
          if ((gp(i)>max(0.,g2(i),g4(i))) .or. &
            (gp(i)<min(0.,g2(i),g4(i))) .or. &
            (gp(i)/=gp(i))) then
            if (gp(i)>max(g1(i),g2(i),g3(i),g4(i))) then
              print*, 'interpolate_linear_ogrid: interpolated value is LARGER than'
              print*, 'interpolate_linear_ogrid: a values at the corner points!'
            elseif (gp(i)<min(g1(i),g2(i),g3(i),g4(i))) then
              print*, 'interpolate_linear_ogrid: interpolated value is smaller than'
              print*, 'interpolate_linear_ogrid: a values at the corner points!'
            elseif (gp(i)/=gp(i)) then
              print*, 'interpolate_linear_ogrid: interpolated value is NaN'
            endif
            print*, 'iproc = ', iproc
            !print*, 'ipar = ', ipar
            print*, 'interpolate_linear_ogrid: xxp=', xxp
            print*, 'interpolate_linear_ogrid: x0, y0, z0=', &
                xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
            print*, 'interpolate_linear_ogrid: i, gp(i)=', i, gp(i)
            print*, 'interpolate_linear_ogrid: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i)
            print*, '------------------'
            call fatal_error('interpolate_linear_ogrid','particle velocity interpolation error')
          endif
        enddo
      endif
!
    endsubroutine interpolate_ogrid_near
!***********************************************************************
    subroutine interpolate_ogrid_near_mom(ivar1,ivar2,xxp,gp,inear_glob)
!
!  Use information from the overlappint curvilinear grid to interpolate the velocity
!  of the particle. Do for all particles within momentum thickness of the
!  surface.
!  Special handling for particles very close to the surface. 
!
!  21-mar-18/Jorgen: Coded
!
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (4) :: inear_glob
!
      real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, f0, f1!g5, g6, g7, g8
      real :: xp0, yp0, zp0
      real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
      integer :: ix0, iy0, iz0, i 
      integer :: ix0_proc, iy0_proc, iz0_proc, proc, ind_proc
!
      intent(in)  :: ivar1, ivar2, xxp, inear_glob
      intent(out) :: gp
!
      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3); proc=inear_glob(4)
      ind_proc = ip_proc_pointer(proc+1)
      if(ind_proc<1) then
         print*, 'ERROR: Pointing to f_array that does not exist'
         print*, 'This can be due to too many processors in parralelization'
      endif
!
!  Check if the grid point interval is really correct.
!
      if ((xglobal_ogrid(ix0)<=xxp(1) .and. xglobal_ogrid(ix0+1)>=xxp(1) .or. nxgrid_ogrid==1) .and. &
          (yglobal_ogrid(iy0)<=xxp(2) .and. yglobal_ogrid(iy0+1)>=xxp(2) .or. nygrid_ogrid==1) .and. &
          (zglobal_ogrid(iz0)<=xxp(3) .and. zglobal_ogrid(iz0+1)>=xxp(3) .or. nzgrid_ogrid==1)) then
        ! Everything okay
      else
        print*, 'interpolate_linear_ogrid: Global interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
        call fatal_error('interpolate_linear_ogrid','point outside of interval for particle interpolation')
        return
      endif
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!  Set it equal to 0 for dimensions having 1 grid points; this will make sure
!  that the interpolation is bilinear for 2D grids.
!
      xp0=0; yp0=0; zp0=0
      if (nxgrid_ogrid/=1) xp0=xxp(1)-xyz0_ogrid(1)
      !if (nxgrid_ogrid/=1) xp0=xxp(1)-xglobal_ogrid(ix0)
      if (nygrid_ogrid/=1) yp0=xxp(2)-yglobal_ogrid(iy0)
      if (nzgrid_ogrid/=1) zp0=xxp(3)-zglobal_ogrid(iz0)
!
!  Calculate derived grid spacing parameters needed for interpolation.
!
      dx1=1./(xglobal_ogrid(ix0+1)-xyz0_ogrid(1))
      !dx1=1./(xglobal_ogrid(ix0+1)-xglobal_ogrid(ix0))
      dy1=1./(yglobal_ogrid(iy0+1)-yglobal_ogrid(iy0))
      if(nzgrid_ogrid/=1) then
        dz1=1./(zglobal_ogrid(iz0+1)-zglobal_ogrid(iz0))
      else 
        dz1=1.
      endif
!
      dxdy1=dx1*dy1; dxdz1=dx1*dz1; dydz1=dy1*dz1
      dxdydz1=dx1*dy1*dz1
!
!  Transform global coordinates to coordinates local to the f_ogrid_procs array
!
      ix0_proc=ix0-nx_ogrid*ip_proc(ind_proc,1)
      iy0_proc=iy0-ny_ogrid*ip_proc(ind_proc,2)
      iz0_proc=iz0-nz_ogrid*ip_proc(ind_proc,3)
!
!  Function values at all corners.
!
      !g1=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc  ,ivar1:ivar2)
      g2=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc  ,ivar1:ivar2)
      !g3=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc  ,ivar1:ivar2)
      g4=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc  ,ivar1:ivar2)
      ! g5=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc+1,ivar1:ivar2)
      ! g6=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc+1,ivar1:ivar2)
      ! g7=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc+1,ivar1:ivar2)
      ! g8=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc+1,ivar1:ivar2)
!
!  Interpolation formula.
!  Linear interpolation along theta(y) and z-direction
!
      if(nzgrid_ogrid/=1) then
        call fatal_error('interpolate_ogrid_near_mom','not implemented in 3D')
      endif
      
      !f0 = g1 + yp0*dy1*(-g1+g3)
      f1 = g2 + yp0*dy1*(-g2+g4)
      
      if(any(f0/=0.)) then
        call fatal_error('interpolate_ogrid_near_mom','interpolated value should be zero at the surface')
      endif
      
      gp = f1*xp0*xp0*dx1*dx1

      !gp = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) + &
      !    xp0*yp0*dxdy1*(g1-g2-g3+g4) + xp0*zp0*dxdz1*(g1-g2-g5+g6) + &
      !    yp0*zp0*dydz1*(g1-g3-g5+g7) + &
      !    xp0*yp0*zp0*dxdydz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
!
!  Do a reality check on the interpolation scheme.
!
      if (lcheck_interpolation) then
        do i=1,ivar2-ivar1+1
          if ((gp(i)>max(0.,g2(i),g4(i))) .or. &
            (gp(i)<min(0.,g2(i),g4(i))) .or. &
            (gp(i)/=gp(i))) then
            if (gp(i)>max(g1(i),g2(i),g3(i),g4(i))) then
              print*, 'interpolate_linear_ogrid: interpolated value is LARGER than'
              print*, 'interpolate_linear_ogrid: a values at the corner points!'
            elseif (gp(i)<min(g1(i),g2(i),g3(i),g4(i))) then
              print*, 'interpolate_linear_ogrid: interpolated value is smaller than'
              print*, 'interpolate_linear_ogrid: a values at the corner points!'
            elseif (gp(i)/=gp(i)) then
              print*, 'interpolate_linear_ogrid: interpolated value is NaN'
            endif
            print*, 'iproc = ', iproc
            !print*, 'ipar = ', ipar
            print*, 'interpolate_linear_ogrid: xxp=', xxp
            print*, 'interpolate_linear_ogrid: x0, y0, z0=', &
                xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
            print*, 'interpolate_linear_ogrid: i, gp(i)=', i, gp(i)
            print*, 'interpolate_linear_ogrid: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i)
            print*, '------------------'
            call fatal_error('interpolate_linear_ogrid','particle velocity interpolation error')
          endif
        enddo
      endif
!
    endsubroutine interpolate_ogrid_near_mom
!***********************************************************************
    subroutine interpolate_pseudoquad(ivar,xxp,gp,inear_glob)
!
!  Use information from the overlappint curvilinear grid to interpolate 
!  flow quantities to the particle, using quadratic inteprolation only
!  for the radial direction.
!
!  31-okt-17/Jorgen: Coded
!
      integer :: ivar
      real, dimension (3) :: xxp
      real :: gp
      integer, dimension (4) :: inear_glob
!
      integer :: ix0, iy0, iz0, ix1, iy1, iz1, ix2
      integer :: ix0_proc, iy0_proc, iz0_proc, proc, ind_proc

      real :: dx10_1, dx21_1, dx20_1, dy_1, dz_1
      real :: dxx0, dxx1, dxx2, dyy0, dyy1, dzz0, dzz1
      real :: g000,g100,g010,g110,g200,g210,g001,g101,g011,g111,g201,g211
      real :: f00, f01, f10, f11, h0, h1

      intent(in)  :: ivar, xxp, inear_glob
      intent(out) :: gp
!
      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3); proc=inear_glob(4)
      ind_proc = ip_proc_pointer(proc+1)
      if(ind_proc<1) then
         print*, 'ERROR: Pointing to f_array that does not exist'
         print*, 'This can be due to too many processors in parralelization'
      endif
!
!  Check if the grid point interval is really correct.
!
      ix1=ix0+1
      iy1=iy0+1
      iz1=iz0+1
      if ((xglobal_ogrid(ix0)<=xxp(1) .and. xglobal_ogrid(ix1)>=xxp(1) .or. nxgrid_ogrid==1) .and. &
          (yglobal_ogrid(iy0)<=xxp(2) .and. yglobal_ogrid(iy1)>=xxp(2) .or. nygrid_ogrid==1) .and. &
          (zglobal_ogrid(iz0)<=xxp(3) .and. zglobal_ogrid(iz1)>=xxp(3) .or. nzgrid_ogrid==1)) then
        ! Everything okay
      else
        print*, 'interpolate_pseudoquad: Global interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
        call fatal_error('interpolate_pseudoquad','point outside of interval for particle interpolation')
        return
      endif
!
!  Redefine closes grid point in radial direciton, if necessary
!
      dxx0=xxp(1)-xglobal_ogrid(ix0)
      dxx1=xxp(1)-xglobal_ogrid(ix1)
      if((abs(dxx0)<abs(dxx1)).and.(xglobal_ogrid(ix0)>xyz0_ogrid(1))) then
        ix1=ix0
        ix0=ix0-1
        dxx0=xxp(1)-xglobal_ogrid(ix0)
        dxx1=xxp(1)-xglobal_ogrid(ix1)
      endif
      ix2=ix0+2

      dx10_1=1./(xglobal_ogrid(ix1)-xglobal_ogrid(ix0))
      dx21_1=1./(xglobal_ogrid(ix2)-xglobal_ogrid(ix1))
      dx20_1=1./(xglobal_ogrid(ix2)-xglobal_ogrid(ix0))
      dy_1=  1./(yglobal_ogrid(iy1)-yglobal_ogrid(iy0))

      dxx2=xxp(1)-xglobal_ogrid(ix2)

      dyy0=xxp(2)-yglobal_ogrid(iy0)
      dyy1=xxp(2)-yglobal_ogrid(iy1)
!
      ix0_proc=ix0-nx_ogrid*ip_proc(ind_proc,1)
      iy0_proc=iy0-ny_ogrid*ip_proc(ind_proc,2)
      iz0_proc=iz0-nz_ogrid*ip_proc(ind_proc,3)

      g000=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc  ,ivar)
      g100=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc  ,ivar)
      g010=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc  ,ivar)
      g110=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc  ,ivar)
      g200=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc  ,iz0_proc  ,ivar)
      g210=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc+1,iz0_proc  ,ivar)
!
!  Simplify if only a 2D-run
!
      if(nzgrid_ogrid==1) then
        f00=g000*dxx1*dxx2*dx10_1*dx20_1+g100*dxx0*dxx2*dx10_1*(-dx21_1)+g200*dxx0*dxx1*dx20_1*dx21_1
        f10=g010*dxx1*dxx2*dx10_1*dx20_1+g110*dxx0*dxx2*dx10_1*(-dx21_1)+g210*dxx0*dxx1*dx20_1*dx21_1
        gp=f00*dyy1*(-dy_1)+f10*dyy0*dy_1
      else
        dzz0=xxp(3)-zglobal_ogrid(iz0)
        dzz1=xxp(3)-zglobal_ogrid(iz1)
        dz_1=1./(zglobal_ogrid(iz1)-zglobal_ogrid(iz0))
!
        g001=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc+1,ivar)
        g101=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc+1,ivar)
        g011=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc+1,ivar)
        g111=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc+1,ivar)
        g201=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc  ,iz0_proc+1,ivar)
        g211=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc+1,iz0_proc+1,ivar)
!
        f00=g000*dxx1*dxx2*dx10_1*dx20_1+g100*dxx0*dxx2*dx10_1*(-dx21_1)+g200*dxx0*dxx1*dx20_1*dx21_1
        f10=g010*dxx1*dxx2*dx10_1*dx20_1+g110*dxx0*dxx2*dx10_1*(-dx21_1)+g210*dxx0*dxx1*dx20_1*dx21_1
        f01=g001*dxx1*dxx2*dx10_1*dx20_1+g101*dxx0*dxx2*dx10_1*(-dx21_1)+g201*dxx0*dxx1*dx20_1*dx21_1
        f11=g011*dxx1*dxx2*dx10_1*dx20_1+g111*dxx0*dxx2*dx10_1*(-dx21_1)+g211*dxx0*dxx1*dx20_1*dx21_1

        h0=f00*dyy1*(-dy_1)+f10*dyy0*dy_1
        h1=f01*dyy1*(-dy_1)+f11*dyy0*dy_1

        gp=h0*dzz1*(-dz_1)+h1*dzz0*dz_1
      endif
!
!  Do a reality check on the interpolation scheme.
!
      if (lcheck_interpolation) then
        if ((gp>max(g000,g100,g010,g110,g200,g210)) .or. &
            (gp<min(g000,g100,g010,g110,g200,g210)) .or. &
            (gp/=gp)) then
          if (gp>max(g000,g100,g010,g110,g200,g210)) then
            print*, 'interpolate_pseudoquad: interpolated value is LARGER than'
            print*, 'interpolate_pseudoquad: a values at the corner points!'
          elseif (gp<min(g000,g100,g010,g110,g200,g210)) then
            print*, 'interpolate_pseudoquad: interpolated value is smaller than'
            print*, 'interpolate_pseudoquad: a values at the corner points!'
          elseif (gp/=gp) then
            print*, 'interpolate_linear_ogrid: interpolated value is NaN'
          endif
          print*, 'iproc = ', iproc
          !print*, 'ipar = ', ipar
          print*, 'interpolate_pseudoquad: xxp=', xxp
          print*, 'interpolate_pseudoquad: x0, y0, z0=', &
              xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
          print*, 'interpolate_pseudoquad: ivar, gp=', ivar, gp
          print*, 'interpolate_pseudoquad: g1...g8=', &
              g000,g100,g010,g110,g200,g210
          print*, '------------------'
          call fatal_error('interpolate_pseudoquad','particle velocity interpolation error')
        endif
      endif
!
    endsubroutine interpolate_pseudoquad
!***********************************************************************
    subroutine interpolate_pseudocubic(ivar,xxp,gp,inear_glob)
!
!  Use information from the overlappint curvilinear grid to interpolate 
!  flow quantities to the particle, using quadratic inteprolation only
!  for the radial direction.
!
!  31-okt-17/Jorgen: Coded
!
      integer :: ivar
      real, dimension (3) :: xxp
      real :: gp
      integer, dimension (4) :: inear_glob
!
      integer :: ix0, iy0, iz0, ix1, iy1, iz1, ix2, ix3
      integer :: ix0_proc, iy0_proc, iz0_proc, proc, ind_proc

      real :: dx10_1, dx21_1, dx20_1, dy_1, dz_1
      real :: dxx0, dxx1, dxx2, dyy0, dyy1, dzz0, dzz1
      real :: g000,g100,g010,g110,g200,g210,g001,g101,g011,g111,g201,g211
      real :: g300,g310,dxx3,dx30_1,dx31_1,dx32_1
      real :: f00, f01, f10, f02,f03, f11, h0, h1

      intent(in)  :: ivar, xxp, inear_glob
      intent(out) :: gp
!
      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3); proc=inear_glob(4)
      ind_proc = ip_proc_pointer(proc+1)
      if(ind_proc<1) then
         print*, 'ERROR: Pointing to f_array that does not exist'
         print*, 'This can be due to too many processors in parralelization'
      endif
!
!  Check if the grid point interval is really correct.
!
      ix1=ix0+1
      iy1=iy0+1
      iz1=iz0+1
      if ((xglobal_ogrid(ix0)<=xxp(1) .and. xglobal_ogrid(ix1)>=xxp(1) .or. nxgrid_ogrid==1) .and. &
          (yglobal_ogrid(iy0)<=xxp(2) .and. yglobal_ogrid(iy1)>=xxp(2) .or. nygrid_ogrid==1) .and. &
          (zglobal_ogrid(iz0)<=xxp(3) .and. zglobal_ogrid(iz1)>=xxp(3) .or. nzgrid_ogrid==1)) then
        ! Everything okay
      else
        print*, 'interpolate_pseudocubic: Global interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
        call fatal_error('interpolate_pseudocubic','point outside of interval for particle interpolation')
        return
      endif
!
!  Redefine closes grid point in radial direciton, if necessary
!
      dxx0=xxp(1)-xglobal_ogrid(ix0)
      dxx1=xxp(1)-xglobal_ogrid(ix1)
      if((abs(dxx0)<abs(dxx1)).and.(xglobal_ogrid(ix0)>xyz0_ogrid(1))) then
        ix1=ix0
        ix0=ix0-1
        dxx0=xxp(1)-xglobal_ogrid(ix0)
        dxx1=xxp(1)-xglobal_ogrid(ix1)
      endif
      ix2=ix0+2
      ix3=ix0+3

      dx10_1=1./(xglobal_ogrid(ix1)-xglobal_ogrid(ix0))
      dx21_1=1./(xglobal_ogrid(ix2)-xglobal_ogrid(ix1))
      dx32_1=1./(xglobal_ogrid(ix3)-xglobal_ogrid(ix2))

      dx20_1=1./(xglobal_ogrid(ix2)-xglobal_ogrid(ix0))
      dx31_1=1./(xglobal_ogrid(ix3)-xglobal_ogrid(ix1))

      dx30_1=1./(xglobal_ogrid(ix3)-xglobal_ogrid(ix0))

      dy_1=  1./(yglobal_ogrid(iy1)-yglobal_ogrid(iy0))

      dxx2=xxp(1)-xglobal_ogrid(ix2)
      dxx3=xxp(1)-xglobal_ogrid(ix3)

      dyy0=xxp(2)-yglobal_ogrid(iy0)
      dyy1=xxp(2)-yglobal_ogrid(iy1)
!
      ix0_proc=ix0-nx_ogrid*ip_proc(ind_proc,1)
      iy0_proc=iy0-ny_ogrid*ip_proc(ind_proc,2)
      iz0_proc=iz0-nz_ogrid*ip_proc(ind_proc,3)

      g000=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc  ,ivar)
      g100=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc  ,ivar)
      g010=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc  ,ivar)
      g110=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc  ,ivar)
      g200=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc  ,iz0_proc  ,ivar)
      g210=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc+1,iz0_proc  ,ivar)
      g300=f_ogrid_procs(ind_proc,ix0_proc+3,iy0_proc  ,iz0_proc  ,ivar)
      g310=f_ogrid_procs(ind_proc,ix0_proc+3,iy0_proc+1,iz0_proc  ,ivar)
!
!  Simplify if only a 2D-run
!
      if(xglobal_ogrid(ix0)<xyz0_ogrid(1)) then
          call fatal_error('interpolate_pseudocubic',&
            'illegal inteprolation point, inside cylinder')
      endif
      if(nzgrid_ogrid==1) then 
        f00=g000*dyy1*(-dy_1)+g010*dyy0*dy_1 !x0
        f01=g100*dyy1*(-dy_1)+g110*dyy0*dy_1 !x1
        f02=g200*dyy1*(-dy_1)+g210*dyy0*dy_1 !x2
        f03=g300*dyy1*(-dy_1)+g310*dyy0*dy_1 !x3

        if(xglobal_ogrid(ix0)==xyz0_ogrid(1)) then
          gp=f00*dxx1*dxx2*(-dx10_1)*(-dx20_1) + &
             f01*dxx0*dxx2*( dx10_1)*(-dx21_1) + &
             f02*dxx0*dxx1*( dx20_1)*( dx21_1)
        else
          gp=f00*dxx1*dxx2*dxx3*(-dx10_1)*(-dx20_1)*(-dx30_1) + &
             f01*dxx0*dxx2*dxx3*( dx10_1)*(-dx21_1)*(-dx31_1) + &
             f02*dxx0*dxx1*dxx3*( dx20_1)*( dx21_1)*(-dx32_1) + &
             f03*dxx0*dxx1*dxx2*( dx30_1)*( dx31_1)*( dx32_1) 
        endif
!         f00=g000*dxx1*dxx2*dx10_1*dx20_1+g100*dxx0*dxx2*dx10_1*(-dx21_1)+g200*dxx0*dxx1*dx20_1*dx21_1
!         f10=g010*dxx1*dxx2*dx10_1*dx20_1+g110*dxx0*dxx2*dx10_1*(-dx21_1)+g210*dxx0*dxx1*dx20_1*dx21_1
!         gp=f00*dyy1*(-dy_1)+f10*dyy0*dy_1
      else
! NOT IMPLEMENTED IN 3D
!         dzz0=xxp(3)-zglobal_ogrid(iz0)
!         dzz1=xxp(3)-zglobal_ogrid(iz1)
!         dz_1=1./(zglobal_ogrid(iz1)-zglobal_ogrid(iz0))
! !
!         g001=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc+1,ivar)
!         g101=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc+1,ivar)
!         g011=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc+1,ivar)
!         g111=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc+1,ivar)
!         g201=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc  ,iz0_proc+1,ivar)
!         g211=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc+1,iz0_proc+1,ivar)
! !
!         f00=g000*dxx1*dxx2*dx10_1*dx20_1+g100*dxx0*dxx2*dx10_1*(-dx21_1)+g200*dxx0*dxx1*dx20_1*dx21_1
!         f10=g010*dxx1*dxx2*dx10_1*dx20_1+g110*dxx0*dxx2*dx10_1*(-dx21_1)+g210*dxx0*dxx1*dx20_1*dx21_1
!         f01=g001*dxx1*dxx2*dx10_1*dx20_1+g101*dxx0*dxx2*dx10_1*(-dx21_1)+g201*dxx0*dxx1*dx20_1*dx21_1
!         f11=g011*dxx1*dxx2*dx10_1*dx20_1+g111*dxx0*dxx2*dx10_1*(-dx21_1)+g211*dxx0*dxx1*dx20_1*dx21_1
! 
!         h0=f00*dyy1*(-dy_1)+f10*dyy0*dy_1
!         h1=f01*dyy1*(-dy_1)+f11*dyy0*dy_1
! 
!         gp=h0*dzz1*(-dz_1)+h1*dzz0*dz_1
      endif
!
!  Do a reality check on the interpolation scheme.
!
      if (lcheck_interpolation) then
        if ((gp>max(g000,g100,g010,g110,g200,g210)) .or. &
            (gp<min(g000,g100,g010,g110,g200,g210)) .or. &
            (gp/=gp)) then
          if (gp>max(g000,g100,g010,g110,g200,g210)) then
            print*, 'interpolate_pseudocubic: interpolated value is LARGER than'
            print*, 'interpolate_pseudocubic: a values at the corner points!'
          elseif (gp<min(g000,g100,g010,g110,g200,g210)) then
            print*, 'interpolate_pseudocubic: interpolated value is smaller than'
            print*, 'interpolate_pseudocubic: a values at the corner points!'
          elseif (gp/=gp) then
            print*, 'interpolate_pseudocubic: interpolated value is NaN'
          endif
          print*, 'iproc = ', iproc
          !print*, 'ipar = ', ipar
          print*, 'interpolate_pseudocubic: xxp=', xxp
          print*, 'interpolate_pseudocubic: x0, y0, z0=', &
              xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
          print*, 'interpolate_pseudocubic: ivar, gp=', ivar, gp
          print*, 'interpolate_pseudocubic: g1...g8=', &
              g000,g100,g010,g110,g200,g210
          print*, '------------------'
          call fatal_error('interpolate_pseudocubic','particle velocity interpolation error')
        endif
      endif
!
    endsubroutine interpolate_pseudocubic
!***********************************************************************
    subroutine interpolate_quad_ogrid(ivar,xxp,gp,inear_glob)
!
!  Use information from the overlappint curvilinear grid to interpolate 
!  flow quantities to the particle, using quadratic inteprolation only
!  for the radial direction.
!
!  31-okt-17/Jorgen: Coded
!
      integer :: ivar
      real, dimension (3) :: xxp
      real :: gp
      integer, dimension (4) :: inear_glob
!
      integer :: ix0, iy0, iz0, ix1, iy1, iz1, ix2
      integer :: ix0_proc, iy0_proc, iz0_proc, proc, ind_proc

      real :: dx10_1, dx21_1, dx20_1, dy10_1, dy21_1, dy20_1!, dz_1
      real :: dxx0, dxx1, dxx2, dyy0, dyy1, dyy2!, dzz0, dzz1
      !real :: g000,g100,g010,g110,g200,g210,g001,g101,g011,g111,g201,g211
      real, dimension(3,3,2) :: gN
      !real :: f00, f01, f10, f11, h0, h1
      real, dimension(3,2) :: fN

      intent(in)  :: ivar, xxp, inear_glob
      intent(out) :: gp
!
      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3); proc=inear_glob(4)
      ind_proc = ip_proc_pointer(proc+1)
      if(ind_proc<1) then
         print*, 'ERROR: Pointing to f_array that does not exist'
         print*, 'This may be due to too many processors in parallel'
      endif
!
!  Check if the grid point interval is really correct.
!
      ix1=ix0+1
      iy1=iy0+1
      iz1=iz0+1
      if ((xglobal_ogrid(ix0)<=xxp(1) .and. xglobal_ogrid(ix1)>=xxp(1) .or. nxgrid_ogrid==1) .and. &
          (yglobal_ogrid(iy0)<=xxp(2) .and. yglobal_ogrid(iy1)>=xxp(2) .or. nygrid_ogrid==1) .and. &
          (zglobal_ogrid(iz0)<=xxp(3) .and. zglobal_ogrid(iz1)>=xxp(3) .or. nzgrid_ogrid==1)) then
        ! Everything okay
      else
        print*, 'interpolate_quad_ogrid: Global interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
        call fatal_error('interpolate_quad_ogrid','point outside of interval for particle interpolation')
        return
      endif
!
!  Redefine closes grid point in radial direciton, if necessary
!
      dxx0=xxp(1)-xglobal_ogrid(ix0)
      dxx1=xxp(1)-xglobal_ogrid(ix1)
      if((abs(dxx0)<abs(dxx1)).and.(xglobal_ogrid(ix0)>xyz0_ogrid(1))) then
        ix1=ix0
        ix0=ix0-1
        dxx0=xxp(1)-xglobal_ogrid(ix0)
        dxx1=xxp(1)-xglobal_ogrid(ix1)
      endif
      ix2=ix0+2
      dxx2=xxp(1)-xglobal_ogrid(ix2)

      dx10_1=1./(xglobal_ogrid(ix1)-xglobal_ogrid(ix0))
      dx21_1=1./(xglobal_ogrid(ix2)-xglobal_ogrid(ix1))
      dx20_1=1./(xglobal_ogrid(ix2)-xglobal_ogrid(ix0))
!
      dyy0=xxp(2)-yglobal_ogrid(iy0)
      dyy1=xxp(2)-yglobal_ogrid(iy1)
      if((abs(dyy0)<abs(dyy1))) then
        iy1=iy0
        iy0=iy0-1
        dyy0=xxp(2)-yglobal_ogrid(iy0)
        dyy1=xxp(2)-yglobal_ogrid(iy1)
      endif
      iy2=iy0+2
      dyy2=xxp(2)-yglobal_ogrid(iy2)

      dy10_1=1./(yglobal_ogrid(iy1)-yglobal_ogrid(iy0))
      dy21_1=1./(yglobal_ogrid(iy2)-yglobal_ogrid(iy1))
      dy20_1=1./(yglobal_ogrid(iy2)-yglobal_ogrid(iy0))
!
      ix0_proc=ix0-nx_ogrid*ip_proc(ind_proc,1)
      iy0_proc=iy0-ny_ogrid*ip_proc(ind_proc,2)
      iz0_proc=iz0-nz_ogrid*ip_proc(ind_proc,3)

      gN(1,1,1)=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc  ,ivar)
      gN(2,1,1)=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc  ,ivar)
      gN(3,1,1)=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc  ,iz0_proc  ,ivar)
      gN(1,2,1)=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc  ,ivar)
      gN(2,2,1)=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc  ,ivar)
      gN(3,2,1)=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc+1,iz0_proc  ,ivar)
      gN(1,3,1)=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+2,iz0_proc  ,ivar)
      gN(2,3,1)=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+2,iz0_proc  ,ivar)
      gN(3,3,1)=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc+2,iz0_proc  ,ivar)
!
!  Simplify if only a 2D-run
!
      if(nzgrid_ogrid==1) then
        fN(1:3,1)=gN(1,:,1)*dxx1*dxx2*dx10_1*dx20_1+gN(2,:,1)*dxx0*dxx2*dx10_1*(-dx21_1)+gN(3,:,1)*dxx0*dxx1*dx20_1*dx21_1
        gp=fN(1,1)*dyy1*dyy2*dy10_1*dy20_1+fN(2,1)*dyy0*dyy2*dy10_1*(-dy21_1)+fN(3,1)*dyy0*dyy1*dy20_1*dy21_1
        !f10=g(1,2,1)*dxx1*dxx2*dx10_1*dx20_1+g(2,2,1)*dxx0*dxx2*dx10_1*(-dx21_1)+g(3,2,1)*dxx0*dxx1*dx20_1*dx21_1
        !f20=g(1,3,1)*dxx1*dxx2*dx10_1*dx20_1+g(2,3,1)*dxx0*dxx2*dx10_1*(-dx21_1)+g(3,3,1)*dxx0*dxx1*dx20_1*dx21_1
      else
        print*, 'ERROR: 3D not implementet'
        !dzz0=xxp(3)-zglobal_ogrid(iz0)
        !dzz1=xxp(3)-zglobal_ogrid(iz1)
        !dz_1=1./(zglobal_ogrid(iz1)-zglobal_ogrid(iz0))
!
        !g001=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc  ,iz0_proc+1,ivar)
        !g101=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc  ,iz0_proc+1,ivar)
        !g011=f_ogrid_procs(ind_proc,ix0_proc  ,iy0_proc+1,iz0_proc+1,ivar)
        !g111=f_ogrid_procs(ind_proc,ix0_proc+1,iy0_proc+1,iz0_proc+1,ivar)
        !g201=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc  ,iz0_proc+1,ivar)
        !g211=f_ogrid_procs(ind_proc,ix0_proc+2,iy0_proc+1,iz0_proc+1,ivar)
!
        !f00=g000*dxx1*dxx2*dx10_1*dx20_1+g100*dxx0*dxx2*dx10_1*(-dx21_1)+g200*dxx0*dxx1*dx20_1*dx21_1
        !f10=g010*dxx1*dxx2*dx10_1*dx20_1+g110*dxx0*dxx2*dx10_1*(-dx21_1)+g210*dxx0*dxx1*dx20_1*dx21_1
        !f01=g001*dxx1*dxx2*dx10_1*dx20_1+g101*dxx0*dxx2*dx10_1*(-dx21_1)+g201*dxx0*dxx1*dx20_1*dx21_1
        !f11=g011*dxx1*dxx2*dx10_1*dx20_1+g111*dxx0*dxx2*dx10_1*(-dx21_1)+g211*dxx0*dxx1*dx20_1*dx21_1

        !h0=f00*dyy1*(-dy_1)+f10*dyy0*dy_1
        !h1=f01*dyy1*(-dy_1)+f11*dyy0*dy_1

        !gp=h0*dzz1*(-dz_1)+h1*dzz0*dz_1
      endif
!
!  Do a reality check on the interpolation scheme.
!
      if (lcheck_interpolation) then
        if ((gp>maxval(gN(:,:,1))) .or. &
            (gp<minval(gN(:,:,1))) .or. &
            (gp/=gp)) then
          if (gp>maxval(gN(:,:,1))) then
            print*, 'interpolate_quad_ogrid: interpolated value is LARGER than'
            print*, 'interpolate_quad_ogrid: a values at the corner points!'
          elseif (gp<minval(gN(:,:,1))) then
            print*, 'interpolate_quad_ogrid: interpolated value is smaller than'
            print*, 'interpolate_quad_ogrid: a values at the corner points!'
          elseif (gp/=gp) then
            print*, 'interpolate_quad_ogrid: interpolated value is NaN'
          endif
          print*, 'iproc = ', iproc
          !print*, 'ipar = ', ipar
          print*, 'interpolate_quad_ogrid: xxp=', xxp
          print*, 'interpolate_quad_ogrid: x0, y0, z0=', &
              xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
          print*, 'interpolate_quad_ogrid: x1, y1=', &
              xglobal_ogrid(ix0+1), yglobal_ogrid(iy0+1)
          print*, 'interpolate_quad_ogrid: x2, y2=', &
              xglobal_ogrid(ix0+2), yglobal_ogrid(iy0+2)
          print*, 'interpolate_quad_ogrid: ivar, gp=', ivar, gp
          print*, 'g(1,:,1)',gN(1,:,1)
          print*, 'g(2,:,1)',gN(2,:,1)
          print*, 'g(3,:,1)',gN(3,:,1)
          print*, '------------------'
          call fatal_error('interpolate_quad_ogrid','particle velocity interpolation error')
        endif
      endif
!
    endsubroutine interpolate_quad_ogrid
!***********************************************************************
!***********************************************************************
    subroutine interpolate_quadratic_spline(farr,ivar1,ivar2,xxp,gp,inear)
!
!  Quadratic spline interpolation of the function g to the point xxp=(xp,yp,zp).
!
!  10-jun-06/anders: coded
!
      real, dimension (:,:,:,:) :: farr
      integer, intent(in) :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (3), intent(in) :: inear
!
      
!TODO
      real, dimension(inear(1)-1:inear(1)+1,inear(2)-1:inear(2)+1,inear(3)-1:inear(3)+1,ivar2-ivar1+1) :: f
      real :: fac_x_m1, fac_x_00, fac_x_p1
      real :: fac_y_m1, fac_y_00, fac_y_p1
      real :: fac_z_m1, fac_z_00, fac_z_p1
      real :: dxp0, dyp0, dzp0
      integer :: ix0, iy0, iz0
!
      intent(in)  :: farr, xxp
      intent(out) :: gp
!TODO
      integer :: i
      f(:,:,:,:)=farr(:,:,:,:)
!
!  Redefine the interpolation point in coordinates relative to nearest grid
!  point and normalize with the cell size.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      dxp0=(xxp(1)-xglobal(ix0))*dx1grid(ix0-nghost)
      dyp0=(xxp(2)-yglobal(iy0))*dy1grid(iy0-nghost)
      dzp0=(xxp(3)-zglobal(iz0))*dz1grid(iz0-nghost)
!
!  Interpolation formulae.
!
      if (dimensionality==0) then
        gp=f(ix0,iy0,iz0,ivar1:ivar2)
      elseif (dimensionality==1) then
        if (nxgrid/=1) then
          gp = 0.5*(0.5-dxp0)**2*f(ix0-1,iy0,iz0,ivar1:ivar2) + &
                  (0.75-dxp0**2)*f(ix0  ,iy0,iz0,ivar1:ivar2) + &
               0.5*(0.5+dxp0)**2*f(ix0+1,iy0,iz0,ivar1:ivar2)
        endif
        if (nygrid/=1) then
          gp = 0.5*(0.5-dyp0)**2*f(ix0,iy0-1,iz0,ivar1:ivar2) + &
                  (0.75-dyp0**2)*f(ix0,iy0  ,iz0,ivar1:ivar2) + &
               0.5*(0.5+dyp0)**2*f(ix0,iy0+1,iz0,ivar1:ivar2)
        endif
        if (nzgrid/=1) then
          gp = 0.5*(0.5-dzp0)**2*f(ix0,iy0,iz0-1,ivar1:ivar2) + &
                  (0.75-dzp0**2)*f(ix0,iy0,iz0  ,ivar1:ivar2) + &
               0.5*(0.5+dzp0)**2*f(ix0,iy0,iz0+1,ivar1:ivar2)
        endif
      elseif (dimensionality==2) then
        if (nxgrid==1) then
          fac_y_m1 = 0.5*(0.5-dyp0)**2
          fac_y_00 = 0.75-dyp0**2
          fac_y_p1 = 0.5*(0.5+dyp0)**2
          fac_z_m1 = 0.5*(0.5-dzp0)**2
          fac_z_00 = 0.75-dzp0**2
          fac_z_p1 = 0.5*(0.5+dzp0)**2
!
          gp= fac_y_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_y_00*( f(ix0,iy0  ,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0  ,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_z_00*( f(ix0,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
              fac_y_p1*( f(ix0,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_y_m1*( f(ix0,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 )
        elseif (nygrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_z_m1 = 0.5*(0.5-dzp0)**2
          fac_z_00 = 0.75-dzp0**2
          fac_z_p1 = 0.5*(0.5+dzp0)**2
!
          gp= fac_x_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_x_00*( f(ix0  ,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0  ,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_z_00*( f(ix0+1,iy0,iz0  ,ivar1:ivar2)*fac_x_p1 + &
                         f(ix0-1,iy0,iz0  ,ivar1:ivar2)*fac_x_m1 ) + &
              fac_x_p1*( f(ix0+1,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0+1,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_x_m1*( f(ix0-1,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0-1,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 )
        elseif (nzgrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_y_m1 = 0.5*(0.5-dyp0)**2
          fac_y_00 = 0.75-dyp0**2
          fac_y_p1 = 0.5*(0.5+dyp0)**2
!
          gp= fac_x_00*fac_y_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_x_00*( f(ix0  ,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0  ,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 ) + &
              fac_y_00*( f(ix0+1,iy0  ,iz0,ivar1:ivar2)*fac_x_p1 + &
                         f(ix0-1,iy0  ,iz0,ivar1:ivar2)*fac_x_m1 ) + &
              fac_x_p1*( f(ix0+1,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0+1,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 ) + &
              fac_x_m1*( f(ix0-1,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0-1,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 )
        endif
      elseif (dimensionality==3) then
        fac_x_m1 = 0.5*(0.5-dxp0)**2
        fac_x_00 = 0.75-dxp0**2
        fac_x_p1 = 0.5*(0.5+dxp0)**2
        fac_y_m1 = 0.5*(0.5-dyp0)**2
        fac_y_00 = 0.75-dyp0**2
        fac_y_p1 = 0.5*(0.5+dyp0)**2
        fac_z_m1 = 0.5*(0.5-dzp0)**2
        fac_z_00 = 0.75-dzp0**2
        fac_z_p1 = 0.5*(0.5+dzp0)**2
!
        gp= fac_x_00*fac_y_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
            fac_x_00*fac_y_00*( f(ix0  ,iy0  ,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0  ,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_z_00*( f(ix0  ,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0  ,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
            fac_y_00*fac_z_00*( f(ix0+1,iy0  ,iz0  ,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0  ,ivar1:ivar2)*fac_x_m1 ) + &
            fac_x_p1*fac_y_p1*( f(ix0+1,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0+1,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_p1*fac_y_m1*( f(ix0+1,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0+1,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_m1*fac_y_p1*( f(ix0-1,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0-1,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_m1*fac_y_m1*( f(ix0-1,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0-1,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_y_p1*( f(ix0  ,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_y_m1*( f(ix0  ,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_y_00*fac_z_p1*( f(ix0+1,iy0  ,iz0+1,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0+1,ivar1:ivar2)*fac_x_m1 ) + &
            fac_y_00*fac_z_m1*( f(ix0+1,iy0  ,iz0-1,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0-1,ivar1:ivar2)*fac_x_m1 ) + &
            fac_z_00*fac_x_p1*( f(ix0+1,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0+1,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
            fac_z_00*fac_x_m1*( f(ix0-1,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0-1,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 )
      endif
!
      if (lcheck_interpolation) then
        do i=1,ivar2-ivar1+1
          if ((gp(i)>maxval(f(:,:,4,i))) .or. &
              (gp(i)<minval(f(:,:,4,i))) .or. &
              (gp(i)/=gp(i))) then
            if (gp(i)>maxval(f(:,:,4,i))) then
              print*, 'interpolate_quadratic_spline: interpolated value is LARGER than'
              print*, 'interpolate_quadratic_spline: values at the corner points!'
            elseif  (gp(i)<minval(f(:,:,4,i))) then
              print*, 'interpolate_quadratic_spline: interpolated value is SMALLER than'
              print*, 'interpolate_quadratic_spline: values at the corner points!'
            elseif (gp(i)/=gp(i)) then
              print*, 'interpolate_quadratic_spline: interpolated value is NaN'
            endif
            print*, 'iproc = ', iproc
            print*, 'dimensionality = ',dimensionality
            print*, 'interpolate_quadratic_spline: xxp=', xxp
            print*, 'interpolate_quadratic_spline: i, gp(i)=', i, gp(i)
            print*, 'Nearest neighbours: xglobal - ', xglobal(ix0-1:ix0+1)
            print*, 'Nearest neighbours: yglobal - ', yglobal(iy0-1:iy0+1)
            !print*, 'interpolate_quadratic_spline: x0, y0, z0=', &
                !xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
            print*, 'interpolate_quadratic_spline: f(ix0-1,iy0-1:iy0+1,4,i)=', f(ix0-1,iy0-1:iy0+1,4,i)
            print*, 'interpolate_quadratic_spline: f(ix0  ,iy0-1:iy0+1,4,i)=', f(ix0  ,iy0-1:iy0+1,4,i)
            print*, 'interpolate_quadratic_spline: f(ix0+1,iy0-1:iy0+1,4,i)=', f(ix0+1,iy0-1:iy0+1,4,i)
            print*, '------------------'
! Commented out the line below for the moment since it is too sensitive when lchemistry
!            call fatal_error('interpolate_quadratic_spline','interpolation error, quadratic spline')
          endif
        enddo
      endif
    endsubroutine interpolate_quadratic_spline
!***********************************************************************
    subroutine interpolate_quadratic_sp_og(farr,ivar1,ivar2,xxp,gp,inear)
!
!  Quadratic spline interpolation of the function g to the point xxp=(xp,yp,zp).
!
!  10-jun-06/anders: coded
!
      real, dimension (:,:,:,:) :: farr
      integer, intent(in) :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (3), intent(in) :: inear
!TODO
      real, dimension(inear(1)-1:inear(1)+1,inear(2)-1:inear(2)+1,inear(3)-1:inear(3)+1,ivar2-ivar1+1) :: f
      real :: fac_x_m1, fac_x_00, fac_x_p1
      real :: fac_y_m1, fac_y_00, fac_y_p1
      real :: fac_z_m1, fac_z_00, fac_z_p1
      real :: dxp0, dyp0, dzp0
      integer :: ix0, iy0, iz0
!
      intent(in)  :: farr, xxp
      intent(out) :: gp
      f(:,:,:,:)=farr(:,:,:,:)
!
!  Redefine the interpolation point in coordinates relative to nearest grid
!  point and normalize with the cell size.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      dxp0=(xxp(1)-xglobal_ogrid(ix0))*dx1global_ogrid(ix0)
      dyp0=(xxp(2)-yglobal_ogrid(iy0))*dy1global_ogrid(iy0)
      dzp0=(xxp(3)-zglobal_ogrid(iz0))*dz1global_ogrid(iz0)
!
!  Interpolation formulae.
!
      if (dimensionality==0) then
        gp=f(ix0,iy0,iz0,ivar1:ivar2)
      elseif (dimensionality==1) then
        if (nxgrid_ogrid/=1) then
          gp = 0.5*(0.5-dxp0)**2*f(ix0-1,iy0,iz0,ivar1:ivar2) + &
                  (0.75-dxp0**2)*f(ix0  ,iy0,iz0,ivar1:ivar2) + &
               0.5*(0.5+dxp0)**2*f(ix0+1,iy0,iz0,ivar1:ivar2)
        endif
        if (nygrid_ogrid/=1) then
          gp = 0.5*(0.5-dyp0)**2*f(ix0,iy0-1,iz0,ivar1:ivar2) + &
                  (0.75-dyp0**2)*f(ix0,iy0  ,iz0,ivar1:ivar2) + &
               0.5*(0.5+dyp0)**2*f(ix0,iy0+1,iz0,ivar1:ivar2)
        endif
        if (nzgrid_ogrid/=1) then
          gp = 0.5*(0.5-dzp0)**2*f(ix0,iy0,iz0-1,ivar1:ivar2) + &
                  (0.75-dzp0**2)*f(ix0,iy0,iz0  ,ivar1:ivar2) + &
               0.5*(0.5+dzp0)**2*f(ix0,iy0,iz0+1,ivar1:ivar2)
        endif
      elseif (dimensionality==2) then
        if (nxgrid_ogrid==1) then
          fac_y_m1 = 0.5*(0.5-dyp0)**2
          fac_y_00 = 0.75-dyp0**2
          fac_y_p1 = 0.5*(0.5+dyp0)**2
          fac_z_m1 = 0.5*(0.5-dzp0)**2
          fac_z_00 = 0.75-dzp0**2
          fac_z_p1 = 0.5*(0.5+dzp0)**2
!
          gp= fac_y_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_y_00*( f(ix0,iy0  ,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0  ,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_z_00*( f(ix0,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
              fac_y_p1*( f(ix0,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_y_m1*( f(ix0,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 )
        elseif (nygrid_ogrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_z_m1 = 0.5*(0.5-dzp0)**2
          fac_z_00 = 0.75-dzp0**2
          fac_z_p1 = 0.5*(0.5+dzp0)**2
!
          gp= fac_x_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_x_00*( f(ix0  ,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0  ,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_z_00*( f(ix0+1,iy0,iz0  ,ivar1:ivar2)*fac_x_p1 + &
                         f(ix0-1,iy0,iz0  ,ivar1:ivar2)*fac_x_m1 ) + &
              fac_x_p1*( f(ix0+1,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0+1,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_x_m1*( f(ix0-1,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0-1,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 )
        elseif (nzgrid_ogrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_y_m1 = 0.5*(0.5-dyp0)**2
          fac_y_00 = 0.75-dyp0**2
          fac_y_p1 = 0.5*(0.5+dyp0)**2
!
          gp= fac_x_00*fac_y_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_x_00*( f(ix0  ,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0  ,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 ) + &
              fac_y_00*( f(ix0+1,iy0  ,iz0,ivar1:ivar2)*fac_x_p1 + &
                         f(ix0-1,iy0  ,iz0,ivar1:ivar2)*fac_x_m1 ) + &
              fac_x_p1*( f(ix0+1,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0+1,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 ) + &
              fac_x_m1*( f(ix0-1,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0-1,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 )
        endif
      elseif (dimensionality==3) then
        fac_x_m1 = 0.5*(0.5-dxp0)**2
        fac_x_00 = 0.75-dxp0**2
        fac_x_p1 = 0.5*(0.5+dxp0)**2
        fac_y_m1 = 0.5*(0.5-dyp0)**2
        fac_y_00 = 0.75-dyp0**2
        fac_y_p1 = 0.5*(0.5+dyp0)**2
        fac_z_m1 = 0.5*(0.5-dzp0)**2
        fac_z_00 = 0.75-dzp0**2
        fac_z_p1 = 0.5*(0.5+dzp0)**2
!
        gp= fac_x_00*fac_y_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
            fac_x_00*fac_y_00*( f(ix0  ,iy0  ,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0  ,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_z_00*( f(ix0  ,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0  ,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
            fac_y_00*fac_z_00*( f(ix0+1,iy0  ,iz0  ,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0  ,ivar1:ivar2)*fac_x_m1 ) + &
            fac_x_p1*fac_y_p1*( f(ix0+1,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0+1,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_p1*fac_y_m1*( f(ix0+1,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0+1,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_m1*fac_y_p1*( f(ix0-1,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0-1,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_m1*fac_y_m1*( f(ix0-1,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0-1,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_y_p1*( f(ix0  ,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_y_m1*( f(ix0  ,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_y_00*fac_z_p1*( f(ix0+1,iy0  ,iz0+1,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0+1,ivar1:ivar2)*fac_x_m1 ) + &
            fac_y_00*fac_z_m1*( f(ix0+1,iy0  ,iz0-1,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0-1,ivar1:ivar2)*fac_x_m1 ) + &
            fac_z_00*fac_x_p1*( f(ix0+1,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0+1,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
            fac_z_00*fac_x_m1*( f(ix0-1,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0-1,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 )
      endif
!
    endsubroutine interpolate_quadratic_sp_og
!***********************************************************************
  subroutine poly_interp_cart(ivar1,ivar2,xyz_ip,gp,id,f_cartesian,order)
    use General, only: polynomial_interpolation
    real, dimension(mx,my,mz,mfarray), intent(in) :: f_cartesian

    real, dimension(3), intent(in) :: xyz_ip
    integer, intent(in) :: id
    integer, intent(in) :: ivar1,ivar2
    integer, intent(in) :: order
    real, dimension(ivar2-ivar1+1), intent(out) :: gp

    real, dimension(order) :: xa1,xa2,ya
    real, dimension(order,order-1) :: up
    real, dimension(ivar2-ivar1+1,order-1) :: gpp
    integer :: ix0,iy0,ivar,j,jj,ii1,ii2

    ix0=cartesian_to_curvilinear(id)%ind_local_neighbour(1)
    iy0=cartesian_to_curvilinear(id)%ind_local_neighbour(2)
    do j=1,order
      jj=-floor(order*0.5)+j-1
      xa1(j)=xglobal(ix0+jj)
      xa2(j)=yglobal(iy0+jj)
    enddo
    ii1=-floor(order*0.5)+1-1
    ii2=-floor(order*0.5)+order-1
    !print*, 'xyz[:,',id,0,']=[',xa1(1),',',xa1(2),',',xa1(3),',',xa1(4),',',xa1(5),']'
    !print*, 'xyz[:,',id,1,']=[',xa2(1),',',xa2(2),',',xa2(3),',',xa2(4),',',xa2(5),']'
    !!print*, 'cartesian:xa2',xa2
    !print*, 'rthp[',id,',:]=[',xyz_ip(1),',',xyz_ip(2),']'
    do ivar=ivar1,ivar2
      do j=1,order
        jj=-floor(order*0.5)+j-1
        ya=f_cartesian(ix0+ii1:ix0+ii2,iy0+jj,4,ivar)
        call polynomial_interpolation(xa1, ya, xyz_ip(1), up(j,:))
      enddo
      ya=up(:,order-1)
      call polynomial_interpolation(xa2, ya, xyz_ip(2), gpp(ivar,:))
    enddo

    gp=gpp(:,order-1)

  endsubroutine poly_interp_cart
!***********************************************************************
  subroutine poly_interp_curv(ivar1,ivar2,xyz_ip,gp,id,order)
    use General, only: polynomial_interpolation

    real, dimension(3), intent(in) :: xyz_ip
    integer, intent(in) :: id
    integer, intent(in) :: ivar1,ivar2
    integer, intent(in) :: order
    real, dimension(ivar2-ivar1+1), intent(out) :: gp

    real, dimension(order) :: xa1,xa2,ya
    real, dimension(order,order-1) :: up
    real, dimension(ivar2-ivar1+1,order-1) :: gpp
    integer :: ix0,iy0,ivar,j,jj,ii1,ii2

    ix0=curvilinear_to_cartesian(id)%ind_local_neighbour(1)
    iy0=curvilinear_to_cartesian(id)%ind_local_neighbour(2)
    do j=1,order
      jj=-floor(order*0.5)+j-1
      xa1(j)=xglobal_ogrid(ix0+jj)
      xa2(j)=yglobal_ogrid(iy0+jj)
    enddo
    ii1=-floor(order*0.5)+1-1
    ii2=-floor(order*0.5)+order-1
    do ivar=ivar1,ivar2
      do j=1,order
        jj=-floor(order*0.5)+j-1
        ya=f_ogrid(ix0+ii1:ix0+ii2,iy0+jj,4,ivar)
        call polynomial_interpolation(xa1, ya, xyz_ip(1), up(j,:))
      enddo
      ya=up(:,order-1)
      call polynomial_interpolation(xa2, ya, xyz_ip(2), gpp(ivar,:))
    enddo

    gp=gpp(:,order-1)

  endsubroutine poly_interp_curv
!***********************************************************************
    real function radius_ogrid(xp,yp) 
!
!  Transform cartesian coordinates xxp to polar coordinates on the cylindrical ogrid
!  Only return r-direction. 
!  Input coordinates xxp can be two- or three-dimensional
!
!  02-aug-17/Jorgen: Coded
!
    real, intent(in) :: xp, yp
!
      radius_ogrid = sqrt((xp-xorigo_ogrid(1))**2 + (yp-xorigo_ogrid(2))**2)
!
    endfunction radius_ogrid
!***********************************************************************
    subroutine get_polar_coords_2D(xxp,rad,theta)
!
!  Transform cartesian coordinates xxp to polar coordinates on the cylindrical ogrid
!  Only return two-dimensional array
!  Input coordinates xxp can be two- or three-dimensional
!
!  02-aug-17/Jorgen: Coded
!
    real, dimension(:), intent(in) :: xxp
    real, intent(out) :: rad, theta
!
    real :: xr, yr
!
      xr=xxp(1)-xorigo_ogrid(1)
      yr=xxp(2)-xorigo_ogrid(2)
      rad=sqrt(xr**2+yr**2)
      theta=atan2(yr,xr)
!
    endsubroutine get_polar_coords_2D
!***********************************************************************
    subroutine get_polar_coords_3D(xxp,rthz)
!
!  Transform cartesian coordinates xxp to polar coordinates on the cylindrical ogrid
!  Input coordinates must be three-dimensional
!
!  02-aug-17/Jorgen: Coded
!
    real, dimension(3), intent(in) :: xxp
    real, dimension(3), intent(out) :: rthz
!
      call get_polar_coords_2D(xxp,rthz(1),rthz(2))
      rthz(3) = xxp(3)
!
    endsubroutine get_polar_coords_3D
!***********************************************************************
    subroutine get_polar_coords_3D_alt(xp,yp,zp,rthz)
!
!  Transform cartesian coordinates xxp to polar coordinates on the cylindrical ogrid
!  Input coordinates must be three-dimensional
!
!  02-aug-17/Jorgen: Coded
!
    real, intent(in) :: xp, yp, zp
    real, dimension(3), intent(out) :: rthz
    real :: xr,yr
!
      xr=xp-xorigo_ogrid(1)
      yr=yp-xorigo_ogrid(2)
      rthz=(/ sqrt(xr**2+yr**2),atan2(yr,xr),zp /)
!
    endsubroutine get_polar_coords_3D_alt
!***********************************************************************
    subroutine adjust_inear_cart(inear,xxp)
!
!  Adjust inear coordinates to guarantee that they point to the CLOSEST point to xxp,
!  not to the bottom left corner of the cell that contains xxp.
!  Necessary for interpolation with asymmetric stencils (e.g., quadratic spline)
!
!  07-sep-17/Jorgen: Coded
!
      integer, dimension(3),intent(inout) :: inear
      real, dimension(3), intent(in) :: xxp
!
      if((xxp(1)-x(inear(1)))>(x(inear(1)+1)-xxp(1))) inear(1) = inear(1)+1
      if((xxp(2)-y(inear(2)))>(y(inear(2)+1)-xxp(2))) inear(2) = inear(2)+1
      if(nzgrid_ogrid>1) then
        if((xxp(3)-z(inear(3)))>(z(inear(3)+1)-xxp(3))) inear(3) = inear(3)+1
      endif
!
    endsubroutine adjust_inear_cart
!***********************************************************************
    subroutine adjust_inear_curv(inear,xxp)
!
!  Adjust inear coordinates to guarantee that they point to the CLOSEST point to xxp,
!  not to the bottom left corner of the cell that contains xxp.
!  Necessary for interpolation with asymmetric stencils (e.g., quadratic spline)
!
!  07-sep-17/Jorgen: Coded
!
      integer, dimension(3),intent(inout) :: inear
      real, dimension(3), intent(in) :: xxp
!
      if((xxp(1)-x_ogrid(inear(1)))>(x_ogrid(inear(1)+1)-xxp(1))) inear(1) = inear(1)+1
      if((xxp(2)-y_ogrid(inear(2)))>(y_ogrid(inear(2)+1)-xxp(2))) inear(2) = inear(2)+1
      if(nzgrid_ogrid>1) then
        if((xxp(3)-z_ogrid(inear(3)))>(z_ogrid(inear(3)+1)-xxp(3))) inear(3) = inear(3)+1
      endif
!
    endsubroutine adjust_inear_curv
!***********************************************************************
    subroutine adjust_inear_cart_glob(inear_glob,xxp)
!
!  Adjust inear GLOBAL coordinates to guarantee that they point to the CLOSEST point to xxp,
!  not to the bottom left corner of the cell that contains xxp.
!  Necessary for interpolation with asymmetric stencils (e.g., quadratic spline)
!
!  14-sep-17/Jorgen: Coded
!
      integer, dimension(3),intent(inout) :: inear_glob
      real, dimension(3), intent(in) :: xxp
!
      if((xxp(1)-xglobal(inear_glob(1)))>(xglobal(inear_glob(1)+1)-xxp(1))) inear_glob(1) = inear_glob(1)+1
      if((xxp(2)-yglobal(inear_glob(2)))>(yglobal(inear_glob(2)+1)-xxp(2))) inear_glob(2) = inear_glob(2)+1
      if(nzgrid>1) then
        if((xxp(3)-zglobal(inear_glob(3)))>(zglobal(inear_glob(3)+1)-xxp(3))) inear_glob(3) = inear_glob(3)+1
      endif
!
    endsubroutine adjust_inear_cart_glob
!***********************************************************************
    subroutine adjust_inear_curv_glob(inear_glob,xxp)
!
!  Adjust inear GLOBAL coordinates to guarantee that they point to the CLOSEST point to xxp,
!  not to the bottom left corner of the cell that contains xxp.
!  Necessary for interpolation with asymmetric stencils (e.g., quadratic spline)
!
!  14-sep-17/Jorgen: Coded
!
      integer, dimension(3),intent(inout) :: inear_glob
      real, dimension(3), intent(in) :: xxp
!
      if((xxp(1)-xglobal_ogrid(inear_glob(1)))>(xglobal_ogrid(inear_glob(1)+1)-xxp(1))) inear_glob(1) = inear_glob(1)+1
      if((xxp(2)-yglobal_ogrid(inear_glob(2)))>(yglobal_ogrid(inear_glob(2)+1)-xxp(2))) inear_glob(2) = inear_glob(2)+1
      if(nzgrid_ogrid>1) then
        if((xxp(3)-zglobal_ogrid(inear_glob(3)))>(zglobal_ogrid(inear_glob(3)+1)-xxp(3))) inear_glob(3) = inear_glob(3)+1
      endif
!
    endsubroutine adjust_inear_curv_glob
!***********************************************************************
    subroutine set_interpolation_limits
!
!  Set interpolation zone for curvilinear to cartesian grid
!  Make sure that no points outside x_ogrid(l2_ogrid) are used
!
      use mpicomm, only: mpibcast_real
      real :: dx_outer, tmp_rad, min_rad, min_tmp_rad
      integer :: ii,jj,i3,j3

        if(lroot) then
          dx_outer = 1./dx1grid_ogrid(nxgrid_ogrid)
          if(interpolation_method==1) then
            r_int_outer=r_ogrid-dx_outer*0.01-dx_outer*interp_shift
          elseif (interpolation_method==2 .or. interpolation_method==3 &
                 .or. interpolation_method==5) then
            r_int_outer=r_ogrid-dx_outer*0.51-dx_outer*interp_shift
            if((xgrid_ogrid(nxgrid_ogrid)-r_int_outer)<(r_int_outer-xgrid_ogrid(nxgrid_ogrid-1))) then
              print*, 'WARNING: An error occured when setting interpolation zone.'
              print*, '         Zone adjusted.'
              print*, 'iproc, r_int_outer first, r_int_outer second',&
                iproc,r_int_outer,r_ogrid-dx_outer*1.01
              r_int_outer=r_ogrid-dx_outer*1.01-dx_outer*interp_shift
            endif
            if(interpolation_method==5) then
              print*, 'WARNING: Polynomal interpolation used, you better know what you are doing!'
            endif
 !         elseif (interpolation_method==4) then
 !           r_int_outer=r_ogrid-dx_outer*1.51-dx_outer*interp_shift
 !           if((xgrid_ogrid(nxgrid_ogrid-1)-r_int_outer)<(r_int_outer-xgrid_ogrid(nxgrid_ogrid-2))) then
 !             print*, 'WARNING: An error occured when setting interpolation zone.'
 !             print*, '         Zone adjusted.'
 !             r_int_outer=r_ogrid-dx_outer*2.01-dx_outer*interp_shift
 !           endif
          elseif (mod(interpolation_method,2)==0) then
            r_int_outer=r_ogrid-dx_outer*((interpolation_method/2-0.5)+0.01)-dx_outer*interp_shift
            if((xgrid_ogrid(nxgrid_ogrid-1)-r_int_outer)<(r_int_outer-xgrid_ogrid(nxgrid_ogrid-2))) then
              print*, 'WARNING: An error occured when setting interpolation zone.'
              print*, '         Zone adjusted.'
              r_int_outer=r_ogrid-dx_outer*((interpolation_method/2)+0.01)-dx_outer*interp_shift
            endif
          else
            call fatal_error('initialize_solid_cells','interpolation method does not exist')
          endif
        endif
!
!  Broadcast the value set for r_int_outer
!
        call mpibcast_real(r_int_outer)
!
!  Set limit of the interpolation zone, r_int_inner
!
        if(interpolation_method<5 .or. interpolation_method>5) then
          min_rad=r_int_outer
          do ii = l1,l2
            do jj = m1,m2
              tmp_rad = radius_ogrid(x(ii),y(jj))

              if(tmp_rad>r_int_outer.and.tmp_rad<(r_int_outer+5*dxmax)) then
                do i3=-3,3
                  do j3=-3,3
                    min_tmp_rad = radius_ogrid(x(ii+i3),y(jj+j3))
                    if(min_tmp_rad<min_rad) min_rad=min_tmp_rad
                  enddo
                enddo
              endif
            enddo
          enddo
          r_int_inner = min(r_int_outer-dxmax*3,min_rad-0.01*dxmax)
        else
          r_int_inner_poly=x_ogrid(l1_ogrid+floor(interpol_order_poly*0.5))
          print*, 'Polynomial integration: r_int_outer, r_int_inner, r_int_inner_poly' &
                  , r_int_outer, r_int_inner, r_int_inner_poly
          print*, 'Polynomial integration: You should no what you are doing...'
        endif


    endsubroutine set_interpolation_limits
!***********************************************************************
    subroutine del2v_etc_ogrid(f,k,del2,graddiv)
!
!  Calculates a number of second derivative expressions of a vector
!  outputs a number of different vector fields.
!  gradcurl is not the vector gradient.
!  Surprisingly, calling derij only if graddiv or curlcurl are present
!  does not speed up the code on Mephisto @ 32x32x64.
!
!  12-sep-01/axel: coded
!  15-mar-07/wlad: added cylindrical coordinates
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray_ogrid) :: f
      real, dimension (nx_ogrid,3,3) :: fjji,fijj
      real, dimension (nx_ogrid,3), optional :: del2,graddiv
      real, dimension (nx_ogrid) :: tmp
      integer :: i,j,k,k1
!
      intent(in) :: f,k
      intent(out) :: del2,graddiv
!
!  calculate f_{i,jj} and f_{j,ji}
!  AJ: graddiv needs diagonal elements from the first tmp (derij only sets
!      off-diagonal elements)
!
      k1=k-1
      do i=1,3
      do j=1,3
        if (present(del2) .or. present(graddiv)) then
          call der2_ogrid(f,k1+i,tmp,  j) 
          fijj(:,i,j)=tmp  ! f_{i,jj}
        endif
        if (present(graddiv)) then
          call derij_ogrid(f,k1+j,tmp,j,i) 
          fjji(:,i,j)=tmp  ! f_{j,ji}
        endif
      enddo
      enddo
!
!  the diagonal terms have not been set in derij; do this now
!  ** They are automatically set above, because derij   **
!  ** doesn't overwrite the value of tmp for i=j!       **
!
!     do j=1,3
!       fjji(:,j,j)=fijj(:,j,j)
!     enddo
!
!  calculate f_{i,jk} for i /= j /= k
!
!
!  del2
!
      if (present(del2)) then
        do i=1,3
          del2(:,i)=fijj(:,i,1)+fijj(:,i,2)+fijj(:,i,3)
        enddo
        !r-component
        call der_ogrid(f,k1+2,tmp,2)
        del2(:,1)=del2(:,1) -(2*tmp+f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,k1+1))*rcyl_mn2_ogrid
        call der_ogrid(f,k1+1,tmp,1)
        del2(:,1)=del2(:,1) + tmp*rcyl_mn1_ogrid
        !phi-component
        call der_ogrid(f,k1+1,tmp,2)
        del2(:,2)=del2(:,2) +(2*tmp-f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,k1+2))*rcyl_mn2_ogrid
        call der_ogrid(f,k1+2,tmp,1)
        del2(:,2)=del2(:,2) + tmp*rcyl_mn1_ogrid
        !z-component
        call der_ogrid(f,k1+3,tmp,1)
        del2(:,3)=del2(:,3) + tmp*rcyl_mn1_ogrid
      endif
!
      if (present(graddiv)) then
        do i=1,3
          graddiv(:,i)=fjji(:,i,1)+fjji(:,i,2)+fjji(:,i,3)
        enddo
        call der_ogrid(f,k1+1,tmp,1)
        graddiv(:,1)=graddiv(:,1)+tmp*rcyl_mn1_ogrid - f(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,k1+1)*rcyl_mn2_ogrid 
        call der_ogrid(f,2,tmp,2)
        graddiv(:,1)=graddiv(:,1)-rcyl_mn1_ogrid*tmp
        call der_ogrid(f,k1+1,tmp,2)
        graddiv(:,2)=graddiv(:,2)+tmp*rcyl_mn1_ogrid
        call der_ogrid(f,k1+1,tmp,3)
        graddiv(:,3)=graddiv(:,3)+tmp*rcyl_mn1_ogrid
      endif
    endsubroutine del2v_etc_ogrid
!***********************************************************************
  subroutine run_tests_ogrid

    real :: velocity
    real :: R2
    real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid) ::  f_pflow
    real, dimension (mx_ogrid, my_ogrid,2,2) ::  df_pflow, df_pflow_ex
    real, dimension (mx_ogrid, my_ogrid,2,2) ::  df2_pflow, df2_pflow_ex
    real, dimension (mx_ogrid, my_ogrid,2,2,2) ::  df2ij_pflow, df2ij_pflow_ex
    real, dimension (2,2) :: df_twonorm, df2_twonorm
    real, dimension (3) :: graddivu_twonorm
    real, dimension (3) :: graddivu2_twonorm
    real, dimension (2,2,2) :: df2ij_twonorm
    real, dimension (nx_ogrid,3)     :: u_vec
    real, dimension (nx_ogrid,3,3) :: uij_tensor
    real, dimension (mx_ogrid,my_ogrid,3) :: graddivu_vec_exact
    real, dimension (mx_ogrid,my_ogrid,3) :: graddivu_vec
    real, dimension (mx_ogrid,my_ogrid,3) :: graddivu_vec_exact2
    integer :: i,j,k
  
    R2 = cylinder_radius**2
    velocity = 1.0
!
!  Set up potential flow
!  No velocity in z-direction, or density variation
!
    do i=1,mx_ogrid
      do j=1,my_ogrid
        do k=1,mz_ogrid
          f_pflow(i,j,k,1) = velocity*(1-R2/(x_ogrid(i)**2))*cos(y_ogrid(j))
          f_pflow(i,j,k,2) =-velocity*(1+R2/(x_ogrid(i)**2))*sin(y_ogrid(j))
          f_pflow(i,j,k,3) =0.
          f_pflow(i,j,k,4) =1.
        enddo
      enddo
    enddo
!
!  Compute first order derivatives
!
    n_ogrid=4
    do m_ogrid=m1_ogrid,m2_ogrid
      call der_ogrid(f_pflow,1,df_pflow(l1_ogrid:l2_ogrid,m_ogrid,1,1),1)
      call der_ogrid(f_pflow,1,df_pflow(l1_ogrid:l2_ogrid,m_ogrid,1,2),2)
      call der_ogrid(f_pflow,2,df_pflow(l1_ogrid:l2_ogrid,m_ogrid,2,1),1)
      call der_ogrid(f_pflow,2,df_pflow(l1_ogrid:l2_ogrid,m_ogrid,2,2),2)
    enddo
!
!  Compute second order derivatives
!
    do m_ogrid=m1_ogrid,m2_ogrid
      call der2_ogrid(f_pflow,1,df2_pflow(l1_ogrid:l2_ogrid,m_ogrid,1,1),1)
      call der2_ogrid(f_pflow,1,df2_pflow(l1_ogrid:l2_ogrid,m_ogrid,1,2),2)
      call der2_ogrid(f_pflow,2,df2_pflow(l1_ogrid:l2_ogrid,m_ogrid,2,1),1)
      call der2_ogrid(f_pflow,2,df2_pflow(l1_ogrid:l2_ogrid,m_ogrid,2,2),2)
    enddo

!
!  Compute mixed derivatives
!
    do m_ogrid=m1_ogrid,m2_ogrid
      call derij_ogrid(f_pflow,1,df2ij_pflow(l1_ogrid:l2_ogrid,m_ogrid,1,1,2),1,2)
      call derij_ogrid(f_pflow,1,df2ij_pflow(l1_ogrid:l2_ogrid,m_ogrid,1,2,1),2,1)
      call derij_ogrid(f_pflow,2,df2ij_pflow(l1_ogrid:l2_ogrid,m_ogrid,2,1,2),1,2)
      call derij_ogrid(f_pflow,2,df2ij_pflow(l1_ogrid:l2_ogrid,m_ogrid,2,2,1),2,1)
    enddo
!
!  Compute grad(div u)
!
    do m_ogrid=m1_ogrid,m2_ogrid
      u_vec(1:nx_ogrid,1)=f_pflow(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,1)
      u_vec(1:nx_ogrid,2)=f_pflow(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,2)
      u_vec(1:nx_ogrid,3)=f_pflow(l1_ogrid:l2_ogrid,m_ogrid,n_ogrid,3)
      call gij_ogrid(f_pflow,iuu,uij_tensor)

!      print*, 'dvr_dr  :', uij_tensor(:,1,1)-df_pflow(l1_ogrid:l2_ogrid,m_ogrid,1,1)
!      print*, 'dvr_dth :', uij_tensor(:,1,2)-df_pflow(l1_ogrid:l2_ogrid,m_ogrid,1,2)
!      print*, 'dvth_dr :', uij_tensor(:,2,1)-df_pflow(l1_ogrid:l2_ogrid,m_ogrid,2,1)
!      print*, 'dvth_dth:', uij_tensor(:,2,2)-df_pflow(l1_ogrid:l2_ogrid,m_ogrid,2,2)

      call gij_etc_ogrid(f_pflow,iuu,u_vec,uij_tensor,GRADDIV=graddivu_vec(l1_ogrid:l2_ogrid,m_ogrid,:))
      print*, 'gij_etc_ogrid:',graddivu_vec(l1_ogrid:l2_ogrid,m_ogrid,1)
      call del2v_etc_ogrid(f_pflow,iuu,GRADDIV=graddivu_vec(l1_ogrid:l2_ogrid,m_ogrid,:))
      print*, 'del2v_etc_ogr:',graddivu_vec(l1_ogrid:l2_ogrid,m_ogrid,1)
    enddo
!
!  Set up exact solutions to derivatives 
!  All values in cylindrical coordinates, hence d/dy = 1/r d/dth, etc.
!
    do i=1,mx_ogrid
      do j=1,my_ogrid
        df_pflow_ex(i,j,1,1) = 2.*velocity*R2*cos(y_ogrid(j))/(x_ogrid(i)**3)
        df_pflow_ex(i,j,1,2) = (velocity*sin(y_ogrid(j))*(R2/(x_ogrid(i)**2)-1.))/x_ogrid(i)
        df_pflow_ex(i,j,2,1) = 2.*velocity*R2*sin(y_ogrid(j))/(x_ogrid(i)**3)
        df_pflow_ex(i,j,2,2) = -(velocity*cos(y_ogrid(j))*(R2/(x_ogrid(i)**2)+1.))/x_ogrid(i)

        
        df2_pflow_ex(i,j,1,1) = -6.*velocity*R2*cos(y_ogrid(j))/(x_ogrid(i)**4)
        df2_pflow_ex(i,j,1,2) = (velocity*cos(y_ogrid(j))*(R2/(x_ogrid(i)**2)-1.))/(x_ogrid(i)**2)
        df2_pflow_ex(i,j,2,1) = -6.*velocity*R2*sin(y_ogrid(j))/(x_ogrid(i)**4)
        df2_pflow_ex(i,j,2,2) = (velocity*sin(y_ogrid(j))*(R2/(x_ogrid(i)**2)+1.))/(x_ogrid(i)**2)

        df2ij_pflow_ex(i,j,1,1,2) = -2.*velocity*R2*sin(y_ogrid(j))/(x_ogrid(i)**4)
        !df2ij_pflow_ex(i,j,1,2,1) = -2.*velocity*R2*sin(y_ogrid(j))/(x_ogrid(i)**4)
        df2ij_pflow_ex(i,j,1,2,1) = df2ij_pflow_ex(i,j,1,1,2)
        df2ij_pflow_ex(i,j,2,1,2) = 2.*velocity*R2*cos(y_ogrid(j))/(x_ogrid(i)**4)
        df2ij_pflow_ex(i,j,2,2,1) = df2ij_pflow_ex(i,j,2,1,2)
        !df2ij_pflow_ex(i,j,2,2,1) = 2.*velocity*R2*cos(y_ogrid(j))/(x_ogrid(i)**4)
        graddivu_vec_exact(i,j,1) = df2_pflow_ex(i,j,1,1)+df2ij_pflow_ex(i,j,2,1,2) + &
                                (1./x_ogrid(i))*(df_pflow_ex(i,j,1,1)-df_pflow_ex(i,j,2,2)) - &
                                (1./(x_ogrid(i)**2))*(f_pflow(i,j,1,1))
        graddivu_vec_exact2(i,j,2) = df2_pflow_ex(i,j,2,2)+df2ij_pflow_ex(i,j,1,1,2) + &
                                (1./x_ogrid(i))*(df_pflow_ex(i,j,1,2))
        graddivu_vec_exact2(i,j,1) = df2_pflow(i,j,1,1)+df2ij_pflow(i,j,2,1,2) + &
                                (1./x_ogrid(i))*(df_pflow(i,j,1,1)-df_pflow(i,j,2,2)) - &
                                (1./(x_ogrid(i)**2))*(f_pflow(i,j,1,1))
        graddivu_vec_exact2(i,j,2) = df2_pflow(i,j,2,2)+df2ij_pflow(i,j,1,1,2) + &
                                (1./x_ogrid(i))*(df_pflow(i,j,1,2)) 
          !f_pflow(i,j,k,1) = velocity*(1-R2/(x_ogrid(i)**2))*cos(y_ogrid(j))
      enddo
    enddo
!
!  Compute two-norms
!
    df_twonorm=0.
    df2_twonorm=0.
    df2ij_twonorm=0.
    graddivu_twonorm=0.
    graddivu2_twonorm=0.
    do i=l1_ogrid,l2_ogrid
      do j=m1_ogrid,m2_ogrid
        df_twonorm(1,1)=df_twonorm(1,1)+(df_pflow_ex(i,j,1,1)-df_pflow(i,j,1,1))**2 &
                        *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        df_twonorm(1,2)=df_twonorm(1,2)+(df_pflow_ex(i,j,1,2)-df_pflow(i,j,1,2))**2 &
                        *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        df_twonorm(2,1)=df_twonorm(2,1)+(df_pflow_ex(i,j,2,1)-df_pflow(i,j,2,1))**2 &
                        *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        df_twonorm(2,2)=df_twonorm(2,2)+(df_pflow_ex(i,j,2,2)-df_pflow(i,j,2,2))**2 &
                            *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        df2_twonorm(1,1)=df2_twonorm(1,1)+(df2_pflow_ex(i,j,1,1)-df2_pflow(i,j,1,1))**2 &
                        *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        df2_twonorm(1,2)=df2_twonorm(1,2)+(df2_pflow_ex(i,j,1,2)-df2_pflow(i,j,1,2))**2 &
                        *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        df2_twonorm(2,1)=df2_twonorm(2,1)+(df2_pflow_ex(i,j,2,1)-df2_pflow(i,j,2,1))**2 &
                        *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        df2_twonorm(2,2)=df2_twonorm(2,2)+(df2_pflow_ex(i,j,2,2)-df2_pflow(i,j,2,2))**2 &
                            *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        df2ij_twonorm(1,1,2)=df2ij_twonorm(1,1,2)+(df2ij_pflow_ex(i,j,1,1,2)-df2ij_pflow(i,j,1,1,2))**2 &
                        *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        df2ij_twonorm(1,2,1)=df2ij_twonorm(1,2,1)+(df2ij_pflow_ex(i,j,1,2,1)-df2ij_pflow(i,j,1,2,1))**2 &
                        *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        df2ij_twonorm(2,1,2)=df2ij_twonorm(2,1,2)+(df2ij_pflow_ex(i,j,2,1,2)-df2ij_pflow(i,j,2,1,2))**2 &
                        *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        df2ij_twonorm(2,2,1)=df2ij_twonorm(2,2,1)+(df2ij_pflow_ex(i,j,2,2,1)-df2ij_pflow(i,j,2,2,1))**2 &
                            *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        graddivu_twonorm(:)=graddivu_twonorm(:)+(graddivu_vec(i,j,:)-graddivu_vec_exact(i,j,:))**2 &
                            *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
        graddivu2_twonorm(:)=graddivu2_twonorm(:)+(graddivu_vec(i,j,:)-graddivu_vec_exact2(i,j,:))**2 &
                            *(1./dx_1_ogrid(i))*(x_ogrid(i)/dy_1_ogrid(j))
      enddo
    enddo
    print*, graddivu_twonorm
    print*, graddivu2_twonorm
    df_twonorm=sqrt(df_twonorm)
    df2_twonorm=sqrt(df2_twonorm)
    df2ij_twonorm=sqrt(df2ij_twonorm)
!
!  Print two-norms
!
    if(SBP) then
      open(10,file='runinfo.dat',status='unknown')
      open(1,file='SBP2norm_df.dat',status='unknown')
      open(2,file='SBP2norm_df2.dat',status='unknown')
      open(3,file='SBP2norm_df2ij.dat',status='unknown')
      write(10,*) '% Summation by parts'
    elseif(BDRY5) then
      open(10,file='runinfo.dat',status='unknown')
      open(1,file='BDRY5norm_df.dat',status='unknown')
      open(2,file='BDRY5norm_df2.dat',status='unknown')
      open(3,file='BDRY5norm_df2ij.dat',status='unknown')
      write(10,*) '# Fifth order boundary closures'
    else
      open(10,file='runinfo.dat',status='unknown')
      open(1,file='NOBDRYnorm_df.dat',status='unknown')
      open(2,file='NOBDRYnorm_df2.dat',status='unknown')
      open(3,file='NOBDRYnorm_df2ij.dat',status='unknown')
      write(10,*) '# No boundary condition given'
    endif
    write(10,*) '# nx_ogrid, ny_ogrid'
    write(10,*)  nx_ogrid,ny_ogrid
    write(*,*)  '# df_twonorm'
    write(*,*)  '# dvrdr, dvrdth'
    write(*,*)  '# dvthdr, dvthdth'
    write(*,*)  df_twonorm(1,1), df_twonorm(1,2)
    write(*,*)  df_twonorm(2,1), df_twonorm(2,2)
    write(*,*)  '# df2_twonorm'
    write(*,*)  '# d2vrdr2, d2vrdth2'
    write(*,*)  '# d2vthdr2, d2vthdth2'
    write(*,*)  df2_twonorm(1,1), df2_twonorm(1,2)
    write(*,*)  df2_twonorm(2,1), df2_twonorm(2,2)
    write(*,*)  '# df2ij_twonorm'
    write(*,*)  '# d2vrdrdth, d2vrdthdr'
    write(*,*)  '# d2vthdrdth, d2vthdthdr'
    write(*,*)  df2ij_twonorm(1,1,2), df2ij_twonorm(1,2,1)
    write(*,*)  df2ij_twonorm(2,1,2), df2ij_twonorm(2,2,1)
    close(1)
    close(2)
    close(3)
    close(10)

    open(111,file='r.dat',status='unknown')
    open(211,file='th.dat',status='unknown')
    open(100,file='x.dat',status='unknown')
    open(11,file='y.dat',status='unknown')
    open(20,file='vr.dat',status='unknown')
    open(21,file='vth.dat',status='unknown')
    open(29,file='dervr_r.dat',status='unknown')
    open(30,file='der2vr_r.dat',status='unknown')
    open(31,file='der2vr_th.dat',status='unknown')
    open(32,file='der2vth_r.dat',status='unknown')
    open(33,file='der2vth_th.dat',status='unknown')
    !open(32,file='der2r_indirect.dat',status='unknown')
    !open(33,file='der2th_indirect.dat',status='unknown')
    write(111,*) x_ogrid(l1_ogrid:l2_ogrid)
    write(211,*) y_ogrid(m1_ogrid:m2_ogrid)
    do i=l1_ogrid,l2_ogrid
      !do j=m1_ogrid,m2_ogrid
        write(100,*) x_ogrid(i)*cos(y_ogrid(m1_ogrid:m2_ogrid))
        write(11,*) x_ogrid(i)*sin(y_ogrid(m1_ogrid:m2_ogrid))
      !enddo
    enddo  
    do i=l1_ogrid,l2_ogrid
      write(20,*) f_pflow(i,m1_ogrid:m2_ogrid,1,1)
      write(21,*) f_pflow(i,m1_ogrid:m2_ogrid,1,2)
      write(29,*) df_pflow(i,m1_ogrid:m2_ogrid,1,1)
      write(30,*) df2_pflow(i,m1_ogrid:m2_ogrid,1,1)
      write(31,*) df2_pflow(i,m1_ogrid:m2_ogrid,1,2)
      write(32,*) df2_pflow(i,m1_ogrid:m2_ogrid,2,1)
      write(33,*) df2_pflow(i,m1_ogrid:m2_ogrid,2,2)
    enddo  
      
    close(111)
    close(211)
    close(100)
    close(11)
    close(20)
    close(21)
    close(29)
    close(30)
    close(31)
    close(32)
    close(33)

  endsubroutine run_tests_ogrid
!***********************************************************************
  subroutine initialize_pade_filter(f_og)
!
!  Initialization of high order padé filtering of solution array.
!  10th order filter requires extension of ghost zones in periodic
!  directions. Extended ghosts zones (halos) allocated here.
!
!  29-nov-17/Jorgen - Coded
!
    real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(in)::  f_og
    
    if(filter_Hsize==0) then 
      print*, 'WARNING: No need for filter halos, ghost zone large enough'
      print*, '         This will not work for parallel runs with current implemtation'
    elseif(filter_Hsize<0) then
      call fatal_error('initialize_pade_filter','Negative filter halo size!')
    elseif(filter_Hsize>nghost) then
      ! Requres a modification of mpi-buffers, not yet implemented
      call fatal_error('initialize_pade_filter','Filter halo too large!')
    endif
    !
    allocate(f_filterH_lowerx(filter_Hsize,my_ogrid,nz_ogrid,mfarray_ogrid))
    allocate(f_filterH_upperx(filter_Hsize,my_ogrid,nz_ogrid,mfarray_ogrid))
    allocate(f_filterH_lowery(mx_ogrid,filter_Hsize,nz_ogrid,mfarray_ogrid))
    allocate(f_filterH_uppery(mx_ogrid,filter_Hsize,nz_ogrid,mfarray_ogrid))
    !
  endsubroutine initialize_pade_filter
!***********************************************************************
  subroutine communicate_filter_zones(f_og,f_Hlox,f_Hupx,f_Hloy,f_Hupy)
    
    use Solid_Cells_Mpicomm, only: initiate_isendrcv_bdry_filter, finalize_isendrcv_bdry_filter
    
    real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid) ::  f_og
    real, dimension (filter_Hsize,my_ogrid,nz_ogrid,mfarray_ogrid) ::  f_Hlox,f_Hupx
    real, dimension (mx_ogrid,filter_Hsize,nz_ogrid,mfarray_ogrid) ::  f_Hloy,f_Hupy
    
    intent(in) :: f_og
    intent(inout) :: f_Hlox,f_Hupx,f_Hloy,f_Hupy
!
!  Communicate additional ghost zones needed for 10th order filter.
!
!  29-nov-17/Jorgen - Coded
!
      if (nprocx*nprocy>1) then
        call initiate_isendrcv_bdry_filter(f_og,filter_Hsize)
        call finalize_isendrcv_bdry_filter(f_Hlox,f_Hupx,f_Hloy,f_Hupy,filter_Hsize)
      else
        call boundconds_y_filter(f_og,f_Hloy,f_Hupy,filter_Hsize)
      endif
  endsubroutine communicate_filter_zones
!***********************************************************************
  subroutine pade_filter(f_og)
    use mpicomm, only: mpibarrier
    use Solid_cells_Mpicomm, only: cyclic_parallel_y,tridag_parallel_x
    use General, only: cyclic, tridag
!
!  high order padé filtering of solution array 
!  10th order on interor points, can choose 6th, 8th or 10th order at cylinder
!  boundary. 
!
!  Filtering will stabilize solution, necessary for certain boundary conditions
!  and grids (e.g., with large stretching)
!  Coefficients from Gaitonde & Visbal (2000)
!
!  WARNING: Only works for 2D serial runs at the moment
!
!  10-nov-17/Jorgen - Coded
!
    real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(inout) ::  f_og
!
    real, dimension(ny_ogrid,mfarray_ogrid) :: by
    real, dimension(nx_ogrid,mfarray_ogrid) :: bx
    real, dimension(ny_ogrid), save :: aWy, aPy, aEy
    real, dimension(nx_ogrid), save :: aWx, aPx, aEx
    integer :: i
    integer, save :: ii,jj
    real, save :: a0, a1, a2, a3, a4, a5
    real, save :: a0_6, a1_6, a2_6, a3_6
    logical :: lfirstcall = .true.
    if(lfirstcall) then
      a0=(193+126*af)/256.
      a1=(105+302*af)/256.
      a2=15*(-1+2*af)/64.
      a3=45*(1-2*af)/512.
      a4=5*(-1+2*af)/256.
      a5=(1-2*af)/512.
      a0_6=11./16.+5.*af/8.
      a1_6=15./32.+17.*af/16.
      a2_6=-3./16.+3.*af/8.
      a3_6=1./32.-af/16.

      aWy=af
      aPy=1.
      aEy=af

      aWx = af
      aPx = 1.
      aEx = af
!  Since we do not filter the point at the boundary
      if(ipx==0) then
        aWx(1) = 0.
        aPx(1) = 1.
        aEx(1) = 0.
      endif
!  Also, do not filter values in the interpolation region 
      if(ipx<nprocx-1) then
        jj=0
      else
        jj=interpol_filter
        aWx(nx_ogrid-jj) = 0.
        aPx(nx_ogrid-jj) = 1.
        aEx(nx_ogrid-jj) = 0.
      endif
      lfirstcall = .false.
    endif
! 
!  Filtering in theta direction
!  10th order on all points
!  Requires special handling near processor bondaries, due to not enough ghost points
!  Note! Only works for serial runs
!
    do i=l1_ogrid,l2_ogrid
      by(3:ny_ogrid-2,:) = a0*f_og(i,m1_ogrid+2:m2_ogrid-2,4,:) &
               + a1*0.5*(f_og(i,m1_ogrid+1:m2_ogrid-3,4,:) + f_og(i,m1_ogrid+3:m2_ogrid-1,4,:)) &
               + a2*0.5*(f_og(i,m1_ogrid  :m2_ogrid-4,4,:) + f_og(i,m1_ogrid+4:m2_ogrid  ,4,:)) &
               + a3*0.5*(f_og(i,m1_ogrid-1:m2_ogrid-5,4,:) + f_og(i,m1_ogrid+5:m2_ogrid+1,4,:)) &
               + a4*0.5*(f_og(i,m1_ogrid-2:m2_ogrid-6,4,:) + f_og(i,m1_ogrid+6:m2_ogrid+2,4,:)) &
               + a5*0.5*(f_og(i,m1_ogrid-3:m2_ogrid-7,4,:) + f_og(i,m1_ogrid+7:m2_ogrid+3,4,:))
!
!  Special handling outside of ordinary ghost points
!
      by(1,:) = a0*f_og(i,m1_ogrid,4,:) &
               + a1*0.5*(f_og(i,m1_ogrid-1,4,:)   + f_og(i,m1_ogrid+1,4,:)) &
               + a2*0.5*(f_og(i,m1_ogrid-2,4,:)   + f_og(i,m1_ogrid+2,4,:)) &
               + a3*0.5*(f_og(i,m1_ogrid-3,4,:)   + f_og(i,m1_ogrid+3,4,:)) &
               + a4*0.5*(f_filterH_lowery(i,2,1,:) + f_og(i,m1_ogrid+4,4,:)) &
               + a5*0.5*(f_filterH_lowery(i,1,1,:) + f_og(i,m1_ogrid+5,4,:))
      by(2,:) = a0*f_og(i,m1_ogrid+1,4,:) &
               + a1*0.5*(f_og(i,m1_ogrid  ,4,:)   + f_og(i,m1_ogrid+2,4,:)) &
               + a2*0.5*(f_og(i,m1_ogrid-1,4,:)   + f_og(i,m1_ogrid+3,4,:)) &
               + a3*0.5*(f_og(i,m1_ogrid-2,4,:)   + f_og(i,m1_ogrid+4,4,:)) &
               + a4*0.5*(f_og(i,m1_ogrid-3,4,:)   + f_og(i,m1_ogrid+5,4,:)) &
               + a5*0.5*(f_filterH_lowery(i,2,1,:) + f_og(i,m1_ogrid+6,4,:))
      by(ny_ogrid-1,:) = a0*f_og(i,m2_ogrid-1,4,:) &
               + a1*0.5*(f_og(i,m2_ogrid-2,4,:) + f_og(i,m2_ogrid  ,4,:)) &
               + a2*0.5*(f_og(i,m2_ogrid-3,4,:) + f_og(i,m2_ogrid+1,4,:)) &
               + a3*0.5*(f_og(i,m2_ogrid-4,4,:) + f_og(i,m2_ogrid+2,4,:)) &
               + a4*0.5*(f_og(i,m2_ogrid-5,4,:) + f_og(i,m2_ogrid+3,4,:)) &
               + a5*0.5*(f_og(i,m2_ogrid-6,4,:) + f_filterH_uppery(i,1,1,:))
      by(ny_ogrid  ,:) = a0*f_og(i,m2_ogrid,4,:) &
               + a1*0.5*(f_og(i,m2_ogrid-1,4,:) + f_og(i,m2_ogrid+1,4,:)) &
               + a2*0.5*(f_og(i,m2_ogrid-2,4,:) + f_og(i,m2_ogrid+2,4,:)) &
               + a3*0.5*(f_og(i,m2_ogrid-3,4,:) + f_og(i,m2_ogrid+3,4,:)) &
               + a4*0.5*(f_og(i,m2_ogrid-4,4,:) + f_filterH_uppery(i,1,1,:)) &
               + a5*0.5*(f_og(i,m2_ogrid-5,4,:) + f_filterH_uppery(i,2,1,:))
      if(nprocy>1) then
        if(.not. lfilter_rhoonly) then
          call cyclic_parallel_y(aWy,aPy,aEy,af,af,by(:,iux),f_og(i,m1_ogrid:m2_ogrid,4,iux),ny_ogrid)
          call cyclic_parallel_y(aWy,aPy,aEy,af,af,by(:,iuy),f_og(i,m1_ogrid:m2_ogrid,4,iuy),ny_ogrid)
          call cyclic_parallel_y(aWy,aPy,aEy,af,af,by(:,iuz),f_og(i,m1_ogrid:m2_ogrid,4,iuz),ny_ogrid)
        endif
        if(lfilter_TT) then
           call cyclic_parallel_y(aWy,aPy,aEy,af,af,by(:,iTT),f_og(i,m1_ogrid:m2_ogrid,4,iTT),ny_ogrid)
        endif
        call cyclic_parallel_y(aWy,aPy,aEy,af,af,by(:,irho),f_og(i,m1_ogrid:m2_ogrid,4,irho),ny_ogrid)
      else
        if(.not. lfilter_rhoonly) then
          call cyclic(aWy,aPy,aEy,af,af,by(:,iux),f_og(i,m1_ogrid:m2_ogrid,4,iux),ny_ogrid)
          call cyclic(aWy,aPy,aEy,af,af,by(:,iuy),f_og(i,m1_ogrid:m2_ogrid,4,iuy),ny_ogrid)
          call cyclic(aWy,aPy,aEy,af,af,by(:,iuz),f_og(i,m1_ogrid:m2_ogrid,4,iuz),ny_ogrid)
        endif
        if(lfilter_TT) then
           call cyclic(aWy,aPy,aEy,af,af,by(:,iTT),f_og(i,m1_ogrid:m2_ogrid,4,iTT),ny_ogrid)
        endif
        call cyclic(aWy,aPy,aEy,af,af,by(:,irho),f_og(i,m1_ogrid:m2_ogrid,4,irho),ny_ogrid)
      endif
    enddo
! 
!  Filtering in radial direction
!  10th order on interior points, 6th order near boundaries on the edge 
!  of cylindrical grid (near interpolation zone).
!  One-sided filter on cylinder surface can be set to 6th, 8th or 10th order
!  Surface point is not filtered, and neither are points in the 'filter'-zone between interpolations
!
    do i=m1_ogrid,m2_ogrid
      if(ipx==0) then
!
!  Special filtering near surface
!
        !call boundary_x_central(bx(1:5,:),f_og,af,i)
        !call boundary_x_6th(bx(1:5,:),f_og,af,i)
        !call boundary_x_8th(bx(1:5,:),f_og,af,i)
        !call boundary_x_10th(bx(1:5,:),f_og,af,i)
        call boundary_x_8_6th(bx(1:5,:),f_og,af,i)
        bx(6:nx_ogrid-2,:) = a0*f_og(l1_ogrid+5:l2_ogrid-2,i,4,:) &
                           + a1*0.5*(f_og(l1_ogrid+4:l2_ogrid-3,i,4,:) + f_og(l1_ogrid+ 6:l2_ogrid-1,i,4,:)) &
                           + a2*0.5*(f_og(l1_ogrid+3:l2_ogrid-4,i,4,:) + f_og(l1_ogrid+ 7:l2_ogrid  ,i,4,:)) &
                           + a3*0.5*(f_og(l1_ogrid+2:l2_ogrid-5,i,4,:) + f_og(l1_ogrid+ 8:l2_ogrid+1,i,4,:)) &
                           + a4*0.5*(f_og(l1_ogrid+1:l2_ogrid-6,i,4,:) + f_og(l1_ogrid+ 9:l2_ogrid+2,i,4,:)) &
                           + a5*0.5*(f_og(l1_ogrid  :l2_ogrid-7,i,4,:) + f_og(l1_ogrid+10:l2_ogrid+3,i,4,:))
      else
!
!  Special handling outside of ordinary ghost points
!
        bx(1,:) = a0*f_og(l1_ogrid,i,4,:) &
                 + a1*0.5*(f_og(l1_ogrid-1,i,4,:)   + f_og(l1_ogrid+1,i,4,:)) &
                 + a2*0.5*(f_og(l1_ogrid-2,i,4,:)   + f_og(l1_ogrid+2,i,4,:)) &
                 + a3*0.5*(f_og(l1_ogrid-3,i,4,:)   + f_og(l1_ogrid+3,i,4,:)) &
                 + a4*0.5*(f_filterH_lowerx(2,i,1,:) + f_og(l1_ogrid+4,i,4,:)) &
                 + a5*0.5*(f_filterH_lowerx(1,i,1,:) + f_og(l1_ogrid+5,i,4,:))
        bx(2,:) = a0*f_og(l1_ogrid+1,i,4,:) &
                 + a1*0.5*(f_og(l1_ogrid  ,i,4,:)    + f_og(l1_ogrid+2,i,4,:)) &
                 + a2*0.5*(f_og(l1_ogrid-1,i,4,:)    + f_og(l1_ogrid+3,i,4,:)) &
                 + a3*0.5*(f_og(l1_ogrid-2,i,4,:)    + f_og(l1_ogrid+4,i,4,:)) &
                 + a4*0.5*(f_og(l1_ogrid-3,i,4,:)    + f_og(l1_ogrid+5,i,4,:)) &
                 + a5*0.5*(f_filterH_lowerx(2,i,1,:) + f_og(l1_ogrid+6,i,4,:))
        bx(3:nx_ogrid-2,:) = a0*f_og(l1_ogrid+2:l2_ogrid-2,i,4,:) &
                 + a1*0.5*(f_og(l1_ogrid+1:l2_ogrid-3,i,4,:) + f_og(l1_ogrid+3:l2_ogrid-1,i,4,:)) &
                 + a2*0.5*(f_og(l1_ogrid  :l2_ogrid-4,i,4,:) + f_og(l1_ogrid+4:l2_ogrid  ,i,4,:)) &
                 + a3*0.5*(f_og(l1_ogrid-1:l2_ogrid-5,i,4,:) + f_og(l1_ogrid+5:l2_ogrid+1,i,4,:)) &
                 + a4*0.5*(f_og(l1_ogrid-2:l2_ogrid-6,i,4,:) + f_og(l1_ogrid+6:l2_ogrid+2,i,4,:)) &
                 + a5*0.5*(f_og(l1_ogrid-3:l2_ogrid-7,i,4,:) + f_og(l1_ogrid+7:l2_ogrid+3,i,4,:))
      endif
      if(ipx<nprocx-1) then
        bx(nx_ogrid-1,:) = a0*f_og(l2_ogrid-1,i,4,:) &
                    + a1*0.5*(f_og(l2_ogrid-2,i,4,:) + f_og(l2_ogrid  ,i,4,:)) &
                    + a2*0.5*(f_og(l2_ogrid-3,i,4,:) + f_og(l2_ogrid+1,i,4,:)) &
                    + a3*0.5*(f_og(l2_ogrid-4,i,4,:) + f_og(l2_ogrid+2,i,4,:)) &
                    + a4*0.5*(f_og(l2_ogrid-5,i,4,:) + f_og(l2_ogrid+3,i,4,:)) &
                    + a5*0.5*(f_og(l2_ogrid-6,i,4,:) + f_filterH_upperx(1,i,1,:))
        bx(nx_ogrid,:) = a0*f_og(l2_ogrid,i,4,:) &
                  + a1*0.5*(f_og(l2_ogrid-1,i,4,:) + f_og(l2_ogrid+1,i,4,:)) &
                  + a2*0.5*(f_og(l2_ogrid-2,i,4,:) + f_og(l2_ogrid+2,i,4,:)) &
                  + a3*0.5*(f_og(l2_ogrid-3,i,4,:) + f_og(l2_ogrid+3,i,4,:)) &
                  + a4*0.5*(f_og(l2_ogrid-4,i,4,:) + f_filterH_upperx(1,i,1,:)) &
                  + a5*0.5*(f_og(l2_ogrid-5,i,4,:) + f_filterH_upperx(2,i,1,:))
      else
        bx(nx_ogrid-1:nx_ogrid,:) = a0_6*f_og(l2_ogrid-1:l2_ogrid,i,4,:) &
                                  + a1_6*0.5*(f_og(l2_ogrid-2:l2_ogrid-1,i,4,:) + f_og(l2_ogrid  :l2_ogrid+1,i,4,:)) &
                                  + a2_6*0.5*(f_og(l2_ogrid-3:l2_ogrid-2,i,4,:) + f_og(l2_ogrid+1:l2_ogrid+2,i,4,:)) &
                                  + a3_6*0.5*(f_og(l2_ogrid-4:l2_ogrid-3,i,4,:) + f_og(l2_ogrid+2:l2_ogrid+3,i,4,:)) 
        bx(nx_ogrid-jj,:) = f_og(l2_ogrid-jj,i,4,:)
      endif
      if(nprocx>1) then
        if(.not. lfilter_rhoonly) then
          call tridag_parallel_x(aWx(1:nx_ogrid-jj),aPx(1:nx_ogrid-jj),aEx(1:nx_ogrid-jj),bx(1:nx_ogrid-jj,iux), &
            f_og(l1_ogrid:l2_ogrid-jj,i,4,iux), nx_ogrid-jj)
          call tridag_parallel_x(aWx(1:nx_ogrid-jj),aPx(1:nx_ogrid-jj),aEx(1:nx_ogrid-jj),bx(1:nx_ogrid-jj,iuy), &
            f_og(l1_ogrid:l2_ogrid-jj,i,4,iuy), nx_ogrid-jj)
          call tridag_parallel_x(aWx(1:nx_ogrid-jj),aPx(1:nx_ogrid-jj),aEx(1:nx_ogrid-jj),bx(1:nx_ogrid-jj,iuz), &
               f_og(l1_ogrid:l2_ogrid-jj,i,4,iuz), nx_ogrid-jj)
        endif
        if(lfilter_TT) then
           call tridag_parallel_x(aWx(1:nx_ogrid-jj),aPx(1:nx_ogrid-jj),aEx(1:nx_ogrid-jj),bx(1:nx_ogrid-jj,iTT), &
                f_og(l1_ogrid:l2_ogrid-jj,i,4,iTT), nx_ogrid-jj)
        endif
        call tridag_parallel_x(aWx(1:nx_ogrid-jj),aPx(1:nx_ogrid-jj),aEx(1:nx_ogrid-jj),bx(1:nx_ogrid-jj,irho), &
          f_og(l1_ogrid:l2_ogrid-jj,i,4,irho), nx_ogrid-jj)
      else
        if(.not. lfilter_rhoonly) then
          call tridag(aWx(1:nx_ogrid-jj),aPx(1:nx_ogrid-jj),aEx(1:nx_ogrid-jj), & 
            bx(1:nx_ogrid-jj,iux),f_og(l1_ogrid:l2_ogrid-jj,i,4,iux))
          call tridag(aWx(1:nx_ogrid-jj),aPx(1:nx_ogrid-jj),aEx(1:nx_ogrid-jj), & 
            bx(1:nx_ogrid-jj,iuy),f_og(l1_ogrid:l2_ogrid-jj,i,4,iuy))
          call tridag(aWx(1:nx_ogrid-jj),aPx(1:nx_ogrid-jj),aEx(1:nx_ogrid-jj), & 
            bx(1:nx_ogrid-jj,iuz),f_og(l1_ogrid:l2_ogrid-jj,i,4,iuz))
        endif
        if(lfilter_TT) then
           call tridag(aWx(1:nx_ogrid-jj),aPx(1:nx_ogrid-jj),aEx(1:nx_ogrid-jj), & 
                bx(1:nx_ogrid-jj,iTT),f_og(l1_ogrid:l2_ogrid-jj,i,4,iTT))
        endif
        call tridag(aWx(1:nx_ogrid-jj),aPx(1:nx_ogrid-jj),aEx(1:nx_ogrid-jj), &
          bx(1:nx_ogrid-jj,irho),f_og(l1_ogrid:l2_ogrid-jj,i,4,irho))
      endif
    enddo
!
  endsubroutine pade_filter
!***********************************************************************
  subroutine boundary_x_10th(bx_bound,f_og,af,i)
!
!  Compute the 10th order filter function for the radial boundary near the surface.
!  Use one-sided differences from Gaitonde & Visbal (2000)
!
!  10-nov-17/Jorgen - Coded
!
    real, dimension(5,mfarray_ogrid), intent(out) :: bx_bound
    real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(in)::  f_og
    real, intent(in) :: af
    integer, intent(in) :: i
    real, dimension(2:5,11), save :: aB_10
    integer :: j
    logical :: lfirstcall = .true.
!
    if(lfirstcall) then
      aB_10(2,1)  = (1   +1022*af)/1024.
      aB_10(2,2)  = (507 + 10*af)/512.
      aB_10(2,3)  = (45  +934*af)/1024.
      aB_10(2,4)  =  15*(-1+2*af)/128.
      aB_10(2,5)  = 105*( 1-2*af)/512.
      aB_10(2,6)  =  63*(-1+2*af)/256.
      aB_10(2,7)  = 105*( 1-2*af)/512.
      aB_10(2,8)  =  15*(-1+2*af)/128.
      aB_10(2,9)  =  45*( 1-2*af)/1024.
      aB_10(2,10) =   5*(-1+2*af)/512.
      aB_10(2,11) =     ( 1-2*af)/1024.

      aB_10(3,1)  =     (-1+2*af)/1024.
      aB_10(3,2)  = (5   +502*af)/512.
      aB_10(3,3)  = (979 + 90*af)/1024.
      aB_10(3,4)  = (15  + 98*af)/128.
      aB_10(3,5)  = 105*(-1+2*af)/512.
      aB_10(3,6)  =  63*( 1-2*af)/256.
      aB_10(3,7)  = 105*(-1+2*af)/512.
      aB_10(3,8)  =  15*( 1-2*af)/128.
      aB_10(3,9)  =  45*(-1+2*af)/1024.
      aB_10(3,10) =   5*( 1-2*af)/512.
      aB_10(3,11) =     (-1+2*af)/1024.

      aB_10(4,1)  =     ( 1-2*af)/1024.
      aB_10(4,2)  =   5*(-1+2*af)/512.
      aB_10(4,3)  = (45  +934*af)/1024.
      aB_10(4,4)  = (113 +30*af)/128.
      aB_10(4,5)  = (105 +302*af)/512.
      aB_10(4,6)  =  63*(-1+2*af)/256.
      aB_10(4,7)  = 105*( 1-2*af)/512.
      aB_10(4,8)  =  15*(-1+2*af)/128.
      aB_10(4,9)  =  45*( 1-2*af)/1024.
      aB_10(4,10) =   5*(-1+2*af)/512.
      aB_10(4,11) =     ( 1-2*af)/1024.

      aB_10(5,1)  =     (-1+2*af)/1024.
      aB_10(5,2)  =   5*( 1-2*af)/512.
      aB_10(5,3)  =  45*(-1+2*af)/1024.
      aB_10(5,4)  = (15  +98*af)/128.
      aB_10(5,5)  = (407 +210*af)/512.
      aB_10(5,6)  = (63  +130*af)/256.
      aB_10(5,7)  = 105*(-1+2*af)/512.
      aB_10(5,8)  =  15*( 1-2*af)/128.
      aB_10(5,9)  =  45*(-1+2*af)/1024.
      aB_10(5,10) =   5*( 1-2*af)/512.
      aB_10(5,11) =     (-1+2*af)/1024.

      lfirstcall = .false.
    endif
!
    do j=2,5
      bx_bound(j,:) = aB_10(j, 1)*f_og(l1_ogrid   ,i,4,:) + &
                      aB_10(j, 2)*f_og(l1_ogrid+ 1,i,4,:) + &
                      aB_10(j, 3)*f_og(l1_ogrid+ 2,i,4,:) + &
                      aB_10(j, 4)*f_og(l1_ogrid+ 3,i,4,:) + &
                      aB_10(j, 5)*f_og(l1_ogrid+ 4,i,4,:) + &
                      aB_10(j, 6)*f_og(l1_ogrid+ 5,i,4,:) + &
                      aB_10(j, 7)*f_og(l1_ogrid+ 6,i,4,:) + &
                      aB_10(j, 8)*f_og(l1_ogrid+ 7,i,4,:) + &
                      aB_10(j, 9)*f_og(l1_ogrid+ 8,i,4,:) + &
                      aB_10(j,10)*f_og(l1_ogrid+ 9,i,4,:) + &
                      aB_10(j,11)*f_og(l1_ogrid+10,i,4,:) 
    enddo
    bx_bound(1,:)=f_og(l1_ogrid,i,4,:)
!
  endsubroutine boundary_x_10th
!***********************************************************************
  subroutine boundary_x_8th(bx_bound,f_og,af,i)
!
!  Compute the 8th order filter function for the radial boundary near the surface.
!  Use central differences point i=5 (with i=1 at the boundary) and 
!  one-sided differences fro points i=2:4 
!  Weights from Gaitonde & Visbal (2000)
!
!  10-nov-17/Jorgen - Coded
!
    real, dimension(5,mfarray_ogrid), intent(out) :: bx_bound
    real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(in)::  f_og
    real, intent(in) :: af
    integer, intent(in) :: i
    real, dimension(2:4,9), save :: aB_8
    real, save :: a0_8, a1_8, a2_8, a3_8, a4_8
    integer :: j
    logical :: lfirstcall = .true.
!
    if(lfirstcall) then
      a0_8=(93+70*af)/128.
      a1_8=(7+18*af)/16.
      a2_8=(-7+14*af)/32.
      a3_8=1/16.-af/8.
      a4_8=-1/128.+af/64.
      aB_8(4,1)  =  1/256. -   af/128.
      aB_8(4,2)  = -1/32.  +   af/16.
      aB_8(4,3)  =  7/64.  +25*af/32.
      aB_8(4,4)  = 25/32.  + 7*af/16.
      aB_8(4,5)  = 35/128. +29*af/64.
      aB_8(4,6)  = -7/32.  + 7*af/16.
      aB_8(4,7)  =  7/64.  - 7*af/32.
      aB_8(4,8)  = -1/32.  +   af/16.
      aB_8(4,9)  =  1/256. -   af/128.

      aB_8(3,1)  = -1/256. +   af/128.
      aB_8(3,2)  =  1/32.  +15*af/16.
      aB_8(3,3)  = 57/64.  + 7*af/32.
      aB_8(3,4)  =  7/32.  + 9*af/16.
      aB_8(3,5)  = 7*(-5 + 10*af)/128.
      aB_8(3,6)  =  7/32.  - 7*af/16.
      aB_8(3,7)  = -7/64.  + 7*af/32.
      aB_8(3,8)  =  1/32.  -   af/16.
      aB_8(3,9)  = -1/256. +   af/128.
      
      aB_8(2,1)  =  1/256. +127*af/128.
      aB_8(2,2)  = 31/32.  +   af/16.
      aB_8(2,3)  =  7/64.  +25*af/32.
      aB_8(2,4)  = -7/32.  + 7*af/16.
      aB_8(2,5)  = 7*( 5 - 10*af)/128.
      aB_8(2,6)  = -7/32.  + 7*af/16.
      aB_8(2,7)  =  7/64.  - 7*af/32.
      aB_8(2,8)  = -1/32.  +   af/16.
      aB_8(2,9)  =  1/256. -   af/128.

      lfirstcall = .false.
    endif

    bx_bound(5,:) = a0_8*f_og(l1_ogrid+4,i,4,:) &
                  + a1_8*0.5*(f_og(l1_ogrid+3,i,4,:) + f_og(l1_ogrid+5,i,4,:)) &
                  + a2_8*0.5*(f_og(l1_ogrid+2,i,4,:) + f_og(l1_ogrid+6,i,4,:)) &
                  + a3_8*0.5*(f_og(l1_ogrid+1,i,4,:) + f_og(l1_ogrid+7,i,4,:)) &
                  + a4_8*0.5*(f_og(l1_ogrid  ,i,4,:) + f_og(l1_ogrid+8,i,4,:)) 
    do j=2,4
      bx_bound(j,:) = aB_8(j,1)*f_og(l1_ogrid  ,i,4,:) + &
                      aB_8(j,2)*f_og(l1_ogrid+1,i,4,:) + &
                      aB_8(j,3)*f_og(l1_ogrid+2,i,4,:) + &
                      aB_8(j,4)*f_og(l1_ogrid+3,i,4,:) + &
                      aB_8(j,5)*f_og(l1_ogrid+4,i,4,:) + &
                      aB_8(j,6)*f_og(l1_ogrid+5,i,4,:) + &
                      aB_8(j,7)*f_og(l1_ogrid+6,i,4,:) + &
                      aB_8(j,8)*f_og(l1_ogrid+7,i,4,:) + &
                      aB_8(j,9)*f_og(l1_ogrid+8,i,4,:)
    enddo
    bx_bound(1,:)=f_og(l1_ogrid,i,4,:)

  endsubroutine boundary_x_8th
!***********************************************************************
  subroutine boundary_x_6th(bx_bound,f_og,af,i)
!
!  Compute the 6th order filter function for the radial boundary near the surface.
!  Use central differences for points i=3:5 (with i=1 at the boundary) and 
!  one-sided differences fro points i=2:3 
!  Weights from Gaitonde & Visbal (2000)
!
!  10-nov-17/Jorgen - Coded
!
    real, dimension(5,mfarray_ogrid), intent(out) :: bx_bound
    real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(in)::  f_og
    real, intent(in) :: af
    integer, intent(in) :: i
    real, dimension(2:3,7), save :: aB_6
    real, save :: a0_6, a1_6, a2_6, a3_6
    logical :: lfirstcall = .true.
!
    if(lfirstcall) then
      a0_6=11./16.+5.*af/8.
      a1_6=15./32.+17.*af/16.
      a2_6=-3./16.+3.*af/8.
      a3_6=1./32.-af/16.

      aB_6(2,1)  =  1/64. +31*af/32.
      aB_6(2,2)  = 29/32. + 3*af/16.
      aB_6(2,3)  = 15/64. +17*af/32.
      aB_6(2,4)  = -5/16. + 5*af/8.
      aB_6(2,5)  = 15/64. -15*af/32.
      aB_6(2,6)  = -3/32. + 3*af/16.
      aB_6(2,7)  =  1/64. -   af/32.

      aB_6(3,1)  = -1/64. +   af/32.
      aB_6(3,2)  =  3/32. +13*af/16.
      aB_6(3,3)  = 49/64. +15*af/32.
      aB_6(3,4)  =  5/16. + 3*af/8.
      aB_6(3,5)  =-15/64. +15*af/32.
      aB_6(3,6)  =  3/32. - 3*af/16.
      aB_6(3,7)  = -1/64. +   af/32.

      lfirstcall = .false.
    endif
    bx_bound(4:5,:) = a0_6*f_og(l1_ogrid+3:l1_ogrid+4,i,4,:) &
                    + a1_6*0.5*(f_og(l1_ogrid+2:l1_ogrid+3,i,4,:) + f_og(l1_ogrid+4:l1_ogrid+5,i,4,:)) &
                    + a2_6*0.5*(f_og(l1_ogrid+1:l1_ogrid+2,i,4,:) + f_og(l1_ogrid+5:l1_ogrid+6,i,4,:)) &
                    + a3_6*0.5*(f_og(l1_ogrid  :l1_ogrid+1,i,4,:) + f_og(l1_ogrid+6:l1_ogrid+7,i,4,:)) 
    bx_bound(3,:) = aB_6(3,1)*f_og(l1_ogrid  ,i,4,:) + &
                    aB_6(3,2)*f_og(l1_ogrid+1,i,4,:) + &
                    aB_6(3,3)*f_og(l1_ogrid+2,i,4,:) + &
                    aB_6(3,4)*f_og(l1_ogrid+3,i,4,:) + &
                    aB_6(3,5)*f_og(l1_ogrid+4,i,4,:) + &
                    aB_6(3,6)*f_og(l1_ogrid+5,i,4,:) + &
                    aB_6(3,7)*f_og(l1_ogrid+6,i,4,:)
    bx_bound(2,:) = aB_6(2,1)*f_og(l1_ogrid  ,i,4,:) + &
                    aB_6(2,2)*f_og(l1_ogrid+1,i,4,:) + &
                    aB_6(2,3)*f_og(l1_ogrid+2,i,4,:) + &
                    aB_6(2,4)*f_og(l1_ogrid+3,i,4,:) + &
                    aB_6(2,5)*f_og(l1_ogrid+4,i,4,:) + &
                    aB_6(2,6)*f_og(l1_ogrid+5,i,4,:) + &
                    aB_6(2,7)*f_og(l1_ogrid+6,i,4,:)
    bx_bound(1,:)=f_og(l1_ogrid,i,4,:)
!
  endsubroutine boundary_x_6th
!***********************************************************************
  subroutine boundary_x_central(bx_bound,f_og,af,i)
!
!  Compute the 6th order filter function for the radial boundary near the surface.
!  Use central differences for points i=3:5 (with i=1 at the boundary) and 
!  one-sided differences fro points i=2:3 
!  Weights from Gaitonde & Visbal (2000)
!
!  10-nov-17/Jorgen - Coded
!
    real, dimension(5,mfarray_ogrid), intent(out) :: bx_bound
    real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(in)::  f_og
    real, intent(in) :: af
    integer, intent(in) :: i
    real, save :: a0_8, a1_8, a2_8, a3_8, a4_8
    real, save :: a0_6, a1_6, a2_6, a3_6
    real, save :: a0_4, a1_4, a2_4
    real, save :: a0_2, a1_2
    logical :: lfirstcall = .true.
!
    if(lfirstcall) then
      a0_8=(93+70*af)/128.
      a1_8=(7+18*af)/16.
      a2_8=(-7+14*af)/32.
      a3_8=1/16.-af/8.
      a4_8=-1/128.+af/64.
!
      a0_6=11./16.+5.*af/8.
      a1_6=15./32.+17.*af/16.
      a2_6=-3./16.+3.*af/8.
      a3_6=1./32.-af/16.

      a0_4=5/8. + 3*af/4.
      a1_4=1/2. + af
      a2_4=-1/8. + af/4.

      a0_2=0.5+af
      a1_2=0.5+af
      lfirstcall = .false.
    endif

    bx_bound(5,:) = a0_8*     f_og(l1_ogrid+4,i,4,:) &
                  + a1_8*0.5*(f_og(l1_ogrid+3,i,4,:) + f_og(l1_ogrid+5,i,4,:)) &
                  + a2_8*0.5*(f_og(l1_ogrid+2,i,4,:) + f_og(l1_ogrid+6,i,4,:)) &
                  + a3_8*0.5*(f_og(l1_ogrid+1,i,4,:) + f_og(l1_ogrid+7,i,4,:)) &
                  + a4_8*0.5*(f_og(l1_ogrid  ,i,4,:) + f_og(l1_ogrid+8,i,4,:)) 
    bx_bound(4,:) = a0_6*     f_og(l1_ogrid+3,i,4,:) &
                  + a1_6*0.5*(f_og(l1_ogrid+2,i,4,:) + f_og(l1_ogrid+4,i,4,:)) &
                  + a2_6*0.5*(f_og(l1_ogrid+1,i,4,:) + f_og(l1_ogrid+5,i,4,:)) &
                  + a3_6*0.5*(f_og(l1_ogrid  ,i,4,:) + f_og(l1_ogrid+6,i,4,:)) 
    bx_bound(3,:) = a0_4*     f_og(l1_ogrid+2,i,4,:) &
                  + a1_4*0.5*(f_og(l1_ogrid+1,i,4,:) + f_og(l1_ogrid+3,i,4,:)) &
                  + a1_4*0.5*(f_og(l1_ogrid  ,i,4,:) + f_og(l1_ogrid+4,i,4,:))
    bx_bound(2,:) = a0_2*     f_og(l1_ogrid+1,i,4,:) &
                  + a1_2*0.5*(f_og(l1_ogrid  ,i,4,:) + f_og(l1_ogrid+1,i,4,:))
    bx_bound(1,:)=f_og(l1_ogrid,i,4,:)
!
  endsubroutine boundary_x_central
!***********************************************************************
  subroutine boundary_x_8_6th(bx_bound,f_og,af,i)
!
!  Compute the 8th order filter function for the radial boundary near the surface.
!  Use central differences point i=5 (with i=1 at the boundary) and 
!  one-sided differences fro points i=2:4 
!  Weights from Gaitonde & Visbal (2000)
!
!  10-nov-17/Jorgen - Coded
!
    real, dimension(5,mfarray_ogrid), intent(out) :: bx_bound
    real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mfarray_ogrid), intent(in)::  f_og
    real, intent(in) :: af
    integer, intent(in) :: i
    real, dimension(3:4,9), save :: aB_8
    real, dimension(7), save :: aB_6
    real, save :: a0_8, a1_8, a2_8, a3_8, a4_8
    logical :: lfirstcall = .true.
    integer :: j
!
    if(lfirstcall) then
      a0_8=(93+70*af)/128.
      a1_8=(7+18*af)/16.
      a2_8=(-7+14*af)/32.
      a3_8=1/16.-af/8.
      a4_8=-1/128.+af/64.
      aB_8(4,1)  =  1/256. -   af/128.
      aB_8(4,2)  = -1/32.  +   af/16.
      aB_8(4,3)  =  7/64.  +25*af/32.
      aB_8(4,4)  = 25/32.  + 7*af/16.
      aB_8(4,5)  = 35/128. +29*af/64.
      aB_8(4,6)  = -7/32.  + 7*af/16.
      aB_8(4,7)  =  7/64.  - 7*af/32.
      aB_8(4,8)  = -1/32.  +   af/16.
      aB_8(4,9)  =  1/256. -   af/128.

      aB_8(3,1)  = -1/256. +   af/128.
      aB_8(3,2)  =  1/32.  +15*af/16.
      aB_8(3,3)  = 57/64.  + 7*af/32.
      aB_8(3,4)  =  7/32.  + 9*af/16.
      aB_8(3,5)  = 7*(-5 + 10*af)/128.
      aB_8(3,6)  =  7/32.  - 7*af/16.
      aB_8(3,7)  = -7/64.  + 7*af/32.
      aB_8(3,8)  =  1/32.  -   af/16.
      aB_8(3,9)  = -1/256. +   af/128.

      aB_6(1)  =  1/64. +31*af/32.
      aB_6(2)  = 29/32. + 3*af/16.
      aB_6(3)  = 15/64. +17*af/32.
      aB_6(4)  = -5/16. + 5*af/8.
      aB_6(5)  = 15/64. -15*af/32.
      aB_6(6)  = -3/32. + 3*af/16.
      aB_6(7)  =  1/64. -   af/32.

      lfirstcall = .false.
    endif

    bx_bound(5,:) = a0_8*f_og(l1_ogrid+4,i,4,:) &
                  + a1_8*0.5*(f_og(l1_ogrid+3,i,4,:) + f_og(l1_ogrid+5,i,4,:)) &
                  + a2_8*0.5*(f_og(l1_ogrid+2,i,4,:) + f_og(l1_ogrid+6,i,4,:)) &
                  + a3_8*0.5*(f_og(l1_ogrid+1,i,4,:) + f_og(l1_ogrid+7,i,4,:)) &
                  + a4_8*0.5*(f_og(l1_ogrid  ,i,4,:) + f_og(l1_ogrid+8,i,4,:)) 
    do j=3,4
      bx_bound(j,:) = aB_8(j,1)*f_og(l1_ogrid  ,i,4,:) + &
                      aB_8(j,2)*f_og(l1_ogrid+1,i,4,:) + &
                      aB_8(j,3)*f_og(l1_ogrid+2,i,4,:) + &
                      aB_8(j,4)*f_og(l1_ogrid+3,i,4,:) + &
                      aB_8(j,5)*f_og(l1_ogrid+4,i,4,:) + &
                      aB_8(j,6)*f_og(l1_ogrid+5,i,4,:) + &
                      aB_8(j,7)*f_og(l1_ogrid+6,i,4,:) + &
                      aB_8(j,8)*f_og(l1_ogrid+7,i,4,:) + &
                      aB_8(j,9)*f_og(l1_ogrid+8,i,4,:)
    enddo
    bx_bound(2,:) = aB_6(1)*f_og(l1_ogrid  ,i,4,:) + &
                    aB_6(2)*f_og(l1_ogrid+1,i,4,:) + &
                    aB_6(3)*f_og(l1_ogrid+2,i,4,:) + &
                    aB_6(4)*f_og(l1_ogrid+3,i,4,:) + &
                    aB_6(5)*f_og(l1_ogrid+4,i,4,:) + &
                    aB_6(6)*f_og(l1_ogrid+5,i,4,:) + &
                    aB_6(7)*f_og(l1_ogrid+6,i,4,:)
    bx_bound(1,:)=f_og(l1_ogrid,i,4,:)

  endsubroutine boundary_x_8_6th
  !***********************************************************************
  subroutine check_cyl_pos(cyl_pos,domstart,domend)
!
! Check if the cylinder is positioned inside the domain
!
    real, dimension(3) :: cyl_pos, domstart, domend

    if ((domstart(1) < cyl_pos(1) .and. cyl_pos(1) < domend(1)) .and. &
         (domstart(2) < cyl_pos(2) .and. cyl_pos(2) < domend(2)) .and. & 
         (domstart(3) < cyl_pos(3) .and. cyl_pos(3) < domend(3))) then
       continue
    else
       call fatal_error('init_solid_cells','Cylinder placed outside domain')
    endif
    
  endsubroutine check_cyl_pos
  !***********************************************************************
  subroutine create_curv_cart_transform(trans_mat)

    real, dimension(my_ogrid,2), intent(out) :: trans_mat

    trans_mat(:,1) = sin(y_ogrid(:))
    trans_mat(:,2) = cos(y_ogrid(:))
    
  endsubroutine create_curv_cart_transform
!***********************************************************************
!
! CURRENTLY NOT USED
!
!***********************************************************************
!  logical function linear_interpolate_cart_HO(farr,ivar1,ivar2,xxp,inear_glob,fp,lcheck,order)
!!
!!  Interpolate the value of f to arbitrary (xp, yp) CARTESIAN coordinate
!!  using the high-order lagrangian interpolation.
!! 
!!  TODO: Extend to 3D
!!  TODO: Extend to arbitrary order
!!
!!  The coefficients are determined by the 2xN grid points surrounding the
!!  interpolation point.
!!  Global coordinates are used for the interpolation, to allow interpolation of 
!!  values outside this processors domain.
!!
!!  26-apr-17/Jorgen: Adapted from linear_interpolate_curv_HO
!!
!      integer :: ivar1, ivar2
!      integer, intent(in) :: order
!      real, dimension (3) :: xxp
!      real, dimension (order,order,2,ivar2-ivar1+1) :: farr
!      real, dimension (ivar2-ivar1+1) :: fp
!      integer, dimension (3) :: inear_glob
!      logical :: lcheck
!!
!      intent(in)  :: farr, ivar1, ivar2, xxp, inear_glob, lcheck
!      intent(out) :: fp
!      integer :: i,ix0,iy0,iz0
!      real, dimension(ivar2-ivar1+1) :: g1,g2
!      real :: xp,yp,l1,l2,l3,l4
!      real, dimension(order,ivar2-ivar1+1) :: g_interp
!      real, dimension(order) :: lagrange
!      real, dimension(order,order) :: dx1,dy1
!      integer :: j,k,l
!!
!!  Determine index value of lowest lying corner point of grid box surrounding
!!  the interpolation point.
!!
!      linear_interpolate_cart_HO= .true.
!!
!      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3)
!!
!!  Check if the grid point interval is really correct.
!!
!      if ((xglobal(ix0)<=xxp(1) .and. xglobal(ix0+1)>=xxp(1) .or. nxgrid==1) .and. &
!          (yglobal(iy0)<=xxp(2) .and. yglobal(iy0+1)>=xxp(2) .or. nygrid==1) .and. &
!          (zglobal(iz0)<=xxp(3) .and. zglobal(iz0+1)>=xxp(3) .or. nzgrid==1)) then
!        ! Everything okay
!      else
!        print*, 'linear_interpolate_cartesian: Global interpolation point does not ' // &
!            'lie within the calculated grid point interval.'
!        print*, 'iproc = ', iproc_world
!        print*, 'mxgrid, xglobal(1), xglobal(mx) = ', mxgrid, xglobal(1), xglobal(mxgrid)
!        print*, 'mygrid, yglobal(1), yglobal(my) = ', mygrid, yglobal(1), yglobal(mygrid)
!        print*, 'mzgrid, zglobal(1), zglobal(mz) = ', mzgrid, zglobal(1), zglobal(mzgrid)
!        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
!        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal(ix0), xglobal(ix0+1)
!        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal(iy0), yglobal(iy0+1)
!        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal(iz0), zglobal(iz0+1)
!        linear_interpolate_cart_HO = .false.
!        return
!      endif
!!
!!  Set up 1D Lagrange basis polynomials in x-direction
!! 
!      k=-floor(order/2.)
!      do i=1,order
!        k=k+1
!        l=-floor(order/2.)
!        do j=1,order
!          l=l+1
!          dx1(i,j)=xglobal(ix0+k)-xglobal(ix0+l)
!        enddo
!        dx1(i,i)=1 ! To avoid division by zero
!      enddo
!      dx1=1./dx1
!      xp=xxp(1)
!      k=-floor(order/2.)
!      do i=1,order
!        k=k+1
!        lagrange(i)=1./(xp-xglobal(ix0+k))
!      enddo
!      k=-floor(order/2.)
!      do i=1,order
!        k=k+1
!        lagrange(:)=lagrange(:)*(xp-xglobal(ix0+k))*dx1(:,i)
!      enddo
!      g_interp=0
!      do i=1,order
!        g_interp(:,ivar1:ivar2)=g_interp(:,ivar1:ivar2)+farr(i,:,1,ivar1:ivar2)*lagrange(i)
!      enddo
!      ! y-dir
!      k=-floor(order/2.)
!      do i=1,order
!        k=k+1
!        l=-floor(order/2.)
!        do j=1,order
!          l=l+1
!          dy1(i,j)=yglobal(iy0+k)-yglobal(iy0+l)
!        enddo
!        dy1(i,i)=1 ! To avoid division by zero
!      enddo
!      dy1=1./dy1
!      yp=xxp(2)
!      k=-floor(order/2.)
!      do i=1,order
!        k=k+1
!        lagrange(i)=1./(yp-yglobal(iy0+k))
!      enddo
!      k=-floor(order/2.)
!      do i=1,order
!        k=k+1
!        lagrange(:)=lagrange(:)*(yp-yglobal(iy0+k))*dy1(:,i)
!      enddo
!      fp=0
!      do i=1,order
!        fp(ivar1:ivar2)=fp(ivar1:ivar2)+g_interp(i,ivar1:ivar2)*lagrange(i)
!      enddo
!!
!!  Do a reality check on the interpolation scheme.
!!
!      if (lcheck) then
!        do i=1,ivar2-ivar1+1
!          if (fp(i)>maxval(farr(:,:,1,i)).and.i/=3) then 
!           l1=(xp-xglobal(ix0+1))/(xglobal(ix0)-xglobal(ix0+1))
!           l2=(xp-xglobal(ix0  ))/(xglobal(ix0+1)-xglobal(ix0))
!           g1=farr(k,k,1,ivar1:ivar2)*l1+farr(k+1,k,1,ivar1:ivar2)*l2
!           g2=farr(k,k+1,1,ivar1:ivar2)*l1+farr(k,k+1,1,ivar1:ivar2)*l2
!           l3=(yp-yglobal(iy0+1))/(yglobal(iy0)-yglobal(iy0+1))
!           l4=(yp-yglobal(iy0  ))/(yglobal(iy0+1)-yglobal(iy0))
!           fp=g1*l3+g2*l4
!          elseif (fp(i)<minval(farr(:,:,1,i)).and.i/=3) then
!           l1=(xp-xglobal(ix0+1))/(xglobal(ix0)-xglobal(ix0+1))
!           l2=(xp-xglobal(ix0  ))/(xglobal(ix0+1)-xglobal(ix0))
!           g1=farr(k,k,1,ivar1:ivar2)*l1+farr(k+1,k,1,ivar1:ivar2)*l2
!           g2=farr(k,k+1,1,ivar1:ivar2)*l1+farr(k,k+1,1,ivar1:ivar2)*l2
!           l3=(yp-yglobal(iy0+1))/(yglobal(iy0)-yglobal(iy0+1))
!           l4=(yp-yglobal(iy0  ))/(yglobal(iy0+1)-yglobal(iy0))
!           fp=g1*l3+g2*l4
!          endif
!          if ((fp(i)>maxval(farr(:,:,1,i)).and.i/=3) .or. (fp(i)<minval(farr(:,:,1,i)).and.i/=3)) then
!            print*, 'linear_interpolate_cart_HO: interpolated value is smaller or larger than'
!            print*, 'linear_interpolate_cart_HO: a values at the corner points, even after linearization!'
!            print*, '------------------'
!            linear_interpolate_cart_HO=.false.
!          endif
!          if (fp(i)/=fp(i)) then
!            print*, 'linear_interpolate_cart_HO: interpolated value is NaN'
!            print*, 'linear_interpolate_cart_HO: xxp=', xxp
!            print*, 'linear_interpolate_cart_HO: x0, y0, z0=', &
!                xglobal(ix0), yglobal(iy0), zglobal(iz0)
!            print*, 'linear_interpolate_cart_HO: i, fp(i)=', i, fp(i)
!            print*, 'linear_interpolate_cart_HO: farr=', farr(:,:,1,i)
!            print*, '------------------'
!            linear_interpolate_cart_HO=.false.
!          endif
!        enddo
!      endif
!!
!  endfunction linear_interpolate_cart_HO
!
!***********************************************************************
!
!!     subroutine send_rcv_all_data(ivar1,ivar2,f_cartesian)
!! 
!! !   Subroutine that exhanges all data in f-arrays, both for curvilinear and cartesian grid,
!! !   between all processors. 
!! !   Very inefficient, but can be useful for testing.
!! !   Only works properly in 2D.
!! !
!! !   30-sep-17/Jorgen: Coded
!! 
!!       use mpicomm, only: mpisend_int, mpisend_real, mpirecv_int, mpirecv_real, mpibcast_real
!!       real, dimension (mx,my,mz,mfarray),intent(in) :: f_cartesian
!!       integer, intent(in) :: ivar1,ivar2
!!       real, dimension (nx_ogrid, ny_ogrid, nz_ogrid,ivar2-ivar1+1) :: fbuf_og
!!       real, dimension (nx,       ny,       nz,      ivar2-ivar1+1) :: fbuf_cg
!!       real, dimension (nxgrid_ogrid, nygrid_ogrid, nzgrid_ogrid,ivar2-ivar1+1) ::  fgrid_ogrid_tmp
!!       real, dimension (nxgrid, nygrid, nzgrid,ivar2-ivar1+1) ::  fgrid_cartesian_tmp
!!       integer, dimension(4) :: nfbuf_og
!!       integer, dimension(4) :: nfbuf_cg
!!       integer, dimension(3) :: ipxyz
!!       integer :: ixdo,ixup,iydo,iyup,izdo,izup
!!       integer :: jx, iproc_recv
!!       integer :: i,j
!! 
!!       nfbuf_og= (/ nxgrid_ogrid, nygrid_ogrid, nzgrid_ogrid,ivar2-ivar1+1/)
!!       nfbuf_cg= (/ nxgrid, nygrid, nzgrid,ivar2-ivar1+1/)
!! !
!!     if (iproc/=root) then
!!       fbuf_og = f_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2)
!!       fbuf_cg = f_cartesian(l1:l2,m1:m2,n1:n2,ivar1:ivar2)
!!       call mpisend_int((/ipx,ipy,ipz/),3,root,111)
!!       do i = 1,nx_ogrid
!!         do j = 1,ny_ogrid
!!           call mpisend_real(fbuf_og(i,j,:,:),nfbuf_og(3:4),root,i*ny_ogrid+j)
!!         enddo
!!       enddo
!!       do i = 1,nx
!!         do j = 1,ny
!!           call mpisend_real(fbuf_cg(i,j,:,:),nfbuf_cg(3:4),root,i*ny_ogrid+j)
!!         enddo
!!       enddo
!!     else
!! !
!! !  The root processor, in turn, receives the data from the others
!! !
!!       do jx=0,ncpus-1
!!         !avoid send-to-self
!!         if (jx/=root) then
!! !
!! !  Formula of the serial processor number:
!! !  iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
!! !  Since for the x-row ipy=ipz=0, this reduces
!! !  to iproc_recv=jx.
!! !
!!           iproc_recv=jx
!!           call mpirecv_int(ipxyz,3,iproc_recv,111)
!!           do i = 1,nx_ogrid
!!             do j = 1,ny_ogrid
!!               call mpirecv_real(fbuf_og(i,j,:,:),nfbuf_og(3:4),iproc_recv,i*ny_ogrid+j)
!!             enddo
!!           enddo
!!           !call mpirecv_real(fbuf_og,nfbuf_og,iproc_recv,111)
!!           do i = 1,nx
!!             do j = 1,ny
!!               call mpirecv_real(fbuf_cg(i,j,:,:),nfbuf_cg(3:4),iproc_recv,i*ny_ogrid+j)
!!             enddo
!!           enddo
!!           !call mpirecv_real(fbuf_cg,nfbuf_cg,iproc_recv,112)
!! !
!!           ixdo=ipxyz(1)*nx_ogrid+1
!!           ixup=(ipxyz(1)+1)*nx_ogrid
!!           iydo=ipxyz(2)*ny_ogrid+1
!!           iyup=(ipxyz(2)+1)*ny_ogrid
!!           izdo=ipxyz(3)*nz_ogrid+1
!!           izup=(ipxyz(3)+1)*nz_ogrid
!! 
!!       fgrid_ogrid_tmp    (ixdo:ixup,iydo:iyup,izdo:izup,ivar1:ivar2) = fbuf_og
!! 
!!           ixdo=ipxyz(1)*nx+1
!!           ixup=(ipxyz(1)+1)*nx
!!           iydo=ipxyz(2)*ny+1
!!           iyup=(ipxyz(2)+1)*ny
!!           izdo=ipxyz(3)*nz+1
!!           izup=(ipxyz(3)+1)*nz
!!       fgrid_cartesian_tmp(ixdo:ixup,iydo:iyup,izdo:izup,ivar1:ivar2) = fbuf_cg
!!       
!!         else
!!       fgrid_ogrid_tmp(1:nx_ogrid,1:ny_ogrid,1:nz_ogrid,ivar1:ivar2) = &
!!           f_ogrid(l1_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2)
!!       fgrid_cartesian_tmp(1:nx,1:ny,1:nz,ivar1:ivar2) = &
!!           f_cartesian(l1:l2,m1:m2,n1:n2,ivar1:ivar2)
!!         endif
!!       enddo
!!     endif
!!     if(iproc==root) then
!!       fgrid_ogrid(nghost+1:mxgrid_ogrid-nghost,nghost+1:mygrid_ogrid-nghost, &
!!         nghost+1:mzgrid_ogrid-nghost,ivar1:ivar2) = fgrid_ogrid_tmp
!!       fgrid_cartesian(nghost+1:mxgrid-nghost,nghost+1:mygrid-nghost, &
!!         nghost+1:mzgrid-nghost,ivar1:ivar2) = fgrid_cartesian_tmp
!! ! Ghosts cells x-direction      
!!       fgrid_ogrid(1:nghost,nghost+1:mygrid_ogrid-nghost, &
!!         nghost+1:mzgrid_ogrid-nghost,ivar1:ivar2) = &
!!         fgrid_ogrid_tmp(nxgrid_ogrid-nghost+1:nxgrid_ogrid,1:nygrid_ogrid, &
!!         1:nzgrid_ogrid,ivar1:ivar2) 
!!       fgrid_ogrid(mxgrid_ogrid-nghost+1:mxgrid_ogrid,nghost+1:mygrid_ogrid-nghost, &
!!         nghost+1:mzgrid_ogrid-nghost,ivar1:ivar2) = &
!!         fgrid_ogrid_tmp(1:nghost,1:nygrid_ogrid, &
!!         1:nzgrid_ogrid,ivar1:ivar2) 
!!       fgrid_cartesian(1:nghost,nghost+1:mygrid-nghost, &
!!         nghost+1:mzgrid-nghost,ivar1:ivar2) = &
!!         fgrid_cartesian_tmp(nxgrid-nghost+1:nxgrid,1:nygrid, &
!!         1:nzgrid,ivar1:ivar2) 
!!       fgrid_cartesian(mxgrid-nghost+1:mxgrid,nghost+1:mygrid-nghost, &
!!         nghost+1:mzgrid-nghost,ivar1:ivar2) = &
!!         fgrid_cartesian_tmp(1:nghost,1:nygrid, &
!!         1:nzgrid,ivar1:ivar2) 
!! ! Ghosts cells y-direction      
!!       fgrid_ogrid(nghost+1:mxgrid_ogrid-nghost,1:nghost, &
!!         nghost+1:mzgrid_ogrid-nghost,ivar1:ivar2) = &
!!         fgrid_ogrid_tmp(1:nxgrid_ogrid,nygrid_ogrid-nghost+1:nygrid_ogrid, &
!!         1:nzgrid_ogrid,ivar1:ivar2) 
!!       fgrid_ogrid(nghost+1:mxgrid_ogrid-nghost,mygrid_ogrid-nghost+1:mygrid_ogrid, &
!!         nghost+1:mzgrid_ogrid-nghost,ivar1:ivar2) = &
!!         fgrid_ogrid_tmp(1:nxgrid_ogrid,1:nghost, &
!!         1:nzgrid_ogrid,ivar1:ivar2) 
!!       fgrid_cartesian(nghost+1:mxgrid-nghost,1:nghost, &
!!         nghost+1:mzgrid-nghost,ivar1:ivar2) = &
!!         fgrid_cartesian_tmp(1:nxgrid,nygrid-nghost+1:nygrid, &
!!         1:nzgrid,ivar1:ivar2) 
!!       fgrid_cartesian(nghost+1:mxgrid-nghost,mygrid-nghost+1:mygrid, &
!!         nghost+1:mzgrid-nghost,ivar1:ivar2) = &
!!         fgrid_cartesian_tmp(1:nxgrid,1:nghost, &
!!         1:nzgrid,ivar1:ivar2) 
!! ! Ghosts cells z-direction (2D runs!)
!!       fgrid_ogrid(nghost+1:mxgrid_ogrid-nghost,nghost+1:mygrid_ogrid-nghost, &
!!         1,ivar1:ivar2) = fgrid_ogrid_tmp(1:nxgrid_ogrid,1:nygrid_ogrid,1,ivar1:ivar2)
!!       fgrid_ogrid(nghost+1:mxgrid_ogrid-nghost,nghost+1:mygrid_ogrid-nghost, &
!!         2,ivar1:ivar2) = fgrid_ogrid_tmp(1:nxgrid_ogrid,1:nygrid_ogrid,1,ivar1:ivar2)
!!       fgrid_ogrid(nghost+1:mxgrid_ogrid-nghost,nghost+1:mygrid_ogrid-nghost, &
!!         3,ivar1:ivar2) = fgrid_ogrid_tmp(1:nxgrid_ogrid,1:nygrid_ogrid,1,ivar1:ivar2)
!!       fgrid_ogrid(nghost+1:mxgrid_ogrid-nghost,nghost+1:mygrid_ogrid-nghost, &
!!         4,ivar1:ivar2) = fgrid_ogrid_tmp(1:nxgrid_ogrid,1:nygrid_ogrid,1,ivar1:ivar2)
!!       fgrid_ogrid(nghost+1:mxgrid_ogrid-nghost,nghost+1:mygrid_ogrid-nghost, &
!!         5,ivar1:ivar2) = fgrid_ogrid_tmp(1:nxgrid_ogrid,1:nygrid_ogrid,1,ivar1:ivar2)
!!       fgrid_ogrid(nghost+1:mxgrid_ogrid-nghost,nghost+1:mygrid_ogrid-nghost, &
!!         6,ivar1:ivar2) = fgrid_ogrid_tmp(1:nxgrid_ogrid,1:nygrid_ogrid,1,ivar1:ivar2)
!!       fgrid_ogrid(nghost+1:mxgrid_ogrid-nghost,nghost+1:mygrid_ogrid-nghost, &
!!         7,ivar1:ivar2) = fgrid_ogrid_tmp(1:nxgrid_ogrid,1:nygrid_ogrid,1,ivar1:ivar2)
!!       fgrid_cartesian(nghost+1:mxgrid-nghost,nghost+1:mygrid-nghost, &
!!         1,ivar1:ivar2) = fgrid_cartesian_tmp(1:nxgrid,1:nygrid,1,ivar1:ivar2)
!!       fgrid_cartesian(nghost+1:mxgrid-nghost,nghost+1:mygrid-nghost, &
!!         2,ivar1:ivar2) = fgrid_cartesian_tmp(1:nxgrid,1:nygrid,1,ivar1:ivar2)
!!       fgrid_cartesian(nghost+1:mxgrid-nghost,nghost+1:mygrid-nghost, &
!!         3,ivar1:ivar2) = fgrid_cartesian_tmp(1:nxgrid,1:nygrid,1,ivar1:ivar2)
!!       fgrid_cartesian(nghost+1:mxgrid-nghost,nghost+1:mygrid-nghost, &
!!         4,ivar1:ivar2) = fgrid_cartesian_tmp(1:nxgrid,1:nygrid,1,ivar1:ivar2)
!!       fgrid_cartesian(nghost+1:mxgrid-nghost,nghost+1:mygrid-nghost, &
!!         5,ivar1:ivar2) = fgrid_cartesian_tmp(1:nxgrid,1:nygrid,1,ivar1:ivar2)
!!       fgrid_cartesian(nghost+1:mxgrid-nghost,nghost+1:mygrid-nghost, &
!!         6,ivar1:ivar2) = fgrid_cartesian_tmp(1:nxgrid,1:nygrid,1,ivar1:ivar2)
!!       fgrid_cartesian(nghost+1:mxgrid-nghost,nghost+1:mygrid-nghost, &
!!         7,ivar1:ivar2) = fgrid_cartesian_tmp(1:nxgrid,1:nygrid,1,ivar1:ivar2)
!!     endif
!! 
!!     nfbuf_og= (/ mxgrid_ogrid, mygrid_ogrid, mzgrid_ogrid,ivar2-ivar1+1/)
!!     nfbuf_cg= (/ mxgrid, mygrid, mzgrid,ivar2-ivar1+1/)
!!     call mpibcast_real(fgrid_ogrid, nfbuf_og, root)
!!     call mpibcast_real(fgrid_cartesian, nfbuf_cg, root)
!!     endsubroutine send_rcv_all_data
!***********************************************************************
!    subroutine initialize_send_ip_points
!!
!! Build arrays of interpolation data on processors that contain data 
!! necessary for interpolation on other processors. 
!!
!! apr-17/Jorgen: Coded
!!
!      use Mpicomm, only: mpirecv_int, mpisend_nonblock_int, mpibarrier, mpiwait
!      use Solid_Cells_Mpicomm, only: finalize_isend_init_interpol
!      integer :: i,iip,npoint
!      integer, dimension(ncpus) :: from_proc_curv_to_cart=0
!      integer, dimension(ncpus) :: from_proc_cart_to_curv=0
!      integer, dimension(:,:,:), allocatable :: ind_from_proc_curv
!      integer, dimension(:,:,:), allocatable :: ind_from_proc_cart
!      integer, dimension(:,:), allocatable :: ip_id_curv_to_cart
!      integer, dimension(:,:), allocatable :: ip_id_cart_to_curv
!      integer, dimension(:), allocatable   :: send_to_curv_to_cart
!      integer, dimension(:), allocatable   :: send_to_cart_to_curv
!      integer, dimension(:,:), allocatable :: send_data_curv_to_cart
!      integer, dimension(:,:), allocatable :: send_data_cart_to_curv
!      integer, dimension(:), allocatable   :: send_id_curv_to_cart
!      integer, dimension(:), allocatable   :: send_id_cart_to_curv
!      integer :: max_from_proc, from_proc
!      integer, dimension(2) :: nelements
!      integer :: size_arr, npoints_requested
!      integer, dimension(:), allocatable   :: tmp_arr1D
!      integer, dimension(:,:), allocatable :: tmp_arr2D
!      integer :: nreq0D,nreq1D,nreq2D
!      integer, dimension(ncpus-1) :: ireq0D,ireq1D
!      integer, dimension(3*(ncpus-1)) :: ireq2D
!      integer :: iter1,iter2
!! TODO: COULD THIS BE MOVED INTO SOLID_CELLS_OGRID_MPICOMM?
!      if(n_ip_curv_to_cart>0) then
!        do i=1,n_ip_curv_to_cart
!          from_proc=curvilinear_to_cartesian(i)%from_proc
!          if(from_proc/=iproc) then
!! Must access from_proc+1 instead of from_proc, to avoid accessing element 0
!            from_proc_curv_to_cart(from_proc+1)=from_proc_curv_to_cart(from_proc+1)+1
!          endif
!        enddo
!      endif
!!
!      max_from_proc=maxval(from_proc_curv_to_cart)
!      if(max_from_proc>0) then
!        allocate(ind_from_proc_curv(ncpus,max_from_proc,3))
!        allocate(ip_id_curv_to_cart(ncpus,max_from_proc))
!        do iip=0,ncpus-1
!          if(from_proc_curv_to_cart(iip+1)>0) then
!            npoint=0
!            do i=1,n_ip_curv_to_cart
!              if(curvilinear_to_cartesian(i)%from_proc==iip) then
!                npoint=npoint+1
!! Must access iip+1 instead of iip, to avoid accessing element 0
!                ind_from_proc_curv(iip+1,npoint,:)=curvilinear_to_cartesian(i)%ind_global_neighbour
!                ip_id_curv_to_cart(iip+1,npoint)=i
!              endif
!            enddo
!          endif
!        enddo
!      endif
!!
!      if(n_ip_cart_to_curv>0) then
!        do i=1,n_ip_cart_to_curv
!          from_proc=cartesian_to_curvilinear(i)%from_proc
!          if(from_proc/=iproc) then
!! Must access from_proc+1 instead of from_proc, to avoid accessing element 0
!            from_proc_cart_to_curv(from_proc+1)=from_proc_cart_to_curv(from_proc+1)+1
!          endif
!        enddo
!      endif
!!
!      max_from_proc=maxval(from_proc_cart_to_curv)
!      if(max_from_proc>0) then
!        allocate(ind_from_proc_cart(ncpus,max_from_proc,3))
!        allocate(ip_id_cart_to_curv(ncpus,max_from_proc))
!        do iip=0,ncpus-1
!         if(from_proc_cart_to_curv(iip+1)>0) then
!            npoint=0
!            do i=1,n_ip_cart_to_curv
!              if(cartesian_to_curvilinear(i)%from_proc==iip) then
!                npoint=npoint+1
!! Must access iip+1 instead of iip, to avoid accessing element 0
!                ind_from_proc_cart(iip+1,npoint,:)=cartesian_to_curvilinear(i)%ind_global_neighbour
!                ip_id_cart_to_curv(iip+1,npoint)=i
!              endif
!            enddo
!          endif
!        enddo
!      endif
!! 
!!  Arrays containing information about which points should be sent by what processor to this
!!  processor has now been created. Now, there should be some communication to let all processors
!!  know which grid points they should SEND and who should RECIEVE them.
!!
!!  Note: Code is repeated twice in stead of being programmed as a function, since some compilers do
!!  not support allocatable arrays as in/out from subroutines/functions
!!  Use som variant of processor number as unique MPI tag (iip,iip+ncpus,etc.) in communication.
!!
!!  Curvilinear to Cartesian
!!
!      nreq0D=0
!      nreq1D=0
!      nreq2D=0
!      do iip=0,ncpus-1
!!  Send number of points requested from each processors, and send what points are requested
!!  if the number of points is larger than zero.
!!  Avoid sending to oneself
!        if(iip/=iproc) then
!          nreq0D=nreq0D+1
!          call mpisend_nonblock_int(from_proc_curv_to_cart(iip+1),iip,iip,ireq0D(nreq0D))
!          if(from_proc_curv_to_cart(iip+1)>0) then
!            nelements=(/ from_proc_curv_to_cart(iip+1),3 /)
!            do i=1,3
!              nreq2D=nreq2D+1
!              call mpisend_nonblock_int(ind_from_proc_curv(iip+1,1:nelements(1),i),nelements(1),iip,200+i,ireq2D(nreq2D))
!            enddo
!            nreq1D=nreq1D+1
!            call mpisend_nonblock_int(ip_id_curv_to_cart(iip+1,1:nelements(1)),nelements(1),iip,iip+2*ncpus,ireq1D(nreq1D))
!          endif
!        endif
!      enddo
!      allocate(send_to_curv_to_cart(0))
!      allocate(send_data_curv_to_cart(0,3))
!      allocate(send_id_curv_to_cart(0))
!      do iip=0,ncpus-1
!!  Recieve data from all processors. If any points are requested, create array of request.
!!  Avoid recieving from oneself
!        if(iip/=iproc) then
!          call mpirecv_int(npoints_requested,iip,iproc)
!!  Allocation/deallocation in a very inefficient manner, but this is only done during pre-processing
!!  so memory effieient code is a priority.
!          if(npoints_requested>0) then
!!  Expand array
!            size_arr=size(send_to_curv_to_cart)
!            allocate(tmp_arr1D(size_arr))
!            tmp_arr1D = send_to_curv_to_cart
!            deallocate(send_to_curv_to_cart)
!            allocate(send_to_curv_to_cart(size_arr+npoints_requested))
!            send_to_curv_to_cart(1:size_arr)=tmp_arr1D
!            deallocate(tmp_arr1D)
!            !
!            send_to_curv_to_cart(size_arr+1:size_arr+npoints_requested)=iip
!            nelements=(/ npoints_requested,3 /)
!!  Expand array
!            allocate(tmp_arr2D(size_arr,3))
!            tmp_arr2D = send_data_curv_to_cart
!            deallocate(send_data_curv_to_cart)
!            allocate(send_data_curv_to_cart(size_arr+npoints_requested,3))
!          
!            send_data_curv_to_cart(1:size_arr,:)=tmp_arr2D
!            deallocate(tmp_arr2D)
!            do i=1,3
!              call mpirecv_int(send_data_curv_to_cart(size_arr+1:size_arr+npoints_requested,i),nelements(1),iip,200+i)
!            enddo
!!  Expand array
!            allocate(tmp_arr1D(size_arr))
!            tmp_arr1D=send_id_curv_to_cart
!            deallocate(send_id_curv_to_cart)
!            allocate(send_id_curv_to_cart(size_arr+npoints_requested))
!            send_id_curv_to_cart(1:size_arr)=tmp_arr1D
!            deallocate(tmp_arr1D)
!            call mpirecv_int(send_id_curv_to_cart(size_arr+1:size_arr+npoints_requested),npoints_requested,iip,iproc+2*ncpus)
!          endif
!        endif
!      enddo
!      do i=1,nreq0D
!        call mpiwait(ireq0D(i))
!      enddo
!      do i=1,nreq1D
!        call mpiwait(ireq1D(i))
!      enddo
!      do i=1,nreq2D
!        call mpiwait(ireq2D(i))
!      enddo
!      call mpibarrier
!      !call finalize_isend_init_interpol(ireq1D,ireq2D,nreq1D,nreq2D)
!!
!!  Cartesian to curvilinear
!!
!      nreq1D=0
!      nreq2D=0
!      do iip=0,ncpus-1
!!  Send number of points requested from each processors, and send what points are requested
!!  if the number of points is larger than zero.
!!  Avoid sending to oneself
!        if(iip/=iproc) then
!          nreq1D=nreq1D+1
!          call mpisend_nonblock_int(from_proc_cart_to_curv(iip+1),iip,iip+3*ncpus,ireq1D(nreq1D))
!          if(from_proc_cart_to_curv(iip+1)>0) then
!            nelements=(/ from_proc_cart_to_curv(iip+1),3 /)
!            nreq2D=nreq2D+2
!            call mpisend_nonblock_int(ind_from_proc_cart(iip+1,1:nelements(1),:),nelements,iip,iip+4*ncpus,ireq2D(nreq2D-1))
!            call mpisend_nonblock_int(ip_id_cart_to_curv(iip+1,1:nelements(1)),nelements(1),iip,iip+5*ncpus,ireq2D(nreq2D))
!          endif
!        endif
!      enddo
!      allocate(send_to_cart_to_curv(0))
!      allocate(send_data_cart_to_curv(0,3))
!      allocate(send_id_cart_to_curv(0))
!      do iip=0,ncpus-1
!!  Recieve data from all processors. If any points are requested, create array of request.
!!  Avoid recieving from oneself
!        if(iip/=iproc) then
!          call mpirecv_int(npoints_requested,iip,iproc+3*ncpus)
!!  Allocation/deallocation in a very inefficient manner, but this is only done during pre-processing
!!  so memory effieient code is a priority.
!          if(npoints_requested>0) then
!!  Expand array
!            size_arr=size(send_to_cart_to_curv)
!            allocate(tmp_arr1D(size_arr))
!            tmp_arr1D = send_to_cart_to_curv
!            deallocate(send_to_cart_to_curv)
!            allocate(send_to_cart_to_curv(size_arr+npoints_requested))
!            send_to_cart_to_curv(1:size_arr)=tmp_arr1D
!            deallocate(tmp_arr1D)
!            !
!            send_to_cart_to_curv(size_arr+1:size_arr+npoints_requested)=iip
!            nelements=(/ npoints_requested,3 /)
!!  Expand array
!            allocate(tmp_arr2D(size_arr,3))
!            tmp_arr2D = send_data_cart_to_curv
!            deallocate(send_data_cart_to_curv)
!            allocate(send_data_cart_to_curv(size_arr+npoints_requested,3))
!            send_data_cart_to_curv(1:size_arr,:)=tmp_arr2D
!            deallocate(tmp_arr2D)
!            call mpirecv_int(send_data_cart_to_curv(size_arr+1:size_arr+npoints_requested,:),nelements,iip,iproc+4*ncpus)
!!  Expand array
!            allocate(tmp_arr1D(size_arr))
!            tmp_arr1D=send_id_cart_to_curv
!            deallocate(send_id_cart_to_curv)
!            allocate(send_id_cart_to_curv(size_arr+npoints_requested))
!            send_id_cart_to_curv(1:size_arr)=tmp_arr1D
!            deallocate(tmp_arr1D)
!            call mpirecv_int(send_id_cart_to_curv(size_arr+1:size_arr+npoints_requested),npoints_requested,iip,iproc+5*ncpus)
!          endif
!        endif
!      enddo
!      call finalize_isend_init_interpol(ireq1D,ireq2D,nreq1D,nreq2D)
!!
!!  Deallocate arrays not not needed later
!!
!      if(allocated(ind_from_proc_curv))  deallocate(ind_from_proc_curv)
!      if(allocated(ind_from_proc_cart))  deallocate(ind_from_proc_cart)
!      if(allocated(ip_id_curv_to_cart))  deallocate(ip_id_curv_to_cart)
!      if(allocated(ip_id_cart_to_curv))  deallocate(ip_id_cart_to_curv)
!!
!!  Translate recieved global indices to local indices and save the module variables for communication 
!!
!      size_arr=size(send_data_curv_to_cart(:,1))
!      allocate(send_curvilinear_to_cartesian(size_arr))
!      do i=1,size_arr
!        send_curvilinear_to_cartesian(i)%send_to_proc=send_to_curv_to_cart(i)
!        send_curvilinear_to_cartesian(i)%ip_id=send_id_curv_to_cart(i)
!        call ind_global_to_local_curv(send_data_curv_to_cart(i,:), &
!            send_curvilinear_to_cartesian(i)%i_near_neighbour,lcheck_init_interpolation)
!      enddo
!      size_arr=size(send_data_cart_to_curv(:,1))
!      allocate(send_cartesian_to_curvilinear(size_arr))
!      do i=1,size_arr
!        send_cartesian_to_curvilinear(i)%send_to_proc=send_to_cart_to_curv(i)
!        send_cartesian_to_curvilinear(i)%ip_id=send_id_cart_to_curv(i)
!        call ind_global_to_local_cart(send_data_cart_to_curv(i,:), &
!            send_cartesian_to_curvilinear(i)%i_near_neighbour,lcheck_init_interpolation)
!      enddo
!!
!!  Set some auxiliary parameters to help with the interpolation communication
!!  Global to module, not to proc
!!
!      size_arr=size(send_data_curv_to_cart(:,1))
!      n_procs_send_curv_to_cart=0
!      if(size_arr>0) then
!        n_procs_send_curv_to_cart=1
!        do i=2,size_arr
!          if(send_curvilinear_to_cartesian(i)%send_to_proc /= &
!              send_curvilinear_to_cartesian(i-1)%send_to_proc) then
!            n_procs_send_curv_to_cart=n_procs_send_curv_to_cart+1
!          endif
!        enddo
!      endif
!      allocate(n_ip_to_proc_curv_to_cart(n_procs_send_curv_to_cart))
!      n_ip_to_proc_curv_to_cart=1
!      do i=2,size_arr
!        if(send_curvilinear_to_cartesian(i)%send_to_proc == &
!            send_curvilinear_to_cartesian(i-1)%send_to_proc) then
!          n_ip_to_proc_curv_to_cart=n_ip_to_proc_curv_to_cart+1
!        endif
!      enddo
!      size_arr=size(send_data_cart_to_curv(:,1))
!      n_procs_send_cart_to_curv=0
!      if(size_arr>0) then
!        n_procs_send_cart_to_curv=1
!        do i=2,size_arr
!          if(send_cartesian_to_curvilinear(i)%send_to_proc /= &
!              send_cartesian_to_curvilinear(i-1)%send_to_proc) then
!            n_procs_send_cart_to_curv=n_procs_send_cart_to_curv+1
!          endif
!        enddo
!      endif
!      allocate(n_ip_to_proc_cart_to_curv(n_procs_send_cart_to_curv))
!      n_ip_to_proc_cart_to_curv=1
!      do i=2,size_arr
!        if(send_cartesian_to_curvilinear(i)%send_to_proc == &
!            send_cartesian_to_curvilinear(i-1)%send_to_proc) then
!          n_ip_to_proc_cart_to_curv=n_ip_to_proc_cart_to_curv+1
!        endif
!      enddo
!!
!      n_procs_recv_curv_to_cart=count(from_proc_curv_to_cart.gt.0)
!      n_procs_recv_cart_to_curv=count(from_proc_cart_to_curv.gt.0)
!      allocate(n_ip_recv_proc_curv_to_cart(n_procs_recv_curv_to_cart))
!      allocate(n_ip_recv_proc_cart_to_curv(n_procs_recv_cart_to_curv))
!      n_ip_recv_proc_curv_to_cart=pack(from_proc_curv_to_cart,from_proc_curv_to_cart.gt.0)
!      n_ip_recv_proc_cart_to_curv=pack(from_proc_cart_to_curv,from_proc_cart_to_curv.gt.0)
!      allocate(procs_recv_curv_to_cart(n_procs_recv_curv_to_cart))
!      allocate(procs_recv_cart_to_curv(n_procs_recv_cart_to_curv))
!      iter1=1
!      iter2=1
!      do iip=0,ncpus-1
!        if(from_proc_cart_to_curv(iip+1)>0) then 
!          procs_recv_cart_to_curv(iter1)=iip
!          iter1=iter1+1
!        endif
!        if(from_proc_curv_to_cart(iip+1)>0) then
!          procs_recv_curv_to_cart(iter2)=iip
!          iter2=iter2+1
!        endif
!      enddo
!      max_send_ip_curv_to_cart=maxval(n_ip_to_proc_curv_to_cart)
!      max_recv_ip_curv_to_cart=maxval(n_ip_recv_proc_curv_to_cart)
!      max_send_ip_cart_to_curv=maxval(n_ip_to_proc_cart_to_curv)
!      max_recv_ip_cart_to_curv=maxval(n_ip_recv_proc_cart_to_curv)
!!
!!  Deallocate arrays
!!
!      deallocate(send_to_curv_to_cart)
!      deallocate(send_to_cart_to_curv)
!      deallocate(send_data_curv_to_cart)
!      deallocate(send_data_cart_to_curv)
!      deallocate(send_id_curv_to_cart)
!      deallocate(send_id_cart_to_curv)
!!
!!  Make sure that all processors complete this initialization before continuing
!!
!      call mpibarrier
!    endsubroutine initialize_send_ip_points
!!***********************************************************************
!  subroutine comm_ip_curv_to_cart_alt(f_cartesian,ivar1,ivar2)
!!
!!  Send and recieve necessary information to perform interpolation from 
!!  the curvilinear to the cartesian grid.
!!  This version of the communication between grids is optimized for the case
!!  where all processors have a part of the ogrid, and not all processors have
!!  parts of the cartesian grid near the ogrid. Hence, we want as much as 
!!  possible to be done with the data before passing it from the curvilinear
!!  to the cartesian grid
!!  
!!  02-okt-17/Jorgen: Coded
!!  
!    use Mpicomm, only: mpisend_nonblock_int,mpisend_nonblock_real,mpirecv_int,mpirecv_real,mpiwait,mpibarrier
!    real, dimension(mx,my,mz,mfarray), intent(inout) :: f_cartesian
!    integer, intent(in) :: ivar1,ivar2
!    integer, dimension(n_procs_send_curv_to_cart) :: ireq1D, ireq5D
!    integer, dimension(5) :: nbuf_farr
!    integer, dimension(max_recv_ip_curv_to_cart) :: id_bufi
!    real, dimension(ivar2-ivar1+1,max_send_ip_curv_to_cart) :: f_bufo
!    real, dimension(ivar2-ivar1+1,max_recv_ip_curv_to_cart) :: f_bufi
!    real, dimension(inter_len,inter_len,inter_len,ivar2-ivar1+1) :: farr
!    integer :: i,j,k,id,ipp
!    integer :: ii1,ii2,jj1,jj2,kk1,kk2
!    integer :: iter, send_to, recv_from
!    integer, dimension(3) :: inear_loc
!    integer, dimension(3) :: inear_glob
!    integer :: ind_send_first, ind_send_last, ind_recv_first, ind_recv_last
!    integer, dimension(max_send_ip_curv_to_cart) :: ip_bufo
!    integer :: ind, ipoly
!    real, dimension(3) :: xyz_ip
!!
!    if(interpolation_method==1) then
!      ii1=0; ii2=1; jj1=0; jj2=1; kk1=0; kk2=1
!    elseif(interpolation_method==2 .or. interpolation_method==3) then
!      ii1=1; ii2=1; jj1=1; jj2=1; kk1=1; kk2=1
! !   elseif(interpolation_method==4) then
! !     ii1=2; ii2=2; jj1=2; jj2=2; kk1=2; kk2=2
!    elseif(interpolation_method==5) then
!      ipoly=floor((interpol_order_poly)*0.5)
!      ii1=ipoly; ii2=ipoly; jj1=ipoly; jj2=ipoly; kk1=ipoly; kk2=ipoly
!    elseif(mod(interpolation_method,2)==0) then
!      ii1=interpolation_method/2; ii2=ii1; jj1=ii1; jj2=ii1; kk1=ii1; kk2=ii1
!    endif
!    nbuf_farr(5)=ivar2-ivar1+1
!!
!!  Before sending data, interpolate. Only send interpolated f-data
!!
!    ind_send_first=1
!    do iter=1,n_procs_send_curv_to_cart
!      ind_send_last=n_ip_to_proc_curv_to_cart(iter)+ind_send_first-1
!      send_to=send_curvilinear_to_cartesian(ind_send_last)%send_to_proc
!      nbuf_farr(5)=ind_send_last-ind_send_first+1
!      do ipp=1,nbuf_farr(5)
!        ind=ind_send_first+ipp-1
!        i=send_curvilinear_to_cartesian(ind)%i_near_neighbour(1)
!        j=send_curvilinear_to_cartesian(ind)%i_near_neighbour(2)
!        k=send_curvilinear_to_cartesian(ind)%i_near_neighbour(3)
!        farr(:,:,:,:)=f_ogrid(i-ii1:i+ii2,j-jj1:j+jj2,k-kk1:k+kk2,ivar1:ivar2)
!!
!!  Need global coordinates for interpolation.
!!
!        inear_glob=send_curvilinear_to_cartesian(ind)%i_global_neighbour
!        xyz_ip=send_curvilinear_to_cartesian(ind)%xyz
!        call interp_point_curv_to_cart_alt(xyz_ip,inear_glob,ivar1,ivar2,farr,f_bufo(:,ipp))
!      enddo
!      ip_bufo(1:nbuf_farr(5)) = send_curvilinear_to_cartesian(ind_send_first:ind_send_last)%ip_id
!      call mpisend_nonblock_int(ip_bufo(1:nbuf_farr(5)),nbuf_farr(5),send_to,send_to,ireq1D(iter))
!      call mpisend_nonblock_real(f_bufo(:,1:nbuf_farr(5)),nbuf_farr,send_to,send_to+ncpus,ireq5D(iter))
!      ind_send_first=ind_send_last+1
!    enddo
!    ind_recv_first=1
!    do iter=1,n_procs_recv_curv_to_cart
!      ind_recv_last=n_ip_recv_proc_curv_to_cart(iter)
!      recv_from=procs_recv_curv_to_cart(iter)
!      nbuf_farr(5)=ind_recv_last-ind_recv_first+1
!      call mpirecv_int(id_bufi(1:nbuf_farr(5)),nbuf_farr(5),recv_from,iproc)
!      call mpirecv_real(f_bufi(:,1:nbuf_farr(5)),nbuf_farr,recv_from,iproc+ncpus)
!      do ipp=1,nbuf_farr(5)
!        call transform_curv_to_cart(f_bufi(:,ipp),f_cartesian,id_bufi(ipp),ivar1,ivar2)
!      enddo
!    enddo
!!
!!  Interpolate remaining points 
!!
!    do id=1,n_ip_curv_to_cart
!      if(curvilinear_to_cartesian(id)%from_proc==iproc) then
!        inear_loc=curvilinear_to_cartesian(id)%ind_local_neighbour
!        farr(:,:,:,:)=f_ogrid(inear_loc(1)-ii1:inear_loc(1)+ii2, &
!          inear_loc(2)-jj1:inear_loc(2)+jj2,inear_loc(3)-kk1:inear_loc(3)+kk2,ivar1:ivar2)
!        call interpolate_point_curv_to_cart(f_cartesian,id,ivar1,ivar2,farr)
!      endif
!    enddo
!!
!!  Finalize nonblocking sends
!!
!    do iter=1,n_procs_send_curv_to_cart
!      call mpiwait(ireq1D(iter))
!      call mpiwait(ireq5D(iter))
!    enddo
!    call mpibarrier
!!
!  endsubroutine comm_ip_curv_to_cart_alt
!!***********************************************************************
!  logical function interp_lagrange4(farr_in,ivar1,ivar2,xxp,inear_glob,fp,lcart_to_curv,lcurv_to_cart,lcheck)
!
!!  Interpolate the value of f to (xp, yp) CURVILINEAR coordinate
!!  using the fourth-order lagrangian interpolation.
! 
!!  TODO: Extend to 3D
!!  TODO: Adjust nearest point to be ACTUALLY NEAREST, not bottom left corner
!!        Needed due to asymetric stencil
!
!!  The coefficients are determined by the 5x5 grid points surrounding the
!!  interpolation point.
!!  Global coordinates are used for the interpolation, to allow interpolation of 
!!  values outside this processors domain.
!
!!  28-sep-17/Jorgen: Coded
!
!      integer :: ivar1, ivar2
!      real, dimension (3) :: xxp
!      real, dimension (5,5,5,ivar2-ivar1+1) :: farr_in
!      !TODO
!      real, dimension (-2:2,-2:2,ivar2-ivar1+1) :: farr
!      real, dimension (ivar2-ivar1+1) :: fp
!      integer, dimension (3) :: inear_glob
!      logical :: lcart_to_curv, lcurv_to_cart, lcheck 
!
!      intent(in)  :: farr_in, ivar1, ivar2, xxp, inear_glob, lcheck, lcart_to_curv, lcurv_to_cart
!      intent(out) :: fp
!
!      real, dimension(-2:2) :: xglob, yglob
!      real :: dx_2_1,dx_20,dx_21,dx_22,dx_1_2,dx_10,dx_11,dx_12,dx0_2,dx0_1, &
!              dx01,dx02,dx1_2,dx1_1,dx10,dx12,dx2_2,dx2_1,dx20,dx21
!      real :: deltax_2,deltax_1,deltax0,deltax1,deltax2
!      real :: dy_2_1,dy_20,dy_21,dy_22,dy_1_2,dy_10,dy_11,dy_12,dy0_2,dy0_1, &
!              dy01,dy02,dy1_2,dy1_1,dy10,dy12,dy2_2,dy2_1,dy20,dy21
!      real :: deltay_2,deltay_1,deltay0,deltay1,deltay2
!      real :: lag_2,lag_1,lag0,lag1,lag2
!      real, dimension(-2:2,ivar2-ivar1+1) :: g_2,g_1,g0,g1,g2
!      real, dimension(-2:2,ivar2-ivar1+1) :: gp
!      integer :: i,ix0,ix1,iy0,iy1
!
!      interp_lagrange4= .true.
!      farr(:,:,:) = farr_in(:,:,3,:)
!
!!  Get grid points
!!
!      if(lcart_to_curv) then
!        xglob(-2:2) = xglobal(inear_glob(1)-2:inear_glob(1)+2)
!        yglob(-2:2) = yglobal(inear_glob(2)-2:inear_glob(2)+2)
!      elseif(lcurv_to_cart) then
!        xglob(-2:2) = xglobal_ogrid(inear_glob(1)-2:inear_glob(1)+2)
!        yglob(-2:2) = yglobal_ogrid(inear_glob(2)-2:inear_glob(2)+2)
!      else
!        print*,'interp_lagrange4:Not interpolated to any specific grid!'
!        interp_lagrange4= .false.
!        return
!      endif
!!
!!  Compute distance from xxp to surrounding grid points
!!  Needed for checking that inear_glob is correct, and for lagrange polynomials
!!
!      deltax_2 = xxp(1)-xglob(-2)
!      deltax_1 = xxp(1)-xglob(-1)
!      deltax0  = xxp(1)-xglob( 0)
!      deltax1  = xxp(1)-xglob( 1)
!      deltax2  = xxp(1)-xglob( 2)
!!
!      deltay_2 = xxp(2)-yglob(-2)
!      deltay_1 = xxp(2)-yglob(-1)
!      deltay0  = xxp(2)-yglob( 0)
!      deltay1  = xxp(2)-yglob( 1)
!      deltay2  = xxp(2)-yglob( 2)
!!
!!  Check that inear_glob actually points to the grid point closest to xxp
!!
!      if(lcheck) then
!        if((abs(deltax0)<=min(abs(deltax_2),abs(deltax_1),abs(deltax1),abs(deltax2))) .or. &
!           (abs(deltay0)<=min(abs(deltay_2),abs(deltay_1),abs(deltay1),abs(deltay2)))) then 
!        ! Everything okay
!        else
!          print*, 'interp_lagrange4: Interpolation point does not lie closest to center grid point.' 
!          print*, 'ix0, iy0, iz0 = ', inear_glob(1:3)
!          print*, 'xp, xglob(-2:2) = ', xxp(1), xglob
!          print*, 'yp, yglob(-2:2) = ', xxp(2), yglob
!          interp_lagrange4= .false.
!          return
!        endif
!      endif
!!
!!  Interpolate in x-direction
!! 
!!  Compute distances
!!
!      dx_2_1=xglob(-2)-xglob(-1)
!      dx_20 =xglob(-2)-xglob( 0)
!      dx_21 =xglob(-2)-xglob( 1)
!      dx_22 =xglob(-2)-xglob( 2)
!
!      dx_1_2=-dx_2_1
!      dx_10 =xglob(-1)-xglob( 0)
!      dx_11 =xglob(-1)-xglob( 1)
!      dx_12 =xglob(-1)-xglob( 2)
!
!      dx0_2 =-dx_20
!      dx0_1 =-dx_10
!      dx01  =xglob( 0)-xglob( 1)
!      dx02  =xglob( 0)-xglob( 2)
!
!      dx1_2 =-dx_21
!      dx1_1 =-dx_11
!      dx10  =-dx01
!      dx12  =xglob( 1)-xglob( 2)
!
!      dx2_2 =-dx_22
!      dx2_1 =-dx_12
!      dx20  =-dx02
!      dx21  =-dx12
!!
!!  Compute products of x-x_k/(x_i-x_k) for (i!=k)
!!
!      lag_2 = (deltax_1/dx_2_1)*(deltax0 /dx_20)*(deltax1/dx_21)*(deltax2/dx_22)
!      lag_1 = (deltax_2/dx_1_2)*(deltax0 /dx_10)*(deltax1/dx_11)*(deltax2/dx_12)
!      lag0  = (deltax_2/dx0_2 )*(deltax_1/dx0_1)*(deltax1/dx01 )*(deltax2/dx02 )
!      lag1  = (deltax_2/dx1_2 )*(deltax_1/dx1_1)*(deltax0/dx10 )*(deltax2/dx12 )
!      lag2  = (deltax_2/dx2_2 )*(deltax_1/dx2_1)*(deltax0/dx20 )*(deltax1/dx21 )
!
!      g_2 = farr(-2,:,:)
!      g_1 = farr(-1,:,:)
!      g0  = farr( 0,:,:)
!      g1  = farr( 1,:,:)
!      g2  = farr( 2,:,:)
!!
!!  Interpolate points in x-direction
!!
!      do i=ivar1,ivar2
!        gp(:,i) = lag_2*g_2(:,i)+lag_1*g_1(:,i)+lag0*g0(:,i)+lag1*g1(:,i)+lag2*g2(:,i)
!      enddo
!      !gp = lag_2*g_2 + lag_1*g_1 + lag0*g0 + lag1*g1 * lag2*g2
!!
!!  Interpolate in y-direction
!! 
!!  Compute distances
!!
!      dy_2_1=yglob(-2)-yglob(-1)
!      dy_20 =yglob(-2)-yglob( 0)
!      dy_21 =yglob(-2)-yglob( 1)
!      dy_22 =yglob(-2)-yglob( 2)
!
!      dy_1_2=-dy_2_1
!      dy_10 =yglob(-1)-yglob( 0)
!      dy_11 =yglob(-1)-yglob( 1)
!      dy_12 =yglob(-1)-yglob( 2)
!
!      dy0_2 =-dy_20
!      dy0_1 =-dy_10
!      dy01  =yglob( 0)-yglob( 1)
!      dy02  =yglob( 0)-yglob( 2)
!
!      dy1_2 =-dy_21
!      dy1_1 =-dy_11
!      dy10  =-dy01
!      dy12  =yglob( 1)-yglob( 2)
!
!      dy2_2 =-dy_22
!      dy2_1 =-dy_12
!      dy20  =-dy02
!      dy21  =-dy12
!!
!!  Compute products of x-x_k/(x_i-x_k) for (i!=k)
!!
!      lag_2 = (deltay_1/dy_2_1)*(deltay0 /dy_20)*(deltay1/dy_21)*(deltay2/dy_22)
!      lag_1 = (deltay_2/dy_1_2)*(deltay0 /dy_10)*(deltay1/dy_11)*(deltay2/dy_12)
!      lag0  = (deltay_2/dy0_2 )*(deltay_1/dy0_1)*(deltay1/dy01 )*(deltay2/dy02 )
!      lag1  = (deltay_2/dy1_2 )*(deltay_1/dy1_1)*(deltay0/dy10 )*(deltay2/dy12 )
!      lag2  = (deltay_2/dy2_2 )*(deltay_1/dy2_1)*(deltay0/dy20 )*(deltay1/dy21 )
!!
!!  Interpolate points in y-direction
!!
!      do i=ivar1,ivar2
!        fp(i) = lag_2*gp(-2,i)+lag_1*gp(-1,i)+lag0*gp(0,i)+lag1*gp(1,i)+lag2*gp(2,i)
!      enddo
!      !fp = lag_2*gp(-2,:)+lag_1*gp(-1,:)+lag0*gp(0,:)+lag1*gp(1,:)+lag2*gp(2,:)
!
!      if (lcheck) then
!        do i=1,ivar2-ivar1+1
!          if ((fp(i)>maxval(farr(:,:,i)).and.i/=3) .or. (fp(i)<minval(farr(:,:,i)).and.i/=3)) then
!!
!!  Compensate for overshoots by linear interpolation
!!
!            ix0=0; ix1=1; iy0=0; iy1=1
!            if(xglob(0)>xxp(1)) then
!              ix0=2; ix1=3
!            else
!              ix0=3; ix1=4
!            endif
!            if(yglob(0)>xxp(2)) then
!              iy0=2; iy1=3
!            else
!              iy0=3; iy1=4
!            endif
!            if(lcart_to_curv) then
!              interp_lagrange4= linear_interpolate_cartesian(farr_in(ix0:ix1,iy0:iy1,3:4,i),i,i,xxp, &
!                                           (/inear_glob(1)+ix0-3,inear_glob(2)+iy0-3,inear_glob(3)/),&
!                                           gp(0,i),lcheck_interpolation)
!            else
!              interp_lagrange4= linear_interpolate_curvilinear(farr_in(ix0:ix1,iy0:iy1,3:4,i),i,i,xxp, &
!                                           (/inear_glob(1)+ix0-3,inear_glob(2)+iy0-3,inear_glob(3)/),&
!                                           gp(0,i),lcheck_interpolation)
!            endif
!!
!            fp(i)=gp(0,i)
!          endif
!          if (fp(i)/=fp(i)) then
!            print*, 'interp_interpolate: interpolated value is NaN'
!            print*, 'interp_interpolate: xxp=', xxp
!            print*, 'interp_interpolate: i, fp(i)=', i, fp(i)
!            print*, '------------------'
!            interp_lagrange4=.false.
!          endif
!        enddo
!      endif
!!
!  endfunction interp_lagrange4
!***********************************************************************
!  logical function HO_interp_curv_loc(farr,ivar1,ivar2,xxp,inear_glob,inear_loc,fp,lcheck,order)
!!
!!  Interpolate the value of f to arbitrary (xp, yp) CURVILINEAR coordinate
!!  using the high-order lagrangian interpolation.
!! 
!!  The coefficients are determined by the 2xN grid points surrounding the
!!  interpolation point.
!! 
!!  TEST VERSION: ONLY APPLICABLE FOR SERIAL RUNS
!!
!!  21-may-17/Jorgen: Adapted from linear_interpolate_curvilinear_HO
!!
!      integer :: ivar1, ivar2
!      integer, intent(in) :: order
!      real, dimension (3) :: xxp
!      real, dimension (ivar2-ivar1+1) :: fp
!      real, dimension (order,order,ivar2-ivar1+1) :: farr
!      integer, dimension (3) :: inear_glob, inear_loc
!      logical :: lcheck
!!
!      intent(in)  :: farr, ivar1, ivar2, xxp, inear_glob, inear_loc, lcheck
!      intent(out) :: fp
!!
!      integer :: i,ix0,iy0,iz0
!      real :: x0,y0,xN,yN,dx,dy,dxN,dyN,xp,yp
!
!      real, dimension(order,ivar2-ivar1+1) :: gp
!      integer :: j
!!
!!  Determine index value of lowest lying corner point of grid box surrounding
!!  the interpolation point.
!!
!      HO_interp_curv_loc=.true.
!!
!      ix0=inear_glob(1); iy0=inear_glob(2); iz0=inear_glob(3)
!!
!!  Check if the grid point interval is really correct.
!!
!      if ((xglobal_ogrid(ix0)<=xxp(1) .and. xglobal_ogrid(ix0+1)>=xxp(1) .or. nxgrid==1) .and. &
!          (yglobal_ogrid(iy0)<=xxp(2) .and. yglobal_ogrid(iy0+1)>=xxp(2) .or. nygrid==1) .and. &
!          (zglobal_ogrid(iz0)<=xxp(3) .and. zglobal_ogrid(iz0+1)>=xxp(3) .or. nzgrid==1)) then
!        ! Everything okay
!      else
!        print*, 'HO_interp_curv_loc: Global interpolation point does not ' // &
!            'lie within the calculated grid point interval.'
!        print*, 'iproc = ', iproc_world
!        print*, 'mxgrid_ogrid, xglobal_ogrid(1), xglobal_ogrid(mxgrid_ogrid) = ', & 
!            mxgrid_ogrid, xglobal_ogrid(1), xglobal_ogrid(mxgrid_ogrid)
!        print*, 'mygrid_ogrid, yglobal_ogrid(1), yglobal_ogrid(mygrid_ogrid) = ', &
!            mygrid_ogrid, yglobal_ogrid(1), yglobal_ogrid(mygrid_ogrid)
!        print*, 'mzgrid_ogrid, zglobal_ogrid(1), zglobal_ogrid(mzgrid_ogrid) = ', & 
!            mzgrid_ogrid, zglobal_ogrid(1), zglobal_ogrid(mzgrid_ogrid)
!        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
!        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
!        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
!        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
!        HO_interp_curv_loc=.false.
!        return
!      endif
!      if ((xglobal_ogrid(ix0)-x_ogrid(inear_loc(1))<10.e-10 ).or. &
!          (yglobal_ogrid(iy0)-y_ogrid(inear_loc(2))<10.e-10 )) then
!        ! Everything okay
!      else
!        print*, 'HO_interp_curv_loc: Global and local interpolation point values do not match' 
!        print*, 'GLOBAL:'
!        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
!        print*, 'xp, xp0, xp1 = ', xxp(1), xglobal_ogrid(ix0), xglobal_ogrid(ix0+1)
!        print*, 'yp, yp0, yp1 = ', xxp(2), yglobal_ogrid(iy0), yglobal_ogrid(iy0+1)
!        print*, 'zp, zp0, zp1 = ', xxp(3), zglobal_ogrid(iz0), zglobal_ogrid(iz0+1)
!        print*, 'LOCAL:'
!        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
!        print*, 'xp, xp0, xp1 = ', xxp(1), x_ogrid(inear_loc(1)), x_ogrid(inear_loc(1)+1)
!        print*, 'yp, yp0, yp1 = ', xxp(2), y_ogrid(inear_loc(2)), y_ogrid(inear_loc(2)+1)
!        print*, 'zp, zp0, zp1 = ', xxp(3), z_ogrid(inear_loc(3)), z_ogrid(inear_loc(3)+1)
!        print*, 'DIFF:'
!        print*, 'xglobal_ogrid(ix0)-x_ogrid(inear_loc(1)) =', xglobal_ogrid(ix0)-x_ogrid(inear_loc(1))
!        print*, 'yglobal_ogrid(iy0)-y_ogrid(inear_loc(2)) =', yglobal_ogrid(iy0)-y_ogrid(inear_loc(2))
!        HO_interp_curv_loc=.false.
!        return
!      endif
!!
!!  Mapping to quadratic area with corners at (0,0) -- (1,1)
!!
!    x0 = xglobal_ogrid(ix0-floor(0.5*order)+1)
!    y0 = yglobal_ogrid(iy0-floor(0.5*order)+1)
!    xN = xglobal_ogrid(ix0+floor(0.5*order))
!    yN = yglobal_ogrid(iy0+floor(0.5*order))
!!
!    dx = xglobal_ogrid(ix0+1)-xglobal_ogrid(ix0)
!    dy = yglobal_ogrid(iy0+1)-yglobal_ogrid(iy0)
!!
!!  Check that gridspacing is correct
!!
!    if((xN-x0 - (order-1)*dx)<10.e-10 .and. (yN-y0 - (order-1)*dy)<10.e-10) then
!      !Do nothing
!    else
!      print*, 'HO_interp_curv_loc: Grid spacing error'
!      print*, 'x0, x1, xN = ', x0,xglobal_ogrid(ix0-floor(0.5*order)+2),xN
!      print*, 'dx, N*dx,xN-x0 = ', dx,(order-1)*dx,xN-x0
!      print*, 'y0, y1, yN = ', y0,yglobal_ogrid(iy0-floor(0.5*order)+2),yN
!      print*, 'dy, N*dy, yN-y0 = ', dy,(order-1)*dy, yN-y0
!    endif
!!
!    dxN = dx*(order-1)
!    dyN = dy*(order-1)
!    xp = (xxp(1)-x0)/dxN
!    yp = (xxp(2)-y0)/dyN
!    x0 = 0.
!    y0 = 0.
!    xN = 1.
!    yN = 1.
!    dx = 1./(order-1)
!    dy = 1./(order-1)
!
!    do i=1,order
!      call lagrange1D(farr(:,i,:),ivar1,ivar2,yp,y0,dy,order,gp(i,:))
!    enddo
!    call lagrange1D(gp,ivar1,ivar2,xp,x0,dx,order,fp) 
!
!!
!!  Do a reality check on the interpolation scheme.
!!
!    if (lcheck) then
!      do i=1,2!ivar2-ivar1+1
!        if (((fp(i)>maxval(farr(:,:,i)).and.i/=3) .or. (fp(i)<minval(farr(:,:,i)).and.i/=3)).and.ix0>floor(nx_ogrid*0.5)) then
!          print*, 'linear_interpolate_curvilinear_HO: interpolated value is smaller or larger than'
!          print*, 'linear_interpolate_curvilinear_HO: a values at the corner points, even after linearization!'
!          print*, 'linear_interpolate_curvilinear: xxp=', xxp
!          print*, 'linear_interpolate_curvilinear: ix0, iy0, iz0=', ix0,iy0,iz0
!          print*, 'linear_interpolate_curvilinear: x0, y0, z0=', &
!              xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
!          print*, 'linear_interpolate_curvilinear: i, fp(i)=', i, fp(i)
!          print*, 'linear_interpolate_curvilinear: farr=', farr(:,:,i)
!          print*, '------------------'
!          print*, 'DETAILS:'
!          do j=1,order
!            print*, 'j,x(j)',ix0-floor(0.5*order)+1+j-1,xglobal_ogrid(ix0-floor(0.5*order)+1+j-1)
!            print*, 'j,y(j)',iy0-floor(0.5*order)+1+j-1,yglobal_ogrid(iy0-floor(0.5*order)+1+j-1)
!            print*, 'farr(j,:,1)',farr(j,:,1)
!            print*, 'farr(j,:,2)',farr(j,:,2)
!          enddo
!          print*, '------------------'
!          HO_interp_curv_loc=.false.
!        endif
!        if (fp(i)/=fp(i)) then
!          print*, 'linear_interpolate_curvilinear: interpolated value is NaN'
!          print*, 'linear_interpolate_curvilinear: xxp=', xxp
!          print*, 'linear_interpolate_curvilinear: x0, y0, z0=', &
!              xglobal_ogrid(ix0), yglobal_ogrid(iy0), zglobal_ogrid(iz0)
!          print*, 'linear_interpolate_curvilinear: i, fp(i)=', i, fp(i)
!          print*, 'linear_interpolate_curvilinear: farr=', farr(:,:,i)
!          print*, '------------------'
!          HO_interp_curv_loc=.false.
!        endif
!      enddo
!    endif
!!
!  endfunction HO_interp_curv_loc
!!
!!***********************************************************************
!  subroutine lagrange1D(farr,ivar1,ivar2,xp,x0,dx,order,fp)
!!
!      integer, intent(in) :: ivar1, ivar2
!      integer, intent(in) :: order
!      real, dimension(order,ivar2-ivar1+1), intent(in)  :: farr
!      real, intent(in) :: xp, x0, dx
!!
!      real, dimension (ivar2-ivar1+1), intent(out) :: fp
!!
!      real, dimension(order) :: l
!      real :: xi
!      integer :: i,j,ivar
!      l = 1.
!      do i=1,order
!        xi = x0+dx*(i-1)
!        do j=1,order
!          if(i/=j) then
!            l(i) = l(i)*(xp-(x0+dx*(j-1)))/(xi-(x0+dx*(j-1)))
!          endif
!        enddo
!      enddo
!      fp=0.
!      do ivar=ivar1,ivar2
!        do i=1,order
!          fp(ivar)=fp(ivar)+l(i)*farr(i,ivar)
!        enddo
!      enddo
!!
!  endsubroutine lagrange1D
!!***********************************************************************
!  subroutine flow_cartesian_to_curvilinear(f_cartesian,f_og)
!
!    use General, only: linear_interpolate
!!
!!  Interpolate all flow variables from cartesian to curvilinear grid
!!  Only need to do this for the radial direction
!!
!!  Find position in (x,y,z)-coordinates from (r,theta,z)-system
!!  Use this to interpolate (linearly) from nearest neighbours
!!  Only works for iux:iuz and scalar values (rho,T,etc.) at present.
!!
!!  NOTE: Does not work for parallell runs
!!        should only be used for testing purposes
!!
!!  16-feb-17/Jorgen: Coded
!!
!    real, dimension (mx,my,mz,mfarray) :: f_cartesian
!    real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f_og
!    real, dimension (3) :: xyz
!    integer :: ivar1,ivar2
!    integer, dimension (3) :: inear
!    real, dimension (mvar) :: gp
!    integer :: i,j,k
!!
!    ivar1=1
!    ivar2=mvar
!    do k=n1_ogrid,n2_ogrid
!      do j=m1_ogrid,m2_ogrid
!        do i=l2_ogrid+1,l2_ogrid+nghost
!          xyz=(/ x_ogrid(i)*cos(y_ogrid(j))+xorigo_ogrid(1), &
!                  x_ogrid(i)*sin(y_ogrid(j))+xorigo_ogrid(2), &
!                  z_ogrid(k) /)
!          call find_near_ind_local_cart(inear,xyz,lcheck_interpolation)
!
!          if ( .not. linear_interpolate(f_cartesian,ivar1,ivar2,xyz,gp,inear,lcheck_interpolation) ) then
!            call fatal_error('linear_interpolate','interpolation from cartesian to curvilinear')
!          endif
!          f_og(i,j,k,iux)=gp(iux)*cos(y_ogrid(j))+gp(iuy)*sin(y_ogrid(j))
!          f_og(i,j,k,iuy)=-gp(iux)*sin(y_ogrid(j))+gp(iuy)*cos(y_ogrid(j))
!          f_og(i,j,k,iuz)=gp(iuz)
!          f_og(i,j,k,irho)=gp(irho)
!        enddo
!      enddo
!    enddo
!
!  endsubroutine flow_cartesian_to_curvilinear
!!***********************************************************************
!  subroutine flow_curvilinear_to_cartesian(f_cartesian)
!!
!!  Interpolate all flow variables from curvilinear to cartesian grid
!!
!!  NOTE: Does not work for parallell runs,
!!        should only be used for testing purposes
!!
!    real, dimension (mx,my,mz,mfarray) :: f_cartesian
!    real, dimension (3) :: rthz
!    integer :: i,j,k
!    integer, dimension (3) :: inear
!    real, dimension (mvar) :: gp
!    integer, parameter :: ivar1=1,ivar2=mvar
!  
!    do k=n1,n2
!      do j=m1,m2
!        do i=l1,l2
!          call get_polar_coords(x(i),y(j),z(k),rthz)
!          if((rthz(1)<=r_int_outer) .and.(rthz(1)>r_int_inner)) then  
!            call find_near_ind_local_curv(inear,rthz,lcheck_interpolation)  
!            if ( .not. linear_interpolate_ogrid(ivar1,ivar2,rthz,gp,inear,lcheck_interpolation) ) then
!              call fatal_error('linear_interpolate_ogrid','interpolation from curvilinear to cartesian')
!            endif
!            f_cartesian(i,j,k,iux)=gp(iux)*cos(rthz(2))-gp(iuy)*sin(rthz(2))
!            f_cartesian(i,j,k,iuy)=gp(iux)*sin(rthz(2))+gp(iuy)*cos(rthz(2))
!            f_cartesian(i,j,k,iuz:ivar2)=gp(iuz:ivar2)
!          endif
!        enddo
!      enddo
!    enddo
!!
!  endsubroutine flow_curvilinear_to_cartesian
!***********************************************************************
endmodule solid_cells
