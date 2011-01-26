! $Id$
!
!  Global variables are defined in this module.
!
!  Only truly global variables are put in Cdata. Module-specific variables
!  should be defined in the relevant modules. Sharing of module-specific
!  variables is not recommended, but may be accomplished using the subroutines
!  in the SharedVariables module.
!
module Cdata
!
  use Cparam
!
  implicit none
!
  public
!
!  Tiny and huge numbers.
!
  real, parameter :: one_real=1.0
  real, parameter :: epsi=5*epsilon(one_real),tini=5*tiny(one_real)
  real, parameter :: huge1=0.2*huge(one_real)
!
!  Default CVS Id.
!
  character (len=120) :: cvsid='no cvsid is given in start.in or run.in!'
!
!  Cartesian coordinate system.
!
  real, dimension (nxgrid) :: xgrid
  real, dimension (nygrid) :: ygrid
  real, dimension (nzgrid) :: zgrid
  real, dimension (mx) :: x,dx_1,dx_tilde,xprim,dVol_x
  real, dimension (my) :: y,dy_1,dy_tilde,yprim,dVol_y
  real, dimension (mz) :: z,dz_1,dz_tilde,zprim,dVol_z
  real, dimension (nx) :: dxyz_2, dxyz_4, dxyz_6
  real :: dx,dy,dz,dxmin,dxmax
  real, dimension (nxgrid) :: kx_fft,kx_fft2
  real, dimension (nygrid) :: ky_fft,ky_fft2
  real, dimension (nzgrid) :: kz_fft,kz_fft2
  real :: kx_ny,ky_ny,kz_ny
  integer :: dimensionality
!
!  Alternative coordinate systems: spherical, cylindric.
!
  character (len=9) :: coord_system='cartesian'
  logical :: lcartesian_coords=.true.
  logical :: lspherical_coords=.false.,lcylindrical_coords=.false.
  logical :: lsphere_in_a_box=.false.,lcylinder_in_a_box=.false.
  logical :: luse_latitude=.false.
  logical :: lcylindrical_gravity=.false.
  logical :: luniform_z_mesh_aspect_ratio=.false.
  real :: drcyl,dsurfxy,dsurfyz,dsurfzx,dvol
  real, dimension (nx) :: r_mn,r1_mn,r2_mn,r2_weight,dvolume,dvolume_1
  real, dimension (my) :: sinth,sin1th,sin2th,costh,cotth,sinth_weight
  real, dimension (my) :: cos1th,tanth
  real, dimension (nygrid) :: sinth_weight_across_proc
  real, dimension (nx) :: rcyl_mn,rcyl_mn1,rcyl_mn2,rcyl_weight
  real, dimension (nx,3) :: dline_1
  real, dimension (nrcyl) :: rcyl  ! used for phi-averages
!
!  Grid parameters.
!
  real, dimension(3) :: coeff_grid=1.0
  real, dimension(3) :: dxi_fact=1.0
  real, dimension(3) :: trans_width=0.0
  real, dimension(3,2) :: xyz_step,xi_step_frac,xi_step_width=1.5
  real :: zeta_grid0=0.
  real :: xbot_slice=0.,xtop_slice=1.
  real :: ybot_slice=0.,ytop_slice=1.
  real :: zbot_slice=0.,ztop_slice=1.
  logical, dimension(3) :: lperi,lshift_origin
  logical, dimension(3) :: lequidist=(/.true.,.true.,.true. /)
  character (len=labellen), dimension(3) :: grid_func='linear'
  real, dimension(0:nprocx) :: procx_bounds
  real, dimension(0:nprocy) :: procy_bounds
  real, dimension(0:nprocz) :: procz_bounds
!
!  Derivative parameters
!
  character (len=labellen) :: der2_type='standard'
!
!  Box dimensions.
!
  real, dimension(3) :: Lxyz,xyz0,xyz1=impossible,xyz_star=(/0.0,0.0,0.0/)
  real, dimension(3) :: Lxyz_loc,xyz0_loc,xyz1_loc
  real :: x0,y0,z0,Lx,Ly,Lz
  real :: r_int=0.,r_ext=impossible   ! for spherical shell problems
  real :: r_int_border=impossible,r_ext_border=impossible
  real :: r_ref=1.,rsmooth=0.,box_volume=1.0
!
!  Time integration parameters.
!
  integer :: nt=10000000, it=1, itorder=3, itsub=0, it_timing=0
  real :: tmax=1e33, tstart=0.0
  real :: max_walltime=0.0  ! in seconds
  double precision :: t
  real :: dt=0.0
  real :: cdt=0.4, cdts=1.0, cdtr=1.0, cdtc=1.0
  real :: cdtv=0.25, cdtv2=0.03, cdtv3=0.01
  real :: cdtvDim, eps_rkf=1e-8, eps_stiff=1e-6
  real :: ddt=0.0
  real :: dt1_last=0.0
  real :: dtmin=1.0e-6, dtmax=1.0e37
  logical :: lini_t_eq_zero=.false.
  real, dimension (nx) :: advec_uu,advec_shear,advec_hall,advec_csn2
  real, dimension (nx) :: advec_cs2,advec_va2,advec_crad2,advec_uud,advec_uun
  real, dimension (nx) :: advec_hypermesh_rho,advec_hypermesh_uu
  real, dimension (nx) :: advec_hypermesh_aa,advec_hypermesh_ss
  real, dimension (nx) :: diffus_nu,diffus_nu2,diffus_nu3
  real, dimension (nx) :: diffus_diffrho,diffus_diffrho3
  real, dimension (nx) :: diffus_eta,diffus_eta2,diffus_eta3
  real, dimension (nx) :: diffus_chi,diffus_chi3
  real, dimension (nx) :: diffus_diffrhon,diffus_diffrhon3
  real, dimension (nx) :: diffus_diffnd,diffus_diffnd3
  real, dimension (nx) :: diffus_pscalar,diffus_pscalar3
  real, dimension (nx) :: diffus_chiral,diffus_cr,diffus_chem
  real, dimension (nx) :: diffus_nud,diffus_nud3
  real, dimension (nx) :: diffus_nun,diffus_nun3
  real, dimension (nx) :: dt1_advec,dt1_diffus,dt1_max
  real, dimension (nx) :: dt1_reac, reac_chem
  real, dimension (3) :: alpha_ts=0.0,beta_ts=0.0,dt_beta_ts=1.0
  logical :: lfirstpoint=.false.,llastpoint=.false.
  logical :: lmaxadvec_sum=.false.,old_cdtv=.false.
  character (len=20) :: timestep_scaling(mvar)='cons_err'
!
!  Use of LSODE to solve the chemistry in a separate step
!  By default, sequential splitting method (1st order)
!  if lsplit_second, Strang splitting procedure (2nd order)
!
  logical :: llsode=.false.
  logical :: lstep1=.true.
  logical :: lchemonly=.false.
  logical :: lsplit_second=.false.
!
!  Input/output of data.
!
  character (len=120) :: datadir='data'
  character (len=120) :: directory='',datadir_snap='',directory_snap=''
  real :: dsnap=100.,d2davg=100.,dvid=0.,dspec=impossible, dsound=0., tsound=0., soundeps=1.e-4
  real :: crash_file_dtmin_factor=-1.0
  integer :: isave=100,ialive=0,isaveglobal=0
  logical :: lread_aux=.false., lwrite_aux=.false., lwrite_dvar=.false.
  logical :: lread_oldsnap=.false., lread_oldsnap_nomag=.false.
  logical :: lread_oldsnap_nopscalar=.false.
  logical :: lread_oldsnap_notestfield=.false.
  logical :: lread_oldsnap_notestscalar=.false.
  logical :: save_lastsnap=.true.
  logical :: noghost_for_isave=.false.
  logical :: lfake_parallel_io=.false.
  logical :: ltec=.false.
  logical :: lformat=.false.
  logical :: lread_less=.false.
!
!  Units (need to be in double precision).
!
  character (len=3) :: unit_system='cgs'
  double precision :: unit_length=impossible,unit_velocity=impossible
  double precision :: unit_density=impossible,unit_temperature=impossible
  double precision :: unit_magnetic=impossible
!
!  Derived units
!
  double precision :: unit_mass,unit_energy,unit_time,unit_flux
  double precision :: k_B,m_u,m_p,m_e,m_H,m_He,eV, &
                      chiH,chiH_,sigmaH_,sigmaSB,kappa_es
  double precision :: c_light=impossible,G_Newton=impossible,hbar=impossible
  real :: mu0=1., mu01=0. !  magnetic permeability [should be in Magnetic]
!
!  Rotation and shear parameters.
!
  real :: Omega=0.0, theta=0.0, qshear=0.0, Sshear=impossible, deltay=0.0
!DM : Omega is now used in the viscosity routine too, for Lambda effect in rotating
! coordinate. This should be taken care of by 'shared variables' if in future
! Omega should be moved from cdata to hydro. 
!
!  Random numbers.
!
  integer, dimension(mseed) :: seed=0
  integer :: nseed=0, seed0=1812
  real, dimension (2) :: fran1,fran2
!
!  Module flags.
!
  logical :: ldensity_nolog=.false., lwrite_stratification=.false.
  logical :: lmpicomm=.false., lforcing_cont=.false.
  logical :: lpostproc=.false.
  logical :: lwrite_slices=.false., lwrite_2daverages=.false.
  logical :: lwrite_sound=.false.
  logical :: lwrite_slice_xy2,lwrite_slice_xy,lwrite_slice_xz,lwrite_slice_yz
  logical :: lwrite_slice_xy3=.false.,lwrite_slice_xy4=.false.
  logical :: lgravx=.false.,lgravy=.false.,lgravz=.false.
  logical :: lgravx_gas=.true.,lgravy_gas=.true.,lgravz_gas=.true.
  logical :: lgravx_dust=.true.,lgravy_dust=.true.,lgravz_dust=.true.
  logical :: lgravr=.false.,lgravr_gas=.false.,lgravr_dust=.false.
  logical :: lwrite_ic=.true.,lnowrite=.false.,lserial_io=.false.
  logical :: lroot=.true.,ldebug=.false.,lfft=.true.
  logical :: lfirst_proc_x=.true.,lfirst_proc_y=.true.,lfirst_proc_z=.true.
  logical :: lfirst_proc_xy=.true.,lfirst_proc_yz=.true.,lfirst_proc_xz=.true.
  logical :: lfirst_proc_xyz=.true.
  logical :: llast_proc_x=.true.,llast_proc_y=.true.,llast_proc_z=.true.
  logical :: llast_proc_xy=.true.,llast_proc_yz=.true.,llast_proc_xz=.true.
  logical :: llast_proc_xyz=.true.
  logical :: lnorth_pole=.false.,lsouth_pole=.false.
  logical :: lpscalar_nolog=.false.
  logical :: lalpm=.false.
  logical :: lradiation_ray=.false.,lradiation_fld=.false.
  logical :: ldustdensity_log=.false.
  logical :: lneutraldensity_nolog=.false.
  logical :: lglobal=.false., lglobal_nolog_density=.false.
  logical :: lvisc_hyper=.false.,lvisc_LES=.false.
  logical :: lvisc_smagorinsky=.false.
  logical :: leos_temperature_ionization=.false.
  logical :: ltemperature_nolog=.false.
  logical :: leos_idealgas=.false., leos_chemistry=.false.
  logical :: leos_ionization=.false.,leos_fixed_ionization=.false.
  logical :: ltestperturb=.false.
  logical :: lweno_transport=.false.
  logical :: lstart=.false., lrun=.false., lreloading=.false.
  logical :: lenergy=.false.
!
!  Variable indices (default zero, set later by relevant physics modules).
!
  integer :: nvar,naux,naux_com
  integer :: ilnrho=0, irho=0
  integer :: irho_b=0, iss_b=0 ! Anelastic auxiliary variables (base state)
  integer :: ipp,irhs=0,iTTold=0
  integer :: ipoly=0
  integer :: ip11=0,ip12=0,ip13=0
  integer :: ip21=0,ip22=0,ip23=0
  integer :: ip31=0,ip32=0,ip33=0
  integer :: ipoly_fr=0
  integer :: iuu=0,iux=0,iuy=0,iuz=0,iss=0
  integer :: iox=0,ioy=0,ioz=0
  integer :: ispecialvar=0
  integer :: iuut=0,iuxt=0,iuyt=0,iuzt=0,ioot=0,ioxt=0,ioyt=0,iozt=0
  integer :: ibbt=0,ibxt=0,ibyt=0,ibzt=0,ijjt=0,ijxt=0,ijyt=0,ijzt=0, &
             ijxb=0, ijxbx=0, ijxby=0, ijxbz=0
  integer :: igg=0,igx=0,igy=0,igz=0,ipotself=0
  integer :: iaa=0,iax=0,iay=0,iaz=0
  integer :: ispx=0,ispy=0,ispz=0
  integer :: ifcr=0,ifcrx=0,ifcry=0,ifcrz=0
  integer :: iaatest=0,iaxtestpq=0,iaztestpq=0,iaxtest=0,iaytest=0,iaztest=0
  integer :: iuutest=0,iuxtestpq=0,iuztestpq=0,ihhtestpq=0
  integer :: ntestscalar=0,ntestfield=0,ntestflow=0
  integer :: icctest=0,icctestpq=0,iug=0
  integer :: iam=0,iamx=0,iamy=0,iamz=0
  integer :: ie=0,iff=0,ifx=0,ify=0,ifz=0,idd=0
  integer :: ivisc_heat=0,ibb=0,ibx=0,iby=0,ibz=0,ijj=0,ijx=0,ijy=0,ijz=0
  integer :: iuxb=0,iugu=0,iugh=0
  integer :: ishock=0,ishock_perp=0
  integer :: iyH=0,ihypvis=0,ihypres=0
  integer :: iecr=0,ismagorinsky
  integer :: iQrad=0,iSrad=0,ikappa=0,ilnTT=0,iTT=0,ikapparho=0
  integer :: iQrad2=0,ikapparho2=0
  integer :: iFrad=0,iFradx=0,iFrady=0,iFradz=0
  integer :: iFrad2=0,iFradx2=0,iFrady2=0,iFradz2=0
  integer :: igpotselfx=0, igpotselfy=0, igpotselfz=0, irhop=0
  integer :: icc=0,ilncc=0,ialpm=0,ietat=0
  integer :: iaphi=0,ibphi=0,ieth=0
  integer, dimension(ndustspec) :: iuud=0,iudx=0,iudy=0,iudz=0
  integer, dimension(ndustspec) :: ilnnd=0, ind=0,imd=0,imi=0
  integer, dimension(nchemspec) :: ichemspec=0
  integer :: ilnrhon=0,irhon=0,iuun=0,iunx=0,iuny=0,iunz=0
  integer :: iglobal_bx_ext=0, iglobal_by_ext=0, iglobal_bz_ext=0
  integer :: iglobal_jx_ext=0, iglobal_jy_ext=0, iglobal_jz_ext=0
  integer :: iglobal_ex_ext=0, iglobal_ey_ext=0, iglobal_ez_ext=0
!
!  Parameters related to message passing.
!
  integer :: ix=-1,iy=-1,iz=-1,iz2=-1,iz3=-1,iz4=-1
  integer :: ix_loc=-1,iy_loc=-1
  integer :: iz_loc=-1,iz2_loc=-1, iz3_loc=-1, iz4_loc=-1
  integer :: iproc,ipx,ipy,ipz,root=0
  logical :: lprocz_slowest=.true.
  integer :: xlneigh,ylneigh,zlneigh ! `lower' processor neighbours
  integer :: xuneigh,yuneigh,zuneigh ! `upper' processor neighbours
  integer :: poleneigh              ! `pole' processor neighbours
  integer :: llcorn,lucorn,uucorn,ulcorn ! (the 4 corners in yz-plane)
!
!  Variables to count the occurance of derivative calls per timestep
!  for optimisation purposes.  To use uncomment the array and
!  set optimise_ders=.true.
!
!debug  integer, dimension(mfarray,8,3,3) :: der_call_count=0 !DERCOUNT
!debug  logical, parameter :: loptimise_ders=.true.             !DERCOUNT
!
!
!  Pencil-related stuff.
!
  integer :: imn,m,n
  integer, dimension (ny*nz) :: mm,nn
  logical, dimension (ny*nz) :: necessary=.false.
  integer, dimension (my,mz) :: imn_array
!
!  Parameters related to the pencil check.
!
  real :: penc0=2.345678 ! `impossible' value -- must not be
                         ! too large, so expressions like
                         ! exp(gamma_m1*p%lnrho) don't cause overflow
  logical :: lpencil_check=.false., lpencil_check_small=.true.
  logical :: lrandom_f_pencil_check=.true.
  logical :: lpencil_init=.false.
  logical :: lpencil_requested_swap=.true., lpencil_diagnos_swap=.false.
  logical :: lpencil_check_diagnos_opti=.false.
  logical :: lpencil_check_at_work=.false.
  integer :: ipencil_swap=0
!
!  Variables related to calculating diagnostic output.
!
  character (len=1) :: comment_char='#'
  integer :: it1=10,it1d=impossible_int
  integer :: nname=0,nnamev=0,nnamexy=0,nnamexz=0,nnamerz=0
  integer :: nnamez=0,nnamey=0,nnamex=0,nnamer=0
  integer :: nname_sound=0, ncoords_sound=0
  integer :: nr_directions=1
  real :: tdiagnos,t1ddiagnos,t2davgfirst
  integer, parameter :: mname=100
  integer, dimension (mname) :: itype_name=0
  real, dimension (mname) :: fweight=0.0
  real, dimension(:)      , allocatable :: fname
  real, dimension(:,:)    , allocatable :: fnamer,fname_sound
  real, dimension(:,:,:)  , allocatable :: fnamex, fnamey, fnamez,fnamexy, fnamexz
  real, dimension(:,:,:,:), allocatable :: fnamerz
  integer, dimension (:,:), allocatable :: sound_coords_list
  real, dimension (nz,nprocz) :: z_allprocs=0.0
  real, dimension (nrcyl,nx) :: phiavg_profile=0.0
  character (len=30), allocatable :: cform(:),cform_sound(:), &
                                     cformxy(:),cformxz(:),cformrz(:), &
                                     cformz(:),cformy(:),cformx(:),cformr(:)
  character (len=30), allocatable :: cname(:),cnamev(:),cname_sound(:), &
                                     cnamexy(:),cnamexz(:),cnamerz(:), &
                                     cnamez(:),cnamey(:),cnamex(:),cnamer(:)
  logical :: lout=.false.,headt=.false.,headtt=.true.,ldt=.true.
  logical :: lfirst=.false.,llast=.false.,ldt_paronly=.false.
  logical :: ldiagnos=.false.,lvideo=.false.,lwrite_prof=.true.,lout_sound=.false.
  logical :: l2davg=.false.,l2davgfirst=.false.
  logical :: l1davg=.false.,l1davgfirst=.false.,l1dphiavg=.false.
  logical :: lwrite_xyaverages=.false.,lwrite_xzaverages=.false.,lwrite_yzaverages=.false., &
             lwrite_phizaverages=.false.
  logical :: lwrite_yaverages=.false.,lwrite_zaverages=.false.,lwrite_phiaverages=.false.
  logical :: ldiagnos_need_zaverages=.false.
  logical :: ltime_integrals=.false.
  character (len=1) :: slice_position='p'
! averaging over smaller box
  logical :: lav_smallx=.false.,loutside_avg=.false.
  real :: xav_max=impossible
  integer :: ixav_max=0
  real :: nVol,nVol1  !  For calculating averages in non-cartesian coordinates
!
! Averages of half the computational box:
! fname_half has two indices, the first contains the quantity averaged
! over northern hemisphere and the second contains the
! quantity averaged over southern hemisphere.
! yequator = [xyz0(2)+0.5*Lxyz(2) assigned in start.f90 and run.f90
! zequator = [xyz0(3)+0.5*Lxyz(3) assigned in start.f90 and run.f90
!
  real :: yequator=0.,zequator=0.
  logical :: lequatory,lequatorz
  integer, parameter :: mname_half=20
  integer, dimension (mname_half) :: itype_name_half=0.
  real, dimension (mname_half,2) :: fname_half
  integer :: name_half_max=0
  character (len=30) :: cname_half(mname_half)
!  Coordinates of the point where some quantities can be printed.
  integer :: lpoint=(l1+l2)/2,mpoint=(m1+m2)/2,npoint=(n1+n2)/2
  integer :: lpoint2=(l1+l2)/4,mpoint2=(m1+m2)/4,npoint2=(n1+n2)/4
!
!  Diagnostic variables (needs to be consistent with reset list in register.90).
!
  integer :: idiag_it=0         ! DIAG_DOC: number of time step
                                ! DIAG_DOC:   \quad(since beginning of job only)
  integer :: idiag_t=0          ! DIAG_DOC: time $t$ \quad(since start.csh)
  integer :: idiag_dt=0         ! DIAG_DOC: time step $\delta t$
  integer :: idiag_walltime=0   ! DIAG_DOC: wall clock time since start of
                                ! DIAG_DOC:   run.x, in seconds
  integer :: idiag_timeperstep=0! DIAG_DOC:
  integer :: idiag_rcylmphi=0   ! PHIAVG_DOC: cylindrical radius
                                ! PHIAVG_DOC: $\varpi = \sqrt{x^2+y^2}$
                                ! PHIAVG_DOC: (useful for debugging
                                ! PHIAVG_DOC:  azimuthal averages)
  integer :: idiag_phimphi=0    ! PHIAVG_DOC: azimuthal angle
                                ! PHIAVG_DOC: $\varphi = \arctan\frac{y}{x}$
                                ! PHIAVG_DOC: (useful for debugging)
  integer :: idiag_zmphi=0      ! PHIAVG_DOC: $z$-coordinate
                                ! PHIAVG_DOC: (useful for debugging)
  integer :: idiag_rmphi=0      ! PHIAVG_DOC: spherical radius
                                ! PHIAVG_DOC: $r=\sqrt{\varpi^2+z^2}$
                                ! PHIAVG_DOC: (useful for debugging)
  integer :: idiag_dtv=0        ! DIAG_DOC:
  integer :: idiag_dtdiffus=0   ! DIAG_DOC:
  integer :: idiag_Rmesh=0      ! DIAG_DOC: $R_{\rm mesh}$
  integer :: idiag_Rmesh3=0     ! DIAG_DOC: $R_{\rm mesh}^{(3)}$
  integer :: idiag_maxadvec=0   ! DIAG_DOC: maxadvec
  integer :: idiag_nu_LES=0     ! DIAG_DOC:
!
!  Emergency brake:
!   When toggled the code will stop at the next convenient point
!   (at the next call to check_emergency_break)
!
  logical :: lemergency_brake=.false.
!
!  Experimential scheme for copying snapshots.
!
  logical :: lcopysnapshots_exp=.false.
!
!  Write snapshots with no ghost cells in the missing direction,
!  can save lots of hard disk space for 2-D runs.
!
  logical :: lwrite_2d=.false.
!
!  Bidiagonal second derivatives; default to true.
!
  logical :: lbidiagonal_derij=.true.
!
!  Constant 'parameters' cannot occur in namelists, so in order to get the
!  now constant module logicals into the lphysics name list...
!  We have some proxies that are used to initialise private local variables
!  called lhydro etc, in the lphysics namelist!
!
  logical, parameter :: lhydro_var=lhydro
  logical, parameter :: ldensity_var=ldensity
  logical, parameter :: lentropy_var=lentropy
  logical, parameter :: ltemperature_var=ltemperature
  logical, parameter :: lshock_var=lshock
  logical, parameter :: lmagnetic_var=lmagnetic
  logical, parameter :: lforcing_var=lforcing
  logical, parameter :: llorenz_gauge_var=llorenz_gauge
  logical, parameter :: ldustvelocity_var=ldustvelocity
  logical, parameter :: ldustdensity_var=ldustdensity
  logical, parameter :: ltestscalar_var=ltestscalar
  logical, parameter :: ltestfield_var=ltestfield
  logical, parameter :: ltestflow_var=ltestflow
  logical, parameter :: linterstellar_var=linterstellar
  logical, parameter :: lcosmicray_var=lcosmicray
  logical, parameter :: lcosmicrayflux_var=lcosmicrayflux
  logical, parameter :: lshear_var=lshear
  logical, parameter :: lpscalar_var=lpscalar
  logical, parameter :: lradiation_var=lradiation
  logical, parameter :: leos_var=leos
  logical, parameter :: lchiral_var=lchiral
  logical, parameter :: lneutralvelocity_var=lneutralvelocity
  logical, parameter :: lneutraldensity_var=lneutraldensity
!
!  Variables related to Fourier spectra and structure functions.
!
  logical :: vel_spec=.false.,mag_spec=.false.,uxj_spec=.false.,vec_spec=.false.
  logical :: j_spec=.false.,jb_spec=.false.,oo_spec=.false.
  logical :: vel_phispec=.false.,mag_phispec=.false.,uxj_phispec=.false.,vec_phispec=.false.
  logical :: uxy_spec=.false., bxy_spec=.false., jxbxy_spec=.false.
  character (LEN=12) :: xy_spec=''
  logical :: EP_spec=.false.
  logical :: ro_spec=.false.,TT_spec=.false.,ss_spec=.false.,cc_spec=.false.,cr_spec=.false.
  logical :: lr_spec=.false.,r2u_spec=.false.,r3u_spec=.false.
  logical :: ou_spec=.false.,ab_spec=.false.,ub_spec=.false.
  logical :: har_spec=.false.,hav_spec=.false.
  logical :: oned=.false.,twod=.false.
  logical :: ab_phispec=.false.,ou_phispec=.false.
  logical :: rhocc_pdf=.false.,cc_pdf=.false.,lncc_pdf=.false.
  logical :: gcc_pdf=.false.,lngcc_pdf=.false.
  logical :: test_nonblocking=.false.,onedall=.false.
  logical :: lsfu=.false.,lsfb=.false.,lsfz1=.false.,lsfz2=.false.
  logical :: lsfflux=.false.
  logical :: lpdfu=.false.,lpdfb=.false.,lpdfz1=.false.,lpdfz2=.false.
!
  ! Auxiliary parameters for boundary conditions:
  real, dimension(mcom) :: fbcx1=0., fbcy1=0., fbcz1=0., fbcz1_1=0., fbcz1_2=0.
  real, dimension(mcom) :: fbcx2=0., fbcy2=0., fbcz2=0., fbcz2_1=0., fbcz2_2=0.
  real, dimension(mcom) :: fbcx1_2=0., fbcx2_2=0.
  ! Auxiliary parameters for distinct use only with bottom or top boundary:
  real, dimension(mcom) :: fbcx_bot=0., fbcx_top=0.
  real, dimension(mcom) :: fbcy_bot=0., fbcy_top=0.
  real, dimension(mcom) :: fbcz_bot=0., fbcz_top=0.
!
  real :: Udrift_bc=0.
  character (len=2*bclen+1), dimension(mcom) :: bcx='p',bcy='p',bcz='p'
  character (len=bclen), dimension(mcom) :: bcx1='',bcx2='', &
                                            bcy1='',bcy2='', &
                                            bcz1='',bcz2=''
  character (len=10), dimension(mfarray) :: varname
  character (len=labellen) :: force_lower_bound='',force_upper_bound=''
!
!  Parameters for freezing boundary zones.
!
  real :: xfreeze_square=impossible,yfreeze_square=impossible
  real :: rfreeze_int=-impossible,rfreeze_ext=-impossible
  real :: wfreeze=0.,wfreeze_int=0.,wfreeze_ext=0.
  real :: wborder=0.,wborder_int=0.,wborder_ext=0.
  real :: tborder=0.
  real :: fshift_int=-1.,fshift_ext=1.
  real, dimension(2) :: border_frac_x=0.0,border_frac_y=0.0,border_frac_z=0.0
  logical :: lfrozen_bcs_x=.false.,lfrozen_bcs_y=.false.,lfrozen_bcs_z=.false.
  logical, dimension(mcom) :: lfrozen_bot_var_x=.false.,lfrozen_top_var_x=.false.
  logical, dimension(mcom) :: lfrozen_bot_var_y=.false.,lfrozen_top_var_y=.false.
  logical, dimension(mcom) :: lfrozen_bot_var_z=.false.,lfrozen_top_var_z=.false.
  logical, dimension(mcom) :: lfreeze_varsquare=.false.
  logical, dimension(mcom) :: lfreeze_varint=.false.,lfreeze_varext=.false.
!
!  A buffer in which to construct an error message.
!
  character (len=255) :: errormsg
!
!  Auxiliary variables.
!
  character (len=10), dimension(maux) :: aux_var
  integer :: aux_count=1
  integer :: mvar_io=0
  integer :: ireac=0
  integer, dimension(nchemspec) :: ireaci=0
!
!  Number of iterations for multigrid solver.
!
  integer :: niter_poisson=30
!
!  Filtering parameters.
!
  integer :: iwig=0,nfilter=0
  real :: awig=1.0
  logical :: lrmwig_rho=.false.,lrmwig_full=.false.,lrmwig_xyaverage=.false.
!
!  Initial conditions.
!
  integer :: init_loops=1
!
!  Implicit advance of the radiative diffusion in the temperature equation.
!
  logical :: lADI=.false.
  real, dimension (mx) :: hcondADI
!
!  Particle-mesh schemes, such as drag force and particle self-gravity,
!  sometimes add force to ghost zones. In that case we must fold the force
!  on the ghost zones onto the other side of the physical domain.
!
  logical :: lfold_df=.false.
!
!  Shift data cube by one grid point each time-step.
!
  logical :: lshift_datacube_x=.false.
!
!  Kinematic flows (computed by nohydro, so must be defined here).
!  [AJ: should probably not be defined here; AB: would hydro.h be better??]
! DM : I suggest we move these to hydro_kinematic as that is going to
! do the kinematic flow.  I have introduced a variable in hydro_kinematic
! called kinematic_flow. If hydro_kinematic is used kinflow is overwritten
! by this variable.
!
  character (len=40) :: kinflow=''
  logical :: lkinflow_as_aux
  real :: ampl_kinflow_x=0., ampl_kinflow_y=0., ampl_kinflow_z=0.
  real :: kx_kinflow=1., ky_kinflow=1., kz_kinflow=1.
  real :: dtphase_kinflow=0.
! (DM) All previous kinematic stuff can go to hydro_kinematic but I am not sure how
! to accomodate the following.
  real :: kinematic_phase=0.
!
!  Scrap yard. Please categorise these variables if you know what they do.
!  Even better: move them to their relevant modules.
!
  real :: ttransient=0.
  real :: b_ell=1., rbound=1.
  real :: grads0=0.   ! (1/c_p)ds/dz
  logical :: lisotropic_advection=.false.
  logical :: lfargo_advection=.false.
  logical :: pretend_lnTT=.false.
  logical :: lmonolithic_io=.false.
  logical :: lrescaling_magnetic=.false.
  logical :: lrescaling_testscalar=.false.
  logical :: lrescaling_testfield=.false.
  logical :: llocal_iso=.false.
!
!  Allow particles modules to do something special with f; like calculating
!  the full vorticity field (see equ.f90:pde)
!
  logical :: lparticles_prepencil_calc
!***********************************************************************
endmodule Cdata
