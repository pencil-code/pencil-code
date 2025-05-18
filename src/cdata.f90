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
!  Default CVS Id.
!
  character (len=linelen) :: cvsid='no cvsid is given in start.in or run.in!'
!
!  Polar grid
!
  integer :: ncoarse=0
  logical :: lcoarse=.false., lcoarse_mn=.false.
  integer, dimension(2) :: mexts=(/-1,-1/)
  integer, dimension(:), allocatable :: nphis
  real, dimension(:), allocatable :: nphis1, nphis2
  real(KIND=rkind8) :: t=0., toutoff=0.
  real :: tslice, eps_rkf=1e-5, eps_stiff=1e-6, eps_rkf0=0.
  real :: dsound=0., tsound=0., soundeps=1.e-4
!
!  Units (need to be in real(KIND=rkind8)).
!
  character (len=3) :: unit_system='cgs'
  logical :: lfix_unit_std=.false.
  real(KIND=rkind8) :: unit_length=impossible,unit_velocity=impossible
  real(KIND=rkind8) :: unit_density=impossible,unit_temperature=impossible
  real(KIND=rkind8) :: unit_magnetic=impossible,unit_entropy=impossible
  real(KIND=rkind8) :: sigma_Thomson=impossible
  real :: mu0=1., mu01=0. !  magnetic permeability [should be in Magnetic]
!
!  Derived units
!
  real(KIND=rkind8) :: unit_mass,unit_energy,unit_time,unit_flux,unit_pressure
  real(KIND=rkind8) :: k_B,m_u,m_p,m_e,m_H,m_He,eV, &
                      chiH,chiH_,sigmaH_,sigmaSB,kappa_es
  real(KIND=rkind8) :: c_light=impossible,G_Newton=impossible,hbar=impossible
!
!  Derived units
!
  real(KIND=rkind8) :: sigmaSB_set=1., c_light_set=1., cp_set=1.
  real(KIND=rkind8) :: k_B_set=1., m_u_set=1.
  integer, dimension(:,:), allocatable :: nexts
  integer, dimension(:,:,:), allocatable :: ninds
  logical :: lfirstpoint=.false.
!
!  Cartesian coordinate system.
!
  real, dimension (nx,3) :: dline_1
  real, dimension (nx) :: dxyz_2, dxyz_4, dxyz_6, dVol
  real, dimension (nx) :: dxmax_pencil,dxmin_pencil
!BEGIN C BINDING

  integer :: l2i=mx-2*nghost+1
  integer :: m2i=my-2*nghost+1
  integer :: n2i=mz-2*nghost+1
  integer :: l2=mx-nghost
  integer :: m2=my-nghost
  integer :: n2=mz-nghost

  real :: dVol_glob
  real, dimension (mx) :: x,dx_1,dx2,dx_tilde,xprim,dVol_x,dVol1_x
  real, dimension (mx) :: dAxy_x, dAxz_x
  real, dimension (my) :: y,dy_1,dy2,dy_tilde,yprim,dVol_y,dVol1_y
  real, dimension (my) :: dAxy_y, dAyz_y
  real, dimension (mz) :: z,dz_1,dz2,dz_tilde,zprim,dVol_z,dVol1_z
  real, dimension (mz) :: dAyz_z, dAxz_z
  real :: dx,dy,dz,dxmin,dxmax
  real, dimension (-nghost:nghost) :: dx2_bound=0., dy2_bound=0., dz2_bound=0.
  real, dimension (nxgrid) :: xgrid, dx1grid, dxtgrid
  real, dimension (nygrid) :: ygrid, dy1grid, dytgrid
  real, dimension (nzgrid) :: zgrid, dz1grid, dztgrid
  real, dimension (mxgrid) :: xglobal
  real, dimension (mygrid) :: yglobal
  real, dimension (mzgrid) :: zglobal
  real :: kx_nyq,ky_nyq,kz_nyq
!
!  Exact coefficients for non-equidistant grid at boundaries
!
  real, dimension(-nghost:nghost,2) :: coeffs_1_x
  real, dimension(-nghost:nghost,2) :: coeffs_1_y
  real, dimension(-nghost:nghost,2) :: coeffs_1_z
!
!  Alternative coordinate systems: spherical, cylindric.
!
  character (len=9) :: coord_system='cartesian'
  logical :: lcartesian_coords=.true.
  logical :: lspherical_coords=.false.,lcylindrical_coords=.false.
  logical :: lpipe_coords=.false.
  logical :: lsphere_in_a_box=.false.,lcylinder_in_a_box=.false.
  logical :: luse_latitude=.false., luse_oldgrid=.true., luse_xyz1=.false.
  logical :: lcylindrical_gravity=.false.
  logical :: luniform_z_mesh_aspect_ratio=.false.
  logical :: lconcurrent=.true.
!
!  Simultaneous foreign code.
!
  integer :: tag_foreign=0
  logical :: lforeign=.false.,lforeign_comm_nblckg=.false.
!
!  Yin-Yang grid.
!
  logical :: lyinyang=.false., lyang=.false., lcutoff_corners=.false.
  character(LEN=labellen) :: cyinyang_intpol_type='bilinear'
  integer :: iyinyang_intpol_type=BILIN
  integer :: nzgrid_eff=nzgrid
  real, dimension(4) :: yy_biquad_weights=impossible
  integer :: nycut=my, nzcut=mz
  real :: rel_dang=0.
!
!  Cubed sphere grid.
!
  logical :: lcubed_sphere=.false.
!
  real :: drcyl,dsurfxy,dsurfyz,dsurfzx
  real, dimension (nx) :: r_mn,r1_mn,r2_mn,r2_weight
  real, dimension (my) :: sinth,sin1th,sin2th,costh,cotth,sinth_weight
  real, dimension (mz) :: sinph,cosph
  real, dimension (my) :: cos1th,tanth
  real, dimension (nygrid) :: sinth_weight_across_proc
  real, dimension (nx) :: rcyl_mn=1.,rcyl_mn1=1.,rcyl_mn2=1.,rcyl_weight
  real, dimension (nx) :: glnCrossSec
  real, dimension (nrcyl) :: rcyl  ! used for phi-averages
  real, dimension (mx) :: x12    ! for slope-limted-diffusion
  real, dimension (my) :: y12    ! for slope-limted-diffusion
  real, dimension (mz) :: z12    ! for slope-limted-diffusion
!
!  Grid parameters.
!
  real, dimension(3) :: coeff_grid=1.0
  real, dimension(3) :: dxi_fact=1.0
  real, dimension(3) :: trans_width=1.0
  real, dimension(3,2) :: xyz_step=1.0,xi_step_frac=1.0,xi_step_width=1.5
  real :: zeta_grid0=0.0
  real :: xbot_slice=0.0,xtop_slice=1.0
  real :: ybot_slice=0.0,ytop_slice=1.0
  real :: zbot_slice=0.0,ztop_slice=1.0
  real :: r_rslice=0.
  integer :: nth_rslice=min(nxgrid,nygrid,nzgrid)/2, nph_rslice=min(nxgrid,nygrid,nzgrid)
  real :: glnCrossSec0=0.0, CrossSec_x1=-1., CrossSec_x2=1., CrossSec_w=.1
  logical, dimension(3) :: lperi=.true., &                                       ! all directions periodic
                           lshift_origin=.false., lshift_origin_lower=.false., & ! don't shift origin
                           lpole=.false., &                                      ! in spherical coords: pole excluded
                           lequidist=.true.                                      ! grid equidistant in every direction
  logical :: lignore_nonequi=.false., lcart_equi=.true.
  character (len=labellen), dimension(3) :: grid_func='linear'
  character (len=labellen) :: pipe_func='error_function'
  integer :: nghost_read_fewer=0
!
! Processor related
!
  real, dimension(0:nprocx) :: procx_bounds
  real, dimension(0:nprocy) :: procy_bounds
  real, dimension(0:nprocz) :: procz_bounds

  integer, dimension(3) :: dim_mask=(/1,2,3/)
!
!  Derivative parameters
!
  character (len=labellen) :: der2_type='standard'
  logical :: lall_onesided=.false.
  character (len=labellen), dimension(3) :: xyz_units='one'
  real, dimension(3) :: Lxyz=impossible,xyz0=-pi,xyz1=impossible,xyz_star=0.0
  real, dimension(3) :: Lxyz_loc,xyz0_loc,xyz1_loc
  real :: r_int=0.,r_ext=impossible   ! for spherical shell problems
  real :: r_int_border=impossible,r_ext_border=impossible
  real :: r_ref=1.,rsmooth=0.,box_volume=1.0
  real :: Area_xy=1., Area_yz=1., Area_xz=1.
  logical :: lfirst=.false.,llast=.false.,ldt_paronly=.false.
  logical :: ldt=.true.
  logical :: lcourant_dt=.true.
  logical :: lupdate_courant_dt=.false.
!
!  Time integration parameters.
!
  real :: dt=0.0
  real :: tmax=1e33, tstart=0.0
  real :: max_walltime=0.0  ! in seconds
  real :: dt_incr=0.0, dt0=0.
  real :: cdt=0.9, cdts=1.0, cdtr=1.0, cdtc=1.0, cdt_poly=1.0
 !real :: cdtv=0.15, cdtv2=0.03, cdtv3=0.01
!AB: 5 autotests failed after having decreased cdtv. I suggest to reassess
!AB: this more carefully and discuss it first in the newsletter.
  real :: cdtv=0.25, cdtv2=0.03, cdtv3=0.01
  real :: cdtsrc=0.2, cdtf=0.9
  real :: ddt=0.0, dtinc=0.5, dtdec=0.5
  real :: dtmin=1.0e-6, dtmax=1.0e37, dt_epsi=1e-7, dt_ratio=1e-5
  real :: nu_sts=0.1
  integer :: permute_sts=0
  integer:: ireset_tstart=2
  integer :: num_substeps = 3
!
!  Parameters related to message passing.
!
!
  integer, dimension(-1:1,-1:1,-1:1) :: neighbors = 0
  integer, dimension(26) :: iproc_comm = -1
  integer :: nproc_comm = 0
  integer :: ix=-1,iy=-1,iy2=-1,iz=-1,iz2=-1,iz3=-1,iz4=-1  !MR: dangerous names  ix -> ix_slice
  integer :: ix_loc=1,iy_loc=1, iy2_loc=1
  integer :: iz_loc=1,iz2_loc=1, iz3_loc=1, iz4_loc=1
  integer :: iproc=0,ipx=0,ipy=0,ipz=0,iproc_world=0,ipatch=0
  logical :: lprocz_slowest=.true.,lzorder=.false.,lmorton_curve=.false.,ltest_bcs=.true.,lcpu_timestep_on_gpu=.false., &
             lsuppress_parallel_reductions=.false.,lread_all_vars_from_device = .false.
  logical :: lac_sparse_autotuning=.false.
  integer :: xlneigh,ylneigh,zlneigh ! `lower' processor neighbours
  integer :: xuneigh,yuneigh,zuneigh ! `upper' processor neighbours
  integer :: poleneigh               ! `pole' processor neighbours
  integer :: nprocx_node=0, nprocy_node=0, nprocz_node=0
!
!  Data for registering of already updated variable ghost zones for only partly
!  updating by the *_after_timestep routines.
!  num_after_timestep: number of such routines; updated_var_ranges: list of already updated
!  variable ranges; ighosts_updated: counter for those, if -1 no registration is performed (default).
!
  integer, parameter :: num_after_timestep=5
  integer, dimension(2,2*num_after_timestep) :: updated_var_ranges=0
  integer :: ighosts_updated=-1
!
!  Box dimensions.
!
  real :: x0, y0, z0, Lx, Ly, Lz, wav1=impossible, wav1z=impossible
!
  logical :: lini_t_eq_zero=.false.
  logical :: lini_t_eq_zero_once=.false.
  real, dimension (nx) :: advec_cs2=0.
  real, dimension (nx) :: maxadvec=0., advec2=0., advec2_hypermesh=0.
  real, dimension (nx) :: maxdiffus=0., maxdiffus2=0., maxdiffus3=0., maxsrc=0.
  real, target, dimension (nx) :: dt1_max
  real, dimension (nx) :: reac_chem, reac_dust
  real    :: trelax_poly, reac_pchem
  real, dimension (5) :: alpha_ts=0.0,beta_ts=0.0,dt_beta_ts=1.0,bhat_ts=0.0,dt_bhat_ts=1.0
  logical :: lfractional_tstep_advance=.false.
  logical :: lfractional_tstep_negative=.true.
  logical :: lmaxadvec_sum=.false.,old_cdtv=.false.,leps_fixed=.true.
  logical :: lmaximal_cdtv=.false., lmaximal_cdt=.false.,lreiterate=.true.
  character (len=20), dimension(mvar) :: timestep_scaling='rel_err'
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
  character (len=fnlen) :: workdir='./', datadir='data', datadir_prestart='data_prestart'
  character (len=fnlen) :: directory='', datadir_snap='', directory_prestart=''
  character (len=fnlen) :: directory_snap='',directory_dist='',directory_collect=''
  character (len=fnlen) :: modify_filename='modify.dat'
  character (len=fmtlen) :: fmt_avgs='e14.5e3'
  logical :: lsnap=.false., lsnap_down=.false., lspec=.false., lspec_start=.false., lspec_at_tplusdt=.false.
  real :: dsnap=100., dsnap_down=0., d1davg=impossible, d2davg=100., dvid=0., dspec=impossible
  real :: dtracers=0., dfixed_points=0.
  real :: crash_file_dtmin_factor=-1.0
  real :: km0EM=0., km1EM=0.
  integer :: farray_smooth_width=6
  integer :: isave=100, ialive=0, isaveglobal=0, nv1_capitalvar=1
  logical :: lwrite_ts_hdf5=.true., lsave=.false.
  logical :: lread_aux=.false., lwrite_aux=.false., lwrite_dvar=.false.
  logical :: lenforce_maux_check=.true., lwrite_avg1d_binary = .false.
  logical :: lread_oldsnap=.false., lwrite_var_anyway=.false., lbackup_snap=.false.
  logical :: lwrite_last_powersnap=.false., lwrite_fsum=.false.
  logical :: lread_oldsnap_rho2lnrho=.false., lread_oldsnap_nomag=.false.
  logical :: lread_oldsnap_lnrho2rho=.false., lread_oldsnap_noshear=.false.
  logical :: lread_oldsnap_nohydro=.false., lread_oldsnap_nohydro_nomu5=.false.
  logical :: lread_oldsnap_onlyA=.false., lread_oldsnap_mskipvar=.false.
  logical :: lread_oldsnap_nohydro_efield=.false., lread_oldsnap_nohydro_ekfield=.false.
  logical :: ldivu_perp=.false.
  logical :: lread_oldsnap_nopscalar=.false.
  logical :: lread_oldsnap_notestfield=.false.
  logical :: lread_oldsnap_notestflow=.false.
  logical :: lread_oldsnap_notestscalar=.false.
  logical :: lread_oldsnap_noisothmhd=.false.
  logical :: lread_oldsnap_nosink=.false.
  logical :: lread_oldsnap_nocoolprof=.false.
  logical :: lnamelist_error=.false., ltolerate_namelist_errors=.false., lparam_nml=.false.
  logical :: lwrite_dim_again=.true., allproc_print=.true.
  logical :: lproc_print=.true.
  logical :: lseparate_persist=.false., ldistribute_persist=.false., lpersist=.true.
  logical :: lomit_add_data=.false.
  logical :: save_lastsnap=.true.
  logical :: noghost_for_isave=.false.
  logical :: ltec=.false.
  logical :: lformat=.false.
  logical :: lread_less=.false., lread_nogrid=.false.
  logical :: lread_global=.true.
  logical :: loutput_varn_at_exact_tsnap=.false.
  logical :: ldirect_access=.false.
  logical :: lread_from_other_prec=.false.       ! works so far only with io_dist!
  integer, dimension(3) :: downsampl=1, firstind=1, ndown=0, ngrid_down
  logical :: ldownsampl=.false., ldownsampling=.false., lrepair_snap=.false., linterpol_on_repair=.false.
  logical :: lastaroth_output=.false.
  character(LEN=fnlen) :: astaroth_dest=''
  integer, dimension(2) :: ivar_omit=(/0,0/)
  logical :: lzaver_on_input=.false.
  logical :: lfatal_num_vector_369=.true.
  logical :: lsmooth_farray=.false.
  logical :: lupdate_cvs=.false.
!
!  Entries related to the scale factor of the universe
!
  logical :: lread_scl_factor_file=.false., lread_scl_factor_file_new=.false.
  real :: scl_factor_target, Hp_target, appa_target, wweos_target
  real :: Hubble=0., ascale=1., sqrt_ascale=1.
  character(LEN=fnlen) :: ascale_type='default'
  integer :: enum_ascale_type = 0
!
! Debugging
!
  integer :: ip=14
!
!  Rotation and shear parameters.
!
  real :: Omega=0.0, theta=0.0, phi=0.0, qshear=0.0, Sshear=0.0, deltay=0.0
!DM : Omega is now used in the viscosity routine too, for Lambda effect in rotating
! coordinate. This should be taken care of by 'shared variables' if in future
! Omega should be moved from cdata to hydro.
!
!  Module flags.
!
  logical :: ldensity_nolog=.false., &
             lreference_state=.false., lfullvar_in_slices=.false., &
             lsubstract_reference_state=.false., ldensity_linearstart=.false.
  logical :: lforcing_cont=.false.


  logical :: lgravx=.false.,lgravy=.false.,lgravz=.false.
  logical :: lgravx_gas=.true.,lgravy_gas=.true.,lgravz_gas=.true.
  logical :: lgravx_dust=.true.,lgravy_dust=.true.,lgravz_dust=.true.
  logical :: lgravr=.false.
  logical :: lwrite_ic=.true.,lnowrite=.false.,lserial_io=.false.
  logical :: lmodify=.false.
  logical :: lroot=.true.,lcaproot=.false.,ldebug=.false.,lfft=.true.
  logical :: lproc_pt=.false., lproc_p2=.false.
  logical :: lfirst_proc_x=.true.,lfirst_proc_y=.true.,lfirst_proc_z=.true.
  logical :: lfirst_proc_xy=.true.,lfirst_proc_yz=.true.,lfirst_proc_xz=.true.
  logical :: lfirst_proc_xyz=.true.
  logical :: llast_proc_x=.true.,llast_proc_y=.true.,llast_proc_z=.true.
  logical :: llast_proc_xy=.true.,llast_proc_yz=.true.,llast_proc_xz=.true.
  logical :: llast_proc_xyz=.true.
  logical :: lnorth_pole=.false.,lsouth_pole=.false.
  logical :: lpscalar_nolog=.false.
  logical :: lalpm=.false., lalpm_alternate=.false.
  logical :: ldustdensity_log=.false.,lmdvar=.false.,ldcore=.false.
  logical :: lneutraldensity_nolog=.false.
  logical :: lvisc_smag=.false.
  logical :: lslope_limit_diff=.false.
  logical :: ltemperature_nolog=.false.
  logical :: ltestperturb=.false.
  logical :: lweno_transport=.false.
  logical :: lstart=.false., lrun=.false., lreloading=.false.
  logical :: ladv_der_as_aux=.false.
  logical :: lghostfold_usebspline = .false.
  logical :: lcooling_ss_mz = .false.
  logical :: lshock_heat = .true.
  real :: density_scale_factor=impossible
!
!  Used together with entropy, turns iss into ilntt (i.e., entropy
!  becomes log temperature). It does the same as using the
!  temperature_idealgas.f90 procedure, but draws on the more available
!  functionality extent in entropy.f90.
!
  logical :: pretend_lnTT=.false.
!END C BINDING
  logical :: lphase=.false.
!
!  Type counters.
!
  integer :: nvar,naux,naux_com,nscratch,nglobal,n_odevars=0
  real, dimension(:), allocatable :: f_ode, df_ode
  logical :: lode=.false.
!
!  Variable indices (default zero, set later by relevant physics modules).
!
  integer :: ilnrho=0, irho=0
  integer :: irho_b=0, iss_b=0 ! Anelastic auxiliary variables (base state)
  integer, dimension(ndustrad) :: iapn=0
  integer :: ipp,irhs=0,iTTold=0,irhsx=0,irhsy=0,irhsz=0
  integer :: ipoly=0
  integer :: ip11=0,ip12=0,ip13=0
  integer :: ip21=0,ip22=0,ip23=0
  integer :: ip31=0,ip32=0,ip33=0
  integer :: ipoly_fr=0
  integer :: iuu=0,iux=0,iuy=0,iuz=0,iss=0,iphiuu=0, ilorentz=0
  integer :: iuu0=0,iu0x=0,iu0y=0,iu0z=0
  integer :: ioo=0, iox=0, ioy=0, ioz=0
  integer :: ivv=0, ivx=0, ivy=0, ivz=0
  integer :: igradu11=0,igradu12=0,igradu13=0
  integer :: igradu21=0,igradu22=0,igradu23=0
  integer :: igradu31=0,igradu32=0,igradu33=0
  integer :: ispecialvar=0, ispecialvar2=0
  integer :: iuut=0,iuxt=0,iuyt=0,iuzt=0,ioot=0,ioxt=0,ioyt=0,iozt=0
  integer :: iuust=0,iuxst=0,iuyst=0,iuzst=0,ioost=0,ioxst=0,ioyst=0,iozst=0
  integer :: ibbt=0,ibxt=0,ibyt=0,ibzt=0,ijjt=0,ijxt=0,ijyt=0,ijzt=0, &
             ijxb=0, ijxbx=0, ijxby=0, ijxbz=0
  integer :: iuxb=0,iuxbx=0,iuxby=0,iuxbz=0
  integer :: iugb=0,iugbx=0,iugby=0,iugbz=0
  integer :: ibgu=0,ibgux=0,ibguy=0,ibguz=0
  integer :: ibdivu=0,ibdivux=0,ibdivuy=0,ibdivuz=0
  integer :: ibxf=0,ibyf=0,ibzf=0,ibbf=0
  integer :: ipotself=0
  integer :: iaa=0,iax=0,iay=0,iaz=0
  integer :: ispx=0,ispy=0,ispz=0
  integer :: ifcr=0,ifcrx=0,ifcry=0,ifcrz=0
  integer :: ihij=0,igij=0
  integer :: ihhT=0,ihhX=0,iggT=0,iggX=0,iStressT=0,iStressX=0,iStress_ij=0
  integer :: ihhTim=0,ihhXim=0,iggTim=0,iggXim=0,iStressTim=0,iStressXim=0
  integer :: iaatest=0,iaztestpq=0,iaxtest=0,iaytest=0,iaztest=0
  integer :: iuutest=0,iuztestpq=0,ihhtestpq=0
  integer :: iqx=0,iqy=0,iqz=0,iqq=0
  integer :: ntestscalar=0,ntestfield=0,ntestflow=0,ntestlnrho=0
  integer :: icctest=0,icctestpq=0,iug=0
  integer :: iam=0,iamx=0,iamy=0,iamz=0
  integer :: ivisc_heat=0,ibb=0,ibx=0,iby=0,ibz=0,ijj=0,ijx=0,ijy=0,ijz=0
  integer :: ieta_planet=0
  integer :: ibb_sph=0, ibb_sphr=0, ibb_spht=0, ibb_sphp=0
  integer :: inusmag=0, ietasmag=0
  integer :: iaak, iaakim, ieek, ieekim
  integer :: iee=0,iex=0,iey=0,iez=0,ialfven=0
  integer :: iFF_diff=0, iFF_diff1=0,  iFF_diff2=0, &
             iFF_div_uu=0, iFF_div_aa=0, iFF_div_ss=0, iFF_div_rho=0, iFF_char_c=0, iFF_heat=0
  integer :: isld_char=0, ivisc_forc=0,ivisc_forcx=0,ivisc_forcy=0,ivisc_forcz=0
  integer :: i_adv_der=0,i_adv_derx=0,i_adv_dery=0,i_adv_derz=0
  integer :: iuxbtest=0,ijxbtest=0,iugutest=0,iughtest=0,iSghtest=0
  integer :: ishock=0,ishock_perp=0
  integer :: iyH=0,ihypvis=0,ihypres=0
  integer :: iecr=0,ismagorinsky,iviscosity=0, inucl=0, inucrate=0
  integer :: isupsat=0
  integer :: iQrad=0,iSrad=0,ilnTT=0,iTT=0,ikapparho=0
  integer :: iKR_Frad=0,iKR_Fradx=0,iKR_Frady=0, iKR_Fradz=0
  integer :: igpotselfx=0, igpotselfy=0, igpotselfz=0
  integer :: icc=0,ilncc=0,ialpm=0,ietat=0
  integer :: iacc=0, issat=0, ittc=0, itauascalar=0
  integer :: iaphi=0,ibphi=0,ieth=0
  integer :: idet = 0
  integer :: iinvgrid=0
  integer :: iguij=0
  integer :: igu11=0,igu12=0,igu13=0
  integer :: igu21=0,igu22=0,igu23=0
  integer :: igu31=0,igu32=0,igu33=0
  integer :: icooling=0, inetheat=0
  integer, dimension(ndustspec) :: iuud=0,iudx=0,iudy=0,iudz=0
  integer, dimension(ndustspec) :: ilnnd=0, ind=0,imd=0,imi=0,idc=0,ilndc=0
  integer, dimension(ndustspec,ndustspec0) :: idcj=0,ilndcj=0
  integer, dimension(nchemspec) :: ichemspec=0
  integer :: ilnrhon=0,irhon=0, irhoe=0, iuun=0,iunx=0,iuny=0,iunz=0
  integer :: iglobal_bx_ext=0, iglobal_by_ext=0, iglobal_bz_ext=0
  integer :: iglobal_ax_ext=0, iglobal_ay_ext=0, iglobal_az_ext=0
  integer, dimension(3) :: iglobal_jext=0, iglobal_eext
  integer :: iglobal_lnrho0=0, iglobal_ss0=0
  integer :: icp=0, igpx=0, igpy=0, iRR=0, iss_run_aver=0
  integer :: iFenth=0, iss_flucz=0, iTT_flucz=0, irho_flucz=0
  integer :: iuu_fluc=0, iuu_flucx=0, iuu_flucy=0, iuu_flucz=0
  integer :: iuu_sph=0, iuu_sphr=0, iuu_spht=0, iuu_sphp=0
  integer :: ics=0, icool_prof=0
!
!  Variables to count the occurance of derivative calls per timestep
!  for optimisation purposes.  To use uncomment the array and
!  set optimise_ders=.true.
!
!debug  integer, dimension(mfarray,8,3,3) :: der_call_count=0 !DERCOUNT
!debug  logical, parameter :: loptimise_ders=.true.             !DERCOUNT
!
!  Pencil-related stuff.
!
  integer :: imn
  integer :: lglob=1
  integer, dimension (ny*nz) :: mm,nn
  logical, dimension (ny*nz) :: necessary=.false.
  integer :: necessary_imn=0
  integer, dimension (my,mz) :: imn_array
  character(LEN=labellen) :: shared_mem_name=''
!
!  Parameters related to the pencils
!
  logical, dimension(npencils) :: lpenc_diagnos   = .false.
  logical, dimension(npencils) :: lpenc_diagnos2d = .false.
  logical, dimension(npencils) :: lpenc_video     = .false.
  logical, dimension(npencils) :: lpenc_requested = .false.
  logical, dimension(npencils) :: lpencil         = .false.
!$omp threadprivate(lpencil)
!
!
!  Parameters related to the pencil check.
!
  real :: penc0=2.345678 ! `impossible' value -- must not be
                         ! too large, so expressions like
                         ! exp(gamma_m1*p%lnrho) don't cause overflow
  logical :: lpencil_check=.false., lpencil_check_small=.true.
  logical :: lpencil_check_no_zeros=.true.
  logical :: lpencil_init=.false.
  logical :: lpencil_requested_swap=.true., lpencil_diagnos_swap=.false.
  logical :: lpencil_check_diagnos_opti=.false.
  logical :: lpencil_check_at_work=.false.
  integer :: ipencil_swap=0
!
!  Variables related to calculating diagnostic output.
!
  character :: comment_char='#'
  integer :: it1=10,it1start=0,it1d=impossible_int,itspec=impossible_int
  integer :: itsnap=impossible_int
  integer :: nname=0,nnamev=0,nnamexy=0,nnamexz=0,nnamerz=0
  integer :: nnamez=0,nnamey=0,nnamex=0,nnamer=0
  integer :: nname_sound=0, ncoords_sound=0
  integer :: nr_directions=1
  integer :: itdiagnos
  real :: tspec,tdiagnos,dtdiagnos,t1ddiagnos,t2davgfirst,eps_rkf_diagnos
  real, dimension (mname) :: fweight=0.0
  integer, dimension(:)   , allocatable :: itype_name
  real, dimension(:)      , allocatable, target :: fname,fname_keep
  real, dimension(:,:)    , allocatable, target :: fnamer,fname_sound
  real, dimension(:,:,:)  , allocatable, target :: fnamex, fnamey, fnamez, fnamexy, fnamexz
  real, dimension(:,:,:,:), allocatable, target :: fnamerz
  integer, dimension(:,:) , allocatable :: sound_coords_list
  integer, dimension(:,:) , allocatable, target :: ncountsz
  character (len=fmtlen), allocatable :: cform(:),cformv(:),cform_sound(:), &
                                         cformxy(:),cformxz(:),cformrz(:), &
                                         cformz(:),cformy(:),cformx(:),cformr(:)
  character (len=fmtlen), allocatable :: cname(:),cnamev(:),cname_sound(:), &
                                         cnamexy(:),cnamexz(:),cnamerz(:), &
                                         cnamez(:),cnamey(:),cnamex(:),cnamer(:)
  integer, dimension(:), allocatable :: inds_max_diags, inds_sum_diags

!END C BINDING
  integer, target :: m,n
  integer :: nt=10000000, it=0, itorder=3, itsub=0, it_timing=0, it_rmv=0
  logical :: ltiming_io=.false.
  logical :: lwrite_slices=.false., lwrite_1daverages=.false., lwrite_2daverages=.false.
  logical :: lwrite_tracers=.false., lwrite_fixed_points=.false.
  logical :: lwrite_sound=.false.
  logical :: lwrite_slice_xy2=.false.,lwrite_slice_xy=.false.,lwrite_slice_xz=.false.,lwrite_slice_yz=.false.
  logical :: lwrite_slice_xy3=.false.,lwrite_slice_xy4=.false.,lwrite_slice_xz2=.false., lwrite_slice_r=.false.

  logical :: lout=.false.,headt=.false.,headtt=.true.,lrmv=.false.
  logical :: ldiagnos=.false.,lvideo=.false.,lwrite_prof=.true.,lout_sound=.false.
  logical :: lrhs_diagnostic_output=.false.
  logical :: ltracers=.false.,lfixed_points=.false.
  logical :: l2davg=.false.,l2davgfirst=.false.
  logical :: l1davg=.false.,l1davgfirst=.false.,l1dphiavg=.false.
  logical :: lwrite_xyaverages=.false.,lwrite_xzaverages=.false.
  logical :: lwrite_yzaverages=.false.,lwrite_phizaverages=.false.
  logical :: lwrite_yaverages=.false.,lwrite_zaverages=.false.
  logical :: lwrite_phiaverages=.false.
  logical :: ldiagnos_need_zaverages=.false.
  logical :: ltime_integrals=.false.
  logical :: lreset_seed=.false.
  logical :: lproper_averages=.false.
  character (len=1) :: slice_position='p'
! averaging over smaller box
  logical :: lav_smallx=.false.,loutside_avg=.false.
  real :: xav_max=impossible
  real :: nVol,nVol1  !  For calculating averages in non-Cartesian coordinates
!
!  Random numbers.
!
  integer, dimension(mseed) :: seed=0, seed2=0
  integer :: nseed=0, seed0=1812, ichannel1=1, ichannel2=1
  real, dimension (2) :: fran1,fran2
  logical :: lseed_global=.true., lseed_procdependent=.false.
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
  integer, dimension (mname_half) :: itype_name_half=0
  real, dimension (mname_half,2) :: fname_half
  integer :: name_half_max=0, ntmax1Dav=500, ntmax2Dav=500
  character (len=30) :: cname_half(mname_half)
!  Radius inside of which diagnostics are calculated for sphere_in_a_box models
  real :: radius_diag=1.0
!  Phase boundaries for ISM
  real :: ssmask1=0.0,ssmask2=0.0
!  Coordinates of the point where some quantities can be printed.
  integer :: lpoint=(mx+1)/2,mpoint=(my+1)/2,npoint=(mz+1)/2
  integer :: lpoint2=(mx+1)/4,mpoint2=(my+1)/4,npoint2=(mz+1)/4
  integer :: iproc_pt=0, iproc_p2=0
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
  integer :: idiag_dtv=0        ! DIAG_DOC: advective timestep as a
                                ! DIAG_DOC:   fraction of the actual one
  integer :: idiag_dtdiffus=0   ! DIAG_DOC: diffusive timestep as a
                                ! DIAG_DOC:   fraction of the actual one
  integer :: idiag_dtdiffus2=0  ! DIAG_DOC: hyperdiffusive (hyper2)  timestep
                                ! DIAG_DOC:   as a fraction of the actual one
  integer :: idiag_dtdiffus3=0  ! DIAG_DOC: hyperdiffusive (hyper3)  timestep
                                ! DIAG_DOC:   as a fraction of the actual one
  integer :: idiag_Rmesh=0      ! DIAG_DOC: $R_{\rm mesh}$
  integer :: idiag_Rmesh3=0     ! DIAG_DOC: $R_{\rm mesh}^{(3)}$
  integer :: idiag_maxadvec=0   ! DIAG_DOC: maxadvec
  integer :: idiag_eps_rkf=0    ! DIAG_DOC: time step accuracy threshold
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
!  Variables related to Fourier spectra and structure functions.
!
  logical :: vel_spec=.false.,mag_spec=.false.,uxj_spec=.false.,vec_spec=.false.
  logical :: j_spec=.false., jb_spec=.false., ja_spec=.false., oo_spec=.false., relvel_spec=.false.
  logical :: vel_phispec=.false.,mag_phispec=.false.,uxj_phispec=.false.,vec_phispec=.false.
  logical :: uxy_spec=.false., bxy_spec=.false., jxbxy_spec=.false.
  character (LEN=labellen*4) :: xy_spec=''
  character (LEN=labellen), dimension(n_xy_specs_max) :: xy_specs=''
  logical :: EP_spec=.false., hEP_spec=.false., nd_spec=.false., ud_spec=.false., abs_u_spec=.false.
  logical :: ro_spec=.false., TT_spec=.false., ss_spec=.false., cc_spec=.false., cr_spec=.false.
  logical :: sp_spec=.false., ssp_spec=.false., sssp_spec=.false., mu_spec=.false.
  logical :: lr_spec=.false., r2u_spec=.false., r3u_spec=.false., oun_spec=.false.
  logical :: np_spec=.false., np_ap_spec=.false., rhop_spec=.false., ele_spec=.false., pot_spec=.false.
  logical :: ux_spec=.false., uy_spec=.false., uz_spec=.false., a0_spec=.false., ucp_spec=.false.
  logical :: ou_spec=.false., ab_spec=.false., azbz_spec=.false., uzs_spec=.false.
  logical :: ub_spec=.false., Lor_spec=.false., EMF_spec=.false., Tra_spec=.false.
  logical :: GWs_spec=.false., GWh_spec=.false., GWm_spec=.false., Str_spec=.false., Stg_spec=.false.
  logical :: Gab_spec=.false., Gan_spec=.false., GBb_spec=.false.
  logical :: GWs_spec_boost=.false., GWh_spec_boost=.false.
  logical :: StT_spec=.false., StX_spec=.false.
  logical :: GWd_spec=.false., GWe_spec=.false., GWf_spec=.false., GWg_spec=.false.
  logical :: SCL_spec=.false., VCT_spec=.false., Tpq_spec=.false., TGW_spec=.false.
  logical :: SCL_spec_boost=.false., VCT_spec_boost=.false.
  logical :: har_spec=.false., hav_spec=.false., bb2_spec=.false., jj2_spec=.false., b2_spec=.false.
  logical :: oned=.false.,twod=.false.
  logical :: ab_phispec=.false.,ou_phispec=.false.
  logical :: rhocc_pdf=.false.,cc_pdf=.false.,lncc_pdf=.false.
  logical :: gcc_pdf=.false., lngcc_pdf=.false., lnspecial_pdf=.false., special_pdf=.false.
  logical :: cosEB_pdf=.false.
  logical :: ang_jb_pdf1d=.false., ang_ub_pdf1d=.false., ang_ou_pdf1d=.false.
  logical :: test_nonblocking=.false.,onedall=.false.
  logical :: lsfu=.false.,lsfb=.false.,lsfz1=.false.,lsfz2=.false.
  logical :: lsfflux=.false.
  logical :: lpdfu=.false.,lpdfb=.false.,lpdfz1=.false.,lpdfz2=.false.
  logical :: ou_omega=.false., cor_uu=.false., ab_kzspec=.false.,ou_kzspec=.false.
  logical :: ou_polar=.false., ab_polar=.false., jb_polar=.false.
  logical :: uut_spec=.false., uut_polar=.false., ouout_spec=.false.,ouout2_spec=.false.,ouout_polar=.false.
  logical :: out_spec=.false., uot_spec=.false.
  logical :: saffman_ub=.false., saffman_mag=.false., saffman_mag_c=.false.
  logical :: saffman_aa=.false., saffman_aa_c=.false., saffman_bb=.false.
  logical :: uu_fft3d=.false., oo_fft3d=.false., bb_fft3d=.false., jj_fft3d=.false.
  logical :: uu_xkyz=.false., oo_xkyz=.false., bb_xkyz=.false., jj_xkyz=.false.
  logical :: uu_kx0z=.false., oo_kx0z=.false., bb_kx0z=.false., jj_kx0z=.false.
  logical :: bb_k00z=.false., ee_k00z=.false., gwT_fft3d=.false.
  logical :: Em_specflux=.false., Hm_specflux=.false., Hc_specflux=.false.
!
  ! Auxiliary parameters for boundary conditions:
  real, dimension(mcom,2) :: fbcx=0., fbcx_2=0.
  real, dimension(mcom,2) :: fbcy=0., fbcy_1=0., fbcy_2=0.
  real, dimension(mcom,2) :: fbcz=0., fbcz_1=0., fbcz_2=0.
  ! Auxiliary parameters for distinct use only with bottom or top boundary:
  real, dimension(mcom) :: fbcx_bot=0., fbcx_top=0.
  real, dimension(mcom) :: fbcy_bot=0., fbcy_top=0.
  real, dimension(mcom) :: fbcz_bot=0., fbcz_top=0.
  ! Switch, if you wanna reset the boundary conditions
  logical :: lreset_boundary_values=.false.
!
  real :: Udrift_bc=0.
  character (len=2*bclen+1), dimension(mcom) :: bcx='p',bcy='p',bcz='p'
  character (len=bclen), dimension(mcom,2) :: bcx12='', bcy12='', bcz12=''
  character (len=labellen), dimension(mfarray) :: varname
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
  real :: theta_lower_border=impossible,wborder_theta_lower=0.
  real :: theta_upper_border=impossible,wborder_theta_upper=0.
  real :: fraction_tborder=1.
  logical :: lmeridional_border_drive=.false.
  real, dimension(2) :: border_frac_x=0.0,border_frac_y=0.0,border_frac_z=0.0, &
                        border_frac_r=0.0
  logical :: lborder_hyper_diff=.true.
  logical :: lfrozen_bcs_x=.false.,lfrozen_bcs_y=.false.,lfrozen_bcs_z=.false.
  logical, dimension(mcom) :: lfrozen_bot_var_x=.false.,lfrozen_top_var_x=.false.
  logical, dimension(mcom) :: lfrozen_bot_var_y=.false.,lfrozen_top_var_y=.false.
  logical, dimension(mcom) :: lfrozen_bot_var_z=.false.,lfrozen_top_var_z=.false.
  logical, dimension(mcom) :: lfreeze_varsquare=.false.
  logical, dimension(mcom) :: lfreeze_varint=.false.,lfreeze_varext=.false.
!
! Parameters for reading data for BCs.
!
  character(LEN=fnlen) :: bc_slc_dir=''
!
!  A buffer in which to construct an error message.
!
  character (len=linelen) :: errormsg
!
!  For mailing from code.
!
  character (len=linelen) :: mailcmd='',mailaddress='',submithost=''
!
  logical :: lstop_on_ioerror=.true.
!
!  Auxiliary variables.
!
  character (len=labellen), dimension(maux) :: aux_var
  integer :: aux_count=1
  integer :: mvar_io=0, mvar_down=-1, maux_down=-1, mskipvar=0
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
!  Particle-mesh schemes, such as drag force and particle self-gravity,
!  sometimes add force to ghost zones. In that case we must fold the force
!  on the ghost zones onto the other side of the physical domain.
!
  logical :: lfold_df=.false.
!
!  Reactive particles need their source terms for particle-fluid interaction
!  distributed over a large number of points to avoid shocks.
!  This means that nodes removed up to 3 nodes from the nearest grid point
!  are affected. This necessitates folding all the ghost zones into the main
!  domain.
!
  logical :: lfold_df_3points=.false.
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
  logical :: lkinflow_as_aux=.false.
  real :: ampl_kinflow_x=0., ampl_kinflow_y=0., ampl_kinflow_z=0.
  real :: kx_kinflow=1., ky_kinflow=1., kz_kinflow=1.
  real :: dtphase_kinflow=0.
!
!  Switch for Galilean-invariant advection for global disks (fargo). Used
!  in connection with special/fargo
!
  logical :: lfargo_advection=.false.
!
!  Switch for running global disks in a corotational frame,
!  and associated corotational radius and omega.
!
  logical :: lcorotational_frame=.false.
  real :: rcorot=1.0, Omega_corot=0.0
!
!  Switch for local isothermal approximation: a hardcoded time-invariant
!  temperature gradient for global disks.
!
  logical :: llocal_iso=.false.
!
!  Switch for using the full speed in less-than-3D simulations
!  when the velocity in the non-existent direction dominates.
!
  logical :: lisotropic_advection=.false.
!
!  Use lreport_undefined_diagnostics=F to suppress reporting of
!  undefined diagnostics. It puts hashes, but sometimes incorrectly so.
!
  logical :: lreport_undefined_diagnostics=.true.
!
!  Scrap yard. Please categorise these variables if you know what they do.
!  Even better: move them to their relevant modules.
!
  real :: ttransient=0.
  real :: b_ell=1., rbound=1.
  real :: grads0=0.   ! (1/c_p)ds/dz
  logical :: lmonolithic_io=.false.
  logical :: lrescaling_magnetic=.false.
  logical :: lrescaling_testscalar=.false.
  logical :: lrescaling_testfield=.false.
!
!  Dynamical diffusion coefficients with fixed mesh Reynolds number.
!
  real :: re_mesh=0.5
  logical :: ldynamical_diffusion=.false.
  logical :: ldyndiff_useumax = .true.
!
!  Background stratification.
!
  logical :: lstratz = .false.
  logical :: lnoghost_strati = .false.
!
!  Inverse timescale for running time average
!
  real :: tau_aver1 = 1.0
!
!  Info whether maux is needed and used on the GPU
!  Index for var is non-zero iff var is used on the GPU
!  The index corresponds to the vertex buffer index on Astaroth
!  Size of mfarray to make sure we can store the handle (for 1 to mvar zero)
!
   integer, dimension(mfarray) :: maux_vtxbuf_index = 0
   integer :: enum_unit_system = 0
!
!  Define and initialize lambda5, so that it can be used to tell whether
!  or not the chiral MHD special module is used.
!
  real :: lambda5 = 0.0
!
!  Variables for concurrency
!
  logical :: lmultithread=.false.
  logical :: l1dphiavg_save, l1davgfirst_save, ldiagnos_save, l2davgfirst_save
  logical :: lout_save, l1davg_save, l2davg_save, lout_sound_save, lvideo_save
  logical :: lchemistry_diag_save

  real(KIND=rkind8) :: t_save,tspec_save
  real :: t1ddiagnos_save,t2davgfirst_save,tslice_save,tsound_save

!$ logical, volatile, dimension(n_helperflags) :: lhelperflags=(/.false.,.false.,.false.,.false./)
!$ logical, volatile, dimension(n_helperflags) :: lmasterflags=(/.false.,.false.,.false.,.false./)
!$ logical, volatile :: lfarray_copied = .false.
  logical, dimension(npencils) :: lpencil_save = .false.
  integer :: num_helper_threads=1, thread_id=1
  integer, dimension(max_threads_possible) :: core_ids
!$ logical, volatile :: lhelper_run=.true., lhelper_perf
!$ logical :: loffload=.false.
!
! threadprivate definitions for OpenMP
!
!$omp threadprivate(inds_max_diags, inds_sum_diags)
!$omp threadprivate(dxyz_2,dxyz_4,dxyz_6,dvol,dxmax_pencil,dxmin_pencil,dline_1,lcoarse_mn, seed, m, n)
!$omp threadprivate(lfirstpoint,thread_id)
!$omp threadprivate(fname,fnamex,fnamey,fnamez,fnamer,fnamexy,fnamexz,fnamerz,fname_keep,fname_sound,ncountsz)
!$omp threadprivate(l1dphiavg, l1davgfirst, l2davgfirst, ldiagnos,lout, l1davg, l2davg, lout_sound, lvideo)
!$omp threadprivate(tspec,tdiagnos,t1ddiagnos,t2davgfirst,tslice,tsound,itdiagnos,dtdiagnos,eps_rkf_diagnos)
!
! For use in offloaded code:
!!$omp declare target(ldensity_nolog,l2,m2,n2)
!
!***********************************************************************
endmodule Cdata
