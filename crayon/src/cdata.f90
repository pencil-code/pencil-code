! $Id: cdata.f90 13579 2010-03-31 11:58:14Z AxelBrandenburg $
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
  real, dimension (mx) :: x,dx_1,dx_tilde,xprim,dVol1
  real, dimension (my) :: y,dy_1,dy_tilde,yprim,dVol2
  real, dimension (mz) :: z,dz_1,dz_tilde,zprim,dVol3
  real, dimension (nx) :: dxyz_2, dxyz_4, dxyz_6
  real :: dx,dy,dz,dxmin,dxmax
  integer :: dimensionality
!
!  Grid parameters.
!
  real, dimension(3,1) :: coeff_grid=1.0
  real, dimension(3,2) :: xyz_step,xi_step_frac,xi_step_width=1.5
  real :: zeta_grid0=0.
  real :: xbot_slice=0.,xtop_slice=1.
  real :: ybot_slice=0.,ytop_slice=1.
  real :: zbot_slice=0.,ztop_slice=1.
  logical, dimension(3) :: lperi,lshift_origin
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
  real :: r_ref=1.,rsmooth=0.,box_volume=1.0
!
!  Time integration parameters.
!
  integer :: nt=10000000,it=1,itorder=3,itsub,it_timing=0
  real :: tmax=1e33, tstart=0.
  real :: max_walltime=0.0  ! in seconds
  double precision :: t
  real :: dt=0.0
  real :: cdt=0.4,cdts=1.,cdtr=1.,cdtc=1.,cdtv=0.25
  real :: eps_rkf=1e-8, eps_stiff=1e-6
  real :: ddt=0.
  real :: dt1_last=0.
  real :: dtmin=1.0e-6,dtmax=1.0e37
  logical :: lini_t_eq_zero=.false.
  real, dimension (nx) :: advec_uu,advec_shear,advec_lnrho
  real, dimension (nx) :: advec_cs2,advec_va2
  real, dimension (nx) :: diffus_nu,diffus_diffrho,diffus_eta,diffus_chi
  real, dimension (nx) :: dt1_advec,dt1_diffus,dt1_max
  real, dimension (3) :: alpha_ts=0.0,beta_ts=0.0,dt_beta_ts=1.0
  logical :: lfirstpoint=.false.,llastpoint=.false.
  logical :: lmaxadvec_sum=.false.,old_cdtv=.false.
  character (len=20) :: timestep_scaling(mvar)='cons_err'
!
!  Input/output of data.
!
  character (len=120) :: datadir='data'
  character (len=120) :: directory='',datadir_snap='',directory_snap=''
  real :: dsnap=100.,d2davg=100.,dvid=0.,dspec=impossible
  real :: crash_file_dtmin_factor=-1.0
  integer :: isave=100,ialive=0,isaveglobal=0
  logical :: lread_aux=.false., lwrite_aux=.false., lwrite_dvar=.false.
  logical :: lread_oldsnap=.false., lread_oldsnap_nomag=.false.
  logical :: save_lastsnap=.true.
  logical :: noghost_for_isave=.false.
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
!
!  Random numbers.
!
  integer, dimension(mseed) :: seed=0
  integer :: nseed=0, seed0=1812
  real, dimension (2) :: fran1,fran2
!
!  Module flags.
!
  logical :: ldensity_nolog=.false.
  logical :: lmpicomm=.false.
  logical :: lpostproc=.false.
  logical :: lwrite_slices=.false., lwrite_2daverages=.false.
  logical :: lwrite_slice_xy2,lwrite_slice_xy,lwrite_slice_xz,lwrite_slice_yz
  logical :: lgravx=.false.,lgravy=.false.,lgravz=.false.
  logical :: lgrav=.false.,lgravx_gas=.true.,lgravy_gas=.true.,lgravz_gas=.true.
  logical :: lwrite_ic=.true.,lnowrite=.false.,lserial_io=.false.
  logical :: lroot=.true.,ldebug=.false.
  logical :: lshear=.false.,lalpm=.false.
  logical :: lglobal=.false., lglobal_nolog_density=.false.
  logical :: lvisc_LES=.false.
  logical :: leos=.false., leos_idealgas=.false.
  logical :: lstart=.false., lrun=.false., lreloading=.false.
!
!  Variable indices (default zero, set later by relevant physics modules).
!
  integer :: nvar,naux,naux_com
  integer :: ilnrho=0, irho=0
  integer :: ipp=0,irhs=0
  integer :: iuu=0,iux=0,iuy=0,iuz=0,iss=0
  integer :: iox=0,ioy=0,ioz=0
  integer :: igg=0,igx=0,igy=0,igz=0,ipotself=0
  integer :: iaa=0,iax=0,iay=0,iaz=0
  integer :: ispx=0,ispy=0,ispz=0
  integer :: iug=0
  integer :: iam=0,iamx=0,iamy=0,iamz=0
  integer :: ie=0,iff=0,ifx=0,ify=0,ifz=0,idd=0
  integer :: ivisc_heat=0,ibb=0,ibx=0,iby=0,ibz=0,ijj=0,ijx=0,ijy=0,ijz=0
  integer :: iuxb=0,ijxb=0,iugu=0,iugh=0
  integer :: ishock=0
  integer :: icc=0,ilncc=0,ialpm=0,ietat=0
  integer :: iaphi=0,ibphi=0
  integer :: iglobal_bx_ext=0, iglobal_by_ext=0, iglobal_bz_ext=0
  integer :: iglobal_jx_ext=0, iglobal_jy_ext=0, iglobal_jz_ext=0
  integer :: iglobal_ex_ext=0, iglobal_ey_ext=0, iglobal_ez_ext=0
!
!  Parameters related to message passing.
!
  integer :: ix=-1,iy=-1,iz=-1,iz2=-1
  integer :: ix_loc=-1,iy_loc=-1
  integer :: iz_loc=-1,iz2_loc=-1
  integer :: iproc,ipx,ipy,ipz,root=0
  logical :: lprocz_slowest=.true.
  integer :: xlneigh,ylneigh,zlneigh ! `lower' processor neighbours
  integer :: xuneigh,yuneigh,zuneigh ! `upper' processor neighbours
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
  integer :: nname=0,nnamev=0,nnamexy=0,nnamexz=0
  integer :: nnamez=0,nnamey=0,nnamex=0
  real :: tdiagnos,t1ddiagnos,t2davgfirst
  integer, parameter :: mname=100,mnamev=100
  integer, dimension (mname) :: itype_name=0
  real, dimension (mname) :: fname=0.0, fweight=0.0
  real, dimension (nz,nprocz) :: z_allprocs=0.0
  real, dimension (:,:,:), allocatable :: fnamex, fnamey, fnamez
  real, dimension(:,:,:), allocatable :: fnamexy, fnamexz
  character (len=30) :: cname(mname),cform(mname)
  character (len=30) :: cnamev(mname)
  character (len=30), allocatable :: cnamexy(:),cformxy(:)
  character (len=30), allocatable :: cnamexz(:),cformxz(:)
  character (len=30), allocatable :: cnamez(:),cformz(:)
  character (len=30), allocatable :: cnamey(:),cformy(:)
  character (len=30), allocatable :: cnamex(:),cformx(:)
  logical :: lout=.false.,headt=.false.,headtt=.true.,ldt=.true.
  logical :: lfirst=.false.,llast=.false.,ldt_paronly=.false.
  logical :: ldiagnos=.false.,lvideo=.false.,lwrite_prof=.true.
  logical :: l2davg=.false.,l2davgfirst=.false.
  logical :: l1davg=.false.,l1davgfirst=.false.
  logical :: lwrite_yaverages=.true.,lwrite_zaverages=.true.
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
  logical, parameter :: lshock_var=lshock
  logical, parameter :: lmagnetic_var=lmagnetic
!
!  Boundary conditions.
!
  real, dimension(mcom) :: fbcx1=0.,fbcy1=0.,fbcz1=0., fbcz1_1=0., fbcz1_2=0.
  real, dimension(mcom) :: fbcx2=0.,fbcy2=0.,fbcz2=0., fbcz2_1=0., fbcz2_2=0.
  real, dimension(mcom) :: fbcx1_2=0.,fbcx2_2=0.
  real :: Udrift_bc=0.
  character (len=2*bclen+1), dimension(mcom) :: bcx='p',bcy='p',bcz='p'
  character (len=bclen), dimension(mcom) :: bcx1='',bcx2='', &
                                            bcy1='',bcy2='', &
                                            bcz1='',bcz2=''
  character (len=10), dimension(mfarray) :: varname
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
!
!  Initial conditions.
!
  integer :: init_loops=1
!
!  Fix SGI reading problem.
!
  logical :: lsgifix=.false.
!
!  Scrap yard. Please categorise these variables if you know what they do.
!  Even better: move them to their relevant modules.
!
  real :: ttransient=0.
  real :: b_ell=1., rbound=1.
  real :: grads0=0.   ! (1/c_p)ds/dz
  logical :: lmonolithic_io=.false.
  logical :: lrescaling_magnetic=.false.
  logical :: test_nonblocking=.false.
!
!***********************************************************************
endmodule Cdata
