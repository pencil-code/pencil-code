! $Id: cdata.f90,v 1.215 2004-04-04 10:36:31 theine Exp $

module Cdata

!!! Global variables

  Use Cparam

  integer :: itorder=3
  real, dimension (mx) :: x
  real, dimension (my) :: y
  real, dimension (mz) :: z
  real, dimension (nrcyl) :: rcyl  ! used for phi-averages
  real, dimension (nx) :: x_mn,y_mn,z_mn,r_mn,rcyl_mn,phi_mn
  real, dimension (nx,3) :: evr    ! spherical unit radius vector
  real, dimension (nx) :: maxadvec2,maxdiffus,maxdss,maxdlnrho

  real, dimension (nx,3,3) :: sij,sdij  ! rate-of-strain tensor

  real, parameter :: pi=3.14159265358979324D0
  real, parameter :: epsi=5*epsilon(1.0),tini=5*tiny(1.0)
  real, dimension(3) :: Lxyz,xyz0,xyz1=impossible
  real :: t,dt=0.
  real :: cdt=0.4,cdtv=0.4,cdts=1.0,cdtr=1.0
  real :: cdtvDim
  real :: dx,dy,dz,dxmin,dxmax,drcyl,dsurfxy,dsurfyz,dsurfzx,dvol
  real :: dsnap=100.,d2davg=100.,dvid=100.,dtmin=1.e-6,dspec=impossible
  real :: r_int=0.,r_ext=impossible   ! for spherical shell problems
  real :: ttransient=0.
  real, dimension (2) :: fran1,fran2

  real, dimension(3) :: border_frac=0

  !  units (need to be in double precision)
  character (len=3) :: unit_system='cgs'
  double precision :: unit_length=1.,unit_velocity=1.,unit_density=1.,unit_temperature=1.
  ! Derived units
  double precision :: unit_mass,unit_energy,unit_time,unit_flux
  
  double precision :: k_B,m_p,m_e,m_H,m_He,eV,hbar, &
                      chiH,chiH_,sigmaH_,sigmaSB,kappa_es

  ! magnetic permeability
  real :: mu0=1., mu01=0.

!ajwm nu moved to viscosity module
!ajwm replaced nu, causes error in forcing to resolve
  real :: nu=0.,cmu,cnu2
  real :: tdiagnos,t2davgfirst
!! not used?  real :: rmean,rrms,rmax,u2m,um2,u2max,divurms,divumax,divu2max
  real :: o2m,om2,oum,epsK_hyper
  real :: UUmax,x0,y0,z0,Lx,Ly,Lz
  real :: grads0=0.   ! (1/c_p)ds/dz
  real :: Omega=0.,qshear=0.,Sshear=impossible
  real :: deltay=0. !(for shear; also used in forcing and output)

  integer, dimension(mseed) :: seed=0
  integer :: nseed=0
  integer :: nvar,naux,iuu=0,iux=0,iuy=0,iuz=0,ilnrho=0,iss=0
  integer :: igg=0,igx=0,igy=0,igz=0
  integer :: iaa=0,iax=0,iay=0,iaz=0
  integer :: ie=0,iff=0,ifx=0,ify=0,ifz=0,idd=0,ishock=0,iyH=0,ihyper=0
  integer :: iecr=0
  integer :: iQrad=0,iSrad=0,ikappa=0,ilnTT=0
  integer :: nt=1000000,it1=10
  integer :: it=1,itsub,ix=-1,iy=-1,iz=-1,iz2=-1
  integer :: ilncc=0
  integer :: iproc,ipx,ipy,ipz,root=0
  integer :: mvar_io=0,dimensionality
  integer :: iinit
  integer, parameter :: ninit=4
  integer, dimension(ndustspec) :: iuud,iudx,iudy,iudz,ind,irhod
  logical, dimension(3) :: lperi,lshift_origin
  character (len=labellen) ::fft_switch='fftpack'
!
!  coordinates of the point where some quantities can be printed
!  for now, these points only apply to the root processor.
!
  integer :: lpoint=(l1+l2)/2,mpoint=(m1+m2)/2,npoint=(n1+n2)/2
!
!  pencil-related stuff
!
  integer :: imn,m,n
  integer, dimension (ny*nz) :: mm,nn
  logical, dimension (ny*nz) :: necessary=.false.
!
!  in this section are all the things related to printing
!
  integer :: nname=0,nnamev=0,nnamez=0,nnamexy=0,nnamerz=0
  integer :: ilabel_max=-1,ilabel_sum=1,ilabel_save=0,ilabel_max_sqrt=-2,ilabel_sum_sqrt=2
  integer :: ilabel_max_dt=-3,ilabel_integrate=3,ilabel_surf=4
  integer :: nr_directions=1
  integer, parameter :: mname=100,mnamev=100,mnamez=20,mnamexy=6,mnamerz=20
  integer, dimension (mname) :: itype_name
  real, dimension (mname) :: fname
  real, dimension (nz,nprocz,mnamez) :: fnamez
  real, dimension (nx,ny,nprocy,mnamexy) :: fnamexy
  real, dimension (nrcyl,0:nz,nprocz,mnamerz) :: fnamerz
  real, dimension (nrcyl,nx) :: phiavg_profile
  real, dimension (nx) :: pomx,pomy,phix,phiy
  character (LEN=30) :: cname(mname),cform(mname)
  character (LEN=30) :: cnamev(mname)
  character (LEN=30) :: cnamexy(mnamexy),cformxy(mnamexy)
  character (LEN=30) :: cnamez(mnamez),cformz(mnamez)
  character (LEN=30) :: cnamerz(mnamerz),cformrz(mnamerz)

  ! other variables (needs to be consistent with reset list in register.90)
  integer :: i_t=0,i_it=0,i_dt=0,i_walltime=0
  integer :: i_rcylmphi=0,i_phimphi=0,i_zmphi=0,i_rmphi=0

  !  initialization of various switches; actual settings depends on the
  !  modules that are linked in (see Makefile.local) and can, in some cases,
  !  be reset also via appropriate namelist entries.

  logical :: lstart=.false., lrun=.false.
  logical :: lhydro=.true., ldensity=.true., lentropy=.false., lmagnetic=.false.
  logical :: lmpicomm=.false., lforcing=.false., lpostproc=.false.
  logical :: lmaxadvec_sum=.false.,old_cdtv=.true.
  logical :: lspecial=.false., lwrite_slices=.false., lwrite_2daverages=.false.
  logical :: lwrite_slice_xy2,lwrite_slice_xy,lwrite_slice_xz,lwrite_slice_yz
  logical :: lgrav=.false., lgravz=.false., lgravr=.false.
  logical :: lout,headt=.false.,headtt=.true.,ldt,lfirst,ldiagnos,lvid
  logical :: l2davg,l2davgfirst
  logical :: lwrite_zaverages=.true.,lwrite_phiaverages=.true.
  logical :: lwrite_ic=.false.,lnowrite=.false.,lserial_io=.false.
  logical :: lroot=.true.,ldebug=.false.,lfft=.true.
  logical :: lshear=.false.,lpscalar=.false.,lviscosity=.false.
  logical :: lradiation=.false.,lradiation_ray=.false.,lradiation_fld=.false.
  logical :: ldustdensity=.false.,ldustvelocity=.false.,linterstellar=.false.
  logical :: lvisc_shock=.false.,lvisc_hyper=.false.
  logical :: lcosmicray=.false.
  logical :: lselfgravity=.false.
  logical :: lmonolithic_io=.false.
  logical :: lionization=.false.,lionization_fixed=.false.

  ! variables to allow modules to share 'precalculated' stuff
  ! when necessary (set in module initialize functions)
  logical :: lneed_sij=.false., lneed_glnrho=.false.
  logical :: lneed_sdij=.false.

  logical :: lfirstpoint=.false.,llastpoint=.false.
  logical :: vel_spec=.false.,mag_spec=.false.,uxj_spec=.false.,vec_spec=.false.
  logical :: ro_spec=.false.,ss_spec=.false.,cc_spec=.false.,cr_spec=.false.
  logical :: ab_spec=.false.,ou_spec=.false.,oned=.false.
  logical :: rhocc_pdf=.false.,cc_pdf=.false.,lncc_pdf=.false.
  logical :: gcc_pdf=.false.,lngcc_pdf=.false.
  logical :: test_nonblocking=.false.,onedall=.false.
  logical :: lsfu=.false.,lsfb=.false.,lsfz1=.false.,lsfz2=.false.
  logical :: lsfflux=.false.
  logical :: lpdfu=.false.,lpdfb=.false.,lpdfz1=.false.,lpdfz2=.false.
!  logical, dimension(mvar + maux) :: lsnap ! flag which variables should be written
                                             ! to the snapshots

  character (len=2*bclen+1), dimension(mvar) :: bcx='p',bcy='p',bcz='p'
  character (len=bclen), dimension(mvar) :: bcx1,bcx2,bcy1,bcy2,bcz1,bcz2
  character (len=120) :: datadir='data' ! default; may be overwritten in
                                        ! Register.initialize()
  character (len=120) :: directory='',datadir_snap='',directory_snap=''
  character (len=120) :: cvsid='[No CVS Id given]'

  character (len=10), dimension(maux) :: aux_var
  integer :: aux_count=1


endmodule Cdata
