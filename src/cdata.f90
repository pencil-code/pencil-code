! $Id: cdata.f90,v 1.142 2003-06-10 19:25:32 mee Exp $

module Cdata

!!! Global variables

  Use Cparam

  integer :: itorder=3
  real, dimension (mx) :: x
  real, dimension (my) :: y
  real, dimension (mz) :: z
  real, dimension (nrcyl) :: rcyl  ! used for phi-averages
  real, dimension (nx) :: x_mn,y_mn,z_mn,r_mn,rcyl_mn
  real, dimension (nx) :: maxadvec2,maxdiffus, maxheating

  real, dimension (nx,3,3) :: sij,sdij  ! rate-of-strain tensor

  real, parameter :: pi=3.14159265358979324D0,epsi=5*epsilon(1.)
  real, dimension(3) :: Lxyz,xyz0,xyz1=impossible
  real :: t,dt=0.,cdt=0.4,cdtv=0.08,ttransient=0.
  real :: dx,dy,dz,dxmin,dxmax,drcyl
  real :: dsnap=100.,dvid=100.,dtmin=1.e-6,dspec=impossible
  real :: tsforce=-10., dtforce=10
  real, dimension (2) :: fran1,fran2

  !  units (need to be in double precision)
  character (len=3) :: unit_system='cgs'
  double precision :: unit_length=1.,unit_velocity=1.,unit_density=1.,unit_temperature=1.
  double precision :: k_B,m_p,m_e,eV,hbar,sigmaH_,sigmaSB,kappa_es

!ajwm nu moved to viscosity module
!ajwm replaced nu, causes error in forcing to resolve
  real :: nu=0.,cmu,cnu2,nud
  real :: tdiagnos,dtu
  real :: rmean,rrms,rmax,u2m,um2,u2max,divurms,divumax,divu2max
  real :: o2m,om2,oum
  real :: UUmax,x0,y0,z0,Lx,Ly,Lz
  real :: grads0=0.   ! (1/c_p)ds/dz
  real :: Omega=0.,qshear=0.,Sshear=impossible
  real :: deltay=0. !(for shear, now also used in forcing and output)

  integer, dimension(mseed) :: seed=0
  integer :: nseed=0
  integer :: nvar,naux,iuu=0,iux=0,iuy=0,iuz=0,ilnrho=0,ient=0
  integer :: iuud=0,iudx=0,iudy=0,iudz=0,ilnrhod=0,igg=0,igx=0,igy=0,igz=0
  integer :: iaa=0,iax=0,iay=0,iaz=0
  integer :: ie=0,iff=0,ifx=0,ify=0,ifz=0,idd=0, ishock=0
  integer :: nt=0,it1=10
  integer :: it,ix=l1,iy=m1,iz=n1,iz2=n2
  integer :: ilncc=0
  integer :: m,n
  integer :: iproc,ipx,ipy,ipz,root=0
  logical, dimension(3) :: lperi
  character (len=labellen) ::fft_switch='fftpack'

!
!  in this section are all the things related to printing
!
  integer :: nname=0,nnamez=0,nnamexy=0,nnamerz=0
  integer :: ilabel_max=-1,ilabel_sum=1,ilabel_save=0,ilabel_max_sqrt=-2,ilabel_sum_sqrt=2,ilabel_integrate=3
  integer :: nr_directions=1
  integer, parameter :: mname=100,mnamez=20,mnamexy=6,mnamerz=6
  integer, dimension (mname) :: itype_name
  real, dimension (mname) :: fname
  real, dimension (nz,nprocz,mnamez) :: fnamez
  real, dimension (nx,ny,nprocy,mnamexy) :: fnamexy
  real, dimension (nrcyl,0:nz,nprocz,mnamerz) :: fnamerz
  real, dimension (nrcyl,nx) :: phiavg_profile
  character (LEN=30) :: cname(mname),cform(mname)
  character (LEN=30) :: cnamexy(mnamexy),cformxy(mnamexy)
  character (LEN=30) :: cnamez(mnamez),cformz(mnamez)
  character (LEN=30) :: cnamerz(mnamerz),cformrz(mnamerz)

  ! other variables (needs to be consistent with reset list in register.90)
  integer :: i_t=0,i_it=0,i_dt=0,i_dtc=0,i_walltime=0

  !  initialization of various switches; actual settings depends on the
  !  modules that are linked in (see Makefile.local) and can, in some cases,
  !  be reset also via appropriate namelist entries.

  logical :: lhydro=.true., ldensity=.true., lentropy=.false., lmagnetic=.false.
  logical :: lmpicomm=.false., lforcing=.false., lpostproc=.false.
  logical :: lgrav=.false., lgravz=.false., lgravr=.false.
  logical :: lout,headt=.false.,headtt=.true.,ldt,lfirst,ldiagnos,lvid
  logical :: lwrite_ic=.false.,lnowrite=.false.,lserial_io=.false.
  logical :: lroot=.true.,ldebug=.false.,lfft=.true.
  logical :: lshear=.false.,lpscalar=.false.,lviscosity=.false.
  logical :: lradiation=.false.,lradiation_ray=.false.,lradiation_fld=.false.
  logical :: ldustdensity=.false.,ldustvelocity=.false.,linterstellar=.false.
  logical :: lselfgravity=.false.

  ! variables to allow modules to share 'precalculated' stuff
  ! when necessary (set in module initialize functions)
  logical :: lneed_sij=.false., lneed_glnrho=.false.
  logical :: lneed_sdij=.false.

  logical :: lfirstpoint
  logical :: vel_spec=.false.,mag_spec=.false.,vec_spec=.false.
  logical :: ab_spec=.false.,ou_spec=.false.,oned=.false.
  logical :: test_nonblocking=.false.
  logical :: lsfu=.false.,lsfb=.false.,lsfz1=.false.,lsfz2=.false.
  logical :: lpdfu=.false.,lpdfb=.false.,lpdfz1=.false.,lpdfz2=.false.
!  logical, dimension(mvar + maux) :: lsnap ! flag which variables should be written
                                             ! to the snapshots

  character (len=2*bclen+1), dimension(mvar) :: bcx,bcy,bcz
  character (len=bclen), dimension(mvar) :: bcx1,bcx2,bcy1,bcy2,bcz1,bcz2
  character (len=120) :: datadir='data' ! default; may be overwritten in
                                        ! Register.initialize()
  character (len=120) :: directory='',directory_snap=''
  character (len=120) :: cvsid='[No CVS Id given]'

endmodule Cdata
