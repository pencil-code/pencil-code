! $Id: cdata.f90,v 1.122 2003-01-24 01:04:51 nilshau Exp $

module Cdata

!!! Global variables

  Use Cparam

  integer :: itorder=3
  real, dimension (mx) :: x
  real, dimension (my) :: y
  real, dimension (mz) :: z
  real, dimension (nx) :: x_mn,y_mn,z_mn,r_mn
  real, dimension (nx) :: maxadvec2,maxdiffus, maxheating

  real, dimension (nx,3,3) :: sij  ! rate-of-strain tensor

  real, parameter :: pi=3.14159265358979324D0,epsi=5*epsilon(1.)
  real, dimension(3) :: Lxyz,xyz0,xyz1=impossible
  real :: t,dt=0.,cdt=0.4,cdtv=0.08,ttransient=0.
  real :: dx,dy,dz,dxmin,dxmax
  real :: dsnap=100.,dvid=100.,dtmin=1.e-6,dspec=impossible
  real :: tsforce=-10., dtforce=1
  real, dimension (2) :: fran1,fran2
!ajwm nu moved to viscosity module
!ajwm replaced nu, causes error in forcing to resolve
  real :: nu=0.,cmu,cnu2
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
  integer :: iaa=0,iax=0,iay=0,iaz=0
  integer :: ie=0,iff=0,ifx=0,ify=0,ifz=0,idd=0, ishock=0
  integer :: nt=0,it1=10
  integer :: it,ix=l1,iy=m1,iz=n1,iz2=n2
  integer :: m,n
  integer :: iproc,ipx,ipy,ipz,root=0
  logical, dimension(3) :: lperi
  character (len=labellen) ::fft_switch='Singleton'

!
!  in this section are all the things related to printing
!
  integer :: nname=0,nnamez=0,nnamexy=0
  integer :: ilabel_max=-1,ilabel_sum=1,ilabel_save=0,ilabel_max_sqrt=-2,ilabel_sum_sqrt=2
  integer :: nr_directions=1
  integer, parameter :: mname=100,mnamez=20,mnamexy=6,mnamerz=6
  integer, dimension (mname) :: itype_name
  real, dimension (mname) :: fname
  real, dimension (nz,nprocz,mnamez) :: fnamez
  real, dimension (nx,ny,nprocy,mnamexy) :: fnamexy
  real, dimension (nx/2,nz,nprocz,mnamerz) :: fnamerz
  character (len=30) :: cname(mname),cform(mname),cnamez(mnamez),cformz(mnamez)
  character (len=30) :: cnamexy(mnamexy),cformxy(mnamexy)

  ! other variables (needs to be consistent with reset list in register.90)
  integer :: i_t=0,i_it=0,i_dt=0,i_dtc=0

  logical :: lhydro=.true., ldensity=.true., lentropy=.false., lmagnetic=.false.
  logical :: lmpicomm=.false., lforcing=.false., lpostproc=.false.
  logical :: lgrav=.false., lgravz=.false., lgravr=.false.
  logical :: lout,headt=.false.,headtt=.true.,ldt,lfirst,ldiagnos,lvid
  logical :: lwrite_ic=.false.,lnowrite=.false.,lserial_io=.false.
  logical :: lroot=.true.,ldebug=.false.,lfft=.true.
  logical :: lshear=.false.,lpscalar=.false.,lradiation=.false., lviscosity=.false.
  logical :: linterstellar=.false.

  ! variables to allow modules to share 'precalculated' stuff
  ! when necessary (set in module initialize functions)
  logical :: lneed_sij=.false., lneed_glnrho=.false.

  logical :: lfirstpoint
  logical :: vel_spec=.false.,mag_spec=.false.,vec_spec=.false.
  logical :: ab_spec=.false.,ou_spec=.false.
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
