! $Id: cdata.f90,v 1.72 2002-07-12 17:40:59 brandenb Exp $

module Cdata

!!! Global variables

  Use Cparam

  integer :: itorder=3
  real, dimension (mx) :: x
  real, dimension (my) :: y
  real, dimension (mz) :: z
  real, dimension (nx) :: x_mn,y_mn,z_mn,r_mn
  real, dimension (nx) :: maxadvec2,maxdiffus
 
  real, parameter :: pi=3.14159265358979323844,epsi=5*epsilon(1.)
  real, dimension(3) :: xyz0,Lxyz
  real :: t,dt=0.,cdt=0.4,cdtv=0.08,dx,dy,dz,dxmin,dxmax
  real :: dsnap=100.,dvid=100.,dtmin=0.
  real :: DD,nu=0.,cmu,cnu2
  real :: tdiagnos,dtu
  real :: rmean,rrms,rmax,u2m,um2,u2max,divurms,divumax,divu2max
  real :: o2m,om2,oum
  real :: UUmax,x0,y0,z0,Lx,Ly,Lz
  real :: grads0=0.   ! (1/c_p)ds/dz

! These are parameters of Entropy, but appear in Boundcond and (worse) in
! wparam (Sub) as well, so they need to be declared here
!AB: Are you sure this cannot be avoided?? It is no longer in Boundcond!!
  real :: hcond0=0,hcond1=0,hcond2=0,whcond=2*epsi
  real :: mpoly0,mpoly1,mpoly2
  real :: cheat,wheat,cool=0.,wcool,Fheat
  integer:: isothtop

  real :: Omega=0.,qshear=0.

  integer, dimension(mseed) :: seed=0
  integer :: nseed
  integer :: nvar,iuu=0,iux=0,iuy=0,iuz=0,ilnrho=0,ient=0
  integer :: iaa=0,iax=0,iay=0,iaz=0
  integer :: nt=0,it1=10
  integer :: it,ix=(mx+1)/2,iy=(my+1)/2,iz=(mz+1)/2
  integer :: ivisc
  integer :: m,n
  integer :: iproc,ipx,ipy,ipz,root=0
  logical, dimension(3) :: lperi
!
!  in this section are all the things related to printing
!
  integer :: nname=0,nnamez=0,nnamexy=0
  integer :: ilabel_max=-1,ilabel_sum=1,ilabel_save=0,ilabel_max_sqrt=-2,ilabel_sum_sqrt=2
  integer, parameter :: mname=100,mnamez=20,mnamexy=5
  integer, dimension (mname) :: itype_name
  real, dimension (mname) :: fname
  real, dimension (nz,nprocz,mnamez) :: fnamez
  real, dimension (nx,ny,nprocy,mnamexy) :: fnamexy
  character (len=30) :: cname(mname),cform(mname),cnamez(mnamez),cformz(mnamez)
  character (len=30) :: cnamexy(mnamexy),cformxy(mnamexy)

  ! other variables (needs to be consistent with reset list in register.90)
  integer :: i_t=0,i_it=0,i_dt=0,i_dtc=0

  logical :: lhydro=.true., ldensity=.true., lentropy=.false., lmagnetic=.false.
  logical :: lmpicomm=.false., lforcing=.false.
  logical :: lgrav=.false., lgravz=.false., lgravr=.false.
  logical :: lout,headt=.true.,headtt=.true.,ldt,lfirst,ldiagnos
  logical :: lwrite_ic=.false.,lnowrite=.false.
  logical :: lroot=.true.,ldebug=.false.
  logical :: lshear=.false.,lpscalar=.false.
  logical :: lfirstpoint

  character (len=2*bclen+1), dimension(mvar) :: bcx,bcy,bcz
  character (len=bclen), dimension(mvar) :: bcx1,bcx2,bcy1,bcy2,bcz1,bcz2
  character (len=12) :: directory
  character (len=120) :: cvsid='[No CVS Id given]'

endmodule Cdata
