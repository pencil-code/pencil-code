! $Id: cdata.f90,v 1.41 2002-06-05 23:45:57 brandenb Exp $

module Cdata

!!! Global variables

  Use Cparam

  real, dimension (mx) :: x
  real, dimension (my) :: y
  real, dimension (mz) :: z
  real, dimension (nx) :: x_mn,y_mn,z_mn,r_mn
!  real, dimension (nx) :: rmn
!
  real, parameter :: pi=3.14159265358979323844,epsi=5*epsilon(1.)
  real, dimension(3) :: xyz0,Lxyz
  real :: t,dt=0.,cdt=0.,cdtv=0.,dx,dy,dz,dxmin,dxmax
  real :: dsnap,dvid,dtmin=0.
  real :: tinit=0.,tdamp=0.,dampu=0.,dampuext=0.,rdamp=1.2,wdamp=0.2
  real :: DD,nu=0.,cmu,cnu2
  real :: tdiagnos,dtu
  real :: rmean,rrms,rmax,u2m,um2,u2max,divurms,divumax,divu2max
  real :: o2m,om2,oum
  real :: UUmax,x0,y0,z0,Lx,Ly,Lz
  real :: z1,z2,ztop
  real :: gravz,ss0,grads0      ! (1/c_p)ds/dz
  real :: urand

! These are parameters of Entropy, but appear in Boundcond and (worse) in
! wparam (Sub) as well, so they need to be declared here
  real :: hcond0,hcond1,hcond2,whcond
  real :: mpoly0,mpoly1,mpoly2
  real :: cheat,wheat,cool,wcool,Fheat
  integer:: isothtop

  integer, dimension (2) :: seed
  integer :: nvar,iuu=0,iux=0,iuy=0,iuz=0,ilnrho=0,ient=0,iaa=0,iax=0,iay=0,iaz=0
  integer :: nt,it1,isave,itorder
  integer :: it,ix,iy,iz
  integer :: ivisc
  integer :: m,n
  integer :: iproc,ipx,ipy,ipz,root=0
  logical, dimension(3) :: lperi

!
!  in this section are all the things related to printing
!
  integer :: nname,nnamez
  integer :: ilabel_max=-1,ilabel_sum=1,ilabel_save=0
  integer, parameter :: mname=100,mnamez=20
  integer, dimension (mname) :: itype_name
  real, dimension (mname) :: fname
  real, dimension (nz,mnamez) :: fnamez
  character (len=30) :: cname(mname),cform(mname),cnamez(mnamez),cformz(mnamez)

  logical :: lhydro=.true., ldensity=.true., lentropy=.false., lmagnetic=.false.
  logical :: lmpicomm=.false., lforcing=.false.
  logical :: lgrav=.false., lgravz=.false., lgravr=.false.
  logical :: lout,headt,headtt,ldt,lfirst,ldiagnos
  logical :: lroot=.true.
  logical :: lfirstpoint

  character (len=2*bclen+1), dimension(mvar) :: bcx,bcy,bcz
  character (len=bclen), dimension(mvar) :: bcx1,bcx2,bcy1,bcy2,bcz1,bcz2
  character (len=12) :: directory

endmodule Cdata
