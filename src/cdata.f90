module Cdata

!!! Global variables

  Use Cparam

  real, dimension (mx) :: x
  real, dimension (my) :: y
  real, dimension (mz) :: z
  real, dimension (nx) :: rmn
  real :: t,dt,dx,dy,dz,dxmin,dxmax
  real :: dsnap,dvid,dforce
  real :: tinit,tdamp,dampu,dampuext,rdamp,wdamp
  real :: cs0,rho0,cs20,gamma,gamma1,force,relhel
  real :: DD,nu,cmu,cnu2,cdiffrho,hcond0,hcond1,hcond2,whcond
  real :: t_diag,rmean,rrms,rmax,urms,umax,u2max,divurms,divumax,divu2max
  real :: orms,omax,o2max,ourms,oumax
  real :: UUmax,viscmax,cdt,cdtv,x0,y0,z0,Lx,Ly,Lz
  real :: z1,z2,z3
  real :: gravz,ss0,grads0      ! (1/c_p)ds/dz
  real :: urand,cheat,wheat,cool,wcool
  real, parameter :: pi=3.14159265358979323844,epsi=5*epsilon(1.)

  integer, dimension (2) :: seed
  integer :: nvar,iuu,iux,iuy,iuz,ilnrho,ient,iaa,iax,iay,iaz
  integer :: iperx,ipery,iperz
  integer :: nt,it1,isave,itorder
  integer :: it,ix,iy,iz
  integer :: ivisc,iforce
  integer :: m,n

  logical :: lmpicomm=.false., lentropy=.false., lmagnetic=.false.
  logical :: lgrav=.false., lgravz=.false., lgravr=.false.
  logical :: lout,headt,headtt,ldt,lfirst
  logical :: lroot=.true.
  logical :: lfirstpoint

  character (LEN=80) :: form1
  character (LEN=2*bclen+1), dimension(mvar) :: bcx,bcy,bcz
  character (LEN=bclen), dimension(mvar) :: bcx1,bcx2,bcy1,bcy2,bcz1,bcz2

endmodule Cdata
