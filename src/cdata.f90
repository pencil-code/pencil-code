module Cdata

!!! Global variables

  Use Cparam

  real, dimension (mx) :: x
  real, dimension (my) :: y
  real, dimension (mz) :: z
  real :: t,dt,dx,dy,dz,dxmin,dxmax
  real :: dsnap,dvid,dforce
  real :: cs0,rho0,cs20,gamma,gamma1,force,relhel
  real :: DD,nu,cmu,cnu2,cdiffrho
  real :: t_diag,rmean,rrms,rmax,urms,umax,u2max,divurms,divumax,divu2max
  real :: orms,omax,o2max,ourms,oumax
  real :: UUmax,cdt,pi,Lx,Ly,Lz
  real :: gravz

  integer :: nvar,iuu,iux,iuy,iuz,ilnrho,ient,iaa,iax,iay,iaz
  integer :: nt,it1,isave,itorder
  integer :: it,ix,iy,iz
  integer :: ivisc,iforce
  integer :: ibc(mvar)
  integer :: m,n
  integer, dimension (2) :: seed

  logical :: lmpicomm=.false., lentropy=.false., lmagnetic=.false.
  logical :: lgrav=.false., lgravz=.false., lgravr=.false.
  logical :: lout,headt,headtt,ldt,lfirst
  logical :: lroot=.true.
  logical :: lfirstpoint

  character (LEN=80) :: form1

endmodule Cdata
