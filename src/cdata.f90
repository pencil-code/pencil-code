module Cdata

!!! Global variables

  Use Cparam

  real, dimension (mx) :: x
  real, dimension (my) :: y
  real, dimension (mz) :: z
  real :: t,dt,dx,dy,dz
  real :: dsnap,dvid,dforce
  real :: cs,cs20,gamma,force,relhel
  real :: DD,nu,cmu,cnu2
  real :: rrms,rmax,urms,umax,divurms,divumax
  real :: orms,omax,ourms,oumax
  real :: UUmax,cdt,pi,Lx,Ly,Lz
  real :: gravz

  integer :: nvar,iuu,iux,iuy,iuz,ilnrho,ient,iaa,iax,iay,iaz
  integer :: nt,it1,isave,itorder
  integer :: it,ix,iy,iz
  integer :: ivisc,iforce
  integer :: ibc(mvar)
  integer :: m,n,im,in
  integer, dimension (2) :: seed

  logical :: lmpicomm=.false., lentropy=.false., lmagnetic=.false.
  logical :: lout,headt,headtt,ldt,lfirst
  logical :: lroot=.true.

  character*80 :: form1

endmodule Cdata
