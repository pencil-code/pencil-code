module Cdata

!!! Global variables

  Use Cparam

  real, dimension (mx) :: x
  real, dimension (my) :: y
  real, dimension (mz) :: z
  real :: t,dt,dx,dy,dz &
         ,dsnap,dvid,dforce &
         ,cs,cs20,gamma,force,relhel &
         ,DD,nu,cmu,cnu2 &
         ,rrms,rmax,urms,umax,divurms,divumax &
         ,orms,omax,ourms,oumax &
         ,UUmax,cdt,pi,Lx,Ly,Lz &
         ,ibc(mvar)

  integer :: nt,it1,isave,itorder
  integer :: it,ix,iy,iz
  integer :: ivisc,iforce
  integer :: m,n,im,in

  logical :: lout,headt,ldt,lfirst

  integer :: iu=1,ilnrho=4
  integer, dimension (2) :: seed

  character*80 :: form1

endmodule Cdata
