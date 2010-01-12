!!!!!!!!!!!!!!!!!!!!!!!
!!!   minmask.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!

!!!  Author: wd (Wolfgang.Dobler@ucalgary.ca)
!!!  Date:   10-Apr-2005
!!!
!!!  Description:
!!!   Using minval(.., MASK=..) seems to fail with ifort 8.0

program Toto

  implicit none

  integer, parameter :: nxgrid=128, nygrid=128, nzgrid=1
  real :: dx=1., dy=1., dz=1.
  real :: dxmin

  dxmin = minval( (/dx,dy,dz,huge(dx)/), &
      MASK=((/nxgrid,nygrid,nzgrid,2/) > 1) )

  print*, 'dxmin = ', dxmin

endprogram Toto

!!! End of file minmask.f90

