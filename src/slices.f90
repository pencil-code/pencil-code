! $Id: slices.f90,v 1.7 2002-11-13 10:08:52 brandenb Exp $

!  This module produces slices for animation purposes

module Slices

  use Cdata

  implicit none

  real, dimension (nx,ny,3) :: uu_xy,uu_xy2,bb_xy,bb_xy2
  real, dimension (nx,ny) :: lnrho_xy,lnrho_xy2,ss_xy,ss_xy2

  real, dimension (nx,nz,3) :: uu_xz,bb_xz
  real, dimension (nx,nz) :: lnrho_xz,ss_xz

  real, dimension (ny,nz,3) :: uu_yz,bb_yz
  real, dimension (ny,nz) :: lnrho_yz,ss_yz

  contains

!***********************************************************************
    subroutine wvid_prepare
!
!  Prepare lvid for writing slices into video file
!  This is useful for visualization of scalar field (or one component
!  of a vector field) on the perifery of a box.
!  Can be visualized in idl using rvid_box.pro
!
!  20-oct-97/axel: coded
!  08-oct-02/tony: increased size of file to handle datadir//'/tvid.dat'
!  13-nov-02/axel: added more fields, use wslice.
!
      use Cdata
!
      real, save :: tvid
      integer, save :: nvid
!
      character (len=4) :: ch
      character (len=130) :: file
!
!  Output vid-data in 'tvid' time intervals
!  Here we read the tvid.dat file every time again.
!  This allows us to change its content on the fly.
!
      file = trim(datadir)//'/tvid.dat'
      call out1 (trim(file),tvid,nvid,dvid,t)
      call out2 (trim(file),tvid,nvid,dvid,t,lvid,ch,.false.)
!
    endsubroutine wvid_prepare
!***********************************************************************
    subroutine wvid(f,path)
!
!  Write slices for animation of scalar field
! (or one component of a vector field) on the perifery of a box.
!  Can be visualized in idl using rvid_box.pro.
!  Data not derived from f (e.g. magnetic field) are prepared in pde.
!
!  13-nov-02/axel: added more fields, use wslice.
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      character (len=*) :: path
      integer :: j
!
!  Velocity field
!
      if (lhydro) then
        do j=1,3
          uu_yz(:,:,j)=f(ix,m1:m2,n1:n2,j+iuu-1)
          uu_xz(:,:,j)=f(l1:l2,iy,n1:n2,j+iuu-1)
          uu_xy(:,:,j)=f(l1:l2,m1:m2,iz,j+iuu-1)
          uu_xy2(:,:,j)=f(l1:l2,m1:m2,iz2,j+iuu-1)
        enddo
        call wslice(path//'ux.yz',uu_yz(:,:,1),ny,nz)
        call wslice(path//'uy.yz',uu_yz(:,:,2),ny,nz)
        call wslice(path//'uz.yz',uu_yz(:,:,3),ny,nz)
        call wslice(path//'ux.xz',uu_xz(:,:,1),nx,nz)
        call wslice(path//'uy.xz',uu_xz(:,:,2),nx,nz)
        call wslice(path//'uz.xz',uu_xz(:,:,3),nx,nz)
        call wslice(path//'ux.xy',uu_xy(:,:,1),nx,ny)
        call wslice(path//'uy.xy',uu_xy(:,:,2),nx,ny)
        call wslice(path//'uz.xy',uu_xy(:,:,3),nx,ny)
        call wslice(path//'ux.Xy',uu_xy2(:,:,1),nx,ny)
        call wslice(path//'uy.Xy',uu_xy2(:,:,2),nx,ny)
        call wslice(path//'uz.Xy',uu_xy2(:,:,3),nx,ny)
      endif
!
!  logarithmic density
!
      if (ldensity) then
        lnrho_yz=f(ix,m1:m2,n1:n2,ilnrho)
        lnrho_xz=f(l1:l2,iy,n1:n2,ilnrho)
        lnrho_xy=f(l1:l2,m1:m2,iz,ilnrho)
        lnrho_xy2=f(l1:l2,m1:m2,iz2,ilnrho)
        call wslice(path//'lnrho.yz',lnrho_yz,ny,nz)
        call wslice(path//'lnrho.xz',lnrho_xz,nx,nz)
        call wslice(path//'lnrho.xy',lnrho_xy,nx,ny)
        call wslice(path//'lnrho.Xy',lnrho_xy2,nx,ny)
      endif
!
!  Entropy
!
      if (lentropy) then
        ss_yz=f(ix,m1:m2,n1:n2,ient)
        ss_xz=f(l1:l2,iy,n1:n2,ient)
        ss_xy=f(l1:l2,m1:m2,iz,ient)
        ss_xy2=f(l1:l2,m1:m2,iz2,ient)
        call wslice(path//'ss.yz',ss_yz,ny,nz)
        call wslice(path//'ss.xz',ss_xz,nx,nz)
        call wslice(path//'ss.xy',ss_xy,nx,ny)
        call wslice(path//'ss.Xy',ss_xy2,nx,ny)
      endif
!
!  Magnetic field
!
      if (lmagnetic) then
        call wslice(path//'bx.yz',bb_yz(:,:,1),ny,nz)
        call wslice(path//'by.yz',bb_yz(:,:,2),ny,nz)
        call wslice(path//'bz.yz',bb_yz(:,:,3),ny,nz)
        call wslice(path//'bx.xz',bb_xz(:,:,1),nx,nz)
        call wslice(path//'by.xz',bb_xz(:,:,2),nx,nz)
        call wslice(path//'bz.xz',bb_xz(:,:,3),nx,nz)
        call wslice(path//'bx.xy',bb_xy(:,:,1),nx,ny)
        call wslice(path//'by.xy',bb_xy(:,:,2),nx,ny)
        call wslice(path//'bz.xy',bb_xy(:,:,3),nx,ny)
        call wslice(path//'bx.Xy',bb_xy2(:,:,1),nx,ny)
        call wslice(path//'by.Xy',bb_xy2(:,:,2),nx,ny)
        call wslice(path//'bz.Xy',bb_xy2(:,:,3),nx,ny)
      endif
!
    endsubroutine wvid
!***********************************************************************
    subroutine wslice(file,a,ndim1,ndim2)
!
!  appending to an existing slice file
!
!  12-nov-02/axel: coded
!
      integer :: ndim1,ndim2
      character (len=*) :: file
      real, dimension (ndim1,ndim2) :: a
!
      open(1,file=file,form='unformatted',position='append')
      write(1) a,t
      close(1)
!
    endsubroutine wslice
!***********************************************************************

endmodule Slices
