! $Id: slices.f90,v 1.32 2003-11-21 01:54:06 brandenb Exp $

!  This module produces slices for animation purposes

module Slices

  use Cdata

  implicit none

  real, dimension (nx,ny,3) :: uu_xy,uu_xy2,uud_xy,uud_xy2,bb_xy,bb_xy2
  real, dimension (nx,ny,3) :: oo_xy,oo_xy2,aa_xy,aa_xy2
  real, dimension (nx,ny) :: lnrho_xy,lnrho_xy2,lnrhod_xy,lnrhod_xy2
  real, dimension (nx,ny) :: divu_xy,divu_xy2,o2_xy,o2_xy2,b2_xy,b2_xy2
  real, dimension (nx,ny) :: ss_xy,ss_xy2,lncc_xy,lncc_xy2
  real, dimension (nx,ny) :: lnTT_xy,lnTT_xy2,yH_xy,yH_xy2,ecr_xy,ecr_xy2
  real, dimension (nx,ny) :: Qrad_xy,Qrad_xy2,shock_xy,shock_xy2
  real, dimension (nx,ny) :: Isurf_xy

  real, dimension (nx,nz,3) :: uu_xz,uud_xz,bb_xz,oo_xz,aa_xz
  real, dimension (nx,nz) :: lnrho_xz,lnrhod_xz,ss_xz,lncc_xz,divu_xz
  real, dimension (nx,nz) :: lnTT_xz,yH_xz,ecr_xz,o2_xz,b2_xz
  real, dimension (nx,nz) :: Qrad_xz,shock_xz

  real, dimension (ny,nz,3) :: uu_yz,uud_yz,bb_yz,oo_yz,aa_yz
  real, dimension (ny,nz) :: lnrho_yz,lnrhod_yz,ss_yz,lncc_yz,divu_yz
  real, dimension (ny,nz) :: lnTT_yz,yH_yz,ecr_yz,o2_yz,b2_yz
  real, dimension (ny,nz) :: Qrad_yz,shock_yz

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
!  18-mar-03/axel: added dust velocity
!
      use Sub
!
      real, save :: tvid
      integer, save :: ifirst,nvid
!
      character (len=4) :: ch
      character (len=130) :: file
!
!  Output vid-data in 'tvid' time intervals
!
      file = trim(datadir)//'/tvid.dat'
      if (ifirst==0) then
        call read_snaptime(trim(file),tvid,nvid,dvid,t)
        ifirst=1
      endif
!
!  This routine sets lvid=T whenever its time to write a slice
!
      call update_snaptime(file,tvid,nvid,dvid,t,lvid,ch,ENUM=.false.)
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
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      character(len=*) :: path
      logical, save :: lfirstloop=.true.
      logical :: lnewfile=.true.
      integer :: inamev
!
!  Loop over slices
!
      do inamev=1,nnamev
      select case (trim(cnamev(inamev)))
!
!  Velocity field
!
      case ('uu')
        uu_yz=f(ix,m1:m2,n1:n2,iux:iuz)
        uu_xz=f(l1:l2,iy,n1:n2,iux:iuz)
        uu_xy=f(l1:l2,m1:m2,iz,iux:iuz)
        uu_xy2=f(l1:l2,m1:m2,iz2,iux:iuz)
        call wslice(path//'ux.yz',uu_yz(:,:,1),x(ix),ny,nz)
        call wslice(path//'uy.yz',uu_yz(:,:,2),x(ix),ny,nz)
        call wslice(path//'uz.yz',uu_yz(:,:,3),x(ix),ny,nz)
        call wslice(path//'ux.xz',uu_xz(:,:,1),y(iy),nx,nz)
        call wslice(path//'uy.xz',uu_xz(:,:,2),y(iy),nx,nz)
        call wslice(path//'uz.xz',uu_xz(:,:,3),y(iy),nx,nz)
        call wslice(path//'ux.xy',uu_xy(:,:,1),z(iz),nx,ny)
        call wslice(path//'uy.xy',uu_xy(:,:,2),z(iz),nx,ny)
        call wslice(path//'uz.xy',uu_xy(:,:,3),z(iz),nx,ny)
        call wslice(path//'ux.Xy',uu_xy2(:,:,1),z(iz2),nx,ny)
        call wslice(path//'uy.Xy',uu_xy2(:,:,2),z(iz2),nx,ny)
        call wslice(path//'uz.Xy',uu_xy2(:,:,3),z(iz2),nx,ny)
!
      case ('oo')
        call wslice(path//'ox.yz',oo_yz(:,:,1),x(ix),ny,nz)
        call wslice(path//'oy.yz',oo_yz(:,:,2),x(ix),ny,nz)
        call wslice(path//'oz.yz',oo_yz(:,:,3),x(ix),ny,nz)
        call wslice(path//'ox.xz',oo_xz(:,:,1),y(iy),nx,nz)
        call wslice(path//'oy.xz',oo_xz(:,:,2),y(iy),nx,nz)
        call wslice(path//'oz.xz',oo_xz(:,:,3),y(iy),nx,nz)
        call wslice(path//'ox.xy',oo_xy(:,:,1),z(iz),nx,ny)
        call wslice(path//'oy.xy',oo_xy(:,:,2),z(iz),nx,ny)
        call wslice(path//'oz.xy',oo_xy(:,:,3),z(iz),nx,ny)
        call wslice(path//'ox.Xy',oo_xy2(:,:,1),z(iz2),nx,ny)
        call wslice(path//'oy.Xy',oo_xy2(:,:,2),z(iz2),nx,ny)
        call wslice(path//'oz.Xy',oo_xy2(:,:,3),z(iz2),nx,ny)
!
      case ('o2')
        call wslice(path//'o2.yz',o2_yz,x(ix),ny,nz)
        call wslice(path//'o2.xz',o2_xz,y(iy),nx,nz)
        call wslice(path//'o2.xy',o2_xy,z(iz),nx,ny)
        call wslice(path//'o2.Xy',o2_xy2,z(iz2),nx,ny)
!
      case ('divu')
        call wslice(path//'divu.yz',divu_yz,x(ix),ny,nz)
        call wslice(path//'divu.xz',divu_xz,y(iy),nx,nz)
        call wslice(path//'divu.xy',divu_xy,z(iz),nx,ny)
        call wslice(path//'divu.Xy',divu_xy2,z(iz2),nx,ny)
!
!  Dust velocity
!
      case ('uud')
        uud_yz=f(ix,m1:m2,n1:n2,iudx:iudz)
        uud_xz=f(l1:l2,iy,n1:n2,iudx:iudz)
        uud_xy=f(l1:l2,m1:m2,iz,iudx:iudz)
        uud_xy2=f(l1:l2,m1:m2,iz2,iudx:iudz)
        call wslice(path//'udx.yz',uud_yz(:,:,1),x(ix),ny,nz)
        call wslice(path//'udy.yz',uud_yz(:,:,2),x(ix),ny,nz)
        call wslice(path//'udz.yz',uud_yz(:,:,3),x(ix),ny,nz)
        call wslice(path//'udx.xz',uud_xz(:,:,1),y(iy),nx,nz)
        call wslice(path//'udy.xz',uud_xz(:,:,2),y(iy),nx,nz)
        call wslice(path//'udz.xz',uud_xz(:,:,3),y(iy),nx,nz)
        call wslice(path//'udx.xy',uud_xy(:,:,1),z(iz),nx,ny)
        call wslice(path//'udy.xy',uud_xy(:,:,2),z(iz),nx,ny)
        call wslice(path//'udz.xy',uud_xy(:,:,3),z(iz),nx,ny)
        call wslice(path//'udx.Xy',uud_xy2(:,:,1),z(iz2),nx,ny)
        call wslice(path//'udy.Xy',uud_xy2(:,:,2),z(iz2),nx,ny)
        call wslice(path//'udz.Xy',uud_xy2(:,:,3),z(iz2),nx,ny)
!
!  Logarithmic density
!
      case ('lnrho')
        lnrho_yz=f(ix,m1:m2,n1:n2,ilnrho)
        lnrho_xz=f(l1:l2,iy,n1:n2,ilnrho)
        lnrho_xy=f(l1:l2,m1:m2,iz,ilnrho)
        lnrho_xy2=f(l1:l2,m1:m2,iz2,ilnrho)
        call wslice(path//'lnrho.yz',lnrho_yz,x(ix),ny,nz)
        call wslice(path//'lnrho.xz',lnrho_xz,y(iy),nx,nz)
        call wslice(path//'lnrho.xy',lnrho_xy,z(iz),nx,ny)
        call wslice(path//'lnrho.Xy',lnrho_xy2,z(iz2),nx,ny)
!
!  Logarithmic dust density
!
      case ('lnrhod')
        lnrhod_yz=f(ix,m1:m2,n1:n2,ilnrhod)
        lnrhod_xz=f(l1:l2,iy,n1:n2,ilnrhod)
        lnrhod_xy=f(l1:l2,m1:m2,iz,ilnrhod)
        lnrhod_xy2=f(l1:l2,m1:m2,iz2,ilnrhod)
        call wslice(path//'lnrhod.yz',lnrhod_yz,x(ix),ny,nz)
        call wslice(path//'lnrhod.xz',lnrhod_xz,y(iy),nx,nz)
        call wslice(path//'lnrhod.xy',lnrhod_xy,z(iz),nx,ny)
        call wslice(path//'lnrhod.Xy',lnrhod_xy2,z(iz2),nx,ny)
!
!  Entropy
!
      case ('ss')
        ss_yz=f(ix,m1:m2,n1:n2,iss)
        ss_xz=f(l1:l2,iy,n1:n2,iss)
        ss_xy=f(l1:l2,m1:m2,iz,iss)
        ss_xy2=f(l1:l2,m1:m2,iz2,iss)
        call wslice(path//'ss.yz',ss_yz,x(ix),ny,nz)
        call wslice(path//'ss.xz',ss_xz,y(iy),nx,nz)
        call wslice(path//'ss.xy',ss_xy,z(iz),nx,ny)
        call wslice(path//'ss.Xy',ss_xy2,z(iz2),nx,ny)
!
!  Shock viscosity
!
      case ('shock')
        shock_yz=f(ix,m1:m2,n1:n2,ishock)
        shock_xz=f(l1:l2,iy,n1:n2,ishock)
        shock_xy=f(l1:l2,m1:m2,iz,ishock)
        shock_xy2=f(l1:l2,m1:m2,iz2,ishock)
        call wslice(path//'shock.yz',shock_yz,x(ix),ny,nz)
        call wslice(path//'shock.xz',shock_xz,y(iy),nx,nz)
        call wslice(path//'shock.xy',shock_xy,z(iz),nx,ny)
        call wslice(path//'shock.Xy',shock_xy2,z(iz2),nx,ny)
!
!  Temperature
!
      case ('lnTT')
        lnTT_yz=f(ix,m1:m2,n1:n2,ilnTT)
        lnTT_xz=f(l1:l2,iy,n1:n2,ilnTT)
        lnTT_xy=f(l1:l2,m1:m2,iz,ilnTT)
        lnTT_xy2=f(l1:l2,m1:m2,iz2,ilnTT)
        call wslice(path//'lnTT.yz',lnTT_yz,x(ix),ny,nz)
        call wslice(path//'lnTT.xz',lnTT_xz,y(iy),nx,nz)
        call wslice(path//'lnTT.xy',lnTT_xy,z(iz),nx,ny)
        call wslice(path//'lnTT.Xy',lnTT_xy2,z(iz2),nx,ny)
!
!  Degree of ionization
!
      case ('yH')
        yH_yz=f(ix,m1:m2,n1:n2,iyH)
        yH_xz=f(l1:l2,iy,n1:n2,iyH)
        yH_xy=f(l1:l2,m1:m2,iz,iyH)
        yH_xy2=f(l1:l2,m1:m2,iz2,iyH)
        call wslice(path//'yH.yz',yH_yz,x(ix),ny,nz)
        call wslice(path//'yH.xz',yH_xz,y(iy),nx,nz)
        call wslice(path//'yH.xy',yH_xy,z(iz),nx,ny)
        call wslice(path//'yH.Xy',yH_xy2,z(iz2),nx,ny)
!
!  Heating rate
!
      case ('Qrad')
        Qrad_yz=f(ix,m1:m2,n1:n2,iQrad)
        Qrad_xz=f(l1:l2,iy,n1:n2,iQrad)
        Qrad_xy=f(l1:l2,m1:m2,iz,iQrad)
        Qrad_xy2=f(l1:l2,m1:m2,iz2,iQrad)
        call wslice(path//'Qrad.yz',Qrad_yz,x(ix),ny,nz)
        call wslice(path//'Qrad.xz',Qrad_xz,y(iy),nx,nz)
        call wslice(path//'Qrad.xy',Qrad_xy,z(iz),nx,ny)
        call wslice(path//'Qrad.Xy',Qrad_xy2,z(iz2),nx,ny)
!
!  Surface intensity
!
      case ('Isurf')
        call wslice(path//'Isurf.xy',Isurf_xy,z(iz2),nx,ny)
!
!  Magnetic field
!
      case ('bb')
        call wslice(path//'bx.yz',bb_yz(:,:,1),x(ix),ny,nz)
        call wslice(path//'by.yz',bb_yz(:,:,2),x(ix),ny,nz)
        call wslice(path//'bz.yz',bb_yz(:,:,3),x(ix),ny,nz)
        call wslice(path//'bx.xz',bb_xz(:,:,1),y(iy),nx,nz)
        call wslice(path//'by.xz',bb_xz(:,:,2),y(iy),nx,nz)
        call wslice(path//'bz.xz',bb_xz(:,:,3),y(iy),nx,nz)
        call wslice(path//'bx.xy',bb_xy(:,:,1),z(iz),nx,ny)
        call wslice(path//'by.xy',bb_xy(:,:,2),z(iz),nx,ny)
        call wslice(path//'bz.xy',bb_xy(:,:,3),z(iz),nx,ny)
        call wslice(path//'bx.Xy',bb_xy2(:,:,1),z(iz2),nx,ny)
        call wslice(path//'by.Xy',bb_xy2(:,:,2),z(iz2),nx,ny)
        call wslice(path//'bz.Xy',bb_xy2(:,:,3),z(iz2),nx,ny)
!
      case ('b2')
        call wslice(path//'b2.yz',b2_yz,x(ix),ny,nz)
        call wslice(path//'b2.xz',b2_xz,y(iy),nx,nz)
        call wslice(path//'b2.xy',b2_xy,z(iz),nx,ny)
        call wslice(path//'b2.Xy',b2_xy2,z(iz2),nx,ny)
!
      case ('aa')
        aa_yz=f(ix,m1:m2,n1:n2,iax:iaz)
        aa_xz=f(l1:l2,iy,n1:n2,iax:iaz)
        aa_xy=f(l1:l2,m1:m2,iz,iax:iaz)
        aa_xy2=f(l1:l2,m1:m2,iz2,iax:iaz)
        call wslice(path//'ax.yz',aa_yz(:,:,1),x(ix),ny,nz)
        call wslice(path//'ay.yz',aa_yz(:,:,2),x(ix),ny,nz)
        call wslice(path//'az.yz',aa_yz(:,:,3),x(ix),ny,nz)
        call wslice(path//'ax.xz',aa_xz(:,:,1),y(iy),nx,nz)
        call wslice(path//'ay.xz',aa_xz(:,:,2),y(iy),nx,nz)
        call wslice(path//'az.xz',aa_xz(:,:,3),y(iy),nx,nz)
        call wslice(path//'ax.xy',aa_xy(:,:,1),z(iz),nx,ny)
        call wslice(path//'ay.xy',aa_xy(:,:,2),z(iz),nx,ny)
        call wslice(path//'az.xy',aa_xy(:,:,3),z(iz),nx,ny)
        call wslice(path//'ax.Xy',aa_xy2(:,:,1),z(iz2),nx,ny)
        call wslice(path//'ay.Xy',aa_xy2(:,:,2),z(iz2),nx,ny)
        call wslice(path//'az.Xy',aa_xy2(:,:,3),z(iz2),nx,ny)
!
!  Passive scalar
!
      case ('lncc')
        lncc_yz=f(ix,m1:m2,n1:n2,ilncc)
        lncc_xz=f(l1:l2,iy,n1:n2,ilncc)
        lncc_xy=f(l1:l2,m1:m2,iz,ilncc)
        lncc_xy2=f(l1:l2,m1:m2,iz2,ilncc)
        call wslice(path//'lncc.yz',lncc_yz,x(ix),ny,nz)
        call wslice(path//'lncc.xz',lncc_xz,y(iy),nx,nz)
        call wslice(path//'lncc.xy',lncc_xy,z(iz),nx,ny)
        call wslice(path//'lncc.Xy',lncc_xy2,z(iz2),nx,ny)
!
!  Cosmic ray energy density
!
      case ('ecr')
        ecr_yz=f(ix,m1:m2,n1:n2,iecr)
        ecr_xz=f(l1:l2,iy,n1:n2,iecr)
        ecr_xy=f(l1:l2,m1:m2,iz,iecr)
        ecr_xy2=f(l1:l2,m1:m2,iz2,iecr)
        call wslice(path//'ecr.yz',ecr_yz,x(ix),ny,nz)
        call wslice(path//'ecr.xz',ecr_xz,y(iy),nx,nz)
        call wslice(path//'ecr.xy',ecr_xy,z(iz),nx,ny)
        call wslice(path//'ecr.Xy',ecr_xy2,z(iz2),nx,ny)
!
      case default
        if (lfirstloop.and.lroot) then
          if (lnewfile) then
            open(1,file='video.err')
            lnewfile=.false.
          else
            open(1,file='video.err',position='append')
          endif
          write(1,*) 'unknown slice: ',trim(cnamev(inamev))
          close(1)
        endif
!
      end select
      enddo

      lfirstloop=.false.
!
    endsubroutine wvid
!***********************************************************************
    subroutine wslice(file,a,pos,ndim1,ndim2)
!
!  appending to an existing slice file
!
!  12-nov-02/axel: coded
!
      integer :: ndim1,ndim2
      character (len=*) :: file
      real, dimension (ndim1,ndim2) :: a
      real, intent(in) :: pos
!
      open(1,file=file,form='unformatted',position='append')
      write(1) a,t,pos
      close(1)
!
    endsubroutine wslice
!***********************************************************************

endmodule Slices
