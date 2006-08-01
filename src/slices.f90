! $Id: slices.f90,v 1.65 2006-08-01 19:02:48 mee Exp $

!  This module produces slices for animation purposes

module Slices

  use Cdata

  implicit none

  private

  public :: wvid, wvid_prepare, setup_slices, wslice

!  Variables for xy slices start here
!!! New slice code reuses the following slice variables.
  real, target, dimension (nx,ny) :: slice_xy = 0., slice_xy2 = 0.
  real, target, dimension (nx,nz) :: slice_xz = 0.
  real, target, dimension (ny,nz) :: slice_yz = 0.

!!! LEGACY SLICE VARIABLES FOLLOW
!  Code variables
  real, public, dimension (nx,ny,3) :: uu_xy,aa_xy,uud_xy,vvp_xy
  real, public, dimension (nx,ny) :: lnrho_xy,ss_xy,cc_xy,lncc_xy
  real, public, dimension (nx,ny) :: XX_chiral_xy,YY_chiral_xy
  real, public, dimension (nx,ny) :: nd_xy
  real, public, dimension (nx,ny,ndustspec) :: md_xy
!  Auxiliary variables  
  real, public, dimension (nx,ny) :: yH_xy,shock_xy,Qrad_xy,ecr_xy
  real, public, dimension (nx,ny,3) :: Frad_xy
!  Derived variables
  real, public, dimension (nx,ny,3) :: oo_xy,bb_xy
  real, public, dimension (nx,ny) :: divu_xy,u2_xy,o2_xy,lnTT_xy,b2_xy,jb_xy,Isurf_xy
  real, public, dimension (nx,ny) :: DQ_chiral_xy,QQ_chiral_xy
  real, public, dimension (nx,ny) :: epsd_xy
  real, public, dimension (nx,ny) :: pp_xy
!  Variables for xy2 slices start here
!  Code variables
  real, public, dimension (nx,ny,3) :: uu_xy2,aa_xy2,uud_xy2,vvp_xy2
  real, public, dimension (nx,ny) :: lnrho_xy2,ss_xy2,cc_xy2,lncc_xy2
  real, public, dimension (nx,ny) :: XX_chiral_xy2,YY_chiral_xy2
  real, public, dimension (nx,ny) :: nd_xy2
  real, public, dimension (nx,ny,ndustspec) :: md_xy2
!  Auxiliary variables  
  real, public, dimension (nx,ny) :: yH_xy2,shock_xy2,Qrad_xy2,ecr_xy2
  real, public, dimension (nx,ny,3) :: Frad_xy2
!  Derived variables
  real, public, dimension (nx,ny,3) :: oo_xy2,bb_xy2
  real, public, dimension (nx,ny) :: divu_xy2,u2_xy2,o2_xy2,lnTT_xy2,b2_xy2,jb_xy2
  real, public, dimension (nx,ny) :: DQ_chiral_xy2,QQ_chiral_xy2
  real, public, dimension (nx,ny) :: epsd_xy2
  real, public, dimension (nx,ny) :: pp_xy2
!  Variables for xz slices start here
!  Code variables
  real, public, dimension (nx,nz,3) :: uu_xz,aa_xz,uud_xz,vvp_xz
  real, public, dimension (nx,nz) :: lnrho_xz,ss_xz,cc_xz,lncc_xz
  real, public, dimension (nx,nz) :: XX_chiral_xz,YY_chiral_xz
  real, public, dimension (nx,nz) :: nd_xz
  real, public, dimension (nx,nz,ndustspec) :: md_xz
!  Auxiliary variables  
  real, public, dimension (nx,nz) :: yH_xz,shock_xz,Qrad_xz,ecr_xz
  real, public, dimension (nx,nz,3) :: Frad_xz
!  Derived variables
  real, public, dimension (nx,nz,3) :: oo_xz,bb_xz
  real, public, dimension (nx,nz) :: divu_xz,u2_xz,o2_xz,lnTT_xz,b2_xz,jb_xz
  real, public, dimension (nx,nz) :: DQ_chiral_xz,QQ_chiral_xz
  real, public, dimension (nx,nz) :: epsd_xz
  real, public, dimension (nx,nz) :: pp_xz
!  Variables for yz slices start here
!  Code variables
  real, public, dimension (ny,nz,3) :: uu_yz,aa_yz,uud_yz,vvp_yz
  real, public, dimension (ny,nz) :: lnrho_yz,ss_yz,cc_yz,lncc_yz
  real, public, dimension (ny,nz) :: XX_chiral_yz,YY_chiral_yz
  real, public, dimension (ny,nz) :: nd_yz
  real, public, dimension (ny,nz,ndustspec) :: md_yz
!  Auxiliary variables  
  real, public, dimension (ny,nz) :: yH_yz,shock_yz,Qrad_yz,ecr_yz
  real, public, dimension (ny,nz,3) :: Frad_yz
!  Derived variables
  real, public, dimension (ny,nz,3) :: oo_yz,bb_yz
  real, public, dimension (ny,nz) :: divu_yz,u2_yz,o2_yz,lnTT_yz,b2_yz,jb_yz
  real, public, dimension (ny,nz) :: DQ_chiral_yz,QQ_chiral_yz
  real, public, dimension (ny,nz) :: epsd_yz
  real, public, dimension (ny,nz) :: pp_yz
!
  real, public :: tvid
  integer, public :: nvid

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
      integer, save :: ifirst=0
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
!  16-aug-04/anders: put slices in the order code-aux-derived
!
      use EquationOfState, only: eoscalc, ilnrho_ss
      use General
      use Messages
      use Particles_main, only: get_slices_particles
      use Interstellar, only: get_slices_interstellar
      use Shock, only: get_slices_shock
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (slice_data) :: slices
      character(len=*) :: path
      character(len=4) :: sindex
      logical, save :: lfirstloop=.true.
      logical :: lnewfile=.true.
      logical :: lslices_legacy=.true.
      integer :: inamev,k,i
      real :: tmpval
      integer :: l

      slices%ix=ix_loc
      slices%iy=iy_loc
      slices%iz=iz_loc
      slices%iz2=iz2_loc
      slices%ready=.false.
      slices%index=0
!
!  Loop over slices
!
      inamev = 1
      do while (inamev <= nnamev)
        slices%xy=>slice_xy
        slices%xz=>slice_xz
        slices%yz=>slice_yz
        slices%xy2=>slice_xy2

        slices%name=trim(cnamev(inamev))
        lslices_legacy=.true.       ! By default assume we're not 
                                    ! using module hooks to get the
                                    ! slice contents
        select case (trim(slices%name))
!
!  Velocity field (code variable)
!
        case ('uu')
          uu_yz=f(ix_loc,m1:m2,n1:n2,iux:iuz)
          uu_xz=f(l1:l2,iy_loc,n1:n2,iux:iuz)
          uu_xy=f(l1:l2,m1:m2,iz_loc,iux:iuz)
          uu_xy2=f(l1:l2,m1:m2,iz2_loc,iux:iuz)
          call wslice(path//'ux.yz',uu_yz(:,:,1),x(ix_loc),ny,nz)
          call wslice(path//'uy.yz',uu_yz(:,:,2),x(ix_loc),ny,nz)
          call wslice(path//'uz.yz',uu_yz(:,:,3),x(ix_loc),ny,nz)
          call wslice(path//'ux.xz',uu_xz(:,:,1),y(iy_loc),nx,nz)
          call wslice(path//'uy.xz',uu_xz(:,:,2),y(iy_loc),nx,nz)
          call wslice(path//'uz.xz',uu_xz(:,:,3),y(iy_loc),nx,nz)
          call wslice(path//'ux.xy',uu_xy(:,:,1),z(iz_loc),nx,ny)
          call wslice(path//'uy.xy',uu_xy(:,:,2),z(iz_loc),nx,ny)
          call wslice(path//'uz.xy',uu_xy(:,:,3),z(iz_loc),nx,ny)
          call wslice(path//'ux.Xy',uu_xy2(:,:,1),z(iz2_loc),nx,ny)
          call wslice(path//'uy.Xy',uu_xy2(:,:,2),z(iz2_loc),nx,ny)
          call wslice(path//'uz.Xy',uu_xy2(:,:,3),z(iz2_loc),nx,ny)
!
!  Logarithmic density (code variable)
!
        case ('lnrho')
          lnrho_yz=f(ix_loc,m1:m2,n1:n2,ilnrho)
          lnrho_xz=f(l1:l2,iy_loc,n1:n2,ilnrho)
          lnrho_xy=f(l1:l2,m1:m2,iz_loc,ilnrho)
          lnrho_xy2=f(l1:l2,m1:m2,iz2_loc,ilnrho)
          call wslice(path//'lnrho.yz',lnrho_yz,x(ix_loc),ny,nz)
          call wslice(path//'lnrho.xz',lnrho_xz,y(iy_loc),nx,nz)
          call wslice(path//'lnrho.xy',lnrho_xy,z(iz_loc),nx,ny)
          call wslice(path//'lnrho.Xy',lnrho_xy2,z(iz2_loc),nx,ny)
! 
!  Entropy (code variable)
!
        case ('ss')
          ss_yz=f(ix_loc,m1:m2,n1:n2,iss)
          ss_xz=f(l1:l2,iy_loc,n1:n2,iss)
          ss_xy=f(l1:l2,m1:m2,iz_loc,iss)
          ss_xy2=f(l1:l2,m1:m2,iz2_loc,iss)
          call wslice(path//'ss.yz',ss_yz,x(ix_loc),ny,nz)
          call wslice(path//'ss.xz',ss_xz,y(iy_loc),nx,nz)
          call wslice(path//'ss.xy',ss_xy,z(iz_loc),nx,ny)
          call wslice(path//'ss.Xy',ss_xy2,z(iz2_loc),nx,ny)
!
!  Magnetic vector potential (code variable)
!
        case ('aa')
          aa_yz=f(ix_loc,m1:m2,n1:n2,iax:iaz)
          aa_xz=f(l1:l2,iy_loc,n1:n2,iax:iaz)
          aa_xy=f(l1:l2,m1:m2,iz_loc,iax:iaz)
          aa_xy2=f(l1:l2,m1:m2,iz2_loc,iax:iaz)
          call wslice(path//'ax.yz',aa_yz(:,:,1),x(ix_loc),ny,nz)
          call wslice(path//'ay.yz',aa_yz(:,:,2),x(ix_loc),ny,nz)
          call wslice(path//'az.yz',aa_yz(:,:,3),x(ix_loc),ny,nz)
          call wslice(path//'ax.xz',aa_xz(:,:,1),y(iy_loc),nx,nz)
          call wslice(path//'ay.xz',aa_xz(:,:,2),y(iy_loc),nx,nz)
          call wslice(path//'az.xz',aa_xz(:,:,3),y(iy_loc),nx,nz)
          call wslice(path//'ax.xy',aa_xy(:,:,1),z(iz_loc),nx,ny)
          call wslice(path//'ay.xy',aa_xy(:,:,2),z(iz_loc),nx,ny)
          call wslice(path//'az.xy',aa_xy(:,:,3),z(iz_loc),nx,ny)
          call wslice(path//'ax.Xy',aa_xy2(:,:,1),z(iz2_loc),nx,ny)
          call wslice(path//'ay.Xy',aa_xy2(:,:,2),z(iz2_loc),nx,ny)
          call wslice(path//'az.Xy',aa_xy2(:,:,3),z(iz2_loc),nx,ny)
!
!  Passive scalar (code variable)
!
        case ('cc')
          if (icc==0) then
            if (lroot) print*,'slices: cannot write cc slice; icc=0'
          else
            cc_yz=f(ix_loc,m1:m2,n1:n2,icc)
            cc_xz=f(l1:l2,iy_loc,n1:n2,icc)
            cc_xy=f(l1:l2,m1:m2,iz_loc,icc)
            cc_xy2=f(l1:l2,m1:m2,iz2_loc,icc)
            call wslice(path//'cc.yz',cc_yz,x(ix_loc),ny,nz)
            call wslice(path//'cc.xz',cc_xz,y(iy_loc),nx,nz)
            call wslice(path//'cc.xy',cc_xy,z(iz_loc),nx,ny)
            call wslice(path//'cc.Xy',cc_xy2,z(iz2_loc),nx,ny)
          endif
!
!  Passive scalar (code variable)
!
        case ('lncc')
          if (ilncc==0) then
            if (lroot) print*,'slices: cannot write lncc slice; ilncc=0'
          else
            lncc_yz=f(ix_loc,m1:m2,n1:n2,ilncc)
            lncc_xz=f(l1:l2,iy_loc,n1:n2,ilncc)
            lncc_xy=f(l1:l2,m1:m2,iz_loc,ilncc)
            lncc_xy2=f(l1:l2,m1:m2,iz2_loc,ilncc)
            call wslice(path//'lncc.yz',lncc_yz,x(ix_loc),ny,nz)
            call wslice(path//'lncc.xz',lncc_xz,y(iy_loc),nx,nz)
            call wslice(path//'lncc.xy',lncc_xy,z(iz_loc),nx,ny)
            call wslice(path//'lncc.Xy',lncc_xy2,z(iz2_loc),nx,ny)
          endif
!
!  Chirality fields: XX (code variable)
!
        case ('XX_chiral')
          XX_chiral_yz=f(ix_loc,m1:m2,n1:n2,iXX_chiral)
          XX_chiral_xz=f(l1:l2,iy_loc,n1:n2,iXX_chiral)
          XX_chiral_xy=f(l1:l2,m1:m2,iz_loc,iXX_chiral)
          XX_chiral_xy2=f(l1:l2,m1:m2,iz2_loc,iXX_chiral)
          call wslice(path//'XX_chiral.yz',XX_chiral_yz,x(ix_loc),ny,nz)
          call wslice(path//'XX_chiral.xz',XX_chiral_xz,y(iy_loc),nx,nz)
          call wslice(path//'XX_chiral.xy',XX_chiral_xy,z(iz_loc),nx,ny)
          call wslice(path//'XX_chiral.Xy',XX_chiral_xy2,z(iz2_loc),nx,ny)
!
!  Chirality fields: YY (code variable)
!
        case ('YY_chiral')
          YY_chiral_yz=f(ix_loc,m1:m2,n1:n2,iYY_chiral)
          YY_chiral_xz=f(l1:l2,iy_loc,n1:n2,iYY_chiral)
          YY_chiral_xy=f(l1:l2,m1:m2,iz_loc,iYY_chiral)
          YY_chiral_xy2=f(l1:l2,m1:m2,iz2_loc,iYY_chiral)
          call wslice(path//'YY_chiral.yz',YY_chiral_yz,x(ix_loc),ny,nz)
          call wslice(path//'YY_chiral.xz',YY_chiral_xz,y(iy_loc),nx,nz)
          call wslice(path//'YY_chiral.xy',YY_chiral_xy,z(iz_loc),nx,ny)
          call wslice(path//'YY_chiral.Xy',YY_chiral_xy2,z(iz2_loc),nx,ny)
!
!  Dust velocity (code variable)
!
        case ('uud')
          do k=1,ndustspec
            call chn(k,sindex)
            if (k == 1) sindex = ''
            uud_yz=f(ix_loc,m1:m2,n1:n2,iudx(k):iudz(k))
            uud_xz=f(l1:l2,iy_loc,n1:n2,iudx(k):iudz(k))
            uud_xy=f(l1:l2,m1:m2,iz_loc,iudx(k):iudz(k))
            uud_xy2=f(l1:l2,m1:m2,iz2_loc,iudx(k):iudz(k))
            call wslice(path//'udx'//trim(sindex)//'.yz', &
                uud_yz(:,:,1),x(ix_loc),ny,nz)
            call wslice(path//'udy'//trim(sindex)//'.yz', &
                uud_yz(:,:,2),x(ix_loc),ny,nz)
            call wslice(path//'udz'//trim(sindex)//'.yz', &
                uud_yz(:,:,3),x(ix_loc),ny,nz)
            call wslice(path//'udx'//trim(sindex)//'.xz', &
                uud_xz(:,:,1),y(iy_loc),nx,nz)
            call wslice(path//'udy'//trim(sindex)//'.xz', &
                uud_xz(:,:,2),y(iy_loc),nx,nz)
            call wslice(path//'udz'//trim(sindex)//'.xz', &
                uud_xz(:,:,3),y(iy_loc),nx,nz)
            call wslice(path//'udx'//trim(sindex)//'.xy', &
                uud_xy(:,:,1),z(iz_loc),nx,ny)
            call wslice(path//'udy'//trim(sindex)//'.xy', &
                uud_xy(:,:,2),z(iz_loc),nx,ny)
            call wslice(path//'udz'//trim(sindex)//'.xy', &
                uud_xy(:,:,3),z(iz_loc),nx,ny)
            call wslice(path//'udx'//trim(sindex)//'.Xy', &
                uud_xy2(:,:,1),z(iz2_loc),nx,ny)
            call wslice(path//'udy'//trim(sindex)//'.Xy', &
                uud_xy2(:,:,2),z(iz2_loc),nx,ny)
            call wslice(path//'udz'//trim(sindex)//'.Xy', &
                uud_xy2(:,:,3),z(iz2_loc),nx,ny)
          enddo
!
!  Dust density (code variable)
!
        case ('nd')
          do k=1,ndustspec
            call chn(k,sindex)
            if (k == 1) sindex = ''
            if (ldustdensity) then
              nd_yz=f(ix_loc,m1:m2,n1:n2,ind(k))
              nd_xz=f(l1:l2,iy_loc,n1:n2,ind(k))
              nd_xy=f(l1:l2,m1:m2,iz_loc,ind(k))
              nd_xy2=f(l1:l2,m1:m2,iz2_loc,ind(k))
              call wslice(path//'nd'//trim(sindex)//'.yz',nd_yz,x(ix_loc),ny,nz)
              call wslice(path//'nd'//trim(sindex)//'.xz',nd_xz,y(iy_loc),nx,nz)
              call wslice(path//'nd'//trim(sindex)//'.xy',nd_xy,z(iz_loc),nx,ny)
              call wslice(path//'nd'//trim(sindex)//'.Xy',nd_xy2,z(iz2_loc),nx,ny)
            else
              if (lroot) call warning('WVID', &
                  "Can't use 'nd' slices with nodustdensity")
            endif
          enddo
!
!  Degree of ionization (auxiliary variable)
!
        case ('yH')
          yH_yz=f(ix_loc,m1:m2,n1:n2,iyH)
          yH_xz=f(l1:l2,iy_loc,n1:n2,iyH)
          yH_xy=f(l1:l2,m1:m2,iz_loc,iyH)
          yH_xy2=f(l1:l2,m1:m2,iz2_loc,iyH)
          call wslice(path//'yH.yz',yH_yz,x(ix_loc),ny,nz)
          call wslice(path//'yH.xz',yH_xz,y(iy_loc),nx,nz)
          call wslice(path//'yH.xy',yH_xy,z(iz_loc),nx,ny)
          call wslice(path//'yH.Xy',yH_xy2,z(iz2_loc),nx,ny)
!
!  Shock viscosity (auxiliary variable)
!
        case ('shock')
          shock_yz=f(ix_loc,m1:m2,n1:n2,ishock)
          shock_xz=f(l1:l2,iy_loc,n1:n2,ishock)
          shock_xy=f(l1:l2,m1:m2,iz_loc,ishock)
          shock_xy2=f(l1:l2,m1:m2,iz2_loc,ishock)
          call wslice(path//'shock.yz',shock_yz,x(ix_loc),ny,nz)
          call wslice(path//'shock.xz',shock_xz,y(iy_loc),nx,nz)
          call wslice(path//'shock.xy',shock_xy,z(iz_loc),nx,ny)
          call wslice(path//'shock.Xy',shock_xy2,z(iz2_loc),nx,ny)
!
!  Heating rate (auxiliary variable)
!
        case ('Qrad')
          Qrad_yz=f(ix_loc,m1:m2,n1:n2,iQrad)
          Qrad_xz=f(l1:l2,iy_loc,n1:n2,iQrad)
          Qrad_xy=f(l1:l2,m1:m2,iz_loc,iQrad)
          Qrad_xy2=f(l1:l2,m1:m2,iz2_loc,iQrad)
          call wslice(path//'Qrad.yz',Qrad_yz,x(ix_loc),ny,nz)
          call wslice(path//'Qrad.xz',Qrad_xz,y(iy_loc),nx,nz)
          call wslice(path//'Qrad.xy',Qrad_xy,z(iz_loc),nx,ny)
          call wslice(path//'Qrad.Xy',Qrad_xy2,z(iz2_loc),nx,ny)
!
!  Radiative Flux (auxiliary variable)
!
        case ('Frad')
          Frad_yz=f(ix_loc,m1:m2,n1:n2,iFradx:iFradz)
          Frad_xz=f(l1:l2,iy_loc,n1:n2,iFradx:iFradz)
          Frad_xy=f(l1:l2,m1:m2,iz_loc,iFradx:iFradz)
          Frad_xy2=f(l1:l2,m1:m2,iz2_loc,iFradx:iFradz)
          call wslice(path//'Fradx.yz',Frad_yz(:,:,1),x(ix_loc),ny,nz)
          call wslice(path//'Frady.yz',Frad_yz(:,:,2),x(ix_loc),ny,nz)
          call wslice(path//'Fradz.yz',Frad_yz(:,:,3),x(ix_loc),ny,nz)
          call wslice(path//'Fradx.xz',Frad_xz(:,:,1),y(iy_loc),nx,nz)
          call wslice(path//'Frady.xz',Frad_xz(:,:,2),y(iy_loc),nx,nz)
          call wslice(path//'Fradz.xz',Frad_xz(:,:,3),y(iy_loc),nx,nz)
          call wslice(path//'Fradx.xy',Frad_xy(:,:,1),z(iz_loc),nx,ny)
          call wslice(path//'Frady.xy',Frad_xy(:,:,2),z(iz_loc),nx,ny)
          call wslice(path//'Fradz.xy',Frad_xy(:,:,3),z(iz_loc),nx,ny)
          call wslice(path//'Fradx.Xy',Frad_xy2(:,:,1),z(iz2_loc),nx,ny)
          call wslice(path//'Frady.Xy',Frad_xy2(:,:,2),z(iz2_loc),nx,ny)
          call wslice(path//'Fradz.Xy',Frad_xy2(:,:,3),z(iz2_loc),nx,ny)
!
!  Cosmic ray energy density (auxiliary variable)
!
        case ('ecr')
          ecr_yz=f(ix_loc,m1:m2,n1:n2,iecr)
          ecr_xz=f(l1:l2,iy_loc,n1:n2,iecr)
          ecr_xy=f(l1:l2,m1:m2,iz_loc,iecr)
          ecr_xy2=f(l1:l2,m1:m2,iz2_loc,iecr)
          call wslice(path//'ecr.yz',ecr_yz,x(ix_loc),ny,nz)
          call wslice(path//'ecr.xz',ecr_xz,y(iy_loc),nx,nz)
          call wslice(path//'ecr.xy',ecr_xy,z(iz_loc),nx,ny)
          call wslice(path//'ecr.Xy',ecr_xy2,z(iz2_loc),nx,ny)
!
!  Divergence of velocity (derived variable)
!
        case ('divu')
          call wslice(path//'divu.yz',divu_yz,x(ix_loc),ny,nz)
          call wslice(path//'divu.xz',divu_xz,y(iy_loc),nx,nz)
          call wslice(path//'divu.xy',divu_xy,z(iz_loc),nx,ny)
          call wslice(path//'divu.Xy',divu_xy2,z(iz2_loc),nx,ny)
!
!  Velocity squared (derived variable)
!
        case ('u2')
          call wslice(path//'u2.yz',u2_yz,x(ix_loc),ny,nz)
          call wslice(path//'u2.xz',u2_xz,y(iy_loc),nx,nz)
          call wslice(path//'u2.xy',u2_xy,z(iz_loc),nx,ny)
          call wslice(path//'u2.Xy',u2_xy2,z(iz2_loc),nx,ny)
!
!  Vorticity (derived variable)
!
        case ('oo')
          call wslice(path//'ox.yz',oo_yz(:,:,1),x(ix_loc),ny,nz)
          call wslice(path//'oy.yz',oo_yz(:,:,2),x(ix_loc),ny,nz)
          call wslice(path//'oz.yz',oo_yz(:,:,3),x(ix_loc),ny,nz)
          call wslice(path//'ox.xz',oo_xz(:,:,1),y(iy_loc),nx,nz)
          call wslice(path//'oy.xz',oo_xz(:,:,2),y(iy_loc),nx,nz)
          call wslice(path//'oz.xz',oo_xz(:,:,3),y(iy_loc),nx,nz)
          call wslice(path//'ox.xy',oo_xy(:,:,1),z(iz_loc),nx,ny)
          call wslice(path//'oy.xy',oo_xy(:,:,2),z(iz_loc),nx,ny)
          call wslice(path//'oz.xy',oo_xy(:,:,3),z(iz_loc),nx,ny)
          call wslice(path//'ox.Xy',oo_xy2(:,:,1),z(iz2_loc),nx,ny)
          call wslice(path//'oy.Xy',oo_xy2(:,:,2),z(iz2_loc),nx,ny)
          call wslice(path//'oz.Xy',oo_xy2(:,:,3),z(iz2_loc),nx,ny)
!
!  Vorticity squared (derived variable)
!
        case ('o2')
          call wslice(path//'o2.yz',o2_yz,x(ix_loc),ny,nz)
          call wslice(path//'o2.xz',o2_xz,y(iy_loc),nx,nz)
          call wslice(path//'o2.xy',o2_xy,z(iz_loc),nx,ny)
          call wslice(path//'o2.Xy',o2_xy2,z(iz2_loc),nx,ny)
!
!  Temperature (derived variable, sometimes code variable)
!
        case ('lnTT')
          if (ilnTT .ne. 0) then        
            lnTT_yz=f(ix_loc,m1:m2,n1:n2,ilnTT)
            lnTT_xz=f(l1:l2,iy_loc,n1:n2,ilnTT)
            lnTT_xy=f(l1:l2,m1:m2,iz_loc,ilnTT)
            lnTT_xy2=f(l1:l2,m1:m2,iz2_loc,ilnTT)
          else
            do m=m1,m2; do n=n1,n2
              call eoscalc(ilnrho_ss,f(ix_loc,m,n,ilnrho),f(ix_loc,m,n,iss),lnTT=tmpval)
              lnTT_yz(m-m1+1,n-n1+1)=tmpval
            enddo; enddo
            do l=l1,l2; do n=n1,n2
              call eoscalc(ilnrho_ss,f(l,iy_loc,n,ilnrho),f(l,iy_loc,n,iss),lnTT=tmpval)
              lnTT_xz(l-l1+1,n-n1+1)=tmpval
            enddo; enddo
            do l=l1,l2; do m=m1,m2
              call eoscalc(ilnrho_ss,f(l,m,iz_loc,ilnrho),f(l,m,iz_loc,iss),lnTT=tmpval)
              lnTT_xy(l-l1+1,m-m1+1)=tmpval
              call eoscalc(ilnrho_ss,f(l,m,iz2_loc,ilnrho),f(l,m,iz2_loc,iss), &
                  lnTT=tmpval)
              lnTT_xy2(l-l1+1,m-m1+1)=tmpval
            enddo; enddo
          endif
          call wslice(path//'lnTT.yz',lnTT_yz,x(ix_loc),ny,nz)
          call wslice(path//'lnTT.xz',lnTT_xz,y(iy_loc),nx,nz)
          call wslice(path//'lnTT.xy',lnTT_xy,z(iz_loc),nx,ny)
          call wslice(path//'lnTT.Xy',lnTT_xy2,z(iz2_loc),nx,ny)
!
!  Pressure (derived variable)
!
        case ('pp')
          do m=m1,m2; do n=n1,n2
            call eoscalc(ilnrho_ss,f(ix_loc,m,n,ilnrho),f(ix_loc,m,n,iss),pp=tmpval)
            pp_yz(m-m1+1,n-n1+1)=tmpval
          enddo; enddo
          do l=l1,l2; do n=n1,n2
            call eoscalc(ilnrho_ss,f(l,iy_loc,n,ilnrho),f(l,iy_loc,n,iss),pp=tmpval)
            pp_xz(l-l1+1,n-n1+1)=tmpval
          enddo; enddo
          do l=l1,l2; do m=m1,m2
            call eoscalc(ilnrho_ss,f(l,m,iz_loc,ilnrho),f(l,m,iz_loc,iss),pp=tmpval)
            pp_xy(l-l1+1,m-m1+1)=tmpval
            call eoscalc(ilnrho_ss,f(l,m,iz2_loc,ilnrho),f(l,m,iz2_loc,iss), &
                pp=tmpval)
            pp_xy2(l-l1+1,m-m1+1)=tmpval
          enddo; enddo
          call wslice(path//'pp.yz',pp_yz,x(ix_loc),ny,nz)
          call wslice(path//'pp.xz',pp_xz,y(iy_loc),nx,nz)
          call wslice(path//'pp.xy',pp_xy,z(iz_loc),nx,ny)
          call wslice(path//'pp.Xy',pp_xy2,z(iz2_loc),nx,ny)
!
!  Magnetic field (derived variable)
!
        case ('bb')
          call wslice(path//'bx.yz',bb_yz(:,:,1),x(ix_loc),ny,nz)
          call wslice(path//'by.yz',bb_yz(:,:,2),x(ix_loc),ny,nz)
          call wslice(path//'bz.yz',bb_yz(:,:,3),x(ix_loc),ny,nz)
          call wslice(path//'bx.xz',bb_xz(:,:,1),y(iy_loc),nx,nz)
          call wslice(path//'by.xz',bb_xz(:,:,2),y(iy_loc),nx,nz)
          call wslice(path//'bz.xz',bb_xz(:,:,3),y(iy_loc),nx,nz)
          call wslice(path//'bx.xy',bb_xy(:,:,1),z(iz_loc),nx,ny)
          call wslice(path//'by.xy',bb_xy(:,:,2),z(iz_loc),nx,ny)
          call wslice(path//'bz.xy',bb_xy(:,:,3),z(iz_loc),nx,ny)
          call wslice(path//'bx.Xy',bb_xy2(:,:,1),z(iz2_loc),nx,ny)
          call wslice(path//'by.Xy',bb_xy2(:,:,2),z(iz2_loc),nx,ny)
          call wslice(path//'bz.Xy',bb_xy2(:,:,3),z(iz2_loc),nx,ny)
!
!  Magnetic field squared (derived variable)
!
        case ('b2')
          call wslice(path//'b2.yz',b2_yz,x(ix_loc),ny,nz)
          call wslice(path//'b2.xz',b2_xz,y(iy_loc),nx,nz)
          call wslice(path//'b2.xy',b2_xy,z(iz_loc),nx,ny)
          call wslice(path//'b2.Xy',b2_xy2,z(iz2_loc),nx,ny)
!
!  Current density (derived variable)
!
        case ('jb')
          call wslice(path//'jb.yz',jb_yz,x(ix_loc),ny,nz)
          call wslice(path//'jb.xz',jb_xz,y(iy_loc),nx,nz)
          call wslice(path//'jb.xy',jb_xy,z(iz_loc),nx,ny)
          call wslice(path//'jb.Xy',jb_xy2,z(iz2_loc),nx,ny)
!
!  chirality fields: DQ (derived variable)
!
        case ('DQ_chiral')
          XX_chiral_yz=f(ix_loc,m1:m2,n1:n2,iXX_chiral)
          XX_chiral_xz=f(l1:l2,iy_loc,n1:n2,iXX_chiral)
          XX_chiral_xy=f(l1:l2,m1:m2,iz_loc,iXX_chiral)
          XX_chiral_xy2=f(l1:l2,m1:m2,iz2_loc,iXX_chiral)
          YY_chiral_yz=f(ix_loc,m1:m2,n1:n2,iYY_chiral)
          YY_chiral_xz=f(l1:l2,iy_loc,n1:n2,iYY_chiral)
          YY_chiral_xy=f(l1:l2,m1:m2,iz_loc,iYY_chiral)
          YY_chiral_xy2=f(l1:l2,m1:m2,iz2_loc,iYY_chiral)
          QQ_chiral_yz=XX_chiral_yz-YY_chiral_yz
          QQ_chiral_xz=XX_chiral_xz-YY_chiral_xz
          QQ_chiral_xy=XX_chiral_xy-YY_chiral_xy
          QQ_chiral_xy2=XX_chiral_xy2-YY_chiral_xy2
          DQ_chiral_yz=QQ_chiral_yz*(1.-QQ_chiral_yz**2)/(1.+QQ_chiral_yz**2)
          DQ_chiral_xz=QQ_chiral_xz*(1.-QQ_chiral_xz**2)/(1.+QQ_chiral_xz**2)
          DQ_chiral_xy=QQ_chiral_xy*(1.-QQ_chiral_xy**2)/(1.+QQ_chiral_xy**2)
          DQ_chiral_xy2=&
              QQ_chiral_xy2*(1.-QQ_chiral_xy2**2)/(1.+QQ_chiral_xy2**2)
          call wslice(path//'DQ_chiral.yz',DQ_chiral_yz,x(ix_loc),ny,nz)
          call wslice(path//'DQ_chiral.xz',DQ_chiral_xz,y(iy_loc),nx,nz)
          call wslice(path//'DQ_chiral.xy',DQ_chiral_xy,z(iz_loc),nx,ny)
          call wslice(path//'DQ_chiral.Xy',DQ_chiral_xy2,z(iz2_loc),nx,ny)
!
!  psi2 - Absolute value of the wave function squared
!
        case ('psi2')
          lnrho_yz=f(ix_loc,m1:m2,n1:n2,1)**2 + f(ix_loc,m1:m2,n1:n2,2)**2 
          lnrho_xz=f(l1:l2,iy_loc,n1:n2,1)**2 + f(l1:l2,iy_loc,n1:n2,2)**2 
          lnrho_xy=f(l1:l2,m1:m2,iz_loc,1)**2 + f(l1:l2,m1:m2,iz_loc,2)**2
          lnrho_xy2=f(l1:l2,m1:m2,iz2_loc,1)**2 + f(l1:l2,m1:m2,iz2_loc,2)**2 
          call wslice(path//'psi2.yz',lnrho_yz,x(ix_loc),ny,nz)
          call wslice(path//'psi2.xz',lnrho_xz,y(iy_loc),nx,nz)
          call wslice(path//'psi2.xy',lnrho_xy,z(iz_loc),nx,ny)
          call wslice(path//'psi2.Xy',lnrho_xy2,z(iz2_loc),nx,ny)
!
!  Dust-to-gas mass ratio (derived variable)
!
        case ('epsd')
          do k=1,ndustspec
            call chn(k,sindex)
            if (k == 1) sindex = ''
            if (ldustdensity_log) then
              nd_yz=exp(f(ix_loc,m1:m2,n1:n2,ind(k)))
              nd_xz=exp(f(l1:l2,iy_loc,n1:n2,ind(k)))
              nd_xy=exp(f(l1:l2,m1:m2,iz_loc,ind(k)))
              nd_xy2=exp(f(l1:l2,m1:m2,iz2_loc,ind(k)))
            else
              nd_yz=f(ix_loc,m1:m2,n1:n2,ind(k))
              nd_xz=f(l1:l2,iy_loc,n1:n2,ind(k))
              nd_xy=f(l1:l2,m1:m2,iz_loc,ind(k))
              nd_xy2=f(l1:l2,m1:m2,iz2_loc,ind(k))
            endif
            epsd_yz=md_yz(:,:,k)*nd_yz/exp(f(ix_loc,m1:m2,n1:n2,ilnrho))
            epsd_xz=md_xz(:,:,k)*nd_xz/exp(f(l1:l2,iy_loc,n1:n2,ilnrho))
            epsd_xy=md_xy(:,:,k)*nd_xy/exp(f(l1:l2,m1:m2,iz_loc,ilnrho))
            epsd_xy2=md_xy2(:,:,k)*nd_xy2/exp(f(l1:l2,m1:m2,iz2_loc,ilnrho))
            call wslice(path//'epsd'//trim(sindex)//'.yz',epsd_yz,x(ix_loc),ny,nz)
            call wslice(path//'epsd'//trim(sindex)//'.xz',epsd_xz,y(iy_loc),nx,nz)
            call wslice(path//'epsd'//trim(sindex)//'.xy',epsd_xy,z(iz_loc),nx,ny)
            call wslice(path//'epsd'//trim(sindex)//'.Xy',epsd_xy2,z(iz2_loc),nx,ny)
          enddo
!
!  Surface intensity (derived variable)
!
        case ('Isurf')
          call wslice(path//'Isurf.xy',Isurf_xy,z(iz2_loc),nx,ny)
!
        case default
!
! Add new slice providing routines here!
!
          lslices_legacy=.false.
          if (lparticles)    call get_slices_particles   (f,slices)
          if (lshock)        call get_slices_shock       (f,slices)
          if (linterstellar) call get_slices_interstellar(f,slices)
!
        endselect

        if (lslices_legacy) then
          inamev=inamev+1
          cycle
        endif

        if (slices%ready) then
          if (slices%index==0) then    ! If this wasn't a multi slice...
            if (associated(slices%yz)) &
              call wslice(path//trim(slices%name)//'.yz',slices%yz, &
                                                     x(slices%ix),ny,nz)
            if (associated(slices%xz)) &
              call wslice(path//trim(slices%name)//'.xz',slices%xz, &
                                                     y(slices%iy),nx,nz)
            if (associated(slices%xy)) &
              call wslice(path//trim(slices%name)//'.xy',slices%xy, &
                                                     z(slices%iz),nx,ny)
            if (associated(slices%xy2)) &
              call wslice(path//trim(slices%name)//'.Xy',slices%xy2, &
                                                     z(slices%iz2),nx,ny)
            inamev=inamev+1
           else 
            call chn(slices%index, sindex)
            if (associated(slices%yz)) &
              call wslice(path//trim(slices%name)//trim(sindex)//'.yz', &
                                       slices%yz, x(slices%ix),ny,nz)
            if (associated(slices%xz)) &
              call wslice(path//trim(slices%name)//trim(sindex)//'.xz', &
                                       slices%xz, y(slices%iy),nx,nz)
            if (associated(slices%xy)) &
              call wslice(path//trim(slices%name)//trim(sindex)//'.xy', &
                                       slices%xy, z(slices%iz),nx,ny)
            if (associated(slices%xy2)) &
              call wslice(path//trim(slices%name)//trim(sindex)//'.Xy', &
                                       slices%xy2, z(slices%iz2),nx,ny)
          endif
        else
          if (slices%index==0) then    ! If this wasn't a multi slice...
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
          else
            slices%index=0
            inamev=inamev+1
          endif
        endif
      enddo
!
!
      lfirstloop=.false.
!
    endsubroutine wvid
!***********************************************************************
    subroutine wslice(filename,a,pos,ndim1,ndim2)
!
!  appending to an existing slice file
!
!  12-nov-02/axel: coded
!  26-jun-06/anders: moved from Slices
!
      use Cdata, only: t, lwrite_slice_xy2, lwrite_slice_xy, lwrite_slice_xz, lwrite_slice_yz
!
      integer :: ndim1,ndim2
      character (len=*) :: filename
      real, dimension (ndim1,ndim2) :: a
      real, intent(in) :: pos
!
!  check whether we want to write a slice on this processor
!
      if ( (lwrite_slice_xy2.and.index(filename,'Xy')>0) .or. &
           (lwrite_slice_xy .and.index(filename,'xy')>0) .or. &
           (lwrite_slice_xz .and.index(filename,'xz')>0) .or. &
           (lwrite_slice_yz .and.index(filename,'yz')>0) ) then
        open(1,file=filename,form='unformatted',position='append')
        write(1) a,t,pos
        close(1)
      endif
!
    endsubroutine wslice
!***********************************************************************
    subroutine setup_slices()
!
!  Determine slice positions and whether slices are to be written on this
!  processor.
!
!  29-may-06/tobi: wrapped code from param_io.f90 into this subroutine
!
      use Cdata
      use Messages, only: fatal_error

!
!  set slice position. The default for slice_position is 'p' for periphery,
!  although setting ix, iy, iz, iz2 by hand will overwrite this.
!  If slice_position is not 'p', then ix, iy, iz, iz2 are overwritten.
!
      if (slice_position=='p' .or. slice_position=='S') then
        ix_loc=l1; iy_loc=m1; iz_loc=n1; iz2_loc=n2
        lwrite_slice_xy2=(ipz==nprocz-1)
        lwrite_slice_xy=(ipz==0)
        lwrite_slice_xz=(ipy==0)
        lwrite_slice_yz=.true.
!
!  slice position when the first meshpoint in z is the equator (sphere)
!  For one z-processor, iz remains n1, but iz2 is set to the middle.
!
!  TH: The text above does not properly describe the code below.
!  TH: A certain processor layout is implied here
!
      elseif (slice_position=='m') then
        ix_loc=(l1+l2)/2; iy_loc=m1; iz_loc=n1; iz2_loc=n2
        if(nprocy==1) then; iy_loc=(m1+m2)/2; endif
        if(nprocz==1) then; iz_loc=(n1+n2)/2; iz2_loc=(iz+n2)/2; endif
        lwrite_slice_xy2=(ipz==nprocz/2)
        lwrite_slice_xy=(ipz==nprocz/2)
        lwrite_slice_xz=(ipy==nprocy/2)
        lwrite_slice_yz=.true.
!
!  slice position when the first meshpoint in z is the equator (sphere)
!  For one z-processor, iz remains n1, but iz2 is set to the middle.
!
!  TH: A certain processor layout is implied here
!
      elseif (slice_position=='e') then
        ix_loc=(l1+l2)/2; iy_loc=m1; iz_loc=n1; iz2_loc=n2
        if(nprocy==1) then; iy_loc=(m1+m2)/2; endif
        if(nprocz==1) then; iz2_loc=(iz+n2)/2; endif
        lwrite_slice_xy2=(ipz==nprocz/4)
        lwrite_slice_xy=(ipz==0)
        lwrite_slice_xz=(ipy==nprocy/2)
        lwrite_slice_yz=.true.
!
!  slice position when the first meshpoint in z is the equator (sphere)
!  For one z-processor, iz remains n1, but iz2 is set to the middle.
!
!  TH: The text above does not properly describe the code below.
!  TH: A certain processor layout is implied here
!
      elseif (slice_position=='c') then
        ix_loc=(l1+l2)/2; iy_loc=m1; iz_loc=n1; iz2_loc=n2
        lwrite_slice_xy2=(ipz==nprocz-1)
        lwrite_slice_xy=(ipz==0)
        lwrite_slice_xz=(ipy==0)
        lwrite_slice_yz=.true.
!
!  periphery of the box, but the other way around
!
      elseif (slice_position=='q') then
        ix_loc=l2; iy_loc=m2; iz_loc=n2; iz2_loc=n1
        lwrite_slice_xy2=(ipz==0)
        lwrite_slice_xy=(ipz==nprocz-1)
        lwrite_slice_xz=(ipy==nprocy-1)
        lwrite_slice_yz=.true.
      else
        if (lroot) then
          call fatal_error('setup_slices', &
                           'No such value for slice_position: '//slice_position)
        endif
      endif
!
!  Overwrite slice postions if any ix,iy,iz,iz2 is greater then Zero
!
      if (ix>0) then
        ix_loc=ix-ipx*nx
        if (ix_loc>=l1.and.ix_loc<=l2) then
          lwrite_slice_yz=.true.
        else
          lwrite_slice_yz=.false.
        endif
      endif

      if (iy>0) then
        iy_loc=iy-ipy*ny
        if (iy_loc>=m1.and.iy_loc<=m2) then
          lwrite_slice_xz=.true.
        else
          lwrite_slice_xz=.false.
        endif
      endif

      if (iz>0) then
        iz_loc=iz-ipz*nz
        if (iz_loc>=n1.and.iz_loc<=n2) then
          lwrite_slice_xy=.true.
        else
          lwrite_slice_xy=.false.
        endif
      endif

      if (iz2>0) then
        iz2_loc=iz2-ipz*nz
        if (iz2_loc>=n1.and.iz2_loc<=n2) then
          lwrite_slice_xy2=.true.
        else
          lwrite_slice_xy2=.false.
        endif
      endif

!
!  write slice position to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/slice_position.dat',STATUS='unknown')
        write(1,'(a)') slice_position
        close(1)
      endif
!  
!  make sure ix_loc,iy_loc,iz_loc,iz2_loc are not outside the boundaries
!
      ix_loc=min(ix_loc,l2); iy_loc=min(iy_loc,m2)
      ix_loc=max(ix_loc,l1); iy_loc=max(iy_loc,m1)
      iz_loc=min(iz_loc,n2); iz2_loc=min(iz2_loc,n2)
      iz_loc=max(iz_loc,n1); iz2_loc=max(iz2_loc,n1)

      if (lroot) then
        write (*,*)'read_runpars: slice_position = '//slice_position
        write (*,'(1x,a,4i4)') &
          'read_runpars: ix,iy,iz,iz2 (video files) =',ix,iy,iz,iz2
      endif

    endsubroutine setup_slices
!***********************************************************************

endmodule Slices
