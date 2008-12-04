! $Id$

!  This module produces slices for animation purposes

module Slices

  use Cdata

  implicit none

  private

  public :: wvid, wvid_prepare, setup_slices, wslice


!  Variables for xy slices start here
!!! New slice code reuses the following slice variables.
  real, target, dimension (nx,ny) :: slice_xy  = 0., slice_xy2 = 0.
  real, target, dimension (nx,nz) :: slice_xz = 0.
  real, target, dimension (ny,nz) :: slice_yz = 0.
!
!  xy3 and xy4 for spherical slices  
!
  real, target, dimension (nx,ny) :: slice_xy3 = 0., slice_xy4 = 0.
!
!!! LEGACY SLICE VARIABLES FOLLOW
!  Code variables
!
! Density Slices
!
  real, public, dimension (nx,ny) :: lnrho_xy
!
! EOS Slices
!
  real, public, dimension (nx,ny) :: ss_xy
  real, public, dimension (nx,ny) :: yH_xy
  real, public, dimension (nx,ny) :: lnTT_xy
  real, public, dimension (nx,ny) :: pp_xy
!
! Pscalar Slices
!
  real, public, dimension (nx,ny) :: cc_xy,lncc_xy
!
! Dust slices
!
  real, public, dimension (nx,ny) :: nd_xy
  real, public, dimension (nx,ny,ndustspec) :: md_xy
  real, public, dimension (nx,ny,3) :: uud_xy,vvp_xy
!
! Chemistry slices
!
  real, public, dimension (nx,ny) :: chemspec_xy
!
!  Cosmic Ray slices
!
  real, public, dimension (nx,ny) :: ecr_xy
!  Derived variables
  real, public, dimension (nx,ny) :: epsd_xy
!  Variables for xy2 slices start here
!  Code variables
  real, public, dimension (nx,ny,3) :: uud_xy2,vvp_xy2
  real, public, dimension (nx,ny) :: lnrho_xy2,ss_xy2,cc_xy2,lncc_xy2
  real, public, dimension (nx,ny) :: nd_xy2
  real, public, dimension (nx,ny,ndustspec) :: md_xy2
  real, public, dimension (nx,ny) :: chemspec_xy2
!  Auxiliary variables
  real, public, dimension (nx,ny) :: yH_xy2,ecr_xy2
!  Derived variables
  real, public, dimension (nx,ny) :: lnTT_xy2
  real, public, dimension (nx,ny) :: epsd_xy2
  real, public, dimension (nx,ny) :: pp_xy2
!  Variables for xy3 and xy4, necessary for 
!  spherical slices, start here
  real, public, dimension (nx,ny) :: lnrho_xy3,lnrho_xy4
  real, public, dimension (nx,ny) :: pp_xy3,pp_xy4
  real, public, dimension (nx,ny) :: ss_xy3,ss_xy4
  real, public, dimension (nx,ny) :: lnTT_xy3,lnTT_xy4
!  Variables for xz slices start here
!  Code variables
  real, public, dimension (nx,nz,3) :: uud_xz,vvp_xz
  real, public, dimension (nx,nz) :: lnrho_xz,ss_xz,cc_xz,lncc_xz
  real, public, dimension (nx,nz) :: nd_xz
  real, public, dimension (nx,nz,ndustspec) :: md_xz
  real, public, dimension (nx,nz) :: chemspec_xz
!  Auxiliary variables
  real, public, dimension (nx,nz) :: yH_xz,ecr_xz
!  Derived variables
  real, public, dimension (nx,nz) :: lnTT_xz
  real, public, dimension (nx,nz) :: epsd_xz
  real, public, dimension (nx,nz) :: pp_xz
!  Variables for yz slices start here
!  Code variables
  real, public, dimension (ny,nz,3) :: uud_yz,vvp_yz
  real, public, dimension (ny,nz) :: lnrho_yz,ss_yz,cc_yz,lncc_yz
  real, public, dimension (ny,nz) :: nd_yz
  real, public, dimension (ny,nz,ndustspec) :: md_yz  
  real, public, dimension (ny,nz) :: chemspec_yz
!  Auxiliary variables
  real, public, dimension (ny,nz) :: yH_yz,ecr_yz
!  Derived variables
  real, public, dimension (ny,nz) :: lnTT_yz
  real, public, dimension (ny,nz) :: epsd_yz
  real, public, dimension (ny,nz) :: pp_yz
!
  real, public :: tvid
  integer, public :: nvid
  real :: tslice
!
  contains
!
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
      character (len=5) :: ch
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
!  Save current time so that the time that is written out in wslice() is not
!  from the next time step
!
      if (lvid) tslice = t
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
!  22-sep-07/axel: changed Xy to xy2, to be compatible with Mac
!
      use EquationOfState, only: eoscalc, ilnrho_ss
      use General
      use Messages
      use Particles_main,  only: get_slices_particles
      use Interstellar,    only: get_slices_interstellar
      use Shock,           only: get_slices_shock
      use Magnetic,        only: get_slices_magnetic
      use Testscalar,      only: get_slices_testscalar
      use Testfield,       only: get_slices_testfield
      use Testflow,        only: get_slices_testflow
      use Hydro,           only: get_slices_hydro
      use Radiation,       only: get_slices_radiation
      use Chiral,          only: get_slices_chiral
      use Special,         only: get_slices_special
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
      character(len=*) :: path
      character(len=5) :: sindex
      logical, save :: lfirstloop=.true.
      logical :: lnewfile=.true.
      logical :: lslices_legacy=.true.
      integer :: inamev,k
      real :: tmpval
      integer :: l

      slices%ix=ix_loc
      slices%iy=iy_loc
      slices%iz=iz_loc
      slices%iz2=iz2_loc
      slices%iz3=iz3_loc
      slices%iz4=iz4_loc
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
        select case (slices%name)
!
!  Logarithmic density (code variable)
!
        case ('lnrho')
          lnrho_yz=f(ix_loc,m1:m2,n1:n2,ilnrho)
          lnrho_xz=f(l1:l2,iy_loc,n1:n2,ilnrho)
          lnrho_xy=f(l1:l2,m1:m2,iz_loc,ilnrho)
          lnrho_xy2=f(l1:l2,m1:m2,iz2_loc,ilnrho)
          lnrho_xy3=f(l1:l2,m1:m2,iz3_loc,ilnrho)
          lnrho_xy4=f(l1:l2,m1:m2,iz4_loc,ilnrho)
          call wslice(path//'lnrho.yz',lnrho_yz,x(ix_loc),ny,nz)
          call wslice(path//'lnrho.xz',lnrho_xz,y(iy_loc),nx,nz)
          call wslice(path//'lnrho.xy',lnrho_xy,z(iz_loc),nx,ny)
          call wslice(path//'lnrho.xy2',lnrho_xy2,z(iz2_loc),nx,ny)
          if (lwrite_slice_xy3) &
              call wslice(path//'lnrho.xy3',lnrho_xy3,z(iz3_loc),nx,ny)
          if (lwrite_slice_xy4) & 
              call wslice(path//'lnrho.xy4',lnrho_xy4,z(iz4_loc),nx,ny)
!
!  Entropy (code variable)
!
        case ('ss')
          ss_yz=f(ix_loc,m1:m2,n1:n2,iss)
          ss_xz=f(l1:l2,iy_loc,n1:n2,iss)
          ss_xy=f(l1:l2,m1:m2,iz_loc,iss)
          ss_xy2=f(l1:l2,m1:m2,iz2_loc,iss)
          ss_xy3=f(l1:l2,m1:m2,iz3_loc,iss)
          ss_xy4=f(l1:l2,m1:m2,iz4_loc,iss)
          call wslice(path//'ss.yz',ss_yz,x(ix_loc),ny,nz)
          call wslice(path//'ss.xz',ss_xz,y(iy_loc),nx,nz)
          call wslice(path//'ss.xy',ss_xy,z(iz_loc),nx,ny)
          call wslice(path//'ss.xy2',ss_xy2,z(iz2_loc),nx,ny)
          if (lwrite_slice_xy3) &
              call wslice(path//'ss.xy3',ss_xy3,z(iz3_loc),nx,ny)
          if (lwrite_slice_xy4) & 
              call wslice(path//'ss.xy4',ss_xy4,z(iz4_loc),nx,ny)
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
            call wslice(path//'cc.xy2',cc_xy2,z(iz2_loc),nx,ny)
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
            call wslice(path//'lncc.xy2',lncc_xy2,z(iz2_loc),nx,ny)
          endif
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
            call wslice(path//'udx'//trim(sindex)//'.xy2', &
                uud_xy2(:,:,1),z(iz2_loc),nx,ny)
            call wslice(path//'udy'//trim(sindex)//'.xy2', &
                uud_xy2(:,:,2),z(iz2_loc),nx,ny)
            call wslice(path//'udz'//trim(sindex)//'.xy2', &
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
              call wslice(path//'nd'//trim(sindex)//'.xy2',nd_xy2,z(iz2_loc),nx,ny)
            else
              if (lroot) call warning('WVID', &
                  "Can't use 'nd' slices with nodustdensity")
            endif
          enddo
!
!  Chemical species mass fractions (code variable)
!
        case ('chemspec')
          do k=1,nchemspec
            call chn(k,sindex)
            if (k == 1) sindex = ''
            if (lchemistry) then
              chemspec_yz=f(ix_loc,m1:m2,n1:n2,ichemspec(k))
              chemspec_xz=f(l1:l2,iy_loc,n1:n2,ichemspec(k))
              chemspec_xy=f(l1:l2,m1:m2,iz_loc,ichemspec(k))
              chemspec_xy2=f(l1:l2,m1:m2,iz2_loc,ichemspec(k))
              call wslice(path//'chemspec'//trim(sindex)//'.yz',chemspec_yz,x(ix_loc),ny,nz)
              call wslice(path//'chemspec'//trim(sindex)//'.xz',chemspec_xz,y(iy_loc),nx,nz)
              call wslice(path//'chemspec'//trim(sindex)//'.xy',chemspec_xy,z(iz_loc),nx,ny)
              call wslice(path//'chemspec'//trim(sindex)//'.xy2',chemspec_xy2,z(iz2_loc),nx,ny)
            else
              if (lroot) call warning('WVID', &
                  "Can't use 'chemspec' slices with nochemistry")
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
          call wslice(path//'yH.xy2',yH_xy2,z(iz2_loc),nx,ny)
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
          call wslice(path//'ecr.xy2',ecr_xy2,z(iz2_loc),nx,ny)
!
!  Temperature (derived variable, sometimes code variable)
!
        case ('lnTT')
          if (ilnTT .ne. 0) then
            lnTT_yz=f(ix_loc,m1:m2,n1:n2,ilnTT)
            lnTT_xz=f(l1:l2,iy_loc,n1:n2,ilnTT)
            lnTT_xy=f(l1:l2,m1:m2,iz_loc,ilnTT)
            lnTT_xy2=f(l1:l2,m1:m2,iz2_loc,ilnTT)
            lnTT_xy3=f(l1:l2,m1:m2,iz3_loc,ilnTT)
            lnTT_xy4=f(l1:l2,m1:m2,iz4_loc,ilnTT)
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
              call eoscalc(ilnrho_ss,f(l,m,iz2_loc,ilnrho),f(l,m,iz2_loc,iss),lnTT=tmpval)
              lnTT_xy2(l-l1+1,m-m1+1)=tmpval
              call eoscalc(ilnrho_ss,f(l,m,iz3_loc,ilnrho),f(l,m,iz3_loc,iss),lnTT=tmpval)
              lnTT_xy3(l-l1+1,m-m1+1)=tmpval
              call eoscalc(ilnrho_ss,f(l,m,iz4_loc,ilnrho),f(l,m,iz4_loc,iss),lnTT=tmpval)
              lnTT_xy4(l-l1+1,m-m1+1)=tmpval
            enddo; enddo
          endif
          call wslice(path//'lnTT.yz',lnTT_yz,x(ix_loc),ny,nz)
          call wslice(path//'lnTT.xz',lnTT_xz,y(iy_loc),nx,nz)
          call wslice(path//'lnTT.xy',lnTT_xy,z(iz_loc),nx,ny)
          call wslice(path//'lnTT.xy2',lnTT_xy2,z(iz2_loc),nx,ny)
          if (lwrite_slice_xy3) &
              call wslice(path//'lnTT.xy3',lnTT_xy3,z(iz3_loc),nx,ny)
          if (lwrite_slice_xy4) & 
              call wslice(path//'lnTT.xy4',lnTT_xy4,z(iz4_loc),nx,ny)
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
            call eoscalc(ilnrho_ss,f(l,m,iz2_loc,ilnrho),f(l,m,iz2_loc,iss),pp=tmpval)
            pp_xy2(l-l1+1,m-m1+1)=tmpval
            call eoscalc(ilnrho_ss,f(l,m,iz3_loc,ilnrho),f(l,m,iz3_loc,iss),pp=tmpval)
            pp_xy3(l-l1+1,m-m1+1)=tmpval
            call eoscalc(ilnrho_ss,f(l,m,iz4_loc,ilnrho),f(l,m,iz4_loc,iss),pp=tmpval)
            pp_xy4(l-l1+1,m-m1+1)=tmpval
          enddo; enddo
          call wslice(path//'pp.yz',pp_yz,x(ix_loc),ny,nz)
          call wslice(path//'pp.xz',pp_xz,y(iy_loc),nx,nz)
          call wslice(path//'pp.xy',pp_xy,z(iz_loc),nx,ny)
          call wslice(path//'pp.xy2',pp_xy2,z(iz2_loc),nx,ny)
          if (lwrite_slice_xy3) &
              call wslice(path//'pp.xy3',pp_xy3,z(iz3_loc),nx,ny)
          if (lwrite_slice_xy4) & 
              call wslice(path//'pp.xy4',pp_xy4,z(iz4_loc),nx,ny)
!
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
            call wslice(path//'epsd'//trim(sindex)//'.xy2',epsd_xy2,z(iz2_loc),nx,ny)
          enddo
!
        case default
!
! Add new slice providing routines here!
!
          lslices_legacy=.false.
          if (lparticles)    call get_slices_particles   (f,slices)
          if (lshock)        call get_slices_shock       (f,slices)
          if (linterstellar) call get_slices_interstellar(f,slices)
          if (lhydro)        call get_slices_hydro       (f,slices)
          if (lmagnetic)     call get_slices_magnetic    (f,slices)
          if (ltestscalar)   call get_slices_testscalar  (f,slices)
          if (ltestfield)    call get_slices_testfield   (f,slices)
          if (ltestflow)     call get_slices_testflow    (f,slices)
          if (lradiation)    call get_slices_radiation   (f,slices)
          if (lchiral)       call get_slices_chiral      (f,slices)
          if (lspecial)      call get_slices_special     (f,slices)
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
              call wslice(path//trim(slices%name)//'.xy2',slices%xy2, &
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
              call wslice(path//trim(slices%name)//trim(sindex)//'.xy2', &
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
          endif
          inamev=inamev+1
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
!  22-sep-07/axel: changed Xy to xy2, to be compatible with Mac
!
      use Cdata, only: t, lwrite_slice_xy ,lwrite_slice_xy2, &
                          lwrite_slice_xy3,lwrite_slice_xy4, &
                          lwrite_slice_xz ,lwrite_slice_yz 
!
      integer :: ndim1,ndim2
      character (len=*) :: filename
      real, dimension (ndim1,ndim2) :: a
      real, intent(in) :: pos
!
!  check whether we want to write a slice on this processor
!
      if ( (lwrite_slice_xy .and.index(filename,'xy' )>0) .or. & 
           (lwrite_slice_xy2.and.index(filename,'xy2')>0) .or. &
           (lwrite_slice_xy3.and.index(filename,'xy3')>0) .or. &
           (lwrite_slice_xy4.and.index(filename,'xy4')>0) .or. &
           (lwrite_slice_xz .and.index(filename,'xz' )>0) .or. &
           (lwrite_slice_yz .and.index(filename,'yz' )>0) ) then
        open(1,file=filename,form='unformatted',position='append')
        write(1) a,tslice,pos
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
!  slice position in middle of the box independ of nprocy,nprocz
!  second horizontal slice is the upper most layer
!
      elseif (slice_position=='m') then
        ix_loc=(l1+l2)/2
        if (mod(nprocy,2)==0) then; iy_loc=m1; else; iy_loc=(m1+m2)/2; endif
        if (mod(nprocz,2)==0) then; iz_loc=n1; else; iz_loc=(n1+n2)/2; endif
        iz2_loc=n2
        lwrite_slice_xy2=(ipz==nprocz-1)
        lwrite_slice_xy=(ipz==nprocz/2)
        lwrite_slice_xz=(ipy==nprocy/2)
        lwrite_slice_yz=.true.
!
!  slice positions for spherical coordinates
!  w is for "wedges" since the outputs are 
!  the midplane (rphi,theta=y(mpoint)) and four 
!  wedges in rtheta (xy)
!
      elseif (slice_position=='w') then
        !midplane slices
        !ix_loc=nxgrid/2+nghost
        iy = nygrid/2+nghost
        !meridional wedges, at 4 different 
        !equally spaced azimuthal locations
        iz =  0*nzgrid/4+1+nghost
        iz2=  1*nzgrid/4+1+nghost
        iz3=  2*nzgrid/4+1+nghost
        iz4=  3*nzgrid/4+1+nghost
        !yz is not needed
        lwrite_slice_yz =.false.        
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
!  slice position similar to periphery (p), but the two z positions
!  can now be given in terms of z (zbot_slice, ztop_slice).
!
      elseif (slice_position=='c') then
        ix_loc=l1; iy_loc=m1
        call zlocation(zbot_slice,iz_loc,lwrite_slice_xy)
        call zlocation(ztop_slice,iz2_loc,lwrite_slice_xy2)
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
!  Spherical admits only position 'sph'. Break if this is not met.
!  Also, turn extra r-theta slices to false in case of 
!  non-spherical coordinates
!
      if (coord_system=='spherical') then
        if (slice_position/='w') &
            call fatal_error("setup_slices",&
            "You are using spherical coordinates. "//&
            "Switch slice_position='w' in run_pars")
      else
        lwrite_slice_xy3=.false.
        lwrite_slice_xy4=.false.
      endif
!
!  Overwrite slice postions if any ix,iy,iz,iz2,iz3,iz4 is greater then Zero
!
      if (ix>0) then
        ix_loc=ix-ipx*nx
        if (ix_loc>=l1.and.ix_loc<=l2) then
          lwrite_slice_yz=.true.
        else
          lwrite_slice_yz=.false.
        endif
      endif
!
      if (iy>0) then
        iy_loc=iy-ipy*ny
        if (iy_loc>=m1.and.iy_loc<=m2) then
          lwrite_slice_xz=.true.
        else
          lwrite_slice_xz=.false.
        endif
      endif
!
      if (iz>0) then
        iz_loc=iz-ipz*nz
        if (iz_loc>=n1.and.iz_loc<=n2) then
          lwrite_slice_xy=.true.
        else
          lwrite_slice_xy=.false.
        endif
      endif
!
      if (iz2>0) then
        iz2_loc=iz2-ipz*nz
        if (iz2_loc>=n1.and.iz2_loc<=n2) then
          lwrite_slice_xy2=.true.
        else
          lwrite_slice_xy2=.false.
        endif
      endif
!
      if (iz3>0) then
        iz3_loc=iz3-ipz*nz
        if (iz3_loc>=n1.and.iz3_loc<=n2) then
          lwrite_slice_xy3=.true.
        else
          lwrite_slice_xy3=.false.
        endif
      endif
!
      if (iz4>0) then
        iz4_loc=iz4-ipz*nz
        if (iz4_loc>=n1.and.iz4_loc<=n2) then
          lwrite_slice_xy4=.true.
        else
          lwrite_slice_xy4=.false.
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
!  write the ipz processor numbers for the two slices
!  The first number (=ipz) is essential, the others just for interest.
!
      if (lwrite_slice_xy.and.ipy==0) then
        open(1,file=trim(directory)//'/zbot_procnum.dat',STATUS='unknown')
        write(1,'(2i5,e12.4)') ipz,iz_loc,z(iz_loc)
        close(1)
      endif
      if (lwrite_slice_xy2.and.ipy==0) then
        open(1,file=trim(directory)//'/ztop_procnum.dat',STATUS='unknown')
        write(1,'(2i5,e12.4)') ipz,iz2_loc,z(iz2_loc)
        close(1)
      endif
      if (lwrite_slice_xy3.and.ipy==0) then
        open(1,file=trim(directory)//'/ztop_procnum.dat',STATUS='unknown')
        write(1,'(2i5,e12.4)') ipz,iz3_loc,z(iz3_loc)
        close(1)
      endif
      if (lwrite_slice_xy4.and.ipy==0) then
        open(1,file=trim(directory)//'/ztop_procnum.dat',STATUS='unknown')
        write(1,'(2i5,e12.4)') ipz,iz4_loc,z(iz4_loc)
        close(1)
      endif
!
!  make sure ix_loc,iy_loc,iz_loc,iz2_loc,iz3_loc,iz4_loc 
!  are not outside the boundaries
!
       ix_loc=min( ix_loc,l2) ;  iy_loc=min( iy_loc,m2)
       ix_loc=max( ix_loc,l1) ;  iy_loc=max( iy_loc,m1)
       iz_loc=min( iz_loc,n2) ; iz2_loc=min(iz2_loc,n2)
       iz_loc=max( iz_loc,n1) ; iz2_loc=max(iz2_loc,n1)
      iz3_loc=min(iz3_loc,n2) ; iz4_loc=min(iz4_loc,n2)
      iz3_loc=max(iz3_loc,n1) ; iz4_loc=max(iz4_loc,n1)
!
      if (lroot) then
        write (*,*)'read_runpars: slice_position = '//slice_position
        write (*,'(1x,a,4i4)') &
          'read_runpars: ix,iy,iz,iz2 (video files) =',&
          ix,iy,iz,iz2
        if (coord_system=='spherical') write (*,'(1x,a,2i4)') &
          'read_runpars: iz3, iz4 (video files) =',iz3,iz4
      endif
!
    endsubroutine setup_slices
!***********************************************************************
    subroutine zlocation(zpos,izpos,lproc)
!
!  if zpos lies within this processor, then lproc=T and zpos=z(izpos).
!  Otherwise lproc=F and izpos=1.
!
!  18-nov-06/axel: coded
!  14-oct-08/ccyang: use half-closed interval and include the top-most plane
!
      use Cdata, only: n1,n2,z
!
      real :: zpos
      integer :: izpos,n
      logical :: lproc
!
!  run through all z positions until we hit the right interval.
!  If the right interval is found, jump out of the loop.
!
      do n=n1,n2
        if (z(n)<=zpos.and.z(n+1)>zpos) then
          izpos=n
          lproc=.true.
          goto 900
        else
        endif
      enddo
!
!  if nothing is found, we set lproc=.false. and
!  and put izpos=1
!
      izpos=1
      lproc=.false.
!
900   continue
!
    endsubroutine zlocation
!***********************************************************************

endmodule Slices
