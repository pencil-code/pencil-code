! $Id$
!
! MODULE_DOC: This module contains Yin-Yang related types and functions
! MODULE_DOC: which are incompatible with FORTRAN 95.
!
!***************************************************************
!
module Yinyang
!
  use Cparam, only: impossible, nygrid, nx, pi
  use Cdata, only: iproc_world, lroot, yy_biquad_weights, iproc, lyang

  implicit none

  include 'yinyang.h'

  type ind_coeffs
    integer, dimension(:,:,:), allocatable :: inds
    real, dimension(:,:,:), allocatable :: coeffs
    real, dimension(:,:,:,:), allocatable :: coeffs2
    real, dimension(:,:,:), allocatable :: pcoors
    integer, dimension(:,:,:), allocatable :: icoors
    integer, dimension(:,:,:), allocatable :: igaps
  end type
!
  integer, dimension(:), allocatable :: overlap_thrngs, overlap_phrngs
! 
  integer :: nasym_th=0,nasym_ph=0
!
!  Indicators for situations in which the interpolation stencil overlaps with the
!  present and the neighboring processor domain: L/R - at left/right domain
!  bound; GAP - in gap; MARG/MARG2 - one/two cell(s) away from domain bound; 
!  NEIGH/NEIGH2 - in domain of neighboring processor one/two cell(s) away from domain bound. 
!                
  integer, parameter :: LGAP=-3, RGAP=3, NOGAP=0, LNEIGH=-4, RNEIGH=4, LMARG=-2, RMARG=2, &
                        LMARG2=-1, RMARG2=1, LNEIGH2=-5, RNEIGH2=5
                        !LMARG2=0, RMARG2=0, LNEIGH2=-5, RNEIGH2=5
!
  contains
!
!**************************************************************************
   subroutine biquad_interp(indcoeffs, ith, iph, f, buffer, i2buf, i3buf)
!
!  Performs biquadratic, bicubic or biquintic interpolation for a pencil at position (ith, iph) of the
!  original strip from values in f-array using the precalculated weights
!  in indcoeffs%coeffs. Result is returned in buffer(i2buf,i3buf)=(ith,iph)
!  or buffer(i2buf,i3buf)=(iph,ith), the latter if transposition is required.
!
!  20-dec-15/MR: coded
! 
      use General, only: transform_spher_cart_yy, notanumber,transform_cart_spher, &
                         transform_thph_yy_other, yy_transform_strip_other,transform_cart_spher_other
      use Cdata, only: iproc_world, itsub, m1, m2, n1, n2, lrun

      type(ind_coeffs),         intent(IN) :: indcoeffs
      integer,                  intent(IN) :: ith, iph, i2buf, i3buf
      real, dimension(:,:,:,:), intent(IN) :: f
      real, dimension(:,:,:,:), intent(OUT):: buffer

      integer :: indthl, indthu, indphl, indphu, igp, igt, nphrng, nthrng, outproc
      real, dimension(:,:,:,:), allocatable :: tmp
      real, dimension(size(buffer,1)) :: fmi, fmx, dist
      logical :: l0
      real, dimension(2,1,1) :: thphprime
      real, dimension(nx,3) :: tmp2

      indthl=indcoeffs%inds(ith,iph,1); 
      if (indthl<-1) return
      indthu=indcoeffs%inds(ith,iph,2)

      indphl=indcoeffs%inds(ith,iph,3); indphu=indcoeffs%inds(ith,iph,4)
      nthrng=indthu-indthl+1; nphrng=indphu-indphl+1

      if (size(buffer,4)==3) then
!
!  Vector field, is transformed to Cartesian basis of *other* grid before being interpolated 
!
        allocate (tmp(size(buffer,1),nthrng,nphrng,3)); tmp=0.
        !call transform_spher_cart_yy(f,indthl,indthu,indphl,indphu,tmp,lyy=.true.)
        call transform_thph_yy_other(f,indthl,indthu,indphl,indphu,tmp)

        do igp=1,nphrng; do igt=1,nthrng
          if (indcoeffs%coeffs2(ith,iph,igt,igp)/=0.) then
          buffer(:,i2buf,i3buf,:) = buffer(:,i2buf,i3buf,:)+indcoeffs%coeffs2(ith,iph,igt,igp)*tmp(:,igt,igp,:)
!if (notanumber(tmp(:,igt,igp,:))) print*, 'NaNs at', iproc_world, 'igt,igp=', igt,igp 
          endif
        enddo; enddo

        ! the following tb used when interpolation is done with the Cartesian components: transformed to spherical ones
        !!call yy_transform_strip_other([indcoeffs%pcoors(ith,iph,1)],[indcoeffs%pcoors(ith,iph,2)],thphprime)
        !!tmp2=buffer(:,i2buf,i3buf,:)
        !!call transform_cart_spher_other(tmp2,thphprime(1,1,1),thphprime(2,1,1))
        !!buffer(:,i2buf,i3buf,:)=tmp2
if (notanumber(buffer(:,i2buf,i3buf,:))) print*, 'indthl,indphl, i2buf,i3buf=', indthl,indphl, i2buf,i3buf
      else
!
!  Scalar field - no transformation.
!
        l0=.true.
        do igp=max(indphl,n1),min(indphu,n2)
          do igt=max(indthl,m1),min(indthu,m2)
            if (indcoeffs%coeffs2(ith,iph,igt-indthl+1,igp-indphl+1)/=0.) then
            buffer(:,i2buf,i3buf,1) = buffer(:,i2buf,i3buf,1) &
                                     +indcoeffs%coeffs2(ith,iph,igt-indthl+1,igp-indphl+1)*f(:,igt,igp,1)
!if (notanumber(f(:,igt,igp,1))) print*, 'NaNs at', iproc_world,  &
!'igt,igp=', igt,igp
            if (l0) then
              l0=.false.
              fmx=f(:,igt,igp,1)
              fmi=fmx
            else
              fmx=max(fmx,f(:,igt,igp,1))
              fmi=min(fmi,f(:,igt,igp,1))
            endif
            endif
          enddo
        enddo
!if (notanumber(buffer(:,i2buf,i3buf,1))) print*, 'NaNs at', iproc_world,  &
!  ',indthl,indphl, i2buf,i3buf,fmx=', indthl,indphl,ith,iph,i2buf,i3buf,fmx
!  ',indthl,indphl, i2buf,i3buf,fmx=', indthl,indphl,ith,iph,i2buf,i3buf,fmx

        return  !!!
        outproc=26
        if (iproc_world==outproc.and.itsub==3) then
          write(88,*) outproc
          write(88,*) indcoeffs%inds(ith,iph,:),buffer(63,i2buf,i3buf,1), &
                      indcoeffs%pcoors(ith,iph,:), indcoeffs%icoors(ith,iph,:),& 
                      f(63,indthl:indthl+nthrng-1,indphl:indphl+nphrng-1,1)
        endif
        return  !!!
!
! experimental: avoids over- and undershoot
!
        where (buffer(:,i2buf,i3buf,1)>fmx) buffer(:,i2buf,i3buf,1)=fmx
        where (buffer(:,i2buf,i3buf,1)<fmi) buffer(:,i2buf,i3buf,1)=fmi
        return
        dist=.035
        where (fmx<0.) 
          fmx=(1.-dist)*fmx
          fmi=(1.+dist)*fmi
        elsewhere
          where (fmi>0.) 
            fmi=(1.-dist)*fmi
            fmx=(1.+dist)*fmx
          elsewhere
            fmi=(1.+dist)*fmi
            fmx=(1.+dist)*fmx
          endwhere
        endwhere
        if (any(buffer(:,i2buf,i3buf,1)>fmx)) then
          if (.not.lyang) print*, 'iproc', iproc, 'overshoot'
        endif
        if (any(buffer(:,i2buf,i3buf,1)<fmi)) then
          if (.not.lyang) print*, 'iproc', iproc, 'undershoot'
        endif
      endif

    endsubroutine biquad_interp
!**************************************************************************
   subroutine bilin_interp(indcoeffs, ith, iph, f, buffer, i2buf, i3buf)
!
!  Performs bilinear interpolation for a pencil at position (ith, iph) of the
!  original strip from values in f-array using the precalculated weights
!  in indcoeffs%coeffs. Result is returned in buffer(i2buf,i3buf)=(ith,iph)
!  or buffer(i2buf,i3buf)=(iph,ith), the latter if transposition is required.
!
!  20-dec-15/MR: coded
! 
      use General, only: transform_spher_cart_yy, notanumber
      use Cdata, only: lrun, nx, lyang, iproc

      type(ind_coeffs),         intent(IN) :: indcoeffs
      integer,                  intent(IN) :: ith, iph, i2buf, i3buf
      real, dimension(:,:,:,:), intent(IN) :: f
      real, dimension(:,:,:,:), intent(OUT):: buffer

      integer :: indth, indph, indx
      real, dimension(:,:,:,:), allocatable :: tmp

      indth=indcoeffs%inds(ith,iph,1)
      indph=indcoeffs%inds(ith,iph,2)

      if (indth<-1 .or. indph<-1) return    ! point at (ith,iph) not caught by calling proc.

      if (size(buffer,4)==3) then
!
!  Vector field, is transformed to Cartesian basis of *other* grid before being interpolated 
!
        allocate (tmp(size(buffer,1),2,2,3)); tmp=0.
        call transform_spher_cart_yy(f,indth-1,indth,indph-1,indph,tmp,lyy=.true.)

        buffer(:,i2buf,i3buf,:) = indcoeffs%coeffs(ith,iph,1)*tmp(:,1,1,:) &
                                 +indcoeffs%coeffs(ith,iph,2)*tmp(:,1,2,:) &
                                 +indcoeffs%coeffs(ith,iph,3)*tmp(:,2,1,:) &
                                 +indcoeffs%coeffs(ith,iph,4)*tmp(:,2,2,:)
if (notanumber(buffer(:,i2buf,i3buf,:))) print*, 'Vector: NaNs at', iproc_world,  &
  ',indth,indph,ith,iph,i2buf,i3buf=', indth,indph,ith,iph,i2buf,i3buf
      else
!
!  Scalar field - no transformation.
!
        buffer(:,i2buf,i3buf,1) = indcoeffs%coeffs(ith,iph,1)*f(:,indth-1,indph-1,1) &
                                 +indcoeffs%coeffs(ith,iph,2)*f(:,indth-1,indph  ,1) &
                                 +indcoeffs%coeffs(ith,iph,3)*f(:,indth  ,indph-1,1) &
                                 +indcoeffs%coeffs(ith,iph,4)*f(:,indth  ,indph  ,1)

if (notanumber(buffer(:,i2buf,i3buf,1))) print*, 'Scalar: NaNs at', iproc_world,  &
  ',indth,indph,ith,iph,i2buf,i3buf=', indth,indph,ith,iph,i2buf,i3buf

if (.false..and.lrun.and..not.lyang) then
  write(iproc+60,*) indth, indph
  indx=nx/2
  write(iproc+60,*) buffer(indx,i2buf,i3buf,1), f(indx,indth-1,indph-1,1), f(indx,indth,indph-1,1), &
                    f(indx,indth-1,indph,1), f(indx,indth,indph,1)
endif
      endif

    endsubroutine bilin_interp
!***********************************************************************
    subroutine prep_biquad_interp(pprime,indth,indph,indcoeffs,ip,jp,ma,me,na,ne)

      use Cdata, only: y, z, dy, dz, lroot

      real, dimension(2),intent(IN)   :: pprime
      integer,           intent(IN)   :: indth,indph,ip,jp,ma,me,na,ne
      type(ind_coeffs),  intent(INOUT):: indcoeffs

      real, dimension(4) :: qth, qph
      integer :: i, igt, igp, nthrng, nphrng
      real, dimension(3) :: gth, gph3, gph4
      real :: dthp, dphp, w1, w2, w1ph, w2ph, w2th
      logical :: lextweights

      lextweights = all(yy_biquad_weights/=impossible)
      if (lextweights) then
        if (abs(yy_biquad_weights(1)+yy_biquad_weights(2) - 1.)>1e-6 .or. &
            abs(yy_biquad_weights(3)+yy_biquad_weights(4) - 1.)>1e-6 ) then 
           if (lroot) print*, 'prep_biquad_interp: Warning - sum of weights in yy_biquad_weights not unity'
        endif
      endif
!
      nthrng=3
      if (indth==ma+1) then
        indcoeffs%inds(ip,jp,1:2)=(/ma,ma+2/)
      elseif (indth==me) then
        indcoeffs%inds(ip,jp,1:2)=(/me-2,me/)
      else
        nthrng=4; indcoeffs%inds(ip,jp,1:2)=(/indth-2,indth+1/)
      endif

      nphrng=3
      if (indph==na+1) then
        indcoeffs%inds(ip,jp,3:4)=(/na,na+2/)
      elseif (indph==ne) then
        indcoeffs%inds(ip,jp,3:4)=(/ne-2,ne/)
      else
        nphrng=4; indcoeffs%inds(ip,jp,3:4)=(/indph-2,indph+1/)
      endif

      do i=1,nthrng
        qth(i) = (pprime(1)-y(i+indcoeffs%inds(ip,jp,1)-1))/dy   ! requires equidistant grid in y
      enddo

      do i=1,nphrng
        qph(i) = (pprime(2)-z(i+indcoeffs%inds(ip,jp,3)-1))/dz   ! requires equidistant grid in z
      enddo

      gth =(/.5*qth(2)*qth(3),-qth(1)*qth(3),.5*qth(1)*qth(2)/)
      gph3=(/.5*qph(2)*qph(3),-qph(1)*qph(3),.5*qph(1)*qph(2)/)
!print*, 'gth,gph3=', sum(gth), sum(gph3)
      do igp=1,3; do igt=1,3
        indcoeffs%coeffs2(ip,jp,igt,igp)=gth(igt)*gph3(igp)
      enddo; enddo
!print*, 'indcoeffs33=', sum(indcoeffs%coeffs2(ip,jp,:,:))

      if (nphrng==4) then

        gph4=(/.5*qph(3)*qph(4),-qph(2)*qph(4),.5*qph(2)*qph(3)/)
!print*, 'gph4=', sum(gph4)
        if (lextweights) then
          w1ph=yy_biquad_weights(3); w2ph=yy_biquad_weights(4)
        else
          w1ph=-qph(3); w2ph=qph(2)
!if (lroot) print*, 'w1ph, w2ph=', w1ph, w2ph
        endif
        indcoeffs%coeffs2(ip,jp,1:3,1:3)=w1ph*indcoeffs%coeffs2(ip,jp,1:3,1:3)

        do igp=1,3; do igt=1,3
          indcoeffs%coeffs2(ip,jp,igt,igp+1)=indcoeffs%coeffs2(ip,jp,igt,igp+1)+w2ph*gth(igt)*gph4(igp)
        enddo; enddo
!print*, 'indcoeffs34=', sum(indcoeffs%coeffs2(ip,jp,:,:))

      endif

      if (nthrng==4) then

        gth=(/.5*qth(3)*qth(4),-qth(2)*qth(4),.5*qth(2)*qth(3)/)
!print*, 'gth4=', sum(gth)
        if (lextweights) then
          w1=yy_biquad_weights(1); w2th=yy_biquad_weights(2)
        else
          w1=-qth(3); w2th=qth(2)
!if (lroot) print*, 'w1th, w2th=', w1, w2th
        endif
        indcoeffs%coeffs2(ip,jp,1:3,:)=w1*indcoeffs%coeffs2(ip,jp,1:3,:)
       
        if (nphrng==4) then; w2=w2th*w1ph; else; w2=w2th; endif

        do igp=1,3; do igt=1,3
          indcoeffs%coeffs2(ip,jp,igt+1,igp)=indcoeffs%coeffs2(ip,jp,igt+1,igp)+w2*gth(igt)*gph3(igp)
        enddo; enddo
!if (nphrng==3) print*, 'indcoeffs43=', sum(indcoeffs%coeffs2(ip,jp,:,:))

        if (nphrng==4) then
          w2=w2th*w2ph
          do igp=1,3; do igt=1,3
            indcoeffs%coeffs2(ip,jp,igt+1,igp+1)=indcoeffs%coeffs2(ip,jp,igt+1,igp+1)+w2*gth(igt)*gph4(igp)
          enddo; enddo
!print*, 'indcoeffs44=', sum(indcoeffs%coeffs2(ip,jp,:,:))
        endif

      endif

    endsubroutine prep_biquad_interp
!***********************************************************************
    subroutine prep_bicub_interp(pprime,indth,indph,indcoeffs,ip,jp,igapt,igapp)
!
!  Establishes data needed for the bicubic interpolation.
!
!  12-jan-2016/MR: coded
!  12-dec-2017/MR: modified for fixing the gap problem
!  10-jan-2018/MR: now on each side of a processor domain three situations
!                  considered in which the interpolation stencil overlaps with
!                  the present and the neighboring processor domain. Then the array
!                  of interpolation coefficients is only partly filled for both.
!                  Situations with asymmetric interpolation detected - warning
!                  launched when present.
!
      use Cdata, only: m1, m2, n1, n2, y, z, dy, dz, &
                       lfirst_proc_y, llast_proc_y, lfirst_proc_z, llast_proc_z

      real, dimension(2),intent(IN)   :: pprime
      integer,           intent(IN)   :: indth,indph,ip,jp,igapt,igapp
      type(ind_coeffs),  intent(INOUT):: indcoeffs

      real, dimension(4) :: qth, qph, gth, gph
      integer :: i, igt, igp

      if (indth==m1+1.and.igapt==NOGAP) then
        indcoeffs%inds(ip,jp,1:2)=(/m1,m1+3/); nasym_th=nasym_th+1
      elseif (indth==m2.and.igapt==NOGAP) then
        indcoeffs%inds(ip,jp,1:2)=(/m2-3,m2/); nasym_th=nasym_th+1
      else
        indcoeffs%inds(ip,jp,1:2)=(/indth-2,indth+1/)
      endif

      if (indph==n1+1.and.igapp==NOGAP) then
        indcoeffs%inds(ip,jp,3:4)=(/n1,n1+3/); nasym_ph=nasym_ph+1
      elseif (indph==n2.and.igapp==NOGAP) then
        indcoeffs%inds(ip,jp,3:4)=(/n2-3,n2/); nasym_ph=nasym_ph+1
      else
        indcoeffs%inds(ip,jp,3:4)=(/indph-2,indph+1/)
      endif

      do i=1,4
        qth(i) = (pprime(1)-y(indcoeffs%inds(ip,jp,1)+i-1))/dy   ! requires equidistant grid in y
        qph(i) = (pprime(2)-z(indcoeffs%inds(ip,jp,3)+i-1))/dz   ! requires equidistant grid in z
      enddo

      gth=0.; gph=0.

      if (igapt>=NOGAP) &
        gth(1)=-.166666*qth(2)*qth(3)*qth(4)
      
      if ((igapt>=NOGAP.and.igapt<=RGAP).or.igapt==LMARG) &
        gth(2)=.5*qth(1)*qth(3)*qth(4)
      
      if ((igapt<=NOGAP.and.igapt>=LGAP).or.igapt==RMARG) &
        gth(3)=-.5*qth(1)*qth(2)*qth(4)

      if (igapt<=NOGAP) &
        gth(4)=.166666*qth(1)*qth(2)*qth(3)

      if (igapp>=NOGAP) &
        gph(1)=-.166666*qph(2)*qph(3)*qph(4)

      if ((igapp>=NOGAP.and.igapp<=RGAP).or.igapp==LMARG) &
        gph(2)=.5*qph(1)*qph(3)*qph(4)

      if ((igapp<=NOGAP.and.igapp>=LGAP).or.igapp==RMARG) &
        gph(3)=-.5*qph(1)*qph(2)*qph(4)

      if (igapp<=NOGAP) &
        gph(4)=.166666*qph(1)*qph(2)*qph(3)

      do igp=1,4; do igt=1,4
        indcoeffs%coeffs2(ip,jp,igt,igp)=gth(igt)*gph(igp)
      enddo; enddo
      
      indcoeffs%pcoors(ip,jp,:)=pprime
      indcoeffs%icoors(ip,jp,:)=(/indth,indph/)
      indcoeffs%igaps(ip,jp,:)=(/igapt,igapp/)
!if (.not.lyang) then
if (.false.) then
if (igapt==NOGAP.and.igapp==NOGAP) then
  if (abs(sum(indcoeffs%coeffs2(ip,jp,1:4,1:4))-1.)>5.e-7) &
    print*, 'coefficient sum /=1: iproc,ip,jp,residual=', iproc,ip,jp,sum(indcoeffs%coeffs(ip,jp,:))-1.
endif
endif

    endsubroutine prep_bicub_interp
!***********************************************************************
    subroutine prep_biquint_interp(pprime,indth,indph,indcoeffs,ip,jp,igapt,igapp)
!
!  Establishes data needed for the biquintic interpolation.
!
!  12-jan-2018/MR: coded. 
!                  On each side of a processor domain there are five situations
!                  possible in which the interpolation stencil overlaps with the
!                  present and the neighboring processor domain. Then the array
!                  of interpolation coefficients is only partly filled for both.
!
      use Cdata, only: m1, m2, n1, n2, y, z, dy, dz, ygrid, zgrid, ny, nz, nghost, ipy, ipz

      real, dimension(2),intent(IN)   :: pprime
      integer,           intent(IN)   :: indth,indph,ip,jp,igapt,igapp
      type(ind_coeffs),  intent(INOUT):: indcoeffs

      real, dimension(6) :: qth, qph, gth, gph
      integer :: i, igt, igp, iyu, izu,lun
      logical, save :: l0=.true.
      logical :: lnasym_th, lnasym_ph

      lnasym_th=.false.; lnasym_ph=.false.; 
      if (igapt==NOGAP.and.indth<=m1+2) then
        indcoeffs%inds(ip,jp,1:2)=(/m1,m1+5/); nasym_th=nasym_th+1; lnasym_th=.true.
      elseif (igapt==NOGAP.and.indth>=m2-1) then
        indcoeffs%inds(ip,jp,1:2)=(/m2-5,m2/); nasym_th=nasym_th+1; lnasym_th=.true.
      else
        indcoeffs%inds(ip,jp,1:2)=(/indth-3,indth+2/)
      endif

      if (igapp==NOGAP.and.indph<=n1+2) then
        indcoeffs%inds(ip,jp,3:4)=(/n1,n1+5/); nasym_ph=nasym_ph+1; lnasym_ph=.true.
      elseif (igapp==NOGAP.and.indph>=n2-1) then
        indcoeffs%inds(ip,jp,3:4)=(/n2-5,n2/); nasym_ph=nasym_ph+1; lnasym_ph=.true.
      else
        indcoeffs%inds(ip,jp,3:4)=(/indph-3,indph+2/)
      endif

      iyu=ipy*ny+indcoeffs%inds(ip,jp,1)-nghost
      izu=ipz*nz+indcoeffs%inds(ip,jp,3)-nghost
      do i=1,6
        qth(i) = (pprime(1)-ygrid(iyu+i-1))/dy   ! requires equidistant grid in y
        qph(i) = (pprime(2)-zgrid(izu+i-1))/dz   ! requires equidistant grid in z
      enddo
!if (.not.lnasym_th.and.(any(qth(4:)>0.).or.any(qth(1:3)<0.))) print*, 'qth<0',ip,indth 
!if (.not.lnasym_ph.and.(any(qph(4:)>0.).or.any(qph(1:3)<0.))) print*, 'qph<0',jp,indph 
      gth=0.; gph=0.
!
! interpolation weights: 1/([-120,24,-12,12,-24,120]*dx^5) = [-0.00833333,0.0416667,-0.0833333,0.0833333,-0.0416667,0.00833333]*dx^-5
!
      if (igapt>=NOGAP) &
        gth(1)=-0.00833333*qth(2)*qth(3)*qth(4)*qth(5)*qth(6)

      if ((igapt>=NOGAP.and.igapt<=RNEIGH).or.igapt==LMARG2) &
        gth(2)=0.0416667*qth(1)*qth(3)*qth(4)*qth(5)*qth(6)

      if ((igapt>=NOGAP.and.igapt<=RGAP).or.igapt==LMARG.or.igapt==LMARG2) &
        gth(3)=-0.0833333*qth(1)*qth(2)*qth(4)*qth(5)*qth(6)

      if ((igapt<=NOGAP.and.igapt>=LGAP).or.igapt==RMARG.or.igapt==RMARG2) &
        gth(4)=0.0833333*qth(1)*qth(2)*qth(3)*qth(5)*qth(6)

      if ((igapt<=NOGAP.and.igapt>=LNEIGH).or.igapt==RMARG2) &
        gth(5)=-0.0416667*qth(1)*qth(2)*qth(3)*qth(4)*qth(6)

      if (igapt<=NOGAP) &
        gth(6)=0.00833333*qth(1)*qth(2)*qth(3)*qth(4)*qth(5)

      if (igapp>=NOGAP) &
        gph(1)=-0.00833333*qph(2)*qph(3)*qph(4)*qph(5)*qph(6)

      if ((igapp>=NOGAP.and.igapp<=RNEIGH).or.igapp==LMARG2) &
        gph(2)=0.0416667*qph(1)*qph(3)*qph(4)*qph(5)*qph(6)

      if ((igapp>=NOGAP.and.igapp<=RGAP).or.igapp==LMARG.or.igapp==LMARG2) &
        gph(3)=-0.0833333*qph(1)*qph(2)*qph(4)*qph(5)*qph(6)

      if ((igapp<=NOGAP.and.igapp>=LGAP).or.igapp==RMARG.or.igapp==RMARG2) &
        gph(4)=0.0833333*qph(1)*qph(2)*qph(3)*qph(5)*qph(6)

      if ((igapp<=NOGAP.and.igapp>=LNEIGH).or.igapp==RMARG2) &
        gph(5)=-0.0416667*qph(1)*qph(2)*qph(3)*qph(4)*qph(6)

      if (igapp<=NOGAP) &
        gph(6)=0.00833333*qph(1)*qph(2)*qph(3)*qph(4)*qph(5)

      do igp=1,6; do igt=1,6
        indcoeffs%coeffs2(ip,jp,igt,igp)=gth(igt)*gph(igp)
      enddo; enddo

if (.false.) then
!if (.not.lyang) then
!if (iproc==0.and.igapt==1.or.iproc==1.and.igapt==-5) then   !.or.iproc==1) then
if (iproc==0.and.igapp==1.or.iproc==3.and.igapp==-5) then   !.or.iproc==1) then
!if (iproc==3.or.iproc==6) then
  !if (igapt==LMARG.or.igapt==RMARG.or.igapt==LNEIGH.or.igapt==RNEIGH) &
  !if (igapt==LNEIGH2.or.igapt==RNEIGH2) &
  !if (igapt==LMARG.or.igapt==LGAP) &
  !if (igapt==-3.or.igapt==-3) then
  if (iproc==0) then
    lun=40
  else
    lun=41
  endif
  !write(lun,'(a,2i2,1x,2(f7.3,1x))') 'iproc, igapt, pprime=', iproc, igapt, pprime
  write(lun,'(a,2i2,1x,2(f7.3,1x))') 'iproc, igapp, pprime=', iproc, igapp, pprime
  !print'(4(i4,","))', indcoeffs%inds(ip,jp,:)
  !print'(6(e11.4,","))', qth
  !print'(6(e11.4,","))', gth
  !print'(6(e11.4,","))', gph
  write(lun,*) '----'
  write(lun,'(6(e11.4,","))') indcoeffs%coeffs2(ip,jp,:,:)
  !endif
  !print'(a,2i2,1x,2(f7.3,1x))', 'iproc, igapp, pprime=', iproc, igapp, pprime
  !if (igapp==LMARG.or.igapp==RMARG.or.igapp==LNEIGH.or.igapp==RNEIGH) &
  !print'(a,2i2,1x,2(f7.3,1x))', 'iproc, igapp, pprime=', iproc, igapp, pprime
endif
endif
if (.false.) then
!if (.true.) then
if (.not.lyang) then
!if (indcoeffs%inds(ip,jp,2)-indcoeffs%inds(ip,jp,1)+1==6.and.indcoeffs%inds(ip,jp,4)-indcoeffs%inds(ip,jp,3)+1==6) then
!if (igapt==NOGAP.and.igapp==NOGAP) then
  if (abs(sum(indcoeffs%coeffs2(ip,jp,:,:))-1.)>1.5e-6) &
    print*, 'coefficient sum /=1: iproc,ip,jp,residual=', iproc,ip,jp,sum(indcoeffs%coeffs2(ip,jp,:,:))-1.
!endif
endif
endif
      indcoeffs%pcoors(ip,jp,:)=pprime
      indcoeffs%igaps(ip,jp,:)=(/igapt,igapp/)
      indcoeffs%icoors(ip,jp,:)=(/indth,indph/)

return
gth=(/-0.00833333*qth(2)*qth(3)*qth(4)*qth(5)*qth(6), &
      0.0416667*qth(1)*qth(3)*qth(4)*qth(5)*qth(6), &
      -0.0833333*qth(1)*qth(2)*qth(4)*qth(5)*qth(6), &
      0.0833333*qth(1)*qth(2)*qth(3)*qth(5)*qth(6), &
      -0.0416667*qth(1)*qth(2)*qth(3)*qth(4)*qth(6), &
      0.00833333*qth(1)*qth(2)*qth(3)*qth(4)*qth(5) /)

gph=(/-0.00833333*qph(2)*qph(3)*qph(4)*qph(5)*qph(6), &
      0.0416667*qph(1)*qph(3)*qph(4)*qph(5)*qph(6), &
      -0.0833333*qph(1)*qph(2)*qph(4)*qph(5)*qph(6), &
      0.0833333*qph(1)*qph(2)*qph(3)*qph(5)*qph(6), &
      -0.0416667*qph(1)*qph(2)*qph(3)*qph(4)*qph(6), &
      0.00833333*qph(1)*qph(2)*qph(3)*qph(4)*qph(5) /)
    endsubroutine prep_biquint_interp
!***********************************************************************
    subroutine prep_quadspline_interp(pprime,indth,indph,indcoeffs,ip,jp,ma,me,na,ne)
!
!  Establishes data needed for the quadratic spline interpolation.
!
!  12-jan-2016/MR: coded
!
      use Cdata, only: y, z, dy, dz, lroot

      real, dimension(2),intent(IN)   :: pprime
      integer,           intent(IN)   :: ip,jp,ma,me,na,ne
      integer,           intent(INOUT):: indth,indph
      type(ind_coeffs),  intent(INOUT):: indcoeffs

      real :: qth, qph, sum
      integer :: igt, igp
      real, dimension(3) :: gth, gph

      if ( pprime(1)-y(indth-1) < y(indth)-pprime(1) ) indth=indth-1

      if (indth<=ma+1) then
        indcoeffs%inds(ip,jp,1:2)=(/ma,ma+2/); nasym_th=nasym_th+1
      elseif (indth==me) then
        indcoeffs%inds(ip,jp,1:2)=(/me-2,me/); nasym_th=nasym_th+1
      else
        indcoeffs%inds(ip,jp,1:2)=(/indth-1,indth+1/)
      endif

      if ( pprime(2)-z(indph-1) < z(indph)-pprime(2) ) indph=indph-1

      if (indph<=na+1) then
        indcoeffs%inds(ip,jp,3:4)=(/na,na+2/); nasym_ph=nasym_ph+1
      elseif (indph==ne) then
        indcoeffs%inds(ip,jp,3:4)=(/ne-2,ne/); nasym_ph=nasym_ph+1
      else
        indcoeffs%inds(ip,jp,3:4)=(/indph-1,indph+1/)
      endif

      qth = (pprime(1)-y(indcoeffs%inds(ip,jp,1)+1))/dy   ! requires equidistant grid in y
      qph = (pprime(2)-z(indcoeffs%inds(ip,jp,3)+1))/dz   ! requires equidistant grid in z
      !qth = (pprime(1)-y(indth))/dy   ! requires equidistant grid in y
      !qph = (pprime(2)-z(indph))/dz   ! requires equidistant grid in z
if (abs(qth)>1..or.abs(qph)>1.) print*, 'iproc_world, sum=', iproc_world,qth,qph

      gth =(/0.5*(0.5-qth)**2,0.75-qth**2,0.5*(0.5+qth)**2/)
      gph =(/0.5*(0.5-qph)**2,0.75-qph**2,0.5*(0.5+qph)**2/)

      sum=0.
      do igp=1,3; do igt=1,3
        indcoeffs%coeffs2(ip,jp,igt,igp)=gth(igt)*gph(igp)
        sum=sum+gth(igt)*gph(igp)
      enddo; enddo
if (abs(sum-1.) > 1.e-5) print*, 'iproc_world, sum=', iproc_world, sum
    endsubroutine prep_quadspline_interp
!***********************************************************************
    function qualify_position_biquin(ind,nl,nu,lrestr_par_low,lrestr_par_up,lrestr_perp) result (qual)
!
! Qualifies the position of a point w.r.t the neighboring processors for quintic
! interpolation (six-points-stencil).
! If the calling proc is the first (lrestr_par_low=T) or last (lrestr_par_up=T)
! proc in the considered direction, a position near the margin of its domain is
! not special, likewise not if it is not the first or last in the perpendicular
! direction (lrestr_perp=F).
!
      integer, intent(IN) :: ind, nl, nu
      logical, intent(IN) :: lrestr_par_low,lrestr_par_up,lrestr_perp
      integer :: qual

      if (ind==nl-2) then       ! in penultimate grid cell of left neighbor's domain
        qual=LNEIGH2
      elseif (ind==nl-1) then   ! in last grid cell of left neighbor's domain
        qual=LNEIGH
      elseif (ind==nl) then     ! in gap to left neighbor's domain
        qual=LGAP
      elseif (ind==nl+1.and..not.lrestr_par_low.and.lrestr_perp) then  ! in first grid cell of calling proc
        qual=LMARG
      elseif (ind==nl+2.and..not.lrestr_par_low.and.lrestr_perp) then  ! in second grid cell of calling proc 
        qual=LMARG2
      elseif (ind==nu-1.and..not.lrestr_par_up.and.lrestr_perp) then   ! in penultimate grid cell of calling proc
        qual=RMARG2
      elseif (ind==nu.and..not.lrestr_par_up.and.lrestr_perp) then     ! in last grid cell of calling proc
        qual=RMARG
      elseif (ind==nu+1) then   ! in gap to right neighbor's domain
        qual=RGAP
      elseif (ind==nu+2) then   ! in first grid cell of right neighbor's domain
        qual=RNEIGH
      elseif (ind==nu+3) then   ! in second grid cell of right neighbor's domain
        qual=RNEIGH2
      else
        qual=NOGAP
      endif

    endfunction qualify_position_biquin
!***********************************************************************
    function qualify_position_bicub(ind,nl,nu,lrestr_par_low,lrestr_par_up,lrestr_perp) result (qual)
!
! Qualifies the position of a point w.r.t the neighboring processors for
! cubic interpolation (four-points-stencil).
!
      integer, intent(IN) :: ind, nl, nu
      logical, intent(IN) :: lrestr_par_low,lrestr_par_up,lrestr_perp
      integer :: qual

      if (ind==nl-1) then
        qual=LNEIGH
      elseif (ind==nl) then
        qual=LGAP 
      elseif (ind==nl+1.and..not.lrestr_par_low.and.lrestr_perp) then
        qual=LMARG
      elseif (ind==nu.and..not.lrestr_par_up.and.lrestr_perp) then
        qual=RMARG
      elseif (ind==nu+1) then
        qual=RGAP 
      elseif (ind==nu+2) then
        qual=RNEIGH
      else
        qual=NOGAP
      endif     

    endfunction qualify_position_bicub
!***********************************************************************
    function qualify_position_bilin(ind,nl,nu) result (qual)
!
! Qualifies the position of a point w.r.t the neighboring processors for
! linear interpolation (two-points-stencil).
!
      integer, intent(IN) :: ind, nl, nu
      integer :: qual

      if (ind==nl) then
        qual=LGAP
      elseif (ind==nu+1) then
        qual=RGAP
      else
        qual=NOGAP
      endif

    endfunction qualify_position_bilin
!***********************************************************************
    function prep_interp(thphprime,indcoeffs,itype,ngap,th_range) result (nok)
!
!  For each of the points in the strip thphprime (with shape 2 x thprime-extent x
!  phprime-extent), arbitrarily positioned in the yz-plane, determine in
!  which cell of the grid y(ma:me) x z(na:ne) it lies, store indices of the
!  cell's upper right corner in indcoeffs%inds and the weights of bilinear
!  interpolation for the four corners in indcoeffs%coeffs. If no cell is found
!  for a point, indcoeffs%inds and indcoeffs%coeffs are set zero.
!  If present, return in th_range the interval in thprime-extent in which
!  interpolation cells could be assigned to any points.
!  Returns number of points in thphprime for which interpolation cell could be
!  found.
!
!  20-dec-15/MR: coded
!  12-dec-17/MR: modifications for fixing the gap problem. Determination of 
!                index ranges for overlap regions added.
!  20-jan-18/MR: biquintic interpolation added: now on each side five situations
!                exist in which two processors contribute to interpoland.
!
      use General, only: find_index_range_hill, notanumber, itoa
      use Cdata,   only: ny, nz, y, z, dy, dz, lfirst_proc_y, lfirst_proc_z, nghost, &
                         llast_proc_y, llast_proc_z, lroot, iproc, lyang, lstart
      use Cparam,  only: m1,m2,n1,n2, BILIN, BIQUAD, BICUB, QUADSPLINE, BIQUIN

      real, dimension(:,:,:),          intent(IN)   :: thphprime
      type(ind_coeffs),                intent(OUT)  :: indcoeffs
      integer,                         intent(IN)   :: itype
      integer, dimension(2), optional, intent(OUT)  :: th_range
      integer,               optional, intent(INOUT):: ngap

      integer :: nok

      integer :: ip, jp, indth, indph, sz1, sz2, ma, na, me, ne, i1, i2, indth_in, indph_in, ind
      real :: dth, dph, dthp, dphp, qth1, qth2, qph1, qph2
      real, parameter :: eps=0.
      logical :: okt, okp
      integer :: igapt, igapp
      character(LEN=240) :: txt

      if (lyang) then
        if (.not.allocated(overlap_thrngs)) then
          if (lfirst_proc_y .or. llast_proc_y) then
            allocate(overlap_thrngs(nz)); overlap_thrngs=0
          endif
        endif
        if (.not.allocated(overlap_phrngs)) then
          if (lfirst_proc_z .or. llast_proc_z) then
            allocate(overlap_phrngs(ny)); overlap_phrngs=0
          endif
        endif
      endif

      sz1=size(thphprime,2); sz2=size(thphprime,3)

      if (allocated(indcoeffs%inds)) deallocate(indcoeffs%inds)
      if (allocated(indcoeffs%coeffs)) deallocate(indcoeffs%coeffs)

      if (itype==BIQUAD.or.itype==BICUB.or.itype==BIQUIN.or.itype==QUADSPLINE) then
        if (allocated(indcoeffs%pcoors)) deallocate(indcoeffs%pcoors,indcoeffs%icoors,indcoeffs%igaps)
        allocate(indcoeffs%inds(sz1,sz2,4))
        allocate(indcoeffs%pcoors(sz1,sz2,2))
        allocate(indcoeffs%icoors(sz1,sz2,2))
        allocate(indcoeffs%igaps(sz1,sz2,2))
      elseif (itype==BILIN) then
        allocate(indcoeffs%inds(sz1,sz2,2))
      else
        if (lroot) print*, &
          'prep_interp: Only bilinear, biquadratic, bicubic, biquintic and quadratic spline interpolations implemented!'
        stop
      endif
!
!  Formally, the lower bound of the index range for interpolation can reach -1
!  (for biquintic interpolation), so -2 indicates "point not caught".
!
      indcoeffs%inds=-2

      nok=0
!
!  In order to catch points which lie in gaps between the domains of procs or in
!  the domain of a neighboring processor, but with interpolation contributions
!  yet from the present processor 
!  the lower and upper cells in either the y or the z directions are included, 
!  with exceptions for the first and last processors.
!
      if (present(th_range)) then        ! this for preparation of z averages (only valid for BILIN interpol)
        ma=m1-1; me=m2+1
      else
        ma=m1; me=m2
        if (lfirst_proc_z.or.llast_proc_z) then

          if (itype==BIQUIN) then
            if (.not.lfirst_proc_y) ma=m1-3
            if (.not.llast_proc_y) me=m2+3
          elseif (itype==BICUB.or.itype==BIQUAD) then
            if (.not.lfirst_proc_y) ma=m1-2
            if (.not.llast_proc_y) me=m2+2
          else     ! bilinear interpolation; quadsplines?
            if (.not.lfirst_proc_y) ma=m1-1
            if (.not.llast_proc_y) me=m2+1
          endif

        endif
      endif
!
!  Likewise for z direction.
!
      if (present(th_range)) then        ! this for preparation of z averages (only valid for BILIN interpol)
        na=n1-1; ne=n2+1
      else
        na=n1; ne=n2
        if (lfirst_proc_y.or.llast_proc_y) then

          if (itype==BIQUIN) then
            if (.not.lfirst_proc_z) na=n1-3
            if (.not.llast_proc_z) ne=n2+3
          elseif (itype==BICUB.or.itype==BIQUAD) then
            if (.not.lfirst_proc_z) na=n1-2
            if (.not.llast_proc_z) ne=n2+2
          else
            if (.not.lfirst_proc_z) na=n1-1
            if (.not.llast_proc_z) ne=n2+1
          endif

        endif
      endif

      do ip=1,sz1
        do jp=1,sz2
          !!if (iproc_world==0) write(iproc_world+300,*) thphprime(:,ip,jp)
!
!  For all points in strip,
!
          okt=.false.; okp=.false.
          !!!if ( y(ma)<=thphprime(1,ip,jp) ) then
          if ( y(ma)-eps*dy<=thphprime(1,ip,jp) ) then
!
!  detect between which y lines it lies.
!
            do indth=ma+1,me
              if ( y(indth)+eps*dy>=thphprime(1,ip,jp) ) then
              !!!if ( y(indth)>=thphprime(1,ip,jp) ) then
                okt=.true.
                if (itype==BIQUIN) then
                  igapt=qualify_position_biquin(indth,m1,m2,lfirst_proc_y,llast_proc_y,lfirst_proc_z.or.llast_proc_z) 
                elseif (itype==BICUB.or.itype==BIQUAD) then
                  igapt=qualify_position_bicub(indth,m1,m2,lfirst_proc_y,llast_proc_y,lfirst_proc_z.or.llast_proc_z) 
                else
                  igapt=qualify_position_bilin(indth,m1,m2)
                endif
                exit
              endif
            enddo
          endif

          if (okt) then
            !!!if ( z(na)<=thphprime(2,ip,jp) ) then
            if ( z(na)-eps*dz<=thphprime(2,ip,jp) ) then
!
!  detect between which z lines it lies.
!
              do indph=na+1,ne
                if ( z(indph)+eps*dz>=thphprime(2,ip,jp) ) then
                !!!if ( z(indph)>=thphprime(2,ip,jp) ) then
                  okp=.true.
                  indcoeffs%inds(ip,jp,1:2) = (/indth,indph/)
                  if (itype==BIQUIN) then
                    igapp=qualify_position_biquin(indph,n1,n2,lfirst_proc_z,llast_proc_z,lfirst_proc_y.or.llast_proc_y)
                  elseif (itype==BICUB.or.itype==BIQUAD) then
                    igapp=qualify_position_bicub(indph,n1,n2,lfirst_proc_z,llast_proc_z,lfirst_proc_y.or.llast_proc_y)
                  else
                    igapp=qualify_position_bilin(indph,n1,n2)
                  endif
!if (.not.lyang.and.(igapt/=0 .or. igapp/=0)) & 
!  print'(a,3(i3,1x),4(f8.3,1x))', 'iproc,inds,gridcoors,pointcoors=', iproc, indth,indph, y(indth), z(indph), thphprime(:,ip,jp)
                  exit
                endif
              enddo
            endif
          endif

          if (okt.and.okp) then
!
!  If both detections positive, determine index ranges for overlap.
!
            if (lyang) then
              indth_in=max(m1,min(m2,indth)); indph_in=max(n1,min(n2,indph))
              ind=indph_in-nghost          
              if (lfirst_proc_y) then
                overlap_thrngs(ind)=max(overlap_thrngs(ind),max(m1,indth_in-1))
              elseif (llast_proc_y) then
                overlap_thrngs(ind)=min(overlap_thrngs(ind),min(m2,indth_in+1))
              endif

              ind=indth_in-nghost          
              if (lfirst_proc_z) then 
                overlap_phrngs(ind)=max(overlap_phrngs(ind),max(n1,indph_in-1))
              elseif (llast_proc_z) then
                overlap_phrngs(ind)=min(overlap_phrngs(ind),min(n2,indph_in+1))
              endif
            endif

            nok=nok+1
            if (present(ngap)) then
              if (igapt/=NOGAP .neqv. igapp/=NOGAP) then
!
! If the position of the point is special w.r.t. exactly one direction, it is
! detected by two processors (finally a division by four is performed by the
! caller!).
!
                ngap=ngap+2
              elseif (igapt/=NOGAP .and. igapp/=NOGAP) then
!
! If the position of the point is special w.r.t. both directions, it would be
! detected by four processors, but this can only happen if in the theta direction
! there are only two processors.
print*, 'DOUBLE GAP: iproc=', iproc
!
                ngap=ngap+1
              endif
            endif
!
!  Calculate interpolation weights.
!
            if (itype==BIQUAD .or. itype==BICUB) then

              if (.not.allocated(indcoeffs%coeffs2)) then
                allocate(indcoeffs%coeffs2(sz1,sz2,4,4))
                indcoeffs%coeffs2=0.
              endif
!if (iproc_world==5) then
!write(24,'(2(i2,1x),2(e13.6,1x))') ip, jp, thphprime(:,ip,jp)
!endif
              if (itype==BIQUAD) then
                call prep_biquad_interp(thphprime(:,ip,jp),indth,indph,indcoeffs,ip,jp,ma,me,na,ne)
              else
                call prep_bicub_interp(thphprime(:,ip,jp),indth,indph,indcoeffs,ip,jp,igapt,igapp)
              endif

            elseif (itype==BIQUIN) then

              if (.not.allocated(indcoeffs%coeffs2)) then
                allocate(indcoeffs%coeffs2(sz1,sz2,6,6))
                indcoeffs%coeffs2=0.
              endif
              call prep_biquint_interp(thphprime(:,ip,jp),indth,indph,indcoeffs,ip,jp,igapt,igapp)

            elseif (itype==QUADSPLINE) then

              if (.not.allocated(indcoeffs%coeffs2)) then
                allocate(indcoeffs%coeffs2(sz1,sz2,3,3))
                indcoeffs%coeffs2=0.
              endif

              call prep_quadspline_interp(thphprime(:,ip,jp),indth,indph,indcoeffs,ip,jp,ma,me,na,ne)

            else      ! bilinear interpolation

              if (.not.allocated(indcoeffs%coeffs)) then
                allocate(indcoeffs%coeffs(sz1,sz2,4))
                indcoeffs%coeffs=0.
              endif

              dthp=thphprime(1,ip,jp)-y(indth-1)
              dphp=thphprime(2,ip,jp)-z(indph-1)

              qth2 = dthp/dy; qth1 = 1.-qth2
              qph2 = dphp/dz; qph1 = 1.-qph2

              indcoeffs%coeffs(ip,jp,:) = (/qth1*qph1,qth1*qph2,qth2*qph1,qth2*qph2/)

!if (igapt/=0 .and. igapp/=0) print*, 'ip,jp, igapt, igapp, thphprime=', ip,jp, igapt, igapp, thphprime(:,ip,jp)
              if (igapt==LGAP) then
                indcoeffs%coeffs(ip,jp,1:2) = 0.
              elseif (igapt==RGAP) then
                indcoeffs%coeffs(ip,jp,3:4) = 0.
              endif
              if (igapp==LGAP) then
                indcoeffs%coeffs(ip,jp,(/1,3/)) = 0.
              elseif (igapp==RGAP) then
                indcoeffs%coeffs(ip,jp,(/2,4/)) = 0.
              endif
              
            endif
          else
            !!if (iproc_world==1) write(iproc_world+100,*) thphprime(:,ip,jp), 9, 9
          endif

        enddo
      enddo

      if (.false.) then
      !if (.not.lyang.and.(nasym_th/=0 .or. nasym_ph/=0)) then
        txt='Warning: in processor '//trim(itoa(iproc))//','
        if (nasym_th/=0) txt=trim(txt)//' '//trim(itoa(nasym_th))
        if (nasym_th/=0.and.nasym_ph/=0) txt=trim(txt)//' and'
        if (nasym_ph/=0) txt=trim(txt)//' '//trim(itoa(nasym_ph))
        txt=trim(txt)//' points subject to asymmetric interpolation in'
        if (nasym_th/=0) txt=trim(txt)//' theta'
        if (nasym_th/=0.and.nasym_ph/=0) txt=trim(txt)//' and'
        if (nasym_ph/=0) txt=trim(txt)//' phi'
        txt=trim(txt)//' direction!'
        print*, trim(txt)
        print*, '        Consider increasing overlap of grids by making dang (more) negative in start.f90!'
      endif

      if (present(th_range)) then
        if (nok>0) then
! 
!  If any points caught, go through all sz2 columns of thphprime and determine
!  for each the range of rows in which points were caught. th_range is the union of all
!  these ranges.
!
          i2=0; i1=sz1+1
          do ip=1,sz2
            th_range=find_index_range_hill(indcoeffs%inds(:,ip,1),(/0,0/))
            if (th_range(1)<i1) i1=th_range(1)
            if (th_range(2)>i2) i2=th_range(2)
          enddo
          th_range=(/i1,i2/)
        else
          th_range=(/0,0/)
        endif
      endif

    endfunction prep_interp
!**************************************************************************
    function in_overlap_mask(indth,indph) result (ok)
!
!  Determines whether point (indth,indph) of Yang grid lies in overlap region.
!
!  12-dec-17/MR: coded
!
      use Cdata, only: lfirst_proc_y, lfirst_proc_z, llast_proc_y, llast_proc_z, nghost, iproc

      integer, intent(IN) :: indth,indph
      logical :: ok

      integer :: ind
      logical, save :: l0=.true.
   
      ok=.true.
      if (lyang) then 
!if (l0.and.iproc==1.and.allocated(overlap_thrngs)) print*, 'overlap_thrngs=', overlap_thrngs
!if (l0.and.iproc==1.and.allocated(overlap_phrngs)) print*, 'overlap_phrngs=', overlap_phrngs
if (l0) l0=.false.

        if (lfirst_proc_y.or.llast_proc_y) then
          ind=indph-nghost
          if (overlap_thrngs(ind)/=0) then
            if (lfirst_proc_y) then
              ok = indth>overlap_thrngs(ind)
            else
              ok = indth<overlap_thrngs(ind)
            endif
          endif
        endif

        if (lfirst_proc_z.or.llast_proc_z) then
          ind=indth-nghost
          if (overlap_phrngs(ind)/=0) then
            if (lfirst_proc_z) then
              ok = indph>overlap_phrngs(ind)
            else
              ok = indth<overlap_phrngs(ind)
            endif
          endif
        endif
      endif
 
    endfunction in_overlap_mask
!**************************************************************************
    subroutine coeffs_to_weights(intcoeffs,indweights)
!
!  Transforms interpolation data w.r.t to a yz strip of (non-grid) points to
!  a representation in which for each (m,n) in the local range (m1,m2)x(n1,n2)
!  it is stored to which phi-coordinate line of the strip (its theta index)
!  a grid point with (m,n) contributes and with which weight. As the use is
!  within z sums (or averages) the contributions are summed up. 
!
!   20-mar-16/MR: coded
!
      use General, only: pos_in_array, itoa
      use Cparam , only: m1,m2,n1,n2
      use Cdata, only: lroot
      !use Messages, only: warning

      type(ind_coeffs), intent(IN) :: intcoeffs
      type(ind_coeffs), intent(OUT):: indweights

      integer :: mm, nn, i, j, indth, indph, dindth, dindph, sz1, sz2, thpos, ith, ithmax
      real :: coeff
      integer, parameter :: nlines_max=nygrid   ! perhaps too small.
      logical :: ltoo_many_lines

       if (allocated(indweights%inds)) deallocate(indweights%inds)
       if (allocated(indweights%coeffs)) deallocate(indweights%coeffs)
       allocate(indweights%inds(m1:m2,n1:n2,nlines_max), &
                indweights%coeffs(m1:m2,n1:n2,nlines_max))

       sz1=size(intcoeffs%inds,1); sz2=size(intcoeffs%inds,2)
       indweights%inds=0

       ltoo_many_lines=.false.; ithmax=0
!
!  Iterate over all (non-ghost) points (mm,nn) in the y-z plane of the local grid.
!
       do nn=n1,n2
         do mm=m1,m2
           ith=0
           do i=1,sz1
             do j=1,sz2

               indth=intcoeffs%inds(i,j,1); indph=intcoeffs%inds(i,j,2)
               if (indth<-1) cycle
!
!  Iterate over all points (i,j) of the strip for which interpolation is (potentially) performed as
!  they lie within the local grid. (indth,indph) refers to the grid point which
!  is employed in interpolation.
!
               dindth=indth-mm; dindph=indph-nn

               if ((dindth==1.or.dindth==0) .and. (dindph==1.or.dindph==0)) then
!
!  Grid point (mm,nn) coincides with one of the four (bilinear interpolation!)
!  grid points employed in interpolation and hence contributes to coordinate line i of strip.
!
                 if (ith==0) then
                   thpos=0
                 else
!
!  Has this point already contributed to a line? Then determine thpos>0.
!
                   thpos=pos_in_array(indth,indweights%inds(mm,nn,:ith))
                 endif

                 if (thpos==0) then
!
!  Not yet contributed -> check for overflow over nlines_max.
!
                   if (ith==nlines_max) then
                     if (.not.ltoo_many_lines) then
                       !call warning(
                       if (lroot) print*, 'coeffs_to_weights: WARNING -- ', &
                                          'More than '//trim(itoa(nlines_max))// &
                                          ' lines contributed to by point at (m,n)=(', &
                                          trim(itoa(mm)),',',trim(itoa(nn)),') in proc '//trim(itoa(iproc_world))
                       ltoo_many_lines=.true.
                     endif
                   else
                     ith=ith+1
                   endif
!
!  If point (mm,nn) contributes to too many lines, ignore this and any further
!  contributions to additional lines. (problematic!)
!
                   if (ltoo_many_lines) cycle
                   thpos=ith
!
!  Store y-index i of line for point (mm,nn), initialize weight.
!
                   indweights%inds(mm,nn,ith)=i
                   indweights%coeffs(mm,nn,ith)=0.
                 endif
!
!  Detect with which "corner" of interpolation cell point (mm,nn) coincides and retrieve interpolation
!  coefficient = weight.
!              
                 if (dindth==1) then
                   coeff=intcoeffs%coeffs(i,j,2-dindph)
                 else
                   coeff=intcoeffs%coeffs(i,j,4-dindph)
                 endif
!
!  Accumulate weights.
!
                 indweights%coeffs(mm,nn,thpos)=indweights%coeffs(mm,nn,thpos)+coeff

               endif
             enddo
           enddo
           ithmax=max(ith,ithmax)    ! not needed
         enddo
       enddo
!print*, 'iproc,ithmax=', iproc_world,ithmax
    endsubroutine coeffs_to_weights
!*******************************************************************
end module
