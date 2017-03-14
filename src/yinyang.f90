! $Id$
!
! MODULE_DOC: This module contains Yin-Yang related types and functions
! MODULE_DOC: which are incompatible with FORTRAN 95.
!
!***************************************************************
!
module Yinyang
!
  use Cdata, only: iproc_world

  implicit none

  include 'yinyang.h'

  type ind_coeffs
    integer, dimension(:,:,:), allocatable :: inds
    real, dimension(:,:,:), allocatable :: coeffs
    real, dimension(:,:,:,:), allocatable :: coeffs2
  end type
!
  contains

!**************************************************************************
   subroutine biquad_interp(indcoeffs, ith, iph, f, buffer, i2buf, i3buf)
!
!  Performs biquadratic interpolation for a pencil at position (ith, iph) of the
!  original strip from values in f-array using the precalculated weights
!  in indcoeffs%coeffs. Result is returned in buffer(i2buf,i3buf)=(ith,iph)
!  or buffer(i2buf,i3buf)=(iph,ith), the latter if transposition is required.
!
!  20-dec-15/MR: coded
! 
      use Cdata, only: y,z
      use General, only: transform_spher_cart_yy, notanumber

      type(ind_coeffs),         intent(IN) :: indcoeffs
      integer,                  intent(IN) :: ith, iph, i2buf, i3buf
      real, dimension(:,:,:,:), intent(IN) :: f
      real, dimension(:,:,:,:), intent(OUT):: buffer

      integer :: indthl, indthu, indphl, indphu, igp, igt, nphrng, nthrng
      real, dimension(:,:,:,:), allocatable :: tmp

      indthl=indcoeffs%inds(ith,iph,1); indthu=indcoeffs%inds(ith,iph,2)
      if (indthl==0) return

      indphl=indcoeffs%inds(ith,iph,3); indphu=indcoeffs%inds(ith,iph,4)
      nthrng=indthu-indthl+1; nphrng=indphu-indphl+1

      if (size(buffer,4)==3) then
!
!  Vector field, is transformed to Cartesian basis of *other* grid before being interpolated 
!
        allocate (tmp(size(buffer,1),nthrng,nphrng,3))
        call transform_spher_cart_yy(f,indthl,indthu,indphl,indphu,tmp,lyy=.true.)

        do igp=1,nphrng; do igt=1,nthrng
          buffer(:,i2buf,i3buf,:) = buffer(:,i2buf,i3buf,:)+indcoeffs%coeffs2(ith,iph,igt,igp)*tmp(:,igt,igp,:)
        enddo; enddo

if (notanumber(buffer(:,i2buf,i3buf,:))) print*, 'indthl,indph, i2buf,i3buf=', indthl,indphl, i2buf,i3buf
      else
!
!  Scalar field - no transformation.
!
        do igp=1,nphrng; do igt=1,nthrng
          buffer(:,i2buf,i3buf,1) = buffer(:,i2buf,i3buf,1) &
                                   +indcoeffs%coeffs2(ith,iph,igt,igp)*f(:,indthl+igt-1,indphl+igp-1,1)
        enddo; enddo
if (iproc_world==5) then
!write(23,'(4(i2,1x),30(e13.6,1x))') ith,iph,nthrng,nphrng,y(indthl:indthu),z(indphl:indphu), f(l1+8,indthl:indthu,indphl:indphu,1), buffer(l1+8,i2buf,i3buf,1)
endif
if (notanumber(buffer(:,i2buf,i3buf,1))) print*, 'NaNs at', iproc_world,  &
  ',indthl,indph, i2buf,i3buf=', indthl,indphl,ith,iph,i2buf,i3buf
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

      type(ind_coeffs),         intent(IN) :: indcoeffs
      integer,                  intent(IN) :: ith, iph, i2buf, i3buf
      real, dimension(:,:,:,:), intent(IN) :: f
      real, dimension(:,:,:,:), intent(OUT):: buffer

      integer :: indth, indph
      real, dimension(:,:,:,:), allocatable :: tmp

      indth=indcoeffs%inds(ith,iph,1)
      indph=indcoeffs%inds(ith,iph,2)

      if (indth==0 .or. indph==0) return

      if (size(buffer,4)==3) then
!
!  Vector field, is transformed to Cartesian basis of *other* grid before being interpolated 
!
        allocate (tmp(size(buffer,1),2,2,3))
        call transform_spher_cart_yy(f,indth-1,indth,indph-1,indph,tmp,lyy=.true.)

        buffer(:,i2buf,i3buf,:) = indcoeffs%coeffs(ith,iph,1)*tmp(:,1,1,:) &
                                 +indcoeffs%coeffs(ith,iph,2)*tmp(:,1,2,:) &
                                 +indcoeffs%coeffs(ith,iph,3)*tmp(:,2,1,:) &
                                 +indcoeffs%coeffs(ith,iph,4)*tmp(:,2,2,:)
if (notanumber(buffer(:,i2buf,i3buf,:))) print*, 'indth,indph, i2buf,i3buf=', indth,indph, i2buf,i3buf
      else
!
!  Scalar field - no transformation.
!
        buffer(:,i2buf,i3buf,1) = indcoeffs%coeffs(ith,iph,1)*f(:,indth-1,indph-1,1) &
                                 +indcoeffs%coeffs(ith,iph,2)*f(:,indth-1,indph  ,1) &
                                 +indcoeffs%coeffs(ith,iph,3)*f(:,indth  ,indph-1,1) &
                                 +indcoeffs%coeffs(ith,iph,4)*f(:,indth  ,indph  ,1)
if (notanumber(buffer(:,i2buf,i3buf,1))) print*, 'NaNs at', iproc_world,  &
  ',indth,indph, i2buf,i3buf=', indth,indph,ith,iph,i2buf,i3buf
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
        w1ph=-qph(3); w2ph=qph(2)
        indcoeffs%coeffs2(ip,jp,1:3,1:3)=w1ph*indcoeffs%coeffs2(ip,jp,1:3,1:3)

        do igp=1,3; do igt=1,3
          indcoeffs%coeffs2(ip,jp,igt,igp+1)=indcoeffs%coeffs2(ip,jp,igt,igp+1)+w2ph*gth(igt)*gph4(igp)
        enddo; enddo
!print*, 'indcoeffs34=', sum(indcoeffs%coeffs2(ip,jp,:,:))

      endif

      if (nthrng==4) then

        gth=(/.5*qth(3)*qth(4),-qth(2)*qth(4),.5*qth(2)*qth(3)/)
!print*, 'gth4=', sum(gth)
        w1=-qth(3); w2th=qth(2)
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
    function prep_interp(thphprime,indcoeffs,itype,th_range) result (nok)
!
!  For each of the points in the strip thphprime (with shape 2 x thprime-extent x
!  phprime-extent), arbitrarily positioned in the yz-plane, determine in
!  which cell of the grid y(ma:me) x z(na:ne) it lies, store indices of the
!  cells upper right corner in indcoeffs%inds and the weights of bilinear
!  interpolation for the four corners in indcoeffs%coeffs. If no cell is found
!  for a point, indcoeffs%inds and indcoeffs%coeffs are set zero.
!  If present, return in th_range the interval in thprime-extent in which
!  interpolation cells could be assigned to any points.
!  Returns number of points in thphprime for which interpolation cell could be
!  found.
!
!  20-dec-15/MR: coded
!
      use General, only: find_index_range_hill
      use Cdata,   only: y, z, lfirst_proc_y, ipy, lfirst_proc_z, ipz, lroot
      use Cparam,  only: m1,m2,n1,n2, BILIN, BIQUAD

      real, dimension(:,:,:),          intent(IN) :: thphprime
      type(ind_coeffs),                intent(OUT):: indcoeffs
      integer,                         intent(IN) :: itype
      integer, dimension(2), optional, intent(OUT):: th_range

      integer :: nok

      integer :: ip, jp, indth, indph, sz1, sz2, ma, na, me, ne, i1, i2
      real :: dth, dph, dthp, dphp, qth1, qth2, qph1, qph2
      logical :: okt, okp

      sz1=size(thphprime,2); sz2=size(thphprime,3)

      if (allocated(indcoeffs%inds)) deallocate(indcoeffs%inds,indcoeffs%coeffs)
      if (itype==BIQUAD) then
        allocate(indcoeffs%inds(sz1,sz2,4))
      elseif (itype==BILIN) then
        allocate(indcoeffs%inds(sz1,sz2,2))
      else
        if (lroot) print*, 'prep_interp: Only bilinear and biquadratic interpolations implemented'
        stop
      endif
      indcoeffs%inds=0

      nok=0
      ma=m1-1; me=m2
!
!  In order to catch points which lie in gaps between the domains of procs,
!  the first lower ghost cell in y direction is included, except for the first
!  proc where it is the first upper, the next proc (ipy=1) then doesn't need a
!  ghost cell.
!
      if (lfirst_proc_y) then
        ma=m1; me=m2+1
      elseif (ipy==1) then
        ma=m1
      endif
!
!  Likewise for z direction.
!
      na=n1-1; ne=n2
      if (lfirst_proc_z) then
        na=n1; ne=n2+1
      elseif (ipz==1) then
        na=n1
      endif

      do ip=1,sz1
        do jp=1,sz2
!
!  For all points in strip,
!
          okt=.false.; okp=.false.
          if ( y(ma)<=thphprime(1,ip,jp) ) then
!
!  detect between which y lines it lies.
!
            do indth=ma+1,me
              if ( y(indth)>=thphprime(1,ip,jp) ) then
                indcoeffs%inds(ip,jp,1) = indth
                okt=.true.
                exit
              endif
            enddo
          endif

          if (okt) then
            if ( z(na)<=thphprime(2,ip,jp) ) then
!
!  detect between which z lines it lies.
!
              do indph=na+1,ne
                if ( z(indph)>=thphprime(2,ip,jp) ) then
                  indcoeffs%inds(ip,jp,2) = indph
                  okp=.true.
                  exit
                endif
              enddo
            endif
            if (.not.okp) indcoeffs%inds(ip,jp,1)=0
          endif

          if (okt.and.okp) then
!
!  If both detections positive, calculate interpolation weights.
!
            nok=nok+1

            if (itype==BIQUAD) then

              if (.not.allocated(indcoeffs%coeffs2)) then
                allocate(indcoeffs%coeffs2(sz1,sz2,4,4))
                indcoeffs%coeffs2=0.
              endif
!if (iproc_world==5) then
!write(24,'(2(i2,1x),2(e13.6,1x))') ip, jp, thphprime(:,ip,jp)
!endif
              call prep_biquad_interp(thphprime(:,ip,jp),indth,indph,indcoeffs,ip,jp,ma,me,na,ne)

            else

              if (.not.allocated(indcoeffs%coeffs)) then
                allocate(indcoeffs%coeffs(sz1,sz2,4))
                indcoeffs%coeffs=0.
              endif

              dth=y(indth)-y(indth-1); dph=z(indph)-z(indph-1)    ! dy, dz !
              dthp=thphprime(1,ip,jp)-y(indth-1)
              dphp=thphprime(2,ip,jp)-z(indph-1)

              qth2 = dthp/dth; qth1 = 1.-qth2
              qph2 = dphp/dph; qph1 = 1.-qph2

              indcoeffs%coeffs(ip,jp,:) = (/qth1*qph1,qth1*qph2,qth2*qph1,qth2*qph2/)
!if (iproc_world==0) &
  !print*, 'iproc_world, coeffs=', iproc_world, indcoeffs%inds(ip,jp,:)
            endif
          endif

        enddo
      enddo

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
      !use Messages, only: warning

      type(ind_coeffs) :: intcoeffs,indweights

      integer :: mm, nn, i, j, indth, indph, dindth, dindph, sz1, sz2, thpos, ith, ithmax
      real :: coeff
      integer, parameter :: nlines_max=128   ! perhaps too small.
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
               if (indth==0) cycle
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
                       print*, 'coeffs_to_weights: WARNING -- ', &
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
                   if (dindph==1) then
                     coeff=intcoeffs%coeffs(i,j,1)
                   else
                     coeff=intcoeffs%coeffs(i,j,2)
                   endif
                 else
                   if (dindph==1) then
                     coeff=intcoeffs%coeffs(i,j,3)
                   else
                     coeff=intcoeffs%coeffs(i,j,4)
                   endif
                 endif
!
!  Accumulate weights.
!
                 indweights%coeffs(mm,nn,thpos)=indweights%coeffs(mm,nn,thpos)+coeff

               endif
             enddo
           enddo
           ithmax=max(ith,ithmax)
         enddo
       enddo
!print*, 'iproc,ithmax=', iproc_world,ithmax

    endsubroutine coeffs_to_weights
!*******************************************************************
end module
