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

  include 'yinyang.h'

  type ind_coeffs
    integer, dimension(:,:,:), allocatable :: inds
    real, dimension(:,:,:), allocatable :: coeffs
  end type

  contains

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
    function prep_bilin_interp(thphprime,indcoeffs,th_range) result (nok)
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
      use Cdata,   only: y, z, lfirst_proc_y, ipy, lfirst_proc_z, ipz
      use Cparam,  only: m1,m2,n1,n2

      real, dimension(:,:,:),          intent(IN) :: thphprime
      type(ind_coeffs),                intent(OUT):: indcoeffs
      integer, dimension(2), optional, intent(OUT):: th_range

      integer :: nok

      integer :: ip, jp, indth, indph, sz1, sz2, ma, na, me, ne, i1, i2
      real :: dth, dph, dthp, dphp, qth1, qth2, qph1, qph2
      logical :: okt, okp

      sz1=size(thphprime,2); sz2=size(thphprime,3)

      if (allocated(indcoeffs%inds)) deallocate(indcoeffs%inds,indcoeffs%coeffs)
      allocate(indcoeffs%inds(sz1,sz2,2))
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

            if (.not.allocated(indcoeffs%coeffs)) then
              allocate(indcoeffs%coeffs(sz1,sz2,4))
              indcoeffs%coeffs=0.
            endif
!
!  If both detections positive, calculate interpolation weights.
!
            nok=nok+1

            dth=y(indth)-y(indth-1); dph=z(indph)-z(indph-1)
            dthp=thphprime(1,ip,jp)-y(indth-1)
            dphp=thphprime(2,ip,jp)-z(indph-1)

            qth2 = dthp/dth; qth1 = 1.-qth2
            qph2 = dphp/dph; qph1 = 1.-qph2

            indcoeffs%coeffs(ip,jp,:) = (/qth1*qph1,qth1*qph2,qth2*qph1,qth2*qph2/)
!if (iproc_world==0) &
  !print*, 'iproc_world, coeffs=', iproc_world, indcoeffs%inds(ip,jp,:)
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

    endfunction prep_bilin_interp
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
!  Has this point already contributed to a line? Then thpos>0.
!
                   thpos=pos_in_array(indth,indweights%inds(mm,nn,:ith))
                 endif

                 if (thpos==0) then
!
!  Not yet contributed -> check for overflow over nlines_max.
!
                   if (ith==nlines_max.and..not.ltoo_many_lines) then
                     !call warning(
                     print*, 'coeffs_to_weights: WARNING -- ', &
                          'More than '//itoa(nlines_max)// &
                          ' points contributed to in proc '//itoa(iproc_world)
                     ltoo_many_lines=.true.
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
