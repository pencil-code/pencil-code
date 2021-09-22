! $Id$
!
!  This module provides routines for obtaining slice data.
!
module Slices_methods
!
  use Cdata
!
  implicit none
!
  interface assign_slices_scal
    module procedure assign_slices_f_scal
    module procedure assign_slices_sep_scal
  endinterface

  interface assign_slices_vec
    module procedure assign_slices_f_vec
    module procedure assign_slices_sep_vec
  endinterface

  interface store_slices
    module procedure store_slices_scal
    module procedure store_slices_vec
  endinterface

  interface process_slices
    module procedure process_slices_func
    module procedure process_slices_fac
  endinterface

  interface alloc_rslice
    module procedure alloc_rslice_scal
    module procedure alloc_rslice_scal_2D
    module procedure alloc_rslice_vec
  endinterface

  public :: store_slices
  public :: process_slices, addto_slices
  public :: assign_slices_vec
  public :: assign_slices_scal, assign_slices_f_scal
  !public :: alloc_slice_buffers
  public :: nullify_slice_pointers
  public :: write_rslice_position
  public :: alloc_rslice, prep_rslice
  public :: exp2d, log2d, abs2d

!  interface alloc_slice_buffers
!    module procedure alloc_slice_buffers_scal
!    module procedure alloc_slice_buffers_vec
!  endinterface
!
private
!
!  For spherical slices in Cartesian geometry:
!
  real, dimension(:), allocatable :: cph_slice,sph_slice,cth_slice,sth_slice
  integer :: ith_min,ith_max,iph_min,iph_max
  real, dimension(:,:,:,:,:), allocatable :: rslice_interp_weights
  integer, dimension(:,:,:),allocatable :: rslice_adjec_corn_inds

  contains
!***********************************************************************
!  include 'slice_methods.f90.removed'
!***********************************************************************
    subroutine assign_slices_sep_scal(slices,xy,xz,yz,xy2,xy3,xy4,xz2,r)
!
!  Assignment of slice pointers according to slice selection switches.
!
!  12-apr-17/MR: coded 
!   8-oct-18/MR: made some pars optional
!
      type (slice_data),                      intent(OUT):: slices
      real, dimension(:,:), target,           intent(IN) :: xy
      real, dimension(:,:), target, optional, intent(IN) :: xz,yz,xy2,xy3,xy4,xz2
      real, dimension(ith_min:ith_max,iph_min:iph_max,-1:0,-1:0,-1:0), target, optional, intent(IN) :: r

      if (lwrite_slice_yz ) slices%yz =>yz
      if (lwrite_slice_xz ) slices%xz =>xz
      if (lwrite_slice_xy ) slices%xy =>xy
      if (lwrite_slice_xy2) slices%xy2=>xy2
      if (lwrite_slice_xy3) slices%xy3=>xy3
      if (lwrite_slice_xy4) slices%xy4=>xy4
      if (lwrite_slice_xz2) slices%xz2=>xz2
      if (present(r).and.lwrite_slice_r) call interp_rslice(slices%r,rslice=r)

      slices%ready=.true.

    endsubroutine assign_slices_sep_scal
!***********************************************************************
    subroutine assign_slices_sep_vec(slices,xy,xz,yz,xy2,xy3,xy4,xz2,r,ncomp)
!
!  Assignment of slice pointers according to slice selection switches.
!
!  12-apr-17/MR: coded 
!
    use General, only: ioptest
!
    type (slice_data)             , intent(INOUT):: slices
    real, dimension(:,:,:), target, intent(IN)   :: xy,xz,yz,xy2,xy3,xy4,xz2
    real, dimension(ith_min:ith_max,iph_min:iph_max,-1:0,-1:0,-1:0,3), target, optional, intent(IN)   :: r
    integer, optional             , intent(IN)   :: ncomp

    integer :: nc

    nc=ioptest(ncomp,3)

    if (slices%index>=nc) then
      slices%ready=.false.
    else
      slices%index=slices%index+1

      if (lwrite_slice_yz ) slices%yz => yz(:,:,slices%index)
      if (lwrite_slice_xz ) slices%xz => xz(:,:,slices%index)
      if (lwrite_slice_xy ) slices%xy => xy(:,:,slices%index)
      if (lwrite_slice_xy2) slices%xy2=>xy2(:,:,slices%index)
      if (lwrite_slice_xy3) slices%xy3=>xy3(:,:,slices%index)
      if (lwrite_slice_xy4) slices%xy4=>xy4(:,:,slices%index)
      if (lwrite_slice_xz2) slices%xz2=>xz2(:,:,slices%index)
      if (present(r).and.lwrite_slice_r) &
        call interp_rslice(slices%r,rslice=r(:,:,:,:,:,slices%index))
      if (slices%index<=nc) slices%ready=.true.
    endif

    endsubroutine assign_slices_sep_vec
!***********************************************************************
    subroutine assign_slices_f_scal(slices,f,ind1,ind2)
!
!  Copying of scalar data from f-array to arrays assigned to the slice pointers 
!  according to slice selection switches.
!
!  12-apr-17/MR: coded 
!  28-may-19/MR: added optional argument ind2
!
      use General, only: ioptest

      type (slice_data)          , intent(OUT):: slices
      real, dimension(mx,my,mz,*), intent(IN) :: f
      integer                    , intent(IN) :: ind1
      integer         , optional , intent(IN) :: ind2

      integer :: ind2_

      ind2_=ioptest(ind2,ind1)
      if (lwrite_slice_yz)  slices%yz =f(ix_loc,m1:m2  ,n1:n2  ,ind1) !:ind2_)
      if (lwrite_slice_xz)  slices%xz =f(l1:l2 ,iy_loc ,n1:n2  ,ind1) !:ind2_)
      if (lwrite_slice_xy)  slices%xy =f(l1:l2 ,m1:m2  ,iz_loc ,ind1) !:ind2_)
      if (lwrite_slice_xy2) slices%xy2=f(l1:l2 ,m1:m2  ,iz2_loc,ind1) !:ind2_) 
      if (lwrite_slice_xy3) slices%xy3=f(l1:l2 ,m1:m2  ,iz3_loc,ind1) !:ind2_)
      if (lwrite_slice_xy4) slices%xy4=f(l1:l2 ,m1:m2  ,iz4_loc,ind1) !:ind2_)
      if (lwrite_slice_xz2) slices%xz2=f(l1:l2 ,iy2_loc,n1:n2  ,ind1) !:ind2_) 
      if (lwrite_slice_r  ) call interp_rslice(slices%r,f,ind1)

      slices%ready=.true.

    endsubroutine assign_slices_f_scal
!***********************************************************************
    subroutine assign_slices_f_vec(slices,f,ind,ncomp)
!
!  Copying of multi-component data from f-array to arrays assigned to the slice pointers 
!  according to slice selection switches. If ncomp is not present the quantity
!  in f(:,:,:,ind:ind+2) is considered a proper vector (to be transformed in
!  Yin-Yang context), otherwise a general multi-component quantity (not to be
!  transformed).
!
!  12-apr-17/MR: coded 
!
    use General, only: transform_thph_yy_other, ioptest

    type (slice_data)       , intent(INOUT):: slices
    real, dimension(:,:,:,:), intent(IN)   :: f
    integer                 , intent(IN)   :: ind
    integer, optional       , intent(IN)   :: ncomp

    real, dimension(:,:,:,:), allocatable, save :: transformed
    integer :: nc

    nc=ioptest(ncomp,3)

    if (slices%index>=nc) then
      slices%ready=.false.
    else
      slices%index=slices%index+1

      if (lwrite_slice_yz) then
        if (lyang.and..not.present(ncomp).and.slices%index>=2) then
!
!  On Yang grid: transform theta and phi components of vector to Yin-grid basis.
!  (phi component is saved in transformed for use in next call.
!
          if (slices%index==2) then
            if (.not.allocated(transformed)) allocate(transformed(1,ny,nz,3))
            call transform_thph_yy_other(f(ix_loc:ix_loc,:,:,ind+1:ind+2),m1,m2,n1,n2,transformed)
          endif
!
!  theta component is used immediately, phi component with next call.
!
          slices%yz=transformed(1,:,:,slices%index-1)
        else
          slices%yz=f(ix_loc,m1:m2,n1:n2,ind-1+slices%index)
        endif

      endif

      if (lwrite_slice_xz) &
        slices%xz =f(l1:l2,iy_loc ,n1:n2  ,ind-1+slices%index)
      if (lwrite_slice_xy) &
        slices%xy =f(l1:l2,m1:m2  ,iz_loc ,ind-1+slices%index)
      if (lwrite_slice_xy2) &
        slices%xy2=f(l1:l2,m1:m2  ,iz2_loc,ind-1+slices%index)
      if (lwrite_slice_xy3) &
        slices%xy3=f(l1:l2,m1:m2  ,iz3_loc,ind-1+slices%index)
      if (lwrite_slice_xy4) &
        slices%xy4=f(l1:l2,m1:m2  ,iz4_loc,ind-1+slices%index)
      if (lwrite_slice_xz2) &
        slices%xz2=f(l1:l2,iy2_loc,n1:n2  ,ind-1+slices%index)
      if (lwrite_slice_r  ) call interp_rslice(slices%r,f,ind-1+slices%index)

      if (slices%index<=nc) slices%ready=.true.
    endif

    endsubroutine assign_slices_f_vec
!***********************************************************************
    subroutine store_slices_scal(pencil,xy,xz,yz,xy2,xy3,xy4,xz2,r)
!
!  Stores scalar data from pencil in auxiliary arrays when needed.
!
!  12-apr-17/MR: coded 
!   8-oct-18/MR: made some pars optional
!
      real, dimension(nx) , intent(IN)           :: pencil
      real, dimension(:,:), intent(OUT)          :: xy
      real, dimension(:,:), intent(OUT),optional :: xz,yz,xy2,xy3,xy4,xz2
      real, dimension(ith_min:ith_max,iph_min:iph_max,-1:0,-1:0,-1:0), intent(OUT), optional :: r

      if (n==iz_loc) xy(:,m-m1+1)=pencil
      if (present(xz) .and.m==iy_loc)  xz(:,n-n1+1)=pencil
      if (present(yz) .and.lwrite_slice_yz) yz(m-m1+1,n-n1+1)=pencil(ix_loc-l1+1)
      if (present(xy2).and.n==iz2_loc) xy2(:,m-m1+1)=pencil
      if (present(xy3).and.n==iz3_loc) xy3(:,m-m1+1)=pencil
      if (present(xy4).and.n==iz4_loc) xy4(:,m-m1+1)=pencil
      if (present(xz2).and.m==iy2_loc) xz2(:,n-n1+1)=pencil
      if (present(r)  .and.lwrite_slice_r) call store_rslice_scal(pencil,r)

    endsubroutine store_slices_scal
!***********************************************************************
    subroutine store_slices_vec(pencil,xy,xz,yz,xy2,xy3,xy4,xz2,r,ncomp)
!
!  Stores multi-component data from pencil in auxiliary arrays when needed.
!
!  12-apr-17/MR: coded 
!
      use General, only: transform_thph_yy, ioptest

      real, dimension(:,:)  , intent(IN) :: pencil
      real, dimension(:,:,:), intent(OUT):: xy,xz,yz,xy2,xy3,xy4,xz2
      real, dimension(ith_min:ith_max,iph_min:iph_max,-1:0,-1:0,-1:0,3), intent(OUT), optional :: r
      integer, optional     , intent(IN) :: ncomp

      real, dimension(1,3) :: tmp

      !nc=ioptest(ncomp,3)
!
!  Note: i[xyz][234]*_loc is the index with respect to array with ghost zones.
!
      if (lwrite_slice_yz) then
        if (lyang.and..not.present(ncomp).and.size(pencil,2)==3) then
!
!  On Yang grid: transform theta and phi components of vector to Yin-grid basis.
!
          call transform_thph_yy(pencil(ix_loc-l1+1:ix_loc-l1+1,:),(/1,1,1/),tmp)
          yz(m-m1+1,n-n1+1,:)=tmp(1,:)
        else
          yz(m-m1+1,n-n1+1,:)=pencil(ix_loc-l1+1,:)
        endif
      endif
!
      if (m==iy_loc ) xz (:,n-n1+1,:)=pencil
      if (n==iz_loc ) xy (:,m-m1+1,:)=pencil
      if (n==iz2_loc) xy2(:,m-m1+1,:)=pencil
      if (n==iz3_loc) xy3(:,m-m1+1,:)=pencil
      if (n==iz4_loc) xy4(:,m-m1+1,:)=pencil
      if (m==iy2_loc) xz2(:,n-n1+1,:)=pencil
      if (present(r).and.lwrite_slice_r) call store_rslice_vec(pencil,r)

    endsubroutine store_slices_vec
!***********************************************************************
    function exp2d(arr) result(res)

    real, dimension(:,:) :: arr
    real, dimension(size(arr,1),size(arr,2)) :: res

    res=exp(arr)
    
    endfunction exp2d
!***********************************************************************
    function log2d(arr) result(res)

    real, dimension(:,:) :: arr
    real, dimension(size(arr,1),size(arr,2)) :: res

    res=alog(arr)

    endfunction log2d
!***********************************************************************
    function abs2d(arr) result(res)

    real, dimension(:,:) :: arr
    real, dimension(size(arr,1),size(arr,2)) :: res

    res=abs(arr)

    endfunction abs2d
!***********************************************************************
    subroutine process_slices_func(slices,func)

      type(slice_data), intent(INOUT):: slices

      interface 
        function func(arr) result(res)
        real, dimension(:,:) :: arr
        real, dimension(size(arr,1),size(arr,2)) :: res
        end function func
      endinterface

      if (lwrite_slice_yz ) slices%yz =func(slices%yz)
      if (lwrite_slice_xz ) slices%xz =func(slices%xz)
      if (lwrite_slice_xy ) slices%xy =func(slices%xy)
      if (lwrite_slice_xy2) slices%xy2=func(slices%xy2)
      if (lwrite_slice_xy3) slices%xy3=func(slices%xy3)
      if (lwrite_slice_xy4) slices%xy4=func(slices%xy4)
      if (lwrite_slice_xz2) slices%xz2=func(slices%xz2)
      if (lwrite_slice_r  ) slices%r  =func(slices%r)

    endsubroutine process_slices_func
!***********************************************************************
    subroutine process_slices_fac(slices,fac)

      type(slice_data), intent(INOUT):: slices
      real,             intent(IN)   :: fac

      if (lwrite_slice_yz ) slices%yz =fac*slices%yz
      if (lwrite_slice_xz ) slices%xz =fac*slices%xz
      if (lwrite_slice_xy ) slices%xy =fac*slices%xy
      if (lwrite_slice_xy2) slices%xy2=fac*slices%xy2
      if (lwrite_slice_xy3) slices%xy3=fac*slices%xy3
      if (lwrite_slice_xy4) slices%xy4=fac*slices%xy4
      if (lwrite_slice_xz2) slices%xz2=fac*slices%xz2
      if (lwrite_slice_r  ) slices%r  =fac*slices%r

    endsubroutine process_slices_fac
!***********************************************************************
    subroutine addto_slices(slices,pencil)
!
!  Adds pencil data to scalar slices according to slice selection switches.
!
!  12-apr-17/MR: coded 
!
      use Messages, only: not_implemented

      type (slice_data),  intent(INOUT):: slices
      real, dimension(nx),intent(IN)   :: pencil

      integer :: i

      if (lwrite_slice_yz ) slices%yz=slices%yz+pencil(ix_loc-l1+1)
      if (lwrite_slice_xz ) then
        do i=1,nz; slices%xz (:,i)=slices%xz (:,i)+pencil; enddo
      endif
      if (lwrite_slice_xy ) then
        do i=1,ny; slices%xy (:,i)=slices%xy (:,i)+pencil; enddo
      endif
      if (lwrite_slice_xy2) then
        do i=1,ny; slices%xy2(:,i)=slices%xy2(:,i)+pencil; enddo
      endif
      if (lwrite_slice_xy3) then
        do i=1,ny; slices%xy3(:,i)=slices%xy3(:,i)+pencil; enddo
      endif
      if (lwrite_slice_xy4) then
        do i=1,ny; slices%xy4(:,i)=slices%xy4(:,i)+pencil; enddo
      endif
      if (lwrite_slice_xz2) then
        do i=1,nz; slices%xz2(:,i)=slices%xz2(:,i)+pencil; enddo
      endif
      if (lwrite_slice_r) call not_implemented('addto_slices','for r-slices')

    endsubroutine addto_slices
!***********************************************************************
    subroutine nullify_slice_pointers(slices)
!
!  Nullifies all pointers in slices struc.
!
!  12-apr-17/MR: coded 
!
      type (slice_data), intent(OUT):: slices

      nullify(slices%yz)
      nullify(slices%xz)
      nullify(slices%xy)
      nullify(slices%xy2)
      nullify(slices%xy3)
      nullify(slices%xy4)
      nullify(slices%xz2)
      nullify(slices%r)

    endsubroutine nullify_slice_pointers
!***********************************************************************
    subroutine prep_rslice

      integer :: ith,iph,ll,mm,nn
      real :: xs,ys,zs,dth,dphi,dth0,dxs,dys,dzs,wsum

      if (allocated(cph_slice)) deallocate(cph_slice,sph_slice,cth_slice,sth_slice)

      allocate(cph_slice(nph_rslice),sph_slice(nph_rslice))
      allocate(cth_slice(nth_rslice),sth_slice(nth_rslice))

      dth=pi/nth_rslice; dphi=2*pi/nph_rslice; dth0=.5*dth
      do iph=1,nph_rslice
        cph_slice(iph)=cos((iph-1)*dphi); sph_slice(iph)=sin((iph-1)*dphi)
      enddo

      ith_min=0; ith_max=0; iph_min=max_int; iph_max=0
      do ith=1,nth_rslice
        cth_slice(ith)=r_rslice*cos(dth0+(ith-1)*dth)
        sth_slice(ith)=r_rslice*sin(dth0+(ith-1)*dth)
        do iph=1,nph_rslice
          xs=sth_slice(ith)*cph_slice(iph)
          ys=sth_slice(ith)*sph_slice(iph)
          zs=cth_slice(ith)
!if (lroot) print*, 'xyzs=', xs, ys, zs
          if ( x(l1-1)<=xs.and.xs<=x(l2) .and. y(m1-1)<=ys.and.ys<=y(m2) .and. &
               z(n1-1)<=zs.and.zs<=z(n2) ) then
            iph_min=min(iph_min,iph); iph_max=max(iph_max,iph)
            if (ith_min==0) ith_min=ith
            ith_max=ith
          endif
        enddo
      enddo
!print*, 'iproc,ith_min,ith_max,iph_min,iph_max=',
!iproc,ith_min,ith_max,iph_min,iph_max

      lwrite_slice_r = (ith_min/=0.and.iph_max/=0)
      if (lwrite_slice_r) then

        if (allocated(rslice_adjec_corn_inds)) &
          deallocate(rslice_adjec_corn_inds,rslice_interp_weights)

        allocate(rslice_adjec_corn_inds(ith_min:ith_max,iph_min:iph_max,3))
        allocate(rslice_interp_weights(ith_min:ith_max,iph_min:iph_max,2,2,2))
        rslice_adjec_corn_inds=0

        do ith=ith_min,ith_max
          do iph=iph_min,iph_max
            xs=sth_slice(ith)*cph_slice(iph)
            ys=sth_slice(ith)*sph_slice(iph)
            zs=cth_slice(ith)
            do nn=n1,n2
              dzs=z(nn)-zs
              if (dzs>=0.) then
                do mm=m1,m2
                  dys=y(mm)-ys
                  if (dys>=0.) then
                    do ll=l1,l2
                      dxs=x(ll)-xs
                      if (dxs>=0.) then
                        rslice_adjec_corn_inds(ith,iph,:)=(/ll,mm,nn/)
                        rslice_interp_weights(ith,iph,:,:,:)= &
                        reshape((/ (/dxs*dys,(dx-dxs)*dys,dxs*(dy-dys),(dx-dxs)*(dy-dys)/)*dzs, &
                                   (/dxs*dys,(dx-dxs)*dys,dxs*(dy-dys),(dx-dxs)*(dy-dys)/)*(dz-dzs) /), &
                        shape(rslice_interp_weights(ith,iph,:,:,:)))/(dx*dy*dz)
                                                                    ! Requires equidistant grid
!wsum=sum(rslice_interp_weights(ith,iph,:,:,:))
!if (abs(wsum-1.)>1e-14) print*, iproc,ith,iph,wsum
                        goto 100
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
100         continue
          enddo
        enddo
!write(20+iproc,'(5(i3,1x))') ((ith, iph, rslice_adjec_corn_inds(ith,iph,:),ith=ith_min,ith_max),iph=iph_min,iph_max)
!write(20+iproc,'(3(i3,1x))') rslice_adjec_corn_inds
!(ith,iph,:), ith=ith_min,ith_max),iph=iph_min,iph_max)
      endif

    endsubroutine prep_rslice
!***********************************************************************
    subroutine alloc_rslice_scal(slice)

      real, dimension(:,:,:,:,:), allocatable :: slice

      if (allocated(slice)) deallocate(slice)
      allocate(slice(ith_min:ith_max,iph_min:iph_max,-1:0,-1:0,-1:0))
      slice=0.

    endsubroutine alloc_rslice_scal
!***********************************************************************
    subroutine alloc_rslice_scal_2D(slice)

      real, dimension(:,:), allocatable :: slice

      if (allocated(slice)) deallocate(slice)
      allocate(slice(ith_min:ith_max,iph_min:iph_max))
      slice=0.

    endsubroutine alloc_rslice_scal_2D
!***********************************************************************
    subroutine alloc_rslice_vec(slice)

      real, dimension(:,:,:,:,:,:), allocatable :: slice

      if (allocated(slice)) deallocate(slice)
      allocate(slice(ith_min:ith_max,iph_min:iph_max,-1:0,-1:0,-1:0,3))
      slice=0.

    endsubroutine alloc_rslice_vec
!***********************************************************************
    subroutine store_rslice_scal(pencil,rslice)

      real, dimension(nx) :: pencil
      real, dimension(ith_min:ith_max,iph_min:iph_max,-1:0,-1:0,-1:0) :: rslice

      integer :: ith, iph, mdif, ndif, lind

      do ith=ith_min,ith_max
        do iph=iph_min,iph_max
          lind=rslice_adjec_corn_inds(ith,iph,1)
          if (lind/=0) then                                     !  proc has point
            mdif=m-rslice_adjec_corn_inds(ith,iph,2)
            if (mdif==-1.or.mdif==0) then                       ! m appropriate
              ndif=n-rslice_adjec_corn_inds(ith,iph,3)
              if (ndif==-1.or.ndif==0) then
                lind=lind-nghost
                rslice(ith,iph,:,mdif,ndif)=pencil((/max(lind-1,1),lind/)) ! take relevant part of pencil
              endif
            endif
          endif
        enddo
      enddo

    endsubroutine store_rslice_scal
!***********************************************************************
    subroutine store_rslice_vec(pencil,rslice)

      real, dimension(nx,3) :: pencil
      real, dimension(ith_min:ith_max,iph_min:iph_max,-1:0,-1:0,-1:0,3) :: rslice

      integer :: ith, iph, mdif, ndif, lind

      do ith=ith_min,ith_max
        do iph=iph_min,iph_max
          lind=rslice_adjec_corn_inds(ith,iph,1)
          if (lind/=0) then                               !  proc has point
            mdif=m-rslice_adjec_corn_inds(ith,iph,2)
            if (mdif==-1.or.mdif==0) then                 ! m appropriate
              ndif=n-rslice_adjec_corn_inds(ith,iph,3)
              if (ndif==-1.or.ndif==0) then               ! n appropriate
                !if (iproc==3) print*, 'ith,iph,inds=', ith,iph,lind,mdif,ndif
                lind=lind-nghost
                rslice(ith,iph,:,mdif,ndif,:)=pencil((/max(lind-1,1),lind/),:) ! take relevant part of pencil
              endif
            endif
          endif
        enddo
      enddo
    
    endsubroutine store_rslice_vec
!***********************************************************************
    subroutine interp_rslice(slice,f,ind,rslice)

      use Messages, only: fatal_error

      real, dimension(ith_min:ith_max,iph_min:iph_max) :: slice
      real, dimension(mx,my,mz,*), optional :: f
      integer, optional :: ind
      real, dimension(ith_min:ith_max,iph_min:iph_max,-1:0,-1:0,-1:0), optional :: rslice

      integer :: ith,iph,lind,mind,nind

      do ith=ith_min,ith_max
        do iph=iph_min,iph_max
          lind=rslice_adjec_corn_inds(ith,iph,1)
          if (lind/=0) then      !  proc has point
            if (present(f)) then
              mind=rslice_adjec_corn_inds(ith,iph,2)
              nind=rslice_adjec_corn_inds(ith,iph,3)
              slice(ith,iph)=sum(f(lind-1:lind,mind-1:mind,nind-1:nind,ind)*rslice_interp_weights(ith,iph,:,:,:))
            elseif (present(rslice)) then 
              slice(ith,iph)=sum(rslice(ith,iph,:,:,:)*rslice_interp_weights(ith,iph,:,:,:))
            else
              call fatal_error('interp_rslice', 'either f and ind or rslice must be provided')
            endif
          endif
        enddo
      enddo
    
    endsubroutine interp_rslice
!***********************************************************************
    subroutine write_rslice_position(unit)

      integer :: unit
      integer :: iph_min_

      if (iph_min==max_int) then
        iph_min_=0
      else
        iph_min_=iph_min
      endif

      write (unit, '(l5,4(1x,i5)," R")') lwrite_slice_r, ith_min,ith_max, iph_min_,iph_max

    endsubroutine write_rslice_position
!***********************************************************************
endmodule Slices_methods
