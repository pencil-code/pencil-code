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

  interface alloc_slice_buffers
    module procedure alloc_slice_buffers_scal
    module procedure alloc_slice_buffers_vec
  endinterface

  interface store_slices
    module procedure store_slices_scal
    module procedure store_slices_vec
  endinterface
!
  contains
!***********************************************************************
    subroutine assign_slices_sep_scal(slices,xy,xz,yz,xy2,xy3,xy4,xz2)
!
!  Assignment of slice pointers according to slice selection switches.
!
!  12-apr-17/MR: coded 
!
      type (slice_data),            intent(OUT):: slices
      real, dimension(:,:), target, intent(IN) :: xy,xz,yz,xy2,xy3,xy4,xz2

      if (lwrite_slice_yz ) slices%yz =>yz
      if (lwrite_slice_xz ) slices%xz =>xz
      if (lwrite_slice_xy ) slices%xy =>xy
      if (lwrite_slice_xy2) slices%xy2=>xy2
      if (lwrite_slice_xy3) slices%xy3=>xy3
      if (lwrite_slice_xy4) slices%xy4=>xy4
      if (lwrite_slice_xz2) slices%xz2=>xz2

      slices%ready=.true.

    endsubroutine assign_slices_sep_scal
!***********************************************************************
    subroutine assign_slices_sep_vec(slices,xy,xz,yz,xy2,xy3,xy4,xz2,ncomp)
!
!  Assignment of slice pointers according to slice selection switches.
!
!  12-apr-17/MR: coded 
!
    use General, only: ioptest
!
    type (slice_data)             , intent(INOUT):: slices
    real, dimension(:,:,:), target, intent(OUT)  :: xy,xz,yz,xy2,xy3,xy4,xz2
    integer, optional             , intent(IN)   :: ncomp

    integer :: nc

    nc=ioptest(ncomp,3)

    if (slices%index>=nc) then
      slices%ready=.false.
    else
      slices%index=slices%index+1

      if (lwrite_slice_yz ) slices%yz =>yz(:,:,slices%index)
      if (lwrite_slice_xz ) slices%xz =>xz(:,:,slices%index)
      if (lwrite_slice_xy ) slices%xy =>xy(:,:,slices%index)
      if (lwrite_slice_xy2) slices%xy2=>xy2(:,:,slices%index)
      if (lwrite_slice_xy3) slices%xy3=>xy3(:,:,slices%index)
      if (lwrite_slice_xy4) slices%xy4=>xy4(:,:,slices%index)
      if (lwrite_slice_xz2) slices%xz2=>xz2(:,:,slices%index)

      if (slices%index<=nc) slices%ready=.true.
    endif

    endsubroutine assign_slices_sep_vec
!***********************************************************************
    subroutine assign_slices_f_scal(slices,f,ind)
!
!  Assignment of slice pointers into f-array according to slice selection switches.
!
!  12-apr-17/MR: coded 
!
      type (slice_data)          , intent(OUT):: slices
      real, dimension(mx,my,mz,*), intent(IN) :: f
      integer                    , intent(IN) :: ind

      if (lwrite_slice_yz) &
        slices%yz=f(ix_loc,m1:m2  ,n1:n2  ,ind)
      if (lwrite_slice_xz) &
        slices%xz =f(l1:l2,iy_loc ,n1:n2  ,ind)
      if (lwrite_slice_xy) &
        slices%xy =f(l1:l2,m1:m2  ,iz_loc ,ind)
      if (lwrite_slice_xy2) &
        slices%xy2=f(l1:l2,m1:m2  ,iz2_loc,ind)
      if (lwrite_slice_xy3) &
        slices%xy3=f(l1:l2,m1:m2  ,iz3_loc,ind)
      if (lwrite_slice_xy4) &
        slices%xy4=f(l1:l2,m1:m2  ,iz4_loc,ind)
      if (lwrite_slice_xz2) &
        slices%xz2=f(l1:l2,iy2_loc,n1:n2  ,ind)

      slices%ready=.true.

    endsubroutine assign_slices_f_scal
!***********************************************************************
    subroutine assign_slices_f_vec(slices,f,ind,ncomp)
!
!  Assignment of slice pointers into f-array according to slice selection
!  switches.
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
        if (lyang.and.slices%index>=2) then
!
!  On Yang grid: transform theta and phi components of vector to Yin-grid basis.
!  (phi component is saved in transformed for use in next call.
!
          if (slices%index==2) then
            if (.not.allocated(transformed)) allocate(transformed(1,ny,nz,2))
            call transform_thph_yy_other(f(ix_loc:ix_loc,m1:m2,n1:n2,ind+1:ind+2),transformed)
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

      if (slices%index<=nc) slices%ready=.true.
    endif

    endsubroutine assign_slices_f_vec
!***********************************************************************
    subroutine store_slices_scal(pencil,xy,xz,yz,xy2,xy3,xy4,xz2)
!
!  Stores data from pencil in auxiliary arrays when needed.
!
!  12-apr-17/MR: coded 
!
      real, dimension(nx) , intent(IN) :: pencil
      real, dimension(:,:), intent(OUT):: xy,xz,yz,xy2,xy3,xy4,xz2

      if (n==iz_loc)  xy(:,m-m1+1)=pencil
      if (m==iy_loc)  xz(:,n-n1+1)=pencil
      if (lwrite_slice_yz) yz(m-m1+1,n-n1+1)=pencil(ix_loc-l1+1)
      if (n==iz2_loc) xy2(:,m-m1+1)=pencil
      if (n==iz3_loc) xy3(:,m-m1+1)=pencil
      if (n==iz4_loc) xy4(:,m-m1+1)=pencil
      if (m==iy2_loc) xz2(:,n-n1+1)=pencil

    endsubroutine store_slices_scal
!***********************************************************************
    subroutine store_slices_vec(pencil,xy,xz,yz,xy2,xy3,xy4,xz2,ncomp)
!
!  Stores data from pencil in auxiliary arrays when needed.
!
!  12-apr-17/MR: coded 
!
      use General, only: transform_thph_yy, ioptest

      real, dimension(:,:)  , intent(IN) :: pencil
      real, dimension(:,:,:), intent(OUT):: xy,xz,yz,xy2,xy3,xy4,xz2
      integer, optional     , intent(IN) :: ncomp

      integer :: nc
      real, dimension(1,3) :: tmp
      
      nc=ioptest(ncomp,3)
!
!  Note: i[xyz][234]*_loc is the index with respect to array with ghost zones.
!
      if (lwrite_slice_yz) then
        if (lyang.and.nc==3) then
!
!  On Yang grid: transform theta and phi components of oo to Yin-grid basis.
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

    endsubroutine store_slices_vec
!***********************************************************************
    subroutine alloc_slice_buffers_scal(xy,xz,yz,xy2,xy3,xy4,xz2)
!
!  Allocation of auxiliary arrays for slices according to slice selection
!  switches.
!
!  12-apr-17/MR: coded 
!
      real,dimension(:,:),allocatable,intent(INOUT):: xy,xy2,xy3,xy4,xz,xz2,yz

      if (lwrite_slice_yz .and..not.allocated(yz)) allocate(yz(ny,nz))
      if (lwrite_slice_xz .and..not.allocated(xz)) allocate(xz(nx,nz))
      if (lwrite_slice_xy .and..not.allocated(xy)) allocate(xy(nx,ny))
      if (lwrite_slice_xy2.and..not.allocated(xy2)) allocate(xy2(nx,ny))
      if (lwrite_slice_xy3.and..not.allocated(xy3)) allocate(xy3(nx,ny))
      if (lwrite_slice_xy4.and..not.allocated(xy4)) allocate(xy4(nx,ny))
      if (lwrite_slice_xz2.and..not.allocated(xz2)) allocate(xz2(nx,nz))

    endsubroutine alloc_slice_buffers_scal
!***********************************************************************
    subroutine alloc_slice_buffers_vec(xy,xz,yz,xy2,xy3,xy4,xz2,ncomp)
!
!  Allocation of auxiliary arrays for slices according to slice selection
!  switches.
!
!  12-apr-17/MR: coded 
!
      use General, only: ioptest
!
      real,dimension(:,:,:),allocatable,intent(INOUT):: xy,xy2,xy3,xy4,xz,xz2,yz
      integer, optional,                intent(IN)   :: ncomp

      integer :: nc

      nc=ioptest(ncomp,3)

      if (lwrite_slice_yz .and..not.allocated(yz)) allocate(yz(ny,nz,nc))
      if (lwrite_slice_xz .and..not.allocated(xz)) allocate(xz(nx,nz,nc))
      if (lwrite_slice_xy .and..not.allocated(xy)) allocate(xy(nx,ny,nc))
      if (lwrite_slice_xy2.and..not.allocated(xy2)) allocate(xy2(nx,ny,nc))            
      if (lwrite_slice_xy3.and..not.allocated(xy3)) allocate(xy3(nx,ny,nc))      
      if (lwrite_slice_xy4.and..not.allocated(xy4)) allocate(xy4(nx,ny,nc))
      if (lwrite_slice_xz2.and..not.allocated(xz2)) allocate(xz2(nx,nz,nc))

    endsubroutine alloc_slice_buffers_vec
!***********************************************************************
    subroutine process_slices(slices,func)
!
!  Assignment of slice pointers according to slice selection switches.
!
!  12-apr-17/MR: coded 
!
      type(slice_data), intent(INOUT):: slices
      character(LEN=40),intent(IN)   :: func

      if (trim(func)=='exp') then
        if (lwrite_slice_yz ) slices%yz =exp(slices%yz)
        if (lwrite_slice_xz ) slices%xz =exp(slices%xz)
        if (lwrite_slice_xy ) slices%xy =exp(slices%xy)
        if (lwrite_slice_xy2) slices%xy2=exp(slices%xy2)
        if (lwrite_slice_xy3) slices%xy3=exp(slices%xy3)
        if (lwrite_slice_xy4) slices%xy4=exp(slices%xy4)
        if (lwrite_slice_xz2) slices%xz2=exp(slices%xz2)
      elseif (trim(func)=='alog') then
        if (lwrite_slice_yz ) slices%yz =alog(slices%yz)
        if (lwrite_slice_xz ) slices%xz =alog(slices%xz)
        if (lwrite_slice_xy ) slices%xy =alog(slices%xy)
        if (lwrite_slice_xy2) slices%xy2=alog(slices%xy2)
        if (lwrite_slice_xy3) slices%xy3=alog(slices%xy3)
        if (lwrite_slice_xy4) slices%xy4=alog(slices%xy4)
        if (lwrite_slice_xz2) slices%xz2=alog(slices%xz2)
      elseif (trim(func)=='abs') then
        if (lwrite_slice_yz ) slices%yz =abs(slices%yz)
        if (lwrite_slice_xz ) slices%xz =abs(slices%xz)
        if (lwrite_slice_xy ) slices%xy =abs(slices%xy)
        if (lwrite_slice_xy2) slices%xy2=abs(slices%xy2)
        if (lwrite_slice_xy3) slices%xy3=abs(slices%xy3)
        if (lwrite_slice_xy4) slices%xy4=abs(slices%xy4)
        if (lwrite_slice_xz2) slices%xz2=abs(slices%xz2)
      endif

    endsubroutine process_slices
!***********************************************************************
    subroutine addto_slices(slices,pencil)
!
!  Adds pencil data to slices according to slice selection switches.
!
!  12-apr-17/MR: coded 
!
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

    endsubroutine addto_slices
!***********************************************************************
endmodule Slices_methods
