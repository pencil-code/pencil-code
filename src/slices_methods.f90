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

!  interface alloc_slice_buffers
!    module procedure alloc_slice_buffers_scal
!    module procedure alloc_slice_buffers_vec
!  endinterface

  contains
!***********************************************************************
!  include 'slice_methods.f90.removed'
!***********************************************************************
    subroutine assign_slices_sep_scal(slices,xy,xz,yz,xy2,xy3,xy4,xz2)
!
!  Assignment of slice pointers according to slice selection switches.
!
!  12-apr-17/MR: coded 
!   8-oct-18/MR: made some pars optional
!
      type (slice_data),                      intent(OUT):: slices
      real, dimension(:,:), target,           intent(IN) :: xy
      real, dimension(:,:), target, optional, intent(IN) :: xz,yz,xy2,xy3,xy4,xz2

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

      if (lwrite_slice_yz ) slices%yz => yz(:,:,slices%index)
      if (lwrite_slice_xz ) slices%xz => xz(:,:,slices%index)
      if (lwrite_slice_xy ) slices%xy => xy(:,:,slices%index)
      if (lwrite_slice_xy2) slices%xy2=>xy2(:,:,slices%index)
      if (lwrite_slice_xy3) slices%xy3=>xy3(:,:,slices%index)
      if (lwrite_slice_xy4) slices%xy4=>xy4(:,:,slices%index)
      if (lwrite_slice_xz2) slices%xz2=>xz2(:,:,slices%index)

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

      if (slices%index<=nc) slices%ready=.true.
    endif

    endsubroutine assign_slices_f_vec
!***********************************************************************
    subroutine store_slices_scal(pencil,xy,xz,yz,xy2,xy3,xy4,xz2)
!
!  Stores scalar data from pencil in auxiliary arrays when needed.
!
!  12-apr-17/MR: coded 
!   8-oct-18/MR: made some pars optional
!
      real, dimension(nx) , intent(IN)           :: pencil
      real, dimension(:,:), intent(OUT)          :: xy
      real, dimension(:,:), intent(OUT),optional :: xz,yz,xy2,xy3,xy4,xz2

      if (n==iz_loc) xy(:,m-m1+1)=pencil
      if (present(xz) .and.m==iy_loc)  xz(:,n-n1+1)=pencil
      if (present(yz) .and.lwrite_slice_yz) yz(m-m1+1,n-n1+1)=pencil(ix_loc-l1+1)
      if (present(xy2).and.n==iz2_loc) xy2(:,m-m1+1)=pencil
      if (present(xy3).and.n==iz3_loc) xy3(:,m-m1+1)=pencil
      if (present(xy4).and.n==iz4_loc) xy4(:,m-m1+1)=pencil
      if (present(xz2).and.m==iy2_loc) xz2(:,n-n1+1)=pencil

    endsubroutine store_slices_scal
!***********************************************************************
    subroutine store_slices_vec(pencil,xy,xz,yz,xy2,xy3,xy4,xz2,ncomp)
!
!  Stores multi-component data from pencil in auxiliary arrays when needed.
!
!  12-apr-17/MR: coded 
!
      use General, only: transform_thph_yy, ioptest

      real, dimension(:,:)  , intent(IN) :: pencil
      real, dimension(:,:,:), intent(OUT):: xy,xz,yz,xy2,xy3,xy4,xz2
      integer, optional     , intent(IN) :: ncomp

      real, dimension(1,3) :: tmp
      
      !nc=ioptest(ncomp,3)
!
!  Note: i[xyz][234]*_loc is the index with respect to array with ghost zones.
!
      if (lwrite_slice_yz) then
        if (lyang.and..not.present(ncomp).and.size(pencil,2)==3) then
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

    endsubroutine process_slices_fac
!***********************************************************************
    subroutine addto_slices(slices,pencil)
!
!  Adds pencil data to scalar slices according to slice selection switches.
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

    endsubroutine nullify_slice_pointers
!***********************************************************************
endmodule Slices_methods
