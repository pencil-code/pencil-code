!
! Type definition for scattered array with time stamp.
!
  type scattered_array
! 
! Dimensions
!
    integer :: dim1,dim2,dim3,dim4
!
! Data
!
    real, dimension(:)      , allocatable :: time
    real, dimension(:,:,:)  , allocatable :: data3
    real, dimension(:,:,:,:), allocatable :: data4
!
! Link to next list element
!
    type (scattered_array), pointer :: next => null()
!
  endtype scattered_array
!
  public :: scattered_array, init_scattered_array, store_scattered_array, get_scattered_array

  interface store_scattered_array
    module procedure store_scattered_array3
    module procedure store_scattered_array4
  endinterface

  interface get_scattered_array
    module procedure get_scattered_array3
    module procedure get_scattered_array4
  endinterface
!
  contains
!***********************************************************************
    subroutine free_scattered_array(arr)
!
!  Comment me.
!
      type (scattered_array), pointer :: arr
      type (scattered_array), pointer :: next
!
      do while (associated(arr))
        next=>arr%next
        if (allocated(arr%time)) deallocate(arr%time)
        if (allocated(arr%data3)) deallocate(arr%data3)
        if (allocated(arr%data4)) deallocate(arr%data4)
        deallocate(arr)
        arr=>next
      enddo
!
      nullify(arr)
!
    endsubroutine free_scattered_array
!***********************************************************************
   subroutine init_scattered_array(arr,dim1,dim2,dim3,dim4,lreset)

     type (scattered_array), pointer :: arr
     logical, optional :: lreset
     integer :: dim1,dim2
     integer, optional :: dim3,dim4
!
     if (loptest(lreset)) then
       call free_scattered_array(arr)
     else
       NULLIFY(arr)
     endif
!
     allocate(arr)

     arr%dim1=dim1
     arr%dim2=dim2
     arr%dim3=ioptest(dim3,0)
     arr%dim4=ioptest(dim4,0)
     arr%time=0.

   endsubroutine init_scattered_array
!***********************************************************************
   subroutine store_scattered_array3(layer,src,dest)

     integer :: layer
     real, dimension(:,:) :: src
     type(scattered_array), pointer :: dest

     integer, parameter :: dim3=50
     integer :: niter, i
     type(scattered_array), pointer :: ptr

     niter=ceiling(real(layer)/dim3)
     ptr=>dest

     do i=1,niter
       if (.not.associated(ptr)) call init_scattered_array(ptr,dest%dim1,dest%dim2)
       if (.not.allocated(ptr%data3)) allocate(ptr%data3(ptr%dim1,ptr%dim2,ptr%dim3))
       if (i==niter) exit
       ptr=>ptr%next
     enddo

     ptr%data3(:,:,layer-(niter-1)*dim3) = src

   endsubroutine store_scattered_array3
!***********************************************************************
   subroutine store_scattered_array4(ind3,layer,src,dest,time)

     use Cdata, only: iproc

     integer :: ind3,layer
     real, dimension(:,:) :: src
     type(scattered_array), pointer :: dest
     real :: time

     integer :: niter, i, ind4
     type(scattered_array), pointer :: ptr, newptr

     niter=ceiling(real(layer)/dest%dim4)
     newptr=>dest; nullify(ptr)

     do i=1,niter
!if (ind3==1.and.iproc==120) print*, 'layer,iter,niter,assoc=', layer,i,niter,associated(newptr)
       if (.not.associated(newptr)) then
         call init_scattered_array(newptr,dest%dim1,dest%dim2,dest%dim3,dest%dim4)
         if (associated(ptr)) ptr%next=>newptr
       endif
       ptr=>newptr
       if (.not.allocated(ptr%data4)) allocate(ptr%data4(ptr%dim1,ptr%dim2,ptr%dim3,ptr%dim4),ptr%time(ptr%dim4))
       if (i==niter) exit
       newptr=>ptr%next
     enddo

     ind4=layer-(niter-1)*ptr%dim4
     ptr%time(ind4)=time
     ptr%data4(:,:,ind3,ind4) = src
!if (ind3==1.and.iproc==120) print*, 'layer,niter,maxval,time,assoc=', layer,niter,maxval(abs(ptr%data4(:,:,ind3,layer-(niter-1)*dim4))),ptr%time,associated(ptr%next)

   endsubroutine store_scattered_array4
!***********************************************************************
   subroutine get_scattered_array3(layer,src,dest)

     integer :: layer
     real, dimension(:,:) :: dest
     type(scattered_array), pointer :: src

     integer, parameter :: dim3=50
     integer :: niter, i
     type(scattered_array), pointer :: ptr

     niter=ceiling(real(layer)/dim3)
     ptr=>src
     do i=1,niter
       if (.not.associated(ptr)) stop
       if (i==niter) exit
       ptr=>ptr%next
     enddo

     dest=ptr%data3(:,:,layer-(niter-1)*dim3)

   endsubroutine get_scattered_array3
!***********************************************************************
   subroutine get_scattered_array4(ivar,layer,src,dest,timediff,ahead)

     use Cdata, only:ldiagnos, iproc
     integer :: ivar,layer
     real, dimension(:,:) :: dest
     real, dimension(:,:), optional :: ahead
     type(scattered_array), pointer :: src
     real :: timediff

     integer :: niter,i,ind4
     type(scattered_array), pointer :: ptr

     niter=ceiling(real(layer)/src%dim4)
     ptr=>src
     do i=1,niter
       if (.not.associated(ptr)) stop
       if (i==niter) exit
       ptr=>ptr%next
     enddo

     ind4=layer-(niter-1)*ptr%dim4
     dest=ptr%data4(:,:,ivar,ind4)
!if (ldiagnos.and.iproc==120) print*, 'in timediff:', ptr%time(ind4), maxval(ptr%data4(:,:,ivar,ind4))

     if (ind4<ptr%dim4) then
       timediff=ptr%time(ind4+1)-ptr%time(ind4)
       if (present(ahead)) ahead=ptr%data4(:,:,ivar,ind4+1)
     elseif (associated(ptr%next)) then
       timediff=ptr%next%time(1)-ptr%time(ind4)
       if (present(ahead)) ahead=ptr%next%data4(:,:,ivar,1)
     else
       timediff=0.
     endif

   endsubroutine get_scattered_array4
