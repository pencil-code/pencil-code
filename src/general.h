!
! Type definition for scattered array with time stamp.
!
  type scattered_array
! 
! Dimensions
!
    integer :: dim1,dim2,dim3,dim4
!
! Link to next list element
!
    type (scattered_array), pointer :: next => null()
!
  endtype scattered_array
!
  public :: scattered_array, init_scattered_array, store_scattered_array, get_scattered_array

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

   endsubroutine init_scattered_array
!***********************************************************************
   subroutine store_scattered_array(ind3,layer,src,dest,time)

     integer :: ind3,layer
     real, dimension(:,:) :: src
     type(scattered_array), pointer :: dest
     real :: time

     print*,'store_scattered_array: not available for f95'
     stop

   endsubroutine store_scattered_array
!***********************************************************************
   subroutine get_scattered_array(ivar,layer,src,dest,timediff,ahead)

     use Cdata, only:ldiagnos
     integer :: ivar,layer
     real, dimension(:,:) :: dest
     real, dimension(:,:), optional :: ahead
     type(scattered_array), pointer :: src
     real :: timediff

     print*,'get_scattered_array: not available for f95'
     stop

   endsubroutine get_scattered_array
!***********************************************************************
