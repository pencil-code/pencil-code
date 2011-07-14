!***********************************************************************
program test
!
  integer, allocatable, dimension (:,:,:,:) :: f
  integer, dimension(1,1,1,1) :: a
!
  allocate(f(1,1,1,2))
!
  f(1,1,1,1)=1
  f(1,1,1,2)=2
!
  print*,'Passing a slice of the array to a subroutine. Output should be 2.'
  call sub(f(:,:,:,2:2))
!
  print*,'Copy the revelant slice beforing passing. Output should be 2.'
  a=f(:,:,:,2:2)
  call sub(a)
!
endprogram
!********************************************************************
subroutine sub(a)
  integer, dimension (1,1,1,1) :: a
  print*, a(1,1,1,1)
endsubroutine sub
!********************************************************************
