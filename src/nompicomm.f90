!!!!!!!!!!!!!!!!!!!!!!!
!!!  nompicomm.f90  !!!
!!!!!!!!!!!!!!!!!!!!!!!

!!!  Module with dummy MPI stuff.
!!!  This allows code to be run on single cpu machine

module Mpicomm

  use Cparam
  use Cdata, only: lroot

  implicit none

  interface mpibcast_real       ! Overload the `mpibcast_real' function
    module procedure mpibcast_real_scl
    module procedure mpibcast_real_arr
  endinterface

  integer, dimension (nx*ny) :: mm,nn
  integer :: ierr,imn
  integer :: iproc,root=0
  integer :: ipy,ipz
  logical, dimension (nx*ny) :: necessary=.false.
!  logical :: lroot=.true.       ! is this the root process?

  character (LEN=12) :: directory

  contains

!***********************************************************************
    subroutine mpicomm_init()
!
      use General
      use Cdata, only: lmpicomm
!
!  sets iproc in order that we write in the correct directory
!
      character (LEN=4) :: chproc
      integer :: m,n
!
!  for single cpu machine, set processor to root value
!
      lmpicomm = .false.
      iproc=root
!
!  produce index-array the the sequence of points to be worked through first
!
      imn=1
      do n=n1,n2
        do m=m1,m2
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
!
!  produce a directory name for the output
!
      call chn(iproc,chproc)
      directory='tmp/proc'//chproc
!
    endsubroutine mpicomm_init
!***********************************************************************
    subroutine initiate_isendrcv_bdry(f)
!
      use Cdata
!
!  for one processor, use periodic boundary conditions
!  but in this dummy routine this is done in finalise_isendrcv_bdry
!
      real, dimension (mx,my,mz,mvar) :: f
      integer :: i,j
!
!  Periodic boundary conditions in x
!
      f( 1:l1-1,:,:,:) = f(l2i:l2,:,:,:)
      f(l2+1:mx,:,:,:) = f(l1:l1i,:,:,:)
!
!  Periodic boundary conditions in y
!
      f(:, 1:m1-1,:,:) = f(:,m2i:m2,:,:)
      f(:,m2+1:my,:,:) = f(:,m1:m1i,:,:)
!
!  Boundary conditions in z
!
      do j=1,mvar
        select case(ibc(j))
        case (0)                ! periodic
          f(:,:, 1:n1-1,j) = f(:,:,n2i:n2,j)
          f(:,:,n2+1:mz,j) = f(:,:,n1:n1i,j)
        case (1)                ! symmetry
          do i=1,nghost
            f(:,:,n1-i,j) = f(:,:,n1+i,j)
            f(:,:,n2+i,j) = f(:,:,n2-i,j)
          enddo
        case (-1)               ! antisymmetry
          f(:,:,n1,j) = 0.      ! ensure boundary values are zero (indep.
                                ! of initial conditions)
          f(:,:,n2,j) = 0.
          do i=1,nghost
            f(:,:,n1-i,j) = -f(:,:,n1+i,j)
            f(:,:,n2+i,j) = -f(:,:,n2-i,j)
          enddo
        case (-2)               ! antisymmetry relative to boundary value
          do i=1,nghost
            f(:,:,n1-i,j) = 2*f(:,:,n1,j)-f(:,:,n1+i,j)
            f(:,:,n2+i,j) = 2*f(:,:,n2,j)-f(:,:,n2-i,j)
          enddo
        case default
          if (lroot) &
               print*,"No such boundary condition ibc = ", ibc(j), " for j=", j
        endselect
      enddo
!
    endsubroutine initiate_isendrcv_bdry
!***********************************************************************
    subroutine finalise_isendrcv_bdry(f)
!
      use General
!
!  do periodic boundary conditions
!
      real, dimension (mx,my,mz,mvar) :: f
      real :: dummy
!
!HPF$ ALIGN (:,:,:,*) WITH tmpl :: f
!HPF$ SHADOW(0,3,3,0) :: f
!
      dummy=f(1,1,1,1)  !(prevent compiler warning "unused variable ...")
!
    endsubroutine finalise_isendrcv_bdry
!***********************************************************************
    subroutine mpibcast_int(ibcast_array,nbcast_array)
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: ibcast_array,dummy
!    
      dummy=ibcast_array
    endsubroutine mpibcast_int
!***********************************************************************
    subroutine mpibcast_real_scl(bcast_array,nbcast_array)
!
      integer :: nbcast_array
      real :: bcast_array,dummy
!
      if (nbcast_array/=1) stop "problem in mpibcast_real_scl"
      dummy=bcast_array
    endsubroutine mpibcast_real_scl
!***********************************************************************
    subroutine mpibcast_real_arr(bcast_array,nbcast_array)
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array,dummy
!
      dummy=bcast_array
    endsubroutine mpibcast_real_arr
!***********************************************************************
    subroutine mpireduce_max(fmax_tmp,fmax,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fmax_tmp, fmax
!
      fmax=fmax_tmp
    endsubroutine mpireduce_max
!***********************************************************************
    subroutine mpireduce_sum(fsum_tmp,fsum,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fsum_tmp,fsum
!
      fsum=fsum_tmp
    endsubroutine mpireduce_sum
!***********************************************************************
    subroutine mpifinalize
    endsubroutine mpifinalize
!***********************************************************************
    subroutine stop_it(msg)
!
!  Print message and stop
!  6-nov-01/wolf: coded
!
      character (LEN=*) :: msg
!      
      if (lroot) write(0,'(A,A)') 'STOPPED: ', msg
      call mpifinalize
      STOP
    endsubroutine stop_it
!***********************************************************************

endmodule Mpicomm
