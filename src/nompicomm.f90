! $Id: nompicomm.f90,v 1.23 2002-05-26 16:42:58 brandenb Exp $

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

  integer, dimension (ny*nz) :: mm,nn
  integer :: ierr,imn
  integer :: iproc,root=0
  integer :: ipx=0,ipy=0,ipz=0
  logical, dimension (nx*ny) :: necessary=.false.

  contains

!***********************************************************************
    subroutine mpicomm_init()
!
      use General
      use Cdata, only: lmpicomm,directory
      use Boundcond
!
!  sets iproc in order that we write in the correct directory
!
      character (len=4) :: chproc
      integer :: m,n
!
!  for single cpu machine, set processor to root value
!
      lmpicomm = .false.
      iproc=root
      lroot = .true.
!
!  produce index-array for the sequence of points to be worked through.
!  Trivial here, since no communication.
!  Need to do the boundary conditions right in the beginning
!
      imn=1
      necessary(imn)=.true.
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
      real :: dummy
!
      dummy=f(1,1,1,1)  !(prevent compiler warning "unused variable ...")
!
    endsubroutine initiate_isendrcv_bdry
!***********************************************************************
    subroutine finalise_isendrcv_bdry(f)
!
      use General
      use Boundcond
!
!  apply boundary conditions
!
      real, dimension (mx,my,mz,mvar) :: f
      character (len=160) :: errmesg
!
!  Boundary conditions in x
!
      call boundconds(f,errmesg)
      if (errmesg /= "") call stop_it(trim(errmesg))
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
      character (len=*) :: msg
!      
      if (lroot) write(0,'(A,A)') 'STOPPED: ', msg
      call mpifinalize
      STOP
    endsubroutine stop_it
!***********************************************************************

endmodule Mpicomm
