! Id: nompicomm.f90,v 1.35 2002/08/16 21:23:48 brandenb Exp $

!!!!!!!!!!!!!!!!!!!!!!!
!!!  nompicomm.f90  !!!
!!!!!!!!!!!!!!!!!!!!!!!

!!!  Module with dummy MPI stuff.
!!!  This allows code to be run on single cpu machine

module Mpicomm

  use Cparam
  use Cdata, only: iproc,ipx,ipy,ipz,root,lroot

  implicit none

  interface mpibcast_real       ! Overload the `mpibcast_real' function
    module procedure mpibcast_real_scl
    module procedure mpibcast_real_arr
  endinterface

  integer, dimension (ny*nz) :: mm,nn
  integer :: ierr,imn
  logical, dimension (ny*nz) :: necessary=.false.

  contains

!***********************************************************************
    subroutine mpicomm_init()
!
!  Before the communication has been completed, the nhost=3 layers next
!  to the processor boundary (m1, m2, n1, or n2) cannot be used yet.
!  In the mean time we can calculate the interior points sufficiently far
!  away from the boundary points. Here we calculate the order in which
!  m and n are executed. At one point, necessary(imn)=.true., which is
!  the moment when all communication must be completed.
!
!   6-jun-02/axel: generalized to allow for ny=1
!
      use General
      use Cdata, only: lmpicomm,directory
!
!  sets iproc in order that we write in the correct directory
!
      character (len=4) :: chproc
      integer :: m,n
!
!  for single cpu machine, set processor to zero
!
      lmpicomm = .false.
      iproc = 0
      lroot = .true.
      ipx = 0
      ipy = 0
      ipz = 0
!
!  produce index-array for the sequence of points to be worked through:
!  first inner box, then boundary zones.
!  Could be somehow simplified here (no communication), but we need to
!  update the ghost zones before using them, so it is best to stick to
!  the same scheme as in Mpicomm.
!
      imn=1
      do n=n1i+1,n2i-1
        do m=m1i+1,m2i-1
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
      if (ip < 10) print*,'NOMPICOMM: setting necessary(',imn,') = .true.' 
      necessary(imn)=.true.
!
!  do the lower stripe in the n-direction
!
      do n=max(n2i,n1+1),n2
        do m=m1i+1,m2i-1
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
!
!  upper stripe in the n-direction
!
      do n=n1,min(n1i,n2)
        do m=m1i+1,m2i-1
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
!
!  left and right hand boxes
!  NOTE: need to have min(m1i,m2) instead of just m1i, and max(m2i,m1+1)
!  instead of just m2i, to make sure the case ny=1 works ok, and
!  also that the same m is not set in both loops.
!
      do n=n1,n2
        do m=m1,min(m1i,m2)
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
        do m=max(m2i,m1+1),m2
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
      use Cparam
!
!  apply boundary conditions
!
      real, dimension (mx,my,mz,mvar) :: f
!
      if (ip==0) print*,'FINALIZE_ISENDRCV_BDRY: f=',f
    endsubroutine finalise_isendrcv_bdry
!***********************************************************************
    subroutine initiate_shearing(f)
!
      real, dimension (mx,my,mz,mvar) :: f
      real :: dummy
!    
      dummy=f(1,1,1,1)
    endsubroutine initiate_shearing
!***********************************************************************
    subroutine finalise_shearing(f)
!
  use Cdata
!
!  for one processor, use periodic boundary conditions
!  but in this dummy routine this is done in finalise_isendrcv_bdry
!
      real, dimension (mx,my,mz,mvar) :: f
      double precision    :: deltay_dy, frak, c1, c2, c3, c4, c5, c6
      integer :: displs
!
!  Periodic boundary conditions in x, with shearing sheat
!
      if (nygrid==1) then !If 2D
         f( 1:l1-1,:,:,:) = f(l2i:l2,:,:,:)
         f(l2+1:mx,:,:,:) = f(l1:l1i,:,:,:)
      else
         deltay_dy=deltay/dy
         displs=int(deltay_dy)
         frak=deltay_dy-displs
         c1 = -          (frak+1.)*frak*(frak-1.)*(frak-2.)*(frak-3.)/120.
         c2 = +(frak+2.)          *frak*(frak-1.)*(frak-2.)*(frak-3.)/24.
         c3 = -(frak+2.)*(frak+1.)     *(frak-1.)*(frak-2.)*(frak-3.)/12.
         c4 = +(frak+2.)*(frak+1.)*frak          *(frak-2.)*(frak-3.)/12.
         c5 = -(frak+2.)*(frak+1.)*frak*(frak-1.)          *(frak-3.)/24.
         c6 = +(frak+2.)*(frak+1.)*frak*(frak-1.)*(frak-2.)          /120.
         f( 1:l1-1,m1:m2,:,:)=c1*cshift(f(l2i:l2,m1:m2,:,:),-displs+2,2) &
                             +c2*cshift(f(l2i:l2,m1:m2,:,:),-displs+1,2) &
                             +c3*cshift(f(l2i:l2,m1:m2,:,:),-displs,2) &
                             +c4*cshift(f(l2i:l2,m1:m2,:,:),-displs-1,2) &
                             +c5*cshift(f(l2i:l2,m1:m2,:,:),-displs-2,2) &
                             +c6*cshift(f(l2i:l2,m1:m2,:,:),-displs-3,2)  
         f(l2+1:mx,m1:m2,:,:)=c1*cshift(f(l1:l1i,m1:m2,:,:),displs-2,2) &
                             +c2*cshift(f(l1:l1i,m1:m2,:,:),displs-1,2) &
                             +c3*cshift(f(l1:l1i,m1:m2,:,:),displs,2) &
                             +c4*cshift(f(l1:l1i,m1:m2,:,:),displs+1,2) &
                             +c5*cshift(f(l1:l1i,m1:m2,:,:),displs+2,2) &
                             +c6*cshift(f(l1:l1i,m1:m2,:,:),displs+3,2) 
      end if
    end subroutine finalise_shearing
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
    subroutine mpibarrier()
    endsubroutine mpibarrier
!***********************************************************************
    subroutine mpifinalize()
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
