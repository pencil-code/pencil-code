! $Id: general.f90,v 1.7 2002-09-26 14:09:04 brandenb Exp $

module General

!  Module with general utility subroutines
!  (Used for example in Sub and Mpicomm)

  implicit none

  interface random_number_wrapper   ! Overload this function
    module procedure random_number_wrapper_1
    module procedure random_number_wrapper_3
  endinterface
!
!  default seed for random numbers
!
  integer(selected_int_kind(9)), save :: iseed=-10

  contains

!***********************************************************************
    subroutine random_number_wrapper_1(a)
!
!  Produces a matrix filled with random numbers calculated
!  with the 'Minimal Standard' random number generator
!
      real, dimension(:) :: a
      integer :: i
!
      do i=1,size(a,1)
        a(i)=ran(iseed)
      enddo
!
    end subroutine random_number_wrapper_1
!***********************************************************************
    subroutine random_number_wrapper_3(a)
!
!  Produces a matrix filled with random numbers calculated
!  with the 'Minimal Standard' random number generator
!
      real, dimension(:,:,:) :: a
      integer :: i,j,k
!
      do i=1,size(a,1)
        do j=1,size(a,2)
          do k=1,size(a,3)
            a(i,j,k)=ran(iseed)
          enddo
        enddo
      enddo
!
    end subroutine random_number_wrapper_3
!***********************************************************************
    function ran0(dummy)
!
!  The 'Minimal Standard' random number generator
!  by Lewis, Goodman and Miller.
!
!  28.08.02/nils: Adapted from Numerical Recipes 
!
      integer,parameter :: ia=16807,im=2147483647,iq=127773,ir=2836, &
           mask=123459876
      real, parameter :: am=1./im
      real :: ran0
      integer :: dummy,k
!
      dummy=ieor(dummy,mask)
      k=dummy/iq
      dummy=ia*(dummy-k*iq)-ir*k
      if (dummy.lt.0) dummy=dummy+im
      ran0=am*dummy
      dummy=ieor(dummy,mask)
      return
!
    end function ran0
!***********************************************************************
    function ran(iseed1)
!
!  From `Numerical Recipes in F90'. Not sure we are allowed to distribute
!  this
!
! 28-aug-02/wolf: Adapted from Numerical Recipes
!
      implicit none
!
      integer(selected_int_kind(9)), intent(inout) :: iseed1
      real :: ran
!
      ! "Minimal" random number generator of Park and Miller combined
      ! with a Marsaglia shift sequence. Returns a uniform random deviate
      ! between 0.0 and 1.0 (exclusive of the endpoint values). This fully
      ! portable, scalar generator has the "traditional" (not Fortran 90)
      ! calling sequence with a random deviate as the returned function
      ! value: call with iseed1 a negative integer to initialize;
      ! thereafter, do not alter iseed1 except to reinitialize. The period
      ! of this generator is about 3.1×10^18.

      real, save :: am
      integer(selected_int_kind(9)), parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
      integer(selected_int_kind(9)), save :: ix=-1,iy=-1,k

      if (iseed1 <= 0 .or. iy < 0) then   ! Initialize.
        am=nearest(1.0,-1.0)/im
        iy=ior(ieor(888889999,abs(iseed1)),1)
        ix=ieor(777755555,abs(iseed1))
        iseed1=abs(iseed1)+1   ! Set iseed1 positive.
      end if
      ix=ieor(ix,ishft(ix,13))   ! Marsaglia shift sequence with period 2^32-1.
      ix=ieor(ix,ishft(ix,-17))
      ix=ieor(ix,ishft(ix,5))
      k=iy/iq   ! Park-Miller sequence by Schrage's method,
      iy=ia*(iy-k*iq)-ir*k   ! period 231 - 2.
      if (iy < 0) iy=iy+im
      ran=am*ior(iand(im,ieor(ix,iy)),1)   ! Combine the two generators with
                                           ! masking to ensure nonzero value.
    end function ran
!***********************************************************************
    subroutine chn(n,ch)
!
      character (len=4) :: ch
      integer :: n
!
      intent(in) :: n
!
!  make a character out of a number
!  take care of numbers that have less than 4 digits
!  30-sep-97/axel: coded
!
      ch='    '
      if (n.lt.0) stop 'chn: lt1'
      if (n.lt.10) then
        write(ch(1:1),'(i1)') n
      elseif (n.lt.100) then
        write(ch(1:2),'(i2)') n
      elseif (n.lt.1000) then
        write(ch(1:3),'(i3)') n
      elseif (n.lt.10000) then
        write(ch(1:4),'(i4)') n
      else
        print*,'n=',n
        stop "chn: n too large"
      endif
!
    endsubroutine chn
!***********************************************************************
    subroutine chk_time(label,time1,time2)
!
      integer :: time1,time2,count_rate
      character (len=*) :: label
!
!  prints time in seconds
!
      call system_clock(count=time2,count_rate=count_rate)
      print*,label,(time2-time1)/real(count_rate)
      time1=time2
!
    endsubroutine chk_time
!***********************************************************************
!
end module General
