! $Id: random_nr.f90,v 1.3 2002-09-04 03:17:18 brandenb Exp $ 

module Random_nr

  use Cparam

  implicit none

  contains

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
    function ran(idum)
!
!  From `Numerical Recipes in F90'. Not sure we are allowed to distribute
!  this
!
! 28-aug-02/wolf: Adapted from Numerical Recipes
!
      implicit none
      integer, parameter :: k4b=selected_int_kind(9)
      integer(k4b), intent(inout) :: idum
      real :: ran

      ! "Minimal" random number generator of Park and Miller combined
      ! with a Marsaglia shift sequence. Returns a uniform random deviate
      ! between 0.0 and 1.0 (exclusive of the endpoint values). This fully
      ! portable, scalar generator has the "traditional" (not Fortran 90)
      ! calling sequence with a random deviate as the returned function
      ! value: call with idum a negative integer to initialize;
      ! thereafter, do not alter idum except to reinitialize. The period
      ! of this generator is about 3.1×10^18.

      integer(k4b), parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
      real, save :: am
      integer(k4b), save :: ix=-1,iy=-1,k
      if (idum <= 0 .or. iy < 0) then   ! Initialize.
        am=nearest(1.0,-1.0)/im
        iy=ior(ieor(888889999,abs(idum)),1)
        ix=ieor(777755555,abs(idum))
        idum=abs(idum)+1   ! Set idum positive.
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
    subroutine min_std(matrix)
!
!  Produces a matrix filled with random numbers calculated
!  with the 'Minimal Standard' random number generator
!      
      real, dimension(mx,my,mz) :: matrix
      integer :: i,j,k
      integer, save :: idum=-10 ! _must_ be < 0 for ran()
!
      do i=1,mx
         do j=1,my
            do k=1,mz
               matrix(i,j,k)=ran(idum)
            enddo
         enddo
      enddo
!
    end subroutine min_std
!***********************************************************************

end module Random_nr
