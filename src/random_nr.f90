! $Id: random_nr.f90,v 1.1 2002-08-26 10:33:49 nilshau Exp $ 

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
      integer dummy,ia,im,iq,ir,MASK
      real ran0,am
      parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836 &
           ,MASK=123459876)
      integer k
!
      dummy=ieor(dummy,MASK)
      k=dummy/iq
      dummy=ia*(dummy-k*iq)-ir*k
      if (dummy.lt.0) dummy=dummy+im
      ran0=am*dummy
      dummy=ieor(dummy,MASK)
      return
!
    end function ran0
!***********************************************************************
    subroutine min_std(matrix)
!
!  Produces a matrix filled with random numbers calculated
!  with the 'Minimal Standard' random number generator
!      
      real, dimension(mx,my,mz) :: matrix
      integer :: i,j,k
!
      do i=1,mx
         do j=1,my
            do k=1,mz
               matrix(i,j,k)=ran0(10)
            enddo
         enddo
      enddo
!
    end subroutine min_std
!***********************************************************************

end module Random_nr
