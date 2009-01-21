! $Id$

module General

  implicit none

  contains

!***********************************************************************
    subroutine tridag(a,b,c,r,u,err)
!
!  solves tridiagonal system
!
!  01-apr-03/tobi: from numerical recipes
!
      implicit none
      real, dimension(:), intent(in) :: a,b,c,r
      real, dimension(:), intent(out) :: u
      real, dimension(size(b)) :: gam
      logical, intent(out), optional :: err
      integer :: n,j
      real :: bet

      if (present(err)) err=.false.
      n=size(b)
      bet=b(1)
      if (bet.eq.0.) then
         print*,'tridag: Error at code stage 1'
         if (present(err)) err=.true.
         return
      endif

      u(1)=r(1)/bet
      do j=2,n
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j-1)*gam(j)
         if (bet.eq.0.) then
            print*,'tridag: Error at code stage 2'
            if (present(err)) err=.true.
            return
         endif
         u(j)=(r(j)-a(j-1)*u(j-1))/bet
      end do
      do j=n-1,1,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
      end do
    end subroutine tridag
!***********************************************************************
    subroutine tridag_double(a,b,c,r,u,err)
!
!  solves tridiagonal system
!
!  01-apr-03/tobi: from numerical recipes
!  11-apr-03/axel: double precision version
!
      implicit none
      double precision, dimension(:), intent(in) :: a,b,c,r
      double precision, dimension(:), intent(out) :: u
      double precision, dimension(size(b)) :: gam
      logical, intent(out), optional :: err
      integer :: n,j
      double precision :: bet

      if (present(err)) err=.false.
      n=size(b)
      bet=b(1)
      if (bet.eq.0.) then
         print*,'tridag_double: Error at code stage 1'
         if (present(err)) err=.true.
         return
      endif

      u(1)=r(1)/bet
      do j=2,n
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j-1)*gam(j)
         if (bet.eq.0.) then
            print*,'tridag_double: Error at code stage 2'
            if (present(err)) err=.true.
            return
         endif
         u(j)=(r(j)-a(j-1)*u(j-1))/bet
      end do
      do j=n-1,1,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
      end do
    end subroutine tridag_double
!***********************************************************************

end module General
