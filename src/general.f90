! $Id: general.f90,v 1.19 2003-04-03 14:56:55 theine Exp $

module General

!  Module with general utility subroutines
!  (Used for example in Sub and Mpicomm)

  use Cparam
!  use Cdata, only: nseed,seed

  implicit none

  interface random_number_wrapper   ! Overload this function
    module procedure random_number_wrapper_1
    module procedure random_number_wrapper_3
  endinterface
!
!  state and default generator of random numbers
!
  integer, save, dimension(mseed) :: rstate=0
  character (len=labellen) :: random_gen='system'

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
      select case(random_gen)
!
      case('system')
        call random_number(a)
      case('min_std')
        do i=1,size(a,1)
          a(i)=ran0(rstate(1))
        enddo
      case default ! 'nr_f90'
        do i=1,size(a,1)
          a(i)=mars_ran()
        enddo
      endselect
!
    endsubroutine random_number_wrapper_1
!***********************************************************************
    subroutine random_number_wrapper_3(a)
!
!  Produces a matrix filled with random numbers calculated
!  with the 'Minimal Standard' random number generator
!
      real, dimension(:,:,:) :: a
      integer :: i,j,k
!
      select case(random_gen)
!
      case('system')
        call random_number(a)
      case('min_std')
      do i=1,size(a,1)
        do j=1,size(a,2)
          do k=1,size(a,3)
            a(i,j,k)=ran0(rstate(1))
          enddo
        enddo
      enddo
      case default ! 'nr_f90'
      do i=1,size(a,1)
        do j=1,size(a,2)
          do k=1,size(a,3)
            a(i,j,k)=mars_ran()
          enddo
        enddo
      enddo
      endselect
!
    endsubroutine random_number_wrapper_3
!***********************************************************************
    subroutine random_seed_wrapper(size,put,get)
!
!  mimics the f90 random_seed routine
!
      real :: dummy
      integer, optional, dimension(:) :: put,get
      integer, optional :: size
      integer :: nseed
!
      select case(random_gen)
!
      case('system')
        call random_seed(SIZE=nseed)
        if(present(size)) size=nseed
        if(present(get))  call random_seed(GET=get(1:nseed))
        if(present(put))  call random_seed(PUT=put(1:nseed))
      case('min_std')
        nseed=1
        if(present(size)) size=nseed
        if(present(get))  get=rstate(1)
        if(present(put))  rstate(1)=min(put(1),-1)  ! ensure rstate(1) < 0
      case default ! 'nr_f90'
        nseed=2
        if(present(size)) size=nseed
        if(present(get))  get=rstate(1:nseed)
        if(present(put)) then
          if (put(2)==0) then   ! state cannot be result from previous
                                ! call, so initialize
            dummy = mars_ran(put(1))
          else
            rstate(1:nseed)=put
          endif
        endif
      endselect
!
    endsubroutine random_seed_wrapper
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
    endfunction ran0

!***********************************************************************
    function mars_ran(init)
!
! 26-sep-02/wolf: Adapted from `Numerical Recipes for F90' ran() routine
!
! "Minimal" random number generator of Park and Miller combined
! with a Marsaglia shift sequence. Returns a uniform random deviate
! between 0.0 and 1.0 (exclusive of the endpoint values).
! Call with (INIT=ival) to initialize.
! The period of this generator is supposed to be about 3.1× 10^18.
!
      implicit none
!
      real :: mars_ran
      real, save :: am             ! will be constant on a given platform
      integer, optional, intent(in) :: init
      integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
      integer :: k,init1=1812   ! default value
!
      if (present(init) .or. rstate(1)==0 .or. rstate(2)<=0) then
        !
        ! initialize
        !
        if (present(init)) init1 = init
        am=nearest(1.0,-1.0)/im
        rstate(1)=ieor(777755555,abs(init1))
        rstate(2)=ior(ieor(888889999,abs(init1)),1)
      endif
      !
      ! Marsaglia shift sequence with period 2^32-1
      !
      rstate(1)=ieor(rstate(1),ishft(rstate(1),13))
      rstate(1)=ieor(rstate(1),ishft(rstate(1),-17))
      rstate(1)=ieor(rstate(1),ishft(rstate(1),5)) 
      !
      ! Park-Miller sequence by Schrage's method, period 2^31-2
      !
      k=rstate(2)/iq
      rstate(2)=ia*(rstate(2)-k*iq)-ir*k
      if (rstate(2) < 0) rstate(2)=rstate(2)+im
      !
      ! combine the two generators with masking to ensure nonzero value
      !
      mars_ran=am*ior(iand(im,ieor(rstate(1),rstate(2))),1)
!
    endfunction mars_ran
!***********************************************************************
    function ran(iseed1)
!
!  (More or less) original routine from `Numerical Recipes in F90'. Not
!  sure we are allowed to distribute this
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
    endfunction ran
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
    subroutine safe_character_assign(dest,src)
!
!  Do character string assignement with check against overflow
!
!  08-oct-02/tony: coded
!  25-oct-02/axel: added tag in output to give name of routine
!
      character (len=*), intent(in):: src
      character (len=*), intent(inout):: dest
      integer :: destLen, srcLen

      destLen = len(dest)
      srcLen = len(src)

      if (destLen<srcLen) then 
         print *, "safe_character_assign: ", &
              "RUNTIME ERROR: FORCED STRING TRUNCATION WHEN ASSIGNING '" & 
               //src//"' to '"//dest//"'"
         dest=src(1:destLen)
      else
         dest=src
      end if

    endsubroutine safe_character_assign
!***********************************************************************
    subroutine safe_character_append(str1,str2)
!
!  08-oct-02/wolf: coded
!
      character (len=*), intent(inout):: str1
      character (len=*), intent(in):: str2
!
      call safe_character_assign(str1, trim(str1) // str2)
!
    endsubroutine safe_character_append
!***********************************************************************
    function spline_derivative(z,f)
!
!  computes derivative of a given function using spline interpolation
!
!  01-apr-03/tobi: originally coded by Aake Nordlund
!
      implicit none
      real, dimension(:) :: z
      real, dimension(:), intent(in) :: f
      real, dimension(size(z)) :: w1,w2,w3
      real, dimension(size(z)) :: d,spline_derivative
      real :: c
      integer :: mz,k

      mz=size(z)

      w1(1)=1./(z(2)-z(1))**2
      w3(1)=-1./(z(3)-z(2))**2
      w2(1)=w1(1)+w3(1)
      d(1)=2.*((f(2)-f(1))/(z(2)-z(1))**3 &
                  -(f(3)-f(2))/(z(3)-z(2))**3)
!
! interior points
!
      w1(2:mz-1)=1./(z(2:mz-1)-z(1:mz-2))
      w3(2:mz-1)=1./(z(3:mz)-z(2:mz-1))
      w2(2:mz-1)=2.*(w1(2:mz-1)+w3(2:mz-1))

      d(2:mz-1)=3.*(w3(2:mz-1)**2*(f(3:mz)-f(2:mz-1)) &
           +w1(2:mz-1)**2*(f(2:mz-1)-f(1:mz-2)))

!
! last point
!
      w1(mz)=1./(z(mz-1)-z(mz-2))**2
      w3(mz)=-1./(z(mz)-z(mz-1))**2
      w2(mz)=w1(mz)+w3(mz)
      d(mz)=2.*((f(mz-1)-f(mz-2))/(z(mz-1)-z(mz-2))**3 &
           -(f(mz)-f(mz-1))/(z(mz)-z(mz-1))**3)
!
! eliminate at first point
!
      c=-w3(1)/w3(2)
      w1(1)=w1(1)+c*w1(2)
      w2(1)=w2(1)+c*w2(2)
      d(1)=d(1)+c*d(2)
      w3(1)=w2(1)
      w2(1)=w1(1)
!
! eliminate at last point
!
      c=-w1(mz)/w1(mz-1)
      w2(mz)=w2(mz)+c*w2(mz-1)
      w3(mz)=w3(mz)+c*w3(mz-1)
      d(mz)=d(mz)+c*d(mz-1)
      w1(mz)=w2(mz)
      w2(mz)=w3(mz)
!
! eliminate subdiagonal
!
      do k=2,mz
         c=-w1(k)/w2(k-1)
         w2(k)=w2(k)+c*w3(k-1)
         d(k)=d(k)+c*d(k-1)
      end do
!
! backsubstitute
!
      d(mz)=d(mz)/w2(mz)
      do k=mz-1,1,-1
         d(k)=(d(k)-w3(k)*d(k+1))/w2(k)
      end do

      spline_derivative=d
    end function spline_derivative
!***********************************************************************
    function spline_integral(z,f,q0)
!
!  computes integral of a given function using spline interpolation
!
!  01-apr-03/tobi: originally coded by Aake Nordlund
!
      implicit none
      real, dimension(:) :: z
      real, dimension(:) :: f
      real, dimension(size(z)) :: df,dz
      real, dimension(size(z)) :: q,spline_integral
      real, optional :: q0
      integer :: mz,k

      mz=size(z)

      q(1)=0.
      if (present(q0)) q(1)=q0
      df=spline_derivative(z,f)
      dz(2:mz)=z(2:mz)-z(1:mz-1)

      q(2:mz)=.5*dz(2:mz)*(f(1:mz-1)+f(2:mz)) &
              +(1./12.)*dz(2:mz)**2*(df(1:mz-1)-df(2:mz))

      do k=2,mz
         q(k)=q(k)+q(k-1)
      end do

      spline_integral=q
    end function spline_integral
!***********************************************************************
    subroutine tridag(a,b,c,r,u,err)
!
!  solves tridiagonal system
!
!  01-apr/tobi: from numerical recipes
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
         print*,'tridag_ser: Error at code stage 1'
         if (present(err)) err=.true.
         return
      endif

      u(1)=r(1)/bet
      do j=2,n
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j-1)*gam(j)
         if (bet.eq.0.) then
            print*,'tridag_ser: Error at code stage 2'
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
end module General
