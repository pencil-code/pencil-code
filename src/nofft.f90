      subroutine fft(a,b,ntot,n,nspan,isn)
!
!  dummy routine for fft.f
!
      real :: a,b
      integer :: ntot,n,nspan,isn
      logical, save :: first=.true.
!
      if(first) then
        print*
        print*,'fft is needed (eg for potential field bc)'
        print*,'you MUST use FFT=fft in Makefile.local'
        print*,'We proceed anyway...'
        if (a==1.) print*,b,ntot,n,nspan,isn
        print*
      endif
      first=.false.
      end
