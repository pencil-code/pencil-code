      subroutine fft(a,b,ntot,n,nspan,isn)
!
!  dummy routine for fft.f
!
        use Cdata, only: lfft
!
      real :: a,b
      integer :: ntot,n,nspan,isn
      logical, save :: first=.true.
!
      if(first) then
        lfft=.false.
        print*
        print*,'fft is needed (eg for potential field bc)'
        print*,'you MUST use FFT=fft in Makefile.local'
        print*,'We proceed anyway...'
        if (a==1.) print*,b,ntot,n,nspan,isn
        print*
      endif
      first=.false.
      end
