!  $Id: nofftpack.f,v 1.6 2003-10-01 13:11:41 theine Exp $
!
!  Dummy routine, to avoid never seeing the compiler warnings from fftpack.
!
!***********************************************************************
      SUBROUTINE CFFTF (N,C,WSAVE)
      DIMENSION C(*),WSAVE(*)
      save iprint
!
!  Dummy routine: put c=0 to prevent further damage
!
      if(iprint.eq.0) print*,'Use FFTPACK=fftpack; ',n,c(1),wsave(1)
      do i=1,n
        c(i)=0.
      enddo
      iprint=1
      if (.false.) print*,wsave(1) ! (keep compiler quiet)
      END
!***********************************************************************
      SUBROUTINE CFFTI (N,WSAVE)
      DIMENSION WSAVE(*)
      save iprint
!
!  Dummy routine: put c=0 to prevent further damage
!
      do i=1,4*n+15
        wsave(i)=0.
      enddo
      if(iprint.eq.0) print*,'Use FFTPACK=fftpack; ',n
      iprint=1
      END
!***********************************************************************
      SUBROUTINE CFFTB (N,C,WSAVE)
      DIMENSION C(*), WSAVE(*)
      save iprint
!
!  Dummy routine: put c=0 to prevent further damage
!
      do i=1,4*n+15
        wsave(i)=0.
      enddo
      if(iprint.eq.0) print*,'Use FFTPACK=fftpack; ',n,c(1),wsave(1)
      iprint=1
      END
!***********************************************************************
      SUBROUTINE COSQB (N,X,WSAVE)
      DIMENSION X(*),WSAVE(*)
      if (.false.) print*,n,x(1),wsave(1) ! (keep compiler quiet)
      END
!***********************************************************************
      SUBROUTINE COSQF (N,X,WSAVE)
      DIMENSION X(*),WSAVE(*)
      if (.false.) print*,n,x(1),wsave(1) ! (keep compiler quiet)
      END
!***********************************************************************
      SUBROUTINE COSQI (N,WSAVE)
      DIMENSION WSAVE(*)
      if (.false.) print*,n,wsave(1) ! (keep compiler quiet)
      END
!***********************************************************************
