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
