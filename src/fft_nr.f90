!
! fft_nr.f
!
! Numerical Recipes (complex-to-complex) fast Fourier transform routine
!
      SUBROUTINE four1(data,nn,isign)

      INTEGER isign,nn
      REAL data(2*nn)

      ! Replaces data(1:2*nn) by its discrete Fourier transform, if
      ! isign is input as 1; or replaces data(1:2*nn) by nn times its
      ! inverse discrete Fourier transform, if isign is input as -1.
      ! data is a complex array of length nn or, equivalently, a real
      ! array of length 2*nn. nn MUST be an integer power of 2 (this is
      ! not checked for!).

      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp ! Double precision for
                                                 ! the trigonometric
                                                 ! recurrences
      n=2*nn
      j=1
      do i=1,n,2                ! Bit-reversal section of the routine
         if (j.gt.i)then
            tempr=data(j)       ! Exchange the two complex numbers
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
         endif
         m=n/2
         do while ((m.ge.2).and.(j.gt.m))
            j=j-m
            m=m/2
         enddo
         j=j+m
      enddo
      mmax=2                    ! Here begins the Danielson-Lanczos
                                ! section of the routine
      do while (n.gt.mmax)      ! Outer loop executed log2 nn times
         istep=2*mmax
         theta=6.28318530717959d0/(isign*mmax) ! Initialize for the
                                               ! trigonometric recurrence
         wpr=-2.d0*sin(0.5d0*theta)**2
         wpi=sin(theta)
         wr=1.d0
         wi=0.d0
         do m=1,mmax,2          ! Here are the two nested inner loops

            do i=m,n,istep
               j=i+mmax         ! This is the Danielson-Lanczos formula:
               tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
               tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
               data(j)=data(i)-tempr
               data(j+1)=data(i+1)-tempi
               data(i)=data(i)+tempr
               data(i+1)=data(i+1)+tempi
            enddo
            wtemp=wr            ! Trigonometric recurrence
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
         enddo
         mmax=istep
      enddo
      return

      END

! End of file
