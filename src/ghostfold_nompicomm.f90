! $Id: ghostfold_nompicomm.f90,v 1.2 2006-07-30 15:05:45 ajohan Exp $
!
!  This module performs some special mpifunctions that 
!  also require the Fourier routines. 
!
module GhostFold

  use Cparam

  implicit none

  private

  public :: fold_df, fold_f

  contains
!***********************************************************************
    subroutine fold_df(df,ivar1,ivar2)
!
!  Fold first ghost zone of df into main part of df.
!
!  15-may-2006/anders: coded
!
      use Cdata
      use Fourier
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: ivar1,ivar2
!
      real, dimension (ny,nz) :: df_tmp_yz
      integer :: ivar
!
!  Fold z-direction first (including first ghost zone in x and y).
!
      if (nzgrid/=1) then
        df(l1-1:l2+1,m1-1:m2+1,n1,ivar1:ivar2)= &
            df(l1-1:l2+1,m1-1:m2+1,n1,ivar1:ivar2) + &
            df(l1-1:l2+1,m1-1:m2+1,n2+1,ivar1:ivar2)
        df(l1-1:l2+1,m1-1:m2+1,n2,ivar1:ivar2)= &
            df(l1-1:l2+1,m1-1:m2+1,n2,ivar1:ivar2) + &
            df(l1-1:l2+1,m1-1:m2+1,n1-1,ivar1:ivar2)
        df(l1-1:l2+1,m1-1:m2+1,n1-1,ivar1:ivar2)=0.0
        df(l1-1:l2+1,m1-1:m2+1,n2+1,ivar1:ivar2)=0.0
      endif
!
!  Then fold y-direction (including first ghost zone in x).
!
      if (nygrid/=1) then
        df(l1-1:l2+1,m1,n1:n2,ivar1:ivar2)= &
            df(l1-1:l2+1,m1,n1:n2,ivar1:ivar2) + &
            df(l1-1:l2+1,m2+1,n1:n2,ivar1:ivar2)
        df(l1-1:l2+1,m2,n1:n2,ivar1:ivar2)= &
            df(l1-1:l2+1,m2,n1:n2,ivar1:ivar2) + &
            df(l1-1:l2+1,m1-1,n1:n2,ivar1:ivar2)
!
!  With shearing boundary conditions we must take care that the information is
!  shifted properly before the final fold.
!
        if (lshear) then
          do ivar=ivar1,ivar2
            df_tmp_yz=df(l1-1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz(df_tmp_yz,-deltay)
            df(l1-1,m1:m2,n1:n2,ivar)=df_tmp_yz
            df_tmp_yz=df(l2+1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz(df_tmp_yz,+deltay)
            df(l2+1,m1:m2,n1:n2,ivar)=df_tmp_yz
          enddo
        endif
        df(l1-1:l2+1,m1-1,n1:n2,ivar1:ivar2)=0.0
        df(l1-1:l2+1,m2+1,n1:n2,ivar1:ivar2)=0.0
      endif
!
!  Finally fold the x-direction.
!
      if (nxgrid/=1) then
        df(l1,m1:m2,n1:n2,ivar1:ivar2)=df(l1,m1:m2,n1:n2,ivar1:ivar2) + &
            df(l2+1,m1:m2,n1:n2,ivar1:ivar2) 
        df(l2,m1:m2,n1:n2,ivar1:ivar2)=df(l2,m1:m2,n1:n2,ivar1:ivar2) + &
            df(l1-1,m1:m2,n1:n2,ivar1:ivar2) 
        df(l1-1,m1:m2,n1:n2,ivar1:ivar2)=0.0
        df(l2+1,m1:m2,n1:n2,ivar1:ivar2)=0.0
      endif
!
    endsubroutine fold_df
!***********************************************************************
    subroutine fold_f(f,ivar1,ivar2)
!
!  Fold first ghost zone of f into main part of f.
!
!  14-jun-2006/anders: adapted
!
      use Cdata
      use Fourier
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: ivar1, ivar2
!
      real, dimension (ny,nz) :: f_tmp_yz
      integer :: ivar
!
!  Fold z-direction first (including first ghost zone in x and y).
!
      if (nzgrid/=1) then
        f(l1-1:l2+1,m1-1:m2+1,n1,ivar1:ivar2)= &
            f(l1-1:l2+1,m1-1:m2+1,n1,  ivar1:ivar2) + &
            f(l1-1:l2+1,m1-1:m2+1,n2+1,ivar1:ivar2)
        f(l1-1:l2+1,m1-1:m2+1,n2,ivar1:ivar2)= &
            f(l1-1:l2+1,m1-1:m2+1,n2,  ivar1:ivar2) + &
            f(l1-1:l2+1,m1-1:m2+1,n1-1,ivar1:ivar2)
        f(l1-1:l2+1,m1-1:m2+1,n1-1,ivar1:ivar2)=0.0
        f(l1-1:l2+1,m1-1:m2+1,n2+1,ivar1:ivar2)=0.0
      endif
!
!  Then fold y-direction (including first ghost zone in x).
!
      if (nygrid/=1) then
        f(l1-1:l2+1,m1,n1:n2,ivar1:ivar2)=f(l1-1:l2+1,m1,n1:n2,ivar1:ivar2) + &
            f(l1-1:l2+1,m2+1,n1:n2,ivar1:ivar2)
        f(l1-1:l2+1,m2,n1:n2,ivar1:ivar2)=f(l1-1:l2+1,m2,n1:n2,ivar1:ivar2) + &
            f(l1-1:l2+1,m1-1,n1:n2,ivar1:ivar2)
!
!  With shearing boundary conditions we must take care that the information is
!  shifted properly before the final fold.
!
        if (lshear) then
          do ivar=ivar1,ivar2
            f_tmp_yz=f(l1-1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz(f_tmp_yz,-deltay)
            f(l1-1,m1:m2,n1:n2,ivar)=f_tmp_yz
            f_tmp_yz=f(l2+1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz(f_tmp_yz,+deltay)
            f(l2+1,m1:m2,n1:n2,ivar)=f_tmp_yz
          enddo
        endif
        f(l1-1:l2+1,m1-1,n1:n2,ivar1:ivar2)=0.0
        f(l1-1:l2+1,m2+1,n1:n2,ivar1:ivar2)=0.0
      endif
!
!  Finally fold the x-direction.
!
      if (nxgrid/=1) then
        f(l1,m1:m2,n1:n2,ivar1:ivar2)=f(l1,m1:m2,n1:n2,ivar1:ivar2) + &
            f(l2+1,m1:m2,n1:n2,ivar1:ivar2) 
        f(l2,m1:m2,n1:n2,ivar1:ivar2)=f(l2,m1:m2,n1:n2,ivar1:ivar2) + &
            f(l1-1,m1:m2,n1:n2,ivar1:ivar2) 
        f(l1-1,m1:m2,n1:n2,ivar1:ivar2)=0.0
        f(l2+1,m1:m2,n1:n2,ivar1:ivar2)=0.0
      endif
!
    endsubroutine fold_f
!***********************************************************************
    subroutine fourier_shift_yz(a_re,shift_y)
!
!  Performs a periodic shift in the y-direction of an entire y-z plane by
!  the amount shift_y. The shift is done in Fourier space for maximum
!  interpolation accuracy.
!
!  19-jul-06/anders: coded
!
      use Cdata, only: ky_fft
      use Fourier
!
      real, dimension (ny,nz) :: a_re
      real :: shift_y
!
      real, dimension(ny,nz) :: a_im
      complex, dimension(ny) :: a_cmplx
      complex, dimension(ny) :: cmplx_shift
      integer :: n
!
      a_im=0.0
      cmplx_shift=exp(cmplx(0.0,-ky_fft*shift_y))
!
!  Transform to Fourier space.
!
      do n=1,nz
        call fourier_transform_other(a_re(:,n),a_im(:,n),+1)
        a_cmplx=cmplx(a_re(:,n),a_im(:,n))
        a_cmplx=a_cmplx*cmplx_shift
        a_re(:,n)=real(a_cmplx)
        a_im(:,n)=aimag(a_cmplx)
      enddo
!
!  Back to real space.
!
      do n=1,nz
        call fourier_transform_other(a_re(:,n),a_im(:,n),-1)
      enddo
!
    endsubroutine fourier_shift_yz
!***********************************************************************
endmodule GhostFold
