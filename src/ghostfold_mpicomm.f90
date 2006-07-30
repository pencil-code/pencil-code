! $Id: ghostfold_mpicomm.f90,v 1.4 2006-07-30 22:36:21 mee Exp $
!
!  This module performs some special mpifunctions that 
!  also require the Fourier routines. 
!
module GhostFold

  use Cparam
  use Messages

  implicit none

  private

  public :: fold_df, fold_f

  contains
!***********************************************************************
    subroutine fold_df(df,ivar1,ivar2)
!
!  Fold first ghost zone of df into main part of df.
!
!  20-may-2006/anders: coded
!
      use Cdata
      use Mpicomm
      use Fourier
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: ivar1, ivar2
!      
      real, dimension (nx+2,ny+2,1,ivar2-ivar1+1) :: df_tmp_xy
      real, dimension (nx+2,1,nz,ivar2-ivar1+1) :: df_tmp_xz
      real, dimension (ny,nz) :: df_tmp_yz
      integer :: nvar_fold, iproc_rcv, ivar
!
      nvar_fold=ivar2-ivar1+1
!  The first ghost zone in the z-direction is folded (the corners will find
!  their proper positions after folding in x and y as well).
      if (nzgrid/=1) then
        if (nprocz==1) then
          df(l1-1:l2+1,m1-1:m2+1,n1,ivar1:ivar2)= &
              df(l1-1:l2+1,m1-1:m2+1,n1,ivar1:ivar2)+ &
              df(l1-1:l2+1,m1-1:m2+1,n2+1,ivar1:ivar2)
          df(l1-1:l2+1,m1-1:m2+1,n2,ivar1:ivar2)= &
              df(l1-1:l2+1,m1-1:m2+1,n2,ivar1:ivar2)+ &
              df(l1-1:l2+1,m1-1:m2+1,n1-1,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xy, &
                  (/nx+2,ny+2,1,nvar_fold/), zlneigh, 10000)
            elseif (iproc_rcv==zuneigh) then
              call mpisend_real(df(l1-1:l2+1,m1-1:m2+1,n2+1:n2+1,ivar1:ivar2), &
                  (/nx+2,ny+2,1,nvar_fold/), zuneigh, 10000)
            endif
            if (iproc==iproc_rcv) df(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2)= &
                df(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2) + df_tmp_xy
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xy, &
                  (/nx+2,ny+2,1,nvar_fold/), zuneigh, 10001)
            elseif (iproc_rcv==zlneigh) then
              call mpisend_real(df(l1-1:l2+1,m1-1:m2+1,n1-1:n1-1,ivar1:ivar2), &
                  (/nx+2,ny+2,1,nvar_fold/), zlneigh, 10001)
            endif
            if (iproc==iproc_rcv) df(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2)= &
                df(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2) + df_tmp_xy
          enddo
        endif
        df(l1-1:l2+1,m1-1:m2+1,n1-1,ivar1:ivar2)=0.0
        df(l1-1:l2+1,m1-1:m2+1,n2+1,ivar1:ivar2)=0.0
      endif
! Then y.
      if (nygrid/=1) then
        if (nprocy==1) then
          df(l1-1:l2+1,m1,n1:n2,ivar1:ivar2)= &
              df(l1-1:l2+1,m1,n1:n2,ivar1:ivar2) + &
              df(l1-1:l2+1,m2+1,n1:n2,ivar1:ivar2)
          df(l1-1:l2+1,m2,n1:n2,ivar1:ivar2)= &
              df(l1-1:l2+1,m2,n1:n2,ivar1:ivar2) + &
              df(l1-1:l2+1,m1-1,n1:n2,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xz, &
                  (/nx+2,1,nz,nvar_fold/), ylneigh, 10002)
            elseif (iproc_rcv==yuneigh) then
              call mpisend_real(df(l1-1:l2+1,m2+1:m2+1,n1:n2,ivar1:ivar2), &
                  (/nx+2,1,nz,nvar_fold/), yuneigh, 10002)
            endif
            if (iproc==iproc_rcv) df(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2)= &
                df(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2) + df_tmp_xz
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xz, &
                  (/nx+2,1,nz,nvar_fold/), yuneigh, 10003)
            elseif (iproc_rcv==ylneigh) then
              call mpisend_real(df(l1-1:l2+1,m1-1:m1-1,n1:n2,ivar1:ivar2), &
                  (/nx+2,1,nz,nvar_fold/), ylneigh, 10003)
            endif
            if (iproc==iproc_rcv) df(l1-1:l2+1,m2:m2,n1:n2,ivar1:ivar2)= &
                df(l1-1:l2+1,m2:m2,n1:n2,ivar1:ivar2) + df_tmp_xz
          enddo
!
        endif
!
!  With shearing boundary conditions we must take care that the information is
!  shifted properly before the final fold.
!
        if (nxgrid>1 .and. lshear) then
          do ivar=ivar1,ivar2
            df_tmp_yz=df(l1-1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz(df_tmp_yz,-deltay)
            df(l1-1,m1:m2,n1:n2,ivar)=df_tmp_yz
            df_tmp_yz=df(l2+1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz(df_tmp_yz,+deltay)
            df(l2+1,m1:m2,n1:n2,ivar)=df_tmp_yz
          enddo
        endif
!
        df(l1-1:l2+1,m1-1,n1:n2,ivar1:ivar2)=0.0
        df(l1-1:l2+1,m2+1,n1:n2,ivar1:ivar2)=0.0
      endif
! Then x (always at the same processor).
      if (nxgrid/=1) then
        df(l1,m1:m2,n1:n2,ivar1:ivar2)=df(l1,m1:m2,n1:n2,ivar1:ivar2) + &
            df(l2+1,m1:m2,n1:n,ivar1:ivar2)
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
!  20-may-2006/anders: coded
!
      use Cdata
      use Mpicomm
      use Fourier
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: ivar1, ivar2
!      
      real, dimension (nx+2,ny+2,1,ivar2-ivar1+1) :: f_tmp_xy
      real, dimension (nx+2,1,nz,ivar2-ivar1+1) :: f_tmp_xz
      real, dimension (ny,nz) :: f_tmp_yz
      integer :: nvar_fold, iproc_rcv, ivar
!
      nvar_fold=ivar2-ivar1+1
!  The first ghost zone in the z-direction is folded (the corners will find
!  their proper positions after folding in x and y as well).
      if (nzgrid/=1) then
        if (nprocz==1) then
          f(l1-1:l2+1,m1-1:m2+1,n1,ivar1:ivar2)= &
              f(l1-1:l2+1,m1-1:m2+1,n1,ivar1:ivar2)+ &
              f(l1-1:l2+1,m1-1:m2+1,n2+1,ivar1:ivar2)
          f(l1-1:l2+1,m1-1:m2+1,n2,ivar1:ivar2)= &
              f(l1-1:l2+1,m1-1:m2+1,n2,ivar1:ivar2)+ &
              f(l1-1:l2+1,m1-1:m2+1,n1-1,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xy, &
                  (/nx+2,ny+2,1,nvar_fold/), zlneigh, 10000)
            elseif (iproc_rcv==zuneigh) then
              call mpisend_real(f(l1-1:l2+1,m1-1:m2+1,n2+1:n2+1,ivar1:ivar2), &
                  (/nx+2,ny+2,1,nvar_fold/), zuneigh, 10000)
            endif
            if (iproc==iproc_rcv) f(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2)= &
                f(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2) + f_tmp_xy
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xy, &
                  (/nx+2,ny+2,1,nvar_fold/), zuneigh, 10001)
            elseif (iproc_rcv==zlneigh) then
              call mpisend_real(f(l1-1:l2+1,m1-1:m2+1,n1-1:n1-1,ivar1:ivar2), &
                  (/nx+2,ny+2,1,nvar_fold/), zlneigh, 10001)
            endif
            if (iproc==iproc_rcv) f(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2)= &
                f(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2) + f_tmp_xy
          enddo
        endif
        f(l1-1:l2+1,m1-1:m2+1,n1-1,ivar1:ivar2)=0.0
        f(l1-1:l2+1,m1-1:m2+1,n2+1,ivar1:ivar2)=0.0
      endif
! Then y.
      if (nygrid/=1) then
        if (nprocy==1) then
          f(l1-1:l2+1,m1,n1:n2,ivar1:ivar2)= &
              f(l1-1:l2+1,m1,n1:n2,ivar1:ivar2) + &
              f(l1-1:l2+1,m2+1,n1:n2,ivar1:ivar2)
          f(l1-1:l2+1,m2,n1:n2,ivar1:ivar2)= &
              f(l1-1:l2+1,m2,n1:n2,ivar1:ivar2) + &
              f(l1-1:l2+1,m1-1,n1:n2,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xz, &
                  (/nx+2,1,nz,nvar_fold/), ylneigh, 10002)
            elseif (iproc_rcv==yuneigh) then
              call mpisend_real(f(l1-1:l2+1,m2+1:m2+1,n1:n2,ivar1:ivar2), &
                  (/nx+2,1,nz,nvar_fold/), yuneigh, 10002)
            endif
            if (iproc==iproc_rcv) f(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2)= &
                f(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2) + f_tmp_xz
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xz, &
                  (/nx+2,1,nz,nvar_fold/), yuneigh, 10003)
            elseif (iproc_rcv==ylneigh) then
              call mpisend_real(f(l1-1:l2+1,m1-1:m1-1,n1:n2,ivar1:ivar2), &
                  (/nx+2,1,nz,nvar_fold/), ylneigh, 10003)
            endif
            if (iproc==iproc_rcv) f(l1-1:l2+1,m2:m2,n1:n2,ivar1:ivar2)= &
                f(l1-1:l2+1,m2:m2,n1:n2,ivar1:ivar2) + f_tmp_xz
          enddo
!
        endif
!
!  With shearing boundary conditions we must take care that the information is
!  shifted properly before the final fold.
!
        if (nxgrid>1 .and. lshear) then
          do ivar=ivar1,ivar2
            f_tmp_yz=f(l1-1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz(f_tmp_yz,-deltay)
            f(l1-1,m1:m2,n1:n2,ivar)=f_tmp_yz
            f_tmp_yz=f(l2+1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz(f_tmp_yz,+deltay)
            f(l2+1,m1:m2,n1:n2,ivar)=f_tmp_yz
          enddo
        endif
!
        f(l1-1:l2+1,m1-1,n1:n2,ivar1:ivar2)=0.0
        f(l1-1:l2+1,m2+1,n1:n2,ivar1:ivar2)=0.0
      endif
! Then x (always at the same processor).
      if (nxgrid/=1) then
        f(l1,m1:m2,n1:n2,ivar1:ivar2)=f(l1,m1:m2,n1:n2,ivar1:ivar2) + &
            f(l2+1,m1:m2,n1:n,ivar1:ivar2)
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
      use Mpicomm
      use Fourier
!
      real, dimension (ny,nz) :: a_re
      real :: shift_y
!
      complex, dimension(nygrid) :: a_cmplx, cmplx_shift
      real, dimension(nygrid,max(nz/nprocy,1)) :: a_re_new, a_im_new
      integer :: n, nz_new, ipy_from, ipy_to, iproc_from, iproc_to
      real, dimension(4*nygrid+15) :: wsavey
!
!  Fourier transform of the subdivided y-interval is done by collecting
!  pencils of length nygrid at each processor. Consider the processor space
!  for a given ipz for nygrid=32, nz=8:
!
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    |               |               |               |               | 8
!    |               |               |               |               | 7
!    |               |               |               |               | 6
!    |               |               |               |               | 5
!  z |     ipy=0     |     ipy=1     |     ipy=2     |     ipy=3     | 4
!    |               |               |               |               | 3
!    |               |               |               |               | 2
!    |               |               |               |               | 1
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                                    y
!
!  This is the resulting processor division that can be used for the Fourier
!  transform:
!
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    |                             ipy=3                             | 8
!    |_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _| 7
!    |                             ipy=2                             | 6
!    |_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _| 5
!  z |                             ipy=1                             | 4
!    |_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _| 3
!    |                             ipy=0                             | 2
!    |                                                               | 1
!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                                    y
!
!  The height at each processor is nz/nprocy.
!
      nz_new=max(nz/nprocy,1)
      if (nprocy/=1) then
        if (nzgrid==1) then
!
!  Degenerate z-direction. Let root processor do the shift of the single nygrid
!  pencil.
!
          do iproc_from=1,ncpus-1
            if (lroot) then
              call mpirecv_real( &
                  a_re_new(iproc_from*ny+1:(iproc_from+1)*ny,1), &
                  ny, iproc_from, 666)
            else
              iproc_to=0
              if (iproc==iproc_from) &
                  call mpisend_real(a_re(:,1), ny, iproc_to, 666)
            endif
          enddo
          if (lroot) a_re_new(1:ny,1)=a_re(:,1)
        else
!
!  Present z-direction. Here nz must be a whole multiple of nprocy (e.g. nz=8,
!  nprocy=4). This constraint could be removed with a bit of work if it should
!  become necessary.
!
          if (modulo(nz,nprocy)/=0) then
            if (lroot) print*, 'fourier_shift_yz: nz must be a whole '// &
                'multiple of nprocy!'
            call fatal_error('fourier_shift_yz','')
          endif
!
          do ipy_from=0,nprocy-1
            iproc_from=ipz*nprocy+ipy_from
            if (ipy/=ipy_from) call mpirecv_real( &
                a_re_new(ipy_from*ny+1:(ipy_from+1)*ny,:), &
                (/ny,nz_new/), iproc_from, 666)
            if (ipy==ipy_from) then
              a_re_new(ipy*ny+1:(ipy+1)*ny,:) = &
                  a_re(:,ipy*nz_new+1:(ipy+1)*nz_new)
              do ipy_to=0,nprocy-1
                iproc_to=ipz*nprocy+ipy_to
                if (ipy/=iproc_to) call mpisend_real( &
                    a_re(:,ipy_to*nz_new+1:(ipy_to+1)*nz_new), &
                    (/ny,nz_new/), iproc_to, 666)
              enddo
            endif
          enddo
        endif
      else
!
!  Only parallelization along z (or not at all).
!
        a_re_new(1:ny,1:nz_new)=a_re(1:ny,1:nz_new)
      endif
!
      a_im_new=0.0
      cmplx_shift=exp(cmplx(0.0,-ky_fft*shift_y))
!
!  Transform to Fourier space.
!
      do n=1,nz_new
        call fourier_transform_other(a_re_new(:,n),a_im_new(:,n),+1)
        a_cmplx=cmplx(a_re_new(:,n),a_im_new(:,n))
        a_cmplx=a_cmplx*cmplx_shift
        a_re_new(:,n)=real(a_cmplx)
        a_im_new(:,n)=aimag(a_cmplx)
      enddo
!
!  Back to real space.
!
      do n=1,nz
        call fourier_transform_other(a_re_new(:,n),a_im_new(:,n),-1)
      enddo
!
!  Reinstate original division of yz-plane.
!
      if (nprocy/=1) then
        if (nzgrid==1) then
!  No z-direction.
          if (.not. lroot) then
            iproc_from=0
            call mpirecv_real(a_re(:,1), ny, iproc_from, 666)
          else
            do iproc_to=1,ncpus-1
              call mpisend_real( &
                  a_re_new(iproc_to*ny+1:(iproc_to+1)*ny,1), &
                  ny, iproc_to, 666)
            enddo
          endif
          if (lroot) a_re(:,1)=a_re_new(1:ny,1)
        else
!  Present z-direction.
          do ipy_from=0,nprocy-1
            iproc_from=ipz*nprocy+ipy_from
            if (ipy/=ipy_from) call mpirecv_real( &
                a_re(:,ipy_from*nz_new+1:(ipy_from+1)*nz_new), &
                (/ny,nz_new/), iproc_from, 666)
            if (ipy==ipy_from) then
              a_re(:,ipy*nz_new+1:(ipy+1)*nz_new)= &
                  a_re_new(ipy*ny+1:(ipy+1)*ny,:)
              do ipy_to=0,nprocy-1
                iproc_to=ipz*nprocy+ipy_to
                if (ipy/=iproc_to) call mpisend_real( &
                    a_re_new(ipy_to*ny+1:(ipy_to+1)*ny,:), &
                    (/ny,nz_new/), iproc_to, 666)
              enddo
            endif
          enddo
        endif
      else
!  Only parallelization along z (or not at all).
        a_re(1:ny,1:nz_new)=a_re_new(1:ny,1:nz_new)
      endif
!
    endsubroutine fourier_shift_yz
!***********************************************************************
endmodule GhostFold
