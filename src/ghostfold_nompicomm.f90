! $Id$
!
!  This module folds ghost zones for a single processor run.
!
!  The module adds the value of f or df in the first ghost zone to the first
!  physical point at the opposite side of the grid. This is useful for example
!  for letting particles assign drag force to the ghost zones and then adding
!  the total assigned drag force back into the physical domain over the
!  periodic boundaries.
!
module GhostFold
!
  use Cdata
  use Cparam
  use Fourier
!
  implicit none
!
  private
!
  public :: fold_df, fold_f
!
  contains
!***********************************************************************
    subroutine fold_df(df,ivar1,ivar2)
!
!  Fold first ghost zone of df into main part of df.
!
!  15-may-2006/anders: coded
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
        if (lshear .and. nxgrid>1) then
          do ivar=ivar1,ivar2
            df_tmp_yz=df(l1-1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(df_tmp_yz,-deltay)
            df(l1-1,m1:m2,n1:n2,ivar)=df_tmp_yz
            df_tmp_yz=df(l2+1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(df_tmp_yz,+deltay)
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
      real, dimension (mx,my,mz,mfarray) :: f
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
        if (nxgrid>1 .and. lshear) then
          do ivar=ivar1,ivar2
            f_tmp_yz=f(l1-1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(f_tmp_yz,-deltay)
            f(l1-1,m1:m2,n1:n2,ivar)=f_tmp_yz
            f_tmp_yz=f(l2+1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(f_tmp_yz,+deltay)
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
endmodule GhostFold
