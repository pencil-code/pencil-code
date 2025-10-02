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
  use Fourier
!
  implicit none
!
  private
!
  public :: fold_df, fold_f,fold_df_3points,reverse_fold_f_3points
  public :: reverse_fold_df_3points
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
    subroutine fold_df_3points(df,ivar1,ivar2)
!
!  Fold whole ghost zone of df into main part of df.
!  This is used for when the gaussian scheme is used for
!  source term distribution, used when reactive particles
!  interact with the flow
!  15-may-2006/anders: coded
!  11-jan-2015/jonas: adapted to fold the whole ghost zone
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: ivar1,ivar2
!
!  Fold z-direction first (including all ghost zones in x and y).
!
      if (nzgrid/=1) then
        df(l1-3:l2+3,m1-3:m2+3,n1:n1+2,ivar1:ivar2)= &
            df(l1-3:l2+3,m1-3:m2+3,n1:n1+2,ivar1:ivar2) + &
            df(l1-3:l2+3,m1-3:m2+3,n2+1:n2+3,ivar1:ivar2)
        df(l1-3:l2+3,m1-3:m2+3,n2-2:n2,ivar1:ivar2)= &
            df(l1-3:l2+3,m1-3:m2+3,n2-2:n2,ivar1:ivar2) + &
            df(l1-3:l2+3,m1-3:m2+3,n1-3:n1-1,ivar1:ivar2)
        df(l1-3:l2+3,m1-3:m2+3,n1-3:n1-1,ivar1:ivar2)=0.0
        df(l1-3:l2+3,m1-3:m2+3,n2+1:n2+3,ivar1:ivar2)=0.0
      endif
!
!  Then fold y-direction (including all ghost zones in x).
!
      if (nygrid/=1) then
        df(l1-3:l2+3,m1:m1+2,n1:n2,ivar1:ivar2)= &
            df(l1-3:l2+3,m1:m1+2,n1:n2,ivar1:ivar2) + &
            df(l1-3:l2+3,m2+1:m2+3,n1:n2,ivar1:ivar2)
        df(l1-3:l2+3,m2-2:m2,n1:n2,ivar1:ivar2)= &
            df(l1-3:l2+3,m2-2:m2,n1:n2,ivar1:ivar2) + &
            df(l1-3:l2+3,m1-3:m1-1,n1:n2,ivar1:ivar2)
        df(l1-3:l2+3,m2+1:m2+3,n1:n2,ivar1:ivar2)=0.0
        df(l1-3:l2+3,m1-3:m1-1,n1:n2,ivar1:ivar2)=0.0
      endif
!
!  Finally fold the x-direction.
!
      if (nxgrid/=1) then
        df(l1:l1+2,m1:m2,n1:n2,ivar1:ivar2)=df(l1:l1+2,m1:m2,n1:n2,ivar1:ivar2) + &
            df(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)
        df(l2-2:l2,m1:m2,n1:n2,ivar1:ivar2)=df(l2-2:l2,m1:m2,n1:n2,ivar1:ivar2) + &
            df(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)
        df(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)=0.0
        df(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)=0.0
      endif
!
    endsubroutine fold_df_3points
!***********************************************************************
subroutine reverse_fold_f_3points(f,ivar1,ivar2)
!
!  Copies the boundary points (l2i:l2) into the ghost zones of the neighbouring
!  node
!  Used for reactive particle runs particle-fluid transfer is diffused before
!  added to the df-array,
!
!
!  20-may-2006/anders: coded
!  11-jan-2015/jonas : adapted for 3 node thick folded zones
!  22-aug-2016/jonas : adapted to reverse
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
!
      integer :: nvar_fold
!
      nvar_fold=ivar2-ivar1+1
!
!  The three ghost zones of the auxiliary array are filled (read: overwritten).
!  by values from the neighbouring domains
!  The ghost zones are communicated in following order (and coordinates:
!  x: (l2i:l2,m1:m2,n1:n2) is sent to (l1-3:l1-1,m1:m2,n1:n2) of the following processor
!     (l1:l1i,m1:m2,n1:n2) is sent to (l2+1:l2+3,m1:m2,n1:n2) of the preceeding processor
!  y: (l1:l2,m2i:m2,n1:n2) is sent to (l1:l2,m1-3:m1-1,n1:n2) of the following processor
!     (l1:l2,m1:m1i,n1:n2) is sent to (l1:l2,m2+1:m2+3,n1:n2) of the preceeding processor
!  z: (l1:l2,m1:m2,n2i:n2) is sent to (l1:l2,m1:m2,n1-3:n1-1) of the following processor
!     (l1:l2,m1:m2,n1:n1i) is sent to (l1:l2,m1:m2,n2+1:n2+3) of the preceeding processor
!
      if (nzgrid/=1) then
! Folding the top ghost zones on the bottom
        f(l1:l2,m1:m2,n1-3:n1-1,ivar1:ivar2)= &
            f(l1:l2,m1:m2,n1-3:n1-1,ivar1:ivar2)+ &
            f(l1:l2,m1:m2,n2i:n2,ivar1:ivar2)
! Folding the bottom ghost zones on the top
        f(l1:l2,m1:m2,n2+1:n2+3,ivar1:ivar2)= &
            f(l1:l2,m1:m2,n2+1:n2+3,ivar1:ivar2)+ &
            f(l1:l2,m1:m2,n1:n1i,ivar1:ivar2)
!        f(l1-3:l2+3,m1-3:m2+3,n1-3:n1-1,ivar1:ivar2)=0.0
!        f(l1-3:l2+3,m1-3:m2+3,n2+1:n2+3,ivar1:ivar2)=0.0
      endif
!
!  Then y.
!
      if (nygrid/=1) then
        f(l1:l2,m1-3:m1-1,n1:n2,ivar1:ivar2)= &
            f(l1:l2,m1-3:m1-1,n1:n2,ivar1:ivar2) + &
            f(l1:l2,m2i:m2,n1:n2,ivar1:ivar2)
        f(l1:l2,m2+1:m2+3,n1:n2,ivar1:ivar2)= &
            f(l1:l2,m2+1:m2+3,n1:n2,ivar1:ivar2) + &
            f(l1:l2,m1:m1i,n1:n2,ivar1:ivar2)
      endif
!
!  Finally x.
!
      if (nxgrid/=1) then
        f(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)= &
            f(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2) + &
            f(l2i:l2,m1:m2,n1:n2,ivar1:ivar2)
        f(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)= &
            f(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2) + &
            f(l1:l1i,m1:m2,n1:n2,ivar1:ivar2)
      endif
!
    endsubroutine reverse_fold_f_3points
!*******************************************************************************
subroutine reverse_fold_df_3points(f,ivar1,ivar2)
!
!  Copies the boundary points (l2i:l2) into the ghost zones of the neighbouring
!  node
!  Used for reactive particle runs particle-fluid transfer is diffused before
!  added to the df-array,
!
!
!  20-may-2006/anders: coded
!  11-jan-2015/jonas : adapted for 3 node thick folded zones
!  22-aug-2016/jonas : adapted to reverse
!
      real, dimension (mx,my,mz,mvar) :: f
      integer :: ivar1, ivar2
!
      integer :: nvar_fold
!
      nvar_fold=ivar2-ivar1+1
!
!  The three ghost zones of the auxiliary array are filled (read: overwritten).
!  by values from the neighbouring domains
!  The ghost zones are communicated in following order (and coordinates:
!  x: (l2i:l2,m1:m2,n1:n2) is sent to (l1-3:l1-1,m1:m2,n1:n2) of the following processor
!     (l1:l1i,m1:m2,n1:n2) is sent to (l2+1:l2+3,m1:m2,n1:n2) of the preceeding processor
!  y: (l1:l2,m2i:m2,n1:n2) is sent to (l1:l2,m1-3:m1-1,n1:n2) of the following processor
!     (l1:l2,m1:m1i,n1:n2) is sent to (l1:l2,m2+1:m2+3,n1:n2) of the preceeding processor
!  z: (l1:l2,m1:m2,n2i:n2) is sent to (l1:l2,m1:m2,n1-3:n1-1) of the following processor
!     (l1:l2,m1:m2,n1:n1i) is sent to (l1:l2,m1:m2,n2+1:n2+3) of the preceeding processor
!
      if (nzgrid/=1) then
! Folding the top ghost zones on the bottom
        f(l1:l2,m1:m2,n1-3:n1-1,ivar1:ivar2)= f(l1:l2,m1:m2,n2i:n2,ivar1:ivar2)
! Folding the bottom ghost zones on the top
        f(l1:l2,m1:m2,n2+1:n2+3,ivar1:ivar2)= f(l1:l2,m1:m2,n1:n1i,ivar1:ivar2)
!        f(l1-3:l2+3,m1-3:m2+3,n1-3:n1-1,ivar1:ivar2)=0.0
!        f(l1-3:l2+3,m1-3:m2+3,n2+1:n2+3,ivar1:ivar2)=0.0
      endif
!
!  Then y.
!
      if (nygrid/=1) then
        f(l1:l2,m1-3:m1-1,n1:n2,ivar1:ivar2)= f(l1:l2,m2i:m2,n1:n2,ivar1:ivar2)
        f(l1:l2,m2+1:m2+3,n1:n2,ivar1:ivar2)= f(l1:l2,m1:m1i,n1:n2,ivar1:ivar2)
      endif
!
!  Finally x.
!
      if (nxgrid/=1) then
        f(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)=f(l2i:l2,m1:m2,n1:n2,ivar1:ivar2)
        f(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)= f(l1:l1i,m1:m2,n1:n2,ivar1:ivar2)
      endif
!
    endsubroutine reverse_fold_df_3points
!*******************************************************************************
endmodule GhostFold
