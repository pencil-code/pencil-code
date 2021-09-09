! $Id$
!
!  This module folds ghost zones for a multiple processor run.
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
  use Messages
  use Mpicomm
!
  implicit none
!
  private
!
  public :: fold_df, fold_f, fold_df_3points,reverse_fold_f_3points
  public :: reverse_fold_df_3points
!
  contains
!***********************************************************************
    subroutine fold_df(df,ivar1,ivar2)
!
!  Fold first ghost zone of df into main part of df.
!
!  20-may-2006/anders: coded
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: ivar1, ivar2
!
      real, dimension (nx+2,ny+2,1,ivar2-ivar1+1) :: df_tmp_xy
      real, dimension (nx+2,1,nz,ivar2-ivar1+1) :: df_tmp_xz
      real, dimension (1,ny,nz,ivar2-ivar1+1) :: df_tmp_yz
      integer :: nvar_fold, iproc_rcv, ivar
      integer :: itag1=10, itag2=11, itag3=12, itag4=13, itag5=14, itag6=15
!
      nvar_fold=ivar2-ivar1+1
!
!  The first ghost zone in the z-direction is folded (the corners will find
!  their proper positions after folding in x and y as well).
!
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
              call mpirecv_real(df_tmp_xy,(/nx+2,ny+2,1,nvar_fold/),zlneigh,itag1)
              df(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2)= &
                  df(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2) + df_tmp_xy
            elseif (iproc_rcv==zuneigh) then
              df_tmp_xy = df(l1-1:l2+1,m1-1:m2+1,n2+1:n2+1,ivar1:ivar2)
              call mpisend_real(df_tmp_xy,(/nx+2,ny+2,1,nvar_fold/),zuneigh,itag1)
            endif
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xy,(/nx+2,ny+2,1,nvar_fold/),zuneigh,itag2)
              df(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2)= &
                  df(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2) + df_tmp_xy
            elseif (iproc_rcv==zlneigh) then
              df_tmp_xy = df(l1-1:l2+1,m1-1:m2+1,n1-1:n1-1,ivar1:ivar2)
              call mpisend_real(df_tmp_xy,(/nx+2,ny+2,1,nvar_fold/),zlneigh,itag2)
            endif
          enddo
        endif
        df(l1-1:l2+1,m1-1:m2+1,n1-1,ivar1:ivar2)=0.0
        df(l1-1:l2+1,m1-1:m2+1,n2+1,ivar1:ivar2)=0.0
      endif
!
!  Then y.
!
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
              call mpirecv_real(df_tmp_xz,(/nx+2,1,nz,nvar_fold/),ylneigh,itag3)
              df(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2)= &
                  df(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2) + df_tmp_xz
            elseif (iproc_rcv==yuneigh) then
              df_tmp_xz = df(l1-1:l2+1,m2+1:m2+1,n1:n2,ivar1:ivar2)
              call mpisend_real(df_tmp_xz,(/nx+2,1,nz,nvar_fold/),yuneigh,itag3)
            endif
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xz,(/nx+2,1,nz,nvar_fold/),yuneigh,itag4)
              df(l1-1:l2+1,m2:m2,n1:n2,ivar1:ivar2)= &
                  df(l1-1:l2+1,m2:m2,n1:n2,ivar1:ivar2) + df_tmp_xz
            elseif (iproc_rcv==ylneigh) then
              df_tmp_xz = df(l1-1:l2+1,m1-1:m1-1,n1:n2,ivar1:ivar2)
              call mpisend_real(df_tmp_xz,(/nx+2,1,nz,nvar_fold/),ylneigh,itag4)
            endif
          enddo
!
        endif
!
!  With shearing boundary conditions we must take care that the information is
!  shifted properly before the final fold.
!
        if (lshear .and. nxgrid > 1 .and. nygrid > 1) call yshift_ghost(df, 1, ivar1, ivar2)
!
        df(l1-1:l2+1,m1-1,n1:n2,ivar1:ivar2)=0.0
        df(l1-1:l2+1,m2+1,n1:n2,ivar1:ivar2)=0.0
      endif
!
!  Finally x.
!
      if (nxgrid/=1) then
        if (nprocx==1) then
          df(l1,m1:m2,n1:n2,ivar1:ivar2)=df(l1,m1:m2,n1:n2,ivar1:ivar2) + &
              df(l2+1,m1:m2,n1:n2,ivar1:ivar2)
          df(l2,m1:m2,n1:n2,ivar1:ivar2)=df(l2,m1:m2,n1:n2,ivar1:ivar2) + &
              df(l1-1,m1:m2,n1:n2,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_yz,(/1,ny,nz,nvar_fold/),xlneigh,itag5)
              df(l1:l1,m1:m2,n1:n2,ivar1:ivar2)= &
                  df(l1:l1,m1:m2,n1:n2,ivar1:ivar2) + df_tmp_yz
            elseif (iproc_rcv==xuneigh) then
              df_tmp_yz = df(l2+1:l2+1,m1:m2,n1:n2,ivar1:ivar2)
              call mpisend_real(df_tmp_yz,(/1,ny,nz,nvar_fold/),xuneigh,itag5)
            endif
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_yz,(/1,ny,nz,nvar_fold/),xuneigh,itag6)
              df(l2:l2,m1:m2,n1:n2,ivar1:ivar2)= &
                  df(l2:l2,m1:m2,n1:n2,ivar1:ivar2) + df_tmp_yz
            elseif (iproc_rcv==xlneigh) then
              df_tmp_yz = df(l1-1:l1-1,m1:m2,n1:n2,ivar1:ivar2)
              call mpisend_real(df_tmp_yz,(/1,ny,nz,nvar_fold/),xlneigh,itag6)
            endif
          enddo
        endif
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
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
!
      real, dimension (nx+2,ny+2,1,ivar2-ivar1+1) :: f_tmp_xy
      real, dimension (nx+2,1,nz,ivar2-ivar1+1)   :: f_tmp_xz
      real, dimension (1,ny,nz,ivar2-ivar1+1) :: f_tmp_yz
      integer :: nvar_fold, iproc_rcv, ivar
      integer :: itag1=10, itag2=11, itag3=12, itag4=13, itag5=14, itag6=15
!
      nvar_fold=ivar2-ivar1+1
!
!  The first ghost zone in the z-direction is folded (the corners will find
!  their proper positions after folding in x and y as well).
!
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
!
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xy,(/nx+2,ny+2,1,nvar_fold/),zlneigh,itag1)
              f(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2)= &
                  f(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2) + f_tmp_xy
            elseif (iproc_rcv==zuneigh) then
              f_tmp_xy = f(l1-1:l2+1,m1-1:m2+1,n2+1:n2+1,ivar1:ivar2)
              call mpisend_real(f_tmp_xy,(/nx+2,ny+2,1,nvar_fold/),zuneigh,itag1)
            endif
!
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xy,(/nx+2,ny+2,1,nvar_fold/),zuneigh,itag2)
              f(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2)= &
                  f(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2) + f_tmp_xy
            elseif (iproc_rcv==zlneigh) then
              f_tmp_xy = f(l1-1:l2+1,m1-1:m2+1,n1-1:n1-1,ivar1:ivar2)
              call mpisend_real(f_tmp_xy,(/nx+2,ny+2,1,nvar_fold/),zlneigh,itag2)
            endif
!
          enddo
        endif
        f(l1-1:l2+1,m1-1:m2+1,n1-1,ivar1:ivar2)=0.0
        f(l1-1:l2+1,m1-1:m2+1,n2+1,ivar1:ivar2)=0.0
      endif
!
!  Then y.
!
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
              call mpirecv_real(f_tmp_xz,(/nx+2,1,nz,nvar_fold/),ylneigh,itag3)
              f(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2)= &
                  f(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2) + f_tmp_xz
            elseif (iproc_rcv==yuneigh) then
              f_tmp_xz = f(l1-1:l2+1,m2+1:m2+1,n1:n2,ivar1:ivar2)
              call mpisend_real(f_tmp_xz,(/nx+2,1,nz,nvar_fold/),yuneigh,itag3)
            endif
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xz,(/nx+2,1,nz,nvar_fold/),yuneigh,itag4)
              f(l1-1:l2+1,m2:m2,n1:n2,ivar1:ivar2)= &
                  f(l1-1:l2+1,m2:m2,n1:n2,ivar1:ivar2) + f_tmp_xz
            elseif (iproc_rcv==ylneigh) then
              f_tmp_xz = f(l1-1:l2+1,m1-1:m1-1,n1:n2,ivar1:ivar2)
              call mpisend_real(f_tmp_xz,(/nx+2,1,nz,nvar_fold/),ylneigh,itag4)
            endif
          enddo
!
        endif
!
!  With shearing boundary conditions we must take care that the information is
!  shifted properly before the final fold.
!
        if (lshear .and. nxgrid > 1 .and. nygrid > 1) call yshift_ghost(f, 1, ivar1, ivar2)
!
        f(l1-1:l2+1,m1-1,n1:n2,ivar1:ivar2)=0.0
        f(l1-1:l2+1,m2+1,n1:n2,ivar1:ivar2)=0.0
      endif
!
!  Finally x.
!
      if (nxgrid/=1) then
        if (nprocx==1) then
          f(l1,m1:m2,n1:n2,ivar1:ivar2)=f(l1,m1:m2,n1:n2,ivar1:ivar2) + &
              f(l2+1,m1:m2,n1:n2,ivar1:ivar2)
          f(l2,m1:m2,n1:n2,ivar1:ivar2)=f(l2,m1:m2,n1:n2,ivar1:ivar2) + &
              f(l1-1,m1:m2,n1:n2,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_yz,(/1,ny,nz,nvar_fold/),xlneigh,itag5)
              f(l1:l1,m1:m2,n1:n2,ivar1:ivar2)= &
                  f(l1:l1,m1:m2,n1:n2,ivar1:ivar2) + f_tmp_yz
            elseif (iproc_rcv==xuneigh) then
              f_tmp_yz = f(l2+1:l2+1,m1:m2,n1:n2,ivar1:ivar2)
              call mpisend_real(f_tmp_yz,(/1,ny,nz,nvar_fold/),xuneigh,itag5)
            endif
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_yz,(/1,ny,nz,nvar_fold/),xuneigh,itag6)
              f(l2:l2,m1:m2,n1:n2,ivar1:ivar2)= &
                  f(l2:l2,m1:m2,n1:n2,ivar1:ivar2) + f_tmp_yz
            elseif (iproc_rcv==xlneigh) then
              f_tmp_yz = f(l1-1:l1-1,m1:m2,n1:n2,ivar1:ivar2)
              call mpisend_real(f_tmp_yz,(/1,ny,nz,nvar_fold/),xlneigh,itag6)
            endif
          enddo
        endif
        f(l1-1,m1:m2,n1:n2,ivar1:ivar2)=0.0
        f(l2+1,m1:m2,n1:n2,ivar1:ivar2)=0.0
      endif
!
    endsubroutine fold_f
!***********************************************************************    
subroutine fold_df_3points(df,ivar1,ivar2)
!
!  Folds the ghost zones of df into main part of df. 
!  Used for particles runs where the gaussian box approach is used,
!  as this influences grid points up to 3 nodes away from the particle
!
!  20-may-2006/anders: coded
!  11-jan-2015/jonas : adapted for 3 node thick folded zones
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: ivar1, ivar2
!
      real, dimension (nx+6,ny+6,3,ivar2-ivar1+1) :: df_tmp_xy
      real, dimension (nx+6,3,nz,ivar2-ivar1+1) :: df_tmp_xz
      real, dimension (3,ny,nz,ivar2-ivar1+1) :: df_tmp_yz
      real, dimension (ny,nz) :: df_tmp_yz_one
      integer :: nvar_fold, iproc_rcv, ivar
      integer :: itag1=10, itag2=11, itag3=12, itag4=13, itag5=14, itag6=15
      integer :: l1m3, l2p3, m1m1, m2p3
!
      nvar_fold=ivar2-ivar1+1
!
!  The first ghost zone in the z-direction is folded (the corners will find
!  their proper positions after folding in x and y as well).
!
      l1m3=l1-3; l2p3=l2+3; 
      m1m1=m1-1; m2p3=m2+3     
      if (nzgrid/=1) then
        if (nprocz==1) then
! Folding the top ghost zones on the bottom
          df(l1m3:l2p3,m1-3:m2+3,n1:n1+2,ivar1:ivar2)= &
              df(l1-3:l2+3,m1-3:m2+3,n1:n1+2,ivar1:ivar2)+ &
              df(l1-3:l2+3,m1-3:m2+3,n2+1:n2+3,ivar1:ivar2)
! Folding the bottom ghost zones on the top
          df(l1-3:l2+3,m1-3:m2+3,n2-2:n2,ivar1:ivar2)= &
              df(l1-3:l2+3,m1-3:m2+3,n2-2:n2,ivar1:ivar2)+ &
              df(l1-3:l2+3,m1-3:m2+3,n1-3:n1-1,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xy,(/nx+6,ny+6,3,nvar_fold/),zlneigh,itag1)
              df(l1-3:l2+3,m1-3:m2+3,n1:n1+2,ivar1:ivar2)= &
                 df(l1-3:l2+3,m1-3:m2+3,n1:n1+2,ivar1:ivar2) + df_tmp_xy
            elseif (iproc_rcv==zuneigh) then
              df_tmp_xy = df(l1-3:l2+3,m1-3:m2+3,n2+1:n2+3,ivar1:ivar2)
              call mpisend_real(df_tmp_xy,(/nx+6,ny+6,3,nvar_fold/),zuneigh,itag1)
            endif
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xy,(/nx+6,ny+6,3,nvar_fold/),zuneigh,itag2)
              df(l1-3:l2+3,m1-3:m2+3,n2-2:n2,ivar1:ivar2)= &
                  df(l1-3:l2+3,m1-3:m2+3,n2-2:n2,ivar1:ivar2) + df_tmp_xy
            elseif (iproc_rcv==zlneigh) then
              df_tmp_xy = df(l1-3:l2+3,m1-3:m2+3,n1-3:n1-1,ivar1:ivar2)
              call mpisend_real(df_tmp_xy,(/nx+6,ny+6,3,nvar_fold/),zlneigh,itag2)
            endif
          enddo
        endif
        df(l1-3:l2+3,m1-3:m2+3,n1-3:n1-1,ivar1:ivar2)=0.0
        df(l1-3:l2+3,m1-3:m2+3,n2+1:n2+3,ivar1:ivar2)=0.0
      endif
!
!  Then y.
!
      if (nygrid/=1) then
        if (nprocy==1) then
          df(l1-3:l2+3,m1:m1+2,n1:n2,ivar1:ivar2)= &
              df(l1-3:l2+3,m1:m1+2,n1:n2,ivar1:ivar2) + &
              df(l1-3:l2+3,m2+1:m2+3,n1:n2,ivar1:ivar2)
          df(l1-3:l2+3,m2-2:m2,n1:n2,ivar1:ivar2)= &
              df(l1-3:l2+3,m2-2:m2,n1:n2,ivar1:ivar2) + &
              df(l1-3:l2+3,m1-3:m1-1,n1:n2,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xz,(/nx+6,3,nz,nvar_fold/),ylneigh,itag3)
              df(l1-3:l2+3,m1:m1+2,n1:n2,ivar1:ivar2)= &
                  df(l1-3:l2+3,m1:m1+2,n1:n2,ivar1:ivar2) + df_tmp_xz
            elseif (iproc_rcv==yuneigh) then
              df_tmp_xz = df(l1-3:l2+3,m2+1:m2+3,n1:n2,ivar1:ivar2)
              call mpisend_real(df_tmp_xz,(/nx+6,3,nz,nvar_fold/),yuneigh,itag3)
            endif
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xz,(/nx+6,3,nz,nvar_fold/),yuneigh,itag4)
              df(l1-3:l2+3,m2-2:m2,n1:n2,ivar1:ivar2)= &
                  df(l1-3:l2+3,m2-2:m2,n1:n2,ivar1:ivar2) + df_tmp_xz
            elseif (iproc_rcv==ylneigh) then
              df_tmp_xz = df(l1-3:l2+3,m1-3:m1-1,n1:n2,ivar1:ivar2)
              call mpisend_real(df_tmp_xz,(/nx+6,3,nz,nvar_fold/),ylneigh,itag4)
            endif
          enddo
        endif
        df(l1-3:l2+3,m2+1:m2+3,n1:n2,ivar1:ivar2)=0.0
        df(l1-3:l2+3,m1-3:m1-1,n1:n2,ivar1:ivar2)=0.0
      endif
!
!  Finally x.
!
      if (nxgrid/=1) then
        if (nprocx==1) then
          df(l1:l1+2,m1:m2,n1:n2,ivar1:ivar2)=&
              df(l1:l1+2,m1:m2,n1:n2,ivar1:ivar2) + &
              df(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)
          df(l2-2:l2,m1:m2,n1:n2,ivar1:ivar2)=&
              df(l2-2:l2,m1:m2,n1:n2,ivar1:ivar2) + &
              df(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_yz,(/3,ny,nz,nvar_fold/),xlneigh,itag5)
              df(l1:l1+2,m1:m2,n1:n2,ivar1:ivar2)= &
                  df(l1:l1+2,m1:m2,n1:n2,ivar1:ivar2) + df_tmp_yz
            elseif (iproc_rcv==xuneigh) then
              df_tmp_yz = df(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)
              call mpisend_real(df_tmp_yz,(/3,ny,nz,nvar_fold/),xuneigh,itag5)
            endif
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_yz,(/3,ny,nz,nvar_fold/),xuneigh,itag6)
              df(l2-2:l2,m1:m2,n1:n2,ivar1:ivar2)= &
                  df(l2-2:l2,m1:m2,n1:n2,ivar1:ivar2) + df_tmp_yz
            elseif (iproc_rcv==xlneigh) then
              df_tmp_yz = df(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)
              call mpisend_real(df_tmp_yz,(/3,ny,nz,nvar_fold/),xlneigh,itag6)
            endif
          enddo
        endif
        df(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)=0.0
        df(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)=0.0
      endif
!
    endsubroutine fold_df_3points
!*******************************************************************************
    subroutine yshift_ghost(f, ng, ivar1, ivar2)
!
!  Shift the ghost zones in y.
!
!  07-apr-16/ccyang: coded
!
!  Input Arguments
!    f
!        4D array of dimensions (mx,my,mz,:).
!    ng
!        Number of ghost zones to be operated on.
!    ivar1, ivar2
!        Range of the components in f to be operated on.
!
        real, dimension(:,:,:,:), intent(inout) :: f
        integer, intent(in) :: ng, ivar1, ivar2
!
        real, dimension(ng,ny,nz) :: work
        integer :: ivar
!
!  Shift the left boundary by -deltay.
!
        first: if (lfirst_proc_x) then
          comp1: do ivar = ivar1, ivar2
            work = f(l1-ng:l1-1,m1:m2,n1:n2,ivar)
            call yshift_block(ng, work, -deltay)
            f(l1-ng:l1-1,m1:m2,n1:n2,ivar) = work
          enddo comp1
        endif first
!
!  Shift the right boundary by +deltay.
!
        last: if (llast_proc_x) then
          comp2: do ivar = ivar1, ivar2
            work = f(l2+1:l2+ng,m1:m2,n1:n2,ivar)
            call yshift_block(ng, work, +deltay)
            f(l2+1:l2+ng,m1:m2,n1:n2,ivar) = work
          enddo comp2
        endif last
!
    endsubroutine yshift_ghost
!*******************************************************************************
    subroutine yshift_block(ng, a, shift)
!
!  Wrapper for shifting a block of data a in y by shift.
!
!  07-apr-16/ccyang: coded.
!
      integer, intent(in) :: ng
      real, dimension(ng,ny,nz), intent(inout) :: a
      real, intent(in) :: shift
!
!  If lghostfold_usebspline is set .true., use B-spline interpolation; use
!  Fourier interpolation, otherwise.
!
      if (lghostfold_usebspline) then
        call yshift_block_bspline(ng, a, shift)
      else
        call yshift_block_fft(ng, a, shift)
      endif
!
    endsubroutine yshift_block
!*******************************************************************************
    subroutine yshift_block_fft(ng, a, shift)
!
!  Shift a block of data a in y by shift using the Fourier interpolation.
!
!  07-apr-16/ccyang: coded.
!
      use Fourier, only: fourier_shift_yz_y
!
      integer, intent(in) :: ng
      real, dimension(ng,ny,nz), intent(inout) :: a
      real, intent(in) :: shift
!
      real, dimension(ny,nz) :: work
      integer :: i
!
      fft: do i = 1, ng
        work = a(i,:,:)
        call fourier_shift_yz_y(work, shift)
        a(i,:,:) = work
      enddo fft
!
    endsubroutine yshift_block_fft
!*******************************************************************************
    subroutine yshift_block_bspline(ng, a, shift)
!
!  Shift a block of data a in y by shift using the B-spline interpolation.
!
!  07-apr-16/ccyang: coded.
!
      use Mpicomm, only: remap_to_pencil_y, unmap_from_pencil_y
      use Sub, only: bspline_precondition, bspline_interpolation, ludcmp
!
      integer, intent(in) :: ng
      real, dimension(ng,ny,nz), intent(inout) :: a
      real, intent(in) :: shift
!
      integer, parameter :: bspline_k = 7
      real, dimension(nygrid,nygrid) :: bspline_ay = 0.0
      integer, dimension(nygrid) :: bspline_iy = 0
      logical :: lfirstcall = .true.
!
      real, dimension(ng,nygrid,nz) :: work
      real, dimension(nygrid) :: penc
      integer :: i, k
      real :: s
!
!  Precondition the B-spline.
!
      precon: if (lfirstcall) then
        call bspline_precondition(nygrid, bspline_k, bspline_ay)
        call ludcmp(bspline_ay, bspline_iy)
        lfirstcall = .false.
      endif precon
!
!  Make the interpolation.
!
      s = shift / dy
      call remap_to_pencil_y(a, work)
      zscan: do k = 1, nz
        xscan: do i = 1, ng
          penc = work(i,:,k)
          call bspline_interpolation(nygrid, bspline_k, penc, bspline_ay, bspline_iy, s)
          work(i,:,k) = penc
        enddo xscan
      enddo zscan
      call unmap_from_pencil_y(work, a)
!
    endsubroutine yshift_block_bspline
!*******************************************************************************
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
      real, dimension (nx,ny,3,ivar2-ivar1+1) :: f_tmp_xy
      real, dimension (nx,3,nz,ivar2-ivar1+1) :: f_tmp_xz
      real, dimension (3,ny,nz,ivar2-ivar1+1) :: f_tmp_yz
      real, dimension (ny,nz) :: f_tmp_yz_one
      integer :: nvar_fold, iproc_rcv, ivar
      integer :: itag1=10, itag2=11, itag3=12, itag4=13, itag5=14, itag6=15
      integer :: l1m1, l1m3, l2p1, l2p3
      integer :: m1m1, m1m3, m2p1, m2p3
      integer :: n1m1, n1m3, n2p1, n2p3
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
!
!
      l1m1=l1-1; l1m3=l1-3; l2p1=l2+1; l2p3=l2+3
      m1m1=m1-1; m1m3=m1-3; m2p1=m2+1; m2p3=m2+3
      n1m1=n1-1; n1m3=n1-3; n2p1=n2+1; n2p3=n2+3
      if (nzgrid/=1) then
        if (nprocz==1) then
! Folding the top ghost zones on the bottom
          f(l1:l2,m1:m2,n1m3:n1m1,ivar1:ivar2)= &
!              f(l1:l2,m1:m2,n1-3:n1-1,ivar1:ivar2)+ &
              f(l1:l2,m1:m2,n2i:n2,ivar1:ivar2)
! Folding the bottom ghost zones on the top
          f(l1:l2,m1:m2,n2p1:n2p3,ivar1:ivar2)= &
!              f(l1:l2,m1:m2,n2+1:n2+3,ivar1:ivar2)+ &
              f(l1:l2,m1:m2,n1:n1i,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xy,(/nx,ny,3,nvar_fold/),zlneigh,itag1)
              f(l1:l2,m1:m2,n1m3:n1m1,ivar1:ivar2)= f_tmp_xy
!                  f(l1:l2,m1:m2,n1-3:n1-1,ivar1:ivar2) + f_tmp_xy
            elseif (iproc_rcv==zuneigh) then
              f_tmp_xy = f(l1:l2,m1:m2,n2i:n2,ivar1:ivar2)
              call mpisend_real(f_tmp_xy,(/nx,ny,3,nvar_fold/),zuneigh,itag1)
            endif
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xy,(/nx,ny,3,nvar_fold/),zuneigh,itag2)
              f(l1:l2,m1:m2,n2+1:n2+3,ivar1:ivar2)= f_tmp_xy 
!                 f(l1:l2,m1:m2,n2+1:n2+3,ivar1:ivar2) + f_tmp_xy 
            elseif (iproc_rcv==zlneigh) then
              f_tmp_xy = f(l1:l2,m1:m2,n1:n1i,ivar1:ivar2)
              call mpisend_real(f_tmp_xy,(/nx,ny,3,nvar_fold/),zlneigh,itag2)
            endif
          enddo
        endif
!        f(l1-3:l2+3,m1-3:m2+3,n1-3:n1-1,ivar1:ivar2)=0.0
!        f(l1-3:l2+3,m1-3:m2+3,n2+1:n2+3,ivar1:ivar2)=0.0
      endif
!
!  Then y.
!
!
      if (nygrid/=1) then
        if (nprocy==1) then
          f(l1:l2,m1m3:m1m1,n1:n2,ivar1:ivar2)= &
!              f(l1:l2,m1-3:m1-1,n1:n2,ivar1:ivar2) + &
              f(l1:l2,m2i:m2,n1:n2,ivar1:ivar2)
          f(l1:l2,m2p1:m2p3,n1:n2,ivar1:ivar2)= &
!              f(l1:l2,m2+1:m2+3,n1:n2,ivar1:ivar2) + &
              f(l1:l2,m1:m1i,n1:n2,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xz,(/nx,3,nz,nvar_fold/),ylneigh,itag3)
              f(l1:l2,m1m3:m1m1,n1:n2,ivar1:ivar2)= &
!                  f(l1:l2,m1-3:m1-1,n1:n2,ivar1:ivar2) + f_tmp_xz
                  f_tmp_xz
            elseif (iproc_rcv==yuneigh) then
              f_tmp_xz = f(l1:l2,m2i:m2,n1:n2,ivar1:ivar2)
              call mpisend_real(f_tmp_xz,(/nx,3,nz,nvar_fold/),yuneigh,itag3)
            endif
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xz,(/nx,3,nz,nvar_fold/),yuneigh,itag4)
              f(l1:l2,m2+1:m2+3,n1:n2,ivar1:ivar2)= &
!                  f(l1:l2,m2+1:m2+3,n1:n2,ivar1:ivar2) + f_tmp_xz
                  f_tmp_xz
            elseif (iproc_rcv==ylneigh) then
              f_tmp_xz = f(l1:l2,m1:m1i,n1:n2,ivar1:ivar2)
              call mpisend_real(f_tmp_xz,(/nx,3,nz,nvar_fold/),ylneigh,itag4)
            endif
          enddo
        endif
      endif
!
!  Finally x.
!
!
      if (nxgrid/=1) then
        if (nprocx==1) then
          f(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)= &
!              f(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2) + &
              f(l2i:l2,m1:m2,n1:n2,ivar1:ivar2)
          f(l2p1:l2p3,m1:m2,n1:n2,ivar1:ivar2)= &
!              f(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2) + &
              f(l1:l1i,m1:m2,n1:n2,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_yz, &
                  (/3,ny,nz,nvar_fold/),xlneigh,itag5)
            elseif (iproc_rcv==xuneigh) then
              call mpisend_real(f(l2i:l2,m1:m2,n1:n2,ivar1:ivar2), &
                  (/3,ny,nz,nvar_fold/),xuneigh,itag5)
            endif
            if (iproc==iproc_rcv) f(l1m3:l1m1,m1:m2,n1:n2,ivar1:ivar2)= &
!                f(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2) + f_tmp_yz
                f_tmp_yz
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_yz, &
                  (/3,ny,nz,nvar_fold/),xuneigh,itag6)
            elseif (iproc_rcv==xlneigh) then
              call mpisend_real(f(l1:l1i,m1:m2,n1:n2,ivar1:ivar2), &
                  (/3,ny,nz,nvar_fold/),xlneigh,itag6)
            endif
            if (iproc==iproc_rcv) f(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)= &
!                f(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2) + f_tmp_yz
                f_tmp_yz
          enddo
        endif
!        f(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)=0.0
!        f(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)=0.0
      endif
!
!      print*, 'inside folding'
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
      real, dimension (nx,ny,3,ivar2-ivar1+1) :: f_tmp_xy
      real, dimension (nx,3,nz,ivar2-ivar1+1) :: f_tmp_xz
      real, dimension (3,ny,nz,ivar2-ivar1+1) :: f_tmp_yz
      real, dimension (ny,nz) :: f_tmp_yz_one
      integer :: nvar_fold, iproc_rcv, ivar
      integer :: itag1=10, itag2=11, itag3=12, itag4=13, itag5=14, itag6=15
      integer :: l1m1, l1m3, l2p1, l2p3
      integer :: m2m1, m2m3, m2p1, m2p3
      integer :: n1m1, n1m3, n2m1, n2m3, n2p1, n2p3
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
      m2m1=m2-1; m2m3=m2-3; m2p1=m2+1; m2p3=m2+3
      n1m1=n1-1; n1m3=n1-3; n2m1=n2-1; n2m3=n2-3; n2p1=n2+1; n2p3=n2+3
      if (nzgrid/=1) then
        if (nprocz==1) then
! Folding the top ghost zones on the bottom
          f(l1:l2,m1:m2,n1m3:n1m1,ivar1:ivar2)= &
!              f(l1:l2,m1:m2,n1-3:n1-1,ivar1:ivar2)+ &
              f(l1:l2,m1:m2,n2i:n2,ivar1:ivar2)
! Folding the bottom ghost zones on the top
          !f(l1:l2,m1:m2,n2+1:n2+3,ivar1:ivar2)= &
          f(l1:l2,m1:m2,n2p1:n2p3,ivar1:ivar2)= &
!              f(l1:l2,m1:m2,n2+1:n2+3,ivar1:ivar2)+ &
              f(l1:l2,m1:m2,n1:n1i,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xy, &
                  (/nx,ny,3,nvar_fold/),zlneigh,itag1)
            elseif (iproc_rcv==zuneigh) then
              call mpisend_real(f(l1:l2,m1:m2,n2i:n2,ivar1:ivar2), &
                  (/nx,ny,3,nvar_fold/),zuneigh,itag1)
            endif
            if (iproc==iproc_rcv) f(l1:l2,m1:m2,n1-3:n1-1,ivar1:ivar2)= f_tmp_xy
!                f(l1:l2,m1:m2,n1-3:n1-1,ivar1:ivar2) + f_tmp_xy
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xy, &
                  (/nx,ny,3,nvar_fold/),zuneigh,itag2)
            elseif (iproc_rcv==zlneigh) then
              call mpisend_real(f(l1:l2,m1:m2,n1:n1i,ivar1:ivar2), &
                  (/nx,ny,3,nvar_fold/),zlneigh,itag2)
            endif
            if (iproc==iproc_rcv) f(l1:l2,m1:m2,n2+1:n2+3,ivar1:ivar2)= f_tmp_xy 
!                f(l1:l2,m1:m2,n2+1:n2+3,ivar1:ivar2) + f_tmp_xy 
          enddo
        endif
!        f(l1-3:l2+3,m1-3:m2+3,n1-3:n1-1,ivar1:ivar2)=0.0
!        f(l1-3:l2+3,m1-3:m2+3,n2+1:n2+3,ivar1:ivar2)=0.0
      endif
!
!  Then y.
!
!
      if (nygrid/=1) then
        if (nprocy==1) then
          f(l1:l2,m1-3:m1-1,n1:n2,ivar1:ivar2)= &
!              f(l1:l2,m1-3:m1-1,n1:n2,ivar1:ivar2) + &
              f(l1:l2,m2i:m2,n1:n2,ivar1:ivar2)
          !f(l1:l2,m2+1:m2+3,n1:n2,ivar1:ivar2)= &
          f(l1:l2,m2p1:m2p3,n1:n2,ivar1:ivar2)= &
!              f(l1:l2,m2+1:m2+3,n1:n2,ivar1:ivar2) + &
              f(l1:l2,m1:m1i,n1:n2,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xz, &
                  (/nx,3,nz,nvar_fold/),ylneigh,itag3)
            elseif (iproc_rcv==yuneigh) then
              call mpisend_real(f(l1:l2,m2i:m2,n1:n2,ivar1:ivar2), &
                  (/nx,3,nz,nvar_fold/),yuneigh,itag3)
            endif
            if (iproc==iproc_rcv) f(l1:l2,m1-3:m1-1,n1:n2,ivar1:ivar2)= &
!                f(l1:l2,m1-3:m1-1,n1:n2,ivar1:ivar2) + f_tmp_xz
                f_tmp_xz
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xz, &
                  (/nx,3,nz,nvar_fold/),yuneigh,itag4)
            elseif (iproc_rcv==ylneigh) then
              call mpisend_real(f(l1:l2,m1:m1i,n1:n2,ivar1:ivar2), &
                  (/nx,3,nz,nvar_fold/),ylneigh,itag4)
            endif
            if (iproc==iproc_rcv) f(l1:l2,m2+1:m2+3,n1:n2,ivar1:ivar2)= &
!                f(l1:l2,m2+1:m2+3,n1:n2,ivar1:ivar2) + f_tmp_xz
                f_tmp_xz
          enddo
        endif
      endif
!
!  Finally x.
!
!
      l1m1=l1-1
      l1m3=l1-3
      l2p1=l2+1
      l2p3=l2+3
      if (nxgrid/=1) then
        if (nprocx==1) then
          f(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)= &
!              f(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2) + &
              f(l2i:l2,m1:m2,n1:n2,ivar1:ivar2)
          f(l2p1:l2p3,m1:m2,n1:n2,ivar1:ivar2)= &
!              f(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2) + &
              f(l1:l1i,m1:m2,n1:n2,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_yz, &
                  (/3,ny,nz,nvar_fold/),xlneigh,itag5)
            elseif (iproc_rcv==xuneigh) then
              call mpisend_real(f(l2i:l2,m1:m2,n1:n2,ivar1:ivar2), &
                  (/3,ny,nz,nvar_fold/),xuneigh,itag5)
            endif
            if (iproc==iproc_rcv) f(l1m3:l1m1,m1:m2,n1:n2,ivar1:ivar2)= &
!                f(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2) + f_tmp_yz
                f_tmp_yz
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_yz, &
                  (/3,ny,nz,nvar_fold/),xuneigh,itag6)
            elseif (iproc_rcv==xlneigh) then
              call mpisend_real(f(l1:l1i,m1:m2,n1:n2,ivar1:ivar2), &
                  (/3,ny,nz,nvar_fold/),xlneigh,itag6)
            endif
            if (iproc==iproc_rcv) f(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)= &
!                f(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2) + f_tmp_yz
                f_tmp_yz
          enddo
        endif
!        f(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)=0.0
!        f(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)=0.0
      endif
!
!      print*, 'inside folding'
!
    endsubroutine reverse_fold_df_3points
!*******************************************************************************
endmodule GhostFold
