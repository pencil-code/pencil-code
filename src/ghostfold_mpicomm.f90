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
  public :: fold_df, fold_f, fold_df_3points
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
      real, dimension (ny,nz) :: df_tmp_yz_one
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
              call mpirecv_real(df_tmp_xy, &
                  (/nx+2,ny+2,1,nvar_fold/),zlneigh,itag1)
            elseif (iproc_rcv==zuneigh) then
              call mpisend_real(df(l1-1:l2+1,m1-1:m2+1,n2+1:n2+1,ivar1:ivar2), &
                  (/nx+2,ny+2,1,nvar_fold/),zuneigh,itag1)
            endif
            if (iproc==iproc_rcv) df(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2)= &
                df(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2) + df_tmp_xy
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xy, &
                  (/nx+2,ny+2,1,nvar_fold/),zuneigh,itag2)
            elseif (iproc_rcv==zlneigh) then
              call mpisend_real(df(l1-1:l2+1,m1-1:m2+1,n1-1:n1-1,ivar1:ivar2), &
                  (/nx+2,ny+2,1,nvar_fold/),zlneigh,itag2)
            endif
            if (iproc==iproc_rcv) df(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2)= &
                df(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2) + df_tmp_xy
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
              call mpirecv_real(df_tmp_xz, &
                  (/nx+2,1,nz,nvar_fold/),ylneigh,itag3)
            elseif (iproc_rcv==yuneigh) then
              call mpisend_real(df(l1-1:l2+1,m2+1:m2+1,n1:n2,ivar1:ivar2), &
                  (/nx+2,1,nz,nvar_fold/),yuneigh,itag3)
            endif
            if (iproc==iproc_rcv) df(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2)= &
                df(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2) + df_tmp_xz
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xz, &
                  (/nx+2,1,nz,nvar_fold/),yuneigh,itag4)
            elseif (iproc_rcv==ylneigh) then
              call mpisend_real(df(l1-1:l2+1,m1-1:m1-1,n1:n2,ivar1:ivar2), &
                  (/nx+2,1,nz,nvar_fold/),ylneigh,itag4)
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
        if (nxgrid>1 .and. lshear .and. lfirst_proc_x) then
          do ivar=ivar1,ivar2
            df_tmp_yz_one=df(l1-1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(df_tmp_yz_one,-deltay)
            df(l1-1,m1:m2,n1:n2,ivar)=df_tmp_yz_one
          enddo
        endif
        if (nxgrid>1 .and. lshear .and. llast_proc_x) then
          do ivar=ivar1,ivar2
            df_tmp_yz_one=df(l2+1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(df_tmp_yz_one,+deltay)
            df(l2+1,m1:m2,n1:n2,ivar)=df_tmp_yz_one
          enddo
        endif
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
              call mpirecv_real(df_tmp_yz, &
                  (/1,ny,nz,nvar_fold/),xlneigh,itag5)
            elseif (iproc_rcv==xuneigh) then
              call mpisend_real(df(l2+1:l2+1,m1:m2,n1:n2,ivar1:ivar2), &
                  (/1,ny,nz,nvar_fold/),xuneigh,itag5)
            endif
            if (iproc==iproc_rcv) df(l1:l1,m1:m2,n1:n2,ivar1:ivar2)= &
                df(l1:l1,m1:m2,n1:n2,ivar1:ivar2) + df_tmp_yz
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_yz, &
                  (/1,ny,nz,nvar_fold/),xuneigh,itag6)
            elseif (iproc_rcv==xlneigh) then
              call mpisend_real(df(l1-1:l1-1,m1:m2,n1:n2,ivar1:ivar2), &
                  (/1,ny,nz,nvar_fold/),xlneigh,itag6)
            endif
            if (iproc==iproc_rcv) df(l2:l2,m1:m2,n1:n2,ivar1:ivar2)= &
                df(l2:l2,m1:m2,n1:n2,ivar1:ivar2) + df_tmp_yz
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
      real, dimension (ny,nz) :: f_tmp_yz_one
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
              call mpirecv_real(f_tmp_xy, &
                  (/nx+2,ny+2,1,nvar_fold/),zlneigh,itag1)
            elseif (iproc_rcv==zuneigh) then
              call mpisend_real(f(l1-1:l2+1,m1-1:m2+1,n2+1:n2+1,ivar1:ivar2), &
                  (/nx+2,ny+2,1,nvar_fold/),zuneigh,itag1)
            endif
!
            if (iproc==iproc_rcv) f(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2)= &
                f(l1-1:l2+1,m1-1:m2+1,n1:n1,ivar1:ivar2) + f_tmp_xy
!
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xy, &
                  (/nx+2,ny+2,1,nvar_fold/),zuneigh,itag2)
            elseif (iproc_rcv==zlneigh) then
              call mpisend_real(f(l1-1:l2+1,m1-1:m2+1,n1-1:n1-1,ivar1:ivar2), &
                  (/nx+2,ny+2,1,nvar_fold/),zlneigh,itag2)
            endif
!
            if (iproc==iproc_rcv) f(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2)= &
                f(l1-1:l2+1,m1-1:m2+1,n2:n2,ivar1:ivar2) + f_tmp_xy
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
              call mpirecv_real(f_tmp_xz, &
                  (/nx+2,1,nz,nvar_fold/),ylneigh,itag3)
            elseif (iproc_rcv==yuneigh) then
              call mpisend_real(f(l1-1:l2+1,m2+1:m2+1,n1:n2,ivar1:ivar2), &
                  (/nx+2,1,nz,nvar_fold/),yuneigh,itag3)
            endif
            if (iproc==iproc_rcv) f(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2)= &
                f(l1-1:l2+1,m1:m1,n1:n2,ivar1:ivar2) + f_tmp_xz
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_xz, &
                  (/nx+2,1,nz,nvar_fold/),yuneigh,itag4)
            elseif (iproc_rcv==ylneigh) then
              call mpisend_real(f(l1-1:l2+1,m1-1:m1-1,n1:n2,ivar1:ivar2), &
                  (/nx+2,1,nz,nvar_fold/),ylneigh,itag4)
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
        if (nxgrid>1 .and. lshear .and. lfirst_proc_x) then
          do ivar=ivar1,ivar2
            f_tmp_yz_one=f(l1-1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(f_tmp_yz_one,-deltay)
            f(l1-1,m1:m2,n1:n2,ivar)=f_tmp_yz_one
          enddo
        endif
        if (nxgrid>1 .and. lshear .and. llast_proc_x) then
          do ivar=ivar1,ivar2
            f_tmp_yz_one=f(l2+1,m1:m2,n1:n2,ivar)
            call fourier_shift_yz_y(f_tmp_yz_one,+deltay)
            f(l2+1,m1:m2,n1:n2,ivar)=f_tmp_yz_one
          enddo
        endif
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
              call mpirecv_real(f_tmp_yz, &
                  (/1,ny,nz,nvar_fold/),xlneigh,itag5)
            elseif (iproc_rcv==xuneigh) then
              call mpisend_real(f(l2+1:l2+1,m1:m2,n1:n2,ivar1:ivar2), &
                  (/1,ny,nz,nvar_fold/),xuneigh,itag5)
            endif
            if (iproc==iproc_rcv) f(l1:l1,m1:m2,n1:n2,ivar1:ivar2)= &
                f(l1:l1,m1:m2,n1:n2,ivar1:ivar2) + f_tmp_yz
            if (iproc==iproc_rcv) then
              call mpirecv_real(f_tmp_yz, &
                  (/1,ny,nz,nvar_fold/),xuneigh,itag6)
            elseif (iproc_rcv==xlneigh) then
              call mpisend_real(f(l1-1:l1-1,m1:m2,n1:n2,ivar1:ivar2), &
                  (/1,ny,nz,nvar_fold/),xlneigh,itag6)
            endif
            if (iproc==iproc_rcv) f(l2:l2,m1:m2,n1:n2,ivar1:ivar2)= &
                f(l2:l2,m1:m2,n1:n2,ivar1:ivar2) + f_tmp_yz
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
!
      nvar_fold=ivar2-ivar1+1
!
!  The first ghost zone in the z-direction is folded (the corners will find
!  their proper positions after folding in x and y as well).
!
      if (nzgrid/=1) then
        if (nprocz==1) then
! Folding the top ghost zones on the bottom
          df(l1-3:l2+3,m1-3:m2+3,n1:n1+2,ivar1:ivar2)= &
              df(l1-3:l2+3,m1-3:m2+3,n1:n1+2,ivar1:ivar2)+ &
              df(l1-3:l2+3,m1-3:m2+3,n2+1:n2+3,ivar1:ivar2)
! Folding the bottom ghost zones on the top
          df(l1-3:l2+3,m1-3:m2+3,n2-2:n2,ivar1:ivar2)= &
              df(l1-3:l2+3,m1-3:m2+3,n2-2:n2,ivar1:ivar2)+ &
              df(l1-3:l2+3,m1-3:m2+3,n1-3:n1-1,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xy, &
                  (/nx+6,ny+6,3,nvar_fold/),zlneigh,itag1)
            elseif (iproc_rcv==zuneigh) then
              call mpisend_real(df(l1-3:l2+3,m1-3:m2+3,n2+1:n2+3,ivar1:ivar2), &
                  (/nx+6,ny+6,3,nvar_fold/),zuneigh,itag1)
            endif
            if (iproc==iproc_rcv) df(l1-3:l2+3,m1-3:m2+3,n1:n1+2,ivar1:ivar2)= &
                df(l1-3:l2+3,m1-3:m2+3,n1:n1+2,ivar1:ivar2) + df_tmp_xy
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xy, &
                  (/nx+6,ny+6,3,nvar_fold/),zuneigh,itag2)
            elseif (iproc_rcv==zlneigh) then
              call mpisend_real(df(l1-3:l2+3,m1-3:m2+3,n1-3:n1-1,ivar1:ivar2), &
                  (/nx+6,ny+6,3,nvar_fold/),zlneigh,itag2)
            endif
            if (iproc==iproc_rcv) df(l1-3:l2+3,m1-3:m2+3,n2-2:n2,ivar1:ivar2)= &
                df(l1-3:l2+3,m1-3:m2+3,n2-2:n2,ivar1:ivar2) + df_tmp_xy 
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
              call mpirecv_real(df_tmp_xz, &
                  (/nx+6,3,nz,nvar_fold/),ylneigh,itag3)
            elseif (iproc_rcv==yuneigh) then
              call mpisend_real(df(l1-3:l2+3,m2+1:m2+3,n1:n2,ivar1:ivar2), &
                  (/nx+6,3,nz,nvar_fold/),yuneigh,itag3)
            endif
            if (iproc==iproc_rcv) df(l1-3:l2+3,m1:m1+2,n1:n2,ivar1:ivar2)= &
                df(l1-3:l2+3,m1:m1+2,n1:n2,ivar1:ivar2) + df_tmp_xz
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_xz, &
                  (/nx+6,3,nz,nvar_fold/),yuneigh,itag4)
            elseif (iproc_rcv==ylneigh) then
              call mpisend_real(df(l1-3:l2+3,m1-3:m1-1,n1:n2,ivar1:ivar2), &
                  (/nx+6,3,nz,nvar_fold/),ylneigh,itag4)
            endif
            if (iproc==iproc_rcv) df(l1-3:l2+3,m2-2:m2,n1:n2,ivar1:ivar2)= &
                df(l1-3:l2+3,m2-2:m2,n1:n2,ivar1:ivar2) + df_tmp_xz
          enddo
        endif
      endif
!
!  Finally x.
!
      if (nxgrid/=1) then
        if (nprocx==1) then
          df(l1:l1+2,m1:m2,n1:n2,ivar1:ivar2)=df(l1:l1+2,m1:m2,n1:n2,ivar1:ivar2) + &
              df(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)
          df(l2-2:l2,m1:m2,n1:n2,ivar1:ivar2)=df(l2-2:l2,m1:m2,n1:n2,ivar1:ivar2) + &
              df(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)
        else
          do iproc_rcv=0,ncpus-1
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_yz, &
                  (/3,ny,nz,nvar_fold/),xlneigh,itag5)
            elseif (iproc_rcv==xuneigh) then
              call mpisend_real(df(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2), &
                  (/3,ny,nz,nvar_fold/),xuneigh,itag5)
            endif
            if (iproc==iproc_rcv) df(l1:l1+2,m1:m2,n1:n2,ivar1:ivar2)= &
                df(l1:l1+2,m1:m2,n1:n2,ivar1:ivar2) + df_tmp_yz
            if (iproc==iproc_rcv) then
              call mpirecv_real(df_tmp_yz, &
                  (/3,ny,nz,nvar_fold/),xuneigh,itag6)
            elseif (iproc_rcv==xlneigh) then
              call mpisend_real(df(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2), &
                  (/3,ny,nz,nvar_fold/),xlneigh,itag6)
            endif
            if (iproc==iproc_rcv) df(l2-2:l2,m1:m2,n1:n2,ivar1:ivar2)= &
                df(l2-2:l2,m1:m2,n1:n2,ivar1:ivar2) + df_tmp_yz
          enddo
        endif
        df(l1-3:l1-1,m1:m2,n1:n2,ivar1:ivar2)=0.0
        df(l2+1:l2+3,m1:m2,n1:n2,ivar1:ivar2)=0.0
      endif
!
    endsubroutine fold_df_3points
!*******************************************************************************
endmodule GhostFold
