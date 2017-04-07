
!***********************************************************************
!  FROM MPICOMM.F90
!***********************************************************************
!  ROUTINES
!    initiate_isendrcv_bdry_ogrid
!    finalize_isendrcv_bdry_ogrid
!    isendrcv_bdry_x_ogrid
!***********************************************************************

module Solid_Cells_Mpicomm

  use Cparam

  implicit none

  include 'solid_cells_mpi.h'

  integer, parameter :: nx_ogrid=nxgrid_ogrid/nprocx,ny_ogrid=nygrid_ogrid/nprocy,nz_ogrid=nzgrid_ogrid/nprocz
  integer, parameter :: mx_ogrid=nx_ogrid+2*nghost,l1_ogrid=1+nghost,l2_ogrid=mx_ogrid-nghost
  integer, parameter :: my_ogrid=ny_ogrid+2*nghost,m1_ogrid=1+nghost,m2_ogrid=my_ogrid-nghost
  integer, parameter :: mz_ogrid=nz_ogrid+2*nghost,n1_ogrid=1+nghost,n2_ogrid=mz_ogrid-nghost

  real, dimension (nghost,ny_ogrid,nz_ogrid,mcom) :: lbufxi,ubufxi,lbufxo,ubufxo
  real, dimension (:,:,:,:), allocatable :: lbufyi,ubufyi,lbufyo,ubufyo
  real, dimension (:,:,:,:), allocatable :: lbufzi,ubufzi,lbufzo,ubufzo
  real, dimension (:,:,:,:), allocatable :: llbufi,lubufi,uubufi,ulbufi
  real, dimension (:,:,:,:), allocatable :: llbufo,lubufo,uubufo,ulbufo
  contains 

!***********************************************************************
    subroutine initiate_isendrcv_bdry_ogrid(f)
!
!  Isend and Irecv boundary values. Called in the beginning of pde.
!  Does not wait for the receives to finish (done in finalize_isendrcv_bdry)
!  leftneigh and rightneigh are initialized by mpicomm_init.
!
!  07-feb-17/Jorgen: adapted from mpicomm.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
!
      integer :: tolowyr,touppyr,tolowys,touppys,tolowzr,touppzr,tolowzs,touppzs ! msg. tags placeholders
      integer :: TOllr,TOulr,TOuur,TOlur,TOlls,TOuls,TOuus,TOlus                 ! placeholder tags
      integer :: ivar1, ivar2, nbufy, nbufz, nbufyz, mxl, comm, bufact, dir
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (ivar2==0) return
!
      mxl=size(f,1)
!
!  Set communication across x-planes.
!
      if (nprocx>1) call isendrcv_bdry_x_ogrid(f)
!
!  Set communication across y-planes.
!  Standard message tags from defaults for surfaces and yz-corners.
!
      tolowyr=tolowy; tolowzr=tolowz; TOlls=TOll; TOllr=TOll
      tolowys=tolowy; tolowzs=tolowz; TOlus=TOlu; TOlur=TOlu
      touppyr=touppy; touppzr=touppz; TOuus=TOuu; TOuur=TOuu
      touppys=touppy; touppzs=touppz; TOuls=TOul; TOulr=TOul
!
!  Allocate and send/receive buffers across y-planes
!
      bufact=mxl*nghost*(ivar2-ivar1+1)

      if (nprocy>1) then
!
!  Internal, periodic exchange y-plane buffers.
!
        lbufyo(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m1i_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2)         !(lower y-zone)
        comm=MPI_COMM_GRID

        nbufy=bufact*bufsizes_yz(INYL,IRCV)
        call MPI_IRECV(lbufyi(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
                       ylneigh,touppyr,comm,irecv_rq_fromlowy,mpierr)

        nbufy=bufact*bufsizes_yz(INYL,ISND)
        call MPI_ISEND(lbufyo(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
                       ylneigh,tolowys,comm,isend_rq_tolowy,mpierr)

        ubufyo(:,:,:,ivar1:ivar2)=f(:,m2i_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(upper y-zone)
        comm=MPI_COMM_GRID
!
        nbufy=bufact*bufsizes_yz(INYU,IRCV)
        call MPI_IRECV(ubufyi(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
                       yuneigh,tolowyr,comm,irecv_rq_fromuppy,mpierr)

        nbufy=bufact*bufsizes_yz(INYU,ISND)
        call MPI_ISEND(ubufyo(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
                       yuneigh,touppys,comm,isend_rq_touppy,mpierr)
      endif
!
!  Set communication across z-planes.
!
      if (nprocz>1) then
!        
        lbufzo(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m2_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2) !lower z-planes
        comm=MPI_COMM_GRID

        nbufz=bufact*bufsizes_yz(INZL,IRCV)
        call MPI_IRECV(lbufzi(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
                       zlneigh,touppzr,comm,irecv_rq_fromlowz,mpierr)

        nbufz=bufact*bufsizes_yz(INZL,ISND)
        call MPI_ISEND(lbufzo(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
                       zlneigh,tolowz,comm,isend_rq_tolowz,mpierr)
!
        ubufzo(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m2_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2) !upper z-planes
        comm=MPI_COMM_GRID

        nbufz=bufact*bufsizes_yz(INZU,IRCV)
        call MPI_IRECV(ubufzi(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
                       zuneigh,tolowzr,comm,irecv_rq_fromuppz,mpierr)

        nbufz=bufact*bufsizes_yz(INZU,ISND)
        call MPI_ISEND(ubufzo(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
                       zuneigh,touppz,comm,isend_rq_touppz,mpierr)

      endif
!
!  The four corners (in counter-clockwise order).
!  (NOTE: this should work even for nprocx>1)
!
      if ((nprocz>1).and.(nprocy>1)) then
!
!  Internal and periodic yz-corner buffers.
!
        bufact=mxl*(ivar2-ivar1+1)
!
!  Lower y, lower z.
!
        llbufo(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m1i_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2)
        comm=MPI_COMM_GRID
!
        nbufyz=bufact*product(bufsizes_yz_corn(:,INLL,IRCV))
        if (llcornr>=0) call MPI_IRECV(llbufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       llcornr,TOuur,comm,irecv_rq_FRll,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INLL,ISND))
        if (llcorns>=0) call MPI_ISEND(llbufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       llcorns,TOlls,comm,isend_rq_TOll,mpierr)
!
!  Upper y, lower z.
!
        ulbufo(:,:,:,ivar1:ivar2)=f(:,m2i_ogrid:m2_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2)
        comm=MPI_COMM_GRID

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUL,IRCV))
        if (ulcornr>=0) call MPI_IRECV(ulbufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       ulcornr,TOlur,comm,irecv_rq_FRul,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUL,ISND))
        if (ulcorns>=0) call MPI_ISEND(ulbufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       ulcorns,TOuls,comm,isend_rq_TOul,mpierr)
!
!  Upper y, upper z.
!
        uubufo(:,:,:,ivar1:ivar2)=f(:,m2i_ogrid:m2_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2)
        comm=MPI_COMM_GRID

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUU,IRCV))
        if (uucornr>=0) call MPI_IRECV(uubufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       uucornr,TOllr,comm,irecv_rq_FRuu,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUU,ISND))
        if (uucorns>=0) call MPI_ISEND(uubufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       uucorns,TOuus,comm,isend_rq_TOuu,mpierr)
!
!  Lower y, upper z.
!
        lubufo(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m1i_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2)
        comm=MPI_COMM_GRID

        nbufyz=bufact*product(bufsizes_yz_corn(:,INLU,IRCV))
        if (lucornr>=0) call MPI_IRECV(lubufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       lucornr,TOulr,comm,irecv_rq_FRlu,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INLU,ISND))
        if (lucorns>=0) call MPI_ISEND(lubufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
                                       lucorns,TOlus,comm,isend_rq_TOlu,mpierr)
!
      endif
!
    endsubroutine initiate_isendrcv_bdry_ogrid
!***********************************************************************
    subroutine finalize_isendrcv_bdry_ogrid(f)
!
!  Make sure the communications initiated with initiate_isendrcv_bdry are
!  finished and insert the just received boundary values.
!  Receive requests do not need to (and on OSF1 cannot) be explicitly
!  freed, since MPI_Wait takes care of this.
!
!  07-feb-17/Jorgen: Adapted from mpicomm.f90
!
      use General, only: transpose_mn, notanumber

      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
!
      integer :: ivar1, ivar2, j
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (ivar2==0) return
!
!  1. wait until data received
!  2. set ghost zones
!  3. wait until send completed, will be overwritten in next time step
!
!  Communication across y-planes (includes periodic bc)
!
      if (nprocy>1) then
        call MPI_WAIT(irecv_rq_fromuppy,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowy,irecv_stat_fl,mpierr)
        do j=ivar1,ivar2
          if (.not. lfirst_proc_y .or. bcy12(j,1)=='p') then  
            f(:,1:m1_ogrid-1,n1_ogrid:n2_ogrid,j)=lbufyi(:,:,:,j)       ! set lower buffer
          endif

          if (.not. llast_proc_y .or. bcy12(j,2)=='p') then
            f(:,m2_ogrid+1:,n1_ogrid:n2_ogrid,j)=ubufyi(:,:,:,j)        ! set upper buffer
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowy,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppy,isend_stat_tu,mpierr)
      endif
!
!  Communication in z (includes periodic bc)
!
      if (nprocz>1) then
        call MPI_WAIT(irecv_rq_fromuppz,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowz,irecv_stat_fl,mpierr)
        do j=ivar1,ivar2
          if (.not. lfirst_proc_z .or. bcz12(j,1)=='p') then
            f(:,m1_ogrid:m2_ogrid,1:n1_ogrid-1,j)=lbufzi(:,:,:,j)       ! set lower buffer
          endif

          if (.not. llast_proc_z .or. bcz12(j,2)=='p') then 
            f(:,m1_ogrid:m2_ogrid,n2_ogrid+1:,j)=ubufzi(:,:,:,j)        ! set upper buffer
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowz,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppz,isend_stat_tu,mpierr)
      endif
!
!  The four yz-corners (in counter-clockwise order)
!
       if (nprocz>1.and.nprocy>1) then

        if (uucornr>=0) call MPI_WAIT(irecv_rq_FRuu,irecv_stat_Fuu,mpierr)
        if (lucornr>=0) call MPI_WAIT(irecv_rq_FRlu,irecv_stat_Flu,mpierr)
        if (llcornr>=0) call MPI_WAIT(irecv_rq_FRll,irecv_stat_Fll,mpierr)
        if (ulcornr>=0) call MPI_WAIT(irecv_rq_FRul,irecv_stat_Ful,mpierr)

        do j=ivar1,ivar2
!
!  Set ll corner
!
          if (.not. lfirst_proc_z .or. bcz12(j,1)=='p') then
            if  (.not.lfirst_proc_y.or.bcy12(j,1)=='p') then    ! inner or periodic proc boundaries
              f(:,1:m1_ogrid-1,1:n1_ogrid-1,j)=llbufi(:,:,:,j)               ! fill lower left corner
            endif
!
!  Set ul corner
!
            if (.not.llast_proc_y .or.bcy12(j,2)=='p') then    ! inner or periodic proc boundaries
                f(:,m2_ogrid+1:,1:n1_ogrid-1,j)=ulbufi(:,:,:,j)                ! fill lower right corner
            endif
          endif
!
!  Set uu corner
!
          if (.not. llast_proc_z .or. bcz12(j,2)=='p') then
            if (.not.llast_proc_y.or.bcy12(j,2)=='p') then    ! inner or periodic proc boundaries
              f(:,m2+1:,n2+1:,j)=uubufi(:,:,:,j)                ! fill upper right corner
            endif
!
!  Set lu corner
!
            if (.not. lfirst_proc_y .or. bcy12(j,1)=='p') then
              f(:,1:m1-1,n2+1:,j)=lubufi(:,:,:,j)                ! fill upper left corner
            endif
          endif
        enddo
        
        if (llcorns>=0) call MPI_WAIT(isend_rq_TOll,isend_stat_Tll,mpierr)
        if (ulcorns>=0) call MPI_WAIT(isend_rq_TOul,isend_stat_Tul,mpierr)
        if (uucorns>=0) call MPI_WAIT(isend_rq_TOuu,isend_stat_Tuu,mpierr)
        if (lucorns>=0) call MPI_WAIT(isend_rq_TOlu,isend_stat_Tlu,mpierr)

      endif
!
!  make sure the other processors don't carry on sending new data
!  which could be mistaken for an earlier time
!
    call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
!
    endsubroutine finalize_isendrcv_bdry_ogrid
!***********************************************************************
    subroutine isendrcv_bdry_x_ogrid(f)
!
!  Isend and Irecv boundary values for x-direction. Sends and receives
!  before continuing to y and z boundaries, as this allows the edges
!  of the grid to be set properly.
!
!  07-feb-17/Jorgen: Adapted from mpicomm.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
!
      integer :: ivar1, ivar2, nbufx, j
!
      ivar1=1; ivar2=min(mcom,size(f,4))
!
      if (nprocx>1) then

        lbufxo(:,:,:,ivar1:ivar2)=f(l1_ogrid:l1i_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(lower x-zone)
        ubufxo(:,:,:,ivar1:ivar2)=f(l2i_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(upper x-zone)
        nbufx=ny_ogrid*nz_ogrid*nghost*(ivar2-ivar1+1)

        call MPI_IRECV(ubufxi(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xuneigh,tolowx,MPI_COMM_GRID,irecv_rq_fromuppx,mpierr)
        call MPI_IRECV(lbufxi(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xlneigh,touppx,MPI_COMM_GRID,irecv_rq_fromlowx,mpierr)
        call MPI_ISEND(lbufxo(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xlneigh,tolowx,MPI_COMM_GRID,isend_rq_tolowx,mpierr)
        call MPI_ISEND(ubufxo(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xuneigh,touppx,MPI_COMM_GRID,isend_rq_touppx,mpierr)
        call MPI_WAIT(irecv_rq_fromuppx,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowx,irecv_stat_fl,mpierr)
!
!  Inner communication or (shear-)periodic boundary conditions in x
!  MR: Communication should only happen under these conditions.
!
        do j=ivar1,ivar2
          if (.not. lfirst_proc_x .or. bcx12(j,1)=='p' .or. &
              (bcx12(j,1)=='she'.and.nygrid==1)) then
            f( 1:l1_ogrid-1,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j)=lbufxi(:,:,:,j)  !!(set lower buffer)
          endif
          if (.not. llast_proc_x .or. bcx12(j,2)=='p' .or. &
              (bcx12(j,2)=='she'.and.nygrid==1)) then
            f(l2_ogrid+1:,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j)=ubufxi(:,:,:,j)  !!(set upper buffer)
          endif
        enddo

        call MPI_WAIT(isend_rq_tolowx,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppx,isend_stat_tu,mpierr)

      endif
!
    endsubroutine isendrcv_bdry_x_ogrid
!***********************************************************************

endmodule Solid_Cells_Mpicomm
