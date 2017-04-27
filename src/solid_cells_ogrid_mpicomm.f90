
!***********************************************************************
!  FROM MPICOMM.F90
!***********************************************************************
!  ROUTINES
!     initiate_isendrcv_bdry_ogrid
!     finalize_isendrcv_bdry_ogrid
!     isendrcv_bdry_x_ogrid
!     initialize_mpicomm_ogrid
!     finalize_isend_init_interpol
!***********************************************************************

module Solid_Cells_Mpicomm

  use Cparam
  use Cdata
  use Mpicomm, only: mpisend_nonblock_real, mpirecv_nonblock_real, mpiwait, mpibarrier

  implicit none

  include 'solid_cells_mpi.h'
  ! Grid parameters
  integer, parameter :: nx_ogrid=nxgrid_ogrid/nprocx,ny_ogrid=nygrid_ogrid/nprocy,nz_ogrid=nzgrid_ogrid/nprocz
  integer, parameter :: mx_ogrid=nx_ogrid+2*nghost,l1_ogrid=1+nghost,l2_ogrid=mx_ogrid-nghost
  integer, parameter :: my_ogrid=ny_ogrid+2*nghost,m1_ogrid=1+nghost,m2_ogrid=my_ogrid-nghost
  integer, parameter :: mz_ogrid=nz_ogrid+2*nghost,n1_ogrid=1+nghost,n2_ogrid=mz_ogrid-nghost
  integer, parameter :: l1i_ogrid=l1_ogrid+nghost-1
  integer, parameter :: l2i_ogrid=mx_ogrid-2*nghost+1
  integer, parameter :: m1i_ogrid=m1_ogrid+nghost-1
  integer, parameter :: m2i_ogrid=my_ogrid-2*nghost+1
  integer, parameter :: n1i_ogrid=n1_ogrid+nghost-1
  integer, parameter :: n2i_ogrid=mz_ogrid-2*nghost+1

  ! MPI variables and parameters
  real, dimension (nghost,ny_ogrid,nz_ogrid,mcom) :: lbufxi_og,ubufxi_og,lbufxo_og,ubufxo_og
  integer :: isend_rq_tolowx,isend_rq_touppx,irecv_rq_fromlowx,irecv_rq_fromuppx

  real, dimension (:,:,:,:), allocatable :: lbufyi_og,ubufyi_og,lbufyo_og,ubufyo_og
  real, dimension (:,:,:,:), allocatable :: lbufzi_og,ubufzi_og,lbufzo_og,ubufzo_og
  real, dimension (:,:,:,:), allocatable :: llbufi_og,lubufi_og,uubufi_og,ulbufi_og
  real, dimension (:,:,:,:), allocatable :: llbufo_og,lubufo_og,uubufo_og,ulbufo_og
  integer, dimension(4,2) :: bufsizes_yz_og
  integer, parameter :: INYU=1, INZU=2, INYL=3, INZL=4
  integer, parameter :: INUU=1, INLU=2, INLL=3, INUL=4
  integer, parameter :: IRCV=1, ISND=2
  integer :: isend_rq_tolowy,isend_rq_touppy,irecv_rq_fromlowy,irecv_rq_fromuppy
  integer :: isend_rq_tolowz,isend_rq_touppz,irecv_rq_fromlowz,irecv_rq_fromuppz
  integer :: isend_rq_TOll,isend_rq_TOul,isend_rq_TOuu,isend_rq_TOlu  !(corners)
  integer :: irecv_rq_FRuu,irecv_rq_FRlu,irecv_rq_FRll,irecv_rq_FRul  !(corners)
 
  logical :: lcorner_yz=.false.
  integer :: llcorn,lucorn,uucorn,ulcorn            ! (the 4 corners in yz-plane)
  integer, parameter :: tolowx=13,touppx=14,tolowy=3,touppy=4,tolowz=5,touppz=6 ! msg. tags
  integer, parameter :: TOll=7,TOul=8,TOuu=9,TOlu=10 ! msg. tags for corners

  contains 

!***********************************************************************
    subroutine initialize_mpicomm_ogrid
!
!  Initialise MPI communication on the overlapping grids. This is run
!  after the standard initialize_mpicomm subroutine.
!
!  19-apr-17/Jorgen: Adapted and modified, from initialize_mpicomm in mpicomm.f90
!
!  Announce myself for pc_run to detect.
!
      if (lroot) print *, 'initialize_mpicomm_ogrid: enabled MPI on overlapping grid'
!
!  Set corners and neighbours. These are the same as those set in mpicomm.f90
!
!  Am I a yz corner?
!
      lcorner_yz=(lfirst_proc_y.or.llast_proc_y).and.(lfirst_proc_z.or.llast_proc_z)
!
!  Set up `lower' and `upper' neighbours, refer to MPI_COMM_GRID.
!
      xlneigh = modulo(ipx-1,nprocx) + ipy*nprocx + ipz*nprocxy
      xuneigh = modulo(ipx+1,nprocx) + ipy*nprocx + ipz*nprocxy
      ylneigh = ipx + modulo(ipy-1,nprocy)*nprocx + ipz*nprocxy
      yuneigh = ipx + modulo(ipy+1,nprocy)*nprocx + ipz*nprocxy
      zlneigh = ipx + ipy*nprocx + modulo(ipz-1,nprocz)*nprocxy
      zuneigh = ipx + ipy*nprocx + modulo(ipz+1,nprocz)*nprocxy
!
!  Set the four corners in the yz-plane (in cyclic order).
!
      llcorn = ipx + modulo(ipy-1,nprocy)*nprocx + modulo(ipz-1,nprocz)*nprocxy
      ulcorn = ipx + modulo(ipy+1,nprocy)*nprocx + modulo(ipz-1,nprocz)*nprocxy
      uucorn = ipx + modulo(ipy+1,nprocy)*nprocx + modulo(ipz+1,nprocz)*nprocxy
      lucorn = ipx + modulo(ipy-1,nprocy)*nprocx + modulo(ipz+1,nprocz)*nprocxy
!
      bufsizes_yz_og(:,IRCV)=(/nz_ogrid,ny_ogrid,nz_ogrid,ny_ogrid/) 
      bufsizes_yz_og(:,ISND)=(/nz_ogrid,ny_ogrid,nz_ogrid,ny_ogrid/) 
!
      allocate( lbufyi_og(mx_ogrid,nghost,bufsizes_yz_og(INYL,IRCV),mcom), &
                ubufyi_og(mx_ogrid,nghost,bufsizes_yz_og(INYU,IRCV),mcom), &
                lbufyo_og(mx_ogrid,nghost,bufsizes_yz_og(INYL,ISND),mcom), &
                ubufyo_og(mx_ogrid,nghost,bufsizes_yz_og(INYU,ISND),mcom), &
                lbufzi_og(mx_ogrid,bufsizes_yz_og(INZL,IRCV),nghost,mcom), &
                ubufzi_og(mx_ogrid,bufsizes_yz_og(INZU,IRCV),nghost,mcom), &
                lbufzo_og(mx_ogrid,bufsizes_yz_og(INZL,ISND),nghost,mcom), &
                ubufzo_og(mx_ogrid,bufsizes_yz_og(INZU,ISND),nghost,mcom)   )
!  
      allocate( llbufi_og(mx_ogrid,nghost,nghost,mcom), &
                llbufo_og(mx_ogrid,nghost,nghost,mcom), &
                lubufi_og(mx_ogrid,nghost,nghost,mcom), &
                lubufo_og(mx_ogrid,nghost,nghost,mcom), &
                ulbufi_og(mx_ogrid,nghost,nghost,mcom), &
                ulbufo_og(mx_ogrid,nghost,nghost,mcom), &
                uubufi_og(mx_ogrid,nghost,nghost,mcom), &
                uubufo_og(mx_ogrid,nghost,nghost,mcom)   )
!
    endsubroutine initialize_mpicomm_ogrid
!***********************************************************************

    subroutine initiate_isendrcv_bdry_ogrid(f)
!
!  Isend and Irecv boundary values. Called in the beginning of pde.
!  Does not wait for the receives to finish (done in finalize_isendrcv_bdry)
!  leftneigh and rightneigh are initialized by mpicomm_init.
!
!  apr-17/Jorgen: adapted from mpicomm.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray) :: f
!
      integer, dimension(4) :: nbuf_y, nbuf_z, nbuf_yz
      integer, parameter :: ivar1=1, ivar2=min(mcom,size(f,4))
      integer, parameter :: sendvar=ivar2-ivar1+1
!
      if (ivar2==0) return
!
!  Set communication across x-planes.
!
      if (nprocx>1) call isendrcv_bdry_x_ogrid(f)
!
!  Allocate and send/receive buffers across y-planes
!
      if (nprocy>1) then
!
!  Internal, periodic exchange y-plane buffers.
!
        lbufyo_og(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m1i_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2)         !(lower y-zone)
        nbuf_y=(/mx_ogrid,nghost,bufsizes_yz_og(INYL,IRCV),sendvar/)
        call mpirecv_nonblock_real(lbufyi_og(:,:,:,ivar1:ivar2),nbuf_y,ylneigh,touppy,irecv_rq_fromlowy)
        call mpisend_nonblock_real(lbufyo_og(:,:,:,ivar1:ivar2),nbuf_y,ylneigh,tolowy,isend_rq_tolowy)

        ubufyo_og(:,:,:,ivar1:ivar2)=f(:,m2i_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(upper y-zone)
        nbuf_y=(/mx_ogrid,nghost,bufsizes_yz_og(INYU,IRCV),sendvar/)
        call mpirecv_nonblock_real(ubufyi_og(:,:,:,ivar1:ivar2),nbuf_y,yuneigh,tolowy,irecv_rq_fromuppy)
        call mpisend_nonblock_real(ubufyo_og(:,:,:,ivar1:ivar2),nbuf_y,yuneigh,touppy,isend_rq_touppy)
      endif
!
!  Set communication across z-planes.
!
      if (nprocz>1) then
!        
        lbufzo_og(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m2_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2) !lower z-planes
        nbuf_z=(/mx_ogrid,bufsizes_yz_og(INZL,IRCV),nghost,sendvar/)
        call mpirecv_nonblock_real(lbufzi_og(:,:,:,ivar1:ivar2),nbuf_z,zlneigh,touppz,irecv_rq_fromlowz)
        call mpisend_nonblock_real(lbufzo_og(:,:,:,ivar1:ivar2),nbuf_z,zlneigh,tolowz,isend_rq_tolowz)
!
        ubufzo_og(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m2_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2) !upper z-planes
        nbuf_z=(/mx_ogrid,bufsizes_yz_og(INZU,IRCV),nghost,sendvar/)
        call mpirecv_nonblock_real(ubufzi_og(:,:,:,ivar1:ivar2),nbuf_z,zuneigh,tolowz,irecv_rq_fromuppz)
        call mpisend_nonblock_real(ubufzo_og(:,:,:,ivar1:ivar2),nbuf_z,zuneigh,touppz,isend_rq_touppz)
      endif
!
!  The four corners (in counter-clockwise order).
!
      if ((nprocz>1).and.(nprocy>1)) then
!
!  Internal and periodic yz-corner buffers.
!
        nbuf_yz=(/mx_ogrid,nghost,nghost,sendvar/)
!
!  Lower y, lower z.
!
        llbufo_og(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m1i_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2)
        if (llcorn>=0) then
          call mpirecv_nonblock_real(llbufi_og(:,:,:,ivar1:ivar2),nbuf_yz,llcorn,TOuu,irecv_rq_FRll)
          call mpisend_nonblock_real(llbufo_og(:,:,:,ivar1:ivar2),nbuf_yz,llcorn,TOll,isend_rq_TOll)
        endif
!
!  Upper y, lower z.
!
        ulbufo_og(:,:,:,ivar1:ivar2)=f(:,m2i_ogrid:m2_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2)
        if (ulcorn>=0) then
          call mpirecv_nonblock_real(ulbufi_og(:,:,:,ivar1:ivar2),nbuf_yz,ulcorn,TOlu,irecv_rq_FRul)
          call mpisend_nonblock_real(ulbufo_og(:,:,:,ivar1:ivar2),nbuf_yz,ulcorn,TOul,isend_rq_TOul)
        endif
!
!  Upper y, upper z.
!
        uubufo_og(:,:,:,ivar1:ivar2)=f(:,m2i_ogrid:m2_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2)
        if (uucorn>=0) then
          call mpirecv_nonblock_real(uubufi_og(:,:,:,ivar1:ivar2),nbuf_yz,uucorn,TOll,irecv_rq_FRuu)
          call mpisend_nonblock_real(uubufo_og(:,:,:,ivar1:ivar2),nbuf_yz,uucorn,TOuu,isend_rq_TOuu)
        endif
!
!  Lower y, upper z.
!
        lubufo_og(:,:,:,ivar1:ivar2)=f(:,m1_ogrid:m1i_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2)
        if (lucorn>=0) then
          call mpirecv_nonblock_real(lubufi_og(:,:,:,ivar1:ivar2),nbuf_yz,lucorn,TOul,irecv_rq_FRlu)
          call mpisend_nonblock_real(lubufo_og(:,:,:,ivar1:ivar2),nbuf_yz,lucorn,TOlu,isend_rq_TOlu)
        endif
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
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mfarray), intent(inout) :: f
      integer, parameter :: ivar1=1, ivar2=min(mcom,size(f,4))
      integer :: j
!
      if (ivar2==0) return
!
!  1. wait until data received
!  2. set ghost zones
!  3. wait until send completed, will be overwritten in next time step
!
!  Communication across y-planes (includes periodic bc)
!  Note: Always periodic BC in theta-dir
!
      if (nprocy>1) then
        call mpiwait(irecv_rq_fromuppy)
        call mpiwait(irecv_rq_fromlowy)
        do j=ivar1,ivar2
          f(:,1:m1_ogrid-1,n1_ogrid:n2_ogrid,j)=lbufyi_og(:,:,:,j)       ! set lower buffer
          f(:,m2_ogrid+1:,n1_ogrid:n2_ogrid,j)=ubufyi_og(:,:,:,j)        ! set upper buffer
        enddo
        call mpiwait(isend_rq_tolowy)
        call mpiwait(isend_rq_touppy)
      endif
!
!  Communication in z (includes periodic bc)
!  Note:Assume same z-periodicity for cartesian and curvilinear grid
!
      if (nprocz>1) then
        call mpiwait(irecv_rq_fromuppz)
        call mpiwait(irecv_rq_fromlowz)
        if (.not. lfirst_proc_z .or. bcz12(j,1)=='p') then            
          do j=ivar1,ivar2
            f(:,m1_ogrid:m2_ogrid,1:n1_ogrid-1,j)=lbufzi_og(:,:,:,j)       ! set lower buffer
          enddo
        endif
        if (.not. llast_proc_z .or. bcz12(j,2)=='p') then 
          do j=ivar1,ivar2
            f(:,m1_ogrid:m2_ogrid,n2_ogrid+1:,j)=ubufzi_og(:,:,:,j)        ! set upper buffer
          enddo
        endif
        call mpiwait(isend_rq_tolowz)
        call mpiwait(isend_rq_touppz)
      endif
!
!  The four yz-corners (in counter-clockwise order)
!
       if (nprocz>1.and.nprocy>1) then
        if (uucorn>=0) call mpiwait(irecv_rq_FRuu)
        if (lucorn>=0) call mpiwait(irecv_rq_FRlu)
        if (llcorn>=0) call mpiwait(irecv_rq_FRll)
        if (ulcorn>=0) call mpiwait(irecv_rq_FRul)
        do j=ivar1,ivar2
          if (.not. lfirst_proc_z .or. bcz12(j,1)=='p') then
            f(:,1:m1_ogrid-1,1:n1_ogrid-1,j)=llbufi_og(:,:,:,j)               ! fill lower left corner (ll)
            f(:,m2_ogrid+1:,1:n1_ogrid-1,j)=ulbufi_og(:,:,:,j)                ! fill lower right corner (ul)
          endif
          if (.not. llast_proc_z .or. bcz12(j,2)=='p') then
            f(:,m2+1:,n2+1:,j)=uubufi_og(:,:,:,j)                             ! fill upper right corner (uu)
            f(:,1:m1-1,n2+1:,j)=lubufi_og(:,:,:,j)                            ! fill upper left corner (lu)
          endif
        enddo
        if (llcorn>=0) call mpiwait(isend_rq_TOll)
        if (ulcorn>=0) call mpiwait(isend_rq_TOul)
        if (uucorn>=0) call mpiwait(isend_rq_TOuu)
        if (lucorn>=0) call mpiwait(isend_rq_TOlu)
      endif
!
!  make sure the other processors don't carry on sending new data
!  which could be mistaken for an earlier time
!
      call mpibarrier
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
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mfarray), intent(inout) :: f
      integer, parameter :: ivar1=1, ivar2=min(mcom,size(f,4))
      integer :: j
      integer, dimension(4), parameter ::  nbuf_x=(/ny_ogrid,nz_ogrid,nghost,ivar2-ivar1+1/)
!
      lbufxo_og(:,:,:,ivar1:ivar2)=f(l1_ogrid:l1i_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(lower x-zone)
      ubufxo_og(:,:,:,ivar1:ivar2)=f(l2i_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(upper x-zone)
      call mpirecv_nonblock_real(ubufxi_og(:,:,:,ivar1:ivar2),nbuf_x,xuneigh,tolowx,irecv_rq_fromuppx)
      call mpirecv_nonblock_real(lbufxi_og(:,:,:,ivar1:ivar2),nbuf_x,xlneigh,touppx,irecv_rq_fromlowx)
      call mpisend_nonblock_real(lbufxo_og(:,:,:,ivar1:ivar2),nbuf_x,xlneigh,tolowx,isend_rq_tolowx)
      call mpisend_nonblock_real(ubufxo_og(:,:,:,ivar1:ivar2),nbuf_x,xuneigh,touppx,isend_rq_touppx)
      call mpiwait(irecv_rq_fromuppx)
      call mpiwait(irecv_rq_fromlowx)
!
!  Inner communication in x
!  Note: Never periodic BC for radial direction
!
      if (.not. lfirst_proc_x) then
        do j=ivar1,ivar2
          f( 1:l1_ogrid-1,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j)=lbufxi_og(:,:,:,j)  !!(set lower buffer)
        enddo
      endif
      if (.not. llast_proc_x) then
        do j=ivar1,ivar2
          f(l2_ogrid+1:,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j)=ubufxi_og(:,:,:,j)  !!(set upper buffer)
        enddo
      endif
      call mpiwait(isend_rq_tolowx)
      call mpiwait(isend_rq_touppx)
!
    endsubroutine isendrcv_bdry_x_ogrid
!***********************************************************************
    subroutine finalize_isend_init_interpol(ireq1D,ireq2D,nreq1D,nreq2D)
!
!  Wait for non-blocking communication in setting up interpolation arrays to finish
!
      integer :: nreq1D, nreq2D
      integer, dimension(nreq1D) :: ireq1D
      integer, dimension(nreq2D) :: ireq2D
      integer :: j

      do j=1,nreq1D
        call mpiwait(ireq1D(j))
      enddo
      do j=1,nreq2D
        call mpiwait(ireq2D(j))
      enddo
!
      call mpibarrier
!
    endsubroutine finalize_isend_init_interpol
!***********************************************************************

endmodule Solid_Cells_Mpicomm
