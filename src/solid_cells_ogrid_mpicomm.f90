
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
  use Mpicomm, only: mpisend_nonblock_real, mpirecv_nonblock_real, mpiwait, mpibarrier, &
                     mpisend_real, mpirecv_real, mpibcast_real

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
  integer :: isend_rq_tolowx,isend_rq_touppx,irecv_rq_fromlowx,irecv_rq_fromuppx

  real, dimension (:,:,:,:), allocatable :: lbufxi_og,ubufxi_og,lbufxo_og,ubufxo_og
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

  ! Padé Filtering
  integer, parameter :: filter_Hsize = 10/2-nghost   ! Also set in solid_cells_ogrid.f90
  real, dimension (:,:,:,:), allocatable :: lbufyi_fi,ubufyi_fi,lbufyo_fi,ubufyo_fi
  real, dimension (:,:,:,:), allocatable :: lbufxi_fi,ubufxi_fi,lbufxo_fi,ubufxo_fi
  integer :: irecv_rq_fromlowy_fi,isend_rq_tolowy_fi,irecv_rq_fromuppy_fi,isend_rq_touppy_fi
  integer :: irecv_rq_fromlowx_fi,isend_rq_tolowx_fi,irecv_rq_fromuppx_fi,isend_rq_touppx_fi

  contains 

!***********************************************************************
    subroutine initialize_mpicomm_ogrid(lfilter_solution)
!
!  Initialise MPI communication on the overlapping grids. This is run
!  after the standard initialize_mpicomm subroutine.
!
!  19-apr-17/Jorgen: Adapted and modified, from initialize_mpicomm in mpicomm.f90
!  29-nov-17/Jorgen: Allocation of filter zones added
!
      logical, intent(in) :: lfilter_solution
!
!  Announce myself for pc_run to detect.
!
      if (lroot) print *, 'initialize_mpicomm_ogrid: enabled MPI on overlapping grid'
!
      if (nprocz>1 .and. lpole(2)) &
        call not_implemented('initialize_mpicomm_ogrid','communication across poles')

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
      if (nprocx>1) &
        allocate( lbufxi_og(nghost,ny_ogrid,nz_ogrid,mcom), &
                  ubufxi_og(nghost,ny_ogrid,nz_ogrid,mcom), &
                  lbufxo_og(nghost,ny_ogrid,nz_ogrid,mcom), &
                  ubufxo_og(nghost,ny_ogrid,nz_ogrid,mcom)   )

      if (nprocy>1) &
        allocate( lbufyi_og(mx_ogrid,nghost,bufsizes_yz_og(INYL,IRCV),mcom), &
                  ubufyi_og(mx_ogrid,nghost,bufsizes_yz_og(INYU,IRCV),mcom), &
                  lbufyo_og(mx_ogrid,nghost,bufsizes_yz_og(INYL,ISND),mcom), &
                  ubufyo_og(mx_ogrid,nghost,bufsizes_yz_og(INYU,ISND),mcom), &
                  lbufzi_og(mx_ogrid,bufsizes_yz_og(INZL,IRCV),nghost,mcom), &
                  ubufzi_og(mx_ogrid,bufsizes_yz_og(INZU,IRCV),nghost,mcom), &
                  lbufzo_og(mx_ogrid,bufsizes_yz_og(INZL,ISND),nghost,mcom), &
                  ubufzo_og(mx_ogrid,bufsizes_yz_og(INZU,ISND),nghost,mcom)   )
!  
      if (nprocy>1.and.nprocz>1) &
        allocate( llbufi_og(mx_ogrid,nghost,nghost,mcom), &
                  llbufo_og(mx_ogrid,nghost,nghost,mcom), &
                  lubufi_og(mx_ogrid,nghost,nghost,mcom), &
                  lubufo_og(mx_ogrid,nghost,nghost,mcom), &
                  ulbufi_og(mx_ogrid,nghost,nghost,mcom), &
                  ulbufo_og(mx_ogrid,nghost,nghost,mcom), &
                  uubufi_og(mx_ogrid,nghost,nghost,mcom), &
                  uubufo_og(mx_ogrid,nghost,nghost,mcom)   )
!
      if (lfilter_solution) then
        if (nprocy>1) &
          allocate( lbufyi_fi(nx_ogrid,filter_Hsize,nz_ogrid,mcom), &
                    ubufyi_fi(nx_ogrid,filter_Hsize,nz_ogrid,mcom), &
                    lbufyo_fi(nx_ogrid,filter_Hsize,nz_ogrid,mcom), &
                    ubufyo_fi(nx_ogrid,filter_Hsize,nz_ogrid,mcom)  )
        if (nprocx>1) &
          allocate( lbufxi_fi(filter_Hsize,ny_ogrid,nz_ogrid,mcom), &
                    ubufxi_fi(filter_Hsize,ny_ogrid,nz_ogrid,mcom), &
                    lbufxo_fi(filter_Hsize,ny_ogrid,nz_ogrid,mcom), &
                    ubufxo_fi(filter_Hsize,ny_ogrid,nz_ogrid,mcom)  )
      endif
!
    endsubroutine initialize_mpicomm_ogrid
!***********************************************************************
    subroutine initiate_isendrcv_bdry_ogrid(f_og)
!
!  Isend and Irecv boundary values. Called in the beginning of pde.
!  Does not wait for the receives to finish (done in finalize_isendrcv_bdry)
!  leftneigh and rightneigh are initialized by mpicomm_init.
!
!  apr-17/Jorgen: adapted from mpicomm.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar) :: f_og
!
      integer, dimension(4) :: nbuf_y, nbuf_z, nbuf_yz
      integer, parameter :: ivar1=1, ivar2=min(mcom,size(f_og,4))
      integer, parameter :: sendvar=ivar2-ivar1+1
!
      if (ivar2==0) return
!
!  Set communication across x-planes.
!
      if (nprocx>1) call isendrcv_bdry_x_ogrid(f_og)
      call mpibarrier
!
!  Allocate and send/receive buffers across y-planes
!
      if (nprocy>1) then
!
!  Internal, periodic exchange y-plane buffers.
!
        lbufyo_og(:,:,:,ivar1:ivar2)=f_og(:,m1_ogrid:m1i_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2)         !(lower y-zone)
        nbuf_y=(/mx_ogrid,nghost,bufsizes_yz_og(INYL,IRCV),sendvar/)
        call mpirecv_nonblock_real(lbufyi_og(:,:,:,ivar1:ivar2),nbuf_y,ylneigh,touppy,irecv_rq_fromlowy)
        call mpisend_nonblock_real(lbufyo_og(:,:,:,ivar1:ivar2),nbuf_y,ylneigh,tolowy,isend_rq_tolowy)

        ubufyo_og(:,:,:,ivar1:ivar2)=f_og(:,m2i_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(upper y-zone)
        nbuf_y=(/mx_ogrid,nghost,bufsizes_yz_og(INYU,IRCV),sendvar/)
        call mpirecv_nonblock_real(ubufyi_og(:,:,:,ivar1:ivar2),nbuf_y,yuneigh,tolowy,irecv_rq_fromuppy)
        call mpisend_nonblock_real(ubufyo_og(:,:,:,ivar1:ivar2),nbuf_y,yuneigh,touppy,isend_rq_touppy)
      endif
!
!  Set communication across z-planes.
!
      if (nprocz>1) then
!        
        lbufzo_og(:,:,:,ivar1:ivar2)=f_og(:,m1_ogrid:m2_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2) !lower z-planes
        nbuf_z=(/mx_ogrid,bufsizes_yz_og(INZL,IRCV),nghost,sendvar/)
        call mpirecv_nonblock_real(lbufzi_og(:,:,:,ivar1:ivar2),nbuf_z,zlneigh,touppz,irecv_rq_fromlowz)
        call mpisend_nonblock_real(lbufzo_og(:,:,:,ivar1:ivar2),nbuf_z,zlneigh,tolowz,isend_rq_tolowz)
!
        ubufzo_og(:,:,:,ivar1:ivar2)=f_og(:,m1_ogrid:m2_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2) !upper z-planes
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
        llbufo_og(:,:,:,ivar1:ivar2)=f_og(:,m1_ogrid:m1i_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2)
        if (llcorn>=0) then
          call mpirecv_nonblock_real(llbufi_og(:,:,:,ivar1:ivar2),nbuf_yz,llcorn,TOuu,irecv_rq_FRll)
          call mpisend_nonblock_real(llbufo_og(:,:,:,ivar1:ivar2),nbuf_yz,llcorn,TOll,isend_rq_TOll)
        endif
!
!  Upper y, lower z.
!
        ulbufo_og(:,:,:,ivar1:ivar2)=f_og(:,m2i_ogrid:m2_ogrid,n1_ogrid:n1i_ogrid,ivar1:ivar2)
        if (ulcorn>=0) then
          call mpirecv_nonblock_real(ulbufi_og(:,:,:,ivar1:ivar2),nbuf_yz,ulcorn,TOlu,irecv_rq_FRul)
          call mpisend_nonblock_real(ulbufo_og(:,:,:,ivar1:ivar2),nbuf_yz,ulcorn,TOul,isend_rq_TOul)
        endif
!
!  Upper y, upper z.
!
        uubufo_og(:,:,:,ivar1:ivar2)=f_og(:,m2i_ogrid:m2_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2)
        if (uucorn>=0) then
          call mpirecv_nonblock_real(uubufi_og(:,:,:,ivar1:ivar2),nbuf_yz,uucorn,TOll,irecv_rq_FRuu)
          call mpisend_nonblock_real(uubufo_og(:,:,:,ivar1:ivar2),nbuf_yz,uucorn,TOuu,isend_rq_TOuu)
        endif
!
!  Lower y, upper z.
!
        lubufo_og(:,:,:,ivar1:ivar2)=f_og(:,m1_ogrid:m1i_ogrid,n2i_ogrid:n2_ogrid,ivar1:ivar2)
        if (lucorn>=0) then
          call mpirecv_nonblock_real(lubufi_og(:,:,:,ivar1:ivar2),nbuf_yz,lucorn,TOul,irecv_rq_FRlu)
          call mpisend_nonblock_real(lubufo_og(:,:,:,ivar1:ivar2),nbuf_yz,lucorn,TOlu,isend_rq_TOlu)
        endif
!
      endif
!
    endsubroutine initiate_isendrcv_bdry_ogrid
!***********************************************************************
    subroutine finalize_isendrcv_bdry_ogrid(f_og)
!
!  Make sure the communications initiated with initiate_isendrcv_bdry are
!  finished and insert the just received boundary values.
!  Receive requests do not need to (and on OSF1 cannot) be explicitly
!  freed, since MPI_Wait takes care of this.
!
!  07-feb-17/Jorgen: Adapted from mpicomm.f90
!
      real, dimension(mx_ogrid,my_ogrid,mz_ogrid,mvar), intent(inout) :: f_og
      integer, parameter :: ivar1=1, ivar2=min(mcom,size(f_og,4))
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
          f_og(:,1:m1_ogrid-1,n1_ogrid:n2_ogrid,j)=lbufyi_og(:,:,:,j)       ! set lower buffer
          f_og(:,m2_ogrid+1:,n1_ogrid:n2_ogrid,j)=ubufyi_og(:,:,:,j)        ! set upper buffer
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
            f_og(:,m1_ogrid:m2_ogrid,1:n1_ogrid-1,j)=lbufzi_og(:,:,:,j)       ! set lower buffer
          enddo
        endif
        if (.not. llast_proc_z .or. bcz12(j,2)=='p') then 
          do j=ivar1,ivar2
            f_og(:,m1_ogrid:m2_ogrid,n2_ogrid+1:,j)=ubufzi_og(:,:,:,j)        ! set upper buffer
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
            f_og(:,1:m1_ogrid-1,1:n1_ogrid-1,j)=llbufi_og(:,:,:,j)               ! fill lower left corner (ll)
            f_og(:,m2_ogrid+1:,1:n1_ogrid-1,j)=ulbufi_og(:,:,:,j)                ! fill lower right corner (ul)
          endif
          if (.not. llast_proc_z .or. bcz12(j,2)=='p') then
            f_og(:,m2+1:,n2+1:,j)=uubufi_og(:,:,:,j)                             ! fill upper right corner (uu)
            f_og(:,1:m1-1,n2+1:,j)=lubufi_og(:,:,:,j)                            ! fill upper left corner (lu)
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
    subroutine isendrcv_bdry_x_ogrid(f_og)
!
!  Isend and Irecv boundary values for x-direction. Sends and receives
!  before continuing to y and z boundaries, as this allows the edges
!  of the grid to be set properly.
!
!  07-feb-17/Jorgen: Adapted from mpicomm.f90
!
      real, dimension (mx_ogrid,my_ogrid,mz_ogrid,mvar), intent(inout) :: f_og
      integer, parameter :: ivar1=1, ivar2=min(mcom,size(f_og,4))
      integer :: j
      integer, dimension(4), parameter ::  nbuf_x=(/nghost,ny_ogrid,nz_ogrid,ivar2-ivar1+1/)
!
      lbufxo_og(:,:,:,ivar1:ivar2)=f_og(l1_ogrid:l1i_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(lower x-zone)
      ubufxo_og(:,:,:,ivar1:ivar2)=f_og(l2i_ogrid:l2_ogrid,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2) !!(upper x-zone)
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
          f_og( 1:l1_ogrid-1,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j)=lbufxi_og(:,:,:,j)  !!(set lower buffer)
        enddo
      endif
      if (.not. llast_proc_x) then
        do j=ivar1,ivar2
          f_og(l2_ogrid+1:,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,j)=ubufxi_og(:,:,:,j)  !!(set upper buffer)
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
    subroutine initiate_isendrcv_bdry_filter(f_og,Hsize)
!
!  Isend and Irecv boundary filter values, called before filter is used.
!  Does not wait for the receives to finish (done in finalize_isendrcv_bdry)
!
!  29-nov-17/Jorgen: adapted from initiate_isendrcv_bdry_ogrid
!  22-jan-18/Jorgen: added filter zone in x-direction 
!
      real, dimension (mx_ogrid, my_ogrid, mz_ogrid,mvar) ::  f_og
      integer :: Hsize
      intent(in) :: f_og,Hsize
!
      integer, dimension(4), save :: nbuf_y, nbuf_x
      integer, parameter :: ivar1=1, ivar2=min(mcom,size(f_og,4))
      logical :: firstcall=.true.
!
      if (ivar2==0) return
      if(firstcall) then
        nbuf_y=(/nx_ogrid,Hsize,nz_ogrid,ivar2-ivar1+1/)
        nbuf_x=(/Hsize,ny_ogrid,nz_ogrid,ivar2-ivar1+1/)
        firstcall=.false.
      endif
!
!  Send/receive buffers across y-planes
!
      lbufyo_fi(:,:,:,:)= &
          f_og(l1_ogrid:l2_ogrid,m1i_ogrid+1:m1i_ogrid+Hsize,n1_ogrid:n2_ogrid,ivar1:ivar2)
      ubufyo_fi(:,:,:,:)= &
          f_og(l1_ogrid:l2_ogrid,m2i_ogrid-Hsize:m2i_ogrid-1,n1_ogrid:n2_ogrid,ivar1:ivar2)
      if (nprocy>1) then
!  Send/recieve lower y-zone
        call mpirecv_nonblock_real(lbufyi_fi,nbuf_y,ylneigh,touppy,irecv_rq_fromlowy_fi)
        call mpisend_nonblock_real(lbufyo_fi,nbuf_y,ylneigh,tolowy,isend_rq_tolowy_fi)
!  Send/recieve upper y-zone
        call mpirecv_nonblock_real(ubufyi_fi,nbuf_y,yuneigh,tolowy,irecv_rq_fromuppy_fi)
        call mpisend_nonblock_real(ubufyo_fi,nbuf_y,yuneigh,touppy,isend_rq_touppy_fi)
      endif
!
!  Send/receive buffers across x-planes
!  Since boundaries are non-periodic, this requires some more if-branching than for
!  y-direction
!
      if (nprocx>1) then
        if(ipx<nprocx-1) then
          ubufxo_fi(:,:,:,:)= &
             f_og(l2i_ogrid-Hsize:l2i_ogrid-1,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2)
!  Send/recieve upper y-zone
          call mpirecv_nonblock_real(ubufxi_fi,nbuf_x,xuneigh,tolowx,irecv_rq_fromuppx_fi)
          call mpisend_nonblock_real(ubufxo_fi,nbuf_x,xuneigh,touppx,isend_rq_touppx_fi)
        endif
        if(ipx>0) then
          lbufxo_fi(:,:,:,:)= &
             f_og(l1i_ogrid+1:l1i_ogrid+Hsize,m1_ogrid:m2_ogrid,n1_ogrid:n2_ogrid,ivar1:ivar2)
!  Send/recieve lower y-zone
          call mpirecv_nonblock_real(lbufxi_fi,nbuf_x,xlneigh,touppx,irecv_rq_fromlowx_fi)
          call mpisend_nonblock_real(lbufxo_fi,nbuf_x,xlneigh,tolowx,isend_rq_tolowx_fi)
        endif
      endif
    endsubroutine initiate_isendrcv_bdry_filter
!***********************************************************************
    subroutine finalize_isendrcv_bdry_filter(f_Hlox,f_Hupx,f_Hloy,f_Hupy,Hsize)
!
!  Make sure the communications initiated with initiate_isendrcv_bdry_filter
!  are finished and insert the just received boundary values.
!
!  29-nov-17/Jorgen: Adapted from finalize_isendrcv_bdry_ogrid
!
      integer, intent(in) :: Hsize
      real, dimension (Hsize,my_ogrid,nz_ogrid,mvar) ::  f_Hlox,f_Hupx
      real, dimension (mx_ogrid,Hsize,nz_ogrid,mvar) ::  f_Hloy,f_Hupy
      intent(inout) :: f_Hlox,f_Hupx,f_Hloy,f_Hupy

      integer, parameter :: ivar1=1, ivar2=min(mcom,size(f_Hloy,4))
      integer :: j
!
      if (ivar2==0) return
!
!  1. wait until data received
!  2. set filter ghost zones
!  3. wait until send completed, will be overwritten in next time step
!
!  Communication across y(theta)-planes, periodic BC 
!
      if (nprocy>1) then
        call mpiwait(irecv_rq_fromuppy_fi)
        call mpiwait(irecv_rq_fromlowy_fi)
        do j=ivar1,ivar2
          f_Hloy(l1_ogrid:l2_ogrid,:,:,j)=lbufyi_fi(1:nx_ogrid,:,:,j)   ! set lower filter halo
          f_Hupy(l1_ogrid:l2_ogrid,:,:,j)=ubufyi_fi(1:nx_ogrid,:,:,j)   ! set upper filter halo
        enddo
        call mpiwait(isend_rq_tolowy_fi)
        call mpiwait(isend_rq_touppy_fi)
      else
        do j=ivar1,ivar2
          f_Hloy(l1_ogrid:l2_ogrid,:,:,j)=ubufyo_fi(1:nx_ogrid,:,:,j)   ! set lower filter halo
          f_Hupy(l1_ogrid:l2_ogrid,:,:,j)=lbufyo_fi(1:nx_ogrid,:,:,j)   ! set upper filter halo
        enddo
      endif
      if (nprocx>1) then
        if(ipx<nprocx-1) then
          call mpiwait(irecv_rq_fromuppx_fi)
          do j=ivar1,ivar2
            f_Hupx(:,m1_ogrid:m2_ogrid,:,j)=ubufxi_fi(:,1:ny_ogrid,:,j)   ! set upper filter halo
          enddo
        endif
        if(ipx>0) then
          call mpiwait(irecv_rq_fromlowx_fi)
          do j=ivar1,ivar2
            f_Hlox(:,m1_ogrid:m2_ogrid,:,j)=lbufxi_fi(:,1:ny_ogrid,:,j)   ! set lower filter halo
          enddo
        endif
        if(ipx<nprocx-1) call mpiwait(isend_rq_touppx_fi)
        if(ipx>0) call mpiwait(isend_rq_tolowx_fi)
      endif
      call mpibarrier
!
    endsubroutine finalize_isendrcv_bdry_filter
!***********************************************************************
    subroutine cyclic_parallel_y(a,b,c,alpha,beta,r,x,n)
!
!  Inversion of a tridiagonal system with periodic BC (alpha and beta
!  coefficients in the left and right corners). Used in the ADI scheme of the
!  implicit_physics module.
!  Note: this subroutine is using twice the tridag one written above by tobi.
!
!  | b1 c1 0  ...       beta | | x1   |   | r1   |
!  | a2 b2 c2 ...            | | x2   |   | r2   |
!  | 0  a3 b3 c3             | | x3   | = | r3   |
!  |          ...            | | ...  |   | ...  |
!  |          an-1 bn-1 cn-1 | | xn-1 |   | rn-1 |
!  | alpha    0    a_n  b_n  | | xn   |   | rn   |
!
!  17.01.18/Jorgen: Adapded from general.f90
!
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: a,b,c,r
      real, dimension(n), intent(inout) :: x
      real, intent(in) :: alpha,beta

      real, dimension(n) :: bb,u,z
      integer :: i
      real    :: gamma,fact
      real, dimension(2) :: recvxNzN
!
      bb=b
      u=0.
      if(ipy==0) then
        gamma=-b(1)
        bb(1)=b(1)-gamma
        u(1)=gamma
        call mpisend_real(gamma,ylneigh,100)
      elseif(ipy==nprocy-1) then
        call mpirecv_real(gamma,yuneigh,100)
        bb(n)=b(n)-alpha*beta/gamma
        u(n)=alpha
      endif
      call tridag_parallel_y(a,bb,c,r,x,n)
      call tridag_parallel_y(a,bb,c,u,z,n)
      if(ipy==0) then
        call mpirecv_real(recvxNzN,2,ylneigh,119)
        fact=(x(1)+beta*recvxNzN(1)/gamma)/(1.+z(1)+beta*recvxNzN(2)/gamma)
      elseif(ipy==nprocy-1) then
        call mpisend_real((/x(n),z(n)/),2,yuneigh,119)
      endif
      call mpibcast_real(fact)
      do i=1,n
        x(i)=x(i)-fact*z(i)
      enddo
!
    endsubroutine cyclic_parallel_y
!***********************************************************************
    subroutine tridag_parallel_y(a,b,c,r,u,n)
!
!  Solves a tridiagonal system, where the system is distributed
!  over many processors.
!
!  17.01.18/Jorgen: Adapded from general.f90
!
!  | b1 c1 0  ...            | | u1   |   | r1   |
!  | a2 b2 c2 ...            | | u2   |   | r2   |
!  | 0  a3 b3 c3             | | u3   | = | r3   |
!  |          ...            | | ...  |   | ...  |
!  |          an-1 bn-1 cn-1 | | un-1 |   | rn-1 |
!  |          0    a_n  b_n  | | un   |   | rn   |
!
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: a,b,c,r
      real, dimension(n), intent(out) :: u

      real, dimension(n) :: gam
      integer :: j
      real :: bet
      real, dimension(3) :: recvBuf
!
      if(ipy==0) then
        bet=b(1)
        u(1)=r(1)/bet
      else
        call mpirecv_real(recvBuf,3,ylneigh,210)
        gam(1)=recvBuf(3)/recvBuf(1)
        bet=b(1)-a(1)*gam(1)
        u(1)=(r(1)-a(1)*recvBuf(2))/bet
      endif

      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      if(ipy<nprocy-1) then
        call mpisend_real((/bet,u(n),c(n)/),3,yuneigh,210)
        call mpirecv_real(recvBuf(1:2),2,yuneigh,211)
        u(n)=u(n)-recvBuf(1)*recvBuf(2)
      endif
!
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
      if(ipy>0) then
        call mpisend_real((/gam(1),u(1)/),2,ylneigh,211)
      endif
!
    endsubroutine tridag_parallel_y
!***********************************************************************
    subroutine tridag_parallel_x(a,b,c,r,u,n)
!
!  Solves a tridiagonal system, where the system is distributed
!  over many processors.
!
!  22.01.18/Jorgen: Adapded from tridag_parallel_y
!
!  | b1 c1 0  ...            | | u1   |   | r1   |
!  | a2 b2 c2 ...            | | u2   |   | r2   |
!  | 0  a3 b3 c3             | | u3   | = | r3   |
!  |          ...            | | ...  |   | ...  |
!  |          an-1 bn-1 cn-1 | | un-1 |   | rn-1 |
!  |          0    a_n  b_n  | | un   |   | rn   |
!
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: a,b,c,r
      real, dimension(n), intent(out) :: u

      real, dimension(n) :: gam
      integer :: j
      real :: bet
      real, dimension(3) :: recvBuf
!
      if(ipx==0) then
        bet=b(1)
        u(1)=r(1)/bet
      else
        call mpirecv_real(recvBuf,3,xlneigh,110)
        gam(1)=recvBuf(3)/recvBuf(1)
        bet=b(1)-a(1)*gam(1)
        u(1)=(r(1)-a(1)*recvBuf(2))/bet
      endif

      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      if(ipx<nprocx-1) then
        call mpisend_real((/bet,u(n),c(n)/),3,xuneigh,110)
        call mpirecv_real(recvBuf(1:2),2,xuneigh,111)
        u(n)=u(n)-recvBuf(1)*recvBuf(2)
      endif
!
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
      if(ipx>0) then
        call mpisend_real((/gam(1),u(1)/),2,xlneigh,111)
      endif
!
    endsubroutine tridag_parallel_x
!***********************************************************************
endmodule Solid_Cells_Mpicomm
