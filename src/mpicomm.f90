! $Id$

!!!!!!!!!!!!!!!!!!!!!
!!!  mpicomm.f90  !!!
!!!!!!!!!!!!!!!!!!!!!

!!! Module with MPI stuff

!  Data layout for each processor (`-' marks ghost points, `+' real
!  points of the processor shown)
!
!         n = mz        - - - - - - - - - . - - - - - - - - -
!             .         - - - - - - - - - . - - - - - - - - -
!             .         - - - - - - - - - . - - - - - - - - -
!             . n2      - - - + + + + + + . + + + + + + - - -
!             . .       - - - + + + + + + . + + + + + + - - -
!             . . n2i   - - - + + + + + + . + + + + + + - - -
!             . .       - - - + + + + + + . + + + + + + - - -
!             . .       - - - + + + + + + . + + + + + + - - -
!             . .       - - - + + + + + + . + + + + + + - - -
!                       . . . . . . . . . . . . . . . . . . .
!             . .       - - - + + + + + + . + + + + + + - - -
!             . .       - - - + + + + + + . + + + + + + - - -
!             . .       - - - + + + + + + . + + + + + + - - -
!             . . n1i   - - - + + + + + + . + + + + + + - - -
!             . .       - - - + + + + + + . + + + + + + - - -
!             . n1      - - - + + + + + + . + + + + + + - - -
!             3         - - - - - - - - - . - - - - - - - - -
!             2         - - - - - - - - - . - - - - - - - - -
!         n = 1         - - - - - - - - - . - - - - - - - - -
!
!                                m1i             m2i
!                             m1. . . . .   . . . . . m2
!               m     = 1 2 3 . . . . . .   . . . . . . . . my
!
!  Thus, the indices for some important regions are
!    ghost zones:
!                        1:nghost (1:m1-1)  and  my-nghost+1:my (m2+1:my)
!    real points:
!                        m1:m2
!    boundary points (which become ghost points for adjacent processors):
!                        m1:m1i  and  m2i:m2
!    inner points for periodic bc (i.e. points where 7-point derivatives are
!    unaffected by ghost information):
!                        m1i+1:m2i-1
!    inner points for general bc (i.e. points where 7-point derivatives are
!    unaffected by ghost information plus boundcond for m1,m2):
!                        m1i+2:m2i-2

module Mpicomm

  use Cparam
  use Cdata, only: iproc,ipx,ipy,ipz,root,lroot, &
                   xlneigh,ylneigh,zlneigh,xuneigh,yuneigh,zuneigh, &
                   llcorn,lucorn,uucorn,ulcorn, &
                   lemergency_brake

  implicit none

  include 'mpicomm.h'

  interface mpirecv_logical
     module procedure mpirecv_logical_scl
     module procedure mpirecv_logical_arr
  endinterface

  interface mpirecv_real
    module procedure mpirecv_real_scl
    module procedure mpirecv_real_arr
    module procedure mpirecv_real_arr2
    module procedure mpirecv_real_arr3
    module procedure mpirecv_real_arr4
  endinterface

  interface mpirecv_int
    module procedure mpirecv_int_scl
    module procedure mpirecv_int_arr
    module procedure mpirecv_int_arr2
  endinterface

  interface mpisend_logical
     module procedure mpisend_logical_scl
     module procedure mpisend_logical_arr
  endinterface

  interface mpisend_real
    module procedure mpisend_real_scl
    module procedure mpisend_real_arr
    module procedure mpisend_real_arr2
    module procedure mpisend_real_arr3
    module procedure mpisend_real_arr4
  endinterface

  interface mpisend_int
    module procedure mpisend_int_scl
    module procedure mpisend_int_arr
    module procedure mpisend_int_arr2
 endinterface

  interface mpibcast_logical
    module procedure mpibcast_logical_scl
    module procedure mpibcast_logical_arr
    module procedure mpibcast_logical_arr2
  endinterface

  interface mpibcast_int
    module procedure mpibcast_int_scl
    module procedure mpibcast_int_arr
  endinterface

  interface mpibcast_real
    module procedure mpibcast_real_scl
    module procedure mpibcast_real_arr
    module procedure mpibcast_real_arr2
    module procedure mpibcast_real_arr3
  endinterface

  interface mpibcast_double
    module procedure mpibcast_double_scl
    module procedure mpibcast_double_arr
  endinterface

  interface mpibcast_char
    module procedure mpibcast_char_scl
    module procedure mpibcast_char_arr
  endinterface

  interface mpiallreduce_sum
    module procedure mpiallreduce_sum_arr2
    module procedure mpiallreduce_sum_arr
    module procedure mpiallreduce_sum_scl
  endinterface

  interface mpiallreduce_sum_int
     module procedure mpiallreduce_sum_int_arr
     module procedure mpiallreduce_sum_int_scl
  endinterface

  interface mpiallreduce_max
    module procedure mpiallreduce_max_arr
    module procedure mpiallreduce_max_scl
  endinterface

  interface mpireduce_max
    module procedure mpireduce_max_arr
    module procedure mpireduce_max_scl
  endinterface

  interface mpireduce_min
    module procedure mpireduce_min_arr
    module procedure mpireduce_min_scl
  endinterface

  interface mpireduce_or
    module procedure mpireduce_or_arr
    module procedure mpireduce_or_scl
  endinterface

  interface mpireduce_and
    module procedure mpireduce_and_arr
    module procedure mpireduce_and_scl
  endinterface

! NOT POSSIBLE BECAUSE OF n-Dimensional array usage
! in equ.f90
!  interface mpireduce_sum
!    module procedure mpireduce_sum_arr
!    module procedure mpireduce_sum_scl
!  endinterface

  interface mpireduce_sum_double
    module procedure mpireduce_sum_double_arr
    module procedure mpireduce_sum_double_scl
  endinterface

  interface mpireduce_sum_int
    module procedure mpireduce_sum_int_arr
    module procedure mpireduce_sum_int_scl
  endinterface

  include 'mpif.h'
!
!  initialize debug parameter for this routine
!
  logical :: ldebug_mpi=.false.
!
!  For f-array processor boundaries
!
  real, dimension (nghost,ny,nz,mcom) :: lbufxi,ubufxi,lbufxo,ubufxo
  real, dimension (mx,nghost,nz,mcom) :: lbufyi,ubufyi,lbufyo,ubufyo
  real, dimension (mx,ny,nghost,mcom) :: lbufzi,ubufzi,lbufzo,ubufzo
  real, dimension (mx,nghost,nghost,mcom) :: llbufi,lubufi,uubufi,ulbufi
  real, dimension (mx,nghost,nghost,mcom) :: llbufo,lubufo,uubufo,ulbufo
!
  real, dimension (nghost,my,mz,mcom) :: fahi,falo,fbhi,fblo,fao,fbo ! For shear
  integer :: nextya, nextyb, lastya, lastyb, displs ! For shear
  integer :: ierr
  integer :: nprocs
!
!  mpi tags
!
  integer :: tolowx=13,touppx=14,tolowy=3,touppy=4,tolowz=5,touppz=6 ! msg. tags
  integer :: TOll=7,TOul=8,TOuu=9,TOlu=10 ! msg. tags for corners
  integer :: io_perm=20,io_succ=21
!
!  mpi tags for radiation
!  the values for those have to differ by a number greater than maxdir=190
!  in order to have unique tags for each boundary and each direction
!
  integer, parameter :: Qtag_zx=300,Qtag_xy=350
  integer, parameter :: tautag_zx=400,tautag_xy=450
  integer, parameter :: Qtag_peri_zx=1000,Qtag_peri_xy=2000
  integer, parameter :: tautag_peri_zx=3000,tautag_peri_xy=4000
!
!  Communicators
!
  integer :: MPI_COMM_ROW,MPI_COMM_COLUMN
!
  integer :: isend_rq_tolowx,isend_rq_touppx,irecv_rq_fromlowx,irecv_rq_fromuppx
  integer :: isend_rq_tolowy,isend_rq_touppy,irecv_rq_fromlowy,irecv_rq_fromuppy
  integer :: isend_rq_tolowz,isend_rq_touppz,irecv_rq_fromlowz,irecv_rq_fromuppz
  integer :: isend_rq_TOll,isend_rq_TOul,isend_rq_TOuu,isend_rq_TOlu  !(corners)
  integer :: irecv_rq_FRuu,irecv_rq_FRlu,irecv_rq_FRll,irecv_rq_FRul  !(corners)
  integer :: isend_rq_tolastya,isend_rq_tonextya, &
             irecv_rq_fromlastya,irecv_rq_fromnextya ! For shear
  integer :: isend_rq_tolastyb,isend_rq_tonextyb, &
             irecv_rq_fromlastyb,irecv_rq_fromnextyb ! For shear

  integer, dimension (MPI_STATUS_SIZE) :: isend_stat_tl,isend_stat_tu
  integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fl,irecv_stat_fu
  integer, dimension (MPI_STATUS_SIZE) :: isend_stat_Tll,isend_stat_Tul, &
                                          isend_stat_Tuu,isend_stat_Tlu
  integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_Fuu,irecv_stat_Flu, &
                                          irecv_stat_Fll,irecv_stat_Ful

  contains

!***********************************************************************
    subroutine mpicomm_init()
!
!  Initialise MPI communication and set up some variables
!  The arrays leftneigh and rghtneigh give the processor numbers
!  to the left and to the right.
!
!
!  20-aug-01/wolf: coded
!  31-aug-01/axel: added to 3-D
!  15-sep-01/axel: adapted from Wolfgang's version
!  21-may-02/axel: communication of corners added
!   6-jun-02/axel: generalized to allow for ny=1
!  23-nov-02/axel: corrected problem with ny=4 or less
!
      use Cdata, only: lmpicomm,lprocz_slowest
!
!  get processor number, number of procs, and whether we are root
!
      lmpicomm = .true.
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, iproc , ierr)
      lroot = (iproc==root)
!
!  consistency checks
!
      if (nprocx /= 1) then
        if (lroot) print*,'WARNING: nprocx > 1 is not yet well tested'
      endif
!
!  check total number of processors
!
      if (nprocs /= nprocx*nprocy*nprocz) then
        if (lroot) then
          print*, 'Compiled with NCPUS = ', ncpus, &
          ', but running on ', nprocs, ' processors'
        endif
        call stop_it('Inconsistency 1')
      endif
!
!  preliminary additions for nprocx \= 1
!  write the value of nprocx*nx at the end for compatibility reasons
!
      if ((nprocx*nx /= nxgrid) .or. &
          (nprocy*ny /= nygrid) .or. &
          (nprocz*nz /= nzgrid)) then
        if (lroot) then
          write(6,'(A,3I4,A,3I4,A)') &
               'nproc[x-z]*n[x-z] = (', &
               nprocy*ny, nprocz*nz, nprocx*nx, &
               ') /= n[xyz]grid= (', &
               nygrid, nzgrid, nxgrid, ")"
        endif
        call stop_it('Inconsistency 2')
      endif
!
!  avoid overlapping ghost zones
!
      if ((nx<nghost) .and. (nxgrid/=1)) &
           call stop_it('Overlapping ghost zones in x-direction: reduce nprocx')
      if ((ny<nghost) .and. (nygrid/=1)) &
           call stop_it('Overlapping ghost zones in y-direction: reduce nprocy')
      if ((nz<nghost) .and. (nzgrid/=1)) &
           call stop_it('Overlapping ghost zones in z-direction: reduce nprocz')
!
!  position on the processor grid
!  x is fastest direction, z slowest (this is the default)
!
      if (lprocz_slowest) then
        ipx = modulo(iproc, nprocx)
        ipy = modulo(iproc/nprocx, nprocy)
        ipz = iproc/(nprocx*nprocy)
       else
        ipx = modulo(iproc, nprocx)
        ipy = iproc/(nprocx*nprocy)
        ipz = modulo(iproc/nprocx, nprocy)
       endif
!
!  set up `lower' and `upper' neighbours
!
      xlneigh = (ipz*nprocx*nprocy+ipy*nprocx+modulo(ipx-1,nprocx))
      xuneigh = (ipz*nprocx*nprocy+ipy*nprocx+modulo(ipx+1,nprocx))
      ylneigh = (ipz*nprocx*nprocy+modulo(ipy-1,nprocy)*nprocx+ipx)
      yuneigh = (ipz*nprocx*nprocy+modulo(ipy+1,nprocy)*nprocx+ipx)
      zlneigh = (modulo(ipz-1,nprocz)*nprocx*nprocy+ipy*nprocx+ipx)
      zuneigh = (modulo(ipz+1,nprocz)*nprocx*nprocy+ipy*nprocx+ipx)
!
!  set the four corners in the yz-plane (in cyclic order)
!
      llcorn=ipx+(modulo(ipy-1,nprocy)+modulo(ipz-1,nprocz)*nprocy)*nprocx
      ulcorn=ipx+(modulo(ipy+1,nprocy)+modulo(ipz-1,nprocz)*nprocy)*nprocx
      uucorn=ipx+(modulo(ipy+1,nprocy)+modulo(ipz+1,nprocz)*nprocy)*nprocx
      lucorn=ipx+(modulo(ipy-1,nprocy)+modulo(ipz+1,nprocz)*nprocy)*nprocx
!
!  this value is not yet the one read in, but the one initialized in cparam.f90
!
!  Print neighbors in counterclockwise order (including the corners),
!  starting with left neighbor.
!  Example with 4x4 processors
!   3 |  0   1   2   3 |  0
!  ---+----------------+---
!  15 | 12  13  14  15 | 12
!  11 |  8   9  10  11 |  8
!   7 |  4   5   6   7 |  4
!   3 |  0   1   2   3 |  0
!  ---+----------------+---
!  15 | 12  13  14  15 | 12
!  should print (3,15,12,13,1,5,4,7) for iproc=0
!
!  print processor numbers and those of their neighbors
!  NOTE: the ip print parameter has *not* yet been read at this point.
!  Therefore it must be invoked by resetting ldebug_mpi appropriately.
!
      if (ldebug_mpi) write(6,'(A,I4,"(",3I4,"): ",8I4)') &
        'mpicomm_init: MPICOMM neighbors ', &
        iproc,ipx,ipy,ipz, &
        ylneigh,llcorn,zlneigh,ulcorn,yuneigh,uucorn,zuneigh,lucorn
!
!  Define MPI communicator MPI_COMM_ROW that includes all processes
!  sharing the same value of ipz. The rank within MPI_COMM_WORLD is
!  given by ipy.
!
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, ipz, ipy, MPI_COMM_ROW, ierr)
!
!  Do the same for columns - all processes sharing the same value of 
!  ipy. The rank within MPI_COMM_WORLD is given by ipz
!
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, ipy, ipz, MPI_COMM_COLUMN, ierr)
!
    endsubroutine mpicomm_init
!***********************************************************************
    subroutine initiate_isendrcv_bdry(f,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values. Called in the beginning of pde.
!  Does not wait for the receives to finish (done in finalize_isendrcv_bdry)
!  leftneigh and rghtneigh are initialized by mpicomm_init
!
!  21-may-02/axel: communication of corners added
!  11-aug-07/axel: communication in the x-direction added
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbufx, nbufy, nbufz, nbufyz
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  Periodic boundary conditions in x
!
      if (nprocx>1) then
        lbufxo(:,:,:,ivar1:ivar2)=f(l1:l1i,m1:m2,n1:n2,ivar1:ivar2) !!(lower x-zone)
        ubufxo(:,:,:,ivar1:ivar2)=f(l2i:l2,m1:m2,n1:n2,ivar1:ivar2) !!(upper x-zone)
        nbufx=ny*nz*nghost*(ivar2-ivar1+1)
        call MPI_IRECV(ubufxi(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xuneigh,tolowx,MPI_COMM_WORLD,irecv_rq_fromuppx,ierr)
        call MPI_IRECV(lbufxi(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xlneigh,touppx,MPI_COMM_WORLD,irecv_rq_fromlowx,ierr)
        call MPI_ISEND(lbufxo(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xlneigh,tolowx,MPI_COMM_WORLD,isend_rq_tolowx,ierr)
        call MPI_ISEND(ubufxo(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xuneigh,touppx,MPI_COMM_WORLD,isend_rq_touppx,ierr)
      endif
!
!  Periodic boundary conditions in y
!
      if (nprocy>1) then
        lbufyo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n1:n2,ivar1:ivar2) !!(lower y-zone)
        ubufyo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n1:n2,ivar1:ivar2) !!(upper y-zone)
        nbufy=mx*nz*nghost*(ivar2-ivar1+1)
        call MPI_IRECV(ubufyi(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
            yuneigh,tolowy,MPI_COMM_WORLD,irecv_rq_fromuppy,ierr)
        call MPI_IRECV(lbufyi(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
            ylneigh,touppy,MPI_COMM_WORLD,irecv_rq_fromlowy,ierr)
        call MPI_ISEND(lbufyo(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
            ylneigh,tolowy,MPI_COMM_WORLD,isend_rq_tolowy,ierr)
        call MPI_ISEND(ubufyo(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
            yuneigh,touppy,MPI_COMM_WORLD,isend_rq_touppy,ierr)
      endif
!
!  Periodic boundary conditions in z
!
      if (nprocz>1) then
        lbufzo(:,:,:,ivar1:ivar2)=f(:,m1:m2,n1:n1i,ivar1:ivar2) !!(lower z-zone)
        ubufzo(:,:,:,ivar1:ivar2)=f(:,m1:m2,n2i:n2,ivar1:ivar2) !!(upper z-zone)
        nbufz=mx*ny*nghost*(ivar2-ivar1+1)
        call MPI_IRECV(ubufzi(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
            zuneigh,tolowz,MPI_COMM_WORLD,irecv_rq_fromuppz,ierr)
        call MPI_IRECV(lbufzi(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
            zlneigh,touppz,MPI_COMM_WORLD,irecv_rq_fromlowz,ierr)
        call MPI_ISEND(lbufzo(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
            zlneigh,tolowz,MPI_COMM_WORLD,isend_rq_tolowz,ierr)
        call MPI_ISEND(ubufzo(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
            zuneigh,touppz,MPI_COMM_WORLD,isend_rq_touppz,ierr)
      endif
!
!  The four corners (in counter-clockwise order)
!  (NOTE: the case nprocx>1 is not yet considered here!)
!
      if (nprocy>1.and.nprocz>1) then
        llbufo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n1:n1i,ivar1:ivar2)
        ulbufo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n1:n1i,ivar1:ivar2)
        uubufo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n2i:n2,ivar1:ivar2)
        lubufo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n2i:n2,ivar1:ivar2)
        nbufyz=mx*nghost*nghost*(ivar2-ivar1+1)
        call MPI_IRECV(uubufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            uucorn,TOll,MPI_COMM_WORLD,irecv_rq_FRuu,ierr)
        call MPI_IRECV(lubufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            lucorn,TOul,MPI_COMM_WORLD,irecv_rq_FRlu,ierr)
        call MPI_IRECV(llbufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            llcorn,TOuu,MPI_COMM_WORLD,irecv_rq_FRll,ierr)
        call MPI_IRECV(ulbufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            ulcorn,TOlu,MPI_COMM_WORLD,irecv_rq_FRul,ierr)
        call MPI_ISEND(llbufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            llcorn,TOll,MPI_COMM_WORLD,isend_rq_TOll,ierr)
        call MPI_ISEND(ulbufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            ulcorn,TOul,MPI_COMM_WORLD,isend_rq_TOul,ierr)
        call MPI_ISEND(uubufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            uucorn,TOuu,MPI_COMM_WORLD,isend_rq_TOuu,ierr)
        call MPI_ISEND(lubufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            lucorn,TOlu,MPI_COMM_WORLD,isend_rq_TOlu,ierr)
      endif
!
!  communication sample
!  (commented out, because compiler does like this for 0-D runs)
!
!     if (ip<7.and.ipy==0.and.ipz==3) &
!       print*,'initiate_isendrcv_bdry: MPICOMM send lu: ',iproc,lubufo(nx/2+4,:,1,2),' to ',lucorn
!
    endsubroutine initiate_isendrcv_bdry
!***********************************************************************
    subroutine finalize_isendrcv_bdry(f,ivar1_opt,ivar2_opt)
!
      use Cdata, only: bcx1,bcx2,bcy1,bcy2,bcz1,bcz2
!
!  Make sure the communications initiated with initiate_isendrcv_bdry are
!  finished and insert the just received boundary values.
!   Receive requests do not need to (and on OSF1 cannot) be explicitly
!  freed, since MPI_Wait takes care of this.
!
!  21-may-02/axel: communication of corners added
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  1. wait until data received
!  2. set ghost zones
!  3. wait until send completed, will be overwritten in next time step
!
!  Communication in x (includes periodic bc)
!
      if (nprocx>1) then
        call MPI_WAIT(irecv_rq_fromuppx,irecv_stat_fu,ierr)
        call MPI_WAIT(irecv_rq_fromlowx,irecv_stat_fl,ierr)
        do j=ivar1,ivar2
          if (ipx/=0 .or. bcx1(j)=='p') then
            f( 1:l1-1,m1:m2,n1:n2,j)=lbufxi(:,:,:,j)  !!(set lower buffer)
          endif
          if (ipx/=nprocx-1 .or. bcx2(j)=='p') then
            f(l2+1:mx,m1:m2,n1:n2,j)=ubufxi(:,:,:,j)  !!(set upper buffer)
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowx,isend_stat_tl,ierr)
        call MPI_WAIT(isend_rq_touppx,isend_stat_tu,ierr)
      endif
!
!  Communication in y (includes periodic bc)
!
      if (nprocy>1) then
        call MPI_WAIT(irecv_rq_fromuppy,irecv_stat_fu,ierr)
        call MPI_WAIT(irecv_rq_fromlowy,irecv_stat_fl,ierr)
        do j=ivar1,ivar2
          if (ipy/=0 .or. bcy1(j)=='p') then
            f(:, 1:m1-1,n1:n2,j)=lbufyi(:,:,:,j)  !!(set lower buffer)
          endif
          if (ipy/=nprocy-1 .or. bcy2(j)=='p') then
            f(:,m2+1:my,n1:n2,j)=ubufyi(:,:,:,j)  !!(set upper buffer)
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowy,isend_stat_tl,ierr)
        call MPI_WAIT(isend_rq_touppy,isend_stat_tu,ierr)
      endif
!
!  Communication in z (includes periodic bc)
!
      if (nprocz>1) then
        call MPI_WAIT(irecv_rq_fromuppz,irecv_stat_fu,ierr)
        call MPI_WAIT(irecv_rq_fromlowz,irecv_stat_fl,ierr)
        do j=ivar1,ivar2
          if (ipz/=0 .or. bcz1(j)=='p') then
            f(:,m1:m2, 1:n1-1,j)=lbufzi(:,:,:,j)  !!(set lower buffer)
          endif
          if (ipz/=nprocz-1 .or. bcz2(j)=='p') then
            f(:,m1:m2,n2+1:mz,j)=ubufzi(:,:,:,j)  !!(set upper buffer)
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowz,isend_stat_tl,ierr)
        call MPI_WAIT(isend_rq_touppz,isend_stat_tu,ierr)
      endif
!
!  The four yz-corners (in counter-clockwise order)
!
      if (nprocy>1.and.nprocz>1) then
        call MPI_WAIT(irecv_rq_FRuu,irecv_stat_Fuu,ierr)
        call MPI_WAIT(irecv_rq_FRlu,irecv_stat_Flu,ierr)
        call MPI_WAIT(irecv_rq_FRll,irecv_stat_Fll,ierr)
        call MPI_WAIT(irecv_rq_FRul,irecv_stat_Ful,ierr)
        do j=ivar1,ivar2
          if (ipz/=0 .or. bcz1(j)=='p') then
            if (ipy/=0 .or. bcy1(j)=='p') then
              f(:, 1:m1-1, 1:n1-1,j)=llbufi(:,:,:,j)  !!(set ll corner)
            endif
            if (ipy/=nprocy-1 .or. bcy2(j)=='p') then
              f(:,m2+1:my, 1:n1-1,j)=ulbufi(:,:,:,j)  !!(set ul corner)
            endif
          endif
          if (ipz/=nprocz-1 .or. bcz2(j)=='p') then
            if (ipy/=nprocy-1 .or. bcy2(j)=='p') then
              f(:,m2+1:my,n2+1:mz,j)=uubufi(:,:,:,j)  !!(set uu corner)
            endif
            if (ipy/=0 .or. bcy1(j)=='p') then
              f(:, 1:m1-1,n2+1:mz,j)=lubufi(:,:,:,j)  !!(set lu corner)
            endif
          endif
        enddo
        call MPI_WAIT(isend_rq_TOll,isend_stat_Tll,ierr)
        call MPI_WAIT(isend_rq_TOul,isend_stat_Tul,ierr)
        call MPI_WAIT(isend_rq_TOuu,isend_stat_Tuu,ierr)
        call MPI_WAIT(isend_rq_TOlu,isend_stat_Tlu,ierr)
      endif
!
!  communication sample
!  (commented out, because compiler does like this for 0-D runs)
!
!     if (ip<7.and.ipy==3.and.ipz==0) &
!       print*,'finalize_isendrcv_bdry: MPICOMM recv ul: ', &
!                       iproc,ulbufi(nx/2+4,:,1,2),' from ',ulcorn
!
!  make sure the other precessors don't carry on sending new data
!  which could be mistaken for an earlier time
!
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!
    endsubroutine finalize_isendrcv_bdry
!***********************************************************************
    subroutine initiate_isendrcv_shockbdry(f,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values. Called in the beginning of pde.
!  Does not wait for the receives to finish (done in finalize_isendrcv_bdry)
!  leftneigh and rghtneigh are initialized by mpicomm_init
!
!  21-may-02/axel: communication of corners added
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbufx, nbufy, nbufz, nbufyz
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  Periodic boundary conditions in x
!
      if (nprocx>1) then
        lbufxo(:,:,:,1:ivar2-ivar1+1)=f( 1:l1-1,m1:m2,n1:n2,ivar1:ivar2) !!(lower x-zone)
        ubufxo(:,:,:,1:ivar2-ivar1+1)=f(l2+1:mx,m1:m2,n1:n2,ivar1:ivar2) !!(upper x-zone)
        nbufx=ny*nz*nghost*(ivar2-ivar1+1)
        call MPI_IRECV(ubufxi,nbufx,MPI_REAL,xuneigh,tolowx,MPI_COMM_WORLD, &
            irecv_rq_fromuppx,ierr)
        call MPI_IRECV(lbufxi,nbufx,MPI_REAL,xlneigh,touppx,MPI_COMM_WORLD, &
            irecv_rq_fromlowx,ierr)
        call MPI_ISEND(lbufxo,nbufx,MPI_REAL,xlneigh,tolowx,MPI_COMM_WORLD, &
            isend_rq_tolowx,ierr)
        call MPI_ISEND(ubufxo,nbufx,MPI_REAL,xuneigh,touppx,MPI_COMM_WORLD, &
            isend_rq_touppx,ierr)
      endif
!
!  Periodic boundary conditions in y
!
      if (nprocy>1) then
        lbufyo(:,:,:,1:ivar2-ivar1+1)=f(:,1:m1-1,n1:n2,ivar1:ivar2) !!(lower y-zone)
        ubufyo(:,:,:,1:ivar2-ivar1+1)=f(:,m2+1:my,n1:n2,ivar1:ivar2) !!(upper y-zone)
        nbufy=mx*nz*nghost*(ivar2-ivar1+1)
        call MPI_IRECV(ubufyi,nbufy,MPI_REAL,yuneigh,tolowy,MPI_COMM_WORLD, &
            irecv_rq_fromuppy,ierr)
        call MPI_IRECV(lbufyi,nbufy,MPI_REAL,ylneigh,touppy,MPI_COMM_WORLD, &
            irecv_rq_fromlowy,ierr)
        call MPI_ISEND(lbufyo,nbufy,MPI_REAL,ylneigh,tolowy,MPI_COMM_WORLD, &
            isend_rq_tolowy,ierr)
        call MPI_ISEND(ubufyo,nbufy,MPI_REAL,yuneigh,touppy,MPI_COMM_WORLD, &
            isend_rq_touppy,ierr)
      endif
!
!  Periodic boundary conditions in z
!
      if (nprocz>1) then
        lbufzo(:,:,:,1:ivar2-ivar1+1)=f(:,m1:m2,1:n1-1,ivar1:ivar2) !!(lower z-zone)
        ubufzo(:,:,:,1:ivar2-ivar1+1)=f(:,m1:m2,n2+1:mz,ivar1:ivar2) !!(upper z-zone)
        nbufz=mx*ny*nghost*(ivar2-ivar1+1)
        call MPI_IRECV(ubufzi,nbufz,MPI_REAL,zuneigh,tolowz,MPI_COMM_WORLD, &
            irecv_rq_fromuppz,ierr)
        call MPI_IRECV(lbufzi,nbufz,MPI_REAL,zlneigh,touppz,MPI_COMM_WORLD, &
            irecv_rq_fromlowz,ierr)
        call MPI_ISEND(lbufzo,nbufz,MPI_REAL,zlneigh,tolowz,MPI_COMM_WORLD, &
            isend_rq_tolowz,ierr)
        call MPI_ISEND(ubufzo,nbufz,MPI_REAL,zuneigh,touppz,MPI_COMM_WORLD, &
            isend_rq_touppz,ierr)
      endif
!
!  The four corners (in counter-clockwise order)
!
      if (nprocy>1.and.nprocz>1) then
        llbufo(:,:,:,1:ivar2-ivar1+1)=f(:,1:m1-1,1:n1-1,ivar1:ivar2)
        ulbufo(:,:,:,1:ivar2-ivar1+1)=f(:,m2+1:my,1:n1-1,ivar1:ivar2)
        uubufo(:,:,:,1:ivar2-ivar1+1)=f(:,m2+1:my,n2+1:mz,ivar1:ivar2)
        lubufo(:,:,:,1:ivar2-ivar1+1)=f(:,1:m1-1,n2+1:mz,ivar1:ivar2)
        nbufyz=mx*nghost*nghost*(ivar2-ivar1+1)
        call MPI_IRECV(uubufi,nbufyz,MPI_REAL,uucorn,TOll,MPI_COMM_WORLD, &
            irecv_rq_FRuu,ierr)
        call MPI_IRECV(lubufi,nbufyz,MPI_REAL,lucorn,TOul,MPI_COMM_WORLD, &
            irecv_rq_FRlu,ierr)
        call MPI_IRECV(llbufi,nbufyz,MPI_REAL,llcorn,TOuu,MPI_COMM_WORLD, &
            irecv_rq_FRll,ierr)
        call MPI_IRECV(ulbufi,nbufyz,MPI_REAL,ulcorn,TOlu,MPI_COMM_WORLD, &
            irecv_rq_FRul,ierr)
        call MPI_ISEND(llbufo,nbufyz,MPI_REAL,llcorn,TOll,MPI_COMM_WORLD, &
            isend_rq_TOll,ierr)
        call MPI_ISEND(ulbufo,nbufyz,MPI_REAL,ulcorn,TOul,MPI_COMM_WORLD, &
            isend_rq_TOul,ierr)
        call MPI_ISEND(uubufo,nbufyz,MPI_REAL,uucorn,TOuu,MPI_COMM_WORLD, &
            isend_rq_TOuu,ierr)
        call MPI_ISEND(lubufo,nbufyz,MPI_REAL,lucorn,TOlu,MPI_COMM_WORLD, &
            isend_rq_TOlu,ierr)
      endif
!
!  communication sample
!  (commented out, because compiler does like this for 0-D runs)
!
!     if (ip<7.and.ipy==0.and.ipz==3) &
!       print*,'initiate_isendrcv_bdry: MPICOMM send lu: ',iproc,lubufo(nx/2+4,:,1,2),' to ',lucorn
!
    endsubroutine initiate_isendrcv_shockbdry
!***********************************************************************
    subroutine finalize_isendrcv_shockbdry(f,ivar1_opt,ivar2_opt)
!
      use Cdata, only: bcx1,bcx2,bcy1,bcy2,bcz1,bcz2
!
!  Make sure the communications initiated with initiate_isendrcv_bdry are
!  finished and insert the just received boundary values.
!   Receive requests do not need to (and on OSF1 cannot) be explicitly
!  freed, since MPI_Wait takes care of this.
!
!  21-may-02/axel: communication of corners added
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  1. wait until data received
!  2. set ghost zones
!  3. wait until send completed, will be overwritten in next time step
!
!  Communication in x (includes periodic bc)
!
      if (nprocx>1) then
        call MPI_WAIT(irecv_rq_fromuppx,irecv_stat_fu,ierr)
        call MPI_WAIT(irecv_rq_fromlowx,irecv_stat_fl,ierr)
        do j=ivar1,ivar2
          if (ipx/=0 .or. bcx1(j)=='p') then
            f(l1:l1i,m1:m2,n1:n2,j)=f(l1:l1i,m1:m2,n1:n2,j)+lbufxi(:,:,:,j-ivar1+1)  !!(set lower buffer)
          endif
          if (ipx/=nprocx-1 .or. bcx2(j)=='p') then
            f(l2i:l2,m1:m2,n1:n2,j)=f(l2i:l2,m1:m2,n1:n2,j)+ubufxi(:,:,:,j-ivar1+1)  !!(set upper buffer)
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowx,isend_stat_tl,ierr)
        call MPI_WAIT(isend_rq_touppx,isend_stat_tu,ierr)
      endif
!
!  Communication in y (includes periodic bc)
!
      if (nprocy>1) then
        call MPI_WAIT(irecv_rq_fromuppy,irecv_stat_fu,ierr)
        call MPI_WAIT(irecv_rq_fromlowy,irecv_stat_fl,ierr)
        do j=ivar1,ivar2
          if (ipy/=0 .or. bcy1(j)=='p') then
            f(:, m1:m1i,n1:n2,j)=f(:, m1:m1i,n1:n2,j) + &
                          lbufyi(:,:,:,j-ivar1+1)  !!(set lower buffer)
          endif
          if (ipy/=nprocy-1 .or. bcy2(j)=='p') then
            f(:,m2i:m2,n1:n2,j)=f(:,m2i:m2,n1:n2,j) + &
                          ubufyi(:,:,:,j-ivar1+1)  !!(set upper buffer)
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowy,isend_stat_tl,ierr)
        call MPI_WAIT(isend_rq_touppy,isend_stat_tu,ierr)
      endif
!
!  Communication in z (includes periodic bc)
!
      if (nprocz>1) then
        call MPI_WAIT(irecv_rq_fromuppz,irecv_stat_fu,ierr)
        call MPI_WAIT(irecv_rq_fromlowz,irecv_stat_fl,ierr)
        do j=ivar1,ivar2
          if (ipz/=0 .or. bcz1(j)=='p') then
            f(:,m1:m2, n1:n1i,j)=f(:,m1:m2, n1:n1i,j) + &
                                 lbufzi(:,:,:,j-ivar1+1)  !!(set lower buffer)
          endif
          if (ipz/=nprocz-1 .or. bcz2(j)=='p') then
            f(:,m1:m2,n2i:n2,j)=f(:,m1:m2,n2i:n2,j) + &
                                 ubufzi(:,:,:,j-ivar1+1)  !!(set upper buffer)
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowz,isend_stat_tl,ierr)
        call MPI_WAIT(isend_rq_touppz,isend_stat_tu,ierr)
      endif
!
!  The four yz-corners (in counter-clockwise order)
!
      if (nprocy>1.and.nprocz>1) then
        call MPI_WAIT(irecv_rq_FRuu,irecv_stat_Fuu,ierr)
        call MPI_WAIT(irecv_rq_FRlu,irecv_stat_Flu,ierr)
        call MPI_WAIT(irecv_rq_FRll,irecv_stat_Fll,ierr)
        call MPI_WAIT(irecv_rq_FRul,irecv_stat_Ful,ierr)
        do j=ivar1,ivar2
          if (ipz/=0 .or. bcz1(j)=='p') then
            if (ipy/=0 .or. bcy1(j)=='p') then
              f(:, m1:m1i, n1:n1i,j)=f(:, m1:m1i, n1:n1i,j) + &
                                   llbufi(:,:,:,j-ivar1+1)  !!(set ll corner)
            endif
            if (ipy/=nprocy-1 .or. bcy2(j)=='p') then
              f(:,m2i:m2, n1:n1i,j)=f(:,m2i:m2, n1:n1i,j) + &
                                   ulbufi(:,:,:,j-ivar1+1)  !!(set ul corner)
            endif
          endif
          if (ipz/=nprocz-1 .or. bcz2(j)=='p') then
            if (ipy/=nprocy-1 .or. bcy2(j)=='p') then
              f(:,m2i:m2,n2i:n2,j)=f(:,m2i:m2,n2i:n2,j) + &
                                   uubufi(:,:,:,j-ivar1+1)  !!(set uu corner)
            endif
            if (ipy/=0 .or. bcy1(j)=='p') then
              f(:, m1:m1i,n2i:n2,j)=f(:, m1:m1i,n2i:n2,j) + &
                                   lubufi(:,:,:,j-ivar1+1)  !!(set lu corner)
            endif
          endif
        enddo
        call MPI_WAIT(isend_rq_TOll,isend_stat_Tll,ierr)
        call MPI_WAIT(isend_rq_TOul,isend_stat_Tul,ierr)
        call MPI_WAIT(isend_rq_TOuu,isend_stat_Tuu,ierr)
        call MPI_WAIT(isend_rq_TOlu,isend_stat_Tlu,ierr)
      endif
!
!  communication sample
!  (commented out, because compiler does like this for 0-D runs)
!
!     if (ip<7.and.ipy==3.and.ipz==0) &
!       print*,'finalize_isendrcv_bdry: MPICOMM recv ul: ', &
!                       iproc,ulbufi(nx/2+4,:,1,2),' from ',ulcorn
!
!  make sure the other precessors don't carry on sending new data
!  which could be mistaken for an earlier time
!
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!
    endsubroutine finalize_isendrcv_shockbdry
!***********************************************************************
   subroutine initiate_shearing(f,ivar1_opt,ivar2_opt)
!
!  Subroutine for shearing sheet boundary conditions
!
!  20-june-02/nils: adapted from pencil_mpi
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      double precision :: deltay_dy, frak, c1, c2, c3, c4, c5, c6
      integer :: ivar1, ivar2, ystep, nbufx_gh
      integer :: tolastya=11, tolastyb=12, tonextya=13, tonextyb=14
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  Sixth order interpolation along the y-direction
!
      deltay_dy=deltay/dy
      displs=int(deltay_dy)
      if (nprocy==1) then
        frak=deltay_dy-displs
        c1 = -          (frak+1.)*frak*(frak-1.)*(frak-2.)*(frak-3.)/120.
        c2 = +(frak+2.)          *frak*(frak-1.)*(frak-2.)*(frak-3.)/24.
        c3 = -(frak+2.)*(frak+1.)     *(frak-1.)*(frak-2.)*(frak-3.)/12.
        c4 = +(frak+2.)*(frak+1.)*frak          *(frak-2.)*(frak-3.)/12.
        c5 = -(frak+2.)*(frak+1.)*frak*(frak-1.)          *(frak-3.)/24.
        c6 = +(frak+2.)*(frak+1.)*frak*(frak-1.)*(frak-2.)          /120.
        f(1:l1-1,m1:m2,:,ivar1:ivar2) = &
             c1*cshift(f(l2i:l2,m1:m2,:,ivar1:ivar2),-displs+2,2) &
            +c2*cshift(f(l2i:l2,m1:m2,:,ivar1:ivar2),-displs+1,2) &
            +c3*cshift(f(l2i:l2,m1:m2,:,ivar1:ivar2),-displs  ,2) &
            +c4*cshift(f(l2i:l2,m1:m2,:,ivar1:ivar2),-displs-1,2) &
            +c5*cshift(f(l2i:l2,m1:m2,:,ivar1:ivar2),-displs-2,2) &
            +c6*cshift(f(l2i:l2,m1:m2,:,ivar1:ivar2),-displs-3,2)
        f(l2+1:mx,m1:m2,:,ivar1:ivar2) = &
             c1*cshift(f(l1:l1i,m1:m2,:,ivar1:ivar2), displs-2,2) &
            +c2*cshift(f(l1:l1i,m1:m2,:,ivar1:ivar2), displs-1,2) &
            +c3*cshift(f(l1:l1i,m1:m2,:,ivar1:ivar2), displs  ,2) &
            +c4*cshift(f(l1:l1i,m1:m2,:,ivar1:ivar2), displs+1,2) &
            +c5*cshift(f(l1:l1i,m1:m2,:,ivar1:ivar2), displs+2,2) &
            +c6*cshift(f(l1:l1i,m1:m2,:,ivar1:ivar2), displs+3,2)
      else
!
!  With more than one CPU in the y-direction it will become necessary to
!  interpolate over data from two different CPUs.  Likewise two different
!  CPUs will require data from this CPU.
!
        ystep = displs/ny
        nextya = ipz*nprocy+modulo(ipy-ystep,nprocy)
        lastya = ipz*nprocy+modulo(ipy-ystep-1,nprocy)
        lastyb = ipz*nprocy+modulo(ipy+ystep,nprocy)
        nextyb = ipz*nprocy+modulo(ipy+ystep+1,nprocy)
        fao(:,:,:,ivar1:ivar2) = f(l1:l1i,:,:,ivar1:ivar2)
        fbo(:,:,:,ivar1:ivar2) = f(l2i:l2,:,:,ivar1:ivar2)
        nbufx_gh=my*mz*nghost*(ivar2-ivar1+1)
        if (lastya/=iproc) then
          call MPI_ISEND(fao(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,lastya, &
              tonextyb,MPI_COMM_WORLD,isend_rq_tolastya,ierr)
        endif
        if (nextyb==iproc) then
          fbhi(:,:,:,ivar1:ivar2)=fao(:,:,:,ivar1:ivar2)
        else
          call MPI_IRECV(fbhi(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,nextyb, &
              tonextyb,MPI_COMM_WORLD,irecv_rq_fromnextyb,ierr)
        endif
        if (nextya/=iproc) then
          call MPI_ISEND(fao(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,nextya, &
              tolastyb,MPI_COMM_WORLD,isend_rq_tonextya,ierr)
        endif
        if (lastyb==iproc) then
          fblo(:,:,:,ivar1:ivar2)=fao(:,:,:,ivar1:ivar2)
        else
          call MPI_IRECV(fblo(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,lastyb, &
              tolastyb,MPI_COMM_WORLD,irecv_rq_fromlastyb,ierr)
        endif
        if (lastyb/=iproc) then
          call MPI_ISEND(fbo(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,lastyb, &
              tonextya,MPI_COMM_WORLD,isend_rq_tolastyb,ierr)
        endif
        if (nextya==iproc) then
          fahi(:,:,:,ivar1:ivar2)=fbo(:,:,:,ivar1:ivar2)
        else
          call MPI_IRECV(fahi(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,nextya, &
              tonextya,MPI_COMM_WORLD,irecv_rq_fromnextya,ierr)
        endif
        if (nextyb/=iproc) then
          call MPI_ISEND(fbo(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,nextyb, &
              tolastya,MPI_COMM_WORLD,isend_rq_tonextyb,ierr)
        endif
        if (lastya==iproc) then
          falo(:,:,:,ivar1:ivar2)=fbo(:,:,:,ivar1:ivar2)
        else
          call MPI_IRECV(falo(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,lastya, &
              tolastya,MPI_COMM_WORLD,irecv_rq_fromlastya,ierr)
        endif
      endif
!
    endsubroutine initiate_shearing
!***********************************************************************
    subroutine finalize_shearing(f,ivar1_opt,ivar2_opt)
!
!  Subroutine for shearing sheet boundary conditions
!
!  20-june-02/nils: adapted from pencil_mpi
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      real, dimension (nghost,2*my-2*nghost,mz,mcom) :: fa, fb
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fal, irecv_stat_fan
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fbl, irecv_stat_fbn
      integer, dimension (MPI_STATUS_SIZE) :: isend_stat_tna, isend_stat_tla
      integer, dimension (MPI_STATUS_SIZE) :: isend_stat_tnb, isend_stat_tlb
      integer :: ivar1, ivar2, m2long, i
      double precision :: deltay_dy, frak, c1, c2, c3, c4, c5, c6
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  Sliding periodic boundary conditions in x
!  ulf:02-mar-02
!
! need to wait till all communication has been recived
!
      if (lastyb/=iproc) call MPI_WAIT(irecv_rq_fromlastyb,irecv_stat_fbl,ierr)
      if (nextyb/=iproc) call MPI_WAIT(irecv_rq_fromnextyb,irecv_stat_fbn,ierr)
      if (lastya/=iproc) call MPI_WAIT(irecv_rq_fromlastya,irecv_stat_fal,ierr)
      if (nextya/=iproc) call MPI_WAIT(irecv_rq_fromnextya,irecv_stat_fan,ierr)
!
! reading communicated information into f
!
      deltay_dy=deltay/dy
      m2long = 2*my-3*nghost
      fa(:,1:m2,:,ivar1:ivar2) = falo(:,1:m2,:,ivar1:ivar2)
      fa(:,m2+1:2*my-2*nghost,:,ivar1:ivar2) = fahi(:,m1:my,:,ivar1:ivar2)
      fb(:,1:m2,:,ivar1:ivar2) = fblo(:,1:m2,:,ivar1:ivar2)
      fb(:,m2+1:2*my-2*nghost,:,ivar1:ivar2) = fbhi(:,m1:my,:,ivar1:ivar2)
      displs = modulo(int(deltay_dy),ny)
      frak = deltay_dy - int(deltay_dy)
      c1 = -          (frak+1.)*frak*(frak-1.)*(frak-2.)*(frak-3.)/120.
      c2 = +(frak+2.)          *frak*(frak-1.)*(frak-2.)*(frak-3.)/24.
      c3 = -(frak+2.)*(frak+1.)     *(frak-1.)*(frak-2.)*(frak-3.)/12.
      c4 = +(frak+2.)*(frak+1.)*frak          *(frak-2.)*(frak-3.)/12.
      c5 = -(frak+2.)*(frak+1.)*frak*(frak-1.)          *(frak-3.)/24.
      c6 = +(frak+2.)*(frak+1.)*frak*(frak-1.)*(frak-2.)          /120.
      f(1:l1-1,m1:m2,:,ivar1:ivar2) = &
           c1*fa(:,m2long-ny-displs+3:m2long-displs+2,:,ivar1:ivar2) &
          +c2*fa(:,m2long-ny-displs+2:m2long-displs+1,:,ivar1:ivar2) &
          +c3*fa(:,m2long-ny-displs+1:m2long-displs-0,:,ivar1:ivar2) &
          +c4*fa(:,m2long-ny-displs-0:m2long-displs-1,:,ivar1:ivar2) &
          +c5*fa(:,m2long-ny-displs-1:m2long-displs-2,:,ivar1:ivar2) &
          +c6*fa(:,m2long-ny-displs-2:m2long-displs-3,:,ivar1:ivar2)
      f(l2+1:mx,m1:m2,:,ivar1:ivar2)= &
           c1*fb(:,m1+displs-2:m2+displs-2,:,ivar1:ivar2) &
          +c2*fb(:,m1+displs-1:m2+displs-1,:,ivar1:ivar2) &
          +c3*fb(:,m1+displs  :m2+displs  ,:,ivar1:ivar2) &
          +c4*fb(:,m1+displs+1:m2+displs+1,:,ivar1:ivar2) &
          +c5*fb(:,m1+displs+2:m2+displs+2,:,ivar1:ivar2) &
          +c6*fb(:,m1+displs+3:m2+displs+3,:,ivar1:ivar2)
!
!  Filling also the x-y corners in order to avoid only zeros at these corners.
!  One should acctually have communicated with an extra processor in order to
!  fill these corners with the right values, but this does not seem to be
!  necessary.
!
!      do i=1,nghost
!        f(1:l1-1 ,i   ,:,ivar1:ivar2)=f(1:l1-1 ,m1+nghost-i,:,ivar1:ivar2)
!        f(1:l1-1 ,m2+i,:,ivar1:ivar2)=f(1:l1-1 ,m2-i+1     ,:,ivar1:ivar2)
!        f(l2+1:mx,i   ,:,ivar1:ivar2)=f(l2+1:mx,m1+nghost-i,:,ivar1:ivar2)
!        f(l2+1:mx,m2+i,:,ivar1:ivar2)=f(l2+1:mx,m2-i+1     ,:,ivar1:ivar2)
!      enddo
!
!  need to wait till buffer is empty before re-using it again
!
      if (nextyb/=iproc) call MPI_WAIT(isend_rq_tonextyb,isend_stat_tnb,ierr)
      if (lastyb/=iproc) call MPI_WAIT(isend_rq_tolastyb,isend_stat_tlb,ierr)
      if (nextya/=iproc) call MPI_WAIT(isend_rq_tonextya,isend_stat_tna,ierr)
      if (lastya/=iproc) call MPI_WAIT(isend_rq_tolastya,isend_stat_tla,ierr)
!
    endsubroutine finalize_shearing
!***********************************************************************
    subroutine radboundary_zx_recv(mrad,idir,Qrecv_zx)
!
!  receive intensities from neighboring processor in y
!
!  11-jul-03/tobi: coded
!  20-jul-05/tobi: use non-blocking MPI calls
!
      integer :: mrad,idir
      real, dimension(mx,mz) :: Qrecv_zx
      integer :: isource
      integer, dimension(MPI_STATUS_SIZE) :: irecv_zx
!
!  Identifier
!
      if (lroot.and.ip<5) print*,'radboundary_zx_recv: ENTER'
!
!  source
!
      if (mrad>0) isource=ylneigh
      if (mrad<0) isource=yuneigh
!
!  actual MPI call
!
      call MPI_RECV(Qrecv_zx,mx*mz,MPI_REAL,isource,Qtag_zx+idir, &
                    MPI_COMM_WORLD,irecv_zx,ierr)
!
    endsubroutine radboundary_zx_recv
!***********************************************************************
    subroutine radboundary_xy_recv(nrad,idir,Qrecv_xy)
!
!  receive intensities from neighboring processor in z
!
!  11-jul-03/tobi: coded
!  20-jul-05/tobi: use non-blocking MPI calls
!
      integer :: nrad,idir
      real, dimension(mx,my) :: Qrecv_xy
      integer :: isource
      integer, dimension(MPI_STATUS_SIZE) :: irecv_xy
!
!  Identifier
!
      if (lroot.and.ip<5) print*,'radboundary_xy_recv: ENTER'
!
!  source
!
      if (nrad>0) isource=zlneigh
      if (nrad<0) isource=zuneigh
!
!  actual MPI call
!
      call MPI_RECV(Qrecv_xy,mx*my,MPI_REAL,isource,Qtag_xy+idir, &
                    MPI_COMM_WORLD,irecv_xy,ierr)
!
    endsubroutine radboundary_xy_recv
!***********************************************************************
    subroutine radboundary_zx_send(mrad,idir,Qsend_zx)
!
!  send intensities to neighboring processor in y
!
!  11-jul-03/tobi: coded
!  20-jul-05/tobi: use non-blocking MPI calls
!
      integer :: mrad,idir
      real, dimension(mx,mz) :: Qsend_zx
      integer :: idest
      integer, dimension(MPI_STATUS_SIZE) :: isend_zx
!
!  Identifier
!
      if (lroot.and.ip<5) print*,'radboundary_zx_send: ENTER'
!
!  destination
!
      if (mrad>0) idest=yuneigh
      if (mrad<0) idest=ylneigh
!
!  actual MPI call
!
      call MPI_SEND(Qsend_zx,mx*mz,MPI_REAL,idest,Qtag_zx+idir, &
                    MPI_COMM_WORLD,isend_zx,ierr)
!
    endsubroutine radboundary_zx_send
!***********************************************************************
    subroutine radboundary_xy_send(nrad,idir,Qsend_xy)
!
!  send intensities to neighboring processor in z
!
!  11-jul-03/tobi: coded
!  20-jul-05/tobi: use non-blocking MPI calls
!
      integer, intent(in) :: nrad,idir
      real, dimension(mx,my) :: Qsend_xy
      integer :: idest
      integer, dimension(MPI_STATUS_SIZE) :: isend_xy
!
!  Identifier
!
      if (lroot.and.ip<5) print*,'radboundary_xy_send: ENTER'
!
!  destination
!
      if (nrad>0) idest=zuneigh
      if (nrad<0) idest=zlneigh
!
!  actual MPI call
!
      call MPI_SEND(Qsend_xy,mx*my,MPI_REAL,idest,Qtag_xy+idir, &
                    MPI_COMM_WORLD,isend_xy,ierr)
!
    endsubroutine radboundary_xy_send
!***********************************************************************
    subroutine radboundary_zx_sendrecv(mrad,idir,Qsend_zx,Qrecv_zx)
!
!  receive intensities from isource and send intensities to idest
!
!  04-aug-03/tobi: coded
!
      integer, intent(in) :: mrad,idir
      real, dimension(mx,mz) :: Qsend_zx,Qrecv_zx
      integer :: idest,isource
      integer, dimension(MPI_STATUS_SIZE) :: isendrecv_zx
!
!  Identifier
!
      if (lroot.and.ip<5) print*,'radboundary_zx_sendrecv: ENTER'
!
!  destination and source id
!
      if (mrad>0) then; idest=yuneigh; isource=ylneigh; endif
      if (mrad<0) then; idest=ylneigh; isource=yuneigh; endif
!
!  actual MPI call
!
      call MPI_SENDRECV(Qsend_zx,mx*mz,MPI_REAL,idest,Qtag_zx+idir, &
                        Qrecv_zx,mx*mz,MPI_REAL,isource,Qtag_zx+idir, &
                        MPI_COMM_WORLD,isendrecv_zx,ierr)

    endsubroutine radboundary_zx_sendrecv
!***********************************************************************
    subroutine radboundary_zx_periodic_ray(Qrad_zx,tau_zx, &
                                           Qrad_zx_all,tau_zx_all)
!
!  Gather all intrinsic optical depths and heating rates into one rank-3 array
!  that is available on each processor.
!
!  19-jul-05/tobi: rewritten
!
      real, dimension(nx,nz), intent(in) :: Qrad_zx,tau_zx
      real, dimension(nx,nz,0:nprocy-1), intent(out) :: Qrad_zx_all,tau_zx_all
!
!  Identifier
!
      if (lroot.and.ip<5) print*,'radboundary_zx_periodic_ray: ENTER'
!
!  actual MPI calls
!
      call MPI_ALLGATHER(tau_zx,nx*nz,MPI_REAL,tau_zx_all,nx*nz,MPI_REAL, &
                         MPI_COMM_ROW,ierr)

      call MPI_ALLGATHER(Qrad_zx,nx*nz,MPI_REAL,Qrad_zx_all,nx*nz,MPI_REAL, &
                         MPI_COMM_ROW,ierr)

    endsubroutine radboundary_zx_periodic_ray
!***********************************************************************
    subroutine mpirecv_logical_scl(bcast_array,nbcast_array,proc_src,tag_id)
!
!  Receive logical scalar from other processor.
!
!  04-sep-06/wlad: coded
!
      integer :: nbcast_array
      logical :: bcast_array
      integer :: proc_src, tag_id
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      call MPI_RECV(bcast_array, nbcast_array, MPI_LOGICAL, proc_src, &
          tag_id, MPI_COMM_WORLD, stat, ierr)
!
    endsubroutine mpirecv_logical_scl
!***********************************************************************
    subroutine mpirecv_logical_arr(bcast_array,nbcast_array,proc_src,tag_id)
!
!  Receive logical array from other processor.
!
!  04-sep-06/anders: coded
!
      integer :: nbcast_array
      logical, dimension(nbcast_array) :: bcast_array
      integer :: proc_src, tag_id
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      call MPI_RECV(bcast_array, nbcast_array, MPI_LOGICAL, proc_src, &
          tag_id, MPI_COMM_WORLD, stat, ierr)
!
    endsubroutine mpirecv_logical_arr
!***********************************************************************
    subroutine mpirecv_real_scl(bcast_array,nbcast_array,proc_src,tag_id)
!
!  Receive real scalar from other processor.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      real :: bcast_array
      integer :: proc_src, tag_id
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
!
      call MPI_RECV(bcast_array, nbcast_array, MPI_REAL, proc_src, &
          tag_id, MPI_COMM_WORLD, stat, ierr)
!
    endsubroutine mpirecv_real_scl
!***********************************************************************
    subroutine mpirecv_real_arr(bcast_array,nbcast_array,proc_src,tag_id)
!
!  Receive real array from other processor.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array
      integer :: proc_src, tag_id
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
!
      call MPI_RECV(bcast_array, nbcast_array, MPI_REAL, proc_src, &
          tag_id, MPI_COMM_WORLD, stat, ierr)
!
    endsubroutine mpirecv_real_arr
!***********************************************************************
    subroutine mpirecv_real_arr2(bcast_array,nbcast_array,proc_src,tag_id)
!
!  Receive real array(:,:) from other processor.
!
!  02-jul-05/anders: coded
!
      integer, dimension(2) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer :: proc_src, tag_id, nbcast
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
!
     nbcast=nbcast_array(1)*nbcast_array(2)
!
      call MPI_RECV(bcast_array, nbcast, MPI_REAL, proc_src, &
          tag_id, MPI_COMM_WORLD, stat, ierr)
!
    endsubroutine mpirecv_real_arr2
!***********************************************************************
    subroutine mpirecv_real_arr3(bcast_array,nbcast_array,proc_src,tag_id)
!
!  Receive real array(:,:,:) from other processor.
!
!  20-may-06/anders: adapted
!
      integer, dimension(3) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2), &
                      nbcast_array(3)) :: bcast_array
      integer :: proc_src, tag_id, nbcast
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
!
     nbcast=nbcast_array(1)*nbcast_array(2)*nbcast_array(3)
!
      call MPI_RECV(bcast_array, nbcast, MPI_REAL, proc_src, &
          tag_id, MPI_COMM_WORLD, stat, ierr)
!
    endsubroutine mpirecv_real_arr3
!***********************************************************************
    subroutine mpirecv_real_arr4(bcast_array,nbcast_array,proc_src,tag_id)
!
!  Receive real array(:,:,:,:) from other processor.
!
!  20-may-06/anders: adapted
!
      integer, dimension(4) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2), &
                      nbcast_array(3),nbcast_array(4)) :: bcast_array
      integer :: proc_src, tag_id, nbcast
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
!
      nbcast=nbcast_array(1)*nbcast_array(2)*nbcast_array(3)*nbcast_array(4)
!
      call MPI_RECV(bcast_array, nbcast, MPI_REAL, proc_src, &
          tag_id, MPI_COMM_WORLD, stat, ierr)
!
    endsubroutine mpirecv_real_arr4
!***********************************************************************
    subroutine mpirecv_int_scl(bcast_array,nbcast_array,proc_src,tag_id)
!
!  Receive integer scalar from other processor.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      integer :: bcast_array
      integer :: proc_src, tag_id
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      call MPI_RECV(bcast_array, nbcast_array, MPI_INTEGER, proc_src, &
          tag_id, MPI_COMM_WORLD, stat, ierr)
!
    endsubroutine mpirecv_int_scl
!***********************************************************************
    subroutine mpirecv_int_arr(bcast_array,nbcast_array,proc_src,tag_id)
!
!  Receive integer array from other processor.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: bcast_array
      integer :: proc_src, tag_id
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      call MPI_RECV(bcast_array, nbcast_array, MPI_INTEGER, proc_src, &
          tag_id, MPI_COMM_WORLD, stat, ierr)
!
    endsubroutine mpirecv_int_arr
!***********************************************************************
    subroutine mpirecv_int_arr2(bcast_array,nbcast_array,proc_src,tag_id)
!
!  Receive 2D integer array from other processor.
!
!  20-fev-08/wlad: adpated from mpirecv_real_arr2
!
      integer, dimension(2) :: nbcast_array
      integer, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer :: proc_src, tag_id, nbcast
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
      nbcast = nbcast_array(1)*nbcast_array(2)
!
      call MPI_RECV(bcast_array, nbcast, MPI_INTEGER, proc_src, &
          tag_id, MPI_COMM_WORLD, stat, ierr)
!
    endsubroutine mpirecv_int_arr2
!***********************************************************************
    subroutine mpisend_logical_scl(bcast_array,nbcast_array,proc_rec,tag_id)
!
!  Send logical scalar to other processor.
!
!  04-sep-06/wlad: coded
!
      integer :: nbcast_array
      logical :: bcast_array
      integer :: proc_rec, tag_id
!
      call MPI_SEND(bcast_array, nbcast_array, MPI_LOGICAL, proc_rec, &
          tag_id, MPI_COMM_WORLD, ierr)
!
    endsubroutine mpisend_logical_scl
!***********************************************************************
    subroutine mpisend_logical_arr(bcast_array,nbcast_array,proc_rec,tag_id)
!
!  Send logical array to other processor.
!
!  04-sep-06/wlad: coded
!
      integer :: nbcast_array
      logical, dimension(nbcast_array) :: bcast_array
      integer :: proc_rec, tag_id
!
      call MPI_SEND(bcast_array, nbcast_array, MPI_LOGICAL, proc_rec, &
          tag_id, MPI_COMM_WORLD,ierr)
!
    endsubroutine mpisend_logical_arr
!***********************************************************************
    subroutine mpisend_real_scl(bcast_array,nbcast_array,proc_rec,tag_id)
!
!  Send real scalar to other processor.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      real :: bcast_array
      integer :: proc_rec, tag_id
!
      call MPI_SEND(bcast_array, nbcast_array, MPI_REAL, proc_rec, &
          tag_id, MPI_COMM_WORLD, ierr)
!
    endsubroutine mpisend_real_scl
!***********************************************************************
    subroutine mpisend_real_arr(bcast_array,nbcast_array,proc_rec,tag_id)
!
!  Send real array to other processor.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array
      integer :: proc_rec, tag_id
!
      call MPI_SEND(bcast_array, nbcast_array, MPI_REAL, proc_rec, &
          tag_id, MPI_COMM_WORLD,ierr)
!
    endsubroutine mpisend_real_arr
!***********************************************************************
    subroutine mpisend_real_arr2(bcast_array,nbcast_array,proc_rec,tag_id)
!
!  Send real array(:,:) to other processor.
!
!  02-jul-05/anders: coded
!
      integer, dimension(2) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer :: proc_rec, tag_id, nbcast
!
      nbcast=nbcast_array(1)*nbcast_array(2)
!
      call MPI_SEND(bcast_array, nbcast, MPI_REAL, proc_rec, &
          tag_id, MPI_COMM_WORLD,ierr)
!
    endsubroutine mpisend_real_arr2
!***********************************************************************
    subroutine mpisend_real_arr3(bcast_array,nbcast_array,proc_rec,tag_id)
!
!  Send real array(:,:,:) to other processor.
!
!  20-may-06/anders: adapted
!
      integer, dimension(3) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2), &
                      nbcast_array(3)) :: bcast_array
      integer :: proc_rec, tag_id, nbcast
!
      nbcast=nbcast_array(1)*nbcast_array(2)*nbcast_array(3)
!
      call MPI_SEND(bcast_array, nbcast, MPI_REAL, proc_rec, &
          tag_id, MPI_COMM_WORLD,ierr)
!
    endsubroutine mpisend_real_arr3
!***********************************************************************
    subroutine mpisend_real_arr4(bcast_array,nbcast_array,proc_rec,tag_id)
!
!  Send real array(:,:,:,:) to other processor.
!
!  20-may-06/anders: adapted
!
      integer, dimension(4) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2), &
                      nbcast_array(3),nbcast_array(4)) :: bcast_array
      integer :: proc_rec, tag_id, nbcast
!
      nbcast=nbcast_array(1)*nbcast_array(2)*nbcast_array(3)*nbcast_array(4)
!
      call MPI_SEND(bcast_array, nbcast, MPI_REAL, proc_rec, &
          tag_id, MPI_COMM_WORLD,ierr)
!
    endsubroutine mpisend_real_arr4
!***********************************************************************
    subroutine mpisend_int_scl(bcast_array,nbcast_array,proc_rec,tag_id)
!
!  Send integer scalar to other processor.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      integer :: bcast_array
      integer :: proc_rec, tag_id
!
      call MPI_SEND(bcast_array, nbcast_array, MPI_INTEGER, proc_rec, &
          tag_id, MPI_COMM_WORLD, ierr)
!
    endsubroutine mpisend_int_scl
!***********************************************************************
    subroutine mpisend_int_arr(bcast_array,nbcast_array,proc_rec,tag_id)
!
!  Send integer array to other processor.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: bcast_array
      integer :: proc_rec, tag_id
!
      call MPI_SEND(bcast_array, nbcast_array, MPI_INTEGER, proc_rec, &
          tag_id, MPI_COMM_WORLD,ierr)
!
    endsubroutine mpisend_int_arr
!***********************************************************************
    subroutine mpisend_int_arr2(bcast_array,nbcast_array,proc_rec,tag_id)
!
!  Send 2d integer array to other processor.
!
!  20-fev-08/wlad: adapted from mpisend_real_arr2
!
      integer, dimension(2) :: nbcast_array
      integer, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer :: proc_rec, tag_id, nbcast
!
      nbcast=nbcast_array(1)*nbcast_array(2)
!
      call MPI_SEND(bcast_array, nbcast, MPI_INTEGER, proc_rec, &
          tag_id, MPI_COMM_WORLD,ierr)
!
    endsubroutine mpisend_int_arr2
!***********************************************************************
    subroutine mpibcast_logical_scl(lbcast_array,nbcast_array,proc)
!
!  Communicate logical scalar between processors
!
      integer :: nbcast_array
      logical :: lbcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(lbcast_array,nbcast_array,MPI_LOGICAL,ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_logical_scl
!***********************************************************************
    subroutine mpibcast_logical_arr(lbcast_array,nbcast_array,proc)
!
!  Communicate logical array between processors
!
      integer :: nbcast_array
      logical, dimension (nbcast_array) :: lbcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(lbcast_array,nbcast_array,MPI_LOGICAL,ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_logical_arr
!***********************************************************************
    subroutine mpibcast_logical_arr2(lbcast_array,nbcast_array,proc)
!
!  Communicate logical array(:,:) to other processor.
!
!  25-may-08/wlad: adapted
!
      integer, dimension(2) :: nbcast_array
      logical, dimension(nbcast_array(1),nbcast_array(2)) :: lbcast_array
      integer, optional :: proc
      integer :: ibcast_proc,nbcast
!
      nbcast=nbcast_array(1)*nbcast_array(2)
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(lbcast_array, nbcast, MPI_LOGICAL, ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_logical_arr2
!***********************************************************************
    subroutine mpibcast_int_scl(ibcast_array,nbcast_array,proc)
!
!  Communicate integer scalar between processors
!
      integer :: nbcast_array
      integer :: ibcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(ibcast_array,1,MPI_INTEGER,ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_int_scl
!***********************************************************************
    subroutine mpibcast_int_arr(ibcast_array,nbcast_array,proc)
!
!  Communicate integer array between processors
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: ibcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(ibcast_array,nbcast_array,MPI_INTEGER,ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_int_arr
!***********************************************************************
    subroutine mpibcast_real_scl(bcast_array,nbcast_array,proc)
!
!  Communicate real scalar between processors
!
      integer :: nbcast_array
      real :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(bcast_array,nbcast_array,MPI_REAL,ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_real_scl
!***********************************************************************
    subroutine mpibcast_real_arr(bcast_array,nbcast_array,proc)
!
!  Communicate real array between processors
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(bcast_array,nbcast_array,MPI_REAL,ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_real_arr
!***********************************************************************
    subroutine mpibcast_real_arr2(bcast_array,nbcast_array,proc)
!
!  Communicate real array(:,:) to other processor.
!
!  25-fev-08/wlad: adapted
!
      integer, dimension(2) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc,nbcast
!
      nbcast=nbcast_array(1)*nbcast_array(2)
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(bcast_array, nbcast, MPI_REAL, ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_real_arr2
!***********************************************************************
    subroutine mpibcast_real_arr3(bcast_array,nb,proc)
!
!  Communicate real array(:,:) to other processor.
!
!  25-fev-08/wlad: adapted
!
      integer, dimension(3) :: nb
      real, dimension(nb(1),nb(2),nb(3)) :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc,nbcast
!
      nbcast=nb(1)*nb(2)*nb(3)
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(bcast_array, nbcast, MPI_REAL, ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_real_arr3
!***********************************************************************
    subroutine mpibcast_double_scl(bcast_array,nbcast_array,proc)
!
!  Communicate real scalar between processors
!
      integer :: nbcast_array
      double precision :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(bcast_array,nbcast_array,MPI_DOUBLE_PRECISION,ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_double_scl
!***********************************************************************
    subroutine mpibcast_double_arr(bcast_array,nbcast_array,proc)
!
!  Communicate real array between processors
!
      integer :: nbcast_array
      double precision, dimension(nbcast_array) :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(bcast_array,nbcast_array,MPI_DOUBLE_PRECISION,ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_double_arr
!***********************************************************************
    subroutine mpibcast_char_scl(cbcast_array,nbcast_array,proc)
!
!  Communicate character scalar between processors
!
      integer :: nbcast_array
      character :: cbcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(cbcast_array,nbcast_array,MPI_CHARACTER,ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_char_scl
!***********************************************************************
    subroutine mpibcast_char_arr(cbcast_array,nbcast_array,proc)
!
!  Communicate character array between processors
!
      integer :: nbcast_array
      character, dimension(nbcast_array) :: cbcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(cbcast_array,nbcast_array,MPI_CHARACTER,ibcast_proc, &
          MPI_COMM_WORLD,ierr)
!
    endsubroutine mpibcast_char_arr
!***********************************************************************
    subroutine mpiallreduce_sum_arr2(fsum_tmp,fsum,nreduce_array,lsumy,lsumz)
!
!  Calculate total sum for each array element and return to 
!  all processors (MPI_COMM_WORLD).
!  Or sum over z and return to the nprocy processors (MPI_COMM_COLUMN). 
!  Or sum over y and return to the nprocz processors (MPI_COMM_ROW).
!
!  23-nov-08/wlad: included the lsumy,lsumz possibilities
!
      integer, dimension(2) :: nreduce_array
      real, dimension(nreduce_array(1),nreduce_array(2)) :: fsum_tmp,fsum
      integer :: nreduce,mpiprocs
      logical, optional :: lsumy,lsumz
!
      nreduce=nreduce_array(1)*nreduce_array(2)
!
      if (present(lsumy)) then
        mpiprocs=MPI_COMM_ROW
      elseif (present(lsumz)) then 
        mpiprocs=MPI_COMM_COLUMN
      else
        mpiprocs=MPI_COMM_WORLD 
      endif
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, nreduce, MPI_REAL, MPI_SUM, &
                      mpiprocs, ierr)
!
    endsubroutine mpiallreduce_sum_arr2
!***********************************************************************
    subroutine mpiallreduce_sum_arr(fsum_tmp,fsum,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fsum_tmp,fsum
!
!  calculate total sum for each array element
!  and return to all processors
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, nreduce, MPI_REAL, MPI_SUM, &
                      MPI_COMM_WORLD, ierr)
!
    endsubroutine mpiallreduce_sum_arr
!***********************************************************************
    subroutine mpiallreduce_sum_scl(fsum_tmp,fsum)
!
      real :: fsum_tmp,fsum
!
!  calculate total sum for each array element
!  and return to all processors
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, 1, MPI_REAL, MPI_SUM, &
                      MPI_COMM_WORLD, ierr)
!
    endsubroutine mpiallreduce_sum_scl
!***********************************************************************
    subroutine mpiallreduce_sum_int_arr(fsum_tmp,fsum,nreduce)
!
      integer :: nreduce
      integer, dimension(nreduce) :: fsum_tmp,fsum
!
!  calculate total sum for each array element
!  and return to all processors
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, nreduce, MPI_INTEGER, MPI_SUM, &
                      MPI_COMM_WORLD, ierr)
!
    endsubroutine mpiallreduce_sum_int_arr
!***********************************************************************
    subroutine mpiallreduce_sum_int_scl(fsum_tmp,fsum)
!
      integer :: fsum_tmp,fsum
!
!  calculate total sum for each array element
!  and return to all processors
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, 1, MPI_INTEGER, MPI_SUM, &
                      MPI_COMM_WORLD, ierr)
!
    endsubroutine mpiallreduce_sum_int_scl
!***********************************************************************
    subroutine mpiallreduce_max_arr(fmax_tmp,fmax,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fmax_tmp,fmax
!
!  calculate total maximum for each array element and return to root
!
      call MPI_ALLREDUCE(fmax_tmp, fmax, nreduce, MPI_REAL, MPI_MAX, &
                      MPI_COMM_WORLD, ierr)
!
    endsubroutine mpiallreduce_max_arr
!***********************************************************************
    subroutine mpiallreduce_max_scl(fmax_tmp,fmax)
!
      real :: fmax_tmp,fmax
!
!  calculate total maximum for each array element and return to root
!
      call MPI_ALLREDUCE(fmax_tmp, fmax, 1, MPI_REAL, MPI_MAX, &
                      MPI_COMM_WORLD, ierr)
!
    endsubroutine mpiallreduce_max_scl
!***********************************************************************
    subroutine mpireduce_max_arr(fmax_tmp,fmax,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fmax_tmp,fmax
!
!  calculate total maximum for each array element and return to root
!
      call MPI_REDUCE(fmax_tmp, fmax, nreduce, MPI_REAL, MPI_MAX, root, &
                      MPI_COMM_WORLD, ierr)
    endsubroutine mpireduce_max_arr
!***********************************************************************
    subroutine mpireduce_max_scl(fmax_tmp,fmax)
!
      real :: fmax_tmp,fmax
!
!  calculate total maximum for each array element and return to root
!
      call MPI_REDUCE(fmax_tmp, fmax, 1, MPI_REAL, MPI_MAX, root, &
                      MPI_COMM_WORLD, ierr)
    endsubroutine mpireduce_max_scl
!***********************************************************************
    subroutine mpireduce_min_arr(fmin_tmp,fmin,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fmin_tmp,fmin
!
!  calculate total maximum for each array element and return to root
!
      call MPI_REDUCE(fmin_tmp, fmin, nreduce, MPI_REAL, MPI_MIN, root, &
                      MPI_COMM_WORLD, ierr)
    endsubroutine mpireduce_min_arr
!***********************************************************************
    subroutine mpireduce_min_scl(fmin_tmp,fmin)
!
      real :: fmin_tmp,fmin
!
!  calculate total minimum for each array element and return to root
!
      call MPI_REDUCE(fmin_tmp, fmin, 1, MPI_REAL, MPI_MIN, root, &
                      MPI_COMM_WORLD, ierr)
    endsubroutine mpireduce_min_scl
!***********************************************************************
    subroutine mpireduce_sum(fsum_tmp,fsum,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fsum_tmp,fsum
!
      intent(in)  :: fsum_tmp,nreduce
      intent(out) :: fsum
!
!  calculate total sum for each array element and return to root
!  Unlike for MPI_MAX, for MPI_SUM cannot handle nprocs=1 correctly!
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, nreduce, MPI_REAL, MPI_SUM, root, &
                        MPI_COMM_WORLD, ierr)
      endif
!
    endsubroutine mpireduce_sum
!***********************************************************************
    subroutine mpireduce_sum_scl(fsum_tmp,fsum)
!
      real :: fsum_tmp,fsum
!
!  calculate total sum for each array element and return to root
!  Unlike for MPI_MAX, for MPI_SUM cannot handle nprocs=1 correctly!
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, 1, MPI_REAL, MPI_SUM, root, &
                        MPI_COMM_WORLD, ierr)
      endif
!
    endsubroutine mpireduce_sum_scl
!***********************************************************************
    subroutine mpireduce_sum_double_arr(dsum_tmp,dsum,nreduce)
!
      integer :: nreduce
      double precision, dimension(nreduce) :: dsum_tmp,dsum
!
!  calculate total sum for each array element and return to root
!  Unlike for MPI_MAX, for MPI_SUM cannot handle nprocs=1 correctly!
!
      if (nprocs==1) then
        dsum=dsum_tmp
      else
        call MPI_REDUCE(dsum_tmp, dsum, nreduce, MPI_DOUBLE_PRECISION, MPI_SUM, root, &
                        MPI_COMM_WORLD, ierr)
      endif
!
    endsubroutine mpireduce_sum_double_arr
!***********************************************************************
    subroutine mpireduce_sum_double_scl(dsum_tmp,dsum)
!
      double precision :: dsum_tmp,dsum
!
!  calculate total sum for each array element and return to root
!  Unlike for MPI_MAX, for MPI_SUM cannot handle nprocs=1 correctly!
!
      if (nprocs==1) then
        dsum=dsum_tmp
      else
        call MPI_REDUCE(dsum_tmp, dsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, root, &
                        MPI_COMM_WORLD, ierr)
      endif
!
    endsubroutine mpireduce_sum_double_scl
!***********************************************************************
    subroutine mpireduce_sum_int_arr(fsum_tmp,fsum,nreduce)
!
!  12-jan-05/anders: coded
!
      integer :: nreduce
      integer, dimension(nreduce) :: fsum_tmp,fsum
!
!  calculate total sum for each array element and return to root
!  Unlike for MPI_MAX, for MPI_SUM cannot handle nprocs=1 correctly!
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, nreduce, MPI_INTEGER, MPI_SUM, root, &
                        MPI_COMM_WORLD, ierr)
      endif
!
    endsubroutine mpireduce_sum_int_arr
!***********************************************************************
    subroutine mpireduce_sum_int_scl(fsum_tmp,fsum)
!
!  16-sep-05/anders: adapted from mpireduce_sum_int
!
      integer :: fsum_tmp,fsum
!
!  Calculate sum and return to root.
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, 1, MPI_INTEGER, MPI_SUM, root, &
                        MPI_COMM_WORLD, ierr)
      endif
!
    endsubroutine mpireduce_sum_int_scl
!***********************************************************************
    subroutine mpireduce_or_arr(flor_tmp,flor,nreduce)
!
!  17-sep-05/anders: coded
!
      integer :: nreduce
      logical, dimension(nreduce) :: flor_tmp, flor
!
!  Calculate logical or over all procs and return to root.
!
      if (nprocs==1) then
        flor=flor_tmp
      else
        call MPI_REDUCE(flor_tmp, flor, nreduce, MPI_LOGICAL, MPI_LOR, root, &
                        MPI_COMM_WORLD, ierr)
      endif
!
    endsubroutine mpireduce_or_arr
!***********************************************************************
    subroutine mpireduce_or_scl(flor_tmp,flor)
!
!  17-sep-05/anders: coded
!
      logical :: flor_tmp, flor
!
!  Calculate logical or over all procs and return to root.
!
      if (nprocs==1) then
        flor=flor_tmp
      else
        call MPI_REDUCE(flor_tmp, flor, 1, MPI_LOGICAL, MPI_LOR, root, &
                        MPI_COMM_WORLD, ierr)
      endif
!
    endsubroutine mpireduce_or_scl
!***********************************************************************
    subroutine mpireduce_and_arr(fland_tmp,fland,nreduce)
!
!  11-mar-09/anders: coded
!
      integer :: nreduce
      logical, dimension(nreduce) :: fland_tmp, fland
!
!  Calculate logical and over all procs and return to root.
!
      if (nprocs==1) then
        fland=fland_tmp
      else
        call MPI_REDUCE(fland_tmp, fland, nreduce, MPI_LOGICAL, MPI_LAND, root,&
                        MPI_COMM_WORLD, ierr)
      endif
!
    endsubroutine mpireduce_and_arr
!***********************************************************************
    subroutine mpireduce_and_scl(fland_tmp,fland)
!
!  17-sep-05/anders: coded
!
      logical :: fland_tmp, fland
!
!  Calculate logical and over all procs and return to root.
!
      if (nprocs==1) then
        fland=fland_tmp
      else
        call MPI_REDUCE(fland_tmp, fland, 1, MPI_LOGICAL, MPI_LAND, root, &
                        MPI_COMM_WORLD, ierr)
      endif
!
    endsubroutine mpireduce_and_scl
!***********************************************************************
    subroutine start_serialize()
!
!  Do block between start_serialize and end_serialize serially in iproc
!  order. root goes first, then sends proc1 permission, waits for succes,
!  then sends proc2 permisssion, waits for success, etc.
!
!  19-nov-02/wolf: coded
!
      integer :: buf
      integer, dimension(MPI_STATUS_SIZE) :: status
!
      buf = 0
      if (.not. lroot) then     ! root starts, others wait for permission
        call MPI_RECV(buf,1,MPI_INTEGER,root,io_perm,MPI_COMM_WORLD,status,ierr)
      endif
!
    endsubroutine start_serialize
!***********************************************************************
    subroutine end_serialize()
!
!  do block between start_serialize and end_serialize serially in iproc order
!  19-nov-02/wolf: coded
!
      integer :: i,buf
      integer, dimension(MPI_STATUS_SIZE) :: status
!
      buf = 0
      if (lroot) then
        do i=1,ncpus-1            ! send permission, wait for success message
          call MPI_SEND(buf,1,MPI_INTEGER,i,io_perm,MPI_COMM_WORLD,ierr)
          call MPI_RECV(buf,1,MPI_INTEGER,i,io_succ,MPI_COMM_WORLD,status,ierr)
        enddo
      else                  ! tell root we're done
        call MPI_SEND(buf,1,MPI_INTEGER,root,io_succ,MPI_COMM_WORLD,ierr)
      endif
!
    endsubroutine end_serialize
!***********************************************************************
    subroutine mpibarrier()
!
!  synchronize nodes
!  23-jul-2002/wolf: coded
!
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!
    endsubroutine mpibarrier
!***********************************************************************
    subroutine mpifinalize()
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_FINALIZE(ierr)
    endsubroutine mpifinalize
!***********************************************************************
    function mpiwtime()
!
      double precision :: mpiwtime
      double precision :: MPI_WTIME   ! definition needed for mpicomm_ to work
!
      mpiwtime = MPI_WTIME()
    endfunction mpiwtime
!***********************************************************************
    function mpiwtick()
!
      double precision :: mpiwtick
      double precision :: MPI_WTICK   ! definition needed for mpicomm_ to work
!
      mpiwtick = MPI_WTICK()
    endfunction mpiwtick
!***********************************************************************
    subroutine touch_file(fname)
!
!  touch file (used for code locking)
!  25-may-03/axel: coded
!  06-mar-07/wolf: moved here from sub.f90, so we can use it below
!
      character (len=*) :: fname
!
      open(1,FILE=fname)
      close(1)
!
    endsubroutine touch_file
!***********************************************************************
    subroutine die_gracefully()
!
!  Stop having shutdown MPI neatly
!  With at least some MPI implementations, this only stops if all
!  processors agree to call stop_it().
!
!  29-jun-05/tony: coded
!
      !
      !  tell the world something went wrong -- mpirun may not propagate
      !  an error status
      !
      if (lroot) call touch_file('ERROR') 

      call mpifinalize
      STOP 1                    ! Return nonzero exit status
    endsubroutine die_gracefully
!***********************************************************************
    subroutine stop_it(msg)
!
!  Print message and stop.
!  With at least some MPI implementations, this only stops if all
!  processors agree to call stop_it(). To stop (collectively) if only one
!  or a few processors find some condition, use stop_it_if_any().
!  6-nov-01/wolf: coded
!
      character (len=*) :: msg
!
      if (lroot) write(0,'(A,A)') 'STOPPED: ', msg

      call die_gracefully()
!
    endsubroutine stop_it
!***********************************************************************
    subroutine stop_it_if_any(stop_flag,msg)
!
!  Conditionally print message and stop.
!  This works unilaterally, i.e. if STOP_FLAG is true on _any_ processor,
!  we will all stop.
!  22-nov-04/wolf: coded
!
      logical :: stop_flag
      character (len=*) :: msg
      logical :: global_stop_flag
!
!  Get global OR of stop_flag and distribute it, so all processors agree
!  on whether to call stop_it():
!
      call MPI_ALLREDUCE(stop_flag,global_stop_flag,1,MPI_LOGICAL, &
                         MPI_LOR,MPI_COMM_WORLD,ierr)
!
      if (global_stop_flag) call stop_it(msg)
!
    endsubroutine stop_it_if_any
!***********************************************************************
    subroutine check_emergency_brake()
!
!  Check the lemergency_brake flag and stop with any provided
!  message if it is set.
!
!  29-jul-06/tony: coded
!
      logical :: global_stop_flag
!
!  Get global OR of lemergency_brake and distribute it, so all
!  processors agree on whether to call stop_it():
!
      call MPI_ALLREDUCE(lemergency_brake,global_stop_flag,1,MPI_LOGICAL, &
                         MPI_LOR,MPI_COMM_WORLD,ierr)
!
      if (global_stop_flag) call stop_it( &
            "Emergency brake activated. Check for error messages above.")
!
    endsubroutine check_emergency_brake
!***********************************************************************
    subroutine transp(a,var)
!
!  Doing the transpose of information distributed on several processors
!  Used for doing FFTs in the y and z directions.
!  This routine is presently restricted to the case nxgrid=nygrid (if var=y)
!  and nygrid=nzgrid (if var=z)
!
!  03-sep-02/nils: coded
!  26-oct-02/axel: comments added
!   6-jun-03/axel: works now also in 2-D (need only nxgrid=nygrid)
!   5-oct-06/tobi: generalized to nxgrid = n*nygrid
!
! TODO: Implement nxgrid = n*nzgrid
!
      real, dimension(nx,ny,nz) :: a
      character :: var
!
      real, dimension(ny,ny,nz) :: send_buf_y, recv_buf_y
      real, dimension(nz,ny,nz) :: send_buf_z, recv_buf_z
      real, dimension(:,:), allocatable :: tmp
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc_y,recvc_y,sendc_z,recvc_z,px
      integer :: ytag=101,ztag=102,partner,ierr
      integer :: m,n,ibox,ix
!
!  Doing x-y transpose if var='y'
!
      if (var=='y') then
!
        if (mod(nxgrid,nygrid)/=0) then
          print*,'transp: nxgrid needs to be an integer multiple of '//&
                 'nygrid for var==y'
          call stop_it('Inconsistency: mod(nxgrid,nygrid)/=0')
        endif
!
!  Allocate temporary scratch array
!
        allocate (tmp(ny,ny))
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.
!
        sendc_y=ny*ny*nz; recvc_y=sendc_y
!
!  Send information to different processors (x-y transpose)
!  Divide x-range in as many intervals as we have processors in y.
!  The index px counts through all of them; partner is the index of the
!  processor we need to communicate with. Thus, px is the ipy of partner,
!  but at the same time the x index of the given block.
!
!  Example: ipy=1, ipz=0, then partner=0,2,3, ..., nprocy-1.
!
!
!        ... |
!          3 |  D  E  F  /
!          2 |  B  C  /  F'
!  ipy=    1 |  A  /  C' E'
!          0 |  /  A' B' D'
!            +--------------
!        px=    0  1  2  3 ..
!

!        ipy
!         ^
!  C D    |
!  A B    |      --> px
!
!  if ipy=1,px=0, then exchange block A with A' on partner=0
!  if ipy=1,px=2, then exchange block C' with C on partner=2
!
!  if nxgrid is an integer multiple of nygrid, we divide the whole domain
!  into nxgrid/nygrid boxes (indexed by ibox below) of unit aspect ratio
!  (grid point wise) and only transpose within those boxes.
!
!  The following communication patterns is kind of self-regulated. It
!  avoids deadlock, because there will always be at least one matching
!  pair of processors; the others will have to wait until their partner
!  posts the corresponding request.
!    Doing send and recv together for a given pair of processors
!  (although in an order that avoids deadlocking) allows us to get away
!  with only send_buf and recv_buf as buffers
!
        do px=0,nprocy-1
          do ibox=0,nxgrid/nygrid-1
            if (px/=ipy) then
              partner=px+ipz*nprocy ! = iproc + (px-ipy)
              if (ip<=6) print*,'transp: MPICOMM: ipy,ipz,px,partner=',ipy,ipz,px,partner
              ix=ibox*nprocy*ny+px*ny
              send_buf_y=a(ix+1:ix+ny,:,:)
              if (px<ipy) then      ! above diagonal: send first, receive then
                call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,ierr)
                call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,ierr)
              elseif (px>ipy) then  ! below diagonal: receive first, send then
                call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,ierr)
                call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,ierr)
              endif
              a(ix+1:ix+ny,:,:)=recv_buf_y
            endif
          enddo
        enddo
!
!  Transposing the received data (x-y transpose)
!  Example:
!
!  |12 13 | 14 15|      | 6  7 | 14 15|      | 3  7 | 11 15|
!  | 8  9 | 10 11|      | 2  3 | 10 11|      | 2  6 | 10 14|
!  |------+------|  ->  |------+------|  ->  |------+------|
!  | 4  5 |  6  7|      | 4  5 | 12 13|      | 1  5 |  9 13|
!  | 0  1 |  2  3|      | 0  1 |  8  9|      | 0  4 |  8 12|
!     original          2x2 blocks         each block
!                       transposed         transposed
!
        do px=0,nprocy-1
          do ibox=0,nxgrid/nygrid-1
            ix=ibox*nprocy*ny+px*ny
            do n=1,nz
              tmp=transpose(a(ix+1:ix+ny,:,n)); a(ix+1:ix+ny,:,n)=tmp
            enddo
          enddo
        enddo
!
!  Deallocate temporary scratch array
!
        deallocate (tmp)

!
!  Doing x-z transpose if var='z'
!
      elseif (var=='z') then

        if (nzgrid/=nxgrid) then
          if (lroot) print*, 'transp: need to have nzgrid=nxgrid for var==z'
          call stop_it('transp: inconsistency - nzgrid/=nxgrid')
        endif
!
!  Calculate the size of buffers.
!  Buffers used for the z-transpose have the same size in z and x.
!
        sendc_z=nz*ny*nz; recvc_z=sendc_z
!
!  Allocate temporary scratch array
!
        allocate (tmp(nz,nz))
!
!  Send information to different processors (x-z transpose)
!  See the discussion above for why we use this communication pattern
        do px=0,nprocz-1
          if (px/=ipz) then
            partner=ipy+px*nprocy ! = iproc + (px-ipz)*nprocy
            send_buf_z=a(px*nz+1:(px+1)*nz,:,:)
            if (px<ipz) then      ! above diagonal: send first, receive then
              call MPI_SEND(send_buf_z,sendc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,ierr)
              call MPI_RECV (recv_buf_z,recvc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,stat,ierr)
            elseif (px>ipz) then  ! below diagonal: receive first, send then
              call MPI_RECV (recv_buf_z,recvc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,stat,ierr)
              call MPI_SSEND(send_buf_z,sendc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,ierr)
            endif
            a(px*nz+1:(px+1)*nz,:,:)=recv_buf_z
          endif
        enddo
!
!  Transposing the received data (x-z transpose)
!
        do px=0,nprocz-1
          do m=1,ny
            tmp=transpose(a(px*nz+1:(px+1)*nz,m,:))
            a(px*nz+1:(px+1)*nz,m,:)=tmp
          enddo
        enddo
!
      else
        if (lroot) print*,'transp: No clue what var=', var, 'is supposed to mean'
      endif
!
!  Synchronize; not strictly necessary, so Axel will prabably remove it..
!
      call mpibarrier()
!
    endsubroutine transp
!***********************************************************************
    subroutine transp_xy(a)
!
!  Doing the transpose of information distributed on several processors.
!  This routine transposes 2D arrays in x and y only.
!
!   6-oct-06/tobi: Adapted from transp
!
! TODO: Implement nxgrid = n*nzgrid
!
      real, dimension(nx,ny), intent(inout) :: a
!
      real, dimension(ny,ny) :: send_buf_y, recv_buf_y, tmp
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc_y,recvc_y,px
      integer :: ytag=101,partner,ierr
      integer :: ibox,iy

      if (mod(nxgrid,nygrid)/=0) then
        print*,'transp: nxgrid needs to be an integer multiple of '//&
               'nygrid for var==y'
        call stop_it('Inconsistency: mod(nxgrid,nygrid)/=0')
      endif

!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.

      sendc_y=ny**2
      recvc_y=ny**2
!
!  Send information to different processors (x-y transpose)
!  Divide x-range in as many intervals as we have processors in y.
!  The index px counts through all of them; partner is the index of the
!  processor we need to communicate with. Thus, px is the ipy of partner,
!  but at the same time the x index of the given block.
!
!  Example: ipy=1, ipz=0, then partner=0,2,3, ..., nprocy-1.
!
!
!        ... |
!          3 |  D  E  F  /
!          2 |  B  C  /  F'
!  ipy=    1 |  A  /  C' E'
!          0 |  /  A' B' D'
!            +--------------
!        px=    0  1  2  3 ..
!

!        ipy
!         ^
!  C D    |
!  A B    |      --> px
!
!  if ipy=1,px=0, then exchange block A with A' on partner=0
!  if ipy=1,px=2, then exchange block C' with C on partner=2
!
!  if nxgrid is an integer multiple of nygrid, we divide the whole domain
!  into nxgrid/nygrid boxes (indexed by ibox below) of unit aspect ratio
!  (grid point wise) and only transpose within those boxes.
!
!  The following communication patterns is kind of self-regulated. It
!  avoids deadlock, because there will always be at least one matching
!  pair of processors; the others will have to wait until their partner
!  posts the corresponding request.
!    Doing send and recv together for a given pair of processors
!  (although in an order that avoids deadlocking) allows us to get away
!  with only send_buf and recv_buf as buffers
!
      do px=0,nprocy-1
        do ibox=0,nxgrid/nygrid-1
          if (px/=ipy) then
            partner=px+ipz*nprocy ! = iproc + (px-ipy)
            iy=(ibox*nprocy+px)*ny
            send_buf_y=a(iy+1:iy+ny,:)
            if (px<ipy) then      ! above diagonal: send first, receive then
              call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,ierr)
              call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,ierr)
            elseif (px>ipy) then  ! below diagonal: receive first, send then
              call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,ierr)
              call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,ierr)
            endif
            a(iy+1:iy+ny,:)=recv_buf_y
          endif
        enddo
      enddo
!
!  Transposing the received data (x-y transpose)
!  Example:
!
!  |12 13 | 14 15|      | 6  7 | 14 15|      | 3  7 | 11 15|
!  | 8  9 | 10 11|      | 2  3 | 10 11|      | 2  6 | 10 14|
!  |------+------|  ->  |------+------|  ->  |------+------|
!  | 4  5 |  6  7|      | 4  5 | 12 13|      | 1  5 |  9 13|
!  | 0  1 |  2  3|      | 0  1 |  8  9|      | 0  4 |  8 12|
!     original          2x2 blocks         each block
!                       transposed         transposed
!
      do px=0,nprocy-1
        do ibox=0,nxgrid/nygrid-1
          iy=(ibox*nprocy+px)*ny
          tmp=transpose(a(iy+1:iy+ny,:)); a(iy+1:iy+ny,:)=tmp
        enddo
      enddo

    endsubroutine transp_xy
!***********************************************************************
    subroutine transp_xy_other(a)
!
!  Doing the transpose of information distributed on several processors.
!  This routine transposes 2D arrays of arbitrary size in x and y only.
!
!   6-oct-06/tobi: Adapted from transp
!
! TODO: Implement nxgrid = n*nzgrid
!
      real, dimension(:,:), intent(inout) :: a
!
      real, dimension(size(a,2),size(a,2)) :: send_buf_y, recv_buf_y, tmp
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc_y,recvc_y,px
      integer :: ytag=101,partner,ierr
      integer :: ibox,iy,nx_other,ny_other
      integer :: nxgrid_other,nygrid_other
!
      nx_other=size(a,1); ny_other=size(a,2)
      nxgrid_other=nx_other
      nygrid_other=ny_other*nprocy   
!
      if (mod(nxgrid_other,nygrid_other)/=0) then
        print*,'transp: nxgrid_other needs to be an integer multiple of '//&
               'nygrid_other for var==y'
        call stop_it('Inconsistency: mod(nxgrid_other,nygrid_other)/=0')
      endif
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.

      sendc_y=ny_other**2
      recvc_y=ny_other**2
!
!  Send information to different processors (x-y transpose)
!  Divide x-range in as many intervals as we have processors in y.
!  The index px counts through all of them; partner is the index of the
!  processor we need to communicate with. Thus, px is the ipy of partner,
!  but at the same time the x index of the given block.
!
!  Example: ipy=1, ipz=0, then partner=0,2,3, ..., nprocy-1.
!
!
!        ... |
!          3 |  D  E  F  /
!          2 |  B  C  /  F'
!  ipy=    1 |  A  /  C' E'
!          0 |  /  A' B' D'
!            +--------------
!        px=    0  1  2  3 ..
!

!        ipy
!         ^
!  C D    |
!  A B    |      --> px
!
!  if ipy=1,px=0, then exchange block A with A' on partner=0
!  if ipy=1,px=2, then exchange block C' with C on partner=2
!
!  if nxgrid is an integer multiple of nygrid, we divide the whole domain
!  into nxgrid/nygrid boxes (indexed by ibox below) of unit aspect ratio
!  (grid point wise) and only transpose within those boxes.
!
!  The following communication patterns is kind of self-regulated. It
!  avoids deadlock, because there will always be at least one matching
!  pair of processors; the others will have to wait until their partner
!  posts the corresponding request.
!    Doing send and recv together for a given pair of processors
!  (although in an order that avoids deadlocking) allows us to get away
!  with only send_buf and recv_buf as buffers
!
      do px=0,nprocy-1
        do ibox=0,nxgrid_other/nygrid_other-1
          if (px/=ipy) then
            partner=px+ipz*nprocy ! = iproc + (px-ipy)
            iy=(ibox*nprocy+px)*ny_other
            send_buf_y=a(iy+1:iy+ny_other,:)
            if (px<ipy) then      ! above diagonal: send first, receive then
              call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,ierr)
              call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,ierr)
            elseif (px>ipy) then  ! below diagonal: receive first, send then
              call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,ierr)
              call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,ierr)
            endif
            a(iy+1:iy+ny_other,:)=recv_buf_y
          endif
        enddo
      enddo
!
!  Transposing the received data (x-y transpose)
!  Example:
!
!  |12 13 | 14 15|      | 6  7 | 14 15|      | 3  7 | 11 15|
!  | 8  9 | 10 11|      | 2  3 | 10 11|      | 2  6 | 10 14|
!  |------+------|  ->  |------+------|  ->  |------+------|
!  | 4  5 |  6  7|      | 4  5 | 12 13|      | 1  5 |  9 13|
!  | 0  1 |  2  3|      | 0  1 |  8  9|      | 0  4 |  8 12|
!     original          2x2 blocks         each block
!                       transposed         transposed
!
      do px=0,nprocy-1
        do ibox=0,nxgrid_other/nygrid_other-1
          iy=(ibox*nprocy+px)*ny_other
          tmp=transpose(a(iy+1:iy+ny_other,:)); a(iy+1:iy+ny_other,:)=tmp
        enddo
      enddo

    endsubroutine transp_xy_other
!***********************************************************************
    subroutine transp_other(a,var)
!
!  Doing the transpose of information distributed on several processors.
!  This routine transposes 3D arrays but is presently restricted to the 
!  case nxgrid=nygrid (if var=y) and nygrid=nzgrid (if var=z)
!
!  08-may-08/wlad: Adapted from transp
!
! TODO: Implement nxgrid = n*nzgrid
!
      real, dimension(:,:,:), intent(inout) :: a
      character :: var
!
      real, dimension(size(a,2),size(a,2),size(a,3)) :: send_buf_y, recv_buf_y
      real, dimension(size(a,3),size(a,2),size(a,3)) :: send_buf_z, recv_buf_z
      real, dimension(:,:), allocatable :: tmp
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc_y,recvc_y,sendc_z,recvc_z,px
      integer :: ytag=101,ztag=202,partner,ierr
      integer :: m,n,ibox,ix,nx_other,ny_other,nz_other
      integer :: nxgrid_other,nygrid_other,nzgrid_other
!
      nx_other=size(a,1); ny_other=size(a,2) ; nz_other=size(a,3)
      nxgrid_other=nx_other
      nygrid_other=ny_other*nprocy   
      nzgrid_other=nz_other*nprocz   
!
      if (var=='y') then
!
        if (mod(nxgrid_other,nygrid_other)/=0) then
          print*,'transp: nxgrid_other needs to be an integer multiple of '//&
               'nygrid_other for var==y'
          call stop_it('Inconsistency: mod(nxgrid_other,nygrid_other)/=0')
        endif
!
!  Allocate temporary scratch array
!
        allocate (tmp(ny_other,ny_other))
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.

        sendc_y=ny_other**2*nz_other ; recvc_y=sendc_y
!
!  Send information to different processors (x-y transpose)
!  Divide x-range in as many intervals as we have processors in y.
!  The index px counts through all of them; partner is the index of the
!  processor we need to communicate with. Thus, px is the ipy of partner,
!  but at the same time the x index of the given block.
!
!  Example: ipy=1, ipz=0, then partner=0,2,3, ..., nprocy-1.
!
!
!        ... |
!          3 |  D  E  F  /
!          2 |  B  C  /  F'
!  ipy=    1 |  A  /  C' E'
!          0 |  /  A' B' D'
!            +--------------
!        px=    0  1  2  3 ..
!

!        ipy
!         ^
!  C D    |
!  A B    |      --> px
!
!  if ipy=1,px=0, then exchange block A with A' on partner=0
!  if ipy=1,px=2, then exchange block C' with C on partner=2
!
!  if nxgrid is an integer multiple of nygrid, we divide the whole domain
!  into nxgrid/nygrid boxes (indexed by ibox below) of unit aspect ratio
!  (grid point wise) and only transpose within those boxes.
!
!  The following communication patterns is kind of self-regulated. It
!  avoids deadlock, because there will always be at least one matching
!  pair of processors; the others will have to wait until their partner
!  posts the corresponding request.
!    Doing send and recv together for a given pair of processors
!  (although in an order that avoids deadlocking) allows us to get away
!  with only send_buf and recv_buf as buffers
!
        do px=0,nprocy-1
          do ibox=0,nxgrid_other/nygrid_other-1
            if (px/=ipy) then
              partner=px+ipz*nprocy ! = iproc + (px-ipy)
              ix=(ibox*nprocy+px)*ny_other
              send_buf_y=a(ix+1:ix+ny_other,:,:)
              if (px<ipy) then      ! above diagonal: send first, receive then
                call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,ierr)
                call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,ierr)
              elseif (px>ipy) then  ! below diagonal: receive first, send then
                call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,ierr)
                call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,ierr)
              endif
              a(ix+1:ix+ny_other,:,:)=recv_buf_y
            endif
          enddo
        enddo
!
!  Transposing the received data (x-y transpose)
!  Example:
!
!  |12 13 | 14 15|      | 6  7 | 14 15|      | 3  7 | 11 15|
!  | 8  9 | 10 11|      | 2  3 | 10 11|      | 2  6 | 10 14|
!  |------+------|  ->  |------+------|  ->  |------+------|
!  | 4  5 |  6  7|      | 4  5 | 12 13|      | 1  5 |  9 13|
!  | 0  1 |  2  3|      | 0  1 |  8  9|      | 0  4 |  8 12|
!     original          2x2 blocks         each block
!                       transposed         transposed
!
        do px=0,nprocy-1
          do ibox=0,nxgrid_other/nygrid_other-1
            ix=(ibox*nprocy+px)*ny_other
            do n=1,nz
              tmp=transpose(a(ix+1:ix+ny_other,:,n)); a(ix+1:ix+ny_other,:,n)=tmp
            enddo
          enddo
        enddo
!
!  Deallocate temporary scratch array
!
        deallocate (tmp)

!
!  Doing x-z transpose if var='z'
!
      elseif (var=='z') then

        if (nzgrid_other/=nxgrid_other) then
          if (lroot) print*, 'transp_other: need to have '//&
          'nzgrid_other=nxgrid_other for var==z'
          call stop_it('transp_other: inconsistency - nzgrid/=nxgrid')
        endif
!
!  Calculate the size of buffers.
!  Buffers used for the z-transpose have the same size in z and x.
!
        sendc_z=nz_other**2*ny_other; recvc_z=sendc_z
!
!  Allocate temporary scratch array
!
        allocate (tmp(nz_other,nz_other))
!
!  Send information to different processors (x-z transpose)
!  See the discussion above for why we use this communication pattern
        do px=0,nprocz-1
          if (px/=ipz) then
            partner=ipy+px*nprocy ! = iproc + (px-ipz)*nprocy
            send_buf_z=a(px*nz_other+1:(px+1)*nz_other,:,:)
            if (px<ipz) then      ! above diagonal: send first, receive then
              call MPI_SEND(send_buf_z,sendc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,ierr)
              call MPI_RECV (recv_buf_z,recvc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,stat,ierr)
            elseif (px>ipz) then  ! below diagonal: receive first, send then
              call MPI_RECV (recv_buf_z,recvc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,stat,ierr)
              call MPI_SSEND(send_buf_z,sendc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,ierr)
            endif
            a(px*nz_other+1:(px+1)*nz_other,:,:)=recv_buf_z
          endif
        enddo
!
!  Transposing the received data (x-z transpose)
!
        do px=0,nprocz-1
          do m=1,ny
            tmp=transpose(a(px*nz_other+1:(px+1)*nz_other,m,:))
            a(px*nz_other+1:(px+1)*nz_other,m,:)=tmp
          enddo
        enddo
!
      else
        if (lroot) print*,'transp_other: No clue what var=', var, &
             'is supposed to mean'
      endif
!
!  Synchronize; not strictly necessary, so Axel will prabably remove it..
!
      call mpibarrier()
!
    endsubroutine transp_other
!***********************************************************************
    subroutine transp_xz(a,b)
!
!  Doing the transpose of information distributed on several processors.
!  This routine transposes 2D arrays in x and z only.
!
!  19-dec-06/anders: Adapted from transp
!
      integer, parameter :: nxt=nx/nprocz
      real, dimension(nx,nz), intent(in) :: a
      real, dimension(nzgrid,nxt), intent (out) :: b
!
      real, dimension(nxt,nz) :: send_buf, recv_buf
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc,recvc,px
      integer :: ztag=101,partner,ierr
!
      if (mod(nxgrid,nprocz)/=0) then
        print*,'transp_xz: nxgrid needs to be an integer multiple of nprocz'
        call stop_it('Inconsistency: mod(nxgrid,nprocz)/=0')
      endif
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.
!
      sendc=nx*nz; recvc=sendc
!
!  Send information to different processors (x-z transpose)
!
      b(ipz*nz+1:(ipz+1)*nz,:)=transpose(a(ipz*nxt+1:(ipz+1)*nxt,:))
      do px=0,nprocz-1
        if (px/=ipz) then
          partner=ipy+px*nprocy ! = iproc + (px-ipz)*nprocy
          send_buf=a(px*nxt+1:(px+1)*nxt,:)
          if (px<ipz) then      ! above diagonal: send first, receive then
            call MPI_SEND(send_buf,sendc,MPI_REAL,partner,ztag,MPI_COMM_WORLD,ierr)
            call MPI_RECV(recv_buf,recvc,MPI_REAL,partner,ztag,MPI_COMM_WORLD,stat,ierr)
          elseif (px>ipz) then  ! below diagonal: receive first, send then
            call MPI_RECV(recv_buf,recvc,MPI_REAL,partner,ztag,MPI_COMM_WORLD,stat,ierr)
            call MPI_SEND(send_buf,sendc,MPI_REAL,partner,ztag,MPI_COMM_WORLD,ierr)
          endif
          b(px*nz+1:(px+1)*nz,:)=transpose(recv_buf)
        endif
      enddo
!
    endsubroutine transp_xz
!***********************************************************************
    subroutine transp_zx(a,b)
!
!  Doing the transpose of information distributed on several processors.
!  This routine transposes 2D arrays in x and z only.
!
!  19-dec-06/anders: Adapted from transp
!
      integer, parameter :: nxt=nx/nprocz
      real, dimension(nzgrid,nxt), intent (in) :: a
      real, dimension(nx,nz), intent(out) :: b
!
      real, dimension(nz,nxt) :: send_buf, recv_buf
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc,recvc,px
      integer :: ztag=101,partner,ierr
!
      if (mod(nxgrid,nprocz)/=0) then
        print*,'transp_xz: nxgrid needs to be an integer multiple of nprocz'
        call stop_it('Inconsistency: mod(nxgrid,nprocz)/=0')
      endif
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.
!
      sendc=nz*nxt; recvc=sendc
!
!  Send information to different processors (x-z transpose)
!
      b(ipz*nxt+1:(ipz+1)*nxt,:)=transpose(a(ipz*nz+1:(ipz+1)*nz,:))
      do px=0,nprocz-1
        if (px/=ipz) then
          partner=ipy+px*nprocy ! = iproc + (px-ipz)*nprocy
          send_buf=a(px*nz+1:(px+1)*nz,:)
          if (px<ipz) then      ! above diagonal: send first, receive then
            call MPI_SEND(send_buf,sendc,MPI_REAL,partner,ztag,MPI_COMM_WORLD,ierr)
            call MPI_RECV(recv_buf,recvc,MPI_REAL,partner,ztag,MPI_COMM_WORLD,stat,ierr)
          elseif (px>ipz) then  ! below diagonal: receive first, send then
            call MPI_RECV(recv_buf,recvc,MPI_REAL,partner,ztag,MPI_COMM_WORLD,stat,ierr)
            call MPI_SEND(send_buf,sendc,MPI_REAL,partner,ztag,MPI_COMM_WORLD,ierr)
          endif
          b(px*nxt+1:(px+1)*nxt,:)=transpose(recv_buf)
        endif
      enddo
!
    endsubroutine transp_zx
!***********************************************************************
    subroutine communicate_bc_aa_pot(f,topbot)
!
!  Helper routine for bc_aa_pot in Magnetic.
!  Needed due to Fourier transforms which only work on (l1:l2,m1:m2)
!
!   8-oct-2006/tobi: Coded
!
      use Cdata, only: iax,iaz

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot

      real, dimension (nx,nghost,nghost+1,3) :: lbufyo,ubufyo,lbufyi,ubufyi
      integer :: nbufy,nn1,nn2

      select case (topbot)
        case ('bot'); nn1=1;  nn2=n1
        case ('top'); nn1=n2; nn2=mz
        case default; call stop_it("communicate_bc_aa_pot: "//topbot//&
                                   " should be either `top' or `bot'")
      end select
!
!  Periodic boundaries in y -- communicate along y if necessary
!
      if (nprocy>1) then

        lbufyo = f(l1:l2, m1:m1i,nn1:nn2,iax:iaz)
        ubufyo = f(l1:l2,m2i:m2 ,nn1:nn2,iax:iaz)

        nbufy=nx*nghost*(nghost+1)*3

        call MPI_IRECV(ubufyi,nbufy,MPI_REAL,yuneigh,tolowy, &
                       MPI_COMM_WORLD,irecv_rq_fromuppy,ierr)
        call MPI_IRECV(lbufyi,nbufy,MPI_REAL,ylneigh,touppy, &
                       MPI_COMM_WORLD,irecv_rq_fromlowy,ierr)
        call MPI_ISEND(lbufyo,nbufy,MPI_REAL,ylneigh,tolowy, &
                       MPI_COMM_WORLD,isend_rq_tolowy,ierr)
        call MPI_ISEND(ubufyo,nbufy,MPI_REAL,yuneigh,touppy, &
                       MPI_COMM_WORLD,isend_rq_touppy,ierr)

        call MPI_WAIT(irecv_rq_fromuppy,irecv_stat_fu,ierr)
        call MPI_WAIT(irecv_rq_fromlowy,irecv_stat_fl,ierr)

        f(l1:l2,   1:m1-1,nn1:nn2,iax:iaz) = lbufyi
        f(l1:l2,m2+1:my  ,nn1:nn2,iax:iaz) = ubufyi

        call MPI_WAIT(isend_rq_tolowy,isend_stat_tl,ierr)
        call MPI_WAIT(isend_rq_touppy,isend_stat_tu,ierr)

      else

        f(l1:l2,   1:m1-1,nn1:nn2,iax:iaz) = f(l1:l2,m2i:m2 ,nn1:nn2,iax:iaz)
        f(l1:l2,m2+1:my  ,nn1:nn2,iax:iaz) = f(l1:l2, m1:m1i,nn1:nn2,iax:iaz)

      endif
!
!  Periodic boundaries in x
!
      f(   1:l1-1,:,nn1:nn2,iax:iaz) = f(l2i:l2 ,:,nn1:nn2,iax:iaz)
      f(l2+1:mx  ,:,nn1:nn2,iax:iaz) = f( l1:l1i,:,nn1:nn2,iax:iaz)

    endsubroutine communicate_bc_aa_pot
!***********************************************************************
    subroutine z2x(a,xi,yj,yproc_no,az)
!
!
      real, dimension(nx,ny,nz), intent(in) :: a
      real, dimension(nzgrid), intent(out) :: az
      real, dimension(nzgrid) :: az_local
      integer, intent(in) :: xi,yj,yproc_no
      integer :: my_iniz,my_finz
!
      az=0.
      az_local=0.
      if (ipy==(yproc_no-1)) then
        my_iniz=ipz*nz+1
        my_finz=(ipz+1)*nz
!       write(*,*)'DM',iproc,ipy,yproc_no,ipz,my_iniz,my_finz,nzgrid
        az_local(my_iniz:my_finz)=a(xi,yj,:)
      else
      az_local=0.
      endif
      call mpireduce_sum(az_local,az,nzgrid)
! maybe we should synchrosize here. 
      call mpibarrier
!
    endsubroutine z2x
!***********************************************************************

endmodule Mpicomm
