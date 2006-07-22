! $Id: mpicomm.f90,v 1.176 2006-07-22 16:01:19 ajohan Exp $

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
  use Cdata, only: iproc,ipx,ipy,ipz,root,lroot

  implicit none

  include 'mpicomm.h'

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
  endinterface

  interface mpibcast_logical
    module procedure mpibcast_logical_scl
    module procedure mpibcast_logical_arr
  endinterface

  interface mpibcast_int
    module procedure mpibcast_int_scl
    module procedure mpibcast_int_arr
  endinterface

  interface mpibcast_real
    module procedure mpibcast_real_scl
    module procedure mpibcast_real_arr
  endinterface

  interface mpibcast_double
    module procedure mpibcast_double_scl
    module procedure mpibcast_double_arr
  endinterface

  interface mpibcast_char
    module procedure mpibcast_char_scl
    module procedure mpibcast_char_arr
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
! For f-array processor boundaries
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
  integer :: tolowy=3,touppy=4,tolowz=5,touppz=6 ! msg. tags
  integer :: TOll=7,TOul=8,TOuu=9,TOlu=10 ! msg. tags for corners
  integer :: io_perm=20,io_succ=21
!  mpi tags for radiation
!  the values for those have to differ by a number greater than maxdir=190
!  in order to have unique tags for each boundary and each direction
  integer, parameter :: Qtag_zx=300,Qtag_xy=350
  integer, parameter :: tautag_zx=400,tautag_xy=450
  integer, parameter :: Qtag_peri_zx=1000,Qtag_peri_xy=2000
  integer, parameter :: tautag_peri_zx=3000,tautag_peri_xy=4000
!
!  Communicators
!
  integer :: MPI_COMM_ROW
!
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

  integer :: ylneigh,zlneigh ! `lower' neighbours
  integer :: yuneigh,zuneigh ! `upper' neighbours
  integer :: llcorn,lucorn,uucorn,ulcorn !!(the 4 corners in yz-plane)

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
      use General
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
      if (nprocx /= 1) &
           call stop_it('Inconsistency: nprocx > 1 not implemented')
      if (nprocs /= nprocy*nprocz) then
        if (lroot) then
          print*, 'Compiled with NCPUS = ', ncpus, &
          ', but running on ', nprocs, ' processors'
        endif
        call stop_it('Inconsistency 1')
      endif
!
      if ((nprocy*ny /= nygrid) .or. &
          (nprocz*nz /= nzgrid)) then
        if (lroot) then
          write(0,'(A,2I4,A,2I4,A)') &
               'nproc[y-z]*n[y-z] = (', &
               nprocy*ny, nprocz*nz, &
               ') /= n[yz]grid= (', &
               nygrid, nzgrid, ")"
        endif
        call stop_it('Inconsistency 2')
      endif
!
      if ((ny<nghost) .and. (nygrid/=1)) &
           call stop_it('Overlapping ghost zones in y-direction: reduce nprocy')
      if ((nz<nghost) .and. (nzgrid/=1)) &
           call stop_it('Overlapping ghost zones in z-direction: reduce nprocz')
!
!  position on the processor grid
!  x is fastest direction, z slowest
!
      if (lprocz_slowest) then
        ipx = 0
        ipy = modulo(iproc, nprocy)
        ipz = iproc/(nprocy)
       else
        ipx = 0
        ipy = iproc/(nprocz)
        ipz = modulo(iproc, nprocz)
       endif
!
!  set up `lower' and `upper' neighbours
!
      ylneigh = (ipz*nprocy+modulo(ipy-1,nprocy))
      yuneigh = (ipz*nprocy+modulo(ipy+1,nprocy))
      zlneigh = (modulo(ipz-1,nprocz)*nprocy+ipy)
      zuneigh = (modulo(ipz+1,nprocz)*nprocy+ipy)
!
!  set the four corners in the yz-plane (in cyclic order)
!
      llcorn=modulo(ipy-1,nprocy)+modulo(ipz-1,nprocz)*nprocy
      ulcorn=modulo(ipy+1,nprocy)+modulo(ipz-1,nprocz)*nprocy
      uucorn=modulo(ipy+1,nprocy)+modulo(ipz+1,nprocz)*nprocy
      lucorn=modulo(ipy-1,nprocy)+modulo(ipz+1,nprocz)*nprocy
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
      if (ip<5) &
           write(*,'(A,I4,"(",2I4,"): ",8I4)') &
           'mpicomm_init: MPICOMM neighbors ', &
           iproc,ipy,ipz, &
           ylneigh,llcorn,zlneigh,ulcorn,yuneigh,uucorn,zuneigh,lucorn
!
      call setup_mm_nn()
!
!  Define MPI communicator MPI_COMM_ROW that includes all processes
!  sharing the same value of ipz. The rank within MPI_COMM_WORLD is
!  given by ipy.
!
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, ipz, ipy, MPI_COMM_ROW, ierr)
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
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbufy, nbufz, nbufyz
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  So far no distribution over x
!
      !(We may never do this after having seen all the trouble involved!!)
!
!  Periodic boundary conditions in y
!
      if (nprocy>1) then
        lbufyo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n1:n2,ivar1:ivar2) !!(lower y-zone)
        ubufyo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n1:n2,ivar1:ivar2) !!(upper y-zone)
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
        lbufzo(:,:,:,ivar1:ivar2)=f(:,m1:m2,n1:n1i,ivar1:ivar2) !!(lower z-zone)
        ubufzo(:,:,:,ivar1:ivar2)=f(:,m1:m2,n2i:n2,ivar1:ivar2) !!(upper z-zone)
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
        llbufo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n1:n1i,ivar1:ivar2)
        ulbufo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n1:n1i,ivar1:ivar2)
        uubufo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n2i:n2,ivar1:ivar2)
        lubufo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n2i:n2,ivar1:ivar2)
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
    endsubroutine initiate_isendrcv_bdry
!***********************************************************************
    subroutine finalize_isendrcv_bdry(f,ivar1_opt,ivar2_opt)
!
      use Cdata, only: bcy1,bcy2,bcz1,bcz2
!
!  Make sure the communications initiated with initiate_isendrcv_bdry are
!  finished and insert the just received boundary values.
!   Receive requests do not need to (and on OSF1 cannot) be explicitly
!  freed, since MPI_Wait takes care of this.
!
!  21-may-02/axel: communication of corners added
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j, nbufy, nbufz, nbufyz
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  1. wait until data received
!  2. set ghost zones
!  3. wait until send completed, will be overwritten in next time step
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbufy, nbufz, nbufyz
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  So far no distribution over x
!
      !(We may never do this after having seen all the trouble involved!!)
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
      use Cdata, only: bcy1,bcy2,bcz1,bcz2
!
!  Make sure the communications initiated with initiate_isendrcv_bdry are
!  finished and insert the just received boundary values.
!   Receive requests do not need to (and on OSF1 cannot) be explicitly
!  freed, since MPI_Wait takes care of this.
!
!  21-may-02/axel: communication of corners added
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j, nbufy, nbufz, nbufyz
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  1. wait until data received
!  2. set ghost zones
!  3. wait until send completed, will be overwritten in next time step
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      do i=1,nghost
        f(1:l1-1 ,i   ,:,ivar1:ivar2)=f(1:l1-1 ,m1+nghost-i,:,ivar1:ivar2)
        f(1:l1-1 ,m2+i,:,ivar1:ivar2)=f(1:l1-1 ,m2-i+1     ,:,ivar1:ivar2)
        f(l2+1:mx,i   ,:,ivar1:ivar2)=f(l2+1:mx,m1+nghost-i,:,ivar1:ivar2)
        f(l2+1:mx,m2+i,:,ivar1:ivar2)=f(l2+1:mx,m2-i+1     ,:,ivar1:ivar2)
      enddo
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
!  Receive real array from other processor.
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
!  Receive real array(:,:) from other processor.
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
!  Receive real array(:,:,:) from other processor.
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
!  Receive real array(:,:,:,:) from other processor.
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
!  Send real scalar to other processor.
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
!  Receive real array from other processor.
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
      call MPI_BCAST(ibcast_array,nbcast_array,MPI_INTEGER,ibcast_proc, &
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
    subroutine mpiallreduce_max_arr(fmax_tmp,fmax,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fmax_tmp,fmax
!
!  calculate total maximum for each array element and return to root
!
      call MPI_ALLREDUCE(fmax_tmp, fmax, nreduce, MPI_REAL, MPI_MAX, &
                      MPI_COMM_WORLD, ierr)
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
      integer :: nreduce
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
    subroutine die_gracefully()
!
!  Stop having shutdown MPI neatly
!  With at least some MPI implementations, this only stops if all
!  processors agree to call stop_it().
!
!  29-jun-05/tony: coded
!
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
      call mpifinalize
      STOP 1                    ! Return nonzero exit status
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
!
      real, dimension(ny,ny,nz) :: send_buf_y, recv_buf_y
      real, dimension(nz,ny,nz) :: send_buf_z, recv_buf_z
      real, dimension(nx,ny,nz) :: a
      real, dimension(nz) :: tmp_z
      real, dimension(ny) :: tmp_y
      integer :: i,j,sendc_y,recvc_y,sendc_z,recvc_z,px
      integer :: ytag=101,ztag=102,partner,ierr
      character :: var
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.
!  Buffers used for the z-transpose have the same size in z and x.
!
      sendc_y=ny*ny*nz; recvc_y=sendc_y
      sendc_z=nz*ny*nz; recvc_z=sendc_z
!
!  Doing x-y transpose if var='y'
!
      if (var=='y') then
        if (nxgrid/=nygrid) then
          print*,'transp: need to have nxgrid=nygrid for var==y'
          call stop_it('Inconsistency: nxgrid/=nygrid')
        endif
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
!  The following communication patterns is kind of self-regulated. It
!  avoids deadlock, because there will always be at least one matching
!  pair of processors; the others will have to wait until their partner
!  posts the corresponding request.
!    Doing send and recv together for a given pair of processors
!  (although in an order that avoids deadlocking) allows us to get away
!  with only send_buf and recv_buf as buffers
!
        do px=0,nprocy-1
          if (px/=ipy) then
            partner=px+ipz*nprocy ! = iproc + (px-ipy)
            if (ip<=6) print*,'transp: MPICOMM: ipy,ipz,px,partner=',ipy,ipz,px,partner
            send_buf_y=a(px*ny+1:(px+1)*ny,:,:)
            if (px<ipy) then      ! above diagonal: send first, receive then
              call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,ierr)
              call MPI_RECV (recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,ierr)
            elseif (px>ipy) then  ! below diagonal: receive first, send then
              call MPI_RECV (recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,ierr)
              call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,ierr)
            endif
            a(px*ny+1:(px+1)*ny,:,:)=recv_buf_y
          endif
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
          do i=1,ny-1
            do j=i+1,ny
              tmp_z=a(i+px*ny,j,:)
              a(i+px*ny,j,:)=a(j+px*ny,i,:)
              a(j+px*ny,i,:)=tmp_z
            enddo
          enddo
        enddo
!
!  Doing x-z transpose if var='z'
!
      elseif (var=='z') then
        if (nzgrid/=nxgrid) then
          if (lroot) print*, 'transp: need to have nzgrid=nxgrid for var==z'
          call stop_it('transp: inconsistency - nzgrid/=nxgrid')
        endif
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
          do i=1,nz-1
            do j=i+1,nz
              tmp_y=a(i+px*nz,:,j)
              a(i+px*nz,:,j)=a(j+px*nz,:,i)
              a(j+px*nz,:,i)=tmp_y
            enddo
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
    subroutine fold_df(df,ivar1,ivar2)
!
!  Fold first ghost zone of df into main part of df.
!
!  20-may-2006/anders: coded
!
      use Cdata
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
    subroutine transform(a1,a2,a3,b1,b2,b3)
!
!  Subroutine to do fourier transform
!  The routine overwrites the input data
!
!  03-sep-02/nils: coded
!  05-nov-02/axel: added normalization factor
!
      real,dimension(nx,ny,nz) :: a1,b1,a2,b2,a3,b3

! Doing the x field
      if (lroot .and. ip<10) print*,'transform: doing fft of x-component'
      call fft(a1,b1, nx*ny*nz, nx, nx,-1) ! x-direction
      call transp(a1,'y')
      call transp(b1,'y')
      call fft(a1,b1, nx*ny*nz, nx, nx,-1) ! y-direction
      call transp(a1,'z')
      call transp(b1,'z')
      call fft(a1,b1, nx*ny*nz, nx, nx,-1) ! z-direction

! Doing the y field
      if (lroot .and. ip<10) print*,'transform: doing fft of y-component'
      call fft(a2,b2, nx*ny*nz, nx, nx,-1) ! x-direction
      call transp(a2,'y')
      call transp(b2,'y')
      call fft(a2,b2, nx*ny*nz, nx, nx,-1) ! y-direction
      call transp(a2,'z')
      call transp(b2,'z')
      call fft(a2,b2, nx*ny*nz, nx, nx,-1) ! z-direction

! Doing the z field
      if (lroot .and. ip<10) print*,'transform: doing fft of z-component'
      call fft(a3,b3, nx*ny*nz, nx, nx,-1) ! x-direction
      call transp(a3,'y')
      call transp(b3,'y')
      call fft(a3,b3, nx*ny*nz, nx, nx,-1) ! y-direction
      call transp(a3,'z')
      call transp(b3,'z')
      call fft(a3,b3, nx*ny*nz, nx, nx,-1) ! z-direction
!
!  Normalize
!
      a1=a1/nwgrid; a2=a2/nwgrid; a3=a3/nwgrid
      b1=b1/nwgrid; b2=b2/nwgrid; b3=b3/nwgrid
      if (lroot .and. ip<10) print*,'transform: fft has finished'
!
    endsubroutine transform
!***********************************************************************
    subroutine transform_i(a_re,a_im)
!
!  Subroutine to do fourier transform
!  The routine overwrites the input data
!
!  22-oct-02/axel+tarek: adapted from transform
!
      real,dimension(nx,ny,nz) :: a_re,a_im

      if (lroot .and. ip<10) print*,'transform_i: doing three FFTs'
      call fft(a_re,a_im, nx*ny*nz, nx, nx,-1)
      call transp(a_re,'y')
      call transp(a_im,'y')
      call fft(a_re,a_im, nx*ny*nz, nx, nx,-1)
      call transp(a_re,'z')
      call transp(a_im,'z')
      call fft(a_re,a_im, nx*ny*nz, nx, nx,-1)
!
!  Normalize
!
      a_re=a_re/nwgrid
      a_im=a_im/nwgrid
      if (lroot .and. ip<10) print*,'transform_i: fft has finished'
!
    endsubroutine transform_i
!***********************************************************************
    subroutine transform_fftpack(a_re,a_im,direction)
!
!  Subroutine to do Fourier transform
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid=nzgrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  27-oct-02/axel: adapted from transform_i, for fftpack
!
      real,dimension(nx,ny,nz) :: a_re,a_im
      complex,dimension(nx) :: ax
      real,dimension(4*nx+15) :: wsavex
      integer :: m,n
      integer,optional :: direction
      logical :: lforward=.true.
!
      if (present(direction)) then
        if (direction==-1) then
          lforward=.false.
        else
          lforward=.true.
        endif
      endif
!
!  Need to initialize cfft only once, because we require nxgrid=nygrid=nzgrid
!
      call cffti(nx,wsavex)
!
      if (lforward) then
        if (lroot .and. ip<10) print*, 'transform_fftpack: doing FFTpack in x'
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nx,ax,wsavex)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
!
!  Transform y-direction.
!      
        if (nygrid/=1) then
          if (nygrid/=nxgrid) then
            if (lroot) print*, 'transform_fftpack: must have nygrid=nxgrid!'
            call stop_it('transform_fftpack')
          endif
          call transp(a_re,'y')
          call transp(a_im,'y')
!
!  The length of the array in the y-direction is nx.
!
          if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in y'
          do n=1,nz; do m=1,ny
            ax=cmplx(a_re(:,m,n),a_im(:,m,n))
            call cfftf(nx,ax,wsavex)
            a_re(:,m,n)=real(ax)
            a_im(:,m,n)=aimag(ax)
          enddo; enddo
        endif
!
!  Transform z-direction.
!      
        if (nzgrid/=1) then
          if (nzgrid/=nxgrid) then
            if (lroot) &
                print*, 'transform_fftpack: must have nzgrid=nxgrid!'
            call stop_it('transform_fftpack')
          endif
          call transp(a_re,'z')
          call transp(a_im,'z')
!
!  The length of the array in the z-direction is also nx.
!
          if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in z'
          do n=1,nz; do m=1,ny
            ax=cmplx(a_re(:,m,n),a_im(:,m,n))
            call cfftf(nx,ax,wsavex)
            a_re(:,m,n)=real(ax)
            a_im(:,m,n)=aimag(ax)
          enddo; enddo
        endif
      else
!
!  Transform from k-space to real space. Works on arrays that have been put in
!  (z,x,y)-order by the parallel Fourier transform.
!
        if (nzgrid/=1) then
!
!  Transform z-direction back.
!      
          if (nzgrid/=nxgrid) then
            if (lroot) print*, 'transform_fftpack: must have nzgrid=nxgrid!'
            call stop_it('transform_fftpack')
          endif
!
          if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in z'
          do n=1,nz; do m=1,ny
            ax=cmplx(a_re(:,m,n),a_im(:,m,n))
            call cfftb(nx,ax,wsavex)
            a_re(:,m,n)=real(ax)
            a_im(:,m,n)=aimag(ax)
          enddo; enddo
        endif
!
!  Transform y-direction back. Must transpose to go from (z,x,y) to (y,x,z).
!
        if (nygrid/=1) then
          if (nygrid/=nxgrid) then
            if (lroot) print*, 'transform_fftpack: must have nygrid=nxgrid!'
            call stop_it('transform_fftpack')
          endif
!
          if (nzgrid/=1) then
            call transp(a_re,'z')
            call transp(a_im,'z')
          endif
!
          if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in y'
          do n=1,nz; do m=1,ny
            ax=cmplx(a_re(:,m,n),a_im(:,m,n))
            call cfftb(nx,ax,wsavex)
            a_re(:,m,n)=real(ax)
            a_im(:,m,n)=aimag(ax)
          enddo; enddo
        endif
!
!  Transform x-direction back. Transpose to go from (y,x,z) to (x,y,z).
!
        if (lroot .and. ip<10) print*, 'transform_fftpack: doing FFTpack in x'
        if (nygrid==1) then
          call transp(a_re,'z')
          call transp(a_im,'z')
        else
          call transp(a_re,'y')
          call transp(a_im,'y')
        endif
        do n=1,nz; do m=1,ny
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nx,ax,wsavex)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/nwgrid
        a_im=a_im/nwgrid
      endif
!
      if (lroot .and. ip<10) print*,'transform_fftpack: fft has finished'
!
    endsubroutine transform_fftpack
!***********************************************************************
    subroutine transform_fftpack_2d(a_re,a_im,dummy)
!
!  Subroutine to do Fourier transform
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid=nzgrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  27-oct-02/axel: adapted from transform_fftpack
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      complex, dimension(nx) :: ax
      real, dimension(4*nx+15) :: wsavex
      integer :: m,n
      integer, optional :: dummy
!
!  check whether nxgrid=nygrid=nzgrid
!
      if (nxgrid/=nygrid .or. (nxgrid/=nzgrid .and. nzgrid/=1)) then
        print*,'transform_fftpack: must have nxgrid=nygrid=nzgrid!'
        call stop_it("must have nxgrid=nygrid=nzgrid!")
      endif
!
!  need to initialize cfft only once, because nxgrid=nygrid=nzgrid
!
      call cffti(nx,wsavex)
!
      if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in x'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        call cfftf(nx,ax,wsavex)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
      call transp(a_re,'z')
      call transp(a_im,'z')
!
!  The length of the array in the z-direction is also nx
!
      if (lroot .and. ip<10) print*,'transform_fftpack: doing FFTpack in z'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        call cfftf(nx,ax,wsavex)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
!
!  Normalize
!
      a_re=a_re/nwgrid
      a_im=a_im/nwgrid
      if (lroot .and. ip<10) print*,'transform_fftpack: fft has finished'
!
      if (NO_WARN) print*,dummy  !(keep compiler quiet)
!
    endsubroutine transform_fftpack_2d
!***********************************************************************
    subroutine transform_fftpack_1d(a_re,a_im)
!
!  Subroutine to do Fourier transform
!  The routine overwrites the input data.
!  This version works currently only when nxgrid=nygrid=nzgrid!
!  The length of the work arrays for ffts is therefore always nx.
!
!  06-feb-03/nils: adapted from transform_fftpack
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      complex, dimension(nx) :: ax
      real, dimension(4*nx+15) :: wsavex
      integer :: m,n
!
!  check whether nxgrid=nygrid=nzgrid
!
      if ( (nygrid/=1.and.nygrid/=nxgrid) .or. &
           (nzgrid/=1.and.nzgrid/=nxgrid) ) then
        if (lroot) &
            print*, 'transform_fftpack_1d: must have nxgrid=nygrid=nzgrid!'
        call stop_it("transform_fftpack_1d: must have nxgrid=nygrid=nzgrid!")
      endif
!
!  need to initialize cfft only once, because nxgrid=nygrid=nzgrid
!
      call cffti(nx,wsavex)
!
      if (lroot .and. ip<10) print*,'transform_fftpack_1d: doing FFTpack in x'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        call cfftf(nx,ax,wsavex)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
!
!  Normalize
!
      a_re=a_re/nxgrid
      a_im=a_im/nxgrid
      if (lroot .and. ip<10) print*,'transform_fftpack_1d: fft has finished'
!
    endsubroutine transform_fftpack_1d
!***********************************************************************
    subroutine transform_fftpack_shear(a_re,a_im,direction)
!
!  Subroutine to do Fourier transform in shearing coordinates.
!  The routine overwrites the input data
!
!  The transform from real space to Fourier space goes as follows:
!
!  (x,y,z)
!     |
!  (y,x,z) - (ky,x,z) - (ky',x,z)
!                            |
!                       (x,ky',z) - (kx,ky',z)
!                                        |
!                                   (z,ky',kx) - (kz,ky',kx)
!
!  Vertical lines refer to a transposition operation, horizontal lines to any
!  other operation. Here ky' refers to a coordinate frame where y has been
!  transformed so that the x-direction is purely periodic (not shear periodic).
!  The transformation from Fourier space to real space takes place like sketched
!  in the diagram, only in the opposite direction (starting lower right).
!
!  For 2-D runs two cases of degeneracy have to be taken into account:
!  (- refers to a missing direction).
!
!  (x,y,-)
!     |
!  (y,x,-) - (ky,x,-) - (ky',x,-)
!                            |
!                       (x, ky',-) - (kx,ky',-)
!
!  (x,-,z) - (kx,-,z)
!                |
!            (z,-,kx) - (kz,-,kx)
!
!  25-may-06/anders: adapted from transform_fftpack
!
      use Cdata, only: pi, Lx, x0, x, deltay, ky_fft
!
      real, dimension (nx,ny,nz) :: a_re,a_im
      integer, optional :: direction
!
      complex, dimension (nxgrid) :: ax
      complex, dimension (nxgrid) :: ay
      complex, dimension (nxgrid) :: az
      real, dimension (4*nxgrid+15) :: wsave
      real :: deltay_x
      integer :: l,m,n
      logical :: lforward=.true.
!
!  if nxgrid/=nygrid/=nzgrid, stop.

      if (nygrid/=nxgrid .and. nygrid /= 1) then
        print*,'transpose_fftpack_shear: need to have nygrid=nxgrid if nygrid /=1.'
        call stop_it('Inconsistency: nygrid/=nxgrid')
      endif
      if (nzgrid/=nxgrid .and. nzgrid /= 1) then
        print*,'transpose_fftpack_shear: need to have nzgrid=nxgrid if nzgrid /=1.'
        call stop_it('Inconsistency: nzgrid/=nxgrid')
      endif
!
      if (present(direction)) then
        if (direction==-1) then
          lforward=.false.
        else
          lforward=.true.
        endif
      endif
!
!  Need to initialize cfft only once, because we require nxgrid=nygrid=nzgrid.
!
      call cffti(nxgrid,wsave)
!
      if (lforward) then
!
!  Transform y-direction. Must start with y, because x is not periodic (yet).
!
        if (nygrid/=1) then
          if (lroot.and.ip<10) &
              print*, 'doing FFTpack in y, direction=',direction
          call transp(a_re,'y')
          call transp(a_im,'y')
          do n=1,nz; do l=1,ny
            ay=cmplx(a_re(:,l,n),a_im(:,l,n))
            call cfftf(nxgrid,ay,wsave)
!  Shift y-coordinate so that x-direction is periodic. This is best done in
!  k-space, by multiplying the complex amplitude by exp[i*ky*deltay(x)].
            deltay_x=-deltay*(x(l+nghost+ipy*ny)-(x0+Lx/2))/Lx
            ay(2:nxgrid)=ay(2:nxgrid)*exp(cmplx(0.0, ky_fft(2:nxgrid)*deltay_x))
            a_re(:,l,n)=real(ay)
            a_im(:,l,n)=aimag(ay)
          enddo; enddo
        endif
!
!  Transform x-direction.
!
        if (lroot.and.ip<10) print*, 'doing FFTpack in x, direction=',direction
        if (nygrid/=1) then
          call transp(a_re,'y')
          call transp(a_im,'y')
        endif
        do m=1,ny; do n=1,nz
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftf(nxgrid,ax,wsave)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
!
!  Transform z-direction.
!
        if (nzgrid/=1) then
          if (lroot.and.ip<10) &
              print*, 'doing FFTpack in z, direction=',direction
          call transp(a_re,'z')
          call transp(a_im,'z')
          do l=1,nz; do m=1,ny
            az=cmplx(a_re(:,m,l),a_im(:,m,l))
            call cfftf(nxgrid,az,wsave)
            a_re(:,m,l)=real(az)
            a_im(:,m,l)=aimag(az)
          enddo; enddo
        endif
      else
!
!  Transform z-direction back.
!
        if (nzgrid/=1) then
          if (lroot.and.ip<10) &
              print*, 'doing FFTpack in z, direction=',direction
          do l=1,nz; do m=1,ny
            az=cmplx(a_re(:,m,l),a_im(:,m,l))
            call cfftb(nxgrid,az,wsave)
            a_re(:,m,l)=real(az)
            a_im(:,m,l)=aimag(az)
          enddo; enddo
        endif
!
!  Transform x-direction back.
!
        if (lroot.and.ip<10) print*, 'doing FFTpack in x, direction=',direction
        if (nzgrid/=1) then
          call transp(a_re,'z')
          call transp(a_im,'z')
        endif
        do m=1,ny; do n=1,nz
          ax=cmplx(a_re(:,m,n),a_im(:,m,n))
          call cfftb(nxgrid,ax,wsave)
          a_re(:,m,n)=real(ax)
          a_im(:,m,n)=aimag(ax)
        enddo; enddo
!
!  Transform y-direction back. Would be more practical to end with the
!  x-direction, but the transformation from y' to y means that we must transform
!  the y-direction last.
!
        if (nygrid/=1) then
          call transp(a_re,'y')
          call transp(a_im,'y')
          if (lroot.and.ip<10) &
              print*, 'doing FFTpack in y, direction=',direction
          do l=1,ny; do n=1,nz
            ay=cmplx(a_re(:,l,n),a_im(:,l,n))
!  Shift y-coordinate back to regular frame (see above).
            deltay_x=-deltay*(x(l+nghost+ipy*ny)-(x0+Lx/2))/Lx
            ay(2:nxgrid)=ay(2:nxgrid)*exp(cmplx(0.0,-ky_fft(2:nxgrid)*deltay_x))
            call cfftb(nxgrid,ay,wsave)
            a_re(:,l,n)=real(ay)
            a_im(:,l,n)=aimag(ay)
          enddo; enddo
          call transp(a_re,'y')  ! Deliver array back in (x,y,z) order.
          call transp(a_im,'y')
        endif
      endif
!
!  Normalize
!
      if (lforward) then
        a_re=a_re/nwgrid
        a_im=a_im/nwgrid
      endif
!
    endsubroutine transform_fftpack_shear
!***********************************************************************
    subroutine transform_nr(a_re,a_im)
!
!  Subroutine to do Fourier transform using Numerical Recipes routine.
!  Note that this routine requires that nx, ny, and nz are powers of 2.
!  The routine overwrites the input data.
!
!  30-oct-02/axel: adapted from transform_fftpack for Numerical Recipes
!
      real, dimension(nx,ny,nz) :: a_re,a_im
      complex, dimension(nx) :: ax
      integer :: m,n
!
!  This Fourier transform would work, but it's very slow!
!  Even the compilation is very slow, so we better get rid of it!
!
      call stop_it("transform_nr: currently disabled!")
!
      if (lroot .and. ip<10) print*,'transform_nr: doing FFT_nr in x'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        !call four1(ax,nx,-1)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
      call transp(a_re,'y')
      call transp(a_im,'y')
!
!  The length of the array in the y-direction is nx
!  (remember: nxgrid=nygrid=nzgrid!)
!
      if (lroot .and. ip<10) print*,'transform_nr: doing FFT_nr in y'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        !call four1(ax,nx,-1)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
      call transp(a_re,'z')
      call transp(a_im,'z')
!
!  The length of the array in the z-direction is also nx
!
      if (lroot .and. ip<10) print*,'transform_nr: doing FFT_nr in z'
      do n=1,nz; do m=1,ny
        ax=cmplx(a_re(:,m,n),a_im(:,m,n))
        !call four1(ax,nx,-1)
        a_re(:,m,n)=real(ax)
        a_im(:,m,n)=aimag(ax)
      enddo; enddo
!
!  Normalize
!
      a_re=a_re/nwgrid
      a_im=a_im/nwgrid
      if (lroot .and. ip<10) print*,'transform_nr: fft has finished'
!
    endsubroutine transform_nr
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
!
      real, dimension (ny,nz) :: a_re
      real :: shift_y
!
      complex, dimension(nygrid) :: ay
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
            call stop_it('fourier_shift_yz')
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
!  Transform to Fourier space.
!
      call cffti(nygrid,wsavey)
      do n=1,nz_new
        ay=cmplx(a_re_new(:,n),0.0)
        call cfftf(nygrid,ay,wsavey)
!  Do the shifting in Fourier space.
        ay(2:nygrid)=ay(2:nygrid)*exp(cmplx(0.0,-ky_fft(2:nygrid)*shift_y))
!
        a_re_new(:,n)=real(ay)
        a_im_new(:,n)=aimag(ay)
      enddo
!
!  Back to real space.
!
      call cffti(nygrid,wsavey)
      do n=1,nz_new
        ay=cmplx(a_re_new(:,n),a_im_new(:,n))
        call cfftb(nygrid,ay,wsavey)
        a_re_new(:,n)=real(ay)/nygrid
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
endmodule Mpicomm
