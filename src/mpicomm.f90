!!!!!!!!!!!!!!!!!!!!!
!!!  mpicomm.f90  !!!
!!!!!!!!!!!!!!!!!!!!!

!!! Module with MPI stuff

! NB: This was previously called mpicommyz.f90 and distributes in two
! directions.

module Mpicomm

  use Cparam
  use Cdata, only: lroot

  implicit none

  interface mpibcast_real		! Overload the `mpibcast_real' function
    module procedure mpibcast_real_scl
    module procedure mpibcast_real_arr
  endinterface

  include 'mpif.h'
 
  integer, parameter :: nbufy=nx*nz*nghost*mvar
  integer, parameter :: nbufz=nx*ny*nghost*mvar

  real, dimension (nx,nghost,nz,mvar) :: lbufyi,ubufyi,lbufyo,ubufyo
  real, dimension (nx,ny,nghost,mvar) :: lbufzi,ubufzi,lbufzo,ubufzo
  integer, dimension (ny*nz) :: mm,nn
  integer :: ierr,imn
  integer :: nprocs,iproc,root=0
  integer :: ipy,ipz
  integer :: tolowy=3,touppy=4,tolowz=5,touppz=6 ! msg. tags
  integer :: isend_rq_tolowy,isend_rq_touppy,irecv_rq_fromlowy,irecv_rq_fromuppy
  integer :: isend_rq_tolowz,isend_rq_touppz,irecv_rq_fromlowz,irecv_rq_fromuppz
  integer, dimension(MPI_STATUS_SIZE) :: isend_stat_tl,isend_stat_tu
  integer, dimension(MPI_STATUS_SIZE) :: irecv_stat_fl,irecv_stat_fu
  integer :: ylneigh,zlneigh ! `lower' neighbours
  integer :: yuneigh,zuneigh ! `upper' neighbours
!! Moved to Cdata to save lots of `use Mpicomm':
!  logical :: lroot              ! is this the root process?
  logical, dimension (ny*nz) :: necessary=.false.
  character directory*12

  contains

!***********************************************************************
    subroutine mpicomm_init()
!
!  Initialise MPI communication and set up some variables
!  The arrays leftneigh and rghtneigh give the processor numbers
!  to the left and to the right.
!
!  20-aug-01/wolf: coded
!  31-aug-01/axel: added to 3-D
!  15-sep-01/axel: adapted from Wolfgang's version
!
      use General
!
      integer :: i,m,n
      character chproc*4
!
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, iproc , ierr)
!
!  consistency checks
!
      if (nprocs /= nprocy*nprocz) then
        if(iproc == root) then
          print*, 'Compiled with NCPUS = ', ncpus, &
          ', but running on ', nprocs, ' processors'
        endif
        STOP 'Inconsistency 1'
      endif
!
      if ((nprocy*ny /= nygrid) .or. &
          (nprocz*nz /= nzgrid)) then
        if(lroot) then
          write(0,'(A,2I4,A,2I4,A)') &
               'nproc[y-z]*n[y-z] = (', &
               nprocy*ny, nprocz*nz, &
               ') /= n[yz]grid= (', &
               nygrid, nzgrid, ")"
        endif
        call stop_all('Inconsistency 2')
      endif
!
!  position on the processor grid
!  x is fastest direction, z slowest
!
      ipy = modulo(iproc, nprocy)
      ipz = iproc/(nprocy)
!
!  set up `lower' and `upper' neighbours
!
      ylneigh = (ipz*nprocy+modulo(ipy-1,nprocy))
      yuneigh = (ipz*nprocy+modulo(ipy+1,nprocy))
      zlneigh = (modulo(ipz-1,nprocz)*nprocy+ipy)
      zuneigh = (modulo(ipz+1,nprocz)*nprocy+ipy)
!
      if (ip<6) then
        write(0,'(I3,": ",A,3(3I3,"  "))') &
             iproc, 'ipy,ipz, {y,z}{l,u}neigh = ', &
             ipy,ipz, &
             ylneigh, zlneigh, &
             yuneigh, zuneigh
      endif
!
!  produce index-array the the sequence of points to be worked through first
!  inner box
!
      imn=1
      do n=n1i+1,n2i-1
        do m=m1i+1,m2i-1
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
      necessary(imn)=.true.
!
!  do the lower stripe in the n-direction
!
      do n=n2i,n2
        do m=m1i+1,m2i-1
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
!
!  upper stripe in the n-direction
!
      do n=n1,n1i
        do m=m1i+1,m2i-1
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
!
!  left and right hand boxes
!
      do n=n1,n2
        do m=m1,m1i
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
        do m=m2i,m2
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
!
!  produce a directory name for the output
!
      call chn(iproc,chproc)
      directory='tmp/proc'//chproc
!
    endsubroutine mpicomm_init
!***********************************************************************
    subroutine initiate_isendrcv_bdry(f)
!
!  Isend and Irecv boundary values. Called in the beginning of pde.
!  Does not wait for the receives to finish (done in finalise_isendrcv_bdry)
!  leftneigh and rghtneigh are initialized by mpicomm_init
!
      real, dimension (mx,my,mz,mvar) :: f
!
!  Periodic boundary conditions in x
!
      f( 1:l1-1,m1:m2,n1:n2,:) = f(l2i:l2,m1:m2,n1:n2,:)
      f(l2+1:mx,m1:m2,n1:n2,:) = f(l1:l1i,m1:m2,n1:n2,:)
!
!  Periodic boundary conditions in y
!
      if (nprocy==0) then
        f(l1:l2, 1:m1-1,n1:n2,:) = f(l1:l2,m2i:m2,n1:n2,:)
        f(l1:l2,m2+1:my,n1:n2,:) = f(l1:l2,m1:m1i,n1:n2,:)
      else
        lbufyo=f(l1:l2,m1:m1i,n1:n2,:)
        ubufyo=f(l1:l2,m2i:m2,n1:n2,:)
        call MPI_ISEND(lbufyo,nbufy,MPI_REAL,ylneigh,tolowy,MPI_COMM_WORLD,isend_rq_tolowy,ierr)
        call MPI_ISEND(ubufyo,nbufy,MPI_REAL,yuneigh,touppy,MPI_COMM_WORLD,isend_rq_touppy,ierr)
        call MPI_IRECV(ubufyi,nbufy,MPI_REAL,yuneigh,tolowy,MPI_COMM_WORLD,irecv_rq_fromuppy,ierr)
        call MPI_IRECV(lbufyi,nbufy,MPI_REAL,ylneigh,touppy,MPI_COMM_WORLD,irecv_rq_fromlowy,ierr)
      endif
!
!  Periodic boundary conditions in z
!
      if (nprocz==0) then
        f(l1:l2,m1:m2, 1:n1-1,:) = f(l1:l2,m1:m2,n2i:n2,:)
        f(l1:l2,m1:m2,n2+1:mz,:) = f(l1:l2,m1:m2,n1:n1i,:)
      else
        lbufzo=f(l1:l2,m1:m2,n1:n1i,:)
        ubufzo=f(l1:l2,m1:m2,n2i:n2,:)
        call MPI_ISEND(lbufzo,nbufz,MPI_REAL,zlneigh,tolowz,MPI_COMM_WORLD,isend_rq_tolowz,ierr)
        call MPI_ISEND(ubufzo,nbufz,MPI_REAL,zuneigh,touppz,MPI_COMM_WORLD,isend_rq_touppz,ierr)
        call MPI_IRECV(ubufzi,nbufz,MPI_REAL,zuneigh,tolowz,MPI_COMM_WORLD,irecv_rq_fromuppz,ierr)
        call MPI_IRECV(lbufzi,nbufz,MPI_REAL,zlneigh,touppz,MPI_COMM_WORLD,irecv_rq_fromlowz,ierr)
      endif
!
    endsubroutine initiate_isendrcv_bdry
!***********************************************************************
    subroutine finalise_isendrcv_bdry(f)
!
!  Make sure the communications initiated with initiate_isendrcv_bdry are
!  finished and insert the just received boundary values.
!    Receive requests do not need to (and on OSF1 cannot) be explicitly
!  freed, since MPI_Wait takes care of this.
!
      real, dimension (mx,my,mz,mvar) :: f
!
!  Periodic boundary conditions in z, combined with communication.
!
      call MPI_WAIT(irecv_rq_fromuppy,irecv_stat_fu,ierr)
      call MPI_WAIT(irecv_rq_fromlowy,irecv_stat_fl,ierr)
      call MPI_WAIT(irecv_rq_fromuppz,irecv_stat_fu,ierr)
      call MPI_WAIT(irecv_rq_fromlowz,irecv_stat_fl,ierr)
!
!  set boundary zones (this cannot be avoided)
!
      f(l1:l2, 1:m1-1,n1:n2,:)=lbufyi
      f(l1:l2,m2+1:my,n1:n2,:)=ubufyi
!
      f(l1:l2,m1:m2, 1:n1-1,:)=lbufzi
      f(l1:l2,m1:m2,n2+1:mz,:)=ubufzi
!
!  need to wait until send is completed before buffer can be overwritten
!
      call MPI_WAIT(isend_rq_tolowy,isend_stat_tl,ierr)
      call MPI_WAIT(isend_rq_touppy,isend_stat_tu,ierr)
      call MPI_WAIT(isend_rq_tolowz,isend_stat_tl,ierr)
      call MPI_WAIT(isend_rq_touppz,isend_stat_tu,ierr)
!
!  make sure the other precessors don't carry on sending new data
!  which could be mistaken for an earlier time
!
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    endsubroutine finalise_isendrcv_bdry
!***********************************************************************
    subroutine mpibcast_int(ibcast_array,nbcast_array)
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: ibcast_array
!
      call MPI_BCAST(ibcast_array,nbcast_array,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    endsubroutine mpibcast_int
!***********************************************************************
    subroutine mpibcast_real_scl(bcast_array,nbcast_array)
!
      integer :: nbcast_array
      real :: bcast_array
!
!  is being used when nbcast_array=1 (eg when dt is being communicated)
!
      call MPI_BCAST(bcast_array,nbcast_array,MPI_REAL,root,MPI_COMM_WORLD,ierr)
    endsubroutine mpibcast_real_scl
!***********************************************************************
    subroutine mpibcast_real_arr(bcast_array,nbcast_array)
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array
!
!  this works for the general case when nbcast_array is not 1
!
      call MPI_BCAST(bcast_array,nbcast_array,MPI_REAL,root,MPI_COMM_WORLD,ierr)
    endsubroutine mpibcast_real_arr
!***********************************************************************
    subroutine mpireduce_max(fmax_tmp,fmax,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fmax_tmp,fmax
!
!  calculate total maximum for each array element and return to root
!
      call MPI_REDUCE(fmax_tmp, fmax, nreduce, MPI_REAL, MPI_MAX, root, &
                      MPI_COMM_WORLD, ierr)
    endsubroutine mpireduce_max
!***********************************************************************
    subroutine mpireduce_sum(fsum_tmp,fsum,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fsum_tmp,fsum
!
!  calculate total sum for each array element and return to root
!
      call MPI_REDUCE(fsum_tmp, fsum, nreduce, MPI_REAL, MPI_SUM, root, &
                      MPI_COMM_WORLD, ierr)
    endsubroutine mpireduce_sum
!***********************************************************************
    subroutine mpifinalize
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_FINALIZE(ierr)
    endsubroutine mpifinalize
!***********************************************************************

endmodule Mpicomm
