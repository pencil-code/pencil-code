! $Id: mpicomm.f90,v 1.68 2002-11-19 20:53:33 dobler Exp $

!!!!!!!!!!!!!!!!!!!!!
!!!  mpicomm.f90  !!!
!!!!!!!!!!!!!!!!!!!!!

!!! Module with MPI stuff

!  Data layout for each processor (`-' marks ghost points, `+' real
!  points of the processor shown)
!
!         n = 1     - - - - - - - . - - - - - - -
!             2     - - - - - - - . - - - - - - -
!             3     - - - - - - - . - - - - - - -
!             . n1  - - - + + + + . + + + + - - -
!             . .   - - - + + + + . + + + + - - -
!             . .   - - - + + + + . + + + + - - -
!             . .   - - - + + + + . + + + + - - -
!                   . . . . . . . . . . . . . . . .
!             . .   - - - + + + + . + + + + - - -
!             . .   - - - + + + + . + + + + - - -
!             . .   - - - + + + + . + + + + - - -
!             . n2  - - - + + + + . + + + + - - -
!             .     - - - - - - - . - - - - - - -
!             .     - - - - - - - . - - - - - - -
!             mz    - - - - - - - . - - - - - - -
!                        
!                         m1. . .   . . . m2
!               m = 1 2 3 . . . .   . . . . . . my
!                           

module Mpicomm

  use Cparam
  use Cdata, only: iproc,ipx,ipy,ipz,root,lroot

  implicit none

  interface mpibcast_real               ! Overload the `mpibcast_real' function
    module procedure mpibcast_real_scl
    module procedure mpibcast_real_arr
  endinterface

  interface mpibcast_logical            ! Overload
    module procedure mpibcast_logical_scl
  endinterface

  include 'mpif.h'

  integer, parameter :: nbufx_gh=my*mz*nghost*mvar ! For shear
  integer, parameter :: nbufy=mx*nz*nghost*mvar
  integer, parameter :: nbufz=mx*ny*nghost*mvar
  integer, parameter :: nbufyz=mx*nghost*nghost*mvar

  real, dimension (mx,nghost,nz,mvar) :: lbufyi,ubufyi,lbufyo,ubufyo
  real, dimension (mx,ny,nghost,mvar) :: lbufzi,ubufzi,lbufzo,ubufzo
  real, dimension (mx,nghost,nghost,mvar) :: llbufi,lubufi,uubufi,ulbufi
  real, dimension (mx,nghost,nghost,mvar) :: llbufo,lubufo,uubufo,ulbufo
  real, dimension (nghost,my,mz,mvar) :: fahi,falo,fbhi,fblo,fao,fbo ! For shear
  integer :: nextya, nextyb, lastya, lastyb, displs ! For shear
  integer, dimension (ny*nz) :: mm,nn
  integer :: ierr,imn
  integer :: nprocs
  integer :: tolowy=3,touppy=4,tolowz=5,touppz=6 ! msg. tags
  integer :: TOll=7,TOul=8,TOuu=9,TOlu=10 ! msg. tags for corners
  integer :: io_perm=20,io_succ=21
  integer :: isend_rq_tolowy,isend_rq_touppy,irecv_rq_fromlowy,irecv_rq_fromuppy
  integer :: isend_rq_tolowz,isend_rq_touppz,irecv_rq_fromlowz,irecv_rq_fromuppz
  integer :: isend_rq_TOll,isend_rq_TOul,isend_rq_TOuu,isend_rq_TOlu  !(corners)
  integer :: irecv_rq_FRuu,irecv_rq_FRlu,irecv_rq_FRll,irecv_rq_FRul  !(corners)
  integer :: isend_rq_tolastya,isend_rq_tonextya, &
             irecv_rq_fromlastya,irecv_rq_fromnextya ! For shear
  integer :: isend_rq_tolastyb,isend_rq_tonextyb, &
             irecv_rq_fromlastyb,irecv_rq_fromnextyb ! For shear
  integer, dimension(MPI_STATUS_SIZE) :: isend_stat_tl,isend_stat_tu
  integer, dimension(MPI_STATUS_SIZE) :: irecv_stat_fl,irecv_stat_fu
  integer, dimension(MPI_STATUS_SIZE) :: isend_stat_Tll,isend_stat_Tul, &
                                         isend_stat_Tuu,isend_stat_Tlu
  integer, dimension(MPI_STATUS_SIZE) :: irecv_stat_Fuu,irecv_stat_Flu, &
                                         irecv_stat_Fll,irecv_stat_Ful
  integer :: ylneigh,zlneigh ! `lower' neighbours
  integer :: yuneigh,zuneigh ! `upper' neighbours
  integer :: llcorn,lucorn,uucorn,ulcorn !!(the 4 corners in yz-plane)
  logical, dimension (ny*nz) :: necessary=.false.

  contains

!***********************************************************************
    subroutine mpicomm_init()
!
!  Initialise MPI communication and set up some variables
!  The arrays leftneigh and rghtneigh give the processor numbers
!  to the left and to the right.
!
!  Before the communication has been completed, the nghost=3 layers next
!  to the processor boundary (m1, m2, n1, or n2) cannot be used yet.
!  In the mean time we can calculate the interior points sufficiently far
!  away from the boundary points. Here we calculate the order in which
!  m and n are executed. At one point, necessary(imn)=.true., which is
!  the moment when all communication must be completed.
!
!  20-aug-01/wolf: coded
!  31-aug-01/axel: added to 3-D
!  15-sep-01/axel: adapted from Wolfgang's version
!  21-may-02/axel: communication of corners added
!   6-jun-02/axel: generalized to allow for ny=1
!
      use General
      use Cdata, only: lmpicomm
!
      integer :: i,m,n
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
        if(lroot) then
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
      ipx = 0
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
           'MPICOMM neighbors ', &
           iproc,ipy,ipz, &
           ylneigh,llcorn,zlneigh,ulcorn,yuneigh,uucorn,zuneigh,lucorn
!
!  produce index-array for the sequence of points to be worked through:
!  first inner box (while communication of ghost values takes place),
!  then boundary zones.
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
      do n=max(n2i,n1+1),n2
        do m=m1i+1,m2i-1
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
!
!  upper stripe in the n-direction
!
      do n=n1,min(n1i,n2)
        do m=m1i+1,m2i-1
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
!
!  left and right hand boxes
!  NOTE: need to have min(m1i,m2) instead of just m1i, and max(m2i,m1)
!  instead of just m2i, to make sure the case ny=1 works ok, and
!  also that the same m is not set in both loops.
!
      do n=n1,n2
        do m=m1,min(m1i,m2)
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
        do m=max(m2i,m1+1),m2
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
!
    endsubroutine mpicomm_init
!***********************************************************************
    subroutine initiate_isendrcv_bdry(f)
!
!  Isend and Irecv boundary values. Called in the beginning of pde.
!  Does not wait for the receives to finish (done in finalise_isendrcv_bdry)
!  leftneigh and rghtneigh are initialized by mpicomm_init
!
!  21-may-02/axel: communication of corners added
!
      real, dimension (mx,my,mz,mvar) :: f
!
!  So far no distribution over x
!
      !(We may never do this after having seen all the trouble involved!!)
!
!  Periodic boundary conditions in y
!
      if (nprocy>1) then
        lbufyo=f(:,m1:m1i,n1:n2,:)  !!(lower y-zone)
        ubufyo=f(:,m2i:m2,n1:n2,:)  !!(upper y-zone)
        call MPI_IRECV(ubufyi,nbufy,MPI_REAL,yuneigh,tolowy,MPI_COMM_WORLD,irecv_rq_fromuppy,ierr)
        call MPI_IRECV(lbufyi,nbufy,MPI_REAL,ylneigh,touppy,MPI_COMM_WORLD,irecv_rq_fromlowy,ierr)
        call MPI_ISEND(lbufyo,nbufy,MPI_REAL,ylneigh,tolowy,MPI_COMM_WORLD,isend_rq_tolowy,ierr)
        call MPI_ISEND(ubufyo,nbufy,MPI_REAL,yuneigh,touppy,MPI_COMM_WORLD,isend_rq_touppy,ierr)
      endif
!
!  Periodic boundary conditions in z
!
      if (nprocz>1) then
        lbufzo=f(:,m1:m2,n1:n1i,:)  !!(lower z-zone)
        ubufzo=f(:,m1:m2,n2i:n2,:)  !!(upper z-zone)
        call MPI_IRECV(ubufzi,nbufz,MPI_REAL,zuneigh,tolowz,MPI_COMM_WORLD,irecv_rq_fromuppz,ierr)
        call MPI_IRECV(lbufzi,nbufz,MPI_REAL,zlneigh,touppz,MPI_COMM_WORLD,irecv_rq_fromlowz,ierr)
        call MPI_ISEND(lbufzo,nbufz,MPI_REAL,zlneigh,tolowz,MPI_COMM_WORLD,isend_rq_tolowz,ierr)
        call MPI_ISEND(ubufzo,nbufz,MPI_REAL,zuneigh,touppz,MPI_COMM_WORLD,isend_rq_touppz,ierr)
      endif
!
!  The four corners (in counter-clockwise order)
!
      if (nprocy>1.and.nprocz>1) then
        llbufo=f(:,m1:m1i,n1:n1i,:)
        ulbufo=f(:,m2i:m2,n1:n1i,:)
        uubufo=f(:,m2i:m2,n2i:n2,:)
        lubufo=f(:,m1:m1i,n2i:n2,:)
        call MPI_IRECV(uubufi,nbufyz,MPI_REAL,uucorn,TOll,MPI_COMM_WORLD,irecv_rq_FRuu,ierr)
        call MPI_IRECV(lubufi,nbufyz,MPI_REAL,lucorn,TOul,MPI_COMM_WORLD,irecv_rq_FRlu,ierr)
        call MPI_IRECV(llbufi,nbufyz,MPI_REAL,llcorn,TOuu,MPI_COMM_WORLD,irecv_rq_FRll,ierr)
        call MPI_IRECV(ulbufi,nbufyz,MPI_REAL,ulcorn,TOlu,MPI_COMM_WORLD,irecv_rq_FRul,ierr)
        call MPI_ISEND(llbufo,nbufyz,MPI_REAL,llcorn,TOll,MPI_COMM_WORLD,isend_rq_TOll,ierr)
        call MPI_ISEND(ulbufo,nbufyz,MPI_REAL,ulcorn,TOul,MPI_COMM_WORLD,isend_rq_TOul,ierr)
        call MPI_ISEND(uubufo,nbufyz,MPI_REAL,uucorn,TOuu,MPI_COMM_WORLD,isend_rq_TOuu,ierr)
        call MPI_ISEND(lubufo,nbufyz,MPI_REAL,lucorn,TOlu,MPI_COMM_WORLD,isend_rq_TOlu,ierr)
      endif
!
!  communication sample
!
      if (ip<7.and.ipy==0.and.ipz==3) &
        print*,'MPICOMM send lu: ',iproc,lubufo(nx/2+4,:,1,2),' to ',lucorn
!
    endsubroutine initiate_isendrcv_bdry
!***********************************************************************
    subroutine finalise_isendrcv_bdry(f)
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
      real, dimension (mx,my,mz,mvar) :: f
      character (len=160) :: errmesg
      integer :: j
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
        do j=1,mvar
           if (ipy /= 0 .OR. bcy1(j)=='p') then
              f(:, 1:m1-1,n1:n2,j)=lbufyi(:,:,:,j)  !!(set lower buffer)
           endif
           if (ipy /= nprocy-1 .OR. bcy2(j)=='p') then 
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
        do j=1,mvar
           if (ipz /= 0 .OR. bcz1(j)=='p') then 
              f(:,m1:m2, 1:n1-1,j)=lbufzi(:,:,:,j)  !!(set lower buffer)
           endif
           if (ipz /= nprocz-1 .OR. bcz2(j)=='p') then 
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
        do j=1,mvar
           if (ipz /= 0 .OR. bcz1(j)=='p') then 
              if (ipy /= 0 .OR. bcy1(j)=='p') then 
                 f(:, 1:m1-1, 1:n1-1,j)=llbufi(:,:,:,j)  !!(set ll corner)
              endif
              if (ipy /= nprocy-1 .OR. bcy2(j)=='p') then
                 f(:,m2+1:my, 1:n1-1,j)=ulbufi(:,:,:,j)  !!(set ul corner)
              endif
           endif
           if (ipz /= nprocz-1 .OR. bcz2(j)=='p') then 
              if (ipy /= nprocy-1 .OR. bcy2(j)=='p') then 
                 f(:,m2+1:my,n2+1:mz,j)=uubufi(:,:,:,j)  !!(set uu corner)
              endif
              if (ipy /= 0 .OR. bcy1(j)=='p') then
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
!
      if (ip<7.and.ipy==3.and.ipz==0) &
        print*,'MPICOMM recv ul: ',iproc,ulbufi(nx/2+4,:,1,2),' from ',ulcorn
!
!  make sure the other precessors don't carry on sending new data
!  which could be mistaken for an earlier time
!
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    endsubroutine finalise_isendrcv_bdry
!***********************************************************************
   subroutine initiate_shearing(f)
!
!  Subroutine for shearing sheet boundary conditions
!
!  20-june-02/nils: adapted from pencil_mpi
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      double precision :: deltay_dy, frak, c1, c2, c3, c4, c5, c6
      integer :: i, ystep
      integer :: tolastya=11, tolastyb=12, tonextya=13, tonextyb=14
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
         f(1:l1-1,m1:m2,:,:)=c1*cshift(f(l2i:l2,m1:m2,:,:),-displs+2,2) &
              +c2*cshift(f(l2i:l2,m1:m2,:,:),-displs+1,2) &
              +c3*cshift(f(l2i:l2,m1:m2,:,:),-displs,2) &
              +c4*cshift(f(l2i:l2,m1:m2,:,:),-displs-1,2) &
              +c5*cshift(f(l2i:l2,m1:m2,:,:),-displs-2,2) &
              +c6*cshift(f(l2i:l2,m1:m2,:,:),-displs-3,2)
         f(l2+1:mx,m1:m2,:,:)=c1*cshift(f(l1:l1i,m1:m2,:,:),displs-2,2) &
              +c2*cshift(f(l1:l1i,m1:m2,:,:),displs-1,2) &
              +c3*cshift(f(l1:l1i,m1:m2,:,:),displs,2) &
              +c4*cshift(f(l1:l1i,m1:m2,:,:),displs+1,2) &
              +c5*cshift(f(l1:l1i,m1:m2,:,:),displs+2,2) &
              +c6*cshift(f(l1:l1i,m1:m2,:,:),displs+3,2)
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
         fao = f(l1:l1i,:,:,:)
         fbo = f(l2i:l2,:,:,:)
         if (lastya/=iproc) then
            call MPI_ISEND(fao,nbufx_gh,MPI_REAL,lastya,tonextyb,MPI_COMM_WORLD,isend_rq_tolastya,ierr)
         endif
         if (nextyb==iproc) then
            fbhi=fao
         else
            call MPI_IRECV(fbhi,nbufx_gh,MPI_REAL,nextyb,tonextyb,MPI_COMM_WORLD,irecv_rq_fromnextyb,ierr)
         endif
         if (nextya/=iproc) then
            call MPI_ISEND(fao,nbufx_gh,MPI_REAL,nextya,tolastyb,MPI_COMM_WORLD,isend_rq_tonextya,ierr)
         endif
         if (lastyb==iproc) then
            fblo=fao
         else
            call MPI_IRECV(fblo,nbufx_gh,MPI_REAL,lastyb,tolastyb,MPI_COMM_WORLD,irecv_rq_fromlastyb,ierr)
         endif
         if (lastyb/=iproc) then
            call MPI_ISEND(fbo,nbufx_gh,MPI_REAL,lastyb,tonextya,MPI_COMM_WORLD,isend_rq_tolastyb,ierr)
         endif
         if (nextya==iproc) then
            fahi=fbo
         else
            call MPI_IRECV(fahi,nbufx_gh,MPI_REAL,nextya,tonextya,MPI_COMM_WORLD,irecv_rq_fromnextya,ierr)
         endif
         if (nextyb/=iproc) then
            call MPI_ISEND(fbo,nbufx_gh,MPI_REAL,nextyb,tolastya,MPI_COMM_WORLD,isend_rq_tonextyb,ierr)
         endif
         if (lastya==iproc) then
            falo=fbo
         else
            call MPI_IRECV(falo,nbufx_gh,MPI_REAL,lastya,tolastya,MPI_COMM_WORLD,irecv_rq_fromlastya,ierr)
         endif
      endif
    endsubroutine initiate_shearing
!***********************************************************************
    subroutine finalise_shearing(f)
!
!  Subroutine for shearing sheet boundary conditions
!
!  20-june-02/nils: adapted from pencil_mpi
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nghost,2*my-2*nghost,mz,mvar) :: fa, fb
      integer, dimension(MPI_STATUS_SIZE) :: irecv_stat_fal, irecv_stat_fan, irecv_stat_fbl, irecv_stat_fbn
      integer, dimension(MPI_STATUS_SIZE) :: isend_stat_tna, isend_stat_tla, isend_stat_tnb, isend_stat_tlb
      integer :: m2long,i
      double precision :: deltay_dy, frak, c1, c2, c3, c4, c5, c6
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
         fa(:,1:m2,:,:) = falo(:,1:m2,:,:)
         fa(:,m2+1:2*my-2*nghost,:,:) = fahi(:,m1:my,:,:)
         fb(:,1:m2,:,:) = fblo(:,1:m2,:,:)
         fb(:,m2+1:2*my-2*nghost,:,:) = fbhi(:,m1:my,:,:)
         displs = modulo(int(deltay_dy),ny)
         frak = deltay_dy - int(deltay_dy)
         c1 = -          (frak+1.)*frak*(frak-1.)*(frak-2.)*(frak-3.)/120.
         c2 = +(frak+2.)          *frak*(frak-1.)*(frak-2.)*(frak-3.)/24.
         c3 = -(frak+2.)*(frak+1.)     *(frak-1.)*(frak-2.)*(frak-3.)/12.
         c4 = +(frak+2.)*(frak+1.)*frak          *(frak-2.)*(frak-3.)/12.
         c5 = -(frak+2.)*(frak+1.)*frak*(frak-1.)          *(frak-3.)/24.
         c6 = +(frak+2.)*(frak+1.)*frak*(frak-1.)*(frak-2.)          /120.
         f(1:l1-1,m1:m2,:,:) = c1*fa(:,m2long-ny-displs+3:m2long-displs+2,:,:) &
              +c2*fa(:,m2long-ny-displs+2:m2long-displs+1,:,:) &
              +c3*fa(:,m2long-ny-displs+1:m2long-displs-0,:,:) &
              +c4*fa(:,m2long-ny-displs-0:m2long-displs-1,:,:) &
              +c5*fa(:,m2long-ny-displs-1:m2long-displs-2,:,:) &
              +c6*fa(:,m2long-ny-displs-2:m2long-displs-3,:,:)
         f(l2+1:mx,m1:m2,:,:) = c1*fb(:,m1+displs-2:m2+displs-2,:,:) &
              +c2*fb(:,m1+displs-1:m2+displs-1,:,:) &
              +c3*fb(:,m1+displs:m2+displs,:,:) &
              +c4*fb(:,m1+displs+1:m2+displs+1,:,:) &
              +c5*fb(:,m1+displs+2:m2+displs+2,:,:) &
              +c6*fb(:,m1+displs+3:m2+displs+3,:,:)
!
!  Filling also the x-y corners in order to avoid only zeros at these corners.
!  One should acctually have communicated with an extra processor in order to
!  fill these corners with the right values, but this does not seem to be
!  necessary.
!
         do i=1,nghost
            f(1:l1-1,i,:,:)=f(1:l1-1,m1+nghost-i,:,:)
            f(1:l1-1,m2+i,:,:)=f(1:l1-1,m2-i+1,:,:)
            f(l2+1:mx,i,:,:)=f(l2+1:mx,m1+nghost-i,:,:)
            f(l2+1:mx,m2+i,:,:)=f(l2+1:mx,m2-i+1,:,:)
         enddo
!
!  need to wait till buffer is empty before re-using it again
!
         if (nextyb/=iproc) call MPI_WAIT(isend_rq_tonextyb,isend_stat_tnb,ierr)
         if (lastyb/=iproc) call MPI_WAIT(isend_rq_tolastyb,isend_stat_tlb,ierr)
         if (nextya/=iproc) call MPI_WAIT(isend_rq_tonextya,isend_stat_tna,ierr)
         if (lastya/=iproc) call MPI_WAIT(isend_rq_tolastya,isend_stat_tla,ierr)
!
       endsubroutine finalise_shearing
!***********************************************************************
    subroutine mpibcast_int(ibcast_array,nbcast_array)
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: ibcast_array
!
      call MPI_BCAST(ibcast_array,nbcast_array,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    endsubroutine mpibcast_int
!***********************************************************************
    subroutine mpibcast_logical_scl(lbcast_array,nbcast_array)
!
      integer :: nbcast_array
      logical :: lbcast_array
!
      call MPI_BCAST(lbcast_array,nbcast_array,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
    endsubroutine mpibcast_logical_scl
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
    subroutine start_serialize()
!
!  Do block between start_serialize and end_serialize serially in iproc
!  order. root goes first, then sends proc1 permission, waits for succes,
!  then sends proc2 permisssion, waits for success, etc. 
!
!  19-nov-02/wolf: coded
!
      integer :: i,buf
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
!call mpibarrier()
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
!
      mpiwtime = MPI_WTIME()
    endfunction mpiwtime
!***********************************************************************
    function mpiwtick()
!
      double precision :: mpiwtick
!
      mpiwtick = MPI_WTICK()
    endfunction mpiwtick
!***********************************************************************
    subroutine stop_it(msg)
!
!  Print message and stop
!  6-nov-01/wolf: coded
!
      character (len=*) :: msg
!
      if (lroot) write(0,'(A,A)') 'STOPPED: ', msg
      call mpifinalize
      STOP
    endsubroutine stop_it
!***********************************************************************
    subroutine transp(a,var)
!
!  Doing the transpose of information distributed on several processors
!  Used for doing FFTs in the y and z directions.
!
!  03-sep-02/nils: coded
!  26-oct-02/axel: comments added
!
      real, dimension(ny,ny,nz) :: send_buf_y, recv_buf_y
      real, dimension(nz,ny,nz) :: send_buf_z, recv_buf_z
      real, dimension(nx,ny,nz) :: a
      real, dimension(nz) :: tmp_z
      real, dimension(ny) :: tmp_y
      integer :: i,j,sendc_y,recvc_y,sendc_z,recvc_z,px
      integer :: ytag=101,ztag=102,partner,ierr
      integer :: isend_rq_y,irecv_rq_y,isend_rq_z,irecv_rq_z
      character :: var
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
!  This routine is presently restricted to the case nxgrid=nygrid=nzgrid
!  Make sure this is the case
!
      if(nxgrid/=nygrid.or.nygrid/=nzgrid) then
        print*,'transp: need to have nxgrid=nygrid=nzgrid'
        call stop_it('Inconsistency: nxgrid/=nygrid.or.nygrid/=nzgrid')
      endif
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
          if(px/=ipy) then
            partner=px+ipz*nprocy ! = iproc + (px-ipy)
            if(ip<=6) print*,'MPICOMM: ipy,ipz,px,partner=',ipy,ipz,px,partner
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
          do i=1,ny
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
!
!  Send information to different processors (x-z transpose)
!  See the discussion above for why we use this communication pattern
        do px=0,nprocz-1
          if(px/=ipz) then
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
          do i=1,nz
            do j=i+1,nz
              tmp_y=a(i+px*nz,:,j)
              a(i+px*nz,:,j)=a(j+px*nz,:,i)
              a(j+px*nz,:,i)=tmp_y
            enddo
          enddo
        enddo
!
      else
        if (lroot) print*,'TRANSP: No clue what var=', var, 'is supposed to mean'
      endif
!
!  Synchronize; not strictly necessary, so Axel will prabably remove it..
!
      call mpibarrier()
!
    endsubroutine transp
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
  if(lroot .AND. ip<10) print*,'doing fft of x-component'
  call fft(a1,b1, nx*ny*nz, nx, nx,-1) ! x-direction
  call transp(a1,'y')
  call transp(b1,'y')
  call fft(a1,b1, nx*ny*nz, nx, nx,-1) ! y-direction
  call transp(a1,'z')
  call transp(b1,'z')
  call fft(a1,b1, nx*ny*nz, nx, nx,-1) ! z-direction

  ! Doing the y field
  if(lroot .AND. ip<10) print*,'doing fft of y-component'
  call fft(a2,b2, nx*ny*nz, nx, nx,-1) ! x-direction
  call transp(a2,'y')
  call transp(b2,'y')
  call fft(a2,b2, nx*ny*nz, nx, nx,-1) ! y-direction
  call transp(a2,'z')
  call transp(b2,'z')
  call fft(a2,b2, nx*ny*nz, nx, nx,-1) ! z-direction

  ! Doing the z field
  if(lroot .AND. ip<10) print*,'doing fft of z-component'
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
  if(lroot .AND. ip<10) print*,'fft has finished'
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

  if(lroot .AND. ip<10) print*,'doing three FFTs'
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
  if(lroot .AND. ip<10) print*,'fft has finished'
!
endsubroutine transform_i
!***********************************************************************
subroutine transform_fftpack(a_re,a_im)
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
  integer :: l,m,n
!
!  check whether nxgrid=nygrid=nzgrid
!
  if(nxgrid/=nygrid .or. nxgrid/=nzgrid) then
    print*,'transform_fftpack: must have nxgrid=nygrid=nzgrid!'
    call stop_it("must have nxgrid=nygrid=nzgrid!")
  endif
!
!  need to initialize cfft only once, because nxgrid=nygrid=nzgrid
!
  call cffti(nx,wsavex)
!
  if(lroot .AND. ip<10) print*,'doing FFTpack in x'
  do n=1,nz
  do m=1,ny
    ax=cmplx(a_re(:,m,n),a_im(:,m,n))
    call cfftf(nx,ax,wsavex)
    a_re(:,m,n)=real(ax)
    a_im(:,m,n)=aimag(ax)
  enddo
  enddo
  call transp(a_re,'y')
  call transp(a_im,'y')
!
!  The length of the array in the y-direction is nx
!  (remember: nxgrid=nygrid=nzgrid!)
!
  if(lroot .AND. ip<10) print*,'doing FFTpack in y'
  do n=1,nz
  do m=1,ny
    ax=cmplx(a_re(:,m,n),a_im(:,m,n))
    call cfftf(nx,ax,wsavex)
    a_re(:,m,n)=real(ax)
    a_im(:,m,n)=aimag(ax)
  enddo
  enddo
  call transp(a_re,'z')
  call transp(a_im,'z')
!
!  The length of the array in the z-direction is also nx
!
  if(lroot .AND. ip<10) print*,'doing FFTpack in z'
  do n=1,nz
  do m=1,ny
    ax=cmplx(a_re(:,m,n),a_im(:,m,n))
    call cfftf(nx,ax,wsavex)
    a_re(:,m,n)=real(ax)
    a_im(:,m,n)=aimag(ax)
  enddo
  enddo
!
!  Normalize
!
  a_re=a_re/nwgrid
  a_im=a_im/nwgrid
  if(lroot .AND. ip<10) print*,'fft has finished'
!
endsubroutine transform_fftpack
!***********************************************************************
subroutine transform_nr(a_re,a_im)
!
!  Subroutine to do Fourier transform using Numerical Recipes routine.
!  Note that this routine requires that nx, ny, and nz are powers of 2.
!  The routine overwrites the input data.
!
!  30-oct-02/axel: adapted from transform_fftpack for Numerical Recipes
!
  real,dimension(nx,ny,nz) :: a_re,a_im
  complex,dimension(nx) :: ax
  integer :: m,n
!
!  This Fourier transform would work, but it's very slow!
!  Even the compilation is very slow, so we better get rid of it!
!
  print*,'fft_nr currently disabled!'
  call stop_it("")
!
  if(lroot .AND. ip<10) print*,'doing FFT_nr in x'
  do n=1,nz
  do m=1,ny
    ax=cmplx(a_re(:,m,n),a_im(:,m,n))
    !call four1(ax,nx,-1)
    a_re(:,m,n)=real(ax)
    a_im(:,m,n)=aimag(ax)
  enddo
  enddo
  call transp(a_re,'y')
  call transp(a_im,'y')
!
!  The length of the array in the y-direction is nx
!  (remember: nxgrid=nygrid=nzgrid!)
!
  if(lroot .AND. ip<10) print*,'doing FFT_nr in y'
  do n=1,nz
  do m=1,ny
    ax=cmplx(a_re(:,m,n),a_im(:,m,n))
    !call four1(ax,nx,-1)
    a_re(:,m,n)=real(ax)
    a_im(:,m,n)=aimag(ax)
  enddo
  enddo
  call transp(a_re,'z')
  call transp(a_im,'z')
!
!  The length of the array in the z-direction is also nx
!
  if(lroot .AND. ip<10) print*,'doing FFT_nr in z'
  do n=1,nz
  do m=1,ny
    ax=cmplx(a_re(:,m,n),a_im(:,m,n))
    !call four1(ax,nx,-1)
    a_re(:,m,n)=real(ax)
    a_im(:,m,n)=aimag(ax)
  enddo
  enddo
!
!  Normalize
!
  a_re=a_re/nwgrid
  a_im=a_im/nwgrid
  if(lroot .AND. ip<10) print*,'fft has finished'
!
endsubroutine transform_nr
!***********************************************************************
endmodule Mpicomm
