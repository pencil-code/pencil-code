program pencil_emul

  implicit none
  include 'mpif.h'
  integer:: mpierr, ncpus,i
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real, parameter :: pi=3.1415
  logical :: lok
  integer :: tag_foreign=1734
  real, dimension(64,64,64,3) :: uu_data
  real, dimension(:,:,:,:), allocatable :: buffer

  integer :: nprocs, iproc, MPI_COMM_PENCIL, tag, iproc_save, nxdel, ind1, ind2, peer
  integer, parameter :: nprocx=1, nprocy=1, nprocz=1, nprocxy=nprocx*nprocy
  integer :: MPI_COMM_XBEAM, ipx, ipy, ipz
  integer, parameter :: nx=16
  INTEGER(KIND=MPI_ADDRESS_KIND) :: iapp
  logical :: flag,lnprocx_mismat
  integer, dimension(-1:0,2) :: xind_rng
  integer :: foreign_name_len
  character(LEN=5) :: foreign_name
  real, dimension(6) :: foreign_extents
  integer, dimension(3) :: nforeign_procs, nforeign_grid
  real :: foreign_dt, dx
  real, dimension(:), allocatable :: xcoors
  real, dimension(nx) :: x
!
  character(LEN=5) :: unit_system
  real :: unit_length, unit_time, unit_BB, unit_T

      call MPI_INIT(mpierr)
      call MPI_COMM_GET_ATTR(MPI_COMM_WORLD, MPI_APPNUM, iapp, flag, mpierr)
!
! Size and rank w.r.t. MPI_COMM_WORLD
!
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, iproc, mpierr)
!
! If mpirun/mpiexec calls also other applications than Pencil:
! Get rank within the set of applications, iapp.
! iapp=0 if there is only one application or Pencil is the first one.
!
!
! New comm MPI_COMM_PENCIL which comprises only the procs of the Pencil
! application. iproc becomes rank in MPI_COMM_PENCIL.
! Attention: If there is more than one application envisaged, Pencil
! needs to be
! compiled
! with FPPFLAGS=-fpp -DMPI_COMM_WORLD=MPI_COMM_PENCIL.
! If there is only one application, iproc is unchanged and
! MPI_COMM_PENCIL=MPI_COMM_WORLD.
!
      iproc_save=iproc
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, iapp, iproc, MPI_COMM_PENCIL,mpierr)
      call MPI_COMM_RANK(MPI_COMM_PENCIL, iproc, mpierr)
      call MPI_COMM_SIZE(MPI_COMM_PENCIL, ncpus, mpierr)
!print*, '#Pencil app,cpus,rank=', iapp,ncpus, iproc
      if (ncpus/=nprocx*nprocy*nprocz) then
        print*, 'Pencil: inconsistent proc numbers!'
        call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
        call MPI_FINALIZE(mpierr)
        stop
      endif

      ipx = modulo(iproc, nprocx)
      ipy = modulo(iproc/nprocx, nprocy)
      ipz = iproc/nprocxy

      call MPI_COMM_SPLIT(MPI_COMM_PENCIL, ipy+nprocy*ipz, ipx, &
                          MPI_COMM_XBEAM, mpierr)
        if (iproc==0) then
!
!  Receive length of name of foreign code.
!
          call MPI_RECV(foreign_name_len,1,MPI_INTEGER,ncpus,tag_foreign,MPI_COMM_WORLD, stat, mpierr)
!
!  Receive name of foreign code.
!
          call MPI_RECV(foreign_name,foreign_name_len,MPI_CHARACTER,ncpus,tag_foreign,MPI_COMM_WORLD, stat, mpierr)
!
!  Receive processor numbers of foreign code.
!
          call MPI_RECV(nforeign_procs,3,MPI_INTEGER,ncpus,tag_foreign,MPI_COMM_WORLD, stat, mpierr)
!
!  Receive gridpoint numbers of foreign code.
!
          call MPI_RECV(nforeign_grid,3,MPI_INTEGER,ncpus,tag_foreign,MPI_COMM_WORLD, stat, mpierr)
!
!  Receive domain extents of foreign code. j loops over r, theta, phi.
!
          call MPI_RECV(foreign_extents,6,MPI_REAL,ncpus,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
!
!  Receive output timestep of foreign code.
!
          call MPI_RECV(foreign_dt,1,MPI_REAL,ncpus,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
!
!  Send unit system name and units.
!
          call MPI_RECV(unit_system,5,MPI_CHARACTER,0,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
          call MPI_RECV(unit_length,1,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
          call MPI_RECV(unit_time,1,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
          call MPI_RECV(unit_BB,1,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
          call MPI_RECV(unit_T,1,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
!
!  Send confirmation flag that setup is acceptable.
!
          lok=.true.
          call MPI_SEND(lok,1,MPI_LOGICAL,ncpus,tag_foreign,MPI_COMM_WORLD,mpierr)
          if (.not.lok) then
            print*, 'Pencil: not ok'
          endif
        endif
!
!  Broadcast nforeign_grid
!
        if (iproc<nprocx) then
          call MPI_BCAST(nforeign_grid,3,MPI_INTEGER,0,MPI_COMM_XBEAM,mpierr)
          allocate(xcoors(nforeign_grid(1)))
        endif
        if (iproc==0) then
!
!  Receive vector of global r-grid points.
!
          call MPI_RECV(xcoors,nforeign_grid(1),MPI_REAL,ncpus,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
!
!print*, 'foreign_name_len, foreign_name=', foreign_name_len, foreign_name
!print*, 'nforeign_procs, nforeign_grid=', nforeign_procs, nforeign_grid
!print*, 'foreign_extents, foreign_dt=', foreign_extents, foreign_dt
!print*, 'xcoors=', xcoors  !(1:10)
!
!  Send number of x-procs to foreign.
!
          call MPI_SEND(nprocx,1,MPI_INTEGER,ncpus,tag_foreign,MPI_COMM_WORLD,mpierr)
!
        endif
!
!  Broadcast xcoors
!
        if (iproc<nprocx) then
          lnprocx_mismat = nprocx>1
          call MPI_BCAST(xcoors,nforeign_grid(1),MPI_REAL,0,MPI_COMM_XBEAM,mpierr)
!
!  Local grid, assume parallization only in x.
!
          dx=0.3/(nprocx*nx)
          x=rangegen(0,nx-1)*dx+dx/2 + iproc*nx*dx + .7

          ind1=find_index(xcoors,x(1))
          if (xcoors(ind1)>x(1).and.ind1>1) ind1=ind1-1 
          ind2=find_index(xcoors,x(nx))
          if (xcoors(ind2)<x(nx).and.ind2<nforeign_grid(1)) ind2=ind2+1 
          xind_rng(-1,:)=(/ind1,ind2/)     
          if (nforeign_procs(1)==1) xind_rng(0,:)=xind_rng(-1,:)      !EULAG case

!print*, 'Pencil: iproc, xind_rng=', iproc, xind_rng(-1,:), xcoors(ind1), xcoors(ind2)

          if (lnprocx_mismat) then
!
!  If number of procs in x-direction is unequal.
!
            tag=tag_foreign+iproc
            !call MPI_SEND(xind_rng,2,MPI_INTEGER,ncpus+iproc,tag,MPI_COMM_WORLD,mpierr)
            !allocate(buffer(xind_rng(1):xind_rng(2),64,64,3))
            !call MPI_RECV(buffer,(xind_rng(2)-xind_rng(1)+1)*64*64*3, &
            !              MPI_REAL,ncpus+iproc,tag,MPI_COMM_WORLD,stat, mpierr)
            ! EULAG case
            call MPI_SEND(xind_rng(-1,:),2,MPI_INTEGER,ncpus,tag,MPI_COMM_WORLD,mpierr)
          else
!
!  Send index range of buddy processors if number of procs in x-direction is equal.
!
            peer=iproc+ncpus
            call MPI_SEND(xind_rng(-1,:),2,MPI_INTEGER,peer,tag_foreign+iproc,MPI_COMM_WORLD,mpierr)
          endif
        endif
      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
      call MPI_FINALIZE(mpierr)
      print*, 'Pencil: successful'
contains
!****************************************************************************
    function rangegen(start,end,step)
!
! Generates vector of integers  start,...,end, analogous to IDL-rangegen.
!
! 5-feb-14/MR: coded
!
      integer, intent(in) :: start, end
      integer, intent(in), optional :: step
      integer, dimension(end-start+1) :: rangegen

      integer :: i

      do i=start,end
        rangegen(i-start+1)=i
      enddo

    endfunction rangegen
!****************************************************************************
    pure integer function find_index(xa, x, dx)
!
!  Returns the index of the element in array xa that is closest to x.
!  with dx present: if distance is bigger than dx/2: return value negative.
!
!  24-feb-13/ccyang: coded
!  14-may-20/MR: new optional argument dx
!               
      real, dimension(:), intent(in) :: xa
      real, intent(in) :: x
      real, intent(in), optional :: dx
!
      integer :: closest(1)
!
      closest = minloc(abs(x - xa))
      find_index = closest(1)

      if (present(dx)) then
        if (abs(x-xa(find_index))>.5*dx) find_index=-find_index
      endif
!
    endfunction find_index
!***********************************************************************
end
