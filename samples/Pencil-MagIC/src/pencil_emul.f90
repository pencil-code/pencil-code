program pencil_emul

  implicit none
  include 'mpif.h'
  integer:: mpierr, ncpus,i
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real :: pi=3.1415, dx
  logical :: lok
  integer :: tag_foreign=1734
  real, dimension(64) :: xcoors
  real, dimension(64,64,64,3) :: uu_data
  real, dimension(:,:,:,:), allocatable :: buffer

  integer :: nprocs, iproc, MPI_COMM_PENCIL, tag, iproc_save
  INTEGER(KIND=MPI_ADDRESS_KIND) :: iapp
  logical :: flag
  integer, dimension(2) :: xind_rng, foreign_name_len
  character(LEN=5) :: foreign_name
  real, dimension(6) :: foreign_extents
  integer, dimension(3) :: nforeign_procs, nforeign_grid
  real :: foreign_dt
!
      call MPI_INIT(mpierr)
      call MPI_COMM_GET_ATTR(MPI_COMM_WORLD, MPI_APPNUM, iapp, flag, mpierr)
print*, '#Pencil app, flag=', iapp, flag
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
print*, '#Pencil-World cpus,rank=', nprocs, iproc_save
print*, '#Pencil cpus,rank=', ncpus, iproc
      if (iproc==0) then
!
!  Receive length of name of foreign code.
!
          call MPI_RECV(foreign_name_len,1,MPI_INTEGER,ncpus,tag_foreign,MPI_COMM_WORLD, stat, mpierr)
!
!  Receive name of foreign code.
!
          call MPI_RECV(foreign_name,5,MPI_CHARACTER,ncpus,tag_foreign,MPI_COMM_WORLD, stat, mpierr)
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
!  Receive output timestep of foreign code (code units).
!
          call MPI_RECV(foreign_dt,1,MPI_REAL,ncpus,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
!
!  Send confirmation flag that setup is acceptable.
!
          call MPI_SEND(lok,1,MPI_LOGICAL,ncpus,tag_foreign,MPI_COMM_WORLD,mpierr)
          if (.not.lok) then
            print*, 'not ok'
          endif
!
!  Receive vector of global r-grid points.
!
          !!!call
          !MPI_RECV(xcoors,64,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
!
        endif
      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
      call MPI_FINALIZE(mpierr)
stop
!
!  Send index range of buddy processors.
!
        tag=tag_foreign+iproc
print*, 'emul: iproc_save, iproc, tag=', iproc_save, iproc, tag
        call MPI_SEND(xind_rng,2,MPI_INTEGER,iproc,tag,MPI_COMM_WORLD,mpierr)
print*, 'emul: xind:rng=', xind_rng
        allocate(buffer(xind_rng(2):xind_rng(1),64,64,3))
        buffer=uu_data(xind_rng(1):xind_rng(2),:,:,:)

        call MPI_RECV(buffer,(xind_rng(2)-xind_rng(1)+1)*64*64*3, &
                      MPI_REAL,iproc,tag,MPI_COMM_WORLD,stat, mpierr)

!print*, 'emul: iapp,iproc,iproc_save=', iapp,iproc,iproc_save
!print*, 'emul: successful'
      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
      call MPI_FINALIZE(mpierr)
stop
end
