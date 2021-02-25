program magic_emul

  include 'mpif.h'
  integer:: mpierr, ncpus=4,i
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real :: pi=3.1415, dx
  logical :: lok
  integer :: tag_foreign=1734
  real, dimension(64) :: xcoors
  real, dimension(64,64,64,3) :: uu_data
  real, dimension(:,:,:,:), allocatable :: buffer

  integer :: nprocs, iproc, iapp, flag, MPI_COMM_MAGIC, tag
  integer, dimension(2) :: xind_rng
!
      call MPI_INIT(mpierr)
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
      call MPI_COMM_GET_ATTR(MPI_COMM_WORLD, MPI_APPNUM, iapp, flag, mpierr)
!
! New comm MPI_COMM_PENCIL which comprises only the procs of the Pencil
! application. iproc becomes rank in MPI_COMM_PENCIL.
! Attention: If there is more than one application envisaged, Pencil needs to be
! compiled
! with FPPFLAGS=-fpp -DMPI_COMM_WORLD=MPI_COMM_PENCIL.
! If there is only one application, iproc is unchanged and
! MPI_COMM_PENCIL=MPI_COMM_WORLD.
!
      iproc_save=iproc
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, iapp, iproc, MPI_COMM_MAGIC, mpierr)
      call MPI_COMM_RANK(MPI_COMM_MAGIC, iproc, mpierr)
      if (iproc==0) then
!
!  Send length of name of foreign code.
!
          call MPI_SEND(5,1,MPI_INTEGER,0,tag_foreign, MPI_COMM_WORLD, mpierr)
!
!  Send name of foreign code.
!
          call MPI_SEND('MagIC',5,MPI_CHARACTER,0,tag_foreign, MPI_COMM_WORLD, mpierr)
!
!  Send processor numbers of foreign code.
!
          call MPI_SEND((/2,1,1/),3,MPI_INTEGER,0,tag_foreign, MPI_COMM_WORLD, mpierr)
!
!  Send gridpoint numbers of foreign code.
!
          call MPI_SEND((/64,64,64/),3,MPI_INTEGER,0,tag_foreign, MPI_COMM_WORLD, mpierr)
!
!  Send domain extents of foreign code. j loops over r, theta, phi.
!
          call MPI_SEND((/.7,1.,0.,pi,0.,2*pi/),6,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,mpierr)
!
!  Send output timestep of foreign code (code units).
!
          call MPI_SEND(1.e-3,1,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,mpierr)
!
!  Receive confirmation flag that setup is acceptable.
!
          call MPI_RECV(lok,1,MPI_LOGICAL,0,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
          if (.not.lok) then
            print*, 'not ok'
          endif
!
!  Send vector of global r-grid points.
!
          dx=0.3/63.
          do i=0,63
            xcoors(i+1)=.7+i*dx
          enddo 
          call MPI_SEND(xcoors,64,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,mpierr)
!
        endif
!
!  Receive index range of buddy processors.
!
        tag=tag_foreign+iproc
print*, 'emul: iproc_save, iproc, tag=', iproc_save, iproc, tag
        call MPI_RECV(xind_rng,2,MPI_INTEGER,iproc,tag,MPI_COMM__WORLD,stat,mpierr)
print*, 'emul: xind:rng=', xind_rng
        allocate(buffer(xind_rng(2):xind_rng(1),64,64,3))
        buffer=uu_data(xind_rng(1):xind_rng(2),:,:,:)

        call MPI_SEND(buffer,(xind_rng(2)-xind_rng(1)+1)*64*64*3, &
                      MPI_REAL,iproc,tag,MPI_COMM_WORLD,mpierr)

!print*, 'emul: iapp,iproc,iproc_save=', iapp,iproc,iproc_save
!print*, 'emul: successful'
      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
      call MPI_FINALIZE(mpierr)
stop
end
