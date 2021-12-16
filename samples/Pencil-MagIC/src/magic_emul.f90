program magic_emul

  implicit none
  include 'mpif.h'
  integer:: mpierr,ncpus,i,iproc,iproc_save,MPI_COMM_MAGIC, &
            nprocs,tag,nprocx_penc, nprocx, ip
  integer, dimension(MPI_STATUS_SIZE) :: stat
  real :: pi=3.1415, dx
  logical :: lok
  integer :: tag_foreign=1734
  real, dimension(64) :: xcoors
  real, dimension(64,64,64,3) :: uu_data

  INTEGER(KIND=MPI_ADDRESS_KIND) :: iapp
  integer, dimension(:,:), allocatable :: xind_rng
  logical :: flag
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
!print*, 'MAGIC: iapp,flag=', iapp, flag
! New comm MPI_COMM_MAGIC which comprises only the procs of the Pencil
! application. iproc becomes rank in MPI_COMM_MAGIC.
! Attention: If there is more than one application envisaged, Pencil needs to be
! compiled
! with FPPFLAGS=-fpp -DMPI_COMM_WORLD=MPI_COMM_MAGIC.
! If there is only one application, iproc is unchanged and
! MPI_COMM_MAGIC=MPI_COMM_WORLD.
!
      iproc_save=iproc
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, iapp, iproc, MPI_COMM_MAGIC, mpierr)
      call MPI_COMM_RANK(MPI_COMM_MAGIC, iproc, mpierr)
      call MPI_COMM_SIZE(MPI_COMM_MAGIC, ncpus, mpierr)
      nprocx=ncpus ! for MagIC only

!print*, '#Magic-World cpus, rank=', nprocs, iproc_save
!print*, '#Magic cpus, rank=', ncpus, iproc
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
          call MPI_SEND((/ncpus,1,1/),3,MPI_INTEGER,0,tag_foreign, MPI_COMM_WORLD, mpierr)
!
!  Send gridpoint numbers of foreign code.
!
          call MPI_SEND((/64,64,64/),3,MPI_INTEGER,0,tag_foreign, MPI_COMM_WORLD, mpierr)
!
!  Send domain extents of foreign code. j loops over r, theta, phi.
!
          call MPI_SEND((/.7,1.,0.,pi,0.,2*pi/),6,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,mpierr)
!
!  Send output timestep of foreign code.
!
          call MPI_SEND(1.e-3,1,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,mpierr)
!
!  Send unit system name and units.
!
          call MPI_SEND('SI   ',5,MPI_CHARACTER,0,tag_foreign,MPI_COMM_WORLD,mpierr)
          call MPI_SEND(1.   ,1,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,mpierr)
          call MPI_SEND(3600.,1,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,mpierr)
          call MPI_SEND(1.e-4,1,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,mpierr)
          call MPI_SEND(1.e6 ,1,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,mpierr)
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
            xcoors(i+1)=i*dx
          enddo 
          xcoors=atan(10.*xcoors)
          xcoors=.3*xcoors/maxval(xcoors)+.7
          call MPI_SEND(xcoors,64,MPI_REAL,0,tag_foreign,MPI_COMM_WORLD,mpierr)
!
!  Receive number of x-procs from Pencil.
!
          call MPI_RECV(nprocx_penc,1,MPI_INTEGER,ncpus,tag_foreign,MPI_COMM_WORLD,stat,mpierr)
!
        endif
!
!  Receive index range of buddy processors.
!
        if (nprocx_penc==nprocx .and. nprocx>1) then     ! non-EULAG case
          call MPI_BCAST(nprocx_penc,1,MPI_INTEGER,0,MPI_COMM_MAGIC,mpierr)
          tag=tag_foreign+iproc
          !call MPI_RECV(xind_rng,2,MPI_INTEGER,iproc,tag,MPI_COMM_WORLD,stat,mpierr)
print*, 'Magic: xind_rng=', xind_rng
        elseif (nprocx==1) then                          ! EULAG case
          if (iproc==0) then
            allocate(xind_rng(0:nprocx_penc-1,2))
            do ip=0,nprocx_penc-1
              tag=tag_foreign+ip
              call MPI_RECV(xind_rng(ip,:),2,MPI_INTEGER,ip,tag,MPI_COMM_WORLD,stat,mpierr)
            enddo
          endif
        endif

        if (nprocx==1.and.iproc==0) then                              ! EULAG case
          do ip=0,nprocx_penc-1
            tag=tag_foreign+ip
            call MPI_SEND(uu_data(:,:,xind_rng(ip,1):xind_rng(ip,2),:),(xind_rng(ip,2)-xind_rng(ip,1)+1)*64*64*3, &
                          MPI_REAL,ip,tag,MPI_COMM_WORLD,mpierr)
          enddo
        endif

!print*, 'Magic: iproc_save, iproc=', iproc_save, iproc
print*, 'Magic: successful'
      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
      call MPI_FINALIZE(mpierr)
stop
end
