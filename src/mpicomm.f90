! $Id$
!
!  This module takes care of MPI communication.
!
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
!
! CPARAM logical, parameter :: lmpicomm = .true.
!
!***************************************************************
!
module Mpicomm
!
  use Cdata
  use General, only: find_proc, ioptest,loptest
  use Yinyang
!
  implicit none
!
  include 'mpif.h'
  include 'mpicomm.h'
!
  interface mpisend_cmplx
    module procedure mpisend_cmplx_arr3
  endinterface
!
  interface mpirecv_cmplx
    module procedure mpirecv_cmplx_arr3
  endinterface
!
!  interface mpigather_and_out
!    module procedure mpigather_and_out_real
!    module procedure mpigather_and_out_cmplx
!  endinterface
!
  integer(kind=MPI_OFFSET_KIND) :: size_of_int = 0, size_of_real = 0, size_of_double = 0
  logical :: lcommunicate_y=.false.
!
!  For f-array processor boundaries
!
  real, dimension (:,:,:,:), allocatable :: lbufxi,ubufxi,lbufxo,ubufxo
  real, dimension (:,:,:,:), allocatable :: lbufyi,ubufyi,lbufyo,ubufyo
  real, dimension (:,:,:,:), allocatable :: lbufzi,ubufzi,lbufzo,ubufzo
  real, dimension (:,:,:,:), allocatable :: llbufi,lubufi,uubufi,ulbufi
  real, dimension (:,:,:,:), allocatable :: llbufo,lubufo,uubufo,ulbufo
  integer, dimension(4,2) :: bufsizes_yz
  integer, dimension(2,4,2) :: bufsizes_yz_corn
  integer, parameter :: INYU=1, INZU=2, INYL=3, INZL=4
  integer, parameter :: INUU=1, INLU=2, INLL=3, INUL=4
  integer, parameter :: IRCV=1, ISND=2
!
  real, dimension (nghost,my-2*nghost,mz,mcom) :: fahi,falo,fbhi,fblo         ! For shear
  real, dimension (nghost,my-2*nghost,mz,mcom) :: fahihi,falolo,fbhihi,fblolo ! For shear
  real, dimension (nghost,my-2*nghost,mz,mcom) :: fao,fbo                     ! For shear
  integer :: ipx_partner, displs                                              ! For shear
  integer :: nextnextya, nextya, lastya, lastlastya ! For shear
  integer :: nextnextyb, nextyb, lastyb, lastlastyb ! For shear
  integer :: llcorn,lucorn,uucorn,ulcorn            ! (the 4 corners in yz-plane)
  integer :: psfcrn,psbcrn,pnfcrn,pnbcrn            ! (the 4 'pole' corners)
  integer :: llcornr, lucornr, ulcornr, uucornr     ! Receiving corner ids
  integer :: llcorns, lucorns, ulcorns, uucorns     ! Sending corner ids
  integer :: nprocs, mpierr
  integer :: serial_level = 0
!
!  mpi tags
!
  integer, parameter :: tolowx=13,touppx=14,tolowy=3,touppy=4,tolowz=5,touppz=6 ! msg. tags
  integer, parameter :: TOll=7,TOul=8,TOuu=9,TOlu=10 ! msg. tags for corners
  !integer :: tolowx=13,touppx=14,tolowy=MPI_ANY_TAG,touppy=MPI_ANY_TAG,tolowz=MPI_ANY_TAG,touppz=MPI_ANY_TAG ! msg. tags
  !integer :: TOll=MPI_ANY_TAG,TOul=MPI_ANY_TAG,TOuu=MPI_ANY_TAG,TOlu=MPI_ANY_TAG ! msg. tags for corners
  integer, parameter :: io_perm=20,io_succ=21
  integer, parameter :: shiftn=22,shifts=23
  integer, parameter :: TOnf=24,TOnb=25,TOsf=26,TOsb=27 ! msg. tags for pole corners
!
!  mpi tags for radiation
!  the values for those have to differ by a number greater than maxdir=190
!  in order to have unique tags for each boundary and each direction
!
  integer, parameter :: Qtag_yz=250, Qtag_zx=300, Qtag_xy=350
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
  integer :: isend_rq_tolastlastya,isend_rq_tonextnextya, &
             irecv_rq_fromlastlastya,irecv_rq_fromnextnextya ! For shear
  integer :: isend_rq_tolastlastyb,isend_rq_tonextnextyb, &
             irecv_rq_fromlastlastyb,irecv_rq_fromnextnextyb ! For shear
!
  integer, dimension (MPI_STATUS_SIZE) :: isend_stat_tl,isend_stat_tu
  integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fl,irecv_stat_fu
  integer, dimension (MPI_STATUS_SIZE) :: isend_stat_Tll,isend_stat_Tul, &
                                          isend_stat_Tuu,isend_stat_Tlu
  integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_Fuu,irecv_stat_Flu, &
                                          irecv_stat_Fll,irecv_stat_Ful
  integer :: REAL_ARR_MAXSIZE
!
!  Data for Yin-Yang communication.
!
  logical :: lcorner_yz=.false.
  type (ind_coeffs), dimension(:), allocatable :: intcoeffs

  integer, parameter :: LEFT=1, MID=2, RIGHT=3, MIDY=MID, NIL=0
  integer :: MIDZ=MID, MIDC=MID
  integer :: len_cornbuf, len_cornstrip_y, len_cornstrip_z, ycornstart, zcornstart
  integer, dimension(3) :: yy_buflens
!
!  Data for communication with foreign code.
!
  type foreign_setup

    character(LEN=labellen) :: name
    integer :: ncpus,tag,root
    integer, dimension(:), allocatable :: recv_req
    real :: dt_out, t_last_recvd
    integer, dimension(2) :: xpeer_rng=-1
    integer, dimension(2) :: ypeer_rng=-1
    integer, dimension(3) :: procnums
    integer, dimension(3) :: dims
    real, dimension(2,3) :: extents
    integer, dimension(3) :: proc_multis
    real, dimension(:), allocatable :: xgrid,ygrid
    integer, dimension(:,:), allocatable :: xind_rng,yind_rng
    character(LEN=5) :: unit_system
    real :: unit_length, unit_time, unit_BB, unit_T, &
            renorm_UU, renorm_t, renorm_L
    real, dimension(mx,2) :: xweights
    real, dimension(my,2) :: yweights
    integer, dimension(mx) :: xinds
    integer, dimension(my) :: yinds
    logical, dimension(3) :: linterpol

  endtype foreign_setup
!
  type(foreign_setup) :: frgn_setup
!
  contains
!
!***********************************************************************
    subroutine mpicomm_init
!
!  Part of the initialization which can be done before input parameters are known.
!  Get processor number, number of procs, and whether we are root.
!
!  20-aug-01/wolf: coded
!  29-jul-2010/anders: separate subroutine
!
      use Syscalls, only: sizeof_real
!
!$    integer :: thread_support
      integer :: iapp=0         ! (KIND=ikind8) ?
      integer :: flag
      integer :: nprocs_penc
!
!$    call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,thread_support,mpierr)
!$    if (thread_support < MPI_THREAD_MULTIPLE) then
!$      print*, "Asked for multiple thread support=',MPI_THREAD_MULTIPLE,' but got: ",thread_support
!$      call die_gracefully
!$    endif
!$    if (.false.) &
      call MPI_INIT(mpierr)
!
! Size and rank w.r.t. MPI_COMM_WORLD
!
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, iproc, mpierr)
!if (lroot) print*, 'Pencil0: nprocs,MPI_COMM_WORLD', nprocs, MPI_COMM_WORLD
!
! If mpirun/mpiexec calls also other applications than Pencil:
! Get rank within the set of applications, iapp.
! iapp=0 if there is only one application or Pencil is the first one.
!
      call MPI_COMM_GET_ATTR(MPI_COMM_WORLD, MPI_APPNUM, iapp, flag, mpierr)
      iapp=0   !!!
!
! New comm MPI_COMM_PENCIL which comprises only the procs of the Pencil
! application. iproc becomes rank in MPI_COMM_PENCIL.
! Attention: If there is more than one application envisaged, Pencil needs to be compiled
! with FPPFLAGS=-fpp -DMPI_COMM_WORLD=MPI_COMM_PENCIL.
! If there is only one application, iproc is unchanged and MPI_COMM_PENCIL=MPI_COMM_WORLD.
!
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, iapp, iproc, MPI_COMM_PENCIL, mpierr)
      call MPI_COMM_SIZE(MPI_COMM_PENCIL, nprocs_penc, mpierr)
      call MPI_COMM_RANK(MPI_COMM_PENCIL, iproc, mpierr)
!
      lroot = (iproc==root)                              ! refers to root of MPI_COMM_PENCIL!
!
      if (iapp>0) call stop_it('in MPMD mode, Pencil must be the first application')

      if (sizeof_real() < 8) then
        mpi_precision = MPI_REAL
        MPI_CMPLX = MPI_COMPLEX
        if (lroot) then
          write (*,*) ""
          write (*,*) "====================================================================="
          write (*,*) "WARNING: using SINGLE PRECISION, you'd better know what you're doing!"
          write (*,*) "====================================================================="
          write (*,*) ""
        endif
        if (linterstellar) call stop_it('using single precision: stop and recompile')
      else
        mpi_precision = MPI_DOUBLE_PRECISION
        MPI_CMPLX = MPI_DOUBLE_COMPLEX
      endif

      call MPI_COMM_DUP(MPI_COMM_PENCIL,MPI_COMM_GRID,mpierr)
      iproc_world=iproc
!
!  Remember the sizes of some MPI elementary types.
!
      call MPI_TYPE_SIZE(MPI_INTEGER, size_of_int, mpierr)
      call stop_it_if_any(mpierr /= MPI_SUCCESS, "unable to find MPI_INTEGER size")
!
      call MPI_TYPE_SIZE(mpi_precision, size_of_real, mpierr)
      call stop_it_if_any(mpierr /= MPI_SUCCESS, "unable to find MPI real size")
!
      call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, size_of_double, mpierr)
      call stop_it_if_any(mpierr /= MPI_SUCCESS, "unable to find MPI_DOUBLE_PRECISION size")
!
!  Check consistency in processor layout.
!
      if (ncpus/=nprocx*nprocy*nprocz) then
        if (lroot) then
          print*, 'Laid out for ncpus (per grid) = ', ncpus, &
                  ', but nprocx*nprocy*nprocz=', nprocx*nprocy*nprocz
        endif
        call stop_it('mpicomm_init')
      endif

      if (ncpus/=nprocs_penc) then
        if (lroot) then
          print*, 'Laid out for ncpus (per grid) = ', ncpus, &
                  ', but called with nprocs_penc=', nprocs_penc
        endif
        call stop_it('mpicomm_init')
      endif
!
!  Check if parallelization and chosen grid numbers make sense.
!
      if ((nprocx>1.and.nxgrid==1).or. &
          (nprocy>1.and.nygrid==1).or. &
          (nprocz>1.and.nzgrid==1)) then
        if (lroot) &
          print*, 'Parallelization in a dimension with ngrid==1 does not work.'
        call stop_it('mpicomm_init')
      endif
      if (mod(nxgrid,nprocx)/=0.or. &
          mod(nygrid,nprocy)/=0.or. &
          mod(nzgrid,nprocz)/=0) then
        if (lroot) &
          print*, 'In each dimension the number of grid points has to be '// &
                  'divisible by the number of processors.'
        call stop_it('mpicomm_init')
      endif
!
!  Avoid overlapping ghost zones.
!
      if ((nx<nghost) .and. (nxgrid/=1)) &
           call stop_it('Overlapping ghost zones in x-direction: reduce nprocx')
      if ((ny<nghost) .and. (nygrid/=1)) &
           call stop_it('Overlapping ghost zones in y-direction: reduce nprocy')
      if ((nz<nghost) .and. (nzgrid/=1)) &
           call stop_it('Overlapping ghost zones in z-direction: reduce nprocz')
!
      call MPI_TYPE_CONTIGUOUS(max_int, mpi_precision, REAL_ARR_MAXSIZE, mpierr)
      call MPI_TYPE_COMMIT(REAL_ARR_MAXSIZE, mpierr)

    endsubroutine mpicomm_init
!***********************************************************************
    subroutine initialize_mpicomm
!
!  Part of the initialization which can only be done after input parameters are known.
!  Initialise MPI communication and set up some variables.
!  The arrays leftneigh and rghtneigh give the processor numbers
!  to the left and to the right.
!
!  20-aug-01/wolf: coded
!  31-aug-01/axel: added to 3-D
!  15-sep-01/axel: adapted from Wolfgang's version
!  21-may-02/axel: communication of corners added
!   6-jun-02/axel: generalized to allow for ny=1
!  23-nov-02/axel: corrected problem with ny=4 or less
!  21-oct-15/fred: adapted periodic y-boundary across the pole
!  20-dec-15/MR: extensions for Yin-Yang grid.
!
!  Announce myself for pc_run to detect.
!
      use General, only: itoa, find_proc, find_proc_coords

      integer :: nprocs_penc, nprocs_foreign

      if (lroot) print *, 'initialize_mpicomm: enabled MPI'
      if (tag_foreign>0) lforeign=.true.
!
!  Check total number of processors assigned to Pencil.
!
      call MPI_COMM_SIZE(MPI_COMM_PENCIL, nprocs_penc, mpierr)
      nprocs_foreign=nprocs-nprocs_penc

      if (.not.( lyinyang.and.nprocs_penc==2*ncpus .or. lcubed_sphere.and.nprocs_penc==6*ncpus .or. &
                 (lforeign.and.(nprocs_foreign>0.or.lstart) .or. .not.lforeign).and.nprocs_penc==ncpus )) then
        if (lroot) then
          if (lyinyang) then
            print*, 'Compiled with 2*ncpus = ', 2*ncpus, &
                ', but running on ', nprocs_penc, ' processors'
          elseif (lcubed_sphere) then
            print*, 'Compiled with 6*ncpus = ', 6*ncpus, &
                ', but running on ', nprocs_penc, ' processors'
          elseif (nprocs_penc/=ncpus) then
            print*, 'Compiled with ncpus = ', ncpus, &
                ', but running on ', nprocs_penc, ' processors'
          elseif (lforeign.and.nprocs_foreign==0) then
            print*, 'number of processors = '//trim(itoa(nprocs))// &
                    ' = number of own processors -> no procs for foreign code left'
          endif
          call stop_it('initialize_mpicomm')
        endif
      endif
!
      if (lyinyang) then

        if (lforeign) then
          print*, 'Yin-Yang grid and foreign companion mutually exclusive!'
          call stop_it('initialize_mpicomm')
        endif

        if (nprocy>5.and..not.lcutoff_corners) then
          if (lroot) print*, 'Processor layout with nprocz>15, nprocy>5 not implemented for Yin-Yang grid.'
          call stop_it('initialize_mpicomm')
        endif
        if (nprocz/=3*nprocy) then
          if (lroot) print*, 'Processor layout with nprocz/=3*nprocy not implemented for Yin-Yang grid.'
          call stop_it('initialize_mpicomm')
        endif

        lyang=iproc>=ncpus             ! if this proc is in first half of all it is in YIN grid otherwise in YANG grid
      endif

      if (lyinyang .or. lcubed_sphere) &
        call MPI_COMM_SPLIT(MPI_COMM_PENCIL, int(iproc/ncpus), mod(iproc,ncpus), MPI_COMM_GRID, mpierr)

!  MPI_COMM_GRID refers to Yin or Yang grid or to PencilCode, when launched as first code (mandatory).
!
      if (lyinyang) then
        if (lyang) then
          iproc=iproc-ncpus
          cyinyang='Yang'
        else
          cyinyang='Yin'
        endif
      endif

!  MPI_COMM_GRID refers to one of the 6 patches of the cubed sphere grid or to PencilCode, when launched as first code (mandatory).
!
      if (lcubed_sphere) then
        iproc=modulo(iproc,ncpus)
        ipatch=int(iproc/ncpus)
      endif
!
!  Position on the processor grid (WITHIN Yin or Yang grid!).
!  x is fastest direction, z slowest (this is the default)
!
      call find_proc_coords(iproc,ipx,ipy,ipz)
      !if (lprocz_slowest) then
      !  ipx = modulo(iproc, nprocx)
      !  ipy = modulo(iproc/nprocx, nprocy)
      !  ipz = iproc/nprocxy
      !else
      !  ipx = modulo(iproc, nprocx)
      !  ipy = iproc/nprocxz
      !  ipz = modulo(iproc/nprocx, nprocz)
      !endif
!
!  Set up flags for leading processors in each possible direction and plane
!
      lfirst_proc_x = (ipx == 0)
      lfirst_proc_y = (ipy == 0)
      lfirst_proc_z = (ipz == 0)
      lfirst_proc_xy = lfirst_proc_x .and. lfirst_proc_y
      lfirst_proc_yz = lfirst_proc_y .and. lfirst_proc_z
      lfirst_proc_xz = lfirst_proc_x .and. lfirst_proc_z
      lfirst_proc_xyz = lfirst_proc_x .and. lfirst_proc_y .and. lfirst_proc_z
!
!  Set up flags for trailing processors in each possible direction and plane
!
      llast_proc_x = (ipx == nprocx-1)
      llast_proc_y = (ipy == nprocy-1)
      llast_proc_z = (ipz == nprocz-1)
      llast_proc_xy = llast_proc_x .and. llast_proc_y
      llast_proc_yz = llast_proc_y .and. llast_proc_z
      llast_proc_xz = llast_proc_x .and. llast_proc_z
      llast_proc_xyz = llast_proc_x .and. llast_proc_y .and. llast_proc_z
!
!  Am I a yz corner?
!
      lcorner_yz=(lfirst_proc_y.or.llast_proc_y).and.(lfirst_proc_z.or.llast_proc_z)
!
!  Set up `lower' and `upper' neighbours, refer to MPI_COMM_GRID.
!
      xlneigh = find_proc(ipx-1,ipy,ipz)  !modulo(ipx-1,nprocx) + ipy*nprocx + ipz*nprocxy
      xuneigh = find_proc(ipx+1,ipy,ipz)  !modulo(ipx+1,nprocx) + ipy*nprocx + ipz*nprocxy
      ylneigh = find_proc(ipx,ipy-1,ipz)  !ipx + modulo(ipy-1,nprocy)*nprocx + ipz*nprocxy
      yuneigh = find_proc(ipx,ipy+1,ipz)  !ipx + modulo(ipy+1,nprocy)*nprocx + ipz*nprocxy
      zlneigh = find_proc(ipx,ipy,ipz-1)  !ipx + ipy*nprocx + modulo(ipz-1,nprocz)*nprocxy
      zuneigh = find_proc(ipx,ipy,ipz+1)  !ipx + ipy*nprocx + modulo(ipz+1,nprocz)*nprocxy
!
!  Set the four corners in the yz-plane (in cyclic order).
!
      llcorn = find_proc(ipx,ipy-1,ipz-1)  !ipx + modulo(ipy-1,nprocy)*nprocx + modulo(ipz-1,nprocz)*nprocxy
      ulcorn = find_proc(ipx,ipy+1,ipz-1)  !ipx + modulo(ipy+1,nprocy)*nprocx + modulo(ipz-1,nprocz)*nprocxy
      uucorn = find_proc(ipx,ipy+1,ipz+1)  !ipx + modulo(ipy+1,nprocy)*nprocx + modulo(ipz+1,nprocz)*nprocxy
      lucorn = find_proc(ipx,ipy-1,ipz+1)  !ipx + modulo(ipy-1,nprocy)*nprocx + modulo(ipz+1,nprocz)*nprocxy
!
!  Overwrite with settings for boundary processors in Yin-Yang grid.
!
      if (lyinyang) call set_yy_neighbors
!
!  Settings for Yin-Yang grid are done later in yyinit.
!
      bufsizes_yz(:,IRCV)=(/nz,ny,nz,ny/)
      bufsizes_yz(:,ISND)=(/nz,ny,nz,ny/)

      bufsizes_yz_corn=nghost

      lcommunicate_y = (nprocz>1 .and. lpole(2))

      if (.not.lyinyang) then
!
!  Allocations for Yin-Yang grid are done later in yyinit.
!
        if (nprocx>1) &
          allocate( lbufxi(nghost,ny,nz,mcom),ubufxi(nghost,ny,nz,mcom), &
                    lbufxo(nghost,ny,nz,mcom),ubufxo(nghost,ny,nz,mcom))

        if (nprocy>1 .or. lcommunicate_y .or. lyinyang) &
          allocate( lbufyi(mx,nghost,bufsizes_yz(INYL,IRCV),mcom),ubufyi(mx,nghost,bufsizes_yz(INYU,IRCV),mcom), &
                    lbufyo(mx,nghost,bufsizes_yz(INYL,ISND),mcom),ubufyo(mx,nghost,bufsizes_yz(INYU,ISND),mcom))

        if (nprocz>1) &
          allocate( lbufzi(mx,bufsizes_yz(INZL,IRCV),nghost,mcom),ubufzi(mx,bufsizes_yz(INZU,IRCV),nghost,mcom), &
                    lbufzo(mx,bufsizes_yz(INZL,ISND),nghost,mcom),ubufzo(mx,bufsizes_yz(INZU,ISND),nghost,mcom))

        if ((nprocy>1 .or. lcommunicate_y .or. lyinyang).and.nprocz>1) &
          allocate( llbufi(mx,bufsizes_yz_corn(1,INLL,IRCV),bufsizes_yz_corn(2,INLL,IRCV),mcom), &
                    llbufo(mx,bufsizes_yz_corn(1,INLL,ISND),bufsizes_yz_corn(2,INLL,ISND),mcom), &
                    lubufi(mx,bufsizes_yz_corn(1,INLU,IRCV),bufsizes_yz_corn(2,INLU,IRCV),mcom), &
                    lubufo(mx,bufsizes_yz_corn(1,INLU,ISND),bufsizes_yz_corn(2,INLU,ISND),mcom), &
                    ulbufi(mx,bufsizes_yz_corn(1,INUL,IRCV),bufsizes_yz_corn(2,INUL,IRCV),mcom), &
                    ulbufo(mx,bufsizes_yz_corn(1,INUL,ISND),bufsizes_yz_corn(2,INUL,ISND),mcom), &
                    uubufi(mx,bufsizes_yz_corn(1,INUU,IRCV),bufsizes_yz_corn(2,INUU,IRCV),mcom), &
                    uubufo(mx,bufsizes_yz_corn(1,INUU,ISND),bufsizes_yz_corn(2,INUU,ISND),mcom)   )
      endif
!
!  Set standard processor ids
!
      llcornr = llcorn; llcorns = llcorn
      ulcornr = ulcorn; ulcorns = ulcorn
      uucornr = uucorn; uucorns = uucorn
      lucornr = lucorn; lucorns = lucorn
!
!  Print neighbors in counterclockwise order (including the corners),
!  starting with left neighbor.
!
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
!  Print processor numbers and those of their neighbors.
!
      if (ip<=7) write(6,'(A,I4,"(",3I4,"): ",8I4)') &
        'initialize_mpicomm: MPICOMM neighbors ', &
        iproc_world,ipx,ipy,ipz, &
        ylneigh,llcorn,zlneigh,ulcorn,yuneigh,uucorn,zuneigh,lucorn
!
!  Define MPI communicators that include all processes sharing the same value
!  of ipx, ipy, or ipz. Refer to MPI_COMM_GRID!
!  The rank within MPI_COMM_GRID is given by a
!  combination of the two other directional processor indices.
!
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipy+nprocy*ipz, ipx, &
                          MPI_COMM_XBEAM, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipx+nprocx*ipz, ipy, &
                          MPI_COMM_YBEAM, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipx+nprocx*ipy, ipz, &
                          MPI_COMM_ZBEAM, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipz, ipx+nprocx*ipy, &
                          MPI_COMM_XYPLANE, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipy, ipx+nprocx*ipz, &
                          MPI_COMM_XZPLANE, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipx, ipy+nprocy*ipz, &
                          MPI_COMM_YZPLANE, mpierr)
!
    endsubroutine initialize_mpicomm
!***********************************************************************
    subroutine create_communicators()

      !used to copy communicators for helper threads
      !this makes MPI communication thread-safe, since all threads have their own communicator

      !For code coupling with f.e. EULAG you have to split MPI_COMM_PENCIL
      !Not done now simply because its simpler this way
      call MPI_COMM_DUP(MPI_COMM_WORLD, MPI_COMM_PENCIL, mpierr)
      call MPI_COMM_DUP(MPI_COMM_PENCIL, MPI_COMM_GRID, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipy+nprocy*ipz, ipx, &
                          MPI_COMM_XBEAM, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipx+nprocx*ipz, ipy, &
                          MPI_COMM_YBEAM, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipx+nprocx*ipy, ipz, &
                          MPI_COMM_ZBEAM, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipz, ipx+nprocx*ipy, &
                          MPI_COMM_XYPLANE, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipy, ipx+nprocx*ipz, &
                          MPI_COMM_XZPLANE, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_GRID, ipx, ipy+nprocy*ipz, &
                          MPI_COMM_YZPLANE, mpierr)
    endsubroutine
!***********************************************************************
    subroutine update_neighbors
!
! Update neighbor processes for communication.
!
! 27-feb-16/ccyang: adapted from particles_mpicomm and
!                   particles_mpicomm_blocks
!  3-mar-16/MR: simplified sorting code
!
      use General, only: find_proc, quick_sort
!
      integer, dimension(:), allocatable :: order
      logical :: left, right
      integer :: dipx, dipy, dipz
      integer :: ipx_rec, ipy_rec, ipz_rec
      integer :: iproc_rec
!
!  Collect the neighbors.
!
      iproc_comm = -1
      nproc_comm = 0
!
      xscan: do dipx = -1, 1
        ipx_rec = ipx + dipx
        left = ipx_rec < 0
        right = ipx_rec > nprocx - 1
        ipx_rec = modulo(ipx_rec, nprocx)  ! assuming periodic boundary conditions
!
        yscan: do dipy = -1, 1
          ipy_rec = ipy + dipy
          shear: if (lshear) then
            if (left) ipy_rec = ipy_rec - ceiling(deltay / Lxyz_loc(2) - 0.5)
            if (right) ipy_rec = ipy_rec + ceiling(deltay / Lxyz_loc(2) - 0.5)
          endif shear
          ipy_rec = modulo(ipy_rec, nprocy)  ! assuming periodic boundary conditions
!
          zscan: do dipz = -1, 1
            ipz_rec = modulo(ipz + dipz, nprocz)  ! assuming periodic boundary conditions
!
            iproc_rec = find_proc(ipx_rec, ipy_rec, ipz_rec)
            neighbors(dipx,dipy,dipz) = iproc_rec
            add: if (iproc_rec /= iproc .and. .not. any(iproc_comm(1:nproc_comm) == iproc_rec)) then
              nproc_comm = nproc_comm + 1
              iproc_comm(nproc_comm) = iproc_rec
            endif add
!
          enddo zscan
!
        enddo yscan
!
      enddo xscan
!
!  Sort the process list.
!
      neighb: if (nproc_comm > 0) then
        allocate (order(nproc_comm))
        call quick_sort(iproc_comm(:nproc_comm), order)
        deallocate (order)
      endif neighb
!
    endsubroutine update_neighbors
!***********************************************************************
    elemental integer function index_to_iproc_comm(iproc_in, mask)
!
!  Converts iproc_in to the index to iproc_comm, returns 0 if iproc_in
!  is iproc itself, and -1 if none of the elements in iproc_comm matches
!  iproc_in.
!
!  iproc_in itself is returned if mask = .false..
!
!  29-feb-16/ccyang: coded.
!   3-mar-16/MR: use loptest
!
      use General, only: binary_search, loptest
!
      integer, intent(in) :: iproc_in
      logical, intent(in), optional :: mask
!
      active: if (loptest(mask,.true.)) then
        nonlocal: if (iproc_in /= iproc) then
!
!  Binary search iproc_comm.
!
          search: if (any(iproc_comm(1:nproc_comm) == iproc_in)) then
            index_to_iproc_comm = binary_search(iproc_in, iproc_comm(1:nproc_comm))
          else
            index_to_iproc_comm = -1
          endif search
        else nonlocal
!
!  The case of iproc_in == iproc.
!
          index_to_iproc_comm = 0
        endif nonlocal
      else active
!
!  Do nothing if mask = .false.
!
        index_to_iproc_comm = iproc_in
      endif active
!
    endfunction index_to_iproc_comm
!***********************************************************************
    subroutine set_cubed_sphere_neighbors
!
!  Cubed mesh
!
!  20-dec-15/MR: coded
!
      real, dimension(:,:,:), allocatable :: gridbuf_midy, gridbuf_midz, &  ! contains grid request of direct neighbour(s)
                                             gridbuf_left, &                !             ~         of left corner neighbour
                                             gridbuf_right                  !             ~         of right corner neighbour

      integer :: patch_neigh_left, patch_neigh_right, patch_neigh_top, patch_neigh_bot
      integer :: ipatch
      logical, save :: lcalled=.false.

      if (lcalled) then
        return
      else
        lcalled=.true.
      endif
!
      if (ipatch<5) then
        patch_neigh_right=modulo(ipatch  ,4)+1
        patch_neigh_left =modulo(ipatch+2,4)+1
        patch_neigh_top  =ZPLUS
        patch_neigh_bot  =ZMINUS
        if (llast_proc_y)  yuneigh=find_proc(ipx,       0,ipz)+(patch_neigh_right-1)*ncpus
        if (lfirst_proc_y) ylneigh=find_proc(ipx,nprocy-1,ipz)+(patch_neigh_left -1)*ncpus
        if (ipatch==1 .or. ipatch==3) then
          if (llast_proc_z)  zuneigh=find_proc(ipx,ipy,       0)+(patch_neigh_top  -1)*ncpus
          if (lfirst_proc_z) zlneigh=find_proc(ipx,ipy,nprocz-1)+(patch_neigh_bot  -1)*ncpus
        else
          if (llast_proc_z)  zuneigh=find_proc(ipx,       0,ipy)+(patch_neigh_top  -1)*ncpus
          if (lfirst_proc_z) zlneigh=find_proc(ipx,nprocy-1,ipy)+(patch_neigh_bot  -1)*ncpus
        endif

print*,'AXEL: patch_neigh_left, patch_neigh_right, patch_neigh_top, patch_neigh_bot=', &
              patch_neigh_left, patch_neigh_right, patch_neigh_top, patch_neigh_bot
      elseif (ipatch==5) then
      elseif (ipatch==6) then
      endif
!
    endsubroutine set_cubed_sphere_neighbors
!***********************************************************************
    subroutine scatter_snapshot(a,f,indvar1,indvar2)
!
!  Scatters a full snapshot (array a) residing in root (w/o ghost zones) to the f-arrays of all ranks
!  (not filling ghost zones there). It's assumed that a contains indvar2-indvar1+1 variables to be stored
!  in f(:,:,:,indvar1:indvar2).
!
!  23-oct-23/MR: coded
!
      real, dimension(:,:,:,:), intent(in):: a
      real, dimension(mx,my,mz,mfarray), intent(out):: f
      integer, intent(in) :: indvar1,indvar2

      integer, dimension(4) :: start_get, start_store
      integer :: type_get, type_store, win
      INTEGER(KIND=MPI_ADDRESS_KIND) :: size

      nvar=indvar2-indvar1+1
      start_get=(/ipx*nx,ipy*ny,ipz*nz,0/)+1
      start_store=(/ipx*nx+1+nghost,ipy*ny+1+nghost,ipz*nz+1+nghost,indvar1/)

      call MPI_Type_create_subarray(4, (/nxgrid,nygrid,nzgrid,nvar/), (/nx,ny,nz,nvar/), &
                                    start_get, MPI_ORDER_FORTRAN, MPI_REAL,type_get,mpierr)
      call MPI_Type_create_subarray(4, (/mx,my,mz,mfarray/), (/nx,ny,nz,nvar/), &
                                    start_store, MPI_ORDER_FORTRAN, MPI_REAL,type_store,mpierr)
      call MPI_Type_commit(type_get,mpierr)
      call MPI_Type_commit(type_store,mpierr)

      size=nwgrid*size_of_real

      if (lroot) then   ! data tb scattered is on root
        call MPI_Win_create(a, size, 1, MPI_INFO_NULL, MPI_COMM_WORLD, win, mpierr)
      else              ! hence other ranks do not declare a window
        call MPI_Win_create(a, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, win, mpierr)
      endif

      call MPI_WIN_FENCE(0, win, mpierr)
      call MPI_Get(f, 1, type_store, root, 0, 1, type_get, win, mpierr)
      call MPI_WIN_FENCE(0, win, mpierr)

      call MPI_TYPE_FREE(type_get,mpierr)
      call MPI_TYPE_FREE(type_store,mpierr)
      call MPI_WIN_FREE(win,mpierr)

    endsubroutine scatter_snapshot
!***********************************************************************
    subroutine yyinit
!
!  For Yin-Yang grid runs: allocate all communication buffers both for
!  initialization of interpolation (gridbuf_*) and real data exchange.
!  Requires set_yy_neighbors to be run first.
!
!  20-dec-15/MR: coded
!
      real, dimension(:,:,:), allocatable :: gridbuf_midy, gridbuf_midz, &! contains grid request of direct neighbour(s)
                                             gridbuf_left, &              !             ~         of left corner neighbour
                                             gridbuf_right                !             ~         of right corner neighbour
      logical, save :: lcalled=.false.

      if (lcalled) then
        return
      else
        lcalled=.true.
      endif
!
      if (lcorner_yz) then
!
!  Procs at yz corners have four outer neighbors: two direct ones for both corner edges
!  and two corner neighbors in the directions away from the yz corner.
!  (A fifth neighbor at the yz corner could be needed, too, for ncpus>75.)
!
        if (nprocy>=6) then
          allocate(intcoeffs(5))
        else
          allocate(intcoeffs(4))
        endif
        MIDZ=4
      else
!
!  Otherwise only three neighbors.
!
        allocate(intcoeffs(3))
        MIDZ=2
      endif
!
!  Set the grid transmission buffer sizes and allocate buffers.
!
      if (lfirst_proc_z.or.llast_proc_z) then
        if (lfirst_proc_z) then
          bufsizes_yz(INZL,ISND)=yy_buflens(MID)
          bufsizes_yz_corn(:,INLL,ISND)=(/yy_buflens(RIGHT),nghost/)
          bufsizes_yz_corn(:,INUL,ISND)=(/yy_buflens(LEFT),nghost/)
        endif
        if (llast_proc_z ) then
          bufsizes_yz(INZU,ISND)=yy_buflens(MID)
          bufsizes_yz_corn(:,INLU,ISND)=(/yy_buflens(LEFT),nghost/)
          bufsizes_yz_corn(:,INUU,ISND)=(/yy_buflens(RIGHT),nghost/)
        endif

        if (.not.((llast_proc_y.and.lfirst_proc_z).or.(llast_proc_z.and.lfirst_proc_y))) &
          allocate(gridbuf_left(2,yy_buflens(LEFT),nghost))

        if (.not.((lfirst_proc_y.and.lfirst_proc_z).or.(llast_proc_z.and.llast_proc_y))) &
          allocate(gridbuf_right(2,yy_buflens(RIGHT),nghost))

        allocate(gridbuf_midy(2,yy_buflens(MID),nghost))
      endif

      if (lfirst_proc_y.or.llast_proc_y) then
        if (lfirst_proc_y) then
          bufsizes_yz(INYL,ISND)=yy_buflens(MID)
          bufsizes_yz_corn(:,INLL,ISND)=(/nghost,yy_buflens(LEFT)/)
          bufsizes_yz_corn(:,INLU,ISND)=(/nghost,yy_buflens(RIGHT)/)
        endif
        if (llast_proc_y ) then
          bufsizes_yz(INYU,ISND)=yy_buflens(MID)
          bufsizes_yz_corn(:,INUU,ISND)=(/nghost,yy_buflens(LEFT)/)
          bufsizes_yz_corn(:,INUL,ISND)=(/nghost,yy_buflens(RIGHT)/)
        endif

        if (.not.((lfirst_proc_y.and.lfirst_proc_z).or.(llast_proc_z.and.llast_proc_y))) &
          allocate(gridbuf_left(2,nghost,yy_buflens(LEFT)))

        if (.not.((llast_proc_y.and.lfirst_proc_z).or.(llast_proc_z.and.lfirst_proc_y))) &
          allocate(gridbuf_right(2,nghost,yy_buflens(RIGHT)))

        allocate(gridbuf_midz(2,nghost,yy_buflens(MID)))

      endif
!
!  Create interpolation data.
!
      call setup_interp_yy(gridbuf_left, gridbuf_midy, gridbuf_midz, gridbuf_right)
!
!  Allocate data communication buffers.
!
      allocate( lbufyi(mx,nghost,bufsizes_yz(INYL,IRCV),mcom),ubufyi(mx,nghost,bufsizes_yz(INYU,IRCV),mcom), &
                lbufzi(mx,bufsizes_yz(INZL,IRCV),nghost,mcom),ubufzi(mx,bufsizes_yz(INZU,IRCV),nghost,mcom))

      if (lfirst_proc_z) then
        allocate(lbufzo(mx,nghost,bufsizes_yz(INZL,ISND),mcom))
      else
        allocate(lbufzo(mx,bufsizes_yz(INZL,ISND),nghost,mcom))
      endif

      if (llast_proc_z) then
        allocate(ubufzo(mx,nghost,bufsizes_yz(INZU,ISND),mcom))
      else
        allocate(ubufzo(mx,bufsizes_yz(INZU,ISND),nghost,mcom))
      endif

      if (lfirst_proc_y.and.ipz>=nprocz/3.and.ipz<2*nprocz/3.) then
        allocate(lbufyo(mx,bufsizes_yz(INYL,ISND),nghost,mcom))
      else
        allocate(lbufyo(mx,nghost,bufsizes_yz(INYL,ISND),mcom))
      endif

      if (llast_proc_y.and.ipz>=nprocz/3.and.ipz<2*nprocz/3.) then
        allocate(ubufyo(mx,bufsizes_yz(INYU,ISND),nghost,mcom))
      else
        allocate(ubufyo(mx,nghost,bufsizes_yz(INYU,ISND),mcom))
      endif

      allocate( llbufi(mx,bufsizes_yz_corn(1,INLL,IRCV),bufsizes_yz_corn(2,INLL,IRCV),mcom), &
                lubufi(mx,bufsizes_yz_corn(1,INLU,IRCV),bufsizes_yz_corn(2,INLU,IRCV),mcom), &
                ulbufi(mx,bufsizes_yz_corn(1,INUL,IRCV),bufsizes_yz_corn(2,INUL,IRCV),mcom), &
                uubufi(mx,bufsizes_yz_corn(1,INUU,IRCV),bufsizes_yz_corn(2,INUU,IRCV),mcom))

      if (lfirst_proc_z.or.lfirst_proc_y.and.ipz>=nprocz/3.and.ipz<2*nprocz/3.) then
        allocate(llbufo(mx,bufsizes_yz_corn(2,INLL,ISND),bufsizes_yz_corn(1,INLL,ISND),mcom))
      else
        allocate(llbufo(mx,bufsizes_yz_corn(1,INLL,ISND),bufsizes_yz_corn(2,INLL,ISND),mcom))
      endif

      if (lfirst_proc_z.or.llast_proc_y.and.ipz>=nprocz/3.and.ipz<2*nprocz/3.) then
        allocate(ulbufo(mx,bufsizes_yz_corn(2,INUL,ISND),bufsizes_yz_corn(1,INUL,ISND),mcom) )
      else
        allocate(ulbufo(mx,bufsizes_yz_corn(1,INUL,ISND),bufsizes_yz_corn(2,INUL,ISND),mcom) )
      endif

      if (llast_proc_z.or.llast_proc_y.and.ipz>=nprocz/3.and.ipz<2*nprocz/3.) then
        allocate(uubufo(mx,bufsizes_yz_corn(2,INUU,ISND),bufsizes_yz_corn(1,INUU,ISND),mcom))
      else
        allocate(uubufo(mx,bufsizes_yz_corn(1,INUU,ISND),bufsizes_yz_corn(2,INUU,ISND),mcom))
      endif

      if (llast_proc_z.or.lfirst_proc_y.and.ipz>=nprocz/3.and.ipz<2*nprocz/3.) then
        allocate(lubufo(mx,bufsizes_yz_corn(2,INLU,ISND),bufsizes_yz_corn(1,INLU,ISND),mcom) )
      else
        allocate(lubufo(mx,bufsizes_yz_corn(1,INLU,ISND),bufsizes_yz_corn(2,INLU,ISND),mcom) )
      endif

      call mpibarrier

    endsubroutine yyinit
!***********************************************************************
    subroutine setup_interp_yy(gridbuf_left, gridbuf_midy, gridbuf_midz, gridbuf_right)
!
!  Initializes interpolation of data to be communicated between the grids:
!  Each proc at a grid boundary tells its neighbors at which points it needs
!  data from them. The neighbors determine the corresponding interpolation data
!  and store them in intcoeffs.
!
!  20-dec-15/MR: coded
!  12-dec-17/MR: adaptations for fixing the gap problem: ngap counts the points
!                which lie either in a y or a z gap.
!
      use General, only: yy_transform_strip, transpose_mn, itoa, reset_triangle, notanumber

      real, dimension(:,:,:) :: gridbuf_midy, gridbuf_midz, gridbuf_left, gridbuf_right
      real, dimension(:,:,:), allocatable :: thphprime_strip_y, thphprime_strip_z, tmp

      integer :: size_corn, size_neigh, lenbuf, lenred
      integer :: nok, noks, noks_all, ngap, ngap_all, nstrip_total
      logical :: stop_flag
      character(LEN=linelen) :: msg
!
      if (cyinyang_intpol_type=='bilinear') then
        iyinyang_intpol_type=BILIN
      elseif (cyinyang_intpol_type=='biquadratic') then
        iyinyang_intpol_type=BIQUAD
      elseif (cyinyang_intpol_type=='bicubic') then
        iyinyang_intpol_type=BICUB
      elseif (cyinyang_intpol_type=='biquintic') then
        iyinyang_intpol_type=BIQUIN
      elseif (cyinyang_intpol_type=='quadspline') then
        iyinyang_intpol_type=QUADSPLINE
      else
        call stop_it('setup_interp_yy: Unknown Yin-Yang interpolation type '//trim(cyinyang_intpol_type))
      endif
!
      noks=0; ngap=0
!
      if (lfirst_proc_z.or.llast_proc_z) then                ! executing proc is at z boundary
!
!  z-boundary procs: communicate only to y-boundary ones -> z strip (nghost x ...) expected.
!
        if (lcorner_yz) then
          lenbuf=len_cornbuf
          lenred=lenbuf-len_cornstrip_z
        else
          lenbuf=my
          lenred=my
        endif

        allocate(thphprime_strip_z(2,nghost,lenbuf))         ! buffer for full y or corner strip
        size_neigh=2*nghost*lenbuf                           ! buffer length for direct neighbors
        size_corn=2*nghost*lenred                            ! buffer length for corner neighbors

        if (lfirst_proc_z) then                              ! lower z boundary

          bufsizes_yz(INZL,IRCV) = lenbuf
          bufsizes_yz_corn(:,INLL,IRCV) = (/lenred,nghost/)  !(/my,nghost/)
          bufsizes_yz_corn(:,INUL,IRCV) = (/lenred,nghost/)  !(/my,nghost/)
!
! transformed coordinates into y strip
!
          if (lcutoff_corners.and.lcorner_yz) then
            if (llast_proc_y) then
              call yy_transform_strip(1,nycut,1,nghost,thphprime_strip_z(:,:,:nycut))
              call yy_transform_strip(nycut+1,my,1,nghost,thphprime_strip_z(:,:,nycut+1:my),iph_shift_=1)
              call yy_transform_strip(nycut+1,my,-nghost+1,0,thphprime_strip_z(:,:,my+1:ycornstart-1),iph_shift_=1)
              call reset_triangle(nghost,2,my,my-(nghost-2),thphprime_strip_z)   ! cutoff overlap tip
            elseif (lfirst_proc_y) then
              call yy_transform_strip(1,my-nycut,mz-len_cornstrip_z+1,mz-len_cornstrip_z+nghost, &
                                      thphprime_strip_z(:,:,:my-nycut),iph_shift_=-1)
              call yy_transform_strip(1,my-nycut,mz-len_cornstrip_z+1-nghost,mz-len_cornstrip_z, &
                                      thphprime_strip_z(:,:,my-nycut+1:2*(my-nycut)),iph_shift_=-1)
              call reset_triangle(nghost,2,1,nghost-1,thphprime_strip_z)       ! cutoff overlap tip
              call yy_transform_strip(my-nycut+1,my,1,nghost,thphprime_strip_z(:,:,2*(my-nycut)+1:ycornstart-1))
            endif
          else
            call yy_transform_strip(1,my,1,nghost,thphprime_strip_z(:,:,:my))
          endif
!
!  At yz corner: create a "cornerstrip" by adding truncated z strip
!
          if (lfirst_proc_y) then
            call yy_transform_strip(1,nghost,mz-len_cornstrip_z+1,mz,thphprime_strip_z(:,:,ycornstart:))   ! ll corner
          elseif (llast_proc_y) then
            call yy_transform_strip(m2+1,my,mz-len_cornstrip_z+1,mz,thphprime_strip_z(:,:,ycornstart:))    ! ul corner
          endif
if (.false.) then
!if (.not.lyang.and.(lfirst_proc_y.or.llast_proc_y)) then
!if (lfirst_proc_y) print*, 'lower left'
!if (llast_proc_y) print*, 'lower right'
!print*, 'iproc, zlneigh, llcorn, ulcorn=', iproc, zlneigh-27, llcorn-27, ulcorn-27
!print*, 'thphprime_strip_z'
!if (lfirst_proc_y) print'(2(f10.4,1x))', thphprime_strip_z
if (lfirst_proc_y) write(100,'(2(f10.4,1x))') thphprime_strip_z
if (llast_proc_y) write(102,'(2(f10.4,1x))') thphprime_strip_z
print*, '-----------------'
endif
!          size(thphprime_strip_z,2), size(thphprime_strip_z,3), sum(thphprime_strip_z(:,:,:57))

          call MPI_IRECV(gridbuf_midy,2*nghost*yy_buflens(MID),mpi_precision, &         ! receive strip from ~
                         zlneigh,zlneigh,MPI_COMM_PENCIL,irecv_rq_fromlowz,mpierr)
                                !MPI_ANY_TAG
          call MPI_ISEND(thphprime_strip_z,size_neigh,mpi_precision, &                   ! send strip to direct neighbor
                         zlneigh,iproc_world,MPI_COMM_PENCIL,isend_rq_tolowz,mpierr)
                                !tolowz
          if (llcorn>=0) then
!
            call MPI_IRECV(gridbuf_right,2*nghost*yy_buflens(RIGHT),mpi_precision, &    ! receive strip from ~
                           llcorn,llcorn,MPI_COMM_PENCIL,irecv_rq_FRll,mpierr)
                                  !MPI_ANY_TAG
!print*, 'proc ', iproc_world, 'receives ', nghost*yy_buflens(RIGHT), 'from ', llcorn
!if (.not.lyang) write(iproc+200,*) iproc, 'sends', size_corn/(2*nghost), 'to ', llcorn
            call MPI_ISEND(thphprime_strip_z(:,:,:lenred),size_corn,mpi_precision, &
!  send strip (without corner part) to right corner neighbor
                           llcorn,iproc_world,MPI_COMM_PENCIL,isend_rq_TOll,mpierr)
                                  !TOll
          endif

          if (ulcorn>=0) then
!
            call MPI_IRECV(gridbuf_left,2*nghost*yy_buflens(LEFT),mpi_precision, &      ! receive strip from ~
                           ulcorn,ulcorn,MPI_COMM_PENCIL,irecv_rq_FRul,mpierr)
                                  !MPI_ANY_TAG
!print*, 'proc ', iproc_world, 'receives ', nghost*yy_buflens(LEFT), 'from ', ulcorn
            call MPI_ISEND(thphprime_strip_z(:,:,:lenred),size_corn,mpi_precision, &
!  send strip (without corner part) to left corner neighbor
                           ulcorn,iproc_world,MPI_COMM_PENCIL,isend_rq_TOul,mpierr)
                                  !TOul
          endif
        endif

        if (llast_proc_z) then                                                     ! upper z boundary

          bufsizes_yz(INZU,IRCV) = lenbuf
          bufsizes_yz_corn(:,INLU,IRCV) = (/lenred,nghost/)  !(/my,nghost/)
          bufsizes_yz_corn(:,INUU,IRCV) = (/lenred,nghost/)  !(/my,nghost/)
!
! transformed coordinates into y strip
!
          if (lcutoff_corners.and.lcorner_yz) then
            if (llast_proc_y) then
              call yy_transform_strip(1,nycut,n2+1,mz,thphprime_strip_z(:,:,:nycut))
              call yy_transform_strip(nycut+1,my,n2+1,mz,thphprime_strip_z(:,:,nycut+1:my),iph_shift_=-1)
              call yy_transform_strip(nycut+1,my,mz+1,mz+nghost,thphprime_strip_z(:,:,my+1:ycornstart-1),iph_shift_=-1)
              call reset_triangle(1,nghost-1,my,my-(nghost-2),thphprime_strip_z)   ! reset overlap tip
            elseif (lfirst_proc_y) then
              call yy_transform_strip(1,my-nycut,len_cornstrip_z-(nghost-1),len_cornstrip_z, &
                                      thphprime_strip_z(:,:,:my-nycut),iph_shift_=1)
              call yy_transform_strip(1,my-nycut,len_cornstrip_z+1,len_cornstrip_z+nghost, &
                                      thphprime_strip_z(:,:,my-nycut+1:2*(my-nycut)),iph_shift_=1)
              call yy_transform_strip(my-nycut+1,my,n2+1,mz,thphprime_strip_z(:,:,2*(my-nycut)+1:ycornstart-1))
              call reset_triangle(1,nghost-1,1,nghost-1,thphprime_strip_z)       ! reset overlap tip
            endif
          else
            call yy_transform_strip(1,my,n2+1,mz,thphprime_strip_z(:,:,:my))
          endif
!
!  At yz corner: create a "cornerstrip" by adding truncated z strip
!
          if (lfirst_proc_y) then
            call yy_transform_strip(1,nghost,1,len_cornstrip_z,thphprime_strip_z(:,:,ycornstart:))     ! lu corner
          elseif (llast_proc_y) then
            call yy_transform_strip(m2+1,my,1,len_cornstrip_z,thphprime_strip_z(:,:,ycornstart:))      ! uu corner
          endif
if (.false.) then
!if (.not.lyang.and.(llast_proc_y.or.lfirst_proc_y)) then
!if (lfirst_proc_y) print*, 'upper left'
!if (llast_proc_y) print*, 'upper right'
!print*, 'thphprime_strip_z'
!print'(2(f10.4,1x))', thphprime_strip_z
if (lfirst_proc_y) write(24,'(2(f10.4,1x))') thphprime_strip_z
if (llast_proc_y) write(26,'(2(f10.4,1x))') thphprime_strip_z
!print*, '-----------------'
endif

!if (lcorner_yz) print*, 'proc ',iproc_world,', sends thphprime_strip to ', zuneigh, &
!size(thphprime_strip_z,2), size(thphprime_strip_z,3), sum(thphprime_strip_z(:,:,:57))
          call MPI_IRECV(gridbuf_midy,2*nghost*yy_buflens(MID),mpi_precision, &         ! receive strip from ~
                         zuneigh,zuneigh,MPI_COMM_PENCIL,irecv_rq_fromuppz,mpierr)
                                !MPI_ANY_TAG
          call MPI_ISEND(thphprime_strip_z,size_neigh,mpi_precision, &                   ! send strip to direct neighbor
                         zuneigh,iproc_world,MPI_COMM_PENCIL,isend_rq_touppz,mpierr)
                                !touppz
          if (lucorn>=0) then

            call MPI_IRECV(gridbuf_left,2*nghost*yy_buflens(LEFT),mpi_precision, &      ! receive strip from ~
                           lucorn,lucorn,MPI_COMM_PENCIL,irecv_rq_FRlu,mpierr)
                                  !MPI_ANY_TAG
!print*, 'proc ', iproc_world, 'receives ', nghost*yy_buflens(LEFT), 'from ', lucorn
            call MPI_ISEND(thphprime_strip_z(:,:,:lenred),size_corn,mpi_precision, &
!  send strip (without corner part) to left corner neighbor
                           lucorn,iproc_world,MPI_COMM_PENCIL,isend_rq_TOlu,mpierr)
                                  !TOlu
          endif

          if (uucorn>=0) then
!
            call MPI_IRECV(gridbuf_right,2*nghost*yy_buflens(RIGHT),mpi_precision, &    ! receive strip from ~
                           uucorn,uucorn,MPI_COMM_PENCIL,irecv_rq_FRuu,mpierr)
                                  !MPI_ANY_TAG
!print*, 'proc ', iproc_world, 'receives ', nghost*yy_buflens(RIGHT), 'from ', uucorn
            call MPI_ISEND(thphprime_strip_z(:,:,:lenred),size_corn,mpi_precision, &              ! send strip to right corner neighbor
                           uucorn,iproc_world,MPI_COMM_PENCIL,isend_rq_TOuu,mpierr)
                                  !TOuu
          endif
        endif
      endif

      if (lfirst_proc_y.or.llast_proc_y) then
!
!  y-boundary processors: communicate to y- or z-boundary ones -> z strip (nghost x ...) or y strip (... x nghost) expected.
!
        if (lcorner_yz) then
          lenbuf=len_cornbuf
          lenred=lenbuf-len_cornstrip_y
        else
          lenbuf=mz
          lenred=mz
        endif

        allocate(thphprime_strip_y(2,nghost,lenbuf))                               ! for communication to y boundary

        size_neigh=2*nghost*lenbuf
        size_corn=2*nghost*lenred

        if (lfirst_proc_y) then                                                    ! lower y boundary

          bufsizes_yz(INYL,IRCV) = lenbuf
          bufsizes_yz_corn(:,INLL,IRCV) = (/nghost,lenred/)   !(/nghost,mz/)
          bufsizes_yz_corn(:,INLU,IRCV) = (/nghost,lenred/)   !(/nghost,mz/)

          if (lcutoff_corners.and.lcorner_yz) then
            if (llast_proc_z) then
              call yy_transform_strip(1,nghost,1,nzcut,thphprime_strip_y(:,:,:nzcut))
              call yy_transform_strip(1,nghost,nzcut+1,mz, &
                                      thphprime_strip_y(:,:,nzcut+1:mz),ith_shift_=1)
              call yy_transform_strip(-nghost+1,0,nzcut+1,mz, &
                                      thphprime_strip_y(:,:,mz+1:zcornstart-1),ith_shift_=1)
              call reset_triangle(nghost,2,mz,mz-(nghost-2),thphprime_strip_y)   ! reset overlap tip
            elseif (lfirst_proc_z) then
              call yy_transform_strip(my-len_cornstrip_y+1,my-len_cornstrip_y+nghost,1,mz-nzcut, &
                                      thphprime_strip_y(:,:,:mz-nzcut),ith_shift_=-1)
              call yy_transform_strip(my-len_cornstrip_y+1-nghost,my-len_cornstrip_y,1,mz-nzcut, &
                                      thphprime_strip_y(:,:,mz-nzcut+1:2*(mz-nzcut)),ith_shift_=-1)
              call yy_transform_strip(1,nghost,mz-nzcut+1,mz,thphprime_strip_y(:,:,2*(mz-nzcut)+1:zcornstart-1))
              call reset_triangle(nghost,2,1,nghost-1,thphprime_strip_y)       ! reset overlap tip
            endif
          else
            call yy_transform_strip(1,nghost,1,mz,thphprime_strip_y(:,:,:mz))          ! full z strip
          endif
!if (lyang.and.lcorner_yz) print*, 'y-ll/lu:', sum(thphprime_strip_y(:,:,:mz))
!
!  At yz corner: create a "cornerstrip"
!
          if (lfirst_proc_z) then
!print*, 'left lower completed: len_cornstrip_y=', len_cornstrip_y
            call yy_transform_strip(my-len_cornstrip_y+1,my,1,nghost,thphprime_strip_y(:,:,zcornstart:))     ! left lower corner completed
          elseif (llast_proc_z) then
!print*, 'left upper completed: len_cornstrip_y=', len_cornstrip_y
            call yy_transform_strip(my-len_cornstrip_y+1,my,n2+1,mz,thphprime_strip_y(:,:,zcornstart:))    ! left upper corner completed
          endif
if (.false.) then
!if (.not.lyang.and.(lfirst_proc_z.or.llast_proc_z)) then
!print*, 'iproc, ylneigh, llcorn, lucorn=', iproc, ylneigh-27, llcorn-27, lucorn-27
!if (lfirst_proc_z) print*, 'left lower'
!if (llast_proc_z) print*, 'left upper'   !: ylneigh, llcorn, lucorn=', ylneigh, llcorn, lucorn
!print*, 'thphprime_strip_y'
if (lfirst_proc_z) write(100,'(2(f10.4,1x))') thphprime_strip_y
if (llast_proc_z) write(24,'(2(f10.4,1x))') thphprime_strip_y
!print*, '-----------------'
endif

          if (ipz>=nprocz/3.and.ipz<2*nprocz/3) then
!
!  Transposition of thphprime_strip_y, can perhaps be avoided.
!
            allocate(tmp(2,nghost,lenbuf))
            tmp=thphprime_strip_y
            deallocate(thphprime_strip_y)
            allocate(thphprime_strip_y(2,lenbuf,nghost))
            call transpose_mn(tmp,thphprime_strip_y)
            deallocate(tmp)
          endif
!if (lcorner_yz) print*, 'proc ',iproc_world,', sends thphprime_strip to ', yuneigh, &
! size(thphprime_strip_y,2), size(thphprime_strip_y,3), sum(thphprime_strip_z(:,:,:57))

!print*, 'proc ', iproc_world, 'receives ', nghost*yy_buflens(MID), 'from ', yuneigh

!print*, 'proc ', iproc_world, 'receives ', nghost*yy_buflens(MID), 'from ', ylneigh
          call MPI_IRECV(gridbuf_midz,2*nghost*yy_buflens(MID),mpi_precision, &         ! receive strip from ~
                         ylneigh,ylneigh,MPI_COMM_PENCIL,irecv_rq_fromlowy,mpierr)
                                !MPI_ANY_TAG
          call MPI_ISEND(thphprime_strip_y,size_neigh,mpi_precision, &                   ! send strip to direct neighbor
                         ylneigh,iproc_world,MPI_COMM_PENCIL,isend_rq_tolowy,mpierr)
                                !tolowy
          if (llcorn>=0) then
!
            call MPI_IRECV(gridbuf_left,2*nghost*yy_buflens(LEFT),mpi_precision, &      ! receive strip from ~
                           llcorn,llcorn,MPI_COMM_PENCIL,irecv_rq_FRll,mpierr)
                                  !MPI_ANY_TAG
!print*, 'proc ', iproc_world, 'receives ', nghost*yy_buflens(LEFT), 'from ', llcorn
            call MPI_ISEND(thphprime_strip_y,size_corn,mpi_precision, &                       ! send strip to left corner neighbor
                           llcorn,iproc_world,MPI_COMM_PENCIL,isend_rq_TOll,mpierr)
                                  !TOll
          endif

          if (lucorn>=0) then
!
            call MPI_IRECV(gridbuf_right,2*nghost*yy_buflens(RIGHT),mpi_precision, &    ! receive strip from ~
                           lucorn,lucorn,MPI_COMM_PENCIL,irecv_rq_FRlu,mpierr)
!print*, 'proc ', iproc_world, 'receives ', nghost*yy_buflens(RIGHT), 'from ', lucorn
            call MPI_ISEND(thphprime_strip_y,size_corn,mpi_precision, &                       ! send strip to right corner neighbor
                           lucorn,iproc_world,MPI_COMM_PENCIL,isend_rq_TOlu,mpierr)
                                  !TOlu
                                  !MPI_ANY_TAG
          endif
        endif

        if (llast_proc_y) then

          bufsizes_yz(INYU,IRCV) = lenbuf
          bufsizes_yz_corn(:,INUL,IRCV) = (/nghost,lenred/)   !(/nghost,mz/)
          bufsizes_yz_corn(:,INUU,IRCV) = (/nghost,lenred/)   !(/nghost,mz/)

          if (lcutoff_corners.and.lcorner_yz) then
            if (llast_proc_z) then
              call yy_transform_strip(m2+1,my,1,nzcut,thphprime_strip_y(:,:,:nzcut))
              call yy_transform_strip(m2+1,my,nzcut+1,mz,thphprime_strip_y(:,:,nzcut+1:mz), &
                                      ith_shift_=-1)
              call yy_transform_strip(my+1,my+nghost,nzcut+1,mz,thphprime_strip_y(:,:,mz+1:zcornstart-1), &
                                      ith_shift_=-1)
              call reset_triangle(1,nghost-1,mz,mz-(nghost-2),thphprime_strip_y)   ! reset overlap tip
            elseif (lfirst_proc_z) then
              call yy_transform_strip(len_cornstrip_y-(nghost-1),len_cornstrip_y,1, &
                                      mz-nzcut,thphprime_strip_y(:,:,:mz-nzcut),ith_shift_=1)
              call yy_transform_strip(len_cornstrip_y+1,len_cornstrip_y+nghost,1, &
                                      mz-nzcut,thphprime_strip_y(:,:,mz-nzcut+1:2*(mz-nzcut)),ith_shift_=1)
              call yy_transform_strip(m2+1,my,mz-nzcut+1,mz,thphprime_strip_y(:,:,2*(mz-nzcut)+1:zcornstart-1))
              call reset_triangle(1,nghost-1,1,nghost-1,thphprime_strip_y)       ! reset overlap tip
            endif
          else
            call yy_transform_strip(m2+1,my,1,mz,thphprime_strip_y(:,:,:mz))         ! full z ghost-strip
          endif
!
!  Create a corner strip.
!
          if (lfirst_proc_z) then
            call yy_transform_strip(1,len_cornstrip_y,1,nghost,thphprime_strip_y(:,:,zcornstart:))   ! right lower corner completed
          elseif (llast_proc_z) then
            call yy_transform_strip(1,len_cornstrip_y,n2+1,mz,thphprime_strip_y(:,:,zcornstart:))    ! right upper corner completed
          endif
if (.false.) then
!if (.not.lyang.and.(lfirst_proc_z.or.llast_proc_z)) then
!if (lfirst_proc_z) print*, 'right lower'
!if (llast_proc_z) print*, 'right upper'  !: yuneigh, ulcorn, uucorn=', yuneigh, ulcorn, uucorn
!print*, 'thphprime_strip_y'
if (lfirst_proc_z) write(102,'(2(f10.4,1x))') thphprime_strip_y
if (llast_proc_z) write(26,'(2(f10.4,1x))') thphprime_strip_y
!print*, '-----------------'
endif

          if (ipz>=nprocz/3.and.ipz<2*nprocz/3) then
!
!  Transposition of thphprime_strip_y, can perhaps be avoided.
!
            allocate(tmp(2,nghost,lenbuf))
            tmp=thphprime_strip_y
            deallocate(thphprime_strip_y)
            allocate(thphprime_strip_y(2,lenbuf,nghost))
            call transpose_mn(tmp,thphprime_strip_y)
            deallocate(tmp)
          endif
!print*, 'at IRECV: proc ,',iproc_world,',  from ', yuneigh, size(gridbuf_midz,1), size(gridbuf_midz,2), size(gridbuf_midz,3)
          call MPI_IRECV(gridbuf_midz,2*nghost*yy_buflens(MID),mpi_precision, &          ! receive strip from ~
                         yuneigh,yuneigh,MPI_COMM_PENCIL,irecv_rq_fromuppy,mpierr)
                                !MPI_ANY_TAG

!print*, 'proc ',iproc_world,', sends thphprime_strip to ', yuneigh, &
!size(thphprime_strip_y,1), size(thphprime_strip_y,2), size(thphprime_strip_y,3)

          call MPI_ISEND(thphprime_strip_y,size_neigh,mpi_precision, &                    ! send strip to direct neighbor
                         yuneigh,iproc_world,MPI_COMM_PENCIL,isend_rq_touppy,mpierr)
                                !touppy
          if (ulcorn>=0) then
! send strip to right corner neighbor

            call MPI_IRECV(gridbuf_right,2*nghost*yy_buflens(RIGHT),mpi_precision, &     ! receive strip from ~
                           ulcorn,ulcorn,MPI_COMM_PENCIL,irecv_rq_FRul,mpierr)
                                  !MPI_ANY_TAG
!print*, 'proc ', iproc_world, 'receives ', nghost*yy_buflens(RIGHT), 'from ', ulcorn
            call MPI_ISEND(thphprime_strip_y,size_corn,mpi_precision, &
                           ulcorn,iproc_world,MPI_COMM_PENCIL,isend_rq_TOul,mpierr)
                                  !TOul
          endif
          if (uucorn>=0) then
            call MPI_IRECV(gridbuf_left,2*nghost*yy_buflens(LEFT),mpi_precision, &       ! receive strip from ~
                           uucorn,uucorn,MPI_COMM_PENCIL,irecv_rq_FRuu,mpierr)
                                  !MPI_ANY_TAG
!print*, 'proc ', iproc_world, 'receives ', nghost*yy_buflens(LEFT), 'from ', uucorn
            call MPI_ISEND(thphprime_strip_y,size_corn,mpi_precision, &                        ! send strip to left corner neighbor
                           uucorn,iproc_world,MPI_COMM_PENCIL,isend_rq_TOuu,mpierr)
                                  !TOuu
          endif
        endif

      endif
!
!  Now finalize all communications and determine the interpolation elements for the transmitted coordinates.
!
      if (lfirst_proc_y) then

        call MPI_WAIT(irecv_rq_fromlowy,irecv_stat_fl,mpierr)
        nok=prep_interp(gridbuf_midz,intcoeffs(MIDZ),iyinyang_intpol_type,ngap); noks=noks+nok
        call MPI_WAIT(isend_rq_tolowy,isend_stat_tl,mpierr)
!if (lyang.and.iproc==15) print*, 'fromlowy: iproc, nok=', iproc, nok, ngap
!if (lyang.and.iproc==18) print*, 'fromlowy: iproc, nok=', iproc, nok, ngap
!if (lyang.and.iproc==6) print*, 'fromlowy: iproc, nok=', iproc, nok, ngap
!if (lyang.and.iproc==9) print*, 'fromlowy: iproc, nok=', iproc, nok, ngap
        if (llcorn>=0) then
          call MPI_WAIT(irecv_rq_FRll,irecv_stat_Fll,mpierr)
          nok=prep_interp(gridbuf_left,intcoeffs(LEFT),iyinyang_intpol_type,ngap); noks=noks+nok
!if (lyang.and.iproc==12) print*, 'fromuppleft: iproc, nok=', iproc, nok, ngap
!if (lyang.and.iproc==3) print*, 'fromuppleft: iproc, nok=', iproc, nok, ngap
!if (iproc_world==6) print*,'ll:', iproc_world,nok
!,maxval(abs(intcoeffs(LEFT)%coeffs)),maxval(intcoeffs(LEFT)%inds),minval(intcoeffs(LEFT)%inds)
          call MPI_WAIT(isend_rq_TOll,isend_stat_Tll,mpierr)
        endif

        if (lucorn>=0) then
          call MPI_WAIT(irecv_rq_FRlu,irecv_stat_Flu,mpierr)
          nok=prep_interp(gridbuf_right,intcoeffs(RIGHT),iyinyang_intpol_type,ngap); noks=noks+nok
!if (lyang.and.iproc==21) print*, 'fromlowleft: iproc, nok=', iproc, nok, ngap
!if (lyang.and.iproc==12) print*, 'fromlowleft: iproc, nok=', iproc, nok, ngap
!maxval(abs(intcoeffs(RIGHT)%coeffs)),maxval(intcoeffs(RIGHT)%inds),minval(intcoeffs(RIGHT)%inds)
          call MPI_WAIT(isend_rq_TOlu,isend_stat_Tlu,mpierr)
        endif

      endif

      if (llast_proc_y) then
!
        call MPI_WAIT(irecv_rq_fromuppy,irecv_stat_fu,mpierr)
        nok=prep_interp(gridbuf_midz,intcoeffs(MIDZ),iyinyang_intpol_type,ngap); noks=noks+nok
! print*,'uppy:',iproc,iproc_world,nok,maxval(abs(intcoeffs(MIDZ)%coeffs)),maxval(intcoeffs(MIDZ)%inds), &
! minval(intcoeffs(MIDZ)%inds)
        call MPI_WAIT(isend_rq_touppy,isend_stat_tu,mpierr)
!if (lyang.and.iproc==17) print*, 'fromuppy: iproc, nok=', iproc, nok, ngap
!if (lyang.and.iproc==20) print*, 'fromuppy: iproc, nok=', iproc, nok, ngap
!if (lyang.and.iproc==11) print*, 'fromuppy: iproc, nok=', iproc, nok, ngap
!if (lyang.and.iproc==8) print*, 'fromuppy: iproc, nok=', iproc, nok, ngap

        if (ulcorn>=0) then
          call MPI_WAIT(irecv_rq_FRul,irecv_stat_Ful,mpierr)
          nok=prep_interp(gridbuf_right,intcoeffs(RIGHT),iyinyang_intpol_type,ngap); noks=noks+nok
!if (lyang.and.iproc==23) print*, 'fromlowright: iproc, nok=', iproc, nok, ngap
!if (lyang.and.iproc==14) print*, 'fromlowright: iproc, nok=', iproc, nok, ngap
! print*,'ul:',iproc,iproc_world,nok,maxval(abs(intcoeffs(RIGHT)%coeffs)),maxval(intcoeffs(RIGHT)%inds), &
! minval(intcoeffs(RIGHT)%inds)
          call MPI_WAIT(isend_rq_TOul,isend_stat_Tul,mpierr)
        endif

        if (uucorn>=0) then
          call MPI_WAIT(irecv_rq_FRuu,irecv_stat_Fuu,mpierr)
          nok=prep_interp(gridbuf_left,intcoeffs(LEFT),iyinyang_intpol_type,ngap); noks=noks+nok
!if (lyang.and.iproc==14) print*, 'fromuppright: iproc, nok=', iproc, nok, ngap
!if (lyang.and.iproc==5) print*, 'fromuppright: iproc, nok=', iproc, nok, ngap
! print*,'uu:',iproc,iproc_world,nok,maxval(abs(intcoeffs(LEFT)%coeffs)),maxval(intcoeffs(LEFT)%inds), &
! minval(intcoeffs(LEFT)%inds)
          call MPI_WAIT(isend_rq_TOuu,isend_stat_Tuu,mpierr)
        endif

      endif

      if (lfirst_proc_z) then
!
        call MPI_WAIT(irecv_rq_fromlowz,irecv_stat_fl,mpierr)
        nok=prep_interp(gridbuf_midy,intcoeffs(MIDY),iyinyang_intpol_type,ngap); noks=noks+nok
! print*,'lowz:',iproc,iproc_world,nok,maxval(abs(intcoeffs(MIDY)%coeffs)),maxval(intcoeffs(MIDY)%inds),&
!  minval(intcoeffs(MIDY)%inds)
if (.not.lyang.and.notanumber(intcoeffs(MIDY)%coeffs)) print*, 'iproc=', iproc

        call MPI_WAIT(isend_rq_tolowz,isend_stat_tl,mpierr)
!
        if (llcorn>=0) then
          call MPI_WAIT(irecv_rq_FRll,irecv_stat_Fll,mpierr)
          nok=prep_interp(gridbuf_right,intcoeffs(RIGHT),iyinyang_intpol_type,ngap); noks=noks+nok
! print*,'ll:',iproc,iproc_world,nok,maxval(abs(intcoeffs(RIGHT)%coeffs)),maxval(intcoeffs(RIGHT)%inds), &
! minval(intcoeffs(RIGHT)%inds)
if (.not.lyang.and.notanumber(intcoeffs(RIGHT)%coeffs)) print*, 'iproc=', iproc
          call MPI_WAIT(isend_rq_TOll,isend_stat_Tll,mpierr)
        endif

        if (ulcorn>=0) then
          call MPI_WAIT(irecv_rq_FRul,irecv_stat_Ful,mpierr)
          nok=prep_interp(gridbuf_left,intcoeffs(LEFT),iyinyang_intpol_type,ngap); noks=noks+nok
! print*,'ul:',iproc,iproc_world,nok,maxval(abs(intcoeffs(LEFT)%coeffs)),maxval(intcoeffs(LEFT)%inds), &
! minval(intcoeffs(LEFT)%inds)
if (.not.lyang.and.notanumber(intcoeffs(LEFT)%coeffs)) print*, 'iproc=', iproc
          call MPI_WAIT(isend_rq_TOul,isend_stat_Tul,mpierr)
        endif

      endif

      if (llast_proc_z) then
!
        call MPI_WAIT(irecv_rq_fromuppz,irecv_stat_fu,mpierr)
        nok=prep_interp(gridbuf_midy,intcoeffs(MIDY),iyinyang_intpol_type,ngap); noks=noks+nok
!  print*, 'proc ,',iproc_world,',  from ', zuneigh, size(gridbuf_midy,1), &
!  size(gridbuf_midy,2), size(gridbuf_midy,3), size(intcoeffs(MIDY)%inds,1), size(intcoeffs(MIDY)%inds,2)
!  print'(22(f8.5,1x))', gridbuf_midy(1,:,:)
!  print'(22(i3,1x))', intcoeffs(MIDY)%inds(:,:,1)
!  print*,'upz:',iproc,iproc_world,nok,maxval(abs(intcoeffs(MIDY)%coeffs)),maxval(intcoeffs(MIDY)%inds), &
!  minval(intcoeffs(MIDY)%inds)
if (.not.lyang.and.notanumber(intcoeffs(MIDY)%coeffs)) print*, 'iproc=', iproc
        call MPI_WAIT(isend_rq_touppz,isend_stat_tu,mpierr)
!
        if (lucorn>=0) then
          call MPI_WAIT(irecv_rq_FRlu,irecv_stat_Flu,mpierr)
          nok=prep_interp(gridbuf_left,intcoeffs(LEFT),iyinyang_intpol_type,ngap); noks=noks+nok
!  print*,'lu:',iproc,iproc_world,nok,maxval(abs(intcoeffs(LEFT)%coeffs)),maxval(intcoeffs(LEFT)%inds), &
!  minval(intcoeffs(LEFT)%inds)
if (.not.lyang.and.notanumber(intcoeffs(LEFT)%coeffs)) print*, 'iproc=', iproc
          call MPI_WAIT(isend_rq_TOlu,isend_stat_Tlu,mpierr)
        endif

        if (uucorn>=0) then
          call MPI_WAIT(irecv_rq_FRuu,irecv_stat_Fuu,mpierr)
          nok=prep_interp(gridbuf_right,intcoeffs(RIGHT),iyinyang_intpol_type,ngap); noks=noks+nok
!  print*,'uu:',iproc,iproc_world,nok,maxval(abs(intcoeffs(RIGHT)%coeffs)),maxval(intcoeffs(RIGHT)%inds), &
!  minval(intcoeffs(RIGHT)%inds)
if (.not.lyang.and.notanumber(intcoeffs(RIGHT)%coeffs)) print*, 'iproc=', iproc
          call MPI_WAIT(isend_rq_TOuu,isend_stat_Tuu,mpierr)
        endif

      endif
!
!  Over each of the two grids, sum up numbers of points (nok) for which interpolation data
!  could be generated. Sum (noks_all) must agree with total number of ghost zone points
!  along boundary of each grid (nstrip_total).
!
      call mpireduce_sum_int(noks, noks_all,comm=MPI_COMM_YZPLANE)
      call mpireduce_sum_int(ngap, ngap_all,comm=MPI_COMM_YZPLANE)

      stop_flag=.false.; msg=''
      if (iproc==root) then
        if (lcutoff_corners) then
          nstrip_total= 2*nghost*((nprocz-1)*mz + (nprocy-1)*my + len_cornstrip_z + len_cornstrip_y)  &
                       -4*(nghost-1)*nghost/2 &                       !??
                       +4*nghost*(my-nycut)-4*(nghost-1)*nghost/2     ! 12 spec for nghost=3!!!
        else
          nstrip_total=2*nghost*(nprocz*mz + nprocy*my - 2*nghost)
        endif
        noks_all=noks_all - ngap_all/4

print*, 'noks_all,ngap_all,nstrip_total=', noks_all,ngap_all,nstrip_total
        if (noks_all/=nstrip_total) then
          msg='setup_interp_yy: '//trim(cyinyang)//' grid: number of caught points '// &
               trim(itoa(noks_all))//' unequal goal '// trim(itoa(nstrip_total)) &
               //'. Reduce dang in start.f90.'
          if (lcutoff_corners) msg=trim(msg)//' Or increase nycut and nzcut.'
          stop_flag=.true.
        endif
      endif
      call stop_it_if_any(stop_flag,msg)
!call  mpibarrier
!call  mpifinalize
!stop
    endsubroutine setup_interp_yy
!***********************************************************************
    subroutine initiate_isendrcv_bdry(f,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values. Called in the beginning of pde.
!  Does not wait for the receives to finish (done in finalize_isendrcv_bdry)
!  leftneigh and rightneigh are initialized by mpicomm_init.
!
!  21-may-02/axel: communication of corners added
!  11-aug-07/axel: communication in the x-direction added
!  22-oct-15/fred: communication for spherical polar coords across poles
!                  added yz corner communication
!  20-dec-15/MR: modified for Yin-Yang grid
!
      use General, only: notanumber

      real, dimension(:,:,:,:), intent(inout):: f
      integer, optional,        intent(in)   :: ivar1_opt, ivar2_opt
!
      integer :: tolowyr,touppyr,tolowys,touppys,tolowzr,touppzr,tolowzs,touppzs ! msg. tags placeholders
      integer :: TOllr,TOulr,TOuur,TOlur,TOlls,TOuls,TOuus,TOlus                 ! placeholder tags
      integer :: ivar1, ivar2, nbufy, nbufz, nbufyz, mxl, comm, bufact, dir,mm,j
!
      ivar1=1; ivar2=min(mcom,size(f,4))

      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
      if (ivar2==0) return
!
      mxl=size(f,1)
!
!  Set communication across x-planes.
!
      if (nprocx>1) call isendrcv_bdry_x(f,ivar1_opt,ivar2_opt)
!
!  Set communication across y-planes.
!  Standard message tags from defaults for surfaces and yz-corners.
!
      tolowyr=tolowy; tolowzr=tolowz; TOlls=TOll; TOllr=TOll
      tolowys=tolowy; tolowzs=tolowz; TOlus=TOlu; TOlur=TOlu
      touppyr=touppy; touppzr=touppz; TOuus=TOuu; TOuur=TOuu
      touppys=touppy; touppzs=touppz; TOuls=TOul; TOulr=TOul
!
!  For spherical polar boundary condition across the pole set up edge and corner
!  neighbours. Assumes the domain is binary communication between z processors.
!  NB nprocz=2*n, n>=1, comms across y-plane parallel in z!
!
      if (lcommunicate_y) then
        poleneigh = find_proc(ipx,     ipy,ipz  +nprocz/2)
        pnbcrn    = find_proc(ipx,       0,ipz-1+nprocz/2)
        pnfcrn    = find_proc(ipx,       0,ipz+1+nprocz/2)
        psfcrn    = find_proc(ipx,nprocy-1,ipz+1+nprocz/2)
        psbcrn    = find_proc(ipx,nprocy-1,ipz-1+nprocz/2)
        !poleneigh = modulo(ipz  +nprocz/2,nprocz)*nprocxy+       ipy*nprocx+ipx
        !pnbcrn    = modulo(ipz-1+nprocz/2,nprocz)*nprocxy+         0*nprocx+ipx !N rev
        !pnfcrn    = modulo(ipz+1+nprocz/2,nprocz)*nprocxy+         0*nprocx+ipx !N fwd
        !psfcrn    = modulo(ipz+1+nprocz/2,nprocz)*nprocxy+(nprocy-1)*nprocx+ipx
        !psbcrn    = modulo(ipz-1+nprocz/2,nprocz)*nprocxy+(nprocy-1)*nprocx+ipx
      endif
!
!  Allocate and send/receive buffers across y-planes
!
      bufact=mxl*nghost*(ivar2-ivar1+1)

      if (nprocy>1 .or. lcommunicate_y .or. lyinyang) then
!
!  Internal, periodic and Yin-Yang-exchange y-plane buffers.
!
        if (lyinyang.and.lfirst_proc_y) then
!
!  Interpolate variables ivar1 ... ivar2 for lower y neighbor in other grid.
!  Result in lbufyo.
!
          call interpolate_yy(f,ivar1,ivar2,lbufyo,MIDZ,iyinyang_intpol_type) !
if (notanumber(lbufyo)) print*, 'lbufyo: iproc=', iproc, iproc_world
          comm=MPI_COMM_PENCIL; touppyr=ylneigh; tolowys=iproc_world   !touppyr=MPI_ANY_TAG
        else
          if (lcommunicate_y.and.lfirst_proc_y) then                        !N-pole
            lbufyo(:,:,:,ivar1:ivar2)=f(:,m1i:m1:-1,n1:n2,ivar1:ivar2)
            tolowys=shiftn; touppyr=shiftn; ylneigh=poleneigh
          else
            lbufyo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n1:n2,ivar1:ivar2)         !(lower y-zone)
          endif
          comm=MPI_COMM_GRID
        endif

        nbufy=bufact*bufsizes_yz(INYL,IRCV)
!if (ldiagnos.and.lfirst.and.iproc_world==0) print*, iproc_world, ' receives', bufsizes_yz(INYL,IRCV), ' from', ylneigh
        call MPI_IRECV(lbufyi(:,:,:,ivar1:ivar2),nbufy,mpi_precision, &
                       ylneigh,touppyr,comm,irecv_rq_fromlowy,mpierr)

        nbufy=bufact*bufsizes_yz(INYL,ISND)
        call MPI_ISEND(lbufyo(:,:,:,ivar1:ivar2),nbufy,mpi_precision, &
                       ylneigh,tolowys,comm,isend_rq_tolowy,mpierr)

        if (lyinyang.and.llast_proc_y) then
!
!  Interpolate variables ivar1 ... ivar2 for upper y neighbor in other grid.
!  Result in ubufyo.
!
          call interpolate_yy(f,ivar1,ivar2,ubufyo,MIDZ,iyinyang_intpol_type)
if (notanumber(ubufyo)) print*, 'ubufyo: iproc=', iproc, iproc_world
          comm=MPI_COMM_PENCIL; tolowyr=yuneigh; touppys=iproc_world  !tolowyr=MPI_ANY_TAG
        else
          if (lcommunicate_y.and.llast_proc_y) then                         !S-pole
!
!  North/south-pole y-plane buffers and message tags. NB N:0/S:pi radians.
!  Swap N(S) lower(upper) buffers between pi-shifted iprocz same y-plane.
!
            ubufyo(:,:,:,ivar1:ivar2)=f(:,m2:m2i:-1,n1:n2,ivar1:ivar2)
            tolowyr=shifts; touppys=shifts; yuneigh=poleneigh
          else
            ubufyo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n1:n2,ivar1:ivar2)         !(upper y-zone)
          endif
          comm=MPI_COMM_GRID
        endif
!
        nbufy=bufact*bufsizes_yz(INYU,IRCV)
        call MPI_IRECV(ubufyi(:,:,:,ivar1:ivar2),nbufy,mpi_precision, &
                       yuneigh,tolowyr,comm,irecv_rq_fromuppy,mpierr)

        nbufy=bufact*bufsizes_yz(INYU,ISND)
        call MPI_ISEND(ubufyo(:,:,:,ivar1:ivar2),nbufy,mpi_precision, &
                       yuneigh,touppys,comm,isend_rq_touppy,mpierr)
      endif
!
!  Set communication across z-planes.
!
      if (nprocz>1) then
!
        if (lyinyang.and.lfirst_proc_z) then
!
!  Interpolate variables ivar1 ... ivar2 for lower z neighbor in other grid.
!  Result in lbufzo.
!
          call interpolate_yy(f,ivar1,ivar2,lbufzo,MIDY,iyinyang_intpol_type)
          comm=MPI_COMM_PENCIL; touppzr=zlneigh; tolowzs=iproc_world  !touppzr=MPI_ANY_TAG
        elseif (lcoarse) then
          do mm=mexts(1),mexts(2)
            lbufzo(:,mm-nghost,:,ivar1:ivar2)=f(:,mm,nexts(mm,1):nexts(mm,1)+(nghost-1)*nphis(mm):nphis(mm),ivar1:ivar2)
          enddo
          if (lfirst_proc_y.and.mexts(2)<m2) then
            lbufzo(:,mexts(2)-nghost+1:,:,ivar1:ivar2)=f(:,mexts(2)+1:m2,n1:n1i,ivar1:ivar2)
          elseif (llast_proc_y.and.mexts(1)>m1) then
            lbufzo(:,:mexts(1)-nghost-1,:,ivar1:ivar2)=f(:,m1:mexts(1)-1,n1:n1i,ivar1:ivar2)
          endif
          comm=MPI_COMM_PENCIL
        else
          lbufzo(:,:,:,ivar1:ivar2)=f(:,m1:m2,n1:n1i,ivar1:ivar2) !lower z-planes
          comm=MPI_COMM_GRID
        endif

        nbufz=bufact*bufsizes_yz(INZL,IRCV)
!if (ldiagnos.and.lfirst.and.iproc_world==0) print*, iproc_world, ' receives', bufsizes_yz(INZL,IRCV), ' from', zlneigh
        call MPI_IRECV(lbufzi(:,:,:,ivar1:ivar2),nbufz,mpi_precision, &
                       zlneigh,touppzr,comm,irecv_rq_fromlowz,mpierr)

        nbufz=bufact*bufsizes_yz(INZL,ISND)
        call MPI_ISEND(lbufzo(:,:,:,ivar1:ivar2),nbufz,mpi_precision, &
                       zlneigh,tolowzs,comm,isend_rq_tolowz,mpierr)
!
        if (lyinyang.and.llast_proc_z) then
!
!  Interpolate variables ivar1 ... ivar2 for upper z neighbor in other grid.
!  Result in ubufyo.
!
          call interpolate_yy(f,ivar1,ivar2,ubufzo,MIDY,iyinyang_intpol_type)
          comm=MPI_COMM_PENCIL; tolowzr=zuneigh; touppzs=iproc_world   !tolowzr=MPI_ANY_TAG
        elseif (lcoarse) then
          do mm=mexts(1),mexts(2)
            ubufzo(:,mm-nghost,:,ivar1:ivar2)=f(:,mm,nexts(mm,2)-(nghost-1)*nphis(mm):nexts(mm,2):nphis(mm),ivar1:ivar2)
          enddo
          if (lfirst_proc_y.and.mexts(2)<m2) then
            ubufzo(:,mexts(2)-nghost+1:,:,ivar1:ivar2)=f(:,mexts(2)+1:m2,n2i:n2,ivar1:ivar2)
          elseif (llast_proc_y.and.mexts(1)>m1) then
            ubufzo(:,:mexts(1)-nghost-1,:,ivar1:ivar2)=f(:,m1:mexts(1)-1,n2i:n2,ivar1:ivar2)
          endif
          comm=MPI_COMM_PENCIL
        else
          ubufzo(:,:,:,ivar1:ivar2)=f(:,m1:m2,n2i:n2,ivar1:ivar2) !upper z-planes
          comm=MPI_COMM_GRID
        endif

        nbufz=bufact*bufsizes_yz(INZU,IRCV)
        call MPI_IRECV(ubufzi(:,:,:,ivar1:ivar2),nbufz,mpi_precision, &
                       zuneigh,tolowzr,comm,irecv_rq_fromuppz,mpierr)

        nbufz=bufact*bufsizes_yz(INZU,ISND)
        call MPI_ISEND(ubufzo(:,:,:,ivar1:ivar2),nbufz,mpi_precision, &
                       zuneigh,touppzs,comm,isend_rq_touppz,mpierr)
!do j=ivar1,ivar2
!if (notanumber(ubufzo(:,:,:,j))) print*, 'ubufzo: iproc,j=', iproc,j
!if (notanumber(lbufzo(:,:,:,:))) print*, 'lbufzo: iproc,j=', iproc
!if (notanumber(ubufyo(:,:,:,j))) print*, 'ubufyo: iproc,j=', iproc,j
!if (notanumber(lbufyo(:,:,:,:))) print*, 'lbufyo: iproc,j=', iproc
!enddo
      endif
!
!  The four corners (in counter-clockwise order).
!  (NOTE: this should work even for nprocx>1)
!
      if (nprocz>1.and.(nprocy>1.or.lyinyang)) then
!
!  Internal and periodic yz-corner buffers.
!
        bufact=mxl*(ivar2-ivar1+1)

!  North/south-pole yz-corner buffers and message tags. NB N:0/S:pi radians.
!  Send N(S) lower(upper) buffers to pi-shifted iprocz +/-1 and receive
!  pi-shifted iprocz -/+1 on the same y-plane.
!
        if (lcommunicate_y) then
          if (lfirst_proc_y) then !N-pole
            if (lcoarse) then
              do mm=min(mexts(2),m1i),m1,-1
                lubufo(:,2*nghost-mm+1,:,ivar1:ivar2)=f(:,mm,nexts(mm,2)-(nghost-1)*nphis(mm):nexts(mm,2):nphis(mm),ivar1:ivar2)
                llbufo(:,2*nghost-mm+1,:,ivar1:ivar2)=f(:,mm,nexts(mm,1):nexts(mm,1)+(nghost-1)*nphis(mm):nphis(mm),ivar1:ivar2)
              enddo
              if (m1i>mexts(2)) then
                lubufo(:,:m1i-mexts(2),:,ivar1:ivar2)= f(:,m1i:mexts(2)+1:-1,n2i:n2,ivar1:ivar2)
                llbufo(:,:m1i-mexts(2),:,ivar1:ivar2)= f(:,m1i:mexts(2)+1:-1,n1:n1i,ivar1:ivar2)
              endif
            else
              lubufo(:,:,:,ivar1:ivar2)=f(:,m1i:m1:-1,n2i:n2,ivar1:ivar2)
              llbufo(:,:,:,ivar1:ivar2)=f(:,m1i:m1:-1,n1:n1i,ivar1:ivar2)
            endif
            lucorns=pnfcrn;lucornr=pnbcrn;llcornr=pnfcrn;llcorns=pnbcrn
            TOlus  =TOnf;  TOulr  =TOnf;  TOuur  =TOnb;  TOlls  =TOnb
          endif
          if (llast_proc_y) then !S-pole
            if (lcoarse) then
              do mm=m2,max(mexts(1),m2i),-1
                ulbufo(:,m2-mm+1,:,ivar1:ivar2)=f(:,mm,nexts(mm,1):nexts(mm,1)+(nghost-1)*nphis(mm):nphis(mm),ivar1:ivar2)
                uubufo(:,m2-mm+1,:,ivar1:ivar2)=f(:,mm,nexts(mm,2)-(nghost-1)*nphis(mm):nexts(mm,2):nphis(mm),ivar1:ivar2)
              enddo
              if (m2i<mexts(1)) then
                ulbufo(:,m2-mexts(1)+2:,:,ivar1:ivar2)=f(:,mexts(1)-1:m2i:-1,n1:n1i,ivar1:ivar2)
                uubufo(:,m2-mexts(1)+2:,:,ivar1:ivar2)=f(:,mexts(1)-1:m2i:-1,n2i:n2,ivar1:ivar2)
              endif
            else
              ulbufo(:,:,:,ivar1:ivar2)=f(:,m2:m2i:-1,n1:n1i,ivar1:ivar2)
              uubufo(:,:,:,ivar1:ivar2)=f(:,m2:m2i:-1,n2i:n2,ivar1:ivar2)
            endif
            uucorns=psfcrn;uucornr=psbcrn;ulcornr=psfcrn;ulcorns=psbcrn
            TOuus  =TOsf;  TOllr  =TOsf;  TOlur  =TOsb;  TOuls  =TOsb
          endif
        endif
!
!  Lower y, lower z.
!
        if (lyinyang.and.(lfirst_proc_z.or.lfirst_proc_y)) then
!
!  "Translation" from (lower y, upper y) x (lower z, upper z) to "LEFT/RIGHT",
!  the latter with respect to looking out of one grid towards the other.
!
          if (lfirst_proc_yz) then
            dir=NIL
          elseif (lfirst_proc_z) then
            dir=RIGHT
          else
            dir=LEFT
          endif
!
!  Interpolate variables ivar1 ... ivar2 for lower-lower corner neighbor in other grid.
!  Result in llbufo.
!
          if (llcorns>=0) call interpolate_yy(f,ivar1,ivar2,llbufo,dir,iyinyang_intpol_type)
          comm=MPI_COMM_PENCIL; TOuur=llcornr; TOlls=iproc_world   !TOuur=MPI_ANY_TAG
        elseif (.not.(lcommunicate_y.and.(lfirst_proc_y.or.llast_proc_y))) then
          llbufo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n1:n1i,ivar1:ivar2)
          comm=MPI_COMM_GRID
        endif

        nbufyz=bufact*product(bufsizes_yz_corn(:,INLL,IRCV))
        if (llcornr>=0) call MPI_IRECV(llbufi(:,:,:,ivar1:ivar2),nbufyz,mpi_precision, &
                                       llcornr,TOuur,comm,irecv_rq_FRll,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INLL,ISND))
        if (llcorns>=0) call MPI_ISEND(llbufo(:,:,:,ivar1:ivar2),nbufyz,mpi_precision, &
                                       llcorns,TOlls,comm,isend_rq_TOll,mpierr)
!
!  Upper y, lower z.
!
        if (lyinyang.and.(lfirst_proc_z.or.llast_proc_y)) then
          if (lfirst_proc_z.and.llast_proc_y) then
            dir=NIL
          elseif (lfirst_proc_z) then
            dir=LEFT
          else
            dir=RIGHT
          endif
! if (iproc_world==44) print*, 'size(intcoeffs)=', size(intcoeffs(RIGHT)%inds,1), &
! size(intcoeffs(RIGHT)%inds,2), size(intcoeffs(RIGHT)%inds,3)
!
!  Interpolate variables ivar1 ... ivar2 for upper-lower corner neighbor in other grid.
!  Result in ulbufo.
!
          if (ulcorns>=0) call interpolate_yy(f,ivar1,ivar2,ulbufo,dir,iyinyang_intpol_type)
          comm=MPI_COMM_PENCIL; TOlur=ulcornr; TOuls=iproc_world   !TOlur=MPI_ANY_TAG
        elseif (.not.(lcommunicate_y.and.(lfirst_proc_y.or.llast_proc_y))) then
          ulbufo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n1:n1i,ivar1:ivar2)
          comm=MPI_COMM_GRID
        endif

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUL,IRCV))
! if (ldiagnos.and.lfirst.and.iproc_world==0) print*, iproc_world, ' receives', bufsizes_yz_corn(:,INUL,IRCV), &
!' from', ulcornr
        if (ulcornr>=0) call MPI_IRECV(ulbufi(:,:,:,ivar1:ivar2),nbufyz,mpi_precision, &
                                       ulcornr,TOlur,comm,irecv_rq_FRul,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUL,ISND))
        if (ulcorns>=0) call MPI_ISEND(ulbufo(:,:,:,ivar1:ivar2),nbufyz,mpi_precision, &
                                       ulcorns,TOuls,comm,isend_rq_TOul,mpierr)
!
!  Upper y, upper z.
!
        if (lyinyang.and.(llast_proc_z.or.llast_proc_y)) then
          if (llast_proc_y.and.llast_proc_z) then
            dir=NIL
          elseif (llast_proc_z) then
            dir=RIGHT
          else
            dir=LEFT
          endif
!
!  Interpolate variables ivar1 ... ivar2 for upper_upper corner neighbor in other grid.
!  Result in uubufo.
!
          if (uucorns>=0) call interpolate_yy(f,ivar1,ivar2,uubufo,dir,iyinyang_intpol_type)
          comm=MPI_COMM_PENCIL; TOllr=uucornr; TOuus=iproc_world   !TOllr=MPI_ANY_TAG
        elseif (.not.(lcommunicate_y.and.(lfirst_proc_y.or.llast_proc_y))) then
          uubufo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n2i:n2,ivar1:ivar2)
          comm=MPI_COMM_GRID
        endif

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUU,IRCV))
        if (uucornr>=0) call MPI_IRECV(uubufi(:,:,:,ivar1:ivar2),nbufyz,mpi_precision, &
                                       uucornr,TOllr,comm,irecv_rq_FRuu,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INUU,ISND))
        if (uucorns>=0) call MPI_ISEND(uubufo(:,:,:,ivar1:ivar2),nbufyz,mpi_precision, &
                                       uucorns,TOuus,comm,isend_rq_TOuu,mpierr)
!
!  Lower y, upper z.
!
        if (lyinyang.and.(llast_proc_z.or.lfirst_proc_y)) then
          if (lfirst_proc_y.and.llast_proc_z) then
            dir=NIL
          elseif (llast_proc_z) then
            dir=LEFT
          else
            dir=RIGHT
          endif
!
!  Interpolate variables ivar1 ... ivar2 for lower_upper corner neighbor in other grid.
!  Result in lubufo.
!
          if (lucorns>=0) call interpolate_yy(f,ivar1,ivar2,lubufo,dir,iyinyang_intpol_type)
          comm=MPI_COMM_PENCIL; TOulr=lucornr; TOlus=iproc_world   !TOulr=MPI_ANY_TAG
        elseif (.not.(lcommunicate_y.and.(lfirst_proc_y.or.llast_proc_y))) then
          lubufo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n2i:n2,ivar1:ivar2)
          comm=MPI_COMM_GRID
        endif

        nbufyz=bufact*product(bufsizes_yz_corn(:,INLU,IRCV))
!if (ldiagnos.and.lfirst.and.iproc_world==0) print*, iproc_world, ' receives', bufsizes_yz_corn(:,INLU,IRCV), &
!' from', lucornr
        if (lucornr>=0) call MPI_IRECV(lubufi(:,:,:,ivar1:ivar2),nbufyz,mpi_precision, &
                                       lucornr,TOulr,comm,irecv_rq_FRlu,mpierr)

        nbufyz=bufact*product(bufsizes_yz_corn(:,INLU,ISND))
        if (lucorns>=0) call MPI_ISEND(lubufo(:,:,:,ivar1:ivar2),nbufyz,mpi_precision, &
                                       lucorns,TOlus,comm,isend_rq_TOlu,mpierr)
!
      endif
!if (itsub>=3.and.it>116) write(78,*) 'after corner comm it,itsub,iproc=', &
!it,itsub,iproc
!
!  communication sample
!  (commented out, because compiler does like this for 0-D runs)
!
!     if (ip<7.and.lfirst_proc_y.and.ipz==3) &
!       print*,'initiate_isendrcv_bdry: MPICOMM send lu: ',iproc,lubufo(nx/2+4,:,1,2),' to ',lucorn
!
    endsubroutine initiate_isendrcv_bdry
!***********************************************************************
    subroutine finalize_isendrcv_bdry(f,ivar1_opt,ivar2_opt)
!
!  Make sure the communications initiated with initiate_isendrcv_bdry are
!  finished and insert the just received boundary values.
!  Receive requests do not need to (and on OSF1 cannot) be explicitly
!  freed, since MPI_Wait takes care of this.
!
!  21-may-02/axel: communication of corners added
!  22-oct-15/fred: communication for spherical polar coords across poles
!                  added yz corner communication
!  20-dec-15/MR: modified for Yin-Yang grid
!
      use General, only: transpose_mn, notanumber, copy_kinked_strip_z, copy_kinked_strip_y, reset_triangle, notanumber

      real, dimension(mx,my,mz,mfarray), intent(inout):: f
      integer, optional,                 intent(in)   :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, j
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
      if (ivar2==0) return
!
!  1. wait until data received
!  2. set ghost zones
!  3. wait until send completed, will be overwritten in next time step
!
!  Communication across y-planes (includes periodic bc)
!
      if (nprocy>1.or.lyinyang) then

        call MPI_WAIT(irecv_rq_fromuppy,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowy,irecv_stat_fl,mpierr)

        do j=ivar1,ivar2

          if (.not. lfirst_proc_y .or. bcy12(j,1)=='p' .or. lyinyang) then
! communication should happen only under these conditions
            if (lfirst_proc_y.and.lyinyang) then
              if (lcorner_yz) then
                if (.not.lcutoff_corners) f(:,:m1-1,:,j)=lbufyi(:,:,:mz,j)    ! fill left vertical ghost strip
if (notanumber(lbufyi(:,:,:mz,j))) print*, 'lbufyi(1:mz): iproc,j=', iproc, iproc_world, j
if (notanumber(lbufyi(:,:,mz+1:,j))) print*, 'lbufyi(mz+1:): iproc,j=', iproc, iproc_world, j
                if (lfirst_proc_z) then
                  if (lcutoff_corners) then
                    call copy_kinked_strip_y(nzcut,my-nycut+2,lbufyi,f,j,-1,.true.)
                    call reset_triangle(1,my-nycut-nghost,1,mz-nzcut-nghost,f(:,:,:,j))
                  endif
                  call transpose_mn( lbufyi(:,:,zcornstart:,j), f(:,my-len_cornstrip_y+1:,:nghost,j))   ! fill lower horizontal ghost strip
                elseif (llast_proc_z) then
                  if (lcutoff_corners) then
                    call copy_kinked_strip_y(nzcut,1,lbufyi,f,j,1,.false.)
                    call reset_triangle(1,my-nycut-nghost,mz,nzcut+1+nghost,f(:,:,:,j))
                  endif
if (notanumber(f(:,:,:,j))) print*, 'f nach lbufyi: iproc,j=', iproc, iproc_world, j
                  call transpose_mn( lbufyi(:,:,zcornstart:,j), f(:,my-len_cornstrip_y+1:,n2+1:,j))   ! fill upper horizontal ghost strip
                endif
              else
                f(:,:m1-1,:,j)=lbufyi(:,:,:mz,j)
              endif
            else
              f(:,1:m1-1,n1:n2,j)=lbufyi(:,:,:,j)                           ! set lower y-ghostzone
            endif
          endif

          if (.not. llast_proc_y .or. bcy12(j,2)=='p' .or.lyinyang) then
            if (llast_proc_y.and.lyinyang) then
              if (lcorner_yz) then
                if (.not.lcutoff_corners) f(:,m2+1:,:,j)=ubufyi(:,:,:mz,j)   ! fill right vertical ghost strip
if (notanumber(ubufyi(:,:,:mz,j))) print*, 'ubufyi(1:mz): iproc,j=', iproc, iproc_world, j
if (notanumber(ubufyi(:,:,mz+1:,j))) print*, 'ubufyi(mz+1:): iproc,j=', iproc, iproc_world, j
                if (lfirst_proc_z) then
                  if (lcutoff_corners) then
                    call copy_kinked_strip_y(nzcut,nycut-nghost,ubufyi,f,j,1,.true.)
                    call reset_triangle(my,nycut+1+nghost,1,mz-nzcut-nghost,f(:,:,:,j))
                  endif
                  call transpose_mn( ubufyi(:,:,zcornstart:,j), f(:,1:len_cornstrip_y,:nghost,j)) ! fill lower horizontal ghost strip
                elseif (llast_proc_z) then
                  if (lcutoff_corners) then
                    call copy_kinked_strip_y(nzcut,m2+1,ubufyi,f,j,-1,.false.)
                    call reset_triangle(my,nycut+1+nghost,mz,nzcut+1+nghost,f(:,:,:,j))
                  endif
                  call transpose_mn( ubufyi(:,:,zcornstart:,j), f(:,1:len_cornstrip_y,n2+1:,j)) ! fill upper horizontal ghost strip
                endif
              else
                f(:,m2+1:,:,j)=ubufyi(:,:,:mz,j)
              endif
!if (notanumber(f(:,:,:,j))) print*, 'ubufyi: iproc,j=', iproc, iproc_world, j
            else
              f(:,m2+1:,n1:n2,j)=ubufyi(:,:,:,j)                           ! set upper y-ghostzone
            endif
          endif
!
!  Spherical y-planes across the poles
!
          if (lcommunicate_y) then
            if (lfirst_proc_y) then
              if (bcy12(j,1)=='pp') f(:, 1:m1-1,n1:n2,j)= lbufyi(:,:,:,j) !N-pole
              if (bcy12(j,1)=='ap') f(:, 1:m1-1,n1:n2,j)=-lbufyi(:,:,:,j) !N-pole
            endif
            if (llast_proc_y) then
              if (bcy12(j,2)=='pp') f(:,m2+1:  ,n1:n2,j)= ubufyi(:,:,:,j) !S-pole
              if (bcy12(j,2)=='ap') f(:,m2+1:  ,n1:n2,j)=-ubufyi(:,:,:,j) !S-pole
            endif
          endif

        enddo

        call MPI_WAIT(isend_rq_tolowy,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppy,isend_stat_tu,mpierr)

      endif
!
!  Communication in z (includes periodic bc)
!
      if (nprocz>1) then

        call MPI_WAIT(irecv_rq_fromuppz,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowz,irecv_stat_fl,mpierr)

        do j=ivar1,ivar2

          if (.not. lfirst_proc_z .or. bcz12(j,1)=='p'.or.lyinyang) then
            if (lfirst_proc_z.and.lyinyang) then
              if (lcorner_yz) then
! fill lower horizontal ghost strip
                if (.not.lcutoff_corners) f(:,:,:n1-1,j)=f(:,:,:n1-1,j)+lbufzi(:,:my,:,j)
                if (lfirst_proc_y) then
                  if (lcutoff_corners) call copy_kinked_strip_z(nycut,mz-nzcut+2,lbufzi,f,j,-1,.true.,ladd_=.true.)
! fill left vertical ghost strip
                  call transpose_mn( lbufzi(:,ycornstart:,:,j), f(:,:nghost,mz-len_cornstrip_z+1:,j),ladd=.true.)
                elseif (llast_proc_y) then
                  if (lcutoff_corners) call copy_kinked_strip_z(nycut,1,lbufzi,f,j,1,.false.,ladd_=.true.)
! fill right vertical ghost strip
                  call transpose_mn( lbufzi(:,ycornstart:,:,j), f(:,m2+1:,mz-len_cornstrip_z+1:,j),ladd=.true.)
                endif
!if (notanumber(f(:,:,:,j))) print*, 'lbufzi: iproc,j=', iproc, iproc_world, j
              else
                f(:,:,:n1-1,j)=lbufzi(:,:my,:,j)
              endif
            else
              f(:,m1:m2,1:n1-1,j)=lbufzi(:,:,:,j)  ! set lower z-ghostzone
            endif
          endif

          if (.not. llast_proc_z .or. bcz12(j,2)=='p'.or.lyinyang) then
            if (llast_proc_z.and.lyinyang) then
              if (lcorner_yz) then
! fill upper horizontal ghost strip
                if (.not.lcutoff_corners) f(:,:,n2+1:,j)=f(:,:,n2+1:,j)+ubufzi(:,:my,:,j)
                if (lfirst_proc_y) then
                  if (lcutoff_corners) call copy_kinked_strip_z(nycut,nzcut-nghost,ubufzi,f,j,1,.true.,ladd_=.true.)
! fill left vertical ghost strip
                  call transpose_mn( ubufzi(:,ycornstart:,:,j), f(:,:nghost,:len_cornstrip_z,j),ladd=.true.)
                elseif (llast_proc_y) then
                  if (lcutoff_corners) call copy_kinked_strip_z(nycut,n2+1,ubufzi,f,j,-1,.false.,ladd_=.true.)
! fill right vertical ghost strip
                  call transpose_mn( ubufzi(:,ycornstart:,:,j), f(:,m2+1:,:len_cornstrip_z,j),ladd=.true.)
                endif
              else
                f(:,:,n2+1:,j)=ubufzi(:,:my,:,j)  ! fill upper horizontal ghost strip
              endif
            else
              f(:,m1:m2,n2+1:,j)=ubufzi(:,:,:,j)  ! set upper z-ghostzone
            endif
          endif

!if (notanumber(ubufzi(:,:,:,j))) print*, 'ubufzi: iproc,j=', iproc, j
!if (notanumber(lbufzi(:,:,:,j))) print*, 'lbufzi: iproc,j=', iproc, j
!if (notanumber(ubufyi(:,:,:,j))) print*, 'ubufyi: iproc,j=', iproc, j
!if (notanumber(lbufyi(:,:,:,j))) print*, 'lbufyi: iproc,j=', iproc, j
        enddo
        call MPI_WAIT(isend_rq_tolowz,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppz,isend_stat_tu,mpierr)

      endif
!
!  The four yz-corners (in counter-clockwise order)
!
       if (nprocz>1.and.(nprocy>1.or.lyinyang)) then

        if (uucornr>=0) call MPI_WAIT(irecv_rq_FRuu,irecv_stat_Fuu,mpierr)
        if (lucornr>=0) call MPI_WAIT(irecv_rq_FRlu,irecv_stat_Flu,mpierr)
        if (llcornr>=0) call MPI_WAIT(irecv_rq_FRll,irecv_stat_Fll,mpierr)
        if (ulcornr>=0) call MPI_WAIT(irecv_rq_FRul,irecv_stat_Ful,mpierr)

        do j=ivar1,ivar2
!
!  Set ll corner
!
          if (.not. lfirst_proc_z .or. bcz12(j,1)=='p' .or. lyinyang) then

            if (.not. lfirst_proc_y .or. bcy12(j,1)=='p' .or. lyinyang) then

              if ((.not.lfirst_proc_z.or.bcz12(j,1)=='p').and. &
                  (.not.lfirst_proc_y.or.bcy12(j,1)=='p')) then    ! inner or periodic proc boundaries
                f(:,1:m1-1,1:n1-1,j)=llbufi(:,:,:,j)               ! fill lower left corner
!if (notanumber(llbufi(:,:,:,j))) print*, 'lbufyi: iproc,j=', iproc, j
              elseif (llcornr>=0.and.lyinyang) then
                if (lfirst_proc_z) then
                  if (llast_proc_y.and.lcutoff_corners) then
                    call copy_kinked_strip_z(nycut,1,llbufi,f,j,1,.false.,ladd_=.true.)
                  else
                    f(:,:,1:n1-1,j)=f(:,:,1:n1-1,j)+llbufi(:,:,:,j)  ! complete lower horizontal strip
                  endif
                elseif (lfirst_proc_y) then
                  if (llast_proc_z.and.lcutoff_corners) then
                    call copy_kinked_strip_y(nzcut,1,llbufi,f,j,1,.false.,ladd_=.true.)
                  else
                    f(:,1:m1-1,:,j)=f(:,1:m1-1,:,j)+llbufi(:,:,:,j)  ! complete left vertical strip
                  endif
                endif
!if (notanumber(f(:,:,:,j))) print*, 'llcorn: iproc,j=', iproc, iproc_world, j
              endif

            endif
!
!  Set ul corner
!
            if (.not. llast_proc_y .or. bcy12(j,2)=='p' .or. lyinyang) then
              if ((.not.lfirst_proc_z.or.bcz12(j,1)=='p').and. &
                  (.not.llast_proc_y .or.bcy12(j,2)=='p')) then    ! inner or periodic proc boundaries
                f(:,m2+1:,1:n1-1,j)=ulbufi(:,:,:,j)                ! fill lower right corner
              elseif (ulcornr>=0.and.lyinyang) then
                if (lfirst_proc_z) then
                  if (lfirst_proc_y.and.lcutoff_corners) then
                    call copy_kinked_strip_z(nycut,mz-nzcut+2,ulbufi,f,j,-1,.true.,ladd_=.true.)
                  else
                    f(:,:,1:n1-1,j)=f(:,:,1:n1-1,j)+ulbufi(:,:,:,j)  ! complete lower horizontal strip
                  endif
                elseif (llast_proc_y) then
                  if (llast_proc_z.and.lcutoff_corners) then
                    call copy_kinked_strip_y(nzcut,m2+1,ulbufi,f,j,-1,.false.,ladd_=.true.)
                  else
                    f(:,m2+1:,:,j)=f(:,m2+1:,:,j)+ulbufi(:,:,:,j)    ! complete right vertical strip
                  endif
                endif
!if (notanumber(f(:,:,:,j))) print*, 'ulcorn: iproc,j=', iproc, iproc_world, j
              endif
            endif

          endif
!
!  Set uu corner
!
          if (.not. llast_proc_z .or. bcz12(j,2)=='p' .or. lyinyang) then

            if (.not. llast_proc_y .or. bcy12(j,2)=='p' .or. lyinyang) then
              if ((.not.llast_proc_z.or.bcz12(j,2)=='p').and. &   ! inner or periodic proc boundaries
                  (.not.llast_proc_y.or.bcy12(j,2)=='p')) then
                f(:,m2+1:,n2+1:,j)=uubufi(:,:,:,j)                ! fill upper right corner
              elseif (uucornr>=0.and.lyinyang) then
                if (llast_proc_z) then
                  if (lfirst_proc_y.and.lcutoff_corners) then
                    call copy_kinked_strip_z(nycut,nzcut-nghost,uubufi,f,j,1,.true.,ladd_=.true.)
                  else
                    f(:,:,n2+1:,j)=f(:,:,n2+1:,j)+uubufi(:,:,:,j)   ! complete upper horizontal strip
                  endif
                elseif (llast_proc_y) then
                  if (lfirst_proc_z.and.lcutoff_corners) then
                    call copy_kinked_strip_y(nzcut,nycut-nghost,uubufi,f,j,1,.true.,ladd_=.true.)
                  else
                    f(:,m2+1:,:,j)=f(:,m2+1:,:,j)+uubufi(:,:,:,j)   ! complete right vertical strip
                  endif
                endif
!if (notanumber(f(:,:,:,j))) print*, 'uucorn: iproc,j=', iproc, iproc_world, j
              endif
            endif
!
!  Set lu corner
!
            if (.not. lfirst_proc_y .or. bcy12(j,1)=='p' .or. lyinyang) then
              if ((.not.llast_proc_z .or.bcz12(j,2)=='p').and. &   ! inner or periodic proc boundaries
                  (.not.lfirst_proc_y.or.bcy12(j,1)=='p')) then
                f(:,1:m1-1,n2+1:,j)=lubufi(:,:,:,j)                ! fill upper left corner
              elseif (lucornr>=0.and.lyinyang) then
                if (llast_proc_z) then
                  if (llast_proc_y.and.lcutoff_corners) then
                    call copy_kinked_strip_z(nycut,n2+1,lubufi,f,j,-1,.false.,ladd_=.true.)
                  else
                    f(:,:,n2+1:,j)=f(:,:,n2+1:,j)+lubufi(:,:,:,j)    ! complete upper horizontal strip
                  endif
                elseif (lfirst_proc_y) then
                  if (lfirst_proc_z.and.lcutoff_corners) then
                    call copy_kinked_strip_y(nzcut,my-nycut+2,lubufi,f,j,-1,.true.,ladd_=.true.)
                  else
                    f(:,1:m1-1,:,j)=f(:,1:m1-1,:,j)+lubufi(:,:,:,j)  ! complete left vertical strip
                  endif
                endif
!if (notanumber(f(:,:,:,j))) print*, 'lucorn: iproc,j=', iproc, iproc_world, j
              endif
            endif

          endif
!
!  Spherical yz-corners across the poles. If lcommunicate_y then periodic in z
!
          if (lcommunicate_y) then

            if (lfirst_proc_y) then
              if (bcy12(j,1)=='pp') then !N-pole periodic
                f(:,:m1-1,:n1-1,j)= lubufi(:,:,:,j)  !(set ll corner)
                f(:,:m1-1,n2+1:,j)= llbufi(:,:,:,j)  !(set lu corner)
              endif
              if (bcy12(j,1)=='ap') then !N-pole antiperiodic
                f(:,:m1-1,:n1-1,j)=-lubufi(:,:,:,j)  !(set ll corner)
                f(:,:m1-1,n2+1:,j)=-llbufi(:,:,:,j)  !(set lu corner)
              endif
            endif

            if (llast_proc_y) then
              if (bcy12(j,2)=='pp') then !S-pole periodic
                f(:,m2+1:,:n1-1,j)= uubufi(:,:,:,j)  !(set ul corner)
                f(:,m2+1:,n2+1:,j)= ulbufi(:,:,:,j)  !(set uu corner)
              endif
              if (bcy12(j,2)=='ap') then !S-pole antiperiodic
                f(:,m2+1:,:n1-1,j)=-uubufi(:,:,:,j)  !(set ul corner)
                f(:,m2+1:,n2+1:,j)=-ulbufi(:,:,:,j)  !(set uu corner)
              endif
            endif
          endif
        enddo

        if (llcorns>=0) call MPI_WAIT(isend_rq_TOll,isend_stat_Tll,mpierr)
        if (ulcorns>=0) call MPI_WAIT(isend_rq_TOul,isend_stat_Tul,mpierr)
        if (uucorns>=0) call MPI_WAIT(isend_rq_TOuu,isend_stat_Tuu,mpierr)
        if (lucorns>=0) call MPI_WAIT(isend_rq_TOlu,isend_stat_Tlu,mpierr)

      endif
!
!  communication sample
!  (commented out, because compiler does like this for 0-D runs)
!
!     if (ip<7.and.ipy==3.and.lfirst_proc_z) &
!       print*,'finalize_isendrcv_bdry: MPICOMM recv ul: ', &
!                       iproc,ulbufi(nx/2+4,:,1,2),' from ',ulcorn
!
!  make sure the other processors don't carry on sending new data
!  which could be mistaken for an earlier time
!
     call mpibarrier
!
    endsubroutine finalize_isendrcv_bdry
!***********************************************************************
    subroutine isendrcv_bdry_x(f,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values for x-direction. Sends and receives
!  before continuing to y and z boundaries, as this allows the edges
!  of the grid to be set properly.
!
!   2-may-09/anders: coded
!
      real, dimension(:,:,:,:), intent(inout) :: f
      integer, intent(in), optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbufx, j
!
      ivar1=1; ivar2=min(mcom,size(f,4))
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
      if (nprocx>1) then

        lbufxo(:,:,:,ivar1:ivar2)=f(l1:l1i,m1:m2,n1:n2,ivar1:ivar2) !!(lower x-zone)
        ubufxo(:,:,:,ivar1:ivar2)=f(l2i:l2,m1:m2,n1:n2,ivar1:ivar2) !!(upper x-zone)
        nbufx=ny*nz*nghost*(ivar2-ivar1+1)

        call MPI_IRECV(ubufxi(:,:,:,ivar1:ivar2),nbufx,mpi_precision, &
            xuneigh,tolowx,MPI_COMM_GRID,irecv_rq_fromuppx,mpierr)
        call MPI_IRECV(lbufxi(:,:,:,ivar1:ivar2),nbufx,mpi_precision, &
            xlneigh,touppx,MPI_COMM_GRID,irecv_rq_fromlowx,mpierr)
        call MPI_ISEND(lbufxo(:,:,:,ivar1:ivar2),nbufx,mpi_precision, &
            xlneigh,tolowx,MPI_COMM_GRID,isend_rq_tolowx,mpierr)
        call MPI_ISEND(ubufxo(:,:,:,ivar1:ivar2),nbufx,mpi_precision, &
            xuneigh,touppx,MPI_COMM_GRID,isend_rq_touppx,mpierr)
        call MPI_WAIT(irecv_rq_fromuppx,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowx,irecv_stat_fl,mpierr)
!
!  Inner communication or (shear-)periodic boundary conditions in x
!  MR: Communication should only happen under these conditions.
!
        do j=ivar1,ivar2
          if (.not. lfirst_proc_x .or. bcx12(j,1)=='p' .or. &
              (bcx12(j,1)=='she'.and.nygrid==1)) then
            f( 1:l1-1,m1:m2,n1:n2,j)=lbufxi(:,:,:,j)  !!(set lower buffer)
          endif
          if (.not. llast_proc_x .or. bcx12(j,2)=='p' .or. &
              (bcx12(j,2)=='she'.and.nygrid==1)) then
            f(l2+1:,m1:m2,n1:n2,j)=ubufxi(:,:,:,j)  !!(set upper buffer)
          endif
        enddo

        call MPI_WAIT(isend_rq_tolowx,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppx,isend_stat_tu,mpierr)

      endif
!
    endsubroutine isendrcv_bdry_x
!***********************************************************************
    subroutine periodic_bdry_x(buffer,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values for x-direction. Sends and receives
!  before continuing to y and z boundaries, as this allows the edges
!  of the grid to be set properly.
!
!   19-jan-23/fred: coded
!
      real, dimension(:,:,:,:), intent(inout) :: buffer
      integer, intent(in), optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbuf, j, ny_loc, nz_loc, nvar_loc
      integer :: isend_rq, irecv_rq
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat, isend_stat
      real, dimension(:,:,:,:), allocatable :: bufo,bufi
      integer :: comm

      comm = MPI_COMM_GRID

!
      if (.not.(lfirst_proc_x.or.llast_proc_x)) return
!
      if (nprocx>1) then
!
        ny_loc=size(buffer,2)-2*nghost
        nz_loc=size(buffer,3)-2*nghost
        ivar1=1; ivar2=min(mcom,size(buffer,4))
        if (present(ivar1_opt)) ivar1=ivar1_opt
        if (present(ivar2_opt)) ivar2=ivar2_opt
        nvar_loc=ivar2-ivar1+1
        allocate(bufo(nghost,ny_loc,nz_loc,nvar_loc),&
                 bufi(nghost,ny_loc,nz_loc,nvar_loc))
!

        nbuf=ny_loc*nz_loc*nghost*nvar_loc
        if (lfirst_proc_x) then
          bufo=buffer(l1:l1i,m1:m1+ny_loc-1,n1:n1+nz_loc-1,ivar1:ivar2) !!(lower x-zone)
          call MPI_ISEND(bufo,nbuf,mpi_precision, &
              xlneigh,tolowx,comm,isend_rq,mpierr)
          call MPI_IRECV(bufi,nbuf,mpi_precision, &
              xlneigh,touppx,comm,irecv_rq,mpierr)
          call MPI_WAIT(irecv_rq,irecv_stat,mpierr)
        else
          bufo=buffer(l2i:l2,m1:m1+ny_loc-1,n1:n1+nz_loc-1,ivar1:ivar2) !!(upper x-zone)
          call MPI_IRECV(bufi,nbuf,mpi_precision, &
              xuneigh,tolowx,comm,irecv_rq,mpierr)
          call MPI_ISEND(bufo,nbuf,mpi_precision, &
              xuneigh,touppx,comm,isend_rq,mpierr)
          call MPI_WAIT(irecv_rq,irecv_stat,mpierr)
        endif
!
!  Inner communication or (shear-)periodic boundary conditions in x
!  MR: Communication should only happen under these conditions.
!
        do j=ivar1,ivar2
          if (.not. lfirst_proc_x .or. bcx12(j,1)=='p' .or. &
              (bcx12(j,1)=='she'.and.nygrid==1)) then
            buffer(1:l1-1,m1:m1+ny_loc-1,n1:n1+nz_loc-1,j)=bufi(:,:,:,j)  !!(set lower buffer)
          endif
          if (.not. llast_proc_x .or. bcx12(j,2)=='p' .or. &
              (bcx12(j,2)=='she'.and.nygrid==1)) then
            buffer(l2+1: ,m1:m1+ny_loc-1,n1:n1+nz_loc-1,j)=bufi(:,:,:,j)  !!(set upper buffer)
          endif
        enddo

        call MPI_WAIT(isend_rq,isend_stat,mpierr)

      endif
!
    endsubroutine periodic_bdry_x
!***********************************************************************
    subroutine periodic_bdry_y(buffer,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values for y-direction.
!
!   19-jan-23/fred: coded
!
      real, dimension(:,:,:,:), intent(inout) :: buffer
      integer, intent(in), optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbuf, j, nx_loc, nz_loc, nvar_loc
      integer :: isend_rq, irecv_rq
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat, isend_stat
      real, dimension(:,:,:,:), allocatable :: bufo,bufi
      integer :: comm
!
      if (.not.(lfirst_proc_y.or.llast_proc_y)) return

      comm = MPI_COMM_GRID
!
      if (nprocy>1) then
!
        nx_loc=size(buffer,1)-2*nghost
        nz_loc=size(buffer,3)-2*nghost
        ivar1=1; ivar2=min(mcom,size(buffer,4))
        if (present(ivar1_opt)) ivar1=ivar1_opt
        if (present(ivar2_opt)) ivar2=ivar2_opt
        nvar_loc=ivar2-ivar1+1
        allocate(bufo(nx_loc,nghost,nz_loc,nvar_loc),&
                 bufi(nx_loc,nghost,nz_loc,nvar_loc))
!
        nbuf=nx_loc*nz_loc*nghost*nvar_loc
        if (lfirst_proc_y) then
          bufo=buffer(l1:l1+nx_loc-1,m1:m1i,n1:n1+nz_loc-1,ivar1:ivar2) !!(lower y-zone)
          call MPI_ISEND(bufo,nbuf,mpi_precision, &
              ylneigh,tolowy,comm,isend_rq,mpierr)
          call MPI_IRECV(bufi,nbuf,mpi_precision, &
              ylneigh,touppy,comm,irecv_rq,mpierr)
          call MPI_WAIT(irecv_rq,irecv_stat,mpierr)
          buffer(l1:l1+nx_loc-1,1:m1-1,n1:n1+nz_loc-1,:)=bufi !!(set lower buffer)
        else
          bufo=buffer(l1:l1+nx_loc-1,m2i:m2,n1:n1+nz_loc-1,ivar1:ivar2) !!(upper y-zone)
          call MPI_IRECV(bufi,nbuf,mpi_precision, &
              yuneigh,tolowy,comm,irecv_rq,mpierr)
          call MPI_ISEND(bufo,nbuf,mpi_precision, &
              yuneigh,touppy,comm,isend_rq,mpierr)
          call MPI_WAIT(irecv_rq,irecv_stat,mpierr)
          buffer(l1:l1+nx_loc-1,m2+1: ,n1:n1+nz_loc-1,:)=bufi  !!(set upper buffer)
        endif
!
        call MPI_WAIT(isend_rq,isend_stat,mpierr)
!
      endif
!
    endsubroutine periodic_bdry_y
!***********************************************************************
    subroutine periodic_bdry_z(buffer,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values for z-direction.
!
!   19-jan-23/fred: coded
!
      real, dimension(:,:,:,:), intent(inout) :: buffer
      integer, intent(in), optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbuf, j, ny_loc, nx_loc, nvar_loc
      integer :: isend_rq, irecv_rq
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat, isend_stat
      real, dimension(:,:,:,:), allocatable :: bufo,bufi
      integer :: comm
!
      if (.not.(lfirst_proc_z.or.llast_proc_z)) return

      comm = MPI_COMM_GRID
!
      if (nprocz>1) then
!
        nx_loc=size(buffer,1)-2*nghost
        ny_loc=size(buffer,2)-2*nghost
        ivar1=1; ivar2=min(mcom,size(buffer,4))
        if (present(ivar1_opt)) ivar1=ivar1_opt
        if (present(ivar2_opt)) ivar2=ivar2_opt
        nvar_loc=ivar2-ivar1+1
        allocate(bufo(nx_loc,ny_loc,nghost,nvar_loc),&
                 bufi(nx_loc,ny_loc,nghost,nvar_loc))
!
        nbuf=nx_loc*ny_loc*nghost*nvar_loc
        if (lfirst_proc_z) then
          bufo=buffer(l1:l1+nx_loc-1,m1:m1+ny_loc-1,n1:n1i,ivar1:ivar2) !!(lower z-zone)
          call MPI_ISEND(bufo,nbuf,mpi_precision, &
              zlneigh,tolowz,comm,isend_rq,mpierr)
          call MPI_IRECV(bufi,nbuf,mpi_precision, &
              zlneigh,touppz,comm,irecv_rq,mpierr)
          call MPI_WAIT(irecv_rq,irecv_stat,mpierr)
          buffer(l1:l1+nx_loc-1,m1:m1+ny_loc-1,1:n1-1,:)=bufi !!(set lower buffer)
        else
          bufo=buffer(l1:l1+nx_loc-1,m1:m1+ny_loc-1,n2i:n2,ivar1:ivar2) !!(upper z-zone)
          call MPI_IRECV(bufi,nbuf,mpi_precision, &
              zuneigh,tolowz,comm,irecv_rq,mpierr)
          call MPI_ISEND(bufo,nbuf,mpi_precision, &
              zuneigh,touppz,comm,isend_rq,mpierr)
          call MPI_WAIT(irecv_rq,irecv_stat,mpierr)
          buffer(l1:l1+nx_loc-1,m1:m1+ny_loc-1,n2+1:,:)=bufi  !!(set upper buffer)
        endif
!
        call MPI_WAIT(isend_rq,isend_stat,mpierr)
!
      endif
!
    endsubroutine periodic_bdry_z
!***********************************************************************
   subroutine initiate_shearing(f,ivar1_opt,ivar2_opt)
!
!  Subroutine for shearing sheet boundary conditions
!
!  27-nov-14/mcnallcp: Now uses 4 shearing neighbours so
!                      y-ghosts are not needed
!  20-june-02/nils: adapted from pencil_mpi
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      real(KIND=rkind8) :: deltay_dy, frac, c1, c2, c3, c4, c5, c6
      integer :: ivar1, ivar2, ystep, nbufx_gh
      integer :: tolastya=11, tolastyb=12, tonextya=13, tonextyb=14
      integer :: tolastlastya=15, tolastlastyb=16, tonextnextya=17, tonextnextyb=18
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  Sixth order interpolation along the y-direction
!
      deltay_dy=deltay/dy
      displs=int(deltay_dy)
      if (nprocx==1 .and. nprocy==1) then
        if (nygrid==1) then ! Periodic boundary conditions.
          f(   1:l1-1,m1:m2,:,ivar1:ivar2) = f(l2i:l2,m1:m2,:,ivar1:ivar2)
          f(l2+1:mx  ,m1:m2,:,ivar1:ivar2) = f(l1:l1i,m1:m2,:,ivar1:ivar2)
        else
          frac=deltay_dy-displs
          c1 = -          (frac+1.)*frac*(frac-1.)*(frac-2.)*(frac-3.)/120.
          c2 = +(frac+2.)          *frac*(frac-1.)*(frac-2.)*(frac-3.)/24.
          c3 = -(frac+2.)*(frac+1.)     *(frac-1.)*(frac-2.)*(frac-3.)/12.
          c4 = +(frac+2.)*(frac+1.)*frac          *(frac-2.)*(frac-3.)/12.
          c5 = -(frac+2.)*(frac+1.)*frac*(frac-1.)          *(frac-3.)/24.
          c6 = +(frac+2.)*(frac+1.)*frac*(frac-1.)*(frac-2.)          /120.
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
        endif
      else
        if (nygrid==1) return ! Periodic boundary conditions already set.
!
!  With more than one CPU in the y-direction it will become necessary to
!  interpolate over data from two different CPUs. Likewise two different
!  CPUs will require data from this CPU.
!
        if (lfirst_proc_x .or. llast_proc_x) then
          ipx_partner=(nprocx-ipx-1)
          ystep = displs/ny
! Terms in these expressions
!     ipz*nprocy*nprocx = the offset to this x-y plane of processors
!     modulo(ipy-ystep,nprocy)*nprocx = the offset to the x-row of the shearing neighbor
!     ipx_partner = the offset to the shearing neighbour within the x-row,
!                   either the low-x end or high-x end of the grid.
!                   Note iproc, ipx, etc are zero-based indexing.
          nextnextya = find_proc(ipx_partner,ipy-ystep+1,ipz)
          nextya     = find_proc(ipx_partner,ipy-ystep  ,ipz)
          lastya     = find_proc(ipx_partner,ipy-ystep-1,ipz)
          lastlastya = find_proc(ipx_partner,ipy-ystep-2,ipz)
!
          lastlastyb = find_proc(ipx_partner,ipy+ystep-1,ipz)
          lastyb     = find_proc(ipx_partner,ipy+ystep  ,ipz)
          nextyb     = find_proc(ipx_partner,ipy+ystep+1,ipz)
          nextnextyb = find_proc(ipx_partner,ipy+ystep+2,ipz)
!
          !!nextnextya = ipz*nprocy*nprocx +modulo(ipy-ystep+1,nprocy)*nprocx + ipx_partner
          !!nextya     = ipz*nprocy*nprocx +modulo(ipy-ystep  ,nprocy)*nprocx + ipx_partner
          !!lastya     = ipz*nprocy*nprocx +modulo(ipy-ystep-1,nprocy)*nprocx + ipx_partner
          !!lastlastya = ipz*nprocy*nprocx +modulo(ipy-ystep-2,nprocy)*nprocx + ipx_partner
!
          !!lastlastyb = ipz*nprocy*nprocx +modulo(ipy+ystep-1,nprocy)*nprocx + ipx_partner
          !!lastyb     = ipz*nprocy*nprocx +modulo(ipy+ystep  ,nprocy)*nprocx + ipx_partner
          !!nextyb     = ipz*nprocy*nprocx +modulo(ipy+ystep+1,nprocy)*nprocx + ipx_partner
          !!nextnextyb = ipz*nprocy*nprocx +modulo(ipy+ystep+2,nprocy)*nprocx + ipx_partner
!
!         The data that gets passed, each set of values goes to 4 places.
!         Only pass active grid points in y, the guard cells are not assumed to be filled yet
          fao(:,:,:,ivar1:ivar2) = f(l1:l1i,m1:m2,:,ivar1:ivar2)
          fbo(:,:,:,ivar1:ivar2) = f(l2i:l2,m1:m2,:,ivar1:ivar2)
!
!         Buffer length for MPI calls, this is the size of the used section of fao and fbo
          nbufx_gh=(my-2*nghost)*mz*nghost*(ivar2-ivar1+1)
!         Here we exchange the fao and fbo data across the shearing boundary.
!         These if statements determinie if we need to copy, or post a MPI send/recieve.
!         Direct copying is done when we discover we are the reciever.
!         the route from send to recieve buffer names is based on values of iproc:
!          nextnextyb -> fao => fahihi
!          nextyb     -> fao => fahi
!          lastyb     -> fao => falo
!          lastlastyb -> fao => falolo
!          nextnextya -> fbo = fahihi
!          nextya     -> fbo = fahi
!          lastya     -> fbo = falo
!          lastlastya -> fbo = falolo
!
!         map of sends
!           lastlastya -> nextnextyb
!           lastya     -> nextyb
!           nextya     -> lastyb
!           nextnextya -> lastlastyb
!
!           lastlastyb -> nextnextya
!           lastyb     -> nextya
!           nextyb     -> lastya
!           nextnextyb -> lastlastya
!
!         Calls to fill the b-side recieve buffers
          if (lastlastya/=iproc) then
            call MPI_ISEND(fao(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,lastlastya, &
                tonextnextyb,MPI_COMM_GRID,isend_rq_tolastlastya,mpierr)
          endif
          if (nextnextyb==iproc) then
            fbhihi(:,:,:,ivar1:ivar2)=fao(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(fbhihi(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,nextnextyb, &
                tonextnextyb,MPI_COMM_GRID,irecv_rq_fromnextnextyb,mpierr)
          endif
!
          if (lastya/=iproc) then
            call MPI_ISEND(fao(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,lastya, &
                tonextyb,MPI_COMM_GRID,isend_rq_tolastya,mpierr)
          endif
          if (nextyb==iproc) then
            fbhi(:,:,:,ivar1:ivar2)=fao(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(fbhi(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,nextyb, &
                tonextyb,MPI_COMM_GRID,irecv_rq_fromnextyb,mpierr)
          endif
!
         if (nextya/=iproc) then
            call MPI_ISEND(fao(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,nextya, &
                tolastyb,MPI_COMM_GRID,isend_rq_tonextya,mpierr)
          endif
          if (lastyb==iproc) then
            fblo(:,:,:,ivar1:ivar2)=fao(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(fblo(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,lastyb, &
                tolastyb,MPI_COMM_GRID,irecv_rq_fromlastyb,mpierr)
          endif
!
          if (nextnextya/=iproc) then
            call MPI_ISEND(fao(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,nextnextya, &
                tolastlastyb,MPI_COMM_GRID,isend_rq_tonextnextya,mpierr)
          endif
          if (lastlastyb==iproc) then
            fblolo(:,:,:,ivar1:ivar2)=fao(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(fblolo(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,lastlastyb, &
                tolastlastyb,MPI_COMM_GRID,irecv_rq_fromlastlastyb,mpierr)
          endif
!         Now fill a-side recieve buffers
          if (lastlastyb/=iproc) then
            call MPI_ISEND(fbo(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,lastlastyb, &
                tonextnextya,MPI_COMM_GRID,isend_rq_tolastlastyb,mpierr)
          endif
          if (nextnextya==iproc) then
            fahihi(:,:,:,ivar1:ivar2)=fbo(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(fahihi(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,nextnextya, &
                tonextnextya,MPI_COMM_GRID,irecv_rq_fromnextnextya,mpierr)
          endif
!
          if (lastyb/=iproc) then
            call MPI_ISEND(fbo(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,lastyb, &
                tonextya,MPI_COMM_GRID,isend_rq_tolastyb,mpierr)
          endif
          if (nextya==iproc) then
            fahi(:,:,:,ivar1:ivar2)=fbo(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(fahi(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,nextya, &
                tonextya,MPI_COMM_GRID,irecv_rq_fromnextya,mpierr)
          endif
!
          if (nextyb/=iproc) then
            call MPI_ISEND(fbo(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,nextyb, &
                tolastya,MPI_COMM_GRID,isend_rq_tonextyb,mpierr)
          endif
          if (lastya==iproc) then
            falo(:,:,:,ivar1:ivar2)=fbo(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(falo(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,lastya, &
                tolastya,MPI_COMM_GRID,irecv_rq_fromlastya,mpierr)
          endif
!
          if (nextnextyb/=iproc) then
            call MPI_ISEND(fbo(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,nextnextyb, &
                tolastlastya,MPI_COMM_GRID,isend_rq_tonextnextyb,mpierr)
          endif
          if (lastlastya==iproc) then
            falolo(:,:,:,ivar1:ivar2)=fbo(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(falolo(:,:,:,ivar1:ivar2),nbufx_gh,mpi_precision,lastlastya, &
                tolastlastya,MPI_COMM_GRID,irecv_rq_fromlastlastya,mpierr)
          endif
!
        endif
      endif
!
    endsubroutine initiate_shearing
!***********************************************************************
    subroutine finalize_shearing(f,ivar1_opt,ivar2_opt)
!
!  Subroutine for shearing sheet boundary conditions
!
!  27-nov-14/mcnallcp: Now uses 4 shearing neighbours so
!                      y-ghosts are not needed
!  20-june-02/nils: adapted from pencil_mpi
!  02-mar-02/ulf: Sliding periodic boundary conditions in x
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
!     fa, fb are the buffers we collect recieved data into
      real, dimension (nghost,4*ny,mz,mcom) :: fa, fb
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fall, irecv_stat_fann
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fal,  irecv_stat_fan
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fbll, irecv_stat_fbnn
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fbl,  irecv_stat_fbn
      integer, dimension (MPI_STATUS_SIZE) :: isend_stat_tnna, isend_stat_tlla
      integer, dimension (MPI_STATUS_SIZE) :: isend_stat_tna,  isend_stat_tla
      integer, dimension (MPI_STATUS_SIZE) :: isend_stat_tnnb, isend_stat_tllb
      integer, dimension (MPI_STATUS_SIZE) :: isend_stat_tnb,  isend_stat_tlb
      integer :: ivar1, ivar2, m2long
      real(KIND=rkind8) :: deltay_dy, frac, c1, c2, c3, c4, c5, c6
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  Some special cases have already finished in initiate_shearing.
!
      if (nygrid/=1 .and. (nprocx>1 .or. nprocy>1) .and. (lfirst_proc_x .or. llast_proc_x)) then
!
!  Need to wait till all communication has been recived.
!
        if (lastlastyb/=iproc) call MPI_WAIT(irecv_rq_fromlastlastyb,irecv_stat_fbll,mpierr)
        if (lastyb/=iproc) call MPI_WAIT(irecv_rq_fromlastyb,irecv_stat_fbl,mpierr)
        if (nextyb/=iproc) call MPI_WAIT(irecv_rq_fromnextyb,irecv_stat_fbn,mpierr)
        if (nextnextyb/=iproc) call MPI_WAIT(irecv_rq_fromnextnextyb,irecv_stat_fbnn,mpierr)
!
        if (lastlastya/=iproc) call MPI_WAIT(irecv_rq_fromlastlastya,irecv_stat_fall,mpierr)
        if (lastya/=iproc) call MPI_WAIT(irecv_rq_fromlastya,irecv_stat_fal,mpierr)
        if (nextya/=iproc) call MPI_WAIT(irecv_rq_fromnextya,irecv_stat_fan,mpierr)
        if (nextnextya/=iproc) call MPI_WAIT(irecv_rq_fromnextnextya,irecv_stat_fann,mpierr)
!
!  Reading communicated information into f.
!
        deltay_dy=deltay/dy
!
        fa(:,1     :  ny,:,ivar1:ivar2) = falolo(:,:,:,ivar1:ivar2)
        fa(:,ny+1  :2*ny,:,ivar1:ivar2) =   falo(:,:,:,ivar1:ivar2)
        fa(:,2*ny+1:3*ny,:,ivar1:ivar2) =   fahi(:,:,:,ivar1:ivar2)
        fa(:,3*ny+1:4*ny,:,ivar1:ivar2) = fahihi(:,:,:,ivar1:ivar2)
!
        fb(:,1     :  ny,:,ivar1:ivar2) = fblolo(:,:,:,ivar1:ivar2)
        fb(:,ny+1  :2*ny,:,ivar1:ivar2) =   fblo(:,:,:,ivar1:ivar2)
        fb(:,2*ny+1:3*ny,:,ivar1:ivar2) =   fbhi(:,:,:,ivar1:ivar2)
        fb(:,3*ny+1:4*ny,:,ivar1:ivar2) = fbhihi(:,:,:,ivar1:ivar2)
!
!       displs is the displacement in cells
        displs = modulo(int(deltay_dy),ny)
        frac = deltay_dy - int(deltay_dy)
        c1 = -          (frac+1.)*frac*(frac-1.)*(frac-2.)*(frac-3.)/120.
        c2 = +(frac+2.)          *frac*(frac-1.)*(frac-2.)*(frac-3.)/24.
        c3 = -(frac+2.)*(frac+1.)     *(frac-1.)*(frac-2.)*(frac-3.)/12.
        c4 = +(frac+2.)*(frac+1.)*frac          *(frac-2.)*(frac-3.)/12.
        c5 = -(frac+2.)*(frac+1.)*frac*(frac-1.)          *(frac-3.)/24.
        c6 = +(frac+2.)*(frac+1.)*frac*(frac-1.)*(frac-2.)          /120.
!
!       These are handy for debugging, and reduce the interpolation order
!        c1=0; c2=0; c5=0; c6=0;
!        c3 = -(frac-1.)
!        c4 = frac
!
!       m2 long is the position of the last value from fahi in fa
!
        m2long = 3*ny
        f(1:l1-1,m1:m2,:,ivar1:ivar2) = &
             c1*fa(:,m2long-ny-displs+3:m2long-displs+2,:,ivar1:ivar2) &
            +c2*fa(:,m2long-ny-displs+2:m2long-displs+1,:,ivar1:ivar2) &
            +c3*fa(:,m2long-ny-displs+1:m2long-displs-0,:,ivar1:ivar2) &
            +c4*fa(:,m2long-ny-displs-0:m2long-displs-1,:,ivar1:ivar2) &
            +c5*fa(:,m2long-ny-displs-1:m2long-displs-2,:,ivar1:ivar2) &
            +c6*fa(:,m2long-ny-displs-2:m2long-displs-3,:,ivar1:ivar2)
!
!       ny+1 is the beginning of the block of interior cell from fblo in fb
        f(l2+1:mx,m1:m2,:,ivar1:ivar2)= &
             c1*fb(:,ny+1+displs-2:2*ny+displs-2,:,ivar1:ivar2) &
            +c2*fb(:,ny+1+displs-1:2*ny+displs-1,:,ivar1:ivar2) &
            +c3*fb(:,ny+1+displs  :2*ny+displs  ,:,ivar1:ivar2) &
            +c4*fb(:,ny+1+displs+1:2*ny+displs+1,:,ivar1:ivar2) &
            +c5*fb(:,ny+1+displs+2:2*ny+displs+2,:,ivar1:ivar2) &
            +c6*fb(:,ny+1+displs+3:2*ny+displs+3,:,ivar1:ivar2)
!
!  Need to wait till buffer is empty before re-using it again.
!
        if (nextnextyb/=iproc) call MPI_WAIT(isend_rq_tonextnextyb,isend_stat_tnnb,mpierr)
        if (nextyb/=iproc)     call MPI_WAIT(isend_rq_tonextyb,isend_stat_tnb,mpierr)
        if (lastyb/=iproc)     call MPI_WAIT(isend_rq_tolastyb,isend_stat_tlb,mpierr)
        if (lastlastyb/=iproc) call MPI_WAIT(isend_rq_tolastlastyb,isend_stat_tllb,mpierr)
!
        if (nextnextya/=iproc) call MPI_WAIT(isend_rq_tonextnextya,isend_stat_tnna,mpierr)
        if (nextya/=iproc)     call MPI_WAIT(isend_rq_tonextya,isend_stat_tna,mpierr)
        if (lastya/=iproc)     call MPI_WAIT(isend_rq_tolastya,isend_stat_tla,mpierr)
        if (lastlastya/=iproc) call MPI_WAIT(isend_rq_tolastlastya,isend_stat_tlla,mpierr)
!
      endif
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
      call MPI_RECV(Qrecv_zx,mx*mz,mpi_precision,isource,Qtag_zx+idir, &
                    MPI_COMM_GRID,irecv_zx,mpierr)
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
      call MPI_RECV(Qrecv_xy,mx*my,mpi_precision,isource,Qtag_xy+idir, &
                    MPI_COMM_GRID,irecv_xy,mpierr)
!
    endsubroutine radboundary_xy_recv
! Felipe: Added below
!***********************************************************************
    subroutine radboundary_yz_recv(lrad,idir,Qrecv_yz)
!
!  receive intensities from neighboring processor in z
!
!  11-jul-03/tobi: coded
!  20-jul-05/tobi: use non-blocking MPI calls
!  14-apr-22/felipe: adapted from radboundary_xy_recv
!
      integer :: lrad,idir
      real, dimension(my,mz) :: Qrecv_yz
      integer :: isource
      integer, dimension(MPI_STATUS_SIZE) :: irecv_yz
!
!  Identifier
!
      if (lroot.and.ip<5) print*,'radboundary_yz_recv: ENTER'
!
!  source
!
      if (lrad>0) isource=xlneigh
      if (lrad<0) isource=xuneigh
!
!  actual MPI call
!
      call MPI_RECV(Qrecv_yz,my*mz,mpi_precision,isource,Qtag_yz+idir, &
                    MPI_COMM_GRID,irecv_yz,mpierr)
!
    endsubroutine radboundary_yz_recv
! Felipe: Added above
!***********************************************************************
    subroutine radboundary_zx_send(mrad,idir,Qsend_zx)
!
!  send intensities to neighboring processor in y
!
!  11-jul-03/tobi: coded
!  20-jul-05/tobi: use non-blocking MPI calls
!   2-jun-15/MR: corrected parameters of MPI_SEND (blocking!)
!
      integer :: mrad,idir
      real, dimension(mx,mz) :: Qsend_zx
      integer :: idest
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
      call MPI_SEND(Qsend_zx,mx*mz,mpi_precision,idest,Qtag_zx+idir, &
                    MPI_COMM_GRID,mpierr)
!
    endsubroutine radboundary_zx_send
!***********************************************************************
    subroutine radboundary_xy_send(nrad,idir,Qsend_xy)
!
!  send intensities to neighboring processor in z
!
!  11-jul-03/tobi: coded
!  20-jul-05/tobi: use non-blocking MPI calls
!   2-jun-15/MR: corrected parameters of MPI_SEND (blocking!)
!
      integer, intent(in) :: nrad,idir
      real, dimension(mx,my) :: Qsend_xy
      integer :: idest
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
      call MPI_SEND(Qsend_xy,mx*my,mpi_precision,idest,Qtag_xy+idir, &
                    MPI_COMM_GRID,mpierr)
!
    endsubroutine radboundary_xy_send
!***********************************************************************
! Felipe: added below
!***********************************************************************
    subroutine radboundary_yz_send(lrad,idir,Qsend_yz)
!
!  send intensities to neighboring processor in x
!
!  11-jul-03/tobi: coded
!  20-jul-05/tobi: use non-blocking MPI calls
!   2-jun-15/MR: corrected parameters of MPI_SEND (blocking!)
!  14-apr-22/felipe: adapted from radboundary_xy_send
!
      integer, intent(in) :: lrad,idir
      real, dimension(my,mz) :: Qsend_yz
      integer :: idest
!
!  Identifier
!
      if (lroot.and.ip<5) print*,'radboundary_yz_send: ENTER'
!
!  destination
!
      if (lrad>0) idest=xuneigh
      if (lrad<0) idest=xlneigh
!
!  actual MPI call
!
      call MPI_SEND(Qsend_yz,my*mz,mpi_precision,idest,Qtag_yz+idir, &
                    MPI_COMM_GRID,mpierr)
!
    endsubroutine radboundary_yz_send
!***********************************************************************
! Felipe: added above
    subroutine radboundary_yz_sendrecv(lrad,idir,Qsend_yz,Qrecv_yz)
!
!  receive intensities from isource and send intensities to idest
!
!  17-nov-14/axel: adapted from radboundary_zx_sendrecv
!
      integer, intent(in) :: lrad,idir
      real, dimension(my,mz) :: Qsend_yz,Qrecv_yz
      integer :: idest,isource
      integer, dimension(MPI_STATUS_SIZE) :: isendrecv_yz
!
!  Identifier
!
      if (lroot.and.ip<5) print*,'radboundary_yz_sendrecv: ENTER'
!
!  destination and source id
!
      if (lrad>0) then; idest=xuneigh; isource=xlneigh; endif
      if (lrad<0) then; idest=xlneigh; isource=xuneigh; endif
!
!  actual MPI call
!
      call MPI_SENDRECV(Qsend_yz,my*mz,mpi_precision,idest,Qtag_yz+idir, &
                        Qrecv_yz,my*mz,mpi_precision,isource,Qtag_yz+idir, &
                        MPI_COMM_GRID,isendrecv_yz,mpierr)
!
    endsubroutine radboundary_yz_sendrecv
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
      call MPI_SENDRECV(Qsend_zx,mx*mz,mpi_precision,idest,Qtag_zx+idir, &
                        Qrecv_zx,mx*mz,mpi_precision,isource,Qtag_zx+idir, &
                        MPI_COMM_GRID,isendrecv_zx,mpierr)
!
    endsubroutine radboundary_zx_sendrecv
!***********************************************************************
    subroutine radboundary_yz_periodic_ray(Qrad_yz,tau_yz, &
                                           Qrad_yz_all,tau_yz_all)
!
!  Gather all intrinsic optical depths and heating rates into one rank-3 array
!  that is available on each processor.
!
!  17-nov-14/axel: adapted from radboundary_zx_periodic_ray
!
      real, dimension(ny,nz), intent(in) :: Qrad_yz,tau_yz
      real, dimension(ny,nz,0:nprocx-1), intent(out) :: Qrad_yz_all,tau_yz_all
!
!  Identifier
!
      if (lroot.and.ip<5) print*,'radboundary_yz_periodic_ray: ENTER'
!
!  actual MPI calls
!
      call MPI_ALLGATHER(tau_yz,ny*nz,mpi_precision,tau_yz_all,ny*nz,mpi_precision, &
                         MPI_COMM_XBEAM,mpierr)
!
      call MPI_ALLGATHER(Qrad_yz,ny*nz,mpi_precision,Qrad_yz_all,ny*nz,mpi_precision, &
                         MPI_COMM_XBEAM,mpierr)
!
    endsubroutine radboundary_yz_periodic_ray
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
      call MPI_ALLGATHER(tau_zx,nx*nz,mpi_precision,tau_zx_all,nx*nz,mpi_precision, &
          MPI_COMM_YBEAM,mpierr)
!
      call MPI_ALLGATHER(Qrad_zx,nx*nz,mpi_precision,Qrad_zx_all,nx*nz,mpi_precision, &
          MPI_COMM_YBEAM,mpierr)
!
    endsubroutine radboundary_zx_periodic_ray
!***********************************************************************
    subroutine mpirecv_char_scl(str,proc_src,tag_id,comm)
!
!  Receive character scalar from other processor.
!
!  04-sep-06/wlad: coded
!
      character(LEN=*), intent(OUT) :: str
      integer,          intent(IN ) :: proc_src, tag_id
      integer, optional,intent(IN)  :: comm

      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      call MPI_RECV(str, len(str), MPI_CHARACTER, proc_src, &
                    tag_id, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
!
    endsubroutine mpirecv_char_scl
!***********************************************************************
    subroutine mpisend_char_scl(str,proc_src,tag_id,comm)
!
!  Send character scalar to other processor.
!
!  04-sep-06/wlad: coded
!
      character(LEN=*) :: str
      integer :: proc_src, tag_id
      integer, optional :: comm
!
      call MPI_SEND(str, len(str), MPI_CHARACTER, proc_src, &
                    tag_id, ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpisend_char_scl
!***********************************************************************
    subroutine mpirecv_logical_scl(bcast_array,proc_src,tag_id,comm)
!
!  Receive logical scalar from other processor.
!
!  04-sep-06/wlad: coded
!
      logical :: bcast_array
      integer :: proc_src, tag_id
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer, optional :: comm
!
      call MPI_RECV(bcast_array, 1, MPI_LOGICAL, proc_src, &
                    tag_id, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
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
      if (nbcast_array == 0) return
!
      call MPI_RECV(bcast_array, nbcast_array, MPI_LOGICAL, proc_src, &
                    tag_id, MPI_COMM_GRID, stat, mpierr)
!
    endsubroutine mpirecv_logical_arr
!***********************************************************************
    subroutine mpirecv_real_scl(bcast_array,proc_src,tag_id,comm,nonblock)
!
!  Receive real scalar from other processor.
!  If present, nonblock is the request id for non-blockin receive.
!
!  02-jul-05/anders: coded
!
      real :: bcast_array
      integer :: proc_src, tag_id
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
!
      integer, optional :: comm, nonblock

      if (present(nonblock)) then
        call MPI_IRECV(bcast_array, 1, mpi_precision, proc_src, &
                      tag_id, ioptest(comm,MPI_COMM_GRID), stat, nonblock, mpierr)
      else
        call MPI_RECV(bcast_array, 1, mpi_precision, proc_src, &
                      tag_id, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
 !
    endsubroutine mpirecv_real_scl
!***********************************************************************
    subroutine mpirecv_real_arr(bcast_array,nbcast_array,proc_src,tag_id,comm,nonblock)
!
!  Receive real array from other processor.
!  If present, nonblock is the request id for non-blockin receive.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array
      integer :: proc_src, tag_id
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer, optional :: comm, nonblock
!
      intent(out) :: bcast_array
!
      if (nbcast_array == 0) return
!
      if (present(nonblock)) then
        call MPI_IRECV(bcast_array, nbcast_array, mpi_precision, proc_src, &
                       tag_id, ioptest(comm,MPI_COMM_GRID), nonblock, mpierr)
      else
        call MPI_RECV(bcast_array, nbcast_array, mpi_precision, proc_src, &
                      tag_id, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine mpirecv_real_arr
!***********************************************************************
    subroutine mpirecv_real_arr2(bcast_array,nbcast_array,proc_src,tag_id,comm)
!
!  Receive real array(:,:) from other processor.
!
!  02-jul-05/anders: coded
!
      integer, dimension(2) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer :: proc_src, tag_id, num_elements
      integer, optional :: comm
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_RECV(bcast_array, num_elements, mpi_precision, proc_src, &
                    tag_id, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
!
    endsubroutine mpirecv_real_arr2
!***********************************************************************
    subroutine mpirecv_real_arr3(bcast_array,nbcast_array,proc_src,tag_id,comm,nonblock)
!
!  Receive real array(:,:,:) from other processor.
!  If present, nonblock is the request id for non-blockin receive.
!
!  20-may-06/anders: adapted
!
      integer, dimension(3) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3)) :: bcast_array
      integer :: proc_src, tag_id, num_elements
      integer, optional :: comm,nonblock
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      if (present(nonblock)) then
        call MPI_IRECV(bcast_array, num_elements, mpi_precision, proc_src, &
                       tag_id, ioptest(comm,MPI_COMM_GRID), nonblock, mpierr)
      else
        call MPI_RECV(bcast_array, num_elements, mpi_precision, proc_src, &
                      tag_id, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine mpirecv_real_arr3
!***********************************************************************
    subroutine mpirecv_cmplx_arr3(bcast_array,nbcast_array,proc_src,tag_id,comm,nonblock)
!
!  Receive complex array(:,:,:) from other processor.
!  If present, nonblock is the request id for non-blockin receive.
!
!  20-may-06/anders: adapted
!
      integer, dimension(3) :: nbcast_array
      complex, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3)) :: bcast_array
      integer :: proc_src, tag_id, num_elements
      integer, optional :: comm,nonblock
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      if (present(nonblock)) then
        call MPI_IRECV(bcast_array, num_elements, MPI_CMPLX, proc_src, &
                       tag_id, ioptest(comm,MPI_COMM_GRID), nonblock, mpierr)
      else
        call MPI_RECV(bcast_array, num_elements, MPI_CMPLX, proc_src, &
                      tag_id, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine mpirecv_cmplx_arr3
!***********************************************************************
    subroutine mpirecv_real_arr4(bcast_array,nbcast_array,proc_src,tag_id,comm,nonblock)
!
!  Receive real array(:,:,:,:) from other processor.
!  If present, nonblock is the request id for non-blockin receive.
!
!  20-may-06/anders: adapted
!
      integer, dimension(4) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3),nbcast_array(4)) :: bcast_array
      integer :: proc_src, tag_id, num_elements
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer, optional :: comm,nonblock
!
      intent(out) :: bcast_array
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      if (present(nonblock)) then
        call MPI_IRECV(bcast_array, num_elements, mpi_precision, proc_src, &
                       tag_id, ioptest(comm,MPI_COMM_GRID), nonblock, mpierr)
      else
        call MPI_RECV(bcast_array, num_elements, mpi_precision, proc_src, &
                      tag_id, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine mpirecv_real_arr4
!***********************************************************************
    subroutine mpirecv_real_arr5(bcast_array,nbcast_array,proc_src,tag_id)
!
!  Receive real array(:,:,:,:) from other processor.
!
!  24-apr-17/Jorgen: adapted
!
      integer, dimension(5) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3),nbcast_array(4),nbcast_array(5)) :: bcast_array
      integer :: proc_src, tag_id, num_elements
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_RECV(bcast_array, num_elements, mpi_precision, proc_src, &
                    tag_id, MPI_COMM_GRID, stat, mpierr)
!
    endsubroutine mpirecv_real_arr5
!***********************************************************************
    subroutine mpirecv_int_scl(bcast_array,proc_src,tag_id,comm,nonblock)
!
!  Receive integer scalar from other processor.
!  If present, nonblock is the request id for non-blockin receive.
!
!  02-jul-05/anders: coded
!
      integer :: bcast_array
      integer :: proc_src, tag_id
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer, optional :: comm, nonblock
!
      if (present(nonblock)) then
        call MPI_IRECV(bcast_array, 1, MPI_INTEGER, proc_src, &
                      tag_id, ioptest(comm,MPI_COMM_GRID), nonblock, mpierr)
      else
        call MPI_RECV(bcast_array, 1, MPI_INTEGER, proc_src, &
                      tag_id, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine mpirecv_int_scl
!***********************************************************************
    subroutine mpirecv_int_arr(bcast_array,nbcast_array,proc_src,tag_id,comm,nonblock)
!
!  Receive integer array from other processor.
!  If present, nonblock is the request id for non-blockin receive.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: bcast_array
      integer :: proc_src, tag_id
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer, optional :: comm, nonblock
!
      if (nbcast_array == 0) return
!
      if (present(nonblock)) then
        call MPI_IRECV(bcast_array, nbcast_array, MPI_INTEGER, proc_src, &
                       tag_id, ioptest(comm,MPI_COMM_GRID), nonblock, mpierr)
      else
        call MPI_RECV(bcast_array, nbcast_array, MPI_INTEGER, proc_src, &
                      tag_id, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
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
      integer :: proc_src, tag_id, num_elements
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      intent(out) :: bcast_array
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_RECV(bcast_array, num_elements, MPI_INTEGER, proc_src, &
                    tag_id, MPI_COMM_GRID, stat, mpierr)
!
    endsubroutine mpirecv_int_arr2
!***********************************************************************
    subroutine mpisend_logical_scl(bcast_array,proc_rec,tag_id,comm)
!
!  Send logical scalar to other processor.
!
!  04-sep-06/wlad: coded
!
      logical :: bcast_array
      integer :: proc_rec, tag_id
      integer, optional :: comm
!
      call MPI_SEND(bcast_array, 1, MPI_LOGICAL, proc_rec, &
                    tag_id, ioptest(comm,MPI_COMM_GRID), mpierr)
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
      if (nbcast_array == 0) return
!
      call MPI_SEND(bcast_array, nbcast_array, MPI_LOGICAL, proc_rec, &
                    tag_id, MPI_COMM_GRID, mpierr)
!
    endsubroutine mpisend_logical_arr
!***********************************************************************
    subroutine mpisend_real_scl(bcast_array,proc_rec,tag_id,comm)
!
!  Send real scalar to other processor.
!
!  02-jul-05/anders: coded
!
      real :: bcast_array
      integer :: proc_rec, tag_id
      integer, optional :: comm
!
      call MPI_SEND(bcast_array, 1, mpi_precision, proc_rec, &
                    tag_id, ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpisend_real_scl
!***********************************************************************
    subroutine mpisend_real_arr(bcast_array,nbcast_array,proc_rec,tag_id,comm)
!
!  Send real array to other processor.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array
      integer :: proc_rec, tag_id
      integer, optional :: comm
!
      if (nbcast_array == 0) return
!
      call MPI_SEND(bcast_array, nbcast_array, mpi_precision, proc_rec, &
                    tag_id, ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpisend_real_arr
!***********************************************************************
    subroutine mpisend_real_arr_assumed(bcast_array,nbcast_array,offset,proc_rec,tag_id,comm)
!
!  Sends nbcast_array elements of real array bcast_array from position offset to other processor.
!  Avoids compilation error when shapes of bcast_array and actual parameter do not agree.
!
!  06-oct-22/MR: coded
!
      real, dimension(*) :: bcast_array
      integer(KIND=ikind8) :: offset
      integer :: nbcast_array,proc_rec,tag_id
      integer, optional :: comm
!
      if (nbcast_array == 0) return
!
      call MPI_SEND(bcast_array(offset), nbcast_array, mpi_precision, proc_rec, &
                    tag_id, ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpisend_real_arr_assumed
!***********************************************************************
    subroutine mpisend_real_arr_huge(array,len_array,partner,tag,comm)
!
!  Allows to communicate arrays with length > max_int.
!  06-oct-22/MR: coded
!
      real, dimension(*) :: array
      integer(KIND=ikind8) :: len_array
      integer :: partner,tag
      integer, optional :: comm
!
      integer :: mult,res

      if (len_array == 0) return

      mult=len_array/max_int
      res =len_array-mult*max_int

      if (mult>0) &
        call MPI_SEND (array, mult, REAL_ARR_MAXSIZE, partner, tag, MPI_COMM_GRID, mpierr)
      if (res>0) &
        call MPI_SEND(array(len_array-res+1), res, mpi_precision, partner, tag+1, &
                     ioptest(comm,MPI_COMM_GRID), mpierr)

    endsubroutine mpisend_real_arr_huge
!***********************************************************************
    subroutine mpirecv_real_arr_huge(array,len_array,partner,tag,comm)
!
!  Allows to communicate arrays with length > max_int.
!  06-oct-22/MR: coded
!
      real, dimension(*) :: array
      integer(KIND=ikind8) :: len_array
      integer :: partner,tag
      integer, optional :: comm
!
      integer :: mult,res
      integer, dimension(MPI_STATUS_SIZE) :: stat

      if (len_array == 0) return

      mult=len_array/max_int
      res =len_array-mult*max_int

      if (mult>0) &
        call MPI_RECV(array, mult, REAL_ARR_MAXSIZE, partner, tag, MPI_COMM_GRID, stat, mpierr)
      if (res>0) &
        call MPI_RECV(array(len_array-res+1), res, mpi_precision, partner, tag+1, &
                     ioptest(comm,MPI_COMM_GRID),stat,mpierr)

    endsubroutine mpirecv_real_arr_huge
!***********************************************************************
    subroutine mpirecv_real_arr_assumed(bcast_array,nbcast_array,offset,proc_rec,tag_id,comm)
!
!  Receives nbcast_array elements of real array bcast_array at position offset from other processor.
!  Avoids compilation error when shapes of bcast_array and actual parameter do not agree.
!
!  06-oct-22/MR: coded
!
      real, dimension(*) :: bcast_array
      integer(KIND=ikind8) :: offset
      integer :: nbcast_array,proc_rec,tag_id
      integer, optional :: comm
!
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      if (nbcast_array == 0) return
!
      call MPI_RECV(bcast_array(offset), nbcast_array, mpi_precision, proc_rec, &
                    tag_id, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
!
    endsubroutine mpirecv_real_arr_assumed
!***********************************************************************
    subroutine mpisend_real_arr2(bcast_array,nbcast_array,proc_rec,tag_id,comm)
!
!  Send real array(:,:) to other processor.
!
!  02-jul-05/anders: coded
!
      integer, dimension(2) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer :: proc_rec, tag_id, num_elements
      integer, optional :: comm
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_SEND(bcast_array, num_elements, mpi_precision, proc_rec, &
                    tag_id, ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpisend_real_arr2
!***********************************************************************
    subroutine mpisend_real_arr3(bcast_array,nbcast_array,proc_rec,tag_id,comm,nonblock)
!
!  Send real array(:,:,:) to other processor.
!  If present, nonblock is the request id for non-blockin send.
!
!  20-may-06/anders: adapted
!
      integer, dimension(3) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3)) :: bcast_array
      integer :: proc_rec, tag_id, num_elements
      integer, optional :: comm,nonblock
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      if (present(nonblock)) then
        call MPI_ISEND(bcast_array, num_elements, mpi_precision, proc_rec, &
                       tag_id,ioptest(comm,MPI_COMM_GRID),nonblock, mpierr)
      else
        call MPI_SEND(bcast_array, num_elements, mpi_precision, proc_rec, &
                      tag_id,ioptest(comm,MPI_COMM_GRID),mpierr)
      endif
!
    endsubroutine mpisend_real_arr3
!***********************************************************************
    subroutine mpisend_cmplx_arr3(bcast_array,nbcast_array,proc_rec,tag_id,comm,nonblock)
!
!  Send real array(:,:,:) to other processor.
!  If present, nonblock is the request id for non-blockin send.
!
!  20-may-06/anders: adapted
!
      integer, dimension(3) :: nbcast_array
      complex, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3)) :: bcast_array
      integer :: proc_rec, tag_id, num_elements
      integer, optional :: comm,nonblock
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      if (present(nonblock)) then
        call MPI_ISEND(bcast_array, num_elements, MPI_CMPLX, proc_rec, &
                       tag_id,ioptest(comm,MPI_COMM_GRID),nonblock, mpierr)
      else
        call MPI_SEND(bcast_array, num_elements, MPI_CMPLX, proc_rec, &
                      tag_id,ioptest(comm,MPI_COMM_GRID),mpierr)
      endif
!
    endsubroutine mpisend_cmplx_arr3
!***********************************************************************
    subroutine mpisend_real_arr4(bcast_array,nbcast_array,proc_rec,tag_id,comm)
!
!  Send real array(:,:,:,:) to other processor.
!
!  20-may-06/anders: adapted
!
      integer, dimension(4) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3),nbcast_array(4)) :: bcast_array
      integer :: proc_rec, tag_id, num_elements
      integer, optional :: comm
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_SEND(bcast_array, num_elements, mpi_precision, proc_rec, &
                    tag_id, ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpisend_real_arr4
!***********************************************************************
    subroutine mpisend_real_arr5(bcast_array,nbcast_array,proc_rec,tag_id)
!
!  Send real array(:,:,:,:) to other processor.
!
!  20-may-06/anders: adapted
!
      integer, dimension(5) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3),nbcast_array(4),nbcast_array(5)) :: bcast_array
      integer :: proc_rec, tag_id, num_elements
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_SEND(bcast_array, num_elements, mpi_precision, proc_rec, &
                    tag_id, MPI_COMM_GRID,mpierr)
!
    endsubroutine mpisend_real_arr5
!***********************************************************************
    subroutine mpisendrecv_int_arr(send_array,sendcnt,proc_dest,sendtag, &
                                   recv_array,proc_src,recvtag,comm)
      integer :: sendcnt
      integer, dimension(sendcnt) :: send_array, recv_array
      integer :: proc_src, proc_dest, sendtag, recvtag
      integer, optional :: comm
      integer, dimension(MPI_STATUS_SIZE) :: stat

      intent(out) :: recv_array
!
      if (sendcnt==0) return

      call MPI_SENDRECV(send_array,sendcnt,MPI_INTEGER,proc_dest,sendtag, &
                        recv_array,sendcnt,MPI_INTEGER,proc_src,recvtag, &
                        ioptest(comm,MPI_COMM_GRID),stat,mpierr)
!
    endsubroutine mpisendrecv_int_arr
!***********************************************************************
    subroutine mpisendrecv_real_scl(send_array,proc_dest,sendtag, &
                                    recv_array,proc_src,recvtag)

      real :: send_array, recv_array
      integer :: proc_src, proc_dest, sendtag, recvtag
      integer, dimension(MPI_STATUS_SIZE) :: stat

      intent(out) :: recv_array

      call MPI_SENDRECV(send_array,1,mpi_precision,proc_dest,sendtag, &
                        recv_array,1,mpi_precision,proc_src,recvtag, &
                        MPI_COMM_GRID,stat,mpierr)

    endsubroutine mpisendrecv_real_scl
!***********************************************************************
    subroutine mpisendrecv_real_arr(send_array,sendcnt,proc_dest,sendtag, &
                                    recv_array,proc_src,recvtag,idir)

      integer :: sendcnt
      real, dimension(sendcnt) :: send_array
      real, dimension(sendcnt) :: recv_array
      integer :: proc_src, proc_dest, sendtag, recvtag
      integer, optional :: idir

      integer, dimension(MPI_STATUS_SIZE) :: stat

      intent(out) :: recv_array
!
      if (sendcnt==0) return

      call MPI_SENDRECV(send_array,sendcnt,mpi_precision,proc_dest,sendtag, &
                        recv_array,sendcnt,mpi_precision,proc_src,recvtag, &
                        mpigetcomm(ioptest(idir)),stat,mpierr)
!
    endsubroutine mpisendrecv_real_arr
!***********************************************************************
    subroutine mpisendrecv_real_arr2(send_array,nbcast_array,proc_dest,sendtag, &
                                     recv_array,proc_src,recvtag,idir)

      integer, dimension(2) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2)) :: send_array
      real, dimension(nbcast_array(1),nbcast_array(2)) :: recv_array
      integer :: proc_src, proc_dest, sendtag, recvtag, num_elements
      integer, dimension(MPI_STATUS_SIZE) :: stat
      intent(out) :: recv_array
      integer, optional :: idir

      if (any(nbcast_array == 0)) return

      num_elements = product(nbcast_array)
      call MPI_SENDRECV(send_array,num_elements,mpi_precision,proc_dest,sendtag, &
                        recv_array,num_elements,mpi_precision,proc_src,recvtag, &
                        mpigetcomm(ioptest(idir)),stat,mpierr)

    endsubroutine mpisendrecv_real_arr2
!***********************************************************************
    subroutine mpisendrecv_real_arr3(send_array,nbcast_array,proc_dest,sendtag, &
                                     recv_array,proc_src,recvtag)

     integer, dimension(3) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3)) :: send_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3)) :: recv_array
      integer :: proc_src, proc_dest, sendtag, recvtag, num_elements
      integer, dimension(MPI_STATUS_SIZE) :: stat
      intent(out) :: recv_array

      if (any(nbcast_array == 0)) return

      num_elements = product(nbcast_array)
      call MPI_SENDRECV(send_array,num_elements,mpi_precision,proc_dest,sendtag, &
                        recv_array,num_elements,mpi_precision,proc_src,recvtag, &
                        MPI_COMM_GRID,stat,mpierr)

    endsubroutine mpisendrecv_real_arr3
!***********************************************************************
    subroutine mpisendrecv_real_arr4(send_array,nbcast_array,proc_dest,sendtag, &
                                     recv_array,proc_src,recvtag)

      integer, dimension(4) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3),nbcast_array(4)) :: send_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3),nbcast_array(4)) :: recv_array
      integer :: proc_src, proc_dest, sendtag, recvtag, num_elements
      integer, dimension(MPI_STATUS_SIZE) :: stat
      intent(out) :: recv_array

      if (any(nbcast_array == 0)) return

      num_elements = product(nbcast_array)
      call MPI_SENDRECV(send_array,num_elements,mpi_precision,proc_dest,sendtag, &
                        recv_array,num_elements,mpi_precision,proc_src,recvtag, &
                        MPI_COMM_GRID,stat,mpierr)

    endsubroutine mpisendrecv_real_arr4
!***********************************************************************
    subroutine mpisend_int_scl(bcast_array,proc_rec,tag_id,comm)
!
!  Send integer scalar to other processor.
!
!  02-jul-05/anders: coded
!
      integer :: bcast_array
      integer :: proc_rec, tag_id
      integer, optional :: comm
!
      call MPI_SEND(bcast_array, 1, MPI_INTEGER, proc_rec, &
                    tag_id, ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpisend_int_scl
!***********************************************************************
    subroutine mpirecv_nonblock_int_scl(bcast_array,proc_src,tag_id,ireq,comm)
!
!  Receive integer scalar from other processor, with non-blocking communication.
!
!  12-dec-14/wlad: adapted
!
      integer :: bcast_array
      integer :: proc_src, tag_id, ireq
      integer, optional :: comm
!
      call MPI_IRECV(bcast_array, 1, MPI_INTEGER, proc_src, &
                     tag_id, ioptest(comm,MPI_COMM_GRID), ireq, mpierr)
!
    endsubroutine mpirecv_nonblock_int_scl
!***********************************************************************
    subroutine mpirecv_nonblock_int_arr(bcast_array,nbcast_array,proc_src,tag_id,ireq)
!
!  Receive integer array from other processor, with non-blocking communication.
!
!  12-dec-14/wlad: adapted
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: bcast_array
      integer :: proc_src, tag_id, ireq

      if (nbcast_array == 0) return
!
      call MPI_IRECV(bcast_array, nbcast_array, MPI_INTEGER, proc_src, &
                     tag_id, MPI_COMM_GRID, ireq, mpierr)
!
    endsubroutine mpirecv_nonblock_int_arr
!***********************************************************************
    subroutine mpirecv_nonblock_int_arr2(bcast_array,nbcast_array,proc_src,tag_id,ireq)
!
!  Receive integer array(:,:) from other processor, with non-blocking communication.
!
!  30-apr-17/Jorgen: adapted
!
      integer, dimension(2) :: nbcast_array
      integer, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer :: proc_src, tag_id, ireq, num_elements

      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_IRECV(bcast_array, nbcast_array, MPI_INTEGER, proc_src, &
                     tag_id, MPI_COMM_GRID, ireq, mpierr)
!
    endsubroutine mpirecv_nonblock_int_arr2
!***********************************************************************
    subroutine mpirecv_nonblock_real_arr(bcast_array,nbcast_array,proc_src,tag_id,ireq)
!
!  Receive real array from other processor, with non-blocking communication.
!
!  12-dec-14/wlad: adapted
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array
      integer :: proc_src, tag_id, ireq
!
      intent(out) :: bcast_array
!
      if (nbcast_array == 0) return
!
      call MPI_IRECV(bcast_array, nbcast_array, mpi_precision, proc_src, &
                     tag_id, MPI_COMM_GRID, ireq, mpierr)
!
    endsubroutine mpirecv_nonblock_real_arr
!***********************************************************************
    subroutine mpirecv_nonblock_real_arr2(bcast_array,nbcast_array,proc_src,tag_id,ireq)
!
!  Receive real array(:,:) from other processor, with non-blocking communication.
!
!  12-dec-14/wlad: adapted
!
      integer, dimension(2) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer :: proc_src, tag_id, ireq, num_elements
!
      intent(out) :: bcast_array

      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_IRECV(bcast_array, num_elements, mpi_precision, proc_src, &
                     tag_id, MPI_COMM_GRID, ireq, mpierr)
!
    endsubroutine mpirecv_nonblock_real_arr2
!***********************************************************************
    subroutine mpirecv_nonblock_real_arr3(bcast_array,nbcast_array,proc_src,tag_id,ireq,comm)
!
!  Receive real array(:,:) from other processor, with non-blocking communication.
!
!  07-jul-17/Jorgen: adapted
!
      integer, dimension(3) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3)) :: bcast_array
      integer :: proc_src, tag_id, ireq, num_elements
      integer, optional :: comm
!
      intent(out) :: bcast_array

      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_IRECV(bcast_array, num_elements, mpi_precision, proc_src, &
                     tag_id, ioptest(comm,MPI_COMM_GRID), ireq, mpierr)
!
    endsubroutine mpirecv_nonblock_real_arr3
!***********************************************************************
    subroutine mpirecv_nonblock_real_arr4(bcast_array,nbcast_array,proc_src,tag_id,ireq,comm)
!
!  Receive real array(:,:,:,:) from other processor, with non-blocking communication.
!
!  12-dec-14/wlad: adapted
!
      integer, dimension(4) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3),nbcast_array(4)) :: bcast_array
      integer :: proc_src, tag_id, ireq, num_elements
      integer, optional :: comm
!
      intent(out) :: bcast_array

      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_IRECV(bcast_array, num_elements, mpi_precision, proc_src, &
                     tag_id, ioptest(comm,MPI_COMM_GRID), ireq, mpierr)
!
    endsubroutine mpirecv_nonblock_real_arr4
!***********************************************************************
    subroutine mpirecv_nonblock_real_arr5(bcast_array,nbcast_array,proc_src,tag_id,ireq)
!
!  Receive real array(:,:,:,:) from other processor, with non-blocking communication.
!
!  24-apr-17/Jorgen: adapted
!
      integer, dimension(5) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3),nbcast_array(4),nbcast_array(5)) :: bcast_array
      integer :: proc_src, tag_id, ireq, num_elements
!
      intent(out) :: bcast_array

      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_IRECV(bcast_array, num_elements, mpi_precision, proc_src, &
                     tag_id, MPI_COMM_GRID, ireq, mpierr)
!
    endsubroutine mpirecv_nonblock_real_arr5
!***********************************************************************
    subroutine mpisend_nonblock_real_arr(bcast_array,nbcast_array,proc_rec,tag_id,ireq)
!
!  Send real array to other processor, with non-blocking communication.
!
!  12-dec-14/wlad: adapted
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array
      integer :: proc_rec, tag_id, ireq
!
      if (nbcast_array == 0) return
!
      call MPI_ISEND(bcast_array, nbcast_array, mpi_precision, proc_rec, &
                     tag_id, MPI_COMM_GRID, ireq, mpierr)
!
    endsubroutine mpisend_nonblock_real_arr
!***********************************************************************
    subroutine mpisend_nonblock_real_arr2(bcast_array,nbcast_array,proc_rec,tag_id,ireq)
!
!  Send real array(:,:,:) to other processor, with non-blocking communication.
!
!  03-okt-17/Jorgen: adapted
!
      integer, dimension(2) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer :: proc_rec, tag_id, ireq, num_elements
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_ISEND(bcast_array, num_elements, mpi_precision, proc_rec, &
                     tag_id, MPI_COMM_GRID,ireq,mpierr)
!
    endsubroutine mpisend_nonblock_real_arr2
!***********************************************************************
    subroutine mpisend_nonblock_real_arr3(bcast_array,nbcast_array,proc_rec,tag_id,ireq)
!
!  Send real array(:,:,:) to other processor, with non-blocking communication.
!
!  03-aug-17/Jorgen: adapted
!
      integer, dimension(3) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3)) :: bcast_array
      integer :: proc_rec, tag_id, ireq, num_elements
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_ISEND(bcast_array, num_elements, mpi_precision, proc_rec, &
                     tag_id, MPI_COMM_GRID,ireq,mpierr)
!
    endsubroutine mpisend_nonblock_real_arr3
!***********************************************************************
    subroutine mpisend_nonblock_real_arr4(bcast_array,nbcast_array,proc_rec,tag_id,ireq)
!
!  Send real array(:,:,:,:) to other processor, with non-blocking communication.
!
!  12-dec-14/wlad: adapted
!
      integer, dimension(4) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3),nbcast_array(4)) :: bcast_array
      integer :: proc_rec, tag_id, ireq, num_elements
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_ISEND(bcast_array, num_elements, mpi_precision, proc_rec, &
                     tag_id, MPI_COMM_GRID,ireq,mpierr)
!
    endsubroutine mpisend_nonblock_real_arr4
!***********************************************************************
    subroutine mpisend_nonblock_real_arr5(bcast_array,nbcast_array,proc_rec,tag_id,ireq)
!
!  Send real array(:,:,:,:) to other processor, with non-blocking communication.
!
!  24-apr-17/Jorgen: adapted
!
      integer, dimension(5) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3),nbcast_array(4),nbcast_array(5)) :: bcast_array
      integer :: proc_rec, tag_id, ireq, num_elements
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_ISEND(bcast_array, num_elements, mpi_precision, proc_rec, &
                     tag_id, MPI_COMM_GRID,ireq,mpierr)
!
    endsubroutine mpisend_nonblock_real_arr5
!***********************************************************************
    subroutine mpisend_nonblock_int_scl(bcast_array,proc_rec,tag_id,ireq)
!
!  Send integer scalar to other processor, with non-blocking communication.
!
!  12-dec-14/wlad: adapted
!
      integer :: bcast_array
      integer :: proc_rec, tag_id, ireq
!
      call MPI_ISEND(bcast_array, 1, MPI_INTEGER, proc_rec, &
                     tag_id, MPI_COMM_GRID, ireq, mpierr)
!
    endsubroutine mpisend_nonblock_int_scl
!***********************************************************************
    subroutine mpisend_nonblock_int_arr(bcast_array,nbcast_array,proc_rec,tag_id,ireq)
!
!  Send integer array to other processor, with non-blocking communication.
!
!  12-dec-14/wlad: adapted
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: bcast_array
      integer :: proc_rec, tag_id, ireq
!
      if (nbcast_array == 0) return
!
      call MPI_ISEND(bcast_array, nbcast_array, MPI_INTEGER, proc_rec, &
          tag_id, MPI_COMM_GRID, ireq, mpierr)
!
    endsubroutine mpisend_nonblock_int_arr
!***********************************************************************
    subroutine mpisend_nonblock_int_arr2(bcast_array,nbcast_array,proc_rec,tag_id,ireq)
!
!  Send 2d integer array to other processor, with non-blocking communication.
!
!  18-apr-17/Jorgen: adapted
!
      integer, dimension(2) :: nbcast_array
      integer, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer :: proc_rec, tag_id, ireq, num_elements
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_ISEND(bcast_array, num_elements, MPI_INTEGER, proc_rec, &
          tag_id, MPI_COMM_GRID, ireq, mpierr)
!
    endsubroutine mpisend_nonblock_int_arr2
!***********************************************************************
    subroutine mpisend_int_arr(bcast_array,nbcast_array,proc_rec,tag_id,comm,nonblock)
!
!  Send integer array to other processor.
!  If present, nonblock is the request id for non-blockin send.
!
!  02-jul-05/anders: coded
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: bcast_array
      integer :: proc_rec, tag_id
      integer, optional :: comm,nonblock
!
      if (nbcast_array == 0) return
!
      if (present(nonblock)) then
        call MPI_ISEND(bcast_array, nbcast_array, MPI_INTEGER, proc_rec, &
                       tag_id, ioptest(comm,MPI_COMM_GRID), nonblock, mpierr)
      else
        call MPI_SEND(bcast_array, nbcast_array, MPI_INTEGER, proc_rec, &
                      tag_id, ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
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
      integer :: proc_rec, tag_id, num_elements
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
      call MPI_SEND(bcast_array, num_elements, MPI_INTEGER, proc_rec, &
                    tag_id, MPI_COMM_GRID,mpierr)
!
    endsubroutine mpisend_int_arr2
!***********************************************************************
    subroutine mpibcast_logical_scl(lbcast_array,proc,comm)
!
!  Communicate logical scalar between processors.
!
      logical :: lbcast_array
      integer, optional :: proc,comm
!
      call MPI_BCAST(lbcast_array,1,MPI_LOGICAL,ioptest(proc,root), &
                     ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_logical_scl
!***********************************************************************
    subroutine mpibcast_logical_arr(lbcast_array,nbcast_array,proc,comm)
!
!  Communicate logical array between processors.
!
      integer :: nbcast_array
      logical, dimension (nbcast_array) :: lbcast_array
      integer, optional :: proc,comm
!
      if (nbcast_array == 0) return
!
      call MPI_BCAST(lbcast_array,nbcast_array,MPI_LOGICAL,ioptest(proc,root), &
                     ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_logical_arr
!***********************************************************************
    subroutine mpibcast_logical_arr2(lbcast_array,nbcast_array,proc,comm)
!
!  Communicate logical array(:,:) to other processor.
!
!  25-may-08/wlad: adapted
!
      integer, dimension(2) :: nbcast_array
      logical, dimension(nbcast_array(1),nbcast_array(2)) :: lbcast_array
      integer, optional :: proc,comm
!
      if (any(nbcast_array == 0)) return
!
      call MPI_BCAST(lbcast_array, product(nbcast_array), MPI_LOGICAL, ioptest(proc,root), &
                     ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_logical_arr2
!***********************************************************************
    subroutine mpibcast_int_scl(ibcast_array,proc,comm)
!
!  Communicate integer scalar between processors.
!
      integer :: ibcast_array
      integer, optional :: proc, comm
      call MPI_BCAST(ibcast_array,1,MPI_INTEGER,ioptest(proc,root), &
                     ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_int_scl
!***********************************************************************
    subroutine mpibcast_int_arr(ibcast_array,nbcast_array,proc,comm)
!
!  Communicate integer array between processors.
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: ibcast_array
      integer, optional :: proc,comm
!
      if (nbcast_array == 0) return
!
      call MPI_BCAST(ibcast_array,nbcast_array,MPI_INTEGER,ioptest(proc,root), &
                     ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_int_arr
!***********************************************************************
    subroutine mpibcast_int_arr2(ibcast_array,nbcast_array,proc,comm)
!
!  Communicate integer array(:,:) to other processors.
!
!  30-apr-17/Jorgen: adapted
!
      integer, dimension(2) :: nbcast_array
      integer, dimension(nbcast_array(1),nbcast_array(2)) :: ibcast_array
      integer, optional :: proc,comm
      integer :: num_elements
!
      if (any(nbcast_array == 0)) return
!
      num_elements = product(nbcast_array)
!
      call MPI_BCAST(ibcast_array,num_elements,MPI_INTEGER,ioptest(proc,root), &
                     ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_int_arr2
!***********************************************************************
    subroutine mpibcast_real_scl(bcast_array,proc,comm)
!
!  Communicate real scalar between processors.
!
      real :: bcast_array
      integer, optional :: proc, comm

      call MPI_BCAST(bcast_array,1,mpi_precision,ioptest(proc,root), &
                     ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_real_scl
!***********************************************************************
    subroutine mpibcast_real_arr(bcast_array,nbcast_array,proc,comm)
!
!  Communicate real array between processors.
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array
      integer, optional :: proc,comm
!
      if (nbcast_array == 0) return
!
      call MPI_BCAST(bcast_array,nbcast_array,mpi_precision,ioptest(proc,root), &
                     ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_real_arr
!***********************************************************************
    subroutine mpibcast_real_arr2(bcast_array,nbcast_array,proc,comm)
!
!  Communicate real array(:,:) to other processor.
!
!  25-feb-08/wlad: adapted
!
      integer, dimension(2) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2)) :: bcast_array
      integer, optional :: proc,comm
      integer :: ibcast_proc
!
      integer :: num_elements
!
      if (any(nbcast_array == 0)) return
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      num_elements = product(nbcast_array)
      call MPI_BCAST(bcast_array, num_elements, mpi_precision, ibcast_proc, &
                     ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_real_arr2
!***********************************************************************
    subroutine mpibcast_real_arr3(bcast_array,nbcast_array,proc)
!
!  Communicate real array(:,:,:) to other processor.
!
!  25-fev-08/wlad: adapted
!
      integer, dimension(3) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3)) :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      integer :: num_elements
!
      if (any(nbcast_array == 0)) return
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      num_elements = product(nbcast_array)
      call MPI_BCAST(bcast_array, num_elements, mpi_precision, ibcast_proc, &
                     MPI_COMM_GRID,mpierr)
!
    endsubroutine mpibcast_real_arr3
!***********************************************************************
    subroutine mpibcast_real_arr4(bcast_array,nbcast_array,proc)
!
!  Communicate real array(:,:,:,:) to other processor.
!
!  21-dec-10/ccyang: adapted
!
      integer, dimension(4) :: nbcast_array
      real, dimension(nbcast_array(1),nbcast_array(2),nbcast_array(3),nbcast_array(4)) :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc, num_elements
!
      if (any(nbcast_array == 0)) return
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      num_elements = product(nbcast_array)
      call MPI_BCAST(bcast_array, num_elements, mpi_precision, ibcast_proc, &
                     MPI_COMM_GRID,mpierr)
!
    endsubroutine mpibcast_real_arr4
!***********************************************************************
    subroutine mpibcast_double_scl(bcast_array,proc,comm)
!
!  Communicate real scalar between processors.
!
      real(KIND=rkind8) :: bcast_array
      integer, optional :: proc,comm
!
      call MPI_BCAST(bcast_array,1,MPI_DOUBLE_PRECISION,ioptest(proc,root), &
                     ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_double_scl
!***********************************************************************
    subroutine mpibcast_double_arr(bcast_array,nbcast_array,proc)
!
!  Communicate real array between processors.
!
      integer :: nbcast_array
      real(KIND=rkind8), dimension(nbcast_array) :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (nbcast_array == 0) return
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(bcast_array,nbcast_array,MPI_DOUBLE_PRECISION,ibcast_proc, &
                     MPI_COMM_GRID,mpierr)
!
    endsubroutine mpibcast_double_arr
!***********************************************************************
    subroutine mpibcast_char_scl(cbcast_array,proc,comm)
!
!  Communicate character scalar between processors.
!
      character(LEN=*) :: cbcast_array
      integer, optional :: proc,comm
!
      call MPI_BCAST(cbcast_array,len(cbcast_array),MPI_CHARACTER,ioptest(proc,root), &
                     ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_char_scl
!***********************************************************************
    subroutine mpibcast_char_arr(cbcast_array,nbcast_array,proc,comm)
!
!  Communicate character array between processors.
!
      integer :: nbcast_array
!
      character(LEN=*), dimension(nbcast_array) :: cbcast_array
      integer, optional :: proc,comm
!
      if (nbcast_array == 0) return
!
      call MPI_BCAST(cbcast_array,len(cbcast_array(1))*nbcast_array,MPI_CHARACTER, &
                     ioptest(proc,root),ioptest(comm,MPI_COMM_GRID),mpierr)
!
    endsubroutine mpibcast_char_arr
!***********************************************************************
    subroutine mpibcast_cmplx_arr(bcast_array,nbcast_array,proc)
!
!  Communicate real array between processors.
!
      integer :: nbcast_array
      complex, dimension(nbcast_array) :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (nbcast_array == 0) return
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(bcast_array,nbcast_array,MPI_CMPLX,ibcast_proc, &
                     MPI_COMM_GRID,mpierr)
!
    endsubroutine mpibcast_cmplx_arr
!***********************************************************************
    subroutine mpibcast_cmplx_arr_dbl(bcast_array,nbcast_array,proc)
!
!  Communicate real array between processors.
!
      integer :: nbcast_array
      complex(KIND=rkind8), dimension(nbcast_array) :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (nbcast_array == 0) return
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(bcast_array,nbcast_array,MPI_DOUBLE_COMPLEX,ibcast_proc, &
                     MPI_COMM_GRID,mpierr)
!
    endsubroutine mpibcast_cmplx_arr_dbl
!***********************************************************************
    subroutine mpibcast_cmplx_arr_sgl(bcast_array,nbcast_array,proc)
!
!  Communicate real array between processors.
!
      integer :: nbcast_array
      complex(KIND=rkind4), dimension(nbcast_array) :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc
!
      if (nbcast_array == 0) return
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(bcast_array,nbcast_array,MPI_COMPLEX,ibcast_proc, &
                     MPI_COMM_GRID,mpierr)
!
    endsubroutine mpibcast_cmplx_arr_sgl
!***********************************************************************
    subroutine mpiscatter_real_arr(src_array,dest_array,proc,comm)
!
!  Scatter real 1D-array between processors of communicator comm from rank proc.
!
      real, dimension(:) :: src_array, dest_array
      integer, optional :: proc,comm

      integer :: count,comm_,np

      comm_=ioptest(comm,MPI_COMM_PENCIL)

      call MPI_COMM_SIZE(comm_,np,mpierr)
      count=size(src_array)/np
      if (size(dest_array)<count) &
        call stop_fatal('mpiscatter_real_arr:'// &
                        'ERROR - sizes of source/destination arrays/number of processors not compatibel')

      call MPI_SCATTER(src_array,count,mpi_precision,dest_array,count,mpi_precision,ioptest(proc,root),comm_,mpierr)

    endsubroutine mpiscatter_real_arr
!***********************************************************************
    subroutine mpiscatter_real_arr2(src_array,dest_array,proc,comm)
!
!  Scatter real 2D-array between processors of communicator comm from rank proc.
!
      use General, only: ioptest, indgen

      real, dimension(:,:) :: src_array
      real, dimension(:,:) :: dest_array
      integer, optional :: proc,comm

      integer :: count1,count2,comm_,np,sizeofreal
      integer, dimension(:), allocatable :: sendcounts,displs
      integer :: block, segment
      integer :: src_sz1, src_sz2, dest_sz1, dest_sz2, locrank, i2

      comm_=ioptest(comm,MPI_COMM_PENCIL)

      dest_sz1=size(dest_array,1); dest_sz2=size(dest_array,2)
      call MPI_COMM_RANK(comm_,locrank,mpierr)

      if (locrank==0) then

        src_sz1=size(src_array,1)
        src_sz2=size(src_array,2)
        count1=src_sz1/dest_sz1
        count2=src_sz2/dest_sz2

        call MPI_COMM_SIZE(comm_,np,mpierr)
        if (np/=count1*count2) &
          call stop_fatal('mpiscatter_real_arr2:'// &
                          'ERROR - sizes of source/destination arrays/number of processors not compatibel')

        allocate(sendcounts(np),displs(np))
        sendcounts=1

        do i2=0,count2-1
          displs(1+i2*count1:(i2+1)*count1)=indgen(count1)-1+count1*dest_sz2*i2
        enddo

      else
        src_sz1=dest_sz1; src_sz2=dest_sz2
      endif

      call MPI_TYPE_CREATE_SUBARRAY(2, (/src_sz1,src_sz2/), (/dest_sz1,dest_sz2/), (/0,0/), &
                                    MPI_ORDER_FORTRAN, mpi_precision, block, mpierr)
      call MPI_TYPE_COMMIT(block,mpierr)
      call MPI_TYPE_SIZE(mpi_precision, sizeofreal, mpierr)
      call MPI_TYPE_CREATE_RESIZED(block, 0, dest_sz1*sizeofreal, segment, mpierr)
      call MPI_TYPE_COMMIT(segment,mpierr)

      call MPI_SCATTERV(src_array,sendcounts,displs,segment,dest_array,dest_sz1*dest_sz2, &
                        mpi_precision,ioptest(proc,root),comm_,mpierr)
      call MPI_TYPE_FREE(block,mpierr)
      call MPI_TYPE_FREE(segment,mpierr)

    endsubroutine mpiscatter_real_arr2
!***********************************************************************
    subroutine mpiallreduce_sum_scl(fsum_tmp,fsum,idir,comm)
!
!  Calculate total sum for each array element and return to all processors.
!
      real :: fsum_tmp,fsum
      integer, optional :: idir,comm
!
      integer :: mpiprocs
!
!  Sum over all processors and return to all (MPI_COMM_GRID).
!  Sum over x beams and return to the ipx=0 processors (MPI_COMM_XBEAM).
!  Sum over y beams and return to the ipy=0 processors (MPI_COMM_YBEAM).
!  Sum over z beams and return to the ipz=0 processors (MPI_COMM_ZBEAM).
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=ioptest(comm,MPI_COMM_GRID)
      endif
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, 1, mpi_precision, MPI_SUM, mpiprocs, mpierr)
!
    endsubroutine mpiallreduce_sum_scl
!***********************************************************************
    subroutine mpiallreduce_sum_arr(fsum_tmp,fsum,nreduce,idir,comm)
!
!  Calculate total sum for each array element and return to all processors.
!
!  3-oct-12/MR: communicator corrected
!
      integer :: nreduce
      real, dimension(nreduce) :: fsum_tmp,fsum
      integer, optional :: idir,comm
!
      integer :: mpiprocs
!
      if (nreduce==0) return
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=MPI_COMM_GRID
      endif
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, nreduce, mpi_precision, MPI_SUM, &
                         mpiprocs, mpierr)
!
    endsubroutine mpiallreduce_sum_arr
!***********************************************************************
    subroutine mpiallreduce_sum_arr2(fsum_tmp,fsum,nreduce,idir,comm)
!
!  Calculate total sum for each array element and return to all processors.
!
!  23-nov-08/wlad: included the idir possibility
!
      integer, dimension(2) :: nreduce
      real, dimension(nreduce(1),nreduce(2)) :: fsum_tmp,fsum
      integer, optional :: idir,comm
!
      integer :: mpiprocs, num_elements
!
      if (any(nreduce==0)) return
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=MPI_COMM_GRID
      endif
!
      num_elements = product(nreduce)
      call MPI_ALLREDUCE(fsum_tmp, fsum, num_elements, mpi_precision, MPI_SUM, &
                         mpiprocs, mpierr)
!
    endsubroutine mpiallreduce_sum_arr2
!***********************************************************************
    subroutine mpiallreduce_sum_arr3(fsum_tmp,fsum,nreduce,idir)
!
!  Calculate total sum for each array element and return to all processors.
!
!  23-nov-08/wlad: included the idir possibility
!
      integer, dimension(3) :: nreduce
      real, dimension(nreduce(1),nreduce(2),nreduce(3)) :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs, num_elements
!
      if (any(nreduce==0)) return
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=MPI_COMM_GRID
      endif
!
      num_elements = product(nreduce)
      call MPI_ALLREDUCE(fsum_tmp, fsum, num_elements, mpi_precision, MPI_SUM, &
                         mpiprocs, mpierr)
!
    endsubroutine mpiallreduce_sum_arr3
!***********************************************************************
    subroutine mpiallreduce_sum_arr4(fsum_tmp,fsum,nreduce,idir)
!
!  Calculate total sum for each array element and return to all processors.
!
!  13-sep-13/MR: derived from mpiallreduce_sum_arr3
!
      integer, dimension(4) :: nreduce
      real, dimension(nreduce(1),nreduce(2),nreduce(3),nreduce(4)) :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs, num_elements
!
      if (any(nreduce==0)) return
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=MPI_COMM_GRID
      endif
!
      num_elements = product(nreduce)
      call MPI_ALLREDUCE(fsum_tmp, fsum, num_elements, mpi_precision, MPI_SUM, &
                         mpiprocs, mpierr)
!
    endsubroutine mpiallreduce_sum_arr4
!***********************************************************************
    subroutine mpiallreduce_sum_arr5(fsum_tmp,fsum,nreduce,idir)
!
!  Calculate total sum for each array element and return to all processors.
!
!  23-apr-14/MR: derived from mpiallreduce_sum_arr4
!
      integer, dimension(5) :: nreduce
      real, dimension(nreduce(1),nreduce(2),nreduce(3),nreduce(4),nreduce(5)) :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs, num_elements
!
      if (any(nreduce==0)) return
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=MPI_COMM_GRID
      endif
!
      num_elements = product(nreduce)
      call MPI_ALLREDUCE(fsum_tmp, fsum, num_elements, mpi_precision, MPI_SUM, &
                         mpiprocs, mpierr)
!
    endsubroutine mpiallreduce_sum_arr5
!***********************************************************************
    subroutine mpiallreduce_sum_arr_inplace(fsum, n)
!
!  Calculate total sum for each array element and return to all
!  processors in place.
!
!  14-nov-20/ccyang: coded
!
      real, dimension(:), intent(inout) :: fsum
      integer, intent(in) :: n
!
      if (n <= 0) return
      call MPI_ALLREDUCE(MPI_IN_PLACE, fsum, n, mpi_precision, MPI_SUM, MPI_COMM_GRID, mpierr)
!
    endsubroutine mpiallreduce_sum_arr_inplace
!***********************************************************************
    subroutine mpiscan_int(num,offset,comm)
!
!  Calculate for each processor offset of local array within a global array.
!
      integer :: num,offset
      integer, optional :: comm
!
      call MPI_SCAN(num, offset, 1, MPI_INTEGER, MPI_SUM, &
                         ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpiscan_int
!***********************************************************************
    subroutine mpiallreduce_sum_int_scl(fsum_tmp,fsum,comm)
!
!  Calculate total sum for each array element and return to all processors.
!
      integer :: fsum_tmp,fsum
      integer, optional :: comm
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, 1, MPI_INTEGER, MPI_SUM, &
                         ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpiallreduce_sum_int_scl
!***********************************************************************
    subroutine mpiallreduce_sum_int_arr(fsum_tmp,fsum,nreduce,idir)
!
!  Calculate total sum for each array element and return to all processors.
!
!  28-Sep-16/MR: introduce optional argument idir.
!
      integer :: nreduce
      integer, dimension(nreduce) :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs
!
      if (nreduce==0) return
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=MPI_COMM_GRID
      endif
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, nreduce, MPI_INTEGER, MPI_SUM, &
                         mpiprocs, mpierr)
!
    endsubroutine mpiallreduce_sum_int_arr
!***********************************************************************
    subroutine mpiallreduce_sum_int_arr_inplace(fsum, n)
!
!  Calculate total sum for each array element and return to all
!  processors in place.
!
!  14-nov-20/ccyang: coded
!
      integer, dimension(:), intent(inout) :: fsum
      integer, intent(in) :: n
!
      if (n <= 0) return
      call MPI_ALLREDUCE(MPI_IN_PLACE, fsum, n, MPI_INTEGER, MPI_SUM, MPI_COMM_GRID, mpierr)
!
    endsubroutine mpiallreduce_sum_int_arr_inplace
!***********************************************************************
    subroutine mpiallreduce_max_scl_sgl(fmax_tmp,fmax,comm)
!
!  Calculate total maximum element and return to all processors.
!
      real(KIND=rkind4) :: fmax_tmp,fmax
      integer, optional :: comm
!
      call MPI_ALLREDUCE(fmax_tmp, fmax, 1, mpi_precision, MPI_MAX, &
                         ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpiallreduce_max_scl_sgl
!***********************************************************************
    subroutine mpiallreduce_max_scl_dbl(fmax_tmp,fmax,comm)
!
!  Calculate total maximum element and return to all processors.
!
      real(KIND=rkind8) :: fmax_tmp,fmax
      integer, optional :: comm
!
      call MPI_ALLREDUCE(fmax_tmp, fmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                         ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpiallreduce_max_scl_dbl
!***********************************************************************
    subroutine mpiallreduce_max_scl_int(imax_tmp,imax,comm)
!
!  Calculate total maximum and return to all processors.
!
      integer :: imax_tmp,imax
      integer, optional :: comm
!
      call MPI_ALLREDUCE(imax_tmp, imax, 1, MPI_INTEGER, MPI_MAX, &
                         ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpiallreduce_max_scl_int
!***********************************************************************
    subroutine mpiallreduce_min_scl_sgl(fmin_tmp,fmin,comm)
!
!  Calculate total minimum and return to all processors.
!
      real(KIND=rkind4) :: fmin_tmp,fmin
      integer, optional :: comm
!
      call MPI_ALLREDUCE(fmin_tmp, fmin, 1, mpi_precision, MPI_MIN, &
                         ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpiallreduce_min_scl_sgl
!***********************************************************************
    subroutine mpiallreduce_min_scl_dbl(fmin_tmp,fmin,comm)
!
!  Calculate total minimum and return to all processors.
!
      real(KIND=rkind8) :: fmin_tmp,fmin
      integer, optional :: comm
!
      call MPI_ALLREDUCE(fmin_tmp, fmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                         ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpiallreduce_min_scl_dbl
!***********************************************************************
    subroutine mpiallreduce_min_scl_int(imin_tmp,imin,comm)
!
!  Calculate total minimum and return to all processors.
!
      integer :: imin_tmp,imin
      integer, optional :: comm
!
      call MPI_ALLREDUCE(imin_tmp, imin, 1, MPI_INTEGER, MPI_MIN, &
                         ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpiallreduce_min_scl_int
!***********************************************************************
    subroutine mpiallreduce_max_arr(fmax_tmp,fmax,nreduce,comm)
!
!  Calculate total maximum for each array element and return to all processors.
!
      integer :: nreduce
      real, dimension(nreduce) :: fmax_tmp,fmax
      integer, optional :: comm
!
      if (nreduce==0) return
!
      call MPI_ALLREDUCE(fmax_tmp, fmax, nreduce, mpi_precision, MPI_MAX, &
                         ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpiallreduce_max_arr
!***********************************************************************
    subroutine mpiallreduce_and_scl(fland_tmp, fland, comm)
!
!  Calculate logical or over all procs and return to all processors.
!
!  14-feb-14/ccyang: coded
!
      logical, intent(in) :: fland_tmp
      logical, intent(out):: fland
      integer, intent(in), optional :: comm
!
      if (nprocs == 1) then
        fland = fland_tmp
      else
        call MPI_ALLREDUCE(fland_tmp, fland, 1, MPI_LOGICAL, MPI_LAND, ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine mpiallreduce_and_scl
!***********************************************************************
    subroutine mpiallreduce_or_scl(flor_tmp, flor, comm)
!
!  Calculate logical or over all procs and return to all processors.
!
!  14-feb-14/ccyang: coded
!
      logical, intent(in) :: flor_tmp
      logical, intent(out):: flor
      integer, intent(in), optional :: comm
!
      if (nprocs == 1) then
        flor = flor_tmp
      else
        call MPI_ALLREDUCE(flor_tmp, flor, 1, MPI_LOGICAL, MPI_LOR, ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine mpiallreduce_or_scl
!***********************************************************************
    subroutine mpiallreduce_or_arr_inplace(lor, n, comm)
!
!  Calculate logical or over all procs and return to all processors in
!  place.
!
!  14-nov-20/ccyang: coded
!
      logical, dimension(:), intent(inout) :: lor
      integer, intent(in) :: n
      integer, intent(in), optional :: comm
!
      if (n <= 0) return
      call MPI_ALLREDUCE(MPI_IN_PLACE, lor, n, MPI_LOGICAL, MPI_LOR, ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpiallreduce_or_arr_inplace
!***********************************************************************
    subroutine mpireduce_max_scl(fmax_tmp,fmax,comm)
!
!  Calculate total maximum for each array element and return to root.
!
      real :: fmax_tmp,fmax
      integer, intent(in), optional :: comm
!
      call MPI_REDUCE(fmax_tmp, fmax, 1, mpi_precision, MPI_MAX, root, &
                      ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpireduce_max_scl
!***********************************************************************
    subroutine mpireduce_max_scl_int(fmax_tmp,fmax,comm)
!
!  Calculate total maximum for each array element and return to root.
!
      integer :: fmax_tmp,fmax
      integer, optional :: comm
!
      call MPI_REDUCE(fmax_tmp, fmax, 1, MPI_INTEGER, MPI_MAX, root, &
                      ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpireduce_max_scl_int
!***********************************************************************
    subroutine mpireduce_max_arr_int(fmax_tmp,fmax,nreduce,comm)
!
!  Calculate total maximum for each array element and return to root.
!
      integer :: nreduce
      integer, dimension(nreduce) :: fmax_tmp,fmax
      integer, optional :: comm
!
      if (nreduce==0) return
!
      call MPI_REDUCE(fmax_tmp, fmax, nreduce, MPI_INTEGER, MPI_MAX, root, &
                      ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpireduce_max_arr_int
!***********************************************************************
    subroutine mpireduce_max_arr(fmax_tmp,fmax,nreduce,comm)
!
!  Calculate total maximum for each array element and return to root.
!
      integer :: nreduce
      real, dimension(nreduce) :: fmax_tmp,fmax
      integer, optional :: comm
!
      if (nreduce==0) return
!
      call MPI_REDUCE(fmax_tmp, fmax, nreduce, mpi_precision, MPI_MAX, root, &
                      ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpireduce_max_arr
!***********************************************************************
    subroutine mpireduce_min_scl(fmin_tmp,fmin,comm)
!
!  Calculate total minimum for each array element and return to root.
!
      real :: fmin_tmp,fmin
      integer, optional :: comm
!
      call MPI_REDUCE(fmin_tmp, fmin, 1, mpi_precision, MPI_MIN, root, &
                      ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpireduce_min_scl
!***********************************************************************
    subroutine mpireduce_min_arr(fmin_tmp,fmin,nreduce,comm)
!
!  Calculate total maximum for each array element and return to root.
!
      integer :: nreduce
      real, dimension(nreduce) :: fmin_tmp,fmin
      integer, optional :: comm
!
      if (nreduce==0) return
!
      call MPI_REDUCE(fmin_tmp, fmin, nreduce, mpi_precision, MPI_MIN, root, &
                      ioptest(comm,MPI_COMM_GRID), mpierr)
!
    endsubroutine mpireduce_min_arr
!***********************************************************************
    subroutine mpireduce_sum_int_scl(fsum_tmp,fsum,comm)
!
!  Calculate sum and return to root.
!
      integer :: fsum_tmp,fsum
      integer, optional :: comm
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, 1, MPI_INTEGER, MPI_SUM, root, &
                        ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine mpireduce_sum_int_scl
!***********************************************************************
    subroutine mpireduce_sum_int_arr(fsum_tmp,fsum,nreduce,comm)
!
!  Calculate total sum for each array element and return to root.
!
      integer :: nreduce
      integer, dimension(nreduce) :: fsum_tmp,fsum
      integer, optional :: comm
!
      if (nreduce==0) return
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, nreduce, MPI_INTEGER, MPI_SUM, root, &
                        ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine mpireduce_sum_int_arr
!***********************************************************************
    subroutine mpireduce_sum_int_arr2(fsum_tmp,fsum,nreduce,comm)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(2) :: nreduce
      integer, dimension(nreduce(1),nreduce(2)) :: fsum_tmp,fsum
      integer, optional :: comm
!
      if (any(nreduce==0)) return
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, product(nreduce), MPI_INTEGER, MPI_SUM, root, &
                        ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine mpireduce_sum_int_arr2
!***********************************************************************
    subroutine mpireduce_sum_int_arr3(fsum_tmp,fsum,nreduce,comm)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(3) :: nreduce
      integer, dimension(nreduce(1),nreduce(2),nreduce(3)) :: fsum_tmp,fsum
      integer, optional :: comm
!
      if (any(nreduce==0)) return
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, product(nreduce), MPI_INTEGER, MPI_SUM, root, &
                        ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine mpireduce_sum_int_arr3
!***********************************************************************
    subroutine mpireduce_sum_int_arr4(fsum_tmp,fsum,nreduce,comm)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(4) :: nreduce
      integer, dimension(nreduce(1),nreduce(2),nreduce(3),nreduce(4)) :: fsum_tmp,fsum
      integer, optional :: comm
!
      if (any(nreduce==0)) return
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, product(nreduce), MPI_INTEGER, MPI_SUM, root, &
                        ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine mpireduce_sum_int_arr4
!***********************************************************************
    subroutine mpireduce_sum_scl(fsum_tmp,fsum,idir)
!
!  Calculate total sum and return to root.
!
      real :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
!
!  Sum over all processors and return to root (MPI_COMM_GRID).
!  Sum over x beams and return to the ipx=0 processors (MPI_COMM_XBEAM).
!  Sum over y beams and return to the ipy=0 processors (MPI_COMM_YBEAM).
!  Sum over z beams and return to the ipz=0 processors (MPI_COMM_ZBEAM).
!
        if (present(idir)) then
          mpiprocs=mpigetcomm(idir)
        else
          mpiprocs=MPI_COMM_GRID
        endif
        call MPI_REDUCE(fsum_tmp, fsum, 1, mpi_precision, MPI_SUM, root, &
                        mpiprocs, mpierr)
      endif
!
    endsubroutine mpireduce_sum_scl
!***********************************************************************
    subroutine mpireduce_sum_arr(fsum_tmp,fsum,nreduce,idir,comm)
!
!  Calculate total sum for each array element and return to root.
!
      integer :: nreduce
      real, dimension(nreduce) :: fsum_tmp,fsum
      integer, optional :: idir,comm
!
      integer :: mpiprocs
!
      intent(in)  :: fsum_tmp,nreduce
      intent(out) :: fsum
!
      if (nreduce==0) return
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        if (present(idir)) then
          mpiprocs=mpigetcomm(idir)
        else
          mpiprocs=ioptest(comm,MPI_COMM_GRID)
        endif
        call MPI_REDUCE(fsum_tmp, fsum, nreduce, mpi_precision, MPI_SUM, root, &
                        mpiprocs, mpierr)
      endif
!
    endsubroutine mpireduce_sum_arr
!***********************************************************************
    subroutine mpireduce_sum_arr2(fsum_tmp,fsum,nreduce,idir,inplace)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(2) :: nreduce
      real, dimension(nreduce(1),nreduce(2)) :: fsum_tmp, fsum
      integer, optional :: idir
      logical, optional :: inplace
!
      integer :: mpiprocs, num_elements
      logical :: inplace_opt
      logical :: MPI_IN_PLACE
!
      intent(in)  :: fsum_tmp,nreduce
      intent(out) :: fsum
!
      if (any(nreduce==0)) return
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        if (present(idir)) then
          mpiprocs=mpigetcomm(idir)
        else
          mpiprocs=MPI_COMM_GRID
        endif
        if (present(inplace)) then
          inplace_opt=inplace
        else
          inplace_opt=.false.
        endif
        num_elements = product(nreduce)
        if (inplace_opt) then
          call MPI_REDUCE(MPI_IN_PLACE, fsum, num_elements, mpi_precision, &
                          MPI_SUM, root, mpiprocs, mpierr)
        else
          call MPI_REDUCE(fsum_tmp, fsum, num_elements, mpi_precision, MPI_SUM, &
                          root, mpiprocs, mpierr)
        endif
      endif
!
    endsubroutine mpireduce_sum_arr2
!***********************************************************************
    subroutine mpireduce_sum_arr3(fsum_tmp,fsum,nreduce,idir)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(3) :: nreduce
      real, dimension(nreduce(1),nreduce(2),nreduce(3)) :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs, num_elements
!
      intent(in)  :: fsum_tmp,nreduce
      intent(out) :: fsum
!
      if (any(nreduce==0)) return
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        if (present(idir)) then
          mpiprocs=mpigetcomm(idir)
        else
          mpiprocs=MPI_COMM_GRID
        endif
        num_elements = product(nreduce)
        call MPI_REDUCE(fsum_tmp, fsum, num_elements, mpi_precision, MPI_SUM, &
                        root, mpiprocs, mpierr)
      endif
!
    endsubroutine mpireduce_sum_arr3
!***********************************************************************
    subroutine mpireduce_sum_arr4(fsum_tmp,fsum,nreduce,idir)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(4) :: nreduce
      real, dimension(nreduce(1),nreduce(2),nreduce(3),nreduce(4)) :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs, num_elements
!
      intent(in)  :: fsum_tmp,nreduce
      intent(out) :: fsum
!
      if (any(nreduce==0)) return
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        if (present(idir)) then
          mpiprocs=mpigetcomm(idir)
        else
          mpiprocs=MPI_COMM_GRID
        endif
        num_elements = product(nreduce)
        call MPI_REDUCE(fsum_tmp, fsum, num_elements, mpi_precision, MPI_SUM, &
                        root, mpiprocs, mpierr)
      endif
!
    endsubroutine mpireduce_sum_arr4
!***********************************************************************
    subroutine mpireduce_or_scl(flor_tmp,flor,comm)
!
!  Calculate logical or over all procs and return to root.
!
!  17-sep-05/anders: coded
!
      logical :: flor_tmp, flor
      integer, optional :: comm
!
      if (nprocs==1) then
        flor=flor_tmp
      else
        call MPI_REDUCE(flor_tmp, flor, 1, MPI_LOGICAL, MPI_LOR, root, &
                        ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine mpireduce_or_scl
!***********************************************************************
    subroutine mpireduce_or_arr(flor_tmp,flor,nreduce,comm)
!
!  Calculate logical or over all procs and return to root.
!
!  17-sep-05/anders: coded
!
      integer :: nreduce
      logical, dimension(nreduce) :: flor_tmp, flor
      integer, optional :: comm
!
      if (nreduce==0) return
!
      if (nprocs==1) then
        flor=flor_tmp
      else
        call MPI_REDUCE(flor_tmp, flor, nreduce, MPI_LOGICAL, MPI_LOR, root, &
                        ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine mpireduce_or_arr
!***********************************************************************
    subroutine mpireduce_and_scl(fland_tmp,fland,comm)
!
!  Calculate logical and over all procs and return to root.
!
!  17-sep-05/anders: coded
!
      logical :: fland_tmp, fland
      integer, optional :: comm
!
      if (nprocs==1) then
        fland=fland_tmp
      else
        call MPI_REDUCE(fland_tmp, fland, 1, MPI_LOGICAL, MPI_LAND, root, &
                        ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine mpireduce_and_scl
!***********************************************************************
    subroutine mpireduce_and_arr(fland_tmp,fland,nreduce,comm)
!
!  Calculate logical and over all procs and return to root.
!
!  11-mar-09/anders: coded
!
      integer :: nreduce
      logical, dimension(nreduce) :: fland_tmp, fland
      integer, optional :: comm
!
      if (nreduce==0) return
!
      if (nprocs==1) then
        fland=fland_tmp
      else
        call MPI_REDUCE(fland_tmp, fland, nreduce, MPI_LOGICAL, MPI_LAND, root,&
                        ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine mpireduce_and_arr
!***********************************************************************
    subroutine start_serialize
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
      serial_level=serial_level+1
      if (serial_level>1) return
!
      buf = 0
      if (.not. lroot) then     ! root starts, others wait for permission
        call MPI_RECV(buf,1,MPI_INTEGER,root,io_perm,MPI_COMM_GRID,status,mpierr)
      endif
!
    endsubroutine start_serialize
!***********************************************************************
    subroutine end_serialize
!
!  Do block between start_serialize and end_serialize serially in iproc order.
!
!  19-nov-02/wolf: coded
!
      integer :: i,buf
      integer, dimension(MPI_STATUS_SIZE) :: status
!
      serial_level=serial_level-1
      if (serial_level>=1) return
      if (serial_level<0) &
          call stop_it('end_serialize: too many end_serialize calls')
!
      buf = 0
      if (lroot) then
        do i=1,ncpus-1            ! send permission, wait for success message
          call MPI_SEND(buf,1,MPI_INTEGER,i,io_perm,MPI_COMM_GRID,mpierr)
          call MPI_RECV(buf,1,MPI_INTEGER,i,io_succ,MPI_COMM_GRID,status,mpierr)
        enddo
      else                  ! tell root we're done
        call MPI_SEND(buf,1,MPI_INTEGER,root,io_succ,MPI_COMM_GRID,mpierr)
      endif
!
    endsubroutine end_serialize
!***********************************************************************
    subroutine mpibarrier(comm)
!
!  Synchronize nodes.
!
!  23-jul-2002/wolf: coded
!
      integer, optional, intent(IN) :: comm
!
      call MPI_BARRIER(ioptest(comm,MPI_COMM_PENCIL), mpierr)
!
    endsubroutine mpibarrier
!***********************************************************************
    subroutine mpifinalize
!
!  Send stop signal to foreign code.
!
      if (lforeign.and.lroot) &
        call MPI_SEND(.true.,1,MPI_LOGICAL,frgn_setup%root,tag_foreign,MPI_COMM_WORLD,mpierr)

      call MPI_TYPE_FREE(REAL_ARR_MAXSIZE, mpierr)

      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
      call MPI_FINALIZE(mpierr)
!
      if (lforeign) then
        if (allocated(frgn_setup%recv_req)) deallocate(frgn_setup%recv_req)
        if (allocated(frgn_setup%xgrid)) deallocate(frgn_setup%xgrid)
        if (allocated(frgn_setup%ygrid)) deallocate(frgn_setup%ygrid)
        if (allocated(frgn_setup%yind_rng)) deallocate(frgn_setup%yind_rng)
      endif
!
    endsubroutine mpifinalize
!***********************************************************************
    subroutine mpiabort
!
!  Aborts the processes in PC (it is documented that in most systems this aborts all MPI processes)
!
    call MPI_ABORT(MPI_COMM_PENCIL,-1,mpierr)
       
    endsubroutine mpiabort
!***********************************************************************
    function mpiwtime()
!
      real(KIND=rkind8) :: mpiwtime
      real(KIND=rkind8) :: MPI_WTIME   ! definition needed for mpicomm_ to work
!
      mpiwtime = MPI_WTIME()
!
    endfunction mpiwtime
!***********************************************************************
    function mpiwtick()
!
      real(KIND=rkind8) :: mpiwtick
      real(KIND=rkind8) :: MPI_WTICK   ! definition needed for mpicomm_ to work
!
      mpiwtick = MPI_WTICK()
!
    endfunction mpiwtick
!***********************************************************************
    subroutine die_gracefully
!
!  Stop having shutdown MPI neatly.
!  With at least some MPI implementations, this only stops if all
!  processors agree to call die_gracefully.
!
!  29-jun-05/tony: coded
!
!  Tell the world something went wrong -- mpirun may not propagate
!  an error status.
!
      use General, only: touch_file

      call touch_file('ERROR')
!
!  mpifinalize may fail when called inside a parallel section need to be looked at further
!
      call mpifinalize
      if (lroot) then
        STOP 1                    ! Return nonzero exit status
      else
        STOP
      endif
!
    endsubroutine die_gracefully
!***********************************************************************
    subroutine die_immediately
!
!  Stop without shutting down MPI.
!  For those MPI implementations, which only finalize when all
!  processors agree to finalize.
!
!  29-jun-05/tony: coded
!
!  Tell the world something went wrong -- mpirun may not propagate
!  an error status.
!
      use General, only: touch_file

      call touch_file('ERROR')
!
      if (lroot) then
        STOP 2                    ! Return nonzero exit status
      else
        STOP
      endif
!
    endsubroutine die_immediately
!***********************************************************************
    subroutine stop_fatal(msg,force)
!
!  Print message and stop. If force, stop without shutting down MPI.
!
!  13-dez-10/Bourdin.KIS: coded
!
      character (len=*) :: msg
      logical, optional :: force
!
      logical :: immediately
!
      immediately = .false.
      if (present (force)) immediately = force
!
      if (lroot .or. immediately) write(*,'(A,A)') 'STOPPED FATAL: ', msg
!
      if (immediately) then
        call die_immediately
      else
        call die_gracefully
      endif
!
    endsubroutine stop_fatal
!***********************************************************************
    subroutine stop_it(msg,code)
!
!  Print message and stop.
!  With at least some MPI implementations, this only stops if all
!  processors agree to call stop_it. To stop (collectively) if only one
!  or a few processors find some condition, use stop_it_if_any.
!
!  6-nov-01/wolf: coded
!  4-nov-11/MR: optional parameter code added
!
      use General, only: itoa
!
      character (len=*) :: msg
      integer, optional :: code
!
      if (lroot) then
        if (present(code)) then
          write(*,'(A,A)') 'STOPPED: ', msg, '. CODE: '//trim(itoa(code))
        elseif (msg /= '') then
          write(*,'(A,A)') 'STOPPED: ', msg
        else
          write(*,'(A)') 'STOPPED'
        endif
      endif
!
      call die_gracefully
!
    endsubroutine stop_it
!***********************************************************************
    subroutine stop_it_if_any(stop_flag,msg,code)
!
!  Conditionally print message and stop.
!  This works unilaterally, i.e. if STOP_FLAG is true on _any_ processor,
!  we will all stop. The error message will be printed together with
!  the MPI rank number, if the message is not empty.
!
!  22-nov-04/wolf: coded
!  15-Feb-2012/Bourdin.KIS: optional parameter code added
!
      logical :: stop_flag
      character (len=*) :: msg
      integer, optional :: code
!
      logical :: global_stop_flag, identical_stop_flag
!
!  Get global OR of stop_flag and distribute it, so all processors agree
!  on whether to call stop_it:
!
      call MPI_ALLREDUCE(stop_flag,global_stop_flag,1,MPI_LOGICAL, &
                         MPI_LOR,MPI_COMM_PENCIL,mpierr)
      call MPI_ALLREDUCE(stop_flag,identical_stop_flag,1,MPI_LOGICAL, &
                         MPI_LAND,MPI_COMM_PENCIL,mpierr)
!
      if (global_stop_flag) then
        if ((.not. lroot) .and. (.not. identical_stop_flag) .and. (msg/='')) &
            write(*,'(A,I8,A,A)') 'RANK ', iproc_world, ' STOPPED: ', msg
        call stop_it(msg, code)
      endif
!
    endsubroutine stop_it_if_any
!***********************************************************************
    subroutine check_emergency_brake
!
!  Check the lemergency_brake flag and stop with any provided
!  message if it is set.
!
!  29-jul-06/tony: coded
!
      logical :: global_stop_flag
!
!  Get global OR of lemergency_brake and distribute it, so all
!  processors agree on whether to call stop_it:
!
      call MPI_ALLREDUCE(lemergency_brake,global_stop_flag,1,MPI_LOGICAL, &
                         MPI_LOR,MPI_COMM_PENCIL,mpierr)
!
      if (global_stop_flag) call stop_it( &
            "Emergency brake activated. Check for error messages above.")
!
    endsubroutine check_emergency_brake
!***********************************************************************
    subroutine transp(a,var, comm, lsync)
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
      integer :: ystag=111,yrtag=112,zstag=113,zrtag=114,partner
      integer :: m,n,ibox,ix
      integer, optional :: comm
      logical, optional :: lsync
!
!  Doing x-y transpose if var='y'
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      if (var=='y') then
!
        if (mod(nxgrid,nygrid)/=0) then
          print*,'transp: nxgrid needs to be an integer multiple of '//&
                 'nygrid for var==y'
          call stop_it_if_any(.true.,'Inconsistency: mod(nxgrid,nygrid)/=0')
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
                call MPI_SEND(send_buf_y,sendc_y,mpi_precision,partner,ystag,MPI_COMM_GRID,mpierr)
                call MPI_RECV(recv_buf_y,recvc_y,mpi_precision,partner,yrtag,MPI_COMM_GRID,stat,mpierr)
              elseif (px>ipy) then  ! below diagonal: receive first, send then
                call MPI_RECV(recv_buf_y,recvc_y,mpi_precision,partner,ystag,MPI_COMM_GRID,stat,mpierr)
                call MPI_SEND(send_buf_y,sendc_y,mpi_precision,partner,yrtag,MPI_COMM_GRID,mpierr)
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
!
        if (nzgrid/=nxgrid) then
          if (lroot) print*, 'transp: need to have nzgrid=nxgrid for var==z'
          call stop_it_if_any(.true.,'transp: inconsistency - nzgrid/=nxgrid')
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
              call MPI_SEND (send_buf_z,sendc_z,mpi_precision,partner,zstag,MPI_COMM_GRID,mpierr)
              call MPI_RECV (recv_buf_z,recvc_z,mpi_precision,partner,zrtag,MPI_COMM_GRID,stat,mpierr)
            elseif (px>ipz) then  ! below diagonal: receive first, send then
              call MPI_RECV (recv_buf_z,recvc_z,mpi_precision,partner,zstag,MPI_COMM_GRID,stat,mpierr)
              call MPI_SEND(send_buf_z,sendc_z,mpi_precision,partner,zrtag,MPI_COMM_GRID,mpierr)
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
      !$omp end single
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
    endsubroutine transp
!***********************************************************************
    subroutine transp_xy(a,comm,lsync)
!
!  Doing the transpose of information distributed on several processors.
!  This routine transposes 2D arrays in x and y only.
!
!   6-oct-06/tobi: Adapted from transp
!
! TODO: Implement nygrid = n*nxgrid
!
      real, dimension(nx,ny), intent(inout) :: a
!
      real, dimension(ny,ny) :: send_buf_y, recv_buf_y, tmp
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc_y,recvc_y,px
      integer :: ytag=101,partner
      integer :: ibox,iy
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      if (nprocx>1) then
        print*,'transp_xy: nprocx must be equal to 1'
        call stop_it_if_any(.true.,'Inconsistency: nprocx>1')
      endif
!
      if (mod(nxgrid,nygrid)/=0) then
        print*,'transp_xy: nxgrid needs to be an integer multiple of nygrid'
        call stop_it_if_any(.true.,'Inconsistency: mod(nxgrid,nygrid)/=0')
      endif
!
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.
!
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
!  avoids deadlocks, because there will always be at least one matching
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
              call MPI_SEND(send_buf_y,sendc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),mpierr)
              call MPI_RECV(recv_buf_y,recvc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),stat,mpierr)
            elseif (px>ipy) then  ! below diagonal: receive first, send then
              call MPI_RECV(recv_buf_y,recvc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),stat,mpierr)
              call MPI_SEND(send_buf_y,sendc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),mpierr)
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
!    original             2x2 blocks           each block
!                         transposed           transposed
!
      do px=0,nprocy-1
        do ibox=0,nxgrid/nygrid-1
          iy=(ibox*nprocy+px)*ny
          tmp=transpose(a(iy+1:iy+ny,:)); a(iy+1:iy+ny,:)=tmp
        enddo
      enddo
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine transp_xy
!***********************************************************************
    subroutine transp_xy_other(a,comm,lsync)
!
!  Doing the transpose of information distributed on several processors.
!  This routine transposes 2D arrays of arbitrary size in x and y only.
!
!   6-oct-06/tobi: Adapted from transp
!
! TODO: Implement nygrid = n*nxgrid
!
      real, dimension(:,:), intent(inout) :: a
!
      real, dimension(size(a,2),size(a,2)) :: send_buf_y, recv_buf_y, tmp
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc_y,recvc_y,px
      integer :: ytag=101,partner
      integer :: ibox,iy,nx_other,ny_other
      integer :: nxgrid_other,nygrid_other
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      nx_other=size(a,1); ny_other=size(a,2)
      nxgrid_other=nx_other
      nygrid_other=ny_other*nprocy
!
      if (mod(nxgrid_other,nygrid_other)/=0) then
        print*,'transp: nxgrid_other needs to be an integer multiple of '//&
               'nygrid_other for var==y'
        call stop_it_if_any(.true.,'Inconsistency: mod(nxgrid_other,nygrid_other)/=0')
      endif
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.
!
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
              call MPI_SEND(send_buf_y,sendc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),mpierr)
              call MPI_RECV(recv_buf_y,recvc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),stat,mpierr)
            elseif (px>ipy) then  ! below diagonal: receive first, send then
              call MPI_RECV(recv_buf_y,recvc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),stat,mpierr)
              call MPI_SEND(send_buf_y,sendc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),mpierr)
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
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine transp_xy_other
!***********************************************************************
    subroutine transp_other(a,var,comm,lsync)
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
      integer :: ytag=101,ztag=202,partner
      integer :: m,n,ibox,ix,nx_other,ny_other,nz_other
      integer :: nxgrid_other,nygrid_other,nzgrid_other
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
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
          call stop_it_if_any(.true.,'Inconsistency: mod(nxgrid_other,nygrid_other)/=0')
        endif
!
!  Allocate temporary scratch array
!
        allocate (tmp(ny_other,ny_other))
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.
!
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
                call MPI_SEND(send_buf_y,sendc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),mpierr)
                call MPI_RECV(recv_buf_y,recvc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),stat,mpierr)
              elseif (px>ipy) then  ! below diagonal: receive first, send then
                call MPI_RECV(recv_buf_y,recvc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),stat,mpierr)
                call MPI_SEND(send_buf_y,sendc_y,mpi_precision,partner,ytag,ioptest(comm,MPI_COMM_GRID),mpierr)
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
!
        if (nzgrid_other/=nxgrid_other) then
          if (lroot) print*, 'transp_other: need to have '//&
          'nzgrid_other=nxgrid_other for var==z'
          call stop_it_if_any(.true.,'transp_other: inconsistency - nzgrid/=nxgrid')
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
              call MPI_SEND(send_buf_z,sendc_z,mpi_precision,partner,ztag,ioptest(comm,MPI_COMM_GRID),mpierr)
              call MPI_RECV (recv_buf_z,recvc_z,mpi_precision,partner,ztag,ioptest(comm,MPI_COMM_GRID),stat,mpierr)
            elseif (px>ipz) then  ! below diagonal: receive first, send then
              call MPI_RECV (recv_buf_z,recvc_z,mpi_precision,partner,ztag,ioptest(comm,MPI_COMM_GRID),stat,mpierr)
              call MPI_SEND(send_buf_z,sendc_z,mpi_precision,partner,ztag,ioptest(comm,MPI_COMM_GRID),mpierr)
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
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine transp_other
!***********************************************************************
    subroutine transp_xz(a,b,comm,lsync)
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
      real, dimension(nxt,nz) :: buf
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc,px
      integer :: ztag=101,partner
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      if (mod(nxgrid,nprocz)/=0) then
        print*,'transp_xz: nxgrid needs to be an integer multiple of nprocz'
        call stop_it_if_any(.true.,'Inconsistency: mod(nxgrid,nprocz)/=0')
      endif
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.
!
      sendc=nxt*nz
!
!  Send information to different processors (x-z transpose)
!
      b(ipz*nz+1:(ipz+1)*nz,:)=transpose(a(ipz*nxt+1:(ipz+1)*nxt,:))
      do px=0,nprocz-1
        if (px/=ipz) then
          partner=ipy+px*nprocy ! = iproc + (px-ipz)*nprocy
          buf=a(px*nxt+1:(px+1)*nxt,:)
          call MPI_SENDRECV_REPLACE(buf,sendc,mpi_precision,partner,ztag,partner,ztag,ioptest(comm,MPI_COMM_GRID),stat,mpierr)
          b(px*nz+1:(px+1)*nz,:)=transpose(buf)
        endif
      enddo
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine transp_xz
!***********************************************************************
    subroutine transp_zx(a,b,comm,lsync)
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
      real, dimension(nz,nxt) :: buf
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc,px
      integer :: ztag=101,partner
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      if (mod(nxgrid,nprocz)/=0) then
        print*,'transp_zx: nxgrid needs to be an integer multiple of nprocz'
        call stop_it_if_any(.true.,'Inconsistency: mod(nxgrid,nprocz)/=0')
      endif
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.
!
      sendc=nz*nxt
!
!  Send information to different processors (x-z transpose)
!
      b(ipz*nxt+1:(ipz+1)*nxt,:)=transpose(a(ipz*nz+1:(ipz+1)*nz,:))
      do px=0,nprocz-1
        if (px/=ipz) then
          partner=ipy+px*nprocy ! = iproc + (px-ipz)*nprocy
          buf=a(px*nz+1:(px+1)*nz,:)
          call MPI_SENDRECV_REPLACE(buf,sendc,mpi_precision,partner,ztag,partner,ztag,ioptest(comm,MPI_COMM_GRID),stat,mpierr)
          b(px*nxt+1:(px+1)*nxt,:)=transpose(buf)
        endif
      enddo
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine transp_zx
!***********************************************************************
    subroutine communicate_vect_field_ghosts(f,topbot,start_index,comm)
!
!  Helper routine for communication of ghost cell values of a vector field.
!  Needed by potential field extrapolations, which only compute nx*ny arrays.
!  Can also be used for synchronization of changed uu values with ghost cells,
!  if the start_index parameter set to iux (default is iax).
!
!   8-oct-2006/tobi: Coded
!  28-dec-2010/Bourdin.KIS: extended to work for any 3D vector field data.
!
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent (in) :: topbot
      integer, intent(in), optional :: start_index
!
      real, dimension (nx,nghost,nghost+1,3) :: lbufyo,ubufyo,lbufyi,ubufyi
      real, dimension (nghost,size(f,2),nghost+1,3) :: lbufxo,ubufxo,lbufxi,ubufxi
      integer :: nbufx,nbufy,nn1,nn2,is,ie
      integer, optional :: comm
!
      is = iax
      if (present (start_index)) is = start_index
      ie = is + 2
!
      nn1=-1; nn2=-1
      select case (topbot)
        case (BOT); nn1=1;  nn2=n1
        case (TOP); nn1=n2; nn2=size(f,3)
        case default; call stop_it_if_any(.true.,"communicate_vect_field_ghosts: topbot should be BOT or TOP")
      end select
!
!  Periodic boundaries in y -- communicate along y if necessary
!
      if (nprocy>1) then
!
        lbufyo = f(l1:l2, m1:m1i,nn1:nn2,is:ie)
        ubufyo = f(l1:l2,m2i:m2 ,nn1:nn2,is:ie)
!
        nbufy=nx*nghost*(nghost+1)*3
!
        call MPI_IRECV(ubufyi,nbufy,mpi_precision,yuneigh,tolowy, &
                       ioptest(comm,MPI_COMM_GRID),irecv_rq_fromuppy,mpierr)
        call MPI_IRECV(lbufyi,nbufy,mpi_precision,ylneigh,touppy, &
                       ioptest(comm,MPI_COMM_GRID),irecv_rq_fromlowy,mpierr)
        call MPI_ISEND(lbufyo,nbufy,mpi_precision,ylneigh,tolowy, &
                       ioptest(comm,MPI_COMM_GRID),isend_rq_tolowy,mpierr)
        call MPI_ISEND(ubufyo,nbufy,mpi_precision,yuneigh,touppy, &
                       ioptest(comm,MPI_COMM_GRID),isend_rq_touppy,mpierr)
!
        call MPI_WAIT(irecv_rq_fromuppy,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowy,irecv_stat_fl,mpierr)
!
        f(l1:l2,   1:m1-1,nn1:nn2,is:ie) = lbufyi
        f(l1:l2,m2+1:    ,nn1:nn2,is:ie) = ubufyi
!
        call MPI_WAIT(isend_rq_tolowy,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppy,isend_stat_tu,mpierr)
!
      else
!
        f(l1:l2,   1:m1-1,nn1:nn2,is:ie) = f(l1:l2,m2i:m2 ,nn1:nn2,is:ie)
        f(l1:l2,m2+1:    ,nn1:nn2,is:ie) = f(l1:l2, m1:m1i,nn1:nn2,is:ie)
!
      endif
!
!  Periodic boundaries in x
!
      if (nprocx>1) then
!
        lbufxo = f( l1:l1i,:,nn1:nn2,is:ie)
        ubufxo = f(l2i:l2 ,:,nn1:nn2,is:ie)
!
        nbufx=nghost*size(f,2)*(nghost+1)*3
!
        call MPI_IRECV(ubufxi,nbufx,mpi_precision,xuneigh,tolowx, &
                       ioptest(comm,MPI_COMM_GRID),irecv_rq_fromuppx,mpierr)
        call MPI_IRECV(lbufxi,nbufx,mpi_precision,xlneigh,touppx, &
                       ioptest(comm,MPI_COMM_GRID),irecv_rq_fromlowx,mpierr)
        call MPI_ISEND(lbufxo,nbufx,mpi_precision,xlneigh,tolowx, &
                       ioptest(comm,MPI_COMM_GRID),isend_rq_tolowx,mpierr)
        call MPI_ISEND(ubufxo,nbufx,mpi_precision,xuneigh,touppx, &
                       ioptest(comm,MPI_COMM_GRID),isend_rq_touppx,mpierr)
!
        call MPI_WAIT(irecv_rq_fromuppx,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowx,irecv_stat_fl,mpierr)
!
        f(   1:l1-1,:,nn1:nn2,is:ie) = lbufxi
        f(l2+1:    ,:,nn1:nn2,is:ie) = ubufxi
!
        call MPI_WAIT(isend_rq_tolowx,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppx,isend_stat_tu,mpierr)
!
      else
!
        f(   1:l1-1,:,nn1:nn2,is:ie) = f(l2i:l2 ,:,nn1:nn2,is:ie)
        f(l2+1:    ,:,nn1:nn2,is:ie) = f( l1:l1i,:,nn1:nn2,is:ie)
!
      endif
!
    endsubroutine communicate_vect_field_ghosts
!***********************************************************************
    subroutine communicate_xy_ghosts(data,comm)
!
!  Helper routine for communication of ghost cells in horizontal direction.
!
!  11-apr-2011/Bourdin.KIS: adapted from communicate_vect_field_ghosts.
!
      real, dimension (:,:), intent (inout) :: data
!
      real, dimension (nx,nghost) :: lbufyo,ubufyo,lbufyi,ubufyi
      real, dimension (nghost,size(data,2)) :: lbufxo,ubufxo,lbufxi,ubufxi
      integer, optional :: comm
      integer :: nbufx,nbufy
!
!  Periodic boundaries in y -- communicate along y if necessary
!
      if (nprocy > 1) then
!
        lbufyo = data(l1:l2, m1:m1i)
        ubufyo = data(l1:l2,m2i:m2 )
!
        nbufy = nx * nghost
!
        call MPI_IRECV (ubufyi, nbufy, mpi_precision, yuneigh, tolowy, &
                       ioptest(comm,MPI_COMM_GRID), irecv_rq_fromuppy, mpierr)
        call MPI_IRECV (lbufyi, nbufy, mpi_precision, ylneigh, touppy, &
                       ioptest(comm,MPI_COMM_GRID), irecv_rq_fromlowy, mpierr)
        call MPI_ISEND (lbufyo, nbufy, mpi_precision, ylneigh, tolowy, &
                       ioptest(comm,MPI_COMM_GRID), isend_rq_tolowy, mpierr)
        call MPI_ISEND (ubufyo, nbufy, mpi_precision, yuneigh, touppy, &
                       ioptest(comm,MPI_COMM_GRID), isend_rq_touppy, mpierr)
!
        call MPI_WAIT (irecv_rq_fromuppy, irecv_stat_fu, mpierr)
        call MPI_WAIT (irecv_rq_fromlowy, irecv_stat_fl, mpierr)
!
        data(l1:l2,   1:m1-1) = lbufyi
        data(l1:l2,m2+1:    ) = ubufyi
!
        call MPI_WAIT (isend_rq_tolowy, isend_stat_tl, mpierr)
        call MPI_WAIT (isend_rq_touppy, isend_stat_tu, mpierr)
!
      else
!
        data(l1:l2,   1:m1-1) = data(l1:l2,m2i:m2 )
        data(l1:l2,m2+1:    ) = data(l1:l2, m1:m1i)
!
      endif
!
!  Periodic boundaries in x
!
      if (nprocx > 1) then
!
        lbufxo = data( l1:l1i,:)
        ubufxo = data(l2i:l2 ,:)
!
        nbufx = nghost * size(data,2)
!
        call MPI_IRECV (ubufxi, nbufx, mpi_precision, xuneigh, tolowx, &
                       ioptest(comm,MPI_COMM_GRID), irecv_rq_fromuppx, mpierr)
        call MPI_IRECV (lbufxi, nbufx, mpi_precision, xlneigh, touppx, &
                       ioptest(comm,MPI_COMM_GRID), irecv_rq_fromlowx, mpierr)
        call MPI_ISEND (lbufxo, nbufx, mpi_precision, xlneigh, tolowx, &
                       ioptest(comm,MPI_COMM_GRID), isend_rq_tolowx, mpierr)
        call MPI_ISEND (ubufxo, nbufx, mpi_precision, xuneigh, touppx, &
                       ioptest(comm,MPI_COMM_GRID), isend_rq_touppx, mpierr)
!
        call MPI_WAIT (irecv_rq_fromuppx, irecv_stat_fu, mpierr)
        call MPI_WAIT (irecv_rq_fromlowx, irecv_stat_fl, mpierr)
!
        data(   1:l1-1,:) = lbufxi
        data(l2+1:    ,:) = ubufxi
!
        call MPI_WAIT (isend_rq_tolowx, isend_stat_tl, mpierr)
        call MPI_WAIT (isend_rq_touppx, isend_stat_tu, mpierr)
!
      else
!
        data(   1:l1-1,:) = data(l2i:l2 ,:)
        data(l2+1:    ,:) = data( l1:l1i,:)
!
      endif
!
    endsubroutine communicate_xy_ghosts
!***********************************************************************
    subroutine fill_zghostzones_3vec(vec,ivar,comm)
!
!  Fills z-direction ghostzones of (mz,3)-array vec depending on the number of
!  processors in z-direction.
!
!  The three components of vec are supposed to be subject to the same
!  z-boundary condiitons like the variables
!  ivar, ivar+1, ivar+2
!
!   18-oct-2009/MR: Coded
!
      real, dimension(mz,3), intent(inout) :: vec
      integer, intent(in)                  :: ivar
!
      integer                    :: nbuf, j
      real, dimension (nghost,3) :: lbufi,ubufi,lbufo,ubufo
      integer, optional :: comm
!
      if (nprocz>1) then
!
        lbufo = vec(n1:n1i,:)                        !!(lower z-zone)
        ubufo = vec(n2i:n2,:)                        !!(upper z-zone)
!
        nbuf=nghost*3
!
        call MPI_IRECV(ubufi,nbuf,mpi_precision, &
                       zuneigh,tolowz,ioptest(comm,MPI_COMM_GRID),irecv_rq_fromuppz,mpierr)
        call MPI_IRECV(lbufi,nbuf,mpi_precision, &
                       zlneigh,touppz,ioptest(comm,MPI_COMM_GRID),irecv_rq_fromlowz,mpierr)
!
        call MPI_ISEND(lbufo,nbuf,mpi_precision, &
                       zlneigh,tolowz,ioptest(comm,MPI_COMM_GRID),isend_rq_tolowz,mpierr)
        call MPI_ISEND(ubufo,nbuf,mpi_precision, &
                       zuneigh,touppz,ioptest(comm,MPI_COMM_GRID),isend_rq_touppz,mpierr)
!
        call MPI_WAIT(irecv_rq_fromuppz,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowz,irecv_stat_fl,mpierr)
!
        do j=1,3
!
          if (.not. lfirst_proc_z .or. bcz12(j-1+ivar,1)=='p') &
            vec(1:n1-1,j)=lbufi(:,j)
!
!  Read from buffer in lower ghostzones.
!
          if (.not. llast_proc_z .or. bcz12(j-1+ivar,2)=='p') &
            vec(n2+1:mz,j)=ubufi(:,j)
!
!  Read from buffer in upper ghostzones.
!
        enddo
!
        call MPI_WAIT(isend_rq_tolowz,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppz,isend_stat_tu,mpierr)
!
      else
!
        do j=1,3
          if ( bcz12(ivar+j-1,1)=='p' ) then
            vec(1   :n1-1     ,j) = vec(n2i:n2 ,j)
            vec(n2+1:n2+nghost,j) = vec(n1 :n1i,j)
          endif
        enddo
!
      endif
!
    endsubroutine fill_zghostzones_3vec
!***********************************************************************
    subroutine sum_xy(in, out)
!
!  Sum up 0D data in the xy-plane and distribute back the sum.
!  This routine needs only to be called from all processors in an xy-plane.
!  Several xy-planes can call this routine at once.
!
!  19-jan-2011/Bourdin.KIS: coded
!
      real, intent(in) :: in
      real, intent(out) :: out
!
      real :: buffer, sum
      integer :: px, py, partner
      integer, parameter :: tag=114
!
      if (lfirst_proc_xy) then
        ! initialize sum with the local data
        sum = in
        ! collect and sum up the remote data
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = find_proc(px,py,ipz)
            if (iproc == partner) cycle
            call mpirecv_real (buffer, partner, tag)
            sum = sum + buffer
          enddo
        enddo
      else
        ! send data to collector
        call mpisend_real (in, find_proc(0,0,ipz), tag)
        sum = 0.0
      endif
!
      ! distribute back the sum
      call distribute_xy (out, sum)
!
    endsubroutine sum_xy
!***********************************************************************
    subroutine distribute_xy_0D(out, in, source_proc,comm)
!
!  This routine distributes a scalar on the source processor
!  to all processors in the xy-plane.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  25-jan-2012/Bourdin.KIS: coded
!
      real, intent(out) :: out
      real, intent(in), optional :: in
      integer, intent(in), optional :: source_proc
!
      integer :: px, py, broadcaster, partner
      integer, parameter :: ytag=115
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer, optional :: comm
!
      broadcaster = find_proc(ioptest(source_proc,0),0,ipz)
!
      if (iproc == broadcaster) then
        ! distribute the data
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = find_proc(px,py,ipz)
            if (iproc == partner) then
              ! data is local
              out = in
            else
              ! send to partner
              call MPI_SEND (in, 1, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
            endif
          enddo
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, 1, mpi_precision, broadcaster, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine distribute_xy_0D
!***********************************************************************
    subroutine distribute_xy_2D(out, in, source_proc, comm)
!
!  This routine divides a large array of 2D data on the source processor
!  and distributes it to all processors in the xy-plane.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  08-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(out) :: out
      real, dimension(:,:), intent(in), optional :: in
      integer, intent(in), optional :: source_proc
      integer, optional :: comm
!
      integer :: bnx, bny ! transfer box sizes
      integer :: px, py, broadcaster, partner, nbox
      integer, parameter :: ytag=115
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      bnx = size (out, 1)
      bny = size (out, 2)
      nbox = bnx*bny
!
      broadcaster = find_proc(ioptest(source_proc,0),0,ipz)
!
      if (iproc == broadcaster) then
        ! distribute the data
        if (bnx * nprocx /= size (in, 1)) &
            call stop_fatal ('distribute_xy_2D: input x dim must be nprocx*output', .true.)
        if (bny * nprocy /= size (in, 2)) &
            call stop_fatal ('distribute_xy_2D: input y dim must be nprocy*output', .true.)
!
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = find_proc(px,py,ipz)
            if (iproc /= partner) then
              ! send to partner
              out = in(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny)
              call MPI_SEND (out, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
            endif
          enddo
        enddo
        ! copy local data
        out = in(ipx*bnx+1:(ipx+1)*bnx,ipy*bny+1:(ipy+1)*bny)
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, mpi_precision, broadcaster, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine distribute_xy_2D
!***********************************************************************
    subroutine distribute_yz_3D(out, in, comm)
!
!  This routine divides a large array of 3D data on the source processor
!  and distributes it to all processors in the yz-plane.
!
!  07-oct-2021/MR: coded
!
      real, dimension(:,:,:), intent(out):: out
      real, dimension(:,:,:), intent(in) :: in
      integer, optional :: comm
!
      integer :: bnx, bny, bnz ! transfer box sizes
      integer :: pz, py, partner, nbox
      integer, parameter :: ytag=115
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      bnx = size (out, 1)
      bny = size (out, 2)
      bnz = size (out, 3)
!
      if (lfirst_proc_yz) then

        if (bnx /= size (in, 1)) &
            call stop_fatal ('distribute_yz_3D: x dim must equal between in and out', .true.)
        if (bny * nprocy /= size (in, 2)) &
            call stop_fatal ('distribute_yz_3D: input y dim must be nprocy*output', .true.)
        if (bnz * nprocz /= size (in, 3)) &
            call stop_fatal ('distribute_yz_3D: input z dim must be nprocz*output', .true.)
!
!  Distribute the data.
!
        nbox = bnx*bny*bnz
        do pz = 0, nprocz-1
          do py = 0, nprocy-1
            if (lprocz_slowest) then
              partner = py + nprocy*pz
            else
              partner = pz + nprocz*py
            endif
            if (partner==0) then
              ! copy local data on broadcaster
              out = in(:,1:bny,1:bnz)
            else
              ! send to partner
              out = in(:,py*bny+1:(py+1)*bny,pz*bnz+1:(pz+1)*bnz)
              call MPI_SEND (out, nbox, mpi_precision, partner, ytag, MPI_COMM_YZPLANE, mpierr)
            endif
          enddo
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, mpi_precision, 0, ytag, MPI_COMM_YZPLANE, stat, mpierr)
      endif
!
    endsubroutine distribute_yz_3D
!***********************************************************************
    subroutine distribute_yz_4D(out, in)
!
!  This routine divides a large array of 4D data on the source processor
!  and distributes it to all processors in the yz-plane.
!
!  07-oct-2021/MR: coded
!
      real, dimension(:,:,:,:), intent(out):: out
      real, dimension(:,:,:,:), intent(in) :: in
!
      integer :: bnx, bny, bnz, bnv ! transfer box sizes
      integer :: pz, py, partner, nbox
      integer, parameter :: ytag=115
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      bnx = size (out, 1)
      bny = size (out, 2)
      bnz = size (out, 3)
      bnv = size (out, 4)
!
      if (lfirst_proc_yz) then

        if (bnx /= size (in, 1)) &
            call stop_fatal ('distribute_yz_3D: x dim must equal between in and out', .true.)
        if (bny * nprocy /= size (in, 2)) &
            call stop_fatal ('distribute_yz_3D: input y dim must be nprocy*output', .true.)
        if (bnz * nprocz /= size (in, 3)) &
            call stop_fatal ('distribute_yz_3D: input z dim must be nprocz*output', .true.)
        if (bnv /= size (in, 4)) &
            call stop_fatal ('distribute_yz_3D: var dim must equal between in and out', .true.)
!
!  Distribute the data.
!
        nbox = bnx*bny*bnz*bnv
        do pz = 0, nprocz-1
          do py = 0, nprocy-1
            if (lprocz_slowest) then
              partner = py + nprocy*pz
            else
              partner = pz + nprocz*py
            endif
            if (partner==0) then
              ! copy local data on broadcaster
              out = in(:,1:bny,1:bnz,:)
            else
              ! send to partner
              out = in(:,py*bny+1:(py+1)*bny,pz*bnz+1:(pz+1)*bnz,:)
              call MPI_SEND (out, nbox, mpi_precision, partner, ytag, MPI_COMM_YZPLANE, mpierr)
            endif
          enddo
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, mpi_precision, 0, ytag, MPI_COMM_YZPLANE, stat, mpierr)
      endif
!
    endsubroutine distribute_yz_4D
!***********************************************************************
    subroutine distribute_xy_3D(out, in, source_proc, comm)
!
!  This routine divides a large array of 3D data on the source processor
!  and distributes it to all processors in the xy-plane.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  08-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(out) :: out
      real, dimension(:,:,:), intent(in), optional :: in
      integer, intent(in), optional :: source_proc
!
      integer :: bnx, bny, bnz ! transfer box sizes
      integer :: px, py, broadcaster, partner, nbox
      integer, parameter :: ytag=115
      integer, optional :: comm
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      bnx = size (out, 1)
      bny = size (out, 2)
      bnz = size (out, 3)
      nbox = bnx*bny*bnz
!
      broadcaster = find_proc(ioptest(source_proc,0),0,ipz)
!
      if (iproc == broadcaster) then
        ! distribute the data
        if (bnx * nprocx /= size (in, 1)) &
            call stop_fatal ('distribute_xy_3D: input x dim must be nprocx*output', .true.)
        if (bny * nprocy /= size (in, 2)) &
            call stop_fatal ('distribute_xy_3D: input y dim must be nprocy*output', .true.)
        if (bnz /= size (in, 3)) &
            call stop_fatal ('distribute_xy_3D: z dim must equal between in and out', .true.)
!
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = find_proc(px,py,ipz)
            if (iproc == partner) then
              ! data is local
              out = in(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:)
            else
              ! send to partner
              call MPI_SEND (in(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:), &
                  nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
            endif
          enddo
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, mpi_precision, broadcaster, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine distribute_xy_3D
!***********************************************************************
    subroutine distribute_xy_4D(out, in, source_proc, comm)
!
!  This routine divides a large array of 4D data on the source processor
!  and distributes it to all processors in the xy-plane.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  08-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(out) :: out
      real, dimension(:,:,:,:), intent(in), optional :: in
      integer, intent(in), optional :: source_proc
      integer, optional :: comm
!
      integer :: bnx, bny, bnz, bna ! transfer box sizes
      integer :: px, py, broadcaster, partner, nbox
      integer, parameter :: ytag=115
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      bnx = size (out, 1)
      bny = size (out, 2)
      bnz = size (out, 3)
      bna = size (out, 4)
      nbox = bnx*bny*bnz*bna
!
      broadcaster = find_proc(ioptest(source_proc,0),0,ipz)
!
      if (iproc == broadcaster) then
        ! distribute the data
        if (bnx * nprocx /= size (in, 1)) &
            call stop_fatal ('distribute_xy_4D: input x dim must be nprocx*output', .true.)
        if (bny * nprocy /= size (in, 2)) &
            call stop_fatal ('distribute_xy_4D: input y dim must be nprocy*output', .true.)
        if (bnz /= size (in, 3)) &
            call stop_fatal ('distribute_xy_4D: z dim must equal between in and out', .true.)
        if (bna /= size (in, 4)) &
            call stop_fatal ('distribute_xy_4D: 4th dim must equal between in and out', .true.)
!
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = find_proc(px,py,ipz)
            if (iproc == partner) then
              ! data is local
              out = in(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:,:)
            else
              ! send to partner
              call MPI_SEND (in(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:,:), &
                  nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
            endif
          enddo
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, mpi_precision, broadcaster, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine distribute_xy_4D
!***********************************************************************
    subroutine collect_xy_0D(in, out, dest_proc, comm)
!
!  Collect 0D data from all processors in the xy-plane
!  and combine it into one large array on one destination processor.
!  'dest_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  08-jan-2011/Bourdin.KIS: coded
!
      real, intent(in) :: in
      real, dimension(:,:), intent(out), optional :: out
      integer, intent(in), optional :: dest_proc
      integer, optional :: comm
!
      integer :: px, py, collector, partner
      integer, parameter :: ytag=116
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real :: buffer
!
      collector = ipz * nprocxy
      if (present (dest_proc)) collector = collector + dest_proc
!
      if (iproc == collector) then
        ! collect the data
        if (nprocx /= size (out, 1)) &
            call stop_fatal ('collect_xy_0D: output x dim must be nprocx', .true.)
        if (nprocy /= size (out, 2)) &
            call stop_fatal ('collect_xy_0D: output y dim must be nprocy', .true.)
!
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = find_proc(px,py,ipz)
            if (iproc == partner) then
              ! data is local
              out(px+1,py+1) = in
            else
              ! receive from partner
              call MPI_RECV (buffer, 1, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
              out(px+1,py+1) = buffer
            endif
          enddo
        enddo
      else
        ! send to collector
        call MPI_SEND (in, 1, mpi_precision, collector, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine collect_xy_0D
!***********************************************************************
    subroutine collect_xy_2D(in, out, dest_proc, comm)
!
!  Collect 2D data from all processors in the xy-plane
!  and combine it into one large array on one destination processor.
!  'dest_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  08-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out), optional :: out
      integer, intent(in), optional :: dest_proc
      integer, optional :: comm
!
      integer :: bnx, bny ! transfer box sizes
      integer :: px, py, collector, partner, nbox, alloc_err
      integer, parameter :: ytag=116
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: buffer
!
      bnx = size (in, 1)
      bny = size (in, 2)
      nbox = bnx*bny
!
      collector = ipz * nprocxy
      if (present (dest_proc)) collector = collector + dest_proc
!
      if (iproc == collector) then
        ! collect the data
        if (bnx * nprocx /= size (out, 1)) &
            call stop_fatal ('collect_xy_2D: output x dim must be nprocx*input', .true.)
        if (bny * nprocy /= size (out, 2)) &
            call stop_fatal ('collect_xy_2D: output y dim must be nprocy*input', .true.)
!
        allocate (buffer(bnx,bny), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('collect_xy_2D: not enough memory for buffer!', .true.)
!
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = find_proc(px,py,ipz)
            if (iproc == partner) then
              ! data is local
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny) = in
            else
              ! receive from partner
              call MPI_RECV (buffer, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny) = buffer
            endif
          enddo
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, nbox, mpi_precision, collector, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine collect_xy_2D
!***********************************************************************
    subroutine collect_xy_3D(in, out, dest_proc, comm)
!
!  Collect 3D data from all processors in the xy-plane
!  and combine it into one large array on one destination processor.
!  'dest_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  08-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out), optional :: out
      integer, intent(in), optional :: dest_proc
      integer, optional :: comm
!
      integer :: bnx, bny, bnz ! transfer box sizes
      integer :: px, py, collector, partner, nbox, alloc_err
      integer, parameter :: ytag=116
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:), allocatable :: buffer
!
      bnx = size (in, 1)
      bny = size (in, 2)
      bnz = size (in, 3)
      nbox = bnx*bny*bnz
!
      collector = ipz * nprocxy
      if (present (dest_proc)) collector = collector + dest_proc
!
      if (iproc == collector) then
        ! collect the data
        if (bnx * nprocx /= size (out, 1)) &
            call stop_fatal ('collect_xy_3D: output x dim must be nprocx*input', .true.)
        if (bny * nprocy /= size (out, 2)) &
            call stop_fatal ('collect_xy_3D: output y dim must be nprocy*input', .true.)
!
        allocate (buffer(bnx,bny,bnz), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('collect_xy_3D: not enough memory for buffer!', .true.)
!
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = find_proc(px,py,ipz)
            if (iproc == partner) then
              ! data is local
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:) = in
            else
              ! receive from partner
              call MPI_RECV (buffer, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:) = buffer
            endif
          enddo
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, nbox, mpi_precision, collector, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine collect_xy_3D
!***********************************************************************
    subroutine collect_xy_4D(in, out, dest_proc, comm)
!
!  Collect 4D data from all processors in the xy-plane
!  and combine it into one large array on one destination processor.
!  'dest_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  08-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out), optional :: out
      integer, intent(in), optional :: dest_proc
      integer, optional :: comm
!
      integer :: bnx, bny, bnz, bna ! transfer box sizes
      integer :: px, py, collector, partner, nbox, alloc_err
      integer, parameter :: ytag=116
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: buffer
!
      bnx = size (in, 1)
      bny = size (in, 2)
      bnz = size (in, 3)
      bna = size (in, 4)
      nbox = bnx*bny*bnz*bna
!
      collector = ipz * nprocxy
      if (present (dest_proc)) collector = collector + dest_proc
!
      if (iproc == collector) then
        ! collect the data
        if (bnx * nprocx /= size (out, 1)) &
            call stop_fatal ('collect_xy_4D: output x dim must be nprocx*input', .true.)
        if (bny * nprocy /= size (out, 2)) &
            call stop_fatal ('collect_xy_4D: output y dim must be nprocy*input', .true.)
        if (bnz /= size (out, 3)) &
            call stop_fatal ('collect_xy_4D: z dim must equal between in and out', .true.)
        if (bna /= size (out, 4)) &
            call stop_fatal ('collect_xy_4D: 4th dim must equal between in and out', .true.)
!
        allocate (buffer(bnx,bny,bnz,bna), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('collect_xy_4D: not enough memory for buffer!', .true.)
!
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = find_proc(px,py,ipz)
            if (iproc == partner) then
              ! data is local
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:,:) = in
            else
              ! receive from partner
              call MPI_RECV (buffer, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:,:) = buffer
            endif
          enddo
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, nbox, mpi_precision, collector, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine collect_xy_4D
!***********************************************************************
    subroutine distribute_z_3D(out, in, source_proc, comm)
!
!  This routine divides a large array of 3D data on the source processor
!  and distributes it to all processors in the z-direction.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding z-direction (Default: 0, equals lfirst_proc_z).
!
!  08-mar-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(out) :: out
      real, dimension(:,:,:), intent(in), optional :: in
      integer, intent(in), optional :: source_proc
      integer, optional :: comm
!
      integer :: bnx, bny, bnz ! transfer box sizes
      integer :: pz, broadcaster, partner, nbox
      integer, parameter :: ytag=117
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      bnx = size (out, 1)
      bny = size (out, 2)
      bnz = size (out, 3)
      nbox = bnx*bny*bnz
!
      broadcaster = find_proc(ipx,ipy,ioptest(source_proc,0))
!
      if (iproc == broadcaster) then
        ! distribute the data
        if (bnx /= size (in, 1)) &
            call stop_fatal ('distribute_z_4D: x dim must be equal between in and out', .true.)
        if (bny /= size (in, 2)) &
            call stop_fatal ('distribute_z_4D: y dim must be equal between in and out', .true.)
        if (bnz * nprocz /= size (in, 3)) &
            call stop_fatal ('distribute_z_4D: input z dim must be nprocz*output', .true.)
!
        do pz = 0, nprocz-1
          partner = find_proc(ipx,ipy,pz)
          if (iproc == partner) then
            ! data is local
            out = in(:,:,pz*bnz+1:(pz+1)*bnz)
          else
            ! send to partner
            call MPI_SEND (in(:,:,pz*bnz+1:(pz+1)*bnz), nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, mpi_precision, broadcaster, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine distribute_z_3D
!***********************************************************************
    subroutine distribute_z_4D(out, in, source_proc, comm)
!
!  This routine divides a large array of 4D data on the source processor
!  and distributes it to all processors in the z-direction.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding z-direction (Default: 0, equals lfirst_proc_z).
!
!  08-mar-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(out) :: out
      real, dimension(:,:,:,:), intent(in), optional :: in
      integer, intent(in), optional :: source_proc
      integer, optional :: comm
!
      integer :: bnx, bny, bnz, bna ! transfer box sizes
      integer :: pz, broadcaster, partner, nbox
      integer, parameter :: ytag=117
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
!
      bnx = size (out, 1)
      bny = size (out, 2)
      bnz = size (out, 3)
      bna = size (out, 4)
      nbox = bnx*bny*bnz*bna
!
      broadcaster = find_proc(ipx,ipy,ioptest(source_proc,0))
!
      if (iproc == broadcaster) then
        ! distribute the data
        if (bnx /= size (in, 1)) &
            call stop_fatal ('distribute_z_4D: x dim must be equal between in and out', .true.)
        if (bny /= size (in, 2)) &
            call stop_fatal ('distribute_z_4D: y dim must be equal between in and out', .true.)
        if (bnz * nprocz /= size (in, 3)) &
            call stop_fatal ('distribute_z_4D: input z dim must be nprocz*output', .true.)
        if (bna /= size (in, 4)) &
            call stop_fatal ('distribute_z_4D: 4th dim must equal between in and out', .true.)
!
        do pz = 0, nprocz-1
          partner = find_proc(ipx,ipy,pz)
          if (iproc == partner) then
            ! data is local
            out = in(:,:,pz*bnz+1:(pz+1)*bnz,:)
          else
            ! send to partner
            call MPI_SEND (in(:,:,pz*bnz+1:(pz+1)*bnz,:), nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, mpi_precision, broadcaster, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine distribute_z_4D
!***********************************************************************
    subroutine collect_z_3D(in, out, dest_proc, comm)
!
!  Collect 3D data from all processors in the z-direction
!  and combine it into one large array on one destination processor.
!  'dest_proc' is the iproc number relative to the first processor
!  in the corresponding z-direction (Default: 0, equals lfirst_proc_z).
!
!  08-mar-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out), optional :: out
      integer, intent(in), optional :: dest_proc
      integer, optional :: comm
!
      integer :: bnx, bny, bnz ! transfer box sizes
      integer :: pz, collector, partner, nbox, alloc_err
      integer, parameter :: ytag=118
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:), allocatable :: buffer
!
      bnx = size (in, 1)
      bny = size (in, 2)
      bnz = size (in, 3)
      nbox = bnx*bny*bnz
!
      collector = ipx + ipy*nprocx
      if (present (dest_proc)) collector = collector + dest_proc*nprocxy
!
      if (iproc == collector) then
        ! collect the data
        if (bnx /= size (out, 1)) &
            call stop_fatal ('collect_z_3D: x dim must equal between in and out', .true.)
        if (bny /= size (out, 2)) &
            call stop_fatal ('collect_z_3D: y dim must equal between in and out', .true.)
        if (bnz * nprocz /= size (out, 3)) &
            call stop_fatal ('collect_z_3D: output z dim must be nprocz*input', .true.)
!
        allocate (buffer(bnx,bny,bnz), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('collect_z_3D: not enough memory for buffer!', .true.)
!
        do pz = 0, nprocz-1
          partner = find_proc(ipx,ipy,pz)
          if (iproc == partner) then
            ! data is local
            out(:,:,pz*bnz+1:(pz+1)*bnz) = in
          else
            ! receive from partner
            call MPI_RECV (buffer, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            out(:,:,pz*bnz+1:(pz+1)*bnz) = buffer
          endif
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, nbox, mpi_precision, collector, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine collect_z_3D
!***********************************************************************
    subroutine collect_z_4D(in, out, dest_proc, comm)
!
!  Collect 4D data from all processors in the z-direction
!  and combine it into one large array on one destination processor.
!  'dest_proc' is the iproc number relative to the first processor
!  in the corresponding z-direction (Default: 0, equals lfirst_proc_z).
!
!  08-mar-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out), optional :: out
      integer, intent(in), optional :: dest_proc
      integer, optional :: comm
!
      integer :: bnx, bny, bnz, bna ! transfer box sizes
      integer :: pz, collector, partner, nbox, alloc_err
      integer, parameter :: ytag=118
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: buffer
!
      bnx = size (in, 1)
      bny = size (in, 2)
      bnz = size (in, 3)
      bna = size (in, 4)
      nbox = bnx*bny*bnz*bna
!
      collector = ipx + ipy*nprocx
      if (present (dest_proc)) collector = collector + dest_proc*nprocxy
!
      if (iproc == collector) then
        ! collect the data
        if (bnx /= size (out, 1)) &
            call stop_fatal ('collect_z_4D: x dim must equal between in and out', .true.)
        if (bny /= size (out, 2)) &
            call stop_fatal ('collect_z_4D: y dim must equal between in and out', .true.)
        if (bnz * nprocz /= size (out, 3)) &
            call stop_fatal ('collect_z_4D: output z dim must be nprocz*input', .true.)
        if (bna /= size (out, 4)) &
            call stop_fatal ('collect_z_4D: 4th dim must equal between in and out', .true.)
!
        allocate (buffer(bnx,bny,bnz,bna), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('collect_z_4D: not enough memory for buffer!', .true.)
!
        do pz = 0, nprocz-1
          partner = find_proc(ipx,ipy,pz)
          if (iproc == partner) then
            ! data is local
            out(:,:,pz*bnz+1:(pz+1)*bnz,:) = in
          else
            ! receive from partner
            call MPI_RECV (buffer, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            out(:,:,pz*bnz+1:(pz+1)*bnz,:) = buffer
          endif
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, nbox, mpi_precision, collector, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine collect_z_4D
!***********************************************************************
    subroutine globalize_xy(in, out, dest_proc, source_pz)
!
!  Globalizes local 4D data first along the x, then along the y-direction to
!  the destination processor. The local data is supposed to include the ghost
!  cells. Inner ghost layers are cut away during the combination of the data.
!  'dest_proc' is the destination iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!  'source_pz' specifies the source pz-layer (Default: ipz).
!
!  23-Apr-2012/Bourdin.KIS: adapted from non-torus-type globalize_xy
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out), optional :: out
      integer, intent(in), optional :: dest_proc, source_pz
!
      integer :: bnx, bny, bnz, bna ! transfer box sizes
      integer :: cnx, cny ! transfer box sizes minus ghosts
      integer :: rny ! y-row box size
      integer :: px, py, pz, collector, partner, alloc_err
      integer(KIND=ikind8) :: nbox, nrow
      integer :: x_add, x_sub, y_add, y_sub
      integer, parameter :: xtag=123, ytag=124
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: buffer, y_row
!
      bnx = size (in, 1)
      bny = size (in, 2)
      bnz = size (in, 3)
      bna = size (in, 4)
      nbox = bnx*bny*bnz*bna
      cnx = bnx - 2*nghost
      cny = bny - 2*nghost
      rny = cny*nprocy + 2*nghost
      nrow = bnx*rny*bnz*bna
!
      collector = ipz * nprocxy
      if (present (dest_proc)) collector = collector + dest_proc
      pz = ioptest(source_pz,ipz)
!
      if (iproc == collector) then
        if (cnx * nprocx + 2*nghost /= size (out, 1)) &
            call stop_fatal ('globalize_xy: output x dim must be nprocx*input minus inner ghosts', .true.)
        if (cny * nprocy + 2*nghost /= size (out, 2)) &
            call stop_fatal ('globalize_xy: output y dim must be nprocy*input minus inner ghosts', .true.)
        if (bnz /= size (out, 3)) &
            call stop_fatal ('globalize_xy: z dim must equal between in and out', .true.)
        if (bna /= size (out, 4)) &
            call stop_fatal ('globalize_xy: 4th dim must equal between in and out', .true.)
      endif
!
      if (lfirst_proc_y .and. (ipz == pz)) then
        ! collect the data of each y-row
        allocate (buffer(bnx,bny,bnz,bna), y_row(bnx,rny,bnz,bna), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('globalize_xy: not enough memory for buffer and y_row!', .true.)
!
        do py = 0, nprocy-1
          partner = find_proc(ipx,py,pz)
          y_add = nghost
          y_sub = nghost
          if (py == 0) y_add = 0
          if (py == nprocy-1) y_sub = 0
          if (iproc == partner) then
            ! data is local
            y_row(:,py*cny+1+y_add:py*cny+my,:,:) = in(:,1+y_add:my,:,:)
          else
            ! receive from y-row partner
            call mpirecv_real_arr_huge(buffer,nbox,partner,ytag)
            ! old version: call MPI_RECV (buffer, nbox, mpi_precision, partner, ytag, MPI_COMM_GRID, stat, mpierr)
            y_row(:,py*cny+1+y_add:py*cny+my-y_sub,:,:) = buffer(:,1+y_add:my-y_sub,:,:)
          endif
        enddo
!
        deallocate (buffer)
!
        if (iproc /= collector) then
          ! send to collector
          call mpisend_real_arr_huge(y_row,nrow,collector,xtag)
          ! old version: call MPI_SEND (y_row, nrow, mpi_precision, collector, xtag, MPI_COMM_GRID, mpierr)
          deallocate (y_row)
        endif
!
      elseif (ipz == pz) then
        ! send to collector of the y-row (lfirst_proc_y)
        partner = find_proc(ipx,0,ipz)
        call mpisend_real_arr_huge(in,nbox,partner,ytag)
        ! old version: call MPI_SEND (in, nbox, mpi_precision, partner, ytag, MPI_COMM_GRID, mpierr)
      endif
!
      if (iproc == collector) then
        ! collect the y-rows into global data
        allocate (buffer(bnx,rny,bnz,bna), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('globalize_xy: not enough memory for buffer!', .true.)
!
        do px = 0, nprocx-1
          partner = find_proc(px,ipy,pz)
          x_add = nghost
          x_sub = nghost
          if (px == 0) x_add = 0
          if (px == nprocx-1) x_sub = 0
          if (iproc == partner) then
            ! y-row is local
            out(px*cnx+1+x_add:px*cnx+mx,:,:,:) = y_row(1+x_add:mx,:,:,:)
            deallocate (y_row)
          else
            ! receive from partner
            call mpirecv_real_arr_huge(buffer,nrow,partner,xtag)
            ! old version: call MPI_RECV (buffer, nrow, mpi_precision, partner, xtag, MPI_COMM_GRID, stat, mpierr)
            out(px*cnx+1+x_add:px*cnx+mx-x_sub,:,:,:) = buffer(1+x_add:mx-x_sub,:,:,:)
          endif
        enddo
!
        deallocate (buffer)
      endif
!
    endsubroutine globalize_xy
!***********************************************************************
    recursive subroutine localize_xy(out, in, source_proc, dest_pz)
!
!  Localizes global 4D data first along the y, then along the x-direction to
!  the destination processor. The global data is supposed to include the outer
!  ghost layers. The returned data will include inner ghost layers.
!  'source_proc' is the source iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!  'dest_pz' specifies the destination pz-layer (Default: ipz).
!
!  23-Apr-2012/Bourdin.KIS: adapted from non-torus-type localize_xy
!
      real, dimension(:,:,:,:), intent(out) :: out
      real, dimension(:,:,:,:), intent(in), optional :: in
      integer, intent(in), optional :: source_proc, dest_pz
!
      integer :: bnx, bny, bnz, bna ! transfer box sizes
      integer :: cnx, cny ! transfer box sizes minus ghosts
      integer :: gnx, gny ! x- and y global data array sizes
      integer :: rnx, rny ! x- and y-row box sizes
      integer :: px, py, pz, broadcaster, partner, alloc_err
      integer(KIND=ikind8) :: nbox,nrow
      integer, parameter :: xtag=125, ytag=126
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: y_row, buffer, extended
!
      bnx = size (out, 1)
      bny = size (out, 2)
      bnz = size (out, 3)
      bna = size (out, 4)
      nbox = bnx*bny*bnz*bna
      cnx = bnx - 2*nghost
      cny = bny - 2*nghost
      gnx = cnx*nprocx
      gny = cny*nprocy
      rnx = gnx + 2*nghost
      rny = gny + 2*nghost
      nrow = bnx*rny*bnz*bna
!
      broadcaster = ipz * nprocxy
      if (present (source_proc)) broadcaster = broadcaster + source_proc
      pz = ipz
      if (present (dest_pz)) pz = dest_pz
!
      if (iproc == broadcaster) then
        if ((gnx == size (in, 1)) .and. (gny == size (in, 2))) then
          ! add outer ghost layers
          allocate (extended(rnx,rny,bnz,bna), stat=alloc_err)
          if (alloc_err > 0) call stop_fatal ('localize_xy: not enough memory for outer ghost cells extension!', .true.)
          extended(nghost+1:nghost+gnx,nghost+1:nghost+gnx,:,:) = in
          if (lperi(1)) then
            extended(1:nghost,:,:,:) = in(gnx-nghost+1:gnx,:,:,:)
            extended(rnx-nghost+1:rnx,:,:,:) = in(1:nghost,:,:,:)
          else
            write (*,*) "WARNING: localize_xy: extending non-periodic outer ghost cells with zeros!"
            extended(1:nghost,:,:,:) = 0.0
            extended(rnx-nghost+1:rnx,:,:,:) = 0.0
          endif
          if (lperi(2)) then
            extended(:,1:nghost,:,:) = in(:,gny-nghost+1:gny,:,:)
            extended(:,rny-nghost+1:rny,:,:) = in(:,1:nghost,:,:)
          else
            write (*,*) "WARNING: localize_xy: extending non-periodic outer ghost cells with zeros!"
            extended(:,1:nghost,:,:) = 0.0
            extended(:,rny-nghost+1:rny,:,:) = 0.0
          endif
          call localize_xy(out, extended, broadcaster - ipz * nprocxy, pz)
          return
        endif
        if (cnx * nprocx + 2*nghost /= size (in, 1)) &
            call stop_fatal ('localize_xy: input x dim must be nprocx*output minus inner ghosts', .true.)
        if (cny * nprocy + 2*nghost /= size (in, 2)) &
            call stop_fatal ('localize_xy: input y dim must be nprocy*output minus inner ghosts', .true.)
        if (bnz /= size (in, 3)) &
            call stop_fatal ('localize_xy: z dim must equal between in and out', .true.)
        if (bna /= size (in, 4)) &
            call stop_fatal ('localize_xy: 4th dim must equal between in and out', .true.)
      endif
!
      if (ipz == pz) then
        allocate (y_row(bnx,rny,bnz,bna), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('localize_xy: not enough memory for y_row!', .true.)
      endif
!
      if (iproc == broadcaster) then
        allocate (buffer(bnx,rny,bnz,bna), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('localize_xy: not enough memory for buffer!', .true.)
        ! distribute the y-rows
        do px = 0, nprocx-1
          partner = find_proc(px,0,pz)
          if (iproc == partner) then
            ! data is local
            y_row = in(px*cnx+1:px*cnx+mx,:,:,:)
          else
            ! send to partner
            buffer = in(px*cnx+1:px*cnx+mx,:,:,:)

            call mpisend_real_arr_huge(buffer,nrow,partner,xtag)
            !old version: call MPI_SEND (buffer, int(nrow), mpi_precision, partner, xtag, MPI_COMM_GRID, mpierr)
          endif
        enddo
        deallocate (buffer)
      endif
!
      if (lfirst_proc_y .and. (ipz == pz)) then
        if (iproc /= broadcaster) then
          ! receive y-row from broadcaster

          call mpirecv_real_arr_huge(y_row,nrow,broadcaster,xtag)
          !old version: call MPI_RECV (y_row, int(nrow), mpi_precision, broadcaster, xtag, MPI_COMM_GRID, stat, mpierr)
        endif
!
        allocate (buffer(bnx,bny,bnz,bna), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('localize_xy: not enough memory for buffer!', .true.)
        ! distribute the data along the y-direction
        do py = 0, nprocy-1
          partner = find_proc(ipx,py,pz)
          if (iproc == partner) then
            ! data is local
            out = y_row(:,py*cny+1:py*cny+my,:,:)
          else
            ! send to partner
            buffer = y_row(:,py*cny+1:py*cny+my,:,:)

            call mpisend_real_arr_huge(buffer,nbox,partner,ytag)
            !old version: call MPI_SEND (buffer, int(nbox), mpi_precision, partner, ytag, MPI_COMM_GRID, mpierr)
          endif
        enddo
        deallocate (buffer)
      elseif (ipz == pz) then
        ! receive local data from y-row partner (lfirst_proc_y)
        partner = find_proc(ipx,0,ipz)

        call mpirecv_real_arr_huge(out,nbox,partner,ytag)
        !old version: call MPI_RECV (out, int(nbox), mpi_precision, partner, ytag, MPI_COMM_GRID, stat, mpierr)
      endif
!
      if (ipz == pz) deallocate (y_row)
!
    endsubroutine localize_xy
!***********************************************************************
    subroutine globalize_z(in, out, dest_proc, comm)
!
!  Globalizes local 1D data in the z-direction to the destination processor.
!  The local data is supposed to include the ghost cells.
!  Inner ghost layers are cut away during the combination of the data.
!  'dest_proc' is the destination ipz-layer number relative to the first
!  processor in the z-direction (Default: 0, equals lfirst_proc_z).
!
!  13-aug-2011/Bourdin.KIS: coded
!
      real, dimension(mz), intent(in) :: in
      real, dimension(mzgrid), intent(out), optional :: out
      integer, intent(in), optional :: dest_proc
      integer, optional :: comm
!
      integer :: pz, z_add, collector, partner, alloc_err
      integer, parameter :: ytag=119
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:), allocatable :: buffer
!
      collector = ipx + ipy * nprocx
      if (present (dest_proc)) collector = collector + dest_proc * nprocxy
!
      if (iproc == collector) then
        ! collect the data
        allocate (buffer(mz), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('globalize_z: not enough memory for buffer!', .true.)
!
        do pz = 0, nprocz-1
          partner = find_proc(ipx,ipy,pz)
          z_add = nghost
          if (pz == 0) z_add = 0
          if (iproc == partner) then
            ! data is local
            out(pz*nz+1+z_add:pz*nz+mz) = in(1+z_add:mz)
          else
            ! receive from partner
            call MPI_RECV (buffer, mz, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            out(pz*nz+1+z_add:pz*nz+mz) = buffer(1+z_add:mz)
          endif
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, mz, mpi_precision, collector, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
    endsubroutine globalize_z
!***********************************************************************
    subroutine localize_z(out, in, source_proc, comm)
!
!  Localizes global 1D data to all processors along the z-direction.
!  The global data is supposed to include the outer ghost layers.
!  The returned data will include inner ghost layers.
!  'source_proc' is the source ipz-layer number relative to the first
!  processor in the z-direction (Default: 0, equals lfirst_proc_z).
!
!  11-Feb-2012/Bourdin.KIS: coded
!
      real, dimension(mz), intent(out) :: out
      real, dimension(mzgrid), intent(in) :: in
      integer, intent(in), optional :: source_proc
      integer, optional :: comm
!
      integer :: pz, broadcaster, partner
      integer, parameter :: ytag=120
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      broadcaster = find_proc(ipx,ipy,ioptest(source_proc,0))
!
      if (iproc == broadcaster) then
        ! collect the data
        do pz = 0, nprocz-1
          partner = find_proc(ipx,ipy,pz)
          if (iproc == partner) then
            ! data is local
            out = in(pz*nz+1:pz*nz+mz)
          else
            ! send to partner
            call MPI_SEND (in(pz*nz+1:pz*nz+mz), mz, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, mz, mpi_precision, broadcaster, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
      endif
!
    endsubroutine localize_z
!***********************************************************************
    subroutine distribute_to_pencil_xy_2D(in, out, broadcaster, comm)
!
!  Distribute data to several processors and reform into pencil shape.
!  This routine divides global data and distributes it in the xy-plane.
!
!  22-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
      integer, intent(in) :: broadcaster
      integer, optional :: comm
!
      integer :: bnx, bny ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ytag=113
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: buffer
!
      if ((nprocx == 1) .and. (nprocy == 1)) then
        out = in
        return
      endif
!
      bnx = size (in, 1)
      bny = size (in, 2) / nprocxy
      nbox = bnx*bny
!
      if (mod (size (in, 2), nprocxy) /= 0) &
          call stop_fatal ('distribute_to_pencil_xy_2D: input y dim must be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
!
      if ((size (out, 1) /= bnx) .or. ((size (out, 2) /= bny))) &
          call stop_fatal ('distribute_to_pencil_xy_2D: output array size mismatch /= bnx,bny', lfirst_proc_xy)
!
      allocate (buffer(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('distribute_to_pencil_xy_2D: not enough memory for buffer!', .true.)
!
      if (iproc == broadcaster) then
        do ibox = 0, nprocxy-1
          partner = find_proc(modulo(ibox,nprocx),ibox/nprocx,ipz)
          if (iproc == partner) then
            ! data is local
            out = in(:,bny*ibox+1:bny*(ibox+1))
          else
            ! send to partner
            buffer = in(:,bny*ibox+1:bny*(ibox+1))
            call MPI_SEND (buffer, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (buffer, nbox, mpi_precision, broadcaster, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
        out = buffer
      endif
!
      deallocate (buffer)
!
    endsubroutine distribute_to_pencil_xy_2D
!***********************************************************************
    subroutine collect_from_pencil_xy_2D(in, out, collector, comm)
!
!  Collect 2D data from several processors and combine into global shape.
!  This routine collects 2D pencil shaped data distributed in the xy-plane.
!
!  22-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
      integer, intent(in) :: collector
      integer, optional :: comm
!
      integer :: bnx, bny ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ytag=114
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: buffer
!
      if ((nprocx == 1) .and. (nprocy == 1)) then
        out = in
        return
      endif
!
      bnx = size (out, 1)
      bny = size (out, 2) / nprocxy
      nbox = bnx*bny
!
      if (mod (size (out, 2), nprocxy) /= 0) &
          call stop_fatal ('collect_from_pencil_xy_2D: output y dim must be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
!
      if ((size (in, 1) /= bnx) .or. ((size (in, 2) /= bny))) &
          call stop_fatal ('collect_from_pencil_xy_2D: input array size mismatch /= bnx,bny', lfirst_proc_xy)
!
      allocate (buffer(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('collect_from_pencil_xy_2D: not enough memory for buffer!', .true.)
!
      if (iproc == collector) then
        do ibox = 0, nprocxy-1
          partner = find_proc(modulo(ibox,nprocx),ibox/nprocx,ipz)
          if (iproc == partner) then
            ! data is local
            out(:,bny*ibox+1:bny*(ibox+1)) = in
          else
            ! receive from partner
            call MPI_RECV (buffer, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            out(:,bny*ibox+1:bny*(ibox+1)) = buffer
          endif
        enddo
      else
        ! send to collector
        buffer = in
        call MPI_SEND (buffer, nbox, mpi_precision, collector, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
      endif
!
      deallocate (buffer)
!
    endsubroutine collect_from_pencil_xy_2D
!***********************************************************************
    subroutine remap_to_pencil_x(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 1D arrays in x only for nprocx>1.
!
!  08-dec-2010/Bourdin.KIS: coded
!
      real, dimension(nx), intent(in) :: in
      real, dimension(nxgrid), intent(out) :: out
!
      integer :: ibox, partner
      integer, parameter :: ltag=102, utag=103
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(nx) :: recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      do ibox = 0, nprocx-1
        partner = find_proc(ibox,ipy,ipz)
        if (iproc == partner) then
          ! data is local
          out(nx*ibox+1:nx*(ibox+1)) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nx, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nx, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nx, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (in, nx, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(nx*ibox+1:nx*(ibox+1)) = recv_buf
        endif
      enddo
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_x
!***********************************************************************
    subroutine unmap_from_pencil_x(in, out)
!
!  Unmaps pencil shaped 1D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocx>1.
!
!  08-dec-2010/Bourdin.KIS: coded
!
      real, dimension(nxgrid), intent(in) :: in
      real, dimension(nx), intent(out) :: out
!
      !$omp workshare
      out = in(nx*ipx+1:nx*(ipx+1))
      !$omp end workshare
!
    endsubroutine unmap_from_pencil_x
!***********************************************************************
    subroutine remap_to_pencil_y_1D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 1D arrays in y only for nprocy>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(ny), intent(in) :: in
      real, dimension(nygrid), intent(out) :: out
!
      integer :: ibox, partner
      integer, parameter :: ltag=102, utag=103
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(ny) :: recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      do ibox = 0, nprocy-1
        partner = find_proc(ipx,ibox,ipz)
        if (iproc == partner) then
          ! data is local
          out(ny*ibox+1:ny*(ibox+1)) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, ny, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, ny, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, ny, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (in, ny, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(ny*ibox+1:ny*(ibox+1)) = recv_buf
        endif
      enddo
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_y_1D
!***********************************************************************
    subroutine remap_to_pencil_y_2D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 2D arrays in y only for nprocy>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(nx,ny), intent(in) :: in
      real, dimension(nx,nygrid), intent(out) :: out
!
      integer :: ibox, partner, nbox
      integer, parameter :: ltag=102, utag=103
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(nx,ny) :: recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      nbox = nx*ny
!
      do ibox = 0, nprocy-1
        partner = find_proc(ipx,ibox,ipz)
        if (iproc == partner) then
          ! data is local
          out(:,ny*ibox+1:ny*(ibox+1)) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (in, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(:,ny*ibox+1:ny*(ibox+1)) = recv_buf
        endif
      enddo
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_y_2D
!***********************************************************************
    subroutine remap_to_pencil_y_3D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 3D arrays in y only for nprocy>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
!
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ltag=102, utag=103
      integer :: inx, inz ! size of the first and third dimension
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(:,:,:), allocatable :: recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      inx = size (in, 1)
      inz = size (in, 3)
      nbox = inx*ny*inz
!
      if (inx /= size (out, 1)) &
          call stop_fatal ('remap_to_pencil_y_3D: first dimension differs for input and output', lfirst_proc_y)
      if (inz /= size (out, 3)) &
          call stop_fatal ('remap_to_pencil_y_3D: third dimension differs for input and output', lfirst_proc_y)
!
      if (size (in, 2) /= ny) &
          call stop_fatal ('remap_to_pencil_y_3D: second dimension of input must be ny', lfirst_proc_y)
      if (size (out, 2) /= nygrid) &
          call stop_fatal ('remap_to_pencil_y_3D: second dimension of output must be nygrid', lfirst_proc_y)
!
      allocate (recv_buf(inx,ny,inz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_y_3D: Could not allocate memory for recv_buf', .true.)
!
      do ibox = 0, nprocy-1
        partner = find_proc(ipx,ibox,ipz)
        if (iproc == partner) then
          ! data is local
          out(:,ny*ibox+1:ny*(ibox+1),:) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (in, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(:,ny*ibox+1:ny*(ibox+1),:) = recv_buf
        endif
      enddo
!
      if (allocated (recv_buf)) deallocate (recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_y_3D
!***********************************************************************
    subroutine remap_to_pencil_y_4D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 4D arrays in y only for nprocy>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
!
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ltag=102, utag=103
      integer :: inx, inz, ina ! size of the first, third, and fourth dimension
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(:,:,:,:), allocatable :: recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      inx = size (in, 1)
      inz = size (in, 3)
      ina = size (in, 4)
      nbox = inx*ny*inz*ina
!
      if (inx /= size (out, 1)) &
          call stop_fatal ('remap_to_pencil_y_4D: first dimension differs for input and output', lfirst_proc_y)
      if (inz /= size (out, 3)) &
          call stop_fatal ('remap_to_pencil_y_4D: third dimension differs for input and output', lfirst_proc_y)
      if (ina /= size (out, 4)) &
          call stop_fatal ('remap_to_pencil_y_4D: fourth dimension differs for input and output', lfirst_proc_y)
!
      if (size (in, 2) /= ny) &
          call stop_fatal ('remap_to_pencil_y_4D: second dimension of input must be ny', lfirst_proc_y)
      if (size (out, 2) /= nygrid) &
          call stop_fatal ('remap_to_pencil_y_4D: second dimension of output must be nygrid', lfirst_proc_y)
!
      allocate (recv_buf(inx,ny,inz,ina), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_y_4D: Could not allocate memory for recv_buf', .true.)
!
      do ibox = 0, nprocy-1
        partner = find_proc(ipx,ibox,ipz)
        if (iproc == partner) then
          ! data is local
          out(:,ny*ibox+1:ny*(ibox+1),:,:) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (in, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(:,ny*ibox+1:ny*(ibox+1),:,:) = recv_buf
        endif
      enddo
!
      deallocate (recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_y_4D
!***********************************************************************
    subroutine unmap_from_pencil_y_1D(in, out)
!
!  Unmaps pencil shaped 1D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocy>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(nygrid), intent(in) :: in
      real, dimension(ny), intent(out) :: out
!
      !$omp workshare
      out = in(ny*ipy+1:ny*(ipy+1))
      !$omp end workshare
!
    endsubroutine unmap_from_pencil_y_1D
!***********************************************************************
    subroutine unmap_from_pencil_y_2D(in, out)
!
!  Unmaps pencil shaped 2D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocy>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(nx,nygrid), intent(in) :: in
      real, dimension(nx,ny), intent(out) :: out
!
      !$omp workshare
      out = in(:,ny*ipy+1:ny*(ipy+1))
      !$omp end workshare
!
    endsubroutine unmap_from_pencil_y_2D
!***********************************************************************
    subroutine unmap_from_pencil_y_3D(in, out)
!
!  Unmaps pencil shaped 3D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocy>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
!
      !$omp workshare
      out = in(:,ny*ipy+1:ny*(ipy+1),:)
      !$omp end workshare
!
    endsubroutine unmap_from_pencil_y_3D
!***********************************************************************
    subroutine unmap_from_pencil_y_4D(in, out)
!
!  Unmaps pencil shaped 4D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocy>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
!
      !$omp workshare
      out = in(:,ny*ipy+1:ny*(ipy+1),:,:)
      !$omp end workshare
!
    endsubroutine unmap_from_pencil_y_4D
!***********************************************************************
    subroutine remap_to_pencil_z_1D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 1D arrays in z only for nprocz>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(nz), intent(in) :: in
      real, dimension(nzgrid), intent(out) :: out
!
      integer :: ibox, partner
      integer, parameter :: ltag=102, utag=103
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(nz) :: recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      do ibox = 0, nprocz-1
        partner = find_proc(ipx,ipy,ibox)
        if (iproc == partner) then
          ! data is local
          out(nz*ibox+1:nz*(ibox+1)) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nz, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nz, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nz, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (in, nz, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(nz*ibox+1:nz*(ibox+1)) = recv_buf
        endif
      enddo
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_z_1D
!***********************************************************************
    subroutine remap_to_pencil_z_2D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 2D arrays in z only for nprocz>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
!
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ltag=102, utag=103
      integer :: ina ! size of the second dimension
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(:,:), allocatable :: recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      ina = size (in, 2)
      nbox = nz*ina
!
      if (size (in, 1) /= nz) &
          call stop_fatal ('remap_to_pencil_z_2D: first dimension of input must be nz', lfirst_proc_z)
      if (size (out, 2) /= nzgrid) &
          call stop_fatal ('remap_to_pencil_z_2D: first dimension of output must be nzgrid', lfirst_proc_y)
      if (ina /= size (out, 2)) &
          call stop_fatal ('remap_to_pencil_z_2D: second dimension differs for input and output', lfirst_proc_z)
!
      ! Allocate memory for large arrays.
      allocate (recv_buf(nz,ina), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_z_2D: Could not allocate memory for recv_buf', .true.)
!
      do ibox = 0, nprocz-1
        partner = find_proc(ipx,ipy,ibox)
        if (iproc == partner) then
          ! data is local
          out(nz*ibox+1:nz*(ibox+1),:) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (in, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(nz*ibox+1:nz*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_z_2D
!***********************************************************************
    subroutine remap_to_pencil_z_3D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 3D arrays in z only for nprocz>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
!
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ltag=102, utag=103
      integer :: inx, iny ! size of the first and third dimension
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(:,:,:), allocatable :: recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      inx = size (in, 1)
      iny = size (in, 2)
      nbox = inx*iny*nz
!
      if (inx /= size (out, 1)) &
          call stop_fatal ('remap_to_pencil_z_3D: first dimension differs for input and output', lfirst_proc_y)
      if (iny /= size (out, 2)) &
          call stop_fatal ('remap_to_pencil_z_3D: second dimension differs for input and output', lfirst_proc_y)
!
      if (size (in, 3) /= nz) &
          call stop_fatal ('remap_to_pencil_z_3D: third dimension of input must be nz', lfirst_proc_y)
      if (size (out, 3) /= nzgrid) &
          call stop_fatal ('remap_to_pencil_z_3D: third dimension of output must be nzgrid', lfirst_proc_y)
!
      allocate (recv_buf(inx,iny,nz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_z_3D: Could not allocate memory for recv_buf', .true.)
!
      do ibox = 0, nprocz-1
        partner = find_proc(ipx,ipy,ibox)
        if (iproc == partner) then
          ! data is local
          out(:,:,nz*ibox+1:nz*(ibox+1)) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (in, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(:,:,nz*ibox+1:nz*(ibox+1)) = recv_buf
        endif
      enddo
!
      deallocate (recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_z_3D
!***********************************************************************
    subroutine remap_to_pencil_z_4D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 4D arrays in z only for nprocz>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
!
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ltag=102, utag=103
      integer :: inx, iny, ina ! size of the first, second, and fourth dimension
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(:,:,:,:), allocatable :: recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      inx = size (in, 1)
      iny = size (in, 2)
      ina = size (in, 4)
      nbox = inx*iny*nz*ina
!
      if (inx /= size (out, 1)) &
          call stop_fatal ('remap_to_pencil_z_4D: first dimension differs for input and output', lfirst_proc_y)
      if (iny /= size (out, 2)) &
          call stop_fatal ('remap_to_pencil_z_4D: second dimension differs for input and output', lfirst_proc_y)
      if (ina /= size (out, 4)) &
          call stop_fatal ('remap_to_pencil_z_4D: fourth dimension differs for input and output', lfirst_proc_y)
!
      if (size (in, 3) /= nz) &
          call stop_fatal ('remap_to_pencil_z_4D: third dimension of input must be nz', lfirst_proc_y)
      if (size (out, 3) /= nzgrid) &
          call stop_fatal ('remap_to_pencil_z_4D: third dimension of output must be nzgrid', lfirst_proc_y)
!
      allocate (recv_buf(inx,iny,nz,ina), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_z_4D: Could not allocate memory for recv_buf', .true.)
!
      do ibox = 0, nprocz-1
        partner = find_proc(ipx,ipy,ibox)
        if (iproc == partner) then
          ! data is local
          out(:,:,nz*ibox+1:nz*(ibox+1),:) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (in, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(:,:,nz*ibox+1:nz*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_z_4D
!***********************************************************************
    subroutine unmap_from_pencil_z_1D(in, out)
!
!  Unmaps pencil shaped 1D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocz>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(nzgrid), intent(in) :: in
      real, dimension(nz), intent(out) :: out
!
      !$omp workshare
      out = in(nz*ipz+1:nz*(ipz+1))
      !$omp end workshare
!
    endsubroutine unmap_from_pencil_z_1D
!***********************************************************************
    subroutine unmap_from_pencil_z_2D(in, out)
!
!  Unmaps pencil shaped 2D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocz>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
!
      !$omp workshare
      out = in(nz*ipz+1:nz*(ipz+1),:)
      !$omp end workshare
!
    endsubroutine unmap_from_pencil_z_2D
!***********************************************************************
    subroutine unmap_from_pencil_z_3D(in, out)
!
!  Unmaps pencil shaped 3D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocz>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
!
      !$omp workshare
      out = in(:,:,nz*ipz+1:nz*(ipz+1))
      !$omp end workshare
!
    endsubroutine unmap_from_pencil_z_3D
!***********************************************************************
    subroutine unmap_from_pencil_z_4D(in, out)
!
!  Unmaps pencil shaped 4D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocz>1.
!
!  13-dec-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
!
      !$omp workshare
      out = in(:,:,nz*ipz+1:nz*(ipz+1),:)
      !$omp end workshare
!
    endsubroutine unmap_from_pencil_z_4D
!***********************************************************************
    subroutine remap_to_pencil_xy_2D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 2D arrays in x and y only for nprocx>1.
!
!   04-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
!
      integer, parameter :: inx=nx, iny=ny
      integer, parameter :: onx=nxgrid, ony=ny/nprocx
      integer, parameter :: bnx=nx, bny=ny/nprocx ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ltag=104, utag=105
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (nprocx == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
!
      nbox = bnx*bny
!
      if (mod (ny, nprocx) /= 0) &
          call stop_fatal ('remap_to_pencil_xy_2D: ny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      if ((size (in, 1) /= inx) .or. ((size (in, 2) /= iny))) &
          call stop_fatal ('remap_to_pencil_xy_2D: input array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (out, 1) /= onx) .or. ((size (out, 2) /= ony))) &
          call stop_fatal ('remap_to_pencil_xy_2D: output array size mismatch /= nxgrid,ny/nprocx', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_xy_2D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_xy_2D: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocx-1
        partner = find_proc(ibox,ipy,ipz)
        if (iproc == partner) then
          ! data is local
          out(bnx*ibox+1:bnx*(ibox+1),:) = in(:,bny*ibox+1:bny*(ibox+1))
        else
          ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1))
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
      !$omp end single

      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_xy_2D
!***********************************************************************
    subroutine remap_to_pencil_xy_2D_other(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 2D arrays in x and y only for nprocx>1.
!
!   04-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
!
      integer :: nnx,nny,inx, iny, onx, ony, bnx, bny
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ltag=104, utag=105
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (nprocx == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif

      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      nnx=size(in,1) ; nny=size(in,2)
      inx=nnx        ; iny=nny
      onx=nprocx*nnx ; ony=nny/nprocx
      bnx=nnx        ; bny=nny/nprocx ! transfer box sizes
!
!
      nbox = bnx*bny
!
      if (mod (nny, nprocx) /= 0) &
          call stop_fatal ('remap_to_pencil_xy_2D_other: nny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      if ((size (in, 1) /= inx) .or. ((size (in, 2) /= iny))) &
          call stop_fatal ('remap_to_pencil_xy_2D_other: input array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (out, 1) /= onx) .or. ((size (out, 2) /= ony))) &
          call stop_fatal ('remap_to_pencil_xy_2D_other: output array size mismatch /= nxgrid,ny/nprocx', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_xy_2D_other: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_xy_2D_other: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocx-1
        partner = find_proc(ibox,ipy,ipz)
        if (iproc == partner) then
          ! data is local
          out(bnx*ibox+1:bnx*(ibox+1),:) = in(:,bny*ibox+1:bny*(ibox+1))
        else
          ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1))
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_xy_2D_other
!***********************************************************************
    subroutine remap_to_pencil_xy_3D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 3D arrays in x and y only for nprocx>1.
!
!  14-jul-2010/Bourdin.KIS: coded
!  13-sep-2014/ccyang: revamped to accommodate either with or without ghost cells.
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
!
      integer, parameter :: bnx = nx, bny = ny / nprocx
      integer, parameter :: ltag = 104, utag = 105
      real, dimension(:,:,:), allocatable :: send_buf, recv_buf
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: inx, iny, inz
      integer :: onx, ony, onz
      integer :: ibox, partner, nbox, alloc_err
      integer :: ngc
      integer, optional :: comm
      logical, optional :: lsync
!
!  No need to remap if nprocx = 1.
!
      if (nprocx == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
!
!  Check the dimensions.
!
      if (mod(ny, nprocx) /= 0) &
        call stop_fatal('remap_to_pencil_xy_3D: ny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      inx = size(in, 1)
      iny = size(in, 2)
      dim: if (inx == nx .and. iny == ny) then
        onx = nxgrid
        ony = ny / nprocx
        ngc = 0
      elseif (inx == mx .and. iny == my) then dim
        onx = mxgrid
        ony = ny / nprocx + 2 * nghost
        ngc = nghost
      else dim
        call stop_fatal('remap_to_pencil_xy_3D: input array size mismatch', lfirst_proc_xy)
      endif dim
!
      inz = size(in, 3)
      onz = size(out, 3)
      nbox = (bnx + 2 * ngc) * (bny + 2 * ngc) * onz
!
      if (size(out,1) /= onx .or. size(out,2) /= ony) &
        call stop_fatal('remap_to_pencil_xy_3D: output array size mismatch', lfirst_proc_xy)
!
      if (inz /= onz) call stop_fatal('remap_to_pencil_xy_3D: sizes differ in the z direction', lfirst_proc_xy)
!
!  Allocate working arrays.
!
      allocate(send_buf(bnx+2*ngc,bny+2*ngc,inz), recv_buf(bnx+2*ngc,bny+2*ngc,inz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal('remap_to_pencil_xy_3D: allocation failed. ', .true.)
!
!  Communicate.
!
      box: do ibox = 0, nprocx - 1
        partner = find_proc(ibox,ipy,ipz)
        local: if (iproc == partner) then  ! data is local
          recv_buf = in(:,bny*ibox+1:bny*(ibox+1)+2*ngc,:)
        else local                         ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1)+2*ngc,:)
          commun: if (iproc > partner) then  ! above diagonal: send first, receive then
            call MPI_SEND(send_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV(recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else commun                        ! below diagonal: receive first, send then
            call MPI_RECV(recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND(send_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif commun
        endif local
        out(ngc+bnx*ibox+1:ngc+bnx*(ibox+1),:,:) = recv_buf(ngc+1:ngc+bnx,:,:)
        ghost: if (ngc > 0) then
          if (ibox == 0) out(1:ngc,:,:) = recv_buf(1:ngc,:,:)
          if (ibox == nprocx - 1) out(onx-ngc+1:onx,:,:) = recv_buf(ngc+bnx+1:bnx+2*ngc,:,:)
        endif ghost
      enddo box
!
!  Deallocate working arrays.
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_xy_3D
!***********************************************************************
    subroutine remap_to_pencil_xy_4D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 4D arrays in x and y only for nprocx>1.
!
!  14-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
!
      integer, parameter :: inx=nx, iny=ny
      integer, parameter :: onx=nxgrid, ony=ny/nprocx
      integer :: inz, ina, onz, ona ! sizes of in and out arrays
      integer, parameter :: bnx=nx, bny=ny/nprocx ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ltag=104, utag=105
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (nprocx == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
!
      inz = size (in, 3)
      ina = size (in, 4)
      onz = size (out, 3)
      ona = size (out, 4)
      nbox = bnx*bny*onz*ona
!
      if (mod (ny, nprocx) /= 0) &
          call stop_fatal ('remap_to_pencil_xy_4D: ny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      if ((size (in, 1) /= inx) .or. ((size (in, 2) /= iny))) &
          call stop_fatal ('remap_to_pencil_xy_4D: input array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (out, 1) /= onx) .or. ((size (out, 2) /= ony))) &
          call stop_fatal ('remap_to_pencil_xy_4D: output array size mismatch /= nxgrid,ny/nprocx', lfirst_proc_xy)
      if (inz /= onz) &
          call stop_fatal ('remap_to_pencil_xy_4D: inz/=onz - sizes differ in the z direction', lfirst_proc_xy)
      if (ina /= ona) &
          call stop_fatal ('remap_to_pencil_xy_4D: ina/=ona - sizes differ in the 4th dimension', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny,onz,ona), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_xy_4D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny,onz,ona), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_xy_4D: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocx-1
        partner = find_proc(ibox,ipy,ipz)
        if (iproc == partner) then
          ! data is local
          out(bnx*ibox+1:bnx*(ibox+1),:,:,:) = in(:,bny*ibox+1:bny*(ibox+1),:,:)
        else
          ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1),:,:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:,:,:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_xy_4D
!***********************************************************************
    subroutine unmap_from_pencil_xy_2D(in, out,comm,lsync)
!
!  Unmaps pencil shaped 2D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocx>1.
!
!   4-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
!
      integer, parameter :: inx=nxgrid, iny=ny/nprocx
      integer, parameter :: onx=nx, ony=ny
      integer, parameter :: bnx=nx, bny=ny/nprocx ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ltag=106, utag=107
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (nprocx == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
!
      nbox = bnx*bny
!
      if (mod (ny, nprocx) /= 0) &
          call stop_fatal ('unmap_from_pencil_xy_2D: ny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      if ((size (in, 1) /= inx) .or. ((size (in, 2) /= iny))) &
          call stop_fatal ('unmap_from_pencil_xy_2D: input array size mismatch /= nxgrid,ny/nprocx', lfirst_proc_xy)
      if ((size (out, 1) /= onx) .or. ((size (out, 2) /= ony))) &
          call stop_fatal ('unmap_from_pencil_xy_2D: output array size mismatch /= nx,ny', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_xy_2D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_xy_2D: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocx-1
        partner = find_proc(ibox,ipy,ipz)
        if (iproc == partner) then
          ! data is local
          out(:,bny*ibox+1:bny*(ibox+1)) = in(bnx*ibox+1:bnx*(ibox+1),:)
        else
          ! communicate with partner
          send_buf = in(bnx*ibox+1:bnx*(ibox+1),:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(:,bny*ibox+1:bny*(ibox+1)) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)

      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine unmap_from_pencil_xy_2D
!***********************************************************************
    subroutine unmap_from_pencil_xy_2D_other(in, out,comm,lsync)
!
!  Unmaps pencil shaped 2D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocx>1.
!
!   4-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
!
      integer :: nnx,nny,inx, iny, onx, ony, bnx, bny, nxgrid_other, nygrid_other
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ltag=106, utag=107
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync

      if (nprocx == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      if (nxgrid/=nygrid) &
           call stop_fatal("unmap_from_pencil_xy_2D_other: this subroutine works only for nxgrid==nygrid",lfirst_proc_xy)
!
      nxgrid_other=size(in,1) ; nygrid_other=nxgrid_other
      nnx=nxgrid_other/nprocx ; nny=nygrid_other/nprocy
      inx=nxgrid_other        ; iny=nny/nprocx
      onx=nnx                 ; ony=nny
      bnx=nnx                 ; bny=nny/nprocx ! transfer box sizes
!
      nbox = bnx*bny
!
      if (mod (nny, nprocx) /= 0) &
          call stop_fatal ('unmap_from_pencil_xy_2D_other: nny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      if ((size (in, 1) /= inx) .or. ((size (in, 2) /= iny))) &
          call stop_fatal ('unmap_from_pencil_xy_2D_other: input array size mismatch /= nxgrid_other,nny/nprocx', lfirst_proc_xy)
      if ((size (out, 1) /= onx) .or. ((size (out, 2) /= ony))) &
          call stop_fatal ('unmap_from_pencil_xy_2D_other: output array size mismatch /= nnx,nny', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_xy_2D_other: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_xy_2D_other: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocx-1
        partner = find_proc(ibox,ipy,ipz)
        if (iproc == partner) then
          ! data is local
          out(:,bny*ibox+1:bny*(ibox+1)) = in(bnx*ibox+1:bnx*(ibox+1),:)
        else
          ! communicate with partner
          send_buf = in(bnx*ibox+1:bnx*(ibox+1),:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(:,bny*ibox+1:bny*(ibox+1)) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine unmap_from_pencil_xy_2D_other
!***********************************************************************
    subroutine unmap_from_pencil_xy_3D(in, out,comm,lsync)
!
!  Unmaps pencil shaped 3D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocx>1.
!
!  14-jul-2010/Bourdin.KIS: coded
!  14-sep-2014/ccyang: revamped to accommodate either with or without ghost cells.
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
!
      integer, parameter :: nypx = ny / nprocx
      integer, parameter :: bnx = nx, bny = nypx
      integer, parameter :: ltag = 106, utag = 107
      real, dimension(:,:,:), allocatable :: send_buf, recv_buf
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: inx, iny, inz
      integer :: onx, ony, onz
      integer :: ibox, partner, nbox, alloc_err
      integer :: ngc
      integer, optional :: comm
      logical, optional :: lsync
!
!  No need to unmap if nprocx = 1.
!
      if (nprocx == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
!
!  Check the dimensions.
!
      if (mod(ny, nprocx) /= 0) &
        call stop_fatal('unmap_from_pencil_xy_3D: ny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      inx = size(in, 1)
      iny = size(in, 2)
      dim: if (inx == nxgrid .and. iny == nypx) then
        onx = nx
        ony = ny
        ngc = 0
      elseif (inx == mxgrid .and. iny == nypx + 2 * nghost) then dim
        onx = mx
        ony = my
        ngc = nghost
      else dim
        call stop_fatal('unmap_from_pencil_xy_3D: input array size mismatch', lfirst_proc_xy)
      endif dim
!
      if (size(out, 1) /= onx .or. size(out, 2) /= ony) &
        call stop_fatal('unmap_from_pencil_xy_3D: output array size mismatch', lfirst_proc_xy)
!
      inz = size(in, 3)
      onz = size(out, 3)
      if (inz /= onz) call stop_fatal('unmap_from_pencil_xy_3D: sizes differ in the z direction', lfirst_proc_xy)
!
      nbox = (bnx + 2 * ngc) * (bny + 2 * ngc) * onz
!
!  Allocate working arrays.
!
      allocate(send_buf(bnx+2*ngc,bny+2*ngc,inz), recv_buf(bnx+2*ngc,bny+2*ngc,inz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal('unmap_from_pencil_xy_3D: allocation failed. ', .true.)
!
!  Communicate.
!
      box: do ibox = 0, nprocx - 1
        partner = find_proc(ibox,ipy,ipz)
        local: if (iproc == partner) then  ! data is local
          recv_buf = in(bnx*ibox+1:bnx*(ibox+1)+2*ngc,:,:)
        else local                         ! communicate with partner
          send_buf = in(bnx*ibox+1:bnx*(ibox+1)+2*ngc,:,:)
          commun: if (iproc > partner) then  ! above diagonal: send first, receive then
            call MPI_SEND(send_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV(recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else commun                        ! below diagonal: receive first, send then
            call MPI_RECV(recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND(send_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif commun
        endif local
        out(:,ngc+bny*ibox+1:ngc+bny*(ibox+1),:) = recv_buf(:,ngc+1:ngc+bny,:)
        ghost: if (ngc > 0) then
          if (ibox == 0) out(:,1:ngc,:) = recv_buf(:,1:ngc,:)
          if (ibox == nprocx - 1) out(:,ony-ngc+1:ony,:) = recv_buf(:,ngc+bny+1:bny+2*ngc,:)
        endif ghost
      enddo box
!
!  Deallocate working arrays.
!
      deallocate(send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine unmap_from_pencil_xy_3D
!***********************************************************************
    subroutine unmap_from_pencil_xy_4D(in, out,comm,lsync)
!
!  Unmaps pencil shaped 4D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocx>1.
!
!  14-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
!
      integer, parameter :: inx=nxgrid, iny=ny/nprocx
      integer, parameter :: onx=nx, ony=ny
      integer :: inz, ina, onz, ona ! sizes of in and out arrays
      integer, parameter :: bnx=nx, bny=ny/nprocx ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ltag=106, utag=107
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (nprocx == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
!
      inz = size (in, 3)
      ina = size (in, 4)
      onz = size (out, 3)
      ona = size (out, 4)
      nbox = bnx*bny*onz*ona
!
      if (mod (ny, nprocx) /= 0) &
          call stop_fatal ('unmap_from_pencil_xy_4D: ny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      if ((size (in, 1) /= inx) .or. ((size (in, 2) /= iny))) &
          call stop_fatal ('unmap_from_pencil_xy_4D: input array size mismatch /= nxgrid,ny/nprocx', lfirst_proc_xy)
      if ((size (out, 1) /= onx) .or. ((size (out, 2) /= ony))) &
          call stop_fatal ('unmap_from_pencil_xy_4D: output array size mismatch /= nx,ny', lfirst_proc_xy)
      if (inz /= onz) &
          call stop_fatal ('unmap_from_pencil_xy_4D: inz/=onz - sizes differ in the z direction', lfirst_proc_xy)
      if (ina /= ona) &
          call stop_fatal ('unmap_from_pencil_xy_4D: ina/=ona - sizes differ in the 4th dimension', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny,onz,ona), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_xy_4D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny,onz,ona), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_xy_4D: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocx-1
        partner = find_proc(ibox,ipy,ipz)
        if (iproc == partner) then
          ! data is local
          out(:,bny*ibox+1:bny*(ibox+1),:,:) = in(bnx*ibox+1:bnx*(ibox+1),:,:,:)
        else
          ! communicate with partner
          send_buf = in(bnx*ibox+1:bnx*(ibox+1),:,:,:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(:,bny*ibox+1:bny*(ibox+1),:,:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine unmap_from_pencil_xy_4D
!***********************************************************************
    subroutine transp_pencil_xy_2D(in, out,comm,lsync)
!
!  Transpose 2D data distributed on several processors.
!  This routine transposes arrays in x and y only.
!  The data must be mapped in pencil shape, especially for nprocx>1.
!
!   4-jul-2010/Bourdin.KIS: coded, adapted parts of transp_xy
!  21-jun-2013/Bourdin.KIS: reworked, parellized MPI communication
!
      use General, only: count_bits
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
!
      integer :: inx, iny, onx, ony ! sizes of in and out arrays
      integer :: bnx, bny, nbox ! destination box sizes and number of elements
      integer :: num_it ! number of iterations for box-data transfer
      integer :: it, ibox, partner, alloc_err
      integer, parameter :: ltag=108, utag=109
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      inx = size (in, 1)
      iny = size (in, 2)
      onx = size (out, 1)
      ony = size (out, 2)
      bnx = onx/nprocxy
      bny = ony
      nbox = bnx*bny
!
      if (mod (onx, nprocxy) /= 0) &
          call stop_fatal ('transp_pencil_xy_2D: onx needs to be an integer multiple of nprocxy', lfirst_proc_xy)
!
      if ((inx /= bny*nprocxy) .or. (iny /= bnx)) &
          call stop_fatal ('transp_pencil_xy_2D: input array has unmatching size', lfirst_proc_xy)
      if ((onx /= bnx*nprocxy) .or. (ony /= bny)) &
          call stop_fatal ('transp_pencil_xy_2D: output array has unmatching size', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('transp_pencil_xy_2D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('transp_pencil_xy_2D: not enough memory for recv_buf!', .true.)
!
      num_it = 2**count_bits (nprocxy-1)
      do it = 0, num_it-1
        ibox = IEOR ((iproc - ipz*nprocxy), it)
        partner = ipz*nprocxy + ibox
        if (ibox >= nprocxy) cycle
        if (iproc == partner) then
          ! data is local
          out(bnx*ibox+1:bnx*(ibox+1),:) = transpose (in(bny*ibox+1:bny*(ibox+1),:))
        else
          ! communicate with partner
          send_buf = transpose (in(bny*ibox+1:bny*(ibox+1),:))
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine transp_pencil_xy_2D
!***********************************************************************
    subroutine transp_pencil_xy_3D(in, out, lghost,comm,lsync)
!
!  Transpose 3D data distributed on several processors.
!  This routine transposes arrays in x and y only.
!  The data must be mapped in pencil shape, especially for nprocx>1.
!
!  14-jul-2010/Bourdin.KIS: coded, adapted parts of transp_xy
!  21-jun-2013/Bourdin.KIS: reworked, parellized MPI communication
!  15-sep-2014/ccyang: revamped to accommodate arrays either with or without ghost cells
!
      use General, only: count_bits
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
      logical, intent(in), optional :: lghost
!
      integer, parameter :: ltag = 108, utag = 109
      real, dimension(:,:,:), allocatable :: send_buf, recv_buf
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: inx, iny, inz, onx, ony, onz ! sizes of in and out arrays
      integer :: bnx, bny, nbox ! destination box sizes and number of elements
      integer :: ibox, partner, alloc_err, iz, ngc
      integer, optional :: comm
      logical, optional :: lsync
!
!  Check if ghost cells are included.
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      ngc = 0
      if (present(lghost)) then
        if (lghost) ngc = nghost
      endif
!
!  Check the dimensions.
!
      inx = size(in, 1)
      iny = size(in, 2)
      inz = size(in, 3)
      onx = size(out, 1)
      ony = size(out, 2)
      onz = size(out, 3)
!
      bnx = (inx - 2 * ngc) / nprocxy
      bny = iny - 2 * ngc
!
      if ((inx /= bnx * nprocxy + 2 * ngc) .or. (iny /= bny + 2 * ngc)) &
        call stop_fatal('transp_pencil_xy_3D: input array has unmatching shape', lfirst_proc_xy)
      if ((onx /= bny * nprocxy + 2 * ngc) .or. (ony /= bnx + 2 * ngc)) &
        call stop_fatal('transp_pencil_xy_3D: output array has unmatching shape', lfirst_proc_xy)
      if (inz /= onz) call stop_fatal('transp_pencil_xy_3D: sizes differ in the z direction', lfirst_proc_xy)
!
      nbox = (bnx + 2 * ngc) * (bny + 2 * ngc) * onz
!
!  Allocate working arrays.
!
      allocate (send_buf(bnx+2*ngc,bny+2*ngc,onz), recv_buf(bnx+2*ngc,bny+2*ngc,onz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal('transp_pencil_xy_3D: allocation failed. ', .true.)
!
!  Communicate.
!
      box: do ibox = 0, nprocxy - 1
        partner = find_proc(modulo(ibox,nprocx),ibox/nprocx,ipz)
        local: if (iproc == partner) then  ! data is local
          recv_buf = in(bnx*ibox+1:bnx*(ibox+1)+2*ngc,:,:)
        else local                        ! communicate with partner
          send_buf = in(bnx*ibox+1:bnx*(ibox+1)+2*ngc,:,:)
          if (iproc > partner) then  ! above diagonal: send first, receive then
            call MPI_SEND(send_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV(recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else
            call MPI_RECV(recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND(send_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
        endif local
        do iz = 1, inz
          out(ngc+bny*ibox+1:ngc+bny*(ibox+1),:,iz) = transpose(recv_buf(:,ngc+1:ngc+bny,iz))
          if (ngc > 0) then
            if (ibox == 0) out(1:ngc,:,iz) = transpose(recv_buf(:,1:ngc,iz))
            if (ibox == nprocxy - 1) out(onx-ngc+1:onx,:,iz) = transpose(recv_buf(:,ngc+bny+1:bny+2*ngc,iz))
          endif
        enddo
      enddo box
!
!  Deallocate working arrays.
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine transp_pencil_xy_3D
!***********************************************************************
    subroutine transp_pencil_xy_4D(in, out,comm,lsync)
!
!  Transpose 4D data distributed on several processors.
!  This routine transposes arrays in x and y only.
!  The data must be mapped in pencil shape, especially for nprocx>1.
!
!  14-jul-2010/Bourdin.KIS: coded, adapted parts of transp_xy
!  21-jun-2013/Bourdin.KIS: reworked, parellized MPI communication
!
      use General, only: count_bits
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
!
      integer :: inx, iny, inz, ina, onx, ony, onz, ona ! sizes of in and out arrays
      integer :: bnx, bny, nbox ! destination box sizes and number of elements
      integer :: num_it ! number of iterations for box-data transfer
      integer :: it, ibox, partner, alloc_err, pos_z, pos_a
      integer, parameter :: ltag=108, utag=109
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
      inx = size (in, 1)
      iny = size (in, 2)
      inz = size (in, 3)
      ina = size (in, 4)
      onx = size (out, 1)
      ony = size (out, 2)
      onz = size (out, 3)
      ona = size (out, 4)
      bnx = onx/nprocxy
      bny = ony
      nbox = bnx*bny*onz*ona
!
      if (mod (onx, nprocxy) /= 0) &
          call stop_fatal ('transp_pencil_xy_4D: onx needs to be an integer multiple of nprocxy', lfirst_proc_xy)
!
      if ((inx /= bny*nprocxy) .or. (iny /= bnx)) &
          call stop_fatal ('transp_pencil_xy_4D: input array has unmatching size', lfirst_proc_xy)
      if ((onx /= bnx*nprocxy) .or. (ony /= bny)) &
          call stop_fatal ('transp_pencil_xy_4D: output array has unmatching size', lfirst_proc_xy)
      if (inz /= onz) &
          call stop_fatal ('transp_pencil_xy_4D: inz/=onz - sizes differ in the z direction', lfirst_proc_xy)
      if (ina /= ona) &
         call stop_fatal ('transp_pencil_xy_4D: ina/=ona - sizes differ in the 4th dimension', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny,onz,ona), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('transp_pencil_xy_4D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny,onz,ona), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('transp_pencil_xy_4D: not enough memory for recv_buf!', .true.)
!
      num_it = 2**count_bits (nprocxy-1)
      do it = 0, num_it-1
        ibox = IEOR ((iproc - ipz*nprocxy), it)
        partner = ipz*nprocxy + ibox
        if (ibox >= nprocxy) cycle
        if (iproc == partner) then
          ! data is local
          do pos_z = 1, onz
            do pos_a = 1, ona
              out(bnx*ibox+1:bnx*(ibox+1),:,pos_z,pos_a) = transpose(in(bny*ibox+1:bny*(ibox+1),:,pos_z,pos_a))
            enddo
          enddo
        else
          ! communicate with partner
          do pos_z = 1, onz
            do pos_a = 1, ona
              send_buf(:,:,pos_z,pos_a) = transpose(in(bny*ibox+1:bny*(ibox+1),:,pos_z,pos_a))
            enddo
          enddo
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, utag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ltag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:,:,:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine transp_pencil_xy_4D
!***********************************************************************
    subroutine remap_to_pencil_yz_3D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into z-pencil shape.
!  This routine remaps 3D arrays in y and z only for nprocz>1.
!
!  27-oct-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
!
      integer, parameter :: iny=ny, inz=nz
      integer, parameter :: ony=ny/nprocz, onz=nzgrid
      integer, parameter :: bny=ny/nprocz, bnz=nz ! transfer box sizes
      integer :: inx, onx ! sizes of in and out arrays
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ytag=110
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (nprocz == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
!
      inx = size (in, 1)
      onx = size (out, 1)
      nbox = onx*bny*bnz
!
      if (mod (ny, nprocz) /= 0) &
          call stop_fatal ('remap_to_pencil_yz_3D: ny needs to be an integer multiple of nprocz', lfirst_proc_yz)
!
      if ((size (in, 2) /= iny) .or. ((size (in, 3) /= inz))) &
          call stop_fatal ('remap_to_pencil_yz_3D: input array size mismatch /= ny,nz', lfirst_proc_yz)
      if ((size (out, 2) /= ony) .or. ((size (out, 3) /= onz))) &
          call stop_fatal ('remap_to_pencil_yz_3D: output array size mismatch /= ny/nprocz,nzgrid', lfirst_proc_yz)
      if (inx /= onx) &
          call stop_fatal ('remap_to_pencil_yz_3D: inx/=onx - sizes differ in the x direction', lfirst_proc_yz)
!
      allocate (send_buf(onx,bny,bnz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_yz_3D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(onx,bny,bnz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_yz_3D: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocz-1
        partner = find_proc(ipx,ipy,ibox)
        if (iproc == partner) then
          ! data is local
          out(:,:,bnz*ibox+1:bnz*(ibox+1)) = in(:,bny*ibox+1:bny*(ibox+1),:)
        else
          ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1),:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
!
          out(:,:,bnz*ibox+1:bnz*(ibox+1)) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_yz_3D
!***********************************************************************
    subroutine remap_to_pencil_yz_4D(in, out,comm,lsync)
!
!  Remaps data distributed on several processors into z-pencil shape.
!  This routine remaps 4D arrays in y and z only for nprocz>1.
!
!  27-oct-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
!
      integer, parameter :: iny=ny, inz=nz
      integer, parameter :: ony=ny/nprocz, onz=nzgrid
      integer, parameter :: bny=ny/nprocz, bnz=nz ! transfer box sizes
      integer :: inx, ina, onx, ona ! sizes of in and out arrays
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ytag=110
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (nprocz == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
!
      inx = size (in, 1)
      ina = size (in, 4)
      onx = size (out, 1)
      ona = size (out, 4)
      nbox = onx*bny*bnz*ona
!
      if (mod (ny, nprocz) /= 0) &
          call stop_fatal ('remap_to_pencil_yz_4D: ny needs to be an integer multiple of nprocz', lfirst_proc_yz)
!
      if ((size (in, 2) /= iny) .or. ((size (in, 3) /= inz))) &
          call stop_fatal ('remap_to_pencil_yz_4D: input array size mismatch /= ny,nz', lfirst_proc_yz)
      if ((size (out, 2) /= ony) .or. ((size (out, 3) /= onz))) &
          call stop_fatal ('remap_to_pencil_yz_4D: output array size mismatch /= ny/nprocz,nzgrid', lfirst_proc_yz)
      if (inx /= onx) &
          call stop_fatal ('remap_to_pencil_yz_4D: inx/=onx - sizes differ in the x direction', lfirst_proc_yz)
      if (ina /= ona) &
          call stop_fatal ('remap_to_pencil_yz_4D: ina/=ona - sizes differ in the 4th dimension', lfirst_proc_yz)
!
      allocate (send_buf(onx,bny,bnz,ona), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_yz_4D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(onx,bny,bnz,ona), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_yz_4D: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocz-1
        partner = find_proc(ipx,ipy,ibox)
        if (iproc == partner) then
          ! data is local
          out(:,:,bnz*ibox+1:bnz*(ibox+1),:) = in(:,bny*ibox+1:bny*(ibox+1),:,:)
        else
          ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1),:,:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(:,:,bnz*ibox+1:bnz*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine remap_to_pencil_yz_4D
!***********************************************************************
    subroutine unmap_from_pencil_yz_3D(in, out,comm,lsync)
!
!  Unmaps z-pencil shaped 3D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocz>1.
!
!  27-oct-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
!
      integer, parameter :: iny=ny/nprocz, inz=nzgrid
      integer, parameter :: ony=ny, onz=nz
      integer :: inx, onx ! sizes of in and out arrays
      integer, parameter :: bny=ny/nprocz, bnz=nz ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ytag=111
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (nprocz == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
!
      inx = size (in, 1)
      onx = size (out, 1)
      nbox = onx*bny*bnz
!
      if (mod (ny, nprocz) /= 0) &
          call stop_fatal ('unmap_from_pencil_yz_3D: ny needs to be an integer multiple of nprocz', lfirst_proc_yz)
!
      if ((size (in, 2) /= iny) .or. ((size (in, 3) /= inz))) &
          call stop_fatal ('unmap_from_pencil_yz_3D: input array size mismatch /= ny/nprocz,nygrid', lfirst_proc_yz)
      if ((size (out, 2) /= ony) .or. ((size (out, 3) /= onz))) &
          call stop_fatal ('unmap_from_pencil_yz_3D: output array size mismatch /= ny,nz', lfirst_proc_yz)
      if (inx /= onx) &
          call stop_fatal ('unmap_from_pencil_yz_3D: inx/=onx - sizes differ in the x direction', lfirst_proc_yz)
!
      allocate (send_buf(onx,bny,bnz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_yz_3D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(onx,bny,bnz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_yz_3D: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocz-1
        partner = find_proc(ipx,ipy,ibox)
        if (iproc == partner) then
          ! data is local
          out(:,bny*ibox+1:bny*(ibox+1),:) = in(:,:,bnz*ibox+1:bnz*(ibox+1))
        else
          ! communicate with partner
          send_buf = in(:,:,bnz*ibox+1:bnz*(ibox+1))
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(:,bny*ibox+1:bny*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine unmap_from_pencil_yz_3D
!***********************************************************************
    subroutine unmap_from_pencil_yz_4D(in, out,comm,lsync)
!
!  Unmaps z-pencil shaped 4D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocz>1.
!
!  27-oct-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
!
      integer, parameter :: iny=ny/nprocz, inz=nzgrid
      integer, parameter :: ony=ny, onz=nz
      integer :: inx, ina, onx, ona ! sizes of in and out arrays
      integer, parameter :: bny=ny/nprocz, bnz=nz ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ytag=111
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: send_buf, recv_buf
      integer, optional :: comm
      logical, optional :: lsync
!
      if (nprocz == 1) then
        !$omp workshare
        out = in
        !$omp end workshare
        return
      endif
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
      !$omp single
!
      inx = size (in, 1)
      ina = size (in, 4)
      onx = size (out, 1)
      ona = size (out, 4)
      nbox = onx*bny*bnz*ona
!
      if (mod (ny, nprocz) /= 0) &
          call stop_fatal ('unmap_from_pencil_yz_4D: ny needs to be an integer multiple of nprocz', lfirst_proc_yz)
!
      if ((size (in, 2) /= iny) .or. ((size (in, 3) /= inz))) &
          call stop_fatal ('unmap_from_pencil_yz_4D: input array size mismatch /= ny/nprocz,nzgrid', lfirst_proc_yz)
      if ((size (out, 2) /= ony) .or. ((size (out, 3) /= onz))) &
          call stop_fatal ('unmap_from_pencil_yz_4D: output array size mismatch /= ny,nz', lfirst_proc_yz)
      if (inx /= onx) &
          call stop_fatal ('unmap_from_pencil_yz_4D: inz/=onz - sizes differ in the x direction', lfirst_proc_yz)
      if (ina /= ona) &
          call stop_fatal ('unmap_from_pencil_yz_4D: ina/=ona - sizes differ in the 4th dimension', lfirst_proc_yz)
!
      allocate (send_buf(onx,bny,bnz,ona), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_yz_4D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(onx,bny,bnz,ona), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_yz_4D: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocz-1
        partner = find_proc(ipx,ipy,ibox)
        if (iproc == partner) then
          ! data is local
          out(:,bny*ibox+1:bny*(ibox+1),:,:) = in(:,:,bnz*ibox+1:bnz*(ibox+1),:)
        else
          ! communicate with partner
          send_buf = in(:,:,bnz*ibox+1:bnz*(ibox+1),:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), stat, mpierr)
            call MPI_SEND (send_buf, nbox, mpi_precision, partner, ytag, ioptest(comm,MPI_COMM_GRID), mpierr)
          endif
          out(:,bny*ibox+1:bny*(ibox+1),:,:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
      !$omp end single
      if (loptest(lsync,.true.)) then
        !$omp barrier
      endif
!
    endsubroutine unmap_from_pencil_yz_4D
!***********************************************************************
    subroutine collect_grid(x, y, z, gx, gy, gz)
!
!  This routine collects the global grid on the root processor.
!
!  05-Jul-2011/PABourdin: coded
!   6-oct-2016/MR: modifications for collection of downsampled grid:
!                  contributions from different processors may have different lengths
!
      use Cparam, only: l1, m1, n1, nghost
!
      real, dimension(:), intent(in) :: x,y,z
      real, dimension(:), intent(out) :: gx,gy,gz
!
      integer :: px, py, pz, mx, my, mz, l2, m2, n2, ie, rank
      integer, parameter :: tag_gx=677, tag_gy=678, tag_gz=679
!
      mx=size(x); my=size(y); mz=size(z)
      l2=mx-nghost; m2=my-nghost; n2=mz-nghost

      if (lroot) then
        ! collect the global x-data from all leading processors in the yz-plane
        gx(1:mx) = x; ie=l2
        if (nprocx > 1) then
          do px = 1, nprocx-1
            rank=find_proc(px,0,0)
            if (ldownsampling) call mpirecv_int(l2,rank,tag_gx)
            call mpirecv_real (gx(ie+1:ie+l2), l2, rank, tag_gx)
            ie=ie+l2-nghost
          enddo
        endif
        ! collect the global y-data from all leading processors in the xz-plane
        gy(1:my) = y; ie=m2
        if (nprocy > 1) then
          do py = 1, nprocy-1
            rank=find_proc(0,py,0)
            if (ldownsampling) call mpirecv_int(m2,rank,tag_gy)
            call mpirecv_real (gy(ie+1:ie+m2), m2, rank, tag_gy)
            ie=ie+m2-nghost
          enddo
        endif
        ! collect the global z-data from all leading processors in the xy-plane
        gz(1:mz) = z; ie=n2
        if (nprocz > 1) then
          do pz = 1, nprocz-1
            rank=find_proc(0,0,pz)
            if (ldownsampling) call mpirecv_int(n2,rank,tag_gz)
            call mpirecv_real (gz(ie+1:ie+n2), n2, rank, tag_gz)
            ie=ie+n2-nghost
          enddo
        endif
      else
        ! leading processors send their local coordinates
        if (lfirst_proc_yz) then
          if (ldownsampling) call mpisend_int(l2,root,tag_gx)
          call mpisend_real (x(l1:), l2, root, tag_gx)
        endif
        if (lfirst_proc_xz) then
          if (ldownsampling) call mpisend_int(m2,root,tag_gy)
          call mpisend_real (y(m1:), m2, root, tag_gy)
        endif
        if (lfirst_proc_xy) then
          if (ldownsampling) call mpisend_int(n2,root,tag_gz)
          call mpisend_real (z(n1:), n2, root, tag_gz)
        endif
      endif
!
    endsubroutine collect_grid
!***********************************************************************
    subroutine distribute_grid(x, y, z, gx, gy, gz)
!
!  This routine distributes the global grid to all processors
!  (concurrent communications).
!
!  16-jan-2023/MR: coded
!
      use General, only: indgen
!
      real, dimension(mx), intent(out) :: x
      real, dimension(my), intent(out) :: y
      real, dimension(mz), intent(out) :: z
      real, dimension(mxgrid), intent(in), optional :: gx
      real, dimension(mygrid), intent(in), optional :: gy
      real, dimension(mzgrid), intent(in), optional :: gz
!
      if (lfirst_proc_yz) then    ! x-beam through root
        if (lroot) then
          call mpiscatterv_real_plain(gx,spread(mx,1,nprocx),(indgen(nprocx)-1)*nx,x,mx,comm=MPI_COMM_XBEAM)
        else
          call mpiscatterv_real_plain(x,(/1/),(/1/),x,mx,comm=MPI_COMM_XBEAM)
        endif
      endif
!
      if (lfirst_proc_xz) then    ! y-beam through root
        if (lroot) then
          call mpiscatterv_real_plain(gy,spread(my,1,nprocy),(indgen(nprocy)-1)*ny,y,my,comm=MPI_COMM_YBEAM)
        else
          call mpiscatterv_real_plain(y,(/1/),(/1/),y,my,comm=MPI_COMM_YBEAM)
        endif
      endif
!
      if (lfirst_proc_xy) then    ! z-beam through root
        if (lroot) then
          call mpiscatterv_real_plain(gz,spread(mz,1,nprocz),(indgen(nprocz)-1)*nz,z,mz,comm=MPI_COMM_ZBEAM)
        else
          call mpiscatterv_real_plain(z,(/1/),(/1/),z,mz,comm=MPI_COMM_ZBEAM)
        endif
      endif

      call mpibcast_real(x,mx,comm=MPI_COMM_YZPLANE)
      call mpibcast_real(y,my,comm=MPI_COMM_XZPLANE)
      call mpibcast_real(z,mz,comm=MPI_COMM_XYPLANE)
!
    endsubroutine distribute_grid
!***********************************************************************
    subroutine y2x(a,xi,zj,zproc_no,ay)
!
!  Load the y dimension of an array in a 1-d array.
!
      real, dimension(nx,ny,nz), intent(in) :: a
      real, dimension(nygrid), intent(out) :: ay
      real, dimension(nygrid) :: ay_local
      integer, intent(in) :: xi,zj,zproc_no
      integer :: my_iniy,my_finy
!
      if (ipz==(zproc_no-1)) then
        my_iniy=ipy*ny+1
        my_finy=(ipy+1)*ny
        ay_local(my_iniy:my_finy)=a(xi,:,zj)
      else
        ay_local=0.
      endif
      call mpireduce_sum(ay_local,ay,nygrid)
!
    endsubroutine y2x
!***********************************************************************
    subroutine z2x(a,xi,yj,yproc_no,az)
!
!  Load the z dimension of an array in a 1-d array.
!
      real, dimension(nx,ny,nz), intent(in) :: a
      real, dimension(nzgrid), intent(out) :: az
      real, dimension(nzgrid) :: az_local
      integer, intent(in) :: xi,yj,yproc_no
      integer :: my_iniz,my_finz
!
      if (ipy==(yproc_no-1)) then
        my_iniz=ipz*nz+1
        my_finz=(ipz+1)*nz
        az_local(my_iniz:my_finz)=a(xi,yj,:)
      else
        az_local=0.
      endif
      call mpireduce_sum(az_local,az,nzgrid)
!
    endsubroutine z2x
!***********************************************************************
    subroutine mpigather_xy( sendbuf, recvbuf, lpz)
!
!  Gathers the chunks of a 2D array from each processor of the z-layer lpz in
!  a big array at the root of the layer. If lpz not present this is done for
!  all layers (not checked).
!
!  Here no parallelization in x allowed.
!
!  18-nov-10/MR: coded
!
      real, dimension(nxgrid,ny)     :: sendbuf   ! nx=nxgrid !
      real, dimension(nxgrid,nygrid) :: recvbuf
      integer, optional, intent(in)  :: lpz
!
      integer :: ncnt
      logical :: cond
      integer :: comm
!
      if (present(lpz)) then
        cond = ipz==lpz
      else
        cond = .true.
      endif
      comm = MPI_COMM_XYPLANE
!
      ncnt = nxgrid*ny
!
      if (cond) &
        call MPI_GATHER(sendbuf, ncnt, mpi_precision, recvbuf, ncnt, mpi_precision, root, comm, mpierr)
!
    endsubroutine mpigather_xy
 !***********************************************************************
    subroutine mpigather_z_1D(sendbuf,recvbuf)
!
!  Gathers the chunks of a 1D array from each processor along a z-beam and broadcasts it.
!
!  25-nov-10/MR: coded
!
      real, dimension(:), intent(in)  :: sendbuf
      real, dimension(:), intent(out) :: recvbuf

      call MPI_GATHER(sendbuf, size(sendbuf), mpi_precision, recvbuf, size(sendbuf), mpi_precision, 0, &
                      MPI_COMM_ZBEAM, mpierr)

    endsubroutine mpigather_z_1D
!***********************************************************************
    subroutine mpigather_z(sendbuf,recvbuf,n1,lproc)
!
!  Gathers the chunks of a 2D array from each processor along a z-beam at
!  position, defined by lproc at root of the beam.
!
!  25-nov-10/MR: coded
!
      integer,                    intent(in)  :: n1
      real, dimension(n1,nz)    , intent(in)  :: sendbuf
      real, dimension(n1,nzgrid), intent(out) :: recvbuf
      integer, optional,          intent(in)  :: lproc
!
      integer :: lpx, lpy
!
      if (present(lproc)) then
        lpy = lproc/nprocx
        lpx = mod(lproc,nprocx)
      else
        lpy=0; lpx=0
      endif
!
      if ( ipx==lpx .and. ipy==lpy ) &
        call MPI_GATHER(sendbuf, n1*nz, mpi_precision, recvbuf, n1*nz, mpi_precision, root, MPI_COMM_ZBEAM, mpierr)
!
    endsubroutine mpigather_z
!***********************************************************************
    subroutine mpigather( sendbuf, recvbuf , comm)
!
!  Gathers the chunks of a 3D array from each processor in a big array at root.
!
!  Here no parallelization in x allowed.
!
!  19-nov-10/MR: coded
!
      real, dimension(nxgrid,ny,nz):: sendbuf   ! nx=nxgrid !
      real, dimension(:,:,:)       :: recvbuf
      integer, optional :: comm
!
      integer :: ncnt, i
!
!  MR: These long integer variables would be necessary for big nxgrid*nygrid,
!      but there is no MPI_GATHERV which would accept a long int shifts argument.
!
      integer(KIND=ikind8) :: nlayer, nshift
      integer(KIND=ikind8), dimension(ncpus) :: shifts
      integer, dimension(ncpus) :: counts
!
      ncnt = nxgrid*ny
!
      if (lroot) then
!
        nlayer = nz*int8(nxgrid)*int8(nygrid)
        if (nlayer>max_int) &
          call stop_it("mpigather: integer overflow in shifts")
        counts = ncnt
!
        shifts(1) = 0
        nshift = nlayer
!
        do i=2,ncpus
!
          if ( mod(i,nprocy)==1 ) then
            shifts(i) = nshift
            nshift = nshift+nlayer
          else
            shifts(i) = shifts(i-1)+ncnt
          endif
!
        enddo
!
      endif
!
      do i=1,nz
        call MPI_GATHERV(sendbuf(1,1,i), ncnt, mpi_precision, recvbuf(1,1,i), counts, int(shifts), &
                         mpi_precision, root, ioptest(comm,MPI_COMM_GRID), mpierr)
      enddo
!
    endsubroutine mpigather
!***********************************************************************
    subroutine mpigather_scl_str(string,string_arr, comm)

      character(LEN=*) :: string
      character(LEN=*), dimension(:) :: string_arr
      integer, optional :: comm

      call MPI_GATHER(string, len(string), MPI_CHARACTER, string_arr, len(string), &
                      MPI_CHARACTER, root, ioptest(comm,MPI_COMM_GRID), mpierr)

    endsubroutine mpigather_scl_str
!***********************************************************************
    logical function get_limits(range, k1g, k2g, ia, ie, is )

      integer, dimension(3) :: range
      integer :: k1g, k2g
      integer :: ia, ie, is

      get_limits=.true.

      if ( range(1) == 0 ) return
!
      get_limits=.false.
      is=range(3)
      if (range(1)>=k1g) then
        ia = range(1)-k1g+1
      else
        ia = mod(is-mod(k1g-range(1),is),is)+1
      endif

      ie = min(k2g,range(2))-k1g+1

    endfunction get_limits
!***********************************************************************
    subroutine mpigather_and_out_real( sendbuf, unit, ltransp, kxrange, kyrange,zrange,comm)
!
!  Transfers the chunks of a 3D array from each processor to root
!  and writes them out in right order.
!
!  Here no parallelization in x allowed.
!
!  22-nov-10/MR: coded
!  06-apr-11/MR: optional parameters kxrange, kyrange, zrange for selective output added
!  03-feb-14/MR: rewritten
!  10-apr-15/MR: corrected for nx/=ny
!  22-apr-15/MR: another correction: order of loops over processors and index ranges exchanged
!
      use General, only: write_full_columns, get_range_no
!
      integer,                             intent(in   ) :: unit
      real,    dimension(:,:,:),           intent(inout) :: sendbuf
      complex, dimension(:,:,:,:),         intent(inout) :: sendbuf_cmplx
      logical,                    optional,intent(in   ) :: ltransp   ! if true, transposition x <-> y
      integer, dimension(3,*),    optional,intent(in   ) :: kxrange, kyrange, zrange
!
      integer, dimension(3,10) :: kxrangel,kyrangel,zrangel
      integer, optional :: comm
!
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer :: ipz, ipy, ipx, ic, ncomp, n1g, n2g, m1g, m2g, l1g, l2g, ig, iy, iz, &
                 irz, iry, irx, iza, ize, izs, iya, iye, iys, ixa, ixe, ixs, nsend, tag, unfilled
      logical :: ltrans, lcomplex
      real,    allocatable :: rowbuf(:)
      complex, allocatable :: rowbuf_cmplx(:)
      integer, dimension(2) :: nprxy,nxy
!
      lcomplex = .false.
      ncomp = 1
      goto 1
!
      entry mpigather_and_out_cmplx( sendbuf_cmplx, unit, ltransp, kxrange, kyrange,zrange, comm)
!
      lcomplex = .true.
      ncomp = size(sendbuf_cmplx,4)
!
1     if (ALWAYS_FALSE) print*,unit
!
      if ( .not.present(ltransp) ) then
        ltrans=.false.
      else
        ltrans=ltransp
      endif
!
      if ( .not.present(kxrange) ) then
        kxrangel = 0
        kxrangel(:,1) = (/1,nxgrid,1/)
      else
        kxrangel = kxrange(:,1:nk_max)
      endif
!
      if ( .not.present(kyrange) ) then
        kyrangel = 0
        kyrangel(:,1)=(/1,nygrid,1/)
      else
        kyrangel = kyrange(:,1:nk_max)
      endif
!
      if ( .not.present(zrange) ) then
        zrangel = 0
        zrangel(:,1) = (/1,nzgrid,1/)
      else
        zrangel = zrange(:,1:nz_max)
      endif
!
      if (lcomplex) then
        allocate( rowbuf_cmplx(nx) )
      else
        allocate( rowbuf(nx) )
      endif
!
      if (ltrans) then
        nxy(1)=ny; nxy(2)=nx
        nprxy(1)=nprocy; nprxy(2)=nprocx
      else
        nxy(1)=nx; nxy(2)=ny
        nprxy(1)=nprocx; nprxy(2)=nprocy
      endif
!
      unfilled=0
!
! loop over all variables for which spectrum has been calculated
!
      do ic=1,ncomp
!
! loop over all ranges of z indices in zrangel
!
        do irz=1,nz_max
!
          n2g=0
!
! loop over all processor array layers of z direction
!
          do ipz=0,nprocz-1
!
! global lower and upper z index bounds for z layer ipz
!
            n1g = n2g+1; n2g = n2g+nz
            if (get_limits( zrangel(:,irz), n1g, n2g, iza, ize, izs )) exit
!
! loop over all z indices in range irz
!
            do iz = iza, ize, izs
!
! loop over all ranges of ky indices in kyrangel
!
              do iry=1,nk_max
!
! loop over all processor array beams in x direction in layer ipz
!
                m2g=0

                do ipy=0,nprxy(2)-1
!
! global lower and upper y index bounds for beam ipy
!
                  m1g=m2g+1; m2g=m2g+nxy(2)
                  if (get_limits( kyrangel(:,iry), m1g, m2g, iya, iye, iys )) exit
                  !if (lroot) print*, 'ipy,ipz,iry,iy*=', ipy,ipz, iry, iya, iye, iys
!
! loop over all ky indices in range iry
!
                  do iy = iya, iye, iys
!
! loop over all ranges of kx indices in kxrangel
!
                    do irx=1,nk_max
!
!
! loop over all processors in beam
!
                      l2g=0
                      do ipx=0,nprxy(1)-1
!
! global processor number
!
                       if (ltrans) then
                          ig = ipz*nprocxy + ipx*nprocx + ipy
                        else
                          ig = ipz*nprocxy + ipy*nprocx + ipx
                        endif
!
! global lower and upper x index bounds for processor ipx
!
                        l1g=l2g+1; l2g=l2g+nxy(1)
                        if (get_limits( kxrangel(:,irx), l1g, l2g, ixa, ixe, ixs )) exit
                        !if (lroot) print*, 'ipx,ipy,ix*=', ipx,ipy,ixa, ixe, ixs
!
                        if (ixe>=ixa) then

                          nsend = get_range_no((/ixa,ixe,ixs/),1)
                          tag = ixa
!
                          if (lroot) then
                            if (ig==0) then
                              if (lcomplex) then
                                if (ltrans) then
                                  rowbuf_cmplx=sendbuf_cmplx(iy,ixa:ixe:ixs,iz,ic)
                                else
                                  rowbuf_cmplx=sendbuf_cmplx(ixa:ixe:ixs,iy,iz,ic)
                                endif
                              else
                                if (ltrans) then
                                  rowbuf=sendbuf(iy,ixa:ixe:ixs,iz)
                                else
                                  rowbuf=sendbuf(ixa:ixe:ixs,iy,iz)
                                endif
                              endif
                            else
                              if (lcomplex) then
                                call MPI_RECV(rowbuf_cmplx, nsend, MPI_CMPLX, ig, tag, ioptest(comm,MPI_COMM_GRID), status, mpierr)
                                !print*, 'iy,irx,ixa:ixe:ixs,kxrangel(:,irx)=', iy,irx,ixa,ixe,ixs , kxrangel(:,irx)
                              else
                                call MPI_RECV(rowbuf, nsend, mpi_precision, ig, tag, ioptest(comm,MPI_COMM_GRID), status, mpierr)
                              endif
                            endif
                            if (lcomplex) then
                              call write_full_columns( 1, rowbuf_cmplx, (/1,nsend,1/), unfilled )
                            else
                              call write_full_columns( 1, rowbuf, (/1,nsend,1/), unfilled )
                            endif
                          else if ( iproc==ig ) then       ! if executing processor is hit by index ig: send to root
!
                            if (lcomplex) then
                              if (ltrans) then
                                call MPI_SEND(sendbuf_cmplx(iy,ixa:ixe:ixs,iz,ic), &
                                              nsend, MPI_CMPLX, root, tag, ioptest(comm,MPI_COMM_GRID), mpierr)
                              else
                                call MPI_SEND(sendbuf_cmplx(ixa:ixe:ixs,iy,iz,ic), &
                                              nsend, MPI_CMPLX, root, tag, ioptest(comm,MPI_COMM_GRID), mpierr)
                              endif
                            else
                              if (ltrans) then
                                call MPI_SEND(sendbuf(iy,ixa:ixe:ixs,iz), &
                                              nsend, mpi_precision, root, tag, ioptest(comm,MPI_COMM_GRID), mpierr)
                              else
                                call MPI_SEND(sendbuf(ixa:ixe:ixs,iy,iz), &
                                              nsend, mpi_precision, root, tag, ioptest(comm,MPI_COMM_GRID), mpierr)
                              endif
                            endif
!
                          endif
                        endif
                      enddo
                    enddo
                    call MPI_BARRIER(ioptest(comm,MPI_COMM_GRID), mpierr)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
!
        if (lroot.and.unfilled>0) then
          write(1,'(a)')
          unfilled=0
        endif
!
      enddo
!
    endsubroutine mpigather_and_out_real
!***********************************************************************
    subroutine mpiwait(bwait)
!
!  Wrapper to go with non-blocking communication functions.
!
!  12-dec-14/wlad: coded
!
      integer, dimension (MPI_STATUS_SIZE) :: stat
      integer :: bwait,ierr
!
      call MPI_WAIT(bwait,stat,ierr)

    endsubroutine mpiwait
!***********************************************************************
    subroutine merge_1d( vec1, vec2, n, type )
!
!  Helper function for mpimerge_1d.
!
!  22-nov-10/MR: coded
!
      integer,            intent(in)    :: n, type
      real, dimension(n), intent(inout) :: vec2
      real, dimension(n), intent(in)    :: vec1
!
    ! merging
      where ((vec2 < 0.) .and. (vec1 >= 0.)) vec2=vec1
!
      if (ALWAYS_FALSE) print *,type
!
    endsubroutine merge_1d
!***********************************************************************
    subroutine mpimerge_1d(vector,nk,idir,comm)
!
!  Merges vectors of processors along idir by filling invalid values (NaN).
!
!  22-nov-10/MR: coded
!
      integer,             intent(in)    :: nk
      real, dimension(nk), intent(inout) :: vector
      integer, optional,   intent(in)    :: idir,comm
!
      integer                            :: mpiprocs,merge
      real, dimension(nk)                :: recvbuf
!
      if (nk==0) return
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=ioptest(comm,MPI_COMM_GRID)
      endif
!
      call MPI_OP_CREATE( merge_1d, .false., merge, mpierr )
      call MPI_REDUCE(vector, recvbuf, nk, mpi_precision, merge, root, mpiprocs, mpierr)
      vector = recvbuf
!
    endsubroutine mpimerge_1d
!***********************************************************************
    integer function mpigetcomm(idir)
!
!  Derives communicator from index idir.
!
!  23-nov-10/MR: coded
!
      integer, intent(in) :: idir
!
      select case(idir)
        case(IXBEAM)
          mpigetcomm=MPI_COMM_XBEAM
        case(IYBEAM)
          mpigetcomm=MPI_COMM_YBEAM
        case(IZBEAM)
          mpigetcomm=MPI_COMM_ZBEAM
        case(IXYPLANE)
          mpigetcomm=MPI_COMM_XYPLANE
        case(IXZPLANE)
          mpigetcomm=MPI_COMM_XZPLANE
        case(IYZPLANE)
          mpigetcomm=MPI_COMM_YZPLANE
        case default
          mpigetcomm=MPI_COMM_GRID
      endselect
!
    endfunction mpigetcomm
!************************************************************************
    logical function report_clean_output( flag, message )
!
!  Generates error message for files (like snapshots) which are distributedly
!  written. Message contains list of processors at which operation failed.
!
!  flag   (IN) : indicates failure of I/O operation at local processor
!  message(OUT): message fragment containing list of processors where
!                operation failed  (only relevant for root)
!  return value: flag for 'synchronize!', identical for all processors
!
!  14-nov-11/MR: coded
!
      use General, only: itoa,safe_character_append,safe_character_prepend
!
      logical,               intent(IN) :: flag
      character (LEN=fnlen), intent(OUT):: message
!
      integer :: mpierr, i, ia, ie, count
      logical, dimension(:), allocatable:: flags
!
      character (LEN=intlen)  :: str
!
      if (lroot) allocate(flags(ncpus))
!
      call MPI_GATHER(flag, 1, MPI_LOGICAL, flags, 1, MPI_LOGICAL, root, MPI_COMM_PENCIL, mpierr)
!
      report_clean_output = .false.
      if (lroot) then
!
        count = 0
        ia = -1; ie = 0
        str = ''; message = ''
!
        if (lroot) then
          do i=1,ncpus
            if ( flags(i) ) then
              if ( i==ia+1 )  then
                ie = i
                count = count+1
              else
                ia = i
                if ( ie>0 ) then
                  str = '-'//trim(itoa(ie))//','
                  ie = 0
                endif
                call safe_character_append(message,trim(itoa(ia))//str)
                str = ','
              endif
            endif
          enddo
        endif
!
        if (count>0) then
          call safe_character_prepend(message,'"at '//trim(itoa(count))//' node(s): ')
          report_clean_output = .true.
        endif
!
      endif
!
! Broadcast flag for 'sychronization necessary'.
!
      call MPI_BCAST(report_clean_output,1,MPI_LOGICAL,root,MPI_COMM_PENCIL,mpierr)
!
    end function report_clean_output
!**************************************************************************
    subroutine set_yy_neighbors
!
!  Determines for the boundary processors of one (Yin or Yang) grid the neighbors
!  (direct or corner) in the other. Determines further the length of the buffers
!  for data exchange between the neighbors and stores them in yy_buflens.
!  Implemented for nprocz/nprocy=3 only!
!
!  At return neighbor ranks refer to MPI_COMM_PENCIL.
!
!  20-dec-15/MR: coded
!
        integer, parameter :: nprocz_rd=nprocz/3
        integer :: lenred

        if (lcutoff_corners) then
          len_cornstrip_y=nycut-1
          len_cornstrip_z=nzcut-1
!
! length of a "cornerstrip" if corner is cut off. Its inclined part has width 2*nghost!
!
          len_cornbuf=max(2*my-nycut+len_cornstrip_z,2*mz-nzcut+len_cornstrip_y)
          ycornstart=2*my-nycut+1; zcornstart=2*mz-nzcut+1
          lenred=max(2*my-nycut,2*mz-nzcut)
        else
          len_cornstrip_y=m2; len_cornstrip_z=n2
          len_cornbuf=my+nz+nghost                                ! length of a "cornerstrip"
          ycornstart=my+1; zcornstart=mz+1
          lenred=max(my,mz)
        endif
!
!  At lower theta boundary (northern hemisphere).
!
        if (lfirst_proc_y) then
          if (ipz<nprocz_rd) then                         ! for first third in phi direction

            yy_buflens=mz                                 ! extent along boundary with ghost zones.
!
!  Direct neighbor.
!
            ylneigh=find_proc(ipx,nprocy-1,ipz+2*nprocz_rd)
!
!  Lower left corner neighbor.
!
            if (lfirst_proc_z) then                       ! not used at yz-corner
              llcorn=-1
            else
              llcorn=find_proc(ipx,nprocy-1,ipz+2*nprocz_rd-1)
            endif
!
!  Lower right corner neighbor.
!
            if (ipz==nprocz_rd-1) then
              lucorn=-1                                   ! last proc in first third = corner opposite
              yy_buflens=(/mz,len_cornbuf,0/)
            else
              lucorn=find_proc(ipx,nprocy-1,ipz+2*nprocz_rd+1)
              if (ipz==nprocz_rd-2) yy_buflens(3)=lenred  ! penultimate proc in first third = cross corner opposite
            endif

          elseif (ipz>=2*nprocz_rd) then                  ! for last third in phi direction

            yy_buflens=mz
!
!  Direct neighbor.
!
            ylneigh=find_proc(ipx,0,5*nprocz_rd-ipz-1)
!
!  Lower left corner neighbor.
!
            if (ipz==2*nprocz_rd) then                    ! first proc in last third = corner opposite
              llcorn=-1
              yy_buflens=(/0,len_cornbuf,mz/)
            else
              llcorn=find_proc(ipx,0,5*nprocz_rd-ipz)
              if (ipz==2*nprocz_rd+1) yy_buflens(1)=lenred ! 2nd proc in last third = cross corner opposite
            endif
!
!  Upper left corner neighbor.
!
            if (llast_proc_z) then                        ! not used at yz-corner
              lucorn=-1
            else
              lucorn=find_proc(ipx,0,5*nprocz_rd-ipz-2)
            endif

          else                                            ! for middle third in phi direction

            yy_buflens=my
!
!  Direct neighbor.
!
            ylneigh=find_proc(ipx,nprocy-1-(ipz-nprocz_rd),nprocz-1)
!
!  Lower left corner neighbor.
!
            if (ipz==nprocz_rd) then                      ! first proc in middle third = corner opposite
              llcorn=-1
              yy_buflens=(/0,len_cornbuf,my/)
            else
              llcorn=find_proc(ipx,nprocy-(ipz-nprocz_rd),nprocz-1)
              if (ipz==nprocz_rd+1) yy_buflens(1)=lenred  ! 2nd proc in middle third = cross corner opposite
            endif
!
!  Upper left corner neighbor.
!
            if (ipz==2*nprocz_rd-1) then                  ! last proc in middle third = corner opposite
              lucorn=-1
              yy_buflens=(/my,len_cornbuf,0/)
            else
              lucorn=find_proc(ipx,nprocy-2-(ipz-nprocz_rd),nprocz-1)
              if (ipz==2*nprocz_rd-2) yy_buflens(3)=lenred  ! penultimate proc in middle third = cross corner opposite
            endif

          endif

          if (.not.lyang) then
!
!  For Yin grid: neighbors are in upper half of overall proc number range.
!
            ylneigh=ylneigh+ncpus
            if (llcorn/=-1) llcorn=llcorn+ncpus
            if (lucorn/=-1) lucorn=lucorn+ncpus
          endif

        endif
!
!  At upper theta boundary (southern hemisphere).
!
        if (llast_proc_y) then
          if (ipz<nprocz_rd) then                              ! for first third in phi direction

            yy_buflens=mz
!
!  Direct neighbor.
!
            yuneigh=find_proc(ipx,nprocy-1,nprocz_rd-1-ipz)
!
!  Upper right neighbor.
!
            if (lfirst_proc_z) then
              ulcorn=-1                                        ! not used at yz-corner
            else
              ulcorn=find_proc(ipx,nprocy-1,nprocz_rd-ipz)
            endif
!
!  Upper left neighbor.
!
            if (ipz==nprocz_rd-1) then                         ! last proc in first third = corner opposite
              uucorn=-1
              yy_buflens=(/0,len_cornbuf,mz/)
            else
              uucorn=find_proc(ipx,nprocy-1,nprocz_rd-2-ipz)
              if (ipz==nprocz_rd-2) yy_buflens(1)=lenred       ! penultimate proc in first third = cross corner opposite
            endif

          elseif (ipz>=2*nprocz_rd) then                       ! for last third in phi direction

            yy_buflens=mz
!
!  Direct neighbor.
!
            yuneigh=find_proc(ipx,0,ipz-2*nprocz_rd)
!
!  Upper right neighbor.
!
            if (ipz==2*nprocz_rd) then                         ! first proc in last third = corner opposite
              ulcorn=-1
              yy_buflens=(/mz,len_cornbuf,0/)
            else
              ulcorn=find_proc(ipx,0,ipz-2*nprocz_rd-1)
              if (ipz==2*nprocz_rd+1) yy_buflens(3)=lenred     ! 2nd proc in last third = cross corner opposite
            endif
!
!  Upper left neighbor.
!
            if (llast_proc_z) then
              uucorn=-1                                       ! not used at yz-corner
            else
              uucorn=find_proc(ipx,0,ipz-2*nprocz_rd+1)
            endif

          else                                                ! middle third in phi direction
!
!  Direct neighbor.
!
            yy_buflens=my

            yuneigh=find_proc(ipx,nprocy-1-(ipz-nprocz_rd),0)
!
!  Upper right neighbor.
!
            if (ipz==nprocz_rd) then                          ! first proc in middle third = corner opposite
              ulcorn=-1
              yy_buflens=(/my,len_cornbuf,0/)
            else
              ulcorn=find_proc(ipx,nprocy-(ipz-nprocz_rd),0)
              if (ipz==nprocz_rd+1) yy_buflens(3)=lenred      ! 2nd proc in middle third = cross corner opposite
            endif
!
!  Upper left neighbor.
!
            if (ipz==2*nprocz_rd-1) then                      ! last proc in middle third = corner opposite
              uucorn=-1
              yy_buflens=(/0,len_cornbuf,my/)
            else
              uucorn=find_proc(ipx,nprocy-2-(ipz-nprocz_rd),0)
              if (ipz==2*nprocz_rd-2) yy_buflens(1)=lenred    ! penultimate proc in middle third = cross corner opposite
            endif

          endif

          if (.not.lyang) then
!
!  For Yin grid: neighbors are in upper half of overall proc number range.
!
            yuneigh=yuneigh+ncpus
            if (ulcorn/=-1) ulcorn=ulcorn+ncpus
            if (uucorn/=-1) uucorn=uucorn+ncpus
          endif

        endif
!
!  At lower phi boundary.
!
        if (lfirst_proc_z) then
!
!  Direct neighbor.
!
          zlneigh=find_proc(ipx,nprocy-1,2*nprocz_rd-1-ipy)
!
!  Lower right neighbor.
!
          if (lfirst_proc_y) then
            llcorn=-1                                       ! not used at yz-corner
          else
            llcorn=find_proc(ipx,nprocy-1,2*nprocz_rd-ipy)
          endif
!
!  Lower left neighbor.
!
          if (llast_proc_y) then
            ulcorn=-1                                       ! not used at yz-corner
          else
            ulcorn=find_proc(ipx,nprocy-1,2*nprocz_rd-2-ipy)
          endif

          yy_buflens=mz

          if (.not.lyang) then
!
!  For Yin grid: neighbors are in upper half of overall proc number range.
!
            zlneigh=zlneigh+ncpus
            if (llcorn/=-1) llcorn=llcorn+ncpus
            if (ulcorn/=-1) ulcorn=ulcorn+ncpus
          endif

        endif
!
!  At upper phi boundary.
!
        if (llast_proc_z) then
!
!  Direct neighbor.
!
          zuneigh=find_proc(ipx,0,2*nprocz_rd-1-ipy)
!
!  Upper left neighbor.
!
          if (lfirst_proc_y) then
            lucorn=-1                                 ! not used at yz-corner
          else
            lucorn=find_proc(ipx,0,2*nprocz_rd-ipy)
          endif
!
!  Upper right neighbor.
!
          if (llast_proc_y) then
            uucorn=-1                                 ! not used at yz-corner
          else
            uucorn=find_proc(ipx,0,2*nprocz_rd-2-ipy)
          endif

          yy_buflens=mz

          if (.not.lyang) then
!
!  For Yin grid: neighbors are in upper half of overall proc number range.
!
            zuneigh=zuneigh+ncpus
            if (lucorn/=-1) lucorn=lucorn+ncpus
            if (uucorn/=-1) uucorn=uucorn+ncpus
          endif

        endif
!if (.not.lyang) then
if (.false.) then
  write(iproc_world+10,*) 'ylneigh,yuneigh=',ylneigh,yuneigh
  write(iproc_world+10,*) 'zlneigh,zuneigh=',zlneigh,zuneigh
  write(iproc_world+10,*) 'llcorn,lucorn=',llcorn,lucorn
  write(iproc_world+10,*) 'ulcorn,uucorn=',ulcorn,uucorn
endif
    endsubroutine set_yy_neighbors
!***********************************************************************
    subroutine interpolate_yy(f,ivar1,ivar2,buffer,pos,type)
!
!  Performs bilinear interpolation for a whole ghost zone strip of the
!  complementary (Yin or Yang) grid and all variables ivar1:ivar2. Selection
!  of the strip (from direct, left or right corner neighbor) by pos. Parameter
!  type is for later use with other interpolation schemes.
!  Result returned in buffer.
!
!  20-dec-15/MR: coded
!
     use General, only: var_is_vec
     use Yinyang, only: bilin_interp, biquad_interp

      real, dimension(:,:,:,:), intent(IN) :: f
      integer,                  intent(IN) :: ivar1, ivar2, pos, type
      real, dimension(:,:,:,:), intent(OUT):: buffer

      integer :: nth, nph, i, j, iv, ive, ibuf, jbuf
      logical :: ltransp

      if (pos==NIL) return

      nth=size(intcoeffs(pos)%inds,1); nph=size(intcoeffs(pos)%inds,2)
      ltransp = size(buffer,3)==nth .and. size(buffer,2)==nph

      buffer=0.
!
!  Loop over all variables.
!
      iv=ivar1
      do while (iv<=ivar2)

        if (var_is_vec(iv)) then
!
!  All components of vectorial variables are interpolated at the same time as transformation is
!  needed.
!
          ive=iv+2
        else
          ive=iv
        endif
!if (.false..and.(iproc_world==52.or.iproc_world==53).and.ive==-1) then
!  if (intcoeffs(pos)%pcoors(15,1,1)/=0.) then
!    print*, intcoeffs(pos)%pcoors(15,1,:)
!    print*,intcoeffs(pos)%coeffs2(15,1,:,:)
!    print*,'----'
!  endif
!endif

        do i=1,nth; do j=1,nph

          if (ltransp) then
            ibuf=j; jbuf=i
          else
            ibuf=i; jbuf=j
          endif
          if (type==BILIN) then
            call bilin_interp(intcoeffs(pos),i,j,f(:,:,:,iv:ive), buffer(:,:,:,iv:ive),ibuf,jbuf)
          elseif (type==BIQUAD .or. type==BICUB .or. type==BIQUIN .or. type==QUADSPLINE) then
            call biquad_interp(intcoeffs(pos),i,j,f(:,:,:,iv:ive), buffer(:,:,:,iv:ive),ibuf,jbuf)
          else
            call stop_it('interpolate_yy: Only bilinear, biquadratic, bicubic and biquintic '// &
                         'interpolations implemented')
          endif

        enddo; enddo
        iv=ive+1

      enddo
!
    endsubroutine interpolate_yy
!***********************************************************************
    subroutine mpiscatterv_real_plain(src,counts,dspls,dest,nlocal,comm)

      real, dimension(:) :: src, dest
      integer, dimension(:) :: counts,dspls
      integer :: nlocal
      integer, optional :: comm

      call MPI_SCATTERV(src,counts,dspls,mpi_precision,dest,nlocal,mpi_precision, &
                        root,ioptest(comm,MPI_COMM_PENCIL),mpierr)

    endsubroutine mpiscatterv_real_plain
!***********************************************************************
    function mpiscatterv_real(nlocal,src,dest) result (lerr)

      integer :: nlocal
      real, dimension(:) :: src, dest
      logical :: lerr

      integer, dimension(ncpus) :: counts,dspls
      integer :: i

      call MPI_GATHER(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,root,MPI_COMM_PENCIL,mpierr)

      if (lroot) then
        lerr=(sum(counts)>size(src))
        if (.not.lerr) then
          dspls(1)=0
          do i=2,ncpus
            dspls(i)=dspls(i-1)+counts(i-1)
          enddo
        endif
      endif

      call mpibcast(lerr)
      if (lerr) then
        return
      else
        call MPI_SCATTERV(src,counts,dspls,mpi_precision,dest,nlocal,mpi_precision,root,MPI_COMM_PENCIL,mpierr)
      endif

    endfunction mpiscatterv_real
!***********************************************************************
    function mpiscatterv_int(nlocal,src,dest) result (lerr)

      integer :: nlocal
      integer, dimension(:) :: src, dest
      logical :: lerr

      integer, dimension(ncpus) :: counts,dspls
      integer :: i

      call MPI_GATHER(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,root,MPI_COMM_PENCIL,mpierr)

      if (lroot) then
        lerr=(sum(counts)>size(src))
        if (.not.lerr) then
          dspls(1)=0
          do i=2,ncpus
            dspls(i)=dspls(i-1)+counts(i-1)
          enddo
        endif
      endif
      call mpibcast(lerr)
      if (lerr) then
        return
      else
        call MPI_SCATTERV(src,counts,dspls,MPI_INTEGER,dest,nlocal,MPI_INTEGER,root,MPI_COMM_PENCIL,mpierr)
      endif

    endfunction mpiscatterv_int
!***********************************************************************
    subroutine initialize_foreign_comm(frgn_buffer)
!
! Initializes communication with foreign code by fetching the foreign setup
! and establishing relationships.
! Partly specific for coupling with MagIC.
!
! 20-oct-21/MR: coded
!
      use General, only: itoa, find_index_range, find_proc_general, get_linterp_weights_1D

      real, dimension(:,:,:,:), allocatable :: frgn_buffer
!
!  Variables for data transfer from foreign code.
!
      integer, dimension(3) :: intbuf
      real, dimension(6) :: floatbuf
      logical :: lok
      character(LEN=128) :: messg
      integer :: ind,j,ll,name_len,nxgrid_foreign,nygrid_foreign,il1,il2,im1,im2,lenx,leny,px,py,tag,peer

      if (lforeign) then

        frgn_setup%root=ncpus
        frgn_setup%tag=tag_foreign
        frgn_setup%t_last_recvd=t
!print*, 'iproc, iproc_world, MPI_COMM_UNIVERSE, MPI_COMM_WORLD=', iproc, &
!        iproc_world, MPI_COMM_UNIVERSE, MPI_COMM_WORLD
        if (lroot) then
!
          lok=.true.; messg=''
!
!  Send communication tag to foreign code.
!
         call mpisend_int(tag_foreign,frgn_setup%root,0,MPI_COMM_WORLD)
!
!  Receive length of name of foreign code.
!
          call mpirecv_int(name_len,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)

          if (name_len<=0) then
            call stop_it('initialize_foreign_comm: length of foreign name <=0 or >')
          elseif (name_len>labellen) then
            call stop_it('initialize_foreign_comm: length of foreign name > labellen ='// &
                         trim(itoa(labellen)))
          endif
!
!  Receive name of foreign code.
!
          call mpirecv_char(frgn_setup%name(1:name_len),frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
          frgn_setup%name=frgn_setup%name(1:name_len)

          if (.not.(trim(frgn_setup%name)=='MagIC'.or.trim(frgn_setup%name)=='EULAG')) &
            call stop_it('initialize_foreign_comm: communication with foreign code "'// &
                         trim(frgn_setup%name)//'" not supported')
!
!  Receive processor numbers of foreign code.
!
          call mpirecv_int(intbuf,3,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
          frgn_setup%procnums=intbuf
          frgn_setup%ncpus=product(intbuf)
          !if (ncpus+frgn_setup%ncpus/=nprocs) &
          !  call fatal_error('initialize_hydro','no of processors '//trim(itoa(nprocs))// &
          !          ' /= no of own + no of foreign processors '//trim(itoa(ncpus+frgn_setup%ncpus)))
!
!	Calculate the multiplicity of processors
!
          frgn_setup%proc_multis = (/nprocx,nprocy,nprocz/)/intbuf
          if ( any(mod((/nprocx,nprocy,nprocz/),intbuf)/=0 ) ) then
            messg="foreign proc numbers don't match;"
            lok=.false.
          endif
!
!  Receive gridpoint numbers of foreign code.
!
          call mpirecv_int(frgn_setup%dims,3,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
!print*, 'Pencil: received frgn_setup%dims=', frgn_setup%dims
          frgn_setup%linterpol = (frgn_setup%dims/=(/nxgrid,nygrid,nzgrid/))
!
!  Receive units of foreign code.
!
          call mpirecv_char(frgn_setup%unit_system,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
          if (trim(frgn_setup%unit_system)/='SI'.and.trim(frgn_setup%unit_system)/='CGS') then
            messg=trim(messg)//' foreign unit system unknown;'
            lok=.false.
          endif
          call mpirecv_real(frgn_setup%unit_length,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
          if (frgn_setup%unit_length<=0.) then
            messg=trim(messg)//' foreign length unit <=0;'
            lok=.false.
          endif
          call mpirecv_real(frgn_setup%unit_time,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
          if (frgn_setup%unit_time<=0.) then
            messg=trim(messg)//' foreign time unit <=0;'
            lok=.false.
          endif
          call mpirecv_real(frgn_setup%unit_BB,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
          if (frgn_setup%unit_BB<=0.) then
            messg=trim(messg)//' foreign B unit <=0;'
            lok=.false.
          endif
          call mpirecv_real(frgn_setup%unit_T,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
          if (frgn_setup%unit_T<=0.) then
            messg=trim(messg)//' foreign temperature unit <=0;'
            lok=.false.
          endif
          if (lok) frgn_setup%renorm_UU=frgn_setup%unit_time/frgn_setup%unit_length  !not the full truth yet
!
!  Receive domain extents of foreign code. j loops over r, theta, phi.
!
          call mpirecv_real(floatbuf,6,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
          ind=1
          do j=1,3
            frgn_setup%extents(:,j)=floatbuf(ind:ind+1)
            if (j==1) then
              floatbuf(1)=floatbuf(1)/floatbuf(2)*xyz1(1)     ! renormalize r-range for spherical coordinates
              floatbuf(2)=xyz1(1)
            endif
!print*, 'floatbuf=', floatbuf(ind:ind+1), xyz0(j),xyz1(j)
            if (any(abs(floatbuf(ind:ind+1)-(/xyz0(j),xyz1(j)/))>1.e-6)) then
               if (j==3) then
                 messg=trim(messg)//"foreign "//trim(coornames(j))//" z-domain extents doesn't match;"
                 lok=.false.
               else
                 frgn_setup%linterpol(j)=.true.
                 print*, "initialize_foreign_comm: Warning -- foreign "//trim(coornames(j))//" domain extent doesn't match"
               endif
            endif
            ind=ind+2
          enddo
!
!  Receive output timestep of foreign code.
!
          call mpirecv_real(frgn_setup%dt_out,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
          if (frgn_setup%dt_out<=0.) then
            messg=trim(messg)//' foreign output step<=0;'
            lok=.false.
          endif
!
!  Send confirmation flag that setup is acceptable.
!
          call mpisend_logical(lok,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
        endif

        call mpibcast_logical(lok,comm=MPI_COMM_PENCIL)

        if (.not.lok) call stop_it('initialize_foreign_comm: '//trim(messg))
        call mpibcast_char(frgn_setup%name,comm=MPI_COMM_PENCIL)
!
!  Broadcast foreign processor numbers and grid size.
!
        call mpibcast_int(frgn_setup%dims,3,comm=MPI_COMM_PENCIL)

        if (lfirst_proc_yz) then
          nxgrid_foreign=frgn_setup%dims(1)
          allocate(frgn_setup%xgrid(nxgrid_foreign))
        endif
        nygrid_foreign=frgn_setup%dims(2)
        allocate(frgn_setup%ygrid(nygrid_foreign))

        if (lroot) then
!
!  Receive vector of foreign global x(r)-grid points.
!
          call mpirecv_real(frgn_setup%xgrid,nxgrid_foreign,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
!
!  Renormalize X-coord to Pencil domain.
!
          frgn_setup%xgrid = frgn_setup%xgrid/frgn_setup%xgrid(nxgrid_foreign)*xyz1(1)
!print*, 'frgn_setup%xgrid=', frgn_setup%xgrid
!
!  Receive vector of foreign global y(theta)-grid points.
!
          call mpirecv_real(frgn_setup%ygrid,nygrid_foreign,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
!
!  Receive normalized time from foreign side
!
          call mpirecv_real(frgn_setup%renorm_t,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
!print*,'PENCIL - dt, dt/tnorm',frgn_setup%dt_out ,frgn_setup%dt_out/frgn_setup%renorm_t
          frgn_setup%renorm_L = frgn_setup%extents(2,1)
          frgn_setup%renorm_UU=frgn_setup%renorm_L/frgn_setup%renorm_t
          frgn_setup%dt_out = frgn_setup%dt_out/frgn_setup%renorm_t
!print*, 'PENCIL - renormUU, t_last_rcv, dt_out', frgn_setup%renorm_UU,frgn_setup%t_last_recvd,frgn_setup%dt_out, t-frgn_setup%t_last_recvd
!
!  Send number of xyz-procs to foreign.
!
          !call mpisend_int((/nprocx,nprocy,nprocz/),frgn_setup%root,tag_foreign,MPI_COMM_WORLD)  !!!
          call mpisend_int(nprocx,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
          call mpisend_int(nprocy,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
          call mpisend_int(nprocz,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
!
!print*, 'frgn_setup%linterpol=', frgn_setup%linterpol
        endif

        call mpibcast_real(frgn_setup%renorm_L,comm=MPI_COMM_PENCIL)
        call mpibcast_real(frgn_setup%renorm_t,comm=MPI_COMM_PENCIL)
        call mpibcast_real(frgn_setup%renorm_UU,comm=MPI_COMM_PENCIL)
        call mpibcast_real(frgn_setup%dt_out,comm=MPI_COMM_PENCIL)
!
        call mpibcast_int(frgn_setup%procnums,3,comm=MPI_COMM_PENCIL)
        call mpibcast_int(frgn_setup%proc_multis,3,comm=MPI_COMM_PENCIL)
        call mpibcast(frgn_setup%linterpol,3,comm=MPI_COMM_PENCIL)
!
!        if (frgn_setup%procnums(1)==1) then   !EULAG case
!print*, 'PENCIL: frgn_setup%proc_multis=', frgn_setup%proc_multis

          !frgn_setup%peer_rng = find_proc_general(ipz/frgn_setup%proc_multis(3), &
          !                                        ipy/frgn_setup%proc_multis(2), &
          !                                        ipx/frgn_setup%proc_multis(1), &
          !                      frgn_setup%procnums(3), frgn_setup%procnums(2), frgn_setup%procnums(1),.true.)
          !peer=frgn_setup%peer_rng(1)
!print*,'PENCIL - peer=', iproc, frgn_setup%peer_rng(1), &
!ipz/frgn_setup%proc_multis(3), ipy/frgn_setup%proc_multis(2), ipx/frgn_setup%proc_multis(1)

!        endif
!
!  Allocate array for index ranges w.r.t. frgn_setup%[xy]grid for each foreign
!  proc. frgn_setup%[xy]ind_rng(-1,:) is overall index range.
!
        allocate(frgn_setup%xind_rng(-1:frgn_setup%procnums(1)-1,2)); frgn_setup%xind_rng=0
        allocate(frgn_setup%yind_rng(-1:frgn_setup%procnums(2)-1,2)); frgn_setup%yind_rng=0
!
        if (lfirst_proc_xz) then             ! on processors of first YBEAM
!
!  Broadcast foreign ygrid to all procs with ipx=ipz=0.
!
          call mpibcast_real(frgn_setup%ygrid,nygrid_foreign,comm=MPI_COMM_YBEAM)
!
!  Determine index range frgn_setup%yind_rng in foreign ygrid which is needed
!  for individual processors in y direction.
!
          call find_index_range(frgn_setup%ygrid,nygrid_foreign,y(1),y(my),im1,im2,lextend=.true.)
          frgn_setup%yind_rng(-1,:)=(/im1,im2/)     ! global y index range of rank iproc
!
!  Ask all foreign processors in first y beam
!  about their share in frgn_setup%yind_rng. No share: receive [0,0].
!
          do py=0,frgn_setup%procnums(2)-1
!
!       Determine the peer using PENCIL (Not EULAG) proc grid conventions.
!
            peer = find_proc_general(0,py,0,frgn_setup%procnums(3),frgn_setup%procnums(2),frgn_setup%procnums(1),.true.)
            call mpisendrecv_int(frgn_setup%yind_rng(-1,:),2,ncpus+peer,tag_foreign+iproc, &
                                 intbuf,ncpus+peer,tag_foreign+peer,MPI_COMM_WORLD)
            frgn_setup%yind_rng(py,:)=intbuf(1:2)

!print*,'PENCIL yind1', iproc, peer,ipy,frgn_setup%yind_rng(py,1),frgn_setup%yind_rng(py,2)
!print*,'PENCIL yind2',iproc,frgn_setup%ypeer_rng(1),frgn_setup%ypeer_rng(2)
            if (frgn_setup%ypeer_rng(1)>=0) then
              if (frgn_setup%yind_rng(py,1)>=frgn_setup%yind_rng(py,2).and.frgn_setup%ypeer_rng(2)<0) then
                frgn_setup%ypeer_rng(2)=py-1
              endif
            else
              if (frgn_setup%yind_rng(py,1)<frgn_setup%yind_rng(py,2)) frgn_setup%ypeer_rng(1)=py
            endif
!print*,'PENCIL yind3',iproc,frgn_setup%ypeer_rng(1),frgn_setup%ypeer_rng(2)
          enddo
          if (frgn_setup%ypeer_rng(2)<0) frgn_setup%ypeer_rng(2)=py-1

!print*,'PENCIL ypeer', iproc,frgn_setup%ypeer_rng
!
          if (frgn_setup%linterpol(2)) &
            call get_linterp_weights_1D(frgn_setup%ygrid,frgn_setup%yind_rng(-1,:),y,frgn_setup%yinds,frgn_setup%yweights)
        endif   !lfirst_proc_xz
!
        if (lfirst_proc_yz) then             ! on processors of first XBEAM, i.e., ipy=ipz=0
!
!  Broadcast frgn_setup%xgrid to all procs of first x-beam.
!
          call mpibcast_real(frgn_setup%xgrid,nxgrid_foreign,comm=MPI_COMM_XBEAM)

!  Determine index range frgn_setup%xind_rng in frgn_setup%xgrid which is needed for individual
!  processors in x direction.
!
          call find_index_range(frgn_setup%xgrid,nxgrid_foreign,x(1),x(mx),il1,il2,lextend=.true.)
!print*,'PENCIL: frgn_setup%xgrid,nxgrid_foreign,x(l1),x(l2),il1,il2=',frgn_setup%xgrid(1), &
!frgn_setup%xgrid(nxgrid_foreign),x(l1),x(l2),il1,il2
!!! GM: PROBABLE INCONSISTENCY IN THE FOLLOWING COMMAND
!!!          if (.not.lfirst_proc_x) il1=il1-1
!!!          if (.not.llast_proc_x) il2=il2+1
!
          frgn_setup%xind_rng(-1,:)=(/il1,il2/)     ! global x index range of rank iproc
!
!              frgn_setup%xind_rng(0,:)=frgn_setup%xind_rng(-1,:)
!              frgn_setup%xpeer_rng=(/0,0/)
!              call mpisend_int(frgn_setup%xind_rng(-1,:),2,peer+ncpus,iproc+tag_foreign,MPI_COMM_WORLD)

!
!  Ask all foreign processors in first x beam
!  about their share in frgn_setup%xind_rng. No share: receive [0,0].
!
          do px=0,frgn_setup%procnums(1)-1
!
!       Determine the peer using EULAG proc grid conventions.
!
            peer = find_proc_general(0,0,px,frgn_setup%procnums(3),frgn_setup%procnums(2),frgn_setup%procnums(1),.true.)
!print*,'PCASA x-send', iproc, peer+ncpus, frgn_setup%xind_rng(-1,:)
            call mpisendrecv_int(frgn_setup%xind_rng(-1,:),2,peer+ncpus,tag_foreign+iproc, &
                                 intbuf,peer+ncpus,tag_foreign+peer,MPI_COMM_WORLD)
            frgn_setup%xind_rng(px,:)=intbuf(1:2)
!
            if (frgn_setup%xpeer_rng(1)>=0) then    ! if start of peer range has already been detected
!
!  If px has no share, px-1 is last of peer range.
!
              if (frgn_setup%xind_rng(px,1)>frgn_setup%xind_rng(px,2).and.frgn_setup%ypeer_rng(2)<0) then
                frgn_setup%xpeer_rng(2)=px-1
              endif
            else                                   ! if start of peer range has not yet been detected
!
!  If px has share, it is first of peer range.
!
              if (frgn_setup%xind_rng(px,1)<=frgn_setup%xind_rng(px,2)) frgn_setup%xpeer_rng(1)=px
            endif
          enddo
          if (frgn_setup%xpeer_rng(2)<0) frgn_setup%xpeer_rng(2)=px-1

          if (frgn_setup%linterpol(1)) &
            call get_linterp_weights_1D(frgn_setup%xgrid,frgn_setup%xind_rng(-1,:),x,frgn_setup%xinds,frgn_setup%xweights)
!print*, 'iproc, frgn_setup%xinds=', iproc, frgn_setup%xinds   !, sum(frgn_setup%xweights,2)
        endif  ! lfirst_proc_yz
!
        call mpibcast_int_arr2(frgn_setup%xind_rng,(/frgn_setup%procnums(1)+1,2/),comm=MPI_COMM_YZPLANE)
        call mpibcast_int_arr(frgn_setup%xpeer_rng,2,comm=MPI_COMM_YZPLANE)
        lenx=frgn_setup%xind_rng(-1,2)-frgn_setup%xind_rng(-1,1)+1

        if (frgn_setup%linterpol(1)) then
          call mpibcast_int(frgn_setup%xinds,mx,comm=MPI_COMM_YZPLANE)
          call mpibcast_real(frgn_setup%xweights,(/mx,2/),comm=MPI_COMM_YZPLANE)
        endif

        call mpibcast_int_arr2(frgn_setup%yind_rng,(/frgn_setup%procnums(2)+1,2/),comm=MPI_COMM_XZPLANE)
        call mpibcast_int_arr(frgn_setup%ypeer_rng,2,comm=MPI_COMM_XZPLANE)
        leny=frgn_setup%yind_rng(-1,2)-frgn_setup%yind_rng(-1,1)+1

        if (frgn_setup%linterpol(2)) then
          call mpibcast_int(frgn_setup%yinds,my,comm=MPI_COMM_XZPLANE)
          call mpibcast_real(frgn_setup%yweights,(/my,2/),comm=MPI_COMM_XZPLANE)
        endif

        if (allocated(frgn_buffer)) deallocate(frgn_buffer)
        allocate(frgn_buffer(lenx,leny,mz,3))
!print*, 'PENCIL0B',iproc,peer, lbound(frgn_buffer),ubound(frgn_buffer)
!print*, 'PENCIL global len*',iproc,peer, lenx, leny, mz

        if (allocated(frgn_setup%recv_req)) deallocate(frgn_setup%recv_req)
        allocate(frgn_setup%recv_req(0:frgn_setup%procnums(1)-1))

      endif    ! if (lforeign)

call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
print*, 'Pencil: after barrier'
!call MPI_FINALIZE(mpierr)
!stop

    endsubroutine initialize_foreign_comm
!***********************************************************************
    subroutine get_foreign_snap_initiate(nvars,frgn_buffer,lnonblock)
!
! Initializes fetching of data snapshot from foreign code.
! Non-blocking mode not yet operational.
!
! 20-oct-21/MR: coded
!
      use General, only: loptest, find_proc_general

      real, dimension(:,:,:,:) :: frgn_buffer
      integer :: nvars
      logical, optional :: lnonblock

      integer :: istart,ixstart,iystart,lenx_loc,leny_loc,px,py,iv,peer
      logical :: lesav !Not needed in the normal run as savings would come normally now. Retained for now, in case is needed again.

!GM: Assume will not save initially
      lesav = .false.

      do px=frgn_setup%xpeer_rng(1),frgn_setup%xpeer_rng(2)

        ixstart=frgn_setup%xind_rng(px,1)-frgn_setup%xind_rng(-1,1)+1!-nghost
        lenx_loc=frgn_setup%xind_rng(px,2)-frgn_setup%xind_rng(px,1)+1!+2*nghost

        do py=frgn_setup%ypeer_rng(1),frgn_setup%ypeer_rng(2)
          peer = find_proc_general(ipz/frgn_setup%proc_multis(3),py,px,frgn_setup%procnums(3), &
                                 frgn_setup%procnums(2),frgn_setup%procnums(1),.true.)

          iystart=frgn_setup%yind_rng(py,1)-frgn_setup%yind_rng(-1,1)+1!-nghost
          leny_loc=frgn_setup%yind_rng(py,2)-frgn_setup%yind_rng(py,1)+1!+2*nghost

          do iv=1,nvars
            if (loptest(lnonblock)) then      ! non-blocking
              call mpirecv_real(frgn_buffer(ixstart:ixstart+lenx_loc-1,iystart:iystart+leny_loc-1,:,iv), &
                               (/lenx_loc,leny_loc,mz/),peer+ncpus,iproc+tag_foreign,MPI_COMM_WORLD,frgn_setup%recv_req(px))
            else                              ! blocking
!print*, 'PENCIL1X', iproc,peer,px,ixstart, ixstart+lenx_loc-1
!print*, 'PENCIL1Y', iproc,peer,py,iystart, iystart+leny_loc-1
!print*, 'PENCIL3', iproc,peer,lenx_loc,leny_loc,mz
!print*,'PENCIL4', iproc, peer,l1, l2,m1,m2
!print*, 'PENCIL1B', iproc, lbound(frgn_buffer(ixstart:ixstart+lenx_loc-1,iystart:iystart+leny_loc-1,:,iv)),ubound(frgn_buffer(ixstart:ixstart+lenx_loc-1,iystart:iystart+leny_loc-1,:,iv))
!print*, 'PENCIL2B', iproc, size(x),size(y),size(z)
              call mpirecv_real(frgn_buffer(ixstart:ixstart+lenx_loc-1,iystart:iystart+leny_loc-1,:,iv), &
                               (/lenx_loc,leny_loc,mz/),peer+ncpus,iproc+tag_foreign,MPI_COMM_WORLD)
            endif
!
          enddo
        enddo
      enddo
!GM: Wait for the saving command
if (iproc==0) then
!  print*, 'PENCIL Before RECEIVED LOGICAL SAVE', lesav
  call mpirecv_logical(lesav,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
!  if (lesav) then
!      call mpirecv_real(time)
!      print*, 'Pencil Time', t, 'Eulag time', time_e
!  endif
endif
!GM: Broadcast saving logical
call mpibcast_logical(lesav,comm=MPI_COMM_PENCIL)
!if (lesav) call write_2daverages() !Not needed for now.
!if (lesav) call output_snap(,nv1=nv1_capitalvar,nv2=m     snap,file=file)
!print *, 'PENCIL - MINMAX W',iproc, minval(frgn_buffer(:,:,:,1)), maxval(frgn_buffer(:,:,:,1))
!print *, 'PENCIL - MINMAX V',iproc, minval(frgn_buffer(:,:,:,2)), maxval(frgn_buffer(:,:,:,2))
!print *, 'PENCIL - MINMAX U',iproc, minval(frgn_buffer(:,:,:,3)), maxval(frgn_buffer(:,:,:,3))

      return !END new version
!
print*, 'PENCIL-BARRIER', iproc
call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
call MPI_FINALIZE(mpierr)
stop
!
! Below old version.
!
      do px=0,frgn_setup%procnums(1)-1
        if (frgn_setup%xind_rng(px,1)>0) then

          istart=frgn_setup%xind_rng(px,1)-frgn_setup%xind_rng(-1,1)+1-nghost
          lenx_loc=frgn_setup%xind_rng(px,2)-frgn_setup%xind_rng(px,1)+1+2*nghost

          if (loptest(lnonblock)) then
            do iv=1,nvars
              call mpirecv_real(frgn_buffer(istart:istart+lenx_loc-1,:,:,iv), &
                                (/lenx_loc,my,mz/),peer,peer-ncpus,MPI_COMM_WORLD,frgn_setup%recv_req(px))
            enddo
          else       ! blocking case
            !peer=frgn_setup%peer_rng(1)
!print*, 'PENCIL recv: iproc,peer,tag=', iproc,peer+ncpus,peer
            do iv=1,nvars
              call mpirecv_real(frgn_buffer(istart:istart+lenx_loc-1,:,:,iv), &
                                (/lenx_loc,my,mz/),peer+ncpus,peer,MPI_COMM_WORLD)
!              if (iv.eq.1) write(100+iproc) frgn_buffer(:,:,:,iv)
!              if (iv.eq.2) write(400+iproc) frgn_buffer(:,:,:,iv)
!              if (iv.eq.3) write(700+iproc) frgn_buffer(:,:,:,iv)
!print *,'PENCIL MINMAX FRGN',iproc, minval(frgn_buffer), maxval(frgn_buffer)
            enddo
!
! If in time loop, receive signal from foreign code to save snapshot for
! synchronous restart capability.
!
            if (it>0) call mpibcast_logical(lsave,frgn_setup%root,MPI_COMM_WORLD)
          endif
        endif
      enddo
!write(200+iproc) frgn_buffer
!if (.not.loptest(lnonblock))then
!else
!  print*, 'PENCIL INIT NON-BLOCK'
!endif
!print*, 'PENCIL get_foreign_snap_initiate: successful', iproc
!print *, 'PENCIL - MIN MAX W',iproc, minval(frgn_buffer(:,:,:,1)), maxval(frgn_buffer(:,:,:,1))
!print *, 'PENCIL - MIN MAX V',iproc, minval(frgn_buffer(:,:,:,2)), maxval(frgn_buffer(:,:,:,2))
!print *, 'PENCIL - MIN MAX U',iproc, minval(frgn_buffer(:,:,:,3)), maxval(frgn_buffer(:,:,:,3))

!call mpibarrier(MPI_COMM_WORLD)
!call MPI_FINALIZE(mpierr)
!stop

    endsubroutine get_foreign_snap_initiate
!***********************************************************************
    subroutine get_foreign_snap_finalize(f,ivar1,ivar2,frgn_buffer,interp_buffer,lnonblock)
!
! Finalizes fetching of data snapshot from foreign code: when nonblocking, call MPI_WAIT;
! interpolate data and scatter across processors.
! Non-blocking mode not yet operational.
!
! 20-oct-21/MR: coded
!
      use General, only: loptest,linear_interpolate_1d, itoa
!      use IO, only: output_snap, output_snap_finalize

      real, dimension(:,:,:,:) :: f,frgn_buffer,interp_buffer
      integer :: ivar1, ivar2, lf1
      logical, optional :: lnonblock
      integer, parameter :: ytag=115
      integer :: istart_y,istart_z,istop_y,istop_z,px,py,pz,partner,ll,nvars
!      logical :: lesav
!      integer, save :: sav_counter = 0
!      character(*), parameter :: filepath = "data/allprocs/"

      !if (mod(frgn_setup%dims(2),nygrid)==0 .and. mod(frgn_setup%dims(3),nzgrid)==0) then
      lf1 = 1
!
! Interpolate/scatter data to array f
!
      if (trim(frgn_setup%name)=='EULAG') then

        if (loptest(lnonblock)) call mpiwait(frgn_setup%recv_req(px))
        if (size(frgn_buffer,1) > size(f,1)) lf1 = 2
!!!        f(:,:,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU  !ORIGINAL SETUP
        if (lfirst_proc_y .and. nprocy==1) then
          if (lfirst_proc_x .and. nprocx==1) then!if only one proc in X
            f(l1:l2,m1:m2,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
!
          else if (.not.(lfirst_proc_x .or. llast_proc_x))then!MID of XBEAM
            f(:,m1:m2,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
!
          else if (lfirst_proc_x) then! on processors of first XBEAM
!
            f(l1:,m1:m2,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
          else if (llast_proc_x) then!on processors of last XBEAM
            f(:l2,m1:m2,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
          endif

        else if (.not.(lfirst_proc_y .or. llast_proc_y)) then!MID of YBEAM
!
          if (lfirst_proc_x .and. nprocx==1) then!if only one proc in X
            f(l1:l2,:,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
!
          else if (.not.(lfirst_proc_x .or. llast_proc_x))then!MID of XBEAM
            f(:,:,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
!
          else if (lfirst_proc_x) then! on processors of first XBEAM
!
            f(l1:,:,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
          else if (llast_proc_x) then!on processors of last XBEAM
            f(:l2,:,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
          endif

!       i
        else if (lfirst_proc_y) then! on processors of first YBEAM
          if (lfirst_proc_x .and. nprocx==1) then !if only one proc in X
            f(l1:l2,m1:,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
!
          else if (.not.(lfirst_proc_x .or. llast_proc_x))then!MID of XBEAM
            f(:,m1:,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
!
          else if (lfirst_proc_x) then! on processors of first XBEAM
!
            f(l1:,m1:,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
          else if (llast_proc_x) then!on processors of last XBEAM
            f(:l2,m1:,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
          endif
!
        else if (llast_proc_y) then!on processors of last YBEAM
          if (lfirst_proc_x .and. nprocx==1) then!if only one proc in X
            f(l1:l2,:m2,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
!
          else if (.not.(lfirst_proc_x .or. llast_proc_x))then!MID of XBEAM
            f(:,:m2,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
          else if (lfirst_proc_x) then! on processors of first XBEAM
!
            f(l1:,:m2,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
          else if (llast_proc_x) then!on processors of last XBEAM
            f(:l2,:m2,:,ivar1:ivar2)=frgn_buffer(lf1:,:,:,:)/frgn_setup%renorm_UU
          endif
        endif
!
!print*,'PENCILFin0', iproc, m1,m2
!print*,'PENCILFin1', iproc, lbound(f),ubound(f)
!print*,'PENCILFin2', iproc, lbound(frgn_buffer),ubound(frgn_buffer)
!print*,'PENCILFin3', iproc, ncpus, nghost+1, mx-nghost, my-nghost, mz-nghost
!print*,'PENCILFin4', iproc, size(frgn_buffer,1) , size(frgn_buffer,2)

!GM: Wait for the saving command/ IF SAVING HERE IS NEEDED
!if (iproc==0) then
!  print*, 'PENCIL Before RECEIVED LOGICAL SAVE', lesav
!  call mpirecv_logical(lesav,frgn_setup%root,tag_foreign,MPI_COMM_WORLD)
!  print*, 'PENCIL RECEIVED LOGICAL SAVE', lesav
!  print*, 'PENCIL TEST savefile, itoa', 'data'//itoa(sav_counter)
!  sav_counter = sav_counter + 1
!  if (lesav) then
    !call mpirecv_real(time)
    !print*, 'Pencil Time', t, 'Eulag time', time_e
!  endif
!endif
!GM: Broadcast saving logical
!call mpibcast_logical(lesav,comm=MPI_COMM_PENCIL)

!if (lesav) then !SAVES FULL DATA
!  call output_snap(f,nv1=1,nv2=mfarray,file='data'//itoa(sav_counter)//'.dat',mode=1)
!  call output_snap_finalize
!endif
goto 125!!!
        !!!f(:,:,:,ivar1:ivar2)=frgn_buffer(lf1:,my::-1,:,:)/frgn_setup%renorm_UU  !For invertion of theta
!
!print*, 'iproc, weights=', iproc, maxval(frgn_setup%xweights), maxval(frgn_setup%yweights)
        if (all(frgn_setup%linterpol)) then
! XYZ-interpolation
        elseif (all(frgn_setup%linterpol(1:2))) then
!
! XY-Interpolation
!
          call interpolate_2d(frgn_buffer(lf1:,:,:,:),f(:,:,:,ivar1:ivar2),frgn_setup%xinds,frgn_setup%yinds, &
                              frgn_setup%xweights,frgn_setup%yweights,1./frgn_setup%renorm_UU)
        elseif (all(frgn_setup%linterpol(2:3))) then
! YZ-Interpolation
        elseif (all(frgn_setup%linterpol((/1,3/)))) then
! XZ-Interpolation
        elseif (frgn_setup%linterpol(1)) then
!
! X-Interpolation
!
          call interpolate_1d(frgn_buffer(lf1:,:,:,:),f(:,:,:,ivar1:ivar2),1,frgn_setup%xinds, &
                              frgn_setup%xweights,1./frgn_setup%renorm_UU)
        elseif (frgn_setup%linterpol(2)) then
!
! Y-Interpolation
!
          call interpolate_1d(frgn_buffer(lf1:,:,:,:),f(:,:,:,ivar1:ivar2),2,frgn_setup%yinds, &
                              frgn_setup%yweights,1./frgn_setup%renorm_UU)
        elseif (frgn_setup%linterpol(3)) then
! Z-Interpolation
        endif
!
125 continue !!!
      elseif (lfirst_proc_yz) then

        if (frgn_setup%proc_multis(1)>1) then
          nvars=ivar2-ivar1+1

          do pz=0,nprocz-1
            istart_z=pz*nz+1; istop_z=istart_z+nz-1

            do py=0,nprocy-1

              if (lprocz_slowest) then
                partner = py + nprocy*pz
              else
                partner = pz + nprocz*py
              endif
              istart_y=py*ny+1; istop_y=istart_y+ny-1

              do ll=l1,l2
                call linear_interpolate_1d(frgn_buffer(:,istart_y:istop_y,istart_z:istop_z,:), &
                                           frgn_setup%xgrid(frgn_setup%xind_rng(-1,1):frgn_setup%xind_rng(-1,2)), &
                                           x(ll),interp_buffer(ll,:,:,:),.true.)
              enddo
              interp_buffer=interp_buffer/frgn_setup%renorm_UU
              call mpisend_real(interp_buffer,(/nx,ny,nz,nvars/),partner,ytag,MPI_COMM_YZPLANE)

              if (size(f,1)==mx) then     ! f has ghosts
                f(l1:l2,m1:m2,n1:n2,ivar1:ivar2)=interp_buffer
              else
                f(:,:,:,ivar1:ivar2)=interp_buffer
              endif
            enddo
          enddo
        endif
          !call mpirecv_real(f(l1:l2,m1:m2,n1:n2,ivar1:ivar2),(/nx,ny,nz,nvars/),0,ytag,MPI_COMM_YZPLANE)
      endif     ! if (name==EULAG)

    endsubroutine get_foreign_snap_finalize
!***********************************************************************
    subroutine interpolate_2d(inbuffer,outbuffer,xinds,yinds,xweights,yweights,scal_)
!
! At the moment only for first two dimensions.
!
      use General, only: roptest

      real, dimension(:,:,:,:) :: inbuffer,outbuffer
      !real, dimension(mx,my,mz,mvar) :: outbuffer
      integer, dimension(:) :: xinds,yinds
      real, dimension(:,:) :: xweights,yweights
      real, optional :: scal_

      real :: yw1,yw2,scal
      integer :: ll,mm,xind1,xind2,yind1,yind2

      scal=roptest(scal_,1.)

      do mm=1,size(outbuffer,2)

        yw1=yweights(mm,1)*scal; yw2=yweights(mm,2)*scal
        yind1=yinds(mm); yind2=yind1+1

        do ll=1,size(outbuffer,1)
          xind1=xinds(ll); xind2=xind1+1
          outbuffer(ll,mm,:,:) =  (xweights(ll,1)*yw1)*inbuffer(xind1,yind1,:,:) &
                                + (xweights(ll,2)*yw1)*inbuffer(xind2,yind1,:,:) &
                                + (xweights(ll,1)*yw2)*inbuffer(xind1,yind2,:,:) &
                                + (xweights(ll,2)*yw2)*inbuffer(xind2,yind2,:,:)
        enddo
      enddo

    endsubroutine interpolate_2d
!***********************************************************************
    subroutine interpolate_1d(inbuffer,outbuffer,dim,inds,weights,scal_)
!
! At the moment only for first two dimensions.
!
      use General, only: roptest

      real, dimension(:,:,:,:) :: inbuffer,outbuffer
      !real, dimension(mx,my,mz,mvar) :: outbuffer
      integer :: dim
      integer, dimension(:) :: inds
      real, dimension(:,:) :: weights
      real, optional :: scal_

      real :: w1,w2,scal
      integer :: ii,ind1,ind2

      scal=roptest(scal_,1.)

      do ii=1,size(outbuffer,dim)

        w1=weights(ii,1)*scal; w2=weights(ii,2)*scal
        ind1=inds(ii); ind2=ind1+1
        select case (dim)
        case (1); outbuffer(ii,:,:,:) = w1*inbuffer(ind1,:,:,:) + w2*inbuffer(ind2,:,:,:)
        case (2); outbuffer(:,ii,:,:) = w1*inbuffer(:,ind1,:,:) + w2*inbuffer(:,ind2,:,:)
        case (3); outbuffer(:,:,ii,:) = w1*inbuffer(:,:,ind1,:) + w2*inbuffer(:,:,ind2,:)
        end select

      enddo

    endsubroutine interpolate_1d
!***********************************************************************
    logical function update_foreign_data(t,t_foreign)
!
! Determines whether it is time to fetch new data from foreign code.
! Returns delivery cadence of foreign code in dt_foreign.
!
! 20-oct-21/MR: coded
!
      real(KIND=rkind8), intent(IN) :: t
      real,              intent(OUT):: t_foreign

      if (t-frgn_setup%t_last_recvd>0.) then
        frgn_setup%t_last_recvd=frgn_setup%t_last_recvd+frgn_setup%dt_out
        update_foreign_data=.true.
        t_foreign = frgn_setup%t_last_recvd
      else
        update_foreign_data=.false.
      endif

    endfunction update_foreign_data
!***********************************************************************
    subroutine set_rslice_communicator

      integer :: rank

      call MPI_COMM_SPLIT(MPI_COMM_PENCIL, lwrite_slice_r, iproc, MPI_COMM_RSLICE, mpierr)
      call MPI_COMM_RANK(MPI_COMM_RSLICE,rank,mpierr)
      if (rank==0) root_rslice=iproc

    endsubroutine set_rslice_communicator
!***********************************************************************
    subroutine mpiallreduce_merge(list,len,comm)
!
!  Merges an integer list across all processors.
!  Returns list of unique elements.
!
!  19-oct-22/MR: coded
!
      use General, only: merge_lists

      integer, dimension(*) :: list
      integer :: len
      integer, optional, intent(in) :: comm

      integer :: maxlen,ip,len_proc
      integer, dimension(:), allocatable :: buffer

      call mpireduce_max(len,maxlen,comm)
      if (lroot) then
        allocate(buffer(maxlen))

        do ip=1,ncpus-1
          call mpirecv_int(len_proc,ip,ip,comm)
          call mpirecv_int(buffer,len_proc,ip,ip,comm)
          call merge_lists(list,len,buffer(1:len_proc))
        enddo

      else
        call mpisend_int(len,root,iproc,comm)
        call mpisend_int(list,len,root,iproc,comm)
      endif

      call mpibcast_int(len,root,comm)
      call mpibcast_int(list,len,root,comm)

    endsubroutine mpiallreduce_merge
!***********************************************************************
  endmodule Mpicomm
