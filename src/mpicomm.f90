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
module Mpicomm
!
  use Cdata
  use Cparam
!
  implicit none
!
  include 'mpicomm.h'
!
  interface mpirecv_logical
     module procedure mpirecv_logical_scl
     module procedure mpirecv_logical_arr
  endinterface
!
  interface mpirecv_real
    module procedure mpirecv_real_scl
    module procedure mpirecv_real_arr
    module procedure mpirecv_real_arr2
    module procedure mpirecv_real_arr3
    module procedure mpirecv_real_arr4
  endinterface
!
  interface mpirecv_int
    module procedure mpirecv_int_scl
    module procedure mpirecv_int_arr
    module procedure mpirecv_int_arr2
  endinterface
!
  interface mpisend_logical
     module procedure mpisend_logical_scl
     module procedure mpisend_logical_arr
  endinterface
!
  interface mpisend_real
    module procedure mpisend_real_scl
    module procedure mpisend_real_arr
    module procedure mpisend_real_arr2
    module procedure mpisend_real_arr3
    module procedure mpisend_real_arr4
  endinterface
!
  interface mpisend_int
    module procedure mpisend_int_scl
    module procedure mpisend_int_arr
    module procedure mpisend_int_arr2
  endinterface
!
  interface mpibcast_logical
    module procedure mpibcast_logical_scl
    module procedure mpibcast_logical_arr
    module procedure mpibcast_logical_arr2
  endinterface
!
  interface mpibcast_int
    module procedure mpibcast_int_scl
    module procedure mpibcast_int_arr
  endinterface
!
  interface mpibcast_real
    module procedure mpibcast_real_scl
    module procedure mpibcast_real_arr
    module procedure mpibcast_real_arr2
    module procedure mpibcast_real_arr3
    module procedure mpibcast_real_arr4
  endinterface
!
  interface mpibcast_double
    module procedure mpibcast_double_scl
    module procedure mpibcast_double_arr
  endinterface
!
  interface mpibcast_char
    module procedure mpibcast_char_scl
    module procedure mpibcast_char_arr
  endinterface
!
  interface mpiallreduce_sum
    module procedure mpiallreduce_sum_scl
    module procedure mpiallreduce_sum_arr
    module procedure mpiallreduce_sum_arr2
    module procedure mpiallreduce_sum_arr3
  endinterface
!
  interface mpiallreduce_sum_int
     module procedure mpiallreduce_sum_int_scl
     module procedure mpiallreduce_sum_int_arr
  endinterface
!
  interface mpiallreduce_max
    module procedure mpiallreduce_max_scl
    module procedure mpiallreduce_max_arr
  endinterface
!
  interface mpireduce_max
    module procedure mpireduce_max_scl
    module procedure mpireduce_max_arr
  endinterface
!
  interface mpireduce_max_int
    module procedure mpireduce_max_scl_int
  endinterface
!
  interface mpireduce_min
    module procedure mpireduce_min_scl
    module procedure mpireduce_min_arr
  endinterface
!
  interface mpireduce_or
    module procedure mpireduce_or_scl
    module procedure mpireduce_or_arr
  endinterface
!
  interface mpireduce_and
    module procedure mpireduce_and_scl
    module procedure mpireduce_and_arr
  endinterface
!
  interface mpireduce_sum_int
    module procedure mpireduce_sum_int_scl
    module procedure mpireduce_sum_int_arr
    module procedure mpireduce_sum_int_arr2
    module procedure mpireduce_sum_int_arr3
    module procedure mpireduce_sum_int_arr4
  endinterface
!
  interface mpireduce_sum
    module procedure mpireduce_sum_scl
    module procedure mpireduce_sum_arr
    module procedure mpireduce_sum_arr2
    module procedure mpireduce_sum_arr3
    module procedure mpireduce_sum_arr4
  endinterface
!
  interface mpireduce_sum_double
    module procedure mpireduce_sum_double_scl
    module procedure mpireduce_sum_double_arr
    module procedure mpireduce_sum_double_arr2
    module procedure mpireduce_sum_double_arr3
    module procedure mpireduce_sum_double_arr4
  endinterface
!
  interface distribute_xy
    module procedure distribute_xy_0D
    module procedure distribute_xy_2D
    module procedure distribute_xy_3D
    module procedure distribute_xy_4D
  endinterface
!
  interface collect_xy
    module procedure collect_xy_0D
    module procedure collect_xy_2D
    module procedure collect_xy_3D
    module procedure collect_xy_4D
  endinterface
!
  interface distribute_z
    module procedure distribute_z_3D
    module procedure distribute_z_4D
  endinterface
!
  interface collect_z
    module procedure collect_z_3D
    module procedure collect_z_4D
  endinterface
!
  interface distribute_to_pencil_xy
    module procedure distribute_to_pencil_xy_2D
  endinterface
!
  interface collect_from_pencil_xy
    module procedure collect_from_pencil_xy_2D
  endinterface
!
  interface remap_to_pencil_y
    module procedure remap_to_pencil_y_1D
    module procedure remap_to_pencil_y_2D
    module procedure remap_to_pencil_y_3D
    module procedure remap_to_pencil_y_4D
  endinterface
!
  interface unmap_from_pencil_y
    module procedure unmap_from_pencil_y_1D
    module procedure unmap_from_pencil_y_2D
    module procedure unmap_from_pencil_y_3D
    module procedure unmap_from_pencil_y_4D
  endinterface
!
  interface remap_to_pencil_z
    module procedure remap_to_pencil_z_1D
    module procedure remap_to_pencil_z_2D
    module procedure remap_to_pencil_z_3D
    module procedure remap_to_pencil_z_4D
  endinterface
!
  interface unmap_from_pencil_z
    module procedure unmap_from_pencil_z_1D
    module procedure unmap_from_pencil_z_2D
    module procedure unmap_from_pencil_z_3D
    module procedure unmap_from_pencil_z_4D
  endinterface
!
  interface remap_to_pencil_xy
    module procedure remap_to_pencil_xy_2D
    module procedure remap_to_pencil_xy_3D
    module procedure remap_to_pencil_xy_4D
  endinterface
!
  interface unmap_from_pencil_xy
    module procedure unmap_from_pencil_xy_2D
    module procedure unmap_from_pencil_xy_3D
    module procedure unmap_from_pencil_xy_4D
  endinterface
!
  interface transp_pencil_xy
    module procedure transp_pencil_xy_2D
    module procedure transp_pencil_xy_3D
    module procedure transp_pencil_xy_4D
  endinterface
!
  interface remap_to_pencil_yz
    module procedure remap_to_pencil_yz_3D
    module procedure remap_to_pencil_yz_4D
  endinterface
!
  interface unmap_from_pencil_yz
    module procedure unmap_from_pencil_yz_3D
    module procedure unmap_from_pencil_yz_4D
  endinterface
!
!  interface mpigather_and_out
!    module procedure mpigather_and_out_real
!    module procedure mpigather_and_out_cmplx
!  endinterface
!
  include 'mpif.h'
!
!  initialize debug parameter for this routine
!
  logical :: ldebug_mpi=.false.
  integer :: mpi_precision
!
!  For f-array processor boundaries
!
  real, dimension (nghost,ny,nz,mcom) :: lbufxi,ubufxi,lbufxo,ubufxo
  real, dimension (nx,nghost,nz,mcom) :: npbufyi,npbufyo,spbufyi,spbufyo
  real, dimension (mx,nghost,nz,mcom) :: lbufyi,ubufyi,lbufyo,ubufyo
  real, dimension (mx,ny,nghost,mcom) :: lbufzi,ubufzi,lbufzo,ubufzo
  real, dimension (mx,nghost,nghost,mcom) :: llbufi,lubufi,uubufi,ulbufi
  real, dimension (mx,nghost,nghost,mcom) :: llbufo,lubufo,uubufo,ulbufo
!
  real, dimension (nghost,my,mz,mcom) :: fahi,falo,fbhi,fblo,fao,fbo ! For shear
  integer :: ipx_partner, nextya, nextyb, lastya, lastyb, displs ! For shear
  integer :: nprocs, mpierr
  integer :: serial_level = 0
!
!  mpi tags
!
  integer :: tolowx=13,touppx=14,tolowy=3,touppy=4,tolowz=5,touppz=6 ! msg. tags
  integer :: TOll=7,TOul=8,TOuu=9,TOlu=10 ! msg. tags for corners
  integer :: io_perm=20,io_succ=21
  integer :: npole_tag=15,spole_tag=16
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
  integer :: MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM
  integer :: MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE
!
  integer :: isend_rq_tolowx,isend_rq_touppx,irecv_rq_fromlowx,irecv_rq_fromuppx
  integer :: isend_rq_spole,isend_rq_npole
  integer :: irecv_rq_spole,irecv_rq_npole
  integer :: isend_rq_tolowy,isend_rq_touppy,irecv_rq_fromlowy,irecv_rq_fromuppy
  integer :: isend_rq_tolowz,isend_rq_touppz,irecv_rq_fromlowz,irecv_rq_fromuppz
  integer :: isend_rq_TOll,isend_rq_TOul,isend_rq_TOuu,isend_rq_TOlu  !(corners)
  integer :: irecv_rq_FRuu,irecv_rq_FRlu,irecv_rq_FRll,irecv_rq_FRul  !(corners)
  integer :: isend_rq_tolastya,isend_rq_tonextya, &
             irecv_rq_fromlastya,irecv_rq_fromnextya ! For shear
  integer :: isend_rq_tolastyb,isend_rq_tonextyb, &
             irecv_rq_fromlastyb,irecv_rq_fromnextyb ! For shear
!
  integer, dimension (MPI_STATUS_SIZE) :: isend_stat_tl,isend_stat_tu
  integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fl,irecv_stat_fu
  integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_np,irecv_stat_sp,&
                                          isend_stat_np,isend_stat_sp
  integer, dimension (MPI_STATUS_SIZE) :: isend_stat_Tll,isend_stat_Tul, &
                                          isend_stat_Tuu,isend_stat_Tlu
  integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_Fuu,irecv_stat_Flu, &
                                          irecv_stat_Fll,irecv_stat_Ful
  integer, dimension (MPI_STATUS_SIZE) :: isend_stat_spole,irecv_stat_spole, &
                                          isend_stat_npole,irecv_stat_npole
!
  contains
!
!***********************************************************************
    subroutine mpicomm_init()
!
!  Get processor number, number of procs, and whether we are root.
!
!  20-aug-01/wolf: coded
!  29-jul-2010/anders: separate subroutine
!
      use Syscalls, only: sizeof_real
!
      lmpicomm = .true.
      call MPI_INIT(mpierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, iproc , mpierr)
      lroot = (iproc==root)
!
      if (sizeof_real() < 8) then
        mpi_precision = MPI_REAL
        if (lroot) write (*,*) "\nWARNING: using SINGLE PRECISION, you'd better know what you do!\n"
      else
        mpi_precision = MPI_DOUBLE_PRECISION
      endif
!
    endsubroutine mpicomm_init
!***********************************************************************
    subroutine initialize_mpicomm()
!
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
!
!  Check consistency in processor layout.
!
      if (ncpus/=nprocx*nprocy*nprocz) then
        if (lroot) then
          print*, 'Compiled with ncpus = ', ncpus, &
              ', but nprocx*nprocy*nprocz=', nprocx*nprocy*nprocz
        endif
        call stop_it('initialize_mpicomm')
      endif
!
!  Check total number of processors.
!
      if (nprocs/=nprocx*nprocy*nprocz) then
        if (lroot) then
          print*, 'Compiled with ncpus = ', ncpus, &
              ', but running on ', nprocs, ' processors'
        endif
        call stop_it('initialize_mpicomm')
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
!  Position on the processor grid.
!  x is fastest direction, z slowest (this is the default)
!
      if (lprocz_slowest) then
        ipx = modulo(iproc, nprocx)
        ipy = modulo(iproc/nprocx, nprocy)
        ipz = iproc/nprocxy
      else
        ipx = modulo(iproc, nprocx)
        ipy = iproc/nprocxy
        ipz = modulo(iproc/nprocx, nprocy)
      endif
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
!  Set up `lower' and `upper' neighbours.
!
      xlneigh = modulo(ipx-1,nprocx) + ipy*nprocx + ipz*nprocxy
      xuneigh = modulo(ipx+1,nprocx) + ipy*nprocx + ipz*nprocxy
      ylneigh = ipx + modulo(ipy-1,nprocy)*nprocx + ipz*nprocxy
      yuneigh = ipx + modulo(ipy+1,nprocy)*nprocx + ipz*nprocxy
      zlneigh = ipx + ipy*nprocx + modulo(ipz-1,nprocz)*nprocxy
      zuneigh = ipx + ipy*nprocx + modulo(ipz+1,nprocz)*nprocxy
!
! For boundary condition across the pole set up pole-neighbours
! This assumes that the domain is equally distributed among the
! processors in the z direction.
!
      poleneigh = modulo(ipz+nprocz/2,nprocz)*nprocxy+ipy*nprocx+ipx
!
!  Set the four corners in the yz-plane (in cyclic order).
!
      llcorn = ipx + modulo(ipy-1,nprocy)*nprocx + modulo(ipz-1,nprocz)*nprocxy
      ulcorn = ipx + modulo(ipy+1,nprocy)*nprocx + modulo(ipz-1,nprocz)*nprocxy
      uucorn = ipx + modulo(ipy+1,nprocy)*nprocx + modulo(ipz+1,nprocz)*nprocxy
      lucorn = ipx + modulo(ipy-1,nprocy)*nprocx + modulo(ipz+1,nprocz)*nprocxy
!
!  This value is not yet the one read in, but the one initialized in cparam.f90.
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
!  NOTE: the ip print parameter has *not* yet been read at this point.
!  Therefore it must be invoked by resetting ldebug_mpi appropriately.
!
      if (ldebug_mpi) write(6,'(A,I4,"(",3I4,"): ",8I4)') &
        'initialize_mpicomm: MPICOMM neighbors ', &
        iproc,ipx,ipy,ipz, &
        ylneigh,llcorn,zlneigh,ulcorn,yuneigh,uucorn,zuneigh,lucorn
!
!  Define MPI communicators that include all processes sharing the same value
!  of ipx, ipy, or ipz. The rank within MPI_COMM_WORLD is given by a
!  combination of the two other directional processor indices.
!
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, ipy+nprocy*ipz, ipx, &
          MPI_COMM_XBEAM, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, ipx+nprocx*ipz, ipy, &
          MPI_COMM_YBEAM, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, ipx+nprocx*ipy, ipz, &
          MPI_COMM_ZBEAM, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, ipz, ipx+nprocx*ipy, &
          MPI_COMM_XYPLANE, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, ipy, ipx+nprocx*ipz, &
          MPI_COMM_XZPLANE, mpierr)
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, ipx, ipy+nprocy*ipz, &
          MPI_COMM_YZPLANE, mpierr)
!
    endsubroutine initialize_mpicomm
!***********************************************************************
    subroutine initiate_isendrcv_bdry(f,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values. Called in the beginning of pde.
!  Does not wait for the receives to finish (done in finalize_isendrcv_bdry)
!  leftneigh and rghtneigh are initialized by mpicomm_init.
!
!  21-may-02/axel: communication of corners added
!  11-aug-07/axel: communication in the x-direction added
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbufy, nbufz, nbufyz
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
      if (ivar2==0) return
!
!  Periodic boundary conditions in x.
!
      if (nprocx>1) call isendrcv_bdry_x(f,ivar1_opt,ivar2_opt)
!
!  Periodic boundary conditions in y.
!
      if (nprocy>1) then
        lbufyo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n1:n2,ivar1:ivar2) !!(lower y-zone)
        ubufyo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n1:n2,ivar1:ivar2) !!(upper y-zone)
        nbufy=mx*nz*nghost*(ivar2-ivar1+1)
        call MPI_IRECV(ubufyi(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
            yuneigh,tolowy,MPI_COMM_WORLD,irecv_rq_fromuppy,mpierr)
        call MPI_IRECV(lbufyi(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
            ylneigh,touppy,MPI_COMM_WORLD,irecv_rq_fromlowy,mpierr)
        call MPI_ISEND(lbufyo(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
            ylneigh,tolowy,MPI_COMM_WORLD,isend_rq_tolowy,mpierr)
        call MPI_ISEND(ubufyo(:,:,:,ivar1:ivar2),nbufy,MPI_REAL, &
            yuneigh,touppy,MPI_COMM_WORLD,isend_rq_touppy,mpierr)
      endif
!
!  Periodic boundary conditions in z.
!
      if (nprocz>1) then
        lbufzo(:,:,:,ivar1:ivar2)=f(:,m1:m2,n1:n1i,ivar1:ivar2) !!(lower z-zone)
        ubufzo(:,:,:,ivar1:ivar2)=f(:,m1:m2,n2i:n2,ivar1:ivar2) !!(upper z-zone)
        nbufz=mx*ny*nghost*(ivar2-ivar1+1)
        call MPI_IRECV(ubufzi(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
            zuneigh,tolowz,MPI_COMM_WORLD,irecv_rq_fromuppz,mpierr)
        call MPI_IRECV(lbufzi(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
            zlneigh,touppz,MPI_COMM_WORLD,irecv_rq_fromlowz,mpierr)
        call MPI_ISEND(lbufzo(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
            zlneigh,tolowz,MPI_COMM_WORLD,isend_rq_tolowz,mpierr)
        call MPI_ISEND(ubufzo(:,:,:,ivar1:ivar2),nbufz,MPI_REAL, &
            zuneigh,touppz,MPI_COMM_WORLD,isend_rq_touppz,mpierr)
        if (lnorth_pole) call isendrcv_bdry_npole(f,ivar1_opt,ivar2_opt)
        if (lsouth_pole) call isendrcv_bdry_spole(f,ivar1_opt,ivar2_opt)
      endif
!
!  The four corners (in counter-clockwise order).
!  (NOTE: this should work even for nprocx>1)
!
      if (nprocy>1.and.nprocz>1) then
        llbufo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n1:n1i,ivar1:ivar2)
        ulbufo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n1:n1i,ivar1:ivar2)
        uubufo(:,:,:,ivar1:ivar2)=f(:,m2i:m2,n2i:n2,ivar1:ivar2)
        lubufo(:,:,:,ivar1:ivar2)=f(:,m1:m1i,n2i:n2,ivar1:ivar2)
        nbufyz=mx*nghost*nghost*(ivar2-ivar1+1)
        call MPI_IRECV(uubufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            uucorn,TOll,MPI_COMM_WORLD,irecv_rq_FRuu,mpierr)
        call MPI_IRECV(lubufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            lucorn,TOul,MPI_COMM_WORLD,irecv_rq_FRlu,mpierr)
        call MPI_IRECV(llbufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            llcorn,TOuu,MPI_COMM_WORLD,irecv_rq_FRll,mpierr)
        call MPI_IRECV(ulbufi(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            ulcorn,TOlu,MPI_COMM_WORLD,irecv_rq_FRul,mpierr)
        call MPI_ISEND(llbufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            llcorn,TOll,MPI_COMM_WORLD,isend_rq_TOll,mpierr)
        call MPI_ISEND(ulbufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            ulcorn,TOul,MPI_COMM_WORLD,isend_rq_TOul,mpierr)
        call MPI_ISEND(uubufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            uucorn,TOuu,MPI_COMM_WORLD,isend_rq_TOuu,mpierr)
        call MPI_ISEND(lubufo(:,:,:,ivar1:ivar2),nbufyz,MPI_REAL, &
            lucorn,TOlu,MPI_COMM_WORLD,isend_rq_TOlu,mpierr)
      endif
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
      if (ivar2==0) return
!
!  1. wait until data received
!  2. set ghost zones
!  3. wait until send completed, will be overwritten in next time step
!
!  Communication in y (includes periodic bc)
!
      if (nprocy>1) then
        call MPI_WAIT(irecv_rq_fromuppy,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowy,irecv_stat_fl,mpierr)
        do j=ivar1,ivar2
          if (.not. lfirst_proc_y .or. bcy1(j)=='p') then
            f(:, 1:m1-1,n1:n2,j)=lbufyi(:,:,:,j)  !!(set lower buffer)
          endif
          if (.not. llast_proc_y .or. bcy2(j)=='p') then
            f(:,m2+1:my,n1:n2,j)=ubufyi(:,:,:,j)  !!(set upper buffer)
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
          if (.not. lfirst_proc_z .or. bcz1(j)=='p') then
            f(:,m1:m2, 1:n1-1,j)=lbufzi(:,:,:,j)  !!(set lower buffer)
          endif
          if (.not. llast_proc_z .or. bcz2(j)=='p') then
            f(:,m1:m2,n2+1:mz,j)=ubufzi(:,:,:,j)  !!(set upper buffer)
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowz,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppz,isend_stat_tu,mpierr)
      endif
!
!  The four yz-corners (in counter-clockwise order)
!
      if (nprocy>1.and.nprocz>1) then
        call MPI_WAIT(irecv_rq_FRuu,irecv_stat_Fuu,mpierr)
        call MPI_WAIT(irecv_rq_FRlu,irecv_stat_Flu,mpierr)
        call MPI_WAIT(irecv_rq_FRll,irecv_stat_Fll,mpierr)
        call MPI_WAIT(irecv_rq_FRul,irecv_stat_Ful,mpierr)
        do j=ivar1,ivar2
          if (.not. lfirst_proc_z .or. bcz1(j)=='p') then
            if (.not. lfirst_proc_y .or. bcy1(j)=='p') then
              f(:, 1:m1-1, 1:n1-1,j)=llbufi(:,:,:,j)  !!(set ll corner)
            endif
            if (.not. llast_proc_y .or. bcy2(j)=='p') then
              f(:,m2+1:my, 1:n1-1,j)=ulbufi(:,:,:,j)  !!(set ul corner)
            endif
          endif
          if (.not. llast_proc_z .or. bcz2(j)=='p') then
            if (.not. llast_proc_y .or. bcy2(j)=='p') then
              f(:,m2+1:my,n2+1:mz,j)=uubufi(:,:,:,j)  !!(set uu corner)
            endif
            if (.not. lfirst_proc_y .or. bcy1(j)=='p') then
              f(:, 1:m1-1,n2+1:mz,j)=lubufi(:,:,:,j)  !!(set lu corner)
            endif
          endif
        enddo
        call MPI_WAIT(isend_rq_TOll,isend_stat_Tll,mpierr)
        call MPI_WAIT(isend_rq_TOul,isend_stat_Tul,mpierr)
        call MPI_WAIT(isend_rq_TOuu,isend_stat_Tuu,mpierr)
        call MPI_WAIT(isend_rq_TOlu,isend_stat_Tlu,mpierr)
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
      call mpibarrier()
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
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbufx, j
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
            xuneigh,tolowx,MPI_COMM_WORLD,irecv_rq_fromuppx,mpierr)
        call MPI_IRECV(lbufxi(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xlneigh,touppx,MPI_COMM_WORLD,irecv_rq_fromlowx,mpierr)
        call MPI_ISEND(lbufxo(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xlneigh,tolowx,MPI_COMM_WORLD,isend_rq_tolowx,mpierr)
        call MPI_ISEND(ubufxo(:,:,:,ivar1:ivar2),nbufx,MPI_REAL, &
            xuneigh,touppx,MPI_COMM_WORLD,isend_rq_touppx,mpierr)
        call MPI_WAIT(irecv_rq_fromuppx,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowx,irecv_stat_fl,mpierr)
        do j=ivar1,ivar2
          if (.not. lfirst_proc_x .or. bcx1(j)=='p' .or. &
              (bcx1(j)=='she'.and.nygrid==1)) then
            f( 1:l1-1,m1:m2,n1:n2,j)=lbufxi(:,:,:,j)  !!(set lower buffer)
          endif
          if (.not. llast_proc_x .or. bcx2(j)=='p' .or. &
              (bcx2(j)=='she'.and.nygrid==1)) then
            f(l2+1:mx,m1:m2,n1:n2,j)=ubufxi(:,:,:,j)  !!(set upper buffer)
          endif
        enddo
        call MPI_WAIT(isend_rq_tolowx,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppx,isend_stat_tu,mpierr)
      endif
!
    endsubroutine isendrcv_bdry_x
!***********************************************************************
    subroutine isendrcv_bdry_npole(f,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values for pole.
!
!   18-june-10/dhruba: aped
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbuf_pole, j
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!
! The following is not a typo, it must be nprocz although the boundary
! is the pole (i.e., along the y direction).
      if (nprocz>1) then
        npbufyo(:,:,:,ivar1:ivar2)=f(l1:l2,m1:m1i,n1:n2,ivar1:ivar2) !!(north pole)
        nbuf_pole=nx*nghost*nz*(ivar2-ivar1+1)
        call MPI_IRECV(npbufyi(:,:,:,ivar1:ivar2),nbuf_pole,MPI_REAL, &
             poleneigh,npole_tag,MPI_COMM_WORLD,irecv_rq_npole,mpierr)
        call MPI_ISEND(npbufyo(:,:,:,ivar1:ivar2),nbuf_pole,MPI_REAL, &
             poleneigh,npole_tag,MPI_COMM_WORLD,isend_rq_npole,mpierr)
        call MPI_WAIT(irecv_rq_npole,irecv_stat_np,mpierr)
        do j=ivar1,ivar2
          if (bcy1(j)=='pp') then
             f(l1:l2,1,n1:n2,j)=npbufyi(:,3,:,j)
             f(l1:l2,2,n1:n2,j)=npbufyi(:,2,:,j)
             f(l1:l2,3,n1:n2,j)=npbufyi(:,1,:,j)
          endif
          if (bcy1(j)=='ap') then
             f(l1:l2,1,n1:n2,j)=-npbufyi(:,3,:,j)
             f(l1:l2,2,n1:n2,j)=-npbufyi(:,2,:,j)
             f(l1:l2,3,n1:n2,j)=-npbufyi(:,1,:,j)
          endif
        enddo
        call MPI_WAIT(isend_rq_npole,isend_stat_np,mpierr)
      endif
!
    endsubroutine isendrcv_bdry_npole
!***********************************************************************
    subroutine isendrcv_bdry_spole(f,ivar1_opt,ivar2_opt)
!
!  Isend and Irecv boundary values for pole.
!
!   18-june-10/dhruba: aped
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      integer :: ivar1, ivar2, nbuf_pole, j
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  The following is not a typo, it must be nprocz although the boundary
!  is the pole (i.e., along the y direction).
!
      if (nprocz>1) then
        spbufyo(:,:,:,ivar1:ivar2)=f(l1:l2,m2i:m2,n1:n2,ivar1:ivar2) !!(south pole)
        nbuf_pole=nx*nghost*nz*(ivar2-ivar1+1)
        call MPI_IRECV(spbufyi(:,:,:,ivar1:ivar2),nbuf_pole,MPI_REAL, &
             poleneigh,spole_tag,MPI_COMM_WORLD,irecv_rq_spole,mpierr)
        call MPI_ISEND(spbufyo(:,:,:,ivar1:ivar2),nbuf_pole,MPI_REAL, &
             poleneigh,spole_tag,MPI_COMM_WORLD,isend_rq_spole,mpierr)
        call MPI_WAIT(irecv_rq_spole,irecv_stat_spole,mpierr)
        do j=ivar1,ivar2
          if (bcy2(j)=='pp') &
              f(l1:l2,m2+1:my,n1:n2,j)=spbufyi(:,:,:,j)
          if (bcy2(j)=='ap') &
              f(l1:l2,m2+1:my,n1:n2,j)=-spbufyi(:,:,:,j)
        enddo
        call MPI_WAIT(isend_rq_spole,isend_stat_spole,mpierr)
      endif
!
    endsubroutine isendrcv_bdry_spole
!***********************************************************************
   subroutine initiate_shearing(f,ivar1_opt,ivar2_opt)
!
!  Subroutine for shearing sheet boundary conditions
!
!  20-june-02/nils: adapted from pencil_mpi
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      double precision :: deltay_dy, frac, c1, c2, c3, c4, c5, c6
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
          if (deltay>=0) then
            nextya=ipz*nprocy*nprocx+modulo(ipy-ystep,nprocy)  *nprocx+ipx_partner
            lastya=ipz*nprocy*nprocx+modulo(ipy-ystep-1,nprocy)*nprocx+ipx_partner
            lastyb=ipz*nprocy*nprocx+modulo(ipy+ystep,nprocy)  *nprocx+ipx_partner
            nextyb=ipz*nprocy*nprocx+modulo(ipy+ystep+1,nprocy)*nprocx+ipx_partner
          else
!
!  The following is probably not quite right (I see some imperfections
!  near the x boundaries.
!
            nextya=ipz*nprocy*nprocx+modulo(ipy-ystep,nprocy)  *nprocx+ipx_partner
            lastya=ipz*nprocy*nprocx+modulo(ipy-ystep-1,nprocy)*nprocx+ipx_partner
            lastyb=ipz*nprocy*nprocx+modulo(ipy+ystep,nprocy)  *nprocx+ipx_partner
            nextyb=ipz*nprocy*nprocx+modulo(ipy+ystep+1,nprocy)*nprocx+ipx_partner
          endif
!
          fao(:,:,:,ivar1:ivar2) = f(l1:l1i,:,:,ivar1:ivar2)
          fbo(:,:,:,ivar1:ivar2) = f(l2i:l2,:,:,ivar1:ivar2)
          nbufx_gh=my*mz*nghost*(ivar2-ivar1+1)
          if (lastya/=iproc) then
            call MPI_ISEND(fao(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,lastya, &
                tonextyb,MPI_COMM_WORLD,isend_rq_tolastya,mpierr)
          endif
          if (nextyb==iproc) then
            fbhi(:,:,:,ivar1:ivar2)=fao(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(fbhi(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,nextyb, &
                tonextyb,MPI_COMM_WORLD,irecv_rq_fromnextyb,mpierr)
          endif
          if (nextya/=iproc) then
            call MPI_ISEND(fao(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,nextya, &
                tolastyb,MPI_COMM_WORLD,isend_rq_tonextya,mpierr)
          endif
          if (lastyb==iproc) then
            fblo(:,:,:,ivar1:ivar2)=fao(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(fblo(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,lastyb, &
                tolastyb,MPI_COMM_WORLD,irecv_rq_fromlastyb,mpierr)
          endif
          if (lastyb/=iproc) then
            call MPI_ISEND(fbo(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,lastyb, &
                tonextya,MPI_COMM_WORLD,isend_rq_tolastyb,mpierr)
          endif
          if (nextya==iproc) then
            fahi(:,:,:,ivar1:ivar2)=fbo(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(fahi(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,nextya, &
                tonextya,MPI_COMM_WORLD,irecv_rq_fromnextya,mpierr)
          endif
          if (nextyb/=iproc) then
            call MPI_ISEND(fbo(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,nextyb, &
                tolastya,MPI_COMM_WORLD,isend_rq_tonextyb,mpierr)
          endif
          if (lastya==iproc) then
            falo(:,:,:,ivar1:ivar2)=fbo(:,:,:,ivar1:ivar2)
          else
            call MPI_IRECV(falo(:,:,:,ivar1:ivar2),nbufx_gh,MPI_REAL,lastya, &
                tolastya,MPI_COMM_WORLD,irecv_rq_fromlastya,mpierr)
          endif
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
!  02-mar-02/ulf: Sliding periodic boundary conditions in x
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer, optional :: ivar1_opt, ivar2_opt
!
      real, dimension (nghost,2*my-2*nghost,mz,mcom) :: fa, fb
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fal, irecv_stat_fan
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fbl, irecv_stat_fbn
      integer, dimension (MPI_STATUS_SIZE) :: isend_stat_tna, isend_stat_tla
      integer, dimension (MPI_STATUS_SIZE) :: isend_stat_tnb, isend_stat_tlb
      integer :: ivar1, ivar2, m2long
      double precision :: deltay_dy, frac, c1, c2, c3, c4, c5, c6
!
      ivar1=1; ivar2=mcom
      if (present(ivar1_opt)) ivar1=ivar1_opt
      if (present(ivar2_opt)) ivar2=ivar2_opt
!
!  Some special cases have already finished in initiate_shearing.
!
      if (nygrid/=1 .and. (nprocx>1 .or. nprocy>1) .and. &
          (lfirst_proc_x .or. llast_proc_x)) then
!
!  Need to wait till all communication has been recived.
!
        if (lastyb/=iproc) &
            call MPI_WAIT(irecv_rq_fromlastyb,irecv_stat_fbl,mpierr)
        if (nextyb/=iproc) &
            call MPI_WAIT(irecv_rq_fromnextyb,irecv_stat_fbn,mpierr)
        if (lastya/=iproc) &
            call MPI_WAIT(irecv_rq_fromlastya,irecv_stat_fal,mpierr)
        if (nextya/=iproc) &
            call MPI_WAIT(irecv_rq_fromnextya,irecv_stat_fan,mpierr)
!
!  Reading communicated information into f.
!
        deltay_dy=deltay/dy
        m2long = 2*my-3*nghost
        fa(:,1:m2,:,ivar1:ivar2) = falo(:,1:m2,:,ivar1:ivar2)
        fa(:,m2+1:2*my-2*nghost,:,ivar1:ivar2) = fahi(:,m1:my,:,ivar1:ivar2)
        fb(:,1:m2,:,ivar1:ivar2) = fblo(:,1:m2,:,ivar1:ivar2)
        fb(:,m2+1:2*my-2*nghost,:,ivar1:ivar2) = fbhi(:,m1:my,:,ivar1:ivar2)
        displs = modulo(int(deltay_dy),ny)
        frac = deltay_dy - int(deltay_dy)
        c1 = -          (frac+1.)*frac*(frac-1.)*(frac-2.)*(frac-3.)/120.
        c2 = +(frac+2.)          *frac*(frac-1.)*(frac-2.)*(frac-3.)/24.
        c3 = -(frac+2.)*(frac+1.)     *(frac-1.)*(frac-2.)*(frac-3.)/12.
        c4 = +(frac+2.)*(frac+1.)*frac          *(frac-2.)*(frac-3.)/12.
        c5 = -(frac+2.)*(frac+1.)*frac*(frac-1.)          *(frac-3.)/24.
        c6 = +(frac+2.)*(frac+1.)*frac*(frac-1.)*(frac-2.)          /120.
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
!  Need to wait till buffer is empty before re-using it again.
!
        if (nextyb/=iproc) call MPI_WAIT(isend_rq_tonextyb,isend_stat_tnb,mpierr)
        if (lastyb/=iproc) call MPI_WAIT(isend_rq_tolastyb,isend_stat_tlb,mpierr)
        if (nextya/=iproc) call MPI_WAIT(isend_rq_tonextya,isend_stat_tna,mpierr)
        if (lastya/=iproc) call MPI_WAIT(isend_rq_tolastya,isend_stat_tla,mpierr)
!
      endif
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
                    MPI_COMM_WORLD,irecv_zx,mpierr)
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
                    MPI_COMM_WORLD,irecv_xy,mpierr)
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
                    MPI_COMM_WORLD,isend_zx,mpierr)
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
                    MPI_COMM_WORLD,isend_xy,mpierr)
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
                        MPI_COMM_WORLD,isendrecv_zx,mpierr)
!
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
          MPI_COMM_YBEAM,mpierr)
!
      call MPI_ALLGATHER(Qrad_zx,nx*nz,MPI_REAL,Qrad_zx_all,nx*nz,MPI_REAL, &
          MPI_COMM_YBEAM,mpierr)
!
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
          tag_id, MPI_COMM_WORLD, stat, mpierr)
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
          tag_id, MPI_COMM_WORLD, stat, mpierr)
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
          tag_id, MPI_COMM_WORLD, stat, mpierr)
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
          tag_id, MPI_COMM_WORLD, stat, mpierr)
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
          tag_id, MPI_COMM_WORLD, stat, mpierr)
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
          tag_id, MPI_COMM_WORLD, stat, mpierr)
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
          tag_id, MPI_COMM_WORLD, stat, mpierr)
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
          tag_id, MPI_COMM_WORLD, stat, mpierr)
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
          tag_id, MPI_COMM_WORLD, stat, mpierr)
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
          tag_id, MPI_COMM_WORLD, stat, mpierr)
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
          tag_id, MPI_COMM_WORLD, mpierr)
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
          tag_id, MPI_COMM_WORLD, mpierr)
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
          tag_id, MPI_COMM_WORLD, mpierr)
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
          tag_id, MPI_COMM_WORLD,mpierr)
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
          tag_id, MPI_COMM_WORLD,mpierr)
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
          tag_id, MPI_COMM_WORLD,mpierr)
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
          tag_id, MPI_COMM_WORLD,mpierr)
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
          tag_id, MPI_COMM_WORLD, mpierr)
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
          tag_id, MPI_COMM_WORLD,mpierr)
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
          tag_id, MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpisend_int_arr2
!***********************************************************************
    subroutine mpibcast_logical_scl(lbcast_array,nbcast_array,proc)
!
!  Communicate logical scalar between processors.
!
      logical :: lbcast_array
      integer, optional :: nbcast_array,proc
      integer :: ibcast_proc,nbcast_arrayl
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      if (present(nbcast_array)) then
        nbcast_arrayl=nbcast_array
      else
        nbcast_arrayl=1
      endif
!
      call MPI_BCAST(lbcast_array,nbcast_arrayl,MPI_LOGICAL,ibcast_proc, &
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_logical_scl
!***********************************************************************
    subroutine mpibcast_logical_arr(lbcast_array,nbcast_array,proc)
!
!  Communicate logical array between processors.
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
          MPI_COMM_WORLD,mpierr)
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
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_logical_arr2
!***********************************************************************
    subroutine mpibcast_int_scl(ibcast_array,nbcast_array,proc)
!
!  Communicate integer scalar between processors.
!
      integer :: ibcast_array
      integer, optional :: nbcast_array,proc
      integer :: ibcast_proc,nbcast_arrayl
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      if (present(nbcast_array)) then
        nbcast_arrayl=nbcast_array
      else
        nbcast_arrayl=1
      endif
!
      call MPI_BCAST(ibcast_array,nbcast_arrayl,MPI_INTEGER,ibcast_proc, &
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_int_scl
!***********************************************************************
    subroutine mpibcast_int_arr(ibcast_array,nbcast_array,proc)
!
!  Communicate integer array between processors.
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
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_int_arr
!***********************************************************************
    subroutine mpibcast_real_scl(bcast_array,nbcast_array,proc)
!
!  Communicate real scalar between processors.
!
      real :: bcast_array
      integer, optional :: nbcast_array,proc
      integer :: ibcast_proc,nbcast_arrayl
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      if (present(nbcast_array)) then
        nbcast_arrayl=nbcast_array
      else
        nbcast_arrayl=1
      endif
!
      call MPI_BCAST(bcast_array,nbcast_arrayl,MPI_REAL,ibcast_proc, &
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_real_scl
!***********************************************************************
    subroutine mpibcast_real_arr(bcast_array,nbcast_array,proc)
!
!  Communicate real array between processors.
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
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_real_arr
!***********************************************************************
    subroutine mpibcast_real_arr2(bcast_array,nbcast_array,proc)
!
!  Communicate real array(:,:) to other processor.
!
!  25-feb-08/wlad: adapted
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
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_real_arr2
!***********************************************************************
    subroutine mpibcast_real_arr3(bcast_array,nb,proc)
!
!  Communicate real array(:,:,:) to other processor.
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
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_real_arr3
!***********************************************************************
    subroutine mpibcast_real_arr4(bcast_array,nb,proc)
!
!  Communicate real array(:,:,:,:) to other processor.
!
!  21-dec-10/ccyang: adapted
!
      integer, dimension(4) :: nb
      real, dimension(nb(1),nb(2),nb(3),nb(4)) :: bcast_array
      integer, optional :: proc
      integer :: ibcast_proc,nbcast
!
      nbcast=nb(1)*nb(2)*nb(3)*nb(4)
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      call MPI_BCAST(bcast_array, nbcast, MPI_REAL, ibcast_proc, &
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_real_arr4
!***********************************************************************
    subroutine mpibcast_double_scl(bcast_array,nbcast_array,proc)
!
!  Communicate real scalar between processors.
!
      double precision :: bcast_array
      integer, optional :: nbcast_array,proc
      integer :: ibcast_proc,nbcast_arrayl
!
      if (present(proc)) then
        ibcast_proc=proc
      else
        ibcast_proc=root
      endif
!
      if (present(nbcast_array)) then
        nbcast_arrayl=nbcast_array
      else
        nbcast_arrayl=1
      endif
!
      call MPI_BCAST(bcast_array,nbcast_arrayl,MPI_DOUBLE_PRECISION,ibcast_proc, &
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_double_scl
!***********************************************************************
    subroutine mpibcast_double_arr(bcast_array,nbcast_array,proc)
!
!  Communicate real array between processors.
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
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_double_arr
!***********************************************************************
    subroutine mpibcast_char_scl(cbcast_array,proc)
!
!  Communicate character scalar between processors.
!
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
      call MPI_BCAST(cbcast_array,1,MPI_CHARACTER,ibcast_proc, &
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_char_scl
!***********************************************************************
    subroutine mpibcast_char_arr(cbcast_array,nbcast_array,proc)
!
!  Communicate character array between processors.
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
          MPI_COMM_WORLD,mpierr)
!
    endsubroutine mpibcast_char_arr
!***********************************************************************
    subroutine mpiallreduce_sum_scl(fsum_tmp,fsum,idir)
!
!  Calculate total sum for each array element and return to all processors.
!
      real :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs
!
!  Sum over all processors and return to root (MPI_COMM_WORLD).
!  Sum over x beams and return to the ipx=0 processors (MPI_COMM_XBEAM).
!  Sum over y beams and return to the ipy=0 processors (MPI_COMM_YBEAM).
!  Sum over z beams and return to the ipz=0 processors (MPI_COMM_ZBEAM).
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=MPI_COMM_WORLD
      endif
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, 1, MPI_REAL, MPI_SUM, mpiprocs, mpierr)
!
    endsubroutine mpiallreduce_sum_scl
!***********************************************************************
    subroutine mpiallreduce_sum_arr(fsum_tmp,fsum,nreduce,idir)
!
!  Calculate total sum for each array element and return to all processors.
!
      integer :: nreduce
      real, dimension(nreduce) :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=MPI_COMM_WORLD
      endif
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, nreduce, MPI_REAL, MPI_SUM, &
          MPI_COMM_WORLD, mpierr)
!
    endsubroutine mpiallreduce_sum_arr
!***********************************************************************
    subroutine mpiallreduce_sum_arr2(fsum_tmp,fsum,nreduce,idir)
!
!  Calculate total sum for each array element and return to all processors.
!
!  23-nov-08/wlad: included the idir possibility
!
      integer, dimension(2) :: nreduce
      real, dimension(nreduce(1),nreduce(2)) :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=MPI_COMM_WORLD
      endif
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, product(nreduce), MPI_REAL, MPI_SUM, &
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
      integer :: mpiprocs
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=MPI_COMM_WORLD
      endif
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, product(nreduce), MPI_REAL, MPI_SUM, &
          mpiprocs, mpierr)
!
    endsubroutine mpiallreduce_sum_arr3
!***********************************************************************
    subroutine mpiallreduce_sum_int_scl(fsum_tmp,fsum)
!
!  Calculate total sum for each array element and return to all processors.
!
      integer :: fsum_tmp,fsum
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, 1, MPI_INTEGER, MPI_SUM, &
          MPI_COMM_WORLD, mpierr)
!
    endsubroutine mpiallreduce_sum_int_scl
!***********************************************************************
    subroutine mpiallreduce_sum_int_arr(fsum_tmp,fsum,nreduce)
!
!  Calculate total sum for each array element and return to all processors.
!
      integer :: nreduce
      integer, dimension(nreduce) :: fsum_tmp,fsum
!
      call MPI_ALLREDUCE(fsum_tmp, fsum, nreduce, MPI_INTEGER, MPI_SUM, &
          MPI_COMM_WORLD, mpierr)
!
    endsubroutine mpiallreduce_sum_int_arr
!***********************************************************************
    subroutine mpiallreduce_max_scl(fmax_tmp,fmax)
!
!  Calculate total maximum for each array element and return to all processors.
!
      real :: fmax_tmp,fmax
!
      call MPI_ALLREDUCE(fmax_tmp, fmax, 1, MPI_REAL, MPI_MAX, &
          MPI_COMM_WORLD, mpierr)
!
    endsubroutine mpiallreduce_max_scl
!***********************************************************************
    subroutine mpiallreduce_max_arr(fmax_tmp,fmax,nreduce)
!
!  Calculate total maximum for each array element and return to all processors.
!
      integer :: nreduce
      real, dimension(nreduce) :: fmax_tmp,fmax
!
      call MPI_ALLREDUCE(fmax_tmp, fmax, nreduce, MPI_REAL, MPI_MAX, &
          MPI_COMM_WORLD, mpierr)
!
    endsubroutine mpiallreduce_max_arr
!***********************************************************************
    subroutine mpireduce_max_scl(fmax_tmp,fmax)
!
!  Calculate total maximum for each array element and return to root.
!
      real :: fmax_tmp,fmax
!
      call MPI_REDUCE(fmax_tmp, fmax, 1, MPI_REAL, MPI_MAX, root, &
          MPI_COMM_WORLD, mpierr)
!
    endsubroutine mpireduce_max_scl
!***********************************************************************
    subroutine mpireduce_max_scl_int(fmax_tmp,fmax)
!
!  Calculate total maximum for each array element and return to root.
!
      integer :: fmax_tmp,fmax
!
      call MPI_REDUCE(fmax_tmp, fmax, 1, MPI_INTEGER, MPI_MAX, root, &
          MPI_COMM_WORLD, mpierr)
!
    endsubroutine mpireduce_max_scl_int
!***********************************************************************
    subroutine mpireduce_max_arr(fmax_tmp,fmax,nreduce)
!
!  Calculate total maximum for each array element and return to root.
!
      integer :: nreduce
      real, dimension(nreduce) :: fmax_tmp,fmax
!
      call MPI_REDUCE(fmax_tmp, fmax, nreduce, MPI_REAL, MPI_MAX, root, &
          MPI_COMM_WORLD, mpierr)
!
    endsubroutine mpireduce_max_arr
!***********************************************************************
    subroutine mpireduce_min_scl(fmin_tmp,fmin)
!
!  Calculate total minimum for each array element and return to root.
!
      real :: fmin_tmp,fmin
!
      call MPI_REDUCE(fmin_tmp, fmin, 1, MPI_REAL, MPI_MIN, root, &
          MPI_COMM_WORLD, mpierr)
!
    endsubroutine mpireduce_min_scl
!***********************************************************************
    subroutine mpireduce_min_arr(fmin_tmp,fmin,nreduce)
!
!  Calculate total maximum for each array element and return to root.
!
      integer :: nreduce
      real, dimension(nreduce) :: fmin_tmp,fmin
!
      call MPI_REDUCE(fmin_tmp, fmin, nreduce, MPI_REAL, MPI_MIN, root, &
          MPI_COMM_WORLD, mpierr)
!
    endsubroutine mpireduce_min_arr
!***********************************************************************
    subroutine mpireduce_sum_int_scl(fsum_tmp,fsum)
!
!  Calculate sum and return to root.
!
      integer :: fsum_tmp,fsum
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, 1, MPI_INTEGER, MPI_SUM, root, &
            MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_sum_int_scl
!***********************************************************************
    subroutine mpireduce_sum_int_arr(fsum_tmp,fsum,nreduce)
!
!  Calculate total sum for each array element and return to root.
!
      integer :: nreduce
      integer, dimension(nreduce) :: fsum_tmp,fsum
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, nreduce, MPI_INTEGER, MPI_SUM, root, &
            MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_sum_int_arr
!***********************************************************************
    subroutine mpireduce_sum_int_arr2(fsum_tmp,fsum,nreduce)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(2) :: nreduce
      integer, dimension(nreduce(1),nreduce(2)) :: fsum_tmp,fsum
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, nreduce, MPI_INTEGER, MPI_SUM, root, &
            MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_sum_int_arr2
!***********************************************************************
    subroutine mpireduce_sum_int_arr3(fsum_tmp,fsum,nreduce)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(3) :: nreduce
      integer, dimension(nreduce(1),nreduce(2),nreduce(3)) :: fsum_tmp,fsum
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, nreduce, MPI_INTEGER, MPI_SUM, root, &
            MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_sum_int_arr3
!***********************************************************************
    subroutine mpireduce_sum_int_arr4(fsum_tmp,fsum,nreduce)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(4) :: nreduce
      integer, dimension(nreduce(1),nreduce(2),nreduce(3),nreduce(4)) :: fsum_tmp,fsum
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        call MPI_REDUCE(fsum_tmp, fsum, nreduce, MPI_INTEGER, MPI_SUM, root, &
            MPI_COMM_WORLD, mpierr)
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
!  Sum over all processors and return to root (MPI_COMM_WORLD).
!  Sum over x beams and return to the ipx=0 processors (MPI_COMM_XBEAM).
!  Sum over y beams and return to the ipy=0 processors (MPI_COMM_YBEAM).
!  Sum over z beams and return to the ipz=0 processors (MPI_COMM_ZBEAM).
!
        if (present(idir)) then
          mpiprocs=mpigetcomm(idir)
        else
          mpiprocs=MPI_COMM_WORLD
        endif
        call MPI_REDUCE(fsum_tmp, fsum, 1, MPI_REAL, MPI_SUM, root, &
            mpiprocs, mpierr)
      endif
!
    endsubroutine mpireduce_sum_scl
!***********************************************************************
    subroutine mpireduce_sum_arr(fsum_tmp,fsum,nreduce,idir)
!
!  Calculate total sum for each array element and return to root.
!
      integer :: nreduce
      real, dimension(nreduce) :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs
!
      intent(in)  :: fsum_tmp,nreduce
      intent(out) :: fsum
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        if (present(idir)) then
          mpiprocs=mpigetcomm(idir)
        else
          mpiprocs=MPI_COMM_WORLD
        endif
        call MPI_REDUCE(fsum_tmp, fsum, nreduce, MPI_REAL, MPI_SUM, root, &
            mpiprocs, mpierr)
      endif
!
    endsubroutine mpireduce_sum_arr
!***********************************************************************
    subroutine mpireduce_sum_arr2(fsum_tmp,fsum,nreduce,idir)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(2) :: nreduce
      real, dimension(nreduce(1),nreduce(2)) :: fsum_tmp,fsum
      integer, optional :: idir
!
      integer :: mpiprocs
!
      intent(in)  :: fsum_tmp,nreduce
      intent(out) :: fsum
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        if (present(idir)) then
          mpiprocs=mpigetcomm(idir)
        else
          mpiprocs=MPI_COMM_WORLD
        endif
        call MPI_REDUCE(fsum_tmp, fsum, product(nreduce), MPI_REAL, MPI_SUM, &
            root, mpiprocs, mpierr)
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
      integer :: mpiprocs
!
      intent(in)  :: fsum_tmp,nreduce
      intent(out) :: fsum
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        if (present(idir)) then
          mpiprocs=mpigetcomm(idir)
        else
          mpiprocs=MPI_COMM_WORLD
        endif
        call MPI_REDUCE(fsum_tmp, fsum, product(nreduce), MPI_REAL, MPI_SUM, &
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
      integer :: mpiprocs
!
      intent(in)  :: fsum_tmp,nreduce
      intent(out) :: fsum
!
      if (nprocs==1) then
        fsum=fsum_tmp
      else
        if (present(idir)) then
          mpiprocs=mpigetcomm(idir)
        else
          mpiprocs=MPI_COMM_WORLD
        endif
        call MPI_REDUCE(fsum_tmp, fsum, product(nreduce), MPI_REAL, MPI_SUM, &
            root, mpiprocs, mpierr)
      endif
!
    endsubroutine mpireduce_sum_arr4
!***********************************************************************
    subroutine mpireduce_sum_double_scl(dsum_tmp,dsum)
!
!  Calculate total sum and return to root.
!
      double precision :: dsum_tmp,dsum
!
      if (nprocs==1) then
        dsum=dsum_tmp
      else
        call MPI_REDUCE(dsum_tmp, dsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            root, MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_sum_double_scl
!***********************************************************************
    subroutine mpireduce_sum_double_arr(dsum_tmp,dsum,nreduce)
!
!  Calculate total sum for each array element and return to root.
!
      integer :: nreduce
      double precision, dimension(nreduce) :: dsum_tmp,dsum
!
      if (nprocs==1) then
        dsum=dsum_tmp
      else
        call MPI_REDUCE(dsum_tmp, dsum, nreduce, MPI_DOUBLE_PRECISION, &
            MPI_SUM, root, MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_sum_double_arr
!***********************************************************************
    subroutine mpireduce_sum_double_arr2(dsum_tmp,dsum,nreduce)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(2) :: nreduce
      double precision, dimension(nreduce(1),nreduce(2)) :: dsum_tmp,dsum
!
      if (nprocs==1) then
        dsum=dsum_tmp
      else
        call MPI_REDUCE(dsum_tmp, dsum, product(nreduce), MPI_DOUBLE_PRECISION,&
            MPI_SUM, root, MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_sum_double_arr2
!***********************************************************************
    subroutine mpireduce_sum_double_arr3(dsum_tmp,dsum,nreduce)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(3) :: nreduce
      double precision, dimension(nreduce(1),nreduce(2),nreduce(3)) :: dsum_tmp,dsum
!
      if (nprocs==1) then
        dsum=dsum_tmp
      else
        call MPI_REDUCE(dsum_tmp, dsum, product(nreduce), MPI_DOUBLE_PRECISION,&
            MPI_SUM, root, MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_sum_double_arr3
!***********************************************************************
    subroutine mpireduce_sum_double_arr4(dsum_tmp,dsum,nreduce)
!
!  Calculate total sum for each array element and return to root.
!
      integer, dimension(4) :: nreduce
      double precision, dimension(nreduce(1),nreduce(2),nreduce(3),nreduce(4)) :: dsum_tmp,dsum
!
      if (nprocs==1) then
        dsum=dsum_tmp
      else
        call MPI_REDUCE(dsum_tmp, dsum, product(nreduce), MPI_DOUBLE_PRECISION,&
            MPI_SUM, root, MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_sum_double_arr4
!***********************************************************************
    subroutine mpireduce_or_scl(flor_tmp,flor)
!
!  Calculate logical or over all procs and return to root.
!
!  17-sep-05/anders: coded
!
      logical :: flor_tmp, flor
!
      if (nprocs==1) then
        flor=flor_tmp
      else
        call MPI_REDUCE(flor_tmp, flor, 1, MPI_LOGICAL, MPI_LOR, root, &
                        MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_or_scl
!***********************************************************************
    subroutine mpireduce_or_arr(flor_tmp,flor,nreduce)
!
!  Calculate logical or over all procs and return to root.
!
!  17-sep-05/anders: coded
!
      integer :: nreduce
      logical, dimension(nreduce) :: flor_tmp, flor
!
      if (nprocs==1) then
        flor=flor_tmp
      else
        call MPI_REDUCE(flor_tmp, flor, nreduce, MPI_LOGICAL, MPI_LOR, root, &
                        MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_or_arr
!***********************************************************************
    subroutine mpireduce_and_scl(fland_tmp,fland)
!
!  Calculate logical and over all procs and return to root.
!
!  17-sep-05/anders: coded
!
      logical :: fland_tmp, fland
!
      if (nprocs==1) then
        fland=fland_tmp
      else
        call MPI_REDUCE(fland_tmp, fland, 1, MPI_LOGICAL, MPI_LAND, root, &
                        MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_and_scl
!***********************************************************************
    subroutine mpireduce_and_arr(fland_tmp,fland,nreduce)
!
!  Calculate logical and over all procs and return to root.
!
!  11-mar-09/anders: coded
!
      integer :: nreduce
      logical, dimension(nreduce) :: fland_tmp, fland
!
      if (nprocs==1) then
        fland=fland_tmp
      else
        call MPI_REDUCE(fland_tmp, fland, nreduce, MPI_LOGICAL, MPI_LAND, root,&
                        MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine mpireduce_and_arr
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
      serial_level=serial_level+1
      if (serial_level>1) return
!
      buf = 0
      if (.not. lroot) then     ! root starts, others wait for permission
        call MPI_RECV(buf,1,MPI_INTEGER,root,io_perm,MPI_COMM_WORLD,status,mpierr)
      endif
!
    endsubroutine start_serialize
!***********************************************************************
    subroutine end_serialize()
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
          call MPI_SEND(buf,1,MPI_INTEGER,i,io_perm,MPI_COMM_WORLD,mpierr)
          call MPI_RECV(buf,1,MPI_INTEGER,i,io_succ,MPI_COMM_WORLD,status,mpierr)
        enddo
      else                  ! tell root we're done
        call MPI_SEND(buf,1,MPI_INTEGER,root,io_succ,MPI_COMM_WORLD,mpierr)
      endif
!
    endsubroutine end_serialize
!***********************************************************************
    subroutine mpibarrier()
!
!  Synchronize nodes.
!
!  23-jul-2002/wolf: coded
!
      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
!
    endsubroutine mpibarrier
!***********************************************************************
    subroutine mpifinalize()
!
      call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
      call MPI_FINALIZE(mpierr)
!
    endsubroutine mpifinalize
!***********************************************************************
    function mpiwtime()
!
      double precision :: mpiwtime
      double precision :: MPI_WTIME   ! definition needed for mpicomm_ to work
!
      mpiwtime = MPI_WTIME()
!
    endfunction mpiwtime
!***********************************************************************
    function mpiwtick()
!
      double precision :: mpiwtick
      double precision :: MPI_WTICK   ! definition needed for mpicomm_ to work
!
      mpiwtick = MPI_WTICK()
!
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
!  processors agree to call die_gracefully().
!
!  29-jun-05/tony: coded
!
!  Tell the world something went wrong -- mpirun may not propagate
!  an error status.
!
      if (lroot) call touch_file('ERROR')
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
    subroutine die_immediately()
!
!  Stop without shuting down MPI
!  For those MPI implementations, which only finalize when all
!  processors agree to finalize.
!
!  29-jun-05/tony: coded
!
!  Tell the world something went wrong -- mpirun may not propagate
!  an error status.
!
      if (lroot) call touch_file('ERROR')
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
        call die_immediately()
      else
        call die_gracefully()
      endif
!
    endsubroutine stop_fatal
!***********************************************************************
    subroutine stop_it(msg,code)
!
!  Print message and stop.
!  With at least some MPI implementations, this only stops if all
!  processors agree to call stop_it(). To stop (collectively) if only one
!  or a few processors find some condition, use stop_it_if_any().
!
!  6-nov-01/wolf: coded
!  4-nov-11/MR: optional parameter code added
!
      use general, only: itoa
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
      call die_gracefully()
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
!  on whether to call stop_it():
!
      call mpibarrier()
      call MPI_ALLREDUCE(stop_flag,global_stop_flag,1,MPI_LOGICAL, &
                         MPI_LOR,MPI_COMM_WORLD,mpierr)
      call MPI_ALLREDUCE(stop_flag,identical_stop_flag,1,MPI_LOGICAL, &
                         MPI_LAND,MPI_COMM_WORLD,mpierr)
!
      if (global_stop_flag) then
        if ((.not. lroot) .and. (.not. identical_stop_flag) .and. (msg/='')) &
            write(*,'(A,I8,A,A)') 'RANK ', iproc, ' STOPPED: ', msg
        if (present (code)) then
          call stop_it(msg, code)
        else
          call stop_it(msg)
        endif
      endif
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
                         MPI_LOR,MPI_COMM_WORLD,mpierr)
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
      integer :: ystag=111,yrtag=112,zstag=113,zrtag=114,partner
      integer :: m,n,ibox,ix
!
!  Doing x-y transpose if var='y'
!
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
                call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ystag,MPI_COMM_WORLD,mpierr)
                call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,yrtag,MPI_COMM_WORLD,stat,mpierr)
              elseif (px>ipy) then  ! below diagonal: receive first, send then
                call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ystag,MPI_COMM_WORLD,stat,mpierr)
                call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,yrtag,MPI_COMM_WORLD,mpierr)
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
              call MPI_SEND (send_buf_z,sendc_z,MPI_REAL,partner,zstag,MPI_COMM_WORLD,mpierr)
              call MPI_RECV (recv_buf_z,recvc_z,MPI_REAL,partner,zrtag,MPI_COMM_WORLD,stat,mpierr)
            elseif (px>ipz) then  ! below diagonal: receive first, send then
              call MPI_RECV (recv_buf_z,recvc_z,MPI_REAL,partner,zstag,MPI_COMM_WORLD,stat,mpierr)
              call MPI_SEND(send_buf_z,sendc_z,MPI_REAL,partner,zrtag,MPI_COMM_WORLD,mpierr)
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
!  Synchronize; not strictly necessary, so Axel will probably remove it..
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
! TODO: Implement nygrid = n*nxgrid
!
      real, dimension(nx,ny), intent(inout) :: a
!
      real, dimension(ny,ny) :: send_buf_y, recv_buf_y, tmp
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc_y,recvc_y,px
      integer :: ytag=101,partner
      integer :: ibox,iy
!
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
              call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,mpierr)
              call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,mpierr)
            elseif (px>ipy) then  ! below diagonal: receive first, send then
              call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,mpierr)
              call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,mpierr)
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
!
    endsubroutine transp_xy
!***********************************************************************
    subroutine transp_xy_other(a)
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
!
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
              call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,mpierr)
              call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,mpierr)
            elseif (px>ipy) then  ! below diagonal: receive first, send then
              call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,mpierr)
              call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,mpierr)
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
!
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
      integer :: ytag=101,ztag=202,partner
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
                call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,mpierr)
                call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,mpierr)
              elseif (px>ipy) then  ! below diagonal: receive first, send then
                call MPI_RECV(recv_buf_y,recvc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,stat,mpierr)
                call MPI_SEND(send_buf_y,sendc_y,MPI_REAL,partner,ytag,MPI_COMM_WORLD,mpierr)
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
              call MPI_SEND(send_buf_z,sendc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,mpierr)
              call MPI_RECV (recv_buf_z,recvc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,stat,mpierr)
            elseif (px>ipz) then  ! below diagonal: receive first, send then
              call MPI_RECV (recv_buf_z,recvc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,stat,mpierr)
              call MPI_SEND(send_buf_z,sendc_z,MPI_REAL,partner,ztag,MPI_COMM_WORLD,mpierr)
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
      real, dimension(nxt,nz) :: buf
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc,px
      integer :: ztag=101,partner
!
      if (mod(nxgrid,nprocz)/=0) then
        print*,'transp_xz: nxgrid needs to be an integer multiple of nprocz'
        call stop_it_if_any(.true.,'Inconsistency: mod(nxgrid,nprocz)/=0')
      endif
!
!  Calculate the size of buffers.
!  Buffers used for the y-transpose have the same size in y and z.
!
      sendc=nx*nz
!
!  Send information to different processors (x-z transpose)
!
      b(ipz*nz+1:(ipz+1)*nz,:)=transpose(a(ipz*nxt+1:(ipz+1)*nxt,:))
      do px=0,nprocz-1
        if (px/=ipz) then
          partner=ipy+px*nprocy ! = iproc + (px-ipz)*nprocy
          buf=a(px*nxt+1:(px+1)*nxt,:)
          call MPI_SENDRECV_REPLACE(buf,sendc,MPI_REAL,partner,ztag,partner,ztag,MPI_COMM_WORLD,stat,mpierr)
          b(px*nz+1:(px+1)*nz,:)=transpose(buf)
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
      real, dimension(nz,nxt) :: buf
      integer, dimension(MPI_STATUS_SIZE) :: stat
      integer :: sendc,px
      integer :: ztag=101,partner
!
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
          call MPI_SENDRECV_REPLACE(buf,sendc,MPI_REAL,partner,ztag,partner,ztag,MPI_COMM_WORLD,stat,mpierr)
          b(px*nxt+1:(px+1)*nxt,:)=transpose(buf)
        endif
      enddo
!
    endsubroutine transp_zx
!***********************************************************************
    subroutine communicate_vect_field_ghosts(f,topbot,start_index)
!
!  Helper routine for communication of ghost cell values of a vector field.
!  Needed by potential field extrapolations, which only compute nx*ny arrays.
!  Can also be used for synchronization of changed uu values with ghost cells,
!  if the start_index parameter set to iux (default is iax).
!
!   8-oct-2006/tobi: Coded
!  28-dec-2010/Bourdin.KIS: extended to work for any 3D vector field data.
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot
      integer, intent(in), optional :: start_index
!
      real, dimension (nx,nghost,nghost+1,3) :: lbufyo,ubufyo,lbufyi,ubufyi
      real, dimension (nghost,my,nghost+1,3) :: lbufxo,ubufxo,lbufxi,ubufxi
      integer :: nbufx,nbufy,nn1,nn2,is,ie
!
      is = iax
      if (present (start_index)) is = start_index
      ie = is + 2
!
      nn1=-1; nn2=-1
      select case (topbot)
        case ('bot'); nn1=1;  nn2=n1
        case ('top'); nn1=n2; nn2=mz
        case default; call stop_it_if_any(.true.,"communicate_vect_field_ghosts: "//topbot//" should be either `top' or `bot'")
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
        call MPI_IRECV(ubufyi,nbufy,MPI_REAL,yuneigh,tolowy, &
                       MPI_COMM_WORLD,irecv_rq_fromuppy,mpierr)
        call MPI_IRECV(lbufyi,nbufy,MPI_REAL,ylneigh,touppy, &
                       MPI_COMM_WORLD,irecv_rq_fromlowy,mpierr)
        call MPI_ISEND(lbufyo,nbufy,MPI_REAL,ylneigh,tolowy, &
                       MPI_COMM_WORLD,isend_rq_tolowy,mpierr)
        call MPI_ISEND(ubufyo,nbufy,MPI_REAL,yuneigh,touppy, &
                       MPI_COMM_WORLD,isend_rq_touppy,mpierr)
!
        call MPI_WAIT(irecv_rq_fromuppy,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowy,irecv_stat_fl,mpierr)
!
        f(l1:l2,   1:m1-1,nn1:nn2,is:ie) = lbufyi
        f(l1:l2,m2+1:my  ,nn1:nn2,is:ie) = ubufyi
!
        call MPI_WAIT(isend_rq_tolowy,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppy,isend_stat_tu,mpierr)
!
      else
!
        f(l1:l2,   1:m1-1,nn1:nn2,is:ie) = f(l1:l2,m2i:m2 ,nn1:nn2,is:ie)
        f(l1:l2,m2+1:my  ,nn1:nn2,is:ie) = f(l1:l2, m1:m1i,nn1:nn2,is:ie)
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
        nbufx=nghost*my*(nghost+1)*3
!
        call MPI_IRECV(ubufxi,nbufx,MPI_REAL,xuneigh,tolowx, &
                       MPI_COMM_WORLD,irecv_rq_fromuppx,mpierr)
        call MPI_IRECV(lbufxi,nbufx,MPI_REAL,xlneigh,touppx, &
                       MPI_COMM_WORLD,irecv_rq_fromlowx,mpierr)
        call MPI_ISEND(lbufxo,nbufx,MPI_REAL,xlneigh,tolowx, &
                       MPI_COMM_WORLD,isend_rq_tolowx,mpierr)
        call MPI_ISEND(ubufxo,nbufx,MPI_REAL,xuneigh,touppx, &
                       MPI_COMM_WORLD,isend_rq_touppx,mpierr)
!
        call MPI_WAIT(irecv_rq_fromuppx,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowx,irecv_stat_fl,mpierr)
!
        f(   1:l1-1,:,nn1:nn2,is:ie) = lbufxi
        f(l2+1:mx  ,:,nn1:nn2,is:ie) = ubufxi
!
        call MPI_WAIT(isend_rq_tolowx,isend_stat_tl,mpierr)
        call MPI_WAIT(isend_rq_touppx,isend_stat_tu,mpierr)
!
      else
!
        f(   1:l1-1,:,nn1:nn2,is:ie) = f(l2i:l2 ,:,nn1:nn2,is:ie)
        f(l2+1:mx  ,:,nn1:nn2,is:ie) = f( l1:l1i,:,nn1:nn2,is:ie)
!
      endif
!
    endsubroutine communicate_vect_field_ghosts
!***********************************************************************
    subroutine communicate_xy_ghosts(data)
!
!  Helper routine for communication of ghost cells in horizontal direction.
!
!  11-apr-2011/Bourdin.KIS: adapted from communicate_vect_field_ghosts.
!
      real, dimension (mx,my), intent (inout) :: data
!
      real, dimension (nx,nghost) :: lbufyo,ubufyo,lbufyi,ubufyi
      real, dimension (nghost,my) :: lbufxo,ubufxo,lbufxi,ubufxi
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
        call MPI_IRECV (ubufyi, nbufy, MPI_REAL, yuneigh, tolowy, &
                       MPI_COMM_WORLD, irecv_rq_fromuppy, mpierr)
        call MPI_IRECV (lbufyi, nbufy, MPI_REAL, ylneigh, touppy, &
                       MPI_COMM_WORLD, irecv_rq_fromlowy, mpierr)
        call MPI_ISEND (lbufyo, nbufy, MPI_REAL, ylneigh, tolowy, &
                       MPI_COMM_WORLD, isend_rq_tolowy, mpierr)
        call MPI_ISEND (ubufyo, nbufy, MPI_REAL, yuneigh, touppy, &
                       MPI_COMM_WORLD, isend_rq_touppy, mpierr)
!
        call MPI_WAIT (irecv_rq_fromuppy, irecv_stat_fu, mpierr)
        call MPI_WAIT (irecv_rq_fromlowy, irecv_stat_fl, mpierr)
!
        data(l1:l2,   1:m1-1) = lbufyi
        data(l1:l2,m2+1:my  ) = ubufyi
!
        call MPI_WAIT (isend_rq_tolowy, isend_stat_tl, mpierr)
        call MPI_WAIT (isend_rq_touppy, isend_stat_tu, mpierr)
!
      else
!
        data(l1:l2,   1:m1-1) = data(l1:l2,m2i:m2 )
        data(l1:l2,m2+1:my  ) = data(l1:l2, m1:m1i)
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
        nbufx = nghost * my
!
        call MPI_IRECV (ubufxi, nbufx, MPI_REAL, xuneigh, tolowx, &
                       MPI_COMM_WORLD, irecv_rq_fromuppx, mpierr)
        call MPI_IRECV (lbufxi, nbufx, MPI_REAL, xlneigh, touppx, &
                       MPI_COMM_WORLD, irecv_rq_fromlowx, mpierr)
        call MPI_ISEND (lbufxo, nbufx, MPI_REAL, xlneigh, tolowx, &
                       MPI_COMM_WORLD, isend_rq_tolowx, mpierr)
        call MPI_ISEND (ubufxo, nbufx, MPI_REAL, xuneigh, touppx, &
                       MPI_COMM_WORLD, isend_rq_touppx, mpierr)
!
        call MPI_WAIT (irecv_rq_fromuppx, irecv_stat_fu, mpierr)
        call MPI_WAIT (irecv_rq_fromlowx, irecv_stat_fl, mpierr)
!
        data(   1:l1-1,:) = lbufxi
        data(l2+1:mx  ,:) = ubufxi
!
        call MPI_WAIT (isend_rq_tolowx, isend_stat_tl, mpierr)
        call MPI_WAIT (isend_rq_touppx, isend_stat_tu, mpierr)
!
      else
!
        data(   1:l1-1,:) = data(l2i:l2 ,:)
        data(l2+1:mx  ,:) = data( l1:l1i,:)
!
      endif
!
    endsubroutine communicate_xy_ghosts
!***********************************************************************
    subroutine fill_zghostzones_3vec(vec,ivar)
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
!
      if (nprocz>1) then
!
        lbufo = vec(n1:n1i,:)                        !!(lower z-zone)
        ubufo = vec(n2i:n2,:)                        !!(upper z-zone)
!
        nbuf=nghost*3
!
        call MPI_IRECV(ubufi,nbuf,MPI_REAL, &
                       zuneigh,tolowz,MPI_COMM_WORLD,irecv_rq_fromuppz,mpierr)
        call MPI_IRECV(lbufi,nbuf,MPI_REAL, &
                       zlneigh,touppz,MPI_COMM_WORLD,irecv_rq_fromlowz,mpierr)
!
        call MPI_ISEND(lbufo,nbuf,MPI_REAL, &
                       zlneigh,tolowz,MPI_COMM_WORLD,isend_rq_tolowz,mpierr)
        call MPI_ISEND(ubufo,nbuf,MPI_REAL, &
                       zuneigh,touppz,MPI_COMM_WORLD,isend_rq_touppz,mpierr)
!
        call MPI_WAIT(irecv_rq_fromuppz,irecv_stat_fu,mpierr)
        call MPI_WAIT(irecv_rq_fromlowz,irecv_stat_fl,mpierr)
!
        do j=1,3
!
          if (.not. lfirst_proc_z .or. bcz1(j-1+ivar)=='p') &
            vec(1:n1-1,j)=lbufi(:,j)
!
!  Read from buffer in lower ghostzones.
!
          if (.not. llast_proc_z .or. bcz2(j-1+ivar)=='p') &
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
          if ( bcz1(ivar+j-1)=='p' ) then
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
!  This routine needs only to be called from all processors a the xy-plane.
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
!
      if (lfirst_proc_xy) then
        ! initialize sum with the local data
        sum = in
        ! collect and sum up the remote data
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) cycle
            call mpirecv_real (buffer, 1, partner, tag)
            sum = sum + buffer
          enddo
        enddo
      else
        ! send data to collector
        call mpisend_real (in, 1, ipz*nprocxy, tag)
        sum = 0.0
      endif
!
      ! distribute back the sum
      call distribute_xy (sum, out)
!
    endsubroutine sum_xy
!***********************************************************************
    subroutine distribute_xy_0D(in, out, source_proc)
!
!  This routine distributes a scalar on the source processor
!  to all processors in the xy-plane.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  25-jan-2012/Bourdin.KIS: coded
!
      real, intent(in) :: in
      real, intent(out) :: out
      integer, intent(in), optional :: source_proc
!
      integer :: px, py, broadcaster, partner
      integer, parameter :: ytag=115
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
!
      broadcaster = ipz * nprocxy
      if (present (source_proc)) broadcaster = broadcaster + source_proc
!
      if (iproc == broadcaster) then
        ! distribute the data
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) then
              ! data is local
              out = in
            else
              ! send to partner
              call MPI_SEND (in, 1, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            endif
          enddo
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, 1, MPI_REAL, broadcaster, ytag, MPI_COMM_WORLD, stat, mpierr)
      endif
!
    endsubroutine distribute_xy_0D
!***********************************************************************
    subroutine distribute_xy_2D(in, out, source_proc)
!
!  This routine divides a large array of 2D data on the source processor
!  and distributes it to all processors in the xy-plane.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  08-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
      integer, intent(in), optional :: source_proc
!
      integer :: bnx, bny ! transfer box sizes
      integer :: px, py, broadcaster, partner, nbox
      integer, parameter :: ytag=115
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
!
      bnx = size (out, 1)
      bny = size (out, 2)
      nbox = bnx*bny
!
      broadcaster = ipz * nprocxy
      if (present (source_proc)) broadcaster = broadcaster + source_proc
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
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) then
              ! data is local
              out = in(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny)
            else
              ! send to partner
              call MPI_SEND (in(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny), &
                  nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            endif
          enddo
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, MPI_REAL, broadcaster, ytag, MPI_COMM_WORLD, stat, mpierr)
      endif
!
    endsubroutine distribute_xy_2D
!***********************************************************************
    subroutine distribute_xy_3D(in, out, source_proc)
!
!  This routine divides a large array of 3D data on the source processor
!  and distributes it to all processors in the xy-plane.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  08-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
      integer, intent(in), optional :: source_proc
!
      integer :: bnx, bny, bnz ! transfer box sizes
      integer :: px, py, broadcaster, partner, nbox
      integer, parameter :: ytag=115
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
!
      bnx = size (out, 1)
      bny = size (out, 2)
      bnz = size (out, 3)
      nbox = bnx*bny*bnz
!
      broadcaster = ipz * nprocxy
      if (present (source_proc)) broadcaster = broadcaster + source_proc
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
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) then
              ! data is local
              out = in(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:)
            else
              ! send to partner
              call MPI_SEND (in(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:), &
                  nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            endif
          enddo
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, MPI_REAL, broadcaster, ytag, MPI_COMM_WORLD, stat, mpierr)
      endif
!
    endsubroutine distribute_xy_3D
!***********************************************************************
    subroutine distribute_xy_4D(in, out, source_proc)
!
!  This routine divides a large array of 4D data on the source processor
!  and distributes it to all processors in the xy-plane.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding xy-plane (Default: 0, equals lfirst_proc_xy).
!
!  08-jan-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
      integer, intent(in), optional :: source_proc
!
      integer :: bnx, bny, bnz, bna ! transfer box sizes
      integer :: px, py, broadcaster, partner, nbox
      integer, parameter :: ytag=115
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
!
      bnx = size (out, 1)
      bny = size (out, 2)
      bnz = size (out, 3)
      bna = size (out, 4)
      nbox = bnx*bny*bnz*bna
!
      broadcaster = ipz * nprocxy
      if (present (source_proc)) broadcaster = broadcaster + source_proc
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
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) then
              ! data is local
              out = in(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:,:)
            else
              ! send to partner
              call MPI_SEND (in(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:,:), &
                  nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            endif
          enddo
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, MPI_REAL, broadcaster, ytag, MPI_COMM_WORLD, stat, mpierr)
      endif
!
    endsubroutine distribute_xy_4D
!***********************************************************************
    subroutine collect_xy_0D(in, out, dest_proc)
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
!
      integer :: px, py, collector, partner
      integer, parameter :: ytag=116
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real :: buffer
!
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
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) then
              ! data is local
              out(px,py) = in
            else
              ! receive from partner
              call MPI_RECV (buffer, 1, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
              out(px,py) = buffer
            endif
          enddo
        enddo
      else
        ! send to collector
        call MPI_SEND (in, 1, MPI_REAL, collector, ytag, MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine collect_xy_0D
!***********************************************************************
    subroutine collect_xy_2D(in, out, dest_proc)
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
!
      integer :: bnx, bny ! transfer box sizes
      integer :: px, py, collector, partner, nbox, alloc_err
      integer, parameter :: ytag=116
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: buffer
!
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
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) then
              ! data is local
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny) = in
            else
              ! receive from partner
              call MPI_RECV (buffer, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny) = buffer
            endif
          enddo
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, nbox, MPI_REAL, collector, ytag, MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine collect_xy_2D
!***********************************************************************
    subroutine collect_xy_3D(in, out, dest_proc)
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
!
      integer :: bnx, bny, bnz ! transfer box sizes
      integer :: px, py, collector, partner, nbox, alloc_err
      integer, parameter :: ytag=116
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:), allocatable :: buffer
!
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
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) then
              ! data is local
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:) = in
            else
              ! receive from partner
              call MPI_RECV (buffer, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:) = buffer
            endif
          enddo
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, nbox, MPI_REAL, collector, ytag, MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine collect_xy_3D
!***********************************************************************
    subroutine collect_xy_4D(in, out, dest_proc)
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
!
      integer :: bnx, bny, bnz, bna ! transfer box sizes
      integer :: px, py, collector, partner, nbox, alloc_err
      integer, parameter :: ytag=116
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: buffer
!
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
            partner = px + py*nprocx + ipz*nprocxy
            if (iproc == partner) then
              ! data is local
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:,:) = in
            else
              ! receive from partner
              call MPI_RECV (buffer, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
              out(px*bnx+1:(px+1)*bnx,py*bny+1:(py+1)*bny,:,:) = buffer
            endif
          enddo
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, nbox, MPI_REAL, collector, ytag, MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine collect_xy_4D
!***********************************************************************
    subroutine distribute_z_3D(in, out, source_proc)
!
!  This routine divides a large array of 3D data on the source processor
!  and distributes it to all processors in the z-direction.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding z-direction (Default: 0, equals lfirst_proc_z).
!
!  08-mar-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
      integer, intent(in), optional :: source_proc
!
      integer :: bnx, bny, bnz ! transfer box sizes
      integer :: pz, broadcaster, partner, nbox
      integer, parameter :: ytag=117
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
!
      bnx = size (out, 1)
      bny = size (out, 2)
      bnz = size (out, 3)
      nbox = bnx*bny*bnz
!
      broadcaster = ipx + ipy*nprocx
      if (present (source_proc)) broadcaster = broadcaster + source_proc*nprocxy
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
          partner = ipx + ipy*nprocx + pz*nprocxy
          if (iproc == partner) then
            ! data is local
            out = in(:,:,pz*bnz+1:(pz+1)*bnz)
          else
            ! send to partner
            call MPI_SEND (in(:,:,pz*bnz+1:(pz+1)*bnz), nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, MPI_REAL, broadcaster, ytag, MPI_COMM_WORLD, stat, mpierr)
      endif
!
    endsubroutine distribute_z_3D
!***********************************************************************
    subroutine distribute_z_4D(in, out, source_proc)
!
!  This routine divides a large array of 4D data on the source processor
!  and distributes it to all processors in the z-direction.
!  'source_proc' is the iproc number relative to the first processor
!  in the corresponding z-direction (Default: 0, equals lfirst_proc_z).
!
!  08-mar-2011/Bourdin.KIS: coded
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
      integer, intent(in), optional :: source_proc
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
      broadcaster = ipx + ipy*nprocx
      if (present (source_proc)) broadcaster = broadcaster + source_proc*nprocxy
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
          partner = ipx + ipy*nprocx + pz*nprocxy
          if (iproc == partner) then
            ! data is local
            out = in(:,:,pz*bnz+1:(pz+1)*bnz,:)
          else
            ! send to partner
            call MPI_SEND (in(:,:,pz*bnz+1:(pz+1)*bnz,:), nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, nbox, MPI_REAL, broadcaster, ytag, MPI_COMM_WORLD, stat, mpierr)
      endif
!
    endsubroutine distribute_z_4D
!***********************************************************************
    subroutine collect_z_3D(in, out, dest_proc)
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
!
      integer :: bnx, bny, bnz ! transfer box sizes
      integer :: pz, collector, partner, nbox, alloc_err
      integer, parameter :: ytag=118
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:), allocatable :: buffer
!
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
          partner = ipx + ipy*nprocx + pz*nprocxy
          if (iproc == partner) then
            ! data is local
            out(:,:,pz*bnz+1:(pz+1)*bnz) = in
          else
            ! receive from partner
            call MPI_RECV (buffer, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            out(:,:,pz*bnz+1:(pz+1)*bnz) = buffer
          endif
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, nbox, MPI_REAL, collector, ytag, MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine collect_z_3D
!***********************************************************************
    subroutine collect_z_4D(in, out, dest_proc)
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
!
      integer :: bnx, bny, bnz, bna ! transfer box sizes
      integer :: pz, collector, partner, nbox, alloc_err
      integer, parameter :: ytag=118
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: buffer
!
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
          partner = ipx + ipy*nprocx + pz*nprocxy
          if (iproc == partner) then
            ! data is local
            out(:,:,pz*bnz+1:(pz+1)*bnz,:) = in
          else
            ! receive from partner
            call MPI_RECV (buffer, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            out(:,:,pz*bnz+1:(pz+1)*bnz,:) = buffer
          endif
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, nbox, MPI_REAL, collector, ytag, MPI_COMM_WORLD, mpierr)
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
      integer :: px, py, pz, collector, partner, nbox, nrow, alloc_err
      integer :: x_add, x_sub, y_add, y_sub
      integer, parameter :: xtag=123, ytag=124
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: buffer, y_row
!
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
      pz = ipz
      if (present (source_pz)) pz = source_pz
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
          partner = ipx + py*nprocx + pz*nprocxy
          y_add = nghost
          y_sub = nghost
          if (py == 0) y_add = 0
          if (py == nprocy-1) y_sub = 0
          if (iproc == partner) then
            ! data is local
            y_row(:,py*cny+1+y_add:py*cny+my,:,:) = in(:,1+y_add:my,:,:)
          else
            ! receive from y-row partner
            call MPI_RECV (buffer, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            y_row(:,py*cny+1+y_add:py*cny+my-y_sub,:,:) = buffer(:,1+y_add:my-y_sub,:,:)
          endif
        enddo
!
        deallocate (buffer)
!
        if (iproc /= collector) then
          ! send to collector
          call MPI_SEND (y_row, nrow, MPI_REAL, collector, xtag, MPI_COMM_WORLD, mpierr)
          deallocate (y_row)
        endif
!
      elseif (ipz == pz) then
        ! send to collector of the y-row (lfirst_proc_y)
        partner = ipx + ipz*nprocxy
        call MPI_SEND (in, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
      endif
!
      if (iproc == collector) then
        ! collect the y-rows into global data
        allocate (buffer(bnx,rny,bnz,bna), stat=alloc_err)
        if (alloc_err > 0) call stop_fatal ('globalize_xy: not enough memory for buffer!', .true.)
!
        do px = 0, nprocx-1
          partner = px + ipy*nprocx + pz*nprocxy
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
            call MPI_RECV (buffer, nrow, MPI_REAL, partner, xtag, MPI_COMM_WORLD, stat, mpierr)
            out(px*cnx+1+x_add:px*cnx+mx-x_sub,:,:,:) = buffer(1+x_add:mx-x_sub,:,:,:)
          endif
        enddo
!
        deallocate (buffer)
      endif
!
    endsubroutine globalize_xy
!***********************************************************************
    subroutine localize_xy(out, in, source_proc, dest_pz)
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
      integer :: rny ! y-row box size
      integer :: px, py, pz, broadcaster, partner, nbox, nrow, alloc_err
      integer, parameter :: xtag=125, ytag=126
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: y_row
!
!
      bnx = size (out, 1)
      bny = size (out, 2)
      bnz = size (out, 3)
      bna = size (out, 4)
      nbox = bnx*bny*bnz*bna
      cnx = bnx - 2*nghost
      cny = bny - 2*nghost
      rny = cny*nprocy + 2*nghost
      nrow = bnx*rny*bnz*bna
!
      broadcaster = ipz * nprocxy
      if (present (source_proc)) broadcaster = broadcaster + source_proc
      pz = ipz
      if (present (dest_pz)) pz = dest_pz
!
      if (iproc == broadcaster) then
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
        if (alloc_err > 0) call stop_fatal ('globalize_xy: not enough memory for y_row!', .true.)
      endif
!
      if (iproc == broadcaster) then
        ! distribute the y-rows
        do px = 0, nprocx-1
          partner = px + pz*nprocxy
          if (iproc == partner) then
            ! data is local
            y_row = in(px*cnx+1:px*cnx+mx,:,:,:)
          else
            ! send to partner
            call MPI_SEND (in(px*cnx+1:px*cnx+mx,:,:,:), &
                nrow, MPI_REAL, partner, xtag, MPI_COMM_WORLD, mpierr)
          endif
        enddo
      endif
!
      if (lfirst_proc_y .and. (ipz == pz)) then
        if (iproc /= broadcaster) then
          ! receive y-row from broadcaster
          call MPI_RECV (y_row, nrow, MPI_REAL, broadcaster, xtag, MPI_COMM_WORLD, stat, mpierr)
        endif
!
        ! distribute the data along the y-direction
        do py = 0, nprocy-1
          partner = ipx + py*nprocx + pz*nprocxy
          if (iproc == partner) then
            ! data is local
            out = y_row(:,py*cny+1:py*cny+my,:,:)
          else
            ! send to partner
            call MPI_SEND (y_row(:,py*cny+1:py*cny+my,:,:), &
                nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
        enddo
      elseif (ipz == pz) then
        ! receive local data from y-row partner (lfirst_proc_y)
        partner = ipx + ipz*nprocxy
        call MPI_RECV (out, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
      endif
!
      if (ipz == pz) deallocate (y_row)
!
    endsubroutine localize_xy
!***********************************************************************
    subroutine globalize_z(in, out, dest_proc)
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
          partner = ipx + ipy*nprocx + pz*nprocxy
          z_add = nghost
          if (pz == 0) z_add = 0
          if (iproc == partner) then
            ! data is local
            out(pz*nz+1+z_add:pz*nz+mz) = in(1+z_add:mz)
          else
            ! receive from partner
            call MPI_RECV (buffer, mz, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            out(pz*nz+1+z_add:pz*nz+mz) = buffer(1+z_add:mz)
          endif
        enddo
!
        deallocate (buffer)
      else
        ! send to collector
        call MPI_SEND (in, mz, MPI_REAL, collector, ytag, MPI_COMM_WORLD, mpierr)
      endif
!
    endsubroutine globalize_z
!***********************************************************************
    subroutine localize_z(out, in, source_proc)
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
      real, dimension(mzgrid), intent(in), optional :: in
      integer, intent(in), optional :: source_proc
!
      integer :: pz, broadcaster, partner
      integer, parameter :: ytag=120
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      broadcaster = ipx + ipy * nprocx
      if (present (source_proc)) broadcaster = broadcaster + source_proc * nprocxy
!
      if (iproc == broadcaster) then
        ! collect the data
        do pz = 0, nprocz-1
          partner = ipx + ipy*nprocx + pz*nprocxy
          if (iproc == partner) then
            ! data is local
            out = in(pz*nz+1:pz*nz+mz)
          else
            ! send to partner
            call MPI_SEND (in(pz*nz+1:pz*nz+mz), mz, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (out, mz, MPI_REAL, broadcaster, ytag, MPI_COMM_WORLD, stat, mpierr)
      endif
!
    endsubroutine localize_z
!***********************************************************************
    subroutine distribute_to_pencil_xy_2D(in, out, broadcaster)
!
!  Distribute data to several processors and reform into pencil shape.
!  This routine divides global data and distributes it in the xy-plane.
!
!  22-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
      integer, intent(in) :: broadcaster
!
      integer :: bnx, bny ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ytag=113
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: buffer
!
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
          partner = ipz*nprocxy + ipy*nprocx + ibox
          if (iproc == partner) then
            ! data is local
            out = in(:,bny*ibox+1:bny*(ibox+1))
          else
            ! send to partner
            buffer = in(:,bny*ibox+1:bny*(ibox+1))
            call MPI_SEND (buffer, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
        enddo
      else
        ! receive from broadcaster
        call MPI_RECV (buffer, nbox, MPI_REAL, broadcaster, ytag, MPI_COMM_WORLD, stat, mpierr)
        out = buffer
      endif
!
      deallocate (buffer)
!
    endsubroutine distribute_to_pencil_xy_2D
!***********************************************************************
    subroutine collect_from_pencil_xy_2D(in, out, collector)
!
!  Collect 2D data from several processors and combine into global shape.
!  This routine collects 2D pencil shaped data distributed in the xy-plane.
!
!  22-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
      integer, intent(in) :: collector
!
      integer :: bnx, bny ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ytag=114
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: buffer
!
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
          partner = ipz*nprocxy + ipy*nprocx + ibox
          if (iproc == partner) then
            ! data is local
            out(:,bny*ibox+1:bny*(ibox+1)) = in
          else
            ! receive from partner
            call MPI_RECV (buffer, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            out(:,bny*ibox+1:bny*(ibox+1)) = buffer
          endif
        enddo
      else
        ! send to collector
        buffer = in
        call MPI_SEND (buffer, nbox, MPI_REAL, collector, ytag, MPI_COMM_WORLD, mpierr)
      endif
!
      deallocate (buffer)
!
    endsubroutine collect_from_pencil_xy_2D
!***********************************************************************
    subroutine remap_to_pencil_x(in, out)
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
      integer, parameter :: xtag=104
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(nx) :: recv_buf
!
!
      do ibox = 0, nprocx-1
        partner = ipz*nprocxy + ipy*nprocx + ibox
        if (iproc == partner) then
          ! data is local
          out(nx*ibox+1:nx*(ibox+1)) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nx, MPI_REAL, partner, xtag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nx, MPI_REAL, partner, xtag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nx, MPI_REAL, partner, xtag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (in, nx, MPI_REAL, partner, xtag, MPI_COMM_WORLD, mpierr)
          endif
          out(nx*ibox+1:nx*(ibox+1)) = recv_buf
        endif
      enddo
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
!
      out = in(nx*ipx+1:nx*(ipx+1))
!
    endsubroutine unmap_from_pencil_x
!***********************************************************************
    subroutine remap_to_pencil_y_1D(in, out)
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
      integer, parameter :: ytag=105
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(ny) :: recv_buf
!
!
      do ibox = 0, nprocy-1
        partner = ipz*nprocxy + ibox*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(ny*ibox+1:ny*(ibox+1)) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, ny, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, ny, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, ny, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (in, ny, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(ny*ibox+1:ny*(ibox+1)) = recv_buf
        endif
      enddo
!
    endsubroutine remap_to_pencil_y_1D
!***********************************************************************
    subroutine remap_to_pencil_y_2D(in, out)
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
      integer, parameter :: ytag=105
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(nx,ny) :: recv_buf
!
!
      nbox = nx*ny
!
      do ibox = 0, nprocy-1
        partner = ipz*nprocxy + ibox*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(:,ny*ibox+1:ny*(ibox+1)) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (in, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,ny*ibox+1:ny*(ibox+1)) = recv_buf
        endif
      enddo
!
    endsubroutine remap_to_pencil_y_2D
!***********************************************************************
    subroutine remap_to_pencil_y_3D(in, out)
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
      integer, parameter :: ytag=105
      integer :: inx, inz ! size of the first and third dimension
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(:,:,:), allocatable :: recv_buf
!
!
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
        partner = ipz*nprocxy + ibox*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(:,ny*ibox+1:ny*(ibox+1),:) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (in, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,ny*ibox+1:ny*(ibox+1),:) = recv_buf
        endif
      enddo
!
      if (allocated (recv_buf)) deallocate (recv_buf)
!
    endsubroutine remap_to_pencil_y_3D
!***********************************************************************
    subroutine remap_to_pencil_y_4D(in, out)
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
      integer, parameter :: ytag=105
      integer :: inx, inz, ina ! size of the first, third, and fourth dimension
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(:,:,:,:), allocatable :: recv_buf
!
!
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
        partner = ipz*nprocxy + ibox*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(:,ny*ibox+1:ny*(ibox+1),:,:) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (in, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,ny*ibox+1:ny*(ibox+1),:,:) = recv_buf
        endif
      enddo
!
      deallocate (recv_buf)
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
!
      out = in(ny*ipy+1:ny*(ipy+1))
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
!
      out = in(:,ny*ipy+1:ny*(ipy+1))
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
!
      out = in(:,ny*ipy+1:ny*(ipy+1),:)
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
!
      out = in(:,ny*ipy+1:ny*(ipy+1),:,:)
!
    endsubroutine unmap_from_pencil_y_4D
!***********************************************************************
    subroutine remap_to_pencil_z_1D(in, out)
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
      integer, parameter :: ztag=106
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(nz) :: recv_buf
!
!
      do ibox = 0, nprocz-1
        partner = ibox*nprocxy + ipy*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(nz*ibox+1:nz*(ibox+1)) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nz, MPI_REAL, partner, ztag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nz, MPI_REAL, partner, ztag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nz, MPI_REAL, partner, ztag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (in, nz, MPI_REAL, partner, ztag, MPI_COMM_WORLD, mpierr)
          endif
          out(nz*ibox+1:nz*(ibox+1)) = recv_buf
        endif
      enddo
!
    endsubroutine remap_to_pencil_z_1D
!***********************************************************************
    subroutine remap_to_pencil_z_2D(in, out)
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
      integer, parameter :: ztag=106
      integer :: ina ! size of the second dimension
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(:,:), allocatable :: recv_buf
!
!
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
        partner = ibox*nprocxy + ipy*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(nz*ibox+1:nz*(ibox+1),:) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (in, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, mpierr)
          endif
          out(nz*ibox+1:nz*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (recv_buf)
!
    endsubroutine remap_to_pencil_z_2D
!***********************************************************************
    subroutine remap_to_pencil_z_3D(in, out)
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
      integer, parameter :: ztag=105
      integer :: inx, iny ! size of the first and third dimension
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(:,:,:), allocatable :: recv_buf
!
!
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
        partner = ibox*nprocxy + ipy*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(:,:,nz*ibox+1:nz*(ibox+1)) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (in, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,:,nz*ibox+1:nz*(ibox+1)) = recv_buf
        endif
      enddo
!
      deallocate (recv_buf)
!
    endsubroutine remap_to_pencil_z_3D
!***********************************************************************
    subroutine remap_to_pencil_z_4D(in, out)
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
      integer, parameter :: ztag=106
      integer :: inx, iny, ina ! size of the first, second, and fourth dimension
      integer, dimension(MPI_STATUS_SIZE) :: stat
      real, dimension(:,:,:,:), allocatable :: recv_buf
!
!
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
        partner = ibox*nprocxy + ipy*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(:,:,nz*ibox+1:nz*(ibox+1),:) = in
        else
          ! communicate with partner
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (in, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (in, nbox, MPI_REAL, partner, ztag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,:,nz*ibox+1:nz*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (recv_buf)
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
!
      out = in(nz*ipz+1:nz*(ipz+1))
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
!
      out = in(nz*ipz+1:nz*(ipz+1),:)
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
!
      out = in(:,:,nz*ipz+1:nz*(ipz+1))
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
!
      out = in(:,:,nz*ipz+1:nz*(ipz+1),:)
!
    endsubroutine unmap_from_pencil_z_4D
!***********************************************************************
    subroutine remap_to_pencil_xy_2D(in, out)
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
      integer, parameter :: ytag=105
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: send_buf, recv_buf
!
!
      if (nprocx == 1) then
        out = in
        return
      endif
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
        partner = ipz*nprocxy + ipy*nprocx + ibox
        if (iproc == partner) then
          ! data is local
          out(bnx*ibox+1:bnx*(ibox+1),:) = in(:,bny*ibox+1:bny*(ibox+1))
        else
          ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1))
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine remap_to_pencil_xy_2D
!***********************************************************************
    subroutine remap_to_pencil_xy_3D(in, out)
!
!  Remaps data distributed on several processors into pencil shape.
!  This routine remaps 3D arrays in x and y only for nprocx>1.
!
!  14-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
!
      integer, parameter :: inx=nx, iny=ny
      integer, parameter :: onx=nxgrid, ony=ny/nprocx
      integer :: inz, onz ! sizes of in and out arrays
      integer, parameter :: bnx=nx, bny=ny/nprocx ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ytag=105
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:), allocatable :: send_buf, recv_buf
!
!
      if (nprocx == 1) then
        out = in
        return
      endif
!
      inz = size (in, 3)
      onz = size (out, 3)
      nbox = bnx*bny*onz
!
      if (mod (ny, nprocx) /= 0) &
          call stop_fatal ('remap_to_pencil_xy_3D: ny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      if ((size (in, 1) /= inx) .or. ((size (in, 2) /= iny))) &
          call stop_fatal ('remap_to_pencil_xy_3D: input array size mismatch /= nx,ny', lfirst_proc_xy)
      if ((size (out, 1) /= onx) .or. ((size (out, 2) /= ony))) &
          call stop_fatal ('remap_to_pencil_xy_3D: output array size mismatch /= nxgrid,ny/nprocx', lfirst_proc_xy)
      if (inz /= onz) &
          call stop_fatal ('remap_to_pencil_xy_3D: inz/=onz - sizes differ in the z direction', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny,onz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_xy_3D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny,onz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('remap_to_pencil_xy_3D: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocx-1
        partner = ipz*nprocxy + ipy*nprocx + ibox
        if (iproc == partner) then
          ! data is local
          out(bnx*ibox+1:bnx*(ibox+1),:,:) = in(:,bny*ibox+1:bny*(ibox+1),:)
        else
          ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1),:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:,:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine remap_to_pencil_xy_3D
!***********************************************************************
    subroutine remap_to_pencil_xy_4D(in, out)
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
      integer, parameter :: ytag=105
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: send_buf, recv_buf
!
!
      if (nprocx == 1) then
        out = in
        return
      endif
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
        partner = ipz*nprocxy + ipy*nprocx + ibox
        if (iproc == partner) then
          ! data is local
          out(bnx*ibox+1:bnx*(ibox+1),:,:,:) = in(:,bny*ibox+1:bny*(ibox+1),:,:)
        else
          ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1),:,:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:,:,:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine remap_to_pencil_xy_4D
!***********************************************************************
    subroutine unmap_from_pencil_xy_2D(in, out)
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
      integer, parameter :: ytag=106
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: send_buf, recv_buf
!
!
      if (nprocx == 1) then
        out = in
        return
      endif
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
        partner = ipz*nprocxy + ipy*nprocx + ibox
        if (iproc == partner) then
          ! data is local
          out(:,bny*ibox+1:bny*(ibox+1)) = in(bnx*ibox+1:bnx*(ibox+1),:)
        else
          ! communicate with partner
          send_buf = in(bnx*ibox+1:bnx*(ibox+1),:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,bny*ibox+1:bny*(ibox+1)) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine unmap_from_pencil_xy_2D
!***********************************************************************
    subroutine unmap_from_pencil_xy_3D(in, out)
!
!  Unmaps pencil shaped 3D data distributed on several processors back to normal shape.
!  This routine is the inverse of the remap function for nprocx>1.
!
!  14-jul-2010/Bourdin.KIS: coded
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
!
      integer, parameter :: inx=nxgrid, iny=ny/nprocx
      integer, parameter :: onx=nx, ony=ny
      integer :: inz, onz ! sizes of in and out arrays
      integer, parameter :: bnx=nx, bny=ny/nprocx ! transfer box sizes
      integer :: ibox, partner, nbox, alloc_err
      integer, parameter :: ytag=106
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:), allocatable :: send_buf, recv_buf
!
!
      if (nprocx == 1) then
        out = in
        return
      endif
!
      inz = size (in, 3)
      onz = size (out, 3)
      nbox = bnx*bny*onz
!
      if (mod (ny, nprocx) /= 0) &
          call stop_fatal ('unmap_from_pencil_xy_3D: ny needs to be an integer multiple of nprocx', lfirst_proc_xy)
!
      if ((size (in, 1) /= inx) .or. ((size (in, 2) /= iny))) &
          call stop_fatal ('unmap_from_pencil_xy_3D: input array size mismatch /= nxgrid,ny/nprocx', lfirst_proc_xy)
      if ((size (out, 1) /= onx) .or. ((size (out, 2) /= ony))) &
          call stop_fatal ('unmap_from_pencil_xy_3D: output array size mismatch /= nx,ny', lfirst_proc_xy)
      if (inz /= onz) &
          call stop_fatal ('unmap_from_pencil_xy_3D: inz/=onz - sizes differ in the z direction', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny,onz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_xy_3D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny,onz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('unmap_from_pencil_xy_3D: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocx-1
        partner = ipz*nprocxy + ipy*nprocx + ibox
        if (iproc == partner) then
          ! data is local
          out(:,bny*ibox+1:bny*(ibox+1),:) = in(bnx*ibox+1:bnx*(ibox+1),:,:)
        else
          ! communicate with partner
          send_buf = in(bnx*ibox+1:bnx*(ibox+1),:,:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,bny*ibox+1:bny*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine unmap_from_pencil_xy_3D
!***********************************************************************
    subroutine unmap_from_pencil_xy_4D(in, out)
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
      integer, parameter :: ytag=106
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: send_buf, recv_buf
!
!
      if (nprocx == 1) then
        out = in
        return
      endif
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
        partner = ipz*nprocxy + ipy*nprocx + ibox
        if (iproc == partner) then
          ! data is local
          out(:,bny*ibox+1:bny*(ibox+1),:,:) = in(bnx*ibox+1:bnx*(ibox+1),:,:,:)
        else
          ! communicate with partner
          send_buf = in(bnx*ibox+1:bnx*(ibox+1),:,:,:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,bny*ibox+1:bny*(ibox+1),:,:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine unmap_from_pencil_xy_4D
!***********************************************************************
    subroutine transp_pencil_xy_2D(in, out)
!
!  Transpose 2D data distributed on several processors.
!  This routine transposes arrays in x and y only.
!  The data must be mapped in pencil shape, especially for nprocx>1.
!
!   4-jul-2010/Bourdin.KIS: coded, adapted parts of transp_xy
!
      real, dimension(:,:), intent(in) :: in
      real, dimension(:,:), intent(out) :: out
!
      integer :: inx, iny, onx, ony ! sizes of in and out arrays
      integer :: bnx, bny, nbox ! destination box sizes and number of elements
      integer :: ibox, partner, alloc_err
      integer, parameter :: ytag=109
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:), allocatable :: send_buf, recv_buf
!
!
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
      do ibox = 0, nprocxy-1
        partner = ipz*nprocxy + ibox
        if (iproc == partner) then
          ! data is local
          out(bnx*ibox+1:bnx*(ibox+1),:) = transpose (in(bny*ibox+1:bny*(ibox+1),:))
        else
          ! communicate with partner
          send_buf = transpose (in(bny*ibox+1:bny*(ibox+1),:))
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine transp_pencil_xy_2D
!***********************************************************************
    subroutine transp_pencil_xy_3D(in, out)
!
!  Transpose 3D data distributed on several processors.
!  This routine transposes arrays in x and y only.
!  The data must be mapped in pencil shape, especially for nprocx>1.
!
!  14-jul-2010/Bourdin.KIS: coded, adapted parts of transp_xy
!
      real, dimension(:,:,:), intent(in) :: in
      real, dimension(:,:,:), intent(out) :: out
!
      integer :: inx, iny, inz, onx, ony, onz ! sizes of in and out arrays
      integer :: bnx, bny, nbox ! destination box sizes and number of elements
      integer :: ibox, partner, alloc_err, pos_z
      integer, parameter :: ytag=109
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:), allocatable :: send_buf, recv_buf
!
!
      inx = size (in, 1)
      iny = size (in, 2)
      inz = size (in, 3)
      onx = size (out, 1)
      ony = size (out, 2)
      onz = size (out, 3)
      bnx = onx/nprocxy
      bny = ony
      nbox = bnx*bny*onz
!
      if (mod (onx, nprocxy) /= 0) &
          call stop_fatal ('transp_pencil_xy_3D: onx needs to be an integer multiple of nprocxy', lfirst_proc_xy)
!
      if ((inx /= bny*nprocxy) .or. (iny /= bnx)) &
          call stop_fatal ('transp_pencil_xy_3D: input array has unmatching size', lfirst_proc_xy)
      if ((onx /= bnx*nprocxy) .or. (ony /= bny)) &
          call stop_fatal ('transp_pencil_xy_3D: output array has unmatching size', lfirst_proc_xy)
      if (inz /= onz) &
          call stop_fatal ('transp_pencil_xy_3D: inz/=onz - sizes differ in the z direction', lfirst_proc_xy)
!
      allocate (send_buf(bnx,bny,onz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('transp_pencil_xy_3D: not enough memory for send_buf!', .true.)
      allocate (recv_buf(bnx,bny,onz), stat=alloc_err)
      if (alloc_err > 0) call stop_fatal ('transp_pencil_xy_3D: not enough memory for recv_buf!', .true.)
!
      do ibox = 0, nprocxy-1
        partner = ipz*nprocxy + ibox
        if (iproc == partner) then
          ! data is local
          do pos_z = 1, onz
            out(bnx*ibox+1:bnx*(ibox+1),:,pos_z) = transpose (in(bny*ibox+1:bny*(ibox+1),:,pos_z))
          enddo
        else
          ! communicate with partner
          do pos_z = 1, onz
            send_buf(:,:,pos_z) = transpose (in(bny*ibox+1:bny*(ibox+1),:,pos_z))
          enddo
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:,:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine transp_pencil_xy_3D
!***********************************************************************
    subroutine transp_pencil_xy_4D(in, out)
!
!  Transpose 4D data distributed on several processors.
!  This routine transposes arrays in x and y only.
!  The data must be mapped in pencil shape, especially for nprocx>1.
!
!  14-jul-2010/Bourdin.KIS: coded, adapted parts of transp_xy
!
      real, dimension(:,:,:,:), intent(in) :: in
      real, dimension(:,:,:,:), intent(out) :: out
!
      integer :: inx, iny, inz, ina, onx, ony, onz, ona ! sizes of in and out arrays
      integer :: bnx, bny, nbox ! destination box sizes and number of elements
      integer :: ibox, partner, alloc_err, pos_z, pos_a
      integer, parameter :: ytag=109
      integer, dimension(MPI_STATUS_SIZE) :: stat
!
      real, dimension(:,:,:,:), allocatable :: send_buf, recv_buf
!
!
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
      do ibox = 0, nprocxy-1
        partner = ipz*nprocxy + ibox
        if (iproc == partner) then
          ! data is local
          do pos_z = 1, onz
            do pos_a = 1, ona
              out(bnx*ibox+1:bnx*(ibox+1),:,pos_z,pos_a) = transpose (in(bny*ibox+1:bny*(ibox+1),:,pos_z,pos_a))
            enddo
          enddo
        else
          ! communicate with partner
          do pos_z = 1, onz
            do pos_a = 1, ona
              send_buf(:,:,pos_z,pos_a) = transpose (in(bny*ibox+1:bny*(ibox+1),:,pos_z,pos_a))
            enddo
          enddo
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(bnx*ibox+1:bnx*(ibox+1),:,:,:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine transp_pencil_xy_4D
!***********************************************************************
    subroutine remap_to_pencil_yz_3D(in, out)
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
!
!
      if (nprocz == 1) then
        out = in
        return
      endif
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
        partner = ibox*nprocxy + ipy*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(:,:,bnz*ibox+1:bnz*(ibox+1)) = in(:,bny*ibox+1:bny*(ibox+1),:)
        else
          ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1),:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,:,bnz*ibox+1:bnz*(ibox+1)) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine remap_to_pencil_yz_3D
!***********************************************************************
    subroutine remap_to_pencil_yz_4D(in, out)
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
!
!
      if (nprocz == 1) then
        out = in
        return
      endif
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
        partner = ibox*nprocxy + ipy*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(:,:,bnz*ibox+1:bnz*(ibox+1),:) = in(:,bny*ibox+1:bny*(ibox+1),:,:)
        else
          ! communicate with partner
          send_buf = in(:,bny*ibox+1:bny*(ibox+1),:,:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,:,bnz*ibox+1:bnz*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine remap_to_pencil_yz_4D
!***********************************************************************
    subroutine unmap_from_pencil_yz_3D(in, out)
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
!
!
      if (nprocz == 1) then
        out = in
        return
      endif
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
        partner = ibox*nprocxy + ipy*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(:,bny*ibox+1:bny*(ibox+1),:) = in(:,:,bnz*ibox+1:bnz*(ibox+1))
        else
          ! communicate with partner
          send_buf = in(:,:,bnz*ibox+1:bnz*(ibox+1))
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,bny*ibox+1:bny*(ibox+1),:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine unmap_from_pencil_yz_3D
!***********************************************************************
    subroutine unmap_from_pencil_yz_4D(in, out)
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
!
!
      if (nprocz == 1) then
        out = in
        return
      endif
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
      do ibox = 0, nprocx-1
        partner = ibox*nprocxy + ipy*nprocx + ipx
        if (iproc == partner) then
          ! data is local
          out(:,bny*ibox+1:bny*(ibox+1),:,:) = in(:,:,bnz*ibox+1:bnz*(ibox+1),:)
        else
          ! communicate with partner
          send_buf = in(:,:,bnz*ibox+1:bnz*(ibox+1),:)
          if (iproc > partner) then ! above diagonal: send first, receive then
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
          else                      ! below diagonal: receive first, send then
            call MPI_RECV (recv_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, stat, mpierr)
            call MPI_SEND (send_buf, nbox, MPI_REAL, partner, ytag, MPI_COMM_WORLD, mpierr)
          endif
          out(:,bny*ibox+1:bny*(ibox+1),:,:) = recv_buf
        endif
      enddo
!
      deallocate (send_buf, recv_buf)
!
    endsubroutine unmap_from_pencil_yz_4D
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
      ay=0.
      ay_local=0.
      if (ipz==(zproc_no-1)) then
        my_iniy=ipy*ny+1
        my_finy=(ipy+1)*ny
        ay_local(my_iniy:my_finy)=a(xi,:,zj)
      else
      ay_local=0.
      endif
      call mpireduce_sum(ay_local,ay,nygrid)
! maybe we should synchrosize here.
      call mpibarrier()
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
      az=0.
      az_local=0.
      if (ipy==(yproc_no-1)) then
        my_iniz=ipz*nz+1
        my_finz=(ipz+1)*nz
        az_local(my_iniz:my_finz)=a(xi,yj,:)
      else
      az_local=0.
      endif
      call mpireduce_sum(az_local,az,nzgrid)
! maybe we should synchrosize here.
      call mpibarrier()
!
    endsubroutine z2x
!***********************************************************************
    subroutine parallel_open(unit,file,form)
!
!  Choose between two reading methods.
!
!  19-nov-10/dhruba.mitra: implemented
!
      integer :: unit
      character (len=*) :: file
      character (len=*), optional :: form
!
      if (lfake_parallel_io) then
        call fake_parallel_open(unit,file,form)
      else
        call true_parallel_open(unit,file,form)
      endif
!
    endsubroutine parallel_open
!***********************************************************************
    subroutine fake_parallel_open(unit,file,form)
!
!  Read a global file.
!
!  18-mar-10/Bourdin.KIS: implemented
!
      integer :: unit
      character (len=*) :: file
      character (len=*), optional :: form
!
      logical :: exists
!
!  Test if file exists.
!
      inquire(FILE=file,exist=exists)
      if (.not. exists) call stop_it('fake_parallel_open: file not found "'//trim(file)//'"')
!
!  Open file.
!
      if (present(form)) then
        open(unit, FILE=file, FORM=form, STATUS='old')
      else
        open(unit, FILE=file, STATUS='old')
      endif
!
    endsubroutine fake_parallel_open
!***********************************************************************
    subroutine true_parallel_open(unit,file,form,recl)
!
!  Read a global file in parallel.
!
!  17-mar-10/Bourdin.KIS: implemented
!
      use Cparam, only: fnlen
      use Syscalls, only: file_size, get_tmp_prefix
!
      integer :: unit
      character (len=*) :: file
      character (len=*), optional :: form
      integer, optional :: recl
!
      logical :: exists
      integer :: bytes, pos
      integer, parameter :: buf_len=fnlen
      character (len=buf_len) :: filename
      character, dimension(:), allocatable :: buffer
      character(len=fnlen) :: get_tmp_prefix_
!
      if (lroot) then
!
!  Test if file exists.
!
        inquire(FILE=file,exist=exists)
        if (.not. exists) call stop_it_if_any(.true., &
            'true_parallel_open: file not found "'//trim(file)//'"')
        bytes=file_size(file)
        if (bytes < 0) call stop_it_if_any(.true., &
            'true_parallel_open: could not determine file size "'//trim(file)//'"')
        if (bytes == 0) call stop_it_if_any(.true., &
            'true_parallel_open: file is empty "'//trim(file)//'"')
      endif
!
!  Catch conditional errors of the MPI root rank.
!
      call stop_it_if_any(.false.,'')
!
!  Broadcast the file size.
!
      call mpibcast_int(bytes, 1)
!
!  Allocate temporary memory.
!
      allocate(buffer(bytes))
      buffer=char(0)
!
      if (lroot) then
!
!  Read file content into buffer.
!
        open(unit, FILE=file, FORM='unformatted', RECL=bytes, ACCESS='direct', STATUS='old')
        read(unit, REC=1) buffer
        close(unit)
      endif
!
!  Broadcast buffer to all MPI ranks.
!
      call mpibcast_char(buffer, bytes)
!
!  Create unique temporary filename.
!
      pos=scan(file, '/')
      do while(pos /= 0)
        file(pos:pos)='_'
        pos=scan(file, '/')
      enddo
      get_tmp_prefix_=get_tmp_prefix()
      write(filename,'(A,A,A,I0)') trim(get_tmp_prefix_), file, '-', iproc
!
!  Write temporary file into local RAM disk (/tmp).
!
!     *** WORK HERE:
!     *** Instead of writing N records with one byte, we want to write 1 record with N bytes...
!     *** Problem is that the record length is not neccessarily in bytes.
!     *** Solution is => we need to devide by the number of bytes per record.
!     *** But, how to determine the number of bytes per record?
!      open(unit, FILE=filename, FORM='unformatted', RECL=bytes/bpr, ACCESS='direct', STATUS='replace')
!      write(unit, REC=1) buffer
!     *** WORK HERE: TEMPORARY REPLACEMENT CODE:
      open(unit, FILE=filename, FORM='formatted', RECL=1, ACCESS='direct', STATUS='replace')
      do pos=1,bytes
        write (unit, '(A)', REC=pos) buffer(pos)
      enddo
      close(unit)
      deallocate(buffer)
!
!  Open temporary file.
!
      if (present(form) .and. present(recl)) then
        open(unit, FILE=filename, FORM=form, RECL=recl, STATUS='old')
      elseif (present(recl)) then
        open(unit, FILE=filename, RECL=recl, STATUS='old')
      elseif (present(form)) then
        open(unit, FILE=filename, FORM=form, STATUS='old')
      else
        open(unit, FILE=filename, STATUS='old')
      endif
!
!  Unit is now reading from RAM and is ready to be used on all ranks in
!  parallel.
!
    endsubroutine true_parallel_open
!***********************************************************************
    subroutine parallel_close(unit)
!
!  Close a file unit opened by parallel_open and remove temporary file.
!
!  17-mar-10/Bourdin.KIS: implemented
!
!
      integer :: unit
!
      if (lfake_parallel_io) then
        call fake_parallel_close(unit)
      else
         call true_parallel_close(unit)
       endif
!
    endsubroutine parallel_close
!***********************************************************************
    subroutine fake_parallel_close(unit)
!
!  Close a file unit opened by parallel_open and remove temporary file.
!
!  17-mar-10/Bourdin.KIS: implemented
!
      integer :: unit
!
      close(unit)
!
    endsubroutine fake_parallel_close
!***********************************************************************
    subroutine true_parallel_close(unit)
!
!  Close a file unit opened by parallel_open and remove temporary file.
!
!  17-mar-10/Bourdin.KIS: implemented
!
      integer :: unit
!
      close(unit,STATUS='delete')
!
    endsubroutine true_parallel_close
!***********************************************************************
    function parallel_count_lines(file)
!
!  Determines in parallel the number of lines in a file.
!
!  Returns:
!  * Integer containing the number of lines in a given file
!  * -1 on error
!
!  23-mar-10/Bourdin.KIS: implemented
!
      use Syscalls, only: count_lines
!
      character(len=*) :: file
      integer :: parallel_count_lines
!
      if (lroot) parallel_count_lines = count_lines(file)
      call mpibcast_int(parallel_count_lines, 1)
!
    endfunction
!***********************************************************************
    function parallel_file_exists(file, delete)
!
!  Determines in parallel if a given file exists.
!  If delete is true, deletes the file.
!
!  Returns:
!  * Integer containing the number of lines in a given file
!  * -1 on error
!
!  23-mar-10/Bourdin.KIS: implemented
!
      use Syscalls, only: file_exists
!
      character(len=*) :: file
      logical :: parallel_file_exists,ldelete
      logical, optional :: delete
!
      if (present(delete)) then
        ldelete=delete
      else
        ldelete=.false.
      endif
!
      ! Let the root node do the dirty work
      if (lroot) parallel_file_exists = file_exists(file,ldelete)
!
      call mpibcast_logical(parallel_file_exists, 1)
!
    endfunction
!***********************************************************************
    subroutine mpigather_xy( sendbuf, recvbuf, lpz )
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
!
      if (present(lpz)) then
        cond = ipz==lpz
      else
        cond = .true.
      endif
!
      ncnt = nxgrid*ny
!
      if (cond) &
        call MPI_GATHER(sendbuf, ncnt, MPI_REAL, recvbuf, ncnt, MPI_REAL, root, MPI_COMM_XYPLANE, mpierr)
!
    endsubroutine mpigather_xy
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
      integer lpx, lpy
!
      if (present(lproc)) then
        lpy = lproc/nprocx
        lpx = mod(lproc,nprocx)
      else
        lpy=0; lpx=0
      endif
!
      if ( ipx==lpx .and. ipy==lpy ) &
        call MPI_GATHER(sendbuf, n1*nz, MPI_REAL, recvbuf, n1*nz, MPI_REAL, root, MPI_COMM_ZBEAM, mpierr)
!
    endsubroutine mpigather_z
!***********************************************************************
    subroutine mpigather( sendbuf, recvbuf )
!
!  Gathers the chunks of a 3D array from each processor in a big array at root.
!
!  Here no parallelization in x allowed.
!
!  19-nov-10/MR: coded
!
      real, dimension(nxgrid,ny,nz)         :: sendbuf   ! nx=nxgrid !
      real, dimension(nxgrid,nygrid,nzgrid) :: recvbuf
!
      integer :: ncnt, nshift, nlayer, i
      integer, dimension(ncpus) :: counts, shifts
!
      ncnt = nxgrid*ny
!
      if (lroot) then
!
        counts = ncnt
        nlayer = nz*nxgrid*nygrid
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
        call MPI_GATHERV(sendbuf(1,1,i), ncnt, MPI_REAL, recvbuf(1,1,i), counts, shifts, &
                     MPI_REAL, root, MPI_COMM_WORLD, mpierr)
      enddo
!
    endsubroutine mpigather
!***********************************************************************
    subroutine mpigather_and_out_real( sendbuf, unit, ltransp, kxrange, kyrange,zrange )
!
!  Transfers the chunks of a 3D array from each processor to root
!  and writes them out in right order.
!
!  Here no parallelization in x allowed.
!
!  22-nov-10/MR: coded
!  06-apr-11/MR: optional parameters kxrange, kyrange, zrange for selective output added
!
      use General, only: write_by_ranges_2d_cmplx, write_by_ranges_2d_real, &
                         write_by_ranges_1d_cmplx, write_by_ranges_1d_real
!
      integer,                               intent(in   ) :: unit, ncomp
      real,    dimension(nxgrid,ny,nz),      intent(inout) :: sendbuf
      complex, dimension(nxgrid,ny,nz,ncomp),intent(inout) :: sendbuf_cmplx
      logical,                      optional,intent(in   ) :: ltransp   ! if true, transposition x <-> y
      integer, dimension(3,*),      optional,intent(in   ) :: kxrange, kyrange, zrange
!
      integer :: i,np,iproca,iproce,nxy,tag
      integer, dimension(3,10) :: kxrangel,kyrangel,zrangel
!
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer :: k,ii,ncompl,ic
      logical :: ltrans, lcomplex
      real,    allocatable :: rowbuf(:)
      complex, allocatable :: rowbuf_cmplx(:)
!
      lcomplex = .false.
      ncompl = 1
      goto 1
!
      entry mpigather_and_out_cmplx( sendbuf_cmplx, ncomp, unit, ltransp, kxrange, kyrange,zrange )
!
      lcomplex = .true.
      ncompl = ncomp
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
        kxrangel = kxrange(:,1:10)
      endif
!
      if ( .not.present(kyrange) ) then
        kyrangel = 0
        kyrangel(:,1)=(/1,nygrid,1/)
      else
        kyrangel = kyrange(:,1:10)
      endif
!
      if ( .not.present(zrange) ) then
        zrangel = 0
        zrangel(:,1) = (/1,nzgrid,1/)
      else
        zrangel = zrange(:,1:10)
      endif
!
      if (ltrans) then
        if (lcomplex) then
          allocate( rowbuf(ny) )
        else
          allocate( rowbuf_cmplx(ny) )
        endif
      endif
!
      nxy = nxgrid*ny
!
      iproca=1; iproce=-1
!
      do ic=1,ncompl
      do np=1,nprocz
!
        iproce = iproce + nprocxy
!
        do k=1,10
          if ( zrangel(1,k) > 0 ) then
            do iz = zrangel(1,k), zrangel(2,k), zrangel(3,k)
              if ( iz >= n1 .and. iz <= n2 ) then
!
                if (ltrans) then
!
                  do ii=1,10
                    if ( kyrangel(1,ii) > 0 ) then
                      do iy = kyrangel(1,ii), kyrangel(2,ii), kyrangel(3,ii)
!
                        if ( lroot .and. np==1 ) then
                          if (lcomplex) then
                            call write_by_ranges_1d_cmplx( 1, sendbuf_cmplx(iy,1,iz,ic), kxrangel )
                          else
                            call write_by_ranges_1d_real( 1, sendbuf(iy,1,iz), kxrangel )
                          endif
                        endif
!
                        do i=iproca,iproce
!
!  overflow possible for large ncpuxy, nz, nxgrid (comment by Anders)
!
                          tag = nprocxy*(iz+1)*iy + nprocxy*iz + i-iproca
!
                          if (lroot) then
                            if (lcomplex) then
                              call MPI_RECV(rowbuf_cmplx, ny, MPI_COMPLEX, i, tag, MPI_COMM_WORLD, status, mpierr)
                              call write_by_ranges_1d_cmplx( 1, rowbuf_cmplx, kxrangel )
                            else
                              call MPI_RECV(rowbuf, ny, MPI_REAL, i, tag, MPI_COMM_WORLD, status, mpierr)
                              call write_by_ranges_1d_real( 1, rowbuf, kxrangel )
                            endif
                          else if ( iproc==i ) then
                            if (lcomplex) then
                              rowbuf_cmplx=sendbuf_cmplx(iy,:,iz,ic)
                              call MPI_SEND(rowbuf_cmplx, ny, MPI_COMPLEX, root, tag, MPI_COMM_WORLD, mpierr)
                            else
                              rowbuf=sendbuf(iy,:,iz)
                              call MPI_SEND(rowbuf, ny, MPI_REAL, root, tag, MPI_COMM_WORLD, mpierr)
                            endif
                          endif
!
                        enddo
!
                      enddo
                    endif
                  enddo
!
                  call mpibarrier()         ! necessary ?
!
                else
!
                  if (lroot .and. np==1 ) then
                    if (lcomplex) then
                      call write_by_ranges_2d_cmplx( 1, sendbuf_cmplx(1,1,iz,ic), kxrangel, kyrangel )
                    else
                      call write_by_ranges_2d_real( 1, sendbuf(1,1,iz), kxrangel, kyrangel )
                    endif
                  endif
!
                  do i=iproca,iproce
!
                    tag = nprocxy*iz+i-iproca
                    if (lroot) then
                      if (lcomplex) then
                        call MPI_RECV(sendbuf_cmplx(1,1,iz,ic), nxy, MPI_COMPLEX, i, tag, MPI_COMM_WORLD, status, mpierr)
                        call write_by_ranges_2d_cmplx( 1, sendbuf_cmplx(1,1,iz,ic), kxrangel, kyrangel )
                      else
                        call MPI_RECV(sendbuf(1,1,iz), nxy, MPI_REAL, i, tag, MPI_COMM_WORLD, status, mpierr)
                        call write_by_ranges_2d_real( 1, sendbuf(1,1,iz), kxrangel, kyrangel )
                      endif
                    else if ( iproc==i ) then
                      if (lcomplex) then
                        call MPI_SEND(sendbuf_cmplx(1,1,iz,ic), nxy, MPI_COMPLEX, root, tag, MPI_COMM_WORLD, mpierr)
                      else
                        call MPI_SEND(sendbuf(1,1,iz), nxy, MPI_REAL, root, tag, MPI_COMM_WORLD, mpierr)
                      endif
                    endif
!
                  enddo
                endif
              endif
            enddo
          endif
        enddo
!
        iproca = iproce+1
!
      enddo
      enddo
!
      if (ltrans) deallocate(rowbuf)
!
    contains
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      integer function update_cnt(fcnt,len,ntot,nl)        ! internal function
!
        integer, intent(in ) :: fcnt,len,ntot
        integer, intent(out) :: nl
!
        integer :: nr
!
        if (fcnt>0) then
          nr=ntot-(len-fcnt)
        else
          nr=ntot
        endif
!
        update_cnt=mod(nr,len)
        nl=nr-update_cnt
!
      endfunction update_cnt
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    endsubroutine mpigather_and_out_real
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
    subroutine mpimerge_1d(vector,nk,idir)
!
!  Merges vectors of processors along idir by filling invalid values (NaN).
!
!  22-nov-10/MR: coded
!
      integer,             intent(in)    :: nk
      real, dimension(nk), intent(inout) :: vector
      integer, optional,   intent(in)    :: idir
!
      integer                            :: mpiprocs,merge
      real, dimension(nk)                :: recvbuf
!
      if (present(idir)) then
        mpiprocs=mpigetcomm(idir)
      else
        mpiprocs=MPI_COMM_WORLD
      endif
!
      call MPI_OP_CREATE( merge_1d, .false., merge, mpierr )
      call MPI_REDUCE(vector, recvbuf, nk, MPI_REAL, merge, root, mpiprocs, mpierr)
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
        case(1)
          mpigetcomm=MPI_COMM_XBEAM
        case(2)
          mpigetcomm=MPI_COMM_YBEAM
        case(3)
          mpigetcomm=MPI_COMM_ZBEAM
        case(12)
          mpigetcomm=MPI_COMM_XYPLANE
        case(13)
          mpigetcomm=MPI_COMM_XZPLANE
        case(23)
          mpigetcomm=MPI_COMM_YZPLANE
        case default
          mpigetcomm=MPI_COMM_WORLD
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
      call MPI_GATHER(flag, 1, MPI_LOGICAL, flags, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, mpierr)
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
        deallocate(flags)
!
        if (count>0) then
          call safe_character_prepend(message,'"at '//trim(itoa(count))//' node(s): ')
          report_clean_output = .true.
        endif
!
      endif
!
      call MPI_BCAST(report_clean_output,1,MPI_LOGICAL,root,MPI_COMM_WORLD,mpierr)
! broadcasts flag for
!
! 'sychronization necessary'
    end function report_clean_output
!**************************************************************************
endmodule Mpicomm
