!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!$Id$
!
  private

  public :: remap_to_pencil_xy_2D_other, unmap_from_pencil_xy_2D_other

  public :: update_neighbors, index_to_iproc_comm

  public :: mpicomm_init, initialize_mpicomm, mpifinalize, yyinit
  public :: mpibarrier
  public :: stop_it, stop_it_if_any
  public :: die_gracefully, die_immediately
  public :: check_emergency_brake

  public :: mpirecv_logical, mpirecv_real, mpirecv_int, mpirecv_char  !, mpirecv_cmplx
  public :: mpisend_logical, mpisend_real, mpiscan_int, mpisend_int, mpisend_char   !, mpisend_cmplx
  public :: mpisendrecv_real, mpisendrecv_int
  public :: mpireduce_sum_int, mpireduce_sum             !, mpireduce_sum_double
  public :: mpireduce_max, mpireduce_max_int, mpireduce_min
  public :: mpiallreduce_max, mpiallreduce_min
  public :: mpiallreduce_sum, mpiallreduce_sum_int
  public :: mpiallreduce_sum_arr, mpiallreduce_sum_arr2
  public :: mpiallreduce_or, mpiallreduce_and
  public :: mpireduce_or, mpireduce_and
  public :: mpibcast, mpibcast_real, mpibcast_logical
  public :: mpibcast_real_arr, mpibcast_cmplx_arr_dbl, mpibcast_cmplx
  public :: mpibcast_double
  public :: mpibcast_int, mpibcast_char, mpireduce_max_scl_int
  public :: mpiscatter
  public :: mpigather_scl_str, mpigather_xy, mpimerge_1d, mpigather_z, &
            mpigather_and_out_cmplx, mpigather_and_out_real
  public :: mpiscatterv
  public :: mpiwtime, mpiwtick

  public :: mpisend_nonblock_real,mpisend_nonblock_int
  public :: mpirecv_nonblock_real,mpirecv_nonblock_int

  public :: start_serialize,end_serialize
  public :: initiate_isendrcv_bdry, finalize_isendrcv_bdry
  public :: isendrcv_bdry_x
  public :: initiate_shearing, finalize_shearing

  public :: transp, transp_xy, transp_xy_other, transp_other
  public :: transp_xz, transp_zx

  public :: communicate_vect_field_ghosts, communicate_xy_ghosts
  public :: fill_zghostzones_3vec

  public :: sum_xy, distribute_xy, collect_xy, distribute_yz
  public :: distribute_z, collect_z
  public :: globalize_xy, localize_xy
  public :: globalize_z, localize_z
  public :: distribute_to_pencil_xy, collect_from_pencil_xy
  public :: remap_to_pencil_x, unmap_from_pencil_x
  public :: remap_to_pencil_y, unmap_from_pencil_y
  public :: remap_to_pencil_z, unmap_from_pencil_z
  public :: remap_to_pencil_xy, unmap_from_pencil_xy, transp_pencil_xy
  public :: remap_to_pencil_yz, unmap_from_pencil_yz
  public :: collect_grid
  public :: set_rslice_communicator,root_rslice
  public :: y2x, z2x

  public :: report_clean_output
  public :: mpiwait
  
! Radiation ray routines
  public :: radboundary_xy_recv, radboundary_xy_send
  public :: radboundary_zx_recv, radboundary_zx_send
  public :: radboundary_yz_sendrecv, radboundary_zx_sendrecv
  public :: radboundary_yz_periodic_ray, radboundary_zx_periodic_ray
  public :: radboundary_yz_recv, radboundary_yz_send

! Foreign application routines.
  public :: initialize_foreign_comm, get_foreign_snap_initiate, get_foreign_snap_finalize, update_foreign_data
! Variables
  public :: ipx, ipy, ipz, lroot, iproc, mpi_precision, nprocs
  public :: lfirst_proc_x, lfirst_proc_y, lfirst_proc_z, lfirst_proc_xy, lfirst_proc_yz, lfirst_proc_xz, lfirst_proc_xyz
  public :: llast_proc_x, llast_proc_y, llast_proc_z, llast_proc_xy, llast_proc_yz, llast_proc_xz, llast_proc_xyz
  public :: MPI_COMM_WORLD, MPI_COMM_GRID, MPI_COMM_PENCIL, MPI_COMM_XYPLANE, MPI_COMM_XZPLANE, MPI_COMM_YZPLANE, &
            MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, MPI_COMM_RSLICE, &
            MPI_INFO_NULL, MPI_ANY_TAG, lyang
  public :: size_int, size_real, size_double
!
  interface mpirecv_logical
    module procedure mpirecv_logical_scl
    module procedure mpirecv_logical_arr
  endinterface
!
  interface mpisend_char
    module procedure mpisend_char_scl
  endinterface
!
  interface mpirecv_char
    module procedure mpirecv_char_scl
  endinterface
!
  interface mpirecv_real
    module procedure mpirecv_real_scl
    module procedure mpirecv_real_arr
    module procedure mpirecv_real_arr2
    module procedure mpirecv_real_arr3
    module procedure mpirecv_real_arr4
    module procedure mpirecv_real_arr5
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
    module procedure mpisend_real_arr5
  endinterface
!
  interface mpisendrecv_int
     module procedure mpisendrecv_int_arr
  endinterface

  interface mpisendrecv_real
    module procedure mpisendrecv_real_scl
    module procedure mpisendrecv_real_arr
    module procedure mpisendrecv_real_arr2
    module procedure mpisendrecv_real_arr3
    module procedure mpisendrecv_real_arr4
  endinterface
!
  interface mpisend_int
    module procedure mpisend_int_scl
    module procedure mpisend_int_arr
    module procedure mpisend_int_arr2
  endinterface
!
  interface mpibcast
    module procedure mpibcast_logical_scl
    module procedure mpibcast_logical_arr
    module procedure mpibcast_logical_arr2
    module procedure mpibcast_int_scl
    module procedure mpibcast_int_arr
    module procedure mpibcast_int_arr2
    module procedure mpibcast_real_scl
    module procedure mpibcast_real_arr
    module procedure mpibcast_real_arr2
    module procedure mpibcast_real_arr3
    module procedure mpibcast_real_arr4
    module procedure mpibcast_cmplx_arr_sgl
    module procedure mpibcast_char_scl
    module procedure mpibcast_char_arr
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
    module procedure mpibcast_int_arr2
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
  interface mpibcast_cmplx
    module procedure mpibcast_cmplx_arr_sgl
  endinterface
!
  interface mpibcast_char
    module procedure mpibcast_char_scl
    module procedure mpibcast_char_arr
  endinterface
!
  interface mpiscatter
    module procedure mpiscatter_real_arr
    module procedure mpiscatter_real_arr2
  endinterface
!
  interface mpiallreduce_sum
    module procedure mpiallreduce_sum_scl
    module procedure mpiallreduce_sum_arr
    module procedure mpiallreduce_sum_arr2
    module procedure mpiallreduce_sum_arr3
    module procedure mpiallreduce_sum_arr4
    module procedure mpiallreduce_sum_arr5
    module procedure mpiallreduce_sum_arr_inplace
  endinterface
!
  interface mpiallreduce_sum_int
    module procedure mpiallreduce_sum_int_scl
    module procedure mpiallreduce_sum_int_arr
    module procedure mpiallreduce_sum_int_arr_inplace
  endinterface
!
  interface mpiallreduce_max
    module procedure mpiallreduce_max_scl_dbl
    module procedure mpiallreduce_max_scl_sgl
    module procedure mpiallreduce_max_scl_int
    module procedure mpiallreduce_max_arr
  endinterface
!
  interface mpiallreduce_min
    module procedure mpiallreduce_min_scl_sgl
    module procedure mpiallreduce_min_scl_dbl
    module procedure mpiallreduce_min_scl_int
  endinterface
!
  interface mpiallreduce_or
    module procedure mpiallreduce_or_scl
    module procedure mpiallreduce_or_arr_inplace !!!
  endinterface
!
  interface mpiallreduce_and
    module procedure mpiallreduce_and_scl
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
  interface distribute_xy
    module procedure distribute_xy_0D
    module procedure distribute_xy_2D
    module procedure distribute_xy_3D
    module procedure distribute_xy_4D
  endinterface
!
  interface distribute_yz
    module procedure distribute_yz_3D
    module procedure distribute_yz_4D
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
  interface mpirecv_nonblock_real
    module procedure mpirecv_nonblock_real_arr
    module procedure mpirecv_nonblock_real_arr2  !!!
    module procedure mpirecv_nonblock_real_arr3
    module procedure mpirecv_nonblock_real_arr4
    module procedure mpirecv_nonblock_real_arr5
  endinterface
!
  interface mpisend_nonblock_real
    module procedure mpisend_nonblock_real_arr
    module procedure mpisend_nonblock_real_arr2
    module procedure mpisend_nonblock_real_arr3
    module procedure mpisend_nonblock_real_arr4
    module procedure mpisend_nonblock_real_arr5
  endinterface
!
  interface mpirecv_nonblock_int
    module procedure mpirecv_nonblock_int_scl
    module procedure mpirecv_nonblock_int_arr
    module procedure mpirecv_nonblock_int_arr2
  endinterface
!
  interface mpisend_nonblock_int
    module procedure mpisend_nonblock_int_scl
    module procedure mpisend_nonblock_int_arr
    module procedure mpisend_nonblock_int_arr2
  endinterface
!
  interface mpiscatterv
    module procedure mpiscatterv_real
    module procedure mpiscatterv_int
  endinterface
!
!  Communicators
!
  integer :: MPI_COMM_GRID, MPI_COMM_PENCIL
  integer :: MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM
  integer :: MPI_COMM_XYPLANE,MPI_COMM_XZPLANE,MPI_COMM_YZPLANE,MPI_COMM_RSLICE
  integer :: root_rslice
!
! for protecting MPI_COMM_WORLD to be redefined by preprocessor
! 
  integer, parameter, public :: MPI_COMM_UNIVERSE=MPI_COMM_WORLD
!
! symbolic constants for beams and planes
!
  integer, parameter, public :: IXBEAM=1, IYBEAM=2, IZBEAM=3, IXYPLANE=12, IXZPLANE=13, IYZPLANE=23

  character(LEN=4), public :: cyinyang=' '
!
  integer :: mpi_precision
!
