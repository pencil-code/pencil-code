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
  public :: mpisendrecv_real
  public :: mpireduce_sum_int, mpireduce_sum             !, mpireduce_sum_double
  public :: mpireduce_max, mpireduce_max_int, mpireduce_min
  public :: mpiallreduce_max, mpiallreduce_min
  public :: mpiallreduce_sum, mpiallreduce_sum_int
  public :: mpiallreduce_sum_arr, mpiallreduce_sum_arr2
  public :: mpiallreduce_or
  public :: mpireduce_or, mpireduce_and
  public :: mpibcast, mpibcast_real, mpibcast_logical
  public :: mpibcast_real_arr, mpibcast_cmplx_arr_dbl, mpibcast_cmplx
  public :: mpibcast_double
  public :: mpibcast_int, mpibcast_char, mpireduce_max_scl_int
  public :: mpigather_scl_str, mpigather_xy, mpimerge_1d, mpigather_z, &
            mpigather_and_out_cmplx, mpigather_and_out_real
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

  public :: sum_xy, distribute_xy, collect_xy
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

  public :: y2x, z2x

  public :: report_clean_output
  public :: mpiwait
  
! Radiation ray routines
  public :: radboundary_xy_recv, radboundary_xy_send
  public :: radboundary_zx_recv, radboundary_zx_send
  public :: radboundary_yz_sendrecv, radboundary_zx_sendrecv
  public :: radboundary_yz_periodic_ray, radboundary_zx_periodic_ray

! Variables
  public :: ipx, ipy, ipz, lroot, iproc, mpi_precision, nprocs
  public :: lfirst_proc_x, lfirst_proc_y, lfirst_proc_z, lfirst_proc_xy, lfirst_proc_yz, lfirst_proc_xz, lfirst_proc_xyz
  public :: llast_proc_x, llast_proc_y, llast_proc_z, llast_proc_xy, llast_proc_yz, llast_proc_xz, llast_proc_xyz
  public :: MPI_COMM_WORLD, MPI_COMM_GRID, MPI_COMM_PENCIL, MPI_COMM_XYPLANE, MPI_COMM_XZPLANE, MPI_COMM_YZPLANE, &
            MPI_COMM_XBEAM,MPI_COMM_YBEAM,MPI_COMM_ZBEAM, &
            MPI_INFO_NULL, MPI_ANY_TAG, lyang
  integer, parameter, public :: IXBEAM=1, IYBEAM=2, IZBEAM=3, IXYPLANE=12, IXZPLANE=13, IYZPLANE=23
!
! for protecting MPI_COMM_WORLD to be redefined by preprocessor
! 
  integer, parameter, public :: MPI_COMM_UNIVERSE=MPI_COMM_WORLD

  character(LEN=4), public :: cyinyang=' '
