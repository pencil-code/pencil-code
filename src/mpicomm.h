!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!$Id$
!
  private

  public :: mpicomm_init, initialize_mpicomm, mpifinalize
  public :: mpibarrier
  public :: stop_it, stop_it_if_any
  public :: die_gracefully, die_immediately
  public :: check_emergency_brake

  public :: mpirecv_logical, mpirecv_real, mpirecv_int
  public :: mpisend_logical, mpisend_real, mpisend_int
  public :: mpireduce_sum_int, mpireduce_sum, mpireduce_sum_double
  public :: mpireduce_max, mpireduce_min
  public :: mpiallreduce_max, mpiallreduce_sum, mpiallreduce_sum_int
  public :: mpiallreduce_sum_arr, mpiallreduce_sum_arr2
  public :: mpireduce_or, mpireduce_and
  public :: mpibcast_real, mpibcast_logical
  public :: mpibcast_real_arr
  public :: mpibcast_double
  public :: mpibcast_int, mpibcast_char, mpireduce_max_scl_int
  public :: mpigather_xy, mpigather_and_out, mpimerge_1d, mpigather_z

  public :: mpiwtime, mpiwtick

  public :: start_serialize,end_serialize
  public :: initiate_isendrcv_bdry, finalize_isendrcv_bdry
  public :: initiate_shearing, finalize_shearing

  public :: transp, transp_xy, transp_xy_other, transp_other
  public :: transp_xz, transp_zx

  public :: communicate_bc_aa_pot, fill_zghostzones_3vec

  public :: distribute_to_pencil_xy, collect_from_pencil_xy
  public :: remap_to_pencil_x, unmap_from_pencil_x
  public :: remap_to_pencil_y, unmap_from_pencil_y
  public :: remap_to_pencil_z, unmap_from_pencil_z
  public :: remap_to_pencil_xy, unmap_from_pencil_xy, transp_pencil_xy
  public :: remap_to_pencil_yz, unmap_from_pencil_yz

  public :: z2x
  public :: MPI_adi_x, MPI_adi_z

  public :: parallel_open, parallel_close
  public :: parallel_file_exists
  public :: parallel_count_lines
  
! Radiation ray routines
  public :: radboundary_xy_recv, radboundary_xy_send
  public :: radboundary_zx_recv, radboundary_zx_send
  public :: radboundary_zx_sendrecv
  public :: radboundary_zx_periodic_ray

! Variables
  public :: ipx, ipy, ipz, lroot, iproc
  public :: lfirst_proc_x, lfirst_proc_y, lfirst_proc_z, lfirst_proc_xy, lfirst_proc_yz, lfirst_proc_xz, lfirst_proc_xyz
  public :: llast_proc_x, llast_proc_y, llast_proc_z, llast_proc_xy, llast_proc_yz, llast_proc_xz, llast_proc_xyz

