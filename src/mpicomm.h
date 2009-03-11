!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!$Id$
!
  private

  public :: mpicomm_init, mpifinalize
  public :: mpibarrier
  public :: stop_it, stop_it_if_any
  public :: die_gracefully
  public :: check_emergency_brake

  public :: mpirecv_logical, mpirecv_real, mpirecv_int
  public :: mpisend_logical, mpisend_real, mpisend_int
  public :: mpireduce_sum, mpireduce_max, mpireduce_min
  public :: mpireduce_sum_scl, mpiallreduce_max, mpiallreduce_sum
  public :: mpireduce_sum_int, mpireduce_sum_double,mpiallreduce_sum_int
  public :: mpireduce_or, mpireduce_and
  public :: mpibcast_real,mpibcast_logical
  public :: mpibcast_real_arr
  public :: mpibcast_double
  public :: mpibcast_int, mpibcast_char

  public :: mpiwtime, mpiwtick

  public :: start_serialize,end_serialize
  public :: initiate_isendrcv_bdry, finalize_isendrcv_bdry
  public :: initiate_isendrcv_shockbdry, finalize_isendrcv_shockbdry
  public :: initiate_shearing, finalize_shearing

  public :: transp,transp_xy,transp_xz,transp_zx
  public :: z2x
  public :: communicate_bc_aa_pot,transp_xy_other,transp_other

! Radiation ray routines
  public :: radboundary_xy_recv,radboundary_xy_send
  public :: radboundary_zx_recv,radboundary_zx_send
  public :: radboundary_zx_sendrecv
  public :: radboundary_zx_periodic_ray

! Variables
  public :: ipx,ipy,ipz,lroot,iproc
