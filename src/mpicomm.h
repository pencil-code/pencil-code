!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
  private

  public :: mpicomm_init, mpifinalize
  public :: mpibarrier
  public :: stop_it, stop_it_if_any
  public :: die_gracefully

  public :: mpirecv_real, mpirecv_int
  public :: mpisend_real, mpisend_int
  public :: mpireduce_sum, mpireduce_max, mpireduce_min
  public :: mpiallreduce_max
  public :: mpireduce_sum_int, mpireduce_sum_double
  public :: mpireduce_or
  public :: mpibcast_real,mpibcast_logical
  public :: mpibcast_double
  public :: mpibcast_int, mpibcast_char

  public :: mpiwtime, mpiwtick, fold_df, fold_f

  public :: start_serialize,end_serialize
  public :: initiate_isendrcv_bdry, finalize_isendrcv_bdry
  public :: initiate_isendrcv_shockbdry, finalize_isendrcv_shockbdry
  public :: initiate_shearing, finalize_shearing

  public :: transp
  public :: transform_fftpack,transform_fftpack_2d, transform_fftpack_1d 
  public :: transform_fftpack_shear
  public :: transform_nr, transform_i, fourier_shift_yz

! Radiation ray routines
  public :: radboundary_xy_recv,radboundary_xy_send
  public :: radboundary_zx_recv,radboundary_zx_send
  public :: radboundary_zx_sendrecv
  public :: radboundary_zx_periodic_ray

! Variables
  public :: ipx,ipy,ipz,lroot,iproc
