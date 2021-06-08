  public :: initialize_deriv, finalize_deriv
  public :: der, der2, der3, der4, der5, der6, derij, der5i1j, der5_single
  public :: der4i2j,der2i2j2k,der3i3j,der3i2j1k,der4i1j1k
  public :: der_pencil, der2_pencil, der6_pencil
  public :: deri_3d_inds
  public :: der_x,der2_x
  public :: der_z,der2_z
  public :: der_upwind1st
  public :: der_onesided_4_slice
  public :: der_onesided_4_slice_other
  public :: der2_minmod
  public :: heatflux_deriv_x
  public :: set_ghosts_for_onesided_ders
  public :: bval_from_neumann, bval_from_3rd, bval_from_4th
!
!debug  integer, parameter :: icount_der   = 1         !DERCOUNT
!debug  integer, parameter :: icount_der2  = 2         !DERCOUNT
!debug  integer, parameter :: icount_der4  = 3         !DERCOUNT
!debug  integer, parameter :: icount_der5  = 4         !DERCOUNT
!debug  integer, parameter :: icount_der6  = 5         !DERCOUNT
!debug  integer, parameter :: icount_derij = 6         !DERCOUNT
!debug  integer, parameter :: icount_der_upwind1st = 7 !DERCOUNT
!debug  integer, parameter :: icount_der_other = 8     !DERCOUNT
!
  interface der                 ! Overload the der function
    module procedure der_main   ! derivative of an 'mvar' variable
    module procedure der_other  ! derivative of another field
    module procedure der_pencil ! derivative of a pencil
  endinterface
!
  interface der2                 ! Overload the der2 function
    module procedure der2_main   ! derivative of an 'mvar' variable
    module procedure der2_other  ! derivative of another field
    module procedure der2_pencil ! derivative of a penci
  endinterface
!
  interface der6                 ! Overload the der6 function
    module procedure der6_main   ! derivative of an 'mvar' variable
    module procedure der6_other  ! derivative of another field
    module procedure der6_pencil ! derivative of a pencil
  endinterface
!
  interface derij                 ! Overload the der function
    module procedure derij_main   ! derivative of an 'mvar' variable
    module procedure derij_other  ! derivative of another field
  endinterface
!
  interface  der_onesided_4_slice                ! Overload the der function
    module procedure  der_onesided_4_slice_main  ! derivative of an 'mvar' variable
    module procedure  der_onesided_4_slice_main_pt
    module procedure  der_onesided_4_slice_other ! derivative of another field
    module procedure  der_onesided_4_slice_other_pt
  endinterface
!
  interface bval_from_neumann
    module procedure bval_from_neumann_scl
    module procedure bval_from_neumann_arr
  endinterface
!
  interface bval_from_3rd
    module procedure bval_from_3rd_scl
    module procedure bval_from_3rd_arr
  endinterface
!
  interface bval_from_4th
    module procedure bval_from_4th_scl
    module procedure bval_from_4th_arr
  endinterface
!

