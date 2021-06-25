! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! PENCILS PROVIDED x_mn; y_mn; z_mn; r_mn; r_mn1
! PENCILS PROVIDED phix; phiy
! PENCILS PROVIDED pomx; pomy
! PENCILS PROVIDED rcyl_mn; rcyl_mn1; phi_mn
! PENCILS PROVIDED evr(3); rr(3); evth(3)
!
!***************************************************************
module Grid
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  private
!
  public :: construct_grid
  public :: pencil_criteria_grid
  public :: pencil_interdep_grid
  public :: calc_pencils_grid
  public :: initialize_grid
  public :: get_grid_mn
  public :: box_vol
  public :: save_grid
  public :: coords_aux
  public :: real_to_index, inverse_grid
  public :: grid_bound_data
  public :: set_coorsys_dimmask
  public :: construct_serial_arrays
  public :: coarsegrid_interp
!
! public :: find_star
  public :: grid_profile

  interface grid_profile
    module procedure grid_profile_0D
    module procedure grid_profile_1D
  endinterface
!
  interface calc_pencils_grid
    module procedure calc_pencils_grid_pencpar
    module procedure calc_pencils_grid_std
  endinterface calc_pencils_grid
!
  integer, parameter :: BOT=1, TOP=2
!
  contains
!***********************************************************************
    subroutine construct_grid(x,y,z,dx,dy,dz)
!
!  Dummy routine.
!
      real, dimension(mx) :: x
      real, dimension(my) :: y
      real, dimension(mz) :: z
      real :: dx,dy,dz
!
    endsubroutine construct_grid
!***********************************************************************
    subroutine set_coorsys_dimmask
!
!  Dummy routine.
!
    endsubroutine set_coorsys_dimmask
!***********************************************************************
    subroutine initialize_grid
!
!  Dummy routine.
!
    endsubroutine initialize_grid
!***********************************************************************
    subroutine coarsegrid_interp(f,ivar1,ivar2)
!
!  Dummy routine.
!
      real, dimension(:,:,:,:) :: f
      integer, optional :: ivar1,ivar2
!
    endsubroutine coarsegrid_interp
!***********************************************************************
    subroutine save_grid(lrestore)
!
!  Dummy routine.
!
      logical, optional :: lrestore
!
    endsubroutine save_grid
!***********************************************************************
    subroutine coords_aux(x,y,z)
!
!  Dummy routine.
!
      real, dimension(:) :: x,y,z
!
    endsubroutine coords_aux
!***********************************************************************
    subroutine box_vol
!
!  Dummy routine.
!
    endsubroutine box_vol
!***********************************************************************
    subroutine pencil_criteria_grid
!
!  Dummy routine.
!
    endsubroutine pencil_criteria_grid
!***********************************************************************
    subroutine pencil_interdep_grid(lpencil_in)
!
!  Dummy routine.
!
      logical, dimension(npencils) :: lpencil_in
!
    endsubroutine pencil_interdep_grid
!***********************************************************************
    subroutine calc_pencils_grid_std(f,p)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
    endsubroutine calc_pencils_grid_std
!***********************************************************************
    subroutine calc_pencils_grid_pencpar(f,p,lpenc_loc)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_grid_pencpar
!***********************************************************************
    subroutine grid_profile_0D(xi,grid_func,g,gder1,gder2,param,dxyz,xistep,delta)
!
!  Dummy routine.
!
      real              :: xi
      character(len=*)  :: grid_func
      real              :: g
      real, optional    :: gder1, gder2
      real, optional    :: param
      real, optional, dimension(3) :: dxyz
      real, optional, dimension(2) :: xistep, delta
!
    endsubroutine grid_profile_0D
!***********************************************************************
    subroutine grid_profile_1D(xi,grid_func,g,gder1,gder2,param,dxyz,xistep,delta)
!
!  Dummy routine.
!
      real, dimension(:)                    :: xi
      character(len=*)                      :: grid_func
      real, dimension(size(xi,1))           :: g
      real, dimension(size(xi,1)), optional :: gder1,gder2
      real, optional                        :: param
      real, optional, dimension(3) :: dxyz
      real, optional, dimension(2) :: xistep,delta
!
    endsubroutine grid_profile_1D
!***********************************************************************
!   function find_star(xi_lo,xi_up,x_lo,x_up,x_star,grid_func) result (xi_star)
!
!  Dummy routine.
!
!     real :: xi_lo,xi_up,x_lo,x_up,x_star
!     character(len=*) :: grid_func
!
      !eal :: xi_star,dxi,tol
      !eal :: g_lo,gder_lo
      !eal :: g_up,gder_up
      !eal :: f   ,fder
      !nteger, parameter :: maxit=1000
      !ogical :: lreturn
      !nteger :: it
!
    !ndfunction find_star
!***********************************************************************
    subroutine real_to_index(n, x, xi)
!
!  Dummy routine.
!
      integer :: n
      real, dimension(n,3) :: x
      real, dimension(n,3) :: xi
!
    endsubroutine real_to_index
!***********************************************************************
    subroutine inverse_grid(dir, x, xi, local)
!
!  Dummy routine.
!
      integer :: dir
      real, dimension(:) :: x
      real, dimension(:) :: xi
      logical, intent(in), optional :: local
!
    endsubroutine inverse_grid
!***********************************************************************
    subroutine construct_serial_arrays
!
!  Dummy routine.
!
    endsubroutine construct_serial_arrays
!***********************************************************************
    subroutine get_grid_mn
!
!  Dummy routine.
!
    endsubroutine get_grid_mn
!***********************************************************************
    subroutine calc_coeffs_1( grid, coeffs )
!
!  Dummy routine.
!
      real, dimension(-2:3) :: grid
      real, dimension(-3:3) :: coeffs
!
  endsubroutine calc_coeffs_1
!***********************************************************************
    subroutine grid_bound_data
!
!  Dummy routine
!
    endsubroutine grid_bound_data
!***********************************************************************
endmodule Grid
