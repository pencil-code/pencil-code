! $Id$
!
!  This module sets the open boundary conditions for density in
!  the vertical direction, scaled by stratificaiton.
!
!  Reference: Simon, J. B., Hawley, J. F., & Beckwith, K. 2011, ApJ, 730, 94
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
!***************************************************************
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
  real, dimension(nghost) :: top = 0.0, bot = 0.0
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  20-jul-15/ccyang: coded.
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  12-aug-15/ccyang: coded
!
      use EquationOfState, only: gamma, cs20
      use Gravity, only: potential
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      real, dimension(mz) :: rho0z
      integer :: k
!
      call keep_compiler_quiet(f)
!
!  Get the vertical density stratification.
!
      do k = 1, mz
        call potential(z=z(k), pot=rho0z(k))
      enddo
      rho0z = exp(-gamma / cs20 * rho0z)
!
!  Find the scale factors needed for the boundary conditions.
!
      bot = rho0z(1:nghost) / rho0z(n1)
      top = rho0z(n2+1:mz) / rho0z(n2)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine special_boundconds(f, bc)
!
!  Sets the boundary conditions in z direction.
!
!  20-jul-15/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      type(boundary_condition), intent(inout) :: bc
!
      integer :: k
!
!  Check the specified name and field.
!
      if (bc%bcname /= 'cps' .or. bc%ivar /= irho) return
!
!  Apply the boundary conditions.
!
      ghost: do k = 1, nghost
        if (bc%location < 0) then
          f(:,:,k,irho) = bot(k) * f(:,:,n1,irho)
        else
          f(:,:,n2+k,irho) = top(k) * f(:,:,n2,irho)
        endif
      enddo ghost
      bc%done = .true.
!
    endsubroutine special_boundconds
!***********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
