! $Id$
!
!  This module provide a way for users to specify custom initial
!  conditions.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
  real, dimension(-nghost:nghost), parameter :: weight = (/1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0/) / 64.0
  real, dimension(3) :: uu_left, uu_right
  real, dimension(3) :: bb_left, bb_right
  real :: rho_left, rho_right
  real :: pp_left, pp_right
!
  namelist /initial_condition_pars/ uu_left, uu_right, bb_left, bb_right, rho_left, rho_right, pp_left, pp_right
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  20-nov-13/ccyang: coded
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  20-nov-13/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f)
!
!  Initializes all the f arrays in one call.  This subroutine is called last.
!
!  21-dec-10/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_all
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  20-nov-13/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: i
!
      comp: do i = 1, 3
        call shocktube(f, iuu+i-1, uu_left(i), uu_right(i))
      enddo comp
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  20-nov-13/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call shocktube(f, ilnrho, log(rho_left), log(rho_right))
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  20-nov-13/ccyang: coded
!
      use EquationOfState, only: gamma_m1
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      energy: if (lthermal_energy) then
        call shocktube(f, ieth, pp_left / gamma_m1, pp_right / gamma_m1)
      else energy
        call fatal_error('initial_condition_ss', 'not yet implemented for this energy module. ')
      endif energy
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  20-nov-13/ccyang: coded
!
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension(:), pointer :: B_ext
      real, dimension(mx) :: penc
      integer :: ierr, i
!
      bfield: if (lbfield) then
        comp: do i = 1, 3
          call shocktube(f, ibb+i-1, bb_left(i), bb_right(i))
        enddo comp
      else bfield
        if (bb_left(1) /= bb_right(1)) call fatal_error('initial_condition_aa', 'discontinuity in bx is not allowed. ')
        call get_shared_variable('B_ext', B_ext, ierr)
        if (ierr /= 0) call fatal_error('initial_condition_aa', 'unable to get B_ext. ')
        B_ext(1) = bb_left(1)
        B_ext(2:3) = 0.0
        f(l1:l2,m1:m2,n1:n2,iax) = 0.0
        where (x < 0.0)
          penc = x * bb_left(3)
        elsewhere
          penc = x * bb_right(3)
        endwhere
        f(l1:l2,m1:m2,n1:n2,iay) = spread(spread(smooth(penc),2,ny),3,nz)
        where (x < 0.0)
          penc = -x * bb_left(2)
        elsewhere
          penc = -x * bb_right(2)
        endwhere
        f(l1:l2,m1:m2,n1:n2,iaz) = spread(spread(smooth(penc),2,ny),3,nz)
      endif bfield
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_pars(unit, iostat)
!
!  Read the namelist initial_condition_pars.
!
!  20-nov-13/ccyang: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit, NML=initial_condition_pars, IOSTAT=iostat)
      else
        read(unit, NML=initial_condition_pars)
      endif
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
!  Write the namelist initial_condition_pars.
!
!  20-nov-13/ccyang: coded
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
    subroutine shocktube(f, ivar, left, right)
!
!  Set up left and right states for variable ivar.
!
!  20-nov-13/ccyang: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: ivar
      real, intent(in) :: left, right
!
      real, dimension(mx) :: penc
!
      where (x < 0.0)
        penc = left
      elsewhere
        penc = right
      endwhere
      f(l1:l2,m1:m2,n1:n2,ivar) = spread(spread(smooth(penc),2,ny),3,nz)
!
    endsubroutine shocktube
!***********************************************************************
    function smooth(penc) result(smoothed)
!
!  Return smoothed pencil.
!
!  20-nov-13/ccyang: coded
!
      real, dimension(mx), intent(in) :: penc
      real, dimension(nx) :: smoothed
!
      integer :: i
!
      smoothed = 0.0
      do i = -nghost, nghost
        smoothed = smoothed + weight(i) * penc(l1+i:l2+i)
      enddo
!
    endfunction smooth 
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
