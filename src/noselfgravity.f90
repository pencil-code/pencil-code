! $Id$
!
!  This module takes care of self gravity.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lselfgravity = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Selfgravity
!
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include 'selfgravity.h'
!
  real :: rhs_poisson_const=0.
  real :: gravitational_const=real(G_Newton_cgs)
  namelist /selfgrav_init_pars/ &
      gravitational_const
!
  contains
!***********************************************************************
    subroutine register_selfgravity
!
!  Register self gravity variables.
!
!  15-may-06/anders+jeff: dummy
!
      use SharedVariables, only: put_shared_variable

      call put_shared_variable('rhs_poisson_const',rhs_poisson_const, caller='register_selfgravity')
      call put_shared_variable('gravitational_const',gravitational_const)

    endsubroutine register_selfgravity
!***********************************************************************
    subroutine initialize_selfgravity(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  15-may-06/anders+jeff: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      pot_spec=.false.
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_selfgravity
!***********************************************************************
    subroutine pencil_criteria_selfgravity
!
!  All pencils that the Selfgravity module depends on are specified here.
!
!  15-may-06/anders+jeff: dummy
!
    endsubroutine pencil_criteria_selfgravity
!***********************************************************************
    subroutine pencil_interdep_selfgravity(lpencil_in)
!
!  Interdependency among pencils from the Selfgravity module is specified here.
!
!  15-may-06/anders+jeff: dummy
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_selfgravity
!***********************************************************************
    subroutine calc_pencils_selfgravity(f,p)
!
!  Calculate Selfgravity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  15-may-06/anders+jeff: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_selfgravity
!***********************************************************************
    subroutine calc_selfpotential(f)
!
!  Calculate the potential of self gravity.
!
!  15-may-06/anders+jeff: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_selfpotential
!***********************************************************************
    subroutine addselfgrav(df,p)
!
!  Add self gravity acceleration on gas.
!
!  15-may-06/anders+jeff: dummy
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine addselfgrav
!***********************************************************************
    subroutine calc_diagnostics_selfgrav(p)
!
      type (pencil_case) :: p

      call keep_compiler_quiet(p)

    endsubroutine calc_diagnostics_selfgrav
!***********************************************************************
    subroutine read_selfgravity_init_pars(iomsg)
!
      use File_io, only: parallel_unit
!
      character(LEN=iomsglen), intent(out) :: iomsg

      integer :: iostat
!
      read(parallel_unit, NML=selfgrav_init_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""
      !TP: Do something ugly since init pars for no module should always be optional
      iomsg=""
!
    endsubroutine read_selfgravity_init_pars
!***********************************************************************
    subroutine write_selfgravity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_selfgravity_init_pars
!***********************************************************************
    subroutine read_selfgravity_run_pars(iomsg)
!
      character(LEN=iomsglen), intent(out) :: iomsg
!
      iomsg=""
!
    endsubroutine read_selfgravity_run_pars
!***********************************************************************
    subroutine write_selfgravity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_selfgravity_run_pars
!***********************************************************************
    subroutine rprint_selfgravity(lreset,lwrite)
!
!  Reads and registers print parameters relevant for gravity advance.
!
!  16-may-06/anders+jeff: dummy
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_selfgravity
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call keep_compiler_quiet(p_par)

    endsubroutine pushpars2c
!***********************************************************************
endmodule Selfgravity
