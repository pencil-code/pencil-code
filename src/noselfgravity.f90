! $Id: noselfgravity.f90,v 1.5 2006-06-30 12:44:58 joishi Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Selfgravity

  use Cdata
  use Cparam
  use Messages

  implicit none

  include 'selfgravity.h'

  real :: rhs_poisson_const=0.

  integer :: idiag_gpoten=0

  contains

!***********************************************************************
    subroutine register_selfgravity()
!
!  Register self gravity variables.
!
!  15-may-06/anders+jeff: dummy
!
      use Cdata
!
      lselfgravity=.false.
!
    endsubroutine register_selfgravity
!***********************************************************************
    subroutine initialize_selfgravity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  15-may-06/anders+jeff: dummy
!
    endsubroutine initialize_selfgravity
!***********************************************************************
    subroutine pencil_criteria_selfgravity()
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
      if (NO_WARN) print*, lpencil_in !(keep compiler quiet)
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
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!      
      if (NO_WARN) print*, f, p !(keep compiler quiet)
!
    endsubroutine calc_pencils_selfgravity
!***********************************************************************
    subroutine calc_selfpotential(f)
!
!  Calculate the potential of self gravity.
!
!  15-may-06/anders+jeff: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if (NO_WARN) print*, f
!
    endsubroutine calc_selfpotential
!***********************************************************************
    subroutine duu_dt_selfgrav(f,df,p)
!
!  Add self gravity acceleration on gas.
!
!  15-may-06/anders+jeff: dummy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      if (NO_WARN) print*, f, df, p !(keep compiler quiet)
!        
    endsubroutine duu_dt_selfgrav
!***********************************************************************
    subroutine read_selfgravity_init_pars(unit,iostat)
!
!  Read self gravity init parameters
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_selfgravity_init_pars
!***********************************************************************
    subroutine write_selfgravity_init_pars(unit)
!
!  Write self gravity init parameters
!
      integer, intent(in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_selfgravity_init_pars
!***********************************************************************
    subroutine read_selfgravity_run_pars(unit,iostat)
!
!  Read self gravity run parameters
!    
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (NO_WARN) print*, unit, iostat
!
    endsubroutine read_selfgravity_run_pars
!***********************************************************************
    subroutine write_selfgravity_run_pars(unit)
!
!  Write self gravity run parameters
!
      integer, intent(in) :: unit
!
      if (NO_WARN) print*, unit
!
    endsubroutine write_selfgravity_run_pars
!***********************************************************************
    subroutine rprint_selfgravity(lreset,lwrite)
!
!  reads and registers print parameters relevant for gravity advance
!  dummy routine
!
!  16-may-06/anders+jeff: adapted
!
      logical :: lreset, lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
      if (NO_WARN) print*, lreset !(to keep compiler quiet)
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_gpoten=',idiag_gpoten
        write(3,*) 'ipotself=0'
      endif
!
    endsubroutine rprint_selfgravity
!***********************************************************************

endmodule Selfgravity
