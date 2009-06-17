! $Id: nopolymer.f90 $
!  This modules deals with all aspects of polymers.
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpolymer = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!
!***************************************************************

module Polymer

  use Cdata
  use Cparam
  use Messages

  implicit none

  include 'record_types.h'
  include 'polymer.h'


  contains

!***********************************************************************
    subroutine register_polymer()
!
!  Initialise variables which should know that we solve for the vector
!  potential: ipp, etc; increase nvar accordingly
!
!  14-Aug-08 : Dhruba 
!
      use Cdata
      use Mpicomm
      use Sub

!
    endsubroutine register_polymer
!***********************************************************************
    subroutine initialize_polymer(f,lstarting)
!
!  Perform any post-parameter-read initialization
!
!  14-aug-08/dhruba: initialize polymer field (dummy at present)
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical :: lstarting
!
      if (NO_WARN) print*, lstarting
!
    endsubroutine initialize_polymer
!***********************************************************************
    subroutine init_pp(f)
!
!  initialise polymer field; called from start.f90
!   14-aug-2008/dhruba: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f

!
    endsubroutine init_pp
!***********************************************************************
    subroutine pencil_criteria_polymer()
!
!   All pencils that the Polymer module depends on are specified here.
!
!  19-11-04/anders: coded
!
      use Cdata

    endsubroutine pencil_criteria_polymer
!***********************************************************************
    subroutine pencil_interdep_polymer(lpencil_in)
!
!  Interdependency among pencils from the Polymer module is specified here.
!
!  18-aug-2008/dhruba: coded
      logical, dimension(npencils) :: lpencil_in
!
    endsubroutine pencil_interdep_polymer
!***********************************************************************
    subroutine calc_pencils_polymer(f,p)
!
!  Calculate Magnetic pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  19-nov-04/anders: coded
!
      use Cdata
      use Sub
      use Deriv

!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
       
    endsubroutine calc_pencils_polymer
!***********************************************************************
    subroutine dpp_dt(f,df,p)
!
!  polymer evolution
!
      real, dimension (mx,my,mz,mfarray) :: f,df
      type (pencil_case) :: p
!
    endsubroutine dpp_dt
!***********************************************************************
    subroutine read_polymer_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!

    endsubroutine read_polymer_init_pars
!***********************************************************************
    subroutine write_polymer_init_pars(unit)
      integer, intent(in) :: unit
!

!
    endsubroutine write_polymer_init_pars
!***********************************************************************
    subroutine read_polymer_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!

    endsubroutine read_polymer_run_pars
!***********************************************************************
    subroutine write_polymer_run_pars(unit)
     integer, intent(in) :: unit
    endsubroutine write_polymer_run_pars

!***********************************************************************
    subroutine get_slices_polymer(f,slices)
!
!  Write slices for animation of polymeric variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices

    endsubroutine get_slices_polymer
!***********************************************************************
    subroutine rprint_polymer(lreset,lwrite)
!
!  reads and registers print parameters relevant for polymer
!
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!

    endsubroutine rprint_polymer
!***********************************************************************
endmodule Polymer
