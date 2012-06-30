! $Id$
!
!  This module contains routines both for delta-correlated
!  and continuous forcing. The fcont pencil is only provided
!  for continuous forcing.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lforcing = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED fcont(3)
!
!***************************************************************
module Forcing
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  logical :: lhydro_forcing=.false.,ltestflow_forcing=.false.
  include 'forcing.h'
!
  integer :: idiag_rufm=0
!
  contains
!***********************************************************************
    subroutine register_forcing()
!
!  add forcing in timestep()
!  11-may-2002/wolf: coded
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine initialize_forcing(lstarting)
!
!  initialize random number generator in processor-dependent fashion
!  see comments in start.f90 for details
!
      logical :: lstarting
!
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_forcing
!***********************************************************************
    subroutine addforce(f)
!
!  add forcing in timestep()
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine addforce
!***********************************************************************
    subroutine calc_lforcing_cont_pars(f)
!
!  precalculate parameters that are new at each timestep,
!  but the same for all pencils
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lforcing_cont_pars
!***********************************************************************
    subroutine pencil_criteria_forcing()
!
!  Dummy routine
!
    endsubroutine pencil_criteria_forcing
!***********************************************************************
    subroutine pencil_interdep_forcing(lpencil_in)
!
!  Dummy routine
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_forcing
!***********************************************************************
    subroutine calc_pencils_forcing(f,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (lpencil(i_fcont)) p%fcont=0.
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_forcing
!***********************************************************************
    subroutine forcing_continuous(df,p)
!
!  dummy routine
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine forcing_continuous
!***********************************************************************
    subroutine forcing_cont(force)
!
      real, dimension (nx,3), intent(out) :: force
!
      call keep_compiler_quiet(force)
!
    endsubroutine forcing_cont
!***********************************************************************
    subroutine read_forcing_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_forcing_init_pars
!***********************************************************************
    subroutine write_forcing_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_forcing_init_pars
!***********************************************************************
    subroutine read_forcing_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_forcing_run_pars
!***********************************************************************
    subroutine write_forcing_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_forcing_run_pars
!***********************************************************************
    subroutine input_persistent_forcing(id,done)
!
!  Read in the stored time of the next SNI
!
      integer :: id
      logical :: done
!
      call keep_compiler_quiet(id)
      call keep_compiler_quiet(done)
!
    endsubroutine input_persistent_forcing
!***********************************************************************
    logical function output_persistent_forcing()
!
!  Writes out the time of the next SNI
!
!  16-nov-11/MR: changed into logical function
!
      output_persistent_forcing = .false.

    endfunction output_persistent_forcing
!***********************************************************************
    subroutine rprint_forcing(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!  26-jan-04/axel: coded
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_rufm=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_noforcing: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rufm',idiag_rufm)
      enddo
!
    endsubroutine rprint_forcing
!***********************************************************************
    subroutine forcing_clean_up
!
!   12-aug-09/dhruba: coded
!   dummy routine.
!
    endsubroutine forcing_clean_up
!***********************************************************************
endmodule Forcing
