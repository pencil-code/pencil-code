! $Id$

!  This module contains routines both for delta-correlated
!  and continuous forcing. The fcont pencil is only provided
!  for continuous forcing.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED fcont(3)
!
!***************************************************************

module Forcing

! Dummy module for hydro/mhd without forcing

  use Cdata
!  use General
  use Sub, only: keep_compiler_quiet
  use Messages

  implicit none

  include 'forcing.h'

  !namelist /forcing_init_pars/ dummy
  !namelist /forcing_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_rufm=0

  contains

!***********************************************************************
    subroutine register_forcing()
!
!  add forcing in timestep()
!  11-may-2002/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      lforcing = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id$")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine initialize_forcing(lstarting)
!
!  initialize random number generator in processor-dependent fashion
!  see comments in start.f90 for details
!
      use Cdata
!
      logical :: lstarting
!
      if (NO_WARN) print*,'lstarting=',lstarting !(to keep compiler quiet)
    endsubroutine initialize_forcing
!***********************************************************************
    subroutine addforce(f)
!
      use Cdata
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
      use Cdata
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      if (NO_WARN) print*,df,p
!
    endsubroutine forcing_continuous
!***********************************************************************
    subroutine read_forcing_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_forcing_init_pars
!***********************************************************************
    subroutine write_forcing_init_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit

    endsubroutine write_forcing_init_pars
!***********************************************************************
    subroutine read_forcing_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_forcing_run_pars
!***********************************************************************
    subroutine write_forcing_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit
    endsubroutine write_forcing_run_pars
!***********************************************************************
    subroutine input_persistent_forcing(id,lun,done)
!
!  Read in the stored time of the next SNI
!
      integer :: id,lun
      logical :: done
!
      if (NO_WARN) print*,id,lun,done
!
    endsubroutine input_persistent_forcing
!***********************************************************************
    subroutine output_persistent_forcing(lun)
!
!  Writes out the time of the next SNI
!
      integer :: lun
!
      if (NO_WARN) print*,lun
!
    endsubroutine output_persistent_forcing
!***********************************************************************
    subroutine rprint_forcing(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!  26-jan-04/axel: coded
!
      use Cdata
      use Sub
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
!  write column where which forcing variable is stored
!
      if (lwr) then
        write(3,*) 'i_rufm=',idiag_rufm
      endif
!
    endsubroutine rprint_forcing
!***********************************************************************

endmodule Forcing
