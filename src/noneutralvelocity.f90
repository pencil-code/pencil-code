! $Id$
!
!  This module takes care of everything related to neutral velocity
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lneutralvelocity = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED uun(3); divun; snij(3,3)
!
!***************************************************************
module NeutralVelocity
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'neutralvelocity.h'
!
  contains
!***********************************************************************
    subroutine register_neutralvelocity()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  18-mar-03/axel: dummy routine
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_neutralvelocity
!***********************************************************************
    subroutine initialize_neutralvelocity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-mar-03/axel: dummy routine
!
    endsubroutine initialize_neutralvelocity
!***********************************************************************
    subroutine init_uun(f)
!
!  Initialise uun; called from start.f90.
!
!  18-mar-03/axel: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_uun
!***********************************************************************
    subroutine pencil_criteria_neutralvelocity()
!
!  All pencils that the Neutralvelocity module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_neutralvelocity
!***********************************************************************
    subroutine pencil_interdep_neutralvelocity(lpencil_in)
!
!  Interdependency among pencils provided by the Neutralvelocity module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_neutralvelocity
!***********************************************************************
    subroutine calc_pencils_neutralvelocity(f,p)
!
!  Calculate Neutralvelocity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  13-nov-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_neutralvelocity
!***********************************************************************
    subroutine duun_dt(f,df,p)
!
!  Velocity evolution.
!  Calculate dud/dt = - ud.gradud - 2Omega x ud + grav + Fvisc
!  No pressure gradient force for neutral!
!
!  18-mar-03/axel: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine duun_dt
!***********************************************************************
    subroutine read_neutralvelocity_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_neutralvelocity_init_pars
!***********************************************************************
    subroutine write_neutralvelocity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_neutralvelocity_init_pars
!***********************************************************************
    subroutine read_neutralvelocity_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_neutralvelocity_run_pars
!***********************************************************************
    subroutine write_neutralvelocity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_neutralvelocity_run_pars
!***********************************************************************
    subroutine rprint_neutralvelocity(lreset,lwrite)
!
!  Reads and registers print parameters relevant for neutral velocity.
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Write column where which neutral velocity variable is stored.
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'iuun=',iuun
        write(3,*) 'iunx=',iunx
        write(3,*) 'iuny=',iuny
        write(3,*) 'iunz=',iunz
      endif
!
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_neutralvelocity
!***********************************************************************
endmodule NeutralVelocity
