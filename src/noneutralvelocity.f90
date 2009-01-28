! $Id$
!  This module takes care of everything related to neutral velocity
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED uun(3); divun; snij(3,3)
!
!***************************************************************

module NeutralVelocity

  use Cparam
  use Messages

  implicit none

  include 'neutralvelocity.h'

  !namelist /neutralvelocity_init_pars/ dummy
  !namelist /neutralvelocity_run_pars/  dummy

  contains

!***********************************************************************
    subroutine register_neutralvelocity()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  18-mar-03/axel: dummy routine
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_neutralvelocity called twice')
      first = .false.
!
      lneutralvelocity = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
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
    subroutine init_uun(f,xx,yy,zz)
!
!  initialise uun; called from start.f90
!
!  18-mar-03/axel: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if (NO_WARN) print*,f,xx,yy,zz  !(keep compiler quiet)
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
      if (NO_WARN) print*, lpencil_in  !(keep compiler quiet)
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
      if (NO_WARN) print*, f, p
!
    endsubroutine calc_pencils_neutralvelocity
!***********************************************************************
    subroutine duun_dt(f,df,p)
!
!  velocity evolution
!  calculate dud/dt = - ud.gradud - 2Omega x ud + grav + Fvisc
!  no pressure gradient force for neutral!
!
!  18-mar-03/axel: dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      if (NO_WARN) print*,f,df,p !(keep compiler quiet)
!
    endsubroutine duun_dt
!***********************************************************************
    subroutine read_neutralvelocity_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_neutralvelocity_init_pars
!***********************************************************************
    subroutine write_neutralvelocity_init_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit

    endsubroutine write_neutralvelocity_init_pars
!***********************************************************************
    subroutine read_neutralvelocity_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_neutralvelocity_run_pars
!***********************************************************************
    subroutine write_neutralvelocity_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit
    endsubroutine write_neutralvelocity_run_pars
!***********************************************************************
    subroutine rprint_neutralvelocity(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
      use Sub
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  write column where which neutral velocity variable is stored
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'iuun=',iuun
        write(3,*) 'iunx=',iunx
        write(3,*) 'iuny=',iuny
        write(3,*) 'iunz=',iunz
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_neutralvelocity
!***********************************************************************

endmodule NeutralVelocity
