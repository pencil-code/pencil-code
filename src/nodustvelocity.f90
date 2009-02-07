! $Id$


!  This module takes care of everything related to velocity

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED uud(3,ndustspec); divud(ndustspec); sdij(3,3,ndustspec)
!
!***************************************************************

module Dustvelocity

  use Cparam
  use Messages

  implicit none

  include 'dustvelocity.h'

  public :: dust_geometry, dimd1, rhods, surfd, mdplus, mdminus
  public :: ad, scolld, ustcst, tausd1, tausd
  public :: unit_md, dust_chemistry, mumon, mmon, mi, md
  character (len=labellen) :: dust_geometry='sphere', dust_chemistry='nothing'
  real, dimension(ndustspec,ndustspec) :: scolld
  real, dimension(nx,ndustspec) :: tausd1
  real, dimension(ndustspec) :: md=1.0, mdplus=0.0, mdminus=0.0, surfd=0.0
  real, dimension(ndustspec) :: mi=0.0, ad=1.0, tausd=0.0
  real :: dimd1=0.0, rhods=0.0, ustcst=0.0, unit_md=0.0
  real :: mumon=0.0, mmon=0.0

  !namelist /dustvelocity_init_pars/ dummy
  !namelist /dustvelocity_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_ud2m=0,idiag_udm2=0,idiag_oudm=0,idiag_od2m=0
  integer :: idiag_udrms=0,idiag_udmax=0,idiag_odrms=0,idiag_odmax=0
  integer :: idiag_udxmz=0,idiag_udymz=0,idiag_udzmz=0,idiag_udmx=0
  integer :: idiag_udmy=0,idiag_udmz=0
  integer :: idiag_udxmxy=0,idiag_udymxy=0,idiag_udzmxy=0
  integer :: idiag_divud2m=0,idiag_epsKd=0

  !! SHOULDN'T REALLY BE PUBLIC!!
  real :: nd0=1.,rhod0=1.
  logical :: ldustcoagulation=.false., ldustcondensation=.false.

  contains

!***********************************************************************
    subroutine register_dustvelocity()
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
      ldustvelocity = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
          "$Id$")
!
    endsubroutine register_dustvelocity
!***********************************************************************
    subroutine initialize_dustvelocity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-mar-03/axel: dummy routine
!
    endsubroutine initialize_dustvelocity
!***********************************************************************
    subroutine copy_bcs_dust
!
!  Copy boundary conditions after first dust species to end of array
!
!  27-feb-04/anders: dummy routine
!
    endsubroutine copy_bcs_dust
!***********************************************************************
    subroutine init_uud(f)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Dustvelocity module, if there was one.
!
!  18-mar-03/axel: dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (NO_WARN) print*,f  !(keep compiler quiet)
!
    endsubroutine init_uud
!***********************************************************************
    subroutine pencil_criteria_dustvelocity()
!
!  All pencils that the Dustvelocity module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_dustvelocity
!***********************************************************************
    subroutine pencil_interdep_dustvelocity(lpencil_in)
!
!  Interdependency among pencils provided by the Dustvelocity module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in  !(keep compiler quiet)
!
    endsubroutine pencil_interdep_dustvelocity
!***********************************************************************
    subroutine calc_pencils_dustvelocity(f,p)
!
!  Calculate Dustvelocity pencils.
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
    endsubroutine calc_pencils_dustvelocity
!***********************************************************************
    subroutine duud_dt(f,df,p)
!
!  velocity evolution
!  calculate dud/dt = - ud.gradud - 2Omega x ud + grav + Fvisc
!  no pressure gradient force for dust!
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
    endsubroutine duud_dt
!***********************************************************************
    subroutine shearingdust(f,df)
!
!  Calculates the shear terms, -uy0*df/dy (shearing sheat approximation)
!
!  6-dec-03/anders: Copied from shearing
!
      use Cparam
      use Deriv
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      if (NO_WARN) print*,f,df !(keep compiler quiet)
!
    end subroutine shearingdust
!***********************************************************************
    subroutine read_dustvelocity_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_dustvelocity_init_pars
!***********************************************************************
    subroutine write_dustvelocity_init_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit

    endsubroutine write_dustvelocity_init_pars
!***********************************************************************
    subroutine read_dustvelocity_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit

    endsubroutine read_dustvelocity_run_pars
!***********************************************************************
    subroutine write_dustvelocity_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit
    endsubroutine write_dustvelocity_run_pars
!***********************************************************************
    subroutine rprint_dustvelocity(lreset,lwrite)
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
!  write column where which dust velocity variable is stored
!
      if (lwr) then
        write(3,*) 'i_ud2m=',idiag_ud2m
        write(3,*) 'i_udm2=',idiag_udm2
        write(3,*) 'i_od2m=',idiag_od2m
        write(3,*) 'i_oudm=',idiag_oudm
        write(3,*) 'i_udrms=',idiag_udrms
        write(3,*) 'i_udmax=',idiag_udmax
        write(3,*) 'i_odrms=',idiag_odrms
        write(3,*) 'i_odmax=',idiag_odmax
        write(3,*) 'i_udmx=',idiag_udmx
        write(3,*) 'i_udmy=',idiag_udmy
        write(3,*) 'i_udmz=',idiag_udmz
        write(3,*) 'i_divud2m=',idiag_divud2m
        write(3,*) 'i_epsKd=',idiag_epsKd
        write(3,*) 'i_udxmz=',idiag_udxmz
        write(3,*) 'i_udymz=',idiag_udymz
        write(3,*) 'i_udzmz=',idiag_udzmz
        write(3,*) 'i_udxmxy=',idiag_udxmxy
        write(3,*) 'i_udymxy=',idiag_udymxy
        write(3,*) 'i_udzmxy=',idiag_udzmxy
        write(3,*) 'nname=',nname
        write(3,*) 'iuud=',iuud
        write(3,*) 'iudx=',iudx
        write(3,*) 'iudy=',iudy
        write(3,*) 'iudz=',iudz
      endif
!
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_dustvelocity
!***********************************************************************

endmodule Dustvelocity
