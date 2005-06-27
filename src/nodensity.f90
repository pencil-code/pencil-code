! $Id: nodensity.f90,v 1.34 2005-06-27 19:47:28 dobler Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldensity = .false.
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED rho,lnrho,rho1,glnrho,uglnrho,sglnrho,del2lnrho,hlnrho
!
!***************************************************************

module Density

  use Cparam
  use EquationOfState, only: cs0,cs20,lnrho0,rho0,lcalc_cp, &
                             gamma,gamma1,cs2top,cs2bot

  implicit none

  include 'density.h'


  !namelist /density_init_pars/ dummy
  !namelist /density_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_rhom=0
  integer :: idiag_rhomin=0,idiag_rhomax=0

  real :: beta_dlnrhodr=0.0, beta_dlnrhodr_scaled=0.0

  contains

!***********************************************************************
    subroutine register_density()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrho; increase nvar accordingly.
!
!   8-jun-02/axel: adapted from density
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_density called twice')
      first = .false.
!
!ajwm      ldensity = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: nodensity.f90,v 1.34 2005-06-27 19:47:28 dobler Exp $")
!
!ajwm Necessary? added incase
!      gamma=1.
!      gamma1=0.
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded 
!
  use Cdata, only: lentropy
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
!
!  in the force-free model, there is no pressure gradient
!  (but allow for the possibility of pressure gradients from entropy)
!
!      if (.not.lentropy) then
!        cs0=0.
!        cs20=0.
!      endif
!
      if (ip == 0) print*,f,lstarting ! keep compiler quiet
    endsubroutine initialize_density
!***********************************************************************
    subroutine init_lnrho(f,xx,yy,zz)
!
!  initialise lnrho; called from start.f90
!
!   7-jun-02/axel: adapted from density
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if(NO_WARN) print*,f,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_lnrho
!***********************************************************************
    subroutine pencil_criteria_density()
! 
!  All pencils that the Density module depends on are specified here.
! 
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_density
!***********************************************************************
    subroutine pencil_interdep_density(lpencil_in)
!
!  Interdependency among pencils from the Density module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in !(keep compiler quiet)
!
    endsubroutine pencil_interdep_density
!***********************************************************************
    subroutine calc_pencils_density(f,p)
!
!  Calculate Density pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f       
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
! rho
      if (lpencil(i_rho)) p%rho=rho0
! lnrho
      if (lpencil(i_lnrho)) p%lnrho=lnrho0
! rho1
      if (lpencil(i_rho1)) p%rho1=1./rho0
! glnrho
      if (lpencil(i_glnrho)) p%glnrho=0.
! uglnrho
      if (lpencil(i_uglnrho)) p%uglnrho=0.
! hlnrho
      if (lpencil(i_hlnrho)) p%hlnrho=0.
!     
      if(NO_WARN) print*,f(1,1,1,1)   !(keep compiler quiet)
!
    endsubroutine calc_pencils_density
!***********************************************************************
    subroutine dlnrho_dt(f,df,p)
!
!  continuity equation, dummy routine
!
!   7-jun-02/axel: adapted from density
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,df,p
!
      if(NO_WARN) print*,f,df,p
!
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine read_density_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
                                                                                                   
      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit
                                                                                                   
    endsubroutine read_density_init_pars
!***********************************************************************
    subroutine write_density_init_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      if (NO_WARN) print*,unit
                                                                                                   
    endsubroutine write_density_init_pars
!***********************************************************************
    subroutine read_density_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
                                                                                                   
      if (present(iostat) .and. (NO_WARN)) print*,iostat
      if (NO_WARN) print*,unit
                                                                                                   
    endsubroutine read_density_run_pars
!***********************************************************************
    subroutine write_density_run_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      if (NO_WARN) print*,unit
    endsubroutine write_density_run_pars
!***********************************************************************
    subroutine rprint_density(lreset,lwrite)
!
!  reads and registers print parameters relevant for compressible part
!
!   8-jun-02/axel: adapted from density
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
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_rhom=0
        idiag_rhomin=0; idiag_rhomax=0
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_rhom=',idiag_rhom
        write(3,*) 'i_rhomin=',idiag_rhomin
        write(3,*) 'i_rhomax=',idiag_rhomax
        write(3,*) 'nname=',nname
        write(3,*) 'ilnrho=',ilnrho
      endif
!
    endsubroutine rprint_density
!***********************************************************************
    subroutine bc_lnrho_temp_z(f,topbot)
!
!  dummy routine
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if(NO_WARN) print*,f,topbot
    endsubroutine bc_lnrho_temp_z
!***********************************************************************
    subroutine bc_lnrho_pressure_z(f,topbot)
!
!  dummy routine
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if(NO_WARN) print*,f,topbot
    endsubroutine bc_lnrho_pressure_z
!***********************************************************************

endmodule Density
