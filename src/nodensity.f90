! $Id: nodensity.f90,v 1.28 2004-10-27 14:52:56 ajohan Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Density

  use Cparam
  use Ionization, only: cs0,cs20,lnrho0,rho0,lcalc_cp,gamma,gamma1,cs2top,cs2bot

  implicit none

  real :: cs2cool=0.
  real :: b_ell=1., rbound=1.
  character (len=labellen) :: initlnrho='nothing', initlnrho2='nothing'

  integer :: dummy           ! We cannot define empty namelists
  namelist /density_init_pars/ dummy
  namelist /density_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: i_ekin=0,i_rhom=0,i_ekintot=0,i_rhomin=0,i_rhomax=0

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
      ldensity = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: nodensity.f90,v 1.28 2004-10-27 14:52:56 ajohan Exp $")
!
!ajwm Necessary? added incase
      gamma=1.
      gamma1=0.
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
      if (.not.lentropy) then
        cs0=0.
        cs20=0.
      endif
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
      if(ip==0) print*,f,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_lnrho
!***********************************************************************
    subroutine calculate_vars_rho(f,lnrho,rho,rho1)
!
!   Calculation of rho1
!   dummy routine
!
!   09-febr-04/bing: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f       
      real, dimension (nx) :: lnrho,rho,rho1

      intent(in) :: f
      intent(out) :: lnrho,rho,rho1
!
!  set rho1=1 (useful for so-called force-free model)
!
      rho1=1.
      rho=1.
      lnrho=0.
!     
      if(ip==0) print*,f(1,1,1,1)   !(keep compiler quiet)
    endsubroutine calculate_vars_rho
!***********************************************************************
    subroutine dlnrho_dt(f,df,uu,divu,lnrho,rho,glnrho,shock,gshock)
!
!  continuity equation, dummy routine
!
!   7-jun-02/axel: adapted from density
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu,glnrho,gshock
      real, dimension (nx) :: lnrho,divu,shock
!
!  will be accessed in noentropy
!
      glnrho=0.
!
      if(ip==0) print*,f,df,uu,divu,lnrho,rho,shock,gshock
!
    endsubroutine dlnrho_dt
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
        i_ekin=0; i_rhom=0; i_ekintot=0; i_rhomin=0; i_rhomax=0
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_ekintot=',i_ekintot
        write(3,*) 'i_ekin=',i_ekin
        write(3,*) 'i_rhom=',i_rhom
        write(3,*) 'i_rhomin=',i_rhomin
        write(3,*) 'i_rhomax=',i_rhomax
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
      if(ip==0) print*,f,topbot
    endsubroutine bc_lnrho_temp_z
!***********************************************************************
    subroutine bc_lnrho_pressure_z(f,topbot)
!
!  dummy routine
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if(ip==0) print*,f,topbot
    endsubroutine bc_lnrho_pressure_z
!***********************************************************************

endmodule Density
