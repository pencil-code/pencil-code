! $Id: nodensity.f90,v 1.18 2003-10-24 13:17:31 dobler Exp $

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

  implicit none

  real :: cs0=1., rho0=1., gamma=1., cs20=1., gamma1=0.
  real :: lnrho0
  logical :: lcalc_cp = .false.
  real :: cs2bot=1., cs2top=1.

  integer :: dummy           ! We cannot define empty namelists
  namelist /density_init_pars/ dummy
  namelist /density_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: i_ekin=0,i_rhom=0,i_ekintot=0

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
           "$Id: nodensity.f90,v 1.18 2003-10-24 13:17:31 dobler Exp $")
!
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded 
!
!  dummy
!
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
    subroutine dlnrho_dt(f,df,uu,glnrho,divu,lnrho,shock,gshock)
!
!  continuity equation, dummy routine
!
!   7-jun-02/axel: adapted from density
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu,glnrho,gshock
      real, dimension (nx) :: lnrho,divu,shock
!
!  will be accessed in noentropy
!
      lnrho=0.
      glnrho=0.
!
      if(ip==0) print*,f,df,uu,divu,shock,gshock
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
        i_ekin=0; i_rhom=0; i_ekintot=0
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_ekintot=',i_ekintot
        write(3,*) 'i_ekin=',i_ekin
        write(3,*) 'i_rhom=',i_rhom
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
