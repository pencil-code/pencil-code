! $Id: nodensity.f90,v 1.6 2002-07-11 11:08:50 nilshau Exp $

module Density

  use Cparam

  implicit none

  real :: cs0=1., rho0=1., gamma=1., cs20=1., gamma1=0.
  real :: lnrho0

  integer :: dummy           ! We cannot define empty namelists
  namelist /density_init_pars/ dummy
  namelist /density_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: i_eth=0,i_ekin=0,i_rhom=0

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
           "$Id: nodensity.f90,v 1.6 2002-07-11 11:08:50 nilshau Exp $")
!
    endsubroutine register_density
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
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_lnrho
!***********************************************************************
    subroutine dlnrho_dt(f,df,uu,glnrho,divu,lnrho)
!
!  continuity equation, dummy routine
!
!   7-jun-02/axel: adapted from density
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu,glnrho
      real, dimension (nx) :: lnrho,divu
!
!  will be accessed in noentropy
!
      lnrho=0.
      glnrho=0.
!
      if(ip==0) print*,f,df,uu,divu
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine rprint_density(lreset)
!
!  reads and registers print parameters relevant for compressible part
!
!   8-jun-02/axel: adapted from density
!
      use Cdata
      use Sub
!
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_eth=0;i_ekin=0;i_rhom=0
      endif
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_eth=',i_eth
      write(3,*) 'i_ekin=',i_ekin
      write(3,*) 'i_rhom=',i_rhom
      write(3,*) 'nname=',nname
      write(3,*) 'ilnrho=',ilnrho
!
    endsubroutine rprint_density
!***********************************************************************
    subroutine bc_lnrho_db(f,topbot)
!
!  dummy routine for density boundary condition
!
!  11-jul-2002/nils: coded
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
!
      if (ip==1) print*,topbot,f(1,1,1,1)  !(to keep compiler quiet)
    end subroutine bc_lnrho_db
!***********************************************************************

endmodule Density
