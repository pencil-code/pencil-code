! $Id: density.f90,v 1.3 2002-06-04 11:08:37 brandenb Exp $

module Density

  use Cparam

  implicit none

  integer :: initlnrho=0
  real :: cs0=1., rho0=1., ampllnrho=1., gamma=5./3., widthlnrho=.1, &
          rho_left=1., rho_right=1., cdiffrho=0., &
          cs20, cs2top, gamma1

  namelist /density_init_pars/ &
       cs0,rho0,ampllnrho,gamma,initlnrho,widthlnrho, &
       rho_left,rho_right

  namelist /density_run_pars/ &
       cs0,rho0,gamma,cdiffrho

  ! other variables (needs to be consistent with reset list below)
  integer :: i_rhom

  contains

!***********************************************************************
    subroutine register_density()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrho; increase nvar accordingly.
!
!   4-jun-02/axel: adapted from hydro
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
      ldensity = .true.
!
      ilnrho = nvar+1           ! indix to access lam
      nvar = nvar+1             ! added 1 variable
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_dnesity:  nvar = ', nvar
        print*, 'ilnrho = ', ilnrho
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$RCSfile: density.f90,v $", &
           "$Revision: 1.3 $", &
           "$Date: 2002-06-04 11:08:37 $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_hydro: nvar > mvar')
      endif
!
    endsubroutine register_density
!***********************************************************************
    subroutine init_lnrho(f,xx,yy,zz)
!
!  initialise lnrho; called from start.f90
!
!  7-nov-2001/wolf: coded
!
      use Cdata
      use Sub
      use Global
      use Gravity
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,prof
!AB   real, dimension (mz) :: stp
!AB   real, dimension (nx,1) :: rmn
!AB   real :: lnrho0
!AB   real :: beta1,lnrhoint,cs2int
!AB   integer :: i
!
!  different initializations of lnrho (called from start).
!  If initrho does't match, f=0 is assumed (default).
!
      select case(initlnrho)
      case(0)
        if (lroot) print*,'uniform lnrho'
        f(:,:,:,ilnrho)=0.
      case(1)
        if (lroot) print*,'lnrho=gravz*zz/cs0^2 (for isothermal/polytropic)'
        f(:,:,:,ilnrho)=(gravz/cs0**2)*zz
      case(2)
        if (lroot) print*,'density jump; rho_left,right=',rho_left,rho_right
        if (lroot) print*,'density jump; widthlnrho=',widthlnrho
        prof=.5*(1.+tanh(zz/widthlnrho))
        f(:,:,:,ilnrho)=alog(rho_left)+alog(rho_left/rho_right)*prof
!
!  sound wave (should be consistent with hydro module)
!
      case(11)
        if (lroot) print*,'x-wave in lnrho; ampllnrho=',ampllnrho
        f(:,:,:,ilnrho)=ampllnrho*sin(xx)
!
!  shock tube test (should be consistent with hydro module)
!  
      case(13)
        if (lroot) print*,'polytopic standing shock'
        prof=.5*(1.+tanh(xx/widthlnrho))
        f(:,:,:,ilnrho)=alog(rho_left)+(alog(rho_right)-alog(rho_left))*prof
!
!  set default
!
      case default
        if (lroot) print*,'Default: initializing lnrho to zero'
      endselect
!
      if(ip==0) print*,prof,yy
    endsubroutine init_lnrho
!***********************************************************************
    subroutine rprint_density(lreset)
!
!  reads and registers print parameters relevant for compressible part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_rhom=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(ip<15) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rhom',i_rhom)
      enddo
!
!  write column where which magnetic variable is stored
!
      open(3,file='tmp/density.pro')
      write(3,*) 'i_rhom=',i_rhom
      write(3,*) 'nname=',nname
      close(3)
!
    endsubroutine rprint_density
!***********************************************************************

endmodule Density
