! $Id: density.f90,v 1.5 2002-06-08 08:01:16 brandenb Exp $

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
           "$Revision: 1.5 $", &
           "$Date: 2002-06-08 08:01:16 $")
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
    subroutine dlnrho_dt(f,df,uu,divu,sij,lnrho,glnrho,rho1)
!
!  continuity equation
!  calculate dlnrho/dt = - u.gradlnrho - divu
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3,3) :: sij
      real, dimension (nx,3) :: uu,glnrho,sglnrho,del2u,graddivu,fvisc
      real, dimension (nx) :: lnrho,divu,uglnrho,glnrho2
      real, dimension (nx) :: murho1,rho1,del2lnrho,rho
      real :: diffrho
      integer :: i
!
!  define lnrho; calculate density gradient and avection term
!
      lnrho=f(l1:l2,m,n,ilnrho)
      call grad(f,ilnrho,glnrho)
      call dot_mn(uu,glnrho,uglnrho)
!
!  continuity equation
!
      df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)-uglnrho-divu
!
!  mass diffusion, in units of dxmin*cs0
!
      if (cdiffrho /= 0.) then
        diffrho=cdiffrho*dxmin*cs0
        call del2(f,ilnrho,del2lnrho)
        call dot2_mn(glnrho,glnrho2)
        df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+diffrho*(del2lnrho+glnrho2)
        maxdiffus=amax1(maxdiffus,diffrho)
      endif
!
!  calculate viscous force for du/dt equation
!
      if(lhydro) then
        rho1=exp(-lnrho)
        if(ivisc==0) then
          if (headtt) print*,'viscous force: nu*del2v'
          call del2v(f,iuu,del2u)
          fvisc=nu*del2u
          maxdiffus=amax1(maxdiffus,nu)
        elseif(ivisc==1) then
          if (headtt) print*,'viscous force: mu/rho*(del2u+graddivu/3)'
          murho1=(nu*rho0)*rho1  !(=mu/rho)
          call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
          do i=1,3
            fvisc(:,i)=murho1*(del2u(:,i)+.333333*graddivu(:,i))
          enddo
          maxdiffus=amax1(maxdiffus,murho1)
        elseif(ivisc==2) then
          if (headtt) print*,'viscous force: nu*(del2u+graddivu/3+2S.glnrho)'
          call del2v_etc(f,iuu,del2u,GRADDIV=graddivu)
          call multmv_mn(sij,glnrho,sglnrho)
          fvisc=2*nu*sglnrho+nu*(del2u+.333333*graddivu)
          maxdiffus=amax1(maxdiffus,nu)
        else
          if (headtt) print*,'no viscous force'
        endif
        df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+fvisc
      endif
!
!  calculate density diagnostics: mean density
!
      if (ldiagnos) then
        rho=exp(f(l1:l2,m,n,ilnrho))
        if (i_rhom/=0) call sum_mn_name(rho,i_rhom)
      endif

    endsubroutine dlnrho_dt
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
      if(lroot.and.ip<14) print*,'run through parse list'
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
