! $Id: density.f90,v 1.11 2002-06-15 19:29:55 dobler Exp $

module Density

  use Cparam

  implicit none

  integer :: initlnrho=0
  real :: cs0=1., rho0=1., ampllnrho=1., gamma=5./3., widthlnrho=.1, &
          rho_left=1., rho_right=1., cdiffrho=0., &
          cs20, cs2top, gamma1, zinfty

  namelist /density_init_pars/ &
       cs0,rho0,ampllnrho,gamma,initlnrho,widthlnrho, &
       rho_left,rho_right,cs2top

  namelist /density_run_pars/ &
       cs0,rho0,gamma,cdiffrho,cs2top

  ! other variables (needs to be consistent with reset list below)
  integer :: i_eth=0,i_ekin=0,i_rhom=0

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
           "$Revision: 1.11 $", &
           "$Date: 2002-06-15 19:29:55 $")
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
      real, dimension (mx,my,mz) :: xx,yy,zz,tmp,p,pot,prof
      real, dimension (mz) :: stp
      real, dimension (nx,1) :: rmn
      real :: lnrho0
      real :: beta1,lnrhoint,cs2int
!
!  different initializations of lnrho (called from start).
!  If initrho does't match, f=0 is assumed (default).
!
      select case(initlnrho)
      case(0)
        if (lroot) print*,'uniform lnrho'
        f(:,:,:,ilnrho)=0.
      case(1)
        if (gamma1==0.) then
          if (lroot) print*,'lnrho=gravz*zz/cs0^2 (for isothermal atmosphere)'
          f(:,:,:,ilnrho)=(gravz/cs0**2)*zz
        else
          !
          !  To maintain continuity with respect to the isothermal case,
          !  one may want to specify cs20, and so zinfty is calculated from that.
          !  On the other hand, for polytropic atmospheres it may be more
          !  physical to specify zinfty, ie the top of the atmosphere.
          !  This is done if zinfty is different from 0.
          !
          print*,'z=',z
          zinfty=-cs20/(gamma1*gravz)
          if (lroot) print*,'rho=(1-z/zinfty)^gamma1; zinfty=',zinfty
          f(:,:,n1:n2,ilnrho)=alog(1.-zz(:,:,n1:n2)/zinfty)/gamma1
        endif
      case(2)
        if (lroot) print*,'density jump; rho_left,right=',rho_left,rho_right
        if (lroot) print*,'density jump; widthlnrho=',widthlnrho
        prof=.5*(1.+tanh(zz/widthlnrho))
        f(:,:,:,ilnrho)=alog(rho_left)+alog(rho_left/rho_right)*prof

      case (3)
        !
        !  hydrostatic density stratification for isentropic atmosphere
        !
        if (lgravz) then
          if (lroot) print*,'vertical density stratification'
          !        f(:,:,:,ilnrho)=-zz
          ! isentropic case:
          !        zmax = -cs20/gamma1/gravz
          !        print*, 'zmax = ', zmax
          !        f(:,:,:,ilnrho) = 1./gamma1 * alog(abs(1-zz/zmax))
          ! linear entropy gradient;
          f(:,:,:,ilnrho) = -grads0*zz &
                            + 1./gamma1*alog( 1 + gamma1*gravz/grads0/cs20 &
                                                  *(1-exp(-grads0*zz)) )
          if (notanumber(f(:,:,:,ilnrho))) then
            STOP "INIT_HYDRO: Imaginary density values"
          endif
        endif
        !
        if (lgravr) then
          if (lroot) print*,'radial density stratification (assumes s=const)'
! nok     call potential(x(l1:l2),y(m),z(n),rmn,pot) ! gravity potential
          call potential(x,y(m),z(n),rmn,pot) ! gravity potential
!          call potential(rr,pot) ! gravity potential
          call output(trim(directory)//'/pot.dat',pot,1)

          ! lnrho at point where cs=cs0 and s=s0 (assuming s0=0)
          if (gamma /= 1) then
            lnrho0 = alog(cs20/gamma)/gamma1
            f(:,:,:,ilnrho) = lnrho0 +  alog(1 - gamma1/cs20*pot) / gamma1
          else                  ! isothermal
            f(:,:,:,ilnrho) = alog(rho0)
          endif
        endif

      case (4)
        !
        !  piecewise polytropic for solar convection stuff
        !
        if (lroot) print*,'piecewise polytropic vertical stratification'
        ! top region
        if (isothtop /= 0) then
          beta1 = 0.
          f(:,:,:,ilnrho) = gamma*gravz/cs20*(zz-ztop)
          ! unstable region
          lnrhoint =  gamma*gravz/cs20*(z2-ztop)
        else
          beta1 = gamma*gravz/(mpoly2+1)
          tmp = 1 + beta1*(zz-ztop)/cs20
          tmp = max(tmp,epsi)  ! ensure arg to log is positive
          f(:,:,:,ilnrho) = mpoly2*alog(tmp)
          ! unstable region
          lnrhoint =  mpoly2*alog(1 + beta1*(z2-ztop)/cs20)
        endif
        ! (lnrho at layer interface z=z2)
        cs2int = cs20 + beta1*(z2-ztop) ! cs2 at layer interface z=z2
        ! NB: beta1 i not dT/dz, but dcs2/dz = (gamma-1)c_pdT/dz
        beta1 = gamma*gravz/(mpoly0+1)
        tmp = 1 + beta1*(zz-z2)/cs2int
        tmp = max(tmp,epsi)  ! ensure arg to log is positive
        tmp = lnrhoint + mpoly0*alog(tmp)
        ! smoothly blend the solutions for the two regions:
        stp = step(z,z2,whcond)
        p = spread(spread(stp,1,mx),2,my)
        f(:,:,:,ilnrho) = p*f(:,:,:,ilnrho)  + (1-p)*tmp
        ! bottom (stable) region
        lnrhoint = lnrhoint + mpoly0*alog(1 + beta1*(z1-z2)/cs2int)
        cs2int = cs2int + beta1*(z1-z2) ! cs2 at layer interface z=z1
        beta1 = gamma*gravz/(mpoly1+1)
        tmp = 1 + beta1*(zz-z1)/cs2int
        tmp = max(tmp,epsi)  ! ensure arg to log is positive
        tmp = lnrhoint + mpoly1*alog(tmp)
        ! smoothly blend the solutions for the two regions:
        stp = step(z,z1,whcond)
        p = spread(spread(stp,1,mx),2,my)
        f(:,:,:,ilnrho) = p*f(:,:,:,ilnrho)  + (1-p)*tmp
        ! Fix origin of log density
        f(:,:,:,ilnrho) = f(:,:,:,ilnrho) + alog(rho0) 

      case(11)
        !
        !  sound wave (should be consistent with hydro module)
        !
        if (lroot) print*,'x-wave in lnrho; ampllnrho=',ampllnrho
        f(:,:,:,ilnrho)=ampllnrho*sin(xx)

      case(13)
        !
        !  shock tube test (should be consistent with hydro module)
        !  
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
      real, dimension (nx) :: murho1,rho1,del2lnrho
      real :: diffrho
      integer :: i
!
!  define lnrho; calculate density gradient and avection term
!
      if (headtt) print*,'solve dlnrho_dt'
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
        i_eth=0;i_ekin=0;i_rhom=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'eth',i_eth)
        call parse_name(iname,cname(iname),cform(iname),'ekin',i_ekin)
        call parse_name(iname,cname(iname),cform(iname),'rhom',i_rhom)
      enddo
!
!  write column where which magnetic variable is stored
!
      open(3,file='tmp/density.pro')
      write(3,*) 'i_eth=',i_eth
      write(3,*) 'i_ekin=',i_ekin
      write(3,*) 'i_rhom=',i_rhom
      write(3,*) 'nname=',nname
      write(3,*) 'ilnrho=',ilnrho
      close(3)
!
    endsubroutine rprint_density
!***********************************************************************

endmodule Density
