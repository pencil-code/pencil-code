! $Id: ionization_fixed.f90,v 1.62 2004-04-19 08:51:02 ajohan Exp $

!
!  Thermodynamics with Fixed ionization fraction
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Ionization

  use Cparam
  use Cdata

  implicit none

  interface eoscalc ! Overload subroutine `eoscalc'
    module procedure eoscalc_farray   ! explicit f implicit m,n
    module procedure eoscalc_point    ! explicit lnrho, ss
    module procedure eoscalc_pencil
  end interface

  interface pressure_gradient ! Overload subroutine `pressure_gradient'
    module procedure pressure_gradient_farray ! explicit f implicit m,n
    module procedure pressure_gradient_point  ! explicit lnrho, ss
  end interface

  interface getentropy ! Overload subroutine `getentropy'
    module procedure getentropy_pencil      ! (dummy routines here --
    module procedure getentropy_point       !  used in noionization.)
  end interface

! integers specifying which independent variables to use in eoscalc
! (only relevant in ionization.f90)
  integer, parameter :: ilnrho_ss=1,ilnrho_ee=2,ilnrho_pp=3,ilnrho_lnTT=4

  ! Constants use in calculation of thermodynamic quantities
  real :: lnTTss,lnTTlnrho,lnTT0

  ! secondary parameters calculated in initialize
  real :: TT_ion,lnTT_ion,TT_ion_,lnTT_ion_
  real :: ss_ion,ee_ion,kappa0,lnchi0,Srad0
  real :: lnrho_H,lnrho_e,lnrho_e_,lnrho_p,lnrho_He
  real :: xHe_term,yH_term,one_yH_term

  !  lionization initialized to .false.
  !  cannot currently be reset to .true. in namelist
  !  because the namelist is now not even read
  real :: yH0=0.,xHe=0.1,xH2=0.

  ! input parameters
  namelist /ionization_init_pars/ yH0,xHe,xH2

  ! run parameters
  namelist /ionization_run_pars/ yH0,xHe,xH2

  ! other variables (needs to be consistent with reset list below)
  integer :: i_yHm=0,i_yHmax=0,i_TTm=0,i_TTmax=0

!ajwm  Moved here from Density.f90
!ajwm  Completely irrelevant to ionization but density and entropy need
!ajwm  reworking to be independent of these things first
  real :: cs0=impossible, rho0=impossible
  real :: cs20=impossible, lnrho0=impossible
  logical :: lcalc_cp = .false.
  real :: gamma=5/3., gamma1
  real :: cp=impossible, cp1=impossible
!ajwm  can't use impossible else it breaks reading param.nml 
  real :: cs2bot=1., cs2top=1. 

  contains

!***********************************************************************
    subroutine register_ionization()
!
!  14-jun-03/axel: adapted from register_ionization
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_ionization: called twice')
      first = .false.
!
      lionization_fixed=.true.
!
      iyH = 0
      ilnTT = 0
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_ionization: ionization nvar = ', nvar
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
          "$Id: ionization_fixed.f90,v 1.62 2004-04-19 08:51:02 ajohan Exp $")
!
!  Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_ionization: naux > maux')
      endif
!
    endsubroutine register_ionization
!*******************************************************************
    subroutine getmu(mu)
!
!  Calculate mean molecular weight of the gas
!
!   12-aug-03/tony: implemented
!   30-mar-04/anders: Added molecular hydrogen to ionization_fixed
!
      real, intent(out) :: mu
!
      mu = (1.+3.97153*xHe)/(1-xH2+xHe)
!
    endsubroutine getmu
!***********************************************************************
    subroutine initialize_ionization()
!
!  Perform any post-parameter-read initialization, e.g. set derived 
!  parameters.
!
!   2-feb-03/axel: adapted from Interstellar module
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real :: mu1yHxHe
!
!  ionization parameters
!  since m_e and chiH, as well as hbar are all very small
!  it is better to divide m_e and chiH separately by hbar.
!
      mu1yHxHe=1.+3.97153*xHe  
      TT_ion=chiH/k_B
      lnTT_ion=log(chiH/k_B)
      TT_ion_=chiH_/k_B
      lnTT_ion_=log(chiH_/k_B)
      lnrho_e=1.5*log((m_e/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_H=1.5*log((m_H/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_p=1.5*log((m_p/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_He=1.5*log((m_He/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_e_=1.5*log((m_e/hbar)*(chiH_/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      ss_ion=k_B/m_H/mu1yHxHe      
      ee_ion=ss_ion*TT_ion
      kappa0=sigmaH_/m_H/mu1yHxHe
      lnchi0=log(kappa0)-log(4.0)
      Srad0=sigmaSB*TT_ion**4/pi
!
      if(lroot) then
        print*,'initialize_ionization: assume cp is not 1, yH0=',yH0
        print*,'initialize_ionization: reference values for ionization'
        print*,'initialize_ionization: TT_ion,lnrho_e,ss_ion=',TT_ion,lnrho_e,ss_ion
      endif
!
      if (yH0>0.) then
        yH_term=yH0*(2*log(yH0)-lnrho_e-lnrho_p)
      elseif (yH0<0.) then
        call stop_it('initialize_ionization: yH0 must not be lower than zero')
      else
        yH_term=0.
      endif
!
      if (yH0<1.) then
        one_yH_term=(1.-yH0)*(log(1.-yH0)-lnrho_H)
      elseif (yH0>1.) then
        call stop_it('initialize_ionization: yH0 must not be greater than one')
      else
        one_yH_term=0.
      endif
!
      if (xHe>0.) then
        xHe_term=xHe*(log(xHe)-lnrho_He)
      elseif (xHe<0.) then
        call stop_it('initialize_ionization: xHe lower than zero makes no sense')
      else
        xHe_term=0.
      endif
!
! Complain if xH2 not between 0 and 0.5
!
      if (xH2 < 0. .or. xH2 > 0.5) &
          call stop_it('initialize_ionization: xH2 must be <= 0.5 and >= 0.0')
!
! Set the reference sound speed (used for noionisation to impossible)
!
      lnTTss=(2./3.)/(1.+yH0+xHe-xH2)/ss_ion
      lnTTlnrho=2./3.
!
      lnTT0=lnTT_ion+(2./3.)*((yH_term+one_yH_term+xHe_term)/ &
          (1+yH0+xHe-xH2)-2.5)
!      
      if(lroot) then
        print*,'initialize_ionization: reference values for ionization'
        print*,'initialize_ionization: TT_ion,ss_ion,kappa0=', &
                TT_ion,ss_ion,kappa0
        print*,'initialize_ionization: lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_=', &
                lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_
      endif
!
!  write scale non-free constants to file; to be read by idl
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro')
        write (1,*) 'TT_ion=',TT_ion
        write (1,*) 'TT_ion_=',TT_ion_
        write (1,*) 'lnrho_e=',lnrho_e
        write (1,*) 'lnrho_H=',lnrho_H
        write (1,*) 'lnrho_p=',lnrho_p
        write (1,*) 'lnrho_He=',lnrho_He
        write (1,*) 'lnrho_e_=',lnrho_e_
        write (1,*) 'ss_ion=',ss_ion
        write (1,*) 'ee_ion=',ee_ion
        write (1,*) 'kappa0=',kappa0
        write (1,*) 'lnchi0=',lnchi0
        write (1,*) 'Srad0=',Srad0
        write (1,*) 'lnTTss=',lnTTss
        write (1,*) 'lnTTlnrho=',lnTTlnrho
        write (1,*) 'lnTT0=',lnTT0
        close (1)
      endif
!
    endsubroutine initialize_ionization
!*******************************************************************
    subroutine rprint_ionization(lreset,lwrite)
!
!  Writes iyH and ilnTT to index.pro file
!
!  14-jun-03/axel: adapted from rprint_radiation
!
      use Cdata
      use Sub
! 
      logical :: lreset
      logical, optional :: lwrite
      integer :: iname
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_yHmax=0
        i_yHm=0
        i_TTmax=0
        i_TTm=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'rprint_ionization: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'yHm',i_yHm)
        call parse_name(iname,cname(iname),cform(iname),'yHmax',i_yHmax)
        call parse_name(iname,cname(iname),cform(iname),'TTm',i_TTm)
        call parse_name(iname,cname(iname),cform(iname),'TTmax',i_TTmax)
      enddo
!
!  write column where which ionization variable is stored
!
      if (present(lwrite)) then
        if (lwrite) then
          write(3,*) 'i_yHmax=',i_yHmax
          write(3,*) 'i_yHm=',i_yHm
          write(3,*) 'i_TTmax=',i_TTmax
          write(3,*) 'i_TTm=',i_TTm
          write(3,*) 'nname=',nname
          write(3,*) 'iyH=',iyH
          write(3,*) 'ilnTT=',ilnTT
        endif
      endif
!
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_ionization
!***********************************************************************
    subroutine ioninit(f)
!
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
!
      if(ip==0) print*,f  !(keep compiler quiet)
!
    endsubroutine ioninit
!***********************************************************************
    subroutine ioncalc(f)
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if(ip==0) print*,f  !(keep compiler quiet)
!
    endsubroutine ioncalc
!***********************************************************************
    subroutine perturb_energy(lnrho,ee,ss,lnTT,yH)
      real, dimension(nx), intent(in) :: lnrho,ee
      real, dimension(nx), intent(out) :: ss,lnTT,yH
      real, dimension(nx) :: TT
!
        yH=yH0
        TT=(2./3.)*TT_ion*(ee/ee_ion-yH0)/(1.+yH0+xHe-xH2)
        lnTT=log(TT)
        ss=(log(TT)-(lnTTlnrho*lnrho)-lnTT0)/lnTTss
!
    end subroutine perturb_energy
!***********************************************************************
    subroutine perturb_mass(lnrho,pp,ss,lnTT,yH)
      real, dimension(nx), intent(in) :: lnrho,pp
      real, dimension(nx), intent(out) :: ss,lnTT,yH
      real, dimension(nx) :: TT
!      real, dimension(nx) :: yH,K
!
        TT= pp / ((1. + yH0 + xHe - xH2) * ss_ion * exp(lnrho))
        lnTT=log(TT)
        ss=(log(TT)-(lnTTlnrho*lnrho)-lnTT0)/lnTTss
        yH=yH0
!        TT= 1.5 * (EE/exp(lnrho)-yH0*ss_ion*TT_ion ) / &
!              ((1.+yH0+xHe-xH2) * ss_ion)
!        K=exp(lnrho_e-lnrho)*(TT/TT_ion)**1.5*exp(-TT_ion/TT)
!        yH=2./(1.+sqrt(1.+4./K))
!        ss=((1.+yH0+xHe-xH2)*(1.5*log(TT/TT_ion)-lnrho+2.5) - yH_term- &
!            one_yH_term-xHe_term)*ss_ion
!
    end subroutine perturb_mass
!***********************************************************************
    subroutine getdensity(EE,TT,yH,rho)

      use Mpicomm, only: stop_it
      
      real, intent(in) :: EE,TT,yH
      real, intent(out) :: rho
      real :: lnrho
print*,'ss_ion,ee_ion,TT_ion',ss_ion,ee_ion,TT_ion
      lnrho = log(EE) - log(1.5*(1.+yH+xHe-xH2)*ss_ion*TT + yH*ee_ion)

      rho=exp(max(lnrho,-15.))
    end subroutine getdensity
!***********************************************************************
    subroutine getentropy_pencil(lnrho,lnTT,ss)
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension(nx), intent(in) :: lnrho,lnTT
      real, dimension(nx), intent(out) :: ss
!
      call stop_it("ionput_pencil: NOT IMPLEMENTED IN IONIZATION_FIXED")
      if (ip==0) print*,lnrho
      if (ip==0) print*,lnTT
      if (ip==0) ss=0.
!
    endsubroutine getentropy_pencil
!***********************************************************************
    subroutine getentropy_point(lnrho,lnTT,ss)
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, intent(in) :: lnrho,lnTT
      real, intent(out) :: ss
!
      call stop_it("ionput_point: NOT IMPLEMENTED IN IONIZATION_FIXED")
      if (ip==0) print*,lnrho
      if (ip==0) print*,lnTT
      if (ip==0) ss=0.
!
    endsubroutine getentropy_point
!***********************************************************************
    subroutine pressure_gradient_farray(f,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx), intent(out) :: cs2,cp1tilde
      real, dimension(nx) :: lnrho,ss,lnTT
!
      lnrho=f(l1:l2,m,n,ilnrho)
      ss=f(l1:l2,m,n,iss)
      lnTT=lnTTss*ss+lnTTlnrho*lnrho+lnTT0
!
      cs2=gamma*(1+yH0+xHe-xH2)*ss_ion*exp(lnTT)
      cp1tilde=(1-gamma1)/(1+yH0+xHe-xH2)/ss_ion
!
    endsubroutine pressure_gradient_farray
!***********************************************************************
    subroutine pressure_gradient_point(lnrho,ss,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      use Cdata
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: cs2,cp1tilde
      real :: lnTT
!
      lnTT=lnTTss*ss+lnTTlnrho*lnrho+lnTT0
!
      cs2=gamma*(1+yH0+xHe-xH2)*ss_ion*exp(lnTT)
      cp1tilde=(1-gamma1)/(1+yH0+xHe-xH2)/ss_ion
!
    endsubroutine pressure_gradient_point
!***********************************************************************
    subroutine temperature_gradient(f,glnrho,gss,glnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx,3), intent(in) :: glnrho,gss
      real, dimension(nx,3), intent(out) :: glnTT
      integer :: j
!
      do j=1,3
        glnTT(:,j)=(2.0/3.0)*(glnrho(:,j)+gss(:,j)/ss_ion/(1+yH0+xHe-xH2))
      enddo
!
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_hessian(f,hlnrho,hss,hlnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally hlnPP and hlnTT
!   hP/rho=cs2*(hlnrho+cp1tilde*hss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx,3,3), intent(in) :: hlnrho,hss
      real, dimension(nx,3,3), intent(out) :: hlnTT
      integer :: i,j
!
      do j=1,3
      do i=1,3
        hlnTT(:,i,j)=(2.0/3.0)*(hlnrho(:,i,j)+hss(:,i,j)/ss_ion/(1+yH0+xHe-xH2))
      enddo
      enddo
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,ss,yH,lnTT,ee,pp,lnchi)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!
      use Cdata
      use Sub
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(out), optional :: lnrho,ss,yH,lnTT
      real, dimension(psize), intent(out), optional :: ee,pp,lnchi
      real, dimension(psize) :: lnrho_,ss_,lnTT_,TT_,yH_
!
      select case (psize)

      case (nx)
        lnrho_=f(l1:l2,m,n,ilnrho)
        ss_=f(l1:l2,m,n,iss)

      case (mx)
        lnrho_=f(:,m,n,ilnrho)
        ss_=f(:,m,n,iss)

      case default
        call stop_it("eoscalc: no such pencil size")

      end select

      lnTT_=lnTTss*ss_+lnTTlnrho*lnrho_+lnTT0
      TT_=exp(lnTT_)
      yH_=yH0
!
      if (present(lnrho)) lnrho=lnrho_
      if (present(ss))    ss=ss_
      if (present(yH))    yH=yH_
      if (present(lnTT))  lnTT=lnTT_
      if (present(ee))    ee=1.5*(1+yH_+xHe-xH2)*ss_ion*TT_+yH_*ss_ion*TT_ion
      if (present(pp))    pp=(1+yH_+xHe-xH2)*exp(lnrho_)*TT_*ss_ion
!
!  Hminus opacity
!
      if (present(lnchi)) then
        lnchi=2*lnrho_-lnrho_e_+1.5*(lnTT_ion_-lnTT_) &
             +TT_ion_/TT_+log(yH_)+log(1-yH_)+lnchi0
      endif
!
      if (ldiagnos.and.psize==nx) then
        if (i_yHmax/=0) call max_mn_name(yH_,i_yHmax)
        if (i_yHm/=0) call sum_mn_name(yH_,i_yHm)
        if (i_TTmax/=0) call max_mn_name(TT_,i_TTmax)
        if (i_TTm/=0) call sum_mn_name(TT_,i_TTm)
      endif
!
    endsubroutine eoscalc_farray
!***********************************************************************
    subroutine eoscalc_point(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!
      use Cdata
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp
      real :: lnrho_,ss_,lnTT_,TT_,rho_,ee_,pp_
!
      select case (ivars)

      case (ilnrho_ss)
        lnrho_ = var1
        ss_    = var2
        lnTT_  = lnTTss*ss_+lnTTlnrho*lnrho_+lnTT0
        TT_    = exp(lnTT_)
        rho_   = exp(lnrho_)
        ee_    = 1.5*(1+yH0+xHe-xH2)*ss_ion*TT_+yH0*ee_ion
        pp_    = (1+yH0+xHe-xH2)*rho_*TT_*ss_ion

      case (ilnrho_ee)
        lnrho_ = var1
        ee_    = var2
        TT_    = (2.0/3.0)*TT_ion*(ee_/ee_ion-yH0)/(1+yH0+xHe-xH2)
        lnTT_  = log(TT_)
        ss_    = (lnTT_-(lnTTlnrho*lnrho_)-lnTT0)/lnTTss
        rho_   = exp(lnrho_)
        pp_    = (1+yH0+xHe-xH2)*rho_*TT_*ss_ion

      case (ilnrho_pp)
        lnrho_ = var1
        pp_    = var2
        rho_   = exp(lnrho_)
        TT_    = pp_/((1+yH0+xHe-xH2)*ss_ion*rho_)
        lnTT_  = log(TT_)
        ss_    = (lnTT_-(lnTTlnrho*lnrho_)-lnTT0)/lnTTss
        ee_    = 1.5*(1+yH0+xHe-xH2)*ss_ion*TT_+yH0*ee_ion

      case (ilnrho_lnTT)
        lnrho_ = var1
        lnTT_  = var2
        ss_    = (lnTT_-lnTTlnrho*lnrho_-lnTT0)/lnTTss

     end select

      if (present(lnrho)) lnrho=lnrho_
      if (present(ss)) ss=ss_
      if (present(yH)) yH=yH0
      if (present(lnTT)) lnTT=lnTT_
      if (present(ee)) ee=ee_
      if (present(pp)) pp=pp_
!
    endsubroutine eoscalc_point
!***********************************************************************
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!
      use Cdata
!
      integer, intent(in) :: ivars
      real, dimension(nx), intent(in) :: var1,var2
      real, dimension(nx), intent(out), optional :: lnrho,ss
      real, dimension(nx), intent(out), optional :: yH,lnTT
      real, dimension(nx), intent(out), optional :: ee,pp
      real, dimension(nx) :: lnrho_,ss_,lnTT_,TT_,rho_,ee_,pp_
!
      select case (ivars)

      case (ilnrho_ss)
        lnrho_ = var1
        ss_    = var2
        lnTT_  = lnTTss*ss_+lnTTlnrho*lnrho_+lnTT0
        TT_    = exp(lnTT_)
        rho_   = exp(lnrho_)
        ee_    = 1.5*(1+yH0+xHe-xH2)*ss_ion*TT_+yH0*ee_ion
        pp_    = (1+yH0+xHe-xH2)*rho_*TT_*ss_ion

      case (ilnrho_ee)
        lnrho_ = var1
        ee_    = var2
        TT_    = (2.0/3.0)*TT_ion*(ee/ee_ion-yH0)/(1+yH0+xHe-xH2)
        lnTT_  = log(TT_)
        ss_    = (lnTT_-(lnTTlnrho*lnrho_)-lnTT0)/lnTTss
        rho_   = exp(lnrho_)
        pp_    = (1+yH0+xHe-xH2)*rho_*TT_*ss_ion

      case (ilnrho_pp)
        lnrho_ = var1
        pp_    = var2
        rho_   = exp(lnrho_)
        TT_    = pp_/((1+yH0+xHe-xH2)*ss_ion*rho_)
        lnTT_  = log(TT_)
        ss_    = (lnTT_-(lnTTlnrho*lnrho_)-lnTT0)/lnTTss
        ee_    = 1.5*(1+yH0+xHe-xH2)*ss_ion*TT_+yH0*ee_ion

      case (ilnrho_lnTT)
        lnrho_ = var1
        lnTT_  = var2
        ss_    = (lnTT_-lnTTlnrho*lnrho_-lnTT0)/lnTTss

     end select

      if (present(lnrho)) lnrho=lnrho_
      if (present(ss)) ss=ss_
      if (present(yH)) yH=yH0
      if (present(lnTT)) lnTT=lnTT_
      if (present(ee)) ee=ee_
      if (present(pp)) pp=pp_
!
    endsubroutine eoscalc_pencil
!***********************************************************************
    subroutine read_ionization_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=ionization_init_pars,ERR=99, IOSTAT=iostat)
      else 
        read(unit,NML=ionization_init_pars,ERR=99) 
      endif

99    return
    endsubroutine read_ionization_init_pars
!***********************************************************************
    subroutine write_ionization_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=ionization_init_pars)
    endsubroutine write_ionization_init_pars
!***********************************************************************
    subroutine read_ionization_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=ionization_run_pars,ERR=99, IOSTAT=iostat)
      else 
        read(unit,NML=ionization_run_pars,ERR=99) 
      endif

99    return
    endsubroutine read_ionization_run_pars
!***********************************************************************
    subroutine write_ionization_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=ionization_run_pars)
    endsubroutine write_ionization_run_pars
!***********************************************************************
    subroutine scale_height_xy(radz0,nrad,f,H_xy)
!
!  calculate characteristic scale height for exponential boundary
!  condition in the radiation module
!
      use Gravity
!
      integer, intent(in) :: radz0,nrad
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(mx,my,radz0), intent(out) :: H_xy
      real, dimension(mx,my,radz0) :: lnrho_xy,ss_xy,TT_xy
!
      if (nrad>0) then
        lnrho_xy=f(:,:,n1-radz0:n1-1,ilnrho)
        ss_xy=f(:,:,n1-radz0:n1-1,iss)
      endif
!
      if (nrad<0) then
        lnrho_xy=f(:,:,n2+1:n2+radz0,ilnrho)
        ss_xy=f(:,:,n2+1:n2+radz0,iss)
      endif
!
      TT_xy=exp(lnTTss*ss_xy+lnTTlnrho*lnrho_xy+lnTT0)
      H_xy=(1.+yH0+xHe-xH2)*ss_ion*TT_xy/gravz
!
    endsubroutine scale_height_xy
!***********************************************************************
    subroutine get_soundspeed(TT,cs2)
!
!  Calculate sound speed for given temperature
!
!  20-Oct-03/tobi: coded
!
      use Mpicomm
!
      real, intent(in)  :: TT
      real, intent(out) :: cs2
!
      call stop_it("get_soundspeed: with ionization, lnrho needs to be known here")
!
      if (ip==0) print*,TT       ! (keep compiler quiet)
      if (ip==0) cs2=0.         ! (keep compiler quiet)
!
    end subroutine get_soundspeed
!***********************************************************************
    subroutine isothermal_entropy(f,T0)
!
!  Isothermal stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!  When entropy is present, this module also initializes entropy.
!
!  Sound speed (and hence Temperature), is
!  initialised to the reference value:
!           sound speed: cs^2_0            from start.in  
!           density: rho0 = exp(lnrho0)
!
!  11-jun-03/tony: extracted from isothermal routine in Density module
!                  to allow isothermal condition for arbitrary density
!  17-oct-03/nils: works also with lionization=T
!  18-oct-03/tobi: distributed across ionization modules
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(inout) :: f
      real, intent(in) :: T0
      real, dimension(nx) :: lnrho,ss
!
      do n=n1,n2
      do m=m1,m2
!
        lnrho=f(l1:l2,m,n,ilnrho)
        ss=ss_ion*((1+yH0+xHe-xH2)*(1.5*log(T0/TT_ion)-lnrho+2.5) & 
                   -yH_term-one_yH_term-xHe_term)
        f(l1:l2,m,n,iss)=ss
!
      enddo
      enddo
!
    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine isothermal_lnrho_ss(f,T0,rho0)
!
!  Isothermal stratification for lnrho and ss (for yH=0!)
!
!  Uses T=T0 everywhere in the box and rho=rho0 in the mid-plane
!
!  Currently only works with grav_profile='linear', but can easily be
!  generalised.
!
!  11-feb-04/anders: Programmed more or less from scratch
!
      use Cdata
      use Gravity, only: grav_profile
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mvar+maux), intent(inout) :: f
      real, intent(in) :: T0,rho0
      real, dimension(nx) :: lnrho,ss,lnTT
!
      if (grav_profile /= 'linear') call stop_it &
          ('isothermal_lnrho_ss: Only implemented for linear gravity profile')
!
!  First calculate hydrostatic density stratification when T=T0
!
      do m=m1,m2
        do n=n1,n2
          f(l1:l2,m,n,ilnrho) = &
              -(Omega*z(n))**2/(2*(1.+xHe-xH2)*ss_ion*T0)+log(rho0)
        enddo
      enddo
!
!  Then calculate entropy as a function of T0 and lnrho
!
      do m=m1,m2
        do n=n1,n2
          lnrho=f(l1:l2,m,n,ilnrho)
          lnTT=log(T0)
          call eoscalc_pencil(ilnrho_lnTT,lnrho,lnTT,ss=ss)
          f(l1:l2,m,n,iss) = ss
        enddo
      enddo
!
    endsubroutine isothermal_lnrho_ss
!***********************************************************************
    subroutine bc_ss_flux(f,topbot,hcond0,hcond1,Fbot, FbotKbot, chi, &
                lmultilayer,lcalc_heatcond_constchi)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!
      use Mpicomm, only: stop_it
      use Cdata
!
      real, intent(in) :: Fbot, FbotKbot, hcond0, hcond1, chi
      logical, intent(in) :: lmultilayer, lcalc_heatcond_constchi
      
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      call stop_it("bc_ss_flux: NOT IMPLEMENTED IN IONIZATION_FIXED")
      if (ip==0) print*,f(1,1,1,1),hcond0,hcond1,Fbot, FbotKbot, chi,lmultilayer,lcalc_heatcond_constchi,topbot
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  23-jun-2003/tony: implemented for lionization_fixed
!
      use Mpicomm, only: stop_it
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      call stop_it("bc_ss_temp_old: NOT IMPLEMENTED IN IONIZATION")
      if (ip==0) print*,f(1,1,1,1),topbot
!
    endsubroutine bc_ss_temp_old
!***********************************************************************
    subroutine bc_ss_temp_x(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      call stop_it("bc_ss_temp_x: NOT IMPLEMENTED IN IONIZATION")
      if (ip==0) print*,f(1,1,1,1),topbot
!
    endsubroutine bc_ss_temp_x
!***********************************************************************
    subroutine bc_ss_temp_y(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      call stop_it("bc_ss_temp_y: NOT IMPLEMENTED IN IONIZATION")
      if (ip==0) print*,f(1,1,1,1),topbot
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      call stop_it("bc_ss_temp_z: NOT IMPLEMENTED IN IONIZATION")
      if (ip==0) print*,f(1,1,1,1),topbot
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      call stop_it("bc_ss_stemp_x: NOT IMPLEMENTED IN IONIZATION")
      if (ip==0) print*,f(1,1,1,1),topbot
!
    endsubroutine bc_ss_stemp_x
!***********************************************************************
    subroutine bc_ss_stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
        call stop_it("bc_ss_stemp_y: NOT IMPLEMENTED IN IONIZATION")
      if (ip==0) print*,f(1,1,1,1),topbot
!
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  3-aug-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      call stop_it("bc_ss_stemp_z: NOT IMPLEMENTED IN IONIZATION_FIXED")
      if (ip==0) print*,f(1,1,1,1),topbot
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_energy(f,topbot)
!
!  boundary condition for entropy
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the entropy module
!
      use Mpicomm, only: stop_it
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my) :: cs2_2d
      integer :: i
!
!  The 'ce' boundary condition for entropy makes the energy constant at
!  the boundaries.
!  This assumes that the density is already set (ie density must register
!  first!)
!
      if (ip==0) print*,f(1,1,1,1),topbot
!
    end subroutine bc_ss_energy
!***********************************************************************
endmodule ionization
