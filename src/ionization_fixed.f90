! $Id: ionization_fixed.f90,v 1.27 2003-10-20 17:21:46 theine Exp $

!  Dummy routine for noionization

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
  use Density

  implicit none

  interface thermodynamics              ! Overload the `thermodynamics' function
    module procedure thermodynamics_pencil   ! explicit f implicit m,n
    module procedure thermodynamics_point    ! explocit lnrho, ss
  end interface

  interface ionget
    module procedure ionget_pencil
    module procedure ionget_point
  end interface

  interface ionput                      ! Overload subroutine ionput
    module procedure ionput_pencil      ! (dummy routines here --
    module procedure ionput_point       !  used in noionization.)
  end interface

  interface perturb_energy              ! Overload subroutine perturb_energy
    module procedure perturb_energy_pencil
    module procedure perturb_energy_point
  end interface

  interface perturb_mass                ! Overload subroutine perturb_energy
    module procedure perturb_mass_pencil
    module procedure perturb_mass_point
  end interface

  ! Constants use in calculation of thermodynamic quantities
  real :: lnTTss,lnTTlnrho,lnTT0
  real :: cs2TT
  real :: eeTT,ee0
  real :: cp1tilde_

  ! secondary parameters calculated in initialize
  real :: TT_ion,TT_ion_,ss_ion,ee_ion,kappa0,Srad0
  real :: lnrho_H,lnrho_e,lnrho_e_,lnrho_p,lnrho_He
  real :: xHe_term,yH_term,one_yH_term

  !  lionization initialized to .false.
  !  cannot currently be reset to .true. in namelist
  !  because the namelist is now not even read
  real :: yH0=.0,xHe=0.1
  logical :: radcalc_test=.false.

  ! input parameters
  namelist /ionization_init_pars/ yH0,xHe,radcalc_test

  ! run parameters
  namelist /ionization_run_pars/ yH0,xHe,radcalc_test

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
      iTT = 0
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_ionization: ionization nvar = ', nvar
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
          "$Id: ionization_fixed.f90,v 1.27 2003-10-20 17:21:46 theine Exp $")
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
!  Calculate average particle mass in the gas relative to
!
!   12-aug-03/tony: implemented
!
      real, intent(out) :: mu
!
      mu=1.+3.97153*xHe
!
! tobi: the real mean molecular weight would be:
!
! mu=(1.+3.97153*xHe)/(1+yH+xHe)
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
      if(headtt) print*,'initialize_ionization: assume cp is not 1, yH0=',yH0
      mu1yHxHe=1.+3.97153*xHe  
      TT_ion=chiH/k_B
      TT_ion_=chiH_/k_B
      lnrho_e=1.5*log((m_e/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_H=1.5*log((m_H/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_p=1.5*log((m_p/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_He=1.5*log((m_He/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_e_=1.5*log((m_e/hbar)*(chiH_/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      ss_ion=k_B/m_H/mu1yHxHe      
      ee_ion=ss_ion*TT_ion
      kappa0=sigmaH_/m_H/mu1yHxHe
      Srad0=sigmaSB*TT_ion**4/pi
!
      if(lroot) then
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
! Set the reference sound speed (used for noionisation to impossible)
!
      lnTTss=(2./3.)/(1.+yH0+xHe)/ss_ion
      lnTTlnrho=2./3.
!
      lnTT0=log(TT_ion)+(2./3.)*((yH_term+one_yH_term+xHe_term)/(1+yH0+xHe)-2.5)
!
      cs2TT=(5./3.)*(1.+yH0+xHe)*ss_ion
      ee0=yH0*ss_ion*TT_ion
      cp1tilde_=(2./5.)/(1.+yH0+xHe)
      eeTT=1.5*(1.+yH0+xHe)*ss_ion
!
      if(lroot) then
        print*,'initialize_ionization: reference values for ionization'
        print*,'initialize_ionization: TT_ion,ss_ion,kappa0=', &
                TT_ion,ss_ion,kappa0
        print*,'initialize_ionization: lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_=', &
                lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_
      endif
!
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
        write (1,*) 'Srad0=',Srad0
        write (1,*) 'lnTTss=',lnTTss
        write (1,*) 'lnTTlnrho=',lnTTlnrho
        write (1,*) 'lnTT0=',lnTT0
      close (1)
!
    endsubroutine initialize_ionization
!*******************************************************************
    subroutine rprint_ionization(lreset)
!
!  Writes iyH and iTT to index.pro file
!
!  14-jun-03/axel: adapted from rprint_radiation
!
      use Cdata
      use Sub
! 
      logical :: lreset
!
!  write column where which ionization variable is stored
!
      write(3,*) 'nname=',nname
      write(3,*) 'iyH=',iyH
      write(3,*) 'iTT=',iTT
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
    subroutine perturb_energy_point(lnrho,ee,ss,TT,yH)
      real,intent(in) :: lnrho,ee
      real, intent(out) :: ss,TT,yH
!      real :: yH,K
!
        yH = yH0
        TT= (ee-ee0) / eeTT
        ss=(log(TT)-(lnTTlnrho*lnrho)-lnTT0)/lnTTss
!        TT= (EE/exp(lnrho)-yH0*ss_ion*TT_ion*2. ) / &
!              (3. * (1.+yH0+xHe) * ss_ion )
!        K=exp(lnrho_e-lnrho)*(TT/TT_ion)**1.5*exp(-TT_ion/TT)
!        yH=2./(1.+sqrt(1.+4./K))
!        ss=((1.+yH0+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5)-yH_term-one_yH_term-xHe_term)*ss_ion
!      
    end subroutine perturb_energy_point
!***********************************************************************
    subroutine perturb_energy_pencil(lnrho,ee,ss,TT,yH)
      real, dimension(nx), intent(in) :: lnrho,ee
      real, dimension(nx), intent(out) :: ss,TT,yH
!      real, dimension(nx) :: yH,K
!
        yH = yH0
        TT= (ee-ee0) / eeTT
        ss=(log(TT)-(lnTTlnrho*lnrho)-lnTT0)/lnTTss

!        TT= 1.5 * (EE/exp(lnrho)-yH0*ss_ion*TT_ion ) / &
!              ((1.+yH0+xHe) * ss_ion)
!        K=exp(lnrho_e-lnrho)*(TT/TT_ion)**1.5*exp(-TT_ion/TT)
!        yH=2./(1.+sqrt(1.+4./K))
!        ss=((1.+yH0+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5)-yH_term-one_yH_term-xHe_term)*ss_ion
!
    end subroutine perturb_energy_pencil
!***********************************************************************
    subroutine perturb_mass_point(lnrho,pp,ss,TT,yH)
      real,intent(in) :: lnrho,pp
      real, intent(out) :: ss,TT,yH
!      real :: yH,K
!
        TT= pp / ((1. + yH0 + xHe) * ss_ion * exp(lnrho))
        ss=(log(TT)-(lnTTlnrho*lnrho)-lnTT0)/lnTTss
        yH=yH0
!        TT= (EE/exp(lnrho)-yH0*ss_ion*TT_ion*2. ) / &
!              (3. * (1.+yH0+xHe) * ss_ion )
!        K=exp(lnrho_e-lnrho)*(TT/TT_ion)**1.5*exp(-TT_ion/TT)
!        yH=2./(1.+sqrt(1.+4./K))
!        ss=((1.+yH0+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5)-yH_term-one_yH_term-xHe_term)*ss_ion
!
    end subroutine perturb_mass_point
!***********************************************************************
    subroutine perturb_mass_pencil(lnrho,pp,ss,TT,yH)
      real, dimension(nx), intent(in) :: lnrho,pp
      real, dimension(nx), intent(out) :: ss,TT,yH
!      real, dimension(nx) :: yH,K
!
        TT= pp / ((1. + yH0 + xHe) * ss_ion * exp(lnrho))
        ss=(log(TT)-(lnTTlnrho*lnrho)-lnTT0)/lnTTss
        yH=yH0
!        TT= 1.5 * (EE/exp(lnrho)-yH0*ss_ion*TT_ion ) / &
!              ((1.+yH0+xHe) * ss_ion)
!        K=exp(lnrho_e-lnrho)*(TT/TT_ion)**1.5*exp(-TT_ion/TT)
!        yH=2./(1.+sqrt(1.+4./K))
!        ss=((1.+yH0+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5)-yH_term-one_yH_term-xHe_term)*ss_ion
!
    end subroutine perturb_mass_pencil
!***********************************************************************
    subroutine getdensity(EE,TT,yH,rho)

      use Mpicomm, only: stop_it
      
      real, intent(in) :: EE,TT,yH
      real, intent(out) :: rho

      rho = EE / ((1.5*(1.+yH+xHe)*TT + yH*TT_ion) * ss_ion) 

    end subroutine getdensity
!***********************************************************************
    subroutine ionget_pencil(f,yH,TT)
!
!  'extract' ionization fraction and temperature from f array pencilwise
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx), intent(out) :: yH,TT
      real, dimension(nx) :: lnrho,ss
!
      lnrho=f(l1:l2,m,n,ilnrho)
      ss=f(l1:l2,m,n,iss)
!
      yH=yH0
      TT=exp(lnTTss*ss+lnTTlnrho*lnrho+lnTT0)
!
    endsubroutine ionget_pencil
!***********************************************************************
    subroutine ionget_point(lnrho,ss,yH,TT)
!
!  'extract' ionization fraction and temperature from f array
!  for an arbitrary point
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: yH,TT
!
      yH=yH0
      TT=exp(lnTTss*ss+lnTTlnrho*lnrho+lnTT0)
!
    endsubroutine ionget_point
!***********************************************************************
    subroutine ionput_pencil(f,yH,TT)
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx), intent(out) :: yH,TT
      real, dimension(nx) :: lnrho,ss
!
      call stop_it("ionput_pencil: NOT IMPLEMENTED IN IONIZATION_FIXED")
!
    endsubroutine ionput_pencil
!***********************************************************************
    subroutine ionput_point(lnrho,ss,yH,TT)
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: yH,TT
!
      call stop_it("ionput_point: NOT IMPLEMENTED IN IONIZATION_FIXED")
!
    endsubroutine ionput_point
!***********************************************************************
    subroutine thermodynamics_pencil(lnrho,ss,yH,TT,cs2,cp1tilde,ee,pp)
!
!  Calculate thermodynamical quantities, cs2, 1/T, and cp1tilde
!  cs2=(dp/drho)_s is the adiabatic sound speed
!  TT1=1/T is the inverse temperature
!  neutral gas: cp1tilde=kB_over_mp/cp=0.4 ("nabla_ad" maybe better name)
!  in general: cp1tilde=dlnPdS/dlnPdlnrho
!
!   2-feb-03/axel: simple example coded
!  15-jun-03/axel: made compatible with current ionization routine
!
      real, dimension(nx), intent(in) :: lnrho,ss,yH,TT
      real, dimension(nx), optional :: cs2,cp1tilde,ee,pp
!
      if (present(cs2))      cs2=cs2TT*TT
      if (present(cp1tilde)) cp1tilde=cp1tilde_
      if (present(ee))       ee=eeTT*TT+ee0
      if (present(pp))       pp=(1.+yH+xHe)*exp(lnrho)*TT*ss_ion
!
    endsubroutine thermodynamics_pencil
!***********************************************************************
    subroutine thermodynamics_point(lnrho,ss,yH,TT,cs2,cp1tilde,ee,pp)
!
!  Calculate thermodynamical quantities, cs2, 1/T, and cp1tilde
!  cs2=(dp/drho)_s is the adiabatic sound speed
!  TT1=1/T is the inverse temperature
!  neutral gas: cp1tilde=kB_over_mp/cp=0.4 ("nabla_ad" maybe better name)
!  in general: cp1tilde=dlnPdS/dlnPdlnrho
!
!   2-feb-03/axel: simple example coded
!  15-jun-03/axel: made compatible with current ionization routine
!
      real, intent(in) :: lnrho,ss,yH,TT
      real, optional :: cs2,cp1tilde,ee,pp
!
      if (present(cs2))      cs2=cs2TT*TT
      if (present(cp1tilde)) cp1tilde=cp1tilde_
      if (present(ee))       ee=eeTT*TT+ee0
      if (present(pp))       pp=(1.+yH+xHe)*exp(lnrho)*TT*ss_ion
!
    endsubroutine thermodynamics_point
!***********************************************************************
    subroutine radcalc(f,kaprho,Srad)
!
!  calculate source function and opacity
!
!  24-mar-03/axel+tobi: coded
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(mx,my,mz), intent(out) :: kaprho,Srad
      real, dimension(mx) :: lnrho,ss,yH,TT
      real :: kx,ky,kz
!
!  test
!
      if(radcalc_test) then
        if(lroot.and.ip<12) print*,'radcalc: put Srad=kaprho=1 (as a test)'
        kx=2*pi/Lx
        ky=2*pi/Ly
        kz=2*pi/Lz
        Srad=1.+.02*spread(spread(cos(kx*x),2,my),3,mz) &
                   *spread(spread(cos(ky*y),1,mx),3,mz) &
                   *spread(spread(cos(kz*z),1,mx),2,my)
        kaprho=2.+spread(spread(cos(2*kx*x),2,my),3,mz) &
                 *spread(spread(cos(2*ky*y),1,mx),3,mz) &
                 *spread(spread(cos(2*kz*z),1,mx),2,my)
        return
      endif
!
!  no test
!
      do n=1,mz
      do m=1,my
!
         lnrho=f(:,m,n,ilnrho)
         ss=f(:,m,n,iss)
         yH=yH0
         TT=exp(lnTTss*ss+lnTTlnrho*lnrho+lnTT0)
!
!  calculate source function
!
         Srad(:,m,n)=Srad0*(TT/TT_ion)**4
!
!  calculate opacity
!
         kaprho(:,m,n)=.25*exp(2.*lnrho-lnrho_e_)*(TT_ion_/TT)**1.5 &
                          *exp(TT_ion_/TT)*yH*(1.-yH)*kappa0
!
      enddo
      enddo
!
    endsubroutine radcalc
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
      H_xy=(1.+yH0+xHe)*ss_ion*TT_xy/gravz
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
        ss=ss_ion*((1+yH0+xHe)*(1.5*log(T0/TT_ion)-lnrho+2.5) & 
                   -yH_term-one_yH_term-xHe_term)
        f(l1:l2,m,n,iss)=ss
!
      enddo
      enddo
!
    endsubroutine isothermal_entropy
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
      use Gravity
!
      real, intent(in) :: Fbot, FbotKbot, hcond0, hcond1, chi
      logical, intent(in) :: lmultilayer, lcalc_heatcond_constchi
      
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my) :: tmp_xy,cs2_xy,rho_xy
      integer :: i
      
!
      call stop_it("bc_ss_flux: NOT IMPLEMENTED IN IONIZATION")
      if(ldebug) print*,'bc_ss_flux: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('strange-bot')
        if(headtt) print*,'bc_ss_flux: hcond0,hcond1=',hcond0,hcond1
        if ((bcz1(ilnrho) /= "a2") .and. (bcz1(ilnrho) /= "a3"))&
             call stop_it("bc_ss_flux: Inconsistent boundary conditions 1.")
        tmp_xy = gamma1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                 * exp(-gamma*f(:,:,n1,iss) &
                       - gamma1*(f(:,:,n1,ilnrho)-lnrho0))
        tmp_xy = Fbot/(hcond0*hcond1) * tmp_xy ! F_heat/(hcond T_0)
        do i=1,nghost
          f(:,:,n1-i,iss) = &
               (2*i*dz*tmp_xy &
                + 2*gamma1*(f(:,:,n1+i,ilnrho)-f(:,:,n1,ilnrho)) &
               )/gamma &
               + f(:,:,n1+i,iss)
        enddo
!
!  bottom boundary
!  ===============
!
      case('bot')
        if (lmultilayer) then
          if(headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
        else
          if(headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0
        endif
!       if(bcz1(ilnrho)/="a2") call stop_it("bc_ss_flux: bad lnrho bc")
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        cs2_xy=cs20*exp(gamma1*(f(:,:,n1,ilnrho)-lnrho0)+gamma*f(:,:,n1,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy. 
!
        if(lcalc_heatcond_constchi) then
          tmp_xy=Fbot/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FbotKbot/cs2_xy
        endif
!
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+gamma1/gamma* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+2*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case('top')
        if(headtt) print*,'bc_ss_flux: hcond0=',hcond0
        if ((bcz2(ilnrho) /= "a2") .and. (bcz2(ilnrho) /= "a3")) &
             call stop_it("bc_ss_flux: Inconsistent boundary conditions 2.")
        tmp_xy = gamma1/cs20 & ! 1/T_0 (i.e. 1/T at boundary)
                 * exp(-gamma*f(:,:,n2,iss) &
                       - gamma1*(f(:,:,n2,ilnrho)-lnrho0))
        tmp_xy = FbotKbot * tmp_xy ! F_heat/(hcond T_0)
        do i=1,nghost
          f(:,:,n2+i,iss) = &
               (-2*i*dz*tmp_xy &
                + 2*gamma1*(f(:,:,n2-i,ilnrho)-f(:,:,n2,ilnrho)) &
               )/gamma &
               + f(:,:,n2-i,iss)
        enddo
      case default
        print*,"bc_ss_flux: invalid argument"
        call stop_it("")
      endselect
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my) :: tmp_xy
      integer :: i
!
      call stop_it("bc_ss_temp_old: NOT IMPLEMENTED IN IONIZATION")
      if(ldebug) print*,'bc_ss_temp_old: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c2' boundary condition (fixed temperature/sound speed) for entropy.
!  This assumes that the density is already set (ie density must register
!  first!)
!  tmp_xy = s(x,y) on the boundary.
!  gamma*s/cp = [ln(cs2/cs20)-(gamma-1)ln(rho/rho0)]
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if ((bcz1(ilnrho) /= "a2") .and. (bcz1(ilnrho) /= "a3")) &
          call stop_it("bc_ss_temp_old: Inconsistent boundary conditions 3.")
        if (lionization.or.lionization_fixed) then
!AB: currently, lionization=.true. regardless of lionization_fixed,
!AB: so ".or.lionization_fixed" is obsolete
           call stop_it("bc_ss_temp_old: NOT IMPLEMENTED FOR IONIZATION CASES")
        else
          if (ldebug) print*, &
                  'bc_ss_temp_old: set bottom temperature: cs2bot=',cs2bot
          if (cs2bot<=0.) &
                print*,'bc_ss_temp_old: cannot have cs2bot<=0'
            tmp_xy = (-gamma1*(f(:,:,n1,ilnrho)-lnrho0) &
                 + alog(cs2bot/cs20)) / gamma
            f(:,:,n1,iss) = tmp_xy
            do i=1,nghost
               f(:,:,n1-i,iss) = 2*tmp_xy - f(:,:,n1+i,iss)
            enddo
         endif
!
!  top boundary
!
      case('top')
        if ((bcz1(ilnrho) /= "a2") .and. (bcz1(ilnrho) /= "a3")) &
          call stop_it("bc_ss_temp_old: Inconsistent boundary conditions 3.")
        if (lionization.or.lionization_fixed) then
           call stop_it("bc_ss_temp_old: NOT IMPLEMENTED FOR IONIZATION CASES")
        else
          if (ldebug) print*, &
                     'bc_ss_temp_old: set top temperature - cs2top=',cs2top
          if (cs2top<=0.) print*, &
                     'bc_ss_temp_old: cannot have cs2top<=0'
  !       if (bcz1(ilnrho) /= "a2") &
  !            call stop_it("BOUNDCONDS: Inconsistent boundary conditions 4.")
          tmp_xy = (-gamma1*(f(:,:,n2,ilnrho)-lnrho0) &
                   + alog(cs2top/cs20)) / gamma
          f(:,:,n2,iss) = tmp_xy
          do i=1,nghost
            f(:,:,n2+i,iss) = 2*tmp_xy - f(:,:,n2-i,iss)
          enddo
        endif
      case default
        print*,"bc_ss_temp_old: invalid argument"
        call stop_it("")
      endselect
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      call stop_it("bc_ss_temp_x: NOT IMPLEMENTED IN IONIZATION")
      if(ldebug) print*,'bc_ss_temp_x: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_x: set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_x: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(l1,:,:,iss) = 0.5*tmp - gamma1/gamma*(f(l1,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
               - gamma1/gamma*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                       'bc_ss_temp_x: set x top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, &
                       'bc_ss_temp_x: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(l2,:,:,iss) = 0.5*tmp - gamma1/gamma*(f(l2,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
               - gamma1/gamma*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
        enddo

      case default
        print*,"bc_ss_temp_x: invalid argument"
        call stop_it("")
      endselect
      

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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      call stop_it("bc_ss_temp_y: NOT IMPLEMENTED IN IONIZATION")
      if(ldebug) print*,'bc_ss_temp_y: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_y: set y bottom temperature - cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_y: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(:,m1,:,iss) = 0.5*tmp - gamma1/gamma*(f(:,m1,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m1-i,:,iss) = -f(:,m1+i,:,iss) + tmp &
               - gamma1/gamma*(f(:,m1+i,:,ilnrho)+f(:,m1-i,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                     'bc_ss_temp_y: set y top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                     'bc_ss_temp_y: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,m2,:,iss) = 0.5*tmp - gamma1/gamma*(f(:,m2,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
               - gamma1/gamma*(f(:,m2-i,:,ilnrho)+f(:,m2+i,:,ilnrho)-2*lnrho0)
        enddo

      case default
        print*,"bc_ss_temp_y: invalid argument"
        call stop_it("")
      endselect
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: tmp
      integer :: i
!
      call stop_it("bc_ss_temp_z: NOT IMPLEMENTED IN IONIZATION")
      if(ldebug) print*,'bc_ss_temp_z: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
       if (lionization.or.lionization_fixed) then
        call stop_it("bc_ss_temp_z: NOT IMPLEMENTED FOR IONISATION CASE")
       else
        if (ldebug) print*, &
                   'bc_ss_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_z: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(:,:,n1,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
               - gamma1/gamma*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
        enddo
     endif
!
!  top boundary
!
      case('top')
       if (lionization.or.lionization_fixed) then
        call stop_it("bc_ss_temp_z: NOT IMPLEMENTED FOR IONISATION CASE")
       else
        if (ldebug) print*, &
                     'bc_ss_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp_z: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,:,n2,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
               - gamma1/gamma*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
        enddo
       endif
      case default
        print*,"bc_ss_temp_z: invalid argument"
        call stop_it("")
      endselect
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i
!
      call stop_it("bc_ss_stemp_x: NOT IMPLEMENTED IN IONIZATION")
      if(ldebug) print*,'bc_ss_stemp_x: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (cs2bot<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2bot<=0'
        do i=1,nghost
          f(l1-i,:,:,iss) = f(l1+i,:,:,iss) &
               + gamma1/gamma*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2top<=0'
        do i=1,nghost
          f(l2+i,:,:,iss) = f(l2-i,:,:,iss) &
               + gamma1/gamma*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
        enddo

      case default
        print*,"bc_ss_stemp_x: invalid argument"
        call stop_it("")
      endselect
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i
!
        call stop_it("bc_ss_stemp_y: NOT IMPLEMENTED IN IONIZATION")
        if(ldebug) print*,'bc_ss_stemp_y: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (cs2bot<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2bot<=0'
        do i=1,nghost
          f(:,m1-i,:,iss) = f(:,m1+i,:,iss) &
               + gamma1/gamma*(f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho))
        enddo
!
!  top boundary
!
      case('top')
        if (cs2top<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2top<=0'
        do i=1,nghost
          f(:,m2+i,:,iss) = f(:,m2-i,:,iss) &
               + gamma1/gamma*(f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho))
        enddo

      case default
        print*,"bc_ss_stemp_y: invalid argument"
        call stop_it("")
      endselect
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i
!
      call stop_it("bc_ss_stemp_z: NOT IMPLEMENTED IN IONIZATION_FIXED")
      if(ldebug) print*,'bc_ss_stemp_z: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case(topbot)
!
!  bottom boundary
!
      case('bot')
        if (lionization.or.lionization_fixed) then
          call stop_it("bc_ss_stemp_z: NOT IMPLEMENTED FOR IONISATION CASE")
        else
          if (cs2bot<=0.) print*, &
                                  'bc_ss_stemp_z: cannot have cs2bot<=0'
          do i=1,nghost
             f(:,:,n1-i,iss) = f(:,:,n1+i,iss) &
                  + gamma1/gamma*(f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho))
          enddo
        endif
!
!  top boundary
!
      case('top')
       if (lionization.or.lionization_fixed) then
        call stop_it("bc_ss_stemp_z: NOT IMPLEMENTED FOR IONISATION CASE")
       else
        if (cs2top<=0.) print*, &
                 'bc_ss_stemp_z: cannot have cs2top<=0'
         do i=1,nghost
           f(:,:,n2+i,iss) = f(:,:,n2-i,iss) &
                + gamma1/gamma*(f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho))
         enddo
        endif
      case default
        print*,"bc_ss_stemp_z: invalid argument"
        call stop_it("")
      endselect
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
    select case(topbot)
!
! Bottom boundary
!
    case('bot')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma1*f(:,:,n1,ilnrho)+gamma*f(:,:,n1,iss))
      do i=1,nghost
         f(:,:,n1-i,iss)=1./gamma*(-gamma1*f(:,:,n1-i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo

!
! Top boundary
!
    case('top')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma1*f(:,:,n2,ilnrho)+gamma*f(:,:,n2,iss))
      do i=1,nghost
         f(:,:,n2+i,iss)=1./gamma*(-gamma1*f(:,:,n2+i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
    case default
      print*,"bc_ss_energy: invalid argument"
      call stop_it("")
    endselect

    end subroutine bc_ss_energy
!***********************************************************************
endmodule ionization
