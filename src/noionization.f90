! $Id: noionization.f90,v 1.33 2003-06-30 09:33:57 brandenb Exp $

!  Dummy routine for noionization

module Ionization

  use Cparam
  use Cdata
  use Density

  implicit none

  interface thermodynamics              ! Overload the `thermodynamics' function
    module procedure thermodynamics_penc     ! explicit f implicit m,n
    module procedure thermodynamics_arbpenc  ! explicit lnrho(nx), ss(nx)
    module procedure thermodynamics_arbpoint ! explocit lnrho, ss
  endinterface

  ! global ionization parameter for yH (here a scalar set to 0)
!ajwm  real :: yyH=0.  shouldn't be used directly outside module

  ! global parameter for perfect gas EOS for either yH=0 or yH=1
  real :: lnTT0,coef_ss,coef_lr,dlnPdlnrho,dlnPdss
  integer :: l0,l3,m0,m3,n0,n3

  ! secondary parameters calculated in initialize
  double precision :: m_H,m_He,mu,twothirds
  double precision :: TT_ion,TT_ion_,chiH,chiH_,ss_ion,kappa0
  double precision :: lnrho_H,lnrho_e,lnrho_e_,lnrho_p,lnrho_He

  !  lionization initialized to .false.
  !  cannot currently be reset to .true. in namelist
  !  because the namelist is now not even read
  logical :: lionization=.false.,lfixed_ionization=.false.
  real :: yH0=impossible,xHe=0.1

  ! input parameters
  integer :: dummy_ni 
  namelist /ionization_init_pars/ dummy_ni 

  ! run parameters
  namelist /ionization_run_pars/  dummy_ni 

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
      if (.not. first) call stop_it('register_ionization called twice')
      first = .false.
!
      iyH = 0
      iTT = 0

      if ((ip<=8) .and. lroot) then
        print*, 'register_ionization: ionization nvar = ', nvar
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: noionization.f90,v 1.33 2003-06-30 09:33:57 brandenb Exp $")
!
!  Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('Register_ionization: naux > maux')
      endif
!
    endsubroutine register_ionization
!***********************************************************************
    subroutine initialize_ionization()
!
!  Perform any post-parameter-read initialization, e.g. set derived 
!  parameters.
!
!   2-feb-03/axel: adapted from Interstellar module
!
      use Cdata
      use General
!
      real :: yH
!
!  ionization parameters
!  since m_e and chiH, as well as hbar are all very small
!  it is better to divide m_e and chiH separately by hbar.
!
      m_H=m_p+m_e
      m_He=3.97153*m_H
      mu=1.+3.97153*xHe  
      chiH=13.6*eV
      chiH_=0.75*eV
      TT_ion=chiH/k_B
      TT_ion_=chiH_/k_B
      lnrho_e=1.5*log((m_e/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu)
      lnrho_H=1.5*log((m_H/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu)
      lnrho_p=1.5*log((m_p/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu)
      lnrho_He=1.5*log((m_He/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu)
      lnrho_e_=1.5*log((m_e/hbar)*(chiH_/hbar)/2./pi)+log(m_H)+log(mu)
      ss_ion=k_B/m_H/mu      ! AKA c_p for noionisation
      kappa0=sigmaH_/m_H/mu
!
!  the following array subscripts may be used to avoid unnecessary
!  calculations in the ghost zones. useful
!  for 1- and 2-dimensional runs with radiation
!
      if (nx>1) then; l0=1; l3=mx; else; l0=l1; l3=l2; endif
      if (ny>1) then; m0=1; m3=my; else; m0=m1; m3=m2; endif
      if (nz>1) then; n0=1; n3=mz; else; n0=n1; n3=n2; endif
!
      if(lroot) then
        print*,'initialize_ionization: reference values for ionization'
        print*,'TT_ion,lnrho_e,ss_ion=',TT_ion,lnrho_e,ss_ion
      endif
!
      yH=amin1(amax1(yH0,1e-5),1.-1e-5)
      lnTT0=log(TT_ion)+(2./3.)*((-1.5*(1.-yH)*log(m_H/m_e) &
                        -1.5*yH*log(m_p/m_e)-1.5*xHe*log(m_He/m_e) &
                        +(1.-yH)*log(1.-yH)+2.*yH*log(yH)+xHe*log(xHe)) &
                        /(1.+yH+xHe)-lnrho_e-2.5)

      dlnPdlnrho=5./3.    ! gamma?
      dlnPdss=(2./3.)/(1.+yH0+xHe)  !(gamma - 1) / (\mu_{effective}/\mu) ?
!
      coef_lr=dlnPdlnrho-1.
      coef_ss=dlnPdss/ss_ion 
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
    subroutine ionset(f,ss,lnrho,yH,TT)
!
!   set degree of ionization and temperature
!   This routine is called from entropy.f90 and operates on a pencil.
!
!   15-jun-03/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (nx), intent(in) :: ss,lnrho
      real, dimension (nx), intent(out) :: yH,TT
!
!  calculate data on pencil only
!
      if(lfixed_ionization) then
        if(headtt) print*,'ionset: assume cp is not 1, yH0=',yH0
        TT=exp(coef_ss*ss+coef_lr*lnrho+lnTT0)
!ajwm but where does the reference density come in to this? 
        yH=yH0
      else
        if(headtt) print*,'ionset: assume cp=1'
        TT=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)/gamma1
!ajwm - NEEDS TO BE DIVIDED BY the effective c_p ...
      endif
!
    endsubroutine ionset
!***********************************************************************
    subroutine output_ionization(lun)
      integer, intent(in) :: lun
      if(ip==0) print*,lun  !(keep compiler quiet)
    endsubroutine output_ionization
!***********************************************************************
    subroutine thermodynamics_penc(f,TT1,cs2,cp1tilde,ee,yH)
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
      use Cdata
      use General
      use Sub
      use Density, only:cs20,lnrho0,gamma
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (nx), intent(out), optional :: cs2,TT1,cp1tilde,ee,yH
      real, dimension (nx) :: TT,lnrho,ss
      logical :: ldummy
!
      lnrho=f(l1:l2,m,n,ilnrho)
      ss=f(l1:l2,m,n,ient)
      if(lfixed_ionization) then
        if(headtt) print*,'thermodynamics: assume cp is not 1, yH0=',yH0
        TT=exp(coef_ss*ss+coef_lr*lnrho+lnTT0)
!ajwm but where does the reference density come in to this? 

        if (present(cs2))      cs2=(1.+yH0+xHe)*ss_ion*TT*dlnPdlnrho
        if (present(cp1tilde)) cp1tilde=dlnPdss/dlnPdlnrho
        if (present(ee))       ee=1.5*(1.+yH0+xHe)*ss_ion*TT+yH0*ss_ion*TT_ion
      else
        if(headtt) print*,'thermodynamics: assume cp=1'
        TT=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)/gamma1
        if (present(cs2))      cs2=gamma1*TT
        if (present(cp1tilde)) cp1tilde=1.  
        if (present(ee))       ee=cs2/(gamma1*gamma)
      endif
      if (present(TT1)) TT1=1./TT
      if (present(yH)) yH=yH0
!
    endsubroutine thermodynamics_penc
!***********************************************************************
    subroutine thermodynamics_arbpenc(lnrho,ss,TT1,cs2,cp1tilde,ee,yH)
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
      use Cdata
      use General
      use Sub
      use Density, only:cs20,lnrho0,gamma
!
      real, dimension (nx), intent(out), optional :: cs2,TT1,cp1tilde,ee,yH
      real, dimension (nx), intent(in) :: lnrho,ss
      real, dimension (nx) :: TT
      logical :: ldummy
!
      if(lfixed_ionization) then
        if(headtt) print*,'thermodynamics: assume cp is not 1, yH0=',yH0
        TT=exp(coef_ss*ss+coef_lr*lnrho+lnTT0)
!ajwm but where does the reference density come in to this? 

        if (present(cs2))      cs2=(1.+yH0+xHe)*ss_ion*TT*dlnPdlnrho
        if (present(cp1tilde)) cp1tilde=dlnPdss/dlnPdlnrho
        if (present(ee))       ee=1.5*(1.+yH0+xHe)*ss_ion*TT+yH0*ss_ion*TT_ion
      else
        if(headtt) print*,'thermodynamics: assume cp=1'
        TT=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)/gamma1
        if (present(cs2))      cs2=gamma1*TT
        if (present(cp1tilde)) cp1tilde=1.  
        if (present(ee))       ee=cs2/(gamma1*gamma)
      endif
      if (present(TT1)) TT1=1./TT
      if (present(yH)) yH=yH0
!
    endsubroutine thermodynamics_arbpenc
!***********************************************************************
    subroutine thermodynamics_arbpoint(lnrho,ss,TT1,cs2,cp1tilde,ee,yH)
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
      use Cdata
      use General
      use Sub
      use Density, only:cs20,lnrho0,gamma
!
      real, intent(out), optional :: cs2,TT1,cp1tilde,ee, yH
      real, intent(in) :: lnrho,ss
      real :: TT
      logical :: ldummy
!
      if(lfixed_ionization) then
        if(headtt) print*,'thermodynamics: assume cp is not 1, yH0=',yH0
        TT=exp(coef_ss*ss+coef_lr*lnrho+lnTT0)
!ajwm but where does the reference density come in to this? 
        if (present(cs2))      cs2=(1.+yH0+xHe)*ss_ion*TT*dlnPdlnrho
        if (present(cp1tilde)) cp1tilde=dlnPdss/dlnPdlnrho
        if (present(ee))       ee=1.5*(1.+yH0+xHe)*ss_ion*TT+yH0*ss_ion*TT_ion
      else
        if(headtt) print*,'thermodynamics: assume cp=1'
        TT=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)/gamma1
        if (present(cs2))      cs2=gamma1*TT
        if (present(cp1tilde)) cp1tilde=1.  
        if (present(ee))       ee=cs2/(gamma1*gamma)
      endif
      if (present(TT1)) TT1=1./TT
      if (present(yH)) yH=yH0
!
    endsubroutine thermodynamics_arbpoint
!***********************************************************************
    subroutine opacity(f,kaprho)
!
!  calculate opacity
!
!  26-jun-03/tobi: coded
!  30-jun-03/tobi: was a function, needed to turn into subroutine
!
      use Cdata
      use Density, only:cs20,lnrho0,gamma
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz) :: kaprho
      real, dimension (mx) :: lnrho,ss,yH,TT,kappa
!
!  (but TT is not currently consistent with fixed ionization)
!
      lnrho(l0:l3)=f(l0:l3,m,n,ilnrho)
      ss(l0:l3)=f(l0:l3,m,n,ient)
      yH(l0:l3)=yH0
      TT(l0:l3)=cs20*exp(gamma1*(lnrho(l0:l3)-lnrho0)+gamma*ss(l0:l3))/gamma1
!
!  opacity: if lkappa_es then take electron scattering opacity only;
!  otherwise use Hminus opacity (but may need to add kappa_es as well).
!  The factor 2 in front of lnrho takes care of the extra rho factor in kaprho
!
      kaprho(l0:l3,m,n)=.25*exp(2.*lnrho(l0:l3)-lnrho_e_)*(TT_ion_/TT(l0:l3))**1.5 & 
            *exp(TT_ion_/TT(l0:l3))*yH(l0:l3)*(1.-yH(l0:l3))*kappa0
!
    endsubroutine opacity
!***********************************************************************
    subroutine sourcefunction(f,Srad)
!
!  calculate opacity
!
!  26-jun-03/tobi: coded
!  30-jun-03/tobi: was a function, needed to turn into subroutine
!
      use Cdata
      use Density, only:cs20,lnrho0,gamma
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz), intent(out) :: Srad
      real, dimension (mx) :: lnrho,ss,TT
!
!  need to calculate within ghost zones
!  (but TT is not currently consistent with fixed ionization)
!
      lnrho(l0:l3)=f(l0:l3,m,n,ilnrho)
      ss(l0:l3)=f(l0:l3,m,n,ient)
      TT(l0:l3)=cs20*exp(gamma1*(lnrho(l0:l3)-lnrho0)+gamma*ss(l0:l3))/gamma1
!
!  final expression
!
      Srad(l0:l3,m,n)=sigmaSB*TT(l0:l3)**4/pi
!
    endsubroutine sourcefunction
!***********************************************************************
    subroutine ioncalc(f)
!
!   calculate degree of ionization and temperature
!
!   13-jun-03/tobi: coded
!
    real, dimension (mx,my,mz,mvar+maux) :: f
    if(ip==0) print*,f(1,1,1,1)  !(keep compiler quiet)
    endsubroutine ioncalc
!***********************************************************************
endmodule ionization
