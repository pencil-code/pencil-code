! $Id: noionization.f90,v 1.19 2003-06-13 21:33:55 theine Exp $

!  Dummy routine for noionization

module Ionization

  use Cparam
  use Cdata
  use Density

  implicit none

  ! global ionization parameter for yH (here a scalar set to 0)
  real :: yyH=0.

  !  secondary parameters calculated in initialize
  double precision :: m_H,m_He,mu
  double precision :: TT_ion,lnrho_ion,ss_ion,chiH
  double precision :: TT_ion_,lnrho_ion_,kappa0,chiH_

  !  lionization initialized to .false.
  !  cannot currently be reset to .true. in namelist
  !  because the namelist is now not even read
  logical :: lionization=.false.,lfixed_ionization=.false.
  real :: yH0=impossible,fHe=0.

  ! input parameters
  integer :: dummy_ni 
  namelist /ionization_init_pars/ dummy_ni 

  ! run parameters
  namelist /ionization_run_pars/  dummy_ni 

  contains

!***********************************************************************
    subroutine register_ionization()
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
!  ionization parameters
!  since m_e and chiH, as well as hbar are all very small
!  it is better to divide m_e and chiH separately by hbar.
!
      m_H=m_p+m_e
      m_He=3.97153*m_H
      mu=1.+3.97153*fHe
      chiH=13.6*eV
      chiH_=0.75*eV
      TT_ion=chiH/k_B
      TT_ion_=chiH_/k_B
      lnrho_ion=1.5*log((m_e/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu)
      lnrho_ion_=1.5*log((m_e/hbar)*(chiH_/hbar)/2./pi)+log(m_H)+log(mu)
      ss_ion=k_B/m_H/mu
      kappa0=sigmaH_/m_H/mu
      if(lroot) then
        print*,'initialize_ionization: reference values for ionization'
        print*,'TT_ion',TT_ion
        print*,'lnrho_ion,exp(lnrho_ion)=',lnrho_ion,exp(lnrho_ion)
        print*,'ss_ion=',ss_ion
      endif

    endsubroutine initialize_ionization
!***********************************************************************
    subroutine ionfrac(f)
      real, dimension (mx,my,mz,mvar), intent(in) :: f
      if(ip==0) print*,f(1,1,1,1)  !(keep compiler quiet)
    endsubroutine ionfrac
!***********************************************************************
    subroutine output_ionization(lun)
      integer, intent(in) :: lun
      if(ip==0) print*,lun  !(keep compiler quiet)
    endsubroutine output_ionization
!***********************************************************************
    subroutine thermodynamics(lnrho,ss,yH,cs2,TT1,cp1tilde, &
      Temperature,InternalEnergy)
!
!  Calculate thermodynamical quantities, cs2, 1/T, and cp1tilde
!  cs2=(dp/drho)_s is the adiabatic sound speed
!  TT1=1/T is the inverse temperature
!  neutral gas: cp1tilde=kB_over_mp/cp=0.4 ("nabla_ad" maybe better name)
!  in general: cp1tilde=dlnPdS/dlnPdlnrho
!
!   2-feb-03/axel: simple example coded
!
      use Cdata
      use General
      use Sub
      use Density, only:cs20, lnrho0,gamma
!
      real, dimension (nx), intent(out), optional :: Temperature
      real, dimension (nx), intent(out), optional :: InternalEnergy
      real, dimension (nx) :: lnrho,ss,yH,cs2,TT1,cp1tilde
      real, dimension (nx) :: TT,dlnPdlnrho,dlnPdss,rho,ee
!
!  calculate cs2, 1/T, and cp1tilde
!  leave this in, in case we may activate it
!
      if(lfixed_ionization) then
        if(headtt) print*,'thermodynamics: assume cp is not 1, yH0=',yH0
        call ioncalc(lnrho,ss,yH,dlnPdlnrho=dlnPdlnrho, &
                                 dlnPdss=dlnPdss, &
                                 TT=TT)
        TT1=1./TT
        cs2=(1.+yH+fHe)*ss_ion*TT*dlnPdlnrho
        cp1tilde=dlnPdss/dlnPdlnrho
        print*,'ss_ion,TT,dlnPdlnrho=',ss_ion,TT,dlnPdlnrho
        if(present(InternalEnergy)) ee=1.5*(1.+yH+fHe)*ss_ion*TT+yH*ss_ion*TT_ion
        if(present(Temperature)) Temperature=TT
        if(present(InternalEnergy)) InternalEnergy=ee
      else
!
!  if ionization turned off, continue assuming cp=1
!  with IONIZATION=noionization this is the only option
!
        if(headtt) print*,'thermodynamics: assume cp=1'
        cs2=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)
        TT1=gamma1/cs2            ! 1/(c_p T) = (gamma-1)/cs^2
        cp1tilde=1.
        if(present(Temperature)) Temperature=cs2/gamma1
        if(present(InternalEnergy)) InternalEnergy=cs2/(gamma*gamma1)
     endif
!
    endsubroutine thermodynamics
!***********************************************************************
    subroutine ioncalc(lnrho,ss,yH,dlnPdlnrho,dlnPdss,TT,kappa)
!
!   calculates thermodynamic quantities under partial ionization
!
!   29-may-03/axel: coded replacement routine for fixed radiation
!
      real, dimension(nx),intent(in)   :: lnrho,ss,yH
      real, dimension(nx), optional    :: dlnPdlnrho,dlnPdss,TT,kappa
                           intent(out) :: dlnPdlnrho,dlnPdss,TT,kappa
      double precision, dimension(nx)  :: lnTT_  ! lnTT_=log(TT/TT_ion)
      double precision :: fHelogfHe
!
!  EOS: p=(1+yH+f)*rho*s0*T
!
      if (present(TT)) then
         fHelogfHe=fHe*log(fHe)
         if (fHe .eq. 0) fHelogfHe=0
         lnTT_=(2./3.)*((ss/ss_ion-1.5*(1.-yH)*log(m_H/m_e) &
                        -1.5*yH*log(m_p/m_e)-1.5*fHe*log(m_He/m_e) &
                        +(1.-yH)*log(1.-yH)+2.*yH*log(yH)+fHelogfHe) &
                        /(1.+yH+fHe)+lnrho-lnrho_ion-2.5)
      endif
!
!  dlnP/dlnrho=gamma
!
      if (present(dlnPdlnrho)) then
         dlnPdlnrho=5./3.
      endif
!
!  dlnP/dss=(gamma-1)/(1.+yH0+fHe)
!
      if (present(dlnPdss)) then
         dlnPdss=(2./3.)/(1.+yH0+fHe)
      endif
!
!  temperature (from lnTT_)
!
      if (present(TT)) TT=exp(lnTT_)*TT_ion
!
!  if kappa is needed, use electron scattering value, kappa_es
!
      if (present(kappa)) then
         kappa=kappa_es
      endif
!
    endsubroutine ioncalc
!***********************************************************************
endmodule ionization
