! $Id: noionization.f90,v 1.26 2003-06-16 13:45:37 mee Exp $

!  Dummy routine for noionization

module Ionization

  use Cparam
  use Cdata
  use Density

  implicit none

  ! global ionization parameter for yH (here a scalar set to 0)
  real :: yyH=0.

  ! global parameter for perfect gas EOS for either yH=0 or yH=1
  real :: lnTT0,coef_ss,coef_lr,dlnPdlnrho,dlnPdss

  ! secondary parameters calculated in initialize
  double precision :: m_H,m_He,mu
  double precision :: TT_ion,lnrho_ion,ss_ion,chiH
!ajwm commented unused quantities
!ajwm   double precision :: TT_ion_,lnrho_ion_,kappa0,chiH_

  !  lionization initialized to .false.
  !  cannot currently be reset to .true. in namelist
  !  because the namelist is now not even read
  logical :: lionization=.false.,lfixed_ionization=.false.
  real :: yH0=impossible,fHe=0.1

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

      !if ((ip<=8) .and. lroot) then
        print*, 'register_ionization: ionization nvar = ', nvar
        print*, 'iyH,iTT = ', iyH,iTT
      !endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: noionization.f90,v 1.26 2003-06-16 13:45:37 mee Exp $")
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
      mu=1.+3.97153*fHe  
      chiH=13.6*eV
!ajwm commented unused quantities
!ajwm      chiH_=0.75*eV
      TT_ion=chiH/k_B
!ajwm      TT_ion_=chiH_/k_B
      lnrho_ion=1.5*log((m_e/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu)
!ajwm      lnrho_ion_=1.5*log((m_e/hbar)*(chiH_/hbar)/2./pi)+log(m_H)+log(mu)
      ss_ion=k_B/m_H/mu      ! AKA c_p for noionisation
!ajwm      kappa0=sigmaH_/m_H/mu
      if(lroot) then
        print*,'initialize_ionization: reference values for ionization'
        print*,'TT_ion,lnrho_ion,ss_ion=',TT_ion,lnrho_ion,ss_ion
      endif
!
      yH=amin1(amax1(yH0,1e-5),1.-1e-5)
      lnTT0=log(TT_ion)+(2./3.)*((-1.5*(1.-yH)*log(m_H/m_e) &
                        -1.5*yH*log(m_p/m_e)-1.5*fHe*log(m_He/m_e) &
                        +(1.-yH)*log(1.-yH)+2.*yH*log(yH)+fHe*log(fHe)) &
                        /(1.+yH+fHe)-lnrho_ion-2.5)
        ! lnTT0 AKA. cs20
      dlnPdlnrho=5./3.    ! gamma?
      dlnPdss=(2./3.)/(1.+yH0+fHe)  !(gamma - 1) / (\mu_{effecive}/\mu) ?
!
      coef_lr=dlnPdlnrho-1.
      coef_ss=dlnPdss/ss_ion  !AKA effective c_p
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
    subroutine thermodynamics(lnrho,ss,TT1,cs2,cp1tilde,ee,yH,TT)
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
      real, dimension (nx), intent(in) :: lnrho,ss
      real, dimension (nx), intent(out) :: cs2,TT1,cp1tilde,ee
      real, dimension (nx), intent(in), optional :: yH,TT
      logical :: ldummy
!
      if (.not. present(yH)) ldummy=.true.
      if (.not. present(TT)) ldummy=.true.
!
      TT1=1./TT
      if(lfixed_ionization) then
        if(headtt) print*,'thermodynamics: assume cp is not 1, yH0=',yH0
        cs2=(1.+yH+fHe)*ss_ion*TT*dlnPdlnrho
        cp1tilde=dlnPdss/dlnPdlnrho
        ee=1.5*(1.+yH+fHe)*ss_ion*TT+yH*ss_ion*TT_ion
      else
        if(headtt) print*,'thermodynamics: assume cp=1'
        cs2=gamma1*TT
        cp1tilde=1.  
        ee=cs2/(gamma1*gamma)
      endif
!
    endsubroutine thermodynamics
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
