! $Id: noionization.f90,v 1.22 2003-06-15 05:39:21 brandenb Exp $

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
           "$Id: noionization.f90,v 1.22 2003-06-15 05:39:21 brandenb Exp $")
!
!  Check we arn't registering too many auxilliary variables
!
      if (nvar > mvar) then
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
    subroutine thermodynamics(lnrho,ss,TT1,cs2,cp1tilde,yH,TT,ee)
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
      use Density, only:cs20,lnrho0,gamma
!
      real, dimension (nx), intent(in) :: lnrho,ss
      real, dimension (nx), intent(out) :: cs2,TT1,cp1tilde
      real, dimension (nx), intent(out), optional :: ee
      real, dimension (nx), intent(in), optional :: yH,TT
      logical :: ldummy
!
      if (.not. present(yH)) ldummy=.true.
      if (.not. present(TT)) ldummy=.true.
!
      cs2=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)
      TT1=gamma1/cs2            ! 1/(c_p T) = (gamma-1)/cs^2
      cp1tilde=1.
!
      if (ldiagnos) then
        if(present(ee)) ee=cs2/(gamma1*gamma)
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
    real, dimension (mx,my,mz,mvar) :: f
    if(ip==0) print*,f(1,1,1,1)  !(keep compiler quiet)
    endsubroutine ioncalc
!***********************************************************************
endmodule ionization
