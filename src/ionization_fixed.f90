! $Id: ionization_fixed.f90,v 1.6 2003-08-13 15:30:07 mee Exp $

!  Dummy routine for noionization

module Ionization

  use Cparam

  implicit none

  interface thermodynamics              ! Overload the `thermodynamics' function
    module procedure thermodynamics_pencil   ! explicit f implicit m,n
    module procedure thermodynamics_point    ! explocit lnrho, ss
  end interface

  interface ioncalc_ss                  ! Overload the 'ioncalc_ss' function
    module procedure ioncalc_ss_penc
    module procedure ioncalc_ss_point
  end interface

  interface ionget
    module procedure ionget_pencil
    module procedure ionget_point
  end interface

  ! Constants use in calculation of thermodynamic quantities
  real :: lnTTss,lnTTlnrho,lnTT0
  real :: cs2TT
  real :: eeTT,ee0
  real :: cp1tilde_

  ! secondary parameters calculated in initialize
  real :: TT_ion,TT_ion_,ss_ion,kappa0
  real :: lnrho_H,lnrho_e,lnrho_e_,lnrho_p,lnrho_He
  real :: yHmin,yHmax
  real :: xHe_term,yH_term,one_yH_term

  !  lionization initialized to .false.
  !  cannot currently be reset to .true. in namelist
  !  because the namelist is now not even read
  logical :: lionization=.true.,lionization_fixed=.true.
  real :: yH0=.0,yHacc=1e-5,xHe=0.1

  ! input parameters
  integer :: dummy_ni 
  namelist /ionization_init_pars/ yH0,yHacc,xHe

  ! run parameters
  namelist /ionization_run_pars/ yH0,yHacc,xHe

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
      iyH = 0
      iTT = 0

      if ((ip<=8) .and. lroot) then
        print*, 'register_ionization: ionization nvar = ', nvar
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
          "$Id: ionization_fixed.f90,v 1.6 2003-08-13 15:30:07 mee Exp $")
!
!  Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_ionization: naux > maux')
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
      use Mpicomm, only: stop_it
!
      real :: mu
!
!  ionization parameters
!  since m_e and chiH, as well as hbar are all very small
!  it is better to divide m_e and chiH separately by hbar.
!
      if(headtt) print*,'initialize_ionization: assume cp is not 1, yH0=',yH0
      mu=1.+3.97153*xHe  
      TT_ion=chiH/k_B
      TT_ion_=chiH_/k_B
      lnrho_e=1.5*log((m_e/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu)
      lnrho_H=1.5*log((m_H/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu)
      lnrho_p=1.5*log((m_p/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu)
      lnrho_He=1.5*log((m_He/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu)
      lnrho_e_=1.5*log((m_e/hbar)*(chiH_/hbar)/2./pi)+log(m_H)+log(mu)
      ss_ion=k_B/m_H/mu      
      kappa0=sigmaH_/m_H/mu
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
      yHmin=0.
      yHmax=1.
!
      if(lroot) then
        print*,'initialize_ionization: reference values for ionization'
        print*,'initialize_ionization: TT_ion,ss_ion,kappa0=', &
                TT_ion,ss_ion,kappa0
        print*,'initialize_ionization: lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_=', &
                lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_
      endif
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
      if(ip==0) print*,f(1,1,1,1)  !(keep compiler quiet)
!
    endsubroutine ioninit
!***********************************************************************
    subroutine ioncalc(f)
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if(ip==0) print*,f(1,1,1,1)  !(keep compiler quiet)
!
    endsubroutine ioncalc
!***********************************************************************
    subroutine ioncalc_ss_point(lnrho,TT,ss)
      real,intent(in) :: lnrho,TT
      real, intent(out) :: ss
      real :: yH,K
!
        K=exp(lnrho_e-lnrho)*(TT/TT_ion)**1.5*exp(-TT_ion/TT)
        yH=2./(1.+sqrt(1.+4./K))
        ss=(1.+yH0+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5)-yH_term-one_yH_term-xHe_term
!
    end subroutine ioncalc_ss_point
!***********************************************************************
    subroutine ioncalc_ss_penc(lnrho,TT,ss)
      real, dimension(nx), intent(in) :: lnrho,TT
      real, dimension(nx), intent(out) :: ss
      real, dimension(nx) :: yH,K
!   
        K=exp(lnrho_e-lnrho)*(TT/TT_ion)**1.5*exp(-TT_ion/TT)
        yH=2./(1.+sqrt(1.+4./K))
        ss=(1.+yH0+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5)-yH_term-one_yH_term-xHe_term
!
    end subroutine ioncalc_ss_penc
!***********************************************************************
    subroutine isothermal_density_ion(pot,tmp)
      real, dimension (nx), intent(in) :: pot
      real, dimension (nx), intent(out) :: tmp
!ajwm WRONG
      tmp=pot/((1.+yH0+xHe)*ss_ion)
    end subroutine isothermal_density_ion
!***********************************************************************
    subroutine ionget_pencil(f,yH,TT)
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(:), intent(out) :: yH,TT
      real, dimension(size(yH)) :: lnrho,ss
!
      yH=yH0
!
      if (size(lnrho)==nx) lnrho=f(l1:l2,m,n,ilnrho)
      if (size(ss)==nx) ss=f(l1:l2,m,n,iss)
!
      if (size(lnrho)==mx) lnrho=f(:,m,n,ilnrho)
      if (size(ss)==mx) ss=f(:,m,n,iss)
!

      TT=exp(lnTTss*ss+lnTTlnrho*lnrho+lnTT0)
!
    endsubroutine ionget_pencil
!***********************************************************************
    subroutine ionget_point(lnrho,ss,yH,TT)
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: yH,TT
!
      yH=yH0
!
      TT=exp(lnTTss*ss+lnTTlnrho*lnrho+lnTT0)
!
    endsubroutine ionget_point
!***********************************************************************
    subroutine thermodynamics_pencil(lnrho,ss,yH,TT,cs2,cp1tilde,ee)
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
      real, dimension(nx), optional :: cs2,cp1tilde,ee
!
      if (present(cs2))      cs2=cs2TT*TT
      if (present(cp1tilde)) cp1tilde=cp1tilde_
      if (present(ee))       ee=eeTT*TT+ee0
!
    endsubroutine thermodynamics_pencil
!***********************************************************************
    subroutine thermodynamics_point(lnrho,ss,yH,TT,cs2,cp1tilde,ee)
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
      real, optional :: cs2,cp1tilde,ee
!
      if (present(cs2))      cs2=cs2TT*TT
      if (present(cp1tilde)) cp1tilde=cp1tilde_
      if (present(ee))       ee=eeTT*TT+ee0
!
    endsubroutine thermodynamics_point
!***********************************************************************
endmodule ionization
