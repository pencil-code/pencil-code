! $Id: noionization.f90,v 1.39 2003-08-02 22:09:36 theine Exp $

!  Dummy routine for noionization

module Ionization

  use Cparam
  use Cdata
  use Density

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

  ! global ionization parameter for yH (here a scalar set to 0)
  !  real :: yyH=0.  shouldn't be used directly outside module so commented out

  ! global parameter for perfect gas EOS for either yH=0 or yH=1
  real :: lnTT0,coef_ss,coef_lr,dlnPdlnrho,dlnPdss

  ! secondary parameters calculated in initialize
  real :: mu,twothirds
  real :: TT_ion,TT_ion_,ss_ion,kappa0
  real :: lnrho_H,lnrho_e,lnrho_e_,lnrho_p,lnrho_He

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
           "$Id: noionization.f90,v 1.39 2003-08-02 22:09:36 theine Exp $")
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
      ss_ion=k_B/m_H/mu      
      kappa0=sigmaH_/m_H/mu

!  the following array subscripts may be used to avoid unnecessary
!  calculations in the ghost zones. useful
!  for 1- and 2-dimensional runs with radiation
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

      dlnPdlnrho=5./3.    
      dlnPdss=(2./3.)/(1.+yH0+xHe) 
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
    subroutine ioninit(f)
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
      if(ip==0) print*,f(1,1,1,1)  !(keep compiler quiet)
    endsubroutine ioninit
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
    subroutine ioncalc_ss_point(lnrho,TT,ss)
      real,intent(in) :: lnrho,TT
      real, intent(out) :: ss

        ss = ((1. + yH0 + xHe) &
                         * (1.5*log(TT/TT_ion)+lnrho_e-lnrho+2.5)  &
                         +1.5*((1.-yH0)*log(m_H/m_e)+yH0*log(m_p/m_e)+xHe*log(m_He/m_e)) &
                         -(1.-yH0)*log(1.-yH0)-2.*yH0*log(yH0)-xHe*log(xHe)) * ss_ion
    end subroutine ioncalc_ss_point
!***********************************************************************
    subroutine ioncalc_ss_penc(lnrho,TT,ss)
      real, dimension(nx), intent(in) :: lnrho,TT
      real, dimension(nx), intent(out) :: ss

         ss = ((1. + yH0 + xHe) &
                         * (1.5*log(TT/TT_ion)+lnrho_e-lnrho+2.5)  &
                         +1.5*((1.-yH0)*log(m_H/m_e)+yH0*log(m_p/m_e)+xHe*log(m_He/m_e)) &
                         -(1.-yH0)*log(1.-yH0)-2.*yH0*log(yH0)-xHe*log(xHe)) * ss_ion
    end subroutine ioncalc_ss_penc
!***********************************************************************
    subroutine isothermal_density_ion(pot,tmp)
      real, dimension (nx), intent(in) :: pot
      real, dimension (nx), intent(out) :: tmp

      tmp=pot/((1+yH0+xHe)*ss_ion)
    end subroutine isothermal_density_ion
!***********************************************************************
    subroutine ionget_pencil(f,yH,TT)
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(:), intent(inout) :: yH,TT
      real, dimension(size(yH)) :: lnrho,ss
!
      yH=yH0
!
      if (size(lnrho)==nx) lnrho=f(l1:l2,m,n,ilnrho)
      if (size(ss)==nx) ss=f(l1:l2,m,n,ient)
!
      if (size(lnrho)==mx) lnrho=f(:,m,n,ilnrho)
      if (size(ss)==mx) ss=f(:,m,n,ient)
!
      if(lfixed_ionization) then
        TT=exp(coef_ss*ss+coef_lr*lnrho+lnTT0)
      else
        TT=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)/gamma1
      endif
!
    endsubroutine ionget_pencil
!***********************************************************************
    subroutine ionget_point(lnrho,ss,yH,TT)
!
      use Cdata
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: yH,TT
      real :: lnTT_
!
      yH=yH0
!
      if(lfixed_ionization) then
        TT=exp(coef_ss*ss+coef_lr*lnrho+lnTT0)
      else
        TT=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)/gamma1
      endif
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
      use Cdata
      use General
      use Sub
      use Density, only:cs20,lnrho0,gamma
!
      real, dimension(nx), intent(in) :: lnrho,ss,yH,TT
      real, dimension(nx), optional :: cs2,cp1tilde,ee
!
      if(lfixed_ionization) then
        if(headtt) print*,'thermodynamics: assume cp is not 1, yH0=',yH0
        if (present(cs2))      cs2=(1.+yH0+xHe)*ss_ion*TT*dlnPdlnrho
        if (present(cp1tilde)) cp1tilde=dlnPdss/dlnPdlnrho
        if (present(ee))       ee=1.5*(1.+yH0+xHe)*ss_ion*TT+yH0*ss_ion*TT_ion
      else
        if(headtt) print*,'thermodynamics: assume cp=1'
        if (present(cs2))      cs2=gamma1*TT
        if (present(cp1tilde)) cp1tilde=1.  
        if (present(ee))       ee=cs2/(gamma1*gamma)
      endif
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
      use Cdata
      use General
      use Sub
      use Density, only:cs20,lnrho0,gamma
!
      real, intent(in) :: lnrho,ss,yH,TT
      real, optional :: cs2,cp1tilde,ee
!
      if(lfixed_ionization) then
        if(headtt) print*,'thermodynamics: assume cp is not 1, yH0=',yH0
        if (present(cs2))      cs2=(1.+yH0+xHe)*ss_ion*TT*dlnPdlnrho
        if (present(cp1tilde)) cp1tilde=dlnPdss/dlnPdlnrho
        if (present(ee))       ee=1.5*(1.+yH0+xHe)*ss_ion*TT+yH0*ss_ion*TT_ion
      else
        if(headtt) print*,'thermodynamics: assume cp=1'
        if (present(cs2))      cs2=gamma1*TT
        if (present(cp1tilde)) cp1tilde=1.  
        if (present(ee))       ee=cs2/(gamma1*gamma)
      endif
!
    endsubroutine thermodynamics_point
!***********************************************************************
endmodule ionization
