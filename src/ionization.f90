! $Id: ionization.f90,v 1.69 2003-08-04 17:09:16 theine Exp $

!  This modules contains the routines for simulation with
!  simple hydrogen ionization.

module Ionization

  use Cparam
  use Cdata
  use Density

  implicit none

  interface thermodynamics              ! Overload subroutine thermodynamics
    module procedure thermodynamics_pencil
    module procedure thermodynamics_point
  endinterface

  interface ioncalc_ss                  ! Overload subroutine ioncalc_ss
    module procedure ioncalc_ss_point
    module procedure ioncalc_ss_pencil
  end interface

  interface ionget                      ! Overload subroutine ionget
    module procedure ionget_pencil
    module procedure ionget_point
  end interface

  !  secondary parameters calculated in initialize
  real, parameter :: twothirds=2./3.
  real :: TT_ion,TT_ion_,ss_ion,kappa0,xHe_term
  real :: lnrho_H,lnrho_e,lnrho_e_,lnrho_p,lnrho_He
  real :: yHmin,yHmax

  !  lionization initialized to .true.
  !  it can be reset to .false. in namelist
  logical :: lionization=.true.,lfixed_ionization=.false.
  real :: yH0=impossible,yHacc=1e-5,xHe=0.1

  ! input parameters
  namelist /ionization_init_pars/ xHe,yH0,yHacc

  ! run parameters
  namelist /ionization_run_pars/ xHe,yH0,yHacc

  contains

!***********************************************************************
    subroutine register_ionization()
!
!   2-feb-03/axel: adapted from Interstellar module
!   13-jun-03/tobi: re-adapted from visc_shock module
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
!  set indices for auxiliary variables
!
      iyH = mvar + naux +1; naux = naux + 1 
      iTT = mvar + naux +1; naux = naux + 1 

      if ((ip<=8) .and. lroot) then
        print*, 'register_ionization: ionization nvar = ', nvar
        print*, 'iyH = ', iyH
        print*, 'iTT = ', iTT
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: ionization.f90,v 1.69 2003-08-04 17:09:16 theine Exp $")
!
!  Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux = ', maux
        call stop_it('register_ionization: naux > maux')
      endif
!
!  Writing files for use with IDL
!
      aux_var(aux_count)=',yh $'
      aux_count=aux_count+1
      if (naux < maux)  aux_var(aux_count)=',TT $'
      if (naux == maux) aux_var(aux_count)=',TT'
      aux_count=aux_count+1
      write(5,*) 'yH = fltarr(mx,my,mz)*one'
      write(5,*) 'TT = fltarr(mx,my,mz)*one'
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
      real :: mu
!
!  ionization parameters
!  since m_e and chiH, as well as hbar are all very small
!  it is better to divide m_e and chiH separately by hbar.
!
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
      if (xHe==0.) then
        xHe_term=0.
      else
        xHe_term=xHe*(log(xHe)-lnrho_He)
      endif
!
      if(lroot) then
        print*,'initialize_ionization: reference values for ionization'
        print*,'TT_ion,ss_ion,kappa0=', &
                TT_ion,ss_ion,kappa0
        print*,'lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_=', &
                lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_
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
    subroutine ioninit(f)
!
!  the ionization fraction has to be set to a value yH0 < yH < yHmax before
!  rtsafe is called for the first time
!
!  12-jul-03/tobi: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
!
      f(:,:,:,iyH)=0.5
!
    endsubroutine ioninit
!***********************************************************************
    subroutine ioncalc(f)
!
!   calculate degree of ionization and temperature
!   This routine is called from equ.f90 and operates on the full 3-D array.
!
!   13-jun-03/tobi: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: lnrho,ss,yH,lnTT_
      integer :: lstart,lstop,mstart,mstop,nstart,nstop
      integer :: l
!
      do n=1,mz
      do m=1,my
      do l=1,mx
         lnrho=f(l,m,n,ilnrho)
         ss=f(l,m,n,ient)
         yH=f(l,m,n,iyH)
         call rtsafe(lnrho,ss,yH)
         f(l,m,n,iyH)=yH
         lnTT_=twothirds*((ss/ss_ion+(1.-yH)*(log(1.-yH)-lnrho_H) &
                           +yH*(2.*log(yH)-lnrho_e-lnrho_p) &
                           +xHe_term)/(1.+yH+xHe) &
                          +lnrho-2.5)
         f(l,m,n,iTT)=exp(lnTT_)*TT_ion
      enddo
      enddo
      enddo
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
      ss=yH*(lnrho_e+lnrho_p-2.*log(yH))+(1.-yH)*(lnrho_H-log(1.-yH))+xHe_term &
         +(1.+yH+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5)
!
    end subroutine ioncalc_ss_point
!***********************************************************************
    subroutine ioncalc_ss_pencil(lnrho,TT,ss)
      real, dimension(nx), intent(in) :: lnrho,TT
      real, dimension(nx), intent(out) :: ss
      real, dimension(nx) :: yH,K
!
      K=exp(lnrho_e-lnrho)*(TT/TT_ion)**1.5*exp(-TT_ion/TT)
      yH=2./(1.+sqrt(1.+4./K))
      ss=yH*(lnrho_e+lnrho_p-2.*log(yH))+(1.-yH)*(lnrho_H-log(1.-yH))+xHe_term &
         +(1.+yH+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5)
!
    end subroutine ioncalc_ss_pencil
!***********************************************************************
    subroutine isothermal_density_ion(pot,tmp)
      real, dimension (nx), intent(in) :: pot
      real, dimension (nx), intent(out) :: tmp
!ajwm WRONG
      tmp=pot/(1+yH0+xHe)/ss_ion
    end subroutine isothermal_density_ion
!***********************************************************************
    subroutine ionget_pencil(f,yH,TT)
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(:), intent(out) :: yH,TT
!
      if (size(yH)==nx) yH=f(l1:l2,m,n,iyH)
      if (size(TT)==nx) TT=f(l1:l2,m,n,iTT)
!
      if (size(yH)==mx) yH=f(:,m,n,iyH)
      if (size(TT)==mx) TT=f(:,m,n,iTT)
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
      yH=0.5
      call rtsafe(lnrho,ss,yH)
      lnTT_=twothirds*((ss/ss_ion+(1.-yH)*(log(1.-yH)-lnrho_H) &
                        +yH*(2.*log(yH)-lnrho_e-lnrho_p) &
                        +xHe_term)/(1.+yH+xHe) &
                       +lnrho-2.5)
      TT=exp(lnTT_)*TT_ion
!
    endsubroutine ionget_point
!***********************************************************************
    subroutine thermodynamics_pencil(lnrho,ss,yH,TT,cs2,cp1tilde,ee)
!
!  Calculate thermodynamical quantities, cs2, 1/T, and cp1tilde
!  cs2=(dp/drho)_s is the adiabatic sound speed
!  neutral gas: cp1tilde=ss_ion/cp=0.4 ("nabla_ad" maybe better name)
!  in general: cp1tilde=dlnPdS/dlnPdlnrho
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!
      use Cdata
      use General
      use Sub
!
      real, dimension(nx), intent(in) :: lnrho,ss,yH,TT
      real, dimension(nx), intent(out), optional :: cs2,cp1tilde,ee
      real, dimension(nx) :: ff,dlnTT_dy,dffdy,dlnTT_dlnrho
      real, dimension(nx) :: dffdlnrho,dydlnrho,dlnPdlnrho
      real, dimension(nx) :: dlnTT_dss,dffdss,dydss,dlnPdss
!
!  calculate 1/T (=TT1), sound speed, and coefficient cp1tilde in
!  the expression (1/rho)*gradp = cs2*(gradlnrho + cp1tilde*gradss)
!
      if (present(cs2).or.present(cp1tilde)) then
         ff=lnrho_e-lnrho+1.5*alog(TT/TT_ion)-TT_ion/TT &
            +log(1.-yH)-2.*log(yH)
         dlnTT_dy=(twothirds*(lnrho_H-lnrho_p-ff-TT_ion/TT)-1)/(1.+yH+xHe)
         dffdy=dlnTT_dy*(1.5+TT_ion/TT)-1./(1.-yH)-2./yH
      endif
      if (present(cs2)) then
         dlnTT_dlnrho=twothirds
         dffdlnrho=twothirds*TT_ion/TT
         dydlnrho=-dffdlnrho/dffdy
         dlnPdlnrho=1.+dydlnrho/(1.+yH+xHe)+dlnTT_dy*dydlnrho+dlnTT_dlnrho
         cs2=(1.+yH+xHe)*ss_ion*TT*dlnPdlnrho
      endif
      if (present(cp1tilde)) then
         dlnTT_dss=twothirds/((1.+yH+xHe)*ss_ion)
         dffdss=(1.+dffdlnrho)/((1.+yH+xHe)*ss_ion)
         dydss=-dffdss/dffdy
         dlnPdss=dydss/(1.+yH+xHe)+dlnTT_dy*dydss+dlnTT_dss
         cp1tilde=dlnPdss/dlnPdlnrho
      endif
!
!  calculate internal energy for calculating thermal energy
!
      if (present(ee)) ee=1.5*(1.+yH+xHe)*ss_ion*TT+yH*ss_ion*TT_ion
!
    endsubroutine thermodynamics_pencil
!***********************************************************************
    subroutine thermodynamics_point(lnrho,ss,yH,TT,cs2,cp1tilde,ee)
!
!  Calculate thermodynamical quantities, cs2, 1/T, and cp1tilde
!  cs2=(dp/drho)_s is the adiabatic sound speed
!  neutral gas: cp1tilde=ss_ion/cp=0.4 ("nabla_ad" maybe better name)
!  in general: cp1tilde=dlnPdS/dlnPdlnrho
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!
      use Cdata
      use General
      use Sub
!
      real, intent(in) ::lnrho,ss,yH,TT
      real, intent(out), optional :: cs2,cp1tilde,ee
      real :: ff,dlnTT_dy,dffdy,dlnTT_dlnrho
      real :: dffdlnrho,dydlnrho,dlnPdlnrho
      real :: dlnTT_dss,dffdss,dydss,dlnPdss
!
!  calculate 1/T (=TT1), sound speed, and coefficient cp1tilde in
!  the expression (1/rho)*gradp = cs2*(gradlnrho + cp1tilde*gradss)
!
      if (present(cs2).or.present(cp1tilde)) then
         ff=lnrho_e-lnrho+1.5*alog(TT/TT_ion)-TT_ion/TT &
            +log(1.-yH)-2.*log(yH)
         dlnTT_dy=(twothirds*(lnrho_H-lnrho_p-ff-TT_ion/TT)-1)/(1.+yH+xHe)
         dffdy=dlnTT_dy*(1.5+TT_ion/TT)-1./(1.-yH)-2./yH
      endif
      if (present(cs2)) then
         dlnTT_dlnrho=twothirds
         dffdlnrho=twothirds*TT_ion/TT
         dydlnrho=-dffdlnrho/dffdy
         dlnPdlnrho=1.+dydlnrho/(1.+yH+xHe)+dlnTT_dy*dydlnrho+dlnTT_dlnrho
         cs2=(1.+yH+xHe)*ss_ion*TT*dlnPdlnrho
      endif
      if (present(cp1tilde)) then
         dlnTT_dss=twothirds/((1.+yH+xHe)*ss_ion)
         dffdss=(1.+dffdlnrho)/((1.+yH+xHe)*ss_ion)
         dydss=-dffdss/dffdy
         dlnPdss=dydss/(1.+yH+xHe)+dlnTT_dy*dydss+dlnTT_dss
         cp1tilde=dlnPdss/dlnPdlnrho
      endif
!
!  calculate internal energy for calculating thermal energy
!
      if (present(ee)) ee=1.5*(1.+yH+xHe)*ss_ion*TT+yH*ss_ion*TT_ion
!
    endsubroutine thermodynamics_point
!***********************************************************************
    subroutine rtsafe(lnrho,ss,yH)
!
!   safe newton raphson algorithm (adapted from NR) !
!   09-apr-03/tobi: changed to subroutine
!
      real, intent (in)    :: lnrho,ss
      real, intent (inout) :: yH
      real                 :: yHmax,dyHold,dyH,fl,fh,f,df,temp
      integer              :: i
      integer, parameter   :: maxit=100

      yHmax=1.-yHacc
      yHmin=yHacc
      dyHold=1.-2.*yHacc
      dyH=dyHold
!
!  return if y is too close to 0 or 1
!
      call saha(yHmin,lnrho,ss,fh,df)
      if (fh.le.0.) then
         yH=yHmin
         return
      endif
      call saha(yHmax,lnrho,ss,fl,df)
      if (fl.ge.0.) then
         yH=yHmax
         return
      endif
!
!  otherwise find root
!
      call saha(yH,lnrho,ss,f,df)
      do i=1,maxit
         if (((yH-yHmin)*df-f)*((yH-yHmax)*df-f).gt.0. &
              .or.abs(2.*f).gt.abs(dyHold*df)) then
            dyHold=dyH
            dyH=.5*(yHmin-yHmax)
            yH=yHmax+dyH
            if (yHmax.eq.yH) return
         else
            dyHold=dyH
            dyH=f/df
            temp=yH
            yH=yH-dyH
            if (temp.eq.yH) return
         endif
         if (abs(dyH).lt.yHacc) return
         call saha(yH,lnrho,ss,f,df)
         if (f.lt.0.) then
            yHmax=yH
         else
            yHmin=yH
         endif
      enddo
      print *,'rtsafe: exceeded maximum iterations',f
!
    endsubroutine rtsafe
!***********************************************************************
    subroutine saha(yH,lnrho,ss,f,df)
!
!   we want to find the root of f
!
!   23-feb-03/tobi: errors fixed
!
      real, intent(in)  :: yH,lnrho,ss
      real, intent(out) :: f,df
      real              :: lnTT_,dlnTT_     ! lnTT_=log(TT/TT_ion)

      lnTT_=twothirds*((ss/ss_ion+(1.-yH)*(log(1.-yH)-lnrho_H) &
                        +yH*(2.*log(yH)-lnrho_e-lnrho_p) &
                        +xHe_term)/(1.+yH+xHe) &
                       +lnrho-2.5)
      f=lnrho_e-lnrho+1.5*lnTT_-exp(-lnTT_)+log(1.-yH)-2.*log(yH)
      dlnTT_=(twothirds*(lnrho_H-lnrho_p-f-exp(-lnTT_))-1)/(1.+yH+xHe)
      df=dlnTT_*(1.5+exp(-lnTT_))-1./(1.-yH)-2./yH
!
    endsubroutine saha
!***********************************************************************
endmodule ionization
