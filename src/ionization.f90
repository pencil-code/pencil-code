! $Id: ionization.f90,v 1.12 2003-02-24 15:43:56 theine Exp $

!  This modules contains the routines for simulation with
!  simple hydrogen ionization.

module Ionization

  use Cparam
  use Cdata
  use Density

  implicit none

  !  secondary parameters calculated in initialize
  real :: TT_ion,lnrho_ion,ss_ion,chiH,m_H

  !  lionization initialized to .true.
  !  it can be reset to .false. in namelist
  logical :: lionization=.true.

  ! input parameters
  integer :: dummy 
  namelist /ionization_init_pars/ lionization

  ! run parameters
  namelist /ionization_run_pars/ lionization

  contains

!***********************************************************************
    subroutine register_ionization()
!
!   2-feb-03/axel: adapted from Interstellar module
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
      if ((ip<=8) .and. lroot) then
        print*, 'Register_ionization'
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: ionization.f90,v 1.12 2003-02-24 15:43:56 theine Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_ionization: nvar > mvar')
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
      logical, save :: first=.true.
      logical :: exist
!
!  ionization parameters
!  since m_e and chiH, as well as hbar are all very small
!  it is better to divide m_e and chiH separately by hbar.
!
      chiH=13.6*eV
      TT_ion=chiH/k_B
      m_H=m_p+m_e
      lnrho_ion=1.5*log((m_e/hbar)*(chiH/hbar)/(2*pi))+log(m_H)
      ss_ion=k_B/m_H
      if(lroot) then
        print*,'initialize_ionization: reference values for ionization'
        print*,'TT_ion',TT_ion
        print*,'lnrho_ion,exp(lnrho_ion)=',lnrho_ion,exp(lnrho_ion)
        print*,'ss_ion=',ss_ion
      endif

    endsubroutine initialize_ionization

!***********************************************************************
    subroutine thermodynamics(lnrho,ss,cs2,TT1,cp1tilde)
!
!  Calculate thermodynamical quantities, cs2, 1/T, and cp1tilde
!  cs2=(dp/drho)_s is the adiabatic sound speed
!  TT1=1/T is the inverse temperature
!  neutral gas: cp1tilde=ss_ion/cp=0.4 ("nabla_ad" maybe better name)
!  in general: cp1tilde=dlnPdS/dlnPdlnrho
!
!   2-feb-03/axel: simple example coded
!
      use Cdata
      use General
      use Sub
!
      real, dimension (nx) :: lnrho,ss,rho1,cs2,TT1,cp1tilde
      real, dimension (nx) :: dlnPdlnrho,dlnPdS,yH,TT,rho,ee,lnTT
      real :: ss0=-5.5542
!
!  calculate cs2, 1/T, and cp1tilde
!
      if(lionization) then
        if(headtt) print*,'thermodynamics: assume cp is not 1'
        call pressure_gradient(lnrho,ss,dlnPdlnrho,dlnPdS,yH,TT)
        TT1=1./TT
!AB: Tobi plans to measure entropy in physical units, and not in kB/mp
!AB: until this is the case, we cannot use T1=1/T, but have to use
!AB: TT1=1/(ss_ion*TT)

!Tobi: entropy is still passed to the ionization routines in units 
!Tobi: of kB/mp, which i think is reasonable.
        TT1=TT1/ss_ion
        cs2=(1.+yH)*ss_ion*TT*dlnPdlnrho
        cp1tilde=dlnPdS/dlnPdlnrho
        if (ldiagnos.and.i_eth/=0) then
          rho=exp(lnrho)
          ee=1.5*(1.+yH)*ss_ion*TT+yH*ss_ion*TT_ion
          call sum_mn_name(rho*ee,i_eth)
        endif
      else
!
!  new case: cp-cv=1, ie cp=2.5
!
        if(headtt) print*,'thermodynamics: assume cp-cv=1 and cp=2.5'
        dlnPdlnrho=gamma
        dlnPdS=gamma1
        lnTT=gamma1*(lnrho+ss-ss0)
        TT=exp(lnTT); TT1=1./TT
        cs2=ss_ion*TT*dlnPdlnrho
        cp1tilde=dlnPdS/dlnPdlnrho
        if (ldiagnos.and.i_eth/=0) then
          rho=exp(lnrho)
          ee=cs2/(gamma*gamma1) !(not correct with ionization)
          call sum_mn_name(rho*ee,i_eth)
        endif
      endif
!
    endsubroutine thermodynamics

!***********************************************************************
    subroutine pressure_gradient(lnrho,ss,dlnPdlnrho,dlnPdS,yH,TT)
!
!   23-feb-03/tobi: rewritten
!
      real, dimension(nx),intent(in)   :: lnrho,ss
      real, dimension(nx), intent(out) :: dlnPdlnrho,dlnPdS,yH,TT
      real, dimension(nx)              :: lnTT_,dlnTT_dy
                                                       ! lnTT_=log(TT/TT_ion)
      real, dimension(nx)              :: f,dfdy,dfdlnrho,dfds
      real, dimension(nx)              :: dydlnrho,dyds
      real, dimension(nx), save        :: yHlast=.5
      integer                          :: i


      do i=1,nx
         yH(i)=rtsafe(lnrho(i),ss(i),yHlast(i))
         yHlast=yH(i)
      enddo

      lnTT_=gamma1*((ss-1.5*(1.-yH)*log(m_H/m_e) &
                        -1.5*yH*log(m_p/m_e)+(1.-yH)*log(1.-yH) &
                        +2.*yH*log(yH))/(1.+yH)+lnrho-lnrho0-2.5)
      f=lnrho_ion-lnrho+1.5*lnTT_-exp(-lnTT_)+log(1.-yH)-2.*log(yH)
      dlnTT_dy=(log(m_H/m_p)-gamma1*(f+exp(-lnTT_))-1.)/(1.+yH)

      dfdy=dlnTT_dy*(1.5+exp(-lnTT_))-1./(1.-yH)-2./yH
      dfdlnrho=gamma1*exp(-lnTT_)
      dfds=(1.+dfdlnrho)/(1.+yH)

      dydlnrho=-dfdlnrho/dfdy
      dyds=-dfds/dfdy

      dlnPdlnrho=1.+dydlnrho/(1.+yH)+dlnTT_dy*dydlnrho+gamma1
      dlnPds=dyds/(1.+yH)+dlnTT_dy*dyds+gamma1/(1.+yH)
      TT=exp(lnTT_)*TT_ion
    endsubroutine pressure_gradient

!***********************************************************************
    function rtsafe(lnrho,ss,yHlast)
!
!   23-feb-03/tobi: minor changes
!
      real, intent (in)    :: lnrho,ss,yHlast
      real                 :: yH0,yH1,dyHold,dyH,fl,fh,f,df,temp,rtsafe
      real, parameter      :: yHacc=1.e-5
      integer              :: i
      integer, parameter   :: maxit=100

      yH1=1.-yHacc
      yH0=yHacc
      dyHold=1.
      dyH=dyHold

      call saha(yH0,lnrho,ss,fh,df)
      if (fh.le.0.) return
      call saha(yH1,lnrho,ss,fl,df)
      if (fl.ge.0.) return
    
      rtsafe=yHlast
      call saha(rtsafe,lnrho,ss,f,df)

      do i=1,maxit
         if (((rtsafe-yH0)*df-f)*((rtsafe-yH1)*df-f).gt.0. &
              .or.abs(2.*f).gt.abs(dyHold*df)) then
            dyHold=dyH
            dyH=.5*(yH0-yH1)
            rtsafe=yH1+dyH
            if (yH1.eq.rtsafe) return
         else
            dyHold=dyH
            dyH=f/df
            temp=rtsafe
            rtsafe=rtsafe-dyH
            if (temp.eq.rtsafe) return
         endif
         if (abs(dyH).lt.yHacc) return
         call saha(rtsafe,lnrho,ss,f,df)
         if (f.lt.0.) then
            yH1=rtsafe
         else
            yH0=rtsafe
         endif
      enddo
      print *,'rtsafe exceeded maximum iterations',f
    endfunction rtsafe

!***********************************************************************
    subroutine saha(yH,lnrho,ss,f,df)
!
!   23-feb-03/tobi: rewritten
!
      real, intent(in)  :: yH,lnrho,ss
      real, intent(out) :: f,df
      real              :: lnTT_,dlnTT_     ! lnTT_=log(TT/TT_ion)

      lnTT_=gamma1*((ss-1.5*(1.-yH)*log(m_H/m_e) &
                        -1.5*yH*log(m_p/m_e)+(1.-yH)*log(1.-yH) &
                        +2.*yH*log(yH))/(1.+yH)+lnrho-lnrho0-2.5)
      f=lnrho_ion-lnrho+1.5*lnTT_-exp(-lnTT_)+log(1.-yH)-2.*log(yH)
      dlnTT_=(log(m_H/m_p)-gamma1*(f+exp(-lnTT_))-1.)/(1.+yH)
      df=dlnTT_*(1.5+exp(-lnTT_))-1./(1.-yH)-2./yH
    endsubroutine saha

endmodule ionization
