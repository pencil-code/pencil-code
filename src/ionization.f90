! $Id: ionization.f90,v 1.10 2003-02-23 10:15:30 theine Exp $

!  This modules contains the routines for simulation with
!  simple hydrogen ionization.

module Ionization

  use Cparam
  use Cdata
  use Density

  implicit none

  !  secondary parameters calculated in initialize
  real :: chiH,lnTT_ion,lnrho_ion,ss_ion,twothirds,m_H

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
           "$Id: ionization.f90,v 1.10 2003-02-23 10:15:30 theine Exp $")
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
!
      chiH=13.6*eV
      lnTT_ion=log(chiH/k_B)
      m_H=m_p+m_e
      lnrho_ion=1.5*log((m_e/hbar)*(chiH/hbar)/(2*pi))+log(m_H)
      !ss_ion=k_B/m_H
      ss_ion=1.
      twothirds=0.6666666667
!
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
!
      real, dimension (nx) :: lnrho,ss,cs2,TT1,cp1tilde
      real, dimension (nx) :: yH,dlnPdlnrho,dlnPdS,lnTT,TT
      real                 :: ss0=-5.5542
!
!  calculate cs2, 1/T, and cp1tilde
!
      if(lionization) then
        if(headtt) print*,'thermodynamics: assume cp is not 1'
        call pressure_gradient(lnrho,ss,yH,dlnPdlnrho,dlnPdS,lnTT)
        TT=exp(lnTT); TT1=1./TT
        cs2=(1.+yH)*ss_ion*TT*dlnPdlnrho
        cp1tilde=dlnPdS/dlnPdlnrho
      else
!
!  if ionization turned off, continue assuming cp=1
!
!       if(headtt) print*,'thermodynamics: assume cp=1'
!
!  old case: cp=1
!  
!       cs2=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)
!       TT1=gamma1/cs2            ! 1/(c_p T) = (gamma-1)/cs^2
!       cp1tilde=1.
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
      endif
!
    endsubroutine thermodynamics

!***********************************************************************
    subroutine pressure_gradient(lnrho,ss,yH,dlnPdlnrho,dlnPdS,lnTT)
!
!   23-feb-03/tobi: coded
!
      use Cdata

      real, dimension(nx),intent(in)   :: lnrho,ss
      real, dimension(nx), intent(out) :: yH,dlnPdlnrho,dlnPdS,lnTT
      real, dimension(nx)              :: lnTT0,dlnTT0dy,f
      real, dimension(nx)              :: dfdy,dfdlnrho,dfds
      real, dimension(nx)              :: dydlnrho,dyds
      integer                          :: i

      do i=1,nx
         yH(i)=rtsafe(lnrho(i),ss(i))
      enddo

      lnTT0=twothirds*((ss/ss_ion-1.5*(1.-yH)*log(m_H/m_e) &
                        -1.5*yH*log(m_p/m_e)+(1.-yH)*log(1.-yH) &
                        +2.*yH*log(yH))/(1.+yH) &
                       +lnrho-lnrho0-2.5)
      f=lnrho_ion-lnrho+1.5*lnTT0-exp(-lnTT0)+log(1.-yH)-2.*log(yH)
      dlnTT0dy=(log(m_H/m_p)-twothirds*(f+exp(-lnTT0))-1.)/(1.+yH)

      dfdy=dlnTT0dy*(1.5+exp(-lnTT0))-1./(1.-yH)-2./yH
      dfdlnrho=twothirds*exp(-lnTT0)
      dfds=(1.+dfdlnrho)/((1.+yH)*ss_ion)

      dydlnrho=-dfdlnrho/dfdy
      dyds=-dfds/dfdy

      dlnPdlnrho=1.+dydlnrho/(1.+yH)+dlnTT0dy*dydlnrho+twothirds
      dlnPds=dyds/(1.+yH)+dlnTT0dy*dyds+twothirds/((1.+yH)*ss_ion)
      lnTT=lnTT0+lnTT_ion
    endsubroutine pressure_gradient

!***********************************************************************
    function rtsafe(lnrho,ss)
!
!   3-feb-03/tobi: coded
!
      real, intent (in)  :: lnrho,ss
      real               :: yH_h,yH_l,dyH_old,dyH,fl,fh,f,df,temp,rtsafe
      real, parameter    :: yH_acc=1.e-5
      real, save         :: yH_last=.5
      integer            :: i
      integer, parameter :: maxit=100
!
      yH_l=yH_acc
      yH_h=1.-yH_acc
      dyH_old=1.
      dyH=dyH_old
!
      call saha(yH_l,lnrho,ss,fl,df)
      if (fl.le.0.) then
         rtsafe=yH_l
         yH_last=rtsafe
         return
      endif
      call saha(yH_h,lnrho,ss,fh,df)
      if (fh.ge.0.) then
         rtsafe=yH_h
         yH_last=rtsafe
         return
      endif
!
      rtsafe=yH_last
      call saha(rtsafe,lnrho,ss,f,df)
!
      do i=1,maxit
         if (((rtsafe-yH_l)*df-f)*((rtsafe-yH_h)*df-f).gt.0. &
              .or.abs(2.*f).gt.abs(dyH_old*df)) then
            dyH_old=dyH
            dyH=.5*(yH_h-yH_l)
            rtsafe=yH_h-dyH
            if (yH_h.eq.rtsafe) then
               yH_last=rtsafe
               return
            endif
         else
            dyH_old=dyH
            dyH=f/df
            temp=rtsafe
            rtsafe=rtsafe-dyH
            if (temp.eq.rtsafe) then
               yH_last=rtsafe
               return
            endif
         endif
         if (abs(dyH).lt.yH_acc) then
            yH_last=rtsafe
            return
         endif
         call saha(rtsafe,lnrho,ss,f,df)
         if (f.lt.0.) then
            yH_h=rtsafe
         else
            yH_l=rtsafe
         endif
      enddo
      print *,'rtsafe exceeded maximum iterations',f
    endfunction rtsafe

!***********************************************************************
    subroutine saha(yH,lnrho,ss,f,df)
!
!   23-feb-03/tobi: coded
!
      use Cdata

      real, intent(in)  :: yH,lnrho,ss
      real, intent(out) :: f,df
      real              :: lnTT0,dlnTT0
!
      lnTT0=twothirds*((ss/ss_ion-1.5*(1.-yH)*log(m_H/m_e) &
                         -1.5*yH*log(m_p/m_e)+(1.-yH)*log(1.-yH) &
                         +2.*yH*log(yH))/(1.+yH) &
                        +lnrho-lnrho0-2.5)
      f=lnrho_ion-lnrho+1.5*lnTT0-exp(-lnTT0)+log(1.-yH)-2.*log(yH)
      dlnTT0=(log(m_H/m_p)-twothirds*(f+exp(-lnTT0))-1.)/(1.+yH)
      df=dlnTT0*(1.5+exp(-lnTT0))-1./(1.-yH)-2./yH
    endsubroutine saha

endmodule ionization
