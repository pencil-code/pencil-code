! $Id: ionization.f90,v 1.121 2003-10-20 17:21:46 theine Exp $

!  This modules contains the routines for simulation with
!  simple hydrogen ionization.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 2
!
!***************************************************************

module Ionization

  use Cparam
  use Cdata
  use Density

  implicit none

  interface thermodynamics              ! Overload subroutine thermodynamics
    module procedure thermodynamics_pencil
    module procedure thermodynamics_point
  endinterface

  interface ionget                      ! Overload subroutine ionget
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

  interface perturb_mass                ! Overload subroutine perturb_mass
    module procedure perturb_mass_pencil
    module procedure perturb_mass_point
  end interface

  !  secondary parameters calculated in initialize
  real :: TT_ion,TT_ion_,ss_ion,ee_ion,kappa0,Srad0,xHe_term
  real :: lnrho_H,lnrho_e,lnrho_e_,lnrho_p,lnrho_He
  integer :: l

  !  lionization initialized to .true.
  !  it can be reset to .false. in namelist
  real :: xHe=0.1
  real :: yHacc=1e-5
  real :: yHmin=2*tiny(yHacc), yHmax=1-2*epsilon(yHacc)
  logical :: radcalc_test=.false.

  ! input parameters
  namelist /ionization_init_pars/ xHe,yHacc,yHmin,yHmax,radcalc_test

  ! run parameters
  namelist /ionization_run_pars/ xHe,yHacc,yHmin,yHmax,radcalc_test

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
      if (.not. first) call stop_it('register_ionization: called twice')
      first = .false.
!
      lionization=.true.
      lionization_fixed=.false.
!
!  set indices for auxiliary variables
!
      iyH = mvar + naux +1; naux = naux + 1 
      iTT = mvar + naux +1; naux = naux + 1 

      if ((ip<=8) .and. lroot) then
        print*, 'register_ionization: ionization nvar = ', nvar
        print*, 'register_ionization: iyH = ', iyH
        print*, 'register_ionization: iTT = ', iTT
      endif
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: ionization.f90,v 1.121 2003-10-20 17:21:46 theine Exp $")
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
!*******************************************************************
    subroutine getmu(mu)
!
!  Calculate average particle mass in the gas relative to
!  Note that the particles density is N = nHI + nHII + ne + nHe
!  = (1-y)*nH + y*nH + y*nH + xHe*nH = (1 + yH + xHe) * nH, where
!  nH is the number of protons per cubic centimeter.
!  The number of particles per mole is therefore 1 + yH + xHe.
!  The mass per mole is M=1.+3.97153*xHe, so the mean molecular weight
!  per particle is M/N = (1.+3.97153*xHe)/(1 + yH + xHe).
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
      use General
      use Mpicomm, only: stop_it
!
      real :: mu1yHxHe
!
!  ionization parameters
!  since m_e and chiH, as well as hbar are all very small
!  it is better to divide m_e and chiH separately by hbar.
!
      mu1yHxHe=1+3.97153*xHe
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
      if (xHe>0) then
        xHe_term=xHe*(log(xHe)-lnrho_He)
      elseif (xHe<0) then
        call stop_it('error (initialize_ionization): xHe lower than zero makes no sense')
      else
        xHe_term=0
      endif
!
      if(lroot) then
        print*,'initialize_ionization: reference values for ionization'
        print*,'initialize_ionization: TT_ion,ss_ion,kappa0=', &
                TT_ion,ss_ion,kappa0
        print*,'initialize_ionization: lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_=', &
                lnrho_e,lnrho_H,lnrho_p,lnrho_He,lnrho_e_
      endif
!
!  write ionization parameters to file; to be read by idl
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
!
      do n=1,mz
      do m=1,my
      do l=1,mx
         lnrho=f(l,m,n,ilnrho)
         ss=f(l,m,n,iss)
         yH=f(l,m,n,iyH)
         call rtsafe('lnrho|ss',lnrho,ss,yHmin,yHmax,yH)
         f(l,m,n,iyH)=yH
         lnTT_=(2./3.)*((ss/ss_ion+(1.-yH)*(log(1.-yH)-lnrho_H) &
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
    subroutine perturb_energy_point(lnrho,ee,ss,TT,yH)
!
!  DOCUMENT ME
!
      use Mpicomm, only: stop_it
      
      real,intent(in) :: lnrho,ee
      real, intent(out) :: ss,TT,yH

      yH=0.5*yHmax
      call rtsafe('lnrho|ee',lnrho,ee,yHmin,yHmax*min(ee/ee_ion,1.0),yH)
      
      TT=(ee-yH*ee_ion)/(1.5*(1.+yH+xHe)*ss_ion)
      ss=ss_ion*((1.+yH+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5) &
                 -yH*(2*log(yH)-lnrho_e-lnrho_p) &
                 -(1.-yH)*(log(1.-yH)-lnrho_H)-xHe_term)
       
    end subroutine perturb_energy_point
!***********************************************************************
    subroutine perturb_energy_pencil(lnrho,ee,ss,TT,yH)
!
!  DOCUMENT ME
!
      use Mpicomm, only: stop_it
      
      real, dimension(nx) ,intent(in) :: lnrho,ee
      real, dimension(nx) ,intent(out) :: ss,TT,yH
      real :: temp

      temp=0.5*min(ee(1)/ee_ion,1.0)

      do l=1,nx 
        temp=min(temp,(1-2*epsilon(yHacc))*ee(l)/ee_ion)
        call rtsafe('lnrho|ee',lnrho(l),ee(l),yHmin,yHmax*min(ee(l)/ee_ion,1.0),temp)
        yH(l)=temp
      enddo

      TT=(ee-yH*ee_ion)/(1.5*(1.+yH+xHe)*ss_ion)
      ss=ss_ion*((1.+yH+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5) &
                 -yH*(2*log(yH)-lnrho_e-lnrho_p) &
                 -(1.-yH)*(log(1.-yH)-lnrho_H)-xHe_term)

    end subroutine perturb_energy_pencil
!***********************************************************************
    subroutine perturb_mass_point(lnrho,pp,ss,TT,yH)
!
!  DOCUMENT ME
!
      use Mpicomm, only: stop_it
      
      real,intent(in) :: lnrho,pp
      real, intent(out) :: ss,TT,yH

      yH=0.5
      call rtsafe('lnrho|pp',lnrho,pp,yHmin,yHmax,yH)
      
      TT=pp/((1.+yH+xHe)*ss_ion*exp(lnrho))
      ss=ss_ion*((1.+yH+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5) &
                 -yH*(2*log(yH)-lnrho_e-lnrho_p) &
                 -(1.-yH)*(log(1.-yH)-lnrho_H)-xHe_term)
       
    end subroutine perturb_mass_point
!***********************************************************************
    subroutine perturb_mass_pencil(lnrho,pp,ss,TT,yH)
!
!  DOCUMENT ME
!

      use Mpicomm, only: stop_it
      
      real, dimension(nx) ,intent(in) :: lnrho,pp
      real, dimension(nx) ,intent(out) :: ss,TT,yH
      real :: temp

      temp=0.5

      do l=1,nx 
        call rtsafe('lnrho|pp',lnrho(l),pp(l),yHmin,yHmax,temp)
        yH(l)=temp
      enddo

      TT=pp/((1.+yH+xHe)*ss_ion*exp(lnrho))
      ss=ss_ion*((1.+yH+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5) &
                 -yH*(2*log(yH)-lnrho_e-lnrho_p) &
                 -(1.-yH)*(log(1.-yH)-lnrho_H)-xHe_term)

    end subroutine perturb_mass_pencil
!***********************************************************************
    subroutine getdensity(EE,TT,yH,rho)
!
!  DOCUMENT ME
!

      use Mpicomm, only: stop_it
      
      real, intent(in) :: EE,TT,yH
      real, intent(out) :: rho

      rho=EE/(1.5*(1.+yH+xHe)*ss_ion*TT+yH*ee_ion)

    end subroutine getdensity
!***********************************************************************
    subroutine ionget_pencil(f,yH,TT)
!
!  extract ionization fraction and temperature from f array pencilwise
!
      use Cdata
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx), intent(out) :: yH,TT
!
      yH=f(l1:l2,m,n,iyH)
      TT=f(l1:l2,m,n,iTT)
!
    endsubroutine ionget_pencil
!***********************************************************************
    subroutine ionget_point(lnrho,ss,yH,TT)
!
!   extract ionization fraction and temperature from f array
!   for an arbitrary point
!
      use Cdata
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: yH,TT
      real :: lnTT_
!
      yH=0.5
      call rtsafe('lnrho|ss',lnrho,ss,yHmin,yHmax,yH)

      lnTT_=(2./3.)*((ss/ss_ion+(1.-yH)*(log(1.-yH)-lnrho_H) &
                        +yH*(2.*log(yH)-lnrho_e-lnrho_p) &
                        +xHe_term)/(1.+yH+xHe) &
                       +lnrho-2.5)
      TT=exp(lnTT_)*TT_ion
!
    endsubroutine ionget_point
!***********************************************************************
    subroutine ionput_pencil(f,yH,TT)
!
!  DOCUMENT ME
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension(nx), intent(out) :: yH,TT
!
      call stop_it("ionput_pencil: NOT IMPLEMENTED IN IONIZATION")
!
    endsubroutine ionput_pencil
!***********************************************************************
    subroutine ionput_point(lnrho,ss,yH,TT)
!
!  DOCUMENT ME
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: yH,TT
!
      call stop_it("ionput_point: NOT IMPLEMENTED IN IONIZATION")
!
    endsubroutine ionput_point
!***********************************************************************
    subroutine thermodynamics_pencil(lnrho,ss,yH,TT,cs2,cp1tilde,ee,pp)
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
      real, dimension(nx), intent(out), optional :: cs2,cp1tilde,ee,pp
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
         dlnTT_dy=((2./3.)*(lnrho_H-lnrho_p-ff-TT_ion/TT)-1)/(1.+yH+xHe)
         dffdy=dlnTT_dy*(1.5+TT_ion/TT)-1./(1.-yH)-2./yH
      endif
      if (present(cs2)) then
         dlnTT_dlnrho=(2./3.)
         dffdlnrho=(2./3.)*TT_ion/TT
         dydlnrho=-dffdlnrho/dffdy
         dlnPdlnrho=1.+dydlnrho/(1.+yH+xHe)+dlnTT_dy*dydlnrho+dlnTT_dlnrho
         cs2=(1.+yH+xHe)*ss_ion*TT*dlnPdlnrho
      endif
      if (present(cp1tilde)) then
         dlnTT_dss=(2./3.)/((1.+yH+xHe)*ss_ion)
         dffdss=(1.+dffdlnrho)/((1.+yH+xHe)*ss_ion)
         dydss=-dffdss/dffdy
         dlnPdss=dydss/(1.+yH+xHe)+dlnTT_dy*dydss+dlnTT_dss
         cp1tilde=dlnPdss/dlnPdlnrho
      endif
!
!  calculate internal energy for calculating thermal energy
!
      if (present(ee)) ee=1.5*(1.+yH+xHe)*ss_ion*TT+yH*ss_ion*TT_ion
      if (present(pp)) pp=(1.+yH+xHe)*exp(lnrho)*TT*ss_ion
!
    endsubroutine thermodynamics_pencil
!***********************************************************************
    subroutine thermodynamics_point(lnrho,ss,yH,TT,cs2,cp1tilde,ee,pp)
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
      real, intent(out), optional :: cs2,cp1tilde,ee,pp
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
         dlnTT_dy=((2./3.)*(lnrho_H-lnrho_p-ff-TT_ion/TT)-1)/(1.+yH+xHe)
         dffdy=dlnTT_dy*(1.5+TT_ion/TT)-1./(1.-yH)-2./yH
      endif
      if (present(cs2)) then
         dlnTT_dlnrho=(2./3.)
         dffdlnrho=(2./3.)*TT_ion/TT
         dydlnrho=-dffdlnrho/dffdy
         dlnPdlnrho=1.+dydlnrho/(1.+yH+xHe)+dlnTT_dy*dydlnrho+dlnTT_dlnrho
         cs2=(1.+yH+xHe)*ss_ion*TT*dlnPdlnrho
      endif
      if (present(cp1tilde)) then
         dlnTT_dss=(2./3.)/((1.+yH+xHe)*ss_ion)
         dffdss=(1.+dffdlnrho)/((1.+yH+xHe)*ss_ion)
         dydss=-dffdss/dffdy
         dlnPdss=dydss/(1.+yH+xHe)+dlnTT_dy*dydss+dlnTT_dss
         cp1tilde=dlnPdss/dlnPdlnrho
      endif
!
!  calculate internal energy for calculating thermal energy
!
      if (present(ee)) ee=1.5*(1.+yH+xHe)*ss_ion*TT+yH*ss_ion*TT_ion
      if (present(pp)) pp=(1.+yH+xHe)*exp(lnrho)*TT*ss_ion
!
    endsubroutine thermodynamics_point
!***********************************************************************
    subroutine rtsafe(variables,variable1,variable2,yHlb,yHub,yH,rterror,rtdebug)
!
!   safe newton raphson algorithm (adapted from NR) !
!   09-apr-03/tobi: changed to subroutine
!
      use Cdata
!
      character(len=*), intent(in)   :: variables
      real, intent(in)               :: variable1,variable2
      real, intent(in)               :: yHlb,yHub
      real, intent(inout)            :: yH
      logical, intent(out), optional :: rterror
      logical, intent(in), optional  :: rtdebug
!
      real               :: dyHold,dyH,yHl,yHh,fl,fh,f,df,temp
      integer            :: i
      integer, parameter :: maxit=1000
!
      if (present(rterror)) rterror=.false.
      if (present(rtdebug).and.rtdebug) print*,'rtsafe: i,yH=',0,yH
!
      yHl=yHlb
      yHh=yHub
      dyH=1
      dyHold=dyH
!
      call saha(variables,variable1,variable2,yH,f,df)
!
      do i=1,maxit
        if (present(rtdebug).and.rtdebug) print*,'rtsafe: i,yH=',i,yH
        if (((yH-yHl)*df-f)*((yH-yHh)*df-f)>0.or.abs(2*f)>abs(dyHold*df)) then
          dyHold=dyH
          dyH=0.5*(yHl-yHh)
          yH=yHh+dyH
          if (yHh==yH) return
        else
          dyHold=dyH
          dyH=f/df
          temp=yH
          yH=yH-dyH
          if (temp==yH) return
        endif
        if (abs(dyH)<yHacc*yH) return
        call saha(variables,variable1,variable2,yH,f,df)
        if (f<0) then
          yHh=yH
        else
          yHl=yH
        endif
      enddo
!
      if (present(rterror)) rterror=.true.
!
    endsubroutine rtsafe
!***********************************************************************
    subroutine saha(variables,variable1,variable2,yH,f,df)
!
!   we want to find the root of f
!
!   23-feb-03/tobi: errors fixed
!
      use Mpicomm, only: stop_it
!
      character(len=*), intent(in) :: variables
      real, intent(in)             :: variable1,variable2,yH
      real, intent(out)            :: f,df
!
      real :: lnrho,ss,ee,pp
      real :: lnTT_,dlnTT_
!
      select case (variables)
      case ('lnrho|ss')
        lnrho=variable1
        ss=variable2
        lnTT_=(2./3.)*((ss/ss_ion+(1.-yH)*(log(1.-yH)-lnrho_H) &
                        +yH*(2.*log(yH)-lnrho_e-lnrho_p) &
                        +xHe_term)/(1.+yH+xHe) &
                       +lnrho-2.5)
      case ('lnrho|ee')
        lnrho=variable1
        ee=variable2
        lnTT_=log(ee/ee_ion-yH)-log(1.5*(1.+yH+xHe))
      case ('lnrho|pp')
        lnrho=variable1
        pp=variable2
        lnTT_=log(pp)-lnrho-log((1.+yH+xHe)*ee_ion)
      case default
        call stop_it("saha: I don't get what the independent variables are.")
      end select
!
      f=lnrho_e-lnrho+1.5*lnTT_-exp(-lnTT_)+log(1.-yH)-2.*log(yH)
      dlnTT_=((2./3.)*(lnrho_H-lnrho_p-f-exp(-lnTT_))-1)/(1.+yH+xHe)
      df=dlnTT_*(1.5+exp(-lnTT_))-1./(1.-yH)-2./yH
!
    endsubroutine saha
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
      real, dimension(mx) :: lnrho,yH,TT
      real :: kx,ky,kz
!
!  test
!
      if(radcalc_test) then
        if(lroot.and.ip<12) print*,'radcalc: using simple profiles for testing'
        !
        ! Periodic profiles
        !
        kx=2*pi/Lx
        ky=2*pi/Ly
        kz=2*pi/Lz
        Srad=1.+.02*spread(spread(cos(kx*x),2,my),3,mz) &
                   *spread(spread(cos(ky*y),1,mx),3,mz) &
                   *spread(spread(cos(kz*z),1,mx),2,my)
        kaprho=2.+spread(spread(cos(2*kx*x),2,my),3,mz) &
                 *spread(spread(cos(2*ky*y),1,mx),3,mz) &
                 *spread(spread(cos(2*kz*z),1,mx),2,my)
        !
        ! kaprho=const and S = tau^n
        !
!        kaprho = 1.e2
!        Srad = spread(spread((z*kaprho(1,1,1))**3,1,mx),2,my)
!
        return
      endif
!
!  no test
!
      do n=1,mz
      do m=1,my
!
         lnrho=f(:,m,n,ilnrho)
         yH=f(:,m,n,iyH)
         TT=f(:,m,n,iTT)
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
      real, dimension(mx,my,radz0) :: yH_xy,TT_xy
!
      if (nrad>0) then
        yH_xy=f(:,:,n1-radz0:n1-1,iyH)
        TT_xy=f(:,:,n1-radz0:n1-1,iTT)
      endif
!
      if (nrad<0) then
        yH_xy=f(:,:,n2+1:n2+radz0,iyH)
        TT_xy=f(:,:,n2+1:n2+radz0,iTT)
      endif
!
      H_xy=(1.+yH_xy+xHe)*ss_ion*TT_xy/gravz
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
      real, dimension(nx) :: lnrho,ss,yH,K,sqrtK,yH_term,one_yH_term
!
      do n=n1,n2
      do m=m1,m2
!
        lnrho=f(l1:l2,m,n,ilnrho)
!
        K=exp(lnrho_e-lnrho-TT_ion/T0)*(T0/TT_ion)**1.5
        sqrtK=sqrt(K)
        yH=2*sqrtK/(sqrtK+sqrt(4+K))
!
        where (yH>0)
          yH_term=yH*(2*log(yH)-lnrho_e-lnrho_p)
        elsewhere
          yH_term=0
        endwhere
!
        where (yH<1)
          one_yH_term=(1-yH)*(log(1-yH)-lnrho_H)
        elsewhere
          one_yH_term=0
        endwhere
!
        ss=ss_ion*((1+yH+xHe)*(1.5*log(T0/TT_ion)-lnrho+2.5) &
                   -yH_term-one_yH_term-xHe_term)
!
        f(l1:l2,m,n,iss)=ss
!
      enddo
      enddo
!
    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine bc_ss_flux(f,topbot,hcond0,hcond1,Fheat,FheatK,chi, &
                lmultilayer,lcalc_heatcond_constchi)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      real, intent(in) :: Fheat, FheatK, hcond0, hcond1, chi
      logical, intent(in) :: lmultilayer, lcalc_heatcond_constchi
      
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my) :: tmp_xy,TT_xy,rho_xy,yH_xy
      integer :: i
      
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  bottom boundary
!  ===============
!
      case('bot')
        if (lmultilayer) then
          if(headtt) print*,'bc_ss_flux: Fbot,hcond=',Fheat,hcond0*hcond1
        else
          if(headtt) print*,'bc_ss_flux: Fbot,hcond=',Fheat,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        TT_xy=f(:,:,n1,iTT)
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy. 
!
        if(lcalc_heatcond_constchi) then
          tmp_xy=Fheat/(rho_xy*chi*TT_xy)
        else
          tmp_xy=FheatK/TT_xy
        endif
!
!  get ionization fraction at bottom boundary
!
        yH_xy=f(:,:,n1,iyH)
!
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+ss_ion*(1+yH_xy+xHe)* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+3*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case('top')
        if (lmultilayer) then
          if(headtt) print*,'bc_ss_flux: Ftop,hcond=',Fheat,hcond0*hcond1
        else
          if(headtt) print*,'bc_ss_flux: Ftop,hcond=',Fheat,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n2,ilnrho))
        TT_xy=f(:,:,n2,iTT)
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy. 
!
        if(lcalc_heatcond_constchi) then
          tmp_xy=Fheat/(rho_xy*chi*TT_xy)
        else
          tmp_xy=FheatK/TT_xy
        endif
!
!  get ionization fraction at top boundary
!
        yH_xy=f(:,:,n2,iyH)
!
!  enforce ds/dz + gamma1/gamma*dlnrho/dz = - gamma1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+ss_ion*(1+yH_xy+xHe)* &
              (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)-3*i*dz*tmp_xy)
        enddo
!
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
!  26-aug-2003/tony: distributed across ionization modules
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
!
!  top boundary
!
      case('top')
        if ((bcz1(ilnrho) /= "a2") .and. (bcz1(ilnrho) /= "a3")) &
          call stop_it("bc_ss_temp_old: Inconsistent boundary conditions 3.")
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
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
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
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
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
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
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
!
!  top boundary
!
      case('top')
        if (ldebug) print*, &
                     'bc_ss_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp_z: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,:,n2,iss) = 0.5*tmp - gamma1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
               - gamma1/gamma*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
        enddo
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
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
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
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
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
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use Cdata
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,nghost) :: lnrho,ss,yH,TT,K,sqrtK,yH_term,one_yH_term
      integer :: i
!
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
        do i=1,nghost
          f(:,:,n1-i,iTT) = f(:,:,n1+i,iTT) 
        enddo
!
        lnrho=f(:,:,1:n1-1,ilnrho)
        TT=f(:,:,1:n1-1,iTT)
!
        K=exp(lnrho_e-lnrho-TT_ion/TT)*(TT/TT_ion)**1.5
        sqrtK=sqrt(K)
        yH=2*sqrtK/(sqrtK+sqrt(4+K))
!
        where (yH>0)
          yH_term=yH*(2*log(yH)-lnrho_e-lnrho_p)
        elsewhere
          yH_term=0
        endwhere
!
        where (yH<1)
          one_yH_term=(1-yH)*(log(1-yH)-lnrho_H)
        elsewhere
          one_yH_term=0
        endwhere
!
        ss=ss_ion*((1+yH+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5) &
                    -yH_term-one_yH_term-xHe_term)
!
        f(:,:,1:n1-1,iyH)=yH
        f(:,:,1:n1-1,iss)=ss
!
!  top boundary
!
      case('top')
        do i=1,nghost
          f(:,:,n2+i,iTT) = f(:,:,n2-i,iTT) 
        enddo
!
        lnrho=f(:,:,n2+1:mz,ilnrho)
        TT=f(:,:,n2+1:mz,iTT)
!
        K=exp(lnrho_e-lnrho-TT_ion/TT)*(TT/TT_ion)**1.5
        sqrtK=sqrt(K)
        yH=2*sqrtK/(sqrtK+sqrt(4+K))
!
        where (yH>0)
          yH_term=yH*(2*log(yH)-lnrho_e-lnrho_p)
        elsewhere
          yH_term=0
        endwhere
!
        where (yH<1)
          one_yH_term=(1-yH)*(log(1-yH)-lnrho_H)
        elsewhere
          one_yH_term=0
        endwhere
!
        ss=ss_ion*((1+yH+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5) &
                    -yH_term-one_yH_term-xHe_term)
!
        f(:,:,n2+1:mz,iyH)=yH
        f(:,:,n2+1:mz,iss)=ss
!
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
!  26-aug-2003/tony: distributed across ionization modules
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
