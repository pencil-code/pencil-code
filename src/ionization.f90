! $Id: ionization.f90,v 1.57 2003-06-30 09:33:57 brandenb Exp $

!  This modules contains the routines for simulation with
!  simple hydrogen ionization.

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

  !  secondary parameters calculated in initialize
  real :: m_H,m_He,mu,twothirds
  real :: TT_ion,TT_ion_,chiH,chiH_,ss_ion,kappa0
  real :: lnrho_H,lnrho_e,lnrho_e_,lnrho_p,lnrho_He
  integer :: l0,l3,m0,m3,n0,n3

  !  lionization initialized to .true.
  !  it can be reset to .false. in namelist
  logical :: lionization=.true.,lfixed_ionization=.false.,output_yH=.false.
  character (len=labellen) :: cionization='hydrogen'
  real :: yH0=impossible,xHe=0.1,lnxHe

  ! input parameters
  namelist /ionization_init_pars/ cionization,output_yH

  ! run parameters
  namelist /ionization_run_pars/  cionization,output_yH

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
           "$Id: ionization.f90,v 1.57 2003-06-30 09:33:57 brandenb Exp $")
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
!  ionization parameters
!  since m_e and chiH, as well as hbar are all very small
!  it is better to divide m_e and chiH separately by hbar.
!
      m_H=m_p+m_e
      m_He=3.97153*m_H
      lnxHe=log(xHe)
      twothirds=2./3.
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
      do n=n0,n3
      do m=m0,m3
      do l=l0,l3
         lnrho=f(l,m,n,ilnrho)
         ss=f(l,m,n,ient)
         yH=f(l,m,n,iyH)
         call rtsafe(lnrho,ss,yH)
         f(l,m,n,iyH)=yH
         lnTT_=twothirds*((ss/ss_ion+(1.-yH)*(log(1.-yH)-lnrho_H) &
                           +yH*(2.*log(yH)-lnrho_e-lnrho_p) &
                           +xHe*(lnxHe-lnrho_He))/(1.+yH+xHe) &
                          +lnrho-2.5)
         f(l,m,n,iTT)=exp(lnTT_)*TT_ion
      enddo
      enddo
      enddo
!
    endsubroutine ioncalc
!***********************************************************************
    subroutine ioncalc_penc(lnrho,ss,TT,yH)
!
!   calculate degree of ionization and temperature
!   This routine is called from equ.f90 and operates on the full 3-D array.
!
!   19-jun-03/tobi: tony
!
      use Cdata
!
      real, dimension(nx), intent(in)  :: lnrho,ss
      real, dimension(nx), intent(inout) :: yH,TT
      real :: yHsave,TTcalc,lnTT_
      integer :: l
!
!  do the loop
!
      do l=1,nx
         yHsave=yH(l)
         call rtsafe(lnrho(l),ss(l),yHsave)
         yH(l)=yHsave
         lnTT_=twothirds*((ss(l)/ss_ion &
                           +(1.-yH(l))*(log(1.-yH(l))-lnrho_H) &
                           +yH(l)*(2.*log(yH(l))-lnrho_e-lnrho_p) &
                           +xHe*(lnxHe-lnrho_He))/(1.+yH(l)+xHe) &
                          +lnrho(l)-2.5)
         TT(l)=exp(lnTT_)*TT_ion
      enddo
!
    endsubroutine ioncalc_penc
!***********************************************************************
    subroutine ioncalc_point(lnrho,ss,TT,yH)
!
!   calculate degree of ionization and temperature
!   This routine is called from equ.f90 and operates on the full 3-D array.
!
!   19-jun-03/tobi: tony
!
      use Cdata
!
      real, intent(in)  :: lnrho,ss
      real, intent(inout) :: yH,TT
      real :: lnTT_
      integer :: l
!
!  do the loop
!
      call rtsafe(lnrho,ss,yH)
      lnTT_=twothirds*((ss/ss_ion+(1.-yH)*(log(1.-yH)-lnrho_H) &
                        +yH*(2.*log(yH)-lnrho_e-lnrho_p) &
                        +xHe*(lnxHe-lnrho_He))/(1.+yH+xHe) &
                       +lnrho-2.5)
      TT=exp(lnTT_)*TT_ion
!
    endsubroutine ioncalc_point
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
!  set data on pencil only
!
      yH=f(l1:l2,m,n,iyH)
      TT=f(l1:l2,m,n,iTT)
!
      if(ip==0) print*,ss,lnrho  !(keep compiler quiet)
    endsubroutine ionset

!***********************************************************************
    subroutine output_ionization(lun)
!
!  Optional output of derived quantities along with VAR-file
!  Called from wsnap
!
!   5-apr-03/axel: coded
!
      use Cdata
!
      integer, intent(in) :: lun
!
!  identifier
!
      !if(headtt.and.ip<10) print*,'output_ionization',yyH(4,4,4)
      !if(output_yH) write(lun) yyH
!
    endsubroutine output_ionization

!***********************************************************************
    subroutine thermodynamics_penc(f,TT1,cs2,cp1tilde,ee,yHout)
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
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (nx), intent(out), optional :: cs2,TT1,cp1tilde,ee,yHout
      real, dimension (nx) :: yH,TT,lnrho,ss
      real, dimension (nx) :: ff,dlnTT_dy,dffdy,dlnTT_dlnrho
      real, dimension (nx) :: dffdlnrho,dydlnrho,dlnPdlnrho
      real, dimension (nx) :: dlnTT_dss,dffdss,dydss,dlnPdss
!
! P...p...p...pick up a Pencil for readability
!
      lnrho=f(l1:l2,m,n,ilnrho)
      ss=f(l1:l2,m,n,ient)
      yH=f(l1:l2,m,n,iyH)
      TT=f(l1:l2,m,n,iTT)
!
!  calculate 1/T (=TT1), sound speed, and coefficient cp1tilde in
!  the expression (1/rho)*gradp = cs2*(gradlnrho + cp1tilde*gradss)
!
      if (present(TT1)) TT1=1./TT
      if (present(yHout)) yHout=yH
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
    endsubroutine thermodynamics_penc
!***********************************************************************
    subroutine thermodynamics_arbpenc(lnrho,ss,TT1,cs2,cp1tilde,ee,yHout)
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
      real, dimension (nx), intent(in) :: lnrho,ss
      real, dimension (nx), intent(out), optional :: cs2,TT1,cp1tilde,ee,yHout
      real, dimension (nx) :: yH,TT
      real, dimension (nx) :: ff,dlnTT_dy,dffdy,dlnTT_dlnrho
      real, dimension (nx) :: dffdlnrho,dydlnrho,dlnPdlnrho
      real, dimension (nx) :: dlnTT_dss,dffdss,dydss,dlnPdss

!  calculate TT, yH for arbitrary pencil
      call ioncalc_penc(lnrho,ss,TT,yH)
!
!  calculate 1/T (=TT1), sound speed, and coefficient cp1tilde in
!  the expression (1/rho)*gradp = cs2*(gradlnrho + cp1tilde*gradss)
!
      if (present(TT1)) TT1=1./TT
      if (present(yHout)) yHout=yH
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
    endsubroutine thermodynamics_arbpenc
!***********************************************************************
    subroutine thermodynamics_arbpoint(lnrho,ss,TT1,cs2,cp1tilde,ee,yHout)
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
      real, intent(in) :: lnrho,ss
      real, intent(out), optional :: cs2,TT1,cp1tilde,ee,yHout
      real :: yH,TT
      real :: ff,dlnTT_dy,dffdy,dlnTT_dlnrho
      real :: dffdlnrho,dydlnrho,dlnPdlnrho
      real :: dlnTT_dss,dffdss,dydss,dlnPdss
!
! P...p...p..pick up a Pencil for readability
!
      call ioncalc_point(lnrho,ss,TT,yH)
!
!  calculate 1/T (=TT1), sound speed, and coefficient cp1tilde in
!  the expression (1/rho)*gradp = cs2*(gradlnrho + cp1tilde*gradss)
!
      if (present(TT1)) TT1=1./TT
      if (present(yHout)) yHout=yH
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
    endsubroutine thermodynamics_arbpoint
!***********************************************************************
    subroutine opacity(f,kaprho)
!
!  calculate opacity
!
!  26-jun-03/tobi: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz) :: kaprho
!
!  The factor 2 in front of lnrho takes care of the extra rho factor in kaprho
!
      kaprho(l0:l3,m,n)=.25*exp(2.*f(l0:l3,m,n,ilnrho)-lnrho_e_) &
                        *(TT_ion_/f(l0:l3,m,n,iTT))**1.5 &
                        *exp(TT_ion_/f(l0:l3,m,n,iTT)) &
                        *f(l0:l3,m,n,iyH)*(1.-f(l0:l3,m,n,iyH))*kappa0
!
    endsubroutine opacity
!***********************************************************************
    subroutine sourcefunction(f,Srad)
!
!  calculate sourcefunction
!
!  26-jun-03/tobi: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz), intent(out) :: Srad
!
      Srad(l0:l3,m,n)=sigmaSB*f(l0:l3,m,n,iTT)**4/pi
!
    endsubroutine sourcefunction
!***********************************************************************
    subroutine rtsafe(lnrho,ss,yH)
!
!   safe newton raphson algorithm (adapted from NR) !
!   09-apr-03/tobi: changed to subroutine
!
      real, intent (in)    :: lnrho,ss
      real, intent (inout) :: yH
      real                 :: yH0,yH1,dyHold,dyH,fl,fh,f,df,temp
      real, parameter      :: yHacc=1e-5
      integer              :: i
      integer, parameter   :: maxit=100

      yH1=1.-yHacc
      yH0=yHacc
      dyHold=1.-2.*yHacc
      dyH=dyHold
!
!  return if y is too close to 0 or 1
!
      call saha(yH0,lnrho,ss,fh,df)
      if (fh.le.0.) then
         yH=yH0
         return
      endif
      call saha(yH1,lnrho,ss,fl,df)
      if (fl.ge.0.) then
         yH=yH1
         return
      endif
!
!  otherwise find root
!
      call saha(yH,lnrho,ss,f,df)
      do i=1,maxit
         if (((yH-yH0)*df-f)*((yH-yH1)*df-f).gt.0. &
              .or.abs(2.*f).gt.abs(dyHold*df)) then
            dyHold=dyH
            dyH=.5*(yH0-yH1)
            yH=yH1+dyH
            if (yH1.eq.yH) return
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
            yH1=yH
         else
            yH0=yH
         endif
      enddo
      print *,'rtsafe: exceeded maximum iterations',f
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
                        +xHe*(lnxHe-lnrho_He))/(1.+yH+xHe) &
                       +lnrho-2.5)
      f=lnrho_e-lnrho+1.5*lnTT_-exp(-lnTT_)+log(1.-yH)-2.*log(yH)
      dlnTT_=(twothirds*(lnrho_H-lnrho_p-f-exp(-lnTT_))-1)/(1.+yH+xHe)
      df=dlnTT_*(1.5+exp(-lnTT_))-1./(1.-yH)-2./yH
    endsubroutine saha

endmodule ionization
