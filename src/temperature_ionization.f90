! $Id: temperature_ionization.f90,v 1.2 2006-04-02 18:29:53 theine Exp $

!  This module takes care of entropy (initial condition
!  and time advance)

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .false.
! CPARAM logical, parameter :: ltemperature = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED glnTT,uglnTT,TT1,yH,mu,pp,cv1,cp1,ee,ss
! PENCILS PROVIDED dlnppdlnTT,dlnppdlnrho,glnpp,cs2,Ma2
!
!***************************************************************
module Entropy

  use Cparam
  use Cdata
  use Messages

  implicit none

  include 'entropy.h'

  real :: radius_lnTT=0.1,ampl_lnTT=0.,widthlnTT=2*epsi
  real :: lnTT_left=1.0,lnTT_right=1.0,lnTT_const=0.0,TT_const=1.0
  real :: kx_lnTT=1.0,ky_lnTT=1.0,kz_lnTT=1.0
  real :: TT_ion,Rgas,rho_H,rho_e,rho_He
  real :: xHe=0.0,mu_0=0.0
  real :: heat_uniform=0.0,Kconst=0.0
  logical :: lpressuregradient_gas=.true.,ladvection_temperature=.true.
  logical :: lupw_lnTT,lcalc_heat_cool=.false.,lheatc_simple=.false.
  character (len=labellen), dimension(ninit) :: initlnTT='nothing'
  character (len=4) :: iinit_str

  ! Delete (or use) me asap!
  real :: hcond0,hcond1,Fbot,FbotKbot,Ftop,Kbot,FtopKtop,chi
  logical :: lmultilayer,lheatc_chiconst

  ! input parameters
  namelist /entropy_init_pars/ &
      initlnTT,radius_lnTT,ampl_lnTT,widthlnTT, &
      lnTT_left,lnTT_right,lnTT_const,TT_const, &
      kx_lnTT,ky_lnTT,kz_lnTT,xHe

  ! run parameters
  namelist /entropy_run_pars/ &
      lupw_lnTT,lpressuregradient_gas,ladvection_temperature, &
      xHe,heat_uniform,Kconst

  ! other variables (needs to be consistent with reset list below)
    integer :: idiag_TTmax=0,idiag_TTmin=0,idiag_TTm=0
    integer :: idiag_yHmax=0,idiag_yHmin=0,idiag_yHm=0
    integer :: idiag_eth=0,idiag_ssm=0,idiag_cv=0,idiag_cp=0
    integer :: idiag_dtchi=0,idiag_dtc

  contains

!***********************************************************************
    subroutine register_entropy()
!
!  initialise variables which should know that we solve an entropy
!  equation: ilnTT, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
      logical, save :: first=.true.
!
      if (.not. first) call fatal_error('register_entropy','module registration called twice')
      first = .false.
!
      ilnTT = nvar+1             ! index to access temperature
      nvar = nvar+1
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_entropy: nvar = ', nvar
        print*, 'register_entropy: ilnTT = ', ilnTT
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: temperature_ionization.f90,v 1.2 2006-04-02 18:29:53 theine Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call fatal_error('register_entropy','nvar > mvar')
      endif
!
!  Put variable name in array
!
      varname(ilnTT) = 'lnTT'
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',lnTT $'
          if (nvar == mvar) write(4,*) ',lnTT'
        else
          write(4,*) ',lnTT $'
        endif
        write(15,*) 'lnTT = fltarr(mx,my,mz)*one'
      endif
!
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
!
!  Check any module dependencies
!
      if (.not.leos_temperature_ionization) then
        call fatal_error('initialize_entropy','EOS/=noeos but'//&
                         'temperature_ionization already include'//&
                         'an EQUATION OF STATE for the fluid')
      endif
!
!  Useful constants for ionization
!
      TT_ion = chiH/k_B
      Rgas = k_B/m_H
      mu_0 = (1 + 4*xHe)
      rho_H = mu_0*m_H*((m_H/hbar)*(chiH/hbar)/(2*pi))**(1.5)
      rho_e = mu_0*m_H*((m_e/hbar)*(chiH/hbar)/(2*pi))**(1.5)
      rho_He = mu_0*m_H*((4*m_H/hbar)*(chiH/hbar)/(2*pi))**(1.5)
!
!  Turn off pressure gradient term and advection for 0-D runs.
!
      if (nxgrid*nygrid*nzgrid==1) then
        lpressuregradient_gas=.false.
        ladvection_temperature=.false.
        print*, 'initialize_entropy: 0-D run, turned off pressure gradient term'
        print*, 'initialize_entropy: 0-D run, turned off advection of entropy'
      endif
!
!  Check whether we want heating/cooling
!
      lcalc_heat_cool = (heat_uniform/=0.0)
!
!  Check whether we want heat conduction
!
      lheatc_simple = (Kconst/=0.0)

      if (NO_WARN) print*,f,lstarting  !(to keep compiler quiet)
!        
      endsubroutine initialize_entropy
!***********************************************************************
    subroutine read_entropy_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=entropy_init_pars,ERR=99, IOSTAT=iostat)
      else 
        read(unit,NML=entropy_init_pars,ERR=99) 
      endif

99    return
    endsubroutine read_entropy_init_pars
!***********************************************************************
    subroutine write_entropy_init_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_init_pars)
    endsubroutine write_entropy_init_pars
!***********************************************************************
    subroutine read_entropy_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=entropy_run_pars,ERR=99, IOSTAT=iostat)
      else 
        read(unit,NML=entropy_run_pars,ERR=99) 
      endif

99    return
    endsubroutine read_entropy_run_pars
!***********************************************************************
    subroutine write_entropy_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=entropy_run_pars)
    endsubroutine write_entropy_run_pars
!!***********************************************************************
    subroutine init_ss(f,xx,yy,zz)
!
!  initialise entropy; called from start.f90
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name]) 
!
      use General, only: chn
      use Sub, only: blob
      use Initcond, only: jump
!
      real, dimension (mx,my,mz,mvar+maux), intent (inout) :: f
      real, dimension (mx,my,mz), intent (in) :: xx,yy,zz
      logical :: lnothing=.true.
!
      do iinit=1,ninit
!
      if (initlnTT(iinit)/='nothing') then
!
      lnothing=.false.

      call chn(iinit,iinit_str)
!
!  select different initial conditions
!
      select case(initlnTT(iinit))

        case('zero', '0'); f(:,:,:,ilnTT) = 0.
        case('const_lnTT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+lnTT_const
        case('const_TT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+log(TT_const)
        case('blob'); call blob(ampl_lnTT,f,ilnTT,radius_lnTT,0.,0.,0.)
        case('xwave'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+ampl_lnTT*sin(kx_lnTT*xx)
        case('ywave'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+ampl_lnTT*sin(ky_lnTT*yy)
        case('zwave'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+ampl_lnTT*sin(kz_lnTT*zz)
        case('xjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'x')
        case('yjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'y')
        case('zjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'z')

        case default
          !
          !  Catch unknown values
          !
          write(unit=errormsg,fmt=*) 'No such value for initss(' &
                           //trim(iinit_str)//'): ',trim(initlnTT(iinit))
          call fatal_error('init_ss',errormsg)

      endselect

      if (lroot) print*,'init_ss: initss(' &
                        //trim(iinit_str)//') = ',trim(initlnTT(iinit))

      endif

      enddo

      if (lnothing.and.lroot) print*,'init_ss: nothing'
!
      if (NO_WARN) print*,xx,yy  !(to keep compiler quiet)
!
    endsubroutine init_ss
!***********************************************************************
    subroutine pencil_criteria_entropy()
! 
!  All pencils that the Entropy module depends on are specified here.
! 
!  20-11-04/anders: coded
!

      if (ldt) lpenc_requested(i_cs2)=.true.

      if (ldensity) then
        lpenc_requested(i_pp)=.true.
        lpenc_requested(i_dlnppdlnTT)=.true.
      endif

      if (lpressuregradient_gas) then
        lpenc_requested(i_pp)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnpp)=.true.
      endif

      if (lviscosity) then
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_TT1)=.true.
      endif

      if (lcalc_heat_cool) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_TT1)=.true.
      endif

      if (ladvection_temperature) lpenc_requested(i_uglnTT)=.true.

      if (lheatc_simple) lpenc_requested(i_cp1)=.true.

!
!  Diagnostics
!
      if (idiag_TTmax/=0) lpenc_diagnos(i_TT1)=.true.
      if (idiag_TTmin/=0) lpenc_diagnos(i_TT1)=.true.
      if (idiag_TTm/=0) lpenc_diagnos(i_TT1)=.true.
      if (idiag_yHmax/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_yHmin/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_yHm/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_eth/=0) then
        lpenc_diagnos(i_rho1)=.true.
        lpenc_diagnos(i_ee)=.true.
      endif
      if (idiag_ssm/=0) lpenc_diagnos(i_ss)=.true.
      if (idiag_cv/=0) lpenc_diagnos(i_cv1)=.true.
      if (idiag_cp/=0) lpenc_diagnos(i_cp1)=.true.
      if (idiag_dtchi/=0) then
        lpenc_diagnos(i_rho1)=.true.
        lpenc_diagnos(i_cv1)=.true.
      endif
      if (idiag_dtchi/=0) lpenc_diagnos(i_cs2)=.true.

    endsubroutine pencil_criteria_entropy
!***********************************************************************
    subroutine pencil_interdep_entropy(lpencil_in)
!       
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_Ma2)) then
        lpencil_in(i_u2)=.true.
        lpencil_in(i_cs2)=.true.
      endif

      if (lpencil_in(i_cs2)) then
        lpencil_in(i_pp)=.true.
        lpencil_in(i_rho1)=.true.
        lpencil_in(i_cv1)=.true.
        lpencil_in(i_TT1)=.true.
        lpencil_in(i_dlnppdlnrho)=.true.
        lpencil_in(i_dlnppdlnTT)=.true.
      endif

      if (lpencil_in(i_glnpp)) then
        lpencil_in(i_dlnppdlnrho)=.true.
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_dlnppdlnTT)=.true.
        lpencil_in(i_glnTT)=.true.
      endif

      if (lpencil_in(i_cv1)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_TT1)=.true.
        lpencil_in(i_mu)=.true.
      endif

      if (lpencil_in(i_pp)) then
        lpencil_in(i_rho1)=.true.
        lpencil_in(i_TT1)=.true.
        lpencil_in(i_mu)=.true.
      endif

      if (lpencil_in(i_dlnppdlnTT)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_TT1)=.true.
      endif

      if (lpencil_in(i_yH)) then
        lpencil_in(i_rho1)=.true.
        lpencil_in(i_TT1)=.true.
      endif

      if (lpencil_in(i_mu)) lpencil_in(i_yH)=.true.

      if (lpencil_in(i_dlnppdlnrho)) lpencil_in(i_yH)=.true.

      if (lpencil_in(i_uglnTT)) lpencil_in(i_glnTT)=.true.

      if (lpencil_in(i_ee)) then
        lpencil_in(i_mu)=.true.
        lpencil_in(i_TT1)=.true.
        lpencil_in(i_yH)=.true.
      endif

      if (lpencil_in(i_ss)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_rho1)=.true.
        lpencil_in(i_TT1)=.true.
      endif

    endsubroutine pencil_interdep_entropy
!***********************************************************************
    subroutine calc_pencils_entropy(f,p)
!       
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use Sub, only: grad,u_dot_gradf

      real, dimension (mx,my,mz,mvar+maux), intent (in) :: f
      type (pencil_case), intent (inout) :: p

      real, dimension (nx) :: tmp,tmp1,tmp2,rhs
      integer :: i

!
!  Temperature gradient
!
      if (lpencil(i_glnTT)) call grad(f,ilnTT,p%glnTT)

!
!  Temperature advection
!
      if (lpencil(i_uglnTT)) then
        call u_dot_gradf(f,ilnTT,p%glnTT,p%uu,p%uglnTT,UPWIND=lupw_lnTT)
      endif

!
!  Inverse temperature
!
      if (lpencil(i_TT1)) then
        p%TT1=exp(-f(l1:l2,m,n,ilnTT))
      endif

!
!  Ionization fraction
!
      if (lpencil(i_yH)) then
        where (TT_ion*p%TT1 < -log(tiny(TT_ion)))
          rhs = rho_e*p%rho1*(p%TT1*TT_ion)**(-1.5)*exp(-TT_ion*p%TT1)
        elsewhere
          rhs = 0
        endwhere
        p%yH = 2*sqrt(rhs)/(sqrt(rhs)+sqrt(4+rhs))
      endif

!
!  Mean molecular weight
!
      if (lpencil(i_mu)) p%mu = mu_0/(1 + p%yH + xHe)

!
!  Pressure
!
      if (lpencil(i_pp)) p%pp = Rgas/(p%mu*p%rho1*p%TT1)

!
!  Inverse specific heat at constant volume (i.e. density)
!
      if (lpencil(i_cv1)) then
        tmp1 = p%yH*(1-p%yH)/((2-p%yH)*(1+p%yH+xHe))
        tmp2 = 1.5 + TT_ion*p%TT1
        p%cv1 = (p%mu/Rgas)/(1.5 + tmp1*tmp2**2)
      endif

!
!  Inverse specific heat at constant pressure
!
      if (lpencil(i_cp1)) then
        tmp1 = p%yH*(1-p%yH)/(2 + xHe*(2-p%yH))
        tmp2 = 2.5+TT_ion*p%TT1
        p%cp1 = (p%mu/Rgas)/(2.5 + tmp1*tmp2**2)
      endif

!
!  Coefficients for the pressure gradient
!
      if (lpencil(i_dlnppdlnTT)) then
        tmp1 = p%yH * (1-p%yH) / ( (2-p%yH) * (1+p%yH+xHe) )
        tmp2 = 1.5+TT_ion*p%TT1
        p%dlnppdlnTT = 1 + tmp1*tmp2
      endif
      if (lpencil(i_dlnppdlnrho)) then
        tmp1 = p%yH * (1-p%yH) / ( (2-p%yH) * (1+p%yH+xHe) )
        p%dlnppdlnrho = 1 - tmp1
      endif

!
!  Logarithmic pressure gradient
!
      if (lpencil(i_glnpp)) then
        do i=1,3
          p%glnpp(:,i) = p%dlnppdlnrho*p%glnrho(:,i) + p%dlnppdlnTT*p%glnTT(:,i)
        enddo
      endif

!
!  Sound speed
!
      if (lpencil(i_cs2)) then
        p%cs2 = p%pp*p%rho1*(p%pp*p%rho1*p%cv1*p%TT1*(p%dlnppdlnTT)**2 &
                              + p%dlnppdlnrho)
      endif

!
!  Mach Speed
!
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
!
!  Energy per unit mass
!
      if (lpencil(i_ee)) p%ee = 1.5*Rgas/(p%mu*p%TT1) + p%yH*(Rgas/mu_0)*TT_ion
!
!  Entropy per unit mass
!  The contributions from each particle species contain the mixing entropy
!
      if (lpencil(i_ss)) then
        tmp = 2.5 - 1.5*log(TT_ion*p%TT1)
        ! Hydrogen
        p%ss = tmp + log(rho_H*p%rho1)
        ! Electrons
        where (p%yH > 0)
          p%ss = p%ss + p%yH*(tmp + log(rho_e*p%rho1) - log(p%yH))
        endwhere
        ! Helium
        if (xHe > 0) then
          p%ss = p%ss + xHe*(tmp + log(rho_He*p%rho1) - log(xHe))
        endif
        ! Overall factor
        p%ss = (Rgas/mu_0)*p%ss
      endif

    endsubroutine calc_pencils_entropy
!**********************************************************************
    subroutine dss_dt(f,df,p)
!
!  calculate right hand side of entropy equation
!  heat condution is currently disabled until old stuff,
!  which in now in calc_heatcond, has been reinstalled.
!
!  17-sep-01/axel: coded
!   9-jun-02/axel: pressure gradient added to du/dt already here
!   2-feb-03/axel: added possibility of ionization
!
      use Viscosity, only: calc_viscous_heat
      use Sub, only: max_mn_name,sum_mn_name,identify_bcs
!
      real, dimension (mx,my,mz,mvar+maux), intent (inout) :: f
      real, dimension (mx,my,mz,mvar), intent (out) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: visc_heat,Hmax
      integer :: j
!
!  Identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dss_dt: SOLVE dss_dt'
      if (headtt) call identify_bcs('lnTT',ilnTT)
!
!  Calculate cs2 in a separate routine
!
      if (headtt) print*,'dss_dt: cs2 =', p%cs2(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
      if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!  Pressure term in momentum equation (setting lpressuregradient_gas to
!  .false. allows suppressing pressure term for test purposes)
!
      if (lhydro.and.lpressuregradient_gas) then
        do j=0,2
          df(l1:l2,m,n,iuu+j) = df(l1:l2,m,n,iuu+j) - p%pp*p%rho1*p%glnpp(:,j+1)
        enddo
      endif
!
!  Advection term
!
      if (ladvection_temperature) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%uglnTT
      endif
!
!  Calculate viscous contribution to temperature
!  (visc_heat has units of energy/mass)
!
      if (lviscosity) then
        call calc_viscous_heat(df,p,visc_heat,Hmax)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + p%cv1*p%TT1*visc_heat
      endif
!
!  Various heating/cooling mechanisms
!
      if (lcalc_heat_cool) call calc_heat_cool(f,df,p)
!
!  Thermal conduction
!
      if (lheatc_simple) call calc_heatcond_simple(f,df,p)
!
!  Need to add left-hand-side of the continuity equation (see manual)
!
      if (ldensity) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + &
          p%pp*p%rho1*p%cv1*p%TT1*p%dlnppdlnTT*df(l1:l2,m,n,ilnrho)
      endif
!
!  Calculate temperature related diagnostics
!
      if (ldiagnos) then
        if (idiag_TTmax/=0) call max_mn_name(1/p%TT1,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-1/p%TT1,idiag_TTmin,lneg=.true.)
        if (idiag_TTm/=0) call sum_mn_name(1/p%TT1,idiag_TTm)
        if (idiag_yHmax/=0) call max_mn_name(p%yH,idiag_yHmax)
        if (idiag_yHmin/=0) call max_mn_name(-p%yH,idiag_yHmin,lneg=.true.)
        if (idiag_yHm/=0) call sum_mn_name(p%yH,idiag_yHm)
        if (idiag_eth/=0) call sum_mn_name(p%ee/p%rho1,idiag_eth)
        if (idiag_ssm/=0) call sum_mn_name(p%ss,idiag_ssm)
        if (idiag_cv/=0) call sum_mn_name(1/p%cv1,idiag_cv)
        if (idiag_cp/=0) call sum_mn_name(1/p%cp1,idiag_cp)
        if (idiag_dtc/=0) then
          call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        endif
      endif

    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_heatcond_simple(f,df,p)

      use Sub, only: max_mn_name,dot,del2

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx) :: chi,gamma,glnTT2,del2lnTT
!
!  Thermal diffusivity
!
      chi = Kconst*p%rho1*p%cp1
!
!  ``gamma''
!
      gamma = p%cv1/p%cp1
!
!  glnTT2 and del2lnTT
!
      call dot(p%glnTT,p%glnTT,glnTT2)
      call del2(f,iss,del2lnTT)
!
!  Add heat conduction to RHS of temperature equation
!
      df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + gamma*chi*(glnTT2 + del2lnTT)
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) then
        diffus_chi=max(diffus_chi,chi*dxyz_2)
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif

    end subroutine calc_heatcond_simple
!***********************************************************************
    subroutine calc_heat_cool(f,df,p)

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx) :: heat
!
!  Initialize
!
      heat=0
!
!  Add spatially uniform heating (usually as a test)
!
      if (heat_uniform/=0.0) heat = heat+heat_uniform
!
!  Add to RHS of temperature equation
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + p%rho1*p%cv1*p%TT1*heat

    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Sub, only: parse_name
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0
        idiag_yHmax=0; idiag_yHmin=0; idiag_yHm=0
        idiag_eth=0; idiag_ssm=0; idiag_cv=0; idiag_cp=0
        idiag_dtchi=0; idiag_dtc=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'yHmax',idiag_yHmax)
        call parse_name(iname,cname(iname),cform(iname),'yHmin',idiag_yHmin)
        call parse_name(iname,cname(iname),cform(iname),'yHm',idiag_yHm)
        call parse_name(iname,cname(iname),cform(iname),'eth',idiag_eth)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'cv',idiag_cv)
        call parse_name(iname,cname(iname),cform(iname),'cp',idiag_cp)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'ilnTT=',ilnTT
        write(3,*) 'iss=',iss
        write(3,*) 'iyH=',iyH
        write(3,*) 'i_TTmax=',idiag_TTmax
        write(3,*) 'i_TTmin=',idiag_TTmin
        write(3,*) 'i_TTm=',idiag_TTm
        write(3,*) 'i_yHmax=',idiag_yHmax
        write(3,*) 'i_yHmin=',idiag_yHmin
        write(3,*) 'i_yHm=',idiag_yHm
        write(3,*) 'i_eth=',idiag_eth
        write(3,*) 'i_ssm=',idiag_ssm
        write(3,*) 'i_cv=',idiag_cv
        write(3,*) 'i_cp=',idiag_cp
        write(3,*) 'i_dtchi=',idiag_dtchi
        write(3,*) 'i_dtc=',idiag_dtc
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************    
endmodule Entropy
