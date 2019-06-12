!$Id$
!  This module can replace the entropy module by using _T_ as dependent
!  variable. For a perfect gas with constant coefficients (no ionization)
!  we have (1-1/gamma) * cp*T = cs02 * exp( (gamma-1)*ln(rho/rho0)-gamma*s/cp )
!
!  At a later point we may want to rename the module Entropy into Energy

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
! PENCILS PROVIDED Ma2; uglnTT; fpres(3)
!
!***************************************************************
module Entropy

  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use EquationOfState, only: mpoly0,mpoly1
  use Interstellar

  implicit none
  public :: ADI_constK, ADI_Kprof


  include 'entropy.h'

  real :: radius_lnTT=0.1,ampl_lnTT=0.,widthlnTT=2*epsi
  real :: lnTT_left=1.0,lnTT_right=1.0,lnTT_const=0.0,TT_const=1.0
  real :: kx_lnTT=1.0,ky_lnTT=1.0,kz_lnTT=1.0
  real :: chi=impossible,heat_uniform=0.0
  real :: zbot=0.0,ztop=0.0
  real :: tau_heat_cor=-1.0,tau_damp_cor=-1.0,zcor=0.0,TT_cor=0.0
  real :: center1_x=0., center1_y=0., center1_z=0.
  real :: r_bcz=0.
  real :: Tbump=0.,Kmin=0.,Kmax=0.
  integer, parameter :: nheatc_max=1
  logical :: lpressuregradient_gas=.true.,ladvection_temperature=.true.
  logical :: lupw_lnTT=.false.
  logical :: lheatc_Kconst=.false.,lheatc_Kprof=.false.,lheatc_Karctan=.false.
  logical :: lheatc_chiconst=.false.,lheatc_chiconst_accurate=.false.
  character (len=labellen) :: iheatcond='nothing'
  logical :: lhcond_global=.false.
  logical :: lviscosity_heat=.true.
  integer :: iglobal_hcond=0
  integer :: iglobal_glhc=0

  character (len=labellen), dimension(ninit) :: initlnTT='nothing'
  character (len=intlen) :: iinit_str

! Delete (or use) me asap!
  real :: hcond0=impossible, hcond1=1.,Fbot,FbotKbot,Ftop,Kbot,FtopKtop
  logical :: lmultilayer=.false.

! input parameters
  namelist /entropy_init_pars/ &
      initlnTT,radius_lnTT,ampl_lnTT,widthlnTT, &
      lnTT_left,lnTT_right,lnTT_const,TT_const, &
      kx_lnTT,ky_lnTT,kz_lnTT,center1_x,center1_y,center1_z, &
      mpoly0,mpoly1,r_bcz, &
      Fbot,Tbump,Kmin,Kmax

! run parameters
  namelist /entropy_run_pars/ &
      lupw_lnTT,lpressuregradient_gas,ladvection_temperature, &
      heat_uniform,chi,iheatcond,tau_heat_cor,tau_damp_cor,zcor,TT_cor, &
      lheatc_chiconst_accurate,hcond0, &
      Tbump,Kmin,Kmax, &
      widthlnTT,mpoly0,mpoly1, &
      lhcond_global,lviscosity_heat, &
      Fbot,Tbump,Kmin,Kmax
!
! other variables (needs to be consistent with reset list below)
  integer :: idiag_TTmax=0,idiag_TTmin=0,idiag_TTm=0
  integer :: idiag_ethm=0,idiag_ssm=0
  integer :: idiag_dtchi=0,idiag_dtc=0
  integer :: idiag_eem=0,idiag_ppm=0,idiag_csm=0
 
  contains

!***********************************************************************
    subroutine register_entropy()
!
!  initialise variables which should know that we solve an entropy
!  equation: ilnTT, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
      use FArrayManager
!
      call farray_register_pde('lnTT',ilnTT)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
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
    subroutine initialize_entropy(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
!
      use Cdata
      use FArrayManager
      use Gravity, only: g0
      use EquationOfState
      use Sub, only: step,der_step
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx)   :: hcond,dhcond
      logical :: lnothing
!
      if (.not. leos) then
         call fatal_error('initialize_entropy','EOS=noeos but temperature_TT requires an EQUATION OF STATE for the fluid')
      endif
!
      call select_eos_variable('TT',ilnTT)
!
!  Check whether we want heat conduction
!
      lheatc_Kconst= .false.
      lheatc_Kprof= .false.
      lheatc_Karctan= .false.
      lheatc_chiconst = .false.
      lnothing = .false.
!
      select case (iheatcond)
        case ('K-const')
          lheatc_Kconst=.true.
          call information('initialize_TT', &
          ' heat conduction: K=cst --> gamma*K/(rho*cp)*lapl(T)')
          if (initlnTT(1).ne.'rad_equil') &
            Fbot=gamma/(gamma-1.)*hcond0*g0/(mpoly0+1.)
        case ('K-profile')
          lheatc_Kprof=.true.
! 
!  TODO..... ailleurs !
!
          hcond1=(mpoly1+1.)/(mpoly0+1.)
          Fbot=gamma/(gamma-1.)*hcond0*g0/(mpoly0+1.)
          call information('initialize_entropy',' heat conduction: K=K(r)')
        case ('K-arctan')
          lheatc_Karctan=.true.         
          call information('initialize_entropy',' heat conduction: arctan profile')
        case ('chi-const')
          lheatc_chiconst=.true.
          call information('initialize_entropy',' heat conduction: constant chi')
        case ('nothing')
          if (lroot .and. (.not. lnothing)) print*,'heat conduction: nothing'
        case default
          if (lroot) then
            write(unit=errormsg,fmt=*)  &
                'No such value iheatcond = ', trim(iheatcond)
            call fatal_error('initialize_entropy',errormsg)
          endif
       endselect
       lnothing=.true.
!
!  compute and store hcond and dhcond if hcond_global=.true.
!
       if (lhcond_global) then
         call farray_register_global("hcond",iglobal_hcond)
         call farray_register_global("glhc",iglobal_glhc)
         do n=n1,n2
         do m=m1,m2
           hcond = 1. + (hcond1-1.)*step(x(l1:l2),r_bcz,-widthlnTT)
           hcond = hcond0*hcond
           dhcond = hcond0*(hcond1-1.)*der_step(x(l1:l2),r_bcz,-widthlnTT)
           f(l1:l2,m,n,iglobal_hcond)=hcond
           f(l1:l2,m,n,iglobal_glhc)=dhcond
         enddo
         enddo
       endif
!
!  A word of warning...
!
      if (lheatc_Kconst .and. hcond0==0.0) then
        call warning('initialize_entropy', 'hcond0 is zero!')
      endif
      if (lheatc_Kprof .and. hcond0==0.0) then
        call warning('initialize_entropy', 'hcond0 is zero!')
      endif
      if (lheatc_chiconst .and. chi==0.0) then
        call warning('initialize_entropy','chi is zero!')
      endif
      if (iheatcond=='nothing') then
        if (hcond0 /= impossible) call warning('initialize_entropy', 'No heat conduction, but hcond0 /= 0')
        if (chi /= impossible) call warning('initialize_entropy', 'No heat conduction, but chi /= 0')
      endif
!
    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f)
!
!  initialise lnTT; called from start.f90
!
!  13-dec-2002/axel+tobi: adapted from init_ss
!
!  initialise entropy; called from start.f90
!  07-nov-2001/wolf: coded 
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!
      use General,  only: itoa
      use Sub,      only: blob
      use Initcond, only: jump
      use InitialCondition, only: initial_condition_ss
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      logical :: lnothing=.true.
      integer :: j
!
      do j=1,ninit
!
        if (initlnTT(j)/='nothing') then
!
          lnothing=.false.
!
          iinit_str=itoa(j)
!
!  select different initial conditions
!
          select case (initlnTT(j))
          case ('zero', '0'); f(:,:,:,ilnTT) = 0.
          case ('const_lnTT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+lnTT_const
          case ('const_TT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+log(TT_const)
          case ('single_polytrope'); call single_polytrope(f)
          case ('rad_equil')
            call rad_equil(f)
            if (ampl_lnTT.ne.0.) then
              print*,'add a bubble with:',ampl_lnTT,radius_lnTT,center1_x,center1_y,center1_z
              call blob(ampl_lnTT,f,ilnTT,radius_lnTT,center1_x,center1_y,center1_z)
              call blob(-ampl_lnTT,f,ilnrho,radius_lnTT,center1_x,center1_y,center1_z)
            endif
          case ('bubble_hs')
!         print*,'init_lnTT: put bubble in hydrostatic equilibrium: radius_lnTT,ampl_lnTT=',radius_lnTT,ampl_lnTT,center1_x,center1_y,center1_z
            call blob(ampl_lnTT,f,ilnTT,radius_lnTT,center1_x,center1_y,center1_z)
            call blob(-ampl_lnTT,f,ilnrho,radius_lnTT,center1_x,center1_y,center1_z)
!
          case default
          !
          !  Catch unknown values
          !
            write(unit=errormsg,fmt=*) 'No such value for init_TT(' &
                //trim(iinit_str)//'): ',trim(initlnTT(j))
            call fatal_error('init_TT',errormsg)

          endselect

          if (lroot) print*,'init_TT: init_TT(' &
              //trim(iinit_str)//') = ',trim(initlnTT(j))
        endif
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_ss(f)
!
      if (lnothing.and.lroot) print*,'init_ss: nothing'
!
    endsubroutine init_ss
!***********************************************************************
    subroutine pencil_criteria_entropy()
!
      use Cdata
!
      if (ldt) lpenc_requested(i_cs2)=.true.
      if (lpressuregradient_gas) lpenc_requested(i_fpres)=.true.
!
      if (lviscosity.and.lviscosity_heat) then
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_visc_heat)=.true.
      endif
!
      if (ldensity) lpenc_requested(i_divu)=.true.
!
      if (ladvection_temperature) lpenc_requested(i_uglnTT)=.true.
!
      if (lheatc_chiconst) lpenc_requested(i_del2lnTT)=.true.
!
      if (lheatc_Kconst) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (lheatc_Kprof) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (lheatc_Karctan) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
!  Diagnostics
!
      if (idiag_TTmax/=0) lpenc_diagnos(i_TT)  =.true.
      if (idiag_TTmin/=0) lpenc_diagnos(i_TT)  =.true.
      if (idiag_TTm/=0)   lpenc_diagnos(i_TT)  =.true.
      if (idiag_ethm/=0) then
                          lpenc_diagnos(i_rho1)=.true.
                          lpenc_diagnos(i_ee)  =.true.
      endif
      if (idiag_ssm/=0)   lpenc_diagnos(i_ss)  =.true.
      if (idiag_dtchi/=0) then
                          lpenc_diagnos(i_rho1)=.true.
                          lpenc_diagnos(i_cv1) =.true.
      endif
      if (idiag_dtchi/=0)  lpenc_diagnos(i_cs2)=.true.
      if (idiag_csm/=0)    lpenc_diagnos(i_cs2)=.true.
      if (idiag_eem/=0)    lpenc_diagnos(i_ee) =.true.
      if (idiag_ppm/=0)    lpenc_diagnos(i_pp) =.true.
!
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
!
      if (lpencil_in(i_uglnTT)) lpencil_in(i_glnTT)=.true.
!
      if (lpencil_in(i_fpres)) then
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_glnTT)=.true.
        lpencil_in(i_TT)=.true.
      endif
!
    endsubroutine pencil_interdep_entropy
!*********************************************************************** 
    subroutine calc_pencils_entropy(f,p)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
! 
      use EquationOfState
      use Sub

      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      type (pencil_case), intent (inout) :: p
      integer :: j
!
!  Mach Speed
!
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
!
!  Temperature advection
!  (Needs to be here because of lupw_lnTT)
!
      if (lpencil(i_uglnTT)) &
        call u_dot_grad(f,ilnTT,p%glnTT,p%uu,p%uglnTT,UPWIND=lupw_lnTT)
!
! fpres
!
      if (lpencil(i_fpres)) then
        do j=1,3
          p%fpres(:,j)=-gamma_m1*gamma1*(p%TT*p%glnrho(:,j) + p%glnTT(:,j))
        enddo
      endif
!
    endsubroutine calc_pencils_entropy
!**********************************************************************
    subroutine dss_dt(f,df,p)
!
!  calculate right hand side of entropy equation
!  heat condution is currently disabled until old stuff,
!  which in now in calc_heatcond, has been reinstalled.
!  DTT/Dt = -gamma_m1*TT*divu + gamma*cp1*rho1*RHS
!
!  13-dec-02/axel+tobi: adapted from entropy
!
      use Cdata
      use Diagnostics
      use Mpicomm
      use Sub
      use Viscosity, only: calc_viscous_heat
      use EquationOfState, only: gamma_m1
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: Hmax=0.
!
      intent(inout) :: f,p
      intent(out) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE dlnTT_dt'
      if (headtt) call identify_bcs('lnTT',ilnTT)
      if (headtt) print*,'dss_dt: lnTT,cs2=', p%lnTT(1), p%cs2(1)
!
!  sound speed squared
!
      if (headtt) print*,'dss_dt: cs20=',p%cs2(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
      if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!  subtract pressure gradient term in momentum equation
!
      if (lhydro.and.lpressuregradient_gas) &
         df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + p%fpres(:,iux:iuz)
!
!  advection term
!
      if (ladvection_temperature) &
         df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%uglnTT
!
!  Calculate viscous contribution to temperature
!
      if (lviscosity.and.lviscosity_heat) call calc_viscous_heat(df,p,Hmax)
!
!  Thermal conduction
!
      if (lheatc_chiconst) call calc_heatcond_constchi(df,p)
      if (lheatc_Kconst)   call calc_heatcond_constK(df,p)
      if (lheatc_Kprof)    call calc_heatcond(f,df,p)
      if (lheatc_Karctan)  call calc_heatcond_arctan(df,p)
!
!  Need to add left-hand-side of the continuity equation (see manual)
!  Check this
      if (ldensity) &
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - gamma_m1*p%TT*p%divu
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        if (idiag_TTmax/=0) call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        if (idiag_TTm/=0)   call sum_mn_name(p%TT,idiag_TTm)
        if (idiag_ethm/=0)   call sum_mn_name(p%ee/p%rho1,idiag_ethm)
        if (idiag_ssm/=0)   call sum_mn_name(p%ss,idiag_ssm)
        if (idiag_dtc/=0) &
          call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        if (idiag_eem/=0) call sum_mn_name(p%ee,idiag_eem)
        if (idiag_ppm/=0) call sum_mn_name(p%pp,idiag_ppm)
        if (idiag_csm/=0) call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine single_polytrope(f)
!
! 04-aug-2007/dintrans: a simple polytrope with index mpoly0
!
      use Cdata
      use Gravity, only: gravz
      use EquationOfState, only: cs20, lnrho0, gamma, gamma_m1, get_cp1, &
                                 cs2bot, cs2top
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: beta, zbot, ztop, cp1, T0, temp
!
!  beta is the (negative) temperature gradient
!  beta = -(g/cp) /[(1-1/gamma)*(m+1)]
!
      call get_cp1(cp1)
      beta=-cp1*gravz/(mpoly0+1.)*gamma/gamma_m1
      ztop=xyz0(3)+Lxyz(3)
      zbot=xyz0(3)
      T0=cs20/gamma_m1
      print*, 'polytrope: mpoly0, beta, T0=', mpoly0, beta, T0
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        temp=T0+beta*(ztop-z(n))
        f(:,m,n,ilnTT)=temp
        f(:,m,n,ilnrho)=lnrho0+mpoly0*log(temp)-mpoly0*log(T0)
      enddo
      cs2bot=gamma_m1*(T0+beta*(ztop-zbot))
      cs2top=cs20
!
    endsubroutine single_polytrope
!***********************************************************************
    subroutine rad_equil(f)
!
! 16-mai-2007/tgastine+dintrans: compute the radiative and hydrostatic 
! equilibria for a given radiative profile (here a hole for the moment)
!
      use Cdata
      use Gravity, only: gravz
      use EquationOfState, only:cs20,lnrho0,cs2top,cs2bot,gamma,gamma_m1
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mz) :: temp,lnrho
      real :: arg,hcond,dtemp,dlnrho
      real :: alp,sig,ecart
      integer :: i
!
      sig=7.
      ecart=0.4
      alp=(Kmax-Kmin)/(pi/2.+atan(sig*ecart**2))
      print*,'sig, ecart, alp=',sig, ecart, alp
      print*,'Kmin, Kmax, Fbot, Tbump=', Kmin, Kmax, Fbot, Tbump
!
! integrate from the top to the bottom
!
      temp(n2)=cs20/gamma_m1
      lnrho(n2)=lnrho0
      f(:,:,n2,ilnTT)=cs20/gamma_m1
      f(:,:,n2,ilnrho)=lnrho0
      do i=n2-1,n1,-1
        arg=sig*(temp(i+1)-Tbump-ecart)*(temp(i+1)-Tbump+ecart)
        hcond=Kmax+alp*(-pi/2.+atan(arg))
        dtemp=Fbot/hcond
        temp(i)=temp(i+1)+dz*dtemp
        dlnrho=2d0*(-gamma/gamma_m1*gravz-dtemp)/(7.d0/6.d0*temp(i)+5.d0/6.d0*temp(i+1))
        lnrho(i)=lnrho(i+1)+dz*dlnrho
        f(:,:,i,ilnTT)=temp(i)
        f(:,:,i,ilnrho)=lnrho(i)
      enddo
! initialize cs2bot by taking into account the new bottom value of temperature
! note: cs2top is already defined in eos_init by assuming cs2top=cs20
!      cs2bot=gamma_m1*temp(n1)
      print*,'cs2top, cs2bot=', cs2top, cs2bot
!
      if (lroot) then
        print*,'--> write the initial setup in data/proc0/setup.dat'
        open(unit=11,file=trim(directory)//'/setup.dat')
        write(11,'(4a14)') 'z','rho','temp','hcond'
        do i=n2,n1,-1
          arg=sig*(temp(i)-Tbump-ecart)*(temp(i)-Tbump+ecart)
          hcond=Kmax+alp*(-pi/2.+atan(arg))
          write(11,'(4e14.5)') z(i),exp(lnrho(i)),temp(i),hcond
        enddo
        close(11)
      endif
!
    endsubroutine rad_equil
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
!  01-mar-07/dintrans: adapted from temperature_ionization
!
!  calculate chi*grad(rho*T*glnTT)/(rho*TT)
!           =chi*(g2.glnTT+g2lnTT) where g2=glnrho+glnTT
!
      use Diagnostics
      use EquationOfState, only: gamma
      use Sub

      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: g2
!
      call dot(p%glnTT+p%glnrho,p%glnTT,g2)
!
!  Add heat conduction to RHS of temperature equation
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chi*(g2 + p%del2lnTT)
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chi*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif

    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_constK(df,p)
!
!  calculate gamma*K/rho/cp*div(grad TT)= gamma*K/rho/cp*del2 TT
!  Be careful! here lnTT means TT...
!
      use Diagnostics
      use EquationOfState, only: gamma
      use Sub

      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: hcond,chix
!
      hcond=hcond0
!
!  Add heat conduction to RHS of temperature equation
!
      chix=p%rho1*hcond*p%cp1
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*p%del2lnTT
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif

    endsubroutine calc_heatcond_constK
!***********************************************************************
    subroutine calc_heatcond_arctan(df,p)
!
! 16-mai-2007/tgastine+dintrans: radiative diffusion with an arctan
!  profile for the conductivity
!  calculate gamma/rho*cp*div(K T*grad lnTT)=
!    gamma*K/rho*cp*(gradlnTT.gradlnTT + del2ln TT + gradlnTT.gradlnK)
!
      use Diagnostics
      use EquationOfState, only: gamma
      use HDF5IO, only: output_profile
      use Sub

      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: arg,hcond,chiT,g2,chix
      real, dimension (nx,3) :: glnhcond=0.
      real :: alp,sig,ecart
!
! todo: must be defined in the namelist (used by rad_equil and _arctan...)
      sig=7.
      ecart=0.4
      alp=(Kmax-Kmin)/(pi/2.+atan(sig*ecart**2))
!
      arg=sig*(p%TT-Tbump-ecart)*(p%TT-Tbump+ecart)
      hcond=Kmax+alp*(-pi/2.+atan(arg))
      chiT=2./hcond*sig*alp*(p%TT-Tbump)/(1.+arg**2)  ! d(ln K)/dT
!
      glnhcond(:,3)=chiT*p%glnTT(:,3)
!      call dot(p%glnTT,p%glnTT,g1)
      call dot(p%glnTT,glnhcond,g2)
!
!  Write out hcond z-profile (during first time step only)
!
      if (m==m1) then
        call output_profile('hcond',(/z(n)/),(/hcond(1)/),'z', lsave_name=(n==n1))
        call output_profile('glnhcond',(/z(n)/),(/glnhcond(1,3)/),'z', lsave_name=(n==n1))
        call output_profile('K_T',(/z(n)/),(/chiT/),'z', lsave_name=(n==n1))
      endif
!
!  Add heat conduction to RHS of temperature equation
!
      chix=p%rho1*hcond*p%cp1
!      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*(g1+g2+p%del2lnTT)
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*(g2+p%del2lnTT)
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif

    endsubroutine calc_heatcond_arctan
!***********************************************************************
    subroutine calc_heatcond(f,df,p)
!
!  12-Mar-2007/dintrans: coded
!  calculate gamma*K/rho*cp*div(T*grad lnTT)= 
!              gamma*K/rho*cp*(gradlnTT.gradln(hcond*TT) + del2ln TT)
!
!
      use Diagnostics
      use EquationOfState, only: gamma
      use Sub

      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: g2,hcond,chix
      real, dimension (nx,3) :: glhc=0.,glnThcond
      integer :: i
      logical :: lwrite_hcond=.true.
      save :: lwrite_hcond
!
      if (lhcond_global) then
        hcond=f(l1:l2,m,n,iglobal_hcond)
        glhc(:,1)=f(l1:l2,m,n,iglobal_glhc)
      else
        hcond = 1. + (hcond1-1.)*step(rcyl_mn,r_bcz,-widthlnTT)
        hcond = hcond0*hcond
        glhc(:,1) = hcond0*(hcond1-1.)*der_step(rcyl_mn,r_bcz,-widthlnTT)
      endif
      if (lroot .and. lwrite_hcond) then
        open(1,file=trim(directory)//'/hcond.dat',position='append')
        write(1,'(3e14.5)') (rcyl_mn(i),hcond(i),glhc(i,1),i=1,nx)
        close(1)
        lwrite_hcond=.false.
      endif
!
      glnThcond = p%glnTT + glhc/spread(hcond,2,3)    ! grad ln(T*hcond)
      call dot(p%glnTT,glnThcond,g2)
!
!  Add heat conduction to RHS of temperature equation
!
      chix=p%rho1*hcond*p%cp1
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*(g2 + p%del2lnTT)
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif

    endsubroutine calc_heatcond
!***********************************************************************
    subroutine read_entropy_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_init_pars, IOSTAT=iostat)
!
    endsubroutine read_entropy_init_pars
!***********************************************************************
    subroutine write_entropy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_init_pars)
!
    endsubroutine write_entropy_init_pars
!***********************************************************************
    subroutine read_entropy_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_run_pars, IOSTAT=iostat)
!
    endsubroutine read_entropy_run_pars
!***********************************************************************
    subroutine write_entropy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_run_pars)
!
    endsubroutine write_entropy_run_pars
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
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
        idiag_ethm=0; idiag_ssm=0
        idiag_dtchi=0; idiag_dtc=0
        idiag_eem=0; idiag_ppm=0; idiag_csm=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'csm',idiag_csm)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        call farray_index_append('nname',nname)
        call farray_index_append('ilnTT',ilnTT)
        call farray_index_append('iss',iss)
        call farray_index_append('i_TTmax',idiag_TTmax)
        call farray_index_append('i_TTmin',idiag_TTmin)
        call farray_index_append('i_TTm',idiag_TTm)
        call farray_index_append('i_ethm',idiag_ethm)
        call farray_index_append('i_ssm',idiag_ssm)
        call farray_index_append('i_dtchi',idiag_dtchi)
        call farray_index_append('i_dtc',idiag_dtc)
        call farray_index_append('i_eem',idiag_eem)
        call farray_index_append('i_ppm',idiag_ppm)
        call farray_index_append('i_csm',idiag_csm)
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine get_slices_entropy(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_entropy
!***********************************************************************
    subroutine ADI_constK(finit,f)
       
      use Cdata
      use Cparam
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top

      implicit none

      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(mx,mz) :: finter,source,rho
      real, dimension(nx)    :: a,b,c
      real, dimension(nz)    :: rhs,work
      real    :: aalpha, bbeta
!
      source=(f(:,4,:,ilnTT)-finit(:,4,:,ilnTT))/dt
      rho=exp(f(:,4,:,ilnrho))
      rho=rho/gamma
!
!  lignes en implicite
!
      do j=n1,n2
        a=-hcond0*dt/(2.*rho(l1:l2,j)*dx**2)
!
        b=1.+hcond0*dt/(rho(l1:l2,j)*dx**2)
!
        c=-hcond0*dt/(2.*rho(l1:l2,j)*dx**2) 
!
        rhs=finit(l1:l2,4,j,ilnTT)+hcond0*dt/(2.*rho(l1:l2,j)*dz**2) &
            *(finit(l1:l2,4,j+1,ilnTT)-2.*finit(l1:l2,4,j,ilnTT)+ &
            finit(l1:l2,4,j-1,ilnTT))+dt/2.*source(l1:l2,j)
!
        aalpha=c(nx)
        bbeta=a(1)
        c(nx)=0.
        a(1)=0.
        call cyclic(a,b,c,aalpha,bbeta,rhs,work,nx)
        finter(l1:l2,j)=work(1:nx)
      enddo
!
      call BC_CT(finter)
!
!  colonnes en implicite
!
      do i=l1,l2
        a=-hcond0*dt/(2.*rho(i,n1:n2)*dz**2)
!
        b=1.+hcond0*dt/(rho(i,n1:n2)*dz**2)
!
        c=-hcond0*dt/(2.*rho(i,n1:n2)*dz**2)
!
        rhs=finter(i,n1:n2)+hcond0*dt/(2.*rho(i,n1:n2)*dx**2) &
           *(finter(i+1,n1:n2)-2.*finter(i,n1:n2)+finter(i-1,n1:n2)) &
           +dt/2.*source(i,n1:n2)
!
        c(nz)=0.
        a(1)=0.
        b(1)=1.
        b(nz)=1.
        c(1)=0.
        a(nz)=0.
        rhs(1)=cs2bot/gamma_m1
        rhs(nx)=cs2top/gamma_m1
        call tridag(a,b,c,rhs,work,nz)
        f(i,4,n1:n2,ilnTT)=work(1:nz)
      enddo
!
      call BC_CT(f(:,4,:,ilnTT))
!
      endsubroutine ADI_constK
!**************************************************************
      subroutine ADI_Kprof(finit,f)
       
      use Cdata
      use Cparam
      use EquationOfState, only: gamma,cs2bot,cs2top

      implicit none

      integer :: i,j
      real    :: aalpha, bbeta
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(mx,mz) :: source,rho,chiprof,dchi,valinter,val
      real, dimension(nx)    :: a,b,c
      real, dimension(nz)    :: rhs,work

      source=(f(:,4,:,ilnTT)-finit(:,4,:,ilnTT))/dt
      rho=exp(f(:,4,:,ilnrho))
      rho=rho/gamma
      call hcond_ADI(f,chiprof,dchi)
!
!  lignes en implicite
!
      do j=n1,n2
       a=-dt/(4d0*rho(l1:l2,j)*dx**2)*(dchi(l1-1:l2-1,j) &
         *(finit(l1-1:l2-1,4,j,ilnTT)-finit(l1:l2,4,j,ilnTT)) &
         +chiprof(l1-1:l2-1,j)+chiprof(l1:l2,j))
!
       b=1d0+dt/(4d0*rho(l1:l2,j)*dx**2)*(dchi(l1:l2,j) &
         *(2d0*finit(l1:l2,4,j,ilnTT)-finit(l1-1:l2-1,4,j,ilnTT) &
         -finit(l1+1:l2+1,4,j,ilnTT))+2d0*chiprof(l1:l2,j) &
       +chiprof(l1+1:l2+1,j)+chiprof(l1-1:l2-1,j))
!
       c=-dt/(4d0*rho(l1:l2,j)*dx**2)*(dchi(l1+1:l2+1,j) &
          *(finit(l1+1:l2+1,4,j,ilnTT)-finit(l1:l2,4,j,ilnTT)) &
       +chiprof(l1:l2,j)+chiprof(l1+1:l2+1,j))
!
       rhs=1d0/(2d0*rho(l1:l2,j)*dz**2)*((chiprof(l1:l2,j+1) &
           +chiprof(l1:l2,j))*(finit(l1:l2,4,j+1,ilnTT)-finit(l1:l2,4,j,ilnTT))&
           -(chiprof(l1:l2,j)+chiprof(l1:l2,j-1)) &
       *(finit(l1:l2,4,j,ilnTT)-finit(l1:l2,4,j-1,ilnTT)))
!
       rhs=rhs+1d0/(2d0*rho(l1:l2,j)*dx**2)*((chiprof(l1+1:l2+1,j) &
         +chiprof(l1:l2,j))*(finit(l1+1:l2+1,4,j,ilnTT)-finit(l1:l2,4,j,ilnTT))&
           -(chiprof(l1:l2,j)+chiprof(l1-1:l2-1,j)) &
           *(finit(l1:l2,4,j,ilnTT)-finit(l1-1:l2-1,4,j,ilnTT)))
!
       aalpha=c(nx)
       bbeta=a(1)
       c(nx)=0d0
       a(1)=0d0
       call cyclic(a,b,c,aalpha,bbeta,rhs,work,nx)
       valinter(l1:l2,j)=work(1:nx)
      enddo
!
!  colonnes en implicite
!
      do i=l1,l2
       a=-dt/(4d0*rho(i,n1:n2)*dz**2)*(dchi(i,n1-1:n2-1) &
         *(finit(i,4,n1-1:n2-1,ilnTT)-finit(i,4,n1:n2,ilnTT))&
         +chiprof(i,n1-1:n2-1)+chiprof(i,n1:n2))
!
       b=1d0+dt/(4d0*rho(i,n1:n2)*dz**2)*(dchi(i,n1:n2)* &
         (2d0*finit(i,4,n1:n2,ilnTT)-finit(i,4,n1-1:n2-1,ilnTT) &
         -finit(i,4,n1+1:n2+1,ilnTT))+2d0*chiprof(i,n1:n2) &
         +chiprof(i,n1+1:n2+1)+chiprof(i,n1-1:n2-1))
!
       c=-dt/(4d0*rho(i,n1:n2)*dz**2)*(dchi(i,n1+1:n2+1) &
         *(finit(i,4,n1+1:n2+1,ilnTT)-finit(i,4,n1:n2,ilnTT))&
         +chiprof(i,n1:n2)+chiprof(i,n1+1:n2+1))
!
       rhs=valinter(i,n1:n2)
!
       c(nz)=0d0
       a(1)=0d0
! Constant flux at the bottom
!       b(1)=-1.
!       c(1)=0.
!       c(1)=1.
!       rhs(1)=0.
! Constant temperature at the top
       b(nx)=1.
       a(nx)=0.
       rhs(nz)=0.
!
       call tridag(a,b,c,rhs,work,nz)
       val(i,n1:n2)=work(1:nz)
      enddo
!
      f(:,4,:,ilnTT)=finit(:,4,:,ilnTT)+dt*val+dt*source
!
!      f(:,:,n1,ilnTT)=cs2bot/(gamma-1d0)
      f(:,:,n2,ilnTT)=cs2top/(gamma-1d0)
!      call BC_CT(f)
      call BC_flux(f,chiprof)
      call hcond_ADI(f,chiprof,dchi)
!
      endsubroutine ADI_Kprof
!**************************************************************
      subroutine BC_CT(f_2d)

      implicit none

      real, dimension(mx,mz) :: f_2d

! z-direction
      f_2d(:,n1-1)=2.*f_2d(:,n1)-f_2d(:,n1+1)
      f_2d(:,n2+1)=2.*f_2d(:,n2)-f_2d(:,n2-1)
! x-direction
      f_2d(1:l1-1,:)=f_2d(l2i:l2,:)
      f_2d(l2+1:mx,:)=f_2d(l1:l1i,:)

      endsubroutine BC_CT
!**************************************************************
      subroutine BC_flux(f,chiprof)

      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz) :: chiprof
      integer :: i

      do i=1,nghost
        f(:,:,n1-i,ilnTT)=f(:,:,n1+i,ilnTT)+2*i*dz*Fbot/chiprof(10,n1+i)
      enddo
      f(:,:,n2+1,ilnTT)=2*f(:,:,n2,ilnTT)-f(:,:,n2-1,ilnTT)
      ! x-direction
      f(1:l1-1,:,:,ilnTT)=f(l2i:l2,:,:,ilnTT)
      f(l2+1:mx,:,:,ilnTT)=f(l1:l1i,:,:,ilnTT)

      endsubroutine BC_flux
!**************************************************************
      subroutine hcond_ADI(f,chiprof,dchi)

      use Sub, only: write_zprof
      real, dimension(mx,my,mz,mfarray) :: f
      real , dimension(mx,mz) :: chiprof, dchi
!      double precision :: chi0
      real, dimension(mx,mz) :: arg
      real :: alp,sig,ecart
!
!      chi0=1d-3
!      chiprof=chi0
!      dchi=0d0
      sig=7.
      ecart=0.4
      alp=(Kmax-Kmin)/(pi/2.+atan(sig*ecart**2))
      arg=sig*(f(:,4,:,ilnTT)-Tbump-ecart)*(f(:,4,:,ilnTT)-Tbump+ecart)
      chiprof=Kmax+alp*(-pi/2.+atan(arg))
      dchi=2d0*alp/(1d0+arg**2)*sig*(f(:,4,:,ilnTT)-Tbump)
!
!      call write_zprof('hcond',chiprof(10,n1:n2))
!      call write_zprof('dhcond/dT',dchi(10,n1:n2))

      endsubroutine hcond_ADI
!**************************************************************
      subroutine tridag(a,b,c,r,u,n)

      integer :: j,n
      integer, parameter :: NMAX=500
      real    :: bet
      real, dimension(n) :: a,b,c,r,u
      real, dimension(NMAX) :: gam
!
      if (b(1).eq.0.) pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if (bet.eq.0.) pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
!
      return
      endsubroutine tridag
!**************************************************************
      subroutine cyclic(a,b,c,alpha,beta,r,x,n)
!
      implicit none
!
      integer :: i,n
      integer, parameter    :: NMAX=500
      real    :: alpha, beta,gamma,fact      
      real, dimension(n)    :: a,b,c,r,x
      real, dimension(NMAX) :: bb,u,z
!
      if (n.le.2)pause 'n too small in cyclic'
      if (n.gt.NMAX)pause 'NMAX too small in cyclic'
      gamma=-b(1)
      bb(1)=b(1)-gamma
      bb(n)=b(n)-alpha*beta/gamma
      do 11 i=2,n-1
        bb(i)=b(i)
11    continue
      call tridag(a,bb,c,r,x,n)
      u(1)=gamma
      u(n)=alpha
      do 12 i=2,n-1
        u(i)=0.
12    continue
      call tridag(a,bb,c,u,z,n)
      fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
      do 13 i=1,n
        x(i)=x(i)-fact*z(i)
13    continue
!
      return
      endsubroutine cyclic
!***************************************************************
    subroutine calc_heatcond_ADI(finit,f)
!
      implicit none
!
      real, dimension(mx,my,mz,mfarray) :: finit,f
!
      call keep_compiler_quiet(finit)
      call keep_compiler_quiet(f)
!
    endsubroutine calc_heatcond_ADI
!***********************************************************************
endmodule Entropy
