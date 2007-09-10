! $Id: temperature_idealgas.f90,v 1.28 2007-09-10 14:16:08 dintrans Exp $
!  This module can replace the entropy module by using lnT or T (with
!  ltemperature_nolog=.true.) as dependent variable. For a perfect gas 
!  with constant coefficients (no ionization) we have:
!  (1-1/gamma) * cp*T = cs02 * exp( (gamma-1)*ln(rho/rho0)-gamma*s/cp )
!
!  Note that to use lnTT as thermal variable, you may rather want to use
!  entropy.f90 with pretend_lnTT=.true. As of March 2007, entropy.f90
!  has way more options and features than temperature_idealgas.f90. 
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
! PENCILS PROVIDED Ma2,uglnTT,fpres
!
!***************************************************************
module Entropy

  use Cparam
  use Cdata
  use Messages
  use Interstellar
  use EquationOfState, only: mpoly0,mpoly1

  implicit none

  public :: calc_heatcond_ADI

  include 'entropy.h'

  interface heatcond_TT ! Overload subroutine `hcond_TT' function
    module procedure heatcond_TT_2d     ! get 2d-arrays (hcond, dhcond)
    module procedure heatcond_TT_point  ! get one value (hcond, dhcond)
  end interface

  real :: radius_lnTT=0.1,ampl_lnTT=0.,widthlnTT=2*epsi
  real :: lnTT_left=1.0,lnTT_right=1.0,lnTT_const=0.0,TT_const=1.0
  real :: kx_lnTT=1.0,ky_lnTT=1.0,kz_lnTT=1.0
  real :: chi=impossible,heat_uniform=0.0,difflnTT_hyper=0.
  real :: zbot=0.0,ztop=0.0
  real :: tau_heat_cor=-1.0,tau_damp_cor=-1.0,zcor=0.0,TT_cor=0.0
  real :: center1_x=0., center1_y=0., center1_z=0.
  real :: r_bcz=0.
! entries for ADI
  real :: Tbump=0.,Kmin=0.,Kmax=0.,hole_slope=0.,hole_width=0.
  real :: hole_alpha ! initialized _after_ the reading
  integer, parameter :: nheatc_max=2
  logical :: lpressuregradient_gas=.true.,ladvection_temperature=.true.
  logical :: lupw_lnTT=.false.,lcalc_heat_cool=.false.,ldiff_hyper=.false.
  logical :: lheatc_Kconst=.false.,lheatc_Kprof=.false.,lheatc_Karctan=.false.
  logical :: lheatc_chiconst=.false.,lheatc_chiconst_accurate=.false.
  logical :: lfreeze_lnTTint=.false.,lfreeze_lnTText=.false.
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  logical :: lhcond_global=.false.
  logical :: lviscosity_heat=.true.
  integer :: iglobal_hcond=0
  integer :: iglobal_glhc=0

  character (len=labellen), dimension(ninit) :: initlnTT='nothing'
  character (len=5) :: iinit_str

! Delete (or use) me asap!
  real :: hcond0=impossible, hcond1=1.,Fbot,FbotKbot,Ftop,Kbot,FtopKtop
  logical :: lmultilayer=.false.

! input parameters
  namelist /entropy_init_pars/ &
      initlnTT,radius_lnTT,ampl_lnTT,widthlnTT, &
      lnTT_left,lnTT_right,lnTT_const,TT_const, &
      kx_lnTT,ky_lnTT,kz_lnTT,center1_x,center1_y,center1_z, &
      mpoly0,mpoly1,r_bcz, &
      Fbot,Tbump,Kmin,Kmax,hole_slope,hole_width

! run parameters
  namelist /entropy_run_pars/ &
      lupw_lnTT,lpressuregradient_gas,ladvection_temperature, &
      heat_uniform,chi,iheatcond,tau_heat_cor,tau_damp_cor,zcor,TT_cor, &
      lheatc_chiconst_accurate,hcond0,lcalc_heat_cool,&
      lfreeze_lnTTint,lfreeze_lnTText,widthlnTT,mpoly0,mpoly1, &
      lhcond_global,lviscosity_heat,difflnTT_hyper, &
      Fbot,Tbump,Kmin,Kmax,hole_slope,hole_width
!
! other variables (needs to be consistent with reset list below)
  integer :: idiag_TTmax=0,idiag_TTmin=0,idiag_TTm=0
  integer :: idiag_yHmax=0,idiag_yHmin=0,idiag_yHm=0
  integer :: idiag_eth=0,idiag_ssm=0,idiag_thcool=0
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
           "$Id: temperature_idealgas.f90,v 1.28 2007-09-10 14:16:08 dintrans Exp $")
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
      use Cdata
      use FArrayManager
      use Gravity, only: g0
      use EquationOfState
      use Sub, only: step,der_step
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx)   :: hcond,dhcond
      logical :: lstarting, lnothing
      type (pencil_case) :: p
      integer :: i
!
      if (.not. leos) then
         call fatal_error('initialize_entropy','EOS=noeos but temperature_idealgas requires an EQUATION OF STATE for the fluid')
      endif
!
      if (ltemperature_nolog) then
        call select_eos_variable('TT',ilnTT)
      else
        call select_eos_variable('lnTT',ilnTT)
      endif
!
!  freeze temperature
!
      if (lfreeze_lnTTint) lfreeze_varint(ilnTT)=.true.
      if (lfreeze_lnTText) lfreeze_varext(ilnTT)=.true.
!
!  Check whether we want heat conduction
!
      lheatc_Kconst= .false.
      lheatc_Kprof= .false.
      lheatc_Karctan= .false.
      lheatc_chiconst = .false.
      ldiff_hyper = .false.
      lnothing = .false.
!
      do i=1,nheatc_max
      select case (iheatcond(i))
        case('K-const')
          lheatc_Kconst=.true.
          if (lroot) call information('initialize_entropy', &
          ' heat conduction: K=cst --> gamma*K/rho/TT/cp*div(T*grad lnTT)')
          if (initlnTT(1).ne.'rad_equil') &
            Fbot=gamma/(gamma-1.)*hcond0*g0/(mpoly0+1.)
        case('K-profile')
          lheatc_Kprof=.true.
! 
!  TODO..... ailleurs !
!
          hcond1=(mpoly1+1.)/(mpoly0+1.)
          Fbot=gamma/(gamma-1.)*hcond0*g0/(mpoly0+1.)
          if (lroot) call information('initialize_entropy',' heat conduction: K=K(r)')
        case('K-arctan')
          lheatc_Karctan=.true.         
          if (lroot) call information('initialize_entropy',' heat conduction: arctan profile')
        case('chi-const')
          lheatc_chiconst=.true.
          if (lroot) call information('initialize_entropy',' heat conduction: constant chi')
        case ('hyper')
          ldiff_hyper=.true.
          if (lroot) call information('initialize_entropy','hyper diffusion')
        case ('nothing')
          if (lroot .and. (.not. lnothing)) print*,'heat conduction: nothing'
        case default
          if (lroot) then
            write(unit=errormsg,fmt=*)  &
                'No such value iheatcond = ', trim(iheatcond(i))
            call fatal_error('initialize_entropy',errormsg)
          endif
       endselect
       enddo
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
!  Some initializations for the ADI setup
!
      if (hole_slope.ne.0.) then
        hole_alpha=(Kmax-Kmin)/(pi/2.+atan(hole_slope*hole_width**2))
        print*,'hole_slope, hole_width, hole_alpha=',hole_slope, &
             hole_width, hole_alpha
        print*,'Kmin, Kmax, Fbot, Tbump=', Kmin, Kmax, Fbot, Tbump
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
      if (iheatcond(1)=='nothing') then
        if (hcond0 /= impossible) call warning('initialize_entropy', 'No heat conduction, but hcond0 /= 0')
        if (chi /= impossible) call warning('initialize_entropy', 'No heat conduction, but chi /= 0')
      endif

    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f,xx,yy,zz)
!
!  initialise lnTT; called from start.f90
!
!  13-dec-2002/axel+tobi: adapted from init_ss
!
!  initialise entropy; called from start.f90
!  07-nov-2001/wolf: coded 
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!
      use General,  only: chn
      use Sub,      only: blob
      use Initcond, only: jump
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
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
      case('single_polytrope'); call single_polytrope(f)
      case('rad_equil')
          call rad_equil(f)
          if (ampl_lnTT.ne.0.) then
            print*,'add a bubble with:',ampl_lnTT,radius_lnTT,center1_x,center1_y,center1_z
            call blob(ampl_lnTT,f,ilnTT,radius_lnTT,center1_x,center1_y,center1_z)
            call blob(-ampl_lnTT,f,ilnrho,radius_lnTT,center1_x,center1_y,center1_z)
          endif
      case('bubble_hs')
!         print*,'init_lnTT: put bubble in hydrostatic equilibrium: radius_lnTT,ampl_lnTT=',radius_lnTT,ampl_lnTT,center1_x,center1_y,center1_z
          call blob(ampl_lnTT,f,ilnTT,radius_lnTT,center1_x,center1_y,center1_z)
          call blob(-ampl_lnTT,f,ilnrho,radius_lnTT,center1_x,center1_y,center1_z)
         !
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
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_visc_heat)=.true.
      endif
!
      if (ldensity) then
         lpenc_requested(i_divu)=.true.
      endif
!
      if (lcalc_heat_cool) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_cv1)=.true.
      endif
!
      if (lheatc_chiconst) then
        lpenc_requested(i_del2lnTT)=.true.
      endif
      if (lheatc_Kconst) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
      if (lheatc_Kprof) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
      if (lheatc_Karctan) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
      if (ldiff_hyper) then
         lpenc_requested(i_del6lnTT)=.true.
      endif

!
      if (ladvection_temperature) lpenc_requested(i_uglnTT)=.true.
!
      !if (lheatc_shock) then
      !   lpenc_requested(i_glnrho)=.true.
      !   lpenc_requested(i_gss)=.true.
      !   lpenc_requested(i_del2lnTT)=.true.
      !   lpenc_requested(i_gshock)=.true.
      !   lpenc_requested(i_shock)=.true.
      !   lpenc_requested(i_glnTT)=.true.
      !endif
!
!  Diagnostics
!
      if (idiag_TTmax/=0) lpenc_diagnos(i_TT)  =.true.
      if (idiag_TTmin/=0) lpenc_diagnos(i_TT)  =.true.
      if (idiag_TTm/=0)   lpenc_diagnos(i_TT)  =.true.
      if (idiag_yHmax/=0) lpenc_diagnos(i_yH)  =.true.
      if (idiag_yHmin/=0) lpenc_diagnos(i_yH)  =.true.
      if (idiag_yHm/=0)   lpenc_diagnos(i_yH)  =.true.
      if (idiag_eth/=0) then
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
      if (idiag_thcool/=0) lpenc_diagnos(i_rho)=.true.
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

      if (lpencil_in(i_uglnTT)) then
        lpencil_in(i_glnTT)=.true.
      endif
!
      if (lpencil_in(i_fpres)) then
        if (ltemperature_nolog) then
          lpencil_in(i_TT)=.true.
        else
          lpencil_in(i_cs2)=.true.
        endif
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_glnTT)=.true.
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
          if (ltemperature_nolog) then
            p%fpres(:,j)=-gamma1*gamma11*(p%TT*p%glnrho(:,j) + p%glnTT(:,j))
          else
            p%fpres(:,j)=-p%cs2*(p%glnrho(:,j) + p%glnTT(:,j))*gamma11
          endif
        enddo
      endif
!
    endsubroutine calc_pencils_entropy
!**********************************************************************
    subroutine dss_dt(f,df,p)
!
!  calculate right hand side of temperature equation
!  heat condution is currently disabled until old stuff,
!  which in now in calc_heatcond, has been reinstalled.
!  lnTT version: DlnTT/Dt = -gamma1*divu + gamma*cp1*rho1*TT1*RHS
!    TT version:   DTT/Dt = -gamma1*TT*divu + gamma*cp1*rho1*RHS
!
!  13-dec-02/axel+tobi: adapted from entropy
!
      use Cdata
      use Mpicomm
      use Sub
      use Global
      use Viscosity, only: calc_viscous_heat
      use EquationOfState, only: gamma1
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: Hmax=0.
      integer :: j,ju
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
!  entropy gradient: needed for advection and pressure gradient
!
      !call grad(f,ilnTT,glnTT)
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
!  advection term and PdV-work
!
      if (ladvection_temperature) &
         df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%uglnTT
!
!  Calculate viscous contribution to temperature
!
      if (lviscosity.and.lviscosity_heat) call calc_viscous_heat(f,df,p,Hmax)
!
!  Various heating conduction contributions
!
      if (lcalc_heat_cool)  call calc_heat_cool(f,df,p)
!
!  Thermal conduction: only chi=cte for the moment
!
      if (lheatc_chiconst) call calc_heatcond_constchi(df,p)
      if (lheatc_Kconst)   call calc_heatcond_constK(df,p)
      if (lheatc_Kprof)    call calc_heatcond(f,df,p)
      if (lheatc_Karctan)  call calc_heatcond_arctan(df,p)
!
!  Hyper diffusion
!
      if (ldiff_hyper) then
         if(headtt) print*,'Hyper diffusion: difflnTT_hyper=',difflnTT_hyper
         !
         df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + difflnTT_hyper*p%del6lnTT
         !
         if (lfirst.and.ldt) diffus_chi=diffus_chi+difflnTT_hyper*dxyz_6
      endif
!
!  Need to add left-hand-side of the continuity equation (see manual)
!  Check this

      if (ldensity) then
        if (ltemperature_nolog) then
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - gamma1*p%TT*p%divu
        else
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - gamma1*p%divu
        endif
      endif
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        if (idiag_TTmax/=0) call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        if (idiag_TTm/=0)   call sum_mn_name(p%TT,idiag_TTm)
        if (idiag_eth/=0)   call sum_mn_name(p%ee/p%rho1,idiag_eth)
        if (idiag_ssm/=0)   call sum_mn_name(p%ss,idiag_ssm)
        if (idiag_dtc/=0) then
          call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        endif
        if (idiag_eem/=0) call sum_mn_name(p%ee,idiag_eem)
        if (idiag_ppm/=0) call sum_mn_name(p%pp,idiag_ppm)
        if (idiag_csm/=0) call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine rad_equil(f)
!
! 16-mai-07/gastine+dintrans: compute the radiative and hydrostatic 
! equilibria for a given radiative profile defined in heatcond_TT.
!
      use Cdata
      use Gravity, only: gravz
      use EquationOfState, only: lnrho0,cs20,cs2top,cs2bot,gamma, &
                                 gamma1,eoscalc,ilnrho_TT
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mz) :: temp,lnrho
      real :: hcond, dhcond, dtemp, dlnrho, ss
      integer :: i
!
      if (.not. ltemperature_nolog) &
        call fatal_error('temperature_idealgas',  &
                         'rad_equil not implemented for lnTT')
      if (lroot) print*,'init_ss: rad_equil for kappa-mechanism pb'
!
! Integrate from top to bottom: z(n2) --> z(n1)
!
      temp(n2)=cs20/gamma1
      lnrho(n2)=lnrho0
      f(:,:,n2,ilnTT)=cs20/gamma1
      f(:,:,n2,ilnrho)=lnrho0
!
! Calculate the n2-1 gridpoint thanks to a 1st order forward Euler scheme
!
      call heatcond_TT(temp(n2),hcond,dhcond)
      dtemp=Fbot/hcond
      temp(n2-1)=temp(n2)+dz*dtemp
      dlnrho=(-gamma/gamma1*gravz-dtemp)/temp(n2)
      lnrho(n2-1)=lnrho(n2)+dz*dlnrho
      f(:,:,n2-1,ilnTT)=temp(n2-1)
      f(:,:,n2-1,ilnrho)=lnrho(n2-1)
!
! Now we use a 2nd order centered scheme for the other gridpoints
!
      do i=n2-1,n1+1,-1
        call heatcond_TT(temp(i),hcond,dhcond)
        dtemp=Fbot/hcond
        temp(i-1)=temp(i+1)+2.*dz*dtemp
        dlnrho=(-gamma/gamma1*gravz-dtemp)/temp(i)
        lnrho(i-1)=lnrho(i+1)+2.*dz*dlnrho
        f(:,:,i-1,ilnTT)=temp(i-1)
        f(:,:,i-1,ilnrho)=lnrho(i-1)
      enddo
!
! Initialize cs2bot by taking into account the new bottom value of temperature
! Note: cs2top=cs20 already defined in eos_idealgas
      cs2bot=gamma1*temp(n1)
      print*,'cs2top, cs2bot=', cs2top, cs2bot
!
      if (lroot) then
        print*,'--> write the initial setup in data/proc0/setup.dat'
        open(unit=11,file=trim(directory)//'/setup.dat')
        write(11,'(5a14)') 'z','rho','temp','ss','hcond'
        do i=n2,n1,-1
          call eoscalc(ilnrho_TT,lnrho(i),temp(i),ss=ss)
          call heatcond_TT(temp(i),hcond,dhcond)
          write(11,'(5e14.5)') z(i),exp(lnrho(i)),temp(i),ss,hcond
        enddo
        close(11)
      endif
!
    endsubroutine rad_equil
!***********************************************************************
    subroutine calc_heat_cool(f,df,p)

      use EquationOfState, only: gamma,gamma1
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: tau,cooling,kappa,a1,a3,pTT
      real :: a2,kappa0,kappa0_cgs
!
!  Initialize
!
      intent(in) :: p
      intent(out) :: df
!
      if (headtt) print*,'enter calc_heat_cool'
!
      kappa0_cgs=2e-6  !cm2/g
      kappa0=kappa0_cgs*unit_density*unit_length
      kappa=kappa0*p%TT**2
!
!  Optical Depth tau=kappa*rho*H
!  If we are using 2D, the pencil value p%rho is actually
!   sigma, the column density, sigma=rho*2*H
!
      if (nzgrid==1) then
         tau = .5*kappa*p%rho
      else
         call fatal_error("calc_heat_cool","opacity not yet implemented for 3D")
      endif
!
! Analytical gray description of Hubeny (1990)
! a1 is the optically thick contribution,
! a3 the optically thin one.
!
      a1=0.375*tau ; a2=0.433013 ; a3=0.25/tau
!
! cooling for Energy: 2*sigmaSB*p%TT**4/(a1+a2+a3)
!
      cooling = 2*sigmaSB*p%rho1*p%TT**4/(a1+a2+a3)
!
!  this cooling has dimension of energy over time
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%cv1*p%TT1*cooling
!
      if (ldiagnos) then
         !cooling power - energy radiated away (luminosity)
         if (idiag_thcool/=0) call sum_lim_mn_name(cooling*p%rho,idiag_thcool,p)
      endif
!
    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
!  01-mar-07/dintrans: adapted from temperature_ionization
!
!  calculate chi*grad(rho*T*glnTT)/(rho*TT)
!           =chi*(g2.glnTT+g2lnTT) where g2=glnrho+glnTT
!
      use Sub, only: max_mn_name,dot,del2,multsv
      use EquationOfState, only: gamma

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
        diffus_chi=max(diffus_chi,gamma*chi*dxyz_2)
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif

    end subroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_constK(df,p)
!
! 
!  calculate gamma*K/rho/TT/cp*div(T*grad lnTT)= 
!            gamma*K/rho/cp*(gradlnTT.gradlnTT + del2ln TT)
!
!
      use Sub, only: max_mn_name,dot,del2,multsv
      use EquationOfState, only: gamma

      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: g2,hcond,chix
!
      hcond=hcond0
!
!  Add heat conduction to RHS of temperature equation
!
      chix=p%rho1*hcond*p%cp1
      if (ltemperature_nolog) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*p%del2lnTT
      else
        call dot(p%glnTT,p%glnTT,g2)
       df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*(g2 + p%del2lnTT)
      endif
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) then
        diffus_chi=max(diffus_chi,gamma*chix*dxyz_2)
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif

    endsubroutine calc_heatcond_constK
!***********************************************************************
    subroutine calc_heatcond_arctan(df,p)
!
! 16-mai-07/gastine+dintrans: radiative diffusion with an arctan
!  profile for the conductivity
!  calculate gamma/rho*cp*div(K T*grad lnTT)=
!    gamma*K/rho*cp*(gradlnTT.gradlnTT + del2ln TT + gradlnTT.gradlnK)
!
      use Sub, only: max_mn_name,dot,del2,multsv,write_zprof
      use EquationOfState, only: gamma

      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: arg,hcond,chiT,g1,g2,chix
      real, dimension (nx,3) :: glnhcond=0.
      real :: alp
!
! todo: must be defined in the namelist (used by rad_equil and _arctan...)
      alp=(Kmax-Kmin)/(pi/2.+atan(hole_slope*hole_width**2))
!
      arg=hole_slope*(p%TT-Tbump-hole_width)*(p%TT-Tbump+hole_width)
      hcond=Kmax+alp*(-pi/2.+atan(arg))
      chiT=2.*p%TT/hcond*hole_slope*alp*(p%TT-Tbump)/(1.+arg**2)  ! d(ln K)/dT
!
      glnhcond(:,3)=chiT*p%glnTT(:,3)
      call dot(p%glnTT,p%glnTT,g1)
      call dot(p%glnTT,glnhcond,g2)
!
!  Write out hcond z-profile (during first time step only)
!
      call write_zprof('hcond',hcond)
      call write_zprof('glnhcond',glnhcond(:,3))
      call write_zprof('K_T',chiT)
!
!  Add heat conduction to RHS of temperature equation
!
      chix=p%rho1*hcond*p%cp1
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*(g1+g2+p%del2lnTT)
!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) then
        diffus_chi=max(diffus_chi,gamma*chix*dxyz_2)
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
      use Sub, only: max_mn_name,dot,del2,multsv,step,der_step
      use EquationOfState, only: gamma

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
        diffus_chi=max(diffus_chi,gamma*chix*dxyz_2)
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif

    endsubroutine calc_heatcond
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
        idiag_eth=0; idiag_ssm=0; idiag_thcool=0
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
        call parse_name(iname,cname(iname),cform(iname),'eth',idiag_eth)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'csm',idiag_csm)
        call parse_name(iname,cname(iname),cform(iname),'thcool',idiag_thcool)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'ilnTT=',ilnTT
        write(3,*) 'iyH=',iyH
        write(3,*) 'iss=',iss
        write(3,*) 'i_TTmax=',idiag_TTmax
        write(3,*) 'i_TTmin=',idiag_TTmin
        write(3,*) 'i_TTm=',idiag_TTm
        write(3,*) 'i_eth=',idiag_eth
        write(3,*) 'i_ssm=',idiag_ssm
        write(3,*) 'i_dtchi=',idiag_dtchi
        write(3,*) 'i_dtc=',idiag_dtc
        write(3,*) 'i_eem=',idiag_eem
        write(3,*) 'i_ppm=',idiag_ppm
        write(3,*) 'i_csm=',idiag_csm
        write(3,*) 'i_thcool=',idiag_thcool
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************
    subroutine single_polytrope(f)
!
! 04-aug-2007/dintrans: a single polytrope with index mpoly0
!
      use Cdata
      use Gravity, only: gravz
      use EquationOfState, only: cs20, lnrho0, gamma, gamma1, get_cp1, &
                                 cs2bot, cs2top
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: beta, zbot, ztop, cp1, T0, temp
!
!  beta is the (negative) temperature gradient
!  beta = -(g/cp) /[(1-1/gamma)*(m+1)]
!
      call get_cp1(cp1)
      beta=-cp1*gravz/(mpoly0+1.)*gamma/gamma1
      ztop=xyz0(3)+Lxyz(3)
      zbot=xyz0(3)
      T0=cs20/gamma1
      print*, 'polytrope: mpoly0, beta, T0=', mpoly0, beta, T0
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        temp=T0+beta*(ztop-z(n))
        if (ltemperature_nolog) then
          f(:,m,n,ilnTT)=temp
        else
          f(:,m,n,ilnTT)=log(temp)
        endif
        f(:,m,n,ilnrho)=lnrho0+mpoly0*log(temp/T0)
      enddo
      cs2bot=gamma1*(T0+beta*(ztop-zbot))
      cs2top=cs20
!
    endsubroutine single_polytrope
!**************************************************************
    subroutine calc_heatcond_ADI(finit,f)
!
!  10-sep-07/gastine+dintrans: wrapper to the two possible ADI subroutines
!  ADI_Kconst: constant radiative conductivity
!  ADI_Kprof: radiative conductivity depends on T, i.e. hcond(T)
!
      use Cparam

      implicit none

      real, dimension(mx,my,mz,mfarray) :: finit,f
!
      if (hcond0 /= impossible) then
        call ADI_Kconst(finit,f)
      else
        call ADI_Kprof(finit,f)
      endif
!
    end subroutine calc_heatcond_ADI
!***********************************************************************
    subroutine ADI_Kconst(finit,f)
!
!  08-Sep-07/gastine+dintrans: coded
!  2-D ADI scheme for the radiative diffusion term (see e.g. 
!  Peaceman & Rachford 1955). Each direction are solved implicitly:
!
!    (1-dt/2*Lambda_x)*T^(n+1/2) = (1+dt/2*Lambda_y)*T^n
!    (1-dt/2*Lambda_y)*T^(n+1)   = (1+dt/2*Lambda_x)*T^(n+1/2)
!
!  where Lambda_x and Lambda_y denote diffusion operators.
!
      use Cdata
      use Cparam
      use EquationOfState, only: gamma, gamma1, cs2bot, cs2top, get_cp1

      implicit none

      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(mx,mz) :: finter,source
      real, dimension(nx)    :: a,b,c,wx
      real, dimension(nz)    :: rhs,work
      real    :: alpha, aalpha, bbeta, cp1, dx_2, dz_2
!
      source=(f(:,4,:,ilnTT)-finit(:,4,:,ilnTT))/dt
      call get_cp1(cp1)
      dx_2=1./dx**2
      dz_2=1./dz**2
!
!  row dealt implicitly
!
      do j=n1,n2
        wx=dt*gamma*hcond0*cp1/exp(f(l1:l2,4,j,ilnrho))
        a=-wx*dx_2/2.
        b=1.+wx*dx_2
        c=a
!
        rhs=finit(l1:l2,4,j,ilnTT)+wx*dz_2/2.*                    &
            (finit(l1:l2,4,j+1,ilnTT)-2.*finit(l1:l2,4,j,ilnTT)+  &
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
!  columns dealt implicitly
!
      do i=l1,l2
        wx=dt*gamma*hcond0*cp1/exp(f(i,4,n1:n2,ilnrho))
        a=-wx*dz_2/2.
        b=1.+wx*dz_2
        c=a
!
        rhs=finter(i,n1:n2)+wx*dx_2/2.*                              &
           (finter(i+1,n1:n2)-2.*finter(i,n1:n2)+finter(i-1,n1:n2))  &
           +dt*source(i,n1:n2)/2.
!
        c(nz)=0.
        a(1)=0.
        b(1)=1.
        b(nz)=1.
        c(1)=0.
        a(nz)=0.
        rhs(1)=cs2bot/gamma1
        rhs(nx)=cs2top/gamma1
        call tridag(a,b,c,rhs,work,nz)
        f(i,4,n1:n2,ilnTT)=work(1:nz)
      enddo
!
      call BC_CT(f(:,4,:,ilnTT))
!
    end subroutine ADI_Kconst
!**************************************************************
    subroutine ADI_Kprof(finit,f)
!
!  10-Sep-07/gastine+dintrans: coded
!  2-D ADI scheme for the radiative diffusion term where the radiative
!  conductivity depends on T (uses heatcond_TT to compute hcond _and_
!  dhcond).
!
      use Cdata
      use Cparam
      use EquationOfState, only: gamma, gamma1, cs2bot, cs2top, get_cp1

      implicit none

      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(mx,mz) :: source,hcond,dhcond,finter,val
      real, dimension(nx)    :: a,b,c, wx
      real, dimension(nz)    :: rhs,work
      real    :: alpha, aalpha, bbeta
      real    :: dx_2, dz_2, cp1

      source=(f(:,4,:,ilnTT)-finit(:,4,:,ilnTT))/dt
      call heatcond_TT(finit(:,4,:,ilnTT),hcond,dhcond)
      call get_cp1(cp1)
      dx_2=1./dx**2
      dz_2=1./dz**2
!
!  lignes en implicite
!
      do j=n1,n2
! a=-dt/2*J_x for i=i-1 (lower diagonal)
       wx=cp1*gamma/exp(f(l1:l2,4,j,ilnrho))
       a=-dt*wx*dx_2/4.*(dhcond(l1-1:l2-1,j)     &
         *(finit(l1-1:l2-1,4,j,ilnTT)-finit(l1:l2,4,j,ilnTT))   &
         +hcond(l1-1:l2-1,j)+hcond(l1:l2,j))
! b=1-dt/2*J_x for i=i (main diagonal)
       b=1.+dt*wx*dx_2/4.*(dhcond(l1:l2,j)       &
         *(2.*finit(l1:l2,4,j,ilnTT)-finit(l1-1:l2-1,4,j,ilnTT) &
         -finit(l1+1:l2+1,4,j,ilnTT))+2.*hcond(l1:l2,j)         &
         +hcond(l1+1:l2+1,j)+hcond(l1-1:l2-1,j))
! c=-dt/2*J_x for i=i+1 (upper diagonal)
       c=-dt*wx*dx_2/4.*(dhcond(l1+1:l2+1,j)     &
          *(finit(l1+1:l2+1,4,j,ilnTT)-finit(l1:l2,4,j,ilnTT))  &
          +hcond(l1:l2,j)+hcond(l1+1:l2+1,j))
! rhs=f_y(T^n) + f_x(T^n) (Eq. 3.6)
! do first f_y(T^n)
       rhs=wx*dz_2/2.*((hcond(l1:l2,j+1)      &
           +hcond(l1:l2,j))*(finit(l1:l2,4,j+1,ilnTT)           &
           -finit(l1:l2,4,j,ilnTT))-(hcond(l1:l2,j)             &
           +hcond(l1:l2,j-1))                                   &
           *(finit(l1:l2,4,j,ilnTT)-finit(l1:l2,4,j-1,ilnTT)))
! then add f_x(T^n)
       rhs=rhs+wx*dx_2/2.*((hcond(l1+1:l2+1,j) &
         +hcond(l1:l2,j))*(finit(l1+1:l2+1,4,j,ilnTT)-finit(l1:l2,4,j,ilnTT))&
           -(hcond(l1:l2,j)+hcond(l1-1:l2-1,j)) &
           *(finit(l1:l2,4,j,ilnTT)-finit(l1-1:l2-1,4,j,ilnTT)))
!
       aalpha=c(nx)
       bbeta=a(1)
       c(nx)=0.
       a(1)=0.
       call cyclic(a,b,c,aalpha,bbeta,rhs,work,nx)
       finter(l1:l2,j)=work(1:nx)
      enddo
!
!  colonnes en implicite
!
      do i=l1,l2
       wx=dt*cp1*gamma*dz_2/exp(f(i,4,n1:n2,ilnrho))
       a=-wx/4.*(dhcond(i,n1-1:n2-1) &
         *(finit(i,4,n1-1:n2-1,ilnTT)-finit(i,4,n1:n2,ilnTT))&
         +hcond(i,n1-1:n2-1)+hcond(i,n1:n2))
!
       b=1.+wx/4.*(dhcond(i,n1:n2)* &
         (2.*finit(i,4,n1:n2,ilnTT)-finit(i,4,n1-1:n2-1,ilnTT) &
         -finit(i,4,n1+1:n2+1,ilnTT))+2.*hcond(i,n1:n2) &
         +hcond(i,n1+1:n2+1)+hcond(i,n1-1:n2-1))
!
       c=-wx/4.*(dhcond(i,n1+1:n2+1) &
         *(finit(i,4,n1+1:n2+1,ilnTT)-finit(i,4,n1:n2,ilnTT))&
         +hcond(i,n1:n2)+hcond(i,n1+1:n2+1))
!
       rhs=finter(i,n1:n2)
!
       c(nz)=0.
       a(1)=0.
! Constant flux at the bottom
!       b(1)=-1.
!       c(1)=0.
!       c(1)=1.
!       rhs(1)=0.
! Constant temperature at the bottom
        b(1)=1.
        c(1)=0.
        rhs(1)=0.
! Constant temperature at the top
       b(nx)=1.
       a(nx)=0.
       rhs(nz)=0.
!
       call tridag(a,b,c,rhs,work,nz)
       val(i,n1:n2)=work(1:nz)
      enddo
!
      f(:,4,:,ilnTT)=finit(:,4,:,ilnTT)+dt*(val+source)
!
! Boris: is it really needed? i.e. filling f(n1,n2) and heatcond_TT
      f(:,:,n1,ilnTT)=cs2bot/gamma1
      f(:,:,n2,ilnTT)=cs2top/gamma1

      call BC_CT(f(:,4,:,ilnTT))
!     call BC_flux(f,hcond)
      call heatcond_TT(f(:,4,:,ilnTT),hcond,dhcond)
!
    end subroutine ADI_Kprof
!**************************************************************
    subroutine heatcond_TT_2d(TT,hcond,dhcond)
!
! 07-Sep-07/gastine: computed the radiative conductivity hcond(T) and
! its derivative dhcond=dhcond(T)/dT when the input TT is a 2d-array.
!
      implicit none

      real, dimension(mx,mz) :: TT, arg, hcond, dhcond
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
! constant hcond for testing purpose
!
!     hcond=hcond0 ; dhcond=0.
!
    end subroutine heatcond_TT_2d
!**************************************************************
    subroutine heatcond_TT_point(TT,hcond,dhcond)
!
! 07-Sep-07/gastine: computed the radiative conductivity hcond(T) and
! its derivative dhcond=dhcond(T)/dT when the input TT is a single
! point.
!
      implicit none

      real :: TT, arg, hcond, dhcond
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
! constant hcond for testing purpose
!
!     hcond=hcond0 ; dhcond=0.
!
    end subroutine heatcond_TT_point
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

    end subroutine BC_CT
!**************************************************************
    subroutine BC_flux(f,hcond)

      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz) :: hcond
      integer :: i

      do i=1,nghost
        f(:,:,n1-i,ilnTT)=f(:,:,n1+i,ilnTT)+2*i*dz*Fbot/hcond(10,n1+i)
      enddo
      f(:,:,n2+1,ilnTT)=2*f(:,:,n2,ilnTT)-f(:,:,n2-1,ilnTT)
      ! x-direction
      f(1:l1-1,:,:,ilnTT)=f(l2i:l2,:,:,ilnTT)
      f(l2+1:mx,:,:,ilnTT)=f(l1:l1i,:,:,ilnTT)

    end subroutine BC_flux
!**************************************************************
    subroutine tridag(a,b,c,r,u,n)

      integer :: j,n
      integer, parameter :: NMAX=500
      real    :: bet
      real, dimension(n) :: a,b,c,r,u
      real, dimension(NMAX) :: gam
!
      if(b(1).eq.0.) pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
!
      return
    end subroutine tridag
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
      if(n.le.2)pause 'n too small in cyclic'
      if(n.gt.NMAX)pause 'NMAX too small in cyclic'
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
    end subroutine cyclic
!***************************************************************
endmodule Entropy
