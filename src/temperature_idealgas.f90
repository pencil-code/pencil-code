! $Id$
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
! PENCILS PROVIDED Ma2; uglnTT; ugTT; fpres(3)
!
!***************************************************************
module Entropy

  use Cdata
  use Cparam
  use EquationOfState, only: mpoly0, mpoly1
  use Interstellar
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include 'entropy.h'

  interface heatcond_TT ! Overload subroutine `hcond_TT' function
    module procedure heatcond_TT_0d  ! get one value (hcond, dhcond)
    module procedure heatcond_TT_1d  ! get 1d-arrays (hcond, dhcond)
    module procedure heatcond_TT_2d  ! get 2d-arrays (hcond, dhcond)
  end interface

  real :: radius_lnTT=0.1,ampl_lnTT=0.,widthlnTT=2*epsi
  real :: lnTT_left=1.0,lnTT_right=1.0,lnTT_const=0.0,TT_const=1.0
  real :: kx_lnTT=1.0,ky_lnTT=1.0,kz_lnTT=1.0
  real :: Kgperp=0.,Kgpara=0.
  real :: chi=impossible,heat_uniform=0.0,difflnTT_hyper=0.
  real :: zbot=0.0,ztop=0.0
  real :: tau_heat_cor=-1.0,tau_damp_cor=-1.0,zcor=0.0,TT_cor=0.0
  real :: center1_x=0., center1_y=0., center1_z=0.
  real :: r_bcz=0.,chi_shock=0.,chi_hyper3=0.
! entries for ADI
  real :: Tbump=0.,Kmin=0.,Kmax=0.,hole_slope=0.,hole_width=0.
  real :: hole_alpha ! initialized _after_ the reading
  integer, parameter :: nheatc_max=2
  logical :: lpressuregradient_gas=.true.,ladvection_temperature=.true.
  logical :: lupw_lnTT=.false.,lcalc_heat_cool=.false.,lheatc_hyper3=.false.
  logical :: lheatc_Kconst=.false.,lheatc_Kprof=.false.,lheatc_Karctan=.false.
  logical :: lheatc_tensordiffusion=.false.
  logical :: lheatc_chiconst=.false.,lheatc_chiconst_accurate=.false.
  logical :: lfreeze_lnTTint=.false.,lfreeze_lnTText=.false.
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  logical :: lhcond_global=.false.
  logical :: lheatc_shock=.false.,lheatc_hyper3_polar=.false.
  logical :: lviscosity_heat=.true.
  integer :: iglobal_hcond=0
  integer :: iglobal_glhc=0
  logical :: linitial_log=.false.

  character (len=labellen), dimension(ninit) :: initlnTT='nothing'
  character (len=5) :: iinit_str

! Delete (or use) me asap!
  real :: hcond0=impossible, hcond1=1.
  real :: Fbot=impossible,FbotKbot,Ftop,Kbot,FtopKtop
  logical :: lmultilayer=.false.
!
! input parameters
  namelist /entropy_init_pars/ &
      initlnTT,radius_lnTT,ampl_lnTT,widthlnTT, &
      lnTT_left,lnTT_right,lnTT_const,TT_const, &
      kx_lnTT,ky_lnTT,kz_lnTT,center1_x,center1_y,center1_z, &
      mpoly0,mpoly1,r_bcz, &
      Fbot,Tbump,Kmin,Kmax,hole_slope,hole_width,ltemperature_nolog,&
      linitial_log,hcond0
!
! run parameters
  namelist /entropy_run_pars/ &
      lupw_lnTT,lpressuregradient_gas,ladvection_temperature, &
      heat_uniform,chi,iheatcond,tau_heat_cor,tau_damp_cor,zcor,TT_cor, &
      lheatc_chiconst_accurate,hcond0,lcalc_heat_cool,&
      lfreeze_lnTTint,lfreeze_lnTText,widthlnTT,mpoly0,mpoly1, &
      lhcond_global,lviscosity_heat,chi_hyper3,chi_shock, &
      Fbot,Tbump,Kmin,Kmax,hole_slope,hole_width,Kgpara,Kgperp
!
! other variables (needs to be consistent with reset list below)
  integer :: idiag_TTmax=0    ! DIAG_DOC: $\max (T)$
  integer :: idiag_gTmax=0    ! DIAG_DOC: $\max (|grad T|)$
  integer :: idiag_TTmin=0    ! DIAG_DOC: $\min (T)$
  integer :: idiag_TTm=0      ! DIAG_DOC: $\left< T \right>$
  integer :: idiag_fradtop=0  ! DIAG_DOC: $<-K{dT\over dz}>_{\text{top}}$ 
                              ! DIAG_DOC: \quad(radiative flux at the top)
  integer :: idiag_yHmax=0,idiag_yHmin=0,idiag_yHm=0
  integer :: idiag_ethm=0,idiag_ssm=0,idiag_thcool=0
  integer :: idiag_eem=0,idiag_ppm=0,idiag_csm=0
  integer :: idiag_dtc=0        ! DIAG_DOC: $\delta t/[c_{\delta t}\,\delta_x
                                ! DIAG_DOC:   /\max c_{\rm s}]$
                                ! DIAG_DOC:   \quad(time step relative to
                                ! DIAG_DOC:   acoustic time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_dtchi=0      ! DIAG_DOC: $\delta t / [c_{\delta t,{\rm v}}\,
                                ! DIAG_DOC:   \delta x^2/\chi_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to time
                                ! DIAG_DOC:   step based on heat conductivity;
                                ! DIAG_DOC:   see \S~\ref{time-step})

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
!  Writing files for use with IDL.
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
      use Gravity, only: g0, gravz
      use EquationOfState, only : cs2bot, cs2top, gamma, gamma_m1, &
                                  select_eos_variable
      use Sub, only: step,der_step
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: hcond, dhcond
      logical :: lstarting, lnothing
      integer :: i, ierr
!
!  Set iTT requal to ilnTT if we are considering non-logarithmic temperature.
!
      if (ltemperature_nolog) iTT=ilnTT
!
      if (.not. leos) then
         call fatal_error('initialize_entropy','EOS=noeos but temperature_idealgas requires an EQUATION OF STATE for the fluid')
      endif
!
      if (ltemperature_nolog) then
        call select_eos_variable('TT',iTT)
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
      lheatc_tensordiffusion=.false.
      lheatc_chiconst = .false.
!
!  Initialize thermal diffusion
!
      lheatc_shock=.false.
      lheatc_hyper3=.false.
      lheatc_hyper3_polar=.false.
!
      lnothing = .false.
!
      do i=1,nheatc_max
      select case (iheatcond(i))
        case ('K-const')
          lheatc_Kconst=.true.
          if (lroot) call information('initialize_entropy', &
          ' heat conduction: K=cst --> gamma*K/rho/TT/cp*div(T*grad lnTT)')
        case ('K-profile')
          lheatc_Kprof=.true.
! 
!  11-Aug-2008/dintrans: better somewhere else?
!
          hcond1=(mpoly1+1.)/(mpoly0+1.)
          Fbot=-gamma/(gamma-1.)*hcond0*g0/(mpoly0+1.)
          if (lroot) call information('initialize_entropy',' heat conduction: K=K(r)')
        case ('K-arctan')
          lheatc_Karctan=.true.         
          if (.not. ltemperature_nolog) &
            call fatal_error('calc_heatcond_arctan','only valid for TT')
          if (lroot) call information('initialize_entropy',' heat conduction: arctan profile')
        case ('chi-const')
          lheatc_chiconst=.true.
          if (lroot) call information('initialize_entropy',' heat conduction: constant chi')
        case ('chi-hyper3')
          lheatc_hyper3=.true.
          if (lroot) call information('initialize_entropy','hyper conductivity')
        case ('hyper3_cyl','hyper3-cyl','hyper3_sph','hyper3-sph')
          lheatc_hyper3_polar=.true.
          if (lroot) call information('initialize_entropy','hyper conductivity: polar coords')
        case ('shock','chi-shock')
          lheatc_shock=.true.
          if (lroot) call information('initialize_entropy','shock conductivity')
        case ('tensor-diffusion')
          lheatc_tensordiffusion=.true.
          if (lroot) print*, 'heat conduction: tensor diffusion'
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
! WL: is this correct? In density, lnothing=.true. is inside the do-loop
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
      if (initlnTT(1).eq.'gaussian') then
! needed when one only works with temperature_idealgas to check the
! radiative diffusion term, i.e. one solves d(TT)/dt=gamma*chi*del2(TT)
! with bcz='cT' (all other modules are down)
        cs2bot=gamma_m1*f(l1,4,n1,ilnTT)
        cs2top=gamma_m1*f(l1,4,n2,ilnTT)
      endif
!
! some tricks regarding Fbot and hcond0 when bcz1='c1' (constant flux)
!
      if (hole_slope.eq.0.) then
       if (bcz1(ilnTT)=='c1' .and. lrun) then
        if (Fbot==impossible .and. hcond0 /= impossible) then
          Fbot=-gamma/gamma_m1*hcond0*gravz/(mpoly0+1.)
          if (lroot) print*, &
                     'initialize_entropy: Calculated Fbot = ', Fbot
        endif
        if (hcond0==impossible .and. Fbot /= impossible) then
          hcond0=-Fbot*gamma_m1/gamma*(mpoly0+1.)/gravz
          if (lroot) print*, &
                     'initialize_entropy: Calculated hcond0 = ', hcond0
        endif
        if (Fbot==impossible .and. hcond0==impossible) &
          call fatal_error('temperature_idealgas',  &
                           'Both Fbot and hcond0 are unknown')
       endif
      endif
!
!  30-nov-2007/dintrans: now hcond0 and Fbot are passed to boundcond()
!  as shared variables
!
      if (initlnTT(1).eq.'rad_equil') then
        ! use hcondADI to dynamically share with boundcond() for the 'c3' BC
        if (lrun) then
          call heatcond_TT_1d(f(:,4,n1,ilnTT), hcondADI)
        else
          hcondADI=spread(Kmax, 1, mx)
        endif
      else
        call put_shared_variable('hcond0', hcond0, ierr)
        if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting hcond0")
      endif
      call put_shared_variable('Fbot', Fbot, ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting Fbot")
!
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
      if (lrun) then
        if (lheatc_hyper3 .and. chi_hyper3==0.0) &
            call fatal_error('initialize_entropy', &
            'Conductivity coefficient chi_hyper3 is zero!')
        if (lheatc_shock .and. chi_shock==0.0) &
            call fatal_error('initialize_entropy', &
            'Conductivity coefficient chi_shock is zero!')
      endif

      if (iheatcond(1)=='nothing') then
        if (hcond0 /= impossible) call warning('initialize_entropy', 'No heat conduction, but hcond0 /= 0')
        if (chi /= impossible) call warning('initialize_entropy', 'No heat conduction, but chi /= 0')
      endif

      call put_shared_variable('lviscosity_heat',lviscosity_heat,ierr)
      if (ierr/=0) call stop_it("initialize_entropy in temperature_idealgas: "//&
           "there was a problem when putting lviscosity_heat")

      call keep_compiler_quiet(lstarting)

    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f)
!
!  initialise lnTT or TT (; called from start.f90
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
      use InitialCondition, only: initial_condition_ss
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, cs20, &
                                 lnrho0, get_cp1
      use Gravity, only: gravz
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
!
      integer :: j
      logical :: lnothing=.true.
      real :: haut, Rgas, cp1, Ttop, alpha, beta, expo
!
      do j=1,ninit
!
        if (initlnTT(j)/='nothing') then
!
          lnothing=.false.

          call chn(j,iinit_str)
!
!  select different initial conditions
!
          select case (initlnTT(j))
          case ('zero', '0'); f(:,:,:,ilnTT) = 0.
!
          case ('const_lnTT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+lnTT_const
!
          case ('const_TT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+log(TT_const)
!
          case ('single_polytrope'); call single_polytrope(f)
!
          case ('gaussian')
            do n=n1,n2
              f(l1:l2,4,n,ilnTT)=exp(-(x(l1:l2)/radius_lnTT)**2)* &
                  exp(-((z(n)-0.5)/radius_lnTT)**2)
            enddo
            cs2bot=gamma_m1*f(l1,4,n1,ilnTT)
            cs2top=gamma_m1*f(l1,4,n2,ilnTT)
!
          case ('rad_equil')
            call rad_equil(f)
!
          case ('blob_hs')
            if (lroot) print*, 'init_lnTT: hydrostatic blob with ', &
                radius_lnTT, ampl_lnTT, center1_x, center1_y, center1_z
            call blob(ampl_lnTT,f,ilnTT,radius_lnTT,center1_x,center1_y,center1_z)
            call blob(-ampl_lnTT,f,ilnrho,radius_lnTT,center1_x,center1_y,center1_z)
!
          case ('isothermal')
            if (lroot) print*, 'init_lnTT: isothermal atmosphere'
            if (ltemperature_nolog) then
              f(:,:,:,iTT)  =cs20/gamma_m1
            else
              f(:,:,:,ilnTT)=log(cs20/gamma_m1)
            endif
            haut=-cs20/gamma/gravz
            do n=n1,n2
              f(:,:,n,ilnrho)=lnrho0+(1.-z(n))/haut
            enddo
!
          case ('hydro_rad')
            if (lroot) print*, 'init_lnTT: hydrostatic+radiative equilibria'
            if (Fbot==impossible .or. hcond0==impossible) &
                call stop_it("initialize_lnTT: Fbot or hcond0 not initialized")
            call get_cp1(cp1)
            Rgas=(1.-1./gamma)/cp1
            Ttop=cs20/gamma_m1
            beta=-Fbot/hcond0
            alpha=Ttop-beta
            expo=-gravz/beta/Rgas
            do n=n1,n2
              if (ltemperature_nolog) then
                f(:,:,n,iTT)  =beta*z(n)+alpha
              else
                f(:,:,n,ilnTT)=log(beta*z(n)+alpha)
              endif
              f(:,:,n,ilnrho)=lnrho0+ &
                  (1.+expo)*log((1.+alpha/beta)/(z(n)+alpha/beta))
            enddo
          !
          case default
          !
          !  Catch unknown values
          !
            write(unit=errormsg,fmt=*) 'No such value for initss(' &
                           //trim(iinit_str)//'): ',trim(initlnTT(j))
            call fatal_error('init_ss',errormsg)

          endselect

          if (lroot) print*,'init_ss: initss(' &
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
      if (ltemperature_nolog.and.linitial_log) f(:,:,:,iTT)=exp(f(:,:,:,ilnTT))
!
    endsubroutine init_ss
!***********************************************************************
    subroutine pencil_criteria_entropy()
!
      use Cdata
!
      if (ldt) lpenc_requested(i_cs2)=.true.
!
      if (lpressuregradient_gas) lpenc_requested(i_fpres)=.true.
!
      if (lviscosity.and.lviscosity_heat) then
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_visc_heat)=.true.
        if (.not.ltemperature_nolog) &
            lpenc_requested(i_TT1)=.true.
      endif
!
      if (ldensity) then
        lpenc_requested(i_divu)=.true.
        if (ltemperature_nolog) lpenc_requested(i_TT)=.true.
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
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (lheatc_Kconst) then
        if (ldensity) lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cp1)=.true.
        if (ltemperature_nolog) then
          lpenc_requested(i_del2TT)=.true.
        else
          lpenc_requested(i_glnTT)=.true.
          lpenc_requested(i_del2lnTT)=.true.
        endif
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
        lpenc_requested(i_gTT)=.true.
        lpenc_requested(i_del2TT)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (lheatc_shock) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_del2lnrho)=.true.
        lpenc_requested(i_gshock)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
!
      if (lheatc_tensordiffusion) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (lheatc_hyper3) lpenc_requested(i_del6lnTT)=.true.
!
      if (ladvection_temperature) then 
        if (ltemperature_nolog) then 
          lpenc_requested(i_ugTT)=.true.
        else
          lpenc_requested(i_uglnTT)=.true.
        endif
      endif
!
!  Diagnostics
!
      if (idiag_TTmax/=0) lpenc_diagnos(i_TT)  =.true.
      if (idiag_gTmax/=0) then
         lpenc_diagnos(i_glnTT) =.true.
         lpenc_diagnos(i_TT) =.true.
      endif
      if (idiag_TTmin/=0) lpenc_diagnos(i_TT)  =.true.
      if (idiag_TTm/=0)   lpenc_diagnos(i_TT)  =.true.
      if (idiag_fradtop/=0) then
        lpenc_diagnos(i_TT) =.true.  ! for hcond computation
        lpenc_diagnos(i_glnTT) =.true.
      endif
      if (idiag_yHmax/=0) lpenc_diagnos(i_yH)  =.true.
      if (idiag_yHmin/=0) lpenc_diagnos(i_yH)  =.true.
      if (idiag_yHm/=0)   lpenc_diagnos(i_yH)  =.true.
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
!
      if (lpencil_in(i_uglnTT)) lpencil_in(i_glnTT)=.true.
      if (lpencil_in(i_ugTT))   lpencil_in(i_gTT)  =.true.
!
      if (lpencil_in(i_fpres)) then
        lpencil_in(i_cs2)=.true.
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
      if (lpencil(i_ugTT)) then
        call u_dot_grad(f,ilnTT,p%gTT,p%uu,p%ugTT,UPWIND=lupw_lnTT)
      endif
!
! fpres
!
      if (lpencil(i_fpres)) then
        do j=1,3
          p%fpres(:,j)=-p%cs2*(p%glnrho(:,j) + p%glnTT(:,j))*gamma_inv
        enddo
      endif
!
    endsubroutine calc_pencils_entropy
!***********************************************************************
    subroutine dss_dt(f,df,p)
!
!  Calculate right hand side of temperature equation
!  heat condution is currently disabled until old stuff,
!  which in now in calc_heatcond, has been reinstalled.
!  lnTT version: DlnTT/Dt = -gamma_m1*divu + gamma*cp1*rho1*TT1*RHS
!    TT version:   DTT/Dt = -gamma_m1*TT*divu + gamma*cp1*rho1*RHS
!
!  13-dec-02/axel+tobi: adapted from entropy
!
      use Cdata
      use Deriv, only: der6
      use Diagnostics
      use EquationOfState, only: gamma_m1,gamma
      use Mpicomm
      use Special, only: special_calc_entropy
      use Sub
      use Viscosity, only: calc_viscous_heat
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: Hmax=0., hcond, thdiff=0.,tmp
      real :: fradtop
      integer :: j
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
      if (ladvection_temperature) then
        if (ltemperature_nolog) then 
          df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   - p%ugTT
        else
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%uglnTT
        endif
      endif
!
!  Calculate viscous contribution to temperature
!
      if (lviscosity.and.lviscosity_heat) call calc_viscous_heat(f,df,p,Hmax)
!
!  Various heating conduction contributions
!
      if (lcalc_heat_cool)  call calc_heat_cool(df,p)
!
!  Thermal conduction: only chi=cte for the moment
!
      if (lheatc_chiconst) call calc_heatcond_constchi(df,p)
      if (lheatc_Kconst)   call calc_heatcond_constK(df,p)
      if (lheatc_Kprof)    call calc_heatcond(f,df,p)
      if (lheatc_Karctan)  call calc_heatcond_arctan(df,p)
      if (lheatc_tensordiffusion) call calc_heatcond_tensor(df,p)
!
!  Hyper diffusion
!
      if (lheatc_hyper3) then
        if (ltemperature_nolog) then 
          thdiff=thdiff+chi_hyper3*p%del6TT
        else
          thdiff=thdiff+chi_hyper3*p%del6lnTT
        endif
        if (lfirst.and.ldt) diffus_chi3=diffus_chi3+chi_hyper3*dxyz_6
        if (headtt) print*,'dss_dt: chi_hyper3=', chi_hyper3
      endif
!
      if (lheatc_hyper3_polar) then
        do j=1,3
          call der6(f,ilnTT,tmp,j,IGNOREDX=.true.)
          if (.not.ltemperature_nolog) tmp=tmp*p%TT1
          thdiff = thdiff + chi_hyper3*pi4_1*tmp*dline_1(:,j)**2
        enddo
        if (lfirst.and.ldt) &
             diffus_chi3=diffus_chi3+chi_hyper3*pi4_1*dxyz_2
        if (headtt) print*,'dss_dt: chi_hyper3=', chi_hyper3
      endif
!
! Shock diffusion 
!
      if (lheatc_shock) call calc_heatcond_shock(df,p)
!
!  Entry possibility for "personal" entries.
!  In that case you'd need to provide your own "special" routine.
!
      if (lspecial) call special_calc_entropy(f,df,p)
!
!  Add thermal diffusion to temperature equation
!
      if (ltemperature_nolog) then
        df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   + thdiff
      else
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + thdiff
      endif
!
!  Information on the timescales
!
      if (lfirst.and.ldt.and.(headtt.or.ldebug)) then
          print*,'dss_dt: max(diffus_chi ) =', maxval(diffus_chi)
          print*,'dss_dt: max(diffus_chi3) =', maxval(diffus_chi3)
      endif
!
!  Need to add left-hand-side of the continuity equation (see manual)
!  Check this
!
      if (ldensity) then
        if (ltemperature_nolog) then
          df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   - gamma_m1*p%TT*p%divu
        else
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - gamma_m1*p%divu
        endif
      endif
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        if (idiag_TTmax/=0) call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_gTmax/=0) then
           call dot2(p%glnTT,tmp)
           call max_mn_name(p%TT*sqrt(tmp),idiag_gTmax)
        endif
        if (idiag_TTmin/=0) call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        if (idiag_TTm/=0)   call sum_mn_name(p%TT,idiag_TTm)
        if (idiag_fradtop/=0.and.n==n2) then
          call heatcond_TT(p%TT, hcond)
          fradtop=sum(-hcond*p%glnTT(:,3))/nx
          call save_name(fradtop, idiag_fradtop)
        endif
        if (idiag_ethm/=0)   call sum_mn_name(p%ee/p%rho1,idiag_ethm)
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
    subroutine calc_heatcond_shock(df,p)
!
!  Adds in shock diffusion to the temperature equation.
!  Ds/Dt = ... + 1/(rho*T) grad(flux), where
!  flux = chi_shock*rho*T*grads
!  (in comments we say chi_shock, but in the code this is "chi_shock*p%shock")
!
!  01-aug-08/wlad: adapted from entropy
!
      use Cdata
      use Diagnostics
      use EquationOfState, only: gamma
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: thdiff,g2,gshockglnTT
!
      intent(in) :: p
      intent(out) :: df
!
      if (headtt) print*,'calc_heatcond_shock: chi_shock=',chi_shock
!
!  calculate terms for shock diffusion
!
      call dot(p%gshock,p%glnTT,gshockglnTT)
      call dot(p%glnrho+p%glnTT,p%glnTT,g2)
!
!  shock entropy diffusivity
!
      if (headtt) print*,'calc_heatcond_shock: chi_shock=',chi_shock
      if (ltemperature_nolog) then
!  WL: can someone explain this equation?      
        thdiff=gamma*chi_shock*(p%shock*(p%del2lnrho+g2)+gshockglnTT)
      else 
!  D(lnTT)/Dt = ... + 1/pp*grad(shock*pp*glnTT)      
        thdiff=chi_shock*(p%shock*(p%del2lnTT+g2)+gshockglnTT)
      endif
!
      df(l1:l2,m,n,ilntt) = df(l1:l2,m,n,ilntt) + thdiff
      if (headtt) print*,'calc_heatcond_shock: added thdiff'
!
      if (lfirst.and.ldt) then
        if (leos_idealgas) then
          diffus_chi=diffus_chi+(gamma*chi_shock*p%shock)*dxyz_2
        else
          diffus_chi=diffus_chi+(chi_shock*p%shock)*dxyz_2
        endif
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif
! 
    endsubroutine calc_heatcond_shock
!***********************************************************************
    subroutine rad_equil(f)
!
! 16-mai-07/gastine+dintrans: compute the radiative and hydrostatic 
! equilibria for a given radiative profile defined in heatcond_TT.
!
      use Cdata
      use Gravity, only: gravz
      use EquationOfState, only: lnrho0,cs20,cs2top,cs2bot,gamma, &
                                 gamma_m1,eoscalc,ilnrho_TT
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mz) :: temp,lnrho
      real :: hcond, dtemp, dlnrho, ss
      integer :: i
!
      if (.not. ltemperature_nolog) &
        call fatal_error('temperature_idealgas',  &
                         'rad_equil not implemented for lnTT')
      if (lroot) print*,'init_ss: rad_equil for kappa-mechanism pb'
!
! Integrate from top to bottom: z(n2) --> z(n1)
!
      temp(n2)=cs20/gamma_m1
      lnrho(n2)=lnrho0
      f(:,:,n2,ilnTT)=cs20/gamma_m1
      f(:,:,n2,ilnrho)=lnrho0
!
! Calculate the n2-1 gridpoint thanks to a 1st order forward Euler scheme
!
      call heatcond_TT_0d(temp(n2), hcond)
      dtemp=Fbot/hcond
      temp(n2-1)=temp(n2)+dz*dtemp
      dlnrho=(-gamma/gamma_m1*gravz-dtemp)/temp(n2)
      lnrho(n2-1)=lnrho(n2)+dz*dlnrho
      f(:,:,n2-1,ilnTT)=temp(n2-1)
      f(:,:,n2-1,ilnrho)=lnrho(n2-1)
!
! Now we use a 2nd order centered scheme for the other gridpoints
!
      do i=n2-1,n1+1,-1
        call heatcond_TT(temp(i), hcond)
        dtemp=Fbot/hcond
        temp(i-1)=temp(i+1)+2.*dz*dtemp
        dlnrho=(-gamma/gamma_m1*gravz-dtemp)/temp(i)
        lnrho(i-1)=lnrho(i+1)+2.*dz*dlnrho
        f(:,:,i-1,ilnTT)=temp(i-1)
        f(:,:,i-1,ilnrho)=lnrho(i-1)
      enddo
!
! Initialize cs2bot by taking into account the new bottom value of temperature
! Note: cs2top=cs20 already defined in eos_idealgas
      cs2bot=gamma_m1*temp(n1)
      print*,'cs2top, cs2bot=', cs2top, cs2bot
!
      if (lroot) then
        print*,'--> write the initial setup in data/proc0/setup.dat'
        open(unit=11,file=trim(directory)//'/setup.dat')
        write(11,'(5a14)') 'z','rho','temp','ss','hcond'
        do i=n2,n1,-1
          call eoscalc(ilnrho_TT,lnrho(i),temp(i),ss=ss)
          call heatcond_TT(temp(i), hcond)
          write(11,'(5e14.5)') z(i),exp(lnrho(i)),temp(i),ss,hcond
        enddo
        close(11)
      endif
!
    endsubroutine rad_equil
!***********************************************************************
    subroutine calc_heat_cool(df,p)
!
      use Diagnostics
      use EquationOfState, only: gamma,gamma_m1
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: tau,cooling,kappa,a1,a3
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
!  Calculate the radiative diffusion term for chi=cte:
!  lnTT version: cp*chi*Div(rho*T*glnTT)/(rho*cv*TT)
!           = gamma*chi*(g2.glnTT+g2lnTT) where g2=glnrho+glnTT
!    TT version: cp*chi*Div(rho*gTT)/(rho*cv)
!           = gamma*chi*(g2.gTT+g2TT) where g2=glnrho
!
      use Diagnostics
      use EquationOfState, only: gamma
      use Sub

      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: g2
!
      if (ltemperature_nolog) then
        call dot(p%glnrho,p%glnTT,g2)
      else
        call dot(p%glnTT+p%glnrho,p%glnTT,g2)
      endif
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
!  Calculate the radiative diffusion term for K=cte:
!  lnTT version: gamma*K/rho/TT/cp*div(T*grad lnTT)
!                =gamma*K/rho/cp*(gradlnTT.gradlnTT + del2ln TT)
!    TT version: gamma*K/rho/cp*del2(TT)=gamma*chi*del2(TT)
!  Note: if ldensity=.false. then rho=1 and chi=K/cp
!
      use Diagnostics
      use EquationOfState, only: gamma
      use Sub

      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case)  :: p
      real, dimension(nx) :: g2, chix
!
!  Add heat conduction to RHS of temperature equation
!
      if (ldensity) then
        chix=p%rho1*hcond0*p%cp1
      else
        chix=hcond0*p%cp1
      endif
!
      if (ltemperature_nolog) then
        df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT)   + gamma*chix*p%del2TT
      else
        call dot(p%glnTT,p%glnTT,g2)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*(g2 + p%del2lnTT)
      endif
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
! 16-mai-07/gastine+dintrans: radiative diffusion with an arctan
!  profile for the conductivity
!  calculate gamma/rho*cp*div(K *grad TT)=
!    gamma*K/rho*cp*(grad LnK.grad TT + del2 TT)
!
      use Diagnostics
      use EquationOfState, only: gamma
      use Sub

      real, dimension(mx,my,mz,mvar) :: df
      real, dimension (nx)   :: hcond, dhcond, g1, chix
      real, dimension (nx,3) :: gLnhcond=0.
      type (pencil_case)     :: p
!
      call heatcond_TT(p%TT, hcond, dhcond)
! must specify the new bottom value of hcond for the 'c1' BC
      if (n == n1) hcond0=hcond(1) 
! grad LnK=grad_T Ln K.grad(TT)
      dhcond=dhcond/hcond
      call multsv(dhcond, p%gTT, gLnhcond)
      call dot(gLnhcond, p%gTT, g1)
!
!  Add heat conduction to RHS of temperature equation
!
      chix=p%rho1*hcond*p%cp1
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + gamma*chix*(g1+p%del2TT)
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
!  12-Mar-07/dintrans: coded
!  calculate gamma*K/rho*cp*div(T*grad lnTT)= 
!              gamma*K/rho*cp*(gradlnTT.gradln(hcond*TT) + del2ln TT)
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
    subroutine calc_heatcond_tensor(df,p)
!
!  Calculates heat conduction parallel and perpendicular (isotropic)
!  to magnetic field lines
!  
!  25-aug-09/bing: moved from dss_dt to here
!
      use Diagnostics
      use EquationOfState, only: gamma
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: cosbgT,gT2,b2
      real, dimension (nx) :: vKpara,vKperp,rhs
!      
      vKpara(:) = Kgpara
      vKperp(:) = Kgperp
!
      call tensor_diffusion_coef(p%glnTT,p%hlnTT,p%bij,p%bb,vKperp,vKpara,rhs,llog=.true.)
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+rhs*p%rho1*gamma*p%cp1
!
      call dot(p%bb,p%glnTT,cosbgT)
      call dot2(p%glnTT,gT2)
      call dot2(p%bb,b2)
!
      where ((gT2.le.tini).or.(b2.le.tini))
         cosbgT=0.
      elsewhere 
         cosbgT=cosbgT/sqrt(gT2*b2)
      endwhere
!
      if (lfirst.and.ldt) then
         diffus_chi=diffus_chi+cosbgT*gamma*Kgpara*p%rho1*p%cp1*dxyz_2
         if (ldiagnos.and.idiag_dtchi/=0) then
            call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
         endif
      endif
!
    endsubroutine calc_heatcond_tensor
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
      use Diagnostics
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
        idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0; idiag_fradtop=0
        idiag_yHmax=0; idiag_yHmin=0; idiag_yHm=0; idiag_gTmax=0
        idiag_ethm=0; idiag_ssm=0; idiag_thcool=0
        idiag_dtchi=0; idiag_dtc=0
        idiag_eem=0; idiag_ppm=0; idiag_csm=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'gTmax',idiag_gTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'fradtop',idiag_fradtop)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
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
        write(3,*) 'i_gTmax=',idiag_gTmax
        write(3,*) 'i_TTmin=',idiag_TTmin
        write(3,*) 'i_TTm=',idiag_TTm
        write(3,*) 'i_fradtop=',idiag_fradtop
        write(3,*) 'i_ethm=',idiag_ethm
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
    subroutine get_slices_entropy(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))

!
!  Temperature.
!
        case ('tt') 
          slices%yz =f(ix_loc,m1:m2,n1:n2,ilnTT)
          slices%xz =f(l1:l2,iy_loc,n1:n2,ilnTT)
          slices%xy =f(l1:l2,m1:m2,iz_loc,ilnTT)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,ilnTT)
          if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,ilnTT)
          if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,ilnTT)
          slices%ready=.true.

      endselect
!
    endsubroutine get_slices_entropy
!***********************************************************************
    subroutine single_polytrope(f)
!
! 04-aug-07/dintrans: a single polytrope with index mpoly0
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
!  gamma*(Rgas/mu)T0 = cs2(ad) = cp*T0*gamma_m1,
!  so T0 = cs20*cp1/gamma_m1
!
      call get_cp1(cp1)
      beta=-cp1*gravz/(mpoly0+1.)*gamma/gamma_m1
      ztop=xyz0(3)+Lxyz(3)
      zbot=xyz0(3)
      T0=cs20*cp1/gamma_m1
      print*, 'polytrope: mpoly0, beta, T0=', mpoly0, beta, T0
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        temp=T0+beta*(ztop-z(n))
        if (ltemperature_nolog) then
          f(:,m,n,iTT)  =temp
        else
          f(:,m,n,ilnTT)=log(temp)
        endif
        f(:,m,n,ilnrho)=lnrho0+mpoly0*log(temp/T0)
      enddo
      cs2bot=gamma_m1*(T0+beta*(ztop-zbot))
      cs2top=cs20
!
    endsubroutine single_polytrope
!***********************************************************************
    subroutine calc_heatcond_ADI(finit,f)
!
!  10-sep-07/gastine+dintrans: wrapper to the two possible ADI subroutines
!  ADI_Kconst: constant radiative conductivity
!  ADI_Kprof: radiative conductivity depends on T, i.e. hcond(T)
!
      use Cparam

      implicit none

      real, dimension(mx,my,mz,mfarray) :: finit, f
!
      if (hcond0 /= impossible) then
        if (nx == 1) then
          call ADI_Kconst_1d(finit,f)
        else
          call ADI_Kconst(finit,f)
!          call toto(finit,f)
        endif
      else
        if (nx == 1) then
          call ADI_Kprof_1d(finit,f)
        else
          call ADI_Kprof(finit,f)
        endif
      endif
!
    endsubroutine calc_heatcond_ADI
!***********************************************************************
    subroutine ADI_Kconst(finit,f)
!
!  08-Sep-07/gastine+dintrans: coded
!  2-D ADI scheme for the radiative diffusion term (see
!  Peaceman & Rachford 1955). Each direction is solved implicitly:
!
!    (1-dt/2*Lambda_x)*T^(n+1/2) = (1+dt/2*Lambda_y)*T^n + source/2
!    (1-dt/2*Lambda_y)*T^(n+1)   = (1+dt/2*Lambda_x)*T^(n+1/2) + source/2
!
!  where Lambda_x and Lambda_y denote diffusion operators and the source
!  term comes from the explicit advance.
!
      use Cdata
      use Cparam
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, get_cp1
      use General, only: tridag
      use Gravity, only: gravz

      implicit none

      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(nx,nz) :: finter, source, rho, TT
      real, dimension(nx)    :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nz)    :: az, bz, cz, wz, rhsz, workz
      real    :: aalpha, bbeta, cp1, dx_2, dz_2, Fbot, tmp_flux
!
      TT=finit(l1:l2,4,n1:n2,ilnTT)
      source=(f(l1:l2,4,n1:n2,ilnTT)-TT)/dt
      if (ldensity) then
        rho=exp(f(l1:l2,4,n1:n2,ilnrho))
      else
        rho=1.
      endif
      call get_cp1(cp1)
      dx_2=1./dx**2
      dz_2=1./dz**2
!
!  row dealt implicitly
!
      do j=1,nz
        wx=dt*gamma*hcond0*cp1/rho(:,j)
        ax=-wx*dx_2/2.
        bx=1.+wx*dx_2
        cx=ax
!
        if (j==1) then
          if (bcz1(ilnTT)=='cT') then
            !constant temperature: T(j-1)=2*T(j)-T(j+1)
            rhsx=TT(:,j)+dt/2.*source(:,j)
          else
            !imposed flux: T(j-1)=T(j+1)-2*dz*tmp_flux with tmp_flux=-Fbot/hcond0
            Fbot=-gamma/(gamma-1.)*hcond0*gravz/(mpoly0+1.)
            tmp_flux=-Fbot/hcond0
            rhsx=TT(:,j)+wx*dz_2/2.*                            &
                (TT(:,j+1)-2.*TT(:,j)+TT(:,j+1)-2.*dz*tmp_flux) &
                +dt/2.*source(:,j)
          endif
        elseif (j==nz) then
          !constant temperature: T(j+1)=2*T(j)-T(j-1)
          rhsx=TT(:,j)+dt/2.*source(:,j)
        else
          rhsx=TT(:,j)+wx*dz_2/2.*                         &
              (TT(:,j+1)-2.*TT(:,j)+TT(:,j-1))             &
               +dt/2.*source(:,j)
        endif
!
! x boundary conditions: periodic
!
        aalpha=cx(nx) ; bbeta=ax(1)
        call cyclic(ax, bx, cx, aalpha, bbeta, rhsx, workx, nx)
        finter(:,j)=workx
      enddo
!
!  columns dealt implicitly
!
      do i=1,nx
        wz=dt*gamma*hcond0*cp1/rho(i,:)
        az=-wz*dz_2/2.
        bz=1.+wz*dz_2
        cz=az
!
        if (i==1) then
          rhsz=finter(i,:)+wz*dx_2/2.*                       &
              (finter(i+1,:)-2.*finter(i,:)+finter(nx,:))    &
              +dt/2.*source(i,:)
        elseif (i==nx) then
          rhsz=finter(i,:)+wz*dx_2/2.*                       &
              (finter(1,:)-2.*finter(i,:)+finter(i-1,:))     &
              +dt/2.*source(i,:)
        else
          rhsz=finter(i,:)+wz*dx_2/2.*                       &
              (finter(i+1,:)-2.*finter(i,:)+finter(i-1,:))   &
              +dt/2.*source(i,:)
        endif
        !
        ! z boundary conditions
        ! Always constant temperature at the top
        !
        bz(nz)=1. ; az(nz)=0.
        rhsz(nz)=cs2top/gamma_m1
        select case (bcz1(ilnTT))
          ! Constant temperature at the bottom
          case ('cT')
            bz(1)=1.  ; cz(1)=0. 
            rhsz(1)=cs2bot/gamma_m1
          ! Constant flux at the bottom
          case ('c1')
!           bz(1)=1.   ; cz(1)=-1
!           rhsz(1)=dz*Fbot/hcond0
! we can use here the second-order relation for the first derivative: 
! (T_{j+1}-T_{j_1})/2dz = dT/dz --> T_{j-1} = T_{j+1} - 2*dz*dT/dz 
! and insert this expression in the difference relation to eliminate T_{j-1}:
! a_{j-1}*T_{j-1} + b_j T_j + c_{j+1}*T_{j+1} = RHS
            cz(1)=cz(1)+az(1)
            rhsz(1)=rhsz(1)-2.*az(1)*dz*Fbot/hcond0
          case default 
           call fatal_error('ADI_Kconst','bcz on TT must be cT or c1')
        endselect
!
        call tridag(az, bz, cz, rhsz, workz)
        f(i+nghost,4,n1:n2,ilnTT)=workz
      enddo
!
    endsubroutine ADI_Kconst
!***********************************************************************
    subroutine ADI_Kprof(finit,f)
!
!  10-Sep-07/gastine+dintrans: coded
!  2-D ADI scheme for the radiative diffusion term where the radiative
!  conductivity depends on T (uses heatcond_TT to compute hcond _and_
!  dhcond). The ADI scheme is of Yakonov's form:
!
!    (1-dt/2*J_x)*lambda = f_x(T^n) + f_y(T^n) + source
!    (1-dt/2*J_y)*beta   = lambda
!    T^(n+1) = T^n + dt*beta
!
!    where J_x and J_y denote Jacobian matrices df/dT.
!
      use Cdata
      use Cparam
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, get_cp1
      use General, only: tridag

      implicit none

      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: finit, f
      real, dimension(mx,mz) :: source, hcond, dhcond, finter, val, TT, rho
      real, dimension(nx)    :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nz)    :: az, bz, cz, wz, rhsz, workz
      real    :: aalpha, bbeta
      real    :: dx_2, dz_2, cp1

      source=(f(:,4,:,ilnTT)-finit(:,4,:,ilnTT))/dt
      call get_cp1(cp1)
      dx_2=1./dx**2
      dz_2=1./dz**2
! BC important not for the x-direction (always periodic) but for 
! the z-direction as we must impose the 'c3' BC at the 2nd-order
! before going in the implicit stuff
      call heatcond_TT_2d(finit(:,4,:,ilnTT), hcond, dhcond)
      call boundary_ADI(finit(:,4,:,ilnTT), hcond(:,n1))
      TT=finit(:,4,:,ilnTT)
      if (ldensity) then
        rho=exp(f(:,4,:,ilnrho))
      else
        rho=1.
      endif
!
!  rows dealt implicitly
!
      do j=n1,n2
       wx=cp1*gamma/rho(l1:l2,j)
! ax=-dt/2*J_x for i=i-1 (lower diagonal)
       ax=-dt*wx*dx_2/4.*(dhcond(l1-1:l2-1,j)    &
         *(TT(l1-1:l2-1,j)-TT(l1:l2,j))          &
         +hcond(l1-1:l2-1,j)+hcond(l1:l2,j))
! bx=1-dt/2*J_x for i=i (main diagonal)
       bx=1.+dt*wx*dx_2/4.*(dhcond(l1:l2,j)      &
         *(2.*TT(l1:l2,j)-TT(l1-1:l2-1,j)        &
         -TT(l1+1:l2+1,j))+2.*hcond(l1:l2,j)     &
         +hcond(l1+1:l2+1,j)+hcond(l1-1:l2-1,j))
! cx=-dt/2*J_x for i=i+1 (upper diagonal)
       cx=-dt*wx*dx_2/4.*(dhcond(l1+1:l2+1,j)    &
          *(TT(l1+1:l2+1,j)-TT(l1:l2,j))         &
          +hcond(l1:l2,j)+hcond(l1+1:l2+1,j))
! rhsx=f_y(T^n) + f_x(T^n) (Eq. 3.6)
! do first f_y(T^n)
       rhsx=wx*dz_2/2.*((hcond(l1:l2,j+1)        &
           +hcond(l1:l2,j))*(TT(l1:l2,j+1)       &
           -TT(l1:l2,j))-(hcond(l1:l2,j)         &
           +hcond(l1:l2,j-1))                    &
           *(TT(l1:l2,j)-TT(l1:l2,j-1)))
! then add f_x(T^n)
       rhsx=rhsx+wx*dx_2/2.*((hcond(l1+1:l2+1,j)         &
         +hcond(l1:l2,j))*(TT(l1+1:l2+1,j)-TT(l1:l2,j))  &
           -(hcond(l1:l2,j)+hcond(l1-1:l2-1,j))          &
           *(TT(l1:l2,j)-TT(l1-1:l2-1,j)))+source(l1:l2,j)
!
! x boundary conditions: periodic
       aalpha=cx(nx) ; bbeta=ax(1)
       call cyclic(ax,bx,cx,aalpha,bbeta,rhsx,workx,nx)
       finter(l1:l2,j)=workx(1:nx)
      enddo
!
!  columns dealt implicitly
!
      do i=l1,l2
       wz=dt*cp1*gamma*dz_2/rho(i,n1:n2)
       az=-wz/4.*(dhcond(i,n1-1:n2-1)   &
         *(TT(i,n1-1:n2-1)-TT(i,n1:n2)) &
         +hcond(i,n1-1:n2-1)+hcond(i,n1:n2))
!
       bz=1.+wz/4.*(dhcond(i,n1:n2)*             &
         (2.*TT(i,n1:n2)-TT(i,n1-1:n2-1)         &
         -TT(i,n1+1:n2+1))+2.*hcond(i,n1:n2)     &
         +hcond(i,n1+1:n2+1)+hcond(i,n1-1:n2-1))
!
       cz=-wz/4.*(dhcond(i,n1+1:n2+1)            &
         *(TT(i,n1+1:n2+1)-TT(i,n1:n2))          &
         +hcond(i,n1:n2)+hcond(i,n1+1:n2+1))
!
       rhsz=finter(i,n1:n2)
!
! z boundary conditions
! Constant temperature at the top: T^(n+1)-T^n=0
       bz(nz)=1. ; az(nz)=0.
       rhsz(nz)=0.
! bottom
       select case (bcz1(ilnTT))
! Constant temperature at the bottom: T^(n+1)-T^n=0
         case ('cT')
          bz(1)=1. ; cz(1)=0.
          rhsz(1)=0.
! Constant flux at the bottom
         case ('c3')
          bz(1)=1. ; cz(1)=-1.
          rhsz(1)=0.
         case default 
          call fatal_error('ADI_Kprof','bcz on TT must be cT or c3')
       endselect
!
       call tridag(az,bz,cz,rhsz,workz)
       val(i,n1:n2)=workz(1:nz)
      enddo
!
      f(:,4,:,ilnTT)=finit(:,4,:,ilnTT)+dt*val
!
! update hcond used for the 'c3' condition in boundcond.f90
!
      call heatcond_TT_1d(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof
!***********************************************************************
    subroutine heatcond_TT_2d(TT, hcond, dhcond)
!
! 07-Sep-07/gastine: computed 2-D radiative conductivity hcond(T) with
! its derivative dhcond=dhcond(T)/dT.
!
      implicit none

      real, dimension(mx,mz) :: TT, arg, hcond
      real, dimension(mx,mz), optional :: dhcond
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
    endsubroutine heatcond_TT_2d
!***********************************************************************
    subroutine heatcond_TT_1d(TT, hcond, dhcond)
!
! 18-Sep-07/dintrans: computed 1-D radiative conductivity 
! hcond(T) with its derivative dhcond=dhcond(T)/dT.
!
      implicit none

      real, dimension(:)           :: TT, hcond
      real, dimension(:), optional :: dhcond
      real, dimension(size(TT,1))  :: arg
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
    endsubroutine heatcond_TT_1d
!***********************************************************************
    subroutine heatcond_TT_0d(TT, hcond, dhcond)
!
! 07-Sep-07/gastine: computed the radiative conductivity hcond(T) 
! with its derivative dhcond=dhcond(T)/dT at a given temperature.
!
      implicit none

      real :: TT, arg, hcond
      real, optional :: dhcond
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
    endsubroutine heatcond_TT_0d
!***********************************************************************
    subroutine boundary_ADI(f_2d, hcond)

! 13-Sep-07/gastine: computed two different types of boundary 
! conditions for the implicit solver:
!     - Always periodic in x-direction
!     - Possibility to choose between 'cT' and 'c3' in z direction
! Note: 'c3' means that the flux is constant at the *bottom only*
!
      implicit none

      real, dimension(mx,mz) :: f_2d
      real, dimension(mx), optional :: hcond

! x-direction: periodic
      f_2d(1:l1-1,:)=f_2d(l2i:l2,:)
      f_2d(l2+1:mx,:)=f_2d(l1:l1i,:)
!
! top boundary condition z=z(n2): always constant temperature
!
      f_2d(:,n2+1)=2.*f_2d(:,n2)-f_2d(:,n2-1)
!
! bottom bondary condition z=z(n1): constant T or imposed flux dT/dz
!
      select case (bcz1(ilnTT))
        case ('cT') ! constant temperature
          f_2d(:,n1-1)=2.*f_2d(:,n1)-f_2d(:,n1+1)
        case ('c3') ! constant flux
          if (.not. present(hcond)) then
            f_2d(:,n1-1)=f_2d(:,n1+1)+2.*dz*Fbot/hcond0
          else 
            f_2d(:,n1-1)=f_2d(:,n1+1)+2.*dz*Fbot/hcond(:)
          endif
      endselect

    endsubroutine boundary_ADI
!***********************************************************************
    subroutine cyclic(a,b,c,alpha,beta,r,x,n)
!
! 08-Sep-07/gastine+dintrans: coded
! inversion of a tridiagonal matrix with periodic BC (alpha and beta
! coefficients); used in the ADI scheme
!
      use General, only: tridag
!
      implicit none
!
      integer :: i,n
      integer, parameter    :: NMAX=600
      real    :: alpha, beta,gamma,fact      
      real, dimension(n)    :: a,b,c,r,x,bb,u,z
!     real, dimension(NMAX) :: bb,u,z
!
      if (n <= 2) call fatal_error('cyclic', 'n too small in cyclic')
      if (n > NMAX) call fatal_error('cyclic', 'NMAX too small in cyclic')
      gamma=-b(1)
      bb(1)=b(1)-gamma
      bb(n)=b(n)-alpha*beta/gamma
      do i=2,n-1
        bb(i)=b(i)
      enddo
      call tridag(a,bb,c,r,x)
      u(1)=gamma
      u(n)=alpha
      do i=2,n-1
        u(i)=0.
      enddo
      call tridag(a,bb,c,u,z)
      fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
      do i=1,n
        x(i)=x(i)-fact*z(i)
      enddo
!
      return
    endsubroutine cyclic
!***********************************************************************
    subroutine ADI_Kconst_1d(finit,f)
!
! 18-sep-07/dintrans: coded
! Implicit Crank Nicolson scheme in 1-D for a constant K (not 
! really an ADI but keep the generic name for commodity).
!
      use Cdata
      use Cparam
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, get_cp1
      use General, only: tridag

      implicit none

      integer :: j, jj
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(mz) :: rho, TT
      real, dimension(nz) :: a, b, c, rhs, work
      real  :: cp1, dz_2, wz
!
!     source=(f(4,4,:,ilnTT)-finit(4,4,:,ilnTT))/dt
!     TT=finit(4,4,:,ilnTT)
      TT=f(4,4,:,ilnTT)
      rho=exp(f(4,4,:,ilnrho))
      call get_cp1(cp1)
      dz_2=1./dz**2
!
      do j=n1,n2
        wz=dt*gamma*hcond0*cp1/rho(j)
        jj=j-nghost
        a(jj)=-wz*dz_2/2.
        b(jj)=1.+wz*dz_2
        c(jj)=a(jj)
!
!       rhs(jj)=TT(j)+wz*dz_2/2.*(TT(j+1)-2.*TT(j)+TT(j-1))+dt*source(j)
        rhs(jj)=TT(j)+wz*dz_2/2.*(TT(j+1)-2.*TT(j)+TT(j-1))
      enddo
! apply the boundary conditions *outside* the j-loop
! Always constant temperature at the top
      b(nz)=1. ; a(nz)=0.
      rhs(nz)=cs2top/gamma_m1
      if (bcz1(ilnTT)=='cT') then
! Constant temperature at the bottom
        b(1)=1. ; c(1)=0. 
        rhs(1)=cs2bot/gamma_m1
      else
! Constant flux at the bottom
        b(1)=1.  ; c(1)=-1.
        rhs(1)=dz*Fbot/hcond0
      endif
      call tridag(a, b, c, rhs, work)
      f(4,4,n1:n2,ilnTT)=work
!
      call keep_compiler_quiet(finit)
!
    endsubroutine ADI_Kconst_1d
!***********************************************************************
    subroutine ADI_Kprof_1d(finit,f)
!
! 18-sep-07/dintrans: coded
! Implicit 1-D case for a temperature-dependent conductivity K(T).
! Not really an ADI but keep the generic name for commodity.
!
      use Cdata
      use Cparam
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, get_cp1
      use General, only: tridag

      implicit none

      integer :: j, jj
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(mz) :: source, rho, TT, hcond, dhcond
      real, dimension(nz) :: a, b, c, rhs, work
      real  :: cp1, dz_2, wz, hcondp, hcondm
!
      source=(f(4,4,:,ilnTT)-finit(4,4,:,ilnTT))/dt
      call get_cp1(cp1)
      dz_2=1./dz**2
      rho=exp(f(4,4,:,ilnrho))
! need to set up the 'c3' BC at the 2nd-order before the implicit stuff
      call heatcond_TT_1d(finit(4,4,:,ilnTT), hcond, dhcond)
      hcondADI=spread(hcond(1), 1, mx)
      call boundary_ADI(finit(:,4,:,ilnTT), hcondADI)
      TT=finit(4,4,:,ilnTT)
!
      do j=n1,n2
        jj=j-nghost
        wz=dt*dz_2*gamma*cp1/rho(j)
        hcondp=hcond(j+1)+hcond(j)
        hcondm=hcond(j)+hcond(j-1)
!
        a(jj)=-wz/4.*(hcondm-dhcond(j-1)*(TT(j)-TT(j-1)))
        b(jj)=1.-wz/4.*(-hcondp-hcondm+dhcond(j)*(TT(j+1)-2.*TT(j)+TT(j-1)))
        c(jj)=-wz/4.*(hcondp+dhcond(j+1)*(TT(j+1)-TT(j)))
!
        rhs(jj)=wz/2.*(hcondp*(TT(j+1)-TT(j))-hcondm*(TT(j)-TT(j-1))) &
          +dt*source(j)
!
! Always constant temperature at the top: T^(n+1)-T^n = 0
        b(nz)=1. ; a(nz)=0.
        rhs(nz)=0.
        if (bcz1(ilnTT)=='cT') then
! Constant temperature at the bottom
          b(1)=1. ; c(1)=0. 
          rhs(1)=0.
        else
! Constant flux at the bottom: d/dz [T^(n+1)-T^n] = 0
          b(1)=1.  ; c(1)=-1.
          rhs(1)=0.
        endif
      enddo
      call tridag(a,b,c,rhs,work)
      f(4,4,n1:n2,ilnTT)=work+TT(n1:n2)
!
! Update the bottom value of hcond used for the 'c3' BC in boundcond
!
      call heatcond_TT_1d(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof_1d
!***********************************************************************
    subroutine toto(finit,f)
!
!  04-sep-2009/dintrans: coded
!  parallel version of the ADI scheme for the K=cte case
!
      use Cdata
      use Cparam
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, get_cp1
      use General, only: tridag
      use Gravity, only: gravz
      use Mpicomm, only: transp_xz, transp_zx, mpisend_real, mpirecv_real

      implicit none

      integer :: i,j,tagr1,tagr2,tags1,tags2
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(nx,nz)  :: finter, source, rho, TT
      real, dimension(nx)     :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nzgrid) :: az, bz, cz, wz, rhsz, workz
      real, dimension(nzgrid,nx/nprocz) :: fintert, rhot, sourcet, wtmp
      real, dimension(nzgrid) :: tmp1, tmp2
      real, dimension(nxgrid) :: tmpx1, tmpx2
      real  :: aalpha, bbeta, cp1, dx_2, dz_2, Fbot
!
      TT=finit(l1:l2,4,n1:n2,ilnTT)
      source=(f(l1:l2,4,n1:n2,ilnTT)-TT)/dt
      if (ldensity) then
        rho=exp(f(l1:l2,4,n1:n2,ilnrho))
      else
        rho=1.
      endif
      call get_cp1(cp1)
      dx_2=1./dx**2
      dz_2=1./dz**2
!
!  row dealt implicitly
!
      tagr1 = 200 + 10*(iproc)   + iproc - 1
      tagr2 = 200 + 10*(iproc)   + iproc
      tags1 = 200 + 10*(iproc-1) + iproc - 2
      tags2 = 200 + 10*(iproc+1) + iproc + 1
      if (iproc==0) then
        call mpisend_real(TT(:,nz), nxgrid, iproc+1, tags2)
        call mpirecv_real(tmpx1, nxgrid, iproc+1, tagr1)
        !constant temperature: T(j-1)=2*T(j)-T(j+1)
        tmpx2 = 2*TT(:,1)-TT(:,2)
      elseif (iproc==nprocz-1) then
        call mpirecv_real(tmpx2, nxgrid, iproc-1, tagr2)
        call mpisend_real(TT(:,1), nxgrid, iproc-1, tags1)
        !constant temperature: T(j+1)=2*T(j)-T(j-1)
        tmpx1 = 2*TT(:,nz)-TT(:,nz-1)
      else
        call mpisend_real(TT(:,1), nxgrid, iproc-1, tags1)
        call mpisend_real(TT(:,nz), nxgrid, iproc+1, tags2)
        call mpirecv_real(tmpx1, nxgrid, iproc+1, tagr1)
        call mpirecv_real(tmpx2, nxgrid, iproc-1, tagr2)
      endif
!
      do j=1,nz
        wx=dt*gamma*hcond0*cp1/rho(:,j)
        ax=-wx*dx_2/2.
        bx=1.+wx*dx_2
        cx=ax
!
        if (j==1) then
          rhsx=TT(:,j)+wx*dz_2/2.*                         &
              (TT(:,j+1)-2.*TT(:,j)+tmpx2)                  &
               +dt/2.*source(:,j)
        elseif (j==nz) then
          rhsx=TT(:,j)+wx*dz_2/2.*                         &
              (tmpx1-2.*TT(:,j)+TT(:,j-1))                  &
               +dt/2.*source(:,j)
        else
          rhsx=TT(:,j)+wx*dz_2/2.*                         &
              (TT(:,j+1)-2.*TT(:,j)+TT(:,j-1))             &
               +dt/2.*source(:,j)
        endif
!
! x boundary conditions: periodic
!
        aalpha=cx(nx) ; bbeta=ax(1)
        call cyclic(ax, bx, cx, aalpha, bbeta, rhsx, workx, nx)
        finter(:,j)=workx
      enddo

! do the transpositions z <--> x
      call transp_xz(finter, fintert)
      call transp_xz(rho, rhot)
      call transp_xz(source, sourcet)

! communicate gridpoints in the x-direction for periodic BC
      tagr1 = 400 + 10*(iproc)   + iproc - 1
      tagr2 = 400 + 10*(iproc)   + iproc
      tags1 = 400 + 10*(iproc-1) + iproc - 2
      tags2 = 400 + 10*(iproc+1) + iproc + 1
      if (iproc==0) then
        call mpisend_real(fintert(:,1), nzgrid, nprocz-1, 111)
        call mpisend_real(fintert(:,nx/nprocz), nzgrid, iproc+1, tags2)
        call mpirecv_real(tmp1, nzgrid, iproc+1, tagr1)
        call mpirecv_real(tmp2, nzgrid, nprocz-1, 112)
      elseif (iproc==nprocz-1) then
        call mpirecv_real(tmp1, nzgrid, 0, 111)
        call mpirecv_real(tmp2, nzgrid, iproc-1, tagr2)
        call mpisend_real(fintert(:,1), nzgrid, iproc-1, tags1)
        call mpisend_real(fintert(:,nx/nprocz), nzgrid, 0, 112)
      else
        call mpisend_real(fintert(:,1), nzgrid, iproc-1, tags1)
        call mpisend_real(fintert(:,nx/nprocz), nzgrid, iproc+1, tags2)
        call mpirecv_real(tmp1, nzgrid, iproc+1, tagr1)
        call mpirecv_real(tmp2, nzgrid, iproc-1, tagr2)
      endif

!
!  columns dealt implicitly
!
      do i=1,nx/nprocz
        wz=dt*gamma*hcond0*cp1/rhot(:,i)
        az=-wz*dz_2/2.
        bz=1.+wz*dz_2
        cz=az
!
        if (i==1) then
          rhsz=fintert(:,i)+wz*dx_2/2.*               &
!             (fintert(:,i+1)-2.*fintert(:,i)+fintert(:,nx/nprocz))   &
              (fintert(:,i+1)-2.*fintert(:,i)+tmp2)   &
              +dt/2.*sourcet(:,i)
        elseif (i==nx/nprocz) then
          rhsz=fintert(:,i)+wz*dx_2/2.*               &
!             (fintert(:,1)-2.*fintert(:,i)+fintert(:,i-1))   &
              (tmp1-2.*fintert(:,i)+fintert(:,i-1))   &
              +dt/2.*sourcet(:,i)
        else
          rhsz=fintert(:,i)+wz*dx_2/2.*                        &
              (fintert(:,i+1)-2.*fintert(:,i)+fintert(:,i-1))  &
              +dt/2.*sourcet(:,i)
        endif
        !
        ! z boundary conditions
        ! Always constant temperature at the top
        !
        bz(nzgrid)=1. ; az(nzgrid)=0.
        rhsz(nzgrid)=cs2top/gamma_m1
        select case (bcz1(ilnTT))
          ! Constant temperature at the bottom
          case ('cT')
            bz(1)=1.  ; cz(1)=0. 
            rhsz(1)=cs2bot/gamma_m1
          ! Constant flux at the bottom
          case ('c1')
!           bz(1)=1.   ; cz(1)=-1
!           rhsz(1)=dz*Fbot/hcond0
! we can use here the second-order relation for the first derivative: 
! (T_{j+1}-T_{j_1})/2dz = dT/dz --> T_{j-1} = T_{j+1} - 2*dz*dT/dz 
! and insert this expression in the difference relation to eliminate T_{j-1}:
! a_{j-1}*T_{j-1} + b_j T_j + c_{j+1}*T_{j+1} = RHS
            cz(1)=cz(1)+az(1)
            rhsz(1)=rhsz(1)-2.*az(1)*dz*Fbot/hcond0
          case default 
           call fatal_error('ADI_Kconst','bcz on TT must be cT or c1')
        endselect
!
        call tridag(az, bz, cz, rhsz, workz)
        wtmp(:,i)=workz
      enddo
      call transp_zx(wtmp, f(l1:l2,4,n1:n2,ilnTT))
!
    endsubroutine toto
!***********************************************************************
    subroutine dss_dt_after_mn(f,df,p)
!
! Do nothing.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(out) :: df
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dss_dt_after_mn
!***********************************************************************
    subroutine calc_pencils_entropy_after_mn(f,p)
!
! Do nothing 
! DM+PC

      use EquationOfState, only: gamma,gamma_m1,cs20,lnrho0,profz_eos
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_entropy_after_mn
!**********************************************************************
endmodule Entropy
