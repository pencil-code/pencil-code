! $Id$

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
! PENCILS PROVIDED Ma2; uglnTT; cvspec(nchemspec)
!
!***************************************************************
module Entropy

  use Cparam
  use Cdata
  use Messages
  use Interstellar

  implicit none

  include 'entropy.h'

  real :: radius_lnTT=0.1,ampl_lnTT=0.,widthlnTT=2*epsi
  real :: lnTT_left=1.0,lnTT_right=1.0,lnTT_const=0.0,TT_const=1.0
  real :: kx_lnTT=1.0,ky_lnTT=1.0,kz_lnTT=1.0
  real :: chi=0.0,heat_uniform=0.0,chi_hyper3=0.0
  real :: zbot=0.0,ztop=0.0
  real :: tau_heat_cor=-1.0,tau_damp_cor=-1.0,zcor=0.0,TT_cor=0.0
  logical :: lpressuregradient_gas=.true.,ladvection_temperature=.true.
  logical :: lupw_lnTT=.false.,lcalc_heat_cool=.false.
  logical :: lheatc_chiconst=.false.,lheatc_chiconst_accurate=.false.
  logical :: lheatc_hyper3=.false.
  logical :: lviscous_heat=.true.
  character (len=labellen), dimension(ninit) :: initlnTT='nothing'
  character (len=5) :: iinit_str

  ! Delete (or use) me asap!
  real :: hcond0,hcond1,Fbot,FbotKbot,Ftop,Kbot,FtopKtop
  logical :: lmultilayer

  ! input parameters
  namelist /entropy_init_pars/ &
      initlnTT,radius_lnTT,ampl_lnTT,widthlnTT, &
      lnTT_left,lnTT_right,lnTT_const,TT_const, &
      kx_lnTT,ky_lnTT,kz_lnTT

  ! run parameters
  namelist /entropy_run_pars/ &
      lupw_lnTT,lpressuregradient_gas,ladvection_temperature, &
      heat_uniform,chi,tau_heat_cor,tau_damp_cor,zcor,TT_cor, &
      lheatc_chiconst_accurate,lheatc_hyper3,chi_hyper3, &
      lviscous_heat

  ! other variables (needs to be consistent with reset list below)
    integer :: idiag_TTmax=0,idiag_TTmin=0,idiag_TTm=0
    integer :: idiag_yHmax=0,idiag_yHmin=0,idiag_yHm=0
    integer :: idiag_eth=0,idiag_ssm=0,idiag_cv=0,idiag_cp=0
    integer :: idiag_dtchi=0,idiag_dtc=0
    integer :: idiag_eem=0,idiag_ppm=0,idiag_csm=0
    integer :: idiag_mum=0

  contains

!***********************************************************************
    subroutine register_entropy()
!
!  initialise variables which should know that we solve an entropy
!  equation: ilnTT, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
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
           "$Id$")
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
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Check any module dependencies
!
      if (.not.leos_temperature_ionization) then
        if (.not.leos_chemistry)  then
        call fatal_error('initialize_entropy','EOS/=noeos but'//&
                         'temperature_ionization already include'//&
                         'an EQUATION OF STATE for the fluid')
        endif
      endif
!
!  Check whether we want heating/cooling
!
      lcalc_heat_cool = (heat_uniform/=0.0.or.tau_heat_cor>0)
!
!  Define bottom and top z positions
!  (TH: This should really be global variables IMHO)
!
      zbot = xyz0(3)
      ztop = xyz0(3) + Lxyz(3)
!
!  Check whether we want heat conduction
!
      lheatc_chiconst = (chi > tiny(chi))
      lheatc_hyper3 = (chi_hyper3 > tiny(chi_hyper3))

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
    subroutine init_ss(f)
!
!  initialise entropy; called from start.f90
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!
      use General, only: chn
      use Sub, only: blob
      use Initcond, only: jump
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
!
      integer :: j
      logical :: lnothing=.true.
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
      select case(initlnTT(j))

        case('zero', '0'); f(:,:,:,ilnTT) = 0.
        case('const_lnTT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+lnTT_const
        case('const_TT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+log(TT_const)
        case('blob'); call blob(ampl_lnTT,f,ilnTT,radius_lnTT,0.,0.,0.)
        case('xwave')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnTT)=f(l1:l2,m,n,ilnTT)+ampl_lnTT*sin(kx_lnTT*x(l1:l2))
          enddo; enddo
        case('ywave')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnTT)=f(l1:l2,m,n,ilnTT)+ampl_lnTT*sin(ky_lnTT*y(m))
          enddo; enddo
        case('zwave')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ilnTT)=f(l1:l2,m,n,ilnTT)+ampl_lnTT*sin(kz_lnTT*z(n))
          enddo; enddo
        case('xjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'x')
        case('yjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'y')
        case('zjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'z')
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

      if (lnothing.and.lroot) print*,'init_ss: nothing'
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
        lpenc_requested(i_gamma1)=.true.
        lpenc_requested(i_delta)=.true.
        lpenc_requested(i_divu)=.true.
      endif

      if (linterstellar) then
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_ee)=.true.
      endif

      if (lpressuregradient_gas) then
        lpenc_requested(i_rho1gpp)=.true.
      endif

      if (lviscosity) then
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_visc_heat)=.true.
      endif

      if (lcalc_heat_cool) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_TT1)=.true.
        if (tau_heat_cor>0) lpenc_requested(i_TT)=.true.
      endif

      if (ladvection_temperature) lpenc_requested(i_uglnTT)=.true.

      if (lheatc_chiconst) then
        lpenc_requested(i_del2lnTT)=.true.
        if (lheatc_chiconst_accurate) then
          lpenc_requested(i_cp1)=.true.
          lpenc_requested(i_gradcp)=.true.
          lpenc_requested(i_gamma)=.true.
        endif
      endif

      if (lheatc_hyper3) lpenc_requested(i_del6lnTT)=.true.

     ! if (lchemistry) then
     !   lpenc_requested(i_DYDt_reac)=.true.
     !   lpenc_requested(i_DYDt_diff)=.true.
     !   if (lheatc_chemistry) then
     !      lpenc_requested(i_lambda)=.true.
     !      lpenc_requested(i_glnlambda)=.true.
     !   endif
     ! endif

!
!  Diagnostics
!
      if (idiag_TTmax/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_TTmin/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_TTm/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_yHmax/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_yHmin/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_yHm/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_eth/=0) then
        lpenc_diagnos(i_rho1)=.true.
        lpenc_diagnos(i_ee)=.true.
      endif
      if (idiag_ssm/=0) lpenc_diagnos(i_ss)=.true.
      if (idiag_cv/=0) lpenc_diagnos(i_cv)=.true.
      if (idiag_cp/=0) lpenc_diagnos(i_cp)=.true.
      if (idiag_dtchi/=0) then
        lpenc_diagnos(i_rho1)=.true.
        lpenc_diagnos(i_cv1)=.true.
      endif
      if (idiag_dtchi/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_csm/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_eem/=0) lpenc_diagnos(i_ee)=.true.
      if (idiag_ppm/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_mum/=0) lpenc_diagnos(i_mu1)=.true.

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

    endsubroutine pencil_interdep_entropy
!***********************************************************************
    subroutine calc_pencils_entropy(f,p)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use Sub, only: u_dot_grad

      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      type (pencil_case), intent (inout) :: p

!
!  Mach Speed
!
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2

!
!  Temperature advection
!  (Needs to be here because of lupw_lnTT)
!
      if (lpencil(i_uglnTT)) then
        call u_dot_grad(f,ilnTT,p%glnTT,p%uu,p%uglnTT,UPWIND=lupw_lnTT)
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
      use Sub, only: max_mn_name,sum_mn_name,identify_bcs,cubic_step
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (mx,my,mz,mvar), intent (out) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: damp
      real, dimension (nx) :: Hmax, sum_DYDt
      real :: prof,react_rate
      integer :: j,k
!
!  Initialize maximum heating to zero
!
      Hmax = 0.0

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
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - p%rho1gpp
      endif
!
!  velocity damping in the coronal heating zone
!
      if (lhydro.and.tau_damp_cor>0) then
        prof = cubic_step(z(n),(ztop+zcor)/2,(ztop-zcor)/2)
        damp = prof*f(l1:l2,m,n,iux:iuz)/tau_damp_cor
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - damp
      endif
!
!  Advection term
!
      if (ladvection_temperature) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%uglnTT
      endif

!
!  Calculate viscous contribution to temperature
!
      if (lviscosity .and. lviscous_heat) call calc_viscous_heat(f,df,p,Hmax)

!
!  Various heating/cooling mechanisms
!
      if (lcalc_heat_cool) call calc_heat_cool(f,df,p)

!
!  Thermal conduction
!
      if (lheatc_chiconst) call calc_heatcond_constchi(df,p)
      if (lheatc_hyper3) call calc_heatcond_hyper3(df,p)

! Natalia: thermal conduction for the chemistry case: lheatc_chemistry=true

    !  if (lheatc_chemistry) call calc_heatcond_chemistry(f,df,p)


!
!  Interstellar radiative cooling and UV heating
!
      if (linterstellar) &
          call calc_heat_cool_interstellar(f,df,p,Hmax)
!
!  Need to add left-hand-side of the continuity equation (see manual)
!


      if (ldensity) then
       if (.not. lchemistry) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%gamma1*p%divu/p%delta
       else 
    !    df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%gamma1*p%divu
    !      sum_DYDt=0.
    !       do k=1,nchemspec
    !         sum_DYDt=sum_DYDt+p%cvspec(:,k)*(p%DYDt_reac(:,k)+p%DYDt_diff(:,k))
    !       enddo
    !    df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - f(l1:l2,m,n,ilnTT)*p%cv1(:)*sum_DYDt(:)
       endif
      endif


!
!  Calculate temperature related diagnostics
!
      if (ldiagnos) then
        if (idiag_TTmax/=0) call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        if (idiag_TTm/=0) call sum_mn_name(p%TT,idiag_TTm)
        if (idiag_yHmax/=0) call max_mn_name(p%yH,idiag_yHmax)
        if (idiag_yHmin/=0) call max_mn_name(-p%yH,idiag_yHmin,lneg=.true.)
        if (idiag_yHm/=0) call sum_mn_name(p%yH,idiag_yHm)
        if (idiag_eth/=0) call sum_mn_name(p%ee/p%rho1,idiag_eth)
        if (idiag_ssm/=0) call sum_mn_name(p%ss,idiag_ssm)
        if (idiag_cv/=0) call sum_mn_name(p%cv,idiag_cv)
        if (idiag_cp/=0) call sum_mn_name(p%cp,idiag_cp)
        if (idiag_dtc/=0) then
          call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        endif
        if (idiag_eem/=0) call sum_mn_name(p%ee,idiag_eem)
        if (idiag_ppm/=0) call sum_mn_name(p%pp,idiag_ppm)
        if (idiag_csm/=0) call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
        if (idiag_mum/=0) call sum_mn_name(1/p%mu1,idiag_mum)
      endif

    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
!  calculate chi*grad(rho*T*glnTT)/(rho*TT)
!           =chi*(g2.glnTT+g2lnTT),
!  where g2=glnrho+glnTT
!
      use Sub, only: max_mn_name,dot,del2,multsv

      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx) :: g2,gamma
      real, dimension (nx,3) :: gradlncp

!
!  g2
!
      if (lheatc_chiconst_accurate) then
        call multsv(p%cp1,p%gradcp,gradlncp)
        call dot(p%glnTT+p%glnrho+gradlncp,p%glnTT,g2)
        gamma = p%gamma
      else
        call dot(p%glnTT+p%glnrho,p%glnTT,g2)
        gamma = 5.0/3.0
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

    end subroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_hyper3(df,p)
!
!
!
      use Sub, only: max_mn_name

      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

!
!  Add heat conduction to RHS of temperature equation
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + chi_hyper3*p%del6lnTT

!
!  check maximum diffusion from thermal diffusion
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+chi_hyper3*dxyz_6
        if (ldiagnos.and.idiag_dtchi/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi,l_dt=.true.)
        endif
      endif

    end subroutine calc_heatcond_hyper3
!***********************************************************************
    subroutine calc_heat_cool(f,df,p)

      use Sub, only: cubic_step

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx) :: heat
      real :: prof
!
!  Initialize
!
      heat=0
!
!  Add spatially uniform heating (usually as a test)
!
      if (heat_uniform/=0.0) heat = heat+heat_uniform
!
!  add "coronal" heating (to simulate a hot corona)
!  assume a linearly increasing reference profile
!  This 1/rho1 business is clumpsy, but so would be obvious alternatives...
!
      if (tau_heat_cor>0) then
        prof = cubic_step(z(n),(ztop+zcor)/2,(ztop-zcor)/2)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)  &
                            - prof*(1-TT_cor/p%TT)/tau_heat_cor
      endif
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
        idiag_eem=0; idiag_ppm=0; idiag_csm=0
        idiag_mum=0
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
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'csm',idiag_csm)
        call parse_name(iname,cname(iname),cform(iname),'mum',idiag_mum)
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
        write(3,*) 'i_yHmax=',idiag_yHmax
        write(3,*) 'i_yHmin=',idiag_yHmin
        write(3,*) 'i_yHm=',idiag_yHm
        write(3,*) 'i_eth=',idiag_eth
        write(3,*) 'i_ssm=',idiag_ssm
        write(3,*) 'i_cv=',idiag_cv
        write(3,*) 'i_cp=',idiag_cp
        write(3,*) 'i_dtchi=',idiag_dtchi
        write(3,*) 'i_dtc=',idiag_dtc
        write(3,*) 'i_eem=',idiag_eem
        write(3,*) 'i_ppm=',idiag_ppm
        write(3,*) 'i_csm=',idiag_csm
        write(3,*) 'i_mum=',idiag_mum
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************

    subroutine calc_heatcond_ADI(finit,f)

      use Cparam

      implicit none

      real, dimension(mx,my,mz,mfarray) :: finit,f

    end subroutine calc_heatcond_ADI
!**************************************************************
endmodule Entropy
