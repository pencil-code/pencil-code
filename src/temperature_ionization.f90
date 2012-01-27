! $Id$
!
!  This module takes care of entropy (initial condition
!  and time advance)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .false.
! CPARAM logical, parameter :: ltemperature = .true.
! CPARAM logical, parameter :: lthermal_energy = .false.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED Ma2; uglnTT; ugTT; cvspec(nchemspec)
!
!***************************************************************
module Entropy
!
  use Cparam
  use Cdata
  use Messages
  use Interstellar
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'entropy.h'
!
  real :: radius_lnTT=0.1,ampl_lnTT=0.,widthlnTT=2*epsi
  real :: lnTT_left=1.0,lnTT_right=1.0,lnTT_const=0.0,TT_const=1.0
  real :: kx_lnTT=1.0,ky_lnTT=1.0,kz_lnTT=1.0
  real :: chi=0.0,heat_uniform=0.0,chi_hyper3=0.0
  real :: zbot=0.0,ztop=0.0
  real :: tau_heat_cor=-1.0,tau_damp_cor=-1.0,zcor=0.0,TT_cor=0.0
  logical, pointer :: lpressuregradient_gas
  logical :: ladvection_temperature=.true.
  logical :: lviscosity_heat=.false.
  logical :: lupw_lnTT=.false.,lcalc_heat_cool=.false.
  logical :: lheatc_chiconst=.false.,lheatc_chiconst_accurate=.false.
  logical :: lheatc_hyper3=.false.
  integer, parameter :: nheatc_max=3
  character (len=labellen), dimension(ninit) :: initlnTT='nothing'
  character (len=labellen), dimension(nheatc_max) :: iheatcond='nothing'
  character (len=intlen) :: iinit_str
!
  namelist /entropy_init_pars/ &
      initlnTT,radius_lnTT,ampl_lnTT,widthlnTT, &
      lnTT_left,lnTT_right,lnTT_const,TT_const, &
      kx_lnTT,ky_lnTT,kz_lnTT,ltemperature_nolog
!
  namelist /entropy_run_pars/ &
      lupw_lnTT,ladvection_temperature, &
      heat_uniform,chi,tau_heat_cor,tau_damp_cor,zcor,TT_cor, &
      lheatc_chiconst_accurate,lheatc_hyper3,chi_hyper3, &
      iheatcond
!
  integer :: idiag_TTmax=0    ! DIAG_DOC: $\max (T)$
  integer :: idiag_TTmin=0    ! DIAG_DOC: $\min (T)$
  integer :: idiag_TTm=0      ! DIAG_DOC: $\left< T \right>$
  integer :: idiag_yHmax=0    ! DIAG_DOC:
  integer :: idiag_yHmin=0    ! DIAG_DOC:
  integer :: idiag_yHm=0      ! DIAG_DOC:
  integer :: idiag_ethm=0     ! DIAG_DOC: $\left< e_{\text{th}}\right> =
                              ! DIAG_DOC:  \left< c_v \rho T \right> $
                              ! DIAG_DOC: \quad(mean thermal energy)
  integer :: idiag_ssm=0      ! DIAG_DOC:
  integer :: idiag_cv=0
  integer :: idiag_cp=0
  integer :: idiag_dtchi=0
  integer :: idiag_dtc=0      ! DIAG_DOC:
  integer :: idiag_eem=0      ! DIAG_DOC: $\left< e \right> =
                              ! DIAG_DOC:  \left< c_v T \right>$
                              ! DIAG_DOC: \quad(mean internal energy)
  integer :: idiag_ppm=0      ! DIAG_DOC:
  integer :: idiag_csm=0
  integer :: idiag_mum=0      ! DIAG_DOC:
  integer :: idiag_ppmax=0    ! DIAG_DOC:
  integer :: idiag_ppmin=0    ! DIAG_DOC:
!
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
      use SharedVariables, only: get_shared_variable
!
      integer :: ierr
!
      call farray_register_pde('lnTT',ilnTT)
!
      call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas,ierr)
      if (ierr/=0) call fatal_error('register_entropy','lpressuregradient_gas')
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
    subroutine initialize_entropy(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
!
      use Mpicomm, only: stop_it
      use SharedVariables, only: put_shared_variable
      use EquationOfState, only : select_eos_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
      integer :: ierr,i
!
!  Set iTT requal to ilnTT if we are considering non-logarithmic temperature.
!
      if (ltemperature_nolog) iTT=ilnTT
!
      if (ltemperature_nolog) then
        call select_eos_variable('TT',iTT)
      else
        call select_eos_variable('lnTT',ilnTT)
      endif
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
!
!  put lviscosity_heat as shared variable for viscosity module
!
      call put_shared_variable('lviscosity_heat',lviscosity_heat,ierr)
      if (ierr/=0) call stop_it("initialize_entropy: "//&
           "there was a problem when putting lviscosity_heat")
!
!  Set iTT equal to ilnTT if we are considering non-logarithmic temperature.
!
      if (ltemperature_nolog) iTT=ilnTT
!
      do i=1,nheatc_max
        select case (iheatcond(i))
        case ('chi-const')
          lheatc_chiconst=.true.
          if (lroot) call information('initialize_entropy', &
              ' heat conduction: constant chi')
        case ('chi-hyper3')
          lheatc_hyper3=.true.
          if (lroot) call information('initialize_entropy','hyper conductivity')
        case ('nothing')
          if (lroot) print*,'heat conduction: nothing'
        case default
          if (lroot) then
            write(unit=errormsg,fmt=*)  &
                'No such value iheatcond = ', trim(iheatcond(i))
            call fatal_error('initialize_entropy',errormsg)
          endif
        endselect
      enddo
!
      if (lheatc_chiconst .and. chi==0.0) then
        call warning('initialize_entropy','chi is zero!')
      endif
      if (lheatc_hyper3 .and. chi_hyper3==0.0) then
        call warning('initialize_entropy', &
            'Conductivity coefficient chi_hyper3 is zero!')
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
      endsubroutine initialize_entropy
!***********************************************************************
    subroutine read_entropy_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=entropy_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=entropy_init_pars,ERR=99)
      endif
!
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
!
      if (present(iostat)) then
        read(unit,NML=entropy_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=entropy_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_entropy_run_pars
!***********************************************************************
    subroutine write_entropy_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=entropy_run_pars)
    endsubroutine write_entropy_run_pars
!!**********************************************************************
    subroutine init_ss(f)
!
!  initialise entropy; called from start.f90
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!
      use General, only: itoa
      use Sub, only: blob
      use Initcond, only: jump, gaunoise
      use InitialCondition, only: initial_condition_ss
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
!
          iinit_str=itoa(j)
!
!  select different initial conditions
!
          select case (initlnTT(j))
!
          case ('zero', '0'); f(:,:,:,ilnTT) = 0.
          case ('const_lnTT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+lnTT_const
          case ('const_TT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+log(TT_const)
          case ('blob'); call blob(ampl_lnTT,f,ilnTT,radius_lnTT,0.,0.,0.)
          case ('xwave')
            do n=n1,n2; do m=m1,m2
              f(l1:l2,m,n,ilnTT)=f(l1:l2,m,n,ilnTT)+ampl_lnTT*sin(kx_lnTT*x(l1:l2))
            enddo; enddo
          case ('ywave')
            do n=n1,n2; do m=m1,m2
              f(l1:l2,m,n,ilnTT)=f(l1:l2,m,n,ilnTT)+ampl_lnTT*sin(ky_lnTT*y(m))
            enddo; enddo
          case ('zwave')
            do n=n1,n2; do m=m1,m2
              f(l1:l2,m,n,ilnTT)=f(l1:l2,m,n,ilnTT)+ampl_lnTT*sin(kz_lnTT*z(n))
            enddo; enddo
          case ('xjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'x')
          case ('yjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'y')
          case ('zjump'); call jump(f,ilnTT,lnTT_left,lnTT_right,widthlnTT,'z')
          case ('gaussian-noise'); call gaunoise(ampl_lnTT,f,ilnTT)
!
          case default
!
!  Catch unknown values
!
            write(unit=errormsg,fmt=*) 'No such value for initss(' &
                //trim(iinit_str)//'): ',trim(initlnTT(j))
            call fatal_error('init_ss',errormsg)
!
          endselect
!
          if (lroot) print*,'init_ss: initss(' &
              //trim(iinit_str)//') = ',trim(initlnTT(j))
!
        endif
!
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_ss(f)
!
!  If unlogarithmic temperature considered, take exp of lnTT resulting from
!  initlnTT
!
      if (ltemperature_nolog) f(:,:,:,iTT)=exp(f(:,:,:,ilnTT))
!
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
!
      if (ldensity) then
       if (.not. lchemistry) then
        lpenc_requested(i_gamma_m1)=.true.
        lpenc_requested(i_delta)=.true.
       endif
        lpenc_requested(i_divu)=.true.
      endif
!
      if (linterstellar) then
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_ee)=.true.
      endif
!
      if (lpressuregradient_gas) then
        lpenc_requested(i_rho1gpp)=.true.
      endif
!
      if (lviscosity) then
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_TT1)=.true.
        if (lviscosity_heat) lpenc_requested(i_visc_heat)=.true.
      endif
!
      if (lcalc_heat_cool) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_TT1)=.true.
        if (tau_heat_cor>0) lpenc_requested(i_TT)=.true.
      endif
!
      if (ladvection_temperature) then
       if (ltemperature_nolog) then
         lpenc_requested(i_ugTT)=.true.
       else
         lpenc_requested(i_uglnTT)=.true.
       endif
      endif
!
      if (lheatc_chiconst) then
        lpenc_requested(i_del2lnTT)=.true.
        if (lheatc_chiconst_accurate) then
          lpenc_requested(i_cp1)=.true.
          lpenc_requested(i_gradcp)=.true.
          lpenc_requested(i_gamma)=.true.
        endif
      endif
!
      if (lheatc_hyper3) lpenc_requested(i_del6lnTT)=.true.
!
      if (ltemperature_nolog) lpenc_requested(i_TT)=.true.
!
    !  if (lchemistry) then
    !     lpenc_requested(i_lnTT)=.true.
    !     lpenc_requested(i_TT)=.true.
    !     lpenc_requested(i_TT1)=.true.
    !     lpenc_requested(i_glnTT)=.true.
    !     lpenc_requested(i_del2lnTT)=.true.
    !  endif
!
!  Diagnostics
!
      if (idiag_TTmax/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_TTmin/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_TTm/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_yHmax/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_yHmin/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_yHm/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_ethm/=0) then
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
      if (idiag_ppmax/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_ppmin/=0) lpenc_diagnos(i_pp)=.true.
      if (idiag_mum/=0) lpenc_diagnos(i_mu1)=.true.
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
      if (lpencil_in(i_uglnTT)) then
        lpencil_in(i_glnTT)=.true.
      endif
      if (lpencil_in(i_ugTT)) then
        lpencil_in(i_gTT)=.true.
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
      use Sub, only: u_dot_grad
!
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
!
      if (lpencil(i_ugTT)) then
        call u_dot_grad(f,iTT,p%gTT,p%uu,p%ugTT,UPWIND=lupw_lnTT)
      endif
!
    endsubroutine calc_pencils_entropy
!***********************************************************************
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
      use Diagnostics, only: max_mn_name,sum_mn_name
      use Special, only: special_calc_entropy
      use Sub, only: cubic_step,identify_bcs
      use Viscosity, only: calc_viscous_heat
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (mx,my,mz,mvar), intent (out) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: damp
      real, dimension (nx) :: Hmax
      real :: prof
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
      if (lpressuregradient_gas) then
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
       if (ltemperature_nolog) then
        df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) - p%ugTT
       else
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%uglnTT
       endif
      endif
!
!  Calculate viscous contribution to temperature
!
      if (lviscosity .and. lviscosity_heat) call calc_viscous_heat(f,df,p,Hmax)
!
!  Various heating/cooling mechanisms
!
      if (lcalc_heat_cool) call calc_heat_cool(df,p)
!
!  Thermal conduction
!
      if (lheatc_chiconst) call calc_heatcond_constchi(df,p)
      if (lheatc_hyper3) call calc_heatcond_hyper3(df,p)
!
! Natalia: thermal conduction for the chemistry case: lheatc_chemistry=true
!
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
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%gamma_m1*p%divu/p%delta
       else
   !    df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%gamma_m1*p%divu
   !      sum_DYDt=0.
   !       do k=1,nchemspec
   !         sum_DYDt=sum_DYDt+p%cvspec(:,k)*(p%DYDt_reac(:,k)+p%DYDt_diff(:,k))
   !       enddo
   !    df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)- &
   !      f(l1:l2,m,n,ilnTT)*p%cv1(:)*sum_DYDt(:)
       endif
      endif
!
      if (lspecial) call special_calc_entropy(f,df,p)
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
        if (idiag_ethm/=0) call sum_mn_name(p%ee/p%rho1,idiag_ethm)
        if (idiag_ssm/=0) call sum_mn_name(p%ss,idiag_ssm)
        if (idiag_cv/=0) call sum_mn_name(p%cv,idiag_cv)
        if (idiag_cp/=0) call sum_mn_name(p%cp,idiag_cp)
        if (idiag_dtc/=0) then
          call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        endif
        if (idiag_eem/=0) call sum_mn_name(p%ee,idiag_eem)
        if (idiag_ppm/=0) call sum_mn_name(p%pp,idiag_ppm)
        if (idiag_ppmax/=0) call max_mn_name(p%pp,idiag_ppmax)
        if (idiag_ppmin/=0) call max_mn_name(-p%pp,idiag_ppmin,lneg=.true.)
        if (idiag_csm/=0) call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
        if (idiag_mum/=0) call sum_mn_name(1/p%mu1,idiag_mum)
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_lentropy_pars(f)
!
!  dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lentropy_pars
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
!  calculate chi*grad(rho*T*glnTT)/(rho*TT)
!           =chi*(g2.glnTT+g2lnTT),
!  where g2=glnrho+glnTT
!
      use Diagnostics, only: max_mn_name
      use Sub!, only: dot,multsv
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
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
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_hyper3(df,p)
!
      use Diagnostics, only: max_mn_name
!
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
!
    endsubroutine calc_heatcond_hyper3
!***********************************************************************
    subroutine calc_heat_cool(df,p)
!
      use Sub, only: cubic_step
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
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
!
    endsubroutine calc_heat_cool
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
        idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0
        idiag_yHmax=0; idiag_yHmin=0; idiag_yHm=0
        idiag_ethm=0; idiag_ssm=0; idiag_cv=0; idiag_cp=0
        idiag_dtchi=0; idiag_dtc=0
        idiag_eem=0; idiag_ppm=0; idiag_csm=0; idiag_ppmax=0; idiag_ppmin=0
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
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'cv',idiag_cv)
        call parse_name(iname,cname(iname),cform(iname),'cp',idiag_cp)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'ppmax',idiag_ppmax)
        call parse_name(iname,cname(iname),cform(iname),'ppmin',idiag_ppmin)
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
        case ('TT')
          slices%yz =exp(f(ix_loc,m1:m2,n1:n2,ilnTT))
          slices%xz =exp(f(l1:l2,iy_loc,n1:n2,ilnTT))
          slices%xy =exp(f(l1:l2,m1:m2,iz_loc,ilnTT))
          slices%xy2=exp(f(l1:l2,m1:m2,iz2_loc,ilnTT))
          if (lwrite_slice_xy3) slices%xy3=exp(f(l1:l2,m1:m2,iz3_loc,ilnTT))
          if (lwrite_slice_xy4) slices%xy4=exp(f(l1:l2,m1:m2,iz4_loc,ilnTT))
          slices%ready=.true.
!  lnTT
        case ('lnTT')
          slices%yz =f(ix_loc,m1:m2,n1:n2,ilnTT)
          slices%xz =f(l1:l2,iy_loc,n1:n2,ilnTT)
          slices%xy =f(l1:l2,m1:m2,iz_loc,ilnTT)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,ilnTT)
          if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,ilnTT)
          if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,ilnTT)
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_entropy
!***********************************************************************
    subroutine fill_farray_pressure(f)
!
!  18-feb-10/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine fill_farray_pressure
!***********************************************************************
    subroutine impose_energy_floor(f)
!
!  Dummy subroutine; may not be necessary for lnTT
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_energy_floor
!***********************************************************************
    subroutine dynamical_thermal_diffusion(umax)
!
!  Dummy subroutine
!
      real, intent(in) :: umax
!
      call keep_compiler_quiet(umax)
      call fatal_error('dynamical_thermal_diffusion', 'not implemented yet')
!
    endsubroutine dynamical_thermal_diffusion
!***********************************************************************
endmodule Entropy
