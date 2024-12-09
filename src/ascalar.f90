! This module is used to solve the equation of supersaturation
! for either the Smoluchowski approach or the swarm model.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lascalar = .true.
!
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 1
!!!! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED acc
! PENCILS PROVIDED gacc(3); ugacc
! PENCILS PROVIDED del2acc
! PENCILS PROVIDED ssat
! PENCILS PROVIDED ttc
! PENCILS PROVIDED gttc(3); ugttc
! PENCILS PROVIDED del2ttc
!
!***************************************************************
module Ascalar
!
  use Cdata
  use Cparam
  use Messages
!
  implicit none
!
  include 'ascalar.h'
!
! Define local variables in this special module, not in "Cdata.f90"
!  integer :: issat=0
!  Init parameters.
!
  real :: acc_const=0., amplacc=0., widthacc=0., ttc_const=0., amplttc=0., widthttc=0.
  real :: z0_acc=0.
  logical :: noascalar=.false., reinitialize_acc=.false.
  character (len=labellen) :: initacc='nothing'
  character (len=labellen) :: initttc='nothing'
  real :: T_env=293.0, qv_env=1.63e-2, Rv_over_Rd_minus_one=0.608, gravity_acceleration=9.81
  real, dimension(1) :: ttc_mean=293.0, acc_mean=1.e-2
  logical :: lbuoyancy=.false., ltauascalar=.false., lttc=.false., lttc_mean=.false.
!
  namelist /ascalar_init_pars/ &
           initacc, acc_const, amplacc, widthacc, z0_acc, &
           initttc, ttc_const, amplttc, widthttc, & 
           T_env, qv_env, lbuoyancy, lttc, lttc_mean
!
!  Run parameters.
!
  real :: thermal_diff=0.0
  real :: ascalar_diff=0.0
  real :: vapor_mixing_ratio_qvs=0.0
  real :: updraft=0.0
  real :: A1=0.0
  real :: latent_heat=2.5e6, cp_constant=1005.0
  real :: Rv=461.5, rhoa=1.06, const1_qvs=2.53e11, const2_qvs=5420.0
  real, target :: ssat0=0.0
  real :: TT_mean=293.25, constTT=293.25
  real, dimension(3) :: gradacc0=(/0.0,0.0,0.0/)
  real, dimension(3) :: gradTT0=(/0.0,0.0,0.0/)
  real, dimension(nx) :: es_T=0.0, qvs_T=0.0
  real, dimension(nx) :: buoyancy=0.0
  logical :: lascalar_sink=.false., Rascalar_sink=.false.,lupdraft=.false., l_T_source=.true.
  logical :: lupw_acc=.false., lcondensation_rate=.false., lconstTT=.false., lTT_mean=.false., lupw_ttc=.false.

  namelist /ascalar_run_pars/ &
      lupw_acc, lascalar_sink, Rascalar_sink, l_T_source, &
      ascalar_diff, gradacc0, lcondensation_rate, vapor_mixing_ratio_qvs, &
      lupdraft, updraft, A1, latent_heat, cp_constant, &
      const1_qvs, const2_qvs, Rv, rhoa, gravity_acceleration, Rv_over_Rd_minus_one, &
      lconstTT, constTT, TT_mean, lTT_mean, thermal_diff, gradTT0
!
!  Diagnostics variables
!
  integer :: idiag_accrms=0, idiag_accmax=0, idiag_accmin=0, idiag_accm=0
  integer :: idiag_ttcrms=0, idiag_ttcmax=0, idiag_ttcmin=0, idiag_ttcm=0
  integer :: idiag_uxaccm=0, idiag_uyaccm=0, idiag_uzaccm=0
  integer :: idiag_tauascalarrms=0, idiag_tauascalarmax=0, idiag_tauascalarmin=0
  integer :: idiag_condensationRaterms=0, idiag_condensationRatemax=0,idiag_condensationRatemin=0,idiag_condensationRatem=0
  integer :: idiag_waterMixingRatiorms=0, idiag_waterMixingRatiomax=0,idiag_waterMixingRatiomin=0,idiag_waterMixingRatiom=0
  integer :: idiag_ssatrms=0, idiag_ssatmax=0, idiag_ssatmin=0, idiag_ssatm=0
  integer :: idiag_buoyancyrms=0, idiag_buoyancym=0, idiag_buoyancymax=0, idiag_buoyancymin=0
  integer :: idiag_esrms=0, idiag_esm=0, idiag_esmax=0, idiag_esmin=0
  integer :: idiag_qvsrms=0, idiag_qvsm=0, idiag_qvsmax=0, idiag_qvsmin=0
  integer :: idiag_ttc_mean=0, idiag_acc_mean=0
  integer :: idiag_accmz=0 ! XYAVG_DOC: $\left<c\right>_{xy}$
!
  contains
!***********************************************************************
    subroutine register_ascalar()
!
!  Initialise the acc variable and increase nvar accordingly
!
!   3-jun-16/xiangyu: adapted from pscalar_nolog
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
      call farray_register_pde('acc', iacc)
!
      call farray_register_auxiliary('ssat', issat)
!
!      call farray_register_auxiliary('ttc', ittc)
      call farray_register_pde('ttc', ittc)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
!      if (lcondensation_rate) then
        call put_shared_variable('ssat0', ssat0, caller='register_ascalar')
!      endif
!
    endsubroutine register_ascalar
!***********************************************************************
    subroutine initialize_ascalar(f)
!
!  Perform any necessary post-parameter read initialization
!  Since the passive scalar is often used for diagnostic purposes
!  one may want to reinitialize it to its initial distribution.
!      
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (lroot) print*, 'Supersaturation routine'
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if (reinitialize_acc) then
        f(:,:,:,iacc)=0.
        f(:,:,:,ittc)=0.
        call init_acc(f)
      endif
!
    endsubroutine initialize_ascalar
!***********************************************************************
    subroutine init_acc(f)
!
!  initialise passive scalar field; called from start.f90
!
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: l
      real, dimension (mz) :: tmp
!
      select case (initacc)
        case ('nothing')
        case ('zero'); f(:,:,:,iacc)=0.0
        case ('constant'); f(:,:,:,iacc)=acc_const
        case ('tanhx')
          do m=m1,m2;do n=n1,n2
             f(:,m,n,iacc)=acc_const+amplacc*tanh(x/widthacc)
          enddo;enddo
        case ('tanhz')
          do l=l1,l2; do m=m1,m2
             f(l,m,:,iacc)=acc_const+amplacc*tanh(z/widthacc)
          enddo;enddo
        case ('tanhz/(1-tanhz)')
          !Such that the mass fraction is a step.
          do l=l1,l2; do m=m1,m2
            tmp = acc_const + amplacc*(1+tanh((z-z0_acc)/widthacc))/2
            if (any(tmp==1)) call fatal_error('init_acc', &
              'specified initial condition leads to infinite values for acc')
            f(l,m,:,iacc) = tmp/(1-tmp)
          enddo;enddo
!
!  Catch unknown values.
!
        case default
          call fatal_error('initacc','initacc value not recognised')
      endselect
!
!  initialise temperature field
      select case (initttc)
        case ('nothing')
        case ('zero'); f(:,:,:,ittc)=0.0
        case ('constant'); f(:,:,:,ittc)=ttc_const
        case ('tanhx')
          do m=m1,m2;do n=n1,n2
             f(:,m,n,ittc)=ttc_const+amplttc*tanh(x/widthttc)
          enddo;enddo
        case ('tanhz')
          do l=l1,l2; do m=m1,m2
             f(l,m,:,ittc)=ttc_const+amplttc*tanh(z/widthttc)
          enddo;enddo
!
!  Catch unknown values.
!
        case default
          call fatal_error('initttc','initttc value not recognised')
      endselect
!
!  modify the density to have lateral pressure equilibrium
!
! XY: These lines cause the short time step, so I commented them out at the moment.
!!print*,'AXEL6', mudry1, muvap1
!       mudry1=1./29.
!       muvap1=1./18.
!       f(:,:,:,ilnrho)=f(:,:,:,ilnrho)-alog((1-f(:,:,:,iacc))*mudry1+f(:,:,:,iacc)*muvap1)
!!XX
!
    endsubroutine init_acc
!***********************************************************************
    subroutine pencil_criteria_ascalar()
!
!  All pencils that the Ascalar module depends on are specified here.
!
      integer :: i
!
!  acc itself always
!
      lpenc_requested(i_acc)=.true.
      lpenc_requested(i_ssat)=.true.
      lpenc_requested(i_ttc)=.true.
!
!  background gradient 
!
! ascalar
!
      do i=1,3
        if (gradacc0(i)/=0.) lpenc_requested(i_uu)=.true.
      enddo
!
!      if (lascalar_sink) then 
      if (lparticles .and. lascalar_sink) then 
        lpenc_requested(i_tauascalar)=.true.
      endif
      if (ascalar_diff/=0.) lpenc_requested(i_del2acc)=.true.
! 
      lpenc_diagnos(i_acc)=.true.
!
! temperature
!
      if (lttc) then
        do i=1,3
          if (gradTT0(i)/=0.) lpenc_requested(i_uu)=.true.
        enddo
      endif
!
      if (thermal_diff/=0.) lpenc_requested(i_del2ttc)=.true.
!
      lpenc_diagnos(i_ttc)=.true.
      if (idiag_accmz/=0) lpenc_diagnos(i_acc)=.true.
!      
! temperature calculated from "temperature_idealgas.f90"
      if (ltemperature) then 
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_cp)=.true.
      endif
!     
    endsubroutine pencil_criteria_ascalar
!***********************************************************************
    subroutine pencil_interdep_ascalar(lpencil_in)
!  Interdependency among pencils provided by the Pscalar module
!  is specified here.
      logical, dimension(npencils) :: lpencil_in
!
      lpencil_in(i_acc)=.true.
      lpencil_in(i_ugacc)=.true.
      if (lpencil_in(i_ugacc)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gacc)=.true.
      endif
!      
      lpencil_in(i_ttc)=.true.
      lpencil_in(i_ugttc)=.true.
      if (lpencil_in(i_ugttc)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gttc)=.true.
      endif
!      
      lpencil_in(i_ssat)=.true.
!
      if (ltauascalar) lpencil_in(i_tauascalar)=.true.
      lpencil_in(i_condensationRate)=.true.
      lpencil_in(i_waterMixingRatio)=.true.
!      
    endsubroutine pencil_interdep_ascalar
!**********************************************************************
    subroutine ascalar_after_boundary(f)

      real, dimension (mx,my,mz,mfarray), intent(IN) :: f

      if (lcondensation_rate.and.lttc.and.lbuoyancy.and.lttc_mean) then
        call calc_ttcmean(f)
        call calc_accmean(f)
      endif

    endsubroutine ascalar_after_boundary
!**********************************************************************
    subroutine calc_pencils_ascalar(f,p)
!
!  Calculate ascalar Pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      integer :: j
!
      intent(in) :: f
      intent(inout) :: p
!
! acc
      if (lpencil(i_acc)) p%acc=f(l1:l2,m,n,iacc)
      if (issat>0) p%ssat=f(l1:l2,m,n,issat)
      if (lpencil(i_ttc)) p%ttc=f(l1:l2,m,n,ittc)
!
!  Compute gacc. Add imposed spatially constant gradient of acc.
!  (This only makes sense for periodic boundary conditions.)
!
      if (lpencil(i_gacc)) then
        call grad(f,iacc,p%gacc)
        do j=1,3
          if (gradacc0(j)/=0.) p%gacc(:,j)=p%gacc(:,j)+gradacc0(j)
        enddo
      endif
! ugacc
      if (lpencil(i_ugacc)) call u_dot_grad(f,iacc,p%gacc,p%uu,p%ugacc,UPWIND=lupw_acc)
! del2acc
      if (lpencil(i_del2acc)) call del2(f,iacc,p%del2acc)
!
!  Compute gttc. Add imposed spatially constant gradient of ttc.
!  (This only makes sense for periodic boundary conditions.)
!
      if (lpencil(i_gttc)) then
        call grad(f,ittc,p%gttc)
        do j=1,3
          if (gradTT0(j)/=0.) p%gttc(:,j)=p%gttc(:,j)+gradTT0(j)
        enddo
      endif
! ugttc
      if (lpencil(i_ugttc)) call u_dot_grad(f,ittc,p%gttc,p%uu,p%ugttc,UPWIND=lupw_ttc)
! del2ttc
      if (lpencil(i_del2ttc)) call del2(f,ittc,p%del2ttc)
!
    endsubroutine calc_pencils_ascalar
!***********************************************************************
    subroutine dacc_dt(f,df,p)
!
!  Active scalar evolution for supersation
!  Calculate dacc/dt=-uu.gacc + supersat_diff*[del2acc + glnrho.gacc].
!
!  27-may-16/xiangyu: adapted from pscalar_nolog
!   4-sep-16/axel: added more diagnostics
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(out) :: df
!
      real, dimension (nx) :: bump, radius_sum, condensation_rate_Cd
!
      character(len=2) :: id
!
!  Identify module and boundary conditions.
!
      if (noascalar) then
        if (headtt.or.ldebug) print*,'not SOLVED: dacc_dt'
      else
        if (headtt.or.ldebug) print*,'SOLVE dacc_dt'
      endif
      if (headtt) then
        write(id,'(i0)')
        call identify_bcs('acc'//trim(id),iacc)
      endif
!
! Initialize auxiliary variables defined in this module and the corresponding pencils
!
      f(:,m,n,issat) = 0.0
!
!  Passive scalar equation.
! ascalar
      df(l1:l2,m,n,iacc)=df(l1:l2,m,n,iacc)-p%ugacc
      if (ascalar_diff/=0.) then
        df(l1:l2,m,n,iacc)=df(l1:l2,m,n,iacc)+ascalar_diff*p%del2acc
        if (lfirst.and.ldt) maxdiffus=max(maxdiffus,ascalar_diff)
      endif
!
! ttc
!
      if (lttc) then
        df(l1:l2,m,n,ittc)=df(l1:l2,m,n,ittc)-p%ugttc
        if (thermal_diff/=0.) then
          df(l1:l2,m,n,ittc)=df(l1:l2,m,n,ittc)+thermal_diff*p%del2ttc
          if (lfirst.and.ldt) maxdiffus=max(maxdiffus,thermal_diff)
        endif
      endif
!
      ! 1-June-16/XY coded 
      if (lascalar_sink) then
        if (Rascalar_sink) then
          bump=A1*p%uu(:,3)
        elseif (lupdraft) then
          bump=A1*updraft
        else
          bump=A1*p%uu(:,3)-f(l1:l2,m,n,iacc)*p%tauascalar
        endif
        df(l1:l2,m,n,iacc)=df(l1:l2,m,n,iacc)+bump
      endif
!
      if (ldustdensity.and.lcondensation_rate) then
!
!  compute sum of particle radii
!
        f(l1:l2,m,n,issat)=f(l1:l2,m,n,issat)+f(l1:l2,m,n,iacc)/vapor_mixing_ratio_qvs-1.
        radius_sum=0.
        condensation_rate_Cd=4.*pi*f(l1:l2,m,n,issat)*p%condensationRate
        df(l1:l2,m,n,iacc)=df(l1:l2,m,n,iacc)-condensation_rate_Cd
       ! if (ltemperature) df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+condensation_rate_Cd/(p%cp*p%TT)
        if (ltemperature) then
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+condensation_rate_Cd/(p%cp*p%TT)
          if (l_T_source) df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)-gravity_acceleration/p%cp*p%uu(:,3)
        endif
      endif
!
! 17-10-18: Xiang-Yu coded. Solve supersaturation by solving equations for the temperature, mixing ratio
!      
      if (lcondensation_rate) then
        df(l1:l2,m,n,iacc)=df(l1:l2,m,n,iacc)-p%condensationRate
        if (lconstTT) then
          es_T=const1_qvs*exp(-const2_qvs/constTT)
          qvs_T=es_T/(Rv*rhoa*constTT)
          ssat0=acc_const/qvs_T(1)-1
        elseif (ltemperature) then
          df(l1:l2,m,n,iTT)=df(l1:l2,m,n,iTT)+p%condensationRate*latent_heat/cp_constant
          if (lbuoyancy) then
            if (lTT_mean) then
              buoyancy = gravity_acceleration*((p%TT+TT_mean-T_env)/(p%TT+TT_mean) &
                        +Rv_over_Rd_minus_one*(p%acc-qv_env)/p%acc-p%waterMixingRatio)
            else
              buoyancy = gravity_acceleration*((p%TT-T_env)/p%TT &
                        +Rv_over_Rd_minus_one*(p%acc-qv_env)/p%acc-p%waterMixingRatio)
            endif
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+buoyancy
          endif
!
          if (lTT_mean) then 
            es_T=const1_qvs*exp(-const2_qvs/(f(l1:l2,m,n,iTT)+TT_mean))
            qvs_T=es_T/(Rv*rhoa*(f(l1:l2,m,n,iTT)+TT_mean))
          else
            es_T=const1_qvs*exp(-const2_qvs/f(l1:l2,m,n,iTT))
            qvs_T=es_T/(Rv*rhoa*f(l1:l2,m,n,iTT))
          endif       
        elseif (lttc) then
          df(l1:l2,m,n,ittc)=df(l1:l2,m,n,ittc)+p%condensationRate*latent_heat/cp_constant
          es_T=const1_qvs*exp(-const2_qvs/f(l1:l2,m,n,ittc))
          qvs_T=es_T/(Rv*rhoa*f(l1:l2,m,n,ittc))
          ssat0=acc_const/((const1_qvs*exp(-const2_qvs/ttc_const))/(Rv*rhoa*ttc_const))-1
          if (lbuoyancy) then
            if (lttc_mean) then
              buoyancy = gravity_acceleration*((p%ttc-ttc_mean(1))/p%ttc+ &
                         Rv_over_Rd_minus_one*(p%acc-acc_mean(1))/p%acc-p%waterMixingRatio)
            else
              buoyancy = gravity_acceleration*((p%ttc-T_env)/p%ttc+ &
                         Rv_over_Rd_minus_one*(p%acc-qv_env)/p%acc-p%waterMixingRatio)
            endif
            df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+buoyancy
          endif
        endif
        f(l1:l2,m,n,issat)=f(l1:l2,m,n,issat)+f(l1:l2,m,n,iacc)/qvs_T-1.
      endif

      call calc_diagnostics_ascalar(p)

    endsubroutine dacc_dt
!***********************************************************************
    subroutine calc_diagnostics_ascalar(p)
!
!  Diagnostics
!
      use Diagnostics
!
      type (pencil_case) :: p
!
      if (ldiagnos) then
        if (idiag_accrms/=0) call sum_mn_name(p%acc**2,idiag_accrms,lsqrt=.true.)
        call max_mn_name(p%acc,idiag_accmax)
        if (idiag_accmin/=0) call max_mn_name(-p%acc,idiag_accmin,lneg=.true.)
        call sum_mn_name(p%acc,idiag_accm)
        if (lttc) then
          if (idiag_ttcrms/=0) call sum_mn_name(p%ttc**2,idiag_ttcrms,lsqrt=.true.)
          call max_mn_name(p%ttc,idiag_ttcmax)
          if (idiag_ttcmin/=0) call max_mn_name(-p%ttc,idiag_ttcmin,lneg=.true.)
          call sum_mn_name(p%ttc,idiag_ttcm)
        endif
        if (idiag_uxaccm/=0) call sum_mn_name(p%uu(:,1)*p%acc,idiag_uxaccm)
        if (idiag_uyaccm/=0) call sum_mn_name(p%uu(:,2)*p%acc,idiag_uyaccm)
        if (idiag_uzaccm/=0) call sum_mn_name(p%uu(:,3)*p%acc,idiag_uzaccm)
        if (ltauascalar) then
          if (idiag_tauascalarrms/=0) call sum_mn_name(p%tauascalar**2,idiag_tauascalarrms,lsqrt=.true.)
          call max_mn_name(p%tauascalar,idiag_tauascalarmax)
          if (idiag_tauascalarmin/=0) call max_mn_name(-p%tauascalar,idiag_tauascalarmin,lneg=.true.)
        endif
        if (idiag_condensationRaterms/=0) &
          call sum_mn_name(p%condensationRate**2,idiag_condensationRaterms,lsqrt=.true.)
        call max_mn_name(p%condensationRate,idiag_condensationRatemax)
        if (idiag_condensationRatemin/=0) call max_mn_name(-p%condensationRate,idiag_condensationRatemin,lneg=.true.)
        call sum_mn_name(p%condensationRate,idiag_condensationRatem)
        if (idiag_waterMixingRatiorms/=0) &
          call sum_mn_name(p%waterMixingRatio**2,idiag_waterMixingRatiorms,lsqrt=.true.)
        call max_mn_name(p%waterMixingRatio,idiag_waterMixingRatiomax)
        if (idiag_waterMixingRatiomin/=0) call max_mn_name(-p%waterMixingRatio,idiag_waterMixingRatiomin,lneg=.true.)
        call sum_mn_name(p%waterMixingRatio,idiag_waterMixingRatiom)
        if (idiag_ssatrms/=0) call sum_mn_name(p%ssat**2,idiag_ssatrms,lsqrt=.true.)
        call max_mn_name(p%ssat,idiag_ssatmax)
        if (idiag_ssatmin/=0) call max_mn_name(-p%ssat,idiag_ssatmin,lneg=.true.)
        call sum_mn_name(p%ssat,idiag_ssatm)
        call max_mn_name(es_T,idiag_esmax)
        if (idiag_esmin/=0) call max_mn_name(-es_T,idiag_esmin,lneg=.true.)
        call sum_mn_name(es_T,idiag_esm)
        if (idiag_esrms/=0) call sum_mn_name(es_T**2,idiag_esrms,lsqrt=.true.)
        call max_mn_name(qvs_T,idiag_qvsmax)
        if (idiag_qvsmin/=0) call max_mn_name(-qvs_T,idiag_qvsmin,lneg=.true.)
        call sum_mn_name(qvs_T,idiag_qvsm)
        if (idiag_qvsrms/=0) call sum_mn_name(qvs_T**2,idiag_qvsrms,lsqrt=.true.)
        if (lbuoyancy) then
          if (idiag_buoyancyrms/=0) call sum_mn_name(buoyancy**2,idiag_buoyancyrms,lsqrt=.true.)
          call sum_mn_name(buoyancy,idiag_buoyancym)
          call max_mn_name(buoyancy,idiag_buoyancymax)
          if (idiag_buoyancymin/=0) call max_mn_name(-buoyancy,idiag_buoyancymin,lneg=.true.)
        endif
!
        if (lcondensation_rate.and.lttc.and.lbuoyancy.and.lttc_mean) then
          call save_name(acc_mean(1),idiag_acc_mean)
          call save_name(ttc_mean(1),idiag_ttc_mean)
        endif
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1.
!
      if (l1davgfirst .or. (ldiagnos .and. ldiagnos_need_zaverages)) then
        call xysum_mn_name_z(p%acc,idiag_accmz)
      endif
!
    endsubroutine calc_diagnostics_ascalar
!***********************************************************************
    subroutine calc_ttcmean(f)
!
!  Calculation of volume averaged mean temperature.
!
!  06-June-18/Xiang-Yu.Li: coded
!
      use Sub, only: global_mean
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
!  Calculate mean of temperature.
!
      call global_mean(f,ittc,ttc_mean)
!
    endsubroutine calc_ttcmean
!***********************************************************************
    subroutine calc_accmean(f)
!
!  Calculation of volume averaged water vapor mixing ratio.
!
!  06-June-18/Xiang-Yu.Li: coded
!
      use Sub, only: global_mean
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
!  Calculate mean of temperature.
!
      call global_mean(f,iacc,acc_mean)
!
    endsubroutine calc_accmean
!***********************************************************************
    subroutine read_ascalar_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=ascalar_init_pars, IOSTAT=iostat)
!
    endsubroutine read_ascalar_init_pars
!***********************************************************************
    subroutine write_ascalar_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=ascalar_init_pars)
!
    endsubroutine write_ascalar_init_pars
!***********************************************************************
    subroutine read_ascalar_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=ascalar_run_pars, IOSTAT=iostat)
!
    endsubroutine read_ascalar_run_pars
!***********************************************************************
    subroutine write_ascalar_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=ascalar_run_pars)
!
    endsubroutine write_ascalar_run_pars
!
!***********************************************************************
    subroutine rprint_ascalar(lreset,lwrite)
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamez
      logical :: lwr
!
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!      
      if (lreset) then
        idiag_accrms=0; idiag_accmax=0; idiag_accmin=0; idiag_accm=0
        idiag_ttcrms=0; idiag_ttcmax=0; idiag_ttcmin=0; idiag_ttcm=0
        idiag_uxaccm=0; idiag_uyaccm=0; idiag_uzaccm=0
        idiag_tauascalarrms=0; idiag_tauascalarmax=0; idiag_tauascalarmin=0
        idiag_condensationRaterms=0; idiag_condensationRatemax=0; idiag_condensationRatemin=0; idiag_condensationRatem=0
        idiag_waterMixingRatiorms=0; idiag_waterMixingRatiomax=0; idiag_waterMixingRatiomin=0; idiag_waterMixingRatiom=0
        idiag_ssatrms=0; idiag_ssatmax=0; idiag_ssatmin=0; idiag_ssatm=0
        idiag_esrms=0; idiag_esm=0; idiag_esmax=0; idiag_esmin=0
        idiag_qvsrms=0; idiag_qvsm=0; idiag_qvsmax=0; idiag_qvsmin=0
        idiag_buoyancyrms=0; idiag_buoyancym=0; idiag_buoyancymax=0; idiag_buoyancymin=0
        idiag_ttc_mean=0; idiag_acc_mean=0
        idiag_accmz=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'accrms',idiag_accrms)
        call parse_name(iname,cname(iname),cform(iname),'accmax',idiag_accmax)
        call parse_name(iname,cname(iname),cform(iname),'accmin',idiag_accmin)
        call parse_name(iname,cname(iname),cform(iname),'accm',idiag_accm)
        call parse_name(iname,cname(iname),cform(iname),'ttcrms',idiag_ttcrms)
        call parse_name(iname,cname(iname),cform(iname),'ttcmax',idiag_ttcmax)
        call parse_name(iname,cname(iname),cform(iname),'ttcmin',idiag_ttcmin)
        call parse_name(iname,cname(iname),cform(iname),'ttcm',idiag_ttcm)
        call parse_name(iname,cname(iname),cform(iname),'uxaccm',idiag_uxaccm)
        call parse_name(iname,cname(iname),cform(iname),'uyaccm',idiag_uyaccm)
        call parse_name(iname,cname(iname),cform(iname),'uzaccm',idiag_uzaccm)
        call parse_name(iname,cname(iname),cform(iname),'tauascalarrms',idiag_tauascalarrms)
        call parse_name(iname,cname(iname),cform(iname),'tauascalarmax',idiag_tauascalarmax)
        call parse_name(iname,cname(iname),cform(iname),'tauascalarmin',idiag_tauascalarmin)
        call parse_name(iname,cname(iname),cform(iname),'condensationRaterms',idiag_condensationRaterms)
        call parse_name(iname,cname(iname),cform(iname),'condensationRatemax',idiag_condensationRatemax)
        call parse_name(iname,cname(iname),cform(iname),'condensationRatemin',idiag_condensationRatemin)
        call parse_name(iname,cname(iname),cform(iname),'condensationRatem',idiag_condensationRatem)
        call parse_name(iname,cname(iname),cform(iname),'waterMixingRatiorms',idiag_waterMixingRatiorms)
        call parse_name(iname,cname(iname),cform(iname),'waterMixingRatiomax',idiag_waterMixingRatiomax)
        call parse_name(iname,cname(iname),cform(iname),'waterMixingRatiomin',idiag_waterMixingRatiomin)
        call parse_name(iname,cname(iname),cform(iname),'waterMixingRatiom',idiag_waterMixingRatiom)
        call parse_name(iname,cname(iname),cform(iname),'ssatrms',idiag_ssatrms)
        call parse_name(iname,cname(iname),cform(iname),'ssatmax',idiag_ssatmax)
        call parse_name(iname,cname(iname),cform(iname),'ssatmin',idiag_ssatmin)
        call parse_name(iname,cname(iname),cform(iname),'ssatm',idiag_ssatm)
        call parse_name(iname,cname(iname),cform(iname),'esrms',idiag_esrms)
        call parse_name(iname,cname(iname),cform(iname),'esmax',idiag_esmax)
        call parse_name(iname,cname(iname),cform(iname),'esmin',idiag_esmin)
        call parse_name(iname,cname(iname),cform(iname),'esm',idiag_esm)
        call parse_name(iname,cname(iname),cform(iname),'qvsrms',idiag_qvsrms)
        call parse_name(iname,cname(iname),cform(iname),'qvsmax',idiag_qvsmax)
        call parse_name(iname,cname(iname),cform(iname),'qvsmin',idiag_qvsmin)
        call parse_name(iname,cname(iname),cform(iname),'qvsm',idiag_qvsm)
        call parse_name(iname,cname(iname),cform(iname),'buoyancyrms',idiag_buoyancyrms)
        call parse_name(iname,cname(iname),cform(iname),'buoyancym',idiag_buoyancym)
        call parse_name(iname,cname(iname),cform(iname),'buoyancymax',idiag_buoyancymax)
        call parse_name(iname,cname(iname),cform(iname),'buoyancymin',idiag_buoyancymin)
        call parse_name(iname,cname(iname),cform(iname),'ttc_mean',idiag_ttc_mean)
        call parse_name(iname,cname(iname),cform(iname),'acc_mean',idiag_acc_mean)
      enddo
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'accmz',idiag_accmz)
      enddo
!
      if (lwr) then 
        call farray_index_append('iacc',iacc)
        call farray_index_append('issat',issat)
        call farray_index_append('ittc',ittc)
      endif
!
    endsubroutine rprint_ascalar 
!***********************************************************************
endmodule Ascalar
