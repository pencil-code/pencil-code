! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  real :: Kpara=0.,Kperp=0.
  real :: cool_RTV=0.,heatamp=0.
  real :: hyper3_chi=0.
  real :: tau_inv_newton=0.,exp_newton=0.
!
  namelist /special_run_pars/ &
      Kpara,Kperp,cool_RTV,hyper3_chi,heatamp,tau_inv_newton, &
      exp_newton
!
! variables for print.in
!
  integer :: idiag_dtchi2=0   ! DIAG_DOC: $\delta t / [c_{\delta t,{\rm v}}\,
                              ! DIAG_DOC:   \delta x^2/\chi_{\rm max}]$
                              ! DIAG_DOC:   \quad(time step relative to time
                              ! DIAG_DOC:   step based on heat conductivity;
                              ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_dtrad=0    ! DIAG_DOC: radiative loss from RTV
  integer :: idiag_dtnewt=0
!
!  variables for video slices:
!
  real, target, dimension (nx,ny) :: spitzer_xy,spitzer_xy2,spitzer_xy3,spitzer_xy4
  real, target, dimension (nx,nz) :: spitzer_xz
  real, target, dimension (ny,nz) :: spitzer_yz
!
  real, target, dimension (nx,ny) :: rtv_xy,rtv_xy2,rtv_xy3,rtv_xy4
  real, target, dimension (nx,nz) :: rtv_xz
  real, target, dimension (ny,nz) :: rtv_yz
!
!  miscellaneous variables 
  real, dimension (mz) :: lnTT_init_prof,lnrho_init_prof
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules and write svn id.
!
!  10-sep-10/bing: coded
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  Called by start.f90 together with lstarting=.true.   and then
!  called by run.f90   together with lstarting=.false.  after reading
!  parameters, but before the time loop.
!
!  13-sep-10/bing: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      logical, intent(in) :: lstarting
      real, dimension (mz) :: ztmp
      character (len=*), parameter :: filename='/strat.dat'
      integer :: lend,unit=12
      real :: dummy=1.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
      inquire(IOLENGTH=lend) dummy
!
      if (.not.lstarting.and.tau_inv_newton/=0) then
        open(unit,file=trim(directory_snap)//filename, &
            form='unformatted',status='unknown',recl=lend*mz)
        read(unit) ztmp
        read(unit) lnrho_init_prof
        read(unit) lnTT_init_prof
        close(unit)
      endif
!
    endsubroutine initialize_special
!***********************************************************************
  subroutine read_special_run_pars(unit,iostat)
!
!  04-sep-10/bing: coded
!
    integer, intent(in) :: unit
    integer, intent(inout), optional :: iostat
!
    if (present(iostat)) then
      read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
    else
      read(unit,NML=special_run_pars,ERR=99)
    endif
!
 99    return
!
  endsubroutine read_special_run_pars
!***********************************************************************
  subroutine write_special_run_pars(unit)
!
!  04-sep-10/bing: coded
!
    integer, intent(in) :: unit
!
    write(unit,NML=special_run_pars)
!
  endsubroutine write_special_run_pars
!***********************************************************************
  subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  04-sep-10/bing: coded
!
    if (Kpara/=0) then
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_bij)=.true.
      lpenc_requested(i_lnTT)=.true.
      lpenc_requested(i_glnTT)=.true.
      lpenc_requested(i_hlnTT)=.true.
      lpenc_requested(i_lnrho)=.true.
      lpenc_requested(i_glnrho)=.true.
    endif
!
    if (cool_RTV/=0) then
      lpenc_requested(i_cp1)=.true.
      lpenc_requested(i_lnTT)=.true.
      lpenc_requested(i_lnrho)=.true.
    endif
!
  endsubroutine pencil_criteria_special
!***********************************************************************
  subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!  04-sep-10/bing: coded
!
    use Diagnostics, only: parse_name
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
      idiag_dtchi2=0.
      idiag_dtrad=0.
      idiag_dtnewt=0
    endif
!
!  iname runs through all possible names that may be listed in print.in
!
    do iname=1,nname
      call parse_name(iname,cname(iname),cform(iname),'dtchi2',idiag_dtchi2)
      call parse_name(iname,cname(iname),cform(iname),'dtrad',idiag_dtrad)
      call parse_name(iname,cname(iname),cform(iname),'dtnewt',idiag_dtnewt)
    enddo
!
!  write column where which variable is stored
!
    if (lwr) then
      write(3,*) 'i_dtchi2=',idiag_dtchi2
      write(3,*) 'i_dtrad=',idiag_dtrad
      write(3,*) 'i_dtnewt=',idiag_dtnewt
    endif
!
  endsubroutine rprint_special
!***********************************************************************
  subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  204-sep-10/bing: coded
!
    real, dimension (mx,my,mz,mfarray) :: f
    type (slice_data) :: slices
!
!  Loop over slices
!
    select case (trim(slices%name))
!
!  Magnetic vector potential (code variable)
!
    case ('spitzer')
      slices%yz =>spitzer_yz
      slices%xz =>spitzer_xz
      slices%xy =>spitzer_xy
      slices%xy2=>spitzer_xy2
      if (lwrite_slice_xy3) slices%xy3=>spitzer_xy3
      if (lwrite_slice_xy4) slices%xy4=>spitzer_xy4
      slices%ready=.true.
!
    case ('rtv')
      slices%yz =>rtv_yz
      slices%xz =>rtv_xz
      slices%xy =>rtv_xy
      slices%xy2=>rtv_xy2
      if (lwrite_slice_xy3) slices%xy3=>rtv_xy3
      if (lwrite_slice_xy4) slices%xy4=>rtv_xy4
      slices%ready=.true.
!
    endselect
!
    call keep_compiler_quiet(f)
!
  endsubroutine get_slices_special
!***********************************************************************
  subroutine special_calc_entropy(f,df,p)
!
! Additional terms to the right hand side of the
! energy equation
!
!  04-sep-10/bing: coded
!
    use Deriv, only: der6
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
    real, dimension (mx,my,mz,mvar), intent(inout) :: df
    type (pencil_case), intent(in) :: p
!
    real, dimension (nx) :: hc,tmp
!
    if (Kpara/=0) call calc_heatcond_spitzer(df,p)
    if (cool_RTV/=0) call calc_heat_cool_RTV(df,p)
    if (heatamp/=0) call calc_artif_heating(df,p)
    if (tau_inv_newton/=0) call calc_heat_cool_newton(df,p)
 !
    if (hyper3_chi/=0) then
      hc(:) = 0.
      call der6(f,ilnTT,tmp,1,IGNOREDX=.true.)
      hc = hc + tmp
      call der6(f,ilnTT,tmp,2,IGNOREDX=.true.)
      hc = hc + tmp
      call der6(f,ilnTT,tmp,3,IGNOREDX=.true.)
      hc = hc + tmp
      hc = hyper3_chi*hc
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + hc
!
!  due to ignoredx hyper3_chi has [1/s]
!
      if (lfirst.and.ldt) diffus_chi3=diffus_chi3  &
          + hyper3_chi
    endif
!
  endsubroutine special_calc_entropy
!***********************************************************************
  subroutine calc_heatcond_spitzer(df,p)
!
!  Calculates heat conduction parallel and perpendicular (isotropic)
!  to magnetic field lines.
!
!  See: Solar MHD; Priest 1982
!
!  04-sep-10/bing: coded
!
    use Diagnostics,     only : max_mn_name
    use EquationOfState, only: gamma, gamma_m1
    use Sub, only: dot2_mn, multsv_mn, tensor_diffusion_coef
!
    real, dimension (mx,my,mz,mvar), intent(inout) :: df
    type (pencil_case), intent(in) :: p
!
    real, dimension (nx,3) :: gvKpara,gvKperp,tmpv,tmpv2
    real, dimension (nx) :: thdiff,chi_spitzer
    real, dimension (nx) :: vKpara,vKperp
    real, dimension (nx) :: gT2_1,gT2,b2,b2_1
!
    integer ::i,j
!
!  Calculate variable diffusion coefficients along pencils.
!
    call dot2_mn(p%bb,b2)
    b2_1=1./max(tini,b2)
!
    vKpara = Kpara * exp(p%lnTT*3.5)
    vKperp = Kperp * b2_1*exp(2.*p%lnrho+0.5*p%lnTT)
!
!  For time step limitiation we have to find the effective heat flux:
!  abs(K) = [K1.delta_ij + (K0-K1).bi.bj].ei.ej
!
    call dot2_mn(p%glnTT,gT2)
    gT2_1=1./max(tini,gT2)
!
    chi_spitzer=0.
    do i=1,3
      do j=1,3
        chi_spitzer=chi_spitzer+ &
            (vKpara-vKperp)*p%bb(:,i)*p%bb(:,j)*p%glnTT(:,i)*p%glnTT(:,j)*b2_1*gT2_1
        if (i==j) chi_spitzer=chi_spitzer+vKperp*p%glnTT(:,i)*p%glnTT(:,j)*gT2_1
      enddo
    enddo
    chi_spitzer = chi_spitzer*exp(-p%lnTT-p%lnrho)
!
!  Calculate gradient of variable diffusion coefficients.
!
    call multsv_mn(3.5*vKpara,p%glnTT,gvKpara)
!
    do i=1,3
      tmpv(:,i)=0.
      do j=1,3
        tmpv(:,i)=tmpv(:,i) + p%bb(:,j)*p%bij(:,j,i)
      enddo
    enddo
    call multsv_mn(2.*b2_1,tmpv,tmpv2)
    tmpv=2.*p%glnrho+0.5*p%glnTT-tmpv2
    call multsv_mn(vKperp,tmpv,gvKperp)
!
!  Calculate diffusion term.
!
    thdiff = 0.
    call tensor_diffusion_coef(p%glnTT,p%hlnTT,p%bij,p%bb,vKperp,vKpara,thdiff,&
        GVKPERP=gvKperp,GVKPARA=gvKpara)
!
    thdiff = thdiff*exp(-p%lnrho-p%lnTT)
!
!  Add to energy equation.
    if (ltemperature) then
      if (ltemperature_nolog) then
        call fatal_error('calc_heatcond_spitzer','not implented for ltemperature_nolog')
      else
        df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + thdiff
      endif
    else if (lentropy) then
      if (pretend_lnTT) then
        df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + thdiff
      else
        call fatal_error('calc_heatcond_spitzer','not implented for lentropy and not pretend_lnTT')
      endif
    endif
!
    if (lfirst.and.ldt) then
      diffus_chi=diffus_chi+gamma*chi_spitzer*dxyz_2
      if (ldiagnos.and.idiag_dtchi2/=0) then
        call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
      endif
    endif
!
  endsubroutine calc_heatcond_spitzer
!***********************************************************************
  subroutine calc_heat_cool_RTV(df,p)
!
!  Computes the radiative loss in the optical thin corona.
!  Zero order estimation for the electron number densities by
!  computing ionization of hydrogen using Saha-equation.
!
!  Electron Temperature should be used for the radiative loss
!  L = n_e * n_H * Q(T_e)
!
!  04-sep-10/bing: coded
!
    use EquationOfState, only: gamma
    use Diagnostics,     only: max_mn_name
!
    real, dimension (mx,my,mz,mvar), intent(inout) :: df
    type (pencil_case), intent(in) :: p
!
    real, dimension (nx) :: lnQ,rtv_cool,lnTT_SI,lnneni
    real, dimension (nx) :: ln_n_K
    real, dimension (nx) :: dE_kB_T
    real :: unit_lnQ
!
    unit_lnQ=3*alog(real(unit_velocity))+&
        5*alog(real(unit_length))+alog(real(unit_density))
    lnTT_SI = p%lnTT + alog(real(unit_temperature))
!
!  calculate ln(ne*ni) :
!  ln(ne*ni) = ln( 1.17*rho^2/(1.34*mp)^2)
!  lnneni = 2*p%lnrho + alog(1.17) - 2*alog(1.34)-2.*alog(real(m_p))
!
    lnneni = 2.*(p%lnrho+61.4412 +alog(real(unit_mass)))
!
!  taking ionization of hydrogen into account using Saha-equation
!
    dE_kB_T = 1.58e5*exp(-p%lnTT)/unit_temperature
!
    ln_n_K = p%lnrho - 1.5*p%lnTT + alog(real(unit_density/unit_temperature**1.5))
!
    ln_n_K = ln_n_K + 12.1313 + dE_kB_T
!
    where (ln_n_K > 15)
      lnneni = lnneni - ln_n_K
    elsewhere (ln_n_K <=15 .and. ln_n_K > -15)
      lnneni = lnneni + 2.*alog((sqrt(1.0+4.0*exp(ln_n_K))-1.0)/2.0)-2.*ln_n_K
    endwhere
!
    lnQ   = get_lnQ(lnTT_SI)
!
    rtv_cool = lnQ-unit_lnQ+lnneni-p%lnTT-p%lnrho
    rtv_cool = gamma*p%cp1*exp(rtv_cool)
!
    rtv_cool = rtv_cool*cool_RTV
!
!     add to temperature equation
!
    if (ltemperature) then
      if (ltemperature_nolog) then
        call fatal_error('calc_heat_cool_RTV','not implemented for ltemperature_nolog')
      else
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)-rtv_cool
      endif
    else if (lentropy) then
      if (pretend_lnTT) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)-rtv_cool
      else
        call fatal_error('calc_heat_cool_RTV','not implemented for lentropy')
      endif
    endif
!
    if (lfirst.and.ldt) then
      dt1_max=max(dt1_max,rtv_cool/cdts)
      if (ldiagnos.and.idiag_dtrad/=0) then
        itype_name(idiag_dtrad)=ilabel_max_dt
        call max_mn_name(rtv_cool/cdts,idiag_dtrad,l_dt=.true.)
      endif
    endif
!
  endsubroutine calc_heat_cool_RTV
!***********************************************************************
    function get_lnQ(lnTT)
!
!  input: lnTT in SI units
!  output: lnP  [p]= W * m^3
!
      real, parameter, dimension (37) :: intlnT = (/ &
          7.74982, 7.9495, 8.18008, 8.39521, 8.71034, 9.24060, 9.67086 &
          , 9.90112, 10.1314, 10.2465, 10.3616, 10.5919, 10.8221, 11.0524 &
          , 11.2827, 11.5129, 11.7432, 11.9734, 12.2037, 12.4340, 12.6642 &
          , 12.8945, 13.1247, 13.3550, 13.5853, 13.8155, 14.0458, 14.2760 &
          , 14.5063, 14.6214, 14.7365, 14.8517, 14.9668, 15.1971, 15.4273 &
          ,  15.6576,  69.0776 /)
      real, parameter, dimension (37) :: intlnQ = (/ &
          -93.9455, -91.1824, -88.5728, -86.1167, -83.8141, -81.6650 &
          , -80.5905, -80.0532, -80.0837, -80.2067, -80.1837, -79.9765 &
          , -79.6694, -79.2857, -79.0938, -79.1322, -79.4776, -79.4776 &
          , -79.3471, -79.2934, -79.5159, -79.6618, -79.4776, -79.3778 &
          , -79.4008, -79.5159, -79.7462, -80.1990, -80.9052, -81.3196 &
          , -81.9874, -82.2023, -82.5093, -82.5477, -82.4172, -82.2637 &
          , -0.66650 /)
!
      real, dimension (nx) :: lnTT,get_lnQ
      real, dimension (nx) :: slope,ordinate
      integer :: i
!
      get_lnQ=-1000.
!
      do i=1,36
        where(lnTT .ge. intlnT(i) .and. lnTT .lt. intlnT(i+1))
          slope=(intlnQ(i+1)-intlnQ(i))/(intlnT(i+1)-intlnT(i))
          ordinate = intlnQ(i) - slope*intlnT(i)
          get_lnQ = slope*lnTT + ordinate
        endwhere
      enddo
!
    endfunction get_lnQ
!***********************************************************************
    subroutine calc_artif_heating(df,p)
!
!  04-sep-10/bing: coded
!
      use EquationOfState, only: gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: heatinput
      real :: z_Mm
      type (pencil_case) :: p
!
! Volumetric heating rate as it can be found
! in the thesis by bingert.
!
! Get height in Mm.
      z_Mm = z(n)*unit_length*1e-6
      !
      ! Compute volumetric heating rate in [W/m^3] as
      ! found in Bingert's thesis.
      heatinput=heatamp*(1e3*exp(-z_Mm/0.2)+1e-4*exp(-z_Mm/10.))
      !
      ! Convert to pencil units.
      heatinput=heatinput/unit_density/unit_velocity**3*unit_length
      !
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)+ &
          p%TT1*p%rho1*gamma*p%cp1*heatinput
!
      if (lfirst.and.ldt) then
        dt1_max=max(dt1_max,heatinput/cdts)
      endif
!
    endsubroutine calc_artif_heating
!***********************************************************************
    subroutine calc_heat_cool_newton(df,p)
!
!  newton cooling dT/dt = -1/tau * (T-T0)
!               dlnT/dt = 1/tau *(1-T0/T-1)
!      where
!        tau = (rho0/rho)^alpha
!
!  13-sep-10/bing: coded
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: lnrho0
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: newton,tau_inv_tmp
!
      if (headtt) print*,'special_calc_entropy: newton cooling',tau_inv_newton
!
!  Get reference temperature
      newton  = exp(lnTT_init_prof(n)-p%lnTT)-1.
!
!  Multiply by density dependend time scale
      tau_inv_tmp = tau_inv_newton * exp(-exp_newton*(lnrho0-p%lnrho))
      newton  = newton * tau_inv_tmp
!
!  Add newton cooling term to entropy
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + newton
!
      if (lfirst.and.ldt) then
        dt1_max=max(dt1_max,tau_inv_tmp/cdts)
        if (ldiagnos.and.idiag_dtnewt/=0) then
          itype_name(idiag_dtnewt)=ilabel_max_dt          
          call max_mn_name(tau_inv_tmp,idiag_dtnewt,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heat_cool_newton
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
