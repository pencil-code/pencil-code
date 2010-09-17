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
  logical :: lgranulation=.false.
  real :: increase_vorticity=15.
!
  namelist /special_run_pars/ &
      Kpara,Kperp,cool_RTV,hyper3_chi,heatamp,tau_inv_newton, &
      exp_newton,lgranulation,increase_vorticity
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
  real, target, dimension (nx,ny) :: rtv_xy,rtv_xy2,rtv_xy3,rtv_xy4
  real, target, dimension (nx,nz) :: rtv_xz
  real, target, dimension (ny,nz) :: rtv_yz
!
!  variables for granulation driver
!
  TYPE point
    real, dimension(2) :: pos
    real, dimension(4) :: data
    type(point),pointer :: next
  end TYPE point
!
  Type(point), pointer :: first
  Type(point), pointer :: previous
  Type(point), pointer :: current
  Type(point), pointer, save :: firstlev
  Type(point), pointer, save :: secondlev
  Type(point), pointer, save :: thirdlev
!
  integer, parameter :: n_gran_level=3
  integer, save, dimension(n_gran_level) :: xrange_arr, yrange_arr
  integer, save, dimension(n_gran_level) :: granr_arr
  real, save, dimension(n_gran_level) :: ampl_arr, lifet_arr
  real, save :: ig, avoid, dxdy2, thresh, pd
  integer, save :: pow
  integer, save, dimension(mseed) :: points_rstate
  real, dimension(nx,ny) :: Ux,Uy
  real, dimension(nx,ny) :: vx,vy,w,avoidarr
!
!  miscellaneous variables
  real, save, dimension (mz) :: lnTT_init_prof,lnrho_init_prof
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
      if (.not.lstarting.and.lgranulation) then
        if (lhydro) then
          call setdrparams()
        else
          call fatal_error &
              ('initialize_special','granulation only works for lhydro=T')
        endif
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
      lpenc_requested(i_rho)=.true.
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
    subroutine special_before_boundary(f)
!
!  Mmodify the f array before the boundaries are communicated.
!
!  13-sep-10/bing: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  Compute photospheric granulation.
      if (lgranulation.and.ipz==0) then
        if (.not. lpencil_check_at_work) call uudriver(f)
      endif
!
    endsubroutine special_before_boundary
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
! !
! !  taking ionization of hydrogen into account using Saha-equation
! !
!     dE_kB_T = 1.58e5*exp(-p%lnTT)/unit_temperature
! !
!     ln_n_K = p%lnrho - 1.5*p%lnTT + alog(real(unit_density/unit_temperature**1.5))
! !
!     ln_n_K = ln_n_K + 12.1313 + dE_kB_T
! !
!     where (ln_n_K > 15)
!       lnneni = lnneni - ln_n_K
!     elsewhere (ln_n_K <=15 .and. ln_n_K > -15)
!       lnneni = lnneni + 2.*alog((sqrt(1.0+4.0*exp(ln_n_K))-1.0)/2.0)-2.*ln_n_K
!     endwhere
!
    lnQ   = get_lnQ(lnTT_SI)
!
    rtv_cool = lnQ-unit_lnQ+lnneni-p%lnTT-p%lnrho
    rtv_cool = gamma*p%cp1*exp(rtv_cool)
!
    rtv_cool = rtv_cool*cool_RTV*(1.-tanh(3e4*(p%rho-1e-4)))/2.
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
    subroutine setdrparams()
!
      real, parameter :: ldif=2.0
      integer :: xrange,yrange
      real :: granr,life_t,ampl
!
! Every granule has 6 values associated with it: data(1-6).
! These contain:
!      x-position, y-position,
!      current amplitude, amplitude at t=t_0,
!      t_0, and life_time.
!
! Gives intergranular area / (granular+intergranular area)
      ig=0.3
!
! Gives exponential power of evolvement. Higher power faster growth/decay.
      pow=2
!
! Fractional difference in granule power
      pd=0.015
!
! Gives average radius of granule + intergranular lane
! (no smaller than 6 grid points across)
!     here radius of granules is 0.8 Mm or bigger (3 times dx)
!
      if (unit_system.eq.'SI') then
        granr=max(0.8*1.e6/unit_length,3*dx,3*dy)
      elseif  (unit_system.eq.'cgs') then
        granr=max(0.8*1.e8/unit_length,3*dx,3*dy)
      else
        call fatal_error("setdrparams","no valid unit system")
        granr=0.
      endif
!
! Fractional distance, where after emergence, no granule can emerge
! whithin this radius.(In order to not 'overproduce' granules).
! This limit is unset when granule reaches t_0.
      avoid=0.8
!
! Lifetime of granule
      life_t=(60.*5./unit_time)
!
      dxdy2=dx**2+dy**2
!
! Typical central velocity of granule(1.5e5 cm/s=1.5km/s)
! Has to be multiplied by the smallest distance, since velocity=ampl/dist
! should now also be dependant on smallest granluar scale resolvable.
!
      if (unit_system.eq.'SI') then
        ampl=sqrt(dxdy2)/granr*0.28e4/unit_velocity
      elseif (unit_system.eq.'cgs') then
        ampl=sqrt(dxdy2)/granr*0.28e6/unit_velocity
      else
        call fatal_error("setdrparams","no valid unit system")
        ampl=0.
      endif
!
! fraction of current amplitude to maximum amplitude to the beginning
! and when the granule disapears
      thresh=0.78
!
      xrange=min(nint(1.5*granr*(1+ig)/dx),nint(nx/2.0)-1)
      yrange=min(nint(1.5*granr*(1+ig)/dy),nint(ny/2.0)-1)
!
      if (lroot) then
        print*,'| solar_corona: settings for granules'
        print*,'-----------------------------------'
        print*,'| radius [Mm]:',granr*unit_length*1e-6
        print*,'| lifetime [min]',life_t*unit_time/60.
        print*,'-----------------------------------'
      endif
!
      granr_arr(1)=granr
      granr_arr(2)=granr*ldif
      granr_arr(3)=granr*ldif*ldif
!
      ampl_arr(1)=ampl
      ampl_arr(2)=ampl/ldif
      ampl_arr(3)=ampl/(ldif*ldif)
!
      lifet_arr(1)=life_t
      lifet_arr(2)=ldif**2*life_t
      lifet_arr(3)=ldif**4*life_t
!
      xrange_arr(1)=xrange
      xrange_arr(2)=min(nint(ldif*xrange),nint(nx/2.-1.))
      xrange_arr(3)=min(nint(ldif*ldif*xrange),nint(nx/2.-1.))
      yrange_arr(1)=yrange
      yrange_arr(2)=min(nint(ldif*yrange),nint(ny/2-1.))
      yrange_arr(3)=min(nint(ldif*ldif*yrange),nint(ny/2-1.))
!
! Don't reset if RELOAD is used
      if (.not.lreloading) then
! Make sure the pointer is defined properly at the first entry
! at the list for each proc and each level.
!
        if (associated(first)) nullify(first)
        if (associated(current)) nullify(current)
        if (associated(previous)) nullify(previous)
        if (associated(firstlev)) nullify(firstlev)
        if (associated(secondlev)) nullify(secondlev)
        if (associated(thirdlev)) nullify(thirdlev)
!
        allocate(firstlev)
        if (associated(firstlev%next)) nullify(firstlev%next)
        allocate(secondlev)
        if (associated(secondlev%next)) nullify(secondlev%next)
        allocate(thirdlev)
        if (associated(thirdlev%next)) nullify(thirdlev%next)
!
        points_rstate(:)=iproc*1000.
!
       endif
!
    endsubroutine setdrparams
!***********************************************************************
    subroutine uudriver(f)
!
! This routine replaces the external computing of granular velocity
! pattern initially written by B. Gudiksen (2004)
!
! It is invoked by setting lgranulation=T in run.in
! additional parameters are
!         Bavoid =0.01 : the magn. field strenght in Tesla at which
!                        no granule is allowed
!         nvod = 5.    : the strength by which the vorticity is
!                        enhanced
!
!  13-sep-10/bing: coded
!
      use General, only: random_seed_wrapper
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer, dimension(mseed) :: global_rstate
!
! Save global random number seed, will be restored after granulation
! is done
      call random_seed_wrapper(GET=global_rstate)
      call random_seed_wrapper(PUT=points_rstate)
!
! Compute granular velocities.
!
      Ux=0.0
      Uy=0.0
!
      call multi_drive3()
!
      call enhance_vorticity()
!      foot point quenching 
!
      f(l1:l2,m1:m2,n1,iux) = Ux
      f(l1:l2,m1:m2,n1,iuy) = Uy
      f(l1:l2,m1:m2,n1,iuz) = 0.
!
      call random_seed_wrapper(GET=points_rstate)
      call random_seed_wrapper(PUT=global_rstate)
!
    endsubroutine uudriver
!***********************************************************************
    subroutine multi_drive3()
!
!  13-sep-10/bing: coded
!
      integer :: level
!
!  loop over the levels and
!  set pointer to the start of the first entry of the n-th level list
!
      do level=1,n_gran_level
        if (associated(first)) nullify(first)
!
        select case (level)
        case (1)
          first => firstlev
          current => firstlev
          if (associated(firstlev%next)) then
            first%next => firstlev%next
            current%next => firstlev%next
          endif
        case (2)
          first => secondlev
          current => secondlev
          if (associated(secondlev%next)) then
            first%next => secondlev%next
            current%next => secondlev%next
          endif
        case (3)
          first => thirdlev
          current => thirdlev
          if (associated(thirdlev%next)) then
            first%next => thirdlev%next
            current%next => thirdlev%next
          endif
        end select
! There is no previous for the first entry
        nullify(previous)
!
! Call driver
        call drive3(level)
!
! Reset and point to the first element which may have changed
        select case (level)
        case (1)
          if (.NOT. associated(firstlev,first)) firstlev => first
          if (.NOT. associated(firstlev%next,first%next)) &
              firstlev%next=>first%next
        case (2)
          if (.NOT. associated(secondlev,first)) secondlev => first
          if (.NOT. associated(secondlev%next,first%next)) &
              secondlev%next=>first%next
        case (3)
          if (.NOT. associated(thirdlev,first)) thirdlev => first
          if (.NOT. associated(thirdlev%next,first%next)) &
            thirdlev%next=>first%next
        end select
!
      enddo  ! level
!
    endsubroutine multi_drive3
!***********************************************************************
    subroutine drive3(level)
!
      integer, intent(in) :: level
!
      call reset_arrays
!
      if (.not.associated(current%next)) then
        call read_points(level)
        if (.not.associated(current%next)) then
          if (itsub==1.and.it==1) call drive_init(level)
          call write_points(level) ! compares to var.dat
          call write_points(level,0) ! compares to VAR0
        endif
      else
        call update_points(level)
        call draw_update(level)
        do
          if (associated(current%next)) then
            call get_next_point
            call draw_update(level)
          else
            exit
          endif
        enddo
        do
          if (minval(avoidarr).eq.1) exit
          call add_point
          call make_new_point(level)
          call draw_update(level)
        enddo
        call reset_pointer
      endif
!
      Ux = Ux + vx
      Uy = Uy + vy
!
    endsubroutine drive3
!***********************************************************************
    subroutine reset_arrays
!
      w(:,:)=0.0
      vx(:,:)=0.0
      vy(:,:)=0.0
      avoidarr(:,:)=0.0
!
    endsubroutine reset_arrays
!***********************************************************************
    subroutine drive_init(level)
!
      use General, only: random_number_wrapper
      use Mpicomm, only: mpisend_real,mpirecv_real, &
          mpisend_logical,mpirecv_logical
!
      integer, intent(in) :: level
      real :: ampl,xrange,yrange,rand
      logical :: lwait_for_points,lneed_points,ltmp
      real, dimension(6) :: tmppoint,tmppoint_recv
      integer :: send_proc,i,j,ipt
!
! Select properties level depending
!
      ampl = ampl_arr(level)
      xrange=xrange_arr(level)
      yrange=yrange_arr(level)
!
! Logicals for the communication of points
      lwait_for_points=.true.
      lneed_points=.true.
!
! Run until every proc has filled his field
!
      do while (lwait_for_points)
! Check if each proc needs a another point
        if (lneed_points) then
!
          call make_new_point(level)
!
! Set randomly some points t0 to the past so they already decay
!
          call random_number_wrapper(rand)
          current%data(3)=t+(rand*2-1)*current%data(4)* &
              (-alog(thresh*ampl/current%data(2)))**(1./pow)
          current%data(1)=current%data(2)*exp(-((t-current%data(3))/current%data(4))**pow)
          call draw_update(level)
!
          tmppoint(1) = current%pos(1) + ipx*nx
          tmppoint(2) = current%pos(2) + ipy*ny
          tmppoint(3:6) = current%data(:)
        else
!
! No new point create set tmp to zero
          tmppoint(:)=0.
        endif
!
! Communicate points to all procs.
!
        do send_proc=0,nprocxy-1    ! send_proc is the reciever.
          if (iproc==send_proc) then
            do i=0,nprocx-1; do j=0,nprocy-1
              ipt=i+nprocx*j
              if (ipt/=send_proc) call mpisend_real(tmppoint,6,ipt,ipt+10*send_proc)
            enddo; enddo
          else
            call mpirecv_real(tmppoint_recv,6,send_proc,iproc+10*send_proc)
!  Check if point received from send_proc is important
            if (sum(tmppoint_recv)/=0.and.loverlapp(tmppoint)) then
              call add_point
              current%pos(1)=tmppoint_recv(1)-ipx*nx
              current%pos(2)=tmppoint_recv(2)-ipy*ny
              current%data(:)=tmppoint_recv(3:6)
              call draw_update(level)
            endif
          endif
        enddo
!
! Check if field is full.
!
        if (minval(avoidarr).eq.1) then
          lwait_for_points=.false.
          lneed_points=.false.
        else
          call add_point
        endif
!
! Communicate to wait for others
! First root collects and then root distributes the result
         if (iproc==0) then
           ltmp=.false.
           do i=0,nprocx-1; do j=0,nprocy-1
             ipt = i+nprocx*j
             if (ipt.ne.0) then
               call mpirecv_logical(ltmp,1,ipt,ipt+222)
               if (ltmp) lwait_for_points=.true.
             endif
           enddo; enddo
         else
           call mpisend_logical(lwait_for_points,1,0,iproc+222)
         endif
!
         if (iproc==0) then
           do i=0,nprocx-1; do j=0,nprocy-1
             ipt = i+nprocx*j
             if (ipt.ne.0) then
               call mpisend_logical(lwait_for_points,1,ipt,ipt+222)
             endif
           enddo; enddo
         else
           call mpirecv_logical(lwait_for_points,1,0,iproc+222)
         endif
      enddo
!
! And reset
!
      call reset_pointer
!
    endsubroutine drive_init
!***********************************************************************
    subroutine read_points(level)
!
      integer, intent(in) :: level
      integer :: unit=12,lend,rn,iostat
      real :: dummy=1.
      real, dimension(6) :: tmppoint
      logical :: ex
!
      character (len=64) :: filename
!
      inquire(IOLENGTH=lend) dummy
!
      write (filename,'("/points_",I1.1,".dat")') level
!
      inquire(file=trim(directory_snap)//trim(filename),exist=ex)
!
      if (ex) then
        open(unit,file=trim(directory_snap)//trim(filename),access="direct",recl=6*lend)
!
        rn=1
        do
          read(unit,iostat=iostat,rec=rn) tmppoint
          if (iostat.eq.0) then
            current%pos(:) =tmppoint(1:2)
            current%data(:)=tmppoint(3:6)
            call draw_update(level)
            call add_point
            rn=rn+1
          else
            nullify(previous%next)
            deallocate(current)
            current => previous
            exit
          endif
        enddo
        close(unit)
!
        if (ip<14) then
          print*,'Proc',iproc,'read',rn-1
          print*,'points for the granulation driver in level',level
        else
          print*,'Read driver points',iproc,level
        endif
!
        call reset_pointer
!
        if (level==n_gran_level) then
          write (filename,'("/seed_",I1.1,".dat")') level
          inquire(file=trim(directory_snap)//trim(filename),exist=ex)
          if (ex) then
            open(unit,file=trim(directory_snap)//trim(filename), &
                status="unknown",access="direct",recl=mseed*lend)
            read(unit,rec=1) points_rstate
            close(unit)
          else
            call fatal_error('read_points','cant find seed list for granules')
          endif
        endif
      else
        print*,'No points lists were found, creating new driver points'
      endif
!
    endsubroutine read_points
!***********************************************************************
    subroutine write_points(level,issnap)
!
!
!  14-sep-10/bing: coded
!
      integer, intent(in) :: level
      integer, optional, intent(in) :: issnap
      integer :: unit=12,lend,rn
      real :: dummy=1.
      real, dimension(6) :: tmppoint
!
      character (len=64) :: filename
!
      inquire(IOLENGTH=lend) dummy
!
      if (present(issnap)) then
        write (filename,'("/points_",I1.1,"_",I3.3,".dat")') level,issnap
      else
        write (filename,'("/points_",I1.1,".dat")') level
      endif
!
      open(unit,file=trim(directory_snap)//trim(filename),access="direct",recl=6*lend)
!
      rn = 1
      do
        tmppoint(1:2)=current%pos
        tmppoint(3:6)=current%data
!
        write(unit,rec=rn) tmppoint
        if (.not.associated(current%next)) exit
        call get_next_point
        rn = rn+1
      enddo
!
      close(unit)
!
! Save seed list for each level. Is needed if levels are spread over 3 procs.
!
      if (present(issnap)) then
        write (filename,'("/seed_",I1.1,"_",I3.3,".dat")') level,issnap
      else
        write (filename,'("/seed_",I1.1,".dat")') level
      endif
      !
      open(unit,file=trim(directory_snap)//trim(filename),access="direct",recl=mseed*lend)
      write(unit,rec=1) points_rstate
      close(unit)
!
      call reset_pointer
!
    endsubroutine write_points
!***********************************************************************
    subroutine make_new_point(level)
!
      use General, only: random_number_wrapper
!
      integer, intent(in) :: level
      integer :: kfind,count,ipos,jpos,i,j
      integer,dimension(nx,ny) :: k
      real :: rand,ampl,life_t
!
      ampl=ampl_arr(level)
      life_t=lifet_arr(level)
      k(:,:)=0; ipos=0; jpos=0
!
      where (avoidarr.eq.0) k=1
!
! Choose and find location of one of them
!
      call random_number_wrapper(rand)
      kfind=int(rand*sum(k))+1
      count=0
      do i=1,nx; do j=1,ny
        if (k(i,j).eq.1) then
          count=count+1
          if (count.eq.kfind) then
            ipos=i
            jpos=j
            exit
          endif
        endif
      enddo; enddo
!
! Create new data for new point
!
      current%pos(1)=ipos
      current%pos(2)=jpos
!
      call random_number_wrapper(rand)
      current%data(2)=ampl*(1+(2*rand-1)*pd)
!
      call random_number_wrapper(rand)
      current%data(4)=life_t*(1+(2*rand-1)/10.)
!
      current%data(3)=t+current%data(4)* &
          (-alog(thresh*ampl/current%data(2)))**(1./pow)
!
      current%data(1)=current%data(2)* &
          exp(-((t-current%data(3))/current%data(4))**pow)
!
    endsubroutine make_new_point
!***********************************************************************
    subroutine draw_update(level)
!
      integer, intent(in) :: level
      real :: xdist,ydist,dist2,dist,wtmp,vv
      integer :: i,ii,j,jj
      real :: dist0,tmp,ampl,granr
      integer :: xrange,yrange
!
      xrange=xrange_arr(level)
      yrange=yrange_arr(level)
      ampl=ampl_arr(level)
      granr=granr_arr(level)
!
! Update weight and velocity for new granule
!
      do jj=int(current%pos(2))-yrange,int(current%pos(2))+yrange
        j = 1+mod(jj-1+ny,ny)
        if (j>=1.and.j<=ny) then
          do ii=int(current%pos(1))-xrange,int(current%pos(1))+xrange
            i = 1+mod(ii-1+nx,nx)
            if (i>=1.and.i<=nx) then
!
              xdist=dx*(ii-current%pos(1))
              ydist=dy*(jj-current%pos(2))
              dist2=max(xdist**2+ydist**2,dxdy2)
              dist=sqrt(dist2)
!
              if (dist.lt.avoid*granr.and.t.lt.current%data(3)) avoidarr(i,j)=1
!
              wtmp=current%data(1)/dist
!
              dist0 = 0.53*granr
              tmp = (dist/dist0)**2
!
              vv=exp(1.)*current%data(1)*tmp*exp(-tmp)
!
              if (wtmp.gt.w(i,j)*(1-ig)) then
                if (wtmp.gt.w(i,j)*(1+ig)) then
                  ! granular area
                  vx(i,j)=vv*xdist/dist
                  vy(i,j)=vv*ydist/dist
                  w(i,j) =wtmp
                else
                  ! intergranular area
                  vx(i,j)=vx(i,j)+vv*xdist/dist
                  vy(i,j)=vy(i,j)+vv*ydist/dist
                  w(i,j) =max(w(i,j),wtmp)
                end if
              endif
              if (w(i,j) .gt. ampl/(granr*(1+ig))) avoidarr(i,j)=1
            endif
          enddo
        endif
      enddo
!
    endsubroutine draw_update
!***********************************************************************
    function loverlapp(tmppoint)
!
      real, dimension(6), intent(in) :: tmppoint
      logical loverlapp
!
      real :: r1,zx,zy,disx,disy
!
      loverlapp =.false.
!
      r1 = sqrt((nx/2.)**2.+(ny/2.)**2.)
      zx = nx *(ipx+ 1/2)
      zy = ny *(ipy+ 1/2)
!
      disx = abs(zx-tmppoint(1))
      disy = abs(zy-tmppoint(2))
!
      if (sqrt(disx**2+disy**2)<r1) then
        loverlapp=.true.
      endif
! check for periodic cases in X-direction
      disx = abs(zx-(tmppoint(1)-nxgrid))
      if (sqrt(disx**2+disy**2)<r1) then
        loverlapp=.true.
      endif
!
      disx = abs(zx-(tmppoint(1)+nxgrid))
      if (sqrt(disx**2+disy**2)<r1) then
        loverlapp=.true.
      endif
! check for periodic cases in Y-direction
      disx = abs(zx-tmppoint(1))
      disy = abs(zy-(tmppoint(2)-nygrid))
      if (sqrt(disx**2+disy**2)<r1) then
        loverlapp=.true.
      endif
!
      disy = abs(zy-(tmppoint(2)+nygrid))
      if (sqrt(disx**2+disy**2)<r1) then
        loverlapp=.true.
      endif
!
    endfunction loverlapp
!***********************************************************************
    subroutine add_point
!
      type(point), pointer :: newpoint
!
      allocate(newpoint)
!
! the next after the last is not defined
      if (associated(newpoint%next)) nullify(newpoint%next)
!
! the actual should have the new point as the next in the list
      current%next => newpoint
!
! the current will be the previous
      previous => current
!
! make the new point as the current one
      current => newpoint
!
    endsubroutine add_point
!***********************************************************************
    subroutine reset_pointer
!
! 14-Sep-10/bing: coded
!
      current => first
      previous => first
!
      if (associated(previous)) nullify(previous)
!
    endsubroutine reset_pointer
!***********************************************************************
    subroutine get_next_point
!
! 14-sep-10/bing: coded
!
      previous => current
      current => current%next
!
    endsubroutine get_next_point
!***********************************************************************
    subroutine enhance_vorticity()
!
      real,dimension(nx,ny) :: wx,wy
      real :: vrms,vtot
!
! Putting sum of velocities back into vx,vy
        vx=Ux
        vy=Uy
!
! Calculating and enhancing rotational part by factor 5
        if (increase_vorticity/=0) then
          call helmholtz(wx,wy)
          !* war vorher 5 ; zum testen auf  50
          ! nvor is now keyword !!!
          vx=(vx+increase_vorticity*wx )
          vy=(vy+increase_vorticity*wy)
        endif
!
! Normalize to given total rms-velocity
        vrms=sqrt(sum(vx**2+vy**2)/(nxgrid*nygrid))+tini
!
        if (unit_system.eq.'SI') then
          vtot=3.*1e3/unit_velocity
        elseif (unit_system.eq.'cgs') then
          vtot=3.*1e5/unit_velocity
        else
          vtot=0.
          call fatal_error('solar_corona','define a valid unit system')
        endif
!
! Reinserting rotationally enhanced velocity field
!
        Ux=vx*vtot/vrms
        Uy=vy*vtot/vrms
!
    endsubroutine enhance_vorticity
!***********************************************************************
    subroutine update_points(level)
!
! Update the amplitude/weight of a point.
!
! 16-sep-10/bing: coded
!
      use Sub, only: notanumber
!
      integer, intent(in) :: level
      real :: ampl
!
      ampl = ampl_arr(level)
!
      do
        if (notanumber(current%data)) &
            call fatal_error('update points','NaN found')
!
! update amplitude
        current%data(1)=current%data(2)* &
            exp(-((t-current%data(3))/current%data(4))**pow)
!
! remove point if amplitude is less than threshold
        if (current%data(1)/ampl.lt.thresh) call remove_point
!
! check if last point is reached
        if (.not. associated(current%next)) exit
!
! if not go to next point
        call get_next_point
      end do
!
      call reset_pointer
!
    endsubroutine update_points
!***********************************************************************
    subroutine remove_point
!
! Remove any pointer from the list.
!
! 12-aug-10/bing: coded
!
      if (associated(current%next)) then
! current is NOT the last one
        if (associated(first,current)) then
! but current is the first point,
! check which level it is
          if (associated(firstlev,first)) firstlev => current%next
          if (associated(secondlev,first)) secondlev => current%next
          if (associated(thirdlev,first)) thirdlev => current%next
          first => current%next
          deallocate(current)
          current => first
          nullify(previous)
        else
! we are in between
          previous%next => current%next
          deallocate(current)
          current => previous%next
        endif
      else
! we are at the end
        deallocate(current)
        current => previous
        nullify(previous)
! BE AWARE THAT PREVIOUS IS NOT POINTING TO THE RIGHT POSITION
      endif
!
    endsubroutine remove_point
!***********************************************************************
    subroutine helmholtz(frx_r,fry_r)
!
! Extracts the rotational part of a 2d vector field
! to increase vorticity of the velocity field.
!
! 16-sep-10/bing: coded
!
      use Fourier, only: fft_xy_parallel
!
      real, dimension(nx,ny) :: kx,ky,k2,filter
      real, dimension(nx,ny) :: fvx_r,fvy_r,fvx_i,fvy_i
      real, dimension(nx,ny) :: frx_r,fry_r,frx_i,fry_i
      real, dimension(nx,ny) :: fdx_r,fdy_r,fdx_i,fdy_i
      real :: k20
!
      fvx_r=vx
      fvx_i=0.
      call fft_xy_parallel(fvx_r,fvx_i)
!
      fvy_r=vy
      fvy_i=0.
      call fft_xy_parallel(fvy_r,fvy_i)
!
! Reference frequency is half the Nyquist frequency.
      k20 = (kx_ny/2.)**2.
!
      kx =spread(kx_fft(ipx*nx+1:(ipx+1)*nx),2,ny)
      ky =spread(ky_fft(ipy*ny+1:(ipy+1)*ny),1,nx)
!
      k2 =kx**2 + ky**2 + tini
!
      frx_r = +ky*(ky*fvx_r - kx*fvy_r)/k2
      frx_i = +ky*(ky*fvx_i - kx*fvy_i)/k2
!
      fry_r = -kx*(ky*fvx_r - kx*fvy_r)/k2
      fry_i = -kx*(ky*fvx_i - kx*fvy_i)/k2
!
      fdx_r = +kx*(kx*fvx_r + ky*fvy_r)/k2
      fdx_i = +kx*(kx*fvx_i + ky*fvy_i)/k2
!
      fdy_r = +ky*(kx*fvx_r + ky*fvy_r)/k2
      fdy_i = +ky*(kx*fvx_i + ky*fvy_i)/k2
!
! Filter out large wave numbers.
      filter = exp(-(k2/k20)**2)
!
      frx_r = frx_r*filter
      frx_i = frx_i*filter
!
      fry_r = fry_r*filter
      fry_i = fry_i*filter
!
      fdx_r = fdx_r*filter
      fdx_i = fdx_i*filter
!
      fdy_r = fdy_r*filter
      fdy_i = fdy_i*filter
!
      call fft_xy_parallel(fdx_r,fdx_i,linv=.true.,lneed_im=.false.)
      vx=fdx_r
      call fft_xy_parallel(fdy_r,fdy_i,linv=.true.,lneed_im=.false.)
      vy=fdy_r
!
      call fft_xy_parallel(frx_r,frx_i,linv=.true.,lneed_im=.false.)
      call fft_xy_parallel(fry_r,fry_i,linv=.true.,lneed_im=.false.)
!
    endsubroutine helmholtz
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
