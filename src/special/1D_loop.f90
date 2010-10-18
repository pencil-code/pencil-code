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
  use Cparam
  use Messages, only: svn_id, fatal_error
  use Sub, only: keep_compiler_quiet,cubic_step
!
  implicit none
!
  include '../special.h'
!
  real :: Kpara=0.,Kperp=0.,cool_RTV=0.
  real :: tau_inv_newton=0.,exp_newton=0.
  real :: init_time=0.
!
  character (len=labellen), dimension(3) :: iheattype='nothing'
  real, dimension(2) :: heat_par_exp=(/0.,1./)
  real, dimension(2) :: heat_par_exp2=(/0.,1./)
  real, dimension(3) :: heat_par_gauss=(/0.,1.,0./)

  namelist /special_run_pars/ &
      Kpara,Kperp,cool_RTV,tau_inv_newton,exp_newton,init_time, &
      iheattype,heat_par_exp,heat_par_exp2,heat_par_gauss
!
! variables for print.in
!
  integer :: idiag_dtchi2=0   ! DIAG_DOC: heatconduction
                                  ! DIAG DOC: in special module
  integer :: idiag_dtrad=0        ! DIAG_DOC: radiative loss from RTV
  integer :: idiag_dtnewt=0
!
!  miscellaneous variables
!
  real, save, dimension (mx) :: lnTT_init_prof,lnrho_init_prof
  integer, save, dimension(mseed) :: nano_seed
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
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
      real, dimension (mx) :: xtmp
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
        read(unit) xtmp
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
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
      endif
!
      if (cool_RTV/=0) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
      endif
!
      if (tau_inv_newton/=0) then
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
    subroutine special_calc_entropy(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  entropy equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      if (Kpara/=0) call calc_heatcond_spitzer(df,p)
      if (cool_RTV/=0) call calc_heat_cool_RTV(df,p)
      if (tau_inv_newton/=0) call calc_heat_cool_newton(df,p)
      if (iheattype(1)/='nothing') call calc_artif_heating(df,p)
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine calc_heatcond_spitzer(df,p)
!
!  Computes Spitzer heat conduction along the 1D loop
!  rhs = Div K T^2.5 Grad(T)
!
!  10-oct-04/bing: coded
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot2
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: chi,glnTT2,rhs
!
      chi=Kpara*exp(p%lnTT*2.5-p%lnrho)* &
          cubic_step(t,init_time,init_time)
!
      call dot2(p%glnTT,glnTT2)
!
      rhs = Kpara * exp(p%lnTT*3.5)*cubic_step(t,init_time,init_time)
      rhs = rhs * (3.5 * glnTT2 + p%del2lnTT)
!
!  Add to energy equation.
      if (ltemperature) then
        if (ltemperature_nolog) then
          call fatal_error('calc_heatcond_spitzer', &
              'not implented for ltemperature_nolog')
        else
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + rhs
        endif
      else if (lentropy) then
        if (pretend_lnTT) then
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + rhs
        else
          call fatal_error('calc_heatcond_spitzer', &
              'not implented for lentropy and not pretend_lnTT')
        endif
      endif
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chi*dxyz_2
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
    lnQ   = get_lnQ(lnTT_SI)
!
    rtv_cool = lnQ-unit_lnQ+lnneni-p%lnTT-p%lnrho
    rtv_cool = gamma*p%cp1*exp(rtv_cool)
!
    rtv_cool = rtv_cool*cool_RTV*(1.-tanh(3e4*(p%rho-1e-4)))/2.
!
    rtv_cool = rtv_cool * cubic_step(t,init_time,init_time)
!
!     add to temperature equation
!
    if (ltemperature) then
      if (ltemperature_nolog) then
        call fatal_error('calc_heat_cool_RTV', &
            'not implemented for ltemperature_nolog')
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
      newton  = exp(lnTT_init_prof(l1:l2)-p%lnTT)-1.
!
!  Multiply by density dependend time scale
      tau_inv_tmp = tau_inv_newton * exp(-exp_newton*(lnrho0-p%lnrho))
!
!  Adjust time scale by the initialization time
      tau_inv_tmp =  tau_inv_tmp * cubic_step(t,init_time,init_time)
!
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
    subroutine calc_artif_heating(df,p)
!
!  Subroutine to calculate intrisic heating.
!  Activated by setting iheattype = exp, exp2 and/or gauss
!  Maximum of 3 different possibible heating types
!  Also set the heating parameters for the chosen heating functions.
!
!  22-sept-10/Tijmen: coded
!
      use EquationOfState, only: gamma
      use Diagnostics, only: max_mn_name
      use General, only: random_number_wrapper,random_seed_wrapper, &
          normal_deviate
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: heatinput,heat_flux
      real, dimension (nx) :: x_Mm,heat_nano,rhs,height
      real, dimension (nx) :: heat_event,heat_event1D
      integer, dimension(mseed) :: global_rstate
      real :: heat_unit
      real :: nano_sigma_t,nano_time,nano_start,nano_sigma_z
      real :: nano_flare_energy,nano_pos_x,nano_pos_z,nano_pos_y
      real :: nano_amplitude
      real, dimension(2) :: event_pos
      type (pencil_case) :: p
      integer :: i
!
      heat_unit= unit_density*unit_velocity**3/unit_length
      x_Mm = x(l1:l2)*unit_length*1e-6
!
! The height in [Mm]
      height = Lxyz(1)/pi*sin((x(l1:l2)-xyz0(1)) / (xyz0(1)+Lxyz(1)) * pi)*unit_length*1e-6
!
      heatinput = 0.
      heat_flux = 0.
!
      do i=1,3
        if (headtt) print*,'iheattype:',iheattype(i)
        select case(iheattype(i))
        case ('nothing')
        case ('one-sided')
          !
          heatinput=heatinput + &
              heat_par_exp(1)*exp(-x(l1:l2)/heat_par_exp(2))/heat_unit
          heatinput=heatinput + &
              heat_par_exp2(1)*exp(-x(l1:l2)/heat_par_exp2(2))/heat_unit
          !
          heat_flux=heat_flux - heat_par_exp(1)*heat_par_exp(2)*1e6
!
          if (headtt) print*,'Flux of exp heating: ',heat_flux(1),' [ W m^(-2)]'
!
        case ('exp')
          ! heat_par_exp(1) should be 530 W/m^3 (amplitude)
          ! heat_par_exp(2) should be 0.3 Mm (scale height)
          !
          heatinput=heatinput + &
              heat_par_exp(1)*exp(-height/heat_par_exp(2))/heat_unit
!
          heat_flux=heat_flux - heat_par_exp(1)*heat_par_exp(2)*1e6
!
          if (headtt) print*,'Flux of exp heating: ',heat_flux(1),' [ W m^(-2)]'
!
        case ('exp2')
          ! A second exponential function
          ! For Sven et al. 2010 set:
          ! heat_par_exp= (10 , 0.2 )
          ! heat_par_exp2= (1e-4 , 4.)
          !
          heatinput=heatinput + &
              heat_par_exp2(1)*exp(-height/heat_par_exp2(2))/heat_unit
!
          heat_flux=heat_flux + heat_par_exp2(1)*heat_par_exp2(2)*1e6

          if (headtt) print*,'Flux for exp2 heating: ', &
              heat_par_exp2(1)*heat_par_exp2(2)*1e-6* &
              (1.-exp(-lz*unit_length*1e-6/heat_par_exp2(2)))
!
        case ('gauss')
          ! heat_par_gauss(1) is Center (z in Mm)
          ! heat_par_gauss(2) is Width (sigma)
          ! heat_par_gauss(3) is the amplitude (Flux)
          !
          heatinput=heatinput + &
              heat_par_gauss(3)*exp(-((x_Mm-heat_par_gauss(1))**2/ &
              (2*heat_par_gauss(2)**2)))/heat_unit
!
        case ('nanof')
          ! simulate nanoflare heating =)
          ! height dependend amplitude und duration
          ! idea: call random numbers , make condition when to flare, when
          ! flare get position
          ! then over time release a flare in shape of gaussian. Also
          ! distribution is gaussian around a point.
          ! gaussian distribution is gained via the dierfc function (double
          ! prec, inverser erf function)
          ! we draw for new flares as long as nano_time is /= 0. when done we
          ! reset nano_time to 0. again.
          !
          ! Later we will implement more options for nanoflaring! =D
          !
          ! SAVE GLOBAL SEED
          ! LOAD NANO_SEED
          call random_seed_wrapper(GET=global_rstate)
          call random_seed_wrapper(PUT=nano_seed)

          if (nano_start .eq. 0.) then
            ! If 0 roll to see if a nanoflare occurs
            call random_number_wrapper(nano_start)
            !
            if (nano_start .gt. 0.95 ) then
              ! 5% chance for a nanoflare to occur, then get the location.
              call normal_deviate(nano_pos_z)
              nano_pos_z=nano_pos_z*lz
              call random_number_wrapper(nano_pos_y)
              call random_number_wrapper(nano_pos_x)
              nano_time=60.
            else
              ! else no nanoflare, reset nano_start to 0. for the next
              ! timestep
              nano_start=0.
            endif
          endif
          !
          if (nano_start .ne. 0.) then
            ! if nano_start is not 0. then there is a nanoflare!
            ! 2nd assumption , nanoflare takes 60 seconds =)
            nano_flare_energy=10.d17 !joules
            nano_sigma_z=0.5
            nano_sigma_t=2.5
!
            nano_amplitude=nano_flare_energy/(pi/2*nano_sigma_t*nano_sigma_z*1.d6 )
!
            heat_nano=nano_amplitude*exp(-((nano_time-5.))**2/( 2*2.**2))* &
                exp(-((x_Mm-nano_pos_z)**2/ (2*0.5**2)))
            nano_time=nano_time-dt*unit_time
!
            heatinput=heatinput + heat_nano/heat_unit
!
            if (nano_time .le. 0.) nano_start=0.
          end if
          !
          !SAVE NANO_SEED
          !RESTORE GLOBAL SEED
          call random_seed_wrapper(GET=nano_seed)
          call random_seed_wrapper(PUT=global_rstate)
!
! amp_t = exp(-((t-nano_time)**2/(2.*nano_dur**2)))
!-> no idea how to properly implement
! spread = exp(-((x_Mm-nano_pos)**2/(2.*nano_spread**2)))
!
! issue! How to put in timelike guassian for a random event starting
! at a random time?
!
        case ('event')
          ! one small point heating event (gaussian to prevent too strong gradients)
          ! one point in time, one point in space!
          if (t*unit_time .gt. 150. .AND. t*unit_time .lt. 1000.) then
            event_pos(1)=7.5
            event_pos(2)=15.
            heat_event=10.*exp(-((250.-t*unit_time))**2/(2*(20.*unit_time)**2))* &
                exp(-((x_Mm-event_pos(1))**2/ (2*0.2**2))) * &
                exp(-((x_Mm-event_pos(2))**2/ (2*0.2**2)))
            heatinput=heatinput + heat_event/heat_unit
          endif
!
        case ('event1D')
          ! one small point heating event (gaussian to prevent too strong gradients)
          ! one point in time, one point in space!
          if (t*unit_time .gt. 300. .AND. t*unit_time .lt. 10000.) then
            if (t*unit_time .gt. 300. .AND. t*unit_time .lt. 301.) &
                print*,'EVENTTTT!!!!!'
            event_pos(1)=10.
            heat_event1D=10.*exp(-((400.-t))**2/( 2*50.**2))* &
                exp(-((x_Mm-event_pos(1))**2/ (2*0.2**2)))
            heatinput=heatinput + heat_event1D/heat_unit
          endif
!
        case default
          if (headtt) call fatal_error('calc_artif_heating', &
              'Please provide correct iheattype')
        endselect
      enddo
      !
      if (headtt) print*,'Total flux for all types:',heat_flux(1)
!
! Add to energy equation
!
      rhs = p%TT1*p%rho1*gamma*p%cp1*heatinput* &
          cubic_step(t,init_time,init_time)
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
!
      if (lfirst.and.ldt) then
        if (ldiagnos.and.idiag_dtnewt/=0) then
          itype_name(idiag_dtnewt)=ilabel_max_dt
          call max_mn_name(rhs/cdts,idiag_dtnewt,l_dt=.true.)
        endif
        dt1_max=max(dt1_max,rhs/cdts)
      endif
!
    endsubroutine calc_artif_heating
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include 'special_dummies.inc'
!********************************************************************
endmodule Special
