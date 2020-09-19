! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
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
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
  real :: Kpara=0.,Kperp=0.,Ksat=0.,Kc=0.,cool_RTV=0.,Kchrom=0.
  real :: exp_RTV=0.,cubic_RTV=0.,tanh_RTV=0.
  real :: tau_inv_newton=0.,exp_newton=0.
  real :: tanh_newton=0.,cubic_newton=0.
  real :: width_newton=0.,width_RTV=0.,hyper3_chi=0.
  real :: init_time=0.,lnTT0_chrom=0.,width_lnTT_chrom=0.
  real :: hcond_grad_iso=0.,chi_aniso=0.
  real :: tau_inv_spitzer=0.
  real :: bullets_x0=0.,bullets_dx=0.
  real :: bullets_t0=0.,bullets_dt=0.
  real :: bullets_h0=0.,hyper3_diffrho=0.
  logical :: lfilter_farray=.false.,lrad_loss=.false.
  real, dimension(mvar) :: filter_strength=0.01
  real :: Ltot,R_hyper3=0.
  real, dimension (nx) :: diffus_chi
!
  character (len=labellen), dimension(3) :: iheattype='nothing'
  character (len=labellen) :: loop_frac='full'
  real, dimension(2) :: heat_par_exp=(/0.,1./)
  real, dimension(2) :: heat_par_exp2=(/0.,1./)
  real, dimension(3) :: heat_par_gauss=(/0.,1.,0./)
  real, dimension(3) :: heat_par_vandoors=(/1.,0.1,50./)
!
  namelist /special_run_pars/ &
      Kpara,Kperp,cool_RTV,tau_inv_newton,exp_newton,init_time, &
      iheattype,heat_par_exp,heat_par_exp2,heat_par_gauss, &
      width_newton,tanh_newton,cubic_newton,Kchrom, &
      lnTT0_chrom,width_lnTT_chrom,width_RTV, &
      exp_RTV,cubic_RTV,tanh_RTV,hcond_grad_iso,Ksat,Kc, &
      tau_inv_spitzer,hyper3_chi,bullets_x0,bullets_dx, &
      bullets_t0,bullets_dt,bullets_h0,hyper3_diffrho, &
      lfilter_farray,filter_strength,loop_frac,&
      heat_par_vandoors,R_hyper3,lrad_loss,chi_aniso
!
! variables for print.in
!
  integer :: idiag_dtchi2=0       ! DIAG_DOC: heatconduction
                                  ! DIAG DOC: in special module
  integer :: idiag_dthyper3=0
  integer :: idiag_dtrad=0        ! DIAG_DOC: radiative loss from RTV
  integer :: idiag_dtspitzer=0    ! DIAG_DOC: Spitzer heat conduction
                                  ! DIAG_DOC: time step
  integer :: idiag_dtnewt=0
  integer :: idiag_qmax=0     ! DIAG_DOC: max of heat flux vector
  integer :: idiag_qrms=0     ! DIAG_DOC: rms of heat flux vector
!
! variables for video.in
!
  integer :: ivid_logQ=0; ivid_rtv=0; ivid_spitzer=0

  real, target, dimension(:,:), allocatable :: rtv_xy, rtv_xy2, rtv_xy3, rtv_xy4
  real, target, dimension(:,:), allocatable :: rtv_xz, rtv_yz, rtv_xz2
  real, target, dimension(:,:), allocatable :: logQ_xy, logQ_xy2, logQ_xy3, logQ_xy4
  real, target, dimension(:,:), allocatable :: logQ_xz, logQ_xz2, logQ_yz
  real, target, dimension(:,:), allocatable :: spitzer_xy, spitzer_xy2, spitzer_xy3, spitzer_xy4
  real, target, dimension(:,:), allocatable :: spitzer_xz, spitzer_xz2, spitzer_yz
!
!  miscellaneous variables
!
  real, save, dimension (mx) :: lnTT_init_prof,lnrho_init_prof
  integer, save, dimension(mseed) :: nano_seed
!
  real :: Kspitzer_para_SI = 2e-11, Kspitzer_para=0.
  real :: Ksaturation_SI = 7e7,Ksaturation=0.,ln_unit_TT=0.
!
  contains
!***********************************************************************
    subroutine register_special()
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  13-sep-10/bing: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: xtmp
      character (len=*), parameter :: filename='/strat.dat'
      integer :: lend,unit=12
      real :: dummy=1.
!
      ln_unit_TT = alog(real(unit_temperature))
!
      Kspitzer_para = Kspitzer_para_SI /unit_density/unit_velocity**3./ &
          unit_length*unit_temperature**(3.5)
!
      Ksaturation = Ksaturation_SI /unit_velocity**3.*unit_temperature**1.5
!
      call keep_compiler_quiet(f)
!
      inquire(IOLENGTH=lend) dummy
!
      if (lrun .and. (tau_inv_newton/=0.)) then
        open(unit,file=trim(directory_snap)//filename, &
            form='unformatted',status='unknown',recl=lend*mx)
        read(unit) xtmp
        read(unit) lnrho_init_prof
        read(unit) lnTT_init_prof
        close(unit)
      endif
!
      select case (loop_frac)
      case ('full')
        Ltot = Lxyz(1) + 2*xyz0(1)
      case ('half')
        Ltot = 2*Lxyz(1) + 2*xyz0(1)
      case default
        call fatal_error('initialize_initial_condition','Wrong loop_frac')
        Ltot = 0.
      endselect
!
      if (ivid_rtv/=0) then
        !call alloc_slice_buffers(rtv_xy,rtv_xz,rtv_yz,rtv_xy2,rtv_xy3,rtv_xy4,rtv_xz2)
        if (lwrite_slice_xy .and..not.allocated(rtv_xy) ) allocate(rtv_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(rtv_xz) ) allocate(rtv_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(rtv_yz) ) allocate(rtv_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(rtv_xy2)) allocate(rtv_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(rtv_xy3)) allocate(rtv_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(rtv_xy4)) allocate(rtv_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(rtv_xz2)) allocate(rtv_xz2(nx,nz))
      endif
      if (ivid_logQ/=0) then
        !call alloc_slice_buffers(logQ_xy,logQ_xz,logQ_yz,logQ_xy2,logQ_xy3,logQ_xy4,logQ_xz2)
        if (lwrite_slice_xy .and..not.allocated(logQ_xy) ) allocate(logQ_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(logQ_xz) ) allocate(logQ_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(logQ_yz) ) allocate(logQ_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(logQ_xy2)) allocate(logQ_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(logQ_xy3)) allocate(logQ_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(logQ_xy4)) allocate(logQ_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(logQ_xz2)) allocate(logQ_xz2(nx,nz))
      endif
      if (ivid_spitzer/=0) then
        if (lwrite_slice_xy .and..not.allocated(spitzer_xy) ) allocate(spitzer_xy (nx,ny))
        if (lwrite_slice_xz .and..not.allocated(spitzer_xz) ) allocate(spitzer_xz (nx,nz))
        if (lwrite_slice_yz .and..not.allocated(spitzer_yz) ) allocate(spitzer_yz (ny,nz))
        if (lwrite_slice_xy2.and..not.allocated(spitzer_xy2)) allocate(spitzer_xy2(nx,ny))
        if (lwrite_slice_xy3.and..not.allocated(spitzer_xy3)) allocate(spitzer_xy3(nx,ny))
        if (lwrite_slice_xy4.and..not.allocated(spitzer_xy4)) allocate(spitzer_xy4(nx,ny))
        if (lwrite_slice_xz2.and..not.allocated(spitzer_xz2)) allocate(spitzer_xz2(nx,nz))
      endif
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  04-sep-10/bing: coded
!
      if (tau_inv_spitzer/=0.) then
        lpenc_requested(i_uglnrho)=.true.
        lpenc_requested(i_divu)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
      endif
!
      if (chi_aniso/=0.) then
         lpenc_requested(i_bij)=.true.
         lpenc_requested(i_bunit)=.true.
         lpenc_requested(i_glnTT)=.true.
         lpenc_requested(i_hlnTT)=.true.
         lpenc_requested(i_glnrho)=.true.
      endif
!
      if (Kpara/=0.) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
!
      if (hcond_grad_iso/=0.) then
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (Kchrom/=0.) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
!
      if (cool_RTV/=0.) then
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_rho)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
      endif
!
      if ((tau_inv_newton/=0.).or.(bullets_h0/=0.)) then
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
      endif
!
      if (iheattype(1)/='nothing') then
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_cp1)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (Ksat /= 0.) then
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if (R_hyper3 /= 0.) lpenc_requested(i_u2)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
      logical, dimension (npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!  04-sep-10/bing: coded
!
      use Diagnostics, only: parse_name
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
        idiag_dtchi2=0.
        idiag_dthyper3=0.
        idiag_dtrad=0.
        idiag_dtnewt=0
        idiag_dtspitzer=0
        idiag_qmax=0
        idiag_qrms=0
        ivid_logQ=0; ivid_rtv=0; ivid_spitzer=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtchi2',idiag_dtchi2)
        call parse_name(iname,cname(iname),cform(iname),'dthyper3',idiag_dthyper3)
        call parse_name(iname,cname(iname),cform(iname),'dtrad',idiag_dtrad)
        call parse_name(iname,cname(iname),cform(iname),'dtnewt',idiag_dtnewt)
        call parse_name(iname,cname(iname),cform(iname),'dtspitzer',idiag_dtspitzer)
        call parse_name(iname,cname(iname),cform(iname),'qmax',idiag_qmax)
        call parse_name(iname,cname(iname),cform(iname),'qrms',idiag_qrms)
      enddo
!
!  check for those quantities for which we want video slices
!
      do iname = 1,nnamev
        call parse_name(iname,cnamev(iname),cformv(iname),'logQ',ivid_logQ)
        call parse_name(iname,cnamev(iname),cformv(iname),'rtv',ivid_rtv)
        call parse_name(iname,cnamev(iname),cformv(iname),'spitzer',ivid_spitzer)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        call farray_index_append('i_dtchi2',idiag_dtchi2)
        call farray_index_append('i_dthyper3',idiag_dthyper3)
        call farray_index_append('i_dtrad',idiag_dtrad)
        call farray_index_append('i_dtnewt',idiag_dtnewt)
        call farray_index_append('i_dtspitzer',idiag_dtspitzer)
        call farray_index_append('i_qmax',idiag_qmax)
        call farray_index_append('i_qrms',idiag_qrms)
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  12-may-11/bingert: coded
!
      use EquationOfState, only: gamma
      use Deriv, only: der6
      use Diagnostics,     only : max_mn_name, sum_mn_name
      use Sub, only: identify_bcs, multsv, dot, del6
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: rhs,tmp,hyper3_coeff,hc
      real, dimension (nx,3) :: K1
      integer :: i,j
!
      intent(in) :: p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
      if (lfilter_farray) call filter_farray(f,df)
!
!  Add hyper diffusion with fixed Reynoldsnmumber
!
      if (R_hyper3 /= 0.) then
        hyper3_coeff = sqrt(p%u2)/dxmax_pencil/R_hyper3
        call del6(f, ilnTT, hc, IGNOREDX=.true.)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + hyper3_coeff * hc
        call del6(f, ilnrho, hc, IGNOREDX=.true.)
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + hyper3_coeff * hc
        call del6(f, iux, hc, IGNOREDX=.true.)
        df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + hyper3_coeff * hc
        call del6(f, iuy, hc, IGNOREDX=.true.)
        df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + hyper3_coeff * hc
        call del6(f, iuz, hc, IGNOREDX=.true.)
        df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + hyper3_coeff * hc
!
        if (lfirst.and.ldt) dt1_max=max(dt1_max,tmp/0.01)
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  26-jun-06/tony: dummy
!
      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
      case ('rtv'); call assign_slices_scal(slices,rtv_xy,rtv_xz,rtv_yz,rtv_xy2, &
                                            rtv_xy3,rtv_xy4,rtv_xz2)
!
      case ('logQ'); call assign_slices_scal(slices,logQ_xy,logQ_xz,logQ_yz,logQ_xy2, &
                                             logQ_xy3,logQ_xy4,logQ_xz2)
!
      case ('spitzer')
        call assign_slices_scal(slices,spitzer_xy,spitzer_xz,spitzer_yz,spitzer_xy2, &
                                spitzer_xy3,spitzer_xy4,spitzer_xz2)
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_,llast)
!
!  10-oct-12/bing: coded
!
      use EquationOfState, only: gamma
!
      logical, intent(in) :: llast
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_
      integer,save :: j=35
      real :: lnTT_SI,dt_step,dt_rad,ln_coeff,one_m_alpha,lnTT_res
      integer :: l
      logical :: notdone
!
      real, parameter, dimension (37) :: intlnT = (/ &
          8.74982, 8.86495, 8.98008, 9.09521, 9.21034, 9.44060, 9.67086 &
          , 9.90112, 10.1314, 10.2465, 10.3616, 10.5919, 10.8221, 11.0524 &
          , 11.2827, 11.5129, 11.7432, 11.9734, 12.2037, 12.4340, 12.6642 &
          , 12.8945, 13.1247, 13.3550, 13.5853, 13.8155, 14.0458, 14.2760 &
          , 14.5063, 14.6214, 14.7365, 14.8517, 14.9668, 15.1971, 15.4273 &
          ,  15.6576,  69.0776 /)
      real, parameter, dimension (36) :: alpha_rad = (/ &
          67.4287,   40.0384,  21.3332, 20.0001,  9.33333,  4.66647 &
          ,  2.33345,-0.566705,-0.199849, 0.199849, 0.899692,  1.33405 &
          ,  1.66611, 0.833237,-0.166805, -1.49977,  0.00000, 0.566656 &
          , 0.233155,-0.966534,-0.633508, 0.800159, 0.433348, -0.0998807 &
          ,-0.499988, -1.00000, -1.96697, -3.06643, -3.60040, -5.80186 &
          , -1.86549, -2.66724,-0.166734, 0.566900, 0.666504,  6.23231 /)
      real, parameter, dimension (36) :: lnchi_rad = (/ &
          -690.935, -448.121, -280.147, -268.022, -169.777, -125.719 &
          , -103.157, -74.4422, -78.1590, -82.2545, -89.5060, -94.1066 &
          , -97.7002, -88.4950, -77.2118, -61.8654, -79.4776, -86.2624 &
          , -82.1925, -67.2755, -71.4930, -89.9794, -85.1652, -78.0439 &
          , -72.6083, -65.7003, -52.1185, -36.4226, -28.6767,  3.51167 &
          , -54.4966, -42.5892, -80.0138, -91.1629, -92.6996, -179.847 /)
!
      if (llast .and. lrad_loss) then
!
        do l=l1,l2
          do m=m1,m2
            do n=n1,n2
              notdone=.true.
              dt_rad = dt*unit_time
              lnTT_SI = f(l,m,n,ilnTT) + ln_unit_TT
              do while (notdone)
!
                if (lnTT_SI > intlnT(j) .and. &
                    lnTT_SI <=  intlnT(j+1) ) then
! we are the propper intervall, compute the factors
!
                  one_m_alpha = (1-alpha_rad(j))
!
                  ln_coeff = f(l,m,n,ilnrho) + alog(real(1.5**2./8.7*(gamma-1.)/0.6/m_p/k_B))+ lnchi_rad(j)
                  ln_coeff = ln_coeff + alog(real(unit_temperature/unit_velocity**2./unit_density/unit_length**6.))
!
! compute resulting temperature in SI units
!
! first check if we would get negative temperatues
                  if (-one_m_alpha*dt_rad*exp(ln_coeff)+exp(one_m_alpha*lnTT_SI) < 0) then
                    lnTT_res = intlnT(j)*0.9
                  else
                    lnTT_res =1./one_m_alpha*alog(  &
                        -one_m_alpha*dt_rad*exp(ln_coeff) &
                        +exp(one_m_alpha*lnTT_SI))
                  endif
!
                  if (lnTT_res >= intlnT(j)) then
                    ! everything is fine
                    f(l,m,n,ilnTT) = lnTT_res - ln_unit_TT
                    notdone = .false.
                  else
                    ! timestep to large, we need to repeat for another intervall
                    dt_step = exp(-ln_coeff)/one_m_alpha*  &
                        (exp(one_m_alpha*lnTT_SI) - exp(one_m_alpha*intlnT(j)))
                    dt_rad = dt_rad - dt_step
                    f(l,m,n,ilnTT)  = intlnT(j)-ln_unit_TT
                    lnTT_SI = f(l,m,n,ilnTT) + ln_unit_TT
                    j = j-1
                    if (j<=0) notdone=.false.
                  endif
                elseif (lnTT_SI == intlnT(1)) then
                  notdone=.false.
                else
                  j = j + sign(1.,lnTT_SI-intlnT(j))
                  if (j <= 0) then
                    j=1
                    notdone=.false.
                  elseif (j >= 37) then
                    j=36
                    notdone=.false.
                  endif
                endif
              enddo  ! while loop
            enddo
          enddo
        enddo
      endif  ! if (llast)
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dt_)
      call keep_compiler_quiet(llast)
!
    endsubroutine  special_after_timestep
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  entropy equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      use Sub, only: del6
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: hc
      integer :: itmp
!
      diffus_chi=0.
!
      if (hcond_grad_iso/=0.) call calc_heatcond_glnTT_iso(df,p)
      if (Kpara/=0.) call calc_heatcond_spitzer(f,df,p)
      if (cool_RTV/=0.) call calc_heat_cool_RTV(df,p)
      if (tau_inv_newton/=0.) call calc_heat_cool_newton(df,p)
      if (iheattype(1)/='nothing') call calc_artif_heating(df,p)
      if (Kchrom/=0.) call calc_heatcond_kchrom(df,p)
      if (bullets_h0/=0.) call calc_bullets_heating(df,p)
      if (chi_aniso/=0.) call calc_heatcond_constchi(df,p)
!
      if (ltemperature) then
        itmp = ilnTT
      else if (lentropy) then
        itmp = iss
      else
        call fatal_error('hyper3_chi special','only for ltemperature')
      endif
!
      if (hyper3_chi /= 0.) then
        call del6(f,itmp,hc,IGNOREDX=.true.)
          df(l1:l2,m,n,itmp) = df(l1:l2,m,n,itmp) + hyper3_chi*hc
!
!  due to ignoredx hyper3_chi has [1/s]
!
          if (lfirst.and.ldt) dt1_max=max(dt1_max,hyper3_chi/0.01)
      endif
!
      if (lfirst.and.ldt) then
         maxdiffus=max(maxdiffus,diffus_chi)
      endif
!
    endsubroutine special_calc_energy
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
      use Sub, only: del6
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: hc
!
      if (hyper3_diffrho /= 0.) then
        if (ldensity_nolog.and.(ilnrho /= irho)) &
            call fatal_error('hyper3_diffrho special','please check')
!
        call del6(f,ilnrho,hc,IGNOREDX=.true.)
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + hyper3_diffrho*hc
!
!  due to ignoredx hyper3_diffrho has [1/s]
!
        if (lfirst.and.ldt) dt1_max=max(dt1_max,hyper3_diffrho/0.01)
      endif
!
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine calc_heatcond_spitzer(f,df,p)
!
!  Computes Spitzer heat conduction along the 1D loop
!  rhs = Div K T^2.5 Grad(T)
!
!  10-oct-04/bing: coded
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot2,dot,cubic_step
      use Slices_methods, only: store_slices
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: chi,glnTT2,rhs,u_spitzer
      real, dimension (nx) :: chi_sat,chi_c,gKc,gKpara,gKsat !,glnTT2_upwind
      real, dimension (nx,3) :: tmpv2,tmpv !,glnTT_upwind
      real :: Ksatb,Kcb
      integer :: i,j
!
      chi=Kpara*exp(p%lnTT*2.5-p%lnrho)*p%cp1
!
      if (init_time /= 0.) chi=chi*cubic_step(real(t),init_time,init_time)
!
!      do i=1,3
!        call der_upwind(f,-p%glnTT,ilnTT,glnTT_upwind(:,i),i)
!      enddo
!      call dot2(glnTT_upwind,glnTT2_upwind)
      call dot2(p%glnTT,glnTT2)
!
      gKpara = 3.5 * glnTT2 !_upwind
!
      if (Ksat /= 0.) then
        Ksatb = Ksat*7.28d7 /unit_velocity**3. * unit_temperature**1.5
        chi_sat =  Ksatb * sqrt(p%TT/max(tini,glnTT2))*p%cp1
        tmpv(:,:)=0.
        do i=1,3
          do j=1,3
            tmpv(:,i)=tmpv(:,i)+p%glnTT(:,j)*p%hlnTT(:,j,i)
          enddo
        enddo
        do i=1,3
          tmpv2(:,i) = p%glnrho(:,1) + 1.5*p%glnTT(:,1) - tmpv(:,1)/max(tini,glnTT2)
        enddo
        call dot(tmpv2,p%glnTT,gKsat)
        where (chi > chi_sat)
          chi = chi_sat
          gKpara = gKsat
        endwhere
      endif
!
! limit the diffusion speed to the speed of light
!
      if (Kc /= 0.) then
        Kcb = Kc * dxmin * c_light * cdtv
        chi_c = Kcb
        call dot(p%glnTT+p%glnrho,p%glnTT,gKc)
        where (chi > chi_c)
          chi = chi_c
          gKpara = gKc
        endwhere
      endif
!
      rhs = gamma * chi * (gKpara + p%del2lnTT)
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
          df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) + rhs/(gamma*p%cp1)
        endif
      else
        call fatal_error('calc_heatcond_spitzer', &
            'not implented for lentropy and not pretend_lnTT')
      endif
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chi*dxyz_2
!
        u_spitzer = 7./2.*gamma*chi*sum(abs(p%glnTT)*dline_1,2)
!
! MR: dt1_max should not be used directly here
!
        dt1_max=max(dt1_max,u_spitzer/cdt)
!
        if (ldiagnos.and.idiag_dtspitzer/=0) then
          call max_mn_name(max(diffus_chi/cdtv,u_spitzer/cdt), &
              idiag_dtspitzer,l_dt=.true.)
        endif
      endif
!
      if (lvideo) then
!
! slices
!
        if (ivid_spitzer/=0) &
          call store_slices(1./chi,spitzer_xy,spitzer_xz,spitzer_yz,spitzer_xy2,spitzer_xy3,spitzer_xy4,spitzer_xz2)
      endif
!
    endsubroutine calc_heatcond_spitzer
!***********************************************************************
    subroutine calc_heatcond_kchrom(df,p)
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: cubic_step, cubic_der_step, dot2,cross,dot
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension(nx) :: K_chrom,glnTT_2
      real, dimension(nx) :: gK_chrom,rhs,chix,tmp
      real, dimension(nx,3) :: tmpv,tmpv2
      integer :: i,j
!
      call dot2(p%glnTT,glnTT_2)
!
      K_chrom=Kchrom*(1.-cubic_step(glnTT_2,lnTT0_chrom,width_lnTT_chrom))
      chix=p%rho1*K_chrom*p%cp1
!
      tmpv(:,1) = p%hlnTT(:,2,3)-p%hlnTT(:,3,2)
      tmpv(:,2) = p%hlnTT(:,3,1)-p%hlnTT(:,1,3)
      tmpv(:,3) = p%hlnTT(:,1,2)-p%hlnTT(:,2,1)
!
      call cross(tmpv,p%glnTT,tmpv2)
!
      tmpv(:,:)=0.
      do i=1,3
        do j=1,3
          tmpv(:,i) = tmpv(:,i) + p%glnTT(:,j)*p%hlnTT(:,i,j)
        enddo
        tmpv(:,i)=tmpv(:,i) + tmpv2(:,i)
      enddo
!
      call dot(2*tmpv,p%glnTT,tmp)
!
      gK_chrom= - Kchrom*cubic_der_step(glnTT_2,lnTT0_chrom,width_lnTT_chrom)
!
      rhs = tmp*gK_chrom + glnTT_2*K_chrom + K_chrom*p%del2lnTT
!
      rhs = p%rho1*p%cp1*rhs
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + rhs
!
!  Check maximum diffusion from thermal diffusion.
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi2/=0.) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_kchrom
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
    use Messages, only: warning
    use Sub, only: cubic_step
    use Slices_methods, only: store_slices
!
    real, dimension (mx,my,mz,mvar), intent(inout) :: df
    type (pencil_case), intent(in) :: p
!
    real, dimension (nx) :: lnQ,rtv_cool,lnTT_SI,lnneni,delta_lnTT
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
    call getlnQ(lnTT_SI,lnQ,delta_lnTT)
!
    rtv_cool = lnQ-unit_lnQ+lnneni-p%lnTT-p%lnrho
    rtv_cool = gamma*p%cp1*exp(rtv_cool)
!
    if (exp_RTV/=0.) then
      call warning('cool_RTV','exp_RTV not yet implemented')
    elseif (tanh_RTV/=0.) then
      rtv_cool=rtv_cool*cool_RTV* &
          0.5*(1.-tanh(width_RTV*(p%lnrho-tanh_RTV)))
!
    elseif (cubic_RTV/=0.) then
      rtv_cool=rtv_cool*cool_RTV* &
          (1.-cubic_step(p%lnrho,cubic_RTV,width_RTV))
!
    endif
!
    if (init_time /= 0.) rtv_cool = rtv_cool * cubic_step(real(t),init_time,init_time)
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
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss)-rtv_cool/(gamma*p%cp1)
      endif
    else
      call fatal_error('calc_heat_cool_RTV','not implemented for lentropy')
    endif
!
    if (lfirst.and.ldt) then
      dt1_max=max(dt1_max,rtv_cool/cdts)
      dt1_max=max(dt1_max,abs(rtv_cool/max(tini,delta_lnTT)))
      if (ldiagnos.and.idiag_dtrad/=0) then
        itype_name(idiag_dtrad)=ilabel_max_dt
        call max_mn_name(rtv_cool/cdts,idiag_dtrad,l_dt=.true.)
        call max_mn_name(abs(rtv_cool/max(tini,delta_lnTT)),idiag_dtrad,l_dt=.true.)
      endif
    endif
!
! fill video slices
!
    if (lvideo) then
      if (ivid_rtv/=0) &
        call store_slices(rtv_cool,rtv_xy,rtv_xz,rtv_yz,rtv_xy2,rtv_xy3,rtv_xy4,rtv_xz2)
      if (ivid_logQ/=0) &
        call store_slices(lnQ*0.43429448,logQ_xy,logQ_xz,logQ_yz,logQ_xy2,logQ_xy3,logQ_xy4,logQ_xz2)
    endif
!
  endsubroutine calc_heat_cool_RTV
!***********************************************************************
  subroutine getlnQ(lnTT,get_lnQ,delta_lnTT)
!
!  input: lnTT in SI units
!  output: lnP  [p]= W * m^3
!
    real, parameter, dimension (36) :: intlnT = (/ &
        8.98008,  9.09521, 9.21034, 9.44060,  9.67086, 9.90112, &
        10.1314, 10.2465, 10.3616, 10.5919, 10.8221, 11.0524, &
        11.2827, 11.5129, 11.7432, 11.9734, 12.2037, 12.4340,&
        12.6642, 12.8945, 13.1247, 13.3550, 13.5853, 13.8155, &
        14.0458, 14.2760, 14.5063, 14.6214, 14.7365, 14.8517, &
        14.9668, 15.1971, 15.4273, 15.6576, 15.8878, 16.1181 /)
    real, parameter, dimension (36) :: intlnQ = (/ &
        -83.9292, -82.8931, -82.4172, -81.2275, -80.5291, -80.0532, &
        -80.1837, -80.2067, -80.1837, -79.9765, -79.6694, -79.2857,&
        -79.0938, -79.1322, -79.4776, -79.4776, -79.3471, -79.2934,&
        -79.5159, -79.6618, -79.4776, -79.3778, -79.4008, -79.5159,&
        -79.7462, -80.1990, -80.9052, -81.3196, -81.9874, -82.2023,&
        -82.5093, -82.5477, -82.4172, -82.2637, -82.1793 , -82.2023 /)
!
    real, dimension (nx), intent(in) :: lnTT
    real, dimension (nx), intent (out) :: get_lnQ,delta_lnTT
    real :: slope,ordinate
    integer :: i,j=18
    logical :: notdone
!
    get_lnQ=-1000.
    delta_lnTT = 100.
!
    do i=1,nx
      notdone=.true.
!
      do while (notdone)
!
        if (lnTT(i) >= intlnT(j) .and. lnTT(i) < intlnT(j+1)) then
!
! define slope and ordinate for linear interpolation
!
          slope=(intlnQ(j+1)-intlnQ(j))/(intlnT(j+1)-intlnT(j))
          ordinate=intlnQ(j) - slope*intlnT(j)
!
          get_lnQ(i) = slope*lnTT(i) + ordinate
          delta_lnTT(i) = intlnT(j+1) - intlnT(j)
          notdone = .false.
        else
          j = j + sign(1.,lnTT(i)-intlnT(j))
          if (j <= 0) then
            j=1
            notdone=.false.
          elseif (j >= 35) then
            j = 35
            !call fatal_error('get_lnQ','lnTT to large')
            notdone=.false.
          endif
        endif
      enddo
    enddo
!
  endsubroutine getlnQ
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
      use Sub, only: cubic_step
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: newton,tau_inv_tmp
!
      if (headtt) print*,'special_calc_energy: newton cooling',tau_inv_newton
!
!  Get reference temperature
      newton  = exp(lnTT_init_prof(l1:l2)-p%lnTT)-1.
!
!  Multiply by density dependend time scale
      if (exp_newton/=0.) then
        tau_inv_tmp = tau_inv_newton * &
            exp(-exp_newton*(lnrho0-p%lnrho))
!
      elseif (tanh_newton/=0.) then
        tau_inv_tmp = tau_inv_newton * &
            0.5*(1+tanh(width_newton*(p%lnrho-tanh_newton)))
!
      elseif (cubic_newton/=0.) then
        tau_inv_tmp = tau_inv_newton * &
            cubic_step(p%lnrho,cubic_newton,width_newton)
!
      endif
!
!  Adjust time scale by the initialization time
      if (init_time /= 0.) tau_inv_tmp =  tau_inv_tmp * cubic_step(real(t),init_time,init_time)
!
      newton  = newton * tau_inv_tmp
!
!  Add newton cooling term to entropy
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + newton
!
      if (lfirst.and.ldt) then
        dt1_max=max(dt1_max,tau_inv_tmp/cdts)
        if (ldiagnos.and.idiag_dtnewt/=0.) then
          call max_mn_name(tau_inv_tmp/cdts,idiag_dtnewt,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heat_cool_newton
!***********************************************************************
    subroutine calc_artif_heating(df,p)
!
!  Subroutine to calculate intrinsic heating.
!  Activated by setting iheattype = exp, exp2 and/or gauss
!  Maximum of 3 different possibible heating types
!  Also set the heating parameters for the chosen heating functions.
!
!  22-sept-10/Tijmen: coded
!
      use EquationOfState, only: gamma
      use Diagnostics, only: max_mn_name
      use General, only: random_number_wrapper,random_seed_wrapper, &
          normal_deviate,notanumber
      use Sub, only: cubic_step, dot2
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: heatinput,doors_heat,doors_a
      real, dimension (nx) :: x_Mm,heat_nano,rhs,height
      real, dimension (nx) :: heat_event,heat_event1D,Hlength
      real, dimension (nx) :: Htemp
      real, dimension (15,2) :: a_arr
      integer, dimension(mseed) :: global_rstate
      real :: heat_unit,heat_flux
      real :: nano_sigma_t,nano_time,nano_start,nano_sigma_z
      real :: nano_flare_energy,nano_pos_x,nano_pos_z,nano_pos_y
      real :: nano_amplitude
      real :: slope
      real, dimension(2) :: event_pos
      type (pencil_case) :: p
      integer :: i,j,k
      logical :: notdone
!
 ! the unit W/m^3 is
      heat_unit= unit_density*unit_velocity**3/unit_length
!
! coordinate along the loop in [Mm]
      x_Mm = x(l1:l2)*unit_length*1e-6
!
! The height in [Mm]
      height = Ltot/pi*sin(x(l1:l2)/Ltot*pi)*unit_length*1e-6
!
      heatinput = 0.
      heat_flux = 0.
!
      do i=1,3
        if (headtt) print*,'iheattype:',iheattype(i)
        select case(iheattype(i))
        case ('nothing')
! do nothing
        case ('one-sided')
          heatinput=heatinput + &
              heat_par_exp(1)*exp(-x(l1:l2)/heat_par_exp(2))/heat_unit
!
          heatinput=heatinput + &
              heat_par_exp2(1)*exp(-x(l1:l2)/heat_par_exp2(2))/heat_unit
!
          heat_flux=heat_flux - heat_par_exp(1)*heat_par_exp(2)*1e6
!
          if (headtt) print*,'Flux of exp heating: ',heat_flux,' [ W m^(-2)]'
!
        case ('exp')
          heatinput=heatinput + &
              heat_par_exp(1)/heat_unit*exp(-height/heat_par_exp(2))
!
!  add the flux at a ref height of 4 Mm
!
          if (headtt) then
            heat_flux = heat_par_exp(1)*heat_par_exp(2)*1e6*exp(-4./heat_par_exp(2))
            write (*,*) 'Exp heating, scale height: ',heat_par_exp(2),' [ Mm ]'
            write (*,*) 'Exp heating, amplitude',heat_par_exp(1)/heat_unit,' [ Wm^-3 ]'
            write (*,*) 'Exp heating, flux at z=4 Mm: ', heat_flux
            write (*,*) 'Exp heating, flux at z=0 Mm: ', heat_par_exp(1)*heat_par_exp(2)*1e6
          endif
!
        case ('exp2')
          ! A second exponential function
          ! For Sven et al. 2010 set:
          ! heat_par_exp= (10 , 0.2 )
          ! heat_par_exp2= (1e-4 , 4.)
          !
          heatinput=heatinput + &
              heat_par_exp2(1)/heat_unit*exp(-height/heat_par_exp2(2))
!
!  add the flux at a ref height of 4 Mm
!
          if (headtt) then
            heat_flux = heat_par_exp2(1)*heat_par_exp2(2)*1e6*exp(-4./heat_par_exp2(2))
            write (*,*) 'Exp heating, scale height: ',heat_par_exp(2),' [ Mm ]'
            write (*,*) 'Exp heating, amplitude',heat_par_exp2(1)/heat_unit,' [ Wm^-3 ]'
            write (*,*) 'Exp heating, flux at z=4 Mm: ', heat_flux
            write (*,*) 'Exp heating, flux at z=0 Mm: ', heat_par_exp2(1)*heat_par_exp2(2)*1e6
          endif
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
!
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
!---------------------------------------------------------------------
        case ('vandoors')
          !heating based on "Van Doorsselaere et al. 2007"
          !pi defined?
!          heat_par_vandoors(1)=1.   !amplitude
          !          heat_par_vandoors(2)=0.1
          !(temp) some amplitude to add a physical heating
          !          heat_par_vandoors(3)=50. Mm,
          !density scale height, later dynamic?
          !  if (headtt) print*,heat_par_vandoors,pi
!
          ! --find a --
          !better in the header of the subroutine?
          a_arr(:,1)=(/0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,2.,2.4,2.8,3.2,3.6,1000. /) !L/pi/H
          a_arr(:,2)=(/0.03426,0.04169,0.04933,0.05717,0.06522,0.07349,0.09068,0.10877,&
              0.12778,0.16853,0.21287,0.26046,0.31068,0.36260,125./) !a
!
           if (notanumber(p%glnrho)) then
   !         print*,p%glnrho
            call fatal_error('p%glnrho','NaN found')
            print*,'-------------------------------'
          endif
!
          call dot2(p%glnrho,Htemp,FAST_SQRT=.true.)
!          Htemp=sqrt(Htemp**2.)
          Hlength=(Htemp+tiny(0D0))**(-1.)/Ltot/pi
!
          j=6
!
          do k=1,nx
            notdone=.true.
            !
            do while (notdone)
              if (Hlength(k) >= a_arr(j,1) .and. Hlength(k)<a_arr(j+1,1)) then
!
                notdone = .false.
              else
                j=j+sign(1.,Hlength(k)-a_arr(j,1))
!
                if (j <= 0) then
                  j=1
                  notdone=.false.
                elseif (j >= 14) then
                  j = 14
                  notdone=.false.
                endif
              endif
            enddo
            slope=(a_arr(j+1,2)-a_arr(j,2))/(a_arr(j+1,1)-a_arr(j,1))
            doors_a(k)=slope*(Hlength(k)-a_arr(j,1))+a_arr(j,2)
          enddo
!
          if (notanumber(doors_a)) then
            print*,Hlength
            call fatal_error('update points','NaN found')
            print*,'-------------------------------'
          endif
!
          where(doors_a < 0)
            doors_a=0D0
          endwhere
          where(doors_a > 0.5)
            doors_a=0.5
          endwhere
!
          doors_heat= heat_par_vandoors(1)* (cos(pi* height/Ltot)**2.+ &
               6.*doors_a*cos(pi*height/Ltot)*cos(3.*pi*height/Ltot))
          where(doors_heat < 0)
            doors_heat=0
          endwhere
!
          do k=0,nx
            if (doors_heat(k) >1. ) then
              print*,doors_heat(k),k
              call fatal_error('doors_heat','toobig!')
            endif
          enddo
!
          heatinput=heatinput + doors_heat!/heat_unit
!print*,max(transpose(heat_par_vandoors(1)*cos(pi* height/Ltot)**2.) )
!
!          print*,max((heat_par_vandoors(1)*cos(pi* height/ L)**2.+6.*heat_par_vandoors(2)*cos(pi*height/L)*cos(3.*pi*height/L))/&
 !             heat_unit)
!---------------------------------------------------------------------
        case default
          if (headtt) call fatal_error('calc_artif_heating', &
              'Please provide correct iheattype')
        endselect
      enddo
      !
      if (headtt) print*,'Total flux for all types:',heat_flux
!
! Add to energy equation
!
      rhs = p%TT1*p%rho1*gamma*p%cp1*heatinput
      if (init_time /=0.) rhs=rhs*cubic_step(real(t),init_time,init_time)
!
      if (ltemperature .and. (.not. ltemperature_nolog)) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
      else if (lentropy .and. (.not. pretend_lnTT)) then
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + rhs/gamma/p%cp1
      else
        call fatal_error('calc_artif_heating', &
            'not for current set of thermodynamic variables')
      endif
!
      if (lfirst.and.ldt) then
        if (ldiagnos.and.idiag_dtnewt/=0) then
!          call max_mn_name(rhs/cdts,idiag_dtnewt,l_dt=.true.)
        endif
        dt1_max=max(dt1_max,rhs/cdts)
      endif
!
    endsubroutine calc_artif_heating
!***********************************************************************
    subroutine calc_heatcond_glnTT_iso(df,p)
!
!  L = Div( Grad(lnT)^2 Grad(T))
!
      use Diagnostics,     only : max_mn_name
      use Sub,             only : dot2,dot,multsv,multmv
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: tmpv
      real, dimension (nx) :: glnT2,glnT_glnr
      real, dimension (nx) :: tmp,rhs,chi
      integer :: i
!
      intent(in) :: p
      intent(out) :: df
!
      call dot2(p%glnTT,glnT2)
      call dot(p%glnTT,p%glnrho,glnT_glnr)
!
      do i=1,3
        tmpv(:,i) = p%glnTT(:,1)*p%hlnTT(:,1,i) + &
                    p%glnTT(:,2)*p%hlnTT(:,2,i) + &
                    p%glnTT(:,3)*p%hlnTT(:,3,i)
      enddo
      call dot(p%glnTT,tmpv,tmp)
!
      chi = glnT2*hcond_grad_iso
!
      rhs = 2*tmp+glnT2*(glnT2+p%del2lnTT+glnT_glnr)
!
      if (ltemperature .and. (.not. ltemperature_nolog)) then
        df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + p%cp1*rhs*gamma*hcond_grad_iso
      else if (lentropy .and. (.not. pretend_lnTT)) then
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss) + rhs*hcond_grad_iso
      else
        call fatal_error('calc_heatcond_glnTT_iso','only for ltemperature')
      endif
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chi*dxyz_2
        if (ldiagnos.and.idiag_dtchi2/=0.) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_glnTT_iso
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
! field-aligned heat conduction that is proportional to rho
! L = Div (b K rho b*Grad(T))
! K = chi [m^2/s] * cV [J/kg/K] = hcond1 [m^4/s^3/K]
! Define chi = K_0/rho
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: rhs, glnrho_b, fdiff
      real, dimension(nx) :: glnTT_H,tmp,hlnTT_Bij,glnTT_b
      real, dimension(nx,3) :: hhh,tmpv
      integer :: i,j,k
!
      call dot(p%glnrho,p%bunit,glnrho_b)
      call dot(p%glnTT,p%bunit,glnTT_b)
!
      do i = 1,3
         hhh(:,i) = 0.
         do j = 1,3
            tmp(:) = 0.
            do k = 1,3
               tmp(:) = tmp(:)-2.*p%bunit(:,k)*p%bij(:,k,j)
            enddo
            hhh(:,i) = hhh(:,i)+p%bunit(:,j)*(p%bij(:,i,j)+p%bunit(:,i)*tmp(:))
         enddo
      enddo
      call dot(hhh,p%glnTT,glnTT_H)
!
      call multmv(p%hlnTT,p%bunit,tmpv)
      call dot(tmpv,p%bunit,hlnTT_Bij)
!
      rhs = (glnTT_H + hlnTT_Bij + (glnrho_b + glnTT_b)*glnTT_b)*chi_aniso*gamma
!
      if (ltemperature) then
         df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
      else
        call fatal_error('calc_heatcond_constchi special','only for ltemperature')
      endif
!
! ANISOTROPIC PART NOT INCLUDED IN THE TIMESTEP CALCULATION
      if (lfirst .and. ldt) then
         advec_cs2 = max(advec_cs2,chi_aniso*maxval(dxyz_2))
         fdiff = gamma * chi_aniso * dxyz_2
         diffus_chi = diffus_chi+fdiff
         if (ldiagnos .and. (idiag_dtchi2 /= 0)) then
            call max_mn_name(fdiff/cdtv,idiag_dtchi2,l_dt=.true.)
         endif
      endif
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine der_upwind(f,uu,k,df,j)
!
!  High order upwind derivative of variable.
!  Useful for advecting.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: df,fac
      integer :: j,k,l
!
      intent(in)  :: f,uu,k,j
      intent(out) :: df
!
      if (nghost /= 3) call fatal_error('der_upwind','only for nghost==3')
      if (lspherical_coords.or.lcylindrical_coords) &
          call fatal_error('der_upwind','NOT IMPLEMENTED for non-cartesian grid')
!
      if (j == 1) then
        if (nxgrid /= 1) then
          fac = 1./6.*dx_1(l1:l2)
          do l=1,nx
            if (uu(l,1) > 0.) then
              df(l) = 11.*f(nghost+l  ,m,n,k)-18.*f(nghost+l-1,m,n,k) &
                     + 9.*f(nghost+l-2,m,n,k)- 2.*f(nghost+l-3,m,n,k)
            else
              df(l) =  2.*f(nghost+l+3,m,n,k)- 9.*f(nghost+l+2,m,n,k) &
                     +18.*f(nghost+l+1,m,n,k)-11.*f(nghost+l  ,m,n,k)
            endif
          enddo
          df = fac*df
        else
          df=0.
          if (ip<=5) print*, 'der_upwind: Degenerate case in x-direction'
        endif
      elseif (j == 2) then
        if (nygrid /= 1) then
          fac = 1./6.*dy_1(m)
          do l=1,nx
            if (uu(l,2) > 0.) then
              df(l) = 11.*f(nghost+l,m  ,n,k)-18.*f(nghost+l,m-1,n,k) &
                     + 9.*f(nghost+l,m-2,n,k)- 2.*f(nghost+l,m-3,n,k)
            else
              df(l) =  2.*f(nghost+l,m+3,n,k)- 9.*f(nghost+l,m+2,n,k) &
                     +18.*f(nghost+l,m+1,n,k)-11.*f(nghost+l,m  ,n,k)
            endif
          enddo
          df = fac*df
        else
          df=0.
          if (ip<=5) print*, 'der_upwind: Degenerate case in y-direction'
        endif
      elseif (j == 3) then
        if (nzgrid /= 1) then
          fac = 1./6.*dz_1(n)
          do l=1,nx
            if (uu(l,3) > 0.) then
              df(l) = 11.*f(nghost+l,m,n  ,k)-18.*f(nghost+l,m,n-1,k) &
                     + 9.*f(nghost+l,m,n-2,k)- 2.*f(nghost+l,m,n-3,k)
            else
              df(l) =  2.*f(nghost+l,m,n+3,k)- 9.*f(nghost+l,m,n+2,k) &
                     +18.*f(nghost+l,m,n+1,k)-11.*f(nghost+l,m,n  ,k)
            endif
          enddo
          df = fac*df
        else
          df=0.
          if (ip<=5) print*, 'der_upwind: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_upwind
!***********************************************************************
    subroutine calc_bullets_heating(df,p)
!
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + bullets_h0* &
          exp( - ((x(l1:l2)-bullets_x0)/bullets_dx)**2.)* &
          exp( - ((t       -bullets_t0)/bullets_dt)**2.)* &
          p%cp1*gamma*exp(-p%lnrho -p%lnTT)
!
    endsubroutine calc_bullets_heating
!***********************************************************************
    subroutine filter_farray(f,df)
!
!  Reduce noise of farray using mesh independend hyper diffusion.
!  Filter strength is equal for all variables up to now and has
!  to be smaller than 1/64 for numerical stability.
!
!  13-jan-12/bing: coded
!
      use Sub, only: del6
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: del6_fj
      integer :: j
!
      if (dt/=0.) then
        do j=1,mvar
          if (filter_strength(j)/=0.) then
            call del6(f,j,del6_fj,IGNOREDX=.true.)
            df(l1:l2,m,n,j) = df(l1:l2,m,n,j) + &
                filter_strength(j)*del6_fj/dt_beta_ts(itsub)
          endif
        enddo
      endif
!
    endsubroutine filter_farray
!***********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
