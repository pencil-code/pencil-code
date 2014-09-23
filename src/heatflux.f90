! $Id$
!
! This module computes the time depentend heat flux equation.
! The non-Fourier form is written
!    tau * dq/dt + q  = F(T,rho ...)
! The analytical solution for constant F is that
! q converges exponential to F with a time scale of tau.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lheatflux = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED qq(3); q2; divq
!
!***************************************************************
!
module Heatflux
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  real :: dummy=0.
!
  namelist /heatflux_init_pars/ &
      dummy
!
  character (len=labellen) :: iheatflux='nothing'
  logical :: lreset_heatflux
  real :: tau_heatflux=0.,saturation_flux=0.,tau1_eighthm=0.
!
  namelist /heatflux_run_pars/ &
      lreset_heatflux,iheatflux,tau_heatflux,saturation_flux,  &
      tau1_eighthm
!
! variables for print.in
!
  integer :: idiag_qmax=0     ! DIAG_DOC: max of heat flux vector
  integer :: idiag_qrms=0     ! DIAG_DOC: rms of heat flux vector
!
!  variables for video slices:
!
  real, target, dimension (nx,ny) :: hflux_xy,hflux_xy2
  real, target, dimension (nx,ny) :: hflux_xy3,hflux_xy4
  real, target, dimension (nx,nz) :: hflux_xz
  real, target, dimension (ny,nz) :: hflux_yz
!
  real, dimension(:), pointer :: B_ext
  integer :: ierr
  real :: tau_inv_spitzer=0.,Kspitzer_para=0.
  real :: saturation_fluxuration=0.,nu_ee,ln_unit_TT
  real :: Ksaturation_SI = 7e7,Ksaturation=0.
  real :: e_m
!
  real, target, dimension (nx,ny) :: divq_xy,divq_xy2,divq_xy3,divq_xy4
  real, target, dimension (nx,nz) :: divq_xz
  real, target, dimension (ny,nz) :: divq_yz
!
  include 'heatflux.h'
!
contains
!***********************************************************************
  subroutine register_heatflux()
!
!  Set up indices for variables in heatflux modules.
!
!  24-apr-13/bing: coded
!
    use FArrayManager, only: farray_register_pde
!
    call farray_register_pde('heatflux',iqq,vector=3)
    iqx=iqq; iqy=iqq+1; iqz=iqq+2
!
    if (lroot) call svn_id( &
        "$Id$")
!
  endsubroutine register_heatflux
!***********************************************************************
  subroutine initialize_heatflux(f,lstarting)
!
!  Called by start.f90 together with lstarting=.true.   and then
!  called by run.f90   together with lstarting=.false.  after reading
!  parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
    use SharedVariables, only: get_shared_variable
!
    real, dimension (mx,my,mz,mfarray) :: f
    logical :: lstarting
    real :: eps0,unit_ampere,e_charge
!
    call keep_compiler_quiet(lstarting)
!
!  Get the external magnetic field if exists.
    if (lmagnetic) then
      call get_shared_variable('B_ext', B_ext, ierr)
      if (ierr /= 0) call fatal_error('calc_hcond_timestep',  &
          'unable to get shared variable B_ext')
    endif
!
!  Set up some important constants
!
    Ksaturation = Ksaturation_SI /unit_velocity**3. * unit_temperature**1.5
!
    eps0 =1./mu0/c_light**2.
!
    unit_ampere = unit_velocity*sqrt(mu0/(4.*pi*1e-7)*unit_density)*unit_length
    e_charge= 1.602176e-19/unit_time/unit_ampere
!
    e_m = e_charge / m_e
!
!  compute the constant for the collision frequency
!  nu_ee = constant
    nu_ee = 4.D0/3. *sqrt(pi) * 20.D0/ sqrt(m_e) * &
        ((e_charge**2./(4.D0*pi*eps0))**2.D0)/(k_B**1.5) * 0.872 /m_p
!
    if (lreset_heatflux) f(:,:,:,iqx:iqz)=0.
!
  endsubroutine initialize_heatflux
!***********************************************************************
  subroutine finalize_heatflux(f,lstarting)
!
!  Called by start.f90 together with lstarting=.true.   and then
!  called by run.f90   together with lstarting=.false.  before exiting.
!
!  14-aug-2011/Bourdin.KIS: coded
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
    logical, intent(in) :: lstarting
!
    call keep_compiler_quiet(f)
    call keep_compiler_quiet(lstarting)
!
  endsubroutine finalize_heatflux
!***********************************************************************
  subroutine init_heatflux(f)
!
!  initialise heatflux condition; called from start.f90
!  06-oct-2003/tony: coded
!
    real, dimension (mx,my,mz,mfarray) :: f
!
    intent(inout) :: f
!
    f(:,:,:,iqx:iqz)=0.
!
  endsubroutine init_heatflux
!***********************************************************************
  subroutine pencil_criteria_heatflux()
!
!  All pencils that this heatflux module depends on are specified here.
!
!  18-07-06/tony: coded
!
    if (iheatflux == 'eighth') then
      lpenc_requested(i_qq)=.true.
      lpenc_requested(i_divq)=.true.
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_b2)=.true.
      lpenc_requested(i_cp1)=.true.
      lpenc_requested(i_lnTT)=.true.
      lpenc_requested(i_glnTT)=.true.
      lpenc_requested(i_lnrho)=.true.
      lpenc_requested(i_uu)=.true.
      lpenc_requested(i_divu)=.true.
      lpenc_requested(i_uij)=.true.
    endif
!
    if (iheatflux == 'spitzer') then
      lpenc_requested(i_qq)=.true.
      lpenc_requested(i_divq)=.true.
      lpenc_requested(i_cp1)=.true.
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_bunit)=.true.
      lpenc_requested(i_b2)=.true.
      lpenc_requested(i_uglnrho)=.true.
      lpenc_requested(i_divu)=.true.
      lpenc_requested(i_lnTT)=.true.
      lpenc_requested(i_glnTT)=.true.
      lpenc_requested(i_lnrho)=.true.
      lpenc_requested(i_glnrho)=.true.
    endif
!
    if (idiag_qmax/=0 .or. idiag_qrms/=0) lpenc_diagnos(i_q2)=.true.
!
  endsubroutine pencil_criteria_heatflux
!***********************************************************************
  subroutine pencil_interdep_heatflux(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
    logical, dimension(npencils), intent(inout) :: lpencil_in
!
    if (lpencil_in(i_q2)) lpencil_in(i_qq)=.true.
!
  endsubroutine pencil_interdep_heatflux
!***********************************************************************
  subroutine calc_pencils_heatflux(f,p)
!
!  Calculate Heatflux pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
    use Sub, only: dot2_mn, div
!
    real, dimension (mx,my,mz,mfarray) :: f
    type (pencil_case) :: p
!
    intent(in) :: f
    intent(inout) :: p
!
    if (lpencil(i_qq)) p%qq=f(l1:l2,m,n,iqx:iqz)
    if (lpencil(i_q2)) call dot2_mn(p%qq,p%q2)
    if (lpencil(i_divq)) call div(f,iqq,p%divq)
!
  endsubroutine calc_pencils_heatflux
!***********************************************************************
  subroutine dheatflux_dt(f,df,p)
!
!  26-apr-13/bing: coded
!
    use Diagnostics, only: max_mn_name,sum_mn_name
    use Sub, only: identify_bcs, del6
!
    real, dimension (mx,my,mz,mfarray) :: f
    real, dimension (mx,my,mz,mvar) :: df
    type (pencil_case) :: p
    real, dimension (nx) :: hc,hyper3_coeff
    integer :: i
!
    intent(in) :: f,p
    intent(inout) :: df
!
    call keep_compiler_quiet(f)
!
    if (headtt.or.ldebug) print*,'dheatflux_dt: SOLVE dheatflux_dt'
!
    if (headtt) then
      call identify_bcs('heatfluxx',iqx)
      call identify_bcs('heatfluxy',iqy)
      call identify_bcs('heatfluxz',iqz)
    endif
!
    select case (iheatflux)
!
    case ('spitzer')
      call non_fourier_spitzer(df,p)
    case ('eighth')
      call eighth_moment_approx(df,p)
    endselect
!
    if (idiag_qmax/=0) call max_mn_name(p%q2,idiag_qmax,lsqrt=.true.)
    if (idiag_qrms/=0) call sum_mn_name(p%q2,idiag_qrms,lsqrt=.true.)
!
  endsubroutine dheatflux_dt
!***********************************************************************
  subroutine read_heatflux_init_pars(unit,iostat)
!
    integer, intent(in) :: unit
    integer, intent(inout), optional :: iostat
!
    if (present(iostat)) then
      read(unit,NML=heatflux_init_pars,ERR=99, IOSTAT=iostat)
    else
      read(unit,NML=heatflux_init_pars,ERR=99)
    endif
!
99  return
!
  endsubroutine read_heatflux_init_pars
!***********************************************************************
  subroutine write_heatflux_init_pars(unit)
!
    integer, intent(in) :: unit
!
    write(unit,NML=heatflux_init_pars)
!
  endsubroutine write_heatflux_init_pars
!***********************************************************************
  subroutine read_heatflux_run_pars(unit,iostat)
!
    integer, intent(in) :: unit
    integer, intent(inout), optional :: iostat
!
    if (present(iostat)) then
      read(unit,NML=heatflux_run_pars,ERR=99, IOSTAT=iostat)
    else
      read(unit,NML=heatflux_run_pars,ERR=99)
    endif
!
99  return
!
  endsubroutine read_heatflux_run_pars
!***********************************************************************
  subroutine write_heatflux_run_pars(unit)
!
    integer, intent(in) :: unit
!
    write(unit,NML=heatflux_run_pars)
!
  endsubroutine write_heatflux_run_pars
!***********************************************************************
  subroutine rprint_heatflux(lreset,lwrite)
!
!  Reads and registers print parameters relevant to heatflux.
!
!  06-oct-03/tony: coded
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
    if (lreset) then
      idiag_qmax=0
      idiag_qrms=0
    endif
!
!  iname runs through all possible names that may be listed in print.in
!
    do iname=1,nname
      call parse_name(iname,cname(iname),cform(iname),'qmax',idiag_qmax)
      call parse_name(iname,cname(iname),cform(iname),'qrms',idiag_qrms)
    enddo
!
!  write column where which variable is stored
!
    if (lwr) then
      write(3,*) 'i_qmax=',idiag_qmax
      write(3,*) 'i_qrms=',idiag_qrms
    endif
!
  endsubroutine rprint_heatflux
!***********************************************************************
  subroutine get_slices_heatflux(f,slices)
!
!  Write slices for animation of Heatflux variables.
!
!  26-jun-06/tony: dummy
!
    real, dimension (mx,my,mz,mfarray) :: f
    type (slice_data) :: slices
!
!  Loop over slices
!
    select case (trim(slices%name))
!
    case ('hflux')
      if (slices%index>=3) then
        slices%ready=.false.
      else
        slices%index=slices%index+1
        slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,iqx-1+slices%index)
        slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,iqx-1+slices%index)
        slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,iqx-1+slices%index)
        slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,iqx-1+slices%index)
        if (lwrite_slice_xy3) &
            slices%xy3=f(l1:l2,m1:m2,iz3_loc,iqx-1+slices%index)
        if (lwrite_slice_xy4) &
              slices%xy4=f(l1:l2,m1:m2,iz4_loc,iqx-1+slices%index)
        if (slices%index<=3) slices%ready=.true.
      endif
!
      case ('divq')
        slices%yz => divq_yz
        slices%xz => divq_xz
        slices%xy => divq_xy
        slices%xy2=> divq_xy2
        if (lwrite_slice_xy3) slices%xy3=> divq_xy3
        if (lwrite_slice_xy4) slices%xy4=> divq_xy4
        slices%ready=.true.
    endselect
!
  endsubroutine get_slices_heatflux
!***********************************************************************
  subroutine non_fourier_spitzer(df,p)
!
!  26-jun-06/tony: dummy
!
    use EquationOfState
    use Sub
!
    real, dimension (mx,my,mz,mvar) :: df
    type (pencil_case) :: p
    real, dimension(nx) :: b2_1,qsat,qabs
    real, dimension(nx) :: rhs,cosgT_b
    real, dimension(nx,3) :: K1,unit_glnTT
    real, dimension(nx,3) :: spitzer_vec
    real, dimension(nx) :: tmp
    integer :: i
!
    b2_1=1./(p%b2+tini)
!
    call multsv(Kspitzer_para*exp(-p%lnrho+3.5*p%lnTT),p%glnTT,K1)
!
    call dot(K1,p%bb,tmp)
    call multsv(-b2_1*tmp,p%bb,spitzer_vec)
!
! Reduce the heat conduction at places of low density or very
! high temperatures
!
    if (saturation_flux/=0.) then
      call dot2(spitzer_vec,qabs,FAST_SQRT=.true.)
      qsat = saturation_flux* exp(1.5*p%lnTT)
!
      qsat = 1./(1./qsat +1./qabs)
      where (qabs > sqrt(tini))
        spitzer_vec(:,1) = spitzer_vec(:,1)*qsat/qabs
        spitzer_vec(:,2) = spitzer_vec(:,2)*qsat/qabs
        spitzer_vec(:,3) = spitzer_vec(:,3)*qsat/qabs
      endwhere
    endif
!
    do i=1,3
      df(l1:l2,m,n,iqq+i-1) = df(l1:l2,m,n,iqq+i-1) + &
          tau_inv_spitzer*(-p%qq(:,i)+spitzer_vec(:,i))  +  &
          p%qq(:,i)*(p%uglnrho + p%divu)
      write (*,*) maxval(abs(df(l1:l2,m,n,iqq+i-1)))
    enddo
!
    call dot(p%qq,p%glnrho,tmp)
!
    rhs = gamma*p%cp1*(p%divq + tmp)*exp(-p%lnTT)
!
    df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - rhs
!
    if (lfirst.and.ldt) then
      call unit_vector(p%glnTT,unit_glnTT)
      call dot(unit_glnTT,p%bunit,cosgT_b)
      rhs = sqrt(Kspitzer_para*exp(2.5*p%lnTT-p%lnrho)* &
          gamma*p%cp1*tau_inv_spitzer*abs(cosgT_b))
      advec_uu = max(advec_uu,rhs/dxmax_pencil)
!
!          if (idiag_dtspitzer/=0) &
!              call max_mn_name(advec_uu/cdt,idiag_dtspitzer,l_dt=.true.)
 !         !
      dt1_max=max(dt1_max,tau_inv_spitzer/cdts)
    endif
!
  endsubroutine non_fourier_spitzer
!***********************************************************************
  subroutine eighth_moment_approx(df,p)
!
!  This heatconduction refers to the eighth moment approximation
!
!  26-jun-06/tony: dummy
!
    use EquationOfState, only: gamma
    use Sub
!
    real, dimension (mx,my,mz,mfarray) :: f
    real, dimension (mx,my,mz,mvar) :: df
    type (pencil_case) :: p
    real, dimension(nx) :: b2_1,rhs
    real, dimension(nx,3) :: K1
    real, dimension(nx) :: dt_1_8th,nu_coll
    real :: coeff,nu_coll_max=1e3
    integer :: i
!
    real, dimension(nx,3,3) ::  qij
    real, dimension(nx,3) :: qgradu,ugradq,qmu,qdivu
    real, dimension(nx,3) :: qxB,tmp_hf
!
    call gij(f,iqx,qij,1)
    call u_dot_grad(f,iqx,qij,p%uu,ugradq)
!
    call u_dot_grad(f,iux,p%uij,p%qq,qgradu)
    call multsv(p%divu,p%qq,qdivu)
!
    call multmv(p%uij,p%qq,qmu)
!
    call cross(p%qq,p%bb,qxB)
!
    tmp_hf =-(ugradq +7./5.*qgradu+7./5.*qdivu+2./5.*qmu) ! - e_m*qxB)
!
    coeff = 0.872*5./2.*(k_B/m_e)*(k_B/m_p)
    call multsv(coeff*exp(2.0*p%lnTT+p%lnrho),p%glnTT,K1)
!
    tmp_hf = tmp_hf - K1
!
    nu_coll = 16./35.*nu_ee*exp(p%lnrho-1.5*p%lnTT)
    where (nu_coll > nu_coll_max)
      nu_coll = nu_coll_max
    endwhere
!
    do i=1,3
      tmp_hf(:,i) = tmp_hf(:,i) - nu_coll*p%qq(:,i)
      !
      df(l1:l2,m,n,iqq+i-1) = df(l1:l2,m,n,iqq+i-1) + tmp_hf(:,i)
    enddo
!
! Add divergence of the flux to the energy equation
!
    rhs = exp(-p%lnrho-p%lnTT)*p%divq
!
   df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - gamma*p%cp1*rhs

   if (lvideo) then
     divq_yz(m-m1+1,n-n1+1)=rhs(ix_loc-l1+1)
     if (m == iy_loc)  divq_xz(:,n-n1+1)= rhs
     if (n == iz_loc)  divq_xy(:,m-m1+1)= rhs
     if (n == iz2_loc) divq_xy2(:,m-m1+1)= rhs
     if (n == iz3_loc) divq_xy3(:,m-m1+1)= rhs
     if (n == iz4_loc) divq_xy4(:,m-m1+1)= rhs
   endif
!
    if (lfirst.and.ldt) then
!      b2_1=1./(p%b2+tini)
      dt_1_8th = nu_coll !+ e_m /sqrt(b2_1)

      dt1_max=max(dt1_max,dt_1_8th/cdts)
      advec_uu = max(advec_uu,advec_uu*7./5.)
!     advec_uu = max(advec_uu,sqrt(coeff*exp(p%lnTT)*gamma*p%cp1)/dxmax_pencil)
    endif
!
  endsubroutine eighth_moment_approx
!***********************************************************************
endmodule Heatflux
