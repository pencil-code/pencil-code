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
  character (len=labellen) :: iheatflux='nothing'
  logical :: lreset_heatflux=.false., lnfs2=.false., ltau_spitzer_va=.false.
  real :: saturation_flux=0.,tau1_eighthm=0.,tau_inv_spitzer=0.
  real :: Ksaturation_SI = 7e7, Ksaturation=0.,Kc=0.
  real :: Kspitzer_para=0.,hyper3_coeff=0., va2max_tau_boris=0.
!
  namelist /heatflux_run_pars/ &
      lreset_heatflux,iheatflux,saturation_flux,  &
      tau1_eighthm,Kspitzer_para,tau_inv_spitzer, &
      hyper3_coeff, lnfs2, ltau_spitzer_va, &
      Kc, va2max_tau_boris
  real, dimension(:), pointer :: B_ext
  real :: nu_ee, e_m
!
!  variables for video slices:
!
  real, target, dimension (:,:), allocatable :: divq_xy,divq_xy2,divq_xy3,divq_xy4
  real, target, dimension (:,:), allocatable :: divq_xz,divq_xz2,divq_yz
!
! Diagnostic variables (need to be consistent with reset list below)
!
  integer :: idiag_dtspitzer=0  ! DIAG_DOC: Spitzer heat conduction time step
  integer :: idiag_dtq=0        ! DIAG_DOC: heatflux time step
  integer :: idiag_dtq2=0       ! DIAG_DOC: heatflux time step due to tau
  integer :: idiag_qmax=0       ! DIAG_DOC: $\max(|\qv|)$
  integer :: idiag_tauqmax=0    ! DIAG_DOC: $\max(|\tau_{\rm Spitzer}|)$
  integer :: idiag_qxmin=0      ! DIAG_DOC: $\min(|q_x|)$
  integer :: idiag_qymin=0      ! DIAG_DOC: $\min(|q_y|)$
  integer :: idiag_qzmin=0      ! DIAG_DOC: $\min(|q_z|)$
  integer :: idiag_qxmax=0      ! DIAG_DOC: $\max(|q_x|)$
  integer :: idiag_qymax=0      ! DIAG_DOC: $\max(|q_y|)$
  integer :: idiag_qzmax=0      ! DIAG_DOC: $\max(|q_z|)$
  integer :: idiag_qrms=0       ! DIAG_DOC: $\sqrt(|\qv|^2)$
  integer :: idiag_qsatmin=0    ! DIAG_DOC: minimum of qsat/qabs
  integer :: idiag_qsatrms=0    ! DIAG_DOC: rms of qsat/abs
  integer :: ivid_divq
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
    call farray_register_pde('qq',iqq,vector=3)
    iqx=iqq; iqy=iqq+1; iqz=iqq+2
!
    if (lroot) call svn_id( &
        "$Id$")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',qq $'
          if (nvar == mvar) write(4,*) ',qq'
        else
          write(4,*) ',qq $'
        endif
        write(15,*) 'qq = fltarr(mx,my,mz,3)*one'
      endif
!
  endsubroutine register_heatflux
!***********************************************************************
  subroutine initialize_heatflux(f)
!
!  Called after reading parameters, but before the time loop.
!
!  07-sept-17/bingert: updated
!
    !use Slices_methods, only: alloc_slice_buffers
    use SharedVariables, only: get_shared_variable
!
    real, dimension (mx,my,mz,mfarray) :: f
    real :: eps0,unit_ampere,e_charge
!
!  Get the external magnetic field if exists.
    if (lmagnetic) &
      call get_shared_variable('B_ext', B_ext, caller='calc_hcond_timestep')
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
    if (ivid_divq/=0) then
      !call alloc_slice_buffers(divq_xy,divq_xz,divq_yz,divq_xy2,divq_xy3,divq_xy4,divq_xz2)
      if (lwrite_slice_xy .and..not.allocated(divq_xy) ) allocate(divq_xy (nx,ny))
      if (lwrite_slice_xz .and..not.allocated(divq_xz) ) allocate(divq_xz (nx,nz))
      if (lwrite_slice_yz .and..not.allocated(divq_yz) ) allocate(divq_yz (ny,nz))
      if (lwrite_slice_xy2.and..not.allocated(divq_xy2)) allocate(divq_xy2(nx,ny))
      if (lwrite_slice_xy3.and..not.allocated(divq_xy3)) allocate(divq_xy3(nx,ny))
      if (lwrite_slice_xy4.and..not.allocated(divq_xy4)) allocate(divq_xy4(nx,ny))
      if (lwrite_slice_xz2.and..not.allocated(divq_xz2)) allocate(divq_xz2(nx,nz))
    endif

  endsubroutine initialize_heatflux
!***********************************************************************
  subroutine finalize_heatflux(f)
!
!  Called right before exiting.
!
!  14-aug-2011/Bourdin.KIS: coded
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    call keep_compiler_quiet(f)
!
  endsubroutine finalize_heatflux
!***********************************************************************
  subroutine init_heatflux(f)
!
!  initialise heatflux condition; called from start.f90
!  07-sept-17/bingert: updated
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
!  07-sept-17/bingert: updated
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
      if (ltau_spitzer_va) lpenc_requested(i_va2)=.true.
    endif
!
    if (iheatflux == 'noadvection-spitzer') then
      lpenc_requested(i_qq)=.true.
      lpenc_requested(i_divq)=.true.
      lpenc_requested(i_cp1)=.true.
      lpenc_requested(i_cv1)=.true.
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_bunit)=.true.
      lpenc_requested(i_b2)=.true.
      lpenc_requested(i_uglnrho)=.true.
      lpenc_requested(i_lnTT)=.true.
      lpenc_requested(i_glnTT)=.true.
      lpenc_requested(i_lnrho)=.true.
      lpenc_requested(i_glnrho)=.true.
      lpenc_requested(i_rho)=.true.
      lpenc_requested(i_rho1)=.true.
      lpenc_requested(i_TT)=.true.
      lpenc_requested(i_TT1)=.true.
      lpenc_requested(i_gamma)=.true.
    endif
!
    if (idiag_qmax/=0 .or. idiag_qrms/=0) then
       lpenc_diagnos(i_q2)=.true.
       lpenc_requested(i_q2)=.true.
       lpenc_requested(i_lnrho)=.true.
    endif
!
    if (idiag_qxmax/=0 .or. idiag_qxmin/=0) lpenc_requested(i_qq)=.true.
    if (idiag_qymax/=0 .or. idiag_qymin/=0) lpenc_requested(i_qq)=.true.
    if (idiag_qzmax/=0 .or. idiag_qzmin/=0) lpenc_requested(i_qq)=.true.
!
  endsubroutine pencil_criteria_heatflux
!***********************************************************************
  subroutine pencil_interdep_heatflux(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  07-sept-17/bingert: updated
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
!  07-sept-17/bingert: updated

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
    real, dimension (nx) :: hc

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
    case ('noadvection-spitzer')
      call noadvection_non_fourier_spitzer(df,p)
    case ('eighth')
      call eighth_moment_approx(df,p)
    endselect
!
    if (hyper3_coeff /= 0.) then
       call del6(f,iqx,hc,IGNOREDX=.true.)
       df(l1:l2,m,n,iqx) = df(l1:l2,m,n,iqx) + hyper3_coeff*hc
       call del6(f,iqy,hc,IGNOREDX=.true.)
       df(l1:l2,m,n,iqy) = df(l1:l2,m,n,iqy) + hyper3_coeff*hc
       call del6(f,iqz,hc,IGNOREDX=.true.)
       df(l1:l2,m,n,iqz) = df(l1:l2,m,n,iqz) + hyper3_coeff*hc
    endif
!
    if (idiag_qmax/=0) call max_mn_name(p%q2,idiag_qmax,lsqrt=.true.)
    if (idiag_qrms/=0) call sum_mn_name(p%q2,idiag_qrms,lsqrt=.true.)
    if (idiag_qxmin/=0) call max_mn_name(-p%qq(:,1),idiag_qxmin,lneg=.true.)
    if (idiag_qymin/=0) call max_mn_name(-p%qq(:,2),idiag_qymin,lneg=.true.)
    if (idiag_qzmin/=0) call max_mn_name(-p%qq(:,3),idiag_qzmin,lneg=.true.)
    if (idiag_qxmax/=0) call max_mn_name(p%qq(:,1),idiag_qxmax)
    if (idiag_qymax/=0) call max_mn_name(p%qq(:,2),idiag_qymax)
    if (idiag_qzmax/=0) call max_mn_name(p%qq(:,3),idiag_qzmax)
!
  endsubroutine dheatflux_dt
!***********************************************************************
    subroutine read_heatflux_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=heatflux_run_pars, IOSTAT=iostat)
!
    endsubroutine read_heatflux_run_pars
!***********************************************************************
    subroutine write_heatflux_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=heatflux_run_pars)
!
    endsubroutine write_heatflux_run_pars
!***********************************************************************
  subroutine rprint_heatflux(lreset,lwrite)
!
!  Reads and registers print parameters relevant to heatflux.
!
!  07-sept-17/bingert: updated
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
    if (lreset) then
      idiag_qmax=0
      idiag_qrms=0
      idiag_dtspitzer=0
      idiag_dtq=0
      idiag_dtq2=0
      idiag_tauqmax=0
      idiag_qsatmin=0
      idiag_qsatrms=0
      ivid_divq=0
    endif
!
!  iname runs through all possible names that may be listed in print.in
!
    do iname=1,nname
      call parse_name(iname,cname(iname),cform(iname),'dtspitzer',idiag_dtspitzer)
      call parse_name(iname,cname(iname),cform(iname),'dtq',idiag_dtq)
      call parse_name(iname,cname(iname),cform(iname),'dtq2',idiag_dtq2)
      call parse_name(iname,cname(iname),cform(iname),'tauqmax',idiag_tauqmax)
      call parse_name(iname,cname(iname),cform(iname),'qmax',idiag_qmax)
      call parse_name(iname,cname(iname),cform(iname),'qrms',idiag_qrms)
      call parse_name(iname,cname(iname),cform(iname),'qsatmin',idiag_qsatmin)
      call parse_name(iname,cname(iname),cform(iname),'qsatrms',idiag_qsatrms)
    enddo
!
!  check for those quantities for which we want video slices
!
    if (lwrite_slices) then 
      where(cnamev=='hflux') cformv='DEFINED'
    endif
    do iname=1,nnamev
      call parse_name(iname,cnamev(iname),cformv(iname),'divq',ivid_divq)
    enddo

    if (lwr) then
      call farray_index_append('iqq',iqq)
      call farray_index_append('iqx',iqx)
      call farray_index_append('iqy',iqy)
      call farray_index_append('iqz',iqz)
      call farray_index_append('i_dtspitzer',idiag_dtspitzer)
      call farray_index_append('i_dtq',idiag_dtq)
      call farray_index_append('i_dtq2',idiag_dtq2)
      call farray_index_append('i_qsatmin',idiag_qsatmin)
      call farray_index_append('i_qsatrms',idiag_qsatrms)
    endif
!
  endsubroutine rprint_heatflux
!***********************************************************************
  subroutine get_slices_heatflux(f,slices)
!
!  Write slices for animation of Heatflux variables.
!
!  07-sept-17/bingert: updated
!
    use Slices_methods, only: assign_slices_scal

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
        if (lwrite_slice_yz) &
          slices%yz=f(ix_loc,m1:m2 ,n1:n2,iqx-1+slices%index) * &
                    exp(-f(ix_loc,m1:m2 ,n1:n2  ,ilnrho))
        if (lwrite_slice_xz) &
          slices%xz=f(l1:l2 ,iy_loc,n1:n2,iqx-1+slices%index) * &
                    exp(-f(l1:l2,iy_loc,n1:n2  ,ilnrho))
        if (lwrite_slice_xy) &
          slices%xy=f(l1:l2 ,m1:m2 ,iz_loc,iqx-1+slices%index) * &
                    exp(-f(l1:l2,m1:m2 ,iz_loc ,ilnrho))
        if (lwrite_slice_xy2) &
          slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,iqx-1+slices%index) * &
                     exp(-f(l1:l2,m1:m2 ,iz2_loc,ilnrho))
        if (lwrite_slice_xy3) &
          slices%xy3=f(l1:l2,m1:m2,iz3_loc,iqx-1+slices%index) * &
                       exp(-f(l1:l2,m1:m2,iz3_loc,ilnrho))
        if (lwrite_slice_xy4) &
          slices%xy4=f(l1:l2,m1:m2,iz4_loc,iqx-1+slices%index) * &
                       exp(-f(l1:l2,m1:m2,iz4_loc,ilnrho))
        if (lwrite_slice_xz2) &
          slices%xz2=f(l1:l2 ,iy_loc,n1:n2  ,iqx-1+slices%index) * &
                     exp(-f(l1:l2,iy2_loc,n1:n2,ilnrho))
        if (slices%index<=3) slices%ready=.true.
      endif
!
    case ('divq')
      call assign_slices_scal(slices,divq_xy,divq_xz,divq_yz,divq_xy2,divq_xy3,divq_xy4,divq_xz2)

    endselect
!
  endsubroutine get_slices_heatflux
!***********************************************************************
  subroutine non_fourier_spitzer(df,p)
!
!  07-sept-17/bingert: updated
!
    use Slices_methods, only: store_slices
    use Diagnostics, only: max_mn_name,sum_mn_name, max_name
    use EquationOfState
    use Sub
!
    real, dimension (mx,my,mz,mvar) :: df
    type (pencil_case) :: p
    real, dimension(nx) :: b2_1,qsat,qabs, dt1_va
    real, dimension(nx) :: rhs,cosgT_b, Kspitzer, K_clight
    real, dimension(nx,3) :: K1,unit_glnTT
    real, dimension(nx,3) :: spitzer_vec
    real, dimension(nx) :: tmp, diffspitz, tau_inv_va
    real, dimension(nx) :: c_spitzer, c_spitzer0
    real :: uplim
    integer :: i
!
! Compute Spizter coefficiant K_0 * T^(5/2)
!
    b2_1=1./(p%b2+tini)
!
    if (lnfs2) then
!
!    for pp=qq/rho there you must divide by rho
!
     Kspitzer=Kspitzer_para*exp(3.5*p%lnTT-p%lnrho)
!
!    Limit the heat flux by the speed of light
!    Kc should be on the order of unity or smaler
!
      if (Kc /= 0.) then
        K_clight = Kc*c_light*dxmin_pencil/(p%cp1*gamma)
!
        where (Kspitzer > K_clight)
          Kspitzer=K_clight
        endwhere
      endif
    else
!
!    For pp=qq*rho there you must a factor of rho
!
     Kspitzer=Kspitzer_para*exp(p%lnrho+3.5*p%lnTT)
!
!    Limit the heat flux by the speed of light
!    Kc should be on the order of unity or smaler
!
      if (Kc /= 0.) then
        K_clight = Kc*c_light*dxmin_pencil*exp(2*p%lnrho)/(p%cp1*gamma)
!
        where (Kspitzer > K_clight)
          Kspitzer=K_clight
        endwhere
      endif
    endif

    call multsv(Kspitzer,p%glnTT,K1)
    call dot(K1,p%bb,tmp)
    call multsv(b2_1*tmp,p%bb,spitzer_vec)
!
!   Reduce the heat conduction at places of low density or very
!   high temperatures
!
    if (saturation_flux/=0.) then
      call dot2(spitzer_vec,qabs,FAST_SQRT=.true.)
!
      if (lnfs2) then
!       for pp=qq/rho, there is no rho factor
        qsat = saturation_flux* exp(1.5*p%lnTT) * Ksaturation
      else
!       for pp=qq*rho, there is an additional rho factor
        qsat = saturation_flux* exp(2.*p%lnrho+1.5*p%lnTT) * Ksaturation
      endif
!
      where (qabs > sqrt(tini))
        qsat = 1./(1./qsat +1./qabs)
        spitzer_vec(:,1) = spitzer_vec(:,1)*qsat/qabs
        spitzer_vec(:,2) = spitzer_vec(:,2)*qsat/qabs
        spitzer_vec(:,3) = spitzer_vec(:,3)*qsat/qabs
      endwhere
      if (ldiagnos) then
!
!   pc_auto-test may digest at maximum 2 digits in the exponent
!
        tmp = qsat/(qabs+sqrt(tini))
        where (tmp > 1d50) tmp = 1d50
        if (idiag_qsatmin/=0) call max_mn_name(-tmp,idiag_qsatmin,lneg=.true.)
        if (idiag_qsatrms/=0) call sum_mn_name(tmp,idiag_qsatrms)
      endif
    endif
!
!
!
    if (ltau_spitzer_va) then
!
!   adjust tau_spitzer to fullfill: propagation speed sqrt(diffspitz/tau)=sqrt(2)*va
!   use tau_inv_spitzer as lower limit for tau_inv -> upper limit for tau
!   use 2 times alfven and advective as upper limit for tau_inv -> lower limit for tau
!
      call unit_vector(p%glnTT,unit_glnTT)
      call dot(unit_glnTT,p%bunit,cosgT_b)
!
!     diffspitz is actually chi * gamma with chi being the usual
!     heat diffusivity, however diffspitz is the one entering the timestep
!     see EQ 24 in the manual.
!
      diffspitz = Kspitzer_para*exp(2.5*p%lnTT-p%lnrho)* &
                 gamma*p%cp1*abs(cosgT_b)
!
      if (va2max_tau_boris /= 0) then
!
!   va2max_tau_boris musst be set to the same value as va2max_boris
!   in magnetic_run_pars
!
        tmp = (1+(p%va2/va2max_tau_boris)**2.)**(-1.0/2.0)
        tau_inv_va = 2.*p%va2*tmp/(diffspitz+sqrt(tini))
        dt1_va=sqrt(p%va2*tmp*dxyz_2)
      else
        tau_inv_va = 2.*p%va2/(diffspitz+sqrt(tini))
        dt1_va=sqrt(p%va2*dxyz_2)
      endif
!
      uplim=max(maxval(dt1_va),maxval(maxadvec))
      where (tau_inv_va > uplim)
        tau_inv_va=uplim
      endwhere
!
      where (tau_inv_va < tau_inv_spitzer)
        tau_inv_va=tau_inv_spitzer
      endwhere
    endif
!
    if (lnfs2) then
!     for pp=qq/rho, it is '+qq*(uglnrho+divu)'
      do i=1,3
        if (ltau_spitzer_va) then
          df(l1:l2,m,n,iqq+i-1) = df(l1:l2,m,n,iqq+i-1) - &
          tau_inv_va*(p%qq(:,i) + spitzer_vec(:,i))  +  &
          p%qq(:,i)*(p%uglnrho + p%divu)
        else
          df(l1:l2,m,n,iqq+i-1) = df(l1:l2,m,n,iqq+i-1) - &
          tau_inv_spitzer*(p%qq(:,i) + spitzer_vec(:,i))  +  &
          p%qq(:,i)*(p%uglnrho + p%divu)
        endif
      enddo
    else
!     for pp=qq*rho, it is '-qq*(uglnrho+divu)'
      do i=1,3
        df(l1:l2,m,n,iqq+i-1) = df(l1:l2,m,n,iqq+i-1) - &
          tau_inv_spitzer*(p%qq(:,i) + spitzer_vec(:,i))  -  &
          p%qq(:,i)*(p%uglnrho + p%divu)
      enddo
    endif
!
! Add to energy equation
!
    call dot(p%qq,p%glnrho,tmp)
!
    if (lnfs2) then
!     for pp=qq/rho, the 1/rho factor is vanishing
!     and there it is '+ tmp'
      rhs = gamma*p%cp1*(p%divq + tmp)*exp(-p%lnTT)    
    else
!     for pp=qq*rho, there is an additional 1/rho factor
!     and there it is '- tmp'
      rhs = gamma*p%cp1*(p%divq - tmp)*exp(-p%lnTT-2*p%lnrho)
    endif
!
!
    if (ltemperature) then
       if (ltemperature_nolog) then
          call fatal_error('non_fourier_spitzer', &
               'not implemented for current set of thermodynamic variables')
       else
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - rhs
       endif
    else if (lentropy.and.pretend_lnTT) then
       call fatal_error('non_fourier_spitzer', &
            'not implemented for current set of thermodynamic variables')
    else if (lthermal_energy .and. ldensity) then
       call fatal_error('non_fourier_spitzer', &
            'not implemented for current set of thermodynamic variables')
    else
       call fatal_error('non_fourier_spitzer', &
            'not implemented for current set of thermodynamic variables')
    endif
!
    if (lfirst.and.ldt) then
      if (ltau_spitzer_va) then
!
!       Define propagation speed c_spitzer
!       c_spitzer0 is upper limit for c_spitzer in the same way as
!       tau_inv_spitzer is the upper limit for tau_inv_va
!
        c_spitzer = sqrt(diffspitz*tau_inv_va)
        c_spitzer0 = sqrt(diffspitz*tau_inv_spitzer)
!
!       The advection time step should include c_spitzer, however, 
!       we have choosen c_spitzer to be sqrt(2) times Alfven speed (va),
!       so we can take this into account by multiplying c_spitzer 
!       by (1-sqrt(2)/2.). To be on the save side we use 0.36 .
!       If tau_inv_va gets so low that it reaches tau_inv_spitzer,
!       c_spitzer becomes larger than sqrt(2)*va. For this cases we need
!       to add the "full" value of c_spitzer to maxadvec.
!       For incomperating all we use:
!
        maxadvec = maxadvec + 0.36*c_spitzer/dxmin_pencil + 0.64*c_spitzer0/dxmin_pencil
!
!       In case tau_inv_va > tau_inv_spitzer is c_spitzer > c_spitzer0 and we get:
!       maxadvec = maxadvec + 0.36*c_spitzer/dxmin_pencil
!       In case tau_inv_va = tau_inv_spitzer is c_spitzer = c_spitzer0 and we get:
!       maxadvec = maxadvec + c_spitzer/dxmin_pencil
!
      else
        call unit_vector(p%glnTT,unit_glnTT)
        call dot(unit_glnTT,p%bunit,cosgT_b)
        diffspitz = Kspitzer_para*exp(2.5*p%lnTT-p%lnrho)* &
                   gamma*p%cp1*abs(cosgT_b)
        c_spitzer = sqrt(diffspitz*tau_inv_spitzer)
        maxadvec = maxadvec + c_spitzer/dxmin_pencil
      endif
!
      if (ldiagnos.and.idiag_dtq/=0) then
        call max_mn_name(c_spitzer/dxmin_pencil/cdt,idiag_dtq,l_dt=.true.)
      endif
!
!     put into dtspitzer, how the time_step would be
!     using the spitzer heatconductivity
!
      if (ldiagnos.and.idiag_dtspitzer/=0) then
        call max_mn_name(diffspitz*dxyz_2/cdtv,idiag_dtspitzer,l_dt=.true.)
      endif
!
!     timestep constraints due to tau directly
!
      if (ltau_spitzer_va) then
        dt1_max=max(dt1_max,maxval(tau_inv_va)/cdts)
      else
        dt1_max=max(dt1_max,tau_inv_spitzer/cdts)
      endif
!
      if (ldiagnos.and.idiag_dtq2/=0) then
        if (ltau_spitzer_va) then
          call max_mn_name(tau_inv_va/cdts,idiag_dtq2,l_dt=.true.)
          if (idiag_tauqmax/=0) then
            call max_mn_name(1./tau_inv_va,idiag_tauqmax)
          endif
        else
          call max_name(tau_inv_spitzer/cdts,idiag_dtq2,l_dt=.true.)
        endif
      endif
!
    endif
!
    if (lvideo.and.lfirst) then
      if (ivid_divq/=0) call store_slices((p%divq - tmp)*exp(-p%lnrho), &
        divq_xy,divq_xz,divq_yz,divq_xy2,divq_xy3,divq_xy4,divq_xz2)
    endif
!
  endsubroutine non_fourier_spitzer
!***********************************************************************
  subroutine eighth_moment_approx(df,p)
!
!  This heatconduction refers to the eighth moment approximation
!
!  07-sept-17/bingert: updated
!
    use Slices_methods, only: store_slices
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

   if (lvideo.and.lfirst) then
     if (ivid_divq/=0) call store_slices(rhs,divq_xy,divq_xz,divq_yz,divq_xy2,divq_xy3,divq_xy4,divq_xz2)
   endif
!
    if (lfirst.and.ldt) then
!      b2_1=1./(p%b2+tini)
      dt_1_8th = nu_coll !+ e_m /sqrt(b2_1)

      dt1_max=max(dt1_max,dt_1_8th/cdts)
!      advec_uu = max(advec_uu,advec_uu*7./5.)
      maxadvec = maxadvec*7./5.
!     advec_uu = max(advec_uu,sqrt(coeff*exp(p%lnTT)*gamma*p%cp1)/dxmax_pencil)
    endif
!
  endsubroutine eighth_moment_approx
!***********************************************************************
  subroutine noadvection_non_fourier_spitzer(df,p)
!
!  19-feb-18/piyali: adapted from non-fourier_spitzer for ionisation equation of
!  state, no advection term and a different free streaming limit
!
    use Slices_methods, only: store_slices
    use Diagnostics, only: max_mn_name,max_name
    use EquationOfState
    use SharedVariables, only: get_shared_variable
    use Sub
!
    real, dimension (mx,my,mz,mvar) :: df
    type (pencil_case) :: p
    real, dimension(nx) :: b2_1,tau1_spitzer_penc
    real, dimension(nx) :: chi_clight,chi_spitzer
    real, dimension(nx) :: rhs,cosgT_b
    real, dimension(nx,3) :: K1,unit_glnTT
    real, dimension(nx,3) :: spitzer_vec
    real, dimension(nx) :: tmp
    real, pointer :: z_cutoff, cool_wid
    integer :: ierr
    integer :: i
!
! Compute Spizter coefficiant K_0 * T^(5/2) * rho where
! we introduced an additional factor rho.
!
    b2_1=1./(p%b2+tini)
!
! Will only work when primary thermodynamic variables are lnrho, lnTT
!
    chi_spitzer=Kspitzer_para*exp(2.5*p%lnTT-p%lnrho)*p%cv1
!
!  Limit heat conduction so that the diffusion speed
!  is smaller than a given fraction of the speed of light (Kc*c_light)
!
    if (Kc /= 0.) then
      chi_clight = Kc * c_light/max(dz_1(n),dx_1(l1:l2))
      where (chi_spitzer > chi_clight)
        chi_spitzer = chi_clight
      endwhere
    endif
    call multsv(chi_spitzer*p%TT*p%rho/p%cv1,p%glnTT,K1)
    call dot(K1,p%bb,tmp)
    call multsv(b2_1*tmp,p%bb,spitzer_vec)
    tau1_spitzer_penc=cdtv**2/(1.0d-6*max(dz_1(n),dx_1(l1:l2)))**2/ &
                      chi_spitzer
    if (lradiation) then
      call get_shared_variable('z_cutoff',&
               z_cutoff,ierr)
      if (ierr/=0) call fatal_error('calc_heatcond_tensor:',&
             'failed to get z_cutoff from radiation_ray')
      call get_shared_variable('cool_wid',&
             cool_wid,ierr)
      if (ierr/=0) call fatal_error('calc_heatcond_tensor:',&
               'failed to get cool_wid from radiation_ray')
      do i=1,3
        df(l1:l2,m,n,iqq+i-1) = df(l1:l2,m,n,iqq+i-1) - &
            tau_inv_spitzer*(p%qq(:,i) + spitzer_vec(:,i))* &
            step(z(n),z_cutoff,cool_wid)
      enddo
    else
      do i=1,3
        df(l1:l2,m,n,iqq+i-1) = df(l1:l2,m,n,iqq+i-1) - &
            tau_inv_spitzer*(p%qq(:,i) + spitzer_vec(:,i)) 
      enddo
    endif
!
! Add to energy equation
!

    rhs = p%cv1*p%divq*p%TT1*p%rho1
!

    if (ltemperature) then
       if (ltemperature_nolog) then
          call fatal_error('non_fourier_spitzer', &
               'not implemented for current set of thermodynamic variables')
       else
          df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - rhs
          if (lfirst.and.ldt) then
            dt1_max=max(dt1_max,maxval(abs(rhs))/cdts)
          endif
       endif
    else if (lentropy.and.pretend_lnTT) then
       call fatal_error('noadvection_non_fourier_spitzer', &
            'not implemented for current set of thermodynamic variables')
    else if (lthermal_energy .and. ldensity) then
       call fatal_error('noadvection_non_fourier_spitzer', &
            'not implemented for current set of thermodynamic variables')
    else
       call fatal_error('noadvection_non_fourier_spitzer', &
            'not implemented for current set of thermodynamic variables')
    endif
!

    if (lfirst.and.ldt) then
      call unit_vector(p%glnTT,unit_glnTT)
      call dot(unit_glnTT,p%bunit,cosgT_b)
      rhs = Kspitzer_para*exp(2.5*p%lnTT-p%lnrho)*p%cv1
!
      if (idiag_dtspitzer/=0) &
           call max_mn_name(rhs/cdtv,idiag_dtspitzer,l_dt=.true.)
!
!
!     timestep constraints due to tau directly
!
      dt1_max=max(dt1_max,tau_inv_spitzer/cdts)
!
      if (ldiagnos.and.idiag_dtq2/=0) then
        call max_name(tau_inv_spitzer/cdts,idiag_dtq2,l_dt=.true.)
      endif
    endif
!
    if (lvideo.and.lfirst) then
      if (ivid_divq/=0) call store_slices(p%divq,divq_xy,divq_xz,divq_yz,divq_xy2,divq_xy3,divq_xy4,divq_xz2)
    endif
!
  endsubroutine noadvection_non_fourier_spitzer
!***********************************************************************
endmodule Heatflux
