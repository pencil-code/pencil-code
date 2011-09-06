! $Id$
!
!  This modules solves the passive scalar advection equation
!  Solves for c, not ln(c).
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpscalar = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED cc; cc1; gcc(3); ugcc; gcc2(3); gcc1(3);
! PENCILS PROVIDED del2cc; hcc(3,3); del6cc; g5cc(3); g5ccglnrho
!
!***************************************************************
module Pscalar
!
  use Cdata
  use Cparam
  use Messages
!
  implicit none
!
  include 'pscalar.h'
!
!  Init parameters.
!
  real, dimension(3) :: gradC0=(/0.0,0.0,0.0/)
  real :: ampllncc=impossible, widthlncc=impossible, lncc_min
  real :: ampllncc2=impossible, radius_lncc=impossible
  real :: kx_lncc=impossible, ky_lncc=impossible,kz_lncc=impossible
  real :: epsilon_lncc=impossible
  real :: amplcc=0.1, widthcc=0.5, cc_min=0.0
  real :: amplcc2=0.0, kx_cc=1.0, ky_cc=1.0, kz_cc=1.0, radius_cc=0.0
  real :: kxx_cc=0.0, kyy_cc=0.0, kzz_cc=0.0
  real :: epsilon_cc=0.0, cc_const=1.0
  real :: zoverh=1.0, hoverr=0.05, powerlr=3.0
  logical :: nopscalar=.false., reinitialize_cc=.false.
  logical :: reinitialize_lncc=.false.
  character (len=labellen) :: initlncc='impossible', initlncc2='impossible'
  character (len=labellen) :: initcc='zero', initcc2='zero'
  character (len=40) :: tensor_pscalar_file
!
  namelist /pscalar_init_pars/ &
      initcc, initcc2,amplcc, amplcc2, kx_cc, ky_cc, kz_cc, radius_cc, &
      epsilon_cc, widthcc, cc_min, cc_const, initlncc, initlncc2, ampllncc, &
      ampllncc2, kx_lncc, ky_lncc, kz_lncc, radius_lncc, epsilon_lncc, &
      widthlncc, kxx_cc, kyy_cc, kzz_cc, hoverr, powerlr, zoverh
!
!  Run parameters.
!
  real :: pscalar_diff=0.0, tensor_pscalar_diff=0.0, soret_diff=0.0
  real :: pscalar_diff_hyper3=0.0, rhoccm=0.0, cc2m=0.0, gcc2m=0.0
  real :: pscalar_sink=0.0, Rpscalar_sink=0.5
  real :: lam_gradC=0.0, om_gradC=0.0, lambda_cc=0.0
  real :: scalaracc=0.0
  real :: LLambda_cc=0.0
  logical :: lpscalar_sink, lgradC_profile=.false., lreactions=.false.
  logical :: lpscalar_diff_simple=.false.
  logical :: lpscalar_per_unitvolume=.false.
  logical :: lpscalar_per_unitvolume_diff=.false.
  logical :: lnotpassive=.false., lupw_cc=.false.
  logical :: lmean_friction_cc=.false.
!
  namelist /pscalar_run_pars/ &
      pscalar_diff, nopscalar, tensor_pscalar_diff, gradC0, soret_diff, &
      pscalar_diff_hyper3, reinitialize_lncc, reinitialize_cc, lpscalar_sink, &
      lmean_friction_cc, LLambda_cc, &
      lpscalar_diff_simple, &
      lpscalar_per_unitvolume, lpscalar_per_unitvolume_diff, &
      pscalar_sink, Rpscalar_sink, lreactions, lambda_cc, lam_gradC, &
      om_gradC, lgradC_profile, lnotpassive, lupw_cc
!
!  Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_rhoccm=0, idiag_ccmax=0, idiag_ccmin=0, idiag_ccm=0
  integer :: idiag_Qrhoccm=0, idiag_Qpsclm=0, idiag_mcct=0
  integer :: idiag_gcc5m=0, idiag_gcc10m=0
  integer :: idiag_ucm=0, idiag_uudcm=0, idiag_Cz2m=0, idiag_Cz4m=0
  integer :: idiag_Crmsm=0, idiag_uxcm=0, idiag_uycm=0, idiag_uzcm=0
  integer :: idiag_cc1m=0, idiag_cc2m=0, idiag_cc3m=0, idiag_cc4m=0
  integer :: idiag_cc5m=0, idiag_cc6m=0, idiag_cc7m=0, idiag_cc8m=0
  integer :: idiag_cc9m=0, idiag_cc10m=0
  integer :: idiag_gcc1m=0, idiag_gcc2m=0, idiag_gcc3m=0, idiag_gcc4m=0
  integer :: idiag_gcc6m=0, idiag_gcc7m=0, idiag_gcc8m=0, idiag_gcc9m=0
  integer :: idiag_ccmx=0, idiag_ccmy=0, idiag_ccmz=0, idiag_ccglnrm=0
  integer :: idiag_uxcmz=0, idiag_uycmz=0, idiag_uzcmz=0
  integer :: idiag_ccmxy=0, idiag_ccmxz=0
!
  contains
!***********************************************************************
    subroutine register_pscalar()
!
!  Initialise variables which should know that we solve for passive
!  scalar: icc; increase nvar accordingly
!
!  6-jul-02/axel: coded
!
      use FArrayManager
!
      lpscalar_nolog = .true.
!
      call farray_register_pde('cc',icc)
      ilncc = 0                 ! needed for idl
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_pscalar
!***********************************************************************
    subroutine initialize_pscalar(f)
!
!  Perform any necessary post-parameter read initialization
!  Since the passive scalar is often used for diagnostic purposes
!  one may want to reinitialize it to its initial distribution.
!
!  24-nov-02/tony: coded
!  20-may-03/axel: reinitialize_cc added
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if (reinitialize_cc) then
        f(:,:,:,icc)=0.
        call init_lncc(f)
      endif
!
      if (lnotpassive) scalaracc=3./5./hoverr**2
    endsubroutine initialize_pscalar
!***********************************************************************
    subroutine init_lncc(f)
!
!  initialise passive scalar field; called from start.f90
!
!   6-jul-2001/axel: coded
!
      use Sub
      use Initcond
      use InitialCondition, only: initial_condition_lncc
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      ! for the time being, keep old name for backward compatibility
      if (initlncc/='impossible') initcc=initlncc
      if (initlncc2/='impossible') initcc2=initlncc2
      if (ampllncc/=impossible) amplcc=ampllncc
      if (ampllncc2/=impossible) amplcc2=ampllncc2
      if (kx_lncc/=impossible) kx_cc=kx_lncc
      if (ky_lncc/=impossible) ky_cc=ky_lncc
      if (kz_lncc/=impossible) kz_cc=kz_lncc
      if (radius_lncc/=impossible) radius_cc=radius_lncc
      if (epsilon_lncc/=impossible) epsilon_cc=epsilon_lncc
      if (widthlncc/=impossible) widthcc=widthlncc
!
      select case (initcc)
        case ('nothing')
        case ('zero'); f(:,:,:,icc)=0.0
        case ('constant'); f(:,:,:,icc)=cc_const
        case ('hatwave-x'); call hatwave(amplcc,f,icc,widthcc,kx=kx_cc)
        case ('hatwave-y'); call hatwave(amplcc,f,icc,widthcc,ky=ky_cc)
        case ('hatwave-z'); call hatwave(amplcc,f,icc,widthcc,kz=kz_cc)
        case ('hat-x'); call hat(amplcc,f,icc,widthcc,kx=kx_cc)
        case ('hat-y'); call hat(amplcc,f,icc,widthcc,ky=ky_cc)
        case ('hat-z'); call hat(amplcc,f,icc,widthcc,kz=kz_cc)
        case ('gaussian-x'); call gaussian(amplcc,f,icc,kx=kx_cc)
        case ('gaussian-y'); call gaussian(amplcc,f,icc,ky=ky_cc)
        case ('gaussian-z'); call gaussian(amplcc,f,icc,kz=kz_cc)
        case ('parabola-x'); call parabola(amplcc,f,icc,kx=kx_cc)
        case ('parabola-y'); call parabola(amplcc,f,icc,ky=ky_cc)
        case ('parabola-z'); call parabola(amplcc,f,icc,kz=kz_cc)
        case ('gaussian-noise'); call gaunoise(amplcc,f,icc,icc)
        case ('wave-x'); call wave(amplcc,f,icc,kx=kx_cc)
        case ('wave-y'); call wave(amplcc,f,icc,ky=ky_cc)
        case ('wave-z'); call wave(amplcc,f,icc,kz=kz_cc)
        case ('linprof-x'); call linprof(amplcc,f,icc,kx=kx_cc)
        case ('linprof-y'); call linprof(amplcc,f,icc,ky=ky_cc)
        case ('linprof-z'); call linprof(amplcc,f,icc,kz=kz_cc)
        case ('propto-ux'); call wave_uu(amplcc,f,icc,kx=kx_cc)
        case ('propto-uy'); call wave_uu(amplcc,f,icc,ky=ky_cc)
        case ('propto-uz'); call wave_uu(amplcc,f,icc,kz=kz_cc)
        case ('cosx_cosy_cosz'); call cosx_cosy_cosz(amplcc,f,icc,kx_cc,ky_cc,kz_cc)
        case ('triquad'); call triquad(amplcc,f,icc,kx_cc,ky_cc,kz_cc, &
            kxx_cc,kyy_cc,kzz_cc)
        case ('semiangmom'); f(:,:,:,icc)=(1-2*powerlr*hoverr**2-1.5*zoverh**2*hoverr**2) &
            *spread(spread(x,2,my),3,mz) &
            +3*zoverh*hoverr*spread(spread(z,1,mx),2,my)
        case ('sound-wave')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,icc)=-amplcc*cos(kx_cc*x(l1:l2))
          enddo; enddo
        case ('tang-discont-z')
          print*,'init_lncc: widthcc=',widthcc
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,icc)=-1.0+2*.5*(1.+tanh(z(n)/widthcc))
          enddo; enddo
        case ('hor-tube'); call htube2(amplcc,f,icc,icc,radius_cc,epsilon_cc)
        case ('jump-x'); call jump(f,icc,cc_const,0.,widthcc,'x')
        case ('jump-x-neg'); call jump(f,icc,0.,cc_const,widthcc,'x')
        case ('jump-y-neg'); call jump(f,icc,0.,cc_const,widthcc,'y')
        case ('jump-z-neg'); call jump(f,icc,0.,cc_const,widthcc,'z')
        case ('jump'); call jump(f,icc,cc_const,0.,widthcc,'z')
        case default; call fatal_error('init_lncc','bad initcc='//trim(initcc))
      endselect
!
!  superimpose something else
!
      select case (initcc2)
        case ('wave-x'); call wave(amplcc2,f,icc,ky=5.)
        case ('constant'); f(:,:,:,icc)=f(:,:,:,icc)+amplcc2
      endselect
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_lncc(f)
!
!  add floor value if cc_min is set
!
      if (cc_min/=0.) then
        if (lroot) print*,'set floor value for cc; cc_min=',cc_min
        f(:,:,:,icc)=max(cc_min,f(:,:,:,icc))
      endif
!
    endsubroutine init_lncc
!***********************************************************************
    subroutine pencil_criteria_pscalar()
!
!  All pencils that the Pscalar module depends on are specified here.
!
!  20-nov-04/anders: coded
!
      integer :: i
!
      lpenc_requested(i_cc)=.true.
      if (.not. nopscalar) lpenc_requested(i_ugcc)=.true.
      if (lpscalar_per_unitvolume) then
        lpenc_requested(i_divu)=.true.
      endif
      if (lpscalar_per_unitvolume_diff) then
        lpenc_requested(i_del2lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gcc)=.true.
      endif
      if (lnotpassive) lpenc_requested(i_cc)=.true.
      if (lpscalar_sink) lpenc_requested(i_rho1)=.true.
      if (pscalar_diff/=0.) then
        lpenc_requested(i_gcc)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
      if (soret_diff/=0.) then
        lpenc_requested(i_cc)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
      if (lreactions) lpenc_requested(i_cc)=.true.
      do i=1,3
        if (gradC0(i)/=0.) lpenc_requested(i_uu)=.true.
      enddo
      if (pscalar_diff/=0.) lpenc_requested(i_del2cc)=.true.
      if (tensor_pscalar_diff/=0.) lpenc_requested(i_hcc)=.true.
      if (pscalar_diff_hyper3/=0.0) then
        lpenc_requested(i_del6cc)=.true.
        lpenc_requested(i_g5ccglnrho)=.true.
      endif
!
      lpenc_diagnos(i_cc)=.true.
      if (idiag_rhoccm/=0 .or. idiag_Cz2m/=0 .or. idiag_Cz4m/=0 .or. &
          idiag_Qrhoccm/=0 .or. idiag_Qpsclm/=0) &
          lpenc_diagnos(i_rho)=.true.
      if (idiag_ucm/=0 .or. idiag_uudcm/=0 .or. idiag_uxcm/=0 .or. &
          idiag_uycm/=0 .or. idiag_uzcm/=0 ) lpenc_diagnos(i_uu)=.true.
      if (idiag_uudcm/=0) lpenc_diagnos(i_ugcc)=.true.
      if (idiag_cc1m/=0 .or. idiag_cc2m/=0 .or. idiag_cc3m/=0 .or. &
          idiag_cc4m/=0 .or. idiag_cc5m/=0 .or. idiag_cc6m/=0 .or. &
          idiag_cc7m/=0 .or. idiag_cc8m/=0 .or. idiag_cc9m/=0 .or. &
          idiag_cc10m/=0) lpenc_diagnos(i_cc1)=.true.
      if (idiag_gcc1m/=0 .or. idiag_gcc2m/=0 .or. idiag_gcc3m/=0 .or. &
          idiag_gcc4m/=0 .or. idiag_gcc5m/=0 .or. idiag_gcc6m/=0 .or. &
          idiag_gcc7m/=0 .or. idiag_gcc8m/=0 .or. idiag_gcc9m/=0 .or. &
          idiag_gcc10m/=0) lpenc_diagnos(i_gcc1)=.true.
      if (idiag_ccglnrm/=0) lpenc_diagnos(i_glnrho)=.true.
!
      if (idiag_ccmxy/=0 .or. idiag_ccmxz/=0) lpenc_diagnos2d(i_cc)=.true.
!
    endsubroutine pencil_criteria_pscalar
!***********************************************************************
    subroutine pencil_interdep_pscalar(lpencil_in)
!
!  Interdependency among pencils provided by the Pscalar module
!  is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_cc1)) lpencil_in(i_cc)=.true.
      if (lpencil_in(i_ugcc)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gcc)=.true.
      endif
      if (lpencil_in(i_g5ccglnrho)) then
        lpencil_in(i_g5cc)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_gcc2)) lpencil_in(i_gcc)=.true.
      if (lpencil_in(i_gcc1)) lpencil_in(i_gcc2)=.true.
      if (tensor_pscalar_diff/=0.) lpencil_in(i_gcc)=.true.
!
    endsubroutine pencil_interdep_pscalar
!**********************************************************************
    subroutine calc_pencils_pscalar(f,p)
!
!  Calculate pscalar Pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
! cc
      if (lpencil(i_cc)) p%cc=f(l1:l2,m,n,icc)
! cc1
      if (lpencil(i_cc1)) p%cc1=1/p%cc
! gcc
      if (lpencil(i_gcc)) call grad(f,icc,p%gcc)
! ugcc
      if (lpencil(i_ugcc)) &
          call u_dot_grad(f,icc,p%gcc,p%uu,p%ugcc,UPWIND=lupw_cc)
! gcc2
      if (lpencil(i_gcc2)) call dot2_mn(p%gcc,p%gcc2)
! gcc1
      if (lpencil(i_gcc1)) p%gcc1=sqrt(p%gcc2)
! del2cc
      if (lpencil(i_del2cc)) call del2(f,icc,p%del2cc)
! hcc
      if (lpencil(i_hcc)) call g2ij(f,icc,p%hcc)
! del6cc
      if (lpencil(i_del6cc)) call del6(f,icc,p%del6cc)
! g5cc
      if (lpencil(i_g5cc)) call grad5(f,icc,p%g5cc)
! g5cc
      if (lpencil(i_g5ccglnrho)) then
        call dot_mn(p%g5cc,p%glnrho,p%g5ccglnrho)
      endif
!
    endsubroutine calc_pencils_pscalar
!***********************************************************************
    subroutine dlncc_dt(f,df,p)
!
!  Passive scalar evolution.
!  Calculate dc/dt=-uu.gcc + pscalar_diff*[del2cc + glnrho.gcc].
!
!  20-may-03/axel: coded
!
      use Diagnostics
      use Special, only: special_calc_pscalar
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: diff_op,diff_op2,bump,cc_xyaver
      real :: lam_gradC_fact=1., om_gradC_fact=1., gradC_fact=1.
      integer :: j
      integer, parameter :: nxy=nxgrid*nygrid
!
      intent(in)  :: f
      intent(out) :: df
!
!  Identify module and boundary conditions.
!
      if (nopscalar) then
        if (headtt.or.ldebug) print*,'not SOLVED: dlncc_dt'
      else
        if (headtt.or.ldebug) print*,'SOLVE dlncc_dt'
      endif
      if (headtt) call identify_bcs('cc',icc)
!
!  Gradient of passive scalar.
!  Allow for possibility to turn off passive scalar
!  without changing file size and recompiling everything.
!
      if (.not. nopscalar) then ! i.e. if (pscalar)
!
!  Passive scalar equation.
!
        df(l1:l2,m,n,icc)=df(l1:l2,m,n,icc) - p%ugcc
!
!  lpscalar_per_unitvolume
!
        if (lpscalar_per_unitvolume) then
          df(l1:l2,m,n,icc)=df(l1:l2,m,n,icc)-p%cc*p%divu
        endif
!
!  Reaction term. Simple Fisher term for now.
!
        if (lreactions) then
          df(l1:l2,m,n,icc)=df(l1:l2,m,n,icc)+lambda_cc*p%cc*(1.0-p%cc)
        endif
!
!  Passive scalar sink.
!
        if (lpscalar_sink) then
          bump=pscalar_sink* &
              exp(-0.5*(x(l1:l2)**2+y(m)**2+z(n)**2)/Rpscalar_sink**2)
          df(l1:l2,m,n,icc)=df(l1:l2,m,n,icc)-bump*f(l1:l2,m,n,icc)
        endif
!
!  Diffusion operator. If lpscalar_per_unitvolume is chosen, use
!  div[kappa*rho*grad(c/rho)] = kappa*(del2c-gradc.gradlnrho-c*del2lnrho).
!  Otherwise, with lpscalar_per_unitvolume and not lpscalar_per_unitvolume_diff
!  use just kappa*del2c.
!
        if (pscalar_diff/=0.) then
          if (headtt) print*,'dlncc_dt: pscalar_diff=',pscalar_diff
          if (lpscalar_per_unitvolume) then
            if (lpscalar_per_unitvolume_diff) then
              call dot_mn(p%glnrho,p%gcc,diff_op)
              diff_op=p%del2cc-diff_op-p%cc*p%del2lnrho
            else
              diff_op=p%del2cc
            endif
          else
            if (lpscalar_diff_simple) then
              diff_op=p%del2cc
            else
              call dot_mn(p%glnrho,p%gcc,diff_op)
              diff_op=diff_op+p%del2cc
            endif
          endif
          df(l1:l2,m,n,icc)=df(l1:l2,m,n,icc)+pscalar_diff*diff_op
        endif
!
!  Hyperdiffusion operator.
!
        if (pscalar_diff_hyper3/=0.) then
          if (headtt) &
              print*,'dlncc_dt: pscalar_diff_hyper3=', pscalar_diff_hyper3
          df(l1:l2,m,n,icc) = df(l1:l2,m,n,icc) + &
              pscalar_diff_hyper3*(p%del6cc+p%g5ccglnrho)
        endif
!
!  Soret diffusion.
!
        if (soret_diff/=0.) then
          if (headtt) print*,'dlncc_dt: soret_diff=',soret_diff
          call dot2_mn(p%glnTT,diff_op2)
          diff_op2=p%cc*(1.-p%cc)*p%TT*(diff_op2+p%del2lnTT)
          df(l1:l2,m,n,icc) = df(l1:l2,m,n,icc) + soret_diff*diff_op2
        endif
!
!  Time-dependent prefactor for the imposed passive scalar gradient.
!
        if (lam_gradC/=0..or.om_gradC/=0.) then
          if (lam_gradC/=0.) lam_gradC_fact=exp(lam_gradC*t)
          if ( om_gradC/=0.)  om_gradC_fact=cos( om_gradC*t)
          gradC_fact=lam_gradC_fact*om_gradC_fact
        endif
!
!  Possibility of an additional z-profile on gradC_fact.
!
        if (lgradC_profile) then
          gradC_fact=gradC_fact*cos(z(n))
        endif
!
!  Add diffusion of imposed spatially constant gradient of c.
!  This makes sense really only for periodic boundary conditions.
!
        do j=1,3
          if (gradC0(j)/=0.) then
            df(l1:l2,m,n,icc)=df(l1:l2,m,n,icc)-gradC0(j)*p%uu(:,j)*gradC_fact
          endif
        enddo
!
!  Tensor diffusion (but keep the isotropic one).
!
        if (tensor_pscalar_diff/=0.) &
            call tensor_diff(df,p,tensor_pscalar_diff)
!
!  Consider here the action of a mean friction term, -LLambda*Cbar.
!  This can be used to compendate for the decay of a horizontally
!  averaged mean concentration and allows thus the determination of
!  the turbulent diffusivity under stationary conditions. Only those
!  results are then comparable with the results of the test-field method.
!
      if (lmean_friction_cc) then
        if (nprocx*nprocy==1) then
          cc_xyaver=sum(f(l1:l2,m1:m2,n,icc))/nxy
          df(l1:l2,m,n,icc)=df(l1:l2,m,n,icc)-LLambda_cc*cc_xyaver
        else
          call fatal_error('pscalar','lmean_friction works only for nprocxy=1')
        endif
      endif
!
!  For the timestep calculation, need maximum diffusion.
!
        if (lfirst.and.ldt) then
          diffus_pscalar =(pscalar_diff+tensor_pscalar_diff)*dxyz_2
          diffus_pscalar3=pscalar_diff_hyper3*dxyz_6
        endif
!
!  Special contributions to this module are called here.
!
        if (lspecial) call special_calc_pscalar(f,df,p)
!
      endif
!
!  Diagnostics.
!
!  output for double and triple correlators (assume z-gradient of cc)
!  <u_k u_j d_j c> = <u_k c uu.gradcc>
!
      if (ldiagnos) then
        if (idiag_Qpsclm/=0)  call sum_mn_name(bump,idiag_Qpsclm)
        if (idiag_Qrhoccm/=0) call sum_mn_name(bump*p%rho*p%cc,idiag_Qrhoccm)
        if (idiag_mcct/=0)    call integrate_mn_name(p%rho*p%cc,idiag_mcct)
        if (idiag_rhoccm/=0)  call sum_mn_name(p%rho*p%cc,idiag_rhoccm)
        if (idiag_ccmax/=0)   call max_mn_name(p%cc,idiag_ccmax)
        if (idiag_ccmin/=0)   call max_mn_name(-p%cc,idiag_ccmin,lneg=.true.)
        if (idiag_uxcm/=0)    call sum_mn_name(p%uu(:,1)*p%cc,idiag_uxcm)
        if (idiag_uycm/=0)    call sum_mn_name(p%uu(:,2)*p%cc,idiag_uycm)
        if (idiag_uzcm/=0)    call sum_mn_name(p%uu(:,3)*p%cc,idiag_uzcm)
        if (lgradC_profile) then
          if (idiag_ucm/=0)   call sum_mn_name(2.*cos(z(n))*p%uu(:,3)*p%cc,idiag_ucm)
        else
          if (idiag_ucm/=0)   call sum_mn_name(p%uu(:,3)*p%cc,idiag_ucm)
        endif
        if (idiag_uudcm/=0)   call sum_mn_name(p%uu(:,3)*p%ugcc,idiag_uudcm)
        if (idiag_Cz2m/=0)    call sum_mn_name(p%rho*p%cc*z(n)**2,idiag_Cz2m)
        if (idiag_Cz4m/=0)    call sum_mn_name(p%rho*p%cc*z(n)**4,idiag_Cz4m)
        if (idiag_Crmsm/=0) &
            call sum_mn_name((p%rho*p%cc)**2,idiag_Crmsm,lsqrt=.true.)
        if (idiag_cc1m/=0)    call sum_mn_name(p%cc1   ,idiag_cc1m)
        if (idiag_cc2m/=0)    call sum_mn_name(p%cc1**2,idiag_cc2m)
        if (idiag_cc3m/=0)    call sum_mn_name(p%cc1**3,idiag_cc3m)
        if (idiag_cc4m/=0)    call sum_mn_name(p%cc1**4,idiag_cc4m)
        if (idiag_cc5m/=0)    call sum_mn_name(p%cc1**5,idiag_cc5m)
        if (idiag_cc6m/=0)    call sum_mn_name(p%cc1**6,idiag_cc6m)
        if (idiag_cc7m/=0)    call sum_mn_name(p%cc1**7,idiag_cc7m)
        if (idiag_cc8m/=0)    call sum_mn_name(p%cc1**8,idiag_cc8m)
        if (idiag_cc9m/=0)    call sum_mn_name(p%cc1**9,idiag_cc9m)
        if (idiag_cc10m/=0)   call sum_mn_name(p%cc1**10,idiag_cc10m)
        if (idiag_gcc1m/=0)   call sum_mn_name(p%gcc1   ,idiag_gcc1m)
        if (idiag_gcc2m/=0)   call sum_mn_name(p%gcc1**2,idiag_gcc2m)
        if (idiag_gcc3m/=0)   call sum_mn_name(p%gcc1**3,idiag_gcc3m)
        if (idiag_gcc4m/=0)   call sum_mn_name(p%gcc1**4,idiag_gcc4m)
        if (idiag_gcc5m/=0)   call sum_mn_name(p%gcc1**5,idiag_gcc5m)
        if (idiag_gcc6m/=0)   call sum_mn_name(p%gcc1**6,idiag_gcc6m)
        if (idiag_gcc7m/=0)   call sum_mn_name(p%gcc1**7,idiag_gcc7m)
        if (idiag_gcc8m/=0)   call sum_mn_name(p%gcc1**8,idiag_gcc8m)
        if (idiag_gcc9m/=0)   call sum_mn_name(p%gcc1**9,idiag_gcc9m)
        if (idiag_gcc10m/=0)  call sum_mn_name(p%gcc1**10,idiag_gcc10m)
        if (idiag_ccglnrm/=0) call sum_mn_name(p%cc*p%glnrho(:,3),idiag_ccglnrm)
      endif
!
      if (l1davgfirst) then
        call xysum_mn_name_z(p%cc,idiag_ccmz)
        if (idiag_ccmy/=0)    call xzsum_mn_name_y(p%cc,idiag_ccmy)
        if (idiag_ccmx/=0)    call yzsum_mn_name_x(p%cc,idiag_ccmx)
        call xysum_mn_name_z(p%uu(:,1)*p%cc,idiag_uxcmz)
        call xysum_mn_name_z(p%uu(:,2)*p%cc,idiag_uycmz)
        call xysum_mn_name_z(p%uu(:,3)*p%cc,idiag_uzcmz)
      endif
!
      if (l2davgfirst) then
        if (idiag_ccmxy/=0)   call zsum_mn_name_xy(p%cc,idiag_ccmxy)
        if (idiag_ccmxz/=0)   call ysum_mn_name_xz(p%cc,idiag_ccmxz)
      endif
!
! AH: notpassive, an angular momentum+gravity workaround
      if (lnotpassive) then
        if (lhydro) then
          df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+(p%cc &
             +(-powerlr*hoverr**2+1.5*zoverh**2*hoverr**2)+ &
             (-1+3*powerlr*hoverr**2-4.5*zoverh**2*hoverr**2)*x(l1:l2))*scalaracc
          df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)+(-hoverr*zoverh &
              -z(n)+3*hoverr*zoverh*x(l1:l2))*scalaracc
        endif
      endif
!
    endsubroutine dlncc_dt
!***********************************************************************
    subroutine read_pscalar_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=pscalar_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=pscalar_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_pscalar_init_pars
!***********************************************************************
    subroutine write_pscalar_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=pscalar_init_pars)
!
    endsubroutine write_pscalar_init_pars
!***********************************************************************
    subroutine read_pscalar_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=pscalar_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=pscalar_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_pscalar_run_pars
!***********************************************************************
    subroutine write_pscalar_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=pscalar_run_pars)
!
    endsubroutine write_pscalar_run_pars
!***********************************************************************
    subroutine rprint_pscalar(lreset,lwrite)
!
!  Reads and registers print parameters relevant for passive scalar.
!
!   6-jul-02/axel: coded
!
      use Diagnostics
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamez, inamey, inamex, inamexy, inamexz
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_rhoccm=0; idiag_ccmax=0; idiag_ccmin=0.; idiag_ccm=0
        idiag_Qrhoccm=0; idiag_Qpsclm=0; idiag_mcct=0
        idiag_ccmz=0; idiag_ccmy=0; idiag_ccmx=0
        idiag_uxcmz=0; idiag_uycmz=0; idiag_uzcmz=0
        idiag_ucm=0; idiag_uudcm=0; idiag_Cz2m=0; idiag_Cz4m=0; idiag_Crmsm=0
        idiag_uxcm=0; idiag_uycm=0; idiag_uzcm=0
        idiag_cc1m=0; idiag_cc2m=0; idiag_cc3m=0; idiag_cc4m=0; idiag_cc5m=0
        idiag_cc6m=0; idiag_cc7m=0; idiag_cc8m=0; idiag_cc9m=0; idiag_cc10m=0
        idiag_gcc1m=0; idiag_gcc2m=0; idiag_gcc3m=0; idiag_gcc4m=0
        idiag_gcc5m=0; idiag_gcc6m=0; idiag_gcc7m=0; idiag_gcc8m=0
        idiag_gcc9m=0; idiag_gcc10m=0; idiag_ccglnrm=0
        idiag_ccmxy=0; idiag_ccmxz=0
      endif
!
!  Check for those quantities that we want to evaluate online.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Qpsclm',idiag_Qpsclm)
        call parse_name(iname,cname(iname),cform(iname),'Qrhoccm',idiag_Qrhoccm)
        call parse_name(iname,cname(iname),cform(iname),'rhoccm',idiag_rhoccm)
        call parse_name(iname,cname(iname),cform(iname),'mcct',idiag_mcct)
        call parse_name(iname,cname(iname),cform(iname),'ccmax',idiag_ccmax)
        call parse_name(iname,cname(iname),cform(iname),'ccmin',idiag_ccmin)
        call parse_name(iname,cname(iname),cform(iname),'ccm',idiag_ccm)
        call parse_name(iname,cname(iname),cform(iname),'ucm',idiag_ucm)
        call parse_name(iname,cname(iname),cform(iname),'uxcm',idiag_uxcm)
        call parse_name(iname,cname(iname),cform(iname),'uycm',idiag_uycm)
        call parse_name(iname,cname(iname),cform(iname),'uzcm',idiag_uzcm)
        call parse_name(iname,cname(iname),cform(iname),'uudcm',idiag_uudcm)
        call parse_name(iname,cname(iname),cform(iname),'Cz2m',idiag_Cz2m)
        call parse_name(iname,cname(iname),cform(iname),'Cz4m',idiag_Cz4m)
        call parse_name(iname,cname(iname),cform(iname),'Crmsm',idiag_Crmsm)
        call parse_name(iname,cname(iname),cform(iname),'cc1m',idiag_cc1m)
        call parse_name(iname,cname(iname),cform(iname),'cc2m',idiag_cc2m)
        call parse_name(iname,cname(iname),cform(iname),'cc3m',idiag_cc3m)
        call parse_name(iname,cname(iname),cform(iname),'cc4m',idiag_cc4m)
        call parse_name(iname,cname(iname),cform(iname),'cc5m',idiag_cc5m)
        call parse_name(iname,cname(iname),cform(iname),'cc6m',idiag_cc6m)
        call parse_name(iname,cname(iname),cform(iname),'cc7m',idiag_cc7m)
        call parse_name(iname,cname(iname),cform(iname),'cc8m',idiag_cc8m)
        call parse_name(iname,cname(iname),cform(iname),'cc9m',idiag_cc9m)
        call parse_name(iname,cname(iname),cform(iname),'cc10m',idiag_cc10m)
        call parse_name(iname,cname(iname),cform(iname),'gcc1m',idiag_gcc1m)
        call parse_name(iname,cname(iname),cform(iname),'gcc2m',idiag_gcc2m)
        call parse_name(iname,cname(iname),cform(iname),'gcc3m',idiag_gcc3m)
        call parse_name(iname,cname(iname),cform(iname),'gcc4m',idiag_gcc4m)
        call parse_name(iname,cname(iname),cform(iname),'gcc5m',idiag_gcc5m)
        call parse_name(iname,cname(iname),cform(iname),'gcc6m',idiag_gcc6m)
        call parse_name(iname,cname(iname),cform(iname),'gcc7m',idiag_gcc7m)
        call parse_name(iname,cname(iname),cform(iname),'gcc8m',idiag_gcc8m)
        call parse_name(iname,cname(iname),cform(iname),'gcc9m',idiag_gcc9m)
        call parse_name(iname,cname(iname),cform(iname),'gcc10m',idiag_gcc10m)
        call parse_name(iname,cname(iname),cform(iname),'ccglnrm',idiag_ccglnrm)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxcmz',idiag_uxcmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uycmz',idiag_uycmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzcmz',idiag_uzcmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ccmz',idiag_ccmz)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'ccmy',idiag_ccmy)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ccmx',idiag_ccmx)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'ccmxz',idiag_ccmxz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'ccmxy',idiag_ccmxy)
      enddo
!
!  Write column where which passive scalar variable is stored.
!
      if (lwr) then
        write(3,*) 'ilncc=0'
        write(3,*) 'icc=',icc
      endif
!
    endsubroutine rprint_pscalar
!***********************************************************************
    subroutine get_slices_pscalar(f,slices)
!
!  Write slices for animation of pscalar variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices.
!
      select case (trim(slices%name))
!
!  Passive scalar.
!
        case ('cc')
          slices%yz =f(ix_loc,m1:m2,n1:n2,icc)
          slices%xz =f(l1:l2,iy_loc,n1:n2,icc)
          slices%xy =f(l1:l2,m1:m2,iz_loc,icc)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,icc)
          if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,icc)
          if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,icc)
          slices%ready=.true.
!
!  Logarithmic passive scalar.
!
        case ('lncc')
          slices%yz =alog(f(ix_loc,m1:m2,n1:n2,icc))
          slices%xz =alog(f(l1:l2,iy_loc,n1:n2,icc))
          slices%xy =alog(f(l1:l2,m1:m2,iz_loc,icc))
          slices%xy2=alog(f(l1:l2,m1:m2,iz2_loc,icc))
          if (lwrite_slice_xy3) slices%xy3=alog(f(l1:l2,m1:m2,iz3_loc,icc))
          if (lwrite_slice_xy4) slices%xy4=alog(f(l1:l2,m1:m2,iz4_loc,icc))
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_pscalar
!***********************************************************************
    subroutine calc_mpscalar
!
!  Calculate mean magnetic field from xy- or z-averages.
!
!  14-apr-03/axel: adaped from calc_mfield
!
      use Diagnostics
!
      logical,save :: first=.true.
      real :: ccm
!
!  Magnetic energy in horizontally averaged field
!  The bxmz and bymz must have been calculated,
!  so they are present on the root processor.
!
      if (idiag_ccm/=0) then
        if (idiag_ccmz==0) then
          if (first) print*
          if (first) print*,"NOTE: to get ccm, ccmz must also be set in xyaver"
          if (first) print*,"      We proceed, but you'll get ccm=0"
          ccm=0.
        else
          ccm=sqrt(sum(fnamez(:,:,idiag_ccmz)**2)/(nz*nprocz))
        endif
        call save_name(ccm,idiag_ccm)
      endif
!
    endsubroutine calc_mpscalar
!***********************************************************************
    subroutine tensor_diff(df,p,tensor_pscalar_diff)
!
!  Reads file.
!
!  11-jul-02/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: tensor_pscalar_diff
!
      real, save, dimension (nx,ny,nz,3) :: bunit,hhh
      real, dimension (nx) :: tmp,scr
      integer :: iy,iz,i,j
      logical, save :: first=.true.
!
!  read H and Bunit arrays and keep them in memory
!
      if (first) then
        open(1,file=trim(directory)//'/bunit.dat',form='unformatted')
        print*,'read bunit.dat with dimension: ',nx,ny,nz,3
        read(1) bunit,hhh
        close(1)
        print*,'read bunit.dat; bunit=',bunit
      endif
!
!  tmp = (Bunit.G)^2 + H.G + Bi*Bj*Gij
!  for details, see tex/mhd/thcond/tensor_der.tex
!
      call dot_mn(bunit,p%gcc,scr)
      call dot_mn(hhh,p%gcc,tmp)
      tmp=tmp+scr**2
!
!  dot with bi*bj
!
      iy=m-m1+1
      iz=n-n1+1
      do j=1,3
      do i=1,3
        tmp=tmp+bunit(:,iy,iz,i)*bunit(:,iy,iz,j)*p%hcc(:,i,j)
      enddo
      enddo
!
!  and add result to the dcc/dt equation
!
      df(l1:l2,m,n,icc)=df(l1:l2,m,n,icc)+tensor_pscalar_diff*tmp
!
      first=.false.
!
    endsubroutine tensor_diff
!***********************************************************************
endmodule Pscalar
