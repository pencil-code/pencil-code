! $Id$
!
!  This modules solves (multiple) passive scalar advection equation(s)
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
! PENCILS PROVIDED cc(npscalar); cc1(npscalar)
! PENCILS PROVIDED gcc(3,npscalar); ugcc(npscalar)
! PENCILS PROVIDED gcc2(npscalar); gcc1(npscalar)
! PENCILS PROVIDED del2cc(npscalar); del6cc(npscalar)
! PENCILS PROVIDED g5cc(3,npscalar); g5ccglnrho(npscalar)
! PENCILS PROVIDED hcc(3,3,npscalar)
!***************************************************************
module Pscalar
!
  use Cdata
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
  real :: cc_left=1., cc_right=0.
  real :: amplcc=0.1, widthcc=0.5, cc_min=0.0
  real :: amplcc2=0.0, kx_cc=1.0, ky_cc=1.0, kz_cc=1.0, radius_cc=0.0
  real :: kxx_cc=0.0, kyy_cc=0.0, kzz_cc=0.0
  real :: epsilon_cc=0.0, cc_const=1.0
  real :: zoverh=1.0, hoverr=0.05, powerlr=3.0
  real :: xjump_mid=0.,yjump_mid=0.,zjump_mid=0.
  logical :: nopscalar=.false., reinitialize_cc=.false.
  logical :: reinitialize_lncc=.false.
  character (len=labellen) :: initlncc='impossible', initlncc2='impossible'
  character (len=labellen) :: initcc='nothing', initcc2='zero'
  character (len=40) :: tensor_pscalar_file
  integer :: ll_sh=-1, mm_sh=-1, n_xprof=0
!
  namelist /pscalar_init_pars/ &
      initcc, initcc2,amplcc, amplcc2, kx_cc, ky_cc, kz_cc, radius_cc, &
      cc_left, cc_right, &
      epsilon_cc, widthcc, cc_min, cc_const, initlncc, initlncc2, ampllncc, &
      ampllncc2, kx_lncc, ky_lncc, kz_lncc, radius_lncc, epsilon_lncc, &
      widthlncc, kxx_cc, kyy_cc, kzz_cc, hoverr, powerlr, zoverh, ll_sh, mm_sh, n_xprof, &
      xjump_mid, yjump_mid, zjump_mid 
!
!  Run parameters.
!
  real :: pscalar_diff=0.0, tensor_pscalar_diff=0.0, soret_diff=0.0
  real :: diffcc_shock = 0.
  real :: pscalar_diff_hyper3=0.0
  real :: rhoccm=0.0, cc2m=0.0, gcc2m=0.0
  real :: pscalar_sink=0.0, Rpscalar_sink=0.5
  real :: lam_gradC=0.0, om_gradC=0.0, lambda_cc=0.0
  real :: scalaracc=0.0
  real :: LLambda_cc=0.0
  logical :: lpscalar_sink=.false., lgradC_profile=.false., lreactions=.false.
  logical :: lpscalar_diff_simple=.false.
  logical :: lpscalar_per_unitvolume=.false.
  logical :: lpscalar_per_unitvolume_diff=.false.
  logical :: lnotpassive=.false., lupw_cc=.false.
  logical :: lmean_friction_cc=.false.
  logical :: lremove_mean=.false.
!
  namelist /pscalar_run_pars/ &
      nopscalar, tensor_pscalar_diff, gradC0, soret_diff, &
      pscalar_diff, diffcc_shock, pscalar_diff_hyper3, &
      reinitialize_lncc, reinitialize_cc, lpscalar_sink, &
      lmean_friction_cc, LLambda_cc, &
      lpscalar_diff_simple, &
      lpscalar_per_unitvolume, lpscalar_per_unitvolume_diff, &
      pscalar_sink, Rpscalar_sink, lreactions, lambda_cc, lam_gradC, &
      om_gradC, lgradC_profile, lnotpassive, lupw_cc, lremove_mean, ll_sh, mm_sh
!
!  Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_rhoccm=0, idiag_ccmax=0, idiag_ccmin=0, idiag_ccm=0
  integer :: idiag_mrclncm=0
  integer :: idiag_Qrhoccm=0, idiag_Qpsclm=0, idiag_mcct=0
  integer :: idiag_gcc5m=0, idiag_gcc10m=0
  integer :: idiag_ucm=0, idiag_uudcm=0, idiag_Cz2m=0, idiag_Cz4m=0
  integer :: idiag_Crmsm=0, idiag_ccrms=0
  integer :: idiag_uxcm=0, idiag_uycm=0, idiag_uzcm=0
  integer :: idiag_cc1m=0, idiag_cc2m=0, idiag_cc3m=0, idiag_cc4m=0
  integer :: idiag_cc5m=0, idiag_cc6m=0, idiag_cc7m=0, idiag_cc8m=0
  integer :: idiag_cc9m=0, idiag_cc10m=0
  integer :: idiag_gcc1m=0, idiag_gcc2m=0, idiag_gcc3m=0, idiag_gcc4m=0
  integer :: idiag_gcc6m=0, idiag_gcc7m=0, idiag_gcc8m=0, idiag_gcc9m=0
  integer :: idiag_cugccm=0, idiag_ccugum=0
  integer :: idiag_ccmx=0, idiag_ccmy=0, idiag_ccmz=0, idiag_ccglnrm=0
  integer :: idiag_uxcmz=0, idiag_uycmz=0, idiag_uzcmz=0, idiag_cc2mz=0
  integer :: idiag_ccmxy=0, idiag_ccmxz=0
  integer :: idiag_cluz_uzlcm=0, idiag_gcguzm=0
!
  real, dimension(:,:), allocatable :: spharm 
  real, dimension(:,:,:,:), allocatable :: bunit,hhh
  real, dimension (nx) :: bump

  contains
!***********************************************************************
    subroutine register_pscalar
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
      call farray_register_pde('cc', icc, array=npscalar)
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
!         17/MR  : added bump for lpscalar_sink defined by spherical harmonics,
!                  specified by ll_sh, mm_sh >=0.
!
      use Sub, only: ylm_other, ylm
      use General, only: yin2yang_coors

      real, dimension (mx,my,mz,mfarray) :: f

      real, dimension(:,:), allocatable :: yz
      integer :: iyz
      real :: sphder
!
!  Print the number of passive scalars.
!
      if (lroot) print*, 'Number of passive scalars = ', npscalar
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if (reinitialize_cc) then
        f(:,:,:,icc)=0.
        call init_lncc(f)
      endif
!
      if (lroot .and. diffcc_shock /= 0.)  &
        print*,'initialize_pscalar: shock diffusion, diffcc_shock = ',diffcc_shock
!
      if (lmean_friction_cc.and.nprocxy/=1) &
        call fatal_error('initialize_pscalar','lmean_friction works only for nprocxy=1')
          
      if (lnotpassive) scalaracc=3./5./hoverr**2
!
      if (tensor_pscalar_diff/=0.) then
        if (.not.allocated(bunit)) allocate(bunit(nx,ny,nz,3),hhh(nx,ny,nz,3))
!
!  read H and Bunit arrays and keep them in memory
!
        open(1,file=trim(directory)//'/bunit.dat',form='unformatted')  !MR: better only root!
        print*,'read bunit.dat with dimension: ',nx,ny,nz,3
        read(1) bunit,hhh
        close(1)
      endif
!
      if (lyinyang) then
        if (lpscalar_sink.and.Rpscalar_sink/=0) &
          call warning('initialize_pscalar','passive scalar sink not correctly implemented on Yin-Yang grid') 
        if (lnotpassive) &
          call warning('initialize_pscalar','lnotpassive=T not correctly implemented on Yin-Yang grid') 
        if (lgradC_profile) &
          call warning('initialize_pscalar','lgradC_profile=T not correctly implemented on Yin-Yang grid') 
      endif

      if (lpscalar_sink.and.pscalar_sink/=0..and.ll_sh>=0.and.mm_sh>=0) then

        allocate(spharm(ny,nz))
        if (lyang) then
          allocate(yz(2,ny*nz))
          call yin2yang_coors(costh(m1:m2),sinth(m1:m2),cosph(n1:n2),sinph(n1:n2),yz)
          iyz=1
          do m=m1,m2
            do n=n1,n2
              spharm(m-m1+1,n-n1+1)=ylm_other(yz(1,iyz),yz(2,iyz),ll_sh,mm_sh,sphder)
              iyz=iyz+1
            enddo
          enddo
        else
          do n=n1,n2
            do m=m1,m2
              spharm(m-m1+1,n-n1+1)=ylm(ll_sh,mm_sh,sphder)
            enddo
          enddo
        endif
 
      endif
!
    endsubroutine initialize_pscalar
!***********************************************************************
    subroutine init_lncc(f)
!
!  initialise passive scalar field; called from start.f90
!
!   6-jul-2001/axel: coded
!  27-jun-2017/MR: added initial condition: spherical harmonic in \theta,\phi,
!                  sinusoidal in r with frequency n_xprof
!
      use Sub
      use General, only: yin2yang_coors
      use Initcond
      use InitialCondition, only: initial_condition_lncc
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx) :: tmpx
      real, dimension(:,:), allocatable :: yz
      integer :: iyz
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
!  Warning for the case of multiple passive scalars.
!
      if (npscalar > 1 .and. .not. linitial_condition) call warning('init_lncc', &
        'only the first species is initialized - use initial_condition facility instead')
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
        case ('wave-z-4') 
           call wave(amplcc/4.,f,icc,kz=1.)
           call wave(amplcc/4.,f,icc,kz=2.)
           call wave(amplcc/4.,f,icc,kz=4.)
           call wave(amplcc/4.,f,icc,kz=6.)
        case ('linprof-x'); call linprof(amplcc,f,icc,kx=kx_cc)
        case ('linprof-y'); call linprof(amplcc,f,icc,ky=ky_cc)
        case ('linprof-z'); call linprof(amplcc,f,icc,kz=kz_cc)
        case ('propto-ux'); call wave_uu(amplcc,f,icc,kx=kx_cc)
        case ('propto-uy'); call wave_uu(amplcc,f,icc,ky=ky_cc)
        case ('propto-uz'); call wave_uu(amplcc,f,icc,kz=kz_cc)
        case ('cosx_cosy_cosz'); call cosx_cosy_cosz(amplcc,f,icc,kx_cc,ky_cc,kz_cc)
        case ('triquad'); call triquad(amplcc,f,icc,kx_cc,ky_cc,kz_cc,kxx_cc,kyy_cc,kzz_cc)
        case ('semiangmom'); f(:,:,:,icc)=(1-2*powerlr*hoverr**2-1.5*zoverh**2*hoverr**2) &
                             *spread(spread(x,2,my),3,mz) + 3*zoverh*hoverr*spread(spread(z,1,mx),2,my)
        case ('sound-wave')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,icc)=-amplcc*cos(kx_cc*x(l1:l2))
          enddo; enddo
       case ('Schur')
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,icc)=-amplcc*sin(kx_cc*x(l1:l2))*sin(ky_cc*y(m))*cos(kz_cc*z(n))
          enddo; enddo
        case ('tang-discont-z')
          print*,'init_lncc: widthcc=',widthcc
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,icc)=-1.0+2*.5*(1.+tanh(z(n)/widthcc))
          enddo; enddo
        case ('hor-tube'); call htube2(amplcc,f,icc,icc,radius_cc,epsilon_cc)
        case ('jump-x'); call jump(f,icc,cc_left,cc_right,widthcc,xjump_mid,yjump_mid,zjump_mid,'x')
        case ('jump-x-neg'); call jump(f,icc,0.,cc_const,widthcc,xjump_mid,yjump_mid,zjump_mid,'x')
        case ('jump-y-neg'); call jump(f,icc,0.,cc_const,widthcc,xjump_mid,yjump_mid,zjump_mid,'y')
        case ('jump-z-neg'); call jump(f,icc,0.,cc_const,widthcc,xjump_mid,yjump_mid,zjump_mid,'z')
        case ('jump'); call jump(f,icc,cc_const,0.,widthcc,xjump_mid,yjump_mid,zjump_mid,'z')
        case('spher-harm')
          if (.not.lspherical_coords) call fatal_error("init_lncc", &
              'initial condition "spher-harm" only meaningful for spherical coordinates'//trim(initcc))
          !tmpx=(x(l1:l2)-xyz0(1))*(x(l1:l2)-xyz1(1)) + (xyz1(1) - 0.5*xyz0(1))*x(l1:l2)         ! S(r)
          tmpx=sin((2.*pi/(Lxyz(1))*n_xprof)*(x(l1:l2)-xyz0(1)))

          if (lyang) then
            allocate(yz(2,ny*nz))
            call yin2yang_coors(costh(m1:m2),sinth(m1:m2),cosph(n1:n2),sinph(n1:n2),yz)
            iyz=1
            do m=m1,m2
              do n=n1,n2
!if (iproc_world==55) print*, 'm,n,yz=', m,n,yz(:,iyz)
                f(l1:l2,m,n,icc) = amplcc*tmpx*ylm_other(yz(1,iyz),yz(2,iyz),ll_sh,mm_sh)
                iyz=iyz+1
              enddo
            enddo
          else
            do n=n1,n2
              do m=m1,m2
                f(l1:l2,m,n,icc) = amplcc*tmpx*ylm(ll_sh,mm_sh)
              enddo
            enddo
          endif

        case default; call fatal_error('init_lncc','no such initcc: '//trim(initcc))
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
    subroutine pencil_criteria_pscalar
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
      if (lpscalar_sink) lpenc_requested(i_cc)=.true.
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
      if (diffcc_shock /= 0.) then
        lpenc_requested(i_gcc) = .true.
        lpenc_requested(i_del2cc) = .true.
        lpenc_requested(i_shock) = .true.
        lpenc_requested(i_gshock) = .true.
      endif
      if (pscalar_diff_hyper3/=0.0) then
        lpenc_requested(i_del6cc)=.true.
        lpenc_requested(i_g5ccglnrho)=.true.
      endif
!
      lpenc_diagnos(i_cc)=.true.
!
      if (idiag_rhoccm/=0 .or. idiag_Cz2m/=0 .or. idiag_Cz4m/=0 .or. &
          idiag_Qrhoccm/=0 .or. idiag_Qpsclm/=0 .or. idiag_mrclncm/=0. ) &
          lpenc_diagnos(i_rho)=.true.
!
      if (idiag_ucm/=0 .or. idiag_uudcm/=0 .or. idiag_uxcm/=0 .or. &
          idiag_uycm/=0 .or. idiag_uzcm/=0 .or. idiag_cluz_uzlcm/=0 ) &
          lpenc_diagnos(i_uu)=.true.
!
      if (idiag_uudcm/=0 .or. idiag_cugccm/=0) lpenc_diagnos(i_ugcc)=.true.
!
      if (idiag_cc1m/=0 .or. idiag_cc2m/=0 .or. idiag_cc3m/=0 .or. &
          idiag_cc4m/=0 .or. idiag_cc5m/=0 .or. idiag_cc6m/=0 .or. &
          idiag_cc7m/=0 .or. idiag_cc8m/=0 .or. idiag_cc9m/=0 .or. &
          idiag_cc10m/=0) lpenc_diagnos(i_cc1)=.true.
!
      if (idiag_gcc1m/=0 .or. idiag_gcc2m/=0 .or. idiag_gcc3m/=0 .or. &
          idiag_gcc4m/=0 .or. idiag_gcc5m/=0 .or. idiag_gcc6m/=0 .or. &
          idiag_gcc7m/=0 .or. idiag_gcc8m/=0 .or. idiag_gcc9m/=0 .or. &
          idiag_gcc10m/=0) lpenc_diagnos(i_gcc1)=.true.
!
      if (idiag_ccglnrm/=0) lpenc_diagnos(i_glnrho)=.true.
      if (idiag_ccugum/=0) lpenc_diagnos(i_ugu)=.true.
      if (idiag_gcguzm/=0) then
       lpenc_diagnos(i_uij)=.true.
       lpenc_diagnos(i_gcc)=.true.
      endif
      if (idiag_cluz_uzlcm/=0) then
        lpenc_diagnos(i_del2u)=.true.
        lpenc_diagnos(i_del2cc)=.true.
      endif
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
!  24-oct-11/ccyang: generalize to multiple scalars
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      real, dimension(nx) :: dot2_tmp
      integer :: i
! cc
      if (lpencil(i_cc)) p%cc=f(l1:l2,m,n,icc:icc+npscalar-1)
! cc1
      if (lpencil(i_cc1)) p%cc1=1./p%cc
! gcc
      if (lpencil(i_gcc)) then
        do i = 1, npscalar
          call grad(f,icc+i-1,p%gcc(:,:,i))
        enddo
      endif
! ugcc
      if (lpencil(i_ugcc)) then
        do i = 1, npscalar
          call u_dot_grad(f,icc+i-1,p%gcc(:,:,i),p%uu,p%ugcc(:,i),UPWIND=lupw_cc)
        enddo
      endif
! gcc2
      if (lpencil(i_gcc2)) then
        do i = 1, npscalar
          call dot2_mn(p%gcc(:,:,i),dot2_tmp)
          p%gcc2(:,i)=dot2_tmp
        enddo
      endif
! gcc1
      if (lpencil(i_gcc1)) p%gcc1=sqrt(p%gcc2)
! del2cc
      if (lpencil(i_del2cc)) then
        do i = 1, npscalar
          call del2(f,icc+i-1,p%del2cc(:,i))
        enddo
      endif
! hcc
      if (lpencil(i_hcc)) then
        do i = 1, npscalar
          call g2ij(f,icc+i-1,p%hcc(:,:,:,i))
        enddo
      endif
! del6cc
      if (lpencil(i_del6cc)) then
        do i = 1, npscalar
          call del6(f,icc+i-1,p%del6cc(:,i))
        enddo
      endif
! g5cc
      if (lpencil(i_g5cc)) then
        do i = 1, npscalar
          call grad5(f,icc+i-1,p%g5cc(:,:,i))
        enddo
      endif
! g5ccglnrho
      if (lpencil(i_g5ccglnrho)) then
        do i = 1, npscalar
          call dot_mn(p%g5cc(:,:,i),p%glnrho,p%g5ccglnrho(:,i))
        enddo
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
      use Special, only: special_calc_pscalar
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: diff_op,diff_op2,diffus_pscalar,diffus_pscalar3
      real :: cc_xyaver
      real :: lam_gradC_fact=1., om_gradC_fact=1., gradC_fact=1.
      integer :: j, k
!
      intent(in)  :: f
      intent(out) :: df
!
      character(len=2) :: id
      integer :: icc2
      real, dimension(nx,3) :: tmp
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) then
        if (nopscalar) then
          print*,'not SOLVED: dlncc_dt'
        else
          print*,'SOLVE dlncc_dt'
        endif
      endif

      if (headtt) then
        do k = 1, npscalar
          write(id,'(i0)') k
          call identify_bcs('cc'//trim(id),icc+k-1)
        enddo
      endif
!
!  Gradient of passive scalar.
!  Allow for possibility to turn off passive scalar
!  without changing file size and recompiling everything.
!
      evolve: if (.not. nopscalar) then ! i.e. if (pscalar)
        icc2=icc+npscalar-1
!
!  Passive scalar equation.
!
        df(l1:l2,m,n,icc:icc2)=df(l1:l2,m,n,icc:icc2)-p%ugcc
!
!  lpscalar_per_unitvolume
!
        if (lpscalar_per_unitvolume) &
          df(l1:l2,m,n,icc:icc2)=df(l1:l2,m,n,icc:icc2)-p%cc*spread(p%divu,2,npscalar)
!
!  Reaction term. Simple Fisher term for now.
!
        if (lreactions) df(l1:l2,m,n,icc:icc2)=df(l1:l2,m,n,icc:icc2)+lambda_cc*p%cc*(1.0-p%cc)
!
!  Passive scalar sink.
!
        if (lpscalar_sink) then
          if (Rpscalar_sink==0) then
            if (ll_sh/=-1.and.mm_sh/=-1) then
              bump=pscalar_sink*spharm(m-m1+1,n-n1+1)
            else
              bump=pscalar_sink
            endif
          else
            bump=pscalar_sink*exp(-0.5*(x(l1:l2)**2+y(m)**2+z(n)**2)/Rpscalar_sink**2)
          endif
          df(l1:l2,m,n,icc:icc2)=df(l1:l2,m,n,icc:icc2)-spread(bump,2,npscalar)*p%cc
        endif
!
!  Diffusion operator. If lpscalar_per_unitvolume is chosen, use
!  div[kappa*rho*grad(c/rho)] = kappa*(del2c-gradc.gradlnrho-c*del2lnrho).
!  Otherwise, with lpscalar_per_unitvolume and not lpscalar_per_unitvolume_diff
!  use just kappa*del2c.
!
        if (pscalar_diff/=0.) then

          if (headtt) print*,'dlncc_dt: pscalar_diff=',pscalar_diff
          do k = 1, npscalar
            if (lpscalar_per_unitvolume) then
              if (lpscalar_per_unitvolume_diff) then
                call dot_mn(p%glnrho,p%gcc(:,:,k),diff_op)
                diff_op=p%del2cc(:,k)-diff_op-p%cc(:,k)*p%del2lnrho
              else
                diff_op=p%del2cc(:,k)
              endif
            else
              if (lpscalar_diff_simple) then
                diff_op=p%del2cc(:,k)
              else
                call dot_mn(p%glnrho,p%gcc(:,:,k),diff_op)
                diff_op=diff_op+p%del2cc(:,k)
              endif
            endif
            df(l1:l2,m,n,icc+k-1)=df(l1:l2,m,n,icc+k-1)+pscalar_diff*diff_op
          enddo

        endif
!
!  Shock Diffusion
!
        if (diffcc_shock /= 0.) then
          do k = 1, npscalar
            call dot_mn(p%gshock, p%gcc(:,:,k), diff_op)
            df(l1:l2,m,n,icc+k-1)=df(l1:l2,m,n,icc+k-1) + diffcc_shock*(p%shock*p%del2cc(:,k) + diff_op)
          enddo
        endif
!
!  Hyperdiffusion operator.
!
        if (pscalar_diff_hyper3/=0.) then
          if (headtt) print*,'dlncc_dt: pscalar_diff_hyper3=', pscalar_diff_hyper3
          df(l1:l2,m,n,icc:icc2)=df(l1:l2,m,n,icc:icc2)+pscalar_diff_hyper3*(p%del6cc+p%g5ccglnrho)
        endif
!
!  Soret diffusion.
!
        if (soret_diff/=0.) then
          if (headtt) print*,'dlncc_dt: soret_diff=',soret_diff
          call dot2_mn(p%glnTT,diff_op)
          diff_op=p%TT*(diff_op+p%del2lnTT)
          do k = 1, npscalar
            diff_op2=p%cc(:,k)*(1.-p%cc(:,k))*diff_op
            df(l1:l2,m,n,icc+k-1)=df(l1:l2,m,n,icc+k-1)+soret_diff*diff_op2
          enddo
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
        if (lgradC_profile) gradC_fact=gradC_fact*cos(z(n))
!
!  Add diffusion of imposed spatially constant gradient of c.
!  This makes sense really only for periodic boundary conditions.
!
        do j=1,3
          if (gradC0(j)/=0.) df(l1:l2,m,n,icc:icc2)=df(l1:l2,m,n,icc:icc2) &
                             -spread(gradC0(j)*p%uu(:,j)*gradC_fact,2,npscalar)
        enddo
!
!  Tensor diffusion (but keep the isotropic one).
!
        if (tensor_pscalar_diff/=0.) call tensor_diff(df,p,tensor_pscalar_diff)
!
!  Consider here the action of a mean friction term, -LLambda*Cbar.
!  This can be used to compendate for the decay of a horizontally
!  averaged mean concentration and allows thus the determination of
!  the turbulent diffusivity under stationary conditions. Only those
!  results are then comparable with the results of the test-field method.
!
        if (lmean_friction_cc) then
          do k = icc, icc2
            cc_xyaver=sum(f(l1:l2,m1:m2,n,k))/nxygrid   !only for nprocxy=1 - tb improved: calc cc_xyaver in before_boundary
            df(l1:l2,m,n,k)=df(l1:l2,m,n,k)-LLambda_cc*cc_xyaver
          enddo
        endif
!
!  For the timestep calculation, need maximum diffusion.
!
        if (lfirst.and.ldt) then
          diffus_pscalar =(pscalar_diff+tensor_pscalar_diff)*dxyz_2
          diffus_pscalar3=pscalar_diff_hyper3*dxyz_6
          maxdiffus=max(maxdiffus,diffus_pscalar)
          maxdiffus3=max(maxdiffus3,diffus_pscalar3)
        endif
!
!  Special contributions to this module are called here.
!
        if (lspecial) call special_calc_pscalar(f,df,p)
!
! AH: notpassive, an angular momentum+gravity workaround
!
        if (lnotpassive.and.lhydro) then
          df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux)+(p%cc(:,1)+(-powerlr*hoverr**2+1.5*zoverh**2*hoverr**2)+ &
                              (-1+3*powerlr*hoverr**2-4.5*zoverh**2*hoverr**2)*x(l1:l2))*scalaracc
          df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz)+(-hoverr*zoverh - z(n)+3*hoverr*zoverh*x(l1:l2))*scalaracc
        endif
!
        call calc_diagnostics_pscalar(p)

      endif evolve

    endsubroutine dlncc_dt
!***********************************************************************
    subroutine calc_diagnostics_pscalar(p)
!
!  Diagnostics (only for the first passive scalar)
!
!  output for double and triple correlators (assume z-gradient of cc)
!  <u_k u_j d_j c> = <u_k c uu.gradcc>
!
      use Diagnostics
      use Sub, only: dot_mn

      type(pencil_case), intent(IN) :: p

      real, dimension (nx) :: gcgu

      if (ldiagnos) then
        call sum_mn_name(bump,idiag_Qpsclm)
        if (idiag_Qrhoccm/=0) call sum_mn_name(bump*p%rho*p%cc(:,1),idiag_Qrhoccm)
        if (idiag_mcct/=0)    call integrate_mn_name(p%rho*p%cc(:,1),idiag_mcct)
        if (idiag_rhoccm/=0)  call sum_mn_name(p%rho*p%cc(:,1),idiag_rhoccm)
        if (idiag_mrclncm/=0) call sum_mn_name(-p%rho*p%cc(:,1)*alog(p%cc(:,1)),idiag_mrclncm)
        call max_mn_name(p%cc(:,1),idiag_ccmax)
        if (idiag_ccmin/=0)   call max_mn_name(-p%cc(:,1),idiag_ccmin,lneg=.true.)
        if (idiag_uxcm/=0)    call sum_mn_name(p%uu(:,1)*p%cc(:,1),idiag_uxcm)
        if (idiag_uycm/=0)    call sum_mn_name(p%uu(:,2)*p%cc(:,1),idiag_uycm)
        if (idiag_uzcm/=0)    call sum_mn_name(p%uu(:,3)*p%cc(:,1),idiag_uzcm)
        if (lgradC_profile) then
          if (idiag_ucm/=0)   call sum_mn_name(2.*cos(z(n))*p%uu(:,3)*p%cc(:,1),idiag_ucm)
        else
          if (idiag_ucm/=0)   call sum_mn_name(p%uu(:,3)*p%cc(:,1),idiag_ucm)
        endif
        if (idiag_uudcm/=0)   call sum_mn_name(p%uu(:,3)*p%ugcc(:,1),idiag_uudcm)
        if (idiag_Cz2m/=0)    call sum_mn_name(p%rho*p%cc(:,1)*z(n)**2,idiag_Cz2m)
        if (idiag_Cz4m/=0)    call sum_mn_name(p%rho*p%cc(:,1)*z(n)**4,idiag_Cz4m)
        if (idiag_Crmsm/=0)   call sum_mn_name((p%rho*p%cc(:,1))**2,idiag_Crmsm,lsqrt=.true.)
        if (idiag_ccrms/=0)   call sum_mn_name(p%cc(:,1)**2,idiag_ccrms,lsqrt=.true.)
        call sum_mn_name(p%cc1(:,1),idiag_cc1m)
        if (idiag_cc2m/=0)    call sum_mn_name(p%cc1(:,1)**2,idiag_cc2m)
        if (idiag_cc3m/=0)    call sum_mn_name(p%cc1(:,1)**3,idiag_cc3m)
        if (idiag_cc4m/=0)    call sum_mn_name(p%cc1(:,1)**4,idiag_cc4m)
        if (idiag_cc5m/=0)    call sum_mn_name(p%cc1(:,1)**5,idiag_cc5m)
        if (idiag_cc6m/=0)    call sum_mn_name(p%cc1(:,1)**6,idiag_cc6m)
        if (idiag_cc7m/=0)    call sum_mn_name(p%cc1(:,1)**7,idiag_cc7m)
        if (idiag_cc8m/=0)    call sum_mn_name(p%cc1(:,1)**8,idiag_cc8m)
        if (idiag_cc9m/=0)    call sum_mn_name(p%cc1(:,1)**9,idiag_cc9m)
        if (idiag_cc10m/=0)   call sum_mn_name(p%cc1(:,1)**10,idiag_cc10m)
        call sum_mn_name(p%gcc1(:,1),idiag_gcc1m)
        if (idiag_gcc2m/=0)   call sum_mn_name(p%gcc1(:,1)**2,idiag_gcc2m)
        if (idiag_gcc3m/=0)   call sum_mn_name(p%gcc1(:,1)**3,idiag_gcc3m)
        if (idiag_gcc4m/=0)   call sum_mn_name(p%gcc1(:,1)**4,idiag_gcc4m)
        if (idiag_gcc5m/=0)   call sum_mn_name(p%gcc1(:,1)**5,idiag_gcc5m)
        if (idiag_gcc6m/=0)   call sum_mn_name(p%gcc1(:,1)**6,idiag_gcc6m)
        if (idiag_gcc7m/=0)   call sum_mn_name(p%gcc1(:,1)**7,idiag_gcc7m)
        if (idiag_gcc8m/=0)   call sum_mn_name(p%gcc1(:,1)**8,idiag_gcc8m)
        if (idiag_gcc9m/=0)   call sum_mn_name(p%gcc1(:,1)**9,idiag_gcc9m)
        if (idiag_gcc10m/=0)  call sum_mn_name(p%gcc1(:,1)**10,idiag_gcc10m)
        if (idiag_ccglnrm/=0) call sum_mn_name(p%cc(:,1)*p%glnrho(:,3),idiag_ccglnrm)
        if (idiag_cugccm/=0)  call sum_mn_name(p%cc(:,1)*p%ugcc(:,1),idiag_cugccm)
        if (idiag_ccugum/=0)  call sum_mn_name(p%cc(:,1)*p%ugu(:,3),idiag_ccugum)
        if (idiag_gcguzm/=0)  then
          call dot_mn(p%gcc(:,:,1),p%uij(:,3,:),gcgu)
          call sum_mn_name(gcgu,idiag_gcguzm)
        endif
        if (idiag_cluz_uzlcm/=0) &
          call sum_mn_name(p%cc(:,1)*p%del2u(:,3)-p%uu(:,3)*p%del2cc(:,1),idiag_cluz_uzlcm)
!
      endif
!
      if (l1davgfirst) then
        call xysum_mn_name_z(p%cc(:,1),idiag_ccmz)
        if (idiag_cc2mz/=0) call xysum_mn_name_z(p%cc(:,1)**2,idiag_cc2mz)
        call xzsum_mn_name_y(p%cc(:,1),idiag_ccmy)
        call yzsum_mn_name_x(p%cc(:,1),idiag_ccmx)
        if (idiag_uxcmz/=0) call xysum_mn_name_z(p%uu(:,1)*p%cc(:,1),idiag_uxcmz)
        if (idiag_uycmz/=0) call xysum_mn_name_z(p%uu(:,2)*p%cc(:,1),idiag_uycmz)
        if (idiag_uzcmz/=0) call xysum_mn_name_z(p%uu(:,3)*p%cc(:,1),idiag_uzcmz)
      endif
!
      if (l2davgfirst) then
        call zsum_mn_name_xy(p%cc(:,1),idiag_ccmxy)
        call ysum_mn_name_xz(p%cc(:,1),idiag_ccmxz)
      endif

    endsubroutine calc_diagnostics_pscalar
!***********************************************************************
    subroutine read_pscalar_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=pscalar_init_pars, IOSTAT=iostat)
!
    endsubroutine read_pscalar_init_pars
!***********************************************************************
    subroutine write_pscalar_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=pscalar_init_pars)
!
    endsubroutine write_pscalar_init_pars
!***********************************************************************
    subroutine read_pscalar_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=pscalar_run_pars, IOSTAT=iostat)
!
    endsubroutine read_pscalar_run_pars
!***********************************************************************
    subroutine write_pscalar_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=pscalar_run_pars)
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
      use FArrayManager, only: farray_index_append
!
      logical :: lreset
      logical, optional :: lwrite
!
      character(len=80) :: fmt
      integer :: iname, inamez, inamey, inamex, inamexy, inamexz
      integer :: i
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
        idiag_mrclncm=0
        idiag_Qrhoccm=0; idiag_Qpsclm=0; idiag_mcct=0
        idiag_ccmz=0; idiag_ccmy=0; idiag_ccmx=0
        idiag_uxcmz=0; idiag_uycmz=0; idiag_uzcmz=0; idiag_cc2mz=0
        idiag_ucm=0; idiag_uudcm=0; idiag_Cz2m=0; idiag_Cz4m=0; idiag_Crmsm=0
        idiag_ccrms=0
        idiag_uxcm=0; idiag_uycm=0; idiag_uzcm=0
        idiag_cc1m=0; idiag_cc2m=0; idiag_cc3m=0; idiag_cc4m=0; idiag_cc5m=0
        idiag_cc6m=0; idiag_cc7m=0; idiag_cc8m=0; idiag_cc9m=0; idiag_cc10m=0
        idiag_gcc1m=0; idiag_gcc2m=0; idiag_gcc3m=0; idiag_gcc4m=0
        idiag_gcc5m=0; idiag_gcc6m=0; idiag_gcc7m=0; idiag_gcc8m=0
        idiag_gcc9m=0; idiag_gcc10m=0; idiag_ccglnrm=0; idiag_cugccm=0
        idiag_ccugum=0; idiag_cluz_uzlcm=0; idiag_gcguzm=0
        idiag_ccmxy=0; idiag_ccmxz=0
      endif
!
!  Check for those quantities that we want to evaluate online.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Qpsclm',idiag_Qpsclm)
        call parse_name(iname,cname(iname),cform(iname),'Qrhoccm',idiag_Qrhoccm)
        call parse_name(iname,cname(iname),cform(iname),'rhoccm',idiag_rhoccm)
        call parse_name(iname,cname(iname),cform(iname),'mrclncm',idiag_mrclncm)
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
        call parse_name(iname,cname(iname),cform(iname),'ccrms',idiag_ccrms)
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
        call parse_name(iname,cname(iname),cform(iname),'cugccm',idiag_cugccm)
        call parse_name(iname,cname(iname),cform(iname),'ccugum',idiag_ccugum)
        call parse_name(iname,cname(iname),cform(iname),'gcguzm',idiag_gcguzm)
        call parse_name(iname,cname(iname),cform(iname),'cluz_uzlcm',idiag_cluz_uzlcm)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uxcmz',idiag_uxcmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uycmz',idiag_uycmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uzcmz',idiag_uzcmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ccmz',idiag_ccmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'cc2mz',idiag_cc2mz)
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
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then 
        where(cnamev=='cc'.or.cnamev=='lncc') cformv='DEFINED'
      endif
!
      if (lwr) then
        call farray_index_append('npscalar', npscalar)
      endif
!
    endsubroutine rprint_pscalar
!***********************************************************************
    subroutine get_slices_pscalar(f,slices)
!
!  Write slices for animation of pscalar variables.
!
!  26-jul-06/tony: coded
!  31-jan-11/ccyang: generalized to multiple scalars
!
      use Slices_methods, only: assign_slices_vec, process_slices, log2d
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type(slice_data), intent(inout) :: slices
      character(LEN=fmtlen) :: sname
!
!  Loop over slices.
!
      sname=trim(slices%name)
!
!  Passive scalar or logarithm of it..
!
      if (sname=='cc'.or.sname=='lncc') then

        call assign_slices_vec(slices,f,icc,npscalar)
!
!  Logarithmic passive scalar.
!
        if (sname=='lncc') call process_slices(slices,log2d)       
!
      endif
!
    endsubroutine get_slices_pscalar
!***********************************************************************
    subroutine pscalar_before_boundary(f) 
!
!  Removes overall means of passive scalars.
!
!  5-dec-11/MR: coded
!
      use Mpicomm, only: mpiallreduce_sum

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f

      real, dimension(npscalar) :: ccm, ccm_tmp
      integer :: i

      if (lremove_mean.and.lrmv) then

        do i=1,npscalar
          ccm_tmp(i) = sum(f(l1:l2,m1:m2,n1:n2,icc+i-1))
        enddo

        call mpiallreduce_sum(ccm_tmp,ccm,npscalar)

        do i=1,npscalar
          f(l1:l2,m1:m2,n1:n2,icc+i-1)=f(l1:l2,m1:m2,n1:n2,icc+i-1)-ccm(i)/nwgrid
        enddo

      endif

    endsubroutine pscalar_before_boundary
!***********************************************************************
    subroutine calc_mpscalar
!
!  Calculate mean magnetic field from xy- or z-averages.
!
!  14-apr-03/axel: adaped from calc_mfield
!
      use Diagnostics
!
      real :: ccm
!
!  RMS value of horizontally averaged concentration of 1st pscalar.
!  ccmz must have been calculated, so it is present on the root processor.
!
      if (idiag_ccm/=0) then
        if (idiag_ccmz==0) then
          if (headt) call warning('calc_mpscalar',"to get ccm, ccmz must also be set in xyaver.in."// &
                                  achar(10)//"                     We proceed, but you'll get ccm=0")
          ccm=0.
        else
          ccm=sqrt(sum(fnamez(:,:,idiag_ccmz)**2)/nzgrid)
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
!  24-oct-11/ccyang: generalized to multiple scalars
!
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real :: tensor_pscalar_diff
!
      real, dimension (nx) :: tmp,scr
      integer :: iy,iz,i,j,k
!
!  tmp = (Bunit.G)^2 + H.G + Bi*Bj*Gij
!  for details, see tex/mhd/thcond/tensor_der.tex
!
      iy=m-m1+1
      iz=n-n1+1
      do k = 1, npscalar
        call dot_mn(bunit(:,iy,iz,:),p%gcc(:,:,k),scr)
        call dot_mn(hhh(:,iy,iz,:),p%gcc(:,:,k),tmp)
        tmp=tmp+scr**2
!
!  dot with bi*bj
!
        do j=1,3
        do i=1,3
          tmp=tmp+bunit(:,iy,iz,i)*bunit(:,iy,iz,j)*p%hcc(:,i,j,k)
        enddo
        enddo
!
!  and add result to the dcc/dt equation
!
        df(l1:l2,m,n,icc+k-1)=df(l1:l2,m,n,icc+k-1)+tensor_pscalar_diff*tmp
      enddo
!
    endsubroutine tensor_diff
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(pscalar_diff,p_par(1))

    endsubroutine pushpars2c
!***********************************************************************
endmodule Pscalar
