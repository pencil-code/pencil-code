! $Id$
!
!  This module supplies a kinematic velocity field.
!  Most of the content of this module was moved with revision r12019
!  by Dhrubaditya Mitra on 5-nov-09 away from nohydro.f90
!  To inspect the revision history of the file before that time,
!  check out nohydro.f90 prior to or at revision r12018.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module.
!
! CPARAM logical, parameter :: lhydro = .false.
! CPARAM logical, parameter :: lhydro_kinematic = .true.
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED oo(3); o2; ou; uij(3,3); uu(3); u2; sij(3,3)
! PENCILS PROVIDED der6u(3)
! PENCILS PROVIDED divu; ugu(3); del2u(3); uij5(3,3); graddivu(3)
! PENCILS PROVIDED uu_advec(3); uuadvec_guu(3)
!***********************************************************************
module Hydro
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'record_types.h'
  include 'hydro.h'
!
  real, dimension (mx,3) :: uumx=0.
  real, dimension (mz,3) :: uumz=0.
  real, dimension (nz,3) :: guumz=0.
  real, dimension (mx,my,3) :: uumxy=0.
  real, dimension (mx,mz,3) :: uumxz=0.
!
  real, dimension(nx) :: profx_kinflow1=1., profx_kinflow2=1., profx_kinflow3=1.
  real, dimension(my) :: profy_kinflow1=1., profy_kinflow2=1., profy_kinflow3=1.
  real, dimension(mz) :: profz_kinflow1=1., profz_kinflow2=1., profz_kinflow3=1.
!
  real :: u_out_kep=0.0
  real :: tphase_kinflow=-1., phase1=0., phase2=0., tsforce=0.
  real :: tsforce_ampl=0., tsforce_wavenumber=0.
  real ::  dtforce=impossible, ampl_random
  real, dimension(3) :: location,location_fixed=(/0.,0.,0./)
  logical :: lupw_uu=.false., lkinflow_as_uudat=.false.
  logical :: lcalc_uumeanz=.false.,lcalc_uumeanxy=.false.
  logical :: lcalc_uumeanx=.false.,lcalc_uumeanxz=.false.
!
  real, allocatable, dimension (:,:) :: KS_k,KS_A,KS_B !or through whole field for each wavenumber?
  real, allocatable, dimension (:) :: KS_omega !or through whole field for each wavenumber?
  integer :: KS_modes = 25
  real, allocatable, dimension (:) :: Zl,dZldr,Pl,dPldtheta
  real :: ampl_fcont_uu=1., random_ampl, random_wavenumber
  logical :: lforcing_cont_uu=.false., lrandom_location=.false.
  logical :: lwrite_random_location=.false., lwrite_random_wavenumber=.false.
  logical :: lrandom_wavenumber=.false., lwrite_random_ampl=.false.
  real, dimension(nx) :: ck_r,ck_rsqr
  integer :: ll_sh=0, mm_sh=0, n_xprof=-1
  integer :: pushpars2c, pushdiags2c  ! should be procedure pointer (F2003)
!
!  init parameters
!  (none)
!
!  run parameters
!
  character (len=labellen) :: kinematic_flow='none'
  real :: ABC_A=1.0, ABC_B=1.0, ABC_C=1.0
  real :: wind_amp=0.,wind_rmin=impossible,wind_step_width=0.
  real :: circ_amp=0.,circ_rmax=0.,circ_step_width=0.
  real :: kx_uukin=1., ky_uukin=1., kz_uukin=1.
  real :: cx_uukin=0., cy_uukin=0., cz_uukin=0.
  real :: phasex_uukin=0., phasey_uukin=0., phasez_uukin=0.
  real :: radial_shear=0.,uphi_at_rzero=0.,uphi_at_rmax=0.,uphi_rmax=1.
  real :: uphi_rbot=1., uphi_rtop=1., uphi_step_width=0.
  real :: gcs_rzero=0.,gcs_psizero=0.
  real :: kinflow_ck_Balpha=0.
  real :: eps_kinflow=0., exp_kinflow=1., omega_kinflow=0., ampl_kinflow=1.
  real :: rp,gamma_dg11=0.4, relhel_uukin=1., chi_uukin=45., del_uukin=0.
  real :: lambda_kinflow=1., zinfty_kinflow=0.
  integer :: kinflow_ck_ell=0, tree_lmax=8, kappa_kinflow=100
  character (len=labellen) :: wind_profile='none'
  logical, target :: lpressuregradient_gas=.false.
  logical :: lkinflow_as_comaux=.false.
  logical :: lrandom_ampl=.false.
!
  namelist /hydro_run_pars/ &
      kinematic_flow,wind_amp,wind_profile,wind_rmin,wind_step_width, &
      circ_rmax,circ_step_width,circ_amp, ABC_A,ABC_B,ABC_C, &
      ampl_kinflow, relhel_uukin, chi_uukin, del_uukin, &
      kx_uukin,ky_uukin,kz_uukin, &
      cx_uukin,cy_uukin,cz_uukin, &
      phasex_uukin, phasey_uukin, phasez_uukin, &
      lrandom_location, lwrite_random_location, &
      lrandom_wavenumber, lwrite_random_wavenumber, &
      location_fixed, dtforce, &
      lwrite_random_ampl, lkinflow_as_comaux, lkinflow_as_uudat, &
      radial_shear, uphi_at_rzero, uphi_rmax, uphi_rbot, uphi_rtop, &
      uphi_step_width,gcs_rzero, &
      gcs_psizero,kinflow_ck_Balpha,kinflow_ck_ell, &
      eps_kinflow,exp_kinflow,omega_kinflow,ampl_kinflow, rp, gamma_dg11, &
      lambda_kinflow, tree_lmax, zinfty_kinflow, kappa_kinflow, &
      ll_sh, mm_sh, n_xprof, lrandom_ampl
!
  integer :: idiag_u2m=0,idiag_um2=0,idiag_oum=0,idiag_o2m=0
  integer :: idiag_uxpt=0,idiag_uypt=0,idiag_uzpt=0
  integer :: idiag_dtu=0,idiag_urms=0,idiag_umax=0,idiag_uzrms=0
  integer :: idiag_uzmax=0,idiag_orms=0,idiag_omax=0
  integer :: idiag_ux2m=0,idiag_uy2m=0,idiag_uz2m=0
  integer :: idiag_uxuym=0,idiag_uxuzm=0,idiag_uyuzm=0,idiag_oumphi=0
  integer :: idiag_ruxm=0,idiag_ruym=0,idiag_ruzm=0,idiag_rumax=0
  integer :: idiag_uxmz=0,idiag_uymz=0,idiag_uzmz=0,idiag_umx=0
  integer :: idiag_umy=0,idiag_umz=0,idiag_uxmxy=0,idiag_uymxy=0,idiag_uzmxy=0
  integer :: idiag_Marms=0,idiag_Mamax=0,idiag_divu2m=0,idiag_epsK=0
  integer :: idiag_urmphi=0,idiag_upmphi=0,idiag_uzmphi=0,idiag_u2mphi=0
  integer :: idiag_phase1=0,idiag_phase2=0
  integer :: idiag_ekintot=0,idiag_ekin=0
  integer :: idiag_divum=0
!
  logical :: lupdate_aux=.false.
  contains
!***********************************************************************
    subroutine register_hydro
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!   6-nov-01/wolf: coded
!
      use Mpicomm, only: lroot
      use SharedVariables, only: put_shared_variable
!
      integer :: ierr
!
!  Identify version number (generated automatically by SVN).
!
      !if (lroot) call svn_id( &
      !    "$Id$")
!
      call put_shared_variable('lpressuregradient_gas',&
          lpressuregradient_gas,caller='register_hydro')
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine initialize_hydro(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!  12-sep-13/MR  : calculation of means added
!  21-sep-13/MR  : adjusted use of calc_pencils_hydro
!  10-jun-17/MR  : added poloidal velocity from a single spherical harmonic,
!                  also oscillating with omega_kinflow
!
      use FArrayManager
      use Sub, only: erfunc, ylm, ylm_other
      use Boundcond, only: update_ghosts
      use General
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: exp_kinflow1, sph, sph_har_der
      real, dimension (nx) :: vel_prof, tmp_mn
      real, dimension (nx,3) :: tmp_nx3
      real, dimension (:,:), allocatable :: yz
      integer :: iyz
!
!  Compute preparatory functions needed to assemble
!  different flow profiles later on in pencil_case.
!
      select case (kinematic_flow)
!
!  Spoke-like differential rotation profile.
!  The minus sign is needed for equatorward acceleration.
!
      case ('Brandt')
        if (lcylindrical_coords) then
          exp_kinflow1=1./exp_kinflow
          profx_kinflow1=+1./(1.+(x(l1:l2)/uphi_rbot)**exp_kinflow)**exp_kinflow1
          profy_kinflow1=-1.
          profz_kinflow1=+1.
        endif
      case ('spoke-like')
        profx_kinflow1=+0.5*(1.+erfunc(((x(l1:l2)-uphi_rbot)/uphi_step_width)))
        profy_kinflow1=-1.5*(5.*cos(y)**2-1.)
        profz_kinflow1=+1.
      case ('spoke-like-NSSL')
        profx_kinflow1=+0.5*(1.+erfunc(((x(l1:l2)-uphi_rbot)/uphi_step_width)))
        profx_kinflow2=+0.5*(1.-erfunc(((x(l1:l2)-uphi_rtop)/uphi_step_width)))
        profx_kinflow3=+0.5*(1.+erfunc(((x(l1:l2)-uphi_rtop)/uphi_step_width)))
        profx_kinflow2=(x(l1:l2)-uphi_rbot)*profx_kinflow1*profx_kinflow2  !(redefined)
        profy_kinflow1=-1.5*(5.*cos(y)**2-1.)
        profy_kinflow2=-1.0*(4.*cos(y)**2-3.)
        profy_kinflow3=-1.0
        profz_kinflow1=+1.
      case ('KS')
        call periodic_KS_setup(-5./3.) !Kolmogorov spec. periodic KS
        !call random_isotropic_KS_setup(-5./3.,1.,(nxgrid)/2.) !old form
        !call random_isotropic_KS_setup_test !Test KS model code with 3 specific modes.
      case ('ck')
        call init_ck
      case default;
        if (lroot .and. (ip < 14)) call information('initialize_hydro','no preparatory profile needed')
      end select
!
! kinflows end here
!
!  Register an extra aux slot for uu if requested (so uu is written
!  to snapshots and can be easily analyzed later). For this to work you
!  must reserve enough auxiliary workspace by setting, for example,
!     ! MAUX CONTRIBUTION 3
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
!
!  After a reload, we need to rewrite index.pro, but the auxiliary
!  arrays are already allocated and must not be allocated again.
!
      if (lkinflow_as_aux.or.lkinflow_as_comaux) then
        if (iuu==0) then
          if (lkinflow_as_comaux) then
            call farray_register_auxiliary('uu',iuu,vector=3,communicated=.true.)
          else
            call farray_register_auxiliary('uu',iuu,vector=3)
          endif
          iux=iuu
          iuy=iuu+1
          iuz=iuu+2
!
!  Possibility to read uu.dat, but currently only for one processor.
!  However, this can be useful when a periodic kinematic flow is to
!  to be used in test-field analyses at subharmic wavenumbers,
!  because then each processor uses the same flow on each pressor.
!  Note, however, that lkinflow_as_uudat=.false. by default.
!
          if (lkinflow_as_uudat) then
            open(1,file='uu.dat',form='unformatted')
            read(1) f(l1:l2,m1:m2,n1:n2,iux:iuz)
            close(1)
            if (ampl_kinflow/=1.) f(l1:l2,m1:m2,n1:n2,iux:iuz) &
                    =ampl_kinflow*f(l1:l2,m1:m2,n1:n2,iux:iuz)
            if (lkinflow_as_comaux) call update_ghosts(f,iux,iuz)
          endif
        else
! set the initial velocity to zero
          if (kinematic_flow/='from-snap') f(:,:,:,iux:iuz) = 0.
          if (lroot .and. (ip<14)) print*, 'initialize_hydro: iuu = ', iuu
          call farray_index_append('iuu',iuu,3)
        endif

        if (kinematic_flow=='spher-harm-poloidal'.or.kinematic_flow=='spher-harm-poloidal-per') then
          if (.not.lspherical_coords) call inevitably_fatal_error("init_uu", &
              '"spher-harm-poloidal" only meaningful for spherical coordinates')
          if (.not.lkinflow_as_aux) call inevitably_fatal_error("init_uu", &
              '"spher-harm-poloidal" requires lkinflow_as_aux=T')
          if (n_xprof==-1) then
            tmp_mn=(x(l1:l2)-xyz0(1))*(x(l1:l2)-xyz1(1))
            vel_prof=tmp_mn/x(l1:l2) + 2.*x(l1:l2)-(xyz0(1)+xyz1(1))     ! f/r + d f/dr 
          else
            tmp_mn=sin((2.*pi/(Lxyz(1))*n_xprof)*(x(l1:l2)-xyz0(1)))
            vel_prof=tmp_mn/x(l1:l2) + (2.*pi/(Lxyz(1))*n_xprof)*cos((2.*pi/(Lxyz(1))*n_xprof)*(x(l1:l2)-xyz0(1)))
          endif
          if (lyang) then
            allocate(yz(2,ny*nz))
            call yin2yang_coors(costh(m1:m2),sinth(m1:m2),cosph(n1:n2),sinph(n1:n2),yz)
            iyz=1
            do m=m1,m2
              do n=n1,n2
                sph=ampl_kinflow*ylm_other(yz(1,iyz),yz(2,iyz),ll_sh,mm_sh,sph_har_der)
                tmp_nx3(:,1)=sph*2.*tmp_mn
                tmp_nx3(:,2)=ampl_kinflow*vel_prof*sph_har_der
                if (mm_sh/=0) then
                  tmp_nx3(:,3) = -sph*vel_prof*mm_sh/sin(yz(1,iyz))*sin(mm_sh*yz(2,iyz))/cos(mm_sh*yz(2,iyz))
                else
                  tmp_nx3(:,3) = 0.
                endif
                call transform_thph_yy( tmp_nx3, (/1,1,1/), f(l1:l2,m,n,iux:iuz), yz(1,iyz), yz(2,iyz) )
                iyz=iyz+1
              enddo
            enddo
          else
            do n=n1,n2
              do m=m1,m2
                sph=ampl_kinflow*ylm(ll_sh,mm_sh,sph_har_der)
                f(l1:l2,m,n,iux) = 2.*tmp_mn*sph
                f(l1:l2,m,n,iuy) = ampl_kinflow*vel_prof*sph_har_der
                if (mm_sh/=0) then
                  f(l1:l2,m,n,iuz) = -vel_prof*sph*mm_sh/sinth(m)*sin(mm_sh*z(n))/cos(mm_sh*z(n))    ! d/d phi / sin(theta)
                else
                  f(l1:l2,m,n,iuz) = 0.
                endif
              enddo
            enddo
          endif
        endif
      endif
!
      call calc_means_hydro(f)
!
    endsubroutine initialize_hydro
!***********************************************************************
    subroutine calc_means_hydro(f)
!
! calculates various means
!
! 14-oct-13/MR: outsourced from initialize_hydro
!
      use Sub, only: finalize_aver

      real, dimension (mx,my,mz,mfarray), intent(IN) :: f
!
      type(pencil_case),dimension(:), allocatable :: p          ! vector as scalar quantities not allocatable
      logical, dimension(:), allocatable :: lpenc_loc
!
      real :: facxy, facyz, facy, facz
      logical :: headtt_save
!
      if (lcalc_uumeanz.or.lcalc_uumeanx.or.lcalc_uumeanxy.or.lcalc_uumeanxz) then

        allocate(p(1),lpenc_loc(npencils))

        facz  = 1./nzgrid
        facy  = 1./nygrid
        facxy = 1./nxygrid
        facyz = 1./nyzgrid

        lpenc_loc = .false.; lpenc_loc(i_uu)=.true.
!
        headtt_save=headtt
        do n=1,mz; do m=1,my
!
          call calc_pencils_hydro(f,p(1),lpenc_loc)
!
          if (lcalc_uumeanz .and. m>=m1 .and. m<=m2) &
            uumz(n,:) = uumz(n,:) + facxy*sum(p(1)%uu,1)
          if (lcalc_uumeanx .and. m>=m1 .and. m<=m2 .and. n>=n1 .and. n<=n2) &
            uumx(l1:l2,:) = uumx(l1:l2,:) + facyz*p(1)%uu
          if (lcalc_uumeanxy .and. n>=n1 .and. n<=n2) &
            uumxy(l1:l2,m,:) = uumxy(l1:l2,m,:) + facz*p(1)%uu
          if (lcalc_uumeanxz .and. m>=m1 .and. m<=m2) &
            uumxz(l1:l2,n,:) = uumxz(l1:l2,n,:) + facy*p(1)%uu
!
          headtt=.false.

        enddo; enddo
!
        headtt=headtt_save

        if (lcalc_uumeanz ) &
          call finalize_aver(nprocxy,12,uumz)
        if (lcalc_uumeanx ) &
          call finalize_aver(nprocyz,23,uumx)
        if (lcalc_uumeanxy) &
          call finalize_aver(nprocz,3,uumxy)
        if (lcalc_uumeanxz) &
          call finalize_aver(nprocy,2,uumxz)
!
      endif
!
    endsubroutine calc_means_hydro
!***********************************************************************
    subroutine init_uu(f)
!
!  Initialise uu.
!
!   7-jun-02/axel: adapted from hydro
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_uu
!***********************************************************************
    subroutine pencil_criteria_hydro
!
!  All pencils that the Hydro module depends on are specified here.
!
!  20-nov-04/anders: coded
!   1-jul-09/axel: added more for kinflow
!
!  pencils for kinflow
!
!AB: i_ugu is not normally required
!--   lpenc_requested(i_ugu)=.true.
!
!DM: The following line with kinflow can be later removed and the variable
!DM: kinematic_flow replaced by kinflow.
!
      kinflow=kinematic_flow
      if (kinflow/='') then
        lpenc_requested(i_uu)=.true.
        if (kinflow=='eddy') then
          lpenc_requested(i_rcyl_mn)=.true.
          lpenc_requested(i_rcyl_mn1)=.true.
        endif
      endif
!
!  Diagnostic pencils.
!
      if (idiag_urms/=0 .or. idiag_umax/=0 .or. idiag_u2m/=0 .or. &
          idiag_um2/=0) lpenc_diagnos(i_u2)=.true.
      if (idiag_orms/=0 .or. idiag_omax/=0 .or. idiag_o2m/=0) &
          lpenc_diagnos(i_o2)=.true.
      if (idiag_oum/=0) lpenc_diagnos(i_ou)=.true.
      if (idiag_divum/=0) lpenc_diagnos(i_divu)=.true.
!
      if (idiag_ekin/=0 .or. idiag_ekintot/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_u2)=.true.
      endif

    endsubroutine pencil_criteria_hydro
!***********************************************************************
    subroutine pencil_interdep_hydro(lpencil_in)
!
!  Interdependency among pencils from the Hydro module is specified here
!
!  20-nov-04/anders: coded
!
      logical, dimension (npencils) :: lpencil_in
!
      if (lpencil_in(i_uglnrho)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_ugrho)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_grho)=.true.
      endif
      if (lpencil_in(i_uij5glnrho)) then
        lpencil_in(i_uij5)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_u2)) lpencil_in(i_uu)=.true.
! oo
      if (lpencil_in(i_ou)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_oo)=.true.
      endif
      if (lpencil_in(i_oo)) lpencil_in(i_uij)=.true.
      if (lpencil_in(i_divu)) lpencil_in(i_uij)=.true.
!  ugu
      if (lpencil_in(i_ugu)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_uij)=.true.
      endif
!
    endsubroutine pencil_interdep_hydro
!***********************************************************************
    subroutine calc_pencils_hydro_std(f,p)
!
! Envelope adjusting calc_pencils_hydro_pencpar to the standard use with
! lpenc_loc=lpencil
!
! 21-sep-13/MR    : coded
!
      real, dimension (mx,my,mz,mfarray),intent(IN) :: f
      type (pencil_case),                intent(OUT):: p
!
      call calc_pencils_hydro_pencpar(f,p,lpencil)
!
      endsubroutine calc_pencils_hydro_std
!***********************************************************************
    subroutine calc_pencils_hydro_pencpar(f,p,lpenc_loc)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   08-nov-04/tony: coded
!   12-sep-13/MR  : optional parameter lpenc_loc added for possibility
!                   to calculate less pencils than in the global setting
!   20-sep-13/MR  : lpenc changed into list of indices in lpencil, penc_inds,
!                   for which pencils are calculated, default: all
!   21-sep-13/MR  : returned to pencil mask as parameter lpenc_loc
!   20-oct-13/MR  : Glen Roberts flow w.r.t. x and z added
!   06-dec-13/MR  : error message for eps_kinflow=0 in Roberts flow IV added
!   15-sep-14/MR  : div for case 'roberts' corrected; div u, du_i/dx_j for case
!                   'roberts_xz' added
      use Diagnostics
      use General
      use Sub
!
      real, dimension (mx,my,mz,mfarray),intent(IN) :: f
      type (pencil_case),                intent(OUT):: p
      logical, dimension(npencils),      intent(IN) :: lpenc_loc
!
      real, dimension(nx) :: kdotxwt, cos_kdotxwt, sin_kdotxwt
      real, dimension(nx) :: local_Omega
      real, dimension(nx) :: wind_prof,div_uprof,der6_uprof
      real, dimension(nx) :: div_vel_prof
      real, dimension(nx) :: vel_prof
      real, dimension(nx) :: tmp_mn, cos1_mn, cos2_mn
      real, dimension(nx) :: rone, argx, pom2
      real, dimension(nx) :: psi1, psi2, psi3, psi4, rho_prof, prof, prof1
      real, dimension(nx) :: random_r, random_p, random_tmp
!      real :: random_r_pt, random_p_pt
      real :: fac, fac2, argy, argz, cxt, cyt, czt, omt, del
      real :: fpara, dfpara, ecost, esint, epst, sin2t, cos2t
      real :: sqrt2, sqrt21k1, eps1=1., WW=0.25, k21
      real :: Balpha, ABC_A1, ABC_B1, ABC_C1
      real :: ro
      real :: xi, slopei, zl1, zlm1, zmax, kappa_kinflow_n, nn_eff
      real :: theta,theta1
      real :: exp_kinflow1,exp_kinflow2
      integer :: modeN, ell, ll, nn, ii, nn_max
!
!  Choose from a list of different flow profiles.
!  Begin with a
!
      lupdate_aux=.true.
!
      select case (kinematic_flow)
!
!constant flow in the x direction.
!
      case ('const-x')
        if (headtt) print*,'const-x'
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=ampl_kinflow*cos(omega_kinflow*t)*exp(eps_kinflow*t)
          p%uu(:,2)=0.
          p%uu(:,3)=0.
        endif
        if (lpenc_loc(i_ugu)) p%ugu=0.
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_del2u)) p%del2u=0.
!
!constant flow in the
!  (ampl_kinflow_x,ampl_kinflow_y,ampl_kinflow_z) direction.
!
      case ('const-xyz')
        if (headtt) print*,'const-xyz'
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=ampl_kinflow_x*cos(omega_kinflow*t)*exp(eps_kinflow*t)
          p%uu(:,2)=ampl_kinflow_y*cos(omega_kinflow*t)*exp(eps_kinflow*t)
          p%uu(:,3)=ampl_kinflow_z*cos(omega_kinflow*t)*exp(eps_kinflow*t)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
! gradient flow, u_x = -lambda x ; u_y = lambda y
!
      case ('grad_xy')
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-ampl_kinflow*lambda_kinflow*x(l1:l2)
          p%uu(:,2)= ampl_kinflow*lambda_kinflow*y(m)
          p%uu(:,3)=0. 
        endif
!
!  ABC-flow
!
      case ('ABC')
        if (headtt) print*,'ABC flow'
! uu
        if (lpenc_loc(i_uu)) then
          ABC_A1=ABC_A*cos(omega_kinflow*t+phasex_uukin)
          ABC_B1=ABC_B*cos(omega_kinflow*t+phasey_uukin)
          ABC_C1=ABC_C*cos(omega_kinflow*t+phasez_uukin)
          p%uu(:,1)=ABC_A1*sin(kz_uukin*z(n))    +ABC_C1*cos(ky_uukin*y(m))
          p%uu(:,2)=ABC_B1*sin(kx_uukin*x(l1:l2))+ABC_A1*cos(kz_uukin*z(n))
          p%uu(:,3)=ABC_C1*sin(ky_uukin*y(m))    +ABC_B1*cos(kx_uukin*x(l1:l2))
        endif
        if (lpenc_loc(i_oo)) then
          p%oo(:,1)=ABC_A1*kz_uukin*sin(kz_uukin*z(n))    +ABC_C1*cos(ky_uukin*y(m))
          p%oo(:,2)=ABC_B1*kx_uukin*sin(kx_uukin*x(l1:l2))+ABC_A1*kz_uukin*cos(kz_uukin*z(n))
          p%oo(:,3)=ABC_C1*ky_uukin*sin(ky_uukin*y(m))    +ABC_B1*kx_uukin*cos(kx_uukin*x(l1:l2))
        endif
! divu
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Willis-flow
!
      case ('Willis')
        if (headtt) print*,'Willis flow'
        fac=2.*one_over_sqrt3
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=fac*sin(ky_uukin*y(m)    )*cos(kz_uukin*z(n)    )
          p%uu(:,2)=fac*sin(kz_uukin*z(n)    )*cos(kx_uukin*x(l1:l2))
          p%uu(:,3)=fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m)    )
        endif
! oo
        if (lpenc_loc(i_oo)) then
          p%oo(:,1)=-fac*(sin(kx_uukin*x(l1:l2))*ky_uukin*sin(ky_uukin*y(m    )) &
                         +cos(kx_uukin*x(l1:l2))*kz_uukin*cos(kz_uukin*z(n    )))
          p%oo(:,2)=-fac*(sin(ky_uukin*y(m    ))*kz_uukin*sin(kz_uukin*z(n    )) &
                         +cos(ky_uukin*y(m    ))*kx_uukin*cos(kx_uukin*x(l1:l2)))
          p%oo(:,3)=-fac*(sin(kz_uukin*z(n    ))*kx_uukin*sin(kx_uukin*x(l1:l2)) &
                         +cos(kz_uukin*z(n    ))*ky_uukin*cos(ky_uukin*y(m    )))
        endif
! divu
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  nocosine or Archontis flow
!
      case ('nocos')
          if (headtt) print*,'nocosine or Archontis flow'
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=ABC_A*sin(kz_uukin*z(n))
          p%uu(:,2)=ABC_B*sin(kx_uukin*x(l1:l2))
          p%uu(:,3)=ABC_C*sin(ky_uukin*y(m))
        endif
! divu
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Roberts I flow with negative helicity
!
      case ('roberts')
        if (headtt) print*,'Glen Roberts flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
! uu
        sqrt2=sqrt(2.)
        if (lpenc_loc(i_uu)) then
          eps1=1.-eps_kinflow
          p%uu(:,1)=+eps1*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,2)=-eps1*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,3)=sqrt2*sin(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
        endif
! divu
        if (lpenc_loc(i_divu)) &
            p%divu= eps1*(kx_uukin-ky_uukin)*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
! uij
        if (lpenc_loc(i_uij)) then
          p%uij(:,1,1)=+eps1*kx_uukin*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uij(:,1,2)=-eps1*ky_uukin*sin(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uij(:,1,3)=+0.
          p%uij(:,2,1)=+eps1*kx_uukin*sin(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uij(:,2,2)=-eps1*ky_uukin*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uij(:,2,3)=+0.
          p%uij(:,3,1)=sqrt2*kx_uukin*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uij(:,3,2)=sqrt2*ky_uukin*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uij(:,3,3)=+0.
        endif
!
      case ('roberts-xz')
        if (headtt) print*,'Glen Roberts flow w.r.t. x and z; kx_uukin,kz_uukin=',kx_uukin,kz_uukin
! uu
        sqrt2=sqrt(2.)
        if (lpenc_loc(i_uu)) then
          eps1=1.-eps_kinflow
          p%uu(:,1)=+eps1 *sin(kx_uukin*x(l1:l2))*cos(kz_uukin*z(n))
          p%uu(:,2)=-sqrt2*sin(kx_uukin*x(l1:l2))*sin(kz_uukin*z(n))
          p%uu(:,3)=-eps1 *cos(kx_uukin*x(l1:l2))*sin(kz_uukin*z(n))
        endif
! divu
        if (lpenc_loc(i_divu) .or. lpenc_loc(i_uij)) &
          call inevitably_fatal_error('hydro_kinematic', 'divu and uij not implemented for "roberts-xz"')

     case ('Roberts-II-xz')
        if (headtt) print*,'Glen Roberts flow II w.r.t. x and z; kx_uukin,kz_uukin=',kx_uukin,kz_uukin
! uu
        fac=ampl_kinflow
        fac2=ampl_kinflow*eps_kinflow
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)= +fac *sin(kx_uukin*x(l1:l2))*cos(kz_uukin*z(n))
          p%uu(:,2)= -fac2*cos(kx_uukin*x(l1:l2))*cos(kz_uukin*z(n))
          p%uu(:,3)= -fac *cos(kx_uukin*x(l1:l2))*sin(kz_uukin*z(n))
        endif
! divu
        if (lpenc_loc(i_divu)) &
            p%divu= fac*(kx_uukin-kz_uukin)*cos(kx_uukin*x(l1:l2))*cos(kz_uukin*z(n))
! uij
        if (lpenc_loc(i_uij)) then
          p%uij(:,1,1)=+fac *kx_uukin*cos(kx_uukin*x(l1:l2))*cos(kz_uukin*z(n))
          p%uij(:,1,2)=+0.
          p%uij(:,1,3)=-fac *kz_uukin*sin(kx_uukin*x(l1:l2))*sin(kz_uukin*z(n))
          p%uij(:,2,1)=+fac2*kx_uukin*sin(kx_uukin*x(l1:l2))*cos(kz_uukin*z(n))
          p%uij(:,2,2)=+0.
          p%uij(:,2,3)=+fac2*kz_uukin*cos(kx_uukin*x(l1:l2))*sin(kz_uukin*z(n))
          p%uij(:,3,1)=+fac *kx_uukin*sin(kx_uukin*x(l1:l2))*sin(kz_uukin*z(n))
          p%uij(:,3,2)=+0.
          p%uij(:,3,3)=-fac *kz_uukin*cos(kx_uukin*x(l1:l2))*cos(kz_uukin*z(n))
        endif
!
! Chandrasekhar-Kendall Flow
!
      case ('ck')
        if (headtt) print*,'Chandrasekhar-Kendall flow'
! uu
        ell=kinflow_ck_ell
        Balpha=kinflow_ck_Balpha
        ck_r = x(l1:l2)
        ck_rsqr = x(l1:l2)*x(l1:l2)
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=ampl_kinflow*Pl(m)*(  &
              (ell*(ell+1)/(Balpha*ck_rsqr))*Zl(l1:l2) &
              -(2./(Balpha*ck_r))*dZldr(l1:l2)  )
          p%uu(:,2)=ampl_kinflow*( &
              dZldr(l1:l2)/Balpha- Zl(l1:l2)/(Balpha*ck_r) &
              )*dPldtheta(m)
          p%uu(:,3)=-ampl_kinflow*Zl(l1:l2)*dPldtheta(m)
        endif
! divu
        if (lpenc_loc(i_divu)) p%divu= 0.
!
!  Roberts I flow with positive helicity.
!  For relhel_uukin=1, we have the normal Roberts flow.
!  For relhel_uukin=0., the flow is purely in closed circles.
!  For relhel_uukin=2., the flow is purely along the z direction.
!
      case ('poshel-roberts')
        if (headtt) print*,'Pos Helicity Roberts flow; chi_uukin=',chi_uukin
        fac =ampl_kinflow*cos(chi_uukin*dtor)*sqrt(2.)
        fac2=ampl_kinflow*sin(chi_uukin*dtor)*2.
        del =                 del_uukin*dtor
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2)    )*sin(ky_uukin*y(m)    )
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)    )*cos(ky_uukin*y(m)    )
          !p%uu(:,3)=fac2*cos(kx_uukin*x(l1:l2)+del)*cos(ky_uukin*y(m)+del)
          p%uu(:,3)=fac2*sin(kx_uukin*x(l1:l2)+del)*sin(ky_uukin*y(m)+del)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Roberts II flow (from Tilgner 2004)
!
      case ('Roberts-II')
        if (headtt) print*,'Roberts-II flow; eps_kinflow=',eps_kinflow
        fac=ampl_kinflow
        fac2=ampl_kinflow*eps_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,2)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,3)=fac2*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
        endif
!
! The following not true for kx_uukin \ne ky_uukin.
!
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  original Roberts III flow
!
      case ('Roberts-III-orig')
        if (headtt) print*,'original Roberts-III flow; eps_kinflow=',eps_kinflow
        fac=ampl_kinflow
        fac2=ampl_kinflow*eps_kinflow*2.
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=fac2*cos(ky_uukin*y(m))*cos(kz_uukin*z(n))
          p%uu(:,2)=+fac*sin(kz_uukin*z(n))
          p%uu(:,3)=+fac*sin(ky_uukin*y(m))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Roberts III flow (from Tilgner 2004)
!
      case ('Roberts-III')
        if (headtt) print*,'Roberts-III flow; eps_kinflow=',eps_kinflow
        fac=ampl_kinflow
        fac2=ampl_kinflow*eps_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,2)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,3)=fac2*(cos(kx_uukin*x(l1:l2))**2-sin(ky_uukin*y(m))**2)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Roberts IIIb flow (from Tilgner 2004)
!
      case ('Roberts-IIIb')
        if (headtt) print*,'Roberts-IIIb flow; eps_kinflow=',eps_kinflow
        fac=ampl_kinflow
        fac2=.5*ampl_kinflow*eps_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,2)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,3)=fac2*(cos(kx_uukin*x(l1:l2))+cos(ky_uukin*y(m)))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Modified Roberts flow (from Tilgner 2004)
!
      case ('Roberts-IV')
        if (headtt) print*,'Roberts-IV flow; eps_kinflow=',eps_kinflow
        if (eps_kinflow==0.) &
          call inevitably_fatal_error('hydro_kinematic','kinflow = "Roberts IV", '//&
                                      'eps_kinflow=0')
        fac=sqrt(2./eps_kinflow)*ampl_kinflow
        fac2=sqrt(eps_kinflow)*ampl_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,2)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,3)=+fac2*sin(kx_uukin*x(l1:l2))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_oo)) then
          p%oo(:,1)=0.
          p%oo(:,2)=-fac2*kx_uukin*cos(kx_uukin*x(l1:l2))
          p%oo(:,3)=2.*fac*kx_uukin*sin(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
        endif
!
!  Modified Roberts flow (from Tilgner 2004)
!
      case ('Roberts-IVb')
        if (headtt) print*,'Roberts-IV flow; eps_kinflow=',eps_kinflow
        fac=sqrt(2./eps_kinflow)*ampl_kinflow
        fac2=sqrt(eps_kinflow)*ampl_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,2)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,3)=+fac2*sin(ky_uukin*y(m))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_oo)) then
          p%oo(:,1)=fac2*kx_uukin*cos(ky_uukin*y(m))
          p%oo(:,2)=0.
          p%oo(:,3)=2.*fac*kx_uukin*sin(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
        endif
!
!  Glen-Roberts flow (positive helicity), alternative version
!
      case ('varhel-roberts')
        if (headtt) print*,'Pos Helicity Roberts flow; eps1=',eps1
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)*eps1
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  1-D Glen-Roberts flow (positive helicity, no y-dependence)
!
      case ('poshel-roberts-1d')
        if (headtt) print*,'Pos Helicity Roberts flow; kx_uukin=',kx_uukin
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*eps1
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*sqrt(2.)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Glen-Roberts flow (x-direction, positive helicity)
!  x -> y
!  y -> z
!  z -> x
!
      case ('xdir-roberts')
        if (headtt) print*,'x-dir Roberts flow; ky_uukin,kz_uukin=',ky_uukin,kz_uukin
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,2)=-fac*cos(ky_uukin*y(m))*sin(kz_uukin*z(n))*eps1
          p%uu(:,3)=+fac*sin(ky_uukin*y(m))*cos(kz_uukin*z(n))*eps1
          p%uu(:,1)=+fac*cos(ky_uukin*y(m))*cos(kz_uukin*z(n))*sqrt(2.)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  z-dependent Roberts flow (positive helicity)
!
      case ('zdep-roberts')
        if (headtt) print*,'z-dependent Roberts flow; kx,ky=',kx_uukin,ky_uukin
        fpara=ampl_kinflow*(quintic_step(z(n),-1.+eps_kinflow,eps_kinflow) &
            -quintic_step(z(n),+1.-eps_kinflow,eps_kinflow))
        dfpara=ampl_kinflow*(quintic_der_step(z(n),-1.+eps_kinflow,eps_kinflow)&
            -quintic_der_step(z(n),+1.-eps_kinflow,eps_kinflow))
!
        sqrt2=sqrt(2.)
        sqrt21k1=1./(sqrt2*kx_uukin)
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m)) &
              -dfpara*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt21k1
          p%uu(:,2)=+sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m)) &
              -dfpara*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*sqrt21k1
          p%uu(:,3)=+fpara*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt2
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with cosinusoidal helicity variation (version 1)
!
      case ('IncohRoberts1')
        if (headtt) print*,'Roberts flow with cosinusoidal helicity;',&
            ' kx_uukin,ky_uukin=',kx_uukin,ky_uukin
! uu
        if (lpenc_loc(i_uu)) then
          fac=ampl_kinflow
          eps1=(1.-eps_kinflow)*cos(omega_kinflow*t)
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*eps1
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*eps1
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with cosinusoidal helicity variation (version 2)
!
      case ('IncohRoberts2')
        if (headtt) print*,'Roberts flow with cosinusoidal helicity;',&
            'kx_uukin,ky_uukin=',kx_uukin,ky_uukin
! uu
        if (lpenc_loc(i_uu)) then
          fac=ampl_kinflow
          eps1=(1.-eps_kinflow)*cos(omega_kinflow*t)
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)*eps1
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with helicity variation and shear (version 1)
!  Must use shear module and set eps_kinflow equal to shear
!
      case ('ShearRoberts1')
        if (headtt) print*,'Roberts flow with cosinusoidal helicity;',&
            'kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        ky_uukin=1.
        kx_uukin=ky_uukin*(mod(.5-eps_kinflow*t,1.D0)-.5)
        if (ip==11.and.m==4.and.n==4) write(21,*) t,kx_uukin
        eps1=cos(omega_kinflow*t)
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*eps1
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*eps1
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with helicity variation and shear (version 2)
!  Must use shear module and set eps_kinflow equal to shear
!
      case ('ShearRoberts2')
        if (headtt) print*,'Roberts flow with cosinusoidal helicity;',&
            'kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        ky_uukin=1.
        kx_uukin=ky_uukin*(mod(.5-eps_kinflow*t,1.D0)-.5)
        if (ip==11.and.m==4.and.n==4) write(21,*) t,kx_uukin
        eps1=cos(omega_kinflow*t)
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)*eps1
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Another planar flow (from XXX 2015)
!
      case ('Another-Planar')
        if (headtt) print*,'Another planar flow; ampl_kinflow=',ampl_kinflow
        fac=ampl_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+fac*cos(ky_uukin*y(m))
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))
          p%uu(:,3)=+fac*(sin(kx_uukin*x(l1:l2))+cos(ky_uukin*y(m)))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_oo)) then
!place holder
          p%oo(:,1)=fac2*kx_uukin*cos(ky_uukin*y(m))
          p%oo(:,2)=0.
          p%oo(:,3)=2.*fac*kx_uukin*sin(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
        endif
!
!  Time-dependent, nearly symmetric flow
!
      case ('Herreman')
        if (headtt) print*,'Herreman flow;',&
            'kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        eps1=ampl_kinflow*eps_kinflow*cos(omega_kinflow*t)
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*sin(ky_uukin*y(m))
          p%uu(:,2)=eps1*sin(kx_uukin*x(l1:l2))
          p%uu(:,3)=+fac*cos(ky_uukin*y(m))+eps1*cos(kx_uukin*x(l1:l2))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Taylor-Green flow
!
      case ('TG')
        if (headtt) print*,'Taylor-Green flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
! uu
        fac=2.*cos(omega_kinflow*t)
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*cos(kz_uukin*z(n))
          p%uu(:,2)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*cos(kz_uukin*z(n))
          p%uu(:,3)=+0.
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_oo)) then
          p%oo(:,1)=-fac*kz_uukin*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*sin(kz_uukin*z(n))
          p%oo(:,2)=-fac*kz_uukin*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sin(kz_uukin*z(n))
          p%oo(:,3)=fac*(kx_uukin+ky_uukin)*sin(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*cos(kz_uukin*z(n))
        endif
!
!  modified Taylor-Green flow
!
      case ('TGmod')
        if (headtt) print*,'modified Taylor-Green flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*cos(kz_uukin*z(n)) &
            +ABC_A*sin(2*kx_uukin*x(l1:l2))*cos(2*kz_uukin*z(n)) &
            +ABC_B*(sin(kx_uukin*x(l1:l2))*cos(3*ky_uukin*y(m))*cos(kz_uukin*z(n)) &
            +5./13.*sin(3*kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*cos(kz_uukin*z(n)))
          p%uu(:,2)=-cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*cos(kz_uukin*z(n)) &
            +ABC_A*sin(2*ky_uukin*y(m))*cos(2*kz_uukin*z(n)) &
            -ABC_B*(cos(3*kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*cos(kz_uukin*z(n)) &
            +5./13.*cos(kx_uukin*x(l1:l2))*sin(3*ky_uukin*y(m))*cos(kz_uukin*z(n)))
          p%uu(:,3)=-ABC_A*(cos(2*kx_uukin*x(l1:l2))+cos(2*ky_uukin*y(m)))*sin(2*kz_uukin*z(n)) &
            +ABC_B*2./13.*(cos(kx_uukin*x(l1:l2))*cos(3*ky_uukin*y(m)) &
                          -cos(3*kx_uukin*x(l1:l2))*cos(ky_uukin*y(m)))*sin(kz_uukin*z(n))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  modified Taylor-Green flow with (x,y,z) -> tilde(y,z,x), i.e., the y direction became z
!
      case ('TGmodY')
        if (headtt) print*,'modified Taylor-Green flow-Y; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,2)=+sin(ky_uukin*y(m))*cos(kz_uukin*z(n))*cos(kx_uukin*x(l1:l2)) &
            +ABC_A*sin(2*ky_uukin*y(m))*cos(2*kx_uukin*x(l1:l2)) &
            +ABC_B*(sin(ky_uukin*y(m))*cos(3*kz_uukin*z(n))*cos(kx_uukin*x(l1:l2)) &
            +5./13.*sin(3*ky_uukin*y(m))*cos(kz_uukin*z(n))*cos(kx_uukin*x(l1:l2)))
          p%uu(:,3)=-cos(ky_uukin*y(m))*sin(kz_uukin*z(n))*cos(kx_uukin*x(l1:l2)) &
            +ABC_A*sin(2*kz_uukin*z(n))*cos(2*kx_uukin*x(l1:l2)) &
            -ABC_B*(cos(3*ky_uukin*y(m))*sin(kz_uukin*z(n))*cos(kx_uukin*x(l1:l2)) &
            +5./13.*cos(ky_uukin*y(m))*sin(3*kz_uukin*z(n))*cos(kx_uukin*x(l1:l2)))
          p%uu(:,1)=-ABC_A*(cos(2*ky_uukin*y(m))+cos(2*kz_uukin*z(n)))*sin(2*kx_uukin*x(l1:l2)) &
            +ABC_B*2./13.*(cos(ky_uukin*y(m))*cos(3*kz_uukin*z(n)) &
                          -cos(3*ky_uukin*y(m))*cos(kz_uukin*z(n)))*sin(kx_uukin*x(l1:l2))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  modified Taylor-Green flow with (y,z,x) -> tilde(z,x,y), i.e., the x direction became z
!
      case ('TGmodX')
        if (headtt) print*,'modified Taylor-Green flow-X; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,3)=+sin(kz_uukin*z(n))*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m)) &
            +ABC_A*sin(2*kz_uukin*z(n))*cos(2*ky_uukin*y(m)) &
            +ABC_B*(sin(kz_uukin*z(n))*cos(3*kx_uukin*x(l1:l2))*cos(ky_uukin*y(m)) &
            +5./13.*sin(3*kz_uukin*z(n))*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m)))
          p%uu(:,1)=-cos(kz_uukin*z(n))*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m)) &
            +ABC_A*sin(2*kx_uukin*x(l1:l2))*cos(2*ky_uukin*y(m)) &
            -ABC_B*(cos(3*kz_uukin*z(n))*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m)) &
            +5./13.*cos(kz_uukin*z(n))*sin(3*kx_uukin*x(l1:l2))*cos(ky_uukin*y(m)))
          p%uu(:,2)=-ABC_A*(cos(2*kz_uukin*z(n))+cos(2*kx_uukin*x(l1:l2)))*sin(2*ky_uukin*y(m)) &
            +ABC_B*2./13.*(cos(kz_uukin*z(n))*cos(3*kx_uukin*x(l1:l2)) &
                          -cos(3*kz_uukin*z(n))*cos(kx_uukin*x(l1:l2)))*sin(ky_uukin*y(m))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  modified Taylor-Green flow with (y,z,x) -> tilde(z,x,y), i.e., the x direction became z
!
      case ('Straining')
        if (headtt) print*,'Straining motion; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
! uu
        if (lpenc_loc(i_uu)) then
          fac=ampl_kinflow
          fac2=-(dimensionality-1)*fac
          argx=kx_uukin*x(l1:l2)+phasex_uukin
          argy=ky_uukin*y(m)+phasey_uukin
          argz=kz_uukin*z(n)+phasez_uukin
          p%uu(:,1)= fac*sin(argx)*cos(argy)*cos(argz)
          p%uu(:,2)= fac*cos(argx)*sin(argy)*cos(argz)
          p%uu(:,3)=fac2*cos(argx)*cos(argy)*sin(argz)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Beltrami-x flow
!
      case ('Beltrami-x')
        if (headtt) print*,'Beltrami-x motion; kx_uukin=',kx_uukin
! uu
        if (lpenc_loc(i_uu)) then
          fac=ampl_kinflow*cos(omega_kinflow*t)
          argx=kx_uukin*x(l1:l2)+phasex_uukin
          p%uu(:,1)=0.
          p%uu(:,2)=fac*sin(argx)*relhel_uukin
          p%uu(:,3)=fac*cos(argx)
        endif
        if (lpenc_loc(i_oo)) then
          p%oo(:,1)=0.
          p%oo(:,2)=fac*kx_uukin*sin(argx)
          p%oo(:,3)=fac*kx_uukin*cos(argx)*relhel_uukin
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Beltrami-y flow
!
      case ('Beltrami-y')
        if (headtt) print*,'Beltrami-y motion; ky_uukin=',ky_uukin
! uu
        if (lpenc_loc(i_uu)) then
          fac=ampl_kinflow*cos(omega_kinflow*t)
          argy=ky_uukin*y(m)+phasey_uukin
          p%uu(:,1)=fac*cos(argy)
          p%uu(:,2)=0.
          p%uu(:,3)=fac*sin(argy)*relhel_uukin
        endif
        if (lpenc_loc(i_oo)) then
          p%oo(:,1)=fac*ky_uukin*cos(argy)*relhel_uukin
          p%oo(:,2)=0.
          p%oo(:,3)=fac*ky_uukin*sin(argy)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Beltrami-z flow
!
      case ('Beltrami-z')
        if (headtt) print*,'Beltrami-z motion; kz_uukin=',kz_uukin
! uu
        if (lpenc_loc(i_uu)) then
          fac=ampl_kinflow*cos(omega_kinflow*t)
          argz=kz_uukin*z(n)+phasez_uukin
          p%uu(:,1)=fac*sin(argz)*relhel_uukin
          p%uu(:,2)=fac*cos(argz)
          p%uu(:,3)=0.
        endif
        if (lpenc_loc(i_oo)) then
          p%oo(:,1)=fac*kz_uukin*sin(argz)
          p%oo(:,2)=fac*kz_uukin*cos(argz)*relhel_uukin
          p%oo(:,3)=0.
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  sinusoidal 1-D flow
!
      case ('sinusoidal')
        if (headtt) print*,'sinusoidal motion; kz_uukin=',kz_uukin
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=ampl_kinflow*sin(kz_uukin*z(n))
        endif
        if (lpenc_loc(i_divu)) p%divu=ampl_kinflow*kz_uukin*cos(kz_uukin*z(n))
!
!  Galloway-Proctor flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!
      case ('Galloway-Proctor')
        if (headtt) print*,'Galloway-Proctor flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        ecost=eps_kinflow*cos(omega_kinflow*t)
        esint=eps_kinflow*sin(omega_kinflow*t)
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*sin(ky_uukin*y(m)    +esint)
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)
          p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+ecost)+cos(ky_uukin*y(m)+esint))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Galloway-Proctor flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!
      case ('Galloway-Proctor-nohel')
        if (headtt) print*,'nonhelical Galloway-Proctor flow;',&
            'kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow*sqrt(1.5)
        fac2=ampl_kinflow*sqrt(6.)
        ecost=eps_kinflow*cos(omega_kinflow*t)
        esint=eps_kinflow*sin(omega_kinflow*t)
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+fac*cos(ky_uukin*y(m)    +esint)
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)
          p%uu(:,3)=-fac2*(sin(kx_uukin*x(l1:l2)+ecost)*cos(ky_uukin*y(m)+esint))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Otani flow, U=curl(psi*zz) + psi*zz, where
!  psi = 2*cos^2t * cosx - 2*csin2t * cosy
!
      case ('Otani')
        if (headtt) print*,'Otani flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=2.*ampl_kinflow
        sin2t=sin(omega_kinflow*t)**2
        cos2t=cos(omega_kinflow*t)**2
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=fac*sin2t*sin(ky_uukin*y(m))
          p%uu(:,2)=fac*cos2t*sin(kx_uukin*x(l1:l2))
          p%uu(:,3)=fac*(cos2t*cos(kx_uukin*x(l1:l2))-sin2t*cos(ky_uukin*y(m)))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_oo)) p%oo=p%uu
!
!  Tilgner flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!
      case ('Tilgner')
        if (headtt) print*,'Tilgner flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow*sqrt(2.)
        epst=eps_kinflow*t
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac* sin(ky_uukin*y(m)         )
          p%uu(:,2)=+fac* sin(kx_uukin*(x(l1:l2)+epst))
          p%uu(:,3)=-fac*(cos(kx_uukin*(x(l1:l2)+epst))+cos(ky_uukin*y(m)))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Tilgner flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!  Here, W in Tilgner's equation is chosen to be 0.25.
!
      case ('Tilgner-orig')
        if (headtt) print*,'original Tilgner flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        epst=eps_kinflow*t
        kx_uukin=2.*pi
        ky_uukin=2.*pi
        sqrt2=sqrt(2.)
        WW=0.25
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+fac*sqrt2*sin(kx_uukin*(x(l1:l2)-epst))*cos(ky_uukin*y(m))
          p%uu(:,2)=-fac*sqrt2*cos(kx_uukin*(x(l1:l2)-epst))*sin(ky_uukin*y(m))
          p%uu(:,3)=+fac*2.*WW*sin(kx_uukin*(x(l1:l2)-epst))*sin(ky_uukin*y(m))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Tree-like flow
!  Define ampl_kinflow > 0 for downflow; therefore the minus sign below.
!
      case ('FatTree')
        if (headtt) print*,'Tree flow; kx_uukin,Lx,Lz=',kx_uukin,Lx,Lz
        zmax=Lxyz(3)*(1.-2./2**tree_lmax)
        nn_max=2**tree_lmax
        fac=-ampl_kinflow*((zinfty_kinflow-z(n))/Lz)**(-1.5)
! uu
        if (lpenc_loc(i_uu)) then
          ll=int(alog(2.*Lxyz(3)/(xyz1(3)-min(z(n),zmax)))/alog(2.))
          nn=2**ll
          zl1 =(xyz1(3)-xyz0(3))*(1.-1./nn)  !(=z_l)
          zlm1=(xyz1(3)-xyz0(3))*(1.-2./nn)  !(=z_{l-1})
          prof=0.
          prof1=0.
          do ii=1,nn
            xi=real(ii)/nn-.5-.5/nn
            slopei=.5*(-1)**ii/nn
            theta=xi-slopei*(zl1-z(n))/(zl1-zlm1)
            argx=x(l1:l2)-Lx*theta
            nn_eff=1./(1.-min(z(n),zmax))
            kappa_kinflow_n=kappa_kinflow*(nn_eff/real(nn_max))**2
            if (ip.le.6) write(*,fmt='(2f10.4)') z(n),theta
            prof=prof+(.5+.5*cos(kx_uukin*argx))**kappa_kinflow_n/nn
            prof1=prof1+(.5+.5*cos(kx_uukin*argx))**kappa_kinflow_n/nn*(-1.)**ii
          enddo
          p%uu(:,1)=fac*prof1*Lx/(2.*Lz)
          p%uu(:,2)=0.
          p%uu(:,3)=fac*prof
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        !if (lpenc_loc(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Tree-like flow
!  Define ampl_kinflow > 0 for downflow; therefore the minus sign below.
!
      case ('Tree')
        if (headtt) print*,'Tree flow; kx_uukin,Lx,Lz=',kx_uukin,Lx,Lz
        zmax=Lxyz(3)*(1.-2./2**tree_lmax)
        fac=-ampl_kinflow*((zinfty_kinflow-z(n))/Lz)**(-1.5)
! uu
        if (lpenc_loc(i_uu)) then
          ll=int(alog(2.*Lxyz(3)/(xyz1(3)-min(z(n),zmax)))/alog(2.))
          nn=2**ll
          zl1 =(xyz1(3)-xyz0(3))*(1.-1./nn)  !(=z_l)
          zlm1=(xyz1(3)-xyz0(3))*(1.-2./nn)  !(=z_{l-1})
          prof=0.
          prof1=0.
          do ii=1,nn
            xi=real(ii)/nn-.5-.5/nn
            slopei=.5*(-1)**ii/nn
            theta=xi-slopei*(zl1-z(n))/(zl1-zlm1)
            argx=x(l1:l2)-Lx*theta
            if (ip.le.6) write(*,fmt='(2f10.4)') z(n),theta
            prof=prof+(.5+.5*cos(kx_uukin*argx))**kappa_kinflow/nn
            prof1=prof1+(.5+.5*cos(kx_uukin*argx))**kappa_kinflow/nn*(-1.)**ii
          enddo
          p%uu(:,1)=fac*prof1*Lx/(2.*Lz)
          p%uu(:,2)=0.
          p%uu(:,3)=fac*prof
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        !if (lpenc_loc(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Fence-like flow
!  Define ampl_kinflow > 0 for downflow; therefore the minus sign below.
!
      case ('Fence')
        if (headtt) print*,'Fence flow; kx_uukin,Lx,Lz=',kx_uukin,Lx,Lz
        zmax=Lxyz(3)*(1.-2./2**tree_lmax)
        if (zinfty_kinflow==0.) then
          fac=-ampl_kinflow
        else
          fac=-ampl_kinflow*((zinfty_kinflow-z(n))/Lz)**(-1.5)
        endif
! uu
        if (lpenc_loc(i_uu)) then
          argx=x(l1:l2)
          prof=(.5+.5*cos(kx_uukin*argx))**kappa_kinflow
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=fac*prof
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        !if (lpenc_loc(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Galloway-Proctor flow with random temporal phase
!
      case ('Galloway-Proctor-RandomTemporalPhase')
        if (headtt) print*,'GP-RandomTemporalPhase; kx,ky=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        if (t>tphase_kinflow) then
          call random_number_wrapper(fran1)
          tphase_kinflow=t+dtphase_kinflow
          phase1=pi*(2*fran1(1)-1.)
          phase2=pi*(2*fran1(2)-1.)
        endif
        ecost=eps_kinflow*cos(omega_kinflow*t+phase1)
        esint=eps_kinflow*sin(omega_kinflow*t+phase2)
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*sin(ky_uukin*y(m)    +esint)
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)
          p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+ecost)+cos(ky_uukin*y(m)+esint))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Galloway-Proctor flow with random phase
!
      case ('Galloway-Proctor-RandomPhase')
        if (headtt) print*,'Galloway-Proctor-RandomPhase; kx,ky=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        if (t>tphase_kinflow) then
          call random_number_wrapper(fran1)
          tphase_kinflow=t+dtphase_kinflow
          phase1=eps_kinflow*pi*(2*fran1(1)-1.)
          phase2=eps_kinflow*pi*(2*fran1(2)-1.)
        endif
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*sin(ky_uukin*y(m)    +phase1)*ky_uukin
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+phase2)*kx_uukin
          p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+phase2)+cos(ky_uukin*y(m)+phase1))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  original Galloway-Proctor flow
!
      case ('Galloway-Proctor-orig')
        if (headtt) print*,'Galloway-Proctor-orig flow; kx_uukin=',kx_uukin
        fac=sqrt(1.5)*ampl_kinflow
        ecost=eps_kinflow*cos(omega_kinflow*t)
        esint=eps_kinflow*sin(omega_kinflow*t)
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+fac*cos(ky_uukin*y(m)    +esint)*ky_uukin
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)*kx_uukin
          p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+ecost)+sin(ky_uukin*y(m)+esint))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Potential flow, u=gradphi, with phi=coskx*X cosky*Y coskz*Z,
!  and X=x-ct, Y=y-ct, Z=z-ct.
!  Possibility of gaussian distributed random amplitudes if lrandom_ampl.
!
      case ('potential')
!
!  Allow for harmonic phase changes.
!
        if (headtt) print*,'potential; ampl_kinflow=', ampl_kinflow
        if (headtt) print*,'potential; ki_uukin=',kx_uukin,ky_uukin,kz_uukin
        random_tmp=random_ampl*ampl_kinflow*cos(omega_kinflow*t)
        cxt=cx_uukin*t
        cyt=cy_uukin*t
        czt=cz_uukin*t
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-random_tmp*kx_uukin*&
              sin(kx_uukin*(x(l1:l2)-cxt))*cos(ky_uukin*(y(m)-cyt))*&
              cos(kz_uukin*(z(n)-czt-phasez_uukin))
          p%uu(:,2)=-random_tmp*ky_uukin*&
              cos(kx_uukin*(x(l1:l2)-cxt))*sin(ky_uukin*(y(m)-cyt))*&
              cos(kz_uukin*(z(n)-czt-phasez_uukin))
          p%uu(:,3)=-random_tmp*kz_uukin*&
              cos(kx_uukin*(x(l1:l2)-cxt))*cos(ky_uukin*(y(m)-cyt))*&
              sin(kz_uukin*(z(n)-czt-phasez_uukin))
        endif
        if (lpenc_loc(i_divu)) p%divu=-random_tmp*(kx_uukin**2+ky_uukin**2+kz_uukin**2) &
            *cos(kx_uukin*x(l1:l2)-cxt)*cos(ky_uukin*y(m)-cyt)*cos(kz_uukin*z(n)-czt)
!
!  2nd Potential flow, u=gradphi, with phi=cos(kx*X+ky*Y+kz*Z),
!  and X=x-ct, Y=y-ct, Z=z-ct.
!
      case ('potential2')
        if (headtt) print*,'2nd potential; ampl_kinflow=',ampl_kinflow
        if (headtt) print*,'2nd potential; ki_uukin=',kx_uukin,ky_uukin,kz_uukin
        fac=ampl_kinflow
        cxt=cx_uukin*t
        cyt=cy_uukin*t
        czt=cz_uukin*t
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*kx_uukin*&
              sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)+phasez_uukin)
          p%uu(:,2)=-fac*ky_uukin*&
              sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)+phasez_uukin)
          p%uu(:,3)=-fac*kz_uukin*&
            sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)+phasez_uukin)
        endif
        if (lpenc_loc(i_divu)) p%divu=-fac*(kx_uukin**2+ky_uukin**2+kz_uukin**2) &
            *cos(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)+phasez_uukin)
!
!  Incompressible 2D, u=curl(yy*psi), with psi=cos(kx*X+ky*Y+kz*Z),
!  and X=x-ct, Y=y-ct, Z=z-ct.
!
      case ('incompressible-2D-xz')
        if (headtt) print*,'incompr, 2D; ampl_kinflow=',ampl_kinflow
        if (headtt) print*,'incompr, 2D; ki_uukin=',kx_uukin,ky_uukin,kz_uukin
        fac=ampl_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+fac*kz_uukin*sin(kx_uukin*x(l1:l2)+kz_uukin*z(n))
          p%uu(:,2)=+0.
          p%uu(:,3)=-fac*kx_uukin*sin(kx_uukin*x(l1:l2)+kz_uukin*z(n))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  1-D Potential flow, u=gradphi, with phi=cos(kx*X+ky*Y+kz*Z),
!  and X=x-ct, Y=y-ct, Z=z-ct.
!
      case ('potentialz')
        if (headtt) print*,'1-D potential; ampl_kinflow,omega_kinflow=',&
            ampl_kinflow,omega_kinflow
        if (headtt) print*,'1-D potential; ki_uukin=',kx_uukin,ky_uukin,kz_uukin
        fac=ampl_kinflow
        omt=omega_kinflow*t
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*kx_uukin*&
              sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)-omt)
          p%uu(:,2)=-fac*ky_uukin*&
              sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)-omt)
          p%uu(:,3)=-fac*kz_uukin*&
              sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)-omt)
        endif
        if (lpenc_loc(i_divu)) p%divu=-fac*(kx_uukin**2+ky_uukin**2+kz_uukin**2) &
            *cos(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)-omt)
!
!  Potential random flow, u=gradphi, with phi=cos(x-x0)*cosy*cosz;
!  assume kx_uukin=ky_uukin=kz_uukin.
!
      case ('potential_random')
        if (headtt) print*,'potential_random; kx_uukin,ampl_kinflow=',ampl_kinflow
        if (headtt) print*,'potential_random; kx_uukin=',kx_uukin,ky_uukin,kz_uukin
        fac=random_ampl*ampl_kinflow
        argx=random_wavenumber*kx_uukin*(x(l1:l2)-location(1))
        argy=random_wavenumber*ky_uukin*(y(m)-location(2))
        argz=random_wavenumber*kz_uukin*(z(n)-location(3))
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*kx_uukin*sin(argx)*cos(argy)*cos(argz)
          p%uu(:,2)=-fac*ky_uukin*cos(argx)*sin(argy)*cos(argz)
          p%uu(:,3)=-fac*kz_uukin*cos(argx)*cos(argy)*sin(argz)
        endif
        if (lpenc_loc(i_divu)) p%divu=fac
!
!  Potential random flow, u=gradphi, with phi=cos(x-x0)*cosy*cosz;
!  assume kx_uukin=ky_uukin=kz_uukin.
!
      case ('potential_ampl_random')
        if (headtt) print*,'potential_random; kx_uukin,ampl_kinflow=',ampl_kinflow
        if (headtt) print*,'potential_random; kx_uukin=',kx_uukin,ky_uukin,kz_uukin
        fac=ampl_kinflow*ampl_random
        !argx=kx_uukin*(x(l1:l2)-location(1))
        !argy=ky_uukin*(y(m)-location(2))
        !argz=kz_uukin*(z(n)-location(3))
        argx=kx_uukin*x(l1:l2)
        argy=ky_uukin*y(m)
        argz=kz_uukin*z(n)
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*kx_uukin*sin(argx)*cos(argy)*cos(argz)
          p%uu(:,2)=-fac*ky_uukin*cos(argx)*sin(argy)*cos(argz)
          p%uu(:,3)=-fac*kz_uukin*cos(argx)*cos(argy)*sin(argz)
          lupdate_aux=.true.
        endif
        if (lpenc_loc(i_divu)) p%divu=fac
!
!  Random flow that is delta correlated in space and time.
!
      case ('delta_correlated')
        if (headtt) print*,'delta_correlated; ampl_kinflow=',ampl_kinflow
        if (lpenc_loc(i_uu)) then
          do ii=1,3
            if (modulo(ii-1,2)==0) then
              call random_number_wrapper(random_r)
              call random_number_wrapper(random_p)
              random_tmp=sqrt(-2*log(random_r))*sin(2*pi*random_p)
            else
              random_tmp=sqrt(-2*log(random_r))*cos(2*pi*random_p)
            endif
            p%uu(:,ii)=one_over_sqrt3*ampl_kinflow*random_tmp
          enddo
        endif
!
!  Random flow that is delta correlated in space, but unchanged in time.
!
      case ('deltax_correlated')
        if (headtt) print*,'delta_correlated_fixed_intime; ampl_kinflow=',ampl_kinflow
        if (lpenc_loc(i_uu)) then
          do ii=1,3
            if (modulo(ii-1,2)==0) then
              seed(1)=-((seed0-1812+1)*10+iproc_world+ncpus*(m+my*n))
              call random_seed_wrapper(PUT=seed)
              call random_number_wrapper(random_r)
              call random_number_wrapper(random_p)
              random_tmp=sqrt(-2*log(random_r))*sin(2*pi*random_p)
            else
              random_tmp=sqrt(-2*log(random_r))*cos(2*pi*random_p)
            endif
            p%uu(:,ii)=one_over_sqrt3*ampl_kinflow*random_tmp
          enddo
        endif
!
!  Convection rolls
!  Stream function: psi_y = cos(kx*x) * cos(kz*z)
!
      case ('rolls')
        if (headtt) print*,'Convection rolls; kx_kinflow,kz_uukin=',kx_kinflow,kz_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+ampl_kinflow*kz_kinflow*cos(kx_kinflow*x(l1:l2))*sin(kz_kinflow*z(n))
          p%uu(:,2)=+0.
          p%uu(:,3)=-ampl_kinflow*kx_kinflow*sin(kx_kinflow*x(l1:l2))*cos(kz_kinflow*z(n))
        endif
! divu
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Convection rolls
!  Stream function: psi_y = sin(kx*x) * sin(kz*z)
!
      case ('rolls2')
        if (headtt) print*,'Convection rolls2; kx_kinflow,kz_uukin=',kx_kinflow,kz_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-ampl_kinflow*kz_kinflow*sin(kx_kinflow*x(l1:l2))*cos(kz_kinflow*z(n))
          p%uu(:,2)=+0.
          p%uu(:,3)=+ampl_kinflow*kx_kinflow*cos(kx_kinflow*x(l1:l2))*sin(kz_kinflow*z(n))
        endif
! divu
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Twist (Yousef & Brandenburg 2003)
!
      case ('Twist')
        if (headtt) print*,'Twist flow; eps_kinflow,kx=',eps_kinflow,kx_uukin
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=eps_kinflow*z(n)*sin(kx_uukin*x(l1:l2))
          p%uu(:,3)=1.+cos(kx_uukin*x(l1:l2))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Eddy (Brandenburg & Zweibel 1994)
!
      case ('eddy')
        if (headtt) print*,'eddy flow; eps_kinflow,kx=',eps_kinflow,kx_uukin
        cos1_mn=max(cos(.5*pi*p%rcyl_mn),0.)
        cos2_mn=cos1_mn**2
        tmp_mn=-.5*pi*p%rcyl_mn1*sin(.5*pi*p%rcyl_mn)*ampl_kinflow* &
            4.*cos2_mn*cos1_mn
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+tmp_mn*y(m)
          p%uu(:,2)=-tmp_mn*x(l1:l2)
          p%uu(:,3)=eps_kinflow*ampl_kinflow*cos2_mn**2
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Shearing wave
!
      case ('ShearingWave')
        if (headtt) print*,'ShearingWave flow; Sshear,eps_kinflow=',Sshear,eps_kinflow
        ky_uukin=1.
        kx_uukin=-ky_uukin*Sshear*t
        k21=1./(kx_uukin**2+ky_uukin**2)
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-ky_uukin*k21*cos(kx_uukin*x(l1:l2)+ky_uukin*y(m))
          p%uu(:,2)=+kx_uukin*k21*cos(kx_uukin*x(l1:l2)+ky_uukin*y(m))
          p%uu(:,3)=eps_kinflow*cos(kx_uukin*x(l1:l2)+ky_uukin*y(m))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Helical Shearing wave
!
      case ('HelicalShearingWave')
        if (headtt) print*,'ShearingWave flow; Sshear,eps_kinflow=',Sshear,eps_kinflow
        ky_uukin=1.
        kx_uukin=-ky_uukin*Sshear*t
        k21=1./(kx_uukin**2+ky_uukin**2)
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*eps1
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*eps1
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
! step function along z
!
      case ('zstep')
        if (headtt) print*,'wind:step function along z'
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=wind_amp*step(z(n),wind_rmin,wind_step_width)
        endif
!
! uniform radial shear with a cutoff at rmax
! uniform radial shear with saturation.
!
      case ('rshear-sat')
        if (headtt) print*,'radial shear saturated'
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
!          omega0=uphi_at_rmax/uphi_rmax
!          shear = arctanh(omega0)/(uphi_rmax-x(l1))
          p%uu(:,3)=ampl_kinflow*x(l1:l2)*sinth(m)*tanh(10.*(x(l1:l2)-x(l1)))
!*tanh(10.*(x(l1:l2)-x(l1)))
!          write(*,*)'DM',m,n,p%uu(:,3)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  U_phi = sin(theta)
!
      case ('Uz=siny')
        if (headtt) print*,'U_phi = sin(theta)',ampl_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=ampl_kinflow*sinth(m)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  U_phi = r*sin(theta)*Omega(r) with Omega(r)=1-r
!
      case ('(1-x)*x*siny')
        if (headtt) print*,'Omega=Omega0*(1-r), Omega0=',ampl_kinflow
        local_Omega=ampl_kinflow*(1.-x(l1:l2))
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=local_Omega*x(l1:l2)*sinth(m)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  U_phi adapted from
!  Ghizaru-Charbonneau-Smolarkiewicz (ApJL 715:L133-L137, 2010)
!
      case ('gcs')
        if (headtt) print*,'gcs:gcs_rzero ',gcs_rzero
        fac=ampl_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          local_Omega=fac*exp(-((x(l1:l2)-xyz1(1))/gcs_rzero)**2- &
              ((pi/2-y(m))/gcs_psizero**2))
          p%uu(:,3)= local_Omega*x(l1:l2)*sinth(m)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Brandt rotation curve (cylindrical geometry)
!
      case ('Brandt')
        if (lcylindrical_coords) then
          if (headtt) print*,'Brandt (cylindrical coords)',ampl_kinflow
! uu
          if (lpenc_loc(i_uu)) then
            local_Omega=ampl_kinflow*profx_kinflow1
            p%uu(:,1)=0.
            p%uu(:,2)=local_Omega*x(l1:l2)
            p%uu(:,3)=0.
          endif
        else
          if (headtt) print*,'Brandt (Cartesian)',ampl_kinflow
            exp_kinflow1=1./exp_kinflow
            exp_kinflow2=.5*exp_kinflow
            pom2=x(l1:l2)**2+y(m)**2
            profx_kinflow1=+1./(1.+(pom2/uphi_rbot**2)**exp_kinflow2)**exp_kinflow1
! uu
          if (lpenc_loc(i_uu)) then
            local_Omega=ampl_kinflow*profx_kinflow1
            p%uu(:,1)=-local_Omega*y(m)
            p%uu(:,2)=+local_Omega*x(l1:l2)
            p%uu(:,3)=0.
          endif
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Spoke-like differential rotation profile
!
      case ('spoke-like')
        if (headtt) print*,'spoke-like ',ampl_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          local_Omega=ampl_kinflow*profx_kinflow1*profy_kinflow1(m)
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=local_Omega*x(l1:l2)*sinth(m)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Spoke-like differential rotation profile with near-surface shear layer
!
      case ('spoke-like-NSSL')
        if (headtt) print*,'spoke-like-NSSL',ampl_kinflow
! uu
        if (lpenc_loc(i_uu)) then
          local_Omega=ampl_kinflow*profx_kinflow1*profy_kinflow1(m) &
                     +ampl_kinflow*profx_kinflow2*profy_kinflow2(m) &
                     +ampl_kinflow*profx_kinflow3*profy_kinflow3(m)
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=local_Omega*x(l1:l2)*sinth(m)
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
!
!  Vertical wind
!
      case ('vertical-wind')
        if (headtt) print*,'vertical-wind along z'
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=wind_amp*z(n)
        endif
        if (lpenc_loc(i_divu)) p%divu=wind_amp
!
!  Vertical wind that goes to zero for z < 0.
!
      case ('vertical-wind-surface')
        if (headtt) print*,'vertical-wind along z'
! uu
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=wind_amp*z(n)*0.5*(1.+erfunc((z(n)/wind_step_width)))
        endif
        if (lpenc_loc(i_divu)) p%divu=wind_amp
!
! Radial wind
!
      case ('radial-wind')
        if (headtt) print*,'Radial wind'
        select case (wind_profile)
        case ('none'); wind_prof=0.;div_uprof=0.
        case ('constant'); wind_prof=1.;div_uprof=0.
        case ('radial-step')
          wind_prof=step(x(l1:l2),wind_rmin,wind_step_width)
          div_uprof=der_step(x(l1:l2),wind_rmin,wind_step_width)
!          der6_uprof=der6_step(x(l1:l2),wind_rmin,wind_step_width)
          der6_uprof=0.
        case default;
          call inevitably_fatal_error('hydro_kinematic', 'kinflow = "radial wind" - '//&
                                      'no such wind profile')
        endselect
!
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=wind_amp*wind_prof
          p%uu(:,2)=0.
          p%uu(:,3)=0.
        endif
        if (lpenc_loc(i_divu)) p%divu=wind_amp*div_uprof
        if (lpenc_loc(i_der6u)) then
          p%der6u(:,1)=wind_amp*der6_uprof
          p%der6u(:,2)= 0.
          p%der6u(:,3)= 0.
        endif
!
!  meridional circulation; psi=.5*(x-x0)*(x-x1)*(y-y0)*(y-y1), so
!  ux=+dpsi/dy=+(x-x0)*(x-x1)*y
!  uy=-dpsi/dx=-x*(y-y0)*(y-y1)
!
      case ('circ_cartesian')
        if (headtt) print*,'just circulation'
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+circ_amp*(x(l1:l2)-xyz0(1))*(x(l1:l2)-xyz1(1))*y(m)
          p%uu(:,2)=-circ_amp*x(l1:l2)*(y(m)-xyz0(2))*(y(m)-xyz1(2))
          p%uu(:,3)=0.
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_der6u)) then
          p%der6u(:,1)=wind_amp*der6_uprof
          p%der6u(:,2)= 0.
          p%der6u(:,3)= 0.
        endif
!
      case ('circ_cartesian_rho1')
        if (headtt) print*,'circulation with 1/rho dependency'
        if (lpenc_loc(i_uu)) then
          rho_prof=(1.33333/(x(l1:l2)+1.13333)-0.97)**1.5
          p%uu(:,1)=+circ_amp*(x(l1:l2)-xyz0(1))*(x(l1:l2)-xyz1(1))*y(m)/rho_prof
          p%uu(:,2)=-circ_amp*x(l1:l2)*(y(m)-xyz0(2))*(y(m)-xyz1(2))/rho_prof
          p%uu(:,3)=0.
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_der6u)) then
          p%der6u(:,1)=wind_amp*der6_uprof
          p%der6u(:,2)= 0.
          p%der6u(:,3)= 0.
        endif
!
!  meridional circulation; psi=.5*(x-x0)*(x-x1)*(y-y0)*(y-y1), so
!  ux=+dpsi/dy=+(x-x0)*(x-x1)*y
!  uy=-dpsi/dx=-x*(y-y0)*(y-y1)
!
      case ('circ_cartesian_xz')
        if (headtt) print*,'just circulation'
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=-(x(l1:l2)-xyz0(1))*(x(l1:l2)-xyz1(1))*z(n)
          p%uu(:,2)=-0.
          p%uu(:,3)=+x(l1:l2)*(z(n)-xyz0(3))*(z(n)-xyz1(3))
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_der6u)) then
          p%der6u(:,1)=wind_amp*der6_uprof
          p%der6u(:,2)= 0.
          p%der6u(:,3)= 0.
        endif
!
!  meridional circulation
!
      case ('circ_spherical')
        if (headtt) print*,'just circulation (spherical)'
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=+(x(l1:l2)-xyz0(1))*(x(l1:l2)-xyz1(1))*y(m)
          p%uu(:,2)=-x(l1:l2)*(y(m)-xyz0(2))*(y(m)-xyz1(2))
          p%uu(:,3)=0.
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_der6u)) then
          p%der6u(:,1)=wind_amp*der6_uprof
          p%der6u(:,2)= 0.
          p%der6u(:,3)= 0.
        endif
!
      case ('mer_flow_dg11')
        if (headtt) print*,'meridional flow as in Dikpati & Gilman 2011 (spherical)'
        if (lpenc_loc(i_uu)) then
          rho_prof=(1./x(l1:l2)-0.97)**1.5
          ro=(1.-0.6)/5.
          psi1=sin(pi*(x(l1:l2)-rp)/(1.-rp))
          psi2=1.-exp(-1.*x(l1:l2)*y(m)**2.0000001)
          psi3=1.-exp(4.*x(l1:l2)*(y(m)-0.5*pi))
          psi4=exp(-(x(l1:l2)-ro)**2/gamma_dg11**2)
          p%uu(:,1)=-(psi1*psi3*psi4*exp(-1.*x(l1:l2)*y(m)**2.0000001) &
               *(1.*2.0000001*x(l1:l2)*y(m)**(2.0000001-1)) &
               +psi1*psi2*psi4*(-4.*x(l1:l2)*exp(4.*x(l1:l2)*(y(m)-0.5*pi)))) &
               /(x(l1:l2)**2*rho_prof*sin(y(m)))
          p%uu(:,2)=(cos(pi*(x(l1:l2)-rp)/(1.-rp))*pi/(1.-rp)*psi2*psi3*psi4 &
               -exp(-1.*x(l1:l2)*y(m)**2.0000001)*(-1.*y(m)**2.0000001)*psi1*psi3*psi4 &
               -exp(4.*x(l1:l2)*(y(m)-0.5*pi))*(4.*(y(m)-0.5*pi))*psi1*psi2*psi4 &
               -2.*(x(l1:l2)-ro)*psi1*psi2*psi3*psi4/gamma_dg11**2)/(sin(y(m))*x(l1:l2)*rho_prof)
          p%uu(:,3)=0.
        endif
        if (lpenc_loc(i_divu)) p%divu=0.
        if (lpenc_loc(i_der6u)) then
          p%der6u(:,1)=wind_amp*der6_uprof
          p%der6u(:,2)= 0.
          p%der6u(:,3)= 0.
        endif
!
! Radial wind+circulation
!
      case ('radial-wind+circ')
        if (headtt) print*,'Radial wind and circulation'
        select case (wind_profile)
        case ('none'); wind_prof=0.;div_uprof=0.
        case ('constant'); wind_prof=1.;div_uprof=0.
        case ('radial-step')
          wind_prof=step(x(l1:l2),wind_rmin,wind_step_width)
          div_uprof=der_step(x(l1:l2),wind_rmin,wind_step_width)
          der6_uprof=0.
        case default;
          call inevitably_fatal_error('hydro_kinematic', 'kinflow="radial_wind-circ" - '//&
                                      'no such wind profile')
        endselect
        rone=xyz0(1)
        theta=y(m)
        theta1=xyz0(2)
        if (lpenc_loc(i_uu)) then
          vel_prof=circ_amp*(1+stepdown(x(l1:l2),circ_rmax,circ_step_width))
          div_vel_prof=-der_step(x(l1:l2),circ_rmax,circ_step_width)
          p%uu(:,1)=vel_prof*(r1_mn**2)*(sin1th(m))*(&
              2*sin(theta-theta1)*cos(theta-theta1)*cos(theta)&
              -sin(theta)*sin(theta-theta1)**2)*&
              (x(l1:l2)-1.)*(x(l1:l2)-rone)**2 + &
              wind_amp*wind_prof
          p%uu(:,2)=-vel_prof*r1_mn*sin1th(m)*(&
              cos(theta)*sin(theta-theta1)**2)*&
              (x(l1:l2)-rone)*(3*x(l1:l2)-rone-2.)
          p%uu(:,3)=0.
        endif
        if (lpenc_loc(i_divu)) p%divu=div_uprof*wind_amp + &
            div_vel_prof*(r1_mn**2)*(sin1th(m))*(&
            2*sin(theta-theta1)*cos(theta-theta1)*cos(theta)&
            -sin(theta)*sin(theta-theta1)**2)*&
            (x(l1:l2)-1.)*(x(l1:l2)-rone)**2
        if (lpenc_loc(i_der6u)) then
          p%der6u(:,1)=wind_amp*der6_uprof
          p%der6u(:,2)= 0.
          p%der6u(:,3)= 0.
        endif
!
!  Radial shear + radial wind + circulation.
!
      case ('rshear-sat+rwind+circ')
        if (headtt) print*,'radial shear, wind, circulation: not complete yet'
!
!  First set wind and cirulation.
!
        select case (wind_profile)
        case ('none'); wind_prof=0.;div_uprof=0.
        case ('constant'); wind_prof=1.;div_uprof=0.
!
        case ('radial-step')
          wind_prof=step(x(l1:l2),wind_rmin,wind_step_width)
          div_uprof=der_step(x(l1:l2),wind_rmin,wind_step_width)
          der6_uprof=0.
        case default;
          call inevitably_fatal_error('hydro_kinematic', 'kinflow="radial-shear+rwind+rcirc" - '//&
                                      'no such wind profile')
        endselect
        rone=xyz0(1)
        theta=y(m)
        theta1=xyz0(2)
        if (lpenc_loc(i_uu)) then
          vel_prof=circ_amp*(1+stepdown(x(l1:l2),circ_rmax,circ_step_width))
          div_vel_prof=-der_step(x(l1:l2),circ_rmax,circ_step_width)
          p%uu(:,1)=vel_prof*(r1_mn**2)*(sin1th(m))*(&
              2*sin(theta-theta1)*cos(theta-theta1)*cos(theta)&
              -sin(theta)*sin(theta-theta1)**2)*&
              (x(l1:l2)-1.)*(x(l1:l2)-rone)**2 + &
              wind_amp*wind_prof
          p%uu(:,2)=-vel_prof*r1_mn*sin1th(m)*(&
              cos(theta)*sin(theta-theta1)**2)*&
              (x(l1:l2)-rone)*(3*x(l1:l2)-rone-2.)
          p%uu(:,3)=x(l1:l2)*sinth(m)*tanh(10.*(x(l1:l2)-x(l1)))
        endif
        if (lpenc_loc(i_divu)) p%divu=div_uprof*wind_amp + &
            div_vel_prof*(r1_mn**2)*(sin1th(m))*(&
            2*sin(theta-theta1)*cos(theta-theta1)*cos(theta)&
            -sin(theta)*sin(theta-theta1)**2)*&
            (x(l1:l2)-1.)*(x(l1:l2)-rone)**2
        if (lpenc_loc(i_der6u)) then
          p%der6u(:,1)=wind_amp*der6_uprof
          p%der6u(:,2)= 0.
          p%der6u(:,3)= 0.
        endif
!
!  KS-flow
!
      case ('KS')
        p%uu=0.
        do modeN=1,KS_modes  ! sum over KS_modes modes
          kdotxwt=KS_k(1,modeN)*x(l1:l2)+(KS_k(2,modeN)*y(m)+KS_k(3,modeN)*z(n))+KS_omega(modeN)*t
          cos_kdotxwt=cos(kdotxwt) ;  sin_kdotxwt=sin(kdotxwt)
          if (lpenc_loc(i_uu)) then
            p%uu(:,1) = p%uu(:,1) + cos_kdotxwt*KS_A(1,modeN) + sin_kdotxwt*KS_B(1,modeN)
            p%uu(:,2) = p%uu(:,2) + cos_kdotxwt*KS_A(2,modeN) + sin_kdotxwt*KS_B(2,modeN)
            p%uu(:,3) = p%uu(:,3) + cos_kdotxwt*KS_A(3,modeN) + sin_kdotxwt*KS_B(3,modeN)
          endif
        enddo
        if (lpenc_loc(i_divu))  p%divu = 0.
!
! flow from snapshot.
!
      case ('from-snap')
        if (lkinflow_as_aux) then
          if (lpenc_loc(i_uu)) p%uu=f(l1:l2,m,n,iux:iuz)
! divu
          if (lpenc_loc(i_divu)) p%divu=0. ! tb implemented
        else
          call inevitably_fatal_error('hydro_kinematic', '"from-snap" requires lkinflow_as_aux=T')
        endif
        lupdate_aux=.false.
!
      case ('Jouve-2008-benchmark-noav')
        if (lpenc_loc(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=ampl_kinflow*x(l1:l2)*sin(y(m))*( &
            -0.011125 + 0.5*(1.0 + erfunc((x(l1:l2)-0.7)/0.02)) &
            *(1.0-0.92-0.2*(cos(y(m)))**2))
        endif
!
! no kinematic flow.
!
      case ('none')
        if (headtt) print*,'kinflow = none'
        if (lpenc_loc(i_uu)) p%uu=0.
! divu
        if (lpenc_loc(i_divu)) p%divu=0.
      case('spher-harm-poloidal')
        if (lpenc_loc(i_uu)) p%uu=f(l1:l2,m,n,iux:iuz)
        lupdate_aux=.false.
      case('spher-harm-poloidal-per')
        if (lpenc_loc(i_uu)) &
          p%uu=f(l1:l2,m,n,iux:iuz)*cos(omega_kinflow*t)
        lupdate_aux=.false.
      case default
        call inevitably_fatal_error('hydro_kinematic', 'kinflow not found')
      end select
!
!  kinflows end here.
!  Note: in some parts above, p%oo and p%divu are already computed.
!  However, when lkinflow_as_comaux then ghost zones are computed
!  and p%uij, p%oo, and p%divu can be computed numerically.
!
      if (lkinflow_as_comaux) then
! uij
        if (lpenc_loc(i_uij)) call gij(f,iuu,p%uij,1)
! oo (=curlu)
        if (lpenc_loc(i_oo)) call curl_mn(p%uij,p%oo,p%uu)
! divu
        if (lpenc_loc(i_divu)) call div_mn(p%uij,p%divu,p%uu)
      endif
! u2
      if (lpenc_loc(i_u2)) call dot2_mn(p%uu,p%u2)
! o2
      if (lpenc_loc(i_o2)) call dot2_mn(p%oo,p%o2)
! ou
      if (lpenc_loc(i_ou)) call dot_mn(p%oo,p%uu,p%ou)
!
!  Calculate maxima and rms values for diagnostic purposes
!
      if (ldiagnos) then
        if (idiag_urms/=0)  call sum_mn_name(p%u2,idiag_urms,lsqrt=.true.)
        if (idiag_orms/=0)  call sum_mn_name(p%o2,idiag_orms,lsqrt=.true.)
        if (idiag_umax/=0)  call max_mn_name(p%u2,idiag_umax,lsqrt=.true.)
        if (idiag_omax/=0)  call max_mn_name(p%o2,idiag_omax,lsqrt=.true.)
        if (idiag_uzrms/=0) &
            call sum_mn_name(p%uu(:,3)**2,idiag_uzrms,lsqrt=.true.)
        if (idiag_uzmax/=0) &
            call max_mn_name(p%uu(:,3)**2,idiag_uzmax,lsqrt=.true.)
        if (idiag_u2m/=0)   call sum_mn_name(p%u2,idiag_u2m)
        if (idiag_um2/=0)   call max_mn_name(p%u2,idiag_um2)
        if (idiag_oum/=0)   call sum_mn_name(p%ou,idiag_oum)
!
        if (idiag_ekin/=0)  call sum_mn_name(.5*p%rho*p%u2,idiag_ekin)
        if (idiag_ekintot/=0) &
            call integrate_mn_name(.5*p%rho*p%u2,idiag_ekintot)
      endif
      if (idiag_divum/=0)  call sum_mn_name(p%divu,idiag_divum)
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_hydro_pencpar
!***********************************************************************
    subroutine hydro_before_boundary(f)
!
!  Dummy routine
!
!   16-dec-10/bing: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine hydro_before_boundary
!***********************************************************************
    subroutine duu_dt(f,df,p)
!
!  velocity evolution, dummy routine
!
!   7-jun-02/axel: adapted from hydro
!
      use Diagnostics
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: advec_uu
      logical, save :: lfirst_aux=.true.
!
      intent(in)  :: df,p
      intent(out) :: f
!
!  uu/dx for timestep (if kinflow is set)
!
      if (kinflow/='') then
        if (lfirst.and.ldt) then
          if (lspherical_coords) then
            advec_uu=abs(p%uu(:,1))*dx_1(l1:l2)+ &
                     abs(p%uu(:,2))*dy_1(  m  )*r1_mn+ &
                     abs(p%uu(:,3))*dz_1(  n  )*r1_mn*sin1th(m)
          elseif (lcylindrical_coords) then
            advec_uu=abs(p%uu(:,1))*dx_1(l1:l2)+ &
                     abs(p%uu(:,2))*dy_1(  m  )*rcyl_mn1+ &
                     abs(p%uu(:,3))*dz_1(  n  )
          else
            advec_uu=abs(p%uu(:,1))*dx_1(l1:l2)+ &
                     abs(p%uu(:,2))*dy_1(  m  )+ &
                     abs(p%uu(:,3))*dz_1(  n  )
          endif
          maxadvec=maxadvec+advec_uu
          if (headtt.or.ldebug) print*, 'duu_dt: max(advec_uu) =', maxval(advec_uu)
        endif
      endif
!
!  Store uu as auxiliary variable in f-array if requested by lkinflow_as_aux.
!  Just neccessary immediately before writing snapshots, but how would we
!  know we are?
!
     if (lpencil(i_uu).and.lkinflow_as_aux.and.(lupdate_aux.or.lfirst_aux)) f(l1:l2,m,n,iux:iuz)=p%uu
     if (.not.lpencil_check_at_work) lfirst_aux=.false.
!
!  Calculate maxima and rms values for diagnostic purposes.
!
      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'duu_dt: diagnostics ...'
!
!  Kinetic field components at one point (=pt).

        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_uxpt/=0) call save_name(p%uu(lpoint-nghost,1),idiag_uxpt)
          if (idiag_uypt/=0) call save_name(p%uu(lpoint-nghost,2),idiag_uypt)
          if (idiag_uzpt/=0) call save_name(p%uu(lpoint-nghost,3),idiag_uzpt)
          if (idiag_phase1/=0) call save_name(phase1,idiag_phase1)
          if (idiag_phase2/=0) call save_name(phase2,idiag_phase2)
        endif
      endif
!
      call keep_compiler_quiet(df)
!
    endsubroutine duu_dt
!***********************************************************************
    subroutine time_integrals_hydro(f,p)
!
!   1-jul-08/axel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine time_integrals_hydro
!***********************************************************************
   subroutine coriolis_cartesian(df,uu,velind)
!
!  Coriolis terms for cartesian geometry.
!
!  30-oct-09/MR: outsourced, parameter velind added
!  checked to be an equivalent change by auto-test conv-slab-noequi, mdwarf
!
      real, dimension (mx,my,mz,mvar), intent(out) :: df
      real, dimension (nx,3),          intent(in)  :: uu
      integer,                         intent(in)  :: velind
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(uu)
      call keep_compiler_quiet(velind)
!
   endsubroutine coriolis_cartesian
!***********************************************************************
    subroutine hydro_after_boundary(f)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine hydro_after_boundary
!***********************************************************************
    subroutine random_isotropic_KS_setup(initpower,kmin,kmax)
!
!  Produces random, isotropic field from energy spectrum following the
!  KS method (Malik and Vassilicos, 1999.)
!
!  More to do; unsatisfactory so far - at least for a steep power-law
!  energy spectrum.
!
!  27-may-05/tony: modified from snod's KS hydro initial
!  03-feb-06/weezy: Attempted rewrite to guarantee periodicity of
!                    KS modes.
!
      use Sub, only: cross, dot2
      use General, only: random_number_wrapper
!
      integer :: modeN
!
      real, dimension (3) :: k_unit
      real, dimension (3) :: ee,e1,e2
      real, dimension (6) :: r
      real :: initpower,kmin,kmax
      real, dimension(KS_modes) :: k,dk,energy,ps
      real :: theta,phi,alpha,beta
      real :: ex,ey,ez,norm,a
!
      allocate(KS_k(3,KS_modes))
      allocate(KS_A(3,KS_modes))
      allocate(KS_B(3,KS_modes))
      allocate(KS_omega(KS_modes))
!
      kmin=2.*pi      !/(1.0*Lxyz(1))
      kmax=128.*pi    !nx*pi
      a=(kmax/kmin)**(1./(KS_modes-1.))
!
!  Loop over all modes.
!
      do modeN=1,KS_modes
!
!  Pick wavenumber.
!
        k=kmin*(a**(modeN-1.))
!
!  Pick 4 random angles for each mode.
!
        call random_number_wrapper(r);
        theta=pi*(r(1) - 0.)
        phi=pi*(2*r(2) - 0.)
        alpha=pi*(2*r(3) - 0.)
        beta=pi*(2*r(4) - 0.)
!
!  Make a random unit vector by rotating fixed vector to random position
!  (alternatively make a random transformation matrix for each k).
!
        k_unit(1)=sin(theta)*cos(phi)
        k_unit(2)=sin(theta)*sin(phi)
        k_unit(3)=cos(theta)
!
        energy=(((k/kmin)**2. +1.)**(-11./6.))*(k**2.) &
            *exp(-0.5*(k/kmax)**2.)
!
!  Make a vector KS_k of length k from the unit vector for each mode.
!
        KS_k(:,modeN)=k*k_unit(:)
        KS_omega(:)=sqrt(energy(:)*(k(:)**3.))
!
!  Construct basis for plane having rr normal to it
!  (bit of code from forcing to construct x', y').
!
      if ((k_unit(2)==0).and.(k_unit(3)==0)) then
          ex=0.; ey=1.; ez=0.
        else
          ex=1.; ey=0.; ez=0.
        endif
        ee = (/ex, ey, ez/)
!
        call cross(k_unit(:),ee,e1)
!  e1: unit vector perp. to KS_k
        call dot2(e1,norm); e1=e1/sqrt(norm)
        call cross(k_unit(:),e1,e2)
!  e2: unit vector perp. to KS_k, e1
        call dot2(e2,norm); e2=e2/sqrt(norm)
!
!  Make two random unit vectors KS_B and KS_A in the constructed plane.
!
        KS_A(:,modeN) = cos(alpha)*e1 + sin(alpha)*e2
        KS_B(:,modeN) = cos(beta)*e1  + sin(beta)*e2
!
!  Make sure dk is set.
!
        call error('random_isotropic_KS_setup', 'Using uninitialized dk')
        dk=0.                     ! to make compiler happy
!
        ps=sqrt(2.*energy*dk)   !/3.0)
!
!  Give KS_A and KS_B length ps.
!
        KS_A(:,modeN)=ps*KS_A(:,modeN)
        KS_B(:,modeN)=ps*KS_B(:,modeN)
!
      enddo
!
!  Form RA = RA x k_unit and RB = RB x k_unit.
!  Note: cannot reuse same vector for input and output.
!
      do modeN=1,KS_modes
        call cross(KS_A(:,modeN),k_unit(:),KS_A(:,modeN))
        call cross(KS_B(:,modeN),k_unit(:),KS_B(:,modeN))
      enddo
!
      call keep_compiler_quiet(initpower)
!
    endsubroutine random_isotropic_KS_setup
!***********************************************************************
    subroutine random_isotropic_KS_setup_test
!
!  Produces random, isotropic field from energy spectrum following the
!  KS method (Malik and Vassilicos, 1999.)
!  This test case only uses 3 very specific modes (useful for comparison
!  with Louise's kinematic dynamo code.
!
!  03-feb-06/weezy: modified from random_isotropic_KS_setup
!
      use Sub, only: cross
      use General, only: random_number_wrapper
!
      integer :: modeN
!
      real, dimension (3,KS_modes) :: k_unit
      real, dimension(KS_modes) :: k,dk,energy,ps
      real :: initpower,kmin,kmax
!
      allocate(KS_k(3,KS_modes))
      allocate(KS_A(3,KS_modes))
      allocate(KS_B(3,KS_modes))
      allocate(KS_omega(KS_modes))
!
      initpower=-5./3.
      kmin=10.88279619
      kmax=23.50952672
!
      KS_k(1,1)=2.00*pi
      KS_k(2,1)=-2.00*pi
      KS_k(3,1)=2.00*pi
!
      KS_k(1,2)=-4.00*pi
      KS_k(2,2)=0.00*pi
      KS_k(3,2)=2.00*pi
!
      KS_k(1,3)=4.00*pi
      KS_k(2,3)=2.00*pi
      KS_k(3,3)=-6.00*pi
!
      KS_k(1,1)=+1; KS_k(2,1)=-1; KS_k(3,1)=1
      KS_k(1,2)=+0; KS_k(2,2)=-2; KS_k(3,2)=1
      KS_k(1,3)=+0; KS_k(2,3)=-0; KS_k(3,3)=1
!
      k(1)=kmin
      k(2)=14.04962946
      k(3)=kmax
!
      do modeN=1,KS_modes
        k_unit(:,modeN)=KS_k(:,modeN)/k(modeN)
      enddo
!
      kmax=k(KS_modes)
      kmin=k(1)
!
      do modeN=1,KS_modes
        if (modeN==1) dk(modeN)=(k(modeN+1)-k(modeN))/2.
        if (modeN>1.and.modeN<KS_modes) &
            dk(modeN)=(k(modeN+1)-k(modeN-1))/2.
        if (modeN==KS_modes) dk(modeN)=(k(modeN)-k(modeN-1))/2.
      enddo
!
      do modeN=1,KS_modes
         energy(modeN)=((k(modeN)**2 +1.)**(-11./6.))*(k(modeN)**2) &
             *exp(-0.5*(k(modeN)/kmax)**2)
      enddo
!
      ps=sqrt(2.*energy*dk)
!
      KS_A(1,1)=1.00/sqrt(2.00)
      KS_A(2,1)=-1.00/sqrt(2.00)
      KS_A(3,1)=0.00
!
      KS_A(1,2)=1.00/sqrt(3.00)
      KS_A(2,2)=1.00/sqrt(3.00)
      KS_A(3,2)=-1.00/sqrt(3.00)
!
      KS_A(1,3)=-1.00/2.00
      KS_A(2,3)=-1.00/2.00
      KS_A(3,3)=1.00/sqrt(2.00)
!
      KS_B(1,3)=1.00/sqrt(2.00)
      KS_B(2,3)=-1.00/sqrt(2.00)
      KS_B(3,3)=0.00
!
      KS_B(1,1)=1.00/sqrt(3.00)
      KS_B(2,1)=1.00/sqrt(3.00)
      KS_B(3,1)=-1.00/sqrt(3.00)
!
      KS_B(1,2)=-1.00/2.00
      KS_B(2,2)=-1.00/2.00
      KS_B(3,2)=1.00/sqrt(2.00)
!
      do modeN=1,KS_modes
        KS_A(:,modeN)=ps(modeN)*KS_A(:,modeN)
        KS_B(:,modeN)=ps(modeN)*KS_B(:,modeN)
      enddo
!
!  Form RA = RA x k_unit and RB = RB x k_unit.
!
       do modeN=1,KS_modes
         call cross(KS_A(:,modeN),k_unit(:,modeN),KS_A(:,modeN))
         call cross(KS_B(:,modeN),k_unit(:,modeN),KS_B(:,modeN))
       enddo
!
    endsubroutine random_isotropic_KS_setup_test
!***********************************************************************
    subroutine periodic_KS_setup(initpower)
!
!  Produces random, isotropic field from energy spectrum following the
!  KS method, however this setup produces periodic velocity field
!  (assuming box at least (-pi,pi)).
!
!  22-jun-2010/abag coded
!
      use Sub
      use General
!
      real, intent(in) :: initpower

      real, dimension(3,KS_modes) :: unit_k,k,A,B,orderK
      real, dimension(KS_modes) :: kk,delk,energy,omega,klengths
      real, dimension(KS_modes) :: ampA, ampB
      real, dimension(3) :: angle,dir_in
      real :: k_option(3,10000),mkunit(10000)
      real :: arg, unity
      real :: bubble, max_box
      real :: j(3),l(3),newa(3),newa2(3)
      integer ::i,s1,num,direction(3)
      logical :: ne
!
!  For the calculation of the velocity field we require
!  KS_k, KS_A/KS_B, KS_omega.
!
      allocate(KS_k(3,KS_modes), KS_A(3,KS_modes),KS_B(3,KS_modes),KS_omega(KS_modes))
!
!  Dummy variable.
!
      num=1
!
!  Space wavenumbers out.
!
      bubble=1.
!
!  Needs adapting for 2D runs.
!
      max_box=min(nxgrid,nygrid,nzgrid)
!
      if (mod(max_box,4.)/=0) print*, 'warning will not be periodic'
!
      print*, 'calculating KS wavenumbers'
!
      do i=1,10000
        call random_number_wrapper(angle)
        if ((angle(1)-0.0 < epsilon(0.0)) .or. &
            (angle(2)-0.0 < epsilon(0.0)) .or. &
            (angle(3)-0.0 < epsilon(0.0))) then
          call random_number_wrapper(angle)
        endif
!
!  Need 4 meshpoints to resolve a wave.
!
        angle=floor((max_box/4)*angle)
        call random_number_wrapper(dir_in)
        direction=nint(dir_in)
        direction=2*direction -1  !positive or negative directions
!
!  A possible orientation provided we havent's already got this length.
!
        k_option(1,i)=direction(1)*2.*pi*angle(1)
        k_option(2,i)=direction(2)*2.*pi*angle(2)
        k_option(3,i)=direction(3)*2.*pi*angle(3)
!
        if (i==1) then
          k_option(1,i)=2.*pi
          k_option(2,i)=0.
          k_option(3,i)=0.
        endif
!
!  Find the length of the current k_option vector.
!
        mkunit(i)=sqrt((k_option(1,i)**2)+(k_option(2,i)**2)+(k_option(3,i)**2))
!
        if (i==1.and.mkunit(i)>0.)then
          k(:,num)=k_option(:,i)
          klengths(num)=mkunit(i)
        endif
!
!  Now we check that the current length is unique (hasn't come before).
!
        if (i>1.and.num<KS_modes)then
          do s1=i-1,1,-1
            if (mkunit(i)>0.0.and. &
                mkunit(i)<=(mkunit(s1)-bubble).or. &
                mkunit(i)>=(mkunit(s1)+bubble)) then
              ne=.true.
            else
              ne=.false.
              exit
            endif
!
!  If length of current k_option is new......
!
            if (s1==1.and.ne) then
              num=num+1
!
!  Load current k_option into k that we keep.
!
              k(:,num)=k_option(:,i)
!
!  Store the length also.
!
              klengths(num)=mkunit(i)
            endif
          enddo
        endif
        if (i==10000.and.num<KS_modes) print*,"Haven't got",KS_modes,"modes!!!!"
      enddo
!
!  The 1 means ascending order.
!
      call KS_order(klengths,KS_modes,1,kk)
!
      do i=1,KS_modes
        do s1=1,KS_modes
          if (kk(i)==klengths(s1))then
            orderK(:,i)=k(:,s1)
          endif
        enddo
      enddo
!
      k=orderK
      do i=1,KS_modes
        unit_k(:,i)=k(:,i)/kk(i)
      enddo
!
      do i=1,KS_modes
!
!  Now we find delk as defined in Malik & Vassilicos' paper.
!
        if (i==1) delk(i)=(kk(i+1)-kk(i))/2.0
        if (i==KS_modes) delk(i)=(kk(i)-kk(i-1))/2.0
        if (i>1.and.i<KS_modes) delk(i)=(kk(i+1)-kk(i-1))/2.0
      enddo
!
!  Now find A&B that are perpendicular to each of our N wave-vectors.
!
      do i=1,KS_modes
!
!  Define "energy" - here we want k^{initpower} in the inertial range.
!
        energy(i)=1.0+(kk(i))**2
        energy(i)=(kk(i)**2)*(energy(i)**((initpower-2.)/2.))
        energy(i)=energy(i)*exp(-0.5*(kk(i)/kk(KS_modes))**2)
!
!  Set the lengths of A& B as defined in Malik & Vassilicos.
!
        ampA(i)=sqrt(2.0*energy(i)*delk(i))
        ampB(i)=ampA(i)
!
        call random_number(newa)
        call random_number(newa2)
        newa=2.0*newa -1.0
        newa2=2.0*newa2 -1.0
        j=newa
        l=newa2
        j=j/(sqrt(sum(newa**2)))
        l=l/(sqrt(sum(newa2**2)))
!
!  Now take the vector product of k with j (->A) and with l (->B).
!
        A(1,i)=(j(2)*k(3,i))-(j(3)*k(2,i))
        A(2,i)=(j(3)*k(1,i))-(j(1)*k(3,i))
        A(3,i)=(j(1)*k(2,i))-(j(2)*k(1,i))
        unity=sqrt((A(1,i)**2)+(A(2,i)**2)+(A(3,i)**2))
        A(:,i)=A(:,i)/unity
        B(1,i)=(l(2)*k(3,i))-(l(3)*k(2,i))
        B(2,i)=(l(3)*k(1,i))-(l(1)*k(3,i))
        B(3,i)=(l(1)*k(2,i))-(l(2)*k(1,i))
        unity=sqrt((B(1,i)**2)+(B(2,i)**2)+(B(3,i)**2))
        B(:,i)=B(:,i)/unity
!
!  Now that we have our unit A's & B's we multiply them by the amplitudes
!  we defined earlier, to create the spectrum.
!
        A(:,i)=ampA(i)*A(:,i)
        B(:,i)=ampB(i)*B(:,i)
      enddo
!
      do i=1,KS_modes
!
!  These are used to define omega - the unsteadyness frequency (co-eff of t).
!
        arg=energy(i)*(kk(i)**3)
        if (arg>0.0)omega(i)=sqrt(arg)
        if (arg==0.0)omega(i)=0.0
      enddo
!
      do i=1,KS_modes
        call cross(A(:,i),unit_k(:,i),KS_A(:,i))
        call cross(B(:,i),unit_k(:,i),KS_B(:,i))
      enddo
!
      KS_omega(:)=omega(:)
      KS_k=k
!
    endsubroutine periodic_KS_setup
!***********************************************************************
    subroutine KS_order(ad_, i_N, i_ord, B)
!
!  Bubble sort algorithm.
!
!  22-feb-10/abag: coded
!
      integer, intent(in) :: i_N, i_ord
      real, intent(in)  :: ad_(i_N)
      real, intent(out) :: B(i_N)
      real :: c=1.
      integer :: n
!
      B = ad_
!
      do while(c>0.0)
        c = -1.
        do n = 1, i_N-1
          if (   (i_ord==1 .and. B(n)>B(n+1)) &
            .or. (i_ord==2 .and. B(n)<B(n+1)) ) then
            c        = B(n)
            B(n)   = B(n+1)
            B(n+1) = c
            c = 1.
          endif
        enddo
      enddo
!
    endsubroutine KS_order
!***********************************************************************
    subroutine input_persistent_hydro(id,done)
!
!  Read in the stored time of the next random phase calculation.
!
!  12-apr-08/axel: adapted from input_persistent_forcing
!
      use IO, only: read_persist
!
      integer :: id
      logical :: done
!
      if (id == id_record_HYDRO_TPHASE) then
        if (read_persist ('HYDRO_TPHASE', tphase_kinflow)) return
        done = .true.
      elseif (id == id_record_HYDRO_PHASE1) then
        if (read_persist ('HYDRO_PHASE1', phase1)) return
        done = .true.
      elseif (id == id_record_HYDRO_PHASE2) then
        if (read_persist ('HYDRO_PHASE2', phase2)) return
        done = .true.
      elseif (id == id_record_HYDRO_LOCATION) then
        if (read_persist ('HYDRO_LOCATION', location)) return
        done = .true.
      elseif (id == id_record_HYDRO_TSFORCE) then
        if (read_persist ('HYDRO_TSFORCE', tsforce)) return
        done = .true.
      elseif (id == id_record_HYDRO_AMPL) then
        if (read_persist ('HYDRO_AMPL', tsforce_ampl)) return
        done = .true.
      elseif (id == id_record_HYDRO_WAVENUMBER) then
        if (read_persist ('HYDRO_WAVENUMBER', tsforce_wavenumber)) return
        done = .true.
      endif
!
      if (lroot) print*,'input_persistent_hydro: ',tphase_kinflow
!
    endsubroutine input_persistent_hydro
!***********************************************************************
    logical function output_persistent_hydro()
!
!  Writes out the time of the next random phase calculation.
!
!  12-apr-08/axel: adapted from output_persistent_forcing
!  16-nov-11/MR: changed into logical function to signal I/O errors, I/O error handling introduced
!  01-Jun-2015/Bourdin.KIS: activated this code by inserting the call in input_/output_persistent
!
      use IO, only: write_persist
!
      output_persistent_hydro = .false.
      if (tphase_kinflow < 0.) return
      if (lroot .and. (ip < 14)) print *,'output_persistent_hydro: ', tphase_kinflow
!
!  Write details.
!
      output_persistent_hydro = .true.
!
      if (write_persist ('HYDRO_TPHASE', id_record_HYDRO_TPHASE, tphase_kinflow)) return
      if (write_persist ('HYDRO_PHASE1', id_record_HYDRO_PHASE1, phase1)) return
      if (write_persist ('HYDRO_PHASE2', id_record_HYDRO_PHASE2, phase2)) return
      if (write_persist ('HYDRO_LOCATION', id_record_HYDRO_LOCATION, location)) return
      if (write_persist ('HYDRO_TSFORCE', id_record_HYDRO_TSFORCE, tsforce)) return
      if (write_persist ('HYDRO_AMPL', id_record_HYDRO_AMPL, tsforce_ampl)) return
      if (write_persist ('HYDRO_WAVENUMBER', id_record_HYDRO_WAVENUMBER, tsforce_wavenumber)) return
!
      output_persistent_hydro = .false.
!
    endfunction output_persistent_hydro
!***********************************************************************
    subroutine read_hydro_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_hydro_init_pars
!***********************************************************************
    subroutine write_hydro_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_hydro_init_pars
!***********************************************************************
    subroutine read_hydro_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=hydro_run_pars, IOSTAT=iostat)
!
    endsubroutine read_hydro_run_pars
!***********************************************************************
    subroutine write_hydro_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=hydro_run_pars)
!
    endsubroutine write_hydro_run_pars
!***********************************************************************
    subroutine rprint_hydro(lreset,lwrite)
!
!  Reads and registers print parameters relevant for hydro part.
!
!   8-jun-02/axel: adapted from hydro
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
!  Reset everything in case of reset
!  (this needs to be consistent with what is defined above!).
!
      if (lreset) then
        idiag_u2m=0; idiag_um2=0; idiag_oum=0; idiag_o2m=0
        idiag_uxpt=0; idiag_uypt=0; idiag_uzpt=0; idiag_dtu=0
        idiag_urms=0; idiag_umax=0; idiag_uzrms=0; idiag_uzmax=0
        idiag_phase1=0; idiag_phase2=0
        idiag_orms=0; idiag_omax=0; idiag_oumphi=0
        idiag_ruxm=0; idiag_ruym=0; idiag_ruzm=0; idiag_rumax=0
        idiag_ux2m=0; idiag_uy2m=0; idiag_uz2m=0
        idiag_uxuym=0; idiag_uxuzm=0; idiag_uyuzm=0
        idiag_umx=0; idiag_umy=0; idiag_umz=0
        idiag_Marms=0; idiag_Mamax=0; idiag_divu2m=0; idiag_epsK=0
        idiag_urmphi=0; idiag_upmphi=0; idiag_uzmphi=0; idiag_u2mphi=0
        idiag_ekin=0; idiag_ekintot=0
        idiag_divum=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_hydro: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ekin',idiag_ekin)
        call parse_name(iname,cname(iname),cform(iname),'ekintot',idiag_ekintot)
        call parse_name(iname,cname(iname),cform(iname),'u2m',idiag_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',idiag_um2)
        call parse_name(iname,cname(iname),cform(iname),'o2m',idiag_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',idiag_oum)
        call parse_name(iname,cname(iname),cform(iname),'dtu',idiag_dtu)
        call parse_name(iname,cname(iname),cform(iname),'urms',idiag_urms)
        call parse_name(iname,cname(iname),cform(iname),'umax',idiag_umax)
        call parse_name(iname,cname(iname),cform(iname),'uzrms',idiag_uzrms)
        call parse_name(iname,cname(iname),cform(iname),'uzmax',idiag_uzmax)
        call parse_name(iname,cname(iname),cform(iname),'ux2m',idiag_ux2m)
        call parse_name(iname,cname(iname),cform(iname),'uy2m',idiag_uy2m)
        call parse_name(iname,cname(iname),cform(iname),'uz2m',idiag_uz2m)
        call parse_name(iname,cname(iname),cform(iname),'uxuym',idiag_uxuym)
        call parse_name(iname,cname(iname),cform(iname),'uxuzm',idiag_uxuzm)
        call parse_name(iname,cname(iname),cform(iname),'uyuzm',idiag_uyuzm)
        call parse_name(iname,cname(iname),cform(iname),'orms',idiag_orms)
        call parse_name(iname,cname(iname),cform(iname),'omax',idiag_omax)
        call parse_name(iname,cname(iname),cform(iname),'divum',idiag_divum)
        call parse_name(iname,cname(iname),cform(iname),'ruxm',idiag_ruxm)
        call parse_name(iname,cname(iname),cform(iname),'ruym',idiag_ruym)
        call parse_name(iname,cname(iname),cform(iname),'ruzm',idiag_ruzm)
        call parse_name(iname,cname(iname),cform(iname),'rumax',idiag_rumax)
        call parse_name(iname,cname(iname),cform(iname),'umx',idiag_umx)
        call parse_name(iname,cname(iname),cform(iname),'umy',idiag_umy)
        call parse_name(iname,cname(iname),cform(iname),'umz',idiag_umz)
        call parse_name(iname,cname(iname),cform(iname),'Marms',idiag_Marms)
        call parse_name(iname,cname(iname),cform(iname),'Mamax',idiag_Mamax)
        call parse_name(iname,cname(iname),cform(iname),'divu2m',idiag_divu2m)
        call parse_name(iname,cname(iname),cform(iname),'epsK',idiag_epsK)
        call parse_name(iname,cname(iname),cform(iname),'uxpt',idiag_uxpt)
        call parse_name(iname,cname(iname),cform(iname),'uypt',idiag_uypt)
        call parse_name(iname,cname(iname),cform(iname),'uzpt',idiag_uzpt)
        call parse_name(iname,cname(iname),cform(iname),'phase1',idiag_phase1)
        call parse_name(iname,cname(iname),cform(iname),'phase2',idiag_phase2)
      enddo
!
!  Write column where which variable is stored.
!
      if (lwr) then
        call farray_index_append('iuu',iuu)
        call farray_index_append('iux',iux)
        call farray_index_append('iuy',iuy)
        call farray_index_append('iuz',iuz)
      endif
!
    endsubroutine rprint_hydro
!***********************************************************************
    subroutine get_slices_hydro(f,slices)
!
!  Write slices for animation of Hydro variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_hydro
!***********************************************************************
    subroutine hydro_after_timestep(f,df,dt_sub)
!
!  Hook for modification of the f and df arrays
!  according to the hydro module, after the       
!  timestep is performed. 
!
!  12-mar-17/wlyra: coded. 
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dt_sub
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dt_sub)

    endsubroutine hydro_after_timestep
!***********************************************************************
    subroutine calc_mflow
!
!  Dummy routine.
!
!  19-jul-03/axel: adapted from hydro
!
    endsubroutine calc_mflow
!***********************************************************************
    subroutine remove_mean_momenta(f)
!
!  Dummy routine.
!
!  32-nov-06/tobi: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine remove_mean_momenta
!***********************************************************************
    subroutine remove_mean_flow(f,indux)
!
!  Dummy.
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      integer,                            intent (in)    :: indux

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(indux)
!
    endsubroutine remove_mean_flow
!***********************************************************************
    subroutine impose_velocity_ceiling(f)
!
!  13-aug-2007/anders: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_velocity_ceiling
!***********************************************************************
    subroutine init_ck
!
!  8-sep-2009/dhruba: coded
!
      integer :: l,m
      real :: Balpha,jl,jlp1,jlm1,LPl,LPlm1
      integer :: ell
!
      print*, 'Initializing variables for Chandrasekhar-Kendall flow'
      if (ldebug) print*, 'Allocating..'
      allocate(Zl(mx),dZldr(mx))
      allocate(Pl(my),dPldtheta(my))
      if (ldebug) print*, 'Allocation done'
      ell=kinflow_ck_ell
      Balpha=kinflow_ck_Balpha
      print*, 'ell=,alpha=',ell,Balpha
!
      do l=1,mx
        call sp_besselj_l(jl,ell,Balpha*x(l))
        call sp_besselj_l(jlp1,ell+1,Balpha*x(l))
        call sp_besselj_l(jlm1,ell-1,Balpha*x(l))
        Zl(l) = jl
        dZldr(l) = ell*jlm1-(ell+1)*jlp1
      enddo
!
      do m=1,my
        call legendre_pl(LPl,ell,y(m))
        call legendre_pl(LPlm1,ell-1,y(m))
        Pl(m) = Lpl
        dPldtheta(m) = -(1/sin(y(m)))*ell*(LPlm1-LPl)
      enddo
!
    endsubroutine init_ck
!***********************************************************************
    subroutine hydro_clean_up
!
!  Deallocate the variables allocated in nohydro.
!
!  8-sep-2009/dhruba: coded
!
      if (ldebug) print*, 'Deallocating some nohydro variables ...'
      if (kinflow=='ck') then
        deallocate(Zl,dZldr)
        deallocate(Pl,dPldtheta)
      elseif (kinflow=='KS') then
        deallocate(KS_k)
        deallocate(KS_A)
        deallocate(KS_B)
        deallocate(KS_omega)
      endif
!
      print*, 'Done.'
!
    endsubroutine hydro_clean_up
!***********************************************************************
    subroutine kinematic_random_phase
!
!  Get a random phase to be used for the whole kinematic velocity field.
!
!  16-feb-2010/dhruba: coded
!
      use General, only: random_number_wrapper
!
      real, dimension(3) :: fran
!
!  Generate random numbers.
!
      if (t>tsforce) then
        if (lrandom_location) then
          call random_number_wrapper(fran)
          location=fran*Lxyz+xyz0
        else
          location=location_fixed
        endif
!
!  Writing down the location.
!
        if (lroot .and. lwrite_random_location) then
          open(1,file=trim(datadir)//'/random_location.dat',status='unknown',position='append')
          write(1,'(4f14.7)') t, location
          close (1)
        endif
!
!  Update next tsforce.
!
        tsforce=t+dtforce
        if (ip<=6) print*,'kinematic_random_phase: location=',location
      endif
!
    endsubroutine kinematic_random_phase
!***********************************************************************
    subroutine kinematic_random_ampl
!
!  Get a random amplitude to be used for the whole kinematic velocity field.
!
!  21-jun-2019/axel: coded
!
      use General, only: random_number_wrapper
!
      real, dimension(2) :: fran
      real :: random_r, random_p
!
!  Generate random numbers.
!
      if (t>=tsforce_ampl) then
        if (lrandom_ampl) then
          call random_number_wrapper(random_r)
          call random_number_wrapper(random_p)
          random_ampl=sqrt(-2*log(random_r))*sin(2*pi*random_p)
        else
          random_ampl=1.
        endif
!
!  Writing down the location.
!
        if (lroot .and. lwrite_random_ampl) then
          open(1,file=trim(datadir)//'/random_ampl.dat',status='unknown',position='append')
          write(1,'(4f14.7)') t, random_ampl
          close (1)
        endif
!
!  Update next tsforce_ampl.
!
        tsforce_ampl=t+dtforce
        if (ip<=6) print*,'kinematic_random_phase: location=',location
      endif
!
    endsubroutine kinematic_random_ampl
!***********************************************************************
    subroutine kinematic_random_wavenumber
!
!  Get a random wavenumber to be used for the whole kinematic velocity field.
!
!  21-jun-2019/axel: coded
!
      use General, only: random_number_wrapper
!
      real :: fran
!
!  Generate random numbers.
!
      if (t>=tsforce_wavenumber) then
        if (lrandom_wavenumber) then
          call random_number_wrapper(fran)
          random_wavenumber=nint(fran*.5*twopi/Lxyz(1)*nxgrid)+1.
        else
          random_wavenumber=1.
        endif
!
!  Writing down the location.
!
        if (lroot .and. lwrite_random_wavenumber) then
          open(1,file=trim(datadir)//'/random_wavenumber.dat',status='unknown',position='append')
          write(1,'(4f14.7)') t, random_wavenumber
          close (1)
        endif
!
!  Update next tsforce_wavenumber.
!
        tsforce_wavenumber=t+dtforce
        if (ip<=6) print*,'kinematic_random_phase: location=',location
      endif
!
    endsubroutine kinematic_random_wavenumber
!***********************************************************************
    subroutine expand_shands_hydro
!
!  Presently dummy, for later use
!
    endsubroutine expand_shands_hydro
!***********************************************************************
    subroutine density_after_timestep(f,df,dtsub)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dtsub
!
    endsubroutine density_after_timestep
!***********************************************************************
    subroutine calc_gradu(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
  
      call fatal_error('calc_gradu','not implemented in hydro_kinematic') 
 
    endsubroutine calc_gradu
!***********************************************************************    
    subroutine update_char_vel_hydro(f)
!
!   25-sep-15/MR+joern: for slope limited diffusion
!
!   calculation of characteristic velocity
!   for slope limited diffusion
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, parameter :: i64_1=1/64.
!
      if (lslope_limit_diff) then
        if (lkinflow_as_aux) then
           f( :mx-1, :my-1, :mz-1,iFF_diff2)=f( :mx-1, :my-1, :mz-1,iFF_diff2) &
           +i64_1 *sum((f( :mx-1, :my-1, :mz-1,iux:iuz) &
                       +f( :mx-1, :my-1,2:mz,  iux:iuz) &
                       +f( :mx-1,2:my ,  :mz-1,iux:iuz) &
                       +f( :mx-1,2:my , 2:mz  ,iux:iuz) &
                       +f(2:mx  , :my-1, :mz-1,iux:iuz) &
                       +f(2:mx  , :my-1,2:mz  ,iux:iuz) &
                       +f(2:mx  ,2:my  , :mz-1,iux:iuz) &
                       +f(2:mx  ,2:my  ,2:mz  ,iux:iuz))**2,4)
        else
          call warning('update_char_vel_hydro','characteristic velocity not implemented for hydro_kinematic')
        endif
      endif

    endsubroutine update_char_vel_hydro
!***********************************************************************
endmodule Hydro
