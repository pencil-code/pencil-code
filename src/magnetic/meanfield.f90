! $Id$
!
!  This modules solves mean-field contributions to both the
!  induction and the momentum equations.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagn_mf = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED mf_EMF(3); mf_EMFdotB; jxb_mf(3); jxbr_mf(3); chiB_mf
! PENCILS PROVIDED mf_qp; mf_Beq21
!
!***************************************************************
module Magnetic_meanfield
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: fatal_error,inevitably_fatal_error,svn_id
  use Magnetic_meanfield_demfdt
!
  implicit none
!
  include 'meanfield.h'
!
!  array for inputting alpha profile
!
  real, dimension (mx,my) :: alpha_input
  real, pointer :: B_ext2, eta
  logical, pointer :: lweyl_gauge
!
  real, dimension (nx,3,3) :: hij
  real, dimension (nx,3,3,3) :: hijk
  real, dimension (nx) :: kf_x, kf_x1
  real, dimension (my) :: kf_y
  real, dimension (mz) :: kf_z
  real, dimension (nz) :: fpatches=0, xpatches=0, ypatches=0
!
!  arrays for the fluctuating alpha effect in axisymmetric spherical coordinates
!
  integer :: channel_magn_mf=1, channel_magn_mf2=2, nr_cell=1, nth_cell=1
  integer, dimension (nx) :: ir_cell=1
  integer, dimension (nx,my) :: itheta_cell=1
  real, dimension (nx) :: tau_cell=1.0
  real, dimension (nx,my) :: rand_alpha=0.
!
! Parameters
!
  character (len=labellen) :: meanfield_kf_profile='const'
  character (len=labellen) :: meanfield_etat_profile='const'
  character (len=labellen) :: meanfield_Beq_profile='const',qp_model='atan'
  character (len=labellen) :: Omega_profile='nothing', alpha_profile='const'
  character (len=labellen) :: EMF_profile='nothing', delta_profile='const'
  character (len=labellen) :: shear_current_profile='nothing',fluc_alp_profile='gaussian'
!
! Input parameters
!
  real :: Omega_ampl=0.0, dummy=0.0
  real :: Calp=0.0, alpha_effect=0.0, alpha_quenching=0.0, delta_effect=0.0, alpha_zz=0.
  real :: gamma_effect=0.0, gamma_quenching=0.0
  real :: chit_quenching=0.0, chi_t0=0.0
  real :: meanfield_etat=0.0, meanfield_etat_height=1., meanfield_pumping=1.
  real :: meanfield_delta_height=10., meanfield_delta_width=.05
  real :: meanfield_Beq=1.0,meanfield_Beq_height=0., meanfield_Beq2_height=0.
  real :: meanfield_etat_exp=1.0, uturb=.1
  real :: alpha_eps=0.0, alpha_pom0=0.0
  real :: x_surface=0., x_surface2=0., z_surface=0., qp_width=impossible, qpx_width=impossible
  real :: alpha_equator=impossible, alpha_equator_gap=0.0, alpha_gap_step=0.0
  real :: alpha_rmin=0.
  real :: alpha_cutoff_up=0.0, alpha_cutoff_down=0.0
  real :: meanfield_qs=0.0, meanfield_qp=0.0, meanfield_qe=0.0, meanfield_qa=0.0
  real :: meanfield_Bs=1.0, meanfield_Bp=1.0, meanfield_Be=1.0, meanfield_Ba=1.0
  real :: meanfield_qp1, meanfield_qs1, meanfield_qe1, meanfield_qa1
  real :: meanfield_etaB=0.0
  real :: x1_alp=0., x2_alp=0., y1_alp=0., y2_alp=0.
  logical :: lOmega_effect=.false., lalpha_Omega_approx=.false.
  logical :: lmeanfield_noalpm=.false., lmeanfield_pumping=.false.
  logical :: lmeanfield_jxb=.false., lmeanfield_jxb_with_vA2=.false.
  logical :: lmeanfield_chitB=.false., lignore_gradB2_inchiB=.false.
  logical :: lchit_with_glnTT=.false., lrho_chit=.true., lchit_Bext2_equil=.false.
  logical :: lturb_temp_diff=.false., lqp_profile=.false., lqpx_profile=.false.
  logical :: lGW_tensor=.false.
!
  namelist /magn_mf_init_pars/ &
      x1_alp, x2_alp, y1_alp, y2_alp
!
! Run parameters
!
  real :: alpha_rmax=0.0, alpha_width=0.0, alpha_width2=0.0, alpha_exp=0.
  real :: meanfield_etat_width=0.0, meanfield_etat_corona=0.0,meanfield_Beq_width=0.0
  real :: meanfield_kf=1.0, meanfield_kf_width=0.0, meanfield_kf_width2=0.0
  real :: meanfield_cs0=1.0
  real :: Omega_rmax=0.0, Omega_rwidth=0.0
  real :: rhs_term_kx=0.0, rhs_term_ampl=0.0
  real :: rhs_term_amplz=0.0, rhs_term_amplphi=0.0
  real :: mf_qJ2=0.0, qp_aniso_factor=1.0
  real :: kx_alpha=1., kx_hij=1., relhel_hij=1., hij_ampl=1.
  real :: GWfac1=1., GWfac2=1., GWfac3=1.
  real :: fluc_alp_m=1.0, sigma_alpha=1.0
  real :: b2_to_u2=0.0, shear_current_sh=0.0
  real :: sigx=0.0, sigz=0.0
  integer :: npatches=1, npatches_actual, seed_magn_mf2=5555
  real, dimension(3) :: alpha_aniso=0.
  real, dimension(3,3) :: alpha_tensor=0., eta_tensor=0.
  real, dimension(ny,3,3) :: alpha_tensor_y=0., eta_tensor_y=0.
  real, dimension(nz,2,2) :: alpha_tensor_z=0., eta_tensor_z=0.
  real, dimension(nx) :: rhs_termz, rhs_termy
  real, dimension(nx) :: etat_x, detat_x, rhs_term
  real, dimension(my) :: etat_y, detat_y
  real, dimension(mz) :: etat_z, detat_z, qp_profile, qp_profder !, qe_profder
  real, dimension(nx) :: qpx_profile, qpx_profder
  logical :: llarge_scale_velocity=.false.
  logical :: lEMF_profile=.false.
  logical :: lalpha_profile_total=.false., lalpha_aniso=.false.
  logical :: ldelta_profile=.false., lalpha_tensor=.false., leta_tensor=.false.
  logical :: lrhs_term=.false., lrhs_term2=.false.
  logical :: lqpcurrent=.false., lNEMPI_correction=.true.
  logical :: lNEMPI_correction_qp_profile=.true.
  logical :: lqp_aniso_factor=.false.
  logical :: lread_alpha_tensor_z=.false., lread_alpha_tensor_z_as_y=.false.
  logical :: lread_eta_tensor_z=.false., lread_eta_tensor_z_as_y=.false.
  logical :: lshear_current_effect=.false., lalphass_disk=.false.
  logical :: ltest_patches=.false., lOmega_effect_meanfield=.false.
  real :: ampluu_kinematic=0.
  real, dimension(:), allocatable :: xcenter, ycenter, zcenter, roty, cy, sy
!
  namelist /magn_mf_run_pars/ &
      Calp, alpha_effect, alpha_quenching, alpha_rmax, alpha_exp, alpha_zz, &
      gamma_effect, gamma_quenching, &
      alpha_eps, alpha_pom0, alpha_width, alpha_width2, alpha_aniso, &
      alpha_tensor, eta_tensor, &
      lalpha_profile_total, lmeanfield_noalpm, alpha_profile, &
      chit_quenching, chi_t0, lqp_profile, lqpx_profile, qp_width, qpx_width, &
      x_surface, x_surface2, z_surface, &
      alpha_rmin, kx_alpha, &
      qp_model, seed_magn_mf2, &
      npatches, sigx, sigz, ltest_patches, lOmega_effect_meanfield, ampluu_kinematic, &
      ldelta_profile, delta_effect, delta_profile, &
      meanfield_etat, meanfield_etat_height, meanfield_etat_profile, &
      meanfield_etat_width, meanfield_etat_exp, meanfield_etat_corona, meanfield_Beq_width, &
      meanfield_kf, meanfield_kf_profile, &
      meanfield_kf_width, meanfield_kf_width2, &
      meanfield_Beq, meanfield_Beq_height, meanfield_Beq2_height, &
      meanfield_Beq_profile, uturb, &
      meanfield_delta_height, meanfield_delta_width, &
      meanfield_cs0, &
      lmeanfield_pumping, meanfield_pumping, &
      lmeanfield_jxb, lmeanfield_jxb_with_vA2, &
      lmeanfield_chitB, lchit_with_glnTT, lrho_chit, lchit_Bext2_equil, &
      lturb_temp_diff, lignore_gradB2_inchiB, &
      meanfield_qs, meanfield_qp, meanfield_qe, meanfield_qa, &
      meanfield_Bs, meanfield_Bp, meanfield_Be, meanfield_Ba, &
      lqpcurrent,mf_qJ2, lNEMPI_correction, lNEMPI_correction_qp_profile, &
      lqp_aniso_factor, qp_aniso_factor, &
      alpha_equator, alpha_equator_gap, alpha_gap_step, &
      alpha_cutoff_up, alpha_cutoff_down, &
      lalpha_Omega_approx, lOmega_effect, Omega_profile, Omega_ampl, &
      llarge_scale_velocity, EMF_profile, lEMF_profile, &
      fluc_alp_profile, fluc_alp_m, sigma_alpha, &
      lshear_current_effect, shear_current_profile, b2_to_u2, shear_current_sh, &
      lrhs_term, lrhs_term2, rhs_term_amplz, rhs_term_amplphi, rhs_term_ampl, &
      Omega_rmax, Omega_rwidth, lread_alpha_tensor_z, lread_eta_tensor_z, &
      lread_alpha_tensor_z_as_y, lread_eta_tensor_z_as_y, &
      x1_alp, x2_alp, y1_alp, y2_alp, &
      lGW_tensor, kx_hij, relhel_hij, hij_ampl, GWfac1, GWfac2, GWfac3
!
! Diagnostic variables (need to be consistent with reset list below)
!
  integer :: idiag_qsm=0        ! DIAG_DOC: $\left<q_p(\overline{B})\right>$
  integer :: idiag_qpm=0        ! DIAG_DOC: $\left<q_p(\overline{B})\right>$
  integer :: idiag_qem=0        ! DIAG_DOC: $\left<q_e(\overline{B})\right>$,
                                ! DIAG_DOC: in the paper referred to as
                                ! DIAG_DOC: $\left<q_g(\overline{B})\right>$
  integer :: idiag_qam=0        ! DIAG_DOC: $\left<q_a(\overline{B})\right>$
  integer :: idiag_alpm=0       ! DIAG_DOC: $\left<\alpha\right>$  !(where is this implemented?)
  integer :: idiag_etatm=0      ! DIAG_DOC: $\left<\eta_{\rm t}\right>$
  integer :: idiag_EMFmz1=0     ! DIAG_DOC: $\left<{\cal E}\right>_{xy}|_x$
  integer :: idiag_EMFmz2=0     ! DIAG_DOC: $\left<{\cal E}\right>_{xy}|_y$
  integer :: idiag_EMFmz3=0     ! DIAG_DOC: $\left<{\cal E}\right>_{xy}|_z$
  integer :: idiag_EMFdotBm=0   ! DIAG_DOC: $\left<{\cal E}\cdot\Bv \right>$
  integer :: idiag_EMFdotB_int=0! DIAG_DOC: $\int{\cal E}\cdot\Bv dV$
  integer :: idiag_alpKjbm=0    ! DIAG_DOC: $\left<\alpha_\mathrm{K}\overline{\Bv}\cdot\overline{\Jv}\right>$
  integer :: idiag_alpKm=0      ! DIAG_DOC: $\left<\alpha_\mathrm{K}\right>$
  integer :: idiag_peffmxz=0    ! YAVG_DOC: $\left<{\cal P}_{\rm eff}\right>_{y}$
  integer :: idiag_alpmxz=0     ! YAVG_DOC: $\left<\alpha\right>_{y}$
!
! xy averaged diagnostics given in xyaver.in
!
  integer :: idiag_qpmz=0       ! XYAVG_DOC: $\left<q_p\right>_{xy}$
!
  contains
!***********************************************************************
    subroutine register_magn_mf()
!
!  No additional variables to be initialized here, but other mean-field
!  modules that are being called from here may need additional variables.
!
!  6-jan-11/axel: adapted from magnetic
!
!  Identify version number.
!
      use SharedVariables, only: put_shared_variable

      if (lroot) call svn_id( &
          "$Id$")
!
!  Register secondary mean-field modules.
!
      if (lmagn_mf_demfdt) call register_magn_mf_demfdt()
!
!  if meanfield theory is invoked, we want to send meanfield_etat to
!  other subroutines
!
      if (lrun) then
        call put_shared_variable('meanfield_etat',meanfield_etat,caller='initialize_magn_mf')
!
!  Quenching parameters: lmeanfield_chitB and chi_t0.
!
        call put_shared_variable('lmeanfield_chitB',lmeanfield_chitB)
        if (lmeanfield_chitB) call put_shared_variable('chi_t0',chi_t0)
!
!  chit_quenching for flux boundary condition
!
        !dummy=meanfield_Beq_profile
        !call put_shared_variable('meanfield_Beq_profile',dummy)
        call put_shared_variable('meanfield_Beq',meanfield_Beq)
        call put_shared_variable('chit_quenching',chit_quenching)
        call put_shared_variable('uturb',uturb)
!
      endif
!
!  share etat profile with viscosity module
!
      if (meanfield_etat/=0.0.and.lviscosity) then
        call put_shared_variable('etat_x',etat_x)
        call put_shared_variable('etat_y',etat_y)
        call put_shared_variable('etat_z',etat_z)
        call put_shared_variable('detat_x',detat_x)
        call put_shared_variable('detat_y',detat_y)
        call put_shared_variable('detat_z',detat_z)
      endif
!
!  Initialize secondary mean-field modules:
!
      if (lmagn_mf_demfdt .or. lalpm .or. lalpm_alternate ) then
        call put_shared_variable('kf_x',kf_x)
        call put_shared_variable('kf_y',kf_y)
        call put_shared_variable('kf_z',kf_z)
        call put_shared_variable('kf_x1',kf_x1)
        call put_shared_variable('etat_x',etat_x)
        call put_shared_variable('etat_y',etat_y)
        call put_shared_variable('etat_z',etat_z)
      endif
!
    endsubroutine register_magn_mf
!***********************************************************************
    subroutine initialize_magn_mf(f)
!
!  Perform any post-parameter-read initialization
!
!  20-may-03/axel: reinitialize_aa added
!
      use Sub, only: erfunc
      use General, only: random_number_wrapper, random_seed_wrapper, random_gen
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: kf_x_tmp, kf_x1_tmp, prof_tmp
      real, dimension (nz) :: rpatches
      real :: kz1, cont_count=0., rpatch_norm
      integer :: ierr, i, j, ipatchz_count, ipatchz, npatchz
      logical :: file_exist
!
!  Get B_ext2 and eta from magnetic module.
!
      call get_shared_variable('B_ext2',B_ext2)
      call get_shared_variable('eta',eta)
!
!  Possibility of giving Calp=alpha_effect/[(eta+etat)*k1]
!
       if (alpha_effect==0. .and. Calp/=0.) then
         kz1=2.*pi/Lxyz(3)
         alpha_effect=Calp*(eta+meanfield_etat)*kz1
         if (lroot) print*,'alpha_effect=',alpha_effect
       endif
!
!  check for alpha profile
!
      if (alpha_profile=='read') then
        print*,'read alpha profile'
        open(1,file='alpha_input.dat',form='unformatted')
        read(1) alpha_input
        close(1)
      endif
!
!  check for (simple-minded) anisotropy
!
      if (any(alpha_aniso/=0.)) then
        lalpha_aniso=.true.
      else
        lalpha_aniso=.false.
      endif
      if (lroot) print*,'lalpha_aniso=',lalpha_aniso
!
!  check for (slightly less-simple minded) anisotropy
!  Put alpha_zz=1 to get alpha_ij = delta_ij - OiOj.
!
      if (any(alpha_tensor/=0.) .or. alpha_zz/=0.) then
        lalpha_tensor=.true.
        if (alpha_zz /= 0.) then
          alpha_tensor(1,1)=1.
          alpha_tensor(2,2)=1.
          alpha_tensor(3,3)=1.
          alpha_tensor_y(:,1,1)=-alpha_zz*cos(y(m1:m2))**2
          alpha_tensor_y(:,2,2)=-alpha_zz*sin(y(m1:m2))**2
          alpha_tensor_y(:,1,2)=+alpha_zz*sin(y(m1:m2))*cos(y(m1:m2))
          alpha_tensor_y(:,2,1)=alpha_tensor_y(:,1,2)
        elseif (lread_alpha_tensor_z.or.lread_alpha_tensor_z_as_y) then
!
!  Read in alpha tensor (z-dependent, and only 2x2)
!
          open(1,file='alpha_tensor_z.dat',form='unformatted')
          read(1) alpha_tensor_z
          close(1)
!
!  Each component can be scaled separately.
!  Conversely, they must all be turned on separately.
!
          do i=1,2; do j=1,2
            alpha_tensor_z(:,i,j)=alpha_tensor(i,j)*alpha_tensor_z(:,i,j)
          enddo; enddo
!
!  if lread_alpha_tensor_z_as_y: swap back (x,y,z) <-- (y,z,x)
!
          if (lread_alpha_tensor_z_as_y) then
            if (ny==nz) then
              alpha_tensor_y(n1:n2,3,3)=alpha_tensor_z(:,1,1)
              alpha_tensor_y(n1:n2,3,1)=alpha_tensor_z(:,1,2)
              alpha_tensor_y(n1:n2,1,3)=alpha_tensor_z(:,2,1)
              alpha_tensor_y(n1:n2,1,1)=alpha_tensor_z(:,2,2)
            else
              call fatal_error('initialize_magn_mf', &
                'lread_alpha_tensor_z_as_y works only when ny=nz')
            endif
          endif
        endif
      else
        lalpha_tensor=.false.
      endif
      if (lroot) print*,'lalpha_tensor=',lalpha_tensor
!
!  check for (slightly less-simple minded) anisotropy
!  Put alpha_zz=1 to get alpha_ij = delta_ij - OiOj.
!
      if (any(eta_tensor/=0.)) then
        leta_tensor=.true.
        if (lread_eta_tensor_z.or.lread_eta_tensor_z_as_y) then
!
!  Read in eta tensor (z-dependent, and only 2x2)
!
          open(1,file='eta_tensor_z.dat',form='unformatted')
          read(1) eta_tensor_z
          close(1)
!
!  Each component can be scaled separately.
!  Conversely, they must all be turned on separately.
!
          do i=1,2; do j=1,2
            eta_tensor_z(:,i,j)=eta_tensor(i,j)*eta_tensor_z(:,i,j)
          enddo; enddo
!
!  if lread_alpha_tensor_z_as_y: swap back (x,y,z) <-- (y,z,x)
!
          if (lread_eta_tensor_z_as_y) then
            if (ny==nz) then
              eta_tensor_y(n1:n2,3,3)=eta_tensor_z(:,1,1)
              eta_tensor_y(n1:n2,3,1)=eta_tensor_z(:,1,2)
              eta_tensor_y(n1:n2,1,3)=eta_tensor_z(:,2,1)
              eta_tensor_y(n1:n2,1,1)=eta_tensor_z(:,2,2)
            else
              call fatal_error('initialize_magn_mf', &
                'lread_eta_tensor_z_as_y works only when ny=nz')
            endif
          endif
        endif
      else
        leta_tensor=.false.
      endif
      if (lroot) print*,'leta_tensor=',leta_tensor
!
!  write profile (uncomment for debugging)
!
!     if (lroot) then
!       do n=n1,n2
!         print*,z(n),eta_z(n)
!       enddo
!     endif
!
!  Possibility of adding a coskx term to the rhs of the dAz/dt equation.
!
      if (lrhs_term) then
        !rhs_term=(rhs_term_ampl/rhs_term_kx)*cos(rhs_term_kx*x(l1:l2))
        !rhs_term=rhs_term_ampl*(.25*x(l1:l2)**2-.0625*x(l1:l2)**4)
        rhs_term=rhs_term_ampl*(1.-x(l1:l2)**2)
        !rhs_term=rhs_term_ampl
      endif
!
      if (lrhs_term2) then
        !rhs_termz=rhs_term_amplz*(-x(l1:l2)**2)/4.      !kk
        !rhs_termy=rhs_term_amplphi*x(l1:l2)/2.          !kk
        rhs_termy=.50*rhs_term_amplphi*(1.-x(l1:l2))*x(l1:l2)
        rhs_termz=.25*rhs_term_amplz*(.75-(1.-.25*x(l1:l2)**2)*x(l1:l2)**2)
      endif
!
!  Compute etat profile and share with other routines.
!  Here we also set the etat_i and getat_i profiles.
!
      if (meanfield_kf/=0.0) then
        select case (meanfield_kf_profile)
        case ('const')
          kf_x1=1./meanfield_kf
          kf_x=meanfield_kf
          kf_y=1.
          kf_z=1.
        case ('sqrt_surface_x1')
          prof_tmp=.5*(1.+erfunc((x(l1:l2)-x_surface2)/meanfield_kf_width2))
          kf_x1_tmp=max(1.-x(l1:l2)/x_surface,tini)/meanfield_kf
          kf_x1=prof_tmp*kf_x1_tmp
          kf_x_tmp=1./kf_x1_tmp
          kf_x=prof_tmp*kf_x_tmp*exp(-(meanfield_kf_width*kf_x_tmp)**2)
          kf_y=1.
          kf_z=1.
        case default;
          call inevitably_fatal_error('initialize_magnetic', &
          'no such meanfield_kf_profile profile')
        endselect
      endif
!
!  Compute etat profile and share with other routines.
!  Here we also set the etat_i and getat_i profiles.
!
      if (meanfield_etat/=0.0) then
        select case (meanfield_etat_profile)
        case ('const')
          etat_x=1.
          etat_y=1.
          etat_z=meanfield_etat
          detat_x=0.
          detat_y=0.
          detat_z=0.
        case ('exp(z/H)')
          etat_x=1.
          etat_y=meanfield_etat
          etat_z=exp(z/meanfield_etat_height)
          detat_x=0.
          detat_y=0.
          detat_z=etat_z/meanfield_etat_height
        case ('surface_x1'); etat_x=0.5 &
          *(1.+erfunc((x(l1:l2)-x_surface2)/meanfield_etat_width))
          etat_y=meanfield_etat
          etat_z=1.
          detat_x=1./(meanfield_etat_width*sqrtpi) &
            *exp(-((x(l1:l2)-x_surface2)/meanfield_etat_width)**2)
          detat_y=0.
          detat_z=0.
        case ('surface_x2'); etat_x=0.25 &
          *(1.-erfunc((x(l1:l2)-x_surface)/meanfield_etat_width)) &
          *(1.+erfunc((x(l1:l2)-x_surface2)/meanfield_etat_width))
          etat_y=meanfield_etat
          etat_z=1.
          detat_x=.5/(meanfield_etat_width*sqrtpi)*( &
             exp(-((x(l1:l2)-x_surface)/meanfield_etat_width)**2) &
            *(1.+erfunc((x(l1:l2)-x_surface2)/meanfield_etat_width)) &
            +(1.+erfunc((x(l1:l2)-x_surface)/meanfield_etat_width)) &
            *exp(-((x(l1:l2)-x_surface2)/meanfield_etat_width)**2))
          detat_y=0.
          detat_z=0.
        case ('sqrt_surface_x1')
          etat_x=sqrt(kf_x1)
          etat_y=meanfield_etat
          etat_z=1.
          if (ip<=6) print*,'etat(sqrt_surface_x1)=',etat_x
          detat_x=1./(meanfield_etat_width*sqrtpi) &
            *exp(-((x(l1:l2)-x_surface2)/meanfield_etat_width)**2)
          detat_y=0.
          detat_z=0.
        case ('sin2y')
          etat_x=meanfield_etat
          etat_y=sin(y)**2
          etat_z=1.
          detat_x=0.
          detat_y=2.*sin(y)*cos(y)
          detat_z=0.
        case ('sin4y')
          etat_x=meanfield_etat
          etat_y=sin(y)**4
          etat_z=1.
          detat_x=0.
          detat_y=4.*sin(y)**3*cos(y)
          detat_z=0.
!
!  general expression, includes previous cases with meanfield_etat_exp=2 and 4
!
        case ('siny**n')
          etat_x=meanfield_etat
          etat_y=sin(y)**meanfield_etat_exp
          etat_z=1.
          detat_x=0.
          detat_y=meanfield_etat_exp*sin(y)**(meanfield_etat_exp-1.)*cos(y)
          detat_z=0.
        case ('Jouve-2008-benchmark')
          etat_x = meanfield_etat*(0.01 + 0.5*(1.-0.01)*(1.0+erfunc((x(l1:l2)-0.7)/0.02)))
          etat_y = 1.
          etat_z = 1.
          detat_x= meanfield_etat*0.5*(1.-0.01)*exp(-((x(l1:l2)-0.7)/0.02)**2)
          detat_y= 0.
          detat_z= 0.
        case ('sinx**2*exp(-z/H)')
          etat_x=sin(x(l1:l2))**2
          etat_y=meanfield_etat
          etat_z=exp(-z/meanfield_etat_height)
          detat_x=2.*sin(x(l1:l2))*cos(x(l1:l2))
          detat_y=0.
          detat_z=-etat_z/meanfield_etat_height
        case ('surface_z'); etat_z=0.5 &
          *(1.-erfunc((z-z_surface)/meanfield_etat_width))
          etat_y=meanfield_etat
          etat_x=1.
          detat_z=-1./(meanfield_etat_width*sqrtpi) &
            *exp(-((z-z_surface)/meanfield_etat_width)**2)
          detat_y=0.
          detat_x=0.
        case ('alphass')
!
!  etat=meanfield_etat*cs2/Omega in spherical coordinates.
!  Here the omega part only. The cs2 part will be added later
!
          etat_x = x(l1:l2)**1.5
          etat_y = max(0.,sin(y))**1.5
          etat_z = meanfield_etat
          detat_x= 0.
          detat_y= 0.
          detat_z= 0.
        case default;
          call inevitably_fatal_error('initialize_magnetic', &
          'no such meanfield_etat_profile: '//trim(meanfield_etat_profile))
        endselect
!
        if (lviscosity) then
          print*,'ipz,z(n),etat_z(n),detat_z(n)'
          do n=n1,n2
            print*,ipz,z(n),etat_z(n),detat_z(n)
          enddo
          print*
        endif
      endif
!
!  define inverse of meanfield_qp
!
      if (meanfield_qp==0.) then
        meanfield_qp1=0.
      else
        meanfield_qp1=1./meanfield_qp
      endif
!
!  define inverse of meanfield_qs
!
      if (meanfield_qs==0.) then
        meanfield_qs1=0.
      else
        meanfield_qs1=1./meanfield_qs
      endif
!
!  define inverse of meanfield_qe
!
      if (meanfield_qe==0.) then
        meanfield_qe1=0.
      else
        meanfield_qe1=1./meanfield_qe
      endif
!
!  define inverse of meanfield_qa
!
      if (meanfield_qa==0.) then
        meanfield_qa1=0.
      else
        meanfield_qa1=1./meanfield_qa
      endif
!
!  define meanfield_qp_profile
!
      if (lqp_profile) then
        qp_profile=0.5*(1.-erfunc((z-z_surface)/qp_width))
        !qe_profile=0.5*(1.-erfunc((z-z_surface)/qe_width))
        qp_profder=-exp(-((z-z_surface)/qp_width)**2)/(qp_width*sqrtpi)
        !qe_profder=-exp(-((z-z_surface)/qe_width)**2)/(qe_width*sqrtpi)
      else
        qp_profile=1.
        qp_profder=0.
      endif
!
!  define meanfield_qpx_profile (x direction)
!
      if (lqpx_profile) then
        qpx_profile=exp(-(x(l1:l2)/qpx_width)**2)
        qpx_profder=-2.*x(l1:l2)/qpx_width**2*exp(-(x(l1:l2)/qpx_width)**2)
      else
        qpx_profile=1.
        qpx_profder=0.
      endif
!
!  Initialize module variables which are parameter dependent
!  wave speed of gauge potential
!
      if (lrun) then
        call get_shared_variable('lweyl_gauge',lweyl_gauge)
        if (lroot) print*,'initialize_magn_mf: lweyl_gauge=',lweyl_gauge
!       if (.not.lweyl_gauge) call get_shared_variable('eta',eta)
      endif
!
      if (lmagn_mf_demfdt .or. lalpm .or. lalpm_alternate) call initialize_magn_mf_demfdt(f)
!
!  compute GW tensor
!
      hij=0.
      hijk=0.
      if (lGW_tensor) then
        hij(:,2,2)=hij_ampl*cos(kx_hij*x(l1:l2))
        hij(:,2,3)=hij_ampl*sin(kx_hij*x(l1:l2))*relhel_hij
        hij(:,3,3)=-hij(:,2,2)
        hij(:,3,2)=+hij(:,2,3)
        hijk(:,2,2,1)=-kx_hij*hij_ampl*sin(kx_hij*x(l1:l2))
        hijk(:,2,3,1)=+kx_hij*hij_ampl*cos(kx_hij*x(l1:l2))*relhel_hij
        hijk(:,3,3,1)=-hijk(:,2,2,1)
        hijk(:,3,2,1)=+hijk(:,2,3,1)
      endif
!
!  thin disk model switch
!
      if (alpha_profile=='alphass' .or. alpha_profile=='fluc-alpha-disk' &
          .or. meanfield_etat_profile=='alphass'.or. &
          shear_current_profile=='alphass') lalphass_disk=.true.
!
!  Determine the coherent cells for fluctuating alpha dynamo
!
      if (alpha_profile=='fluc-alpha-disk') call calc_fluc_alp_cells
!
!  Patches allows one to place small dynamos in the domain. Their distribution
!  may be uniform in space if the file npatchz.txt does not exist.
!
      if (alpha_profile=='patches') then
        inquire(file='npatchz.txt',exist=file_exist)
        if (file_exist) then
          open(1,file='npatchz.txt')
          do n=1,nz
            read(1,*) i,fpatches(n),xpatches(n)
          enddo
          close(1)
        endif
!
        allocate(xcenter(npatches))
        allocate(ycenter(npatches))
        allocate(zcenter(npatches))
        allocate(roty(npatches))
        allocate(cy(npatches))
        allocate(sy(npatches))
!
!  Use channel=channel_magn_mf2=2 for all these patches.
!
        if (random_gen/='nr_f90') call fatal_error('initialize_magn_mf', &
          'must use random_gen=nr_f90, because we need channel=2!')
        seed2=0; seed2(1)=seed_magn_mf2
        call random_seed_wrapper(PUT=seed2,CHANNEL=channel_magn_mf2)
        call random_number_wrapper(xcenter,CHANNEL=channel_magn_mf2)
        call random_number_wrapper(ycenter,CHANNEL=channel_magn_mf2)
        call random_number_wrapper(zcenter,CHANNEL=channel_magn_mf2)
        call random_number_wrapper(rpatches,CHANNEL=channel_magn_mf2)
!
!  Assume full x-range if the file npatchz.txt does not exist.
!
        if (file_exist) then
          ipatchz_count=0
          do n=1,nz
!
!  number of patches per z value
!
            rpatch_norm=nz/sum(rpatches)
            npatchz=nint(rpatch_norm*rpatches(n)*fpatches(n)*npatches)
            cont_count=cont_count+2.*rpatches(n)*fpatches(n)*npatches
!
!  Allow for multiple patches; loop over all of them.
!  Print a line of there are too many for each patch (for now).
!  The actual number of patches (npatches_actual) may be smaller than npatches.
!
            do ipatchz=1,npatchz
              ipatchz_count=ipatchz_count+1
              if (ipatchz_count<=npatches) then
                xcenter(ipatchz_count)=-xpatches(n)+2.*xpatches(n)*xcenter(ipatchz_count)
                ycenter(ipatchz_count)=-ypatches(n)+2.*ypatches(n)*ycenter(ipatchz_count)
                zcenter(ipatchz_count)=z(n+nghost)
                print*,n,z(n+nghost),npatchz,ipatchz_count,cont_count,2.*rpatches(n)*fpatches(n)*npatches
              else
                print*,'too many patches: ipatchz_count,cont_count=',ipatchz_count,cont_count
              endif
            enddo
          enddo
          npatches_actual=ipatchz_count
        else
          xcenter=xyz0(1)+(xyz1(1)-xyz0(1))*xcenter 
          ycenter=xyz0(2)+(xyz1(2)-xyz0(2))*ycenter 
          zcenter=xyz0(3)+(xyz1(3)-xyz0(3))*zcenter 
          npatches_actual=npatches
        endif
!
!  Assume random orientation angles for all patches.
!
        call random_number_wrapper(roty,CHANNEL=channel_magn_mf2)
        cy=cos(360.*roty*dtor)
        sy=sin(360.*roty*dtor)
      endif
    endsubroutine initialize_magn_mf
!***********************************************************************
    subroutine init_aa_mf(f)
!
!  Initialise mean-field related magnetic field; called from magnetic.f90
!  At the moment, no own initial conditions are allowed, but we need this
!  to call secondary modules
!
!   6-jan-2011/axel: adapted from magnetic
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Initialize secondary mean-field modules.
!
      if (lmagn_mf_demfdt) call init_aa_mf_demfdt(f)
!
    endsubroutine init_aa_mf
!***********************************************************************
    subroutine pencil_criteria_magn_mf()
!
!   All pencils that the magnetic mean-field module depends on
!   are specified here.
!
!  28-jul-10/axel: adapted from magnetic
!
      lpenc_requested(i_bb)=.true.
!
      if (meanfield_etat/=0.0.or.ietat/=0) &
          lpenc_requested(i_del2a)=.true.
!
!  In mean-field theory, with variable etat, need divA for resistive gauge.
!
      if (meanfield_etat_profile/='const') then
        lpenc_requested(i_diva)=.true.
      endif
!
!  In spherical coordinates, we also need grad(divA)
!
      if (lspherical_coords) then
        lpenc_requested(i_graddiva)=.true.
        lpenc_requested(i_x_mn)=lOmega_effect
      endif
!
!  For mean-field modelling in momentum equation:
!
      if (lmeanfield_jxb) then
        lpenc_requested(i_b2)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_jxbr_mf)=.true.
        if (lmeanfield_jxb_with_vA2) lpenc_requested(i_va2)=.true.
      endif
!
!  For mean-field modelling in entropy equation:
!
      if (lmeanfield_chitB) then
        lpenc_requested(i_b2)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_chiB_mf)=.true.
        if (lrho_chit) lpenc_requested(i_rho1)=.true.
        if (lturb_temp_diff) then
          lpenc_requested(i_cp)=.true.
          lpenc_requested(i_glnTT)=.true.
          lpenc_requested(i_del2lnTT)=.true.
        else
          lpenc_requested(i_gss)=.true.
          lpenc_requested(i_del2ss)=.true.
        endif
      endif
!
!  Turbulent diffusivity
!
      if (meanfield_etat/=0.0 .or. ietat/=0 .or. &
          alpha_effect/=0.0 .or. delta_effect/=0.0 .or. &
          gamma_effect/=0.0 .or. lread_alpha_tensor_z .or. &
          lread_eta_tensor_z .or. lread_alpha_tensor_z_as_y .or. &
          lread_eta_tensor_z_as_y) &
          lpenc_requested(i_mf_EMF)=.true.
      if (delta_effect/=0.0) lpenc_requested(i_oxJ)=.true.
!
!  J is needed when the eta tensor is invoked.
!  (But there should be other such requests that have apparently
!  not yet caused pencil checks to fail yet...)
!
      if (leta_tensor) then
        lpenc_requested(i_jj)=.true.
      endif
!
      if (idiag_EMFdotBm/=0.or.idiag_EMFdotB_int/=0) lpenc_diagnos(i_mf_EMFdotB)=.true.
!
! If large_scale_velocity is included we need this pencil
!
      if (llarge_scale_velocity) lpenc_requested(i_uxb)=.true.
!
!  Pencil criteria for secondary modules
!
      if (lmagn_mf_demfdt) call pencil_criteria_magn_mf_demfdt()
!
!  Nempi model C requires currents
!
      if (lqpcurrent) then
        lpenc_requested(i_jj)=.true.
        lpenc_requested(i_j2)=.true.
        lpenc_requested(i_jij)=.true.
      endif
!
      if (lNEMPI_correction) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_b2)=.true.
      endif
!
      if (lGW_tensor) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_mf_EMF)=.true.
      endif
!
!  If idiag_alpKjbm/=0, we need j.b.
!
      if (idiag_alpKjbm/=0) lpenc_requested(i_jb)=.true.
!
!  thin-disk dynamo models
!
      if (lalphass_disk) lpenc_requested(i_cs2)=.true.
!
      if (meanfield_Beq_profile=='alphass'.or. &
          meanfield_Beq_profile=='alphass-algebraic'.or. &
          meanfield_Beq_profile=='fluc-alpha-disk'.or.&
          shear_current_profile=='alphass') lpenc_requested(i_pp)=.true.
!
    endsubroutine pencil_criteria_magn_mf
!***********************************************************************
    subroutine pencil_interdep_magn_mf(lpencil_in)
!
!  Interdependency among pencils from the Magnetic module is specified here.
!
!  28-jul-10/axel: adapted from magnetic
!
      logical, dimension(npencils) :: lpencil_in
!
!  exa
!
      if (lpencil_in(i_exa)) then
        lpencil_in(i_aa)=.true.
        lpencil_in(i_mf_EMF)=.true.
      endif
!
!  mf_EMFdotB
!
      if (lpencil_in(i_mf_EMFdotB)) then
        lpencil_in(i_mf_EMF)=.true.
        lpencil_in(i_bb)=.true.
      endif
!
!  mf_EMFdotB
!
      if (lpencil_in(i_mf_EMF)) then
        if (lspherical_coords) then
          lpencil_in(i_jj)=.true.
          lpencil_in(i_graddivA)=.true.
        endif
        lpencil_in(i_b2)=.true.
!
!  compute del2a regardless of gauge
!
        if (meanfield_etat/=0.0 .or. ietat/=0) then
!         if (lweyl_gauge) then
!           lpencil_in(i_jj)=.true.
!         else
            lpencil_in(i_del2a)=.true.
!         endif
        endif
      endif
!
!  oxJ effect
!
        if (lpencil_in(i_oxJ)) lpencil_in(i_jj)=.true.
!
!  ??
!
!     if (lpencil_in(i_del2A)) then
!       if (lspherical_coords) then
!       endif
!     endif
!
!  Mean-field Lorentz force: jxb_mf
!
      if (lpencil_in(i_jxbr_mf)) lpencil_in(i_jxb_mf)=.true.
      if (lpencil_in(i_jxb_mf)) lpencil_in(i_jxb)=.true.
!
!  Quenching prescription using plasma beta
!
      if (meanfield_Beq_profile=='alphass') lpenc_requested(i_b2)=.true.
!
!  Shear-current effect
!
      if (lshear_current_effect) lpencil_in(i_jj)=.true.
!
      if (shear_current_profile=='alphass') lpencil_in(i_b2)=.true.
!
!  Pencil criteria for secondary modules
!
      if (lmagn_mf_demfdt) call pencil_interdep_magn_mf_demfdt(lpencil_in)
!
    endsubroutine pencil_interdep_magn_mf
!***********************************************************************
    subroutine calc_pencils_magn_mf(f,p)
!
!  Calculate Magnetic mean-field pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  28-jul-10/axel: adapted from magnetic
!
      use Sub
      use General, only: bessj
      use Diagnostics, only: sum_mn_name, xysum_mn_name_z, ysum_mn_name_xz
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (nx) :: alpha_total, rr, pp, pp0, spiral
      real, dimension (nx) :: meanfield_etat_tmp
      real, dimension (nx) :: alpha_tmp, alpha_quenching_tmp, delta_tmp
      real, dimension (nx) :: kf_tmp, EMF_prof, alpm, prefact, Beq21
      real, dimension (nx) :: meanfield_qs_func, meanfield_qp_func, meanfield_qa_func
      real, dimension (nx) :: meanfield_qe_func, meanfield_qs_der, meanfield_qa_der
      real, dimension (nx) :: meanfield_qp_der, meanfield_qe_der, BiBk_Bki
      real, dimension (nx) :: meanfield_Bs21, meanfield_Bp21, meanfield_Be21, meanfield_Ba21
      real, dimension (nx) :: meanfield_etaB2, g2, chit_prof !, quench_chiB
      real, dimension (nx) :: oneQbeta02, oneQbeta2, XXj_BBj, BjBzj
      real, dimension (nx) :: shear_current_sh_tmp, disk_height, z_over_h
      real, dimension (nx,3) :: Bk_Bki, exa_meanfield, glnchit_prof, glnchit, XXj
      real, dimension (nx,3) :: meanfield_getat_tmp, getat_cross_B_tmp, B2glnrho, glnchit2
      real, dimension (nx) :: r2, x1, y1, z1, x11, y11, z11
      real :: kx,fact
      integer :: i, j, k, nn, l, ipatch
!
      intent(inout) :: f,p
!
      save :: chit_prof, glnchit_prof
!
!  mean-field Lorentz force
!
      if (lpencil(i_jxbr_mf)) then
!
!  The following 9 lines have not been used for any publications so far.
!
!        if (lmeanfield_jxb_with_vA2) then
!          call inevitably_fatal_error('calc_pencils_magnetic', &
!            'lmeanfield_jxb_with_vA2 should be checked')
!--
!         meanfield_urms21=1./(3.*meanfield_kf*meanfield_etat)**2
!         meanfield_qs_func=meanfield_qs*(1.-2*pi_1*atan(p%vA2*meanfield_urms21))
!         meanfield_qp_func=meanfield_qp*(1.-2*pi_1*atan(p%vA2*meanfield_urms21))
!         meanfield_qe_func=meanfield_qe*(1.-2*pi_1*atan(p%b2*meanfield_Be21))
!         meanfield_qs_der=2*pi_1*meanfield_qs/(1.+(p%vA2*meanfield_urms21)**2)
!         meanfield_qp_der=2*pi_1*meanfield_qp/(1.+(p%vA2*meanfield_urms21)**2)
!         meanfield_qe_der=2*pi_1*meanfield_qe*meanfield_Be21/(1.+(p%b2*meanfield_Be21)**2)
!         call multsv_mn(meanfield_qs_func,p%jxb,tmp_jxb); p%jxb_mf=tmp_jxb
!--
!
!  The following (not lmeanfield_jxb_with_vA2) has been used for the
!  various publications so far.
!
!        else
          select case (meanfield_Beq_profile)
          case ('exp(z/H)');
!           H = -2* density scale height, old version
            Beq21=1./meanfield_Beq**2*exp(-2*z(n)/meanfield_Beq_height)
!
!  Beq ~ exp(-z/2H), where H=Hrho is the density scale height.
!  Give here meanfield_Beq2_height, which is the scale height for Beq^2.
!
          case ('exp(z/2H)');
!           H = density scale height, more straightforward
            Beq21=1./meanfield_Beq**2*exp(z(n)/meanfield_Beq2_height)
!
!  Beq profile for constant uturb
!
          case ('uturbconst');
!
!  allows to take into account a shifted equilibrium solution
!  caused by the additional pressure contribution
!
            Beq21=mu01*p%rho1/(uturb**2)
!
!  Beq21 becomes large in the corona
!
          case ('uturb_surface');
            Beq21=mu01*p%rho1/(uturb**2)* &
              0.5*(1.-erfunc((z(n)-z_surface)/meanfield_Beq_width))
!
!  thin disk models; meanfield_Beq=alphaSS
!
          case ('alphass');
            Beq21=1.0/max(1e-10, 3.*meanfield_Beq**2*mu0*p%pp-p%b2)
          case ('alphass-algebraic')
            Beq21=mu01/( 3*meanfield_Beq**2*(1.-meanfield_Beq)*p%pp )
          case ('fluc-alpha-disk')
            Beq21=(sqrt(fluc_alp_m)*sigma_alpha-1.)*4./meanfield_Beq/mu0/p%pp
!
!  default
!
          case default;
            Beq21=1./meanfield_Beq**2
          endselect
          p%mf_Beq21=Beq21
!
!  Here, p%b2*meanfield_Bp21 = beta^2/beta_p^2, etc.
!
          meanfield_Bs21=Beq21/meanfield_Bs**2
          meanfield_Bp21=Beq21/meanfield_Bp**2
          meanfield_Be21=Beq21/meanfield_Be**2
          meanfield_Ba21=Beq21/meanfield_Ba**2
!
!  Compute qp(B^2), qs(B^2), qe(B^2), and their derivatives for 2 different parameterizations
!  qp=qp0/(1+B^2*meanfield_Bp21)
!  qp'=-qp0/(1+B^2*meanfield_Bp21)^2*meanfield_Bp21=-qp^2*meanfield_Bp21/qp0
!
          if(qp_model=='rational') then
            meanfield_qp_func=meanfield_qp/(1.+p%b2*meanfield_Bp21)*qp_profile(n)*qpx_profile
            meanfield_qs_func=meanfield_qs/(1.+p%b2*meanfield_Bs21)
            meanfield_qe_func=meanfield_qe/(1.+p%b2*meanfield_Be21)
            meanfield_qa_func=meanfield_qa/(1.+p%b2*meanfield_Ba21)
            meanfield_qp_der=-meanfield_qp_func**2*meanfield_Bp21*meanfield_qp1
            meanfield_qs_der=-meanfield_qs_func**2*meanfield_Bs21*meanfield_qs1
            meanfield_qe_der=-meanfield_qe_func**2*meanfield_Be21*meanfield_qe1
            meanfield_qa_der=-meanfield_qa_func**2*meanfield_Ba21*meanfield_qa1
          else
            meanfield_qp_func=meanfield_qp*(1.-2*pi_1*atan(p%b2*meanfield_Bp21))*qp_profile(n)
            meanfield_qs_func=meanfield_qs*(1.-2*pi_1*atan(p%b2*meanfield_Bs21))
            meanfield_qe_func=meanfield_qe*(1.-2*pi_1*atan(p%b2*meanfield_Be21))
            meanfield_qp_der=-2*pi_1*meanfield_qp*meanfield_Bp21/(1.+(p%b2*meanfield_Bp21)**2)
            meanfield_qs_der=-2*pi_1*meanfield_qs*meanfield_Bs21/(1.+(p%b2*meanfield_Bs21)**2)
            meanfield_qe_der=-2*pi_1*meanfield_qe*meanfield_Be21/(1.+(p%b2*meanfield_Be21)**2)
          endif
          p%mf_qp=meanfield_qp_func
!
!  Add (1/2)*grad[qp*B^2]. This initializes p%jxb_mf.
!  Note: p%jij is not J_i,j; omit extra term for the time being.
!  Also: "p%b2*meanfield_qp_der" is the same as dqp/dlnbeta^2.
!
          call multmv_transp(p%bij,p%bb,Bk_Bki) !=1/2 grad B^2
          if (lqpcurrent) then
            call multsv_mn(mu01*(1+p%j2*mf_qJ2)* &
              (meanfield_qp_func+p%b2*meanfield_qp_der), &
              Bk_Bki,p%jxb_mf)
          else
            call multsv_mn(mu01*(meanfield_qp_func+p%b2*meanfield_qp_der),Bk_Bki,p%jxb_mf)
            if (lNEMPI_correction) then
              call multsv_mn_add(-.5*mu01*p%b2**2*meanfield_qp_der,p%glnrho,p%jxb_mf)
            endif
!
!  Allow for qp profile
!
            if (lNEMPI_correction_qp_profile) then
              p%jxb_mf(:,3)=p%jxb_mf(:,3)+.5*mu01*qp_profder(n)*meanfield_qp_func*p%b2
            endif
          endif
!
!  Allow for simple-minded anisotropy
!
          if (lqp_aniso_factor) p%jxb_mf(:,3)=p%jxb_mf(:,3)*qp_aniso_factor
!
!  Add -B.grad[qs*B_i]. This term does not promote instability.
!
          call multsv_mn_add(-meanfield_qs_func,p%jxb+mu01*Bk_Bki,p%jxb_mf)
          call dot(Bk_Bki,p%bb,BiBk_Bki)
          call multsv_mn_add(-2.*mu01*meanfield_qs_der*BiBk_Bki,p%bb,p%jxb_mf)
!
!  Add e_z*grad(qe*B^2). This has not yet been found to promote instability.
!  Added below (in commented out form, how I think it should be if used)
!
          p%jxb_mf(:,3)=p%jxb_mf(:,3)+2.*mu01*(meanfield_qe_der*p%b2+meanfield_qe_func)*Bk_Bki(:,3)
          !p%jxb_mf(:,3)=p%jxb_mf(:,3)+mu01*(meanfield_qe_der*p%b2+meanfield_qe_func)*Bk_Bki(:,3) &
          !                        -.5*mu01*p%b2**2*meanfield_qe_der,p%glnrho(:,3) &
          !                        +.5*mu01*qe_profder(n)*meanfield_qe_func*p%b2
!
!  Add div[qa*Bz*(Bi*gj+Bj*gi)]
!  Begin by initializing auxiliary vector X_j (see notes)
!
          XXj=2.*Bk_Bki
          call multsv_mn_add(-mu01*p%rho1*p%b2,p%grho,XXj)
!
          if (qp_model=='rational') then
!
!  Do dqa/dB2 * Bz * (BiX3+z3B.X) term
!
            call dot(XXj,p%bb,XXj_BBj)
            p%jxb_mf(:,3)=p%jxb_mf(:,3)+meanfield_qa_der*p%bb(:,3)*XXj_BBj
            call multsv_mn_add(meanfield_qa_der*p%bb(:,3)*XXj(:,3),p%bb,p%jxb_mf)
!
!  Do qa*[Bz,z*(Bi+Bj*Bz,j*zi)] term
!
            call dot(p%bb,p%bij(:,3,:),BjBzj)
            p%jxb_mf(:,3)=p%jxb_mf(:,3)+mu01*meanfield_qa_func*p%bij(:,3,3)*BjBzj
            call multsv_mn_add(mu01*meanfield_qa_func*p%bij(:,3,3),p%bb,p%jxb_mf)
!
!  Do qa*Bz*B_i,z
!
            call multsv_mn_add(mu01*meanfield_qa_func*p%bb(:,3),p%bij(:,:,3),p%jxb_mf)
!
          endif
!
        call multsv_mn(p%rho1,p%jxb_mf,p%jxbr_mf)
      endif
!
!  Needed quantities for a thin disk model
!
      if (lalphass_disk) then
        disk_height = sqrt(p%cs2)/sqrt(1.0)*x(l1:l2)**1.5*max(0.01,sin(y(m)))**1.5
        z_over_h = x(l1:l2)*cos(y(m))/disk_height
      endif
!
!  compute alpha effect for EMF
!
!  mf_EMF for alpha effect dynamos
!
      if (lpencil(i_mf_EMF)) then
!
!  compute alpha profile (alpha_tmp)
!
        if (nxgrid/=1) kx=2*pi/Lx
        select case (alpha_profile)
        case ('const'); alpha_tmp=1.
        case ('box')
          if (y(m) >= y1_alp .and. y(m) <= y2_alp) then
            where(x(l1:l2)>=x1_alp .and. x(l1:l2)<=x2_alp)
              alpha_tmp=1.
            elsewhere
              alpha_tmp=0.
            endwhere
          else
            alpha_tmp=0.
          endif
        case ('coskx'); alpha_tmp=sqrt2*cos(kx_alpha*x(l1:l2))
        case ('siny'); alpha_tmp=sin(y(m))
        case ('sinz'); alpha_tmp=sin(z(n))
        case ('cos(z/2)'); alpha_tmp=cos(.5*z(n))
        case ('cos(z/2)_with_halo'); alpha_tmp=max(cos(.5*z(n)),0.)
        case ('patches')
          alpha_tmp=0.
          do ipatch=1, npatches_actual
            x1=2.*atan(tan(.5*(x(l1:l2)-xcenter(ipatch))))
            y1=2.*atan(tan(.5*(y(m    )-ycenter(ipatch))))
            z1=2.*atan(tan(.5*(z(n    )-zcenter(ipatch))))
            x11=+cy(ipatch)*x1+sy(ipatch)*z1
            z11=-sy(ipatch)*x1+cy(ipatch)*z1
            y11=            y1
            r2=.5*(x11/sigx)**2+.5*(y11/sigx)**2+.5*(z11/sigz)**2
            alpha_tmp=alpha_tmp+z11*exp(-r2)
            if (iux/=0) then
              if (ltest_patches) then
                p%uu(:,1)=p%uu(:,1)+ampluu_kinematic*    exp(-r2)
                p%uu(:,2)=p%uu(:,2)+ampluu_kinematic*x11*exp(-r2)
                p%uu(:,3)=p%uu(:,3)+ampluu_kinematic*z11*exp(-r2)
              else
                p%uu(:,1)=p%uu(:,1)-ampluu_kinematic*y11*exp(-r2)
                p%uu(:,2)=p%uu(:,2)+ampluu_kinematic*x11*exp(-r2)
!if (n>40 .and. n<128) print*,'AXEL1: z(n),p%uu(128,2)=',z(n),maxval(p%uu(:,2))
!print*,'AXEL1: z(n),p%uu(128,2)=',z(n),minval(p%uu(:,2)),maxval(p%uu(:,2))
!print*,'AXEL1: n,p%uu(128,2)=',n,maxval(p%uu(:,2))
              endif
            endif
          enddo
        case ('sphere')
          if (lspherical_coords) then
            rr=x(l1:l2)
          elseif (lcylindrical_coords) then
            rr=sqrt(x(l1:l2)**2+z(n)**2)
          else
            rr=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
          endif
          alpha_tmp=.5*(1.-erfunc((rr-alpha_rmax)/alpha_width))
        case ('z')
          if (alpha_rmax/=0.) then
            rr=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
            pp=.5*(1.-erfunc((rr-alpha_rmax)/alpha_width))
            alpha_tmp=z(n)*pp
          else
            alpha_tmp=z(n)
          endif
        case ('z/r')
          rr=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
          if (alpha_rmax/=0.) then
            pp=.5*(1.-erfunc((rr-alpha_rmax)/alpha_width))
            alpha_tmp=z(n)*pp/rr
          else
            alpha_tmp=z(n)/rr
          endif
        case ('1-tanhz'); alpha_tmp=.5*(1.-tanh(z(n)/alpha_width))
        case ('z/H'); alpha_tmp=z(n)/xyz1(3)
        case ('z/H_0'); alpha_tmp=z(n)/xyz1(3); if (z(n)==xyz1(3)) alpha_tmp=0.
        case ('y/H'); alpha_tmp=y(m)/xyz1(3)
        case ('J0x'); do l=l1,l2; alpha_tmp(l-l1+1)=bessj(0,k1bessel0*x(l)); enddo
        case ('J0x_2nd'); do l=l1,l2; alpha_tmp(l-l1+1)=bessj(0,k2bessel0*x(l)); enddo
        case ('cosy'); alpha_tmp=cos(y(m))
        case ('cosy*sin2y'); alpha_tmp=cos(y(m))*sin(y(m))**2
        case ('cosy*sin4y'); alpha_tmp=cos(y(m))*sin(y(m))**4
        case ('cosy*sin6y'); alpha_tmp=cos(y(m))*sin(y(m))**6
        case ('cosy*sin8y'); alpha_tmp=cos(y(m))*sin(y(m))**8
        case ('cosy*siny**n'); alpha_tmp=cos(y(m))*sin(y(m))**alpha_exp
        case ('surface_x*cosy'); alpha_tmp=0.5*(1.-erfunc((x(l1:l2)-x_surface)/alpha_width))*cos(y(m))
        case ('surface_x2*cosy'); alpha_tmp=0.25 &
          *(1.-erfunc((x(l1:l2)-x_surface)/alpha_width)) &
          *(1.+erfunc((x(l1:l2)-x_surface2)/alpha_width2))*cos(y(m))
        case ('sqrt_surface_x2*cosy')
          !alpha_tmp=sqrt(kf_x)*cos(y(m))
          alpha_tmp=kf_x1*cos(y(m))
          if (ip<=6.and.m==m1) print*,'alpha(sqrt_surface_x2)=',alpha_tmp
        case ('y*(1+eps*sinx)'); alpha_tmp=y(m)*(1.+alpha_eps*sin(kx*x(l1:l2)))
        case ('step-nhemi'); alpha_tmp=-tanh((y(m)-pi/2)/alpha_gap_step)
        case ('stepy'); alpha_tmp=-tanh((y(m)-yequator)/alpha_gap_step)
        case ('stepz'); alpha_tmp=-tanh((z(n)-zequator)/alpha_gap_step)
        case ('ystep-xcutoff')
           alpha_tmp=-tanh((y(m)-pi/2)/alpha_gap_step)&
             *(1+stepdown(x(l1:l2),alpha_rmax,alpha_width))
        case ('cosy*sin**n-xcutoff')
           alpha_tmp=(cos(y(m))*sin(y(m))**alpha_exp)&
             *(1+stepdown(x(l1:l2),alpha_rmax,alpha_width))
        case ('step-drop'); alpha_tmp=(1. &
                -step(y(m),pi/2.-alpha_equator_gap,alpha_gap_step) &
                -step(y(m),pi/2.+alpha_equator_gap,alpha_gap_step) &
                -step(alpha_cutoff_up,y(m),alpha_gap_step) &
                +step(y(m),alpha_cutoff_down,alpha_gap_step))
        case ('surface_z'); alpha_tmp=0.5*(1.-erfunc((z(n)-z_surface)/alpha_width))
        case ('above_x0'); alpha_tmp=0.5*(1.+erfunc((x(l1:l2)-alpha_rmin)/alpha_width))
        case ('z/H+surface_z'); alpha_tmp=(z(n)/z_surface)*0.5*(1.-erfunc((z(n)-z_surface)/alpha_width))
          if (headtt) print*,'alpha_profile=z/H+surface_z: z_surface,alpha_width=',z_surface,alpha_width
!
!  Galactic alpha effect profile
!
        case ('z/H*exp(-z2/2h2)')
          if (lcylindrical_coords) then
            if (headtt) print*,'z/H*exp(-z2/2h2) profile (cylindrical coords)',alpha_rmax
            if (alpha_rmax/=0.) then
              rr=x(l1:l2)
              pp=.5*(1.-erfunc((rr-alpha_rmax)/alpha_width2))
              alpha_tmp=z(n)/alpha_width*exp(-.5*(z(n)/alpha_width)**2)*pp
              !alpha_tmp=exp(-.5*(z(n)/alpha_width)**2)*pp
              if (alpha_pom0/=0.) alpha_tmp=alpha_tmp/(1.+(rr/alpha_pom0)**4)**.25
            else
              alpha_tmp=z(n)/alpha_width*exp(-.5*(z(n)/alpha_width)**2)
            endif
          else
            if (headtt) print*,'z/H*exp(-z2/2h2) profile (Cartesian)',alpha_rmax
            if (alpha_rmax/=0.) then
              rr=sqrt(x(l1:l2)**2+y(m)**2)
              pp=.5*(1.-erfunc((rr-alpha_rmax)/alpha_width2))
              alpha_tmp=z(n)/alpha_width*exp(-.5*(z(n)/alpha_width)**2)*pp
              if (alpha_pom0/=0.) alpha_tmp=alpha_tmp/(1.+(rr/alpha_pom0)**4)**.25
            else
              alpha_tmp=z(n)/alpha_width*exp(-.5*(z(n)/alpha_width)**2)
            endif
          endif
!
!  Galactic alpha effect profile with spiral
!
        case ('z/H*exp(-z2/2h2)*spiral')
          rr=sqrt(x(l1:l2)**2+y(m)**2)
          pp=atan2(y(m),x(l1:l2))
          pp0=2.*alog(rr)
          spiral=sin((pp-pp0))**2
          alpha_tmp=z(n)/alpha_width*exp(-.5*(z(n)/alpha_width)**2)*spiral
!
!  Shakura-Sunyaev disk, and alphaSS=alpha_effect
!
        case ('alphass')
          alpha_tmp = z_over_h*sqrt(p%cs2)*tanh(5.*(abs(z_over_h)-0.5))*0.5*(1.0+tanh(4.-2.*abs(z_over_h)))
!
!  alpha=alpha_effect^(1/2)*cs*rand.
!  Since we multiply alpha_effect later, here we divide by sqrt(alpha_effect) first.
!
        case ('fluc-alpha-disk')
          alpha_tmp=sqrt(p%cs2)/sqrt(alpha_effect)*rand_alpha(:,m)*0.5*(1.0+tanh(4.-2.*abs(z_over_h)))
        case ('z/H*erfunc(H-z)'); alpha_tmp=z(n)/xyz1(3)*erfunc((xyz1(3)-abs(z(n)))/alpha_width)
        case ('read'); alpha_tmp=alpha_input(l1:l2,m)
        case ('Jouve-2008-benchmark')
          alpha_tmp=(3.*sqrt(3.)/4.)*(sin(y(m)))**2*cos(y(m))*(1.0+erfunc((x(l1:l2)-0.7)/0.02))
        case ('nothing');
          call inevitably_fatal_error('calc_pencils_magnetic', &
            'alpha_profile="nothing" has been renamed to "const", please update your run.in')
        case default;
          call inevitably_fatal_error('calc_pencils_magnetic', &
            'alpha_profile no such alpha profile')
        endselect
!
!  compute delta effect profile (delta_tmp)
!
        select case (delta_profile)
        case ('const'); delta_tmp=1.
        case ('cos(z/2)_with_halo'); delta_tmp=max(cos(.5*z(n)),0.)
        case ('sincos(z/2)_with_halo'); delta_tmp=max(cos(.5*z(n)),0.)*sin(.5*z(n))
        case ('cos(theta)'); delta_tmp=cos(y(m))
        case ('sin(theta)'); delta_tmp=sin(y(m))
        case ('cosy*sin4y'); delta_tmp=cos(y(m))*sin(y(m))**4
        case ('tanhz'); delta_tmp=.5*(1.-tanh((abs(z(n))-meanfield_delta_height)/meanfield_delta_width))
        case default;
          call inevitably_fatal_error('calc_pencils_magnetic', &
            'delta_profile no such delta profile')
        endselect
!
!  Compute diffusion term.
!  This sets the meanfield_etat_tmp term at each time step.
!
        if (meanfield_etat/=0.0) then
          meanfield_etat_tmp=etat_x*etat_y(m)*etat_z(n)
          meanfield_getat_tmp(:,1)=detat_x*etat_y(m)*etat_z(n)
          meanfield_getat_tmp(:,2)=etat_x*detat_y(m)*etat_z(n)
          meanfield_getat_tmp(:,3)=etat_x*etat_y(m)*detat_z(n)
!
!  1/r correction for spherical symmetry (assuming axisymmetric profiles,
!  i.e.  meanfield_getat_tmp(:,3)=0.)
!
          if (lspherical_coords) then
            meanfield_getat_tmp(:,2)=r1_mn*meanfield_getat_tmp(:,2)
          endif
!
!  Magnetic etat quenching (contribution to pumping currently ignored)
!
          if (meanfield_etaB/=0.0) then
            meanfield_etaB2=meanfield_etaB**2
            meanfield_etat_tmp=meanfield_etat_tmp/sqrt(1.+p%b2/meanfield_etaB2)
!           call quench_simple(p%b2,meanfield_etaB2,meanfield_etat_tmp)
          endif
!
!  For accretion disks, etat~cs2/Omega
!
          if (meanfield_etat_profile=='alphass') then
            meanfield_etat_tmp=meanfield_etat_tmp*p%cs2  * (1.+z_over_h**2.)*3./8
            if (meanfield_etat_corona>0.0) then
              where (abs(z_over_h)>meanfield_etat_width) &
                meanfield_etat_tmp = meanfield_etat_corona*x(l1:l2)**meanfield_etat_exp
            endif
          endif
        endif
!
!  Here we initialize alpha_total.
!  Allow for the possibility of dynamical alpha.
!  By default, lmeanfield_noalpm=F, so normally treat
!  alpha_tmp as a general profile in front of the total alpha.
!  Invoke lmeanfield_noalpm to treat magnetic alpha separately
!  (which makes physically more sense!)
!
        if (lalpm .and. .not. lmeanfield_noalpm) then
          if (lalpha_profile_total) then
            alpha_total=(alpha_effect+f(l1:l2,m,n,ialpm))*alpha_tmp
            if (headtt) print*,'use alp=(alpK+alpM)*profile'
          else
            alpha_total=alpha_effect*alpha_tmp+f(l1:l2,m,n,ialpm)
            if (headtt) print*,'use alp=alpK*profile+alpM'
          endif
!
!  Possibility of *alternate* dynamical alpha.
!  Here we initialize alpha_total.
!
        elseif (lalpm_alternate .and. .not. lmeanfield_noalpm) then
          kf_tmp=kf_x*kf_y(m)*kf_z(n)
          prefact=meanfield_etat_tmp*(kf_tmp/meanfield_Beq)**2
          alpm=prefact*(f(l1:l2,m,n,ialpm)-p%ab)
          if (lalpha_profile_total) then
            alpha_total=(alpha_effect+alpm)*alpha_tmp
          else
            alpha_total=alpha_effect*alpha_tmp+alpm
          endif
        else
          alpha_total=alpha_effect*alpha_tmp
        endif
!
!  Possibility of conventional alpha quenching (rescales alpha_total).
!  Initialize EMF with alpha_total*bb.
!  Here we initialize p%mf_EMF.
!  NOTE: the following with alpha_quenching*sqrt(kf_x) was a hardwired fix
!  and and is now dealt with using meanfield_Beq_profile='sqrt(kf_x)'.
!  NOTE2: the calculation of Beq21 is here repeated, which should be avoided.
!
        if (alpha_quenching/=0.0) then
          select case (meanfield_Beq_profile)
          case ('uturbconst');
            Beq21=mu01*p%rho1/(uturb**2)
          case ('sqrt(kf_x)');
            Beq21=sqrt(kf_x)/meanfield_Beq**2
          case ('alphass');
            Beq21=1.0/max(1e-10, 3.*meanfield_Beq**2*mu0*p%pp-p%b2)
          case ('alphass-algebraic')
            Beq21=mu01/( 3*meanfield_Beq**2*(1.-meanfield_Beq)*p%pp )
          case ('fluc-alpha-disk')
            Beq21=(sqrt(fluc_alp_m)*sigma_alpha-1.)*4./meanfield_Beq/mu0/p%pp
          case default;
            Beq21=1./meanfield_Beq**2
          endselect
          alpha_quenching_tmp=alpha_quenching
          alpha_total=alpha_total/(1.+alpha_quenching_tmp*p%b2*Beq21)
        endif
!
!  Apply alpha effect; allow for anisotropy
!
        if (lalpha_aniso) then
          do j=1,3
            p%mf_EMF(:,j)=alpha_total*alpha_aniso(j)*p%bb(:,j)
          enddo
        elseif (lalpha_tensor) then
          do i=1,3
            p%mf_EMF(:,i)=0.
            do j=1,3
              if (alpha_zz /= 0.) then
                p%mf_EMF(:,i)=p%mf_EMF(:,i)+alpha_total &
                  *(alpha_tensor(i,j)+alpha_tensor_y(m-m1+1,i,j))*p%bb(:,j)
              elseif (lread_alpha_tensor_z_as_y) then
                p%mf_EMF(:,i)=p%mf_EMF(:,i)+ &
                  alpha_tensor_y(m-m1+1,i,j)*p%bb(:,j)
              elseif (lread_alpha_tensor_z) then
                if (i<=2.and.j<=2) then
                  !p%mf_EMF(:,i)=p%mf_EMF(:,i)+alpha_total &
!AB: working with alpha_total is questionable in this case
                  p%mf_EMF(:,i)=p%mf_EMF(:,i)+1. &
                      *alpha_tensor_z(n-n1+1,i,j)*p%bb(:,j)
                endif
              else
                p%mf_EMF(:,i)=p%mf_EMF(:,i)+alpha_total*alpha_tensor(i,j)*p%bb(:,j)
              endif
            enddo
          enddo
        else
          call multsv_mn(alpha_total,p%bb,p%mf_EMF)
        endif
!
!  Optionally, run the alpha-Omega approximation, i.e. apply
!  alpha effect only to the toroidal component. Since p%mf_EMF
!  was initialized only in the previous line, we can just set
!  the r and theta components to zero (in spherical coordinates).
!  In Cartesian coordinates, we want the x and z components to vanish,
!  so we keep only the y-components.
!
        if (lalpha_Omega_approx) then
          if (lspherical_coords) then
            p%mf_EMF(:,1:2)=0.
            call fatal_error("calc_pencils_magn_mf: ", &
                "lalpha_Omega_approx not implemented for this case")
          else
            p%mf_EMF(:,1)=0.
            p%mf_EMF(:,3)=0.
          endif
        endif
!
!  Apply eta tensor, but subtract part from etat for stability reasons.
!  In other words, any constant meanfield_etat should have formally no
!  effect, but is always needed because the pure Weyl gauge is unstable.
!
        if (leta_tensor) then
          if (lread_eta_tensor_z_as_y) then
            do i=1,3
              do j=1,3
                p%mf_EMF(:,i)=p%mf_EMF(:,i)-eta_tensor_y(m-m1+1,i,j)*p%jj(:,j)
              enddo
              p%mf_EMF(:,i)=p%mf_EMF(:,i)+meanfield_etat*p%jj(:,i)
            enddo
          else
            do i=1,2
              do j=1,2
                p%mf_EMF(:,i)=p%mf_EMF(:,i)-eta_tensor_z(n-n1+1,i,j)*p%jj(:,j)
              enddo
              p%mf_EMF(:,i)=p%mf_EMF(:,i)+meanfield_etat*p%jj(:,i)
            enddo
          endif
        endif
!
!  Add possible delta x J effect and turbulent diffusion to EMF.
!
        if (ldelta_profile) then
          p%mf_EMF=p%mf_EMF+delta_effect*p%oxJ
        else
          if (lspherical_coords) then
! assuming 1D in theta, delta is only delta_theta, J is only J_theta
            p%mf_EMF(:,1)=p%mf_EMF(:,1)+delta_effect*delta_tmp*p%jj(:,3)
            p%mf_EMF(:,3)=p%mf_EMF(:,3)-delta_effect*delta_tmp*p%jj(:,1)
          else
            p%mf_EMF(:,1)=p%mf_EMF(:,1)-delta_effect*delta_tmp*p%jj(:,2)
            p%mf_EMF(:,2)=p%mf_EMF(:,2)+delta_effect*delta_tmp*p%jj(:,1)
          endif
        endif
!
!  apply pumping effect: EMF=...-.5*grad(etat) x B
!
        if (lmeanfield_pumping .and. (meanfield_etat /= 0.0)) then
          fact=.5*meanfield_pumping
          call cross_mn(meanfield_getat_tmp, p%bb, getat_cross_B_tmp)
          p%mf_EMF=p%mf_EMF-fact*getat_cross_B_tmp
        endif
!
!  apply pumping effect in spherical coordinates
!
        if (gamma_effect/=0.) then
          fact=gamma_effect
          p%mf_EMF(:,2)=p%mf_EMF(:,2)-fact*p%bb(:,3)
          p%mf_EMF(:,3)=p%mf_EMF(:,3)+fact*p%bb(:,2)
        endif
!
!  Apply diffusion term: simple in Weyl gauge, which is not the default!
!  In diffusive gauge, add (divA) grad(etat) term.
!
        if (meanfield_etat /= 0.0) then
          if (lweyl_gauge) then
            call multsv_mn_add(-meanfield_etat_tmp,p%jj,p%mf_EMF)
          else
            call multsv_mn_add(meanfield_etat_tmp,p%del2a,p%mf_EMF)
            call multsv_mn_add(p%diva,meanfield_getat_tmp,p%mf_EMF)
          endif
        endif
!
!  Allow for possibility of variable etat from the f-array.
!
        if (ietat/=0) then
          call multsv_mn_add(-f(l1:l2,m,n,ietat),p%jj,p%mf_EMF)
        endif
!
!  Possibility of adding contribution from large-scale velocity.
!
        if (llarge_scale_velocity) p%mf_EMF=p%mf_EMF+p%uxb
!
!  Possibility of shear-current effect (ref: Zhou & Blackman2021).
!  Note that only the off-diagonal terms in eta_ij are included,
!  i.e., not including turbulent diffusion, but including delta
!  effect, and the shear-current effect itself.
!
        if (lshear_current_effect) then
          select case (shear_current_profile)
          case ('alphass')
!
!  thin-disk model. meanfield_Beq=alphaSS
!
            shear_current_sh_tmp = shear_current_sh*0.5*(1.0+tanh(4.-2.*abs(z_over_h))) &
                                     /( 1.0 + 3.0*p%b2/7.0/meanfield_Beq/mu0/p%pp )
            !
            p%mf_EMF(:,1) = p%mf_EMF(:,1) - meanfield_etat_tmp *          &
                     shear_current_sh_tmp * (0.63-0.26*b2_to_u2) * p%jj(:,3)
            p%mf_EMF(:,3) = p%mf_EMF(:,3) - meanfield_etat_tmp *          &
                   shear_current_sh_tmp * (-0.03-0.124*b2_to_u2) * p%jj(:,1)
          case default
            call fatal_error('calc_pencils_magn_mf','Invalid shear_current_profile!')
          endselect
        endif
!
!  Possibility of turning EMF to zero in a certain region.
!
        if (lEMF_profile) then
          select case (EMF_profile)
          case ('xcutoff');
            EMF_prof= 1+stepdown(x(l1:l2),alpha_rmax,alpha_width)
          case ('surface_z');
            EMF_prof=0.5*(1.-erfunc((z(n)-z_surface)/alpha_width))
          case ('nothing');
          call inevitably_fatal_error('calc_pencils_magnetic', &
            'lEMF_profile=T, but no profile selected !')
          endselect
            p%mf_EMF(:,1)=p%mf_EMF(:,1)*EMF_prof(:)
            p%mf_EMF(:,2)=p%mf_EMF(:,2)*EMF_prof(:)
            p%mf_EMF(:,3)=p%mf_EMF(:,3)*EMF_prof(:)
        endif
      endif
!
!  compute GW part,
!  Use also indices m and n, which are normally used to address positions in the f-array
!  when underneath the mn-loop, but this is not the case 
!
      if (lGW_tensor) then
        do j=1,3
        do k=1,3
        do i=1,3
        do nn=1,3
          p%mf_EMF(:,j)=p%mf_EMF(:,j)+GWfac1*levi_civita(k,nn,j)*hij(:,i,nn)*p%bij(:,k,i) &
                                     -GWfac2*levi_civita(k,nn,i)*hij(:,j,nn)*p%bij(:,k,i) &
                                     -GWfac3*levi_civita(k,nn,i)*hijk(:,j,nn,i)*p%bb(:,k)
        enddo
        enddo
        enddo
        enddo
      endif
!
!  Evaluate exa, but do this at the end after mf_EMF has been
!  fully assembled.
!
      if (lpencil(i_exa)) then
        if (lmagn_mf) then
          call cross_mn(-p%mf_EMF,p%aa,exa_meanfield)
          p%exa=p%exa+exa_meanfield
        endif
      endif
!
!  EMFdotB
!
      if (lpencil(i_mf_EMFdotB)) call dot_mn(p%mf_EMF,p%bb,p%mf_EMFdotB)
!
!  Compute chiB_term in specific entropy equation.
!  Ds/Dt = ... (1/rho*T)*div[rho*T*chit(B)*grads]
!        = ... chit(B)*{[gradln(rho*T)+gradln(chit(B))].gs+del2s}
!        = ... chit(B)*[gradln(rho*T).grads+del2s] + chit(B)*gradln(chit(B)).gs
!
      if (lpencil(i_chiB_mf)) then
        select case (meanfield_Beq_profile)
        case ('uturbconst');
          Beq21=mu01*p%rho1/(uturb**2)
        case default;
          Beq21=1./meanfield_Beq**2
        endselect
!
!  Choice between 2 versions:
!
        oneQbeta2=1.+chit_quenching*p%b2*Beq21
!
!  Currently no additional profile, so zero gradient.
!
        chit_prof=1.
        glnchit_prof=0.
!
!  Add contribution from gradient of chiB*B^2/Beq^2 term.
!  The B2glnrho term is critical for the magneto-thermodiffusive instability.
!
        call multsv_mn(p%b2,p%glnrho,B2glnrho)
        if (lignore_gradB2_inchiB) then
          call multsv_mn(-chit_quenching*Beq21/oneQbeta2,-B2glnrho,glnchit)
        else
          call multmv_transp(p%bij,p%bb,Bk_Bki) !=1/2 grad B^2
          call multsv_mn(-chit_quenching*Beq21/oneQbeta2,2.*Bk_Bki-B2glnrho,glnchit)
        endif
!
!  In the presence of a uniform imposed field, chit or calKt are no
!  longer constant. To prevent this, we compensate by a 1+Q*B0^2 factor.
!  This leads to an additional contribution in the glnchit expression.
!
        if (lchit_Bext2_equil) then
          oneQbeta02=1.+chit_quenching*B_ext2*Beq21
          call multsv_mn(-chit_quenching*B_ext2*Beq21/oneQbeta02,p%glnrho,glnchit2)
          glnchit=glnchit+glnchit2
        endif
!
!       if (lchit_with_glnTT) then
!         call dot(p%glnrho+p%glnTT+glnchit_prof+glnchit,p%gss,g2)
!         p%chiB_mf=chi_t0*quench_chiB*(g2+p%del2ss)
!       else
!
!  The lrho_chit option corresponds to solving the equation
!  Ds/Dt = (calKt/rho)*(glncalKt.grads + del2s). Otherwise we have
!  Ds/Dt = chit*[(glnrho+glncalKt).grads + del2s].
!  This corresponds to turbulent diffusion with entropy gradient.
!
          if (lrho_chit) then
            call dot(glnchit_prof+glnchit,p%gss,g2)
            if (lchit_Bext2_equil) then
              p%chiB_mf=p%rho1*chi_t0*oneQbeta02/oneQbeta2*(g2+p%del2ss)
            else
              p%chiB_mf=p%rho1*chi_t0/oneQbeta2*(g2+p%del2ss)
            endif
!
!  turbulent diffusion with temperature gradient
!
          elseif (lturb_temp_diff) then
            call dot(p%glnrho+p%glnTT+glnchit_prof+glnchit,p%glnTT,g2)
            if (lchit_Bext2_equil) then
              p%chiB_mf=p%cp*chi_t0*oneQbeta02/oneQbeta2*(g2+p%del2lnTT)
            else
              p%chiB_mf=p%cp*chi_t0/oneQbeta2*(g2+p%del2lnTT)
            endif
!
!  Turbulent diffusion with entropy gradient and coefficient
!  proportional to chi_t*rho, without temperature factor.
!
          else
            call dot(p%glnrho+glnchit_prof+glnchit,p%gss,g2)
            if (lchit_Bext2_equil) then
              p%chiB_mf=chi_t0*oneQbeta02/oneQbeta2*(g2+p%del2ss)
            else
              p%chiB_mf=chi_t0/oneQbeta2*(g2+p%del2ss)
            endif
          endif
!       endif
!?? end of wrong indentation??
!?? end of chit and chiB apparently (!AB)
      endif
!
!  Calculate diagnostics.
!
      if (ldiagnos) then
        if (idiag_qsm/=0) call sum_mn_name(meanfield_qs_func,idiag_qsm)
        if (idiag_qpm/=0) call sum_mn_name(meanfield_qp_func,idiag_qpm)
        if (idiag_qem/=0) call sum_mn_name(meanfield_qe_func,idiag_qem)
        if (idiag_alpKm/=0) call sum_mn_name(alpha_total,idiag_alpKm)
        if (idiag_alpKjbm/=0) call sum_mn_name(alpha_total*p%jb,idiag_alpKjbm)
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1.
!
      if (l1davgfirst .or. (ldiagnos .and. ldiagnos_need_zaverages)) then
        call xysum_mn_name_z(meanfield_qp_func,idiag_qpmz)
      endif
!
!  2-D averages.
!  y-averaged alpha effect
!
      if (l2davgfirst) then
        if (idiag_alpmxz/=0)  then
          call ysum_mn_name_xz(alpha_total,idiag_alpmxz)
        endif
      endif
!
    endsubroutine calc_pencils_magn_mf
!***********************************************************************
    subroutine daa_dt_meanfield(f,df,p)
!
!  add mean-field evolution to magnetic field.
!
!  27-jul-10/axel: coded
!
      use Sub, only: cross_mn
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(inout) :: df,f
!
      real, dimension(nx) :: diffus_eta
      real, dimension(nx,3) :: tmpv
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'daa_dt_meanfield: SOLVE'
!
!  Add jxb/rho to momentum equation.
!
      if (lhydro) then
        if (lmeanfield_jxb) df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+p%jxbr_mf
      endif
!
!  Add div(chit(B)*rho*T*grads) to entropy equation
!AB: the following is not correct, so disable it for now!
!
!      if (lentropy) then
!        if (lmeanfield_chitB) df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+p%chiB_mf
!      endif
!
!  Multiply resistivity by Nyquist scale, for resistive time-step.
!  We include possible contribution from meanfield_etat, which is however
!  only invoked in mean field models.
!  Allow for variable etat (mean field theory).
!
      if (lfirst.and.ldt) then
        diffus_eta=meanfield_etat*dxyz_2
        if (headtt.or.ldebug) &
          print*, 'daa_dt_meanfield: max(diffus_eta)  =', maxval(diffus_eta)
        maxdiffus=max(maxdiffus,diffus_eta)
      endif
!
!  Alpha effect.
!  Additional terms if Mean Field Theory is included.
!
      if (lmagn_mf .and. &
        (meanfield_etat/=0.0 .or. ietat/=0 .or. &
        alpha_effect/=0.0 .or. delta_effect/=0.0) .or. &
        gamma_effect/=0.0 .or. lread_alpha_tensor_z .or. &
        lread_eta_tensor_z .or. lread_alpha_tensor_z_as_y .or. &
        lread_eta_tensor_z_as_y .or. lGW_tensor) then
!
!  Apply p%mf_EMF only if .not.lmagn_mf_demfdt; otherwise postpone.
!
        if (.not.lmagn_mf_demfdt) then
          df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+p%mf_EMF
        endif
!
!  Apply Omega effect. This is currently applied directly to A
!  and is therefore not part of the pencil p%mf_EMF.
!
        if (lOmega_effect) call Omega_effect(f,df,p)
      endif
!
!  Impose a By(x)=B0*sinkx field by adding a dAz/dt= ... + (B/tau*k)*coskx term.
!
      if (lrhs_term) df(l1:l2,m,n,iaz)=df(l1:l2,m,n,iaz)+rhs_term
!
!  2-component external emf
!
      if (lrhs_term2) then
        df(l1:l2,m,n,iay)=df(l1:l2,m,n,iay)+rhs_termy
        df(l1:l2,m,n,iaz)=df(l1:l2,m,n,iaz)+rhs_termz
      endif
!
!  Apply here the Omega effect when the mean flow is only calculated in
!  this mean-field module. In that case, it was not previously applied.
!
      if (lOmega_effect_meanfield) then
        call cross_mn(p%uu,p%bb,tmpv)
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+tmpv
!if (n>90 .and. n<128) print*,'AXEL2: n,p%uu(128,2)=',n,p%uu(128,2)
!print*,'AXEL2: z(n),p%uu(128,2)=',z(n),minval(p%uu(:,2)),maxval(p%uu(:,2))
      endif
!
!  Time-advance of secondary mean-field modules.
!
      if (lmagn_mf_demfdt) call demf_dt_meanfield(f,df,p)
!
      call calc_diagnostics_meanfield(f,p)
!
    endsubroutine daa_dt_meanfield
!***********************************************************************
    subroutine calc_diagnostics_meanfield(f,p)
!
!  Calculate diagnostic quantities.
!  Diagnostic output for mean field dynamos.
!
      use Diagnostics

      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (ldiagnos) then
        if (idiag_EMFdotBm/=0) call sum_mn_name(p%mf_EMFdotB,idiag_EMFdotBm)
        if (idiag_EMFdotB_int/=0) call integrate_mn_name(p%mf_EMFdotB,idiag_EMFdotB_int)
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1.
!
      if (l1davgfirst .or. (ldiagnos .and. ldiagnos_need_zaverages)) then
!
!  Calculate magnetic helicity flux (ExA contribution).
!
        if (idiag_EMFmz1/=0) call xysum_mn_name_z(p%mf_EMF(:,1),idiag_EMFmz1)
        if (idiag_EMFmz2/=0) call xysum_mn_name_z(p%mf_EMF(:,2),idiag_EMFmz2)
        if (idiag_EMFmz3/=0) call xysum_mn_name_z(p%mf_EMF(:,3),idiag_EMFmz3)
      endif
!
!  2-D averages.
!  y-averaged effective magnetic pressure,
!  only makes sense for the 'uturbconst' profile
!
      if (l2davgfirst) then
        if (idiag_peffmxz/=0) call ysum_mn_name_xz(.5*(1.-p%mf_qp)*p%b2*p%mf_Beq21,idiag_peffmxz)
      endif
!  Note that this does not necessarily happen with ldiagnos=.true.
!
!     if (l2davgfirst) then
!       if (idiag_Ezmxz/=0) call ysum_mn_name_xz(p%uxb(:,3),idiag_Ezmxz)
!     else
!
!  We may need to calculate bxmxy without calculating bmx. The following
!  if condition was messing up calculation of bmxy_rms
!
!       if (ldiagnos) then
!         if (idiag_bxmxy/=0) call zsum_mn_name_xy(p%bb(:,1),idiag_bxmxy)
!       endif
!     endif
!
!  Diagnostics of secondary mean-field modules.
!
      if (lmagn_mf_demfdt) call calc_diagnostics_dt_meanfield(f)
!
    endsubroutine  calc_diagnostics_meanfield
!***********************************************************************
    subroutine meanfield_chitB(rho,b2,quench)
!
!  Calculate magnetic properties needed for z boundary conditions.
!  This routine contails calls to more specialized routines.
!
!   8-jun-13/axel: coded
!
      real, dimension (nx) :: rho,b2,Beq21,quench
!
!  compute Beq21 = 1/Beq^2
!
      select case (meanfield_Beq_profile)
      case ('uturbconst');
        Beq21=mu01/(rho*uturb**2)
      case default;
        Beq21=1./meanfield_Beq**2
      endselect
!
!  compute chit_quenching
!
      quench=1./(1.+chit_quenching*b2*Beq21)
!
    endsubroutine meanfield_chitB
!***********************************************************************
    subroutine Omega_effect(f,df,p)
!
!  Omega effect coded (normally used in context of mean-field theory)
!  Can do uniform shear (0,Sx,0), and the cosx*cosz profile (solar CZ).
!  In most cases the Omega effect can be modeled using hydro_kinematic,
!  but this is not possible when the flow varies in a direction that
!  is not a coordinate direction, e.g. when U=(0,Sx,0) and A=A(z,t).
!  In such cases the Omega effect must be rewritten in terms of
!  velocity gradients operating on A, i.e. (gradU)^T.A, instead of UxB.
!
!  30-apr-05/axel: coded
!
      use Sub, only: stepdown
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real kx
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  use gauge transformation, uxB = -Ay*grad(Uy) + gradient-term
!
      select case (Omega_profile)
      case ('nothing'); print*,'Omega_profile=nothing'
      case ('(0,Sx,0)')
        if (headtt) print*,'Omega_effect: uniform shear in x, S=',Omega_ampl
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)-Omega_ampl*f(l1:l2,m,n,iay)
      case ('(0,0,Sx)')
        if (headtt) print*,'Omega_effect: uniform shear in x, S=',Omega_ampl
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)-Omega_ampl*f(l1:l2,m,n,iaz)
      case ('(0,0,pomSx)')
        if (headtt) print*,'Omega_effect: uniform shear in radius, S=',Omega_ampl
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)-Omega_ampl*f(l1:l2,m,n,iaz)*p%x_mn
      case ('(Sz,0,0)')
        if (headtt) print*,'Omega_effect: uniform shear in z, S=',Omega_ampl
        df(l1:l2,m,n,iaz)=df(l1:l2,m,n,iaz)-Omega_ampl*f(l1:l2,m,n,iax)
        if (lhydro) df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-Omega_ampl*f(l1:l2,m,n,iuz)
      case ('(0,cosx*cosz,0)')
        if (headtt) print*,'Omega_effect: solar shear, S=',Omega_ampl
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)+Omega_ampl*f(l1:l2,m,n,iay) &
            *sin(x(l1:l2))*cos(z(n))
        df(l1:l2,m,n,iaz)=df(l1:l2,m,n,iaz)+Omega_ampl*f(l1:l2,m,n,iay) &
            *cos(x(l1:l2))*sin(z(n))
      case ('(0,0,cosx)')
        kx=2*pi/Lx
        if (headtt) print*,'Omega_effect: (0,0,cosx), S,kx=',Omega_ampl,kx
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)+Omega_ampl*f(l1:l2,m,n,iaz) &
            *kx*sin(kx*x(l1:l2))
      case ('(0,0,siny)')
        if (headtt) print*,'Omega_effect: (0,0,siny), Omega_ampl=',Omega_ampl
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)+Omega_ampl*f(l1:l2,m,n,iaz) &
            *sin(y(m))
      case('rcutoff_sin_theta')
        if (headtt) print*,'Omega_effect: r cutoff sin(theta), Omega_ampl=',Omega_ampl
!        df(l1:l2,m,n,iax)=Omega_ampl*df(l1:l2,m,n,iax)
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)-Omega_ampl*p%bb(:,2) &
            *sin(y(m))*x(l1:l2)*(1+stepdown(x(l1:l2),Omega_rmax,Omega_rwidth))
        df(l1:l2,m,n,iay)=df(l1:l2,m,n,iay)+Omega_ampl*p%bb(:,1) &
            *sin(y(m))*x(l1:l2)*(1+stepdown(x(l1:l2),Omega_rmax,Omega_rwidth))
      case default; print*,'Omega_profile=unknown'
      endselect
!
    endsubroutine Omega_effect
!***********************************************************************
    subroutine read_magn_mf_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=magn_mf_init_pars, IOSTAT=iostat)
!
!  read namelist for secondary modules in mean-field theory (if invoked)
!
      if (lmagn_mf_demfdt) call read_magn_mf_demfdt_init_pars(iostat)
!
    endsubroutine read_magn_mf_init_pars
!***********************************************************************
    subroutine write_magn_mf_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=magn_mf_init_pars)
!
!  write namelist for secondary modules in mean-field theory (if invoked)
!
      if (lmagn_mf_demfdt) call write_magn_mf_demfdt_init_pars(unit)
!
    endsubroutine write_magn_mf_init_pars
!***********************************************************************
    subroutine read_magn_mf_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=magn_mf_run_pars, IOSTAT=iostat)
!
!  read namelist for secondary modules in mean-field theory (if invoked)
!
      if (lmagn_mf_demfdt) call read_magn_mf_demfdt_run_pars(iostat)
!
    endsubroutine read_magn_mf_run_pars
!***********************************************************************
    subroutine write_magn_mf_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=magn_mf_run_pars)
!
!  write namelist for secondary modules in mean-field theory (if invoked)
!
      if (lmagn_mf_demfdt) call write_magn_mf_demfdt_run_pars(unit)
!
    endsubroutine write_magn_mf_run_pars
!***********************************************************************
    subroutine rprint_magn_mf(lreset,lwrite)
!
!  Reads and registers print parameters relevant for magnetic fields.
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics, only: parse_name
!
      integer :: iname,inamey,inamez,ixy,ixz,irz,inamer
!      integer :: inamex
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of RELOAD.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_qsm=0; idiag_qpm=0; idiag_qem=0;
        idiag_EMFmz1=0; idiag_EMFmz2=0; idiag_EMFmz3=0; idiag_peffmxz=0
        idiag_qpmz=0; idiag_alpmxz=0
        idiag_alpKm=0; idiag_alpKjbm=0
      endif
!
!  Check for those quantities that we want to evaluate online.
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'qsm',idiag_qsm)
        call parse_name(iname,cname(iname),cform(iname),'qpm',idiag_qpm)
        call parse_name(iname,cname(iname),cform(iname),'qem',idiag_qem)
        call parse_name(iname,cname(iname),cform(iname),'alpKm',idiag_alpKm)
        call parse_name(iname,cname(iname),cform(iname),'alpKjbm',idiag_alpKjbm)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
!     do inamex=1,nnamex
!    enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'EMFmz1',idiag_EMFmz1)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'EMFmz2',idiag_EMFmz2)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'EMFmz3',idiag_EMFmz3)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'qpmz',idiag_qpmz)
      enddo
!
!  Check for those quantities for which we want y-averages.
!
      do ixz=1,nnamexz
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'peffmxz',idiag_peffmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'alpmxz',idiag_alpmxz)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do ixy=1,nnamexy
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do irz=1,nnamerz
      enddo
!
!  Check for those quantities for which we want phiz-averages.
!
      do inamer=1,nnamer
      enddo
!
!  Write column, idiag_XYZ, where our variable XYZ is stored.
!
      if (lwr) then
      endif
!
    endsubroutine rprint_magn_mf
!***********************************************************************
    subroutine pc_aasb_const_alpha(f,topbot,j)
!
!  Perfect-conductor BC for mean field in 1D (only z dependent) with constant alpha and etaT,
!  IF A IS CONSIDERED AS B. 
!  Only implemented at lower z boundary.
!
!  2-jan-17/MR: coded
!
    real, dimension(:,:,:,:) :: f
    integer :: j
    integer :: topbot

    real :: g0, D
    integer :: k

      if (topbot==BOT) then

        if (j/=iax) return

        g0=-147./(60.*dz); D=(alpha_effect**2+g0**2)*60.*dz
        k=n1
        f(:,:,n1,iax) = ( 360.*(-g0*f(:,:,k+1,iax)-alpha_effect*f(:,:,k+1,iay)) &
                         -450.*(-g0*f(:,:,k+2,iax)-alpha_effect*f(:,:,k+2,iay)) &
                         +400.*(-g0*f(:,:,k+3,iax)-alpha_effect*f(:,:,k+3,iay)) &
                         -225.*(-g0*f(:,:,k+4,iax)-alpha_effect*f(:,:,k+4,iay)) &
                         + 72.*(-g0*f(:,:,k+5,iax)-alpha_effect*f(:,:,k+5,iay)) &
                         - 10.*(-g0*f(:,:,k+6,iax)-alpha_effect*f(:,:,k+6,iay))  )/D

        f(:,:,n1,iay) = ( 360.*(-g0*f(:,:,k+1,iay)+alpha_effect*f(:,:,k+1,iax)) &
                         -450.*(-g0*f(:,:,k+2,iay)+alpha_effect*f(:,:,k+2,iax)) &
                         +400.*(-g0*f(:,:,k+3,iay)+alpha_effect*f(:,:,k+3,iax)) &
                         -225.*(-g0*f(:,:,k+4,iay)+alpha_effect*f(:,:,k+4,iax)) &
                         + 72.*(-g0*f(:,:,k+5,iay)+alpha_effect*f(:,:,k+5,iax)) &
                         - 10.*(-g0*f(:,:,k+6,iay)+alpha_effect*f(:,:,k+6,iax))  )/D
      else 
        call fatal_error('pc_aasb_const_alpha','topbot should be BOT or TOP')
      endif
      
    endsubroutine pc_aasb_const_alpha
!***********************************************************************
    subroutine meanfield_after_boundary(f)
!
!  precalculate parameters that are new at each timestep,
!  but the same for all pencils
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      if (alpha_profile=='fluc-alpha-disk') call refresh_fluc_alpha_phase
!
      call keep_compiler_quiet(f)
!
    endsubroutine meanfield_after_boundary
!***********************************************************************
    subroutine calc_fluc_alp_cells
!
!  calculate the location-dependent length and time scales for alpha
!  fluctuations in the Shakura-Sunyaev accretion disk model.
!
      use Gravity, only: g0
!
      integer :: i,j
      real :: r0,s0,pos_tmp,s_tmp
      real :: epsi
!
!  determine alpha cells iteratively:
!  pos(1)=xyz0(1), pos(i+1)=pos(i)+s(i),
!  s(i)= m * sqrt(alphass) * h0 * (r/r0)**exp
!  m=fluc_alp_m>1.0 determines how much larger is the coherent structure than the turbulence scale
!  alphass=alpha_effect is the Shakura-Sunyaev alpha parameter
!  h0=cs0/Omega0 is the disk scale height at r=r0
!  exp=alpha_exp determines how turbulence length scale with R
!
      r0=1.0
      s0 = fluc_alp_m * sqrt(alpha_effect) * meanfield_cs0/sqrt(g0/r0**3)
      pos_tmp = r0
      nr_cell=0   !  total number of cells in the r direction
      nth_cell=0  !  maximal number of cells in the theta direction
      i=1
      do while (pos_tmp<=xyz1(1))
        s_tmp = s0 * (pos_tmp/r0)**alpha_exp
        pos_tmp = pos_tmp + s_tmp
        epsi = s_tmp/(pos_tmp-0.5*s_tmp)
        nr_cell=nr_cell+1
        nth_cell=max(nth_cell,ceiling(xyz1(2)/epsi))
!
!  the grid at (l,m)=(i,j) belongs to the ir_cell(i)-th alpha cell in the radial direction,
!  and the itheta_cell(j)-th cell in the polar direction.
!
        do while (i<=nx .and. x(l1+i-1)<pos_tmp)
          ir_cell(i) = nr_cell
          tau_cell(i) = fluc_alp_m/sqrt(g0/(pos_tmp - 0.5*s_tmp)**3)
          itheta_cell(i,:) = ceiling((0.001+y(:))/epsi)
          !
          i=i+1
        enddo
      enddo
!
    endsubroutine calc_fluc_alp_cells
!***********************************************************************
    subroutine refresh_fluc_alpha_phase
!
!  For each alpha cell, update its value every tau_cell time
!
      use General, only: random_number_wrapper
!
      integer :: ix,iy,i,j
      real, dimension(:,:,:), allocatable :: rand1,rand2,update_chance
!
!  Skipping it=1 to avoid failing the pencil check
!
      if (it>1 .and. lfirst) then
        allocate( update_chance (nr_cell,nth_cell,1) )
        call random_number_wrapper(update_chance,CHANNEL=channel_magn_mf)
!
        select case (fluc_alp_profile)
        case ('gaussian')
          allocate( rand1 (nr_cell,nth_cell,1) )
          allocate( rand2 (nr_cell,nth_cell,1) )
!
          call random_number_wrapper(rand1,CHANNEL=channel_magn_mf)
          call random_number_wrapper(rand2,CHANNEL=channel_magn_mf)
!
          do ix=1,nx
            i = ir_cell(ix)
            do iy=1,my
              j = itheta_cell(ix,iy)
!
!  Use Box-Muller transform for a Gaussian ditribution N(0,1) for rand_alpha
!
              if (update_chance(i,j,1)<=dt/tau_cell(ix)) &
                  rand_alpha(ix,iy)=sqrt(-2.0*log(rand1(i,j,1)))*cos(2.0*pi*rand2(i,j,1))
            enddo
          enddo
          rand_alpha=sigma_alpha*rand_alpha
          where (rand_alpha>1.) rand_alpha=1.
          where (rand_alpha<-1.) rand_alpha=-1.
          deallocate(rand1,rand2)
!
        case ('uniform')
          allocate( rand1 (nr_cell,nth_cell,1) )
          call random_number_wrapper(rand1,CHANNEL=channel_magn_mf)
!
          do ix=1,nx
            i = ir_cell(ix)
            do iy=1,my
              j = itheta_cell(ix,iy)
              if (update_chance(i,j,1)<=dt/tau_cell(ix)) rand_alpha(ix,iy)=2.*rand1(i,j,1)-1.
            enddo
          enddo
          deallocate(rand1)
        endselect
!
        deallocate(update_chance)
      endif
!
    endsubroutine refresh_fluc_alpha_phase
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General, only: string_to_enum

    integer, parameter :: n_pars=1100
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(etat_x,p_par(1)) ! (nx)
    call copy_addr(etat_y,p_par(2)) ! (my)
    call copy_addr(etat_z,p_par(3)) ! (mz)
    call copy_addr(detat_x,p_par(4)) ! (nx)
    call copy_addr(detat_y,p_par(5)) ! (my)
    call copy_addr(detat_z,p_par(6)) ! (mz)

    endsubroutine
!***********************************************************************
endmodule Magnetic_meanfield
