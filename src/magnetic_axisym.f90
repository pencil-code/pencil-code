! $Id$
!
!  This modules deals with all aspects of magnetic fields; if no
!  magnetic fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  magnetically relevant subroutines listed in here.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagnetic = .true.
!
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED aa(3); a2; aij(3,3); bb(3); bbb(3); ab; uxb(3)
! PENCILS PROVIDED b2; bij(3,3); del2a(3); graddiva(3); jj(3)
! PENCILS PROVIDED j2; jb; va2; jxb(3); jxbr(3); jxbr2; ub; uxb(3); uxb2
! PENCILS PROVIDED uxj(3); beta1; uga(3); djuidjbi; jo
! PENCILS PROVIDED ujxb; oxu(3); oxuxb(3); jxbxb(3); jxbrxb(3)
! PENCILS PROVIDED glnrhoxb(3); del4a(3); del6a(3); oxj(3); diva
! PENCILS PROVIDED jij(3,3); sj; ss12; mf_EMF(3); mf_EMFdotB
! PENCILS PROVIDED cosjb,jparallel;jperp
! PENCILS PROVIDED cosub
!
!***************************************************************
module Magnetic
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'record_types.h'
  include 'magnetic.h'
!
! Slice precalculation buffers
!
  real, target, dimension (nx,ny,3) :: bb_xy,jj_xy
  real, target, dimension (nx,ny,3) :: bb_xy2,jj_xy2
  real, target, dimension (nx,ny,3) :: bb_xy3,jj_xy3
  real, target, dimension (nx,ny,3) :: bb_xy4,jj_xy4
  real, target, dimension (nx,nz,3) :: bb_xz,jj_xz
  real, target, dimension (ny,nz,3) :: bb_yz,jj_yz
!
  real, target, dimension (nx,ny) :: b2_xy,jb_xy
  real, target, dimension (nx,ny) :: b2_xy2,jb_xy2
  real, target, dimension (nx,ny) :: b2_xy3,jb_xy3
  real, target, dimension (nx,ny) :: b2_xy4,jb_xy4
  real, target, dimension (ny,nz) :: b2_yz,jb_yz
  real, target, dimension (nx,nz) :: b2_xz,jb_xz
!
  real, target, dimension (nx,ny) :: beta1_xy
  real, target, dimension (nx,ny) :: beta1_xy2
  real, target, dimension (nx,ny) :: beta1_xy3
  real, target, dimension (nx,ny) :: beta1_xy4
  real, target, dimension (ny,nz) :: beta1_yz
  real, target, dimension (nx,nz) :: beta1_xz
!
  real, dimension (mx,my) :: alpha_input
!
! Parameters
!
  integer, parameter :: nresi_max=4
!
  real, dimension (ninit) :: amplaa=0.0,kx_aa=1.,ky_aa=1.,kz_aa=1.
  real, dimension (ninit) :: phasex_aa=0., phasey_aa=0., phasez_aa=0.
  character (len=labellen), dimension(ninit) :: initaa='nothing'
  character (len=labellen) :: borderaa='nothing'
  character (len=labellen), dimension(nresi_max) :: iresistivity=''
  character (len=labellen) :: alpha_xprofile='const',alpha_yprofile='const'
  character (len=labellen) :: Omega_xprofile='const',Omega_yprofile='const'
  character (len=labellen) :: fring_profile='tanh'
  ! input parameters
  complex, dimension(3) :: coefaa=(/0.,0.,0./), coefbb=(/0.,0.,0./)
  real, dimension(3) :: B_ext=(/0.,0.,0./),B1_ext,B_ext_tmp,eta_aniso_hyper3
  real, dimension(3) :: axisr1=(/0,0,1/),dispr1=(/0.,0.5,0./)
  real, dimension(3) :: axisr2=(/1,0,0/),dispr2=(/0.,-0.5,0./)
  real, dimension(3) :: axisr3=(/1,0,0/),dispr3=(/0.,-0.5,0./)
!
!  profile functions
!
  real, dimension(nx) :: alpha_x,dalpha_x,Omega_x,dOmega_x
  real, dimension(my) :: alpha_y,dalpha_y,Omega_y,dOmega_y,Omega_tmp0
!
  real, target :: zmode=1. !(temporary)
!
!  profile parameters
!
  real :: alpha_x1=0.,alpha_dx1=0.
  real :: Omega_x1=0.,Omega_dx1=0.
  real :: Omega_yc1=0.,Omega_yc2=0.
!
  real :: fring1=0.,Iring1=0.,Rring1=1.,wr1=0.3
  real :: fring2=0.,Iring2=0.,Rring2=1.,wr2=0.3
  real :: fring3=0.,Iring3=0.,Rring3=1.,wr3=0.3
  real :: radius=.1,epsilonaa=1e-2,widthaa=.5,x0aa=0.,z0aa=0.
  real :: by_left=0.,by_right=0.,bz_left=0.,bz_right=0.
  real :: ABC_A=1.,ABC_B=1.,ABC_C=1.
  real :: bthresh=0.,bthresh_per_brms=0.,brms=0.,bthresh_scl=1.
  real :: eta_shock=0.
  real :: rhomin_jxb=0.,va2max_jxb=0.
  real :: omega_Bz_ext=0.
  real :: mu_r=-0.5 !(still needed for backwards compatibility)
  real :: mu_ext_pot=-0.5,inclaa=0.
  real :: mu012=.5 !(=1/2mu0)
  real :: rescale_aa=0.
  real :: ampl_B0=0.,D_smag=0.17,B_ext21,B_ext11
  real :: Omega_ampl=0.
  real :: rmode=1.,rm_int=0.,rm_ext=0.
  real :: alpha_effect=0.,alpha_quenching=0.,delta_effect=0.,meanfield_etat=0.
  real :: alpha_eps=0.
  real :: alpha_equator=impossible,alpha_equator_gap=0.,alpha_gap_step=0.
  real :: alpha_cutoff_up=0.,alpha_cutoff_down=0.
  real :: meanfield_Qs=1.,meanfield_Qp=1.
  real :: meanfield_Bs=1.,meanfield_Bp=1.
  real :: meanfield_kf=1.,meanfield_etaB=0.
  real :: displacement_gun=0.
  real :: pertamplaa=0.
  real :: initpower_aa=0.,cutoff_aa=0.,brms_target=1.,rescaling_fraction=1.
  real :: phase_beltrami=0., ampl_beltrami=0.
  real :: bmz=0, bmz_beltrami_phase=0.
  real :: taareset=0.,daareset=0.
  real :: center1_x=0., center1_y=0., center1_z=0.
  real :: fluxtube_border_width=impossible
  integer :: nbvec,nbvecmax=nx*ny*nz/4,va2power_jxb=5
  integer :: N_modes_aa=1, naareset
  integer :: nrings=2
  logical :: lpress_equil=.false., lpress_equil_via_ss=.false.
  logical :: llorentzforce=.true.,linduction=.true.
  logical :: lalpha_phi_equation=.true.
  logical :: lresi_eta_const=.false.
  logical :: lresi_etaSS=.false.
  logical :: lresi_hyper2=.false.
  logical :: lresi_hyper3=.false.
  logical :: lresi_hyper3_polar=.false.
  logical :: lresi_hyper3_strict=.false.
  logical :: lresi_zdep=.false., lresi_dust=.false., lresi_rdep=.false.
  logical :: lresi_hyper3_aniso=.false.
  logical :: lresi_eta_shock=.false.
  logical :: lresi_eta_shock_perp=.false.
  logical :: lresi_shell=.false.
  logical :: lresi_smagorinsky=.false.
  logical :: lresi_smagorinsky_cross=.false.
  logical, target, dimension (3) :: lfrozen_bb_bot=(/.false.,.false.,.false./)
  logical, target, dimension (3) :: lfrozen_bb_top=(/.false.,.false.,.false./)
  logical :: reinitialize_aa=.false., lohmic_heat=.true., lneutralion_heat=.true.
  logical :: lB_ext_pot=.false.
  logical :: lforce_free_test=.false.
  logical :: lmeanfield_theory=.false.,lOmega_effect=.false.
  logical :: lmeanfield_noalpm=.false.
  logical :: lmeanfield_jxb=.false.,lmeanfield_jxb_with_vA2=.false.
  logical :: lgauss=.false.
  logical :: lbb_as_aux=.false.,ljj_as_aux=.false.
  logical :: lbbt_as_aux=.false.,ljjt_as_aux=.false.
  logical :: lbext_curvilinear=.true., lcheck_positive_va2=.false.
  logical :: lreset_aa=.false.
  character (len=labellen) :: pertaa='zero'
!
  namelist /magnetic_init_pars/ &
       B_ext, lohmic_heat, &
       fring1,Iring1,Rring1,wr1,axisr1,dispr1, &
       fring2,Iring2,Rring2,wr2,axisr2,dispr2, &
       fring3,Iring3,Rring3,wr3,axisr3,dispr3, &
       fring_profile,nrings, &
       radius,epsilonaa,x0aa,z0aa,widthaa, &
       by_left,by_right,bz_left,bz_right, &
       initaa,amplaa,kx_aa,ky_aa,kz_aa,coefaa,coefbb, &
       phasex_aa,phasey_aa,phasez_aa, &
       inclaa,lpress_equil,lpress_equil_via_ss,mu_r, &
       mu_ext_pot,lB_ext_pot,lforce_free_test, &
       ampl_B0,initpower_aa,cutoff_aa,N_modes_aa, &
       rmode,zmode,rm_int,rm_ext,lgauss,lcheck_positive_va2, &
       lbb_as_aux,ljj_as_aux,lbext_curvilinear, &
       lbbt_as_aux,ljjt_as_aux, lneutralion_heat, &
       center1_x,center1_y,center1_z,fluxtube_border_width, &
       va2max_jxb,va2power_jxb
!
  ! run parameters
  real :: eta=0.,eta1=0.,eta_hyper2=0.,eta_hyper3=0.,height_eta=0.,eta_out=0.
  real :: meanfield_molecular_eta=0.
  real :: eta_int=0.,eta_ext=0.,wresistivity=.01
  real :: tau_aa_exterior=0.
  real :: sigma_ratio=1.,eta_width=0.,eta_z0=1.
  real :: alphaSSm=0.
  real :: k1_ff=1.,ampl_ff=1.,swirl=1.
  real :: k1x_ff=1.,k1y_ff=1.,k1z_ff=1.
  real :: inertial_length=0.,linertial_2
  real :: forcing_continuous_aa_phasefact=1.
  real :: forcing_continuous_aa_amplfact=1., ampl_fcont_aa=1.
  real :: LLambda_aa=0.
  real, dimension(nx,my) :: eta_r
  real, dimension(nx,my,3) :: geta_r
  real, dimension(mz) :: coskz,sinkz,eta_z
  real, dimension(mz,3) :: geta_z
  logical :: lfreeze_aint=.false.,lfreeze_aext=.false.
  logical :: lweyl_gauge=.false.
  logical :: lupw_aa=.false.
  logical :: lforcing_cont_aa=.false.
  logical :: lelectron_inertia=.false.
  logical :: lkinematic=.false.
  logical :: luse_Bext_in_b2=.false.
  logical :: lmean_friction=.false.
  character (len=labellen) :: zdep_profile='fs'
  character (len=labellen) :: rdep_profile='schnack89'
  character (len=labellen) :: iforcing_continuous_aa='fixed_swirl'
!
  namelist /magnetic_run_pars/ &
       eta,eta1,eta_hyper2,eta_hyper3,B_ext,omega_Bz_ext, &
       lmeanfield_theory,alpha_effect,alpha_quenching,delta_effect, &
       alpha_eps, &
       alpha_xprofile,alpha_x1,alpha_dx1, &
       alpha_yprofile, &
       lalpha_phi_equation, &
       lOmega_effect,Omega_ampl, &
       Omega_xprofile,Omega_x1,Omega_dx1, &
       Omega_yprofile,Omega_yc1,Omega_yc2, &
       lmeanfield_noalpm, &
       meanfield_etat, lohmic_heat, &
       lmeanfield_jxb,lmeanfield_jxb_with_vA2, &
       meanfield_Qs, meanfield_Qp, &
       meanfield_Bs, meanfield_Bp, &
       meanfield_kf,meanfield_etaB, &
       alpha_equator,alpha_equator_gap,alpha_gap_step,&
       alpha_cutoff_up,alpha_cutoff_down,&
       height_eta,eta_out,tau_aa_exterior, &
       kx_aa,ky_aa,kz_aa,ABC_A,ABC_B,ABC_C, &
       lforcing_cont_aa,iforcing_continuous_aa, &
       forcing_continuous_aa_phasefact, &
       forcing_continuous_aa_amplfact, &
       k1_ff,ampl_ff,swirl,radius, &
       k1x_ff,k1y_ff,k1z_ff,lcheck_positive_va2, &
       lmean_friction,LLambda_aa, &
       bthresh,bthresh_per_brms, &
       iresistivity,lweyl_gauge,lupw_aa, &
       alphaSSm, &
       eta_int,eta_ext,eta_shock,wresistivity, &
       rhomin_jxb,va2max_jxb,va2power_jxb,llorentzforce,linduction, &
       reinitialize_aa,rescale_aa,lB_ext_pot, &
       displacement_gun, &
       pertaa,pertamplaa,D_smag,brms_target,rescaling_fraction, &
       lfreeze_aint,lfreeze_aext, &
       sigma_ratio,zdep_profile,eta_width,eta_z0, &
       borderaa,eta_aniso_hyper3, &
       lbext_curvilinear, &
       lbb_as_aux,ljj_as_aux,lkinematic, &
       lbbt_as_aux,ljjt_as_aux, &
       lneutralion_heat, lreset_aa, daareset, &
       luse_Bext_in_b2, ampl_fcont_aa
!
  ! diagnostic variables (need to be consistent with reset list below)
  integer :: idiag_aphi2m=0     ! DIAG_DOC: $\left<A_\phi^2\right>$
  integer :: idiag_bphi2m=0     ! DIAG_DOC: $\left<B_\phi^2\right>$
  integer :: idiag_bpol2m=0     ! DIAG_DOC: $\left<B_p^2\right>$
  !
  !  plus all the others currently in nomagnetic
  !
  integer :: idiag_b2m=0,idiag_bm2=0,idiag_j2m=0,idiag_jm2=0,idiag_abm=0
  integer :: idiag_jbm=0,idiag_epsM=0,idiag_vArms=0,idiag_vAmax=0
  integer :: idiag_brms=0,idiag_bmax=0,idiag_jrms=0,idiag_jmax=0
  integer :: idiag_bx2m=0, idiag_by2m=0, idiag_bz2m=0,idiag_bmz=0
  integer :: idiag_bxmz=0,idiag_bymz=0,idiag_bzmz=0,idiag_bmx=0,idiag_bmy=0
  integer :: idiag_bxmxy=0,idiag_bymxy=0,idiag_bzmxy=0
  integer :: idiag_uxbm=0,idiag_oxuxbm=0,idiag_jxbxbm=0,idiag_uxDxuxbm=0
  integer :: idiag_b2mphi=0
  integer :: idiag_bmxy_rms=0
  integer :: idiag_bsinphz=0
  integer :: idiag_bcosphz=0
!
  contains
!
!***********************************************************************
    subroutine register_magnetic()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaphi and bphi, etc; increase nvar accordingly
!
!  1-may-02/wolf: coded
!
      use FArrayManager
!
      call farray_register_pde('aphi',iaphi)
      call farray_register_pde('bphi',ibphi)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',aphi,bphi $'
          if (nvar == mvar) write(4,*) ',aphi,bphi'
        else
          write(4,*) ',aphi,bphi $'
        endif
        write(15,*) 'aphi = fltarr(mx,my,mz)*one'
        write(15,*) 'bphi = fltarr(mx,my,mz)*one'
      endif
!
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic(f,lstarting)
!
!  Perform any post-parameter-read initialization
!
!  24-nov-02/tony: dummy routine - nothing to do at present
!  20-may-03/axel: reinitialize_aa added
!
      use BorderProfiles, only: request_border_driving
      use FArrayManager
      use SharedVariables
      use Sub, only: erfunc
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
      integer :: i,ierr
!
!  Precalculate 1/mu (moved here from register.f90)
!
      mu01=1./mu0
      mu012=.5*mu01
!
!  Precalculate eta if 1/eta (==eta1) is given instead
!
      if (eta1/=0.) then
        eta=1./eta1
      endif
!
!  calculate B_ext21 = 1/B_ext**2 and the unit vector B1_ext = B_ext/|B_ext|
!
      B_ext21=B_ext(1)**2+B_ext(2)**2
      if (B_ext21/=0.) then
        B_ext21=1./B_ext21
      else
        B_ext21=1.
      endif
      B_ext11=sqrt(B_ext21)
      B1_ext=B_ext*B_ext11
!
!  rescale magnetic field by a factor reinitialize_aa
!
      if (reinitialize_aa) then
        f(:,:,:,iax:iaz)=rescale_aa*f(:,:,:,iax:iaz)
      endif
!
!  set lrescaling_magnetic=T if linit_aa=T
!
      if (lreset_aa) then
        lrescaling_magnetic=.true.
      endif
!
      if (lfreeze_aint) lfreeze_varint(iax:iaz) = .true.
      if (lfreeze_aext) lfreeze_varext(iax:iaz) = .true.
!
!  Initialize resistivity.
!
      if (iresistivity(1)=='') iresistivity(1)='eta-const'  ! default
      lresi_eta_const=.false.
      lresi_hyper2=.false.
      lresi_hyper3=.false.
      lresi_hyper3_polar=.false.
      lresi_hyper3_strict=.false.
      lresi_hyper3_aniso=.false.
      lresi_eta_shock=.false.
      lresi_eta_shock_perp=.false.
      lresi_smagorinsky=.false.
      lresi_smagorinsky_cross=.false.
!
      do i=1,nresi_max
        select case (iresistivity(i))
        case ('eta-const')
          if (lroot) print*, 'resistivity: constant eta'
          lresi_eta_const=.true.
        case ('etaSS')
          if (lroot) print*, 'resistivity: etaSS (Shakura-Sunyaev)'
          lresi_etaSS=.true.
        case ('hyper2')
          if (lroot) print*, 'resistivity: hyper2'
          lresi_hyper2=.true.
        case ('hyper3')
          if (lroot) print*, 'resistivity: hyper3'
          lresi_hyper3=.true.
        case ('hyper3_cyl','hyper3-cyl','hyper3_sph','hyper3-sph')
          if (lroot) print*, 'resistivity: hyper3 curvilinear'
          lresi_hyper3_polar=.true.
        case ('hyper3_strict')
          if (lroot) print*, 'resistivity: strict hyper3 with positive definite heating rate'
          lresi_hyper3_strict=.true.
        case ('rdep')
          if (lroot) print*, 'resistivity: r-dependent'
          lresi_rdep=.true.
          call eta_rdep(eta_r,geta_r,rdep_profile)
        case ('zdep')
          if (lroot) print*, 'resistivity: z-dependent'
          lresi_zdep=.true.
          call eta_zdep(eta_z,geta_z,zdep_profile)
        case ('dust')
          if (lroot) print*, 'resistivity: depending on dust density'
          lresi_dust=.true.
        case ('hyper3-aniso')
          if (lroot) print*, 'resistivity: hyper3_aniso'
          lresi_hyper3_aniso=.true.
        case ('shell')
          if (lroot) print*, 'resistivity: shell'
          lresi_shell=.true.
        case ('shock','eta-shock')
          if (lroot) print*, 'resistivity: shock'
          lresi_eta_shock=.true.
          if (.not. lshock) &
              call fatal_error('initialize_magnetic', &
              'shock resistivity, but module setting SHOCK=noshock')
        case ('shock-perp')
          if (lroot) print*, 'resistivity: shock_perp'
          lresi_eta_shock_perp=.true.
          if (.not. lshock) &
              call fatal_error('initialize_magnetic', &
              'shock resistivity, but module setting SHOCK=noshock')
        case ('smagorinsky')
          if (lroot) print*, 'resistivity: smagorinsky'
          lresi_smagorinsky=.true.
        case ('smagorinsky-cross')
          if (lroot) print*, 'resistivity: smagorinsky_cross'
          lresi_smagorinsky_cross=.true.
        case ('none')
          ! do nothing
        case ('')
          ! do nothing
        case default
          if (lroot) print*, 'No such such value for iresistivity(',i,'): ', &
              trim(iresistivity(i))
          call fatal_error('initialize_magnetic','')
        endselect
      enddo
!
!  If we're timestepping, die or warn if the the resistivity coefficient that
!  corresponds to the chosen resistivity type is not set.
!
      if (lrun) then
        if (lresi_eta_const.and.(eta==0.0.and.meanfield_etat==0.0)) &
            call warning('initialize_magnetic', &
            'Resistivity coefficient eta is zero!')
        if (lresi_hyper2.and.eta_hyper2==0.0) &
            call fatal_error('initialize_magnetic', &
            'Resistivity coefficient eta_hyper2 is zero!')
        if (lresi_hyper3.and.eta_hyper3==0.0) &
            call fatal_error('initialize_magnetic', &
            'Resistivity coefficient eta_hyper3 is zero!')
        if (lresi_hyper3_polar.and.eta_hyper3==0.0) &
             call fatal_error('initialize_magnetic', &
            'Resistivity coefficient eta_hyper3 is zero!')
        if (lresi_hyper3_strict.and.eta_hyper3==0.0) &
            call fatal_error('initialize_magnetic', &
            'Resistivity coefficient eta_hyper3 is zero!')
        if ( (lresi_hyper3_aniso) .and.  &
             ((eta_aniso_hyper3(1)==0. .and. nxgrid/=1 ).or. &
              (eta_aniso_hyper3(2)==0. .and. nygrid/=1 ).or. &
              (eta_aniso_hyper3(3)==0. .and. nzgrid/=1 )) ) &
            call fatal_error('initialize_magnetic', &
            'A resistivity coefficient of eta_aniso_hyper3 is zero!')
        if (lresi_eta_shock.and.eta_shock==0.0) &
            call fatal_error('initialize_magnetic', &
            'Resistivity coefficient eta_shock is zero!')
        if (lresi_eta_shock_perp.and.eta_shock==0.0) &
            call fatal_error('initialize_magnetic', &
            'Resistivity coefficient eta_shock is zero!')
!        if (lentropy .and. lohmic_heat .and. .not. lresi_eta_const) &
!            call fatal_error('initialize_magnetic', &
!            'Resistivity heating only works with regular resistivity!')
        if (lresi_hyper2.and.lresi_hyper3) &
            call warning('initialize_magnetic', &
            '4th & 6th order hyperdiffusion are both set. ' // &
            'Timestep is currently only sensitive to fourth order.')
      endif
!
!  alpha profile
!  =============
!
      select case (alpha_xprofile)
      case ('const'); alpha_x=1.; dalpha_x=0.
      case ('erfx'); alpha_x=.5*(1.+erfunc((x(l1:l2)-alpha_x1)/alpha_dx1))
        dalpha_x=exp(-((x(l1:l2)-alpha_x1)/alpha_dx1)**2)/(alpha_dx1*sqrtpi)
      case default
        if (lroot) print*,'No such such value for alpha_yprofile'
        call fatal_error('initialize_magnetic','')
      endselect
!
!  y direction
!
      select case (alpha_yprofile)
      case ('const'); alpha_y=1.; dalpha_y=0.
      case ('cosy'); alpha_y=cos(y); dalpha_y=-sin(y)
      case ('cosy*sin2y'); alpha_y=1.5*sqrt(3.)*cos(y)*sin(y)**2
        dalpha_y=1.5*sqrt(3.)*(2.-3.*sin(y)**2)*sin(y)
      case ('read')
        print*,'read alpha profile'
        open(1,file='alpha_input.dat',form='unformatted')
        read(1) alpha_y
!--     read(1) alpha_input
        close(1)
      case default
        if (lroot) print*,'No such such value for alpha_yprofile'
        call fatal_error('initialize_magnetic','')
      endselect
!
!  Omega effect and its gradient
!  =============================
!
      if (lOmega_effect) then
        select case (Omega_xprofile)
        case ('const'); Omega_x=1.; dOmega_x=0.
        case ('Omega=x'); Omega_x=x(l1:l2); dOmega_x=1.
        case ('erfx'); Omega_x=.5*(1.+erfunc((x(l1:l2)-Omega_x1)/Omega_dx1))
          dOmega_x=exp(-((x(l1:l2)-Omega_x1)/Omega_dx1)**2)/(Omega_dx1*sqrtpi)
        case default
          if (lroot) print*,'No such such value for Omega_xprofile'
          call fatal_error('initialize_magnetic','')
        endselect
!
!  y direction; compute Omega_y and dOmega_y = (dOmega_y/dy)
!
        select case (Omega_yprofile)
        case ('const'); Omega_y=1.; dOmega_y=0.
        case ('c1+c2*cos2y'); Omega_y=Omega_yc1+Omega_yc2*cos(y)**2
          dOmega_y=-2.*Omega_yc2*cos(y)*sin(y)
        case default
          if (lroot) print*,'No such such value for Omega_yprofile'
          call fatal_error('initialize_magnetic','')
        endselect
      endif
!
!  if meanfield theory is invoked, we want to send meanfield_etat to
!  other subroutines
!
      call put_shared_variable('lmeanfield_theory',lmeanfield_theory,ierr)
      if (lmeanfield_theory) then
        call put_shared_variable('meanfield_etat',meanfield_etat,ierr)
        call put_shared_variable('eta',eta,ierr)
      endif
!
!  Tell the BorderProfiles module if we intend to use border driving, so
!  that the module can request the right pencils.
!
      if (borderaa/='nothing') call request_border_driving(borderaa)
!
!  Register an extra aux slot for bb if requested (so bb and jj are written
!  to snapshots and can be easily analyzed later). For this to work you
!  must reserve enough auxiliary workspace by setting, for example,
!     ! MAUX CONTRIBUTION 6
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
!
!  After a reload, we need to rewrite index.pro, but the auxiliary
!  arrays are already allocated and must not be allocated again.
!
      if (lbb_as_aux) then
        if (ibb==0) then
          call farray_register_auxiliary('bb',ibb,vector=3)
          ibx=ibb
          iby=ibb+1
          ibz=ibb+2
        endif
        if (ibb/=0.and.lroot) then
          print*, 'initialize_magnetic: ibb = ', ibb
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'ibb=',ibb
          write(3,*) 'ibx=',ibx
          write(3,*) 'iby=',iby
          write(3,*) 'ibz=',ibz
          close(3)
        endif
      endif
!
!  do the same for jj (current density)
!
      if (ljj_as_aux) then
        if (ijj==0) then
          call farray_register_auxiliary('jj',ijj,vector=3)
          ijx=ijj
          ijy=ijj+1
          ijz=ijj+2
        endif
        if (ijj/=0.and.lroot) then
          print*, 'initialize_magnetic: ijj = ', ijj
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'ijj=',ijj
          write(3,*) 'ijx=',ijx
          write(3,*) 'ijy=',ijy
          write(3,*) 'ijz=',ijz
          close(3)
        endif
      endif
!
!  After a reload, we need to rewrite index.pro, but the auxiliary
!  arrays are already allocated and must not be allocated again.
!
      if (lbbt_as_aux) then
        if (ibbt==0) then
          call farray_register_auxiliary('bbt',ibbt,vector=3)
          ibxt=ibbt
          ibyt=ibbt+1
          ibzt=ibbt+2
        endif
        if (ibbt/=0.and.lroot) then
          print*, 'initialize_velocity: ibbt = ', ibbt
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'ibbt=',ibbt
          write(3,*) 'ibxt=',ibxt
          write(3,*) 'ibyt=',ibyt
          write(3,*) 'ibzt=',ibzt
          close(3)
        endif
      endif
!
      if (ljjt_as_aux) then
        if (ijjt==0) then
          call farray_register_auxiliary('jjt',ijjt,vector=3)
          ijxt=ijjt
          ijyt=ijjt+1
          ijzt=ijjt+2
        endif
        if (ijjt/=0.and.lroot) then
          print*, 'initialize_velocity: ijjt = ', ijjt
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'ijjt=',ijjt
          write(3,*) 'ijxt=',ijxt
          write(3,*) 'ijyt=',ijyt
          write(3,*) 'ijzt=',ijzt
          close(3)
        endif
      endif
!
      if (any(initaa=='Alfven-zconst')) then
        call put_shared_variable('zmode',zmode,ierr)
        if (ierr/=0) call fatal_error('initialize_magnetic',&
             'there was a problem when sharing zmode')
      endif
!
      call put_shared_variable('lfrozen_bb_bot',lfrozen_bb_bot,ierr)
      if (ierr/=0) call fatal_error('initialize_magnetic',&
           'there was a problem when sharing lfrozen_bb_bot')
      call put_shared_variable('lfrozen_bb_top',lfrozen_bb_top,ierr)
      if (ierr/=0) call fatal_error('initialize_magnetic',&
           'there was a problem when sharing lfrozen_bb_top')
!
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_magnetic
!***********************************************************************
    subroutine init_aa(f)
!
!  initialise magnetic field; called from start.f90
!  AB: maybe we should here call different routines (such as rings)
!  AB: and others, instead of accummulating all this in a huge routine.
!  We have an init parameter (initaa) to stear magnetic i.c. independently.
!
!   7-nov-2001/wolf: coded
!
      use EquationOfState
      use FArrayManager
      use Gravity, only: gravz, z1, z2
      use Initcond
      use InitialCondition, only: initial_condition_aa
      use Mpicomm
      use SharedVariables
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (mz) :: tmp
      real, dimension (nx,3) :: bb
      real, dimension (nx) :: b2,fact
      real :: beq2
      integer :: j
!
      do j=1,ninit
!
        select case (initaa(j))
        case ('nothing'); if (lroot .and. j==1) print*,'init_aa: nothing'
        case ('zero', '0'); f(:,:,:,iaphi) = 0.; f(:,:,:,ibphi) = 0.
        case ('gaussian-noise'); call gaunoise(amplaa(j),f,iaphi); call gaunoise(amplaa(j),f,ibphi)
        case ('dipolar-field'); f(:,:,:,ibphi) = 0.; call phi_siny_over_r2(f,iaphi)
!
        case default
!
!  Catch unknown values
!
          call fatal_error('init_aa', &
              'init_aa value "' // trim(initaa(j)) // '" not recognised')
!
        endselect
!
!  End loop over initial conditions
!
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_aa(f)
!
!  allow for pressure equilibrium (for isothermal tube)
!  assume that ghost zones have already been set.
!  corrected expression below for gamma /= 1 case.
!  The beq2 expression for 2*mu0*p is not general yet.
!
      if (lpress_equil.or.lpress_equil_via_ss) then
        if (lroot) print*,'init_aa: adjust lnrho to have pressure equilib; cs0=',cs0
        do n=n1,n2
        do m=m1,m2
          call curl(f,iaa,bb)
          call dot2_mn(bb,b2)
          if (gamma==1.) then
            f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)-b2/(2.*cs0**2)
          else
            beq2=2.*rho0*cs0**2
            fact=max(1e-6,1.-b2/beq2)
            if (lentropy.and.lpress_equil_via_ss) then
              f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss)+fact/gamma
            else
              f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)+fact/gamma_m1
            endif
          endif
        enddo
        enddo
      endif
!
    endsubroutine init_aa
!***********************************************************************
    subroutine pencil_criteria_magnetic()
!
!   All pencils that the Magnetic module depends on are specified here.
!
!  19-nov-04/anders: coded
!
      use Mpicomm, only: stop_it
!
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_uxb)=.true.
!
      if (dvid/=0.0) lpenc_video(i_b2)=.true.
!
!  jj pencil always needed when in Weyl gauge
!
      if ( height_eta/=0. .or. ip<=4 .or. &
          (lweyl_gauge) .or. (lspherical_coords) ) &
!  WL: but doesn't seem to be needed for the cylindrical case
          lpenc_requested(i_jj)=.true.
      if ((eta/=0..or.meanfield_etat/=0.).and. &
          (.not.lweyl_gauge)) lpenc_requested(i_del2a)=.true.
      if (dvid/=0.) lpenc_video(i_jb)=.true.
      if (lresi_eta_const .or. lresi_shell .or. &
          lresi_eta_shock .or. lresi_smagorinsky .or. &
          lresi_zdep .or. lresi_rdep .or. &
          lresi_smagorinsky_cross) lpenc_requested(i_del2a)=.true.
      if (lresi_eta_shock) then
        lpenc_requested(i_shock)=.true.
        if (lweyl_gauge) then
          lpenc_requested(i_jj)=.true.
        else
          lpenc_requested(i_gshock)=.true.
          lpenc_requested(i_diva)=.true.
        endif
      endif
      if (lresi_shell) then
        lpenc_requested(i_r_mn)=.true.
        lpenc_requested(i_evr)=.true.
      endif
      if (lresi_eta_shock_perp) then
        lpenc_requested(i_shock_perp)=.true.
        if (lweyl_gauge) then
          lpenc_requested(i_jj)=.true.
        else
          lpenc_requested(i_gshock_perp)=.true.
          lpenc_requested(i_diva)=.true.
        endif
      endif
      if (lupw_aa) then
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_aij)=.true.
      endif
      if (lresi_shell .or. lresi_zdep .or. lresi_rdep) lpenc_requested(i_diva)=.true.
      if (lresi_smagorinsky_cross) lpenc_requested(i_jo)=.true.
      if (lresi_hyper2) lpenc_requested(i_del4a)=.true.
      if (lresi_hyper3) lpenc_requested(i_del6a)=.true.
!WL: for the cylindrical case, lpencil_check says graddiva is not needed
      if (lspherical_coords) lpenc_requested(i_graddiva)=.true.
      if (lentropy .or. lresi_smagorinsky .or. ltemperature) then
        lpenc_requested(i_j2)=.true.
      endif
      if (lresi_dust) lpenc_requested(i_rhop)=.true.
      if (lentropy .or. ltemperature .or. ldt) lpenc_requested(i_rho1)=.true.
      if (lentropy .or. ltemperature) lpenc_requested(i_TT1)=.true.
      if (ltemperature) lpenc_requested(i_cv1)=.true.
!
!  for mean-field modelling
!
      if (lmeanfield_jxb) then
        lpenc_requested(i_va2)=.true.
        lpenc_requested(i_rho1)=.true.
      endif
!
!  ambipolar diffusion
!
      if (lhydro .and. llorentzforce) &
          lpenc_requested(i_jxbr)=.true.
      if (lresi_smagorinsky_cross .or. delta_effect/=0.) &
          lpenc_requested(i_oo)=.true.
      if (lmeanfield_theory) then
        if (meanfield_etat/=0. .or. alpha_effect/=0. .or. delta_effect/=0.) &
            lpenc_requested(i_mf_EMF)=.true.
        if (delta_effect/=0.) lpenc_requested(i_oxj)=.true.
      endif
!
    endsubroutine pencil_criteria_magnetic
!***********************************************************************
    subroutine pencil_interdep_magnetic(lpencil_in)
!
!  Interdependency among pencils from the Magnetic module is specified here.
!
!  19-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_jparallel).or.lpencil_in(i_jperp)) then
        lpencil_in(i_cosjb)=.true.
        lpencil_in(i_jxb)=.true.
      endif
      if (lpencil_in(i_cosjb)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_j2)=.true.
        lpencil_in(i_jb)=.true.
      endif
      if (lpencil_in(i_a2)) lpencil_in(i_aa)=.true.
      if (lpencil_in(i_ab)) then
        lpencil_in(i_aa)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_va2)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_rho1)=.true.
      endif
      if (lpencil_in(i_j2)) lpencil_in(i_jj)=.true.
      if (lpencil_in(i_uxj)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_jj)=.true.
      endif
      if (lpencil_in(i_jb)) then
        lpencil_in(i_bb)=.true.
        lpencil_in(i_jj)=.true.
      endif
      if (lpencil_in(i_jxbr) .and. va2max_jxb>0) lpencil_in(i_va2)=.true.
      if (lpencil_in(i_jxbr)) then
        lpencil_in(i_jxb)=.true.
        lpencil_in(i_rho1)=.true.
      endif
      if (lpencil_in(i_jxb)) then
        lpencil_in(i_jj)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_uxb2)) lpencil_in(i_uxb)=.true.
      if (lpencil_in(i_uxb)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_cosub)) then
        lpencil_in(i_ub)=.true.
        lpencil_in(i_u2)=.true.
        lpencil_in(i_b2)=.true.
      endif
      if (lpencil_in(i_ub)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_beta1)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_pp)=.true.
      endif
      if (lpencil_in(i_b2)) lpencil_in(i_bb)=.true.
      if (lpencil_in(i_jj)) lpencil_in(i_bij)=.true.
      if (lpencil_in(i_bb)) then
        if (.not.lcartesian_coords) lpencil_in(i_aa)=.true.
        lpencil_in(i_aij)=.true.
      endif
      if (lpencil_in(i_djuidjbi)) then
        lpencil_in(i_uij)=.true.
        lpencil_in(i_bij)=.true.
      endif
      if (lpencil_in(i_jo)) then
        lpencil_in(i_oo)=.true.
        lpencil_in(i_jj)=.true.
      endif
      if (lpencil_in(i_ujxb)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_jxb)=.true.
      endif
      if (lpencil_in(i_oxu)) then
        lpencil_in(i_oo)=.true.
        lpencil_in(i_uu)=.true.
      endif
      if (lpencil_in(i_oxuxb)) then
        lpencil_in(i_oxu)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_jxbxb)) then
        lpencil_in(i_jxb)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_jxbrxb)) then
        lpencil_in(i_jxbr)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_glnrhoxb)) then
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_oxj)) then
        lpencil_in(i_oo)=.true.
        lpencil_in(i_jj)=.true.
      endif
      if (lpencil_in(i_jij)) lpencil_in(i_bij)=.true.
      if (lpencil_in(i_sj)) then
        lpencil_in(i_sij)=.true.
        lpencil_in(i_jij)=.true.
      endif
      if (lpencil_in(i_ss12)) lpencil_in(i_sij)=.true.
      if (lpencil_in(i_mf_EMFdotB)) then
        lpencil_in(i_mf_EMF)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_mf_EMF)) then
        if (lspherical_coords) then
          lpencil_in(i_jj)=.true.
          lpencil_in(i_graddivA)=.true.
        endif
        lpencil_in(i_b2)=.true.
        lpencil_in(i_bb)=.true.
        if (delta_effect/=0.) lpencil_in(i_oxJ)=.true.
        if (meanfield_etat/=0.) then
          if (lweyl_gauge) then
            lpencil_in(i_jj)=.true.
          else
            lpencil_in(i_del2a)=.true.
          endif
        endif
      endif
      if (lpencil_in(i_del2A)) then
        if (lspherical_coords) then
!WL: for the cylindrical case, lpencil_check says these pencils are not needed
          lpencil_in(i_jj)=.true.
          lpencil_in(i_graddivA)=.true.
        endif
      endif
      if (lpencil_in(i_uga)) then
        lpencil_in(i_aij)=.true.
        lpencil_in(i_uu)=.true.
      endif
!
      if (lpencil_in(i_ss12)) lpencil_in(i_sj)=.true.
!  Pencils bij, del2a and graddiva come in a bundle.
!     if ( lpencil_in(i_bij) .and. &
!         (lpencil_in(i_del2a).or.lpencil_in(i_graddiva)) ) then
!       lpencil_in(i_del2a)=.false.
!       lpencil_in(i_graddiva)=.false.
!     endif
!     if (lpencil_in(i_del2a) .and. &
!         (lpencil_in(i_bij).or.lpencil_in(i_graddiva)) ) then
!       lpencil_in(i_bij)=.false.
!       lpencil_in(i_graddiva)=.false.
!     endif
!
    endsubroutine pencil_interdep_magnetic
!***********************************************************************
    subroutine calc_pencils_magnetic(f,p)
!
!  Calculate Magnetic pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  19-nov-04/anders: coded
!
      use Sub
      use Deriv
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: bb_ext,bb_ext_pot
      real, dimension (nx) :: rho1_jxb,alpha_total
      real, dimension (nx) :: alpha_tmp
      real, dimension (nx) :: jcrossb2
      real, dimension (nx) :: meanfield_Qs_func, meanfield_Qp_func
      real, dimension (nx) :: meanfield_Qs_der, meanfield_Qp_der, BiBk_Bki
      real, dimension (nx) :: meanfield_Bs21, meanfield_Bp21
      real, dimension (nx) :: meanfield_Bs2, meanfield_Bp2
      real, dimension (nx) :: meanfield_urms21,meanfield_etaB2
      real, dimension (nx,3) :: Bk_Bki
      real :: B2_ext,c,s,kx
      integer :: i,j,ix
!
      intent(inout) :: f,p
! aa
      if (lpencil(i_aa)) p%aa=f(l1:l2,m,n,iax:iaz)
! a2
      if (lpencil(i_a2)) call dot2_mn(p%aa,p%a2)
! aij
      if (lpencil(i_aij)) call gij(f,iaa,p%aij,1)
! diva
      if (lpencil(i_diva)) call div_mn(p%aij,p%diva,p%aa)
! bb
      if (lpencil(i_bb)) then
        call curl_mn(p%aij,p%bb,p%aa)
!
!  save field before adding imposed field (for diagnostics)
!
        p%bbb=p%bb
        B2_ext=B_ext(1)**2+B_ext(2)**2+B_ext(3)**2
!
!  allow external field to precess about z-axis
!  with frequency omega_Bz_ext
!
        if (B2_ext/=0.) then
          if (lbext_curvilinear.or.lcartesian_coords) then
!
!  luse_curvilinear_bext is default. The B_ext the user defines in
!  magnetic_init_pars respects the coordinate system of preference
!  which means that B_ext=(0.,1.,0.) is an azimuthal field in cylindrical
!  coordinates and a polar one in spherical.
!
            if (omega_Bz_ext==0.) then
              B_ext_tmp=B_ext
            elseif (omega_Bz_ext/=0.) then
              c=cos(omega_Bz_ext*t)
              s=sin(omega_Bz_ext*t)
              B_ext_tmp(1)=B_ext(1)*c-B_ext(2)*s
              B_ext_tmp(2)=B_ext(1)*s+B_ext(2)*c
              B_ext_tmp(3)=B_ext(3)
            endif
          else if (lcylindrical_coords) then
            if (omega_Bz_ext/=0.) &
                 call fatal_error("calc_pencils_magnetic",&
                  "precession of the external field not "//&
                  "implemented for cylindrical coordinates")
!
!  transform b_ext to other coordinate systems
!
            B_ext_tmp(1)=  B_ext(1)*cos(y(m)) + B_ext(2)*sin(y(m))
            B_ext_tmp(2)= -B_ext(1)*sin(y(m)) + B_ext(2)*cos(y(m))
            B_ext_tmp(3)=  B_ext(3)
          else if (lspherical_coords) then
            if (omega_Bz_ext/=0.) &
                 call fatal_error("calc_pencils_magnetic",&
                  "precession of the external field not "//&
                  "implemented for spherical coordinates")
            B_ext_tmp(1)= B_ext(1)*sinth(m)*cos(z(n)) + B_ext(2)*sinth(m)*sin(z(n)) + B_ext(3)*costh(m)
            B_ext_tmp(2)= B_ext(1)*costh(m)*cos(z(n)) + B_ext(2)*costh(m)*sin(z(n)) - B_ext(3)*sinth(m)
            B_ext_tmp(3)=-B_ext(1)         *sin(z(n)) + B_ext(2)         *cos(z(n))
          endif
!
!  add the external field
!
          if (B_ext_tmp(1)/=0.) p%bb(:,1)=p%bb(:,1)+B_ext_tmp(1)
          if (B_ext_tmp(2)/=0.) p%bb(:,2)=p%bb(:,2)+B_ext_tmp(2)
          if (B_ext_tmp(3)/=0.) p%bb(:,3)=p%bb(:,3)+B_ext_tmp(3)
          if (headtt) print*,'calc_pencils_magnetic: B_ext=',B_ext
          if (headtt) print*,'calc_pencils_magnetic: B_ext_tmp=',B_ext_tmp
        endif
!
!  add the external potential field
!
        if (lB_ext_pot) then
!          call get_global(bb_ext_pot,m,n,'B_ext_pot')
!          p%bb=p%bb+bb_ext_pot
        endif
!
!  add external B-field.
!
        if (iglobal_bx_ext/=0) p%bb(:,1)=p%bb(:,1)+f(l1:l2,m,n,iglobal_bx_ext)
        if (iglobal_by_ext/=0) p%bb(:,2)=p%bb(:,2)+f(l1:l2,m,n,iglobal_by_ext)
        if (iglobal_bz_ext/=0) p%bb(:,3)=p%bb(:,3)+f(l1:l2,m,n,iglobal_bz_ext)
      endif
!
! b2 (default is that B_ext is not included), but this can be changed
! by setting luse_Bext_in_b2=.true.
!
      if (luse_Bext_in_b2) then
        if (lpencil(i_b2)) call dot2_mn(p%bb,p%b2)
      else
        if (lpencil(i_b2)) call dot2_mn(p%bbb,p%b2)
      endif
! ab
      if (lpencil(i_ab)) call dot_mn(p%aa,p%bbb,p%ab)
! uxb
      if (lpencil(i_uxb)) then
        call cross_mn(p%uu,p%bb,p%uxb)
!  add external e-field.
        if (iglobal_ex_ext/=0) p%uxb(:,1)=p%uxb(:,1)+f(l1:l2,m,n,iglobal_ex_ext)
        if (iglobal_ey_ext/=0) p%uxb(:,2)=p%uxb(:,2)+f(l1:l2,m,n,iglobal_ey_ext)
        if (iglobal_ez_ext/=0) p%uxb(:,3)=p%uxb(:,3)+f(l1:l2,m,n,iglobal_ez_ext)
      endif
! uga
! DM : this requires later attention
      if (lpencil(i_uga)) then
        if (.not.lcartesian_coords) then
          call warning("calc_pencils_magnetic","u_dot_grad A not implemented for non-cartesian coordinates")
        else
          call u_dot_grad(f,iaa,p%aij,p%uu,p%uga,UPWIND=lupw_aa)
        endif
      endif
!
!  bij, del2a, graddiva
!  For non-cartesian coordinates jj is always required for del2a=graddiva-jj
!
      if (lpencil(i_bij) .or. lpencil(i_del2a) .or. lpencil(i_graddiva) .or. &
          lpencil(i_jj) ) then
        if (lcartesian_coords) then
          call gij_etc(f,iaa,p%aa,p%aij,p%bij,p%del2a,p%graddiva)
          if (.not. lpencil(i_bij)) p%bij=0.0      ! Avoid warnings from pencil
          if (.not. lpencil(i_del2A)) p%del2A=0.0  ! consistency check...
          if (.not. lpencil(i_graddiva)) p%graddiva=0.0
          if (lpencil(i_jj)) call curl_mn(p%bij,p%jj,p%bb)
        else
          call gij_etc(f,iaa,p%aa,p%aij,p%bij,GRADDIV=p%graddiva)
          call curl_mn(p%bij,p%jj,p%bb)
          if (lpencil(i_del2a)) p%del2a=p%graddiva-p%jj
!           if (lpencil(i_del2a)) call del2v(f,iaa,p%del2a,p%aij,p%aa)
        endif
      endif
! jj
      if (lpencil(i_jj)) then
        p%jj=mu01*p%jj
!
!  add external j-field.
!
        if (iglobal_jx_ext/=0) p%jj(:,1)=p%jj(:,1)+f(l1:l2,m,n,iglobal_jx_ext)
        if (iglobal_jy_ext/=0) p%jj(:,2)=p%jj(:,2)+f(l1:l2,m,n,iglobal_jy_ext)
        if (iglobal_jz_ext/=0) p%jj(:,3)=p%jj(:,3)+f(l1:l2,m,n,iglobal_jz_ext)
      endif
! j2
      if (lpencil(i_j2)) call dot2_mn(p%jj,p%j2)
! jb
      if (lpencil(i_jb)) call dot_mn(p%jj,p%bbb,p%jb)
! va2
      if (lpencil(i_va2)) then
        p%va2=p%b2*mu01*p%rho1
        if (lcheck_positive_va2 .and. minval(p%va2)<0.0) then
          print*, 'calc_pencils_magnetic: Alfven speed is imaginary!'
          print*, 'calc_pencils_magnetic: it, itsub, iproc=', it, itsub, iproc
          print*, 'calc_pencils_magnetic: m, y(m), n, z(n)=', m, y(m), n, z(n)
          p%va2=abs(p%va2)
        endif
      endif
! jxb
      if (lpencil(i_jxb)) call cross_mn(p%jj,p%bb,p%jxb)
! jxbr
      if (lpencil(i_jxbr)) rho1_jxb=p%rho1
! cosjb
      if (lpencil(i_cosjb)) then
        do ix=1,nx
          if ((abs(p%j2(ix))<=tini).or.(abs(p%b2(ix))<=tini))then
            p%cosjb(ix)=0.
          else
            p%cosjb(ix)=p%jb(ix)/sqrt(p%j2(ix)*p%b2(ix))
          endif
        enddo
        if (lpencil_check) then
        ! map penc0 value back to interval [-1,1]
          p%cosjb = modulo(p%cosjb + 1.0, 2.0) - 1
        endif
      endif
! jparallel and jperp
      if (lpencil(i_jparallel).or.lpencil(i_jperp)) then
        p%jparallel=sqrt(p%j2)*p%cosjb
        call dot2_mn(p%jxb,jcrossb2)
        do ix=1,nx
          if ((abs(p%j2(ix))<=tini).or.(abs(p%b2(ix))<=tini))then
            p%jperp=0
          else
            p%jperp=sqrt(jcrossb2(ix))/sqrt(p%b2(ix))
          endif
        enddo
!        if (lpencil_check) then
!          sinjb=sqrt(1-(modulo(p%cosjb + 1.0, 2.0) - 1)**2)
!        else
!          sinjb=sqrt(1-p%cosjb**2)
!        endif
!        p%jperp=sqrt(p%j2)*sinjb
      endif
! jxbr
      if (lpencil(i_jxbr)) then
        rho1_jxb=p%rho1
!
!  Set rhomin_jxb>0 in order to limit the jxb term at very low densities.
!  Set va2max_jxb>0 in order to limit the jxb term at very high Alfven speeds.
!  Set va2power_jxb to an integer value in order to specify the power of the
!  limiting term,
!
        if (rhomin_jxb>0) rho1_jxb=min(rho1_jxb,1/rhomin_jxb)
        if (va2max_jxb>0) then
          rho1_jxb = rho1_jxb &
                   * (1+(p%va2/va2max_jxb)**va2power_jxb)**(-1.0/va2power_jxb)
        endif
        if (lmeanfield_jxb) then
          if (lmeanfield_jxb_with_vA2) then
            meanfield_urms21=1./(3.*meanfield_kf*meanfield_etat)**2
            meanfield_Qs_func=1.+(meanfield_Qs-1.)*(1.-2*pi_1*atan(p%vA2*meanfield_urms21))
            meanfield_Qp_func=1.+(meanfield_Qp-1.)*(1.-2*pi_1*atan(p%vA2*meanfield_urms21))
            meanfield_Qs_der=-2*pi_1*(meanfield_Qs-1.)/(1.+(p%vA2*meanfield_urms21)**2)
            meanfield_Qp_der=-2*pi_1*(meanfield_Qp-1.)/(1.+(p%vA2*meanfield_urms21)**2)
            call multsv_mn(meanfield_Qs_func,p%jxb,p%jxb)
!           call multmv_transp(p%bij,p%bb,Bk_Bki)
            !call multsv_mn_add(meanfield_Qs_func-meanfield_Qp_func-p%b2*meanfield_Qp_der,Bk_Bki,p%jxb)
            !call multsv_mn_add(meanfield_Qs_func-meanfield_Qp_func-p%vA2*meanfield_urms21/meanfield_Bp2*meanfield_Qp_der,Bk_Bki,p%jxb)
!           call multsv_mn_add(meanfield_Qs_func-meanfield_Qp_func,Bk_Bki,p%jxb)
            !call dot(Bk_Bki,p%bb,BiBk_Bki)
!           call multsv_mn_add(2*meanfield_Qp_der*BiBk_Bki*p%rho1*meanfield_urms21,p%bb,p%jxb)
          else
            meanfield_Bs21=1./meanfield_Bs**2
            meanfield_Bp21=1./meanfield_Bp**2
            meanfield_Qs_func=1.+(meanfield_Qs-1.)*(1.-2*pi_1*atan(p%b2*meanfield_Bs21))
            meanfield_Qp_func=1.+(meanfield_Qp-1.)*(1.-2*pi_1*atan(p%b2*meanfield_Bp21))
            meanfield_Qs_der=-2*pi_1*(meanfield_Qs-1.)*meanfield_Bs21/(1.+(p%b2*meanfield_Bs21)**2)
            meanfield_Qp_der=-2*pi_1*(meanfield_Qp-1.)*meanfield_Bp21/(1.+(p%b2*meanfield_Bp21)**2)
            call multsv_mn(meanfield_Qs_func,p%jxb,p%jxb)
            call multmv_transp(p%bij,p%bb,Bk_Bki)
            call multsv_mn_add(meanfield_Qs_func-meanfield_Qp_func-p%b2*meanfield_Qp_der,Bk_Bki,p%jxb)
            call dot(Bk_Bki,p%bb,BiBk_Bki)
            call multsv_mn_add(2*meanfield_Qs_der*BiBk_Bki,p%bb,p%jxb)
          endif
        endif
        call multsv_mn(rho1_jxb,p%jxb,p%jxbr)
      endif
! jxbr2
      if (lpencil(i_jxbr2)) call dot2_mn(p%jxbr,p%jxbr2)
! ub
      if (lpencil(i_ub)) call dot_mn(p%uu,p%bb,p%ub)
! cosub
      if (lpencil(i_cosub)) then
        do ix=1,nx
          if ((abs(p%u2(ix))<=tini).or.(abs(p%b2(ix))<=tini)) then
            p%cosub(ix)=0.
          else
            p%cosub(ix)=p%ub(ix)/sqrt(p%u2(ix)*p%b2(ix))
          endif
        enddo
        if (lpencil_check) then
        ! map penc0 value back to interval [-1,1]
          p%cosub = modulo(p%cosub + 1.0, 2.0) - 1
        endif
      endif
! uxb2
      if (lpencil(i_uxb2)) call dot2_mn(p%uxb,p%uxb2)
! uxj
      if (lpencil(i_uxj)) call cross_mn(p%uu,p%jj,p%uxj)
! beta1
      if (lpencil(i_beta1)) p%beta1=0.5*p%b2*mu01/p%pp
! djuidjbi
      if (lpencil(i_djuidjbi)) call multmm_sc(p%uij,p%bij,p%djuidjbi)
! jo
      if (lpencil(i_jo)) call dot(p%jj,p%oo,p%jo)
! ujxb
      if (lpencil(i_ujxb)) call dot_mn(p%uu,p%jxb,p%ujxb)
! oxu
      if (lpencil(i_oxu)) call cross_mn(p%oo,p%uu,p%oxu)
! oxuxb
      if (lpencil(i_oxuxb)) call cross_mn(p%oxu,p%bb,p%oxuxb)
! jxbxb
      if (lpencil(i_jxbxb)) call cross_mn(p%jxb,p%bb,p%jxbxb)
! jxbrxb
      if (lpencil(i_jxbrxb)) call cross_mn(p%jxbr,p%bb,p%jxbrxb)
! glnrhoxb
      if (lpencil(i_glnrhoxb)) call cross_mn(p%glnrho,p%bb,p%glnrhoxb)
! del4a
      if (lpencil(i_del4a)) call del4v(f,iaa,p%del4a)
! del6a
      if (lpencil(i_del6a)) call del6v(f,iaa,p%del6a)
! oxj
      if (lpencil(i_oxj)) call cross_mn(p%oo,p%jj,p%oxJ)
! jij
      if (lpencil(i_jij)) then
        do j=1,3
          do i=1,3
            p%jij(:,i,j)=.5*(p%bij(:,i,j)+p%bij(:,j,i))
          enddo
        enddo
      endif
! sj
      if (lpencil(i_sj)) call multmm_sc(p%sij,p%jij,p%sj)
! ss12
      if (lpencil(i_ss12)) p%ss12=sqrt(abs(p%sj))
!
! mf_EMF
! needed if a mean field (mf) model is calculated
!
      if (lpencil(i_mf_EMF)) then
        kx=2*pi/Lx
!
!  possibility of dynamical alpha
!
        if (lalpm.and..not.lmeanfield_noalpm) then
          alpha_total=alpha_effect*alpha_tmp+f(l1:l2,m,n,ialpm)
        else
          alpha_total=alpha_effect*alpha_tmp
        endif
!
!  possibility of conventional alpha quenching (rescales alpha_total)
!  initialize EMF with alpha_total*bb
!
        if (alpha_quenching/=0.) alpha_total=alpha_total/(1.+alpha_quenching*p%b2)
        call multsv_mn(alpha_total,p%bb,p%mf_EMF)
!
!  add possible delta x J effect and turbulent diffusion to EMF
!
        if (delta_effect/=0.) p%mf_EMF=p%mf_EMF+delta_effect*p%oxJ
        if (meanfield_etat/=0.) then
          if (lweyl_gauge) then
            if (meanfield_etaB/=0.) then
              meanfield_etaB2=meanfield_etaB**2
              call multsv_mn_add(meanfield_etat/sqrt(1.+p%b2/meanfield_etaB2),p%jj,p%mf_EMF)
            else
              p%mf_EMF=p%mf_EMF-meanfield_etat*p%jj
            endif
          else
            p%mf_EMF=p%mf_EMF+meanfield_etat*p%del2a
          endif
!
!  allow for possibility of variable etat
!
          if (ietat/=0) then
            call multsv_mn_add(-f(l1:l2,m,n,ietat),p%jj,p%mf_EMF)
          endif
        endif
      endif
      if (lpencil(i_mf_EMFdotB)) call dot_mn(p%mf_EMF,p%bb,p%mf_EMFdotB)
!
!  Store bb in auxiliary variable if requested.
!  Just neccessary immediately before writing snapshots, but how would we
!  know we are?
!
     if (lbb_as_aux) f(l1:l2,m,n,ibx:ibz)=p%bb
     if (ljj_as_aux) f(l1:l2,m,n,ijx:ijz)=p%jj
!
    endsubroutine calc_pencils_magnetic
!***********************************************************************
    subroutine daa_dt(f,df,p)
!
!  magnetic field evolution
!
!  calculate dA/dt=uxB+3/2 Omega_0 A_y x_dir -eta mu_0 J
!  for mean field calculations one can also add dA/dt=...+alpha*bb+delta*WXJ
!  add jxb/rho to momentum equation
!  add eta mu_0 j2/rho to entropy equation
!
!   7-oct-09/axel: adapted from magnetic
!
      use Deriv, only: der6
      use Diagnostics
      use EquationOfState, only: eoscalc,gamma_m1
      use Io, only: output_pencil
      use Mpicomm, only: stop_it
      use Special, only: special_calc_magnetic
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: bpol
      real, dimension (nx) :: aphi,bphi,d2aphi,d2bphi,bpol2
      real, dimension (nx) :: alpha_tmp,alpha_tmp2,Omega_tmp
      real :: eta_tot
      integer :: i
      integer, parameter :: nxy=nxgrid*nygrid
!
      intent(inout)  :: f,p
      intent(inout)  :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'daa_dt: SOLVE'
      if (headtt) then
        call identify_bcs('Aphi',iaphi)
        call identify_bcs('Bphi',ibphi)
      endif
!
!  define aphi and bphi, and calculate bpol
!
      aphi=f(l1:l2,m,n,iaphi)
      bphi=f(l1:l2,m,n,ibphi)
      call curl_horizontal(f,iaphi,Bpol)
!
!  Diffusion operator
!
      call del2(f,iaphi,d2aphi)
      call del2(f,ibphi,d2bphi)
!
!  calculate D2a = del2a - 1/pomega^2 and same for b
!
      d2aphi=d2aphi-aphi*r2_mn*sin2th(m)
      d2bphi=d2bphi-bphi*r2_mn*sin2th(m)
!
!  add diffusion term to the right-hand side.
!  Use total eta, eta_tot, which includes microscopic and meanfield eta.
!
      eta_tot=eta+meanfield_etat
      df(l1:l2,m,n,iaphi)=df(l1:l2,m,n,iaphi)+eta_tot*d2aphi
      df(l1:l2,m,n,ibphi)=df(l1:l2,m,n,ibphi)+eta_tot*d2bphi
!
!  add alpha effect, note that j=-D2a
!
      alpha_tmp=alpha_effect*alpha_x*alpha_y(m)
     alpha_tmp2=alpha_effect*(dalpha_x*Bpol(:,2)-r1_mn*dalpha_y(m)*Bpol(:,1))
!
!  Add alpha effect. At the moment this ignores the grad(alpha) terms
!
      df(l1:l2,m,n,iaphi)=df(l1:l2,m,n,iaphi)+alpha_tmp*bphi
      if (lalpha_phi_equation) then
        df(l1:l2,m,n,ibphi)=df(l1:l2,m,n,ibphi)-alpha_tmp*d2aphi+alpha_tmp2
      endif
!
!  differential rotation, need poloidal field for this, and
!  add pomega*Bpol.grad(Omega) to the dBphi/dt equation.
!  Note that grad_theta(Omega)=r1_mn*dOmega_y, but r1_mn cancels with r_mn.
!
      if (lOmega_effect) then
 !!      Omega_tmp=sinth(m)*(r_mn*Bpol(:,1)*dOmega_x+Bpol(:,2)*dOmega_y(m))
       Omega_tmp=sinth(m)*(r_mn*Bpol(:,1)*dOmega_x*Omega_y(m)+2*Bpol(:,2)*dOmega_y(m)*Omega_x)
        df(l1:l2,m,n,ibphi)=df(l1:l2,m,n,ibphi)+Omega_ampl*Omega_tmp
      endif
!
!  allow for special routines
!
      if (lspecial) call special_calc_magnetic(f,df,p)
!
!  Multiply resistivity by Nyquist scale, for resistive time-step.
!  We include possible contribution from meanfield_etat, which is however
!  only invoked in mean field models.
!  Allow for variable etat (mean field theory)
!
      if (lfirst.and.ldt) then
        diffus_eta=diffus_eta+eta_tot
        diffus_eta=diffus_eta*dxyz_2
      endif
!
      if (headtt.or.ldebug) then
        print*, 'daa_dt: max(diffus_eta)  =', maxval(diffus_eta)
      endif
!
!     if (linduction) &
!        df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + uxb_upw + fres
!     endif
!
!  Alpha effect
!  additional terms if Mean Field Theory is included
!
!     if (lmeanfield_theory.and. &
!       (meanfield_etat/=0..or.alpha_effect/=0..or.delta_effect/=0.)) then
!       df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+p%mf_EMF
!       if (lOmega_effect) call Omega_effect(f,df)
!     endif
!
      if (ldiagnos) then
        if (idiag_aphi2m/=0) call sum_mn_name(aphi**2,idiag_aphi2m)
        if (idiag_bphi2m/=0) call sum_mn_name(bphi**2,idiag_bphi2m)
        if (idiag_bpol2m/=0) then
          call dot2(Bpol,Bpol2)
          call sum_mn_name(Bpol2,idiag_bpol2m)
        endif
      endif
!
    endsubroutine daa_dt
!***********************************************************************
    subroutine time_integrals_magnetic(f,p)
!
!  Calculate time_integrals within each pencil (as long as each
!  pencil case p still contains the current data). This routine
!  is now being called at the end of equ.
!
!  28-jun-07/axel+mreinhard: coded
!  24-jun-08/axel: moved call to this routine to the individual pde routines
!   1-jul-08/axel: moved this part to magnetic
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(inout) :: f
      intent(in) :: p
!
      if (ibbt/=0) f(l1:l2,m,n,ibxt:ibzt)=f(l1:l2,m,n,ibxt:ibzt)+dt*p%bb
      if (ijjt/=0) f(l1:l2,m,n,ijxt:ijzt)=f(l1:l2,m,n,ijxt:ijzt)+dt*p%jj
!
    endsubroutine time_integrals_magnetic
!***********************************************************************
    subroutine df_diagnos_magnetic(df,p)
!
!  calculate diagnostics that involves df
!  Here we calculate <du/dt x b> and <u x db/dt>.
!  The latter is calculated as <divu dai/dt> -  <uji daj/dt>
!  This is used in dynamo theory for checking the minimal tau approximation.
!
!  10-oct-06/axel: coded
!
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: uudot,aadot,udotxb,B1_gradu
      real, dimension (nx) :: B1dot_udotxb,B1dot_uxbdot,B1dot_aadot,uxbdot2
!
      intent(in)  :: df, p
!
    endsubroutine df_diagnos_magnetic
!***********************************************************************
    subroutine set_border_magnetic(f,df,p)
!
!  Calculates the driving term for the border profile
!  of the aa variable.
!
!  28-jul-06/wlad: coded
!
      use BorderProfiles, only: border_driving,set_border_initcond
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx,3) :: f_target
      integer :: ju,j
!
      select case (borderaa)
!
      case ('zero','0')
        f_target=0.
!
      case ('initial-condition')
        do j=1,3
          ju=j+iaa-1
          call set_border_initcond(f,ju,f_target(:,j))
        enddo
!
      case ('nothing')
        if (lroot.and.ip<=5) &
             print*,"set_border_magnetic: borderaa='nothing'"
!
      case default
        write(unit=errormsg,fmt=*) &
             'set_border_magnetic: No such value for borderaa: ', &
             trim(borderaa)
        call fatal_error('set_border_magnetic',errormsg)
      endselect
!
      if (borderaa /= 'nothing') then
        do j=1,3
          ju=j+iaa-1
          call border_driving(f,df,p,f_target(:,j),ju)
        enddo
      endif
!
    endsubroutine set_border_magnetic
!***********************************************************************
    subroutine eta_shell(p,eta_mn,geta)
!
!  24-nov-03/dave: coded
!  23-jun-09/axel: generalized to lcylinder_in_a_box
!
      use Sub, only: step, der_step
      use Mpicomm, only: stop_it
!
      type (pencil_case) :: p
      real, dimension (nx) :: eta_mn
      real, dimension (nx) :: prof,eta_r
      real, dimension (nx,3) :: geta
      real :: d_int,d_ext
!
      eta_r=0.
!
      if (eta_int > 0.) then
        d_int = eta_int - eta
      else
        d_int = 0.
      endif
      if (eta_ext > 0.) then
        d_ext = eta_ext - eta
      else
        d_ext = 0.
      endif
!
!  calculate steps in resistivity
!  make this dependent on the geometry used
!
!  (i) lcylinder_in_a_box
!
      if (lcylinder_in_a_box.or.lcylindrical_coords) then
        prof=step(p%rcyl_mn,r_int,wresistivity)
        eta_mn=d_int*(1-prof)
        prof=step(p%rcyl_mn,r_ext,wresistivity)
        eta_mn=eta+eta_mn+d_ext*prof
!
!     calculate radial derivative of steps and gradient of eta
!
        prof=der_step(p%rcyl_mn,r_int,wresistivity)
        eta_r=-d_int*prof
        prof=der_step(p%rcyl_mn,r_ext,wresistivity)
        eta_r=eta_r+d_ext*prof
        geta=p%evr*spread(eta_r,2,3)
!
!  (ii) lsphere_in_a_box
!
      elseif (lsphere_in_a_box.or.lspherical_coords) then
        prof=step(p%r_mn,r_int,wresistivity)
        eta_mn=d_int*(1-prof)
        prof=step(p%r_mn,r_ext,wresistivity)
        eta_mn=eta+eta_mn+d_ext*prof
!
!     calculate radial derivative of steps and gradient of eta
!
        prof=der_step(p%r_mn,r_int,wresistivity)
        eta_r=-d_int*prof
        prof=der_step(p%r_mn,r_ext,wresistivity)
        eta_r=eta_r+d_ext*prof
        geta=p%evr*spread(eta_r,2,3)
!
!  (iii) other cases are not implemented yet
!
      else
        call stop_it("eta_shell works only for spheres or cylinders")
      endif
!
    endsubroutine eta_shell
!***********************************************************************
    subroutine calc_bthresh()
!
!  calculate bthresh from brms, give warnings if there are problems
!
!   6-aug-03/axel: coded
!
!  give warning if brms is not set in prints.in
!
      if (idiag_brms==0) then
        if (lroot.and.lfirstpoint) then
          print*,'calc_bthresh: need to set brms in print.in to get bthresh'
        endif
      endif
!
!  if nvec exceeds nbvecmax (=1/4) of points per processor, then begin to
!  increase scaling factor on bthresh. These settings will stay in place
!  until the next restart
!
      if (nbvec>nbvecmax.and.lfirstpoint) then
        print*,'calc_bthresh: processor ',iproc,': bthresh_scl,nbvec,nbvecmax=', &
                                                   bthresh_scl,nbvec,nbvecmax
        bthresh_scl=bthresh_scl*1.2
      endif
!
!  calculate bthresh as a certain fraction of brms
!
      bthresh=bthresh_scl*bthresh_per_brms*brms
!
    endsubroutine calc_bthresh
!***********************************************************************
    subroutine rescaling_magnetic(f)
!
!  Rescale magnetic field by factor rescale_aa,
!
!  22-feb-05/axel: coded
!  10-feb-09/petri: adapted from testfield
!
      use Sub, only: update_snaptime, read_snaptime
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=fnlen) :: file
      logical :: lmagnetic_out
      integer,save :: ifirst=0
!
      intent(inout) :: f
!
!  Reinitialize aa periodically if requested
!
      if (lreset_aa) then
        file=trim(datadir)//'/treset_aa.dat'
        if (ifirst==0) then
          call read_snaptime(trim(file),taareset,naareset,daareset,t)
          if (taareset==0 .or. taareset < t-daareset) then
            taareset=t+daareset
          endif
          ifirst=1
        endif
!
!  Rescale when the time has come
!  (Note that lmagnetic_out and ch are not used here)
!
        if (t >= taareset) then
          f(:,:,:,iax:iaz)=rescale_aa*f(:,:,:,iax:iaz)
          call update_snaptime(file,taareset,naareset,daareset,t,lmagnetic_out)
        endif
      endif
!
    endsubroutine rescaling_magnetic
!***********************************************************************
    subroutine calc_tau_aa_exterior(f,df)
!
!  magnetic field relaxation to zero on time scale tau_aa_exterior within
!  exterior region. For the time being this means z > zgrav.
!
!  29-jul-02/axel: coded
!
      use Gravity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: scl
      integer :: j
!
      intent(in) :: f
      intent(out) :: df
!
      if (headtt) print*,'calc_tau_aa_exterior: tau=',tau_aa_exterior
      if (z(n)>zgrav) then
        scl=1./tau_aa_exterior
        do j=iax,iaz
          df(l1:l2,m,n,j)=df(l1:l2,m,n,j)-scl*f(l1:l2,m,n,j)
        enddo
      endif
!
    endsubroutine calc_tau_aa_exterior
!***********************************************************************
    subroutine read_magnetic_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=magnetic_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=magnetic_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_magnetic_init_pars
!***********************************************************************
    subroutine write_magnetic_init_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=magnetic_init_pars)
!
    endsubroutine write_magnetic_init_pars
!***********************************************************************
    subroutine read_magnetic_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=magnetic_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=magnetic_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_magnetic_run_pars
!***********************************************************************
    subroutine write_magnetic_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=magnetic_run_pars)
!
    endsubroutine write_magnetic_run_pars
!***********************************************************************
    subroutine get_slices_magnetic(f,slices)
!
!  Write slices for animation of Magnetic variables.
!
!  26-jul-06/tony: coded
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
        case ('aa')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =f(ix_loc,m1:m2 ,n1:n2  ,iax-1+slices%index)
            slices%xz =f(l1:l2 ,iy_loc,n1:n2  ,iax-1+slices%index)
            slices%xy =f(l1:l2 ,m1:m2 ,iz_loc ,iax-1+slices%index)
            slices%xy2=f(l1:l2 ,m1:m2 ,iz2_loc,iax-1+slices%index)
            if (lwrite_slice_xy3) &
                 slices%xy3=f(l1:l2,m1:m2,iz3_loc,iax-1+slices%index)
            if (lwrite_slice_xy4) &
                 slices%xy4=f(l1:l2,m1:m2,iz4_loc,iax-1+slices%index)
            if (slices%index<=3) slices%ready=.true.
          endif
!
!  Magnetic field (derived variable)
!
        case ('bb')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =>bb_yz(:,:,slices%index)
            slices%xz =>bb_xz(:,:,slices%index)
            slices%xy =>bb_xy(:,:,slices%index)
            slices%xy2=>bb_xy2(:,:,slices%index)
            if (lwrite_slice_xy3) slices%xy3=>bb_xy3(:,:,slices%index)
            if (lwrite_slice_xy4) slices%xy4=>bb_xy4(:,:,slices%index)
            if (slices%index<=3) slices%ready=.true.
          endif
!
!  Magnetic field (derived variable)
!
        case ('jj')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =>jj_yz(:,:,slices%index)
            slices%xz =>jj_xz(:,:,slices%index)
            slices%xy =>jj_xy(:,:,slices%index)
            slices%xy2=>jj_xy2(:,:,slices%index)
            if (lwrite_slice_xy3) slices%xy3=>jj_xy3(:,:,slices%index)
            if (lwrite_slice_xy4) slices%xy4=>jj_xy4(:,:,slices%index)
            if (slices%index<=3) slices%ready=.true.
          endif
!
!  Magnetic field squared (derived variable)
!
        case ('b2')
          slices%yz =>b2_yz
          slices%xz =>b2_xz
          slices%xy =>b2_xy
          slices%xy2=>b2_xy2
          if (lwrite_slice_xy3) slices%xy3=>b2_xy3
          if (lwrite_slice_xy4) slices%xy4=>b2_xy4
          slices%ready=.true.
!
!  Current density (derived variable)
!
        case ('jb')
          slices%yz =>jb_yz
          slices%xz =>jb_xz
          slices%xy =>jb_xy
          slices%xy2=>jb_xy2
          if (lwrite_slice_xy3) slices%xy3=>jb_xy3
          if (lwrite_slice_xy4) slices%xy4=>jb_xy4
          slices%ready=.true.
!
!  Plasma beta
!
       case ('beta1')
          slices%yz =>beta1_yz
          slices%xz =>beta1_xz
          slices%xy =>beta1_xy
          slices%xy2=>beta1_xy2
          if (lwrite_slice_xy3) slices%xy3=>beta1_xy3
          if (lwrite_slice_xy4) slices%xy4=>beta1_xy4
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_magnetic
!***********************************************************************
    subroutine calc_mfield
!
!  calculate mean magnetic field from xy- or z-averages
!
!  19-jun-02/axel: moved from print to here
!   9-nov-02/axel: corrected bxmy(m,j); it used bzmy instead!
!
      use Mpicomm
      use Sub
!
!  For vector output (of bb vectors) we need brms
!  on all processors. It suffices to have this for times when lout=.true.,
!  but we need to broadcast the result to all procs.
!
!  calculate brms (this requires that brms is set in print.in)
!  broadcast result to other processors
!
!  The following calculation involving spatial averages
!
    endsubroutine calc_mfield
!***********************************************************************
    subroutine alfven_x(ampl,f,iuu,iaa,ilnrho,kx)
!
!  Alfven wave propagating in the x-direction
!
!  ux = +sink(x-vA*t)
!  Az = -cosk(x-vA*t)*sqrt(rho*mu0)/k
!
!  Alfven and slow magnetosonic are the same here and both incompressible, and
!  a fast magnetosonic (compressible) wave is also excited, but decoupled.
!
!  satisfies the four equations
!  dlnrho/dt = -ux'
!  dux/dt = -cs2*(lnrho)'
!  duy/dt = B0*By'  ==>  duy/dt = -B0*Az''
!  dBy/dt = B0*uy'  ==>  dAz/dt = -B0*ux
!
!   8-nov-03/axel: coded
!  29-apr-03/axel: added sqrt(rho*mu0)/k factor
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rho,ampl_Az
      real :: ampl,kx,ampl_lr,ampl_ux,ampl_uy
      integer :: iuu,iaa,ilnrho
!
!  Amplitude factors
!
      ampl_lr=+0.
      ampl_ux=+0.
      ampl_uy=+ampl
!
!  ux and Ay.
!  Don't overwrite the density, just add to the log of it.
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,ilnrho)=ampl_lr*(sin(kx*x(l1:l2))+f(l1:l2,m,n,ilnrho))
        f(l1:l2,m,n,iuu+0 )=ampl_ux*sin(kx*x(l1:l2))
        f(l1:l2,m,n,iuu+1 )=ampl_uy*sin(kx*x(l1:l2))
        rho=exp(f(l1:l2,m,n,ilnrho))
        ampl_Az=-ampl*sqrt(rho*mu0)/kx
        f(l1:l2,m,n,iaa+2 )=ampl_Az*cos(kx*x(l1:l2))
      enddo; enddo
!
    endsubroutine alfven_x
!***********************************************************************
    subroutine alfven_y(ampl,f,iuu,iaa,ky,mu0)
!
!  Alfven wave propagating in the y-direction; can be used in 2-d runs.
!  ux = cos(ky-ot), for B0y=1 and rho=1.
!  Az = sin(ky-ot), ie Bx=-cos(ky-ot)
!
!  [wd nov-2006: There should be a 1/ky in the aa term here and in
!  alfven_x, I think]
!
!  satisfies the equations
!  dux/dt = Bx'  ==>  dux/dt = -Az''
!  dBx/dt = ux'  ==>  dAz/dt = -ux.
!
!  06-dec-06/wolf: adapted from alfven_z
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,ky,mu0
      integer :: iuu,iaa
!
!  ux and Az
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+0) = +ampl*cos(ky*y(m))
        f(l1:l2,m,n,iaa+2) = -ampl*sin(ky*y(m))*sqrt(mu0)/ky
      enddo; enddo
!
    endsubroutine alfven_y
!***********************************************************************
    subroutine alfven_z(ampl,f,iuu,iaa,kz,mu0)
!
!  Alfven wave propagating in the z-direction
!  ux = cos(kz-ot), for B0z=1 and rho=1.
!  Ay = sin(kz-ot), ie Bx=-cos(kz-ot)
!
!  satisfies the equations
!  dux/dt = Bx'  ==>  dux/dt = -Ay''
!  dBx/dt = ux'  ==>  dAy/dt = -ux.
!
!  18-aug-02/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kz,mu0
      integer :: iuu,iaa
!
!  ux and Ay
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+0)=+ampl*cos(kz*z(n))
        f(l1:l2,m,n,iaa+1)=+ampl*sin(kz*z(n))*sqrt(mu0)
      enddo; enddo
!
    endsubroutine alfven_z
!***********************************************************************
    subroutine alfven_xy(ampl,f,iuu,iaa,kx,ky)
!
!  Alfven wave propagating in the xy-direction; can be used in 2-d runs.
!  uz = cos(kx*x+ky*y-ot), for B0=(1,1,0) and rho=1.
!  Ax = sin(kx*x+ky*y-ot),
!  Ay = sin(kx*x+ky*y-ot),
!
!  16-jun-07/axel: adapted from alfven_y
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx,ky,om
      real, parameter :: mu0=1.
      integer :: iuu,iaa
!
!  set ux, Ax, and Ay
!
      om=B_ext(1)*kx+B_ext(2)*ky
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+2)=+ampl*cos(kx*x(l1:l2)+ky*y(m))
        f(l1:l2,m,n,iaa+0)=+ampl*sin(kx*x(l1:l2)+ky*y(m))*sqrt(mu0)/om*B_ext(2)
        f(l1:l2,m,n,iaa+1)=-ampl*sin(kx*x(l1:l2)+ky*y(m))*sqrt(mu0)/om*B_ext(1)
      enddo; enddo
!
    endsubroutine alfven_xy
!***********************************************************************
    subroutine alfven_xz(ampl,f,iuu,iaa,kx,kz)
!
!  Alfven wave propagating in the xz-direction; can be used in 2-d runs.
!  uz = cos(kx*x+kz*z-ot), for B0=(1,1,0) and rho=1.
!  Ax = sin(kx*x+kz*z-ot),
!  Az = sin(kx*x+kz*z-ot),
!
!  16-jun-07/axel: adapted from alfven_xy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kx,kz,om
      real, parameter :: mu0=1.
      integer :: iuu,iaa
!
!  set ux, Ax, and Az
!
      om=B_ext(1)*kx+B_ext(3)*kz
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+2)=+ampl*cos(kx*x(l1:l2)+kz*z(n))
        f(l1:l2,m,n,iaa+0)=+ampl*sin(kx*x(l1:l2)+kz*z(n))*sqrt(mu0)/om*B_ext(2)
        f(l1:l2,m,n,iaa+2)=-ampl*sin(kx*x(l1:l2)+kz*z(n))*sqrt(mu0)/om*B_ext(1)
      enddo; enddo
!
    endsubroutine alfven_xz
!***********************************************************************
    subroutine alfvenz_rot(ampl,f,iuu,iaa,kz,O)
!
!  Alfven wave propagating in the z-direction (with Coriolis force)
!  ux = cos(kz-ot), for B0z=1 and rho=1.
!  Ay = sin(kz-ot), ie Bx=-cos(kz-ot)
!
!  satisfies the equations
!  dux/dt - 2Omega*uy = -Ay''
!  duy/dt + 2Omega*ux = +Ax''
!  dAx/dt = +uy
!  dAy/dt = -ux
!
!  18-aug-02/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kz,O,fac
      integer :: iuu,iaa
!
!  ux, uy, Ax and Ay
!
      if (lroot) print*,'alfvenz_rot: Alfven wave with rotation; O,kz=',O,kz
      fac=-O+sqrt(O**2+kz**2)
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+0)=-ampl*sin(kz*z(n))*fac/kz
        f(l1:l2,m,n,iuu+1)=-ampl*cos(kz*z(n))*fac/kz
        f(l1:l2,m,n,iaa+0)=+ampl*sin(kz*z(n))/kz
        f(l1:l2,m,n,iaa+1)=+ampl*cos(kz*z(n))/kz
      enddo; enddo
!
    endsubroutine alfvenz_rot
!***********************************************************************
    subroutine alfvenz_rot_shear(ampl,f,iuu,iaa,kz,OO)
!
!  Alfven wave propagating in the z-direction (with Coriolis force and shear)
!
!  satisfies the equations
!  dux/dt - 2*Omega*uy = -Ay''
!  duy/dt + (2-q)*Omega*ux = +Ax''
!  dAx/dt = q*Omega*Ay + uy
!  dAy/dt = -ux
!
!  Assume B0=rho0=mu0=1
!
!  28-june-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl,kz,OO
      complex :: fac
      integer :: iuu,iaa
!
!  ux, uy, Ax and Ay
!
      if (lroot) print*,'alfvenz_rot_shear: '// &
          'Alfven wave with rotation and shear; OO,kz=',OO,kz
      fac=cmplx(OO-sqrt(16*kz**2+OO**2),0.)
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+0)=f(l1:l2,m,n,iuu+0)+ampl*fac/(4*kz)*sin(kz*z(n))
        f(l1:l2,m,n,iuu+1)=f(l1:l2,m,n,iuu+1)+ampl*real(exp(cmplx(0,z(n)*kz))* &
            fac*sqrt(2*kz**2+OO*fac)/(sqrt(2.)*kz*(-6*OO-fac)))
        f(l1:l2,m,n,iaa+0)=ampl*sin(kz*z(n))/kz
        f(l1:l2,m,n,iaa+1)=-ampl*2*sqrt(2.)*aimag(exp(cmplx(0,z(n)*kz))* &
            sqrt(2*kz**2+OO*fac)/(-6*OO-fac)/(cmplx(0,kz)))
      enddo; enddo
!
    endsubroutine alfvenz_rot_shear
!***********************************************************************
    subroutine fluxrings(ampl,f,ivar1,ivar2,profile)
!
!  Magnetic flux rings. Constructed from a canonical ring which is the
!  rotated and translated:
!    AA(xxx) = D*AA0(D^(-1)*(xxx-xxx_disp)) ,
!  where AA0(xxx) is the canonical ring and D the rotation matrix
!  corresponding to a rotation by phi around z, followed by a
!  rotation by theta around y.
!  The array was already initialized to zero before calling this
!  routine.
!  Optional argument `profile' allows to choose a different profile (see
!  norm_ring())
!
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: tmpv
      real, dimension (nx) :: xx1,yy1,zz1
      real, dimension(3) :: axis,disp
      real :: ampl,phi,theta,ct,st,cp,sp
      real :: fring,Iring,R0,width
      integer :: i,ivar,ivar1,ivar2,ivar3
      character (len=*), optional :: profile
      character (len=labellen) :: prof
!
      if (present(profile)) then
        prof = profile
      else
        prof = 'tanh'
      endif
!
!  fix ivar3=ivar1 (for now)
!
      ivar3=ivar1
!
!  initialize each ring
!
      if (any((/fring1,fring2,Iring1,Iring2/) /= 0.)) then
        ! fringX is the magnetic flux, IringX the current
        if (lroot) then
          print*, 'fluxrings: Initialising magnetic flux rings'
        endif
        do i=1,nrings
          if (i==1) then
            fring = fring1      ! magnetic flux along ring
            Iring = Iring1      ! current along ring (for twisted flux tube)
            R0    = Rring1      ! radius of ring
            width = wr1         ! ring thickness
            axis  = axisr1      ! orientation
            disp  = dispr1      ! position
            ivar  = ivar1
          elseif (i==2) then
            fring = fring2
            Iring = Iring2
            R0    = Rring2
            width = wr2
            axis  = axisr2
            disp  = dispr2
            ivar  = ivar2
          elseif (i==3) then
            fring = fring3
            Iring = Iring3
            R0    = Rring3
            width = wr3
            axis  = axisr3
            disp  = dispr3
            ivar  = ivar3
          else
            call stop_it('fluxrings: nrings is too big')
          endif
          phi   = atan2(axis(2),axis(1)+epsi)
          theta = atan2(sqrt(axis(1)**2+axis(2)**2)+epsi,axis(3))
          ct = cos(theta); st = sin(theta)
          cp = cos(phi)  ; sp = sin(phi)
          ! Calculate D^(-1)*(xxx-disp)
          do n=n1,n2; do m=m1,m2
            xx1= ct*cp*(x(l1:l2)-disp(1))+ct*sp*(y(m)-disp(2))-st*(z(n)-disp(3))
            yy1=-   sp*(x(l1:l2)-disp(1))+   cp*(y(m)-disp(2))
            zz1= st*cp*(x(l1:l2)-disp(1))+st*sp*(y(m)-disp(2))+ct*(z(n)-disp(3))
            call norm_ring(xx1,yy1,zz1,fring,Iring,R0,width,tmpv,PROFILE=prof)
            ! calculate D*tmpv
            f(l1:l2,m,n,ivar  ) = f(l1:l2,m,n,ivar  ) + ampl*( &
                 + ct*cp*tmpv(:,1) - sp*tmpv(:,2) + st*cp*tmpv(:,3))
            f(l1:l2,m,n,ivar+1) = f(l1:l2,m,n,ivar+1) + ampl*( &
                 + ct*sp*tmpv(:,1) + cp*tmpv(:,2) + st*sp*tmpv(:,3))
            f(l1:l2,m,n,ivar+2) = f(l1:l2,m,n,ivar+2) + ampl*( &
                 - st   *tmpv(:,1)                + ct   *tmpv(:,3))
          enddo; enddo
        enddo
      endif
      if (lroot) print*, 'fluxrings: Magnetic flux rings initialized'
!
    endsubroutine fluxrings
!***********************************************************************
    subroutine norm_ring(xx1,yy1,zz1,fring,Iring,R0,width,vv,profile)
!
!  Generate vector potential for a flux ring of magnetic flux FRING,
!  current Iring (not correctly normalized), radius R0 and thickness
!  WIDTH in normal orientation (lying in the x-y plane, centred at (0,0,0)).
!
!   1-may-02/wolf: coded (see the manual, Section C.3)
!   7-jun-09/axel: added gaussian and constant (or box) profiles
!
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (nx,3) :: vv
      real, dimension (nx) :: xx1,yy1,zz1,phi,tmp
      real :: fring,Iring,R0,width
      character (len=*) :: profile
!
!  magnetic ring, define r-R
!
      tmp = sqrt(xx1**2+yy1**2)-R0
!
!  choice of different profile functions
!
      select case (profile)
!
!  gaussian profile, exp(-.5*(x/w)^2)/(sqrt(2*pi)*eps),
!  so its derivative is .5*(1.+erf(-x/(sqrt(2)*eps))
!
      case ('gaussian')
        vv(:,3) = - fring * .5*(1.+erfunc(tmp/(sqrt(2.)*width))) &
                          * exp(-.5*(zz1/width)**2)/(sqrt(2.*pi)*width)
!
!  tanh profile, so the delta function is approximated by 1/cosh^2.
!  The name tanh is misleading, because the actual B frofile is
!  1./cosh^2, but this is harder to write.
!
      case ('tanh')
        vv(:,3) = - fring * 0.5*(1+tanh(tmp/width)) &
                          * 0.5/width/cosh(zz1/width)**2
!
!  constant profile, so the delta function is approximated by the function
!  delta(x) = 1/2w, if -w < x < w.
!
      case ('const')
        vv(:,3) = - fring * 0.5*(1.+max(-1.,min(tmp/width,1.))) &
                          * 0.25/width*(1.-sign(1.,abs(zz1)-width))
!
!  there is no default option here
!
      case default
        call stop_it('norm_ring: No such fluxtube profile')
      endselect
!
!  current ring (to twist the B-lines)
!
      tmp = width - sqrt(tmp**2 + zz1**2)
      tmp = Iring*0.5*(1+tanh(tmp/width))     ! Now the A_phi component
      phi = atan2(yy1,xx1)
      vv(:,1) = - tmp*sin(phi)
      vv(:,2) =   tmp*cos(phi)
!
    endsubroutine norm_ring
!***********************************************************************
    subroutine torus_test(ampl,f)
!
!  Initial field concentrated along torus well inside the computational
!  domain.
!  Implements the same field for cartesian and spherical cordinates.
!  The field is of mixed parity (bb_pol symmetric, bb_tor antisymmetric)
!  and the relative contributions of toroidal and poloidal field are
!  determined by
!    ampl(1) -- bb_pol (through aa_tor)
!    ampl(3) -- bb_tor (through aa_pol)
!  Uses x_max as reference radius.
!
!   05-may-2008/wolf: coded
!
      real :: ampl
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx) :: xxi2,ee
      real, dimension (nx) :: costh,sinth,cosphi,sinphi,ss,rr,aar,aap
      real :: radius,width,r_cent
!
      intent(in)    :: ampl
      intent(inout) :: f
!
      radius = xyz1(1)
      width  = 0.1*radius
      r_cent = 0.6*radius
!
      if (lspherical_coords) then
        do n=n1,n2; do m=m1,m2
          xxi2 = (x(l1:l2)*sin(y(m))-r_cent)**2 + x(l1:l2)**2*cos(y(m))**2
          ee = ampl * exp(-0.5 * xxi2 / width**2)
          f(l1:l2,m,n,iax) = f(l1:l2,m,n,iax) + ee * x(l1:l2)*cos(y(m))
          f(l1:l2,m,n,iaz) = f(l1:l2,m,n,iaz) + ee
        enddo; enddo
      else
        do n=n1,n2; do m=m1,m2
          xxi2 = (sqrt(x(l1:l2)**2+y(m)**2) - r_cent)**2 + z(n)**2
          ee = ampl * exp(-0.5 * xxi2 / width**2)
          aar = z(n) * ee
          aap = ee
          ss = sqrt(x(l1:l2)**2+y(m)**2)
          rr = sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
          ss = max(ss, tini)
          rr = max(rr, tini)
          costh = z(n)/rr
          sinth = ss/rr
          cosphi   = x(l1:l2)/ss
          sinphi   = y(m)/ss
          f(l1:l2,m,n,iax) = f(l1:l2,m,n,iax) + aar*sinth*cosphi - aap*sinphi
          f(l1:l2,m,n,iay) = f(l1:l2,m,n,iay) + aar*sinth*sinphi + aap*cosphi
          f(l1:l2,m,n,iaz) = f(l1:l2,m,n,iaz) + aar*costh
        enddo; enddo
      endif
!
    endsubroutine torus_test
!***********************************************************************
    subroutine force_free_jet(mu)
!
!  Force free magnetic field configuration for jet simulations
!  with a fixed accretion disc at the bottom boundary.
!
!  The input parameter mu specifies the radial dependency of
!  the magnetic field in the disc.
!
!  Solves the laplace equation in cylindrical coordinates for the
!  phi-component of the vector potential. A_r and A_z are taken to
!  be zero.
!
!    nabla**2 A_phi - A_phi / r**2 = 0
!
!  For the desired boundary condition in the accretion disc
!
!    B_r=B0*r**(mu-1)  (z == 0)
!
!  the solution is
!
!    A_phi = Hypergeometric2F1( (1-mu)/2, (2+mu)/2, 2, xi**2 )
!            *xi*(r**2+z**2)**(mu/2)
!
!  where xi = sqrt(r**2/(r**2+z**2))
!
!
!  30-may-04/tobi: coded
!
      use Sub, only: hypergeometric2F1,gamma_function
      use Deriv, only: der
      use IO, only: output
!
      real, intent(in) :: mu
      real :: xi2,A_phi
      real :: r2
      real :: B1r_,B1z_,B1
      real, parameter :: tol=10*epsilon(1.0)
      integer :: l
      real, dimension(mx,my,mz) :: Ax_ext,Ay_ext
      real, dimension(nx,3) :: bb_ext_pot
      real, dimension(nx) :: bb_x,bb_y,bb_z
!
!  calculate un-normalized |B| at r=r_ref and z=0 for later normalization
!
      if (lroot.and.ip<=5) print*,'FORCE_FREE_JET: calculating normalization'
!
      B1r_=sin(pi*mu/2)*gamma_function(   abs(mu) /2) / &
                        gamma_function((1+abs(mu))/2)
!
      B1z_=cos(pi*mu/2)*gamma_function((1+abs(mu))/2) / &
                        gamma_function((2+abs(mu))/2)
!
      B1=sqrt(4/pi)*r_ref**(mu-1)*sqrt(B1r_**2+B1z_**2)
!
!  calculate external vector potential
!
      if (lroot) print*,'FORCE_FREE_JET: calculating external vector potential'
!
      if (lforce_free_test) then
!
        if (lroot) print*,'FORCE_FREE_JET: using analytic solution for mu=-1'
        do l=1,mx; do m=1,my; do n=1,mz
          Ax_ext=-2*y(m)*(1-z(n)/sqrt(x(l)**2+y(m)**2+z(n)**2))/(x(l)**2+y(m)**2)/B1
          Ay_ext= 2*x(l)*(1-z(n)/sqrt(x(l)**2+y(m)**2+z(n)**2))/(x(l)**2+y(m)**2)/B1
        enddo; enddo; enddo
!
      else
!
        do l=1,mx; do m=1,my; do n=1,mz
          r2=x(l)**2+y(m)**2
          xi2=r2/(r2+z(n)**2)
          A_phi=hypergeometric2F1((1-mu)/2,(2+mu)/2,2.0,xi2,tol) &
               *sqrt(xi2)*sqrt(r2+z(n)**2)**mu/B1
!
          Ax_ext(l,m,n)=-y(m)*A_phi/sqrt(r2)
          Ay_ext(l,m,n)= x(l)*A_phi/sqrt(r2)
        enddo; enddo; enddo
!
      endif
!
!  calculate external magnetic field
!
      if (lroot.and.ip<=5) &
        print*,'FORCE_FREE_JET: calculating the external magnetic field'
!
      do n=n1,n2
      do m=m1,m2
!        call der(Ay_ext,bb_x,3)
!        bb_ext_pot(:,1)=-bb_x
!        call der(Ax_ext,bb_y,3)
!        bb_ext_pot(:,2)= bb_y
!        call der(Ay_ext,bb_z,1)
!        bb_ext_pot(:,3)= bb_z
!        call der(Ax_ext,bb_z,2)
!        bb_ext_pot(:,3)=bb_ext_pot(:,3)-bb_z
!        call set_global(bb_ext_pot,m,n,'B_ext_pot',nx)
      enddo
      enddo
!
      if (ip<=5) then
        call output(trim(directory)//'/Ax_ext.dat',Ax_ext,1)
        call output(trim(directory)//'/Ay_ext.dat',Ay_ext,1)
      endif
!
    endsubroutine force_free_jet
!***********************************************************************
    subroutine piecew_dipole_aa(ampl,inclaa,f,ivar)
!
!  A field that is vertical uniform for r<R_int, inclined dipolar for
!  r>R_ext, and potential in the shell R_int<r<R_ext.
!  This mimics a neutron star just after the Meissner effect forced the
!  internal field to become vertical (aligned with rotation axis).
!
!  AMPL represents mu/4 pi, where  mu = 1/2 Int rr jj dV  is the
!  magnetic moment of the external dipole field.
!  INCLAA is the inclination of the dipolar field.
!
!  Pencilized in order to minimize memory consumption with all the
!  auxiliary variables used.
!
!  23-jul-05/wolf:coded
!
      real, intent(inout), dimension (mx,my,mz,mfarray) :: f
      real, intent(in) :: ampl,inclaa
      real, dimension (nx) :: r_1_mn,r_2_mn,sigma0,sigma1, r_mn
      real :: fact
      real, dimension(2) :: beta(0:1)
      real, dimension(2,3) :: a(0:1,1:3),b(0:1,1:3)
      integer :: ivar
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
        r_1_mn = 1./max(r_mn,tini)
        r_2_mn = 1./max(r_mn**2,tini)
!
        fact = ampl
        ! beta = [beta_1^0, beta_1^1] combines coefficients for m=0, m=1
        beta =  fact * (/ cos(inclaa), -sin(inclaa)/sqrt(2.) /)
        ! a and b for m=0, m=1 (index 1) and interior, shell, exterior index 2)
        a(0,:) = (/ 1./r_ext**3, 1./r_ext**3                  , 0. /) * beta(0)
        a(1,:) = (/ 0.         , 1./(r_ext**3-r_int**3)       , 0. /) * beta(1)
        !
        b(0,:) = (/ 0.         , 0.                           , 1. /) * beta(0)
        b(1,:) = (/ 0.         , -r_int**3/(r_ext**3-r_int**3), 1. /) * beta(1)
        !
        ! The following could be coded much clearer using elsewhere, but
        ! that is not in F90 (and pgf90 doesn't support F95)
        ! r_in < r < r_ext
        sigma0 = a(0,2)*r_mn + b(0,2)*r_2_mn
        sigma1 = a(1,2)*r_mn + b(1,2)*r_2_mn
        where(r_mn>r_ext) ! r > r_ext
          sigma0 = a(0,3)*r_mn + b(0,3)*r_2_mn
          sigma1 = a(1,3)*r_mn + b(1,3)*r_2_mn
        endwhere
        where(r_mn<r_int) ! r < r_int
          sigma0 = a(0,1)*r_mn + b(0,1)*r_2_mn
          sigma1 = a(1,1)*r_mn + b(1,1)*r_2_mn
        endwhere
        sigma1 = sigma1*sqrt(2.)
        f(l1:l2,m,n,ivar+0) = -sigma0*y(m)*r_1_mn
        f(l1:l2,m,n,ivar+1) =  sigma0*x(l1:l2)*r_1_mn + sigma1*z(n)*r_1_mn
        f(l1:l2,m,n,ivar+2) =                     - sigma1*y(m)*r_1_mn
      enddo
!
    endsubroutine piecew_dipole_aa
!***********************************************************************
    subroutine geo_benchmark_B(f)
!
!  30-june-04/grs: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx) :: theta_mn,ar,atheta,aphi,r_mn,phi_mn
      real :: C_int,C_ext,A_int,A_ext
      integer :: j
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
        theta_mn=acos(spread(z(n),1,nx)/r_mn)
        phi_mn=atan2(spread(y(m),1,nx),x(l1:l2))
!
! calculate ax,ay,az (via ar,atheta,aphi) inside shell (& leave zero outside shell)
!
        do j=1,ninit
           select case (initaa(j))
           case ('geo-benchmark-case1')
              if (lroot .and. imn==1) print*, 'geo_benchmark_B: geo-benchmark-case1'
              C_int=-( -1./63.*r_int**4 + 11./84.*r_int**3*r_ext            &
                     + 317./1050.*r_int**2*r_ext**2                         &
                     - 1./5.*r_int**2*r_ext**2*log(r_int) )
              C_ext=-( -1./63.*r_ext**9 + 11./84.*r_ext**8*r_int            &
                     + 317./1050.*r_ext**7*r_int**2                         &
                     - 1./5.*r_ext**7*r_int**2*log(r_ext) )
              A_int=5./2.*(r_ext-r_int)
              A_ext=5./8.*(r_ext**4-r_int**4)
!
              where (r_mn < r_int)
                ar=C_int*ampl_B0*80.*2.*(3.*sin(theta_mn)**2-2.)*r_mn
                atheta=3.*C_int*ampl_B0*80.*sin(2.*theta_mn)*r_mn
                aphi=ampl_B0*A_int*r_mn*sin(theta_mn)
              endwhere
!
              where (r_mn <= r_ext .and. r_mn >= r_int)
                ar=ampl_B0*80.*2.*(3.*sin(theta_mn)**2-2.)*                 &
                   (   1./36.*r_mn**5 - 1./12.*(r_int+r_ext)*r_mn**4        &
                     + 1./14.*(r_int**2+4.*r_int*r_ext+r_ext**2)*r_mn**3    &
                     - 1./3.*(r_int**2*r_ext+r_int*r_ext**2)*r_mn**2        &
                     - 1./25.*r_int**2*r_ext**2*r_mn                        &
                     + 1./5.*r_int**2*r_ext**2*r_mn*log(r_mn) )
                atheta=-ampl_B0*80.*sin(2.*theta_mn)*                       &
                   (   7./36.*r_mn**5 - 1./2.*(r_int+r_ext)*r_mn**4         &
                     + 5./14.*(r_int**2+4.*r_int*r_ext+r_ext**2)*r_mn**3    &
                     - 4./3.*(r_int**2*r_ext+r_int*r_ext**2)*r_mn**2        &
                     + 2./25.*r_int**2*r_ext**2*r_mn                        &
                     + 3./5.*r_int**2*r_ext**2*r_mn*log(r_mn) )
                aphi=ampl_B0*5./8.*sin(theta_mn)*                           &
                   ( 4.*r_ext*r_mn - 3.*r_mn**2 - r_int**4/r_mn**2 )
              endwhere
!
              where (r_mn > r_ext)
                ar=C_ext*ampl_B0*80.*2.*(3.*sin(theta_mn)**2-2.)/r_mn**4
                atheta=-2.*C_ext*ampl_B0*80.*sin(2.*theta_mn)/r_mn**4
                aphi=ampl_B0*A_ext/r_mn**2*sin(theta_mn)
              endwhere
!
          ! debug checks -- look at a pencil near the centre...
              if (ip<=4 .and. imn==(ny+1)*nz/2) then
                 print*,'r_int,r_ext',r_int,r_ext
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(r_mn), imn, iproc:', &
                      iproc, imn, minval(r_mn), maxval(r_mn)
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(theta_mn), imn, iproc:', &
                      iproc, imn, minval(theta_mn), maxval(theta_mn)
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(phi_mn), imn, iproc:', &
                      iproc, imn, minval(phi_mn), maxval(phi_mn)
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(ar), imn, iproc:', &
                      iproc, imn, minval(ar), maxval(ar)
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(atheta), imn, iproc:', &
                      iproc, imn, minval(atheta), maxval(atheta)
                 write(*,'(a45,2i6,2f15.7)') &
                      'geo_benchmark_B: minmax(aphi), imn, iproc:', &
                      iproc, imn, minval(aphi), maxval(aphi)
              endif
!
           case ('geo-benchmark-case2')
              if (lroot .and. imn==1) print*, 'geo_benchmark_B: geo-benchmark-case2 not yet coded.'
!
           case default
              if (lroot .and. imn==1) print*,'geo_benchmark_B: case not defined!'
              call stop_it("")
           endselect
        enddo
        f(l1:l2,m,n,iax)=sin(theta_mn)*cos(phi_mn)*ar + cos(theta_mn)*cos(phi_mn)*atheta - sin(phi_mn)*aphi
        f(l1:l2,m,n,iay)=sin(theta_mn)*sin(phi_mn)*ar + cos(theta_mn)*sin(phi_mn)*atheta + cos(phi_mn)*aphi
        f(l1:l2,m,n,iaz)=cos(theta_mn)*ar - sin(theta_mn)*atheta
     enddo
!
     if (ip<=14) then
        print*,'geo_benchmark_B: minmax(ax) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iax)),maxval(f(l1:l2,m1:m2,n1:n2,iax))
        print*,'geo_benchmark_B: minmax(ay) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iay)),maxval(f(l1:l2,m1:m2,n1:n2,iay))
        print*,'geo_benchmark_B: minmax(az) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iaz)),maxval(f(l1:l2,m1:m2,n1:n2,iaz))
     endif
!
    endsubroutine geo_benchmark_B
!***********************************************************************
    subroutine eta_rdep(eta_r,geta_r,rdep_profile)
!
!   2-jul-2009/koen: creates an r-dependent resistivity for RFP studies
!
      real, dimension(nx,my) :: eta_r,r2,rmax2,gradr_eta_r
      real, dimension(nx,my,3)  :: geta_r
      character (len=labellen) :: rdep_profile
      integer :: i,j
!
      intent(out) :: eta_r,geta_r
!
!  adapted from etazdep
!  Note the usage of mixed array lengths (nx and my)
!
      select case (rdep_profile)
      case ('schnack89')
        do i=1,nx
        do j=1,my
          r2(i,j)=x(i+l1-1)**2+y(j)**2
          rmax2=1. !value should come from start.in?
        enddo
        enddo
!
!  define eta_r: resistivity profile from Y.L. Ho, S.C. Prager &
!              D.D. Schnack, Phys rev letters vol 62 nr 13 1989
!  and define gradr_eta_r: 1/r *d_r(eta_r))
!
        eta_r = eta*(1+9*(r2/rmax2)**15)**2
        gradr_eta_r= 540*eta*(1+9*(r2/rmax2)**15)*(r2/rmax2)**14/rmax2**0.5
!
!  gradient
!
        do i=1,nx
          geta_r(i,:,1) = x(i+l1-1)*gradr_eta_r(i,:)
        enddo
        do i=1,my
          geta_r(:,i,2) = y(i)*gradr_eta_r(:,i)
        enddo
        geta_r(:,:,3) = 0.
!
!  possibility of constant eta_r (as a test)
!
      case ('constant')
        eta_r=eta
!
!        gradient
        geta_r(:,:,1) = 0.
        geta_r(:,:,2) = 0.
        geta_r(:,:,3) = 0.
      endselect
!
    endsubroutine eta_rdep
!***********************************************************************
    subroutine eta_zdep(eta_z,geta_z,zdep_profile)
!
!  creates a z-dependent resistivity for protoplanetary disk studies
!
!  12-jul-2005/joishi: coded
!
      use General, only:erfcc
!
      real, dimension(mz) :: eta_z,z2
      real, dimension(mz,3) :: geta_z
      character (len=labellen) :: zdep_profile
!      integer :: i
!
      intent(out) :: eta_z,geta_z
!
      select case (zdep_profile)
        case ('fs')
          z2 = z**2.
!  resistivity profile from Fleming & Stone (ApJ 585:908-920)
          eta_z = eta*exp(-z2/2.+sigma_ratio*erfcc(abs(z))/4.)
!
! its gradient:
          geta_z(:,1) = 0.
          geta_z(:,2) = 0.
          geta_z(:,3) = eta_z*(-z-sign(1.,z)*sigma_ratio*exp(-z2)/(2.*sqrt(pi)))
!
        case ('tanh')
!  default to spread gradient over ~5 grid cells.
           if (eta_width == 0.) eta_width = 5.*dz
           eta_z = eta*0.5*(tanh((z + eta_z0)/eta_width) &
             - tanh((z - eta_z0)/eta_width))
!
! its gradient:
           geta_z(:,1) = 0.
           geta_z(:,2) = 0.
           geta_z(:,3) = -eta/(2.*eta_width) * ((tanh((z + eta_z0)/eta_width))**2. &
             - (tanh((z - eta_z0)/eta_width))**2.)
!
      endselect
!
    endsubroutine eta_zdep
!************************************************************************
    subroutine bb_unitvec_shock(f,bb_hat)
!
!  Compute unit vector along the magnetic field.
!  Accurate to 2nd order.
!  Tries to avoid division by zero.
!  Taken from http://nuclear.llnl.gov/CNP/apt/apt/aptvunb.html.
!  If anybody knows a more accurate way of doing this, please modify.
!
!  16-aug-06/tobi: coded
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,3), intent (out) :: bb_hat
!
      !Tobi: Not sure about this value
      real, parameter :: tol=1e-11
!
      real, dimension (mx,3) :: bb,bb2
      real, dimension (mx) :: bb_len,aerr2
      real :: fac
      integer :: j
!
!  Compute magnetic field from vector potential.
!
      bb=0.
!
      if (nxgrid/=1) then
        fac = 1/(2*dx)
        bb(l1-2:l2+2,3) = bb(l1-2:l2+2,3) + fac*( f(l1-1:l2+3,m  ,n  ,iay)   &
                                                - f(l1-3:l2+1,m  ,n  ,iay) )
        bb(l1-2:l2+2,2) = bb(l1-2:l2+2,2) - fac*( f(l1-1:l2+3,m  ,n  ,iaz)   &
                                                - f(l1-3:l2+1,m  ,n  ,iaz) )
      endif
!
      if (nygrid/=1) then
        fac = 1/(2*dy)
        bb(l1-2:l2+2,1) = bb(l1-2:l2+2,1) + fac*( f(l1-2:l2+2,m+1,n  ,iaz)   &
                                                - f(l1-2:l2+2,m-1,n  ,iaz) )
        bb(l1-2:l2+2,3) = bb(l1-2:l2+2,3) - fac*( f(l1-2:l2+2,m+1,n  ,iax)   &
                                                - f(l1-2:l2+2,m-1,n  ,iax) )
      endif
!
      if (nzgrid/=1) then
        fac = 1/(2*dz)
        bb(l1-2:l2+2,2) = bb(l1-2:l2+2,2) + fac*( f(l1-2:l2+2,m  ,n+1,iax)   &
                                                - f(l1-2:l2+2,m  ,n-1,iax) )
        bb(l1-2:l2+2,1) = bb(l1-2:l2+2,1) - fac*( f(l1-2:l2+2,m  ,n+1,iay)   &
                                                - f(l1-2:l2+2,m  ,n-1,iay) )
      endif
!
!  Add external magnetic field.
!
      do j=1,3; bb(:,j) = bb(:,j) + B_ext(j); enddo
!
!  Truncate small components to zero.
!
      bb2 = bb**2
!
      aerr2 = tol**2 * max(sum(bb2,2),1.)
!
      do j=1,3
        where (bb2(:,j) < aerr2)
          bb_hat(:,j) = 0.
        elsewhere
          bb_hat(:,j) = bb(:,j)
        endwhere
      enddo
!
!  Get unit vector.
!
      bb_len = sqrt(sum(bb_hat**2,2))
!
      do j=1,3; bb_hat(:,j) = bb_hat(:,j)/(bb_len+tini); enddo
!
    endsubroutine bb_unitvec_shock
!***********************************************************************
    subroutine input_persistent_magnetic(id,lun,done)
!
!  Read in the stored phase and amplitude for the
!  correction of the Beltrami wave forcing
!
!   5-apr-08/axel: adapted from input_persistent_forcing
!
      integer :: id,lun
      logical :: done
!
      if (id==id_record_MAGNETIC_PHASE) then
        read (lun) phase_beltrami
        done=.true.
        if (lroot) print*, 'input_persistent_magnetic: ', phase_beltrami
      elseif (id==id_record_MAGNETIC_AMPL) then
        read (lun) ampl_beltrami
        done=.true.
        if (lroot) print*, 'input_persistent_magnetic: ', ampl_beltrami
      endif
!
    endsubroutine input_persistent_magnetic
!***********************************************************************
    logical function output_persistent_magnetic(lun)
!
!  Write the stored phase and amplitude for the
!  correction of the Beltrami wave forcing
!
!   5-apr-08/axel: adapted from output_persistent_forcing
!  16-nov-11/MR: changed into logical function to signal I/O errors, I/O error handling introduced
!
      integer :: lun
!
      integer :: iostat
!
      if (lroot.and.ip<14) then
        if (phase_beltrami>=0.) &
            print*,'output_persistent_magnetic: ', &
              phase_beltrami,ampl_beltrami
      endif
!
!  write details
!
      output_persistent_magnetic = .true.
!
      write (lun,IOSTAT=iostat) id_record_MAGNETIC_PHASE
      if (outlog(lun,'write id_record_MAGNETIC_PHASE')) return
      write (lun,IOSTAT=iostat) phase_beltrami
      if (outlog(lun,'write phase_beltrami')) return
      write (lun,IOSTAT=iostat) id_record_MAGNETIC_AMPL
      if (outlog(lun,'write id_record_MAGNETIC_AMPL')) return
      write (lun,IOSTAT=iostat) ampl_beltrami
      if (outlog(lun,'write ampl_beltrami')) return
!
      output_persistent_magnetic = .false.
!
    endfunction output_persistent_magnetic
!***********************************************************************
    subroutine rprint_magnetic(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Diagnostics
!
      integer :: iname,inamex,inamey,inamez,ixy,ixz,irz,inamer,iname_half
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_aphi2m=0; idiag_bphi2m=0; idiag_bpol2m=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'aphi2m',idiag_aphi2m)
        call parse_name(iname,cname(iname),cform(iname),'bphi2m',idiag_bphi2m)
        call parse_name(iname,cname(iname),cform(iname),'bpol2m',idiag_bpol2m)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'nnamexz=',nnamexz
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'iaphi=',iaphi
        write(3,*) 'ibphi=',ibphi
      endif
!
    endsubroutine rprint_magnetic
!***********************************************************************
endmodule Magnetic
