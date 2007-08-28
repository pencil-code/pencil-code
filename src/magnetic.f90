! $Id: magnetic.f90,v 1.453 2007-08-28 01:07:12 wlyra Exp $
!  This modules deals with all aspects of magnetic fields; if no
!  magnetic fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  magnetically relevant subroutines listed in here.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lmagnetic = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED aa,a2,aij,bb,bbb,ab,uxb,b2,bij,del2a,graddiva,jj
! PENCILS PROVIDED j2,jb,va2,jxb,jxbr,ub,uxb,uxb2,uxj,beta,uga
! PENCILS PROVIDED djuidjbi,jo,ujxb,oxu,oxuxb,jxbxb,jxbrxb
! PENCILS PROVIDED glnrhoxb,del4a,del6a,oxj,diva,jij,sj,ss12
! PENCILS PROVIDED mf_EMF, mf_EMFdotB
!
!***************************************************************

module Magnetic

  use Cdata
  use Cparam
  use Messages

  implicit none

  include 'magnetic.h'
!
! Slice precalculation buffers
!
  real, target, dimension (nx,ny,3) :: bb_xy,jj_xy
  real, target, dimension (nx,ny,3) :: bb_xy2,jj_xy2
  real, target, dimension (nx,nz,3) :: bb_xz,jj_xz
  real, target, dimension (ny,nz,3) :: bb_yz,jj_yz
!
  real, target, dimension (nx,ny) :: b2_xy,jb_xy
  real, target, dimension (nx,ny) :: b2_xy2,jb_xy2
  real, target, dimension (ny,nz) :: b2_yz,jb_yz
  real, target, dimension (nx,nz) :: b2_xz,jb_xz
!
  real, dimension (mx,my) :: alpha_input
!
! Parameters
!
  integer, parameter :: nresi_max=4
!
  real, dimension (ninit) :: amplaa=0.0,kx_aa=1.,ky_aa=1.,kz_aa=1., phasey_aa=0.
  character (len=labellen), dimension(ninit) :: initaa='nothing'
  character (len=labellen) :: borderaa='nothing'
  character (len=labellen), dimension(nresi_max) :: iresistivity=''
  character (len=labellen) :: Omega_profile='nothing',alpha_profile='nothing'
  ! input parameters
  complex, dimension(3) :: coefaa=(/0.,0.,0./), coefbb=(/0.,0.,0./)
  real, dimension(3) :: B_ext=(/0.,0.,0./),B1_ext,B_ext_tmp,eta_aniso_hyper3
  real, dimension(3) :: axisr1=(/0,0,1/),dispr1=(/0.,0.5,0./)
  real, dimension(3) :: axisr2=(/1,0,0/),dispr2=(/0.,-0.5,0./)
  real, dimension(nx,3) :: uxbb !(temporary)
  real, target :: zmode=1. !(temporary)
  real :: fring1=0.,Iring1=0.,Rring1=1.,wr1=0.3
  real :: fring2=0.,Iring2=0.,Rring2=1.,wr2=0.3
  real :: radius=.1,epsilonaa=1e-2,widthaa=.5,x0aa=0.,z0aa=0.
  real :: by_left=0.,by_right=0.,bz_left=0.,bz_right=0.
  real :: ABC_A=1.,ABC_B=1.,ABC_C=1.
  real :: bthresh=0.,bthresh_per_brms=0.,brms=0.,bthresh_scl=1.
  real :: eta_shock=0.
  real :: rhomin_jxb=0.,va2max_jxb=0.
  real :: omega_Bz_ext=0.
  real :: mu_r=-0.5 !(still needed for backwards compatibility)
  real :: mu_ext_pot=-0.5,inclaa=0.
  real :: rescale_aa=1.
  real :: ampl_B0=0.,D_smag=0.17,B_ext21,B_ext11
  real :: Omega_ampl=0.
  real :: rmode=1.,rm_int=0.,rm_ext=0.
  real :: nu_ni=0.,nu_ni1,hall_term=0.
  real :: alpha_effect=0.,alpha_quenching=0.,delta_effect=0.,meanfield_etat=0.
  real :: displacement_gun=0.
  real :: pertamplaa=0., beta_const=1.0
  real :: initpower_aa=0.,cutoff_aa=0.,brms_target=1.,rescaling_fraction=1.
  integer :: nbvec,nbvecmax=nx*ny*nz/4,va2power_jxb=5
  integer :: N_modes_aa=1
  integer :: iglobal_bx_ext=0, iglobal_by_ext=0, iglobal_bz_ext=0
  integer :: iglobal_jx_ext=0, iglobal_jy_ext=0, iglobal_jz_ext=0
  integer :: iglobal_ex_ext=0, iglobal_ey_ext=0, iglobal_ez_ext=0
  logical :: lpress_equil=.false., lpress_equil_via_ss=.false.
  logical :: llorentzforce=.true.,linduction=.true.
  logical :: lresi_eta_const=.false.
  logical :: lresi_etaSS=.false.
  logical :: lresi_hyper2=.false.
  logical :: lresi_hyper3=.false.
  logical :: lresi_hyper3_strict=.false.
  logical :: lresi_zdep=.false.
  logical :: lresi_hyper3_aniso=.false.
  logical :: lresi_eta_shock=.false.
  logical :: lresi_eta_shock_perp=.false.
  logical :: lresi_shell=.false.
  logical :: lresi_smagorinsky=.false.
  logical :: lresi_smagorinsky_cross=.false.
  logical, dimension (3) :: lfrozen_bb_bot=(/.false.,.false.,.false./)
  logical, dimension (3) :: lfrozen_bb_top=(/.false.,.false.,.false./)
  logical :: reinitalize_aa=.false., lohmic_heat=.true.
  logical :: lB_ext_pot=.false.
  logical :: lforce_free_test=.false.
  logical :: lmeanfield_theory=.false.,lOmega_effect=.false.
  logical :: lmeanfield_noalpm=.false.
  logical :: lgauss=.false.
  logical :: lbb_as_aux=.false.,ljj_as_aux=.false.
  character (len=labellen) :: pertaa='zero'

  namelist /magnetic_init_pars/ &
       B_ext, lohmic_heat, &
       fring1,Iring1,Rring1,wr1,axisr1,dispr1, &
       fring2,Iring2,Rring2,wr2,axisr2,dispr2, &
       radius,epsilonaa,x0aa,z0aa,widthaa, &
       by_left,by_right,bz_left,bz_right, &
       initaa,amplaa,kx_aa,ky_aa,kz_aa,phasey_aa,coefaa,coefbb, &
       inclaa,lpress_equil,lpress_equil_via_ss,mu_r, &
       mu_ext_pot,lB_ext_pot,lforce_free_test, &
       ampl_B0,initpower_aa,cutoff_aa,N_modes_aa, &
       rmode,zmode,rm_int,rm_ext,lgauss, &
       lbb_as_aux,ljj_as_aux,beta_const

  ! run parameters
  real :: eta=0.,eta_hyper2=0.,eta_hyper3=0.,height_eta=0.,eta_out=0.
  real :: eta_int=0.,eta_ext=0.,wresistivity=.01
  real :: tau_aa_exterior=0.
  real :: sigma_ratio=1.,eta_width=0.,eta_z0=1.
  real :: alphaSSm=0.
  real :: k1_ff=1.,ampl_ff=1.,swirl=1.
  real :: k1x_ff=1.,k1y_ff=1.,k1z_ff=1.
  real :: inertial_length=0.,linertial_2
  real, dimension(mz) :: eta_z
  real, dimension(mz,3) :: geta_z
  logical :: lfreeze_aint=.false.,lfreeze_aext=.false.
  logical :: lweyl_gauge=.false.
  logical :: lupw_aa=.false.
  logical :: lforcing_continuous_aa=.false.
  logical :: lelectron_inertia=.false.
  character (len=labellen) :: zdep_profile='fs'
  character (len=labellen) :: iforcing_continuous_aa='fixed_swirl'

  namelist /magnetic_run_pars/ &
       eta,eta_hyper2,eta_hyper3,B_ext,omega_Bz_ext,nu_ni,hall_term, &
       lmeanfield_theory,alpha_effect,alpha_quenching,delta_effect, &
       lmeanfield_noalpm,alpha_profile, &
       meanfield_etat, lohmic_heat, &
       height_eta,eta_out,tau_aa_exterior, &
       kx_aa,ky_aa,kz_aa,phasey_aa,ABC_A,ABC_B,ABC_C, &
       lforcing_continuous_aa,iforcing_continuous_aa, &
       k1_ff,ampl_ff,swirl,radius, &
       k1x_ff,k1y_ff,k1z_ff, &
       bthresh,bthresh_per_brms, &
       iresistivity,lweyl_gauge,lupw_aa, &
       alphaSSm, &
       eta_int,eta_ext,eta_shock,wresistivity, &
       rhomin_jxb,va2max_jxb,va2power_jxb,llorentzforce,linduction, &
       reinitalize_aa,rescale_aa,lB_ext_pot, &
       displacement_gun, &
       pertaa,pertamplaa,D_smag,brms_target,rescaling_fraction, &
       lOmega_effect,Omega_profile,Omega_ampl,lfreeze_aint,lfreeze_aext, &
       sigma_ratio,zdep_profile,eta_width,eta_z0, &
       borderaa,eta_aniso_hyper3, &
       lelectron_inertia,inertial_length, &
       lbb_as_aux,ljj_as_aux

  ! diagnostic variables (need to be consistent with reset list below)
  integer :: idiag_b2m=0        ! DIAG_DOC: $\left<\Bv^2\right>$
  integer :: idiag_bm2=0        ! DIAG_DOC: $\max(\Bv^2)$
  integer :: idiag_j2m=0        ! DIAG_DOC: $\left<\jv^2\right>$
  integer :: idiag_jm2=0        ! DIAG_DOC: $\max(\jv^2)$
  integer :: idiag_abm=0        ! DIAG_DOC: $\left<\Av\cdot\Bv\right>$
  integer :: idiag_jbm=0        ! DIAG_DOC: $\left<\jv\cdot\Bv\right>$
  integer :: idiag_ubm=0        ! DIAG_DOC: $\left<\uv\cdot\Bv\right>$
  integer :: idiag_epsM=0       ! DIAG_DOC: $\left<2\eta\mu_0\jv^2\right>$
  integer :: idiag_bxpt=0       ! DIAG_DOC: $B_x(x_0,y_0,z_0,t)$
  integer :: idiag_bypt=0       ! DIAG_DOC: $B_y(x_0,y_0,z_0,t)$
  integer :: idiag_bzpt=0       ! DIAG_DOC: $B_z(x_0,y_0,z_0,t)$
  integer :: idiag_epsM_LES=0   ! DIAG_DOC:
  integer :: idiag_aybym2=0     ! DIAG_DOC:
  integer :: idiag_exaym2=0     ! DIAG_DOC:
  integer :: idiag_exjm2=0      ! DIAG_DOC:
  integer :: idiag_brms=0       ! DIAG_DOC: $\left<\Bv^2\right>^{1/2}$
  integer :: idiag_bmax=0       ! DIAG_DOC: $\max(|\Bv|)$
  integer :: idiag_jrms=0       ! DIAG_DOC: $\left<\jv^2\right>^{1/2}$
  integer :: idiag_jmax=0       ! DIAG_DOC: $\max(|\jv|)$
  integer :: idiag_vArms=0      ! DIAG_DOC: $\left<\Bv^2/\varrho\right>^{1/2}$
  integer :: idiag_vAmax=0      ! DIAG_DOC: $\max(\Bv^2/\varrho)^{1/2}$
  integer :: idiag_dtb=0        ! DIAG_DOC: $\delta t / [c_{\delta t}\,\delta x
                                ! DIAG_DOC:   /v_{\rm A,max}]$
                                ! DIAG_DOC:   \quad(time step relative to
                                ! DIAG_DOC:   Alfv{\'e}n time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_dteta=0      ! DIAG_DOC: $\delta t/[c_{\delta t,{\rm v}}\,
                                ! DIAG_DOC:   \delta x^2/\eta_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to
                                ! DIAG_DOC:   resistive time step;
                                ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_arms=0       ! DIAG_DOC:
  integer :: idiag_amax=0       ! DIAG_DOC:
  integer :: idiag_beta1m=0     ! DIAG_DOC: $\left<\Bv^2/(2\mu_0 p)\right>$
                                ! DIAG_DOC:   \quad(mean inverse plasma beta)
  integer :: idiag_beta1max=0   ! DIAG_DOC: $\max[\Bv^2/(2\mu_0 p)]$
                                ! DIAG_DOC:   \quad(maximum inverse plasma beta)
  integer :: idiag_bxm=0        ! DIAG_DOC:
  integer :: idiag_bym=0        ! DIAG_DOC:
  integer :: idiag_bzm=0        ! DIAG_DOC:
  integer :: idiag_bx2m=0       ! DIAG_DOC:
  integer :: idiag_by2m=0       ! DIAG_DOC:
  integer :: idiag_bz2m=0       ! DIAG_DOC:
  integer :: idiag_bxbym=0      ! DIAG_DOC:
  integer :: idiag_bxbzm=0      ! DIAG_DOC:
  integer :: idiag_bybzm=0      ! DIAG_DOC:
  integer :: idiag_djuidjbim=0  ! DIAG_DOC:
  integer :: idiag_bxbymz=0     ! DIAG_DOC:
  integer :: idiag_bxbzmz=0     ! DIAG_DOC:
  integer :: idiag_bybzmz=0     ! DIAG_DOC:
  integer :: idiag_b2mz=0       ! DIAG_DOC:
  integer :: idiag_bxmz=0       ! DIAG_DOC:
  integer :: idiag_bymz=0       ! DIAG_DOC:
  integer :: idiag_bzmz=0       ! DIAG_DOC:
  integer :: idiag_bmx=0        ! DIAG_DOC: $\left<\left<\Bv\right>_{yz}^2
                                ! DIAG_DOC:   \right>^{1/2}$
                                ! DIAG_DOC:   \quad(energy of $yz$-averaged
                                ! DIAG_DOC:   mean field)
  integer :: idiag_bmy=0        ! DIAG_DOC: $\left<\left<\Bv\right>_{xz}^2
                                ! DIAG_DOC:   \right>^{1/2}$
                                ! DIAG_DOC:   \quad(energy of $xz$-averaged
                                ! DIAG_DOC:   mean field)
  integer :: idiag_bmz=0        ! DIAG_DOC: $\left<\left<\Bv\right>_{xy}^2
                                ! DIAG_DOC:   \right>^{1/2}$
                                ! DIAG_DOC:   \quad(energy of $xy$-averaged
                                ! DIAG_DOC:   mean field)
  integer :: idiag_bx2mz=0      ! DIAG_DOC:
  integer :: idiag_by2mz=0      ! DIAG_DOC:
  integer :: idiag_bz2mz=0      ! DIAG_DOC:
  integer :: idiag_bxmxy=0      ! DIAG_DOC:
  integer :: idiag_bymxy=0      ! DIAG_DOC:
  integer :: idiag_bzmxy=0      ! DIAG_DOC:
  integer :: idiag_bxmxz=0      ! DIAG_DOC:
  integer :: idiag_bymxz=0      ! DIAG_DOC:
  integer :: idiag_bzmxz=0      ! DIAG_DOC:
  integer :: idiag_uxbm=0       ! DIAG_DOC:
  integer :: idiag_oxuxbm=0     ! DIAG_DOC:
  integer :: idiag_jxbxbm=0     ! DIAG_DOC:
  integer :: idiag_gpxbm=0      ! DIAG_DOC:
  integer :: idiag_uxDxuxbm=0   ! DIAG_DOC:
  integer :: idiag_jbmphi=0     ! DIAG_DOC:
  integer :: idiag_b3b21m=0     ! DIAG_DOC:
  integer :: idiag_b1b32m=0     ! DIAG_DOC:
  integer :: idiag_b2b13m=0     ! DIAG_DOC:
  integer :: idiag_EMFdotBm=0   ! DIAG_DOC:
  integer :: idiag_udotxbm=0    ! DIAG_DOC:
  integer :: idiag_uxbdotm=0    ! DIAG_DOC:
  integer :: idiag_uxbmx=0      ! DIAG_DOC:
  integer :: idiag_uxbmy=0      ! DIAG_DOC:
  integer :: idiag_uxbmz=0      ! DIAG_DOC:
  integer :: idiag_uxjm=0       ! DIAG_DOC:
  integer :: idiag_brmphi=0     ! DIAG_DOC:
  integer :: idiag_bpmphi=0     ! DIAG_DOC:
  integer :: idiag_bzmphi=0     ! DIAG_DOC:
  integer :: idiag_b2mphi=0     ! DIAG_DOC:
  integer :: idiag_uxbrmphi=0   ! DIAG_DOC:
  integer :: idiag_uxbpmphi=0   ! DIAG_DOC:
  integer :: idiag_uxbzmphi=0   ! DIAG_DOC:
  integer :: idiag_ujxbm=0      ! DIAG_DOC:
  integer :: idiag_jxbrmphi=0   ! DIAG_DOC:
  integer :: idiag_jxbpmphi=0   ! DIAG_DOC:
  integer :: idiag_jxbzmphi=0   ! DIAG_DOC:
  integer :: idiag_armphi=0     ! DIAG_DOC:
  integer :: idiag_apmphi=0     ! DIAG_DOC:
  integer :: idiag_azmphi=0     ! DIAG_DOC:
  integer :: idiag_uxBrms=0     ! DIAG_DOC:
  integer :: idiag_Bresrms=0    ! DIAG_DOC:
  integer :: idiag_Rmrms=0      ! DIAG_DOC:
  integer :: idiag_jfm=0        ! DIAG_DOC:
  integer :: idiag_brbpmr=0     ! DIAG_DOC:
  integer :: idiag_vA2m=0       ! DIAG_DOC:
  integer :: idiag_b2mr=0       ! DIAG_DOC:
  integer :: idiag_brmr=0       ! DIAG_DOC:
  integer :: idiag_bpmr=0       ! DIAG_DOC:
  integer :: idiag_bzmr=0       ! DIAG_DOC:
  integer :: idiag_armr=0       ! DIAG_DOC:
  integer :: idiag_apmr=0       ! DIAG_DOC:
  integer :: idiag_azmr=0       ! DIAG_DOC:
  integer :: idiag_bxmx=0       ! DIAG_DOC:
  integer :: idiag_bymy=0       ! DIAG_DOC:
  integer :: idiag_mflux_x=0    ! DIAG_DOC:
  integer :: idiag_mflux_y=0    ! DIAG_DOC:
  integer :: idiag_mflux_z=0    ! DIAG_DOC:
  integer :: idiag_bmxy_rms=0   ! DIAG_DOC: $\sqrt{[\left<b_x\right>_z(x,y)]^2 + 
                                ! DIAG_DOC: [\left<b_x\right>_z(x,y)]^2 +
                                ! DIAG_DOC: [\left<b_x>_z(x,y)\right>]^2} $ 
  contains

!***********************************************************************
    subroutine register_magnetic()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaa, etc; increase nvar accordingly
!
!  1-may-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_aa called twice')
      first = .false.
!
      iaa = nvar+1              ! indices to access aa
      iax = iaa
      iay = iaa+1
      iaz = iaa+2
      nvar = nvar+3             ! added 3 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_magnetic: nvar = ', nvar
        print*, 'register_magnetic: iaa,iax,iay,iaz = ', iaa,iax,iay,iaz
      endif
!
!  Put variable names in array
!
      varname(iax) = 'ax'
      varname(iay) = 'ay'
      varname(iaz) = 'az'
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: magnetic.f90,v 1.453 2007-08-28 01:07:12 wlyra Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_magnetic: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',aa $'
          if (nvar == mvar) write(4,*) ',aa'
        else
          write(4,*) ',aa $'
        endif
        write(15,*) 'aa = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic(f,lstarting)
!
!  Perform any post-parameter-read initialization
!
!  24-nov-02/tony: dummy routine - nothing to do at present
!  20-may-03/axel: reinitalize_aa added
!
      use Cdata
      use Messages, only: fatal_error
      use FArrayManager
      use SharedVariables
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
      integer :: i,ierr
!
!  Precalculate 1/mu (moved here from register.f90)
!
      mu01=1./mu0
!
!  Precalculate 1/inertial_length^2
!
      if (inertial_length /= 0.) then
        linertial_2 = inertial_length**(-2)
      else
        linertial_2 = 0.
        ! make sure not to use this value by checking that
        ! (inertial_length /= 0.)...
      endif
!
!  Precalculate 1/nu_ni
!
      if (nu_ni /= 0.) then
        nu_ni1=1./nu_ni
      else
        nu_ni1=0.
      endif
!
!  calculate B_ext21 = 1/B_ext**2 and the unit vector B1_ext = B_ext/|B_ext|
!
      B_ext21=B_ext(1)**2+B_ext(2)**2+B_ext(3)**2
      if (B_ext21/=0.) then
        B_ext21=1./B_ext21
      else
        B_ext21=1.
      endif
      B_ext11=sqrt(B_ext21)
      B1_ext=B_ext*B_ext11
!
!  set to zero and then rescale the magnetic field
!  (in future, could call something like init_aa_simple)
!
      if (reinitalize_aa) then
        f(:,:,:,iax:iaz)=rescale_aa*f(:,:,:,iax:iaz)
        call pert_aa(f)
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
      lresi_hyper3_strict=.false.
      lresi_hyper3_aniso=.false.
      lresi_eta_shock=.false.
      lresi_eta_shock_perp=.false.
      lresi_smagorinsky=.false.
      lresi_smagorinsky_cross=.false.

      do i=1,nresi_max
        select case (iresistivity(i))
        case ('eta-const')
          if (lroot) print*, 'resistivity: constant eta'
          lresi_eta_const=.true.
        case ('etaSS')
          if (lroot) print*, 'resistivity: etaSS (Shakura-Sunyaev)'
          lresi_etaSS=.true.
        case('hyper2')
          if (lroot) print*, 'resistivity: hyper2'
          lresi_hyper2=.true.
        case('hyper3')
          if (lroot) print*, 'resistivity: hyper3'
          lresi_hyper3=.true.
        case('hyper3_strict')
          if (lroot) print*, 'resistivity: strict hyper3 with positive definite heating rate'
          lresi_hyper3_strict=.true.
        case('zdep')
          if (lroot) print*, 'resistivity: z-dependent'
          lresi_zdep=.true.
          call eta_zdep(eta_z,geta_z,zdep_profile)
        case('hyper3-aniso')
          if (lroot) print*, 'resistivity: hyper3_aniso'
          lresi_hyper3_aniso=.true.
        case('shell')
          if (lroot) print*, 'resistivity: shell'
          lresi_shell=.true.
        case ('shock')
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
        if (lresi_eta_const.and.eta==0.0) &
            call warning('initialize_magnetic', &
            'Resistivity coefficient eta is zero!')
        if (lresi_hyper2.and.eta_hyper2==0.0) &
            call fatal_error('initialize_magnetic', &
            'Resistivity coefficient eta_hyper2 is zero!')
        if (lresi_hyper3.and.eta_hyper3==0.0) &
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
!  check for alpha profile
!
      if (alpha_profile=='read') then
        print*,'read alpha profile'
        open(1,file='alpha_input.dat',form='unformatted')
        read(1) alpha_input
        close(1)
      endif
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
!  Register global variables
!
      if (any(initaa=='hydrostatic_magnetic')) then
        call farray_register_global('iglobal_by_ext',iglobal_by_ext)
        call farray_register_global('iglobal_jx_ext',iglobal_jx_ext)
      endif
!
      if (any(initaa=='Alfven-zconst')) then    
        call put_shared_variable('zmode',zmode,ierr)
        if (ierr/=0) call fatal_error('initialize_magnetic',&
             'there was a problem when sharing zmode')
      endif
!
      if (NO_WARN) print*, lstarting
!
    endsubroutine initialize_magnetic
!***********************************************************************
    subroutine init_aa(f,xx,yy,zz)
!
!  initialise magnetic field; called from start.f90
!  AB: maybe we should here call different routines (such as rings)
!  AB: and others, instead of accummulating all this in a huge routine.
!  We have an init parameter (initaa) to stear magnetic i.c. independently.
!
!   7-nov-2001/wolf: coded
!
      use Cdata
      use EquationOfState
      use FArrayManager
      use Gravity, only: gravz
      use Initcond
      use Mpicomm
      use SharedVariables
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz,tmp,prof
!
      real, dimension (nx,3) :: bb
      real, dimension (nx) :: b2,fact,rho
      real :: beq2, scaleH
      real, pointer :: nu_epicycle
      integer :: j, ierr
!
      do j=1,ninit

        select case(initaa(j))

        case('nothing'); if(lroot .and. j==1) print*,'init_uu: nothing'
        case('zero', '0'); f(:,:,:,iax:iaz) = 0.
        case('rescale'); f(:,:,:,iax:iaz)=amplaa(j)*f(:,:,:,iax:iaz)
        case('mode'); call modev(amplaa(j),coefaa,f,iaa,kx_aa(j),ky_aa(j),kz_aa(j),xx,yy,zz)
        case('modeb'); call modeb(amplaa(j),coefbb,f,iaa,kx_aa(j),ky_aa(j),kz_aa(j),xx,yy,zz)
        case('const_lou'); call const_lou(amplaa(j),f,iaa,xx,yy,zz)
        case('power_randomphase')
          call power_randomphase(amplaa(j),initpower_aa,cutoff_aa,f,iax,iaz)
        case('random-isotropic-KS')
          call random_isotropic_KS(amplaa(j),initpower_aa,cutoff_aa,f,iax,iaz,N_modes_aa)
        case('gaussian-noise'); call gaunoise(amplaa(j),f,iax,iaz)
        case('gaussian-noise-rprof')
          tmp=sqrt(xx**2+yy**2+zz**2)
          call gaunoise_rprof(amplaa(j),tmp,prof,f,iax,iaz)
        case('Beltrami-x', '11'); call beltrami(amplaa(j),f,iaa,KX=kx_aa(j))
        case('Beltrami-y', '12'); call beltrami(amplaa(j),f,iaa,KY=ky_aa(j))
        case('Beltrami-z', '1');  call beltrami(amplaa(j),f,iaa,KZ=kz_aa(j))
        case('propto-ux'); call wave_uu(amplaa(j),f,iaa,kx=kx_aa(j))
        case('propto-uy'); call wave_uu(amplaa(j),f,iaa,ky=ky_aa(j))
        case('propto-uz'); call wave_uu(amplaa(j),f,iaa,kz=kz_aa(j))
        case('diffrot'); call diffrot(amplaa(j),f,iay,xx,yy,zz)
        case('hor-tube'); call htube(amplaa(j),f,iax,iaz,xx,yy,zz,radius,epsilonaa)
        case('hor-fluxlayer'); call hfluxlayer(amplaa(j),f,iaa,xx,yy,zz,z0aa,widthaa)
        case('ver-fluxlayer'); call vfluxlayer(amplaa(j),f,iaa,xx,yy,zz,x0aa,widthaa)
        case('mag-support'); call magsupport(amplaa(j),f,zz,gravz,cs0,rho0)
        case('pattern-xy'); call vecpatternxy(amplaa(j),f,iaa)
        case('arcade-x'); call arcade_x(amplaa(j),f,iaa,xx,yy,zz,kx_aa(j),kz_aa(j))
        case('halfcos-Bx'); call halfcos_x(amplaa(j),f,iaa,xx,yy,zz)
        case('uniform-Bx'); call uniform_x(amplaa(j),f,iaa,xx,yy,zz)
        case('uniform-By'); call uniform_y(amplaa(j),f,iaa,xx,yy,zz)
        case('uniform-Bz'); call uniform_z(amplaa(j),f,iaa,xx,yy,zz)
        case('uniform-Bphi'); call uniform_phi(amplaa(j),f,iaa,xx,yy,zz)
        case('Bz(x)', '3'); call vfield(amplaa(j),f,iaa,xx)
        case('vfield2'); call vfield2(amplaa(j),f,iaa,xx)
        case('vecpatternxy'); call vecpatternxy(amplaa(j),f,iaa)
        case('xjump'); call bjump(f,iaa,by_left,by_right,bz_left,bz_right,widthaa,'x')
        case('fluxrings', '4'); call fluxrings(amplaa(j),f,iaa,xx,yy,zz)
        case('sinxsinz'); call sinxsinz(amplaa(j),f,iaa,kx_aa(j),ky_aa(j),kz_aa(j))
        case('sin2xsin2y'); call sin2x_sin2y_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case('cosxcosy'); call cosx_cosy_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case('sinxsiny'); call sinx_siny_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case('xsiny'); call x_siny_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case('x1siny'); call x1_siny_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.,phasey_aa(j))
        case('sinxcosz'); call sinx_siny_cosz(amplaa(j),f,iay,kx_aa(j),ky_aa(j),kz_aa(j))
        case('sinycosz'); call cosx_siny_cosz(amplaa(j),f,iax,kx_aa(j),ky_aa(j),0.)
        case('cosysinz'); call cosy_sinz(amplaa(j),f,iax,ky_aa(j),kz_aa(j))
        case('x3cosycosz'); call x3_cosy_cosz(amplaa(j),f,iax,ky_aa(j),kz_aa(j))
        case('Ax=cosysinz'); call cosy_sinz(amplaa(j),f,iax,ky_aa(j),kz_aa(j))
        case('magnetogram'); call mdi_init(f)
        case('cosxcoscosy'); call cosx_coscosy_cosz(amplaa(j),f,iaz,kx_aa(j),ky_aa(j),0.)
        case('crazy', '5'); call crazy(amplaa(j),f,iaa)
        case('sinwave-x'); call sinwave(amplaa(j),f,iaa,kx=kx_aa(j))
        case('linear-zx'); f(:,:,:,iay)=-.5*amplaa(j)*zz**2/Lxyz(3)
        case('Alfven-x'); call alfven_x(amplaa(j),f,iuu,iaa,ilnrho,xx,kx_aa(j))
        case('Alfven-y'); call alfven_y(amplaa(j),f,iuu,iaa,yy,ky_aa(j),mu0)
        case('Alfven-z'); call alfven_z(amplaa(j),f,iuu,iaa,zz,kz_aa(j),mu0)
        case('Alfven-xy'); call alfven_xy(amplaa(j),f,iuu,iaa,xx,yy,kx_aa(j),ky_aa(j))
        case('Alfven-xz'); call alfven_xz(amplaa(j),f,iuu,iaa,xx,zz,kx_aa(j),kz_aa(j))
        case('Alfven-rphi'); call alfven_rphi(amplaa(j),f,xx,yy,rmode)
        case('Alfven-zconst'); call alfven_zconst(f,xx,yy)
        case('Alfven-rz'); call alfven_rz(amplaa(j),f,xx,yy,rmode)
        case('Alfvenz-rot'); call alfvenz_rot(amplaa(j),f,iuu,iaa,zz,kz_aa(j),Omega)
        case('Alfvenz-rot-shear'); call alfvenz_rot_shear(amplaa(j),f,iuu,iaa,zz,kz_aa(j),Omega)
        case('piecewise-dipole'); call piecew_dipole_aa (amplaa(j),inclaa,f,iaa,xx,yy,zz)
        case('tony-nohel')
          f(:,:,:,iay) = amplaa(j)/kz_aa(j)*cos(kz_aa(j)*2.*pi/Lz*zz)
        case('tony-nohel-yz')
          f(:,:,:,iay) = amplaa(j)/kx_aa(j)*sin(kx_aa(j)*2.*pi/Lx*xx)
        case('tony-hel-xy')
          f(:,:,:,iax) = amplaa(j)/kz_aa(j)*sin(kz_aa(j)*2.*pi/Lz*zz)
          f(:,:,:,iay) = amplaa(j)/kz_aa(j)*cos(kz_aa(j)*2.*pi/Lz*zz)
        case('tony-hel-yz')
          f(:,:,:,iay) = amplaa(j)/kx_aa(j)*sin(kx_aa(j)*2.*pi/Lx*xx)
          f(:,:,:,iaz) = amplaa(j)/kx_aa(j)*cos(kx_aa(j)*2.*pi/Lx*xx)
        case('force-free-jet')
          lB_ext_pot=.true.
          call force_free_jet(mu_ext_pot,xx,yy,zz)
        case('Alfven-circ-x')
          !
          !  circularly polarised Alfven wave in x direction
          !
          if (lroot) print*,'init_aa: circular Alfven wave -> x'
          f(:,:,:,iay) = amplaa(j)/kx_aa(j)*sin(kx_aa(j)*xx)
          f(:,:,:,iaz) = amplaa(j)/kx_aa(j)*cos(kx_aa(j)*xx)
        case('geo-benchmark-case1','geo-benchmark-case2'); call geo_benchmark_B(f)
        case('hydrostatic_magnetic')
          scaleH=1.0*sqrt(1+1/beta_const)
          print*, 'init_aa: hydrostatic_magnetic: scaleH=', scaleH
          do m=m1,m2; do n=n1,n2
            if (ldensity_nolog) then
              f(l1:l2,m,n,ilnrho) = rho0*exp(-0.5*(z(n)/scaleH)**2)
              rho=f(l1:l2,m,n,ilnrho)
            else
              f(l1:l2,m,n,ilnrho) = alog(rho0)-0.5*(z(n)/scaleH)**2
              rho=exp(f(l1:l2,m,n,ilnrho))
            endif
            f(l1:l2,m,n,iglobal_by_ext)= &
                sqrt(2*mu0*1.0**2*rho/beta_const)
            f(l1:l2,m,n,iglobal_jx_ext)= &
                sqrt(0.5*mu0*1.0**2*beta_const)*sqrt(rho)*z(n)/scaleH**2
          enddo; enddo
        case('hydrostatic_disk')
          call get_shared_variable('nu_epicycle',nu_epicycle,ierr)
          if (ierr/=0) then
            if (lroot) print*, 'init_aa: '// &
                'there was a problem when getting nu_epicycle!'
            call fatal_error('init_aa','')
          endif
          if (lroot) &
            print*, 'init_aa: hydrostatic_disk: '// &
                'fetched shared variable nu_epicycle=', nu_epicycle
          ! This assumes an isothermal equation of state
          scaleH = (cs0/nu_epicycle)*sqrt(1+1/beta_const)
          print*, 'init_aa: hydrostatic_disk: scaleH=', scaleH
          do n=n1,n2; do m=m1,m2
            if (ldensity_nolog) then
              f(l1:l2,m,n,ilnrho) = rho0*exp(-0.5*(z(n)/scaleH)**2)
            else
              f(l1:l2,m,n,ilnrho) = alog(rho0)-0.5*(z(n)/scaleH)**2
            endif
            f(l1:l2,m,n,iax) = scaleH*sqrt(2*mu0*pi*rho0*cs20/beta_const) &
                              *erfunc(0.5*z(n)/scaleH)
          enddo; enddo
        case('hydrostatic_disk-tanh')
! requires gravz_profile='tanh' and zref=2*scaleH to work!        
          scaleH = cs20*(1+1/beta_const)/1.0
          print*, 'init_aa: hydrostatic_disk-tanh: scaleH=', scaleH
          do n=n1,n2; do m=m1,m2
            if (ldensity_nolog) then
              f(l1:l2,m,n,ilnrho) = rho0*exp(-2*alog(cosh(z(n)/(2*scaleH))))
            else
              f(l1:l2,m,n,ilnrho) = alog(rho0)-2*alog(cosh(z(n)/(2*scaleH)))
            endif
            f(l1:l2,m,n,iax) = sqrt(2*mu0*rho0*cs20/beta_const) * &
                              4*scaleH*atan(tanh(z(n)/(4*scaleH)))
          enddo; enddo
        case('hydrostatic_disk_sinx')
          call get_shared_variable('nu_epicycle',nu_epicycle,ierr)
          if (ierr/=0) then
            if (lroot) print*, 'init_aa: '// &
                'there was a problem when getting nu_epicycle!'
            call fatal_error('init_aa','')
          endif
          if (lroot) &
            print*, 'init_aa: hydrostatic_disk_sinx: '// &
                'fetched shared variable nu_epicycle=', nu_epicycle
          ! This assumes an isothermal equation of state
          scaleH = (cs0/nu_epicycle)*sqrt(1+1/beta_const)
          print*, 'init_aa: hydrostatic_disk_sinx: scaleH=', scaleH
          do n=n1,n2; do m=m1,m2
            if (ldensity_nolog) then
              f(l1:l2,m,n,ilnrho) = rho0*exp(-0.5*(z(n)/scaleH)**2)
            else
              f(l1:l2,m,n,ilnrho) = alog(rho0)-0.5*(z(n)/scaleH)**2
            endif
            f(l1:l2,m,n,iax) = 2*scaleH*sqrt(2*mu0*pi*rho0*cs20/beta_const) &
                              *sin(kx_aa(j)*x(l1:l2))*erfunc(0.5*z(n)/scaleH)
          enddo; enddo

        case default
          !
          !  Catch unknown values
          !
          if (lroot) &
              print*, 'init_aa: No such value for initaa: ', trim(initaa(j))
          call stop_it("")

        endselect
        !
        !  End loop over initial conditions
        !
      enddo
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
              f(l1:l2,m,n,ilnrho)=f(l1:l2,m,n,ilnrho)+fact/gamma1
            endif
          endif
        enddo
        enddo
      endif
!
    endsubroutine init_aa
!***********************************************************************
    subroutine pert_aa(f)
!
!   perturb magnetic field when reading old NON-magnetic snapshot
!   called from run.f90; this uses a lot of memory and should be modified.
!
!   30-july-2004/dave: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      xx=spread(spread(x,2,my),3,mz)
      yy=spread(spread(y,1,mx),3,mz)
      zz=spread(spread(z,1,mx),2,my)
      initaa=pertaa
      amplaa=pertamplaa
      call init_aa(f,xx,yy,zz)
!
    endsubroutine pert_aa
!***********************************************************************
    subroutine pencil_criteria_magnetic()
!
!   All pencils that the Magnetic module depends on are specified here.
!
!  19-11-04/anders: coded
!
      use Cdata
!
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_uxb)=.true.
!
      if (dvid/=0.0) lpenc_video(i_b2)=.true.
!
!  jj pencil always needed when in Weyl gauge
!
      if ( (hall_term/=0. .and. ldt) .or. height_eta/=0. .or. ip<=4 .or. &
          (lweyl_gauge) .or. (.not.lcartesian_coords) ) &
          lpenc_requested(i_jj)=.true.
      if (eta/=0..and.(.not.lweyl_gauge)) lpenc_requested(i_del2a)=.true.
      if (dvid/=0.) lpenc_video(i_jb)=.true.
      if (lresi_eta_const .or. lresi_shell .or. &
          lresi_eta_shock .or. lresi_smagorinsky .or. &
          lresi_zdep .or. &
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
      if (lresi_shell .or. lresi_zdep) lpenc_requested(i_diva)=.true.
      if (lresi_smagorinsky_cross) lpenc_requested(i_jo)=.true.
      if (lresi_hyper2) lpenc_requested(i_del4a)=.true.
      if (lresi_hyper3) lpenc_requested(i_del6a)=.true.
      if (lspherical_coords) lpenc_requested(i_graddiva)=.true.
      if (lentropy .or. lresi_smagorinsky .or. ltemperature) then
        lpenc_requested(i_j2)=.true.
      endif

      if ((borderaa=='Alfven-rz').or.(borderaa=='Alfven-zconst')) then
        lpenc_requested(i_rcyl_mn) =.true.
        lpenc_requested(i_rcyl_mn1)=.true.
        lpenc_requested(i_phix)    =.true.
        lpenc_requested(i_phiy)    =.true.
      endif

      if (borderaa=='toroidal') &
           lpenc_requested(i_rcyl_mn) =.true.

      if (lentropy .or. ltemperature .or. ldt) lpenc_requested(i_rho1)=.true.
      if (lentropy .or. ltemperature) lpenc_requested(i_TT1)=.true.
      if (ltemperature) lpenc_requested(i_cv1)=.true.
      if (nu_ni/=0.) lpenc_requested(i_va2)=.true.
      if (hall_term/=0.) lpenc_requested(i_jxb)=.true.
      if ((lhydro .and. llorentzforce) .or. nu_ni/=0.) &
          lpenc_requested(i_jxbr)=.true.
      if (lresi_smagorinsky_cross .or. delta_effect/=0.) &
          lpenc_requested(i_oo)=.true.
      if (nu_ni/=0.) lpenc_requested(i_va2)=.true.
      if (lmeanfield_theory) then
        if (alpha_effect/=0. .or. delta_effect/=0.) lpenc_requested(i_mf_EMF)=.true.
        if (delta_effect/=0.) lpenc_requested(i_oxj)=.true.
      endif
      if (nu_ni/=0.) lpenc_diagnos(i_jxbrxb)=.true.
!
      if (     idiag_brmphi/=0  .or. idiag_uxbrmphi/=0 .or. idiag_jxbrmphi/=0 &
          .or. idiag_armphi/=0  .or. idiag_brmr/=0     .or. idiag_armr/=0 ) then
        lpenc_diagnos(i_pomx)=.true.
        lpenc_diagnos(i_pomy)=.true.
      endif
!
      if (     idiag_bpmphi/=0  .or. idiag_uxbpmphi/=0 .or. idiag_jxbpmphi/=0 &
          .or. idiag_bpmr/=0    .or. idiag_brbpmr/=0   .or. idiag_apmphi/=0  &
          .or. idiag_apmr/=0 ) then
        lpenc_diagnos(i_phix)=.true.
        lpenc_diagnos(i_phiy)=.true.
      endif
!
      if (idiag_armr/=0 .or. idiag_apmr/=0 .or. idiag_azmr/=0) &
           lpenc_diagnos(i_aa)=.true.
!
      if (idiag_armphi/=0 .or. idiag_apmphi/=0 .or. idiag_azmphi/=0) &
           lpenc_diagnos2d(i_aa)=.true.
!
      if (idiag_aybym2/=0 .or. idiag_exaym2/=0) lpenc_diagnos(i_aa)=.true.
      if (idiag_arms/=0 .or. idiag_amax/=0) lpenc_diagnos(i_a2)=.true.
      if (idiag_abm/=0) lpenc_diagnos(i_ab)=.true.
      if (idiag_djuidjbim/=0 .or. idiag_b3b21m/=0 .or. &
          idiag_b1b32m/=0 .or.  idiag_b2b13m/=0) &
          lpenc_diagnos(i_bij)=.true.
      if (idiag_j2m/=0 .or. idiag_jm2/=0 .or. idiag_jrms/=0 .or. &
          idiag_jmax/=0 .or. idiag_epsM/=0 .or. idiag_epsM_LES/=0) &
          lpenc_diagnos(i_j2)=.true.
      if (idiag_jbm/=0) lpenc_diagnos(i_jb)=.true.
      if (idiag_jbmphi/=0) lpenc_diagnos2d(i_jb)=.true.
      if (idiag_vArms/=0 .or. idiag_vAmax/=0 .or. idiag_vA2m/=0) lpenc_diagnos(i_va2)=.true.
      if (idiag_ubm/=0) lpenc_diagnos(i_ub)=.true.
      if (idiag_djuidjbim/=0 .or. idiag_uxDxuxbm/=0) lpenc_diagnos(i_uij)=.true.
      if (idiag_uxjm/=0) lpenc_diagnos(i_uxj)=.true.
      if (idiag_vArms/=0 .or. idiag_vAmax/=0) lpenc_diagnos(i_va2)=.true.
      if (idiag_uxBrms/=0 .or. idiag_Rmrms/=0) lpenc_diagnos(i_uxb2)=.true.
      if (idiag_beta1m/=0 .or. idiag_beta1max/=0) lpenc_diagnos(i_beta)=.true.
      if (idiag_djuidjbim/=0) lpenc_diagnos(i_djuidjbi)=.true.
      if (idiag_ujxbm/=0) lpenc_diagnos(i_ujxb)=.true.
      if (idiag_gpxbm/=0) lpenc_diagnos(i_glnrhoxb)=.true.
      if (idiag_jxbxbm/=0) lpenc_diagnos(i_jxbxb)=.true.
      if (idiag_oxuxbm/=0) lpenc_diagnos(i_oxuxb)=.true.
      if (idiag_exaym2/=0 .or. idiag_exjm2/=0) lpenc_diagnos(i_jj)=.true.
      if (idiag_b2m/=0 .or. idiag_bm2/=0 .or. idiag_brms/=0 .or. &
          idiag_bmax/=0) lpenc_diagnos(i_b2)=.true.
!
!  pencils for meanfield dynamo diagnostics
!
      if (idiag_EMFdotBm/=0) lpenc_diagnos(i_mf_EMFdotB)=.true.
!
      if (lisotropic_advection) lpenc_requested(i_va2)=.true.
!
    endsubroutine pencil_criteria_magnetic
!***********************************************************************
    subroutine pencil_interdep_magnetic(lpencil_in)
!
!  Interdependency among pencils from the Magnetic module is specified here.
!
!  19-11-04/anders: coded
!
      use Cdata
!
      logical, dimension(npencils) :: lpencil_in
!
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
      if (lpencil_in(i_ub)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_bb)=.true.
      endif
      if (lpencil_in(i_beta)) then
        lpencil_in(i_b2)=.true.
        lpencil_in(i_pp)=.true.
      endif
      if (lpencil_in(i_b2)) lpencil_in(i_bb)=.true.
      if (lpencil_in(i_jj)) lpencil_in(i_bij)=.true.
      if (lpencil_in(i_bb)) then
        if (lspherical_coords) lpencil_in(i_aa)=.true.
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
        lpencil_in(i_b2)=.true.
        lpencil_in(i_bb)=.true.
        if (delta_effect/=0.) lpencil_in(i_oxJ)=.true.
        if (meanfield_etat/=0.) lpencil_in(i_jj)=.true.
      endif
      if (lpencil_in(i_del2A)) then
        if (.not.lcartesian_coords) then
          lpencil_in(i_jj)=.true.
          lpencil_in(i_graddivA)=.true.
        endif
      endif
      if (lpencil_in(i_uga)) then
        lpencil_in(i_aij)=.true.
        lpencil_in(i_uu)=.true.
      endif
!XX
!     if (lspherical_coords) lpencil_in(i_aa)=.true.
!XX
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
      use Cdata
      use Sub
      use Deriv
      use Global, only: get_global
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: bb_ext,bb_ext_pot,ee_ext,jj_ext
      real, dimension (nx) :: rho1_jxb,alpha_total
      real, dimension (nx) :: alpha_tmp
      real :: B2_ext,c,s
      integer :: i,j
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
          if (omega_Bz_ext==0.) then
            B_ext_tmp=B_ext
          elseif (omega_Bz_ext/=0.) then
            c=cos(omega_Bz_ext*t)
            s=sin(omega_Bz_ext*t)
            B_ext_tmp(1)=B_ext(1)*c-B_ext(2)*s
            B_ext_tmp(2)=B_ext(1)*s+B_ext(2)*c
            B_ext_tmp(3)=B_ext(3)
          endif
!  add the external field
          if (B_ext_tmp(1)/=0.) p%bb(:,1)=p%bb(:,1)+B_ext_tmp(1)
          if (B_ext_tmp(2)/=0.) p%bb(:,2)=p%bb(:,2)+B_ext_tmp(2)
          if (B_ext_tmp(3)/=0.) p%bb(:,3)=p%bb(:,3)+B_ext_tmp(3)
          if (headtt) print*,'calc_pencils_magnetic: B_ext=',B_ext
          if (headtt) print*,'calc_pencils_magnetic: B_ext_tmp=',B_ext_tmp
        endif
!  add the external potential field
        if (lB_ext_pot) then
          call get_global(bb_ext_pot,m,n,'B_ext_pot')
          p%bb=p%bb+bb_ext_pot
        endif
!  add external B-field.
        if (iglobal_bx_ext/=0) p%bb(:,1)=p%bb(:,1)+f(l1:l2,m,n,iglobal_bx_ext)
        if (iglobal_by_ext/=0) p%bb(:,2)=p%bb(:,2)+f(l1:l2,m,n,iglobal_by_ext)
        if (iglobal_bz_ext/=0) p%bb(:,3)=p%bb(:,3)+f(l1:l2,m,n,iglobal_bz_ext)
      endif
! ab
      if (lpencil(i_ab)) call dot_mn(p%aa,p%bbb,p%ab)
! uxb
      if (lpencil(i_uxb)) then
        call cross_mn(p%uu,p%bb,p%uxb)
        call cross_mn(p%uu,p%bbb,uxbb)
!  add external e-field.
        if (iglobal_ex_ext/=0) p%uxb(:,1)=p%uxb(:,1)+f(l1:l2,m,n,iglobal_ex_ext)
        if (iglobal_ey_ext/=0) p%uxb(:,2)=p%uxb(:,2)+f(l1:l2,m,n,iglobal_ey_ext)
        if (iglobal_ez_ext/=0) p%uxb(:,3)=p%uxb(:,3)+f(l1:l2,m,n,iglobal_ez_ext)
      endif
! uga
! DM : this requires later attention
      if (lpencil(i_uga)) then
        if(lspherical_coords) then
          call warning("calc_pencils_magnetic","u_dot_grad A not implemented for spherical coordinates") 
        else
          call u_dot_grad(f,iaa,p%aij,p%uu,p%uga,UPWIND=lupw_aa)
        endif
      endif
! b2
      if (lpencil(i_b2)) call dot2_mn(p%bb,p%b2)
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
!  add external j-field.
        if (iglobal_jx_ext/=0) p%jj(:,1)=p%jj(:,1)+f(l1:l2,m,n,iglobal_jx_ext)
        if (iglobal_jy_ext/=0) p%jj(:,2)=p%jj(:,2)+f(l1:l2,m,n,iglobal_jy_ext)
        if (iglobal_jz_ext/=0) p%jj(:,3)=p%jj(:,3)+f(l1:l2,m,n,iglobal_jz_ext)
      endif
! j2
      if (lpencil(i_j2)) call dot2_mn(p%jj,p%j2)
! jb
      if (lpencil(i_jb)) call dot_mn(p%jj,p%bbb,p%jb)
! va2
      if (lpencil(i_va2)) p%va2=p%b2*mu01*p%rho1
! jxb
      if (lpencil(i_jxb)) call cross_mn(p%jj,p%bb,p%jxb)
! jxbr
      if (lpencil(i_jxbr)) then
        rho1_jxb=p%rho1
!  set rhomin_jxb>0 in order to limit the jxb term at very low densities.
!  set va2max_jxb>0 in order to limit the jxb term at very high Alfven speeds.
!  set va2power_jxb to an integer value in order to specify the power
!  of the limiting term,
        if (rhomin_jxb>0) rho1_jxb=min(rho1_jxb,1/rhomin_jxb)
        if (va2max_jxb>0) then
          rho1_jxb = rho1_jxb &
                   * (1+(p%va2/va2max_jxb)**va2power_jxb)**(-1.0/va2power_jxb)
        endif
        call multsv_mn(rho1_jxb,p%jxb,p%jxbr)
      endif
! ub
      if (lpencil(i_ub)) call dot_mn(p%uu,p%bb,p%ub)
! uxb2
      if (lpencil(i_uxb2)) call dot2_mn(p%uxb,p%uxb2)
! uxj
      if (lpencil(i_uxj)) call cross_mn(p%uu,p%jj,p%uxj)
! beta
      if (lpencil(i_beta)) p%beta=0.5*p%b2/p%pp
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
        select case(alpha_profile)
        case('nothing'); alpha_tmp=1.
        case('siny'); alpha_tmp=sin(y(m))
        case('cosy'); alpha_tmp=cos(y(m))
        case('read'); alpha_tmp=alpha_input(l1:l2,m)
        endselect
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
        if (meanfield_etat/=0.) p%mf_EMF=p%mf_EMF-meanfield_etat*p%jj
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
!  22-nov-01/nils: coded
!   1-may-02/wolf: adapted for pencil_modular
!  17-jun-03/ulf:  added bx^2, by^2 and bz^2 as separate diagnostics
!   8-aug-03/axel: introduced B_ext21=1./B_ext**2, and set =1 to avoid div. by 0
!  12-aug-03/christer: added alpha effect (alpha in the equation above)
!  26-may-04/axel: ambipolar diffusion added
!  18-jun-04/axel: Hall term added
!
      use Cdata
      use Deriv, only: der6
      use EquationOfState, only: eoscalc,gamma1
      use Io, only: output_pencil
      use Mpicomm, only: stop_it
      use Special, only: special_calc_magnetic
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: geta,uxDxuxb,fres,uxb_upw,tmp2
      real, dimension (nx) :: uxb_dotB0,oxuxb_dotB0,jxbxb_dotB0,uxDxuxb_dotB0
      real, dimension (nx) :: gpxb_dotB0,uxj_dotB0,b3b21,b1b32,b2b13,sign_jo,rho1_jxb
      real, dimension (nx) :: B1dot_glnrhoxb
      real, dimension (nx) :: eta_mn,eta_smag,etatotal,fres2,etaSS,penc
      real :: tmp,eta_out1,OmegaSS=1.
      integer :: i,j,k
!
      intent(in)     :: f,p
      intent(inout)  :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'daa_dt: SOLVE'
      if (headtt) then
        call identify_bcs('Ax',iax)
        call identify_bcs('Ay',iay)
        call identify_bcs('Az',iaz)
      endif
!
!  add jxb/rho to momentum equation
!
      if (lhydro) then
        if (llorentzforce) df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+p%jxbr
      endif
!
!  Restivivity term
!
!  Because of gauge invariance, we can add the gradient of an arbitrary scalar
!  field Phi to the induction equation without changing the magnetic field,
!    dA/dt = u x B - eta j + grad(Phi).
!
!  If lweyl_gauge=T, we choose Phi = const. and solve
!    dA/dt = u x B - eta j.
!  Else, if lweyl_gauge=F, we make the gauge choice Phi = eta div(A)
!  and thus solve
!    dA/dt = u x B + eta laplace(A) + div(A) grad(eta).
!
!  Note: lweyl_gauge=T is so far only implemented for shock resistivity.
!
      fres=0.0
      etatotal=0.0
!
      if (lresi_eta_const) then
        if (lweyl_gauge) then
          fres = fres - eta*p%jj
        else
          fres = fres + eta*p%del2a
        endif
        etatotal=etatotal+eta
      endif
!
!  Shakura-Sunyaev type resistivity (mainly just as a demo to show
!  how resistivity can be made depend on temperature.
!  Since etaSS is nonuniform, we use this contribution only for -etaSS*JJ
!  and keep the constant piece with +eta*del2A. (The divA term is eliminated
!  by a suitable gauge transformation.) A sample run is checked in under
!  pencil-runs/1d-tests/bdecay
!
      if (lresi_etaSS) then
        etaSS=alphaSSm*p%cs2/OmegaSS
        do j=1,3
          fres(:,j)=fres(:,j)+eta*p%del2a(:,j)-etaSS*p%jj(:,j)
        enddo
        etatotal=etaSS+eta
      endif
!
      if (lresi_hyper2) then
        fres=fres+eta_hyper2*p%del4a
      endif
!
      if (lresi_hyper3) then
        fres=fres+eta_hyper3*p%del6a
      endif
!
      if (lresi_hyper3_strict) then
        fres=fres+eta_hyper3*f(l1:l2,m,n,ihypres:ihypres+2)
      endif
!
      if (lresi_zdep) then 
        do j=1,3
          fres(:,j)=fres(:,j)+eta_z(n)*p%del2a(:,j)+geta_z(n,j)*p%diva
        enddo
        etatotal=etatotal + eta_z(n)
      endif
!
      if (lresi_hyper3_aniso) then
         call del6fjv(f,eta_aniso_hyper3,iaa,tmp2)
         fres=fres+tmp2
      endif
!
      if (lresi_shell) then
        call eta_shell(p,eta_mn,geta)
        do j=1,3
          fres(:,j)=fres(:,j)+eta_mn*p%del2a(:,j)+geta(:,j)*p%diva
        enddo
        etatotal=etatotal+eta_mn
      endif
!
      if (lresi_eta_shock) then
        if (lweyl_gauge) then
          do i=1,3
            fres(:,i) = fres(:,i) - eta_shock*p%shock*p%jj(:,i)
          enddo
        else
          do i=1,3
            fres(:,i) = fres(:,i) &
                      + eta_shock*(p%shock*p%del2a(:,i)+p%diva*p%gshock(:,i))
          enddo
        endif
        etatotal=etatotal+eta_shock*p%shock
      endif
!
      if (lresi_eta_shock_perp) then
        if (lweyl_gauge) then
          do i=1,3
            fres(:,i) = fres(:,i) - eta_shock*p%shock_perp*p%jj(:,i)
          enddo
        else
          do i=1,3
            fres(:,i) = fres(:,i) &
                      + eta_shock*(p%shock_perp*p%del2a(:,i) &
                                  +p%diva*p%gshock_perp(:,i))
          enddo
        endif
        etatotal=etatotal+eta_shock*p%shock_perp
      endif
!
      if (lresi_smagorinsky) then
        eta_smag=(D_smag*dxmax)**2.*sqrt(p%j2)
        call multsv(eta_smag+eta,p%del2a,fres)
        etatotal=etatotal+eta_smag+eta
      endif
!
      if (lresi_smagorinsky_cross) then
        sign_jo=1.
        do i=1,nx
          if (p%jo(i) .lt. 0) sign_jo(i)=-1.
        enddo
        eta_smag=(D_smag*dxmax)**2.*sign_jo*sqrt(p%jo*sign_jo)
        call multsv(eta_smag+eta,p%del2a,fres)
        etatotal=eta_smag+eta
      endif
!
      if (headtt) print*,'daa_dt: iresistivity=',iresistivity
!
!  add eta mu_0 j2/rho to entropy or temperature equation
!
      if (lentropy .and. lohmic_heat) then
        df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) &
                          + etatotal*mu0*p%j2*p%rho1*p%TT1
      endif

      if (ltemperature .and. lohmic_heat) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) &
                            + etatotal*mu0*p%j2*p%rho1*p%cv1*p%TT1
      endif
!
!  Switch off diffusion in boundary slice if requested by boundconds
!
!  Only need to do this on bottommost (topmost) processors
!  and in bottommost (topmost) pencils
!
      do j=1,3
        if (lfrozen_bb_bot(j).and.ipz==0       .and.n==n1) fres(:,j)=0.
        if (lfrozen_bb_top(j).and.ipz==nprocz-1.and.n==n2) fres(:,j)=0.
      enddo
!
!  Induction equation
!
      if (.not.lupw_aa) then
        df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + p%uxb + fres
      else
!
!  Use upwinding for the advection term.
!
!  We only do upwinding for advection-like terms, u_i f_k,j,
!  for which i=j. This means that for instance in the evolution
!  equation for A_x,
!
!  d(A_x)/dt + u_y A_x,y + u_z A_x,z = u_y A_y,x + u_z A_z,x
!
!  we only do upwinding for the advection-type terms on the
!  left hand side.
!
        if (lupw_aa.and.headtt) then
          print *,'calc_pencils_magnetic: upwinding advection term. '//&
                  'Not well tested; use at own risk!'; endif
!  Add Lorentz force that results from the external field.
!  Note: For now, this only works for uniform external fields.
        uxb_upw(:,1) = p%uu(:,2)*B_ext(3) - p%uu(:,3)*B_ext(2)
        uxb_upw(:,2) = p%uu(:,3)*B_ext(1) - p%uu(:,1)*B_ext(3)
        uxb_upw(:,3) = p%uu(:,1)*B_ext(2) - p%uu(:,2)*B_ext(1)
!  Add u_k A_k,j and `upwinded' advection term
!  Note: this works currently only in cartesian geometry!
        do j=1,3; do k=1,3
          if (k/=j) then
            uxb_upw(:,j) = uxb_upw(:,j) + p%uu(:,k)*(p%aij(:,k,j)-p%aij(:,j,k))
            call der6(f,iaa+j-1,penc,k,upwind=.true.)
            uxb_upw(:,j) = uxb_upw(:,j) + abs(p%uu(:,k))*penc
          endif
        enddo; enddo
!  Full right hand side of the induction equation
        df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + uxb_upw + fres
      endif
!
!  Ambipolar diffusion in the strong coupling approximation
!
      if (nu_ni/=0.) then
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+nu_ni1*p%jxbrxb
        etatotal=etatotal+nu_ni1*p%va2
      endif
!
!  Hall term
!
      if (hall_term/=0.) then
        if (headtt) print*,'daa_dt: hall_term=',hall_term
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-hall_term*p%jxb
        if (lfirst.and.ldt) then
          advec_hall=abs(p%uu(:,1)-hall_term*p%jj(:,1))*dx_1(l1:l2)+ &
                     abs(p%uu(:,2)-hall_term*p%jj(:,2))*dy_1(  m  )+ &
                     abs(p%uu(:,3)-hall_term*p%jj(:,3))*dz_1(  n  )
        endif
        if (headtt.or.ldebug) print*,'daa_dt: max(advec_hall) =',&
                                     maxval(advec_hall)
      endif
!
!  Alpha effect
!  additional terms if Mean Field Theory is included
!
      if (lmeanfield_theory.and.(alpha_effect/=0..or.delta_effect/=0.)) then
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+p%mf_EMF
        if (lOmega_effect) call Omega_effect(f,df)
      endif
!
!  Possibility of adding extra diffusivity in some halo of given geometry:
!  Note that eta_out is total eta in halo (not eta_out+eta)
!
      if (height_eta/=0.) then
        if (headtt) print*,'daa_dt: height_eta,eta_out=',height_eta,eta_out
        tmp=(z(n)/height_eta)**2
        eta_out1=eta_out*(1.-exp(-tmp**5/max(1.-tmp,1e-5)))-eta
        df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)-(eta_out1*mu0)*p%jj
      endif
!
!  add possibility of forcing that is not delta-correlated in time
!
      if (lforcing_continuous_aa) call forcing_continuous(df,p)
!
!  possibility of relaxation of A in exterior region
!
      if (tau_aa_exterior/=0.) call calc_tau_aa_exterior(f,df)
!
!  ``va^2/dx^2'' and ``eta/dx^2'' for timestep
!  in the diffusive timestep, we include possible contribution from
!  meanfield_etat, which is however only invoked in mean field models
!  Consider advective timestep only when lhydro=T.
!
      if (lfirst.and.ldt) then
        if (lhydro) then
          rho1_jxb=p%rho1
          if (rhomin_jxb>0) rho1_jxb=min(rho1_jxb,1/rhomin_jxb)
          if (va2max_jxb>0) then
            rho1_jxb = rho1_jxb &
                     * (1+(p%va2/va2max_jxb)**va2power_jxb)**(-1.0/va2power_jxb)
          endif
          if (lspherical_coords) then 
            advec_va2=((p%bb(:,1)*dx_1(l1:l2))**2+ &
                       (p%bb(:,2)*dy_1(  m  )*r1_mn)**2+ &
                       (p%bb(:,3)*dz_1(  n  )*r1_mn*sin1th(m))**2)*mu01*rho1_jxb
          elseif (lcylindrical_coords) then 
            advec_va2=((p%bb(:,1)*dx_1(l1:l2))**2+ &
                       (p%bb(:,2)*dy_1(  m  )*rcyl_mn1)**2+ &
                       (p%bb(:,3)*dz_1(  n  ))**2)*mu01*rho1_jxb
          else
            advec_va2=((p%bb(:,1)*dx_1(l1:l2))**2+ &
                       (p%bb(:,2)*dy_1(  m  ))**2+ &
                       (p%bb(:,3)*dz_1(  n  ))**2)*mu01*rho1_jxb
        endif
      endif
!
!WL: don't know if this is correct, but it's the only way I can make
!    some 1D and 2D samples work when the non-existent direction has the
!    largest velocity (like a 2D rz slice of a Keplerian disk that rotates
!    on the phi direction)
!    Please check
!
        if (lisotropic_advection) then
          if (lfirst.and.ldt) then
            if ((nxgrid==1).or.(nygrid==1).or.(nzgrid==1)) &
                 advec_va2=sqrt(p%va2*dxyz_2)
          endif
        endif

!
!  resistive time step considerations
!
        if (lresi_hyper3_aniso) then
           diffus_eta=eta_aniso_hyper3(1)*dx_1(l1:l2)**6 + &
                eta_aniso_hyper3(2)*dy_1(m)**6 + &
                eta_aniso_hyper3(3)*dz_1(n)**6 + &
                etatotal*dxyz_2
        elseif (lresi_hyper3) then
           diffus_eta=eta_hyper3*dxyz_6 + (etatotal+meanfield_etat)*dxyz_2
        elseif (lresi_hyper2) then
           diffus_eta=eta_hyper2*dxyz_4 + (etatotal+meanfield_etat)*dxyz_2
        else
           diffus_eta=(etatotal+meanfield_etat)*dxyz_2
        endif
        if (ldiagnos.and.idiag_dteta/=0) then
          call max_mn_name(diffus_eta/cdtv,idiag_dteta,l_dt=.true.)
        endif
      endif
      if (headtt.or.ldebug) then
        print*,'daa_dt: max(advec_va2) =',maxval(advec_va2)
        print*,'daa_dt: max(diffus_eta) =',maxval(diffus_eta)
      endif
!
!  Special contributions to this module are called here
!
      if (lspecial) call special_calc_magnetic(f,df,p)
!
!  Apply border profiles
!
      if (lborder_profiles) call set_border_magnetic(f,df,p)
!
!  Calculate diagnostic quantities
!
      if (ldiagnos) then
 
        if (idiag_beta1m/=0) call sum_mn_name(p%beta,idiag_beta1m)
        if (idiag_beta1max/=0) call max_mn_name(p%beta,idiag_beta1max)

        if (idiag_b2m/=0) call sum_mn_name(p%b2,idiag_b2m)
        if (idiag_bm2/=0) call max_mn_name(p%b2,idiag_bm2)
        if (idiag_brms/=0) call sum_mn_name(p%b2,idiag_brms,lsqrt=.true.)
        if (idiag_bmax/=0) call max_mn_name(p%b2,idiag_bmax,lsqrt=.true.)
        if (idiag_aybym2/=0) &
            call sum_mn_name(2*p%aa(:,2)*p%bb(:,2),idiag_aybym2)
        if (idiag_abm/=0) call sum_mn_name(p%ab,idiag_abm)
        if (idiag_ubm/=0) call sum_mn_name(p%ub,idiag_ubm)
        if (idiag_bxm/=0) call sum_mn_name(p%bbb(:,1),idiag_bxm)
        if (idiag_bym/=0) call sum_mn_name(p%bbb(:,2),idiag_bym)
        if (idiag_bzm/=0) call sum_mn_name(p%bbb(:,3),idiag_bzm)
        if (idiag_bx2m/=0) call sum_mn_name(p%bbb(:,1)**2,idiag_bx2m)
        if (idiag_by2m/=0) call sum_mn_name(p%bbb(:,2)**2,idiag_by2m)
        if (idiag_bz2m/=0) call sum_mn_name(p%bbb(:,3)**2,idiag_bz2m)
        if (idiag_bxbym/=0) call sum_mn_name(p%bbb(:,1)*p%bbb(:,2),idiag_bxbym)
        if (idiag_bxbzm/=0) call sum_mn_name(p%bbb(:,1)*p%bbb(:,3),idiag_bxbzm)
        if (idiag_bybzm/=0) call sum_mn_name(p%bbb(:,2)*p%bbb(:,3),idiag_bybzm)

        if (idiag_djuidjbim/=0) call sum_mn_name(p%djuidjbi,idiag_djuidjbim)
!
!  magnetic field components at one point (=pt)
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_bxpt/=0) call save_name(p%bb(lpoint-nghost,1),idiag_bxpt)
          if (idiag_bypt/=0) call save_name(p%bb(lpoint-nghost,2),idiag_bypt)
          if (idiag_bzpt/=0) call save_name(p%bb(lpoint-nghost,3),idiag_bzpt)
        endif
!
!  v_A = |B|/sqrt(rho); in units where mu_0=1
!
        if (idiag_vA2m/=0)   call sum_mn_name(p%va2,idiag_vA2m)
        if (idiag_vArms/=0) call sum_mn_name(p%va2,idiag_vArms,lsqrt=.true.)
        if (idiag_vAmax/=0) call max_mn_name(p%va2,idiag_vAmax,lsqrt=.true.)
        if (idiag_dtb/=0) &
            call max_mn_name(sqrt(advec_va2)/cdt,idiag_dtb,l_dt=.true.)
!
! <J.B>
!
        if (idiag_jbm/=0) call sum_mn_name(p%jb,idiag_jbm)
        if (idiag_j2m/=0) call sum_mn_name(p%j2,idiag_j2m)
        if (idiag_jm2/=0) call max_mn_name(p%j2,idiag_jm2)
        if (idiag_jrms/=0) call sum_mn_name(p%j2,idiag_jrms,lsqrt=.true.)
        if (idiag_jmax/=0) call max_mn_name(p%j2,idiag_jmax,lsqrt=.true.)
        if (idiag_epsM_LES/=0) call sum_mn_name(eta_smag*p%j2,idiag_epsM_LES)
!
!  Not correct for hyperresistivity:
!
        if (idiag_epsM/=0) call sum_mn_name(eta*p%j2,idiag_epsM)
!
! <A^2> and A^2|max
!
        if (idiag_arms/=0) call sum_mn_name(p%a2,idiag_arms,lsqrt=.true.)
        if (idiag_amax/=0) call max_mn_name(p%a2,idiag_amax,lsqrt=.true.)
!
!  calculate surface integral <2ExA>*dS
!
        if (idiag_exaym2/=0) call helflux(p%aa,p%uxb,p%jj)
!
!  calculate surface integral <2ExJ>*dS
!
        if (idiag_exjm2/=0) call curflux(p%uxb,p%jj)
!
!  calculate emf for alpha effect (for imposed field)
!  Note that uxbm means <EMF.B0>/B0^2, so it gives already alpha=EMF/B0.
!
        if (idiag_uxbm/=0 .or. idiag_uxbmx/=0 .or. idiag_uxbmy/=0 &
            .or. idiag_uxbmz/=0) then
          uxb_dotB0=B_ext(1)*p%uxb(:,1)+B_ext(2)*p%uxb(:,2)+B_ext(3)*p%uxb(:,3)
          uxb_dotB0=uxb_dotB0*B_ext21
          if (idiag_uxbm/=0) call sum_mn_name(uxb_dotB0,idiag_uxbm)
          if (idiag_uxbmx/=0) call sum_mn_name(uxbb(:,1),idiag_uxbmx)
          if (idiag_uxbmy/=0) call sum_mn_name(uxbb(:,2),idiag_uxbmy)
          if (idiag_uxbmz/=0) call sum_mn_name(uxbb(:,3),idiag_uxbmz)
        endif
!
!  calculate <uxj>.B0/B0^2
!
        if (idiag_uxjm/=0) then
          uxj_dotB0=B_ext(1)*p%uxj(:,1)+B_ext(2)*p%uxj(:,2)+B_ext(3)*p%uxj(:,3)
          uxj_dotB0=uxj_dotB0*B_ext21
          call sum_mn_name(uxj_dotB0,idiag_uxjm)
        endif
!
!  calculate <u x B>_rms, <resistive terms>_rms, <ratio ~ Rm>_rms
!
        if (idiag_uxBrms/=0) call sum_mn_name(p%uxb2,idiag_uxBrms,lsqrt=.true.)
        if (idiag_Bresrms/=0 .or. idiag_Rmrms/=0) then
          call dot2_mn(fres,fres2)
          if (idiag_Bresrms/=0) &
              call sum_mn_name(fres2,idiag_Bresrms,lsqrt=.true.)
          if (idiag_Rmrms/=0) &
              call sum_mn_name(p%uxb2/fres2,idiag_Rmrms,lsqrt=.true.)
        endif
!
!  calculate <u.(jxb)>
!
        if (idiag_ujxbm/=0) call sum_mn_name(p%ujxb,idiag_ujxbm)
!
!  magnetic triple correlation term (for imposed field)
!
        if (idiag_jxbxbm/=0) then
          jxbxb_dotB0=B_ext(1)*p%jxbxb(:,1)+B_ext(2)*p%jxbxb(:,2)+B_ext(3)*p%jxbxb(:,3)
          jxbxb_dotB0=jxbxb_dotB0*B_ext21
          call sum_mn_name(jxbxb_dotB0,idiag_jxbxbm)
        endif
!
!  triple correlation from Reynolds tensor (for imposed field)
!
        if (idiag_oxuxbm/=0) then
          oxuxb_dotB0=B_ext(1)*p%oxuxb(:,1)+B_ext(2)*p%oxuxb(:,2)+B_ext(3)*p%oxuxb(:,3)
          oxuxb_dotB0=oxuxb_dotB0*B_ext21
          call sum_mn_name(oxuxb_dotB0,idiag_oxuxbm)
        endif
!
!  triple correlation from pressure gradient (for imposed field)
!  (assume cs2=1, and that no entropy evolution is included)
!  This is ok for all applications currently under consideration.
!
        if (idiag_gpxbm/=0) then
          call dot_mn_sv(B1_ext,p%glnrhoxb,B1dot_glnrhoxb)
          call sum_mn_name(B1dot_glnrhoxb,idiag_gpxbm)
        endif
!
!  < u x curl(uxB) > = < E_i u_{j,j} - E_j u_{j,i} >
!   ( < E_1 u2,2 + E1 u3,3 - E2 u2,1 - E3 u3,1 >
!     < E_2 u1,1 + E2 u3,3 - E1 u2,1 - E3 u3,2 >
!     < E_3 u1,1 + E3 u2,2 - E1 u3,1 - E2 u2,3 > )
!
        if (idiag_uxDxuxbm/=0) then
          uxDxuxb(:,1)=p%uxb(:,1)*(p%uij(:,2,2)+p%uij(:,3,3))-p%uxb(:,2)*p%uij(:,2,1)-p%uxb(:,3)*p%uij(:,3,1)
          uxDxuxb(:,2)=p%uxb(:,2)*(p%uij(:,1,1)+p%uij(:,3,3))-p%uxb(:,1)*p%uij(:,1,2)-p%uxb(:,3)*p%uij(:,3,2)
          uxDxuxb(:,3)=p%uxb(:,3)*(p%uij(:,1,1)+p%uij(:,2,2))-p%uxb(:,1)*p%uij(:,1,3)-p%uxb(:,2)*p%uij(:,2,3)
          uxDxuxb_dotB0=B_ext(1)*uxDxuxb(:,1)+B_ext(2)*uxDxuxb(:,2)+B_ext(3)*uxDxuxb(:,3)
          uxDxuxb_dotB0=uxDxuxb_dotB0*B_ext21
          call sum_mn_name(uxDxuxb_dotB0,idiag_uxDxuxbm)
        endif
!
!  alpM11=<b3*b2,1>
!
        if (idiag_b3b21m/=0) then
          b3b21=p%bb(:,3)*p%bij(:,2,1)
          call sum_mn_name(b3b21,idiag_b3b21m)
        endif
!
!  alpM22=<b1*b3,2>
!
        if (idiag_b1b32m/=0) then
          b1b32=p%bb(:,1)*p%bij(:,3,2)
          call sum_mn_name(b1b32,idiag_b1b32m)
        endif
!
!  alpM33=<b2*b1,3>
!
        if (idiag_b2b13m/=0) then
          b2b13=p%bb(:,2)*p%bij(:,1,3)
          call sum_mn_name(b2b13,idiag_b2b13m)
        endif
!
!  diagnostic output for mean field dynamos
!
        if (idiag_EMFdotBm/=0) call sum_mn_name(p%mf_EMFdotB,idiag_EMFdotBm)
      endif ! endif (ldiagnos)
!                                                                               
!  1d-averages. Happens at every it1d timesteps, NOT at every it1               
!             
      if (l1ddiagnos .or. (ldiagnos .and. ldiagnos_need_zaverages)) then
        if (idiag_bxmz/=0)   call xysum_mn_name_z(p%bb(:,1),idiag_bxmz)
        if (idiag_bymz/=0)   call xysum_mn_name_z(p%bb(:,2),idiag_bymz)
        if (idiag_bzmz/=0)   call xysum_mn_name_z(p%bb(:,3),idiag_bzmz)
        if (idiag_bx2mz/=0)  call xysum_mn_name_z(p%bb(:,1)**2,idiag_bx2mz)
        if (idiag_by2mz/=0)  call xysum_mn_name_z(p%bb(:,2)**2,idiag_by2mz)
        if (idiag_bz2mz/=0)  call xysum_mn_name_z(p%bb(:,3)**2,idiag_bz2mz)
        if (idiag_bxbymz/=0) &
            call xysum_mn_name_z(p%bbb(:,1)*p%bbb(:,2),idiag_bxbymz)
        if (idiag_bxbzmz/=0) &
            call xysum_mn_name_z(p%bbb(:,1)*p%bbb(:,3),idiag_bxbzmz)
        if (idiag_bybzmz/=0) &
            call xysum_mn_name_z(p%bbb(:,2)*p%bbb(:,3),idiag_bybzmz)
        if (idiag_b2mz/=0)   call xysum_mn_name_z(p%b2,idiag_b2mz)
!  phi-z averages
        if (idiag_b2mr/=0)   call phizsum_mn_name_r(p%b2,idiag_b2mr)
        if (idiag_brmr/=0)   &
             call phizsum_mn_name_r(p%bb(:,1)*p%pomx+p%bb(:,2)*p%pomy,idiag_brmr)
        if (idiag_bpmr/=0)   &
             call phizsum_mn_name_r(p%bb(:,1)*p%phix+p%bb(:,2)*p%phiy,idiag_bpmr)
        if (idiag_bzmr/=0)   &
             call phizsum_mn_name_r(p%bb(:,3),idiag_bzmr)
        if (idiag_armr/=0)   &
             call phizsum_mn_name_r(p%aa(:,1)*p%pomx+p%aa(:,2)*p%pomy,idiag_armr)
        if (idiag_apmr/=0)   &
             call phizsum_mn_name_r(p%aa(:,1)*p%phix+p%aa(:,2)*p%phiy,idiag_apmr)
        if (idiag_azmr/=0)   &
             call phizsum_mn_name_r(p%aa(:,3),idiag_azmr)
        if (idiag_bxmx/=0) call yzsum_mn_name_x(p%bb(:,1),idiag_bxmx)
        if (idiag_bymy/=0) call xzsum_mn_name_y(p%bb(:,2),idiag_bymy)
        if (idiag_mflux_x/=0) &
          call yzintegrate_mn_name_x(p%bb(:,1),idiag_mflux_x)
        if (idiag_mflux_y/=0) &
          call xzintegrate_mn_name_y(p%bb(:,2),idiag_mflux_y)
        if (idiag_mflux_z/=0) &
          call xyintegrate_mn_name_z(p%bb(:,3),idiag_mflux_z)
      endif
!
!  phi-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        call phisum_mn_name_rz(p%bb(:,1)*p%pomx+p%bb(:,2)*p%pomy,idiag_brmphi)
        call phisum_mn_name_rz(p%bb(:,1)*p%phix+p%bb(:,2)*p%phiy,idiag_bpmphi)
        call phisum_mn_name_rz(p%bb(:,3),idiag_bzmphi)
        call phisum_mn_name_rz(p%b2,idiag_b2mphi)
        if (idiag_jbmphi/=0) call phisum_mn_name_rz(p%jb,idiag_jbmphi)
        if (any((/idiag_uxbrmphi,idiag_uxbpmphi,idiag_uxbzmphi/) /= 0)) then
          call phisum_mn_name_rz(p%uxb(:,1)*p%pomx+p%uxb(:,2)*p%pomy,idiag_uxbrmphi)
          call phisum_mn_name_rz(p%uxb(:,1)*p%phix+p%uxb(:,2)*p%phiy,idiag_uxbpmphi)
          call phisum_mn_name_rz(p%uxb(:,3)                     ,idiag_uxbzmphi)
        endif
        if (any((/idiag_jxbrmphi,idiag_jxbpmphi,idiag_jxbzmphi/) /= 0)) then
          call phisum_mn_name_rz(p%jxb(:,1)*p%pomx+p%jxb(:,2)*p%pomy,idiag_jxbrmphi)
          call phisum_mn_name_rz(p%jxb(:,1)*p%phix+p%jxb(:,2)*p%phiy,idiag_jxbpmphi)
          call phisum_mn_name_rz(p%jxb(:,3)                         ,idiag_jxbzmphi)
        endif
        if (any((/idiag_armphi,idiag_apmphi,idiag_azmphi/) /= 0)) then
          call phisum_mn_name_rz(p%aa(:,1)*p%pomx+p%aa(:,2)*p%pomy,idiag_armphi)
          call phisum_mn_name_rz(p%aa(:,1)*p%phix+p%aa(:,2)*p%phiy,idiag_apmphi)
          call phisum_mn_name_rz(p%aa(:,3)                        ,idiag_azmphi)
        endif
        if (idiag_bxmxy/=0) call zsum_mn_name_xy(p%bb(:,1),idiag_bxmxy)
        if (idiag_bymxy/=0) call zsum_mn_name_xy(p%bb(:,2),idiag_bymxy)
        if (idiag_bzmxy/=0) call zsum_mn_name_xy(p%bb(:,3),idiag_bzmxy)
        if (idiag_bxmxz/=0) call ysum_mn_name_xz(p%bb(:,1),idiag_bxmxz)
        if (idiag_bymxz/=0) call ysum_mn_name_xz(p%bb(:,2),idiag_bymxz)
        if (idiag_bzmxz/=0) call ysum_mn_name_xz(p%bb(:,3),idiag_bzmxz)
      else
!
!  idiag_bxmxy and idiag_bymxy also need to be calculated when
!  ldiagnos and idiag_bmx and/or idiag_bmy, so
!
        if (ldiagnos .and. (idiag_bmx/=0 .or. idiag_bmy/=0)) then
          if (idiag_bxmxy/=0) call zsum_mn_name_xy(p%bb(:,1),idiag_bxmxy)
          if (idiag_bymxy/=0) call zsum_mn_name_xy(p%bb(:,2),idiag_bymxy)
          if (idiag_bzmxy/=0) call zsum_mn_name_xy(p%bb(:,3),idiag_bzmxy)
        endif
      endif
!
!  debug output
!
      if (headtt .and. lfirst .and. ip<=4) then
        call output_pencil(trim(directory)//'/aa.dat',p%aa,3)
        call output_pencil(trim(directory)//'/bb.dat',p%bb,3)
        call output_pencil(trim(directory)//'/jj.dat',p%jj,3)
        call output_pencil(trim(directory)//'/del2A.dat',p%del2a,3)
        call output_pencil(trim(directory)//'/JxBr.dat',p%jxbr,3)
        call output_pencil(trim(directory)//'/JxB.dat',p%jxb,3)
        call output_pencil(trim(directory)//'/df.dat',df(l1:l2,m,n,:),mvar)
      endif
!
!  write B-slices for output in wvid in run.f90
!  Note: ix is the index with respect to array with ghost zones.
!
      if (lvid.and.lfirst) then
        do j=1,3
          bb_yz(m-m1+1,n-n1+1,j)=p%bb(ix_loc-l1+1,j)
          if (m==iy_loc)  bb_xz(:,n-n1+1,j)=p%bb(:,j)
          if (n==iz_loc)  bb_xy(:,m-m1+1,j)=p%bb(:,j)
          if (n==iz2_loc) bb_xy2(:,m-m1+1,j)=p%bb(:,j)
        enddo
        do j=1,3
          jj_yz(m-m1+1,n-n1+1,j)=p%jj(ix_loc-l1+1,j)
          if (m==iy_loc)  jj_xz(:,n-n1+1,j)=p%jj(:,j)
          if (n==iz_loc)  jj_xy(:,m-m1+1,j)=p%jj(:,j)
          if (n==iz2_loc) jj_xy2(:,m-m1+1,j)=p%jj(:,j)
        enddo
        b2_yz(m-m1+1,n-n1+1)=p%b2(ix_loc-l1+1)
        if (m==iy_loc)  b2_xz(:,n-n1+1)=p%b2
        if (n==iz_loc)  b2_xy(:,m-m1+1)=p%b2
        if (n==iz2_loc) b2_xy2(:,m-m1+1)=p%b2
        jb_yz(m-m1+1,n-n1+1)=p%jb(ix_loc-l1+1)
        if (m==iy_loc)  jb_xz(:,n-n1+1)=p%jb
        if (n==iz_loc)  jb_xy(:,m-m1+1)=p%jb
        if (n==iz2_loc) jb_xy2(:,m-m1+1)=p%jb
        if (bthresh_per_brms/=0) call calc_bthresh
        call vecout(41,trim(directory)//'/bvec',p%bb,bthresh,nbvec)
      endif
!
    endsubroutine daa_dt
!***********************************************************************
    subroutine df_diagnos_magnetic(f,df,p)
!
!  calculate diagnostics that involves df
!  Here we calculate <du/dt x b> and <u x db/dt>.
!  The latter is calculated as <divu dai/dt> -  <uji daj/dt>
!  This is used in dynamo theory for checking the minimal tau approximation.
!
!  10-oct-06/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: uudot,aadot,udotxb,B1_gradu
      real, dimension (nx) :: B1dot_udotxb,B1dot_uxbdot,B1dot_aadot,uxbdot2
!
      intent(in)  :: f, df, p
!
!  this routine is only called when ldiagnos=T
!  start with <du/dt x b>
!
      if (idiag_udotxbm/=0) then
        uudot=df(l1:l2,m,n,iux:iuz)
        call cross_mn(uudot,p%bb,udotxb)
        call dot_mn_sv(B1_ext,udotxb,B1dot_udotxb)
        call sum_mn_name(B1dot_udotxb,idiag_udotxbm)
      endif
!
!  next, do <divu dai/dt> -  <uji daj/dt>
!
      if (idiag_uxbdotm/=0) then
        aadot=df(l1:l2,m,n,iax:iaz)
        call dot_mn_sv(B1_ext,aadot,B1dot_aadot)
        call dot_mn_sm(B1_ext,p%uij,B1_gradu)
        call dot_mn(B1_gradu,aadot,uxbdot2)
        B1dot_uxbdot=p%divu*B1dot_aadot-uxbdot2
        call sum_mn_name(B1dot_uxbdot,idiag_uxbdotm)
      endif
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
      use Cdata
      use BorderProfiles, only: border_driving
      use Mpicomm, only: stop_it
      use Gravity, only: qgshear
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx,3) :: f_target
      real :: kr,kr1,Aphi,B0
      integer :: ju,j,i
!
      select case(borderaa)
!
      case('zero','0')
         f_target=0.
!
      case('toroidal')
        if (any(initaa=='uniform-Bphi')) &
             call stop_it("borderaa: this border profile "//&
             "is to be used with alfven_rphi (sinusoidal azimuthal wave) only")
        kr = 2*pi*rmode/(r_ext-r_int)
        kr1 = 1./kr
        do i=1,nx
          if ( ((p%rcyl_mn(i).ge.r_int).and.(p%rcyl_mn(i).le.r_int+2*wborder_int)).or.&
               ((p%rcyl_mn(i).ge.r_ext-2*wborder_ext).and.(p%rcyl_mn(i).le.r_ext))) then
            f_target(i,1) = 0.
            f_target(i,2) = 0.
            f_target(i,3) = -amplaa(1)*kr1*sin(kr*(p%rcyl_mn(i)-r_int))
          endif
        enddo
!
      case('Alfven-rz')
         kr = 2*pi*rmode/(r_ext-r_int)
         kr1 = 1./kr
         do i=1,nx
           if ( ((p%rcyl_mn(i).ge.r_int).and.(p%rcyl_mn(i).le.r_int+2*wborder_int)).or.&
                ((p%rcyl_mn(i).ge.r_ext-2*wborder_ext).and.(p%rcyl_mn(i).le.r_ext))) then
!
             Aphi =  amplaa(1)*kr1 * sin(kr*(p%rcyl_mn(i)-r_int)) + &
                  amplaa(1)*kr1**2*p%rcyl_mn1(i)*cos(kr*(p%rcyl_mn(i)-r_int))
!
             f_target(:,1) = Aphi * p%phix(i)
             f_target(:,2) = Aphi * p%phiy(i)
             f_target(:,3) = 0.
           endif
         enddo
!
      case('Alfven-zconst')
        B0=Lxyz(3)/(2*zmode*pi)
        do i=1,nx
          if ( ((p%rcyl_mn(i).ge.r_int).and.(p%rcyl_mn(i).le.r_int+2*wborder_int)).or.&
               ((p%rcyl_mn(i).ge.r_ext-2*wborder_ext).and.(p%rcyl_mn(i).le.r_ext))) then
!            
            if ((qgshear.eq.1.5).and.(rsmooth.eq.0.)) then
              !default unsmoothed keplerian
              Aphi=2*B0*sqrt(p%rcyl_mn1(i))
            else
              Aphi=B0/(p%rcyl_mn(i)*(2-qgshear))*(p%rcyl_mn(i)**2+rsmooth**2)**(1-qgshear/2.)
            endif
!
            f_target(:,1) = Aphi * p%phix(i)
            f_target(:,2) = Aphi * p%phiy(i)
            f_target(:,3) = 0.
          endif
        enddo
!
      case('nothing')
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
!   24-nov-03/dave: coded
!
      use Cdata
      use Sub, only: step, der_step
!
      type (pencil_case) :: p
      real, dimension (nx) :: eta_mn
      real, dimension (nx) :: prof,eta_r
      real, dimension (nx,3) :: geta
      real :: d_int=0.,d_ext=0.
!
      eta_r=0.
!
      if (eta_int > 0.) d_int=eta_int-eta
      if (eta_ext > 0.) d_ext=eta_ext-eta
!
!     calculate steps in resistivity
!
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
    endsubroutine eta_shell
!***********************************************************************
    subroutine calc_bthresh()
!
!  calculate bthresh from brms, give warnings if there are problems
!
!   6-aug-03/axel: coded
!
      use Cdata
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
    subroutine rescaling(f)
!
!  This routine could be turned into a wrapper routine later on,
!  if we want to do dynamic rescaling also on other quantities.
!
!  22-feb-05/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: scl
      integer :: j
!
      intent(inout) :: f
!
!  do rescaling only if brms is finite.
!  Note: we rely here on the brms that is update every it1 timesteps.
!  This may not always be sufficient.
!
      if (brms/=0) then
        scl=1.+rescaling_fraction*(brms_target/brms-1.)
        if (headtt) print*,'rescaling: scl=',scl
        do j=iax,iaz
          do n=n1,n2
            f(l1:l2,m1:m2,n,j)=scl*f(l1:l2,m1:m2,n,j)
          enddo
        enddo
      endif
!
    endsubroutine rescaling
!***********************************************************************
    subroutine calc_tau_aa_exterior(f,df)
!
!  magnetic field relaxation to zero on time scale tau_aa_exterior within
!  exterior region. For the time being this means z > zgrav.
!
!  29-jul-02/axel: coded
!
      use Cdata
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
    subroutine Omega_effect(f,df)
!
!  Omega effect coded (normally used in context of mean field theory)
!  Can do uniform shear (0,Sx,0), and the cosx*cosz profile (solar CZ).
!
!  30-apr-05/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      intent(in) :: f
      intent(inout) :: df
!
!  use gauge transformation, uxB = -Ay*grad(Uy) + gradient-term
!
      select case(Omega_profile)
      case('nothing'); print*,'Omega_profile=nothing'
      case('(0,Sx,0)')
        if (headtt) print*,'Omega_effect: uniform shear in x, S=',Omega_ampl
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)-Omega_ampl*f(l1:l2,m,n,iay)
      case('(Sz,0,0)')
        if (headtt) print*,'Omega_effect: uniform shear in z, S=',Omega_ampl
        df(l1:l2,m,n,iaz)=df(l1:l2,m,n,iaz)-Omega_ampl*f(l1:l2,m,n,iax)
        if (lhydro) df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-Omega_ampl*f(l1:l2,m,n,iuz)
      case('(0,cosx*cosz,0)')
        if (headtt) print*,'Omega_effect: solar shear, S=',Omega_ampl
        df(l1:l2,m,n,iax)=df(l1:l2,m,n,iax)+Omega_ampl*f(l1:l2,m,n,iay) &
            *sin(x(l1:l2))*cos(z(n))
        df(l1:l2,m,n,iaz)=df(l1:l2,m,n,iaz)+Omega_ampl*f(l1:l2,m,n,iay) &
            *cos(x(l1:l2))*sin(z(n))
      case default; print*,'Omega_profile=unknown'
      endselect
!
    endsubroutine Omega_effect
!***********************************************************************
    subroutine helflux(aa,uxb,jj)
!
!  magnetic helicity flux (preliminary)
!
!  14-aug-03/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (nx,3), intent(in) :: aa,uxb,jj
      real, dimension (nx,3) :: ee
      real, dimension (nx) :: FHx,FHz
      real :: FH
!
      ee=eta*jj-uxb
!
!  calculate magnetic helicity flux in the X and Z directions
!
      FHx=-2*ee(:,3)*aa(:,2)*dsurfyz
      FHz=+2*ee(:,1)*aa(:,2)*dsurfxy
!
!  sum up contribution per pencil
!  and then stuff result into surf_mn_name for summing up all processors.
!
      FH=FHx(nx)-FHx(1)
      if (ipz==0       .and.n==n1) FH=FH-sum(FHz)
      if (ipz==nprocz-1.and.n==n2) FH=FH+sum(FHz)
      call surf_mn_name(FH,idiag_exaym2)
!
    endsubroutine helflux
!***********************************************************************
    subroutine curflux(uxb,jj)
!
!  current helicity flux (preliminary)
!
!  27-nov-03/axel: adapted from helflux
!
      use Cdata
      use Sub
!
      real, dimension (nx,3), intent(in) :: uxb,jj
      real, dimension (nx,3) :: ee
      real, dimension (nx) :: FCx,FCz
      real :: FC
!
      ee=eta*jj-uxb
!
!  calculate current helicity flux in the X and Z directions
!  Could speed up by only calculating here boundary points!
!
      FCx=2*(ee(:,2)*jj(:,3)-ee(:,3)*jj(:,2))*dsurfyz
      FCz=2*(ee(:,1)*jj(:,2)-ee(:,2)*jj(:,1))*dsurfxy
!
!  sum up contribution per pencil
!  and then stuff result into surf_mn_name for summing up all processors.
!
      FC=FCx(nx)-FCx(1)
      if (ipz==0       .and.n==n1) FC=FC-sum(FCz)
      if (ipz==nprocz-1.and.n==n2) FC=FC+sum(FCz)
      call surf_mn_name(FC,idiag_exjm2)
!
    endsubroutine curflux
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
    subroutine forcing_continuous(df,p)
!
!  add a continuous forcing term (here currently only for localized rotors)
!
!  21-jan-07/axel: adapted from hydro
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: forcing_rhs,fxb
      real, dimension (nx) :: jf,phi
      real, dimension (mx), save :: phix,sinx,cosx
      real, dimension (my), save :: phiy,siny,cosy
      real, dimension (mz), save :: phiz,sinz,cosz
      real, save :: R2,R12
      type (pencil_case) :: p
      integer, save :: ifirst
      integer :: j,jfff,ifff
      real :: fact
!
!  at the first step, the sin and cos functions are calculated for all
!  x,y,z points and are then saved and used for all subsequent steps
!  and pencils
!
      if(ip<=6) print*,'forcing_continuous: ifirst=',ifirst
      if (ifirst==0) then
        if (iforcing_continuous_aa=='fixed_swirl') then
          if (lroot) print*,'forcing_continuous: fixed_swirl; swirl=',swirl
          R2=radius**2
          R12=1./R2
          phix=exp(-R12*x**2)
          phiy=exp(-R12*y**2)
          phiz=exp(-R12*z**2)
        elseif (iforcing_continuous_aa=='cosxcosz') then
          cosx=cos(k1x_ff*x)
          cosz=cos(k1z_ff*z)
        elseif (iforcing_continuous_aa=='RobertsFlow') then
          if (lroot) print*,'forcing_continuous: RobertsFlow'
          sinx=sin(k1_ff*x); cosx=cos(k1_ff*x)
          siny=sin(k1_ff*y); cosy=cos(k1_ff*y)
        endif
      endif
      ifirst=ifirst+1
      if(ip<=6) print*,'forcing_continuous: dt, ifirst=',dt,ifirst
!
!  calculate forcing
!
      if (iforcing_continuous_aa=='fixed_swirl') then
        fact=ampl_ff
        phi=2.*R12*fact*phix(l1:l2)*phiy(m)*phiz(n)
        forcing_rhs(:,1)=(-swirl*y(m    )+2.*x(l1:l2)*z(n))*phi
        forcing_rhs(:,2)=(+swirl*x(l1:l2)+2.*y(m    )*z(n))*phi
        forcing_rhs(:,3)=(R2-x(l1:l2)**2-y(m)**2)*2.*R12*phi
      elseif (iforcing_continuous_aa=='cosxcosz') then
        fact=ampl_ff
        forcing_rhs(:,1)=0.
        forcing_rhs(:,2)=fact*cosx(l1:l2)*cosz(n)
        forcing_rhs(:,3)=0.
      elseif (iforcing_continuous_aa=='RobertsFlow') then
        fact=ampl_ff
        forcing_rhs(:,1)=-fact*cosx(l1:l2)*siny(m)
        forcing_rhs(:,2)=+fact*sinx(l1:l2)*cosy(m)
        forcing_rhs(:,3)=+fact*cosx(l1:l2)*cosy(m)*sqrt(2.)
      endif
!
!  apply forcing in uncurled induction equation
!
      ifff=iax
      do j=1,3
        jfff=j+ifff-1
        df(l1:l2,m,n,jfff)=df(l1:l2,m,n,jfff)+forcing_rhs(:,j)
      enddo
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_jfm/=0) then
          call dot_mn(p%jj,forcing_rhs,jf)
          call sum_mn_name(jf,idiag_jfm)
        endif
      endif
!
    endsubroutine forcing_continuous
!***********************************************************************
    subroutine rprint_magnetic(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
      use Sub
!
      integer :: iname,inamex,inamey,inamez,ixy,ixz,irz,inamer
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
        idiag_b2m=0; idiag_bm2=0; idiag_j2m=0; idiag_jm2=0; idiag_abm=0
        idiag_jbm=0; idiag_ubm=0; idiag_epsM=0
        idiag_bxpt=0; idiag_bypt=0; idiag_bzpt=0; idiag_epsM_LES=0
        idiag_aybym2=0; idiag_exaym2=0; idiag_exjm2=0
        idiag_brms=0; idiag_bmax=0; idiag_jrms=0; idiag_jmax=0; idiag_vArms=0
        idiag_vAmax=0; idiag_dtb=0; idiag_arms=0; idiag_amax=0
        idiag_beta1m=0; idiag_beta1max=0
        idiag_bxm=0; idiag_bym=0; idiag_bzm=0
        idiag_bx2m=0; idiag_by2m=0; idiag_bz2m=0
        idiag_bxbymz=0; idiag_bxbzmz=0; idiag_bybzmz=0; idiag_b2mz=0
        idiag_bxbym=0; idiag_bxbzm=0; idiag_bybzm=0; idiag_djuidjbim=0
        idiag_bxmz=0; idiag_bymz=0; idiag_bzmz=0; idiag_bmx=0; idiag_bmy=0
        idiag_bx2mz=0; idiag_by2mz=0; idiag_bz2mz=0
        idiag_bmz=0; idiag_bxmxy=0; idiag_bymxy=0; idiag_bzmxy=0
        idiag_uxbm=0; idiag_oxuxbm=0; idiag_jxbxbm=0.; idiag_gpxbm=0.
        idiag_uxDxuxbm=0.; idiag_uxbmx=0; idiag_uxbmy=0; idiag_uxbmz=0
        idiag_uxjm=0; idiag_ujxbm=0
        idiag_b3b21m=0; idiag_b1b32m=0; idiag_b2b13m=0
        idiag_EMFdotBm=0
        idiag_udotxbm=0; idiag_uxbdotm=0
        idiag_brmphi=0; idiag_bpmphi=0; idiag_bzmphi=0; idiag_b2mphi=0
        idiag_jbmphi=0; idiag_uxbrmphi=0; idiag_uxbpmphi=0; idiag_uxbzmphi=0
        idiag_jxbrmphi=0; idiag_jxbpmphi=0; idiag_jxbzmphi=0
        idiag_armphi=0; idiag_apmphi=0; idiag_azmphi=0
        idiag_dteta=0; idiag_uxBrms=0; idiag_Bresrms=0; idiag_Rmrms=0
        idiag_jfm=0; idiag_brbpmr=0; idiag_va2m=0
        idiag_b2mr=0; idiag_brmr=0; idiag_bpmr=0; idiag_bzmr=0
        idiag_armr=0; idiag_apmr=0; idiag_azmr=0
        idiag_bxmx=0; idiag_bymy=0
        idiag_mflux_x=0; idiag_mflux_y=0; idiag_mflux_z=0
        idiag_bmxy_rms=0
!
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dteta',idiag_dteta)
        call parse_name(iname,cname(iname),cform(iname),'aybym2',idiag_aybym2)
        call parse_name(iname,cname(iname),cform(iname),'exaym2',idiag_exaym2)
        call parse_name(iname,cname(iname),cform(iname),'exjm2',idiag_exjm2)
        call parse_name(iname,cname(iname),cform(iname),'abm',idiag_abm)
        call parse_name(iname,cname(iname),cform(iname),'jbm',idiag_jbm)
        call parse_name(iname,cname(iname),cform(iname),'ubm',idiag_ubm)
        call parse_name(iname,cname(iname),cform(iname),'b2m',idiag_b2m)
        call parse_name(iname,cname(iname),cform(iname),'bm2',idiag_bm2)
        call parse_name(iname,cname(iname),cform(iname),'j2m',idiag_j2m)
        call parse_name(iname,cname(iname),cform(iname),'jm2',idiag_jm2)
        call parse_name(iname,cname(iname),cform(iname),'epsM',idiag_epsM)
        call parse_name(iname,cname(iname),cform(iname),&
            'epsM_LES',idiag_epsM_LES)
        call parse_name(iname,cname(iname),cform(iname),'brms',idiag_brms)
        call parse_name(iname,cname(iname),cform(iname),'bmax',idiag_bmax)
        call parse_name(iname,cname(iname),cform(iname),'jrms',idiag_jrms)
        call parse_name(iname,cname(iname),cform(iname),'jmax',idiag_jmax)
        call parse_name(iname,cname(iname),cform(iname),'arms',idiag_arms)
        call parse_name(iname,cname(iname),cform(iname),'amax',idiag_amax)
        call parse_name(iname,cname(iname),cform(iname),'vArms',idiag_vArms)
        call parse_name(iname,cname(iname),cform(iname),'vAmax',idiag_vAmax)
        call parse_name(iname,cname(iname),cform(iname),'vA2m',idiag_vA2m)
        call parse_name(iname,cname(iname),cform(iname),&
            'beta1m',idiag_beta1m)
        call parse_name(iname,cname(iname),cform(iname),&
            'beta1max',idiag_beta1max)
        call parse_name(iname,cname(iname),cform(iname),'dtb',idiag_dtb)
        call parse_name(iname,cname(iname),cform(iname),'bxm',idiag_bxm)
        call parse_name(iname,cname(iname),cform(iname),'bym',idiag_bym)
        call parse_name(iname,cname(iname),cform(iname),'bzm',idiag_bzm)
        call parse_name(iname,cname(iname),cform(iname),'bx2m',idiag_bx2m)
        call parse_name(iname,cname(iname),cform(iname),'by2m',idiag_by2m)
        call parse_name(iname,cname(iname),cform(iname),'bz2m',idiag_bz2m)
        call parse_name(iname,cname(iname),cform(iname),'bxbym',idiag_bxbym)
        call parse_name(iname,cname(iname),cform(iname),'bxbzm',idiag_bxbzm)
        call parse_name(iname,cname(iname),cform(iname),'bybzm',idiag_bybzm)
        call parse_name(iname,cname(iname),cform(iname),&
            'djuidjbim',idiag_djuidjbim)
        call parse_name(iname,cname(iname),cform(iname),'uxbm',idiag_uxbm)
        call parse_name(iname,cname(iname),cform(iname),'uxbmx',idiag_uxbmx)
        call parse_name(iname,cname(iname),cform(iname),'uxbmy',idiag_uxbmy)
        call parse_name(iname,cname(iname),cform(iname),'uxbmz',idiag_uxbmz)
        call parse_name(iname,cname(iname),cform(iname),'uxjm',idiag_uxjm)
        call parse_name(iname,cname(iname),cform(iname),'ujxbm',idiag_ujxbm)
        call parse_name(iname,cname(iname),cform(iname),'jxbxbm',idiag_jxbxbm)
        call parse_name(iname,cname(iname),cform(iname),'oxuxbm',idiag_oxuxbm)
        call parse_name(iname,cname(iname),cform(iname),'gpxbm',idiag_gpxbm)
        call parse_name(iname,cname(iname),cform(iname),&
            'uxDxuxbm',idiag_uxDxuxbm)
        call parse_name(iname,cname(iname),cform(iname),'b3b21m',idiag_b3b21m)
        call parse_name(iname,cname(iname),cform(iname),'b1b32m',idiag_b1b32m)
        call parse_name(iname,cname(iname),cform(iname),'b2b13m',idiag_b2b13m)
        call parse_name(iname,cname(iname),cform(iname),'EMFdotBm',idiag_EMFdotBm)
        call parse_name(iname,cname(iname),cform(iname),'udotxbm',idiag_udotxbm)
        call parse_name(iname,cname(iname),cform(iname),'uxbdotm',idiag_uxbdotm)
        call parse_name(iname,cname(iname),cform(iname),'bmx',idiag_bmx)
        call parse_name(iname,cname(iname),cform(iname),'bmy',idiag_bmy)
        call parse_name(iname,cname(iname),cform(iname),'bmz',idiag_bmz)
        call parse_name(iname,cname(iname),cform(iname),'bxpt',idiag_bxpt)
        call parse_name(iname,cname(iname),cform(iname),'bypt',idiag_bypt)
        call parse_name(iname,cname(iname),cform(iname),'bzpt',idiag_bzpt)
        call parse_name(iname,cname(iname),cform(iname),'uxBrms',idiag_uxBrms)
        call parse_name(iname,cname(iname),cform(iname),'Bresrms',idiag_Bresrms)
        call parse_name(iname,cname(iname),cform(iname),'Rmrms',idiag_Rmrms)
        call parse_name(iname,cname(iname),cform(iname),'jfm',idiag_jfm)
        call parse_name(iname,cname(iname),cform(iname),'bmxy_rms',idiag_bmxy_rms)
!
      enddo
!
! Currently need to force zaverage calculation at every lout step for
! bmx and bmy.
!
      if ((idiag_bmx+idiag_bmy)>0) ldiagnos_need_zaverages=.true.
!
!  check for those quantities for which we want xy-averages
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'bxmx',idiag_bxmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex), &
                        'mflux_x',idiag_mflux_x)
      enddo
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'bymy',idiag_bymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey), &
                        'mflux_y',idiag_mflux_y)
      enddo
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bxmz',idiag_bxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bymz',idiag_bymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bzmz',idiag_bzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bx2mz',idiag_bx2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'by2mz',idiag_by2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bz2mz',idiag_bz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bxbymz',idiag_bxbymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bxbzmz',idiag_bxbzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bybzmz',idiag_bybzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'b2mz',idiag_b2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez), &
                        'mflux_z',idiag_mflux_z)
      enddo
!
!  check for those quantities for which we want y-averages
!
      do ixz=1,nnamexz
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bxmxz',idiag_bxmxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bymxz',idiag_bymxz)
        call parse_name(ixz,cnamexz(ixz),cformxz(ixz),'bzmxz',idiag_bzmxz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bxmxy',idiag_bxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bymxy',idiag_bymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'bzmxy',idiag_bzmxy)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'brmphi'  ,idiag_brmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'bpmphi'  ,idiag_bpmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'bzmphi'  ,idiag_bzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'b2mphi'  ,idiag_b2mphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'jbmphi'  ,idiag_jbmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uxbrmphi',idiag_uxbrmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uxbpmphi',idiag_uxbpmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'uxbzmphi',idiag_uxbzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'jxbrmphi',idiag_jxbrmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'jxbpmphi',idiag_jxbpmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'jxbzmphi',idiag_jxbzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'armphi'  ,idiag_armphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'apmphi'  ,idiag_apmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'azmphi'  ,idiag_azmphi)

      enddo
!
!  check for those quantities for which we want phiz-averages
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'brmr',  idiag_brmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'bpmr',  idiag_bpmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'bzmr',  idiag_bzmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'armr',  idiag_armr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'apmr',  idiag_apmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'azmr',  idiag_azmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'b2mr',  idiag_b2mr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'brbpmr',idiag_brbpmr)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'i_dteta=',idiag_dteta
        write(3,*) 'i_aybym2=',idiag_aybym2
        write(3,*) 'i_exaym2=',idiag_exaym2
        write(3,*) 'i_exjm2=',idiag_exjm2
        write(3,*) 'i_abm=',idiag_abm
        write(3,*) 'i_jbm=',idiag_jbm
        write(3,*) 'i_ubm=',idiag_ubm
        write(3,*) 'i_b2m=',idiag_b2m
        write(3,*) 'i_bm2=',idiag_bm2
        write(3,*) 'i_j2m=',idiag_j2m
        write(3,*) 'i_jm2=',idiag_jm2
        write(3,*) 'i_epsM=',idiag_epsM
        write(3,*) 'i_epsM_LES=',idiag_epsM_LES
        write(3,*) 'i_brms=',idiag_brms
        write(3,*) 'i_bmax=',idiag_bmax
        write(3,*) 'i_jrms=',idiag_jrms
        write(3,*) 'i_jmax=',idiag_jmax
        write(3,*) 'i_arms=',idiag_arms
        write(3,*) 'i_amax=',idiag_amax
        write(3,*) 'i_vArms=',idiag_vArms
        write(3,*) 'i_vAmax=',idiag_vAmax
        write(3,*) 'i_vA2m=',idiag_vA2m
        write(3,*) 'i_beta1m=',idiag_beta1m
        write(3,*) 'i_beta1max=',idiag_beta1max
        write(3,*) 'i_dtb=',idiag_dtb
        write(3,*) 'i_bxm=',idiag_bxm
        write(3,*) 'i_bym=',idiag_bym
        write(3,*) 'i_bzm=',idiag_bzm
        write(3,*) 'i_bx2m=',idiag_bx2m
        write(3,*) 'i_by2m=',idiag_by2m
        write(3,*) 'i_bz2m=',idiag_bz2m
        write(3,*) 'i_bxbym=',idiag_bxbym
        write(3,*) 'i_bxbzm=',idiag_bxbzm
        write(3,*) 'i_bybzm=',idiag_bybzm
        write(3,*) 'i_djuidjbim=',idiag_djuidjbim
        write(3,*) 'i_uxbm=',idiag_uxbm
        write(3,*) 'i_uxbmx=',idiag_uxbmx
        write(3,*) 'i_uxbmy=',idiag_uxbmy
        write(3,*) 'i_uxbmz=',idiag_uxbmz
        write(3,*) 'i_uxjm=',idiag_uxjm
        write(3,*) 'i_ujxbm=',idiag_ujxbm
        write(3,*) 'i_oxuxbm=',idiag_oxuxbm
        write(3,*) 'i_jxbxbm=',idiag_jxbxbm
        write(3,*) 'i_gpxbm=',idiag_gpxbm
        write(3,*) 'i_uxDxuxbm=',idiag_uxDxuxbm
        write(3,*) 'i_b3b21m=',idiag_b3b21m
        write(3,*) 'i_b1b32m=',idiag_b1b32m
        write(3,*) 'i_b2b13m=',idiag_b2b13m
        write(3,*) 'i_EMFdotBm=',idiag_EMFdotBm
        write(3,*) 'i_udotxbm=',idiag_udotxbm
        write(3,*) 'i_uxbdotm=',idiag_uxbdotm
        write(3,*) 'i_bxmz=',idiag_bxmz
        write(3,*) 'i_bymz=',idiag_bymz
        write(3,*) 'i_bzmz=',idiag_bzmz
        write(3,*) 'i_bx2mz=',idiag_bxmz
        write(3,*) 'i_by2mz=',idiag_bymz
        write(3,*) 'i_bz2mz=',idiag_bzmz
        write(3,*) 'i_bxbymz=',idiag_bxbymz
        write(3,*) 'i_b2mz=',idiag_b2mz
        write(3,*) 'i_bmx=',idiag_bmx
        write(3,*) 'i_bmy=',idiag_bmy
        write(3,*) 'i_bmz=',idiag_bmz
        write(3,*) 'i_bxpt=',idiag_bxpt
        write(3,*) 'i_bypt=',idiag_bypt
        write(3,*) 'i_bzpt=',idiag_bzpt
        write(3,*) 'i_bxmxy=',idiag_bxmxy
        write(3,*) 'i_bymxy=',idiag_bymxy
        write(3,*) 'i_bzmxy=',idiag_bzmxy
        write(3,*) 'i_bxmxz=',idiag_bxmxz
        write(3,*) 'i_bymxz=',idiag_bymxz
        write(3,*) 'i_bzmxz=',idiag_bzmxz
        write(3,*) 'i_brmphi=',idiag_brmphi
        write(3,*) 'i_bpmphi=',idiag_bpmphi
        write(3,*) 'i_bzmphi=',idiag_bzmphi
        write(3,*) 'i_b2mphi=',idiag_b2mphi
        write(3,*) 'i_armphi=',idiag_armphi
        write(3,*) 'i_apmphi=',idiag_apmphi
        write(3,*) 'i_azmphi=',idiag_azmphi
        write(3,*) 'i_brmr=',idiag_brmr
        write(3,*) 'i_bpmr=',idiag_bpmr
        write(3,*) 'i_bzmr=',idiag_bzmr
        write(3,*) 'i_armr=',idiag_armr
        write(3,*) 'i_apmr=',idiag_apmr
        write(3,*) 'i_azmr=',idiag_azmr
        write(3,*) 'i_b2mr=',idiag_b2mr
        write(3,*) 'i_jbmphi=',idiag_jbmphi
        write(3,*) 'i_uxBrms=',idiag_uxBrms
        write(3,*) 'i_Bresrms=',idiag_Bresrms
        write(3,*) 'i_Rmrms=',idiag_Rmrms
        write(3,*) 'i_jfm=',idiag_jfm
        write(3,*) 'i_bxmx=',idiag_bxmx
        write(3,*) 'i_bymy=',idiag_bymy
        write(3,*) 'i_mflux_x=',idiag_mflux_x
        write(3,*) 'i_mflux_y=',idiag_mflux_y
        write(3,*) 'i_mflux_z=',idiag_mflux_z
        write(3,*) 'nname=',nname
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'nnamexz=',nnamexz
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'iaa=',iaa
        write(3,*) 'iax=',iax
        write(3,*) 'iay=',iay
        write(3,*) 'iaz=',iaz
        write(3,*) 'ihypres=',ihypres
        write(3,*) 'i_bmxy_rms=',idiag_bmxy_rms
      endif
!
    endsubroutine rprint_magnetic
!***********************************************************************
    subroutine get_slices_magnetic(f,slices)
!
!  Write slices for animation of magnetic variables.
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
          if (slices%index >= 3) then
            slices%ready = .false.
          else
            slices%yz=f(slices%ix,m1:m2    ,n1:n2,iax+slices%index)
            slices%xz=f(l1:l2    ,slices%iy,n1:n2,iax+slices%index)
            slices%xy=f(l1:l2    ,m1:m2    ,slices%iz,iax+slices%index)
            slices%xy2=f(l1:l2    ,m1:m2    ,slices%iz2,iax+slices%index)
            slices%index = slices%index+1
            if (slices%index < 3) slices%ready = .true.
          endif
!
!  Magnetic field (derived variable)
!
        case ('bb')
          if (slices%index >= 3) then
            slices%ready = .false.
          else
            slices%index = slices%index+1
            slices%yz=>bb_yz(:,:,slices%index)
            slices%xz=>bb_xz(:,:,slices%index)
            slices%xy=>bb_xy(:,:,slices%index)
            slices%xy2=>bb_xy2(:,:,slices%index)
            if (slices%index < 3) slices%ready = .true.
          endif
!
!  Magnetic field (derived variable)
!
        case ('jj')
          if (slices%index >= 3) then
            slices%ready = .false.
          else
            slices%index = slices%index+1
            slices%yz=>jj_yz(:,:,slices%index)
            slices%xz=>jj_xz(:,:,slices%index)
            slices%xy=>jj_xy(:,:,slices%index)
            slices%xy2=>jj_xy2(:,:,slices%index)
            if (slices%index < 3) slices%ready = .true.
          endif
!
!  Magnetic field squared (derived variable)
!
        case ('b2')
          slices%yz=>b2_yz
          slices%xz=>b2_xz
          slices%xy=>b2_xy
          slices%xy2=>b2_xy2
          slices%ready = .true.
!
!  Current density (derived variable)
!
        case ('jb')
          slices%yz=>jb_yz
          slices%xz=>jb_xz
          slices%xy=>jb_xy
          slices%xy2=>jb_xy2
          slices%ready = .true.
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
      use Cdata
      use Mpicomm
      use Sub
!
      logical,save :: first=.true.
      real, dimension(nx) :: bymx,bzmx
      real, dimension(ny,nprocy) :: bxmy,bzmy
      real :: bmx,bmy,bmz,bxmxy,bymxy,bzmxy,bmxy_rms
      integer :: l,j
!
!  For vector output (of bb vectors) we need brms
!  on all processors. It suffices to have this for times when lout=.true.,
!  but we need to broadcast the result to all procs.
!
!  calculate brms (this requires that brms is set in print.in)
!  broadcast result to other processors
!
      if (idiag_brms/=0) then
        if (iproc==0) brms=fname(idiag_brms)
        call mpibcast_real(brms,1)
      endif

      if (.not.lroot) return
!
!  Magnetic energy in vertically averaged field!
!  The bymxy and bzmxy must have been calculated,
!  so they are present on the root processor.
!
      if (idiag_bmx/=0) then
        if (idiag_bymxy==0.or.idiag_bzmxy==0) then
          if (first) print*,"calc_mfield:                  WARNING"
          if (first) print*, &
                  "calc_mfield: NOTE: to get bmx, bymxy and bzmxy must also be set in zaver"
          if (first) print*, &
                  "calc_mfield:       We proceed, but you'll get bmx=0"
          bmx=0.
        else
          do l=1,nx
            bymx(l)=sum(fnamexy(l,:,:,idiag_bymxy))/(ny*nprocy)
            bzmx(l)=sum(fnamexy(l,:,:,idiag_bzmxy))/(ny*nprocy)
          enddo
          bmx=sqrt(sum(bymx**2+bzmx**2)/nx)
        endif
        call save_name(bmx,idiag_bmx)
      endif
!  similarly for bmy
!
      if (idiag_bmy/=0) then
        if (idiag_bxmxy==0.or.idiag_bzmxy==0) then
          if (first) print*,"calc_mfield:                  WARNING"
          if (first) print*, &
                  "calc_mfield: NOTE: to get bmy, bxmxy and bzmxy must also be set in zaver"
          if (first) print*, &
                  "calc_mfield:       We proceed, but you'll get bmy=0"
          bmy=0.
        else
          do j=1,nprocy
          do m=1,ny
            bxmy(m,j)=sum(fnamexy(:,m,j,idiag_bxmxy))/nx
            bzmy(m,j)=sum(fnamexy(:,m,j,idiag_bzmxy))/nx
          enddo
          enddo
          bmy=sqrt(sum(bxmy**2+bzmy**2)/(ny*nprocy))
        endif
        call save_name(bmy,idiag_bmy)
      endif
!
!  Magnetic energy in horizontally averaged field
!  The bxmz and bymz must have been calculated,
!  so they are present on the root processor.
!
      if (idiag_bmz/=0) then
        if (idiag_bxmz==0.or.idiag_bymz==0) then
          if (first) print*,"calc_mfield:                  WARNING"
          if (first) print*, &
                  "calc_mfield: NOTE: to get bmz, bxmz and bymz must also be set in xyaver"
          if (first) print*, &
                  "calc_mfield:       This may be because we renamed zaver.in into xyaver.in"
          if (first) print*, &
                  "calc_mfield:       We proceed, but you'll get bmz=0"
          bmz=0.
        else
          bmz=sqrt(sum(fnamez(:,:,idiag_bxmz)**2+fnamez(:,:,idiag_bymz)**2)/(nz*nprocz))
        endif
        call save_name(bmz,idiag_bmz)
      endif
!
!  Magnetic energy in z averaged field 
!  The bxmxy, bymxy and bzmxy must have been calculated,
!  so they are present on the root processor.
!
      if (idiag_bmxy_rms/=0) then
        if (idiag_bxmxy==0.or.idiag_bymxy==0.or.idiag_bzmxy==0) then
          if (first) print*,"calc_mfield:                  WARNING"
          if (first) print*, &
                  "calc_mfield: NOTE: to get bmxy_rms, bxmxy, bymxy and bzmxy must also be set in zaver"
          if (first) print*, &
                  "calc_mfield:       We proceed, but you'll get bmxy_rms=0"
          bmxy_rms=0.
        else
          bmxy_rms=0.
          do l=1,nx
            do m=1,ny
              bxmxy=sum(fnamexy(l,m,:,idiag_bxmxy))/nprocy
              bymxy=sum(fnamexy(l,m,:,idiag_bymxy))/nprocy
              bzmxy=sum(fnamexy(l,m,:,idiag_bzmxy))/nprocy
              bmxy_rms = bmxy_rms+bxmxy**2+bymxy**2+bzmxy**2
            enddo
          enddo
          bmxy_rms = bmxy_rms/(nx*ny)
          bmxy_rms = sqrt(bmxy_rms)
        endif
        call save_name(bmxy_rms,idiag_bmxy_rms)
      endif

      first = .false.
    endsubroutine calc_mfield
!***********************************************************************
    subroutine alfven_x(ampl,f,iuu,iaa,ilnrho,xx,kx)
!
!  Alfven wave propagating in the x-direction
!  ux = +sin(kx-ot), for B0x=1 and rho=1.
!  Az = -cos(kx-ot), ie By = sin(kx-ot)
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
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: xx
      real :: ampl,kx
      integer :: iuu,iaa,ilnrho
!
!  ux and Ay.
!  Don't overwrite the density, just add to the log of it.
!
      f(:,:,:,ilnrho)=ampl*sin(kx*xx)+f(:,:,:,ilnrho)
      f(:,:,:,iuu+0)=+ampl*sin(kx*xx)
      f(:,:,:,iuu+1)=+ampl*sin(kx*xx)
      f(:,:,:,iaa+2)=-ampl*cos(kx*xx)
!
    endsubroutine alfven_x
!***********************************************************************
    subroutine alfven_y(ampl,f,iuu,iaa,yy,ky,mu0)
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
      real, dimension (mx,my,mz)         :: yy
      real                               :: ampl,ky,mu0
      integer                            :: iuu,iaa
!
!  ux and Az
!
      f(:,:,:,iuu+0) = +ampl*cos(ky*yy)
      f(:,:,:,iaa+2) = -ampl*sin(ky*yy)*sqrt(mu0)/ky
!
    endsubroutine alfven_y
!***********************************************************************
    subroutine alfven_z(ampl,f,iuu,iaa,zz,kz,mu0)
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
      real, dimension (mx,my,mz) :: zz
      real :: ampl,kz,mu0
      integer :: iuu,iaa
!
!  ux and Ay
!
      f(:,:,:,iuu+0)=+ampl*cos(kz*zz)
      f(:,:,:,iaa+1)=+ampl*sin(kz*zz)*sqrt(mu0)
!
    endsubroutine alfven_z
!***********************************************************************
    subroutine alfven_xy(ampl,f,iuu,iaa,xx,yy,kx,ky)
!
!  Alfven wave propagating in the xy-direction; can be used in 2-d runs.
!  uz = cos(kx*x+ky*y-ot), for B0=(1,1,0) and rho=1.
!  Ax = sin(kx*x+ky*y-ot),
!  Ay = sin(kx*x+ky*y-ot),
!
!  16-jun-07/axel: adapted from alfven_y
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz)         :: xx,yy
      real                               :: ampl,kx,ky,mu0=1.,om
      integer                            :: iuu,iaa
!
!  set ux, Ax, and Ay
!
      om=B_ext(1)*kx+B_ext(2)*ky
      f(:,:,:,iuu+2)=+ampl*cos(kx*xx+ky*yy)
      f(:,:,:,iaa+0)=+ampl*sin(kx*xx+ky*yy)*sqrt(mu0)/om*B_ext(2)
      f(:,:,:,iaa+1)=-ampl*sin(kx*xx+ky*yy)*sqrt(mu0)/om*B_ext(1)
!
    endsubroutine alfven_xy
!***********************************************************************
    subroutine alfven_xz(ampl,f,iuu,iaa,xx,zz,kx,kz)
!
!  Alfven wave propagating in the xz-direction; can be used in 2-d runs.
!  uz = cos(kx*x+kz*z-ot), for B0=(1,1,0) and rho=1.
!  Ax = sin(kx*x+kz*z-ot),
!  Az = sin(kx*x+kz*z-ot),
!
!  16-jun-07/axel: adapted from alfven_xy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz)         :: xx,zz
      real                               :: ampl,kx,kz,mu0=1.,om
      integer                            :: iuu,iaa
!
!  set ux, Ax, and Az
!
      om=B_ext(1)*kx+B_ext(3)*kz
      f(:,:,:,iuu+2)=+ampl*cos(kx*xx+kz*zz)
      f(:,:,:,iaa+0)=+ampl*sin(kx*xx+kz*zz)*sqrt(mu0)/om*B_ext(2)
      f(:,:,:,iaa+2)=-ampl*sin(kx*xx+kz*zz)*sqrt(mu0)/om*B_ext(1)
!
    endsubroutine alfven_xz
!***********************************************************************
    subroutine alfven_rphi(B0,f,xx,yy,mode)
!
!  Alfven wave propagating on radial direction with
!  field pointing to the phi direction.
!
!  Bphi = B0 cos(k r) ==> Az = -1/k B0 sin(k r)
!
!  04-oct-06/wlad: coded
!
      use Cdata
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz) :: xx,yy,rrcyl
      real :: B0,kr,mode
!
      kr = 2*pi*mode/(r_ext-r_int)
      rrcyl = sqrt(xx**2 + yy**2)
!
      f(:,:,:,iaz) =  -B0/kr*sin(kr*(rrcyl-r_int))
!
    endsubroutine alfven_rphi
!***********************************************************************
    subroutine alfven_zconst(f,xx,yy)
!
!  Radially variable field pointing in the z direction
!  4 Balbus Hawley wavelengths in the vertical direction
!
!  Bz=Lz/(8pi)*Omega      ==> Aphi = Lz/(8pi) Omega*r/(2-q)
!
!  The smoothed case should be general, since it reduces 
!  to the non-smoothed for r0_pot=0.
!
!  B=C*(r2+r02)^-q ==> Aphi=C/(r*(2-q))*(r2+r02)^(1-q/2)
!
!  04-oct-06/wlad: coded
!
      use Cdata
      use Gravity, only: qgshear,r0_pot
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz) :: xx,yy,rrcyl,Aphi
      real :: B0
!
      B0=Lxyz(3)/(2*zmode*pi)
      rrcyl=sqrt(xx**2+yy**2)
      Aphi=B0/(rrcyl*(2-qgshear))*(rrcyl**2+r0_pot**2)**(1-qgshear/2.)
!
      f(:,:,:,iax) =  -Aphi*yy/rrcyl
      f(:,:,:,iay) =   Aphi*xx/rrcyl
!
    endsubroutine alfven_zconst
!***********************************************************************
    subroutine alfven_rz(B0,f,xx,yy,mode)
!
!  Alfven wave propagating on radial direction with
!  field pointing to the z direction.
!
!  Bz = B0 cos(k r) ==> Aphi = B0/k sin(k r) + B0/(k2*r)*cos(k r)
!
!  04-oct-06/wlad: coded
!
      use Cdata
      use Mpicomm,only:stop_it
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz) :: xx,yy,rrcyl,Aphi
      real :: B0,kr,mode,k1,const
!
      if (headtt) print*,'radial alfven wave propagating on z direction'
      if (.not.lcylinder_in_a_box) &
           call stop_it("alfven_rz, this initial condition works only for embedded cylinders")
!
! Choose between cases. The non-smoothed case is singular in r=0 for the potential, thus
! can only be used if freezing is used with a non-zero r_int. The smoothed case
! has a linear component that prevents singularities and a exponential that prevents
! the field from growing large in the outer disk
!
      rrcyl = max(sqrt(xx**2 + yy**2),tini)
      if (r_int.gt.0.) then
         if (lroot) print*,'freezing is being used, ok to use singular potentials'
         if (lroot) print*,'Bz=B0cos(k.r) ==> Aphi=B0/k*sin(k.r)+B0/(k^2*r)*cos(k r)'
         kr = 2*pi*mode/(r_ext-r_int)
         Aphi =  B0/kr * sin(kr*(rrcyl-r_int)) + &
              B0/(kr**2*rrcyl)*cos(kr*(rrcyl-r_int))
      else   
         if (lroot) print*,'Softened magnetic field in the center'
         if (mode .lt. 5) call stop_it("put more wavelengths in the field")
          kr = 2*pi*mode/r_ext
          k1 = 1. !not tested for other values
          const=B0*exp(1.)*k1/cos(kr/k1)
          Aphi=const/kr*rrcyl*exp(-k1*rrcyl)*sin(kr*rrcyl)
      endif
!
      f(:,:,:,iax) = Aphi * (-yy/rrcyl)
      f(:,:,:,iay) = Aphi * ( xx/rrcyl)
!
    endsubroutine alfven_rz
!***********************************************************************
    subroutine alfvenz_rot(ampl,f,iuu,iaa,zz,kz,O)
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
      use Cdata, only: lroot
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: zz
      real :: ampl,kz,O,fac
      integer :: iuu,iaa
!
!  ux, uy, Ax and Ay
!
      if (lroot) print*,'alfvenz_rot: Alfven wave with rotation; O,kz=',O,kz
      fac=-O+sqrt(O**2+kz**2)
      f(:,:,:,iuu+0)=-ampl*sin(kz*zz)*fac/kz
      f(:,:,:,iuu+1)=-ampl*cos(kz*zz)*fac/kz
      f(:,:,:,iaa+0)=+ampl*sin(kz*zz)/kz
      f(:,:,:,iaa+1)=+ampl*cos(kz*zz)/kz
!
    endsubroutine alfvenz_rot
!***********************************************************************
    subroutine alfvenz_rot_shear(ampl,f,iuu,iaa,zz,kz,O)
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
      use Cdata, only: lroot
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: zz
      real :: ampl,kz,O
      complex :: fac
      integer :: iuu,iaa
!
!  ux, uy, Ax and Ay
!
      if (lroot) print*,'alfvenz_rot_shear: '// &
          'Alfven wave with rotation and shear; O,kz=',O,kz
      fac=cmplx(O-sqrt(16*kz**2+O**2),0.)
      f(:,:,:,iuu+0)=f(:,:,:,iuu+0) + ampl*fac/(4*kz)*sin(kz*zz)
      f(:,:,:,iuu+1)=f(:,:,:,iuu+1) + ampl*real(exp(cmplx(0,zz*kz))* &
          fac*sqrt(2*kz**2+O*fac)/(sqrt(2.)*kz*(-6*O-fac)))
      f(:,:,:,iaa+0)=ampl*sin(kz*zz)/kz
      f(:,:,:,iaa+1)=-ampl*2*sqrt(2.)*aimag(exp(cmplx(0,zz*kz))* &
          sqrt(2*kz**2+O*fac)/(-6*O-fac)/(cmplx(0,kz)))
!
    endsubroutine alfvenz_rot_shear
!***********************************************************************
    subroutine fluxrings(ampl,f,ivar,xx,yy,zz,profile)
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
      use Cdata, only: lroot, epsi
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,3)    :: tmpv
      real, dimension (mx,my,mz)      :: xx,yy,zz,xx1,yy1,zz1
      real, dimension(3) :: axis,disp
      real    :: ampl,phi,theta,ct,st,cp,sp
      real    :: fring,Iring,R0,width
      integer :: i,ivar
      character (len=*), optional :: profile
      character (len=labellen) :: prof
!
      if (present(profile)) then
        prof = profile
      else
        prof = 'tanh'
      endif

      if (any((/fring1,fring2,Iring1,Iring2/) /= 0.)) then
        ! fringX is the magnetic flux, IringX the current
        if (lroot) then
          print*, 'fluxrings: Initialising magnetic flux rings'
        endif
        do i=1,2
          if (i==1) then
            fring = fring1      ! magnetic flux along ring
            Iring = Iring1      ! current along ring (for twisted flux tube)
            R0    = Rring1      ! radius of ring
            width = wr1         ! ring thickness
            axis  = axisr1 ! orientation
            disp  = dispr1    ! position
          else
            fring = fring2
            Iring = Iring2
            R0    = Rring2
            width = wr2
            axis  = axisr2
            disp  = dispr2
          endif
          phi   = atan2(axis(2),axis(1)+epsi)
          theta = atan2(sqrt(axis(1)**2+axis(2)**2)+epsi,axis(3))
          ct = cos(theta); st = sin(theta)
          cp = cos(phi)  ; sp = sin(phi)
          ! Calculate D^(-1)*(xxx-disp)
          xx1 =  ct*cp*(xx-disp(1)) + ct*sp*(yy-disp(2)) - st*(zz-disp(3))
          yy1 = -   sp*(xx-disp(1)) +    cp*(yy-disp(2))
          zz1 =  st*cp*(xx-disp(1)) + st*sp*(yy-disp(2)) + ct*(zz-disp(3))
          call norm_ring(xx1,yy1,zz1,fring,Iring,R0,width,tmpv,PROFILE=prof)
          ! calculate D*tmpv
          f(:,:,:,ivar  ) = f(:,:,:,ivar  ) + ampl*( &
               + ct*cp*tmpv(:,:,:,1) - sp*tmpv(:,:,:,2) + st*cp*tmpv(:,:,:,3))
          f(:,:,:,ivar+1) = f(:,:,:,ivar+1) + ampl*( &
               + ct*sp*tmpv(:,:,:,1) + cp*tmpv(:,:,:,2) + st*sp*tmpv(:,:,:,3))
          f(:,:,:,ivar+2) = f(:,:,:,ivar+2) + ampl*( &
               - st   *tmpv(:,:,:,1)                    + ct   *tmpv(:,:,:,3))
        enddo
      endif
      if (lroot) print*, 'fluxrings: Magnetic flux rings initialized'
!
    endsubroutine fluxrings
!***********************************************************************
    subroutine norm_ring(xx,yy,zz,fring,Iring,R0,width,vv,profile)
!
!  Generate vector potential for a flux ring of magnetic flux FRING,
!  current Iring (not correctly normalized), radius R0 and thickness
!  WIDTH in normal orientation (lying in the x-y plane, centred at (0,0,0)).
!
!   1-may-02/wolf: coded
!
      use Cdata, only: mx,my,mz
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,3) :: vv
      real, dimension (mx,my,mz)   :: xx,yy,zz,phi,tmp
      real :: fring,Iring,R0,width
      character (len=*) :: profile
!
      vv = 0.
!
!  magnetic ring
!
      tmp = sqrt(xx**2+yy**2)-R0

      select case(profile)

      case('tanh')
        vv(:,:,:,3) = - fring * 0.5*(1+tanh(tmp/width)) &
                              * 0.5/width/cosh(zz/width)**2

      case default
        call stop_it('norm_ring: No such fluxtube profile')
      endselect
!
!  current ring (to twist the B-lines)
!
!      tmp = tmp**2 + zz**2 + width**2  ! need periodic analog of this
      tmp = width - sqrt(tmp**2 + zz**2)
      tmp = Iring*0.5*(1+tanh(tmp/width))     ! Now the A_phi component
      phi = atan2(yy,xx)
      vv(:,:,:,1) = - tmp*sin(phi)
      vv(:,:,:,2) =   tmp*cos(phi)
!
    endsubroutine norm_ring
!***********************************************************************
    subroutine force_free_jet(mu,xx,yy,zz)
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
      use Cdata, only: x,y,z,lroot,directory,ip,m,n,pi,r_ref
      use Sub, only: hypergeometric2F1,gamma_function
      use Global, only: set_global
      use Deriv, only: der
      use IO, only: output

      real, intent(in) :: mu
      real, dimension(mx,my,mz), intent(in) :: xx,yy,zz
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

      B1r_=sin(pi*mu/2)*gamma_function(   abs(mu) /2) / &
                        gamma_function((1+abs(mu))/2)

      B1z_=cos(pi*mu/2)*gamma_function((1+abs(mu))/2) / &
                        gamma_function((2+abs(mu))/2)

      B1=sqrt(4/pi)*r_ref**(mu-1)*sqrt(B1r_**2+B1z_**2)
!
!  calculate external vector potential
!
      if (lroot) print*,'FORCE_FREE_JET: calculating external vector potential'

      if (lforce_free_test) then

        if (lroot) print*,'FORCE_FREE_JET: using analytic solution for mu=-1'
        Ax_ext=-2*yy*(1-zz/sqrt(xx**2+yy**2+zz**2))/(xx**2+yy**2)/B1
        Ay_ext= 2*xx*(1-zz/sqrt(xx**2+yy**2+zz**2))/(xx**2+yy**2)/B1

      else

        do l=1,mx
        do m=1,my
        do n=1,mz

          r2=x(l)**2+y(m)**2
          xi2=r2/(r2+z(n)**2)
          A_phi=hypergeometric2F1((1-mu)/2,(2+mu)/2,2.0,xi2,tol) &
               *sqrt(xi2)*sqrt(r2+z(n)**2)**mu/B1

          Ax_ext(l,m,n)=-y(m)*A_phi/sqrt(r2)
          Ay_ext(l,m,n)= x(l)*A_phi/sqrt(r2)

        enddo
        enddo
        enddo

      endif

!
!  calculate external magnetic field
!
      if (lroot.and.ip<=5) &
        print*,'FORCE_FREE_JET: calculating the external magnetic field'

      do n=n1,n2
      do m=m1,m2
        call der(Ay_ext,bb_x,3)
        bb_ext_pot(:,1)=-bb_x
        call der(Ax_ext,bb_y,3)
        bb_ext_pot(:,2)= bb_y
        call der(Ay_ext,bb_z,1)
        bb_ext_pot(:,3)= bb_z
        call der(Ax_ext,bb_z,2)
        bb_ext_pot(:,3)=bb_ext_pot(:,3)-bb_z
        call set_global(bb_ext_pot,m,n,'B_ext_pot',nx)
      enddo
      enddo

      if (ip<=5) then
        call output(trim(directory)//'/Ax_ext.dat',Ax_ext,1)
        call output(trim(directory)//'/Ay_ext.dat',Ay_ext,1)
      endif

    endsubroutine force_free_jet
!***********************************************************************
    subroutine piecew_dipole_aa(ampl,inclaa,f,ivar,xx,yy,zz)
!
!  A field that is vertical uniform for r<R_int, inclined dipolar for
!  r>R_ext, and potential in the shell R_int<r<R_ext.
!  This mimics a neutron star just after the Meissner effect forced the
!  internal field to become vertical (aligned with rotation axis).
!
!  AMPL represents mu/4 pi, where  mu = 1/2 Int rr × jj dV  is the
!  magnetic moment of the external dipole field.
!  INCLAA is the inclination of the dipolar field.
!
!  Pencilized in order to minimize memory consumption with all the
!  auxiliary variables used.
!
!  23-jul-05/wolf:coded
!
      use Cdata, only: pi,tini, &
                       imn,m,n,mm,nn, &
                       r_int,r_ext
!
      real, intent(inout), dimension (mx,my,mz,mfarray) :: f
      real, intent(in), dimension (mx,my,mz) :: xx,yy,zz
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
      if (NO_WARN) print*, xx(1,1,1),yy(1,1,1),zz(1,1,1) !(keep compiler quiet)
!
    endsubroutine piecew_dipole_aa
!***********************************************************************
    subroutine geo_benchmark_B(f)
!
!  30-june-04/grs: coded
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx) :: theta_mn,ar,atheta,aphi,r_mn,phi_mn
      real :: C_int,C_ext,A_int,A_ext
      integer :: j

      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
        theta_mn=acos(spread(z(n),1,nx)/r_mn)
        phi_mn=atan2(spread(y(m),1,nx),x(l1:l2))

 ! calculate ax,ay,az (via ar,atheta,aphi) inside shell (& leave zero outside shell)

        do j=1,ninit
           select case(initaa(j))
           case('geo-benchmark-case1')
              if (lroot .and. imn==1) print*, 'geo_benchmark_B: geo-benchmark-case1'
              C_int=-( -1./63.*r_int**4 + 11./84.*r_int**3*r_ext            &
                     + 317./1050.*r_int**2*r_ext**2                         &
                     - 1./5.*r_int**2*r_ext**2*log(r_int) )
              C_ext=-( -1./63.*r_ext**9 + 11./84.*r_ext**8*r_int            &
                     + 317./1050.*r_ext**7*r_int**2                         &
                     - 1./5.*r_ext**7*r_int**2*log(r_ext) )
              A_int=5./2.*(r_ext-r_int)
              A_ext=5./8.*(r_ext**4-r_int**4)

              where (r_mn < r_int)
                ar=C_int*ampl_B0*80.*2.*(3.*sin(theta_mn)**2-2.)*r_mn
                atheta=3.*C_int*ampl_B0*80.*sin(2.*theta_mn)*r_mn
                aphi=ampl_B0*A_int*r_mn*sin(theta_mn)
              endwhere

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

              where (r_mn > r_ext)
                ar=C_ext*ampl_B0*80.*2.*(3.*sin(theta_mn)**2-2.)/r_mn**4
                atheta=-2.*C_ext*ampl_B0*80.*sin(2.*theta_mn)/r_mn**4
                aphi=ampl_B0*A_ext/r_mn**2*sin(theta_mn)
              endwhere

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

           case('geo-benchmark-case2')
              if (lroot .and. imn==1) print*, 'geo_benchmark_B: geo-benchmark-case2 not yet coded.'

           case default
              if (lroot .and. imn==1) print*,'geo_benchmark_B: case not defined!'
              call stop_it("")
           endselect
        enddo
        f(l1:l2,m,n,iax)=sin(theta_mn)*cos(phi_mn)*ar + cos(theta_mn)*cos(phi_mn)*atheta - sin(phi_mn)*aphi
        f(l1:l2,m,n,iay)=sin(theta_mn)*sin(phi_mn)*ar + cos(theta_mn)*sin(phi_mn)*atheta + cos(phi_mn)*aphi
        f(l1:l2,m,n,iaz)=cos(theta_mn)*ar - sin(theta_mn)*atheta
     enddo


     if (ip<=14) then
        print*,'geo_benchmark_B: minmax(ax) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iax)),maxval(f(l1:l2,m1:m2,n1:n2,iax))
        print*,'geo_benchmark_B: minmax(ay) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iay)),maxval(f(l1:l2,m1:m2,n1:n2,iay))
        print*,'geo_benchmark_B: minmax(az) on iproc:', iproc, minval(f(l1:l2,m1:m2,n1:n2,iaz)),maxval(f(l1:l2,m1:m2,n1:n2,iaz))
     endif

    endsubroutine geo_benchmark_B

!***********************************************************************
    subroutine bc_frozen_in_bb(topbot,j)
!
!  Set flags to indicate that magnetic flux is frozen-in at the
!  z boundary. The implementation occurs in daa_dt where magnetic
!  diffusion is switched off in that layer.
!
      use Cdata
!
      character (len=3) :: topbot
      integer :: j
!
      select case(topbot)
      case('bot')               ! bottom boundary
        lfrozen_bb_bot(j-iax+1) = .true.    ! set flag
      case('top')               ! top boundary
        lfrozen_bb_top(j-iax+1) = .true.    ! set flag
      case default
        print*, "bc_frozen_in_bb: ", topbot, " should be `top' or `bot'"
      endselect
!
    endsubroutine bc_frozen_in_bb
!***********************************************************************
    subroutine bc_aa_pot3(f,topbot)
!
!  Pontential field boundary condition
!
!  11-oct-06/wolf: Adapted from Tobi's bc_aa_pot2
!
      use Fourier, only: fourier_transform_xy_xy
      use Mpicomm, only: communicate_bc_aa_pot

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot

      real, dimension (nx,ny,iax:iaz) :: aa_re,aa_im
      real, dimension (nx,ny) :: kx,ky,kappa,kappa1,exp_fact
      real, dimension (nx,ny) :: tmp_re,tmp_im
      real, dimension (nx,ny) :: fac
      real    :: delta_z
      integer :: i,j
!
!  Get local wave numbers
!
      kx = spread(kx_fft(ipx*nx+1:ipx*nx+nx),2,ny)
      ky = spread(ky_fft(ipy*ny+1:ipy*ny+ny),1,nx)
!
!  Calculate 1/k^2, zero mean
!
      kappa = sqrt(kx**2 + ky**2)
      where (kappa > 0)
        kappa1 = 1/kappa
      elsewhere
        kappa1 = 0
      endwhere
!
!  Check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  Potential field condition at the bottom
!
      case('bot')

        do j=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          delta_z  = z(n1+j) - z(n1-j)
          exp_fact = exp(-kappa*delta_z)
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
          do i=iax,iaz
            tmp_re = f(l1:l2,m1:m2,n1+j,i)
            tmp_im = 0.0
            call fourier_transform_xy_xy(tmp_re,tmp_im)
            aa_re(:,:,i) = tmp_re*exp_fact
            aa_im(:,:,i) = tmp_im*exp_fact
          enddo

         ! Transform back
          do i=iax,iaz
            tmp_re = aa_re(:,:,i)
            tmp_im = aa_im(:,:,i)
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            f(l1:l2,m1:m2,n1-j,i) = tmp_re
          enddo

        enddo
!
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
        call communicate_bc_aa_pot(f,topbot)
!
!  Potential field condition at the top
!
      case('top')

        do j=1,nghost
!
! Calculate delta_z based on z(), not on dz to improve behavior for
! non-equidistant grid (still not really correct, but could be OK)
!
          delta_z  = z(n2+j) - z(n2-j)
          exp_fact = exp(-kappa*delta_z)
!
!  Determine potential field in ghost zones
!
          !  Fourier transforms of x- and y-components on the boundary
          do i=iax,iaz
            tmp_re = f(l1:l2,m1:m2,n2-j,i)
            tmp_im = 0.0
            call fourier_transform_xy_xy(tmp_re,tmp_im)
            aa_re(:,:,i) = tmp_re*exp_fact
            aa_im(:,:,i) = tmp_im*exp_fact
          enddo

          ! Transform back
          do i=iax,iaz
            tmp_re = aa_re(:,:,i)
            tmp_im = aa_im(:,:,i)
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            f(l1:l2,m1:m2,n2+j,i) = tmp_re
          enddo

        enddo
!
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
        call communicate_bc_aa_pot(f,topbot)

      case default

        if (lroot) print*,"bc_aa_pot2: invalid argument"

      endselect

    endsubroutine bc_aa_pot3
!***********************************************************************
    subroutine bc_aa_pot2(f,topbot)
!
!  Pontential field boundary condition
!
!  10-oct-06/tobi: Coded
!
      use Fourier, only: fourier_transform_xy_xy
      use Mpicomm, only: communicate_bc_aa_pot

      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      character (len=3), intent (in) :: topbot

      real, dimension (nx,ny,iax:iaz) :: aa_re,aa_im
      real, dimension (nx,ny) :: kx,ky,kappa,kappa1
      real, dimension (nx,ny) :: tmp_re,tmp_im
      real, dimension (nx,ny) :: fac
      integer :: i,j
!
!  Get local wave numbers
!
      kx = spread(kx_fft(ipx*nx+1:ipx*nx+nx),2,ny)
      ky = spread(ky_fft(ipy*ny+1:ipy*ny+ny),1,nx)
!
!  Calculate 1/k^2, zero mean
!
      kappa = sqrt(kx**2 + ky**2)
      where (kappa > 0)
        kappa1 = 1/kappa
      elsewhere
        kappa1 = 0
      endwhere
!
!  Check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  Pontential field condition at the bottom
!
      case('bot')
!
!  Fourier transforms of x- and y-components on the boundary
!
        do i=iax,iaz
          tmp_re = f(l1:l2,m1:m2,n1,i)
          tmp_im = 0.0
          call fourier_transform_xy_xy(tmp_re,tmp_im)
          aa_re(:,:,i) = tmp_re
          aa_im(:,:,i) = tmp_im
        enddo
!
!  Determine potential field in ghost zones
!
        do j=1,nghost
          fac = exp(-j*kappa*dz)
          do i=iax,iaz
            tmp_re = fac*aa_re(:,:,i)
            tmp_im = fac*aa_im(:,:,i)
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            f(l1:l2,m1:m2,n1-j,i) = tmp_re
          enddo
        enddo
!
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
        call communicate_bc_aa_pot(f,topbot)
!
!  Pontential field condition at the top
!
      case('top')
!
!  Fourier transforms of x- and y-components on the boundary
!
        do i=iax,iaz
          tmp_re = f(l1:l2,m1:m2,n2,i)
          tmp_im = 0.0
          call fourier_transform_xy_xy(tmp_re,tmp_im)
          aa_re(:,:,i) = tmp_re
          aa_im(:,:,i) = tmp_im
        enddo
!
!  Determine potential field in ghost zones
!
        do j=1,nghost
          fac = exp(-j*kappa*dz)
          do i=iax,iaz
            tmp_re = fac*aa_re(:,:,i)
            tmp_im = fac*aa_im(:,:,i)
            call fourier_transform_xy_xy(tmp_re,tmp_im,linv=.true.)
            f(l1:l2,m1:m2,n2+j,i) = tmp_re
          enddo
        enddo
!
!  The vector potential needs to be known outside of (l1:l2,m1:m2) as well
!
        call communicate_bc_aa_pot(f,topbot)

      case default

        if (lroot) print*,"bc_aa_pot2: invalid argument"

      endselect

    endsubroutine bc_aa_pot2
!***********************************************************************
      subroutine bc_aa_pot(f,topbot)
!
!  Potential field boundary condition for magnetic vector potential at
!  bottom or top boundary (in z).
!
!  14-jun-2002/axel: adapted from similar
!   8-jul-2002/axel: introduced topbot argument
!
      use Cdata
      use Mpicomm, only: stop_it,communicate_bc_aa_pot
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,ny) :: f2,f3
      real, dimension (nx,ny,nghost+1) :: fz
      integer :: j
!
!  pontential field condition
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case(topbot)
!
!  pontential field condition at the bottom
!
      case('bot')
        if (headtt) print*,'bc_aa_pot: pot-field bdry cond at bottom'
        if (mod(nxgrid,nygrid)/=0) &
             call stop_it("bc_aa_pot: pot-field doesn't work "//&
                          "with mod(nxgrid,nygrid)/=1")
        do j=0,1
          f2=f(l1:l2,m1:m2,n1+1,iax+j)
          f3=f(l1:l2,m1:m2,n1+2,iax+j)
          call potential_field(fz,f2,f3,-1)
          f(l1:l2,m1:m2,1:n1,iax+j)=fz
        enddo
        !
        f2=f(l1:l2,m1:m2,n1,iax)
        f3=f(l1:l2,m1:m2,n1,iay)
        call potentdiv(fz,f2,f3,-1)
        f(l1:l2,m1:m2,1:n1,iaz)=-fz
        call communicate_bc_aa_pot(f,topbot)
!
!  pontential field condition at the top
!
      case('top')
        if (headtt) print*,'bc_aa_pot: pot-field bdry cond at top'
        if (mod(nxgrid,nygrid)/=0) &
             call stop_it("bc_aa_pot: pot-field doesn't work "//&
                          "with mod(nxgrid,nygrid)/=1")
        do j=0,1
          f2=f(l1:l2,m1:m2,n2-1,iax+j)
          f3=f(l1:l2,m1:m2,n2-2,iax+j)
          call potential_field(fz,f2,f3,+1)
          f(l1:l2,m1:m2,n2:mz,iax+j)=fz
        enddo
        !
        f2=f(l1:l2,m1:m2,n2,iax)
        f3=f(l1:l2,m1:m2,n2,iay)
        call potentdiv(fz,f2,f3,+1)
        f(l1:l2,m1:m2,n2:mz,iaz)=-fz
        call communicate_bc_aa_pot(f,topbot)
      case default
        if (lroot) print*,"bc_aa_pot: invalid argument"
      endselect
!
      endsubroutine bc_aa_pot
!***********************************************************************
      subroutine potential_field(fz,f2,f3,irev)
!
!  solves the potential field boundary condition;
!  fz is the boundary layer, and f2 and f3 are the next layers inwards.
!  The condition is the same on the two sides.
!
!  20-jan-00/axel+wolf: coded
!  22-mar-00/axel: corrected sign (it is the same on both sides)
!  29-sep-06/axel: removed multiple calls, removed normalization, non-para
!
      use Cdata
      use Fourier
!
      real, dimension (nx,ny) :: fac,kk,f1r,f1i,g1r,g1i,f2,f2r,f2i,f3,f3r,f3i
      real, dimension (nx,ny,nghost+1) :: fz
      real, dimension (nx) :: kx
      real, dimension (nygrid) :: ky
      real :: delz
      integer :: i,irev
!
!  initialize workspace
!
      f2r=f2; f2i=0
      f3r=f3; f3i=0
!
!  Transform; real and imaginary parts
!
      call fourier_transform_xy_xy(f2r,f2i)
      call fourier_transform_xy_xy(f3r,f3i)
!
!  define wave vector
!
      kx=cshift((/(i-(nx-1)/2,i=0,nx-1)/),+(nx-1)/2)*2*pi/Lx
      ky=cshift((/(i-(nygrid-1)/2,i=0,nygrid-1)/),+(nygrid-1)/2)*2*pi/Ly
!
!  calculate 1/k^2, zero mean
!
      kk=sqrt(spread(kx**2,2,ny)+spread(ky(ipy*ny+1:(ipy+1)*ny)**2,1,nx))
!
!  one-sided derivative
!
      fac=1./(3.+2.*dz*kk)
      f1r=fac*(4.*f2r-f3r)
      f1i=fac*(4.*f2i-f3i)
!
!  set ghost zones
!
      do i=0,nghost
        delz=i*dz
        fac=exp(-kk*delz)
        g1r=fac*f1r
        g1i=fac*f1i
!
!  Transform back
!
        call fourier_transform_xy_xy(g1r,g1i,linv=.true.)
!
!  reverse order if irev=-1 (if we are at the bottom)
!
        if (irev==+1) fz(:,:,       i+1) = g1r
        if (irev==-1) fz(:,:,nghost-i+1) = g1r
      enddo
!
    endsubroutine potential_field
!***********************************************************************
      subroutine potentdiv(fz,f2,f3,irev)
!
!  solves the divA=0 for potential field boundary condition;
!  f2 and f3 correspond to Ax and Ay (input) and fz corresponds to Ax (out)
!  In principle we could save some ffts, by combining with the potential
!  subroutine above, but this is now easier
!
!  22-mar-02/axel: coded
!  29-sep-06/axel: removed multiple calls, removed normalization, non-para
!   7-oct-06/axel: corrected sign for irev==+1.
!
      use Cdata
      use Fourier
!
      real, dimension (nx,ny) :: fac,kk,kkkx,kkky,f1r,f1i,g1r,g1i,f2,f2r,f2i,f3,f3r,f3i
      real, dimension (nx,ny,nghost+1) :: fz
      real, dimension (nx) :: kx
      real, dimension (nygrid) :: ky
      real :: delz
      integer :: i,irev
!
      f2r=f2; f2i=0
      f3r=f3; f3i=0
!
!  Transform
!
      call fourier_transform_xy_xy(f2r,f2i)
      call fourier_transform_xy_xy(f3r,f3i)
!
!  define wave vector
!
      kx=cshift((/(i-nx/2,i=0,nx-1)/),+nx/2)*2*pi/Lx
      ky=cshift((/(i-nygrid/2,i=0,nygrid-1)/),+nygrid/2)*2*pi/Ly
!
!  calculate 1/k^2, zero mean
!
      kk=sqrt(spread(kx**2,2,ny)+spread(ky(ipy*ny+1:(ipy+1)*ny)**2,1,nx))
      kkkx=spread(kx,2,ny)
      kkky=spread(ky(ipy*ny+1:(ipy+1)*ny),1,nx)
!
!  calculate 1/kk
!
      kk(1,1)=1.
      fac=1./kk
      fac(1,1)=0.
!
      f1r=fac*(-kkkx*f2i-kkky*f3i)
      f1i=fac*(+kkkx*f2r+kkky*f3r)
!
!  set ghost zones
!
      do i=0,nghost
        delz=i*dz
        fac=exp(-kk*delz)
        g1r=fac*f1r
        g1i=fac*f1i
!
!  Transform back
!
        call fourier_transform_xy_xy(g1r,g1i,linv=.true.)
!
!  reverse order if irev=-1 (if we are at the bottom)
!  but reverse sign if irev=+1 (if we are at the top)
!
        if (irev==+1) fz(:,:,       i+1) = -g1r
        if (irev==-1) fz(:,:,nghost-i+1) = +g1r
      enddo
!
    endsubroutine potentdiv
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
        case('fs')
          z2 = z**2.
!  resistivity profile from Fleming & Stone (ApJ 585:908-920)
          eta_z = eta*exp(-z2/2.+sigma_ratio*erfcc(abs(z))/4.)
!
! its gradient: 
          geta_z(:,1) = 0.
          geta_z(:,2) = 0.
          geta_z(:,3) = eta_z*(-z-sign(1.,z)*sigma_ratio*exp(-z2)/(2.*sqrt(pi)))
!
        case('tanh')
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
      use Cdata

      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,3), intent (out) :: bb_hat

      !Tobi: Not sure about this value
      real, parameter :: tol=1e-11

      real, dimension (mx,3) :: bb,bb2
      real, dimension (mx) :: bb_len,aerr2
      real :: fac
      integer :: j

!
!  Compute magnetic field from vector potential.
!
      bb=0.

      if (nxgrid/=1) then
        fac = 1/(2*dx)
        bb(l1-2:l2+2,3) = bb(l1-2:l2+2,3) + fac*( f(l1-1:l2+3,m  ,n  ,iay)   &
                                                - f(l1-3:l2+1,m  ,n  ,iay) )
        bb(l1-2:l2+2,2) = bb(l1-2:l2+2,2) - fac*( f(l1-1:l2+3,m  ,n  ,iaz)   &
                                                - f(l1-3:l2+1,m  ,n  ,iaz) )
      endif

      if (nygrid/=1) then
        fac = 1/(2*dy)
        bb(l1-2:l2+2,1) = bb(l1-2:l2+2,1) + fac*( f(l1-2:l2+2,m+1,n  ,iaz)   &
                                                - f(l1-2:l2+2,m-1,n  ,iaz) )
        bb(l1-2:l2+2,3) = bb(l1-2:l2+2,3) - fac*( f(l1-2:l2+2,m+1,n  ,iax)   &
                                                - f(l1-2:l2+2,m-1,n  ,iax) )
      endif

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

      aerr2 = tol**2 * max(sum(bb2,2),1.)

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

      do j=1,3; bb_hat(:,j) = bb_hat(:,j)/(bb_len+tini); enddo

    endsubroutine bb_unitvec_shock
!***********************************************************************
endmodule Magnetic
