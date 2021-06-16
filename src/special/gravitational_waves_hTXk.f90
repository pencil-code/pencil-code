! $Id$
!
!  This module provide a way for users to specify custom
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc.
!
!  The module provides a set of standard hooks into the Pencil-Code and
!  currently allows the following customizations:
!
!  Description                                     | Relevant function call
!  ---------------------------------------------------------------------------
!  Special variable registration                   | register_special
!    (pre parameter read)                          |
!  Special variable initialization                 | initialize_special
!    (post parameter read)                         |
!  Special variable finalization                   | finalize_special
!    (deallocation, etc.)                          |
!                                                  |
!  Special initial condition                       | init_special
!   this is called last so may be used to modify   |
!   the mvar variables declared by this module     |
!   or optionally modify any of the other f array  |
!   variables.  The latter, however, should be     |
!   avoided where ever possible.                   |
!                                                  |
!  Special term in the mass (density) equation     | special_calc_density
!  Special term in the momentum (hydro) equation   | special_calc_hydro
!  Special term in the energy equation             | special_calc_energy
!  Special term in the induction (magnetic)        | special_calc_magnetic
!     equation                                     |
!                                                  |
!  Special equation                                | dspecial_dt
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MAUX CONTRIBUTION 18
!
! PENCILS PROVIDED stress_ij(6)
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!   lspecial = .true.
! to enable use of special hooks.
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code SVN repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module
! selections to say something like:
!
!   SPECIAL=special/geo_kws
!
! Where geo_kws it replaced by the filename of your new module
! upto and not including the .f90
!
module Special
!
  use Cparam
  use Cdata
  use Initcond
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error, warning
!
  implicit none
!
  include '../special.h'
!
! Declare index of new variables in f array (if any).
!
  character (len=labellen) :: inithij='nothing'
  character (len=labellen) :: initgij='nothing'
  character (len=labellen) :: ctrace_factor='1/3'
  character (len=labellen) :: cstress_prefactor='6'
  character (len=labellen) :: fourthird_in_stress='4/3'
  character (len=labellen) :: cc_light='1'
  character (len=labellen) :: aux_stress='stress'
  real :: amplhij=0., amplgij=0., dummy=0.
  real :: trace_factor=0., stress_prefactor, fourthird_factor, EGWpref
  real :: nscale_factor_conformal=1., tshift=0.
  logical :: lno_transverse_part=.false., lgamma_factor=.false.
  logical :: lswitch_sign_e_X=.true., lswitch_symmetric=.false., ldebug_print=.false.
  logical :: lStress_as_aux=.true., lreynolds=.false., lkinGW=.true.
  logical :: lelectmag=.false.
  logical :: lggTX_as_aux=.true., lhhTX_as_aux=.true.
  logical :: lremove_mean_hij=.false., lremove_mean_gij=.false.
  logical :: GWs_spec_complex=.true. !(fixed for now)
  logical :: lreal_space_hTX_as_aux=.false., lreal_space_gTX_as_aux=.false.
  logical :: linflation=.false., lreheating_GW=.false.
  logical :: lonly_mag=.false.
  logical :: lstress=.true., lstress_ramp=.false., ldelkt=.false.
  logical :: lnonlinear_source=.false., lnonlinear_Tpq_trans=.true.
  real, dimension(3,3) :: ij_table
  real :: c_light2=1., delk=0., tdelk=0., tstress_ramp=0.
!
  real, dimension (:,:,:,:), allocatable :: Tpq_re, Tpq_im
  real, dimension (:,:,:,:), allocatable :: nonlinear_Tpq_re, nonlinear_Tpq_im
  real :: kscale_factor, tau_stress_comp=0., exp_stress_comp=0.
  real :: tau_stress_kick=0., tnext_stress_kick=1., fac_stress_kick=2., accum_stress_kick=1.
  real :: nonlinear_source_fact=0.
!
! input parameters
  namelist /special_init_pars/ &
    ctrace_factor, cstress_prefactor, fourthird_in_stress, lno_transverse_part, &
    inithij, initgij, amplhij, amplgij, lStress_as_aux, lgamma_factor, &
    lreal_space_hTX_as_aux, lreal_space_gTX_as_aux, &
    lelectmag, lggTX_as_aux, lhhTX_as_aux, linflation, lreheating_GW, lonly_mag
!
! run parameters
  namelist /special_run_pars/ &
    ctrace_factor, cstress_prefactor, fourthird_in_stress, lno_transverse_part, &
    ldebug_print, lswitch_sign_e_X, lswitch_symmetric, lStress_as_aux, &
    nscale_factor_conformal, tshift, cc_light, lgamma_factor, &
    lStress_as_aux, lkinGW, aux_stress, tau_stress_comp, exp_stress_comp, &
    lelectmag, tau_stress_kick, fac_stress_kick, delk, tdelk, ldelkt, &
    lreal_space_hTX_as_aux, lreal_space_gTX_as_aux, &
    lggTX_as_aux, lhhTX_as_aux, lremove_mean_hij, lremove_mean_gij, &
    lstress, lstress_ramp, tstress_ramp, linflation, lreheating_GW, &
    lnonlinear_source, lnonlinear_Tpq_trans, nonlinear_source_fact
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_STrept=0      ! DIAG_DOC: $Re S_{T}(k_1,k_1,k_1,t)$
  integer :: idiag_STimpt=0      ! DIAG_DOC: $Im S_{T}(k_1,k_1,k_1,t)$
  integer :: idiag_SXrept=0      ! DIAG_DOC: $Re S_{X}(k_1,k_1,k_1,t)$
  integer :: idiag_SXimpt=0      ! DIAG_DOC: $Im S_{X}(k_1,k_1,k_1,t)$
  integer :: idiag_hTrept=0      ! DIAG_DOC: $Re h_{T}(k_1,k_1,k_1,t)$
  integer :: idiag_hTimpt=0      ! DIAG_DOC: $Im h_{T}(k_1,k_1,k_1,t)$
  integer :: idiag_hXrept=0      ! DIAG_DOC: $Re h_{X}(k_1,k_1,k_1,t)$
  integer :: idiag_hXimpt=0      ! DIAG_DOC: $Im h_{X}(k_1,k_1,k_1,t)$
  integer :: idiag_gTrept=0      ! DIAG_DOC: $Re h_{T}(k_1,k_1,k_1,t)$
  integer :: idiag_gTimpt=0      ! DIAG_DOC: $Im h_{T}(k_1,k_1,k_1,t)$
  integer :: idiag_gXrept=0      ! DIAG_DOC: $Re h_{X}(k_1,k_1,k_1,t)$
  integer :: idiag_gXimpt=0      ! DIAG_DOC: $Im h_{X}(k_1,k_1,k_1,t)$
  integer :: idiag_STrep2=0      ! DIAG_DOC: $Re S_{T}(k_2,k_2,k_2,t)$
  integer :: idiag_STimp2=0      ! DIAG_DOC: $Im S_{T}(k_2,k_2,k_2,t)$
  integer :: idiag_SXrep2=0      ! DIAG_DOC: $Re S_{X}(k_2,k_2,k_2,t)$
  integer :: idiag_SXimp2=0      ! DIAG_DOC: $Im S_{X}(k_2,k_2,k_2,t)$
  integer :: idiag_hTrep2=0      ! DIAG_DOC: $Re h_{T}(k_2,k_2,k_2,t)$
  integer :: idiag_hTimp2=0      ! DIAG_DOC: $Im h_{T}(k_2,k_2,k_2,t)$
  integer :: idiag_hXrep2=0      ! DIAG_DOC: $Re h_{X}(k_2,k_2,k_2,t)$
  integer :: idiag_hXimp2=0      ! DIAG_DOC: $Im h_{X}(k_2,k_2,k_2,t)$
  integer :: idiag_gTrep2=0      ! DIAG_DOC: $Re g_{T}(k_2,k_2,k_2,t)$
  integer :: idiag_gTimp2=0      ! DIAG_DOC: $Im g_{T}(k_2,k_2,k_2,t)$
  integer :: idiag_gXrep2=0      ! DIAG_DOC: $Re g_{X}(k_2,k_2,k_2,t)$
  integer :: idiag_gXimp2=0      ! DIAG_DOC: $Im g_{X}(k_2,k_2,k_2,t)$
  integer :: idiag_g11pt=0       ! DIAG_DOC: $g_{11}(x_1,y_1,z_1,t)$
  integer :: idiag_g22pt=0       ! DIAG_DOC: $g_{22}(x_1,y_1,z_1,t)$
  integer :: idiag_g33pt=0       ! DIAG_DOC: $g_{33}(x_1,y_1,z_1,t)$
  integer :: idiag_g12pt=0       ! DIAG_DOC: $g_{12}(x_1,y_1,z_1,t)$
  integer :: idiag_g23pt=0       ! DIAG_DOC: $g_{23}(x_1,y_1,z_1,t)$
  integer :: idiag_g31pt=0       ! DIAG_DOC: $g_{31}(x_1,y_1,z_1,t)$
  integer :: idiag_hhTpt=0       ! DIAG_DOC: $h_{T}(x_1,y_1,z_1,t)$
  integer :: idiag_hhXpt=0       ! DIAG_DOC: $h_{X}(x_1,y_1,z_1,t)$
  integer :: idiag_ggTpt=0       ! DIAG_DOC: $\dot{h}_{T}(x_1,y_1,z_1,t)$
  integer :: idiag_ggXpt=0       ! DIAG_DOC: $\dot{h}_{X}(x_1,y_1,z_1,t)$
  integer :: idiag_hhTp2=0       ! DIAG_DOC: $h_{T}(x_1,y_1,z_1,t)$
  integer :: idiag_hhXp2=0       ! DIAG_DOC: $h_{X}(x_1,y_1,z_1,t)$
  integer :: idiag_ggTp2=0       ! DIAG_DOC: $\dot{h}_{T}(x_1,y_1,z_1,t)$
  integer :: idiag_ggXp2=0       ! DIAG_DOC: $\dot{h}_{X}(x_1,y_1,z_1,t)$
  integer :: idiag_hrms=0        ! DIAG_DOC: $\bra{h_T^2+h_X^2}^{1/2}$
  integer :: idiag_EEGW=0        ! DIAG_DOC: $\bra{g_T^2+g_X^2}\,c^2/(32\pi G)$
  integer :: idiag_gg2m=0        ! DIAG_DOC: $\bra{g_T^2+g_X^2}$
  integer :: idiag_Stgm=0        ! DIAG_DOC: $\bra{S_Tg_T+S_Xg_X}$
  integer :: idiag_hhT2m=0       ! DIAG_DOC: $\bra{h_T^2}$
  integer :: idiag_hhX2m=0       ! DIAG_DOC: $\bra{h_X^2}$
  integer :: idiag_hhTXm=0       ! DIAG_DOC: $\bra{h_T h_X}$
  integer :: idiag_ggT2m=0       ! DIAG_DOC: $\bra{g_T^2}$
  integer :: idiag_ggX2m=0       ! DIAG_DOC: $\bra{g_X^2}$
  integer :: idiag_ggTXm=0       ! DIAG_DOC: $\bra{g_T g_X}$
!
  integer :: ihhT_realspace, ihhX_realspace
  integer :: iggT_realspace, iggX_realspace
  integer, parameter :: nk=nxgrid/2
  type, public :: GWspectra
    real, dimension(nk) :: GWs   ,GWh   ,GWm   ,Str   ,Stg
    real, dimension(nk) :: GWshel,GWhhel,GWmhel,Strhel,Stghel
    real, dimension(nk) :: SCL, VCT, Tpq, TGW
    complex, dimension(nx) :: complex_Str_T, complex_Str_X
  endtype GWspectra

  type(GWspectra) :: spectra

  contains
!***********************************************************************
    subroutine register_special
!
!  Set up indices for variables in special modules.
!
!   3-aug-17/axel: coded
!  28-mar-21/axel:allowed for lreal_space_hTX_as_aux
!
      use Sub, only: register_report_aux
      use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Register ggT and ggX as auxiliary arrays
!  May want to do this only when Fourier transform is enabled.
!
      if (lggTX_as_aux) then
        call farray_register_auxiliary('ggT',iggT)
        call farray_register_auxiliary('ggX',iggX)
        call farray_register_auxiliary('ggTim',iggTim)
        call farray_register_auxiliary('ggXim',iggXim)
      endif
!
      if (lhhTX_as_aux) then
        call farray_register_auxiliary('hhT',ihhT)
        call farray_register_auxiliary('hhX',ihhX)
        call farray_register_auxiliary('hhTim',ihhTim)
        call farray_register_auxiliary('hhXim',ihhXim)
      endif
!
      if (lStress_as_aux) then
        call farray_register_auxiliary('StT',iStressT)
        call farray_register_auxiliary('StX',iStressX)
        call farray_register_auxiliary('StTim',iStressTim)
        call farray_register_auxiliary('StXim',iStressXim)
        call farray_register_auxiliary('Str',iStress_ij,array=6)
      endif
!
!  To get hT and hX in real space, invoke lreal_space_hTX_as_aux
!
      if (lreal_space_hTX_as_aux) then
        call farray_register_auxiliary('hhT_realspace',ihhT_realspace)
        call farray_register_auxiliary('hhX_realspace',ihhX_realspace)
      endif
!
!  To get hT and hX in real space, invoke lreal_space_gTX_as_aux
!
      if (lreal_space_gTX_as_aux) then
        call farray_register_auxiliary('ggT_realspace',iggT_realspace)
        call farray_register_auxiliary('ggX_realspace',iggX_realspace)
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      use EquationOfState, only: cs0
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: stat, i
!
!  set index table
!
      ij_table(1,1)=1
      ij_table(2,2)=2
      ij_table(3,3)=3
      ij_table(1,2)=4
      ij_table(2,3)=5
      ij_table(3,1)=6
      ij_table(2,1)=4
      ij_table(3,2)=5
      ij_table(1,3)=6
!
!  determine trace factor
!
      select case (ctrace_factor)
        case ('0'); trace_factor=.0
        case ('1/2'); trace_factor=.5
        case ('1/3'); trace_factor=onethird
        case default
          call fatal_error("initialize_special: No such value for ctrace_factor:" &
              ,trim(ctrace_factor))
      endselect
!
!  Determine fourthird_in_stress. This factor is normally 4/3, but it can be
!  set to unity in case we want to pretend that the kinematic Beltrami field
!  has the same prefactor as a magnetic one.
!
      select case (fourthird_in_stress)
        case ('1'); fourthird_factor=1.
        case ('4/3'); fourthird_factor=fourthird
        case default
          call fatal_error("initialize_special: No such value for fourthird_in_stress:" &
              ,trim(fourthird_in_stress))
      endselect
!
!  determine stress_prefactor and GW energy prefactor,
!  which is EGWpref=.5*16.*pi/stress_prefactor**2
!  At the moment, only the case stress_prefactor=6 is to be used.
!  The other cases are kept for backward compatibility.
!
      select case (cstress_prefactor)
        case ('1'); stress_prefactor=1.; EGWpref=8.*pi
        case ('6'); stress_prefactor=6.; EGWpref=1./6.
        case ('24'); stress_prefactor=24.; EGWpref=1./6.
        case ('6old'); stress_prefactor=6.; EGWpref=1./(32.*pi)
        case ('16pi'); stress_prefactor=16.*pi; EGWpref=1./(32.*pi)
        case ('16pi_corr'); stress_prefactor=16.*pi; EGWpref=1./(16.*pi)
        case ('16piG/c^2'); stress_prefactor=16.*pi*G_Newton_cgs/c_light_cgs**2;
          EGWpref=c_light_cgs**2/(32.*pi*G_Newton_cgs)
        case default
          call fatal_error("initialize_special: No such value for ctrace_factor:" &
              ,trim(ctrace_factor))
      endselect
      if (headt) print*,'stress_prefactor=',stress_prefactor
      if (headt) print*,'EGWpref=',EGWpref
!
!  set speed of light
!
      select case (cc_light)
        case ('1'); c_light2=1.
        case ('cgs'); c_light2=c_light_cgs**2
        case default
          call fatal_error("initialize_special: No such value for cc_light:" &
              ,trim(ctrace_factor))
      endselect
      if (headt) print*,'c_light2=',c_light2
!
!  Determine whether Reynolds stress needs to be computed:
!  Compute Reynolds stress when lhydro or lhydro_kinematic,
!  provided lkinGW=T
!
      lreynolds=(lhydro.or.lhydro_kinematic).and.lkinGW
!
!  give a warning if cs0**2=1
!
!     if (cs0==1.) call fatal_error('gravitational_waves_hij6', &
!         'cs0 should probably not be unity')
!
      if (.not.allocated(Tpq_re)) allocate(Tpq_re(nx,ny,nz,6),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for Tpq_re')
!
      if (.not.allocated(Tpq_im)) allocate(Tpq_im(nx,ny,nz,6),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for Tpq_im')
!
!  Allocate memory for nonlinear source
!
      if (lnonlinear_source) then
        if (.not.allocated(nonlinear_Tpq_re)) allocate(nonlinear_Tpq_re(nx,ny,nz,6),stat=stat)
        if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for nonlinear_Tpq_re')
!
!  Need imaginary part only if lnonlinear_Tpq_trans=T
!
        if (lnonlinear_Tpq_trans) then
          if (.not.allocated(nonlinear_Tpq_im)) allocate(nonlinear_Tpq_im(nx,ny,nz,6),stat=stat)
          if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for nonlinear_Tpq_im')
        endif
      endif
!
!  calculate kscale_factor (for later binning)
!
      kscale_factor=2*pi/Lx
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine finalize_special(f)
!
!  Called right before exiting.
!
!  14-aug-2011/Bourdin.KIS: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine finalize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  initial condition for hij
!
      select case (inithij)
        case ('nothing')
          if (lroot) print*,'init_special: nothing'
          f(:,:,:,ihhT:ihhXim)=0.
        case default
          call fatal_error("init_special: No such value for inithij:" &
              ,trim(inithij))
      endselect
!
!  initial condition for gij
!
      select case (initgij)
        case ('nothing')
          if (lroot) print*,'init_special: nothing'
          f(:,:,:,iggT:iggXim)=0.
        case default
          call fatal_error("init_special: No such value for initgij:" &
              ,trim(initgij))
      endselect
!
      if (lStress_as_aux) then
        f(:,:,:,iStressT)=0.
        f(:,:,:,iStressX)=0.
        f(:,:,:,iStressTim)=0.
        f(:,:,:,iStressXim)=0.
      endif
!
      if (lreal_space_gTX_as_aux) then
        f(:,:,:,iggT_realspace)=0.
        f(:,:,:,iggX_realspace)=0.
      endif
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!   1-apr-06/axel: coded
!
!  Velocity field needed for Reynolds stress
!
      if (lreynolds) then
        lpenc_requested(i_uu)=.true.
        lpenc_requested(i_rho)=.true.
        if (trace_factor/=0..or.lgamma_factor) lpenc_requested(i_u2)=.true.
      endif
!
!  Magnetic field needed for Maxwell stress
!
      if (lmagnetic) then
        lpenc_requested(i_bb)=.true.
        if (trace_factor/=0.) lpenc_requested(i_b2)=.true.
      endif
!
!  Electric field needed for Maxwell stress
!
      if (lelectmag) then
        lpenc_requested(i_el)=.true.
        if (trace_factor/=0.) lpenc_requested(i_e2)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-jul-06/tony: coded
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-aug-17/axel: coded
!   7-jun-18/axel: included 4/3*rho factor
!
      use Deriv, only: derij
      use Sub, only: dot2_mn
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: prefactor, EEE2
      real, dimension (nx,3) :: EEEE
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
      integer :: i, j, ij
      real :: fact
!
!  Construct stress tensor; notice opposite signs for u and b.
!
      if (lgamma_factor) then
        prefactor=fourthird_factor/(1.-p%u2)
      else
        prefactor=fourthird_factor
      endif
!
!  Construct stress tensor; notice opposite signs for u and b.
!
      p%stress_ij=0.0
      if (lstress) then
        do j=1,3
        do i=1,j
          ij=ij_table(i,j)
          if (lonly_mag) then
            if (lmagnetic) p%stress_ij(:,ij)=p%stress_ij(:,ij)-p%bb(:,i)*p%bb(:,j)
          else
            if (lreynolds) p%stress_ij(:,ij)=p%stress_ij(:,ij)+p%uu(:,i)*p%uu(:,j)*prefactor*p%rho
            if (lmagnetic) p%stress_ij(:,ij)=p%stress_ij(:,ij)-p%bb(:,i)*p%bb(:,j)
            if (lelectmag) p%stress_ij(:,ij)=p%stress_ij(:,ij)-p%el(:,i)*p%el(:,j)
          endif
!
!  Remove trace.
!
          if (i==j) then
            if (lonly_mag) then  
              if (lmagnetic) p%stress_ij(:,ij)=p%stress_ij(:,ij)+trace_factor*p%b2
            else
              if (lreynolds) p%stress_ij(:,ij)=p%stress_ij(:,ij)-trace_factor*p%u2*prefactor*p%rho
              if (lmagnetic) p%stress_ij(:,ij)=p%stress_ij(:,ij)+trace_factor*p%b2
              if (lelectmag) p%stress_ij(:,ij)=p%stress_ij(:,ij)+trace_factor*p%e2
            endif
          endif
        enddo
        enddo
!
!  Possibility of gradually ramping up the stress on time scale tstress_ramp.
!
        if (lstress_ramp) then
          fact=min(real(t-tstart)/tstress_ramp, 1.)
          p%stress_ij(:,:)=p%stress_ij(:,:)*fact
        endif
      endif
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  Several precalculated Pencils of information are passed for
!  efficiency.
!
!  This routine computes various diagnostic quantities, but those
!  would be more easily done further below when sign_switch is known.
!  But this only affects the helicities. The values of EEGW and hrms are
!  however correct and agree with those of gravitational_waves_hij6.
!
!  06-oct-03/tony: coded
!  07-feb-18/axel: added nscale_factor=0 (no expansion), =.5 (radiation era)
!
      use Diagnostics
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: scale_factor, stress_prefactor2, sign_switch=0, fac_stress_comp
      type (pencil_case) :: p
!
      integer :: ij
!
      intent(in) :: p
      intent(inout) :: f,df
!
!  Identify module and boundary conditions.
!
      if (lfirst) then
        if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
!  Compute scale factor.
!  Note: to prevent division by zero, it is best to put tstart=1. in start.in.
!  If that is not done, one can put here tshift=1., for example.
!  If that is not the case either, we put scale_factor=1.
!  At the next timestep, this will no longer be a problem.
!
      if (lreheating_GW) then
        scale_factor=.25*(t+1.)**2
      else
        if (t+tshift==0.) then
          scale_factor=1.
        else
          scale_factor=(t+tshift)**nscale_factor_conformal
        endif
      endif
      stress_prefactor2=stress_prefactor/scale_factor
!
!  Possibilty to compensate against the decaying stress in decaying turbulence.
!
      if (tau_stress_comp>0.) then
        fac_stress_comp=(1.+(t-tstart)/tau_stress_comp)**exp_stress_comp
        stress_prefactor2=stress_prefactor2*fac_stress_comp
      endif
!
!  Possibilty to introduce later kicks by a factor fac_stress_kick.
!
      if (tau_stress_kick>0.) then
        if (t>=tnext_stress_kick) then
          tnext_stress_kick=tnext_stress_kick+tau_stress_kick
          accum_stress_kick=accum_stress_kick+fac_stress_kick
        endif
        stress_prefactor2=stress_prefactor2*accum_stress_kick
      endif
!
!  Assemble rhs of GW equations.
!
      do ij=1,6
        f(l1:l2,m,n,iStress_ij+ij-1)=stress_prefactor2*p%stress_ij(:,ij)
      enddo
!
!  diagnostics
!
       if (ldiagnos) then
         if (lggTX_as_aux) then
           if (idiag_EEGW/=0) call sum_mn_name((f(l1:l2,m,n,iggT)**2+f(l1:l2,m,n,iggTim)**2 &
                                               +f(l1:l2,m,n,iggX)**2+f(l1:l2,m,n,iggXim)**2 &
                                               )*nwgrid*EGWpref,idiag_EEGW)
           if (idiag_gg2m/=0) call sum_mn_name((f(l1:l2,m,n,iggT)**2+f(l1:l2,m,n,iggTim)**2 &
                                               +f(l1:l2,m,n,iggX)**2+f(l1:l2,m,n,iggXim)**2 &
                                               )*nwgrid,idiag_gg2m)
           if (idiag_Stgm/=0) call sum_mn_name((f(l1:l2,m,n,iStressT  )*f(l1:l2,m,n,iggT  ) &
                                               +f(l1:l2,m,n,iStressTim)*f(l1:l2,m,n,iggTim) &
                                               +f(l1:l2,m,n,iStressX  )*f(l1:l2,m,n,iggX  ) &
                                               +f(l1:l2,m,n,iStressXim)*f(l1:l2,m,n,iggXim) &
                                               )*nwgrid,idiag_Stgm)
           if (idiag_ggT2m/=0) call sum_mn_name((f(l1:l2,m,n,iggT)**2+f(l1:l2,m,n,iggTim)**2 &
                                                )*nwgrid,idiag_ggT2m)
           if (idiag_ggX2m/=0) call sum_mn_name((f(l1:l2,m,n,iggX)**2+f(l1:l2,m,n,iggXim)**2 &
                                                )*nwgrid,idiag_ggX2m)
           if (idiag_ggTXm/=0) call sum_mn_name((f(l1:l2,m,n,iggT  )*f(l1:l2,m,n,iggXim) &
                                                -f(l1:l2,m,n,iggTim)*f(l1:l2,m,n,iggX  ) &
                                                )*nwgrid*sign_switch,idiag_ggTXm)
         endif
         if (lhhTX_as_aux) then
           if (idiag_hrms/=0) call sum_mn_name((f(l1:l2,m,n,ihhT)**2+f(l1:l2,m,n,ihhTim)**2 &
                                               +f(l1:l2,m,n,ihhX)**2+f(l1:l2,m,n,ihhXim)**2 &
                                               )*nwgrid,idiag_hrms,lsqrt=.true.)
           if (idiag_hhT2m/=0) call sum_mn_name((f(l1:l2,m,n,ihhT)**2+f(l1:l2,m,n,ihhTim)**2 &
                                                )*nwgrid,idiag_hhT2m)
           if (idiag_hhX2m/=0) call sum_mn_name((f(l1:l2,m,n,ihhX)**2+f(l1:l2,m,n,ihhXim)**2 &
                                                )*nwgrid,idiag_hhX2m)
           if (idiag_hhTXm/=0) call sum_mn_name((f(l1:l2,m,n,ihhT  )*f(l1:l2,m,n,ihhXim) &
                                                -f(l1:l2,m,n,ihhTim)*f(l1:l2,m,n,ihhX  ) &
                                                )*nwgrid*sign_switch,idiag_hhTXm)
         endif
!
         if (lproc_pt.and.m==mpoint.and.n==npoint) then
           if (idiag_STrept/=0) call save_name(f(lpoint,m,n,iStressT  ),idiag_STrept)
           if (idiag_STimpt/=0) call save_name(f(lpoint,m,n,iStressTim),idiag_STimpt)
           if (idiag_SXrept/=0) call save_name(f(lpoint,m,n,iStressX  ),idiag_SXrept)
           if (idiag_SXimpt/=0) call save_name(f(lpoint,m,n,iStressXim),idiag_SXimpt)
           if (idiag_hTrept/=0) call save_name(f(lpoint,m,n,ihhT  ),idiag_hTrept)
           if (idiag_hTimpt/=0) call save_name(f(lpoint,m,n,ihhTim),idiag_hTimpt)
           if (idiag_hXrept/=0) call save_name(f(lpoint,m,n,ihhX  ),idiag_hXrept)
           if (idiag_hXimpt/=0) call save_name(f(lpoint,m,n,ihhXim),idiag_hXimpt)
           if (idiag_gTrept/=0) call save_name(f(lpoint,m,n,iggT  ),idiag_gTrept)
           if (idiag_gTimpt/=0) call save_name(f(lpoint,m,n,iggTim),idiag_gTimpt)
           if (idiag_gXrept/=0) call save_name(f(lpoint,m,n,iggX  ),idiag_gXrept)
           if (idiag_gXimpt/=0) call save_name(f(lpoint,m,n,iggXim),idiag_gXimpt)
           if (idiag_g11pt/=0) call save_name(f(lpoint,m,n,igij+1-1),idiag_g11pt)
           if (idiag_g22pt/=0) call save_name(f(lpoint,m,n,igij+2-1),idiag_g22pt)
           if (idiag_g33pt/=0) call save_name(f(lpoint,m,n,igij+3-1),idiag_g33pt)
           if (idiag_g12pt/=0) call save_name(f(lpoint,m,n,igij+4-1),idiag_g12pt)
           if (idiag_g23pt/=0) call save_name(f(lpoint,m,n,igij+5-1),idiag_g23pt)
           if (idiag_g31pt/=0) call save_name(f(lpoint,m,n,igij+6-1),idiag_g31pt)
           if (lhhTX_as_aux) then
             if (idiag_hhTpt/=0) call save_name(f(lpoint,m,n,ihhT),idiag_hhTpt)
             if (idiag_hhXpt/=0) call save_name(f(lpoint,m,n,ihhX),idiag_hhXpt)
           endif
           if (lreal_space_hTX_as_aux) then
             if (idiag_hhTpt/=0) call save_name(f(lpoint,m,n,ihhT_realspace),idiag_hhTpt)
             if (idiag_hhXpt/=0) call save_name(f(lpoint,m,n,ihhX_realspace),idiag_hhXpt)
             if (idiag_hhTp2/=0) call save_name(f(lpoint2,m,n,ihhT_realspace),idiag_hhTp2)
             if (idiag_hhXp2/=0) call save_name(f(lpoint2,m,n,ihhX_realspace),idiag_hhXp2)
           endif
           if (lggTX_as_aux) then
             if (idiag_ggTpt/=0) call save_name(f(lpoint,m,n,iggT),idiag_ggTpt)
             if (idiag_ggXpt/=0) call save_name(f(lpoint,m,n,iggX),idiag_ggXpt)
           endif
           if (lreal_space_gTX_as_aux) then
             if (idiag_ggTpt/=0) call save_name(f(lpoint,m,n,iggT_realspace),idiag_ggTpt)
             if (idiag_ggXpt/=0) call save_name(f(lpoint,m,n,iggX_realspace),idiag_ggXpt)
             if (idiag_ggTp2/=0) call save_name(f(lpoint2,m,n,iggT_realspace),idiag_ggTp2)
             if (idiag_ggXp2/=0) call save_name(f(lpoint2,m,n,iggX_realspace),idiag_ggXp2)
           endif
         endif
!
         if (lproc_p2.and.m==mpoint2.and.n==npoint2) then
           if (idiag_STrep2/=0) call save_name(f(lpoint2,m,n,iStressT  ),idiag_STrep2)
           if (idiag_STimp2/=0) call save_name(f(lpoint2,m,n,iStressTim),idiag_STimp2)
           if (idiag_SXrep2/=0) call save_name(f(lpoint2,m,n,iStressX  ),idiag_SXrep2)
           if (idiag_SXimp2/=0) call save_name(f(lpoint2,m,n,iStressXim),idiag_SXimp2)
           if (idiag_hTrep2/=0) call save_name(f(lpoint2,m,n,ihhT  ),idiag_hTrep2)
           if (idiag_hTimp2/=0) call save_name(f(lpoint2,m,n,ihhTim),idiag_hTimp2)
           if (idiag_hXrep2/=0) call save_name(f(lpoint2,m,n,ihhX  ),idiag_hXrep2)
           if (idiag_hXimp2/=0) call save_name(f(lpoint2,m,n,ihhXim),idiag_hXimp2)
           if (idiag_gTrep2/=0) call save_name(f(lpoint2,m,n,iggT  ),idiag_gTrep2)
           if (idiag_gTimp2/=0) call save_name(f(lpoint2,m,n,iggTim),idiag_gTimp2)
           if (idiag_gXrep2/=0) call save_name(f(lpoint2,m,n,iggX  ),idiag_gXrep2)
           if (idiag_gXimp2/=0) call save_name(f(lpoint2,m,n,iggXim),idiag_gXimp2)
           if (lhhTX_as_aux) then
             if (idiag_hhTp2/=0) call save_name(f(lpoint2,m,n,ihhT),idiag_hhTp2)
             if (idiag_hhXp2/=0) call save_name(f(lpoint2,m,n,ihhX),idiag_hhXp2)
           endif
           if (lggTX_as_aux) then
             if (idiag_ggTp2/=0) call save_name(f(lpoint2,m,n,iggT),idiag_ggTp2)
             if (idiag_ggXp2/=0) call save_name(f(lpoint2,m,n,iggX),idiag_ggXp2)
           endif
         endif
       endif
      else
        if (headtt.or.ldebug) print*,'dspecial_dt: DONT SOLVE dspecial_dt'
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
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
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  13-may-18/axel: added remove_mean_value for hij and gij
!
      use Sub, only: remove_mean_value
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  07-aug-17/axel: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_,llast)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_
      logical, intent(in) :: llast
!
!  Compute the transverse part of the stress tensor by going into Fourier space.
!
      if (lfirst) call compute_gT_and_gX_from_gij(f,'St')
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dt_)
!
    endsubroutine special_after_timestep
!***********************************************************************
    subroutine make_spectra(f)
!
!  16-oct-19/MR: carved out from special_calc_spectra
!
      use Fourier, only: kx_fft, ky_fft, kz_fft
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: ikx, iky, ikz, q, p, pq, ik
      real :: k1, k2, k3, ksqr,one_over_k2,one_over_k4,sign_switch
      real :: k1mNy, k2mNy, k3mNy, SCL_re, SCL_im
      real, dimension(3) :: VCT_re, VCT_im, kvec
!
      spectra%GWs=0.; spectra%GWshel=0.
      spectra%GWh=0.; spectra%GWhhel=0.
      spectra%GWm=0.; spectra%GWmhel=0.
      spectra%Str=0.; spectra%Strhel=0.
      spectra%Stg=0.; spectra%Stghel=0.
      spectra%SCL=0.; spectra%VCT=0.; spectra%Tpq=0.
      spectra%TGW=0.
!
!  Define negative Nyquist wavenumbers if lswitch_symmetric
!
      if (lswitch_symmetric) then
        k1mNy=kx_fft(nxgrid/2+1)
        k2mNy=ky_fft(nygrid/2+1)
        k3mNy=kz_fft(nzgrid/2+1)
      endif
!
!  Loop over all positions in k-space.
!
      do ikz=1,nz
        do iky=1,ny
          do ikx=1,nx
!
            k1=kx_fft(ikx+ipx*nx)
            k2=ky_fft(iky+ipy*ny)
            k3=kz_fft(ikz+ipz*nz)
!
            ksqr=k1**2+k2**2+k3**2
!
            if (lroot.and.ikx==1.and.iky==1.and.ikz==1) then
              one_over_k2=0.
            else
              one_over_k2=1./ksqr
            endif
            one_over_k4=one_over_k2**2
!
!  possibility of swapping the sign
!
            sign_switch=1.
            if (lswitch_sign_e_X) then
              if (k3<0.) then
                sign_switch=-1.
              elseif (k3==0.) then
                if (k2<0.) then
                  sign_switch=-1.
                elseif (k2==0.) then
                  if (k1<0.) sign_switch=-1.
                endif
              endif
            endif
!
!  Put sign_switch to zero for the negative Nyquist values, because they
!  don't exist for the corresponding positive values and cause asymmetry.
!
            if (lswitch_symmetric) then
              if (k1==k1mNy .or.k2==k2mNy .or.  k3==k3mNy) sign_switch=0.
            endif
!
!  Sum up energy and helicity spectra. Divide by kscale_factor to have integers
!  for the Fortran index ik. Note, however, that 
!
!  set k vector
!
            kvec(1)=k1
            kvec(2)=k2
            kvec(3)=k3
!
            if (SCL_spec) then
              SCL_re=0.
              SCL_im=0.
              do q=1,3
              do p=1,3
                pq=ij_table(p,q)
                SCL_re=SCL_re-1.5*kvec(p)*kvec(q)*one_over_k4*Tpq_re(ikx,iky,ikz,pq)
                SCL_im=SCL_im-1.5*kvec(p)*kvec(q)*one_over_k4*Tpq_im(ikx,iky,ikz,pq)
              enddo
              enddo
            endif
!
!  V_i = -i*ki (2/3) S - (4/3) i*kj/k^4 Tij
!      = -i*ki (2/3) (S'+iS") -  2*i*kj/k^2 (Tij'+iTik")
!
            if (VCT_spec) then
              do q=1,3
                VCT_re(q)=+twothird*kvec(q)*SCL_im
                VCT_im(q)=-twothird*kvec(q)*SCL_re
                do p=1,3
                  pq=ij_table(p,q)
                  VCT_re(q)=VCT_re(q)+2.*kvec(p)*one_over_k2*Tpq_im(ikx,iky,ikz,pq)
                  VCT_im(q)=VCT_im(q)-2.*kvec(p)*one_over_k2*Tpq_re(ikx,iky,ikz,pq)
                enddo
              enddo
            endif
!
            ik=1+nint(sqrt(ksqr)/kscale_factor)
!
!  Debug output
!
            if (ldebug_print) then
              if (ik <= 5) write(*,1000) iproc,ik,k1,k2,k3,f(nghost+ikx,nghost+iky,nghost+ikz,iggX  )
              1000 format(2i5,1p,4e11.2)
            endif
!
            if (ik <= nk) then
!
!  Gravitational wave energy spectrum computed from hdot (=g)
!
              if (GWs_spec) then
                spectra%GWs(ik)=spectra%GWs(ik) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iggX  )**2 &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iggXim)**2 &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iggT  )**2 &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iggTim)**2
                spectra%GWshel(ik)=spectra%GWshel(ik)+2*sign_switch*( &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iggXim) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,iggT  ) &
                   -f(nghost+ikx,nghost+iky,nghost+ikz,iggX  ) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,iggTim) )
              endif
!
              if (GWs_spec_complex) then
                if (k2==0. .and. k3==0.) then
                  spectra%complex_Str_T(ikx)=cmplx(f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  ), &
                                                   f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim))
                  spectra%complex_Str_X(ikx)=cmplx(f(nghost+ikx,nghost+iky,nghost+ikz,ihhX  ), &
                                                   f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim))
                else
                  spectra%complex_Str_T(ikx)=0.
                  spectra%complex_Str_X(ikx)=0.
                endif
              endif
!
!  Gravitational wave strain spectrum computed from h
!
              if (GWh_spec) then
                spectra%GWh(ik)=spectra%GWh(ik) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,ihhX  )**2 &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim)**2 &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  )**2 &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim)**2
                spectra%GWhhel(ik)=spectra%GWhhel(ik)+2*sign_switch*( &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  ) &
                   -f(nghost+ikx,nghost+iky,nghost+ikz,ihhX  ) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim) )
              endif
!
!  Gravitational wave mixed spectrum computed from h and g
!
              if (GWm_spec) then
                spectra%GWm(ik)=spectra%GWm(ik) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,ihhX)  *f(nghost+ikx,nghost+iky,nghost+ikz,iggX  ) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim)*f(nghost+ikx,nghost+iky,nghost+ikz,iggXim) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,ihhT)  *f(nghost+ikx,nghost+iky,nghost+ikz,iggT  ) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim)*f(nghost+ikx,nghost+iky,nghost+ikz,iggTim)
                spectra%GWmhel(ik)=spectra%GWmhel(ik)-sign_switch*( &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,iggT  ) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iggXim) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  ) &
                   -f(nghost+ikx,nghost+iky,nghost+ikz,ihhX  ) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,iggTim) &
                   -f(nghost+ikx,nghost+iky,nghost+ikz,iggX  ) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim) )
              endif
!
!  Gravitational wave production spectrum for TTgT
!
              if (Stg_spec) then
                spectra%Stg(ik)=spectra%Stg(ik) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iStressT)  *f(nghost+ikx,nghost+iky,nghost+ikz,iggT  ) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iStressTim)*f(nghost+ikx,nghost+iky,nghost+ikz,iggTim) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iStressX)  *f(nghost+ikx,nghost+iky,nghost+ikz,iggX  ) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iStressXim)*f(nghost+ikx,nghost+iky,nghost+ikz,iggXim)
                spectra%Stghel(ik)=spectra%Stghel(ik)-sign_switch*( &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iStressXim) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,iggT  ) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iStressXim) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  ) &
                   -f(nghost+ikx,nghost+iky,nghost+ikz,iStressX  ) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,iggTim) &
                   -f(nghost+ikx,nghost+iky,nghost+ikz,iStressX  ) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim) )
              endif
!
!  Stress spectrum computed from Str
!  ?not used currently
!
              if (Str_spec) then
                spectra%Str(ik)=spectra%Str(ik) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iStressX  )**2 &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iStressXim)**2 &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iStressT  )**2 &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iStressTim)**2
                spectra%Strhel(ik)=spectra%Strhel(ik)+2*sign_switch*( &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iStressXim) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,iStressT  ) &
                   -f(nghost+ikx,nghost+iky,nghost+ikz,iStressX  ) &
                   *f(nghost+ikx,nghost+iky,nghost+ikz,iStressTim) )
              endif
!
!  Stress spectrum computed from the scalar mode, SCL.
!
              if (SCL_spec) then
                spectra%SCL(ik)=spectra%SCL(ik)+.5*(SCL_re**2+SCL_im**2)
              endif
!
!  Spectrum computed from the vector modes...
!
              if (VCT_spec) then
                do q=1,3
                  spectra%VCT(ik)=spectra%VCT(ik)+.5*(VCT_re(q)**2+VCT_im(q)**2)
                enddo
              endif
!
!  Spectrum computed from the total unprojected stress
!
              if (Tpq_spec) then
                do q=1,3
                do p=1,3
                  pq=ij_table(p,q)
                  spectra%Tpq(ik)=spectra%Tpq(ik)+.5*(Tpq_re(ikx,iky,ikz,pq)**2+Tpq_im(ikx,iky,ikz,pq)**2)
                enddo
                enddo
              endif
!
! Added for nonlinear GW memory effect
!
            if (TGW_spec) then
              if (.not.lnonlinear_Tpq_trans.and.lroot) print*,'WARNING: TGW_spec incorrect; nonlinear_Tpq is still in real space'
              do q=1,3
              do p=1,3
                pq=ij_table(p,q)
                spectra%TGW(ik)=spectra%TGW(ik)+.5*(nonlinear_Tpq_re(ikx,iky,ikz,pq)**2+nonlinear_Tpq_im(ikx,iky,ikz,pq)**2)
              enddo
              enddo
            endif
            
          endif   
        enddo
      enddo
    enddo

    endsubroutine make_spectra
!***********************************************************************
    subroutine special_calc_spectra(f,spectrum,spectrum_hel,lfirstcall,kind)
!
!  Calculates GW spectra. For use with a single special module.
!
!  16-oct-19/MR: carved out from compute_gT_and_gX_from_gij
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:) :: spectrum,spectrum_hel
      logical :: lfirstcall
      character(LEN=3) :: kind

      if (lfirstcall) then
        call make_spectra(f)
        lfirstcall=.false.
      endif

      select case(kind)
      case ('GWs'); spectrum=spectra%GWs; spectrum_hel=spectra%GWshel
      case ('GWh'); spectrum=spectra%GWh; spectrum_hel=spectra%GWhhel
      case ('GWm'); spectrum=spectra%GWm; spectrum_hel=spectra%GWmhel
      case ('Str'); spectrum=spectra%Str; spectrum_hel=spectra%Strhel 
      case ('Stg'); spectrum=spectra%Stg; spectrum_hel=spectra%Stghel 
      case ('SCL'); spectrum=spectra%SCL; spectrum_hel=0. 
      case ('VCT'); spectrum=spectra%VCT; spectrum_hel=0. 
      case ('Tpq'); spectrum=spectra%Tpq; spectrum_hel=0. 
      case ('TGW'); spectrum=spectra%TGW; spectrum_hel=0.
      case ('StT'); spectrum=real(spectra%complex_Str_T)
                    spectrum_hel=aimag(spectra%complex_Str_T)
      case ('StX'); spectrum=real(spectra%complex_Str_X)
                    spectrum_hel=aimag(spectra%complex_Str_X)
      case default; if (lroot) call warning('special_calc_spectra', &
                      'kind of spectrum "'//kind//'" not implemented')
      endselect

    endsubroutine special_calc_spectra
!***********************************************************************
    subroutine special_calc_spectra_byte(f,spectrum,spectrum_hel,lfirstcall,kind,len)
!
!  Calculates GW spectra. For use with multiple special modules.
!
!  16-oct-19/MR: carved out from compute_gT_and_gX_from_gij
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nk) :: spectrum,spectrum_hel
      logical :: lfirstcall
      integer(KIND=ikind1), dimension(3) :: kind
      integer :: len

      character(LEN=3) :: kindstr

      if (lfirstcall) then
        call make_spectra(f)
        lfirstcall=.false.
      endif

      kindstr=char(kind(1))//char(kind(2))//char(kind(3))

      select case(kindstr)
      case ('GWs'); spectrum=spectra%GWs; spectrum_hel=spectra%GWshel
      case ('GWh'); spectrum=spectra%GWh; spectrum_hel=spectra%GWhhel
      case ('GWm'); spectrum=spectra%GWm; spectrum_hel=spectra%GWmhel
      case ('Str'); spectrum=spectra%Str; spectrum_hel=spectra%Strhel 
      case ('Stg'); spectrum=spectra%Stg; spectrum_hel=spectra%Stghel 
      case ('SCL'); spectrum=spectra%SCL; spectrum_hel=0. 
      case ('VCT'); spectrum=spectra%VCT; spectrum_hel=0. 
      case ('Tpq'); spectrum=spectra%Tpq; spectrum_hel=0. 
      case ('TGW'); spectrum=spectra%TGW; spectrum_hel=0.
      case default; if (lroot) call warning('special_calc_spectra', &
                      'kind of spectrum "'//kindstr//'" not implemented')
      endselect

    endsubroutine special_calc_spectra_byte
!***********************************************************************
    subroutine compute_gT_and_gX_from_gij(f,label)
!
!  Compute the transverse part of the stress tensor by going into Fourier space.
!
!  07-aug-17/axel: coded
!
      use Fourier, only: fourier_transform, fft_xyz_parallel, kx_fft, ky_fft, kz_fft
      use SharedVariables, only: put_shared_variable
!
      real, dimension (:,:,:), allocatable :: S_T_re, S_T_im, S_X_re, S_X_im, g2T_re, g2T_im, g2X_re, g2X_im
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (6) :: Pij=0., e_T, e_X, Sij_re, Sij_im, delij=0.
      real, dimension (:,:,:,:,:), allocatable :: Hijkre, Hijkim
      real, dimension (3) :: e1, e2, kvec
      integer :: i,j,p,q,ik,ikx,iky,ikz,stat,ij,pq,ip,jq,jStress_ij
      real :: fact, delkt
      real :: ksqr, one_over_k2, k1, k2, k3, k1sqr, k2sqr, k3sqr
      real :: hhTre, hhTim, hhXre, hhXim, coefAre, coefAim
      real :: ggTre, ggTim, ggXre, ggXim, coefBre, coefBim
      real :: cosot, sinot, sinot_minus, om12, om, om1, om2
      intent(inout) :: f
      character (len=2) :: label
      logical :: lsign_om2
!
!  Check that the relevant arrays are registered
!
      if (.not.lStress_as_aux.and.label=='St') call fatal_error('compute_gT_and_gX_from_gij','lStress_as_aux must be true')
!
!  For testing purposes, if lno_transverse_part=T, we would not need to
!  compute the Fourier transform, so we would skip the rest.
!
!  Allocate memory for arrays.
!
      allocate(S_T_re(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_T_re')
!
      allocate(S_T_im(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_T_im')
!
      allocate(S_X_re(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_X_re')
!
      allocate(S_X_im(nx,ny,nz),stat=stat)
      if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for S_X_im')
!
!  Allocate 18 chunks of memory for nonlinear source
!
      if (lnonlinear_source) then
        allocate(Hijkre(nx,ny,nz,3,6),stat=stat)
        if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for Hijkre')
!
        allocate(Hijkim(nx,ny,nz,3,6),stat=stat)
        if (stat>0) call fatal_error('compute_gT_and_gX_from_gij','Could not allocate memory for Hijkim')
      endif
!
!  set delta_ij
!
      delij(1:3)=1.
      delij(4:6)=0.
!
!  Set k^2 array. Note that in Fourier space, kz is the fastest index and has
!  the full nx extent (which, currently, must be equal to nxgrid).
!  But call it one_over_k2.
!  Added for computation of Hijk
!  Begin by setting nonlinear_Tpq_re=0.
!
      if (lnonlinear_source) then
        nonlinear_Tpq_re=0.
        do ikz=1,nz
          do iky=1,ny
            do ikx=1,nx
!
!  compute e_T and e_X; determine first preferred direction,
!  which is a component with the smallest component by modulus.
!  Note: this is duplicated code and any changes here need to be done also elsewhere.
!
              k1=kx_fft(ikx+ipx*nx)
              k2=ky_fft(iky+ipy*ny)
              k3=kz_fft(ikz+ipz*nz)
              kvec(1)=k1
              kvec(2)=k2
              kvec(3)=k3
              k1sqr=k1**2
              k2sqr=k2**2
              k3sqr=k3**2
              ksqr=k1sqr+k2sqr+k3sqr
!
              if (lroot.and.ikx==1.and.iky==1.and.ikz==1) then
                e1=0.
                e2=0.
              else
!
                if(abs(k1)<abs(k2)) then
                  if(abs(k1)<abs(k3)) then !(k1 is pref dir)
                    e1=(/0.,-k3,+k2/)
                    e2=(/k2sqr+k3sqr,-k2*k1,-k3*k1/)
                  else !(k3 is pref dir)
                    e1=(/k2,-k1,0./)
                    e2=(/k1*k3,k2*k3,-(k1sqr+k2sqr)/)
                  endif
                else !(k2 smaller than k1)
                  if(abs(k2)<abs(k3)) then !(k2 is pref dir)
                    e1=(/-k3,0.,+k1/)
                    e2=(/+k1*k2,-(k1sqr+k3sqr),+k3*k2/)
                  else !(k3 is pref dir)
                    e1=(/k2,-k1,0./)
                    e2=(/k1*k3,k2*k3,-(k1sqr+k2sqr)/)
                  endif
                endif
                e1=e1/sqrt(e1(1)**2+e1(2)**2+e1(3)**2)
                e2=e2/sqrt(e2(1)**2+e2(2)**2+e2(3)**2)
              endif
!
!  compute e_T and e_X
!
              do j=1,3
              do i=1,3
                ij=ij_table(i,j)
                e_T(ij)=e1(i)*e1(j)-e2(i)*e2(j)
                e_X(ij)=e1(i)*e2(j)+e2(i)*e1(j)
              enddo
              enddo
!
!  possibility of swapping the sign of e_X
!
              if (lswitch_sign_e_X) then
                if (k3<0.) then
                  e_X=-e_X
                elseif (k3==0.) then
                  if (k2<0.) then
                    e_X=-e_X
                  elseif (k2==0.) then
                    if (k1<0.) then
                      e_X=-e_X
                    endif
                  endif
                endif
              endif
!
!  Compute exact solution for hT, hX, gT, and gX in Fourier space.
!
              hhTre=f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  )
              hhXre=f(nghost+ikx,nghost+iky,nghost+ikz,ihhX  )
              hhTim=f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim)
              hhXim=f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim)
!
! compute Hijk for nonlinear GW memory effect (still in Fourier space)
!             
              do i=1,3 
              do j=1,6 
                Hijkim(ikx,iky,ikz,i,j)=+kvec(i)*(e_T(j)*hhTre+e_X(j)*hhXre) 
                Hijkre(ikx,iky,ikz,i,j)=-kvec(i)*(e_T(j)*hhTim+e_X(j)*hhXim)
              enddo
              enddo
!
!  end of ikx, iky, and ikz loops
!
            enddo
          enddo
        enddo
!
!  Go back with Hijkre into real space:
!
        do i=1,3
        do j=1,6
          call fft_xyz_parallel(Hijkre(:,:,:,i,j),Hijkim(:,:,:,i,j))
        enddo
        enddo
!
!  Now we compute nonlinear_Tpq_re in real space, so the imaginary
!  part must be zero and is not used.
!
        do i=1,3
        do j=1,3
          ij=ij_table(i,j)
          do p=1,3
          do q=1,3
            pq=ij_table(p,q)
            nonlinear_Tpq_re(:,:,:,pq)=nonlinear_Tpq_re(:,:,:,pq) &
                +Hijkre(:,:,:,p,ij)*Hijkre(:,:,:,q,ij)
! -             +Hijkre(:,:,:,p,ij)*Hijkre(:,:,:,q,ij) &
! -             -Hijkim(:,:,:,p,ij)*Hijkim(:,:,:,q,ij)
          enddo
          enddo
        enddo
        enddo
!
! end of if condition for nonlinear_source
!
      endif
!
!  Assemble stress, Tpq, and transform to Fourier space.
!  Add nonlinear source here, before transforming
!
      if (label=='St') then
        Tpq_re(:,:,:,:)=f(l1:l2,m1:m2,n1:n2,iStress_ij:iStress_ij+5)
        Tpq_im=0.0
        if (lnonlinear_source) then
          if (nonlinear_source_fact/=0.) nonlinear_Tpq_re=nonlinear_Tpq_re*nonlinear_source_fact
          if (lnonlinear_Tpq_trans) then
            nonlinear_Tpq_im=0.
            call fft_xyz_parallel(nonlinear_Tpq_re(:,:,:,:),nonlinear_Tpq_im(:,:,:,:))
          else
            Tpq_re(:,:,:,:)=Tpq_re(:,:,:,:)+nonlinear_Tpq_re(:,:,:,:)
          endif
        endif
        call fft_xyz_parallel(Tpq_re(:,:,:,:),Tpq_im(:,:,:,:))
      endif
!
!  determine time-dependent delkt
!
      delkt=delk
      if (ldelkt) then
        if (t>tdelk) delkt=0.
      endif
!
!  Set ST=SX=0 and reset all spectra.
!
      S_T_re=0. ; S_T_im=0.
      S_X_re=0. ; S_X_im=0.
!
!  P11, P22, P33, P12, P23, P31
!
      do ikz=1,nz
        do iky=1,ny
          do ikx=1,nx
!
!  compute e_T and e_X; determine first preferred direction,
!  which is a component with the smallest component by modulus.
!
            k1=kx_fft(ikx+ipx*nx)
            k2=ky_fft(iky+ipy*ny)
            k3=kz_fft(ikz+ipz*nz)
            k1sqr=k1**2
            k2sqr=k2**2
            k3sqr=k3**2
            ksqr=k1sqr+k2sqr+k3sqr
!
!  find two vectors e1 and e2 to compute e_T and e_X
!
            if (lroot.and.ikx==1.and.iky==1.and.ikz==1) then
              e1=0.
              e2=0.
              Pij(1)=0.
              Pij(2)=0.
              Pij(3)=0.
              Pij(4)=0.
              Pij(5)=0.
              Pij(6)=0.
              om=0.
            else
!
!  compute omega (but assume c=1), omega*t, etc.
!
              one_over_k2=1./ksqr
              if (linflation) then
                om2=4.*ksqr-2./t**2
                lsign_om2=(om2 >= 0.)
                om=sqrt(abs(om2))
              elseif (lreheating_GW) then
                om2=ksqr-2./(t+1.)**2
                lsign_om2=(om2 >= 0.)
                om=sqrt(abs(om2))
              else
                if (delkt/=0.) then
                  om=sqrt(ksqr+delkt**2)
                else
                  om=sqrt(ksqr)
                endif
                lsign_om2=.true.
              endif
!
              if(abs(k1)<abs(k2)) then
                if(abs(k1)<abs(k3)) then !(k1 is pref dir)
                  e1=(/0.,-k3,+k2/)
                  e2=(/k2sqr+k3sqr,-k2*k1,-k3*k1/)
                else !(k3 is pref dir)
                  e1=(/k2,-k1,0./)
                  e2=(/k1*k3,k2*k3,-(k1sqr+k2sqr)/)
                endif
              else !(k2 smaller than k1)
                if(abs(k2)<abs(k3)) then !(k2 is pref dir)
                  e1=(/-k3,0.,+k1/)
                  e2=(/+k1*k2,-(k1sqr+k3sqr),+k3*k2/)
                else !(k3 is pref dir)
                  e1=(/k2,-k1,0./)
                  e2=(/k1*k3,k2*k3,-(k1sqr+k2sqr)/)
                endif
              endif
              e1=e1/sqrt(e1(1)**2+e1(2)**2+e1(3)**2)
              e2=e2/sqrt(e2(1)**2+e2(2)**2+e2(3)**2)
              Pij(1)=1.-k1sqr*one_over_k2
              Pij(2)=1.-k2sqr*one_over_k2
              Pij(3)=1.-k3sqr*one_over_k2
              Pij(4)=-k1*k2*one_over_k2
              Pij(5)=-k2*k3*one_over_k2
              Pij(6)=-k3*k1*one_over_k2
            endif
!
!  compute e_T and e_X
!
            do j=1,3
            do i=1,3
              ij=ij_table(i,j)
              e_T(ij)=e1(i)*e1(j)-e2(i)*e2(j)
              e_X(ij)=e1(i)*e2(j)+e2(i)*e1(j)
            enddo
            enddo
!
!  possibility of swapping the sign of e_X
!
            if (lswitch_sign_e_X) then
              if (k3<0.) then
                e_X=-e_X
              elseif (k3==0.) then
                if (k2<0.) then
                  e_X=-e_X
                elseif (k2==0.) then
                  if (k1<0.) then
                    e_X=-e_X
                  endif
                endif
              endif
            endif
!
!  Traceless-tansverse projection:
!  Sij = (Pip*Pjq-.5*Pij*Ppq)*Tpq
!  MR: Tpq should not be used if (label/='St')
!
            Sij_re=0.
            Sij_im=0.
            do j=1,3
            do i=1,j
            do q=1,3
            do p=1,3
              ij=ij_table(i,j)
              pq=ij_table(p,q)
              ip=ij_table(i,p)
              jq=ij_table(j,q)
              Sij_re(ij)=Sij_re(ij)+(Pij(ip)*Pij(jq)-.5*Pij(ij)*Pij(pq))*Tpq_re(ikx,iky,ikz,pq)
              Sij_im(ij)=Sij_im(ij)+(Pij(ip)*Pij(jq)-.5*Pij(ij)*Pij(pq))*Tpq_im(ikx,iky,ikz,pq)
              if (lnonlinear_source.and.lnonlinear_Tpq_trans) then
                Sij_re(ij)=Sij_re(ij)+(Pij(ip)*Pij(jq)-.5*Pij(ij)*Pij(pq))*nonlinear_Tpq_re(ikx,iky,ikz,pq)
                Sij_im(ij)=Sij_im(ij)+(Pij(ip)*Pij(jq)-.5*Pij(ij)*Pij(pq))*nonlinear_Tpq_im(ikx,iky,ikz,pq)
              endif
            enddo
            enddo
            enddo
            enddo
!
!  Compute S_T and S_X. Loop over all i and j for simplicity to avoid
!  special treatments of diagonal and off-diagonal terms,
!  AB: it looks like one could reuse part of Tpq_re and Tpq_im for S_T_re, ...,
!  AB: S_X_im to save memory
!
            do j=1,3
            do i=1,3
              ij=ij_table(i,j)
              S_T_re(ikx,iky,ikz)=S_T_re(ikx,iky,ikz)+.5*e_T(ij)*Sij_re(ij)
              S_T_im(ikx,iky,ikz)=S_T_im(ikx,iky,ikz)+.5*e_T(ij)*Sij_im(ij)
              S_X_re(ikx,iky,ikz)=S_X_re(ikx,iky,ikz)+.5*e_X(ij)*Sij_re(ij)
              S_X_im(ikx,iky,ikz)=S_X_im(ikx,iky,ikz)+.5*e_X(ij)*Sij_im(ij)
            enddo
            enddo
!
!  Compute exact solution for hT, hX, gT, and gX in Fourier space.
!
            hhTre=f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  )
            hhXre=f(nghost+ikx,nghost+iky,nghost+ikz,ihhX  )
            hhTim=f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim)
            hhXim=f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim)
!
            ggTre=f(nghost+ikx,nghost+iky,nghost+ikz,iggT  )
            ggXre=f(nghost+ikx,nghost+iky,nghost+ikz,iggX  )
            ggTim=f(nghost+ikx,nghost+iky,nghost+ikz,iggTim)
            ggXim=f(nghost+ikx,nghost+iky,nghost+ikz,iggXim)
!
!  compute cos(om*dt) and sin(om*dt) to get from one timestep to the next.
!
            if (om/=0.) then

              om1=1./om
              om12=om1**2
!
!  check whether om^2 is positive. If om^2 is positive, we have the standard
!  rotation matrix, whose third element is negative (sinot_minus), but
!  if om^2 is negative, we have cosh and sinh, always with a plus sign.
!
              if (lsign_om2) then
                cosot=cos(om*dt)
                sinot=sin(om*dt)
                sinot_minus=-sinot
              else
                cosot=cosh(om*dt)
                sinot=sinh(om*dt)
                sinot_minus=+sinot
              endif
!
!  Solve wave equation for hT and gT from one timestep to the next.
!
              coefAre=(hhTre-om12*S_T_re(ikx,iky,ikz))
              coefAim=(hhTim-om12*S_T_im(ikx,iky,ikz))
              coefBre=ggTre*om1
              coefBim=ggTim*om1
              f(nghost+ikx,nghost+iky,nghost+ikz,ihhT  )=coefAre*cosot+coefBre*sinot+om12*S_T_re(ikx,iky,ikz)
              f(nghost+ikx,nghost+iky,nghost+ikz,ihhTim)=coefAim*cosot+coefBim*sinot+om12*S_T_im(ikx,iky,ikz)
              f(nghost+ikx,nghost+iky,nghost+ikz,iggT  )=coefBre*cosot*om+coefAre*om*sinot_minus
              f(nghost+ikx,nghost+iky,nghost+ikz,iggTim)=coefBim*cosot*om+coefAim*om*sinot_minus
!
!  Debug output
!
              if (ldebug_print) then
                if (nint(k1)==2.and.nint(k2)==0.and.nint(k3)==0) then
                  print*,'AXEL0: ',coefAre,coefBre,hhTre,ggTre
                endif
              endif
!
!  Solve wave equation for hX and gX from one timestep to the next.
!
              coefAre=(hhXre-om12*S_X_re(ikx,iky,ikz))
              coefAim=(hhXim-om12*S_X_im(ikx,iky,ikz))
              coefBre=ggXre*om1
              coefBim=ggXim*om1
              f(nghost+ikx,nghost+iky,nghost+ikz,ihhX  )=coefAre*cosot+coefBre*sinot+om12*S_X_re(ikx,iky,ikz)
              f(nghost+ikx,nghost+iky,nghost+ikz,ihhXim)=coefAim*cosot+coefBim*sinot+om12*S_X_im(ikx,iky,ikz)
              f(nghost+ikx,nghost+iky,nghost+ikz,iggX  )=coefBre*cosot*om+coefAre*om*sinot_minus
              f(nghost+ikx,nghost+iky,nghost+ikz,iggXim)=coefBim*cosot*om+coefAim*om*sinot_minus

            else
!
!  Set origin to zero. It is given by (1,1,1) on root processor.
!
              f(nghost+1,nghost+1,nghost+1,ihhT  ) = 0.
              f(nghost+1,nghost+1,nghost+1,ihhTim) = 0.
              f(nghost+1,nghost+1,nghost+1,iggT  ) = 0.
              f(nghost+1,nghost+1,nghost+1,iggTim) = 0.
!
              f(nghost+1,nghost+1,nghost+1,ihhX  ) = 0.
              f(nghost+1,nghost+1,nghost+1,ihhXim) = 0.
              f(nghost+1,nghost+1,nghost+1,iggX  ) = 0.
              f(nghost+1,nghost+1,nghost+1,iggXim) = 0.

            endif
!
!  Set stress components in f-array.
!
            f(nghost+ikx,nghost+iky,nghost+ikz,iStressT  )=S_T_re(ikx,iky,ikz)
            f(nghost+ikx,nghost+iky,nghost+ikz,iStressTim)=S_T_im(ikx,iky,ikz)
            f(nghost+ikx,nghost+iky,nghost+ikz,iStressX  )=S_X_re(ikx,iky,ikz)
            f(nghost+ikx,nghost+iky,nghost+ikz,iStressXim)=S_X_im(ikx,iky,ikz)
!
!  end of ikx, iky, and ikz loops
!
          enddo
        enddo
      enddo
!
!  back to real space: hTX
!
      if (lreal_space_hTX_as_aux) then
        S_T_re=f(l1:l2,m1:m2,n1:n2,ihhT  )
        S_X_re=f(l1:l2,m1:m2,n1:n2,ihhX  )
        S_T_im=f(l1:l2,m1:m2,n1:n2,ihhTim)
        S_X_im=f(l1:l2,m1:m2,n1:n2,ihhXim)
        call fft_xyz_parallel(S_T_re,S_T_im,linv=.true.)
        call fft_xyz_parallel(S_X_re,S_X_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ihhT_realspace)=S_T_re
        f(l1:l2,m1:m2,n1:n2,ihhX_realspace)=S_X_re
      endif
!
!  back to real space: gTX
!
      if (lreal_space_gTX_as_aux) then
        S_T_re=f(l1:l2,m1:m2,n1:n2,iggT  )
        S_X_re=f(l1:l2,m1:m2,n1:n2,iggX  )
        S_T_im=f(l1:l2,m1:m2,n1:n2,iggTim)
        S_X_im=f(l1:l2,m1:m2,n1:n2,iggXim)
        call fft_xyz_parallel(S_T_re,S_T_im,linv=.true.)
        call fft_xyz_parallel(S_X_re,S_X_im,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,iggT_realspace)=S_T_re
        f(l1:l2,m1:m2,n1:n2,iggX_realspace)=S_X_re
      endif
!
    endsubroutine compute_gT_and_gX_from_gij
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics
!!      use FArrayManager, only: farray_index_append
!
      integer :: iname
      logical :: lreset,lwrite
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
        idiag_STrept=0; idiag_STimpt=0; idiag_SXrept=0; idiag_SXimpt=0
        idiag_hTrept=0; idiag_hTimpt=0; idiag_hXrept=0; idiag_hXimpt=0
        idiag_gTrept=0; idiag_gTimpt=0; idiag_gXrept=0; idiag_gXimpt=0
        idiag_STrep2=0; idiag_STimp2=0; idiag_SXrep2=0; idiag_SXimp2=0
        idiag_hTrep2=0; idiag_hTimp2=0; idiag_hXrep2=0; idiag_hXimp2=0
        idiag_gTrep2=0; idiag_gTimp2=0; idiag_gXrep2=0; idiag_gXimp2=0
        idiag_g11pt=0; idiag_g22pt=0; idiag_g33pt=0
        idiag_g12pt=0; idiag_g23pt=0; idiag_g31pt=0
        idiag_hhTpt=0; idiag_hhXpt=0; idiag_ggTpt=0; idiag_ggXpt=0
        idiag_hhTp2=0; idiag_hhXp2=0; idiag_ggTp2=0; idiag_ggXp2=0
        idiag_hhT2m=0; idiag_hhX2m=0; idiag_hhTXm=0; idiag_hrms=0
        idiag_ggT2m=0; idiag_ggX2m=0; idiag_ggTXm=0; idiag_gg2m=0
        idiag_Stgm=0 ; idiag_EEGW=0
        cformv=''
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'STrept',idiag_STrept)
        call parse_name(iname,cname(iname),cform(iname),'STimpt',idiag_STimpt)
        call parse_name(iname,cname(iname),cform(iname),'SXrept',idiag_SXrept)
        call parse_name(iname,cname(iname),cform(iname),'SXimpt',idiag_SXimpt)
        call parse_name(iname,cname(iname),cform(iname),'hTrept',idiag_hTrept)
        call parse_name(iname,cname(iname),cform(iname),'hTimpt',idiag_hTimpt)
        call parse_name(iname,cname(iname),cform(iname),'hXrept',idiag_hXrept)
        call parse_name(iname,cname(iname),cform(iname),'hXimpt',idiag_hXimpt)
        call parse_name(iname,cname(iname),cform(iname),'gTrept',idiag_gTrept)
        call parse_name(iname,cname(iname),cform(iname),'gTimpt',idiag_gTimpt)
        call parse_name(iname,cname(iname),cform(iname),'gXrept',idiag_gXrept)
        call parse_name(iname,cname(iname),cform(iname),'gXimpt',idiag_gXimpt)
        call parse_name(iname,cname(iname),cform(iname),'g11pt',idiag_g11pt)
        call parse_name(iname,cname(iname),cform(iname),'g22pt',idiag_g22pt)
        call parse_name(iname,cname(iname),cform(iname),'g33pt',idiag_g33pt)
        call parse_name(iname,cname(iname),cform(iname),'g12pt',idiag_g12pt)
        call parse_name(iname,cname(iname),cform(iname),'g23pt',idiag_g23pt)
        call parse_name(iname,cname(iname),cform(iname),'g31pt',idiag_g31pt)
        call parse_name(iname,cname(iname),cform(iname),'STrep2',idiag_STrep2)
        call parse_name(iname,cname(iname),cform(iname),'STimp2',idiag_STimp2)
        call parse_name(iname,cname(iname),cform(iname),'SXrep2',idiag_SXrep2)
        call parse_name(iname,cname(iname),cform(iname),'SXimp2',idiag_SXimp2)
        call parse_name(iname,cname(iname),cform(iname),'hTrep2',idiag_hTrep2)
        call parse_name(iname,cname(iname),cform(iname),'hTimp2',idiag_hTimp2)
        call parse_name(iname,cname(iname),cform(iname),'hXrep2',idiag_hXrep2)
        call parse_name(iname,cname(iname),cform(iname),'hXimp2',idiag_hXimp2)
        call parse_name(iname,cname(iname),cform(iname),'gTrep2',idiag_gTrep2)
        call parse_name(iname,cname(iname),cform(iname),'gTimp2',idiag_gTimp2)
        call parse_name(iname,cname(iname),cform(iname),'gXrep2',idiag_gXrep2)
        call parse_name(iname,cname(iname),cform(iname),'gXimp2',idiag_gXimp2)
        if (lhhTX_as_aux) then
          call parse_name(iname,cname(iname),cform(iname),'hhTpt',idiag_hhTpt)
          call parse_name(iname,cname(iname),cform(iname),'hhXpt',idiag_hhXpt)
          call parse_name(iname,cname(iname),cform(iname),'hhTp2',idiag_hhTp2)
          call parse_name(iname,cname(iname),cform(iname),'hhXp2',idiag_hhXp2)
          call parse_name(iname,cname(iname),cform(iname),'hrms',idiag_hrms)
          call parse_name(iname,cname(iname),cform(iname),'hhT2m',idiag_hhT2m)
          call parse_name(iname,cname(iname),cform(iname),'hhX2m',idiag_hhX2m)
          call parse_name(iname,cname(iname),cform(iname),'hhTXm',idiag_hhTXm)
        endif
        if (lggTX_as_aux) then
          call parse_name(iname,cname(iname),cform(iname),'ggTpt',idiag_ggTpt)
          call parse_name(iname,cname(iname),cform(iname),'ggXpt',idiag_ggXpt)
          call parse_name(iname,cname(iname),cform(iname),'ggTp2',idiag_ggTp2)
          call parse_name(iname,cname(iname),cform(iname),'ggXp2',idiag_ggXp2)
          call parse_name(iname,cname(iname),cform(iname),'EEGW',idiag_EEGW)
          call parse_name(iname,cname(iname),cform(iname),'gg2m',idiag_gg2m)
          call parse_name(iname,cname(iname),cform(iname),'Stgm',idiag_Stgm)
          call parse_name(iname,cname(iname),cform(iname),'ggT2m',idiag_ggT2m)
          call parse_name(iname,cname(iname),cform(iname),'ggX2m',idiag_ggX2m)
          call parse_name(iname,cname(iname),cform(iname),'ggTXm',idiag_ggTXm)
        endif
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='hhT'.or.cnamev=='hhX'.or.cnamev=='ggT'.or.cnamev=='ggX'.or. &
              cnamev=='hhTre'.or.cnamev=='hhTim'.or. &
              cnamev=='StTre'.or.cnamev=='StTim' &
             ) cformv='DEFINED'
      endif
!
!!!  write column where which magnetic variable is stored
!!      if (lwrite) then
!!        call farray_index_append('i_SPECIAL_DIAGNOSTIC',i_SPECIAL_DIAGNOSTIC)
!!      endif
!!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
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
!  hhT
!
        case ('hhT')
          if (lreal_space_hTX_as_aux) then
            call assign_slices_scal(slices,f,ihhT_realspace)
          else
            call assign_slices_scal(slices,f,ihhT)
          endif
!
!  hhTre
!
        case ('hhTre')
          call assign_slices_scal(slices,f,ihhT)
!
!  hhTim
!
        case ('hhTim')
          call assign_slices_scal(slices,f,ihhTim)
!
!  hhX
!
        case ('hhX')
          if (lreal_space_hTX_as_aux) then
            call assign_slices_scal(slices,f,ihhX_realspace)
          else
            call assign_slices_scal(slices,f,ihhX)
          endif
!
!  ggT
!
        case ('ggT')
          if (lreal_space_hTX_as_aux) then
            call assign_slices_scal(slices,f,iggT_realspace)
          else
            call assign_slices_scal(slices,f,iggT)
          endif
!
!  ggX
!
        case ('ggX')
          if (lreal_space_hTX_as_aux) then
            call assign_slices_scal(slices,f,iggX_realspace)
          else
            call assign_slices_scal(slices,f,iggX)
          endif
!
!  StTre
!
        case ('StTre')
          call assign_slices_scal(slices,f,iStressT)
!
!  StTim
!
        case ('StTim')
          call assign_slices_scal(slices,f,iStressTim)
!
      endselect
!
!  The following is just a comment to remind ourselves how
!  the remaining 3 offdiagonal terms are being accessed.
!
      !ij_table(1,2)=4
      !ij_table(2,3)=5
      !ij_table(3,1)=6
      !ij_table(2,1)=4
      !ij_table(3,2)=5
      !ij_table(1,3)=6
!
    endsubroutine get_slices_special
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
