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
!  Special variable registration                   | register_magnetic
!    (pre parameter read)                          |
!  Special variable initialization                 | initialize_magnetic
!    (post parameter read)                         |
!  Special variable finalization                   | finalize_magnetic
!    (deallocation, etc.)                          |
!                                                  |
!  Special initial condition                       | init_magnetic
!   this is called last so may be used to modify   |
!   the mvar variables declared by this module     |
!   or optionally modify any of the other f array  |
!   variables.  The latter, however, should be     |
!   avoided where ever possible.                   |
!                                                  |
!  Special term in the mass (density) equation     | magnetic_calc_density
!  Special term in the momentum (hydro) equation   | magnetic_calc_hydro
!  Special term in the energy equation             | magnetic_calc_energy
!  Special term in the induction (magnetic)        | magnetic_calc_magnetic
!     equation                                     |
!                                                  |
!  Special equation                                | dmagnetic_dt
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of magnetic_dummies.inc) the number of f array
! variables and auxiliary variables added by this module.
! Need 3*2 for Ak, 3*2 for Ek, and then 3+3 for E and B, so 3*6=18.
! However, the space for E and B needs to be allocated in cparam.local.
! Need to allocate 12 basic chunks if full 3-D fields are used,
! or only 4 chunks if lpolarization_basis=T.
!
! CPARAM logical, parameter :: lmagnetic = .true.
! CPARAM logical, parameter :: lbfield = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED bb(3); bbb(3); bij(3,3); jxbr(3); ss12; uxb(3); jj(3)
! PENCILS PROVIDED el(3); e2; a2; b2; bf2
! PENCILS PROVIDED aa(3); diva; del2a(3); aij(3,3); bunit(3); va2
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!   lmagnetic = .true.
! to enable use of magnetic hooks.
!
! The rest of this file may be used as a template for your own
! magnetic module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/magnetic directory.  This path has
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
!   SPECIAL=magnetic/geo_kws
!
! Where geo_kws it replaced by the filename of your new module
! upto and not including the .f90
!
module Magnetic
!
  use Cparam
  use Cdata
  use Initcond
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error, warning
!
  implicit none
!
  include '../magnetic.h'
  !include 'record_types.h'
!
  complex, dimension (:,:,:,:), allocatable :: epol
  real, dimension(3) :: B_ext = 0.0
  real, dimension(3) :: B_ext_inv=(/0.0,0.0,0.0/)
  real, dimension (mz,3) :: aamz
  real :: inertial_length=0.,linertial_2
  real :: ux_const=0., ampl_uy=0.
  logical :: lelectron_inertia=.false.
  logical :: lcalc_aameanz=.false., lcalc_aamean=.false.
  logical :: reinitialize_aa=.false.
  logical, dimension(7) :: lresi_dep=.false. 
  logical :: lcoulomb=.false.
  integer :: pushpars2c, pushdiags2c  ! should be procedure pointer (F2003)
  integer :: iLam=0
!
  integer :: idiag_axmz=0,idiag_aymz=0
  integer :: idiag_bxmz=0,idiag_bymz=0
  integer :: idiag_bsinphz=0,idiag_bcosphz=0
!
! Declare index of new variables in f array (if any).
!
  character (len=labellen) :: conductivity='const'
  character (len=labellen) :: initaak='nothing'
  character (len=labellen) :: initeek='nothing'
  character (len=labellen) :: cc_light='1'
  real :: alpha_inflation=0., beta_inflation=0., eps_quench=.3, nquench=10.
  real :: sigma=0., sigma_t1=0., sigma_t2=0., t1_sigma=0., t2_sigma=0.
  logical :: laa_as_aux=.false., lbb_as_aux=.true., lee_as_aux=.true., ljj_as_aux=.false.
  logical :: lemf=.false., linflation=.false., lreheating=.false., ldebug_print=.false.
  real, dimension(3) :: aakre, aakim, eekre, eekim
  real :: c_light2=1.
  real :: alpha2_inflation, beta1_inflation, beta2_inflation, kscale_factor
  real :: amplaa=1e-4, initpower_aa=0.0, initpower2_aa=-11./3., cutoff_aa=0.0, ncutoff_aa=1.
  real :: kpeak_aa=10., kgaussian_aa=0.
  real :: relhel_aa=1., ksign=1.
  real :: k1hel=0., k2hel=1.
  real :: mu012=0.5 !(=1/2mu0)
  real :: phase_aa=0.
  real :: rescale_aa=0.0
  integer :: kx_aa=0, ky_aa=0, kz_aa=0
  logical :: lpolarization_basis=.false., lswitch_sign_e2=.false., lminus_mode=.false.
  logical :: lscale_tobox=.true., lskip_projection_aa=.false.
  logical :: lalpha_inflation, lbeta_inflation, lquench_inflation=.false.
  logical :: lback_to_unscaled=.true.
!
! input parameters
  namelist /magnetic_init_pars/ &
    linflation, lreheating, alpha_inflation, beta_inflation, &
    lpolarization_basis, lswitch_sign_e2, lminus_mode, initaak, initeek, &
    ux_const, ampl_uy, rescale_aa, &
    laa_as_aux, lbb_as_aux, lee_as_aux, ljj_as_aux, &
    conductivity, sigma, sigma_t1, sigma_t2, t1_sigma, t2_sigma, &
    amplaa, initpower_aa, initpower2_aa, cutoff_aa, ncutoff_aa, kpeak_aa, &
    lscale_tobox, kgaussian_aa, lskip_projection_aa, relhel_aa, &
    lemf, B_ext, &
    k1hel, k2hel, &
    kx_aa, ky_aa, kz_aa, phase_aa
!
! run parameters
  namelist /magnetic_run_pars/ &
    linflation, lreheating, alpha_inflation, beta_inflation, &
    lswitch_sign_e2, lminus_mode, ldebug_print, initaak, initeek, &
    reinitialize_aa, rescale_aa, cc_light, &
    ksign, lemf, B_ext, lback_to_unscaled, &
    conductivity, sigma, sigma_t1, sigma_t2, t1_sigma, t2_sigma, &
    lquench_inflation, eps_quench, nquench
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_aa2m=0      ! DIAG_DOC: $\langle A^2\rangle$
  integer :: idiag_ee2m=0      ! DIAG_DOC: $\langle E^2\rangle$
  integer :: idiag_EEEM=0      ! DIAG_DOC: $\langle(E^2+B^2)/2\rangle$
  integer :: idiag_akxpt=0     ! DIAG_DOC: $Akx^{pt}$
  integer :: idiag_ekxpt=0     ! DIAG_DOC: $Ekx^{pt}$
  integer :: idiag_sigma=0     ! DIAG_DOC: $\sigma$
  integer :: idiag_emag=0      ! DIAG_DOC: $\int_V{1\over2\mu_0}\Bv^2\, dV$
  integer :: idiag_bmax=0      ! DIAG_DOC: $\max(|\Bv|)$
  integer :: idiag_brms=0      ! DIAG_DOC: $\left<\Bv^2\right>^{1/2}$
  integer :: idiag_arms=0      ! DIAG_DOC: $\left<\Av^2\right>^{1/2}$
  integer :: idiag_erms=0      ! DIAG_DOC: $\left<\Ev^2\right>^{1/2}$
  integer :: idiag_bfrms=0     ! DIAG_DOC: $\left<{\Bv'}^2\right>^{1/2}$
!
  integer :: iakx, iaky, iakz, iakxim, iakyim, iakzim
  integer :: iekx, ieky, iekz, iekxim, iekyim, iekzim
  integer, parameter :: nk=nxgrid/2
  type :: magspectra
    real, dimension(nk) :: mag   ,ele
    real, dimension(nk) :: maghel,elehel
  endtype magspectra

  type(magspectra) :: spectra

  contains
!***********************************************************************
    subroutine register_magnetic
!
!  Set up indices for variables in magnetic modules.
!  Need 2x2x3+3=5*3=15 chunks.
!
!  30-mar-21/axel: adapted from gravitational_waves_hTXk.f90
!
      use Sub, only: register_report_aux
      use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Register aak and eek as auxiliary arrays
!  May want to do this only when Fourier transform is enabled.
!  Need 4 chunks (real and imaginary parts of A and E if
!  lpolarization_basis=T and 3*4=12 otherwise for 3 components.
!  In the former case, define iakz=iaky, etc.
!
      if (lpolarization_basis) then
        call farray_register_auxiliary('aak',iaak)
        call farray_register_auxiliary('aakim',iaakim)
        call farray_register_auxiliary('eek',ieek)
        call farray_register_auxiliary('eekim',ieekim)
        iakx  =iaak  ; iakz  =iaak
        iakxim=iaakim; iakzim=iaakim
        iekx  =ieek  ; iekz  =ieek
        iekxim=ieekim; iekzim=ieekim
      else
        call farray_register_auxiliary('aak',iaak,vector=3)
        call farray_register_auxiliary('aakim',iaakim,vector=3)
        call farray_register_auxiliary('eek',ieek,vector=3)
        call farray_register_auxiliary('eekim',ieekim,vector=3)
        iakx  =iaak  ; iaky  =iaak  +1; iakz  =iaak  +2
        iakxim=iaakim; iakyim=iaakim+1; iakzim=iaakim+2 
        iekx  =ieek  ; ieky  =ieek  +1; iekz  =ieek  +2
        iekxim=ieekim; iekyim=ieekim+1; iekzim=ieekim+2 
      endif
!
!  register B array as aux
!
      if (lee_as_aux) call farray_register_auxiliary('ee',iee,vector=3)
      if (lbb_as_aux) call farray_register_auxiliary('bb',ibb,vector=3)
      if (ljj_as_aux) call farray_register_auxiliary('jj',ijj,vector=3)
      if (laa_as_aux) call farray_register_auxiliary('aa',iaa,vector=3)
      iex=iee; iey=iee+1; iez=iee+2
      iax=iaa; iay=iaa+1; iaz=iaa+2
      ibx=ibb; iby=ibb+1; ibz=ibb+2
      ijx=ijj; ijy=ijj+1; ijz=ijj+2
!
    endsubroutine register_magnetic
!***********************************************************************
    subroutine initialize_magnetic(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      use Fourier, only: kx_fft, ky_fft, kz_fft
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (3) :: e1, e2, kvec
      real :: ksqr, k1, k2, k3, k1sqr, k2sqr, k3sqr
      integer :: stat, i, j, ikx, iky, ikz
      complex :: fact=cmplx(0.,-1./sqrt(2.))
!
!  set speed of light
!
      select case (cc_light)
        case ('1'); c_light2=1.
        case ('cgs'); c_light2=c_light_cgs**2
        case default
          call fatal_error("initialize_magnetic","no such cc_light: "//trim(cc_light))
      endselect
      if (headt) print*,'c_light2=',c_light2
!
!  calculate kscale_factor (for later binning)
!
      kscale_factor=2*pi/Lx
!
!  Compute alpha*(alpha+1) from Sharma+17 paper, and 2beta*(2beta+1).
!  Also compute -2*beta. Note the minus sign, so f'/f=beta1_inflation/(t+1).
!
      if (beta_inflation/=0.) then
        beta1_inflation=-2.*beta_inflation
        beta2_inflation=+2.*beta_inflation*(2.*beta_inflation+1.)
        lbeta_inflation=.true.
      endif
!
      if (alpha_inflation/=0.) then
        alpha2_inflation=alpha_inflation*(alpha_inflation+1.)
        lalpha_inflation=.true.
      endif
!
      if (lalpha_inflation.and.lbeta_inflation) &
        call fatal_error("initialize_magnetic","only alpha or beta /= 0")
!
!  Possibility of computing polarization basis here:
!
      if (lpolarization_basis) then
        allocate(epol(nx,ny,nz,3),stat=stat)
        if (stat>0) call fatal_error('initialize_magnetic','could not allocate epol')
        do ikz=1,nz
          do iky=1,ny
            do ikx=1,nx
!
!  Collect k vector and compute k^2 at each point.
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
!  Possibility of swapping the sign of e2 in halfspace of k.
!
              if (lswitch_sign_e2) then
                if (k3<0.) then
                  e2=-e2
                elseif (k3==0.) then
                  if (k2<0.) then
                    e2=-e2
                  elseif (k2==0.) then
                    if (k1<0.) then
                      e2=-e2
                    endif
                  endif
                endif
              endif
!
!  compute epol=(-i/sqrt(2))*(e1+i*e2)
!
              do j=1,3
                epol(ikx,iky,ikz,j)=fact*cmplx(e1(j),e2(j))
              enddo
!
!  end of ikx, iky, and ikz loops
!
            enddo
          enddo
        enddo
      endif
!
!  Possibility of switching the sign of the polarization.
!
      if (lminus_mode) epol=conjg(epol)
!
!  Reinitialize magnetic field using a small selection of perturbations
!  that were mostly also available as initial conditions.
!
      if (reinitialize_aa) then
        select case (initaak)
        case ('rescale')
          f(:,:,:,iaak:iakzim)=rescale_aa*f(:,:,:,iaak:iakzim)
          f(:,:,:,ieek:iekzim)=rescale_aa*f(:,:,:,ieek:iekzim)
        case default
        endselect
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_magnetic
!***********************************************************************
    subroutine init_aa(f)
!
!  initialise magnetic condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Fourier, only: kx_fft, ky_fft, kz_fft
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, parameter :: lvectorpotential=.true.
      real :: k1, k2, k3, ksqr, k
      integer :: ikx, iky, ikz
!
      intent(inout) :: f
!
!  initial condition for aak
!
      f(:,:,:,iaak:iakzim)=0.
      f(:,:,:,ieek:iekzim)=0.
      if (lee_as_aux) f(:,:,:,iee :iee   +2)=0.
      if (laa_as_aux) f(:,:,:,iaa :iaa   +2)=0.
      if (lbb_as_aux) f(:,:,:,ibb :ibb   +2)=0.
      if (ljj_as_aux) f(:,:,:,ijj :ijj   +2)=0.
      select case (initaak)
        case ('nothing')
          if (lroot) print*,'init_magnetic: nothing'
          f(:,:,:,iaak:iakzim)=0.
        case ('single')
          if (lroot) print*,'init_magnetic for A: single'
          f(:,:,:,iaak:iakzim)=0.
          if (lroot) f(l1+1,m1+1,n1+1,iaak)=1.
        case ('Beltrami-general')
               !call beltramik_general(amplaa,f,iaak,kx_aa,ky_aa,kz_aa,phase_aa)
        case ('power_randomphase_hel')
          call power_randomphase_hel(amplaa,initpower_aa,initpower2_aa, &
            cutoff_aa,ncutoff_aa,kpeak_aa,f,iaak,iakz,relhel_aa,kgaussian_aa, &
            lskip_projection_aa, lvectorpotential, &
            lscale_tobox=lscale_tobox, k1hel=k1hel, k2hel=k2hel, &
            lremain_in_fourier=.true.)
        case ('Alfven-x'); call alfvenk_x(amplaa,f,iuu,iaak,kx_aa)
        case default
          call fatal_error("init_magnetic","no such initaak: "//trim(initaak))
      endselect
!
!  initial condition for eek
!
      select case (initeek)
        case ('nothing')
          if (lroot) print*,'init_magnetic: nothing'
          f(:,:,:,ieek:iekzim)=0.
        case ('single')
          if (lroot) print*,'init_magnetic for E: single'
          f(:,:,:,ieek:iekzim)=0.
          if (lroot) f(l1+1,m1+1,n1+1,ieek)=1.
!
!  dA/deta = -E
!  d^2A/deta^2 = -dE/deta = -(k^2-alpha*(alpha+1), or
!  d^2A/deta^2 = -dE/deta = -(k^2-2beta*(2beta+1)
!
        case ('ikA')
          do ikz=1,nz
            do iky=1,ny
              do ikx=1,nx
!
!  collect k vector and compute k^2 at each point
!
                k1=kx_fft(ikx+ipx*nx)
                k2=ky_fft(iky+ipy*ny)
                k3=kz_fft(ikz+ipz*nz)
                ksqr=k1**2+k2**2+k3**2
                k=sqrt(ksqr)
!
!  Compute Ek = -ik*Ak
!  Re Ek = -k*Im(A)
!  Im Ek = +k*Re(A)
!
                f(nghost+ikx,nghost+iky,nghost+ikz,ieek  :iekz)=-k* &
                f(nghost+ikx,nghost+iky,nghost+ikz,iaakim:iakzim)
                f(nghost+ikx,nghost+iky,nghost+ikz,ieekim:iekzim)=+k* &
                f(nghost+ikx,nghost+iky,nghost+ikz,iaak  :iakz)
              enddo
            enddo
          enddo
        case ('powerlaw??>')
          !call powerlaw ...
        case default
          call fatal_error("init_magnetic","no such initeek: "//trim(initeek))
      endselect
!
    endsubroutine init_aa
!***********************************************************************
    subroutine pencil_criteria_magnetic
!
!  All pencils that this magnetic module depends on are specified here.
!
!   1-apr-06/axel: coded
!
      if (idiag_brms/=0 .or. idiag_EEEM/=0 .or. idiag_bfrms/=0 .or. &
          idiag_bmax/=0 .or. idiag_emag/=0 ) &
          lpenc_diagnos(i_b2)=.true.
      if (idiag_erms/=0 .or. idiag_EEEM/=0) lpenc_requested(i_e2)=.true.
      if (idiag_arms/=0) lpenc_diagnos(i_a2)=.true.
      if (idiag_erms/=0) lpenc_diagnos(i_e2)=.true.
      if (lpenc_diagnos(i_a2)) lpenc_requested(i_aa)=.true.
      if (lpenc_diagnos(i_b2)) lpenc_requested(i_bb)=.true.
      if (lpenc_diagnos(i_e2)) lpenc_requested(i_el)=.true.
!
    endsubroutine pencil_criteria_magnetic
!***********************************************************************
    subroutine pencil_interdep_magnetic(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-jul-06/tony: coded
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_magnetic
!***********************************************************************
    subroutine daa_dt(f,df,p)
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
!  This routine computes various diagnostic quantities.
!
!  06-oct-03/tony: coded
!  07-feb-18/axel: added nscale_factor=0 (no expansion), =.5 (radiation era)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(inout) :: f,df
!
!  Identify module and boundary conditions.
!
      if (lfirst) then
        if (headtt.or.ldebug) print*,'daa_dt: SOLVE daa_dt using aak and eek in 1st substep'

        call calc_diagnostics_magnetic(f,p)

      else
        if (headtt.or.ldebug) print*,'daa_dt: DONT SOLVE aa_dt in other substeps'
      endif
!
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(f,df)

    endsubroutine daa_dt
!***********************************************************************
    subroutine calc_diagnostics_magnetic(f,p)
!
      use Diagnostics

      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
!  diagnostics
!
      if (ldiagnos) then
        if (lpolarization_basis) then
          if (idiag_ee2m/=0) call sum_mn_name(( &
              f(l1:l2,m,n,iekx  )**2+ &
              f(l1:l2,m,n,iekxim)**2 &
                                              )*nwgrid,idiag_ee2m)
          if (idiag_aa2m/=0) call sum_mn_name(( &
              f(l1:l2,m,n,iakx  )**2+ &
              f(l1:l2,m,n,iakxim)**2 &
                                              )*nwgrid,idiag_aa2m)
        else
          if (idiag_ee2m/=0) call sum_mn_name(( &
              f(l1:l2,m,n,iekx  )**2+f(l1:l2,m,n,ieky  )**2+f(l1:l2,m,n,iekz  )**2+ &
              f(l1:l2,m,n,iekxim)**2+f(l1:l2,m,n,iekyim)**2+f(l1:l2,m,n,iekzim)**2 &
                                              )*nwgrid,idiag_ee2m)
          if (idiag_aa2m/=0) call sum_mn_name(( &
              f(l1:l2,m,n,iakx  )**2+f(l1:l2,m,n,iaky  )**2+f(l1:l2,m,n,iakz  )**2+ &
              f(l1:l2,m,n,iakxim)**2+f(l1:l2,m,n,iakyim)**2+f(l1:l2,m,n,iakzim)**2 &
                                              )*nwgrid,idiag_aa2m)
        endif

        if (idiag_emag/=0) call integrate_mn_name(mu012*p%b2,idiag_emag)
        call max_mn_name(p%b2,idiag_bmax,lsqrt=.true.)
        call sum_mn_name(p%b2,idiag_brms,lsqrt=.true.)
        call sum_mn_name(p%a2,idiag_arms,lsqrt=.true.)
        call sum_mn_name(p%e2,idiag_erms,lsqrt=.true.)
        if (idiag_EEEM/=0) call sum_mn_name(.5*(p%e2+p%b2),idiag_EEEM)
        call sum_mn_name(p%bf2,idiag_bfrms,lsqrt=.true.)
!
        if (lproc_pt.and.m==mpoint.and.n==npoint) then
          if (idiag_akxpt/=0) call save_name(f(lpoint,m,n,iakx),idiag_akxpt)
          if (idiag_ekxpt/=0) call save_name(f(lpoint,m,n,iekx),idiag_ekxpt)
          if (idiag_sigma/=0) call save_name(sigma,idiag_sigma)
        endif
      endif
!
      call keep_compiler_quiet(p)
!
    endsubroutine calc_diagnostics_magnetic
!***********************************************************************
    subroutine time_integrals_magnetic(f,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine time_integrals_magnetic
!***********************************************************************
    subroutine df_diagnos_magnetic(df,p)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)  ::  df, p
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine df_diagnos_magnetic
!***********************************************************************
    subroutine rescaling_magnetic(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine rescaling_magnetic
!***********************************************************************
    subroutine read_magnetic_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=magnetic_init_pars, IOSTAT=iostat)
!
    endsubroutine read_magnetic_init_pars
!***********************************************************************
    subroutine write_magnetic_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=magnetic_init_pars)
!
    endsubroutine write_magnetic_init_pars
!***********************************************************************
    subroutine read_magnetic_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=magnetic_run_pars, IOSTAT=iostat)
!
    endsubroutine read_magnetic_run_pars
!***********************************************************************
    subroutine write_magnetic_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=magnetic_run_pars)
!
    endsubroutine write_magnetic_run_pars
!***********************************************************************
    subroutine magnetic_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    endsubroutine magnetic_before_boundary
!***********************************************************************
    subroutine calc_pencils_magnetic_std(f,p)
!
!  Standard version (_std): global variable lpencil contains information about needed pencils.
!
      real, dimension (mx,my,mz,mfarray), intent(inout):: f
      type (pencil_case),                 intent(out)  :: p
!
      call calc_pencils_magnetic_pencpar(f,p,lpencil)
!
    endsubroutine calc_pencils_magnetic_std
!***********************************************************************
    subroutine calc_pencils_magnetic_pencpar(f,p,lpenc_loc)
!
!  Calculate Magnetic pencils.
!  Most basic pencils should come first, as others may depend on them.
! 
!  Version with formal parameter lpencil_loc instead of global lpencil for cases
!  in which not necessarily all generally needed pencil are to be calculated.
!
!  19-nov-04/anders: coded
!  18-jun-13/axel: b2 now includes B_ext by default (luse_Bext_in_b2=T is kept)
!  20-jun-16/fred: added derivative tensor option and streamlined gij_etc
!
      use Sub
      use EquationOfState, only: rho0
!
      real, dimension (mx,my,mz,mfarray), intent(inout):: f
      type (pencil_case),                 intent(out)  :: p
      logical, dimension(:),              intent(in)   :: lpenc_loc
!
      real, dimension (nx,3) :: tmp ! currently unused: bb_ext_pot
      real, dimension (nx) :: rho1_jxb, quench, StokesI_ncr, tmp1
  !   real :: c,s
  !   integer :: i,j !,ix
      integer :: j
! ee
      if (lpenc_loc(i_el)) p%el=f(l1:l2,m,n,iex:iez)
! aa
      if (lpenc_loc(i_aa)) p%aa=f(l1:l2,m,n,iax:iaz)
! bb
      if (lpenc_loc(i_bb)) p%bb=f(l1:l2,m,n,ibx:ibz)
! jj
      if (lpenc_loc(i_jj)) p%jj=f(l1:l2,m,n,ijx:ijz)
! bbb
      if (lpenc_loc(i_bbb)) p%bbb=p%bb
! e2
      if (lpenc_loc(i_e2)) call dot2_mn(p%el,p%e2)
! a2
      if (lpenc_loc(i_a2)) call dot2_mn(p%aa,p%a2)
! b2
      if (lpenc_loc(i_b2)) call dot2_mn(p%bb,p%b2)
      if (lpenc_loc(i_bf2)) call dot2_mn(p%bbb,p%bf2)
!
    endsubroutine calc_pencils_magnetic_pencpar
!***********************************************************************
    subroutine update_char_vel_magnetic(f)
!
! Dummy 
!
!  25-sep-15/MR+joern: coded
!
      real, dimension(mx,my,mz,mfarray), intent(INOUT) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine update_char_vel_magnetic
!***********************************************************************
    subroutine magnetic_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  07-aug-17/axel: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f

      if (lfirst) then
!
!  Choice of conductivity profiles.
!
        select case (conductivity)
          case ('const')
            if (headtt.or.ldebug) print*,'sigma=const=',sigma
!
!  Time-dependent profile for sigma.
!
          case ('t-dep')
            if (t<=t1_sigma) then
              sigma=sigma_t1
            elseif (t<=t2_sigma) then
              sigma=sigma_t1+(sigma_t2-sigma_t1)*(t-t1_sigma)/(t2_sigma-t1_sigma)
            else
              sigma=sigma_t2
            endif
!
          case default; call fatal_error('daa_dt','no such conductivity: '//trim(conductivity))
        endselect
      endif
!
    endsubroutine magnetic_after_boundary
!***********************************************************************
    subroutine make_spectra(f)
!
!  31-mar-21/axel: adapted from gravitational_waves_hTXk.f90
!
      use Fourier, only: kx_fft, ky_fft, kz_fft
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: ikx, iky, ikz, q, p, ik
      real :: k1, k2, k3, ksqr
!
      spectra%mag=0.; spectra%maghel=0.
      spectra%ele=0.; spectra%elehel=0.
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
            ik=1+nint(sqrt(ksqr)/kscale_factor)
!
!  Debug output.
!
            if (ldebug_print) then
              if (ik <= 5) write(*,1000) iproc,ik,k1,k2,k3,f(nghost+ikx,nghost+iky,nghost+ikz,iggX  )
              1000 format(2i5,1p,4e11.2)
            endif
!
            if (ik <= nk) then
!
!  Electromagnetic wave spectrum computed from aak.
!
              if (ab_spec) then
                if (lpolarization_basis) then
                  spectra%mag(ik)=spectra%mag(ik) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaak  +0)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaakim+0)**2
                  spectra%maghel(ik)=spectra%maghel(ik) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaak  +0)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaakim+0)**2
                else
                  spectra%mag(ik)=spectra%mag(ik) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaak  +0)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaakim+0)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaak  +1)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaakim+1)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaak  +2)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaakim+2)**2
                  spectra%maghel(ik)=spectra%maghel(ik)+2*( &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaak  +0) &
                    *(f(nghost+ikx,nghost+iky,nghost+ikz,iaakim+1)*(+k3) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaakim+2)*(-k2)) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaak  +1) &
                    *(f(nghost+ikx,nghost+iky,nghost+ikz,iaakim+2)*(+k1) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaakim+0)*(-k3)) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaak  +2) &
                    *(f(nghost+ikx,nghost+iky,nghost+ikz,iaakim+0)*(+k2) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,iaakim+1)*(-k1)) &
                                                          )
                endif
              endif
!
!  Electric field spectrum computed from eek
!
              if (ele_spec) then
                if (lpolarization_basis) then
                  spectra%ele(ik)=spectra%ele(ik) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieek  +0)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieekim+0)**2
                  spectra%elehel(ik)=spectra%elehel(ik) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieek  +0)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieekim+0)**2
                else
                  spectra%ele(ik)=spectra%ele(ik) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieek  +0)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieekim+0)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieek  +1)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieekim+1)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieek  +2)**2 &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieekim+2)**2
                  spectra%elehel(ik)=spectra%elehel(ik)+2*( &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieek  +0) &
                    *(f(nghost+ikx,nghost+iky,nghost+ikz,ieekim+1)*(+k3) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieekim+2)*(-k2)) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieek  +1) &
                    *(f(nghost+ikx,nghost+iky,nghost+ikz,ieekim+2)*(+k1) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieekim+0)*(-k3)) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieek  +2) &
                    *(f(nghost+ikx,nghost+iky,nghost+ikz,ieekim+0)*(+k2) &
                     +f(nghost+ikx,nghost+iky,nghost+ikz,ieekim+1)*(-k1)) &
                                                          )
                endif
              endif
            endif
          enddo
        enddo
      enddo

    endsubroutine make_spectra
!***********************************************************************
    subroutine magnetic_calc_spectra(f,spectrum,spectrum_hel,lfirstcall,kind)
!
!  Calculates magnetic spectra. For use with a single magnetic module.
!
!  16-apr-21/axel: adapted from gravitational_waves_hTXk.f90
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:) :: spectrum,spectrum_hel
      logical :: lfirstcall
      character(LEN=3) :: kind

      if (lfirstcall) then
        call make_spectra(f)
        lfirstcall=.false.
      endif
!
      select case(kind)
      case ('mag'); spectrum=spectra%mag; spectrum_hel=spectra%maghel
      case ('ele'); spectrum=spectra%ele; spectrum_hel=spectra%elehel
      case default; call warning('magnetic_calc_spectra', &
                      'kind of spectrum "'//kind//'" not implemented')
      endselect
!
    endsubroutine magnetic_calc_spectra
!***********************************************************************
    subroutine magnetic_calc_spectra_byte(f,spectrum,spectrum_hel,lfirstcall,kind,len)
!
!  Calculates magnetic and electric spectra. For use with multiple magnetic modules.
!
!  30-mar-21/axel: adapted from gravitational_waves_hTXk.f90
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
      case ('mag'); spectrum=spectra%mag; spectrum_hel=spectra%maghel
      case ('ele'); spectrum=spectra%ele; spectrum_hel=spectra%elehel
      case default; call warning('magnetic_calc_spectra', &
                      'kind of spectrum "'//kindstr//'" not implemented')
      endselect

    endsubroutine magnetic_calc_spectra_byte
!***********************************************************************
    subroutine compute_bb_from_aak_and_eek(f)
!
!  Compute the transverse part of the stress tensor by going into Fourier space.
!  This is currently only called in magnetic_after_timestep, and by that time,
!  some real-space variables haven't been computed yet. Could fix this by
!  calling it also once after start (and from start).
!
!  07-aug-17/axel: coded
!
      use Fourier, only: fft_xyz_parallel, kx_fft, ky_fft, kz_fft
      use Sub, only: cross_mn
!
      real, dimension (:,:,:,:), allocatable :: bbkre, bbkim
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu, bb, uxb
      real, dimension (3) :: coefAre, coefAim, coefBre, coefBim
      integer :: i,j,p,q,ik,ikx,iky,ikz,stat
      complex, dimension (3) :: Acomplex, Ecomplex, Acomplex_new, Ecomplex_new
      complex :: discrim, det1, lam1, lam2, explam1t, explam2t
      complex :: cosotA, cosotE, sinotA, sinotE
      real :: discrim2, sigmaeff, quench
      real :: ksqr, ksqr_eff, k, k1, k2, k3, fact, kdotEMF, ascl, finv
      intent(inout) :: f
!
!  For testing purposes, if lno_transverse_part=T, we would not need to
!  compute the Fourier transform, so we would skip the rest.
!
!  Allocate memory for arrays.
!
      allocate(bbkre(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('compute_bb_from_aak_and_eek','Could not allocate bbkre')
!
      allocate(bbkim(nx,ny,nz,3),stat=stat)
      if (stat>0) call fatal_error('compute_bb_from_aak_and_eek','Could not allocate bbkim')
!
!  Compute electromotive force in real space, and go then to Fourier space.
!  We use the temporary variables bbkre,bbkim for this.
!
      if (lemf) then
        do m=m1,m2
        do n=n1,n2
          uu=f(l1:l2,m,n,iux:iuz)
          bb=f(l1:l2,m,n,ibx:ibz)
          call cross_mn(uu,bb,uxb)
          bbkre(:,m-m1+1,n-n1+1,:)=uxb
          bbkim(:,m-m1+1,n-n1+1,:)=0.
        enddo
        enddo
        call fft_xyz_parallel(bbkre,bbkim)
      endif
!
!  possibility of quenching factor to make magnetogenesis turn off smoothly
!
      if (lquench_inflation) then
        quench=1./(1.+((t+1.)/(2.-eps_quench))**nquench)
      else
        quench=1.
      endif
!
!  Set k^2 array. Note that in Fourier space, kz is the fastest index and has
!  the full nx extent (which, currently, must be equal to nxgrid).
!
      do ikz=1,nz
        do iky=1,ny
          do ikx=1,nx
!
!  collect k vector and compute k^2 at each point
!
            k1=kx_fft(ikx+ipx*nx)
            k2=ky_fft(iky+ipy*ny)
            k3=kz_fft(ikz+ipz*nz)
            ksqr=k1**2+k2**2+k3**2
            k=sqrt(ksqr)
!
!  Define effective squared wavenumber, which can be negative when t is small.
!  With lbeta_inflation, we have f"/f = beta2_inflation/(t+1.)**2.
!
            if (linflation) then
              if (lpolarization_basis) then
                ksqr_eff=ksqr-2.*ksign*k/t
              else
                ksqr_eff=ksqr-alpha2_inflation/t**2
              endif
            elseif (lreheating) then
              if (lalpha_inflation) ksqr_eff=ksqr-alpha2_inflation*2./(t+1.)**2
!
!  With lbeta_inflation, we have f"/f = beta2_inflation/(t+1.)**2, and
!  f'/f = beta1_inflation/(t+1.). Note that f'/f itself is multiplied
!  by another factor of 2. Remember: beta1_inflation = -2*beta.
!  The factor ksign should normally always be =1 to get the fasted growing mode,
!  regardless of whether or not lminus_mode=T or not (default: lminus_mode=F).
!  To reproduce the OF21 paper, however, we put ksign=2.5 (and beta=3).
!
              if (lbeta_inflation) then
                if (lpolarization_basis) then
                  ksqr_eff= ksqr-quench*(beta2_inflation/(t+1.)**2 &
                           -2.*ksign*k*beta1_inflation/(t+1.))
                else
                  ksqr_eff=ksqr-quench*beta2_inflation/(t+1.)**2
                endif
              endif
            else
              ksqr_eff=ksqr
            endif
!
!  Compute eigenvalues; see implementation documentation.
!  Prevent singularity of discrim2 by resetting it to tini if zero.
!
            if (ksqr_eff/=0.) then
              sigmaeff=sigma
              discrim2=sigmaeff**2-4.*ksqr_eff
              if (discrim2==0.) discrim2=tini
              discrim=sqrt(cmplx(discrim2,0.))
              lam1=.5*(-sigmaeff+discrim)
              lam2=.5*(-sigmaeff-discrim)
!
!  Compute exact solution for hT, hX, gT, and gX in Fourier space.
!
              aakre=f(nghost+ikx,nghost+iky,nghost+ikz,iaak  :iakz)
              aakim=f(nghost+ikx,nghost+iky,nghost+ikz,iaakim:iakzim)
!
              eekre=f(nghost+ikx,nghost+iky,nghost+ikz,ieek  :iekz)
              eekim=f(nghost+ikx,nghost+iky,nghost+ikz,ieekim:iekzim)
!
!  Prepare electromotive force EMF, but first project out solenoidal part
!
              if (lemf) then
!
!  Do projection (in 3-D) of (bbkre,bbkim). Use bbk as dummy name.
!  Real part of (ux, uy, uz) -> vx, vy, vz
!  (kk.uu)/k2, vi = ui - ki kj uj
!
                kdotEMF=k1*bbkre(ikx,iky,ikz,1) &
                       +k2*bbkre(ikx,iky,ikz,2) &
                       +k3*bbkre(ikx,iky,ikz,3)
                bbkre(ikx,iky,ikz,1)=bbkre(ikx,iky,ikz,1)-k1*kdotEMF
                bbkre(ikx,iky,ikz,2)=bbkre(ikx,iky,ikz,2)-k2*kdotEMF
                bbkre(ikx,iky,ikz,3)=bbkre(ikx,iky,ikz,3)-k3*kdotEMF
!
!  Do the same for the imaginary part
!
                kdotEMF=k1*bbkim(ikx,iky,ikz,1) &
                       +k2*bbkim(ikx,iky,ikz,2) &
                       +k3*bbkim(ikx,iky,ikz,3)
                bbkim(ikx,iky,ikz,1)=bbkim(ikx,iky,ikz,1)-k1*kdotEMF
                bbkim(ikx,iky,ikz,2)=bbkim(ikx,iky,ikz,2)-k2*kdotEMF
                bbkim(ikx,iky,ikz,3)=bbkim(ikx,iky,ikz,3)-k3*kdotEMF
!
!  Prepare Atilde = A - (sig/k^2)*EMF. Define fact=sig/k^2.
!
                fact=sigmaeff/ksqr_eff
                bbkre(ikx,iky,ikz,:)=fact*bbkre(ikx,iky,ikz,:)
                bbkim(ikx,iky,ikz,:)=fact*bbkim(ikx,iky,ikz,:)
                aakre=aakre-bbkre(ikx,iky,ikz,:)
                aakim=aakim-bbkim(ikx,iky,ikz,:)
              endif
!
              do j=1,3
                Acomplex(j)=cmplx(aakre(j),aakim(j))
                Ecomplex(j)=cmplx(eekre(j),eekim(j))
              enddo
!
!  compute cos(om*dt) and sin(om*dt) to get from one timestep to the next.
!
              explam1t=exp(lam1*dt)
              explam2t=exp(lam2*dt)
              det1=1./discrim
              cosotA=det1*(lam1*explam2t-lam2*explam1t)
              cosotE=det1*(lam1*explam1t-lam2*explam2t)
              sinotA=det1*(     explam2t-     explam1t)
              sinotE=-sinotA*lam1*lam2
!
!  Solve wave equation for hT and gT from one timestep to the next.
!
              Acomplex_new=cosotA*Acomplex+sinotA*Ecomplex
              Ecomplex_new=sinotE*Acomplex+cosotE*Ecomplex
              f(nghost+ikx,nghost+iky,nghost+ikz,iaak  :iaak  +2)= real(Acomplex_new)
              f(nghost+ikx,nghost+iky,nghost+ikz,iaakim:iaakim+2)=aimag(Acomplex_new)
              f(nghost+ikx,nghost+iky,nghost+ikz,ieek  :ieek  +2)= real(Ecomplex_new)
              f(nghost+ikx,nghost+iky,nghost+ikz,ieekim:ieekim+2)=aimag(Ecomplex_new)
!
!  Add back electromotive force to Ak.
!  After this, bbkre and bbrim can be used for other purposes.
!
              if (lemf) then
                f(nghost+ikx,nghost+iky,nghost+ikz,iaak  :iaak  +2)= &
                f(nghost+ikx,nghost+iky,nghost+ikz,iaak  :iaak  +2)+bbkre(ikx,iky,ikz,:)
                f(nghost+ikx,nghost+iky,nghost+ikz,iaakim:iaakim+2)= &
                f(nghost+ikx,nghost+iky,nghost+ikz,iaakim:iaakim+2)+bbkim(ikx,iky,ikz,:)
              endif
!
!  Debug output: check matrix elements in vacuum limit.
!
              if (ldebug_print.and.lroot) then
                if (ikx==2.and.iky==1.and.ikz==1) then
                  print*,'AXEL: Acomplex_new=',Acomplex_new
                  print*,'AXEL: Ecomplex_new=',Ecomplex_new
                  print*,'AXEL: cosotA=',cosotA,cos(sqrt(ksqr_eff)*dt)
                  print*,'AXEL: cosotE=',cosotE,cos(sqrt(ksqr_eff)*dt)
                  print*,'AXEL: sinotA=',sinotA,-sin(sqrt(ksqr_eff)*dt)/sqrt(ksqr_eff)
                  print*,'AXEL: sinotE=',sinotE,+sin(sqrt(ksqr_eff)*dt)*sqrt(ksqr_eff)
                  print*,'AXEL: discrim=',discrim
                endif
              endif
!
            else
!
!  Set origin to zero. It is given by (1,1,1) on root processor.
!
              f(nghost+1,nghost+1,nghost+1,iaak  :iakz  )=0.
              f(nghost+1,nghost+1,nghost+1,iaakim:iakzim)=0.
              f(nghost+1,nghost+1,nghost+1,ieek  :iekz  )=0.
              f(nghost+1,nghost+1,nghost+1,ieekim:iekzim)=0.
            endif
!
!  end of ikx, iky, and ikz loops
!
          enddo
        enddo
      enddo
!
!  compute 1/f if lback_to_unscaled=T
!  Remember that f=a^(-beta), so finv=a^beta (positive exponent!)
!
        if (lback_to_unscaled) then
          ascl=.25*(t+1.)**2
          finv=ascl**beta_inflation
        else
          finv=1.
        endif
!
!  ee back to real space, use the names bbkre and bbkim for ee.
!  In real space, divide by f (if lback_to_unscaled), or multiply by finv.
!  The full expression for E is E=-dA/dt=-d/dt(A/f)=calE+calA*f'/f,
!  where calE=-dcalA/dt, and calA=A/f. Note that f'/f=beta1_inflation/(t+1).
!  Apply the second term later, when calA has been computed in real space.
!  Note also that beta1_inflation=-2*beta has been defined with minus sign.
!  Also, since the stress has now an f^2 factor, we undo the finv multiplication
!  and just keep the inner derivative term on E, so E => E+A*f'/f, and B => B.
!
      if (lee_as_aux) then
        if (lpolarization_basis) then
          do j=1,3
            bbkre(:,:,:,j)=f(l1:l2,m1:m2,n1:n2,ieek  )*real(epol(:,:,:,j)) &
                          -f(l1:l2,m1:m2,n1:n2,ieekim)*aimag(epol(:,:,:,j))
            bbkim(:,:,:,j)=f(l1:l2,m1:m2,n1:n2,ieekim)*real(epol(:,:,:,j)) &
                          +f(l1:l2,m1:m2,n1:n2,ieek  )*aimag(epol(:,:,:,j))
          enddo
        else
          bbkre=f(l1:l2,m1:m2,n1:n2,ieek  :ieek  +2)
          bbkim=f(l1:l2,m1:m2,n1:n2,ieekim:ieekim+2)
        endif
        call fft_xyz_parallel(bbkre,bbkim,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,iee:iee+2)=bbkre
      endif
!
!  aa back to real space, use the names bbkre and bbkim for aa.
!
      if (laa_as_aux) then
        if (lpolarization_basis) then
          do j=1,3
            bbkre(:,:,:,j)=f(l1:l2,m1:m2,n1:n2,iaak  )*real(epol(:,:,:,j)) &
                          -f(l1:l2,m1:m2,n1:n2,iaakim)*aimag(epol(:,:,:,j))
            bbkim(:,:,:,j)=f(l1:l2,m1:m2,n1:n2,iaakim)*real(epol(:,:,:,j)) &
                          +f(l1:l2,m1:m2,n1:n2,iaak  )*aimag(epol(:,:,:,j))
          enddo
        else
          bbkre=f(l1:l2,m1:m2,n1:n2,iaak  :iaak  +2)
          bbkim=f(l1:l2,m1:m2,n1:n2,iaakim:iaakim+2)
        endif
        call fft_xyz_parallel(bbkre,bbkim,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,iaa:iaa+2)=bbkre
        if (lback_to_unscaled) then
          if (lee_as_aux) then
            f(l1:l2,m1:m2,n1:n2,iee:iee+2)=f(l1:l2,m1:m2,n1:n2,iee:iee+2) &
                 +(beta1_inflation/(t+1.))*f(l1:l2,m1:m2,n1:n2,iaa:iaa+2)
          endif
        endif
      endif
!
!  jj back to real space, use the names bbkre and bbkim for jj.
!
      if (ljj_as_aux) then
        do ikz=1,nz
          do iky=1,ny
            do ikx=1,nx
              k1=kx_fft(ikx+ipx*nx)
              k2=ky_fft(iky+ipy*ny)
              k3=kz_fft(ikz+ipz*nz)
              ksqr=k1**2+k2**2+k3**2
              if (lpolarization_basis) then
                do j=1,3
                  bbkre(ikx,iky,ikz,j)=ksqr*( &
                    f(nghost+ikx,nghost+iky,nghost+ikz,iaak  )*real(epol(ikx,iky,ikz,j)) &
                   -f(nghost+ikx,nghost+iky,nghost+ikz,iaakim)*aimag(epol(ikx,iky,ikz,j)))
                  bbkim(ikx,iky,ikz,j)=ksqr*( &
                    f(nghost+ikx,nghost+iky,nghost+ikz,iaakim)*real(epol(ikx,iky,ikz,j)) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iaak  )*aimag(epol(ikx,iky,ikz,j)))
                enddo
              else
                bbkre(ikx,iky,ikz,:)=ksqr*f(nghost+ikx,nghost+iky,nghost+ikz,iaak  :iaak  +2)
                bbkim(ikx,iky,ikz,:)=ksqr*f(nghost+ikx,nghost+iky,nghost+ikz,iaakim:iaakim+2)
              endif
            enddo
          enddo
        enddo
        call fft_xyz_parallel(bbkre,bbkim,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ijj:ijj+2)=bbkre
      endif
!
!  Back to real space, but first compute B in Fourier space as
!  Bk = ik x Ak, so Bk'+iBk'' = ik*(Ak'+iAk'') = ik*Ak'-k*Ak''.
!  Therefore Re(Bk) = -k x Im(Ak), Im(Bk) = +k x Re(Ak)
!
      if (lbb_as_aux) then
        do ikz=1,nz
          do iky=1,ny
            do ikx=1,nx
              k1=kx_fft(ikx+ipx*nx)
              k2=ky_fft(iky+ipy*ny)
              k3=kz_fft(ikz+ipz*nz)
              if (lpolarization_basis) then
                k=sqrt(k1**2+k2**2+k3**2)
                do j=1,3
                  bbkre(ikx,iky,ikz,j)=+k*( &
                    f(nghost+ikx,nghost+iky,nghost+ikz,iaak  )*real(epol(ikx,iky,ikz,j)) &
                   -f(nghost+ikx,nghost+iky,nghost+ikz,iaakim)*aimag(epol(ikx,iky,ikz,j)))
                  bbkim(ikx,iky,ikz,j)=+k*( &
                    f(nghost+ikx,nghost+iky,nghost+ikz,iaakim)*real(epol(ikx,iky,ikz,j)) &
                   +f(nghost+ikx,nghost+iky,nghost+ikz,iaak  )*aimag(epol(ikx,iky,ikz,j)))
                enddo
              else
!
!  standard method; non-polarization based:
!
                aakre=f(nghost+ikx,nghost+iky,nghost+ikz,iaak  :iaak  +2)
                aakim=f(nghost+ikx,nghost+iky,nghost+ikz,iaakim:iaakim+2)
!
!  Re(Bk) = -k x Im(Ak)
!
                bbkre(ikx,iky,ikz,1)=-k2*aakim(3)+k3*aakim(2)
                bbkre(ikx,iky,ikz,2)=-k3*aakim(1)+k1*aakim(3)
                bbkre(ikx,iky,ikz,3)=-k1*aakim(2)+k2*aakim(1)
!
!  Im(Bk) = +k x Re(Ak)
!
                bbkim(ikx,iky,ikz,1)=+k2*aakre(3)+k3*aakre(2)
                bbkim(ikx,iky,ikz,2)=+k3*aakre(1)+k1*aakre(3)
                bbkim(ikx,iky,ikz,3)=+k1*aakre(2)+k2*aakre(1)
              endif
            enddo
          enddo
        enddo
        call fft_xyz_parallel(bbkre,bbkim,linv=.true.)
        f(l1:l2,m1:m2,n1:n2,ibb:ibb+2)=bbkre
!
!  Add external (imposed) field B_ext, if nonvanishing.
!
        if (any(B_ext/=0.)) then
          forall(j = 1:3, B_ext(j) /= 0.0) &
              f(l1:l2,m1:m2,n1:n2,ibb-1+j)=f(l1:l2,m1:m2,n1:n2,ibb-1+j)+B_ext(j)
          if (headtt) print *, 'calc_pencils_magnetic_pencpar: B_ext = ', B_ext
        endif
      endif
!
    endsubroutine compute_bb_from_aak_and_eek
!***********************************************************************
    subroutine rprint_magnetic(lreset,lwrite)
!
!  Reads and registers print parameters relevant to magnetic.
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
        idiag_ee2m=0; idiag_aa2m=0
        idiag_akxpt=0
        idiag_ekxpt=0
        idiag_emag=0
        idiag_bmax=0; idiag_brms=0; idiag_arms=0; idiag_erms=0
        idiag_EEEM=0; idiag_bfrms=0
        idiag_sigma=0
        cformv=''
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ee2m',idiag_ee2m)
        call parse_name(iname,cname(iname),cform(iname),'aa2m',idiag_aa2m)
        call parse_name(iname,cname(iname),cform(iname),'akxpt',idiag_akxpt)
        call parse_name(iname,cname(iname),cform(iname),'ekxpt',idiag_ekxpt)
        call parse_name(iname,cname(iname),cform(iname),'emag',idiag_emag)
        call parse_name(iname,cname(iname),cform(iname),'bmax',idiag_bmax)
        call parse_name(iname,cname(iname),cform(iname),'brms',idiag_brms)
        call parse_name(iname,cname(iname),cform(iname),'arms',idiag_arms)
        call parse_name(iname,cname(iname),cform(iname),'erms',idiag_erms)
        call parse_name(iname,cname(iname),cform(iname),'EEEM',idiag_EEEM)
        call parse_name(iname,cname(iname),cform(iname),'bfrms',idiag_bfrms)
        call parse_name(iname,cname(iname),cform(iname),'sigma',idiag_sigma)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='bb'.or.cnamev=='ee'.or.cnamev=='jj'.or. &
              cnamev=='h22'.or.cnamev=='h33'.or.cnamev=='h23') cformv='DEFINED'
      endif
!
!!!  write column where which magnetic variable is stored
!!      if (lwrite) then
!!        call farray_index_append('i_SPECIAL_DIAGNOSTIC',i_SPECIAL_DIAGNOSTIC)
!!      endif
!!
    endsubroutine rprint_magnetic
!***********************************************************************
    subroutine get_slices_magnetic(f,slices)
!
!  Write slices for animation of magnetic variables.
!
!  11-apr-21/axel: adapted from gravitational_waves_hTXk.f90
!
      use Slices_methods, only: assign_slices_vec
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  bb
!
        case ('bb')
          if (lbb_as_aux) call assign_slices_vec(slices,f,ibb)
!
      endselect
!
    endsubroutine get_slices_magnetic
!***********************************************************************
    subroutine bdry_magnetic(f,quench,task)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (nx) :: quench
      character (len=*), intent(in) :: task

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(quench)
      call keep_compiler_quiet(task)
!
      call fatal_error('bdry_magnetic','not to be called w/o B-field')
!
    endsubroutine bdry_magnetic
!***********************************************************************
    subroutine calc_mfield
!
!  Dummy routine
!
    endsubroutine calc_mfield
!***********************************************************************
    subroutine bb_unitvec_shock(f,bb_hat)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      real, dimension (mx,3), intent (out) :: bb_hat
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bb_hat)
!
    endsubroutine bb_unitvec_shock
!***********************************************************************
    subroutine input_persistent_magnetic(id,done)
!
!  Dummy routine
!
      integer, optional :: id
      logical, optional :: done
!
      if (present (id)) call keep_compiler_quiet(id)
      if (present (done)) call keep_compiler_quiet(done)
!
    endsubroutine input_persistent_magnetic
!***********************************************************************
    logical function output_persistent_magnetic()
!
!  Dummy routine
!
      output_persistent_magnetic = .false.
!
    endfunction output_persistent_magnetic
!***********************************************************************
    subroutine dynamical_resistivity(uc)
!
!  dummy
!
      real, intent(in) :: uc
!
      call keep_compiler_quiet(uc)
!
    endsubroutine dynamical_resistivity
!***********************************************************************
    subroutine split_update_magnetic(f)
!
!  dummy
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine split_update_magnetic
!***********************************************************************
    subroutine magnetic_after_timestep(f,df,dtsub)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dtsub
!
!  Compute the transverse part of the stress tensor by going into Fourier space.
!
      if (lfirst) call compute_bb_from_aak_and_eek(f)
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dtsub)
!
    endsubroutine magnetic_after_timestep
!***********************************************************************
    subroutine expand_shands_magnetic
!
!  Dummy
!
    endsubroutine expand_shands_magnetic
!***********************************************************************
    subroutine get_bext(B_ext_out)
!
!  Dummy
!
      real, dimension(3), intent(out) :: B_ext_out
!
      B_ext_out = 0.0
!
    endsubroutine get_bext
!***********************************************************************
    subroutine alfvenk_x(ampl,f,iuu,iaak,kx)
!
!  Alfven wave propagating in the x-direction
!
!  uy = +sink(x-vA*t)
!  Az = -cosk(x-vA*t)*sqrt(rho*mu0)/k
!
!  Alfven and slow magnetosonic are the same here and both incompressible, and
!  a fast magnetosonic (compressible) wave is also excited, but decoupled.
!
!  satisfies the four equations
!  duy/dt = B0*By'  ==>  duy/dt = -B0*Az''
!  dBy/dt = B0*uy'  ==>  dAz/dt = -B0*ux
!
!   8-nov-03/axel: coded
!  29-apr-03/axel: added sqrt(rho*mu0)/k factor
!   7-aug-17/axel: added sqrt(.75) for lrelativistic_eos=T
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rho,ampl_Az
      real :: ampl
      integer :: iuu,iaak,kx
!
!  Amplitude factors
!
   !  ampl_uy=+ampl
!
!  ux and Ay.
!  Don't overwrite the density, just add to the log of it.
!
      do n=n1,n2; do m=m1,m2
        f(l1:l2,m,n,iuu+0)=ux_const
        f(l1:l2,m,n,iuu+1)=ampl_uy*sin(kx*x(l1:l2))
      enddo; enddo
!
!  preparations; j is the first index of the array for the complex part.
!
!
      if (kx>0.and.kx<=nx-1) then
        if (ipx==0) then
          f(l1+kx,m1,n1,iaak+2)=+.5*ampl
        endif
!
        if (ipx==nprocx-1) then
          f(l2+1-kx,m1,n1,iaak+2)=+.5*ampl
        endif
      endif
      if (lroot) print*,'alfven_x: kx, ampl=',kx, ampl
!
    endsubroutine alfvenk_x
!***********************************************************************
    real function beltrami_phase()
!
!  Dummy routine
!
      use Mpicomm, only: mpibcast_real
!
    endfunction beltrami_phase
!***********************************************************************
endmodule Magnetic
