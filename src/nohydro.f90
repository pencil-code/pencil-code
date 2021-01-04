! $Id$
!
! MODULE_DOC: no variable $\uv$: useful for kinematic dynamo runs.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lhydro = .false.
! CPARAM logical, parameter :: lhydro_kinematic = .false.
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED uu(3); u2; oo(3); ou;  oxu(3); uij(3,3); sij(3,3); sij2
! PENCILS PROVIDED divu; uij5(3,3); graddivu(3); ugu(3); ogu(3)
! PENCILS PROVIDED del2u(3), curlo(3), uu_advec(3)
!
!***************************************************************
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
  real, dimension (mz,3) :: uumzg=0.
  real, dimension (nz,3) :: guumz=0.
  real, dimension (mx,my,3) :: uumxy=0.
  real, dimension (mx,mz,3) :: uumxz=0.
!
  real :: u_out_kep=0.0
  logical, target :: lpressuregradient_gas=.false.
  logical :: lcalc_uumean=.false.,lupw_uu=.false.
  logical :: lcalc_uumeanx=.false.,lcalc_uumeanz=.false.
  logical :: lcalc_uumeanxy=.false.,lcalc_uumeanxz=.false.
!
  real, allocatable, dimension (:,:) :: KS_k,KS_A,KS_B !or through whole field for each wavenumber?
  real, allocatable, dimension (:) :: KS_omega !or through whole field for each wavenumber?
  integer :: KS_modes = 3
  real, allocatable, dimension (:) :: Zl,dZldr,Pl,dPldtheta
  real :: ampl_fcont_uu=1., w_sld_cs=0.0
  logical :: lforcing_cont_uu=.false.
  integer :: pushpars2c, pushdiags2c  ! should be procedure pointer (F2003)
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
  integer :: idiag_ekintot=0, idiag_ekin=0
!
  contains
!***********************************************************************
    subroutine register_hydro
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use Mpicomm, only: lroot
      use SharedVariables, only: put_shared_variable
!
      integer :: ierr
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Share lpressuregradient_gas so Entropy module knows whether to apply
!  pressure gradient or not.
!
      call put_shared_variable('lpressuregradient_gas',lpressuregradient_gas,ierr)
      if (ierr/=0) call fatal_error('register_hydro','there was a problem sharing lpressuregradient_gas')
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine initialize_hydro(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (kinflow=='KS') then
!        call random_isotropic_KS_setup(-5./3.,1.,(nxgrid)/2.)
!
!  Use constant values for testing KS model code with 3
!  specific modes.
!
        call random_isotropic_KS_setup_test
      endif
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
      if (lkinflow_as_aux) then
        if (iuu==0) then
          call farray_register_auxiliary('uu',iuu,vector=3)
          iux=iuu
          iuy=iuu+1
          iuz=iuu+2
        else
          if (lroot .and. (ip<14)) print*, 'initialize_hydro: iuu = ', iuu
          call farray_index_append('iuu',iuu,3)
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_hydro
!***********************************************************************
      subroutine calc_means_hydro(f)
!
!  dummy routine
!
!  14-oct-13/MR: coded
!
        real, dimension (mx,my,mz,mfarray), intent(IN) :: f
        call keep_compiler_quiet(f)
!
      endsubroutine calc_means_hydro
!***********************************************************************
    subroutine init_uu(f)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Hydro module, if there was one.
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
      if (kinflow/='') then
        lpenc_requested(i_uu)=.true.
        if (kinflow=='eddy') then
          lpenc_requested(i_rcyl_mn)=.true.
          lpenc_requested(i_rcyl_mn1)=.true.
        endif
      endif
!
!  disgnostic pencils
!
      if (idiag_urms/=0 .or. idiag_umax/=0 .or. idiag_u2m/=0 .or. &
          idiag_um2/=0) lpenc_diagnos(i_u2)=.true.
      if (idiag_oum/=0) lpenc_diagnos(i_ou)=.true.
!
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
!ajwm May be overkill... Perhaps only needed for certain kinflow?
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
      if (lpencil_in(i_ou) .or. lpencil_in(i_oxu)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_oo)=.true.
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
!
      use Diagnostics, only: sum_mn_name, max_mn_name, integrate_mn_name
      use General, only: random_number_wrapper
      use Sub, only: quintic_step, quintic_der_step, dot_mn, dot2_mn
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(:) :: lpenc_loc
      integer :: kk
!
      intent(in) :: f, lpenc_loc
      intent(inout) :: p
! uu
      if (lpenc_loc(i_uu)) p%uu=0.0
! u2
      if (lpenc_loc(i_u2)) p%u2=0.0
! oo
      if (lpenc_loc(i_oo)) p%oo=0.0
! ou
      if (lpenc_loc(i_ou)) p%ou=0.0
! oxu
      if (lpenc_loc(i_oxu)) p%oxu=0.0
! uij
      if (lpenc_loc(i_uij)) p%uij=0.0
! sij
      if (lpenc_loc(i_sij)) p%sij=0.0
! sij2
      if (lpenc_loc(i_sij2)) p%sij2=0.0
! divu
      if (lpenc_loc(i_divu)) p%divu=0.0
! uij5
      if (lpenc_loc(i_uij5)) p%uij5=0.0
! graddivu
      if (lpenc_loc(i_graddivu)) p%graddivu=0.0
! ugu
      if (lpenc_loc(i_ugu)) p%ugu=0.0
! ogu
      if (lpenc_loc(i_ogu)) p%ogu=0.0
! del2u
      if (lpenc_loc(i_del2u)) p%del2u=0.0
! curlo
      if (lpenc_loc(i_curlo)) p%curlo=0.0
!
!  Calculate maxima and rms values for diagnostic purposes
!
      if (ldiagnos) then
        if (idiag_urms/=0)  call sum_mn_name(p%u2,idiag_urms,lsqrt=.true.)
        if (idiag_umax/=0)  call max_mn_name(p%u2,idiag_umax,lsqrt=.true.)
        if (idiag_uzrms/=0) &
            call sum_mn_name(p%uu(:,3)**2,idiag_uzrms,lsqrt=.true.)
        if (idiag_uzmax/=0) &
            call max_mn_name(p%uu(:,3)**2,idiag_uzmax,lsqrt=.true.)
        if (idiag_u2m/=0)   call sum_mn_name(p%u2,idiag_u2m)
        if (idiag_um2/=0)   call max_mn_name(p%u2,idiag_um2)
!
        if (idiag_ekin/=0)  call sum_mn_name(.5*p%rho*p%u2,idiag_ekin)
        if (idiag_ekintot/=0) &
            call integrate_mn_name(.5*p%rho*p%u2,idiag_ekintot)
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_hydro_pencpar
!***********************************************************************
    subroutine hydro_before_boundary(f)
!
!  Actions to take before boundary conditions are set, dummy routine.
!
!   17-dec-2010/Bourdin.KIS: coded
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
      use Diagnostics, only: sum_mn_name, save_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)  :: df,p
      intent(out) :: f
!
!  Calculate maxima and rms values for diagnostic purposes
!
      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'duu_dt: diagnostics ...'
        if (idiag_oum/=0) call sum_mn_name(p%ou,idiag_oum)
        !if (idiag_orms/=0) call sum_mn_name(p%o2,idiag_orms,lsqrt=.true.)
        !if (idiag_omax/=0) call max_mn_name(p%o2,idiag_omax,lsqrt=.true.)
!
!  kinetic field components at one point (=pt)
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_uxpt/=0) call save_name(p%uu(lpoint-nghost,1),idiag_uxpt)
          if (idiag_uypt/=0) call save_name(p%uu(lpoint-nghost,2),idiag_uypt)
          if (idiag_uzpt/=0) call save_name(p%uu(lpoint-nghost,3),idiag_uzpt)
        endif
      endif
!
      call keep_compiler_quiet(f,df)
!
    endsubroutine duu_dt
!***********************************************************************
    subroutine calc_diagnostics_hydro(f,p)

      real, dimension(:,:,:,:) :: f
      type(pencil_case), intent(in) :: p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)

    endsubroutine calc_diagnostics_hydro
!******************************************************************************
    subroutine df_diagnos_hydro(df,p)

      use Diagnostics, only: sum_mn_name

      type(pencil_case), intent(in) :: p
      real, dimension(:,:,:,:) :: df

      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)

    endsubroutine df_diagnos_hydro
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
!  coriolis terms for cartesian geometry
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
!  dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
!    Slope limited diffusion: update characteristic speed
!    Not staggered yet
!
     if (lslope_limit_diff .and. llast) then
           f(:,:,:,isld_char)=0.
     endif
!
    endsubroutine hydro_after_boundary
!***********************************************************************
    subroutine random_isotropic_KS_setup_tony(initpower,kmin,kmax)
!
!   produces random, isotropic field from energy spectrum following the
!   KS method (Malik and Vassilicos, 1999.)
!
!   more to do; unsatisfactory so far - at least for a steep power-law
!   energy spectrum
!
!   27-may-05/tony: modified from snod's KS hydro initial
!   03-feb-06/weezy: Tony's code doesn't appear to have the
!                    correct periodicity.
!                    renamed from random_isotropic_KS_setup
!
    use Sub, only: cross
    use General, only: random_number_wrapper
!
    integer :: modeN
!
    real, dimension (3) :: k_unit
    real, dimension (3) :: e1,e2
    real, dimension (6) :: r
    real, dimension (3) ::j,l  !get rid of this - these replace ee,ee1
    real :: initpower,kmin,kmax
    real, dimension (KS_modes) :: k,dk,energy,ps
    real :: theta,phi,alpha,beta
    real :: a,mkunit
    real :: newthet,newphi  !get rid of this line if there's no change
!
    allocate(KS_k(3,KS_modes))
    allocate(KS_A(3,KS_modes))
    allocate(KS_B(3,KS_modes))
    allocate(KS_omega(KS_modes))
!
!    minlen=Lxyz(1)/(nx-1)
!    kmax=2.*pi/minlen
!    KS_modes=int(0.5*(nx-1))
!    hh=Lxyz(1)/(nx-1)
!    pta=(nx)**(1.0/(nx-1))
!    do modeN=1,KS_modes
!       ggt=(kkmax-kkmin)/(KS_modes-1)
!       ggt=(kkmax/kkmin)**(1./(KS_modes-1))
!        k(modeN)=kmin+(ggt*(modeN-1))
!        k(modeN)=(modeN+3)*2*pi/Lxyz(1)
!       k(modeN)=kkmin*(ggt**(modeN-1)
!    enddo
!
!    do modeN=1,KS_modes
!       if (modeN==1)delk(modeN)=(k(modeN+1)-K(modeN))
!       if (modeN==KS_modes)delk(modeN)=(k(modeN)-k(modeN-1))
!       if (modeN>1.and.modeN<KS_modes)delk(modeN)=(k(modeN+1)-k(modeN-2))/2.0
!    enddo
!          mk=(k2*k2)*((1.0 + (k2/(bk_min*bk_min)))**(0.5*initpower-2.0))
!
!  set kmin
!
       kmin=2.*pi      !/(1.0*Lxyz(1))
!       kmin=kmin*2.*pi
       kmax=128.*pi    !nx*pi
       a=(kmax/kmin)**(1./(KS_modes-1.))
!
!
    do modeN=1,KS_modes
!
!  pick wavenumber
!
!       k=modeN*kmin
      k=kmin*(a**(modeN-1.))
!
!  calculate dk
!
!       print *,kmin,kmax,k
!       dk=1.0*kmin
!
      if (modeN==1)&
              dk=kmin*(a-1.)/2.
      if (modeN>1.and.modeN<KS_modes) &
              dk=(a**(modeN-2.))*kmin*((a**2.) -1.)/2.
      if (modeN==KS_modes) &
              dk=(a**(KS_modes -2.))*kmin*(a -1.)/2.
!
       call random_number_wrapper(r)
       theta=r(1)*pi
       phi=r(2)*2.0*pi
       alpha=r(3)*pi
       beta=r(4)*2.0*pi
       newthet=r(5)*pi
       newphi=r(6)*2.0*pi
!
       k_unit(1)=sin(theta)*cos(phi)
       k_unit(2)=sin(theta)*sin(phi)
       k_unit(3)=cos(theta)
!
       j(1)=sin(alpha)*cos(beta)
       j(2)=sin(alpha)*sin(beta)
       j(3)=cos(alpha)
!
       l(1)=sin(newthet)*cos(newphi)
       l(2)=sin(newthet)*sin(newphi)
       l(3)=cos(newthet)
!
       KS_k(:,modeN)=k*k_unit(:)
!
       call cross(KS_k(:,modeN),j,e1)
       call cross(KS_k(:,modeN),l,e2)
!
!  Make e1 & e2 unit vectors so that we can later make them
!  the correct lengths
!
       mkunit=sqrt(e1(1)**2+e1(2)**2+e1(3)**2)
       e1=e1/mkunit
!
       mkunit=sqrt(e2(1)**2+e2(2)**2+e2(3)**2)
       e2=e2/mkunit
!
!        energy=(((k/1.)**2. +1.)**(-11./6.))*(k**2.) &
!                            *exp(-0.5*(k/kmax)**2.)
!  The energy above is how this code has it. i
!  I've changed the divisor of k.
       energy=(((k/kmin)**2. +1.)**(-11./6.))*(k**2.) &
                       *exp(-0.5*(k/kmax)**2.)
       energy=1.*energy
       ps=sqrt(2.*energy*dk)   !/3.0)
!
       KS_A(:,modeN) = ps*e1
       KS_B(:,modeN) = ps*e2
!
    enddo
!
!   form RA = RA x k_unit and RB = RB x k_unit
!
    do modeN=1,KS_modes
      call cross(KS_A(:,modeN),k_unit(:),KS_A(:,modeN))
      call cross(KS_B(:,modeN),k_unit(:),KS_B(:,modeN))
    enddo
!
    call keep_compiler_quiet(initpower)
!
    endsubroutine random_isotropic_KS_setup_tony
!***********************************************************************
    subroutine random_isotropic_KS_setup(initpower,kmin,kmax)
!
!   produces random, isotropic field from energy spectrum following the
!   KS method (Malik and Vassilicos, 1999.)
!
!   more to do; unsatisfactory so far - at least for a steep power-law
!   energy spectrum
!
!   27-may-05/tony: modified from snod's KS hydro initial
!   03-feb-06/weezy: Attempted rewrite to guarantee periodicity of
!                    KS modes.
!
    use Sub, only: cross, dot2
    use General, only: random_number_wrapper
!
    integer :: modeN
!
    real, dimension (3) :: k_unit
    real, dimension (3) :: ee,e1,e2
!    real, dimension (4) :: r
    real, dimension (6) :: r
    real :: initpower,kmin,kmax
    real, dimension (KS_modes) :: k,dk,energy,ps
    real :: theta,phi,alpha,beta
    real :: ex,ey,ez,norm,a
!
    allocate(KS_k(3,KS_modes))
    allocate(KS_A(3,KS_modes))
    allocate(KS_B(3,KS_modes))
    allocate(KS_omega(KS_modes))
!
    kmin=2.*pi      !/(1.0*Lxyz(1))
!    kmin=kmin*2.*pi
    kmax=128.*pi    !nx*pi
    a=(kmax/kmin)**(1./(KS_modes-1.))
!
!
    do modeN=1,KS_modes
!
!  pick wavenumber
!
!      k=modeN*kmin
      k=kmin*(a**(modeN-1.))
!
!weezy need to investigate if this is still needed
!weezy !
!weezy !  calculate dk
!weezy !
!weezy       print *,kmin,kmax,k
!weezy       dk=1.0*kmin
!
!weezy       if (modeN==1)dk=kmin*(a-1.)/2.
!weezy       if (modeN>1.and.modeN<KS_modes)dk=(a**(modeN-2.))*kmin*((a**2.) -1.)/2.
!weezy       if (modeN==KS_modes)dk=(a**(KS_modes -2.))*kmin*(a -1.)/2.
!
!
!  pick 4 random angles for each mode
!
      call random_number_wrapper(r);
      theta=pi*(r(1) - 0.)
      phi=pi*(2*r(2) - 0.)
      alpha=pi*(2*r(3) - 0.)
      beta=pi*(2*r(4) - 0.)
!
!  random phase?
!      call random_number_wrapper(r); gamma(modeN)=pi*(2*r - 0.)
!
!  make a random unit vector by rotating fixed vector to random position
!  (alternatively make a random transformation matrix for each k)
!
      k_unit(1)=sin(theta)*cos(phi)
      k_unit(2)=sin(theta)*sin(phi)
      k_unit(3)=cos(theta)
!
      energy=(((k/kmin)**2. +1.)**(-11./6.))*(k**2.) &
                       *exp(-0.5*(k/kmax)**2.)
!      energy=(((k/1.)**2. +1.)**(-11./6.))*(k**2.) &
!                       *exp(-0.5*(k/kmax)**2.)
!
!  make a vector KS_k of length k from the unit vector for each mode
!
      KS_k(:,modeN)=k*k_unit(:)
!      KS_omega(modeN)=k**(2./3.)
      KS_omega(:)=sqrt(energy(:)*(k(:)**3.))
!
!  construct basis for plane having rr normal to it
!  (bit of code from forcing to construct x', y')
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
!  make two random unit vectors KS_B and KS_A in the constructed plane
!
      KS_A(:,modeN) = cos(alpha)*e1 + sin(alpha)*e2
      KS_B(:,modeN) = cos(beta)*e1  + sin(beta)*e2
!
!  define the power spectrum (ps=sqrt(2.*power_spectrum(k)*delta_k/3.))
!
!      ps=(k**(initpower/2.))*sqrt(dk*2./3.)
!  The factor of 2 just after the sqrt may need to be 2./3.
!
!
!  With the `weezey' stuff above commented out, dk is currently used, but
!  never set, so we better abort
!
      call error('random_isotropic_KS_setup', 'Using uninitialized dk')
      dk=0.                     ! to make compiler happy
!
      ps=sqrt(2.*energy*dk)   !/3.0)
!
!  give KS_A and KS_B length ps
!
      KS_A(:,modeN)=ps*KS_A(:,modeN)
      KS_B(:,modeN)=ps*KS_B(:,modeN)
!
    enddo
!
!  form RA = RA x k_unit and RB = RB x k_unit
!  Note: cannot reuse same vector for input and output
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
!   produces random, isotropic field from energy spectrum following the
!   KS method (Malik and Vassilicos, 1999.)
!   This test case only uses 3 very specific modes (useful for comparison
!   with Louise's kinematic dynamo code.
!
!   03-feb-06/weezy: modified from random_isotropic_KS_setup
!
    use Sub, only: cross
!
    integer :: modeN
!
    real, dimension (3,KS_modes) :: k_unit
    real, dimension (KS_modes) :: k,dk,energy,ps
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
!-----------------------------
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
!-----------------------------
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
!   form RA = RA x k_unit and RB = RB x k_unit
!
!
     do modeN=1,KS_modes
       call cross(KS_A(:,modeN),k_unit(:,modeN),KS_A(:,modeN))
       call cross(KS_B(:,modeN),k_unit(:,modeN),KS_B(:,modeN))
     enddo
!
    endsubroutine random_isotropic_KS_setup_test
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
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_hydro_run_pars
!***********************************************************************
    subroutine write_hydro_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_hydro_run_pars
!***********************************************************************
    subroutine input_persistent_hydro(id,done)
!
      integer, optional :: id
      logical, optional :: done
!
      if (present (id)) call keep_compiler_quiet(id)
      if (present (done)) call keep_compiler_quiet(done)
!
    endsubroutine input_persistent_hydro
!***********************************************************************
    logical function output_persistent_hydro()
!
      output_persistent_hydro = .false.
!
    endfunction output_persistent_hydro
!***********************************************************************
    subroutine rprint_hydro(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
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
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_u2m=0; idiag_um2=0; idiag_oum=0; idiag_o2m=0
        idiag_uxpt=0; idiag_uypt=0; idiag_uzpt=0; idiag_dtu=0
        idiag_urms=0; idiag_umax=0; idiag_uzrms=0; idiag_uzmax=0
        idiag_orms=0; idiag_omax=0; idiag_oumphi=0
        idiag_ruxm=0; idiag_ruym=0; idiag_ruzm=0; idiag_rumax=0
        idiag_ux2m=0; idiag_uy2m=0; idiag_uz2m=0
        idiag_uxuym=0; idiag_uxuzm=0; idiag_uyuzm=0
        idiag_umx=0; idiag_umy=0; idiag_umz=0
        idiag_Marms=0; idiag_Mamax=0; idiag_divu2m=0; idiag_epsK=0
        idiag_urmphi=0; idiag_upmphi=0; idiag_uzmphi=0; idiag_u2mphi=0
        idiag_ekin=0; idiag_ekintot=0
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
      enddo
!
!  write column where which hydro variable is stored
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
    subroutine calc_mflow
!
!  dummy routine
!
!  19-jul-03/axel: adapted from hydro
!
    endsubroutine calc_mflow
!***********************************************************************
    subroutine remove_mean_momenta(f)
!
!  dummy routine
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
    subroutine hydro_clean_up
!
!  Deallocate the variables allocated in nohydro
!
!  8-sep-2009/dhruba: coded
!
      print*, 'Deallocating some nohydro variables ...'
      if (kinflow=='KS') then
         deallocate(KS_k)
         deallocate(KS_A)
         deallocate(KS_B)
         deallocate(KS_omega)
       endif
      print*, 'Done.'
!
    endsubroutine hydro_clean_up
!***********************************************************************
    subroutine kinematic_random_phase
!
!  This is a dummy routine.
!
!  16-feb-2010/bing:
!
      call fatal_error('kinematic_random_phase', &
          'Use HYDRO=hydro_kinematic in Makefile.local instead')
!
    endsubroutine kinematic_random_phase
!***********************************************************************
    subroutine kinematic_random_ampl
!
!  This is a dummy routine.
!
!  26-jun-2019/axel: coded
!
      call fatal_error('kinematic_random_ampl', &
          'Use HYDRO=hydro_kinematic in Makefile.local instead')
!
    endsubroutine kinematic_random_ampl
!***********************************************************************
    subroutine kinematic_random_wavenumber
!
!  This is a dummy routine.
!
!  26-jun-2019/axel: coded
!
      call fatal_error('kinematic_random_wavenumber', &
          'Use HYDRO=hydro_kinematic in Makefile.local instead')
!
    endsubroutine kinematic_random_wavenumber
!***********************************************************************
    subroutine find_umax(f,umax)
!
!  Dummy subroutine
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, intent(out) :: umax
!
      call keep_compiler_quiet(f)
      umax = 0.
!
    endsubroutine find_umax
!***********************************************************************
    subroutine expand_shands_hydro
!
!  Dummy
!
    endsubroutine expand_shands_hydro
!***********************************************************************
    subroutine hydro_after_timestep(f,df,dtsub)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dtsub
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dtsub)

    endsubroutine hydro_after_timestep
!***********************************************************************
    subroutine update_char_vel_hydro(f)
!
!  Dummy
!
!   25-sep-15/MR+joern: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f

      call keep_compiler_quiet(f)

    endsubroutine update_char_vel_hydro
!***********************************************************************
    subroutine calc_gradu(f)
!
! Dummy 
!
    real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)

    endsubroutine calc_gradu
!***********************************************************************
endmodule Hydro

