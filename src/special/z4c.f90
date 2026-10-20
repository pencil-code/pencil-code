! $Id$
!
!  Vacuum Einstein equations in the Z4c formulation (Bernuzzi & Hilditch
!  2010, PRD 81, 084003; Hilditch et al. 2013, PRD 88, 084057), written in
!  the form used by GR-Athena++ (Daszuta et al. 2021, ApJS 257, 25, eqs. 8-13)
!  and AthenaK (Zhu et al. 2025, ApJS 278, 52).
!
!  Evolved variables: chi, gt_ij, Khat, At_ij, Theta, Gam^i, alpha, beta^i
!  (22 PDE variables). Gauge: Bona-Masso lapse (1+log, harmonic or frozen)
!  and integrated Gamma-driver shift (or zero shift).
!
!  Initial data: flat space, linearized gravitational wave and gauge wave in
!  an arbitrary periodic direction (apples-with-apples tests, Alcubierre et
!  al. 2004, CQG 21, 589), random noise for the robust stability test, and
!  a single Brill-Lindquist puncture with pre-collapsed lapse.
!  A sinusoidal initial shift (ampl_beta) can be added to test the shift terms.
!  For the wave tests the error with respect to the exact solution is
!  available as diagnostics (gerrij, alperr).
!
!  Vacuum only, Cartesian coordinates, CPU only.
!
!  13-oct-26: coded
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 22
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED gr_sqrtg; gr_alp; gr_dalp(3); gr_beta(3); gr_dbeta(3,3)
! PENCILS PROVIDED gr_gamu(6); gr_dgam(6,3); gr_Kij(6)
!
!***************************************************************
!
module Special
!
  use Cdata
  use Quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
!  Names of the symmetric index pairs in the storage order of ij_table
!  (i11,i22,i33,i12,i23,i13 in cparam).
!
!
! Input parameters
!
  character(len=labellen) :: initz4c='flat'
  real :: ampl_wave=1e-8, ampl_noise=0.
  real, dimension(3) :: ampl_beta=0., punc_pos=0.
  real :: punc_mass=1.
  integer, dimension(3) :: nkwave=(/1,0,0/)
!
  namelist /special_init_pars/ initz4c, ampl_wave, ampl_noise, ampl_beta, nkwave, &
      punc_mass, punc_pos
!
! Run parameters
!
  character(len=labellen) :: lapse_type='1+log', shift_type='gamma-driver'
  integer :: enum_lapse_type=0
  real :: kappa1=0.02, kappa2=0., shift_mu=1., shift_eta=2., hyper3_mesh=0.
  real :: chi_floor=1e-6, r_excise_diag=0., harmonic_f=1.
  logical :: lenforce_algebraic=.true., lupw_beta=.false., lsommerfeld=.false.
!
  namelist /special_run_pars/ lapse_type, shift_type, kappa1, kappa2, &
      shift_mu, shift_eta, hyper3_mesh, chi_floor, lenforce_algebraic, lupw_beta, &
      r_excise_diag, lsommerfeld, harmonic_f

  logical, target :: lgr_metric=.true.
!
! f-array indices
!
  integer :: ichi=0, iKh=0, ith=0, ialp=0
  integer, dimension(6) :: igt=0, iAt=0
  integer, parameter :: nz4c=22   ! number of Z4c variables, ichi:ibeta(3)
  integer, dimension(3) :: iGam=0, ibeta=0
!
!
!
  integer :: iTij=0
  logical, pointer :: lgr_backreaction
!
  logical :: lshift=.true.
  real, dimension(3) :: kvec=0., nhat=0., ppol=0., qpol=0.
  real :: kmag=0.
!
!  Geometric quantities shared by dspecial_dt and calc_diagnostics_special,
!  computed in calc_pencils_special.
!
  type InternalPencils
    real, dimension(nx) :: chi, Kh, th, alp, Ktr, AA, trR, ham, divb
    real, dimension(nx,6) :: gt, At, gtu, Rij, DDalp, Atu
    real, dimension(nx,3) :: Gam, beta, dchi, dKh, dth, dalp, Gamd
    real, dimension(nx,3,3) :: dbeta, dGam, Aud, Suu
    real, dimension(nx,6,3) :: dgt, dAt
    real, dimension(nx,6,3,3) :: d2gt
    real, dimension(nx,3,3,3) :: d2beta
    real, dimension(nx,3,6) :: Gam1, Gam2
    real, dimension(nx) :: E, Strace
    real, dimension(nx,3) :: Sv
    real, dimension(nx,6) :: Sij
  endtype InternalPencils
!
  type(InternalPencils) :: q
  !$omp threadprivate(q)
!
! Diagnostic variables (needs to be consistent with reset list below)
!
  integer :: idiag_hamrms=0    ! DIAG_DOC: $\left<{\cal H}^2\right>^{1/2}$
  integer :: idiag_hammax=0    ! DIAG_DOC: $\max|{\cal H}|$
  integer :: idiag_hamxrms=0   ! DIAG_DOC: $\left<{\cal H}^2\right>^{1/2}$ for $r>r_{\rm excise}$
  integer :: idiag_momrms=0    ! DIAG_DOC: $\left<\gamma^{ij}{\cal M}_i{\cal M}_j\right>^{1/2}$
  integer :: idiag_thrms=0     ! DIAG_DOC: $\left<\Theta^2\right>^{1/2}$
  integer :: idiag_Gconrms=0   ! DIAG_DOC: $\left<(\tilde\Gamma^i-\tilde\Gamma^i_d)^2\right>^{1/2}$
  integer :: idiag_alpm=0      ! DIAG_DOC: $\left<\alpha\right>$
  integer :: idiag_alpmin=0    ! DIAG_DOC: $\min\alpha$
  integer :: idiag_alpmax=0    ! DIAG_DOC: $\max\alpha$
  integer :: idiag_chimin=0    ! DIAG_DOC: $\min\chi$
  integer :: idiag_Khrms=0     ! DIAG_DOC: $\left<\hat K^2\right>^{1/2}$
  integer :: idiag_betarms=0   ! DIAG_DOC: $\left<\beta^2\right>^{1/2}$
  integer :: idiag_alperr=0    ! DIAG_DOC: $\left<|\alpha-\alpha_{\rm exact}|\right>$
  integer :: idiag_gerrmax=0   ! DIAG_DOC: $\max|\gamma_{ij}-\gamma_{ij}^{\rm exact}|$
  integer, dimension(6) :: idiag_gerr=0 ! DIAG_DOC: $\left<|\gamma_{ij}-\gamma_{ij}^{\rm exact}|\right>$
  integer, dimension(6) :: idiag_gmx=0  ! YZAVG_DOC: $\left<\gamma_{ij}\right>_{yz}$
  integer :: idiag_alpmx=0     ! YZAVG_DOC: $\left<\alpha\right>_{yz}$
  integer :: idiag_hammx=0     ! YZAVG_DOC: $\left<{\cal H}\right>_{yz}$
!
  contains
!***********************************************************************
    subroutine register_special
!
!  Register the Z4c variables. The registration order makes all of them
!  consecutive in the f array.
!
!  13-oct-26/Claude: coded
!
      use FArrayManager
!
      integer :: i
!
      call farray_register_pde('chi',ichi)
      do i=1,6
        call farray_register_pde('gt'//csym(i),igt(i))
      enddo
      call farray_register_pde('Kh',iKh)
      do i=1,6
        call farray_register_pde('At'//csym(i),iAt(i))
      enddo
      call farray_register_pde('theta',ith)
      do i=1,3
        call farray_register_pde('Gam'//cvec(i),iGam(i))
      enddo
      call farray_register_pde('alp',ialp)
      do i=1,3
        call farray_register_pde('beta'//cvec(i),ibeta(i))
      enddo
!
      call svn_id( &
           "$Id$")

!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  13-oct-26/Claude: coded
!
      use FArrayManager, only: farray_index_by_name
      use Sub, only: cross
      use SharedVariables, only: get_shared_variable
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      integer :: i, imin, nfound
      real, dimension(3) :: e
!
      if (.not.lcartesian_coords) call fatal_error('initialize_special', &
          'z4c only works in Cartesian coordinates')
!
      select case (lapse_type)
      case ('1+log','harmonic','frozen')
      case default
        call fatal_error('initialize_special','no such lapse_type: '//trim(lapse_type))
      endselect
!
      select case (shift_type)
      case ('gamma-driver'); lshift=.true.
      case ('zero'); lshift=.false.
      case default
        call fatal_error('initialize_special','no such shift_type: '//trim(shift_type))
      endselect
!
      if (hyper3_mesh/=0. .and. nghost<3) call fatal_error('initialize_special', &
          'Hyper diffusion needs nghost>=3')
!
!  Matter sources, if another module provides them.
!
      iTij=farray_index_by_name('Tij')
      call get_shared_variable('lgr_backreaction',lgr_backreaction,default_val=.false.)
!
!  Wave vector and polarization of the wave initial data:
!  plus polarization e_ij = p_i p_j - q_i q_j with p, q orthogonal to k.
!
      kvec=2*pi*nkwave/Lxyz
      kmag=sqrt(sum(kvec**2))
      if (kmag>0.) then
        nhat=kvec/kmag
        imin=3
        do i=2,1,-1
          if (abs(nhat(i))<abs(nhat(imin))) imin=i
        enddo
        e=0.; e(imin)=1.
        call cross(nhat,e,ppol); ppol=ppol/sqrt(sum(ppol**2))
        call cross(nhat,ppol,qpol)
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  Set the ADM initial data and convert them to Z4c variables.
!
!  13-oct-26/Claude: coded
!
      use General, only: random_number_wrapper
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      real, dimension(mx,6) :: g, Kij, gt, At
      real, dimension(mx,3) :: Gam
      real, dimension(mx) :: alp, chi, Kh
      real, dimension(mx,my,mz) :: noise
      integer :: mm, nn, i
!
      do nn=1,mz; do mm=1,my
        call exact_adm(x,y(mm),z(nn),real(0.),g,Kij,alp,Gam)
        call adm_to_z4c(g,Kij,chi,gt,Kh,At)
        f(:,mm,nn,ichi)=chi
        f(:,mm,nn,iKh)=Kh
        f(:,mm,nn,ith)=0.
        f(:,mm,nn,ialp)=alp
        do i=1,6
          f(:,mm,nn,igt(i))=gt(:,i)
          f(:,mm,nn,iAt(i))=At(:,i)
        enddo
        do i=1,3
          f(:,mm,nn,iGam(i))=Gam(:,i)
          f(:,mm,nn,ibeta(i))=ampl_beta(i)*sin(kvec(1)*x+kvec(2)*y(mm)+kvec(3)*z(nn))
        enddo
      enddo; enddo
!
!  Robust stability test: random noise in all variables.
!
      if (ampl_noise/=0.) then
        do i=ichi,ibeta(3)
          call random_number_wrapper(noise)
          f(:,:,:,i)=f(:,:,:,i)+ampl_noise*(2*noise-1)
        enddo
      endif
!
    endsubroutine init_special
!***********************************************************************
    subroutine exact_adm(xx,yy,zz,tt,g,Kij,alp,Gam)
!
!  ADM variables of the initial data / exact solutions along a line
!  of points xx at fixed (yy,zz); Gam is the exact conformal connection
!  function -d_j gt^ij.
!
!  13-oct-26/Claude: coded
!
      real, dimension(:) :: xx
      real :: yy, zz, tt
      real, dimension(:,:) :: g, Kij, Gam
      real, dimension(:) :: alp
!
      real, dimension(size(xx)) :: ph, hh, dhdph, cosph
      integer :: i, j
!
      g=0.; Kij=0.; Gam=0.; alp=1.
      g(:,i11)=1.; g(:,i22)=1.; g(:,i33)=1.
!
      if (initz4c=='flat') return
      if (initz4c/='puncture') then
        ph=kvec(1)*xx+kvec(2)*yy+kvec(3)*zz-kmag*tt
        cosph=cos(ph)
      endif
!
      select case (initz4c)
!
!  Schwarzschild puncture, gamma_ij = psi^4 delta_ij with psi = 1 + M/(2r),
!  K_ij = 0 and pre-collapsed lapse alpha = psi^-2 (not a solution in time).
!
      case ('puncture')
        hh=1+0.5*punc_mass/sqrt((xx-punc_pos(1))**2+(yy-punc_pos(2))**2+(zz-punc_pos(3))**2)
        alp=hh**(-2)
        do i=1,3
          g(:,ij_table(i,i))=hh**4
        enddo
!
!  Linearized wave, ds^2 = -dt^2 + (delta_ij + H e_ij) dx^i dx^j with
!  H = A sin(k.x - |k| t). Here det(gamma) = 1 - H^2.
!
      case ('linear_wave')
        hh=ampl_wave*sin(ph)
        dhdph=ampl_wave*cosph
        do j=1,3; do i=1,j
          g(:,ij_table(i,j))=g(:,ij_table(i,j))+hh*(ppol(i)*ppol(j)-qpol(i)*qpol(j))
          Kij(:,ij_table(i,j))=0.5*kmag*dhdph*(ppol(i)*ppol(j)-qpol(i)*qpol(j))
        enddo; enddo
        do i=1,3
          Gam(:,i)=(2./3.)*kmag*nhat(i)*(1-hh**2)**(-2./3.)*hh*dhdph
        enddo
!
!  Gauge wave, ds^2 = H(-dt^2 + dn^2) + (transverse) with
!  H = 1 - A sin(k.x - |k| t). Here det(gamma) = H.
!
      case ('gauge_wave')
        hh=1-ampl_wave*sin(ph)
        dhdph=-ampl_wave*cosph
        alp=sqrt(hh)
        do j=1,3; do i=1,j
          g(:,ij_table(i,j))=g(:,ij_table(i,j))+(hh-1)*nhat(i)*nhat(j)
          Kij(:,ij_table(i,j))=0.5*kmag*dhdph/alp*nhat(i)*nhat(j)
        enddo; enddo
        do i=1,3
          Gam(:,i)=(2./3.)*kmag*nhat(i)*hh**(-5./3.)*dhdph
        enddo
!
      case default
        call fatal_error('exact_adm','no such initz4c: '//trim(initz4c))
      endselect
!
    endsubroutine exact_adm
!***********************************************************************
    subroutine adm_to_z4c(g,Kij,chi,gt,Kh,At)
!
!  13-oct-26/Claude: coded
!
      use Sub, only: sym_inverse, contract
!
      real, dimension(:,:) :: g, Kij, gt, At
      real, dimension(:) :: chi, Kh
!
      real, dimension(size(chi),6) :: gu
      real, dimension(size(chi)) :: detg
      integer :: i
!
      call sym_inverse(g,gu,detg)
      chi=detg**(-1./3.)
      Kh=contract(gu,Kij)
      do i=1,6
        gt(:,i)=chi*g(:,i)
        At(:,i)=chi*(Kij(:,i)-g(:,i)*Kh/3.)
      enddo
!
    endsubroutine adm_to_z4c
!***********************************************************************
    subroutine pencil_criteria_special
!
!  Metric pencils from special/z4c for the energy equation in curved spacetime.
!
      if (ldensity) then
        lpenc_requested(i_gr_sqrtg)=.true.
        lpenc_requested(i_gr_alp)=.true.
        lpenc_requested(i_gr_dalp)=.true.
        lpenc_requested(i_gr_beta)=.true.
        lpenc_requested(i_gr_gamu)=.true.
        lpenc_requested(i_gr_Kij)=.true.
      endif
!
!  Metric pencils from special/z4c for the momentum equation in curved spacetime.
!
      if (lhydro) then
        lpenc_requested(i_gr_sqrtg)=.true.
        lpenc_requested(i_gr_alp)=.true.
        lpenc_requested(i_gr_dalp)=.true.
        lpenc_requested(i_gr_beta)=.true.
        lpenc_requested(i_gr_dbeta)=.true.
        lpenc_requested(i_gr_gamu)=.true.
        lpenc_requested(i_gr_dgam)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Enforce the algebraic constraints det(gt)=1 and gt^ij At_ij=0 at every
!  substep, and keep chi above chi_floor.
!
!  13-oct-26/Claude
!
      use Sub, only: sym_inverse, contract
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      real, dimension(nx,6) :: gt, gtu
      real, dimension(nx) :: detg, trA
      integer :: i, mm, nn
!
      if (.not.lenforce_algebraic) return
!
      do nn=n1,n2; do mm=m1,m2
        do i=1,6
          gt(:,i)=f(l1:l2,mm,nn,igt(i))
        enddo
        call sym_inverse(gt,gtu,detg)
        detg=detg**(-1./3.)
        do i=1,6
          gt(:,i)=gt(:,i)*detg
          f(l1:l2,mm,nn,igt(i))=gt(:,i)
          gtu(:,i)=gtu(:,i)/detg
        enddo
        trA=contract(gtu,f(l1:l2,mm,nn,iAt))
        do i=1,6
          f(l1:l2,mm,nn,iAt(i))=f(l1:l2,mm,nn,iAt(i))-gt(:,i)*trA/3.
        enddo
        f(l1:l2,mm,nn,ichi)=max(f(l1:l2,mm,nn,ichi),chi_floor)
      enddo; enddo
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Geometric quantities shared by the right-hand side and the diagnostics.
!
!  13-oct-26/Claude: coded
!
      use Sub, only: grad, g2ij, sym_inverse, contract,covariant_g2ij
      use Sub, only: lower_first_index,raise_second_index
!
      real, contiguous, dimension(:,:,:,:) :: f
      type(pencil_case) :: p
!
      real, dimension(nx,6) :: DDchi
      real, dimension(nx) :: detg, lapchi, dchi2, dchida, rhs
      real, dimension(nx,3,3) :: Smix
      real, dimension(nx,6) :: gam
      integer :: i, j, k, l, mm, ij
!
!  Variables and their derivatives.
!
      q%chi=f(l1:l2,m,n,ichi); q%Kh=f(l1:l2,m,n,iKh); q%th=f(l1:l2,m,n,ith)
      q%alp=f(l1:l2,m,n,ialp)
      do i=1,6
        q%gt(:,i)=f(l1:l2,m,n,igt(i)); q%At(:,i)=f(l1:l2,m,n,iAt(i))
      enddo
      do i=1,3
        q%Gam(:,i)=f(l1:l2,m,n,iGam(i)); q%beta(:,i)=f(l1:l2,m,n,ibeta(i))
      enddo
!
      call grad(f,ichi,q%dchi); 
      call grad(f,ialp,q%dalp); 
      call grad(f,iKh,q%dKh)
      call grad(f,ith,q%dth)
      call g2ij(f,igt(1),q%d2gt)
      call grad(f,igt(1),q%dgt)
      call grad(f,iAt(1),q%dAt)
      call grad(f,iGam(1),q%dGam)
      if (lshift) then
        call g2ij(f,ibeta(1),q%d2beta)
        call grad(f,ibeta(1),q%dbeta)
        q%divb=q%dbeta(:,1,1)+q%dbeta(:,2,2)+q%dbeta(:,3,3)
      endif
!
!  Conformal Christoffel symbols: Gam1(k,ij)=Gamt_kij, Gam2(k,ij)=Gamt^k_ij,
!  and the contracted ones Gamd^k computed from the metric.
!
      do k=1,3; do j=1,3; do i=1,j
        q%Gam1(:,k,ij_table(i,j))=0.5*(q%dgt(:,ij_table(j,k),i)+q%dgt(:,ij_table(i,k),j)-q%dgt(:,ij_table(i,j),k))
      enddo; enddo; enddo

      call sym_inverse(q%gt,q%gtu,detg)

      do k=1,3; do ij=1,6
        q%Gam2(:,k,ij)=q%gtu(:,ij_table(k,1))*q%Gam1(:,1,ij)+q%gtu(:,ij_table(k,2))*q%Gam1(:,2,ij) &
                      +q%gtu(:,ij_table(k,3))*q%Gam1(:,3,ij)
      enddo; enddo
      do k=1,3
        q%Gamd(:,k)=contract(q%gtu,q%Gam2(:,k,:))
      enddo

      call covariant_g2ij(f,ichi,q%dchi,q%Gam2,DDchi)
      lapchi=contract(q%gtu,DDchi)
      dchi2=0.; dchida=0.
      do j=1,3; do i=1,3
        dchi2=dchi2+q%gtu(:,ij_table(i,j))*q%dchi(:,i)*q%dchi(:,j)
        dchida=dchida+q%gtu(:,ij_table(i,j))*q%dchi(:,i)*q%dalp(:,j)
      enddo; enddo
!
!  Ricci tensor R_ij = Rt_ij + Rchi_ij 
!
      do j=1,3; do i=1,j

        !Rt_ij
        ij=ij_table(i,j)
        rhs=0.
        do l=1,3; do mm=1,3
          rhs=rhs-0.5*q%gtu(:,ij_table(l,mm))*q%d2gt(:,ij,l,mm)
        enddo; enddo
        do k=1,3
          rhs=rhs+0.5*(q%gt(:,ij_table(k,i))*q%dGam(:,k,j)+q%gt(:,ij_table(k,j))*q%dGam(:,k,i)) &
                 +0.5*q%Gamd(:,k)*(q%Gam1(:,i,ij_table(j,k))+q%Gam1(:,j,ij_table(i,k)))
          do l=1,3; do mm=1,3
            rhs=rhs+q%gtu(:,ij_table(l,mm))*(q%Gam2(:,k,ij_table(l,i))*q%Gam1(:,j,ij_table(k,mm)) &
                                        +q%Gam2(:,k,ij_table(l,j))*q%Gam1(:,i,ij_table(k,mm)) &
                                        +q%Gam2(:,k,ij_table(i,mm))*q%Gam1(:,k,ij_table(l,j)))
          enddo; enddo
        enddo

        !Rchi_ij
        rhs=rhs+0.5/q%chi*(DDchi(:,ij)+q%gt(:,ij)*lapchi)&
            -0.25/q%chi**2*q%dchi(:,i)*q%dchi(:,j)-0.75/q%chi**2*q%gt(:,ij)*dchi2
        q%Rij(:,ij)=rhs
      enddo; enddo

!
!  Physical second covariant derivative of the lapse.
!
      call covariant_g2ij(f,ialp,q%dalp,q%Gam2,q%DDalp)
      do j=1,3; do i=1,j
        ij=ij_table(i,j)
        q%DDalp(:,ij)=q%DDalp(:,ij)+0.5/q%chi*(q%dalp(:,i)*q%dchi(:,j)+q%dalp(:,j)*q%dchi(:,i) &
                                               -q%gt(:,ij)*dchida)
      enddo; enddo

      q%trR=q%chi*contract(q%gtu,q%Rij)
!
!  Raised extrinsic curvature: Aud(i,j)=At^i_j, Atu(ij)=At^ij.
!
      do j=1,3; do i=1,3
        q%Aud(:,i,j)=q%gtu(:,ij_table(i,1))*q%At(:,ij_table(1,j))+q%gtu(:,ij_table(i,2))*q%At(:,ij_table(2,j)) &
                    +q%gtu(:,ij_table(i,3))*q%At(:,ij_table(3,j))
      enddo; enddo
      do j=1,3; do i=1,j
        q%Atu(:,ij_table(i,j))=q%Aud(:,i,1)*q%gtu(:,ij_table(1,j))+q%Aud(:,i,2)*q%gtu(:,ij_table(2,j)) &
                          +q%Aud(:,i,3)*q%gtu(:,ij_table(3,j))
      enddo; enddo
      q%AA=0.
      do j=1,3; do i=1,3
        q%AA=q%AA+q%Aud(:,i,j)*q%Aud(:,j,i)
      enddo; enddo
      q%Ktr=q%Kh+2*q%th
      q%ham=q%trR-q%AA+(2./3.)*q%Ktr**2
!
!  Metric pencils for the matter equations (hydro and density with lgr_metric):
!  sqrt(gamma), alpha, d_i alpha, beta^i, d_j beta^i (gr_dbeta(:,i,j)),
!  gamma^ij, d_k gamma_ij (gr_dgam(:,ij,k)) and K_ij, symmetric pairs in ij_table order.
!  Without shift evolution, beta is zero.
!
      if (lpencil(i_gr_sqrtg)) p%gr_sqrtg=sqrt(detg/q%chi**3)
      if (lpencil(i_gr_alp)) p%gr_alp=q%alp
      if (lpencil(i_gr_dalp)) p%gr_dalp=q%dalp
      if (lpencil(i_gr_beta)) then
        if (lshift) then
          p%gr_beta=q%beta
        else
          p%gr_beta=0.
        endif
      endif
      if (lpencil(i_gr_dbeta)) then
        if (lshift) then
          p%gr_dbeta=q%dbeta
        else
          p%gr_dbeta=0.
        endif
      endif
      if (lpencil(i_gr_gamu)) p%gr_gamu=q%gtu*spread(q%chi,2,6)
      if (lpencil(i_gr_dgam)) then
        do k=1,3; do ij=1,6
          p%gr_dgam(:,ij,k)=(q%dgt(:,ij,k)-q%gt(:,ij)*q%dchi(:,k)/q%chi)/q%chi
        enddo; enddo
      endif
      if (lpencil(i_gr_Kij)) then
        do ij=1,6
          p%gr_Kij(:,ij)=(q%At(:,ij)+q%gt(:,ij)*q%Ktr/3.)/q%chi
        enddo
      endif
!
!  Matter sources; S = gamma^ij S_ij = chi gt^ij S_ij.
!
      if (lgr_backreaction) then
        do i=1,3; do j=1,3
          Smix(:,j,i)=(f(l1:l2,m,n,iTij+3*(i-1)+j-1)+p%gr_beta(:,j)*f(l1:l2,m,n,iux+i-1))/p%gr_alp
        enddo; enddo
        do i=1,6
          gam(:,i)=q%gt(:,i)/q%chi
        enddo

        call raise_second_index(Smix,p%gr_gamu,q%Suu)
        call lower_first_index(Smix,gam,q%Sij)
        do i=1,6
          q%Sij(:,i) = q%Sij(:,i)/p%gr_sqrtg 
        enddo

        q%E=f(l1:l2,m,n,irho)/p%gr_sqrtg
        do i=1,3
          q%Sv(:,i)=f(l1:l2,m,n,iux-1+i)/p%gr_sqrtg
        enddo
        q%Strace=contract(p%gr_gamu,q%Sij)
        q%ham=q%ham-16*pi*q%E
      endif
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  Right-hand side of the Z4c equations, GR-Athena++ eqs. (8)-(13) and (22)
!  without matter.
!
!  13-oct-26/Claude: coded
!
      use Deriv, only: der6
      use Sub, only: u_dot_grad, contract
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(mx,my,mz,mvar) :: df
      type(pencil_case) :: p
!
      real, dimension(nx,6) :: Xij
      real, dimension(nx,nz4c) :: df_in
      real, dimension(nx) :: lapalp, trX, rhs, tmp, speed
      logical, dimension(nx) :: lbound
      integer :: i, j, k, ij
!
!  Points on non-periodic outer boundaries get the Sommerfeld condition
!  instead of the evolution equations, see sommerfeld below.
!
      lbound=.false.
      if (lsommerfeld) then
        lbound=sommerfeld_mask()
        if (any(lbound)) call save_df_z4c(df,df_in)
      endif
!
      lapalp=q%chi*contract(q%gtu,q%DDalp)
!
!  chi
!
      rhs=(2./3.)*q%chi*q%alp*q%Ktr
      if (lshift) then
        rhs=rhs-(2./3.)*q%chi*q%divb
        call u_dot_grad(f,ichi,q%dchi,q%beta,rhs,upwind=lupw_beta,ladd=.true.)
      endif
      df(l1:l2,m,n,ichi)=df(l1:l2,m,n,ichi)+rhs
!
!  gt_ij and At_ij
!
      do j=1,3; do i=1,j
        ij=ij_table(i,j)
        Xij(:,ij)=-q%DDalp(:,ij)+q%alp*q%Rij(:,ij)
        if (lgr_backreaction) Xij(:,ij)=Xij(:,ij)-8*pi*q%alp*q%Sij(:,ij)
      enddo; enddo
      trX=contract(q%gtu,Xij)
      do j=1,3; do i=1,j
        ij=ij_table(i,j)
        rhs=-2*q%alp*q%At(:,ij)
        if (lshift) then
          rhs=rhs-(2./3.)*q%gt(:,ij)*q%divb
          call u_dot_grad(f,igt(ij),q%dgt(:,ij,:),q%beta,rhs,upwind=lupw_beta,ladd=.true.)
          do k=1,3
            rhs=rhs+q%gt(:,ij_table(i,k))*q%dbeta(:,k,j)+q%gt(:,ij_table(j,k))*q%dbeta(:,k,i)
          enddo
        endif
        df(l1:l2,m,n,igt(ij))=df(l1:l2,m,n,igt(ij))+rhs
!
        rhs=q%chi*(Xij(:,ij)-q%gt(:,ij)*trX/3.)+q%alp*q%Ktr*q%At(:,ij)
        do k=1,3
          rhs=rhs-2*q%alp*q%At(:,ij_table(i,k))*q%Aud(:,k,j)
        enddo
        if (lshift) then
          rhs=rhs-(2./3.)*q%At(:,ij)*q%divb
          call u_dot_grad(f,iAt(ij),q%dAt(:,ij,:),q%beta,rhs,upwind=lupw_beta,ladd=.true.)
          do k=1,3
            rhs=rhs+q%At(:,ij_table(i,k))*q%dbeta(:,k,j)+q%At(:,ij_table(j,k))*q%dbeta(:,k,i)
          enddo
        endif
        df(l1:l2,m,n,iAt(ij))=df(l1:l2,m,n,iAt(ij))+rhs
      enddo; enddo
!
!  Khat and Theta
!
      rhs=-lapalp+q%alp*(q%AA+q%Ktr**2/3.)+q%alp*kappa1*(1-kappa2)*q%th
      if (lgr_backreaction) rhs=rhs+4*pi*q%alp*(q%Strace+q%E)
      if (lshift) call u_dot_grad(f,iKh,q%dKh,q%beta,rhs,upwind=lupw_beta,ladd=.true.)
      df(l1:l2,m,n,iKh)=df(l1:l2,m,n,iKh)+rhs
!
      rhs=0.5*q%alp*(q%ham-2*kappa1*(2+kappa2)*q%th)
      if (lshift) call u_dot_grad(f,ith,q%dth,q%beta,rhs,upwind=lupw_beta,ladd=.true.)
      df(l1:l2,m,n,ith)=df(l1:l2,m,n,ith)+rhs
!
!  Gam^i
!
      do i=1,3
        rhs=2*q%alp*(contract(q%Atu,q%Gam2(:,i,:))-kappa1*(q%Gam(:,i)-q%Gamd(:,i)))
        do j=1,3
          rhs=rhs-2*q%Atu(:,ij_table(i,j))*(q%dalp(:,j)+1.5*q%alp/q%chi*q%dchi(:,j)) &
                 -(2./3.)*q%alp*q%gtu(:,ij_table(i,j))*(2*q%dKh(:,j)+q%dth(:,j))
          if (lgr_backreaction) rhs=rhs-16*pi*q%alp*q%gtu(:,ij_table(i,j))*q%Sv(:,j)
        enddo
        if (lshift) then
          rhs=rhs+(2./3.)*q%Gamd(:,i)*q%divb
          call u_dot_grad(f,iGam(i),q%dGam(:,i,:),q%beta,rhs,upwind=lupw_beta,ladd=.true.)
          do j=1,3
            do k=1,3
              rhs=rhs+q%gtu(:,ij_table(j,k))*q%d2beta(:,i,j,k)+(1./3.)*q%gtu(:,ij_table(i,j))*q%d2beta(:,k,k,j)
            enddo
            rhs=rhs-q%Gamd(:,j)*q%dbeta(:,i,j)
          enddo
        endif
        df(l1:l2,m,n,iGam(i))=df(l1:l2,m,n,iGam(i))+rhs
      enddo
!
!  Gauge: Bona-Masso lapse and integrated Gamma-driver shift.
!
      select case (lapse_type)
      case ('1+log');    rhs=-2*q%alp*q%Kh
      case ('harmonic'); rhs=-harmonic_f*q%alp**2*q%Kh
      case ('frozen');   rhs=0.
      endselect
      if (lshift) call u_dot_grad(f,ialp,q%dalp,q%beta,rhs,upwind=lupw_beta,ladd=.true.)
      df(l1:l2,m,n,ialp)=df(l1:l2,m,n,ialp)+rhs
!
      if (lshift) then
        do i=1,3
          rhs=shift_mu*q%Gam(:,i)-shift_eta*q%beta(:,i)
          call u_dot_grad(f,ibeta(i),q%dbeta(:,i,:),q%beta,rhs,upwind=lupw_beta,ladd=.true.)
          df(l1:l2,m,n,ibeta(i))=df(l1:l2,m,n,ibeta(i))+rhs
        enddo
      endif
!
!  sixth-order hyperdiffusion, sigma/64 dx^5 d^6/dx^6.
!
!  The variables are visited by name (not by a loop over ichi:ibeta(3)) so
!  that the GPU transpiler can map each df access to its field.
!
      if (hyper3_mesh/=0.) then
        call hyper_diffus(f,df,ichi)
        call hyper_diffus(f,df,iKh)
        call hyper_diffus(f,df,ith)
        call hyper_diffus(f,df,ialp)
        do ij=1,6
          call hyper_diffus(f,df,igt(ij))
          call hyper_diffus(f,df,iAt(ij))
        enddo
        do i=1,3
          call hyper_diffus(f,df,iGam(i))
          if (lshift) call hyper_diffus(f,df,ibeta(i))
        enddo
      endif
!
      if (any(lbound)) call sommerfeld(f,df,df_in,lbound)
!
!  Time step: light speed and gauge speed along each direction.
!
      if (lupdate_courant_dt) then
        select case (lapse_type)
        case ('1+log'); tmp=max(1.,sqrt(2/q%alp))
        case ('harmonic'); tmp=max(1.,sqrt(harmonic_f))
        case default;   tmp=1.
        endselect
        speed=0.
        do j=1,3
          speed=speed+(abs(q%beta(:,j))+q%alp*tmp*sqrt(q%chi*q%gtu(:,ij_table(j,j))))*dline_1(:,j)
        enddo
        maxadvec=max(maxadvec,speed)
      endif
!
      call calc_diagnostics_special(f,p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!  Metric/source part of the right-hand side of sqrt(gamma)*S_i for lgr_metric=T
!  (e.g. Banyuls et al. 1997, ApJ 476, 221; Font 2008, Living Rev. Rel. 11, 7):
!
!    d_t(sg S_i) + d_j F^j_i = sg (-E d_i alpha + S_k d_i beta^k + alpha/2 S^jk d_i gamma_jk)
!
!  The metric quantities are pencils of special/z4c; S^j_i is recovered from
!  the stored flux F^j_i. Flux is added in hydro
!
!
!  19-oct-26/TP: coded
!
      use Sub, only: div_tensor

      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p

      real, dimension(nx) :: EE, rhs, tmp
      integer :: i, j, k
!
      do i=1,3
        rhs=-f(l1:l2,m,n,irho)*p%gr_dalp(:,i)
        do k=1,3
          rhs=rhs+f(l1:l2,m,n,iuu-1+k)*p%gr_dbeta(:,k,i)
          do j=1,3
            rhs=rhs+.5*p%gr_alp*q%Suu(:,j,k)*p%gr_dgam(:,ij_table(j,k),i)
          enddo
        enddo
        df(l1:l2,m,n,iuu+i-1) = df(l1:l2,m,n,iuu+i-1)+ rhs
      enddo
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!  Metric/source part of the right-hand side of sqrt(gamma)*E (in irho) for lgr_metric=T
!  (e.g. Banyuls et al. 1997, ApJ 476, 221; Font 2008, Living Rev. Rel. 11, 7):
!
!    d_t(sg E) + d_j F_E^j = sg (alpha S^jk K_jk - S^j d_j alpha)
!
!  The metric quantities are pencils of special/z4c. The fluxes F_E^j and
!  T^j_i are computed by hydro (hydro_after_boundary_gr); S^j_i is recovered
!  from T^j_i = sg (alpha S^j_i - beta^j S_i).
!
!  14-oct-26/Claude: moved from hydro
!  19-oct-26/TP: moved to z4c and skipped unnecessary multiplications by sqrt(gamma)
!
      use Deriv, only: der
      use Sub, only: raise_index
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type (pencil_case), intent(in) :: p

      real, dimension(nx) :: rhs
!
      real, dimension(nx,3) :: Su
      real, dimension(nx) :: tmp,divFgrE
      integer :: i, j, k
!
      call raise_index(p%gr_gamu,f(l1:l2,m,n,iux:iuz),Su)
!
      rhs=0.
      do j=1,3
        rhs=rhs-Su(:,j)*p%gr_dalp(:,j)
        do k=1,3
          rhs=rhs+p%gr_alp*q%Suu(:,j,k)*p%gr_Kij(:,ij_table(j,k))
        enddo
      enddo
      df(l1:l2,m,n,irho) = df(l1:l2,m,n,irho) + rhs
    endsubroutine special_calc_density
!***********************************************************************
!
    function sommerfeld_mask()
!
!  Points of the current pencil that lie on a non-periodic outer boundary.
!
!  13-oct-26/Claude: coded
!
      logical, dimension(nx) :: sommerfeld_mask
!
      sommerfeld_mask=.false.
!  The x boundary points are selected pointwise (not as elements 1 and nx of
!  the pencil) so that the GPU transpiler can map the mask.
!
      if (.not.lperi(1)) then
        if (lfirst_proc_x) sommerfeld_mask=sommerfeld_mask .or. x(l1:l2)==x(l1)
        if (llast_proc_x) sommerfeld_mask=sommerfeld_mask .or. x(l1:l2)==x(l2)
      endif
      if (.not.lperi(2)) then
        if ((lfirst_proc_y.and.m==m1).or.(llast_proc_y.and.m==m2)) sommerfeld_mask=.true.
      endif
      if (.not.lperi(3)) then
        if ((lfirst_proc_z.and.n==n1).or.(llast_proc_z.and.n==n2)) sommerfeld_mask=.true.
      endif
!
    endfunction sommerfeld_mask
!***********************************************************************
    subroutine hyper_diffus(f,df,k)
!
!  Sixth-order hyperdiffusion, sigma/64 dx^5 d^6/dx^6, of variable k.
!  Called Kreiss-Oliger dissipation by other NR codes.
!
!  20-oct-26/TP: coded
!
      use Deriv, only: der6
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(mx,my,mz,mvar) :: df
      integer :: k
!
      real, dimension(nx) :: tmp
      integer :: j
!
      do j=1,3
        call der6(f,k,tmp,j,ignoredx=.true.)
        df(l1:l2,m,n,k)=df(l1:l2,m,n,k)+(hyper3_mesh/64.)*dline_1(:,j)*tmp
      enddo
!
    endsubroutine hyper_diffus
!***********************************************************************
    subroutine save_df_z4c(df,df_in)
!
!  Copy of df of the Z4c variables, in the order chi, gt(1:6), Kh, At(1:6),
!  theta, Gam(1:3), alp, beta(1:3).
!
!  15-oct-26/Claude: coded
!
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx,nz4c) :: df_in
!
      integer :: i, ij
!
      df_in(:,1)=df(l1:l2,m,n,ichi)
      df_in(:,8)=df(l1:l2,m,n,iKh)
      df_in(:,15)=df(l1:l2,m,n,ith)
      df_in(:,19)=df(l1:l2,m,n,ialp)
      do ij=1,6
        df_in(:,1+ij)=df(l1:l2,m,n,igt(ij))
        df_in(:,8+ij)=df(l1:l2,m,n,iAt(ij))
      enddo
      do i=1,3
        df_in(:,15+i)=df(l1:l2,m,n,iGam(i))
        df_in(:,19+i)=df(l1:l2,m,n,ibeta(i))
      enddo
!
    endsubroutine save_df_z4c
!***********************************************************************
    subroutine sommerfeld(f,df,df_in,lbound)
!
!  Radiative (Sommerfeld) condition d_t u = -(x^i/r) d_i u - (u-u0)/r with
!  unit speed and flat-space asymptotic values u0, replacing the right-hand
!  side at the boundary points lbound. df_in is df before dspecial_dt added
!  to it (see save_df_z4c for its order). u0=1 for chi, alp and the diagonal
!  of gt (ij_table puts 11, 22, 33 first).
!
!  13-oct-26/Claude: coded
!  15-oct-26/TP: variables visited by name for the GPU transpiler
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx,nz4c) :: df_in
      logical, dimension(nx) :: lbound
!
      real :: u0
      integer :: i, ij
!
      call sommerfeld_var(f,df,ichi,df_in(:,1),1.,lbound)
      call sommerfeld_var(f,df,iKh,df_in(:,8),0.,lbound)
      call sommerfeld_var(f,df,ith,df_in(:,15),0.,lbound)
      call sommerfeld_var(f,df,ialp,df_in(:,19),1.,lbound)
      do ij=1,6
        u0=0.
        if (ij<=3) u0=1.
        call sommerfeld_var(f,df,igt(ij),df_in(:,1+ij),u0,lbound)
        call sommerfeld_var(f,df,iAt(ij),df_in(:,8+ij),0.,lbound)
      enddo
      do i=1,3
        call sommerfeld_var(f,df,iGam(i),df_in(:,15+i),0.,lbound)
        call sommerfeld_var(f,df,ibeta(i),df_in(:,19+i),0.,lbound)
      enddo
!
    endsubroutine sommerfeld
!***********************************************************************
    subroutine sommerfeld_var(f,df,k,dfk_in,u0,lbound)
!
!  Sommerfeld condition for variable k, see sommerfeld.
!
!  15-oct-26/Claude: split off from sommerfeld
!
      use Sub, only: grad
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(mx,my,mz,mvar) :: df
      integer :: k
      real, dimension(nx) :: dfk_in
      real :: u0
      logical, dimension(nx) :: lbound
!
      real, dimension(nx,3) :: gu
      real, dimension(nx) :: rr
!
      rr=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
      call grad(f,k,gu)
      where (lbound) df(l1:l2,m,n,k)=dfk_in &
          -(x(l1:l2)*gu(:,1)+y(m)*gu(:,2)+z(n)*gu(:,3))/rr-(f(l1:l2,m,n,k)-u0)/rr
!
    endsubroutine sommerfeld_var
!***********************************************************************
    subroutine calc_diagnostics_special(f,p)
!
!  13-oct-26/Claude: coded
!
      use Diagnostics
!
      real, contiguous, dimension(:,:,:,:) :: f
      type(pencil_case) :: p
!
      real, dimension(nx,6) :: gex, Kex
      real, dimension(nx,3) :: mom, Gamex
      real, dimension(nx) :: tmp
      integer :: i, j, k, l
!
      if (ldiagnos) then
        call sum_mn_name(q%ham**2,idiag_hamrms,lsqrt=.true.)
        call max_mn_name(abs(q%ham),idiag_hammax)
        if (idiag_hamxrms/=0) then
          tmp=sqrt((x(l1:l2)-punc_pos(1))**2+(y(m)-punc_pos(2))**2+(z(n)-punc_pos(3))**2)
          call sum_mn_name(merge(q%ham**2,0.,tmp>r_excise_diag),idiag_hamxrms,lsqrt=.true.)
        endif
        if (idiag_momrms/=0) then
!
!  Momentum constraint M_i = Dt_j At^j_i - 3/(2 chi) At^j_i d_j chi - 2/3 d_i K.
!
          do i=1,3
            mom(:,i)=-(2./3.)*(q%dKh(:,i)+2*q%dth(:,i))
            if (lgr_backreaction) mom(:,i)=mom(:,i)-8*pi*q%Sv(:,i)
            do j=1,3
              mom(:,i)=mom(:,i)-1.5/q%chi*q%Aud(:,j,i)*q%dchi(:,j)
              do k=1,3
                tmp=q%dAt(:,ij_table(k,i),j)
                do l=1,3
                  tmp=tmp-q%Gam2(:,l,ij_table(j,k))*q%At(:,ij_table(l,i))-q%Gam2(:,l,ij_table(j,i))*q%At(:,ij_table(k,l))
                enddo
                mom(:,i)=mom(:,i)+q%gtu(:,ij_table(j,k))*tmp
              enddo
            enddo
          enddo
          tmp=0.
          do j=1,3; do i=1,3
            tmp=tmp+q%chi*q%gtu(:,ij_table(i,j))*mom(:,i)*mom(:,j)
          enddo; enddo
          call sum_mn_name(tmp,idiag_momrms,lsqrt=.true.)
        endif
        call sum_mn_name(q%th**2,idiag_thrms,lsqrt=.true.)
        call sum_mn_name(sum((q%Gam-q%Gamd)**2,2),idiag_Gconrms,lsqrt=.true.)
        call sum_mn_name(q%alp,idiag_alpm)
        call max_mn_name(-q%alp,idiag_alpmin,lneg=.true.)
        call max_mn_name(q%alp,idiag_alpmax)
        call max_mn_name(-q%chi,idiag_chimin,lneg=.true.)
        call sum_mn_name(q%Kh**2,idiag_Khrms,lsqrt=.true.)
        call sum_mn_name(sum(q%beta**2,2),idiag_betarms,lsqrt=.true.)
!
!  Error with respect to the exact solution (wave tests).
!
        if (any(idiag_gerr/=0).or.idiag_gerrmax/=0.or.idiag_alperr/=0) then
          call exact_adm(x(l1:l2),y(m),z(n),real(t),gex,Kex,tmp,Gamex)
          call sum_mn_name(abs(q%alp-tmp),idiag_alperr)
          do i=1,6
            gex(:,i)=q%gt(:,i)/q%chi-gex(:,i)
            call sum_mn_name(abs(gex(:,i)),idiag_gerr(i))
          enddo
          call max_mn_name(maxval(abs(gex),2),idiag_gerrmax)
        endif
      endif
!
      if (l1davgfirst) then
        do i=1,6
          call yzsum_mn_name_x(q%gt(:,i)/q%chi,idiag_gmx(i))
        enddo
        call yzsum_mn_name_x(q%alp,idiag_alpmx)
        call yzsum_mn_name_x(q%ham,idiag_hammx)
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_diagnostics_special
!***********************************************************************
    subroutine read_special_init_pars(iomsg)
!
      use File_io, only: parallel_unit
      use SharedVariables, only: put_shared_variable
!
      character(len=*), intent(out) :: iomsg
      integer :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""
!TP: Bit naughty to put it here but has to come before register_hydro
      call put_shared_variable('lgr_metric',lgr_metric,caller='read_special_init_pars')
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
    subroutine read_special_run_pars(iomsg)
!
      use File_io, only: parallel_unit
!
      character(len=*), intent(out) :: iomsg
      integer :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""
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
    subroutine rprint_special(lreset,lwrite)
!
!  13-oct-26/Claude: coded
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, i
!
      if (lreset) then
        idiag_hamrms=0; idiag_hammax=0; idiag_hamxrms=0; idiag_momrms=0; idiag_thrms=0
        idiag_Gconrms=0; idiag_alpm=0; idiag_alpmin=0; idiag_alpmax=0
        idiag_chimin=0; idiag_Khrms=0; idiag_betarms=0; idiag_alperr=0
        idiag_gerrmax=0; idiag_gerr=0; idiag_gmx=0; idiag_alpmx=0; idiag_hammx=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'hamrms',idiag_hamrms)
        call parse_name(iname,cname(iname),cform(iname),'hammax',idiag_hammax)
        call parse_name(iname,cname(iname),cform(iname),'hamxrms',idiag_hamxrms)
        call parse_name(iname,cname(iname),cform(iname),'momrms',idiag_momrms)
        call parse_name(iname,cname(iname),cform(iname),'thrms',idiag_thrms)
        call parse_name(iname,cname(iname),cform(iname),'Gconrms',idiag_Gconrms)
        call parse_name(iname,cname(iname),cform(iname),'alpm',idiag_alpm)
        call parse_name(iname,cname(iname),cform(iname),'alpmin',idiag_alpmin)
        call parse_name(iname,cname(iname),cform(iname),'alpmax',idiag_alpmax)
        call parse_name(iname,cname(iname),cform(iname),'chimin',idiag_chimin)
        call parse_name(iname,cname(iname),cform(iname),'Khrms',idiag_Khrms)
        call parse_name(iname,cname(iname),cform(iname),'betarms',idiag_betarms)
        call parse_name(iname,cname(iname),cform(iname),'alperr',idiag_alperr)
        call parse_name(iname,cname(iname),cform(iname),'gerrmax',idiag_gerrmax)
        do i=1,6
          call parse_name(iname,cname(iname),cform(iname),'gerr'//csym(i),idiag_gerr(i))
        enddo
      enddo
!
      do inamex=1,nnamex
        do i=1,6
          call parse_name(inamex,cnamex(inamex),cformx(inamex),'g'//csym(i)//'mx',idiag_gmx(i))
        enddo
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'alpmx',idiag_alpmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'hammx',idiag_hammx)
      enddo
!
      call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=100
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(kappa1,p_par(1))
    call copy_addr(kappa2,p_par(2))
    call copy_addr(shift_mu,p_par(3))
    call copy_addr(shift_eta,p_par(4))
    call copy_addr(hyper3_mesh,p_par(5))
    call copy_addr(harmonic_f,p_par(6))
    call copy_addr(lupw_beta,p_par(7)) ! bool
    call copy_addr(lsommerfeld,p_par(8)) ! bool
    call copy_addr(lgr_metric,p_par(9)) ! bool
    call copy_addr(ichi,p_par(10)) ! int
    call copy_addr(ikh,p_par(11)) ! int
    call copy_addr(ith,p_par(12)) ! int
    call copy_addr(ialp,p_par(13)) ! int
    call copy_addr(lshift,p_par(16)) ! bool
    call copy_addr(igt,p_par(17)) ! int (6)
    call copy_addr(iat,p_par(18)) ! int (6)
    call copy_addr(igam,p_par(19)) ! int3
    call copy_addr(ibeta,p_par(20)) ! int3
    call string_to_enum(enum_lapse_type,lapse_type)
    call copy_addr(enum_lapse_type,p_par(23)) ! int
    call copy_addr(lenforce_algebraic,p_par(24)) ! bool
    call copy_addr(chi_floor,p_par(25))
    call copy_addr(iTij,p_par(26)) ! int

  endsubroutine pushpars2c
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
!
endmodule Special
