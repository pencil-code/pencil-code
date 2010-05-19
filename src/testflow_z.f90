! $Id$
!
!  This modules deals with all aspects of testflow fields; if no
!  testflow fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testflow relevant subroutines listed in here.
!
!  Note: this routine requires that MVAR and MAUX contributions
!  together with njtestflow are set correctly in the cparam.local file.
!  njtestflow must be set at the end of the file such that
!  4*(njtestflow+1)=MVAR.
!
!  Example:
!  ! MVAR CONTRIBUTION 28
!  ! MAUX CONTRIBUTION 28
!  integer, parameter :: njtestflow=6
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestflow = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Testflow
!
  use Cparam
  use Messages
  use Cdata, only : mname
!
  implicit none
!
  include 'testflow.h'
!
  interface insert                ! Overload the `insert' function
    module procedure insert_array
    module procedure insert_array_mult
  endinterface
!
  integer, parameter :: njtestflow=6  !obsolete, but needed for compilation (MR)
!
! Slice precalculation buffers
!
  real, target, dimension (nx,ny,3) :: uu11_xy, &
                                       uu11_xy2
  real, target, dimension (nx,nz,3) :: uu11_xz
  real, target, dimension (ny,nz,3) :: uu11_yz
!
!  cosine and sine function for setting test fields and analysis
!
  real, dimension(mz) :: cz,sz,k1cz,k1sz,zq2
!
  real, parameter :: cs2test=1., cs2test1=1./cs2test
!
  logical :: ltestflow_upw_lnrho=.false., ltestflow_upw_uu=.false.
!
  integer :: nuuinit
  real :: tuuinit=0.,duuinit=0.
!
  ! input parameters
!
  character (len=labellen), dimension(ninit*njtestflow) :: inituutest='nothing'
  real, dimension (ninit*njtestflow) :: ampluutest=0.
!
  real :: kx_uutest=1.
!
  logical :: zextent=.true.,    &
             lugu_as_aux=.false.
!
  namelist /testflow_init_pars/ zextent,    &    ! ??
                                inituutest, &    ! name(s) of initialization kind (up to five different to be overlayed),
                                                 ! legal values: 'zero', 'gaussian-noise-<jtest>', 'sinwave-x-<jtest>', 'user'
                                ampluutest, &    ! amplitude(s) of initial condition for testflow solution (up to five different)
                                kx_uutest,  &    ! wavenumber of the initial condition (if sinusoidal)
                                lugu_as_aux      ! ??
  ! run parameters
!
  real :: nutest=0.,    &
          nutest1=0.,   &
          wamp=1.,      &
          ktestflow=1., &
          ktestflow1=1.
!
  logical :: reinitialize_uutest=.false., &
             lsoca_testflow=.true., &
             lkinem_testflow=.false., &
             lburgers_testflow=.false.
!
  character (len=labellen) :: itestflow='U11-U21'
!
  namelist /testflow_run_pars/ reinitialize_uutest, &    ! flag for reinitialization of testflow solution
                               zextent,             &    ! ??
                               lsoca_testflow,      &    ! flag for SOCA
                               lkinem_testflow,     &    ! flag for kinematic calculation
                               lburgers_testflow,   &    ! flag for disconnecting enthalpy(pressure) from velocity
                               nutest,              &    ! viscosity in testflow equations
                               nutest1,             &    ! reciprocal viscosity
                               itestflow,           &    ! name of used testflow set, legal values
                                                         ! 'W11-W22', 'quadratic', 'quadratic+G', 'U=0'
                               ktestflow,           &    ! wavenumber for testflow (if sinusoidal)
                               lugu_as_aux,         &    ! ??
                               duuinit,             &    !
                               wamp                      ! amplitude of the testflows
!
  ! diagnostic variables (needs to be consistent with reset list below)
!
  integer :: idiag_gal=0                            ! DIAG_DOC: GAL-coefficients, couple  $\overline F$ and $\overline U$
  integer :: idiag_nu=0                             ! DIAG_DOC: $\nu$-tensor,     couples $\overline F$ and $\partial \overline U/\partial z$
  integer :: idiag_aklam=0                          ! DIAG_DOC: $\lambda$-tensor, couples $\overline F$ and $\overline W$
  integer :: idiag_chi=0                            ! DIAG_DOC: $\chi$-vector,    couples $\overline F$ and ${\overline G}_z$
  integer :: idiag_kappa=0                          ! DIAG_DOC: $\kappa$-vector,  couples $\overline Q$ and $\overline U$
  integer :: idiag_psi=0                            ! DIAG_DOC: $\psi$-vector,    couples $\overline Q$ and $\overline W$
  integer :: idiag_pi=0                             ! DIAG_DOC: $\pi$-scalar,     couples $\overline Q$ and ${\overline G}_z$
!
  integer, dimension(2,2) :: idiag_galij=0, &       ! DIAG_DOC: $\alpha_{K,ij}$
                             idiag_aklamij=0, &     ! DIAG_DOC: $\lambda_{ij}$
                             idiag_nuij=0           ! DIAG_DOC: $\nu_{ij}$
!
  integer, dimension(2) :: idiag_chii=0, &          ! DIAG_DOC: $\chi_i$
                           idiag_kappai=0, &        ! DIAG_DOC: $\kappa_i$
                           idiag_psii=0             ! DIAG_DOC: $\psi_i$
!
  integer, dimension(3,njtestflow) :: idiag_Fipq=0      ! DIAG_DOC: ${\cal F}_i^{pq}$
  integer, dimension(njtestflow)   :: idiag_Qpq=0       ! DIAG_DOC: ${\cal Q}^{pq}$
  integer, dimension(0:njtestflow) :: idiag_upqrms=0    ! DIAG_DOC: $\left<{u^{pq}}^2\right>$
  integer, dimension(0:njtestflow) :: idiag_hpqrms=0    ! DIAG_DOC: $\left<{h^{pq}}^2\right>$
!
  integer :: idiag_ux0mz=0                          ! DIAG_DOC: $\left<u_{x}\right>_{xy}$
  integer :: idiag_uy0mz=0                          ! DIAG_DOC: $\left<u_{y}\right>_{xy}$
  integer :: idiag_uz0mz=0                          ! DIAG_DOC: $\left<u_{z}\right>_{xy}$
!
  integer :: nname_old
  integer, dimension(mname) :: idiag_map

  contains
!
!***********************************************************************
    subroutine register_testflow()          ! -> Standardschnittstelle
!
!  Initialise variables which should know that we solve for the vector
!  potential: iuutest, etc; increase nvar accordingly
!
!   3-jun-05/axel: adapted from register_magnetic
!
      use Cdata
      use Mpicomm
      use Sub
!
      integer :: j
!
!  Set first index of test flow in f
!
      iuutest=nvar+1
!
      ntestflow=4*(njtestflow+1)
!
      nvar=nvar+ntestflow
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_testflow: nvar = ', nvar
        print*, 'register_testflow: iuutest = ', iuutest
      endif
!
!  Put variable names in array
!
      do j=iuutest,nvar
        varname(j) = 'uutest'
      enddo
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_testflow: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',uutest $'
          if (nvar == mvar) write(4,*) ',uutest'
        else
          write(4,*) ',uutest $'
        endif
        write(15,*) 'uutest = fltarr(mx,my,mz,ntestflow)*one'
      endif
!
    endsubroutine register_testflow
!***********************************************************************
    subroutine initialize_testflow(f)           ! -> Standardschnittstelle
!
!  Perform any post-parameter-read initialization
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use FArrayManager
      use density, only: lcalc_glnrhomean, lupw_lnrho
      use Hydro, only: lcalc_uumean, lupw_uu
      use Forcing, only:ltestflow_forcing, lhydro_forcing
      use Viscosity, only:getnu
      use Mpicomm, only:stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
!
! Set nutest from the value uf nu in Hydro, if not specified here
!
      if ( nutest==0. .and. nutest1==0. ) &
        call getnu(nutest)
!
!  Precalculate nutest if 1/nutest (==nutest1) is given instead
!
      if (nutest1/=0.) then
        nutest=1./nutest1
      else if ( nutest==0. ) then
        !!call stop_it('duutest_dt: Error - nutest==0!!!')
      endif
!
!  set to zero and then rescale the testflow
!  (in future, could call something like init_uu_simple)
!
     if ( itestflow=='quadratic' ) then
       zq2 = 0.5*z*z
     else
!
!  set cosine and sine function for setting test fields and analysis
!
        cz=cos(ktestflow*z)
        sz=sin(ktestflow*z)
!
!  Also calculate the inverse of the testflow wavenumber, but only if different from zero
!
        if (ktestflow==0.) then
          ktestflow1=1.
        else
          ktestflow1=1./ktestflow
        endif
!
!  cosine and sine functions multiplied with k1
!
        k1cz=ktestflow1*cos(ktestflow*z)
        k1sz=ktestflow1*sin(ktestflow*z)
!
     endif
!
! ensure that hydro calculates mean velocity
!
     lcalc_uumean = .true.
!
! ensure that density calculates gradient of mean logdensity
!
     lcalc_glnrhomean = .true.
!
! ensure that for calculation of terms of the form v.grad(lnrho) upwinding is used as in density.f90
!
     ltestflow_upw_lnrho = lupw_lnrho
!
! ensure that for calculation of terms of the form (u_1.grad)u_2 upwinding is used as in hydro.f90
!
     ltestflow_upw_uu = lupw_uu
!
! ensure that forcing adds (discontinuous) force in f
!
     ltestflow_forcing = .true.
!
     if ( lkinem_testflow ) &
       lhydro_forcing = .false.
!
!  Register an extra aux slot for ugu if requested (so ugu is written
!  to snapshots and can be easily analyzed later). For this to work you
!  must reserve enough auxiliary workspace by setting, for example,
!     ! MAUX CONTRIBUTION 9
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
!
!  After a reload, we need to rewrite index.pro, but the auxiliary
!  arrays are already allocated and must not be allocated again.
!
      if (lugu_as_aux) then
        if (iugu==0) then
          call farray_register_auxiliary('ugu',iugu,vector=3*njtestflow)
        endif
        if (iugu/=0.and.lroot) then
          print*, 'initialize_magnetic: iugu = ', iugu
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'iugu=',iugu
          close(3)
        endif
      endif
!
!  write testflow information to a file (for convenient post-processing)
!
      if (lroot) then
!
        open(1,file=trim(datadir)//'/testflow_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'zextent=',merge(1,0,zextent)
        write(1,'(a,i1)') 'lsoca_testflow='  ,merge(1,0,lsoca_testflow)
        write(1,'(3a)') "itestflow='",trim(itestflow)//"'"
        write(1,'(a,f5.2)') 'ktestflow=',ktestflow
        close(1)
!
      endif
!
    endsubroutine initialize_testflow
!***********************************************************************
    subroutine init_uutest(f)          ! -> Standardschnittstelle
!
!  initialise testflow-solutions; called from start.f90
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use Mpicomm
      use Initcond
      use Sub
      use InitialCondition, only: initial_condition_uutest
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: j,ihhtest,jtest,offset
      character(len=labellen) selector
!
      jtest=1
!
      ihhtest=iuutest+3
!
!  initialize to zero
!
!AXEL f(:,:,:,iuutest:iuutest+ntestflow-1)=0.
!
      do j=1,ninit*njtest
!
!  loop over different choices of initial conditions
!
        if ( inituutest(j)(1:15).eq.'gaussian-noise-' ) then
!
          read( inituutest(j)(16:16), '(i1)' ) jtest
          selector='gaussian-noise'
!
        else if ( inituutest(j)(1:10).eq.'sinwave-x-' ) then
!
          read( inituutest(j)(11:11), '(i1)' ) jtest
          selector='sinwave-x'
        else if ( inituutest(j)(1:9).eq.'constant-' ) then
!
          read( inituutest(j)(10:10), '(i1)' ) jtest
          selector='constant'

        else
          selector = inituutest(j)
        endif
!
        offset = 4*jtest
!
        select case (selector)
!
        case ('zero'); f(:,:,:,iuutest:iuutest+ntestflow-1)=0.
        case ('constant'); f(:,:,:,iuutest+offset:ihhtest+offset)=ampluutest(j)
        case ('gaussian-noise'); call gaunoise(ampluutest(j),f,iuutest+offset,ihhtest+offset )
        case ('sinwave-x')
          call sinwave(ampluutest(j),f,iuutest+offset,kx=kx_uutest)
          call sinwave(ampluutest(j),f,ihhtest+offset,kx=kx_uutest)
!
        case ('nothing'); !(do nothing)
!
        case ('user')
!
!  Interface for user's own initial condition
!
          call initial_condition_uutest(f)
!
        case default
        !
        !  Catch unknown values
        !
          if (lroot) print*, 'init_uutest: check inituutest: ', trim(inituutest(j))
          call stop_it("")
!
        endselect
!
      enddo
!
    endsubroutine init_uutest
!***********************************************************************
    subroutine pencil_criteria_testflow()          ! -> Standardschnittstelle
!
!   All pencils that the Testflow module depends on are specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
      lpenc_requested(i_uu)=.true.
      lpenc_requested(i_glnrho)=.true.
      lpenc_requested(i_uij)=.true.
      lpenc_requested(i_sij)=.true.
!
    endsubroutine pencil_criteria_testflow
!***********************************************************************
    subroutine pencil_interdep_testflow(lpencil_in)         ! -> Standardschnittstelle
!
!  Interdependency among pencils from the Testflow module is specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
      logical, dimension(npencils) :: lpencil_in
!
    endsubroutine pencil_interdep_testflow
!***********************************************************************
    subroutine read_testflow_init_pars(unit,iostat)     ! -> Standardschnittstelle
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=testflow_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testflow_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_testflow_init_pars
!***********************************************************************
    subroutine write_testflow_init_pars(unit)         ! -> Standardschnittstelle
!
      integer, intent(in) :: unit
!
      write(unit,NML=testflow_init_pars)
!
    endsubroutine write_testflow_init_pars
!***********************************************************************
    subroutine read_testflow_run_pars(unit,iostat)    ! -> Standardschnittstelle
!
      integer, intent(in)              :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=testflow_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testflow_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_testflow_run_pars
!***********************************************************************
    subroutine write_testflow_run_pars(unit)          ! -> Standardschnittstelle
      integer, intent(in) :: unit
!
      write(unit,NML=testflow_run_pars)
!
    endsubroutine write_testflow_run_pars
!***********************************************************************
    subroutine duutest_dt(f,df,p)                    ! -> Standardschnittstelle
!
!  testflow evolution:
!
!  calculate du^(pq)/dt = -( u^(0) .grad u^(pq) + u^(pq).grad u ) + < u^(0) .grad u^(pq) + u^(pq).grad u >
!        [alternatively:  -( u^(pq).grad u^(0)  + u.grad u^(pq) ) + < u^(pq).grad u^(0)  + u.grad u^(pq) >]
!                                                 - grad h^(pq)
!                         -U^(pq).grad u - u.grad U^(pq)
!                         + viscous terms
!                         + Lorentz force (added in testfield)
!
!  and dh^(pq)/dt = -( u^(0).grad h^(pq) + u^(pq).grad h ) + < u^(0).grad h^(pq) + u^(pq).grad h >
!  [alternatively:  -( u^(pq).grad h^(0) + u.grad h^(pq) ) + < u^(pq).grad h^(0) + u.grad h^(pq) > ]
!                                       - cs2*div u^(pq)
!                   -U.grad h -u.grad H
!
!
!  12-mar-08/axel: coded
!  24-jun-08/MR: modified
!
      use Cdata
      use Diagnostics
      use Sub
      use Mpicomm, only: stop_it
      use density, only: glnrhomz
      use Hydro, only: uumz, guumz, traceless_strain, coriolis_cartesian, &
                       lforcing_cont_uu, ampl_fcont_uu
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)    :: f,p
      intent(inout) :: df
!
! local variables
!
      real, dimension (nx,3) :: uutest,uufluct,del2utest,graddivutest,ghtest,ghfluct,U0testgu,ugU0test
      real, dimension (nx,3,3) :: uijfluct
      real, dimension (nx,3) :: U0test,gU0test                      !!!MR: überflüssige x-Ausdehnung bei 0-Größen
      real, dimension (nx)   :: hhtest,U0ghtest,upq2,divutest, divufluct
      real :: gH0test
!
      logical :: ltestflow_out
!
      integer :: jtest,j,i,i3,i4, nl
      integer,save :: ifirst=0
      integer :: iuxtest, iuytest, iuztest, ihhtest
!
      character :: cjtest
      character (len=5) :: ch
      character (len=130) :: file

      if (n==4.and.m==4) then
      !!print*, 'uumz(1):', uumz(:,1)
      !!print*, 'guumz(1):', guumz(:,1)
      endif
!
      if ( iuutest.eq.0 ) return

      U0test=0.; gU0test=0.; gH0test=0.
!
!  identify module
!
      if (headtt.or.ldebug) print*,'duutest_dt: SOLVE'
!
      iuxtest=iuutest
      iuytest=iuutest+1
      iuztest=iuutest+2
      ihhtest=iuutest+3
!
      if (lkinem_testflow) then
!
        uufluct = f(l1:l2,n,m,iuxtest:iuztest)           ! ufluct = u0
        call grad(f,ihhtest,ghfluct)                     ! ghfluct = grad(h0)
!
      else
!
!  calculate uufluct=U-Umean out of the main run
!
        do j=1,3
          uufluct(:,j)=p%uu(:,j)-uumz(n-n1+1,j)
        enddo
!
! calculate ghfluct=gh-ghmean
!
        do j=1,3
          ghfluct(:,j)=p%glnrho(:,j)-glnrhomz(n-n1+1,j)
        enddo

      endif
!
!  do each of the njtestflow test flows at a time
!  jtest=0  refers to primary turbulence (u^(0), h^(0))
!
      do jtest=0,njtestflow
!
!  identify boundary conditions
!
        if (headtt) then
!
          write(cjtest,'(i1)') jtest
          call identify_bcs('uxtest-'//cjtest,iuxtest)
          call identify_bcs('uytest-'//cjtest,iuytest)
          call identify_bcs('uztest-'//cjtest,iuztest)
          call identify_bcs('hhtest-'//cjtest,ihhtest)
!
        endif
!
!  velocity vector and enthalpy
!
        uutest=f(l1:l2,m,n,iuxtest:iuztest)	! testflow solution # jtest
        hhtest=f(l1:l2,m,n,ihhtest)
        !!if ( jtest.eq.1 ) print*, 'uutest, hhtest:', minval(uutest),maxval(uutest), minval(hhtest),maxval(hhtest)
!
!  velocity gradient matrix and div u term
!
        call div(f,iuxtest,divutest)
!
!  gradient of (pseudo) enthalpy
!
        call grad(f,ihhtest,ghtest)
!
        if (jtest.gt.0) then
!
          select case (itestflow)                       ! get testflow U^(pq), grad U^(pq), grad H^(pq)
!
            case ('W11-W22')
              call set_U0test_W11_W22(U0test,gU0test,jtest)
              gH0test=0.
!
            case ('onlyconstant')
!
              call set_U0test_onlyconstant(U0test,gU0test,jtest)
              gH0test=0.
!
            case ('onlylinear')
!
              call set_U0test_onlylinear(U0test,gU0test,jtest)
              gH0test=0.
!
            case ('quadratic')
              call set_U0test_quadratic(U0test,gU0test,jtest)
              gH0test=0.
!
            case ('quasiperiodic')
              call set_U0test_quasiperiodic(U0test,gU0test,jtest)
              gH0test=0.
!
            case ('quadratic+G')
!
              if ( jtest.le.6 ) then
                call set_U0test_quadratic(U0test,gU0test,jtest)
                gH0test=0.
              else
                U0test=0.; gU0test=0.; gH0test=wamp
              endif
!
            case ('U=0')
              U0test=0.; gU0test=0.; gH0test=wamp

            case default
              call fatal_error('duutest_dt','undefined itestflow value "'//itestflow//'"')
!
          endselect
!
        endif ! otherwise = for primary turbulence, equal to zero
!
!  rhs of continuity equation (nonlinear terms already in df, if needed!),
!  dh^pq/dt = n.l.Terms - cs^2*div u^pq -u.gradH^pq - U^pq.gradh
!
        if ( .not.lburgers_testflow ) &
          df(l1:l2,m,n,ihhtest)=df(l1:l2,m,n,ihhtest) - cs2test*divutest
!
!       testflow inhomogeneity
!
        if (jtest.gt.0) then
!
          call dot_mn(U0test,ghfluct,U0ghtest)          ! MR: checked by numbers
          !!call u_dot_grad(f,iuutest+3, ghfluct, U0test, U0ghtest,UPWIND=ltestflow_upw_lnrho)   ! Utest.grad(h)!!!MR: noch falsch, da unter ilnrho nicht hhfluct!!
                          !!!ilnrho
          !!if ( jtest.eq.2 .and. ( maxval(ghfluct(:,2)-U0ghtest(:))/=0. .or. minval(ghfluct(:,2)-U0ghtest(:))/=0. ) ) &
            !!print*, 'U0ghtest, n, m:', n, m, ghfluct(:,2)-U0ghtest(:)
!
          if ( .not.lburgers_testflow ) &
            df(l1:l2,m,n,ihhtest)=df(l1:l2,m,n,ihhtest)-uufluct(:,3)*gH0test-U0ghtest

        endif
!
!  rhs of momentum equation (nonlinear terms already in df, if needed!),
!  du^pq/dt = n.l.Terms - grad h^pq -u.gradU^pq - U^pq.gradu
!
        if ( .not.lburgers_testflow ) &
          df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) - ghtest
!
!       testflow inhomogeneity
!
        if ( jtest.gt.0 ) then
!
          if (lkinem_testflow) then
!
            call gij(f,iuxtest,uijfluct,1)
!
          else
!
            uijfluct = p%uij
!
            do i=1,3                                        ! determine fluctuating part of uij
              uijfluct(:,i,3) = uijfluct(:,i,3) - guumz(n-n1+1,i)
            enddo
!
          endif
!
          call h_dot_grad(U0test,uijfluct,uufluct,U0testgu)
                                          !this parameter without meaning in cartesian geometry
!
          call multsv(uufluct(:,3),gU0test,ugU0test)
!
          df(l1:l2,m,n,iuxtest:iuztest) = df(l1:l2,m,n,iuxtest:iuztest) - U0testgu-ugU0test
!
!         viscous contributions
!
          if (nutest/=0.) then
!
            if (lkinem_testflow) then

              call div_mn( uijfluct, divufluct, uufluct )
              call traceless_strain( uijfluct, divufluct, uijfluct, uufluct )
                                                         !this parameter without meaning in cartesian geometry
            else

              uijfluct = p%sij                                               ! uijfluct now contains sij !!!

              nl = n-n1+1
              do i=1,2                                                       ! determine fluctuating part
                uijfluct(:,i,i) = uijfluct(:,i,i) + guumz(nl,3)/3.           ! ~
              enddo
!
              uijfluct(:,3,3) = uijfluct(:,3,3) - 2.*guumz(nl,3)/3.          ! ~
              uijfluct(:,1,3) = uijfluct(:,1,3) - 0.5*guumz(nl,1)            ! ~
              uijfluct(:,3,1) = uijfluct(:,1,3)                              ! ~
              uijfluct(:,2,3) = uijfluct(:,2,3) - 0.5*guumz(nl,2)            ! ~
              uijfluct(:,3,2) = uijfluct(:,2,3)                              ! of sij
!
            endif
!
            U0testgu = uijfluct(:,:,3)*gH0test                               ! S(u).grad(H^T)
!
            call multsv(0.5*ghfluct(:,3),gU0test,ugU0test)             ! S(U^T).grad(h)        MR: checked by numbers
            ugU0test(:,3) = ugU0test(:,3) + 0.5*ghfluct(:,3)*gU0test(:,3)    !  ~
!
            call multsv(-(1./3.)*gU0test(:,3),ghfluct,ugU0test,.true.)       !  ~
!
            if ( .not.lburgers_testflow ) &
              df(l1:l2,m,n,iuxtest:iuztest) =   df(l1:l2,m,n,iuxtest:iuztest) &
                                              + 2.*nutest*cs2test1*(U0testgu+ugU0test)
          endif

        endif
!
!  add linear part of diffusion term nu*(del2u^pq + divu^pq/3)
!
        if (nutest/=0.) then
!
          call del2v_etc(f,iuxtest,DEL2=del2utest,GRADDIV=graddivutest)
          df(l1:l2,m,n,iuxtest:iuztest) =   df(l1:l2,m,n,iuxtest:iuztest) &
                                          + nutest*(del2utest+(1./3.)*graddivutest)
        endif
!
! add continuous forcing
!
        if ( jtest==0.and.lforcing_cont_uu) &                                 ! requires that forcing is fluctuation
          df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest)+&
                                        ampl_fcont_uu*p%fcont
!
! add Coriolis/centrifugal force         (MR: precession still missing!!!)
!
        call coriolis_cartesian(df,uutest,iuxtest)
!
!  check for testflow timestep
!
        if (lfirst.and.ldt) &
          advec_uu=max(advec_uu,abs(uutest(:,1))*dx_1(l1:l2)+ &
                                abs(uutest(:,2))*dy_1(  m  )+ &
                                abs(uutest(:,3))*dz_1(  n  ))
        iuxtest=iuxtest+4
        iuztest=iuztest+4
        ihhtest=ihhtest+4
!
      enddo
!
!  diffusive time step, just take the max of diffus_nu (if existent)
!  and whatever is calculated here. Check also for testsound timestep.
!
      if (lfirst.and.ldt) then
        advec_cs2=max(advec_cs2,cs2test*dxyz_2)
        diffus_nu=max(diffus_nu,nutest*dxyz_2)
      endif
!
      if (idiag_ux0mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,iuutest  ),idiag_ux0mz)!MR: only testflow # 0
      if (idiag_uy0mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,iuutest+1),idiag_uy0mz)
      if (idiag_uz0mz/=0) call xysum_mn_name_z(f(l1:l2,m,n,iuutest+2),idiag_uz0mz)
!
!  rms values of small scales fields upq in response to the test fields Upq
!  Obviously idiag_u0rms and idiag_u12rms cannot both be invoked!
!  Needs modification!
!
      if (ldiagnos) then
!
        iuxtest = iuutest
        iuztest = iuutest+2
        ihhtest = iuutest+3
!
        do jtest=0,njtestflow
!
          if (idiag_upqrms(jtest)/=0) then
!
            call dot2(f(l1:l2,m,n,iuxtest:iuztest),upq2)
            !!if (jtest==1) print*, 'upq,jtest:', jtest,minval(upq2),maxval(upq2)
!
            call sum_mn_name(upq2,idiag_upqrms(jtest),lsqrt=.true.)

          endif
!
          if (idiag_hpqrms(jtest)/=0) then
!
            upq2=f(l1:l2,m,n,ihhtest)*f(l1:l2,m,n,ihhtest)
            call sum_mn_name(upq2,idiag_hpqrms(jtest),lsqrt=.true.)

          endif

          iuxtest = iuxtest+4
          iuztest = iuztest+4
          ihhtest = ihhtest+4
!
        enddo
      endif
!
!  write utest-slices for output in wvid in run.f90
!  Note: ix is the index with respect to array with ghost zones.
!
      if (lvideo.and.lfirst) then
        do j=1,3
          uu11_yz(m-m1+1,n-n1+1,j)             = f(ix_loc-l1+1,m,n,iuutest+j-1)        !MR: only testflow # 0
          if (m==iy_loc)  uu11_xz(:,n-n1+1,j)  = f(l1:l2      ,m,n,iuutest+j-1)
          if (n==iz_loc)  uu11_xy(:,m-m1+1,j)  = f(l1:l2      ,m,n,iuutest+j-1)
          if (n==iz2_loc) uu11_xy2(:,m-m1+1,j) = f(l1:l2      ,m,n,iuutest+j-1)
        enddo
      endif
!
! initialize uutest periodically if requested
!
      if (reinitialize_uutest) then
!
        file=trim(datadir)//'/tinit_uutest.dat'
!
        if (ifirst==0) then
!
          call read_snaptime(trim(file),tuuinit,nuuinit,duuinit,t)
!
          if (tuuinit==0 .or. tuuinit < t-duuinit) then
            tuuinit=t+duuinit
          endif
          ifirst=1
!
        endif
!
        if (t >= tuuinit) then

          call init_uutest(f) !MR: ist nicht das gemeint???
          !!!call initialize_testflow(f,reinitialize_uutest)
          call update_snaptime(file,tuuinit,nuuinit,duuinit,t,ltestflow_out,ch,.false.)
!
        endif
      endif
!
    endsubroutine duutest_dt
!***********************************************************************
    subroutine get_slices_testflow(f,slices)          ! -> Standardschnittstelle
!
!  Write slices for animation of velocity variables.
!
!  12-sep-09/axel: adapted from the corresponding magnetic routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Velocity field
!
        case ('uu11')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =>uu11_yz(:,:,slices%index)
            slices%xz =>uu11_xz(:,:,slices%index)
            slices%xy =>uu11_xy(:,:,slices%index)
            slices%xy2=>uu11_xy2(:,:,slices%index)
            if (slices%index<=3) slices%ready=.true.
          endif
      endselect
!
    endsubroutine get_slices_testflow
!***********************************************************************
    subroutine calc_ltestflow_nonlin_terms(f,df)          ! -> Standardschnittstelle (als do_prepencilstep oder so)
!
!  calculates < -u0.grad(u0) + (2nu/cs^2)*grad(h0).Sij(u0) >, < u0.gradh0 >,
!             < -u0.grad(utest) - utest.grad(u) + (2nu/cs^2)*( grad(h0).Sij(utest) + grad(htest).Sij(u) ) >,
!             <  u0.grad(htest) + utest.grad(h) >
!
!  which is needed when lsoca_testflow=.false., resp.
!  this is done prior to the pencil loop (calc_fluct=false)
!
!  calculates ( -u0.grad(u0) + (2nu/cs^2)*grad(h0).Sij(u0) )', ( u0.gradh0 )',
!             ( -u0.grad(utest) - utest.grad(u) ) + (2nu/cs^2)*( grad(h0).Sij(utest) + grad(htest).Sij(u) )',
!             (  u0.grad(htest) + utest.grad(h) )'
!
!  which is needed when lsoca_unl=.false.,  resp.
!  this is done inside the pencil loop
!
!  15-mar-08/axel: coded
!  24-jun-08/MR: modified
!  01-sep-09/MR: further processed
!
      use Cdata
      use Sub
      use Diagnostics
      use density, only: glnrhomz
      use Hydro, only: uumz, traceless_strain
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension (mx,my,mz,mfarray), intent(in)    :: f
      real, dimension (mx,my,mz,mvar),    intent(inout) :: df
!
      real, dimension (3,0:njtestflow) :: unltestm, unltestm1      ! unltestm1, hnltestm1 local fields for MPI correspondence
      real, dimension (0:njtestflow) :: hnltestm, hnltestm1
!
      real, dimension (nx,3)   :: uufluct,uutest, uu0, ghfluct, ghtest, gh0, sghtest, unltest
      real, dimension (nx,3,3) :: sijtest,uijtest,sij0,uij0
      real, dimension (nx)     :: divutest,hnltest, divufluct

      integer :: jtest,i,j,ju,k,iii, jtesto
      integer :: iuxtest, iuztest, ihhtest
      logical :: headtt_save
      real :: fac, gal1, gal2, aklam1, aklam2

      fac=1./(nxgrid*nygrid)
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt

      lfirstpoint=.true.
!
zloop:do n=n1,n2

        unltestm=0.
        hnltestm=0.
!
yloop:  do m=m1,m2
!
          iuxtest=iuutest
          iuztest=iuutest+2
          ihhtest=iuutest+3
!
          if (lkinem_testflow) then
!
            uufluct = f(l1:l2,n,m,iuxtest:iuztest)           ! ufluct = u0
            call grad(f,ihhtest,ghfluct)                     ! ghfluct = grad(h0)
!
          else
!
!  calculate uufluct=U-Umean out of the main run
!
            uufluct=f(l1:l2,m,n,iux:iux+2)
!
!!            if ( n.eq.5.and.m.eq.5 ) print*,'nl:uufluct:', maxval(uufluct), minval(uufluct)
            do j=1,3
              uufluct(:,j)=uufluct(:,j)-uumz(n-n1+1,j)
            enddo
!
! calculate ghfluct=gh-ghmean out of the main run
!
            call grad(f,ilnrho,ghfluct)
            do j=1,3
              ghfluct(:,j)=ghfluct(:,j)-glnrhomz(n-n1+1,j)
            enddo
!
          endif
!
!  do each of the njtestflow (=2, 4, or 6) test flows at a time
!
testloop: do jtest=0,njtestflow                           ! jtest=0 : primary turbulence, testflow=0
!
!  velocity vector and enthalpy gradient of testflow solution
!
            uutest=f(l1:l2,m,n,iuxtest:iuztest)
            call grad(f,ihhtest,ghtest)                 ! grad(htest)
!
!  velocity gradient matrix and velocity divergence of testflow solution
!
            call gij(f,iuxtest,uijtest,1)               !!p
            call div_mn(uijtest,divutest,uutest)        !!p
!
!  calculate traceless strain tensor sijtest from uijtest

            call traceless_strain( uijtest, divutest, sijtest, uutest )
                                                               ! this parameter is without meaning in Cartesian geometry
            if ( jtest.eq.0 ) then    ! primary turbulence

              gh0  = ghtest           ! save primary turbulence
              uu0  = uutest
              uij0 = uijtest
              sij0 = sijtest
!
!  u.grad(u) term and u.grad(h) term
!
              !!call multmv(uijtest,uutest,unltest)
              call u_dot_grad(f,iuxtest,uijtest,uutest,unltest,UPWIND=ltestflow_upw_uu)     ! (u0.grad)(u0)

              !!call dot_mn(uutest,ghtest,hnltest)                                          ! u0.grad(h0)
              if ( .not.lburgers_testflow ) &
                call u_dot_grad(f,ihhtest,ghtest,uutest,hnltest,UPWIND=ltestflow_upw_lnrho)   ! u0.grad(h0)

            else
!
              !!call multmv(uijtest,uufluct,unltest)                                        ! (u.grad)(utest)
              call u_dot_grad(f,iuxtest,uijtest, uufluct,unltest,UPWIND=ltestflow_upw_uu)
!
              !!call multmv(uij0,uutest,unltest,.true.)                                     ! (utest.grad)(u0)
              call u_dot_grad(f,iuutest,uij0,uutest,unltest,ltestflow_upw_uu,.true.)
!
              if ( .not.lburgers_testflow ) then
!
                !!call dot_mn(uufluct,ghtest ,hnltest)                                        ! u.grad(htest)
                call u_dot_grad(f,ihhtest,ghtest,uufluct,hnltest,ltestflow_upw_lnrho)         ! u.grad(htest)
!
                !!call dot_mn(uutest,gh0,hnltest,.true.)                                      ! utest.grad(h0)
                call u_dot_grad(f,iuutest+3,gh0,uutest,hnltest,ltestflow_upw_lnrho,.true.)    ! utest.grad(h0)
!
              endif
!
            endif
!
!  add nonlinear part of diffusion term nu*2*S.grad(h)/cs2
!
            if (nutest/=0.) then
!
              if ( jtest.eq.0 ) then
!
                call multmv(sijtest,ghtest,sghtest)                            ! S(u0).grad(h0)
!
              else

                call multmv(sijtest,ghfluct,sghtest)                           ! S(utest).grad(h)
                call multmv(sij0   ,ghtest ,sghtest,.true.)                    ! S(u0).grad(htest)
!
              endif
!
              if ( .not.lburgers_testflow ) &
                unltest = unltest - (2.*nutest*cs2test1)*sghtest
!
            endif
!
            if ( jtest.eq.0 .or. .not.lsoca_testflow ) then

              df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest)-unltest      ! nonlinear parts stored in df
!
              if ( .not.lburgers_testflow ) &
                df(l1:l2,m,n,ihhtest        )=df(l1:l2,m,n,ihhtest        )-hnltest

            endif

            if ( jtest.eq.0 .or. .not.lsoca_testflow .or. ldiagnos ) then  ! means of nonlinear terms always needed for 0th testflow,
                                                                           ! if not SOCA for all testflows, except from that if the
                                                                           ! the turbulence coefficients have to be calculated
              do j=1,3
                unltestm(j,jtest)=unltestm(j,jtest)+sum(unltest(:,j))      ! sums of nonlinear parts (here sum sums over x extent only)
              enddo

              hnltestm(jtest)=hnltestm(jtest)+sum(hnltest)

            endif
!
            iuxtest=iuxtest+4
            iuztest=iuztest+4
            ihhtest=ihhtest+4
!
          enddo testloop
!
!  reset headtt
!
          headtt=.false.           !???MR: richtige Stelle?
!
        enddo yloop
!
!  do communication for arrays of size 3*njtestflow and njtestflow, resp.
!
        if ( .not.lsoca_testflow .or. ldiagnos ) then         ! see above
          jtesto=njtestflow
        else
          jtesto=0
        endif
!
        if (nprocy>1) then

          call mpiallreduce_sum(unltestm,unltestm1,(/3,jtesto+1/),idir=2)
          unltestm(:,0:jtesto) = unltestm1(:,0:jtesto)
!
          call mpiallreduce_sum(hnltestm,hnltestm1,jtesto+1,idir=2)
          hnltestm(0:jtesto) = hnltestm1(0:jtesto)
!
        endif
!
        unltestm(:,0:jtesto)=fac*unltestm(:,0:jtesto)         ! means of nonlinear parts
        hnltestm(0:jtesto)=fac*hnltestm(0:jtesto)
!
!  means are completely determined -> calculation of the  f l u c t u a t i o n s  of the nonlinear parts
!
        if ( lsoca_testflow ) then
          jtesto=0
        else
          jtesto= njtestflow
        endif
!
        do jtest=0,jtesto
!
          iuxtest=iuutest+4*jtest
          ihhtest=iuxtest+3

          do j=1,3
!
            ju = iuxtest+j-1
!
            df(l1:l2,m1:m2,n,ju)=df(l1:l2,m1:m2,n,ju)+unltestm(j,jtest)
!
          enddo

          df(l1:l2,m1:m2,n,ihhtest)=df(l1:l2,m1:m2,n,ihhtest)+hnltestm(jtest)

        enddo
!
        !!print*, 'unltestm, hnltestm:', minval(unltestm),maxval(unltestm), minval(hnltestm),maxval(hnltestm)
        !!if (ldiagnos) print*, 'unltestm, hnltestm:', z(n), unltestm(1,1), unltestm(2,1), unltestm(1,2), unltestm(2,2)

        if (ldiagnos) call calc_coeffcients(n,unltestm,hnltestm)
!
        lfirstpoint=.false.
!
      enddo zloop
!
! remove gal-parts from Fipq
!
      if (ldiagnos) then
!
        do k=1,2

          gal1 = get_from_fname(idiag_galij(k,1))
          gal2 = get_from_fname(idiag_galij(k,2))

          do n=n1,n2
!
            select case (itestflow)

            case ('W11-W22')
!
              if (idiag_aklamij(k,1)/=0) &
                call surf_mn_name(  wamp*(-cz(n)*gal2*k1sz(n)+sz(n)*gal2*k1cz(n)),      idiag_aklamij(k,1) )

              if (idiag_nuij(k,2)/=0) &
                call surf_mn_name( ( wamp*(k1sz(n)*gal2*k1sz(n)+k1cz(n)*gal2*k1cz(n))), idiag_nuij(k,2)    )
!
              if (idiag_aklamij(k,2)/=0) &
                call surf_mn_name(  wamp*(+cz(n)*gal1*k1sz(n)-sz(n)*gal1*k1cz(n)),      idiag_aklamij(k,2) )
!
              if (idiag_nuij(k,1)/=0) &
                call surf_mn_name( wamp*(+k1sz(n)*gal1*k1sz(n)+k1cz(n)*gal1*k1cz(n)), idiag_nuij(k,1)    )

            case ('quadratic')
!
              if (idiag_aklamij(k,1)/=0) &
                  call surf_mn_name( -wamp*gal2*z(n), idiag_aklamij(k,1) )
!
              if (idiag_aklamij(k,2)/=0) &
                  call surf_mn_name(  wamp*gal1*z(n), idiag_aklamij(k,2) )
            case default
!
            end select
!
          enddo
!
          if ( itestflow=='quadratic') then

            aklam1 = get_from_fname(idiag_aklamij(k,1))
            aklam2 = get_from_fname(idiag_aklamij(k,2))
!
            do n=n1,n2
!
              if (idiag_nuij(k,1)/=0) &
                call surf_mn_name( wamp*(-gal1*zq2(n)+aklam2*z(n)), idiag_nuij(k,1) )
!
              if (idiag_nuij(k,2)/=0) &
                call surf_mn_name( wamp*(-gal2*zq2(n)-aklam1*z(n)), idiag_nuij(k,2) )
!
            enddo
!
          endif
!
        enddo
!
      endif
!
      headtt=headtt_save
!
    endsubroutine calc_ltestflow_nonlin_terms
!***********************************************************************
    subroutine calc_coeffcients(indz,Fipq,Qipq)
!
    use Diagnostics
    use Cdata
!
    implicit none
!
    integer, intent(in) :: indz
!
    real, dimension (3,0:njtestflow) :: Fipq            ! double index pq in single one (second) subsumed
    real, dimension (0:njtestflow) :: Qipq
!
    real, dimension (2,2) :: aklam
    integer :: i,j,i3,i4,i5,i6,k
!
      Fipq=Fipq/(wamp*nzgrid)
      Qipq=Qipq/(wamp*nzgrid)                               ! factor nzgrid for averaging over z
!
      !!print*,'z,Fipq=', z(indz), Fipq(1:2,1), Fipq(1:2,2)
!
      do j=1,njtestflow
!
        do i=1,3
          if (idiag_Fipq(i,j)/=0) call surf_mn_name( Fipq(i,j), idiag_Fipq(i,j) )   ! surf_mn_name because only simple addition needed
        enddo
!
        if (idiag_Qpq(j)/=0) call surf_mn_name( Qipq(j), idiag_Qpq(j) )
!
      enddo
!
!  calculate gal, aka-lambda and nu tensors
!
      if ( idiag_aklam/=0 .and. njtestflow.lt.4 ) then
!
        print*,'Warning: Number of testflows njtestflow=',njtestflow,' is too small for AKA-Lambda tensor to be calculated'
        idiag_aklam=0
!
      endif
!
      if ( idiag_nu/=0 .and. njtestflow.lt.6 ) then
!
        print*,'Warning: Number of testflows njtestflow=',njtestflow,' is too small for viscosity tensor to be calculated'
        idiag_nu=0
!
      endif

      select case (itestflow)

        case ('W11-W22')
!
          if (idiag_gal/=0) then
!
            do i=1,2
            do j=1,2
              call surf_mn_name(Fipq(i,j+4),idiag_galij(i,j))          ! \gal_{ij}
            enddo
            enddo
!
          endif
!
          do k=1,2
!
            if (idiag_aklamij(k,1)/=0) &
              call surf_mn_name( -cz(indz)*Fipq(k,1)-sz(indz)*Fipq(k,2),     idiag_aklamij(k,1) )
!
            if (idiag_aklamij(k,2)/=0) &
              call surf_mn_name( -cz(indz)*Fipq(k,3)-sz(indz)*Fipq(k,4),     idiag_aklamij(k,2) )
!
            if (idiag_nuij(k,1)/=0) &
              call surf_mn_name( -k1sz(indz)*Fipq(k,3)+k1cz(indz)*Fipq(k,4), idiag_nuij(k,1)    )
!
            if (idiag_nuij(k,2)/=0) &
              call surf_mn_name(  k1sz(indz)*Fipq(k,1)-k1cz(indz)*Fipq(k,2), idiag_nuij(k,2)    )
!
          enddo
          !!print*, 'indz=', indz, -k1sz(indz)*Fipq(1,3)+k1cz(indz)*Fipq(1,4), k1sz(indz)*Fipq(2,1)-k1cz(indz)*Fipq(2,2), &
          !!                       -k1sz(indz)*Fipq(2,3)+k1cz(indz)*Fipq(2,4), k1sz(indz)*Fipq(1,1)-k1cz(indz)*Fipq(1,2)
!
        case ('quadratic')
!
          if (idiag_gal/=0) then
!
            do i=1,2
            do j=1,2
              call surf_mn_name(Fipq(i,j),idiag_galij(i,j))          ! \gal_{ij}
            enddo
            enddo
!
!!            if ( indz>=n1+2 .and. indz<=n2-2 ) then
            if (idiag_aklam/=0) then
!
              i3=3
              i4=4

              aklam(1,1) =  Fipq(1,i4)
              aklam(1,2) = -Fipq(1,i3)
              aklam(2,1) =  Fipq(2,i4)
              aklam(2,2) = -Fipq(2,i3)
              !!print*, 'lam=', idiag_gal, idiag_aklam, idiag_nu

              do i=1,2
              do j=1,2
                call surf_mn_name( aklam(i,j), idiag_aklamij(i,j)  )  ! \aklam_{ij}
              enddo
              enddo
!
              if (idiag_nu/=0) then

                i5=5
                i6=6

                call surf_mn_name( Fipq(1,i5), idiag_nuij(1,1) ) ! \nu_{11}
                call surf_mn_name( Fipq(1,i6), idiag_nuij(1,2) ) ! \nu_{12}
                call surf_mn_name( Fipq(2,i5), idiag_nuij(2,1) ) ! \nu_{21}
                call surf_mn_name( Fipq(2,i6), idiag_nuij(2,2) ) ! \nu_{22}
!
                !!print*, 'nu11,22=', z(indz), Fipq(1,i5), Fipq(2,i6)
              endif
            endif
!!            endif
          endif

          if ( idiag_kappa/=0) then
!
            do i=1,2
              call surf_mn_name( Qipq(i), idiag_kappa+i-1 )        ! \kappa_i
            enddo
!
            if ( idiag_psi/=0) then
!
              do i=1,2
                call surf_mn_name( -Qipq(i+2)+Qipq(i)*z(indz), idiag_psi+i-1 )        ! \psi_i
              enddo
            endif
!
          endif
!
        case ('onlyconstant')
!
          do i=1,2
            do j=1,2
              call surf_mn_name(Fipq(i,j),idiag_galij(i,j))          ! \gal_{ij}
            enddo
          enddo
!
        case ('onlylinear')
!
          if (idiag_aklam/=0) then
!
            i3=3
            i4=4

            aklam(1,1) =  Fipq(1,i4)
            aklam(1,2) = -Fipq(1,i3)
            aklam(2,1) =  Fipq(2,i4)
            aklam(2,2) = -Fipq(2,i3)
            !!print*, 'lam=', idiag_gal, idiag_aklam, idiag_nu

            do i=1,2
              do j=1,2
                call surf_mn_name( aklam(i,j), idiag_aklamij(i,j)  )  ! \aklam_{ij}
              enddo
            enddo
!
          endif
        case ('quasi-periodic')
!
        case ('W=0')

        case default
              call fatal_error('duutest_dt','undefined itestflow value')
      endselect
!
!  print warning if aklam12 and aklam21 are needed, but njtestflow is too small
!
!      if ((idiag_aklamij(1,2)/=0.or.idiag_aklamij(2,2)/=0 &
!         .or.idiag_nuij(1,2)/=0.or.idiag_nuij(2,2)/=0).and.njtestflow<=2) then
!        call stop_it('njtestflow is too small if aklam12, aklam22, nu12, or nu22 are needed')
!      else
!        if (idiag_aklamij(1,2)/=0) call sum_name(+cz(iz)*Fipq(1,3)+sz(indz)*Fipq(1,4),idiag_aklamij(1,2))
!
    endsubroutine calc_coeffcients
!***********************************************************************
    subroutine set_uutest(uutest,jtest)
!
!  set testflow
!
!   3-jun-05/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (nx,3) :: uutest
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: uutest
!
!  set uutest for each of the 9 cases
!
      select case (jtest)
!
      case (1); uutest(:,1)=cz(n); uutest(:,2)=0.; uutest(:,3)=0.
      case (2); uutest(:,1)=sz(n); uutest(:,2)=0.; uutest(:,3)=0.
      case (3); uutest(:,1)=0.   ; uutest(:,2)=0.; uutest(:,3)=0.
      case default; uutest=0.
!
      endselect
!
    endsubroutine set_uutest
!***********************************************************************
    subroutine set_uutest_U11_U21 (uutest,jtest)
!
!  set testflow
!
!   3-jun-05/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (nx,3) :: uutest
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: uutest
!
!  set uutest for each of the 9 cases
!
      select case (jtest)
      case (1); uutest(:,1)=cz(n); uutest(:,2)=0.; uutest(:,3)=0.
      case (2); uutest(:,1)=sz(n); uutest(:,2)=0.; uutest(:,3)=0.
      case default; uutest=0.
      endselect
!
    endsubroutine set_uutest_U11_U21
!***********************************************************************
    subroutine set_U0test_W11_W22 (U0test,gU0test,jtest)
!
!  set testflow
!
!  23-mar-08/axel: adapted from testflow_z
!
      use Cdata
      use Sub
!
      real, dimension (nx,3), intent(out) :: U0test,gU0test
      integer,                intent(in)  :: jtest
!
!  set U0test and gU0test for each of the various cases
!
      select case (jtest)
!
      case (1)
        U0test(:,1)=0.; U0test(:,2)=-wamp*k1sz(n); U0test(:,3)=0.
        gU0test(:,1)=0.; gU0test(:,2)=-wamp*cz(n); gU0test(:,3)=0.
!
      case (2)
        U0test(:,1)=0.; U0test(:,2)=+wamp*k1cz(n); U0test(:,3)=0.
        gU0test(:,1)=0.; gU0test(:,2)=-wamp*sz(n); gU0test(:,3)=0.
!
      case (3)
        U0test(:,1)=+wamp*k1sz(n); U0test(:,2)=0.; U0test(:,3)=0.
        gU0test(:,1)=+wamp*cz(n); gU0test(:,2)=0.; gU0test(:,3)=0.
!
      case (4)
        U0test(:,1)=-wamp*k1cz(n); U0test(:,2)=0.; U0test(:,3)=0.
        gU0test(:,1)=+wamp*sz(n); gU0test(:,2)=0.; gU0test(:,3)=0.
!
      case (5)
        U0test(:,1)=wamp; U0test(:,2)=0.; U0test(:,3)=0.
        gU0test=0.
!
      case (6)
        U0test(:,1)=0.; U0test(:,2)=wamp; U0test(:,3)=0.
        gU0test=0.
!
     case default
        U0test=0.;gU0test=0.
!
      endselect
!
    endsubroutine set_U0test_W11_W22
!***********************************************************************
    subroutine set_U0test_quasiperiodic(U0test,gU0test,jtest)
!
!  set testflow
!
!  27-aug-09/MR: adapted from set_U0test_W11_W22
!
      use Cdata
      use Sub
!
      real, dimension (nx,3), intent(out) :: U0test,gU0test
      integer,                 intent(in) :: jtest
!
!  set U0test and gU0test for each of the testflow indices jtest (which stands for p and q)
!  order convention: jtest=1,2,3,4,5,6 -> (p,q) = (1,1), (2,1), (1,2), (2,2), (1,3), (2,3)
!  p indicates component which is non-zero, q-1 indicates order of the quadratic
!
      integer, parameter :: P1Q1=1, P2Q1=2, P1Q2=3, P2Q2=4, P1Q3=5, P2Q3=6, G=7
!
      select case (jtest)
!
      case (P1Q1)
        U0test(:,1)=wamp; U0test(:,2)=0.; U0test(:,3)=0.
        gU0test=0.
!
      case (P2Q1)
        U0test(:,1)=0.; U0test(:,2)=wamp; U0test(:,3)=0.
        gU0test=0.
!
      case (P1Q2)
!        U0test(:,1)=wamp*(1.+z(n)-zc(n)); U0test(:,2)=0.; U0test(:,3)=0.
!       gU0test(:,1)=wamp-...; gU0test(:,2)=0.; gU0test(:,3)=0.
!
      case (P2Q2)
!        U0test(:,1)=0.; U0test(:,2)=wamp*(1.+z(n)-zc(n)); U0test(:,3)=0.
!        gU0test(:,1)=0.; gU0test(:,2)=wamp-...; gU0test(:,3)=0.
!
      case (P1Q3)
!        U0test(:,1)=wamp*(1.-z(n)+zc(n)); U0test(:,2)=0.; U0test(:,3)=0.
!        gU0test(:,1)=-wamp+...; gU0test(:,2)=0.; gU0test(:,3)=0.
!
      case (P2Q3)
!        U0test(:,1)=0.; U0test(:,2)=wamp*(1.-z(n)+zc(n)); U0test(:,3)=0.
!        gU0test(:,1)=0.; gU0test(:,2)=-wamp+...; gU0test(:,3)=0.
!
      case default
        U0test=0.;gU0test=0.
!
      endselect
!
    endsubroutine set_U0test_quasiperiodic
!***********************************************************************
    subroutine set_U0test_onlyconstant(U0test,gU0test,jtest)
!
      use Cdata
      use Sub
!
      real, dimension (nx,3), intent(out) :: U0test,gU0test
      integer,                intent(in)  :: jtest
!
!  set U0test and gU0test for each of the testflow indices jtest (which stands for p and q)
!  order convention: jtest=1,2,3,4,5,6 -> (p,q) = (1,1), (2,1), (1,2), (2,2), (1,3), (2,3)
!  p indicates component which is non-zero, q-1 indicates degree of the polynom
!
      integer, parameter :: P1Q1=1, P2Q1=2, P1Q2=3, P2Q2=4, P1Q3=5, P2Q3=6, G=7

      select case (jtest)
!
      case (P1Q1)
        U0test(:,1)=wamp; U0test(:,2)=0.; U0test(:,3)=0.
        gU0test=0.
!
      case (P2Q1)
        U0test(:,1)=0.; U0test(:,2)=wamp; U0test(:,3)=0.
        gU0test=0.

      case default
        U0test=0.;gU0test=0.
!
      end select
!
    endsubroutine set_U0test_onlyconstant
!***********************************************************************
    subroutine set_U0test_onlylinear(U0test,gU0test,jtest)
!
      use Cdata
      use Sub
!
      real, dimension (nx,3), intent(out) :: U0test,gU0test
      integer,                intent(in)  :: jtest
!
!  set U0test and gU0test for each of the testflow indices jtest (which stands for p and q)
!  order convention: jtest=1,2,3,4,5,6 -> (p,q) = (1,1), (2,1), (1,2), (2,2), (1,3), (2,3)
!  p indicates component which is non-zero, q-1 indicates degree of the polynom
!
      integer, parameter :: P1Q1=1, P2Q1=2, P1Q2=3, P2Q2=4, P1Q3=5, P2Q3=6, G=7

      select case (jtest)
!
      case (P1Q2)
        U0test(:,1)=wamp*z(n); U0test(:,2)=0.; U0test(:,3)=0.
        gU0test(:,1)=wamp; gU0test(:,2)=0.; gU0test(:,3)=0.
!
      case (P2Q2)
        U0test(:,1)=0.; U0test(:,2)=wamp*z(n); U0test(:,3)=0.
        gU0test(:,1)=0.; gU0test(:,2)=wamp; gU0test(:,3)=0.

      case default
        U0test=0.;gU0test=0.
!
      end select
!
    endsubroutine set_U0test_onlylinear
!***********************************************************************
    subroutine set_U0test_quadratic(U0test,gU0test,jtest)
!
!  set testflow
!
!  14-sep-09/MR: adapted from set_U0test_W11_W22
!
      use Cdata
      use Sub
!
      real, dimension (nx,3), intent(out) :: U0test,gU0test
      integer,                intent(in)  :: jtest
!
!  set U0test and gU0test for each of the testflow indices jtest (which stands for p and q)
!  order convention: jtest=1,2,3,4,5,6 -> (p,q) = (1,1), (2,1), (1,2), (2,2), (1,3), (2,3)
!  p indicates component which is non-zero, q-1 indicates degree of the polynom
!
      integer, parameter :: P1Q1=1, P2Q1=2, P1Q2=3, P2Q2=4, P1Q3=5, P2Q3=6, G=7
!
      select case (jtest)
!
      case (P1Q1)
        U0test(:,1)=wamp; U0test(:,2)=0.; U0test(:,3)=0.
        gU0test=0.
!
      case (P2Q1)
        U0test(:,1)=0.; U0test(:,2)=wamp; U0test(:,3)=0.
        gU0test=0.
!
      case (P1Q2)
        U0test(:,1)=wamp*z(n); U0test(:,2)=0.; U0test(:,3)=0.
        gU0test(:,1)=wamp; gU0test(:,2)=0.; gU0test(:,3)=0.
!
      case (P2Q2)
        U0test(:,1)=0.; U0test(:,2)=wamp*z(n); U0test(:,3)=0.
        gU0test(:,1)=0.; gU0test(:,2)=wamp; gU0test(:,3)=0.
!
      case (P1Q3)
        U0test(:,1)=wamp*zq2(n); U0test(:,2)=0.; U0test(:,3)=0.
        gU0test(:,1)=wamp*z(n); gU0test(:,2)=0.; gU0test(:,3)=0.
!
      case (P2Q3)
        U0test(:,1)=0.; U0test(:,2)=wamp*zq2(n); U0test(:,3)=0.
        gU0test(:,1)=0.; gU0test(:,2)=wamp*z(n); gU0test(:,3)=0.
!
      case default
        U0test=0.;gU0test=0.
!
      endselect
!
    endsubroutine set_U0test_quadratic
!***********************************************************************
    subroutine rprint_testflow(lreset,lwrite)                   ! -> Standardschnittstelle
!
!  reads and registers print parameters relevant for testflow fields
!
!   3-jun-05/axel: adapted from rprint_magnetic
!
      use Cdata
      use Diagnostics
!
      logical           :: lreset
      logical, optional :: lwrite
!
      integer           :: iname,i,j,p,ifound,ifoundold
      character         :: cind
      character(len=2)  :: cind2
      character(len=20) :: name
      logical           :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
!
        idiag_Fipq=0; idiag_Qpq=0
!
        idiag_gal=0; idiag_galij=0
        idiag_aklam=0; idiag_aklamij=0
        idiag_nu=0 ; idiag_nuij=0
!
        idiag_ux0mz=0; idiag_uy0mz=0; idiag_uz0mz=0;
        idiag_upqrms=0; idiag_hpqrms=0
!
      endif
!
!  check for those quantities that we want to evaluate online
!
      ifoundold=-1
!
      do iname=1,nname
!

        ifound = fparse_name(iname,cname(iname),cform(iname),'gal',idiag_gal)
        ifound = ifound + fparse_name(iname,cname(iname),cform(iname),'aklam',idiag_aklam)
        ifound = ifound + fparse_name(iname,cname(iname),cform(iname),'nu',idiag_nu)
        ifound = ifound + fparse_name(iname,cname(iname),cform(iname),'kappa',idiag_kappa)
        ifound = ifound + fparse_name(iname,cname(iname),cform(iname),'xi',idiag_chi)
!
        do i=1,2
        do j=1,2
!
          cind2 = gen_2ind(i,j)

          ifound = ifound + fparse_name(iname,cname(iname),cform(iname),'gal'//cind2,idiag_galij(i,j))
          ifound = ifound + fparse_name(iname,cname(iname),cform(iname),'aklam'//cind2,idiag_aklamij(i,j))
          ifound = ifound + fparse_name(iname,cname(iname),cform(iname),'nu'//cind2,idiag_nuij(i,j))
!
        enddo
        enddo

        p=1
        do j=0,njtestflow
!
          if (j.eq.0) then
            name = '0rms'
          else
!
            cind2 = gen_2ind(p,floor((j+1)/2.))
            name = cind2//'rms'
!
            p=3-p

          endif
!
          ifound = ifound + fparse_name(iname,cname(iname),cform(iname),'u'//name,idiag_upqrms(j))
          ifound = ifound + fparse_name(iname,cname(iname),cform(iname),'h'//name,idiag_hpqrms(j))

          if ( j.gt.0 ) then
!
            do i=1,3
              ifound = ifound + fparse_name(iname,cname(iname),cform(iname),'F'//cind2,idiag_Fipq(i,j))
            enddo
!
            write(cind,'(i1)') j
            ifound = ifound + fparse_name(iname,cname(iname),cform(iname),'Q'//cind,idiag_Qpq(j))

          endif

        enddo
!
        idiag_map(iname) = iname       ! initialising the index mapping vector for cname und cform

        if ( ifound==0 .and. ifoundold>0 ) &
          print*, 'testflow_z, rprint_testflow: Warning - diagnostic ouput for line',iname, &
                  ' and beyond will be messed up!!!'
!
        if ( ifound/=0 ) &
          ifoundold=ifound
!
      enddo

      nname_old = nname
!
      call update_diag( idiag_gal, idiag_galij, 'gal' )
      call update_diag( idiag_aklam, idiag_aklamij, 'aklam' )
      call update_diag( idiag_nu, idiag_nuij, 'nu' )

      if ( idiag_gal/=-1 )   call correct_inds(idiag_galij)
      if ( idiag_aklam/=-1 ) call correct_inds(idiag_aklamij)
      if ( idiag_nu/=-1 )    call correct_inds(idiag_nuij)
!
      if ( idiag_nu /= 0 ) &
        idiag_aklam = -1
!
      if ( idiag_aklam /= 0 ) &
        idiag_gal = -1
!
!  check for those quantities for which we want xy-averages
!
      do iname=1,nnamez
!
        call parse_name(iname,cnamez(iname),cformz(iname),'ux0mz',idiag_ux0mz)
        call parse_name(iname,cnamez(iname),cformz(iname),'uy0mz',idiag_uy0mz)
        call parse_name(iname,cnamez(iname),cformz(iname),'uz0mz',idiag_uz0mz)
!
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
!
        do i=1,2
        do j=1,2
!
          cind2 = gen_2ind(i,j)
          write(3,*) 'idiag_gal'//cind2//'=',idiag_galij(i,j)
          write(3,*) 'idiag_aklam'//cind2//'=',idiag_aklamij(i,j)
          write(3,*) 'idiag_nu'//cind2//'=',idiag_nuij(i,j)
!
        enddo
        enddo

        do j=1,njtestflow
!
          do i=1,3

            cind2 = gen_2ind(i,j)
            write(3,*) 'idiag_Fipq'//cind2//'=',idiag_Fipq(i,j)
!
          enddo
!
          write(cind,'(i1)') j
          write(3,*) 'idiag_Qpq'//cind//'=',idiag_Qpq(j)
!
        enddo
!
        p=1
        do j=0,njtestflow
!
          if (j.eq.0) then
            name = '0rms='
          else
!
            name = gen_2ind(p,floor((j+1)/2.))//'rms='
            p=3-p
!
          endif
!
          write(3,*) 'idiag_u'//name,idiag_upqrms(j)
          write(3,*) 'idiag_h'//name,idiag_hpqrms(j)
!
        enddo
!
        write(3,*) 'idiag_ux0mz=',idiag_ux0mz
        write(3,*) 'idiag_uy0mz=',idiag_uy0mz
        write(3,*) 'idiag_uz0mz=',idiag_uz0mz
!
        write(3,*) 'iuutest=',iuutest
!
        write(3,*) 'ntestflow=',ntestflow
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'nnamexy=',nnamexy
      endif
!
    endsubroutine rprint_testflow
!***********************************************************************
    character(len=2) function gen_2ind(i,j)
!
      integer :: i,j
      intent(in) :: i,j
!
      character(len=20) c1, c2
      write(c1,fmt='(i19)') i
      write(c2,fmt='(i19)') j

      gen_2ind = trim(adjustl(c1))//trim(adjustl(c2))

    endfunction gen_2ind
!***********************************************************************
    subroutine mark_del_elems( carray, leng, indices )

      integer, dimension(2,2) :: indices
      integer                 :: leng
!
      character*(*), dimension(*) :: carray
!
      intent(inout) :: indices, carray, leng
!
      integer :: i,j,ind
!
      do i=1,2
      do j=1,2
        if ( indices(i,j) /= 0 ) then
!
          ind = indices(i,j)

          call del_elem( carray, idiag_map(ind), leng )
!
          idiag_map(ind) = 0
          if ( ind.lt.nname_old ) &
            idiag_map(ind+1:nname_old) = idiag_map(ind+1:nname_old)-1

        endif
      enddo
      enddo
!
    endsubroutine mark_del_elems
!***********************************************************************
    subroutine del_elem( carray, index, leng )

      integer, intent(in)          :: index
      integer, intent(inout)       :: leng
      character*(*), dimension(*), intent(inout) :: carray
!
      if ( index>0.and.index<=leng ) then
        carray(index:leng-1) = carray(index+1:leng)
        leng = leng-1
      endif
!
    endsubroutine del_elem
!***********************************************************************
    subroutine insert_array( carray, cinsert, leng_insert, index, leng )

      integer                     :: index, leng, leng_insert, i
      character*(*), dimension(*) :: carray, cinsert
!
      intent(in)    :: index, cinsert, leng_insert
      intent(inout) :: leng, carray
!
      if ( index>0.and.index<=leng+1 ) then
!
        do i=leng,index,-1
          carray(i+leng_insert) = carray(i)
        enddo
!
        carray(index:index+leng_insert-1) = cinsert(1:leng_insert)
!
        leng = leng+leng_insert

      endif
!
    endsubroutine insert_array
!***********************************************************************
    subroutine insert_array_mult( carray, cinsert, mult, index, leng )

      integer                     :: index, leng, mult, i
      character*(*), dimension(*) :: carray
      character*(*)               :: cinsert
!
      intent(in)    :: index, cinsert, mult
      intent(inout) :: leng, carray
!
      if ( index>0.and.index<=leng+1 ) then
!
        do i=leng,index,-1
          carray(i+mult) = carray(i)
        enddo
!
        carray(index:index+mult-1) = cinsert
!
        leng = leng+mult
!
      endif
    endsubroutine insert_array_mult
!***********************************************************************
    subroutine update_diag( mat_diag_ind, mat_diag_indices, cdiagname )
!
      use Cdata
      use Diagnostics
!
      integer                 :: mat_diag_ind
      integer, dimension(2,2) :: mat_diag_indices
!
      intent(inout) :: mat_diag_ind, mat_diag_indices
!
      character*(*), intent(in) :: cdiagname
!
      character(len=30) cformat
!
      integer :: iname, i, j, nname_form, indx
!
      if ( mat_diag_ind.eq.0 ) then                           ! if complete tensor was not already inquired
!
        do iname=1,nname
        do i=1,2
        do j=1,2
!
          if ( mat_diag_indices(i,j)/=0 ) then                 ! if any of the elements is inquired
!
            mat_diag_ind = -2                                  ! flag the whole tensor to be calculated
            return
!
          endif
!
        enddo
        enddo
        enddo
!
      else
!
        call mark_del_elems( cname, nname, mat_diag_indices )     ! delete inquiries for individual elements and mark them in idiag_map by 0
!
        indx = idiag_map(mat_diag_ind)
        cformat = cname(indx)
        cformat = cformat(index(cformat,'('):)

        nname_form=nname

        call del_elem( cname, indx, nname )                                                         ! replace inquiry for the whole tensor by inquiries for
        call insert( cname, (/cdiagname//'11'//cformat,cdiagname//'12'//cformat, &
                              cdiagname//'21'//cformat,cdiagname//'22'//cformat/), 4, indx, nname ) ! all elements
        do i=1,2
          do j=1,2
!
            mat_diag_indices(i,j) = indx
            indx = indx+1
!
          enddo
        enddo

        if ( mat_diag_ind.lt.nname_old ) &
          idiag_map(mat_diag_ind+1:nname_old) = idiag_map(mat_diag_ind+1:nname_old)+3
!
        cformat = cform(indx)
        call insert( cform, cformat, 3, indx+1, nname_form )
!
        mat_diag_ind = -1
!
      endif
!
    endsubroutine update_diag
!***********************************************************************
    subroutine correct_inds(mat_diag_indices)
!
      integer, dimension(2,2) :: mat_diag_indices
!
      integer i,j,ind
!
      do i=1,2
        do j=1,2
!
          ind = mat_diag_indices(i,j)
!
          if (ind.gt.0) &
            mat_diag_indices(i,j) = idiag_map(ind)
!
        enddo
      enddo
!
    endsubroutine correct_inds
!***********************************************************************
      SUBROUTINE ISORTP(N,RA,IP)
!
      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(IN)  :: RA(*)
      INTEGER, INTENT(OUT) :: IP(*)

      INTEGER RRA,L,IR,I,J,IPA,J2,K
!
      IP(1) = 1
      IF (N.LE.1) GOTO 1001
      L = N/2+1
      IR = N
!
      DO I=2,N
         IP(I) = I
      ENDDO
!
   21 CONTINUE
         IF (IR.GT.1) THEN
            IF (L.GT.1) THEN
               L = L-1
               J = L
            ELSE
               IPA = IP(1)
               IP(1) = IP(IR)
               IP(IR) = IPA
               IR = IR-1
               J = 1
            END IF
            IPA = IP(J)
            RRA = RA(IPA)
   22       CONTINUE
               J2 = J+J
               IF (J2.LE.IR) THEN
                  K = J2
                  IF (K.LT.IR) THEN
                     IF ( RA(IP(K)).LT.RA(IP(K+1)) ) K = K+1
                  END IF
                  IF ( RRA.LT.RA(IP(K)) ) THEN
                     IP(J) = IP(K)
                     J = K
                  ELSE
                     GOTO 23
                  END IF
                  GOTO 22
               END IF
   23          IP(J) = IPA
            GOTO 21
         END IF

 1001 RETURN
 endsubroutine isortp
!***********************************************************************
endmodule Testflow
