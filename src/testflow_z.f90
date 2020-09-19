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
!  integer, parameter :: njtestflow=8
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
  use Diagnostics              !!, only : idiag_map
!
  implicit none
!
  include 'testflow.h'
!
  interface insert                 ! Overload the 'insert' function
    module procedure insert_array
    module procedure insert_array_mult
  endinterface
!
  interface gen_ind                ! Overload the 'gen_ind' function
    module procedure gen_1ind
    module procedure gen_2ind
  endinterface
!
  interface mark_del_elems         ! Overload the 'mark_del_elems' function
    module procedure mark_del_elems1
    module procedure mark_del_elems2
  endinterface

  interface correct_inds           ! Overload the 'correct_inds' function
    module procedure correct_0inds
    module procedure correct_1inds
    module procedure correct_2inds
  endinterface


  interface update_diag            ! Overload the 'update_diag' function
    module procedure update_diag1
    module procedure update_diag2
  endinterface

  !!integer, parameter :: njtestflow=6
  integer :: njtestflow_loc
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
  real :: cs2test=1., cs2test1=1.
!
  logical :: ltestflow_upw_lnrho=.false., ltestflow_upw_uu=.false.
!
  integer :: nuuinit

  integer, dimension(2) :: izrange
  real, dimension(2) :: zrange, ca, cl, cq, cqa, cqq
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
  real :: nutest=0.,     &
          nutest1=0.,    &
          wamp=1.,       &
          ktestflow=1.,  &
          ktestflow1=1., &
          valid_zrange=1.
!
  logical :: reinitialize_uutest=.false., &
             lsoca_testflow=.true., &
             lkinem_testflow=.false., &
             lburgers_testflow=.false., &
             lprescribed_velocity=.false., &
             lremove_mean_momenta_testflow=.false., &
             lremove_mean_enthalpy_z=.false., &
             lshear_as_param=.false.
!
  character (len=labellen) :: itestflow='W11-W22'
!
  namelist /testflow_run_pars/ reinitialize_uutest, &    ! flag for reinitialization of testflow solution
                               zextent,             &    ! ??
                               lsoca_testflow,      &    ! flag for SOCA
                               lkinem_testflow,     &    ! flag for kinematic calculation
                               lburgers_testflow,   &    ! flag for disconnecting enthalpy(pressure) from velocity
                               lprescribed_velocity,&    ! flag for prescribed velocity, prescription via p%fcont, only effective if lkinem_testflow=.true.
                               lremove_mean_momenta_testflow, &  ! flag for removing mean momenta in 0-solution
                               lremove_mean_enthalpy_z,& ! flag for removing xy-mean of enthalpy in 0-solution
                               lshear_as_param,     &    ! flag for treating of shear as a problem parameter, instead of a mean flow
                               nutest,              &    ! viscosity in testflow equations
                               nutest1,             &    ! reciprocal viscosity
                               itestflow,           &    ! name of used testflow set, legal values
                                                         ! 'W11-W22', 'quadratic', 'quadratic+G', 'none'
                               ktestflow,           &    ! wavenumber for testflow (if sinusoidal)
                               lugu_as_aux,         &    ! ??
                               duuinit,             &    !
                               wamp,                &    ! amplitude of the testflows
                               valid_zrange              ! relative size of z-range from which coefficients are calculated
!
  ! diagnostic variables (needs to be consistent with reset list below)
!
  integer :: idiag_gal=0                            ! DIAG_DOC: GAL-coefficients,     couple  $\overline F$ and $\overline U$
  integer :: idiag_aklam=0                          ! DIAG_DOC: AKA-$\lambda$-tensor, couples $\overline F$ and $\overline W = \nabla\times{\overline U}$
  integer :: idiag_gamma=0                          ! DIAG_DOC: $\gamma$-vector,      couples $\overline F$ and $\nabla\cdot{\overline U}$
  integer :: idiag_nu=0                             ! DIAG_DOC: $\nu$-tensor,         couples $\overline F$ and $\partial^2 {\overline U}/\partial z^2$
  integer :: idiag_zeta=0                           ! DIAG_DOC: $\zeta$-vector,       couples $\overline F$ and ${\overline G}_z = \nabla_z {\overline H}$
  integer :: idiag_xi=0                             ! DIAG_DOC: $\xi$-vector,         couples $\overline F$ and $\partial^2 {\overline H}/\partial z^2$

  integer :: idiag_aklamQ=0                         ! DIAG_DOC: $aklam^Q$-vector,     couples $\overline Q$ and $\overline W$
  integer :: idiag_gammaQ=0                         ! DIAG_DOC: $\gamma^Q$-scalar,    couples $\overline Q$ and $\nabla\cdot{\overline U}=dU_z/dz$
  integer :: idiag_nuQ=0                            ! DIAG_DOC: $\nu^Q$-vector,       couples $\overline Q$ and $\partial^2 \overline U/\partial z^2$
  integer :: idiag_zetaQ=0                          ! DIAG_DOC: $\zeta^Q$-scalar,      couples $\overline Q$ and ${\overline G}_z$
  integer :: idiag_xiQ=0                            ! DIAG_DOC: $\xi^Q$-scalar,        couples $\overline Q$ and $\partial^2 {\overline H}/\partial z^2$
!
  integer, dimension(3,3) :: idiag_galij=0          ! DIAG_DOC:
  integer, dimension(3,2) :: idiag_aklamij=0        ! DIAG_DOC: $\alpha_{K,ij}$
  integer, dimension(3)   :: idiag_gammai=0         ! DIAG_DOC: $\gamma_i$
  integer, dimension(3,3) :: idiag_nuij=0           ! DIAG_DOC: $\nu_{ij}$
!
  integer, dimension(3) :: idiag_zetai=0, &         ! DIAG_DOC: $\zeta_i$
                           idiag_xii=0, &           ! DIAG_DOC: $\xi_i$
                           idiag_nuQi=0             ! DIAG_DOC: $\nu^Q_i$
  integer, dimension(2) :: idiag_aklamQi=0          ! DIAG_DOC: $aklam^Q_i$
!
  integer, dimension(3,njtestflow) :: idiag_Fipq=0      ! DIAG_DOC: ${\cal F}_i^{pq}$
  integer, dimension(njtestflow)   :: idiag_Qpq=0       ! DIAG_DOC: ${\cal Q}^{pq}$
  integer, dimension(0:njtestflow) :: idiag_upqrms=0    ! DIAG_DOC: $\left<{u^{pq}}^2\right>$
  integer, dimension(0:njtestflow) :: idiag_hpqrms=0    ! DIAG_DOC: $\left<{h^{pq}}^2\right>$
!
  integer :: idiag_ux0mz=0                              ! DIAG_DOC: $\left<u_{x}\right>_{xy}$
  integer :: idiag_uy0mz=0                              ! DIAG_DOC: $\left<u_{y}\right>_{xy}$
  integer :: idiag_uz0mz=0                              ! DIAG_DOC: $\left<u_{z}\right>_{xy}$
!
  integer :: nname_old
!
  integer :: idiag_map(100)
!
  contains
!
!***********************************************************************
    subroutine register_testflow()          ! -> Default interface
!
!  Initialise variables which should know that we solve for the vector
!  potential: iuutest, etc; increase nvar accordingly
!
!   3-jun-05/axel: adapted from register_magnetic
!
      use Cdata
      use Mpicomm, only: stop_it
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
    subroutine initialize_testflow(f)           ! -> Default interface
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
      use sub, only:zlocation
      use EquationOfState, only : cs0
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: range, cn

      integer :: izpos
      logical :: lproc
!
      if (itestflow=='W11-W22') then
        njtestflow_loc=njtestflow            ! remove unnecessary memory consumption, if possible
      else
        njtestflow_loc=njtestflow
      endif
!
! Set nutest from the value uf nu in Hydro, if not specified here
!
      if ( nutest==0. .and. nutest1==0. ) &
        call getnu(nu_input=nutest)
!
!  Precalculate nutest if 1/nutest (==nutest1) is given instead
!
      if (nutest1/=0.) then
        nutest=1./nutest1
      else if ( nutest==0. ) then
        !!call stop_it('duutest_dt: Error - nutest==0!!!')
      endif
!
      if ( cs0 /= 0. ) then
        cs2test = cs0*cs0; cs2test1=1./cs2test
      endif
!
!  set to zero and then rescale the testflow
!  (in future, could call something like init_uu_simple)
!
     if ( itestflow=='quadratic' ) then
       zq2 = 0.5*z*z
     elseif ( itestflow=='W11-W22' ) then
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
! ensure that module forcing adds (discontinuous) force in f
!
     ltestflow_forcing = .true.
!
     if (lkinem_testflow) &
       lhydro_forcing = .false. ! ???

! ensure that lprescribed_velocity finds continuous forcing

     if ( lprescribed_velocity .and. .not.lforcing_cont ) &
       lprescribed_velocity=.false.

! ensure that lprescribed_velocity finds lkinem_testflow=.true.

     if (lprescribed_velocity) &
       lkinem_testflow=.true.

     if (lburgers_testflow) &
       lremove_mean_enthalpy_z=.false.
!
! calculate z-index bounds from valid_zrange for actual processor
!
     izrange(1)=n1
     izrange(2)=n2

     if ( itestflow=='quasi-periodic' ) then

       valid_zrange = max(0.3,valid_zrange)
       valid_zrange = min(1.0,valid_zrange)

       if ( valid_zrange<1. ) then

         range = Lz*valid_zrange

         zrange(1) = xyz0(3) + 0.5*(Lz-range)        ! lower and upper z-limit of valid range
         zrange(2) = zrange(1) + range

         izrange(1) = n1-1
         izrange(2) = n2+1

         if ( zrange(1)<=z(n1) ) then
           zrange(1) = z(n1)-.01
         else

           call zlocation(zrange(1),izrange(1),lproc)
           if (.not.lproc) izrange(1)=n2+1

           cn = (xyz0(3)-zrange(1))**2
           ca(1) = -xyz0(3)*zrange(1)**2/cn           ! coefficients of periodized linear function
           cl(1) = (xyz0(3)**2+zrange(1)**2)/cn
           cq(1) = -xyz0(3)/cn

           cqa(1) =  xyz0(3)*zrange(1)
           cqq(1) = -zrange(1)/(xyz0(3)-zrange(1))      ! coefficients of periodized quadratic function

         endif
!
         if ( zrange(2)>=z(n2) ) then
           zrange(2) = z(n2)+.01
         else

           call zlocation(zrange(2),izrange(2),lproc)
           if (.not.lproc) izrange(2)=n1-1

           cn = (xyz1(3)-zrange(2))**2
           ca(2) = -xyz1(3)*zrange(2)**2/cn           ! coefficients of periodized linear function
           cl(2) = (xyz1(3)**2+zrange(2)**2)/cn
           cq(2) = -xyz1(3)/cn

           cqa(2) =  xyz1(3)*zrange(2)                ! coefficients of periodized quadratic function
           cqq(2) = -zrange(2)/(xyz1(3)-zrange(2))

         endif

 !!        print*, 'zrange, ca, cl, cq:', xyz0(3), xyz1(3), zrange, ca, cl, cq
        else
         zrange(1) = z(n1)-.01
         zrange(2) = z(n2)+.01
       endif
     endif

    if ( itestflow=='W11-W22' .and. idiag_gal/=0 ) then
      idiag_gal=0
      idiag_galij=0
    endif

    if ( idiag_aklam/=0 .and. njtestflow<4 ) then

      print*,'Warning: Number of testflows njtestflow=',njtestflow,' is too small for AKA-Lambda tensor to be calculated'
      idiag_aklam=0

    endif

    if ( idiag_nu/=0 .and. ( itestflow=='W11-W22' .and. njtestflow<4 .or. itestflow=='quadratic' .and. njtestflow<6) ) then

      print*,'Warning: Number of testflows njtestflow=',njtestflow,' is too small for viscosity tensor to be calculated'
      idiag_nu=0

    endif

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
        if (iugutest==0) &
          call farray_register_auxiliary('ugu',iugutest,vector=3*njtestflow)
        if (iugutest/=0) then
          if (lroot) print*, 'initialize_testflow: iugu = ', iugutest
          call farray_index_append('iugutest',iugutest)
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
    subroutine init_uutest(f)          ! -> Default interface
!
!  initialise testflow-solutions; called from start.f90
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use Mpicomm, only:stop_it
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
      do j=1,ninit*njtestflow
!
!  loop over different choices of initial conditions
!
        if ( inituutest(j)(1:15)=='gaussian-noise-' ) then
!
          read( inituutest(j)(16:16), '(i1)' ) jtest
          selector='gaussian-noise'
!
        else if ( inituutest(j)(1:10)=='sinwave-x-' ) then
!
          read( inituutest(j)(11:11), '(i1)' ) jtest
          selector='sinwave-x'
        else if ( inituutest(j)(1:9)=='constant-' ) then
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
    subroutine pencil_criteria_testflow()          ! -> Default interface
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
    subroutine pencil_interdep_testflow(lpencil_in)         ! -> Default interface
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
    subroutine read_testflow_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testflow_init_pars, IOSTAT=iostat)
!
    endsubroutine read_testflow_init_pars
!***********************************************************************
    subroutine write_testflow_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testflow_init_pars)
!
    endsubroutine write_testflow_init_pars
!***********************************************************************
    subroutine read_testflow_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testflow_run_pars, IOSTAT=iostat)
!
    endsubroutine read_testflow_run_pars
!***********************************************************************
    subroutine write_testflow_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testflow_run_pars)
!
!      ! remove unnecessary memory consumption
!      if (itestflow=='W11-W22') njtestflow=4
!
    endsubroutine write_testflow_run_pars
!***********************************************************************
    subroutine duutest_dt(f,df,p)                    ! -> Default interface
!
!  testflow evolution:
!
!  calculate du^(pq)/dt = -( u^(0) .grad u^(pq) + u^(pq).grad u ) + < u^(0) .grad u^(pq) + u^(pq).grad u >
!        [alternatively:  -( u^(pq).grad u^(0)  + u.grad u^(pq) ) + < u^(pq).grad u^(0)  + u.grad u^(pq) >]
!                                                 - grad h^(pq)
!                         -U^(pq).grad u - u.grad U^(pq)
!                         + viscous terms
!
!  and dh^(pq)/dt = -( u^(0).grad h^(pq) + u^(pq).grad h ) + < u^(0).grad h^(pq) + u^(pq).grad h >
!  [alternatively:  -( u^(pq).grad h^(0) + u.grad h^(pq) ) + < u^(pq).grad h^(0) + u.grad h^(pq) > ]
!                                       - cs2*div u^(pq)
!                   -U.grad h -u.grad H
!
!
!  12-mar-08/axel: coded
!  24-jun-08/MR: modified
!  15-feb-13/MR: handling of test cases with mean enthalpy added
!
      use Cdata
      use Sub
      use Mpicomm, only: stop_it
      use density, only: glnrhomz
      use Hydro, only: uumz, guumz, coriolis_cartesian, ampl_fcont_uu
      use Forcing, only: forcing_cont
      use Shear, only: shear_variables
      use Deriv, only: der5_single
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
      real, dimension (nx,3) :: U0test,gU0test                      !!!MR: needless x-dimension for 0-quantities
      real, dimension (nx)   :: hhtest,U0ghtest,gH0test,upq2,divutest,divufluct,diffus_nu,advec_uu
!
      logical :: ltestflow_out
!
      integer :: jtest,j,i,i3,i4, nl
      logical, save :: lfirst_call=.true.
      integer :: iuxtest, iuytest, iuztest, ihhtest
!
      character :: cjtest
      character (len=fnlen) :: file

      !if (ldiagnos .and. n==4.and.m==4) then
      !print*, 'uumz(1):', uumz(:,1)
      !print*, 'guumz(1):', guumz(:,1)
      !endif
!
      if ( iuutest==0 ) return

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
        uufluct = f(l1:l2,m,n,iuxtest:iuztest)           ! ufluct = u0
        if (.not.lburgers_testflow) &
          call grad(f,ihhtest,ghfluct)                   ! ghfluct = grad(h0)
!
      else
!
!  calculate uufluct=U-Umean out of the main run
!
        do j=1,3
          uufluct(:,j)=p%uu(:,j)-uumz(n,j)
        enddo
!
! calculate ghfluct=gh-ghmean
!
        ghfluct=p%glnrho
        ghfluct(:,3)=ghfluct(:,3)-glnrhomz(n-n1+1)

      endif
!
!  do each of the njtestflow test flows at a time
!  jtest=0  refers to primary turbulence (u^(0), h^(0))
!
      advec_uu=0.
      do jtest=0,njtestflow_loc
!
!  identify boundary conditions
!
        if (headtt) then
!
          write(cjtest,'(i1)') jtest
          call identify_bcs('uxtest-'//cjtest,iuxtest)
          call identify_bcs('uytest-'//cjtest,iuytest)
          call identify_bcs('uztest-'//cjtest,iuztest)
          if (.not.lburgers_testflow) &
            call identify_bcs('hhtest-'//cjtest,ihhtest)
!
        endif
!
!  velocity vector and enthalpy
!
        if ( jtest>0 .or. .not.lprescribed_velocity ) then

          uutest=f(l1:l2,m,n,iuxtest:iuztest)! testflow solution # jtest
          if (.not.lburgers_testflow) hhtest=f(l1:l2,m,n,ihhtest)
          !!if ( jtest==1 ) print*, 'uutest, hhtest:', minval(uutest),maxval(uutest), minval(hhtest),maxval(hhtest)
!
!  div u term
!
          call div(f,iuxtest,divutest)
!
!  gradient of (pseudo) enthalpy
!
          if (.not.lburgers_testflow) call grad(f,ihhtest,ghtest)
!
        endif

        if (jtest>0) then
!
          select case (itestflow)                       ! get testflow U^(pq), grad U^(pq), grad H^(pq)
!
            case ('W11-W22')
              call set_U0test_W11_W22(U0test,gU0test,gH0test,jtest)
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

              !!if ( (jtest==3 .or. jtest==4) .and. m==5 ) &
                 !!print*, jtest,z(n), U0test(1,jtest-2),gU0test(1,jtest-2)

            case ('quadratic+G')
!
              if ( jtest<=6 ) then
                call set_U0test_quadratic(U0test,gU0test,jtest)
                gH0test=0.
              else
                U0test=0.; gU0test=0.; gH0test=wamp
              endif
!
            case ('none')
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
        if ( .not.lburgers_testflow .and. (jtest>0 .or. .not.lprescribed_velocity) ) &
          df(l1:l2,m,n,ihhtest)=df(l1:l2,m,n,ihhtest) - cs2test*divutest
!
!       testflow inhomogeneity
!
        if (jtest>0 .and. .not.lburgers_testflow) then
!
          if (lkinem_testflow) then
            !!call dot_mn(U0test,ghfluct,U0ghtest)          ! MR: checked by numbers
            call u_dot_grad(f,iuutest+3, ghfluct, U0test, U0ghtest,UPWIND=ltestflow_upw_lnrho)   ! Utest.grad(h0)
          else
            call u_dot_grad(f,ilnrho, ghfluct, U0test, U0ghtest,UPWIND=ltestflow_upw_lnrho)      ! Utest.grad(h)
            if (ltestflow_upw_lnrho) &
              U0ghtest = U0ghtest - abs(U0test(:,3))* &     ! in slot ilnrho there is hh, not hhfluct
                         der5_single(glnrhomz,n-n1+1,dz_1(n1:n2))/(60.*dz_1(n-n1+1)**5)
          endif
!
          if ( .not.lburgers_testflow) &
            df(l1:l2,m,n,ihhtest)=df(l1:l2,m,n,ihhtest)-uufluct(:,3)*gH0test-U0ghtest
!
        endif
!
!  rhs of momentum equation (nonlinear terms already in df, if needed!),
!  du^pq/dt = n.l.Terms - grad h^pq -u.gradU^pq - U^pq.gradu
!
        if ( .not.lburgers_testflow .and. (jtest>0 .or. .not.lprescribed_velocity) ) &
          df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) - ghtest
!
!       testflow inhomogeneity
!
        if ( jtest>0 ) then
!
          if (lkinem_testflow) then
!
            call gij(f,iuutest,uijfluct,1)
!
          else
!
            uijfluct = p%uij
!
            do i=1,3                                        ! determine fluctuating part of uij
              uijfluct(:,i,3) = uijfluct(:,i,3) - guumz(n,i)
            enddo
!
          endif
!
          call h_dot_grad(U0test,uijfluct,uufluct,U0testgu) ! no upwinding?
                                          !this parameter without meaning in cartesian geometry
!
          call multsv(uufluct(:,3),gU0test,ugU0test)
!
          df(l1:l2,m,n,iuxtest:iuztest) = df(l1:l2,m,n,iuxtest:iuztest) - ugU0test - U0testgu
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

              do i=1,2                                                       ! determine fluctuating part
                uijfluct(:,i,i) = uijfluct(:,i,i) + guumz(n,3)/3.            ! ~
              enddo
!
              uijfluct(:,3,3) = uijfluct(:,3,3) - 2.*guumz(n,3)/3.           ! ~
              uijfluct(:,1,3) = uijfluct(:,1,3) - 0.5*guumz(n,1)             ! ~
              uijfluct(:,3,1) = uijfluct(:,1,3)                              ! ~
              uijfluct(:,2,3) = uijfluct(:,2,3) - 0.5*guumz(n,2)             ! ~
              uijfluct(:,3,2) = uijfluct(:,2,3)                              ! of sij
!
            endif
!
            if ( .not.lburgers_testflow ) then
!
              do i=1,3                                                       !
                U0testgu(:,i) = uijfluct(:,i,3)*gH0test                      ! S(u).grad(H^T)
              enddo
!
              call multsv(0.5*ghfluct(:,3),gU0test,ugU0test)                 ! S(U^T).grad(h)        MR: checked by numbers
              ugU0test(:,3) = ugU0test(:,3) + 0.5*ghfluct(:,3)*gU0test(:,3)  !  ~
              call multsv_mn_add(-(1./3.)*gU0test(:,3),ghfluct,ugU0test)     !  ~
!
              df(l1:l2,m,n,iuxtest:iuztest) =   df(l1:l2,m,n,iuxtest:iuztest) &
                                              + (2.*nutest*cs2test1)*(U0testgu+ugU0test)
            endif

          endif

        endif
!  add linear part of diffusion term nu*(del2u^pq + divu^pq/3)
!
        if (nutest/=0.) then
!
          if ( jtest>0 .or. .not.lprescribed_velocity ) then

            call del2v_etc(f,iuxtest,DEL2=del2utest,GRADDIV=graddivutest)
            if ( .false..and. m==10 .and. n==10 .and.jtest==1) then
              print'(14(e12.6,","))', del2utest(:,1)
              print'(14(e12.6,","))',f(l1:l2,m,n,iuxtest)
              print*, '================'
            endif
            df(l1:l2,m,n,iuxtest:iuztest) =   df(l1:l2,m,n,iuxtest:iuztest) &
                                            + nutest*(del2utest+(1./3.)*graddivutest)
          endif

        endif
!
! add continuous forcing
!
        if ( jtest==0 .and. lforcing_cont ) then                      ! requires that forcing is 'fluctuation': <force>=0

          if (.not.lprescribed_velocity) &
            df(l1:l2,m,n,iuxtest:iuztest) =   df(l1:l2,m,n,iuxtest:iuztest) &
                                            + ampl_fcont_uu*p%fcont
        endif
!
! add Coriolis/centrifugal force         (MR: precession still missing!!!)
!
        if ( jtest>0 .or. .not.lprescribed_velocity ) then
          call coriolis_cartesian(df,uutest,iuxtest)
!
! add Shear if considered a parameter
!
        if ( lshear .and. lshear_as_param ) &
          call shear_variables(f,df,ntestflow,iuutest,4,.true.)
!
!  check for testflow timestep
!
          if (lfirst.and.ldt) &
            advec_uu=max(advec_uu,sum(abs(uutest)*dline_1,2))
        endif

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
        maxadvec=maxadvec+advec_uu
        diffus_nu=nutest*dxyz_2
        maxdiffus=max(maxdiffus,diffus_nu)
      endif
!
      call xysum_mn_name_z(f(l1:l2,m,n,iuutest  ),idiag_ux0mz)!MR: only testflow # 0
      call xysum_mn_name_z(f(l1:l2,m,n,iuutest+1),idiag_uy0mz)
      call xysum_mn_name_z(f(l1:l2,m,n,iuutest+2),idiag_uz0mz)
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
        do jtest=0,njtestflow_loc
!
          if (idiag_upqrms(jtest)/=0) then
!
            call dot2(f(l1:l2,m,n,iuxtest:iuztest),upq2)
            !!if (m==5 .and. n==5.and.jtest==1) &
              !!print*, 'upq,jtest,iuxtest,iuztest:', &
                  !!jtest,iuxtest,iuztest,minval(f(l1:l2,m,n,iuxtest:iuztest)),maxval(f(l1:l2,m,n,iuxtest:iuztest))
!
            call sum_mn_name(upq2,idiag_upqrms(jtest),lsqrt=.true.)

          endif
!
          if (.not.lburgers_testflow.and.idiag_hpqrms(jtest)/=0) then
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
        if (lfirst_call) then
!
          call read_snaptime(trim(file),tuuinit,nuuinit,duuinit,t)
!
          if (tuuinit==0 .or. tuuinit < t-duuinit) then
            tuuinit=t+duuinit
          endif
          lfirst_call=.false.
!
        endif
!
        if (t >= tuuinit) then

          call init_uutest(f) !MR: isn't this meant???
          !!!call initialize_testflow(f,reinitialize_uutest)
          call update_snaptime(file,tuuinit,nuuinit,duuinit,t,ltestflow_out)
!
        endif
      endif
!
    endsubroutine duutest_dt
!***********************************************************************
    subroutine get_slices_testflow(f,slices)          ! -> Default interface
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
    subroutine testflow_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!   15-dec-10/MR: adapted from density
!
      use Hydro, only: remove_mean_momenta
      use Cdata
      use Mpicomm

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      if (lremove_mean_momenta_testflow) &
        call remove_mean_momenta(f,iuutest)

      if (lremove_mean_enthalpy_z) &
        call remove_mean_enthalpy_z(f)
!
    endsubroutine testflow_before_boundary
!***********************************************************************
    subroutine remove_mean_enthalpy_z(f)
!
!  removes xy average of enthaply in 0 solution.
!
!   15-dec-10/MR: coded
!
      use Cdata
      use Mpicomm, only: mpiallreduce_sum

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f

      real    :: hm, hm1, fac
      integer :: ihhtest,i,j

      fac = 1./(nxgrid*nygrid)

      ihhtest = iuutest+3
!
      do i=n1,n2

        hm = 0.
        do j=m1,m2
          hm = hm + sum(f(:,j,i,ihhtest))
        enddo

        if (nprocy>1) then                       ! no x parallelization allowed
          call mpiallreduce_sum(hm,hm1,idir=12)
          hm=hm1
        endif

        f(:,:,i,ihhtest) = f(:,:,i,ihhtest) - fac*hm

      enddo
!
    endsubroutine remove_mean_enthalpy_z
!***********************************************************************
    subroutine calc_ltestflow_nonlin_terms(f,df)          ! -> Default interface (as do_prepencilstep or so)
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
      use density, only: glnrhomz
      use Hydro, only: uumz
      use Mpicomm, only: mpiallreduce_sum
      use Forcing, only:forcing_cont
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar),    intent(inout) :: df
!
      real, dimension (3,0:njtestflow) :: unltestm, unltestm1      ! unltestm1, hnltestm1 local fields for MPI correspondence
      real, dimension (  0:njtestflow) :: hnltestm, hnltestm1
!
      real, dimension (nx,3)   :: uufluct,uutest, uu0, ghfluct, ghtest, gh0, sghtest, unltest, force
      real, dimension (nx,3,3) :: sijtest,uijtest,sij0,uij0
      real, dimension (nx)     :: divutest,hnltest, divufluct

      integer :: jtest,i,j,ju,k,iii,jtestu,jtesto
      integer :: iuxtest, iuztest,ihhtest
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
            if (lprescribed_velocity) then
!
              call forcing_cont(uufluct)
              f(l1:l2,m,n,iuxtest:iuztest) = uufluct

              ghfluct = 0.
              if (.not.lburgers_testflow) f(l1:l2,m,n,ihhtest) = 0.

            else

              uufluct = f(l1:l2,m,n,iuxtest:iuztest)                      ! ufluct = u0
              if (.not.lburgers_testflow) call grad(f,ihhtest,ghfluct)    ! ghfluct = grad(h0)

            endif

          else
!
!  calculate uufluct=U-Umean out of the main run
!
            uufluct=f(l1:l2,m,n,iux:iux+2)
!
            do j=1,3
              uufluct(:,j)=uufluct(:,j)-uumz(n,j)
            enddo
!
! calculate ghfluct=gh-ghmean out of the main run
!
            call grad(f,ilnrho,ghfluct)
            ghfluct(:,3)=ghfluct(:,3)-glnrhomz(n-n1+1)
!
          endif
!
!  do each of the njtestflow (=2, 4, or 6) test flows at a time
!
testloop: do jtest=0,njtestflow_loc                           ! jtest=0 : primary turbulence, testflow=0
!
!  velocity vector and enthalpy gradient of testflow solution
!
            uutest=f(l1:l2,m,n,iuxtest:iuztest)
            if (.not.lburgers_testflow) call grad(f,ihhtest,ghtest)                 ! grad(htest)
!
!  velocity gradient matrix and velocity divergence of testflow solution
!
            call gij(f,iuxtest,uijtest,1)               ! TBC
            call div_mn(uijtest,divutest,uutest)        ! TBC
!
!  calculate traceless strain tensor sijtest from uijtest

            call traceless_strain( uijtest, divutest, sijtest, uutest )
                                                               ! this parameter is without meaning in Cartesian geometry
            if ( jtest==0 ) then      ! primary turbulence

              gh0  = ghtest           ! save primary turbulence
              uu0  = uutest
              uij0 = uijtest
              sij0 = sijtest
!
              if (.not.lprescribed_velocity) then
!  u.grad(u) term and u.grad(h) term
!
                !!call multmv(uijtest,uutest,unltest)
                call u_dot_grad(f,iuxtest,uijtest,uutest,unltest,UPWIND=ltestflow_upw_uu)       ! (u0.grad)(u0)

                !!call dot_mn(uutest,ghtest,hnltest)
                if ( .not.lburgers_testflow ) &
                  call u_dot_grad(f,ihhtest,ghtest,uutest,hnltest,UPWIND=ltestflow_upw_lnrho)   ! u0.grad(h0)

              endif
!
            else
!
              !!call multmv(uijtest,uufluct,unltest)                                        ! (u.grad)(utest)
              call u_dot_grad(f,iuxtest,uijtest, uufluct,unltest,UPWIND=ltestflow_upw_uu)
!
!  initialize unltest and hnltest (unless lburgers_testflow)
!
              !!call multmv(uij0,uutest,unltest,.true.)                                     ! (utest.grad)(u0)
              call u_dot_grad(f,iuutest,uij0,uutest,unltest,ltestflow_upw_uu,.true.)
!
              if ( .not.lburgers_testflow ) then
!
                !!call dot_mn(uufluct,ghtest ,hnltest)
                call u_dot_grad(f,ihhtest,ghtest,uufluct,hnltest,ltestflow_upw_lnrho)         ! u.grad(htest)
!
                !!call dot_mn(uutest,gh0,hnltest,.true.)
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
              if ( jtest==0 ) then
!
                if (.not.lprescribed_velocity) &
                  call multmv(sijtest,ghtest,sghtest)                          ! S(u0).grad(h0)

              else

                call multmv(sijtest,ghfluct,sghtest)                           ! S(utest).grad(h)
                call multmv(sij0   ,ghtest ,sghtest,.true.)                    ! S(u0).grad(htest)
!
              endif
!
!  part of dissipative term
!
              if ( .not.lburgers_testflow .and. (jtest>0 .or. .not.lprescribed_velocity) ) &
                unltest = unltest - (2.*nutest*cs2test1)*sghtest
!
            endif
!
!  unless SOCA or other things, incrementing df with unltest
!
            if ( jtest==0.and..not.lprescribed_velocity .or. jtest>0.and..not.lsoca_testflow ) then

              df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest)-unltest      ! nonlinear parts stored in df
!
              if ( .not.lburgers_testflow ) &
                df(l1:l2,m,n,ihhtest)=df(l1:l2,m,n,ihhtest)-hnltest

            endif
!
!  Means of nonlinear terms always needed for 0th testflow (if not prescribed),
!  if not SOCA for all testflows, except from that if the
!  the turbulence coefficients have to be calculated.
!
            if ( (jtest==0.and..not.lprescribed_velocity) .or. &
                 (jtest>0 .and.(.not.lsoca_testflow .or. ldiagnos) ) ) then
!
!  sums of nonlinear parts (here sum sums over x extent only)
!
              do j=1,3
                unltestm(j,jtest)=unltestm(j,jtest)+sum(unltest(:,j))
              enddo

              if ( .not.lburgers_testflow ) &
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
          headtt=.false.           ! TBC
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

        if (nprocy>1) then

          !!if (ldiagnos) then
          !!print*, 'iproc,n=', iproc, n
          !!print*, hnltestm(0)
          !!endif

          call mpiallreduce_sum(unltestm,unltestm1,(/3,jtesto+1/),idir=2)
          unltestm(:,0:jtesto) = unltestm1(:,0:jtesto)
!
          call mpiallreduce_sum(hnltestm,hnltestm1,jtesto+1,idir=2)
          hnltestm(0:jtesto) = hnltestm1(0:jtesto)
!
        endif
!
        unltestm(:,0:jtesto)=fac*unltestm(:,0:jtesto)         ! means of nonlinear parts
        hnltestm(  0:jtesto)=fac*hnltestm(0:jtesto)
!
!  means are completely determined -> calculation of the  *fluctuations* of the nonlinear parts
!
        if (lprescribed_velocity) then
          jtestu=1
        else
          jtestu=0
        endif

        if ( lsoca_testflow ) then
          jtesto=0
        else
          jtesto=njtestflow
        endif
!
        do jtest=jtestu,jtesto
!
          iuxtest=iuutest+4*jtest
          ihhtest=iuxtest+3

          do j=1,3
!
            ju = iuxtest+j-1
!
!  subtract unltestm (mean flow)
!
            df(l1:l2,m1:m2,n,ju)=df(l1:l2,m1:m2,n,ju)+unltestm(j,jtest)
!
          enddo

          if ( .not.lburgers_testflow ) &
            df(l1:l2,m1:m2,n,ihhtest)=df(l1:l2,m1:m2,n,ihhtest)+hnltestm(jtest)

        enddo
!
        !!print*, 'unltestm, hnltestm:', minval(unltestm),maxval(unltestm), minval(hnltestm),maxval(hnltestm)
        !!if (ldiagnos) print*, 'hnltestm:', hnltestm     !!z(n), unltestm(1,1), unltestm(2,1), unltestm(1,2), unltestm(2,2)

        if (ldiagnos) call calc_coefficients(n,unltestm,hnltestm)
!
        lfirstpoint=.false.
!
      enddo zloop
!
! remove gal-parts from Fipq
!
      if ( ldiagnos .and. itestflow/='W11-W22' ) then
!
        do k=1,2

          gal1 = get_from_fname(idiag_galij(k,1))
          gal2 = get_from_fname(idiag_galij(k,2))

          do n=izrange(1),izrange(2)
!
            select case (itestflow)

            case('quadratic','quasi-periodic')
!
                call surf_mn_name( -wamp*gal2*z(n), idiag_aklamij(k,1) )
                call surf_mn_name(  wamp*gal1*z(n), idiag_aklamij(k,2) )

            case default
!
            end select
!
          enddo
!
          if ( itestflow=='quadratic' .or. itestflow=='quasi-periodic' ) then

            aklam1 = get_from_fname(idiag_aklamij(k,1))
            aklam2 = get_from_fname(idiag_aklamij(k,2))
!
            do n=izrange(1),izrange(2)
!
              call surf_mn_name( wamp*(-gal1*zq2(n)+aklam2*z(n)), idiag_nuij(k,1) )
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
    subroutine calc_coefficients(indz,Fipq,Qipq)
!
!  15-feb-13/MR: diagnostics completed
!
    use Cdata
!
    integer, intent(in) :: indz
!
    real, dimension (3,0:njtestflow):: Fipq            ! double index pq in single one (second) subsumed
    real, dimension (0:njtestflow)  :: Qipq
!
    real, dimension (2,2) :: aklam
    integer :: i,j,i3,i4,i5,i6,k
!
      Fipq=Fipq/(-wamp*valid_zrange*nzgrid*nprocy)                               ! as at call Fipq, Qipq have inverted sign yet
      Qipq=Qipq/(-wamp*valid_zrange*nzgrid*nprocy)                               ! factor nzgrid for averaging over z
!
      do j=1,njtestflow_loc
!
        do i=1,3
          call surf_mn_name( Fipq(i,j), idiag_Fipq(i,j) ) ! surf_mn_name because only simple addition needed
        enddo
!
        if (.not.lburgers_testflow) call surf_mn_name( Qipq(j), idiag_Qpq(j) )
!
      enddo
!
      select case (itestflow)

        case ('W11-W22')
!
          do k=1,3
!
!  calculate aka-lambda and nu tensors
!
            if ( idiag_aklamij(k,1)/=0 ) &
              call surf_mn_name(    cz(indz)*Fipq(k,1)+  sz(indz)*Fipq(k,2), idiag_aklamij(k,1) )
            if ( idiag_aklamij(k,2)/=0 ) &
              call surf_mn_name(    cz(indz)*Fipq(k,3)+  sz(indz)*Fipq(k,4), idiag_aklamij(k,2) )
!
            if ( idiag_nuij(k,1)/=0 ) &
              call surf_mn_name( -k1sz(indz)*Fipq(k,3)+k1cz(indz)*Fipq(k,4), idiag_nuij(k,1) )
            if ( idiag_nuij(k,2)/=0 ) &
              call surf_mn_name(  k1sz(indz)*Fipq(k,1)-k1cz(indz)*Fipq(k,2), idiag_nuij(k,2) )
            if ( idiag_nuij(k,3)/=0 ) &
              call surf_mn_name( -k1sz(indz)*Fipq(k,5)+k1cz(indz)*Fipq(k,6), idiag_nuij(k,3) )
!
            if ( idiag_gammai(k)/=0 ) &
              call surf_mn_name( cz(indz)*Fipq(k,5)+sz(indz)*Fipq(k,6), idiag_gammai(k) )
!
          enddo
!
          !!print*, 'indz=', indz, -k1sz(indz)*Fipq(1,3)+k1cz(indz)*Fipq(1,4), k1sz(indz)*Fipq(2,1)-k1cz(indz)*Fipq(2,2), &
          !!                       -k1sz(indz)*Fipq(2,3)+k1cz(indz)*Fipq(2,4), k1sz(indz)*Fipq(1,1)-k1cz(indz)*Fipq(1,2)
!
!  calculate aklamQ and nuQ vectors
!
          if (.not.lburgers_testflow) then
!
            if ( idiag_aklamQi(1)/=0 ) &
              call surf_mn_name( Qipq(1)*cz(indz)+Qipq(2)*sz(indz), idiag_aklamQi(1) )       ! \aklamQ_1
            if ( idiag_aklamQi(2)/=0 ) &
              call surf_mn_name( Qipq(3)*cz(indz)+Qipq(4)*sz(indz), idiag_aklamQi(2) )       ! \aklamQ_2
!
            if ( idiag_nuQi(1)/=0 ) &
              call surf_mn_name( -Qipq(3)*k1sz(indz)+Qipq(4)*k1cz(indz), idiag_nuQi(1) )      ! \nuQ_1
            if ( idiag_nuQi(2)/=0 ) &
              call surf_mn_name(  Qipq(1)*k1sz(indz)-Qipq(2)*k1cz(indz), idiag_nuQi(2) )      ! \nuQ_2
            if ( idiag_nuQi(3)/=0 ) &
              call surf_mn_name( -Qipq(5)*k1sz(indz)+Qipq(6)*k1cz(indz), idiag_nuQi(3) )      ! \nuQ_3
!
!  calculate gammaQ-scalar
!
            if ( idiag_gammaQ/=0 ) &
              call surf_mn_name( Qipq(5)*cz(indz)+Qipq(6)*sz(indz), idiag_gammaQ )           ! \gamma^Q
!
            if ( njtestflow>6 ) then

              do k=1,3
                if ( idiag_zetai(k)/=0 ) &
                  call surf_mn_name( -Fipq(k,7)*  cz(indz)-Fipq(k,8)*  sz(indz), idiag_zetai(k) )    ! \zeta_k
                if ( idiag_xii(k)/=0 ) &
                  call surf_mn_name(  Fipq(k,7)*k1sz(indz)-Fipq(k,8)*k1cz(indz), idiag_xii(k) )      ! \xi_k
              enddo
!
              if ( idiag_zetaQ/=0 ) &
                call surf_mn_name( -Qipq(7)*  cz(indz)-Qipq(8)*  sz(indz), idiag_zetaQ )    ! \zetaQ
              if ( idiag_xiQ/=0 ) &
                call surf_mn_name(  Qipq(7)*k1sz(indz)-Qipq(8)*k1cz(indz), idiag_xiQ )      ! \xiQ

            endif

          endif
!
        case('quadratic','quasi-periodic')
!
          if (idiag_gal/=0) then
!
            do i=1,2
            do j=1,2
              call surf_mn_name(Fipq(i,j),idiag_galij(i,j))          ! \gal_{ij}
            enddo
            enddo
!
            if ( indz>=izrange(1) .and. indz<=izrange(2) ) then
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
            endif
          endif

          if (.not.lburgers_testflow) then
!
            if ( idiag_nuQ/=0) then
!
              do i=1,2
                call surf_mn_name( -Qipq(i+2)+Qipq(i)*z(indz), idiag_nuQ+i-1 ) ! \nuQ_i   TBC
              enddo
!
            endif
!
          endif
!
        case ('onlyconstant')
!
          do i=1,2
            do j=1,2
              call surf_mn_name(Fipq(i,j),idiag_galij(i,j))                   ! \gal_{ij}
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
!
        case ('none')
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
    endsubroutine calc_coefficients
!***********************************************************************
    subroutine set_uutest(uutest,jtest)
!
!  set testflow
!
!   3-jun-05/axel: coded
!
      use Cdata

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
      case default; uutest=0.
!
      endselect
!
    endsubroutine set_uutest_U11_U21
!***********************************************************************
    subroutine set_U0test_W11_W22 (U0test,gU0test,gH0test,jtest)
!
!  set testflow
!
!  23-mar-08/axel: adapted from testflow_z
!  15-feb-13/MR: parameter gH0test/test cases 7,8 added
!
      use Cdata

      real, dimension (nx,3), intent(out) :: U0test,gU0test
      real, dimension (nx),   intent(out) :: gH0test
      integer,                intent(in)  :: jtest
!
!  set U0test and gU0test for each of the various cases
!
      gH0test=0.

      select case (jtest)
!
      case (1)
        U0test (:,1)=0.; U0test (:,2)=-wamp*k1sz(n); U0test (:,3)=0.
        gU0test(:,1)=0.; gU0test(:,2)=-wamp*cz(n);   gU0test(:,3)=0.
!
      case (2)
        U0test (:,1)=0.; U0test (:,2)=+wamp*k1cz(n); U0test (:,3)=0.
        gU0test(:,1)=0.; gU0test(:,2)=-wamp*sz(n);   gU0test(:,3)=0.
!
      case (3)
        U0test (:,1)=+wamp*k1sz(n); U0test (:,2)=0.; U0test( :,3)=0.
        gU0test(:,1)=+wamp*cz(n);   gU0test(:,2)=0.; gU0test(:,3)=0.
!
      case (4)
        U0test (:,1)=-wamp*k1cz(n); U0test (:,2)=0.; U0test (:,3)=0.
        gU0test(:,1)=+wamp*sz(n);   gU0test(:,2)=0.; gU0test(:,3)=0.
!
      case (5)
        U0test (:,1)=0.; U0test (:,2)=0.; U0test (:,3)=+wamp*k1sz(n)
        gU0test(:,1)=0.; gU0test(:,2)=0.; gU0test(:,3)=+wamp*cz(n)
!
      case (6)
        U0test (:,1)=0.; U0test (:,2)=0.; U0test (:,3)=-wamp*k1cz(n)
        gU0test(:,1)=0.; gU0test(:,2)=0.; gU0test(:,3)=+wamp*sz(n)

      case (7)
        U0test=0.; gU0test=0.
        gH0test=-wamp*cz(n)
!
      case (8)
        U0test=0.; gU0test=0.
        gH0test=-wamp*sz(n)

      case default
        U0test=0.; gU0test=0.

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
      integer :: comp, indx

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
      case(P1Q2,P2Q2)
!        U0test(:,1)=wamp*(1.+z(n)-zc(n)); U0test(:,2)=0.; U0test(:,3)=0.
!       gU0test(:,1)=wamp-...; gU0test(:,2)=0.; gU0test(:,3)=0.
!
        comp=jtest-P1Q2+1

        if ( z(n)<zrange(1) ) then
          indx=1
        else if ( z(n)>zrange(2) ) then
          indx=2
        else
          indx=0
        endif

!      case (P2Q2)
!        U0test(:,1)=0.; U0test(:,2)=wamp*(1.+z(n)-zc(n)); U0test(:,3)=0.
!        gU0test(:,1)=0.; gU0test(:,2)=wamp-...; gU0test(:,3)=0.
!
        if (indx>0) then
          U0test (:,comp)=wamp*(ca(indx)+cl(indx)*z(n)+   cq(indx)*z(n)*z(n))
          gU0test(:,comp)=wamp*(         cl(indx)     +2.*cq(indx)*z(n)     )
        else
          U0test(:,comp)=wamp*z(n)
          gU0test(:,comp)=wamp
        endif

!      case (P1Q3)
!        U0test(:,1)=wamp*(1.-z(n)+zc(n)); U0test(:,2)=0.; U0test(:,3)=0.
!        gU0test(:,1)=-wamp+...; gU0test(:,2)=0.; gU0test(:,3)=0.
!
        !!print*, comp,indx,zu,z(n),zo,U0test(:,comp),gU0test(:,comp)

!      case (P2Q3)
!        U0test(:,1)=0.; U0test(:,2)=wamp*(1.-z(n)+zc(n)); U0test(:,3)=0.
!        gU0test(:,1)=0.; gU0test(:,2)=-wamp+...; gU0test(:,3)=0.
!
        U0test(:,3-comp)=0.; U0test(:,3-comp)=0.
        gU0test(:,3-comp)=0.; gU0test(:,3-comp)=0.

      case(P1Q3,P2Q3)

        comp=jtest-P1Q3+1

        if ( z(n)<zrange(1) ) then
          indx=1
        else if ( z(n)>zrange(2) ) then
          indx=2
        else
          indx=0
        endif

        if (indx==1) then

          U0test (:,comp)=0.5*wamp*(cqa(indx)+cqq(indx)*(xyz0(3)-z(n))**2)
          gU0test(:,comp)=   -wamp*           cqq(indx)*(xyz0(3)-z(n))

        else if (indx==2) then

          U0test (:,comp)=0.5*wamp*(cqa(indx)+cqq(indx)*(xyz1(3)-z(n))**2)
          gU0test(:,comp)=   -wamp*           cqq(indx)*(xyz1(3)-z(n))

        else
          U0test (:,comp)=0.5*wamp*z(n)**2
          gU0test(:,comp)=    wamp*z(n)
        endif

        !!print*, comp,indx,zu,z(n),zo,U0test(:,comp),gU0test(:,comp)

        U0test(:,3-comp)=0.; U0test(:,3-comp)=0.
        gU0test(:,3-comp)=0.; gU0test(:,3-comp)=0.

        U0test=0.;gU0test=0.!!!!

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
    subroutine rprint_testflow(lreset,lwrite)                   ! -> Default interface
!
!  reads and registers print parameters relevant for testflow fields
!
!   3-jun-05/axel: adapted from rprint_magnetic
!
      use Cdata
      use FArrayManager, only: farray_index_append
      use General, only: loptest
!
      logical           :: lreset
      logical, optional :: lwrite
!
      integer           :: iname,i,j,p,ifound,ifoundold,ijmax
      character         :: cind
      character(len=2)  :: cind2
      character(len=1)  :: cind1
      character(len=20) :: name
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
!
        idiag_Fipq=0; idiag_Qpq=0
!
        idiag_gal=0  ; idiag_galij=0
        idiag_aklam=0; idiag_aklamij=0
        idiag_gamma=0; idiag_gammai=0
        idiag_nu=0   ; idiag_nuij=0
!
        idiag_zeta=0; idiag_zetai=0
        idiag_xi=0; idiag_xii=0

        idiag_aklamQ=0; idiag_aklamQi=0
        idiag_nuQ=0; idiag_nuQi=0

        idiag_zetaQ=0; idiag_xiQ=0; idiag_gammaQ=0

        idiag_ux0mz=0; idiag_uy0mz=0; idiag_uz0mz=0
        idiag_upqrms=0; idiag_hpqrms=0
!
      endif
!
!  check for those quantities that we want to evaluate online
!
      ifoundold=0
!
      do iname=1,nname
!
        ifound=0
        if ( itestflow/='W11-W22' ) &
          ifound = fparse_name(iname,cname(iname),'gal',idiag_gal,cform(iname))

        ifound = ifound + fparse_name(iname,cname(iname),'aklam',idiag_aklam,cform(iname))
        ifound = ifound + fparse_name(iname,cname(iname),'gamma',idiag_gamma,cform(iname))
        ifound = ifound + fparse_name(iname,cname(iname),'nu',idiag_nu,cform(iname))
        ifound = ifound + fparse_name(iname,cname(iname),'zeta',idiag_zeta,cform(iname))
        ifound = ifound + fparse_name(iname,cname(iname),'xi',idiag_xi,cform(iname))

        ifound = ifound + fparse_name(iname,cname(iname),'aklamQ',idiag_aklamQ,cform(iname))
        ifound = ifound + fparse_name(iname,cname(iname),'gammaQ',idiag_gammaQ,cform(iname))
        ifound = ifound + fparse_name(iname,cname(iname),'nuQ',idiag_nuQ,cform(iname))
        ifound = ifound + fparse_name(iname,cname(iname),'zetaQ',idiag_zetaQ,cform(iname))
        ifound = ifound + fparse_name(iname,cname(iname),'xiQ',idiag_xiQ,cform(iname))
!
        do i=1,3
          do j=1,3
!
            cind2 = gen_ind(i,j)

            if ( itestflow/='W11-W22' ) &
              ifound = ifound + fparse_name(iname,cname(iname),'gal'//cind2,idiag_galij(i,j),cform(iname))

            if (j<3) ifound = ifound + fparse_name(iname,cname(iname),'aklam'//cind2,idiag_aklamij(i,j),cform(iname))
            ifound = ifound + fparse_name(iname,cname(iname),'nu'//cind2,idiag_nuij(i,j),cform(iname))
!
          enddo

          cind1 = gen_ind(i)
          ifound = ifound + fparse_name(iname,cname(iname),'gamma'//cind1,idiag_gammai(i),cform(iname))
          ifound = ifound + fparse_name(iname,cname(iname),'zeta'//cind1,idiag_zetai(i),cform(iname))
          ifound = ifound + fparse_name(iname,cname(iname),'xi'//cind1,idiag_xii(i),cform(iname))

          if (i<3) ifound = ifound + fparse_name(iname,cname(iname),'aklamQ'//cind1,idiag_aklamQi(i),cform(iname))
          ifound = ifound + fparse_name(iname,cname(iname),'nuQ'//cind1,idiag_nuQi(i),cform(iname))

        enddo

        p=1
        do j=0,njtestflow
!
          if (j==0) then
            name = '0rms'
          else
!
            cind2 = gen_ind(p,floor((j+1)/2.))
            name = cind2//'rms'
!
            p=3-p

          endif
!
          ifound = ifound + fparse_name(iname,cname(iname),'u'//name,idiag_upqrms(j),cform(iname))
          ifound = ifound + fparse_name(iname,cname(iname),'h'//name,idiag_hpqrms(j),cform(iname))

          if ( j>0 ) then
!
            do i=1,3
              ifound = ifound + fparse_name(iname,cname(iname),'F'//cind2,idiag_Fipq(i,j),cform(iname))
            enddo
!
            cind = gen_ind(j)
            ifound = ifound + fparse_name(iname,cname(iname),'Q'//cind,idiag_Qpq(j),cform(iname))

          endif

        enddo
!
        idiag_map(iname) = iname       ! initialising the index mapping vector for cname und cform

        if ( ifound==0 .and. ifoundold>0 ) &
          print*, 'testflow_z, rprint_testflow: Warning - diagnostic ouput for line',iname, &
                  ' and beyond will be messed up!!!'
!
        if ( ifound/=0 ) ifoundold=ifound
!
      enddo
      nname_old = nname
!
      if (njtestflow==4) then
        ijmax=2
      else
        ijmax=3
      endif

      call update_diag( idiag_gal,   idiag_galij,   'gal'  , ijmax)
      call update_diag( idiag_aklam, idiag_aklamij, 'aklam', ijmax)
      call update_diag( idiag_nu,    idiag_nuij,    'nu'   , ijmax)
      call update_diag( idiag_gamma, idiag_gammai,  'gamma', ijmax)
      call update_diag( idiag_zeta,  idiag_zetai,   'zeta' , ijmax)
      call update_diag( idiag_xi,    idiag_xii,     'xi'   , ijmax)

      call update_diag( idiag_aklamQ,idiag_aklamQi, 'aklamQ', ijmax)
      call update_diag( idiag_nuQ,   idiag_nuQi,    'nuQ'   , ijmax)

      if ( idiag_gal/=-1   ) call correct_inds(idiag_galij, ijmax)
      if ( idiag_aklam/=-1 ) call correct_inds(idiag_aklamij, ijmax)
      if ( idiag_nu/=-1    ) call correct_inds(idiag_nuij, ijmax)
!
      if ( idiag_gamma/=-1 ) call correct_inds(idiag_gammai, ijmax)
      if ( idiag_zeta/=-1  ) call correct_inds(idiag_zetai, ijmax)
      if ( idiag_xi/=-1    ) call correct_inds(idiag_xii, ijmax)

      if ( idiag_aklamQ/=-1) call correct_inds(idiag_aklamQi)
      if ( idiag_nuQ/=-1   ) call correct_inds(idiag_nuQi, ijmax)

      if ( idiag_gammaQ/=-1) call correct_inds(idiag_gammaQ)
      if ( idiag_zetaQ/=-1 ) call correct_inds(idiag_zetaQ)
      if ( idiag_xiQ/=-1   ) call correct_inds(idiag_xiQ)

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
      if (loptest(lwrite)) then
        call farray_index_append('iuutest',iuutest)
        call farray_index_append('ntestflow',ntestflow)
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
    character(len=1) function gen_1ind(i)
!
      integer :: i
!
      character(len=20) c
      write(c,fmt='(i19)') i

      gen_1ind = trim(adjustl(c))

    endfunction gen_1ind
!***********************************************************************
    subroutine mark_del_elems2( carray, leng, indices )

      integer, dimension(:,:) :: indices
      integer                 :: leng
!
      character*(*), dimension(*) :: carray
!
      intent(inout) :: indices, carray, leng
!
      integer :: i,j,ind
!
      do i=1,size(indices,1)
      do j=1,size(indices,2)

        ind = indices(i,j)
        if ( ind /= 0 ) then
!
          call del_elem( carray, idiag_map(ind), leng )
!
          idiag_map(ind) = 0
          if ( ind<nname_old ) &
            idiag_map(ind+1:nname_old) = idiag_map(ind+1:nname_old)-1

        endif
      enddo
      enddo
!
    endsubroutine mark_del_elems2
!***********************************************************************
    subroutine mark_del_elems1( carray, leng, indices )

      integer, dimension(:) :: indices
      integer               :: leng
!
      character*(*), dimension(*) :: carray
!
      intent(inout) :: indices, carray, leng
!
      integer :: i,ind
!
      do i=1,size(indices)

        ind = indices(i)
        if ( ind /= 0 ) then
!
          call del_elem( carray, idiag_map(ind), leng )
!
          idiag_map(ind) = 0
          if ( ind<nname_old ) &
            idiag_map(ind+1:nname_old) = idiag_map(ind+1:nname_old)-1

        endif

      enddo
!
    endsubroutine mark_del_elems1
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
    subroutine insert_array( carray, cinsert, index, leng )
!
! 15-feb-2013/MR: parameter leng_insert removed (now derived from cinsert)
!
      integer                     :: index, leng, leng_insert, i
      character*(*), dimension(*) :: carray
      character*(*), dimension(:) :: cinsert
!
      intent(in)    :: index, cinsert
      intent(inout) :: leng, carray
!
      if ( index>0.and.index<=leng+1 ) then
!
        leng_insert = size(cinsert)
        do i=leng,index,-1
          carray(i+leng_insert) = carray(i)
        enddo
!
        carray(index:index+leng_insert-1) = cinsert
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
    subroutine update_diag2( mat_diag_ind, mat_diag_indices, cdiagname, ijmax )
!
! 15-feb-2013/MR: optional parameter ijmax for upper bound
!                 of mat_diag_indices added; variable size
!                 of mat_diag_indices handled
!
      use Cdata
      use Diagnostics
!
      integer                 :: mat_diag_ind, ijmax
      integer, dimension(:,:) :: mat_diag_indices
!
      intent(in)    :: ijmax
      intent(inout) :: mat_diag_ind, mat_diag_indices
!
      character*(*), intent(in) :: cdiagname
!
      character(len=30) cformat
!
      integer :: iname, i, j, nname_form, indx, matsize1, matsize2, addsize
      character cind1, cind2
      character(len=60), dimension(:), allocatable :: entries
!
      matsize1=min(size(mat_diag_indices,1),ijmax)
      matsize2=min(size(mat_diag_indices,2),ijmax)

      if ( mat_diag_ind==0 ) then                                 ! if complete tensor was not already inquired
!
        do iname=1,nname
          do i=1,matsize1
            do j=1,matsize2
!
              if ( mat_diag_indices(i,j)/=0 ) then                ! if any of the elements is inquired
!
                mat_diag_ind = -2                                 ! flag the whole tensor to be calculated
                return
!
              endif
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

        call del_elem( cname, indx, nname )                       ! replace inquiry for the whole tensor by inquiries for

        allocate(entries(matsize1*matsize2))
        ind=1
        do i=1,matsize1
          cind1=gen_ind(i)
          do j=1,matsize2
            cind2=gen_ind(j)
            entries(ind) = cdiagname//cind1//cind2//cformat
            ind=ind+1
          enddo
        enddo

        call insert( cname,entries,indx,nname)  ! by inquiries for all elements
        deallocate(entries)

        do i=1,matsize1
          do j=1,matsize2
!
            mat_diag_indices(i,j) = indx
            indx = indx+1
!
          enddo
        enddo

        addsize=matsize1*matsize2-1
        if ( mat_diag_ind<nname_old ) &
          idiag_map(mat_diag_ind+1:nname_old) = idiag_map(mat_diag_ind+1:nname_old)+addsize
!
        cformat = cform(indx)
        call insert( cform, cformat, addsize, indx+1, nname_form )
!
        mat_diag_ind = -1
!
      endif
!
    endsubroutine update_diag2
!***********************************************************************
    subroutine update_diag1( mat_diag_ind, mat_diag_indices, cdiagname, ijmax )
!
! 15-feb-2013/MR: optional parameter ijmax for upper bound
!                 of mat_diag_indices added; variable size
!                 of mat_diag_indices handled
!
      use Cdata
      use Diagnostics
!
      integer               :: mat_diag_ind, ijmax
      integer, dimension(:) :: mat_diag_indices
!
      intent(in)    :: ijmax
      intent(inout) :: mat_diag_ind, mat_diag_indices
!
      character*(*), intent(in) :: cdiagname
!
      character(len=30) cformat
!
      integer :: iname, i, nname_form, indx, matsize
!
      matsize=min(size(mat_diag_indices),ijmax)
!
      if ( mat_diag_ind==0) then             ! if complete tensor was not already inquired
!
        do iname=1,nname
        do i=1,matsize
!
          if ( mat_diag_indices(i)/=0 ) then ! if any of the elements is inquired
!
            mat_diag_ind = -2                ! flag the whole tensor to be calculated
            return
!
          endif
!
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
        if (matsize==2) then
          call insert( cname, (/cdiagname//'1'//cformat,cdiagname//'2'//cformat/), indx, nname )    ! all elements
        else
          call insert( cname, (/cdiagname//'1'//cformat,cdiagname//'2'//cformat,cdiagname//'3'//cformat/), indx, nname )   ! all elements
        endif

        do i=1,matsize
!
          mat_diag_indices(i) = indx
          indx = indx+1
!
        enddo
!
        if ( mat_diag_ind<nname_old ) &
          idiag_map(mat_diag_ind+1:nname_old) = idiag_map(mat_diag_ind+1:nname_old)+matsize-1
!
        cformat = cform(indx)
        call insert( cform, cformat, matsize-1, indx+1, nname_form )
!
        mat_diag_ind = -1
!
      endif
!
    endsubroutine update_diag1
!***********************************************************************
    subroutine correct_2inds(mat_diag_indices, ijmax)
!
! 15-feb-2013/MR: optional parameter ijmax for upper bounds
!                 of mat_diag_indices added
!
      integer, dimension(:,:) :: mat_diag_indices
      integer, optional :: ijmax
!
      integer i,j,ind,ijm
!
      if (.not.present(ijmax)) then
        ijm=3
      else
        ijm=ijmax
      endif
!
      do i=1,min(size(mat_diag_indices,1),ijm)
        do j=1,min(size(mat_diag_indices,2),ijm)
!
          ind = mat_diag_indices(i,j)
!
          if (ind>0) &
            mat_diag_indices(i,j) = idiag_map(ind)
!
        enddo
      enddo
!
    endsubroutine correct_2inds
!***********************************************************************
    subroutine correct_1inds(mat_diag_indices, ijmax)
!
! 15-feb-2013/MR: optional parameter ijmax for upper bound
!                 of mat_diag_indices added
!
      integer, dimension(:) :: mat_diag_indices
      integer, optional :: ijmax
!
      integer i,ind,ijm
!
      if (.not.present(ijmax)) then
        ijm=3
      else
        ijm=ijmax
      endif
!
      do i=1,min(size(mat_diag_indices),ijm)
!
        ind = mat_diag_indices(i)
!
        if (ind>0) &
          mat_diag_indices(i) = idiag_map(ind)
!
      enddo
!
    endsubroutine correct_1inds
!***********************************************************************
    subroutine correct_0inds(mat_diag_indices)
!
! 15-feb-2013/MR: adapted frome correct_2inds for
!                 scalar quantities
!
      integer :: mat_diag_indices
!
      if (mat_diag_indices>0) &
        mat_diag_indices = idiag_map(mat_diag_indices)
!
    endsubroutine correct_0inds
!***********************************************************************
      SUBROUTINE ISORTP(N,RA,IP)
!
      INTEGER, INTENT(IN)  :: N
      INTEGER, INTENT(IN)  :: RA(*)
      INTEGER, INTENT(OUT) :: IP(*)

      INTEGER RRA,L,IR,I,J,IPA,J2,K
!
      IP(1) = 1
      IF (N<=1) GOTO 1001
      L = N/2+1
      IR = N
!
      DO I=2,N
         IP(I) = I
      ENDDO
!
   21 CONTINUE
         IF (IR>1) THEN
            IF (L>1) THEN
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
               IF (J2<=IR) THEN
                  K = J2
                  IF (K<IR) THEN
                     IF ( RA(IP(K))<RA(IP(K+1)) ) K = K+1
                  END IF
                  IF ( RRA<RA(IP(K)) ) THEN
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
