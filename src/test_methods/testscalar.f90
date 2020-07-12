! $Id$
!
!  This modules deals with all aspects of testscalar fields; if no
!  testscalar fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testscalar relevant subroutines listed in here.
!
!  Note: this routine requires that MVAR and MAUX contributions
!  together with njtestscalar are set correctly in the cparam.local file.
!  njtestscalar must be set at the end of the file such that njtestscalar=MVAR.
!
!  Example:
!  ! MVAR CONTRIBUTION 2
!  ! MAUX CONTRIBUTION 2
!  integer, parameter :: njtestscalar=2
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestscalar = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Testscalar
!
  use Cparam
  use Messages
!
  implicit none
!
  include '../testscalar.h'
!
!  cosine and sine function for setting test fields and analysis
!
  real, dimension(nx) :: cx,sx
  real, dimension(my) :: cy,sy
  real, dimension(mz) :: cz,sz
!
  character (len=labellen), dimension(ninit) :: initcctest='nothing'
  real, dimension (ninit) :: kx_cctest=1.,ky_cctest=1.,kz_cctest=1.
  real, dimension (ninit) :: phasex_cctest=0.,phasez_cctest=0.
  real, dimension (ninit) :: amplcctest=0.

  ! input parameters
  real, dimension (nx) :: ccc
  real :: tccinit=0.,dccinit=0.,tccinit_previous=0.
  logical :: reinitialize_cctest=.false.
  logical :: zextent=.true.,lsoca_ug=.false.,lset_cctest2=.false.
  logical :: lug_as_aux=.false.,linit_cctest=.false.
  logical :: lignore_ugtestm=.false.
  character (len=labellen) :: itestscalar='G1-G2'
  real :: ktestscalar=1., ktestscalar1=1.
  real :: kxtestscalar=1., kxtestscalar1=1.
  real :: kytestscalar=1., kytestscalar1=1.
  real :: lam_testscalar=0.,om_testscalar=0.,delta_testscalar=0.
  real :: delta_testscalar_next=0.
  integer, parameter :: mtestscalar=njtestscalar
  integer, parameter :: jtestz1=1,jtestz2=2,jtestx1=3,jtestx2=4,jtesty1=5,jtesty2=6
  integer :: nccinit
  real :: camp=1.,camp1=1.
!
  namelist /testscalar_init_pars/ &
       zextent,initcctest, &
       amplcctest,kx_cctest,ky_cctest,kz_cctest, &
       phasex_cctest,phasez_cctest, &
       lug_as_aux

  ! run parameters
  real :: kappatest=0.,kappatest1=0.
  real, dimension(njtestscalar) :: rescale_cctest=0.
  logical :: ltestscalar_ugc=.false.
  logical :: ltestscalar_newz=.true.
  logical :: ltestscalar_newx=.false.
  logical :: ltestscalar_newy=.false.
  logical :: ltestscalar_per_unitvolume=.false.
  namelist /testscalar_run_pars/ &
       reinitialize_cctest,zextent,lsoca_ug, &
       lset_cctest2,kappatest,kappatest1,itestscalar, &
       kxtestscalar,kytestscalar,ktestscalar, &
       lam_testscalar,om_testscalar,delta_testscalar, &
       ltestscalar_newx,ltestscalar_newz, &
       ltestscalar_per_unitvolume, &
       lug_as_aux,lignore_ugtestm, &
       dccinit,linit_cctest,camp, &
       rescale_cctest

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_gam11=0      ! DIAG_DOC: $\gamma_{1}^{(1)}$
  integer :: idiag_gam12=0      ! DIAG_DOC: $\gamma_{2}^{(1)}$
  integer :: idiag_gam13=0      ! DIAG_DOC: $\gamma_{3}^{(1)}$
  integer :: idiag_gam21=0      ! DIAG_DOC: $\gamma_{1}^{(2)}$
  integer :: idiag_gam22=0      ! DIAG_DOC: $\gamma_{2}^{(2)}$
  integer :: idiag_gam23=0      ! DIAG_DOC: $\gamma_{3}^{(2)}$
  integer :: idiag_gam31=0      ! DIAG_DOC: $\gamma_{1}^{(3)}$
  integer :: idiag_gam32=0      ! DIAG_DOC: $\gamma_{2}^{(3)}$
  integer :: idiag_gam33=0      ! DIAG_DOC: $\gamma_{3}^{(3)}$
  integer :: idiag_kap11=0      ! DIAG_DOC: $\kappa_{11}$
  integer :: idiag_kap21=0      ! DIAG_DOC: $\kappa_{21}$
  integer :: idiag_kap31=0      ! DIAG_DOC: $\kappa_{31}$
  integer :: idiag_kap12=0      ! DIAG_DOC: $\kappa_{12}$
  integer :: idiag_kap22=0      ! DIAG_DOC: $\kappa_{22}$
  integer :: idiag_kap32=0      ! DIAG_DOC: $\kappa_{32}$
  integer :: idiag_kap13=0      ! DIAG_DOC: $\kappa_{13}$
  integer :: idiag_kap23=0      ! DIAG_DOC: $\kappa_{23}$
  integer :: idiag_kap33=0      ! DIAG_DOC: $\kappa_{33}$
  integer :: idiag_gam11z=0     ! DIAG_DOC: $\gamma_{1}^{(1)}(z,t)$
  integer :: idiag_gam12z=0     ! DIAG_DOC: $\gamma_{2}^{(1)}(z,t)$
  integer :: idiag_gam13z=0     ! DIAG_DOC: $\gamma_{3}^{(1)}(z,t)$
  integer :: idiag_gam21z=0     ! DIAG_DOC: $\gamma_{1}^{(2)}(z,t)$
  integer :: idiag_gam22z=0     ! DIAG_DOC: $\gamma_{2}^{(2)}(z,t)$
  integer :: idiag_gam23z=0     ! DIAG_DOC: $\gamma_{3}^{(2)}(z,t)$
  integer :: idiag_gam31z=0     ! DIAG_DOC: $\gamma_{1}^{(3)}(z,t)$
  integer :: idiag_gam32z=0     ! DIAG_DOC: $\gamma_{2}^{(3)}(z,t)$
  integer :: idiag_gam33z=0     ! DIAG_DOC: $\gamma_{3}^{(3)}(z,t)$
  integer :: idiag_kap11z=0     ! DIAG_DOC: $\kappa_{11}(z,t)$
  integer :: idiag_kap21z=0     ! DIAG_DOC: $\kappa_{21}(z,t)$
  integer :: idiag_kap31z=0     ! DIAG_DOC: $\kappa_{31}(z,t)$
  integer :: idiag_kap12z=0     ! DIAG_DOC: $\kappa_{12}(z,t)$
  integer :: idiag_kap22z=0     ! DIAG_DOC: $\kappa_{22}(z,t)$
  integer :: idiag_kap32z=0     ! DIAG_DOC: $\kappa_{32}(z,t)$
  integer :: idiag_kap13z=0     ! DIAG_DOC: $\kappa_{13}(z,t)$
  integer :: idiag_kap23z=0     ! DIAG_DOC: $\kappa_{23}(z,t)$
  integer :: idiag_kap33z=0     ! DIAG_DOC: $\kappa_{33}(z,t)$
  integer :: idiag_mgam33=0     ! DIAG_DOC: $\tilde\gamma_{33}$
  integer :: idiag_mkap33=0     ! DIAG_DOC: $\tilde\kappa_{33}$
  integer :: idiag_ngam33=0     ! DIAG_DOC: $\hat\gamma_{33}$
  integer :: idiag_nkap33=0     ! DIAG_DOC: $\hat\kappa_{33}$
  integer :: idiag_c1rms=0      ! DIAG_DOC: $\left<c_{1}^2\right>^{1/2}$
  integer :: idiag_c2rms=0      ! DIAG_DOC: $\left<c_{2}^2\right>^{1/2}$
  integer :: idiag_c3rms=0      ! DIAG_DOC: $\left<c_{3}^2\right>^{1/2}$
  integer :: idiag_c4rms=0      ! DIAG_DOC: $\left<c_{4}^2\right>^{1/2}$
  integer :: idiag_c5rms=0      ! DIAG_DOC: $\left<c_{5}^2\right>^{1/2}$
  integer :: idiag_c6rms=0      ! DIAG_DOC: $\left<c_{6}^2\right>^{1/2}$
  integer :: idiag_c1pt=0       ! DIAG_DOC: $c^{1}$
  integer :: idiag_c2pt=0       ! DIAG_DOC: $c^{2}$
  integer :: idiag_c3pt=0       ! DIAG_DOC: $c^{3}$
  integer :: idiag_c4pt=0       ! DIAG_DOC: $c^{4}$
  integer :: idiag_c5pt=0       ! DIAG_DOC: $c^{5}$
  integer :: idiag_c6pt=0       ! DIAG_DOC: $c^{6}$
  integer :: idiag_F11z=0       ! DIAG_DOC: ${\cal F}_1^{1}$
  integer :: idiag_F21z=0       ! DIAG_DOC: ${\cal F}_2^{1}$
  integer :: idiag_F31z=0       ! DIAG_DOC: ${\cal F}_3^{1}$
  integer :: idiag_F12z=0       ! DIAG_DOC: ${\cal F}_1^{2}$
  integer :: idiag_F22z=0       ! DIAG_DOC: ${\cal F}_2^{2}$
  integer :: idiag_F32z=0       ! DIAG_DOC: ${\cal F}_3^{2}$
!
!  arrays for horizontally averaged uc
!
  real, dimension (nz,mtestscalar) :: ugtestm
  real, dimension (nx,mtestscalar) :: ugtestmx
  real, dimension (ny,mtestscalar) :: ugtestmy

  contains

!***********************************************************************
    subroutine register_testscalar()
!
!  Initialise variables which should know that we solve for the vector
!  potential: icctest, etc; increase nvar accordingly
!
!  26-nov-08/axel: adapted from testfield_z.f90
!  17-apr-09/axel: added y column of kappa tensor
!
      use Cdata
      use Mpicomm, only: stop_it
!
      integer :: j
!
!  Set first and last index of text field
!  Note: iaxtest, iaytest, and iaztest are initialized to the first test field.
!  These values are used in this form in start, but later overwritten.
!
      icctest=nvar+1
      icctestpq=icctest+njtestscalar-1
      ntestscalar=mtestscalar
      nvar=nvar+ntestscalar
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_testscalar: nvar = ', nvar
        print*, 'register_testscalar: icctest = ', icctest
      endif
!
!  Put variable names in array
!
      do j=icctest,nvar
        varname(j) = 'cctest'
      enddo
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_testscalar: nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',cctest $'
          if (nvar == mvar) write(4,*) ',cctest'
        else
          write(4,*) ',cctest $'
        endif
        write(15,*) 'cctest = fltarr(mx,my,mz,ntestscalar)*one'
      endif
!
    endsubroutine register_testscalar
!***********************************************************************
    subroutine initialize_testscalar(f)
!
!  Perform any post-parameter-read initialization
!
!  26-nov-08/axel: adapted from testfield_z.f90
!  27-dec-08/axel: extended to x-dependent mean fields
!
      use Cdata
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx) :: xtestscalar
      real, dimension(my) :: ytestscalar
      real, dimension(mz) :: ztestscalar
      real :: kxtestscalar_effective,kytestscalar_effective
      real :: ktestscalar_effective
      integer :: jtest, jcctest
!
!  Precalculate kxappatest if 1/kappatest (==kappatest1) is given instead
!
      if (kappatest1/=0.) then
        kappatest=1./kappatest1
      endif
!
!  set cosine and sine function for setting test fields and analysis
!
!  Choice of using rescaled x-array or original x-array
!  Define ktestscalar_effective to deal with boxes bigger than 2pi.
!
      if (ltestscalar_newx) then
        kxtestscalar_effective=kxtestscalar*(2.*pi/Lx)
        xtestscalar=kxtestscalar_effective*(x(l1:l2)-x0)-pi
      else
        kxtestscalar_effective=kxtestscalar
        xtestscalar=x(l1:l2)
      endif
      cx=cos(kxtestscalar*xtestscalar)
      sx=sin(kxtestscalar*xtestscalar)
!
!  Choice of using rescaled y-array or original y-array
!  Define kytestscalar_effective to deal with boxes bigger than 2pi.
!
      if (ltestscalar_newy) then
        kytestscalar_effective=kytestscalar*(2.*pi/Ly)
        ytestscalar=kytestscalar_effective*(y-y0)-pi
      else
        kytestscalar_effective=kytestscalar
        ytestscalar=y
      endif
      cy=cos(kytestscalar*ytestscalar)
      sy=sin(kytestscalar*ytestscalar)
!
!  Choice of using rescaled z-array or original z-array
!  Define kztestscalar_effective to deal with boxes bigger than 2pi.
!
      if (ltestscalar_newz) then
        ktestscalar_effective=ktestscalar*(2.*pi/Lz)
        ztestscalar=ktestscalar_effective*(z-z0)-pi
      else
        ktestscalar_effective=ktestscalar
        ztestscalar=z
      endif
      cz=cos(ktestscalar_effective*ztestscalar)
      sz=sin(ktestscalar_effective*ztestscalar)
!
!  Optionally, one can determine the phase in the actual field
!  and modify the following calculations in testscalar_after_boundary.
!
!  debug output
!
      if (lroot.and.ip<14) then
        print*,'x0=',x0
        print*,'cx=',cx
        print*,'sx=',sx
        print*,'y0=',y0
        print*,'cy=',cy
        print*,'sy=',sy
        print*,'z0=',z0
        print*,'cz=',cz
        print*,'sz=',sz
      endif
!
!  Also calculate its inverse, but only if different from zero
!  first kx
!
      if (kxtestscalar==0) then
        kxtestscalar1=1.
      else
        kxtestscalar1=1./kxtestscalar_effective
      endif
!
      if (ktestscalar==0) then
        ktestscalar1=1.
      else
        ktestscalar1=1./ktestscalar_effective
      endif
!
!  next ky
!
      if (kytestscalar==0) then
        kytestscalar1=1.
      else
        kytestscalar1=1./kytestscalar_effective
      endif
!
!  finally kz (still called k)
!
      if (ktestscalar==0) then
        ktestscalar1=1.
      else
        ktestscalar1=1./ktestscalar_effective
      endif
!
!  calculate inverse testscalar amplitude (unless it is set to zero)
!
      if (camp==0.) then
        camp1=1.
      else
        camp1=1./camp
      endif
!
!  set to zero and then rescale the testscalar
!  (in future, could call something like init_cc_simple)
!
      if (reinitialize_cctest) then
        do jtest=1,njtestscalar
          jcctest=icctest+(jtest-1)
          f(:,:,:,jcctest)=rescale_cctest(jtest)*f(:,:,:,jcctest)
        enddo
      endif
!
!  set lrescaling_testscalar=T if linit_cctest=T
!
      if (linit_cctest) then
        lrescaling_testscalar=.true.
      endif
!
!  Register an extra aux slot for uc if requested (so uc is written
!  to snapshots and can be easily analyzed later). For this to work you
!  must reserve enough auxiliary workspace by setting, for example,
!     ! MAUX CONTRIBUTION 9
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
!
!  After a reload, we need to rewrite index.pro, but the auxiliary
!  arrays are already allocated and must not be allocated again.
!
      if (lug_as_aux) then
        if (iug==0) then
          call farray_register_auxiliary('ug',iug,array=njtestscalar)
        else
          if (lroot) print*, 'initialize_testscalar: iug = ', iug
          call farray_index_append('iug',iug)
        endif
      endif
!
!  Write testscalar information to a file (for convenient post-processing).
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testscalar_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'zextent=',merge(1,0,zextent)
        write(1,'(a,i1)') 'lsoca_ug='  ,merge(1,0,lsoca_ug)
        write(1,'(3a)') "itestscalar='",trim(itestscalar)//"'"
        write(1,'(a,f5.2)') 'kxtestscalar=',kxtestscalar
        write(1,'(a,f5.2)') 'kytestscalar=',kytestscalar
        write(1,'(a,f5.2)') 'ktestscalar=',ktestscalar
        write(1,'(a,f7.4)') 'lam_testscalar=',lam_testscalar
        write(1,'(a,f7.4)') 'om_testscalar=', om_testscalar
        write(1,'(a,f7.4)') 'delta_testscalar=',delta_testscalar
        close(1)
      endif
!
    endsubroutine initialize_testscalar
!***********************************************************************
    subroutine init_cctest(f)
!
!  initialise testscalar; called from start.f90
!
!  26-nov-08/axel: adapted from testfield_z.f90
!
      use Cdata
      use Mpicomm, only: stop_it
      use InitialCondition, only: initial_condition_cctest
!
      real, dimension (mx,my,mz,mfarray), intent(OUT) :: f
      integer :: j
!
      do j=1,ninit
        select case (initcctest(j))

        case ('zero'); f(:,:,:,icctest:icctest+ntestscalar-1)=0.
        case ('nothing'); !(do nothing)

        case default
        !
        !  Catch unknown values
        !
          if (lroot) print*, 'init_cctest: check initcctest: ', trim(initcctest(j))
          call stop_it("")

        endselect
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_cctest(f)
!
    endsubroutine init_cctest
!***********************************************************************
    subroutine pencil_criteria_testscalar()
!
!  All pencils that the Testscalar module depends on are specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
!  Pencils requested for computing the rhs.
!
      lpenc_requested(i_uu)=.true.
      if (ltestscalar_per_unitvolume) lpenc_requested(i_divu)=.true.
!
!  Diagnostics pencils.
!
      lpenc_diagnos(i_divu)=.true.
!
    endsubroutine pencil_criteria_testscalar
!***********************************************************************
    subroutine pencil_interdep_testscalar(lpencil_in)
!
!  Interdependency among pencils from the Testscalar module is specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
      logical, dimension(npencils) :: lpencil_in
!
    endsubroutine pencil_interdep_testscalar
!***********************************************************************
    subroutine read_testscalar_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testscalar_init_pars, IOSTAT=iostat)
!
    endsubroutine read_testscalar_init_pars
!***********************************************************************
    subroutine write_testscalar_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testscalar_init_pars)
!
    endsubroutine write_testscalar_init_pars
!***********************************************************************
    subroutine read_testscalar_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testscalar_run_pars, IOSTAT=iostat)
!
    endsubroutine read_testscalar_run_pars
!***********************************************************************
    subroutine write_testscalar_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testscalar_run_pars)
!
    endsubroutine write_testscalar_run_pars
!***********************************************************************
    subroutine dcctest_dt(f,df,p)
!
!  testscalar evolution:
!
!  calculate dc^(p)/dt=-U.g^(p)-u.G^(p)-u.g+<u.g>]+kappa*del2c^(p),
!  where p=1,2 and g=grad c
!
!  26-nov-08/axel: adapted from testfield_z.f90
!  27-dec-08/axel: extended to x-dependent mean fields
!  17-apr-09/axel: added y column of kappa tensor
!  27-sep-13/MR  : changes due to ugtestm(mz,...  -->  ugtestm(nz,...,
!                                ugtestmy(my,...  -->  ugtestmy(ny,...
!
      use Cdata
      use Diagnostics
      use Sub
      use Hydro, only: uumz,lcalc_uumeanz
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx) :: cc,cctest,C0test,del2ctest,ug
      real, dimension (nx,3) :: ggtest, G0test=0.,uctest
      real, dimension (nx) :: ugtest,dctest,dugtest
      real, dimension (nx,3,njtestscalar) :: Fipq,Gipq,Hipq
      real, dimension (nx,njtestscalar) :: cpq
      real, dimension (nx,3) :: uufluct
      integer :: jcctest,jtest,j,nl,ml,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6
      logical,save :: ltest_ug=.false.
      real, dimension (nx) :: diffus_eta
!
      intent(in)    :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dcctest_dt: SOLVE'
      if (headtt) then
        if (icctest /= 0) call identify_bcs('cctest',icctest)
      endif
!
!  calculate uufluct=U-Umean
!
      if (lcalc_uumeanz) then
        do j=1,3
          uufluct(:,j)=p%uu(:,j)-uumz(n,j)
        enddo
      else
        uufluct=p%uu
      endif
!
!  multiply by exponential factor if lam_testscalar is different from zero
!  Keep camp1=1 for oscillatory fields.
!
      if (lam_testscalar/=0..or.om_testscalar/=0..or.delta_testscalar/=0.) then
        if (lam_testscalar/=0.) then
          tccinit_previous=tccinit-dccinit
          camp=exp(lam_testscalar*(t-tccinit_previous))
          camp1=1./camp
        endif
        if (om_testscalar/=0.) then
          camp=cos(om_testscalar*t)
          camp1=1.
        endif
        if (delta_testscalar/=0.) then
          if (t>=delta_testscalar_next.and.t<(delta_testscalar_next+dt)) then
            camp=1./(dt*delta_testscalar)
            delta_testscalar_next=t+delta_testscalar
          else
            camp=0.
          endif
          camp1=1.
        endif
      endif
!
!  do each of the 2+2+2 test scalars at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further down in the file.
!
      nl=n-n1+1; ml=m-m1+1
!
      do jtest=1,njtestscalar
        jcctest=icctest+(jtest-1)
        call del2(f,jcctest,del2ctest)
        call grad(f,jcctest,ggtest)
        select case (itestscalar)
          case ('G1-G2'); call set_ggtest_G1_G2(C0test,G0test,jtest)
          case ('G1-G2+const'); call set_ggtest_G1_G2_const(C0test,G0test,jtest)
          case ('G=0') !(don't do anything)
        case default
          call fatal_error('dcctest_dt','undefined itestscalar value')
        endselect
!
!  Prepare all the SOCA terms
!  Assume no mean flow (for now)
!
        call dot_mn(uufluct,G0test,ug)
        if (ltestscalar_per_unitvolume) ug=ug+p%divu*C0test
!
!  add SOCA terms
!
        df(l1:l2,m,n,jcctest)=df(l1:l2,m,n,jcctest)-ug+kappatest*del2ctest
!
!  compute and apply non-soca terms
!
        if (.not.lsoca_ug) then
!
!  use f-array for uc (if space has been allocated for this) and
!  if we don't test (i.e. if ltest_uc=.false.)
!
          if (iug/=0.and..not.ltest_ug) then
            ugtest=f(l1:l2,m,n,iug+(jtest-1))
          else
            call dot_mn(uufluct,ggtest,ugtest)
          endif
!
!  subtract average flux, unless we ignore the <ug> term (lignore_ugtestm=T)
!
          if (lignore_ugtestm) then
            dugtest=ugtest
          else
            if (jtest>=jtestz1 .and. jtest<=jtestz2) then
              dugtest=ugtest-ugtestm(nl,jtest)
            elseif (jtest>=jtestx1 .and. jtest<=jtestx2) then
              dugtest=ugtest-ugtestmx(:,jtest)
            elseif (jtest>=jtesty1 .and. jtest<=jtesty2) then
              dugtest=ugtest-ugtestmy(ml,jtest)
            endif
          endif
!
!  add to the right-hand side of the equation
!
          df(l1:l2,m,n,jcctest)=df(l1:l2,m,n,jcctest)-dugtest
        endif
!
!  Calculate gamma and kappa; begin by calculating dctest, uctest, and ugtest.
!
        if (ldiagnos.or.l1davgfirst) then
          cctest=f(l1:l2,m,n,jcctest)
          dctest=p%divu*cctest
          do j=1,3
            uctest(:,j)=uufluct(:,j)*cctest
          enddo
          call dot_mn(uufluct,ggtest,ugtest)
!
!  Prepare flux terms for diagnostics.
!
          cpq(:,jtest)=cctest
          Fipq(:,:,jtest)=uctest*camp1
          Gipq(:,3,jtest)=ugtest*camp1
          Hipq(:,3,jtest)=dctest*camp1
        endif
      enddo
!
!  Evaluate diffusive time step; just take the max of diffus_eta (if existent)
!  and whatever is calculated here.
!
      if (lfirst.and.ldt) then
        diffus_eta=kappatest*dxyz_2
        maxdiffus=max(maxdiffus,diffus_eta)
      endif
!
!  in the following block, we have already swapped the 4-6 entries with 7-9
!  The g95 compiler doesn't like to see an index that is out of bounds,
!  so prevent this warning by writing i3=3 and i4=4
!
      if (ldiagnos) then
!       if (idiag_F11x/=0) call yzsum_mn_name_x(Fipq(:,1,i1),idiag_F11x)
!       if (idiag_F21x/=0) call yzsum_mn_name_x(Fipq(:,2,i1),idiag_F21x)
!       if (idiag_F31x/=0) call yzsum_mn_name_x(Fipq(:,3,i1),idiag_F31x)
!       if (idiag_F12x/=0) call yzsum_mn_name_x(Fipq(:,1,i2),idiag_F12x)
!       if (idiag_F22x/=0) call yzsum_mn_name_x(Fipq(:,2,i2),idiag_F22x)
!       if (idiag_F32x/=0) call yzsum_mn_name_x(Fipq(:,3,i2),idiag_F32x)
        if (idiag_F11z/=0) call xysum_mn_name_z(Fipq(:,1,i1),idiag_F11z)
        if (idiag_F21z/=0) call xysum_mn_name_z(Fipq(:,2,i1),idiag_F21z)
        if (idiag_F31z/=0) call xysum_mn_name_z(Fipq(:,3,i1),idiag_F31z)
        if (idiag_F12z/=0) call xysum_mn_name_z(Fipq(:,1,i2),idiag_F12z)
        if (idiag_F22z/=0) call xysum_mn_name_z(Fipq(:,2,i2),idiag_F22z)
        if (idiag_F32z/=0) call xysum_mn_name_z(Fipq(:,3,i2),idiag_F32z)
!
!  Check whether njtestscalar is large enough.
!
        if (idiag_kap13/=0.or.idiag_kap23/=0.or.idiag_kap33/=0.or.&
            idiag_kap13z/=0.or.idiag_kap23z/=0.or.idiag_kap33z/=0) then
          if (njtestscalar<2) call stop_it('dcctest_dt: njtestscalar < 2 is insufficient')
        endif
        if (idiag_kap11/=0.or.idiag_kap21/=0.or.idiag_kap31/=0.or.&
            idiag_kap11z/=0.or.idiag_kap21z/=0.or.idiag_kap31z/=0) then
          if (njtestscalar<4) call stop_it('dcctest_dt: njtestscalar < 4 is insufficient')
        endif
        if (idiag_kap12/=0.or.idiag_kap22/=0.or.idiag_kap32/=0.or.&
            idiag_kap12z/=0.or.idiag_kap22z/=0.or.idiag_kap32z/=0) then
          if (njtestscalar<6) call stop_it('dcctest_dt: njtestscalar < 6 is insufficient')
        endif
!
!  First consider results from  z-dependent test fields.
!
        if (idiag_mgam33/=0) call sum_mn_name(  +cz(n)*Gipq(:,3,i1)+sz(n)*Gipq(:,3,i2) ,idiag_mgam33)
        if (idiag_mkap33/=0) call sum_mn_name(-(-sz(n)*Gipq(:,3,i1)+cz(n)*Gipq(:,3,i2))*ktestscalar1,idiag_mkap33)
        if (idiag_ngam33/=0) call sum_mn_name(  +cz(n)*Hipq(:,3,i1)+sz(n)*Hipq(:,3,i2) ,idiag_ngam33)
        if (idiag_nkap33/=0) call sum_mn_name(-(-sz(n)*Hipq(:,3,i1)+cz(n)*Hipq(:,3,i2))*ktestscalar1,idiag_nkap33)
!
!  Now do remaining kappa terms.
!
        if (idiag_kap11 /=0) call   sum_mn_name  (-(+cx(:)*Fipq(:,1,i3)+sx(:)*Fipq(:,1,i4)),idiag_kap11)
        if (idiag_kap11z/=0) call xysum_mn_name_z(-(+cx(:)*Fipq(:,1,i3)+sx(:)*Fipq(:,1,i4)),idiag_kap11z)
        if (idiag_kap21 /=0) call   sum_mn_name  (-(+cx(:)*Fipq(:,2,i3)+sx(:)*Fipq(:,2,i4)),idiag_kap21)
        if (idiag_kap21z/=0) call xysum_mn_name_z(-(+cx(:)*Fipq(:,2,i3)+sx(:)*Fipq(:,2,i4)),idiag_kap21z)
        if (idiag_kap31 /=0) call   sum_mn_name  (-(+cx(:)*Fipq(:,3,i3)+sx(:)*Fipq(:,3,i4)),idiag_kap31)
        if (idiag_kap31z/=0) call xysum_mn_name_z(-(+cx(:)*Fipq(:,3,i3)+sx(:)*Fipq(:,3,i4)),idiag_kap31z)
        if (idiag_kap12 /=0) call   sum_mn_name  (-(+cy(m)*Fipq(:,1,i5)+sy(m)*Fipq(:,1,i6)),idiag_kap12)
        if (idiag_kap12z/=0) call xysum_mn_name_z(-(+cy(m)*Fipq(:,1,i5)+sy(m)*Fipq(:,1,i6)),idiag_kap12z)
        if (idiag_kap22 /=0) call   sum_mn_name  (-(+cy(m)*Fipq(:,2,i5)+sy(m)*Fipq(:,2,i6)),idiag_kap22)
        if (idiag_kap22z/=0) call xysum_mn_name_z(-(+cy(m)*Fipq(:,2,i5)+sy(m)*Fipq(:,2,i6)),idiag_kap22z)
        if (idiag_kap32 /=0) call   sum_mn_name  (-(+cy(m)*Fipq(:,3,i5)+sy(m)*Fipq(:,3,i6)),idiag_kap32)
        if (idiag_kap32z/=0) call xysum_mn_name_z(-(+cy(m)*Fipq(:,3,i5)+sy(m)*Fipq(:,3,i6)),idiag_kap32z)
        if (idiag_kap13 /=0) call   sum_mn_name  (-(+cz(n)*Fipq(:,1,i1)+sz(n)*Fipq(:,1,i2)),idiag_kap13)
        if (idiag_kap13z/=0) call xysum_mn_name_z(-(+cz(n)*Fipq(:,1,i1)+sz(n)*Fipq(:,1,i2)),idiag_kap13z)
        if (idiag_kap23 /=0) call   sum_mn_name  (-(+cz(n)*Fipq(:,2,i1)+sz(n)*Fipq(:,2,i2)),idiag_kap23)
        if (idiag_kap23z/=0) call xysum_mn_name_z(-(+cz(n)*Fipq(:,2,i1)+sz(n)*Fipq(:,2,i2)),idiag_kap23z)
        if (idiag_kap33 /=0) call   sum_mn_name  (-(+cz(n)*Fipq(:,3,i1)+sz(n)*Fipq(:,3,i2)),idiag_kap33)
        if (idiag_kap33z/=0) call xysum_mn_name_z(-(+cz(n)*Fipq(:,3,i1)+sz(n)*Fipq(:,3,i2)),idiag_kap33z)
!
!  Finally do remaining gamma terms (pumping effect).
!
        if (idiag_gam11 /=0) call   sum_mn_name  (-(-sz(n)*Fipq(:,1,i3)+cz(n)*Fipq(:,1,i4))*ktestscalar1,idiag_gam11)
        if (idiag_gam11z/=0) call xysum_mn_name_z(-(-sz(n)*Fipq(:,1,i3)+cz(n)*Fipq(:,1,i4))*ktestscalar1,idiag_gam11z)
        if (idiag_gam21 /=0) call   sum_mn_name  (-(-sz(n)*Fipq(:,2,i3)+cz(n)*Fipq(:,2,i4))*ktestscalar1,idiag_gam21)
        if (idiag_gam21z/=0) call xysum_mn_name_z(-(-sz(n)*Fipq(:,2,i3)+cz(n)*Fipq(:,2,i4))*ktestscalar1,idiag_gam21z)
        if (idiag_gam31 /=0) call   sum_mn_name  (-(-sz(n)*Fipq(:,3,i3)+cz(n)*Fipq(:,3,i4))*ktestscalar1,idiag_gam31)
        if (idiag_gam31z/=0) call xysum_mn_name_z(-(-sz(n)*Fipq(:,3,i3)+cz(n)*Fipq(:,3,i4))*ktestscalar1,idiag_gam31z)
        if (idiag_gam12 /=0) call   sum_mn_name  (-(-sz(n)*Fipq(:,1,i5)+cz(n)*Fipq(:,1,i6))*ktestscalar1,idiag_gam12)
        if (idiag_gam12z/=0) call xysum_mn_name_z(-(-sz(n)*Fipq(:,1,i5)+cz(n)*Fipq(:,1,i6))*ktestscalar1,idiag_gam12z)
        if (idiag_gam22 /=0) call   sum_mn_name  (-(-sz(n)*Fipq(:,2,i5)+cz(n)*Fipq(:,2,i6))*ktestscalar1,idiag_gam22)
        if (idiag_gam22z/=0) call xysum_mn_name_z(-(-sz(n)*Fipq(:,2,i5)+cz(n)*Fipq(:,2,i6))*ktestscalar1,idiag_gam22z)
        if (idiag_gam32 /=0) call   sum_mn_name  (-(-sz(n)*Fipq(:,3,i5)+cz(n)*Fipq(:,3,i6))*ktestscalar1,idiag_gam32)
        if (idiag_gam32z/=0) call xysum_mn_name_z(-(-sz(n)*Fipq(:,3,i5)+cz(n)*Fipq(:,3,i6))*ktestscalar1,idiag_gam32z)
        if (idiag_gam13 /=0) call   sum_mn_name  (-(-sz(n)*Fipq(:,1,i1)+cz(n)*Fipq(:,1,i2))*ktestscalar1,idiag_gam13)
        if (idiag_gam13z/=0) call xysum_mn_name_z(-(-sz(n)*Fipq(:,1,i1)+cz(n)*Fipq(:,1,i2))*ktestscalar1,idiag_gam13z)
        if (idiag_gam23 /=0) call   sum_mn_name  (-(-sz(n)*Fipq(:,2,i1)+cz(n)*Fipq(:,2,i2))*ktestscalar1,idiag_gam23)
        if (idiag_gam23z/=0) call xysum_mn_name_z(-(-sz(n)*Fipq(:,2,i1)+cz(n)*Fipq(:,2,i2))*ktestscalar1,idiag_gam23z)
        if (idiag_gam33 /=0) call   sum_mn_name  (-(-sz(n)*Fipq(:,3,i1)+cz(n)*Fipq(:,3,i2))*ktestscalar,idiag_gam33)
        if (idiag_gam33z/=0) call xysum_mn_name_z(-(-sz(n)*Fipq(:,3,i1)+cz(n)*Fipq(:,3,i2))*ktestscalar,idiag_gam33z)
!
!  Extract values at one point.
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_c1pt/=0) call save_name(cpq(lpoint-nghost,i1),idiag_c1pt)
          if (idiag_c2pt/=0) call save_name(cpq(lpoint-nghost,i2),idiag_c2pt)
          if (idiag_c3pt/=0) call save_name(cpq(lpoint-nghost,i3),idiag_c3pt)
          if (idiag_c4pt/=0) call save_name(cpq(lpoint-nghost,i4),idiag_c4pt)
          if (idiag_c5pt/=0) call save_name(cpq(lpoint-nghost,i5),idiag_c5pt)
          if (idiag_c6pt/=0) call save_name(cpq(lpoint-nghost,i6),idiag_c6pt)
        endif
!
!  RMS values of small scales fields cpq in response to the test fields Bpq.
!
        if (idiag_c1rms/=0) call sum_mn_name(cpq(:,i1)**2,idiag_c1rms,lsqrt=.true.)
        if (idiag_c2rms/=0) call sum_mn_name(cpq(:,i2)**2,idiag_c2rms,lsqrt=.true.)
        if (idiag_c3rms/=0) call sum_mn_name(cpq(:,i3)**2,idiag_c3rms,lsqrt=.true.)
        if (idiag_c4rms/=0) call sum_mn_name(cpq(:,i4)**2,idiag_c4rms,lsqrt=.true.)
        if (idiag_c5rms/=0) call sum_mn_name(cpq(:,i5)**2,idiag_c5rms,lsqrt=.true.)
        if (idiag_c6rms/=0) call sum_mn_name(cpq(:,i6)**2,idiag_c6rms,lsqrt=.true.)
!
      endif
!
    endsubroutine dcctest_dt
!***********************************************************************
    subroutine get_slices_testscalar(f,slices)
! 
!  Write slices for animation.
! 
!  26-nov-08/axel: adapted from testfield_z.f90
! 
      use Cdata, only: icctest
      use Slices_methods, only: assign_slices_vec
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
! 
!  Loop over slices
! 
      select case (trim(slices%name))
!
!  Testfield slice
!
        case ('cctest')
          call assign_slices_vec(slices,f,icctest,6)

      endselect
!
    endsubroutine get_slices_testscalar
!***********************************************************************
    subroutine testscalar_after_boundary(f)
!
!  Calculate either <ug> or, if ltestscalar_per_unitvolume=T, <ug>+<du>,
!  where g=gradc and d=divu. This is needed when lsoca_ug=F.
!
!  26-nov-08/axel: adapted from testfield_z.f90
!  27-dec-08/axel: also calculate yz-averages
!   8-aug-10/axel: renamed calc_ltestscalar_pars => testscalar_after_boundary
!  27-sep-13/MR  : changes due to ugtestm(mz,...  -->  ugtestm(nz,... ,
!                                ugtestmy(my,...  -->  ugtestmy(ny,... ;
!                  restricted calculation of hydro pencils; simplified communication.
!
      use Cdata
      use Sub
      use Hydro, only: calc_pencils_hydro
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (nx) :: cctest,ugtest
      real, dimension (nx,3) :: ggtest
      integer :: jcctest,jtest,nl,ml
      logical :: headtt_save
      real :: fac_xy,fac_yz,fac_xz
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac_xy=1./nxygrid
      fac_yz=1./nyzgrid
      fac_xz=1./nxzgrid
!
!  Do each of the 2+2 test fields at a time,
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
!  ** Start with z-dependent mean fields **
!
      lpenc_loc=.false.; lpenc_loc(i_uu)=.true.
      if (ltestscalar_per_unitvolume) lpenc_loc(i_divu)=.true.
!
      do jtest=jtestz1,jtestz2
        jcctest=icctest+(jtest-1)
        if (lsoca_ug) then
          ugtestm(:,jtest)=0.
        else
          do n=n1,n2
            nl=n-n1+1
            ugtestm(nl,jtest)=0.
            do m=m1,m2
              cctest=f(l1:l2,m,n,jcctest)
              call calc_pencils_hydro(f,p,lpenc_loc)
              call grad(f,jcctest,ggtest)
              call dot_mn(p%uu,ggtest,ugtest)
!
!  Add divu*c term if ltestscalar_per_unitvolume=T.
!
              if (ltestscalar_per_unitvolume) ugtest=ugtest+p%divu*cctest
!
!  Put <ug> or <ug-dc> into auxiliary array and compute its average.
!
              if (iug/=0) f(l1:l2,m,n,iug+(jtest-1))=ugtest
              ugtestm(nl,jtest)=ugtestm(nl,jtest)+fac_xy*sum(ugtest)
              headtt=.false.
            enddo
          enddo
        endif
      enddo
!
!  ** Now do x-dependent mean fields **
!
      if (njtestscalar>=jtestx2) then
        do jtest=jtestx1,jtestx2
          jcctest=icctest+(jtest-1)
          if (lsoca_ug) then
            ugtestmx(:,jtest)=0.
          else
            ugtestmx(:,jtest)=0.
            do n=n1,n2
              do m=m1,m2
                cctest=f(l1:l2,m,n,jcctest)
                call calc_pencils_hydro(f,p,lpenc_loc)
                call grad(f,jcctest,ggtest)
                call dot_mn(p%uu,ggtest,ugtest)
!
!  Add divu*c term if ltestscalar_per_unitvolume=T.
!
              if (ltestscalar_per_unitvolume) ugtest=ugtest+p%divu*cctest
!
!  Put <ug> or <ug-dc> into auxiliary array and compute its average.
!
                if (iug/=0) f(l1:l2,m,n,iug+(jtest-1))=ugtest
                ugtestmx(:,jtest)=ugtestmx(:,jtest)+fac_yz*ugtest
                headtt=.false.
              enddo
            enddo
          endif
        enddo
      endif
!
!  ** Finally do y-dependent mean fields **
!
      if (njtestscalar>=jtesty2) then
        do jtest=jtesty1,jtesty2
          jcctest=icctest+(jtest-1)
          if (lsoca_ug) then
            ugtestmy(:,jtest)=0.
          else
            ugtestmy(:,jtest)=0.
            do m=m1,m2
              ml=m-m1+1
              do n=n1,n2
                cctest=f(l1:l2,m,n,jcctest)
                call calc_pencils_hydro(f,p,lpenc_loc)
                call grad(f,jcctest,ggtest)
                call dot_mn(p%uu,ggtest,ugtest)
!
!  Add divu*c term if ltestscalar_per_unitvolume=T.
!
              if (ltestscalar_per_unitvolume) ugtest=ugtest+p%divu*cctest
!
!  Put <ug> or <ug-dc> into auxiliary array and compute its average.
!
                if (iug/=0) f(l1:l2,m,n,iug+(jtest-1))=ugtest
                ugtestmy(ml,jtest)=ugtestmy(ml,jtest)+fac_xz*sum(ugtest)
                headtt=.false.
              enddo
            enddo
          endif
        enddo
      endif
!
!  do communication
!  Do this first for z-dependent mean fields
!
      call finalize_aver(nprocxy,12,ugtestm)
!
!  Next, do this for x-dependent mean fields
!
      call finalize_aver(nprocyz,23,ugtestmx)
!
!  Finally, do this for y-dependent mean fields
!
      call finalize_aver(nprocxz,13,ugtestmy)
!
!  reset headtt
!
      headtt=headtt_save
!
    endsubroutine testscalar_after_boundary
!***********************************************************************
    subroutine rescaling_testscalar(f)
!
!  Rescale testscalar by factor rescale_cctest(jtest),
!  which could be different for different testscalar
!
!  26-nov-08/axel: adapted from testfield_z.f90
!
      use Cdata
      use Sub, only: read_snaptime, update_snaptime
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=fnlen) :: file
      logical :: ltestscalar_out
      integer,save :: ifirst=0
      integer :: j,jtest
!
      intent(inout) :: f
!
! reinitialize cctest periodically if requested
!
      if (linit_cctest) then
        file=trim(datadir)//'/tinit_cctest.dat'
        if (ifirst==0) then
          call read_snaptime(trim(file),tccinit,nccinit,dccinit,t)
          if (tccinit==0 .or. tccinit < t-dccinit) then
            tccinit=t+dccinit
          endif
          ifirst=1
        endif
!
!  Do only one xy plane at a time (for cache efficiency)
!
        if (t >= tccinit) then
          do jtest=1,njtestscalar
            j=icctest+(jtest-1)
            f(l1:l2,m1:m2,n1:n2,j)=rescale_cctest(jtest)*f(l1:l2,m1:m2,n1:n2,j)
          enddo
          call update_snaptime(file,tccinit,nccinit,dccinit,t,ltestscalar_out)
        endif
      endif
!
    endsubroutine rescaling_testscalar
!***********************************************************************
    subroutine set_ggtest_G1_G2(C0test,G0test,jtest)
!
!  set testscalar
!
!  26-nov-08/axel: adapted from testfield_z.f90
!  27-dec-08/axel: extended to x-dependent mean fields
!
      use Cdata
!
      real, dimension (nx,3) :: G0test
      real, dimension (nx) :: C0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: C0test,G0test
!
!  set G0test for each of the 2+2 cases
!
      select case (jtest)
      case (1); C0test=+camp*ktestscalar1*sz(n); G0test(:,1)=0.; G0test(:,2)=0.; G0test(:,3)=camp*cz(n)
      case (2); C0test=-camp*ktestscalar1*cz(n); G0test(:,1)=0.; G0test(:,2)=0.; G0test(:,3)=camp*sz(n)
      case (3); C0test=+camp*ktestscalar1*sx(:); G0test(:,1)=camp*cx(:); G0test(:,2)=0.; G0test(:,3)=0.
      case (4); C0test=-camp*ktestscalar1*cx(:); G0test(:,1)=camp*sx(:); G0test(:,2)=0.; G0test(:,3)=0.
      case (5); C0test=+camp*ktestscalar1*sy(m); G0test(:,1)=0.; G0test(:,2)=camp*cy(m); G0test(:,3)=0.
      case (6); C0test=-camp*ktestscalar1*cy(m); G0test(:,1)=0.; G0test(:,2)=camp*sy(m); G0test(:,3)=0.
      case default; C0test(:)=0.; G0test(:,:)=0.
      endselect
!
    endsubroutine set_ggtest_G1_G2
!***********************************************************************
    subroutine set_ggtest_G1_G2_const(C0test,G0test,jtest)
!
!  set testscalar
!
!  26-nov-08/axel: adapted from testfield_z.f90
!  27-dec-08/axel: extended to x-dependent mean fields
!
      use Cdata
!
      real, dimension (nx,3) :: G0test
      real, dimension (nx) :: C0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: C0test,G0test
!
!  set G0test for each of the 2+2 cases
!
      select case (jtest)
      case (1); G0test(:,1)=0.; G0test(:,2)=0.; G0test(:,3)=camp*cz(n)
      case (2); G0test(:,1)=0.; G0test(:,2)=0.; G0test(:,3)=camp*sz(n)
      case (3); G0test(:,1)=0.; G0test(:,2)=0.; G0test(:,3)=camp
      case (4); G0test(:,1)=camp*cx(:); G0test(:,2)=0.; G0test(:,3)=0.
      case (5); G0test(:,1)=camp*sx(:); G0test(:,2)=0.; G0test(:,3)=0.
      case (6); G0test(:,1)=camp      ; G0test(:,2)=0.; G0test(:,3)=0.
      case (7); G0test(:,1)=0.; G0test(:,2)=camp*cy(m); G0test(:,3)=0.
      case (8); G0test(:,1)=0.; G0test(:,2)=camp*sy(m); G0test(:,3)=0.
      case (9); G0test(:,1)=0.; G0test(:,2)=camp      ; G0test(:,3)=0.
      case default; C0test(:)=0.; G0test(:,:)=0.
      endselect
!
    endsubroutine set_ggtest_G1_G2_const
!***********************************************************************
    subroutine rprint_testscalar(lreset,lwrite)
!
!  reads and registers print parameters relevant for testscalar fields
!
!  26-nov-08/axel: adapted from testfield_z.f90
!  27-dec-08/axel: included kap11, kap21, and kap31
!
      use Cdata
      use Diagnostics
      use FArrayManager, only: farray_index_append
      use General, only: loptest
!
      integer :: iname,inamez
      logical :: lreset
      logical, optional :: lwrite
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_F11z=0; idiag_F21z=0; idiag_F31z=0
        idiag_F12z=0; idiag_F22z=0; idiag_F32z=0
        idiag_kap11=0; idiag_kap21=0; idiag_kap31=0
        idiag_kap12=0; idiag_kap22=0; idiag_kap32=0
        idiag_kap13=0; idiag_kap23=0; idiag_kap33=0
        idiag_kap11z=0; idiag_kap21z=0; idiag_kap31z=0
        idiag_kap12z=0; idiag_kap22z=0; idiag_kap32z=0
        idiag_kap13z=0; idiag_kap23z=0; idiag_kap33z=0
        idiag_mkap33=0; idiag_nkap33=0
        idiag_gam11=0; idiag_gam21=0; idiag_gam31=0
        idiag_gam12=0; idiag_gam22=0; idiag_gam32=0
        idiag_gam13=0; idiag_gam23=0; idiag_gam33=0
        idiag_gam11z=0; idiag_gam21z=0; idiag_gam31z=0
        idiag_gam12z=0; idiag_gam22z=0; idiag_gam32z=0
        idiag_gam13z=0; idiag_gam23z=0; idiag_gam33z=0
        idiag_c1rms=0; idiag_c2rms=0
        idiag_c1pt=0; idiag_c2pt=0
      endif
!
!  check for those quantities that we want to evaluate online
! 
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'kap11',idiag_kap11)
        call parse_name(iname,cname(iname),cform(iname),'kap21',idiag_kap21)
        call parse_name(iname,cname(iname),cform(iname),'kap31',idiag_kap31)
        call parse_name(iname,cname(iname),cform(iname),'kap12',idiag_kap12)
        call parse_name(iname,cname(iname),cform(iname),'kap22',idiag_kap22)
        call parse_name(iname,cname(iname),cform(iname),'kap32',idiag_kap32)
        call parse_name(iname,cname(iname),cform(iname),'kap13',idiag_kap13)
        call parse_name(iname,cname(iname),cform(iname),'kap23',idiag_kap23)
        call parse_name(iname,cname(iname),cform(iname),'kap33',idiag_kap33)
        call parse_name(iname,cname(iname),cform(iname),'mgam33',idiag_mgam33)
        call parse_name(iname,cname(iname),cform(iname),'mkap33',idiag_mkap33)
        call parse_name(iname,cname(iname),cform(iname),'ngam33',idiag_ngam33)
        call parse_name(iname,cname(iname),cform(iname),'nkap33',idiag_nkap33)
        call parse_name(iname,cname(iname),cform(iname),'gam11',idiag_gam11)
        call parse_name(iname,cname(iname),cform(iname),'gam21',idiag_gam21)
        call parse_name(iname,cname(iname),cform(iname),'gam31',idiag_gam31)
        call parse_name(iname,cname(iname),cform(iname),'gam12',idiag_gam12)
        call parse_name(iname,cname(iname),cform(iname),'gam22',idiag_gam22)
        call parse_name(iname,cname(iname),cform(iname),'gam32',idiag_gam32)
        call parse_name(iname,cname(iname),cform(iname),'gam13',idiag_gam13)
        call parse_name(iname,cname(iname),cform(iname),'gam23',idiag_gam23)
        call parse_name(iname,cname(iname),cform(iname),'gam33',idiag_gam33)
        call parse_name(iname,cname(iname),cform(iname),'c1rms',idiag_c1rms)
        call parse_name(iname,cname(iname),cform(iname),'c2rms',idiag_c2rms)
        call parse_name(iname,cname(iname),cform(iname),'c3rms',idiag_c3rms)
        call parse_name(iname,cname(iname),cform(iname),'c4rms',idiag_c4rms)
        call parse_name(iname,cname(iname),cform(iname),'c5rms',idiag_c5rms)
        call parse_name(iname,cname(iname),cform(iname),'c6rms',idiag_c6rms)
        call parse_name(iname,cname(iname),cform(iname),'c1pt',idiag_c1pt)
        call parse_name(iname,cname(iname),cform(iname),'c2pt',idiag_c2pt)
        call parse_name(iname,cname(iname),cform(iname),'c3pt',idiag_c3pt)
        call parse_name(iname,cname(iname),cform(iname),'c4pt',idiag_c4pt)
        call parse_name(iname,cname(iname),cform(iname),'c5pt',idiag_c5pt)
        call parse_name(iname,cname(iname),cform(iname),'c6pt',idiag_c6pt)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F11z',idiag_F11z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F21z',idiag_F21z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F31z',idiag_F31z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F12z',idiag_F12z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F22z',idiag_F22z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F32z',idiag_F32z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kap11z',idiag_kap11z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kap21z',idiag_kap21z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kap31z',idiag_kap31z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kap12z',idiag_kap12z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kap22z',idiag_kap22z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kap32z',idiag_kap32z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kap13z',idiag_kap13z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kap23z',idiag_kap23z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kap33z',idiag_kap33z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gam11z',idiag_gam11z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gam21z',idiag_gam21z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gam31z',idiag_gam31z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gam12z',idiag_gam12z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gam22z',idiag_gam22z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gam32z',idiag_gam32z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gam13z',idiag_gam13z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gam23z',idiag_gam23z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gam33z',idiag_gam33z)
      enddo
!
!  check for those quantities for which we want video slices
!
      where(cnamev=='cctest') cformv='DEFINED'
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (loptest(lwrite)) then
        call farray_index_append('icctest',icctest)
        call farray_index_append('ntestscalar',ntestscalar)
      endif
!
    endsubroutine rprint_testscalar
!***********************************************************************
endmodule Testscalar
