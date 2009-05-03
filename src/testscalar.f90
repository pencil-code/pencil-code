! $Id$

!  This modules deals with all aspects of testscalar fields; if no
!  testscalar fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testscalar relevant subroutines listed in here.

!  Note: this routine requires that MVAR and MAUX contributions
!  together with njtest are set correctly in the cparam.local file.
!  njtest must be set at the end of the file such that njtest=MVAR.
!
!  Example:
!  ! MVAR CONTRIBUTION 2
!  ! MAUX CONTRIBUTION 2
!  integer, parameter :: njtest=2

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

  use Cparam
  use Messages

  implicit none

  include 'testscalar.h'
!
! Slice precalculation buffers
!
  real, target, dimension (nx,ny,1) :: cc1_xy
  real, target, dimension (nx,ny,1) :: cc1_xy2
  real, target, dimension (nx,nz,1) :: cc1_xz
  real, target, dimension (ny,nz,1) :: cc1_yz
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
  real :: lam_testscalar=0.,om_testscalar=0.,delta_testscalar=0.
  real :: delta_testscalar_next=0.
  integer, parameter :: mtestscalar=njtest
  integer :: jtestz1=1,jtestz2=2,jtestx1=3,jtestx2=4,jtesty1=5,jtesty2=6
  integer :: nccinit
  real :: camp=1.,camp1=1.
  namelist /testscalar_init_pars/ &
       zextent,initcctest, &
       amplcctest,kx_cctest,ky_cctest,kz_cctest, &
       phasex_cctest,phasez_cctest, &
       lug_as_aux

  ! run parameters
  real :: kappatest=0.,kappatest1=0.
  real, dimension(njtest) :: rescale_cctest=0.
  logical :: ltestscalar_newz=.true.,leta_rank2=.false.
  logical :: ltestscalar_newx=.false.
  logical :: ltestscalar_newy=.false.
  namelist /testscalar_run_pars/ &
       reinitialize_cctest,zextent,lsoca_ug, &
       lset_cctest2,kappatest,kappatest1,itestscalar,ktestscalar, &
       lam_testscalar,om_testscalar,delta_testscalar, &
       ltestscalar_newx,ltestscalar_newz, &
       lug_as_aux,lignore_ugtestm, &
       dccinit,linit_cctest,camp, &
       rescale_cctest

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_kap11=0      ! DIAG_DOC: $\kappa_{11}$
  integer :: idiag_kap21=0      ! DIAG_DOC: $\kappa_{21}$
  integer :: idiag_kap31=0      ! DIAG_DOC: $\kappa_{31}$
  integer :: idiag_kap12=0      ! DIAG_DOC: $\kappa_{12}$
  integer :: idiag_kap22=0      ! DIAG_DOC: $\kappa_{22}$
  integer :: idiag_kap32=0      ! DIAG_DOC: $\kappa_{32}$
  integer :: idiag_kap13=0      ! DIAG_DOC: $\kappa_{13}$
  integer :: idiag_kap23=0      ! DIAG_DOC: $\kappa_{23}$
  integer :: idiag_kap33=0      ! DIAG_DOC: $\kappa_{33}$
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
  real, dimension (mz,mtestscalar) :: ugtestm
  real, dimension (nx,mtestscalar) :: ugtestmx
  real, dimension (my,mtestscalar) :: ugtestmy

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
      use Mpicomm
      use Sub
!
      integer :: j
!
!  Set first and last index of text field
!  Note: iaxtest, iaytest, and iaztest are initialized to the first test field.
!  These values are used in this form in start, but later overwritten.
!
      icctest=nvar+1
      icctestpq=icctest+(njtest-1)
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
      do j=1,ntestscalar
        varname(j) = 'cctest'
      enddo
!
!  identify version number
!
      if (lroot) call cvs_id( &
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
      real :: ktestscalar_effective
      integer :: jtest, jcctest
!
!  Precalculate kappatest if 1/kappatest (==kappatest1) is given instead
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
        ktestscalar_effective=ktestscalar*(2.*pi/Lx)
        xtestscalar=ktestscalar_effective*(x(l1:l2)-x0)-pi
      else
        ktestscalar_effective=ktestscalar
        xtestscalar=x(l1:l2)
      endif
      cx=cos(ktestscalar*xtestscalar)
      sx=sin(ktestscalar*xtestscalar)
!
!  Choice of using rescaled y-array or original y-array
!  Define ktestscalar_effective to deal with boxes bigger than 2pi.
!
      if (ltestscalar_newy) then
        ktestscalar_effective=ktestscalar*(2.*pi/Ly)
        ytestscalar=ktestscalar_effective*(y-y0)-pi
      else
        ktestscalar_effective=ktestscalar
        ytestscalar=y
      endif
      cy=cos(ktestscalar*ytestscalar)
      sy=sin(ktestscalar*ytestscalar)
!
!  Choice of using rescaled z-array or original z-array
!  Define ktestscalar_effective to deal with boxes bigger than 2pi.
!
      if (ltestscalar_newz) then
        ktestscalar_effective=ktestscalar*(2.*pi/Lz)
        ztestscalar=ktestscalar_effective*(z-z0)-pi
      else
        ktestscalar_effective=ktestscalar
        ztestscalar=z
      endif
      cz=cos(ktestscalar*ztestscalar)
      sz=sin(ktestscalar*ztestscalar)
!
!  Optionally, one can determine the phase in the actual field
!  and modify the following calculations in calc_ltestscalar_pars.
!
!  debug output
!
      if (lroot.and.ip<14) then
        print*,'cz=',cz
        print*,'sz=',sz
      endif
!
!  Also calculate its inverse, but only if different from zero
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
        do jtest=1,njtest
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
          call farray_register_auxiliary('ug',iug,vector=njtest)
        endif
        if (iug/=0.and.lroot) then
          print*, 'initialize_magnetic: iug = ', iug
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'iug=',iug
          close(3)
        endif
      endif
!
!  write testscalar information to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testscalar_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'zextent=',merge(1,0,zextent)
        write(1,'(a,i1)') 'lsoca_ug='  ,merge(1,0,lsoca_ug)
        write(1,'(3a)') "itestscalar='",trim(itestscalar)//"'"
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
      use Mpicomm
      use Initcond
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: cc
      real, dimension (nx) :: b2,fact
      real :: beq2
      integer :: j
!
      do j=1,ninit

      select case(initcctest(j))

      case('zero'); f(:,:,:,icctest:icctest+ntestscalar-1)=0.
      case('nothing'); !(do nothing)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'init_cctest: check initcctest: ', trim(initcctest(j))
        call stop_it("")

      endselect
      enddo
!
    endsubroutine init_cctest
!***********************************************************************
    subroutine pencil_criteria_testscalar()
!
!   All pencils that the Testscalar module depends on are specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
      lpenc_requested(i_uu)=.true.
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
    subroutine read_testscalar_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=testscalar_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testscalar_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_testscalar_init_pars
!***********************************************************************
    subroutine write_testscalar_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=testscalar_init_pars)

    endsubroutine write_testscalar_init_pars
!***********************************************************************
    subroutine read_testscalar_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=testscalar_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testscalar_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_testscalar_run_pars
!***********************************************************************
    subroutine write_testscalar_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=testscalar_run_pars)

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
!
      use Cdata
      use Diagnostics
      use Sub
      use Hydro, only: uumz,lcalc_uumean
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx) :: cc,cctest,del2ctest,ug
      real, dimension (nx,3) :: ggtest, G0test=0.,uctest
      real, dimension (nx) :: ugtest,dugtest
      real, dimension (nx,3,njtest) :: Fipq
      real, dimension (nx,njtest) :: cpq
      real, dimension (nx,3) :: uufluct
      integer :: jcctest,jtest,jfnamez,j,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6
      logical,save :: ltest_ug=.false.
!
      intent(in)     :: f,p
      intent(inout)  :: df
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
      if (lcalc_uumean) then
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
          if (t.ge.delta_testscalar_next.and.t.lt.(delta_testscalar_next+dt)) then
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
      do jtest=1,njtest
        jcctest=icctest+(jtest-1)
        call del2(f,jcctest,del2ctest)
        call grad(f,jcctest,ggtest)
        select case(itestscalar)
          case('G1-G2'); call set_ggtest_G1_G2(G0test,jtest)
          case('G=0') !(don't do anything)
        case default
          call fatal_error('dcctest_dt','undefined itestscalar value')
        endselect
!
        call dot_mn(uufluct,G0test,ug)
        if (lsoca_ug) then
          df(l1:l2,m,n,jcctest)=df(l1:l2,m,n,jcctest) &
            +ug+kappatest*del2ctest
        else
!
!  use f-array for uc (if space has been allocated for this) and
!  if we don't test (i.e. if ltest_uc=.false.)
!
          if (iug/=0.and..not.ltest_ug) then
            ugtest=f(l1:l2,m,n,iug+(jtest-1))
          else
            call dot_mn(p%uu,ggtest,ugtest)
          endif
!
!  subtract average flux, unless we ignore the <ug> term (lignore_ugtestm=T)
!
          if (lignore_ugtestm) then
            dugtest=ugtest
          else
            if (jtest>=jtestz1 .and. jtest<=jtestz2) then
              dugtest=ugtest-ugtestm(n,jtest)
            elseif (jtest>=jtestx1 .and. jtest<=jtestx2) then
              dugtest=ugtest-ugtestmx(:,jtest)
            elseif (jtest>=jtesty1 .and. jtest<=jtesty2) then
              dugtest=ugtest-ugtestmy(m,jtest)
            endif
          endif
!
!  advance test field equation
!
          df(l1:l2,m,n,jcctest)=df(l1:l2,m,n,jcctest) &
            +ug+kappatest*del2ctest+dugtest
        endif
!
!  calculate kappa, begin by calculating uctest (not ugtest!)
!
        if (ldiagnos.or.l1ddiagnos) then
          cctest=f(l1:l2,m,n,jcctest)
          do j=1,3
            uctest(:,j)=p%uu(:,j)*cctest
          enddo
        endif
        cpq(:,jtest)=cctest
        Fipq(:,:,jtest)=uctest*camp1
      enddo
!
!  diffusive time step, just take the max of diffus_eta (if existent)
!  and whatever is calculated here
!
      if (lfirst.and.ldt) then
        diffus_eta=max(diffus_eta,kappatest*dxyz_2)
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
!  averages of kappa
!
        if (idiag_kap11/=0) call sum_mn_name(+cx(:)*Fipq(:,1,i3)+sx(:)*Fipq(:,1,i4),idiag_kap11)
        if (idiag_kap21/=0) call sum_mn_name(+cx(:)*Fipq(:,2,i3)+sx(:)*Fipq(:,2,i4),idiag_kap21)
        if (idiag_kap31/=0) call sum_mn_name(+cx(:)*Fipq(:,3,i3)+sx(:)*Fipq(:,3,i4),idiag_kap31)
        if (idiag_kap12/=0) call sum_mn_name(+cy(m)*Fipq(:,1,i5)+sy(m)*Fipq(:,1,i6),idiag_kap12)
        if (idiag_kap22/=0) call sum_mn_name(+cy(m)*Fipq(:,2,i5)+sy(m)*Fipq(:,2,i6),idiag_kap22)
        if (idiag_kap32/=0) call sum_mn_name(+cy(m)*Fipq(:,3,i5)+sy(m)*Fipq(:,3,i6),idiag_kap32)
        if (idiag_kap13/=0) call sum_mn_name(+cz(n)*Fipq(:,1,i1)+sz(n)*Fipq(:,1,i2),idiag_kap13)
        if (idiag_kap23/=0) call sum_mn_name(+cz(n)*Fipq(:,2,i1)+sz(n)*Fipq(:,2,i2),idiag_kap23)
        if (idiag_kap33/=0) call sum_mn_name(+cz(n)*Fipq(:,3,i1)+sz(n)*Fipq(:,3,i2),idiag_kap33)
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
!  rms values of small scales fields cpq in response to the test fields Bpq
!  Obviously idiag_b0rms and idiag_b12rms cannot both be invoked!
!  Needs modification!
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
!  write B-slices for output in wvid in run.f90
!  Note: ix is the index with respect to array with ghost zones.
! 
      if (lvideo.and.lfirst) then
        cc1_yz(m-m1+1,n-n1+1,1)=cpq(ix_loc-l1+1,1)
        if (m==iy_loc)  cc1_xz(:,n-n1+1,1)=cpq(:,1)
        if (n==iz_loc)  cc1_xy(:,m-m1+1,1)=cpq(:,1)
        if (n==iz2_loc) cc1_xy2(:,m-m1+1,1)=cpq(:,1)
      endif
!
    endsubroutine dcctest_dt
!***********************************************************************
    subroutine get_slices_testscalar(f,slices)
! 
!  Write slices for animation of magnetic variables.
! 
!  26-nov-08/axel: adapted from testfield_z.f90
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
        case ('cc1')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =>cc1_yz(:,:,slices%index)
            slices%xz =>cc1_xz(:,:,slices%index)
            slices%xy =>cc1_xy(:,:,slices%index)
            slices%xy2=>cc1_xy2(:,:,slices%index)
            if (slices%index<3) slices%ready=.true.
          endif
      endselect
!
    endsubroutine get_slices_testscalar
!***********************************************************************
    subroutine calc_ltestscalar_pars(f)
!
!  calculate <ug>, which is needed when lsoca_ug=.false.
!
!  26-nov-08/axel: adapted from testfield_z.f90
!  27-dec-08/axel: also calculate yz-averages
!
      use Cdata
      use Sub
      use Hydro, only: calc_pencils_hydro
      use Mpicomm, only: mpireduce_sum, mpibcast_real, mpibcast_real_arr
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: c,s
!
      real, dimension (nx,nprocy,nprocz,njtest) :: ugtestmx1=0.,ugtestmx1_tmp=0.
      real, dimension (ny,nprocy,njtest) :: ugtestmy1=0.,ugtestmy1_tmp=0.
      real, dimension (nz,nprocz,njtest) :: ugtestm1=0.,ugtestm1_tmp=0.
!
      real, dimension (nx) :: cctest,ugtest
      real, dimension (nx,3) :: ggtest
      integer :: jcctest,jtest,j,jug,jpy,jpz
      integer :: nxy=nxgrid*nygrid,nyz=nygrid*nzgrid,nxz=nxgrid*nzgrid
      logical :: headtt_save
      real :: fac_xy,fac_yz,fac_xz
      type (pencil_case) :: p
!
      intent(inout) :: f
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac_xy=1./nxy
      fac_yz=1./nyz
      fac_xz=1./nxz
!
!  do each of the 2+2 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!  Start with z-dependent mean fields
!
      do jtest=jtestz1,jtestz2
        jcctest=icctest+(jtest-1)
        if (lsoca_ug) then
          ugtestm(:,jtest)=0.
        else
          do n=n1,n2
            ugtestm(n,jtest)=0.
            do m=m1,m2
              cctest=f(l1:l2,m,n,jcctest)
              call calc_pencils_hydro(f,p)
              call grad(f,jcctest,ggtest)
              call dot_mn(p%uu,ggtest,ugtest)
              jug=iug+(jtest-1)
              if (iug/=0) f(l1:l2,m,n,jug)=ugtest
              ugtestm(n,jtest)=ugtestm(n,jtest)+fac_xy*sum(ugtest)
              headtt=.false.
            enddo
            ugtestm1(n-n1+1,ipz+1,jtest)=ugtestm(n,jtest)
          enddo
        endif
      enddo
!
!  Now do x-dependent mean fields
!
      do jtest=jtestx1,jtestx2
        jcctest=icctest+(jtest-1)
        if (lsoca_ug) then
          ugtestmx(:,jtest)=0.
        else
          ugtestmx(:,jtest)=0.
          do n=n1,n2
            do m=m1,m2
              cctest=f(l1:l2,m,n,jcctest)
              call calc_pencils_hydro(f,p)
              call grad(f,jcctest,ggtest)
              call dot_mn(p%uu,ggtest,ugtest)
              jug=iug+(jtest-1)
              if (iug/=0) f(l1:l2,m,n,jug)=ugtest
              ugtestmx(:,jtest)=ugtestmx(:,jtest)+fac_yz*ugtest
              headtt=.false.
            enddo
          enddo
          ugtestmx1(:,ipy+1,ipz+1,jtest)=ugtestmx(:,jtest)
        endif
      enddo
!
!  Finally do y-dependent mean fields
!
      do jtest=jtesty1,jtesty2
        jcctest=icctest+(jtest-1)
        if (lsoca_ug) then
          ugtestmy(:,jtest)=0.
        else
          ugtestmy(:,jtest)=0.
          do m=m1,m2
            do n=n1,n2
              cctest=f(l1:l2,m,n,jcctest)
              call calc_pencils_hydro(f,p)
              call grad(f,jcctest,ggtest)
              call dot_mn(p%uu,ggtest,ugtest)
              jug=iug+(jtest-1)
              if (iug/=0) f(l1:l2,m,n,jug)=ugtest
              ugtestmy(m,jtest)=ugtestmy(m,jtest)+fac_xz*sum(ugtest)
              headtt=.false.
            enddo
            ugtestmy1(m-m1+1,ipy+1,jtest)=ugtestmy(m,jtest)
          enddo
        endif
      enddo
!
!  do communication for array of size nz*nprocz*3*njtest
!  Do this first for z-dependent mean fields
!
      if (nprocy>1) then
        call mpireduce_sum(ugtestm1,ugtestm1_tmp,nz*nprocz*njtest)
        call mpibcast_real_arr(ugtestm1_tmp,nz*nprocz*njtest)
        do jtest=1,njtest
          do n=n1,n2
            ugtestm(n,jtest)=ugtestm1_tmp(n-n1+1,ipz+1,jtest)
          enddo
        enddo
      endif
!
!  Next, do this for x-dependent mean fields
!
      if (ncpus>1) then
        call mpireduce_sum(ugtestmx1,ugtestmx1_tmp,nx*nprocy*nprocz*njtest)
        call mpibcast_real_arr(ugtestmx1_tmp,nx*nprocy*nprocz*njtest)
        do jtest=1,njtest
          ugtestmx(:,jtest)=0.
          do jpz=1,nprocz
            do jpy=1,nprocy
              ugtestmx(:,jtest)=ugtestmx(:,jtest)+ugtestmx1_tmp(:,jpy,jpz,jtest)
            enddo
          enddo
        enddo
      endif
!
!  Finally, do this for y-dependent mean fields
!
      if (nprocz>1) then
        call mpireduce_sum(ugtestmy1,ugtestmy1_tmp,ny*nprocy*njtest)
        call mpibcast_real_arr(ugtestmy1_tmp,ny*nprocy*njtest)
        do jtest=1,njtest
          do m=m1,m2
            ugtestmy(m,jtest)=ugtestmy1_tmp(m-m1+1,ipy+1,jtest)
          enddo
        enddo
      endif
!
!  reset headtt
!
      headtt=headtt_save
!
    endsubroutine calc_ltestscalar_pars
!***********************************************************************
    subroutine rescaling_testscalar(f)
!
!  Rescale testscalar by factor rescale_cctest(jtest),
!  which could be different for different testscalar
!
!  26-nov-08/axel: adapted from testfield_z.f90
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=130) :: file
      character (len=5) :: ch
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
          do jtest=1,njtest
            j=icctest+(jtest-1)
            do n=n1,n2
              f(l1:l2,m1:m2,n,j)=rescale_cctest(jtest)*f(l1:l2,m1:m2,n,j)
            enddo
          enddo
          call update_snaptime(file,tccinit,nccinit,dccinit,t, &
            ltestscalar_out,ch,.false.)
        endif
      endif
!
    endsubroutine rescaling_testscalar
!***********************************************************************
    subroutine set_ggtest_G1_G2(G0test,jtest)
!
!  set testscalar
!
!  26-nov-08/axel: adapted from testfield_z.f90
!  27-dec-08/axel: extended to x-dependent mean fields
!
      use Cdata
!
      real, dimension (nx,3) :: G0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: G0test
!
!  set G0test for each of the 2+2 cases
!
      select case(jtest)
      case(1); G0test(:,1)=0.; G0test(:,2)=0.; G0test(:,3)=camp*cz(n)
      case(2); G0test(:,1)=0.; G0test(:,2)=0.; G0test(:,3)=camp*sz(n)
      case(3); G0test(:,1)=camp*cx(:); G0test(:,2)=0.; G0test(:,3)=0.
      case(4); G0test(:,1)=camp*sx(:); G0test(:,2)=0.; G0test(:,3)=0.
      case(5); G0test(:,1)=0.; G0test(:,2)=camp*cy(m); G0test(:,3)=0.
      case(6); G0test(:,1)=0.; G0test(:,2)=camp*sy(m); G0test(:,3)=0.
      case default; G0test(:,:)=0.
      endselect
!
    endsubroutine set_ggtest_G1_G2
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
!
      integer :: iname,inamez,inamexz
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
        idiag_F11z=0; idiag_F21z=0; idiag_F31z=0
        idiag_F12z=0; idiag_F22z=0; idiag_F32z=0
        idiag_kap11=0; idiag_kap21=0; idiag_kap31=0
        idiag_kap12=0; idiag_kap22=0; idiag_kap32=0
        idiag_kap13=0; idiag_kap23=0; idiag_kap33=0
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
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'idiag_kap11=',idiag_kap11
        write(3,*) 'idiag_kap21=',idiag_kap21
        write(3,*) 'idiag_kap31=',idiag_kap31
        write(3,*) 'idiag_kap12=',idiag_kap12
        write(3,*) 'idiag_kap22=',idiag_kap22
        write(3,*) 'idiag_kap32=',idiag_kap32
        write(3,*) 'idiag_kap13=',idiag_kap13
        write(3,*) 'idiag_kap23=',idiag_kap23
        write(3,*) 'idiag_kap33=',idiag_kap33
        write(3,*) 'idiag_c1rms=',idiag_c1rms
        write(3,*) 'idiag_c2rms=',idiag_c2rms
        write(3,*) 'idiag_c3rms=',idiag_c3rms
        write(3,*) 'idiag_c4rms=',idiag_c4rms
        write(3,*) 'idiag_c5rms=',idiag_c5rms
        write(3,*) 'idiag_c6rms=',idiag_c6rms
        write(3,*) 'idiag_c1pt=',idiag_c1pt
        write(3,*) 'idiag_c2pt=',idiag_c2pt
        write(3,*) 'idiag_c3pt=',idiag_c3pt
        write(3,*) 'idiag_c4pt=',idiag_c4pt
        write(3,*) 'idiag_c5pt=',idiag_c5pt
        write(3,*) 'idiag_c6pt=',idiag_c6pt
        write(3,*) 'idiag_F11z=',idiag_F11z
        write(3,*) 'idiag_F21z=',idiag_F21z
        write(3,*) 'idiag_F31z=',idiag_F31z
        write(3,*) 'idiag_F12z=',idiag_F12z
        write(3,*) 'idiag_F22z=',idiag_F22z
        write(3,*) 'idiag_F32z=',idiag_F32z
        write(3,*) 'icctest=',icctest
        write(3,*) 'ntestscalar=',ntestscalar
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'nnamexz=',nnamexz
      endif
!
    endsubroutine rprint_testscalar

endmodule Testscalar
