! $Id$

!  Test-field module for axisymmetric turbulence.
!  No special settings in cparam.local are needed.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestfield = .true.
! CPARAM integer, parameter :: njtest=3
!
! MVAR CONTRIBUTION 9
! MAUX CONTRIBUTION 9
!
! CPARAM logical, parameter :: ltestfield = .true.
! CPARAM logical, parameter :: ltestfield_xy = .false.
! CPARAM logical, parameter :: ltestfield_z = .false.
! CPARAM logical, parameter :: ltestfield_xz  = .false.
!***************************************************************

module Testfield

  use Cdata
  use Messages

  implicit none

  include 'testfield.h'
!
! Slice precalculation buffers
!
  real, target, dimension (:,:,:), allocatable :: bb1_xy,bb1_xy2,bb1_xy3,bb1_xy4
  real, target, dimension (:,:,:), allocatable :: bb1_xz,bb1_xz2,bb1_yz
  real, target, dimension (:,:,:,:,:,:), allocatable :: bb1_r
!
!  cosine and sine function for setting test fields and analysis
!
  real, dimension(nx) :: sx,cx
  real, dimension(my) :: sy,cy
  real, dimension(mz) :: cz,sz,c2z,csz,s2z,c2kz,s2kz,kyBx1,kyBz1,bampz,bampz1
  real :: phase_testfield=.0
!
  character (len=labellen), dimension(ninit) :: initaatest='nothing'
  real, dimension (ninit) :: kx_aatest=1.,ky_aatest=1.,kz_aatest=1.
  real, dimension (ninit) :: phasex_aatest=0.,phasez_aatest=0.
  real, dimension (ninit) :: amplaatest=0.
  integer, dimension (njtest) :: nuxb=0

  ! input parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./)
  real :: taainit=0.,daainit=0.,taainit_previous=0.
  logical :: reinitialize_aatest=.false.
  logical :: zextent=.true.,lsoca=.false.,lsoca_jxb=.true.,lset_bbtest2=.false.
  logical :: luxb_as_aux=.false.,ljxb_as_aux=.false.,linit_aatest=.false.
  logical :: lignore_uxbtestm=.false., lphase_adjust=.false.
  character (len=labellen) :: itestfield='linear'
  real :: ktestfield=1., kxtestfield=1., kytestfield=1.
  real :: ktestfield1=1., kytestfield1, ky1, kz1
  real :: ztestfield_offset=0.
  real :: lin_testfield=0.,lam_testfield=0.,om_testfield=0.,delta_testfield=0.
  real :: delta_testfield_next=0., delta_testfield_time=0.
  integer :: naainit
  real :: bamp=1.,bamp1=1.,bamp12=1.
  namelist /testfield_init_pars/ &
       B_ext,zextent,initaatest, &
       amplaatest,kx_aatest,ky_aatest,kz_aatest, &
       phasex_aatest,phasez_aatest, &
       luxb_as_aux,ljxb_as_aux

  ! run parameters
  real :: etatest=0.,etatest1=0.
  real :: tau_aatest=0.,tau1_aatest=0.
  real :: ampl_fcont_aatest=1.
  real, dimension(njtest) :: rescale_aatest=0.
  logical :: ltestfield_newz=.true.,leta_rank2=.true.
  logical :: ltestfield_taver=.false.
  logical :: llorentzforce_testfield=.false.
  logical :: lforcing_cont_aatest=.false.
  logical :: ltestfield_artifric=.false.
  namelist /testfield_run_pars/ &
       B_ext,reinitialize_aatest,zextent,lsoca,lsoca_jxb, &
       lset_bbtest2,etatest,etatest1,itestfield, &
       ktestfield,kxtestfield,kytestfield, &
       ztestfield_offset, &
       lin_testfield,lam_testfield,om_testfield,delta_testfield, &
       ltestfield_newz,leta_rank2,lphase_adjust,phase_testfield, &
       ltestfield_taver,llorentzforce_testfield, &
       luxb_as_aux,ljxb_as_aux,lignore_uxbtestm, &
       lforcing_cont_aatest,ampl_fcont_aatest, &
       daainit,linit_aatest,bamp, &
       rescale_aatest,tau_aatest

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_alpPERP=0    ! DIAG_DOC: $\alpha_\perp$
  integer :: idiag_alpPARA=0    ! DIAG_DOC: $\alpha_\perp$
  integer :: idiag_gam=0        ! DIAG_DOC: $\gamma$
  integer :: idiag_betPERP=0    ! DIAG_DOC: $\beta_\perp$
  integer :: idiag_betPARA=0    ! DIAG_DOC: $\beta_\perp$
  integer :: idiag_del=0        ! DIAG_DOC: $\delta$
  integer :: idiag_kapPERP=0    ! DIAG_DOC: $\kappa_\perp$
  integer :: idiag_kapPARA=0    ! DIAG_DOC: $\kappa_\perp$
  integer :: idiag_mu=0         ! DIAG_DOC: $\mu$
  integer :: idiag_alpPERPz=0   ! DIAG_DOC: $\alpha_\perp(z)$
  integer :: idiag_alpPARAz=0   ! DIAG_DOC: $\alpha_\perp(z)$
  integer :: idiag_gamz=0       ! DIAG_DOC: $\gamma(z)$
  integer :: idiag_betPERPz=0   ! DIAG_DOC: $\beta_\perp(z)$
  integer :: idiag_betPARAz=0   ! DIAG_DOC: $\beta_\perp(z)$
  integer :: idiag_delz=0       ! DIAG_DOC: $\delta(z)$
  integer :: idiag_kapPERPz=0   ! DIAG_DOC: $\kappa_\perp(z)$
  integer :: idiag_kapPARAz=0   ! DIAG_DOC: $\kappa_\perp(z)$
  integer :: idiag_muz=0        ! DIAG_DOC: $\mu(z)$
  integer :: idiag_bx1pt=0      ! DIAG_DOC: $b_x^{1}$
  integer :: idiag_bx2pt=0      ! DIAG_DOC: $b_x^{2}$
  integer :: idiag_bx3pt=0      ! DIAG_DOC: $b_x^{3}$
  integer :: idiag_b1rms=0      ! DIAG_DOC: $\left<b_{1}^2\right>^{1/2}$
  integer :: idiag_b2rms=0      ! DIAG_DOC: $\left<b_{2}^2\right>^{1/2}$
  integer :: idiag_b3rms=0      ! DIAG_DOC: $\left<b_{3}^2\right>^{1/2}$
  integer :: ivid_bb1=0
  integer :: idiag_bcosphz=0
  integer :: idiag_bsinphz=0
!
!  arrays for horizontally averaged uxb and jxb
!
  real, dimension (nz,3,njtest) :: uxbtestm,jxbtestm
  real, dimension (nx,3,njtest) :: Eipq,bpq

  contains
!***********************************************************************
    subroutine register_testfield()
!
!  Initialise variables which should know that we solve for the vector
!  potential: iaatest, etc; increase nvar accordingly
!
!   3-jun-05/axel: adapted from register_magnetic
!
      use Cdata
      use FArrayManager
      use Mpicomm
      use Sub
!
!  Register test field.
!
      ntestfield = 3*njtest
      call farray_register_pde('aatest',iaatest,array=ntestfield)
      call farray_index_append('ntestfield',ntestfield)
!
!  Set first and last index of text field
!  Note: iaxtest, iaytest, and iaztest are initialized to the first test field.
!  These values are used in this form in start, but later overwritten.
!
      iaxtest=iaatest
      iaytest=iaatest+1
      iaztest=iaatest+2
      iaztestpq=iaatest+ntestfield-1
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
          if (nvar < mvar) write(4,*) ',aatest $'
          if (nvar == mvar) write(4,*) ',aatest'
        else
          write(4,*) ',aatest $'
        endif
        write(15,*) 'aatest = fltarr(mx,my,mz,ntestfield)*one'
      endif
!
    endsubroutine register_testfield
!***********************************************************************
    subroutine initialize_testfield(f)
!
!  Perform any post-parameter-read initialization
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use FArrayManager
      use Slices_methods, only: alloc_slice_buffers
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(mz) :: ztestfield, c, s
      real :: ktestfield_effective
      integer :: jtest
!
!  Precalculate etatest if 1/etatest (==etatest1) is given instead
!
      if (etatest1 /= 0.) etatest=1./etatest1
      if (lroot) print*,'initialize_testfield: etatest=',etatest
!
!  set cosine and sine function for setting test fields and analysis
!  Choice of using rescaled z-array or original z-array
!  Define ktestfield_effective to deal with boxes bigger than 2pi.
!
      if (ltestfield_newz) then
        ktestfield_effective=ktestfield*(2.*pi/Lz)
        ztestfield=ktestfield_effective*(z-z0)-pi
      else
        ktestfield_effective=ktestfield
        ztestfield=z*ktestfield_effective
      endif
      cz=cos(ztestfield)
      sz=sin(ztestfield)
      c2kz=cos(2*ztestfield)
      s2kz=sin(2*ztestfield)
!
!  define sinky and cosky, called sy and cy
!
      sy=sin(kytestfield*y)
      cy=cos(kytestfield*y)
!
!  define sinkx
!
      sx=sin(kxtestfield*x(l1:l2))
      cx=cos(kxtestfield*x(l1:l2))
!
!  calculate cosz*sinz, cos^2, and sinz^2, to take moments with
!  of alpij and etaij. This is useful if there is a mean Beltrami field
!  in the main calculation (lmagnetic=.true.) with phase zero.
!  Optionally, one can determine the phase in the actual field
!  and modify the following calculations in testfield_after_boundary.
!  They should also be with respect to k*z, not just z.
!
      c=cz
      s=sz
      c2z=c**2
      s2z=s**2
      csz=c*s
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
      if (ktestfield==0) then
        ktestfield1=1.
      else
        ktestfield1=1./ktestfield_effective
      endif
      kz1=ktestfield1
!
!  Calculate inverse kytestfield, but only if different from zero
!  Introduce kyBz1 as abbreviation of 1/(kytestfield*Btestz)
!
      if (kytestfield==0) then
        kytestfield1=1.
      else
        kytestfield1=1./kytestfield
      endif
      kyBx1=kytestfield1/(z(n)+ztestfield_offset)
      kyBz1=kytestfield1/(z(n)+ztestfield_offset)
      ky1=kytestfield1
!
!  calculate inverse testfield amplitude (unless it is set to zero)
!
      if (bamp==0.) then
        bamp1=1.
        bamp12=1.
      else
        bamp1=1./bamp
        bamp12=1./bamp**2
      endif
!
      bampz=bamp*(ztestfield_offset+z)
!
!  avoid dividing by zero
!
      bampz1=1.
      where (abs(bampz)>=(.1*dz))
        bampz1=1./bampz
      endwhere
!
!  set to zero and then rescale the testfield
!  (in future, could call something like init_aa_simple)
!
      if (reinitialize_aatest) then
        do jtest=1,njtest
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2
          f(:,:,:,iaxtest:iaztest)=rescale_aatest(jtest)*f(:,:,:,iaxtest:iaztest)
        enddo
      endif
!
!  set lrescaling_testfield=T if linit_aatest=T
!
      if (linit_aatest) then
        lrescaling_testfield=.true.
      endif
!
!  check for possibility of artificial friction force
!
      if (tau_aatest/=0.) then
        ltestfield_artifric=.true.
        tau1_aatest=1./tau_aatest
        if (lroot) print*,'initialize_testfield: tau1_aatest=',tau1_aatest
      endif

!  Register an extra aux slot for uxb if requested (so uxb is written
!  to snapshots and can be easily analyzed later). For this to work you
!  must reserve enough auxiliary workspace by setting, for example,
!     ! MAUX CONTRIBUTION 9
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
!
!  After a reload, we need to rewrite index.pro, but the auxiliary
!  arrays are already allocated and must not be allocated again.
!
      if (luxb_as_aux) then
        if (iuxbtest==0) then
          call farray_register_auxiliary('uxb',iuxbtest,vector=3*njtest)
        else
          if (lroot) print*, 'initialize_testfield: iuxbtest = ', iuxbtest
          call farray_index_append('iuxbtest',iuxbtest)
        endif
      endif
!
!  possibility of using jxb as auxiliary array (is intended to be
!  used in connection with testflow method)
!
      if (ljxb_as_aux) then
        if (ijxbtest==0) then
          call farray_register_auxiliary('jxb',ijxbtest,vector=3*njtest)
        else
          if (lroot) print*, 'initialize_testfield: ijxbtest = ', ijxbtest
          call farray_index_append('ijxbtest',ijxbtest)
        endif
      endif
!
      if (ivid_bb1/=0) &
        call alloc_slice_buffers(bb1_xy,bb1_xz,bb1_yz,bb1_xy2,bb1_xy3,bb1_xy4,bb1_xz2,bb1_r)
!
!  write testfield information to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testfield_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'zextent=',merge(1,0,zextent)
        write(1,'(a,i1)') 'lsoca='  ,merge(1,0,lsoca)
        write(1,'(a,i1)') 'lsoca_jxb='  ,merge(1,0,lsoca_jxb)
        write(1,'(3a)') "itestfield='",trim(itestfield)//"'"
        write(1,'(a,f5.2)') 'ktestfield=',ktestfield
        write(1,'(a,f5.2)') 'kytestfield=',kytestfield
        write(1,'(a,f7.4)') 'lin_testfield=',lin_testfield
        write(1,'(a,f7.4)') 'lam_testfield=',lam_testfield
        write(1,'(a,f7.4)') 'om_testfield=', om_testfield
        write(1,'(a,f7.4)') 'delta_testfield=',delta_testfield
        close(1)
      endif
!
    endsubroutine initialize_testfield
!***********************************************************************
    subroutine init_aatest(f)
!
!  initialise testfield; called from start.f90
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use Mpicomm
      use Initcond
      use Sub
      use InitialCondition, only: initial_condition_aatest
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      do j=1,ninit

      select case (initaatest(j))

      case ('zero'); f(:,:,:,iaatest:iaatest+3*njtest-1)=0.
      case ('gaussian-noise-1'); call gaunoise(amplaatest(j),f,iaxtest+0,iaztest+0)
      case ('gaussian-noise-2'); call gaunoise(amplaatest(j),f,iaxtest+3,iaztest+3)
      case ('gaussian-noise-3'); call gaunoise(amplaatest(j),f,iaxtest+6,iaztest+6)
      case ('gaussian-noise-4'); call gaunoise(amplaatest(j),f,iaxtest+9,iaztest+9)
      case ('gaussian-noise-5'); call gaunoise(amplaatest(j),f,iaxtest+12,iaztest+12)
      case ('sinwave-x-1'); call sinwave(amplaatest(j),f,iaxtest+0+1,kx=kx_aatest(j))
      case ('sinwave-x-2'); call sinwave(amplaatest(j),f,iaxtest+3+1,kx=kx_aatest(j))
      case ('sinwave-x-3'); call sinwave(amplaatest(j),f,iaxtest+6+1,kx=kx_aatest(j))
      case ('Beltrami-x-1'); call beltrami(amplaatest(j),f,iaxtest+0,kx=-kx_aatest(j),phase=phasex_aatest(j))
      case ('Beltrami-z-1'); call beltrami(amplaatest(j),f,iaxtest+0,kz=-kz_aatest(j),phase=phasez_aatest(j))
      case ('Beltrami-z-2'); call beltrami(amplaatest(j),f,iaxtest+3,kz=-kz_aatest(j),phase=phasez_aatest(j))
      case ('Beltrami-z-3'); call beltrami(amplaatest(j),f,iaxtest+6,kz=-kz_aatest(j),phase=phasez_aatest(j))
      case ('Beltrami-z-4'); call beltrami(amplaatest(j),f,iaxtest+9,kz=-kz_aatest(j),phase=phasez_aatest(j))
      case ('Beltrami-z-5'); call beltrami(amplaatest(j),f,iaxtest+12,kz=-kz_aatest(j),phase=phasez_aatest(j))
      case ('nothing'); !(do nothing)

      case default
        call fatal_error("init_aatest","no such initaatest: "//trim(initaatest(j)))
      endselect
      enddo
!
!  Interface for user's own subroutine
!
      if (linitial_condition) call initial_condition_aatest(f)
!
    endsubroutine init_aatest
!***********************************************************************
    subroutine pencil_criteria_testfield()
!
!   All pencils that the Testfield module depends on are specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
      lpenc_requested(i_uu)=.true.
      if (lforcing_cont_aatest) lpenc_requested(i_fcont)=.true.
!
    endsubroutine pencil_criteria_testfield
!***********************************************************************
    subroutine pencil_interdep_testfield(lpencil_in)
!
!  Interdependency among pencils from the Testfield module is specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use General, only: keep_compiler_quiet
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_testfield
!***********************************************************************
    subroutine read_testfield_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testfield_init_pars, IOSTAT=iostat)
!
    endsubroutine read_testfield_init_pars
!***********************************************************************
    subroutine write_testfield_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testfield_init_pars)
!
    endsubroutine write_testfield_init_pars
!***********************************************************************
    subroutine read_testfield_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=testfield_run_pars, IOSTAT=iostat)
!
    endsubroutine read_testfield_run_pars
!***********************************************************************
    subroutine write_testfield_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=testfield_run_pars)
!
    endsubroutine write_testfield_run_pars
!***********************************************************************
    subroutine daatest_dt(f,df,p)
!
!  testfield evolution:
!
!  calculate da^(pq)/dt=Uxb^(pq)+uxB^(pq)+uxb-<uxb>+eta*del2A^(pq),
!    where p=1,2 and q=1 (if B11-B21) and optionally q=2 (if B11-B22)
!
!  also calculate corresponding Lorentz force in connection with
!  testflow method
!
!   3-jun-05/axel: coded
!  16-mar-08/axel: Lorentz force added for testfield method
!  25-jan-09/axel: added Maxwell stress tensor calculation
!
      use Hydro, only: uumz,lcalc_uumeanz
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx,3) :: uxB,B0test=0,bbtest
      real, dimension (nx,3) :: uxbtest,duxbtest,jxbtest,djxbrtest,eetest
      real, dimension (nx,3) :: J0test=0,jxB0rtest,J0xbrtest
      real, dimension (nx,3,3,njtest) :: Mijpq
      real, dimension (nx,3,njtest) :: jpq
      real, dimension (nx,3) :: del2Atest,uufluct
      real, dimension (nx,3) :: del2Atest2,graddivatest,aatest,jjtest,jxbrtest
      real, dimension (nx,3,3) :: aijtest
      real, dimension (nx) :: jbpq,Epq2,s2kzDF1,s2kzDF2,unity=1.
      integer :: jtest,j,nl
      logical,save :: ltest_uxb=.false.,ltest_jxb=.false.
      real, dimension(nx) :: diffus_eta
!
      intent(in)     :: f,p
      intent(inout)  :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'daatest_dt: SOLVE'
      if (headtt) then
        if (iaxtest /= 0) call identify_bcs('Axtest',iaxtest)
        if (iaytest /= 0) call identify_bcs('Aytest',iaytest)
        if (iaztest /= 0) call identify_bcs('Aztest',iaztest)
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
!  Multiply by exponential factor if lam_testfield is different from zero.
!  Allow also for linearly increasing testfields
!  Keep bamp1=1 for oscillatory fields.
!
      if (lam_testfield/=0..or.lin_testfield/=0. .or. &
          om_testfield/=0..or.delta_testfield/=0.) then
        if (lam_testfield/=0.) then
          taainit_previous=taainit-daainit
          bamp=exp(lam_testfield*(t-taainit_previous))
          bamp1=1./bamp
        endif
        if (lin_testfield/=0.) then
          taainit_previous=taainit-daainit
          bamp=lin_testfield*(t-taainit_previous-daainit/2.)
          bamp1=1./bamp
        endif
        if (om_testfield/=0.) then
          bamp=cos(om_testfield*t)
          bamp1=1.
        endif
        if (delta_testfield/=0.) then
          if (t>=delta_testfield_next.and.t<=(delta_testfield_next+dt)) then
            delta_testfield_time=t
            delta_testfield_next=t+delta_testfield
          endif
          if (t>=delta_testfield_time-.1*dt.and.t<=delta_testfield_time+.1*dt) then
            bamp=1.
          else
            bamp=0.
          endif
          bamp1=1.
        endif
      endif
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further down in the file.
!
      nl=n-n1+1
      do jtest=1,njtest
        iaxtest=iaatest+3*(jtest-1)
        iaztest=iaxtest+2
        call del2v(f,iaxtest,del2Atest)
        select case (itestfield)
          case ('linear'); call set_bbtest_linear(B0test,jtest)
          case ('linear_old'); call set_bbtest_linear_old(B0test,jtest)
          case ('sxsysz'); call set_bbtest_sxsysz(B0test,jtest)
          case ('sinkz'); call set_bbtest_sinkz(B0test,jtest)
          case ('B=0') !(dont do anything)
        case default
          call fatal_error('daatest_dt','no such itestfield: '//trim(itestfield))
        endselect
!
!  add an external field, if present
!
        if (B_ext(1)/=0.) B0test(:,1)=B0test(:,1)+B_ext(1)
        if (B_ext(2)/=0.) B0test(:,2)=B0test(:,2)+B_ext(2)
        if (B_ext(3)/=0.) B0test(:,3)=B0test(:,3)+B_ext(3)
!
        call cross_mn(uufluct,B0test,uxB)
        if (lsoca) then
          df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+uxB+etatest*del2Atest
        else
!
!  use f-array for uxb (if space has been allocated for this) and
!  if we don't test (i.e. if ltest_uxb=.false.)
!
          if (iuxbtest/=0.and..not.ltest_uxb) then
            uxbtest=f(l1:l2,m,n,iuxbtest+3*(jtest-1):iuxbtest+3*jtest-1)
          else
            aatest=f(l1:l2,m,n,iaxtest:iaztest)
            call gij(f,iaxtest,aijtest,1)
            call curl_mn(aijtest,bbtest,aatest)
            call cross_mn(p%uu,bbtest,uxbtest)
          endif
!
!  subtract average emf, unless we ignore the <uxb> term (lignore_uxbtestm=T)
!
          if (lignore_uxbtestm) then
            duxbtest=uxbtest
          else
            do j=1,3
              duxbtest(:,j)=uxbtest(:,j)-uxbtestm(n,j,jtest)
            enddo
          endif
!
!  advance test field equation
!
          df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
            +uxB+etatest*del2Atest+duxbtest
        endif
!
!  add possibility of forcing that is not delta-correlated in time
!
        if (lforcing_cont_aatest) &
          df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
                                        +ampl_fcont_aatest*p%fcont(:,:,1)  ! first forcing added
!
!  add possibility of artificial friction
!
        if (ltestfield_artifric) &
          df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
              -tau1_aatest*f(l1:l2,m,n,iaxtest:iaztest)
!
!  calculate alpha, begin by calculating uxbtest (if not already done above)
!
        if ((ldiagnos.or.l1davgfirst).and. (lsoca.or.ltest_uxb.or. &
           idiag_b1rms/=0.or.idiag_b2rms/=0.or.idiag_b3rms/=0)) then
          aatest=f(l1:l2,m,n,iaxtest:iaztest)
          call gij(f,iaxtest,aijtest,1)
          call curl_mn(aijtest,bbtest,aatest)
          call cross_mn(p%uu,bbtest,uxbtest)
        endif
        bpq(:,:,jtest)=bbtest
        Eipq(:,:,jtest)=uxbtest*bamp1
        if (ldiagnos) jpq(:,:,jtest)=jjtest
      enddo
!
!  diffusive time step, just take the max of diffus_eta (if existent)
!  and whatever is calculated here
!
      if (lfirst.and.ldt) then
        diffus_eta=etatest*dxyz_2
        maxdiffus=max(maxdiffus,diffus_eta)
      endif

      call calc_diagnostics_testfield(p)
!
    endsubroutine daatest_dt
!***********************************************************************
    subroutine calc_diagnostics_testfield(p)
!
      use Diagnostics
      use Slices_methods, only: store_slices
      use Sub, only:dot2
!
      type (pencil_case) :: p
!
      real, dimension (nx) :: bpq2
      integer :: i1=1, i2=2, i3=3, i4=4
!
!  in the following block, we have already swapped the 4-6 entries with 7-9
!  The g95 compiler doesn't like to see an index that is out of bounds,
!  so prevent this warning by writing i3=3 and i4=4
!
      if (ldiagnos) then
!
!  averages of alpha and eta
!  Don't subtract E0 field if iE=0
!
        select case (itestfield)
!
!  test fields linear in z, where the mean field depends onm x and y
!
        case ('linear')
!
!  Volume averages first.
!
          if (idiag_gam    /=0) call sum_mn_name(+2.*sx*sy(m)*(Eipq(:,1,1)-Eipq(:,2,1)),idiag_gam    )
          if (idiag_alpPERP/=0) call sum_mn_name(-2.*sx*sy(m)*(Eipq(:,1,1)+Eipq(:,2,1)),idiag_alpPERP)
          if (idiag_alpPARA/=0) call sum_mn_name(-4.*sx*sy(m)* Eipq(:,3,1)             ,idiag_alpPARA)
!
          if (idiag_mu     /=0) call sum_mn_name(-4.*sy(m)*(sx*(Eipq(:,2,2)+.5*bampz(n)*(Eipq(:,1,1)-Eipq(:,2,1))) &
            +cx*bampz1(n)*Eipq(:,2,3)),idiag_mu     )
          if (idiag_betPERP/=0) call sum_mn_name(-2.*sy(m)*(sx*(Eipq(:,2,2)+.5*bampz(n)*(Eipq(:,1,1)-Eipq(:,2,1))) &
            -cx*bampz1(n)*Eipq(:,2,3)),idiag_betPERP)
          if (idiag_betPARA/=0) call sum_mn_name(+4.*cy(m)*                                                        &
             sx*bampz1(n)*Eipq(:,3,2) ,idiag_betPARA)
!
          if (idiag_del    /=0) call sum_mn_name(+2.*sy(m)*(sx*(Eipq(:,1,2)-.5*bampz(n)*(Eipq(:,1,1)+Eipq(:,2,1))) &
            -cx*bampz1(n)*Eipq(:,1,3)),idiag_del    )
          if (idiag_kapPERP/=0) call sum_mn_name(-4.*sy(m)*(sx*(Eipq(:,1,2)-.5*bampz(n)*(Eipq(:,1,1)+Eipq(:,2,1))) &
            +cx*bampz1(n)*Eipq(:,1,3)),idiag_kapPERP)
          if (idiag_kapPARA/=0) call sum_mn_name(-4.*sy(m)* sx*(Eipq(:,3,3)   -bampz(n)* Eipq(:,3,1)             ) &
                                      ,idiag_kapPARA)
!
!  Same, but for z-dependent output (xy-averages)
!
          call xysum_mn_name_z(+2.*sx*sy(m)*(Eipq(:,1,1)-Eipq(:,2,1)),idiag_gamz    )
          call xysum_mn_name_z(-2.*sx*sy(m)*(Eipq(:,1,1)+Eipq(:,2,1)),idiag_alpPERPz)
          call xysum_mn_name_z(-4.*sx*sy(m)* Eipq(:,3,1)             ,idiag_alpPARAz)
!
          call xysum_mn_name_z(-4.*sy(m)*(sx*(Eipq(:,2,2)+.5*bampz(n)*(Eipq(:,1,1)-Eipq(:,2,1))) &
            +cx*bampz1(n)*Eipq(:,2,3)),idiag_muz     )
          call xysum_mn_name_z(-2.*sy(m)*(sx*(Eipq(:,2,2)+.5*bampz(n)*(Eipq(:,1,1)-Eipq(:,2,1))) &
            -cx*bampz1(n)*Eipq(:,2,3)),idiag_betPERPz)
          call xysum_mn_name_z(+4.*cy(m)*                                                        &
             sx*bampz1(n)*Eipq(:,3,2) ,idiag_betPARAz)
!
          call xysum_mn_name_z(+2.*sy(m)*(sx*(Eipq(:,1,2)-.5*bampz(n)*(Eipq(:,1,1)+Eipq(:,2,1))) &
            -cx*bampz1(n)*Eipq(:,1,3)),idiag_delz    )
          call xysum_mn_name_z(-4.*sy(m)*(sx*(Eipq(:,1,2)-.5*bampz(n)*(Eipq(:,1,1)+Eipq(:,2,1))) &
            +cx*bampz1(n)*Eipq(:,1,3)),idiag_kapPERPz)
          call xysum_mn_name_z(-4.*sy(m)* sx*(Eipq(:,3,3)   -bampz(n)* Eipq(:,3,1)             ) &
                                      ,idiag_kapPARAz)
!
!  test fields linear in z, but mean field depends only on y
!
        case ('linear_old')
          if (idiag_gam    /=0) call sum_mn_name(+.5*(Eipq(:,1,1)-Eipq(:,2,1)),idiag_gam    )
          if (idiag_alpPERP/=0) call sum_mn_name(-.5*(Eipq(:,1,1)+Eipq(:,2,1)),idiag_alpPERP)
          if (idiag_alpPARA/=0) call sum_mn_name(-    Eipq(:,3,1)             ,idiag_alpPARA)
!
          if (idiag_mu     /=0) call sum_mn_name(-2.*(sy(m)*Eipq(:,2,2)+kyBz1(n)*cy(m)*Eipq(:,1,3)),idiag_mu     )
          if (idiag_betPERP/=0) call sum_mn_name(-    sy(m)*Eipq(:,2,2)-kyBz1(n)*cy(m)*Eipq(:,1,3),idiag_betPERP)
          if (idiag_betPARA/=0) call sum_mn_name(+2.*                   kyBx1(n)*cy(m)*Eipq(:,3,3),idiag_betPARA)
!
          if (idiag_del    /=0) call sum_mn_name(+    sy(m)*Eipq(:,1,2)-kyBz1(n)*cy(m)*Eipq(:,2,3),idiag_del    )
          if (idiag_kapPERP/=0) call sum_mn_name(-2.*(sy(m)*Eipq(:,1,2)+kyBz1(n)*cy(m)*Eipq(:,2,3)),idiag_kapPERP)
          if (idiag_kapPARA/=0) call sum_mn_name(-2.* sy(m)*Eipq(:,3,3)                            ,idiag_kapPARA)
!
!  test fields sinusoidal also in the z direction
!
        case ('sinkz')
          if (idiag_gam    /=0) call sum_mn_name(+2*sy(m)*sz(n)*(Eipq(:,1,1)-Eipq(:,2,1)),idiag_gam    )
          if (idiag_alpPERP/=0) call sum_mn_name(-2*sy(m)*sz(n)*(Eipq(:,1,1)+Eipq(:,2,1)),idiag_alpPERP)
          if (idiag_alpPARA/=0) call sum_mn_name(-4*sy(m)*sz(n)* Eipq(:,3,1)             ,idiag_alpPARA)
!
          if (idiag_mu     /=0) call sum_mn_name(-4*(kz1*sy(m)*cz(n)*Eipq(:,2,2)-ky1*cy(m)*sz(n)*Eipq(:,1,3)),idiag_mu     )
          if (idiag_betPERP/=0) call sum_mn_name(-2*(kz1*sy(m)*cz(n)*Eipq(:,2,2)+ky1*cy(m)*sz(n)*Eipq(:,1,3)),idiag_betPERP)
          if (idiag_betPARA/=0) call sum_mn_name(+4*                             ky1*cy(m)*sz(n)*Eipq(:,3,3) ,idiag_betPARA)
!
          if (idiag_del    /=0) call sum_mn_name(+2*(kz1*sy(m)*cz(n)*Eipq(:,1,2)-ky1*cy(m)*sz(n)*Eipq(:,2,3)),idiag_del    )
          if (idiag_kapPERP/=0) call sum_mn_name(-2*(kz1*sy(m)*cz(n)*Eipq(:,1,2)+ky1*cy(m)*sz(n)*Eipq(:,2,3)),idiag_kapPERP)
          if (idiag_kapPARA/=0) call sum_mn_name(-2* kz1*sy(m)*cz(n)*Eipq(:,3,3)                             ,idiag_kapPARA)
!
!  Test fields sinusoidal in all three directions.
!
        case ('sxsysz')
!
!  Volume averages first.
!
          if (idiag_gam    /=0) call sum_mn_name(+4*sx*sy(m)*sz(n)*(Eipq(:,1,1)-Eipq(:,2,1)),idiag_gam    )
          if (idiag_alpPERP/=0) call sum_mn_name(-4*sx*sy(m)*sz(n)*(Eipq(:,1,1)+Eipq(:,2,1)),idiag_alpPERP)
          if (idiag_alpPARA/=0) call sum_mn_name(-8*sx*sy(m)*sz(n)* Eipq(:,3,1)             ,idiag_alpPARA)
!
          if (idiag_mu     /=0) call sum_mn_name(-8*sx*(kz1*sy(m)*cz(n)*Eipq(:,2,2)-ky1*cy(m)*sz(n)*Eipq(:,1,3)),idiag_mu     )
          if (idiag_betPERP/=0) call sum_mn_name(-4*sx*(kz1*sy(m)*cz(n)*Eipq(:,2,2)+ky1*cy(m)*sz(n)*Eipq(:,1,3)),idiag_betPERP)
          if (idiag_betPARA/=0) call sum_mn_name(+8*sx*                             ky1*cy(m)*sz(n)*Eipq(:,3,2) ,idiag_betPARA)
!
          if (idiag_del    /=0) call sum_mn_name(+4*sx*(kz1*sy(m)*cz(n)*Eipq(:,1,2)-ky1*cy(m)*sz(n)*Eipq(:,2,3)),idiag_del    )
          if (idiag_kapPERP/=0) call sum_mn_name(-8*sx*(kz1*sy(m)*cz(n)*Eipq(:,1,2)+ky1*cy(m)*sz(n)*Eipq(:,2,3)),idiag_kapPERP)
          if (idiag_kapPARA/=0) call sum_mn_name(-8*sx* kz1*sy(m)*cz(n)*Eipq(:,3,3)                             ,idiag_kapPARA)
!
!  xy-averages next.
!
          call xysum_mn_name_z(+4*sx*sy(m)*sz(n)*(Eipq(:,1,1)-Eipq(:,2,1)),idiag_gamz    )
          call xysum_mn_name_z(-4*sx*sy(m)*sz(n)*(Eipq(:,1,1)+Eipq(:,2,1)),idiag_alpPERPz)
          call xysum_mn_name_z(-8*sx*sy(m)*sz(n)* Eipq(:,3,1)             ,idiag_alpPARAz)
!
          call xysum_mn_name_z(-8*sx*(kz1*sy(m)*cz(n)*Eipq(:,2,2)-ky1*cy(m)*sz(n)*Eipq(:,1,3)),idiag_muz     )
          call xysum_mn_name_z(-4*sx*(kz1*sy(m)*cz(n)*Eipq(:,2,2)+ky1*cy(m)*sz(n)*Eipq(:,1,3)),idiag_betPERPz)
          call xysum_mn_name_z(+8*sx*                             ky1*cy(m)*sz(n)*Eipq(:,3,2) ,idiag_betPARAz)
!
          call xysum_mn_name_z(+4*sx*(kz1*sy(m)*cz(n)*Eipq(:,1,2)-ky1*cy(m)*sz(n)*Eipq(:,2,3)),idiag_delz    )
          call xysum_mn_name_z(-8*sx*(kz1*sy(m)*cz(n)*Eipq(:,1,2)+ky1*cy(m)*sz(n)*Eipq(:,2,3)),idiag_kapPERPz)
          call xysum_mn_name_z(-8*sx* kz1*sy(m)*cz(n)*Eipq(:,3,3)                             ,idiag_kapPARAz)
!
        case default
          call fatal_error('daatest_dt','no such itestfield: '//trim(itestfield))
        endselect
!
!  diagnostics for single points
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          call save_name(bpq(lpoint-nghost,1,i1),idiag_bx1pt)
          call save_name(bpq(lpoint-nghost,1,i2),idiag_bx2pt)
        endif
!
!  rms values of small scales fields bpq in response to the test fields Bpq
!
        if (idiag_b1rms/=0) then
          call dot2(bpq(:,:,1),bpq2)
          call sum_mn_name(bpq2,idiag_b1rms,lsqrt=.true.)
        endif
!
        if (idiag_b2rms/=0) then
          call dot2(bpq(:,:,i2),bpq2)
          call sum_mn_name(bpq2,idiag_b2rms,lsqrt=.true.)
        endif
!
        if (idiag_b3rms/=0) then
          call dot2(bpq(:,:,i3),bpq2)
          call sum_mn_name(bpq2,idiag_b3rms,lsqrt=.true.)
        endif
!
      endif
!
!  write B-slices for output in wvid in run.f90
!
      if (lvideo.and.lfirst.and.ivid_bb1/=0) &
        call store_slices(bpq(:,:,1),bb1_xy,bb1_xz,bb1_yz,bb1_xy2,bb1_xy3,bb1_xy4,bb1_xz2,bb1_r)
!
    endsubroutine calc_diagnostics_testfield
!***********************************************************************
    subroutine get_slices_testfield(f,slices)
!
!  Write slices for animation of magnetic variables.
!
!  12-sep-09/axel: adapted from the corresponding magnetic routine
!
      use General, only: keep_compiler_quiet
      use Slices_methods, only: assign_slices_vec
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Magnetic field
!
      case ('bb1')
        call assign_slices_vec(slices,bb1_xy,bb1_xz,bb1_yz,bb1_xy2,bb1_xy3,bb1_xy4,bb1_xz2,bb1_r)

      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine get_slices_testfield
!***********************************************************************
    subroutine testfield_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!    4-oct-18/axel+nishant: adapted from testflow
!
      use General, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine testfield_before_boundary
!***********************************************************************
    subroutine testfield_after_boundary(f)
!
!  calculate <uxb>, which is needed when lsoca=.false.
!
!  21-jan-06/axel: coded
!
      use Cdata
      use Sub
      use Hydro, only: calc_pencils_hydro
      use Mpicomm, only: mpiallreduce_sum, mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mz) :: c,s
!
      real, dimension (nz,nprocz,3,njtest) :: uxbtestm1=0.,uxbtestm1_tmp=0.
      real, dimension (nz,nprocz,3,njtest) :: jxbtestm1=0.,jxbtestm1_tmp=0.
!
      real, dimension (nx,3,3) :: aijtest,bijtest
      real, dimension (nx,3) :: aatest,bbtest,jjtest,uxbtest,jxbtest
      real, dimension (nx,3) :: del2Atest2,graddivatest
      integer :: jtest,j,juxb,jjxb,nl
      logical :: headtt_save
      real :: fac, bcosphz, bsinphz, fac1=0., fac2=1.
!
      intent(inout) :: f
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac=1./nxygrid
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
      do jtest=1,njtest
        iaxtest=iaatest+3*(jtest-1)
        iaztest=iaxtest+2
        if (lsoca) then
          uxbtestm(:,:,jtest)=0.
        else
!
!  prepare coefficients for possible time integration of uxb,
!  But do this only when we are on the last time step
!
          if (ltestfield_taver) then
            if (llast) then
              fac2=1./(nuxb(jtest)+1.)
              fac1=1.-fac2
            endif
          endif
!
!  calculate uxb, and put it into f-array
!  If ltestfield_taver is turned on, the time evolution becomes unstable,
!  so that does not currently work.
!
          do n=n1,n2
            nl=n-n1+1
            uxbtestm(nl,:,jtest)=0.
            do m=m1,m2
              aatest=f(l1:l2,m,n,iaxtest:iaztest)
              call gij(f,iaxtest,aijtest,1)
              call curl_mn(aijtest,bbtest,aatest)
              call cross_mn(f(l1:l2,m,n,iux:iuz),bbtest,uxbtest)
              juxb=iuxbtest+3*(jtest-1)
              if (ltestfield_taver) then
                if (llast) then
                  if (iuxbtest/=0) f(l1:l2,m,n,juxb:juxb+2)= &
                          fac1*f(l1:l2,m,n,juxb:juxb+2)+fac2*uxbtest
                endif
              else
                if (iuxbtest/=0) f(l1:l2,m,n,juxb:juxb+2)=uxbtest
              endif
              uxbtestm(nl,:,jtest)=uxbtestm(nl,:,jtest)+fac*sum(uxbtest,1)
              headtt=.false.
            enddo
            uxbtestm1(:,ipz+1,:,:)=uxbtestm
          enddo
!
!  update counter, but only when we are on the last substep
!
          if (ltestfield_taver) then
            if (llast) then
              nuxb(jtest)=nuxb(jtest)+1
            endif
          endif
        endif
      enddo
!
!  do communication for array of size nz*nprocz*3*njtest
!
      if (nprocxy>1) then
        call mpiallreduce_sum(uxbtestm1,uxbtestm1_tmp,(/nz,nprocz,3,njtest/),idir=12)
        uxbtestm=uxbtestm1_tmp(:,ipz+1,:,:)
      endif
!
!  Do the same for jxb; do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
      do jtest=1,njtest
        iaxtest=iaatest+3*(jtest-1)
        iaztest=iaxtest+2
        if (lsoca_jxb) then
          jxbtestm(:,:,jtest)=0.
        else
          do n=n1,n2
            nl=n-n1+1
            jxbtestm(nl,:,jtest)=0.
            do m=m1,m2
              aatest=f(l1:l2,m,n,iaxtest:iaztest)
              call gij(f,iaxtest,aijtest,1)
              call gij_etc(f,iaxtest,aatest,aijtest,bijtest,del2Atest2,graddivatest)
              call curl_mn(aijtest,bbtest,aatest)
              call curl_mn(bijtest,jjtest,bbtest)
              call cross_mn(jjtest,bbtest,jxbtest)
              jjxb=ijxbtest+3*(jtest-1)
              if (ijxbtest/=0) f(l1:l2,m,n,jjxb:jjxb+2)=jxbtest
              jxbtestm(nl,:,jtest)=jxbtestm(nl,:,jtest)+fac*sum(jxbtest,1)
              headtt=.false.
            enddo
            jxbtestm1(:,ipz+1,:,:)=jxbtestm
          enddo
        endif
      enddo
!
!  do communication for array of size nz*nprocz*3*njtest
!
      if (nprocxy>1) then
        call mpiallreduce_sum(jxbtestm1,jxbtestm1_tmp,(/nz,nprocz,3,njtest/),idir=12)
        jxbtestm=jxbtestm1_tmp(:,ipz+1,:,:)
      endif
!
!  calculate cosz*sinz, cos^2, and sinz^2, to take moments with
!  of alpij and etaij. This is useful if there is a mean Beltrami field
!  in the main calculation (lmagnetic=.true.) with phase zero.
!  Here we modify the calculations depending on the phase of the
!  actual field.
!
!  Calculate phase_testfield (for Beltrami fields)
!
      if (lphase_adjust) then
        if (lroot) then
          if (idiag_bcosphz/=0.and.idiag_bsinphz/=0) then
            bcosphz=fname(idiag_bcosphz)
            bsinphz=fname(idiag_bsinphz)
            phase_testfield=atan2(bsinphz,bcosphz)
          else
            call fatal_error('testfield_after_boundary', &
            'need bcosphz, bsinphz in print.in for lphase_adjust=T')
          endif
        endif
        call mpibcast_real(phase_testfield)
        c=cos(z+phase_testfield)
        s=sin(z+phase_testfield)
        c2z=c**2
        s2z=s**2
        csz=c*s
      endif
      if (ip<13) print*,'iproc,phase_testfield=',iproc,phase_testfield
!
!  reset headtt
!
      headtt=headtt_save
!
    endsubroutine testfield_after_boundary
!***********************************************************************
    subroutine rescaling_testfield(f)
!
!  Rescale testfield by factor rescale_aatest(jtest),
!  which could be different for different testfields
!
!  18-may-08/axel: rewrite from rescaling as used in magnetic
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=fnlen) :: file
      logical :: ltestfield_out
      logical, save :: lfirst_call=.true.
      integer :: j,jtest
!
      intent(inout) :: f
!
!  reinitialize aatest periodically if requested
!
      if (linit_aatest) then
        file=trim(datadir)//'/tinit_aatest.dat'
        if (lfirst_call) then
          call read_snaptime(trim(file),taainit,naainit,daainit,t)
          if (taainit==0 .or. taainit < t-daainit) then
            taainit=t+daainit
          endif
          lfirst_call=.false.
        endif
!
!  Do only one xy plane at a time (for cache efficiency)
!  Also: reset nuxb=0, which is used for time-averaged testfields
!  Do this for the full nuxb array.
!
        if (t >= taainit) then
          if (ltestfield_taver) nuxb=0
          do jtest=1,njtest
            iaxtest=iaatest+3*(jtest-1)
            iaztest=iaxtest+2
            do j=iaxtest,iaztest
              do n=n1,n2
                f(l1:l2,m1:m2,n,j)=rescale_aatest(jtest)*f(l1:l2,m1:m2,n,j)
              enddo
            enddo
          enddo
          call update_snaptime(file,taainit,naainit,daainit,t,ltestfield_out)
        endif
      endif
!
    endsubroutine rescaling_testfield
!***********************************************************************
    subroutine set_bbtest_linear (B0test,jtest)
!
!  set testfield
!
!  29-jun-10/axel: adapted from set_bbtest_linear_old
!
      use Cdata
!
      real, dimension (nx,3) :: B0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: B0test
!
!  set B0test for each of the 9 cases
!    B^T0 = (1, 1, 1) sx*sy
!    B^T1 = (z, 0, 0) sx*sy
!    B^T2 = (0, 0, z) sx*sy
!
      select case (jtest)
      case (1); B0test(:,1)=bamp*sx*sy(m); B0test(:,2)=B0test(:,1); B0test(:,3)=B0test(:,1)
      case (2); B0test(:,1)=bampz(n)*sx*sy(m); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bampz(n)*sx*sy(m)
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_linear
!***********************************************************************
    subroutine set_bbtest_linear_old (B0test,jtest)
!
!  set testfield
!
!  12-feb-10/axel: adapted from testfield_z
!  29-jun-10/axel: renamed to set_bbtest_linear_old
!
      use Cdata
!
      real, dimension (nx,3) :: B0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: B0test
!
!  set B0test for each of the 9 cases
!
      select case (jtest)
      case (1); B0test(:,1)=bamp; B0test(:,2)=bamp; B0test(:,3)=bamp
      case (2); B0test(:,1)=bamp*sy(m)*z(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*sy(m)*z(n)
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_linear_old
!***********************************************************************
    subroutine set_bbtest_sinkz (B0test,jtest)
!
!  set testfield
!
!  12-feb-10/axel: adapted from testfield_z
!
      use Cdata
!
      real, dimension (nx,3) :: B0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: B0test
!
!  set B0test for each of the 9 cases
!
      select case (jtest)
      case (1); B0test(:,1)=bamp*sy(m)*sz(n); B0test(:,2)=B0test(:,1); B0test(:,3)=B0test(:,1)
      case (2); B0test(:,1)=bamp*sy(m)*sz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*sy(m)*sz(n)
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_sinkz
!***********************************************************************
    subroutine set_bbtest_sxsysz (B0test,jtest)
!
!  set testfield
!
!  12-feb-10/axel: adapted from testfield_z
!
      use Cdata
!
      real, dimension (nx,3) :: B0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: B0test
!
!  set B0test for each of the 9 cases
!
      select case (jtest)
      case (1); B0test(:,1)=bamp*sx*sy(m)*sz(n); B0test(:,2)=B0test(:,1); B0test(:,3)=B0test(:,1)
      case (2); B0test(:,1)=bamp*sx*sy(m)*sz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*sx*sy(m)*sz(n)
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_sxsysz
!***********************************************************************
    subroutine rprint_testfield(lreset,lwrite)
!
!  reads and registers print parameters relevant for testfield fields
!
!   3-jun-05/axel: adapted from rprint_magnetic
!
      use Cdata
      use Diagnostics
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
        idiag_alpPERP=0; idiag_alpPARA=0; idiag_gam=0
        idiag_betPERP=0; idiag_betPARA=0; idiag_del=0
        idiag_kapPERP=0; idiag_kapPARA=0; idiag_mu=0
        idiag_alpPERPz=0; idiag_alpPARAz=0; idiag_gamz=0
        idiag_betPERPz=0; idiag_betPARAz=0; idiag_delz=0
        idiag_kapPERPz=0; idiag_kapPARAz=0; idiag_muz=0
        idiag_b1rms=0; idiag_b2rms=0; idiag_b3rms=0
        idiag_bx1pt=0; idiag_bx2pt=0; idiag_bx3pt=0
        ivid_bb1=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'alpPERP',idiag_alpPERP)
        call parse_name(iname,cname(iname),cform(iname),'alpPARA',idiag_alpPARA)
        call parse_name(iname,cname(iname),cform(iname),'gam',idiag_gam)
        call parse_name(iname,cname(iname),cform(iname),'betPERP',idiag_betPERP)
        call parse_name(iname,cname(iname),cform(iname),'betPARA',idiag_betPARA)
        call parse_name(iname,cname(iname),cform(iname),'del',idiag_del)
        call parse_name(iname,cname(iname),cform(iname),'kapPERP',idiag_kapPERP)
        call parse_name(iname,cname(iname),cform(iname),'kapPARA',idiag_kapPARA)
        call parse_name(iname,cname(iname),cform(iname),'mu',idiag_mu)
        call parse_name(iname,cname(iname),cform(iname),'bx1pt',idiag_bx1pt)
        call parse_name(iname,cname(iname),cform(iname),'bx2pt',idiag_bx2pt)
        call parse_name(iname,cname(iname),cform(iname),'bx3pt',idiag_bx3pt)
        call parse_name(iname,cname(iname),cform(iname),'b1rms',idiag_b1rms)
        call parse_name(iname,cname(iname),cform(iname),'b2rms',idiag_b2rms)
        call parse_name(iname,cname(iname),cform(iname),'b3rms',idiag_b3rms)
        call parse_name(iname,cname(iname),cform(iname),'bcosphz',idiag_bcosphz)
        call parse_name(iname,cname(iname),cform(iname),'bsinphz',idiag_bsinphz)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alpPERPz',idiag_alpPERPz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'alpPARAz',idiag_alpPARAz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'gamz',idiag_gamz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'betPERPz',idiag_betPERPz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'betPARAz',idiag_betPARAz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'delz',idiag_delz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kapPERPz',idiag_kapPERPz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'kapPARAz',idiag_kapPARAz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'muz',idiag_muz)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then 
        do iname=1,nnamev
          call parse_name(iname,cnamev(iname),cformv(iname),'bb1',ivid_bb1)
        enddo
      endif
!
    endsubroutine rprint_testfield
!**********************************************************************************
endmodule Testfield
