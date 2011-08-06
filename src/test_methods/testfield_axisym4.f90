! $Id$

!  Test-field module for axisymmetric turbulence.
!  No special settings in cparam.local are needed.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestfield = .true.
!
! MVAR CONTRIBUTION 12
! MAUX CONTRIBUTION 12
!
!***************************************************************

module Testfield

  use Cparam
  use Messages

  implicit none

  include '../testfield.h'
!
! Slice precalculation buffers
!
  real, target, dimension (nx,ny,3) :: bb1_xy
  real, target, dimension (nx,ny,3) :: bb1_xy2
  real, target, dimension (nx,nz,3) :: bb1_xz
  real, target, dimension (ny,nz,3) :: bb1_yz
!
!  cosine and sine function for setting test fields and analysis
!
  real, dimension(nx) :: sx,cx
  real, dimension(my) :: sy,cy
  real, dimension(mz) :: cz,sz,c2z,csz,s2z,c2kz,s2kz,kxBx1,kyBx1,kxBz1,kyBz1,zmask
  integer, parameter :: njtest=4
  real :: phase_testfield=.0
!
  character (len=labellen), dimension(ninit) :: initaatest='nothing'
  real, dimension (ninit) :: kx_aatest=1.,ky_aatest=1.,kz_aatest=1.
  real, dimension (ninit) :: phasex_aatest=0.,phasez_aatest=0.
  real, dimension (ninit) :: amplaatest=0.
  integer, dimension (njtest) :: nuxb=0

  ! input parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./)
  real, dimension(2) :: testfield_zaver_range=(/-max_real,max_real/)
  real :: taainit=0.,daainit=0.,taainit_previous=0.
  logical :: reinitialize_aatest=.false.
  logical :: zextent=.true.,lsoca=.false.,lsoca_jxb=.true.,lset_bbtest2=.false.
  logical :: luxb_as_aux=.false.,ljxb_as_aux=.false.,linit_aatest=.false.
  logical :: lignore_uxbtestm=.false., lphase_adjust=.false.
  character (len=labellen) :: itestfield='linear'
  real :: ktestfield=1., kxtestfield=1., kytestfield=1.
  real :: ktestfield1=1., kxtestfield1, kytestfield1, kx1, ky1, kz1
  real :: ztestfield_offset=0.
  real :: lin_testfield=0.,lam_testfield=0.,om_testfield=0.,delta_testfield=0.
  real :: delta_testfield_next=0., delta_testfield_time=0.
  integer, parameter :: mtestfield=3*njtest
  integer :: naainit
  real :: bamp=1.,bamp1=1.,bamp12=1.
  namelist /testfield_init_pars/ &
       B_ext,testfield_zaver_range,zextent,initaatest, &
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
       B_ext,testfield_zaver_range, &
       reinitialize_aatest,zextent,lsoca,lsoca_jxb, &
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
  integer :: idiag_betPERP2=0   ! DIAG_DOC: $\beta_\perp^{(2)}$
  integer :: idiag_betPARA=0    ! DIAG_DOC: $\beta_\perp$
  integer :: idiag_del=0        ! DIAG_DOC: $\delta$
  integer :: idiag_del2=0       ! DIAG_DOC: $\delta^{(2)}$
  integer :: idiag_kapPERP=0    ! DIAG_DOC: $\kappa_\perp$
  integer :: idiag_kapPERP2=0   ! DIAG_DOC: $\kappa_\perp^{(2)}$
  integer :: idiag_kapPARA=0    ! DIAG_DOC: $\kappa_\perp$
  integer :: idiag_mu=0         ! DIAG_DOC: $\mu$
  integer :: idiag_mu2=0        ! DIAG_DOC: $\mu^{(2)}$
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
!
!  arrays for horizontally averaged uxb and jxb
!
  real, dimension (mz,3,mtestfield/3) :: uxbtestm,jxbtestm

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
      use Mpicomm
      use Sub
!
      integer :: j
!
!  Set first and last index of text field
!  Note: iaxtest, iaytest, and iaztest are initialized to the first test field.
!  These values are used in this form in start, but later overwritten.
!
      iaatest=nvar+1
      iaxtest=iaatest
      iaytest=iaatest+1
      iaztest=iaatest+2
      ntestfield=mtestfield
      nvar=nvar+ntestfield
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_testfield: nvar = ', nvar
        print*, 'register_testfield: iaatest = ', iaatest
      endif
!
!  Put variable names in array
!
      do j=iaatest,nvar
        varname(j) = 'aatest'
      enddo
      iaztestpq=nvar
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           "$Id$")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('register_testfield: nvar > mvar')
      endif
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
    subroutine initialize_testfield(f,lstarting)
!
!  Perform any post-parameter-read initialization
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(mz) :: ztestfield, c, s
      real :: ktestfield_effective
      logical, intent(in) :: lstarting
      integer :: jtest
!
!  Precalculate etatest if 1/etatest (==etatest1) is given instead
!
      if (etatest1/=0.) then
        etatest=1./etatest1
      endif
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
!  define sinkx and coskx, called sx and cx
!
      sx=sin(kxtestfield*x(l1:l2))
      cx=cos(kxtestfield*x(l1:l2))
!
!  define sinky and cosky, called sy and cy
!
      sy=sin(kytestfield*y)
      cy=cos(kytestfield*y)
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
!  Compute mask for z-averaging where z is in testfield_zaver_range.
!  Normalize such that the average over the full domain
!  gives still unity.
!
      where (z>=testfield_zaver_range(1) .and. z<=testfield_zaver_range(2))
        zmask=1.
      elsewhere
        zmask=0.
      endwhere
      testfield_zaver_range(1)=max(testfield_zaver_range(1),xyz0(3))
      testfield_zaver_range(2)=min(testfield_zaver_range(2),xyz1(3))
      zmask=zmask*Lxyz(3)/(testfield_zaver_range(2)-testfield_zaver_range(1))
!
!  debug output
!
      if (lroot.and.ip<14) then
        print*,'cz=',cz
        print*,'sz=',sz
        print*,'zmask=',zmask
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
!  Calculate inverse kxtestfield, but only if different from zero
!  Introduce kxBz1 as abbreviation of 1/(kxtestfield*Btestz)
!
      if (kxtestfield==0) then
        kxtestfield1=1.
      else
        kxtestfield1=1./kxtestfield
      endif
      kxBx1=kxtestfield1/(z(n)+ztestfield_offset)
      kxBz1=kxtestfield1/(z(n)+ztestfield_offset)
      kx1=kxtestfield1
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
!     ! MAUX CONTRIBUTION 12
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
!
!  After a reload, we need to rewrite index.pro, but the auxiliary
!  arrays are already allocated and must not be allocated again.
!
      if (luxb_as_aux) then
        if (iuxb==0) then
          call farray_register_auxiliary('uxb',iuxb,vector=3*njtest)
        endif
        if (iuxb/=0.and.lroot) then
          print*, 'initialize_magnetic: iuxb = ', iuxb
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'iuxb=',iuxb
          close(3)
        endif
      endif
!
!  possibility of using jxb as auxiliary array (is intended to be
!  used in connection with testflow method)
!
      if (ljxb_as_aux) then
        if (ijxb==0) then
          call farray_register_auxiliary('jxb',ijxb,vector=3*njtest)
        endif
        if (ijxb/=0.and.lroot) then
          print*, 'initialize_magnetic: ijxb = ', ijxb
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'ijxb=',ijxb
          close(3)
        endif
      endif
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

      case ('zero'); f(:,:,:,iaatest:iaatest+ntestfield-1)=0.
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
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'init_aatest: check initaatest: ', trim(initaatest(j))
        call stop_it("")

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
      use Sub, only: keep_compiler_quiet
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_testfield
!***********************************************************************
    subroutine read_testfield_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=testfield_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testfield_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_testfield_init_pars
!***********************************************************************
    subroutine write_testfield_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=testfield_init_pars)

    endsubroutine write_testfield_init_pars
!***********************************************************************
    subroutine read_testfield_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=testfield_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testfield_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_testfield_run_pars
!***********************************************************************
    subroutine write_testfield_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=testfield_run_pars)

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
      use Cdata
      use Diagnostics
      use Hydro, only: uumz,lcalc_uumean
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx,3) :: uxB,B0test=0,bbtest
      real, dimension (nx,3) :: uxbtest,duxbtest,jxbtest,djxbrtest,eetest
      real, dimension (nx,3) :: J0test=0,jxB0rtest,J0xbrtest
      real, dimension (nx,3,3,njtest) :: Mijpq
      real, dimension (nx,3,njtest) :: Eipq,bpq,jpq
      real, dimension (nx,3) :: del2Atest,uufluct
      real, dimension (nx,3) :: del2Atest2,graddivatest,aatest,jjtest,jxbrtest
      real, dimension (nx,3,3) :: aijtest,bijtest,Mijtest
      real, dimension (nx) :: jbpq,bpq2,Epq2,s2kzDF1,s2kzDF2,unity=1.
      integer :: jtest,j, i1=1, i2=2, i3=3, i4=4, iuxtest, iuytest, iuztest
      logical,save :: ltest_uxb=.false.,ltest_jxb=.false.
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
      if (lcalc_uumean) then
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
      do jtest=1,njtest
        iaxtest=iaatest+3*(jtest-1)
        iaztest=iaxtest+2
        call del2v(f,iaxtest,del2Atest)
        select case (itestfield)
          case ('linear'); call set_bbtest_linear(B0test,jtest)
          case ('sxsysz'); call set_bbtest_sxsysz(B0test,jtest)
          case ('sinkz'); call set_bbtest_sinkz(B0test,jtest)
          case ('B=0') !(dont do anything)
        case default
          call fatal_error('daatest_dt','undefined itestfield value')
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
          df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
            +uxB+etatest*del2Atest
        else
!
!  use f-array for uxb (if space has been allocated for this) and
!  if we don't test (i.e. if ltest_uxb=.false.)
!
          if (iuxb/=0.and..not.ltest_uxb) then
            uxbtest=f(l1:l2,m,n,iuxb+3*(jtest-1):iuxb+3*jtest-1)
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
            duxbtest(:,:)=uxbtest(:,:)
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
              +ampl_fcont_aatest*p%fcont
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
        diffus_eta=max(diffus_eta,etatest*dxyz_2)
      endif
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
!  test fields linear in z
!
        case ('linear')
          call fatal_error('testfield_axisym4','dont use this part')
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
          call fatal_error('testfield_axisym4','dont use this part')
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
          if (idiag_gam    /=0) call sum_mn_name(-4*sx &
            *sy(m)*(cz(n)*Eipq(:,2,1)+sz(n)*Eipq(:,2,2)) &
            *zmask(n),idiag_gam    )
          if (idiag_alpPERP/=0) call sum_mn_name(-4*sx &
            *sy(m)*(cz(n)*Eipq(:,1,1)+sz(n)*Eipq(:,1,2)) &
            *zmask(n),idiag_alpPERP)
          if (idiag_alpPARA/=0) call sum_mn_name(-4*sx &
            *sy(m)*(cz(n)*Eipq(:,3,3)+sz(n)*Eipq(:,3,4)) &
            *zmask(n),idiag_alpPARA)
!
          if (idiag_mu     /=0) call sum_mn_name(+4*sx*( &
            -kz1*sy(m)*(-sz(n)*Eipq(:,2,1)+cz(n)*Eipq(:,2,2)) &
            +ky1*cy(m)*(+cz(n)*Eipq(:,1,3)+sz(n)*Eipq(:,1,4)) &
            )*zmask(n),idiag_mu)
          if (idiag_mu2    /=0) call sum_mn_name(+4*( &
            -kz1*sx*sy(m)*(-sz(n)*Eipq(:,2,1)+cz(n)*Eipq(:,2,2)) &
            -kx1*cx*sy(m)*(+cz(n)*Eipq(:,2,3)+sz(n)*Eipq(:,2,4)) &
            )*zmask(n),idiag_mu2)
          if (idiag_betPERP/=0) call sum_mn_name(-2*sx*( &
            +kz1*sy(m)*(-sz(n)*Eipq(:,2,1)+cz(n)*Eipq(:,2,2)) &
            +ky1*cy(m)*(+cz(n)*Eipq(:,1,3)+sz(n)*Eipq(:,1,4)) &
            )*zmask(n),idiag_betPERP)
          if (idiag_betPERP2/=0) call sum_mn_name(-2*( &
            +kz1*sx*sy(m)*(-sz(n)*Eipq(:,2,1)+cz(n)*Eipq(:,2,2)) &
            -kx1*cx*sy(m)*(+cz(n)*Eipq(:,2,3)+sz(n)*Eipq(:,2,4)) &
            )*zmask(n),idiag_betPERP2)
!
          if (idiag_del    /=0) call sum_mn_name(+2*sx*( &
            +kz1*sy(m)*(-sz(n)*Eipq(:,1,1)+cz(n)*Eipq(:,1,2)) &
            -ky1*cy(m)*(+cz(n)*Eipq(:,2,3)+sz(n)*Eipq(:,2,4)) &
            )*zmask(n),idiag_del)
          if (idiag_del2   /=0) call sum_mn_name(+2*( &
            +kz1*sx*sy(m)*(-sz(n)*Eipq(:,1,1)+cz(n)*Eipq(:,1,2)) &
            -kx1*cx*sy(m)*(+cz(n)*Eipq(:,1,3)+sz(n)*Eipq(:,1,4)) &
            )*zmask(n),idiag_del2)
          if (idiag_kapPERP/=0) call sum_mn_name(-4*sx*( &
            +kz1*sy(m)*(-sz(n)*Eipq(:,1,1)+cz(n)*Eipq(:,1,2)) &
            +ky1*cy(m)*(+cz(n)*Eipq(:,2,3)+sz(n)*Eipq(:,2,4)) &
            )*zmask(n),idiag_kapPERP)
          if (idiag_kapPERP2/=0) call sum_mn_name(-4*( &
            +kz1*sx*sy(m)*(-sz(n)*Eipq(:,1,1)+cz(n)*Eipq(:,1,2)) &
            +kx1*cx*sy(m)*(+cz(n)*Eipq(:,1,3)+sz(n)*Eipq(:,1,4)) &
            )*zmask(n),idiag_kapPERP2)
!
          if (idiag_betPARA/=0) call sum_mn_name(+4*sx* &
             ky1*cy(m)*(+cz(n)*Eipq(:,3,1)+sz(n)*Eipq(:,3,2)) &
             *zmask(n),idiag_betPARA)
          if (idiag_kapPARA/=0) call sum_mn_name(-4*sx* &
             kz1*sy(m)*(-sz(n)*Eipq(:,3,3)+cz(n)*Eipq(:,3,4)) &
             *zmask(n),idiag_kapPARA)
!
!  xy-averages next.
!
          if (idiag_gamz    /=0) call xysum_mn_name_z(-4*sx*sy(m)*(cz(n)*Eipq(:,2,1)+sz(n)*Eipq(:,2,2)),idiag_gamz    )
          if (idiag_alpPERPz/=0) call xysum_mn_name_z(-4*sx*sy(m)*(cz(n)*Eipq(:,1,1)+sz(n)*Eipq(:,1,2)),idiag_alpPERPz)
          if (idiag_alpPARAz/=0) call xysum_mn_name_z(-4*sx*sy(m)*(cz(n)*Eipq(:,3,3)+sz(n)*Eipq(:,3,4)),idiag_alpPARAz)
!
          if (idiag_muz     /=0) call xysum_mn_name_z(+4*sx*(-kz1*sy(m)*(-sz(n)*Eipq(:,2,1)+cz(n)*Eipq(:,2,2)) &
                                                             +ky1*cy(m)*(+cz(n)*Eipq(:,1,3)+sz(n)*Eipq(:,1,4))),idiag_muz)
          if (idiag_betPERPz/=0) call xysum_mn_name_z(-2*sx*(+kz1*sy(m)*(-sz(n)*Eipq(:,2,1)+cz(n)*Eipq(:,2,2)) &
                                                             +ky1*cy(m)*(+cz(n)*Eipq(:,1,3)+sz(n)*Eipq(:,1,4))),idiag_betPERPz)
!
          if (idiag_delz    /=0) call xysum_mn_name_z(+2*sx*(+kz1*sy(m)*(-sz(n)*Eipq(:,1,1)+cz(n)*Eipq(:,1,2)) &
                                                             -ky1*cy(m)*(+cz(n)*Eipq(:,2,3)+sz(n)*Eipq(:,2,4))),idiag_delz)
          if (idiag_kapPERPz/=0) call xysum_mn_name_z(-4*sx*(+kz1*sy(m)*(-sz(n)*Eipq(:,1,1)+cz(n)*Eipq(:,1,2)) &
                                                             +ky1*cy(m)*(+cz(n)*Eipq(:,2,3)+sz(n)*Eipq(:,2,4))),idiag_kapPERPz)
!
          if (idiag_betPARAz/=0) call xysum_mn_name_z(+4*sx*  ky1*cy(m)*(+cz(n)*Eipq(:,3,1)+sz(n)*Eipq(:,3,2)) ,idiag_betPARAz)
          if (idiag_kapPARAz/=0) call xysum_mn_name_z(-4*sx*  kz1*sy(m)*(-sz(n)*Eipq(:,3,3)+cz(n)*Eipq(:,3,4)) ,idiag_kapPARAz)
!
        case default
          call fatal_error('daatest_dt','undefined itestfield value')
        endselect
!
!  diagnostics for single points
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_bx1pt/=0) call save_name(bpq(lpoint-nghost,1,i1),idiag_bx1pt)
          if (idiag_bx2pt/=0) call save_name(bpq(lpoint-nghost,1,i2),idiag_bx2pt)
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
!  Note: ix is the index with respect to array with ghost zones.
! 
      if (lvideo.and.lfirst) then
        do j=1,3
          bb1_yz(m-m1+1,n-n1+1,j)=bpq(ix_loc-l1+1,j,1)
          if (m==iy_loc)  bb1_xz(:,n-n1+1,j)=bpq(:,j,1)
          if (n==iz_loc)  bb1_xy(:,m-m1+1,j)=bpq(:,j,1)
          if (n==iz2_loc) bb1_xy2(:,m-m1+1,j)=bpq(:,j,1)
        enddo
      endif
!
    endsubroutine daatest_dt
!***********************************************************************
    subroutine get_slices_testfield(f,slices)
! 
!  Write slices for animation of magnetic variables.
! 
!  12-sep-09/axel: adapted from the corresponding magnetic routine
! 
      use Sub, only: keep_compiler_quiet
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
        if (slices%index>=3) then
          slices%ready=.false.
        else
          slices%index=slices%index+1
          slices%yz =>bb1_yz(:,:,slices%index)
          slices%xz =>bb1_xz(:,:,slices%index)
          slices%xy =>bb1_xy(:,:,slices%index)
          slices%xy2=>bb1_xy2(:,:,slices%index)
          if (slices%index<=3) slices%ready=.true.
        endif
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine get_slices_testfield
!***********************************************************************
    subroutine testfield_after_boundary(f,p)
!
!  calculate <uxb>, which is needed when lsoca=.false.
!
!  21-jan-06/axel: coded
!
      use Cdata
      use Sub
      use Hydro, only: calc_pencils_hydro
      use Magnetic, only: idiag_bcosphz, idiag_bsinphz
      use Mpicomm, only: mpireduce_sum, mpibcast_real, mpibcast_real_arr
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension (mz) :: c,s
!
      real, dimension (nz,nprocz,3,njtest) :: uxbtestm1=0.,uxbtestm1_tmp=0.
      real, dimension (nz,nprocz,3,njtest) :: jxbtestm1=0.,jxbtestm1_tmp=0.
!
      real, dimension (nx,3,3) :: aijtest,bijtest
      real, dimension (nx,3) :: aatest,bbtest,jjtest,uxbtest,jxbtest
      real, dimension (nx,3) :: del2Atest2,graddivatest
      integer :: jtest,j,nxy=nxgrid*nygrid,juxb,jjxb
      logical :: headtt_save
      real :: fac, bcosphz, bsinphz, fac1=0., fac2=1.
!
      intent(inout) :: f
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac=1./nxy
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
            uxbtestm(n,:,jtest)=0.
            do m=m1,m2
              aatest=f(l1:l2,m,n,iaxtest:iaztest)
              call calc_pencils_hydro(f,p)
              call gij(f,iaxtest,aijtest,1)
              call curl_mn(aijtest,bbtest,aatest)
              call cross_mn(p%uu,bbtest,uxbtest)
              juxb=iuxb+3*(jtest-1)
              if (ltestfield_taver) then
                if (llast) then
                  if (iuxb/=0) f(l1:l2,m,n,juxb:juxb+2)= &
                          fac1*f(l1:l2,m,n,juxb:juxb+2)+fac2*uxbtest
                endif
              else
                if (iuxb/=0) f(l1:l2,m,n,juxb:juxb+2)=uxbtest
              endif
              do j=1,3
                uxbtestm(n,j,jtest)=uxbtestm(n,j,jtest)+fac*sum(uxbtest(:,j))
              enddo
              headtt=.false.
            enddo
            do j=1,3
              uxbtestm1(n-n1+1,ipz+1,j,jtest)=uxbtestm(n,j,jtest)
            enddo
          enddo
!
!  Update counter, but only when we are on the last substep.
!
          if (ltestfield_taver) then
            if (llast) then
              nuxb(jtest)=nuxb(jtest)+1
            endif
          endif
        endif
      enddo
!
!  Do communication for array of size nz*nprocz*3*njtest.
!  The factor 3 is because of the 3 components of each vector field.
!
      if (nprocy>1) then
        call mpireduce_sum(uxbtestm1,uxbtestm1_tmp,(/nz,nprocz,3,njtest/))
        call mpibcast_real_arr(uxbtestm1_tmp,nz*nprocz*3*njtest)
        do jtest=1,njtest
          do n=n1,n2
            do j=1,3
              uxbtestm(n,j,jtest)=uxbtestm1_tmp(n-n1+1,ipz+1,j,jtest)
            enddo
          enddo
        enddo
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
            jxbtestm(n,:,jtest)=0.
            do m=m1,m2
              aatest=f(l1:l2,m,n,iaxtest:iaztest)
              call gij(f,iaxtest,aijtest,1)
              call gij_etc(f,iaxtest,aatest,aijtest,bijtest,del2Atest2,graddivatest)
              call curl_mn(aijtest,bbtest,aatest)
              call curl_mn(bijtest,jjtest,bbtest)
              call cross_mn(jjtest,bbtest,jxbtest)
              jjxb=ijxb+3*(jtest-1)
              if (ijxb/=0) f(l1:l2,m,n,jjxb:jjxb+2)=jxbtest
              do j=1,3
                jxbtestm(n,j,jtest)=jxbtestm(n,j,jtest)+fac*sum(jxbtest(:,j))
              enddo
              headtt=.false.
            enddo
            do j=1,3
              jxbtestm1(n-n1+1,ipz+1,j,jtest)=jxbtestm(n,j,jtest)
            enddo
          enddo
        endif
      enddo
!
!  do communication for array of size nz*nprocz*3*njtest
!
      if (nprocy>1) then
        call mpireduce_sum(jxbtestm1,jxbtestm1_tmp,(/nz,nprocz,3,njtest/))
        call mpibcast_real_arr(jxbtestm1_tmp,nz*nprocz*3*njtest)
        do jtest=1,njtest
          do n=n1,n2
            do j=1,3
              jxbtestm(n,j,jtest)=jxbtestm1_tmp(n-n1+1,ipz+1,j,jtest)
            enddo
          enddo
        enddo
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
        call mpibcast_real(phase_testfield,1)
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
      integer,save :: ifirst=0
      integer :: j,jtest
!
      intent(inout) :: f
!
!  reinitialize aatest periodically if requested
!
      if (linit_aatest) then
        file=trim(datadir)//'/tinit_aatest.dat'
        if (ifirst==0) then
          call read_snaptime(trim(file),taainit,naainit,daainit,t)
          if (taainit==0 .or. taainit < t-daainit) then
            taainit=t+daainit
          endif
          ifirst=1
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
      case (1); B0test(:,1)=bamp; B0test(:,2)=bamp; B0test(:,3)=bamp
      case (2); B0test(:,1)=bamp*sy(m)*z(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*sy(m)*z(n)
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_linear
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
!  16-jun-11/axel: adapted from testfield_axisym2
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
      case (1); B0test(:,1)=bamp*sx*sy(m)*cz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=bamp*sx*sy(m)*sz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*sx*sy(m)*cz(n)
      case (4); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*sx*sy(m)*sz(n)
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
!
      integer :: iname,inamez
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
        idiag_alpPERP=0; idiag_alpPARA=0; idiag_gam=0
        idiag_betPERP=0; idiag_betPERP2=0; idiag_betPARA=0; idiag_del=0; idiag_del2=0
        idiag_kapPERP=0; idiag_kapPERP2=0; idiag_kapPARA=0; idiag_mu=0; idiag_mu2=0
        idiag_alpPERPz=0; idiag_alpPARAz=0; idiag_gamz=0
        idiag_betPERPz=0; idiag_betPARAz=0; idiag_delz=0
        idiag_kapPERPz=0; idiag_kapPARAz=0; idiag_muz=0
        idiag_b1rms=0; idiag_b2rms=0; idiag_b3rms=0
        idiag_bx1pt=0; idiag_bx2pt=0; idiag_bx3pt=0
      endif
!
!  check for those quantities that we want to evaluate online
! 
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'alpPERP',idiag_alpPERP)
        call parse_name(iname,cname(iname),cform(iname),'alpPARA',idiag_alpPARA)
        call parse_name(iname,cname(iname),cform(iname),'gam',idiag_gam)
        call parse_name(iname,cname(iname),cform(iname),'betPERP',idiag_betPERP)
        call parse_name(iname,cname(iname),cform(iname),'betPERP2',idiag_betPERP2)
        call parse_name(iname,cname(iname),cform(iname),'betPARA',idiag_betPARA)
        call parse_name(iname,cname(iname),cform(iname),'del',idiag_del)
        call parse_name(iname,cname(iname),cform(iname),'del2',idiag_del2)
        call parse_name(iname,cname(iname),cform(iname),'kapPERP',idiag_kapPERP)
        call parse_name(iname,cname(iname),cform(iname),'kapPERP2',idiag_kapPERP2)
        call parse_name(iname,cname(iname),cform(iname),'kapPARA',idiag_kapPARA)
        call parse_name(iname,cname(iname),cform(iname),'mu',idiag_mu)
        call parse_name(iname,cname(iname),cform(iname),'mu2',idiag_mu2)
        call parse_name(iname,cname(iname),cform(iname),'bx1pt',idiag_bx1pt)
        call parse_name(iname,cname(iname),cform(iname),'bx2pt',idiag_bx2pt)
        call parse_name(iname,cname(iname),cform(iname),'bx3pt',idiag_bx3pt)
        call parse_name(iname,cname(iname),cform(iname),'b1rms',idiag_b1rms)
        call parse_name(iname,cname(iname),cform(iname),'b2rms',idiag_b2rms)
        call parse_name(iname,cname(iname),cform(iname),'b3rms',idiag_b3rms)
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

!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'iaatest=',iaatest
        write(3,*) 'ntestfield=',ntestfield
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'nnamexz=',nnamexz
      endif
!
    endsubroutine rprint_testfield

endmodule Testfield
