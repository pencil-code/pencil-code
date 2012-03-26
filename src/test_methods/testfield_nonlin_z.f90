! $Id$

!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.

!  Note: this routine requires that MVAR and MAUX contributions
!  together with njtest are set correctly in the cparam.local file.
!  njtest must be set at the end of the file such that 6*njtest=MVAR.
!
!  Example:
!  ! MVAR CONTRIBUTION 24
!  ! MAUX CONTRIBUTION 24
!  integer, parameter :: njtest=5

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestfield = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
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
  real, target, dimension (nx,ny,3) :: bb11_xy
  real, target, dimension (nx,ny,3) :: bb11_xy2
  real, target, dimension (nx,nz,3) :: bb11_xz
  real, target, dimension (ny,nz,3) :: bb11_yz
!
!  cosine and sine function for setting test fields and analysis
!
  real, dimension(mz) :: cz,sz,c2z,csz,s2z,c2kz,s2kz
  real :: phase_testfield=.0
!
  character (len=labellen), dimension(ninit) :: initaatest='nothing'
  real, dimension (ninit) :: kx_aatest=1.,ky_aatest=1.,kz_aatest=1.
  real, dimension (ninit) :: phasex_aatest=0.,phasez_aatest=0.
  real, dimension (ninit) :: amplaatest=0.
  integer :: iuxtest=0,iuytest=0,iuztest=0,iuztestpq=0
  integer :: iu0xtest=0,iu0ztest=0
  integer, dimension (njtest) :: nuxb=0
  integer :: iE0=0

  ! input parameters
  real, dimension(3) :: Btest_ext=(/0.,0.,0./)
  real :: taainit=0.,daainit=0.,taainit_previous=0.
  logical :: reinitialize_aatest=.false.
  logical :: reinitialize_from_mainrun=.false.
  logical :: zextent=.true.,lsoca=.false.,lsoca_jxb=.true.,lset_bbtest2=.false.
  logical :: luxb_as_aux=.false.,ljxb_as_aux=.false.,linit_aatest=.false.
  logical :: lignore_uxbtestm=.false., lignore_jxbtestm=.false., lphase_adjust=.false.
  character (len=labellen) :: itestfield='B11-B21',itestfield_method='(i)'
  real :: ktestfield=1., ktestfield1=1.
  real :: lin_testfield=0.,lam_testfield=0.,om_testfield=0.,delta_testfield=0.
  real :: delta_testfield_next=0., delta_testfield_time=0.
  integer, parameter :: mtestfield=6*njtest
  integer :: naainit
  real :: bamp=1.,bamp1=1.,bamp12=1.
  namelist /testfield_init_pars/ &
       Btest_ext,zextent,initaatest, &
       amplaatest,kx_aatest,ky_aatest,kz_aatest, &
       phasex_aatest,phasez_aatest, &
       luxb_as_aux,ljxb_as_aux

  ! run parameters
  real :: etatest=0.,etatest1=0.,nutest=0.,nutest1=0.
  real :: ampl_fcont_aatest=1.,ampl_fcont_uutest=1.
  real, dimension(njtest) :: rescale_aatest=0.,rescale_uutest=0.
  logical :: ltestfield_newz=.true.,leta_rank2=.true.
  logical :: lforcing_cont_aatest=.false.,lforcing_cont_uutest=.false.
  namelist /testfield_run_pars/ &
       reinitialize_aatest,reinitialize_from_mainrun, &
       Btest_ext,zextent,lsoca,lsoca_jxb, &
       lset_bbtest2,itestfield,ktestfield,itestfield_method, &
       etatest,etatest1,nutest,nutest1, &
       lin_testfield,lam_testfield,om_testfield,delta_testfield, &
       ltestfield_newz,leta_rank2,lphase_adjust,phase_testfield, &
       luxb_as_aux,ljxb_as_aux,lignore_uxbtestm,lignore_jxbtestm, &
       lforcing_cont_aatest,ampl_fcont_aatest, &
       lforcing_cont_uutest,ampl_fcont_uutest, &
       daainit,linit_aatest,bamp, &
       rescale_aatest,rescale_uutest

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_alp11=0      ! DIAG_DOC: $\alpha_{11}$
  integer :: idiag_alp21=0      ! DIAG_DOC: $\alpha_{21}$
  integer :: idiag_alp31=0      ! DIAG_DOC: $\alpha_{31}$
  integer :: idiag_alp12=0      ! DIAG_DOC: $\alpha_{12}$
  integer :: idiag_alp22=0      ! DIAG_DOC: $\alpha_{22}$
  integer :: idiag_alp32=0      ! DIAG_DOC: $\alpha_{32}$
  integer :: idiag_eta11=0      ! DIAG_DOC: $\eta_{11}k$
  integer :: idiag_eta21=0      ! DIAG_DOC: $\eta_{21}k$
  integer :: idiag_eta12=0      ! DIAG_DOC: $\eta_{12}k$
  integer :: idiag_eta22=0      ! DIAG_DOC: $\eta_{22}k$
  integer :: idiag_alpK=0       ! DIAG_DOC: $\alpha^K$
  integer :: idiag_alpM=0       ! DIAG_DOC: $\alpha^M$
  integer :: idiag_alpMK=0      ! DIAG_DOC: $\alpha^{MK}$
  integer :: idiag_phi11=0      ! DIAG_DOC: $\phi_{11}$
  integer :: idiag_phi21=0      ! DIAG_DOC: $\phi_{21}$
  integer :: idiag_phi12=0      ! DIAG_DOC: $\phi_{12}$
  integer :: idiag_phi22=0      ! DIAG_DOC: $\phi_{22}$
  integer :: idiag_phi32=0      ! DIAG_DOC: $\phi_{32}$
  integer :: idiag_psi11=0      ! DIAG_DOC: $\psi_{11}k$
  integer :: idiag_psi21=0      ! DIAG_DOC: $\psi_{21}k$
  integer :: idiag_psi12=0      ! DIAG_DOC: $\psi_{12}k$
  integer :: idiag_psi22=0      ! DIAG_DOC: $\psi_{22}k$
  integer :: idiag_phiK=0       ! DIAG_DOC: $\phi^K$
  integer :: idiag_phiM=0       ! DIAG_DOC: $\phi^M$
  integer :: idiag_phiMK=0      ! DIAG_DOC: $\phi^{MK}$
  integer :: idiag_alp11cc=0    ! DIAG_DOC: $\alpha_{11}\cos^2 kz$
  integer :: idiag_alp21sc=0    ! DIAG_DOC: $\alpha_{21}\sin kz\cos kz$
  integer :: idiag_alp12cs=0    ! DIAG_DOC: $\alpha_{12}\cos kz\sin kz$
  integer :: idiag_alp22ss=0    ! DIAG_DOC: $\alpha_{22}\sin^2 kz$
  integer :: idiag_eta11cc=0    ! DIAG_DOC: $\eta_{11}\cos^2 kz$
  integer :: idiag_eta21sc=0    ! DIAG_DOC: $\eta_{21}\sin kz\cos kz$
  integer :: idiag_eta12cs=0    ! DIAG_DOC: $\eta_{12}\cos kz\sin kz$
  integer :: idiag_eta22ss=0    ! DIAG_DOC: $\eta_{22}\sin^2 kz$
  integer :: idiag_s2kzDFm=0    ! DIAG_DOC: $\left<\sin2kz\nabla\cdot F\right>$
  integer :: idiag_M11=0        ! DIAG_DOC: ${\cal M}_{11}$
  integer :: idiag_M22=0        ! DIAG_DOC: ${\cal M}_{22}$
  integer :: idiag_M33=0        ! DIAG_DOC: ${\cal M}_{33}$
  integer :: idiag_M11cc=0      ! DIAG_DOC: ${\cal M}_{11}\cos^2 kz$
  integer :: idiag_M11ss=0      ! DIAG_DOC: ${\cal M}_{11}\sin^2 kz$
  integer :: idiag_M22cc=0      ! DIAG_DOC: ${\cal M}_{22}\cos^2 kz$
  integer :: idiag_M22ss=0      ! DIAG_DOC: ${\cal M}_{22}\sin^2 kz$
  integer :: idiag_M12cs=0      ! DIAG_DOC: ${\cal M}_{12}\cos kz\sin kz$
  integer :: idiag_bx11pt=0     ! DIAG_DOC: $b_x^{11}$
  integer :: idiag_bx21pt=0     ! DIAG_DOC: $b_x^{21}$
  integer :: idiag_bx12pt=0     ! DIAG_DOC: $b_x^{12}$
  integer :: idiag_bx22pt=0     ! DIAG_DOC: $b_x^{22}$
  integer :: idiag_bx0pt=0      ! DIAG_DOC: $b_x^{0}$
  integer :: idiag_by11pt=0     ! DIAG_DOC: $b_y^{11}$
  integer :: idiag_by21pt=0     ! DIAG_DOC: $b_y^{21}$
  integer :: idiag_by12pt=0     ! DIAG_DOC: $b_y^{12}$
  integer :: idiag_by22pt=0     ! DIAG_DOC: $b_y^{22}$
  integer :: idiag_by0pt=0      ! DIAG_DOC: $b_y^{0}$
  integer :: idiag_u11rms=0     ! DIAG_DOC: $\left<u_{11}^2\right>^{1/2}$
  integer :: idiag_u21rms=0     ! DIAG_DOC: $\left<u_{21}^2\right>^{1/2}$
  integer :: idiag_u12rms=0     ! DIAG_DOC: $\left<u_{12}^2\right>^{1/2}$
  integer :: idiag_u22rms=0     ! DIAG_DOC: $\left<u_{22}^2\right>^{1/2}$
  integer :: idiag_j11rms=0     ! DIAG_DOC: $\left<j_{11}^2\right>^{1/2}$
  integer :: idiag_b11rms=0     ! DIAG_DOC: $\left<b_{11}^2\right>^{1/2}$
  integer :: idiag_b21rms=0     ! DIAG_DOC: $\left<b_{21}^2\right>^{1/2}$
  integer :: idiag_b12rms=0     ! DIAG_DOC: $\left<b_{12}^2\right>^{1/2}$
  integer :: idiag_b22rms=0     ! DIAG_DOC: $\left<b_{22}^2\right>^{1/2}$
  integer :: idiag_ux0m=0       ! DIAG_DOC: $\left<u_{0_x}\right>$
  integer :: idiag_uy0m=0       ! DIAG_DOC: $\left<u_{0_y}\right>$
  integer :: idiag_ux11m=0      ! DIAG_DOC: $\left<u_{11_x}\right>$
  integer :: idiag_uy11m=0      ! DIAG_DOC: $\left<u_{11_y}\right>$
  integer :: idiag_u0rms=0      ! DIAG_DOC: $\left<u_{0}^2\right>^{1/2}$
  integer :: idiag_b0rms=0      ! DIAG_DOC: $\left<b_{0}^2\right>^{1/2}$
  integer :: idiag_jb0m=0       ! DIAG_DOC: $\left<jb_{0}\right>$
  integer :: idiag_E11rms=0     ! DIAG_DOC: $\left<{\cal E}_{11}^2\right>^{1/2}$
  integer :: idiag_E21rms=0     ! DIAG_DOC: $\left<{\cal E}_{21}^2\right>^{1/2}$
  integer :: idiag_E12rms=0     ! DIAG_DOC: $\left<{\cal E}_{12}^2\right>^{1/2}$
  integer :: idiag_E22rms=0     ! DIAG_DOC: $\left<{\cal E}_{22}^2\right>^{1/2}$
  integer :: idiag_E0rms=0      ! DIAG_DOC: $\left<{\cal E}_{0}^2\right>^{1/2}$
  integer :: idiag_Ex11pt=0     ! DIAG_DOC: ${\cal E}_x^{11}$
  integer :: idiag_Ex21pt=0     ! DIAG_DOC: ${\cal E}_x^{21}$
  integer :: idiag_Ex12pt=0     ! DIAG_DOC: ${\cal E}_x^{12}$
  integer :: idiag_Ex22pt=0     ! DIAG_DOC: ${\cal E}_x^{22}$
  integer :: idiag_Ex0pt=0      ! DIAG_DOC: ${\cal E}_x^{0}$
  integer :: idiag_Ey11pt=0     ! DIAG_DOC: ${\cal E}_y^{11}$
  integer :: idiag_Ey21pt=0     ! DIAG_DOC: ${\cal E}_y^{21}$
  integer :: idiag_Ey12pt=0     ! DIAG_DOC: ${\cal E}_y^{12}$
  integer :: idiag_Ey22pt=0     ! DIAG_DOC: ${\cal E}_y^{22}$
  integer :: idiag_Ey0pt=0      ! DIAG_DOC: ${\cal E}_y^{0}$
  integer :: idiag_bamp=0       ! DIAG_DOC: bamp
  integer :: idiag_E111z=0      ! DIAG_DOC: ${\cal E}_1^{11}$
  integer :: idiag_E211z=0      ! DIAG_DOC: ${\cal E}_2^{11}$
  integer :: idiag_E311z=0      ! DIAG_DOC: ${\cal E}_3^{11}$
  integer :: idiag_E121z=0      ! DIAG_DOC: ${\cal E}_1^{21}$
  integer :: idiag_E221z=0      ! DIAG_DOC: ${\cal E}_2^{21}$
  integer :: idiag_E321z=0      ! DIAG_DOC: ${\cal E}_3^{21}$
  integer :: idiag_E112z=0      ! DIAG_DOC: ${\cal E}_1^{12}$
  integer :: idiag_E212z=0      ! DIAG_DOC: ${\cal E}_2^{12}$
  integer :: idiag_E312z=0      ! DIAG_DOC: ${\cal E}_3^{12}$
  integer :: idiag_E122z=0      ! DIAG_DOC: ${\cal E}_1^{22}$
  integer :: idiag_E222z=0      ! DIAG_DOC: ${\cal E}_2^{22}$
  integer :: idiag_E322z=0      ! DIAG_DOC: ${\cal E}_3^{22}$
  integer :: idiag_E10z=0       ! DIAG_DOC: ${\cal E}_1^{0}$
  integer :: idiag_E20z=0       ! DIAG_DOC: ${\cal E}_2^{0}$
  integer :: idiag_E30z=0       ! DIAG_DOC: ${\cal E}_3^{0}$
  integer :: idiag_EBpq=0       ! DIAG_DOC: ${\cal E}\cdot\Bv^{pq}$
  integer :: idiag_E0Um=0       ! DIAG_DOC: ${\cal E}^0\cdot\Uv$
  integer :: idiag_E0Wm=0       ! DIAG_DOC: ${\cal E}^0\cdot\Wv$
  integer :: idiag_bx0mz=0      ! DIAG_DOC: $\left<b_{x}\right>_{xy}$
  integer :: idiag_by0mz=0      ! DIAG_DOC: $\left<b_{y}\right>_{xy}$
  integer :: idiag_bz0mz=0      ! DIAG_DOC: $\left<b_{z}\right>_{xy}$
  integer :: idiag_M11z=0       ! DIAG_DOC: $\left<{\cal M}_{11}\right>_{xy}$
  integer :: idiag_M22z=0       ! DIAG_DOC: $\left<{\cal M}_{22}\right>_{xy}$
  integer :: idiag_M33z=0       ! DIAG_DOC: $\left<{\cal M}_{33}\right>_{xy}$
!
!  arrays for horizontally averaged uxb and jxb
!
  real, dimension (mz,3,mtestfield/6) :: uxbtestmz,jxbtestmz

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
      iaztestpq=iaatest+3*njtest-1
!
!  Allocate mtestfield slots; the first half is used for aatest
!  and the second for uutest.
!
      iuutest=nvar+1+mtestfield/2
      iuxtest=iuutest
      iuytest=iuutest+1
      iuztest=iuutest+2
      iuztestpq=iuutest+3*njtest-1
!
!  set ntestfield and nvar
!
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
      do j=1,ntestfield/2
        varname(iaatest-1+j) = 'aatest'
      enddo
!
      do j=ntestfield/2+1,ntestfield
        varname(iaatest-1+j) = 'uutest'
      enddo
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
!
!  Precalculate nutest if 1/nutest (==nutest1) is given instead
!
      if (nutest1/=0.) then
        nutest=1./nutest1
      endif
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
!  calculate iE0
!
      if (.not.lstarting) then
        select case (itestfield)
          case ('B11-B22'); iE0=0
        case default
          call fatal_error('initialize_testfield','undefined itestfield value')
        endselect
      endif
!
!  Override iE0 if njtest is big enough.
!  This method of identifying the location of the reference field is not very elegant.
!  Don't rescale the reference field, so put rescale_aatest=rescale_uutest=1.
!
      if (njtest==5) then
        iE0=5
        rescale_aatest(iE0)=1.
        rescale_uutest(iE0)=1.
      endif
!
!  set to zero and then rescale the testfield
!  (in future, could call something like init_aa_simple)
!
      if (reinitialize_aatest) then
        do jtest=1,njtest
          iaxtest=iaatest+3*(jtest-1); iaztest=iaxtest+2
          iuxtest=iuutest+3*(jtest-1); iuztest=iuxtest+2
          f(:,:,:,iaxtest:iaztest)=rescale_aatest(jtest)*f(:,:,:,iaxtest:iaztest)
          f(:,:,:,iuxtest:iuztest)=rescale_uutest(jtest)*f(:,:,:,iuxtest:iuztest)
        enddo
      endif
!
!  set lrescaling_testfield=T if linit_aatest=T
!
      if (linit_aatest) then
        lrescaling_testfield=.true.
      endif
!
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
        write(1,'(3a)') "itestfield_method='",trim(itestfield_method)//"'"
        write(1,'(a,f5.2)') 'ktestfield=',ktestfield
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
!  27-nov-09/axel: adapted from init_aatest in testfield_z
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
      lpenc_diagnos(i_bbb)=.true.
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
!  and du^(pq)/dt=Jxb^(pq)+jxB^(pq)+jxb-<jxb>+eta*del2U^(pq),
!    where p=1,2 and q=1 (if B11-B21) and optionally q=2 (if B11-B22),
!  and  da^(0)/dt=uxb-<uxb>+Eext+eta*del2A^(0),
!  with du^(0)/dt=jxb-<jxb>+Fext+eta*del2U^(0),
!
!  also calculate corresponding Lorentz force in connection with
!  testflow method
!
!   3-jun-05/axel: coded
!  16-mar-08/axel: Lorentz force added for testfield method
!  27-nov-09/axel: adapted from testfield_z, and added velocity equation
!  25-jan-09/axel: added Maxwell stress tensor calculation
!
      use Cdata
      use Diagnostics
      use Hydro, only: uumz,lcalc_uumean
      use Magnetic, only: aamz,bbmz,jjmz,lcalc_aamean,B_ext_inv
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx,3) :: uxB,bbtest,B0_imposed,uum,umxbtest
      real, dimension (nx,3) :: B0test=0,J0test=0
      real, dimension (nx,3) :: duxbtest,djxbtest,eetest
      real, dimension (nx,3) :: uxbtest,uxbtestK,uxbtestM,uxbtestMK
      real, dimension (nx,3) :: jxbtest,jxbtestK,jxbtestM,jxbtestMK
      real, dimension (nx,3,3,njtest) :: Mijpq
      real, dimension (nx,3,njtest) :: Eipq,upq
      real, dimension (nx,3,njtest) :: Fipq,bpq,jpq
      real, dimension (nx) :: alpK,alpM,alpMK
      real, dimension (nx) :: phiK,phiM,phiMK
      real, dimension (nx,3) :: del2Atest,uufluct,bbfluct,jjfluct
      real, dimension (nx,3) :: graddivAtest,aatest,jjtest
      real, dimension (nx,3) :: jxbrtest,jxbtest1,jxbtest2
      real, dimension (nx,3) :: del2Utest,uutest
      real, dimension (nx,3) :: u0ref,b0ref,j0ref
      real, dimension (nx,3,3) :: aijtest,bijtest,Mijtest
      real, dimension (nx) :: jbpq,upq2,jpq2,bpq2,Epq2,s2kzDF1,s2kzDF2,unity=1.
      integer :: jtest,j, i1=1, i2=2, i3=3, i4=4
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
!  loop over all fields, but do it backwards,
!  so we compute the zero field first
!
      do jtest=njtest,1,-1
        iaxtest=iaatest+3*(jtest-1); iaztest=iaxtest+2
        iuxtest=iuutest+3*(jtest-1); iuztest=iuxtest+2
!
!  compute uutest, bbtest, etc
!
        aatest=f(l1:l2,m,n,iaxtest:iaztest)
        uutest=f(l1:l2,m,n,iuxtest:iuztest)
!       call del2v(f,iaxtest,del2Atest)
        call del2v(f,iuxtest,del2Utest)
        call gij(f,iaxtest,aijtest,1)
        call gij_etc(f,iaxtest,aatest,aijtest,bijtest,del2Atest,graddivatest)
        call curl_mn(aijtest,bbtest,aatest)
        call curl_mn(bijtest,jjtest,bbtest)
!-test- call del2v_etc(f,iaxtest,CURLCURL=jjtest)
!
!  Get u0ref, b0ref, and j0ref (if iE0=5).
!  Also compute u0 x b0 and j0 x b0, and put into corresponding array.
!  They continue to exist throughout the jtest loop.
!
        if (jtest==iE0) then
          u0ref=uutest
          b0ref=bbtest
          j0ref=jjtest
!          call cross_mn(u0ref,b0ref,uxbtest)
!          call cross_mn(j0ref,b0ref,jxbtest)
        endif
!
!  do diffusion terms
!
        df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
          +etatest*del2Atest
        df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) &
          +nutest*del2Utest
!
!  With imposed field, calculate uutest x B0 and jjtest x B0 terms.
!  This applies to all terms, including the reference fields.
!
        do j=1,3
          B0_imposed(:,j)=Btest_ext(j)
        enddo
        call cross_mn(uutest,B0_imposed,uxbtest)
        call cross_mn(jjtest,B0_imposed,jxbtest)
        df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+uxbtest
        df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest)+jxbtest
!
!  add Ubar x b^0 and Ubar x b^T terms
!
        if (lcalc_uumean) then
          do j=1,3
            uum(:,j)=uumz(n,j)
          enddo
          call cross_mn(uum,bbtest,umxbtest)
          df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+umxbtest
        endif
!
!  possibility of non-soca terms
!
        if (.not.lsoca) then
!
!  subtract average emf, unless we ignore the <uxb> term (lignore_uxbtestm=T)
!
          if (iuxb/=0.and..not.ltest_uxb) then
            uxbtest=f(l1:l2,m,n,iuxb+3*(jtest-1):iuxb+3*jtest-1)
            if (lignore_uxbtestm) then
              duxbtest(:,:)=uxbtest(:,:)
            else
              do j=1,3
                duxbtest(:,j)=uxbtest(:,j)-uxbtestmz(n,j,jtest)
              enddo
            endif
            df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
                +duxbtest
          endif
!
!  subtract average jxb, unless we ignore the <jxb> term (lignore_jxbtestm=T)
!
          if (ijxb/=0.and..not.ltest_jxb) then
            jxbtest=f(l1:l2,m,n,ijxb+3*(jtest-1):ijxb+3*jtest-1)
            if (lignore_jxbtestm) then
              djxbtest(:,:)=jxbtest(:,:)
            else
              do j=1,3
                djxbtest(:,j)=jxbtest(:,j)-jxbtestmz(n,j,jtest)
              enddo
            endif
            df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) &
              +djxbtest
          endif
!
!  end of .not.lsoca
!
        endif
!
!  Do things differently for the reference fields.
!
        if (jtest==iE0) then
!
!  Add possibility of forcing that is not delta-correlated in time.
!
          if (lforcing_cont_aatest) &
            df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
              +ampl_fcont_aatest*p%fcont
          if (lforcing_cont_uutest) &
            df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) &
              +ampl_fcont_uutest*p%fcont
        else
!
!  Calculate uufluct=U-Umean.
!-  Note that uumz has dimensions mz*3, not nz*3.
!
          if (lcalc_uumean) then
            do j=1,3
              uufluct(:,j)=p%uu(:,j)-uumz(n,j)
            enddo
          else
            uufluct=p%uu
          endif
!
!  Calculate bbfluct=B-Bmean and jjfluct=J-Jmean.
!-  Note that, unlike uumz, bbmz and jjmz have dimensions nz*3.
!
          if (lcalc_aamean) then
            do j=1,3
              bbfluct(:,j)=p%bbb(:,j)-bbmz(n-n1+1,j)
              jjfluct(:,j)=p%jj(:,j)-jjmz(n-n1+1,j)
            enddo
          else
            bbfluct=p%bbb
            jjfluct=p%jj
          endif
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further down in the file.
!
          select case (itestfield)
          case ('B11-B22')
            call set_bbtest_B11_B22(B0test,jtest)
            call set_J0test_B11_B22(J0test,jtest)
          case default
            call fatal_error('daatest_dt','undefined itestfield value')
          endselect
!
!  u x B^T
!
          call cross_mn(uufluct,B0test,uxB)
!
!  jxB^T + J^Txb
!
          call cross_mn(jjfluct,B0test,jxbtest1)
          call cross_mn(J0test,bbfluct,jxbtest2)
          !call multsv_mn(p%rho1,jxbtest1+jxbtest2,jxbrtest)
          jxbrtest=jxbtest1+jxbtest2
          df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+uxB
          df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest)+jxbrtest
!
!  advance test field equation
!
        endif
!
!  calculate alpha, begin by calculating uxbtest (if not already done above)
!
        if ((ldiagnos.or.l1davgfirst).and. &
          (lsoca.or.ltest_uxb.or.idiag_b0rms/=0.or. &
           idiag_j11rms/=0.or.idiag_b11rms/=0.or.idiag_b21rms/=0.or. &
           idiag_b12rms/=0.or.idiag_b22rms/=0.or. &
           idiag_s2kzDFm/=0.or. &
           idiag_M11cc/=0.or.idiag_M11ss/=0.or. &
           idiag_M22cc/=0.or.idiag_M22ss/=0.or. &
           idiag_M12cs/=0.or. &
           idiag_M11/=0.or.idiag_M22/=0.or.idiag_M33/=0.or. &
           idiag_M11z/=0.or.idiag_M22z/=0.or.idiag_M33z/=0)) then
!
!  uutest and bbtest were only needed for the non-SOCA parts that
!  are calculated elsewhere. Therefore, for diagnostics, we need
!  to recalculate them here.
!
          uutest=f(l1:l2,m,n,iuxtest:iuztest)
          aatest=f(l1:l2,m,n,iaxtest:iaztest)
          call gij(f,iaxtest,aijtest,1)
          call curl_mn(aijtest,bbtest,aatest)
!
!  Settings for quadratic quantities (magnetic stress)
!
          if (idiag_M11cc/=0.or.idiag_M11ss/=0.or. &
              idiag_M22cc/=0.or.idiag_M22ss/=0.or. &
              idiag_M12cs/=0.or. &
              idiag_M11/=0.or.idiag_M22/=0.or.idiag_M33/=0.or. &
              idiag_M11z/=0.or.idiag_M22z/=0.or.idiag_M33z/=0) then
            call dyadic2(bbtest,Mijtest)
            Mijpq(:,:,:,jtest)=Mijtest*bamp12
          endif
        endif
        bpq(:,:,jtest)=bbtest
        upq(:,:,jtest)=uutest
!
!  Restore uxbtest and jxbtest from f-array, and compute uxbtestK for alpK
!  computation for comparison. Do the same for jxb.
!
        uxbtest=f(l1:l2,m,n,iuxb+3*(jtest-1):iuxb+3*jtest-1)
        jxbtest=f(l1:l2,m,n,ijxb+3*(jtest-1):ijxb+3*jtest-1)
!
!  evaluate different contributions to <uxb> and <jxb>
!
        Eipq(:,:,jtest)=uxbtest*bamp1
        Fipq(:,:,jtest)=jxbtest*bamp1
        if (ldiagnos.and.(idiag_jb0m/=0.or.idiag_j11rms/=0)) &
            jpq(:,:,jtest)=jjtest
!
!  enddo loop for jtest
!
      enddo
!
!  compute kinetic, magnetic, and magneto-kinetic contributions
!
      if (any(B_ext_inv/=0.)) then
        call cross_mn(p%uu,p%bbb-b0ref,uxbtestK)
        call cross_mn(p%uu-u0ref,p%bbb,uxbtestM)
        call cross_mn(p%uu-u0ref,p%bbb-b0ref,uxbtestMK)
        call cross_mn(p%jj,p%bbb-b0ref,jxbtestK)
        call cross_mn(p%jj-j0ref,p%bbb,jxbtestM)
        call cross_mn(p%jj-j0ref,p%bbb-b0ref,jxbtestMK)
        call dot(B_ext_inv,uxbtestK,alpK)
        call dot(B_ext_inv,uxbtestM,alpM)
        call dot(B_ext_inv,uxbtestMK,alpMK)
        call dot(B_ext_inv,jxbtestK,phiK)
        call dot(B_ext_inv,jxbtestM,phiM)
        call dot(B_ext_inv,jxbtestMK,phiMK)
      endif
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
        call xysum_mn_name_z(bpq(:,1,iE0),idiag_bx0mz)
        call xysum_mn_name_z(bpq(:,2,iE0),idiag_by0mz)
        call xysum_mn_name_z(bpq(:,3,iE0),idiag_bz0mz)
        call xysum_mn_name_z(Eipq(:,1,i1),idiag_E111z)
        call xysum_mn_name_z(Eipq(:,2,i1),idiag_E211z)
        call xysum_mn_name_z(Eipq(:,3,i1),idiag_E311z)
        call xysum_mn_name_z(Eipq(:,1,i2),idiag_E121z)
        call xysum_mn_name_z(Eipq(:,2,i2),idiag_E221z)
        call xysum_mn_name_z(Eipq(:,3,i2),idiag_E321z)
        call xysum_mn_name_z(Eipq(:,1,i3),idiag_E112z)
        call xysum_mn_name_z(Eipq(:,2,i3),idiag_E212z)
        call xysum_mn_name_z(Eipq(:,3,i3),idiag_E312z)
        call xysum_mn_name_z(Eipq(:,1,i4),idiag_E122z)
        call xysum_mn_name_z(Eipq(:,2,i4),idiag_E222z)
        call xysum_mn_name_z(Eipq(:,3,i4),idiag_E322z)
        call xysum_mn_name_z(Eipq(:,1,iE0),idiag_E10z)
        call xysum_mn_name_z(Eipq(:,2,iE0),idiag_E20z)
        call xysum_mn_name_z(Eipq(:,3,iE0),idiag_E30z)
        call xysum_mn_name_z(Mijpq(:,1,1,i1),idiag_M11z)
        call xysum_mn_name_z(Mijpq(:,2,2,i1),idiag_M22z)
        call xysum_mn_name_z(Mijpq(:,3,3,i1),idiag_M33z)
!
!  averages of alpha and eta
!
        if (idiag_alp11/=0) call sum_mn_name(+cz(n)*Eipq(:,1,1)+sz(n)*Eipq(:,1,i2),idiag_alp11)
        if (idiag_alp21/=0) call sum_mn_name(+cz(n)*Eipq(:,2,1)+sz(n)*Eipq(:,2,i2),idiag_alp21)
        if (idiag_alp31/=0) call sum_mn_name(+cz(n)*Eipq(:,3,1)+sz(n)*Eipq(:,3,i2),idiag_alp31)
        if (idiag_eta12/=0) call sum_mn_name(-(-sz(n)*Eipq(:,1,i1)+cz(n)*Eipq(:,1,i2))*ktestfield1,idiag_eta12)
        if (idiag_eta22/=0) call sum_mn_name(-(-sz(n)*Eipq(:,2,i1)+cz(n)*Eipq(:,2,i2))*ktestfield1,idiag_eta22)
!
!  same for jxb
!
        if (idiag_phi11/=0) call sum_mn_name(+cz(n)*Fipq(:,1,1)+sz(n)*Fipq(:,1,i2),idiag_phi11)
        if (idiag_phi21/=0) call sum_mn_name(+cz(n)*Fipq(:,2,1)+sz(n)*Fipq(:,2,i2),idiag_phi21)
        if (idiag_psi12/=0) call sum_mn_name(-(-sz(n)*Fipq(:,1,i1)+cz(n)*Fipq(:,1,i2))*ktestfield1,idiag_psi12)
        if (idiag_psi22/=0) call sum_mn_name(-(-sz(n)*Fipq(:,2,i1)+cz(n)*Fipq(:,2,i2))*ktestfield1,idiag_psi22)
!
!  weighted averages alpha and eta
!  Still need to do this for iE0 /= 0 case.
!
        if (idiag_alp11cc/=0) call sum_mn_name(c2z(n)*(+cz(n)*Eipq(:,1,1)+sz(n)*Eipq(:,1,i2)),idiag_alp11cc)
        if (idiag_alp21sc/=0) call sum_mn_name(csz(n)*(+cz(n)*Eipq(:,2,1)+sz(n)*Eipq(:,2,i2)),idiag_alp21sc)
        if (idiag_eta12cs/=0) call sum_mn_name(-csz(n)*(-sz(n)*Eipq(:,1,i1)+cz(n)*Eipq(:,1,i2))*ktestfield1,idiag_eta12cs)
        if (idiag_eta22ss/=0) call sum_mn_name(-s2z(n)*(-sz(n)*Eipq(:,2,i1)+cz(n)*Eipq(:,2,i2))*ktestfield1,idiag_eta22ss)
!
!  Compute kinetic, magnetic, and magneto-kinetic contributions with
!  imposed-field method
!
        if (idiag_alpK/=0) call sum_mn_name(alpK,idiag_alpK)
        if (idiag_alpM/=0) call sum_mn_name(alpM,idiag_alpM)
        if (idiag_alpMK/=0) call sum_mn_name(alpMK,idiag_alpMK)
        if (idiag_phiK/=0) call sum_mn_name(phiK,idiag_phiK)
        if (idiag_phiM/=0) call sum_mn_name(phiM,idiag_phiM)
        if (idiag_phiMK/=0) call sum_mn_name(phiMK,idiag_phiMK)
!
!  Divergence of current helicity flux
!
        if (idiag_s2kzDFm/=0) then
          eetest=etatest*jjtest-duxbtest
          s2kzDF1=2*ktestfield*c2kz(n)*(&
            bijtest(:,1,3)*eetest(:,1)+&
            bijtest(:,2,3)*eetest(:,2)+&
            bijtest(:,3,3)*eetest(:,3))
          s2kzDF2=4*ktestfield**2*s2kz(n)*(&
            bbtest(:,1)*eetest(:,1)+&
            bbtest(:,2)*eetest(:,2))
          call sum_mn_name(s2kzDF1-s2kzDF2,idiag_s2kzDFm)
        endif
!
!  Maxwell tensor and its weighted averages
!
        if (idiag_M11/=0)   call sum_mn_name(       Mijpq(:,1,1,i1),idiag_M11)
        if (idiag_M22/=0)   call sum_mn_name(       Mijpq(:,2,2,i1),idiag_M22)
        if (idiag_M33/=0)   call sum_mn_name(       Mijpq(:,3,3,i1),idiag_M33)
        if (idiag_M11cc/=0) call sum_mn_name(c2z(n)*Mijpq(:,1,1,i1),idiag_M11cc)
        if (idiag_M11ss/=0) call sum_mn_name(s2z(n)*Mijpq(:,1,1,i1),idiag_M11ss)
        if (idiag_M22cc/=0) call sum_mn_name(c2z(n)*Mijpq(:,2,2,i1),idiag_M22cc)
        if (idiag_M22ss/=0) call sum_mn_name(s2z(n)*Mijpq(:,2,2,i1),idiag_M22ss)
        if (idiag_M12cs/=0) call sum_mn_name(cz(n)*sz(n)*Mijpq(:,1,2,i1),idiag_M12cs)
!
!  Projection of EMF from testfield against testfield itself
!
        if (idiag_EBpq/=0) call sum_mn_name(cz(n)*Eipq(:,1,1) &
                                           +sz(n)*Eipq(:,2,1),idiag_EBpq)
!
!  print warning if alp12 and alp12 are needed, but njtest is too small XX
!
        if (njtest<=2 .and. &
            (idiag_alp12/=0.or.idiag_alp22/=0.or.idiag_alp32/=0.or. &
            (leta_rank2.and.(idiag_eta11/=0.or.idiag_eta21/=0)).or. &
            (.not.leta_rank2.and.(idiag_eta12/=0.or.idiag_eta22/=0)).or. &
            (leta_rank2.and.(idiag_eta11cc/=0.or.idiag_eta21sc/=0)))) then
          call stop_it('njtest is too small if alpi2 or etai2 for i=1,2,3 are needed')
        else
!
!  Remaining coefficients
!
          if (idiag_alp12/=0) call sum_mn_name(+cz(n)*Eipq(:,1,i3)+sz(n)*Eipq(:,1,i4),idiag_alp12)
          if (idiag_alp22/=0) call sum_mn_name(+cz(n)*Eipq(:,2,i3)+sz(n)*Eipq(:,2,i4),idiag_alp22)
          if (idiag_alp32/=0) call sum_mn_name(+cz(n)*Eipq(:,3,i3)+sz(n)*Eipq(:,3,i4),idiag_alp32)
          if (idiag_alp12cs/=0) call sum_mn_name(csz(n)*(+cz(n)*Eipq(:,1,i3)+sz(n)*Eipq(:,1,i4)),idiag_alp12cs)
          if (idiag_alp22ss/=0) call sum_mn_name(s2z(n)*(+cz(n)*Eipq(:,2,i3)+sz(n)*Eipq(:,2,i4)),idiag_alp22ss)
          if (idiag_eta11/=0) call sum_mn_name((-sz(n)*Eipq(:,1,i3)+cz(n)*Eipq(:,1,i4))*ktestfield1,idiag_eta11)
          if (idiag_eta21/=0) call sum_mn_name((-sz(n)*Eipq(:,2,i3)+cz(n)*Eipq(:,2,i4))*ktestfield1,idiag_eta21)
          if (idiag_eta11cc/=0) call sum_mn_name(c2z(n)*(-sz(n)*Eipq(:,1,i3)+cz(n)*Eipq(:,1,i4))*ktestfield1,idiag_eta11cc)
          if (idiag_eta21sc/=0) call sum_mn_name(csz(n)*(-sz(n)*Eipq(:,2,i3)+cz(n)*Eipq(:,2,i4))*ktestfield1,idiag_eta21sc)
!
!  same for jxb
!
          if (idiag_phi12/=0) call sum_mn_name(+cz(n)*Fipq(:,1,i3)+sz(n)*Fipq(:,1,i4),idiag_phi12)
          if (idiag_phi22/=0) call sum_mn_name(+cz(n)*Fipq(:,2,i3)+sz(n)*Fipq(:,2,i4),idiag_phi22)
          if (idiag_phi32/=0) call sum_mn_name(+cz(n)*Fipq(:,3,i3)+sz(n)*Fipq(:,3,i4),idiag_phi32)
          if (idiag_psi11/=0) call sum_mn_name((-sz(n)*Fipq(:,1,i3)+cz(n)*Fipq(:,1,i4))*ktestfield1,idiag_psi11)
          if (idiag_psi21/=0) call sum_mn_name((-sz(n)*Fipq(:,2,i3)+cz(n)*Fipq(:,2,i4))*ktestfield1,idiag_psi21)
        endif
!
!  Volume-averaged dot products of mean emf and velocity and of mean emf and vorticity
!
        if (iE0/=0) then
          if (idiag_E0Um/=0) call sum_mn_name(uxbtestmz(n,1,iE0)*p%uu(:,1) &
                                             +uxbtestmz(n,2,iE0)*p%uu(:,2) &
                                             +uxbtestmz(n,3,iE0)*p%uu(:,3),idiag_E0Um)
          if (idiag_E0Wm/=0) call sum_mn_name(uxbtestmz(n,1,iE0)*p%oo(:,1) &
                                             +uxbtestmz(n,2,iE0)*p%oo(:,2) &
                                             +uxbtestmz(n,3,iE0)*p%oo(:,3),idiag_E0Wm)
        endif
!
!  diagnostics for delta function driving, but doesn't seem to work
!
        if (idiag_bamp/=0) call sum_mn_name(bamp*unity,idiag_bamp)
!
!  diagnostics for single points
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_bx0pt/=0)  call save_name(bpq(lpoint-nghost,1,iE0),idiag_bx0pt)
          if (idiag_bx11pt/=0) call save_name(bpq(lpoint-nghost,1,i1),idiag_bx11pt)
          if (idiag_bx21pt/=0) call save_name(bpq(lpoint-nghost,1,i2),idiag_bx21pt)
          if (idiag_bx12pt/=0) call save_name(bpq(lpoint-nghost,1,i3),idiag_bx12pt)
          if (idiag_bx22pt/=0) call save_name(bpq(lpoint-nghost,1,i4),idiag_bx22pt)
          if (idiag_by0pt/=0)  call save_name(bpq(lpoint-nghost,2,iE0),idiag_by0pt)
          if (idiag_by11pt/=0) call save_name(bpq(lpoint-nghost,2,i1),idiag_by11pt)
          if (idiag_by21pt/=0) call save_name(bpq(lpoint-nghost,2,i2),idiag_by21pt)
          if (idiag_by12pt/=0) call save_name(bpq(lpoint-nghost,2,i3),idiag_by12pt)
          if (idiag_by22pt/=0) call save_name(bpq(lpoint-nghost,2,i4),idiag_by22pt)
          if (idiag_Ex0pt/=0)  call save_name(Eipq(lpoint-nghost,1,iE0),idiag_Ex0pt)
          if (idiag_Ex11pt/=0) call save_name(Eipq(lpoint-nghost,1,i1),idiag_Ex11pt)
          if (idiag_Ex21pt/=0) call save_name(Eipq(lpoint-nghost,1,i2),idiag_Ex21pt)
          if (idiag_Ex12pt/=0) call save_name(Eipq(lpoint-nghost,1,i3),idiag_Ex12pt)
          if (idiag_Ex22pt/=0) call save_name(Eipq(lpoint-nghost,1,i4),idiag_Ex22pt)
          if (idiag_Ey0pt/=0)  call save_name(Eipq(lpoint-nghost,2,iE0),idiag_Ey0pt)
!         if (idiag_bamp/=0)   call save_name(bamp,idiag_bamp)
          if (idiag_Ey11pt/=0) call save_name(Eipq(lpoint-nghost,2,i1),idiag_Ey11pt)
          if (idiag_Ey21pt/=0) call save_name(Eipq(lpoint-nghost,2,i2),idiag_Ey21pt)
          if (idiag_Ey12pt/=0) call save_name(Eipq(lpoint-nghost,2,i3),idiag_Ey12pt)
          if (idiag_Ey22pt/=0) call save_name(Eipq(lpoint-nghost,2,i4),idiag_Ey22pt)
        endif
!
!  rms values of small scales fields bpq in response to the test fields Bpq
!  Obviously idiag_b0rms and idiag_b12rms cannot both be invoked!
!  Needs modification!
!
        if (idiag_jb0m/=0) then
          call dot(jpq(:,:,iE0),bpq(:,:,iE0),jbpq)
          call sum_mn_name(jbpq,idiag_jb0m)
        endif
!
        if (idiag_ux0m/=0) call sum_mn_name(upq(:,1,iE0),idiag_ux0m)
        if (idiag_uy0m/=0) call sum_mn_name(upq(:,2,iE0),idiag_uy0m)
        if (idiag_ux11m/=0) call sum_mn_name(upq(:,1,i1),idiag_ux11m)
        if (idiag_uy11m/=0) call sum_mn_name(upq(:,2,i1),idiag_uy11m)
!
        if (idiag_u0rms/=0) then
          call dot2(upq(:,:,iE0),bpq2)
          call sum_mn_name(bpq2,idiag_u0rms,lsqrt=.true.)
        endif
!
        if (idiag_b0rms/=0) then
          call dot2(bpq(:,:,iE0),bpq2)
          call sum_mn_name(bpq2,idiag_b0rms,lsqrt=.true.)
        endif
!
        if (idiag_u11rms/=0) then
          call dot2(upq(:,:,1),upq2)
          call sum_mn_name(upq2,idiag_u11rms,lsqrt=.true.)
        endif
!
        if (idiag_u21rms/=0) then
          call dot2(upq(:,:,i2),upq2)
          call sum_mn_name(upq2,idiag_u21rms,lsqrt=.true.)
        endif
!
        if (idiag_u12rms/=0) then
          call dot2(upq(:,:,i3),upq2)
          call sum_mn_name(upq2,idiag_u12rms,lsqrt=.true.)
        endif
!
        if (idiag_u22rms/=0) then
          call dot2(upq(:,:,i4),upq2)
          call sum_mn_name(upq2,idiag_u22rms,lsqrt=.true.)
        endif
!
        if (idiag_b11rms/=0) then
          call dot2(bpq(:,:,1),bpq2)
          call sum_mn_name(bpq2,idiag_b11rms,lsqrt=.true.)
        endif
!
        if (idiag_j11rms/=0) then
          call dot2(jpq(:,:,1),jpq2)
          call sum_mn_name(jpq2,idiag_j11rms,lsqrt=.true.)
          !-test- call sum_mn_name(jpq(:,1,1)**2,idiag_j11rms,lsqrt=.true.)
        endif
!
        if (idiag_b21rms/=0) then
          call dot2(bpq(:,:,i2),bpq2)
          call sum_mn_name(bpq2,idiag_b21rms,lsqrt=.true.)
        endif
!
        if (idiag_b12rms/=0) then
          call dot2(bpq(:,:,i3),bpq2)
          call sum_mn_name(bpq2,idiag_b12rms,lsqrt=.true.)
        endif
!
        if (idiag_b22rms/=0) then
          call dot2(bpq(:,:,i4),bpq2)
          call sum_mn_name(bpq2,idiag_b22rms,lsqrt=.true.)
        endif
!
        if (idiag_E0rms/=0) then
          call dot2(Eipq(:,:,iE0),Epq2)
          call sum_mn_name(Epq2,idiag_E0rms,lsqrt=.true.)
        endif
!
        if (idiag_E11rms/=0) then
          call dot2(Eipq(:,:,i1),Epq2)
          call sum_mn_name(Epq2,idiag_E11rms,lsqrt=.true.)
        endif
!
        if (idiag_E21rms/=0) then
          call dot2(Eipq(:,:,i2),Epq2)
          call sum_mn_name(Epq2,idiag_E21rms,lsqrt=.true.)
        endif
!
        if (idiag_E12rms/=0) then
          call dot2(Eipq(:,:,i3),Epq2)
          call sum_mn_name(Epq2,idiag_E12rms,lsqrt=.true.)
        endif
!
        if (idiag_E22rms/=0) then
          call dot2(Eipq(:,:,i4),Epq2)
          call sum_mn_name(Epq2,idiag_E22rms,lsqrt=.true.)
        endif
!
      endif
!
!  write B-slices for output in wvid in run.f90
!  Note: ix is the index with respect to array with ghost zones.
!
      if (lvideo.and.lfirst) then
        do j=1,3
          bb11_yz(m-m1+1,n-n1+1,j)=bpq(ix_loc-l1+1,j,1)
          if (m==iy_loc)  bb11_xz(:,n-n1+1,j)=bpq(:,j,1)
          if (n==iz_loc)  bb11_xy(:,m-m1+1,j)=bpq(:,j,1)
          if (n==iz2_loc) bb11_xy2(:,m-m1+1,j)=bpq(:,j,1)
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
        case ('bb11')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =>bb11_yz(:,:,slices%index)
            slices%xz =>bb11_xz(:,:,slices%index)
            slices%xy =>bb11_xy(:,:,slices%index)
            slices%xy2=>bb11_xy2(:,:,slices%index)
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
!  Calculate <uxb^T> + <u^Txb>, which is needed when lsoca=.false.
!  Also calculate <jxb^T> + <j^Txb>, which is needed when lsoca_jxb=.false.
!
!  30-nov-09/axel: adapted from testfield_z.f90
!
      use Cdata
      use Sub
      use Hydro, only: calc_pencils_hydro,uumz,lcalc_uumean
      use Magnetic, only: calc_pencils_magnetic, idiag_bcosphz, idiag_bsinphz, &
        aamz,bbmz,jjmz,lcalc_aamean
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
      real, dimension (nx,3) :: aatest,bbtest,jjtest,uutest,uxbtest,jxbtest
      real, dimension (nx,3) :: uxbtest1,jxbtest1
      real, dimension (nx,3) :: uxbtest2,jxbtest2
      real, dimension (nx,3) :: del2Atest2,graddivatest
      real, dimension (nx,3) :: u0ref,b0ref,j0ref
      real, dimension (nx,3) :: uufluct,bbfluct,jjfluct
      integer :: jtest,j,nxy=nxgrid*nygrid,juxb,jjxb
      logical :: headtt_save
      real :: fac, bcosphz, bsinphz
!
      intent(inout) :: f
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac=1./nxy
!
!  Stop if iE0 is too small.
!
      if (iE0/=5) then
        call fatal_error('testfield_after_boundary','need ntest=5 for u0ref')
      endif
!
!  initialize counter for mean fields
!
      uxbtestmz=0.
      jxbtestmz=0.
!
!  Start mn loop
!
      do n=n1,n2
      do m=m1,m2
!
!  Begin by getting/computing fields from main run.
!
        call calc_pencils_hydro(f,p)
        call calc_pencils_magnetic(f,p)
!
!  Calculate uufluct=U-Umean.
!-  Note that uumz has dimensions mz*3, not nz*3.
!
        if (lcalc_uumean) then
          do j=1,3
            uufluct(:,j)=p%uu(:,j)-uumz(n,j)
          enddo
        else
          uufluct=p%uu
        endif
!
!  Calculate bbfluct=B-Bmean and jjfluct=J-Jmean.
!-  Note that, unlike uumz, bbmz and jjmz have dimensions nz*3.
!
        if (lcalc_aamean) then
          do j=1,3
            bbfluct(:,j)=p%bbb(:,j)-bbmz(n-n1+1,j)
            jjfluct(:,j)=p%jj(:,j)-jjmz(n-n1+1,j)
          enddo
        else
          bbfluct=p%bbb
          jjfluct=p%jj
        endif
!
!  Count jtest backward, so we have already access to the reference fields.
!
        do jtest=njtest,1,-1
!
!  Compute test fields aatest, bbtest, jjtest, and uutest,
!  and set bbref, jjref, uuref if jtest=iE0 (first point in loop)
!
          iaxtest=iaatest+3*(jtest-1); iaztest=iaxtest+2
          iuxtest=iuutest+3*(jtest-1); iuztest=iuxtest+2
          aatest=f(l1:l2,m,n,iaxtest:iaztest)
          uutest=f(l1:l2,m,n,iuxtest:iuztest)
!
!  calculate uxb, and put it into f-array
!
          call gij(f,iaxtest,aijtest,1)
          call gij_etc(f,iaxtest,aatest,aijtest,bijtest,del2Atest2,graddivatest)
          call curl_mn(aijtest,bbtest,aatest)
          call curl_mn(bijtest,jjtest,bbtest)
!
!  Get u0ref, b0ref, and j0ref (if iE0=5).
!  Also compute u0 x b0 and j0 x b0, and put into corresponding array.
!  They continue to exist throughout the jtest loop.
!
          if (jtest==iE0) then
            u0ref=uutest
            b0ref=bbtest
            j0ref=jjtest
            call cross_mn(u0ref,b0ref,uxbtest)
            call cross_mn(j0ref,b0ref,jxbtest)
          else
!
!  Calculate uxb and jxb, and put it into f-array, depending on whether we use
!  Testfield Method (i) or (ii), or mixed ones of either (iii) or (iv).
!
            select case (itestfield_method)
            case ('ju', '(i)')
              call cross_mn(jjfluct,bbtest,jxbtest1)
              call cross_mn(jjtest,b0ref,jxbtest2)
              call cross_mn(uufluct,bbtest,uxbtest1)
              call cross_mn(uutest,b0ref,uxbtest2)
            case ('bb', '(ii)')
              call cross_mn(j0ref,bbtest,jxbtest1)
              call cross_mn(jjtest,bbfluct,jxbtest2)
              call cross_mn(u0ref,bbtest,uxbtest1)
              call cross_mn(uutest,bbfluct,uxbtest2)
            case ('bu', '(iii)')
              call cross_mn(j0ref,bbtest,jxbtest1)
              call cross_mn(jjtest,bbfluct,jxbtest2)
              call cross_mn(uufluct,bbtest,uxbtest1)
              call cross_mn(uutest,b0ref,uxbtest2)
            case ('jb', '(iv)')
              call cross_mn(jjfluct,bbtest,jxbtest1)
              call cross_mn(jjtest,b0ref,jxbtest2)
              call cross_mn(u0ref,bbtest,uxbtest1)
              call cross_mn(uutest,bbfluct,uxbtest2)
            case default
              call fatal_error('testfield_after_boundary','??itestfield_method')
            endselect
            uxbtest=uxbtest1+uxbtest2
            jxbtest=jxbtest1+jxbtest2
          endif
          juxb=iuxb+3*(jtest-1)
          jjxb=ijxb+3*(jtest-1)
          if (iuxb/=0) f(l1:l2,m,n,juxb:juxb+2)=uxbtest
          if (ijxb/=0) f(l1:l2,m,n,jjxb:jjxb+2)=jxbtest
!
!  Add corresponding contribution into averaged arrays, uxbtestmz, jxbtestmz
!
          do j=1,3
            uxbtestmz(n,j,jtest)=uxbtestmz(n,j,jtest)+fac*sum(uxbtest(:,j))
            jxbtestmz(n,j,jtest)=jxbtestmz(n,j,jtest)+fac*sum(jxbtest(:,j))
          enddo
          headtt=.false.
!
!  finish jtest and mn loops
!
        enddo
      enddo
      enddo
!
!  Set uxbtestm1 and jxbtestm1 for communication.
!
      do jtest=njtest,1,-1
        do j=1,3
          do n=n1,n2
            uxbtestm1(n-n1+1,ipz+1,j,jtest)=uxbtestmz(n,j,jtest)
            jxbtestm1(n-n1+1,ipz+1,j,jtest)=jxbtestmz(n,j,jtest)
          enddo
        enddo
      enddo
!
!  do communication for array of size nz*nprocz*3*njtest
!
      if (nprocy>1) then
        call mpireduce_sum(uxbtestm1,uxbtestm1_tmp,(/nz,nprocz,3,njtest/))
        call mpireduce_sum(jxbtestm1,jxbtestm1_tmp,(/nz,nprocz,3,njtest/))
        call mpibcast_real_arr(uxbtestm1_tmp,nz*nprocz*3*njtest)
        call mpibcast_real_arr(jxbtestm1_tmp,nz*nprocz*3*njtest)
        do jtest=1,njtest
          do n=n1,n2
            do j=1,3
              uxbtestmz(n,j,jtest)=uxbtestm1_tmp(n-n1+1,ipz+1,j,jtest)
              jxbtestmz(n,j,jtest)=jxbtestm1_tmp(n-n1+1,ipz+1,j,jtest)
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
      if (ip<9) print*,'iproc,phase_testfield=',iproc,phase_testfield
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
      use Hydro, only: uumz,lcalc_uumean
      use Magnetic, only: aamz,bbmz,jjmz,lcalc_aamean
!
      real, dimension (mx,my,mz,mfarray) :: f
      character (len=fnlen) :: file
      logical :: ltestfield_out
      integer,save :: ifirst=0
      integer :: j,jtest,jaatest,juutest,jaa,juu
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
          do jtest=1,njtest
            iaxtest=iaatest+3*(jtest-1); iaztest=iaxtest+2
            iuxtest=iuutest+3*(jtest-1); iuztest=iuxtest+2
            do n=n1,n2
              do j=iaxtest,iaztest
                f(l1:l2,m1:m2,n,j)=rescale_aatest(jtest)*f(l1:l2,m1:m2,n,j)
              enddo
              do j=iuxtest,iuztest
                f(l1:l2,m1:m2,n,j)=rescale_uutest(jtest)*f(l1:l2,m1:m2,n,j)
              enddo
            enddo
          enddo
!
!  Reinitialize reference fields with fluctuations of main run.
!
          if (reinitialize_from_mainrun) then
            if (lcalc_aamean.and.lcalc_uumean) then
              jtest=iE0
              iaxtest=iaatest+3*(jtest-1)
              iuxtest=iuutest+3*(jtest-1)
              do n=n1,n2
                do j=1,3
                  jaatest=iaxtest+j-1;  jaa=iax+j-1
                  juutest=iuxtest+j-1;  juu=iux+j-1
                  f(l1:l2,m1:m2,n,jaatest)=f(l1:l2,m1:m2,n,jaa)-aamz(n,j)
                  f(l1:l2,m1:m2,n,juutest)=f(l1:l2,m1:m2,n,juu)-uumz(n,j)
                enddo
              enddo
            else
              call fatal_error('rescaling_testfield', &
                  'need lcalc_aamean.and.lcalc_uumean')
            endif
          endif
!
!  Update next time for rescaling
!
          call update_snaptime(file,taainit,naainit,daainit,t,ltestfield_out)
        endif
      endif
!
    endsubroutine rescaling_testfield
!***********************************************************************
    subroutine set_bbtest_Beltrami(B0test,jtest)
!
!  set testfield
!
!  29-mar-08/axel: coded
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
      case (1); B0test(:,1)=bamp*cz(n); B0test(:,2)=bamp*sz(n); B0test(:,3)=0.
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_Beltrami
!***********************************************************************
    subroutine set_bbtest_B11_B21(B0test,jtest)
!
!  set testfield
!
!   3-jun-05/axel: coded
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
      case (1); B0test(:,1)=bamp*cz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=bamp*sz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_B11_B21
!***********************************************************************
    subroutine set_J0test_B11_B21(J0test,jtest)
!
!  set testfield
!
!   3-jun-05/axel: coded
!
      use Cdata
!
      real, dimension (nx,3) :: J0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: J0test
!
!  set J0test for each of the 9 cases
!
      select case (jtest)
      case (1); J0test(:,1)=0.; J0test(:,2)=-bamp*ktestfield*sz(n); J0test(:,3)=0.
      case (2); J0test(:,1)=0.; J0test(:,2)=+bamp*ktestfield*cz(n); J0test(:,3)=0.
      case default; J0test(:,:)=0.
      endselect
!
    endsubroutine set_J0test_B11_B21
!***********************************************************************
    subroutine set_J0test_B11_B22(J0test,jtest)
!
!  set testfield
!
!   3-jun-05/axel: coded
!
      use Cdata
!
      real, dimension (nx,3) :: J0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: J0test
!
!  set J0test for each of the 9 cases
!
      select case (jtest)
      case (1); J0test(:,1)=0.; J0test(:,2)=-bamp*ktestfield*sz(n); J0test(:,3)=0.
      case (2); J0test(:,1)=0.; J0test(:,2)=+bamp*ktestfield*cz(n); J0test(:,3)=0.
      case (3); J0test(:,1)=+bamp*ktestfield*sz(n); J0test(:,2)=0.; J0test(:,3)=0.
      case (4); J0test(:,1)=-bamp*ktestfield*cz(n); J0test(:,2)=0.; J0test(:,3)=0.
      case default; J0test(:,:)=0.
      endselect
!
    endsubroutine set_J0test_B11_B22
!***********************************************************************
    subroutine set_bbtest_B11_B22 (B0test,jtest)
!
!  set testfield
!
!   3-jun-05/axel: coded
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
      case (1); B0test(:,1)=bamp*cz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=bamp*sz(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=bamp*cz(n); B0test(:,3)=0.
      case (4); B0test(:,1)=0.; B0test(:,2)=bamp*sz(n); B0test(:,3)=0.
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_B11_B22
!***********************************************************************
    subroutine set_bbtest_B11_B22_lin (B0test,jtest)
!
!  set testfield
!
!  25-Mar-09/axel: adapted from set_bbtest_B11_B22
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
      case (1); B0test(:,1)=bamp     ; B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=bamp*z(n); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=bamp     ; B0test(:,3)=0.
      case (4); B0test(:,1)=0.; B0test(:,2)=bamp*z(n); B0test(:,3)=0.
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_B11_B22_lin
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
        idiag_bx0mz=0; idiag_by0mz=0; idiag_bz0mz=0
        idiag_E111z=0; idiag_E211z=0; idiag_E311z=0
        idiag_E121z=0; idiag_E221z=0; idiag_E321z=0
        idiag_E10z=0; idiag_E20z=0; idiag_E30z=0
        idiag_EBpq=0; idiag_E0Um=0; idiag_E0Wm=0
        idiag_alp11=0; idiag_alp21=0; idiag_alp31=0
        idiag_alp12=0; idiag_alp22=0; idiag_alp32=0
        idiag_eta11=0; idiag_eta21=0
        idiag_eta12=0; idiag_eta22=0
        idiag_phi11=0; idiag_phi21=0
        idiag_phi12=0; idiag_phi22=0; idiag_phi32=0
        idiag_psi11=0; idiag_psi21=0
        idiag_psi12=0; idiag_psi22=0
        idiag_alp11cc=0; idiag_alp21sc=0; idiag_alp12cs=0; idiag_alp22ss=0
        idiag_eta11cc=0; idiag_eta21sc=0; idiag_eta12cs=0; idiag_eta22ss=0
        idiag_alpK=0; idiag_alpM=0; idiag_alpMK=0
        idiag_phiK=0; idiag_phiM=0; idiag_phiMK=0
        idiag_s2kzDFm=0
        idiag_M11=0; idiag_M22=0; idiag_M33=0
        idiag_M11cc=0; idiag_M11ss=0; idiag_M22cc=0; idiag_M22ss=0
        idiag_M12cs=0
        idiag_M11z=0; idiag_M22z=0; idiag_M33z=0; idiag_bamp=0
        idiag_jb0m=0; idiag_u0rms=0; idiag_b0rms=0; idiag_E0rms=0
        idiag_ux0m=0; idiag_uy0m=0
        idiag_ux11m=0; idiag_uy11m=0
        idiag_u11rms=0; idiag_u21rms=0; idiag_u12rms=0; idiag_u22rms=0
        idiag_j11rms=0; idiag_b11rms=0; idiag_b21rms=0; idiag_b12rms=0; idiag_b22rms=0
        idiag_E11rms=0; idiag_E21rms=0; idiag_E12rms=0; idiag_E22rms=0
        idiag_bx0pt=0; idiag_bx11pt=0; idiag_bx21pt=0; idiag_bx12pt=0; idiag_bx22pt=0
        idiag_by0pt=0; idiag_by11pt=0; idiag_by21pt=0; idiag_by12pt=0; idiag_by22pt=0
        idiag_Ex0pt=0; idiag_Ex11pt=0; idiag_Ex21pt=0; idiag_Ex12pt=0; idiag_Ex22pt=0
        idiag_Ey0pt=0; idiag_Ey11pt=0; idiag_Ey21pt=0; idiag_Ey12pt=0; idiag_Ey22pt=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'alp11',idiag_alp11)
        call parse_name(iname,cname(iname),cform(iname),'alp21',idiag_alp21)
        call parse_name(iname,cname(iname),cform(iname),'alp31',idiag_alp31)
        call parse_name(iname,cname(iname),cform(iname),'alp12',idiag_alp12)
        call parse_name(iname,cname(iname),cform(iname),'alp22',idiag_alp22)
        call parse_name(iname,cname(iname),cform(iname),'alp32',idiag_alp32)
        call parse_name(iname,cname(iname),cform(iname),'eta11',idiag_eta11)
        call parse_name(iname,cname(iname),cform(iname),'eta21',idiag_eta21)
        call parse_name(iname,cname(iname),cform(iname),'eta12',idiag_eta12)
        call parse_name(iname,cname(iname),cform(iname),'eta22',idiag_eta22)
        call parse_name(iname,cname(iname),cform(iname),'alpK',idiag_alpK)
        call parse_name(iname,cname(iname),cform(iname),'alpM',idiag_alpM)
        call parse_name(iname,cname(iname),cform(iname),'alpMK',idiag_alpMK)
        call parse_name(iname,cname(iname),cform(iname),'phi11',idiag_phi11)
        call parse_name(iname,cname(iname),cform(iname),'phi21',idiag_phi21)
        call parse_name(iname,cname(iname),cform(iname),'phi12',idiag_phi12)
        call parse_name(iname,cname(iname),cform(iname),'phi22',idiag_phi22)
        call parse_name(iname,cname(iname),cform(iname),'phi32',idiag_phi32)
        call parse_name(iname,cname(iname),cform(iname),'psi11',idiag_psi11)
        call parse_name(iname,cname(iname),cform(iname),'psi21',idiag_psi21)
        call parse_name(iname,cname(iname),cform(iname),'psi12',idiag_psi12)
        call parse_name(iname,cname(iname),cform(iname),'psi22',idiag_psi22)
        call parse_name(iname,cname(iname),cform(iname),'phiK',idiag_phiK)
        call parse_name(iname,cname(iname),cform(iname),'phiM',idiag_phiM)
        call parse_name(iname,cname(iname),cform(iname),'phiMK',idiag_phiMK)
        call parse_name(iname,cname(iname),cform(iname),'alp11cc',idiag_alp11cc)
        call parse_name(iname,cname(iname),cform(iname),'alp21sc',idiag_alp21sc)
        call parse_name(iname,cname(iname),cform(iname),'alp12cs',idiag_alp12cs)
        call parse_name(iname,cname(iname),cform(iname),'alp22ss',idiag_alp22ss)
        call parse_name(iname,cname(iname),cform(iname),'eta11cc',idiag_eta11cc)
        call parse_name(iname,cname(iname),cform(iname),'eta21sc',idiag_eta21sc)
        call parse_name(iname,cname(iname),cform(iname),'eta12cs',idiag_eta12cs)
        call parse_name(iname,cname(iname),cform(iname),'eta22ss',idiag_eta22ss)
        call parse_name(iname,cname(iname),cform(iname),'s2kzDFm',idiag_s2kzDFm)
        call parse_name(iname,cname(iname),cform(iname),'M11',idiag_M11)
        call parse_name(iname,cname(iname),cform(iname),'M22',idiag_M22)
        call parse_name(iname,cname(iname),cform(iname),'M33',idiag_M33)
        call parse_name(iname,cname(iname),cform(iname),'M11cc',idiag_M11cc)
        call parse_name(iname,cname(iname),cform(iname),'M11ss',idiag_M11ss)
        call parse_name(iname,cname(iname),cform(iname),'M22cc',idiag_M22cc)
        call parse_name(iname,cname(iname),cform(iname),'M22ss',idiag_M22ss)
        call parse_name(iname,cname(iname),cform(iname),'M12cs',idiag_M12cs)
        call parse_name(iname,cname(iname),cform(iname),'bx11pt',idiag_bx11pt)
        call parse_name(iname,cname(iname),cform(iname),'bx21pt',idiag_bx21pt)
        call parse_name(iname,cname(iname),cform(iname),'bx12pt',idiag_bx12pt)
        call parse_name(iname,cname(iname),cform(iname),'bx22pt',idiag_bx22pt)
        call parse_name(iname,cname(iname),cform(iname),'bx0pt',idiag_bx0pt)
        call parse_name(iname,cname(iname),cform(iname),'by11pt',idiag_by11pt)
        call parse_name(iname,cname(iname),cform(iname),'by21pt',idiag_by21pt)
        call parse_name(iname,cname(iname),cform(iname),'by12pt',idiag_by12pt)
        call parse_name(iname,cname(iname),cform(iname),'by22pt',idiag_by22pt)
        call parse_name(iname,cname(iname),cform(iname),'by0pt',idiag_by0pt)
        call parse_name(iname,cname(iname),cform(iname),'Ex11pt',idiag_Ex11pt)
        call parse_name(iname,cname(iname),cform(iname),'Ex21pt',idiag_Ex21pt)
        call parse_name(iname,cname(iname),cform(iname),'Ex12pt',idiag_Ex12pt)
        call parse_name(iname,cname(iname),cform(iname),'Ex22pt',idiag_Ex22pt)
        call parse_name(iname,cname(iname),cform(iname),'Ex0pt',idiag_Ex0pt)
        call parse_name(iname,cname(iname),cform(iname),'Ey11pt',idiag_Ey11pt)
        call parse_name(iname,cname(iname),cform(iname),'Ey21pt',idiag_Ey21pt)
        call parse_name(iname,cname(iname),cform(iname),'Ey12pt',idiag_Ey12pt)
        call parse_name(iname,cname(iname),cform(iname),'Ey22pt',idiag_Ey22pt)
        call parse_name(iname,cname(iname),cform(iname),'Ey0pt',idiag_Ey0pt)
        call parse_name(iname,cname(iname),cform(iname),'u11rms',idiag_u11rms)
        call parse_name(iname,cname(iname),cform(iname),'u21rms',idiag_u21rms)
        call parse_name(iname,cname(iname),cform(iname),'u12rms',idiag_u12rms)
        call parse_name(iname,cname(iname),cform(iname),'u22rms',idiag_u22rms)
        call parse_name(iname,cname(iname),cform(iname),'j11rms',idiag_j11rms)
        call parse_name(iname,cname(iname),cform(iname),'b11rms',idiag_b11rms)
        call parse_name(iname,cname(iname),cform(iname),'b21rms',idiag_b21rms)
        call parse_name(iname,cname(iname),cform(iname),'b12rms',idiag_b12rms)
        call parse_name(iname,cname(iname),cform(iname),'b22rms',idiag_b22rms)
        call parse_name(iname,cname(iname),cform(iname),'jb0m',idiag_jb0m)
        call parse_name(iname,cname(iname),cform(iname),'bamp',idiag_bamp)
        call parse_name(iname,cname(iname),cform(iname),'ux0m',idiag_ux0m)
        call parse_name(iname,cname(iname),cform(iname),'uy0m',idiag_uy0m)
        call parse_name(iname,cname(iname),cform(iname),'ux11m',idiag_ux11m)
        call parse_name(iname,cname(iname),cform(iname),'uy11m',idiag_uy11m)
        call parse_name(iname,cname(iname),cform(iname),'u0rms',idiag_u0rms)
        call parse_name(iname,cname(iname),cform(iname),'b0rms',idiag_b0rms)
        call parse_name(iname,cname(iname),cform(iname),'E11rms',idiag_E11rms)
        call parse_name(iname,cname(iname),cform(iname),'E21rms',idiag_E21rms)
        call parse_name(iname,cname(iname),cform(iname),'E12rms',idiag_E12rms)
        call parse_name(iname,cname(iname),cform(iname),'E22rms',idiag_E22rms)
        call parse_name(iname,cname(iname),cform(iname),'E0rms',idiag_E0rms)
        call parse_name(iname,cname(iname),cform(iname),'EBpq',idiag_EBpq)
        call parse_name(iname,cname(iname),cform(iname),'E0Um',idiag_E0Um)
        call parse_name(iname,cname(iname),cform(iname),'E0Wm',idiag_E0Wm)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bx0mz',idiag_bx0mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'by0mz',idiag_by0mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'bz0mz',idiag_bz0mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E111z',idiag_E111z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E211z',idiag_E211z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E311z',idiag_E311z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E121z',idiag_E121z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E221z',idiag_E221z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E321z',idiag_E321z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E112z',idiag_E112z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E212z',idiag_E212z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E312z',idiag_E312z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E122z',idiag_E122z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E222z',idiag_E222z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E322z',idiag_E322z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E10z',idiag_E10z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E20z',idiag_E20z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'E30z',idiag_E30z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'M11z',idiag_M11z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'M22z',idiag_M22z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'M33z',idiag_M33z)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'idiag_alp11=',idiag_alp11
        write(3,*) 'idiag_alp21=',idiag_alp21
        write(3,*) 'idiag_alp31=',idiag_alp31
        write(3,*) 'idiag_alp12=',idiag_alp12
        write(3,*) 'idiag_alp22=',idiag_alp22
        write(3,*) 'idiag_alp32=',idiag_alp32
        write(3,*) 'idiag_eta11=',idiag_eta11
        write(3,*) 'idiag_eta21=',idiag_eta21
        write(3,*) 'idiag_eta12=',idiag_eta12
        write(3,*) 'idiag_eta22=',idiag_eta22
        write(3,*) 'idiag_alpK=',idiag_alpK
        write(3,*) 'idiag_alpM=',idiag_alpM
        write(3,*) 'idiag_alpMK=',idiag_alpMK
        write(3,*) 'idiag_phi11=',idiag_phi11
        write(3,*) 'idiag_phi21=',idiag_phi21
        write(3,*) 'idiag_phi12=',idiag_phi12
        write(3,*) 'idiag_phi22=',idiag_phi22
        write(3,*) 'idiag_phi32=',idiag_phi32
        write(3,*) 'idiag_psi11=',idiag_psi11
        write(3,*) 'idiag_psi21=',idiag_psi21
        write(3,*) 'idiag_psi12=',idiag_psi12
        write(3,*) 'idiag_psi22=',idiag_psi22
        write(3,*) 'idiag_phiK=',idiag_phiK
        write(3,*) 'idiag_phiM=',idiag_phiM
        write(3,*) 'idiag_phiMK=',idiag_phiMK
        write(3,*) 'idiag_alp11cc=',idiag_alp11cc
        write(3,*) 'idiag_alp21sc=',idiag_alp21sc
        write(3,*) 'idiag_alp12cs=',idiag_alp12cs
        write(3,*) 'idiag_alp22ss=',idiag_alp22ss
        write(3,*) 'idiag_eta11cc=',idiag_eta11cc
        write(3,*) 'idiag_eta21sc=',idiag_eta21sc
        write(3,*) 'idiag_eta12cs=',idiag_eta12cs
        write(3,*) 'idiag_eta22ss=',idiag_eta22ss
        write(3,*) 'idiag_s2kzDFm=',idiag_s2kzDFm
        write(3,*) 'idiag_M11=',idiag_M11
        write(3,*) 'idiag_M22=',idiag_M22
        write(3,*) 'idiag_M33=',idiag_M33
        write(3,*) 'idiag_M11cc=',idiag_M11cc
        write(3,*) 'idiag_M11ss=',idiag_M11ss
        write(3,*) 'idiag_M22cc=',idiag_M22cc
        write(3,*) 'idiag_M22ss=',idiag_M22ss
        write(3,*) 'idiag_M12cs=',idiag_M12cs
        write(3,*) 'idiag_bx11pt=',idiag_bx11pt
        write(3,*) 'idiag_bx21pt=',idiag_bx21pt
        write(3,*) 'idiag_bx12pt=',idiag_bx12pt
        write(3,*) 'idiag_bx22pt=',idiag_bx22pt
        write(3,*) 'idiag_bx0pt=',idiag_bx0pt
        write(3,*) 'idiag_by11pt=',idiag_by11pt
        write(3,*) 'idiag_by21pt=',idiag_by21pt
        write(3,*) 'idiag_by12pt=',idiag_by12pt
        write(3,*) 'idiag_by22pt=',idiag_by22pt
        write(3,*) 'idiag_by0pt=',idiag_by0pt
        write(3,*) 'idiag_Ex11pt=',idiag_Ex11pt
        write(3,*) 'idiag_Ex21pt=',idiag_Ex21pt
        write(3,*) 'idiag_Ex12pt=',idiag_Ex12pt
        write(3,*) 'idiag_Ex22pt=',idiag_Ex22pt
        write(3,*) 'idiag_Ex0pt=',idiag_Ex0pt
        write(3,*) 'idiag_Ey11pt=',idiag_Ey11pt
        write(3,*) 'idiag_Ey21pt=',idiag_Ey21pt
        write(3,*) 'idiag_Ey12pt=',idiag_Ey12pt
        write(3,*) 'idiag_Ey22pt=',idiag_Ey22pt
        write(3,*) 'idiag_Ey0pt=',idiag_Ey0pt
        write(3,*) 'idiag_u11rms=',idiag_u11rms
        write(3,*) 'idiag_u21rms=',idiag_u21rms
        write(3,*) 'idiag_u12rms=',idiag_u12rms
        write(3,*) 'idiag_u22rms=',idiag_u22rms
        write(3,*) 'idiag_j11rms=',idiag_j11rms
        write(3,*) 'idiag_b11rms=',idiag_b11rms
        write(3,*) 'idiag_b21rms=',idiag_b21rms
        write(3,*) 'idiag_b12rms=',idiag_b12rms
        write(3,*) 'idiag_b22rms=',idiag_b22rms
        write(3,*) 'idiag_jb0m=',idiag_jb0m
        write(3,*) 'idiag_ux0m=',idiag_ux0m
        write(3,*) 'idiag_uy0m=',idiag_uy0m
        write(3,*) 'idiag_ux11m=',idiag_ux11m
        write(3,*) 'idiag_uy11m=',idiag_uy11m
        write(3,*) 'idiag_u0rms=',idiag_u0rms
        write(3,*) 'idiag_b0rms=',idiag_b0rms
        write(3,*) 'idiag_E11rms=',idiag_E11rms
        write(3,*) 'idiag_E21rms=',idiag_E21rms
        write(3,*) 'idiag_E12rms=',idiag_E12rms
        write(3,*) 'idiag_E22rms=',idiag_E22rms
        write(3,*) 'idiag_E0rms=',idiag_E0rms
        write(3,*) 'idiag_bx0mz=',idiag_bx0mz
        write(3,*) 'idiag_by0mz=',idiag_by0mz
        write(3,*) 'idiag_bz0mz=',idiag_bz0mz
        write(3,*) 'idiag_bamp=',idiag_bamp
        write(3,*) 'idiag_E111z=',idiag_E111z
        write(3,*) 'idiag_E211z=',idiag_E211z
        write(3,*) 'idiag_E311z=',idiag_E311z
        write(3,*) 'idiag_E121z=',idiag_E121z
        write(3,*) 'idiag_E221z=',idiag_E221z
        write(3,*) 'idiag_E321z=',idiag_E321z
        write(3,*) 'idiag_E112z=',idiag_E112z
        write(3,*) 'idiag_E212z=',idiag_E212z
        write(3,*) 'idiag_E312z=',idiag_E312z
        write(3,*) 'idiag_E122z=',idiag_E122z
        write(3,*) 'idiag_E222z=',idiag_E222z
        write(3,*) 'idiag_E322z=',idiag_E322z
        write(3,*) 'idiag_E10z=',idiag_E10z
        write(3,*) 'idiag_E20z=',idiag_E20z
        write(3,*) 'idiag_E30z=',idiag_E30z
        write(3,*) 'idiag_M11z=',idiag_M11z
        write(3,*) 'idiag_M22z=',idiag_M22z
        write(3,*) 'idiag_M33z=',idiag_M33z
        write(3,*) 'idiag_EBpq=',idiag_EBpq
        write(3,*) 'idiag_E0Um=',idiag_E0Um
        write(3,*) 'idiag_E0Wm=',idiag_E0Wm
        write(3,*) 'iaatest=',iaatest
        write(3,*) 'iuutest=',iuutest
        write(3,*) 'ntestfield=',ntestfield
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'nnamexz=',nnamexz
      endif
!
    endsubroutine rprint_testfield

endmodule Testfield
