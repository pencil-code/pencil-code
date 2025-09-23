! $Id$

!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.

!  Note: this routine requires that MVAR and MAUX contributions
!  together with njtest are set correctly in the cparam.local file.
!  njtest must be set at the end of the file such that 3*njtest=MVAR.
!
!  Example:
!  ! MVAR CONTRIBUTION 12
!  ! MAUX CONTRIBUTION 12
!  integer, parameter :: njtest=4

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltestfield = .true.
! CPARAM logical, parameter :: ltestfield_x = .true.
! CPARAM logical, parameter :: ltestfield_z = .false.
! CPARAM logical, parameter :: ltestfield_xy = .false.
! CPARAM logical, parameter :: ltestfield_xz  = .false.
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
  real, dimension(nx) :: cx,sx,xx,x2
!
  character (len=labellen), dimension(ninit) :: initaatest='nothing'
  real :: kx_aatest=1.,ky_aatest=1.,kz_aatest=1.
  real, dimension (ninit) :: amplaatest=0.
  integer :: iE0=0

  ! input parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./)
  real :: taainit=0.,daainit=0.
  logical :: reinitialize_aatest=.false.
  logical :: zextent=.true.,lsoca=.false.,lsoca_jxb=.true.,lset_bbtest2=.false.
  logical :: luxb_as_aux=.false.,ljxb_as_aux=.false.,linit_aatest=.false.
  logical :: lignore_uxbtestm=.false.
  character (len=labellen) :: itestfield='B11-B21'
  real :: ktestfield=1., ktestfield1=1.
  integer :: naainit
  real :: bamp=1.
  namelist /testfield_init_pars/ &
       B_ext,zextent,initaatest, &
       amplaatest,kx_aatest,ky_aatest,kz_aatest, &
       luxb_as_aux,ljxb_as_aux

  ! run parameters
  real :: etatest=0.,etatest1=0.
  real, dimension(njtest) :: rescale_aatest=0.
  logical :: ltestfield_newx=.true.,leta_rank2=.true.
  namelist /testfield_run_pars/ &
       B_ext,reinitialize_aatest,zextent,lsoca,lsoca_jxb, &
       lset_bbtest2,etatest,etatest1,itestfield,ktestfield, &
       ltestfield_newx, &
       luxb_as_aux,ljxb_as_aux,lignore_uxbtestm, &
       daainit,linit_aatest,bamp, &
       rescale_aatest

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
!
!  sin and cos weighted averages
!
  integer :: idiag_alp11cc=0    ! DIAG_DOC: $\alpha_{11}\cos^2 kx$
  integer :: idiag_alp21sc=0    ! DIAG_DOC: $\alpha_{21}\sin kx\cos kx$
  integer :: idiag_alp12cs=0    ! DIAG_DOC: $\alpha_{12}\cos kx\sin kx$
  integer :: idiag_alp22ss=0    ! DIAG_DOC: $\alpha_{22}\sin^2 kx$
  integer :: idiag_eta11cc=0    ! DIAG_DOC: $\eta_{11}\cos^2 kx$
  integer :: idiag_eta21sc=0    ! DIAG_DOC: $\eta_{21}\sin kx\cos kx$
  integer :: idiag_eta12cs=0    ! DIAG_DOC: $\eta_{12}\cos kx\sin kx$
  integer :: idiag_eta22ss=0    ! DIAG_DOC: $\eta_{22}\sin^2 kx$
!
!  x-weighted averages
!
  integer :: idiag_alp11_x=0      ! DIAG_DOC: $\alpha_{11}x$
  integer :: idiag_alp21_x=0      ! DIAG_DOC: $\alpha_{21}x$
  integer :: idiag_alp12_x=0      ! DIAG_DOC: $\alpha_{12}x$
  integer :: idiag_alp22_x=0      ! DIAG_DOC: $\alpha_{22}x$
  integer :: idiag_eta11_x=0      ! DIAG_DOC: $\eta_{11}kx$
  integer :: idiag_eta21_x=0      ! DIAG_DOC: $\eta_{21}kx$
  integer :: idiag_eta12_x=0      ! DIAG_DOC: $\eta_{12}kx$
  integer :: idiag_eta22_x=0      ! DIAG_DOC: $\eta_{22}kx$
!
!  x^2-weighted averages
!
  integer :: idiag_alp11_x2=0      ! DIAG_DOC: $\alpha_{11}x^2$
  integer :: idiag_alp21_x2=0      ! DIAG_DOC: $\alpha_{21}x^2$
  integer :: idiag_alp12_x2=0      ! DIAG_DOC: $\alpha_{12}x^2$
  integer :: idiag_alp22_x2=0      ! DIAG_DOC: $\alpha_{22}x^2$
  integer :: idiag_eta11_x2=0      ! DIAG_DOC: $\eta_{11}kx^2$
  integer :: idiag_eta21_x2=0      ! DIAG_DOC: $\eta_{21}kx^2$
  integer :: idiag_eta12_x2=0      ! DIAG_DOC: $\eta_{12}kx^2$
  integer :: idiag_eta22_x2=0      ! DIAG_DOC: $\eta_{22}kx^2$
!
!  other quantities
!
  integer :: idiag_b11rms=0     ! DIAG_DOC: $\left<b_{11}^2\right>^{1/2}$
  integer :: idiag_b21rms=0     ! DIAG_DOC: $\left<b_{21}^2\right>^{1/2}$
  integer :: idiag_b12rms=0     ! DIAG_DOC: $\left<b_{12}^2\right>^{1/2}$
  integer :: idiag_b22rms=0     ! DIAG_DOC: $\left<b_{22}^2\right>^{1/2}$
  integer :: idiag_b0rms=0      ! DIAG_DOC: $\left<b_{0}^2\right>^{1/2}$
  integer :: idiag_E11rms=0     ! DIAG_DOC: $\left<{\cal E}_{11}^2\right>^{1/2}$
  integer :: idiag_E21rms=0     ! DIAG_DOC: $\left<{\cal E}_{21}^2\right>^{1/2}$
  integer :: idiag_E12rms=0     ! DIAG_DOC: $\left<{\cal E}_{12}^2\right>^{1/2}$
  integer :: idiag_E22rms=0     ! DIAG_DOC: $\left<{\cal E}_{22}^2\right>^{1/2}$
  integer :: idiag_E0rms=0      ! DIAG_DOC: $\left<{\cal E}_{0}^2\right>^{1/2}$
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
  integer :: idiag_bx0mz=0      ! DIAG_DOC: $\left<b_{x}\right>_{xy}$
  integer :: idiag_by0mz=0      ! DIAG_DOC: $\left<b_{y}\right>_{xy}$
  integer :: idiag_bz0mz=0      ! DIAG_DOC: $\left<b_{z}\right>_{xy}$
  integer :: idiag_alp11x=0     ! DIAG_DOC: $\alpha_{11}(x,t)$
  integer :: idiag_alp21x=0     ! DIAG_DOC: $\alpha_{21}(x,t)$
  integer :: idiag_alp12x=0     ! DIAG_DOC: $\alpha_{12}(x,t)$
  integer :: idiag_alp22x=0     ! DIAG_DOC: $\alpha_{22}(x,t)$
  integer :: idiag_eta11x=0     ! DIAG_DOC: $\eta_{11}(x,t)$
  integer :: idiag_eta21x=0     ! DIAG_DOC: $\eta_{21}(x,t)$
  integer :: idiag_eta12x=0     ! DIAG_DOC: $\eta_{12}(x,t)$
  integer :: idiag_eta22x=0     ! DIAG_DOC: $\eta_{22}(x,t)$
!
!  arrays for horizontally averaged uxb and jxb
!
  real, dimension (nx,3,njtest) :: uxbtestm,jxbtestm

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
!  Here always ltestfield=T
!
      iaxtest=iaatest
      iaytest=iaatest+1
      iaztest=iaatest+2
      iaztestpq=iaatest+ntestfield-1
!
!  identify version number
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
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension(nx) :: xtestfield
      integer :: jtest
!
!  Precalculate etatest if 1/etatest (==etatest1) is given instead
!
      if (etatest1/=0.) then
        etatest=1./etatest1
      endif
!
!  set cosine and sine function for setting test fields and analysis
!  Choice of using rescaled x-array or original x-array
!
      if (ltestfield_newx) then
        xtestfield=2.*pi*(x(l1:l2)-x0)/Lx-pi
      else
        xtestfield=x(l1:l2)
      endif
      cx=cos(ktestfield*xtestfield)
      sx=sin(ktestfield*xtestfield)
      x2=xtestfield**2
      xx=xtestfield
!
!  Also calculate its inverse, but only if different from zero
!
      if (ktestfield==0) then
        ktestfield1=1.
      else
        ktestfield1=1./ktestfield
      endif
!
!  calculate iE0
!
      if (lrun) then
        select case (itestfield)
        case ('Beltrami'); iE0=1
        case ('B11-B21+B=0'); iE0=3
        case ('B11-B22+B=0'); iE0=5
        case ('B11-B21'); iE0=0
        case ('B11-B22'); iE0=0
        case ('B=0') !(dont do anything)
        case default
          call fatal_error('initialize_testfield','undefined itestfield value')
        endselect
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
!  write testfield information to a file (for convenient post-processing)
!
      if (lroot) then
        open(1,file=trim(datadir)//'/testfield_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'zextent=',merge(1,0,zextent)
        write(1,'(a,i1)') 'lsoca='  ,merge(1,0,lsoca)
        write(1,'(a,i1)') 'lsoca_jxb='  ,merge(1,0,lsoca_jxb)
        write(1,'(3a)') "itestfield='",trim(itestfield)//"'"
        write(1,'(a,f5.2)') 'ktestfield=',ktestfield
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
      case ('sinwave-x-1'); call sinwave(amplaatest(j),f,iaxtest+0+1,kx=kx_aatest)
      case ('sinwave-x-2'); call sinwave(amplaatest(j),f,iaxtest+3+1,kx=kx_aatest)
      case ('sinwave-x-3'); call sinwave(amplaatest(j),f,iaxtest+6+1,kx=kx_aatest)
      case ('Beltrami-z-1'); call beltrami(amplaatest(j),f,iaxtest+0,kz=kz_aatest)
      case ('Beltrami-z-3'); call beltrami(amplaatest(j),f,iaxtest+6,kz=kz_aatest)
      case ('Beltrami-z-5'); call beltrami(amplaatest(j),f,iaxtest+12,kz=kz_aatest)
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
!
    endsubroutine pencil_criteria_testfield
!***********************************************************************
    subroutine pencil_interdep_testfield(lpencil_in)
!
!  Interdependency among pencils from the Testfield module is specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
      logical, dimension(npencils) :: lpencil_in
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
!  19-sep-13/MR  : corrected uumz -> uumx
!
      use Cdata
      use Diagnostics
      use Hydro, only: uumx,lcalc_uumeanx
      use Mpicomm, only: stop_it
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx,3) :: uxB,B0test=0,bbtest
      real, dimension (nx,3) :: uxbtest,duxbtest,jxbtest,djxbrtest
      real, dimension (nx,3) :: J0test,jxB0rtest,J0xbrtest
      real, dimension (nx,3,njtest) :: Eipq,bpq
      real, dimension (nx,3) :: del2Atest,uufluct
      real, dimension (nx,3) :: del2Atest2,graddivatest,aatest,jjtest,jxbrtest
      real, dimension (nx,3,3) :: aijtest,bijtest
      real, dimension (nx) :: bpq2,Epq2,diffus_eta
      integer :: jtest, j, iuxtest, iuytest, iuztest
      integer :: i1=1, i2=2, i3=3, i4=4
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
      if (lcalc_uumeanx) then
        do j=1,3
          !uufluct(:,j)=p%uu(:,j)-uumx(:,j)
!AB: quick fix
          uufluct(:,j)=p%uu(:,j)-uumx(l1:l2,j)
        enddo
      else
        uufluct=p%uu
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
          case ('Beltrami'); call set_bbtest_Beltrami(B0test,jtest)
          case ('B11-B21+B=0'); call set_bbtest_B11_B21(B0test,jtest)
          case ('B11-B22+B=0'); call set_bbtest_B11_B22(B0test,jtest)
          case ('B11-B21'); call set_bbtest_B11_B21(B0test,jtest)
          case ('B11-B22'); call set_bbtest_B11_B22(B0test,jtest)
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
!  add diffusion
!
        df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest) &
            +etatest*del2Atest
!
        call cross_mn(uufluct,B0test,uxB)
        df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+uxB
!
        if (.not.lsoca) then
!
!  Add non-SOCA terms:
!  use f-array for uxb (if space has been allocated for this) and
!  if we don't test (i.e. if ltest_uxb=.false.)
!
          if (iuxbtest/=0.and..not.ltest_uxb) then
            uxbtest=f(l1:l2,m,n,iuxbtest+3*(jtest-1):iuxbtest+3*jtest-1)
          else
            call curl(f,iaxtest,bbtest)
            call cross_mn(p%uu,bbtest,uxbtest)
          endif
!
!  subtract average emf, unless we ignore the <uxb> term (lignore_uxbtestm=T)
!
          if (lignore_uxbtestm) then
            duxbtest=uxbtest
          else
            duxbtest=uxbtest-uxbtestm(:,:,jtest)
          endif
!
!  advance test field equation
!
          df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+duxbtest
        endif
!
!  Calculate Lorentz force
!
        if (ltestflow) then
          iuxtest=iuutest+4*(jtest-1)                   !!!MR: only correct if jtest refers to number of testflow
          iuytest=iuxtest+1 !(even though its not used)
          iuztest=iuxtest+2
          aatest=f(l1:l2,m,n,iaxtest:iaztest)
          call gij(f,iaxtest,aijtest,1)
          call gij_etc(f,iaxtest,aatest,aijtest,bijtest,del2Atest2,graddivatest)
!
!  calculate jpq x bpq
!
          call curl_mn(aijtest,bbtest,aatest)
          call curl_mn(bijtest,jjtest,bbtest)
!
!  calculate jpq x B0pq
!
          call cross_mn(jjtest,B0test,jxB0rtest)
!
!  calculate J0pq x bpq
!
          select case (itestfield)
!           case ('B11-B21+B=0'); call set_J0test(J0test,jtest)
            case ('B11-B21'); call set_J0test_B11_B21(J0test,jtest)
!           case ('B11-B22'); call set_J0test_B11_B22(J0test,jtest)
            case ('B=0') !(dont do anything)
          case default
            call fatal_error('daatest_dt','undefined itestfield value')
          endselect
          call cross_mn(J0test,bbtest,J0xbrtest)
!
!  add them all together
!
        if (lsoca_jxb) then
          df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) &
            +jxB0rtest+J0xbrtest
        else
!
!  use f-array for uxb (if space has been allocated for this) and
!  if we don't test (i.e. if ltest_jxb=.false.)
!
          if (ijxbtest/=0.and..not.ltest_jxb) then
            jxbtest=f(l1:l2,m,n,ijxbtest+3*(jtest-1):ijxbtest+3*jtest-1)
          else
            call cross_mn(jjtest,bbtest,jxbrtest)
          endif
!
!  subtract average jxb
!
          djxbrtest=jxbtest-jxbtestm(:,:,jtest)
          df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) &
            +jxB0rtest+J0xbrtest+djxbrtest
        endif
        endif
!
!  calculate alpha, begin by calculating uxbtest (if not already done above)
!
        if ((ldiagnos.or.l1davgfirst).and. &
          ((lsoca.or.iuxbtest/=0).and.(.not.ltest_uxb))) then
          call curl(f,iaxtest,bbtest)
          call cross_mn(p%uu,bbtest,uxbtest)
        endif
        bpq(:,:,jtest)=bbtest
        Eipq(:,:,jtest)=uxbtest/bamp
      enddo
!
!  diffusive time step, just take the max of diffus_eta (if existent)
!  and whatever is calculated here
!
      if (lfirst.and.ldt) then
        diffus_eta=etatest*dxyz_2
        maxdiffus=max(maxdiffus,diffus_eta)
      endif
!
!  in the following block, we have already swapped the 4-6 entries with 7-9
!  The g95 compiler doesn't like to see an index that is out of bounds,
!  so prevent this warning by writing i3=3 and i4=4
!  Apparently, 1,2 really means 2,3, but it seems to work ok.
!
      if (ldiagnos) then
          call yzsum_mn_name_x(cx*Eipq(:,2,i1)+sx*Eipq(:,2,i2),idiag_alp11x)
          call yzsum_mn_name_x(cx*Eipq(:,3,i1)+sx*Eipq(:,3,i2),idiag_alp21x)
          call yzsum_mn_name_x(cx*Eipq(:,2,i3)+sx*Eipq(:,2,i4),idiag_alp12x)
          call yzsum_mn_name_x(cx*Eipq(:,3,i3)+sx*Eipq(:,3,i4),idiag_alp22x)
          call yzsum_mn_name_x(-(sx*Eipq(:,2,i3)-cx*Eipq(:,2,i4))*ktestfield1,idiag_eta11x)
          call yzsum_mn_name_x(+(sx*Eipq(:,2,i1)-cx*Eipq(:,2,i2))*ktestfield1,idiag_eta12x)
          call yzsum_mn_name_x(-(sx*Eipq(:,3,i3)-cx*Eipq(:,3,i4))*ktestfield1,idiag_eta21x)
          call yzsum_mn_name_x(+(sx*Eipq(:,3,i1)-cx*Eipq(:,3,i2))*ktestfield1,idiag_eta22x)
        if (idiag_bx0mz/=0) call yzsum_mn_name_x(bpq(:,1,i3),idiag_bx0mz)
        if (idiag_by0mz/=0) call yzsum_mn_name_x(bpq(:,2,i3),idiag_by0mz)
        if (idiag_bz0mz/=0) call yzsum_mn_name_x(bpq(:,3,i3),idiag_bz0mz)
        if (idiag_E111z/=0) call yzsum_mn_name_x(Eipq(:,1,1),idiag_E111z)
        if (idiag_E211z/=0) call yzsum_mn_name_x(Eipq(:,2,1),idiag_E211z)
        if (idiag_E121z/=0) call yzsum_mn_name_x(Eipq(:,1,2),idiag_E121z)
        if (idiag_E221z/=0) call yzsum_mn_name_x(Eipq(:,2,2),idiag_E221z)
        if (idiag_E112z/=0) call yzsum_mn_name_x(Eipq(:,1,i3),idiag_E112z)
        if (idiag_E212z/=0) call yzsum_mn_name_x(Eipq(:,2,i3),idiag_E212z)
        if (idiag_E122z/=0) call yzsum_mn_name_x(Eipq(:,1,i4),idiag_E122z)
        if (idiag_E222z/=0) call yzsum_mn_name_x(Eipq(:,2,i4),idiag_E222z)
        if (idiag_E10z/=0) call yzsum_mn_name_x(Eipq(:,1,iE0),idiag_E10z)
        if (idiag_E20z/=0) call yzsum_mn_name_x(Eipq(:,2,iE0),idiag_E20z)
!
!  averages of alpha and eta
!
        if (idiag_alp11/=0) call sum_mn_name(+cx*Eipq(:,2,1)+sx*Eipq(:,2,2),idiag_alp11)
        if (idiag_alp21/=0) call sum_mn_name(+cx*Eipq(:,3,1)+sx*Eipq(:,3,2),idiag_alp21)
        if (idiag_eta11/=0) call sum_mn_name((-sx*Eipq(:,2,i3)+cx*Eipq(:,2,i4))*ktestfield1,idiag_eta11)
        if (idiag_eta21/=0) call sum_mn_name((-sx*Eipq(:,3,i3)+cx*Eipq(:,3,i4))*ktestfield1,idiag_eta21)
!
!  weighted averages alpha and eta
!
        if (idiag_alp11cc/=0) call sum_mn_name(cx**2*(+cx*Eipq(:,2,1)+sx*Eipq(:,2,2)),idiag_alp11cc)
        if (idiag_alp21sc/=0) call sum_mn_name(sx*cx*(+cx*Eipq(:,3,1)+sx*Eipq(:,3,2)),idiag_alp21sc)
        if (idiag_eta11cc/=0) call sum_mn_name(cx**2*(-sx*Eipq(:,2,i3)+cx*Eipq(:,2,i4))*ktestfield1,idiag_eta11cc)
        if (idiag_eta21sc/=0) call sum_mn_name(sx*cx*(-sx*Eipq(:,3,i3)+cx*Eipq(:,3,i4))*ktestfield1,idiag_eta21sc)
!
!  x-weighted averages alpha and eta
!
        if (idiag_alp11_x/=0) call sum_mn_name(xx*(+cx*Eipq(:,2,1)+sx*Eipq(:,2,2)),idiag_alp11_x)
        if (idiag_alp21_x/=0) call sum_mn_name(xx*(+cx*Eipq(:,3,1)+sx*Eipq(:,3,2)),idiag_alp21_x)
        if (idiag_eta11_x/=0) call sum_mn_name(xx*(-sx*Eipq(:,2,i3)+cx*Eipq(:,2,i4))*ktestfield1,idiag_eta11_x)
        if (idiag_eta21_x/=0) call sum_mn_name(xx*(-sx*Eipq(:,3,i3)+cx*Eipq(:,3,i4))*ktestfield1,idiag_eta21_x)
!
!  x^2-weighted averages alpha and eta
!
        if (idiag_alp11_x2/=0) call sum_mn_name(x2*(+cx*Eipq(:,2,1)+sx*Eipq(:,2,2)),idiag_alp11_x2)
        if (idiag_alp21_x2/=0) call sum_mn_name(x2*(+cx*Eipq(:,3,1)+sx*Eipq(:,3,2)),idiag_alp21_x2)
        if (idiag_eta11_x2/=0) call sum_mn_name(x2*(-sx*Eipq(:,2,i3)+cx*Eipq(:,2,i4))*ktestfield1,idiag_eta11_x2)
        if (idiag_eta21_x2/=0) call sum_mn_name(x2*(-sx*Eipq(:,3,i3)+cx*Eipq(:,3,i4))*ktestfield1,idiag_eta21_x2)
!
!  Projection of EMF from testfield against testfield itself
!
        if (idiag_EBpq/=0) call sum_mn_name(cx*Eipq(:,2,1) &
                                           +sx*Eipq(:,3,1),idiag_EBpq)
!
!  print warning if alp12 and alp12 are needed, but njtest is too small XX
!
        if ((idiag_alp12/=0.or.idiag_alp22/=0 &
         .or.idiag_eta12/=0.or.idiag_eta22/=0 &
         .or.idiag_alp12cs/=0.or.idiag_alp22ss/=0  &
         .or.idiag_eta12cs/=0.or.idiag_eta22ss/=0 &
         .or.idiag_alp12x/=0.or.idiag_alp22x/=0 &
         .or.idiag_eta12x/=0.or.idiag_eta22x/=0 &
        ).and.njtest<=2) then
          call stop_it('njtest is too small if alpi2 or etai2 for i=1,2,3 are needed')
        else
          if (idiag_alp12/=0) call sum_mn_name(+cx*Eipq(:,2,i3)+sx*Eipq(:,2,i4),idiag_alp12)
          if (idiag_alp22/=0) call sum_mn_name(+cx*Eipq(:,3,i3)+sx*Eipq(:,3,i4),idiag_alp22)
          if (idiag_alp12cs/=0) call sum_mn_name(cx*sx*(+cx*Eipq(:,2,i3)+sx*Eipq(:,2,i4)),idiag_alp12cs)
          if (idiag_alp22ss/=0) call sum_mn_name(sx**2*(+cx*Eipq(:,3,i3)+sx*Eipq(:,3,i4)),idiag_alp22ss)
          if (idiag_eta12/=0) call sum_mn_name(-(-sx*Eipq(:,2,i1)+cx*Eipq(:,2,i2))*ktestfield1,idiag_eta12)
          if (idiag_eta22/=0) call sum_mn_name(-(-sx*Eipq(:,3,i1)+cx*Eipq(:,3,i2))*ktestfield1,idiag_eta22)
          if (idiag_eta12cs/=0) call sum_mn_name(-cx*sx*(-sx*Eipq(:,2,i1)+cx*Eipq(:,2,i2))*ktestfield1,idiag_eta12cs)
          if (idiag_eta22ss/=0) call sum_mn_name(-sx**2*(-sx*Eipq(:,3,i1)+cx*Eipq(:,3,i2))*ktestfield1,idiag_eta22ss)
          if (idiag_alp12_x/=0) call sum_mn_name(xx*(+cx*Eipq(:,2,i3)+sx*Eipq(:,2,i4)),idiag_alp12_x)
          if (idiag_alp22_x/=0) call sum_mn_name(xx*(+cx*Eipq(:,3,i3)+sx*Eipq(:,3,i4)),idiag_alp22_x)
          if (idiag_eta12_x/=0) call sum_mn_name(xx*(-sx*Eipq(:,2,i1)+cx*Eipq(:,2,i2))*ktestfield1,idiag_eta12_x)
          if (idiag_eta22_x/=0) call sum_mn_name(xx*(-sx*Eipq(:,3,i1)+cx*Eipq(:,3,i2))*ktestfield1,idiag_eta22_x)
          if (idiag_alp12_x2/=0) call sum_mn_name(x2*(+cx*Eipq(:,2,i3)+sx*Eipq(:,2,i4)),idiag_alp12_x2)
          if (idiag_alp22_x2/=0) call sum_mn_name(x2*(+cx*Eipq(:,3,i3)+sx*Eipq(:,3,i4)),idiag_alp22_x2)
          if (idiag_eta12_x2/=0) call sum_mn_name(x2*(-sx*Eipq(:,2,i1)+cx*Eipq(:,2,i2))*ktestfield1,idiag_eta12_x2)
          if (idiag_eta22_x2/=0) call sum_mn_name(x2*(-sx*Eipq(:,3,i1)+cx*Eipq(:,3,i2))*ktestfield1,idiag_eta22_x2)
        endif
!
!  rms values of small scales fields bpq in response to the test fields Bpq
!  Obviously idiag_b0rms and idiag_b12rms cannot both be invoked!
!  Needs modification!
!
        if (idiag_b0rms/=0) then
          call dot2(bpq(:,:,iE0),bpq2)
          call sum_mn_name(bpq2,idiag_b0rms,lsqrt=.true.)
        endif
!
        if (idiag_b11rms/=0) then
          call dot2(bpq(:,:,1),bpq2)
          call sum_mn_name(bpq2,idiag_b11rms,lsqrt=.true.)
        endif
!
        if (idiag_b21rms/=0) then
          call dot2(bpq(:,:,2),bpq2)
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
    endsubroutine get_slices_testfield
!***********************************************************************
    subroutine testfield_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!    4-oct-18/axel+nishant: adapted from testflow
!
      use General, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine testfield_before_boundary
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
      use Mpicomm, only: mpireduce_sum, mpibcast_real, mpibcast_real_arr
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (nx,3,3) :: aijtest,bijtest
      real, dimension (nx,3) :: aatest,bbtest,jjtest,uxbtest,jxbtest
      real, dimension (nx,3) :: del2Atest2,graddivatest
      integer :: jtest,j,juxb,jjxb
      logical :: headtt_save
      real :: fac
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac=1./nyzgrid
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
      uxbtestm=0.
  
      if (.not.lsoca) then
        lpenc_loc = .false.; lpenc_loc(i_uu)=.true.

        do jtest=1,njtest
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2
          juxb=iuxbtest+3*(jtest-1)
          do n=n1,n2
            do m=m1,m2
              call calc_pencils_hydro(f,p,lpenc_loc)
              call curl(f,iaxtest,bbtest)
              call cross_mn(p%uu,bbtest,uxbtest)
              if (iuxbtest/=0) f(l1:l2,m,n,juxb:juxb+2)=uxbtest
              uxbtestm(:,:,jtest)=uxbtestm(:,:,jtest)+fac*uxbtest
              headtt=.false.
            enddo
          enddo
        enddo
!
!  do communication
!
        call finalize_aver(nprocyz,23,uxbtestm)
      endif
!
!  Do the same for jxb; do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
      jxbtestm=0.
      if (.not.lsoca_jxb) then
        do jtest=1,njtest
          iaxtest=iaatest+3*(jtest-1)
          iaztest=iaxtest+2
          jjxb=ijxbtest+3*(jtest-1)
          do n=n1,n2
            do m=m1,m2
              aatest=f(l1:l2,m,n,iaxtest:iaztest)
              call gij(f,iaxtest,aijtest,1)
              call gij_etc(f,iaxtest,aatest,aijtest,bijtest,del2Atest2,graddivatest)
              call curl_mn(aijtest,bbtest,aatest)
              call curl_mn(bijtest,jjtest,bbtest)
              call cross_mn(jjtest,bbtest,jxbtest)
              if (ijxbtest/=0) f(l1:l2,m,n,jjxb:jjxb+2)=jxbtest
              jxbtestm(:,:,jtest)=jxbtestm(:,:,jtest)+fac*jxbtest
              headtt=.false.
            enddo
          enddo
        enddo
!
!  do communication 
!
        call finalize_aver(nprocyz,23,jxbtestm)
      endif
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
! reinitialize aatest periodically if requested
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
!
        if (t >= taainit) then
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
      case (1); B0test(:,1)=cx; B0test(:,2)=sx; B0test(:,3)=0.
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
      case (1); B0test(:,1)=cx; B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=sx; B0test(:,2)=0.; B0test(:,3)=0.
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
      case (1); J0test(:,1)=0.; J0test(:,2)=-ktestfield*sx; J0test(:,3)=0.
      case (2); J0test(:,1)=0.; J0test(:,2)=+ktestfield*cx; J0test(:,3)=0.
      case default; J0test(:,:)=0.
      endselect
!
    endsubroutine set_J0test_B11_B21
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
      case (1); B0test(:,1)=0.; B0test(:,2)=bamp*cx; B0test(:,3)=0.
      case (2); B0test(:,1)=0.; B0test(:,2)=bamp*sx; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*cx
      case (4); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*sx
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_B11_B22
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
      integer :: iname,inamex,inamez,inamexz
      logical :: lreset
      logical, optional :: lwrite
!
!  reset everything in case of RELOAD
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_bx0mz=0; idiag_by0mz=0; idiag_bz0mz=0
        idiag_E111z=0; idiag_E211z=0; idiag_E311z=0
        idiag_E121z=0; idiag_E221z=0; idiag_E321z=0
        idiag_alp11x=0; idiag_alp21x=0; idiag_alp12x=0; idiag_alp22x=0;
        idiag_eta11x=0; idiag_eta21x=0; idiag_eta12x=0; idiag_eta22x=0
        idiag_E10z=0; idiag_E20z=0; idiag_E30z=0; idiag_EBpq=0
        idiag_alp11=0; idiag_alp21=0; idiag_alp31=0
        idiag_alp12=0; idiag_alp22=0; idiag_alp32=0
        idiag_eta11=0; idiag_eta21=0
        idiag_eta12=0; idiag_eta22=0
        idiag_alp11cc=0; idiag_alp21sc=0; idiag_alp12cs=0; idiag_alp22ss=0
        idiag_eta11cc=0; idiag_eta21sc=0; idiag_eta12cs=0; idiag_eta22ss=0
        idiag_alp11_x=0; idiag_alp21_x=0; idiag_alp12_x=0; idiag_alp22_x=0
        idiag_eta11_x=0; idiag_eta21_x=0; idiag_eta12_x=0; idiag_eta22_x=0
        idiag_alp11_x2=0; idiag_alp21_x2=0; idiag_alp12_x2=0; idiag_alp22_x2=0
        idiag_eta11_x2=0; idiag_eta21_x2=0; idiag_eta12_x2=0; idiag_eta22_x2=0
        idiag_b0rms=0; idiag_E0rms=0
        idiag_b11rms=0; idiag_b21rms=0; idiag_b12rms=0; idiag_b22rms=0
        idiag_E11rms=0; idiag_E21rms=0; idiag_E12rms=0; idiag_E22rms=0
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
        call parse_name(iname,cname(iname),cform(iname),'alp11cc',idiag_alp11cc)
        call parse_name(iname,cname(iname),cform(iname),'alp21sc',idiag_alp21sc)
        call parse_name(iname,cname(iname),cform(iname),'alp12cs',idiag_alp12cs)
        call parse_name(iname,cname(iname),cform(iname),'alp22ss',idiag_alp22ss)
        call parse_name(iname,cname(iname),cform(iname),'eta11cc',idiag_eta11cc)
        call parse_name(iname,cname(iname),cform(iname),'eta21sc',idiag_eta21sc)
        call parse_name(iname,cname(iname),cform(iname),'eta12cs',idiag_eta12cs)
        call parse_name(iname,cname(iname),cform(iname),'eta22ss',idiag_eta22ss)
        call parse_name(iname,cname(iname),cform(iname),'alp11_x',idiag_alp11_x)
        call parse_name(iname,cname(iname),cform(iname),'alp21_x',idiag_alp21_x)
        call parse_name(iname,cname(iname),cform(iname),'alp12_x',idiag_alp12_x)
        call parse_name(iname,cname(iname),cform(iname),'alp22_x',idiag_alp22_x)
        call parse_name(iname,cname(iname),cform(iname),'eta11_x',idiag_eta11_x)
        call parse_name(iname,cname(iname),cform(iname),'eta21_x',idiag_eta21_x)
        call parse_name(iname,cname(iname),cform(iname),'eta12_x',idiag_eta12_x)
        call parse_name(iname,cname(iname),cform(iname),'eta22_x',idiag_eta22_x)
        call parse_name(iname,cname(iname),cform(iname),'alp11_x2',idiag_alp11_x2)
        call parse_name(iname,cname(iname),cform(iname),'alp21_x2',idiag_alp21_x2)
        call parse_name(iname,cname(iname),cform(iname),'alp12_x2',idiag_alp12_x2)
        call parse_name(iname,cname(iname),cform(iname),'alp22_x2',idiag_alp22_x2)
        call parse_name(iname,cname(iname),cform(iname),'eta11_x2',idiag_eta11_x2)
        call parse_name(iname,cname(iname),cform(iname),'eta21_x2',idiag_eta21_x2)
        call parse_name(iname,cname(iname),cform(iname),'eta12_x2',idiag_eta12_x2)
        call parse_name(iname,cname(iname),cform(iname),'eta22_x2',idiag_eta22_x2)
        call parse_name(iname,cname(iname),cform(iname),'b11rms',idiag_b11rms)
        call parse_name(iname,cname(iname),cform(iname),'b21rms',idiag_b21rms)
        call parse_name(iname,cname(iname),cform(iname),'b12rms',idiag_b12rms)
        call parse_name(iname,cname(iname),cform(iname),'b22rms',idiag_b22rms)
        call parse_name(iname,cname(iname),cform(iname),'b0rms',idiag_b0rms)
        call parse_name(iname,cname(iname),cform(iname),'E11rms',idiag_E11rms)
        call parse_name(iname,cname(iname),cform(iname),'E21rms',idiag_E21rms)
        call parse_name(iname,cname(iname),cform(iname),'E12rms',idiag_E12rms)
        call parse_name(iname,cname(iname),cform(iname),'E22rms',idiag_E22rms)
        call parse_name(iname,cname(iname),cform(iname),'E0rms',idiag_E0rms)
        call parse_name(iname,cname(iname),cform(iname),'EBpq',idiag_EBpq)
      enddo
!
!  check for those quantities for which we want yz-averages
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'alp11x',idiag_alp11x)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'alp21x',idiag_alp21x)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'alp12x',idiag_alp12x)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'alp22x',idiag_alp22x)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'eta11x',idiag_eta11x)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'eta21x',idiag_eta21x)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'eta12x',idiag_eta12x)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'eta22x',idiag_eta22x)
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
      enddo
!
    endsubroutine rprint_testfield

endmodule Testfield
