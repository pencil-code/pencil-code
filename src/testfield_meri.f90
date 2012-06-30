! $Id$
!***************************************************************
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

  include 'testfield.h'
!
! Slice precalculation buffers
!
  real, target, dimension (nx,ny,3) :: bb11_xy
  real, target, dimension (nx,ny,3) :: bb11_xy2
  real, target, dimension (nx,nz,3) :: bb11_xz
  real, target, dimension (ny,nz,3) :: bb11_yz
!
! Define the EMF and the Test fields here
!
  real, dimension (nx,3,njtest) :: Eipq,bpq,jpq
  real, dimension (nx,3,3) :: atilde,alpha
  real, dimension (nx,3,3,2) :: btilde,beta
!
!  spherical bessel and legendre function for setting test fields and analysis
!
  real, dimension(mx) :: j0r,n0r,dj0dr,dn0dr,atilde_denom1,btilde_denom1
  real, dimension(my) :: P1,dP1_dtheta

!  real, dimension(:,:), pointer :: geta_z
!  real, dimension(:), pointer :: eta_z
  real :: phase_testfield=.0
!
  character (len=labellen), dimension(ninit) :: initaatest='zero'
!  real, dimension (ninit) :: kx_aatest=1.,ky_aatest=1.,kz_aatest=1.
!  real, dimension (ninit) :: phasex_aatest=0.,phasez_aatest=0.
  real, dimension (ninit) :: amplaatest=0.
  integer, dimension (njtest) :: nuxb=0
  integer :: iE0=0

  ! input parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./)
  real :: taainit=0.,daainit=0.,taainit_previous=0.
  logical :: reinitialize_aatest=.false.
  logical :: lsoca=.false.,lsoca_jxb=.true.,lset_bbtest2=.false.
  logical :: luxb_as_aux=.false.,ljxb_as_aux=.false.,linit_aatest=.false.
  logical :: lignore_uxbtestm=.false., lphase_adjust=.false.
  character (len=labellen) :: itestfield='jzero-pzero'
  real :: krtf=1.,krtf1=1.
  real :: khtf=1.,khtf1=1.
  real, dimension(nx) :: xtf
  real, dimension(ny) :: ytf,csec
!  real :: lin_testfield=0.,lam_testfield=0.,om_testfield=0.,delta_testfield=0.
!  real :: delta_testfield_next=0., delta_testfield_time=0.
  integer, parameter :: mtestfield=3*njtest
  integer :: i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9
  integer :: naainit
  real :: bamp=1.,bamp1=1.,bamp12=1.
  namelist /testfield_init_pars/ &
       B_ext,initaatest, &
       amplaatest,&
       luxb_as_aux,ljxb_as_aux

  ! run parameters
  real :: etatest=0.,etatest1=0.
  real :: tau_aatest=0.,tau1_aatest=0.
  real :: ampl_fcont_aatest=1.
  real, dimension(njtest) :: rescale_aatest=0.
  logical :: ltestfield_taver=.false.
  logical :: ltestfield_linear=.false.
  logical :: llorentzforce_testfield=.false.
  logical :: lforcing_cont_aatest=.false.
  logical :: ltestfield_artifric=.false.
  logical :: ltestfield_profile_eta_z=.false.
  namelist /testfield_run_pars/ &
       B_ext,reinitialize_aatest,lsoca,lsoca_jxb, &
       lset_bbtest2,etatest,etatest1,itestfield,&
       ltestfield_taver,llorentzforce_testfield, &
       ltestfield_profile_eta_z, &
       luxb_as_aux,ljxb_as_aux,lignore_uxbtestm, &
       lforcing_cont_aatest,ampl_fcont_aatest, &
       daainit,linit_aatest,bamp, &
       rescale_aatest,tau_aatest
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_E11xy=0      ! DIAG_DOC: $E_{11xy}$
  integer :: idiag_E12xy=0      ! DIAG_DOC: $E_{12xy}$
  integer :: idiag_E13xy=0      ! DIAG_DOC: $E_{13xy}$
  integer :: idiag_E21xy=0      ! DIAG_DOC: $E_{21xy}$
  integer :: idiag_E22xy=0      ! DIAG_DOC: $E_{22xy}$
  integer :: idiag_E23xy=0      ! DIAG_DOC: $E_{23xy}$
  integer :: idiag_E31xy=0      ! DIAG_DOC: $E_{31xy}$
  integer :: idiag_E32xy=0      ! DIAG_DOC: $E_{32xy}$
  integer :: idiag_E33xy=0      ! DIAG_DOC: $E_{33xy}$
  integer :: idiag_E41xy=0      ! DIAG_DOC: $E_{41xy}$
  integer :: idiag_E42xy=0      ! DIAG_DOC: $E_{42xy}$
  integer :: idiag_E43xy=0      ! DIAG_DOC: $E_{43xy}$
  integer :: idiag_E51xy=0      ! DIAG_DOC: $E_{51xy}$
  integer :: idiag_E52xy=0      ! DIAG_DOC: $E_{52xy}$
  integer :: idiag_E53xy=0      ! DIAG_DOC: $E_{53xy}$
  integer :: idiag_E61xy=0      ! DIAG_DOC: $E_{61xy}$
  integer :: idiag_E62xy=0      ! DIAG_DOC: $E_{62xy}$
  integer :: idiag_E63xy=0      ! DIAG_DOC: $E_{63xy}$
  integer :: idiag_E71xy=0      ! DIAG_DOC: $E_{71xy}$
  integer :: idiag_E72xy=0      ! DIAG_DOC: $E_{72xy}$
  integer :: idiag_E73xy=0      ! DIAG_DOC: $E_{73xy}$
  integer :: idiag_E81xy=0      ! DIAG_DOC: $E_{81}$
  integer :: idiag_E82xy=0      ! DIAG_DOC: $E_{82}$
  integer :: idiag_E83xy=0      ! DIAG_DOC: $E_{83}$
  integer :: idiag_E91xy=0      ! DIAG_DOC: $E_{91}$
  integer :: idiag_E92xy=0      ! DIAG_DOC: $E_{92}$
  integer :: idiag_E93xy=0      ! DIAG_DOC: $E_{93}$
! The 27 coefficients in blocks of 3 
  integer :: idiag_a11xy=0      ! DIAG_DOC: $\alpha_{11}$
  integer :: idiag_a12xy=0      ! DIAG_DOC: $\alpha_{12}$
  integer :: idiag_a13xy=0      ! DIAG_DOC: $\alpha_{13}$
  integer :: idiag_a21xy=0      ! DIAG_DOC: $\alpha_{21}$
  integer :: idiag_a22xy=0      ! DIAG_DOC: $\alpha_{22}$
  integer :: idiag_a23xy=0      ! DIAG_DOC: $\alpha_{23}$
  integer :: idiag_a31xy=0      ! DIAG_DOC: $\alpha_{31}$
  integer :: idiag_a32xy=0      ! DIAG_DOC: $\alpha_{32}$
  integer :: idiag_a33xy=0      ! DIAG_DOC: $\alpha_{33}$
!
  integer :: idiag_b111xy=0      ! DIAG_DOC: $\b_{111}$ 
  integer :: idiag_b121xy=0      ! DIAG_DOC: $\b_{121}$ 
  integer :: idiag_b131xy=0      ! DIAG_DOC: $\b_{131}$ 
  integer :: idiag_b211xy=0      ! DIAG_DOC: $\b_{211}$ 
  integer :: idiag_b221xy=0      ! DIAG_DOC: $\b_{221}$ 
  integer :: idiag_b231xy=0      ! DIAG_DOC: $\b_{231}$ 
  integer :: idiag_b311xy=0      ! DIAG_DOC: $\b_{311}$ 
  integer :: idiag_b321xy=0      ! DIAG_DOC: $\b_{321}$ 
  integer :: idiag_b331xy=0      ! DIAG_DOC: $\b_{331}$ 
!
  integer :: idiag_b112xy=0      ! DIAG_DOC: $\b_{112}$ 
  integer :: idiag_b122xy=0      ! DIAG_DOC: $\b_{122}$ 
  integer :: idiag_b132xy=0      ! DIAG_DOC: $\b_{132}$ 
  integer :: idiag_b212xy=0      ! DIAG_DOC: $\b_{212}$ 
  integer :: idiag_b222xy=0      ! DIAG_DOC: $\b_{222}$ 
  integer :: idiag_b232xy=0      ! DIAG_DOC: $\b_{232}$ 
  integer :: idiag_b312xy=0      ! DIAG_DOC: $\b_{312}$ 
  integer :: idiag_b322xy=0      ! DIAG_DOC: $\b_{322}$ 
  integer :: idiag_b332xy=0      ! DIAG_DOC: $\b_{332}$ 
!
!  arrays for azimuthally averaged uxb and jxb
!
  real, dimension (mx,my,3,mtestfield/3) :: uxbtestm,jxbtestm

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
!  Set first and last index of test field
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
      use Mpicomm, only : stop_it
      use FArrayManager
      use SharedVariables, only : get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, intent(in) :: lstarting
      integer :: jtest, ierr
! 
! Stop the code if we are not in spherical coordinates
!
      if (.not.lspherical_coords) & 
        call stop_it('initialize_testfiled: testfield_meri works only in spherical coordinates') 
!
!  Precalculate etatest if 1/etatest (==etatest1) is given instead
!
      if (etatest1/=0.) then
        etatest=1./etatest1
      endif
      if (lroot) print*,'initialize_testfield: etatest=',etatest
!
!  set cosine and sine function for setting test fields and analysis
!
      j0r= sin(x)/x
      n0r= -cos(x)/x
      dj0dr= cos(x)/x - sin(x)/(x*x)
      dn0dr= -sin(x)/x + cos(x)/(x*x)
      P1 = cos(y)
      dP1_dtheta = -sin(y)
      atilde_denom1=j0r*dn0dr-n0r*dj0dr
      btilde_denom1=dn0dr*j0r-dj0dr*n0r
!
!  debug output
!
      if (lroot.and.ip<14) then
        print*,'j0r=',j0r
        print*,'n0r',n0r
        print*,'P1',P1
      endif
!
!  calculate inverse testfield amplitude (unless it is set to zero)
!
!REVISIT
      if (bamp==0.) then
        bamp1=1.
        bamp12=1.
      else
        bamp1=1./bamp
        bamp12=1./bamp**2
      endif
!
!  calculate iE0; set ltestfield_linear in specific cases
!
      ltestfield_linear=.false.
      if (.not.lstarting) then
        select case (itestfield)
        case ('j0-P1'); iE0=0
        case ('SRSRC07'); iE0=0
        case('harmonic') 
          xtf=krtf*(x(l1:l2)-x(l1))/Lx
          if (krtf /= 0.) krtf1=Lx/krtf
          ytf=khtf*y(m1:m2);
          csec=1./sin(ytf)
          if (khtf /= 0.) khtf1=1./khtf
        case default
          call fatal_error('initialize_testfield','undefined itestfield value')
        endselect
      endif
!REVISIT
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
        write(1,'(a,i1)') 'lsoca='  ,merge(1,0,lsoca)
        write(1,'(a,i1)') 'lsoca_jxb='  ,merge(1,0,lsoca_jxb)
        write(1,'(3a)') "itestfield='",trim(itestfield)//"'"
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
      use Hydro, only: uumxy,lcalc_uumeanxy
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
!      real, dimension (nx,3,njtest) :: Eipq,bpq,jpq
      real, dimension (nx,3) :: del2Atest,uufluct
      real, dimension (nx,3) :: del2Atest2,graddiv_atest,aatest,jjtest,jxbrtest
      real, dimension (nx,3,3) :: aijtest,bijtest,Mijtest
      real, dimension (nx) :: jbpq,bpq2,Epq2,s2kzDF1,s2kzDF2,divatest,unity=1.,&
                              temp
      integer :: jtest, j, jaatest, iuxtest, iuytest, iuztest
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
      if (lcalc_uumeanxy) then
        do j=1,3
          uufluct(:,j)=p%uu(:,j)-uumxy(l1:l2,m,j)
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
!
! aatest
!
        aatest=f(l1:l2,m,n,iaxtest:iaztest)
!
! aijtest
!
        call gij(f,iaxtest,aijtest,1)
!
! bbtest
! 
       call curl_mn(aijtest,bbtest,aatest)
!
! bijtest,gradfiv_atest
! 
       call gij_etc(f,iaxtest,aatest,aijtest,bijtest,GRADDIV=graddiv_atest)
!
! jjtest
!
       call curl_mn(bijtest,jjtest,bbtest)
!
! del2A
!
        del2Atest=graddiv_atest-jjtest
!
! Select which testfield to use.
!
        select case (itestfield)
          case ('j0-P1'); call set_bbtest_j0_P1(B0test,jtest)
!
! Simplest test fields not obeying solenoidal condition from 
! Table 1 of Schrinner et al. (2007)
!
          case ('SRSRC07'); call set_bbtest_srsrc07(B0test,jtest)
          case ('harmonic'); call set_bbtest_harmonic(B0test,jtest)
          case ('B=0'); B0test=0.
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
!  Compute u=U-Ubar
!
        call cross_mn(uufluct,B0test,uxB)
        df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+uxB
!
!  Add non-SOCA terms:
!  Use f-array for uxb (if space has been allocated for this) and
!  if we don't test (i.e. if ltest_uxb=.false.).
!
        if (.not.lsoca) then
           call cross_mn(p%uu,bbtest,uxbtest)
!
!  subtract average emf, unless we ignore the <uxb> term (lignore_uxbtestm=T)
!
           if (lignore_uxbtestm) then
              duxbtest(:,:)=uxbtest(:,:)
           else
              do j=1,3
                 duxbtest(:,j)=uxbtest(:,j)-uxbtestm(l1:l2,m,j,jtest)
              enddo
           endif
!
!  add to advance test-field equation
!
           df(l1:l2,m,n,iaxtest:iaztest)=df(l1:l2,m,n,iaxtest:iaztest)+duxbtest
        endif
!
!  calculate alpha, begin by calculating uxbtest (if not already done above)
!
        if ((ldiagnos.or.l1davgfirst).and. &
          (lsoca.or.ltest_uxb)) then
          call cross_mn(p%uu,bbtest,uxbtest)
       endif
       bpq(:,:,jtest)=bbtest
       Eipq(:,:,jtest)=uxbtest*bamp1
    enddo
!
!  diffusive time step, just take the max of diffus_eta (if existent)
!  and whatever is calculated here
!
!DM check if the following is correct in spherical coordinates
      if (lfirst.and.ldt) then
        diffus_eta=max(diffus_eta,etatest*dxyz_2)
      endif
!
!  in the following block, we have already swapped the 4-6 entries with 7-9
!  The g95 compiler doesn't like to see an index that is out of bounds,
!  so prevent this warning by writing i3=3 and i4=4
!
      if (ldiagnos) then
        if (idiag_E11xy/=0) call zsum_mn_name_xy(Eipq(:,1,i1),idiag_E11xy)
        if (idiag_E12xy/=0) call zsum_mn_name_xy(Eipq(:,2,i1),idiag_E12xy)
        if (idiag_E13xy/=0) call zsum_mn_name_xy(Eipq(:,3,i1),idiag_E13xy)
        if (idiag_E21xy/=0) call zsum_mn_name_xy(Eipq(:,1,i2),idiag_E21xy)
        if (idiag_E22xy/=0) call zsum_mn_name_xy(Eipq(:,2,i2),idiag_E22xy)
        if (idiag_E23xy/=0) call zsum_mn_name_xy(Eipq(:,3,i2),idiag_E23xy)
        if (idiag_E31xy/=0) call zsum_mn_name_xy(Eipq(:,1,i3),idiag_E31xy)
        if (idiag_E32xy/=0) call zsum_mn_name_xy(Eipq(:,2,i3),idiag_E32xy)
        if (idiag_E33xy/=0) call zsum_mn_name_xy(Eipq(:,3,i3),idiag_E33xy)
        if (idiag_E41xy/=0) call zsum_mn_name_xy(Eipq(:,1,i4),idiag_E41xy)
        if (idiag_E42xy/=0) call zsum_mn_name_xy(Eipq(:,2,i4),idiag_E42xy)
        if (idiag_E43xy/=0) call zsum_mn_name_xy(Eipq(:,3,i4),idiag_E43xy)
        if (idiag_E51xy/=0) call zsum_mn_name_xy(Eipq(:,1,i5),idiag_E51xy)
        if (idiag_E52xy/=0) call zsum_mn_name_xy(Eipq(:,2,i5),idiag_E52xy)
        if (idiag_E53xy/=0) call zsum_mn_name_xy(Eipq(:,3,i5),idiag_E53xy)
        if (idiag_E61xy/=0) call zsum_mn_name_xy(Eipq(:,1,i6),idiag_E61xy)
        if (idiag_E62xy/=0) call zsum_mn_name_xy(Eipq(:,2,i6),idiag_E62xy)
        if (idiag_E63xy/=0) call zsum_mn_name_xy(Eipq(:,3,i6),idiag_E63xy)
        if (idiag_E71xy/=0) call zsum_mn_name_xy(Eipq(:,1,i7),idiag_E71xy)
        if (idiag_E72xy/=0) call zsum_mn_name_xy(Eipq(:,2,i7),idiag_E72xy)
        if (idiag_E73xy/=0) call zsum_mn_name_xy(Eipq(:,3,i7),idiag_E73xy)
        if (idiag_E81xy/=0) call zsum_mn_name_xy(Eipq(:,1,i8),idiag_E81xy)
        if (idiag_E82xy/=0) call zsum_mn_name_xy(Eipq(:,2,i8),idiag_E82xy)
        if (idiag_E83xy/=0) call zsum_mn_name_xy(Eipq(:,3,i8),idiag_E83xy)
        if (idiag_E91xy/=0) call zsum_mn_name_xy(Eipq(:,1,i9),idiag_E91xy)
        if (idiag_E92xy/=0) call zsum_mn_name_xy(Eipq(:,2,i9),idiag_E92xy)
        if (idiag_E93xy/=0) call zsum_mn_name_xy(Eipq(:,3,i9),idiag_E93xy)
!
! Invert the testfield equations here to get the transport coeffecients,
! in terms of the testfields and the calculated EMF. 
! ( \tilde{a} (Schrinner et al 2007 ArXiv:astro-ph/0609752
!  http://arxiv.org/abs/astro-ph/0609752
!  see also notes in tex/notes/testfield/spherical.tex )
!
        call invert_testfield_eqn
!
! Now calculate the (a,b) from tilde (a,b) 
!
        call get_ab_from_tildeab
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
    subroutine invert_testfield_eqn
! 
! Invert the testfield equations to get the 'tilde'-(a,b) in terms of
! EMF and testfields. For different choice of testfield different 
! subroutines are called. 
! 
!  dhruba+piyali: 
! 
      select case (itestfield)
      case ('j0-P1'); call invert_bbtest_j0_P1
      case ('SRSRC07'); call invert_bbtest_srsrc07
      case ('harmonic'); call invert_bbtest_harmonic
      case ('B=0'); call fatal_error('invert_testfield_eqn','cannot invert for zero testfield')
      case default
        call fatal_error('invert_testfield_eqn','undefined itestfield value')
      endselect
!
    endsubroutine invert_testfield_eqn
!***********************************************************************
    subroutine invert_bbtest_srsrc07
! 
! Inversion for the testfield in Schrinner '07 paper.
! 
!  dhruba+piyali: 
! 
      use Cdata
      integer :: ivec
!
      do ivec=1,3
! For the testfield (1,0,0)
        atilde(:,ivec,1) = Eipq(:,ivec,i1)
! For the testfield (0,1,0)
        atilde(:,ivec,2) = Eipq(:,ivec,i2)
! For the testfield (0,0,1)
        atilde(:,ivec,3) = Eipq(:,ivec,i3)
! For the testfield (r,0,0)
        do ix=1,nx
          btilde(ix,ivec,1,1) = Eipq(ix,ivec,i4) - x(ix+3)*atilde(ix,ivec,1)
! For the testfield (0,r,0)
          btilde(ix,ivec,2,1) = Eipq(ix,ivec,i5) - x(ix+3)*atilde(ix,ivec,2)
! For the testfield (0,0,r)
          btilde(ix,ivec,3,1) = Eipq(ix,ivec,i6) - x(ix+3)*atilde(ix,ivec,3)
! For the testfield (theta,0,0)
          btilde(ix,ivec,1,2) = x(ix+3)*(Eipq(ix,ivec,i7) - y(m)*atilde(ix,ivec,1))
! For the testfield (0,theta,0)
          btilde(ix,ivec,2,2) = x(ix+3)*(Eipq(ix,ivec,i8) - y(m)*atilde(ix,ivec,2))
! For the testfield (0,0,theta)
          btilde(ix,ivec,3,2) = x(ix+3)*(Eipq(ix,ivec,i9) - y(m)*atilde(ix,ivec,3))
        enddo
      enddo
!
    endsubroutine invert_bbtest_srsrc07
!***********************************************************************
    subroutine invert_bbtest_harmonic
! 
! Inversion for the harmonic testfield in r and \theta.
! 
!  dhruba+piyali: 
! 
      use Cdata
      integer :: ivec
!
      do ivec=1,3
        atilde(:,ivec,1) = Eipq(:,ivec,i1)*cos(xtf)+&
                           Eipq(:,ivec,i2)*sin(xtf)
        atilde(:,ivec,2) = Eipq(:,ivec,i3)*cos(xtf)+&
                           Eipq(:,ivec,i4)*sin(xtf)
        atilde(:,ivec,3) = Eipq(:,ivec,i5)*cos(xtf)+&
                           Eipq(:,ivec,i6)*sin(xtf)
        do ix=1,nx
          btilde(ix,ivec,1,1) =krtf1*(-Eipq(ix,ivec,i1)*sin(xtf(ix))+&
                                       Eipq(ix,ivec,i2)*cos(xtf(ix)))
          btilde(ix,ivec,2,1) =krtf1*(-Eipq(ix,ivec,i3)*sin(xtf(ix))+&
                                       Eipq(ix,ivec,i4)*cos(xtf(ix)))
          btilde(ix,ivec,3,1) =krtf1*(-Eipq(ix,ivec,i5)*sin(xtf(ix))+&
                                       Eipq(ix,ivec,i6)*cos(xtf(ix)))
          btilde(ix,ivec,1,2) = khtf1*x(ix)*csec(m)*(atilde(ix,ivec,1)*&
                                cos(ytf(m))-Eipq(ix,ivec,i7))
          btilde(ix,ivec,2,2) = khtf1*x(ix)*csec(m)*(atilde(ix,ivec,2)*&
                                cos(ytf(m))-Eipq(ix,ivec,i8))
          btilde(ix,ivec,3,2) = khtf1*x(ix)*csec(m)*(atilde(ix,ivec,3)*&
                                cos(ytf(m))-Eipq(ix,ivec,i9))
        enddo
      enddo
!
    endsubroutine invert_bbtest_harmonic
!***********************************************************************
    subroutine invert_bbtest_j0_P1
      use Cdata
! 
! Inversion for the testfield for spherical bessel and legendre.
! 
!  dhruba+piyali: 
! 
      integer :: ivec
!
      call fatal_error('invert_bbtest_j0_P1','not coded yet')
! get inspiration from below
!          temp=(dn0dr(l1:l2)*Eipq(:,1,i1)-dj0dr(l1:l2)*Eipq(:,1,i2))/(atilde_denom1(l1:l2))
!          if (idiag_a11xy/=0) call zsum_mn_name_xy(temp,idiag_a11xy)
!          temp=(dn0dr(l1:l2)*Eipq(:,2,i1)-dj0dr(l1:l2)*Eipq(:,2,i2))/(atilde_denom1(l1:l2))
!          if (idiag_a21xy/=0) call zsum_mn_name_xy(temp,idiag_a21xy)
!          temp=(dn0dr(l1:l2)*Eipq(:,3,i1)-dj0dr(l1:l2)*Eipq(:,3,i2))/(atilde_denom1(l1:l2))
!          if (idiag_a31xy/=0) call zsum_mn_name_xy(temp,idiag_a31xy)
! \tilde{b}
!          temp=(n0r(l1:l2)*Eipq(:,1,i1)-j0r(l1:l2)*Eipq(:,1,i2))/(btilde_denom1(l1:l2))
!          if (idiag_b111xy/=0) call zsum_mn_name_xy(temp,idiag_b111xy)
!          temp=(n0r(l1:l2)*Eipq(:,2,i1)-j0r(l1:l2)*Eipq(:,2,i2))/(btilde_denom1(l1:l2))
!          if (idiag_b211xy/=0) call zsum_mn_name_xy(temp,idiag_b211xy)
!          temp=(n0r(l1:l2)*Eipq(:,3,i1)-j0r(l1:l2)*Eipq(:,3,i2))/(btilde_denom1(l1:l2))
!          if (idiag_b311xy/=0) call zsum_mn_name_xy(temp,idiag_b311xy)
!        endif
!
    endsubroutine invert_bbtest_j0_P1
!***********************************************************************
    subroutine get_ab_from_tildeab
      use Cdata
! 
! Get the a and b (alpha and beta in our notation) from the tilde (a,b)
! Eq. (16) page 6, Schrinner 2007 (Arxiv version)
! 
!  dhruba+piyali: 
!
!
      do ix=1, nx 
        alpha(ix,:,1) = atilde(ix,:,1) - btilde(ix,:,2,2)/x(ix+3)
        alpha(ix,:,2) = atilde(ix,:,2) - btilde(ix,:,1,2)/x(ix+3)
      enddo
      alpha(:,:,3) = atilde(:,:,3) 
      beta = btilde

    endsubroutine get_ab_from_tildeab
!***********************************************************************
    subroutine get_slices_testfield(f,slices)
! 
!  Write slices for animation of magnetic variables.
! 
!  12-sep-09/axel: adapted from the corresponding magnetic routine
! 
      use General, only: keep_compiler_quiet
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
!  calculate <uxb>, which is needed when lsoca=.false.
!
!  21-jan-06/axel: coded
!
      use Cdata
      use Sub
      use Hydro, only: uumxy,lcalc_uumeanxy
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (mx,my,3) :: uxbtestm_temp,jxbtestm_temp
!
      real, dimension (nx,3,3) :: aijtest,bijtest
      real, dimension (nx,3) :: aatest,bbtest,jjtest,uxbtest,jxbtest,uu,uufluct
      real, dimension (nx,3) :: del2Atest2,graddiv_atest
      integer :: jtest,j,juxb,jjxb
      logical :: headtt_save
      real :: fac
!
      intent(inout) :: f
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac=1./nzgrid
!
!  do each of the 9 test fields at a time
!  but exclude redundancies, e.g. if the averaged field lacks x extent.
!  Note: the same block of lines occurs again further up in the file.
!
      do jtest=1,njtest
        iaxtest=iaatest+3*(jtest-1)
        iaztest=iaxtest+2
        uxbtestm(:,:,:,jtest) = 0.
        jxbtestm(:,:,:,jtest)=0.
        if (.not.lsoca) then
!
          do m=m1,m2
            do n=n1,n2
              aatest=f(l1:l2,m,n,iaxtest:iaztest)
              call gij(f,iaxtest,aijtest,1)
              call curl_mn(aijtest,bbtest,aatest)
!
!  Subtract the mean flow from uu.
!  This does not matter for the calculation of <Ubar x b>, which
!  is zero, but for Ubar x b, if it is saved in the f array.
!
!                 uu=f(l1:l2,m,n,iux:iuz)
! 
              uu=f(l1:l2,m,n,iux:iuz)
              if (lcalc_uumeanxy) then
                do j=1,3
                  uufluct(:,j)=uu(:,j)-uumxy(l1:l2,m,j)
                enddo
              else
                uufluct=uu
              endif
              call cross_mn(uufluct,bbtest,uxbtest)
              juxb=iuxb+3*(jtest-1)
              if (iuxb/=0) f(l1:l2,m,n,juxb:juxb+2)=uxbtest
              do j=1,3
                uxbtestm_temp(l1:l2,m,j)=uxbtestm_temp(l1:l2,m,j)+fac*uxbtest(:,jtest)
              enddo
              if (.not.lsoca_jxb) then
                call gij_etc(f,iaxtest,aatest,aijtest,bijtest,GRADDIV=graddiv_atest)
                call curl_mn(bijtest,jjtest,bbtest)
                call cross_mn(jjtest,bbtest,jxbtest)
                jjxb=ijxb+3*(jtest-1)
                if (ijxb/=0) f(l1:l2,m,n,jjxb:jjxb+2)=jxbtest
                do j=1,3
                  jxbtestm_temp(l1:l2,m,j)=jxbtestm_temp(l1:l2,m,j)+fac*jxbtest(:,jtest)
                enddo
              endif
            enddo
          enddo
!
!  do communication in the z (or phi) direction
!
          if (nprocz>1) then
            call mpiallreduce_sum(uxbtestm_temp,uxbtestm(:,:,:,jtest),(/mx,my,3/),idir=3)
            if (.not.lsoca_jxb) &
                call mpiallreduce_sum(jxbtestm_temp,jxbtestm(:,:,:,jtest),(/mx,my,3/),idir=3)
          endif
!
        endif
      enddo
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
    subroutine set_bbtest_j0_P1(B0test,jtest)
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
      case (1); B0test(:,1)=bamp*j0r(l1:l2); B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=bamp*n0r(l1:l2); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=bamp*P1(m); B0test(:,2)=0.; B0test(:,3)=0.
!
      case (4); B0test(:,1)=0.; B0test(:,2)=bamp*j0r(l1:l2); B0test(:,3)=0.
      case (5); B0test(:,1)=0.; B0test(:,2)=bamp*n0r(l1:l2); B0test(:,3)=0.
      case (6); B0test(:,1)=0.; B0test(:,2)=bamp*P1(m); B0test(:,3)=0.
!
      case (7); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*j0r(l1:l2)
      case (8); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*n0r(l1:l2)
      case (9); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*P1(m)
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_j0_P1
!***********************************************************************
    subroutine set_bbtest_srsrc07(B0test,jtest)
!
!  set testfield
!
!   25-nov-10/piyali: copied set_bbtest_j0_P1 and modified according
!                     to Table.~1 of Schrinner et al. (2007)
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
      case (1); B0test(:,1)=bamp; B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=0.; B0test(:,2)=bamp; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp
!
      case (4); B0test(:,1)=bamp*x(l1:l2); B0test(:,2)=0.; B0test(:,3)=0.
      case (5); B0test(:,1)=0.; B0test(:,2)=bamp*x(l1:l2); B0test(:,3)=0.
      case (6); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*x(l1:l2)
!
      case (7); B0test(:,1)=bamp*y(m); B0test(:,2)=0.; B0test(:,3)=0.
      case (8); B0test(:,1)=0.; B0test(:,2)=bamp*y(m); B0test(:,3)=0.
      case (9); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*y(m)
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_srsrc07
!***********************************************************************
    subroutine set_bbtest_harmonic(B0test,jtest)
!
!  set testfield
!
!   25-nov-05/piyali: copied set_bbtest_j0_P1 and modified according
!                     to Table.~1 of Schrinner et al. (2007)
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
      case (1); B0test(:,1)=bamp*cos(xtf); B0test(:,2)=0.; B0test(:,3)=0.
      case (2); B0test(:,1)=bamp*sin(xtf); B0test(:,2)=0.; B0test(:,3)=0.
      case (3); B0test(:,1)=0.; B0test(:,2)=bamp*cos(xtf); B0test(:,3)=0.
      case (4); B0test(:,1)=0.; B0test(:,2)=bamp*sin(xtf); B0test(:,3)=0.
      case (5); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*cos(xtf)
      case (6); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*sin(xtf)
      case (7); B0test(:,1)=bamp*cos(ytf(m)); B0test(:,2)=0.; B0test(:,3)=0.
      case (8); B0test(:,1)=0.; B0test(:,2)=bamp*cos(ytf(m)); B0test(:,3)=0.
      case (9); B0test(:,1)=0.; B0test(:,2)=0.; B0test(:,3)=bamp*cos(ytf(m))
      case default; B0test(:,:)=0.
      endselect
!
    endsubroutine set_bbtest_harmonic
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
!DM work here
      endif
!
!  check for those quantities that we want to evaluate online
! 
      do iname=1,nname
!DM work here
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
!DM maybe nothing here
      enddo
!
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
