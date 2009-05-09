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
! CPARAM logical, parameter :: ltestflow = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Testflow

  use Cparam
  use Messages

  implicit none

  include 'testflow.h'
!
! Slice precalculation buffers
!
  real, target, dimension (nx,ny,3) :: uu11_xy
  real, target, dimension (nx,ny,3) :: uu11_xy2
  real, target, dimension (nx,nz,3) :: uu11_xz
  real, target, dimension (ny,nz,3) :: uu11_yz
!
!  cosine and sine function for setting test fields and analysis
!
  real, dimension(mz) :: cz,sz,k1cz,k1sz
!
  character (len=labellen), dimension(ninit) :: inituutest='nothing'
  real, dimension (ninit) :: ampluutest=0.
  real, parameter :: cs2test=1., cs2test1=1./cs2test

  ! input parameters
  real, dimension(3) :: B_ext=(/0.,0.,0./)
  real, dimension (nx,3) :: bbb
  real :: kx_uutest=1.,ky_uutest=1.,kz_uutest=1.
  real :: tuuinit=0.,duuinit=0.
  logical :: reinitialize_uutest=.false.
  logical :: zextent=.true.,lsoca_testflow=.true.
  logical :: lset_bbtest2=.false.,lset_U0test=.true.,lignore_ugutestm=.false.
  logical :: lugu_as_aux=.false.,linit_uutest=.false.
  logical :: lforcing_cont_uutest=.false.
  character (len=labellen) :: itestfield='B11-B21'
  real :: ktestfield=1., ktestfield1=1.
  integer, parameter :: ntestflow=4*njtest
  integer :: nuuinit
  real :: wamp=1.
  namelist /testflow_init_pars/ &
       B_ext,zextent,inituutest, &
       ampluutest,kx_uutest,ky_uutest,kz_uutest, &
       lugu_as_aux

  ! run parameters
  real :: nutest=0.,nutest1=0.
  namelist /testflow_run_pars/ &
       B_ext,reinitialize_uutest,zextent,lsoca_testflow, &
       lset_bbtest2,lset_U0test,lignore_ugutestm, &
       nutest,nutest1,itestfield,ktestfield, &
       lugu_as_aux,duuinit,linit_uutest,wamp, &
       lforcing_cont_uutest

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_lam11=0      ! DIAG_DOC: $\lambda_{11}$
  integer :: idiag_lam21=0      ! DIAG_DOC: $\lambda_{21}$
  integer :: idiag_lam12=0      ! DIAG_DOC: $\lambda_{12}$
  integer :: idiag_lam22=0      ! DIAG_DOC: $\lambda_{22}$
  integer :: idiag_nu11=0       ! DIAG_DOC: $\nu_{113}k$
  integer :: idiag_nu21=0       ! DIAG_DOC: $\nu_{213}k$
  integer :: idiag_nu12=0       ! DIAG_DOC: $\nu_{123}k$
  integer :: idiag_nu22=0       ! DIAG_DOC: $\nu_{223}k$
  integer :: idiag_u0rms=0      ! DIAG_DOC: $\left<u_{0}^2\right>$
  integer :: idiag_u11rms=0     ! DIAG_DOC: $\left<u_{11}^2\right>$
  integer :: idiag_u21rms=0     ! DIAG_DOC: $\left<u_{21}^2\right>$
  integer :: idiag_u12rms=0     ! DIAG_DOC: $\left<u_{12}^2\right>$
  integer :: idiag_u22rms=0     ! DIAG_DOC: $\left<u_{22}^2\right>$
  integer :: idiag_F111z=0      ! DIAG_DOC: ${\cal F}_1^{11}$
  integer :: idiag_F211z=0      ! DIAG_DOC: ${\cal F}_2^{11}$
  integer :: idiag_F311z=0      ! DIAG_DOC: ${\cal F}_3^{11}$
  integer :: idiag_F121z=0      ! DIAG_DOC: ${\cal F}_1^{21}$
  integer :: idiag_F221z=0      ! DIAG_DOC: ${\cal F}_2^{21}$
  integer :: idiag_F321z=0      ! DIAG_DOC: ${\cal F}_3^{21}$
  integer :: idiag_F112z=0      ! DIAG_DOC: ${\cal F}_1^{12}$
  integer :: idiag_F212z=0      ! DIAG_DOC: ${\cal F}_2^{12}$
  integer :: idiag_F312z=0      ! DIAG_DOC: ${\cal F}_3^{12}$
  integer :: idiag_F122z=0      ! DIAG_DOC: ${\cal F}_1^{22}$
  integer :: idiag_F222z=0      ! DIAG_DOC: ${\cal F}_2^{22}$
  integer :: idiag_F322z=0      ! DIAG_DOC: ${\cal F}_3^{22}$
  integer :: idiag_F10z=0       ! DIAG_DOC: ${\cal F}_1^{0}$
  integer :: idiag_F20z=0       ! DIAG_DOC: ${\cal F}_2^{0}$
  integer :: idiag_F30z=0       ! DIAG_DOC: ${\cal F}_3^{0}$
  integer :: idiag_ux0mz=0      ! DIAG_DOC: $\left<u_{x}\right>_{xy}$
  integer :: idiag_uy0mz=0      ! DIAG_DOC: $\left<u_{y}\right>_{xy}$
  integer :: idiag_uz0mz=0      ! DIAG_DOC: $\left<u_{z}\right>_{xy}$

  real, dimension (mz,3,0:njtest) :: unltestm
  real, dimension (mz,0:njtest) :: hnltestm

  contains

!***********************************************************************
    subroutine register_testflow()
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
!  Set first and last index of text field
!  Note that iuxtest, ..., ihhtest are being overwritten later,
!  and only iuutest stays fixed.
!
      iuutest=nvar+1
      iuxtest=iuutest
      iuytest=iuutest+1
      iuztest=iuutest+2
      ihhtest=iuutest+3
!     iuxtestpq=iuutest+4*(njtest-1)
!     iuztestpq=iuxtestpq+2
!     ihhtestpq=iuxtestpq+3
      nvar=nvar+ntestflow
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_testflow: nvar = ', nvar
        print*, 'register_testflow: iuutest = ', iuutest
      endif
!
!  Put variable names in array
!
      do j=1,ntestflow
        varname(j) = 'uutest'
      enddo
!
!  Identify version number.
!
      if (lroot) call cvs_id( &
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
    subroutine initialize_testflow(f)
!
!  Perform any post-parameter-read initialization
!
!   2-jun-05/axel: adapted from magnetic
!
      use Cdata
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Precalculate nutest if 1/nutest (==nutest1) is given instead
!
      if (nutest1/=0.) then
        nutest=1./nutest1
      endif
!
!  set to zero and then rescale the testflow
!  (in future, could call something like init_uu_simple)
!
      if (reinitialize_uutest) then
        f(:,:,:,iuutest:iuutest+ntestflow-1)=0.
      endif
!
!  set cosine and sine function for setting test fields and analysis
!
      cz=cos(ktestfield*z)
      sz=sin(ktestfield*z)
!
!  Also calculate its inverse, but only if different from zero
!
      if (ktestfield==0) then
        ktestfield1=1.
      else
        ktestfield1=1./ktestfield
      endif
!
!  cosine and sine functions multiplied with k1
!
      k1cz=ktestfield1*cos(ktestfield*z)
      k1sz=ktestfield1*sin(ktestfield*z)
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
          call farray_register_auxiliary('ugu',iugu,vector=3*njtest)
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
        open(1,file=trim(datadir)//'/testflow_info.dat',STATUS='unknown')
        write(1,'(a,i1)') 'zextent=',merge(1,0,zextent)
        write(1,'(a,i1)') 'lsoca_testflow='  ,merge(1,0,lsoca_testflow)
        write(1,'(3a)') "itestfield='",trim(itestfield)//"'"
        write(1,'(a,f5.2)') 'ktestfield=',ktestfield
        close(1)
      endif
!
    endsubroutine initialize_testflow
!***********************************************************************
    subroutine init_uutest(f)
!
!  initialise testflow; called from start.f90
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
!     real, dimension (nx,3) :: bb
      real, dimension (nx) :: b2,fact
      real :: beq2
      integer :: j
!
      do j=1,ninit

      select case(inituutest(j))

      case('zero'); f(:,:,:,iuutest:iuutest+ntestflow-1)=0.
      case('gaussian-noise-1'); call gaunoise(ampluutest(j),f,iuxtest+0,iuztest+0)
      case('gaussian-noise-2'); call gaunoise(ampluutest(j),f,iuxtest+4,iuztest+4)
      case('gaussian-noise-3'); call gaunoise(ampluutest(j),f,iuxtest+8,iuztest+8)
      case('sinwave-x-1')
        call sinwave(+ampluutest(j),f,ihhtest+0,kx=kx_uutest)
        call sinwave(+ampluutest(j),f,iuxtest+0,kx=kx_uutest)
      case('sinwave-x-2')
        call sinwave(+ampluutest(j),f,iuxtest+4,kx=kx_uutest)
        call sinwave(+ampluutest(j),f,ihhtest+4,kx=kx_uutest)
      case('nothing'); !(do nothing)

      case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'init_uutest: check inituutest: ', trim(inituutest(j))
        call stop_it("")

      endselect
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_uutest(f)
!
    endsubroutine init_uutest
!***********************************************************************
    subroutine pencil_criteria_testflow()
!
!   All pencils that the Testflow module depends on are specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
      lpenc_requested(i_uu)=.true.
      lpenc_requested(i_lnrho)=.true.
!
    endsubroutine pencil_criteria_testflow
!***********************************************************************
    subroutine pencil_interdep_testflow(lpencil_in)
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
    subroutine read_testflow_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=testflow_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testflow_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_testflow_init_pars
!***********************************************************************
    subroutine write_testflow_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=testflow_init_pars)

    endsubroutine write_testflow_init_pars
!***********************************************************************
    subroutine read_testflow_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=testflow_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=testflow_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_testflow_run_pars
!***********************************************************************
    subroutine write_testflow_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=testflow_run_pars)

    endsubroutine write_testflow_run_pars
!***********************************************************************
    subroutine duutest_dt(f,df,p)
!
!  testflow evolution:
!
!  calculate du^(pq)/dt = -( u^(0) .grad u^(pq) + u^(pq).grad u ) + < u^(0) .grad u^(pq) + u^(pq).grad u >
!        [alternatively:  -( u^(pq).grad u^(0)  + u.grad u^(pq) ) + < u^(pq).grad u^(0)  + u.grad u^(pq) >]
!                                                 - grad h^(pq)
!                         -U^(pq).grad u - u.grad U^(pq)
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

      use Cdata
      use Diagnostics
      use Sub
      use Hydro, only: uumz,lcalc_uumean
      use Mpicomm, only: stop_it
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
      real, dimension (nx,3) :: uutest,uufluct, del2utest,graddivutest,ghtest,ghfluct
      real, dimension (nx,3) :: U0test=0,gU0test=0,U0gutest,ugU0test,gH0test=0
      real, dimension (nx)   :: hhtest,H0test=0,ugH0test,U0ghtest,upq2,divutest,help
      real, dimension (nx,3,0:njtest) :: Fipq,upq
      real, dimension (nx,0:njtest) :: Qipq,hpq
      real, dimension (nx,3,3) :: uijtest

      logical :: ltestflow_out

      integer :: jtest,jfnamez,j,i,i3,i4
      integer,save :: ifirst=0
      
      character (len=5) :: ch
      character (len=130) :: file
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'duutest_dt: SOLVE'
      if (headtt) then
        if (iuxtest /= 0) call identify_bcs('uxtest',iuxtest)
        if (iuytest /= 0) call identify_bcs('uytest',iuytest)
        if (iuztest /= 0) call identify_bcs('uztest',iuztest)
        if (ihhtest /= 0) call identify_bcs('hhtest',ihhtest)
      endif
!
!  calculate uufluct=U-Umean
!
      uufluct=p%uu
      ghfluct=p%glnrho
      
      if (lcalc_uumean) then

        do j=1,3
          uufluct(:,j)=uufluct(:,j)-uumz(n,j)
        enddo
!
!  rhomz is not currently calculated....
!
!!!  ghfluct=ghfluct-alog(rhomz(n))  !MR: noch falsch!!!

      endif
!
!  do each of the njtest test flow at a time            
!  jtest=0  refers to primary turbulence (u^(0), h^(0))
!
      do jtest=0,njtest
        iuxtest=iuutest+4*jtest
        iuztest=iuxtest+2
        ihhtest=iuxtest+3
!
!  velocity vector and enthalpy
!
        uutest=f(l1:l2,m,n,iuxtest:iuztest)
        hhtest=f(l1:l2,m,n,ihhtest)

        upq(:,:,jtest)=uutest
        hpq(:  ,jtest)=hhtest
!
!  velocity gradient matrix and div u term
!                                                               
        call gij(f,iuxtest,uijtest,1)
        call div_mn(uijtest,divutest,uutest)                    !MR: better determine divutest without determining uijtest
!
!  gradient of (pseudo) enthalpy
!
        call grad(f,ihhtest,ghtest)

        if (jtest.gt.0) then
          select case(itestfield)                       ! get testfield U^(pq), gradU^(pq), H^(pq), gradH^(pq)
            case('W11-W22'); call set_U0test_W11_W22(U0test,gU0test,H0test,gH0test,jtest) 
            case('W=0'); U0test=0; gU0test=0; H0test=0.;gH0test=0.
            case default
              call fatal_error('duutest_dt','undefined itestfield value')
          endselect
        endif
!
!  rhs of continuity equation (nonlinear terms already in df!), dh^pq/dt = n.l.Terms - cs^2*div u^pq -u.gradH^pq - U^pq.gradh,
!
        df(l1:l2,m,n,ihhtest)=df(l1:l2,m,n,ihhtest) - cs2test*divutest
!
!  testfield inhomogeneity
!
        if (jtest.gt.0) then

          if (lset_U0test) then
            call dot_mn(uufluct,gH0test,ugH0test)
            call dot_mn(U0test,ghfluct,U0ghtest)
            df(l1:l2,m,n,ihhtest)=df(l1:l2,m,n,ihhtest)-ugH0test-U0ghtest
          endif

        endif
!
!  rhs of momentum equation (nonlinear terms already in df!)
!
        df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) - ghtest
!
!  testfield inhomogeneity
!
        if ( jtest.gt.0 ) then
          if (lset_U0test) then
!
            call h_dot_grad(U0test,p%uij,uufluct,U0gutest)
            call multsv(uufluct(:,3),gU0test,ugU0test,ladd=.true.) 

            df(l1:l2,m,n,iuxtest:iuztest) = df(l1:l2,m,n,iuxtest:iuztest) - U0gutest-ugU0test
!
            call multmv(p%sij,gH0test,U0gutest)
!
            call multsv(ghfluct(:,3),gU0test,ugU0test)
            call multsv(-(1./3.)*gU0test(:,3),ghfluct,ugU0test,ladd=.true.)
            ghfluct(:,3)=0.       !MR: from here on ghfluct no longer valid 
            call dot_mn(gU0test,ghfluct,help)
            ugU0test(:,3) = ugU0test(:,3) + help

            df(l1:l2,m,n,iuxtest:iuztest) = df(l1:l2,m,n,iuxtest:iuztest) + 2.*nutest*cs2test1*(U0gutest+ugU0test)
!
        endif
      endif
!
!  add linear part of diffusion term nu*(del2u^pq + divu^pq/3)
!
        if (nutest/=0.) then

          call del2v_etc(f,iuxtest,DEL2=del2utest,GRADDIV=graddivutest)
          df(l1:l2,m,n,iuxtest:iuztest)=df(l1:l2,m,n,iuxtest:iuztest) &
            +nutest*(del2utest+(1./3.)*graddivutest)

        endif
!
!  calculate mean force and mean source
!  (this will not be ok on multiproceesor jobs)
!
        do j=1,3
          Fipq(:,j,jtest)=unltestm(n,j,jtest)/wamp
        enddo
        Qipq(:,  jtest)=hnltestm(n,  jtest)/wamp
!
!  check for testflow timestep
!
        if (lfirst.and.ldt) then
          advec_uu=max(advec_uu,abs(uutest(:,1))*dx_1(l1:l2)+ &
                                abs(uutest(:,2))*dy_1(  m  )+ &
                                abs(uutest(:,3))*dz_1(  n  ))
        endif
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
!  in the following block, we have already swapped the 4-6 entries with 7-9
!  The g95 compiler doesn't like to see an index that is out of bounds,
!  so prevent this warning by writing i3=3 and i4=4
!
      i3=3
      i4=4
      if (ldiagnos) then
        if (idiag_ux0mz/=0) call xysum_mn_name_z(upq(:,1,i3),idiag_ux0mz)
        if (idiag_uy0mz/=0) call xysum_mn_name_z(upq(:,2,i3),idiag_uy0mz)
        if (idiag_uz0mz/=0) call xysum_mn_name_z(upq(:,3,i3),idiag_uz0mz)
        if (idiag_F111z/=0) call xysum_mn_name_z(Fipq(:,1,1),idiag_F111z)
        if (idiag_F211z/=0) call xysum_mn_name_z(Fipq(:,2,1),idiag_F211z)
        if (idiag_F311z/=0) call xysum_mn_name_z(Fipq(:,3,1),idiag_F311z)
        if (idiag_F121z/=0) call xysum_mn_name_z(Fipq(:,1,2),idiag_F121z)
        if (idiag_F221z/=0) call xysum_mn_name_z(Fipq(:,2,2),idiag_F221z)
        if (idiag_F321z/=0) call xysum_mn_name_z(Fipq(:,3,2),idiag_F321z)
        if (idiag_F112z/=0) call xysum_mn_name_z(Fipq(:,1,i3),idiag_F112z)
        if (idiag_F212z/=0) call xysum_mn_name_z(Fipq(:,2,i3),idiag_F212z)
        if (idiag_F312z/=0) call xysum_mn_name_z(Fipq(:,3,i3),idiag_F312z)
        if (idiag_F122z/=0) call xysum_mn_name_z(Fipq(:,1,i4),idiag_F122z)
        if (idiag_F222z/=0) call xysum_mn_name_z(Fipq(:,2,i4),idiag_F222z)
        if (idiag_F322z/=0) call xysum_mn_name_z(Fipq(:,3,i4),idiag_F322z)
        if (idiag_F10z/=0) call xysum_mn_name_z(Fipq(:,1,i3),idiag_F10z)
        if (idiag_F20z/=0) call xysum_mn_name_z(Fipq(:,2,i3),idiag_F20z)
        if (idiag_F30z/=0) call xysum_mn_name_z(Fipq(:,3,i3),idiag_F30z)
!
!  lambda and nu tensors
!
        if (idiag_lam11/=0) call sum_mn_name(+cz(n)*Fipq(:,1,1)+sz(n)*Fipq(:,1,2),idiag_lam11)
        if (idiag_lam21/=0) call sum_mn_name(+cz(n)*Fipq(:,2,1)+sz(n)*Fipq(:,2,2),idiag_lam21)
        if (idiag_nu11/=0) call sum_mn_name((-sz(n)*Fipq(:,1,1)+cz(n)*Fipq(:,1,2))*ktestfield1,idiag_nu11)
        if (idiag_nu21/=0) call sum_mn_name((-sz(n)*Fipq(:,2,1)+cz(n)*Fipq(:,2,2))*ktestfield1,idiag_nu21)
!
!  print warning if lam12 and lam21 are needed, but njtest is too small
!
        if ((idiag_lam12/=0.or.idiag_lam22/=0 &
         .or.idiag_nu12/=0.or.idiag_nu22/=0).and.njtest<=2) then
          call stop_it('njtest is too small if lam12, lam22, nu12, or nu22 are needed')
        else
          if (idiag_lam12/=0) call sum_mn_name(+cz(n)*Fipq(:,1,i3)+sz(n)*Fipq(:,1,i4),idiag_lam12)
          if (idiag_lam22/=0) call sum_mn_name(+cz(n)*Fipq(:,2,i3)+sz(n)*Fipq(:,2,i4),idiag_lam22)
          if (idiag_nu12/=0) call sum_mn_name((-sz(n)*Fipq(:,1,i3)+cz(n)*Fipq(:,1,i4))*ktestfield1,idiag_nu12)
          if (idiag_nu22/=0) call sum_mn_name((-sz(n)*Fipq(:,2,i3)+cz(n)*Fipq(:,2,i4))*ktestfield1,idiag_nu22)
        endif
!
!  rms values of small scales fields upq in response to the test fields Upq
!  Obviously idiag_u0rms and idiag_u12rms cannot both be invoked!
!  Needs modification!
!
        if (idiag_u0rms/=0) then
          call dot2(upq(:,:,i3),upq2)
          call sum_mn_name(upq2,idiag_u0rms,lsqrt=.true.)
        endif
!
        if (idiag_u11rms/=0) then
          call dot2(upq(:,:,1),upq2)
          call sum_mn_name(upq2,idiag_u11rms,lsqrt=.true.)
        endif
!
        if (idiag_u21rms/=0) then
          call dot2(upq(:,:,2),upq2)
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
      endif
!
!  write utest-slices for output in wvid in run.f90
!  Note: ix is the index with respect to array with ghost zones.
! 
      if (lvideo.and.lfirst) then
        do j=1,3
          uu11_yz(m-m1+1,n-n1+1,j)=upq(ix_loc-l1+1,j,1)
          if (m==iy_loc)  uu11_xz(:,n-n1+1,j)=upq(:,j,1)
          if (n==iz_loc)  uu11_xy(:,m-m1+1,j)=upq(:,j,1)
          if (n==iz2_loc) uu11_xy2(:,m-m1+1,j)=upq(:,j,1)
        enddo
      endif
!
! initialize uutest periodically if requested
!
      if (linit_uutest) then
         file=trim(datadir)//'/tinit_uutest.dat'
         if (ifirst==0) then
            call read_snaptime(trim(file),tuuinit,nuuinit,duuinit,t)
            if (tuuinit==0 .or. tuuinit < t-duuinit) then
              tuuinit=t+duuinit
            endif
            ifirst=1
         endif
!
         if (t >= tuuinit) then
            reinitialize_uutest=.true.
            call initialize_testflow(f)
            reinitialize_uutest=.false.
            call update_snaptime(file,tuuinit,nuuinit,duuinit,t,ltestflow_out,ch,.false.)
         endif
      endif
!
    endsubroutine duutest_dt
!***********************************************************************
    subroutine get_slices_testflow(f,slices)
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
        case ('uu11')
          if (slices%index>=3) then
            slices%ready=.false.
          else
            slices%index=slices%index+1
            slices%yz =>uu11_yz(:,:,slices%index)
            slices%xz =>uu11_xz(:,:,slices%index)
            slices%xy =>uu11_xy(:,:,slices%index)
            slices%xy2=>uu11_xy2(:,:,slices%index)
            if (slices%index<3) slices%ready=.true.
          endif
      endselect
!
    endsubroutine get_slices_testflow
!***********************************************************************
    subroutine calc_ltestflow_nonlin_terms(f,df)
!
!  calculates <-u0.gradu0 + (2nu/cs^2)*gradh0.Sij0>, <u0.gradh0>,
!             <-u0.grad utest - utest.grad u + (2nu/cs^2)*(grad h0.Sijtest + grad htest.Sij)>, 
!             <u0.grad htest + utest.gradh>
!  which is needed when lsoca_testflow=.false., lsoca_testflow=.false.,  resp.
!  this is done prior to the pencil loop (calc_fluct=false)
!
!  calculates (-u0.gradu0 + (2nu/cs^2)*gradh0.Sij0)', (u0.gradh0)',
!             (-u0.grad utest - utest.grad u) + (2nu/cs^2)*(grad h0.Sijtest + grad htest.Sij)', 
!             (u0.grad htest + utest.gradh)'

!  which is needed when lsoca_unl=.false.,  resp. 
!  this is done inside the pencil loop 
!
!  15-mar-08/axel: coded
!  24-jun-08/MR: modified
!
      use Cdata
      use Sub
      use Hydro, only: uumz,lcalc_uumean
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f 
      real, dimension (mx,my,mz,mvar) :: df
!
      intent(inout) :: df
      intent(in)    :: f
  
      real, dimension (nz,nprocz,3,0:njtest) :: unltestm1
      real, dimension (nz,nprocz,0:njtest  ) :: hnltestm1
  
      real, dimension (nz,nprocz,3,0:njtest) :: unltestm1_tmp  !MR: richtig? ?Axel
      real, dimension (nz,nprocz,0:njtest  ) :: hnltestm1_tmp
 
      real, dimension (nx,3)   :: uufluct,uutest, uu0, ghfluct, ghtest, gh0, sghtest, unltest
      real, dimension (nx,3,3) :: sijtest,uijtest,sij0,uij0
      real, dimension (nx)     :: hhfluct,divutest,hnltest
 
      integer :: jtest,i,j,nxy=nxgrid*nygrid,ju,jugu,iux0
      logical :: headtt_save
      real :: fac
!
!  In this routine we will reset headtt after the first pencil,
!  so we need to reset it afterwards.
!
      headtt_save=headtt
      fac=1./nxy

!  do each of the njtest (=2, 4, or 9) test flows at a time
!  but exclude redundancies, e.g. if the averaged flows lacks x extent.
!
      do jtest=0,njtest
    
        iuxtest=iuutest+4*jtest
        iuztest=iuxtest+2
        ihhtest=iuxtest+3
  
        do n=n1,n2
  
          unltestm(n,:,jtest)=0.
          hnltestm(n,jtest)=0.

          do m=m1,m2
!
!  calculate uufluct=U-Umean
!
            uufluct=f(l1:l2,m,n,iux:iux+2)
            hhfluct=f(l1:l2,m,n,iux+3)

            if (lcalc_uumean) then
  
              do j=1,3
                uufluct(:,j)=uufluct(:,j)-uumz(n,j)
              enddo

!!!           hhfluct=hhfluct-alog(rhomz(n))     !MR: noch falsch!!!

            endif
!
!  velocity vector and enthalpy gradient
!
            uutest=f(l1:l2,m,n,iuxtest:iuztest)
            call grad(f,ihhtest,ghtest)
!
!  velocity gradient matrix and velocity divergence
!  
            call gij(f,iuxtest,uijtest,1)
            call div_mn(uijtest,divutest,uutest)
!
!  calculate stress tensor sijtest

            do j=1,3
  
              do i=j,3
                sijtest(:,i,j)=.5*(uijtest(:,i,j)+uijtest(:,j,i))
                if ( i.ne.j ) sijtest(:,j,i)=sijtest(:,i,j)
              enddo

              sijtest(:,j,j)=sijtest(:,j,j)-(1./3.)*divutest

            enddo

            if ( jtest.eq.0 ) then    ! primary turbulence
  
              gh0  = ghtest           ! save primary turbulence
              uu0  = uutest
              sij0 = sijtest
              iux0 = iuxtest

!  u.gradu term and u.gradh term

              call u_dot_grad(f,iuxtest,uijtest,uutest,unltest)
              call dot_mn(uutest,ghtest,hnltest)

            elseif ( .not.lsoca_testflow ) then

              call u_dot_grad(f,iuxtest,uijtest,uufluct,unltest)
              call u_dot_grad(f,iux0   ,uij0   ,uutest ,unltest,ladd=.true.)

              call dot_mn(uu0   ,ghtest ,hnltest)
              call dot_mn(uutest,ghfluct,hnltest,ladd=.true.)

            endif
!
!  add nonlinear part of diffusion term nu*2*Sgradh/cs2
!
            if (nutest/=0.) then

!  calculate stress tensor sij from uij

              if ( jtest.eq.0 ) then

                call multmv(sijtest,ghtest,sghtest) 

              elseif ( .not.lsoca_testflow ) then
 
                call multmv(sijtest,ghfluct,sghtest)
                call multmv(sij0   ,ghtest ,sghtest, ladd=.true.)
 
              endif

              unltest = -unltest + 2.*nutest*cs2test1*sghtest   
 
            endif

            if ( .not.lsoca_testflow ) then

              df(l1:l2,m,n,jugu:jugu+2)=df(l1:l2,m,n,jugu:jugu+2)+unltest
              df(l1:l2,m,n,jugu+3     )=df(l1:l2,m,n,jugu+3     )+hnltest

            endif
 
            do j=1,3
              unltestm(n,j,jtest)=unltestm(n,j,jtest)+fac*sum(unltest(:,j))
            enddo

            hnltestm(n,jtest)=hnltestm(n,jtest)+fac*sum(hnltest)

            headtt=.false.
 
          enddo

          do j=1,3
            unltestm1(n-n1+1,ipz+1,j,jtest)=unltestm(n,j,jtest)
          enddo

          hnltestm1(n-n1+1,ipz+1,jtest)=hnltestm(n,jtest)

        enddo
      enddo
!
!  do communication for arrays of size nz*nprocz*3*njtest and nz*nprocz*njtest, resp.
!
    if (nprocy>1) then
 
      call mpireduce_sum(unltestm1,unltestm1_tmp,nz*nprocz*3*njtest)  !MR: allreduce?
!!!   call mpibcast_real(unltestm1_tmp,nz*nprocz*3*njtest)

      call mpireduce_sum(hnltestm1,hnltestm1_tmp,nz*nprocz*njtest)  !MR: allreduce?
!!!   call mpibcast_real(hnltestm1_tmp,nz*nprocz*njtest)

      do jtest=0,njtest
        do n=n1,n2
  
          do j=1,3
            unltestm(n,j,jtest)=unltestm1_tmp(n-n1+1,ipz+1,j,jtest)
          enddo

          hnltestm(n,jtest)=hnltestm1_tmp(n-n1+1,ipz+1,jtest)

        enddo
      enddo
    endif
!
!  start with zero
!
    do jtest=0,njtest

      iuxtest=iuutest+4*jtest
      ihhtest=iuxtest+3
 
      do n=n1,n2
        do m=m1,m2

          if ( jtest.eq.0 .or. .not.lsoca_testflow ) then
            
            do j=1,3
              ju = iuxtest+j-1
              df(l1:l2,m,n,ju)=df(l1:l2,m,n,ju)-unltestm(n,j,jtest)
            enddo
            df(l1:l2,m,n,ihhtest)=df(l1:l2,m,n,ihhtest)-hnltestm(n,jtest)

          endif
        enddo
      enddo
    enddo
!
!  reset headtt
!
print*,'unltestm,hnltestm=',unltestm,hnltestm
    headtt=headtt_save
!
    endsubroutine calc_ltestflow_nonlin_terms
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
      select case(jtest)
      case(1); uutest(:,1)=cz(n); uutest(:,2)=0.; uutest(:,3)=0.
      case(2); uutest(:,1)=sz(n); uutest(:,2)=0.; uutest(:,3)=0.
      case(3); uutest(:,1)=0.   ; uutest(:,2)=0.; uutest(:,3)=0.
      case default; uutest(:,:)=0.
      endselect
!
    endsubroutine set_uutest
!***********************************************************************
    subroutine set_uutest_B11_B21 (uutest,jtest)
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
      select case(jtest)
      case(1); uutest(:,1)=cz(n); uutest(:,2)=0.; uutest(:,3)=0.
      case(2); uutest(:,1)=sz(n); uutest(:,2)=0.; uutest(:,3)=0.
      case default; uutest(:,:)=0.
      endselect
!
    endsubroutine set_uutest_B11_B21
!***********************************************************************
    subroutine set_U0test_W11_W22 (U0test,gU0test,H0test,gH0test,jtest)
!
!  set testflow
!
!  23-mar-08/axel: adapted from testflow_z
!
      use Cdata
      use Sub
!
      real, dimension (nx,3) :: U0test,gU0test
      real, dimension (nx) :: H0test, gH0test
      integer :: jtest
!
      intent(in)  :: jtest
      intent(out) :: U0test,gU0test,H0test,gH0test
!
!  set U0test and gU0test for each of the various cases
!
      select case(jtest)
      case(1)
        U0test(:,1)=0.; U0test(:,2)=-wamp*k1sz(n); U0test(:,3)=0.
        gU0test(:,1)=0.; gU0test(:,2)=-wamp*cz(n); gU0test(:,3)=0.
      case(2)
        U0test(:,1)=0.; U0test(:,2)=+wamp*k1cz(n); U0test(:,3)=0.
        gU0test(:,1)=0.; gU0test(:,2)=-wamp*sz(n); gU0test(:,3)=0.
      case(3)
        U0test(:,1)=+wamp*k1sz(n); U0test(:,2)=0.; U0test(:,3)=0.
        gU0test(:,1)=+wamp*cz(n); gU0test(:,2)=0.; gU0test(:,3)=0.
      case(4)
        U0test(:,1)=-wamp*k1cz(n); U0test(:,2)=0.; U0test(:,3)=0.
        gU0test(:,1)=+wamp*sz(n); gU0test(:,2)=0.; gU0test(:,3)=0.
      case default; U0test(:,:)=0.;gU0test=0.
      endselect
      
      H0test=0.; gH0test=0.;
!
    endsubroutine set_U0test_W11_W22
!***********************************************************************
    subroutine rprint_testflow(lreset,lwrite)
!
!  reads and registers print parameters relevant for testflow fields
!
!   3-jun-05/axel: adapted from rprint_magnetic
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
        idiag_ux0mz=0; idiag_uy0mz=0; idiag_uz0mz=0
        idiag_F111z=0; idiag_F211z=0; idiag_F311z=0
        idiag_F121z=0; idiag_F221z=0; idiag_F321z=0
        idiag_F10z=0; idiag_F20z=0; idiag_F30z=0
        idiag_lam11=0; idiag_lam21=0; idiag_lam12=0; idiag_lam22=0
        idiag_nu11=0; idiag_nu21=0; idiag_nu12=0; idiag_nu22=0
        idiag_u11rms=0; idiag_u21rms=0; idiag_u12rms=0; idiag_u22rms=0; idiag_u0rms=0
      endif
!
!  check for those quantities that we want to evaluate online
! 
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'lam11',idiag_lam11)
        call parse_name(iname,cname(iname),cform(iname),'lam21',idiag_lam21)
        call parse_name(iname,cname(iname),cform(iname),'lam12',idiag_lam12)
        call parse_name(iname,cname(iname),cform(iname),'lam22',idiag_lam22)
        call parse_name(iname,cname(iname),cform(iname),'nu11',idiag_nu11)
        call parse_name(iname,cname(iname),cform(iname),'nu21',idiag_nu21)
        call parse_name(iname,cname(iname),cform(iname),'nu12',idiag_nu12)
        call parse_name(iname,cname(iname),cform(iname),'nu22',idiag_nu22)
        call parse_name(iname,cname(iname),cform(iname),'u11rms',idiag_u11rms)
        call parse_name(iname,cname(iname),cform(iname),'u21rms',idiag_u21rms)
        call parse_name(iname,cname(iname),cform(iname),'u12rms',idiag_u12rms)
        call parse_name(iname,cname(iname),cform(iname),'u22rms',idiag_u22rms)
        call parse_name(iname,cname(iname),cform(iname),'u0rms',idiag_u0rms)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ux0mz',idiag_ux0mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uy0mz',idiag_uy0mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uz0mz',idiag_uz0mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F111z',idiag_F111z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F211z',idiag_F211z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F311z',idiag_F311z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F121z',idiag_F121z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F221z',idiag_F221z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F321z',idiag_F321z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F112z',idiag_F112z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F212z',idiag_F212z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F312z',idiag_F312z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F122z',idiag_F122z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F222z',idiag_F222z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F322z',idiag_F322z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F10z',idiag_F10z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F20z',idiag_F20z)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'F30z',idiag_F30z)
      enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
      if (lwr) then
        write(3,*) 'idiag_lam11=',idiag_lam11
        write(3,*) 'idiag_lam21=',idiag_lam21
        write(3,*) 'idiag_lam12=',idiag_lam12
        write(3,*) 'idiag_lam22=',idiag_lam22
        write(3,*) 'idiag_nu11=',idiag_nu11
        write(3,*) 'idiag_nu21=',idiag_nu21
        write(3,*) 'idiag_nu12=',idiag_nu12
        write(3,*) 'idiag_nu22=',idiag_nu22
        write(3,*) 'idiag_u0rms=',idiag_u0rms
        write(3,*) 'idiag_u11rms=',idiag_u11rms
        write(3,*) 'idiag_u21rms=',idiag_u21rms
        write(3,*) 'idiag_u12rms=',idiag_u12rms
        write(3,*) 'idiag_u22rms=',idiag_u22rms
        write(3,*) 'idiag_ux0mz=',idiag_ux0mz
        write(3,*) 'idiag_uy0mz=',idiag_uy0mz
        write(3,*) 'idiag_uz0mz=',idiag_uz0mz
        write(3,*) 'idiag_F111z=',idiag_F111z
        write(3,*) 'idiag_F211z=',idiag_F211z
        write(3,*) 'idiag_F311z=',idiag_F311z
        write(3,*) 'idiag_F121z=',idiag_F121z
        write(3,*) 'idiag_F221z=',idiag_F221z
        write(3,*) 'idiag_F321z=',idiag_F321z
        write(3,*) 'idiag_F112z=',idiag_F112z
        write(3,*) 'idiag_F212z=',idiag_F212z
        write(3,*) 'idiag_F312z=',idiag_F312z
        write(3,*) 'idiag_F122z=',idiag_F122z
        write(3,*) 'idiag_F222z=',idiag_F222z
        write(3,*) 'idiag_F322z=',idiag_F322z
        write(3,*) 'idiag_F10z=',idiag_F10z
        write(3,*) 'idiag_F20z=',idiag_F20z
        write(3,*) 'idiag_F30z=',idiag_F30z
        write(3,*) 'iuutest=',iuutest
!       write(3,*) 'iuxtestpq=',iuxtestpq
!       write(3,*) 'iuztestpq=',iuztestpq
        write(3,*) 'ntestflow=',ntestflow
        write(3,*) 'nnamez=',nnamez
        write(3,*) 'nnamexy=',nnamexy
        write(3,*) 'nnamexz=',nnamexz
      endif
!
    endsubroutine rprint_testflow

endmodule Testflow
