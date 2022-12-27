! $Id$
!
!  This modules solves two reactive scalar advection equations
!  This is used for modeling the spatial evolution of left and
!  right handed aminoacids.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lchiral = .true.
!
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Chiral
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'chiral.h'
!
  integer :: iXX_chiral=0,iYY_chiral=0
  character (len=labellen) :: initXX_chiral='zero',initYY_chiral='zero'
  logical :: llorentzforceEP=.false.
  logical :: linitialize_aa_from_EP=.false.
  logical :: linitialize_aa_from_EP_alpgrad=.false.
  logical :: linitialize_aa_from_EP_betgrad=.false.
  real :: tinitialize_aa_from_EP=0.
  real :: amplXX_chiral=.1, widthXX_chiral=.5
  real :: amplYY_chiral=.1, widthYY_chiral=.5
  real :: kx_XX_chiral=1.,ky_XX_chiral=1.,kz_XX_chiral=1.,radiusXX_chiral=0.
  real :: kx_YY_chiral=1.,ky_YY_chiral=1.,kz_YY_chiral=1.,radiusYY_chiral=0.
  real :: xposXX_chiral=0.,yposXX_chiral=0.,zposXX_chiral=0.
  real :: xposYY_chiral=0.,yposYY_chiral=0.,zposYY_chiral=0.
!
  namelist /chiral_init_pars/ &
       initXX_chiral,amplXX_chiral,kx_XX_chiral,ky_XX_chiral,kz_XX_chiral, &
       initYY_chiral,amplYY_chiral,kx_YY_chiral,ky_YY_chiral,kz_YY_chiral, &
       radiusXX_chiral,widthXX_chiral, &
       radiusYY_chiral,widthYY_chiral, &
       xposXX_chiral,yposXX_chiral,zposXX_chiral, &
       xposYY_chiral,yposYY_chiral,zposYY_chiral
!
  real :: chiral_diffXX=impossible, chiral_diff=0., chiral_crossinhibition=1.,chiral_fidelity=1.
  real :: chiral_fishernu=0., chiral_fisherK=1. !fishers equation growth rate, carrying capac.
  real :: chiral_fisherR=0., chiral_fisherR2=0. !(reinfections, second model corresponds to SIRS)
  real, dimension(3) :: gradX0=(/0.0,0.0,0.0/), gradY0=(/0.0,0.0,0.0/)
  character (len=labellen) :: chiral_reaction='BAHN_model'
  logical :: limposed_gradient=.false.
  logical :: lupw_chiral=.false.
!
  namelist /chiral_run_pars/ &
       chiral_diffXX, chiral_diff, chiral_crossinhibition, chiral_fidelity, &
       chiral_reaction, limposed_gradient, gradX0, gradY0, &
       chiral_fishernu, chiral_fisherK, chiral_fisherR, chiral_fisherR2, &
       lupw_chiral, llorentzforceEP, &
       linitialize_aa_from_EP,tinitialize_aa_from_EP, &
       linitialize_aa_from_EP_alpgrad,linitialize_aa_from_EP_betgrad
!
  integer :: idiag_XX_chiralmax=0, idiag_XX_chiralm=0
  integer :: idiag_YY_chiralmax=0, idiag_YY_chiralm=0
  integer :: idiag_QQm_chiral=0, idiag_QQ21m_chiral=0, idiag_QQ21QQm_chiral=0
  integer :: idiag_brmsEP=0, idiag_bmaxEP=0
  integer :: idiag_jrmsEP=0, idiag_jmaxEP=0
  integer :: idiag_jbmEP=0
!
  contains
!***********************************************************************
    subroutine register_chiral()
!
!  Initialise variables which should know that we solve for passive
!  scalar: iXX_chiral and iYY_chiral; increase nvar accordingly
!
!  28-may-04/axel: adapted from pscalar
!
      use FArrayManager
!
      call farray_register_pde('XX_chiral',iXX_chiral)
      call farray_register_pde('YY_chiral',iYY_chiral)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_chiral
!***********************************************************************
    subroutine initialize_chiral(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  28-may-04/axel: adapted from pscalar
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  default: chiral_diffXX=chiral_diff
!
      if (chiral_diffXX==impossible) chiral_diffXX=chiral_diff
      if (headt) print*,'chiral_diffXX,chiral_diff=',chiral_diffXX,chiral_diff
!
!  set f to zero and then call the same initial condition
!  that was used in start.csh
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_chiral
!***********************************************************************
    subroutine init_chiral(f)
!
!  initialise passive scalar field; called from start.f90
!
!  28-may-04/axel: adapted from pscalar
!
      use Mpicomm
      use Sub
      use Initcond
      use InitialCondition, only: initial_condition_chiral
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  check first for initXX_chiral
!
      select case (initXX_chiral)
        case ('zero'); f(:,:,:,iXX_chiral)=0.
        case ('const'); f(:,:,:,iXX_chiral)=amplXX_chiral
        case ('blob'); call blob(amplXX_chiral,f,iXX_chiral,radiusXX_chiral,xposXX_chiral,yposXX_chiral,zposXX_chiral)
        case ('hat-x'); call hat(amplXX_chiral,f,iXX_chiral,widthXX_chiral,kx=kx_XX_chiral)
        case ('hat-y'); call hat(amplXX_chiral,f,iXX_chiral,widthXX_chiral,ky=ky_XX_chiral)
        case ('hat-z'); call hat(amplXX_chiral,f,iXX_chiral,widthXX_chiral,kz=kz_XX_chiral)
        case ('gaussian-x'); call gaussian(amplXX_chiral,f,iXX_chiral,kx=kx_XX_chiral)
        case ('gaussian-y'); call gaussian(amplXX_chiral,f,iXX_chiral,ky=ky_XX_chiral)
        case ('gaussian-z'); call gaussian(amplXX_chiral,f,iXX_chiral,kz=kz_XX_chiral)
        case ('gaussian-noise'); call gaunoise(amplXX_chiral,f,iXX_chiral)
        case ('positive-noise'); call posnoise(amplXX_chiral,f,iXX_chiral)
        case ('wave-x'); call wave(amplXX_chiral,f,iXX_chiral,kx=kx_XX_chiral)
        case ('wave-y'); call wave(amplXX_chiral,f,iXX_chiral,ky=ky_XX_chiral)
        case ('wave-z'); call wave(amplXX_chiral,f,iXX_chiral,kz=kz_XX_chiral)
        case ('cos(x-cosz)'); call cosxz_cosz(amplXX_chiral,f,iXX_chiral,kx_XX_chiral,kz_XX_chiral)
        case ('cosy_sinz'); call cosy_sinz(amplXX_chiral,f,iXX_chiral,ky_XX_chiral,kz_XX_chiral)
        case ('cosx_cosz'); call cosx_cosz(amplXX_chiral,f,iXX_chiral,kx_XX_chiral,kz_XX_chiral)
        case ('cosx_cosy_cosz'); call cosx_cosy_cosz(amplXX_chiral,f,iXX_chiral,kx_XX_chiral,ky_XX_chiral,kz_XX_chiral)
        case ('cosx_siny_cosz'); call cosx_siny_cosz(amplXX_chiral,f,iXX_chiral,kx_XX_chiral,ky_XX_chiral,kz_XX_chiral)
        case default; call stop_it('init_chiral: bad init_chiral='//trim(initXX_chiral))
      endselect
!
!  check next for initYY_chiral
!
      select case (initYY_chiral)
        case ('zero'); f(:,:,:,iYY_chiral)=0.
        case ('const'); f(:,:,:,iYY_chiral)=amplYY_chiral
        case ('blob'); call blob(amplYY_chiral,f,iYY_chiral,radiusYY_chiral,xposYY_chiral,yposYY_chiral,zposYY_chiral)
        case ('hat-x'); call hat(amplYY_chiral,f,iYY_chiral,widthYY_chiral,kx=kx_YY_chiral)
        case ('hat-y'); call hat(amplYY_chiral,f,iYY_chiral,widthYY_chiral,ky=ky_YY_chiral)
        case ('hat-z'); call hat(amplYY_chiral,f,iYY_chiral,widthYY_chiral,kz=kz_YY_chiral)
        case ('gaussian-x'); call gaussian(amplYY_chiral,f,iYY_chiral,kx=kx_YY_chiral)
        case ('gaussian-y'); call gaussian(amplYY_chiral,f,iYY_chiral,ky=ky_YY_chiral)
        case ('gaussian-z'); call gaussian(amplYY_chiral,f,iYY_chiral,kz=kz_YY_chiral)
        case ('gaussian-noise'); call gaunoise(amplYY_chiral,f,iYY_chiral)
        case ('positive-noise'); call posnoise(amplYY_chiral,f,iYY_chiral)
        case ('wave-x'); call wave(amplYY_chiral,f,iYY_chiral,kx=kx_YY_chiral)
        case ('wave-y'); call wave(amplYY_chiral,f,iYY_chiral,ky=ky_YY_chiral)
        case ('wave-z'); call wave(amplYY_chiral,f,iYY_chiral,kz=kz_YY_chiral)
        case ('cos(y-sinz)'); call cosyz_sinz(amplYY_chiral,f,iYY_chiral,ky_YY_chiral,kz_YY_chiral)
        case ('cosy_sinz'); call cosy_sinz(amplYY_chiral,f,iYY_chiral,ky_YY_chiral,kz_YY_chiral)
        case ('cosy_cosz'); call cosy_cosz(amplYY_chiral,f,iYY_chiral,ky_YY_chiral,kz_YY_chiral)
        case ('cosx_cosy_cosz'); call cosx_cosy_cosz(amplYY_chiral,f,iYY_chiral,kx_YY_chiral,ky_YY_chiral,kz_YY_chiral)
        case ('cosx_siny_cosz'); call cosx_siny_cosz(amplYY_chiral,f,iYY_chiral,kx_YY_chiral,ky_YY_chiral,kz_YY_chiral)
        case ('chiral_list'); call chiral_list(amplYY_chiral,f,iYY_chiral)
        case default; call stop_it('init_chiral: bad init_chiral='//trim(initYY_chiral))
      endselect
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_chiral(f)
!
    endsubroutine init_chiral
!***********************************************************************
    subroutine chiral_list(fact,f,i)
!
!  Read intial conditions from file
!
!  13-apr-20/axel: coded
!
      integer :: nrow, irow, ix, iy, i, l, m, n=1
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ampl, fact
!
      open(1,file='chiral_list.dat')
      read(1,*) nrow
      do irow=1,nrow
        read(1,*) ix,iy,ampl
        if (ix/nx==ipx.and.iy/ny==ipy) then
          l=l1+mod(ix,nx)
          m=m1+mod(iy,ny)
          n=n1
          print*,'AXEL: ',ix,iy,ipx,ipy,l,m
          f(l,m,n,i)=f(l,m,n,i)+fact*ampl
        endif
      enddo
      close(1)
!
    endsubroutine chiral_list
!***********************************************************************
    subroutine pencil_criteria_chiral()
!
!  All pencils that the Chiral module depends on are specified here.
!
!  21-nov-04/anders: coded
!
      lpenc_requested(i_uu)=.true.
      if (llorentzforceEP) lpenc_requested(i_rho1)=.true.
!
    endsubroutine pencil_criteria_chiral
!***********************************************************************
    subroutine pencil_interdep_chiral(lpencil_in)
!
!  Interdependency among pencils provided by the Chiral module
!  is specified here.
!
!  21-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_chiral
!***********************************************************************
    subroutine calc_pencils_chiral(f,p)
!
!  Calculate Chiral pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  21-11-04/anders: coded
!
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_chiral
!***********************************************************************
    subroutine dXY_chiral_dt(f,df,p)
!
!  passive scalar evolution
!  calculate chirality equations in reduced form; see q-bio/0401036
!
!  28-may-04/axel: adapted from pscalar
!   1-jul-09/axel: included gradX, gradY, and allowed for different reactions
!
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx,3,3) :: gXXij_chiral,gYYij_chiral
      real, dimension (nx,3) :: gXX_chiral,gYY_chiral,bbEP,jjEP,jxbEP
      real, dimension (nx) :: bbEP2,jjEP2,jbEP
      real, dimension (nx) :: XX_chiral,ugXX_chiral,del2XX_chiral,dXX_chiral
      real, dimension (nx) :: YY_chiral,ugYY_chiral,del2YY_chiral,dYY_chiral
      real, dimension (nx) :: RRXX_chiral,XX2_chiral
      real, dimension (nx) :: RRYY_chiral,YY2_chiral
      real, dimension (nx) :: RR21_chiral
      real, dimension (nx) :: QQ_chiral,QQ21_chiral,QQ21QQ_chiral
      real, dimension (nx) :: diffus_chiral
      real :: pp,qq,lamchiral
      integer :: j
!
      intent(in)  :: p
      intent(inout) :: f
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE dXY_dt'
      if (headtt) call identify_bcs('XX_chiral',iXX_chiral)
      if (headtt) call identify_bcs('YY_chiral',iYY_chiral)
!
!  gradient of passive scalar
!
      call grad(f,iXX_chiral,gXX_chiral)
      call grad(f,iYY_chiral,gYY_chiral)
!
!  Add diffusion of imposed spatially constant gradient of X or Y.
!  This makes sense mainly for periodic boundary conditions.
!
      if (limposed_gradient) then
        do j=1,3
          gXX_chiral(:,j)=gXX_chiral(:,j)+gradX0(j)
          gYY_chiral(:,j)=gYY_chiral(:,j)+gradY0(j)
        enddo
      endif
!
!  advection term
!
      !call dot_mn(p%uu,gXX_chiral,ugXX_chiral)
      !call dot_mn(p%uu,gYY_chiral,ugYY_chiral)
      call u_dot_grad(f,iXX_chiral,gXX_chiral,p%uu,ugXX_chiral,UPWIND=lupw_chiral)
      call u_dot_grad(f,iYY_chiral,gYY_chiral,p%uu,ugYY_chiral,UPWIND=lupw_chiral)
      df(l1:l2,m,n,iXX_chiral)=df(l1:l2,m,n,iXX_chiral)-ugXX_chiral
      df(l1:l2,m,n,iYY_chiral)=df(l1:l2,m,n,iYY_chiral)-ugYY_chiral
!
!  diffusion term
!
      call del2(f,iXX_chiral,del2XX_chiral)
      call del2(f,iYY_chiral,del2YY_chiral)
      df(l1:l2,m,n,iXX_chiral)=df(l1:l2,m,n,iXX_chiral)+chiral_diffXX*del2XX_chiral
      df(l1:l2,m,n,iYY_chiral)=df(l1:l2,m,n,iYY_chiral)+chiral_diff  *del2YY_chiral
!
!  For Euler Potentials, possibility to add Lorentz force
!
      if (llorentzforceEP) then
        if (lhydro) then
          call cross(gXX_chiral,gYY_chiral,bbEP)
          do j=1,3
            jjEP(:,j)=gXX_chiral(:,j)*del2YY_chiral &
                     -gYY_chiral(:,j)*del2XX_chiral
          enddo
          call cross(jjEP,bbEP,jxbEP)
          if (ldensity) call multsv_mn(p%rho1,jxbEP,jxbEP)
          df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+jxbEP
        endif
      endif
!
!  selection of different reaction terms
!
      select case (chiral_reaction)
!
!  BAHN model
!
      case ('BAHN_model')
!
!  reaction terms
!  X^2/Rtilde^2 - X*R
!  Y^2/Rtilde^2 - Y*R
!
!  for finite crossinhibition (=kI/kS) and finite fidelity (=f) we have
!  R --> RX=X+Y*kI/kS, R --> RY=Y+X*kI/kS, and
!  X2tilde=X^2/2RX, Y2tilde=Y^2/2RY.
!
      XX_chiral=f(l1:l2,m,n,iXX_chiral)
      YY_chiral=f(l1:l2,m,n,iYY_chiral)
      RRXX_chiral=XX_chiral+YY_chiral*chiral_crossinhibition
      RRYY_chiral=YY_chiral+XX_chiral*chiral_crossinhibition
!
!  abbreviations for quadratic quantities
!
      XX2_chiral=.5*XX_chiral**2/max(RRXX_chiral, tini)
      YY2_chiral=.5*YY_chiral**2/max(RRYY_chiral, tini)
      RR21_chiral=1./max(XX2_chiral+YY2_chiral, tini)
!
!  fidelity factor
!
      pp=.5*(1.+chiral_fidelity)
      qq=.5*(1.-chiral_fidelity)
!
!  final reaction equation
!
      dXX_chiral=(pp*XX2_chiral+qq*YY2_chiral)*RR21_chiral-XX_chiral*RRXX_chiral
      dYY_chiral=(pp*YY2_chiral+qq*XX2_chiral)*RR21_chiral-YY_chiral*RRYY_chiral
      df(l1:l2,m,n,iXX_chiral)=df(l1:l2,m,n,iXX_chiral)+dXX_chiral
      df(l1:l2,m,n,iYY_chiral)=df(l1:l2,m,n,iYY_chiral)+dYY_chiral
!
      case('fisher')
      if (headtt) print*,"chiral_reaction='fishers equation'"
      if (headtt) print*,"growth rate=", chiral_fishernu,"carrying capacity=",chiral_fisherK
      XX_chiral=f(l1:l2,m,n,iXX_chiral)
      YY_chiral=f(l1:l2,m,n,iYY_chiral)
      dXX_chiral=(1.-XX_chiral/chiral_fisherK)*XX_chiral*chiral_fishernu
      dYY_chiral=(1.-YY_chiral/chiral_fisherK)*YY_chiral*chiral_fishernu
      df(l1:l2,m,n,iXX_chiral)=df(l1:l2,m,n,iXX_chiral)+dXX_chiral
      df(l1:l2,m,n,iYY_chiral)=df(l1:l2,m,n,iYY_chiral)+dYY_chiral
!
      case('SIR')
      if (headtt) print*,"chiral_reaction='SIR equation'"
      if (headtt) print*,"growth rate=", chiral_fishernu,"carrying capacity=",chiral_fisherK
      XX_chiral=f(l1:l2,m,n,iXX_chiral)
      YY_chiral=f(l1:l2,m,n,iYY_chiral)
      dXX_chiral=-chiral_fishernu*XX_chiral*YY_chiral &
        +chiral_fisherR2*(1.-(XX_chiral+YY_chiral))
      dYY_chiral=+chiral_fishernu*XX_chiral*YY_chiral-chiral_fisherK*YY_chiral &
        +chiral_fisherR*(1.-(XX_chiral+YY_chiral))
      df(l1:l2,m,n,iXX_chiral)=df(l1:l2,m,n,iXX_chiral)+dXX_chiral
      df(l1:l2,m,n,iYY_chiral)=df(l1:l2,m,n,iYY_chiral)+dYY_chiral
!
      case ('nothing')
        if (lroot.and.ip<=5) print*,"chiral_reaction='nothing'"
!
      case default
        write(unit=errormsg,fmt=*) &
             'dXY_chiral_dt: No such value for chiral_reaction: ', &
             trim(chiral_reaction)
        call fatal_error('chiral_reaction',errormsg)
      endselect
!
!  For the timestep calculation, need maximum diffusion
!
      if (lfirst.and.ldt) then
        diffus_chiral=max(chiral_diffXX,chiral_diff)*dxyz_2
        if (headtt.or.ldebug) print*,'dXY_chiral_dt: max(diffus_chiral) =', &
                                      maxval(diffus_chiral)
        maxdiffus=max(maxdiffus,diffus_chiral)
      endif
!
!  diagnostics
!
!  output for double and triple correlators (assume z-gradient of cc)
!  <u_k u_j d_j c> = <u_k c uu.gradXX_chiral>
!
      if (ldiagnos) then
        if (idiag_XX_chiralmax/=0) &
            call max_mn_name(XX_chiral,idiag_XX_chiralmax)
        if (idiag_YY_chiralmax/=0) &
            call max_mn_name(YY_chiral,idiag_YY_chiralmax)
        if (idiag_XX_chiralm/=0) call sum_mn_name(XX_chiral,idiag_XX_chiralm)
        if (idiag_YY_chiralm/=0) call sum_mn_name(YY_chiral,idiag_YY_chiralm)
!
!  extra diagnostics
!
        lamchiral=2.*chiral_fidelity-1.
        QQ_chiral=XX_chiral-YY_chiral
        QQ21_chiral=1.-QQ_chiral**2
        QQ21QQ_chiral=(lamchiral-QQ_chiral**2)/(1.+QQ_chiral**2)*QQ_chiral
        if (idiag_QQm_chiral/=0) call sum_mn_name(QQ_chiral,idiag_QQm_chiral)
        if (idiag_QQ21m_chiral/=0) &
            call sum_mn_name(QQ21_chiral,idiag_QQ21m_chiral)
        if (idiag_QQ21QQm_chiral/=0) &
            call sum_mn_name(QQ21QQ_chiral,idiag_QQ21QQm_chiral)
!
!  Calculate BB = gXX_chiral x gYY_chiral
!
        if (idiag_brmsEP/=0.or.idiag_bmaxEP/=0) then
          call cross(gXX_chiral,gYY_chiral,bbEP)
          call dot2(bbEP,bbEP2)
          if (idiag_brmsEP/=0) &
              call sum_mn_name(bbEP2,idiag_brmsEP,lsqrt=.true.)
          if (idiag_bmaxEP/=0) &
              call max_mn_name(bbEP2,idiag_bmaxEP,lsqrt=.true.)
        endif
!
!  Calculate Ji = Xi*del2Y - Yi*del2X + Xij Yj - Yij Xj
!  and then the max and rms values of that
!
        if (idiag_jrmsEP/=0.or.idiag_jmaxEP/=0.or.idiag_jbmEP/=0) then
          do j=1,3
            jjEP(:,j)=gXX_chiral(:,j)*del2YY_chiral &
                     -gYY_chiral(:,j)*del2XX_chiral
          enddo
          call g2ij(f,iXX_chiral,gXXij_chiral)
          call g2ij(f,iYY_chiral,gYYij_chiral)
          call multmv_mn(+gXXij_chiral,gYY_chiral,jjEP,ladd=.true.)
          call multmv_mn(-gYYij_chiral,gXX_chiral,jjEP,ladd=.true.)
          call dot2(jjEP,jjEP2)
          if (idiag_jrmsEP/=0) &
              call sum_mn_name(jjEP2,idiag_jrmsEP,lsqrt=.true.)
          if (idiag_jmaxEP/=0) &
              call max_mn_name(jjEP2,idiag_jmaxEP,lsqrt=.true.)
!
!  For J.B, we need B, but only if not already calculated.
!
          if (idiag_jbmEP/=0) then
            if (idiag_brmsEP/=0.or.idiag_bmaxEP/=0) then
              call cross(gXX_chiral,gYY_chiral,bbEP)
            endif
            call dot(jjEP,bbEP,jbEP)
            call sum_mn_name(jbEP,idiag_jbmEP)
          endif
        endif
      endif
!
    endsubroutine dXY_chiral_dt
!***********************************************************************
    subroutine chiral_before_boundary(f)
!
!  initialize aa from Euler potentials (EP), but only once when t=0
!
!   4-jul-09/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: gXX_chiral,gYY_chiral
      integer :: j
!
      intent(inout) :: f
!
!  Possibility to initialize magnetic vector potential in terms
!  of Euler potentials X and Y.
!
      if (linitialize_aa_from_EP) then
        if (t>=tinitialize_aa_from_EP) then
          linitialize_aa_from_EP=.false.
          if (lmagnetic) then
            if (lroot) print*,'initialize_aa_from_EP at t=',t
            do n=n1,n2
            do m=m1,m2
              call grad(f,iXX_chiral,gXX_chiral)
              call grad(f,iYY_chiral,gYY_chiral)
              if (limposed_gradient) then
                do j=1,3
                  gXX_chiral(:,j)=gXX_chiral(:,j)+gradX0(j)
                  gYY_chiral(:,j)=gYY_chiral(:,j)+gradY0(j)
                enddo
              endif
              do j=1,3
                if (linitialize_aa_from_EP_alpgrad) then
                  f(l1:l2,m,n,j+iaa-1)=f(l1:l2,m,n,iXX_chiral)*gYY_chiral(:,j)
                elseif (linitialize_aa_from_EP_betgrad) then
                  f(l1:l2,m,n,j+iaa-1)=-f(l1:l2,m,n,iYY_chiral)*gXX_chiral(:,j)
                else
                  f(l1:l2,m,n,j+iaa-1)=.5*( &
                    f(l1:l2,m,n,iXX_chiral)*gYY_chiral(:,j) &
                   -f(l1:l2,m,n,iYY_chiral)*gXX_chiral(:,j))
                endif
              enddo
            enddo
            enddo
          endif
        endif
      endif
!
    endsubroutine chiral_before_boundary
!***********************************************************************
    subroutine read_chiral_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=chiral_init_pars, IOSTAT=iostat)
!
    endsubroutine read_chiral_init_pars
!***********************************************************************
    subroutine write_chiral_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=chiral_init_pars)
!
    endsubroutine write_chiral_init_pars
!***********************************************************************
    subroutine read_chiral_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=chiral_run_pars, IOSTAT=iostat)
!
    endsubroutine read_chiral_run_pars
!***********************************************************************
    subroutine write_chiral_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=chiral_run_pars)
!
    endsubroutine write_chiral_run_pars
!***********************************************************************
    subroutine rprint_chiral(lreset,lwrite)
!
!  reads and registers print parameters relevant for chirality
!
!  28-may-04/axel: adapted from pscalar
!
      use Diagnostics
!
      integer :: iname
      logical :: lreset
      logical, optional :: lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_XX_chiralmax=0; idiag_XX_chiralm=0
        idiag_YY_chiralmax=0; idiag_YY_chiralm=0
        idiag_QQm_chiral=0; idiag_QQ21m_chiral=0; idiag_QQ21QQm_chiral=0
        idiag_brmsEP=0; idiag_bmaxEP=0
        idiag_jrmsEP=0; idiag_jmaxEP=0
        idiag_jbmEP=0
      endif
!
!  check for those quantities that we want to evaluate online
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),&
            'XXm',idiag_XX_chiralm)
        call parse_name(iname,cname(iname),cform(iname),&
            'YYm',idiag_YY_chiralm)
        call parse_name(iname,cname(iname),cform(iname),&
            'XXmax',idiag_XX_chiralmax)
        call parse_name(iname,cname(iname),cform(iname),&
            'YYmax',idiag_YY_chiralmax)
        call parse_name(iname,cname(iname),cform(iname),&
            'QQm',idiag_QQm_chiral)
        call parse_name(iname,cname(iname),cform(iname),&
            'QQ21m',idiag_QQ21m_chiral)
        call parse_name(iname,cname(iname),cform(iname),&
            'QQ21QQm',idiag_QQ21QQm_chiral)
        call parse_name(iname,cname(iname),cform(iname),&
            'brmsEP',idiag_brmsEP)
        call parse_name(iname,cname(iname),cform(iname),&
            'bmaxEP',idiag_bmaxEP)
        call parse_name(iname,cname(iname),cform(iname),&
            'jrmsEP',idiag_jrmsEP)
        call parse_name(iname,cname(iname),cform(iname),&
            'jmaxEP',idiag_jmaxEP)
        call parse_name(iname,cname(iname),cform(iname),&
            'jbmEP',idiag_jbmEP)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='XX_chiral'.or.cnamev=='YY_chiral'.or.cnamev=='DQ_chiral') &
          cformv='DEFINED'
      endif
!
    endsubroutine rprint_chiral
!***********************************************************************
    subroutine get_slices_chiral(f,slices)
!
!  Write slices for animation of Chiral variables.
!
!  26-jul-06/tony: coded
!
      use Slices_methods, only:assign_slices_scal

      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Chirality fields: XX
!
        case ('XX_chiral'); call assign_slices_scal(slices,f,iXX_chiral)
!
!  Chirality fields: YY
!
        case ('YY_chiral'); call assign_slices_scal(slices,f,iYY_chiral)
!
!  Chirality fields: DQ
!
        case ('DQ_chiral')
          if (lwrite_slice_yz ) &
            call QQ_chiral(f(ix_loc,m1:m2,n1:n2,iXX_chiral)-f(ix_loc,m1:m2,n1:n2,iYY_chiral),slices%yz)
          if (lwrite_slice_xz ) &
            call QQ_chiral(f(l1:l2,iy_loc,n1:n2,iXX_chiral)-f(l1:l2,iy_loc,n1:n2,iYY_chiral),slices%xz)
          if (lwrite_slice_xy ) &
            call QQ_chiral(f(l1:l2,m1:m2,iz_loc,iXX_chiral)-f(l1:l2,m1:m2,iz_loc,iYY_chiral),slices%xy)
          if (lwrite_slice_xy2) &
            call QQ_chiral(f(l1:l2,m1:m2,iz2_loc,iXX_chiral)-f(l1:l2,m1:m2,iz2_loc,iYY_chiral),slices%xy2)
          if (lwrite_slice_xy3) &
            call QQ_chiral(f(l1:l2,m1:m2,iz3_loc,iXX_chiral)-f(l1:l2,m1:m2,iz3_loc,iYY_chiral),slices%xy3)
          if (lwrite_slice_xy4) &
            call QQ_chiral(f(l1:l2,m1:m2,iz4_loc,iXX_chiral)-f(l1:l2,m1:m2,iz4_loc,iYY_chiral),slices%xy4)
          if (lwrite_slice_xz2) &
            call QQ_chiral(f(l1:l2,iy2_loc,n1:n2,iXX_chiral)-f(l1:l2,iy2_loc,n1:n2,iYY_chiral),slices%xz2)
          slices%ready=.true.
!
      endselect
!
contains
      subroutine QQ_chiral(source,dest)
!
        real, dimension(:,:) :: source,dest
!
        dest=source**2 
        dest=source*(1.-dest)/(1.+dest)
!
      endsubroutine QQ_chiral

    endsubroutine get_slices_chiral
!***********************************************************************
endmodule Chiral
