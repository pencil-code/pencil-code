! $Id$

!  This modules solves two reactive scalar advection equations
!  This is used for modeling the spatial evolution of left and
!  right handed aminoacids.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Chiral

  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include 'chiral.h'

  integer :: iXX_chiral=0,iYY_chiral=0

  character (len=labellen) :: initXX_chiral='zero',initYY_chiral='zero'

  ! input parameters
  real :: amplXX_chiral=.1, widthXX_chiral=.5
  real :: amplYY_chiral=.1, widthYY_chiral=.5
  real :: kx_XX_chiral=1.,ky_XX_chiral=1.,kz_XX_chiral=1.,radiusXX_chiral=0.
  real :: kx_YY_chiral=1.,ky_YY_chiral=1.,kz_YY_chiral=1.,radiusYY_chiral=0.
  real :: xposXX_chiral=0.,yposXX_chiral=0.,zposXX_chiral=0.
  real :: xposYY_chiral=0.,yposYY_chiral=0.,zposYY_chiral=0.

  namelist /chiral_init_pars/ &
       initXX_chiral,amplXX_chiral,kx_XX_chiral,ky_XX_chiral,kz_XX_chiral, &
       initYY_chiral,amplYY_chiral,kx_YY_chiral,ky_YY_chiral,kz_YY_chiral, &
       radiusXX_chiral,widthXX_chiral, &
       radiusYY_chiral,widthYY_chiral, &
       xposXX_chiral,yposXX_chiral,zposXX_chiral, &
       xposYY_chiral,yposYY_chiral,zposYY_chiral

  ! run parameters
  real :: chiral_diff=0., chiral_crossinhibition=1.,chiral_fidelity=1.

  namelist /chiral_run_pars/ &
       chiral_diff,chiral_crossinhibition,chiral_fidelity

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_XX_chiralmax=0, idiag_XX_chiralm=0
  integer :: idiag_YY_chiralmax=0, idiag_YY_chiralm=0
  integer :: idiag_QQm_chiral=0, idiag_QQ21m_chiral=0, idiag_QQ21QQm_chiral=0

  contains

!***********************************************************************
    subroutine register_chiral()
!
!  Initialise variables which should know that we solve for passive
!  scalar: iXX_chiral and iYY_chiral; increase nvar accordingly
!
!  28-may-04/axel: adapted from pscalar
!
      use Cdata
      use FArrayManager
!
      lchiral = .true.
!
      call farray_register_pde('XX_chiral',iXX_chiral)
      call farray_register_pde('YY_chiral',iYY_chiral)
!
!  Identify version number.
!
      if (lroot) call cvs_id( &
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
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if (NO_WARN) print*,'f=',f
    endsubroutine initialize_chiral
!***********************************************************************
    subroutine init_chiral(f)
!
!  initialise passive scalar field; called from start.f90
!
!  28-may-04/axel: adapted from pscalar
!
      use Cdata
      use Mpicomm
      use Sub
      use Initcond
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  check first for initXX_chiral
!
      select case(initXX_chiral)
        case('zero'); f(:,:,:,iXX_chiral)=0.
        case('const'); f(:,:,:,iXX_chiral)=amplXX_chiral
        case('blob'); call blob(amplXX_chiral,f,iXX_chiral,radiusXX_chiral,xposXX_chiral,yposXX_chiral,zposXX_chiral)
        case('hat-x'); call hat(amplXX_chiral,f,iXX_chiral,widthXX_chiral,kx=kx_XX_chiral)
        case('hat-y'); call hat(amplXX_chiral,f,iXX_chiral,widthXX_chiral,ky=ky_XX_chiral)
        case('hat-z'); call hat(amplXX_chiral,f,iXX_chiral,widthXX_chiral,kz=kz_XX_chiral)
        case('gaussian-x'); call gaussian(amplXX_chiral,f,iXX_chiral,kx=kx_XX_chiral)
        case('gaussian-y'); call gaussian(amplXX_chiral,f,iXX_chiral,ky=ky_XX_chiral)
        case('gaussian-z'); call gaussian(amplXX_chiral,f,iXX_chiral,kz=kz_XX_chiral)
        case('positive-noise'); call posnoise(amplXX_chiral,f,iXX_chiral)
        case('wave-x'); call wave(amplXX_chiral,f,iXX_chiral,kx=kx_XX_chiral)
        case('wave-y'); call wave(amplXX_chiral,f,iXX_chiral,ky=ky_XX_chiral)
        case('wave-z'); call wave(amplXX_chiral,f,iXX_chiral,kz=kz_XX_chiral)
        case('cosx_cosy_cosz'); call cosx_cosy_cosz(amplXX_chiral,f,iXX_chiral,kx_XX_chiral,ky_XX_chiral,kz_XX_chiral)
        case default; call stop_it('init_chiral: bad init_chiral='//trim(initXX_chiral))
      endselect
!
!  check next for initYY_chiral
!
      select case(initYY_chiral)
        case('zero'); f(:,:,:,iYY_chiral)=0.
        case('const'); f(:,:,:,iYY_chiral)=amplYY_chiral
        case('blob'); call blob(amplYY_chiral,f,iYY_chiral,radiusYY_chiral,xposYY_chiral,yposYY_chiral,zposYY_chiral)
        case('hat-x'); call hat(amplYY_chiral,f,iYY_chiral,widthYY_chiral,kx=kx_YY_chiral)
        case('hat-y'); call hat(amplYY_chiral,f,iYY_chiral,widthYY_chiral,ky=ky_YY_chiral)
        case('hat-z'); call hat(amplYY_chiral,f,iYY_chiral,widthYY_chiral,kz=kz_YY_chiral)
        case('gaussian-x'); call gaussian(amplYY_chiral,f,iYY_chiral,kx=kx_YY_chiral)
        case('gaussian-y'); call gaussian(amplYY_chiral,f,iYY_chiral,ky=ky_YY_chiral)
        case('gaussian-z'); call gaussian(amplYY_chiral,f,iYY_chiral,kz=kz_YY_chiral)
        case('positive-noise'); call posnoise(amplYY_chiral,f,iYY_chiral)
        case('wave-x'); call wave(amplYY_chiral,f,iYY_chiral,kx=kx_YY_chiral)
        case('wave-y'); call wave(amplYY_chiral,f,iYY_chiral,ky=ky_YY_chiral)
        case('wave-z'); call wave(amplYY_chiral,f,iYY_chiral,kz=kz_YY_chiral)
        case('cosx_cosy_cosz'); call cosx_cosy_cosz(amplYY_chiral,f,iYY_chiral,kx_YY_chiral,ky_YY_chiral,kz_YY_chiral)
        case default; call stop_it('init_chiral: bad init_chiral='//trim(initYY_chiral))
      endselect
!
    endsubroutine init_chiral
!***********************************************************************
    subroutine pencil_criteria_chiral()
!
!  All pencils that the Chiral module depends on are specified here.
!
!  21-11-04/anders: coded
!
      use Cdata
!
      lpenc_requested(i_uu)=.true.
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
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      real, dimension (nx,3) :: gXX_chiral,gYY_chiral
      real, dimension (nx) :: XX_chiral,ugXX_chiral,del2XX_chiral,dXX_chiral
      real, dimension (nx) :: YY_chiral,ugYY_chiral,del2YY_chiral,dYY_chiral
      real, dimension (nx) :: RRXX_chiral,XX2_chiral
      real, dimension (nx) :: RRYY_chiral,YY2_chiral
      real, dimension (nx) :: RR21_chiral
      real, dimension (nx) :: QQ_chiral,QQ21_chiral,QQ21QQ_chiral
      real :: pp,qq,lamchiral
      integer :: j
!
      intent(in)  :: f,p
      intent(out) :: df
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
      call dot_mn(p%uu,gXX_chiral,ugXX_chiral)
      call dot_mn(p%uu,gYY_chiral,ugYY_chiral)
!
!  advection term
!
      if (lhydro) df(l1:l2,m,n,iXX_chiral)=df(l1:l2,m,n,iXX_chiral)-ugXX_chiral
      if (lhydro) df(l1:l2,m,n,iYY_chiral)=df(l1:l2,m,n,iYY_chiral)-ugYY_chiral
!
!  diffusion term
!
      call del2(f,iXX_chiral,del2XX_chiral)
      call del2(f,iYY_chiral,del2YY_chiral)
      df(l1:l2,m,n,iXX_chiral)=df(l1:l2,m,n,iXX_chiral)+chiral_diff*del2XX_chiral
      df(l1:l2,m,n,iYY_chiral)=df(l1:l2,m,n,iYY_chiral)+chiral_diff*del2YY_chiral
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
!  For the timestep calculation, need maximum diffusion
!
      if (lfirst.and.ldt) diffus_chiral=diffus_chiral+chiral_diff*dxyz_2
      if (headtt.or.ldebug) print*,'dXY_chiral_dt: max(diffus_chiral) =', &
                                    maxval(diffus_chiral)
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
      endif
!
    endsubroutine dXY_chiral_dt
!***********************************************************************
    subroutine read_chiral_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=chiral_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=chiral_init_pars,ERR=99)
      endif


99    return
    endsubroutine read_chiral_init_pars
!***********************************************************************
    subroutine write_chiral_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=chiral_init_pars)

    endsubroutine write_chiral_init_pars
!***********************************************************************
    subroutine read_chiral_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=chiral_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=chiral_run_pars,ERR=99)
      endif


99    return
    endsubroutine read_chiral_run_pars
!***********************************************************************
    subroutine write_chiral_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=chiral_run_pars)

    endsubroutine write_chiral_run_pars
!***********************************************************************
    subroutine rprint_chiral(lreset,lwrite)
!
!  reads and registers print parameters relevant for chirality
!
!  28-may-04/axel: adapted from pscalar
!
      use Sub
!
      integer :: iname,inamez
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_XX_chiralmax=0; idiag_XX_chiralm=0
        idiag_YY_chiralmax=0; idiag_YY_chiralm=0
        idiag_QQm_chiral=0; idiag_QQ21m_chiral=0; idiag_QQ21QQm_chiral=0
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
      enddo
!
!  write column where which chiral variable is stored
!
      if (lwr) then
        write(3,*) 'i_XX_chiralm=',idiag_XX_chiralm
        write(3,*) 'i_YY_chiralm=',idiag_YY_chiralm
        write(3,*) 'i_XX_chiralmax=',idiag_XX_chiralmax
        write(3,*) 'i_YY_chiralmax=',idiag_YY_chiralmax
        write(3,*) 'i_QQm_chiral=',idiag_QQm_chiral
        write(3,*) 'i_QQ21m_chiral=',idiag_QQ21m_chiral
        write(3,*) 'i_QQ21QQm_chiral=',idiag_QQ21QQm_chiral
        write(3,*) 'iXX_chiral=',iXX_chiral
        write(3,*) 'iYY_chiral=',iYY_chiral
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
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
      real, dimension (nx,ny) :: QQ_chiral_xy
      real, dimension (nx,ny) :: QQ_chiral_xy2
      real, dimension (nx,nz) :: QQ_chiral_xz
      real, dimension (ny,nz) :: QQ_chiral_yz
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Chirality fields: XX
!
        case ('XX_chiral')
          slices%yz=f(ix_loc,m1:m2,n1:n2,iXX_chiral)
          slices%xz=f(l1:l2,iy_loc,n1:n2,iXX_chiral)
          slices%xy=f(l1:l2,m1:m2,iz_loc,iXX_chiral)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,iXX_chiral)
          slices%ready = .true.
!
!  Chirality fields: YY
!
        case ('YY_chiral')
          slices%yz=f(ix_loc,m1:m2,n1:n2,iYY_chiral)
          slices%xz=f(l1:l2,iy_loc,n1:n2,iYY_chiral)
          slices%xy=f(l1:l2,m1:m2,iz_loc,iYY_chiral)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,iYY_chiral)
          slices%ready = .true.
!
!  Chirality fields: DQ
!
        case ('DQ_chiral')
          QQ_chiral_yz=f(ix_loc,m1:m2,n1:n2,iXX_chiral)-f(ix_loc,m1:m2,n1:n2,iYY_chiral)
          QQ_chiral_xz=f(l1:l2,iy_loc,n1:n2,iXX_chiral)-f(l1:l2,iy_loc,n1:n2,iYY_chiral)
          QQ_chiral_xy=f(l1:l2,m1:m2,iz_loc,iXX_chiral)-f(l1:l2,m1:m2,iz_loc,iYY_chiral)
          QQ_chiral_xy2=f(l1:l2,m1:m2,iz2_loc,iXX_chiral)-f(l1:l2,m1:m2,iz2_loc,iYY_chiral)
          slices%yz=QQ_chiral_yz*(1.-QQ_chiral_yz**2)/(1.+QQ_chiral_yz**2)
          slices%xz=QQ_chiral_xz*(1.-QQ_chiral_xz**2)/(1.+QQ_chiral_xz**2)
          slices%xy=QQ_chiral_xy*(1.-QQ_chiral_xy**2)/(1.+QQ_chiral_xy**2)
          slices%xy2=&
              QQ_chiral_xy2*(1.-QQ_chiral_xy2**2)/(1.+QQ_chiral_xy2**2)
          slices%ready = .true.
!
      endselect
!
    endsubroutine get_slices_chiral
!***********************************************************************
endmodule Chiral

