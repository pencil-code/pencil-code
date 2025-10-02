! $Id$
!
!  This module takes care of everything related to velocity
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lneutralvelocity = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED divun; visc_heatn; un2; unij(3,3); uun(3); snij(3,3); snij2
! PENCILS PROVIDED ungun(3); del2un(3); del6un(3); graddivun(3); advec_uun; advec_csn2
!
!***************************************************************
module NeutralVelocity
!
  use Cdata
  use Messages
!
  implicit none
!
  include 'neutralvelocity.h'
!
!  Init parameters.
!
  real :: ampl_unx=0.0, ampl_uny=0.0, ampl_unz=0.0
  real :: kx_uun=1., ky_uun=1., kz_uun=1.
  real, dimension (ninit) :: ampluun=0.0
  character (len=labellen), dimension(ninit) :: inituun='nothing'
  character (len=labellen) :: borderuun='nothing'
  real, dimension(3) :: uun_const=(/0.,0.,0./)
  logical :: lcoriolis_force=.true., lcentrifugal_force=.false.
  logical :: ladvection_velocity=.true.,lpressuregradient=.true.
  logical :: lviscneutral=.true.,lelectron_pressure=.false.
  real :: colldrag=0,electron_pressure=1
  real :: nun=0.,csn0=0.,csn20,nun_hyper3=0.
  real :: rnoise_int=impossible,rnoise_ext=impossible
  real :: uun_right=0., uun_left=0., widthuun=.1
  real, dimension (nx,3,3) :: unij5
  character (len=labellen),dimension(ninit) :: iviscn=''
!
  namelist /neutralvelocity_init_pars/ &
      ampluun, ampl_unx, ampl_uny, ampl_unz, &
      inituun, uun_const, Omega, lcoriolis_force, lcentrifugal_force, &
      ladvection_velocity,colldrag,csn0,kx_uun,ky_uun,kz_uun,&
      uun_right, uun_left, widthuun, &
      rnoise_int,rnoise_ext
!
!  Run parameters.
!
  logical :: lupw_uun=.false.
  logical :: lfreeze_unint=.false.,lfreeze_unext=.false.
!
  namelist /neutralvelocity_run_pars/ &
       Omega,theta, lupw_uun, &
       borderuun, lfreeze_unint, lpressuregradient, &
       lfreeze_unext,lcoriolis_force,lcentrifugal_force,ladvection_velocity, &
       colldrag,nun,lviscneutral,iviscn,nun,csn0,nun_hyper3, &
       lelectron_pressure,electron_pressure
!
!  Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_un2m=0,idiag_unm2=0
  integer :: idiag_unxpt=0,idiag_unypt=0,idiag_unzpt=0,idiag_dtcn=0
  integer :: idiag_dtun=0,idiag_unrms=0,idiag_unmax=0,idiag_unzrms=0
  integer :: idiag_unzrmaxs=0
  integer :: idiag_unxmax=0,idiag_unymax=0,idiag_unzmax=0
  integer :: idiag_unxm=0,idiag_unym=0,idiag_unzm=0
  integer :: idiag_unx2m=0,idiag_uny2m=0,idiag_unz2m=0
  integer :: idiag_unx2mz=0,idiag_uny2mz=0,idiag_unz2mz=0
  integer :: idiag_unx2my=0,idiag_uny2my=0,idiag_unz2my=0
  integer :: idiag_unx2mx=0,idiag_uny2mx=0,idiag_unz2mx=0
  integer :: idiag_unxunym=0,idiag_unxunzm=0,idiag_unyunzm=0
  integer :: idiag_unxunymz=0,idiag_unxunzmz=0,idiag_unyunzmz=0,idiag_rnunxunymz=0
  integer :: idiag_unxunymy=0,idiag_unxunzmy=0,idiag_unyunzmy=0
  integer :: idiag_unxunymx=0,idiag_unxunzmx=0,idiag_unyunzmx=0
  integer :: idiag_unxmz=0,idiag_unymz=0,idiag_unzmz=0,idiag_unmx=0,idiag_unmy=0
  integer :: idiag_unxmy=0,idiag_unymy=0,idiag_unzmy=0,idiag_un2mz=0
  integer :: idiag_unmz=0,idiag_unxmxy=0,idiag_unymxy=0,idiag_unzmxy=0
  integer :: idiag_unxmx=0,idiag_unymx=0,idiag_unzmx=0
  integer :: idiag_unrmphi=0,idiag_unpmphi=0,idiag_unzmphi=0,idiag_un2mphi=0
  integer :: idiag_neutralangmom=0
  integer :: idiag_un2mr=0,idiag_unrunpmr=0
  integer :: idiag_unrmr=0,idiag_unpmr=0,idiag_unzmr=0
  integer :: idiag_divunm=0, idiag_pndivunm=0, idiag_dtnun=0
  integer :: idiag_fricneut=0, idiag_fricions=0
  integer :: idiag_epsKn=0       ! DIAG_DOC: $\left<2\nu_n\varrho_n\Strain_n^2\right>$
!
! Auxiliaries
!
  real, dimension(nx) :: cions_rhon,cneut_rho,diffus_nun
  integer, dimension(ninit) :: enum_iviscn = 0
  integer :: enum_borderuun = 0

  contains
!***********************************************************************
    subroutine register_neutralvelocity()
!
!  Initialise variables which should know that we solve the neutralvelocity
!  equations: iuun, etc; increase nvar accordingly.
!
!  28-feb-07/wlad: adapted
!
      use FArrayManager, only: farray_register_pde
!
      if (.not.lcartesian_coords) call not_implemented('register_neutralvelocity','non-Cartesian')
!
!  Indices to access uun.
!
      call farray_register_pde('uun',iuun,vector=3)
      iunx = iuun; iuny = iuun+1; iunz = iuun+2
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id$")
!
!  Writing files for use with IDL.
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',uun $'
          if (nvar == mvar) write(4,*) ',uun'
        else
          write(4,*) ',uun $'
        endif
        write(15,*) 'uun = fltarr(mx,my,mz,3)*one'
      endif
!
    endsubroutine register_neutralvelocity
!***********************************************************************
    subroutine initialize_neutralvelocity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  28-feb-07/wlad: adapted
!
      use BorderProfiles, only: request_border_driving
!
! Check any module dependencies
!
      if (.not.leos) call fatal_error('initialize_neutralvelocity', &
                                      'EOS=noeos but an EQUATION OF STATE for the fluid is required')
!
!  set freezing arrays
!
      if (lfreeze_unint) lfreeze_varint(iunx:iunz) = .true.
      if (lfreeze_unext) lfreeze_varext(iunx:iunz) = .true.
!
!  csn20
!
      csn20=csn0**2
!
!  Turn off advection for 0-D runs.
!
      if (dimensionality==0) then
        ladvection_velocity=.false.
        if (lroot) print*, 'initialize_neutralvelocity: 0-D run, turned off advection of velocity'
      endif
!
      if (ladvection_velocity.and.lupw_uun) &
        call warning('initialize_neutralvelocity','upwinding advection term not well tested')

      select case (borderuun)
      case ('zero','0','constant')
!
!  Tell the BorderProfiles module if we intend to use border driving, so
!  that the modules can request the right pencils.
!
        call request_border_driving(borderuun)
      case ('initial-condition')
        call fatal_error("initialize_neutralvelocity","borderuu = 'initial-condition': set mcount/ncount")
      case ('nothing')
        if (lroot.and.ip<=5) print*,"initialize_neutralvelocity: borderuu='nothing'"
      case default
        call fatal_error('initialize_neutralvelocity','no such borderuun: '//trim(borderuun))
      endselect
!
!  Turn off neutral viscosity if zero viscosity
!
      if ((nun == 0.).and.(nun_hyper3==0.)) lviscneutral=.false.
      if (lroot) print*,'initialize_neutralvelocity: lviscneutral,nun,nun_hyper3=', &
                        lviscneutral,nun,nun_hyper3
!
      endsubroutine initialize_neutralvelocity
!***********************************************************************
    subroutine read_neutralvelocity_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=neutralvelocity_init_pars, IOSTAT=iostat)
!
    endsubroutine read_neutralvelocity_init_pars
!***********************************************************************
    subroutine write_neutralvelocity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=neutralvelocity_init_pars)
!
    endsubroutine write_neutralvelocity_init_pars
!***********************************************************************
    subroutine read_neutralvelocity_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=neutralvelocity_run_pars, IOSTAT=iostat)
!
    endsubroutine read_neutralvelocity_run_pars
!***********************************************************************
    subroutine write_neutralvelocity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=neutralvelocity_run_pars)
!
    endsubroutine write_neutralvelocity_run_pars
!***********************************************************************
    subroutine init_uun(f)
!
!  initialise uun and lnrhon; called from start.f90
!
!  28-feb-07/wlad: adapted
!
      use Initcond
      use InitialCondition, only: initial_condition_uun
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: prof
      integer :: j,i
!
!  inituun corresponds to different initializations of uun (called from start).
!
      do j=1,ninit

        select case (inituun(j))

        case ('nothing'); if (lroot .and. j==1) print*,'init_uun: nothing'
        case ('zero', '0')
          if (lroot) print*,'init_uu: zero velocity'
          ! Ensure really is zero, as may have used lread_oldsnap
          f(:,:,:,iunx:iunz)=0.
        case ('const_uun'); do i=1,3; f(:,:,:,iuun+i-1) = uun_const(i); enddo
        case ('gaussian-noise'); call gaunoise(ampluun(j),f,iunx,iunz)
        case ('gaussian-noise-x'); call gaunoise(ampluun(j),f,iunx)
        case ('gaussian-noise-y'); call gaunoise(ampluun(j),f,iuny)
        case ('gaussian-noise-z'); call gaunoise(ampluun(j),f,iunz)
        case ('gaussian-noise-xy'); call gaunoise(ampluun(j),f,iunx,iuny)
        case ('sinwave-x'); call sinwave(ampluun(j),f,iunx,kx=kx_uun)
        case ('sinwave-y'); call sinwave(ampluun(j),f,iuny,ky=ky_uun)
        case ('sinwave-z'); call sinwave(ampluun(j),f,iunz,kz=kz_uun)
        case ('gaussian-noise-rprof')
          call gaunoise_rprof(ampluun(j),f,iunx,iunz,rnoise_int,rnoise_ext)
        case ('follow-ions'); f(:,:,:,iunx:iunz)=f(:,:,:,iux:iuz)
!
        case ('shock-tube')
!
!  shock tube test (should be consistent with density module)
!
          if (lroot) print*,'init_uun: polytopic standing shock'
          do n=n1,n2; do m=m1,m2
            prof=.5*(1.+tanh(x(l1:l2)/widthuun))
            f(l1:l2,m,n,iunx)=uun_left+(uun_right-uun_left)*prof
          enddo; enddo
!
        case ('tanhx')
!
!  Burgers shock
!
          if (lroot) print*,'init_uun: Burgers shock'
          prof=-ampluun(j)*tanh(.5*x(l1:l2)/widthuun)
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,iunx)=prof
          enddo; enddo
!
        case default
          call fatal_error("init_uun",'No such value for inituu: '//trim(inituun(j)))

        endselect
!
!  End loop over initial conditions
!
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_uun(f)
!
    endsubroutine init_uun
!***********************************************************************
    subroutine pencil_criteria_neutralvelocity()
!
!  All pencils that the Neutralvelocity module depends on are specified here.
!
!  28-feb-07/wlad: adapted
!
      if (ladvection_velocity) lpenc_requested(i_ungun)=.true.
      if (ldt) lpenc_requested(i_uun)=.true.
!
      lpenc_requested(i_uu)   =.true.
      lpenc_requested(i_rho)  =.true.
      lpenc_requested(i_rhon) =.true.
      lpenc_requested(i_rho1) =.true.
      lpenc_requested(i_rhon1)=.true.
      lpenc_requested(i_alpha)=.true.
      lpenc_requested(i_zeta) =.true.
!
      if (any(iviscn=='nun-const')) then
        !lpenc_requested(i_snij)=.true.
        !lpenc_requested(i_glnrhon)=.true.
        lpenc_requested(i_del2un)=.true.
        lpenc_requested(i_snglnrhon)=.true.
        lpenc_requested(i_graddivun)=.true.
      endif
      if (any(iviscn=='hyper3_nun-const')) then
        lpenc_requested(i_del6un)=.true.
        lpenc_requested(i_glnrhon)=.true.
      endif
      if (any(iviscn=='rhon_nun-const')) then
        lpenc_requested(i_del2un)=.true.
        lpenc_requested(i_graddivun)=.true.
      endif
      !if ( lneutraldensity.and.                                    &
      !     any((iviscn=='shock').or.(iviscn=='nun-shock'))) then
      !  lpenc_requested(i_graddivun)=.true.
      !  lpenc_requested(i_shock)=.true.
      !  lpenc_requested(i_gshock)=.true.
      !  lpenc_requested(i_divun)=.true.
      !  lpenc_requested(i_glnrhon)=.true.
      !endif
!
      if (lpressuregradient) lpenc_requested(i_glnrhon)=.true.

      if (lhydro.and.lelectron_pressure) lpenc_requested(i_fpres)=.true.
!
!  diagnostic pencils
!
      lpenc_diagnos(i_uun)=.true.
      if (idiag_un2mphi/=0) lpenc_diagnos2d(i_un2)=.true.
      if (idiag_unrms/=0 .or. idiag_unmax/=0 .or.  &
          idiag_un2m/=0 .or. idiag_unm2/=0 .or. idiag_un2mz/=0) &
        lpenc_diagnos(i_un2)=.true.
      if (idiag_unrmr/=0 .or. idiag_unrunpmr/=0) then
        lpenc_diagnos(i_pomx)=.true.
        lpenc_diagnos(i_pomy)=.true.
      endif
      if (idiag_divunm/=0 .or. idiag_pndivunm/=0) lpenc_diagnos(i_divun)=.true.
!
      if (idiag_unpmr/=0 .or. idiag_unrunpmr/=0) then
        lpenc_diagnos(i_phix)=.true.
        lpenc_diagnos(i_phiy)=.true.
      endif
!
      if (idiag_epsKn/=0) then
        lpenc_diagnos(i_rhon)=.true.
        lpenc_diagnos(i_snij2)=.true.
        lpenc_diagnos(i_visc_heatn)=.true.
      endif
!
    endsubroutine pencil_criteria_neutralvelocity
!***********************************************************************
    subroutine pencil_interdep_neutralvelocity(lpencil_in)
!
!  Interdependency among pencils from the Neutralvelocity module
!   is specified here.
!
!  28-feb-07/wlad: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_un2)) lpencil_in(i_uun)=.true.
      if (lpencil_in(i_divun)) lpencil_in(i_unij)=.true.
      if (lpencil_in(i_snij)) then
        lpencil_in(i_unij)=.true.
        lpencil_in(i_divun)=.true.
      endif
      if (lpencil_in(i_ungun)) then
        lpencil_in(i_uun)=.true.
        lpencil_in(i_unij)=.true.
      endif
!
      if (lpencil_in(i_snij)) then
        if (any(iviscn=='nun-const')) then
          lpencil_in(i_unij)=.true.
          lpencil_in(i_divun)=.true.
        endif
      endif
!
      if (lpencil_in(i_snij2)) lpencil_in(i_snij)=.true.
!
    endsubroutine pencil_interdep_neutralvelocity
!***********************************************************************
    subroutine calc_pencils_neutralvelocity(f,p)
!
!  Calculate Neutralvelocity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  28-feb-07/wlad: adapted
!
      use General, only: notanumber
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      integer :: i, j
!
      intent(in) :: f
      intent(inout) :: p
! uun
      if (lpencil(i_uun)) p%uun=f(l1:l2,m,n,iunx:iunz)
! un2
      if (lpencil(i_un2)) call dot2_mn(p%uun,p%un2)
! unij
      if (lpencil(i_unij)) call gij(f,iuun,p%unij,1)
! divun
      if (lpencil(i_divun)) call div_mn(p%unij,p%divun,p%uun)
! snij
      if (lpencil(i_snij)) then
        do j=1,3
          do i=1,3
            p%snij(:,i,j)=.5*(p%unij(:,i,j)+p%unij(:,j,i))
          enddo
          p%snij(:,j,j)=p%snij(:,j,j)-onethird*p%divun
        enddo
      endif
! snij2
      if (lpencil(i_snij2)) call multm2_sym_mn(p%snij,p%snij2)
! ungun
      if (lpencil(i_ungun)) &
        call u_dot_grad(f,iuun,p%unij,p%uun,p%ungun,UPWIND=lupw_uun)
! del6un
      if (lpencil(i_del6un)) call del6v(f,iuun,p%del6un)
! del2un
! graddivun
      if (lpencil(i_del2un)) then
        if (lpencil(i_graddivun)) then
          call del2v_etc(f,iuun,DEL2=p%del2un,GRADDIV=p%graddivun)
        else
          call del2v(f,iuun,p%del2un)
        endif
      else
        if (lpencil(i_graddivun)) call del2v_etc(f,iuun,GRADDIV=p%graddivun)
      endif

!
      if (lupdate_courant_dt.and.dimensionality>0) then
!
!  ``uun/dx'' for timestep
!
        p%advec_uun=sum(abs(p%uun)*dline_1,2)
        if (notanumber(p%advec_uun)) print*, 'p%advec_uun  =',p%advec_uun
!
! csn2/dx^2 for timestep
! have to include a selection of equation of state...
!
        p%advec_csn2=csn20*dxyz_2
        if (notanumber(p%advec_csn2)) print*, 'advec_csn2 =',p%advec_csn2
!
        if (headtt.or.ldebug) then
          print*,'calc_pencils_neutralvelocity: max(advec_csn2) =',maxval(p%advec_csn2)
          print*,'calc_pencils_neutralvelocity: max(advec_uun) =',maxval(p%advec_uun)
        endif
      endif
!
    endsubroutine calc_pencils_neutralvelocity
!***********************************************************************
    subroutine duun_dt(f,df,p)
!
!  velocity evolution
!  calculate dun/dt = - un.gradun - 2Omega x un + Fpres + grav + Fvisc
!
!  28-feb-07/wlad: adapted
!
      use Diagnostics, only: max_mn_name
      use General, only: notanumber
      use Sub, only: identify_bcs
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)    :: f
      intent(inout) :: df,p
!
      real, dimension (nx) :: ionization,recombination
!
      real :: c2,s2
      integer :: j
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'duun_dt: SOLVE'
      if (headtt) then
        call identify_bcs('unx',iunx)
        call identify_bcs('uny',iuny)
        call identify_bcs('unz',iunz)
      endif
!
!  advection term, -u.gradu
!
      if (ladvection_velocity) df(l1:l2,m,n,iunx:iunz) = df(l1:l2,m,n,iunx:iunz) - p%ungun
!
!  Coriolis force, -2*Omega x u (unless lprecession=T)
!  Omega=(-sin_theta, 0, cos_theta)
!  theta corresponds to latitude, but to have the box located on the
!  right hand side of the sphere (grav still pointing dowward and then
!  Omega to the left), one should choose Omega=-90 for the equator,
!  for example.
!
      if (Omega/=0.) then
        if (theta==0) then
          if (lcoriolis_force) then
            if (headtt) print*,'duun_dt: add Coriolis force; Omega=',Omega
            c2=2*Omega
            df(l1:l2,m,n,iunx)=df(l1:l2,m,n,iunx)+c2*p%uun(:,2)
            df(l1:l2,m,n,iuny)=df(l1:l2,m,n,iuny)-c2*p%uun(:,1)
!
!  add centrifugal force (doing this with periodic boundary
!  conditions in x and y would not be compatible, so it is
!  therefore usually ignored in those cases!)
!
          endif
          if (lcentrifugal_force) then
            if (headtt) print*,'duun_dt: add Centrifugal force; Omega=',Omega
            df(l1:l2,m,n,iunx)=df(l1:l2,m,n,iunx)+x(l1:l2)*Omega**2
            df(l1:l2,m,n,iuny)=df(l1:l2,m,n,iuny)+y(  m  )*Omega**2
          endif
        else
!
!  add Coriolis force with an angle (defined such that theta=-60,
!  for example, would correspond to 30 degrees latitude).
!  Omega=(sin(theta), 0, cos(theta)).
!
          if (lcoriolis_force) then
            if (headtt) &
                print*,'duun_dt: Coriolis force; Omega, theta=', Omega, theta
            c2=2*Omega*cos(theta*pi/180.)
            s2=2*Omega*sin(theta*pi/180.)
            df(l1:l2,m,n,iunx)=df(l1:l2,m,n,iunx)+c2*p%uun(:,2)
            df(l1:l2,m,n,iuny)=df(l1:l2,m,n,iuny)-c2*p%uun(:,1)+s2*p%uun(:,3)
            df(l1:l2,m,n,iunz)=df(l1:l2,m,n,iunz)              -s2*p%uun(:,2)
          endif
        endif
      endif
!
!  Neutral-ion collision, ionization and recombination
!  Remember that cneut*p%rho enters below, so another p%rho is factored out.
!
      ionization=p%zeta*p%rho1
      recombination=p%alpha*p%rho*p%rhon1
      cions_rhon=(colldrag+ionization)*p%rhon
      cneut_rho=(colldrag+recombination)*p%rho
!
      do j=1,3
!
! neutrals gain momentum by recombination
!
        df(l1:l2,m,n,iuun-1+j)=df(l1:l2,m,n,iuun-1+j) + cneut_rho*(p%uu(:,j)-p%uun(:,j))
!
! ions gain momentum by ionization and electron pressure
!
        if (lhydro) then
          df(l1:l2,m,n,iuu-1+j)=df(l1:l2,m,n,iuu-1+j) - cions_rhon*(p%uu(:,j)-p%uun(:,j))
!
! add electron pressure to the ions if needed
! This adds to the already entered contribution from noentropy.f90
! df(l1:l2,m,n,iux:iuz)=df(l1:l2,m,n,iux:iuz)+p%fpres
! and thus implies altogether a factor of 2, which is correct.
!
          if (lelectron_pressure) df(l1:l2,m,n,iuu-1+j)=df(l1:l2,m,n,iuu-1+j)+electron_pressure*p%fpres(:,j)
        endif
!
      enddo

      if (dimensionality>0) then
!
! calculate viscous force on neutrals
!
        if (lviscneutral) call calc_viscous_force_neutral(f,df,p)
!
! add pressure gradient on neutrals
!
        if (lpressuregradient) df(l1:l2,m,n,iunx:iunz)=df(l1:l2,m,n,iunx:iunz)-csn20*p%glnrhon
!
      endif
!
      if (lupdate_courant_dt.and.dimensionality>0) then
        advec2=advec2+p%advec_csn2
        maxadvec=maxadvec+p%advec_uun
      endif
!
!  Apply border profiles
!
      if (lborder_profiles) call set_border_neutralvelocity(f,df,p)
!
      call calc_diagnostics_neutralvel(p)

    endsubroutine duun_dt
!***********************************************************************
    subroutine calc_diagnostics_neutralvel(p)
!
!  Calculate maxima and rms values for diagnostic purposes
!
      use Diagnostics
      use Sub, only: dot_mn

      type (pencil_case) :: p

      real, dimension (nx) :: udelu_neut, udelu_ions

      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'duun_dt: Calculate maxima and rms values...'
        call sum_mn_name(p%un2,idiag_unrms,lsqrt=.true.)
        call max_mn_name(p%un2,idiag_unmax,lsqrt=.true.)
        if (idiag_unzrms/=0) call sum_mn_name(p%uun(:,3)**2,idiag_unzrms,lsqrt=.true.)
        if (idiag_unzrmaxs/=0) call max_mn_name(p%uun(:,3)**2,idiag_unzrmaxs,lsqrt=.true.)
        call max_mn_name(p%uun(:,1),idiag_unxmax)
        call max_mn_name(p%uun(:,2),idiag_unymax)
        call max_mn_name(p%uun(:,3),idiag_unzmax)
        call sum_mn_name(p%un2,idiag_un2m)
        call sum_mn_name(p%divun,idiag_divunm)
        if (idiag_pndivunm/=0) call sum_mn_name(csn20*p%rhon*p%divun,idiag_pndivunm)
        if (idiag_fricneut/=0) then
          call dot_mn(p%uu-p%uun,p%uun,udelu_neut)
          call sum_mn_name(cneut_rho*p%rhon*udelu_neut,idiag_fricneut)
        endif
        if (idiag_fricions/=0) then
          call dot_mn(p%uu-p%uun,p%uu,udelu_ions)
          call sum_mn_name(-cions_rhon*p%rho*udelu_ions,idiag_fricions)
        endif
        call max_mn_name(p%un2,idiag_unm2)
        call sum_mn_name(p%uun(:,1),idiag_unxm)
        call sum_mn_name(p%uun(:,2),idiag_unym)
        call sum_mn_name(p%uun(:,3),idiag_unzm)
        if (idiag_unx2m/=0)    call sum_mn_name(p%uun(:,1)**2,idiag_unx2m)
        if (idiag_uny2m/=0)    call sum_mn_name(p%uun(:,2)**2,idiag_uny2m)
        if (idiag_unz2m/=0)    call sum_mn_name(p%uun(:,3)**2,idiag_unz2m)
        if (idiag_unxunym/=0)   call sum_mn_name(p%uun(:,1)*p%uun(:,2),idiag_unxunym)
        if (idiag_unxunzm/=0)   call sum_mn_name(p%uun(:,1)*p%uun(:,3),idiag_unxunzm)
        if (idiag_unyunzm/=0)   call sum_mn_name(p%uun(:,2)*p%uun(:,3),idiag_unyunzm)
        !if (idiag_dunxdzma/=0) call sum_mn_name(abs(p%uij(:,1,3)),idiag_dunxdzma)
        !if (idiag_dunydzma/=0) call sum_mn_name(abs(p%uij(:,2,3)),idiag_dunydzma)
        !if (idiag_neutralangmom/=0) &
        !     call sum_lim_mn_name(p%rhon*(p%uun(:,2)*x(l1:l2)-p%uun(:,1)*y(m)),&
        !     idiag_neutralangmom,p)
        if (idiag_epsKn/=0) call sum_mn_name(p%visc_heatn*p%rhon,idiag_epsKn)
!
!  kinetic field components at one point (=pt)
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          call save_name(p%uun(lpoint-nghost,1),idiag_unxpt)
          call save_name(p%uun(lpoint-nghost,2),idiag_unypt)
          call save_name(p%uun(lpoint-nghost,3),idiag_unzpt)
        endif
!
        if (ldt) then
          if (lviscneutral) then
            if (idiag_dtnun/=0) call max_mn_name(diffus_nun/cdtv,idiag_dtnun,l_dt=.true.)
          endif

          if (dimensionality>0) then
            if (idiag_dtun/=0) call max_mn_name(p%advec_uun/cdt,idiag_dtun,l_dt=.true.)
            if (idiag_dtcn/=0) call max_mn_name(sqrt(p%advec_csn2)/cdt,idiag_dtcn,l_dt=.true.)
          endif
        endif
!
!  this doesn't need to be as frequent (check later)
!
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1davgfirst) then
        call xysum_mn_name_z(p%uun(:,1),idiag_unxmz)
        call xysum_mn_name_z(p%uun(:,2),idiag_unymz)
        call xysum_mn_name_z(p%uun(:,3),idiag_unzmz)
        call xzsum_mn_name_y(p%uun(:,1),idiag_unxmy)
        call xzsum_mn_name_y(p%uun(:,2),idiag_unymy)
        call xzsum_mn_name_y(p%uun(:,3),idiag_unzmy)
        call yzsum_mn_name_x(p%uun(:,1),idiag_unxmx)
        call yzsum_mn_name_x(p%uun(:,2),idiag_unymx)
        call yzsum_mn_name_x(p%uun(:,3),idiag_unzmx)
        if (idiag_unx2mz/=0) call xysum_mn_name_z(p%uun(:,1)**2,idiag_unx2mz)
        if (idiag_uny2mz/=0) call xysum_mn_name_z(p%uun(:,2)**2,idiag_uny2mz)
        if (idiag_unz2mz/=0) call xysum_mn_name_z(p%uun(:,3)**2,idiag_unz2mz)
        if (idiag_unx2my/=0) call xzsum_mn_name_y(p%uun(:,1)**2,idiag_unx2my)
        if (idiag_uny2my/=0) call xzsum_mn_name_y(p%uun(:,2)**2,idiag_uny2my)
        if (idiag_unz2my/=0) call xzsum_mn_name_y(p%uun(:,3)**2,idiag_unz2my)
        if (idiag_unx2mx/=0) call yzsum_mn_name_x(p%uun(:,1)**2,idiag_unx2mx)
        if (idiag_uny2mx/=0) call yzsum_mn_name_x(p%uun(:,2)**2,idiag_uny2mx)
        if (idiag_unz2mx/=0) call yzsum_mn_name_x(p%uun(:,3)**2,idiag_unz2mx)
        if (idiag_unxunymz/=0) call xysum_mn_name_z(p%uun(:,1)*p%uun(:,2),idiag_unxunymz)
        if (idiag_unxunzmz/=0) call xysum_mn_name_z(p%uun(:,1)*p%uun(:,3),idiag_unxunzmz)
        if (idiag_unyunzmz/=0) call xysum_mn_name_z(p%uun(:,2)*p%uun(:,3),idiag_unyunzmz)
        if (idiag_rnunxunymz/=0) call xysum_mn_name_z(p%rhon*p%uun(:,1)*p%uun(:,2),idiag_rnunxunymz)
        if (idiag_unxunymy/=0) call xzsum_mn_name_y(p%uun(:,1)*p%uun(:,2),idiag_unxunymy)
        if (idiag_unxunzmy/=0) call xzsum_mn_name_y(p%uun(:,1)*p%uun(:,3),idiag_unxunzmy)
        if (idiag_unyunzmy/=0) call xzsum_mn_name_y(p%uun(:,2)*p%uun(:,3),idiag_unyunzmy)
        if (idiag_unxunymx/=0) call yzsum_mn_name_x(p%uun(:,1)*p%uun(:,2),idiag_unxunymx)
        if (idiag_unxunzmx/=0) call yzsum_mn_name_x(p%uun(:,1)*p%uun(:,3),idiag_unxunzmx)
        if (idiag_unyunzmx/=0) call yzsum_mn_name_x(p%uun(:,2)*p%uun(:,3),idiag_unyunzmx)
!  phi-z averages
        call phizsum_mn_name_r(p%un2,idiag_un2mr)
        if (idiag_unrmr/=0) call phizsum_mn_name_r(p%uun(:,1)*p%pomx+p%uun(:,2)*p%pomy,idiag_unrmr)
        if (idiag_unpmr/=0) call phizsum_mn_name_r(p%uun(:,1)*p%phix+p%uun(:,2)*p%phiy,idiag_unpmr)
        call phizsum_mn_name_r(p%uun(:,3),idiag_unzmr)
      endif
!
!  2-D averages.
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        if (idiag_unrmphi/=0) call phisum_mn_name_rz(p%uun(:,1)*p%pomx+p%uun(:,2)*p%pomy,idiag_unrmphi)
        if (idiag_unpmphi/=0) call phisum_mn_name_rz(p%uun(:,1)*p%phix+p%uun(:,2)*p%phiy,idiag_unpmphi)
        if (idiag_unzmphi/=0) call phisum_mn_name_rz(p%uun(:,3),idiag_unzmphi)
        call phisum_mn_name_rz(p%un2,idiag_un2mphi)
        call zsum_mn_name_xy(p%uun(:,1),idiag_unxmxy)
        call zsum_mn_name_xy(p%uun(:,2),idiag_unymxy)
        call zsum_mn_name_xy(p%uun(:,3),idiag_unzmxy)
        call zsum_mn_name_xy(p%un2,idiag_un2mz)
      endif
!
    endsubroutine calc_diagnostics_neutralvel
!***********************************************************************
    subroutine set_border_neutralvelocity(f,df,p)
!
!  Calculates the driving term for the border profile
!  of the uu variable.
!
!  28-jul-06/wlad: coded
!
      use BorderProfiles, only: border_driving
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      real, dimension(mx,my,mz,mvar) :: df

      real, dimension(nx,3) :: f_target
      integer :: ju,j
!
! MR: the following comment apparently belongs to somewhere else:
! these tmps and where's are needed because these square roots
! go negative in the frozen inner disc if the sound speed is big enough
! (like a corona, no neutralvelocitystatic equilibrium)
!
      select case (borderuun)
      case ('zero','0')
        f_target=0.
        do j=1,3
          ju=j+iuun-1
          call border_driving(f,df,p,f_target(:,j),ju)
        enddo
      case ('initial-condition')
        !f_target=f(l1:l2,mcount,ncount,iunx:iunz)
        do j=1,3
          ju=j+iuun-1
          call border_driving(f,df,p,f_target(:,j),ju)
        enddo
      case ('constant')
        do j=1,3
          f_target(:,j) = uun_const(j)
        enddo
        do j=1,3
          ju=j+iuun-1
          call border_driving(f,df,p,f_target(:,j),ju)
        enddo
      case ('nothing')
      endselect
!
!
    endsubroutine set_border_neutralvelocity
!***********************************************************************
    subroutine calc_viscous_force_neutral(f,df,p)
!
!  calculate viscous force term for right hand side of momentum equation
!
!  28-feb-07/wlad: coded
!
      use Deriv, only: der6
      use Diagnostics
      use General, only: itoa
      use Sub, only: multmv, gij
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
      intent(inout) :: df

      real, dimension(nx,3) :: fvisc,unij5glnrhon
      real, dimension(nx) :: munrhon1,tmp,diffus_nun3
      real, dimension (nx,3,3) :: unij5
      integer :: i,j,jj,ju

!
      fvisc=0.
      diffus_nun=0.
      diffus_nun3=0.
!
!  Initialize p%visc_heatn to zero.
!
      if (lpencil(i_visc_heatn)) p%visc_heatn=0.
!
!  Loop over different species.
!
      do j=1,ninit
         select case (iviscn(j))
!
         case ('rhon_nun-const')
!
!  viscous force: mu/rho*(del2u+graddivu/3)
!  -- the correct expression for rho*nu=const
!
            if (headtt) print*,'Viscous force (neutral):  mun/rhon*(del2un+graddivun/3)'
            munrhon1=nun*p%rhon1  !(=mun/rhon)
            do i=1,3
               fvisc(:,i)=fvisc(:,i)+munrhon1*(p%del2un(:,i)+onethird*p%graddivun(:,i))
            enddo
            if (lpencil(i_visc_heatn)) p%visc_heatn=p%visc_heatn + 2*nun*p%snij2*p%rhon1
            if (lupdate_courant_dt) diffus_nun=diffus_nun+munrhon1*dxyz_2

         case ('nun-const')
!
!  viscous force: nu*(del2u+graddivu/3+2S.glnrho)
!  -- the correct expression for nu=const
!
            if (headtt) &
                 print*,'Viscous force (neutral): nun*(del2un+graddivun/3+2Sn.glnrhon)'
            if (lneutraldensity) then
               fvisc=fvisc+2*nun*p%snglnrhon+nun*(p%del2un+onethird*p%graddivun)
            else
               fvisc=fvisc+nun*(p%del2un+onethird*p%graddivun)
            endif
!
            if (lpencil(i_visc_heatn)) p%visc_heatn=p%visc_heatn + 2*nun*p%snij2
            if (lupdate_courant_dt) diffus_nun=diffus_nun+nun*dxyz_2
!
         case ('hyper3_nun-const')
!
!  Viscous force: nun*(del6un+Sn.glnrhon), where Sn_ij=d^5 un_i/dx_j^5
!
            call gij(f,iuun,unij5,5)
            call multmv(unij5,p%glnrhon,unij5glnrhon)
!
            if (headtt) print*, 'Viscous force (neutral): nun*(del6un+Sn.glnrhon)'
            fvisc = fvisc + nun_hyper3*(p%del6un+unij5glnrhon)
            if (lupdate_courant_dt) diffus_nun3=diffus_nun3+nun_hyper3*dxyz_6
!
         case ('hyper3-cyl','hyper3_cyl','hyper3-sph','hyper3_sph')
!
!  Viscous force: anysotropic hyperviscosity
!
           do jj=1,3
             ju=jj+iuun-1
             do i=1,3
               call der6(f,ju,tmp,i,IGNOREDX=.true.)
               fvisc(:,jj)=fvisc(:,jj)+nun_hyper3*pi4_1*tmp*dline_1(:,i)**2
             enddo
             if (lupdate_courant_dt) diffus_nun3=diffus_nun3+nun_hyper3*pi4_1*dxmin_pencil**4
           enddo
!
     !    case ('shock','nun-shock')
     !      if (.not. lshock) &
     !           call stop_it('calc_viscous_force_neutral: shock viscosity'// &
     !           ' but module setting SHOCK=noshock')
     !      if (nun_shock==0.) &
     !           call fatal_error('calc_viscous_force_neutral:', &
     !           'Viscosity coefficient nun_shock is zero!')
     !      if (lneutraldensity) then
     !     !tobi: The following only works with operator overloading for pencils
     !     !      (see sub.f90). Commented out because it seems to be slower.
     !     !tmp=nu_shock*(p%shock*(p%divu*p%glnrho+p%graddivu)+p%divu*p%gshock)
     !        call multsv(p%divun,p%glnrhon,tmp2)
     !        tmp=tmp2 + p%graddivun
     !        call multsv(nun_shock*p%shock,tmp,tmp2)
     !        call multsv_add(tmp2,nun_shock*p%divun,p%gshockn,tmp)
     !        fvisc=fvisc+tmp
     !        if (lupdate_courant_dt) diffus_total=diffus_total+(nu_shock*p%shock)
     !      endif
            !
         case ('')
            ! do nothing
         case default
            call fatal_error('calc_viscous_force_neutral','no such iviscn('// &
                             trim(itoa(i))//'): '//trim(iviscn(i)))
         endselect
      enddo
!
      if (lupdate_courant_dt) then
        maxdiffus=max(maxdiffus,diffus_nun)
        maxdiffus3=max(maxdiffus3,diffus_nun3)
      endif
!
! Add viscosity to the equation of motion
!
     df(l1:l2,m,n,iunx:iunz) = df(l1:l2,m,n,iunx:iunz) + fvisc
!
    endsubroutine calc_viscous_force_neutral
!***********************************************************************
    subroutine rprint_neutralvelocity(lreset,lwrite)
!
!  reads and registers print parameters relevant for neutralvelocity part
!
!  28-feb-07/wlad: adapted
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
!
      integer :: iname,inamez,inamey,inamex,ixy,irz,inamer
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
        idiag_un2m=0; idiag_unm2=0
        idiag_unxpt=0; idiag_unypt=0; idiag_unzpt=0; idiag_dtun=0
        idiag_dtnun=0; idiag_dtcn=0
        idiag_unrms=0; idiag_unmax=0; idiag_unzrms=0; idiag_unzrmaxs=0
        idiag_unxmax=0; idiag_unymax=0; idiag_unzmax=0
        idiag_unxm=0; idiag_unym=0; idiag_unzm=0
        idiag_unx2m=0; idiag_uny2m=0; idiag_unz2m=0
        idiag_unxunym=0; idiag_unxunzm=0; idiag_unyunzm=0
        idiag_unxunymz=0; idiag_unxunzmz=0; idiag_unyunzmz=0; idiag_unxunymz=0
        idiag_unmx=0; idiag_unmy=0; idiag_unmz=0
        idiag_unrmphi=0; idiag_unpmphi=0; idiag_unzmphi=0; idiag_un2mphi=0
        idiag_unxmy=0; idiag_unymy=0; idiag_unzmy=0
        idiag_unx2my=0; idiag_uny2my=0; idiag_unz2my=0
        idiag_unxunymy=0; idiag_unxunzmy=0; idiag_unyunzmy=0
        idiag_neutralangmom=0;
        idiag_unrunpmr=0; idiag_divunm=0; idiag_pndivunm=0
        idiag_fricneut=0; idiag_fricions=0
        idiag_un2mr=0; idiag_unrmr=0; idiag_unpmr=0; idiag_unzmr=0
        idiag_epsKn=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_neutralvelocity: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'un2m',idiag_un2m)
        call parse_name(iname,cname(iname),cform(iname),'unm2',idiag_unm2)
        call parse_name(iname,cname(iname),cform(iname),'dtun',idiag_dtun)
        call parse_name(iname,cname(iname),cform(iname),'dtcn',idiag_dtcn)
        call parse_name(iname,cname(iname),cform(iname),'dtnun',idiag_dtnun)
        call parse_name(iname,cname(iname),cform(iname),'divunm',idiag_divunm)
        call parse_name(iname,cname(iname),cform(iname),'pndivunm',idiag_pndivunm)
        call parse_name(iname,cname(iname),cform(iname),'fricneut',idiag_fricneut)
        call parse_name(iname,cname(iname),cform(iname),'fricions',idiag_fricions)
        call parse_name(iname,cname(iname),cform(iname),'unrms',idiag_unrms)
        call parse_name(iname,cname(iname),cform(iname),'unmax',idiag_unmax)
        call parse_name(iname,cname(iname),cform(iname),'unxmax',idiag_unxmax)
        call parse_name(iname,cname(iname),cform(iname),'unymax',idiag_unymax)
        call parse_name(iname,cname(iname),cform(iname),'unzmax',idiag_unzmax)
        call parse_name(iname,cname(iname),cform(iname),'unzrms',idiag_unzrms)
        call parse_name(iname,cname(iname),cform(iname),'unzrmaxs',idiag_unzrmaxs)
        call parse_name(iname,cname(iname),cform(iname),'unxm',idiag_unxm)
        call parse_name(iname,cname(iname),cform(iname),'unym',idiag_unym)
        call parse_name(iname,cname(iname),cform(iname),'unzm',idiag_unzm)
        call parse_name(iname,cname(iname),cform(iname),'unx2m',idiag_unx2m)
        call parse_name(iname,cname(iname),cform(iname),'uny2m',idiag_uny2m)
        call parse_name(iname,cname(iname),cform(iname),'unz2m',idiag_unz2m)
        call parse_name(iname,cname(iname),cform(iname),'unxunym',idiag_unxunym)
        call parse_name(iname,cname(iname),cform(iname),'unxunzm',idiag_unxunzm)
        call parse_name(iname,cname(iname),cform(iname),'unyunzm',idiag_unyunzm)
        call parse_name(iname,cname(iname),cform(iname),'unmx',idiag_unmx)
        call parse_name(iname,cname(iname),cform(iname),'unmy',idiag_unmy)
        call parse_name(iname,cname(iname),cform(iname),'unmz',idiag_unmz)
        call parse_name(iname,cname(iname),cform(iname),'unxpt',idiag_unxpt)
        call parse_name(iname,cname(iname),cform(iname),'unypt',idiag_unypt)
        call parse_name(iname,cname(iname),cform(iname),'unzpt',idiag_unzpt)
        call parse_name(iname,cname(iname),cform(iname),'neutralangmom',idiag_neutralangmom)
        call parse_name(iname,cname(iname),cform(iname),'epsKn',idiag_epsKn)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'unxmz',idiag_unxmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'unymz',idiag_unymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'unzmz',idiag_unzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'unx2mz',idiag_unx2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'uny2mz',idiag_uny2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'unz2mz',idiag_unz2mz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'unxunymz',idiag_unxunymz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'unxunzmz',idiag_unxunzmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'unyunzmz',idiag_unyunzmz)
      enddo
!
!  check for those quantities for which we want xz-averages
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'unxmy',idiag_unxmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'unymy',idiag_unymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'unzmy',idiag_unzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'unx2my',idiag_unx2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'uny2my',idiag_uny2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'unz2my',idiag_unz2my)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'unxunymy',idiag_unxunymy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'unxunzmy',idiag_unxunzmy)
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'unyunzmy',idiag_unyunzmy)
      enddo
!
!  check for those quantities for which we want yz-averages
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'unxmx',idiag_unxmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'unymx',idiag_unymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'unzmx',idiag_unzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'unx2mx',idiag_unx2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'uny2mx',idiag_uny2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'unz2mx',idiag_unz2mx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'unxunymx',idiag_unxunymx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'unxunzmx',idiag_unxunzmx)
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'unyunzmx',idiag_unyunzmx)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do ixy=1,nnamexy
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'unxmxy',idiag_unxmxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'unymxy',idiag_unymxy)
        call parse_name(ixy,cnamexy(ixy),cformxy(ixy),'unzmxy',idiag_unzmxy)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'unrmphi',idiag_unrmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'unpmphi',idiag_unpmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'unzmphi',idiag_unzmphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'un2mphi',idiag_un2mphi)
      enddo
!
!  check for those quantities for which we want phiz-averages
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'unrmr',  idiag_unrmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'unpmr',  idiag_unpmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'unzmr',  idiag_unzmr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'un2mr',  idiag_un2mr)
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'unrunpmr',idiag_unrunpmr)
      enddo
!
!  write column where which neutralvelocity variable is stored
!
      if (lwr) then
        call farray_index_append('iuun',iuun)
        call farray_index_append('iunx',iunx)
        call farray_index_append('iuny',iuny)
        call farray_index_append('iunz',iunz)
      endif
!
    endsubroutine rprint_neutralvelocity
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=20
    integer(KIND=ikind8), dimension(n_pars) :: p_par
    integer :: i

    call copy_addr(lcoriolis_force,p_par(1)) ! bool
    call copy_addr(lcentrifugal_force,p_par(2)) ! bool
    call copy_addr(ladvection_velocity,p_par(3)) ! bool
    call copy_addr(lpressuregradient,p_par(4)) ! bool
    call copy_addr(lviscneutral,p_par(5)) ! bool
    call copy_addr(lelectron_pressure,p_par(6)) ! bool
    call copy_addr(colldrag,p_par(7))
    call copy_addr(electron_pressure,p_par(8))
    call copy_addr(nun,p_par(9))
    call copy_addr(csn20,p_par(10))
    call copy_addr(nun_hyper3,p_par(11))
    call copy_addr(lupw_uun,p_par(12)) ! bool
    call copy_addr(idiag_fricneut,p_par(13)) ! int
    call copy_addr(idiag_fricions,p_par(14)) ! int
    call copy_addr(uun_const,p_par(15)) ! real3
    do i = 1,ninit; call string_to_enum(enum_iviscn(i),iviscn(i)); enddo
    call copy_addr(enum_iviscn,p_par(16)) ! int (ninit)
    call string_to_enum(enum_borderuun,borderuun)
    call copy_addr(enum_borderuun,p_par(17)) ! int

    endsubroutine pushpars2c
!***********************************************************************
endmodule Neutralvelocity
