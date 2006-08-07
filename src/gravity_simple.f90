! $Id: gravity_simple.f90,v 1.10 2006-08-07 10:45:03 ajohan Exp $

!
!  This module takes care of simple types of gravity, i.e. where
!    gx=gx(x) or gy=gy(y) or gz=gz(z)
!  Here the gravity master pencils gravx_xpencil, gravy_ypencil and
!  gravz_zpencil only need to be calculated once, and then these can 
!  simply be added to the equations of motion again and again.
!

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED gg
!
!***************************************************************

module Gravity

  use Cdata
  use Cparam
  use Messages

  implicit none

  include 'gravity.h'

  interface potential
    module procedure potential_global
    module procedure potential_penc
    module procedure potential_point
  endinterface
!
!  Parameters used throughout entire module
!
  real, dimension(nx) :: gravx_xpencil=0.,potx_xpencil=0.
  real, dimension(ny) :: gravy_ypencil=0.,poty_ypencil=0.
  real, dimension(nz) :: gravz_zpencil=0.,potz_zpencil=0.
  real :: gravx=0.,gravy=0.,gravz=0.
  real :: kx_gg=1.,ky_gg=1.,kz_gg=1.,grav_const=1.,reduced_top=1.
  real :: xgrav=impossible,ygrav=impossible,zgrav=impossible
  real :: xinfty=0.,yinfty=0.,zinfty=0.
  real :: dgravx=0.,pot_ratio=1.
  real :: z1=0.,z2=1.,zref=0.
  real :: nu_epicycle=1.,nu_epicycle2=1.
  real :: r0_pot=0.    ! peak radius for smoothed potential
  real :: g_A, g_C
  real, parameter :: g_A_cgs=4.4e-9, g_C_cgs=1.7e-9
  double precision :: g_B, g_D
  double precision, parameter :: g_B_cgs=6.172D20, g_D_cgs=3.086D21
  integer :: n_pot=10  ! exponent for smoothed potential
  character (len=labellen) :: gravx_profile='zero',gravy_profile='zero'
  character (len=labellen) :: gravz_profile='zero'
!
!  Parameters used by other modules (only defined for other gravities) 
!
  logical :: lnumerical_equilibrium=.false.
  real :: g0=0.
  real :: lnrho_bot=0.,lnrho_top=0.,ss_bot=0.,ss_top=0.
  character (len=labellen) :: grav_profile='const'

  namelist /grav_init_pars/ &
       gravx_profile,gravy_profile,gravz_profile,gravx,gravy,gravz, &
       xgrav,ygrav,zgrav,kx_gg,ky_gg,kz_gg,dgravx,nu_epicycle,pot_ratio, &
       z1,z2,zref,lnrho_bot,lnrho_top,ss_bot,ss_top, &
       lgravx_gas,lgravx_dust,lgravy_gas,lgravy_dust,lgravz_gas,lgravz_dust

  namelist /grav_run_pars/ &
       gravx_profile,gravy_profile,gravz_profile,gravx,gravy,gravz, &
       xgrav,ygrav,zgrav,kx_gg,ky_gg,kz_gg,dgravx,nu_epicycle,pot_ratio, &
       lgravx_gas,lgravx_dust,lgravy_gas,lgravy_dust,lgravz_gas,lgravz_dust, &
       zref

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_curlggrms=0,idiag_curlggmax=0,idiag_divggrms=0
  integer :: idiag_divggmax=0

  contains

!***********************************************************************
    subroutine register_gravity()
!
!  Initialise gravity variables (currently none)
!
!  12-nov-04/anders: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_gravity: called twice')
      first = .false.
!
!  Identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: gravity_simple.f90,v 1.10 2006-08-07 10:45:03 ajohan Exp $")
!
!  Set lgrav and lgravz (the latter for backwards compatibility)
!
      lgrav=.true.
      lgravz=.true.
!
    endsubroutine register_gravity
!***********************************************************************
    subroutine initialize_gravity()
!
!  Calculate master pencils for gravity. These are put into gravity pencils
!  in the subroutine calc_pencils_grav.
!
!  12-nov-04/anders: coded, copied init conds from grav_x, grav_y and grav_y.
!
      use Mpicomm, only: stop_it
      use Sub, only: notanumber, sine_step
!
      real :: ztop, prof  
!
!  Set up physical units.
!
      if (unit_system=='cgs') then
          g_A = g_A_cgs/unit_velocity*unit_time
          g_B = g_B_cgs/unit_length
          g_C = g_C_cgs/unit_velocity*unit_time
          g_D = g_D_cgs/unit_length
      else if (unit_system=='SI') then
        call stop_it('initialize_gravity: SI unit conversions not inplemented')
      endif
!
!  Different x-gravity profiles
!
      select case (gravx_profile)

      case('zero')
        if (lroot) print*,'initialize_gravity: no x-gravity'
        gravx_xpencil=0.
        lgravx_gas=.false.
        lgravx_dust=.false.

      case('const')
        if (lroot) print*,'initialize_gravity: constant x-grav=',gravx
        gravx_xpencil=gravx
        potx_xpencil=-gravx*(x(l1:l2)-xinfty)
!
!  tanh profile
!  for isothermal EOS, we have 0=-cs2*dlnrho+gravx
!  pot_ratio gives the resulting ratio in the density.
!
      case('tanh')
        if (dgravx==0.) call stop_it("initialize_gravity: dgravx=0 not OK")
        if (lroot) print*,'initialize_gravity: tanh x-grav, gravx=',gravx
        if (lroot) print*,'initialize_gravity: xgrav,dgravx=',xgrav,dgravx
        gravx=-alog(pot_ratio)/dgravx
        gravx_xpencil=gravx*.5/cosh((x(l1:l2)-xgrav)/dgravx)**2
        potx_xpencil=-gravx*.5*(1.+tanh((x(l1:l2)-xgrav)/dgravx))*dgravx

      case('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal x-grav, gravx=',gravx
        gravx_xpencil = -gravx*sin(kx_gg*x(l1:l2))

      case('kepler')
        if (lroot) print*,'initialize_gravity: kepler x-grav, gravx=',gravx
        gravx_xpencil = -gravx/x(l1:l2)**2

      case default
        if (lroot) print*, &
            'initialize_gravity: unknown gravx_profile ', gravx_profile
        call stop_it('initialize_gravity: chosen gravx_profile not valid')

      endselect
!
!  Different y-gravity profiles
!
      select case (gravy_profile)

      case('zero')
        if (lroot) print*,'initialize_gravity: no y-gravity'
        gravy_ypencil=0.
        lgravy_gas=.false.
        lgravy_dust=.false.

      case('const')
        if (lroot) print*,'initialize_gravity: constant y-grav=', gravy
        gravy_ypencil=gravy
        poty_ypencil=-gravy*(y(m1:m2)-yinfty)

      case('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal y-grav, gravy=', gravy
        gravy_ypencil = -gravy*sin(ky_gg*y(m1:m2))

      case('kepler')
        if (lroot) print*,'initialize_gravity: kepler gravy, y-grav=', gravy
        gravy_ypencil = -gravy/y(m1:m2)**2

      case default
        if (lroot) print*, &
            'initialize_gravity: unknown gravy_profile ', gravy_profile
        call stop_it('initialize_gravity: chosen gravy_profile not valid')

      endselect
!
!  Different z-gravity profiles
!
      select case (gravz_profile)

      case('zero')
        if (lroot) print*,'initialize_gravity: no z-gravity'
        gravz_zpencil=0.
        lgravz_gas=.false.
        lgravz_dust=.false.

      case('const')
        if (lroot) print*,'initialize_gravity: constant gravz=', gravz
        gravz_zpencil=gravz
        potz_zpencil=-gravz*(z(n1:n2)-zinfty)

      case('const_zero')  !  Const. gravity acc. (but zero for z>zgrav)
        if (headtt) print*,'initialize_gravity: const_zero gravz=', gravz
        if (zgrav==impossible .and. lroot) &
            print*,'initialize_gravity: zgrav is not set!'
        do n=n1,n2
          if (z(n)<=zgrav) gravz_zpencil(n-3) = gravz
        enddo

      case('linear')      !  Linear gravity profile (for accretion discs)
        nu_epicycle2=nu_epicycle**2
        if(lroot) print*,'initialize_gravity: linear z-grav, nu=', nu_epicycle
        gravz_zpencil = -nu_epicycle2*z(n1:n2)
        potz_zpencil=0.5*nu_epicycle2*(z(n1:n2)-zinfty)**2

      case('sinusoidal')
        if (lroot) print*,'initialize_gravity: sinusoidal z-grav, gravz=', gravz
        gravz_zpencil = -gravz*sin(kz_gg*z(n1:n2))
        potz_zpencil = -gravz/kz_gg*cos(kz_gg*z(n1:n2))

      case('kepler')
        if (lroot) print*,'initialize_gravity: kepler z-grav, gravz=', gravz
        gravz_zpencil = -gravz/z(n1:n2)**2

      case('Ferriere')
!
!  gravity profile from K. Ferriere, ApJ 497, 759, 1998, eq (34)
!   at solar radius.  (for interstellar runs)
!
!  nb: 331.5 is conversion factor: 10^-9 cm/s^2 -> kpc/Gyr^2)  (/= 321.1 ?!?)
!AB: These numbers should be inserted in the appropriate unuts.
!AB: As it is now, it can never make much sense.
        gravz_zpencil = -(g_A*z(n1:n2)/sqrt(z(n1:n2)**2+g_B**2) + g_C*z(n1:n2)/g_D)

      case('reduced_top') 
        if (headtt) print*,'duu_dt_grav: reduced, gravz=',gravz
        if (zgrav==impossible.and.lroot) print*,'zgrav is not set!'
        ztop = xyz0(3)+Lxyz(3)
        prof = sine_step(z(n),(zgrav+ztop)/2,(ztop-zgrav)/2)
        gravz_zpencil = (1 - prof*(1-reduced_top))*gravz

      case default
        if (lroot) print*, &
            'initialize_gravity: unknown gravz_profile ', gravz_profile
        call stop_it('initialize_gravity: chosen gravz_profile not valid')

      endselect
!
!  Sanity check
!
      if (notanumber(gravx_xpencil)) then
        call stop_it('initialize_gravity: found NaN or +/-Inf in gravx_xpencil')
      endif
      if (notanumber(gravy_ypencil)) then
        call stop_it('initialize_gravity: found NaN or +/-Inf in gravy_ypencil')
      endif
      if (notanumber(gravz_zpencil)) then
        call stop_it('initialize_gravity: found NaN or +/-Inf in gravz_zpencil')
      endif
!
    endsubroutine initialize_gravity
!***********************************************************************
    subroutine init_gg(f,xx,yy,zz)
!
!  Initialise gravity; called from start.f90
!
!  12-nov-04/anders: coded
! 
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
!  Don't do anything
!
      if (NO_WARN) print*,f,xx,yy,zz !(keep compiler quiet)
!
    endsubroutine init_gg
!***********************************************************************
    subroutine pencil_criteria_gravity()
! 
!  All pencils that the Gravity module depends on are specified here.
! 
!  20-11-04/anders: coded
!
!
      lpenc_requested(i_gg)=.true.
!
    endsubroutine pencil_criteria_gravity
!***********************************************************************
    subroutine pencil_interdep_gravity(lpencil_in)
!
!  Interdependency among pencils from the Gravity module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (NO_WARN) print*, lpencil_in !(keep compiler quiet)
!
    endsubroutine pencil_interdep_gravity
!***********************************************************************
    subroutine calc_pencils_gravity(f,p)
!
!  Calculate Gravity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  12-nov-04/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!      
      intent(in) :: f
      intent(inout) :: p
!
      if (lpencil(i_gg)) then
        p%gg(:,1) = gravx_xpencil
        p%gg(:,2) = gravy_ypencil(m-nghost)
        p%gg(:,3) = gravz_zpencil(n-nghost)
      endif
!
      if (NO_WARN) print*, f
!
    endsubroutine calc_pencils_gravity
!***********************************************************************
    subroutine duu_dt_grav(f,df,p)
!
!  Add gravitational acceleration to gas and dust
!
!  12-nov-04/anders: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      integer :: k
!
      intent(in) :: f,p
      intent(out) :: df
!
!  Add gravity acceleration on gas
!
      if (lhydro) then
        if (lgravx_gas) df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + p%gg(:,1)
        if (lgravy_gas) df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + p%gg(:,2)
        if (lgravz_gas) df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + p%gg(:,3)
      endif
!
!  Add gravity acceleration on dust
!
      if (ldustvelocity) then
        do k=1,ndustspec
          if (lgravx_dust) &
              df(l1:l2,m,n,iudx(k)) = df(l1:l2,m,n,iudx(k)) + p%gg(:,1)
          if (lgravy_dust) &
              df(l1:l2,m,n,iudy(k)) = df(l1:l2,m,n,iudy(k)) + p%gg(:,2)
          if (lgravz_dust) &
              df(l1:l2,m,n,iudz(k)) = df(l1:l2,m,n,iudz(k)) + p%gg(:,3)
        enddo
      endif
!
      if (NO_WARN) print*,f,p !(keep compiler quiet)
!        
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential_global(xx,yy,zz,pot,pot0)
!
!  Calculates gravity potential globally
!  DEPRECATED (use loop and potential_penc instead to avoid global xx,yy,zz)
!
!  13-nov-04/anders: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz) :: xx,yy,zz, pot
      real, optional :: pot0
!
      call stop_it('potential_global: this subroutine has been '// & 
          'deprecated for gravity_simple')
!
      if (NO_WARN) print*,xx(1,1,1)+yy(1,1,1)+zz(1,1,1), &
          pot(1,1,1),pot0  !(keep compiler quiet)
!
    endsubroutine potential_global
!***********************************************************************
    subroutine potential_penc(xmn,ymn,zmn,pot,pot0,grav,rmn)
!
!  Calculates gravity potential on a pencil
!
!  13-nov-04/anders: coded
!
      real, dimension (nx) :: pot
      real, optional :: ymn,zmn,pot0
      real, optional, dimension (nx) :: xmn,rmn
      real, optional, dimension (nx,3) :: grav
!
      intent(in) :: xmn,ymn,zmn,rmn
      intent(out) :: pot
!
!  Calculate potential from master pencils defined in initialize_gravity
!
      pot = potx_xpencil + poty_ypencil(m-nghost) + potz_zpencil(n-nghost)
!
      if (NO_WARN) print*, xmn, ymn, zmn, pot0, grav, rmn
!
    endsubroutine potential_penc
!***********************************************************************
    subroutine potential_point(x,y,z,r, pot,pot0, grav)
!
!  Calculates gravity potential in one point
!
!  13-nov-04/anders: coded
!
      use Mpicomm, only: stop_it
!
      real :: pot
      real, optional :: x,y,z,r
      real, optional :: pot0,grav
!
      call stop_it('potential_point: Not implemented for gravity_simple')
!
      if(NO_WARN) print*,x,y,z,r,pot,pot0,grav     !(to keep compiler quiet)
!        
    endsubroutine potential_point
!***********************************************************************
    subroutine read_gravity_init_pars(unit,iostat)
!
!  Read gravity init parameters
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=grav_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=grav_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_gravity_init_pars
!***********************************************************************
    subroutine write_gravity_init_pars(unit)
!
!  Write gravity init parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=grav_init_pars)
!
    endsubroutine write_gravity_init_pars
!***********************************************************************
    subroutine read_gravity_run_pars(unit,iostat)
!
!  Read gravity run parameters
!    
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=grav_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=grav_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_gravity_run_pars
!***********************************************************************
    subroutine write_gravity_run_pars(unit)
!
!  Write gravity run parameters
!
      integer, intent(in) :: unit
!
      write(unit,NML=grav_run_pars)
!
    endsubroutine write_gravity_run_pars
!***********************************************************************
    subroutine rprint_gravity(lreset,lwrite)
!
!  reads and registers print parameters relevant for gravity advance
!  dummy routine
!
!  12-jun-04/axel: adapted from grav_z
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!  idl needs this even if everything is zero
!
      if (lwr) then
        write(3,*) 'i_curlggrms=',idiag_curlggrms
        write(3,*) 'i_curlggmax=',idiag_curlggmax
        write(3,*) 'i_divggrms=',idiag_divggrms
        write(3,*) 'i_divggmax=',idiag_divggmax
        write(3,*) 'igg=',igg
        write(3,*) 'igx=',igx
        write(3,*) 'igy=',igy
        write(3,*) 'igz=',igz
      endif
!
      if(NO_WARN) print*,lreset  !(to keep compiler quiet)
!        
    endsubroutine rprint_gravity
!***********************************************************************

endmodule Gravity
