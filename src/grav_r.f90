! $Id: grav_r.f90,v 1.64 2004-06-18 12:19:03 mcmillan Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Gravity

!
!  Radial gravity
!

  use Cparam

  implicit none

  interface potential
    module procedure potential_global
    module procedure potential_penc
    module procedure potential_point
  endinterface

  ! coefficients for potential
  real, dimension (5) :: cpot = (/ 0., 0., 0., 0., 0. /)
  real :: nu_epicycle=1.
  real :: lnrho_bot,lnrho_top,ss_bot,ss_top
  real :: grav_const=1.
  real :: g0=0.
  real :: r0_pot=0.    ! peak radius for smoothed potential
  integer :: n_pot=10  ! exponent for smoothed potential

  character (len=labellen) :: ipotential='zero'

  ! variables for compatibility with grav_z (used by Entropy and Density):
  real :: z1,z2,zref,zgrav,gravz,zinfty
  character (len=labellen) :: grav_profile='const'

  namelist /grav_init_pars/ ipotential,g0,r0_pot,n_pot

  namelist /grav_run_pars/  ipotential,g0,r0_pot,n_pot

  ! other variables (needs to be consistent with reset list below)
  integer :: i_curlggrms=0,i_curlggmax=0,i_divggrms=0,i_divggmax=0

  contains

!***********************************************************************
    subroutine register_gravity()
!
!  initialise gravity flags
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_grav: called twice')
      first = .false.
!
!  identify version number
!
      if (lroot) call cvs_id("$Id: grav_r.f90,v 1.64 2004-06-18 12:19:03 mcmillan Exp $")
!
      lgrav = .true.
      lgravz = .false.
      lgravr = .true.
!
    endsubroutine register_gravity
!***********************************************************************
    subroutine initialize_gravity()
!
!  Set up cpot according to the value of ipotential, and initialize the
!  global variable gg (gravity field).
!  Needed by both start.f90 and run.f90
!
!  16-jul-02/wolf: coded
!  22-nov-02/tony: renamed
!
!ajwm - need to figure out how to call this from start.f90
      use Cdata
      use Sub, only: poly, calc_unitvects_sphere
      use Mpicomm
      use Global
!
      real, dimension (nx,3) :: gg_mn
      real, dimension (nx) :: g_r

      logical, save :: first=.true.
      logical :: lpade=.true. ! set to false for 1/r potential

      !ajwm - should this be done on RELOAD too??
      if (first) then
!
!  for lpade=.true. set coefficients for potential (coefficients a0, a2, a3,
!  b2, b3) for the rational approximation
!
!              a_0   +   a_2 r^2 + a_3 r^3
!    Phi(r) = ---------------------------------------
!               1    +   b_2 r^2 + b_3 r^3 + a_3 r^4
!
      select case(ipotential)

        case ('zero')           ! zero potential
          if (lroot) print*, 'initialize_gravity: zero gravity potential'
          cpot = 0.

        case ('solar')          ! solar case
          if (lroot) print*, 'initialize_gravity: solar gravity potential'
          cpot = (/ 5.088, -4.344, 61.36, 10.91, -13.93 /)

        case ('M5-dwarf')       ! M5 dwarf
          if (lroot) print*, 'initialize_gravity: M5 dwarf gravity potential'
          cpot = (/ 2.3401, 0.44219, 2.5952, 1.5986, 0.20851 /)

        case ('M2-sgiant')       ! M super giant
          if (lroot) print*, 'M super giant gravity potential'
          cpot = (/ 1.100, 0.660, 2.800, 1.400, 0.100 /)

        case ('A7-star')       ! Ap star 
          if (lroot) print*, 'A star gravity potential'
          cpot = (/ 4.080, -3.444, 15.2000, 11.2000, -12.1000 /)

        case ('simple')         ! simple potential for tests
          if (lroot) print*, 'initialize_gravity: very simple gravity potential'
          cpot =  (/ 1., 0., 0., 1., 0. /)

        case ('simple-2')       ! another simple potential for tests
          if (lroot) print*, 'initialize_gravity: simple gravity potential'
          cpot =  (/ 1., 1., 0., 1., 1. /)

        case ('smoothed-newton')
          if (lroot) print*,'initialize_gravity: smoothed 1/r potential'
          lpade=.false.

        ! geodynamo
        case ('geo-kws-approx')     ! approx. 1/r potential between r=.5 and r=1
          if (lroot) print*, 'initialize_gravity: approximate 1/r potential'
          cpot = (/ 0., 2.2679, 0., 0., 1.1697 /)
        case ('geo-benchmark')      ! for geodynamo benchmark runs
          if (lroot) print*, 'initialize_gravity: gravity linear in radius'
          cpot = (/ 0., .5, 0., 0., 0. /)
        case ('geo-kws')
          if (lroot) print*, 'initialize_gravity: '//&
                             'smoothed 1/r potential in spherical shell'
          if (r0_pot < epsi) print*, 'WARNING: grav_r: r0_pot is too small.'//&
                                     'Can be set in grav_r namelists.'
          lpade=.false.
        ! end geodynamo

        case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'initialize_gravity: '//&
                           'No such value for ipotential: ', trim(ipotential)
        call stop_it("")
        
      endselect
!
!  initialize gg, so we can later retrieve gravity via get_global
!
      do imn=1,ny*nz

        n=nn(imn)
        m=mm(imn)
!
        call calc_unitvects_sphere()
!
        if (lpade) then

          g_r = - r_mn * poly( (/ 2*(cpot(1)*cpot(4)-cpot(2)), &
                                  3*(cpot(1)*cpot(5)-cpot(3)), &
                                  4*cpot(1)*cpot(3), &
                                  cpot(5)*cpot(2)-cpot(3)*cpot(4), &
                                  2*cpot(2)*cpot(3), &
                                  cpot(3)**2  /), r_mn) &
                       / poly( (/ 1., 0., cpot(4), cpot(5), cpot(3) /), r_mn)**2

        else

          ! smoothed 1/r potential in a spherical shell
          g_r=-g0*r_mn**(n_pot-1)*(r_mn**n_pot+r0_pot**n_pot)**(-1.0/n_pot-1.0)

        endif
!
        gg_mn = evr*spread(g_r,2,3)
        call set_global(gg_mn,m,n,'gg')

      enddo
!
      endif

    endsubroutine initialize_gravity
!***********************************************************************
    subroutine init_gg(f,xx,yy,zz)
!
!  initialise gravity; called from start.f90
!  10-jan-02/wolf: coded
!  24-nov-02/tony: renamed from init_grav for consistancy (i.e. init_[variable name])
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
! Not doing anything (this might change if we decide to save gg to a file)
!
      if(ip==0) print*,f,xx,yy,zz  !(to keep compiler quiet)
    endsubroutine init_gg
!***********************************************************************
    subroutine duu_dt_grav(f,df,uu,rho1)
!
!  add duu/dt according to gravity
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use Sub
      use IO
      use Global
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: gg_mn,uu
      real, dimension (nx) :: g_r,rho1
!
!  evr is the radial unit vector
!
if (.false.) then               ! switch between the two methods for timing
!if (.true.) then               ! switch between the two methods for timing
!
       g_r = - r_mn * poly( (/ 2*(cpot(1)*cpot(4)-cpot(2)), &
                              3*(cpot(1)*cpot(5)-cpot(3)), &
                              4*cpot(1)*cpot(3), &
                              cpot(5)*cpot(2)-cpot(3)*cpot(4), &
                              2*cpot(2)*cpot(3), &
                              cpot(3)**2  /), r_mn) &
                   / poly( (/ 1., 0., cpot(4), cpot(5), cpot(3) /), r_mn)**2
!      g_r = - r_mn**2 * poly( (/ 3., 0., 1. /), r_mn) &
!                   / poly( (/ 1., 0., 1., 1. /), r_mn)**2
!
! dgm: radial unit vector evr now calculated in subroutine 
! calc_unitvects_sphere()
!
      gg_mn = evr*spread(g_r,2,3)
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + gg_mn
else
      call get_global(gg_mn,m,n,'gg')
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + gg_mn
endif
!
! if (headtt) call output_pencil(trim(datadir)//'/proc0/gg0.dat',gg_mn,3)
!
      if(ip==0) print*,f,uu,rho1  !(to keep compiler quiet)
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential_global(xx,yy,zz, pot,pot0)
!    subroutine potential(rr, pot)
!
!  gravity potential; version called by init_hydro, which operates on
!  full global coordinate arrays
!
!  16-jul-02/wolf: coded
!
      use Cdata, only: mx,my,mz
      use Sub, only: poly

      real, dimension (mx,my,mz) :: xx,yy,zz, rr, pot
      real, optional :: pot0           ! potential ar r=0
!
!  remove this if you are sure rr is already calculated elsewhere      
!
      rr=sqrt(xx**2+yy**2+zz**2)

      select case (ipotential)

      case ('geo-kws','smoothed-newton')
        pot=-g0*(rr**n_pot+r0_pot**n_pot)**(-1.0/n_pot)
        if (present(pot0)) pot0=-1.0/r0_pot

      case default
        pot = - poly((/cpot(1), 0., cpot(2), cpot(3)/), rr) &
                / poly((/1., 0., cpot(4), cpot(5), cpot(3)/), rr)
        if (present(pot0)) pot0=-cpot(1)

      endselect

    endsubroutine potential_global
!***********************************************************************
    subroutine potential_penc(xmn,ymn,zmn,rmn, pot,pot0, grav)
!
!  Gravity potential along one pencil
!
!  21-jan-02/wolf: coded
!
      use Cdata, only: nx
      use Sub, only: poly
      use Mpicomm, only: stop_it

      real, dimension (nx) :: rad, pot
      real, optional :: ymn,zmn,pot0
      real, optional, dimension (nx) :: xmn,rmn
      real, optional, dimension (nx,3) :: grav
      
      if (present(rmn)) then
        rad = rmn
      else
        if (present(xmn) .and. present(ymn) .and. present(zmn)) then
          rad = sqrt(xmn**2+ymn**2+zmn**2)
        else
          call stop_it("POTENTIAL_PENC: Need to specify either x,y,z or r.")
        endif
      endif

      select case (ipotential)

      case ('geo-kws','smoothed-newton')
        pot=-g0*(rmn**n_pot+r0_pot**n_pot)**(-1.0/n_pot)
        if (present(pot0)) pot0=-1.0/r0_pot

      case default
        pot = - poly((/cpot(1), 0., cpot(2), cpot(3)/), rmn) &
                / poly((/1., 0., cpot(4), cpot(5), cpot(3)/), rmn)
        if (present(pot0)) pot0=-cpot(1)

      endselect

      if(ip==0) print*,present(grav)

    endsubroutine potential_penc
!***********************************************************************
    subroutine potential_point(x,y,z,r, pot,pot0, grav)
!
!  Gravity potential in one point
!
!  20-dec-03/wolf: coded
!
      use Sub, only: poly
      use Mpicomm, only: stop_it

      real :: pot,rad
      real, optional :: x,y,z,r
      real, optional :: pot0,grav

      if (present(r)) then
        rad = r
      else
        if (present(x) .and. present(y) .and. present(z)) then
          rad = sqrt(x**2+y**2+z**2)
        else
          call stop_it("Need to specify either x,y,z or r in potential_point()")
        endif
      endif

      select case (ipotential)

      case ('geo-kws','smoothed-newton')
        pot=-g0*(rad**n_pot+r0_pot**n_pot)**(-1.0/n_pot)
        if (present(pot0)) pot0=-1.0/r0_pot

      case default
        pot = - poly((/cpot(1), 0., cpot(2), cpot(3)/), rad) &
                / poly((/1., 0., cpot(4), cpot(5), cpot(3)/), rad)
        if (present(pot0)) pot0=-cpot(1)

      endselect

      if(ip==0) print*,present(grav)

    endsubroutine potential_point
!***********************************************************************
    subroutine rprint_gravity(lreset,lwrite)
!
!  reads and registers print parameters relevant for gravity advance
!  dummy routine
!
!  26-apr-03/axel: coded
!
      use Cdata
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  write column, i_XYZ, where our variable XYZ is stored
!  idl needs this even if everything is zero
!
      if (lwr) then
        write(3,*) 'i_curlggrms=',i_curlggrms
        write(3,*) 'i_curlggmax=',i_curlggmax
        write(3,*) 'i_divggrms=',i_divggrms
        write(3,*) 'i_divggmax=',i_divggmax
        write(3,*) 'igg=',igg
        write(3,*) 'igx=',igx
        write(3,*) 'igy=',igy
        write(3,*) 'igz=',igz
      endif
!
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_gravity
!***********************************************************************

endmodule Gravity
