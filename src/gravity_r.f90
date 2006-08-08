! $Id: gravity_r.f90,v 1.1 2006-08-08 15:02:16 ajohan Exp $

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
  use Messages

  implicit none

  include 'gravity.h'

  interface potential
    module procedure potential_global
    module procedure potential_penc
    module procedure potential_point
  endinterface

  ! coefficients for potential
  real, dimension(nx) :: gravx_pencil=0.,gravy_pencil=0.,gravz_pencil=0.
  real, dimension (5) :: cpot = (/ 0., 0., 0., 0., 0. /)
  real :: nu_epicycle=1.
  real :: lnrho_bot,lnrho_top,ss_bot,ss_top
  real :: grav_const=1.,reduced_top=1.
  real :: g0=0.
  real :: r0_pot=0.    ! peak radius for smoothed potential
  integer :: n_pot=10  ! exponent for smoothed potential

  character (len=labellen) :: ipotential='zero'

  ! variables for compatibility with grav_z (used by Entropy and Density):
  real :: z1,z2,zref,zgrav,gravz,zinfty
  character (len=labellen) :: grav_profile='const'
  logical :: lnumerical_equilibrium=.false.

  namelist /grav_init_pars/ ipotential,g0,r0_pot,n_pot,lnumerical_equilibrium
  
  namelist /grav_run_pars/  ipotential,g0,r0_pot,n_pot,lnumerical_equilibrium

  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_curlggrms=0,idiag_curlggmax=0,idiag_divggrms=0
  integer :: idiag_divggmax=0

  contains

!***********************************************************************
    subroutine register_gravity()
!
!  initialise gravity flags
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_grav: called twice')
      first = .false.
!
!  identify version number
!
      if (lroot) call cvs_id("$Id: gravity_r.f90,v 1.1 2006-08-08 15:02:16 ajohan Exp $")
!
      lgrav =.true.
      lgravr=.true.
      lgravr_gas =.true.
      lgravr_dust=.true.
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
      use Cdata
      use Sub, only: poly, calc_unitvects_sphere
      use Mpicomm
      use Global
!
      real, dimension (nx,3) :: gg_mn=0.
      real, dimension (nx) :: g_r,rr_mn

      logical, save :: first=.true.
      logical :: lpade=.true. ! set to false for 1/r potential

      integer :: i

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

        case ('A0-star')       ! A0 star 
          if (lroot) print*, 'A0 star gravity potential'
!          cpot = (/ 4.7446,  -1.9456,  0.6884,  4.8007, 1.79877 /)
          cpot = (/ 4.3641,  -1.5612,  0.4841, 4.0678, 1.2548 /)
 
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
          cpot = (/ 0., -.5, 0., 0., 0. /)
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
      if (lnumerical_equilibrium) then

        if (lroot) then
          print*,'inititialize_gravity: numerical exact equilibrium -- gravity'
          print*,'                      will be calculated in density module'
        endif

      else

        do n=n1,n2
        do m=m1,m2
!
           if (lcylindrical) then
              rr_mn=sqrt(x(l1:l2)**2+y(m)**2)+tini
           else   
              rr_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)+tini
           endif
!
          if (lpade) then

            g_r = - rr_mn * poly( (/ 2*(cpot(1)*cpot(4)-cpot(2)), &
                                     3*(cpot(1)*cpot(5)-cpot(3)), &
                                      4*cpot(1)*cpot(3), &
                                        cpot(5)*cpot(2)-cpot(3)*cpot(4), &
                                      2*cpot(2)*cpot(3), &
                                        cpot(3)**2  /), rr_mn) &
                         / poly( (/ 1., 0., cpot(4), cpot(5), &
                                            cpot(3) /), rr_mn)**2

          else

             ! smoothed 1/r potential in a spherical shell
             g_r=-g0*rr_mn**(n_pot-1) &
                  *(rr_mn**n_pot+r0_pot**n_pot)**(-1./n_pot-1.)
!
          endif
!
          gg_mn(:,1) = x(l1:l2)/rr_mn*g_r
          gg_mn(:,2) = y(  m  )/rr_mn*g_r
          gg_mn(:,3) = z(  n  )/rr_mn*g_r
          if (lcylindrical) gg_mn(:,3)=0.
         
          call set_global(gg_mn,m,n,'gg',nx)

          enddo
        enddo
!
      endif

      endif
!
    endsubroutine initialize_gravity
!***********************************************************************
    subroutine read_gravity_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=grav_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=grav_init_pars,ERR=99)
      endif
                                                                                                                                                                                                
99    return
    endsubroutine read_gravity_init_pars
!***********************************************************************
    subroutine write_gravity_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=grav_init_pars)

    endsubroutine write_gravity_init_pars
!***********************************************************************
    subroutine read_gravity_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=grav_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=grav_run_pars,ERR=99)
      endif
                                                                                                                                                                                                
99    return
    endsubroutine read_gravity_run_pars
!***********************************************************************
    subroutine write_gravity_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=grav_run_pars)

    endsubroutine write_gravity_run_pars
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
      if (NO_WARN) print*,f,xx,yy,zz  !(to keep compiler quiet)
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
      intent(in) :: f,p
!
      if (NO_WARN) print*, f, p  !(keep compiler quiet)
!
    endsubroutine calc_pencils_gravity
!***********************************************************************
    subroutine duu_dt_grav(f,df,p)
!
!  add duu/dt according to gravity
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use Sub
      use Global
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: gg_mn
      real, dimension (nx) :: g_r
!
!  evr is the radial unit vector
!
!if (.false.) then               ! switch between the two methods for timing
!if (.true.) then               ! switch between the two methods for timing
!
!       g_r = - r_mn * poly( (/ 2*(cpot(1)*cpot(4)-cpot(2)), &
!                              3*(cpot(1)*cpot(5)-cpot(3)), &
!                              4*cpot(1)*cpot(3), &
!                              cpot(5)*cpot(2)-cpot(3)*cpot(4), &
!                              2*cpot(2)*cpot(3), &
!                              cpot(3)**2  /), r_mn) &
!                   / poly( (/ 1., 0., cpot(4), cpot(5), cpot(3) /), r_mn)**2
!      g_r = - r_mn**2 * poly( (/ 3., 0., 1. /), r_mn) &
!                   / poly( (/ 1., 0., 1., 1. /), r_mn)**2
!
! dgm: radial unit vector evr now calculated in subroutine 
! calc_unitvects_sphere()
!
!      gg_mn = evr*spread(g_r,2,3)
!      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + gg_mn
!else

      call get_global(gg_mn,m,n,'gg')
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + gg_mn
!
!
!endif
!
! if (headtt) call output_pencil(trim(datadir)//'/proc0/gg0.dat',gg_mn,3)
!
      if (NO_WARN) print*,f,p !(to keep compiler quiet)
!        
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
      use Cdata, only: mx,my,mz,dx
      use Sub, only: poly

      real, dimension (mx,my,mz) :: xx,yy,zz, pot
      real, optional :: pot0           ! potential at r=0

      real, dimension (mx,my,mz) :: rr
!
!  remove this if you are sure rr is already calculated elsewhere      
!
      rr=sqrt(xx**2+yy**2+zz**2)

      select case (ipotential)

      case ('geo-kws','smoothed-newton')
        pot = -g0*(rr**n_pot+r0_pot**n_pot)**(-1.0/n_pot)
        if (present(pot0)) pot0=-g0/r0_pot

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
      use Cdata, only: nx,dx
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
        if (present(pot0)) pot0=-g0/r0_pot

      case default
        pot = - poly((/cpot(1), 0., cpot(2), cpot(3)/), rad) &
                / poly((/1., 0., cpot(4), cpot(5), cpot(3)/), rad)
        if (present(pot0)) pot0=-cpot(1)

      endselect

      if (present(grav)) call stop_it("POTENTIAL_PENC: Argument grav"//&
                                      "not implemented")

    endsubroutine potential_penc
!***********************************************************************
    subroutine potential_point(x,y,z,r, pot,pot0, grav)
!
!  Gravity potential in one point
!
!  20-dec-03/wolf: coded
!
      use Cdata, only: dx
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
        if (present(pot0)) pot0=-g0/r0_pot

      case default
        pot = - poly((/cpot(1), 0., cpot(2), cpot(3)/), rad) &
                / poly((/1., 0., cpot(4), cpot(5), cpot(3)/), rad)
        if (present(pot0)) pot0=-cpot(1)

      endselect

      if (present(grav)) call stop_it("POTENTIAL_PENC: Argument grav"//&
                                      "not implemented")

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
      use Sub
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
      if (NO_WARN) print*,lreset  !(to keep compiler quiet)
!
    endsubroutine rprint_gravity
!***********************************************************************
endmodule Gravity
