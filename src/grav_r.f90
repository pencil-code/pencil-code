! $Id: grav_r.f90,v 1.21 2002-07-16 21:35:22 dobler Exp $

module Gravity

!
!  Radial gravity
!

  use Cparam

  implicit none

  interface potential
    module procedure potential_global
    module procedure potential_penc
  endinterface


  ! coefficients for potential
  real, dimension (5) :: cpot = (/ 0., 0., 0., 0., 0. /)

  character (len=labellen) :: ipotential

  ! variables for compatibility with grav_z (used by Entropy and Density):
  real :: z1,z2,zref,gravz,zinfty
  character (len=labellen) :: grav_profile='const'

  namelist /grav_init_pars/ ipotential

  namelist /grav_run_pars/  ipotential



  contains

!***********************************************************************
    subroutine register_grav()
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
      if (.not. first) call stop_it('register_grav called twice')
      first = .false.
!
!  identify version number
!
      if (lroot) call cvs_id("$Id: grav_r.f90,v 1.21 2002-07-16 21:35:22 dobler Exp $")
!
      lgrav = .true.
      lgravz = .false.
      lgravr = .true.
!
    endsubroutine register_grav
!***********************************************************************
    subroutine init_grav(f,xx,yy,zz)
!
!  initialise gravity; called from start.f90
!  10-jan-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
! Not doing anything (this might change if we decide to store gg)
!
!
    endsubroutine init_grav
!***********************************************************************
    subroutine setup_grav()
!
!  Set up some variables for gravity, in particular set cpot according to
!  the value of ipotential.
!  Needed by both start.f90 and run.f90
!
!  16-jul-02/wolf: coded
!
      use Cdata
      use Mpicomm
!
!  set coefficients for potential (coefficients a0, a2, a3, b2, b3)
!  for the rational approximation
!
!              a_0   +   a_2 r^2 + a_3 r^3
!    Phi(r) = ---------------------------------------
!               1    +   b_2 r^2 + b_3 r^3 + a_3 r^4
!
      select case(ipotential)

        case ('zero')           ! zero potential
          cpot = 0.

        case ('solar')          ! solar case
          cpot = (/ 5.088, -4.344, 61.36, 10.91, -13.93 /)

        case ('M5-dwarf')       ! M5 dwarf
          cpot = (/ 2.3401, 0.44219, 2.5952, 1.5986, 0.20851 /)

        case ('simple')         ! simple potential for tests
          cpot =  (/ 1., 0., 0., 1., 0. /)

        case ('simple-2')       ! another simple potential for tests
          cpot =  (/ 1., 1., 0., 1., 1. /)

        case default
        !
        !  Catch unknown values
        !
        if (lroot) print*, 'No such value for ipotential: ', trim(ipotential)
        call stop_it("")

      endselect
!
    endsubroutine setup_grav
!***********************************************************************
    subroutine duu_dt_grav(f,df)
!
!  add duu/dt according to gravity
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use Sub
!      use Global
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: evr,gg
      real, dimension (nx) :: g_r
!
!  evr is the radial unit vector
!
      evr(:,1) = x_mn
      evr(:,2) = y_mn
      evr(:,3) = z_mn
      evr = evr / spread(r_mn+epsi,2,3)
      g_r = - r_mn * poly( (/ 2*(cpot(1)*cpot(4)-cpot(2)), &
                              3*(cpot(1)*cpot(5)-cpot(3)), &
                              4*cpot(1)*cpot(3), &
                              cpot(5)*cpot(2)-cpot(3)*cpot(4), &
                              2*cpot(2)*cpot(3), &
                              cpot(3)**2  /), r_mn) &
                   / poly( (/ 1., 0., cpot(4), cpot(5), cpot(3) /), r_mn)**2
!      g_r = - r_mn**2 * poly( (/ 3., 0., 1. /), r_mn) &
!                   / poly( (/ 1., 0., 1., 1. /), r_mn)**2
      gg = evr*spread(g_r,2,3)
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + gg
!
if (headt .and. lfirst) call output_pencil('tmp/proc0/gg.dat',gg,3)
!

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
      use Cdata, only: mx,my,mz
      use Sub, only: poly
!
      real, dimension (mx,my,mz) :: xx,yy,zz, rr, pot
      real, optional :: pot0
!
!  remove this if you are sure rr is already calculated elsewhere      
!
      rr=sqrt(xx**2+yy**2+zz**2)

      pot = - poly((/cpot(1), 0., cpot(2), cpot(3)/), rr) &
              / poly((/1., 0., cpot(4), cpot(5), cpot(3)/), rr)
!
      if (present(pot0)) then
        pot0 = cpot(1)            ! potential at r=0
      endif
!
    endsubroutine potential_global
!***********************************************************************
    subroutine potential_penc(xmn,ymn,zmn, pot,pot0, grav,rmn)
!
!  gravity potential. The grav/rmn stuff is currently not used
!
!  21-jan-02/wolf: coded
!
      use Cdata, only: nx,ny,nz
      use Sub, only: poly
!
      real, dimension (nx) :: xmn, pot
      real :: ymn,zmn
      real, optional :: pot0
      real, optional, dimension (nx) :: rmn
      real, optional, dimension (nx,3) :: grav
!      
      pot = - poly((/cpot(1), 0., cpot(2), cpot(3)/), rmn) &
              / poly((/1., 0., cpot(4), cpot(5), cpot(3)/), rmn)
!
      if (present(pot0)) then
        pot0 = cpot(1)            ! potential at r=0
      endif
!
    endsubroutine potential_penc
!***********************************************************************

endmodule Gravity
