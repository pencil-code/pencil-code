module Grid

  implicit none

  interface grid_profile        ! Overload the grid_profile' subroutine
    module procedure grid_profile_point
    module procedure grid_profile_1d
  endinterface

  contains
!***********************************************************************
    subroutine construct_grid(x,y,z,dx,dy,dz,x00,y00,z00)
!
!  Constructs a non-equidistant grid x(xi) of an equidistant grid xi with
!  grid spacing dxi=1. For grid_func='linear' this is equivalent to an
!  equidistant grid.
!
!  dx_1 and dx_tilde are the coefficients that enter the formulae for the
!  1st and 2nd derivative:
!
!    ``df/dx'' = ``df/dxi'' * dx_1
!    ``d2f/dx2'' = ``df2/dxi2'' * dx_1**2 + dx_tilde * ``df/dxi'' * dx_1
!
!  These coefficients are also very useful when adapting grid dependend stuff
!  such as the timestep. A simple substitution
!    1./dx -> dx_1
!  should suffice in most cases.
!
!  Currently, the only non-equidistant grid-function is sinh. Its inflection
!  point in each direction is specified by xyz_star.
!  (Suggestions for a better name than ``xyz_star'' are welcome.)
!
!
!  25-jun-04/tobi+wolf: coded
!
      use Cdata, only: nx,ny,nz
      use Cdata, only: mx,my,mz
      use Cdata, only: Lx,Ly,Lz
      use Cdata, only: dx_1,dy_1,dz_1
      use Cdata, only: dx_tilde,dy_tilde,dz_tilde
      use Cdata, only: x0,y0,z0
      use Cdata, only: nxgrid,nygrid,nzgrid
      use Cdata, only: ipx,ipy,ipz
      use Cdata, only: nghost,coeff_grid,grid_func
      use Cdata, only: lperi,lshift_origin,xyz_star,lequidist
      use Mpicomm, only: stop_it

      real, dimension(mx), intent(out) :: x
      real, dimension(my), intent(out) :: y
      real, dimension(mz), intent(out) :: z
      real, intent(out) :: dx,dy,dz,x00,y00,z00

      real :: xi1lo,xi1up,g1lo,g1up
      real :: xi2lo,xi2up,g2lo,g2up
      real :: xi3lo,xi3up,g3lo,g3up
      real :: xi1star,xi2star,xi3star

      real, dimension(mx) :: g1,g1der1,g1der2,xi1,xprim,xprim2
      real, dimension(my) :: g2,g2der1,g2der2,xi2,yprim,yprim2
      real, dimension(mz) :: g3,g3der1,g3der2,xi3,zprim,zprim2

      real :: a,dummy
      integer :: i
      logical :: err


      lequidist=(grid_func=='linear')

      if (lperi(1)) then
        dx=Lx/nxgrid
        x00=x0+.5*dx
      else
        dx=Lx/(nxgrid-1)
        x00=x0
        if (lshift_origin(1)) x00=x0+.5*dx
      endif
      if (lperi(2)) then
        dy=Ly/nygrid
        y00=y0+.5*dy
      else
        dy=Ly/(nygrid-1)
        y00=y0
        if (lshift_origin(2)) y00=y0+.5*dy
      endif
      if (lperi(3)) then
        dz=Lz/nzgrid
        z00=z0+.5*dz
      else
        dz=Lz/(nzgrid-1)
        z00=z0
        if (lshift_origin(3)) z00=z0+.5*dz
      endif

      do i=1,mx; xi1(i)=i-nghost-1+ipx*nx; enddo
      do i=1,my; xi2(i)=i-nghost-1+ipy*ny; enddo
      do i=1,mz; xi3(i)=i-nghost-1+ipz*nz; enddo
!
!  The following is correct for periodic and non-periodic case
!
      xi1lo=0; xi1up=nxgrid-merge(0,1,lperi(1))
      xi2lo=0; xi2up=nygrid-merge(0,1,lperi(2))
      xi3lo=0; xi3up=nzgrid-merge(0,1,lperi(3))
!
!  Construct nonequidistant grid
!
!  x coordinate
!
      if (nxgrid == 1) then
        x = x00
        ! hopefully, we will only ever multiply by the following quantities:
        xprim = 0.
        xprim2 = 0.
        dx_1 = 0.
        dx_tilde = 0.
      else
        call grid_profile(dummy,grid_func(1),dummy,err=err)
        if (err) call &
             stop_it("CONSTRUCT_GRID: unknown grid_func "//grid_func(1))

        a=coeff_grid(1,1)*dx
        xi1star=find_star(a*xi1lo,a*xi1up,x00,x00+Lx,xyz_star(1),grid_func(1))/a

        call grid_profile(a*(xi1  -xi1star),grid_func(1),g1,g1der1,g1der2)
        call grid_profile(a*(xi1lo-xi1star),grid_func(1),g1lo)
        call grid_profile(a*(xi1up-xi1star),grid_func(1),g1up)

        x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
        xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
        xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)

        dx_1=1./xprim
        dx_tilde=-xprim2/xprim**2
      endif
!
!  y coordinate
!
      if (nygrid == 1) then
        y = y00
        ! hopefully, we will only ever multiply by the following quantities:
        yprim = 0.
        yprim2 = 0.
        dy_1 = 0.
        dy_tilde = 0.
      else
        call grid_profile(dummy,grid_func(2),dummy,err=err)
        if (err) &
             call stop_it("CONSTRUCT_GRID: unknown grid_func "//grid_func(2))

        a=coeff_grid(2,1)*dy
        xi2star=find_star(a*xi2lo,a*xi2up,y00,y00+Ly,xyz_star(2),grid_func(2))/a

        call grid_profile(a*(xi2  -xi2star),grid_func(2),g2,g2der1,g2der2)
        call grid_profile(a*(xi2lo-xi2star),grid_func(2),g2lo)
        call grid_profile(a*(xi2up-xi2star),grid_func(2),g2up)

        y     =y00+Ly*(g2  -  g2lo)/(g2up-g2lo)
        yprim =    Ly*(g2der1*a   )/(g2up-g2lo)
        yprim2=    Ly*(g2der2*a**2)/(g2up-g2lo)

        dy_1=1./yprim
        dy_tilde=-yprim2/yprim**2
      endif
!
!  z coordinate
!
      if (nzgrid == 1) then
        z = z00
        ! hopefully, we will only ever multiply by the following quantities:
        zprim = 0.
        zprim2 = 0.
        dz_1 = 0.
        dz_tilde = 0.
      else
        call grid_profile(dummy,grid_func(3),dummy,err=err)
        if (err) &
             call stop_it("CONSTRUCT_GRID: unknown grid_func "//grid_func(3))

        a=coeff_grid(3,1)*dz
        xi3star=find_star(a*xi3lo,a*xi3up,z00,z00+Lz,xyz_star(3),grid_func(3))/a

        call grid_profile(a*(xi3  -xi3star),grid_func(3),g3,g3der1,g3der2)
        call grid_profile(a*(xi3lo-xi3star),grid_func(3),g3lo)
        call grid_profile(a*(xi3up-xi3star),grid_func(3),g3up)

        z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
        zprim =    Lz*(g3der1*a   )/(g3up-g3lo)
        zprim2=    Lz*(g3der2*a**2)/(g3up-g3lo)

        dz_1=1./zprim
        dz_tilde=-zprim2/zprim**2
      endif

    endsubroutine construct_grid
!***********************************************************************
    subroutine grid_profile_point(xi,grid_func,g,gder1,gder2,err)
!
!  Specify the functional form of the grid profile function g
!  and calculate g,g',g''.
!  Must be in sync with grid_profile_2.
!  Much nicer as one `elemental subroutine', but some of our compilers
!  don't speak F95 fluently
!
!  25-jun-04/tobi+wolf: coded
!
      real              :: xi
      character(len=*)  :: grid_func
      real              :: g
      real, optional    :: gder1,gder2
      logical, optional :: err
!
      intent(in)  :: xi,grid_func
      intent(out) :: g, gder1,gder2,err

      if (present(err)) err=.false.

      select case (grid_func)

      case ('linear')
        g=xi
        if (present(gder1)) gder1=1.0
        if (present(gder2)) gder2=0.0

      case ('sinh')
        g=sinh(xi)
        if (present(gder1)) gder1=cosh(xi)
        if (present(gder2)) gder2=sinh(xi)

      case default
        if (present(err)) err=.true.

      endselect

    endsubroutine grid_profile_point
!***********************************************************************
    subroutine grid_profile_1d(xi,grid_func,g,gder1,gder2,err)
!
!  Same as grid_profile_1 for 1d arrays as arguments
!
!  25-jun-04/tobi+wolf: coded
!
      real, dimension(:)                    :: xi
      character(len=*)                      :: grid_func
      real, dimension(size(xi,1))           :: g
      real, dimension(size(xi,1)), optional :: gder1,gder2
      logical, optional                     :: err
!
      intent(in)  :: xi,grid_func
      intent(out) :: g, gder1,gder2,err

      if (present(err)) err=.false.

      select case (grid_func)

      case ('linear')
        g=xi
        if (present(gder1)) gder1=1.0
        if (present(gder2)) gder2=0.0

      case ('sinh')
        g=sinh(xi)
        if (present(gder1)) gder1=cosh(xi)
        if (present(gder2)) gder2=sinh(xi)

      case default
        if (present(err)) err=.true.

      endselect

    endsubroutine grid_profile_1d
!***********************************************************************
    function find_star(xi_lo,xi_up,x_lo,x_up,x_star,grid_func) result (xi_star)
!
!  Finds the xi that corresponds to the inflection point of the grid-function
!  by means of a newton-raphson root-finding algorithm.
!
!
!  25-jun-04/tobi+wolf: coded
!
      use Cdata, only: epsi
      use Mpicomm, only: stop_it

      real, intent(in) :: xi_lo,xi_up,x_lo,x_up,x_star
      character(len=*), intent(in) :: grid_func

      real :: xi_star,dxi,tol
      real :: g_lo,gder_lo
      real :: g_up,gder_up
      real :: f   ,fder
      integer, parameter :: maxit=100
      logical :: lreturn
      integer :: it


      if (xi_lo>=xi_up) &
           call stop_it("FIND_STAR: xi1 >= xi2 -- this should not happen")

      tol=epsi*(xi_up-xi_lo)
      xi_star= (xi_up+xi_lo)/2

      lreturn=.false.

      do it=1,maxit

        call grid_profile(xi_lo-xi_star,grid_func,g_lo,gder_lo)
        call grid_profile(xi_up-xi_star,grid_func,g_up,gder_up)

        f   =-(x_up-x_star)*g_lo   +(x_lo-x_star)*g_up
        fder= (x_up-x_star)*gder_lo-(x_lo-x_star)*gder_up

        dxi=f/fder
        xi_star=xi_star-dxi

        if (lreturn) return

        if (abs(dxi)<tol) lreturn=.true.
        
      enddo

      call stop_it("FIND_STAR: maximum number of iterations exceeded")

    endfunction find_star
!***********************************************************************
endmodule Grid
