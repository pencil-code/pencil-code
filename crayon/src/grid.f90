! $Id: grid.f90 12795 2010-01-03 14:03:57Z ajohan@strw.leidenuniv.nl $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! PENCILS PROVIDED x_mn; y_mn; z_mn; r_mn; r_mn1
! PENCILS PROVIDED phix; phiy
! PENCILS PROVIDED pomx; pomy
! PENCILS PROVIDED rcyl_mn; rcyl_mn1; phi_mn
! PENCILS PROVIDED evr(3); rr(3); evth(3)
!
!***************************************************************
module Grid
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  private
!
  public :: construct_grid
  public :: pencil_interdep_grid
  public :: calc_pencils_grid
!
  interface grid_profile        ! Overload the grid_profile' subroutine
    module procedure grid_profile_point
    module procedure grid_profile_1d
  endinterface
!
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
!  25-jun-04/tobi+wolf: coded
!
      real, dimension(mx), intent(out) :: x
      real, dimension(my), intent(out) :: y
      real, dimension(mz), intent(out) :: z
      real, intent(out) :: dx,dy,dz,x00,y00,z00
!
      real :: xi1lo,xi1up,g1lo,g1up
      real :: xi2lo,xi2up,g2lo,g2up
      real :: xi3lo,xi3up,g3lo,g3up
      real :: dxmin_x,dxmax_x,dxmin_y,dxmax_y,dxmin_z,dxmax_z
      real :: xi1star,xi2star,xi3star
!
      real, dimension(mx) :: g1,g1der1,g1der2,xi1,xprim2
      real, dimension(my) :: g2,g2der1,g2der2,xi2,yprim2
      real, dimension(mz) :: g3,g3der1,g3der2,xi3,zprim2
!
      real, dimension(0:2*nprocx+1) :: g1proc
      real, dimension(0:2*nprocy+1) :: g2proc
      real, dimension(0:2*nprocz+1) :: g3proc
!
      real :: a,dummy1=0.,dummy2=0.
      integer :: i
!
      if (lperi(1)) then
        dx=Lx/nxgrid
        x00=x0+.5*dx
        if (lshift_origin(1)) x00=x0+dx
      else
        dx=Lx/(nxgrid-1)
        x00=x0
        if (lshift_origin(1)) x00=x0+.5*dx
      endif
      if (lperi(2)) then
        dy=Ly/nygrid
        y00=y0+.5*dy
        if (lshift_origin(2)) y00=y0+dy
      else
        dy=Ly/(nygrid-1)
        y00=y0
        if (lshift_origin(2)) y00=y0+.5*dy
      endif
      if (lperi(3)) then
        dz=Lz/nzgrid
        z00=z0+.5*dz
        if (lshift_origin(3)) z00=z0+dz
      else
        dz=Lz/(nzgrid-1)
        z00=z0
        if (lshift_origin(3)) z00=z0+.5*dz
      endif
!
!  produce index arrays xi1, xi2, and xi3
!
      do i=1,mx; xi1(i)=i-nghost-1+ipx*nx; enddo
      do i=1,my; xi2(i)=i-nghost-1+ipy*ny; enddo
      do i=1,mz; xi3(i)=i-nghost-1+ipz*nz; enddo
!
!  The following is correct for periodic and non-periodic case
!
      xi1lo=0.; xi1up=nxgrid-merge(0.,1.,lperi(1))
      xi2lo=0.; xi2up=nygrid-merge(0.,1.,lperi(2))
      xi3lo=0.; xi3up=nzgrid-merge(0.,1.,lperi(3))
!
!  Construct nonequidistant grid
!
!  x coordinate
!
      if (nxgrid==1) then
        x = x00
        ! hopefully, we will only ever multiply by the following quantities:
        xprim = 0.
        xprim2 = 0.
        dx_1 = 0.
        dx_tilde = 0.
        g1proc=x00
      else
        ! Test whether grid function is valid
        call grid_profile(dummy1,dummy2)
!
        a=coeff_grid(1,1)*dx
        xi1star=find_star(a*xi1lo,a*xi1up,x00,x00+Lx,xyz_star(1))/a
        call grid_profile(a*(xi1  -xi1star),g1,g1der1,g1der2)
        call grid_profile(a*(xi1lo-xi1star),g1lo)
        call grid_profile(a*(xi1up-xi1star),g1up)
!
        x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
        xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
        xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
!
        dx_1=1./xprim
        dx_tilde=-xprim2/xprim**2
      endif
!
!  y coordinate
!
      if (nygrid==1) then
        y = y00
        ! hopefully, we will only ever multiply by the following quantities:
        yprim = 0.
        yprim2 = 0.
        dy_1 = 0.
        dy_tilde = 0.
        g2proc=y00
      else
!
        a=coeff_grid(2,1)*dy
        xi2star=find_star(a*xi2lo,a*xi2up,y00,y00+Ly,xyz_star(2))/a
        call grid_profile(a*(xi2  -xi2star),g2,g2der1,g2der2)
        call grid_profile(a*(xi2lo-xi2star),g2lo)
        call grid_profile(a*(xi2up-xi2star),g2up)
!
        y     =y00+Ly*(g2  -  g2lo)/(g2up-g2lo)
        yprim =    Ly*(g2der1*a   )/(g2up-g2lo)
        yprim2=    Ly*(g2der2*a**2)/(g2up-g2lo)
!
        dy_1=1./yprim
        dy_tilde=-yprim2/yprim**2
!
      endif
!
!  z coordinate
!
      if (nzgrid==1) then
        z = z00
        ! hopefully, we will only ever multiply by the following quantities:
        zprim = 0.
        zprim2 = 0.
        dz_1 = 0.
        dz_tilde = 0.
        g3proc=z00
      else
        a=coeff_grid(3,1)*dz
        xi3star=find_star(a*xi3lo,a*xi3up,z00,z00+Lz,xyz_star(3))/a
        call grid_profile(a*(xi3  -xi3star),g3,g3der1,g3der2)
        call grid_profile(a*(xi3lo-xi3star),g3lo)
        call grid_profile(a*(xi3up-xi3star),g3up)

        z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
        zprim =    Lz*(g3der1*a   )/(g3up-g3lo)
        zprim2=    Lz*(g3der2*a**2)/(g3up-g3lo)
!       
        dz_1=1./zprim
        dz_tilde=-zprim2/zprim**2
!
      endif
!
!  determine global minimum and maximum of grid spacing in any direction
!
      if (nxgrid <= 1) then
        dxmin_x = dx
        dxmax_x = dx
      endif
      !
      if (nygrid <= 1) then
        dxmin_y = dy
        dxmax_y = dy
      endif
      !
      if (nzgrid <= 1) then
        dxmin_z = dz
        dxmax_z = dz
      endif

      dxmin = minval( (/dxmin_x, dxmin_y, dxmin_z, huge(dx)/), &
                MASK=((/nxgrid, nygrid, nzgrid, 2/) > 1) )

      dxmax = maxval( (/dxmax_x, dxmax_y, dxmax_z, epsilon(dx)/), &
                MASK=((/nxgrid, nygrid, nzgrid, 2/) > 1) )
!
    endsubroutine construct_grid
!***********************************************************************
    subroutine pencil_interdep_grid(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  15-nov-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_rcyl_mn1)) lpencil_in(i_rcyl_mn)=.true.
      if (lpencil_in(i_evr)) lpencil_in(i_r_mn)=.true.
      if (lpencil_in(i_evth)) then
         lpencil_in(i_pomx)=.true.
         lpencil_in(i_pomy)=.true.
         lpencil_in(i_rcyl_mn)=.true.
         lpencil_in(i_r_mn1)=.true.
      endif
      if (  lpencil_in(i_pomx) &
       .or. lpencil_in(i_pomy) &
       .or. lpencil_in(i_phix) &
       .or. lpencil_in(i_phiy)) then
        lpencil_in(i_rcyl_mn1)=.true.
      endif
!
      if (lpencil_in(i_rr).or.lpencil_in(i_evr)) then
        lpencil_in(i_x_mn)=.true.
        lpencil_in(i_y_mn)=.true.
        lpencil_in(i_z_mn)=.true.
        if (lpencil_in(i_evr)) lpencil_in(i_r_mn)=.true.
      endif

!
    endsubroutine pencil_interdep_grid
!***********************************************************************
    subroutine calc_pencils_grid(f,p)
!
!  Calculate Grid/geometry related pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   15-nov-06/tony: coded
!   27-aug-07/wlad: generalized for cyl. and sph. coordinates
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
!coordinates vectors
      if (lpencil(i_x_mn))     p%x_mn    = x(l1:l2)
      if (lpencil(i_y_mn))     p%y_mn    = spread(y(m),1,nx)
      if (lpencil(i_z_mn))     p%z_mn    = spread(z(n),1,nx)
!spherical distance  
      if (lpencil(i_r_mn))     p%r_mn    = sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
!cylindrical distance (pomega)
      if (lpencil(i_rcyl_mn))  p%rcyl_mn = sqrt(x(l1:l2)**2+y(m)**2)
!azimuthal angle (phi)
      if (lpencil(i_phi_mn))   p%phi_mn  = atan2(y(m),x(l1:l2))
!inverse cylindrical distance 1/pomega
      if (lpencil(i_rcyl_mn1)) p%rcyl_mn1=1./max(p%rcyl_mn,tini)
!inverse spherical distance 1/r
      if (lpencil(i_r_mn1))    p%r_mn1   =1./max(p%r_mn,tini)
!pomega unit vectors: pomx=cos(phi) and pomy=sin(phi) where phi=azimuthal angle
      if (lpencil(i_pomx))     p%pomx    = x(l1:l2)*p%rcyl_mn1
      if (lpencil(i_pomy))     p%pomy    = y(  m  )*p%rcyl_mn1
!phi unit vectors
      if (lpencil(i_phix))     p%phix    =-y(  m  )*p%rcyl_mn1
      if (lpencil(i_phiy))     p%phiy    = x(l1:l2)*p%rcyl_mn1
!
!  set position vector
!
      if (lpencil(i_rr)) then
        p%rr(:,1)=p%x_mn
        p%rr(:,2)=p%y_mn
        p%rr(:,3)=p%z_mn
      endif
!
!  evr is the radial unit vector
!
      if (lpencil(i_evr)) then
        p%evr(:,1) = p%x_mn
        p%evr(:,2) = p%y_mn
        p%evr(:,3) = p%z_mn
        p%evr = p%evr / spread(p%r_mn+tini,2,3)
      endif
!
!  evth is the latitudinal unit vector
!
      if (lpencil(i_evth)) then
        p%evth(:,1) = -z(n)*p%r_mn1*p%pomx
        p%evth(:,2) = -z(n)*p%r_mn1*p%pomy
        p%evth(:,3) = p%rcyl_mn*p%r_mn1
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_grid
!***********************************************************************
    subroutine grid_profile_point(xi,g,gder1,gder2)
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
      real              :: g
      real, optional    :: gder1,gder2
!
      intent(in)  :: xi
      intent(out) :: g,gder1,gder2
!
      g=xi
      if (present(gder1)) gder1=1.0
      if (present(gder2)) gder2=0.0
!
    endsubroutine grid_profile_point
!***********************************************************************
    subroutine grid_profile_1d(xi,g,gder1,gder2)
!
!  Same as grid_profile_1 for 1d arrays as arguments
!
!  25-jun-04/tobi+wolf: coded
!
      real, dimension(:)                    :: xi
      real, dimension(size(xi,1))           :: g
      real, dimension(size(xi,1)), optional :: gder1,gder2
!
      intent(in)  :: xi
      intent(out) :: g, gder1,gder2
!
      g=xi
      if (present(gder1)) gder1=1.0
      if (present(gder2)) gder2=0.0
!
    endsubroutine grid_profile_1d
!***********************************************************************
    function find_star(xi_lo,xi_up,x_lo,x_up,x_star) result (xi_star)
!
!  Finds the xi that corresponds to the inflection point of the grid-function
!  by means of a newton-raphson root-finding algorithm.
!
!  25-jun-04/tobi+wolf: coded
!
      real, intent(in) :: xi_lo,xi_up,x_lo,x_up,x_star
!
      real :: xi_star,dxi,tol
      real :: g_lo,gder_lo
      real :: g_up,gder_up
      real :: f   ,fder
      integer, parameter :: maxit=1000
      logical :: lreturn
      integer :: it
!
      if (xi_lo>=xi_up) &
           call fatal_error('find_star','xi1 >= xi2 -- this should not happen')
!
      tol=epsi*(xi_up-xi_lo)
      xi_star= (xi_up+xi_lo)/2
!
      lreturn=.false.
!
      do it=1,maxit
!
        call grid_profile(xi_lo-xi_star,g_lo,gder_lo)
        call grid_profile(xi_up-xi_star,g_up,gder_up)
!
        f   =-(x_up-x_star)*g_lo   +(x_lo-x_star)*g_up
        fder= (x_up-x_star)*gder_lo-(x_lo-x_star)*gder_up
!
        dxi=f/fder
        xi_star=xi_star-dxi
!
        if (lreturn) return
!
        if (abs(dxi)<tol) lreturn=.true.
!
      enddo
!
      call fatal_error('find_star','maximum number of iterations exceeded')
!
    endfunction find_star
!***********************************************************************
endmodule Grid
