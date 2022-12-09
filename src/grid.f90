! $Id$
!
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
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  private
!
  public :: construct_grid
  public :: pencil_criteria_grid
  public :: pencil_interdep_grid
  public :: calc_pencils_grid
  public :: initialize_grid
  public :: get_grid_mn
  public :: box_vol
  public :: save_grid
  public :: coords_aux
  public :: real_to_index, inverse_grid
  public :: grid_bound_data
  public :: set_coorsys_dimmask
  public :: construct_serial_arrays
  public :: coarsegrid_interp
!
  public :: find_star
  public :: calc_bound_coeffs
  public :: grid_profile

  interface grid_profile
    module procedure grid_profile_0D
    module procedure grid_profile_1D
  endinterface
!
  interface calc_pencils_grid
    module procedure calc_pencils_grid_pencpar
    module procedure calc_pencils_grid_std
  endinterface calc_pencils_grid
!
  integer, parameter :: BOT=1, TOP=2
!
  contains
!***********************************************************************
    !subroutine construct_grid(x,y,z,dx,dy,dz,x00,y00,z00)
    subroutine construct_grid(x,y,z,dx,dy,dz)
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
!  25-jun-04/tobi+wolf: coded
!   3-mar-15/MR: calculation of d[xyz]2_bound added: contain twice the distances of
!                three neighbouring points from the boundary point
!   9-mar-17/MR: removed unneeded use of optional parameters in calls to grid_profile
!   6-mar-18/MR: moved call construct_serial_arrays here from initialize_grid
!
      use Cdata, only: xprim, yprim, zprim, dx_1, dx_tilde, dy_1, dy_tilde, dz_1, dz_tilde

      real, dimension(mx), intent(out) :: x
      real, dimension(my), intent(out) :: y
      real, dimension(mz), intent(out) :: z
      real, intent(out) :: dx,dy,dz
!
      real :: x00,y00,z00
      real :: xi1lo,xi1up,g1lo,g1up
      real :: xi2lo,xi2up,g2lo,g2up
      real :: xi3lo,xi3up,g3lo,g3up
      real, dimension(3,2) :: xi_step
      real, dimension(3,3) :: dxyz_step
      real :: xi1star,xi2star,xi3star,bound_prim1,bound_prim2
!
      real, dimension(mx) :: g1,g1der1,g1der2,xi1,xprim2
      real, dimension(my) :: g2,g2der1,g2der2,xi2,yprim2
      real, dimension(mz) :: g3,g3der1,g3der2,xi3,zprim2
!
      real, dimension(0:nprocx) :: xi1proc, g1proc
      real, dimension(0:nprocy) :: xi2proc, g2proc
      real, dimension(0:nprocz) :: xi3proc, g3proc
!
      real :: a,b,c
      integer :: i
!
      lequidist=(grid_func=='linear')
!
!  Abbreviations
!
      x0 = xyz0(1)
      y0 = xyz0(2)
      z0 = xyz0(3)
      Lx = Lxyz(1)
      Ly = Lxyz(2)
      Lz = Lxyz(3)
!
!  Set the lower boundary and the grid size.
!
      x00 = x0
      y00 = y0
      z00 = z0
!
      dx = Lx / merge(nxgrid, max(nxgrid-1,1), lperi(1))
      dy = Ly / merge(nygrid, max(nygrid-1,1), (lperi(2).or.lpole(2)))
      dz = Lz / merge(nzgrid, max(nzgrid-1,1), lperi(3))
!
!  Shift the lower boundary if requested, but only for periodic directions.
!
      if (lshift_origin(1)) x00 = x0 + 0.5 * dx
      if (lshift_origin(2)) y00 = y0 + 0.5 * dy
      if (lshift_origin(3)) z00 = z0 + 0.5 * dz
!
!  Shift the lower boundary if requested, but only for periodic directions.
!  Contrary to the upper case (lshift_origin)
!
      if (lshift_origin_lower(1)) x00 = x0 - 0.5 * dx
      if (lshift_origin_lower(2)) y00 = y0 - 0.5 * dy
      if (lshift_origin_lower(3)) z00 = z0 - 0.5 * dz
!
!  Produce index arrays xi1, xi2, and xi3:
!    xi = 0, 1, 2, ..., N-1 for non-periodic grid
!    xi = 0.5, 1.5, 2.5, ..., N-0.5 for periodic grid
!
      do i=1,mx; xi1(i)=i-nghost-1+ipx*nx; enddo
      do i=1,my; xi2(i)=i-nghost-1+ipy*ny; enddo
      do i=1,mz; xi3(i)=i-nghost-1+ipz*nz; enddo
!
      if (lperi(1)) xi1 = xi1 + 0.5
      if (lperi(2).or.lpole(2)) xi2 = xi2 + 0.5
      if (lperi(3)) xi3 = xi3 + 0.5
!
!  Produce index arrays for processor boundaries, which are needed for
!  particle migration (see redist_particles_bounds). The select cases
!  should use these arrays to set g{2,3}proc using the grid function.
!
      xi1proc = (/ (real(i * nx), i = 0, nprocx) /)
      xi1p: if (lactive_dimension(1) .and. .not. lperi(1)) then
        xi1proc = xi1proc - 0.5
        xi1proc(0) = 0.0
        xi1proc(nprocx) = real(nxgrid - 1)
      endif xi1p
!
      xi2proc = (/ (real(i * ny), i = 0, nprocy) /)
      xi2p: if (lactive_dimension(2) .and. .not. (lperi(2) .or. lpole(2))) then
        xi2proc = xi2proc - 0.5
        xi2proc(0) = 0.0
        xi2proc(nprocy) = real(nygrid - 1)
      endif xi2p
!
      xi3proc = (/ (real(i * nz), i = 0, nprocz) /)
      xi3p: if (lactive_dimension(3) .and. .not. lperi(3)) then
        xi3proc = xi3proc - 0.5
        xi3proc(0) = 0.0
        xi3proc(nprocz) = real(nzgrid - 1)
      endif xi3p
!
!  The following is correct for periodic and non-periodic case
!    Periodic: x(xi=0) = x0 and x(xi=N) = x1
!    Non-periodic: x(xi=0) = x0 and x(xi=N-1) = x1
!
      xi1lo=0.; xi1up=nxgrid-merge(0.,1.,lperi(1))
      xi2lo=0.; xi2up=nygrid-merge(0.,1.,(lperi(2).or.lpole(2)))
      xi3lo=0.; xi3up=nzgrid-merge(0.,1.,lperi(3))
!
      if (lpole(2)) xyz_star(2) = pi/2.
!
!  Construct nonequidistant grid
!
!  x coordinate
!
      if (nxgrid==1) then
        x = x00 + 0.5 * dx
        ! hopefully, we will only ever multiply by the following quantities:
        ! [PABourdin] We should find a way to have valid grid functions as:
        ! xprim = dx  ;  dx_1 = 1/dx
        ! together with a degenerated-case flag: ldegenerated_x = .true.
        xprim = 0.
        xprim2 = 0.
        dx_1 = 0.
        dx_tilde = 0.
        g1proc(0) = x00
        g1proc(1) = x00 + Lx
      else
!
        select case (grid_func(1))
!
        case ('linear','sinh')
          a=coeff_grid(1)*dx
          xi1star=find_star(a*xi1lo,a*xi1up,x00,x00+Lx,xyz_star(1),grid_func(1))/a
          call grid_profile(a*(xi1  -xi1star),grid_func(1),g1,g1der1,g1der2)
          call grid_profile(a*(xi1lo-xi1star),grid_func(1),g1lo)
          call grid_profile(a*(xi1up-xi1star),grid_func(1),g1up)
          call grid_profile(a * (xi1proc - xi1star), grid_func(1), g1proc)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
          g1proc = x00 + Lx * (g1proc - g1lo) / (g1up - g1lo)
!
        case ('cos','tanh')
          ! Approximately equidistant at the boundaries, linear in the middle
          a=pi*dx/trans_width(1)
          b=dxi_fact(1)
          xi1star = xi1lo + (xi1up-xi1lo) / (1.0 + (Lx/(xyz_star(1)-x00) - 1.0) / b)
          call grid_profile(a*(xi1  -xi1star),grid_func(1),g1,g1der1,g1der2,param=b)
          call grid_profile(a*(xi1lo-xi1star),grid_func(1),g1lo,param=b)
          call grid_profile(a*(xi1up-xi1star),grid_func(1),g1up,param=b)
          call grid_profile(a * (xi1proc - xi1star), grid_func(1), g1proc, param=b)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
          g1proc = x00 + Lx * (g1proc - g1lo) / (g1up - g1lo)
!
        case ('arsinh')
          ! Approximately linear in infinity, linear in the middle
          a=pi*dx/trans_width(1)
          b=dxi_fact(1)
          xi1star = xi1lo + (xi1up-xi1lo) / (1.0 + (Lx/(xyz_star(1)-x00) - 1.0) / b)
          call grid_profile(a*(xi1  -xi1star),grid_func(1),g1,g1der1,g1der2,param=b)
          call grid_profile(a*(xi1lo-xi1star),grid_func(1),g1lo,param=b)
          call grid_profile(a*(xi1up-xi1star),grid_func(1),g1up,param=b)
          call grid_profile(a * (xi1proc - xi1star), grid_func(1), g1proc, param=b)
!
          ! Slope should be 1 at the minimum:
          b = minval (g1der1)
          if (b < 1.0) then
            g1der1 = g1der1 + 1.0 - b
            g1     = g1   + (1.0 - b) * a * (xi1   - xi1star)
            g1lo   = g1lo + (1.0 - b) * a * (xi1lo - xi1star)
            g1up   = g1up + (1.0 - b) * a * (xi1up - xi1star)
            g1proc = g1proc + (1.0 - b) * a * (xi1proc - xi1star)
          endif
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
          g1proc = x00 + Lx * (g1proc - g1lo) / (g1up - g1lo)
!
        case ('step-linear')
!
          xi_step(1,1)=xi_step_frac(1,1)*(nxgrid-1.0)
          xi_step(1,2)=xi_step_frac(1,2)*(nxgrid-1.0)
          dxyz_step(1,1)=(xyz_step(1,1)-x00)/(xi_step(1,1)-0.0)
          dxyz_step(1,2)=(xyz_step(1,2)-xyz_step(1,1))/ &
                                (xi_step(1,2)-xi_step(1,1))
          dxyz_step(1,3)=(x00+Lx-xyz_step(1,2))/(nxgrid-1.0-xi_step(1,2))
!
          call grid_profile(xi1,grid_func(1),g1,g1der1,g1der2, &
           dxyz=dxyz_step(1,:),xistep=xi_step(1,:),delta=xi_step_width(1,:))
          call grid_profile(xi1lo,grid_func(1),g1lo, &
           dxyz=dxyz_step(1,:),xistep=xi_step(1,:),delta=xi_step_width(1,:))
          call grid_profile(xi1proc, grid_func(1), g1proc, dxyz=dxyz_step(1,:), xistep=xi_step(1,:), delta=xi_step_width(1,:))
!
          x     = x00 + g1-g1lo
          xprim = g1der1
          xprim2= g1der2
          g1proc = x00 + g1proc - g1lo
!
        case ('duct')
          a = pi/(max(nxgrid-1,1))
          call grid_profile(a*xi1  -pi/2,grid_func(1),g1,g1der1,g1der2)
          call grid_profile(a*xi1lo-pi/2,grid_func(1),g1lo)
          call grid_profile(a*xi1up-pi/2,grid_func(1),g1up)
          call grid_profile(a * xi1proc - 0.5 * pi, grid_func(1), g1proc)
!
          x     =x00+Lx*(g1-g1lo)/2
          xprim =    Lx*(g1der1*a   )/2
          xprim2=    Lx*(g1der2*a**2)/2
          g1proc = x00 + 0.5 * Lx * (g1proc - g1lo)
!
          if (lfirst_proc_x) then
            bound_prim1=x(l1+1)-x(l1)
            do i=1,nghost
              x(l1-i)=x(l1)-i*bound_prim1
              xprim(1:l1)=bound_prim1
            enddo
          endif
          if (llast_proc_x) then
            bound_prim2=x(l2)-x(l2-1)
            do i=1,nghost
              x(l2+i)=x(l2)+i*bound_prim2
              xprim(l2:mx)=bound_prim2
            enddo
          endif
!
! half-duct profile : like the duct above but the grid are closely spaced
! at the outer boundary.
!
        case ('half-duct')
          a =-pi/(2.*max(nxgrid-1,1))
          call grid_profile(pi/2.+a*xi1,grid_func(1),g1,g1der1,g1der2)
          call grid_profile(pi/2.+a*xi1lo,grid_func(1),g1lo)
          call grid_profile(pi/2.+a*xi1up,grid_func(1),g1up)
          call grid_profile(0.5 * pi + a * xi1proc, grid_func(1), g1proc)
!
          x     =x00+Lx*g1
          xprim =    Lx*(g1der1*a   )
          xprim2=    Lx*(g1der2*a**2)
          g1proc = x00 + Lx * g1proc
!
          if (lfirst_proc_x) then
            bound_prim1=x(l1+1)-x(l1)
            do i=1,nghost
              x(l1-i)=x(l1)-i*bound_prim1
              xprim(1:l1)=bound_prim1
            enddo
          endif
          if (llast_proc_x) then
            bound_prim2=x(l2)-x(l2-1)
            do i=1,nghost
              x(l2+i)=x(l2)+i*bound_prim2
              xprim(l2:mx)=bound_prim2
            enddo
          endif
!
        case ('log')
          ! Grid distance increases logarithmically
          ! .i.e., grid spacing increases linearly
          ! d[log x] = const ==> dx = const*x
          a= log(xyz1(1)/xyz0(1))/(xi1up-xi1lo)
          b= .5*(xi1up+xi1lo-log(xyz1(1)*xyz0(1))/a)
!
          call grid_profile(a*(xi1-b)  ,grid_func(1),g1,g1der1,g1der2)
          call grid_profile(a*(xi1lo-b),grid_func(1),g1lo)
          call grid_profile(a*(xi1up-b),grid_func(1),g1up)
          call grid_profile(a * (xi1proc - b), grid_func(1), g1proc)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
          g1proc = x00 + Lx * (g1proc - g1lo) / (g1up - g1lo)
!
        case ('power-law')
          ! Grid distance increases as a power-law set by c=coeff_grid(1)
          ! i.e., grid distance increases as an inverse power-law
          ! d[x^c] = const ==> dx = const/x^(c-1)
          if (lroot) then
            print*,"Constructing power-law grid. It will "
            print*,"generate a grid obeying the constrain "
            print*,""
            print*,"  d[x^c] = const ==> dx = const/x^(c-1)"
            print*,""
            print*,"where c is the value of coeff_grid."
          endif
!
          if (coeff_grid(1) == 0.) &
               call fatal_error('construct_grid:', 'Cannot create '//&
               'a grid for a power-law of zero. Use "log" instead.')
          c= coeff_grid(1)
          a= (xyz1(1)**c-xyz0(1)**c)/(xi1up-xi1lo)
          b= .5*(xi1up+xi1lo-(xyz1(1)**c+xyz0(1)**c)/a)
!
          call grid_profile(a*(xi1-b)  ,grid_func(1),g1,g1der1,g1der2,param=c)
          call grid_profile(a*(xi1lo-b),grid_func(1),g1lo,param=c)
          call grid_profile(a*(xi1up-b),grid_func(1),g1up,param=c)
          call grid_profile(a * (xi1proc - b), grid_func(1), g1proc, param=c)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
          g1proc = x00 + Lx * (g1proc - g1lo) / (g1up - g1lo)
!
        case ('squared')
          ! Grid distance increases linearily
          a = real(max(nxgrid, 1))
          !b=-max(nxgrid,1)/10  ! 23-nov-20/ccyang: is it correct to use integer arithmetic like this?
          b = 0.0
          call grid_profile(a*(xi1  -b),grid_func(1),g1,g1der1,g1der2)
          call grid_profile(a*(xi1lo-b),grid_func(1),g1lo)
          call grid_profile(a*(xi1up-b),grid_func(1),g1up)
          call grid_profile(a * (xi1proc - b), grid_func(1), g1proc)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
          g1proc = x00 + Lx * (g1proc - g1lo) / (g1up - g1lo)
!
          if (lfirst_proc_x) then
            bound_prim1=x(l1+1)-x(l1)
            do i=1,nghost
              x(l1-i)=x(l1)-i*bound_prim1
              xprim(1:l1)=bound_prim1
            enddo
          endif
          if (llast_proc_x) then
            bound_prim2=x(l2)-x(l2-1)
            do i=1,nghost
              x(l2+i)=x(l2)+i*bound_prim2
              xprim(l2:mx)=bound_prim2
            enddo
          endif
!
        case ('frozensphere')
          ! Just like sinh, except set dx constant below a certain radius, and
          ! constant for top ghost points.
          a=coeff_grid(1)*dx
          xi1star=find_star(a*xi1lo,a*xi1up,x00,x00+Lx,0.8*xyz_star(1),grid_func(1))/a
          call grid_profile(a*(xi1  -xi1star),grid_func(1),g1,g1der1,g1der2)
          call grid_profile(a*(xi1lo-xi1star),grid_func(1),g1lo)
          call grid_profile(a*(xi1up-xi1star),grid_func(1),g1up)
          call grid_profile(a * (xi1proc - xi1star), grid_func(1), g1proc)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
          g1proc = x00 + Lx * (g1proc - g1lo) / (g1up - g1lo)
!
          if (llast_proc_x) then
            bound_prim2=x(l2-2)-x(l2-3)
            do i=1,nghost+2
              x(l2-2+i)=x(l2-2)+i*bound_prim2
              xprim(l2-2:mx)=bound_prim2
            enddo
          endif
!
        case default
          call fatal_error('construct_grid', &
                           'No such x grid function - '//trim(grid_func(1)))
        endselect
!
!  Fitting a polynomial function to the grid function to set the xprim2 zero at the boundary,
!  this makes sure that the second derivative is zero, if you choose a2 for example.
!  Still in the testing phase, that why it is commented out.
!
!        if ( .not.lequidist(1) ) then
!          if (xprim2(l1) .ne. 0) then
!            print*, 'use polynomia extension in x bottom'
!            if (lfirst_proc_x) then
!              a3=xprim2(l1+1)/6.
!              a1=xprim(l1+1) - xprim2(l1+1)/2.
!             a0=x(l1+1) - xprim(l1+1) + xprim2(l1+1)/3.
!             x(l1)=a0
!             xprim(l1) = a1
!             xprim2(l1) = 0.
!           endif
!          endif
!          if (xprim2(l2) .ne. 0) then
!            print*, 'use polynomia extension in x top'
!            if (llast_proc_x) then
!              a3=-xprim2(l2-1)/6.
!              a1=xprim(l2-1) + xprim2(l2-1)/2.
!              a0=x(l2-1) +xprim(l2-1) + xprim2(l2-1)/3.
!              x(l2)=a0
!              xprim(l2) = a1
!              xprim2(l2) = 0
!            endif
!          endif
!        endif
!
        dx2=xprim**2
        dx_1=1./xprim
        dx_tilde=-xprim2/dx2
!
      endif
!
!  y coordinate
!
      if (nygrid==1) then
        y = y00 + 0.5 * dy
        ! hopefully, we will only ever multiply by the following quantities:
        ! [PABourdin] We should find a way to have valid grid functions as:
        ! xprim = dx  ;  dx_1 = 1/dx
        ! together with a degenerated-case flag: ldegenerated_x = .true.
        yprim = 0.
        yprim2 = 0.
        dy_1 = 0.
        dy_tilde = 0.
        g2proc(0) = y00
        g2proc(1) = y00 + Ly
      else
!
        select case (grid_func(2))
!
        case ('linear','sinh')
!
          a=coeff_grid(2)*dy
          xi2star=find_star(a*xi2lo,a*xi2up,y00,y00+Ly,xyz_star(2),grid_func(2))/a
          call grid_profile(a*(xi2  -xi2star),grid_func(2),g2,g2der1,g2der2)
          call grid_profile(a*(xi2lo-xi2star),grid_func(2),g2lo)
          call grid_profile(a*(xi2up-xi2star),grid_func(2),g2up)
          call grid_profile(a * (xi2proc - xi2star), grid_func(2), g2proc)
!
          y     =y00+Ly*(g2  -  g2lo)/(g2up-g2lo)
          yprim =    Ly*(g2der1*a   )/(g2up-g2lo)
          yprim2=    Ly*(g2der2*a**2)/(g2up-g2lo)
          g2proc = y00 + Ly * (g2proc - g2lo) / (g2up - g2lo)
!
        case ('cos','tanh')
          ! Approximately equidistant at the boundaries, linear in the middle
          a=pi*dy/trans_width(2)
          b=dxi_fact(2)
          xi2star = xi2lo + (xi2up-xi2lo) / (1.0 + (Ly/(xyz_star(2)-y00) - 1.0) / b)
          call grid_profile(a*(xi2  -xi2star),grid_func(2),g2,g2der1,g2der2,param=b)
          call grid_profile(a*(xi2lo-xi2star),grid_func(2),g2lo,param=b)
          call grid_profile(a*(xi2up-xi2star),grid_func(2),g2up,param=b)
          call grid_profile(a * (xi2proc - xi2star), grid_func(2), g2proc, param=b)
!
          y     =y00+Ly*(g2  -  g2lo)/(g2up-g2lo)
          yprim =    Ly*(g2der1*a   )/(g2up-g2lo)
          yprim2=    Ly*(g2der2*a**2)/(g2up-g2lo)
          g2proc = y00 + Ly * (g2proc - g2lo) / (g2up - g2lo)
!
        case ('arsinh')
          ! Approximately linear in infinity, linear in the middle
          a=pi*dy/trans_width(2)
          b=dxi_fact(2)
          xi2star = xi2lo + (xi2up-xi2lo) / (1.0 + (Ly/(xyz_star(2)-y00) - 1.0) / b)
          call grid_profile(a*(xi2  -xi2star),grid_func(2),g2,g2der1,g2der2,param=b)
          call grid_profile(a*(xi2lo-xi2star),grid_func(2),g2lo,param=b)
          call grid_profile(a*(xi2up-xi2star),grid_func(2),g2up,param=b)
          call grid_profile(a * (xi2proc - xi2star), grid_func(2), g2proc, param=b)
!
          ! Slope should be 1 at the minimum:
          b = minval (g2der1)
          if (b < 1.0) then
            g2der1 = g2der1 + 1.0 - b
            g2     = g2   + (1.0 - b) * a * (xi2   - xi2star)
            g2lo   = g2lo + (1.0 - b) * a * (xi2lo - xi2star)
            g2up   = g2up + (1.0 - b) * a * (xi2up - xi2star)
            g2proc = g2proc + (1.0 - b) * a * (xi2proc - xi2star)
          endif
!
          y     =y00+Ly*(g2  -  g2lo)/(g2up-g2lo)
          yprim =    Ly*(g2der1*a   )/(g2up-g2lo)
          yprim2=    Ly*(g2der2*a**2)/(g2up-g2lo)
          g2proc = y00 + Ly * (g2proc - g2lo) / (g2up - g2lo)
!
        case ('duct')
!
          a = pi/max(nygrid-1, 1)
          call grid_profile(a*xi2  -pi/2,grid_func(2),g2,g2der1,g2der2)
          call grid_profile(a*xi2lo-pi/2,grid_func(2),g2lo)
          call grid_profile(a*xi2up-pi/2,grid_func(2),g2up)
          call grid_profile(a * xi2proc - 0.5 * pi, grid_func(2), g2proc)
!
          y     =y00+Ly*(g2-g2lo)/2
          yprim =    Ly*(g2der1*a   )/2
          yprim2=    Ly*(g2der2*a**2)/2
          g2proc = y00 + 0.5 * Ly * (g2proc - g2lo)
!
          if (lfirst_proc_y) then
            bound_prim1=y(m1+1)-y(m1)
            do i=1,nghost
              y(m1-i)=y(m1)-i*bound_prim1
              yprim(1:m1)=bound_prim1
            enddo
          endif
          if (llast_proc_y) then
            bound_prim2=y(m2)-y(m2-1)
            do i=1,nghost
              y(m2+i)=y(m2)+i*bound_prim2
              yprim(m2:my)=bound_prim2
            enddo
          endif
!
        case ('step-linear')
!
          xi_step(2,1)=xi_step_frac(2,1)*merge(nygrid+0.,nygrid-1.,lpole(2))
          xi_step(2,2)=xi_step_frac(2,2)*merge(nygrid+0.,nygrid-1.,lpole(2))
          dxyz_step(2,1)=(xyz_step(2,1)-y00)/(xi_step(2,1)-0.0)
          dxyz_step(2,2)=(xyz_step(2,2)-xyz_step(2,1))/ &
              (xi_step(2,2)-xi_step(2,1))
          dxyz_step(2,3)=(y00+Ly-xyz_step(2,2))/&
              (merge(nygrid+0.,nygrid-1.,lpole(2))-xi_step(2,2))
!
          if (lpole(2)) then 
            do i=1,my
              xi2(i)=0.5-nghost+0.5+(ipy*ny+i-1)*&
                (nygrid+2.*(nghost-0.5)-1)/(nygrid+2.*nghost-1)
            enddo
            xi2proc = (/ (real(1 - nghost) + (real(nghost + i * ny) - 0.5) * real(mygrid - 2) &
                                                                           / real(mygrid - 1), i = 0, nprocy) /)
          endif
!
          call grid_profile(xi2,grid_func(2),g2,g2der1,g2der2,dxyz= &
              dxyz_step(2,:),xistep=xi_step(2,:),delta=xi_step_width(2,:))
          call grid_profile(xi2lo,grid_func(2),g2lo,dxyz= &
              dxyz_step(2,:),xistep=xi_step(2,:),delta=xi_step_width(2,:))
          call grid_profile(xi2proc, grid_func(2), g2proc, dxyz=dxyz_step(2,:), xistep=xi_step(2,:), delta=xi_step_width(2,:))
!
          if (lpole(2)) then
            call grid_profile(xi2up,grid_func(2),g2up,dxyz= &
                dxyz_step(2,:),xistep=xi_step(2,:),delta=xi_step_width(2,:))
            y = y00 + g2 -0.5*(g2lo+g2up-pi)
            g2proc = y00 + g2proc - 0.5 * (g2lo + g2up - pi)
          else
!if(iproc<2) print*, 'iproc, y00, g2-g2lo=', iproc, y00, g2-g2lo
            y = y00 + g2-g2lo
            g2proc = y00 + g2proc - g2lo
          endif
!
          yprim = g2der1
          yprim2= g2der2
!
        case ('trough')
          ! Grid distance is almost equidistant at boundaries and then decreases
          ! at the middle
          a=0.05
          b=40
          c=160-40
!
          call grid_profile(xi2, grid_func(2),g2,g2der1,g2der2,param=a,xistep=(/b,c/),delta=(/0.5,0.1/))
          call grid_profile(xi2lo,grid_func(2),g2lo,param=a,xistep=(/b,c/),delta=(/0.5,0.1/))
          call grid_profile(xi2up,grid_func(2),g2up,param=a,xistep=(/b,c/),delta=(/0.5,0.1/))
          call grid_profile(xi2proc, grid_func(2), g2proc, param=a, xistep=(/b,c/), delta=(/0.5,0.1/))
!
          y     =y00+Ly*(g2  -  g2lo)/(g2up-g2lo)
          yprim =    Ly*g2der1/(g2up-g2lo)
          yprim2=    Ly*g2der2/(g2up-g2lo)
          g2proc = y00 + Ly * (g2proc - g2lo) / (g2up - g2lo)
!
        case default
          call fatal_error('construct_grid', &
                           'No such y grid function - '//trim(grid_func(2)))
!
        endselect
!
! Added parts for spherical coordinates and cylindrical coordinates.
! From now on dy = d\theta but dy_1 = 1/rd\theta and similarly for \phi.
! corresponding r and rsin\theta factors for equ.f90 (where CFL timesteps
! are estimated) are removed.
!
        if (lpole(2)) then                       !apply grid symmetry across the poles 
          if (lfirst_proc_y) then
            y(     1:nghost) = -y(     m1i:m1:-1)
            yprim( 1:nghost) =  yprim( m1i:m1:-1)
            yprim2(1:nghost) = -yprim2(m1i:m1:-1)
          endif
          if (llast_proc_y) then
            y(     m2+1:m2+nghost) = 2.*pi-y(     m2:m2i:-1)
            yprim( m2+1:m2+nghost) =       yprim( m2:m2i:-1)
            yprim2(m2+1:m2+nghost) =      -yprim2(m2:m2i:-1)
          endif
        endif
!
        dy2=yprim**2
        dy_1=1./yprim
        dy_tilde=-yprim2/dy2
!
      endif
!
!  z coordinate
!
      if (nzgrid==1) then
        z = z00 + 0.5 * dz
        ! hopefully, we will only ever multiply by the following quantities:
        ! [PABourdin] We should find a way to have valid grid functions as:
        ! xprim = dx  ;  dx_1 = 1/dx
        ! together with a degenerated-case flag: ldegenerated_x = .true.
        zprim = 0.
        zprim2 = 0.
        dz_1 = 0.
        dz_tilde = 0.
        g3proc(0) = z00
        g3proc(1) = z00 + Lz
      else
!
        select case (grid_func(3))
!
        case ('linear','sinh')
!
          a=coeff_grid(3)*dz
          xi3star=find_star(a*xi3lo,a*xi3up,z00,z00+Lz,xyz_star(3),grid_func(3))/a
          call grid_profile(a*(xi3  -xi3star),grid_func(3),g3,g3der1,g3der2)
          call grid_profile(a*(xi3lo-xi3star),grid_func(3),g3lo)
          call grid_profile(a*(xi3up-xi3star),grid_func(3),g3up)
          call grid_profile(a * (xi3proc - xi3star), grid_func(3), g3proc)
!
          z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
          zprim =    Lz*(g3der1*a   )/(g3up-g3lo)
          zprim2=    Lz*(g3der2*a**2)/(g3up-g3lo)
          g3proc = z00 + Lz * (g3proc - g3lo) / (g3up - g3lo)
!
        case ('cos','tanh')
          ! Approximately equidistant at the boundaries, linear in the middle
          a=pi*dz/trans_width(3)
          b=dxi_fact(3)
          xi3star = xi3lo + (xi3up-xi3lo) / (1.0 + (Lz/(xyz_star(3)-z00) - 1.0) / b)
          call grid_profile(a*(xi3  -xi3star),grid_func(3),g3,g3der1,g3der2,param=b)
          call grid_profile(a*(xi3lo-xi3star),grid_func(3),g3lo,param=b)
          call grid_profile(a*(xi3up-xi3star),grid_func(3),g3up,param=b)
          call grid_profile(a * (xi3proc - xi3star), grid_func(3), g3proc, param=b)
!
          z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
          zprim =    Lz*(g3der1*a   )/(g3up-g3lo)
          zprim2=    Lz*(g3der2*a**2)/(g3up-g3lo)
          g3proc = z00 + Lz * (g3proc - g3lo) / (g3up - g3lo)
!
        case ('arsinh')
          ! Approximately linear in infinity, linear in the middle
          a=pi*dz/trans_width(3)
          b=dxi_fact(3)
          xi3star = xi3lo + (xi3up-xi3lo) / (1.0 + (Lz/(xyz_star(3)-z00) - 1.0) / b)
          call grid_profile(a*(xi3  -xi3star),grid_func(3),g3,g3der1,g3der2,param=b)
          call grid_profile(a*(xi3lo-xi3star),grid_func(3),g3lo,param=b)
          call grid_profile(a*(xi3up-xi3star),grid_func(3),g3up,param=b)
          call grid_profile(a * (xi3proc - xi3star), grid_func(3), g3proc, param=b)
!
          ! Slope should be 1 at the minimum:
          b = minval (g3der1)
          if ((grid_func(3) == 'arsinh') .and. (b < 1.0)) then
            g3der1 = g3der1 + 1.0 - b
            g3     = g3   + (1.0 - b) * a * (xi3   - xi3star)
            g3lo   = g3lo + (1.0 - b) * a * (xi3lo - xi3star)
            g3up   = g3up + (1.0 - b) * a * (xi3up - xi3star)
            g3proc = g3proc + (1.0 - b) * a * (xi3proc - xi3star)
          endif
!
          z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
          zprim =    Lz*(g3der1*a   )/(g3up-g3lo)
          zprim2=    Lz*(g3der2*a**2)/(g3up-g3lo)
          g3proc = z00 + Lz * (g3proc - g3lo) / (g3up - g3lo)
!
        case ('step-linear')
!
          xi_step(3,1)=xi_step_frac(3,1)*(nzgrid-1.0)
          xi_step(3,2)=xi_step_frac(3,2)*(nzgrid-1.0)
          dxyz_step(3,1)=(xyz_step(3,1)-z00)/(xi_step(3,1)-0.0)
          dxyz_step(3,2)=(xyz_step(3,2)-xyz_step(3,1))/ &
                                (xi_step(3,2)-xi_step(3,1))
          dxyz_step(3,3)=(z00+Lz-xyz_step(3,2))/(nzgrid-1.0-xi_step(3,2))
!
          call grid_profile(xi3,grid_func(3),g3,g3der1,g3der2, &
           dxyz=dxyz_step(3,:),xistep=xi_step(3,:),delta=xi_step_width(3,:))
          call grid_profile(xi3lo,grid_func(3),g3lo, &
           dxyz=dxyz_step(3,:),xistep=xi_step(3,:),delta=xi_step_width(3,:))
          call grid_profile(xi3proc, grid_func(3), g3proc, dxyz=dxyz_step(3,:), xistep=xi_step(3,:), delta=xi_step_width(3,:))
!
          z     = z00 + g3-g3lo
          zprim = g3der1
          zprim2= g3der2
          g3proc = z00 + g3proc - g3lo
!
        case ('linear+log')
          ! Grid distance is equidistant first and then increases logarithmically
          ! .i.e., grid spacing increases linearly
          ! d[log z] = const ==> dz = const*z
          a=4.0/float(nzgrid)
          b=200
          c=5
!
          call grid_profile(a*(xi3-b)  ,grid_func(3),g3,g3der1,g3der2,param=c)
          call grid_profile(a*(xi3lo-b),grid_func(3),g3lo,param=c)
          call grid_profile(a*(xi3up-b),grid_func(3),g3up,param=c)
          call grid_profile(a * (xi3proc - b), grid_func(3), g3proc, param=c)
!
          z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
          zprim =    Lz*(g3der1*a   )/(g3up-g3lo)
          zprim2=    Lz*(g3der2*a**2)/(g3up-g3lo)
          g3proc = z00 + Lz * (g3proc - g3lo) / (g3up - g3lo)
!
        case ('squared')
          ! Grid distance increases linearily
          a=max(nzgrid,1)
          !b=-max(nzgrid,1)/10  ! 23-nov-20/ccyang: is it correct to use integer arithmetic like this?
          b = 0.0
          call grid_profile(a*(xi3  -b),grid_func(3),g3,g3der1,g3der2)
          call grid_profile(a*(xi3lo-b),grid_func(3),g3lo)
          call grid_profile(a*(xi3up-b),grid_func(3),g3up)
          call grid_profile(a * (xi3proc - b), grid_func(3), g3proc)
!
          z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
          zprim =    Lz*(g3der1*a   )/(g3up-g3lo)
          zprim2=    Lz*(g3der2*a**2)/(g3up-g3lo)
          g3proc = z00 + Lz * (g3proc - g3lo) / (g3up - g3lo)
!
          if (lfirst_proc_z) then
            bound_prim1=z(n1+1)-z(n1)
            do i=1,nghost
              z(n1-i)=z(n1)-i*bound_prim1
              zprim(1:n1)=bound_prim1
            enddo
          endif
          if (llast_proc_z) then
            bound_prim2=z(n2)-z(n2-1)
            do i=1,nghost
              z(n2+i)=z(n2)+i*bound_prim2
              zprim(n2:mz)=bound_prim2
            enddo
          endif
!
        case ('trough')
          ! Grid distance is almost equidistant at boundaries and then decreases
          ! at the middle
          a=0.05
          b=47
          c=192-47
!
          call grid_profile(xi3  ,grid_func(3),g3,g3der1,g3der2,param=a,xistep=(/b,c/),delta=(/0.5,0.1/))
          call grid_profile(xi3lo,grid_func(3),g3lo,param=a,xistep=(/b,c/),delta=(/0.5,0.1/))
          call grid_profile(xi3up,grid_func(3),g3up,param=a,xistep=(/b,c/),delta=(/0.5,0.1/))
          call grid_profile(xi3proc, grid_func(3), g3proc, param=a, xistep=(/b,c/), delta=(/0.5,0.1/))
!
          z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
          zprim =    Lz*g3der1/(g3up-g3lo)
          zprim2=    Lz*g3der2/(g3up-g3lo)
          g3proc = z00 + Lz * (g3proc - g3lo) / (g3up - g3lo)
!
        case default
          call fatal_error('construct_grid', &
                           'No such z grid function - '//trim(grid_func(3)))
        endselect
!
        dz2=zprim**2
        dz_1=1./zprim
        dz_tilde=-zprim2/dz2
!
      endif
!
!  Record the processor boundary.
!
      procx_bounds = g1proc
      procy_bounds = g2proc
      procz_bounds = g3proc
!
      call grid_bound_data
!
!  Fargo (orbital advection acceleration) is implemented for polar coordinates only.
!  Die otherwise.       
!
      if (lfargo_advection.and.coord_system=='cartesian') then
        if (lroot) then
          print*,""
          print*,"Fargo advection is only implemented for"
          print*,"polar coordinates. Switch"
          print*," coord_system='cylindric' or 'spherical'"
          print*,"in init_pars of start.in if you"
          print*,"want to use the fargo algorithm"
          print*,""
        endif
        call fatal_error("construct_grid","")
      endif
!
!  Set equator if possible.
!
      if (abs(xyz0(2)+xyz1(2))-pi<1e-3) lequatory=.true.
!
!  Set the serial grid arrays, that contain the coordinate values
!  from all processors.
!
!!!      call construct_serial_arrays  !!! creates trouble
!
    endsubroutine construct_grid
!***********************************************************************
    subroutine set_coorsys_dimmask
!
!   Sets switches for the different coordinate systems.
!
!  18-dec-15/MR: outsourced from initialize_grid
!  16-jan-17/MR: added initialization of dimensionality mask.
!
      lcartesian_coords=.false.
      lspherical_coords=.false.
      lcylindrical_coords=.false.
      lpipe_coords=.false.
!
      if (coord_system=='cartesian') then
        lcartesian_coords=.true.
!
!  Introduce new names (spherical_coords), in addition to the old ones.
!
      elseif (    coord_system=='spherical' &
              .or.coord_system=='spherical_coords') then
        lspherical_coords=.true.
!
!  Introduce new names (cylindrical_coords), in addition to the old ones.
!
      elseif (    coord_system=='cylindric' &
              .or.coord_system=='cylindrical_coords') then
        lcylindrical_coords=.true.
      else if (coord_system=='pipeflows') then
        lpipe_coords=.true.
      elseif (coord_system=='Lobachevskii') then
        call fatal_error('set_coorsys_dimmask', &
                         'Lobachevskii ccordinates not implemented')
      endif
!
!  Initialize dimensionality mask.
!
      if (nxgrid==1) then
        if (nygrid==1) then
          dim_mask(1)=3
        else
          dim_mask(1:2)=(/2,3/)
        endif
      else
        if (nygrid==1) dim_mask(1:2)=(/1,3/)
      endif

    endsubroutine set_coorsys_dimmask
!***********************************************************************
    subroutine initialize_grid
!
!  Coordinate-related issues: nonuniform meshes, different coordinate systems
!
!  20-jul-10/wlad: moved here from register
!  3-mar-14/MR: outsourced calculation of box_volume into box_vol
!  29-sep-14/MR: outsourced calculation of auxiliary quantities for curvilinear
!                coordinates into coords_aux; set coordinate switches at the
!                beginning
!   9-jun-15/MR: calculation of Area_* added
!  18-dec-15/MR: outsourced setting of switches to set_coords_switches; added 
!                initialization of Yin-Yang grid; added dimensionality mask:
!                lists the indices of the non-degenerate directions in the first 
!                dimensionality elements of dim_mask 
!  10-oct-17/MR: avoided communication in calculation of r_int and r_ext
!  10-jan-17/MR: moved call construct_serial_arrays to beginning
!
      use Sub, only: remove_prof
      use Mpicomm
      use IO, only: lcollective_IO
      use General, only: indgen, itoa, find_proc
!
      real :: fact, dxmin_x, dxmin_y, dxmin_z, dxmax_x, dxmax_y, dxmax_z
      integer :: xj,yj,zj,itheta,nphi,na,ne,mm,nn,nni
      integer, dimension(3) :: idzl,idzr
!
      call construct_serial_arrays
!
!  For curvilinear coordinate systems, calculate auxiliary quantities as, e.g., for spherical coordinates 1/r, cot(theta)/r, etc.
!
      call coords_aux(x,y,z)
!
!  determine global minimum and maximum of grid spacing in any direction
!
      if (lequidist(1) .or. nxgrid <= 1) then
        dxmin_x = dx
        dxmax_x = dx
      else
        dxmin_x = minval(xprim(l1:l2))
        dxmax_x = maxval(xprim(l1:l2))
      endif
!
      if (lequidist(2) .or. nygrid <= 1) then
        if (lspherical_coords.or.lcylindrical_coords) then
          dxmin_y = dy*minval(x(l1:l2))
          dxmax_y = dy*maxval(x(l1:l2))
        else
          dxmin_y = dy                 ! possibly incorrect for pipe coordinates
          dxmax_y = dy
        endif
      else
        dxmin_y = minval(yprim(m1:m2))
        dxmax_y = maxval(yprim(m1:m2))
      endif
!
      if (lequidist(3) .or. nzgrid <= 1) then
        if (lspherical_coords) then
          dxmin_z = dz*minval(x(l1:l2))*minval(sinth(m1:m2))
          dxmax_z = dz*maxval(x(l1:l2))*maxval(sinth(m1:m2))
        else
          dxmin_z = dz                ! possibly incorrect for pipe coordinates
          dxmax_z = dz
        endif
      else
        dxmin_z = minval(zprim(n1:n2))
        dxmax_z = maxval(zprim(n1:n2))
      endif
!
!  Find minimum/maximum grid spacing. Note that
!    minval( (/dxmin_x,dxmin_y,dxmin_z/), MASK=((/nxgrid,nygrid,nzgrid/) > 1) )
!  will be undefined if all n[xyz]grid==1, so we have to add the fourth
!  component with a test that is always true
!
      dxmin = minval( (/dxmin_x, dxmin_y, dxmin_z, huge(dx)/), &
                MASK=((/nxgrid, nygrid, nzgrid, 2/) > 1) )

      call mpiallreduce_min(dxmin,dxmin_x)
      dxmin=dxmin_x
!
      if (dxmin == 0) &
        call fatal_error ("initialize_grid", "check Lx,Ly,Lz: is one of them 0?", .true.)
!
      dxmax = maxval( (/dxmax_x, dxmax_y, dxmax_z, epsilon(dx)/), &
                MASK=((/nxgrid, nygrid, nzgrid, 2/) > 1) )

      call mpiallreduce_max(dxmax,dxmax_x)
      dxmax=dxmax_x
!       
!  Grid spacing. For non-equidistant grid or non-Cartesian coordinates
!  the grid spacing is calculated in the (m,n) loop.
!
      if (lcartesian_coords .and. all(lequidist)) then
!
!  FAG replaced old_cdtv flag with more general coordinate independent lmaximal   
!        if (old_cdtv) then
!          dxyz_2 = max(dx_1(l1:l2)**2,dy_1(m1)**2,dz_1(n1)**2)
!        else
          dline_1(:,1)=dx_1(l1:l2)
          dline_1(:,2)=dy_1(m1)
          dline_1(:,3)=dz_1(n1)
!
          if (lmaximal_cdtv) then
            dxyz_2 = max(dline_1(:,1)**2, dline_1(:,2)**2, dline_1(:,3)**2)
            dxyz_4 = max(dline_1(:,1)**4, dline_1(:,2)**4, dline_1(:,3)**4)
            dxyz_6 = max(dline_1(:,1)**6, dline_1(:,2)**6, dline_1(:,3)**6)
          else
            dxyz_2 = dline_1(:,1)**2 + dline_1(:,2)**2 + dline_1(:,3)**2
            dxyz_4 = dline_1(:,1)**4 + dline_1(:,2)**4 + dline_1(:,3)**4
            dxyz_6 = dline_1(:,1)**6 + dline_1(:,2)**6 + dline_1(:,3)**6
          endif
        !  dxyz_2 = dline_1(:,1)**2+dline_1(:,2)**2+dline_1(:,3)**2
        !  dxyz_4 = dline_1(:,1)**4+dline_1(:,2)**4+dline_1(:,3)**4
        !  dxyz_6 = dline_1(:,1)**6+dline_1(:,2)**6+dline_1(:,3)**6
!        endif
!
!  Fill pencil with maximum gridspacing. Will be overwritten
!  during the mn loop in the non-equidistant or non-Cartesian case.
!
        dxmax_pencil = dxmax
        dxmin_pencil = dxmin
!
      endif
!
! Box volume
!
      call box_vol
!
!  Initialize Yin-Yang grid.
!
      if (lyinyang) call yyinit
!
!  Volume element and area of coordinate surfaces.
!  Note that in the area a factor depending only on the coordinate x_i which defines the surface by x_i=const. is dropped.
!
      Area_xy=1.; Area_yz=1.; Area_xz=1.
!
      if (lcartesian_coords) then
!
!  x-extent
!
        if (nxgrid/=1) then
          dVol_x=xprim
          Area_xy=Area_xy*Lxyz(1)
          Area_xz=Area_xz*Lxyz(1)
        else
          dVol_x=1.
        endif
!
!  y-extent
!
        if (nygrid/=1) then
          dVol_y=yprim
          Area_xy=Area_xy*Lxyz(2)
          Area_yz=Area_yz*Lxyz(2)
        else
          dVol_y=1.
        endif
!
!  z-extent
!
        if (nzgrid/=1) then
          dVol_z=zprim
          Area_yz=Area_yz*Lxyz(3)
          Area_xz=Area_xz*Lxyz(3)
        else
          dVol_z=1.
        endif
!
!  single value volume element dVol applicable only to Cartesian equidistant grid
!
        if (all(lequidist)) dVol=dVol_x(l1:l2)*dVol_y(m1)*dVol_z(n1)
!
!  Spherical coordinate system
!
      elseif (lspherical_coords) then
!
!  Volume element - it is wrong for spherical, since
!  sinth also changes with y-position.
!
!  WL: Is this comment above still relevant?
!
!  Split up volume differential as (dr) * (r*dtheta) * (r*sinth*dphi)
!  and assume that sinth=1 if there is no theta extent.
!  This should always give a volume of 4pi/3*(r2^3-r1^3) for constant integrand
!  r extent:
!
        if (nxgrid/=1) then
          dVol_x=x**2*xprim
          Area_xy=Area_xy*1./2.*(xyz1(1)**2-xyz0(1)**2)
          Area_xz=Area_xz*1./2.*(xyz1(1)**2-xyz0(1)**2)
        else
!          dVol_x=1./3.*(xyz1(1)**3-xyz0(1)**3)   !???
!          Area_xy=Area_xy*1./2.*(xyz1(1)**2-xyz0(1)**2)    !???
!          Area_xz=Area_xz*1./2.*(xyz1(1)**2-xyz0(1)**2)
          dVol_x=1./3.
          Area_xy=Area_xy*1./2.
          Area_xz=Area_xz*1./2.

        endif
!
!  Theta extent (if non-radially symmetric)
!
        if (nygrid/=1) then
          dVol_y=sinth*yprim
          Area_xy=Area_xy*Lxyz(2)
          Area_yz=Area_yz*(cos(xyz0(2))-cos(xyz1(2)))
        else
          dVol_y=2.
          Area_xy=Area_xy*pi
          Area_yz=Area_yz*2.
        endif
!
!  phi extent (if non-axisymmetry)
!
        if (nzgrid/=1) then
          dVol_z=zprim
          Area_xz=Area_xz*Lxyz(3)
          Area_yz=Area_yz*Lxyz(3)
        else
          dVol_z=2.*pi
          Area_xz=Area_xz*2.*pi
          Area_yz=Area_yz*2.*pi
        endif
!
!  weighted coordinates for integration purposes
!  Need to modify for 2-D and 1-D cases!
!AB: for now, allow only if nxgrid>1. Dhruba, please check
!
        r2_weight=x(l1:l2)**2
        sinth_weight=sinth
        if (nxgrid>1) then
          do itheta=1,nygrid
            sinth_weight_across_proc(itheta)=sin(xyz0(2)+dy*itheta)
          enddo
        else
          sinth_weight_across_proc=1.
        endif
!
!  Calculate the volume of the box, for spherical coordinates.
!  MR: Why needed? Box_volume should be the same (but more accurate).
!      nVol, nVol1 never used!
!
        nVol=0.
        do xj=l1,l2
          do yj=m1,m2
            do zj=n1,n2
              nVol=nVol+x(xj)*x(xj)*sinth(yj)
            enddo
          enddo
        enddo
        nVol1=1./nVol
!
!  Trapezoidal rule
!
        if (lfirst_proc_x) r2_weight( 1)=.5*r2_weight( 1)
        if (llast_proc_x ) r2_weight(nx)=.5*r2_weight(nx)
!
        if (lfirst_proc_y) sinth_weight(m1)=.5*sinth_weight(m1)
        if (llast_proc_y ) sinth_weight(m2)=.5*sinth_weight(m2)
        sinth_weight_across_proc(1     )=0.5*sinth_weight_across_proc(1)
        sinth_weight_across_proc(nygrid)=0.5*sinth_weight_across_proc(nygrid)
!
!  End of coord_system=='spherical_coords' query.
!
      elseif (lcylindrical_coords) then
!
!  Volume element.
!
        if (nxgrid/=1) then
          dVol_x=x*xprim
          Area_xy=Area_xy*1./2.*(xyz1(1)**2-xyz0(1)**2)
          Area_xz=Area_xz*Lxyz(1)
        else
!          dVol_x=1./2.*(xyz1(1)**2-xyz0(1)**2)   !???
          dVol_x=1./2.
          Area_xy=Area_xy*dVol_x(1)
          Area_xz=Area_xz*Lxyz(1)
        endif
!
!  theta extent (non-cylindrically symmetric)
!
        if (nygrid/=1) then
          dVol_y=yprim
          Area_xy=Area_xy*Lxyz(2)
          Area_yz=Area_yz*Lxyz(2)
        else
          dVol_y=2.*pi
          Area_xy=Area_xy*2.*pi
          Area_yz=Area_yz*2.*pi
        endif
!
!  z extent (vertically extended)
!
        if (nzgrid/=1) then
          dVol_z=zprim
          Area_xz=Area_xz*Lxyz(3)
          Area_yz=Area_yz*Lxyz(3)
        else
          dVol_z=1.
        endif
!
!  Trapezoidal rule
!
        rcyl_weight=rcyl_mn
        if (lfirst_proc_x) rcyl_weight( 1)=.5*rcyl_weight( 1)
        if (llast_proc_x ) rcyl_weight(nx)=.5*rcyl_weight(nx)
!
!  Pipe coordinates (for hydraulic applications)
!
      elseif (lpipe_coords) then
!
!  Compute profile function.
!
        select case (pipe_func)
!
        case ('error_function')
          call warning('initialize_grid','calculation of Area_* not implemented')
          fact=.5/CrossSec_w**2
          glnCrossSec=glnCrossSec0*(-exp(-fact*(x(l1:l2)-CrossSec_x1)**2) &
                                    +exp(-fact*(x(l1:l2)-CrossSec_x2)**2))
!
        case default
          call fatal_error('initialize_grid', &
                           'no such pipe function - '//trim(pipe_func))
        endselect
!
!  Lobachevskii space
!
!      else
!
!  Stop if no existing coordinate system is specified
!
      else
        call fatal_error('initialize_grid', &
                         'no such coordinate system: '//coord_system)
      endif
!
!  Inverse volume elements
!
      dVol1_x = dVol_x
      dVol1_y = dVol_y
      dVol1_z = dVol_z
!
!  Avoid 1/0 at axis or in origin.
!
      if (lspherical_coords) then
        if (dVol1_x(l1) == 0.) dVol1_x(l1)=.5*x(l1+1)**2*xprim(l1+1)
        if (dVol1_y(m1) == 0.) dVol1_y(m1)=.5*sinth(m1+1)*yprim(m1+1)
        if (dVol1_y(m2) == 0.) dVol1_y(m2)=.5*sinth(m2-1)*yprim(m2-1)
      elseif (lcylindrical_coords) then
        if (dVol1_x(l1) == 0.) dVol1_x(l1)=.5*x(l1+1)*xprim(l1+1)
      endif

      dVol1_x = 1./dVol1_x
      dVol1_y = 1./dVol1_y
      dVol1_z = 1./dVol1_z
!
      if (lspherical_coords.or.lcylindrical_coords) then
!
!  Define inner and outer radii for non-cartesian coords.
!  If the user did not specify them yet (in start.in),
!  these are the first point of the first x-processor,
!  and the last point of the last x-processor.
!  r_int and r_ext can be taken from the domain's x-extent as periodic BC can't 
!  appear for r.
!
        if (r_int == 0)         r_int=xyz0(1)
        if (r_ext ==impossible) r_ext=xyz1(1)
        
        if (lroot.and.ip<14) print*,'initialize_grid, r_int,r_ext=',r_int,r_ext
      endif
!
!  For a non-periodic mesh, multiply boundary points by 1/2.
!  Do it for each direction in turn.
!  If a direction has no extent, it is automatically periodic
!  and the corresponding step is therefore not called.
!
      if (.not.lperi(1)) then
        if (lfirst_proc_x) dVol_x(l1)=.5*dVol_x(l1)
        if (llast_proc_x ) dVol_x(l2)=.5*dVol_x(l2)
      endif
!
      if (.not.lperi(2)) then
        if (lfirst_proc_y) dVol_y(m1)=.5*dVol_y(m1)
        if (llast_proc_y ) dVol_y(m2)=.5*dVol_y(m2)
      endif
!
      if (.not.lperi(3)) then
        if (lfirst_proc_z) dVol_z(n1)=.5*dVol_z(n1)
        if (llast_proc_z ) dVol_z(n2)=.5*dVol_z(n2)
      endif
!
! Check max possible ncoarse.
!
      if (nzgrid==1) then       ! switch off coarsening if no z-extent
        ncoarse=0
      elseif (ncoarse>nz/nghost) then
        ncoarse=nz/nghost
        call warning('initialize_grid','there are jumped-over processors due to grid coarsening'// &
                     ' -> ncoarse reduced to floor(nz/nghost)='//trim(itoa(ncoarse)))
      endif
!
!  Set lcoarse for coarsening grid near poles.
!
      lcoarse=lspherical_coords .and. lpole(2) .and. ncoarse>1 .and. nprocy>1
      if (.not.lcoarse) ncoarse=1
   
      if (lcoarse.and.(lfirst_proc_y.or.llast_proc_y)) then

!MR: TB generalized to more than two procs, in which coarsening happens.

        if (lfirst_proc_y) then
          mexts=(/m1,m1+ncoarse-2/)
        elseif (llast_proc_y) then
          mexts=(/m2-ncoarse+2,m2/)
        endif
!
!MR: missing - nprocy=1 case: two m intervals in one proc!
!
        allocate(nphis(mexts(1):mexts(2)), nphis1(mexts(1):mexts(2)), &
                 nphis2(mexts(1):mexts(2)))
        allocate(nexts(mexts(1):mexts(2),2))
        allocate(ninds(-nghost:nghost,mexts(1):mexts(2),n1:n2))
        ninds=0

        do mm=mexts(1),mexts(2)

          if (lfirst_proc_y) then
            nphi = ceiling(ncoarse/(mm-m1+1.))
          else
            nphi = ceiling(ncoarse/(m2-mm+1.))
          endif

          if (mod(nzgrid,nphi)/=0) &
            call warning('initialize_grid','coarsened grid non-periodic: nphi(m='// &
                         trim(itoa(mm))//')='//trim(itoa(nphi)))

          nphis(mm)=nphi
          nphis1(mm)=1./nphi
          nphis2(mm)=1./(nphi*nphi)
          na=n1 + mod((nphi-mod(nz,nphi))*ipz,nphi)

          do nn=na,n2,nphi
            ninds(0,mm,nn)=nn

            ninds(-1:-nghost:-1,mm,nn)=nn-indgen(nghost)*nphi
            where(ninds(-nghost:-1,mm,nn)<n1) &
              ninds(-nghost:-1,mm,nn)=(ninds(-nghost:-1,mm,nn)-n1+1)/nphi-1+n1

            ninds(+1:+nghost,mm,nn)=nn+indgen(nghost)*nphi

            where(ninds(+1:+nghost,mm,nn)>n2) &
              ninds(+1:+nghost,mm,nn)=(ninds(+1:+nghost,mm,nn)-n2-1)/nphi+1+n2

            ne=nn

            do nni=nn+1,min(n2,nn+nphi-1)
              idzl=nni-nn+(indgen(nghost)-1)*nphi
              idzr=nn+nphi-nni+(indgen(nghost)-1)*nphi
              ninds(:,mm,nni)=(/idzl(nghost:1:-1),-nn,idzr/)
            enddo

          enddo

          do nni=n1,na-1

            idzl=nphis(mm)-na+nni+(indgen(nghost)-1)*nphis(mm)
            idzr=na-nni+(indgen(nghost)-1)*nphis(mm)
            ninds(:,mm,nni)=(/idzl(nghost:1:-1),0,idzr/)

          enddo

          nexts(mm,:)=(/na,ne/)
        enddo

        if (lfirst_proc_x.and.lfirst_proc_z) then
          print*, 'On processors '//trim(itoa(iproc))//' ...(iprocy='//trim(itoa(ipy))// &
                  ')... '//trim(itoa(find_proc(nprocx-1,ipy,nprocz-1)))// &
                  ' the grid is coarsened for m = '// &
                  trim(itoa(mexts(1)))//' ... '//trim(itoa(mexts(2)))
          print*, 'with coarsening factors:'
          print'(30(1x,i2))', nphis
        endif
!write(iproc+30,*) 'nexts=', nexts(:,:)
!write(iproc+30,*) 'nphis=', nphis(:)
!write(iproc+30,'(7(i4,1x))')  (ninds(:,mm,:), mm=mexts(1),mexts(2))
      else
        lcoarse=.false.    ! now processor-dependent!
      endif
!
!  Print the value for which output is being produced.
!  (Have so far only bothered about single processor output.)
!  This is now done on all the processors.
!
      lpoint=min(max(l1,lpoint),l2)
      mpoint=min(max(m1,mpoint),m2)
      npoint=min(max(n1,npoint),n2)
      lpoint2=min(max(l1,lpoint2),l2)
      mpoint2=min(max(m1,mpoint2),m2)
      npoint2=min(max(n1,npoint2),n2)
!
!  Clean up profile files.
!
      if (lroot.or..not.lcollective_IO) call remove_prof('z')
      lwrite_prof=.true.

      if (lslope_limit_diff) then
        if (lroot) print*,'initialize_grid: Set up half grid x12, y12, z12'
        call generate_halfgrid(x12,y12,z12)
      endif
!
    endsubroutine initialize_grid
!***********************************************************************
    subroutine coarsegrid_interp(f,ivar1,ivar2)

      use General, only: ioptest

      real, dimension(:,:,:,:) :: f
      integer, optional :: ivar1,ivar2

      integer :: mm,ll,nn,coarse_neigh,iv,iv1,iv2
      integer, dimension(3), parameter :: left_inds=(/-3,-2,-1/), right_inds=(/1,2,3/)
      integer, dimension(6) :: neighs,dists
      real, dimension(6) :: ws
      real :: fac

      iv1=ioptest(ivar1,1)
      iv2=ioptest(ivar2,size(f,4))

      do mm=mexts(1),mexts(2)

        do nn=n1,n2
          if (ninds(0,mm,nn)<=0) then
!
!  Interpolation neighbors and distances:
!
            if (ninds(0,mm,nn)==0) then
              coarse_neigh=nexts(mm,1)             ! right coarsegrid neighbor
              neighs=(/ninds(left_inds,mm,coarse_neigh),coarse_neigh,ninds(right_inds(1:2),mm,coarse_neigh)/)
            else
              coarse_neigh=-ninds(0,mm,nn)         ! left coarsegrid neighbor
              neighs=(/ninds(left_inds(2:3),mm,coarse_neigh),coarse_neigh,ninds(right_inds,mm,coarse_neigh)/)
            endif

            dists=(/ninds(left_inds,mm,nn),-ninds(right_inds,mm,nn)/)

!  should be precalculated
            fac=1./nphis(mm)**5
            ws(1)=-0.0083333333333333*(dists(2)*dists(3)*dists(4)*dists(5)*dists(6))*fac
            ws(2)= 0.041666666666667 *(dists(1)*dists(3)*dists(4)*dists(5)*dists(6))*fac
            ws(3)=-0.083333333333333 *(dists(1)*dists(2)*dists(4)*dists(5)*dists(6))*fac
            ws(4)= 0.083333333333333 *(dists(1)*dists(2)*dists(3)*dists(5)*dists(6))*fac
            ws(5)=-0.041666666666667 *(dists(1)*dists(2)*dists(3)*dists(4)*dists(6))*fac
            ws(6)= 0.0083333333333333*(dists(1)*dists(2)*dists(3)*dists(4)*dists(5))*fac
if (abs(sum(ws)-1.)>1e-7) write(iproc+40,'(6(e12.5,1x), e12.5)') ws, sum(ws)
            do ll=l1,l2
              do iv=iv1,iv2
                f(ll,mm,nn,iv) = sum(ws*f(ll,mm,neighs,iv))
!if (ll==0.and.iv==5) then
!  write(iproc+50,'(i3,7(e14.7,1x))') nn, f(ll,mm,nn,iv), sum(ws*f(ll,mm,neighs,iv)), f(ll,mm,neighs,iv)
!endif
              enddo
            enddo
          endif
        enddo

      enddo

    endsubroutine coarsegrid_interp
!***********************************************************************
    subroutine quintic_interp(nn,a,ninds,dc)

      integer :: nn
      real, dimension(:,:) :: a
      integer, dimension(:) :: ninds
      real :: dc
! interpolation weights: 1/([-120,24,-12,12,-24,120]*dx^5) =
! [-0.00833333,0.0416667,-0.0833333,0.0833333,-0.0416667,0.00833333]*dx^-5

      integer :: iv,izu
      real, dimension(6), parameter :: coefs = &
            (/-0.00833333,0.0416667,-0.0833333,0.0833333,-0.0416667,0.00833333/)

      izu=ipz*nz+nghost

 !     dels=zgrid(inds)-zgrid(nn)
      do iv=1,mvar
        a(nn,iv) = sum(a(ninds,iv)*coefs)/dc**(-5)
      enddo

    endsubroutine quintic_interp
!***********************************************************************
    subroutine save_grid(lrestore)
!
!  Saves grid into local statics (needed for downsampled output)
!
!  6-mar-14/MR: coded
!
      use General, only : loptest
!
      logical, optional :: lrestore
!
      real, dimension(mx), save :: xs,dx_1s,dx_tildes
      real, dimension(my), save :: ys,dy_1s,dy_tildes
      real, dimension(mz), save :: zs,dz_1s,dz_tildes
      real, dimension(nx), save :: r_mns,r1_mns,r2_mns
      real, dimension(my), save :: sinths,sin1ths,sin2ths,cosths, &
                                   cotths,cos1ths,tanths
      real, dimension(mz), save :: sinphs, cosphs
      real, dimension(nx), save :: rcyl_mns,rcyl_mn1s,rcyl_mn2s

      real, save :: dxs, dys, dzs

      logical, save :: lfirst=.true.

      if (loptest(lrestore)) then
        if (lfirst) then
          call fatal_error('save_grid','first call must have lrestore=F')
        else
          dx=dxs; dy=dys; dz=dzs
          x=xs; y=ys; z=zs
          dx_1=dx_1s; dy_1=dy_1s; dz_1=dz_1s
          dx_tilde=dx_tildes; dy_tilde=dy_tildes; dz_tilde=dz_tildes

          if (lspherical_coords) then
            r_mn=r_mns; r1_mn=r1_mns; r2_mn=r2_mns
            sinth=sinths; sin1th=sin1ths; sin2th=sin2ths; costh=cosths
            cotth=cotths; cos1th=cos1ths; tanth=tanths
            sinphs=sinph; cosphs=cosph
          elseif (lcylindrical_coords) then
            rcyl_mn=rcyl_mns; rcyl_mn1=rcyl_mn1s; rcyl_mn2=rcyl_mn2s
          endif
        endif
      elseif (lfirst) then
        lfirst=.false.
        dxs=dx; dys=dy; dzs=dz
        xs=x; ys=y; zs=z
        dx_1s=dx_1; dy_1s=dy_1; dz_1s=dz_1
        dx_tildes=dx_tilde; dy_tildes=dy_tilde; dz_tildes=dz_tilde

        if (lspherical_coords) then
          r_mns=r_mn; r1_mns=r1_mn; r2_mns=r2_mn
          sinths=sinth; sin1ths=sin1th; sin2ths=sin2th; cosths=costh
          cotths=cotth; cos1ths=cos1th; tanths=tanth
          sinph=sinphs; cosph=cosphs
        elseif (lcylindrical_coords) then
          rcyl_mns=rcyl_mn; rcyl_mn1s=rcyl_mn1; rcyl_mn2s=rcyl_mn2
        endif
      endif

    endsubroutine save_grid
!***********************************************************************
    subroutine coords_aux(x,y,z)
!
!  6-mar-14/MR: outsourced from initialize_grid
!
      real, dimension(:) :: x,y,z
!
      real, dimension(size(y)) :: lat
      real, parameter :: sinth_min=1e-5,costh_min=1e-5
!
      if (lspherical_coords) then
!
! For spherical coordinates
!
        r_mn=x(l1:l2)
        if (x(l1)==0.) then
          r1_mn(2:)=1./x(l1+1:l2)
          r1_mn(1)=0.
        else
          r1_mn=1./x(l1:l2)
        endif
        r2_mn=r1_mn**2
!
!  Calculate sin(theta). Make sure that sinth=1 if there is no y extent,
!  regardless of the value of y. This is needed for correct integrations.
!
        if (ny==1) then
          sinth=1.
        else
          sinth=sin(y)
        endif
!
!  Calculate cos(theta) via latitude, which allows us to ensure
!  that sin(lat(midpoint)) = 0 exactly.
!
        if (luse_latitude) then
          lat=pi/2-y
          costh=sin(lat)
        else
          costh=cos(y)
        endif
!
!  Calculate 1/sin(theta). To avoid the axis we check that sinth
!  is always larger than a minmal value, sinth_min. The problem occurs
!  on theta=pi, because the theta range is normally only specified
!  with no more than 6 digits, e.g. theta = 0., 3.14159.
!
        where(abs(sinth)>sinth_min)
          sin1th=1./sinth
        elsewhere
          sin1th=0.
        endwhere
        sin2th=sin1th**2
!
!  Calculate cot(theta).
!
        cotth=costh*sin1th
!
!  Calculate 1/cos(theta). To avoid the equator we check that costh
!  is always larger than a minmal value, costh_min. The problem occurs
!  on theta=pi/2., because the theta range is normally only specified
!  with no more than 6 digits, e.g. theta = 0., 3.14159.
!
        where(abs(costh)>costh_min)
          cos1th=1./costh
        elsewhere
          cos1th=0.
        endwhere
!
!  Calculate tan(theta).
!
        tanth=sinth*cos1th
!
! also calculate sin(phi) and cos(phi) useful to calculate spherical harmonic
! decomposition.
!
        if (nzgrid.gt.1) then
          sinph=sin(z)
          cosph=cos(z)
        else
           sinph=0.;cosph=1
        endif
!
      elseif (lcylindrical_coords) then
!
!  Note: for consistency with spherical, 1/rcyl should really be rcyl1_mn,
!  not rcyl_mn1.
!
        rcyl_mn=x(l1:l2)
        if (x(l1)==0.) then
          rcyl_mn1(2:)=1./x(l1+1:l2)
          rcyl_mn1(1)=0.
        else
          rcyl_mn1=1./x(l1:l2)
        endif
        rcyl_mn2=rcyl_mn1**2
!
      endif
!
    endsubroutine coords_aux
!***********************************************************************
    subroutine box_vol
!
! calculates box volume
!
! 3-mar-14/MR: outsourced from initialize_grid
! 6-mar-14/MR: changed into subroutine setting global variable box_volume
!
      box_volume=1.
      if (lcartesian_coords) then
!
!  x,y,z-extent
!
        if (nxgrid/=1) box_volume = box_volume*Lxyz(1)
        if (nygrid/=1) box_volume = box_volume*Lxyz(2)
        if (nzgrid/=1) box_volume = box_volume*Lxyz(3)
!
!  Spherical coordinate system
!
      elseif (lspherical_coords) then 
!
!  Assume that sinth=1 if there is no theta extent.
!  This should always give a volume of 4pi/3*(r2^3-r1^3) for constant integrand
!
!  r extent
!
        if (nxgrid/=1) &
          box_volume = box_volume*1./3.*(xyz1(1)**3-xyz0(1)**3)
!
!  theta extent (if non-radially symmetric)
!
        if (nygrid/=1) then
          box_volume = box_volume*(-(cos(xyz1(2))  -cos(xyz0(2))))
        else
          box_volume = box_volume*2.
        endif
!
!  phi extent (if non-axisymmetry)
!
        if (nzgrid/=1) then
          box_volume = box_volume*Lxyz(3)
        else
          box_volume = box_volume*2.*pi
        endif
!
      elseif (lcylindrical_coords) then
!
!  Box volume and volume element.
!
        if (nxgrid/=1) &
          box_volume = box_volume*.5*(xyz1(1)**2-xyz0(1)**2)
!
!  theta extent (non-cylindrically symmetric)
!
        if (nygrid/=1) then
          box_volume = box_volume*Lxyz(2)
        else
          box_volume = box_volume*2.*pi
        endif
!
!  z extent (vertically extended)
!
        if (nzgrid/=1) box_volume = box_volume*Lxyz(3)
!
      endif
!
!  Volume calculation for pipe coordinates missing!
!
    endsubroutine box_vol
!***********************************************************************
    subroutine pencil_criteria_grid
!
!  All pencils that this special module depends on are specified here.
!
!  15-nov-06/tony: coded
!
      if (any(lfreeze_varext).or.any(lfreeze_varint)) then
        if (lcylinder_in_a_box.or.lcylindrical_coords) then
          lpenc_requested(i_rcyl_mn)=.true.
        else
          lpenc_requested(i_r_mn)=.true.
        endif
      endif
!
    endsubroutine pencil_criteria_grid
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
      if (lpencil_in(i_r_mn1)) lpencil_in(i_r_mn)=.true.
      if (lpencil_in(i_evr)) lpencil_in(i_r_mn)=.true.
      if (lpencil_in(i_evth).or.lpencil_in(i_evr)) then
         lpencil_in(i_pomx)=.true.
         lpencil_in(i_pomy)=.true.
         lpencil_in(i_rcyl_mn)=.true.
         lpencil_in(i_r_mn1)=.true.
      endif
      if (  lpencil_in(i_pomx) &
       .or. lpencil_in(i_pomy) &
       .or. lpencil_in(i_phix) &
       .or. lpencil_in(i_phiy)) then
        if (lcartesian_coords) then
          lpencil_in(i_rcyl_mn1)=.true.
        endif
      endif
      if (lspherical_coords.and.lpencil_in(i_phi_mn)) then
        lpencil_in(i_x_mn)=.true.
        lpencil_in(i_y_mn)=.true.
      endif
!
      if (lpencil_in(i_rr)) then
        lpencil_in(i_x_mn)=.true.
        lpencil_in(i_y_mn)=.true.
        lpencil_in(i_z_mn)=.true.
      endif
!
    endsubroutine pencil_interdep_grid
!***********************************************************************
    subroutine calc_pencils_grid_std(f,p)
!
! Envelope adjusting calc_pencils_hydro_pencpar to the standard use with
! lpenc_loc=lpencil
!
! 10-oct-17/MR: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call calc_pencils_grid_pencpar(f,p,lpencil)
!
    endsubroutine calc_pencils_grid_std
!***********************************************************************
    subroutine calc_pencils_grid_pencpar(f,p,lpenc_loc)
!
!  Calculate Grid/geometry related pencils. Uses arbitrary pencil mask lpenc_loc.
!  Most basic pencils should come first, as others may depend on them.
!
!   15-nov-06/tony: coded
!   27-aug-07/wlad: generalized for cyl. and sph. coordinates
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      intent(in) :: f
      intent(inout) :: p
      intent(in) :: lpenc_loc
!
      if (lcartesian_coords) then
!coordinates vectors
        if (lpenc_loc(i_x_mn))     p%x_mn    = x(l1:l2)
        if (lpenc_loc(i_y_mn))     p%y_mn    = spread(y(m),1,nx)
        if (lpenc_loc(i_z_mn))     p%z_mn    = spread(z(n),1,nx)
!spherical distance
        if (lpenc_loc(i_r_mn))     p%r_mn    = sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
!cylindrical distance (pomega)
        if (lpenc_loc(i_rcyl_mn))  p%rcyl_mn = sqrt(x(l1:l2)**2+y(m)**2)
!azimuthal angle (phi)
        if (lpenc_loc(i_phi_mn))   p%phi_mn  = atan2(y(m),x(l1:l2))
!inverse cylindrical distance 1/pomega
        if (lpenc_loc(i_rcyl_mn1)) p%rcyl_mn1=1./max(p%rcyl_mn,tini)
!inverse spherical distance 1/r
        if (lpenc_loc(i_r_mn1))    p%r_mn1   =1./max(p%r_mn,tini)
!pomega unit vectors: pomx=cos(phi) and pomy=sin(phi) where phi=azimuthal angle
        if (lpenc_loc(i_pomx))     p%pomx    = x(l1:l2)*p%rcyl_mn1
        if (lpenc_loc(i_pomy))     p%pomy    = y(  m  )*p%rcyl_mn1
!phi unit vectors
        if (lpenc_loc(i_phix))     p%phix    =-y(  m  )*p%rcyl_mn1
        if (lpenc_loc(i_phiy))     p%phiy    = x(l1:l2)*p%rcyl_mn1
!
      elseif (lcylindrical_coords) then
        if (lpenc_loc(i_x_mn))     p%x_mn    = x(l1:l2)*cos(y(m))
        if (lpenc_loc(i_y_mn))     p%y_mn    = x(l1:l2)*sin(y(m))
        if (lpenc_loc(i_z_mn))     p%z_mn    = spread(z(n),1,nx)
        if (lpenc_loc(i_r_mn))     p%r_mn    = sqrt(x(l1:l2)**2+z(n)**2)
        if (lpenc_loc(i_rcyl_mn))  p%rcyl_mn = x(l1:l2)
        if (lpenc_loc(i_phi_mn))   p%phi_mn  = spread(y(m),1,nx)
        if (lpenc_loc(i_rcyl_mn1)) p%rcyl_mn1=1./max(p%rcyl_mn,tini)
        if (lpenc_loc(i_r_mn1))    p%r_mn1   =1./max(p%r_mn,tini)
        if (lpenc_loc(i_pomx))     p%pomx    = 1.
        if (lpenc_loc(i_pomy))     p%pomy    = 0.
        if (lpenc_loc(i_phix))     p%phix    = 0.
        if (lpenc_loc(i_phiy))     p%phiy    = 1.
      elseif (lspherical_coords) then
        if (lpenc_loc(i_x_mn))     p%x_mn    = x(l1:l2)*sin(y(m))*cos(z(n))
        if (lpenc_loc(i_y_mn))     p%y_mn    = x(l1:l2)*sin(y(m))*sin(z(n))
        if (lpenc_loc(i_z_mn))     p%z_mn    = x(l1:l2)*cos(y(m))
        if (lpenc_loc(i_r_mn))     p%r_mn    = x(l1:l2)
        if (lpenc_loc(i_rcyl_mn))  p%rcyl_mn = x(l1:l2)*sin(y(m))
        if (lpenc_loc(i_phi_mn))   p%phi_mn  = spread(z(n),1,nx)
        if (lpenc_loc(i_rcyl_mn1)) p%rcyl_mn1=1./max(p%rcyl_mn,tini)
        if (lpenc_loc(i_r_mn1))    p%r_mn1   =1./max(p%r_mn,tini)
        if (lpenc_loc(i_pomx).or.lpenc_loc(i_pomy).or.&
            lpenc_loc(i_phix).or.lpenc_loc(i_phiy)) &
            call fatal_error('calc_pencils_grid', &
                'pomx, pomy, phix and phix not implemented for '// &
                'spherical polars')
      endif
!
!  set position vector
!
      if (lpenc_loc(i_rr)) then
        if (lcartesian_coords) then
          p%rr(:,1)=p%x_mn
          p%rr(:,2)=p%y_mn
          p%rr(:,3)=p%z_mn
         else
           call fatal_error('calc_pencils_grid', &
               'position vector not implemented for '//&
               'non-cartesian coordinates')
         endif
      endif
!
!  evr is the radial unit vector
!
      if (lpenc_loc(i_evr)) then
        if (lcartesian_coords) then
          p%evr(:,1) = p%rcyl_mn*p%r_mn1*p%pomx
          p%evr(:,2) = p%rcyl_mn*p%r_mn1*p%pomy
          p%evr(:,3) = z(n)*p%r_mn1
        else
          call fatal_error('calc_pencils_grid', &
              'radial unit vector not implemented for '//&
              'non-cartesian coordinates')
        endif
      endif
!
!  evth is the co-latitudinal unit vector
!
      if (lpenc_loc(i_evth)) then
        if (lcartesian_coords) then
          p%evth(:,1) = z(n)*p%r_mn1*p%pomx
          p%evth(:,2) = z(n)*p%r_mn1*p%pomy
          p%evth(:,3) = -p%rcyl_mn*p%r_mn1
        else
          call fatal_error('calc_pencils_grid', &
              'co-latitudinal unit vector not implemented for '//&
              'non-cartesian coordinates')
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_grid_pencpar
!***********************************************************************
    subroutine grid_profile_0D(xi,grid_func,g,gder1,gder2,param,dxyz,xistep,delta)
!
!  Scalar wrapper for the "elemental" subroutine 'grid_profile_1D'.
!
!  14-dec-10/Bourdin.KIS: coded
!
      real              :: xi
      character(len=*)  :: grid_func
      real              :: g
      real, optional    :: gder1, gder2
      real, optional    :: param
      real, optional, dimension(3) :: dxyz
      real, optional, dimension(2) :: xistep, delta
!
      intent(in)  :: xi, grid_func, param, dxyz, xistep, delta
      intent(out) :: g, gder1, gder2
!
      real, dimension(1) :: tmp_g, tmp_gder1, tmp_gder2
!
      call grid_profile ((/xi/), grid_func, tmp_g, tmp_gder1, tmp_gder2, param, dxyz, xistep, delta)
!
      g = tmp_g(1)
      if (present (gder1)) gder1 = tmp_gder1(1)
      if (present (gder2)) gder2 = tmp_gder2(1)
!
    endsubroutine grid_profile_0D
!***********************************************************************
    subroutine grid_profile_1D(xi,grid_func,g,gder1,gder2,param,dxyz,xistep,delta)
!
!  Specify the functional form of the grid profile function g
!  and calculate g,g',g''.
!
!  25-jun-04/tobi+wolf: coded
!   9-mar-17/MR: blocked used of unpresent optionals
!
      real, dimension(:)                    :: xi
      character(len=*)                      :: grid_func
      real, dimension(size(xi,1))           :: g
      real, dimension(size(xi,1)), optional :: gder1,gder2
      real, optional                        :: param
      real, optional, dimension(3) :: dxyz
      real, optional, dimension(2) :: xistep,delta
      real :: m
!
      intent(in)  :: xi,grid_func,param,dxyz,xistep,delta
      intent(out) :: g,gder1,gder2
!
      select case (grid_func)
!
      case ('linear')
        g=xi
        if (present(gder1)) gder1=1.0
        if (present(gder2)) gder2=0.0
!
      case ('sinh', 'sinh2')
        g=sinh(xi)
        if (present(gder1)) gder1=cosh(xi)
        if (present(gder2)) gder2=sinh(xi)
!
      case ('cos')
        ! Cos grid:
        ! Approximately equidistant at the boundaries, linear in the middle
        if (.not. present (param)) &
            call fatal_error ('grid_profile', "'cos' needs its parameter.")
!
        ! Transition goes from slope 1 (xi < pi/2) to slope param (xi > pi/2).
        m = 0.5 * (param - 1)
        where (xi <= -pi/2.0)
          g = xi + pi/2.0 + 1
        elsewhere (xi >= pi/2.0)
          g = param * (xi - pi/2.0) + pi/2.0 + m * pi/2.0 + 1 + pi/2.0 * (1 + m)
        elsewhere
          g = xi + m * (xi - cos (xi)) + 1 + pi/2.0 * (1 + m)
        endwhere
        if (present (gder1)) then
          where (xi <= -pi/2.0)
            gder1 = 1.0
          elsewhere (xi >= pi/2.0)
            gder1 = param
          elsewhere
            gder1 = 1 + 0.5 * (m - 1) * (1 + sin (xi))
          endwhere
        endif
        if (present (gder2)) then
          where ((xi <= -pi/2.0) .or. (xi >= pi/2.0))
            gder2 = 0.0
          elsewhere
            gder2 = 0.5 * (m - 1) * cos (xi)
          endwhere
        endif
!
      case ('arsinh')
        ! Area sinh grid:
        ! Approximately equidistant at the boundaries, linear in the middle
        if (.not. present (param)) &
            call fatal_error ('grid_profile', "'arsinh' needs its parameter.")
!
        ! ('asinh' is not available in F95, therefore it is replaced by 'ln'.)
        g = xi + param * (xi * log (xi + sqrt (xi**2 + 1)) - sqrt (xi**2 + 1))
        if (present (gder1)) gder1 = 1.0 + param * log (xi + sqrt (xi**2 + 1))
        if (present (gder2)) gder2 = param / (sqrt (xi**2 + 1))
!
      case ('tanh')
        ! Tanh grid:
        ! Approximately equidistant at the boundaries, linear in the middle
        if (.not. present (param)) &
            call fatal_error ('grid_profile', "'tanh' needs its parameter.")
!
        ! Transition goes from slope 1 (xi -> -oo) to slope param (xi -> +oo).
        m = 0.5 * (param - 1)
        g = xi * (m + 1) + m * log (cosh (xi))
        if (present (gder1)) gder1 = m * (1 + tanh (xi)) + 1
        if (present (gder2)) gder2 = m * (1 - tanh (xi)**2)
!
      case ('duct')
        g=sin(xi)
        if (present(gder1)) gder1= cos(xi)
        if (present(gder2)) gder2=-sin(xi)
!
      case ('half-duct')
! duct, but only on one boundary:
! Points are much denser near the boundaries than in the middle
        g=cos(xi)
        if (present(gder1)) gder1=-sin(xi)
        if (present(gder2)) gder2=-cos(xi)
!
      case ('squared')
        ! Grid distance increases linearily
        g=0.5*xi**2
        if (present(gder1)) gder1= xi
        if (present(gder2)) gder2= 0.
!
      case ('log')
        ! Grid distance increases logarithmically
        g=exp(xi)
        if (present(gder1)) gder1=  g
        if (present(gder2)) gder2=  g
!
      case ('linear+log')
        ! Grid distance is equidistant near lower boundary and then
        ! increases logarithmically
        if (present(param)) then
            g=(1+xi)*0.5*(1.-tanh(param*xi))+exp(xi)*0.5*(1.+tanh(param*xi))
          if (present(gder1)) then
              gder1=  0.5*((1.-tanh(param*xi))+ &
                (1+xi)*param*(tanh(param*xi)**2-1.0)+ &
                exp(xi)*(1.+param+tanh(param*xi)-param*tanh(param*xi)**2))
          endif
          if (present(gder2)) then
              gder2 = param*(tanh(param*xi)**2-1.) + &
                param**2*(1+xi)*(1.-tanh(param*xi)**2)*tanh(param*xi) + &
                0.5*exp(xi)*(1.+2*param+(1.-2*param**2)*tanh(param*xi)- &
                2*param*tanh(param*xi)**2+2*param**2*tanh(param*xi)**3)
          endif
        else
          g=exp(xi)
          if (present(gder1)) gder1=  g
          if (present(gder2)) gder2=  g
        endif
!
      case ('power-law')
        ! Grid distance increases according to a power-law
        if (.not. present (param)) &
            call fatal_error ('grid_profile', "'power-law' needs its parameter.")
        g=xi**(1./param)
        if (present(gder1)) gder1=1./param*              xi**(1/param-1)
        if (present(gder2)) gder2=1./param*(1./param-1.)*xi**(1/param-2)
!
      case ('trough')
        ! Grid distance larger and linear near boundaries and smaller and linear in
        ! middle. Needs xistep(2), delta(2) and param parameters specified in
        ! start.in. E.g.,for solar atmosphere can use delta=0.4,0.5, xistar=100,300,
        ! param=0.02 on a grid of mzgrid=646
        if (.not. present (param)) &
          call fatal_error ('grid_profile', "'trough' needs its parameter.")
        g=(-alog(exp(param*(xi-xistep(1)))+exp(-param*(xi-xistep(1))))/param+ &
            alog(exp(param*(xi-xistep(2)))+exp(-param*(xi-xistep(2))))/param+   &
            (2.+delta(1))*xi)/(2.+delta(1))
        if (present(gder1)) then 
          gder1=(2.+delta(1)-tanh(param*(xi-xistep(1)))+ &
                tanh(param*(xi-xistep(2))))/(2.+delta(1))
        endif
        if (present(gder2)) then 
          gder2=(-param*(1.-(tanh(param*(xi-xistep(1))))**2)+ &
                 param*(1.-(tanh(param*(xi-xistep(2))))**2))/ &
                (2.+delta(1))
        endif
!
      case ('frozensphere')
        ! Just like sinh, except set dx constant below a certain radius.
        m = 4.
        where (xi<0)
          g = m*xi
        elsewhere
          g=sinh(xi)
        endwhere
        if (present(gder1)) then
          where (xi<0)
            gder1 = m
          elsewhere
            gder1=cosh(xi)
          endwhere
        endif
        if (present(gder2)) then
          where (xi<0)
            gder2 = 0.
          elsewhere
            gder2=sinh(xi)
          endwhere
        endif
!
      case ('step-linear')
        if (.not. (present(dxyz) .and. present(xistep) .and. present(delta))) &
            call fatal_error('grid_profile',"'step-linear' needs its parameters.")
        if (xistep(1)/=0.0) then
          g=                                                                    &
           dxyz(1)*0.5*(xi-delta(1)*log(cosh(dble((xi-xistep(1))/delta(1))))) + &
           dxyz(2)*0.5*(   delta(1)*log(cosh(dble((xi-xistep(1))/delta(1))))  - &
                           delta(2)*log(cosh(dble((xi-xistep(2))/delta(2)))) )+ &
           dxyz(3)*0.5*(xi+delta(2)*log(cosh(dble((xi-xistep(2))/delta(2)))) )
!
          if (present(gder1)) then
            gder1=                                                      &
              dxyz(1)*0.5*(1.0-tanh(dble((xi-xistep(1))/delta(1))) ) +  &
              dxyz(2)*0.5*(    tanh(dble((xi-xistep(1))/delta(1)))   -  &
                               tanh(dble((xi-xistep(2))/delta(2))) ) +  &
              dxyz(3)*0.5*(1.0+tanh(dble((xi-xistep(2))/delta(2))) )
!
          endif
          if (present(gder2)) then
            gder2=  &
                + 0.5/delta(1)*(dxyz(2)-dxyz(1))/cosh(dble((xi-xistep(1))/delta(1)))**2 &
                + 0.5/delta(2)*(dxyz(3)-dxyz(2))/cosh(dble((xi-xistep(2))/delta(2)))**2
          endif
        else ! if xistep(1)=0.0, only a single step is needed
          g=                                                                    &
           dxyz(2)*0.5*(xi-delta(2)*log(cosh(dble((xi-xistep(2))/delta(2)))) )+ &
           dxyz(3)*0.5*(xi+delta(2)*log(cosh(dble((xi-xistep(2))/delta(2)))) )
!
          if (present(gder1)) then
            gder1=                                                      &
              dxyz(2)*0.5*(1.0-tanh(dble((xi-xistep(2))/delta(2))) ) +  &
              dxyz(3)*0.5*(1.0+tanh(dble((xi-xistep(2))/delta(2))) )
!
          endif
          if (present(gder2)) then
            gder2=  &
                + 0.5/delta(2)*(dxyz(3)-dxyz(2))/cosh(dble((xi-xistep(2))/delta(2)))**2
          endif
        endif
!
      case default
        call fatal_error('grid_profile','grid function not implemented: '//trim(grid_func))
!
      endselect
!
    endsubroutine grid_profile_1D
!***********************************************************************
    function find_star(xi_lo,xi_up,x_lo,x_up,x_star,grid_func) result (xi_star)
!
!  Finds the xi that corresponds to the inflection point of the grid-function
!  by means of a newton-raphson root-finding algorithm.
!
!  25-jun-04/tobi+wolf: coded
!
      real, intent(in) :: xi_lo,xi_up,x_lo,x_up,x_star
      character(len=*), intent(in) :: grid_func
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
        call grid_profile(xi_lo-xi_star,grid_func,g_lo,gder_lo)
        call grid_profile(xi_up-xi_star,grid_func,g_up,gder_up)
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
    subroutine real_to_index(n, x, xi)
!
!  Transforms coordinates in real space to those in index space.
!
!  10-sep-15/ccyang: coded.
!
      integer, intent(in) :: n
      real, dimension(n,3), intent(in) :: x
      real, dimension(n,3), intent(out) :: xi
!
      real, parameter :: ngp1 = nghost + 1
      integer :: i
!
!  Work on each direction.
!
      nonzero: if (n > 0) then
        dir: do i = 1, 3
          if (lactive_dimension(i)) then
            call inverse_grid(i, x(:,i), xi(:,i), local=.true.)
          else
            xi(:,i) = ngp1
          endif
        enddo dir
      endif nonzero
!
    endsubroutine real_to_index
!***********************************************************************
    subroutine inverse_grid(dir, x, xi, local)
!
!  Transform the x coordinates in real space to the xi coordinates in
!  index space in dir direction, where dir = 1, 2, or, 3.
!  If local is present and .true., the index space is with respect to
!  the local grid.
!
!  24-dec-14/ccyang: coded.
!
      use General, only: arcsinh
!
      integer, intent(in) :: dir
      real, dimension(:), intent(in) :: x
      real, dimension(:), intent(out) :: xi
      logical, intent(in), optional :: local
!
      integer, dimension(size(xi)) :: inear
      character(len=linelen) :: msg
      logical :: loc
      integer :: shift, inear_max
      real :: h, a, b, c, xiup
!
!  Sanity check.
!
      if (any(lshift_origin) .or. any(lshift_origin_lower)) &
          call fatal_error('inverse_grid', 'lshift_origin and lshift_origin_lower are not supported. ')
!
!  Global or local index space?
!
      loc = .false.
      if (present(local)) loc = local
!
!  Check the direction (assuming xilo = 0).
!
      h = 0.0
      xiup = 0.0
      shift = 0
      inear_max = nghost
      ckdir: select case (dir)
      case (1) ckdir
        h = dx
        xiup = real(nxgrid - merge(0, 1, lperi(1)))
        if (loc) shift = nx * ipx
        if (loc) inear_max = nghost + nx
      case (2) ckdir
        h = dy
        xiup = real(nygrid - merge(0, 1, (lperi(2) .or. lpole(2))))
        if (loc) shift = ny * ipy
        if (loc) inear_max = nghost + ny
      case (3) ckdir
        h = dz
        xiup = real(nzgrid - merge(0, 1, lperi(3)))
        if (loc) shift = nz * ipz
        if (loc) inear_max = nghost + nz
      case default ckdir
        write(msg,*) 'unknown direction dir = ', dir
        call fatal_error('inverse_grid', trim(msg))
      endselect ckdir
!
!  Make the inversion according to the grid function.
!
      func: select case (grid_func(dir))
!
      case ('linear') func
        xi = (x - xyz0(dir)) / h
!
      case ('sinh') func
        a = coeff_grid(dir) * Lxyz(dir)
        b = sinh(a)
        c = cosh(a) - 1.0
        a = (xyz_star(dir) - xyz0(dir)) / Lxyz(dir)
        a = a * b / sqrt(1.0 + 2.0 * a * (1.0 - a) * c)
        b = (sqrt(1.0 + a * a) * b - a * c) / Lxyz(dir)
        xi = (arcsinh(a) + arcsinh(b * (x - xyz0(dir)) - a)) / (coeff_grid(dir) * h)
!
      case ("log") func
        xi = xiup / log(xyz1(dir) / xyz0(dir)) * log(x / xyz0(dir))
!
      case default func
        call fatal_error('inverse_grid in grid.f90', &
          'unknown grid function ' // trim(grid_func(dir)))
!
      endselect func
!
!  Shift to match the global index space.
!
      if (lperi(dir)) then
        xi = xi + real(nghost) + 0.5
      else
        xi = xi + real(nghost + 1)
      endif
!
!  Convert to the local index space if requested.
!
      getloc: if (loc) then
        xi = xi - real(shift)
        inear = nint(xi)
        where ((x < xyz1(dir) .and. inear > inear_max) .or. &
               (x < xyz0(dir) .and. inear > nghost)) xi = nearest(xi, -1.0)
      endif getloc
!
    endsubroutine inverse_grid
!***********************************************************************
    subroutine symmetrize_grid(grid,ngrid,idim)

      integer                :: ngrid,idim
      real, dimension(ngrid) :: grid

      if ( lequidist(idim) .and. xyz0(idim)+xyz1(idim)<epsi ) then
        if ( mod(ngrid,4)==0 ) then
           grid(ngrid:3*ngrid/4+1:-1)=xyz1(idim)-grid(ngrid/2+1:3*ngrid/4)
           grid(1:ngrid/2)=-grid(ngrid:ngrid/2+1:-1)
        else
           grid(1:ngrid/2)=-grid(ngrid:ngrid/2+2:-1)
           grid(1:ngrid/2+1)=0.
        endif
      endif

    endsubroutine symmetrize_grid
!***********************************************************************
    subroutine construct_serial_arrays(lprecise_symmetry)
!
!  The arrays xyz are local only, yet sometimes the serial array is
!  needed. Construct here the serial arrays out of the local ones,
!  but gathering them processor-wise and broadcasting the constructed
!  array. This is only done in start time, so legibility (3 near-copies
!  of the same code) is preferred over code-reusability (one general
!  piece of code called three times).
!
!  19-oct-10/wlad: coded
!  08-may-12/ccyang: include dx_1, dx_tilde, ... arrays
!  25-feb-13/ccyang: construct global coordinates including ghost cells.
!
      use Mpicomm, only: mpisend_real,mpirecv_real,mpibcast_real, mpiallreduce_sum_int, MPI_COMM_WORLD
      use General, only: loptest
!
      logical, optional :: lprecise_symmetry

      real, dimension(nx) :: xrecv, x1recv, x2recv
      real, dimension(ny) :: yrecv, y1recv, y2recv
      real, dimension(nz) :: zrecv, z1recv, z2recv
      integer :: jx,jy,jz,iup,ido,iproc_recv
      integer :: iproc_first, iproc_last
!
      xrecv=0.; yrecv=0.; zrecv=0.
      x1recv=0.; y1recv=0.; z1recv=0.
      x2recv=0.; y2recv=0.; z2recv=0.
!
!  Serial x array
!
      if (iproc/=root) then
!
!  All processors of the same row (ipx,ipy or ipz)
!  send their array values to the root.
!
        if ((ipy==0).and.(ipz==0)) then
          call mpisend_real(x(l1:l2),nx,root,111)
          call mpisend_real(dx_1(l1:l2),nx,root,112)
          call mpisend_real(dx_tilde(l1:l2),nx,root,113)
        endif
      else
!
!  The root processor, in turn, receives the data from the others
!
        do jx=0,nprocx-1
          !avoid send-to-self
          if (jx/=root) then
!
!  Formula of the serial processor number:
!  iproc=ipx+nprocx*ipy+nprocx*nprocy*ipz
!  Since for the x-row ipy=ipz=0, this reduces
!  to iproc_recv=jx.
!
            iproc_recv=jx
            call mpirecv_real(xrecv,nx,iproc_recv,111)
            call mpirecv_real(x1recv,nx,iproc_recv,112)
            call mpirecv_real(x2recv,nx,iproc_recv,113)
!
            ido=jx    *nx + 1
            iup=(jx+1)*nx
            xgrid(ido:iup)=xrecv
            dx1grid(ido:iup)=x1recv
            dxtgrid(ido:iup)=x2recv
          else
            !the root just copies its value to the serial array
            xgrid(1:nx)=x(l1:l2)
            dx1grid(1:nx)=dx_1(l1:l2)
            dxtgrid(1:nx)=dx_tilde(l1:l2)
          endif
        enddo
      endif
!
!  Serial array constructed. Broadcast the result. Repeat the
!  procedure for y and z arrays.
!
      if (loptest(lprecise_symmetry)) call symmetrize_grid(xgrid,nxgrid,1)

      call mpibcast_real(xgrid,nxgrid,comm=MPI_COMM_WORLD)
      call mpibcast_real(dx1grid,nxgrid,comm=MPI_COMM_WORLD)
      call mpibcast_real(dxtgrid,nxgrid,comm=MPI_COMM_WORLD)
!
!  Serial y-array
!
      if (iproc/=root) then
        if (ipx==0.and.ipz==0) then
          call mpisend_real(y(m1:m2),ny,root,221)
          call mpisend_real(dy_1(m1:m2),ny,root,222)
          call mpisend_real(dy_tilde(m1:m2),ny,root,223)
        endif
      else
        do jy=0,nprocy-1
          if (jy/=root) then
            iproc_recv=nprocx*jy
            call mpirecv_real(yrecv,ny,iproc_recv,221)
            call mpirecv_real(y1recv,ny,iproc_recv,222)
            call mpirecv_real(y2recv,ny,iproc_recv,223)
            ido=jy    *ny + 1
            iup=(jy+1)*ny
            ygrid(ido:iup)=yrecv
            dy1grid(ido:iup)=y1recv
            dytgrid(ido:iup)=y2recv
          else
            ygrid(1:ny)=y(m1:m2)
            dy1grid(1:ny)=dy_1(m1:m2)
            dytgrid(1:ny)=dy_tilde(m1:m2)
          endif
        enddo
      endif

      if (loptest(lprecise_symmetry)) call symmetrize_grid(ygrid,nygrid,2)

      call mpibcast_real(ygrid,nygrid)
      call mpibcast_real(dy1grid,nygrid)
      call mpibcast_real(dytgrid,nygrid)
!
!  Serial z-array
!
      if (iproc/=root) then
        if (ipx==0.and.ipy==0) then
          call mpisend_real(z(n1:n2),nz,root,331)
          call mpisend_real(dz_1(n1:n2),nz,root,332)
          call mpisend_real(dz_tilde(n1:n2),nz,root,333)
        endif
      else
        do jz=0,nprocz-1
          if (jz/=root) then
            iproc_recv=nprocx*nprocy*jz
            call mpirecv_real(zrecv,nz,iproc_recv,331)
            call mpirecv_real(z1recv,nz,iproc_recv,332)
            call mpirecv_real(z2recv,nz,iproc_recv,333)
            ido=jz    *nz + 1
            iup=(jz+1)*nz
            zgrid(ido:iup)=zrecv
            dz1grid(ido:iup)=z1recv
            dztgrid(ido:iup)=z2recv
          else
            zgrid(1:nz)=z(n1:n2)
            dz1grid(1:nz)=dz_1(n1:n2)
            dztgrid(1:nz)=dz_tilde(n1:n2)
          endif
        enddo
      endif

      if (loptest(lprecise_symmetry)) call symmetrize_grid(zgrid,nzgrid,3)

      call mpibcast_real(zgrid,nzgrid)
      call mpibcast_real(dz1grid,nzgrid)
      call mpibcast_real(dztgrid,nzgrid)
!
!  Check the first and last processors.
!
      iup = 0
      if (lfirst_proc_xyz) iup = iproc
      ido = 0
      if (llast_proc_xyz) ido = iproc
      call mpiallreduce_sum_int(iup, iproc_first)
      call mpiallreduce_sum_int(ido, iproc_last)
!
!  Communicate the ghost cells.
!
      xglobal(nghost+1:mxgrid-nghost) = xgrid
      yglobal(nghost+1:mygrid-nghost) = ygrid
      zglobal(nghost+1:mzgrid-nghost) = zgrid
!
      xglobal(1:nghost) = x(1:nghost)
      yglobal(1:nghost) = y(1:nghost)
      zglobal(1:nghost) = z(1:nghost)
!
      xglobal(mxgrid-nghost+1:mxgrid) = x(mx-nghost+1:mx)
      yglobal(mygrid-nghost+1:mygrid) = y(my-nghost+1:my)
      zglobal(mzgrid-nghost+1:mzgrid) = z(mz-nghost+1:mz)
!
      call mpibcast_real(xglobal(1:nghost), nghost, iproc_first)
      call mpibcast_real(yglobal(1:nghost), nghost, iproc_first)
      call mpibcast_real(zglobal(1:nghost), nghost, iproc_first)
!
      call mpibcast_real(xglobal(mxgrid-nghost+1:mxgrid), nghost, iproc_last)
      call mpibcast_real(yglobal(mygrid-nghost+1:mygrid), nghost, iproc_last)
      call mpibcast_real(zglobal(mzgrid-nghost+1:mzgrid), nghost, iproc_last)
!
    endsubroutine construct_serial_arrays
!***********************************************************************
    subroutine get_grid_mn
!
!  Gets the geometry of the pencil at each (m,n) in the mn-loop.
!
!  03-jul-13/ccyang: extracted from Equ.
!
!      obsolete: if (old_cdtv) then
!       The following is only kept for backwards compatibility. Will be deleted in the future.
!        dxyz_2 = max(dx_1(l1:l2)**2, dy_1(m)**2, dz_1(n)**2)
!
!      else obsolete
!
        if (lspherical_coords) then
          dline_1(:,1) = dx_1(l1:l2)
          dline_1(:,2) = r1_mn * dy_1(m)
          dline_1(:,3) = r1_mn * sin1th(m) * dz_1(n)
          if (lcoarse_mn) dline_1(:,3) = dline_1(:,3)*nphis1(m)
        else if (lcylindrical_coords) then
          dline_1(:,1) = dx_1(l1:l2)
          dline_1(:,2) = rcyl_mn1 * dy_1(m)
          dline_1(:,3) = dz_1(n)
        else if (lcartesian_coords) then
          dline_1(:,1) = dx_1(l1:l2)
          dline_1(:,2) = dy_1(m)
          dline_1(:,3) = dz_1(n)
        else if (lpipe_coords) then
          dline_1(:,1) = dx_1(l1:l2)
          dline_1(:,2) = dy_1(m)
          dline_1(:,3) = dz_1(n)
        endif
!
        dxmax_pencil = 0.
        if (nxgrid /= 1) dxmax_pencil =     1.0 / dline_1(:,1)
        if (nygrid /= 1) dxmax_pencil = max(1.0 / dline_1(:,2), dxmax_pencil)
        if (nzgrid /= 1) dxmax_pencil = max(1.0 / dline_1(:,3), dxmax_pencil)
!
        dxmin_pencil = 0.
        if (nxgrid /= 1) dxmin_pencil =     1.0 / dline_1(:,1)
        if (nygrid /= 1) dxmin_pencil = min(1.0 / dline_1(:,2), dxmin_pencil)
        if (nzgrid /= 1) dxmin_pencil = min(1.0 / dline_1(:,3), dxmin_pencil)
!
        if (lmaximal_cdtv) then
          dxyz_2 = max(dline_1(:,1)**2, dline_1(:,2)**2, dline_1(:,3)**2)
          dxyz_4 = max(dline_1(:,1)**4, dline_1(:,2)**4, dline_1(:,3)**4)
          dxyz_6 = max(dline_1(:,1)**6, dline_1(:,2)**6, dline_1(:,3)**6)
        else
          dxyz_2 = dline_1(:,1)**2 + dline_1(:,2)**2 + dline_1(:,3)**2
          dxyz_4 = dline_1(:,1)**4 + dline_1(:,2)**4 + dline_1(:,3)**4
          dxyz_6 = dline_1(:,1)**6 + dline_1(:,2)**6 + dline_1(:,3)**6
        endif
!
        dVol = dVol_x(l1:l2)*dVol_y(m)*dVol_z(n)
!
!      endif obsolete
!
    endsubroutine get_grid_mn
!***********************************************************************
    subroutine calc_bound_coeffs(coors,coeffs)
!
!  Calculates the coefficients of the 6th order difference formula for the first
!  derivative at the boundary points.
!  The grid is provided in form of the coordinate vector coors.
!
!  26-mar-15/MR: extracted from deriv_alt.
!
  use Deriv, only: calc_coeffs_1
!
      real, dimension(:),                intent(IN ) :: coors
      real, dimension(-nghost:nghost,2), intent(OUT) :: coeffs
!
      real, dimension(-nghost+1:nghost) :: dc
      integer                           :: sc
!
      sc=size(coors)

      dc( 1:nghost)  = coors(nghost+2:2*nghost+1)-coors(nghost+1:2*nghost)
      dc(-nghost+1:0)= dc(nghost:1:-1)
!print*, 'DX,Y,Z=', dx,dy,dz
      call calc_coeffs_1(dc,coeffs(:,BOT))
      dc(-nghost+1:0)= coors(sc-2*nghost+1:sc-nghost)-coors(sc-2*nghost:sc-nghost-1)
      dc( 1:nghost)  = dc(0:-nghost+1:-1)

      call calc_coeffs_1(dc,coeffs(:,TOP))

    endsubroutine calc_bound_coeffs
!***********************************************************************
    subroutine grid_bound_data
!
!  Assign local boundary data after the grid is constructed.
!
!  08-apr-15/MR: coded
!  23-nov-20/ccyang: added xyz[01]_loc and Lxyz_loc
!
      xyz0_loc(1) = procx_bounds(ipx)
      xyz1_loc(1) = procx_bounds(ipx+1)
!
      xyz0_loc(2) = procy_bounds(ipy)
      xyz1_loc(2) = procy_bounds(ipy+1)
!
      xyz0_loc(3) = procz_bounds(ipz)
      xyz1_loc(3) = procz_bounds(ipz+1)
!
      Lxyz_loc = xyz1_loc - xyz0_loc
!
      if (nxgrid>1) then
        if (lfirst_proc_x) &
          dx2_bound(-1:-nghost:-1)= 2.*(x(l1+1:l1+nghost)-x(l1))
        if (llast_proc_x) &
          dx2_bound(nghost:1:-1)  = 2.*(x(l2)-x(l2-nghost:l2-1))
!
        call calc_bound_coeffs(x,coeffs_1_x)
      endif

      if (nygrid>1) then
        if (lfirst_proc_y) &
          dy2_bound(-1:-nghost:-1)= 2.*(y(m1+1:m1+nghost)-y(m1))
        if (llast_proc_y) &
          dy2_bound(nghost:1:-1)  = 2.*(y(m2)-y(m2-nghost:m2-1))
!
        call calc_bound_coeffs(y,coeffs_1_y)
      endif

      if (nzgrid>1) then
        if (lfirst_proc_z) &
          dz2_bound(-1:-nghost:-1)= 2.*(z(n1+1:n1+nghost)-z(n1))
        if (llast_proc_z) &
          dz2_bound(nghost:1:-1)  = 2.*(z(n2)-z(n2-nghost:n2-1))
!
        call calc_bound_coeffs(z,coeffs_1_z)
      endif
!
    endsubroutine grid_bound_data
!***********************************************************************
    subroutine generate_halfgrid(x12,y12,z12)
!
! x[l1:l2]+0.5/dx_1[l1:l2]-0.25*dx_tilde[l1:l2]/dx_1[l1:l2]^2
!
      real, dimension (mx), intent(out) :: x12
      real, dimension (my), intent(out) :: y12
      real, dimension (mz), intent(out) :: z12
!
      if (nxgrid == 1) then
        x12 = x
      else
        x12 = x + 0.5/dx_1 - 0.25*dx_tilde/dx_1**2
      endif
!
      if (nygrid == 1) then
        y12 = y
      else
        y12 = y + 0.5/dy_1 - 0.25*dy_tilde/dy_1**2
      endif
!
      if (nzgrid == 1) then
        z12 = z
      else
        z12 = z + 0.5/dz_1 - 0.25*dz_tilde/dz_1**2
      endif
!
    endsubroutine generate_halfgrid
!***********************************************************************
endmodule Grid
