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
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
!  Special definition for degenerate directions
!
  logical :: ldegenerate_directions=.false.
!
  private
!
  public :: construct_grid
  public :: pencil_criteria_grid
  public :: pencil_interdep_grid
  public :: calc_pencils_grid
  public :: initialize_grid
!
  interface grid_profile
    module procedure grid_profile_0D
    module procedure grid_profile_1D
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
      real, dimension(3,2) :: xi_step
      real, dimension(3,3) :: dxyz_step
      real :: dxmin_x,dxmax_x,dxmin_y,dxmax_y,dxmin_z,dxmax_z
      real :: xi1star,xi2star,xi3star,bound_prim1,bound_prim2
!
      real, dimension(mx) :: g1,g1der1,g1der2,xi1,xprim2
      real, dimension(my) :: g2,g2der1,g2der2,xi2,yprim2
      real, dimension(mz) :: g3,g3der1,g3der2,xi3,zprim2
!
      real, dimension(0:2*nprocx+1) :: xi1proc,g1proc
      real, dimension(0:2*nprocy+1) :: xi2proc,g2proc
      real, dimension(0:2*nprocz+1) :: xi3proc,g3proc
!
      real :: a,b,c
      integer :: i
!
      lequidist=(grid_func=='linear')
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
!  Produce index arrays for processor boundaries, which are needed for
!  particle migration (see redist_particles_bounds). The select cases
!  should use these arrays to set g{2,3}proc using the grid function.
!
      if (lparticles .or. lsolid_cells) then
        do i=0,nprocx
          xi1proc(2*i)  =i*nx-1
          xi1proc(2*i+1)=i*nx
        enddo
        do i=0,nprocy
          xi2proc(2*i)  =i*ny-1
          xi2proc(2*i+1)=i*ny
        enddo
        do i=0,nprocz
          xi3proc(2*i)  =i*nz-1
          xi3proc(2*i+1)=i*nz
        enddo
      endif
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
!
        select case (grid_func(1))
!
        case ('linear','sinh')
          a=coeff_grid(1)*dx
          xi1star=find_star(a*xi1lo,a*xi1up,x00,x00+Lx,xyz_star(1),grid_func(1))/a
          call grid_profile(a*(xi1  -xi1star),grid_func(1),g1,g1der1,g1der2)
          call grid_profile(a*(xi1lo-xi1star),grid_func(1),g1lo)
          call grid_profile(a*(xi1up-xi1star),grid_func(1),g1up)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
!
          if (lparticles .or. lsolid_cells) then
            call grid_profile(a*(xi1proc-xi1star),grid_func(1),g1proc)
            g1proc=x00+Lx*(g1proc  -  g1lo)/(g1up-g1lo)
          endif
!
        case ('cos','tanh')
          ! Approximately equidistant at the boundaries, linear in the middle
          a=pi*dx/trans_width(1)
          b=dxi_fact(1)
          xi1star = xi1lo + (xi1up-xi1lo) / (1.0 + (Lx/(xyz_star(1)-x00) - 1.0) / b)
          call grid_profile(a*(xi1  -xi1star),grid_func(1),g1,g1der1,g1der2,param=b)
          call grid_profile(a*(xi1lo-xi1star),grid_func(1),g1lo,param=b)
          call grid_profile(a*(xi1up-xi1star),grid_func(1),g1up,param=b)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
!
          if (lparticles .or. lsolid_cells) then
            call grid_profile(a*(xi1proc-xi1star),grid_func(1),g1proc,param=b)
            g1proc=x00+Lx*(g1proc  -  g1lo)/(g1up-g1lo)
          endif
!
        case ('arsinh')
          ! Approximately linear in infinity, linear in the middle
          a=pi*dx/trans_width(1)
          b=dxi_fact(1)
          xi1star = xi1lo + (xi1up-xi1lo) / (1.0 + (Lx/(xyz_star(1)-x00) - 1.0) / b)
          call grid_profile(a*(xi1  -xi1star),grid_func(1),g1,g1der1,g1der2,param=b)
          call grid_profile(a*(xi1lo-xi1star),grid_func(1),g1lo,param=b)
          call grid_profile(a*(xi1up-xi1star),grid_func(1),g1up,param=b)
!
          ! Slope should be 1 at the minimum:
          b = minval (g1der1)
          if (b < 1.0) then
            g1der1 = g1der1 + 1.0 - b
            g1     = g1   + (1.0 - b) * a * (xi1   - xi1star)
            g1lo   = g1lo + (1.0 - b) * a * (xi1lo - xi1star)
            g1up   = g1up + (1.0 - b) * a * (xi1up - xi1star)
          endif
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
!
          if (lparticles .or. lsolid_cells) then
            b=dxi_fact(3)
            call grid_profile(a*(xi1proc-xi1star),grid_func(1),g1proc,param=b)
            g1proc=x00+Lx*(g1proc  -  g1lo)/(g1up-g1lo)
          endif
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
!
          x     = x00 + g1-g1lo
          xprim = g1der1
          xprim2= g1der2
!
          if (lparticles .or. lsolid_cells) then
            call grid_profile(xi1proc,grid_func(1),g1proc, &
              dxyz=dxyz_step(1,:),xistep=xi_step(1,:),delta=xi_step_width(1,:))
            g1proc=x00+g1proc-g1lo
          endif
!
        case ('duct')
          a = pi/(max(nxgrid-1,1))
          call grid_profile(a*xi1  -pi/2,grid_func(1),g1,g1der1,g1der2)
          call grid_profile(a*xi1lo-pi/2,grid_func(1),g1lo)
          call grid_profile(a*xi1up-pi/2,grid_func(1),g1up)
!
          x     =x00+Lx*(g1-g1lo)/2
          xprim =    Lx*(g1der1*a   )/2
          xprim2=    Lx*(g1der2*a**2)/2
!
          if (lparticles .or. lsolid_cells) then
            call grid_profile(a*xi1proc-pi/2,grid_func(1),g1proc)
            g1proc=x00+Lx*(g1proc-g1lo)/2
            g1proc(0)=g1proc(1)-x(l1+1)+x(l1)
            g1proc(2*nprocx+1)=g1proc(2*nprocx)+x(l2)-x(l2-1)
          endif
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
!
          x     =x00+Lx*g1
          xprim =    Lx*(g1der1*a   )
          xprim2=    Lx*(g1der2*a**2)
!
          if (lparticles .or. lsolid_cells) then
             call fatal_error('construct_grid: non-equidistant grid', &
                  'half-duct not implemented for particles.')
!            call grid_profile(a*xi1proc-pi/2,grid_func(1),g1proc)
!            g1proc=x00+Lx*(g1proc-g1lo)/2
!            g1proc(0)=g1proc(1)-x(l1+1)+x(l1)
!            g1proc(2*nprocx+1)=g1proc(2*nprocx)+x(l2)-x(l2-1)
          endif
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
          call grid_profile(a*(xi1-b)  ,grid_func(1),g1,g1der1,g1der2,param=c)
          call grid_profile(a*(xi1lo-b),grid_func(1),g1lo,param=c)
          call grid_profile(a*(xi1up-b),grid_func(1),g1up,param=c)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
!
        case ('power-law')
          ! Grid distance increases as a power-law set by c=coeff_grid(1)
          ! i.e., grid distance increases as an inverse power-law
          ! d[x^c] = const ==> dx = const/x^(c-1)
          if (coeff_grid(1) == 0.) &
               call fatal_error('construct_grid:', 'Cannot create '//&
               'a grid for a power-law of zero. Please check.')
          c= coeff_grid(1)
          a= (xyz1(1)**c-xyz0(1)**c)/(xi1up-xi1lo)
          b= .5*(xi1up+xi1lo-(xyz1(1)**c+xyz0(1)**c)/a)
!
          call grid_profile(a*(xi1-b)  ,grid_func(1),g1,g1der1,g1der2,param=c)
          call grid_profile(a*(xi1lo-b),grid_func(1),g1lo,param=c)
          call grid_profile(a*(xi1up-b),grid_func(1),g1up,param=c)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)

        case ('squared')
          ! Grid distance increases linearily
          a=max(nxgrid,1)
          b=-max(nxgrid,1)/10
          !b=0.
          call grid_profile(a*(xi1  -b),grid_func(1),g1,g1der1,g1der2)
          call grid_profile(a*(xi1lo-b),grid_func(1),g1lo)
          call grid_profile(a*(xi1up-b),grid_func(1),g1up)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
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
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
!
          if (llast_proc_x) then
            bound_prim2=x(l2-2)-x(l2-3)
            do i=1,nghost+2
              x(l2-2+i)=x(l2-2)+i*bound_prim2
              xprim(l2-2:mx)=bound_prim2
            enddo
          endif
        case default
          call fatal_error('construct_grid', &
                           'No such x grid function - '//grid_func(1))
        endselect
!
        dx_1=1./xprim
        dx_tilde=-xprim2/xprim**2
!
!DM should this be xprim**3 ?
!        dx_tilde=-xprim2/xprim**3
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
        select case (grid_func(2))
!
        case ('linear','sinh')
!
          a=coeff_grid(2)*dy
          xi2star=find_star(a*xi2lo,a*xi2up,y00,y00+Ly,xyz_star(2),grid_func(2))/a
          call grid_profile(a*(xi2  -xi2star),grid_func(2),g2,g2der1,g2der2)
          call grid_profile(a*(xi2lo-xi2star),grid_func(2),g2lo)
          call grid_profile(a*(xi2up-xi2star),grid_func(2),g2up)
!
          y     =y00+Ly*(g2  -  g2lo)/(g2up-g2lo)
          yprim =    Ly*(g2der1*a   )/(g2up-g2lo)
          yprim2=    Ly*(g2der2*a**2)/(g2up-g2lo)
!
          if (lparticles .or. lsolid_cells) then
            call grid_profile(a*(xi2proc-xi2star),grid_func(2),g2proc)
            g2proc=y00+Ly*(g2proc  -  g2lo)/(g2up-g2lo)
          endif
!
        case ('cos','tanh')
          ! Approximately equidistant at the boundaries, linear in the middle
          a=pi*dy/trans_width(2)
          b=dxi_fact(2)
          xi2star = xi2lo + (xi2up-xi2lo) / (1.0 + (Ly/(xyz_star(2)-y00) - 1.0) / b)
          call grid_profile(a*(xi2  -xi2star),grid_func(2),g2,g2der1,g2der2,param=b)
          call grid_profile(a*(xi2lo-xi2star),grid_func(2),g2lo,param=b)
          call grid_profile(a*(xi2up-xi2star),grid_func(2),g2up,param=b)
!
          y     =y00+Ly*(g2  -  g2lo)/(g2up-g2lo)
          yprim =    Ly*(g2der1*a   )/(g2up-g2lo)
          yprim2=    Ly*(g2der2*a**2)/(g2up-g2lo)
!
          if (lparticles .or. lsolid_cells) then
            call grid_profile(a*(xi2proc-xi2star),grid_func(2),g2proc,param=b)
            g2proc=y00+Ly*(g2proc  -  g2lo)/(g2up-g2lo)
          endif
!
        case ('arsinh')
          ! Approximately linear in infinity, linear in the middle
          a=pi*dy/trans_width(2)
          b=dxi_fact(2)
          xi2star = xi2lo + (xi2up-xi2lo) / (1.0 + (Ly/(xyz_star(2)-y00) - 1.0) / b)
          call grid_profile(a*(xi2  -xi2star),grid_func(2),g2,g2der1,g2der2,param=b)
          call grid_profile(a*(xi2lo-xi2star),grid_func(2),g2lo,param=b)
          call grid_profile(a*(xi2up-xi2star),grid_func(2),g2up,param=b)
!
          ! Slope should be 1 at the minimum:
          b = minval (g2der1)
          if (b < 1.0) then
            g2der1 = g2der1 + 1.0 - b
            g2     = g2   + (1.0 - b) * a * (xi2   - xi2star)
            g2lo   = g2lo + (1.0 - b) * a * (xi2lo - xi2star)
            g2up   = g2up + (1.0 - b) * a * (xi2up - xi2star)
          endif
!
          y     =y00+Ly*(g2  -  g2lo)/(g2up-g2lo)
          yprim =    Ly*(g2der1*a   )/(g2up-g2lo)
          yprim2=    Ly*(g2der2*a**2)/(g2up-g2lo)
!
          if (lparticles .or. lsolid_cells) then
            b=dxi_fact(3)
            call grid_profile(a*(xi2proc-xi2star),grid_func(2),g2proc,param=b)
            g2proc=y00+Ly*(g2proc  -  g2lo)/(g2up-g2lo)
          endif
!
        case ('duct')
!
          a = pi/max(nygrid-1, 1)
          call grid_profile(a*xi2  -pi/2,grid_func(2),g2,g2der1,g2der2)
          call grid_profile(a*xi2lo-pi/2,grid_func(2),g2lo)
          call grid_profile(a*xi2up-pi/2,grid_func(2),g2up)
!
          y     =y00+Ly*(g2-g2lo)/2
          yprim =    Ly*(g2der1*a   )/2
          yprim2=    Ly*(g2der2*a**2)/2
!
          if (lparticles .or. lsolid_cells) then
            call grid_profile(a*xi2proc-pi/2,grid_func(2),g2proc)
            g2proc=y00+Ly*(g2proc-g2lo)/2
            g2proc(0)=g2proc(1)-y(m1+1)+y(m1)
            g2proc(2*nprocy+1)=g2proc(2*nprocy)+y(m2)-y(m2-1)
          endif
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
!
        case ('step-linear')
!
          xi_step(2,1)=xi_step_frac(2,1)*(nygrid-1.0)
          xi_step(2,2)=xi_step_frac(2,2)*(nygrid-1.0)
          dxyz_step(2,1)=(xyz_step(2,1)-y00)/(xi_step(2,1)-0.0)
          dxyz_step(2,2)=(xyz_step(2,2)-xyz_step(2,1))/ &
                                (xi_step(2,2)-xi_step(2,1))
          dxyz_step(2,3)=(y00+Ly-xyz_step(2,2))/(nygrid-1.0-xi_step(2,2))
!
          call grid_profile(xi2,grid_func(2),g2,g2der1,g2der2, &
           dxyz=dxyz_step(2,:),xistep=xi_step(2,:),delta=xi_step_width(2,:))
          call grid_profile(xi2lo,grid_func(2),g2lo, &
           dxyz=dxyz_step(2,:),xistep=xi_step(2,:),delta=xi_step_width(2,:))
          y     = y00 + g2-g2lo
          yprim = g2der1
          yprim2= g2der2
!
          if (lparticles .or. lsolid_cells) then
            call grid_profile(xi2proc,grid_func(2),g2proc, &
              dxyz=dxyz_step(2,:),xistep=xi_step(2,:),delta=xi_step_width(2,:))
            g2proc=y00+g2proc-g2lo
          endif
!
        case default
          call fatal_error('construct_grid', &
                           'No such y grid function - '//grid_func(2))
!
        endselect
!
! Added parts for spherical coordinates and cylindrical coordinates.
! From now on dy = d\theta but dy_1 = 1/rd\theta and similarly for \phi.
! corresponding r and rsin\theta factors for equ.f90 (where CFL timesteps
! are estimated) are removed.
        dy_1=1./yprim
        dy_tilde=-yprim2/yprim**2
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
!
          z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
          zprim =    Lz*(g3der1*a   )/(g3up-g3lo)
          zprim2=    Lz*(g3der2*a**2)/(g3up-g3lo)
!
          if (lparticles .or. lsolid_cells) then
            call grid_profile(a*(xi3proc-xi3star),grid_func(3),g3proc)
            g3proc=z00+Lz*(g3proc-g3lo)/(g3up-g3lo)
          endif
!
        case ('cos','tanh')
          ! Approximately equidistant at the boundaries, linear in the middle
          a=pi*dz/trans_width(3)
          b=dxi_fact(3)
          xi3star = xi3lo + (xi3up-xi3lo) / (1.0 + (Lz/(xyz_star(3)-z00) - 1.0) / b)
          call grid_profile(a*(xi3  -xi3star),grid_func(3),g3,g3der1,g3der2,param=b)
          call grid_profile(a*(xi3lo-xi3star),grid_func(3),g3lo,param=b)
          call grid_profile(a*(xi3up-xi3star),grid_func(3),g3up,param=b)
!
          z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
          zprim =    Lz*(g3der1*a   )/(g3up-g3lo)
          zprim2=    Lz*(g3der2*a**2)/(g3up-g3lo)
!
          if (lparticles .or. lsolid_cells) then
            call grid_profile(a*(xi3proc-xi3star),grid_func(3),g3proc,param=b)
            g3proc=z00+Lz*(g3proc  -  g3lo)/(g3up-g3lo)
          endif
!
        case ('arsinh')
          ! Approximately linear in infinity, linear in the middle
          a=pi*dz/trans_width(3)
          b=dxi_fact(3)
          xi3star = xi3lo + (xi3up-xi3lo) / (1.0 + (Lz/(xyz_star(3)-z00) - 1.0) / b)
          call grid_profile(a*(xi3  -xi3star),grid_func(3),g3,g3der1,g3der2,param=b)
          call grid_profile(a*(xi3lo-xi3star),grid_func(3),g3lo,param=b)
          call grid_profile(a*(xi3up-xi3star),grid_func(3),g3up,param=b)
!
          ! Slope should be 1 at the minimum:
          b = minval (g3der1)
          if ((grid_func(3) == 'arsinh') .and. (b < 1.0)) then
            g3der1 = g3der1 + 1.0 - b
            g3     = g3   + (1.0 - b) * a * (xi3   - xi3star)
            g3lo   = g3lo + (1.0 - b) * a * (xi3lo - xi3star)
            g3up   = g3up + (1.0 - b) * a * (xi3up - xi3star)
          endif
!
          z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
          zprim =    Lz*(g3der1*a   )/(g3up-g3lo)
          zprim2=    Lz*(g3der2*a**2)/(g3up-g3lo)
!
          if (lparticles .or. lsolid_cells) then
            b=dxi_fact(3)
            call grid_profile(a*(xi3proc-xi3star),grid_func(3),g3proc,param=b)
            g3proc=z00+Lz*(g3proc  -  g3lo)/(g3up-g3lo)
          endif
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
          z     = z00 + g3-g3lo
          zprim = g3der1
          zprim2= g3der2
!
          if (lparticles .or. lsolid_cells) then
            call grid_profile(xi3proc,grid_func(3),g3proc, &
              dxyz=dxyz_step(3,:),xistep=xi_step(3,:),delta=xi_step_width(3,:))
            g3proc=z00+g3proc-g3lo
          endif
!
        case default
          call fatal_error('construct_grid', &
                           'No such z grid function - '//grid_func(3))
        endselect
!
        dz_1=1./zprim
        dz_tilde=-zprim2/zprim**2
      endif
!
!  Compute averages across processor boundaries to calculate the physical
!  boundaries
!
      if (lparticles .or. lsolid_cells) then
        do i=0,nprocx
          procx_bounds(i)=(g1proc(2*i)+g1proc(2*i+1))*0.5
        enddo
        do i=0,nprocy
          procy_bounds(i)=(g2proc(2*i)+g2proc(2*i+1))*0.5
        enddo
        do i=0,nprocz
          procz_bounds(i)=(g3proc(2*i)+g3proc(2*i+1))*0.5
        enddo
      endif
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
        dxmin_y = dy
        if (lspherical_coords) dxmin_y = dy*minval(x(l1:l2))
        if (lcylindrical_coords) dxmin_y = dy*minval(x(l1:l2))
        dxmax_y = dy
        if (lspherical_coords) dxmax_y = dy*maxval(x(l1:l2))
        if (lcylindrical_coords) dxmax_y = dy*maxval(x(l1:l2))
      else
        dxmin_y = minval(yprim(m1:m2))
        dxmax_y = maxval(yprim(m1:m2))
      endif
      !
      if (lequidist(3) .or. nzgrid <= 1) then
        dxmin_z = dz
        dxmax_z = dz
        if (lspherical_coords) dxmin_z = dz*minval(x(l1:l2))*minval(sinth(m1:m2))
        if (lspherical_coords) dxmax_z = dz*maxval(x(l1:l2))*maxval(sinth(m1:m2))
      else
        dxmin_z = minval(zprim(n1:n2))
        dxmax_z = maxval(zprim(n1:n2))
      endif
!
      dxmin = minval( (/dxmin_x, dxmin_y, dxmin_z, huge(dx)/), &
                MASK=((/nxgrid, nygrid, nzgrid, 2/) > 1) )
!
      dxmax = maxval( (/dxmax_x, dxmax_y, dxmax_z, epsilon(dx)/), &
                MASK=((/nxgrid, nygrid, nzgrid, 2/) > 1) )
!
! Identify whether you are a processor at the pole or not
!
      if (nprocz>0) then
        if ((coord_system=='spherical').and. &
            (y(m1)==xyz0(2)).and.(y(m1)<pi/4.)) lnorth_pole=.true.
        if ((coord_system=='spherical').and. &
            (y(m2)==xyz1(2)).and.(y(m2)>3*pi/4.)) lsouth_pole=.true.
      endif
!
    endsubroutine construct_grid
!***********************************************************************
    subroutine initialize_grid
!
!  Coordinate-related issues: nonuniform meshes, different coordinate systems
!
!  20-jul-10/wlad: moved here from register
!  21-apr-11/axel: insert special ldegenerate_directions option for degenerate directions
!
      use Sub, only: remove_zprof
      use Mpicomm
!
      real, dimension(my) :: lat
      real, dimension (nz,nprocz) :: z_allprocs_tmp
      real :: sinth_min=1e-5,costh_min=1e-5
      integer :: xj,yj,zj,itheta
!
! WL: Axel, I coded serial arrays xgrid, ygrid, zgrid, that maybe 
! make this z_allprocs obsolete. Can you check if you can use zgrid 
! instead of z_allprocs? (19-oct-2010)
!
!  Set z_allprocs, which contains the z values from all processors
!  ignore the ghost zones.
!
        z_allprocs(:,ipz+1)=z(n1:n2)
!
!  Communicate z_allprocs over all processors (if there are more than 1)
!  the final result is only present on the root processor.
!
      if (nprocz>1) then
        z_allprocs_tmp=z_allprocs
        call mpireduce_sum(z_allprocs_tmp,z_allprocs,(/nz,nprocz/))
      endif
!
!  For spherical coordinate system, calculate 1/r, cot(theta)/r, etc
!  Introduce new names (spherical_coords), in addition to the old ones.
!
      if (coord_system=='cartesian') then
        lcartesian_coords=.true.
        lspherical_coords=.false.
        lcylindrical_coords=.false.
!
!  Box volume and volume element.
!  x-extent
!
        box_volume=1.
        if (nxgrid/=1) then
          box_volume = box_volume*Lxyz(1)
          dVol_x=xprim
        else
          dVol_x=1.
        endif
!
!  y-extent
!
        if (nygrid/=1) then
          box_volume = box_volume*Lxyz(2)
          dVol_y=yprim
        else
          dVol_y=1.
        endif
!
!  z-extent
!
        if (nzgrid/=1) then
          box_volume = box_volume*Lxyz(3)
          dVol_z=zprim
        else
          dVol_z=1.
        endif
!
!  Spherical coordinate system
!
      elseif (coord_system=='spherical' &
        .or.coord_system=='spherical_coords') then
        lcartesian_coords=.false.
        lspherical_coords=.true.
        lcylindrical_coords=.false.
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
!  Calculate 1/cos(theta). To avoid the axis we check that costh
!  is always larger than a minmal value, costh_min. The problem occurs
!  on theta=pi, because the theta range is normally only specified
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
!  Box volume and volume element - it is wrong for spherical, since
!  sinth also changes with y-position.
!
!  WL: Is this comment above still relevant?
!
!  Split up volume differential as (dr) * (r*dtheta) * (r*sinth*dphi)
!  and assume that sinth=1 if there is no theta extent.
!  This should always give a volume of 4pi/3*(r2^3-r1^3) for constant integrand
!  r extent:
!
        box_volume=1.
        if (nxgrid/=1) then
          box_volume = box_volume*1./3.*(xyz1(1)**3-xyz0(1)**3)
          dVol_x=x**2*xprim
        else
!         
          if (ldegenerate_directions) then
            dVol_x=1.
          else
            dVol_x=1./3.*(xyz1(1)**3-xyz0(1)**3)
          endif
        endif
!
!  Theta extent (if non-radially symmetric)
!
        if (nygrid/=1) then
          box_volume = box_volume*(-(cos(xyz1(2))  -cos(xyz0(2))))
          dVol_y=sinth*yprim
        else
          if (ldegenerate_directions) then
            dVol_y=1.
          else
            box_volume = box_volume*2.
            dVol_y=2.
          endif
        endif
!
!  phi extent (if non-axisymmetry)
!
        if (nzgrid/=1) then
          box_volume = box_volume*Lxyz(3)
          dVol_z=zprim
        else
          if (ldegenerate_directions) then
            dVol_z=1.
          else
            box_volume = box_volume*2.*pi
            dVol_z=2.*pi
          endif
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
        endif
!
!  Calculate the volume of the box, for non-cartesian coordinates.
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
        if (lfirst_proc_x)  r2_weight( 1)=.5*r2_weight( 1)
        if (llast_proc_x) r2_weight(nx)=.5*r2_weight(nx)
!
        if (lfirst_proc_y)  sinth_weight(m1)=.5*sinth_weight(m1)
        if (llast_proc_y) sinth_weight(m2)=.5*sinth_weight(m2)
        sinth_weight_across_proc(1)=0.5*sinth_weight_across_proc(1)
        sinth_weight_across_proc(nygrid)=0.5*sinth_weight_across_proc(nygrid)
!
!  End of coord_system=='spherical_coords' query.
!  Introduce new names (cylindrical_coords), in addition to the old ones.
!
      elseif (coord_system=='cylindric' &
          .or.coord_system=='cylindrical_coords') then
        lcartesian_coords=.false.
        lspherical_coords=.false.
        lcylindrical_coords=.true.
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
!  Box volume and volume element.
!
        box_volume=1.
        if (nxgrid/=1) then
          box_volume = box_volume*.5*(xyz1(1)**2-xyz0(1)**2)
          dVol_x=x*xprim
        else
          if (ldegenerate_directions) then
            dVol_x=1.
          else
            dVol_x=1./2.*(xyz1(1)**2-xyz0(1)**2)
          endif
        endif
!
!  theta extent (non-cylindrically symmetric)
!
        if (nygrid/=1) then
          box_volume = box_volume*Lxyz(2)
          dVol_y=yprim
        else
          if (ldegenerate_directions) then
            dVol_y=1.
          else
            box_volume = box_volume*2.*pi
            dVol_y=2.*pi
          endif
        endif
!
!  z extent (vertically extended)
!
        if (nzgrid/=1) then
          box_volume = box_volume*Lxyz(3)
          dVol_z=zprim
        else
          dVol_z=1.
        endif
!
!  Trapezoidal rule
!
        rcyl_weight=rcyl_mn
        if (lfirst_proc_x)  rcyl_weight( 1)=.5*rcyl_weight( 1)
        if (llast_proc_x) rcyl_weight(nx)=.5*rcyl_weight(nx)
!
!  Lobachevskii space
!
      elseif (coord_system=='Lobachevskii') then
        lcartesian_coords=.false.
        lspherical_coords=.false.
        lcylindrical_coords=.false.
!
      endif
!
!  Inverse volume elements
!      
      dVol1_x = 1./dVol_x
      dVol1_y = 1./dVol_y
      dVol1_z = 1./dVol_z
!
!  Define inner and outer radii for non-cartesian coords.
!  If the user did not specify them yet (in start.in),
!  these are the first point of the first x-processor,
!  and the last point of the last x-processor.
!
      if (lspherical_coords.or.lcylindrical_coords) then
!
        if (nprocx/=1) then
!
!  The root (iproc=0) has by default the first value of x
!
          if (lroot) then
            if (r_int==0) r_int=x(l1)
!
!  The root should also receive the value of r_ext from
!  from the last x-processor (which is simply nprocx-1
!  for iprocy=0 and iprocz=0) for broadcasting.
!
            if (r_ext==impossible) &
                 call mpirecv_real(r_ext,1,nprocx-1,111)
          endif
!
!  The last x-processor knows the value of r_ext, and sends
!  it to root, for broadcasting.
!
          if ((r_ext==impossible).and.&
               (llast_proc_x.and.lfirst_proc_y.and.lfirst_proc_z)) then
            r_ext=x(l2)
            call mpisend_real(r_ext,1,0,111)
          endif
!
!  Broadcast the values of r_int and r_ext
!
          call mpibcast_real(r_int,1)
          call mpibcast_real(r_ext,1)
        else
!
!  Serial-x. Just get the local grid values.
!
          if (r_int == 0)         r_int=x(l1)
          if (r_ext ==impossible) r_ext=x(l2)
        endif
        if (lroot) print*,'initialize_modules, r_int,r_ext=',r_int,r_ext
      endif
!
!  For a non-periodic mesh, multiply boundary points by 1/2.
!  Do it for each direction in turn.
!  If a direction has no extent, it is automatically periodic
!  and the corresponding step is therefore not called.
!
      if (.not.lperi(1)) then
        if (lfirst_proc_x) dVol_x(l1)=.5*dVol_x(l1)
        if (llast_proc_x) dVol_x(l2)=.5*dVol_x(l2)
      endif
!
      if (.not.lperi(2)) then
        if (lfirst_proc_y) dVol_y(m1)=.5*dVol_y(m1)
        if (llast_proc_y)  dVol_y(m2)=.5*dVol_y(m2)
      endif
!
      if (.not.lperi(3)) then
        if (lfirst_proc_z) dVol_z(n1)=.5*dVol_z(n1)
        if (llast_proc_z)  dVol_z(n2)=.5*dVol_z(n2)
      endif
!
!  Print the value for which output is being produced.
!  (Have so far only bothered about single processor output.)
!
      if (lroot) then
        lpoint=min(max(l1,lpoint),l2)
        mpoint=min(max(m1,mpoint),m2)
        npoint=min(max(n1,npoint),n2)
        lpoint2=min(max(l1,lpoint2),l2)
        mpoint2=min(max(m1,mpoint2),m2)
        npoint2=min(max(n1,npoint2),n2)
        print*,'(x,y,z)(point)=',x(lpoint),y(mpoint),z(npoint)
        print*,'(x,y,z)(point2)=',x(lpoint2),y(mpoint2),z(npoint2)
      endif
!
!  Clean up profile files.
!
      call remove_zprof()
      lwrite_prof=.true.
!
!  Set the the serial grid arrays, that contain the coordinate values 
!  from all processors.
!
      call construct_serial_arrays()
!
    endsubroutine initialize_grid
!***********************************************************************
    subroutine pencil_criteria_grid()
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
      if (lcartesian_coords) then
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
      elseif (lcylindrical_coords) then
        if (lpencil(i_x_mn))     p%x_mn    = x(l1:l2)*cos(y(m))
        if (lpencil(i_y_mn))     p%y_mn    = x(l1:l2)*sin(y(m))
        if (lpencil(i_z_mn))     p%z_mn    = spread(z(n),1,nx)
        if (lpencil(i_r_mn))     p%r_mn    = sqrt(x(l1:l2)**2+z(n)**2)
        if (lpencil(i_rcyl_mn))  p%rcyl_mn = x(l1:l2)
        if (lpencil(i_phi_mn))   p%phi_mn  = spread(y(m),1,nx)
        if (lpencil(i_rcyl_mn1)) p%rcyl_mn1=1./max(p%rcyl_mn,tini)
        if (lpencil(i_r_mn1))    p%r_mn1   =1./max(p%r_mn,tini)
        if (lpencil(i_pomx))     p%pomx    = 1.
        if (lpencil(i_pomy))     p%pomy    = 0.
        if (lpencil(i_phix))     p%phix    = 0.
        if (lpencil(i_phiy))     p%phiy    = 1.
      elseif (lspherical_coords) then
        if (lpencil(i_x_mn))     p%x_mn    = x(l1:l2)*sin(y(m))*cos(z(n))
        if (lpencil(i_y_mn))     p%y_mn    = x(l1:l2)*sin(y(m))*sin(z(n))
        if (lpencil(i_z_mn))     p%z_mn    = x(l1:l2)*cos(y(m))
        if (lpencil(i_r_mn))     p%r_mn    = x(l1:l2)
        if (lpencil(i_rcyl_mn))  p%rcyl_mn = x(l1:l2)*sin(y(m))
        if (lpencil(i_phi_mn))   p%phi_mn  = spread(z(n),1,nx)
        if (lpencil(i_rcyl_mn1)) p%rcyl_mn1=1./max(p%rcyl_mn,tini)
        if (lpencil(i_r_mn1))    p%r_mn1   =1./max(p%r_mn,tini)
        if (lpencil(i_pomx).or.lpencil(i_pomy).or.&
            lpencil(i_phix).or.lpencil(i_phiy)) &
            call fatal_error('calc_pencils_grid', &
                'pomx, pomy, phix and phix not implemented for '// &
                'spherical polars')
      endif
!
!  set position vector
!
      if (lpencil(i_rr)) then
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
      if (lpencil(i_evr)) then
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
      if (lpencil(i_evth)) then
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
    endsubroutine calc_pencils_grid
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
      real, dimension(1) :: tmp_xi, tmp_g, tmp_gder1, tmp_gder2
      real :: tmp_param
      real, dimension(3) :: tmp_dxyz
      real, dimension(2) :: tmp_xistep, tmp_delta
!
!
      tmp_param = 0.0
      if (present (param)) tmp_param = param
      tmp_dxyz = 0.0
      if (present (dxyz)) tmp_dxyz = dxyz
      tmp_xistep = 0.0
      if (present (xistep)) tmp_xistep = xistep
      tmp_delta = 0.0
      if (present (delta)) tmp_delta = delta
!
      tmp_xi(:) = xi
!
      call grid_profile (tmp_xi, grid_func, tmp_g, tmp_gder1, tmp_gder2, tmp_param, tmp_dxyz, tmp_xistep, tmp_delta)
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
      case ('sinh')
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
      case ('power-law')
        ! Grid distance increases according to a power-law
        g=xi**(1./param)
        if (present(gder1)) gder1=1./param*              xi**(1/param-1)
        if (present(gder2)) gder2=1./param*(1./param-1.)*xi**(1/param-2)
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
        g=                                                                     &
         dxyz(1)*0.5*(xi - delta(1)*log(cosh(dble((xi-xistep(1))/delta(1))))) +&
         dxyz(2)*0.5*( delta(1)*log(cosh(dble((xi-xistep(1))/delta(1)))) -     &
                           delta(2)*log(cosh(dble((xi-xistep(2))/delta(2)))) )+&
         dxyz(3)*0.5*( xi + delta(2)*log(cosh(dble((xi-xistep(2))/delta(2)))) )
!
        if (present(gder1)) then
         gder1=                                                           &
            dxyz(1)*0.5*( 1.0 - tanh(dble((xi-xistep(1))/delta(1))) )  +  &
            dxyz(2)*0.5*( tanh(dble((xi-xistep(1))/delta(1)))  -          &
                              tanh(dble((xi-xistep(2))/delta(2))) )    +  &
            dxyz(3)*0.5*( 1.0 + tanh(dble((xi-xistep(2))/delta(2))) )
!
        endif
        if (present(gder2)) then
         gder2=                                                               &
          dxyz(1)*0.5* (-1.0)/delta(1)/cosh(dble((xi-xistep(1))/delta(1)))**2 +&
          dxyz(2)*0.5*(( 1.0)/delta(1)/cosh(dble((xi-xistep(1))/delta(1)))**2 -&
                       (-1.0)/delta(2)/cosh(dble((xi-xistep(2))/delta(2)))**2)+&
          dxyz(3)*0.5* ( 1.0)/delta(2)/cosh(dble((xi-xistep(2))/delta(2)))**2
!
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
    subroutine construct_serial_arrays
!
!  The arrays xyz are local only, yet sometimes the serial array is 
!  needed. Construct here the serial arrays out of the local ones, 
!  but gathering them processor-wise and broadcasting the constructed 
!  array. This is only done in start time, so legibility (3 near-copies 
!  of the same code) is preferred over code-reusability (one general 
!  piece of code called three times). 
!
!  19-oct-10/wlad: coded
!
      use Mpicomm, only: mpisend_real,mpirecv_real,mpibcast_real
!
      real, dimension(nx) :: xrecv,dx1recv,dVol1_x_recv
      real, dimension(ny) :: yrecv,dy1recv,dVol1_y_recv
      real, dimension(nz) :: zrecv,dz1recv,dVol1_z_recv
      integer :: jx,jy,jz,iup,ido,iproc_recv
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
          call mpisend_real(dx_1(l1:l2),nx,root,119)
          call mpisend_real(dVol1_x(l1:l2),nx,root,118)
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
            call mpirecv_real(dx1recv,nx,iproc_recv,119)
            call mpirecv_real(dVol1_x_recv,nx,iproc_recv,118)
!
            ido=jx    *nx + 1 
            iup=(jx+1)*nx
            xgrid(ido:iup)=xrecv
            dxgrid_1(ido:iup)=dx1recv
            dxvolume_1(ido:iup)=dVol1_x_recv
          else
            !the root just copies its value to the serial array
            xgrid(1:nx)=x(l1:l2)
            dxgrid_1(1:nx)=dx_1(l1:l2)
            dxvolume_1(1:nx)=dVol1_x(l1:l2)
          endif
        enddo
      endif
!
!  Serial array constructed. Broadcast the result. Repeat the 
!  procedure for y and z arrays.
!
      call mpibcast_real(xgrid,nxgrid)
      call mpibcast_real(dxgrid_1,nxgrid)
      call mpibcast_real(dxvolume_1,nxgrid)
!
!  Serial y-array
!
      if (iproc/=root) then
        if (ipx==0.and.ipz==0) then 
          call mpisend_real(y(m1:m2),ny,root,222)
          call mpisend_real(dy_1(m1:m2),ny,root,229)
          call mpisend_real(dVol1_y(m1:m2),ny,root,228)
        endif
      else
        do jy=0,nprocy-1
          if (jy/=root) then
            iproc_recv=nprocx*jy
            call mpirecv_real(yrecv,ny,iproc_recv,222)
            call mpirecv_real(dy1recv,ny,iproc_recv,229)
            call mpirecv_real(dVol1_y_recv,ny,iproc_recv,228)
            ido=jy    *ny + 1
            iup=(jy+1)*ny
            ygrid(ido:iup)=yrecv
            dygrid_1(ido:iup)=dy1recv
            dyvolume_1(ido:iup)=dVol1_y_recv
          else
            ygrid(1:ny)=y(m1:m2)
            dygrid_1(1:ny)=dy_1(m1:m2)
            dyvolume_1(1:ny)=dVol1_y(m1:m2)
          endif
        enddo
      endif
      call mpibcast_real(ygrid,nygrid)
      call mpibcast_real(dygrid_1,nygrid)
      call mpibcast_real(dyvolume_1,nygrid)
!
!  Serial z-array
!
      if (iproc/=root) then
        if (ipx==0.and.ipy==0) then 
          call mpisend_real(z(n1:n2),nz,root,333)
          call mpisend_real(dz_1(n1:n2),nz,root,339)
          call mpisend_real(dVol1_z(n1:n2),nz,root,338)
        endif
      else
        do jz=0,nprocz-1
          if (jz/=root) then
            iproc_recv=nprocx*nprocy*jz
            call mpirecv_real(zrecv,nz,iproc_recv,333)
            call mpirecv_real(dz1recv,nz,iproc_recv,339)
            call mpirecv_real(dvol1_z_recv,nz,iproc_recv,338)
            ido=jz    *nz + 1
            iup=(jz+1)*nz
            zgrid(ido:iup)=zrecv
            dzgrid_1(ido:iup)=dz1recv
            dzvolume_1(ido:iup)=dVol1_z_recv
          else
            zgrid(1:nz)=z(n1:n2)
            dzgrid_1(1:nz)=dz_1(n1:n2)
            dzvolume_1(1:nz)=dVol1_z(n1:n2)
          endif
        enddo
      endif
      call mpibcast_real(zgrid,nzgrid)
      call mpibcast_real(dzgrid_1,nzgrid)
      call mpibcast_real(dzvolume_1,nzgrid)
!
    endsubroutine construct_serial_arrays
!***********************************************************************
endmodule Grid
