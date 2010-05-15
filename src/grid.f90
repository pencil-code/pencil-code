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
  private
!
  public :: construct_grid
  public :: pencil_criteria_grid
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
      real :: a,b,dummy1=0.,dummy2=0.
      integer :: i
      logical :: err
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
      if (lparticles) then
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
        ! Test whether grid function is valid
        call grid_profile(dummy1,grid_func(1),dummy2,err=err)
        if (err) call &
            fatal_error('construct_grid','unknown grid_func '//grid_func(1))
!
        select case (grid_func(1))
!
        case ('linear','sinh')
          a=coeff_grid(1,1)*dx
          xi1star=find_star(a*xi1lo,a*xi1up,x00,x00+Lx,xyz_star(1),grid_func(1))/a
          call grid_profile(a*(xi1  -xi1star),grid_func(1),g1,g1der1,g1der2)
          call grid_profile(a*(xi1lo-xi1star),grid_func(1),g1lo)
          call grid_profile(a*(xi1up-xi1star),grid_func(1),g1up)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
!
          if (lparticles) then
            call grid_profile(a*(xi1proc-xi1star),grid_func(1),g1proc)
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
          if (lparticles) then
            g1proc=x00+g1proc-g1lo
            call grid_profile(xi1proc,grid_func(1),g1proc, &
              dxyz=dxyz_step(2,:),xistep=xi_step(2,:),delta=xi_step_width(2,:))
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
          if (lparticles) then
            g1proc=x00+Lx*(g1proc-g1lo)/2
            call grid_profile(a*xi1proc-pi/2,grid_func(1),g1proc)
            g1proc(0)=g1proc(1)-x(l1+1)+x(l1)
            g1proc(2*nprocx+1)=g1proc(2*nprocx)+x(l2)-x(l2-1)
          endif
!
          if (ipx==0) then
            bound_prim1=x(l1+1)-x(l1)
            do i=1,nghost
              x(l1-i)=x(l1)-i*bound_prim1
              xprim(1:l1)=bound_prim1
            enddo
          endif
          if (ipx==nprocx-1) then
            bound_prim2=x(l2)-x(l2-1)
            do i=1,nghost
              x(l2+i)=x(l2)+i*bound_prim2
              xprim(l2:mx)=bound_prim2
            enddo
          endif
!
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
          if (ipx==0) then
            bound_prim1=x(l1+1)-x(l1)
            do i=1,nghost
              x(l1-i)=x(l1)-i*bound_prim1
              xprim(1:l1)=bound_prim1
            enddo
          endif
          if (ipx==nprocx-1) then
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
          a=coeff_grid(1,1)*dx
          xi1star=find_star(a*xi1lo,a*xi1up,x00,x00+Lx,0.8*xyz_star(1),grid_func(1))/a
          call grid_profile(a*(xi1  -xi1star),grid_func(1),g1,g1der1,g1der2)
          call grid_profile(a*(xi1lo-xi1star),grid_func(1),g1lo)
          call grid_profile(a*(xi1up-xi1star),grid_func(1),g1up)
!
          x     =x00+Lx*(g1  -  g1lo)/(g1up-g1lo)
          xprim =    Lx*(g1der1*a   )/(g1up-g1lo)
          xprim2=    Lx*(g1der2*a**2)/(g1up-g1lo)
!
          if (ipx==nprocx-1) then
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
        ! Test whether grid function is valid
        call grid_profile(dummy1,grid_func(2),dummy2,err=err)
        if (err) &
            call fatal_error('construct_grid','unknown grid_func '//grid_func(2))
!
        select case (grid_func(2))
!
        case ('linear','sinh')
!
          a=coeff_grid(2,1)*dy
          xi2star=find_star(a*xi2lo,a*xi2up,y00,y00+Ly,xyz_star(2),grid_func(2))/a
          call grid_profile(a*(xi2  -xi2star),grid_func(2),g2,g2der1,g2der2)
          call grid_profile(a*(xi2lo-xi2star),grid_func(2),g2lo)
          call grid_profile(a*(xi2up-xi2star),grid_func(2),g2up)
!
          y     =y00+Ly*(g2  -  g2lo)/(g2up-g2lo)
          yprim =    Ly*(g2der1*a   )/(g2up-g2lo)
          yprim2=    Ly*(g2der2*a**2)/(g2up-g2lo)
!
          if (lparticles) then
            call grid_profile(a*(xi2proc-xi2star),grid_func(2),g2proc)
            g2proc=y00+Ly*(g2proc  -  g2lo)/(g2up-g2lo)
          endif
!
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
          if (lparticles) then
            call grid_profile(a*xi2proc-pi/2,grid_func(2),g2proc)
            g2proc=y00+Ly*(g2proc-g2lo)/2
            g2proc(0)=g2proc(1)-y(m1+1)+y(m1)
            g2proc(2*nprocy+1)=g2proc(2*nprocy)+y(m2)-y(m2-1)
          endif
!
          if (ipy==0) then
            bound_prim1=y(m1+1)-y(m1)
            do i=1,nghost
              y(m1-i)=y(m1)-i*bound_prim1
              yprim(1:m1)=bound_prim1
            enddo
          endif
          if (ipy==nprocy-1) then
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
          if (lparticles) then
            g2proc=y00+g2proc-g2lo
            call grid_profile(xi2proc,grid_func(2),g2proc, &
              dxyz=dxyz_step(2,:),xistep=xi_step(2,:),delta=xi_step_width(2,:))
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
        ! Test whether grid function is valid
        call grid_profile(dummy1,grid_func(3),dummy2,ERR=err)
        if (err) &
            call fatal_error('construct_grid','unknown grid_func '//grid_func(3))
!
        select case (grid_func(3))
!
        case ('linear','sinh')
!
          a=coeff_grid(3,1)*dz
          xi3star=find_star(a*xi3lo,a*xi3up,z00,z00+Lz,xyz_star(3),grid_func(3))/a
          call grid_profile(a*(xi3  -xi3star),grid_func(3),g3,g3der1,g3der2)
          call grid_profile(a*(xi3lo-xi3star),grid_func(3),g3lo)
          call grid_profile(a*(xi3up-xi3star),grid_func(3),g3up)
!
          z     =z00+Lz*(g3  -  g3lo)/(g3up-g3lo)
          zprim =    Lz*(g3der1*a   )/(g3up-g3lo)
          zprim2=    Lz*(g3der2*a**2)/(g3up-g3lo)
!
          if (lparticles) then
            g3proc=z00+Lz*(g3proc-g3lo)/(g3up-g3lo)
            call grid_profile(a*(xi3proc-xi3star),grid_func(3),g3proc)
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
          if (lparticles) then
            g3proc=z00+g3proc-g3lo
            call grid_profile(xi3proc,grid_func(2),g3proc, &
              dxyz=dxyz_step(3,:),xistep=xi_step(3,:),delta=xi_step_width(3,:))
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
      if (lparticles) then
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
    endsubroutine construct_grid
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
        if (lcartesian_coords) then
          lpencil_in(i_rcyl_mn1)=.true.
        endif
      endif
      if (lspherical_coords.and.lpencil_in(i_phi_mn)) then
        lpencil_in(i_x_mn)=.true.
        lpencil_in(i_y_mn)=.true.
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
      if (lcartesian_coords) then
!coordinates vectors
        if (lpencil(i_x_mn))     p%x_mn    = x(l1:l2)
        if (lpencil(i_y_mn))     p%y_mn    = spread(y(m),1,nx)
        if (lpencil(i_z_mn))     p%z_mn    = spread(z(n),1,nx)
!spherical distance
        if (lpencil(i_r_mn))     p%r_mn    = sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
        if (lpencil(i_r_mn).and.lspherical_coords) p%r_mn=x(l1:l2)
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
          p%evr(:,1) = p%x_mn
          p%evr(:,2) = p%y_mn
          p%evr(:,3) = p%z_mn
          p%evr = p%evr / spread(p%r_mn+tini,2,3)
        else
          call fatal_error('calc_pencils_grid', &
              'radial unit vector not implemented for '//&
              'non-cartesian coordinates')
        endif
      endif
!
!  evth is the latitudinal unit vector
!
      if (lpencil(i_evth)) then
        if (lcartesian_coords) then
          p%evth(:,1) = -z(n)*p%r_mn1*p%pomx
          p%evth(:,2) = -z(n)*p%r_mn1*p%pomy
          p%evth(:,3) = p%rcyl_mn*p%r_mn1
        else
          call fatal_error('calc_pencils_grid', &
              'latitudinal unit vector not implemented for '//&
              'non-cartesian coordinates')
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_grid
!***********************************************************************
    subroutine grid_profile_point(xi,grid_func,g,gder1,gder2,err, &
                                  dxyz,xistep,delta)
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
      real,optional,dimension(3) :: dxyz
      real,optional,dimension(2) :: xistep,delta
      real a1
!
      intent(in)  :: xi,grid_func,dxyz,xistep,delta
      intent(out) :: g,gder1,gder2,err
!
      if (present(err)) err=.false.
!
      select case (grid_func)
!
      case ('linear')
        ! Equidistant grid
        g=xi
        if (present(gder1)) gder1=1.0
        if (present(gder2)) gder2=0.0
!
      case ('sinh')
        ! Sinh grid:
        ! Approximately equidistant near the middle, but approximately
        ! exponential near the boundaries
        g=sinh(xi)
        if (present(gder1)) gder1=cosh(xi)
        if (present(gder2)) gder2=sinh(xi)
!
      case ('duct')
        ! Chebyshev-type grid in all Cartesian directions:
        ! Points are much denser near the boundaries than in the middle
        g=sin(xi)
        if (present(gder1)) gder1= cos(xi)
        if (present(gder2)) gder2=-sin(xi)
!
      case ('squared')
        ! Grid distance increases linearily
        g=0.5*xi**2
        if (present(gder1)) gder1= xi
        if (present(gder2)) gder2= 0.
!
      case ('frozensphere')
        ! Just like sinh, except set dx constant below a certain radius.
        a1 = 4.
        if (xi<0) then
          g = a1*xi
        else
          g=sinh(xi)
        endif
        if (present(gder1)) then
          if (xi<0) then
            gder1 = a1
          else
            gder1=cosh(xi)
          endif
        endif
        if (present(gder2)) then
          if (xi<0) then
            gder2 = 0.
          else
            gder2=sinh(xi)
          endif
        endif
!
      case ('step-linear')
       ! [Document me!
       ! This is certainly _not_ stepwise linear as the name would suggest]
       if (present(dxyz) .and. present(xistep) .and. present(delta)) then
        g=                                                                    &
         dxyz(1)*0.5*(xi-delta(1)*log(cosh(dble((xi-xistep(1))/delta(1))))) + &
         dxyz(2)*0.5*(delta(1)*log(cosh(dble((xi-xistep(1))/delta(1)))) -     &
                         delta(2)*log(cosh(dble((xi-xistep(2))/delta(2))))) + &
         dxyz(3)*0.5*(xi+delta(2)*log(cosh(dble((xi-xistep(2))/delta(2)))))
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
         gder2=                                                                &
          dxyz(1)*0.5*(-1.0)/delta(1)/cosh(dble((xi-xistep(1))/delta(1)))**2 + &
          dxyz(2)*0.5*((1.0)/delta(1)/cosh(dble((xi-xistep(1))/delta(1)))**2 - &
                      (-1.0)/delta(2)/cosh(dble((xi-xistep(2))/delta(2)))**2)+ &
          dxyz(3)*0.5*( 1.0)/delta(2)/cosh(dble((xi-xistep(2))/delta(2)))**2
!
        endif
       endif
!
      case default
        if (present(err)) err=.true.
!
      endselect
!
    endsubroutine grid_profile_point
!***********************************************************************
    subroutine grid_profile_1d(xi,grid_func,g,gder1,gder2,err, &
                               dxyz,xistep,delta)
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
      real, optional, dimension(3) :: dxyz
      real, optional, dimension(2) :: xistep,delta
      real a1
!
      intent(in)  :: xi,grid_func,dxyz,xistep,delta
      intent(out) :: g, gder1,gder2,err
!
      if (present(err)) err=.false.
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
      case ('duct')
        g=sin(xi)
        if (present(gder1)) gder1= cos(xi)
        if (present(gder2)) gder2=-sin(xi)
!
      case ('squared')
        ! Grid distance increases linearily
        g=0.5*xi**2
        if (present(gder1)) gder1= xi
        if (present(gder2)) gder2= 0.
!
      case ('frozensphere')
        ! Just like sinh, except set dx constant below a certain radius.
        a1 = 4.
        where (xi<0)
          g = a1*xi
        elsewhere
          g=sinh(xi)
        endwhere
        if (present(gder1)) then
          where (xi<0)
            gder1 = a1
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
       if (present(dxyz) .and. present(xistep) .and. present(delta)) then
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
       endif
!
      case default
        if (present(err)) err=.true.
!
      endselect
!
    endsubroutine grid_profile_1d
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
endmodule Grid
