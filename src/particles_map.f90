! $Id$
!
!  This module contains subroutines for mapping particles on the mesh.
!  Different domain decompositions have different versions of this module.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_blocks = .false.
!
!***************************************************************
module Particles_map
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_mpicomm
!
  implicit none
!
  include 'particles_map.h'
!
  interface interp_field_pencil_wrap
    module procedure interp_field_pencil_0
    module procedure interp_field_pencil_1
  endinterface
!
  interface interpolate_linear
    module procedure interpolate_linear_range
    module procedure interpolate_linear_scalar
  endinterface
!
  contains
!***********************************************************************
    subroutine initialize_particles_map()
!
!  Perform any post-parameter-read initialization.
!
!  29-mar-16/ccyang: coded.
!
!  Note: Currently, this subroutine is called after modules
!    Particles_mpicomm and Particles.
!
      integer :: i
      real :: total_gab_weights
!
      if (particles_module .ne. "lagrangian") then
!
!  Check the particle-mesh interpolation method.
!
      pm: select case (particle_mesh)
      case ('ngp', 'NGP') pm
!       Nearest-Grid-Point
        lparticlemesh_cic = .false.
        lparticlemesh_tsc = .false.
        lparticlemesh_gab = .false.
        if (lroot) print *, 'initialize_particles_map: selected nearest-grid-point for particle-mesh method. '
      case ('cic', 'CIC') pm
!       Cloud-In-Cell
        lparticlemesh_cic = .true.
        lparticlemesh_tsc = .false.
        lparticlemesh_gab = .false.
        if (lroot) print *, 'initialize_particles_map: selected cloud-in-cell for particle-mesh method. '
      case ('tsc', 'TSC') pm
!       Triangular-Shaped-Cloud
        lparticlemesh_cic = .false.
        lparticlemesh_tsc = .true.
        lparticlemesh_gab = .false.
        if (lroot) print *, 'initialize_particles_map: selected triangular-shaped-cloud for particle-mesh method. '
      case ('gab', 'GAB') pm
!       Gaussian box
        lparticlemesh_cic = .false.
        lparticlemesh_tsc = .false.
        lparticlemesh_gab = .true.
        if (lroot) print *, 'initialize_particles_map: selected gaussian-box for particle-mesh method. '

!NILS: Must figure out what to do with GAB!!!!!!!

        do i = 1,4
          gab_weights(i) = exp(-(real(i)-1.)**2/(gab_width**2))
        enddo
        total_gab_weights = sum(gab_weights)+sum(gab_weights(2:4))
        gab_weights = gab_weights/total_gab_weights
        if (lroot) print *,'The number of cells representing one standard deviation is: ', gab_width
        if (lroot) print *,'The weights for the gaussian box are: ', gab_weights
      case ('') pm
!       Let the logical switches decide.
!       TSC assignment/interpolation overwrites CIC in case they are both set.
        switch: if (lparticlemesh_tsc) then
          lparticlemesh_cic = .false.
          particle_mesh = 'tsc'
        elseif (lparticlemesh_cic) then switch
          particle_mesh = 'cic'
        else switch
          particle_mesh = 'ngp'
        endif switch
        if (lroot) print *, 'initialize_particles_map: particle_mesh = ' // trim(particle_mesh)
      case default pm
        call fatal_error('initialize_particles_map', 'unknown particle-mesh type ' // trim(particle_mesh))
      endselect pm

    endif

      if (lparticlemesh_gab) then
        lfold_df_3points=.true.
      endif
!
    endsubroutine initialize_particles_map
!***********************************************************************
    subroutine interpolate_linear_range(f,ivar1,ivar2,xxp,gp,inear,iblock,ipar)
!
!  Interpolate the value of g to arbitrary (xp, yp, zp) coordinate
!  using the linear interpolation formula
!
!    g(x,y,z) = A*x*y*z + B*x*y + C*x*z + D*y*z + E*x + F*y + G*z + H .
!
!  The coefficients are determined by the 8 grid points surrounding the
!  interpolation point.
!
!  30-dec-04/anders: coded
!
      use Solid_Cells
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      real, dimension (mvar) :: f_tmp
      integer, dimension (3) :: inear
      integer :: iblock, ipar
!
      real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, g5, g6, g7, g8
      real :: xp0, yp0, zp0
      real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
      integer :: ivar, i, ix0, iy0, iz0, icyl=1
      logical :: lfirstcall=.true.
!
      intent(in)  :: f, xxp, ivar1, inear
      intent(out) :: gp
!
!  Determine index value of lowest lying corner point of grid box surrunding
!  the interpolation point.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      if ( (x(ix0)>xxp(1)) .and. nxgrid/=1) ix0=ix0-1
      if ( (y(iy0)>xxp(2)) .and. nygrid/=1) iy0=iy0-1
      if ( (z(iz0)>xxp(3)) .and. nzgrid/=1) iz0=iz0-1
!
!  Check if the grid point interval is really correct.
!
      if ((x(ix0)<=xxp(1) .and. x(ix0+1)>=xxp(1)) .or. nxgrid==1) then
      ! x-dirn ok
      else
        print*, 'interpolate_linear: Interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'ipar = ', ipar
        print*, 'mx, x(1), x(mx) = ', mx, x(1), x(mx)
        print*, 'xp0, xp, xp1 = ', x(ix0), xxp(1), x(ix0+1)
        call fatal_error('interpolate_linear','point outside of x-interval')
      endif 
      if ((y(iy0)<=xxp(2) .and. y(iy0+1)>=xxp(2)) .or. nygrid==1) then
      ! y-dirn ok
      else
        print*, 'interpolate_linear: Interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'ipar = ', ipar
        print*, 'my, y(1), y(my) = ', my, y(1), y(my)
        print*, 'yp0, yp, yp1 = ', y(iy0), xxp(2), y(iy0+1)
        call fatal_error('interpolate_linear','point outside of y-interval')
      endif 
      if ((z(iz0)<=xxp(3) .and. z(iz0+1)>=xxp(3)) .or. nzgrid==1) then
      ! z-dirn ok
      else
        print*, 'interpolate_linear: Interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'ipar = ', ipar
        print*, 'mz, z(1), z(mz) = ', mz, z(1), z(mz)
        print*, 'zp0, zp, zp1 = ', z(iz0), xxp(3), z(iz0+1)
        call fatal_error('interpolate_linear','point outside of z-interval')
      endif 
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!  Set it equal to 0 for dimensions having 1 grid points; this will make sure
!  that the interpolation is bilinear for 2D grids.
!
      xp0=0; yp0=0; zp0=0
      if (nxgrid/=1) xp0=xxp(1)-x(ix0)
      if (nygrid/=1) yp0=xxp(2)-y(iy0)
      if (nzgrid/=1) zp0=xxp(3)-z(iz0)
!
!  Calculate derived grid spacing parameters needed for interpolation.
!  For an equidistant grid we only need to do this at the first call.
!
      if (lequidist(1)) then
        if (lfirstcall) dx1=dx_1(ix0) !1/dx
      else
        dx1=dx_1(ix0)
      endif
!
      if (lequidist(2)) then
        if (lfirstcall) dy1=dy_1(iy0)
      else
        dy1=dy_1(iy0)
      endif
!
      if (lequidist(3)) then
        if (lfirstcall) dz1=dz_1(iz0)
      else
        dz1=dz_1(iz0)
      endif
!
      if ( (.not. all(lequidist)) .or. lfirstcall) then
        dxdy1=dx1*dy1; dxdz1=dx1*dz1; dydz1=dy1*dz1
        dxdydz1=dx1*dy1*dz1
      endif
!
!  Function values at all corners.
!
      g1=f(ix0  ,iy0  ,iz0  ,ivar1:ivar2)
      g2=f(ix0+1,iy0  ,iz0  ,ivar1:ivar2)
      g3=f(ix0  ,iy0+1,iz0  ,ivar1:ivar2)
      g4=f(ix0+1,iy0+1,iz0  ,ivar1:ivar2)
      g5=f(ix0  ,iy0  ,iz0+1,ivar1:ivar2)
      g6=f(ix0+1,iy0  ,iz0+1,ivar1:ivar2)
      g7=f(ix0  ,iy0+1,iz0+1,ivar1:ivar2)
      g8=f(ix0+1,iy0+1,iz0+1,ivar1:ivar2)
!
!  Interpolation formula.
!
      gp = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) + &
          xp0*yp0*dxdy1*(g1-g2-g3+g4) + xp0*zp0*dxdz1*(g1-g2-g5+g6) + &
          yp0*zp0*dydz1*(g1-g3-g5+g7) + &
          xp0*yp0*zp0*dxdydz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
!
!  If we have solid geometry we might want some special treatment very close
!  to the surface of the solid geometry
!
      if (lsolid_cells.and.(.not.lsolid_ogrid)) then
        do ivar=ivar1,ivar2
          f_tmp(ivar)=gp(ivar-ivar1+1)
        enddo
        call close_interpolation(f,ix0,iy0,iz0,icyl,xxp,&
            f_tmp,.false.)
        do ivar=ivar1,ivar2
          gp(ivar-ivar1+1)=f_tmp(ivar)
        enddo
      endif
!
!  Do a reality check on the interpolation scheme.
!
      if (linterp_reality_check) then
        do i=1,ivar2-ivar1+1
          if (gp(i)>max(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
            print*, 'interpolate_linear: interpolated value is LARGER than'
            print*, 'interpolate_linear: a values at the corner points!'
            print*, 'interpolate_linear: ipar, xxp=', ipar, xxp
            print*, 'interpolate_linear: x0, y0, z0=', &
                x(ix0), y(iy0), z(iz0)
            print*, 'interpolate_linear: i, gp(i)=', i, gp(i)
            print*, 'interpolate_linear: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
          endif
          if (gp(i)<min(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
            print*, 'interpolate_linear: interpolated value is smaller than'
            print*, 'interpolate_linear: a values at the corner points!'
            print*, 'interpolate_linear: xxp=', xxp
            print*, 'interpolate_linear: x0, y0, z0=', &
                x(ix0), y(iy0), z(iz0)
            print*, 'interpolate_linear: i, gp(i)=', i, gp(i)
            print*, 'interpolate_linear: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
          endif
        enddo
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
      call keep_compiler_quiet(iblock)
!
    endsubroutine interpolate_linear_range
!***********************************************************************
    subroutine interpolate_linear_scalar(f,ivar,xxp,gp,inear,iblock,ipar)
!
!  Interpolate the value of g to arbitrary (xp, yp, zp) coordinate
!  using the linear interpolation formula
!
!    g(x,y,z) = A*x*y*z + B*x*y + C*x*z + D*y*z + E*x + F*y + G*z + H .
!
!  The coefficients are determined by the 8 grid points surrounding the
!  interpolation point.
!
!  30-dec-04/anders: coded
!
      use Solid_Cells
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (3) :: xxp
      real, dimension (mvar) :: f_tmp
      integer, dimension (3) :: inear
      integer :: iblock, ipar
!
      real :: gp, g1, g2, g3, g4, g5, g6, g7, g8
      real :: xp0, yp0, zp0
      real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
      integer :: ivar, ix0, iy0, iz0, icyl=1
      logical :: lfirstcall=.true.
!
      intent(in)  :: f, xxp, ivar, inear
      intent(out) :: gp
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      if ( (x(ix0)>xxp(1)) .and. nxgrid/=1) ix0=ix0-1
      if ( (y(iy0)>xxp(2)) .and. nygrid/=1) iy0=iy0-1
      if ( (z(iz0)>xxp(3)) .and. nzgrid/=1) iz0=iz0-1
!
!  Check if the grid point interval is really correct.
!
      if (((x(ix0)<=xxp(1) .and. x(ix0+1)>=xxp(1)) .or. nxgrid==1) .and. &
          ((y(iy0)<=xxp(2) .and. y(iy0+1)>=xxp(2)) .or. nygrid==1) .and. &
          ((z(iz0)<=xxp(3) .and. z(iz0+1)>=xxp(3)) .or. nzgrid==1)) then
        ! Everything okay
      else
        print*, 'interpolate_linear_scalar: Interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'ipar = ', ipar
        print*, 'mx, x(1), x(mx) = ', mx, x(1), x(mx)
        print*, 'my, y(1), y(my) = ', my, y(1), y(my)
        print*, 'mz, z(1), z(mz) = ', mz, z(1), z(mz)
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), x(ix0), x(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), y(iy0), y(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), z(iz0), z(iz0+1)
        call fatal_error('interpolate_linear_scalar','point outside of interval')
      endif
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!  Set it equal to 0 for dimensions having 1 grid points; this will make sure
!  that the interpolation is bilinear for 2D grids.
!
      xp0=0; yp0=0; zp0=0
      if (nxgrid/=1) xp0=xxp(1)-x(ix0)
      if (nygrid/=1) yp0=xxp(2)-y(iy0)
      if (nzgrid/=1) zp0=xxp(3)-z(iz0)
!
!  Calculate derived grid spacing parameters needed for interpolation.
!  For an equidistant grid we only need to do this at the first call.
!
      if (lequidist(1)) then
        if (lfirstcall) dx1=dx_1(ix0) !1/dx
      else
        dx1=dx_1(ix0)
      endif
!
      if (lequidist(2)) then
        if (lfirstcall) dy1=dy_1(iy0)
      else
        dy1=dy_1(iy0)
      endif
!
      if (lequidist(3)) then
        if (lfirstcall) dz1=dz_1(iz0)
      else
        dz1=dz_1(iz0)
      endif
!
      if ( (.not. all(lequidist)) .or. lfirstcall) then
        dxdy1=dx1*dy1; dxdz1=dx1*dz1; dydz1=dy1*dz1
        dxdydz1=dx1*dy1*dz1
      endif
!
!  Function values at all corners.
!
      g1=f(ix0  ,iy0  ,iz0  ,ivar)
      g2=f(ix0+1,iy0  ,iz0  ,ivar)
      g3=f(ix0  ,iy0+1,iz0  ,ivar)
      g4=f(ix0+1,iy0+1,iz0  ,ivar)
      g5=f(ix0  ,iy0  ,iz0+1,ivar)
      g6=f(ix0+1,iy0  ,iz0+1,ivar)
      g7=f(ix0  ,iy0+1,iz0+1,ivar)
      g8=f(ix0+1,iy0+1,iz0+1,ivar)
!
!  Interpolation formula.
!
      gp = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) + &
          xp0*yp0*dxdy1*(g1-g2-g3+g4) + xp0*zp0*dxdz1*(g1-g2-g5+g6) + &
          yp0*zp0*dydz1*(g1-g3-g5+g7) + &
          xp0*yp0*zp0*dxdydz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
!
!  If we have solid geometry we might want some special treatment very close
!  to the surface of the solid geometry
!
      if (lsolid_cells.and.(.not.lsolid_ogrid)) then
        f_tmp(ivar)=gp
        call close_interpolation(f,ix0,iy0,iz0,icyl,xxp,&
            f_tmp,.false.)
        gp=f_tmp(ivar)
      endif
!
!  Do a reality check on the interpolation scheme.
!
      if (linterp_reality_check) then
        if (gp>max(g1,g2,g3,g4,g5,g6,g7,g8)) then
          print*, 'interpolate_linear_scalar: interpolated value is LARGER than'
          print*, 'interpolate_linear_scalar: a values at the corner points!'
          print*, 'interpolate_linear_scalar: ipar, xxp=', ipar, xxp
          print*, 'interpolate_linear_scalar: x0, y0, z0=', &
              x(ix0), y(iy0), z(iz0)
          print*, 'interpolate_linear_scalar: gp=', gp
          print*, 'interpolate_linear_scalar: g1...g8=', &
              g1, g2, g3, g4, g5, g6, g7, g8
          print*, '------------------'
        endif
        if (gp<min(g1,g2,g3,g4,g5,g6,g7,g8)) then
          print*, 'interpolate_linear_scalar: interpolated value is smaller than'
          print*, 'interpolate_linear_scalar: a values at the corner points!'
          print*, 'interpolate_linear_scalar: xxp=', xxp
          print*, 'interpolate_linear_scalar: x0, y0, z0=', &
              x(ix0), y(iy0), z(iz0)
          print*, 'interpolate_linear_scalar: gp=', gp
          print*, 'interpolate_linear_scalar: g1...g8=', &
              g1, g2, g3, g4, g5, g6, g7, g8
          print*, '------------------'
        endif
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
      call keep_compiler_quiet(iblock)
!
    endsubroutine interpolate_linear_scalar
!***********************************************************************
    subroutine interpolate_quadratic(f,ivar1,ivar2,xxp,gp,inear,iblock,ipar)
!
!  Quadratic interpolation of g to arbitrary (xp, yp, zp) coordinate
!  using the biquadratic interpolation function
!
!    g(x,y,z) = (1+x+x^2)*(1+y+y^2)*(1+z+z^2)
!
!  The coefficients (9, one for each unique term) are determined by the 9
!  grid points surrounding the interpolation point.
!
!  The interpolation matrix M is defined through the relation
!    M#c = g
!  Here c are the coefficients and g is the value of the function at the grid
!  points. An equidistant grid has the following value of M:
!
!    invmat(:,1)=(/ 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00/)
!    invmat(:,2)=(/ 0.00, 0.00, 0.00,-0.50, 0.00, 0.50, 0.00, 0.00, 0.00/)
!    invmat(:,3)=(/ 0.00, 0.00, 0.00, 0.50,-1.00, 0.50, 0.00, 0.00, 0.00/)
!    invmat(:,4)=(/ 0.00,-0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00/)
!    invmat(:,5)=(/ 0.00, 0.50, 0.00, 0.00,-1.00, 0.00, 0.00, 0.50, 0.00/)
!    invmat(:,6)=(/ 0.25, 0.00,-0.25, 0.00, 0.00, 0.00,-0.25, 0.00, 0.25/)
!    invmat(:,7)=(/-0.25, 0.50,-0.25, 0.00, 0.00, 0.00, 0.25,-0.50, 0.25/)
!    invmat(:,8)=(/-0.25, 0.00, 0.25, 0.50, 0.00,-0.50,-0.25, 0.00, 0.25/)
!    invmat(:,9)=(/ 0.25,-0.50, 0.25,-0.50, 1.00,-0.50, 0.25,-0.50, 0.25/)
!
!    invmat(:,1)=invmat(:,1)
!    invmat(:,2)=invmat(:,2)/dx
!    invmat(:,3)=invmat(:,3)/dx**2
!    invmat(:,4)=invmat(:,4)/dz
!    invmat(:,5)=invmat(:,5)/dz**2
!    invmat(:,6)=invmat(:,6)/(dx*dz)
!    invmat(:,7)=invmat(:,7)/(dx**2*dz)
!    invmat(:,8)=invmat(:,8)/(dx*dz**2)
!    invmat(:,9)=invmat(:,9)/(dx**2*dz**2)
!
!  Space coordinates are defined such that the nearest grid point is at (0,0).
!  The grid points are counted from lower left:
!
!    7  8  9
!    4  5  6
!    1  2  3
!
!  The nearest grid point has index number 5.
!
!  09-jun-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (3) :: inear
      integer :: iblock, ipar
!
      real, dimension (9,ivar2-ivar1+1) :: cc
      real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, g5, g6, g7, g8, g9
      real :: dxp, dzp
      real, save :: dx1, dx2, dz1, dz2
      real, save :: dx1dz1, dx2dz1, dx1dz2, dx2dz2
      integer :: ix0, iy0, iz0
      logical, save :: lfirstcall=.true.
!
      intent(in)  :: f, xxp, ivar1, inear
      intent(out) :: gp
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
!
!  Not implemented in y-direction yet (but is easy to generalise).
!
      if (nygrid/=1) then
        if (lroot) print*, 'interpolate_quadratic: not implemented in y'
        call fatal_error('interpolate_quadratic','')
      endif
!
!  A few values that only need to be calculated once for equidistant grids.
!
      if (lequidist(1)) then
        if (lfirstcall) then
          dx1=1/dx; dx2=1/dx**2
        endif
      else
        dx1=dx_1(ix0); dx2=dx1**2
      endif
!
      if (lequidist(3)) then
        if (lfirstcall) then
          dz1=1/dz; dz2=1/dz**2
        endif
      else
        dz1=dz_1(iz0); dz2=dz1**2
      endif
!
      if (lequidist(1).and.lequidist(3)) then
        if (lfirstcall) then
          dx1dz1=1/(dx*dz)
          dx2dz1=1/(dx**2*dz); dx1dz2=1/(dx*dz**2); dx2dz2=1/(dx**2*dz**2)
        endif
      else
        dx1dz1=dx1*dz1
        dx2dz1=dx2*dz1; dx1dz2=dx1*dz1; dx2dz2=dx2*dz2
      endif
!
!  Define function values at the grid points.
!
      g1=f(ix0-1,iy0,iz0-1,ivar1:ivar2)
      g2=f(ix0  ,iy0,iz0-1,ivar1:ivar2)
      g3=f(ix0+1,iy0,iz0-1,ivar1:ivar2)
      g4=f(ix0-1,iy0,iz0  ,ivar1:ivar2)
      g5=f(ix0  ,iy0,iz0  ,ivar1:ivar2)
      g6=f(ix0+1,iy0,iz0  ,ivar1:ivar2)
      g7=f(ix0-1,iy0,iz0+1,ivar1:ivar2)
      g8=f(ix0  ,iy0,iz0+1,ivar1:ivar2)
      g9=f(ix0+1,iy0,iz0+1,ivar1:ivar2)
!
!  Calculate the coefficients of the interpolation formula (see introduction).
!
      cc(1,:)=                                g5
      cc(2,:)=dx1   *0.5 *(             -g4     +  g6           )
      cc(3,:)=dx2   *0.5 *(              g4-2*g5+  g6           )
      cc(4,:)=dz1   *0.5 *(     -g2                     +  g8   )
      cc(5,:)=dz2   *0.5 *(      g2        -2*g5        +  g8   )
      cc(6,:)=dx1dz1*0.25*( g1     -g3               -g7     +g9)
      cc(7,:)=dx2dz1*0.25*(-g1+2*g2-g3               +g7-2*g8+g9)
      cc(8,:)=dx1dz2*0.25*(-g1     +g3+2*g4     -2*g6-g7     +g9)
      cc(9,:)=dx2dz2*0.25*( g1-2*g2+g3-2*g4+4*g5-2*g6+g7-2*g8+g9)
!
!  Calculate the value of the interpolation function at the point (dxp,dzp).
!
      dxp=xxp(1)-x(ix0)
      dzp=xxp(3)-z(iz0)
!
      gp = cc(1,:)            + cc(2,:)*dxp        + cc(3,:)*dxp**2        + &
           cc(4,:)*dzp        + cc(5,:)*dzp**2     + cc(6,:)*dxp*dzp       + &
           cc(7,:)*dxp**2*dzp + cc(8,:)*dxp*dzp**2 + cc(9,:)*dxp**2*dzp**2
!
      call keep_compiler_quiet(ipar,iblock)
!
    endsubroutine interpolate_quadratic
!***********************************************************************
    subroutine interpolate_quadratic_spline(f,ivar1,ivar2,xxp,gp,inear,iblock,ipar)
!
!  Quadratic spline interpolation of the function g to the point xxp=(xp,yp,zp).
!
!  10-jun-06/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (3) :: inear
      integer :: iblock, ipar
!
      real :: fac_x_m1, fac_x_00, fac_x_p1
      real :: fac_y_m1, fac_y_00, fac_y_p1
      real :: fac_z_m1, fac_z_00, fac_z_p1
      real :: dxp0, dyp0, dzp0
      integer :: ix0, iy0, iz0
!
      intent(in)  :: f, xxp, ivar1, inear
      intent(out) :: gp
!
!  Redefine the interpolation point in coordinates relative to nearest grid
!  point and normalize with the cell size.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      dxp0=(xxp(1)-x(ix0))*dx_1(ix0)
      dyp0=(xxp(2)-y(iy0))*dy_1(iy0)
      dzp0=(xxp(3)-z(iz0))*dz_1(iz0)
!
!  Interpolation formulae.
!
      if (dimensionality==0) then
        gp=f(ix0,iy0,iz0,ivar1:ivar2)
      elseif (dimensionality==1) then
        if (nxgrid/=1) then
          gp = 0.5*(0.5-dxp0)**2*f(ix0-1,iy0,iz0,ivar1:ivar2) + &
                  (0.75-dxp0**2)*f(ix0  ,iy0,iz0,ivar1:ivar2) + &
               0.5*(0.5+dxp0)**2*f(ix0+1,iy0,iz0,ivar1:ivar2)
        endif
        if (nygrid/=1) then
          gp = 0.5*(0.5-dyp0)**2*f(ix0,iy0-1,iz0,ivar1:ivar2) + &
                  (0.75-dyp0**2)*f(ix0,iy0  ,iz0,ivar1:ivar2) + &
               0.5*(0.5+dyp0)**2*f(ix0,iy0+1,iz0,ivar1:ivar2)
        endif
        if (nzgrid/=1) then
          gp = 0.5*(0.5-dzp0)**2*f(ix0,iy0,iz0-1,ivar1:ivar2) + &
                  (0.75-dzp0**2)*f(ix0,iy0,iz0  ,ivar1:ivar2) + &
               0.5*(0.5+dzp0)**2*f(ix0,iy0,iz0+1,ivar1:ivar2)
        endif
      elseif (dimensionality==2) then
        if (nxgrid==1) then
          fac_y_m1 = 0.5*(0.5-dyp0)**2
          fac_y_00 = 0.75-dyp0**2
          fac_y_p1 = 0.5*(0.5+dyp0)**2
          fac_z_m1 = 0.5*(0.5-dzp0)**2
          fac_z_00 = 0.75-dzp0**2
          fac_z_p1 = 0.5*(0.5+dzp0)**2
!
          gp= fac_y_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_y_00*( f(ix0,iy0  ,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0  ,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_z_00*( f(ix0,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
              fac_y_p1*( f(ix0,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_y_m1*( f(ix0,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 )
        elseif (nygrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_z_m1 = 0.5*(0.5-dzp0)**2
          fac_z_00 = 0.75-dzp0**2
          fac_z_p1 = 0.5*(0.5+dzp0)**2
!
          gp= fac_x_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_x_00*( f(ix0  ,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0  ,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_z_00*( f(ix0+1,iy0,iz0  ,ivar1:ivar2)*fac_x_p1 + &
                         f(ix0-1,iy0,iz0  ,ivar1:ivar2)*fac_x_m1 ) + &
              fac_x_p1*( f(ix0+1,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0+1,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_x_m1*( f(ix0-1,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0-1,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 )
        elseif (nzgrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_y_m1 = 0.5*(0.5-dyp0)**2
          fac_y_00 = 0.75-dyp0**2
          fac_y_p1 = 0.5*(0.5+dyp0)**2
!
          gp= fac_x_00*fac_y_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_x_00*( f(ix0  ,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0  ,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 ) + &
              fac_y_00*( f(ix0+1,iy0  ,iz0,ivar1:ivar2)*fac_x_p1 + &
                         f(ix0-1,iy0  ,iz0,ivar1:ivar2)*fac_x_m1 ) + &
              fac_x_p1*( f(ix0+1,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0+1,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 ) + &
              fac_x_m1*( f(ix0-1,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0-1,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 )
        endif
      elseif (dimensionality==3) then
        fac_x_m1 = 0.5*(0.5-dxp0)**2
        fac_x_00 = 0.75-dxp0**2
        fac_x_p1 = 0.5*(0.5+dxp0)**2
        fac_y_m1 = 0.5*(0.5-dyp0)**2
        fac_y_00 = 0.75-dyp0**2
        fac_y_p1 = 0.5*(0.5+dyp0)**2
        fac_z_m1 = 0.5*(0.5-dzp0)**2
        fac_z_00 = 0.75-dzp0**2
        fac_z_p1 = 0.5*(0.5+dzp0)**2
!
        gp= fac_x_00*fac_y_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
            fac_x_00*fac_y_00*( f(ix0  ,iy0  ,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0  ,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_z_00*( f(ix0  ,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0  ,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
            fac_y_00*fac_z_00*( f(ix0+1,iy0  ,iz0  ,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0  ,ivar1:ivar2)*fac_x_m1 ) + &
            fac_x_p1*fac_y_p1*( f(ix0+1,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0+1,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_p1*fac_y_m1*( f(ix0+1,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0+1,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_m1*fac_y_p1*( f(ix0-1,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0-1,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_m1*fac_y_m1*( f(ix0-1,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0-1,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_y_p1*( f(ix0  ,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_y_m1*( f(ix0  ,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_y_00*fac_z_p1*( f(ix0+1,iy0  ,iz0+1,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0+1,ivar1:ivar2)*fac_x_m1 ) + &
            fac_y_00*fac_z_m1*( f(ix0+1,iy0  ,iz0-1,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0-1,ivar1:ivar2)*fac_x_m1 ) + &
            fac_z_00*fac_x_p1*( f(ix0+1,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0+1,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
            fac_z_00*fac_x_m1*( f(ix0-1,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0-1,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 )
      endif
!
      call keep_compiler_quiet(ipar,iblock)
!
    endsubroutine interpolate_quadratic_spline
!***********************************************************************
    subroutine interpolate_fourth(f,ivar1,ivar2,xxp,gp,inear,iblock,ipar)
!
!  Interpolate using 4th order polynomials in x and z. The weight function
!  is so complicated that a direct analytical solution gives little insight.
!  Instead we derive the weight function through the following steps:
!
!    f_i = M_ij*C_j
!
!  Here f_i is the value of the quantity at the grid cells, M_ij is the
!  interpolation matrix and C_j the interpolation parameters to be found.
!  We order the grid cells from (-2,-2) to (+2,+2) in x and z, incrementing
!  along x as the main direction. The coefficients are sorted as
!
!    C_{ 1- 5} = 1  , x    , x^2    , x^3     , x^4
!    C_{ 6-10} = z  , x*z  , x^2*z  , x^3*z  , x^4*z
!    C_{11-15} = z^2, x*z^2, x^2*z^2, x^3*z^2, x^4*z^2
!    C_{16-20} = z^3, x*z^3, x^2*z^3, x^3*z^3, x^4*z^3
!    C_{21-25} = z^4, x*z^4, x^2*z^4, x^3*z^4, x^4*z^4
!
!  The interpolated value is found through:
!
!    f = C_j*s_j = M^{-1}_ij*f_j*s_i = M^{-1}_ij*s_i*f_j = W_j*f_j
!
!  Here s is the separation vector s=(1,dx,dx^2,dx^3,dx^4,dz,dx*dz,...).
!
!  The weight function W = transpose(M^{-1}#s) can then be used directly for
!  the back-reaction as well.
!
!  18-oct-14/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (3) :: inear
      integer :: iblock, ipar
!
      real, dimension (25,25) :: invmat, invmatt
      real, dimension (25,ivar2-ivar1+1) :: gg
      real, dimension (25) :: sepvec, weight
      real :: dxp, dzp, dxp2, dzp2, dxp3, dzp3, dxp4, dzp4
      real, save :: dx1, dx2, dx3, dx4
      real, save :: dz1, dz2, dz3, dz4
      real, save :: dx1dz1, dx2dz1, dx3dz1, dx4dz1
      real, save :: dx1dz2, dx2dz2, dx3dz2, dx4dz2
      real, save :: dx1dz3, dx2dz3, dx3dz3, dx4dz3
      real, save :: dx1dz4, dx2dz4, dx3dz4, dx4dz4
      integer :: ix0, iy0, iz0, ipoint, ivar
      logical, save :: lfirstcall=.true.
!
      intent(in)  :: f, xxp, ivar1, inear
      intent(out) :: gp
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
!
!  Not implemented in y-direction yet (but is easy to generalise).
!
      if (nygrid/=1) then
        if (lroot) print*, 'interpolate_fourth: not implemented in y'
        call fatal_error('interpolate_quadratic','')
      endif
!
!  Precalculate values that only need to be calculated once.
!
      if (lfirstcall) then
        dx1=1/dx; dx2=1/dx**2; dx3=1/dx**3; dx4=1/dx**4
        dz1=1/dz; dz2=1/dz**2; dz3=1/dz**3; dz4=1/dz**4
        dx1dz1=dx1*dz1; dx2dz1=dx2*dz1; dx3dz1=dx3*dz1; dx4dz1=dx4*dz1
        dx1dz2=dx1*dz2; dx2dz2=dx2*dz2; dx3dz2=dx3*dz2; dx4dz2=dx4*dz2
        dx1dz3=dx1*dz3; dx2dz3=dx2*dz3; dx3dz3=dx3*dz3; dx4dz3=dx4*dz3
        dx1dz4=dx1*dz4; dx2dz4=dx2*dz4; dx3dz4=dx3*dz4; dx4dz4=dx4*dz4
!
!  The inverse of the interpolation matrix is input analytically. Note that
!  zeros here denote zeros in the actual matrix, to avoid infinities.
!
        invmat(1,:)=(/0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0/)
        invmat(2,:)=(/0,0,0,0,0,0,0,0,0,0,12,-1,0,1,-11,0,0,0,0,0,0,0,0,0,0/)
        invmat(3,:)=(/0,0,0,0,0,0,0,0,0,0,-24,1,0,1,-24,0,0,0,0,0,0,0,0,0,0/)
        invmat(4,:)=(/0,0,0,0,0,0,0,0,0,0,-12,6,0,-6,11,0,0,0,0,0,0,0,0,0,0/)
        invmat(5,:)=(/0,0,0,0,0,0,0,0,0,0,24,-6,4,-6,24,0,0,0,0,0,0,0,0,0,0/)
        invmat(6,:)=(/0,0,12,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-11,0,0/)
        invmat(7,:)=(/144,-18,0,18,-143,-18,2,0,-2,18,0,0,0,0,0,17,-2,0,2,-17,-144, 18,0,-17,143/)
        invmat(8,:)=(/-287,18,-9,18,-287,36,-2,1,-2,36,0,0,0,0,0,-35,2,-1,2,-35,287, -18,9,-18,287/)
        invmat(9,:)=(/-144,72,0,-72,144,18,-9,0,9,-18,0,0,0,0,0,-18,9,0,-9,18,144, -72,0,71,-143/)
        invmat(10,:)=(/288,-72,48,-72,288,-36,9,-6,9,-36,0,0,0,0,0,36,-9,6,-9,36, -288,72,-47,72,-288/)
        invmat(11,:)=(/0,0,-24,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-24,0,0/)
        invmat(12,:)=(/-287,36,0,-35,287,18,-2,0,2,-18,-9,1,0,-1,9,18,-2,0,2,-18, -287,36,0,-35,287/)
        invmat(13,:)=(/575,-36,19,-36,575,-36,2,-1,2, -36,19,-1,0,-1,19,-36,2,-1,2,-36,575,-36,19,-36,575/)
        invmat(14,:)=(/288,-144,0,143,-288,-18,9,0,-9,18,9,-4,0,4,-9,-18,9,0,-9,18,288, -144,0,143,-288/)
        invmat(15,:)=(/-576,144,-96,144,-576,36,-9,6,-9, 36,-19,4,-3,4,-19,36,-9,6,-9,36,-576,144,-96,144,-576/)
        invmat(16,:)=(/0,0,-12,0,0,0,0,6,0,0,0,0,0,0,0,0,0,-6,0,0,0,0,11,0,0/)
        invmat(17,:)=(/-144,18,0,-18,144,72,-9,0,9,-72,0,0,0,0,0,-71,9,0,-9,71,144, -18,0,18,-143/)
        invmat(18,:)=(/288,-18,9,-18,288,-144,9,-4,9,-144,0,0,0,0,0,143,-9,4,-9,143, -288,18,-9,18,-288/)
        invmat(19,:)=(/144,-72,0,72,-144,-72,36,0,-36,72,0,0,0,0,0,72,-36,0,36,-72,-144, 72,0,-72,143/)
        invmat(20,:)=(/-288,72,-48,72,-288,144,-36,24,-36,144,0,0,0,0,0,-144,36,-24,36,-144,288,-72,47,-72,288/)
        invmat(21,:)=(/0,0,24,0,0,0,0,-6,0,0,0,0,4,0,0,0,0,-6,0,0,0,0,24,0,0/)
        invmat(22,:)=(/288,-36,0,36,-288,-72,9,0,-9,72,48,-6,0,6,-47,-72,9,0,-9,72,288, -36,0,36,-288/)
        invmat(23,:)=(/-576,36,-19,36,-576,144,-9,4,-9,144,-96,6,-3,6,-96,144,-9,4,-9,144,-576,36,-19,36,-576/)
        invmat(24,:)=(/-288,144,0,-144,288,72,-36,0,36,-72,-48,24,0,-24,47,72,-36,0,36,-72,-288,144,0,-144,288/)
        invmat(25,:)=(/576,-144,96,-144,576,-144,36,-24,36,-144,96,-24,16,-24,96,-144,36,-24,36,-144,576,-144,96,-144,576/)
        where (invmat/=0.0) invmat=1/invmat
        invmatt=transpose(invmat)
      endif
!
!  Define function values at the grid points.
!
      ipoint=1
      do iz=-2,+2; do ix=-2,+2
        gg(ipoint,:)=f(ix0+ix,iy0,iz0+iz,ivar1:ivar2)
        ipoint=ipoint+1
      enddo; enddo
!
!  Calculate the separation of the interpolation point from the central grid
!  cell (and various powers of it).
!
      dxp=xxp(1)-x(ix0)
      dzp=xxp(3)-z(iz0)
      dxp2=dxp**2
      dzp2=dzp**2
      dxp3=dxp**3
      dzp3=dzp**3
      dxp4=dxp**4
      dzp4=dzp**4
!
!  Calculate the elements of the separation vector (see discussion above).
!
      sepvec( 1)=1.0
      sepvec( 2)=dxp*dx1
      sepvec( 3)=dxp2*dx2
      sepvec( 4)=dxp3*dx3
      sepvec( 5)=dxp4*dx4
      sepvec( 6)=dzp*dz1
      sepvec( 7)=dxp*dzp*dx1dz1
      sepvec( 8)=dxp2*dzp*dx2dz1
      sepvec( 9)=dxp3*dzp*dx3dz1
      sepvec(10)=dxp4*dzp*dx4dz1
      sepvec(11)=dzp2*dz2
      sepvec(12)=dxp*dzp2*dx1dz2
      sepvec(13)=dxp2*dzp2*dx2dz2
      sepvec(14)=dxp3*dzp2*dx3dz2
      sepvec(15)=dxp4*dzp2*dx4dz2
      sepvec(16)=dzp3*dz3
      sepvec(17)=dxp*dzp3*dx1dz3
      sepvec(18)=dxp2*dzp3*dx2dz3
      sepvec(19)=dxp3*dzp3*dx3dz3
      sepvec(20)=dxp4*dzp3*dx4dz3
      sepvec(21)=dzp4*dz4
      sepvec(22)=dxp*dzp4*dx1dz4
      sepvec(23)=dxp2*dzp4*dx2dz4
      sepvec(24)=dxp3*dzp4*dx3dz4
      sepvec(25)=dxp4*dzp4*dx4dz4
!
!  Calculate the weight vector.
!
      weight=matmul(invmatt,sepvec)
!
!  Finally use the weight vector to calculate the interpolated value as a
!  simple sum over the involved grid cells.
!
      do ivar=1,ivar2-ivar1+1
        gp(ivar) = sum(weight*gg(:,ivar))
      enddo
!
      call keep_compiler_quiet(ipar,iblock)
!
    endsubroutine interpolate_fourth
!***********************************************************************
    subroutine map_nearest_grid(fp,ineargrid,k1_opt,k2_opt)
!
!  Find index (ix0, iy0, iz0) of nearest grid point of all or some of the
!  particles.
!
!  23-jan-05/anders: coded
!  08-jul-08/kapelrud: support for non-equidistant grids
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, optional :: k1_opt, k2_opt
!
      double precision, save :: dx1, dy1, dz1
      integer :: k, k1, k2, ix0, iy0, iz0
      logical, save :: lfirstcall=.true.
      real :: t_sp   ! t in single precision for backwards compatibility
!
      intent(in)  :: fp
      intent(out) :: ineargrid
!
      t_sp = t
!
      if (present(k1_opt)) then
        k1=k1_opt
      else
        k1=1
      endif
!
      if (present(k2_opt)) then
        k2=k2_opt
      else
        k2=npar_loc
      endif
!
!  Default values in case of missing directions.
!
      ix0=nghost+1; iy0=nghost+1; iz0=nghost+1
!
      if (lfirstcall) then
        dx1=dx_1(l1); dy1=dy_1(m1); dz1=dz_1(n1)
        lfirstcall=.false.
      endif
!
      do k=k1,k2
!
!  Find nearest grid point in x-direction.
!  Find nearest grid point by bisection if the grid is not equidistant.
!
        if (nxgrid/=1) then
          if (lequidist(1)) then
            ix0 = nint((fp(k,ixp)-x(1))*dx1) + 1
          else
            call find_index_by_bisection(fp(k,ixp),x,ix0)
          endif
        endif
!
!  Find nearest grid point in y-direction.
!
        if (nygrid/=1) then
          if (lequidist(2)) then
            iy0 = nint((fp(k,iyp)-y(1))*dy1) + 1
          else
            call find_index_by_bisection(fp(k,iyp),y,iy0)
          endif
        endif
!
!  Find nearest grid point in z-direction.
!
        if (nzgrid/=1) then
          if (lequidist(3)) then
            iz0 = nint((fp(k,izp)-z(1))*dz1) + 1
          else
            call find_index_by_bisection(fp(k,izp),z,iz0)
          endif
        endif
!
        ineargrid(k,1)=ix0; ineargrid(k,2)=iy0; ineargrid(k,3)=iz0
!
!  Round off errors may put a particle closer to a ghost point than to a
!  physical point. Either stop the code with a fatal error or fix problem
!  by forcing the nearest grid point to be a physical point.
!
          if (ineargrid(k,1)<=l1-1.or.ineargrid(k,1)>=l2+1.or. &
              ineargrid(k,2)<=m1-1.or.ineargrid(k,2)>=m2+1.or. &
              ineargrid(k,3)<=n1-1.or.ineargrid(k,3)>=n2+1) then
            if (lcheck_exact_frontier) then
              if (ineargrid(k,1)<=l1-1) then
                ineargrid(k,1)=l1
              elseif (ineargrid(k,1)>=l2+1) then
                ineargrid(k,1)=l2
              endif
              if (ineargrid(k,2)<=m1-1) then
                ineargrid(k,2)=m1
              elseif (ineargrid(k,2)>=m2+1) then
                ineargrid(k,2)=m2
              endif
              if (ineargrid(k,3)<=n1-1) then
                ineargrid(k,3)=n1
              elseif (ineargrid(k,3)>=n2+1) then
                ineargrid(k,3)=n2
              endif
            else
              print*, 'map_nearest_grid: particle must never be closer to a '//&
                      'ghost point than'
              print*, '                  to a physical point.'
              print*, '                  Consider using double precision to '//&
                      'avoid this problem'
              print*, '                  or set lcheck_exact_frontier in '// &
                      '&particles_run_pars.'
              print*, 'Information about what went wrong:'
              print*, '----------------------------------'
              print*, 'it, itsub, t=', it, itsub, t_sp
              print*, 'iproc  =', iproc
              print*, 'ipar, k     =', ipar(k), k
              print*, 'xxp         =', fp(k,ixp:izp)
              if (ivpx/=0) print*, 'vvp         =', fp(k,ivpx:ivpz)
              print*, 'ineargrid   =', ineargrid(k,:)
              print*, 'l1, m1, n1  =', l1, m1, n1
              print*, 'l2, m2, n2  =', l2, m2, n2
              print*, 'x1, y1, z1  =', x(l1), y(m1), z(n1)
              print*, 'x2, y2, z2  =', x(l2), y(m2), z(n2)
              call fatal_error_local('map_nearest_grid','')
            endif
        endif
      enddo
!
    endsubroutine map_nearest_grid
!***********************************************************************
    subroutine sort_particles_imn(fp,ineargrid,ipar,dfp,f)
!
!  Sort the particles so that they appear in the same order as the (m,n) loop.
!
!  20-apr-06/anders: coded
!
      use General, only: safe_character_assign
!
      real, dimension (mx,my,mz,mfarray),optional :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
!
      real :: t_sp   ! t in single precision for backwards compatibility
      real, dimension (mparray) :: fp_tmp
      real, dimension (mpvar) :: dfp_tmp
      integer, dimension (3) :: ineargrid_tmp
      integer, dimension (ny*nz) :: kk
      integer :: ilmn_par_tmp, ipark_sorted_tmp, ipar_tmp
      integer, dimension (mpar_loc) :: ilmn_par, ipark_sorted
      integer :: j, k, ix0, iy0, iz0, ih, lun
      integer, save :: ncount=-1, isorttype=3
      integer, parameter :: nshellsort=21
      integer, dimension (nshellsort), parameter :: &
          hshellsort=(/ 14057, 9371, 6247, 4177, 2777, 1861, 1237, 823, 557, &
                        367, 251, 163, 109, 73, 37, 19, 11, 7, 5, 3, 1 /)
      logical, save :: lrunningsort=.false.
      character (len=fnlen) :: filename
!
      intent(inout)  :: fp, ineargrid, ipar, dfp
!
      t_sp = t
!
!  Determine beginning and ending index of particles in pencil (m,n).
!
      call particle_pencil_index(ineargrid)
!
!  Choose sort algorithm.
!  WARNING - choosing the wrong one might make the code unnecessarily slow.
!
      if ( lstart .or. &
           (lshear .and. Omega/=0.0) .and. (nxgrid>1 .and. nygrid>1) ) then
        isorttype=3
        lrunningsort=.false.
      else
        isorttype=1
        lrunningsort=.true.
      endif
      ncount=0
!
!  Calculate integer value to sort after.
!
      do k=1,npar_loc
        ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
        ilmn_par(k)=imn_array(iy0,iz0)!-1)*ny*nz+ix0
        ipark_sorted(k)=k
      enddo
!
!  Sort using either straight insertion (1), shell sorting (2) or counting
!  sort (3).
!
      select case (isorttype)
!  Straight insertion.
      case (1)
        do k=2,npar_loc
!
          j=k
!
          do while ( ilmn_par(k)<ilmn_par(j-1) )
            j=j-1
            if (j==1) exit
          enddo
!
          if (j/=k) then
            ncount=ncount+k-j
!
            ilmn_par_tmp=ilmn_par(k)
            ilmn_par(j+1:k)=ilmn_par(j:k-1)
            ilmn_par(j)=ilmn_par_tmp
            ipark_sorted_tmp=ipark_sorted(k)
            ipark_sorted(j+1:k)=ipark_sorted(j:k-1)
            ipark_sorted(j)=ipark_sorted_tmp
!  Sort particle data on the fly (practical for almost ordered arrays).
            if (lrunningsort) then
              fp_tmp=fp(k,:)
              fp(j+1:k,:)=fp(j:k-1,:)
              fp(j,:)=fp_tmp
!
              if (present(dfp)) then
                dfp_tmp=dfp(k,:)
                dfp(j+1:k,:)=dfp(j:k-1,:)
                dfp(j,:)=dfp_tmp
              endif
!
              ineargrid_tmp=ineargrid(k,:)
              ineargrid(j+1:k,:)=ineargrid(j:k-1,:)
              ineargrid(j,:)=ineargrid_tmp
!
              ipar_tmp=ipar(k)
              ipar(j+1:k)=ipar(j:k-1)
              ipar(j)=ipar_tmp
!
            endif
          endif
        enddo
!  Shell sort.
      case (2)
!
        do ih=1,21
          do k=1+hshellsort(ih),npar_loc
            ilmn_par_tmp=ilmn_par(k)
            ipark_sorted_tmp=ipark_sorted(k)
            j=k
            do while (ilmn_par(j-hshellsort(ih)) > ilmn_par_tmp)
              ncount=ncount+1
              ilmn_par(j)=ilmn_par(j-hshellsort(ih))
              ilmn_par(j-hshellsort(ih))=ilmn_par_tmp
              ipark_sorted(j)=ipark_sorted(j-hshellsort(ih))
              ipark_sorted(j-hshellsort(ih))=ipark_sorted_tmp
              j=j-hshellsort(ih)
              if (j-hshellsort(ih)<1) exit
            enddo
          enddo
        enddo
!  Counting sort.
      case (3)
        kk=k1_imn
        do k=1,npar_loc
          ipark_sorted(kk(ilmn_par(k)))=k
          kk(ilmn_par(k))=kk(ilmn_par(k))+1
        enddo
        ncount=npar_loc
      endselect
!
!  Sort particle data according to sorting index.
!
      if (lrunningsort .and. isorttype/=1) then
        if (lroot) print*, 'sort_particles_imn: lrunningsort is only '// &
             'allowed for straight insertion sort.'
        call fatal_error('sort_particles_imn','')
      endif
!
      if ( (.not. lrunningsort) .and. (ncount/=0) ) then
        fp(1:npar_loc,:)=fp(ipark_sorted(1:npar_loc),:)
        if (present(dfp)) dfp(1:npar_loc,:)=dfp(ipark_sorted(1:npar_loc),:)
        ineargrid(1:npar_loc,:)=ineargrid(ipark_sorted(1:npar_loc),:)
        ipar(1:npar_loc)=ipar(ipark_sorted(1:npar_loc))
      endif
!
      if (lroot.and.ldiagnos) then
        call safe_character_assign(filename,trim(datadir)//'/sort_particles.dat')
        lun=1
        open(lun,file=trim(filename),action='write',position='append')
        write(lun,'(A15,f7.3)') '------------ t=', t_sp
        write(lun,'(A40,3i9,l9)')  'iproc, ncount, isorttype, lrunningsort=', &
             iproc, ncount, isorttype, lrunningsort
        close (lun)
      endif
!
      if (ip<=8) print '(A,i4,i8,i4,l4)', &
           'sort_particles_imn: iproc, ncount, isorttype, lrunningsort=', &
           iproc, ncount, isorttype, lrunningsort
!
!  Possible to randomize particles inside each pencil. This screws with the
!  pencil consistency check, so we turn it off when the test is running.
!
      if (lrandom_particle_pencils .and. (.not.lpencil_check_at_work)) then
        if (present(dfp)) then
          call random_particle_pencils(fp,ineargrid,ipar,dfp)
        else
          call random_particle_pencils(fp,ineargrid,ipar)
        endif
      endif
!
    endsubroutine sort_particles_imn
!***********************************************************************
    subroutine boundcond_neighbour_list
!
! Copy the number of neighbours to the boundary points of the 
! neighbour list
!
!      nlist(1-neighbourx:0) = nlist()
!      allocate(nlist(1-neighbourx:nx+neighbourx,1-neighboury:nx+neighboury, &
!        1-neighbourz:nx+neighbourz,Nneighbour+1))

!
    endsubroutine boundcond_neighbour_list
!***********************************************************************
    subroutine invert_ineargrid_list(fp,ineargrid,ipar,dfp,f)
!
!  Update the neighbour list for all the particles
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
!
      integer :: k,ix0,ix,im,in,imn
      integer, dimension(nx) :: pinl, pin_cell
      integer, dimension(0:nx-1) :: cuml_pinl
      integer, dimension(npar_loc) :: lpark_sorted
!
! Go through the whole array fp to distribute them in cells
!
      f(:,:,:,iinvgrid:iinvgrid+1) = 0
      do imn=1,ny*nz
        im=nn(imn)
        in=mm(imn)
        if (npar_imn(imn)>=2) then
 !         write(*,*) 'DHRUBA, in here ? ',imn
          pinl = 0
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            write(*,*) 'ix0 and k!',ix0,k
            pinl(ix0)=pinl(ix0)+1
          enddo
          cuml_pinl(0)=0
          do ix=1,nx-1
            cuml_pinl(ix)=cuml_pinl(ix-1)+pinl(ix)
          enddo
 !         write(*,*)'DHRUBA, cuml',cuml_pinl
          pin_cell = 0
          do k=k1_imn(imn),k2_imn(imn)
            ix0=ineargrid(k,1)
            lpark_sorted(k-k1_imn(imn)+1) = k1_imn(imn)+cuml_pinl(ix0-1)+pin_cell(ix0)
            pin_cell(ix0) = pin_cell(ix0)+1
          enddo
          fp(k1_imn(imn):k2_imn(imn),:)=fp(lpark_sorted(1:npar_imn(imn)),:)
          if (present(dfp)) &
            dfp(k1_imn(imn):k2_imn(imn),:)=dfp(lpark_sorted(1:npar_imn(imn)),:)
          ineargrid(k1_imn(imn):k2_imn(imn),:)=ineargrid(lpark_sorted(1:npar_imn(imn)),:)
          ipar(k1_imn(imn):k2_imn(imn))=ipar(lpark_sorted(1:npar_imn(imn)))
          do ix=1,nx
            if (pin_cell(ix) .gt. 0) then
              f(ix+nghost,im,in,iinvgrid)= k1_imn(imn)+cuml_pinl(ix-1)
              f(ix+nghost,im,in,iinvgrid+1)=k1_imn(imn)+cuml_pinl(ix-1)+pin_cell(ix)-1
!              write(*,*) cuml_pinl(ix-1),k1_imn(imn),k2_imn(imn),f(ix+nghost,im,in,iinvgrid),f(ix+nghost,im,in,iinvgrid+1)
            else
              f(ix+nghost,im,in,iinvgrid)= 0
              f(ix+nghost,im,in,iinvgrid+1)= 0
            endif
          enddo
        elseif (npar_imn(imn).eq.1) then
          k=k1_imn(imn)
          ix0=ineargrid(k,1)
          f(ix0+nghost,iy,iz,iinvgrid)=k
          f(ix0+nghost,iy,iz,iinvgrid+1)=k
        endif
      enddo   !loop over imn ends
!      do ix=l1,l2;do iy=m1,m2;do iz=n1,n2;
!        write(*,*) 'ix,iy,iz,iinvgrid,f',ix,iy,iz,iinvgrid,f(ix,iy,iz,iinvgrid)
        !if(f(ix,iy,iz,iinvgrid+1).ne.0.) write(*,*) 'ix,iy,iz,iinvgrid,f',ix,iy,iz,iinvgrid+1,f(ix,iy,iz,iinvgrid+1)
!      enddo;enddo;enddo
!
    endsubroutine invert_ineargrid_list
!***********************************************************************
    subroutine random_particle_pencils(fp,ineargrid,ipar,dfp)
!
!  Randomize particles within each pencil to avoid low index particles
!  always being considered first.
!
!  Slows down simulation by around 10%.
!
!  12-nov-09/anders: coded
!
      use General, only: random_number_wrapper
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
!
      real, dimension (mparray) :: fp_swap
      real, dimension (mpvar) :: dfp_swap
      real :: r
      integer, dimension (3) :: ineargrid_swap
      integer :: ipar_swap, imn, k, kswap
!
      intent (inout) :: fp, ineargrid, ipar, dfp
!
      do imn=1,ny*nz
        if (npar_imn(imn)>=2) then
          do k=k1_imn(imn),k2_imn(imn)
            call random_number_wrapper(r)
            kswap=k1_imn(imn)+floor(r*npar_imn(imn))
            if (kswap/=k) then
              fp_swap=fp(kswap,:)
              ineargrid_swap=ineargrid(kswap,:)
              ipar_swap=ipar(kswap)
              if (present(dfp)) dfp_swap=dfp(kswap,:)
              fp(kswap,:)=fp(k,:)
              ineargrid(kswap,:)=ineargrid(k,:)
              ipar(kswap)=ipar(k)
              if (present(dfp)) dfp(kswap,:)=dfp(k,:)
              fp(k,:)=fp_swap
              ineargrid(k,:)=ineargrid_swap
              ipar(k)=ipar_swap
              if (present(dfp)) dfp(k,:)=dfp_swap
            endif
          enddo
        endif
      enddo
!
    endsubroutine random_particle_pencils
!***********************************************************************
    subroutine map_xxp_grid(f,fp,ineargrid,lmapsink_opt)
!
!  Map the particles as a continuous density field on the grid.
!
!  27-nov-05/anders: coded
!
      use GhostFold,     only: fold_f
      use Particles_sub, only: get_rhopswarm
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      logical, intent(in), optional :: lmapsink_opt
!
      real, dimension(nx) :: rhop_swarm_mn
      real :: weight0, weight, weight_x, weight_y, weight_z
      integer :: ik, k, ix0, iy0, iz0, ixx, iyy, izz
      integer :: ixx0, ixx1, iyy0, iyy1, izz0, izz1, irhopm, iapn_ik
      logical :: lsink, lmapsink, max_iapn_use
!
!  Possible to map sink particles by temporarily switching irhop to irhops.
!
      if (present(lmapsink_opt)) then
        lmapsink=lmapsink_opt
        if (lmapsink) then
          irhopm=irhop
          irhop=ipotself
        endif
      else
        lmapsink=.false.
      endif
!
!  Calculate the number of particles in each grid cell.
!
      if (inp/=0 .and. (.not. lnocalc_np) .and. (.not. lmapsink)) then
        f(:,:,:,inp)=0.0
        do k=1,npar_loc
          ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
          f(ix0,iy0,iz0,inp) = f(ix0,iy0,iz0,inp) + 1.0
        enddo
      endif
!
!  Calculate the number of particles in each grid cell, for each particle size.
!
      max_iapn_use=.false.
      loop_ap: do ik=1,ndustrad
        if ( iapn(ik) /=0 ) then
          f(:,:,:,iapn(ik))=0.0
          max_iapn_use=.true.
        else
          max_iapn_use=.false.
        endif
      enddo loop_ap
!
      calculate__np_ap_spec: if ( max_iapn_use ) then
        if (lparticlemesh_cic) then
!
!  Cloud In Cell (CIC) scheme.
!
          do k=1,npar_loc
!
            call map_xxp_grid__loop_ap(fp,k,iapn_ik)
!
            ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
            ixx0=ix0; iyy0=iy0; izz0=iz0
            ixx1=ix0; iyy1=iy0; izz1=iz0
            if ( (x(ix0)>fp(k,ixp)) .and. nxgrid/=1) ixx0=ixx0-1
            if ( (y(iy0)>fp(k,iyp)) .and. nygrid/=1) iyy0=iyy0-1
            if ( (z(iz0)>fp(k,izp)) .and. nzgrid/=1) izz0=izz0-1
            if (nxgrid/=1) ixx1=ixx0+1
            if (nygrid/=1) iyy1=iyy0+1
            if (nzgrid/=1) izz1=izz0+1
!
            weight0=1.0
!
            do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
!
              weight=weight0
              if (nxgrid/=1) then
                weight_x = 1.0-abs(fp(k,ixp)-x(ixx))*dx_1(ixx)
                weight=weight*weight_x
              endif
              if (nygrid/=1) then
                weight_y = 1.0-abs(fp(k,iyp)-y(iyy))*dy_1(iyy)
                weight=weight*weight_y
              endif
              if (nzgrid/=1) then
                weight_z = 1.0-abs(fp(k,izp)-z(izz))*dz_1(izz)
                weight=weight*weight_z
              endif
!
              f(ixx,iyy,izz,iapn_ik)=f(ixx,iyy,izz,iapn_ik) + weight
!
            enddo; enddo; enddo
          enddo
!
!  Triangular Shaped Cloud (TSC) scheme.
!
        elseif (lparticlemesh_tsc) then
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
          do k=1,npar_loc
!
            call map_xxp_grid__loop_ap(fp,k,iapn_ik)
!
            ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
            if (nxgrid/=1) then
              ixx0=ix0-1; ixx1=ix0+1
            else
              ixx0=ix0  ; ixx1=ix0
            endif
            if (nygrid/=1) then
              iyy0=iy0-1; iyy1=iy0+1
            else
              iyy0=iy0  ; iyy1=iy0
            endif
            if (nzgrid/=1) then
              izz0=iz0-1; izz1=iz0+1
            else
              izz0=iz0  ; izz1=iz0
            endif
!
            weight0=1.0
!
!  The nearest grid point is influenced differently than the left and right
!  neighbours are. A particle that is situated exactly on a grid point gives
!  3/4 contribution to that grid point and 1/8 to each of the neighbours.
!
            do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
              if ( ((ixx-ix0)==-1) .or. ((ixx-ix0)==+1) ) then
                weight_x = 1.125 - 1.5* abs(fp(k,ixp)-x(ixx))*dx_1(ixx) + &
                                   0.5*(abs(fp(k,ixp)-x(ixx))*dx_1(ixx))**2
              else
                if (nxgrid/=1) &
                    weight_x = 0.75  -   ((fp(k,ixp)-x(ixx))*dx_1(ixx))**2
              endif
              if ( ((iyy-iy0)==-1) .or. ((iyy-iy0)==+1) ) then
                weight_y = 1.125 - 1.5* abs(fp(k,iyp)-y(iyy))*dy_1(iyy) + &
                                   0.5*(abs(fp(k,iyp)-y(iyy))*dy_1(iyy))**2
              else
                if (nygrid/=1) &
                    weight_y = 0.75  -   ((fp(k,iyp)-y(iyy))*dy_1(iyy))**2
              endif
              if ( ((izz-iz0)==-1) .or. ((izz-iz0)==+1) ) then
                weight_z = 1.125 - 1.5* abs(fp(k,izp)-z(izz))*dz_1(izz) + &
                                   0.5*(abs(fp(k,izp)-z(izz))*dz_1(izz))**2
              else
                if (nzgrid/=1) &
                    weight_z = 0.75  -   ((fp(k,izp)-z(izz))*dz_1(izz))**2
              endif
!
              weight=weight0
!
              if (nxgrid/=1) weight=weight*weight_x
              if (nygrid/=1) weight=weight*weight_y
              if (nzgrid/=1) weight=weight*weight_z
!
              f(ixx,iyy,izz,iapn_ik)=f(ixx,iyy,izz,iapn_ik) + weight
!
            enddo; enddo; enddo
          enddo
!
!  Nearest Grid Point (NGP) method.
!
        else
          do k=1,npar_loc
!
            call map_xxp_grid__loop_ap(fp,k,iapn_ik)
!
            ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
!
            f(ix0,iy0,iz0,iapn_ik)=f(ix0,iy0,iz0,iapn_ik) + 1.0
!
          enddo
        endif
!
!  Fold first ghost zone of f.
!
        if (lparticlemesh_cic.or.lparticlemesh_tsc) then
          do ik=1,ndustrad
            call fold_f(f,iapn(ik),iapn(ik))
          enddo
        endif
      end if calculate__np_ap_spec
!
!  Calculate the smooth number of particles in each grid cell. Three methods are
!  implemented for assigning a particle to the mesh (see Hockney & Eastwood):
!
!    0. NGP (Nearest Grid Point)
!       The entire effect of the particle goes to the nearest grid point.
!    1. C (Cloud In Cell)
!       The particle has a region of influence with the size of a grid cell.
!       This is equivalent to a first order (spline) interpolation scheme.
!    2. TSC (Triangular Shaped Cloud)
!       The particle is spread over a length of two grid cells, but with
!       a density that falls linearly outwards.
!       This is equivalent to a second order spline interpolation scheme.
!
      if ( (irhop/=0 .and. (.not. lnocalc_rhop)) .or. lmapsink .or. ipeh/=0) then
        f(:,:,:,irhop)=0.0
        if (ipeh/=0) f(:,:,:,ipeh)=0.0
        if (lparticlemesh_cic) then
!
!  Cloud In Cell (CIC) scheme.
!
          do k=1,npar_loc
            lsink=.false.
            if (lparticles_sink) then
              if (fp(k,iaps)>0.0) lsink=.true.
            endif
            if (lmapsink) lsink=.not.lsink
            if (.not.lsink) then
              ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
              ixx0=ix0; iyy0=iy0; izz0=iz0
              ixx1=ix0; iyy1=iy0; izz1=iz0
              if ( (x(ix0)>fp(k,ixp)) .and. nxgrid/=1) ixx0=ixx0-1
              if ( (y(iy0)>fp(k,iyp)) .and. nygrid/=1) iyy0=iyy0-1
              if ( (z(iz0)>fp(k,izp)) .and. nzgrid/=1) izz0=izz0-1
              if (nxgrid/=1) ixx1=ixx0+1
              if (nygrid/=1) iyy1=iyy0+1
              if (nzgrid/=1) izz1=izz0+1
!
              if (lparticles_density) then
                weight0=fp(k,irhopswarm)
              elseif (lparticles_radius.and.lparticles_number) then
                weight0=four_pi_rhopmat_over_three*fp(k,iap)**3*fp(k,inpswarm)
              elseif (lparticles_radius) then
                weight0=four_pi_rhopmat_over_three*fp(k,iap)**3*np_swarm
              elseif (lparticles_number) then
                weight0=mpmat*fp(k,inpswarm)
              else
                weight0=1.0
              endif
!
              do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
!
                weight=weight0
                if (nxgrid/=1) then
                  weight_x = 1.0-abs(fp(k,ixp)-x(ixx))*dx_1(ixx)
                  weight=weight*weight_x
                endif
                if (nygrid/=1) then
                  weight_y = 1.0-abs(fp(k,iyp)-y(iyy))*dy_1(iyy)
                  weight=weight*weight_y
                endif
                if (nzgrid/=1) then
                  weight_z = 1.0-abs(fp(k,izp)-z(izz))*dz_1(izz)
                  weight=weight*weight_z
                endif
!
                f(ixx,iyy,izz,irhop)=f(ixx,iyy,izz,irhop) + weight
                if (lparticles_lyapunov) then
                   if (ibbf/=0) then
                        f(ixx,iyy,izz,ibxf:ibzf)=f(ixx,iyy,izz,ibxf:ibzf)+weight*fp(k,ibpx:ibpz) 
                   endif
                endif
!
                if (ipeh/=0 .and. (lcylindrical_coords .or. lspherical_coords)) then
                  if (lparticles_number) then
                    weight = fp(k,inpswarm)*(fp(k,iap)**2.0)
                  else
                    weight = np_swarm*(fp(k,iap)**2.0)
                  endif
                  if (nxgrid/=1) weight=weight*weight_x
                  if (nygrid/=1) weight=weight*weight_y
                  if (nzgrid/=1) weight=weight*weight_z
                  f(ixx,iyy,izz,ipeh) = f(ixx,iyy,izz,ipeh)+weight
                endif
!
              enddo; enddo; enddo
            endif
          enddo
!
!  Triangular Shaped Cloud (TSC) scheme.
!
        elseif (lparticlemesh_tsc) then
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
          do k=1,npar_loc
            lsink=.false.
            if (lparticles_sink) then
              if (fp(k,iaps)>0.0) lsink=.true.
            endif
            if (lmapsink) lsink=.not.lsink
            if (.not.lsink) then
              ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
              if (nxgrid/=1) then
                ixx0=ix0-1; ixx1=ix0+1
              else
                ixx0=ix0  ; ixx1=ix0
              endif
              if (nygrid/=1) then
                iyy0=iy0-1; iyy1=iy0+1
              else
                iyy0=iy0  ; iyy1=iy0
              endif
              if (nzgrid/=1) then
                izz0=iz0-1; izz1=iz0+1
              else
                izz0=iz0  ; izz1=iz0
              endif
!
              if (lparticles_density) then
                weight0=fp(k,irhopswarm)
              elseif (lparticles_radius.and.lparticles_number) then
                weight0=four_pi_rhopmat_over_three*fp(k,iap)**3*fp(k,inpswarm)
              elseif (lparticles_radius) then
                weight0=four_pi_rhopmat_over_three*fp(k,iap)**3*np_swarm
              elseif (lparticles_number) then
                weight0=mpmat*fp(k,inpswarm)
              else
                weight0=1.0
              endif
!
!  The nearest grid point is influenced differently than the left and right
!  neighbours are. A particle that is situated exactly on a grid point gives
!  3/4 contribution to that grid point and 1/8 to each of the neighbours.
!
              do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
                if ( ((ixx-ix0)==-1) .or. ((ixx-ix0)==+1) ) then
                  weight_x = 1.125 - 1.5* abs(fp(k,ixp)-x(ixx))*dx_1(ixx) + &
                                     0.5*(abs(fp(k,ixp)-x(ixx))*dx_1(ixx))**2
                else
                  if (nxgrid/=1) &
                      weight_x = 0.75  -   ((fp(k,ixp)-x(ixx))*dx_1(ixx))**2
                endif
                if ( ((iyy-iy0)==-1) .or. ((iyy-iy0)==+1) ) then
                  weight_y = 1.125 - 1.5* abs(fp(k,iyp)-y(iyy))*dy_1(iyy) + &
                                     0.5*(abs(fp(k,iyp)-y(iyy))*dy_1(iyy))**2
                else
                  if (nygrid/=1) &
                      weight_y = 0.75  -   ((fp(k,iyp)-y(iyy))*dy_1(iyy))**2
                endif
                if ( ((izz-iz0)==-1) .or. ((izz-iz0)==+1) ) then
                  weight_z = 1.125 - 1.5* abs(fp(k,izp)-z(izz))*dz_1(izz) + &
                                     0.5*(abs(fp(k,izp)-z(izz))*dz_1(izz))**2
                else
                  if (nzgrid/=1) &
                      weight_z = 0.75  -   ((fp(k,izp)-z(izz))*dz_1(izz))**2
                endif
!
                weight=weight0
!
                if (nxgrid/=1) weight=weight*weight_x
                if (nygrid/=1) weight=weight*weight_y
                if (nzgrid/=1) weight=weight*weight_z
!
                f(ixx,iyy,izz,irhop)=f(ixx,iyy,izz,irhop) + weight
                if (lparticles_lyapunov) then
                   if (ibbf/=0) then
                        f(ixx,iyy,izz,ibxf:ibzf)=f(ixx,iyy,izz,ibxf:ibzf)+weight*fp(k,ibpx:ibpz) 
                   endif
                endif
!
                if (ipeh/=0 .and. (lcylindrical_coords .or. lspherical_coords)) then
                  if (lparticles_number) then
                    weight = fp(k,inpswarm)*(fp(k,iap)**2.0)
                  else
                    weight = np_swarm*(fp(k,iap)**2.0)
                  endif
                  if (nxgrid/=1) weight=weight*weight_x
                  if (nygrid/=1) weight=weight*weight_y
                  if (nzgrid/=1) weight=weight*weight_z
                  f(ixx,iyy,izz,ipeh) = f(ixx,iyy,izz,ipeh)+weight
                endif
!
              enddo; enddo; enddo
            endif
          enddo
!
!  Nearest Grid Point (NGP) method.
!
        else
          if (lparticles_radius.or.lparticles_number.or.lparticles_density) then
            do k=1,npar_loc
              lsink=.false.
              if (lparticles_sink) then
                if (fp(k,iaps)>0.0) lsink=.true.
              endif
              if (lmapsink) lsink=.not.lsink
              if (.not.lsink) then
                ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
!
                if (lparticles_density) then
                  weight0=fp(k,irhopswarm)
                elseif (lparticles_radius.and.lparticles_number) then
                  weight0=four_pi_rhopmat_over_three*fp(k,iap)**3*fp(k,inpswarm)
                elseif (lparticles_radius) then
                  weight0=four_pi_rhopmat_over_three*fp(k,iap)**3*np_swarm
                elseif (lparticles_number) then
                  weight0=mpmat*fp(k,inpswarm)
                endif
!
                f(ix0,iy0,iz0,irhop)=f(ix0,iy0,iz0,irhop) + weight0
                if (lparticles_lyapunov) then
                   if (ibbf/=0) then
                        f(ixx,iyy,izz,ibxf:ibzf)=f(ixx,iyy,izz,ibxf:ibzf)+weight0*fp(k,ibpx:ibpz) 
                   endif
                endif
!
                if (ipeh/=0 .and. (lcylindrical_coords .or. lspherical_coords)) then
                  if (lparticles_number) then
                    weight0 = fp(k,inpswarm)*(fp(k,iap)**2.0)
                  else
                    weight0 = np_swarm*(fp(k,iap)**2.0)
                  endif
                  f(ixx,iyy,izz,ipeh) = f(ixx,iyy,izz,ipeh)+weight0
                endif
!
              endif
            enddo
          else
            f(l1:l2,m1:m2,n1:n2,irhop)=f(l1:l2,m1:m2,n1:n2,inp)
            if (lparticles_lyapunov) then
              if (ibbf/=0) then
                 f(ixx,iyy,izz,ibxf:ibzf)=f(ixx,iyy,izz,ibxf:ibzf)+fp(k,ibpx:ibpz) 
              endif
            endif
          endif
        endif
!
!  Fold first ghost zone of f.
!
        if (lparticlemesh_cic.or.lparticlemesh_tsc) then 
          call fold_f(f,irhop,irhop)
          if (lparticles_lyapunov.and.ibbf/=0) then
            call fold_f(f,ibxf,ibxf)
            call fold_f(f,ibyf,ibyf)
            call fold_f(f,ibzf,ibzf)
          endif
        endif
        if (.not.(lparticles_radius.or.lparticles_number.or. &
            lparticles_density)) then
          do m=m1,m2; do n=n1,n2
            call get_rhopswarm(mp_swarm,fp,1,m,n,rhop_swarm_mn)
            f(l1:l2,m,n,irhop)=rhop_swarm_mn*f(l1:l2,m,n,irhop)
          enddo; enddo
        endif
!        call sharpen_tsc_density(f)
      endif
!
!  Restore normal particle index if mapping sink particles.
!
      if (lmapsink) irhop=irhopm
!
    endsubroutine map_xxp_grid
!***********************************************************************
    subroutine map_xxp_grid__loop_ap(fp,k,index)
!
!  Find the radius bin of particle(k):
!
!  10-mar-20/christer: coded
!
      use SharedVariables, only: get_shared_variable
      integer, intent(in) :: k
      integer, intent(out) :: index
      integer :: ik,ierr
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, pointer :: ap0(:)
!
      call get_shared_variable('ap0',ap0,ierr)
      loop_ap__cic:do ik=1,ndustrad
        if ((fp(k,iap) > ap0(ik)*0.99) .and. (fp(k,iap) < ap0(ik)*1.01)) then
          index = iapn(ik)
          exit loop_ap__cic
        endif
      enddo loop_ap__cic
!
    endsubroutine map_xxp_grid__loop_ap
!***********************************************************************
    subroutine map_vvp_grid(f,fp,ineargrid)
!
!  Map the particle velocities as vector field on the grid.
!
!  07-oct-08/anders: coded
!
      use GhostFold,     only: fold_f
      use Particles_sub, only: get_rhopswarm
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
      real, dimension(nx) :: rhop_swarm_mn
      real :: weight, weight0, weight_x, weight_y, weight_z
      integer :: ivp, k, ix0, iy0, iz0, ixx, iyy, izz
      integer :: ixx0, ixx1, iyy0, iyy1, izz0, izz1
!
!  Calculate the smooth velocity field of particles in each grid cell. Three
!  methods are implemented for assigning a particle to the mesh (see Hockney &
!  Eastwood):
!
!    0. NGP (Nearest Grid Point)
!       The entire effect of the particle goes to the nearest grid point.
!    1. CIC (Cloud In Cell)
!       The particle has a region of influence with the size of a grid cell.
!       This is equivalent to a first order (spline) interpolation scheme.
!    2. TSC (Triangular Shaped Cloud)
!       The particle is spread over a length of two grid cells, but with
!       a density that falls linearly outwards.
!       This is equivalent to a second order spline interpolation scheme.
!
      if (iupx/=0) then
        do ivp=0,2
          f(:,:,:,iupx+ivp)=0.0
          if (lparticlemesh_cic) then
!
!  Cloud In Cell (CIC) scheme.
!
            do k=1,npar_loc
                ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
                ixx0=ix0; iyy0=iy0; izz0=iz0
                ixx1=ix0; iyy1=iy0; izz1=iz0
                if ( (x(ix0)>fp(k,ixp)) .and. nxgrid/=1) ixx0=ixx0-1
                if ( (y(iy0)>fp(k,iyp)) .and. nygrid/=1) iyy0=iyy0-1
                if ( (z(iz0)>fp(k,izp)) .and. nzgrid/=1) izz0=izz0-1
                if (nxgrid/=1) ixx1=ixx0+1
                if (nygrid/=1) iyy1=iyy0+1
                if (nzgrid/=1) izz1=izz0+1
                call get_rhopswarm(mp_swarm, fp, k, ix0, iy0, iz0, weight0)
                do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
                  weight = weight0
                  if (nxgrid/=1) &
                       weight=weight*( 1.0-abs(fp(k,ixp)-x(ixx))*dx_1(ixx) )
                  if (nygrid/=1) &
                       weight=weight*( 1.0-abs(fp(k,iyp)-y(iyy))*dy_1(iyy) )
                  if (nzgrid/=1) &
                       weight=weight*( 1.0-abs(fp(k,izp)-z(izz))*dz_1(izz) )
                  f(ixx,iyy,izz,iupx+ivp)=f(ixx,iyy,izz,iupx+ivp) + &
                      weight*fp(k,ivpx+ivp)
                enddo; enddo; enddo
            enddo
!
!  Triangular Shaped Cloud (TSC) scheme.
!
          elseif (lparticlemesh_tsc) then
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
            do k=1,npar_loc
                ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
                if (nxgrid/=1) then
                  ixx0=ix0-1; ixx1=ix0+1
                else
                  ixx0=ix0  ; ixx1=ix0
                endif
                if (nygrid/=1) then
                  iyy0=iy0-1; iyy1=iy0+1
                else
                  iyy0=iy0  ; iyy1=iy0
                endif
                if (nzgrid/=1) then
                  izz0=iz0-1; izz1=iz0+1
                else
                  izz0=iz0  ; izz1=iz0
                endif
!
                call get_rhopswarm(mp_swarm, fp, k, ix0, iy0, iz0, weight0)
!
!  The nearest grid point is influenced differently than the left and right
!  neighbours are. A particle that is situated exactly on a grid point gives
!  3/4 contribution to that grid point and 1/8 to each of the neighbours.
!
                do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
                  if ( ((ixx-ix0)==-1) .or. ((ixx-ix0)==+1) ) then
                    weight_x = 1.125 - 1.5* abs(fp(k,ixp)-x(ixx))*dx_1(ixx) + &
                                       0.5*(abs(fp(k,ixp)-x(ixx))*dx_1(ixx))**2
                  else
                    if (nxgrid/=1) &
                         weight_x = 0.75  -   ((fp(k,ixp)-x(ixx))*dx_1(ixx))**2
                  endif
                  if ( ((iyy-iy0)==-1) .or. ((iyy-iy0)==+1) ) then
                    weight_y = 1.125 - 1.5* abs(fp(k,iyp)-y(iyy))*dy_1(iyy) + &
                                       0.5*(abs(fp(k,iyp)-y(iyy))*dy_1(iyy))**2
                  else
                    if (nygrid/=1) &
                         weight_y = 0.75  -   ((fp(k,iyp)-y(iyy))*dy_1(iyy))**2
                  endif
                  if ( ((izz-iz0)==-1) .or. ((izz-iz0)==+1) ) then
                    weight_z = 1.125 - 1.5* abs(fp(k,izp)-z(izz))*dz_1(izz) + &
                                       0.5*(abs(fp(k,izp)-z(izz))*dz_1(izz))**2
                  else
                    if (nzgrid/=1) &
                         weight_z = 0.75  -   ((fp(k,izp)-z(izz))*dz_1(izz))**2
                  endif
!
                  weight = weight0
                  if (nxgrid/=1) weight=weight*weight_x
                  if (nygrid/=1) weight=weight*weight_y
                  if (nzgrid/=1) weight=weight*weight_z
                  f(ixx,iyy,izz,iupx+ivp)=f(ixx,iyy,izz,iupx+ivp) + &
                      weight*fp(k,ivpx+ivp)
                enddo; enddo; enddo
            enddo
          else
!
!  Nearest Grid Point (NGP) method.
!
            do k=1,npar_loc
              ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
              call get_rhopswarm(mp_swarm, fp, k, ix0, iy0, iz0, weight)
              f(ix0,iy0,iz0,iupx+ivp) = f(ix0,iy0,iz0,iupx+ivp) + weight * fp(k,ivpx+ivp)
            enddo
          endif
!
!  Fold first ghost zone of f.
!
          if (lparticlemesh_cic.or.lparticlemesh_tsc) &
              call fold_f(f,iupx+ivp,iupx+ivp)
!
!  Normalize the assigned momentum by the particle density in the grid cell.
!
          do m=m1,m2 ; do n=n1,n2
            where (f(l1:l2,m,n,irhop)/=0.0)
              f(l1:l2,m,n,iupx+ivp) = f(l1:l2,m,n,iupx+ivp) / f(l1:l2,m,n,irhop)
            endwhere
          enddo;enddo
        enddo
!
      endif
!
    endsubroutine map_vvp_grid
!***********************************************************************
    subroutine particle_pencil_index(ineargrid)
!
!  Calculate the beginning and ending index of particles in a pencil.
!
!  24-apr-06/anders: coded
!
      integer, dimension (mpar_loc,3) :: ineargrid
!
      integer :: k, iy0, iz0
!
      intent(in)  :: ineargrid
!
      npar_imn=0
!
!  Calculate the number of particles in each pencil.
!
      do k=1,npar_loc
        iy0=ineargrid(k,2); iz0=ineargrid(k,3)
        npar_imn(imn_array(iy0,iz0))=npar_imn(imn_array(iy0,iz0))+1
      enddo
!
!  Calculate beginning and ending particle index for each pencil.
!
      k=0
      do imn=1,ny*nz
        if (npar_imn(imn)/=0) then
          k1_imn(imn)=k + 1
          k2_imn(imn)=k1_imn(imn) + npar_imn(imn) - 1
          k=k+npar_imn(imn)
        else
          k1_imn(imn)=0
          k2_imn(imn)=0
        endif
      enddo
!
    endsubroutine particle_pencil_index
!***********************************************************************
    subroutine shepherd_neighbour_pencil(fp,ineargrid,kshepherd,kneighbour)
!
!  Create a shepherd/neighbour list of particles in the pencil.
!
!  24-oct-05/anders: coded
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (nx) :: kshepherd
      integer, dimension (:) :: kneighbour
!
      integer :: k, ix0
!
      kshepherd=0
      if (imn==1) kneighbour=0
!
      if (npar_imn(imn)/=0) then
        do k=k1_imn(imn),k2_imn(imn)
          ix0=ineargrid(k,1)
          kneighbour(k)=kshepherd(ix0-nghost)
          kshepherd(ix0-nghost)=k
        enddo
      endif
!
      call keep_compiler_quiet(fp)
!
    endsubroutine shepherd_neighbour_pencil
!***********************************************************************
    subroutine shepherd_neighbour_block(fp,ineargrid,kshepherdb,kneighbour, &
        iblock)
!
!  Create a shepherd/neighbour list of particles in the block.
!
!  17-nov-09/anders: dummy
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (nxb,nyb,nzb) :: kshepherdb
      integer, dimension (:) :: kneighbour
      integer :: iblock
!
      intent (in) :: fp, ineargrid, iblock
      intent (out) :: kshepherdb, kneighbour
!
      call fatal_error('shepherd_neighbour_block', &
          'only implemented for block domain decomposition')
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(kshepherdb(1,1,1))
      call keep_compiler_quiet(kneighbour)
      call keep_compiler_quiet(iblock)
!
    endsubroutine shepherd_neighbour_block
!***********************************************************************
    subroutine interpolation_consistency_check()
!
!  Check that all interpolation requirements are satisfied:
!
      use Particles_cdata
!
      if (interp%luu .and. (.not. &
          (lhydro.or.(lhydro_kinematic.and.lkinflow_as_aux)))) then
        call warning('initialize_particles','interpolated uu '// &
          'is set to zero because there is no Hydro module!')
      endif
!
      if (interp%loo .and. (.not.lparticles_spin)) then
        call fatal_error('initialize_particles','interpolation of oo '// &
          'impossible without the Particles_lift module!')
      endif
!
      if (interp%lTT .and. ((.not.lentropy).and.(.not.ltemperature))) then
        call fatal_error('initialize_particles','interpolation of TT '//&
          'impossible without the Entropy (or temperature_idealgas) module!')
      endif
!
!  allow for kinematic flows as well
!
      if (.not.lhydro_kinematic) then
        if (interp%lrho .and. (.not.ldensity) .and. (.not.lanelastic)) then
          call fatal_error('initialize_particles','interpolation of rho '// &
            'impossible without the Density module!')
        endif
      endif
!
    endsubroutine interpolation_consistency_check
!***********************************************************************
    subroutine interpolate_quantities(f,fp,p,ineargrid)
!
!  Interpolate the needed sub-grid quantities according to preselected
!  interpolation policies.
!
!  28-jul-08/kapelrud: coded
!
      use Particles_cdata
!
      real,dimension(mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer,dimension(mpar_loc,3) :: ineargrid
      type (pencil_case) :: p
!
      intent(in) :: f,fp,ineargrid
!
      integer :: k1,k2,k
!
      k1=k1_imn(imn); k2=k2_imn(imn)
!
!  Flow velocity:
!
      if (interp%luu) then
        allocate(interp_uu(k1:k2,3))
        if (.not.allocated(interp_uu)) then
          print*,'interpolate_quantities: unable to allocate '// &
                 'sufficient memory for interp_uu'
          call fatal_error('interpolate_quantities','')
        endif
!
!  For kinematic flow, we must have the flow in the f-array
!
        if (lhydro.or.(lhydro_kinematic.and.lkinflow_as_aux)) then
          call interp_field_pencil_wrap(f,iux,iuz,fp,ineargrid,interp_uu, &
              interp%pol_uu)
        else
          interp_uu=0.0
        endif
      endif
!
!  Flow vorticity:
!
      if (interp%loo) then
        allocate(interp_oo(k1:k2,3))
        if (.not.allocated(interp_oo)) then
          print*,'interpolate_quantities: unable to allocate '// &
                 'sufficient memory for interp_oo'
          call fatal_error('interpolate_quantities','')
        endif
        call interp_field_pencil_wrap(f,iox,ioz,fp,ineargrid,interp_oo, &
          interp%pol_oo)
      endif
!
!  Temperature:
!
      if (interp%lTT) then
        allocate(interp_TT(k1:k2))
        if (.not.allocated(interp_TT)) then
          print*,'interpolate_quantities: unable to allocate '// &
                 'sufficient memory for interp_TT'
          call fatal_error('interpolate_quantities','')
        endif
!
!  Determine what quantity to interpolate (temperature or entropy)
!
        if (lentropy) then
          call interp_field_pencil_wrap(f,iss,iss,fp,ineargrid,interp_TT, &
            interp%pol_TT)
        elseif (ltemperature) then
          call interp_field_pencil_wrap(f,ilnTT,ilnTT,fp,ineargrid,interp_TT, &
               interp%pol_TT)
        endif
!
        if ( (lentropy.and.pretend_lnTT)   .or. &
           (ltemperature.and.(.not.ltemperature_nolog)) ) then
          interp_TT=exp(interp_TT)
        elseif (lentropy.and.(.not.pretend_lnTT)) then
          call fatal_error('interpolate_quantities','enable flag '//&
            'pretend_lnTT in init_pars to be able to interpolate'// &
            ' the temperature if using the regular temperature module!')
        endif
      endif
!
!  Temperature gradient:
!
      if (interp%lgradTT) then
        allocate(interp_gradTT(k1:k2,3))
        if (.not.allocated(interp_gradTT)) then
          print*,'interpolate_quantities: unable to allocate '// &
                 'sufficient memory for interp_gradTT'
          call fatal_error('interpolate_quantities','')
        endif
!
        if (npar_imn(imn) /= 0) then
           do k=k1,k2
              interp_gradTT(k,:)=p%gTT(ineargrid(k,1)-nghost,:)
           enddo
        endif
      endif
!
!  Density:
!
      if (interp%lrho) then
        allocate(interp_rho(k1:k2))
        if (.not.allocated(interp_rho)) then
          print*,'interpolate_quantities: unable to allocate '// &
                 'sufficient memory for interp_rho'
          call fatal_error('interpolate_quantities','')
        endif
        if (ldensity) then
          call interp_field_pencil_wrap(f,ilnrho,ilnrho,fp,ineargrid,&
            interp_rho,interp%pol_rho)
        else
          interp_rho=1. !(should be made general)
        endif
        if (.not.ldensity_nolog) then
          interp_rho=exp(interp_rho)
        endif
      endif
!
!  Pressure:
!
      if (interp%lpp) then
!!$        call fatal_error('interpolate_quantities',&
!!$            'Check that interpolation of pressure is properly implemented!')
        allocate(interp_pp(k1:k2))
        if (.not.allocated(interp_pp)) then
          print*,'interpolate_quantities: unable to allocate '// &
                 'sufficient memory for interp_pp'
          call fatal_error('interpolate_quantities','')
        endif
        if (npar_imn(imn) /= 0) then
          do k=k1,k2
            interp_pp(k)=p%pp(ineargrid(k,1)-nghost)
          enddo
        endif
      endif
!
!  Viscosity:
!
      if (interp%lnu) then
        allocate(interp_nu(k1:k2))
        if (.not.allocated(interp_nu)) then
          print*,'interpolate_quantities: unable to allocate '// &
                 'sufficient memory for interp_nu'
          call fatal_error('interpolate_quantities','')
        endif
        if (npar_imn(imn) /= 0) then
          do k=k1,k2
            interp_nu(k)=p%nu(ineargrid(k,1)-nghost)
          enddo
        endif
      endif
!
!  Species:
!
      if (interp%lspecies) then
        allocate(interp_species(k1:k2,nchemspec))
        if (.not.allocated(interp_species)) then
          print*,'interpolate_quantities: unable to allocate '// &
                 'sufficient memory for interp_species'
          call fatal_error('interpolate_quantities','')
        endif
        call interp_field_pencil_wrap(f,ichemspec(1),ichemspec(nchemspec),&
            fp,ineargrid,interp_species,interp%pol_species)
      endif
!
    endsubroutine interpolate_quantities
!***********************************************************************
    subroutine cleanup_interpolated_quantities
!
!  Deallocate memory from particle pencil interpolation variables
!
!  28-jul-08/kapelrud: coded
!
      use Particles_cdata
!
      if (allocated(interp_uu)) deallocate(interp_uu)
      if (allocated(interp_oo)) deallocate(interp_oo)
      if (allocated(interp_TT)) deallocate(interp_TT)
      if (allocated(interp_rho)) deallocate(interp_rho)
      if (allocated(interp_gradTT)) deallocate(interp_gradTT)
      if (allocated(interp_nu)) deallocate(interp_nu)
      if (allocated(interp_pp)) deallocate(interp_pp)
      if (allocated(interp_species)) deallocate(interp_species)
!
    endsubroutine cleanup_interpolated_quantities
!***********************************************************************
    subroutine interp_field_pencil_0(f,i1,i2,fp,ineargrid,vec,policy)
!
!  Overloaded interpolation wrapper for scalar fields.
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i1,i2
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
      real, dimension(:) :: vec
      integer :: policy
!
      intent(in) :: f,i1,i2,fp,ineargrid, policy
      intent(inout) :: vec
!
      if(lsolid_ogrid) then
        call interp_field_pencil_ogrid(f,i1,i2,fp,ineargrid,1,vec)
      else
        call interp_field_pencil(f,i1,i2,fp,ineargrid,1,vec,policy)
      endif
!
    endsubroutine interp_field_pencil_0
!***********************************************************************
    subroutine interp_field_pencil_1(f,i1,i2,fp,ineargrid,vec,policy)
!
!  Overloaded interpolation wrapper for vector fields.
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i1,i2
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
      real, dimension(:,:) :: vec
      integer :: policy
!
      intent(in) :: f,i1,i2,fp,ineargrid, policy
      intent(inout) :: vec
!
      if(lsolid_ogrid) then
        call interp_field_pencil_ogrid(f,i1,i2,fp,ineargrid,ubound(vec,2),vec)
      else
         call interp_field_pencil(f,i1,i2,fp,ineargrid,ubound(vec,2),vec,policy)
      endif
!
    endsubroutine interp_field_pencil_1
!***********************************************************************
    subroutine interp_field_pencil(f,i1,i2,fp,ineargrid,uvec2,vec,policy)
!
!  Interpolate stream field to all sub grid particle positions in the
!  current pencil.
!  i1 & i2 sets the component to interpolate; ex.: iux:iuz, or iss:iss
!
!  16-jul-08/kapelrud: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i1,i2
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
      integer :: uvec2,policy
      real, dimension(k1_imn(imn):k2_imn(imn),uvec2) :: vec
!
      intent(in) :: f,i1,i2,fp,ineargrid, policy
      intent(inout) :: vec
!
      integer :: k
!
      if (npar_imn(imn)/=0) then
        do k=k1_imn(imn),k2_imn(imn)
          select case (policy)
          case (cic)
            call interpolate_linear( &
                f,i1,i2,fp(k,ixp:izp),vec(k,:),ineargrid(k,:),0,ipar(k) )
          case (tsc)
            if (linterpolate_spline) then
              call interpolate_quadratic_spline( &
                   f,i1,i2,fp(k,ixp:izp),vec(k,:),ineargrid(k,:),0,ipar(k) )
            else
              call interpolate_quadratic( &
                   f,i1,i2,fp(k,ixp:izp),vec(k,:),ineargrid(k,:),0,ipar(k) )
            endif
          case (ngp)
            vec(k,:)=f(ineargrid(k,1),ineargrid(k,2),ineargrid(k,3),i1:i2)
          endselect
        enddo
      endif
!
    endsubroutine interp_field_pencil
!***********************************************************************
    subroutine interp_field_pencil_ogrid(f,i1,i2,fp,ineargrid,uvec2,vec)
!
!  Interpolate stream field to all sub grid particle positions in the
!  current pencil.
!  i1 & i2 sets the component to interpolate; ex.: iux:iuz, or iss:iss
!
!  16-jul-08/kapelrud: coded
!
      use Solid_Cells, only: interpolate_particles_ogrid, map_nearest_grid_ogrid, &
                             r_int_outer, xorigo_ogrid
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i1,i2
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
      integer :: uvec2
      real, dimension(k1_imn(imn):k2_imn(imn),uvec2) :: vec
!
      real, dimension(k1_imn(imn):k2_imn(imn)) :: rp2
      real :: r2_ogrid
      integer, dimension(4) :: inear_ogrid
      real, dimension(3) :: rthz

      intent(in) :: f,i1,i2,fp,ineargrid
      intent(inout) :: vec
!
      integer :: k
      ! TODO
      real, dimension(k1_imn(imn):k2_imn(imn),uvec2) :: vec2

      r2_ogrid=r_int_outer*r_int_outer
!
      if (npar_imn(imn)/=0) then
        rp2 = (fp(k1_imn(imn):k2_imn(imn),ixp)-xorigo_ogrid(1))**2 &
             +(fp(k1_imn(imn):k2_imn(imn),iyp)-xorigo_ogrid(2))**2
        do k=k1_imn(imn),k2_imn(imn)
          if(rp2(k)>r2_ogrid) then
!
!  NILS: Commented out the line below  since we should enter here even
!  NILS: when iTT=0 (e.g. for interpolation of velocity).
!              if(iTT > 0 .or. i1<=iTT) then
                 call interpolate_linear( &
                      f,i1,i2,fp(k,ixp:izp),vec(k,:),ineargrid(k,:),0,ipar(k) )
 !             endif
           else
              call map_nearest_grid_ogrid(fp(k,ixp:izp),inear_ogrid,rthz)
              call interpolate_particles_ogrid(i1,i2,rthz,vec(k,:),inear_ogrid)
           endif
         
        enddo
     endif
!
    endsubroutine interp_field_pencil_ogrid
!***********************************************************************
    subroutine sort_particles_iblock(fp,ineargrid,ipar,dfp)
!
!  Sort the particles so that they appear in order of the global brick index.
!  That is, sorted first by processor number and then by local brick index.
!
!  16-nov-09/anders: dummy
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
!
      call fatal_error('sort_particles_iblock', &
          'only implemented for block domain decomposition')
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(ipar)
      if (present(dfp)) call keep_compiler_quiet(dfp)
!
    endsubroutine sort_particles_iblock
!***********************************************************************
    subroutine fill_blocks_with_bricks(a,ab,marray,ivar)
!
!  Fill adopted blocks with bricks from the f-array.
!
!  04-nov-09/anders: dummy
!
      integer :: marray
      real, dimension (mx,my,mz,marray) :: a
      real, dimension (mxb,myb,mzb,marray,0:nblockmax-1) :: ab
      integer :: ivar
!
      call keep_compiler_quiet(a)
      call keep_compiler_quiet(ab(1,1,1,1,0))
      call keep_compiler_quiet(marray)
      call keep_compiler_quiet(ivar)
!
    endsubroutine fill_blocks_with_bricks
!***********************************************************************
    subroutine fill_bricks_with_blocks(a,ab,marray,ivar)
!
!  Fill adopted blocks with bricks from the f-array.
!
!  04-nov-09/anders: dummy
!
      integer :: marray
      real, dimension (mx,my,mz,marray) :: a
      real, dimension (mxb,myb,mzb,marray,0:nblockmax-1) :: ab
      integer :: ivar
!
      call keep_compiler_quiet(a)
      call keep_compiler_quiet(ab(1,1,1,1,0))
      call keep_compiler_quiet(marray)
      call keep_compiler_quiet(ivar)
!
    endsubroutine fill_bricks_with_blocks
!***********************************************************************
    subroutine shepherd_neighbour_pencil3d(fp,ineargrid_c,kshepherd_c,kneighbour_c)
!
!  20-July-2010: coded AlexHubbard
!
!  Create a shepherd/neighbour list of particles in the pencil.
!  On collisional grid
!  Adapted from particles_map
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid_c
      integer, dimension (:,:,:) :: kshepherd_c
      integer, dimension (mpar_loc) :: kneighbour_c
!
      integer :: k, ix0, iy0, iz0
!
      kshepherd_c=0
      kneighbour_c=0
!
! kshepherd contains highest k particle for each collisional grid point
! kneighbour contains next highest k particle on the same collisional grid
! point as that of particle k.
!
      do k=1,npar_loc
        ix0=ineargrid_c(k,1); iy0=ineargrid_c(k,2);iz0=ineargrid_c(k,3)
        kneighbour_c(k)=kshepherd_c(ix0, iy0, iz0)
        kshepherd_c(ix0, iy0, iz0)=k
      enddo
!
      call keep_compiler_quiet(fp)
!
    endsubroutine shepherd_neighbour_pencil3d
!***********************************************************************
endmodule Particles_map
