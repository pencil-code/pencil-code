module Grid

  implicit none

  contains
!***********************************************************************
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
      real, intent(out) :: dx,dy,dz

      real :: x00,y00,z00

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
!  set x,y,z arrays
!
        do i=1,mx; x(i)=x00+(i-nghost-1+ipx*nx)*dx; enddo
        do i=1,my; y(i)=y00+(i-nghost-1+ipy*ny)*dy; enddo
        do i=1,mz; z(i)=z00+(i-nghost-1+ipz*nz)*dz; enddo


    endsubroutine construct_grid
!***********************************************************************
endmodule Grid
