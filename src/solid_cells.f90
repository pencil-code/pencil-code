! $Id$
!
!  This module add solid (as in no-fluid) cells in the domain.
!  This can be used e.g. in order to simulate a cylinder in a cross flow.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lsolid_cells = .true.
!
!***************************************************************
module Solid_Cells
!
  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!  
  include 'solid_cells.h'
!
  integer, parameter            :: max_items=5
  integer                       :: ncylinders=0,nrectangles,nspheres=0,dummy
  integer                       :: nobjects, nlong, nlat
  integer                       :: nforcepoints=200
  real, dimension(max_items)    :: cylinder_radius
  real, dimension(max_items)    :: cylinder_temp=703.0
  real, dimension(max_items)    :: cylinder_xpos,cylinder_ypos,cylinder_zpos
  real, dimension(max_items)    :: sphere_radius
  real, dimension(max_items)    :: sphere_temp=703.0
  real, dimension(max_items)    :: sphere_xpos,sphere_ypos,sphere_zpos
  integer, dimension(mx,my,mz,4):: ba,ba_shift
  real :: skin_depth=0, init_uu=0, ampl_noise=0, object_skin=0
  character (len=labellen), dimension(ninit) :: initsolid_cells='nothing'
  character (len=labellen) :: interpolation_method='staircase'
  logical :: lclose_interpolation=.false., lclose_linear=.false.
  logical :: lnointerception=.false.
  real                          :: rhosum
  integer                       :: irhocount
  real                          :: theta_shift=1e-2
  real                          :: limit_close_linear=0.5
!
  type solid_object
    character(len=10) :: form
    real :: r,T
    real, dimension(3) :: x
  end type solid_object
!
  type(solid_object), dimension(max_items) :: objects
!
  namelist /solid_cells_init_pars/ &
       cylinder_temp, ncylinders, cylinder_radius, cylinder_xpos, &
       cylinder_ypos, cylinder_zpos, initsolid_cells, skin_depth, init_uu, &
       ampl_noise,interpolation_method, nforcepoints,object_skin,&
       lclose_interpolation,lclose_linear,limit_close_linear,lnointerception,&
       nspheres,sphere_radius,sphere_xpos,sphere_ypos,sphere_zpos,sphere_temp
!
  namelist /solid_cells_run_pars/  &
       interpolation_method,object_skin,lclose_interpolation,lclose_linear,&
       limit_close_linear,lnointerception
!
!  Diagnostic variables (need to be consistent with reset list below).
!
  integer :: idiag_c_dragx=0       ! DIAG_DOC: 
  integer :: idiag_c_dragy=0       ! DIAG_DOC: 
  integer :: idiag_c_dragz=0       ! DIAG_DOC: 
!
  integer, allocatable :: fpnearestgrid(:,:,:)
  real, allocatable    :: c_dragx(:), c_dragy(:), c_dragz(:)
!
  contains
!***********************************************************************
    subroutine initialize_solid_cells
!
!  Define the geometry of the solids.
!  There might be many separate solid objects of different geometries (currently
!  only cylinders and spheres are implemented however).
!
!  19-nov-2008/nils: coded
!  28-sep-2010/nils: added spheres
!  nov-2010/kragset: updated allocations related to drag calculations
!
      integer :: icyl,isph      
!
!  Define the geometry of the solid object.
!  For more complex geometries (i.e. for objects different than cylinders,
!  spheres or rectangles) this shold probably be included as a geometry.local 
!  file such that
!  one can define complex geometries on a case to case basis.
!  Alternatively one will here end up with a terribly long series
!  of case checks. 
!
!  Loop over all cylinders
!
      do icyl=1,ncylinders
        if (cylinder_radius(icyl)>0) then
          objects(icyl)%r    = cylinder_radius(icyl)
          objects(icyl)%x(1) = cylinder_xpos(icyl)
          objects(icyl)%x(2) = cylinder_ypos(icyl)
          objects(icyl)%x(3) = cylinder_zpos(icyl)
          objects(icyl)%T    = cylinder_temp(icyl)
          objects(icyl)%form = 'cylinder'
        else
          call fatal_error('initialize_solid_cells',&
               'All cylinders must have non-zero radii!')
        endif
      enddo
!
!  Loop over all spheres
!
      do isph=1,nspheres
        if (sphere_radius(isph)>0) then
          objects(isph+ncylinders)%r   = sphere_radius(isph)
          objects(isph+ncylinders)%x(1)= sphere_xpos(isph)
          objects(isph+ncylinders)%x(2)= sphere_ypos(isph)
          objects(isph+ncylinders)%x(3)= sphere_zpos(isph)
          objects(isph+ncylinders)%T   = sphere_temp(isph)
          objects(isph+ncylinders)%form= 'sphere'
        else
          call fatal_error('initialize_solid_cells',&
               'All spheres must have non-zero radii!')
        endif
      enddo
      nobjects=ncylinders+nspheres      
!
!  In order to avoid problems with how namelists are written not all
!  slots of an array which should be stored in the namelist can be zero
!
      if (nspheres==0) then
        sphere_radius(1)=impossible
        sphere_xpos(1)=impossible
        sphere_ypos(1)=impossible
        sphere_zpos(1)=impossible
        sphere_temp(1)=impossible
      else
!
! Find number of lines of longitude and latitude such that nforcepoints=nlong*(nlat+1)
! and nlat=nlong/2-1
!
        nlong=int(sqrt(2.*nforcepoints))
        nlat=int(.5*nlong)-1
        if (nlong*(nlat+1)/=nforcepoints) print*, "Warning: 2*nforcepoints should be square"
      end if
      if (ncylinders==0) then
        cylinder_radius(1)=impossible
        cylinder_xpos(1)=impossible
        cylinder_ypos(1)=impossible
        cylinder_zpos(1)=impossible
        cylinder_temp(1)=impossible
      end if
!
!  Prepare the solid geometry
!
      call find_solid_cell_boundaries
      call calculate_shift_matrix
!
!
! Find nearest grid point of the "forcepoints" on all cylinders. Arrays
! are only allocated if c_dragx, c_dragy or c_dragz is set in print.in. 
! This must be changed if drag force is required for other purposes, 
! e.g. if solid object is allowed to move.
!
      if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
          idiag_c_dragz /= 0) then 
        allocate(c_dragx(nobjects))
        allocate(c_dragy(nobjects))
        allocate(c_dragz(nobjects))
        allocate(fpnearestgrid(nobjects,nforcepoints,3))
        call fp_nearest_grid
        rhosum    = 0.0
        irhocount = 0
      end if
!
    endsubroutine initialize_solid_cells
!***********************************************************************
    subroutine init_solid_cells(f)
!
!  Initial conditions for cases where we have solid structures in the domain.
!  Typically the flow field is set such that we have no-slip conditions
!  at the solid structure surface.
! 
!  28-nov-2008/nils: coded
!
      use Cdata
      use Sub
      use Initcond
      use InitialCondition, only: initial_condition_solid_cells
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: a2,rr2, wall_smoothing,rr2_low,rr2_high,shiftx,shifty
      real :: wall_smoothing_temp,xr,yr
      integer i,j,k,cyl,jj,icyl
!
      do jj=1,ninit
      select case (initsolid_cells(jj))
!
!  This overrides any initial conditions set in the Hydro module.
!
      case ('nothing')
        if (lroot) print*,'init_solid_cells: nothing'
      case ('cylinderstream_x')
!  Stream functions for flow around a cylinder as initial condition. 
        call gaunoise(ampl_noise,f,iux,iuz)
        f(:,:,:,iux)=f(:,:,:,iux)+init_uu
        shiftx=0
        do i=l1,l2
        do j=m1,m2
        do k=n1,n2
!
!  Loop over all cylinders
!
          do icyl=1,ncylinders
            a2 = objects(icyl)%r**2
            xr=x(i)-objects(icyl)%x(1)
            if (objects(icyl)%x(2) /= 0) then
              print*,'When using cylinderstream_x all cylinders must have'
              print*,'zero offset in y-direction!'
              call fatal_error('init_solid_cells:','')
            endif
            yr=y(j)
            rr2 = xr**2+yr**2
            if (rr2 > a2) then
              do cyl=0,100
                if (cyl==0) then
                  wall_smoothing=1-exp(-(rr2-a2)/skin_depth**2)
                  f(i,j,k,iuy) = f(i,j,k,iuy)-init_uu*&
                       2*xr*yr*a2/rr2**2*wall_smoothing
                  f(i,j,k,iux) = f(i,j,k,iux)+init_uu*&
                       (0. - a2/rr2 + 2*yr**2*a2/rr2**2)&
                       *wall_smoothing
                  if (ilnTT /= 0) then
                    wall_smoothing_temp=1-exp(-(rr2-a2)/(sqrt(a2))**2)
                    f(i,j,k,ilnTT) = wall_smoothing_temp*f(i,j,k,ilnTT)&
                         +objects(icyl)%T*(1-wall_smoothing_temp)
                    f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)&
                         *f(l2,m2,n2,ilnTT)/f(i,j,k,ilnTT)
                  endif
                else
                  shifty=cyl*Lxyz(2)
                  rr2_low =(xr+shiftx)**2+(yr+shifty)**2
                  rr2_high=(xr-shiftx)**2+(yr-shifty)**2
                  f(i,j,k,iux) = f(i,j,k,iux)+init_uu*( &
                       +2*(yr-shifty)**2*a2/rr2_high**2-a2/rr2_high&
                       +2*(yr+shifty)**2*a2/rr2_low**2 -a2/rr2_low)
                  f(i,j,k,iuy) = f(i,j,k,iuy)-init_uu*( &
                       +2*(xr-shiftx)*(yr-shifty)&
                       *a2/rr2_high**2&
                       +2*(xr+shiftx)*(yr+shifty)&
                       *a2/rr2_low**2)
                endif
              enddo
            else
              if (ilnTT /= 0) then
                f(i,j,k,ilnTT) = objects(icyl)%T
                f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)&
                     *f(l2,m2,n2,ilnTT)/objects(icyl)%T
              endif
            endif
          enddo
        enddo
        enddo
        enddo
      case ('cylinderstream_y')
!  Stream functions for flow around a cylinder as initial condition.
        call gaunoise(ampl_noise,f,iux,iuz)
        f(:,:,:,iuy)=f(:,:,:,iuy)+init_uu
        shifty=0
        do i=l1,l2
        do j=m1,m2
        do k=n1,n2
          do icyl=1,ncylinders
            a2 = objects(icyl)%r**2
            yr=y(j)-objects(icyl)%x(2)
            if (objects(icyl)%x(1) /= 0) then
              print*,'When using cylinderstream_y all cylinders must have'
              print*,'zero offset in x-direction!'
              call fatal_error('init_solid_cells:','')
            endif
            xr=x(i)
            rr2 = xr**2+yr**2
            if (rr2 > a2) then
              do cyl=0,100
                if (cyl==0) then
                  wall_smoothing=1-exp(-(rr2-a2)/skin_depth**2)
                  f(i,j,k,iux) = f(i,j,k,iux)-init_uu*&
                       2*xr*yr*a2/rr2**2*wall_smoothing
                  f(i,j,k,iuy) = f(i,j,k,iuy)+init_uu*&
                       (0. - a2/rr2 + 2*xr**2*a2/rr2**2)&
                       *wall_smoothing
                  if (ilnTT /= 0) then
                    wall_smoothing_temp=1-exp(-(rr2-a2)/(sqrt(a2))**2)
                    f(i,j,k,ilnTT) = wall_smoothing_temp*f(i,j,k,ilnTT)&
                         +objects(icyl)%T*(1-wall_smoothing_temp)
                    f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)&
                         *f(l2,m2,n2,ilnTT)/f(i,j,k,ilnTT)
                  endif
                else
                  shiftx=cyl*Lxyz(1)
                  rr2_low =(xr+shiftx)**2+(yr+shifty)**2
                  rr2_high=(xr-shiftx)**2+(yr-shifty)**2
                  f(i,j,k,iuy) = f(i,j,k,iuy)+init_uu*( &
                       +2*(xr-shiftx)**2*a2/rr2_high**2-a2/rr2_high&
                       +2*(xr+shiftx)**2*a2/rr2_low**2 -a2/rr2_low)
                  f(i,j,k,iux) = f(i,j,k,iux)-init_uu*( &
                       +2*(xr-shiftx)*(y(j)-shifty)&
                       *a2/rr2_high**2&
                       +2*(xr+shiftx)*(y(j)+shifty)&
                       *a2/rr2_low**2)
                endif
              enddo
            else
              if (ilnTT /= 0) then
                f(i,j,k,ilnTT) = objects(icyl)%T
                f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)&
                     *f(l2,m2,n2,ilnTT)/objects(icyl)%T
              endif
            endif
          enddo
        enddo
        enddo
        enddo
if (llast_proc_y) f(:,m2-5:m2,:,iux)=0
      case default
!
!  Catch unknown values
!
        if (lroot) print*,'No such value for init_solid_cells:',&
             trim(initsolid_cells(jj))
        call fatal_error('init_solid_cells','')
      endselect
    enddo
!
!  Interface for user's own initial condition
!
    if (linitial_condition) call initial_condition_solid_cells(f)
!
    endsubroutine init_solid_cells
!***********************************************************************  
  subroutine fp_nearest_grid
!
!  Find coordinates for nearest grid point of all the 
!  "forcepoints" (fp) for each object (assume object with axis
!  parallel to the z direction. Assign values to fpnearestgrid.
!
!  mar-2009/kragset: coded
!  nov-2010/kragset: updated to include spheres
!
!
    integer           :: iobj,iforcepoint, ipoint, inearest, icoord(8,3)
    integer           :: ilong,ilat
    integer           :: ixl, iyl, izl, ixu, iyu, izu, ju, jl, jm
    real              :: robj, xobj, yobj, zobj,fpx, fpy, fpz
    real              :: dx1, dy1, dz1, longitude, latitude
    real              :: dist_to_fp2(8), dist_to_cent2(8), twopi,dlong,dlat
    logical           :: interiorpoint
    character(len=10) :: objectform
!
    dx1=1/dx
    dy1=1/dy
    dz1=1/dz
!
    twopi=2.*pi
!
!  Loop over all objects 
    do iobj=1,nobjects
      robj = objects(iobj)%r
      xobj = objects(iobj)%x(1)
      yobj = objects(iobj)%x(2)
      objectform = objects(iobj)%form
      if (objectform == 'cylinder') then 
        zobj = z(n1)
        dlong = twopi/nforcepoints
      else if (objectform == 'sphere') then
        zobj  = objects(iobj)%x(3)
        dlong = twopi/nlong
        dlat  = pi/(nlat+1)
      else
        print*, "Warning: Subroutine fp_nearest_grid not implemented ", &
            "for this objectform."
      end if

!
!  Loop over all forcepoints on each object, iobj
      do iforcepoint=1,nforcepoints
!        
!  Marking whether fp is within this processor's domain or not 
        interiorpoint = .true.
!
!  Fp coordinates
!  Shifting the location of the forcpoints in the theta direction
!  in order to avoid problems with autotesting
        if (objectform == 'cylinder') then
          longitude = (iforcepoint-theta_shift)*dlong
          fpx = xobj - robj * sin(longitude)
          fpy = yobj - robj * cos(longitude)
          fpz = z(n1)
        elseif (objectform == 'sphere') then
!  Note definition of lines of longitude: ilong = [0,..,nlong-1]
          ilong  = mod(iforcepoint-1,nlong)
!  Note definition of lines of latitude: ilat  = [0,..,nlat]
          ilat   = int((iforcepoint-1)/nlong)
          longitude = (ilong+.5-theta_shift)*dlong
          latitude  = (ilat+.5)*dlat
          fpx = xobj - robj*sin(longitude)*sin(latitude)
          fpy = yobj - robj*cos(longitude)*sin(latitude)
          fpz = zobj + robj*cos(latitude)
        end if
!
!  Find nearest grid point in x-direction
!          
        if (nxgrid/=1) then
          if (fpx >= x(l1-1) .and. fpx <= x(l2+1)) then
            if (lequidist(1)) then
              ixl = int((fpx-x(1))*dx1) + 1
              ixu = ixl+1
            else
!
!  Find nearest grid point by bisection if grid is not equidistant
!
              ju=l2+1; jl=l1-1
              do while((ju-jl)>1)
                jm=(ju+jl)/2
                if (fpx > x(jm)) then
                  jl=jm
                else
                  ju=jm
                endif
              enddo
              ixl=jl
              ixu=ju
            endif
          else
            interiorpoint=.false.
          endif
        else
          print*,"WARNING: Solid cells need nxgrid > 1."
        endif
!
!  Find nearest grid point in y-direction
!          
        if (nygrid/=1) then
          if (fpy >= y(m1-1) .and. fpy <= y(m2+1)) then
            if (lequidist(2)) then
              iyl = int((fpy-y(1))*dy1) + 1
              iyu = iyl+1
            else
!
!  Find nearest grid point by bisection if grid is not equidistant
!
              ju=m2; jl=m1
              do while((ju-jl)>1)
                jm=(ju+jl)/2
                if (fpy > y(jm)) then
                  jl=jm
                else
                  ju=jm
                endif
              enddo
              iyl=jl
              iyu=ju
            endif
          else
            interiorpoint=.false.
          endif
        else
          print*,"WARNING: Solid cells need nygrid > 1."
        endif
!
!  Find nearest grid point in z-direction
!          
        if (nzgrid/=1) then
          if (fpz >= z(n1-1) .and. fpz <= z(n2+1)) then
            if (lequidist(3)) then
              izl = int((fpz-z(1))*dz1) + 1
              izu = izl+1
            else
!
!  Find nearest grid point by bisection if grid is not equidistant
!
              ju=n2; jl=n1
              do while((ju-jl)>1)
                jm=(ju+jl)/2
                if (fpz > z(jm)) then
                  jl=jm
                else
                  ju=jm
                endif
              enddo
              izl=jl
              izu=ju
            endif
          else
            interiorpoint=.false.
          endif
        else
!  z direction is irrelevant when in 2D
          izl=n1
          izu=n1
        endif
!
!  Now, we have the upper and lower (x,y,z)-coordinates: 
!  ixl, ixu, iyl, iyu, izl, izu,
!  i.e. the eight corners of the grid cell containing the forcepoint (fp).
!  Decide which ones are outside the object, and which one of these
!  is the closest one to fp:
!
!  Check if fp is within this processor's local domain
        if (interiorpoint) then
          dist_to_fp2(1) = (x(ixl)-fpx)**2+(y(iyl)-fpy)**2+(z(izl)-fpz)**2 
          dist_to_fp2(2) = (x(ixu)-fpx)**2+(y(iyl)-fpy)**2+(z(izl)-fpz)**2 
          dist_to_fp2(3) = (x(ixu)-fpx)**2+(y(iyu)-fpy)**2+(z(izl)-fpz)**2 
          dist_to_fp2(4) = (x(ixl)-fpx)**2+(y(iyu)-fpy)**2+(z(izl)-fpz)**2 
          dist_to_fp2(5) = (x(ixl)-fpx)**2+(y(iyl)-fpy)**2+(z(izu)-fpz)**2 
          dist_to_fp2(6) = (x(ixu)-fpx)**2+(y(iyl)-fpy)**2+(z(izu)-fpz)**2 
          dist_to_fp2(7) = (x(ixu)-fpx)**2+(y(iyu)-fpy)**2+(z(izu)-fpz)**2 
          dist_to_fp2(8) = (x(ixl)-fpx)**2+(y(iyu)-fpy)**2+(z(izu)-fpz)**2 
          icoord(1,:) = (/ixl,iyl,izl/)
          icoord(2,:) = (/ixu,iyl,izl/)
          icoord(3,:) = (/ixu,iyu,izl/)
          icoord(4,:) = (/ixl,iyu,izl/)
          icoord(5,:) = (/ixl,iyl,izu/)
          icoord(6,:) = (/ixu,iyl,izu/)
          icoord(7,:) = (/ixu,iyu,izu/)
          icoord(8,:) = (/ixl,iyu,izu/)
          inearest=0
          do ipoint=1,8 
!  Test if we are in a fluid cell, i.e.
!  that mod(ba(ix,iy,iz,1),10) = 0 
            if (mod(ba(icoord(ipoint,1),icoord(ipoint,2),icoord(ipoint,3),1),10)&
                == 0 .and. inearest == 0) then
              inearest=ipoint
            else if ( &
                mod(ba(icoord(ipoint,1),icoord(ipoint,2),icoord(ipoint,3),1),10)&
                == 0 ) then
              if (dist_to_fp2(ipoint) <= dist_to_fp2(inearest)) then
                inearest=ipoint
              endif
            endif
          enddo
!          
!  Coordinates of nearest grid point. Zero if outside local domain.
          if (inearest > 0) then
            fpnearestgrid(iobj,iforcepoint,:) = icoord(inearest,:)
          else
            print*, "WARNING: Could not find fpnearestgrid!"
          endif

        else ! fp is outside local domain and fpnearestgrid shouldn't exist
          fpnearestgrid(iobj,iforcepoint,:) = 0
        endif
      enddo
    enddo
!
  endsubroutine fp_nearest_grid
!***********************************************************************  
  subroutine dsolid_dt(f,df,p)
!
!  Find pressure and stress in all the forcepoints (fp) positioned on 
!  object surface, based on values in nearest grid point.
!
!  mar-2009/kragset: coded
!  okt-2009/kragset: updated to include multiple objects
!  nov-2010/kragset: updated to include spheres
!
    use viscosity, only: getnu
!    
    real, dimension (mx,my,mz,mfarray), intent(in):: f
    real, dimension (mx,my,mz,mvar), intent(in)   :: df
    type (pencil_case), intent(in)                :: p

    real    :: fp_pressure
    real    :: fp_stress(3,3)
    integer :: iobj, ifp, ix0, iy0, iz0, i, ilong, ilat
    real    :: nu, twonu, longitude, latitude, dlong, dlat, robj
    real    :: force_x, force_y, force_z
    real    :: twopi, nvec(3), surfaceelement,surfacecoeff
    character(len=10) :: objectform
!
    if (ldiagnos) then
!      
!  Reset cumulating quantities before calculations in first pencil
!
      if (imn == 1) then
        c_dragx=0.
        c_dragy=0.
        c_dragz=0.
        rhosum=0
        irhocount=0
      endif
!
      if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 &
          .or. idiag_c_dragz /= 0) then 
        call getnu(nu)
        twopi=2.*pi
        twonu=2.*nu
!
        do iobj=1,nobjects
          robj = objects(iobj)%r
          objectform = objects(iobj)%form
          if (objectform=='cylinder') then
            dlong = twopi/nforcepoints
            surfaceelement = dlong
          else if (objectform=='sphere') then
            dlong = twopi/nlong
            dlat  = pi/(nlat+1)
            surfacecoeff = 2.*dlong*dlat/pi
          else
            print*, "Warning: Subroutine dsolid_dt not implemented ", &
                "for this objectform."
          end if
          do ifp=1,nforcepoints
            iy0=fpnearestgrid(iobj,ifp,2)
            iz0=fpnearestgrid(iobj,ifp,3)
            if (objectform=='cylinder') iz0=n 
!
!  Test: Use this pencil for force calculation?
!
            if (iy0 == m .and. iz0 == n) then
              ix0=fpnearestgrid(iobj,ifp,1)
              ! Test: ix0 in local domain?
              if (ix0 >= l1 .and. ix0 <= l2) then
!
!  Acquire pressure and stress from grid point (ix0,iy0,iz0).
!  Shifting the location of the forcpoints in the thetal direction
!  in order to avoid problems with autotesting
!
                fp_pressure=p%pp(ix0-nghost)
                fp_stress(:,:)=twonu*p%rho(ix0-nghost)*p%sij(ix0-nghost,:,:)
!                
                if (objectform=='cylinder') then
                  longitude = (ifp-theta_shift)*dlong
                  nvec(1) = -sin(longitude)
                  nvec(2) = -cos(longitude)
                  nvec(3) = 0
                elseif (objectform == 'sphere') then
                  ilong  = mod(ifp-1,nlong)
                  ilat   = int((ifp-1)/nlong)
                  longitude = (ilong+.5-theta_shift)*dlong
                  latitude  = (ilat+.5)*dlat
                  nvec(1) = -sin(longitude)*sin(latitude)
                  nvec(2) = -cos(longitude)*sin(latitude)
                  nvec(3) = cos(latitude)
                  surfaceelement = surfacecoeff*sin(latitude)
                end if
!
!  Force in x direction
!
                force_x = (-fp_pressure*nvec(1) &
                    + fp_stress(1,1)*nvec(1) &
                    + fp_stress(1,2)*nvec(2) & 
                    + fp_stress(1,3)*nvec(3)) * surfaceelement
!                
!  Force in y direction
!
                force_y = (-fp_pressure*nvec(2) &
                    + fp_stress(2,1)*nvec(1) &
                    + fp_stress(2,2)*nvec(2) & 
                    + fp_stress(2,3)*nvec(3)) * surfaceelement
!                
!  Force in z direction
!
                force_z = (-fp_pressure*nvec(3) &
                    + fp_stress(3,1)*nvec(1) &
                    + fp_stress(3,2)*nvec(2) & 
                    + fp_stress(3,3)*nvec(3)) * surfaceelement
!
                c_dragx(iobj) = c_dragx(iobj) + force_x
                c_dragy(iobj) = c_dragy(iobj) + force_y
                c_dragz(iobj) = c_dragz(iobj) + force_z
              endif
            endif
          enddo
        enddo
!
!  Calculate average density of the domain, solid cell regions excluded:
!
        do i=l1,l2
          if (mod(ba(i,m,n,1),10) == 0) then
            rhosum = rhosum + p%rho(i-nghost)
            irhocount = irhocount+1
          endif
        enddo
      endif
    endif
!
    call keep_compiler_quiet(df,f)
!
  endsubroutine dsolid_dt
!***********************************************************************  
  subroutine dsolid_dt_integrate
!
!  Calculate drag- and lift-coefficients for solid cell objects
!  by integrating fluid force on object surface. 
!
!  mar-2009/kragset: coded
!  okt-2009/kragset: updated to include multiple objects
!  nov-2010/kragset: updated to include spheres
!
    use mpicomm
    use general

    real    :: rhosum_all, c_dragx_all(nobjects), c_dragy_all(nobjects)
    real    :: c_dragz_all(nobjects)
    integer :: irhocount_all,iobj
    real    :: norm, refrho0
    character*50  :: numberstring
    character*500 :: solid_cell_drag
!    
    if (ldiagnos) then
      if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 &
          .or. idiag_c_dragz /= 0) then 
!
!  Collect and sum rhosum, irhocount, c_dragx, c_dragz, and c_dragy.
        call mpireduce_sum(rhosum,rhosum_all)
        call mpireduce_sum_int(irhocount,irhocount_all)
        call mpireduce_sum(c_dragx,c_dragx_all,nobjects)
        call mpireduce_sum(c_dragy,c_dragy_all,nobjects)
        call mpireduce_sum(c_dragz,c_dragz_all,nobjects)
!        
        if (lroot) then          
          refrho0 = rhosum_all / irhocount_all
!  Normalizing factor. Additional factors was included in subroutine dsolid_dt.
          norm = 1. / (refrho0*init_uu**2)
!          
          c_dragx = c_dragx_all * norm
          c_dragy = c_dragy_all * norm
          c_dragz = c_dragz_all * norm
!          
!  Write drag coefficients for all objects
!  (may need to expand solid_cell_drag to more
!  characters if large number of objects).
! 
          open(unit=81,file='data/dragcoeffs.dat',position='APPEND')
          write(solid_cell_drag,84) it-1, t
          do iobj=1,nobjects
            write(numberstring,82) c_dragx(iobj), c_dragy(iobj),c_dragz(iobj)
            call safe_character_append(solid_cell_drag,numberstring)
          enddo
          write(81,*) trim(solid_cell_drag)
          close(81)
84        format(1I8,1F15.8)
82        format(3F15.8)
        endif
      endif
      if (idiag_c_dragx /= 0) fname(idiag_c_dragx)=c_dragx(1)
      if (idiag_c_dragy /= 0) fname(idiag_c_dragy)=c_dragy(1)
      if (idiag_c_dragz /= 0) fname(idiag_c_dragz)=c_dragz(1)
    endif
 !   
  endsubroutine dsolid_dt_integrate
 !***********************************************************************  
  subroutine rprint_solid_cells(lreset,lwrite)
!
!  Reads and registers print parameters relevant for solid cells
!
!   mar-2009/kragset: coded
!   nov-2010/kragset: generalized to include drag in z-direction
!
!    
    use cdata
    use sub
    use diagnostics
!    
    integer :: iname
    logical :: lreset,lwr
    logical, optional :: lwrite
!
    lwr = .false.
    if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset
! 
    if (lreset) then
      idiag_c_dragx=0 
      idiag_c_dragy=0
      idiag_c_dragz=0
    endif
!
!  check for those quantities that we want to evaluate online
!
    do iname=1,nname
      call parse_name(iname,cname(iname),cform(iname),'c_dragx',idiag_c_dragx)
      call parse_name(iname,cname(iname),cform(iname),'c_dragy',idiag_c_dragy)
      call parse_name(iname,cname(iname),cform(iname),'c_dragz',idiag_c_dragz)
    enddo
!
!  write column, idiag_XYZ, where our variable XYZ is stored
!
    if (lwr) then
!
    endif
!    
  endsubroutine rprint_solid_cells
!***********************************************************************  
    subroutine update_solid_cells(f)
!
!  Set the boundary values of the solid area such that we get a 
!  correct fluid-solid interface.
!
!  19-nov-2008/nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,k,idir,xind,yind,zind,iobj
      
      real :: z_obj, y_obj, x_obj, r_obj, r_new, r_point, sin_theta, cos_theta
      real :: xmirror, ymirror, zmirror, phi, dr, cos_phi, sin_phi
      integer :: lower_i, upper_i, lower_j, upper_j, ii, jj, kk
      integer :: lower_k, upper_k, ndims
      logical :: bax, bay, baz
      real :: gpp
      real, dimension(3) :: xxp
      character(len=10) :: form
!
!  Find ghost points based on the mirror interpolation method
!
      if (interpolation_method=='mirror') then
        do i=l1,l2
        do j=m1,m2
        do k=n1,n2
          bax=(ba(i,j,k,1) /= 0).and.(ba(i,j,k,1)/=9).and.(ba(i,j,k,1)/=10)
          bay=(ba(i,j,k,2) /= 0).and.(ba(i,j,k,2)/=9).and.(ba(i,j,k,2)/=10)
          if (form=='sphere') then
            baz=(ba(i,j,k,3) /= 0).and.(ba(i,j,k,3)/=9).and.(ba(i,j,k,3)/=10)
          else
            baz=.false.
          endif
!
!  Check if we are in a point which must be interpolated, i.e. we are inside
!  a solid geometry AND we are not more than three grid points from the 
!  closest solid-fluid interface
!
          if (bax.or.bay.or.baz) then
!
!  Find x, y and z values of mirror point
!
            iobj=ba(i,j,k,4)
            x_obj=objects(iobj)%x(1)
            y_obj=objects(iobj)%x(2)
            z_obj=objects(iobj)%x(3)
            r_obj=objects(iobj)%r
            form=objects(iobj)%form
            if (form=='cylinder') then
              r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2)
              r_new=r_obj+(r_obj-r_point)
              sin_theta=(y(j)-y_obj)/r_point
              cos_theta=(x(i)-x_obj)/r_point
              xmirror=cos_theta*r_new+x_obj
              ymirror=sin_theta*r_new+y_obj 
              zmirror=0.
            elseif (form=='sphere') then
              r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2+(z(k)-z_obj)**2)
              r_new=r_obj+(r_obj-r_point)
              sin_theta=(y(j)-y_obj)/r_point
              cos_theta=(x(i)-x_obj)/r_point
              cos_phi  =(z(k)-z_obj)/r_point
              sin_phi  =sqrt(1-cos_phi**2)
              xmirror=sin_phi*cos_theta*r_new+x_obj
              ymirror=sin_phi*sin_theta*r_new+y_obj            
              zmirror=cos_phi*r_new+z_obj
!
! Check if mirror point is inside domain
!
              if (xmirror<xyz0(1) .and. lperi(1)) xmirror=xmirror+Lxyz(1)
              if (ymirror<xyz0(2) .and. lperi(2)) ymirror=ymirror+Lxyz(2)
              if (zmirror<xyz0(3) .and. lperi(3)) zmirror=zmirror+Lxyz(3)
              if (xmirror>xyz1(1) .and. lperi(1)) xmirror=xmirror-Lxyz(1)
              if (ymirror>xyz1(2) .and. lperi(2)) ymirror=ymirror-Lxyz(2)
              if (zmirror>xyz1(3) .and. lperi(3)) zmirror=zmirror-Lxyz(3)
             endif
!
!  Check that we are indeed inside the solid geometry
!
            if (r_point>r_obj) then
              print*,'i,j,k=',i,j,k
              print*,'x(i),x_obj=',x(i),x_obj
              print*,'y(j),y_obj=',y(j),y_obj
              print*,'z(k),z_obj=',z(k),z_obj
              print*,'r_point,r_new,r_obj=',r_point,r_new,r_obj
              call fatal_error('update_solid_cells:','r_point>r_obj')
            endif
!
!  Find i, j and k indeces for points to be used during interpolation 
!
            lower_i=0
            upper_i=0
            do ii=1,mx
              if (x(ii)>xmirror) then
                lower_i=ii-1
                upper_i=ii
                exit
              endif
            enddo
!
            lower_j=0
            upper_j=0
            do jj=1,my
              if (y(jj)>ymirror) then
                lower_j=jj-1
                upper_j=jj
                exit
              endif
            enddo
!
            if (form=='sphere') then
              lower_k=0
              upper_k=0
              do kk=1,mz
                if (z(kk)>zmirror) then
                  lower_k=kk-1
                  upper_k=kk
                  exit
                endif
              enddo
            else
              lower_k=k
              upper_k=k
            endif
!
!  Issue with domain borders: A mirror point can be outside a
!  processor's local domain (including ghost points). Some sort
!  communication has to be implemented!
!
            if (lower_i == 0 .or. upper_i == 0) then
              call fatal_error('update_solid_cells:','lower_i==0 or upper_i==0')
            endif
            if (lower_j == 0 .or. upper_j == 0) then
              call fatal_error('update_solid_cells:','lower_j==0 or upper_j==0')
            endif            
            if (form=='sphere') then
              if (lower_k == 0 .or. upper_k == 0) then
                call fatal_error('update_solid_cells:','lower_k==0 or upper_k==0')
              endif
            endif
!
!  First we use interpolations to find the value of the mirror point.
!  Then we use the interpolated value to find the value of the ghost point
!  by empoying either Dirichlet or Neuman boundary conditions.
!
            if (objects(iobj)%form=='cylinder') then
              ndims=2
            elseif (objects(iobj)%form=='sphere') then
              ndims=3
            endif
            call interpolate_mirror_point(f,phi,iux,lower_i,upper_i,lower_j,&
                upper_j,lower_k,upper_k,iobj,xmirror,ymirror,zmirror,ndims)
            f(i,j,k,iux)=-phi
            call interpolate_mirror_point(f,phi,iuy,lower_i,upper_i,lower_j,&
                upper_j,lower_k,upper_k,iobj,xmirror,ymirror,zmirror,ndims)
            f(i,j,k,iuy)=-phi
            call interpolate_mirror_point(f,phi,iuz,lower_i,upper_i,lower_j,&
                upper_j,lower_k,upper_k,iobj,xmirror,ymirror,zmirror,ndims)
            f(i,j,k,iuz)=-phi
            if (ilnrho>0) then
              call interpolate_mirror_point(f,phi,ilnrho,lower_i,upper_i,&
                  lower_j,upper_j,lower_k,upper_k,iobj,xmirror,ymirror,&
                  zmirror,ndims)
              f(i,j,k,ilnrho)=phi
            endif
            if (ilnTT>0) then
              call interpolate_mirror_point(f,phi,ilnTT,lower_i,upper_i,&
                  lower_j,upper_j,lower_k,upper_k,iobj,xmirror,ymirror,&
                  zmirror,ndims)
              f(i,j,k,ilnTT)=2*objects(iobj)%T-phi
            endif
          else
!
!  For fluid points very close to the solid surface the value of the point
!  is set from interpolation between the value at the closest grid line
!  and the value at the solid surface.
!
            if (lclose_linear) then
              if (ba(i,j,k,1)==10) then
                iobj=ba(i,j,k,4)
                x_obj=objects(iobj)%x(1)
                y_obj=objects(iobj)%x(2)
                r_obj=objects(iobj)%r
                r_point=sqrt(((x(i)-x_obj)**2+(y(j)-y_obj)**2))
                dr=r_point-r_obj
                if ((dr > 0) .and. (dr<dxmin*limit_close_linear)) then
                  xxp=(/x(i),y(j),z(k)/)
                  call close_interpolation(f,i,j,k,iobj,iux,xxp,gpp,.true.)
                  f(i,j,k,iux)=gpp
                  call close_interpolation(f,i,j,k,iobj,iuy,xxp,gpp,.true.)
                  f(i,j,k,iuy)=gpp
                  call close_interpolation(f,i,j,k,iobj,iuz,xxp,gpp,.true.)
                  f(i,j,k,iuz)=gpp
                endif
              endif
            endif
          endif
        enddo
        enddo
        enddo
!
!  Find ghost points based on the staircase interpolation method
!
      elseif (interpolation_method=='staircase') then
        do i=l1,l2
        do j=m1,m2
        do k=n1,n2
          do idir=1,3
            if (ba_shift(i,j,k,idir)/=0) then
              xind=i
              yind=j
              zind=k
              if (idir==1) then
                xind=i-ba_shift(i,j,k,idir)
              elseif (idir==2) then
                yind=j-ba_shift(i,j,k,idir)
              elseif (idir==3) then
                zind=k-ba_shift(i,j,k,idir)
              else
                print*,'No such idir!...exiting!'
                stop
              endif
!                
!  Only update the solid cell "ghost points" if all indeces are non-zero.
!  In this way we might loose the innermost "ghost point" if the processor
!  border is two grid cells inside the solid structure, but this will 
!  probably just have a very minor effect.
!
              if (xind/=0 .and. yind/=0 .and. zind/=0) then
                iobj=ba_shift(i,j,k,4)
                f(i,j,k,iux:iuz)=-f(xind,yind,zind,iux:iuz)
                if (ilnrho>0) f(i,j,k,ilnrho) = f(xind,yind,zind,ilnrho)
                if (ilnTT>0) f(i,j,k,ilnTT) = &
                    2*objects(iobj)%T-f(xind,yind,zind,ilnTT)
              endif
            endif
          enddo
        enddo
        enddo
        enddo
      endif
!
    endsubroutine update_solid_cells
!***********************************************************************  
    subroutine interpolate_mirror_point(f,phi,ivar,lower_i,upper_i,lower_j,&
        upper_j,lower_k,upper_k,iobj,xmirror,ymirror,zmirror,ndims)
!
!  Interpolate value in a mirror point from the eight corner values
!
!  23-dec-2008/nils: coded
!  22-apr-2009/nils: added special treatment close to the solid surface
!
      use General, only: linear_interpolate
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: iobj,ndims
      integer, intent(in) :: lower_i,upper_i,lower_j,upper_j,ivar
      integer, intent(in) :: lower_k,upper_k
      real,    intent(in) :: xmirror,ymirror,zmirror
      real,    intent(out):: phi
!
      real, dimension(3) :: xxp
      real, dimension(1) :: gp
      integer, dimension(3) :: inear
!
      xxp=(/xmirror,ymirror,zmirror/)
      inear=(/lower_i,lower_j,lower_k/)
      call linear_interpolate(f,ivar,ivar,xxp,gp,inear,.false.)
      phi=gp(1)
!
!  If the mirror point is very close to the surface of the object 
!  some special treatment is required.
!
      if (lclose_interpolation .and. ivar < 4) then        
        call close_interpolation(f,lower_i,lower_j,lower_k,iobj,ivar,xxp,&
            phi,.false.)
      endif
!
    endsubroutine interpolate_mirror_point
!***********************************************************************  
    subroutine close_interpolation(f,ix0_,iy0_,iz0_,iobj,ivar1,xxp,gpp,&
        fluid_point)
!
!  20-mar-2009/nils: coded
!  
!  If fluid_point=.true. this means that we are handling a grid point.
!  For fluid points very close to the solid surface the value of the point
!  is found from interpolation between the value at the closest grid line "g"
!  and the value at the solid surface "s". 
!  
!  If fluid_point=.false. we are handling a point which is NOT a grid point,
!  i.e. the point we are interested in are somewhere between the grid points.
!  This situation typically appears for mirror points used 
!  to find the value of the ghost points INSIDE the solid geometry, or when
!  a particle is very close to the surface.
!  If fluid_point=.false. the routine check if any of the neighbouring
!  grid points, used for interpolation, are inside the solid geometry.
!  If so the routine use the value at the surface
!  of the solid geometry together with the interpolated value at the nearest
!  grid line in the direction away from the solid geometry to set a value
!  at a grid point which is very close to the solid geometry.
!  
!  The interpolation point is called "p".
!  Define a straight line "l" which pass through the point "p" and is normal 
!  to the surface of the object.
!  The point where "l" cross the surface of the object is called "s"
!  The point "g" is the point where "l" cross the first grid plane (line in 2D)
!  outside of "p".
!  The point "object" is the center point of the solid object.
!  The variable "bv" give the x,y and z values of the neighbouring grid points
!  to "p".
!  
!  
!---------------------------------------------------------------------------
! Point   Description               Global coord. sys.    Local coord. syst.
!                                                         (origo in center 
!                                                         of object)
!---------------------------------------------------------------------------
! p       The interpolated point    p_global(3)           p_local(3)
! s       "l" cross surface         -                     -
! g       "l" cross grid plane      g_global(3)           -
! object  Object center             o_global(3)           -
! bv      Border value              bordervalue(3*2)      -
!---------------------------------------------------------------------------
!
      use General, only: linear_interpolate
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: ix0_,iy0_,iz0_,ivar1
      integer :: ix0,iy0,iz0
      real, intent(inout) :: gpp
      real :: rs,verylarge=1e9,varval,rint1,rint2,fint,rps,rintp,phi
      real :: r_pg,r_sg,rl,rlmin,xmirror,ymirror,zmirror
      integer :: ix1,iy1,iz1,min,dir,lower_i,lower_j,lower_k,ndims,ndir
      integer :: upper_i,upper_j,upper_k,vardir1,vardir2
      real, dimension(3), intent(in) :: xxp
      real, dimension(2,2,2) :: rij
      real, dimension(3,2) :: bordervalue
      integer, dimension(3,2) :: borderindex
      integer, dimension(6) :: constdir_arr, vardir_arr, topbot_arr
      real, dimension(3) :: xyint, p_global, p_local, o_global
      real, dimension(3) :: ngrids, g_global
      real :: xtemp,r,R1,Rsmall,xs,ys,zs,rp,dist
      integer :: constdir,vardir,topbot_tmp,dirconst,dirvar,iobj,counter,topbot
      integer :: maxcounter
      real :: x1,x2,f1,f2,rij_min,rij_max,inputvalue,smallx
      logical, intent(in) :: fluid_point
      real :: fintx, finty,fint_ur,fint_ut,drp,dri,f2x,f2y,f1x,f1y
      real, save :: urp,utp
      logical :: quadratic=.true.
      real, dimension(1) :: gp
      integer, dimension(3) :: inear
!
!  Check if we really want this special treatment close to the fluid-solid 
!  interface
!
      if ((.not. fluid_point .and. lclose_interpolation) &
          .or. ( fluid_point .and. lclose_linear)) then
!
!  This subrutine is not working (and should never be used) with other
!  variables than the velocity.
!
        if (ivar1 > iuz) call fatal_error('close_interpolation',&
            'This subroutine should never be called for anything but velocity!')
!
!  Define some help variables
!
        o_global=objects(iobj)%x
        rs=objects(iobj)%r
        p_global=xxp
!
!  Find the corner points of the grid cell we are in
!
        if (fluid_point) then
          smallx=dx*1e-5
          iz0=iz0_
          if (p_global(1) < o_global(1)) then
            ix0=ix0_-1
            p_global(1)=p_global(1)-smallx
          else
            ix0=ix0_
            p_global(1)=p_global(1)+smallx
          endif
          if (p_global(2) < o_global(2)) then
            iy0=iy0_-1
            p_global(2)=p_global(2)-smallx
          else
            iy0=iy0_
            p_global(2)=p_global(2)+smallx
          endif
          if (p_global(3) < o_global(3)) then
            iz0=iz0_-1
            p_global(3)=p_global(3)-smallx
          else
            iz0=iz0_
            p_global(3)=p_global(3)+smallx
          endif
        else
          ix0=ix0_
          iy0=iy0_
          iz0=iz0_
        endif
        ix1=ix0+1
        iy1=iy0+1
        iz1=iz0+1
!
!  Find distance from corner points to the object center
!
        rij(1,1,1)=sqrt((x(ix0)-o_global(1))**2+(y(iy0)-o_global(2))**2&
            +(z(iz0)-o_global(3))**2)
        rij(1,2,1)=sqrt((x(ix0)-o_global(1))**2+(y(iy1)-o_global(2))**2&
            +(z(iz0)-o_global(3))**2)
        rij(2,1,1)=sqrt((x(ix1)-o_global(1))**2+(y(iy0)-o_global(2))**2&
            +(z(iz0)-o_global(3))**2)
        rij(2,2,1)=sqrt((x(ix1)-o_global(1))**2+(y(iy1)-o_global(2))**2&
            +(z(iz0)-o_global(3))**2) 
        rij(1,1,2)=sqrt((x(ix0)-o_global(1))**2+(y(iy0)-o_global(2))**2&
            +(z(iz1)-o_global(3))**2)
        rij(1,2,2)=sqrt((x(ix0)-o_global(1))**2+(y(iy1)-o_global(2))**2&
            +(z(iz1)-o_global(3))**2)
        rij(2,1,2)=sqrt((x(ix1)-o_global(1))**2+(y(iy0)-o_global(2))**2&
            +(z(iz1)-o_global(3))**2)
        rij(2,2,2)=sqrt((x(ix1)-o_global(1))**2+(y(iy1)-o_global(2))**2&
            +(z(iz1)-o_global(3))**2) 
!
!  Check if we want special treatment
!
        if ((minval(rij) < rs) .or. fluid_point) then
!
!  Put help variables into arrays
!
          bordervalue(1,1)=x(ix0)
          bordervalue(2,1)=y(iy0)
          bordervalue(3,1)=z(iz0)
          bordervalue(1,2)=x(ix1)
          bordervalue(2,2)=y(iy1)
          bordervalue(3,2)=z(iz1)
          borderindex(1,1)=ix0
          borderindex(2,1)=iy0
          borderindex(3,1)=iz0
          borderindex(1,2)=ix1
          borderindex(2,2)=iy1
          borderindex(3,2)=iz1
          if (objects(iobj)%form=='cylinder') then
            constdir_arr=(/2,2,1,1,0,0/)
            vardir_arr  =(/1,1,2,2,0,0/)
            topbot_arr  =(/2,1,2,1,0,0/)
          elseif (objects(iobj)%form=='sphere') then
            constdir_arr=(/3,3,2,2,1,1/)
            vardir_arr  =(/1,1,1,1,3,3/)
            topbot_arr  =(/2,1,2,1,2,1/)
          endif
          R1=verylarge
          Rsmall=verylarge/2.0
!
!  Find the x, y and z coordinates of "p" in a coordiante system with origin
!  in the center of the object
!
          p_local=p_global-o_global
!
!  Find the distance from the object center to "p"
!
          if (objects(iobj)%form=='cylinder') then
            rp=sqrt(p_local(1)**2+p_local(2)**2)
            p_local(3)=0
          elseif (objects(iobj)%form=='sphere') then
            rp=sqrt(p_local(1)**2+p_local(2)**2+p_local(3)**2)
          endif
!
!  Currently there are two implementations for the close surface treatment.
!  The first one is more general and works both for spheres and cylinders,
!  but it does not have quadratic interpolation of radial velocities
!  implemented yet. Furthermore it does not handle particles correctly
!  for the cases where one of the corner points of "gridplane" is insed
!  the solid geometry.
!
          if (objects(iobj)%form=='sphere' .or. .not. quadratic) then
            if (lparticles .and. lclose_interpolation) then
              print*,'Using close interpolation together with particles '
              print*,'does not yet work with spherical solid geometries.'
              print*,'You should either:'
              print*,' -Turn of lclose_interpolation'
              print*,' -Use cylinders instead of spheres'
              print*,' -Implement what is missing in order to solve the problem'
              call fatal_error('close_interpolation','')
            endif
!
!  Find which grid line is the closest one in the direction
!  away from the object surface
!
!  Check the distance to all the six (four in 2D) possible 
!  surface planes (lines in 2D). 
!  Define a straight line "l" which pass through the point "p" and is normal 
!  to the surface of the object.
!  Pick the grid plane, OUTSIDE the point "p", that the line "l" cross first,
!  call this plane "gridplane". 
!
            ngrids(1)=nxgrid
            ngrids(2)=nygrid
            ngrids(3)=nzgrid
            ndir=3
            ndims=0
            rlmin=verylarge
            do dir=1,ndir
! 
!  Loop through all directions "dir" and find which grid plane, outside of "p",
!  that is crossed by the line "l" first. Lets call this plane "gridplane". 
!  If "dir" for the first crossing is 1 then "gridplane" is the 'yz' plane 
!  and so on.....
!
              if (ngrids(dir) .gt. 1) then
                ndims=ndims+1
                do topbot_tmp=1,2
!
!  Find the variable "rl" such that
!  rl*p_local(dir)+o_global(dir)=bordervalue(dir)  
!
                  rl=(bordervalue(dir,topbot_tmp)-o_global(dir))/p_local(dir)
!
!  If "rl" is larger than unity the line "l" cross the grid plane
!  outside of "p". If in addition "rl" is smaller than the smallest "rl" so
!  far (rlmin) this state must be stored since it might be the plane
!  ("gridplane") which we are looking for. In addition we must
!  find the distance, r, from the center of the object
!  to the point where the normal cross the grid line.
!  The point where "l" cross "gridplane" is called "g".
!
                  if (rl > 1.0 .and. rl<rlmin) then
                    constdir=dir
                    vardir1=mod(constdir+0,3)+1
                    vardir2=mod(constdir+1,3)+1
                    rlmin=rl
                    topbot=topbot_tmp
                    r=rl*sqrt(p_local(1)**2+p_local(2)**2+p_local(3)**2)
                    g_global(constdir)=bordervalue(constdir,topbot)
                    g_global(vardir1)=rl*p_local(vardir1)+o_global(vardir1)
                    g_global(vardir2)=rl*p_local(vardir2)+o_global(vardir2)
                  endif
                enddo
              endif
            enddo
!
!  The direction of the normal to "gridplane", called "N_grid", is either
!  in the x, y or z direction since the grid is Cartesian. See table below
!  for values of constdir, vardir1, and vardir2 for different "N_grid" directions.
!
!  N_grid  constdir  vardir1  vardir2 
!--------------------------------------------------
!     x        1        2        3    
!     y        2        3        1
!     z        3        1        2
!
!  Check that we have found a valid distance
!
            if (rlmin==verylarge) then
              print*,'fluid_point=',fluid_point
              print*,'lclose_interpolation=',lclose_interpolation
              print*,'lclose_linear=',lclose_linear
              print*,'o_global=',o_global
              print*,'ix0,iy0,iz0=',ix0,iy0,iz0
              print*,'ix1,iy1,iz1=',ix1,iy1,iz1
              print*,'p_global=',p_global
              print*,'r,rs,rp=',r,rs,rp
              print*,'rij,rs=',rij,rs
              print*,'x(ix0),p_global(1),x(ix1)=',x(ix0),p_global(1),x(ix1)
              print*,'y(iy0),p_global(2),y(iy1)=',y(iy0),p_global(2),y(iy1)
              print*,'dirvar,dirconst,topbot,iz0=',dirvar,dirconst,topbot,iz0
              call fatal_error('close_interpolation',&
                  'A valid radius is not found!')            
            endif
!
!  Depending on the value of constdir the indeces of the corner points 
!  specifying "gridplane" is given by lower_i, upper_i, lower_j ........
!  The physical position of "g" is stored in xmirror, ymirror and zmirror.
!
            if (constdir==1) then
              lower_i=borderindex(constdir,topbot)
              upper_i=borderindex(constdir,topbot)
              lower_j=borderindex(vardir1,1)
              upper_j=borderindex(vardir1,2)
              lower_k=borderindex(vardir2,1)
              upper_k=borderindex(vardir2,2)
              xmirror=g_global(constdir)
              ymirror=g_global(vardir1)
              zmirror=g_global(vardir2)
            elseif (constdir==2) then
              lower_j=borderindex(constdir,topbot)
              upper_j=borderindex(constdir,topbot)
              lower_k=borderindex(vardir1,1)
              upper_k=borderindex(vardir1,2)
              lower_i=borderindex(vardir2,1)
              upper_i=borderindex(vardir2,2)
              ymirror=g_global(constdir)
              xmirror=g_global(vardir2)
              zmirror=g_global(vardir1)
            elseif (constdir==3) then
              lower_k=borderindex(constdir,topbot)
              upper_k=borderindex(constdir,topbot)
              lower_i=borderindex(vardir1,1)
              upper_i=borderindex(vardir1,2)
              lower_j=borderindex(vardir2,1)
              upper_j=borderindex(vardir2,2)
              zmirror=g_global(constdir)
              xmirror=g_global(vardir1)
              ymirror=g_global(vardir2)
            endif
!
!  Let "phi" be the physical value of the variable "ivar1" on "gridplane".
!  The value of "phi" is then found by interpolation between the four corner
!  points of "gridplane".
!
            inear=(/lower_i,lower_j,lower_k/)
            call linear_interpolate(f,ivar1,ivar1,g_global,gp,inear,.false.)
            phi=gp(1)
!
!  Now we know the value associated with the variable "ivar1" in the point "g".
!  Furthermore we know the value value associated the "ivar1" in point "s"
!  on the object surface to be zero for any of the velocities and equal to 
!  the solid temperature for the temperature.
!  By interpolation it is now straight forward to find the value also in "p".
!  First find the distance, "r_pg", from "p" to "g" to be used as weight, 
!  then normalize by the full distance, "r_sg", between "s" and "g". 
!
            r_pg=r-rp
            r_sg=r-rs
            gpp=(1-r_pg/r_sg)*phi
          else
            maxcounter=6
            if (objects(iobj)%form=='cylinder') maxcounter=4
            do counter=1,maxcounter
              constdir=constdir_arr(counter)
              vardir=vardir_arr(counter)
              topbot_tmp=topbot_arr(counter)
!
!  Find the position, xtemp, in the variable direction
!  where the normal cross the grid line
!
              xtemp=(p_local(vardir)/(p_local(constdir)+tini))&
                  *(bordervalue(constdir,topbot_tmp)-objects(iobj)%x(constdir))
!
!  Find the distance, r, from the center of the object
!  to the point where the normal cross the grid line
!
              if (abs(xtemp) > verylarge) then
                r=verylarge*2
              else
                r=sqrt(xtemp**2+(bordervalue(constdir,topbot_tmp)&
                    -objects(iobj)%x(constdir))**2)
              endif
!
!  Check if the point xtemp is outside the object,
!  outside the point "p" and that it cross the grid
!  line within this grid cell
!
              if ((r > rs) .and. (r > rp) &
                  .and.(xtemp+objects(iobj)%x(vardir)>=bordervalue(vardir,1))&
                  .and.(xtemp+objects(iobj)%x(vardir)<=bordervalue(vardir,2)))then
                R1=r
              else
                R1=verylarge
              endif
!
!  If we have a new all time low (in radius) then go on....
!
              if (R1 < Rsmall) then
                Rsmall=R1                
                xyint(vardir)=xtemp+objects(iobj)%x(vardir)
                xyint(constdir)=bordervalue(constdir,topbot_tmp)
                dirconst=constdir
                dirvar=vardir
                topbot=topbot_tmp
                if (constdir == 2) then
                  rij_min=rij(1,topbot_tmp,1)
                  rij_max=rij(2,topbot_tmp,1)
                else
                  rij_min=rij(topbot_tmp,1,1)
                  rij_max=rij(topbot_tmp,2,1)
                endif
                inputvalue=bordervalue(constdir,topbot_tmp)&
                    -objects(iobj)%x(constdir)
              endif
            enddo
!
!  Check that we have found a valid distance
!
            if (Rsmall==verylarge/2.0) then
              print*,'fluid_point=',fluid_point
              print*,'lclose_interpolation=',lclose_interpolation
              print*,'lclose_linear=',lclose_linear
              print*,'o_global=',o_global
              print*,'ix0,iy0,iz0=',ix0,iy0,iz0
              print*,'ix1,iy1,iz1=',ix1,iy1,iz1
              print*,'p_global=',p_global
              print*,'xtemp=',xtemp
              print*,'r,rs,rp=',r,rs,rp
              print*,'R1,Rsmall=',R1,Rsmall
              print*,'rij,rs=',rij,rs
              print*,'x(ix0),p_global(1),x(ix1)=',x(ix0),p_global(1),x(ix1)
              print*,'y(iy0),p_global(2),y(iy1)=',y(iy0),p_global(2),y(iy1)
              print*,'dirvar,dirconst,topbot,iz0=',dirvar,dirconst,topbot,iz0
              call fatal_error('close_interpolation',&
                  'A valid radius is not found!')            
            endif
!
!  Check if the endpoints in the variable direction are
!  outside the objects. If they are not then define the endpoints
!  as where the grid line cross the object surface.
!  Find the variable value at the endpoints.
!
            min=1
            if (dirconst == 2) then
              varval=f(borderindex(dirvar,1),borderindex(dirconst,topbot),&
                  iz0,ivar1)
            else
              varval=f(borderindex(dirconst,topbot),borderindex(dirvar,1),&
                  iz0,ivar1)
            endif
            call find_point(rij_min,rs,varval,inputvalue,x1,&
                bordervalue(dirvar,1),bordervalue(dirvar,2),&
                min,f1,objects(iobj)%x(dirvar))
!
! If we want quadratic interpolation of the radial velocity we
! must find both the interploated x and y velocity in order to 
! do interpolations for the radial and theta directions.
!
            if (quadratic .and. ivar1==iux) then
              if (dirconst == 2) then
                varval=f(borderindex(dirvar,1),borderindex(dirconst,topbot),&
                    iz0,iuy)
              else
                varval=f(borderindex(dirconst,topbot),borderindex(dirvar,1),&
                    iz0,iuy)
              endif
              call find_point(rij_min,rs,varval,inputvalue,x1,&
                  bordervalue(dirvar,1),bordervalue(dirvar,2),&
                  min,f1y,objects(iobj)%x(dirvar))
              f1x=f1
            endif
!
            min=0
            if (dirconst == 2) then
              varval=f(borderindex(dirvar,2),borderindex(dirconst,topbot),&
                  iz0,ivar1)
            else
              varval=f(borderindex(dirconst,topbot),borderindex(dirvar,2),&
                  iz0,ivar1)
            endif
            call find_point(rij_max,rs,varval,inputvalue,x2,&
                bordervalue(dirvar,1),bordervalue(dirvar,2),&
                min,f2,objects(iobj)%x(dirvar))
!
! If we want quadratic interpolation of the radial velocity we
! must find both the interploated x and y velocity in order to 
! do interpolations for the radial and theta directions.
!
            if (quadratic .and. ivar1==iux) then
              if (dirconst == 2) then
                varval=f(borderindex(dirvar,2),borderindex(dirconst,topbot),&
                    iz0,iuy)
              else
                varval=f(borderindex(dirconst,topbot),borderindex(dirvar,2),&
                    iz0,iuy)
              endif
              call find_point(rij_max,rs,varval,inputvalue,x2,&
                  bordervalue(dirvar,1),bordervalue(dirvar,2),&
                  min,f2y,objects(iobj)%x(dirvar))
              f2x=f2
            endif
!
!  Find the interpolation values between the two endpoints of
!  the line and the normal from the objects.
!
            rint1=xyint(dirvar)-x1
            rint2=x2-xyint(dirvar)
!
            if (quadratic .and. (ivar1 /= iuz)) then
              if (ivar1==iux) then
                fintx=(rint1*f2x+rint2*f1x)/(x2-x1)
                finty=(rint1*f2y+rint2*f1y)/(x2-x1)
                fint_ur    =fintx*p_local(1)/rs+finty*p_local(2)/rs
                fint_ut=finty*p_local(1)/rs-fintx*p_local(2)/rs
                drp=rp-rs
                dri=Rsmall-rs
                urp=(drp/dri)**2*fint_ur
                utp=(drp/dri)*fint_ut
                gpp=urp*p_local(1)/rs-utp*p_local(2)/rs
              elseif (ivar1==iuy) then
                gpp=urp*p_local(2)/rs+utp*p_local(1)/rs
              else
                call fatal_error('close_interpolation',&
                    'Yor ivar1 is not correct!') 
              endif
            else
!
!  Find the interpolated value on the line
!
              fint=(rint1*f2+rint2*f1)/(x2-x1)
!
!  Find the weigthing factors for the point on the line
!  and the point on the object surface.
!          
              rps=rp-rs
              rintp=Rsmall-rp
!
!  Perform the final interpolation
!
              gpp=(rps*fint+rintp*0)/(Rsmall-rs)
            endif
          endif
        endif
      endif
!
    endsubroutine close_interpolation
!***********************************************************************  
    function in_solid_cell(part_pos,part_rad)
!
!  Check if the position px,py,pz is within a colid cell
!
!  02-dec-2008/nils: coded
!
      logical :: in_solid_cell
      real, dimension(3) :: obj_pos, part_pos
      real :: obj_rad,distance2,part_rad,rad_part
      integer :: iobj, i, ndims
!
      in_solid_cell=.false.
!
      do iobj=1,nobjects
        obj_rad=objects(iobj)%r
        obj_pos=objects(iobj)%x(1:3)
        distance2=0
!
!  Loop only over the number of dimensions required
!
        ndims=2
        if (objects(iobj)%form=='sphere') ndims=3
        do i=1,ndims
          distance2=distance2+(obj_pos(i)-part_pos(i))**2
        enddo
!
!  Check if we want to include interception or not
!
        if (lnointerception) then
          rad_part=0
        else
          rad_part=part_rad
        endif
!
!  The object_skin is the closest a particle can get to the solid 
!  cell before it is captured (this variable is normally zero).
!
        if (sqrt(distance2)<obj_rad+rad_part+object_skin) then
          in_solid_cell=.true.
        endif
      enddo
!
    endfunction in_solid_cell
!***********************************************************************  
    subroutine freeze_solid_cells(df)
!
!  If we are in a solid cell (or in a cell where the value of the variables are
!  found from interpolation) set df=0 for all variables
!
!  19-nov-2008/nils: coded
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: i
!
      do i=l1,l2
        if (&
            (ba(i,m,n,1)/=0).or.&
            (ba(i,m,n,2)/=0).or.&
            (ba(i,m,n,3)/=0)) then
!
!  If this is a fluid point which has to be interpolated because it is very
!  close to the solid geometry (i.e. ba(i,m,n,1) == 10) then only the 
!  velocity components should be frozen.
!
          if (ba(i,m,n,1) == 10) then
            df(i,m,n,iux:iuz)=0
          else
            df(i,m,n,:)=0
          endif
        endif
      enddo
!
    endsubroutine freeze_solid_cells
!***********************************************************************
    subroutine read_solid_cells_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=solid_cells_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=solid_cells_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_solid_cells_init_pars
!***********************************************************************
    subroutine read_solid_cells_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=solid_cells_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=solid_cells_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_solid_cells_run_pars
!***********************************************************************
    subroutine write_solid_cells_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=solid_cells_init_pars)

    endsubroutine write_solid_cells_init_pars
!***********************************************************************
    subroutine write_solid_cells_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=solid_cells_run_pars)

    endsubroutine write_solid_cells_run_pars
!***********************************************************************
    subroutine find_solid_cell_boundaries
!
!  Find the boundaries of the geometries such that we can set the
!  ghost points inside the solid geometry in order to achieve the
!  correct no-slip boundaries.
!  
!  Store data in the ba array.
!  If ba(ip,jp,kp,1)= 0 we are in a fluid cell (i.e. NOT inside a solid geometry)
!  If ba(ip,jp,kp,1)=10 we are in a fluid cell which are so close to the 
!                       surface of the solid geometry that we set the value of
!                       this point by interpolating between the value at the 
!                       solid surface and the interpolated value at the first
!                       grid line crossed by the normal to the solid surface.
!  If ba(ip,jp,kp,1)= 9 we are inside a solid geometry, but far from the boundary
!  If ba(ip,jp,kp,1)=-1 we are inside a solid geometry, and the point at ip+1
!                       is outside the geometry. 
!  If ba(ip,jp,kp,1)=-3 we are inside a solid geometry, and the point at ip+3
!                       is outside the geometry. 
!  If ba(ip,jp,kp,2)=-3 we are inside a solid geometry, and the point at jp+3
!                       is outside the geometry. 
!  If ba(ip,jp,kp,2)=11 we are inside a solid geometry, either close to or far 
!                       from the boundary, but the position (ip,jp,kp) is a ghost
!                       point at the current processor.
!
!  The number stored in ba(ip,jp,kp,4) is the number of the object
!
!  19-nov-2008/nils: coded
!
      integer :: i,j,k,iobj,cw
      real :: x2,y2,z2,xval_p,xval_m,yval_p,yval_m, zval_p,zval_m
      real :: dr,r_point,x_obj,y_obj,z_obj,r_obj
      character(len=10) :: form
!
!  Initialize ba
!
      ba=0
!
!  Loop over all objects 
!
      do iobj=1,nobjects
        x_obj=objects(iobj)%x(1)
        y_obj=objects(iobj)%x(2)
        z_obj=objects(iobj)%x(3)
        r_obj=objects(iobj)%r
        form=objects(iobj)%form
!
!  First we look in x-direction
!
        do k=n1,n2
        do j=m1,m2
!
!  Check if we are inside the object for y(j) and z(k) (i.e. if x2>0)
!  This depens on the form of the solid geometry
!
          if (form=='cylinder') then
            x2&
                =objects(iobj)%r**2&
                -(y(j)-objects(iobj)%x(2))**2
          else if (form=='sphere') then
            x2&
                =objects(iobj)%r**2&
                -(y(j)-objects(iobj)%x(2))**2&
                -(z(k)-objects(iobj)%x(3))**2
          else
            call fatal_error('find_solid_cell_boundaries','No such form!')
          endif
          if (x2>0) then
!
!  Find upper and lower x-values for the surface of the object for y(j) and z(k)
!
            xval_p=objects(iobj)%x(1)+sqrt(x2)
            xval_m=objects(iobj)%x(1)-sqrt(x2)
            do i=l1,l2
              if (x(i)<xval_p .and. x(i)>xval_m) then
                !
                if (x(i+1)>xval_p) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=-1
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==1) ba(i,j,k,1)=-1
                  endif
                endif
                !
                if (x(i+2)>xval_p .and. x(i+1)<xval_p) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=-2
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==1) ba(i,j,k,1)=-2
                  endif
                endif
                !
                if (x(i+3)>xval_p .and. x(i+2)<xval_p) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=-3
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==1) ba(i,j,k,1)=-3
                  endif
                endif
                !
                if (x(i-1)<xval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=1
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-1) ba(i,j,k,1)=1
                  endif
                endif
                !
                if (x(i-2)<xval_m .and. x(i-1)>xval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=2
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-1) ba(i,j,k,1)=2
                  endif
                endif
                !
                if (x(i-3)<xval_m .and. x(i-2)>xval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=3
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-1) ba(i,j,k,1)=3
                  endif
                endif
                !
                if (ba(i,j,k,1)==0) then
                  ba(i,j,k,1)=9
                  ba(i,j,k,4)=iobj
                endif
                !
              endif
            enddo
          endif
        enddo
        enddo
!
!  Then we look in y-direction
!
        do k=n1,n2
        do i=l1,l2
!
!  Check if we are inside the object for x(i) (i.e. if y2>0)
!  This depens on the form of the solid geometry
!
          if (form=='cylinder') then
            y2&
                =objects(iobj)%r**2&
                -(x(i)-objects(iobj)%x(1))**2
          else if (form=='sphere') then
            y2&
                =objects(iobj)%r**2&
                -(x(i)-objects(iobj)%x(1))**2&
                -(z(k)-objects(iobj)%x(3))**2
          else
            call fatal_error('find_solid_cell_boundaries','No such form!')
          endif
          if (y2>0) then
!
!  Find upper and lower y-values for the surface of the object for x(i)
!
            yval_p=objects(iobj)%x(2)+sqrt(y2)
            yval_m=objects(iobj)%x(2)-sqrt(y2)            
            do j=m1,m2
              if (y(j)<yval_p .and. y(j)>yval_m) then
                if (y(j+1)>yval_p) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,2)=-1
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==2) ba(i,j,k,2)=-1                  
                  endif
                endif
!
                if (y(j+2)>yval_p .and. y(j+1)<yval_p) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,2)=-2
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==2) ba(i,j,k,2)=-2                  
                  endif
                endif
!
                if (y(j+3)>yval_p .and. y(j+2)<yval_p) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,2)=-3
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==2) ba(i,j,k,2)=-3                  
                  endif
                endif
!
                if (y(j-1)<yval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,2)=1
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-2) ba(i,j,k,2)=1                  
                  endif
                endif
!
                if (y(j-2)<yval_m .and. y(j-1)>yval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,2)=2
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-2) ba(i,j,k,2)=2                  
                  endif
                endif
!
                if (y(j-3)<yval_m .and. y(j-2)>yval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,2)=3
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-2) ba(i,j,k,2)=3                  
                  endif
                endif
!
                if (ba(i,j,k,2)==0) then
                  ba(i,j,k,2)=9
                  ba(i,j,k,4)=iobj
                endif
              endif
            enddo
          endif
        enddo
        enddo
!
!  If form='sphere' we must also look in the z-direction
!
        if (form .ne. 'cylinder') then
        do i=l1,l2
        do j=m1,m2
!
!  Check if we are inside the object for y(j) and x(i) (i.e. if z2>0)
!
          if (form=='cylinder') then
            call fatal_error('find_solid_cell_boundaries',&
                'no cylinders when variable z')
          else if (form=='sphere') then
            z2&
                =objects(iobj)%r**2&
                -(y(j)-objects(iobj)%x(2))**2&
                -(x(i)-objects(iobj)%x(1))**2
          else
            call fatal_error('find_solid_cell_boundaries','No such form!')
          endif
          if (z2>0) then
!
!  Find upper and lower x-values for the surface of the object for y(j) and z(k)
!
            zval_p=objects(iobj)%x(3)+sqrt(z2)
            zval_m=objects(iobj)%x(3)-sqrt(z2)            
            do k=n1,n2
              if (z(k)<zval_p .and. z(k)>zval_m) then
                !
                if (z(k+1)>zval_p) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=-1
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==1) ba(i,j,k,1)=-1
                  endif
                endif
                !
                if (z(k+2)>zval_p .and. z(k+1)<zval_p) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=-2
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==1) ba(i,j,k,1)=-2
                  endif
                endif
                !
                if (z(k+3)>zval_p .and. z(k+2)<zval_p) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=-3
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==1) ba(i,j,k,1)=-3
                  endif
                endif
                !
                if (z(k-1)<zval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=1
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-1) ba(i,j,k,1)=1
                  endif
                endif
                !
                if (z(k-2)<zval_m .and. z(k-1)>zval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=2
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-1) ba(i,j,k,1)=2
                  endif
                endif
                !
                if (z(k-3)<zval_m .and. z(k-2)>zval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,1)=3
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-1) ba(i,j,k,1)=3
                  endif
                endif
                !
                if (ba(i,j,k,1)==0) then
                  ba(i,j,k,1)=9
                  ba(i,j,k,4)=iobj
                endif
                !
              endif
            enddo
          endif
        enddo
        enddo
        endif
!
!  If we want to interpolate points which are very close to the solid surface
!  these points have to be "marked" for later use.
!
        if (lclose_linear) then
!
!  Loop over all points
!
          do i=l1,l2
            do j=m1,m2
              do k=n1,n2
                if (form=='cylinder') then 
                  r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2)
                elseif (form=='sphere') then
                  r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2+(z(k)-z_obj)**2)
                else
                  call fatal_error('find_solid_cell_boundaries','No such form!')
                endif
                dr=r_point-r_obj
                if ((dr >= 0) .and. (dr<limit_close_linear*dxmin)) then
                  ba(i,j,k,1)=10
                  ba(i,j,k,4)=iobj
                endif
              enddo
            enddo
          enddo
        endif
!
!  Fill ba array also for ghost points - need only know whether
!  we are actually inside object (then ba = 11), not how close we are to 
!  the border.
!
! Lower and upper ghost points in z direction
!
        do i=1,mx
        do j=1,my
        do k=1,nghost
          if (form=='cylinder') then 
            ba(i,j,k,1:3)=11
            ba(i,j,k,4)=iobj
            ba(i,j,mz-nghost+k,1:3)=11
            ba(i,j,mz-nghost+k,4)=iobj
          elseif (form=='sphere') then
            r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2+(z(k)-z_obj)**2)
            !  Lower (left) ghost points
            if (r_point < r_obj) then
              ba(i,j,k,1:3)=11
              ba(i,j,k,4)=iobj
            endif
            !  Upper (right) ghost points
            if (r_point < r_obj) then
              ba(i,j,mz-nghost+k,1:3)=11
              ba(i,j,mz-nghost+k,4)=iobj
            endif
          else
            call fatal_error('find_solid_cell_boundaries','No such form!')
          endif
        enddo
        enddo
        enddo
!
!  Lower and upper ghost points in y direction
!
        do j=1,nghost
        do k=1,mz
        do i=1,mx
            !  Lower ghost points
            if (form=='cylinder') then 
              r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2)
            elseif (form=='sphere') then
              r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2+(z(k)-z_obj)**2)
            else
              call fatal_error('find_solid_cell_boundaries','No such form!')
            endif
            if (r_point < r_obj) then
              ba(i,j,k,1:3)=11
              ba(i,j,k,4)=iobj
            endif
            !  Upper ghost points
            if (form=='cylinder') then 
              r_point=sqrt((x(i)-x_obj)**2+(y(my-nghost+j)-y_obj)**2)
            elseif (form=='sphere') then
              r_point=sqrt((x(i)-x_obj)**2+(y(my-nghost+j)-y_obj)**2&
                  +(z(k)-z_obj)**2)
            else
              call fatal_error('find_solid_cell_boundaries','No such form!')
            endif
            if (r_point < r_obj) then
              ba(i,my-nghost+j,k,1:3)=11
              ba(i,my-nghost+j,k,4)=iobj
            endif
        enddo
        enddo
        enddo
!
! Lower and upper ghost points in x direction
!
        do k=1,mz
        do j=1,my
        do i=1,nghost
          !  Lower (left) ghost points
          if (form=='cylinder') then 
            r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2)
          elseif (form=='sphere') then
            r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2+(z(k)-z_obj)**2)
          else
            call fatal_error('find_solid_cell_boundaries','No such form!')
          endif
          if (r_point < r_obj) then
            ba(i,j,k,1:3)=11
            ba(i,j,k,4)=iobj
          endif
          !  Upper (right) ghost points
          if (form=='cylinder') then 
            r_point=sqrt((x(mx-nghost+i)-x_obj)**2+(y(j)-y_obj)**2)
          elseif (form=='sphere') then
            r_point=sqrt((x(mx-nghost+i)-x_obj)**2+(y(j)-y_obj)**2&
                +(z(k)-z_obj)**2)
          else
            call fatal_error('find_solid_cell_boundaries','No such form!')
          endif

          if (r_point < r_obj) then
            ba(mx-nghost+i,j,k,1:3)=11
            ba(mx-nghost+i,j,k,4)=iobj
          endif
        enddo
        enddo
        enddo
!
! Finalize loop over all objects
!
      enddo
!
    endsubroutine find_solid_cell_boundaries
!***********************************************************************
    subroutine calculate_shift_matrix
!
!  Set up the shift matrix
!
!  19-nov-2008/nils: coded
!
      integer :: i,j,k,idir
      integer :: sgn
!
      ba_shift=0
!
      do i=l1,l2
      do j=m1,m2
      do k=n1,n2
        do idir=1,3
!
!  If ba is non-zero find the shift matrix
!
          if (ba(i,j,k,idir)/=0 .and. ba(i,j,k,idir)/=9) then
            sgn=-ba(i,j,k,idir)/abs(ba(i,j,k,idir))
            ba_shift(i,j,k,idir)=2*ba(i,j,k,idir)+sgn
            ba_shift(i,j,k,4)=ba(i,j,k,4)
          endif
        enddo
      enddo
      enddo
      enddo
!
    endsubroutine calculate_shift_matrix
!***********************************************************************  
    subroutine find_closest_wall(i,j,k,iobj,cw)
!
!  Find the direction of the closest wall for given grid point and object
!
!  28-nov-2008/nils: coded
!
      integer :: i,j,k,cw,iobj
      real :: xval_p,xval_m,yval_p,yval_m,x2,y2,z2,minval,dist
      real :: zval_p,zval_m
!
      if (objects(iobj)%form == 'cylinder') then
        x2=objects(iobj)%r**2-(y(j)-objects(iobj)%x(2))**2
        y2=objects(iobj)%r**2-(x(i)-objects(iobj)%x(1))**2
      elseif (objects(iobj)%form == 'sphere') then
        x2=objects(iobj)%r**2&
            -(y(j)-objects(iobj)%x(2))**2&
            -(z(k)-objects(iobj)%x(3))**2
        y2=objects(iobj)%r**2&
            -(x(i)-objects(iobj)%x(1))**2&
            -(z(k)-objects(iobj)%x(3))**2
        z2=objects(iobj)%r**2&
            -(x(i)-objects(iobj)%x(1))**2&
            -(y(j)-objects(iobj)%x(2))**2        
        zval_p=objects(iobj)%x(3)+sqrt(z2)
        zval_m=objects(iobj)%x(3)-sqrt(z2)            
      endif
      xval_p=objects(iobj)%x(1)+sqrt(x2)
      xval_m=objects(iobj)%x(1)-sqrt(x2)            
      yval_p=objects(iobj)%x(2)+sqrt(y2)
      yval_m=objects(iobj)%x(2)-sqrt(y2)            
!
      minval=impossible
      cw=0
!
      dist=xval_p-x(i)
      if (dist<minval) then
        minval=dist
        cw=1
      endif
!
      dist=x(i)-xval_m
      if (dist<minval) then
        minval=dist
        cw=-1
      endif
!
      dist=yval_p-y(j)
      if (dist<minval) then
        minval=dist
        cw=2
      endif
!
      dist=y(j)-yval_m
      if (dist<minval) then
        minval=dist
        cw=-2
      endif
!
      if (objects(iobj)%form == 'sphere') then
        dist=zval_p-z(k)
        if (dist<minval) then
          minval=dist
          cw=3
        endif
!
        dist=z(k)-zval_m
        if (dist<minval) then
          minval=dist
          cw=-3
        endif
      endif
!
      call keep_compiler_quiet(k)
!
    endsubroutine find_closest_wall
!***********************************************************************  
    function ba_defined(i,j,k)
!
!  28-nov-2008/nils: coded
!
!  Check if ba for the point of interest has been defined for another direction.
!  This is only interesting if interpolation_method=='staircase',
!  otherwise this function always return .false.
!
      integer, intent(in) :: i,j,k
      logical :: lba1=.true.,lba2=.true.
      logical :: ba_defined
!
      if (interpolation_method=='staircase') then
        if (ba(i,j,k,1)==0 .or. ba(i,j,k,1)==9) then
          lba1=.false.
        endif
!
        if (ba(i,j,k,2)==0 .or. ba(i,j,k,2)==9) then
          lba2=.false.
        endif
!
        if (lba1 .or. lba2) then
          ba_defined=.true.
        else
          ba_defined=.false.
        endif
      else
        ba_defined=.false.
      endif
!
    endfunction ba_defined
!***********************************************************************  
    subroutine find_point(rij,rs,f,yin,xout,xmin,xmax,min,fout,x0)
!
!  20-mar-2009/nils: coded
!
!  Check if a grid line has any of it ends inside a solid cell - if so
!  find the point where the grid line enters the solid cell.
!
      integer, intent(in) :: min
      real, intent(in) :: xmin,xmax,rij,rs,f,yin,x0
      real, intent(out) :: fout,xout
      real :: xvar,xout0
!
      if (min == 1) then
        xvar=xmin
      else
        xvar=xmax
      endif
!
      if (rij > rs) then
        xout=xvar
        fout=f
      else
        xout0=sqrt(rs**2-yin**2)
        xout=xout0+x0
        if ((xout > xmax) .or. (xout < xmin)) then
          xout=x0-xout0
        endif
        fout=0
      endif
!
    endsubroutine find_point
!***********************************************************************  
    subroutine pencil_criteria_solid_cells()
!
!  All pencils that the Solid_Cells module depends on are specified here.
!
!  mar-2009/kragset: coded
!
      use Cdata
!
!  Request p and sij-pencils here
!  Request rho-pencil
      lpenc_requested(i_pp)=.true.
      lpenc_requested(i_sij)=.true.
      lpenc_requested(i_rho)=.true.
!
    endsubroutine pencil_criteria_solid_cells
!***********************************************************************  
    subroutine solid_cells_clean_up
!
!  Deallocate the variables allocated in solid_cells
!
!  7-oct-2010/dhruba: aped from hydro_kinematic
!
      print*, 'Deallocating some solid_cells variables ...'
      deallocate(fpnearestgrid)
      deallocate(c_dragx)
      deallocate(c_dragy)
      deallocate(c_dragz)
      print*, '..Done.'
!
    endsubroutine solid_cells_clean_up
!***********************************************************************
endmodule Solid_Cells
