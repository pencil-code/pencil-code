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
  logical :: lnointerception=.false.,lcheck_ba=.false.
  logical :: lclose_quad_rad_inter=.true.
  real                          :: rhosum,flow_dir,T0
  integer                       :: irhocount
  real                          :: theta_shift=1e-2
  real                          :: limit_close_linear=0.5
  real                          :: ineargridshift=1
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
       nspheres,sphere_radius,sphere_xpos,sphere_ypos,sphere_zpos,sphere_temp,&
       lclose_quad_rad_inter,ineargridshift
!
  namelist /solid_cells_run_pars/  &
       interpolation_method,object_skin,lclose_interpolation,lclose_linear,&
       limit_close_linear,lnointerception,nforcepoints,lcheck_ba,&
       lclose_quad_rad_inter,ineargridshift
!
!  Diagnostic variables (need to be consistent with reset list below).
!
  integer :: idiag_c_dragx=0       ! DIAG_DOC:
  integer :: idiag_c_dragy=0       ! DIAG_DOC:
  integer :: idiag_c_dragz=0       ! DIAG_DOC:
  integer :: idiag_c_dragx_p=0       ! DIAG_DOC:
  integer :: idiag_c_dragy_p=0       ! DIAG_DOC:
  integer :: idiag_c_dragz_p=0       ! DIAG_DOC:
  integer :: idiag_Nusselt=0            ! DIAG_DOC:
!
  integer, allocatable :: fpnearestgrid(:,:,:)
  real, allocatable    :: c_dragx(:), c_dragy(:), c_dragz(:), Nusselt(:)
  real, allocatable    :: c_dragx_p(:), c_dragy_p(:), c_dragz_p(:)
!
  contains
!***********************************************************************
    subroutine initialize_solid_cells(f)
!
!  Define the geometry of the solids.
!  There might be many separate solid objects of different geometries (currently
!  only cylinders and spheres are implemented however).
!
!  19-nov-2008/nils: coded
!  28-sep-2010/nils: added spheres
!  nov-2010/kragset: updated allocations related to drag calculations
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: icyl,isph,i
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
        if (nlong*(nlat+1)/=nforcepoints) then
          print*, "Warning: 2*nforcepoints should be square"
          print*,'nforcepoints=',nforcepoints
          print*,'nlong,nlat=',nlong,nlat
        endif
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
      call find_solid_cell_boundaries(f)
      call calculate_shift_matrix
!
!
! Find nearest grid point of the "forcepoints" on all cylinders. Arrays
! are only allocated if c_dragx, c_dragy or c_dragz is set in print.in.
! This must be changed if drag force is required for other purposes,
! e.g. if solid object is allowed to move.
!
      if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
          idiag_c_dragz /= 0 .or. idiag_Nusselt /= 0 .or. &
          idiag_c_dragx_p /= 0 .or. idiag_c_dragy_p /= 0 .or. &
          idiag_c_dragz_p /= 0 ) then
        allocate(fpnearestgrid(nobjects,nforcepoints,3))
        call fp_nearest_grid
        rhosum    = 0.0
        irhocount = 0
      endif
      if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
          idiag_c_dragz /= 0 .or. idiag_c_dragx_p /= 0 .or. &
          idiag_c_dragy_p /= 0 .or. idiag_c_dragz_p /= 0) then
        allocate(c_dragx(nobjects))
        allocate(c_dragy(nobjects))
        allocate(c_dragz(nobjects))
        allocate(c_dragx_p(nobjects))
        allocate(c_dragy_p(nobjects))
        allocate(c_dragz_p(nobjects))
      endif
      if (idiag_Nusselt /= 0) allocate(Nusselt(nobjects))
!
! Try to find flow direction
!
      flow_dir=0
      if (fbcx1(1) > 0) flow_dir= 1
      if (fbcx2(1) < 0) flow_dir=-1
      if (fbcy1(2) > 0) flow_dir= 2
      if (fbcy2(2) < 0) flow_dir=-2
      if (fbcz1(3) > 0) flow_dir= 3
      if (fbcz2(3) < 0) flow_dir=-3
      if (flow_dir > 0) then
        if (lroot) then
          print*,'By using fbc[x,y,z] I found the flow direction to be in the ',&
              flow_dir,' direction.'
        endif
      else
        do i=1,3
          if (lperi(i)) then
            if (.not. lperi(mod(i,3)+1) .and. .not. lperi(mod(i+1,3)+1)) then
              flow_dir=i
              if (lroot) then
              print*,'By using lperi I found the flow direction to be in the ',&
                  flow_dir,' direction.'
              endif
            endif
          endif
        enddo
        if (flow_dir == 0) then
          call fatal_error('initialize_solid_cells',&
              'I was not able to determine the flow direction!')
        endif
      endif
!
! Find inlet temperature
!
      if (ilnTT /= 0) then
        if (flow_dir== 1) T0=fbcx1(ilnTT)
        if (flow_dir==-1) T0=fbcx2(ilnTT)
        if (flow_dir== 2) T0=fbcy1(ilnTT)
        if (flow_dir==-2) T0=fbcy2(ilnTT)
        if (flow_dir== 3) T0=fbcz1(ilnTT)
        if (flow_dir==-3) T0=fbcz2(ilnTT)
        if (.not. ltemperature_nolog) T0=exp(T0)
      endif
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
          do icyl=1,nobjects
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
!  Assume a minimum radius for the forcepoints
      robj = robj+dxmin*ineargridshift
      else
        print*, "Warning: Subroutine fp_nearest_grid not implemented ", &
            "for this objectform."
      end if
!
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
          dist_to_cent2(1) = (x(ixl)-xobj)**2+(y(iyl)-yobj)**2+(z(izl)-zobj)**2
          dist_to_cent2(2) = (x(ixu)-xobj)**2+(y(iyl)-yobj)**2+(z(izl)-zobj)**2
          dist_to_cent2(3) = (x(ixu)-xobj)**2+(y(iyu)-yobj)**2+(z(izl)-zobj)**2
          dist_to_cent2(4) = (x(ixl)-xobj)**2+(y(iyu)-yobj)**2+(z(izl)-zobj)**2
          dist_to_cent2(5) = (x(ixl)-xobj)**2+(y(iyl)-yobj)**2+(z(izu)-zobj)**2
          dist_to_cent2(6) = (x(ixu)-xobj)**2+(y(iyl)-yobj)**2+(z(izu)-zobj)**2
          dist_to_cent2(7) = (x(ixu)-xobj)**2+(y(iyu)-yobj)**2+(z(izu)-zobj)**2
          dist_to_cent2(8) = (x(ixl)-xobj)**2+(y(iyu)-yobj)**2+(z(izu)-zobj)**2
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
!  that forcepoints are outside robj.
            if (dist_to_cent2(ipoint) .gt. robj**2 .and. inearest == 0) then
              inearest=ipoint
            else if (dist_to_cent2(ipoint) .gt. robj**2) then
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
!
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
    use Sub, only: dot2, dot
!
    real, dimension (mx,my,mz,mfarray), intent(in):: f
    real, dimension (mx,my,mz,mvar), intent(in)   :: df
    type (pencil_case), intent(in)                :: p
!
    real    :: fp_pressure, fp_gradT
    real    :: fp_stress(3,3)
    integer :: iobj, ifp, ix0, iy0, iz0, i, ilong, ilat
    real    :: longitude, latitude, dlong, dlat, robj, rforce
    real, dimension(nx) :: nu, twonu
    real    :: force_x, force_y, force_z, loc_Nus
    real    :: twopi, nvec(3), surfaceelement,surfacecoeff
    real    :: deltaT,Tobj, drag_norm, nusselt_norm
    character(len=10) :: objectform
!
    if (ldiagnos) then
      if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
          idiag_c_dragz /= 0 .or. idiag_Nusselt /= 0 .or. &
          idiag_c_dragx_p /= 0 .or. idiag_c_dragy_p /= 0 .or. &
          idiag_c_dragz_p /= 0) then
!
!  Reset cumulating quantities before calculations in first pencil
!
        if (imn == 1) then
          if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
              idiag_c_dragz /= 0) then
            c_dragx=0.
            c_dragy=0.
            c_dragz=0.
          endif
          if (idiag_Nusselt /= 0) Nusselt=0.
          if (idiag_c_dragx_p /= 0 .or. idiag_c_dragy_p /= 0 .or. &
              idiag_c_dragz_p /= 0) then
            c_dragx_p=0.
            c_dragy_p=0.
            c_dragz_p=0.
          endif
          rhosum=0
          irhocount=0
        endif
!
        call getnu(nu_pencil=nu,p=p)
        twopi=2.*pi
        twonu=2.*nu
!
        do iobj=1,nobjects
          robj = objects(iobj)%r
!  Integrating at radius rforce (for spheres)
          rforce = robj+dxmin*ineargridshift
          objectform = objects(iobj)%form
          if (objectform=='cylinder') then
            dlong = twopi/nforcepoints
            surfaceelement = dlong*rforce
            drag_norm=1./(2*robj)
          else if (objectform=='sphere') then
            dlong = twopi/nlong
            dlat  = pi/(nlat+1)
!  Surface term, normalized by the squared radius of the object.
!  Additional normalizing factors can be found in subroutine
!  dsolid_dt_integrate.
            surfacecoeff = dlong*dlat*rforce**2
            drag_norm=1./(pi*robj**2)
            nusselt_norm=1./(4*pi*robj**2)
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
!  Shifting the location of the forcpoints in the theta direction
!  in order to avoid problems with autotesting
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
                else
                  call fatal_error('dsolid_dt','No such objectform!')
                  call keep_compiler_quiet(nvec)
                end if
!
! Find force in x,y and z direction
!
                if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
                    idiag_c_dragz /= 0 .or. idiag_c_dragx_p /= 0 .or. &
                    idiag_c_dragy_p /= 0 .or. idiag_c_dragz_p /= 0) then
                  fp_pressure=p%pp(ix0-nghost)
                  fp_stress(:,:)=twonu(ix0-nghost)*p%rho(ix0-nghost)&
                      *p%sij(ix0-nghost,:,:)
!
!  Force in x-,y-, and z-directions
!
                  force_x = (-fp_pressure*nvec(1) &
                      + fp_stress(1,1)*nvec(1) &
                      + fp_stress(1,2)*nvec(2) &
                      + fp_stress(1,3)*nvec(3)) * surfaceelement
!
                  force_y = (-fp_pressure*nvec(2) &
                      + fp_stress(2,1)*nvec(1) &
                      + fp_stress(2,2)*nvec(2) &
                      + fp_stress(2,3)*nvec(3)) * surfaceelement
!
                  force_z = (-fp_pressure*nvec(3) &
                      + fp_stress(3,1)*nvec(1) &
                      + fp_stress(3,2)*nvec(2) &
                      + fp_stress(3,3)*nvec(3)) * surfaceelement
!
                end if
!
!  Local Nusselt number
!
                if (idiag_Nusselt /= 0) then
                  call dot(p%gTT(ix0-nghost,:),-nvec,fp_gradT)
                  Tobj=objects(iobj)%T
                  if (.not. ltemperature_nolog) Tobj=exp(Tobj)
                  deltaT=Tobj-T0
                  loc_Nus=fp_gradT*robj*2./deltaT * surfaceelement
                endif
!
                if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
                    idiag_c_dragz /= 0 .or. idiag_c_dragx_p /= 0 .or. &
                    idiag_c_dragy_p /= 0 .or. idiag_c_dragz_p /= 0) then
                  c_dragx(iobj) = c_dragx(iobj) + force_x * drag_norm
                  c_dragy(iobj) = c_dragy(iobj) + force_y * drag_norm
                  c_dragz(iobj) = c_dragz(iobj) + force_z * drag_norm
                  c_dragx_p(iobj) = c_dragx_p(iobj) + &
                      -fp_pressure*nvec(1) * drag_norm
                  c_dragy_p(iobj) = c_dragy_p(iobj) + &
                      -fp_pressure*nvec(2) * drag_norm
                  c_dragz_p(iobj) = c_dragz_p(iobj) + &
                      -fp_pressure*nvec(3) * drag_norm
                endif
                if (idiag_Nusselt /= 0) Nusselt(iobj) = Nusselt(iobj) &
                    + loc_Nus * nusselt_norm
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
!
    real    :: rhosum_all, c_dragx_all(nobjects), c_dragy_all(nobjects)
    real    :: c_dragz_all(nobjects), Nusselt_all(nobjects)
    real    :: c_dragx_p_all(nobjects),c_dragy_p_all(nobjects)
    real    :: c_dragz_p_all(nobjects)
    integer :: irhocount_all,iobj
    real    :: norm, refrho0
    character*100  :: numberstring
    character*500 :: solid_cell_drag
!
    if (ldiagnos) then
      if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 &
          .or. idiag_c_dragz /= 0 .or. idiag_Nusselt /= 0 &
          .or. idiag_c_dragx_p /= 0 .or. idiag_c_dragy_p /= 0 &
          .or. idiag_c_dragz_p /= 0) then
!
!  Collect and sum rhosum, irhocount, c_dragx, c_dragz, and c_dragy.
        call mpireduce_sum(rhosum,rhosum_all)
        call mpireduce_sum_int(irhocount,irhocount_all)
        if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
            idiag_c_dragz /= 0 .or. idiag_c_dragx_p /= 0 .or. &
            idiag_c_dragy_p /= 0 .or. idiag_c_dragz_p /= 0) then
          call mpireduce_sum(c_dragx,c_dragx_all,nobjects)
          call mpireduce_sum(c_dragy,c_dragy_all,nobjects)
          call mpireduce_sum(c_dragz,c_dragz_all,nobjects)
          call mpireduce_sum(c_dragx_p,c_dragx_p_all,nobjects)
          call mpireduce_sum(c_dragy_p,c_dragy_p_all,nobjects)
          call mpireduce_sum(c_dragz_p,c_dragz_p_all,nobjects)
        endif
        if (idiag_Nusselt /= 0) call mpireduce_sum(Nusselt,Nusselt_all,nobjects)
        if () then
        endif
!
        if (lroot) then
          refrho0 = rhosum_all / irhocount_all
!
!  Find drag and lift
!
          if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
              idiag_c_dragz /= 0 .or. idiag_c_dragx_p /= 0 .or. &
              idiag_c_dragy_p /= 0 .or. idiag_c_dragz_p /= 0) then
!  Normalizing factor. Additional factors was included in subroutine dsolid_dt.
            norm = 2. / (refrho0*init_uu**2)
            c_dragx = c_dragx_all * norm
            c_dragy = c_dragy_all * norm
            c_dragz = c_dragz_all * norm
            c_dragx_p = c_dragx_p_all * norm
            c_dragy_p = c_dragy_p_all * norm
            c_dragz_p = c_dragz_p_all * norm
!
!  Write drag coefficients for all objects
!  (may need to expand solid_cell_drag to more
!  characters if large number of objects).
!
            open(unit=81,file='data/dragcoeffs.dat',position='APPEND')
            write(solid_cell_drag,84) it-1, t
            do iobj=1,nobjects
              write(numberstring,82) c_dragx(iobj), c_dragy(iobj),c_dragz(iobj), &
                  c_dragx_p(iobj), c_dragy_p(iobj),c_dragz_p(iobj)
              call safe_character_append(solid_cell_drag,numberstring)
            enddo
            write(81,*) trim(solid_cell_drag)
            close(81)
84          format(1I8,1F15.8)
82          format(6F15.8)
          endif
!
!  Find Nusselt number
!
          if (idiag_Nusselt /= 0) then
            Nusselt = Nusselt_all
          endif
        endif
      endif
      if (idiag_c_dragx /= 0) fname(idiag_c_dragx)=c_dragx(1)
      if (idiag_c_dragy /= 0) fname(idiag_c_dragy)=c_dragy(1)
      if (idiag_c_dragz /= 0) fname(idiag_c_dragz)=c_dragz(1)
      if (idiag_c_dragx_p /= 0) fname(idiag_c_dragx_p)=c_dragx_p(1)
      if (idiag_c_dragy_p /= 0) fname(idiag_c_dragy_p)=c_dragy_p(1)
      if (idiag_c_dragz_p /= 0) fname(idiag_c_dragz_p)=c_dragz_p(1)
      if (idiag_Nusselt /= 0) fname(idiag_Nusselt)=Nusselt(1)
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
      idiag_c_dragx_p=0
      idiag_c_dragy_p=0
      idiag_c_dragz_p=0
      idiag_Nusselt=0
    endif
!
!  check for those quantities that we want to evaluate online
!
    do iname=1,nname
      call parse_name(iname,cname(iname),cform(iname),'c_dragx',idiag_c_dragx)
      call parse_name(iname,cname(iname),cform(iname),'c_dragy',idiag_c_dragy)
      call parse_name(iname,cname(iname),cform(iname),'c_dragz',idiag_c_dragz)
      call parse_name(iname,cname(iname),cform(iname),'c_dragx_p',idiag_c_dragx_p)
      call parse_name(iname,cname(iname),cform(iname),'c_dragy_p',idiag_c_dragy_p)
      call parse_name(iname,cname(iname),cform(iname),'c_dragz_p',idiag_c_dragz_p)
      call parse_name(iname,cname(iname),cform(iname),'Nusselt',idiag_Nusselt)
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
      real, dimension (mvar) :: f_tmp
      integer :: i,j,k,idir,xind,yind,zind,iobj

      real :: z_obj, y_obj, x_obj, r_obj, r_new, r_point, sin_theta, cos_theta
      real :: xmirror, ymirror, zmirror, dr
      integer :: lower_i, upper_i, lower_j, upper_j, ii, jj, kk
      integer :: lower_k, upper_k, ndims
      logical :: bax, bay, baz, lnew_interpolation_method
      real, dimension(3) :: xxp, phi
      character(len=10) :: form
!
!  Find ghost points based on the mirror interpolation method
!
      if (interpolation_method=='mirror') then
!
!  For fluid points very close to the solid surface the value of the point
!  is set from interpolation between the value at the closest grid line
!  and the value at the solid surface.
!
        do i=l1,l2
        do j=m1,m2
        do k=n1,n2
          if (lclose_linear) then
            if (ba(i,j,k,1)==10) then
              iobj=ba(i,j,k,4)
!
!  Check if we will use the old or the new interpolation method
!
              if (objects(iobj)%form=='sphere' .or. &
                  .not. lclose_quad_rad_inter) then
                lnew_interpolation_method=.true.
              else
                lnew_interpolation_method=.false.
              endif
!
              x_obj=objects(iobj)%x(1)
              y_obj=objects(iobj)%x(2)
              z_obj=objects(iobj)%x(3)
              r_obj=objects(iobj)%r
              if (objects(iobj)%form=='cylinder') then
                r_point=sqrt(((x(i)-x_obj)**2+(y(j)-y_obj)**2))
              elseif (objects(iobj)%form=='sphere') then
                r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2+(z(k)-z_obj)**2)
              endif
              dr=r_point-r_obj
              if ((dr > 0) .and. (dr<dxmin*limit_close_linear)) then
                xxp=(/x(i),y(j),z(k)/)
                if (lnew_interpolation_method) then
                  call close_interpolation(f,i,j,k,iobj,xxp,f_tmp,.true.,&
                      lnew_interpolation_method)
                  f(i,j,k,iux:iuz)=f_tmp(iux:iuz)
                  if (ilnTT > 0) f(i,j,k,ilnTT)=f_tmp(ilnTT)
                else
                  call close_interpolation(f,i,j,k,iobj,xxp,f_tmp,.true.,&
                      lnew_interpolation_method)
                  f(i,j,k,iux:iuz)=f_tmp(iux:iuz)
                  if (ilnTT > 0) then
                    f(i,j,k,ilnTT)=f_tmp(ilnTT)
                  endif
                endif
              endif
            endif
          endif
        enddo
        enddo
        enddo
!
!  Find ghost points based on the mirror interpolation method
!
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
              zmirror=z(k)
            elseif (form=='sphere') then
              r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2+(z(k)-z_obj)**2)
              r_new=r_obj+(r_obj-r_point)
              xmirror=(x(i)-x_obj)*r_new/r_point+x_obj
              ymirror=(y(j)-y_obj)*r_new/r_point+y_obj
              zmirror=(z(k)-z_obj)*r_new/r_point+z_obj
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
!  Check if we will use the old or the new interpolation method
!
             if (objects(iobj)%form=='sphere' .or. &
                 .not. lclose_quad_rad_inter) then
               lnew_interpolation_method=.true.
             else
               lnew_interpolation_method=.false.
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
            if (lnew_interpolation_method) then
              call interpolate_mirror_point_new(f,f_tmp,lower_i,upper_i,lower_j,&
                  upper_j,lower_k,upper_k,iobj,xmirror,ymirror,zmirror,ndims,&
                  lnew_interpolation_method)
              f(i,j,k,1:mvar)=f_tmp
            else
              call interpolate_mirror_point(f,phi,iux,lower_i,upper_i,lower_j,&
                  upper_j,lower_k,upper_k,iobj,xmirror,ymirror,zmirror,ndims,&
                  lnew_interpolation_method)
              f(i,j,k,iux)=-phi(1)
              call interpolate_mirror_point(f,phi,iuy,lower_i,upper_i,lower_j,&
                  upper_j,lower_k,upper_k,iobj,xmirror,ymirror,zmirror,ndims,&
                  lnew_interpolation_method)
              f(i,j,k,iuy)=-phi(1)
              call interpolate_mirror_point(f,phi,iuz,lower_i,upper_i,lower_j,&
                  upper_j,lower_k,upper_k,iobj,xmirror,ymirror,zmirror,ndims,&
                  lnew_interpolation_method)
              f(i,j,k,iuz)=-phi(1)
              if (ilnrho>0) then
                call interpolate_mirror_point(f,phi,ilnrho,lower_i,upper_i,&
                    lower_j,upper_j,lower_k,upper_k,iobj,xmirror,ymirror,&
                    zmirror,ndims,lnew_interpolation_method)
                f(i,j,k,ilnrho)=phi(1)
              endif
              if (ilnTT>0) then
                call interpolate_mirror_point(f,phi,ilnTT,lower_i,upper_i,&
                    lower_j,upper_j,lower_k,upper_k,iobj,xmirror,ymirror,&
                    zmirror,ndims,lnew_interpolation_method)
                f(i,j,k,ilnTT)=2*objects(iobj)%T-phi(1)
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
                call fatal_error('update_solid_cells','No such idir!')
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
    subroutine interpolate_mirror_point(f,phi_,ivar,lower_i,upper_i,lower_j,&
        upper_j,lower_k,upper_k,iobj,xmirror,ymirror,zmirror,ndims,&
        lnew_interpolation_method)
!
!  Interpolate value in a mirror point from the eight corner values
!
!  23-dec-2008/nils: coded
!  22-apr-2009/nils: added special treatment close to the solid surface
!
      use General, only: linear_interpolate
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mvar) :: f_tmp
      integer, intent(in) :: iobj,ndims
      integer, intent(in) :: lower_i,upper_i,lower_j,upper_j,ivar
      integer, intent(in) :: lower_k,upper_k
      real,    intent(in) :: xmirror,ymirror,zmirror
      real,dimension(3), intent(out):: phi_
      logical, intent(in) :: lnew_interpolation_method
!
      real, dimension(3) :: xxp
      real, dimension(3) :: gp
      integer, dimension(3) :: inear
!
      call keep_compiler_quiet(upper_i)
      call keep_compiler_quiet(upper_j)
      call keep_compiler_quiet(upper_k)
      call keep_compiler_quiet(ndims)
!
      xxp=(/xmirror,ymirror,zmirror/)
      inear=(/lower_i,lower_j,lower_k/)
      if (lnew_interpolation_method .and. ivar==iux) then
        call linear_interpolate(f,iux,iuz,xxp,gp,inear,.false.)
        phi_=gp
      else
        call linear_interpolate(f,ivar,ivar,xxp,gp(1),inear,.false.)
        phi_(1)=gp(1)
      endif
!
!  If the mirror point is very close to the surface of the object
!  some special treatment is required.
!
      if (lclose_interpolation .and. (ivar < 4 .or. ivar==ilnTT)) then
        if (lnew_interpolation_method .and. ivar==iux) then
          f_tmp(iux:iuz)=phi_
          call close_interpolation(f,lower_i,lower_j,lower_k,iobj,xxp,&
              f_tmp,.false.,lnew_interpolation_method)
          phi_=f_tmp(iux:iuz)
        else
          f_tmp(ivar)=phi_(1)
          call close_interpolation(f,lower_i,lower_j,lower_k,iobj,xxp,&
              f_tmp,.false.,lnew_interpolation_method)
          phi_(1)=f_tmp(ivar)
        endif
      endif
!
    endsubroutine interpolate_mirror_point
!***********************************************************************
    subroutine interpolate_mirror_point_new(f,f_tmp,lower_i,upper_i,lower_j,&
        upper_j,lower_k,upper_k,iobj,xmirror,ymirror,zmirror,ndims,&
        lnew_interpolation_method)
!
!  Interpolate value in a mirror point from the eight corner values
!
!  23-dec-2008/nils: coded
!  22-apr-2009/nils: added special treatment close to the solid surface
!
      use General, only: linear_interpolate
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mvar), intent(out) :: f_tmp
      integer, intent(in) :: iobj,ndims
      integer, intent(in) :: lower_i,upper_i,lower_j,upper_j
      integer, intent(in) :: lower_k,upper_k
      real,    intent(in) :: xmirror,ymirror,zmirror
      logical, intent(in) :: lnew_interpolation_method
!
      real, dimension(3) :: xxp
      integer, dimension(3) :: inear
      real :: rho_ghost_point,rho_fluid_point,T_fluid_point,T_ghost_point
!
      call keep_compiler_quiet(upper_i)
      call keep_compiler_quiet(upper_j)
      call keep_compiler_quiet(upper_k)
      call keep_compiler_quiet(ndims)
!
!  Linear interploation (use the eight corner points of the cube), same
! method as used in the particle module
!
      xxp=(/xmirror,ymirror,zmirror/)
      inear=(/lower_i,lower_j,lower_k/)
      call linear_interpolate(f,iux,mvar,xxp,f_tmp,inear,.false.)
!
!  If the mirror point is very close to the surface of the object
!  some special treatment is required.
!
      if (lclose_interpolation) then
        call close_interpolation(f,lower_i,lower_j,lower_k,iobj,xxp,&
            f_tmp,.false.,lnew_interpolation_method)
      endif
!
!  For the temperature boundaries being antisymmetric relative to the
!  object temperature are used
!
      if (ilnTT>0) then
        T_fluid_point=f_tmp(ilnTT)
        if (.not. ltemperature_nolog) T_fluid_point=exp(T_fluid_point)
        T_ghost_point=2*objects(iobj)%T-T_fluid_point
        if (.not. ltemperature_nolog) then
          f_tmp(ilnTT)=log(T_ghost_point)
        else
          f_tmp(ilnTT)=T_ghost_point
        endif
      endif
!
!  The gradient of the pressure should be zero at the fluid-solid intereface.
!  For isothermal cases this means that the density gradient is zero
!
      rho_fluid_point=f_tmp(ilnrho)
      if (.not. ldensity_nolog) rho_fluid_point=exp(rho_fluid_point)
      if (ilnTT>0) then
        rho_ghost_point=rho_fluid_point*T_fluid_point/T_ghost_point
      else
        rho_ghost_point=rho_fluid_point
      endif
      if (.not. ldensity_nolog) then
        f_tmp(ilnrho)=log(rho_ghost_point)
      else
        f_tmp(ilnrho)=rho_ghost_point
      endif
!
!  Antisymmetric boundaries are used for the velocity vector
!
      f_tmp(1:3)=-f_tmp(1:3)
!
    endsubroutine interpolate_mirror_point_new
!***********************************************************************
    subroutine close_interpolation(f,ix0_,iy0_,iz0_,iobj,xxp,f_tmp,&
        fluid_point,lnew_interpolation_method)
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
! bv      Border value              cornervalue(3*2)      -
!---------------------------------------------------------------------------
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      logical, intent(in) :: lnew_interpolation_method
      integer, intent(in) :: ix0_,iy0_,iz0_, iobj
      real, dimension(mvar), intent(inout) :: f_tmp
      real, dimension(3), intent(in) :: xxp
      logical, intent(in) :: fluid_point
!
      real, dimension(3) :: p_local,p_global,o_global, gpp
      real :: rs, rp
      real, dimension(2,2,2) :: rij
      real, dimension(3,2) :: cornervalue
      integer, dimension(3,2) :: cornerindex
      integer :: itop_bot,jtop_bot,ktop_bot
!
!  Check if we really want special treatment close to the fluid-solid
!  interface
!
      if ((.not. fluid_point .and. lclose_interpolation) &
          .or. ( fluid_point .and. lclose_linear)) then
!
!  Define some help variables
!
        o_global=objects(iobj)%x
        rs=objects(iobj)%r
        p_global=xxp
!
!  Find the corner points of the grid cell we are in based on one of the
!  corner points (given by ix0_,iy0_,iz0_)
!
        call find_corner_points(fluid_point,cornervalue,cornerindex,&
            ix0_,iy0_,iz0_,p_global,o_global)
!
!  Find the distance rij from all eight corner points to the object center
!
          do itop_bot=1,2
          do jtop_bot=1,2
          do ktop_bot=1,2
            rij(itop_bot,jtop_bot,ktop_bot)=sqrt(&
                (cornervalue(1,itop_bot)-o_global(1))**2+&
                (cornervalue(2,jtop_bot)-o_global(2))**2+&
                (cornervalue(3,ktop_bot)-o_global(3))**2)
          enddo
          enddo
          enddo
!
!  We want special treatment if at least one of the corner points are
!  inside the solid geometry, or if this is a fluid point.
!
        if ((minval(rij) < rs) .or. fluid_point) then
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
!  The new one is more general and works both for spheres and cylinders,
!  but it does not handle particles correctly
!  for the cases where one of the corner points of "gridplane" is inside
!  the solid geometry.
!
          if (lnew_interpolation_method) then
            call close_inter_new(f,f_tmp,p_local,p_global,o_global,rs,rp,&
                cornervalue,cornerindex, fluid_point,iobj)
          else
            call close_inter_old(f,gpp, rij, o_global, p_global, fluid_point,&
                iobj, cornervalue, &
                cornerindex,p_local, iux, rs, rp)
            f_tmp(iux)=gpp(1)
            call close_inter_old(f,gpp, rij, o_global, p_global, fluid_point,&
                iobj, cornervalue, &
                cornerindex,p_local, iuy, rs, rp)
            f_tmp(iuy)=gpp(1)
            call close_inter_old(f,gpp, rij, o_global, p_global, fluid_point,&
                iobj, cornervalue, &
                cornerindex,p_local, iuz, rs, rp)
            f_tmp(iuz)=gpp(1)
            if (ilnTT>0) then
              call close_inter_old(f,gpp, rij, o_global, p_global, fluid_point,&
                  iobj, cornervalue, &
                  cornerindex,p_local, ilnTT, rs, rp)
              f_tmp(ilnTT)=gpp(1)
            endif
          endif
        endif
      endif
!
    endsubroutine close_interpolation
!***********************************************************************
    subroutine close_inter_new(f,f_tmp,p_local,p_global,o_global,rs,rp,&
        cornervalue,cornerindex, fluid_point,iobj)
!
      use General, only: linear_interpolate
      use Sub
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mvar) :: fvar,f_tmp
      real, dimension(3) :: o_global, p_global, p_local, g_global
      integer, dimension(3) :: ngrids, inear
      real :: rp,rs, verylarge=1e9, rlmin, rl, r, r_pg, r_sg, r_sp, surf_val
      integer :: ndir, ndims, dir, vardir1,vardir2, constdir, topbot_tmp
      integer :: topbot, iobj
      real, dimension(3,2) :: cornervalue
      integer, dimension(3,2) :: cornerindex
      logical :: fluid_point
      integer :: lower_i, lower_j, lower_k, upper_i, upper_j, upper_k
      real :: xmirror, ymirror, zmirror, phi, theta
      real,  dimension(3) :: nr_hat, ntheta_hat, nphi_hat
      real :: vg_r, vg_phi, vg_theta
      real :: vp_r, vp_phi, vp_theta
!
      intent(out) :: f_tmp
      intent(in) ::  iobj
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
!  rl*p_local(dir)+o_global(dir)=cornervalue(dir)
!
            rl=(cornervalue(dir,topbot_tmp)-o_global(dir))/p_local(dir)
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
              g_global(constdir)=cornervalue(constdir,topbot)
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
        print*,'cornerindex=',cornerindex
        print*,'cornervalue=',cornervalue
        print*,'p_global=',p_global
        print*,'r,rs,rp=',r,rs,rp
        print*,'rs=',rs
        call fatal_error('close_interpolation',&
            'A valid radius is not found!')
      endif
!
!  Due to roundoff errors the value of g_global might end up outside the
!  grid cell - in that case put it back in where it belongs
!
      if (&
          (cornervalue(1,1)<=g_global(1).and.&
          cornervalue(1,2)>=g_global(1).or.nxgrid==1).and.&
          (cornervalue(2,1)<=g_global(2).and.&
          cornervalue(2,2)>=g_global(2).or.nygrid==1).and.&
          (cornervalue(3,1)<=g_global(3).and.&
          cornervalue(3,2)>=g_global(3).or.nzgrid==1)) then
        ! Everything okay
      else
        if (g_global(1)>cornervalue(1,2)) g_global(1)=cornervalue(1,2)
        if (g_global(1)<cornervalue(1,1)) g_global(1)=cornervalue(1,1)
        if (g_global(2)>cornervalue(2,2)) g_global(2)=cornervalue(2,2)
        if (g_global(2)<cornervalue(2,1)) g_global(2)=cornervalue(2,1)
        if (g_global(3)>cornervalue(3,2)) g_global(3)=cornervalue(3,2)
        if (g_global(3)<cornervalue(3,1)) g_global(3)=cornervalue(3,1)
      endif
!
!  Depending on the value of constdir the indeces of the corner points
!  specifying "gridplane" is given by lower_i, upper_i, lower_j ........
!  The physical position of "g" is stored in xmirror, ymirror and zmirror.
!
      if (constdir==1) then
        lower_i=cornerindex(constdir,topbot)
        upper_i=cornerindex(constdir,topbot)
        lower_j=cornerindex(vardir1,1)
        upper_j=cornerindex(vardir1,2)
        lower_k=cornerindex(vardir2,1)
        upper_k=cornerindex(vardir2,2)
        xmirror=g_global(constdir)
        ymirror=g_global(vardir1)
        zmirror=g_global(vardir2)
      elseif (constdir==2) then
        lower_j=cornerindex(constdir,topbot)
        upper_j=cornerindex(constdir,topbot)
        lower_k=cornerindex(vardir1,1)
        upper_k=cornerindex(vardir1,2)
        lower_i=cornerindex(vardir2,1)
        upper_i=cornerindex(vardir2,2)
        ymirror=g_global(constdir)
        xmirror=g_global(vardir2)
        zmirror=g_global(vardir1)
      elseif (constdir==3) then
        lower_k=cornerindex(constdir,topbot)
        upper_k=cornerindex(constdir,topbot)
        lower_i=cornerindex(vardir1,1)
        upper_i=cornerindex(vardir1,2)
        lower_j=cornerindex(vardir2,1)
        upper_j=cornerindex(vardir2,2)
        zmirror=g_global(constdir)
        xmirror=g_global(vardir1)
        ymirror=g_global(vardir2)
      endif
!
!  Let "fvar" be the physical value of a give variable on "gridplane".
!  The value of "fvar" is then found by interpolation between the four corner
!  points of "gridplane".
!
      inear=(/lower_i,lower_j,lower_k/)
      call linear_interpolate(f,1,mvar,g_global,fvar,inear,.false.)
!
!  Now we know the value associated with any variable in the point "g",
!  given by "fvar".
!  Furthermore we know the value associated with any variable in point "s"
!  on the object surface to be zero for any of the velocities and equal to
!  the solid temperature for the temperature. This value is given by "surf_val".
!  By interpolation it is now straight forward to find the value also in "p".
!  First find the distance, "r_pg", from "p" to "g" to be used as weight,
!  then normalize by the full distance, "r_sg", between "s" and "g".
!
      r_pg=r-rp
      r_sg=r-rs
      r_sp=rp-rs
!
!  Find temperature
!
      if (ilnTT>0) then
        if (.not. ltemperature_nolog) then
          call fatal_error('close_inter_new',&
              'Due to interpolation it is not correc to use lnTT!')
        endif
        surf_val=objects(iobj)%T
        f_tmp(ilnTT)=(fvar(ilnTT)*r_sp+surf_val*r_pg)/r_sg
      endif
!
!  Find velocity.
!  If lclose_quad_rad_inter = true we must find the
!  velocities in the r, theta and phi
!  directions at the point "g". This will then be used to set up a
!  linear interpolation for v_theta and v_phi and a quadratic interpolation
!  for v_r.
!
      surf_val=0
      if (lclose_quad_rad_inter) then
!
!  The unity vector "nr_hat" is normal to the solid surface, while
!  "nphi_hat" and "ntheta_hat" are the unit vectors in the two angular
!  directions. The angle "theta" is zero in the positive x-direction,
!  while "phi" is zero in the positive z-direction.
!
        phi=acos(p_local(3)/rp)
        theta=atan(p_local(2)/p_local(1))
        if (p_local(2) < 0) then
          if (theta > 0) then
            theta=theta+pi
          else
            theta=theta+2*pi
          endif
        else
          if (theta<0) then
            theta=theta+pi
          endif
        endif
!
        nr_hat    =(/cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)/)
        nphi_hat  =(/-cos(phi)*cos(theta),-cos(phi)*sin(theta),sin(phi)/)
        ntheta_hat=(/-sin(theta),cos(theta),0./)
!
!  Having found the unit vectors in the r, theta and phi directions we can now
!  find the velocities in the same three directions at point "g".
!
        call dot(nr_hat    ,fvar,vg_r)
        call dot(nphi_hat  ,fvar,vg_phi)
        call dot(ntheta_hat,fvar,vg_theta)
!
!  Now it is time to use linear and quadratic interpolation to find the
!  velocities in point "p".
!
        vp_phi  =(vg_phi  *r_sp+surf_val*r_pg)/r_sg
        vp_theta=(vg_theta*r_sp+surf_val*r_pg)/r_sg
        vp_r    =(vg_r    *(r_sp/r_sg)**2)
!
!  Finally the velocities found in the spherical coordinate system can
!  now be transfered back to the cartesian coordinate system.
!
        f_tmp(iux:iuz)=vp_r*nr_hat+vp_theta*ntheta_hat+vp_phi*nphi_hat
      else
        f_tmp(iux:iuz)=(fvar(iux:iuz)*r_sp+surf_val*r_pg)/r_sg
      endif
!
    end subroutine close_inter_new
!***********************************************************************
    subroutine find_corner_points(fluid_point,cornervalue,cornerindex,&
        ix0_,iy0_,iz0_,p_global,o_global)
!
!  8-dec-10: coded (nils)
!
!  Based on one of the corner points this routine find all corner points
!  of the fluid cell inwhich we are.
!  Furthermore; if we are at a fluid point p_global is shifted slightly
!  inside the domain.
!
      logical, intent(in) :: fluid_point
      integer, intent(in) :: ix0_,iy0_,iz0_
      real, dimension(3), intent(inout) :: p_global
      real, dimension(3), intent(in) :: o_global
      real, dimension(3,2), intent(out) :: cornervalue
      integer, dimension(3,2), intent(out) :: cornerindex
      real :: smallx
      integer :: ix0,iy0,iz0,ix1,iy1,iz1
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
!  Put help variables into arrays
!
          cornervalue(1,1)=x(ix0)
          cornervalue(2,1)=y(iy0)
          cornervalue(3,1)=z(iz0)
          cornervalue(1,2)=x(ix1)
          cornervalue(2,2)=y(iy1)
          cornervalue(3,2)=z(iz1)
          cornerindex(1,1)=ix0
          cornerindex(2,1)=iy0
          cornerindex(3,1)=iz0
          cornerindex(1,2)=ix1
          cornerindex(2,2)=iy1
          cornerindex(3,2)=iz1
!
      end subroutine find_corner_points
!***********************************************************************
    subroutine close_inter_old(f,gpp, rij, o_global, p_global, fluid_point,&
        iobj, cornervalue, cornerindex,&
        p_local, ivar1, rs, rp)
!
!  7-dec-2010/nils: moved from close_interpolation
!
!  This is the old version of the close interpolation which only work
!  for cylindrical geometries.
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension(2,2,2) :: rij
      integer, dimension(6) :: constdir_arr, vardir_arr, topbot_arr
      integer :: vardir, constdir
      real, dimension(3) :: xyint, gpp, o_global, p_global, p_local
      logical :: fluid_point
      integer :: iobj, maxcounter, counter, topbot_tmp
      real :: R1, verylarge, Rsmall, xtemp, r, rp, rs, rij_min
      real :: rij_max, inputvalue, varval, x1, f1, f1y, f1x, x2, f2, f2y, f2x
      real :: rint1, rint2, fintx, fint, finty, rps, fint_ur, rintp
      real, save :: urp, utp
      real :: surf_val, fint_ut, drp, dri
      real, dimension(3,2) :: cornervalue
      integer, dimension(3,2) :: cornerindex
      integer :: dirconst, dirvar, topbot, min, ivar1
!
      if (objects(iobj)%form=='cylinder') then
        constdir_arr=(/2,2,1,1,0,0/)
        vardir_arr  =(/1,1,2,2,0,0/)
        topbot_arr  =(/2,1,2,1,0,0/)
      elseif (objects(iobj)%form=='sphere') then
        constdir_arr=(/3,3,2,2,1,1/)
        vardir_arr  =(/1,1,1,1,3,3/)
        topbot_arr  =(/2,1,2,1,2,1/)
      endif
      verylarge=1e9
      R1=verylarge
      Rsmall=verylarge/2.0
!
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
            *(cornervalue(constdir,topbot_tmp)-objects(iobj)%x(constdir))
!
!  Find the distance, r, from the center of the object
!  to the point where the normal cross the grid line
!
        if (abs(xtemp) > verylarge) then
          r=verylarge*2
        else
          r=sqrt(xtemp**2+(cornervalue(constdir,topbot_tmp)&
              -objects(iobj)%x(constdir))**2)
        endif
!
!  Check if the point xtemp is outside the object,
!  outside the point "p" and that it cross the grid
!  line within this grid cell
!
        if ((r > rs) .and. (r > rp) &
            .and.(xtemp+objects(iobj)%x(vardir)>=cornervalue(vardir,1))&
            .and.(xtemp+objects(iobj)%x(vardir)<=cornervalue(vardir,2)))then
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
          xyint(constdir)=cornervalue(constdir,topbot_tmp)
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
          inputvalue=cornervalue(constdir,topbot_tmp)&
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
        print*,'cornerindex=',cornerindex
        print*,'cornervalue=',cornervalue
        print*,'p_global=',p_global
        print*,'xtemp=',xtemp
        print*,'r,rs,rp=',r,rs,rp
        print*,'R1,Rsmall=',R1,Rsmall
        print*,'rij,rs=',rij,rs
        print*,'dirvar,dirconst,topbot=',dirvar,dirconst,topbot
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
        varval=f(cornerindex(dirvar,1),cornerindex(dirconst,topbot),&
            cornerindex(3,1),ivar1)
      else
        varval=f(cornerindex(dirconst,topbot),cornerindex(dirvar,1),&
            cornerindex(3,1),ivar1)
      endif
      call find_point(rij_min,rs,varval,inputvalue,x1,&
          cornervalue(dirvar,1),cornervalue(dirvar,2),&
          min,f1,objects(iobj)%x(dirvar))
!
! If we want quadratic interpolation of the radial velocity we
! must find both the interploated x and y velocity in order to
! do interpolations for the radial and theta directions.
!
      if (lclose_quad_rad_inter .and. ivar1==iux) then
        if (dirconst == 2) then
          varval=f(cornerindex(dirvar,1),cornerindex(dirconst,topbot),&
              cornerindex(3,1),iuy)
        else
          varval=f(cornerindex(dirconst,topbot),cornerindex(dirvar,1),&
              cornerindex(3,1),iuy)
        endif
        call find_point(rij_min,rs,varval,inputvalue,x1,&
            cornervalue(dirvar,1),cornervalue(dirvar,2),&
            min,f1y,objects(iobj)%x(dirvar))
        f1x=f1
      endif
!
      min=0
      if (dirconst == 2) then
        varval=f(cornerindex(dirvar,2),cornerindex(dirconst,topbot),&
            cornerindex(3,1),ivar1)
      else
        varval=f(cornerindex(dirconst,topbot),cornerindex(dirvar,2),&
            cornerindex(3,1),ivar1)
      endif
      call find_point(rij_max,rs,varval,inputvalue,x2,&
          cornervalue(dirvar,1),cornervalue(dirvar,2),&
          min,f2,objects(iobj)%x(dirvar))
!
! If we want quadratic interpolation of the radial velocity we
! must find both the interploated x and y velocity in order to
! do interpolations for the radial and theta directions.
!
      if (lclose_quad_rad_inter .and. ivar1==iux) then
        if (dirconst == 2) then
          varval=f(cornerindex(dirvar,2),cornerindex(dirconst,topbot),&
              cornerindex(3,1),iuy)
        else
          varval=f(cornerindex(dirconst,topbot),cornerindex(dirvar,2),&
              cornerindex(3,1),iuy)
        endif
        call find_point(rij_max,rs,varval,inputvalue,x2,&
            cornervalue(dirvar,1),cornervalue(dirvar,2),&
            min,f2y,objects(iobj)%x(dirvar))
        f2x=f2
      else
        ! To keep the compiler quiet
        f1x=0
        f2x=0
      endif
!
!  Find the interpolation values between the two endpoints of
!  the line and the normal from the objects.
!
      rint1=xyint(dirvar)-x1
      rint2=x2-xyint(dirvar)
!
      if (lclose_quad_rad_inter .and. (ivar1 /= iuz) &
          .and. (ivar1 /= ilnTT)) then
        if (ivar1==iux) then
          fintx=(rint1*f2x+rint2*f1x)/(x2-x1)
          finty=(rint1*f2y+rint2*f1y)/(x2-x1)
          fint_ur    =fintx*p_local(1)/rs+finty*p_local(2)/rs
          fint_ut=finty*p_local(1)/rs-fintx*p_local(2)/rs
          drp=rp-rs
          dri=Rsmall-rs
          urp=(drp/dri)**2*fint_ur
          utp=(drp/dri)*fint_ut
          gpp(1)=urp*p_local(1)/rs-utp*p_local(2)/rs
        elseif (ivar1==iuy) then
          gpp(1)=urp*p_local(2)/rs+utp*p_local(1)/rs
        else
          call fatal_error('close_interpolation',&
              'Your ivar1 is not correct!')
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
        surf_val=0
        if (ivar1==ilnTT) surf_val=objects(iobj)%T
        gpp(1)=(rps*fint+rintp*surf_val)/(Rsmall-rs)
      endif
!
    end subroutine close_inter_old
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
!
      if (present(iostat)) then
        read(unit,NML=solid_cells_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=solid_cells_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_solid_cells_init_pars
!***********************************************************************
    subroutine read_solid_cells_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=solid_cells_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=solid_cells_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_solid_cells_run_pars
!***********************************************************************
    subroutine write_solid_cells_init_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=solid_cells_init_pars)
!
    endsubroutine write_solid_cells_init_pars
!***********************************************************************
    subroutine write_solid_cells_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=solid_cells_run_pars)
!
    endsubroutine write_solid_cells_run_pars
!***********************************************************************
    subroutine find_solid_cell_boundaries(f)
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
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
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
                    ba(i,j,k,3)=-1
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==1) ba(i,j,k,3)=-1
                  endif
                endif
                !
                if (z(k+2)>zval_p .and. z(k+1)<zval_p) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,3)=-2
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==1) ba(i,j,k,3)=-2
                  endif
                endif
                !
                if (z(k+3)>zval_p .and. z(k+2)<zval_p) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,3)=-3
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==1) ba(i,j,k,3)=-3
                  endif
                endif
                !
                if (z(k-1)<zval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,3)=1
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-1) ba(i,j,k,3)=1
                  endif
                endif
                !
                if (z(k-2)<zval_m .and. z(k-1)>zval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,3)=2
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-1) ba(i,j,k,3)=2
                  endif
                endif
                !
                if (z(k-3)<zval_m .and. z(k-2)>zval_m) then
                  if (.not. ba_defined(i,j,k)) then
                    ba(i,j,k,3)=3
                    ba(i,j,k,4)=iobj
                  else
                    call find_closest_wall(i,j,k,iobj,cw)
                    if (cw==-1) ba(i,j,k,3)=3
                  endif
                endif
                !
                if (ba(i,j,k,3)==0) then
                  ba(i,j,k,3)=9
                  ba(i,j,k,4)=iobj
                endif
                !
              endif
            enddo
          endif
        enddo
        enddo
        else
!
!  If the object is a cylinder then every point inside the cylinder will
!  be infinetly far from the surface in the z-direction.
!
          do i=l1,l2
            do j=m1,m2
              do k=n1,n2
                if ((ba(i,j,k,1).ne.0) .and. (ba(i,j,k,1).ne.10)) then
                  ba(i,j,k,3)=9
                endif
              enddo
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
              r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2)
            elseif (form=='sphere') then
              r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2&
                  +(z(mz-nghost+k)-z_obj)**2)
            else
              call fatal_error('find_solid_cell_boundaries','No such form!')
            endif
            if (r_point < r_obj) then
              ba(i,j,mz-nghost+k,1:3)=11
              ba(i,j,mz-nghost+k,4)=iobj
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
!
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
!  Set zero value of all variables inside the solid geometry far from
!  all interfaces. This is done for more easy interpretation of postprocessing.
!
      if (it==1) then
        do iobj=1,nobjects
          do i=1,mx
          do j=1,my
          do k=1,mz
            if (ba(i,j,k,1)==9 .and. ba(i,j,k,2)==9 .and. ba(i,j,k,3)==9) then
              f(i,j,k,iux:iuz)=0
            endif
          enddo
          enddo
          enddo
        enddo
      endif
!
!  Check that a fluid point is really outside a solid geometry
!
      if (lcheck_ba) then
        do iobj=1,nobjects
          x_obj=objects(iobj)%x(1)
          y_obj=objects(iobj)%x(2)
          z_obj=objects(iobj)%x(3)
          r_obj=objects(iobj)%r
          form=objects(iobj)%form
          do i=1,mx
            do j=1,my
              do k=1,mz
                if (form=='cylinder') then
                  r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2)
                elseif (form=='sphere') then
                  r_point=sqrt((x(i)-x_obj)**2+(y(j)-y_obj)**2+(z(k)-z_obj)**2)
                else
                  call fatal_error('find_solid_cell_boundaries','No such form!')
                endif
                if (r_point > r_obj) then
                  if ((ba(i,j,k,1) .ne. 0 ) .and. (ba(i,j,k,1) .ne. 10)) then
                    print*,'i,j,k=',i,j,k
                    print*,'ba(i,j,k,1)=',ba(i,j,k,1)
                    print*,'r_point,r_obj=',r_point,r_obj
                    print*,'x(i),y(j),z(k)=',x(i),y(j),z(k)
                    call fatal_error('find_solid_cell_boundaries',&
                        'Point marked as fluid point but seems not to be...')
                  endif
                else
                  if ((ba(i,j,k,1)==0).or.(ba(i,j,k,1)==10).or.&
                      (ba(i,j,k,2)==0).or.(ba(i,j,k,2)==10).or.&
                      (ba(i,j,k,3)==0).or.(ba(i,j,k,3)==10))then
                    print*,'i,j,k=',i,j,k
                    print*,'ba(i,j,k,1)=',ba(i,j,k,1)
                    print*,'ba(i,j,k,2)=',ba(i,j,k,2)
                    print*,'ba(i,j,k,3)=',ba(i,j,k,3)
                    print*,'r_point,r_obj=',r_point,r_obj
                    print*,'x(i),y(j),z(k)=',x(i),y(j),z(k)
                    call fatal_error('find_solid_cell_boundaries',&
                        'Point marked as a solid point but seems not to be...')
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
      endif
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
          if (ba(i,j,k,idir)/=0 .and. ba(i,j,k,idir)/=9.) then
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
      if (idiag_Nusselt /= 0) lpenc_requested(i_gTT)=.true.
      if (idiag_Nusselt /= 0) lpenc_requested(i_tcond)=.true.
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
      deallocate(c_dragx_p)
      deallocate(c_dragy_p)
      deallocate(c_dragz_p)
      deallocate(Nusselt)
      print*, '..Done.'
!
    endsubroutine solid_cells_clean_up
!***********************************************************************
endmodule Solid_Cells
