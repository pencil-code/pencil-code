! $Id: solid_cells.f90 20898 2013-08-19 11:05:32Z nils.e.haugen@gmail.com $
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
  use General, only: keep_compiler_quiet
  use Messages
  use Chemistry
  use EquationOfState
!
  implicit none
!
  include 'solid_cells.h'
!
  integer, parameter            :: max_items=5
  integer                       :: ncylinders=0,nspheres=0,dummy
  integer                       :: nobjects,nlong,nlat,nforcepoints=200,nsvar=0
  integer                       :: ixs=0,iys=0,izs=0,ius=0,ivs=0,iws=0,iRs=0,iTs=0
  integer                       :: irhocount, osci_dir=1
  real, dimension(max_items)    :: cylinder_radius, sphere_radius
  real, dimension(max_items)    :: cylinder_temp=330.0, sphere_temp=330.0
  real, dimension(max_items)    :: cylinder_xpos, cylinder_ypos
  real, dimension(max_items)    :: cylinder_xvel=0.0, cylinder_yvel=0.0
  real, dimension(max_items)    :: cylinder_theta=0.0
  real, dimension(max_items)    :: sphere_xpos, sphere_ypos, sphere_zpos
  real, dimension(max_items)    :: sphere_xvel=0.0,sphere_yvel=0.0,sphere_zvel=0.0
  real, dimension(max_items)    :: sphere_theta=0.0,sphere_phi=0.0
  real, dimension(max_items,8)  :: dfs,fs
  integer, dimension(mx,my,mz,4):: ba,ba_shift
  character (len=labellen), dimension(ninit) :: initsolid_cells='nothing'
  character (len=10), dimension(8) :: svarname
  logical :: lnointerception=.false.,lcheck_ba=.false.
  logical :: lerror_norm=.false.,lNusselt_output=.false.,locdensity_error=.false.
  logical :: locchemspec_error=.false.,loutput_local_reaction_rate=.false.
  logical :: lset_flow_dir=.false.,lpos_advance=.false.,lradius_advance=.false.
  logical :: lsecondorder_rho=.true., lsecondorder_chem=.true.
  logical :: lstefan_flow=.true.
  real    :: rhosum,flow_dir,T0,flow_dir_set,theta_shift=0.0
  real    :: skin_depth=0,init_uu=0,ampl_noise=0,object_skin=0
  real    :: limit_close_linear=0.7, ineargridshift=1.0
  real    :: osci_A=0.5,osci_f=4.0,osci_t, solid_reactions_intro_time=0.0
!
! For chemestry
!
  integer :: isO2, ichemsO2, isCO2, ichemsCO2, isCO, ichemsCO, isN2, ichemsN2
  logical :: lsO2, lsCO2, lsN2, lsCO
  real    :: Rgas_unit_sys, pressure0=10.13e5, solid_ds=0.0, solid_dt=0.0
  real    :: srho=1.08   ![g/cm^3]
  real    :: scp=1.465e7    ![erg/g/K]
  real, dimension(nchemspec) :: Mspecies
  logical :: new_Stefan=.false.
  logical :: lfull_diff_ib=.false.
  logical :: lsimple_diff_ib=.true.
!
  type solid_object
    character(len=10) :: form
    real :: r,T
    real, dimension(3) :: x, vel
    real, dimension(2) :: rot
  end type solid_object
!
  type(solid_object), dimension(max_items) :: objects
!
  namelist /solid_cells_init_pars/ &
       ncylinders,cylinder_radius,cylinder_xpos,cylinder_ypos,cylinder_xvel, &
       cylinder_yvel,cylinder_temp,cylinder_theta, &
       nspheres,sphere_radius,sphere_xpos,sphere_ypos,sphere_zpos,sphere_xvel, &
       sphere_yvel,sphere_zvel,sphere_theta,sphere_phi,sphere_temp, &
       ineargridshift,lset_flow_dir,flow_dir_set,nforcepoints,object_skin, &
       initsolid_cells,skin_depth,init_uu,ampl_noise,limit_close_linear, &
       pressure0, solid_reactions_intro_time,new_Stefan
!
  namelist /solid_cells_run_pars/  &
       lnointerception,lcheck_ba,lerror_norm,lpos_advance,lradius_advance, &
       lNusselt_output,osci_A,osci_f,osci_dir,osci_t, lsecondorder_rho, &
       lsecondorder_chem, solid_reactions_intro_time, lstefan_flow, locdensity_error, &
       locchemspec_error,loutput_local_reaction_rate,new_Stefan
!
!  Diagnostic variables (need to be consistent with reset list below).
!
  integer :: idiag_c_dragx=0       ! DIAG_DOC:
  integer :: idiag_c_dragy=0       ! DIAG_DOC:
  integer :: idiag_c_dragz=0       ! DIAG_DOC:
  integer :: idiag_c_dragx_p=0     ! DIAG_DOC:
  integer :: idiag_c_dragy_p=0     ! DIAG_DOC:
  integer :: idiag_c_dragz_p=0     ! DIAG_DOC:
  integer :: idiag_Nusselt=0       ! DIAG_DOC:
!
  integer, allocatable :: fpnearestgrid(:,:,:)
  real, allocatable    :: c_dragx(:), c_dragy(:), c_dragz(:), Nusselt(:)
  real, allocatable    :: c_dragx_p(:), c_dragy_p(:), c_dragz_p(:)
  real, allocatable    :: vs_normal(:), heat_cond(:), char_consumption(:)
!
  contains
!******************************************************************************************
    subroutine register_solid_cells()
!
!  Set up indices for access to the fs and dfs arrays
!
    integer :: isvar, kchem
!
!    if (lroot) call svn_id("$Id: solid_cells.f90 zhuangzhenya@126.com $")
!
!  Indices for solid_cells position.
!
    if (lpos_advance) then
      ixs=nsvar+1
      svarname(nsvar+1)='ixs'
      iys=nsvar+2
      svarname(nsvar+2)='iys'
      izs=nsvar+3
      svarname(nsvar+3)='izs'
!  Indices for solid velocity.
      ius=nsvar+4
      svarname(nsvar+4)='ius'
      ivs=nsvar+5
      svarname(nsvar+5)='ivs'
      iws=nsvar+6
      svarname(nsvar+6)='iws'
!  Increase nsvar accordingly.
      nsvar=nsvar+6
    endif
!
!  Increase nsvar for radius when needed
!
    if (lradius_advance) then
      iRs=nsvar+1
      svarname(nsvar+1)='iRs'
      nsvar=nsvar+1
    endif
!
!  Increase nsvar for temperature when needed.
!
    iTs=nsvar+1
    svarname(nsvar+1)='iTs'
    nsvar=nsvar+1
!
    if (lroot .and. nsvar>0) then
      open(3,file=trim(datadir)//'/svarname.dat',status='replace')
      do isvar=1,nsvar
        write(3,"(i4,2x,a)") isvar, svarname(isvar)
      enddo
      close(3)
    endif
!
!  Set variables for solid chemistry.
!
    if (lchemistry) then
      do kchem=1, nchemspec
        Mspecies(kchem)=species_constants(kchem,1)
      enddo
!
      call find_species_index('O2', isO2, ichemsO2, lsO2 )
      call find_species_index('CO2',isCO2,ichemsCO2,lsCO2)
      call find_species_index('CO', isCO, ichemsCO, lsCO )
      call find_species_index('N2', isN2, ichemsN2, lsN2 )
    else
      call fatal_error('register_solid_cells','This is for chemistry only')
    endif
!
    end subroutine register_solid_cells
!*************************************************************************
    subroutine initialize_solid_cells(f)
!
!  Define the geometry of the solids.
!  Currently only cylinders and spheres are implemented.
!
      use General, only: safe_character_assign
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: icyl, isph, i, iobj, j, k, iv
      logical :: file_exists, dfile_exists
      character(len=fnlen) :: directory_snap='/proc0'
      character(len=fnlen) :: file_directory
!
!  Loop over all cylinders
!
      do icyl=1,ncylinders
        if (cylinder_radius(icyl)>0) then
          objects(icyl)%r    = cylinder_radius(icyl)
          objects(icyl)%x(1) = cylinder_xpos(icyl)
          objects(icyl)%x(2) = cylinder_ypos(icyl)
          objects(icyl)%x(3) = 0.0
          objects(icyl)%T    = cylinder_temp(icyl)
          objects(icyl)%form = 'cylinder'
          objects(icyl)%rot(1)=cylinder_theta(icyl)
          objects(icyl)%rot(2)=0.0
          objects(icyl)%vel(1)=cylinder_xvel(icyl)
          objects(icyl)%vel(2)=cylinder_yvel(icyl)
          objects(icyl)%vel(3)=0.0
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
          objects(isph+ncylinders)%r    = sphere_radius(isph)
          objects(isph+ncylinders)%x(1) = sphere_xpos(isph)
          objects(isph+ncylinders)%x(2) = sphere_ypos(isph)
          objects(isph+ncylinders)%x(3) = sphere_zpos(isph)
          objects(isph+ncylinders)%T    = sphere_temp(isph)
          objects(isph+ncylinders)%form = 'sphere'
          objects(isph+ncylinders)%rot(1)=sphere_theta(isph)
          objects(isph+ncylinders)%rot(2)=sphere_phi(isph)
          objects(isph+ncylinders)%vel(1)=sphere_xvel(icyl)
          objects(isph+ncylinders)%vel(2)=sphere_yvel(icyl)
          objects(isph+ncylinders)%vel(3)=sphere_zvel(icyl)
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
        sphere_xvel(1)=impossible
        sphere_yvel(1)=impossible
        sphere_zvel(1)=impossible
        sphere_theta(1)=impossible
        sphere_phi(1)=impossible
      else
!
! Find number of lines of longitude and latitude such that
! nforcepoints=nlong*(nlat+1) and nlat=nlong/2-1
!
        nlong=int(sqrt(2.0*nforcepoints))
        nlat=int(0.5*nlong)-1
        if (nlong*(nlat+1)/=nforcepoints) then
          print*, "Warning: 2*nforcepoints should be square"
          print*,'nforcepoints=',nforcepoints
          print*,'nlong,nlat=',nlong,nlat
        endif
      endif
      if (ncylinders==0) then
        cylinder_radius(1)=impossible
        cylinder_xpos(1)=impossible
        cylinder_ypos(1)=impossible
        cylinder_temp(1)=impossible
        cylinder_xvel(1)=impossible
        cylinder_yvel(1)=impossible
        cylinder_theta(1)=impossible
      endif
!
! Search if there is a snapshot for solid objects.
!
      inquire(file=trim(datadir)//trim(directory_snap)//'/svar.dat',exist=file_exists)
      inquire(file=trim(datadir)//trim(directory_snap)//'/dsvar.dat',exist=dfile_exists)
      if (lroot .and. file_exists .and. dfile_exists) &
        print*, 'Find data/svar.dat and data/dsvar.dat'
      if (file_exists .and. dfile_exists) then
        call safe_character_assign(file_directory,trim(datadir)//trim(directory_snap))
        call read_snapshot_solid_cells(file_directory)
      else
        if (lpos_advance) then
          do iobj=1, nobjects
            fs(iobj,ixs)=objects(iobj)%x(1)
            fs(iobj,iys)=objects(iobj)%x(2)
            fs(iobj,izs)=objects(iobj)%x(3)
            fs(iobj,ius)=objects(iobj)%vel(1)
            fs(iobj,ivs)=objects(iobj)%vel(2)
            fs(iobj,iws)=objects(iobj)%vel(3)
          enddo
        endif
        fs(1:nobjects,iTs)=objects(1:nobjects)%T
        if (lradius_advance) fs(1:nobjects,iRs)=objects(1:nobjects)%r
      endif
!
      if (lradius_advance)  allocate(vs_normal(nobjects))
      allocate(heat_cond(nobjects))
      allocate(char_consumption(nobjects))
!
!  Prepare the solid geometry
!
      call find_solid_cell_boundaries(f)
      call calculate_shift_matrix
!
! Find nearest grid point of the "forcepoints" on all cylinders. Arrays
! are only allocated if c_dragx, c_dragy or c_dragz is set in print.in.
! This must be changed if drag force is required for other purposes,
! e.g. if solid object is allowed to move.
!
      if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
          idiag_c_dragz /= 0 .or. idiag_Nusselt /= 0 .or. &
          idiag_c_dragx_p /= 0 .or. idiag_c_dragy_p /= 0 .or. &
          idiag_c_dragz_p /= 0 .or. lchemistry) then
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
      if (fbcx(1,1) > 0) flow_dir= 1
      if (fbcx(1,2) < 0) flow_dir=-1
      if (fbcy(2,1) > 0) flow_dir= 2
      if (fbcy(2,2) < 0) flow_dir=-2
      if (fbcz(3,1) > 0) flow_dir= 3
      if (fbcz(3,2) < 0) flow_dir=-3
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
        if (lset_flow_dir) flow_dir = flow_dir_set
        if (flow_dir == 0) then
          if (lperi(1).and.lperi(2).and.lperi(3)) then
            print*,'the medium is quiescent'
          else
            call fatal_error('initialize_solid_cells','no flow direction!')
          endif
        endif
      endif
!
! Find inlet temperature
!
      if (ilnTT /= 0) then
        if (flow_dir== 1) T0=fbcx(ilnTT,1)
        if (flow_dir==-1) T0=fbcx(ilnTT,2)
        if (flow_dir== 2) T0=fbcy(ilnTT,1)
        if (flow_dir==-2) T0=fbcy(ilnTT,2)
        if (flow_dir== 3) T0=fbcz(ilnTT,1)
        if (flow_dir==-3) T0=fbcz(ilnTT,2)
        if (.not. ltemperature_nolog) T0=exp(T0)
      endif
!
      if (lroot) then
        print*,'initialize_solid_cells'
        do iobj=1, nobjects
          print*,'objects(',iobj,')%r=',objects(iobj)%r
          print*,'objects(',iobj,')%T=',objects(iobj)%T
        enddo
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
      use Initcond, only: gaunoise
      use InitialCondition, only: initial_condition_solid_cells
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: a2,rr2, wall_smoothing,rr2_low,rr2_high,shiftx,shifty
      real :: wall_smoothing_temp,xr,yr,mu1_,sumspecies
      integer :: i,j,k,cyl,jj,icyl,kchem
!
      do jj=1,ninit
      select case (initsolid_cells(jj))
!
!  This overrides any initial conditions set in the Hydro module.
!
      case ('nothing')
        if (lroot) print*,'init_solid_cells: nothing'
      case ('cylinder')
!  Initial condition for cyilinder in quiescent medium
        call gaunoise(ampl_noise,f,iux,iuz)
        shiftx=0
        if (lroot) print*,'init_solid_cells: cylinder'
        do i=l1,l2
        do j=m1,m2
        do k=n1,n2
!
!  Loop over all cylinders
!
          do icyl=1,nobjects
            a2 = objects(icyl)%r**2
            xr=x(i)-objects(icyl)%x(1)
            yr=y(j)-objects(icyl)%x(2)
            rr2 = xr**2+yr**2
            if (rr2 < a2) then
              f(i,j,k,ilnTT) = objects(icyl)%T
              f(i,j,k,isO2)  = 0.5
              f(i,j,k,isCO2) = 0.5
            endif
          enddo
        enddo
        enddo
        enddo
      case ('cylinder_combustion')
!  Initial condition for cyilinder burning in quiescent medium
        call gaunoise(ampl_noise,f,iux,iuz)
        shiftx=0
        if (lroot) print*,'init_solid_cells: cylinder_combustion'
        do i=l1,l2
        do j=m1,m2
        do k=n1,n2
!
!  Loop over all cylinders
!
          do icyl=1,nobjects
            a2 = objects(icyl)%r**2
            xr=x(i)-objects(icyl)%x(1)
            yr=y(j)-objects(icyl)%x(2)
            rr2 = xr**2+yr**2
            if (rr2 > a2) then
              wall_smoothing_temp=1-exp(-(rr2-a2)/(sqrt(a2))**2)
              f(i,j,k,ilnTT) = log(wall_smoothing_temp*exp(f(i,j,k,ilnTT))&
                         +objects(icyl)%T*(1-wall_smoothing_temp))
              if (solid_reactions_intro_time==0) then
                f(i,j,k,isCO) = wall_smoothing_temp*f(i,j,k,isCO)
                f(i,j,k,isO2) = wall_smoothing_temp*f(i,j,k,isO2)
                f(i,j,k,isCO2) = wall_smoothing_temp*f(i,j,k,isCO2)+0.2678*(1.0-wall_smoothing_temp)
                f(i,j,k,isN2) = wall_smoothing_temp*f(i,j,k,isN2)+0.7322*(1.0-wall_smoothing_temp)
              endif
            else
              f(i,j,k,ilnTT) = log(objects(icyl)%T)
              f(i,j,k,isCO)=0.0
              f(i,j,k,isO2)=0.0
              f(i,j,k,isCO2)=0.2678
              f(i,j,k,isN2)=0.7322
            endif
            sumspecies=0.0; mu1_=0.0
            do kchem=1,nchemspec
              sumspecies=sumspecies+f(i,j,k,ichemspec(kchem))
            enddo
            do kchem=1,nchemspec
              f(i,j,k,ichemspec(kchem))=f(i,j,k,ichemspec(kchem))/sumspecies
              mu1_=mu1_+f(i,j,k,ichemspec(kchem))/Mspecies(kchem)
            enddo
            if (unit_system == 'cgs') then
              Rgas_unit_sys=k_B_cgs/m_u_cgs
            endif
            f(i,j,k,ilnrho)=log(pressure0)-log(Rgas_unit_sys)-f(i,j,k,ilnTT)-log(mu1_)
          enddo
        enddo
        enddo
        enddo
      case ('cylinderstream_x')
!  Stream functions for flow around a cylinder as initial condition.
        call gaunoise(ampl_noise,f,iux,iuz)
        f(:,:,:,iux)=f(:,:,:,iux)+init_uu
        shiftx=0
        if (lroot) print*,'init_solid_cells: cylinderstream_x'
        do i=l1,l2
        do j=m1,m2
        do k=n1,n2
!
!  Loop over all cylinders
!
          do icyl=1,nobjects
            a2=objects(icyl)%r**2
            xr=x(i)-objects(icyl)%x(1)
            if (objects(icyl)%x(2) /= 0) then
              print*,'When using cylinderstream_x all cylinders must have'
              print*,'zero offset in y-direction!'
              call fatal_error('init_solid_cells:','')
            endif
            yr=y(j)
            rr2=xr**2+yr**2
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
                         +objects(icyl)%T*(1.0-wall_smoothing_temp)
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
        if (lroot) print*,'init_solid_cells: cylinderstream_y'
        shifty=0
        do i=l1,l2
        do j=m1,m2
        do k=n1,n2
          do icyl=1,ncylinders
            a2=objects(icyl)%r**2
            yr=y(j)-objects(icyl)%x(2)
            if (objects(icyl)%x(1) /= 0) then
!              write(*,*) 'DM',objects(icyl)%x(1),icyl
              print*,'When using cylinderstream_y all cylinders must have'
              print*,'zero offset in x-direction!'
              call fatal_error('init_solid_cells:','')
            endif
            xr=x(i)
            rr2=xr**2+yr**2
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
                f(i,j,k,ilnTT)=objects(icyl)%T
                f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)*f(l2,m2,n2,ilnTT)/objects(icyl)%T
              endif
            endif
          enddo
        enddo
        enddo
        enddo
        if (llast_proc_y) f(:,m2-5:m2,:,iux)=0
!
      case ('cylinder_combustion_x')
        if (lroot) print*,'init_solid_cells: cylinder_combustion_x'
        do i=l1,l2
        do j=m1,m2
        do k=n1,n2
!
!  Loop over all cylinders
!
          do icyl=1,nobjects
            a2=objects(icyl)%r**2
            xr=x(i)-objects(icyl)%x(1)
            if (objects(icyl)%x(2) /= 0) then
              print*,'When using cylinder_combustion_x all cylinders must have'
              print*,'zero offset in y-direction!'
              call fatal_error('init_solid_cells:','')
            endif
            yr=y(j)
            rr2=xr**2+yr**2
            if (rr2 > a2) then
              do cyl=0,100
                if (cyl==0) then
                  wall_smoothing=1-exp(-(rr2-a2)/skin_depth**2)
                  f(i,j,k,iux) = f(i,j,k,iux)-init_uu*&
                       2*xr*yr*a2/rr2**2*wall_smoothing
                  f(i,j,k,iuy) = f(i,j,k,iuy)+init_uu*&
                       (0. - a2/rr2 + 2*xr**2*a2/rr2**2)&
                       *wall_smoothing
                  wall_smoothing_temp=1-exp(-(rr2-a2)/(sqrt(a2))**2)
                  f(i,j,k,ilnTT) = log(wall_smoothing_temp*exp(f(i,j,k,ilnTT))&
                         +objects(icyl)%T*(1-wall_smoothing_temp))
                  if (solid_reactions_intro_time==0) then
                    f(i,j,k,isCO) = wall_smoothing_temp*f(i,j,k,isCO)
                    f(i,j,k,isO2) = wall_smoothing_temp*f(i,j,k,isO2)
                    f(i,j,k,isCO2) = wall_smoothing_temp*f(i,j,k,isCO2)+0.2678*(1.0-wall_smoothing_temp)
                    f(i,j,k,isN2) = wall_smoothing_temp*f(i,j,k,isN2)+0.7322*(1.0-wall_smoothing_temp)
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
              f(i,j,k,ilnTT) = log(objects(icyl)%T)
              f(i,j,k,isCO)=0.0
              f(i,j,k,isO2)=0.0
              f(i,j,k,isCO2)=0.2678
              f(i,j,k,isN2)=0.7322
            endif
            sumspecies=0.0; mu1_=0.0
            do kchem=1,nchemspec
              sumspecies=sumspecies+f(i,j,k,ichemspec(kchem))
            enddo
            do kchem=1,nchemspec
              f(i,j,k,ichemspec(kchem))=f(i,j,k,ichemspec(kchem))/sumspecies
              mu1_=mu1_+f(i,j,k,ichemspec(kchem))/Mspecies(kchem)
            enddo
            if (unit_system == 'cgs') then
              Rgas_unit_sys=k_B_cgs/m_u_cgs
            endif
            f(i,j,k,ilnrho)=log(pressure0)-log(Rgas_unit_sys)-f(i,j,k,ilnTT)-log(mu1_)
          enddo
        enddo
        enddo
        enddo
!
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
    integer           :: iobj, iforcepoint, ipoint, inearest, icoord(8,3)
    integer           :: ilong,ilat
    integer           :: ixl, iyl, izl, ixu, iyu, izu, ju, jl, jm
    real              :: robj, xobj, yobj, zobj,fpx, fpy, fpz
    real              :: dx1, dy1, dz1, longitude, latitude
    real              :: dist_to_fp2(8), dist_to_cent2(8), twopi,dlong,dlat
    logical           :: interiorpoint
    character(len=10) :: objectform
!
    dx1=1.0/dx
    dy1=1.0/dy
    dz1=1.0/dz
!
    twopi=2.0*pi
!
!  Loop over all objects
!
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
      endif
!
!  Loop over all forcepoints on each object, iobj
!
      do iforcepoint=1,nforcepoints
!
!  Marking whether fp is within this processor's domain or not
!
        interiorpoint = .true.
!
!  Fp coordinates. Shifting the location of the forcpoints in the
!  theta direction in order to avoid problems with autotesting
!   
        if (objectform == 'cylinder') then
          longitude = (iforcepoint-theta_shift)*dlong
          fpx = xobj - robj * sin(longitude)
          fpy = yobj - robj * cos(longitude)
          fpz = z(n1)
        elseif (objectform == 'sphere') then
!  Note definition of lines of longitude: ilong = [0,..,nlong-1]
          ilong = mod(iforcepoint-1,nlong)
!  Note definition of lines of latitude: ilat  = [0,..,nlat]
          ilat  = int((iforcepoint-1)/nlong)
          longitude = (ilong+0.5-theta_shift)*dlong
          latitude  = (ilat+0.5)*dlat
          fpx = xobj - robj*sin(longitude)*sin(latitude)
          fpy = yobj - robj*cos(longitude)*sin(latitude)
          fpz = zobj + robj*cos(latitude)
        endif
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
            if (dist_to_cent2(ipoint) > robj**2 .and. inearest == 0) then
              inearest=ipoint
            else if (dist_to_cent2(ipoint) > robj**2) then
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
    real, dimension (mx,my,mz,mfarray), intent(in):: f
    real, dimension (mx,my,mz,mvar), intent(in)   :: df
    type (pencil_case), intent(in)                :: p
    integer :: i
!
    if (ldiagnos) then
      if (idiag_c_dragx/=0 .or. idiag_c_dragy/=0 .or. idiag_c_dragz/=0 &
        .or. idiag_Nusselt/=0 .or. idiag_c_dragx_p/=0 .or. idiag_c_dragy_p/=0 &
        .or. idiag_c_dragz_p/=0) then
!
        call output_solid_cells(f,df,p)
!
!  Calculate average density of the domain, solid cell regions excluded:
!
        do i=l1,l2
          if (mod(ba(i,m,n,1),10)==0) then
            rhosum=rhosum+p%rho(i-nghost)
            irhocount=irhocount+1
          endif
        enddo
      endif
    endif
!
    call calc_solid_cells_chemistry(f,df,p)
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
    use General, only: safe_character_assign, safe_character_append
    use Mpicomm, only: mpireduce_sum, mpireduce_sum_int, mpibcast_real
!
    real    :: rhosum_all, c_dragx_all(nobjects), c_dragy_all(nobjects)
    real    :: c_dragz_all(nobjects), Nusselt_all(nobjects)
    real    :: c_dragx_p_all(nobjects), c_dragy_p_all(nobjects)
    real    :: c_dragz_p_all(nobjects), heat_all(nobjects)
    real    :: char_consumption_all(nobjects)
    integer :: irhocount_all, iobj
    real    :: norm, refrho0
    character (len=fnlen) :: file1, file2, file3
    character(len=100) :: numberstring, time_string
    character(len=200) :: solid_cell_drag, solid_cell_time
!
    call mpireduce_sum(heat_cond,heat_all,nobjects)
    call mpireduce_sum(char_consumption,char_consumption_all,nobjects)
    heat_cond(1:nobjects)=heat_all(1:nobjects)
    char_consumption(1:nobjects)=char_consumption_all(1:nobjects)
    call mpibcast_real(heat_cond,nobjects)
    call mpibcast_real(char_consumption,nobjects)
!
    dfs(1:nobjects,iTs)=dfs(1:nobjects,iTs)+heat_cond(1:nobjects) &
       /(pi*objects(1:nobjects)%r**2)/srho/scp
    if (lradius_advance) then
      dfs(1:nobjects,iRs)=dfs(1:nobjects,iRs)+char_consumption(1:nobjects)/ &
         (2.0*pi*objects(1:nobjects)%r)/srho
      vs_normal(1:nobjects)=dfs(1:nobjects,iRs)
    endif
!
    if (ldiagnos) then
      if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 &
          .or. idiag_c_dragz /= 0 .or. idiag_Nusselt /= 0 &
          .or. idiag_c_dragx_p /= 0 .or. idiag_c_dragy_p /= 0 &
          .or. idiag_c_dragz_p /= 0) then
!
!  Collect and sum rhosum, irhocount, c_dragx, c_dragz, and c_dragy.
!
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
            norm = 2.0 / (refrho0*init_uu**2)
            if (init_uu==0) norm = 1.0
            c_dragx = c_dragx_all * norm
            c_dragy = c_dragy_all * norm
            c_dragz = c_dragz_all * norm
            c_dragx_p = c_dragx_p_all * norm
            c_dragy_p = c_dragy_p_all * norm
            c_dragz_p = c_dragz_p_all * norm
!
!  Write drag coefficients for all objects (may need to expand solid_cell_drag to more
!  characters if large number of objects).
!
            call safe_character_assign(file1,trim(datadir)//'/dragcoeffs.dat')
            open(unit=81, file=file1, position='APPEND')
            write(solid_cell_drag,"(1I8,1F15.8)") it-1, t
            do iobj=1,nobjects
              write(numberstring,"(4F15.8)") c_dragx(iobj),c_dragx_p(iobj),&
                  c_dragy(iobj),c_dragy_p(iobj)
              call safe_character_append(solid_cell_drag,numberstring)
            enddo
            write(81,*) trim(solid_cell_drag)
            close(81)
!
!  Write solid data.
!
            call safe_character_assign(file2,trim(datadir)//'/timeadvance.dat')
            open(unit=83,file=file2,position='APPEND')
            write(solid_cell_time,"(1F15.8)") t
            do iobj=1,nobjects
              write(time_string,"(2F15.8)") objects(iobj)%T, objects(iobj)%r
              if (lpos_advance) then
                if (osci_dir==1) then
                  write(time_string,"(2F15.8)") objects(iobj)%x(1),objects(iobj)%vel(1)
                elseif (osci_dir==2) then
                  write(time_string,"(2F15.8)") objects(iobj)%x(2),objects(iobj)%vel(2)
                elseif (osci_dir==3) then
                  write(time_string,"(2F15.8)") objects(iobj)%x(3),objects(iobj)%vel(3)
                endif
              endif
            enddo
            call safe_character_append(solid_cell_time,time_string)
            write(83,*) trim(solid_cell_time)
            close(83)
!
!  Write char consumption rate.
!
            call safe_character_assign(file3,trim(datadir)//'/char_reaction.dat')
            open(unit=85,file=file3,position='APPEND')
            do iobj=1, nobjects
              write (85,"(4F15.8)") real(it-1), t, real(iobj), char_consumption(iobj)
            enddo
            close(85)
!
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
    use Diagnostics, only: parse_name
    use Sub
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mvar) :: f_tmp
      integer :: i,j,k,iobj,kchem,iv
      real :: z_obj,y_obj,x_obj,r_obj,r_new,r_point,dr
      real :: xghost, yghost, zghost, sumval
      logical :: bax, bay, baz, sax, say, saz
      real, dimension(3) :: xxp
      character(len=10) :: form
!
!  Find ghost points based on the mirror interpolation method
!
      do i=l1,l2
      do j=m1,m2
      do k=n1,n2
        bax=(ba(i,j,k,1)/=0).and.(ba(i,j,k,1)/=9).and.(ba(i,j,k,1)/=10)
        bay=(ba(i,j,k,2)/=0).and.(ba(i,j,k,2)/=9).and.(ba(i,j,k,2)/=10)
        if (form=='sphere') then
          baz=(ba(i,j,k,3)/=0).and.(ba(i,j,k,3)/=9).and.(ba(i,j,k,3)/=10)
        else
          baz=.false.
        endif
!
!  Check if we are in a point which must be interpolated, i.e. we are inside
!  a solid geometry AND we are not more than three grid points from the
!  closest solid-fluid interface
!
        if (bax .or. bay .or. baz) then
!  Find x, y and z values of ghost point
          xghost=x(i); yghost=y(j); zghost=z(k)
          iobj=ba(i,j,k,4)
          call interpolate_point(f,f_tmp,iobj,xghost,yghost,zghost)
          f(i,j,k,1:mvar)=f_tmp(1:mvar)
        endif
      enddo
      enddo
      enddo
!stop

!
!  Renew the flow field variables inside the solid object when these variables
!  change with time.
!
      if (lpos_advance .or. lradius_advance) then
        do i=l1,l2
        do j=m1,m2
        do k=n1,n2
          sax=(ba(i,j,k,1)==9); say=(ba(i,j,k,2)==9)
          if (form=='sphere') then
            saz=(ba(i,j,k,3)==9)
          else
            saz=.true.
          endif
          if (sax .and. say .and. saz) then
            iobj=ba(i,j,k,4)
            f(i,j,k,iux)=objects(iobj)%vel(1)
            f(i,j,k,iuy)=objects(iobj)%vel(2)
            f(i,j,k,iuz)=objects(iobj)%vel(3)
            if (ilnTT>0) then
              if (ltemperature_nolog) then
                f(i,j,k,ilnTT)=objects(iobj)%T
                f(i,j,k,ilnrho)=f(l2,m2,n2,ilnTT)/objects(iobj)%T
              else
                f(i,j,k,ilnTT)=log(objects(iobj)%T)
                f(i,j,k,ilnrho)=exp(f(l2,m2,n2,ilnTT))/objects(iobj)%T
              endif
              if (ldensity_nolog) then
                f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)*f(i,j,k,ilnrho)
              else
                f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)+log(f(i,j,k,ilnrho))
              endif
            else
              f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)
            endif
          endif
        enddo
        enddo
        enddo
      endif
!
    end subroutine update_solid_cells
!***********************************************************************
    subroutine interpolate_point(f,f_tmp,iobj,xghost,yghost,zghost)
!
!  Interpolate value in a mirror point from the eight corner values
!
      use General, only: linear_interpolate, notanumber
      use Messages, only: fatal_error
      use Sub
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mvar), intent(out) :: f_tmp
      integer, intent(in) :: iobj
      real,    intent(in) :: xghost, yghost, zghost
!
      real, dimension(3) :: o_global, ffp, xxp, ibvel, ibp
      real, dimension(mvar) :: f_mir0, f_mir1, f_ib
      real, dimension(nchemspec) :: ibchem, chem_fp0, chem_fp1, reac_rate
      real :: rs, rp1, r_sp1, rp0, r_sp0, Diff_ib, char_reac, T_fp0, T_fp1,Stefan_flow
      real :: x_obj, y_obj, z_obj, Tobj, xmir1, ymir1, zmir1
      real :: rho_gp, T_gp, rg, r_sg, xmir0, ymir0, zmir0, rho_fp0, rho_fp1
      real :: sumspecies,x_ib,y_ib,z_ib,rho_ib, mu1_fp0, mu1_fp1, mu1_gp
      real :: mu_fp0, mu_fp1, rp, coe_a, coe_b, coe_c
      integer :: lower_i, upper_i, lower_j, upper_j, lower_k, upper_k
      integer :: lower_i0,upper_i0,lower_j0,upper_j0,lower_k0,upper_k0
      integer :: lower_i1,upper_i1,lower_j1,upper_j1,lower_k1,upper_k1
      integer :: kchem
      integer, dimension(3) :: inear0, inear1, inear_ib
      character(len=10) :: form
      real    :: mu1_ib,mu_ib
      integer, dimension(3)      :: nearest_grid_index
      integer                    :: ix0,iy0,iz0
      real                       :: T_ib
      real, dimension(nchemspec) :: chem_ib
      real, dimension(nchemspec) :: Diff_ib_ks
!
      x_obj=objects(iobj)%x(1)
      y_obj=objects(iobj)%x(2)
      z_obj=objects(iobj)%x(3)
      rs=objects(iobj)%r
      Tobj=objects(iobj)%T
      form =objects(iobj)%form
      if (form=='cylinder') then
        rg=sqrt((xghost-x_obj)**2+(yghost-y_obj)**2)
        rp0=2.0*rs-rg
!
!  The distance of the probe points need to satisfy two conditions:
!  1. The first point must have four neighbouring fluid point, untill
!     we find another proper way to calculate it;
!  2. The second point cannot be 3 grid away from the boundary, which
!     is a constraint for the moving boundary case.
!
        if ((rp0-rs) <= (0.5*sqrt(2.0)*dxmin)) then
          rp0=rp0+0.5*sqrt(2.0)*dxmin
          rp1=rp0+sqrt(2.0)*dxmin
        elseif ((rp0-rs) >= (3.0*dxmin)) then
          rp1=3.0*dxmin
          rp0=rp1-sqrt(2.0)*dxmin
        else
          rp1=rp0+sqrt(2.0)*dxmin
        endif
        
        !
        xmir0=(xghost-x_obj)/rg*rp0+x_obj
        ymir0=(yghost-y_obj)/rg*rp0+y_obj
        zmir0=zghost
        !
        xmir1=(xghost-x_obj)/rg*rp1+x_obj
        ymir1=(yghost-y_obj)/rg*rp1+y_obj
        zmir1=zghost
        !
        x_ib=(xghost-x_obj)/rg*rs+x_obj
        y_ib=(yghost-y_obj)/rg*rs+y_obj
        z_ib=zghost
!
      elseif (form=='sphere') then
        rg=sqrt((xghost-x_obj)**2+(yghost-y_obj)**2+(zghost-z_obj)**2)
        rp0=2.0*rs-rg
        if ((rp0-rs) <= (0.5*sqrt(2.0)*dxmin)) then
          rp0=rp0+0.5*sqrt(2.0)*dxmin
          rp1=rp0+sqrt(2.0)*dxmin
        elseif ((rp0-rs) >= (3.0*dxmin)) then
          rp1=3.0*dxmin
          rp0=rp1-sqrt(2.0)*dxmin
        else
          rp1=rp0+sqrt(2.0)*dxmin
        endif
        !
        xmir0=(xghost-x_obj)*rp0/rg+x_obj
        ymir0=(yghost-y_obj)*rp0/rg+y_obj
        zmir0=(zghost-z_obj)*rp0/rg+z_obj
        !
        xmir1=(xghost-x_obj)/rg*rp1+x_obj
        ymir1=(yghost-y_obj)/rg*rp1+y_obj
        zmir1=(zghost-z_obj)/rg*rp1+z_obj
        !
        x_ib=(xghost-x_obj)/rg*rs+x_obj
        y_ib=(yghost-y_obj)/rg*rs+y_obj
        z_ib=zghost
!
! Check if mirror point is inside domain
!
        if (xmir0<xyz0(1) .and. lperi(1)) xmir0=xmir0+Lxyz(1)
        if (ymir0<xyz0(2) .and. lperi(2)) ymir0=ymir0+Lxyz(2)
        if (zmir0<xyz0(3) .and. lperi(3)) zmir0=zmir0+Lxyz(3)
        if (xmir0>xyz1(1) .and. lperi(1)) xmir0=xmir0-Lxyz(1)
        if (ymir0>xyz1(2) .and. lperi(2)) ymir0=ymir0-Lxyz(2)
        if (zmir0>xyz1(3) .and. lperi(3)) zmir0=zmir0-Lxyz(3)
!
        if (xmir1<xyz0(1) .and. lperi(1)) xmir1=xmir1+Lxyz(1)
        if (ymir1<xyz0(2) .and. lperi(2)) ymir1=ymir1+Lxyz(2)
        if (zmir1<xyz0(3) .and. lperi(3)) zmir1=zmir1+Lxyz(3)
        if (xmir1>xyz1(1) .and. lperi(1)) xmir1=xmir1-Lxyz(1)
        if (ymir1>xyz1(2) .and. lperi(2)) ymir1=ymir1-Lxyz(2)
        if (zmir1>xyz1(3) .and. lperi(3)) zmir1=zmir1-Lxyz(3)
!
        if (x_ib<xyz0(1) .and. lperi(1)) x_ib=x_ib+Lxyz(1)
        if (y_ib<xyz0(2) .and. lperi(2)) y_ib=y_ib+Lxyz(2)
        if (z_ib<xyz0(3) .and. lperi(3)) z_ib=z_ib+Lxyz(3)
        if (x_ib>xyz1(1) .and. lperi(1)) x_ib=x_ib-Lxyz(1)
        if (y_ib>xyz1(2) .and. lperi(2)) y_ib=y_ib-Lxyz(2)
        if (z_ib>xyz1(3) .and. lperi(3)) z_ib=z_ib-Lxyz(3)
      endif
!
!  Check that we are indeed inside the solid geometry
!
      if (rg>rs) then
        print*,'x(i),x_obj=',xghost,x_obj
        print*,'y(j),y_obj=',yghost,y_obj
        print*,'z(k),z_obj=',zghost,z_obj
        print*,'rg,r_obj,rp0,rp1=',rg,rs,rp0,rp1
        call fatal_error('interpolate_point:','rg>rs')
      endif
!
!  Find i, j and k indeces for points to be used during interpolation
!
      xxp=(/xmir0,ymir0,zmir0/); ffp=(/xmir1,ymir1,zmir1/)
      ibp=(/x_ib,y_ib,z_ib/)
!
      call find_near_indeces(lower_i,upper_i,lower_j,upper_j,lower_k,&
           upper_k,x,y,z,xxp)
      inear0=(/lower_i,lower_j,lower_k /)
!
      call find_near_indeces(lower_i1,upper_i1,lower_j1,upper_j1,lower_k1,&
           upper_k1,x,y,z,ffp)
      inear1=(/lower_i1,lower_j1,lower_k1/)
!
      call find_near_indeces(lower_i0,upper_i0,lower_j0,upper_j0,lower_k0,&
           upper_k0,x,y,z,ibp)
      inear_ib=(/lower_i0,lower_j0,lower_k0/)
!
      if (lower_i==0 .or. upper_i==0 .or. lower_i1==0 .or. upper_i1==0 &
        .or. lower_i0==0 .or. upper_i0==0) then
        call fatal_error('interpolate_point:','lower_i==0 or upper_i==0')
      endif
      if (lower_j==0 .or. upper_j==0 .or. lower_j1==0 .or. upper_j1==0 &
        .or. lower_j0==0 .or. upper_j0==0) then
        call fatal_error('interpolate_point:','lower_j==0 or upper_j==0')
      endif
      if (form=='sphere') then
        if (lower_k==0 .or. upper_k==0 .or. lower_k1==0 .or. upper_k1==0 &
          .or. lower_k0==0 .or. upper_k0==0) then
          call fatal_error('interpolate_point:','lower_k==0 or upper_k==0')
        endif
      endif
!
      o_global=objects(iobj)%x
      r_sp0=rp0-rs; r_sp1=rp1-rs; r_sg=rs-rg
!
      if (.not. linear_interpolate(f,iux,mvar,xxp,f_mir0,inear0,.false.))&
          call fatal_error('linear_interpolate','')
      if (.not. linear_interpolate(f,iux,mvar,ffp,f_mir1,inear1,.false.))&
          call fatal_error('linear_interpolate','')
      if (lfull_diff_ib) then
!
!  Find out the nearest grid point for every IB point
        call ib_nearest_grid(x_ib,y_ib,z_ib,x_obj,y_obj,z_obj,rs,nearest_grid_index)
        ix0=nearest_grid_index(1);iy0=nearest_grid_index(2);iz0=nearest_grid_index(3)
        f_ib(:)=f(ix0,iy0,iz0,iux:mvar)
      else
        if (.not. linear_interpolate(f,iux,mvar,ibp,f_ib,inear_ib,.false.))&
          call fatal_error('linear_interpolate','')
      endif
!
!  For the temperature boundaries being antisymmetric relative to the
!  object temperature are used
!
      if (ilnTT>0) then
        T_fp0=f_mir0(ilnTT); T_fp1=f_mir1(ilnTT)
        if (.not. ltemperature_nolog) then
          T_fp0=exp(T_fp0); T_fp1=exp(T_fp1)
        endif
        T_gp=((r_sg+r_sp0)*Tobj-r_sg*T_fp0)/r_sp0
        if (.not. ltemperature_nolog) then
          f_tmp(ilnTT)=log(T_gp)
        else
          f_tmp(ilnTT)=T_gp
        endif
      endif
!
!  Calculate the chemistry related variables.
!
      if (.not. ldensity_nolog) then
        rho_fp0=exp(f_mir0(ilnrho)); rho_fp1=exp(f_mir1(ilnrho))
        rho_ib=exp(f_ib(ilnrho))
      else
        rho_fp0=f_mir0(ilnrho); rho_fp1=f_mir1(ilnrho); rho_ib=f_ib(ilnrho)
      endif
!
      chem_ib(:)=f_ib(ichemspec(1:nchemspec))
      if (.not. ltemperature_nolog) then
        T_ib=exp(f_ib(ilnTT))
      else
        T_ib=f_ib(ilnTT)
      endif
!
      chem_fp0(1:nchemspec)=f_mir0(ichemspec(1:nchemspec))
      chem_fp1(1:nchemspec)=f_mir1(ichemspec(1:nchemspec))
      ibchem(1:nchemspec)=f_ib(ichemspec(1:nchemspec))*rho_ib
!
      call calc_reaction_rate(f,iobj,ibchem,reac_rate,char_reac)
!
      if (new_Stefan) then
! The species boundary conditions is based on rho*D*d_Yk/d_n+M*sum(dot{m}_i/M_i)+dot{m}_k=0
        mu1_ib=0.0
        do kchem=1, nchemspec
          mu1_ib=mu1_ib+f_ib(ichemspec(kchem))/Mspecies(kchem)
        enddo
        mu_ib=1.0/mu1_ib
!
        Stefan_flow=0.0
        do kchem=1, nchemspec
          Stefan_flow=Stefan_flow+reac_rate(kchem)/Mspecies(kchem)
        enddo
        Stefan_flow=-Stefan_flow*mu_ib
        char_reac=Stefan_flow
      endif
!
!  Calcualte diffusion coefficient at the boundary intersection point.
!
      if (lsimple_diff_ib) Diff_ib=2.58e-4/rho_ib*exp(0.7*log(Tobj/298.0))
      if (lfull_diff_ib) call calc_Diff_ib(rho_ib,T_ib,chem_ib, Diff_ib_ks)
      do kchem=1, nchemspec
        if (lfull_diff_ib) Diff_ib=Diff_ib_ks(kchem)
        if (lsecondorder_chem) then
          coe_b=(-chem_fp0(kchem)-reac_rate(kchem)/char_reac &
               +(chem_fp1(kchem)-chem_fp0(kchem))/(r_sp1**2-r_sp0**2)*r_sp0**2) &
               /(-r_sp0*r_sp1/(r_sp0+r_sp1)+rho_ib*Diff_ib/char_reac)
!         coe_b=-reac_rate(kchem)/rho_ib*Diff_ib
          coe_a=-(reac_rate(kchem)+rho_ib*Diff_ib*coe_b)/char_reac
!         coe_a=chem_fp0(kchem)+r_sp0*reac_rate(kchem)/rho_ib*Diff_ib+r_sp0**2/ &
!              (r_sp0**2-r_sp1**2)*(chem_fp0(kchem)-chem_fp1(kchem))+r_sp0**2/(r_sp0+r_sp1)* &
!              reac_rate(kchem)/rho_ib*Diff_ib
          coe_c=(chem_fp1(kchem)-chem_fp0(kchem)-coe_b*(r_sp1-r_sp0))/(r_sp1**2-r_sp0**2)
!         coe_c=(chem_fp0(kchem)-chem_fp1(kchem))/(r_sp0**2-r_sp1**2)+ &
!               reac_rate(kchem)/rho_ib*Diff_ib/(r_sp0+r_sp1)
          f_tmp(ichemspec(kchem))=coe_a-coe_b*r_sg+coe_c*r_sg**2
        else
          f_tmp(ichemspec(kchem))= &
           ((Diff_ib*rho_ib+r_sg*char_reac)*chem_fp0(kchem)+(r_sg+r_sp0)*reac_rate(kchem)) &
           /(Diff_ib*rho_ib-r_sp0*char_reac)
        endif
      enddo
!
      mu1_gp=0.0; mu1_fp0=0.0; mu1_fp1=0.0
      do kchem=1, nchemspec
        mu1_fp0=mu1_fp0+f_mir0(ichemspec(kchem))/Mspecies(kchem)
        mu1_fp1=mu1_fp1+f_mir1(ichemspec(kchem))/Mspecies(kchem)
        mu1_gp=mu1_gp+f_tmp(ichemspec(kchem))/Mspecies(kchem)
      enddo
      mu_fp0=1.0/mu1_fp0; mu_fp1=1.0/mu1_fp1
!
      if (lsecondorder_rho) then
        rho_gp=(rho_fp0*T_fp0*mu1_fp0+(rho_fp1*T_fp1*mu1_fp1-rho_fp0*T_fp0*mu1_fp0)*(r_sg**2-r_sp0**2) &
            /(r_sp1**2-r_sp0**2))/mu1_gp/T_gp
      else
        rho_gp=rho_fp0*T_fp0*mu1_fp0/T_gp/mu1_gp
      endif
!
!        if (notanumber(rho_gp) .or. rho_gp<0.0) then
!          print*,'rho_gp=',rho_gp
!          print*,'T_gp, T_fp0, T_fp1 =',T_gp, T_fp0, T_fp1 
!          print*,'mu1_gp, mu1_fp0, mu1_fp1=',mu1_gp, mu1_fp0, mu1_fp1
!          print*,'rho_fp0, rho_fp1=',rho_fp0, rho_fp1
!          print*,'r_sp1, r_sg, r_sp0=',r_sp1, r_sg, r_sp0
!          print*,'xghost,yghost=',xghost,yghost
!          print*,'f_tmp(ichemspec(:))=',f_tmp(ichemspec(:))
!          print*,'reac_rate(:)=',reac_rate(:)
!          call fatal_error('interpolate_point','ghost point density is unvalid')
!        endif
!
      if (ldensity_nolog) then
        f_tmp(ilnrho)=rho_gp
      else
        f_tmp(ilnrho)=log(rho_gp)
      endif
!
!  The boundary velocity should be calculated for every ghost point, since the
!  particle shrinking process in lradius_advance make the velocity at IB change
!  according to the location.
!
      if (lstefan_flow) then
        call calc_boundary_velocity(f,iobj,rp0,xxp,rho_ib,char_reac,ibvel)
      else
        ibvel=0.0
      endif
      f_tmp(1:3)=((r_sg+r_sp0)*ibvel-r_sg*f_mir0(1:3))/r_sp0
!
    end subroutine interpolate_point
!***************************************************************************
    subroutine find_near_indeces(lower_i,upper_i,lower_j,upper_j,&
        lower_k,upper_k,x,y,z,ppp)
!
!  Find i, j and k indeces for all neighbouring grid points
!
      integer :: ii,jj,kk
      integer, intent(out) :: lower_i,upper_i,lower_j,upper_j,lower_k,upper_k
      real, intent(in), dimension(mx) :: x
      real, intent(in), dimension(my) :: y
      real, intent(in), dimension(mz) :: z
      real, intent(in), dimension(3)  :: ppp
!
      lower_i=0
      upper_i=0
      do ii=1,mx
        if (x(ii)>ppp(1)) then
          lower_i=ii-1
          upper_i=ii
          exit
        endif
      enddo
!
      lower_j=0
      upper_j=0
      do jj=1,my
        if (y(jj)>ppp(2)) then
          lower_j=jj-1
          upper_j=jj
          exit
        endif
      enddo
!
      if (nzgrid==1) then
        lower_k=n1
        upper_k=n1
      else
        lower_k=0
        upper_k=0
        do kk=1,mz
          if (z(kk)>ppp(3)) then
            lower_k=kk-1
            upper_k=kk
            exit
          endif
        enddo
      endif
!
    end subroutine find_near_indeces
!***********************************************************************
    subroutine find_unit_vectors(p_local,rp,iobj,nr_hat,nphi_hat,ntheta_hat)
!
!  The unity vector "nr_hat" is normal to the solid surface, while
!  "nphi_hat" and "ntheta_hat" are the unit vectors in the two angular
!  directions. The angle "theta" is zero in the positive x-direction,
!  while "phi" is zero in the positive z-direction.
!
      real, dimension(3) :: p_local, nr_hat, ntheta_hat, nphi_hat
      integer :: iobj
      real :: phi, theta, rp
!
      intent(in) :: p_local,rp,iobj
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
      if (objects(iobj)%form=='cylinder') then
        nr_hat    =(/cos(theta),sin(theta),0.0/)
        nphi_hat  =(/0.0, 0.0, 1.0/)
        ntheta_hat=(/-sin(theta),cos(theta),0.0/)
      else
        nr_hat    =(/cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)/)
        nphi_hat  =(/-cos(phi)*cos(theta),-cos(phi)*sin(theta),sin(phi)/)
        ntheta_hat=(/-sin(theta),cos(theta),0.0/)
      endif
!
    end subroutine find_unit_vectors
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
    function in_solid_cell(part_pos,part_rad)
!
!  Check if the position px,py,pz is within a solid cell
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
        if ((ba(i,m,n,1)/=0) .or. (ba(i,m,n,2)/=0) .or. &
            (ba(i,m,n,3)/=0)) then
!
!  If this is a fluid point which has to be interpolated because it is very
!  close to the solid geometry (i.e. ba(i,m,n,1) == 10) then only the
!  temperature and the velocity components should be frozen.
!
          if (ba(i,m,n,1) == 10) then
            df(i,m,n,iux:iuz)=0
            if (ilnTT>0) df(i,m,n,ilnTT)=0
          else
            df(i,m,n,:)=0
          endif
        endif
      enddo
!
    endsubroutine freeze_solid_cells
!***********************************************************************
    subroutine read_solid_cells_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      !read(parallel_unit, NML=solid_cells_init_pars, IOSTAT=iostat)
      iostat = 0
      read(parallel_unit, NML=solid_cells_init_pars)
!
    endsubroutine read_solid_cells_init_pars
!***********************************************************************
    subroutine read_solid_cells_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=solid_cells_run_pars, IOSTAT=iostat)
!
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
!                       surface of the solid geometry that we must set the
!                       value of this point by some special method
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
!  This depends on the form of the solid geometry
!
          if (form=='cylinder') then
            x2=objects(iobj)%r**2-(y(j)-objects(iobj)%x(2))**2
          else if (form=='sphere') then
            x2=objects(iobj)%r**2-(y(j)-objects(iobj)%x(2))**2&
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
        if (form /= 'cylinder') then
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
                if ((ba(i,j,k,1)/=0) .and. (ba(i,j,k,1)/=10)) then
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
        if ( .false. ) then
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
                  if ((ba(i,j,k,1) /= 0 ) .and. (ba(i,j,k,1) /= 10)) then
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
    end subroutine find_solid_cell_boundaries
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
    end subroutine calculate_shift_matrix
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
    end subroutine find_closest_wall
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
      ba_defined=.false.
!
    end function ba_defined
!***********************************************************************
    subroutine find_point(rij,rs,f,yin,xout,xmin,xmax,min,fout,x0,surf_val)
!
!  20-mar-2009/nils: coded
!
!  Check if a grid line has any of it ends inside a solid cell - if so
!  find the point where the grid line enters the solid cell.
!
      integer, intent(in) :: min
      real, intent(in) :: xmin,xmax,rij,rs,f,yin,x0,surf_val
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
        fout=surf_val
      endif
!
    end subroutine find_point
!***********************************************************************
    subroutine pencil_criteria_solid_cells()
!
!  All pencils that the Solid_Cells module depends on are specified here.
!
!  mar-2009/kragset: coded
!
!  Request p and sij-pencils here
!  Request rho-pencil
      lpenc_requested(i_pp)=.true.
      lpenc_requested(i_sij)=.true.
      lpenc_requested(i_rho)=.true.
      lpenc_requested(i_gTT)=.true.
      lpenc_requested(i_TT)=.true.
      lpenc_requested(i_lambda)=.true.
      if (lchemistry) lpenc_requested(i_hhk_full)=.true.
!
    end subroutine pencil_criteria_solid_cells
!***********************************************************************
    subroutine solid_cells_clean_up
!
!  Deallocate the variables allocated in solid_cells
!
!  7-oct-2010/dhruba: adeped from hydro_kinematic
!  21-jul-2011/bing: fixed, only deallocate variable if allocted
!
      print*, 'Deallocating some solid_cells variables ...'
      if (allocated(fpnearestgrid)) deallocate(fpnearestgrid)
      if (allocated(c_dragx)) deallocate(c_dragx)
      if (allocated(c_dragy)) deallocate(c_dragy)
      if (allocated(c_dragz)) deallocate(c_dragz)
      if (allocated(c_dragx_p)) deallocate(c_dragx_p)
      if (allocated(c_dragy_p)) deallocate(c_dragy_p)
      if (allocated(c_dragz_p)) deallocate(c_dragz_p)
      if (allocated(Nusselt)) deallocate(Nusselt)
      print*, '..Done.'
!
    end subroutine solid_cells_clean_up
!***********************************************************************
    subroutine linear_interpolate_quadratic(gp,A_corners,x_corners,xxp)
!
!  Interpolate the value of g to arbitrary (xp, yp, zp) coordinate
!  using the linear interpolation formula
!
!    g(x,y,z) = A*x*y*z + B*x*y + C*x*z + D*y*z + E*x + F*y + G*z + H .
!
!  The coefficients are determined by the 8 grid points surrounding the
!  interpolation point.
!
!  21-may-31/NILS: Adapted form the version in general.f90
!
      use Cdata
!
      real, dimension(2,2,2,3) :: x_corners, A_corners
      real, dimension (3) :: xxp, A_p, gp
      real, dimension (3) :: g1, g2, g3, g4, g5, g6, g7, g8
      real :: xp0, yp0, zp0, drr
      real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
      integer :: i
      logical :: lfirstcall=.true.
!
      intent(in)  :: A_corners,x_corners,xxp
      intent(out) :: gp
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!  Set it equal to 0 for dimensions having 1 grid points; this will make sure
!  that the interpolation is bilinear for 2D grids.
!
      xp0=0; yp0=0; zp0=0
      if (nxgrid/=1) xp0=xxp(1)-x_corners(1,1,1,1)
      if (nygrid/=1) yp0=xxp(2)-x_corners(1,1,1,2)
      if (nzgrid/=1) zp0=xxp(3)-x_corners(1,1,1,3)
!
!  Inverse grid spacing
!
      if ( (.not. all(lequidist)) .or. lfirstcall) then
        dx1=0
        dy1=0
        dz1=0
        if (nxgrid/=1) dx1=1/(x_corners(2,1,1,1)-x_corners(1,1,1,1))
        if (nygrid/=1) dy1=1/(x_corners(1,2,1,2)-x_corners(1,1,1,2))
        if (nzgrid/=1) dz1=1/(x_corners(1,1,2,3)-x_corners(1,1,1,3))
        dxdy1=dx1*dy1; dxdz1=dx1*dz1; dydz1=dy1*dz1
        dxdydz1=dx1*dy1*dz1
      endif
!
!  Function values at all corners.
!
      g1=A_corners(1,1,1,1:3)
      g2=A_corners(2,1,1,1:3)
      g3=A_corners(1,2,1,1:3)
      g4=A_corners(2,2,1,1:3)
      g5=A_corners(1,1,2,1:3)
      g6=A_corners(2,1,2,1:3)
      g7=A_corners(1,2,2,1:3)
      g8=A_corners(2,2,2,1:3)
!
!  Interpolation formula.
!
      gp = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) + &
          xp0*yp0*dxdy1*(g1-g2-g3+g4) + xp0*zp0*dxdz1*(g1-g2-g5+g6) + &
          yp0*zp0*dydz1*(g1-g3-g5+g7) + &
          xp0*yp0*zp0*dxdydz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
!
      if (lfirstcall) lfirstcall=.false.
!
    end subroutine linear_interpolate_quadratic
!***********************************************************************
    subroutine r_theta_phi_velocity_in_point(f,g_global,inear,iobj,o_global,&
        rs,r_sg,v_r,v_theta,v_phi)
!
!  Find values of the velocity in the r, phi and theta directions for
!  any given position (does not have to be a grid point).
!
      USE sub
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension(3) :: o_global,g_global
      integer, dimension(3) :: inear
      real :: rs, r_sg, drr,rpp
      integer :: i,j,k,iobj
      real,  dimension(3) :: xc,A_g,xc_local
      real,  dimension(3) :: nrc_hat, nthetac_hat, nphic_hat
      real, dimension(2,2,2,3) :: x_corners, A_corners
      real :: v_r, v_phi, v_theta
      real :: vc_r, vc_phi, vc_theta

!
!  For all the 8 corner points the velocity in the r, phi and theta
!  directions must be found. These are then used to found a scaling
!  coefficient "A" for all three directions.
!
      do k=inear(3),inear(3)+1
      do j=inear(2),inear(2)+1
      do i=inear(1),inear(1)+1
        xc=(/x(i),y(j),z(k)/)
        xc_local=xc-o_global
        if (objects(iobj)%form=='cylinder') then
          rpp=sqrt((x(i)-o_global(1))**2+(y(j)-o_global(2))**2)
        elseif (objects(iobj)%form=='sphere') then
          rpp=sqrt((x(i)-o_global(1))**2+(y(j)-o_global(2))**2+&
              (z(k)-o_global(3))**2)
        endif
        drr=rpp-rs
        call find_unit_vectors(xc_local,rpp,iobj,nrc_hat,nphic_hat,nthetac_hat)
        call dot(nrc_hat    ,f(i,j,k,iux:iuz),vc_r)
        call dot(nphic_hat  ,f(i,j,k,iux:iuz),vc_phi)
        call dot(nthetac_hat,f(i,j,k,iux:iuz),vc_theta)
        A_corners(i-inear(1)+1,j-inear(2)+1,k-inear(3)+1,1)=vc_r/drr**2
        A_corners(i-inear(1)+1,j-inear(2)+1,k-inear(3)+1,2)=vc_phi/drr
        A_corners(i-inear(1)+1,j-inear(2)+1,k-inear(3)+1,3)=vc_theta/drr
        x_corners(i-inear(1)+1,j-inear(2)+1,k-inear(3)+1,:)=xc
      enddo
      enddo
      enddo
!
!  Having found the scaling component for all 8 corner points and for all three
!  velocity components (stored in "A_corners")
!  we can now find the interpolated velocity of the
!  three scaling components in point "g".
!
      call linear_interpolate_quadratic(A_g,A_corners,x_corners,g_global)
!
!  Since we now know "A_g", the value of the three scaling components in
!  the point "g" we can use "A_g" to find the three velocity components of
!  interest in "g".
!
      v_r    =A_g(1)*r_sg**2
      v_phi  =A_g(2)*r_sg
      v_theta=A_g(3)*r_sg
!
    end subroutine r_theta_phi_velocity_in_point
!***********************************************************************
    subroutine solid_cells_timestep_first(f)
!
!  Setup dfs in the beginning of each itsub.
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
      if (lfirst) then
        dfs(1:nobjects,1:nsvar)=0.0
      else
        dfs(1:nobjects,1:nsvar)=alpha_ts(itsub)*dfs(1:nobjects,1:nsvar)
      endif
!
    end subroutine solid_cells_timestep_first
!***********************************************************************
    subroutine solid_cells_timestep_second(f,int_dt,int_ds)
!
!  Time evolution of solid_cells variables.
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: iobj
      real :: int_dt, int_ds
!
      solid_dt=int_dt; solid_ds=int_ds
!
      fs(1:nobjects,iTs)=fs(1:nobjects,iTs)+solid_dt*dfs(1:nobjects,iTs)
!
      if (lpos_advance) call advance_solid_cells_position()
!
      if (lpos_advance) then
        do iobj=1, nobjects
          objects(iobj)%x(1)=fs(iobj,ixs)
          objects(iobj)%x(2)=fs(iobj,iys)
          objects(iobj)%x(3)=fs(iobj,izs)
          objects(iobj)%vel(1)=fs(iobj,ius)
          objects(iobj)%vel(2)=fs(iobj,ivs)
          objects(iobj)%vel(3)=fs(iobj,iws)
        enddo
      endif
!
      if (lradius_advance) objects(1:nobjects)%r=fs(1:nobjects,iRs)
!
!  Reprepare the solid geometry
!
      if (lpos_advance .or. lradius_advance) then
        call find_solid_cell_boundaries(f)
        call calculate_shift_matrix
        call fp_nearest_grid
      endif
!
      objects(1:nobjects)%T=fs(1:nobjects,iTs)
      heat_cond=0.0 
      char_consumption=0.0
!
    endsubroutine solid_cells_timestep_second
!***********************************************************************
   subroutine advance_solid_cells_position()
!
!  Time evolution of solid_cells position and velocity.
!
     integer :: iobj
!
     do iobj=1, nobjects
       fs(iobj,ius)=0.0; fs(iobj,ixs)=0.0
       fs(iobj,ivs)=0.0; fs(iobj,iys)=0.0
       fs(iobj,iws)=0.0; fs(iobj,izs)=0.0
!
       if (osci_dir==1) then
         fs(iobj,ius)=2.0*pi*osci_f*osci_A*cos(2.0*pi*osci_f*(t-osci_t))
         fs(iobj,ixs)=osci_A*sin(2.0*pi*osci_f*(t-osci_t))
       endif
!
       if (osci_dir==2) then
         fs(iobj,ivs)=2.0*pi*osci_f*osci_A*cos(2.0*pi*osci_f*(t-osci_t))
         fs(iobj,iys)=osci_A*sin(2.0*pi*osci_f*(t-osci_t))
       endif
!
       if (osci_dir==3) then
         fs(iobj,iws)=2.0*pi*osci_f*osci_A*cos(2.0*pi*osci_f*(t-osci_t))
         fs(iobj,izs)=osci_A*sin(2.0*pi*osci_f*(t-osci_t))
       endif
!
     enddo
!
   end subroutine advance_solid_cells_position
!*************************************************************************
   subroutine calc_solid_cells_chemistry(f,df,p)
!
!  Time evolution of solid_cells temperature.
!
     use Sub, only: dot
!
     real, dimension(mx,my,mz,mfarray), intent(in) :: f
     real, dimension (mx,my,mz,mvar), intent(in)   :: df
     type (pencil_case), intent(in)                :: p
!
     real :: xobj, yobj, zobj, robj, Tobj, rforce, twopi, latitude
     real :: surfacecoeff, surfaceelement, dlong, dlat, longitude
     real :: normal_gradT ,reac_heat, fp_heat, char_reac
     real, dimension(3) :: nvec
     real, dimension(nchemspec) :: fpchem, reac_rate
     integer, dimension(3) :: inear
     integer :: iobj, ifp, ilong, ilat, ix0, iy0, iz0, kchem
     character(len=10) :: objectform
!
     twopi=2.0*pi
!
     do iobj=1, nobjects
       robj=objects(iobj)%r
       xobj=objects(iobj)%x(1)
       yobj=objects(iobj)%x(2)
       objectform=objects(iobj)%form
       Tobj=objects(iobj)%T
       rforce=robj+dxmin*ineargridshift
       if (objectform=='cylinder') then
         zobj=z(n1)
         dlong=twopi/nforcepoints
         surfaceelement=dlong*rforce
       elseif (objectform=='sphere') then
         zobj=objects(iobj)%x(3)
         dlong=twopi/nlong
         dlat=pi/(nlat+1)
         surfacecoeff=dlong*dlat*rforce**2
       else
         print*,"Warning: Subroutine advance_solid_cells not implemented ",&
             "for this objectform."
       endif
       do ifp=1, nforcepoints
         iy0=fpnearestgrid(iobj,ifp,2)
         iz0=fpnearestgrid(iobj,ifp,3)
         if (objectform=='cylinder') iz0=n
!
         if (iy0 == m .and. iz0 == n) then
           ix0=fpnearestgrid(iobj,ifp,1)
           if (ix0 >= l1 .and. ix0 <= l2) then
             if (objectform=='cylinder') then
               longitude=(ifp-theta_shift)*dlong
               nvec(1)=-sin(longitude)
               nvec(2)=-cos(longitude)
               nvec(3)=0
             elseif (objectform=='sphere') then
               ilong=mod(ifp-1,nlong)
               ilat =int((ifp-1)/nlong)
               longitude=(ilong+0.5-theta_shift)*dlong
               latitude=(ilat+0.5)*dlat
               nvec(1)=-sin(longitude)*sin(latitude)
               nvec(2)=-cos(longitude)*sin(latitude)
               nvec(3)=cos(latitude)
               surfaceelement=surfacecoeff*sin(latitude)
             else
               call fatal_error('advance_solid_cells','No such objectform!')
               call keep_compiler_quiet(nvec)
             endif
!  Conductivity heat transfer.
             call dot(p%gTT(ix0-nghost,:), nvec, normal_gradT)
!  Reaction heat.
             fpchem(1:nchemspec)=f(ix0,iy0,iz0,ichemspec(1):ichemspec(nchemspec))
             if (ldensity_nolog) then
               fpchem(1:nchemspec)=fpchem(1:nchemspec)*f(ix0,iy0,iz0,ilnrho)
             else
               fpchem(1:nchemspec)=fpchem(1:nchemspec)*exp(f(ix0,iy0,iz0,ilnrho))
             endif
             call calc_reaction_rate(f,iobj,fpchem,reac_rate,char_reac)
             reac_heat=0.0
             do kchem=1,nchemspec
               reac_heat=reac_heat+reac_rate(kchem)*p%hhk_full(ix0-nghost,kchem)
             enddo
!  Add radiation heat.
             fp_heat=p%lambda(ix0-nghost)*normal_gradT-reac_heat & 
                 -sigmaSB_cgs*(Tobj**4-T0**4)
!  Sum heat.
             heat_cond(iobj)=heat_cond(iobj)+fp_heat*surfaceelement
!
!  Calculate char consumption.
!
             char_consumption(iobj)=char_consumption(iobj)+ &
                  char_reac*surfaceelement
           endif
         endif
       enddo
     enddo
!
   end subroutine calc_solid_cells_chemistry
!***********************************************************************
   subroutine calc_boundary_velocity(f,iobj,rp,xxp,rho_loc,char_reac,ibvel)
!
!  This subroutine is designed to calculate the velocity at the IB points. 
!  Since the rotation, translation and shrink is involved, it becomes more
!  complicated to determined the velocity values at the IB points.
!
     real, dimension(mx,my,mz,mfarray), intent(in) :: f
     integer, intent(in) :: iobj
     real, intent(in) :: rp, rho_loc
     real, dimension(3), intent(in)  :: xxp
     real, dimension(3), intent(out) :: ibvel
     real, intent(in) :: char_reac
!
     real :: vs_theta, vs_phi, vs_r, mu1_
     real, dimension(3) :: o_global, xxp_local, nr_hat, nphi_hat, ntheta_hat
     integer :: kchem
!
     ibvel=0.0
!
     o_global=objects(iobj)%x
     xxp_local=xxp-o_global
     call find_unit_vectors(xxp_local,rp,iobj,nr_hat,nphi_hat,ntheta_hat)
!
!  Boundary velocity due to solid rotation.
!
     vs_theta=objects(iobj)%rot(1)
     vs_phi=objects(iobj)%rot(2)
!
!  Boundary velocity due to Stefan flow.
!
     vs_r=-char_reac/rho_loc
!
     if (lradius_advance) vs_r=vs_r+vs_normal(iobj)
!
     ibvel=vs_r*nr_hat+vs_theta*ntheta_hat+vs_phi*nphi_hat
!
!  Boundary velocity due to solid transform.
!
     if (lpos_advance) then
       ibvel=ibvel+objects(iobj)%vel
     endif
!
   end subroutine calc_boundary_velocity
!**************************************************************************
   subroutine calc_reaction_rate(f,iobj,ibchem,reac_rate,char_reac)
!
!  This is designed for calculating the heterogenerous reactions on the
!  char particle surface. Three heterogeneous reactions are included:
! ---------------------------------------------------------------------
! | Number |     Reactions     |    B_n    |  alpha_n  |  E_an(J/mol) |
! ---------------------------------------------------------------------
! |   (1)  |   2C + O2 => 2CO  |   1.97e7  |    0.0    |   197997.91  |
! ---------------------------------------------------------------------
! |   (2)  |   C + CO2 => 2CO  |  1.291e5  |    0.0    |  191022.464  |
! ---------------------------------------------------------------------
!
! ---------------------------------------------------------------------
! | Number |     Reactions     |    B_n    |  alpha_n  |  E_an(J/mol) |
! ---------------------------------------------------------------------
! |   (1)  |   2C + O2 => 2CO  |   3.007e5  |    0.0    |   149370  |
! ---------------------------------------------------------------------
! |   (2)  |   C + CO2 => 2CO  |   4.016e8  |    0.0    |   247670  |
! ---------------------------------------------------------------------
!
!  The reaction rate is expressed using Arrhenius equation:
!     k_i=B_n(i)*T**alpha_n(i)*exp(-E_an(i)/Rcal/T_loc)
!
     real, dimension(mx,my,mz,mfarray), intent(in) :: f
     real, dimension(nchemspec), intent(in) :: ibchem
     integer, intent(in) :: iobj
     real, intent(out), optional :: char_reac
     real, dimension(nchemspec), intent(out) :: reac_rate
!
     real :: Rcal, Rcal1, T_loc, T_loc1
     real, dimension(2) :: alpha, prod, vreact_p, kf, B_n, E_an
     integer :: kchem, i
!
     alpha=0.0
     B_n=(/1.97e9,1.291e7/)  ! cm/s
     E_an=(/47301.701,45635.267/)  ! cal/mol
!    B_n=(/3.007e7,4.016e10/)  ! cm/s
!    E_an=(/35684.493,59168.363/)  ! cal/mol
!
     if (unit_system == 'cgs') then
       Rgas_unit_sys=k_B_cgs/m_u_cgs
     endif
     Rcal=Rgas_unit_sys/4.14*1e-7
     Rcal1=1.0/Rcal
     T_loc=objects(iobj)%T
     T_loc1=1.0/T_loc
!
     prod(1)=ibchem(ichemsO2)/Mspecies(ichemsO2)
     prod(2)=ibchem(ichemsCO2)/Mspecies(ichemsCO2)
!
!  The mass fraction of O2 and CO2 could be negative because of large
!  negative values at the ghost points. If a large negative value appears
!  at the ghost point species value, it means that the concentration of this
!  species has been zero at the boundary location.
!
     if (prod(1)<0.0) then
       prod(1)=0.0
     endif
     if (prod(2)<0.0) then
       prod(2)=0.0
     endif
!
!  Calculate the forward reaction rate. The reactions are irrevisible,
!  so the backward reaction rate is not calculated.
!
     do i=1,2
       kf(i)=log(B_n(i))+alpha(i)*log(T_loc)-E_an(i)*Rcal1*T_loc1
       vreact_p(i)=prod(i)*exp(kf(i))
     enddo
!
!  Because the three reactions are irrevisible, vreact=vreact_p.
!  Calculate production rate for all species k (called \dot(\omega)_k
!  in the chemkin manual)
!
     reac_rate(ichemsO2)=-vreact_p(1)*Mspecies(ichemsO2)
     reac_rate(ichemsCO)=(2.0*vreact_p(1)+2.0*vreact_p(2))*Mspecies(ichemsCO)
     reac_rate(ichemsCO2)=-vreact_p(2)*Mspecies(ichemsCO2)
     reac_rate(ichemsN2)=0.0
!
     reac_rate=reac_rate*unit_time
!
     if (present(char_reac)) then
       char_reac=-(2.0*vreact_p(1)+vreact_p(2))*12.0*unit_time
     endif
!
!  Gradually turn solid reactions on
!
     if (t < solid_reactions_intro_time) then
       reac_rate=reac_rate*t/solid_reactions_intro_time
       if (present(char_reac)) then
         char_reac=char_reac*t/solid_reactions_intro_time
       endif
     endif
!
   end subroutine calc_reaction_rate
!**************************************************************************
   subroutine read_snapshot_solid_cells(chsnap)
!
     character(len=*) :: chsnap
     character(len=fnlen) :: filename, dfilename
     intent(in) :: chsnap
!
     filename =trim(chsnap)//'/svar.dat'
     dfilename=trim(chsnap)//'/dsvar.dat'
!
     open(1,FILE=filename,FORM='formatted')
     read(1,*) fs(1:nobjects,1:nsvar)
     if (ip<=8) print*, 'input_particles: read ', filename
     close(1)
!
     open(2,FILE=dfilename,FORM='formatted')
     read(2,*) dfs(1:nobjects,1:nsvar)
     if (ip<=8) print*, 'input_particles: read ', dfilename
     close(2)
!
     if (lpos_advance) then
       objects(1:nobjects)%x(1)=fs(1:nobjects,ixs)
       objects(1:nobjects)%x(2)=fs(1:nobjects,iys)
       objects(1:nobjects)%x(3)=fs(1:nobjects,izs)
       objects(1:nobjects)%vel(1)=fs(1:nobjects,ius)
       objects(1:nobjects)%vel(2)=fs(1:nobjects,ivs)
       objects(1:nobjects)%vel(3)=fs(1:nobjects,iws)
     endif
     objects(1:nobjects)%T=fs(1:nobjects,iTs)
     if (lradius_advance) objects(1:nobjects)%r=fs(1:nobjects,iRs)
!
   end subroutine read_snapshot_solid_cells
!************************************************************************
   subroutine write_snapshot_solid_cells(chsnap,f)
!
     character(len=*), intent(in) :: chsnap
     character(len=fnlen) :: filename, dfilename
     real, dimension(mx,my,mz,mfarray) :: f
!
     filename =trim(chsnap)//'/svar.dat'
     dfilename=trim(chsnap)//'/dsvar.dat'
!
     call wsnap_solid_cells(filename,f,fs)
     call wsnap_solid_cells(dfilename,f,dfs)
!
   end subroutine write_snapshot_solid_cells
!************************************************************************
   subroutine wsnap_solid_cells(chsnap,f,fs)
!
     use IO, only : lun_output
!
     character(len=*), intent(in) :: chsnap
     real, dimension(mx,my,mz,mfarray) :: f
     real, dimension(:,:),intent(in) :: fs
!
     open(lun_output,FILE=chsnap)
     write(lun_output,*) fs(1:nobjects,1:nsvar)
     close(lun_output)
!
   end subroutine wsnap_solid_cells
!************************************************************************
   subroutine output_solid_cells(f,df,p)
!   
     use Sub, only: dot
     use Viscosity, only: getnu
     use General, only: safe_character_assign,safe_character_append
     use General, only: linear_interpolate
     use General, only: itoa
!
     real, dimension (mx,my,mz,mfarray), intent(in):: f
     real, dimension (mx,my,mz,mvar), intent(in)   :: df
     type (pencil_case), intent(in)                :: p
!
     real    :: fp_pressure, fp_gradT
     real    :: fp_stress(3,3)
     integer :: iobj, ifp, ix0, iy0, iz0, i, ilong, ilat
     integer :: lower_i, lower_j, lower_k, upper_i, upper_j, upper_k
     real    :: longitude, latitude, dlong, dlat, robj, rforce
     real, dimension(nx) :: nu, twonu
     real    :: force_x, force_y, force_z, loc_Nus, loc_density
     real    :: twopi, nvec(3), surfaceelement, surfacecoeff
     real    :: deltaT, Tobj, drag_norm, nusselt_norm, char_reac, prodout(2)
     real    :: fpx, fpy, fpz, xobj, yobj, zobj, theta, phi, c_press
     real, dimension(3) :: fp_location, fp_vel, vel_error
     real, dimension(4) :: chemspec_error
     real, dimension(1) :: rho_error,t_error
     real, dimension(nchemspec) :: fpchem, reac_rate
     integer, dimension(3) :: inear_error
     character (len=fnlen) :: file1, file2, file3
     character(len=10) :: objectform
     character(len=60) :: Nus_string
     character(len=100):: error_string, reac_string
     character(len=120):: solid_cell_Nus
     character(len=150):: solid_cell_error, solid_cell_reac
     character(len=10) :: chproc
!
!  Reset cumulating quantities before calculations in first pencil
!
     if (imn == 1) then
       if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
           idiag_c_dragz /= 0 .or. idiag_c_dragx_p /= 0 .or. &
           idiag_c_dragy_p /= 0 .or. idiag_c_dragz_p /= 0) then
         c_dragx=0.0; c_dragy=0.0; c_dragz=0.0
         c_dragx_p=0.0; c_dragy_p=0.0; c_dragz_p=0.0
       endif
       if (idiag_Nusselt /= 0) Nusselt=0.0
       rhosum=0
       irhocount=0
     endif
!
     call getnu(nu_pencil=nu,p=p)
     twopi=2.0*pi
     twonu=2.0*nu
!
     do iobj=1,nobjects
       robj = objects(iobj)%r
       xobj = objects(iobj)%x(1)
       yobj = objects(iobj)%x(2)
       objectform = objects(iobj)%form
       rforce = robj+dxmin*ineargridshift
       if (objectform=='cylinder') then
         zobj = z(n1)
         dlong = twopi/nforcepoints
         surfaceelement = dlong*rforce
         drag_norm=1.0/(2.0*rforce)
         nusselt_norm=1.0/(twopi*rforce)
       elseif (objectform=='sphere') then
         zobj  = objects(iobj)%x(3)
         dlong = twopi/nlong
         dlat  = pi/(nlat+1)
         surfacecoeff = dlong*dlat*rforce**2
         drag_norm=1.0/(pi*rforce**2)
         nusselt_norm=1.0/(4*pi*rforce**2)
       else
         print*, "Warning: Subroutine dsolid_dt not implemented ", &
             "for this objectform."
       endif
!
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
               theta=(longitude/pi-0.5)*180.0
               nvec(1) = -sin(longitude)
               nvec(2) = -cos(longitude)
               nvec(3) = 0
               fpx=xobj+robj*nvec(1)
               fpy=yobj+robj*nvec(2)
               fpz=z(n1)
             elseif (objectform == 'sphere') then
               ilong=mod(ifp-1,nlong)
               ilat =int((ifp-1)/nlong)
               longitude=(ilong+0.5-theta_shift)*dlong
               theta=(longitude/pi-0.5)*180.0
               latitude=(ilat+0.5)*dlat
               phi=latitude/pi*180.0
               nvec(1)=-sin(longitude)*sin(latitude)
               nvec(2)=-cos(longitude)*sin(latitude)
               nvec(3)=cos(latitude)
               surfaceelement=surfacecoeff*sin(latitude)
               fpx=xobj+robj*nvec(1)
               fpy=yobj+robj*nvec(2)
               fpz=zobj+robj*nvec(3)
             else
               call fatal_error('dsolid_dt','No such objectform!')
               call keep_compiler_quiet(nvec)
             endif
             fp_location=(/fpx, fpy, fpz/)
!
! Find force in x,y and z direction
!
             if (idiag_c_dragx /= 0 .or. idiag_c_dragy /= 0 .or. &
                 idiag_c_dragz /= 0 .or. idiag_c_dragx_p /= 0 .or. &
                 idiag_c_dragy_p /= 0 .or. idiag_c_dragz_p /= 0) then
!
!  Calculate pressure on the immersed boundary in a bilinear way.
!
               fp_pressure=p%pp(ix0-nghost)
               loc_density=p%rho(ix0-nghost)
               fp_stress(:,:)=twonu(ix0-nghost)*p%rho(ix0-nghost)*p%sij(ix0-nghost,:,:)
               fp_vel(1:3)=f(ix0,m,n,iux:iuz)
!
!  Force in x-,y-, and z-directions
!
               force_x = (-fp_pressure*nvec(1) &
                   +fp_stress(1,1)*nvec(1)+fp_stress(1,2)*nvec(2)+fp_stress(1,3)*nvec(3) &
                   ) * surfaceelement
!
               force_y = (-fp_pressure*nvec(2) &
                   +fp_stress(2,1)*nvec(1)+fp_stress(2,2)*nvec(2)+fp_stress(2,3)*nvec(3) &
                   ) * surfaceelement
!
               force_z = (-fp_pressure*nvec(3) &
                   +fp_stress(3,1)*nvec(1)+fp_stress(3,2)*nvec(2)+fp_stress(3,3)*nvec(3) &
                   ) * surfaceelement
!
               if (lpos_advance) then
                 force_x=force_x+(fp_vel(1)*nvec(1)+fp_vel(2)*nvec(2)+fp_vel(3)*nvec(3)) &
                   *loc_density*fp_vel(1)*surfaceelement
!
                 force_y=force_y+(fp_vel(1)*nvec(1)+fp_vel(2)*nvec(2)+fp_vel(3)*nvec(3)) &
                   *loc_density*fp_vel(2)*surfaceelement
!
                 force_z=force_z+(fp_vel(1)*nvec(1)+fp_vel(2)*nvec(2)+fp_vel(3)*nvec(3)) &
                   *loc_density*fp_vel(3)*surfaceelement
               endif
!
               c_dragx(iobj) = c_dragx(iobj) + force_x * drag_norm
               c_dragy(iobj) = c_dragy(iobj) + force_y * drag_norm
               c_dragz(iobj) = c_dragz(iobj) + force_z * drag_norm
!
               c_dragx_p(iobj)=c_dragx_p(iobj)-fp_pressure*nvec(1)*drag_norm*surfaceelement
               c_dragy_p(iobj)=c_dragy_p(iobj)-fp_pressure*nvec(2)*drag_norm*surfaceelement
               c_dragz_p(iobj)=c_dragz_p(iobj)-fp_pressure*nvec(3)*drag_norm*surfaceelement
!
               c_press=(fp_pressure-pressure0)/(0.5*init_uu**2)
             endif
!
!  Local Nusselt number
!
             if (idiag_Nusselt /= 0) then
               call dot(p%gTT(ix0-nghost,:),-nvec,fp_gradT)
               Tobj=objects(iobj)%T
               if (.not. ltemperature_nolog) Tobj=exp(Tobj)
               deltaT=Tobj-T0
               loc_Nus=fp_gradT*robj*2.0/deltaT
               Nusselt(iobj)=Nusselt(iobj)+loc_Nus/nforcepoints
             endif
!
             chproc=itoa(iproc)
!
             if(lNusselt_output) then
               call safe_character_assign(file1,trim(datadir)//'/local_Nusselt.dat'//chproc)
               open(unit=87,file=file1,position='APPEND')
               write(solid_cell_Nus,98) it-1, t
               if (ilnTT>0) then
                 write(Nus_string,96)  theta, c_press, loc_Nus
               else
                 write(Nus_string,"(2F15.8)")  theta, c_press
               endif
               call safe_character_append(solid_cell_Nus,Nus_string)
               write(87,*) trim(solid_cell_Nus)
               close(87)
             endif
!
             if (lerror_norm .or. locdensity_error .or. locchemspec_error) then
               call find_near_indeces(lower_i,upper_i,lower_j,upper_j,&
                                      lower_k,upper_k,x,y,z,fp_location)
               inear_error=(/lower_i,lower_j,lower_k/)
             endif
!
             if (it==nt) then
               if (lerror_norm) then
                 call safe_character_assign(file2,trim(datadir)//'/error_norm.dat'//chproc)
                 open(unit=89,file=file2,position='APPEND')

                 if(.not. linear_interpolate(f,iux,iuz,fp_location,vel_error,inear_error,.false.))&
                   call fatal_error('linear_interpolate','')
                 if(.not. linear_interpolate(f,ilnTT,ilnTT,fp_location,t_error,inear_error,.false.))&
                   call fatal_error('linear_interpolate','')
                 write(solid_cell_error,"(1I8)") it-1
                 if (ltemperature_nolog) then
                   write(error_string,94) theta,vel_error(1),vel_error(2),t_error
                 else
                   write(error_string,94) theta,vel_error(1),vel_error(2),exp(t_error)
                 endif
                 call safe_character_append(solid_cell_error,error_string)
                 write(89,*) trim(solid_cell_error)
                 close(89)
               endif
             endif ! finalize if "it==nt"
!
! Output Local density
             if (it==nt) then
               if(locdensity_error) then
                 open(unit=86,file='data/locdensity_error.dat'//chproc,position='APPEND')
                 if(.not. linear_interpolate(f,ilnrho,ilnrho,fp_location,rho_error,inear_error,.false.))&
                   call fatal_error('linear_interpolate','')
                 if (ldensity_nolog) then
                   write(86,*) theta, rho_error
                 else
                   write(86,*) theta, exp(rho_error)
                 endif
                 close(86)
               endif
             endif ! Finalize if "(it==nt)"
!
! Output Local species mass fraction
             if (it==nt) then
               if(locchemspec_error) then
                 open(unit=88,file='data/locchemspec_error.dat'//chproc,position='APPEND')
                 if(.not. linear_interpolate(f,ichemspec(1),ichemspec(nchemspec),&
                    fp_location,chemspec_error,inear_error,.false.))&
                   call fatal_error('linear_interpolate','')
                 write(88,'(1X,5e16.8)') theta, chemspec_error(:)
                 close(88)
               endif
             endif ! Finalize if "(it==nt)"
!
! Output Local reaction rate
             if (it==nt) then
               if (loutput_local_reaction_rate) then
                 fpchem(1:nchemspec)=f(ix0,iy0,iz0,ichemspec(1):ichemspec(nchemspec))*loc_density
                 call calc_reaction_rate(f,iobj,fpchem,reac_rate,char_reac)
!
                 call safe_character_assign(file3,trim(datadir)//'/local_reaction.dat'//chproc)
                 open(unit=91,file=file3,position='APPEND')
!                write(solid_cell_reac,"(1I8)") it-1
!                write(reac_string,"(5F12.7)") theta,reac_rate(ichemsO2),reac_rate(ichemsCO),&
!                 reac_rate(ichemsCO2),char_reac
                 write(91,"(5F12.7)") theta,reac_rate(ichemsO2),reac_rate(ichemsCO),&
                  reac_rate(ichemsCO2),char_reac
!                call safe_character_append(solid_cell_reac,reac_string)
!                write(91,*) trim(solid_cell_reac)
                 close(91)
               endif ! Finalize if "loutput_local_reaction_rate"
             endif ! Finalize if "it==nt"
!
98  format(1I8,1F15.8)
96  format(3F15.8)
94  format(4F15.8)
!
           endif
         endif
       enddo
     enddo
!
     call keep_compiler_quiet(df,f)
!   
   end subroutine output_solid_cells
!************************************************************************
    subroutine close_interpolation(f)
!
!  Dummy.
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    end subroutine close_interpolation
!************************************************************************
   subroutine update_solid_cells_pencil(f)
!
!  Dummy
!
     real, dimension (mx,my,mz,mfarray) :: f
!
     call keep_compiler_quiet(f)
!
   endsubroutine update_solid_cells_pencil
!***********************************************************************
  subroutine ib_nearest_grid(ibx,iby,ibz,xobj,yobj,zobj,robj,ibnearestgrid)
!
!  Find coordinates for nearest grid point of a given IB point.
!
!  April-2016/Chaoli: Adapted from fp_nearest_grid.
!
    real, intent(in)  :: robj, xobj, yobj, zobj
    real, intent(in)  :: ibx, iby, ibz
    integer           :: ipoint, inearest, icoord(8,3)
    integer           :: ixl, iyl, izl, ixu, iyu, izu, ju, jl, jm
    real              :: dx1, dy1, dz1
    real              :: dist_to_ib2(8), dist_to_cent2(8)
    logical           :: interiorpoint
    integer, dimension(3), intent(out)  :: ibnearestgrid
!
    dx1=1.0/dx
    dy1=1.0/dy
    dz1=1.0/dz
!
    interiorpoint = .true.
!
!  Find nearest grid point in x-direction
!
      if (nxgrid/=1) then
        if (ibx >= x(l1-1) .and. ibx <= x(l2+1)) then
          if (lequidist(1)) then
            ixl = int((ibx-x(1))*dx1) + 1
            ixu = ixl+1
          else
!
!  Find nearest grid point by bisection if grid is not equidistant
!
            ju=l2+1; jl=l1-1
            do while((ju-jl)>1)
              jm=(ju+jl)/2
              if (ibx > x(jm)) then
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
        if (iby >= y(m1-1) .and. iby <= y(m2+1)) then
          if (lequidist(2)) then
            iyl = int((iby-y(1))*dy1) + 1
            iyu = iyl+1
          else
!
!  Find nearest grid point by bisection if grid is not equidistant
!
            ju=m2; jl=m1
            do while((ju-jl)>1)
              jm=(ju+jl)/2
              if (iby > y(jm)) then
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
        if (ibz >= z(n1-1) .and. ibz <= z(n2+1)) then
          if (lequidist(3)) then
            izl = int((ibz-z(1))*dz1) + 1
            izu = izl+1
          else
!
!  Find nearest grid point by bisection if grid is not equidistant
!
            ju=n2; jl=n1
            do while((ju-jl)>1)
              jm=(ju+jl)/2
              if (ibz > z(jm)) then
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
!  i.e. the eight corners of the grid cell containing the IB point (ib).
!  Decide which ones are outside the object, and which one of these
!  is the closest one to ib:
!
!  Check if ib is within this processor's local domain
        if (interiorpoint) then
          dist_to_ib2(1) = (x(ixl)-ibx)**2+(y(iyl)-iby)**2+(z(izl)-ibz)**2
          dist_to_ib2(2) = (x(ixu)-ibx)**2+(y(iyl)-iby)**2+(z(izl)-ibz)**2
          dist_to_ib2(3) = (x(ixu)-ibx)**2+(y(iyu)-iby)**2+(z(izl)-ibz)**2
          dist_to_ib2(4) = (x(ixl)-ibx)**2+(y(iyu)-iby)**2+(z(izl)-ibz)**2
          dist_to_ib2(5) = (x(ixl)-ibx)**2+(y(iyl)-iby)**2+(z(izu)-ibz)**2
          dist_to_ib2(6) = (x(ixu)-ibx)**2+(y(iyl)-iby)**2+(z(izu)-ibz)**2
          dist_to_ib2(7) = (x(ixu)-ibx)**2+(y(iyu)-iby)**2+(z(izu)-ibz)**2
          dist_to_ib2(8) = (x(ixl)-ibx)**2+(y(iyu)-iby)**2+(z(izu)-ibz)**2
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
            if (dist_to_cent2(ipoint) > robj**2 .and. inearest == 0) then
              inearest=ipoint
            else if (dist_to_cent2(ipoint) > robj**2) then
              if (dist_to_ib2(ipoint) <= dist_to_ib2(inearest)) then
                inearest=ipoint
              endif
            endif
          enddo
!
!  Coordinates of nearest grid point. Zero if outside local domain.
          if (inearest > 0) then
            ibnearestgrid(:) = icoord(inearest,:)
          else
            print*, "WARNING: Could not find ibnearestgrid!"
          endif
!
        else ! ib is outside local domain and ibnearestgrid shouldn't exist
          ibnearestgrid(:) = 0
        endif
!
  endsubroutine ib_nearest_grid
!***********************************************************************
  subroutine calc_Diff_ib(rho_ib,T_ib,chem_ib, Diff_ib_ks)
!
!   August-2016/Chaoli: Adapted from calc_pencils_chemistry in module chemistry.
!

    real, intent(in)  :: rho_ib,T_ib
    real, dimension(nchemspec), intent(in)  :: chem_ib
    real, dimension(nchemspec), intent(out) :: Diff_ib_ks
    real  :: lambda_ib
    real, dimension(nchemspec,nchemspec)  :: Bin_Diff_coef_ib
    real  :: mu1_full_ib
    real,dimension(nchemspec)  :: XX_ib
    real :: tmp_sum, tmp_sum2
    integer  :: k, j

    call get_mu1_full_ib(chem_ib,mu1_full_ib)
    call get_XX_ib(chem_ib,mu1_full_ib,XX_ib)
    call calc_diff_visc_coef_ib(T_ib,rho_ib,mu1_full_ib,Bin_Diff_coef_ib)

    do k = 1,nchemspec
      tmp_sum = 0.
      tmp_sum2 = 0.
      do j = 1,nchemspec
        if (Mspecies(k) > 0.) then
          if (j /= k) then
            tmp_sum = tmp_sum &
                      +XX_ib(j)/Bin_Diff_coef_ib(j,k)
            tmp_sum2 = tmp_sum2 &
                      +XX_ib(j)*Mspecies(j)
!
          endif
        endif
      enddo
      Diff_ib_ks(k) = mu1_full_ib*tmp_sum2/tmp_sum
    enddo

  endsubroutine calc_Diff_ib
!***********************************************************************
    subroutine calc_diff_visc_coef_ib(T_ib,rho_ib,mu1_full_ib,Bin_Diff_coef_ib)
!
!  August-2016/Chaoli: Adapted from calc_pencils_chemistry in module chemistry.
!
!  Calculation of the binary diffusion coefficients and the species viscosities.
!
      use EquationOfState, only: tran_data

      real, intent(in)  :: T_ib,rho_ib,mu1_full_ib
      real :: Omega_kl, prefactor
      real :: lnTjk
      integer :: k, j
      real :: eps_jk, sigma_jk, m_jk, delta_jk, delta_st
      real :: Na=6.022E23, tmp_local, delta_jk_star
      real, dimension(nchemspec,nchemspec), intent(out) :: Bin_Diff_coef_ib
!
!  Find binary diffusion coefficients
!
      tmp_local = 3./16.*sqrt(2.*k_B_cgs**3/pi)

      prefactor = tmp_local*sqrt(T_ib) &
              *unit_length**3/(Rgas_unit_sys*rho_ib)
!
!  Do non-simplified binary diffusion coefficient
!
      do k = 1,nchemspec
        do j = k,nchemspec
!  Account for the difference between eq. 5-4 and 5-31 in the Chemkin theory
!  manual
!
           if (j /= k) then
             eps_jk = sqrt(tran_data(j,2)*tran_data(k,2))
             sigma_jk = 0.5*(tran_data(j,3)+tran_data(k,3))*1e-8
             m_jk = (Mspecies(j)*Mspecies(k)) &
                    /(Mspecies(j)+Mspecies(k))/Na
             delta_jk = 0.5*tran_data(j,4)*tran_data(k,4)*1e-18*1e-18
           else
             eps_jk = tran_data(j,2)
             sigma_jk = tran_data(j,3)*1e-8
             m_jk = Mspecies(j)/(2*Na)
             delta_jk = 0.5*(tran_data(j,4)*1e-18)*(tran_data(j,4)*1e-18)
           endif
!
           lnTjk = log(T_ib/eps_jk)
!
           Omega_kl = &
                 1./(6.96945701E-1+3.39628861E-1*lnTjk &
                 +1.32575555E-2*lnTjk**2 &
                 -3.41509659E-2*lnTjk**3 &
                 +7.71359429E-3*lnTjk**4 &
                 +6.16106168E-4*lnTjk**5 &
                 -3.27101257E-4*lnTjk**6 &
                 +2.51567029E-5*lnTjk**7)
           delta_jk_star = delta_jk/(eps_jk*k_B_cgs*sigma_jk**3)
!
           Omega_kl = Omega_kl &
                      +0.19*delta_jk_star*delta_jk_star/(T_ib/eps_jk)
            if (j /= k) then
               Bin_Diff_coef_ib(k,j) = prefactor/mu1_full_ib &
                        /(sqrt(m_jk)*sigma_jk*sigma_jk*Omega_kl)
            else
               Bin_Diff_coef_ib(k,j) = prefactor &
                        /(sqrt(m_jk)*sigma_jk*sigma_jk*Omega_kl)*Mspecies(k)
!
            endif

          enddo
        enddo
!

        do k = 1,nchemspec
          do j = 1,k-1
            Bin_Diff_coef_ib(k,j) = Bin_Diff_coef_ib(j,k)
          enddo
        enddo
!
    endsubroutine calc_diff_visc_coef_ib
!***********************************************************************
   subroutine get_mu1_full_ib(chem_ib,mu1_full_tmp)
!
!  Calculate mean molecular weight
!
      real, intent(in), dimension (nchemspec) :: chem_ib
      real, intent(out) :: mu1_full_tmp
      integer :: k
!
!  Mean molecular weight
!
      mu1_full_tmp=0.
      do k=1,nchemspec
        if (Mspecies(k)>0.) then
          mu1_full_tmp= &
          mu1_full_tmp+unit_mass*chem_ib(k) &
                  /Mspecies(k)
        endif
      enddo
!
    endsubroutine get_mu1_full_ib
!***********************************************************************
   subroutine get_XX_ib(chem_ib,mu1_full_ib,XX_ib)
!
!  Calculate mean molecular weight
!
      real, intent(in), dimension (nchemspec) :: chem_ib
      real, intent(in) :: mu1_full_ib
      real, intent(out), dimension (nchemspec) :: XX_ib
!
      XX_ib(:) = chem_ib(:)*unit_mass &
                      /(Mspecies(:)*mu1_full_ib)
!
    endsubroutine get_XX_ib
!***********************************************************************
endmodule Solid_Cells
