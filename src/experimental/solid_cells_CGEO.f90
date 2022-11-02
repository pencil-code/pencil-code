! $Id: solid_cells_stl.f90 2018-05-03 mao_chaoli@zju.edu.cn $
! $Id: zjulk@zju.edu.cn, jintai@zju.edu.cn $
!  This module add solid (as in no-fluid) cells in the domain.
!  This can be used e.g. in order to simulate a cylinder in a cross flow.
!  Now, apart from the circle geometry, a 2D geometry file similar to STL
!  can be used to represent the solid cells.Details can be found in the paper 
!  "A ghost-cell immersed boundary method for the simulations of heat
!  transfer in compressible flows under different boundary conditions
!  Part-II: Complex geometries"
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
!
      implicit none
!
      include 'solid_cells_stl.h'
!
      integer, parameter  :: maxterm=5 ! The max number of solid cells
      integer             :: nobjects
      character(len=32), dimension(maxterm) ::geometryfile='1stl' 
      integer, parameter :: DIM_STL = 10000 ! increased as needed.
      real, dimension(2,2,DIM_STL,maxterm) :: vertex
      real, dimension(3,DIM_STL,maxterm)   :: norm_line
      real, dimension(2,DIM_STL,maxterm)   :: centr_line 
      real, dimension(3,DIM_STL,maxterm)   :: loc_coor_utvec
      real, dimension(DIM_STL,maxterm)     :: lineelement
      real, dimension(maxterm)             :: sum_line
      integer, dimension(maxterm)          :: n_lines 
!  Reverse_normal_value is used to reverse the normal vector as needed,
      integer :: reverse_normal_value=1.0
!  TOL_STL is used to ignore very small line,
      real    :: TOL_STL=1e-23
!  SCALE_STL is used for transformation of unit_system.	
!  TX_STL... are used for transmation of all the line vertex as a whole.
      real    :: SCALE_STL=1.0, TX_STL=0.0, TY_STL=0.0
!
      real, dimension(maxterm)    :: XMIN_STL, XMAX_STL
      real, dimension(maxterm)    :: YMIN_STL, YMAX_STL
!
      integer, dimension(mx,my,mz,4):: ba
      integer, allocatable :: fpnearestgrid(:,:,:)
      real, allocatable    :: c_dragx(:), c_dragy(:), c_dragz(:), Nusselt(:)
      real, allocatable    :: c_dragx_p(:), c_dragy_p(:), c_dragz_p(:)
!
!  Diagnostic variables (need to be consistent with reset list below).
      integer :: idiag_c_dragx=0       ! DIAG_DOC:
      integer :: idiag_c_dragy=0       ! DIAG_DOC:
      integer :: idiag_c_dragz=0       ! DIAG_DOC:
      integer :: idiag_c_dragx_p=0     ! DIAG_DOC:
      integer :: idiag_c_dragy_p=0     ! DIAG_DOC:
      integer :: idiag_c_dragz_p=0     ! DIAG_DOC:
      integer :: idiag_Nusselt=0       ! DIAG_DOC:
!
!     real, dimension(maxterm)    :: drag_norm
!     real, dimension(maxterm)    :: nusselt_norm
      real, dimension(maxterm)    :: charac_len  ! for local nusselt number
!
      integer :: flow_dir
      logical :: lset_flow_dir=.false.
      integer :: flow_dir_set
      real    :: rhosum
      integer :: irhocount
      real    :: T0
      character(len=12)  :: bc_thermo_type='Dirichlet'
      character(len=12)  :: bc_vel_type='no-slip'
      character(len=12)  :: interpolation_method='IDW'
      logical :: lLoc_Nusselt_output=.false., lLoc_Temp_output=.false.
      logical :: lLoc_density_output=.false.
      logical :: lLoc_c_press_output=.false.
      logical :: lerror_norm=.false.
! 
      real,dimension(3)  :: solid_vel=(/ 0.0,0.0,0.0/)
      real      :: solid_temp=330.0
      real      :: solid_flux=30.0
      real      :: solid_tempflux=333.0
      real      :: solid_xpos=0.0, solid_ypos=0.0, solid_zpos=0.0
      real      :: solid_theta=0.0, solid_phi=0.0
      real      :: coefa=1.0,coefb=1.0
  
      real      :: init_uu=0
      character (len=labellen), dimension(ninit) :: initsolid_cells='nothing'    
!
  namelist /solid_cells_init_pars/ &
       nobjects,reverse_normal_value,TOL_STL, &
       SCALE_STL,TX_STL,TY_STL,charac_len, &
       lset_flow_dir,flow_dir_set,bc_thermo_type, &
       bc_vel_type,init_uu,solid_xpos,solid_ypos,solid_zpos, &
       initsolid_cells,geometryfile
!
  namelist /solid_cells_run_pars/  &
       interpolation_method,lLoc_Nusselt_output,lLoc_Temp_output, &
       lLoc_density_output,lLoc_c_press_output,solid_vel, &
       solid_temp,solid_flux, &
       solid_tempflux,solid_theta,solid_phi,coefa,coefb, &
       lerror_norm
!
      contains
!***********************************************************************
      subroutine initialize_solid_cells(f)
!
!  Define the geometry of the solid object.
!
!  19-nov-2008/nils: coded
!  28-sep-2010/nils: added spheres
!  nov-2010/kragset: updated allocations related to drag calculations
!  20-sep-2015/chaoli: adapted to read stl geometry file
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer  :: iobj
      integer  :: i
!
!  Prepare the solid geometry and fill ba matrix
!
      call get_stl_data
      call find_solid_cell_boundary(f)

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
        do iobj=1, nobjects
          allocate(fpnearestgrid(nobjects,n_lines(iobj),3))
        enddo
        allocate(c_dragx(nobjects))
        allocate(c_dragy(nobjects))
        allocate(c_dragz(nobjects))
        allocate(c_dragx_p(nobjects))
        allocate(c_dragy_p(nobjects))
        allocate(c_dragz_p(nobjects))
        call fp_nearest_grid
        rhosum    = 0.0
        irhocount = 0
!
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
          call fatal_error('initialize_solid_cells',&
              'I was not able to determine the flow direction!')
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
      endsubroutine initialize_solid_cells
!***********************************************************************
      subroutine get_stl_data
!
!     13-09-2015//: Coded by Chaoli.
!     This subroutine is used to read line vertices and 
!     normal vectors from a STL file
!
      implicit none
!
      integer :: iobj  ! iobj is used for loop of the number of stl files.
      integer :: ignored_lines 
      logical :: lstl_file_exist,keep_reading,lignore_current_line
      real    :: v1x,v1y
      real    :: v2x,v2y
      real    :: x12,y12
      real    :: d12
      real    :: n1,n2,norm
      real    :: abs_trans
      character(len=32) ::test_char, buff_char 
      real    :: abstrans
      real, dimension(2)  :: vec_12
!
! Start reading each STL FILE
!
      ignored_lines = 0
!
!  Loop over all solid objects
!
      do iobj=1,nobjects
        if (lroot) print*, 'Reading geometry from '//trim(geometryfile(iobj))//' ...'
!
        inquire(file=trim(geometryfile(iobj)),EXIST=lstl_file_exist)
        if (.not.lstl_file_exist) then
          if (lroot) then
            print*, 'Stl file'// trim(geometryfile(iobj))//'missing'
          endif
          call fatal_error('get_stl_data','Stl file does not exist!')
        endif
!
!  Open geometry.stl ASCII FILE
!
        open(unit=333, file=trim(geometryfile(iobj)), FORM='FORMATTED', STATUS='OLD')
!
        if (lroot) print*,'STL file opened. Starting reading data...'
!
        keep_reading = .TRUE.
!
!  Initialize sum_line
        sum_line(iobj)=0.0
!  Initialize n_lines
        n_lines(iobj)=0
!
        do while(keep_reading)
!
          read(333,*) test_char
!         if (lroot) print *,'test_char=',test_char
          if (trim(test_char) == 'facet') then
!
            backspace(333)
            lignore_current_line = .FALSE.
!
            read(333,*) buff_char,buff_char,n1,n2 ! Read normal vector
            read(333,*)
!
!  Read the corordinates of the three vertexs of every line.
            read(333,*) buff_char, v1x,v1y
            read(333,*) buff_char, v2x,v2y
!
!  Reverse normal vector if needed (this will switch fluid and blocked cells)
            n1 = n1 * reverse_normal_value
            n2 = n2 * reverse_normal_value
!
            x12 = v2x - v1x
            y12 = v2y - v1y
            vec_12=(/x12,y12/)
!
            d12 = sqrt(x12**2+y12**2)
!
            if (d12>TOL_STL) then
              lignore_current_line=.false.
            else
              lignore_current_line=.true.
            endif
!
            norm = sqrt(n1**2+n2**2)
!
            if (norm>TOL_STL) then
              n1 = n1 /norm
              n2 = n2 /norm
            else
              lignore_current_line = .TRUE.  ! Ignore facets with zero normal vector
            endif
!
            if (lignore_current_line) then
              ignored_lines = ignored_lines + 1
            else  
!			   
!  Save vertex coordinates for valid facets and perform translation
              n_lines(iobj) = n_lines(iobj) + 1
!
              norm_line(1,n_lines(iobj),iobj) = n1
              norm_line(2,n_lines(iobj),iobj) = n2
!  For the slip velocity boundary condition 
              norm_line(3,n_lines(iobj),iobj) = 0

              loc_coor_utvec(1,n_lines(iobj),iobj)=x12/d12
              loc_coor_utvec(2,n_lines(iobj),iobj)=y12/d12
!  For the slip velocity boundary condition 			  
              loc_coor_utvec(3,n_lines(iobj),iobj)=0
!
              vertex(1,1,n_lines(iobj),iobj) = SCALE_STL*v1x + TX_STL
              vertex(1,2,n_lines(iobj),iobj) = SCALE_STL*v1y + TY_STL

              vertex(2,1,n_lines(iobj),iobj) = SCALE_STL*v2x + TX_STL
              vertex(2,2,n_lines(iobj),iobj) = SCALE_STL*v2y + TY_STL

              centr_line(1,n_lines(iobj),iobj) = (vertex(1,1,n_lines(iobj),iobj) + &
                                                  vertex(2,1,n_lines(iobj),iobj))/2.0
              centr_line(2,n_lines(iobj),iobj) = (vertex(1,2,n_lines(iobj),iobj) + &
                                                  vertex(2,2,n_lines(iobj),iobj))/2.0

!  Calculate the length of the line element and the total length of the solid  									  
              lineelement(n_lines(iobj),iobj)=d12
              sum_line(iobj)=sum_line(iobj)+lineelement(n_lines(iobj),iobj)
!
            endif
!
          elseif (trim(test_char) == 'endsolid') then
            keep_reading = .FALSE.
          endif
!
        enddo ! finish reading
!		 
       if (lroot)  print*,'Done reading file.'
!
        close(333)
!		 
!  To obtain the coordinates of the northest, southest, westest, eastest, 
!  topest and bottomest of the stl geometry, for the sake of logical determination. 
        XMIN_STL(iobj) = minval(vertex(:,1,1:n_lines(iobj),iobj))
        XMAX_STL(iobj) = maxval(vertex(:,1,1:n_lines(iobj),iobj))
        YMIN_STL(iobj) = minval(vertex(:,2,1:n_lines(iobj),iobj))
        YMAX_STL(iobj) = maxval(vertex(:,2,1:n_lines(iobj),iobj))

        if (lroot) then
          print*,'STL file(s) successfully read.'
          print*,'STL file number',iobj
          print*,'Total number of facets read =',n_lines(iobj) + ignored_lines
          print*,'Number of valid facets      =',n_lines(iobj)
          print*,'Number of ignored facets    =',ignored_lines
          print*,'Geometry range from stl file:'
          if (SCALE_STL/=1.0) then
            print*,'Scaling factor:',SCALE_STL
          endif
          abstrans = abs(TX_STL)+abs(TY_STL)
          if (abstrans>TOL_STL) then
            print*,'Translation vector (X,Y,Z)=',TX_STL,TY_STL
          endif
          print*,'x-range = ', XMIN_STL(iobj),XMAX_STL(iobj)
          print*,'y-range = ', YMIN_STL(iobj),YMAX_STL(iobj)
          print*,' '
        endif
      enddo  ! finalize the loop of nobjects

      return

      end subroutine get_stl_data
!***********************************************************************
      subroutine find_solid_cell_boundary(f)
!
!  Find the boundaries of the geometries such that we can set the
!  ghost points inside the solid geometry in order to achieve the
!  correct no-slip boundaries.
!
!  Store data in the ba array.
!  If ba(ip,jp,kp,1)= 0 we are in a fluid cell (i.e. iline not inside a solid geometry)
!  If ba(ip,jp,kp,1)=10 we are in a fluid cell which are so close to the
!                       surface of the solid geometry that we must set the
!                       value of this point by some special method.(can't.)
!  If ba(ip,jp,kp,1)= 9 we are inside a solid geometry, but far from the boundary
!  If ba(ip,jp,kp,1)=-1 we are inside a solid geometry, and the point at ip+1
!                       is outside the geometry.
!  If ba(ip,jp,kp,1)=-3 we are inside a solid geometry, and the point at ip+3
!                       is outside the geometry.
!  If ba(ip,jp,kp,2)=-3 we are inside a solid geometry, and the point at jp+3
!                       is outside the geometry.
!  If ba(ip,jp,kp,2)=11 we are inside a solid geometry, either close to or far
!                       from the boundary, but the position (ip,jp,kp) is a ghost
!                       point at the current processor. (can't.)
!
!  The number stored in ba(ip,jp,kp,4) is the number of the object
!
!  19-nov-2008/nils: coded
!  14-09-2015/chaoli: adapted to include stl geometry files.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: iline
      integer :: i,j,k
      integer :: iobj
!
!  Distance between a given grid point and the centroid of the nearest line
      real    :: p_centr_dist
!
      real    :: dist_standard
!
!  Recorder of the index of the nearest line to the given grid point
      integer :: iline_nearest
!
!  Direction vector from the given grid point to centroid of nearest line
      real, dimension(3)    :: centr_p_v
      real                  :: deter_dot_product
!
!  Initialize ba
!
      ba=0
!
!  Loop over all grid points, including the ghost points for now
!
!  These files are used to output the geometry of solid, debug
!     open(222, file='data/solid_point', form='formatted')
!     open(223, file='data/fluid_point', form='formatted')
!     open(224, file='data/nearest_line_centr', form='formatted')
!
      do i=l1, l2
      do j=m1, m2
      do k=n1, n2
!
!  Initialize dist_standard
      dist_standard=impossible
! Loop over all objects
        do iobj=1, nobjects
!
          if ((x(i)>=(XMIN_STL(iobj)-2.0*dxmin)) .and. (x(i)<=(XMAX_STL(iobj)+2.0*dxmin)) &
                                                                            .and. &
              (y(j)>=(YMIN_STL(iobj)-2.0*dxmin)) .and. (y(j)<=(YMAX_STL(iobj)+2.0*dxmin))) then
!
!  Loop over all facets
!
!  These files are used for the sake of debug.
!           open(220, file='data/line_centr', form='formatted')
!           open(221, file='data/line_normal', form='formatted')
!
            do iline=1, n_lines(iobj)
!             if (it==1) write(220,*) centr_line(1,iline,iobj),centr_line(2,iline,iobj)
!             if (it==1) write(221,*) norm_line(1,iline,iobj),norm_line(2,iline,iobj)
              p_centr_dist=sqrt((x(i)-centr_line(1,iline,iobj))**2 + &
                                (y(j)-centr_line(2,iline,iobj))**2)
!
!  Find the closest line to a given grid point
              if (p_centr_dist<dist_standard) then
                iline_nearest=iline
                dist_standard=p_centr_dist
              endif
!
            enddo ! end loop of lines
!
! For the sake of debug.
!           close(220)
!           close(221)
!
!  To determine whether the grid point is inside stl geometry or mot.
!  Using the dot product between the direction vector from the centroid of 
!  the nearest line to the grid point and normal vector of the nearest line.
!
            centr_p_v(1)=x(i)-centr_line(1,iline_nearest,iobj)
            centr_p_v(2)=y(j)-centr_line(2,iline_nearest,iobj)
            deter_dot_product=centr_p_v(1)*norm_line(1,iline_nearest,iobj) + &
                              centr_p_v(2)*norm_line(2,iline_nearest,iobj)
!
            if (deter_dot_product<0) then ! inside solid stl geometry
              ba(i,j,k,1:3)=9
              ba(i,j,k,4)=iobj
!             write(222, *) x(i),y(j)
            else
              ba(i,j,k,1:3)=0
!             if (it==1) write(223, *) x(i),y(j)
!             if (it==1) write(224, *) centr_line(1,iline_nearest,iobj),centr_line(2,iline_nearest,iobj)
            endif
!
!  End the if to determin whether a grid point is within range
!  (stl_xmin,stl_xmax),(stl_ymin,stl_ymax) and (stl_zmin,stl_zmax)
          endif
!
!  Finalize loop over all solid cells
        enddo 
!		
      enddo
      enddo
      enddo  ! end loop of grid points on the local processor.
!
!  For the sake of debug
!     close(222)
!     close(223)
!     close(224)
!
!  Find the ghost points which are those grid points that has at least 
!  one neighbor fluid point in the east/west, north/south or top/bottom 
!  direction.
!
!  Loop all the grid points on the local processor to mark ghost points
!
! This file coatains coordinates of ghost points, for debug
      open(225, file='data/ghost_point', form='formatted')
!
      do i=l1, l2
      do j=m1, m2
      do k=n1, n2
!
        if (ba(i,j,k,1)==9 .and. ba(i,j,k,2)==9 .and. &
            ba(i,j,k,3)==9) then
!
!  First look x direction 
!
!         x+ direction
          if (ba(i+1,j,k,1)==0) then
            ba(i,j,k,1)=-1  ! the first layer of ghost points.
          elseif (ba(i+2,j,k,1)==0) then
            ba(i,j,k,1)=-2  ! the second layer of ghost points.
          elseif (ba(i+3,j,k,1)==0) then
            ba(i,j,k,1)=-3  ! the third layer of ghost points.
          endif
!
!         x- direction		
          if (ba(i-1,j,k,1)==0) then
            ba(i,j,k,1)=1  ! the first layer of ghost points.
          elseif (ba(i-2,j,k,1)==0) then
            ba(i,j,k,1)=2  ! the second layer of ghost points.
          elseif (ba(i-3,j,k,1)==0) then
            ba(i,j,k,1)=3  ! the third layer of ghost points.
          endif
!
!  Secondly look y direction
!
!         y+ direction
          if (ba(i,j+1,k,2)==0) then
            ba(i,j,k,2)=-1  ! the first layer of ghost points.
          elseif (ba(i,j+2,k,2)==0) then
            ba(i,j,k,2)=-2  ! the second layer of ghost points.
          elseif (ba(i,j+3,k,2)==0) then
            ba(i,j,k,2)=-3  ! the third layer of ghost points.
          endif
!
!          y- direction
          if (ba(i,j-1,k,2)==0) then
            ba(i,j,k,2)=1  ! the first layer of ghost points.
          elseif (ba(i,j-2,k,2)==0) then
            ba(i,j,k,2)=2  ! the second layer of ghost points.
          elseif (ba(i,j-3,k,2)==0) then
            ba(i,j,k,2)=3  ! the third layer of ghost points.
          endif
!
!  End if that the determanation whether a grid point is within 
!  solid boundary
        endif   
!
! for debug
        if (((ba(i,j,k,1)/=0) .and. (ba(i,j,k,1)/=9)) .or. &
           ((ba(i,j,k,2)/=0) .and. (ba(i,j,k,2)/=9))) then
          if (it==1) write(225, *) x(i),y(j)
        endif
!
      enddo
      enddo
      enddo ! end loop of grid points on the local processor.
!
!  For the sake of debug
      close(225)
!
!  Set zero value of all variables inside the solid geometry far from
!  all interfaces. This is done for more easy interpretation of postprocessing.
!
      if (.not.lreloading) then
          do i=l1,l2
          do j=m1,m2
          do k=n1,n2
            if (ba(i,j,k,1)==9 .and. ba(i,j,k,2)==9) then
              f(i,j,k,iux:iuz)=0.0
!
!  Set the temperature of the object to points inside solid geometry
!  far from all interfaces. By Mcl.
!
              if (ilnTT/=0) then
                if (ltemperature_nolog) then
                  f(i,j,k, iTT)=solid_temp
                else
                  f(i,j,k, ilnTT)=log(solid_temp)
                endif
              endif
!
! This is for debug.
!           elseif (ba(i,j,k,1)/=0 .or. ba(i,j,k,2)/=0) then
!             f(i,j,k,iux:iuz)=0.2
!
            endif
          enddo
          enddo
          enddo
      endif
!
      end subroutine find_solid_cell_boundary
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
      real :: wall_smoothing
      real :: wall_smoothing_temp
      integer :: jj
!
      do jj=1,ninit
        select case (initsolid_cells(jj))
!
!  This overrides any initial conditions set in the Hydro module.
!
        case ('nothing')
          if (lroot) print*,'init_solid_cells: nothing'
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
!  "forcepoints" (fp) for solid object. Here, we use the centroid of the line
!  as forceppoint.
!
!  mar-2009/kragset: coded
!  nov-2010/kragset: updated to include spheres
!  sep-2015/chaoli: adapted to be used for stl geometry
!
      integer           :: iline
      integer           :: ipoint, inearest, icoord(8,3)
      integer           :: ixl, iyl, izl, ixu, iyu, izu, ju, jl, jm
      real              :: fpx, fpy, fpz
      real              :: dx1, dy1, dz1
      real              :: dist_to_fp2(8)
      logical           :: interiorpoint
      integer           :: i, j, k
      logical           :: lfluid_point
      integer           :: iobj
!
      integer              :: xindex, yindex
!
      dx1=1/dx
      dy1=1/dy
      dz1=1/dz
!
!  Loop over all forcepoints on each object, iobj
!  Loop over all objects
!
!     open(444, file='data/fp_nearest', form='formatted')
!
      do iobj=1, nobjects
      do iline=1,n_lines(iobj)
!
!  Marking whether fp is within this processor's domain or not
        interiorpoint = .true.
!
!  Find fp coordinates. We use the centroid of the line as forcepoint.
        fpx = centr_line(1,iline,iobj)
        fpy = centr_line(2,iline,iobj)
        fpz = z(n1)
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
          dist_to_fp2(3) = (x(ixl)-fpx)**2+(y(iyu)-fpy)**2+(z(izl)-fpz)**2
          dist_to_fp2(4) = (x(ixu)-fpx)**2+(y(iyu)-fpy)**2+(z(izl)-fpz)**2
          dist_to_fp2(5) = (x(ixl)-fpx)**2+(y(iyl)-fpy)**2+(z(izu)-fpz)**2
          dist_to_fp2(6) = (x(ixu)-fpx)**2+(y(iyl)-fpy)**2+(z(izu)-fpz)**2
          dist_to_fp2(7) = (x(ixl)-fpx)**2+(y(iyu)-fpy)**2+(z(izu)-fpz)**2
          dist_to_fp2(8) = (x(ixu)-fpx)**2+(y(iyu)-fpy)**2+(z(izu)-fpz)**2
          icoord(1,:) = (/ixl,iyl,izl/)
          icoord(2,:) = (/ixu,iyl,izl/)
          icoord(3,:) = (/ixl,iyu,izl/)
          icoord(4,:) = (/ixu,iyu,izl/)
          icoord(5,:) = (/ixl,iyl,izu/)
          icoord(6,:) = (/ixu,iyl,izu/)
          icoord(7,:) = (/ixl,iyu,izu/)
          icoord(8,:) = (/ixu,iyu,izu/)
          inearest=0
          ipoint=0
          do k=izl, izu
          do j=iyl, iyu
          do i=ixl, ixu
!  Test if we are in a fluid cell.
            lfluid_point=((ba(i,j,k,1)==0).and. &
                          (ba(i,j,k,2)==0).and.(ba(i,j,k,3)==0))
            ipoint=ipoint+1
            if (lfluid_point.and. inearest == 0) then
              inearest=ipoint
            else if (lfluid_point) then
              if (dist_to_fp2(ipoint) <= dist_to_fp2(inearest)) then
                inearest=ipoint
              endif
            endif
          enddo
          enddo
          enddo
!
!  Coordinates of nearest grid point. Zero if outside local domain.
          if (inearest > 0) then
            fpnearestgrid(iobj,iline,:)=icoord(inearest,:)
          else
            if (lroot) print*, "WARNING: Could not find fpnearestgrid!"
          endif
! 
! for debug
          xindex=fpnearestgrid(iobj,iline,1)
          yindex=fpnearestgrid(iobj,iline,2)
!         write(444,*) x(xindex),y(yindex)
!
        else ! fp is outside local domain and fpnearestgrid shouldn't exist
          fpnearestgrid(iobj,iline,:) = 0
        endif
      enddo ! Finalize loop over all lines
      enddo ! Finalize loop over all objects
!     close(444)
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
!  sep-2015/chaoli: adapted to include stl geometry
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
!  sep-2015/chaoli: adapted to include stl geometry
!
      use General, only: safe_character_append
      use Mpicomm, only: mpireduce_sum, mpireduce_sum_int, mpibcast_real
!
      real    :: rhosum_all, c_dragx_all(nobjects), c_dragy_all(nobjects)
      real    :: c_dragz_all(nobjects), Nusselt_all(nobjects)
      real    :: c_dragx_p_all(nobjects), c_dragy_p_all(nobjects)
      real    :: c_dragz_p_all(nobjects)
      integer :: irhocount_all
      real    :: norm, refrho0
      character(len=100) :: numberstring
      character(len=500) :: solid_cell_drag
      integer            :: iobj
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
              open(unit=81,file='data/dragcoeffs.dat',position='APPEND')
              write(solid_cell_drag,84) it-1, t
              do iobj=1,nobjects
                write(numberstring,82) c_dragx(iobj),c_dragx_p(iobj), c_dragy(iobj)
                call safe_character_append(solid_cell_drag,numberstring)
              enddo
              write(81,*) trim(solid_cell_drag)
              close(81)
84            format(1I8,1F15.8)
82            format(4F15.8)
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
!  sep-2015/chaoli: Adapted to be used for stl geometry
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mvar) :: f_tmp
      integer :: i,j,k
      real :: z_obj,y_obj,x_obj,r_obj,r_new,r_point,dr
      real :: xghost, yghost, zghost
      logical :: bax, bay, baz
      real, dimension(3) :: xxp
      real :: xib, yib, zib
      real :: xmirror, ymirror, zmirror
      integer  :: ix0, iy0, iz0
      integer  :: line_index
      real     :: dist_mir_gp
      integer  :: iobj
      integer  :: lower_i, upper_i
      integer  :: lower_j, upper_j
      integer  :: lower_k, upper_k
      logical  :: lfind_proj_point=.false.
!
!  For debug
!     real, dimension(3)  :: g_ib
!     real, dimension(3)  :: line_v
!     real                :: check_dot
!
!  For debug
!     real, dimension(3)  :: g_mir
!     real, dimension(3)  :: line_v_mir
!     real                :: check_dot_mir
!     real                :: dist_g_origin
!     real                :: xmirror_a, ymirror_a, zmirror_a
!     real                :: check_dot_a
!
!  Loop over all grid points belonging to local processor
!
!  These two files are used for debug
       open(226, file='data/nearest_point', form='formatted')
       open(227, file='data/proj_point', form='formatted')
!      open(229, file='data/mir_point', form='formatted')
!      open(230, file='data/interpolate_message', form='formatted')
!      open(231, file='data/surroundgrid_point', form='formatted')
      do i=l1,l2
      do j=m1,m2
      do k=n1,n2
!
        bax=((ba(i,j,k,1)/=0).and.(ba(i,j,k,1)/=9))
        bay=((ba(i,j,k,2)/=0).and.(ba(i,j,k,2)/=9))
!
!  Check if we are in a point which must be interpolated, i.e. we are inside
!  a solid geometry AND we are not more than three grid points from the
!  closest solid-fluid interface
!
        if (bax .or. bay) then
!  Find which object this ghost point belongs to
          iobj=ba(i,j,k,4)
          if (iobj==0) call fatal_error('update_solid_cells','iobj=0')
!  Find x, y and z values of ghost point
          xghost=x(i); yghost=y(j); zghost=z(k)
!
! For debug
!         dist_g_origin=sqrt((xghost)**2+(yghost)**2+(zghost)**2)
!         xmirror_a=(0.05-dist_g_origin)/dist_g_origin*xghost
!         ymirror_a=(0.05-dist_g_origin)/dist_g_origin*yghost
!         zmirror_a=z(k)
! 
!  First we find the coordinates of projective point(i.e.,IB point) of 
!  ghost point on the nearest line
          call find_IB_point(xghost,yghost,zghost,xib,yib,zib, &
                             line_index,iobj,lfind_proj_point)
!
!  For debug
!         g_ib(1)=xghost-xib;g_ib(2)=yghost-yib;g_ib(3)=zghost-zib
!         line_v(1)=vertex(2,1,line_index,iobj)-vertex(1,1,line_index,iobj)
!         line_v(2)=vertex(2,2,line_index,iobj)-vertex(1,2,line_index,iobj)
!         line_v(3)=vertex(2,3,line_index,iobj)-vertex(1,3,line_index,iobj)
!         check_dot=g_ib(1)*line_v(1)+g_ib(2)*line_v(2)+g_ib(3)*line_v(3)
!         if (it==1) then
!           open(unit=666,file='data/check_dot',position='APPEND')
!           write(666,*) check_dot
!           close(666)
!         endif
!
!  Output for debug
          if (it==1) then
            if (lfind_proj_point) then 
              write(227,*) xib,yib
            else
              write(226,*) xib,yib
            endif
          endif
!
!  Coordinates of the moirror points
          xmirror=2.0*xib-xghost
          ymirror=2.0*yib-yghost
          zmirror=z(k)
!
!  For debug
!         g_mir(1)=xghost-xmirror;g_mir(2)=yghost-ymirror;g_mir(3)=zghost-zmirror
!         line_v_mir(1)=vertex(2,1,line_index,iobj)-vertex(1,1,line_index,iobj)
!         line_v_mir(2)=vertex(2,2,line_index,iobj)-vertex(1,2,line_index,iobj)
!         line_v_mir(3)=vertex(2,3,line_index,iobj)-vertex(1,3,line_index,iobj)
!         check_dot_mir=g_mir(1)*line_v_mir(1)+g_mir(2)*line_v_mir(2)+g_mir(3)*line_v_mir(3)
!         check_dot_a=(xghost-xmirror_a)*line_v_mir(1)+(yghost-ymirror_a)*line_v_mir(2) &
!                     +(zghost-zmirror_a)*line_v_mir(3)
!         if (it==1) then
!           open(unit=667,file='data/check_dot_mir',position='APPEND')
!           write(667,*) check_dot_mir
!           write(667,*) check_dot_a
!           close(667)
!         endif
!
!         if (it==1) then
!           write(229,*) xmirror,ymirror
!         endif
!
!  Distance between mirror and ghost point
          dist_mir_gp=sqrt((xmirror-xghost)**2+ &
                           (ymirror-yghost)**2)
!
!  Check if mirror point lies within local processor domain
          if ((xmirror>=x(1).and.xmirror<=x(mx)).and. &
              (ymirror>=y(1).and.ymirror<=y(my)).and. &
	      (zmirror>=z(1).and.zmirror<=z(mz))) then
! Everything OK
          else 
            print*,'mirror point is outside local domain'
            print*, 'iproc_loc = ', iproc
            print*, 'mx, x(1), x(mx) = ', mx, x(1), x(mx)
            print*, 'my, y(1), y(my) = ', my, y(1), y(my)
            print*, 'mz, z(1), z(mz) = ', mz, z(1), z(mz)
            print*, 'xmirr=', xmirror
            print*, 'ymirr=', ymirror
            print*, 'zmirr=', zmirror
!
            call fatal_error('update_solid_cells','mp outside local domain')
          endif
!
!  Find i, j and k indeces for points to be used during interpolation
          xxp=(/xmirror, ymirror, zmirror /) 
          call find_near_indeces(lower_i,upper_i,lower_j,upper_j,lower_k, &
                                 upper_k,x,y,z,xxp)
!  
!  Check if mirror point really lies within the range (x(lower_i),x(upper_i)),
!  (y(lower_j),y(upper_j)) and (z(lower_k),z(upper_k)).
          ix0=lower_i; iy0=lower_j; iz0=lower_k	
!
!         if (it==1) then
!           write(231,*) x(ix0),y(iy0)
!           write(231,*) x(ix0+1),y(iy0)
!           write(231,*) x(ix0),y(iy0+1)
!           write(231,*) x(ix0+1),y(iy0+1)
!         endif
!
          if ((x(ix0)<=xxp(1) .and. x(ix0+1)>=xxp(1) &
		       .or. nxgrid==1) .and. &
              (y(iy0)<=xxp(2) .and. y(iy0+1)>=xxp(2) &
                       .or. nygrid==1) .and. &
              (z(iz0)<=xxp(3) .and. z(iz0+1)>=xxp(3) &
                       .or. nzgrid==1)) then
!  Everything okay
          else
            print*, 'update_solid_cells: Interpolation point does not ' // &
                    'lie within the calculated grid point interval.'
            print*, 'iproc = ', iproc
            print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
            print*, 'xp, xp0, xp1 = ', xxp(1), x(ix0), x(ix0+1)
            print*, 'yp, yp0, yp1 = ', xxp(2), y(iy0), y(iy0+1)
            print*, 'zp, zp0, zp1 = ', xxp(3), z(iz0), z(iz0+1)
!
            call fatal_error('update_solid_cells','')
          endif
!
!  First we use interpolations to find the value of the mirror point.
!  Then we use the interpolated value to find the value of the ghost point
!  by empoying either Dirichlet or Neuman boundary conditions.
!
!         open(230, file='data/interpolate_message', form='formatted')
          if (interpolation_method=='IDW') then
!           if (it==1) write(230,*) 'IDW interpolate method is used.'
            call interpolate_idw(f,f_tmp,lower_i,lower_j,lower_k, &
			         xmirror,ymirror,zmirror,line_index,dist_mir_gp,iobj)
            f(i,j,k,1:mvar)=f_tmp
          elseif (interpolation_method=='LINEAR') then
!           if (it==1) write(230,*) 'LINEAR interpolate method is used.'
            call interpolate_linear(f,f_tmp,lower_i,lower_j,lower_k,xmirror,ymirror, &
                                    zmirror,xib,yib,zib,line_index,iobj,dist_mir_gp)
            f(i,j,k,1:mvar)=f_tmp
          elseif (interpolation_method=='MATRIX') then
!           if (it==1) write(230,*) 'MATRIX interpolate method is used.'
            call interpolate_matrix(f,f_tmp,lower_i,lower_j,lower_k, &
			            xmirror,ymirror,zmirror,xghost,yghost,zghost, &
                                    line_index,dist_mir_gp,iobj)
            f(i,j,k,1:mvar)=f_tmp
          elseif (interpolation_method==' ') then
!           if (it==1) write(230,*) 'NO interpolate method is used.'
!           do nothing. debug to see if we can identify the right interface
          else
!           if (it==1) write(230,*) 'unknown interpolation_method is used.'
            call fatal_error('update_solid_cells','unknown interpolation_method')
          endif
!
! For debug 
!         if (it==1) then
!           open (229, file='data/interpolation_message', position='APPEND')
!           write(229,*) 'after interpolation_method'
!         endif
!         close(229)
!
        endif ! end if this is a ghost point
!
      enddo
      enddo
      enddo
!
!  For debug
      close(226)
      close(227)
!     close(229)
!     close(230)
!     close(231)
!
      endsubroutine update_solid_cells
!***********************************************************************
      subroutine update_solid_cells_pencil(f)
!
!  Set the boundary values of the solid area such that we get a
!  correct fluid-solid interface.
!
!  30-mar-15/Jorgen+nils: coded
!  17-jun-15/chaoli: do nothing
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      call keep_compiler_quiet(f)
!
      endsubroutine update_solid_cells_pencil
!***********************************************************************
      subroutine find_IB_point(xa,yb,zc,IB_x,IB_y,IB_z,line_index,iobj,lfind_proj_point)
!
!  For a given ghost point, find the according IB point. IB points are defined as
!  the normal intercepts between the normal vector starting from the given ghost point 
!  normal to solid-fluid interface and the solid boundary. 
!  Firstly, the nearest vertex to the ghost point is found. Secondly,  find the normal
!  intercept between the ghost point and the facets that share the nearest vertex.    
!
!  15-09-2015/chaoli: coded
!
      real, intent(in)    :: xa,yb,zc
      integer, intent(in) :: iobj
      integer :: iline
      integer, intent(out) :: line_index
      real    :: gp_v1_dist,gp_v2_dist,gp_v3_dist,gp_v_dist
!
!  Coordinates of the nearest vertex of the three/two ones of a line
      real, dimension(2)  :: nearest_of_three
      real    :: dist_standard
!
!  Coordinates of the nearest vertex to the given ghost point
      real, dimension(2)  :: nearest_vertex
!
      logical, dimension(n_lines(iobj)) :: lnearest_line
!
      integer             :: ii, jj
      integer             :: n_nearest_line
!
      real, dimension(2)  :: gp_v1_vec
      real                :: line_n_norm
      real                :: IB_x_temp,IB_y_temp,IB_z_temp
      real, intent(inout) :: IB_x,IB_y,IB_z
!
      real                :: gp_ib_dist
      logical             :: lon_line
      logical, optional   :: lfind_proj_point

      integer, parameter  :: N_MAX_FACET=20 ! should be modified as needed
!
!  A recorder to store the index iline of the nearest line
!  20 can be modified as needed
      integer, dimension(5)  :: nearest_line_index 
!
!  For debug
      real,    dimension(2)  :: v1_1, v1_2, vj_1, vj_2
      integer                :: iline_1, iline_j
      logical                :: lshare_1, lshare_2, lshare_3, lshare_4
!
!  Initialize dist_standard
      dist_standard=impossible
!
!  Initialize the logical matrix lnearest_line to be .true.
!  for all facets
      lnearest_line(:)=.true.
!
!  Loop over all the lines
      do iline=1, n_lines(iobj)
!
!  Find the distance between ghost point and the three/two verteices
        gp_v1_dist=sqrt((xa-vertex(1,1,iline,iobj))**2+ &
                        (yb-vertex(1,2,iline,iobj))**2)
        gp_v2_dist=sqrt((xa-vertex(2,1,iline,iobj))**2+ &
                        (yb-vertex(2,2,iline,iobj))**2)
! 
        gp_v_dist=min(gp_v1_dist, gp_v2_dist)
!
!  Determine which one of the two verteices is the nearest one
        if (gp_v_dist==gp_v1_dist) then
	  nearest_of_three(:)=vertex(1,1:2,iline,iobj)
        elseif (gp_v_dist==gp_v2_dist) then
	  nearest_of_three(:)=vertex(2,1:2,iline,iobj)
        endif
!
!  Find the nearest vertex to a given ghost
        if (gp_v_dist<dist_standard) then
          nearest_vertex(:)=nearest_of_three(1:2)
          dist_standard=gp_v_dist
!
!  For the situation where the nearest vertex is IB point
          line_index=iline
!
          if (iline>=2) then
            lnearest_line(1:iline-1)=.false.
          endif
!
        elseif (gp_v_dist>dist_standard) then
          lnearest_line(iline)=.false.
        else
!         do nothing
        endif ! end if "gp_v_dist<dist_standard"
!
!{  Initialize the number of nearest facets. 
!       n_nearest_line=0
!
!  The total number of facets that are nearest to a given point
!       do ii=1, iline
!
!         if (lnearest_line(ii)) then
!           n_nearest_line=n_nearest_line+1
!           nearest_line_index(n_nearest_line)=ii
!         endif
!
!       enddo
!
!  N_MAX_FACET is MAX number of facets around a common vertex.
!       if (n_nearest_line>=N_MAX_FACET) exit }(more efficient??)
!
      enddo ! end loop of lines
!
!  Assign the coordinates of the nearest vertex to IB point, temperorarily
!  If we fail to find a projection as IB point, then this is the situation.
      IB_x=nearest_vertex(1)
      IB_y=nearest_vertex(2)
      IB_z=zc
!
!  Initialize the number of nearest facets.
      n_nearest_line=0
!
!  The total number of facets surrounding the nearest vertex
!  that are closest to a given ghost point
      do ii=1, n_lines(iobj)
!
        if (lnearest_line(ii)) then
          n_nearest_line=n_nearest_line+1
          nearest_line_index(n_nearest_line)=ii
        endif
!
      enddo
!
!     if (lroot) print*,'number of nearest line element', n_nearest_line
!  Check if n_nearest_line equal to 0
      if (n_nearest_line==0) &
        call fatal_error ('find_IB_point','n_nearest_line=0')
!  Debug
      if (n_nearest_line/=2) &
        call fatal_error ('find_IB_point','n_nearest_line/=2')
!
!  Add a check here that facets with lnearest_line==.true. 
!  share a common vertex. Debug.
      do jj=1, n_nearest_line
        iline_j=nearest_line_index(jj)
        iline_1=nearest_line_index(1)        
        v1_1=vertex(1,1:2,iline_1,iobj)
        v1_2=vertex(2,1:2,iline_1,iobj)
        vj_1=vertex(1,1:2,iline_j,iobj)
        vj_2=vertex(2,1:2,iline_j,iobj)
        if (jj==1) then
          cycle
        endif
        lshare_1=((v1_1(1)==vj_1(1)).and.(v1_1(2)==vj_1(2)))
        lshare_2=((v1_1(1)==vj_2(1)).and.(v1_1(2)==vj_2(2)))
        lshare_3=((v1_2(1)==vj_1(1)).and.(v1_2(2)==vj_1(2)))
        lshare_4=((v1_2(1)==vj_2(1)).and.(v1_2(2)==vj_2(2)))

        if (lshare_1 .or. lshare_2 .or. &
            lshare_3 .or. lshare_4) then
          exit
        endif

      enddo
      if (jj>n_nearest_line) call fatal_error('find_IB_point','Not sharing a vertex')
! For debug
!     open (unit=555,file='data/sharing_vertex',position='APPEND')
!     if (jj>n_nearest_line) write(555,*) 'Not sharing a vertex'
!     close(555)
!
!  Initialize dist_standard for the loop of nearest facets
      dist_standard=impossible
!
!  Initialize lfind_proj_point
!     lfind_proj_point=.false.
!
!  Find the normal intercepts(IB points),with the min distance to gp.
!
!  This file is for debug
!     open(228, file='data/rproj_point', form='formatted')
!
      do jj=1, n_nearest_line
        iline=nearest_line_index(jj)
!
        if (lnearest_line(iline)) then
!
! For debug
!         if (it==1) then
!           open (unit=556,file='data/lnearest_line',position='APPEND')
!           write(556,*) vertex(1,1:2,iline,iobj), iline
!           write(556,*) vertex(2,1:2,iline,iobj), iline
!           close(556)
!         endif 
!
!  Find coordinates of the projection
          call find_projection(xa,yb,zc,iline,IB_x_temp,IB_y_temp, &
                               IB_z_temp,iobj)
!
!  Output for the sake of debug
!         write(228,*) IB_x_temp,IB_y_temp
!
!  Check if above IB point lies within line
          call is_point_on_line(IB_x_temp,IB_y_temp,IB_z_temp,&
                                iline,lon_line,iobj)
!
! For debug
!         if (it==1) then
!           open (229, file='data/find_proj_message', position='APPEND')
!           write(229,*) 'after is_point_on_line'
!         endif
!         close(229)
!
          if (lon_line) then
            if (present(lfind_proj_point)) lfind_proj_point=.true.
!
! For debug
!           open (unit=666,file='data/find_ib_message',position='APPEND')
!           if (it==1) write(666,*) 'Find IB point for solid point surrounding mirror point'
!           if (it==1) write(666,*) lon_line
!           close(666)
!
!  Distance between IB point and ghost point. |a|cos(theta)=(a dot b)/|b|
!
            gp_v1_vec(1)=vertex(1,1,iline,iobj)-xa
            gp_v1_vec(2)=vertex(1,2,iline,iobj)-yb
            line_n_norm=sqrt(norm_line(1,iline,iobj)**2 + &
			     norm_line(2,iline,iobj)**2)
!
            gp_ib_dist=abs((gp_v1_vec(1)*norm_line(1,iline,iobj)+ &
			    gp_v1_vec(2)*norm_line(2,iline,iobj))/ &
                           (line_n_norm))

            if (gp_ib_dist<dist_standard) then
              IB_x=IB_x_temp
              IB_y=IB_y_temp
              IB_z=IB_z_temp
!
              line_index=iline
!
              dist_standard=gp_ib_dist
            endif
!
          endif ! end if lon_line
!
        else
!
          if (lroot) print*,'lnearest_line should be true', &
			    'for range jj=1 to n_nearest_line', &
                            'as has been assigned above'					  
          call fatal_error('find_IB_point','lnearest_line should be true')
!
        endif ! end "if" lnearest_line=T
!
      enddo ! end loop of nearest lines that sharing the nearest vertex
!
! For debug
!     if (it==1) then
!       open (230, file='data/find_proj_message', position='APPEND')
!       write(230,*) 'after loop of nearest lines that sharing the nearest vertex'
!     endif
!     close(230)
!
!     close(228)
!
      end subroutine find_IB_point
!***********************************************************************
      subroutine find_IB_point_mir(xa,yb,zc,IB_x,IB_y,IB_z,line_index,iobj,lfind_proj_point)
!
!  For a given ghost point, find the according IB point. IB points are defined as
!  the normal intercepts between the normal vector starting from the given ghost point 
!  normal to solid-fluid interface and the solid boundary. 
!  Firstly, the nearest vertex to the ghost point is found. Secondly,  find the normal
!  intercept between the ghost point and the facets that share the nearest vertex.    
!
!  15-09-2015/chaoli: coded
!
      real, intent(in)    :: xa,yb,zc
      integer, intent(in) :: iobj
      integer :: iline
      integer, intent(out) :: line_index
      real    :: gp_v1_dist,gp_v2_dist,gp_v3_dist,gp_v_dist
!
!  Coordinates of the nearest vertex of the three/two ones of a line
      real, dimension(2)  :: nearest_of_three
      real    :: dist_standard
!
!  Coordinates of the nearest vertex to the given ghost point
      real, dimension(2)  :: nearest_vertex
!
      logical, dimension(n_lines(iobj)) :: lnearest_line
!
      integer             :: ii, jj
      integer             :: n_nearest_line
!
      real, dimension(2)  :: gp_v1_vec
      real                :: line_n_norm
      real                :: IB_x_temp,IB_y_temp,IB_z_temp
      real, intent(inout) :: IB_x,IB_y,IB_z
!
      real                :: gp_ib_dist
      logical             :: lon_line
      logical, optional   :: lfind_proj_point

      integer, parameter  :: N_MAX_FACET=20 ! should be modified as needed
!
!  A recorder to store the index iline of the nearest line
!  20 can be modified as needed
      integer, dimension(5)  :: nearest_line_index 
!
!  For debug
      real,    dimension(2)  :: v1_1, v1_2, vj_1, vj_2
      integer                :: iline_1, iline_j
      logical                :: lshare_1, lshare_2, lshare_3, lshare_4
!
! For debug
!        if (it==1) then
!          open (229, file='data/find_IB_mir_message', position='APPEND')
!          write(229,*) 'before find_IB_mir'
!        endif
!        close(229)
!
!  Initialize dist_standard
      dist_standard=impossible
!
!  Initialize the logical matrix lnearest_line to be .true.
!  for all facets
      lnearest_line(:)=.true.
!
!  Loop over all the lines
      do iline=1, n_lines(iobj)
!
!  Find the distance between ghost point and the three/two verteices
        gp_v1_dist=sqrt((xa-vertex(1,1,iline,iobj))**2+ &
                        (yb-vertex(1,2,iline,iobj))**2)
        gp_v2_dist=sqrt((xa-vertex(2,1,iline,iobj))**2+ &
                        (yb-vertex(2,2,iline,iobj))**2)
! 
        gp_v_dist=min(gp_v1_dist, gp_v2_dist)
!
! For debug
!        if (it==1) then
!          open (232, file='data/find_IB_mir_message3', position='APPEND')
!          write(232,*) gp_v_dist
!        endif
!        close(232)
!
!  Determine which one of the two verteices is the nearest one
        if (gp_v_dist==gp_v1_dist) then
	  nearest_of_three(:)=vertex(1,1:2,iline,iobj)
        elseif (gp_v_dist==gp_v2_dist) then
	  nearest_of_three(:)=vertex(2,1:2,iline,iobj)
        endif
!
!  Find the nearest vertex to a given ghost
        if (gp_v_dist<dist_standard) then
          nearest_vertex(:)=nearest_of_three(1:2)
          dist_standard=gp_v_dist
!
!  For the situation where the nearest vertex is IB point
          line_index=iline
!
          if (iline>=2) then
            lnearest_line(1:iline-1)=.false.
          endif
!
        elseif (gp_v_dist>dist_standard) then
          lnearest_line(iline)=.false.
        else
!         do nothing
        endif ! end if "gp_v_dist<dist_standard"
!
!{  Initialize the number of nearest facets. 
!       n_nearest_line=0
!
!  The total number of facets that are nearest to a given point
!       do ii=1, iline
!
!         if (lnearest_line(ii)) then
!           n_nearest_line=n_nearest_line+1
!           nearest_line_index(n_nearest_line)=ii
!         endif
!
!       enddo
!
!  N_MAX_FACET is MAX number of facets around a common vertex.
!       if (n_nearest_line>=N_MAX_FACET) exit }(more efficient??)
!
      enddo ! end loop of all lines
!
!  Assign the coordinates of the nearest vertex to IB point, temperorarily
!  If we fail to find a projection as IB point, then this is the situation.
      IB_x=nearest_vertex(1)
      IB_y=nearest_vertex(2)
      IB_z=zc
!
!  Initialize the number of nearest facets.
      n_nearest_line=0
!
!  The total number of facets surrounding the nearest vertex
!  that are closest to a given ghost point
      do ii=1, n_lines(iobj)
!
        if (lnearest_line(ii)) then
          n_nearest_line=n_nearest_line+1
          nearest_line_index(n_nearest_line)=ii
        endif
!
      enddo
!
!     if (lroot) print*,'number of nearest line element', n_nearest_line
!  Check if n_nearest_line equal to 0
      if (n_nearest_line==0) &
        call fatal_error ('find_IB_point_mir','n_nearest_line=0')
!
! For debug
!        if (it==1) then
!          open (231, file='data/find_IB_mir_message2', position='APPEND')
!          write(231,*) n_nearest_line
!        endif
!        close(231)
!
!  Debug
      if (n_nearest_line/=2) &
        call fatal_error ('find_IB_point_mir','n_nearest_line/=2')
!
!  Add a check here that facets with lnearest_line==.true. 
!  share a common vertex. Debug.
      do jj=1, n_nearest_line
        iline_j=nearest_line_index(jj)
        iline_1=nearest_line_index(1)        
        v1_1=vertex(1,1:2,iline_1,iobj)
        v1_2=vertex(2,1:2,iline_1,iobj)
        vj_1=vertex(1,1:2,iline_j,iobj)
        vj_2=vertex(2,1:2,iline_j,iobj)
        if (jj==1) then
          cycle
        endif
        lshare_1=((v1_1(1)==vj_1(1)).and.(v1_1(2)==vj_1(2)))
        lshare_2=((v1_1(1)==vj_2(1)).and.(v1_1(2)==vj_2(2)))
        lshare_3=((v1_2(1)==vj_1(1)).and.(v1_2(2)==vj_1(2)))
        lshare_4=((v1_2(1)==vj_2(1)).and.(v1_2(2)==vj_2(2)))

        if (lshare_1 .or. lshare_2 .or. &
            lshare_3 .or. lshare_4) then
          exit
        endif

      enddo
!
! For debug
         if (it==1) then
           open (230, file='data/find_IB_mir_message1', position='APPEND')
           write(230,*) 'after find_IB_mir'
         endif
         close(230)
!
      if (jj>n_nearest_line) call fatal_error('find_IB_point_mir','Not sharing a vertex')
!
! For debug
!     open (unit=555,file='data/sharing_vertex',position='APPEND')
!     if (jj>n_nearest_line) write(555,*) 'Not sharing a vertex'
!     close(555)
!
!  Initialize dist_standard for the loop of nearest facets
      dist_standard=impossible
!
!  Initialize lfind_proj_point
!     lfind_proj_point=.false.
!
!  Find the normal intercepts(IB points),with the min distance to gp.
!
!  This file is for debug
!     open(228, file='data/rproj_point', form='formatted')
!
      do jj=1, n_nearest_line
        iline=nearest_line_index(jj)
!
        if (lnearest_line(iline)) then
!
! For debug
!         if (it==1) then
!           open (unit=556,file='data/lnearest_line',position='APPEND')
!           write(556,*) vertex(1,1:2,iline,iobj), iline
!           write(556,*) vertex(2,1:2,iline,iobj), iline
!           close(556)
!         endif 
!
!  Find coordinates of the projection
          call find_projection(xa,yb,zc,iline,IB_x_temp,IB_y_temp, &
                               IB_z_temp,iobj)
!
!  Output for the sake of debug
!         write(228,*) IB_x_temp,IB_y_temp
!
!  Check if above IB point lies within line
          call is_point_on_line(IB_x_temp,IB_y_temp,IB_z_temp,&
                                iline,lon_line,iobj)
!
! For debug
!         if (it==1) then
!           open (229, file='data/find_proj_message', position='APPEND')
!           write(229,*) 'after is_point_on_line'
!         endif
!         close(229)
!
          if (lon_line) then
            if (present(lfind_proj_point)) lfind_proj_point=.true.
!
! For debug
!           open (unit=666,file='data/find_ib_message',position='APPEND')
!           if (it==1) write(666,*) 'Find IB point for solid point surrounding mirror point'
!           if (it==1) write(666,*) lon_line
!           close(666)
!
!  Distance between IB point and ghost point. |a|cos(theta)=(a dot b)/|b|
!
            gp_v1_vec(1)=vertex(1,1,iline,iobj)-xa
            gp_v1_vec(2)=vertex(1,2,iline,iobj)-yb
            line_n_norm=sqrt(norm_line(1,iline,iobj)**2 + &
			     norm_line(2,iline,iobj)**2)
!
            gp_ib_dist=abs((gp_v1_vec(1)*norm_line(1,iline,iobj)+ &
			    gp_v1_vec(2)*norm_line(2,iline,iobj))/ &
                           (line_n_norm))

            if (gp_ib_dist<dist_standard) then
              IB_x=IB_x_temp
              IB_y=IB_y_temp
              IB_z=IB_z_temp
!
              line_index=iline
!
              dist_standard=gp_ib_dist
            endif
!
          endif ! end if lon_line
!
        else
!
          if (lroot) print*,'lnearest_line should be true', &
			    'for range jj=1 to n_nearest_line', &
                            'as has been assigned above'					  
          call fatal_error('find_IB_point_mir','lnearest_line should be true')
!
        endif ! end "if" lnearest_line=T
!
      enddo ! end loop of nearest lines that sharing the nearest vertex
!
! For debug
!     if (it==1) then
!       open (230, file='data/find_proj_message', position='APPEND')
!       write(230,*) 'after loop of nearest lines that sharing the nearest vertex'
!     endif
!     close(230)
!
!     close(228)
!
      end subroutine find_IB_point_mir
!***********************************************************************
      subroutine find_projection(xa,yb,zc,iline,x_proj,y_proj,z_proj,iobj)
!
! 17-09-2015/chaoli: coded.
!
! This subroutine is to find the projective point on the first set of triangle 
! facets for a given ghost point.
!
      real, intent(in)        :: xa,yb,zc
      integer, intent(in)     :: iline
      integer, intent(in)     :: iobj
      real, intent(out)       :: x_proj,y_proj,z_proj

      real, dimension(2)      :: gp_v_vec
      real                    :: line_norm2
      real                    :: v1x, v1y
      real                    :: v2x, v2y
!
! Calculate the the aux coeeficient
      v1x=vertex(1,1,iline,iobj)
      v1y=vertex(1,2,iline,iobj)
      v2x=vertex(2,1,iline,iobj)
      v2y=vertex(2,2,iline,iobj)
      gp_v_vec(1)=xa-v1x
      gp_v_vec(2)=yb-v1y

      line_norm2=(v2x-v1x)**2+(v2y-v1y)**2
!
!  Find the coordinates of IB point temporarily
!
      x_proj=v1x+(gp_v_vec(1)*(v2x-v1x)**2+gp_v_vec(2)*(v2y-v1y)* &
                              (v2x-v1x))/line_norm2
      y_proj=v1y+(gp_v_vec(1)*(v2x-v1x)*(v2y-v1y)+gp_v_vec(2)* &
                              (v2y-v1y)**2)/line_norm2
      z_proj=zc

      end subroutine find_projection
!***********************************************************************
      subroutine is_point_on_line(px,py,pz,line,lon_line,iobj)
!
!   This subroutine is to determine whether the projection is on 
!   the line. 
!
!  15-09-2015/chaoli: coded.
!
      integer :: line,vv
      integer, intent(in)  :: iobj
      real, intent(in)     :: px,py,pz
      real :: dot_check
      logical, intent(out) :: lon_line
      real :: vx,vy
      real :: minval_D
      real, dimension(2) :: D
      real, dimension(2) :: v1_IB_vec
      real, dimension(2) :: v12_vec
!
!  If point p is very close to one of the line vertices,
!  consider that p belongs to line, and return.
!
      do vv = 1,2
        vx = vertex(vv,1,line,iobj)
        vy = vertex(vv,2,line,iobj)
        D(vv) = sqrt((px - vx)**2 + (py - vy)**2)
      enddo
!
      minval_D = minval(D)
!
      if (minval_D < TOL_STL) then
        lon_line = .true.
        return
      endif
!
!  line direction vector (normalized)
      v12_vec(1) = vertex(2,1,line,iobj)-vertex(1,1,line,iobj)
      v12_vec(2) = vertex(2,2,line,iobj)-vertex(1,2,line,iobj)
!
!  vector pointing from vertex1 to IB point
      v1_IB_vec(1) =px-vertex(1,1,line,iobj) 
      v1_IB_vec(2) =py-vertex(1,2,line,iobj)	
!
      dot_check=v12_vec(1)*v1_IB_vec(1)+v12_vec(2)*v1_IB_vec(2)		
!
      if (dot_check>=-TOL_STL) then
!
! For Debug
!
!       if (it==1) then
!         open(777, file='data/on_line_dot_check',position='APPEND')
!         write(777,*) dot_check
!       endif
!
        lon_line=.true.
        return
      endif

      end subroutine is_point_on_line
!***********************************************************************
      subroutine interpolate_idw(f,f_tmp,lower_i,lower_j, &
                                 lower_k,xmirror,ymirror, &
                                 zmirror,line_index,dist_mir_gp,iobj)
!
!  Interpolate value in a mirror point from the four or eight corner values
!  the inverse distance weighting interpolation(IDW) is implemented.
!  Chaudhuri et al. 2011
!
!  14-apr-2015/Yujuan: coded
!  21-sep-2015/Chaoli: adapted for the use of stl geometry
!
      use Sub 
!
      real, dimension (mx,my,mz,mfarray), intent(in)  :: f
      real, dimension (mvar),             intent(out) :: f_tmp
      integer, intent(in) :: lower_i,lower_j,lower_k
      real,    intent(in) :: xmirror,ymirror,zmirror
      real                :: dis_grid
      real, intent(in)    :: dist_mir_gp 
!
      integer                 :: ix0, iy0, iz0
      real, dimension(3)      :: xxp
      real                    :: vm_n,vm_t,vs_n
      real                    :: coe_total, ratio
      logical                 :: lfluid_point
      integer                 :: ii, jj
!
      real, dimension(2,2,mvar)   :: fvar
      real, dimension(2,2)        :: dij, alfa, coe
     
      real,    dimension(3,2)     :: cornervalue
      integer, dimension(3,2)     :: cornerindex
!
      integer, intent(in)         :: line_index
!
      real                        :: rho_ghost_point,rho_mirror_point
      real                        :: T_mirror_point,T_ghost_point
!
      integer                     :: i,j
!
      real, parameter             :: eps_1=10.0e-10
      real, dimension(3)          :: norm_vec, tan_vec
      integer, intent(in)         :: iobj
!
!  Prepare help parameters
!
      xxp=(/xmirror,ymirror,zmirror/)
!
      coe_total=0.0; f_tmp=0.0 
!
      ix0=lower_i
      iy0=lower_j
      iz0=lower_k
!
!  Find 4 or 8 corner points.
!
      call find_corner_points(.false.,cornervalue,cornerindex, &
	                           ix0,iy0,iz0)
!
!  Find the length of diagonal line of the grid cell that surrounds the mirror point
      dis_grid=sqrt((cornervalue(1,1)-cornervalue(1,2))**2+ &
	            (cornervalue(2,1)-cornervalue(2,2))**2)
!
!  Find the distance rij from all 4 or 8 corner points to the mirror point
!
      do i=1,2; do j=1,2
!
        ii=cornerindex(1,i)
        jj=cornerindex(2,j)

        dij(i,j)=sqrt((cornervalue(1,i)-xxp(1))**2+(cornervalue(2,j)-xxp(2))**2)
!
        fvar(i,j,:)=f(ii,jj,n1,1:mvar)
!
        lfluid_point=((ba(ii,jj,n1,1)==0).and. &
                      (ba(ii,jj,n1,2)==0).and.(ba(ii,jj,n1,3)==0))
        if (lfluid_point) then
          alfa(i,j)=1.0
        else
          alfa(i,j)=0.0
        endif
!
        coe(i,j)=alfa(i,j)*1.0/(dij(i,j)**2)
        coe_total=coe_total+coe(i,j)
!
      enddo; enddo
!
! Calculate the value of mirror point
! If the mirror point is very close to a grid point, special treatment is needed!
!
      do i=1,2; do j=1,2
        ii=cornerindex(1,i)
        jj=cornerindex(2,j)
        lfluid_point=((ba(ii,jj,n1,1)==0).and. &
		      (ba(ii,jj,n1,2)==0).and. &
                      (ba(ii,jj,n1,3)==0))
	ratio=dij(i,j)/dis_grid
        if ((ratio.le.eps_1) .and. (lfluid_point)) then
! f_tmp is to store the variable value on mirror point
          f_tmp(1:mvar)=fvar(i,j,1:mvar)
          exit
        endif
        f_tmp(1:mvar)=f_tmp(1:mvar)+(coe(i,j)/coe_total)*fvar(i,j,1:mvar)
      enddo; enddo
!
! Boundary conditions is implemented.
!
! Non_slip boundary condition for velocity.
!
      if (bc_vel_type=='no-slip') then
        f_tmp(1:3)=2.0*solid_vel(1:3)-f_tmp(1:3)
      elseif (bc_vel_type=='slip') then
        norm_vec(1:3)=norm_line(1:3,line_index,iobj)
        tan_vec(1:3)=loc_coor_utvec(1:3,line_index,iobj)
        call dot(norm_vec(1:3),f_tmp(iux:iuz),vm_n)
        call dot(tan_vec(1:3),f_tmp(iux:iuz),vm_t)
        call dot(norm_vec(1:3),solid_vel,vs_n)

        f_tmp(iux:iuz)=(2.0*vs_n-vm_n)*norm_vec(1:3)+vm_t*tan_vec(1:3)
      endif
!
! BC for temperature.
!
      if (ilnTT>0) then
        T_mirror_point=f_tmp(ilnTT)
        if (.not. ltemperature_nolog) T_mirror_point=exp(T_mirror_point)
        if (bc_thermo_type=='Dirichlet') then
          T_ghost_point=2*solid_temp-T_mirror_point
        elseif (bc_thermo_type=='Neumann') then
          T_ghost_point=dist_mir_gp*solid_flux+T_mirror_point
        elseif (bc_thermo_type=='Robin') then
          T_ghost_point=4.0*dist_mir_gp*solid_tempflux/ &
		        (2.0*coefa+2.0*coefb*dist_mir_gp)+ &
                        (2.0*coefa-2.0*coefb*dist_mir_gp)/ &
                        (2.0*coefa+2.0*coefb*dist_mir_gp)*T_mirror_point
        endif
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
      rho_mirror_point=f_tmp(ilnrho)
      if (.not. ldensity_nolog) rho_mirror_point=exp(rho_mirror_point)
      if (ilnTT>0) then
        rho_ghost_point=rho_mirror_point*T_mirror_point/T_ghost_point
      else
        rho_ghost_point=rho_mirror_point
      endif
      if (.not. ldensity_nolog) then
        f_tmp(ilnrho)=log(rho_ghost_point)
      else
        f_tmp(ilnrho)=rho_ghost_point
      endif
!
      endsubroutine interpolate_idw
!***********************************************************************
      subroutine interpolate_linear(f,f_tmp,lower_i,lower_j, &
                                    lower_k,xmirror,ymirror, &
                                    zmirror,xib,yib,zib,line_index,iobj, &
                                    dist_mir_gp)
!
!  Interpolate value in a mirror point from the four or eight corner values
!  linear interpolation method is implemented.
!
      use Sub
      use General, only: linear_interpolate
!
      real, dimension (mx,my,mz,mfarray), intent(in)  :: f
      real, dimension (mvar),             intent(out) :: f_tmp
      real, dimension (mvar)                         :: fvar
      integer, intent(in) :: lower_i,lower_j,lower_k
      real,    intent(in) :: xmirror,ymirror,zmirror
      real,    intent(in) :: xib,yib,zib
!
      integer                 :: ix0, iy0, iz0
      integer                 :: ix1, iy1, iz1
      real, dimension(3)      :: xxm, xxib,xxg
      integer, dimension(3)   :: inear, inearg
!
      real                        :: fpvar_T
      real                        :: fpvar_rho
     
      real,    dimension(3,2)     :: cornervalue
      integer, dimension(3,2)     :: cornerindex
      logical                 :: lfluid_point1, lfluid_point2
      logical                 :: lfluid_point3, lfluid_point4
      real                    :: dist_ib_mir, dist_mir_gri
      real                    :: dist_ib_gri
!
      integer, intent(in)         :: line_index
      integer, intent(in)         :: iobj
      real                        :: vm_n,vm_t,vs_n
      real                        :: rho_ghost_point,rho_mirror_point
      real                        :: T_mirror_point,T_ghost_point
      real, dimension(3)          :: norm_vec, tan_vec
      real, intent(in)            :: dist_mir_gp
!
!  Prepare help parameters
!
      ix0=lower_i
      iy0=lower_j
      iz0=lower_k
      ix1=ix0+1
      iy1=iy0+1
      iz1=iz0+1
! 
      xxm=(/xmirror,ymirror,zmirror/)
      xxib=(/xib,yib,zib/)
      inear=(/lower_i,lower_j,lower_k/)
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
! To check if there is a solid point among the 4 or 8 grid points surrounding a mirror point
      lfluid_point1=((ba(ix0,iy0,n1,1)==0).and. &
                     (ba(ix0,iy0,n1,2)==0).and.(ba(ix0,iy0,n1,3)==0))
      lfluid_point2=((ba(ix1,iy0,n1,1)==0).and. &
                     (ba(ix1,iy0,n1,2)==0).and.(ba(ix1,iy0,n1,3)==0))
      lfluid_point3=((ba(ix0,iy1,n1,1)==0).and. &
                     (ba(ix0,iy1,n1,2)==0).and.(ba(ix0,iy1,n1,3)==0))
      lfluid_point4=((ba(ix1,iy1,n1,1)==0).and. &
                     (ba(ix1,iy1,n1,2)==0).and.(ba(ix1,iy1,n1,3)==0))
!
      if (lfluid_point1 .and. lfluid_point2 .and. &
          lfluid_point3 .and. lfluid_point4) then
!
! f_tmp is to store the variable value on mirror point
        if (.not. linear_interpolate(f,iux,mvar,xxm,f_tmp,inear,.false.))&
          call fatal_error('linear_interpolate','interpolate_linear')
!
! There is at least one solid point
      else
!  Find the closest grid plane or line
          call find_closest_grid_plane(xxg,inearg,xxm,xxib,cornervalue,cornerindex)
! 
! Calculate the variable value on closest grid plane or line
! fvar is to store the variable value on closest grid line point
          if (.not. linear_interpolate(f,iux,mvar,xxg,fvar,inearg,.false.))&
          call fatal_error('linear_interpolate','interpolate_linear grid point')
! Calculate several aux distance
          dist_ib_mir=sqrt((xxm(1)-xxib(1))**2+(xxm(2)-xxib(2))**2)
          dist_ib_gri=sqrt((xxg(1)-xxib(1))**2+(xxg(2)-xxib(2))**2)
          dist_mir_gri=sqrt((xxm(1)-xxg(1))**2+(xxm(2)-xxg(2))**2)
! Find the velocity value on mirror point
          f_tmp(iux:iuz)=(fvar(iux:iuz)*dist_ib_mir+solid_vel*dist_mir_gri)/dist_ib_gri
! Find temperature value on mirror point
          if (ilnTT>0) then
            fpvar_T=fvar(ilnTT)
            if (.not. ltemperature_nolog) fpvar_T=exp(fpvar_T)
            if (bc_thermo_type=='Dirichlet') then
              f_tmp(ilnTT)=(fpvar_T*dist_ib_mir+solid_temp*dist_mir_gri)/dist_ib_gri
            elseif (bc_thermo_type=='Neumann') then
              f_tmp(ilnTT)=solid_flux*dist_mir_gri+fpvar_T
            elseif (bc_thermo_type=='Robin') then
              f_tmp(ilnTT)=(dist_mir_gri*solid_tempflux+ &
                            (coefb*dist_ib_mir+coefa)*fpvar_T)/(coefb*dist_ib_gri+coefa)
            endif
          endif
! Find density value on mirror point
          fpvar_rho=fvar(ilnrho)
          if (.not. ldensity_nolog) then
            fpvar_rho=exp(fpvar_rho)
          endif
          if (ilnTT>0) then
            f_tmp(ilnrho)=fpvar_rho*fpvar_T/f_tmp(ilnTT) ! IDEA ONE
          else
            f_tmp(ilnrho)=fpvar_rho ! IDEA ONE
          endif
          if (.not. ldensity_nolog) then
            f_tmp(ilnrho)=log(f_tmp(ilnrho))
          endif
        endif
!
! Boundary conditions is implemented.
!
! Non_slip boundary condition for velocity.
!
      if (bc_vel_type=='no-slip') then
        f_tmp(1:3)=2.0*solid_vel(1:3)-f_tmp(1:3)
      elseif (bc_vel_type=='slip') then
        norm_vec(1:3)=norm_line(1:3,line_index,iobj)
        tan_vec(1:3)=loc_coor_utvec(1:3,line_index,iobj)
        call dot(norm_vec(1:3),f_tmp(iux:iuz),vm_n)
        call dot(tan_vec(1:3),f_tmp(iux:iuz),vm_t)
        call dot(norm_vec(1:3),solid_vel,vs_n)

        f_tmp(iux:iuz)=(2.0*vs_n-vm_n)*norm_vec(1:3)+vm_t*tan_vec(1:3)
      endif
!
! BC for temperature.
!
      if (ilnTT>0) then
        T_mirror_point=f_tmp(ilnTT)
        if (.not. ltemperature_nolog) T_mirror_point=exp(T_mirror_point)
        if (bc_thermo_type=='Dirichlet') then
          T_ghost_point=2*solid_temp-T_mirror_point
        elseif (bc_thermo_type=='Neumann') then
          T_ghost_point=dist_mir_gp*solid_flux+T_mirror_point
        elseif (bc_thermo_type=='Robin') then
          T_ghost_point=4.0*dist_mir_gp*solid_tempflux/ &
		        (2.0*coefa+2.0*coefb*dist_mir_gp)+ &
                        (2.0*coefa-2.0*coefb*dist_mir_gp)/ &
                        (2.0*coefa+2.0*coefb*dist_mir_gp)*T_mirror_point
        endif
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
      rho_mirror_point=f_tmp(ilnrho)
      if (.not. ldensity_nolog) rho_mirror_point=exp(rho_mirror_point)
      if (ilnTT>0) then
        rho_ghost_point=rho_mirror_point*T_mirror_point/T_ghost_point
      else
        rho_ghost_point=rho_mirror_point
      endif
      if (.not. ldensity_nolog) then
        f_tmp(ilnrho)=log(rho_ghost_point)
      else
        f_tmp(ilnrho)=rho_ghost_point
      endif
!
      endsubroutine interpolate_linear
!***********************************************************************
      subroutine interpolate_matrix(f,f_tmp,lower_i,lower_j, &
                                    lower_k,xmirror,ymirror, &
                                    zmirror,xghost,yghost,zghost, &
                                    line_index,dist_mir_gp,iobj)
! 
! If interpolate point is inside the geometry, then it will be 
! substituted by the corresponding boundary intercept point
! Bilinear interpolation is used for 2D case
! Trilinear interpolation is used for 3D case
! Only non-slip boundary condion for velocity can be implemented.
! Cindy Merlin et al.(2013)
!
! Oct-20-2015//: Coded by Yujuan
! Oct-20-2015//: Adapted by chaoli for the use of stl geometry.
!
      use Messages, only: fatal_error
! 
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mvar), intent(out)            :: f_tmp
      real, intent(in)    :: xmirror, ymirror, zmirror
      integer, intent(in) :: lower_i, lower_j, lower_k
!
      real, dimension(3,2)    :: cornervalue
      integer, dimension(3,2) :: cornerindex
      real, dimension(3)      :: xxp
! 
      real, dimension(3)      :: norm_vec, tan_vec     
      integer                 :: i, j, k, counter, ivar, nrow
      integer                 :: ix0, iy0, iz0
      real                    :: xib, yib, zib

      real                    :: rho_ghost_point,rho_mirror_point
      real                    :: T_mirror_point,T_ghost_point
!
      real, allocatable, dimension(:,:,:) :: xyz_matrix
      real, allocatable, dimension(:,:)   :: coe_matrix
      real, allocatable, dimension(:)     :: var_matrix
      real, allocatable, dimension(:,:)   :: fvar
      real, allocatable, dimension(:,:)   :: inpoint
!
      logical    :: lfluid_point
      integer    :: ii, jj
      real       :: xm, ym, zm
!
      integer               :: ib1, ib2
      real                  :: minval
!
!  line_index is to record the index of a lineelement on which the ghost point has 
!  its IB point. the according mirror point is xxp. This is line_index is used to implemented 
!  slip velocity boundary condition.
      integer, intent(in)    :: line_index, iobj
      real,    intent(in)    :: dist_mir_gp
!
!  To record the index of an lineelement on which a non-fluid grid point has its 
!  IB point
      integer    :: line_index_nf
      integer    :: iobj_nf
!
      integer    :: iline, iline_nearest
      real       :: mp_centr_dist, dist_standard
!
      real, intent(in)      :: xghost, yghost, zghost 
      real                  :: xmirror1, ymirror1, zmirror1
      real                  :: dist_mir_gp1
      real, dimension(3)    :: xxp1
      integer  :: lower_i1, upper_i1
      integer  :: lower_j1, upper_j1
      integer  :: lower_k1, upper_k1
      integer  :: ix01, iy01, iz01
      real, dimension(3,2)    :: cornervalue1
      integer, dimension(3,2) :: cornerindex1
      integer                 :: counter1
      real, allocatable, dimension(:,:,:) :: xyz_matrix1
      real, allocatable, dimension(:,:)   :: coe_matrix1
      real, allocatable, dimension(:)     :: var_matrix1
      real, allocatable, dimension(:,:)   :: fvar1
      real, dimension(mvar)               :: f_tmp1
      real                                :: T_mirror_point1
      real                                :: rho_mirror_point1
      logical                             :: lfluid_point1
!
!  Calculate coordinates of another one mirror point
      dist_mir_gp1=dist_mir_gp+sqrt(2.0)*dxmin
      xmirror1=xghost+(xmirror-xghost)*dist_mir_gp1/dist_mir_gp
      ymirror1=yghost+(ymirror-yghost)*dist_mir_gp1/dist_mir_gp
      zmirror1=zghost
!
! For debug
!     open(unit=336,file='data/2nd_mirror_point_coor',position='APPEND')
!     if (it==1) write(336,*) xmirror1,ymirror1
!     close(336)
!
! For debug
!     open(unit=337,file='data/1st_mirror_point_coor',position='APPEND')
!     if (it==1) write(337,*) xmirror,ymirror
!     close(337)
!
      xxp=(/xmirror,ymirror,zmirror/)
      xxp1=(/xmirror1,ymirror1,zmirror1/)
!
      call find_near_indeces(lower_i1,upper_i1,lower_j1,upper_j1,lower_k1, &
                                 upper_k1,x,y,z,xxp1)
!
!  Find 4 points, now only for 2D case.
!
      ix0=lower_i; iy0=lower_j; iz0=lower_k
      call find_corner_points(.false.,cornervalue,cornerindex,ix0,iy0, &
                              iz0)
      ix01=lower_i1; iy01=lower_j1; iz01=lower_k1
      call find_corner_points(.false.,cornervalue1,cornerindex1,ix01,iy01, &
                              iz01)
!
!  Now only for 2D case
      nrow=4
!
      if (.not. allocated(xyz_matrix)) allocate(xyz_matrix(2,nrow,nrow))
      if (.not. allocated(coe_matrix)) allocate(coe_matrix(nrow,mvar))
      if (.not. allocated(var_matrix)) allocate(var_matrix(nrow))
      if (.not. allocated(fvar))       allocate(fvar(nrow,mvar))
      if (.not. allocated(inpoint))    allocate(inpoint(nrow,3))
!
      if (.not. allocated(xyz_matrix1)) allocate(xyz_matrix1(2,nrow,nrow))
      if (.not. allocated(coe_matrix1)) allocate(coe_matrix1(nrow,mvar))
      if (.not. allocated(var_matrix1)) allocate(var_matrix1(nrow))
      if (.not. allocated(fvar1))       allocate(fvar1(nrow,mvar))
!
!  Counter is used to record the index of surrounding(mir point) fluid points
!
      counter=0
!
!  Open a file to save the surrounding grid points which is not a solid point, for debug
!      open(unit=330,file='data/surroundfluid_point',position='APPEND')
!      open(unit=332,file='data/xyz_matrix1',position='APPEND')
!      open(unit=333,file='data/xyz_matrix2',position='APPEND')
!      open(unit=334,file='data/var_matrix',position='APPEND')
!      open(unit=331,file='data/message_file',position='APPEND')
!
!
! For debug 
!     if (it==1) then
!       open (231, file='data/interpolation_matrix_message', position='APPEND')
!       write(231,*) 'before loop'
!     endif
!     close(231)
!  For the first mirror point
      do i=1,2
      do j=1,2
       counter=counter+1

       ii=cornerindex(1,i)
       jj=cornerindex(2,j)
       lfluid_point=((ba(ii,jj,n1,1)==0).and. &
                     (ba(ii,jj,n1,2)==0).and. &
                     (ba(ii,jj,n1,3)==0))
!
       if (lfluid_point) then
! 
! for debug
!       if (it==1) write(330,*) cornervalue(1,i), cornervalue(2,j)
!
        inpoint(counter,:)=(/cornervalue(1,i),cornervalue(2,j),z(n1)/)
        xyz_matrix(:,counter,1)=inpoint(counter,1)*inpoint(counter,2)
        xyz_matrix(:,counter,2)=inpoint(counter,1)
        xyz_matrix(:,counter,3)=inpoint(counter,2)
        xyz_matrix(:,counter,4)=1.0
        fvar(counter,1:mvar)=f(ii,jj,n1,1:mvar)
!
!  Boundary condtion for pressure(rho*T) is Neumann type
!  For isothermo flow, there is no need to multiply Temperature
        if (ilnTT>0) fvar(counter,ilnrho)=f(ii,jj,n1,ilnrho)*f(ii,jj,n1,ilnTT)
       else
!
! For debug
!        open(unit=338,file='data/1st_mirror_point',position='APPEND')
!        if (it==1) write(338,*) 'there is a solid point surrounding 1st mirror point'
!        close(338)
!  If not a fluid point, find the IB point on fluid-solid interface instead.
         xm=cornervalue(1,i)
         ym=cornervalue(2,j)
         zm=z(n1)
!
!       if (it==1) write(330,*) xa, yb  ! for debug
!
         iobj_nf=iobj
!        if (it==1) write(330,*) iobj_nf  ! for debug
!
! For debug 
!        if (it==1) then
!          open (233, file='data/find_IB_message', position='APPEND')
!          write(233,*) 'before call find_IB_point_mir'
!        endif
!        close(233)
!
         call find_IB_point_mir(xm,ym,zm,xib,yib,zib,line_index_nf,iobj_nf)
!
! For debug 
!        if (it==1) then
!          open (234, file='data/find_IB_message2', position='APPEND')
!          write(234,*) 'after call find_IB_point_mir'
!        endif
!        close(234)
!
! Alternatively(to replace find_IB_point_mir).
!        dist_standard=impossible
!        do iline=1, n_lines(iobj)
!          mp_centr_dist=sqrt((xm-centr_line(1,iline,iobj))**2 + &
!                             (ym-centr_line(2,iline,iobj))**2)
!          if (mp_centr_dist<dist_standard) then
!               iline_nearest=iline
!               dist_standard=mp_centr_dist
!          endif
!        enddo
!        xib=centr_line(1,iline_nearest,iobj)
!        yib=centr_line(2,iline_nearest,iobj)
!        zib=z(n1)
!
! For debug 
!        if (it==1) then
!          open (235, file='data/find_IB_message3', position='APPEND')
!          write(235,*) 'after find_IB_point'
!        endif
!        close(235)
!
!  For debug
!        if (it==1) write(330,*) xib, yib
!
        inpoint(counter,:)=(/xib,yib,zib/)
!
!  Prepare for Dirichlet boundary condition
!
        xyz_matrix(1,counter,1)=inpoint(counter,1)*inpoint(counter,2)
        xyz_matrix(1,counter,2)=inpoint(counter,1)
        xyz_matrix(1,counter,3)=inpoint(counter,2)
        xyz_matrix(1,counter,4)=1.0
!
        fvar(counter,iux:iuz)=solid_vel(1:3)
        if (ilnTT>0) then
          if (bc_thermo_type=='Dirichlet') fvar(counter,ilnTT)=solid_temp
        endif
!
! Prepare for Neumann boundary condition 
!
        norm_vec(1:3)=norm_line(1:3,line_index_nf,iobj_nf)
        tan_vec(1:3)=loc_coor_utvec(1:3,line_index_nf,iobj_nf)
!
        xyz_matrix(2,counter,1)=norm_vec(1)*inpoint(counter,2)+norm_vec(2)*inpoint(counter,1)
        xyz_matrix(2,counter,2)=norm_vec(1)
        xyz_matrix(2,counter,3)=norm_vec(2)
        xyz_matrix(2,counter,4)=0.0
!
!  Boundary condtion for pressure(rho*T) is Neumann type, That is to set rho*T=0
!  at the immersed boundary.
        fvar(counter,ilnrho)=0.0
        if (ilnTT>0) then
          if (bc_thermo_type=='Neumann') fvar(counter,ilnTT)=-solid_flux
        endif
!
       endif
!
!       if (it==1) then
!          write(331,*) 'now we are here.'
!          write(332,*) xyz_matrix(1,counter,1)
!          write(332,*) xyz_matrix(1,counter,2)
!          write(332,*) xyz_matrix(1,counter,3)
!          write(332,*) xyz_matrix(1,counter,4)
!          write(333,*) xyz_matrix(2,counter,1)
!          write(333,*) xyz_matrix(2,counter,2)
!          write(333,*) xyz_matrix(2,counter,3)
!          write(333,*) xyz_matrix(2,counter,4)
!          write(334,*) fvar(counter,iux)
!          write(334,*) fvar(counter,iuy)
!          write(334,*) fvar(counter,iuz)
!          write(334,*) fvar(counter,ilnrho)
!        endif
!
      enddo
      enddo
!
! For debug 
!     if (it==1) then
!       open (232, file='data/interpolation_matrix2_message', position='APPEND')
!       write(232,*) 'after loop'
!     endif
!     close(232)
!
!      close(330)
!      close(331)
!      close(332)
!      close(333)
!      close(334)
!
      counter1=0
!  For the second mirror point
      do i=1,2
      do j=1,2
       counter1=counter1+1

       ii=cornerindex1(1,i)
       jj=cornerindex1(2,j)
       lfluid_point1=((ba(ii,jj,n1,1)==0).and. &
                     (ba(ii,jj,n1,2)==0).and. &
                     (ba(ii,jj,n1,3)==0))
!
       if (lfluid_point1) then
!
        xyz_matrix1(:,counter1,1)=cornervalue1(1,i)*cornervalue1(2,j)
        xyz_matrix1(:,counter1,2)=cornervalue1(1,i)
        xyz_matrix1(:,counter1,3)=cornervalue1(2,j)
        xyz_matrix1(:,counter1,4)=1.0
        fvar1(counter1,1:mvar)=f(ii,jj,n1,1:mvar)
!
!  Boundary condtion for pressure(rho*T) is Neumann type
!  For isothermo flow, there is no need to multiply Temperature
        if (ilnTT>0) fvar1(counter1,ilnrho)=f(ii,jj,n1,ilnrho)*f(ii,jj,n1,ilnTT)
!
       else
!
         call fatal_error('interpolate_matrix','there are solid points surrounding the 2nd mirror')
!
!        open(unit=335,file='data/2nd_mirror_point',position='APPEND')
!        if (it==1) write(335,*) 'warning: there is a solid point surrounding 2nd mirror point'
!        close(335)
!  If not a fluid point, find the IB point on fluid-solid interface instead.
!        xm=cornervalue1(1,i)
!        ym=cornervalue1(2,j)
!        zm=z(n1)

!        dist_standard=impossible
!        do iline=1, n_lines(iobj)
!          mp_centr_dist=sqrt((xm-centr_line(1,iline,iobj))**2 + &
!                             (ym-centr_line(2,iline,iobj))**2)
!          if (mp_centr_dist<dist_standard) then
!               iline_nearest=iline
!               dist_standard=mp_centr_dist
!          endif
!        enddo
!        xib=centr_line(1,iline_nearest,iobj)
!        yib=centr_line(2,iline_nearest,iobj)
!        zib=z(n1)
!
!  Prepare for Dirichlet boundary condition
!
!        xyz_matrix1(1,counter1,1)=xib*yib
!        xyz_matrix1(1,counter1,2)=xib
!        xyz_matrix1(1,counter1,3)=yib
!        xyz_matrix1(1,counter1,4)=1.0
!
!        fvar1(counter1,iux:iuz)=solid_vel(1:3)
!        if (ilnTT>0) then
!          if (bc_thermo_type=='Dirichlet') fvar1(counter1,ilnTT)=solid_temp
!        endif
!
! Prepare for Neumann boundary condition 
!
!        norm_vec(1:3)=norm_line(1:3,iline_nearest,iobj)
!        tan_vec(1:3)=loc_coor_utvec(1:3,iline_nearest,iobj)
!
!        xyz_matrix1(2,counter1,1)=norm_vec(1)*yib+norm_vec(2)*xib
!        xyz_matrix1(2,counter1,2)=norm_vec(1)
!        xyz_matrix1(2,counter1,3)=norm_vec(2)
!        xyz_matrix1(2,counter1,4)=0.0
!
!  Boundary condtion for pressure(rho*T) is Neumann type, That is to set rho*T=0
!  at the immersed boundary.
!        fvar1(counter1,ilnrho)=0.0
!        if (ilnTT>0) then
!          if (bc_thermo_type=='Neumann') fvar1(counter1,ilnTT)=-solid_flux
!        endif
!
       endif ! finalize if(lfluid_point1).
!
      enddo
      enddo
!
!
!  Determine the velocity coefficient matrix
!
!     open(unit=331,file='data/message_file',position='APPEND')
!
! For debug
!     open(unit=339,file='data/interpolate_m_m_file',position='APPEND')
!     if (it==1) write(339,*) 'interpolate_matrix is called'
!     close (339)
!
! For debug 
!     if (it==1) then
!       open (234, file='data/find_IB_message2', position='APPEND')
!       write(234,*) 'before velocity coefficient'
!     endif
!     close(234)
!
      do ivar=iux, iuz
        var_matrix(:)=fvar(:,ivar)
        call calc_matrix(xyz_matrix(1,:,:),var_matrix,coe_matrix(:,ivar))
!
!  There is no need to use the second mirror point for velocity
        var_matrix1(:)=fvar1(:,ivar)
        call calc_matrix(xyz_matrix1(1,:,:),var_matrix1,coe_matrix1(:,ivar))
!
!       if (it==1) write(331,*) 'velocity coefficient ready.'
      enddo
!
! For debug 
!     if (it==1) then
!       open (235, file='data/find_IB_message3', position='APPEND')
!       write(235,*) 'after velocity coefficient'
!     endif
!     close(235)
!
!     close(331)
! 
!
! For debug 
!    if (it==1) then
!      open (234, file='data/find_IB_message2', position='APPEND')
!      write(234,*) 'before temperature coefficient'
!    endif
!    close(234)
!  Determine the temperature coefficient matrix
! 
!     open(unit=331,file='data/message_file',position='APPEND')
      if (ilnTT>0) then
        if (.not. ltemperature_nolog) then
          call fatal_error('interpolate_matrix',&
              'Due to interpolation it is not correct to use lnTT!')
        endif
        var_matrix(:)=fvar(:,ilnTT)
        var_matrix1(:)=fvar1(:,ilnTT)
        if (bc_thermo_type=='Dirichlet') then
          call calc_matrix(xyz_matrix(1,:,:),var_matrix,coe_matrix(:,ilnTT))
          call calc_matrix(xyz_matrix1(1,:,:),var_matrix1,coe_matrix1(:,ilnTT))
        endif
        if (bc_thermo_type=='Neumann') then
          call calc_matrix(xyz_matrix(2,:,:),var_matrix,coe_matrix(:,ilnTT))
!
!  The second mirror point
          call calc_matrix(xyz_matrix1(2,:,:),var_matrix1,coe_matrix1(:,ilnTT))
        endif
!
! For debug 
!    if (it==1) then
!      open (235, file='data/find_IB_message3', position='APPEND')
!      write(235,*) 'after temperature coefficient'
!    endif
!    close(235)
!
!       if (it==1) write(331,*) 'temperature coefficient ready.'
      endif
!     close(331)
!
!  Determine the density coefficient matrix
!
!     open(unit=331,file='data/message_file',position='APPEND')
!
! For debug 
!    if (it==1) then
!      open (236, file='data/find_IB_message4', position='APPEND')
!      write(236,*) 'before density coefficient'
!    endif
!    close(236)
!
      if (.not. ldensity_nolog) then
        call fatal_error('interpolate_matrix',&
              'Due to interpolation it is not correct to use lnrho!')
      else
        var_matrix(:)=fvar(:,ilnrho)
        call calc_matrix(xyz_matrix(2,:,:),var_matrix,coe_matrix(:,ilnrho))
!
        var_matrix1(:)=fvar1(:,ilnrho)
        call calc_matrix(xyz_matrix1(2,:,:),var_matrix1,coe_matrix1(:,ilnrho))
!
!       if (it==1) write(331,*) 'density coefficient ready.'
      endif
!
! For debug 
!    if (it==1) then
!      open (237, file='data/find_IB_message5', position='APPEND')
!      write(237,*) 'after density coefficient'
!    endif
!    close(237)
!    close(331)
!
!  For debug
!      do ib1=1, nrow
!        do ib2=1, mvar
!          minval=abs(coe_matrix(ib1,ib2))
!          if (minval<=1e-9) then
!            call fatal_error('interpolate_matrix','zero in coe_matrix')
!          endif
!        enddo
!      enddo
!
!     if (lroot) print*,'coe_matrix is ready'
! 
!  Calculate final results on mirror point, now only for 2D case.
!
      f_tmp(1:mvar)=coe_matrix(1,1:mvar)*xxp(1)*xxp(2)+coe_matrix(2,1:mvar)*xxp(1)+ &
                    coe_matrix(3,1:mvar)*xxp(2)+coe_matrix(4,1:mvar)
!
      f_tmp1(1:mvar)=coe_matrix1(1,1:mvar)*xxp1(1)*xxp1(2)+coe_matrix1(2,1:mvar)*xxp1(1)+ &
                     coe_matrix1(3,1:mvar)*xxp1(2)+coe_matrix1(4,1:mvar)
!
! Boundary conditions is implemented.
!
! Non-slip boundary for velocity, for now.
!
      f_tmp(1:3)=2.0*solid_vel(1:3)-f_tmp(1:3)
!
! Boundary conditions for temperature.
!
      if (ilnTT>0) then
        T_mirror_point=f_tmp(ilnTT);T_mirror_point1=f_tmp1(ilnTT)
        if (.not. ltemperature_nolog) then
          T_mirror_point=exp(T_mirror_point)
          T_mirror_point1=exp(T_mirror_point1)
        endif
        if (bc_thermo_type=='Dirichlet') then
          T_ghost_point=2*solid_temp-T_mirror_point
        elseif (bc_thermo_type=='Neumann') then
          T_ghost_point=0.5*(T_mirror_point+T_mirror_point1)+ &
                        0.5*dist_mir_gp*solid_flux+ &
                        0.5*dist_mir_gp1/(sqrt(2.0)*dxmin)* &
                        (T_mirror_point-T_mirror_point1)
!         T_ghost_point=dist_mir_gp*solid_flux+T_mirror_point
        else
          call fatal_error('interpolate_matrix',&
              'now only for Dirichlet and Neumann boundary condition')
        endif
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
      rho_mirror_point=f_tmp(ilnrho);rho_mirror_point1=f_tmp1(ilnrho)
      if (.not. ldensity_nolog) then
        rho_mirror_point=exp(rho_mirror_point)
        rho_mirror_point1=exp(rho_mirror_point1)
      endif
      if (ilnTT>0) then
!
!  For non-isothermal case rho_mirror_point is rho_mirror_point*T_mirror_point in fact.
!
!       rho_ghost_point=rho_mirror_point/T_ghost_point
        rho_ghost_point=(0.5*(rho_mirror_point+rho_mirror_point1)+ &
                         0.5*dist_mir_gp1/(sqrt(2.0)*dxmin)* &
                        (rho_mirror_point-rho_mirror_point1))/T_ghost_point
      else
!
! For isothermal case the same Temperature is multiplied and divided.
! rho_mirror_point is rho_mirror_point in fact.
!
        rho_ghost_point=0.5*(rho_mirror_point+rho_mirror_point1)+ &
                        0.5*dist_mir_gp1/(sqrt(2.0)*dxmin)* &
                        (rho_mirror_point-rho_mirror_point1)
      endif
      if (.not. ldensity_nolog) then
        f_tmp(ilnrho)=log(rho_ghost_point)
      else
        f_tmp(ilnrho)=rho_ghost_point
      endif
!
      if (allocated(xyz_matrix)) deallocate(xyz_matrix)
      if (allocated(coe_matrix)) deallocate(coe_matrix)
      if (allocated(var_matrix)) deallocate(var_matrix)
      if (allocated(fvar))       deallocate(fvar)
      if (allocated(inpoint))    deallocate(inpoint)
      if (allocated(xyz_matrix1)) deallocate(xyz_matrix1)
      if (allocated(coe_matrix1)) deallocate(coe_matrix1)
      if (allocated(var_matrix1)) deallocate(var_matrix1)
      if (allocated(fvar1))       deallocate(fvar1)
!
!  For debug
!     open(unit=331,file='data/message_file',position='APPEND')
!     if (it==1) write(331,*) 'interpolate_matrix successful.'
!     close(331)
!
      end subroutine interpolate_matrix
!***********************************************************************
     subroutine calc_matrix(xyz_matrix,var_matrix,coe_matrix)
!
!  20-Oct-2015//: Yujuan coded
!
      real, dimension(:,:), intent(in)    :: xyz_matrix
      real, dimension(:),   intent(in)    :: var_matrix
      real, dimension(:),   intent(out)   :: coe_matrix
!
      real, allocatable, dimension(:,:) :: ab_matrix, xyz_up
      real, allocatable, dimension(:)   :: temp, var_up
      real :: column_max, coef
      integer :: i, j, k, id_max, nrow
!
      integer               :: ib1, ib2
      real                  :: minval
!
      nrow=size(xyz_matrix(:,1))
!
      if (.not. allocated(ab_matrix)) allocate(ab_matrix(nrow,nrow+1))
      if (.not. allocated(temp))      allocate(temp(nrow+1))
      if (.not. allocated(xyz_up))    allocate(xyz_up(nrow,nrow))
      if (.not. allocated(var_up))    allocate(var_up(nrow))
! 
      ab_matrix(1:nrow,1:nrow)=xyz_matrix
      ab_matrix(1:nrow,nrow+1)=var_matrix
!
!     open(unit=342,file='data/ab_matrix',position='APPEND')
!     do i=1, nrow
!       if (it==1) write(342,*) ab_matrix(i,:)
!     enddo
!
!     open(unit=335,file='data/message1_file',position='APPEND')
!     open(unit=336,file='data/message2_file',position='APPEND')
!
!  For debug
!      do ib1=1, nrow
!        do ib2=1, nrow+1
!          minval=abs(ab_matrix(ib1,ib2))
!          if (minval<=1e-9) then
!            call fatal_error('interpolate_matrix','zero in ab_matrix')
!          endif
!        enddo
!      enddo
!
      do k=1,nrow-1
        column_max=abs(ab_matrix(k,k))
        id_max=k
!
!        if (it==1) write(335,*) 'we are here now'  ! for debug.
!
        do i=k+1,nrow
          if (abs(ab_matrix(i,k))>column_max) then
            column_max=abs(ab_matrix(i,k))
            id_max=i
          endif
        enddo
        temp=ab_matrix(id_max,:); ab_matrix(id_max,:)=ab_matrix(k,:)
        ab_matrix(k,:)=temp
        do i=k+1,nrow 
          if (abs(ab_matrix(k,k))<1e-9) call fatal_error('calc_matrix','zero be divided')
          coef=ab_matrix(i,k)/ab_matrix(k,k)
          ab_matrix(i,:)=ab_matrix(i,:)-coef*ab_matrix(k,:)
        enddo
      enddo
      xyz_up=ab_matrix(1:nrow,1:nrow)
      var_up=ab_matrix(:,nrow+1)
!
!     if (it==1) write(336,*) 'are we here now'
!
!     close(335)
!     close(336)
!
! For debug
!     open(unit=340,file='data/interpolate_c_m_file',position='APPEND')
!     open(unit=341,file='data/xyz_up_n_n',position='APPEND')
!     if (it==1) write(340,*) 'calc_matrix is called'
!     if (it==1) then
!       if (abs(xyz_up(nrow,nrow))<1e-9) write(341,*) xyz_up(nrow,nrow)
!     endif
!     close (340)
!     close (341)
!
     if (abs(xyz_up(nrow,nrow))<1e-9) call fatal_error('calc_matrix','zero be divided')
!
      coe_matrix(nrow)=var_up(nrow)/xyz_up(nrow,nrow)
      do i=nrow-1,1,-1
        coe_matrix(i)=var_up(i)
        do j=i+1,nrow
          coe_matrix(i)=coe_matrix(i)-xyz_up(i,j)*coe_matrix(j)
        enddo
        if (abs(xyz_up(i,i))<1e-9) call fatal_error('calc_matrix','zero be divided')
        coe_matrix(i)=coe_matrix(i)/xyz_up(i,i)
      enddo
!
      if (allocated(ab_matrix)) deallocate(ab_matrix)
      if (allocated(temp))      deallocate(temp)
      if (allocated(xyz_up))    deallocate(xyz_up)
      if (allocated(var_up))    deallocate(var_up)  
!
      end subroutine calc_matrix
!***********************************************************************
      subroutine find_closest_grid_plane(xxg,inearg,xxm,xxib,cornervalue,cornerindex)
!
! Find g_global based on the gridplane which is first crossed by the normal
! from the surface.
!
      use Sub
!
      real, dimension(3), intent(in) :: xxm, xxib
      real, dimension(3) :: p_local, o_global
      real :: verylarge=1e9, r0min, r0
      integer :: ndir, dir, vardir1, constdir, topbot_tmp
      integer :: topbot
      real, dimension(3,2) :: cornervalue
      integer, dimension(3,2) :: cornerindex
      integer :: lower_i, lower_j, lower_k, upper_i, upper_j, upper_k
!
      real, dimension(3), intent(out) :: xxg
      integer, dimension(3), intent(out) :: inearg
!
!  Find which grid line is the closest one in the direction
!  away from the object surface
!
!  Check the distance to all the six (four in 2D) possible surface planes (lines in 2D).
!  Define a straight line "l" which pass through the point "p" and is normal
!  to the surface of the object.
!  Pick the grid plane, OUTSIDE the point "p", that the line "l" cross first,
!  call this plane "gridplane".
!
      ndir=2
      r0min=verylarge
      p_local=xxm-xxib
      o_global=xxib
      do dir=1,ndir
!
!  Loop through all directions "dir" and find which grid plane, outside of "p",
!  that is crossed by the line "l" first. Lets call this plane "gridplane".
!  If "dir" for the first crossing is 1 then "gridplane" is the 'yz' plane
!  and so on.....
!
          do topbot_tmp=1,2
!
!  Find the variable "rl" such that
!  rl*p_local(dir)+o_global(dir)=cornervalue(dir)
!
            r0=(cornervalue(dir,topbot_tmp)-o_global(dir))/p_local(dir)
!
!  If "rl" is larger than unity the line "l" cross the grid plane
!  outside of "p". If in addition "rl" is smaller than the smallest "rl" so
!  far (rlmin) this state must be stored since it might be the plane
!  ("gridplane") which we are looking for. In addition we must
!  find the distance, rg, from the center of the object
!  to the point where the normal cross the "gridplane".
!  The point where "l" cross "gridplane" is called "g".
!
            if (r0 > 0.0 .and. r0<r0min) then
              constdir=dir
              vardir1=mod(constdir+0,2)+1
              r0min=r0
              topbot=topbot_tmp
              xxg(constdir)=cornervalue(constdir,topbot)
              xxg(vardir1)=r0*p_local(vardir1)+o_global(vardir1)
            endif
          enddo
      enddo
      xxg(3)=xxib(3)
!
!  The direction of the normal to "gridplane", called "N_grid", is either
!  in the x or y direction since the grid is 2d-Cartesian. See table below
!  for values of constdir, vardir1, and vardir2 for different "N_grid" directions.
!
!  N_grid  constdir  vardir1 
!------------------------------
!     x        1        2      
!     y        2        1      
!
!  Check that we have found a valid distance
!
      if (r0min==verylarge) then
        print*,'o_global=',o_global
        print*,'cornerindex=',cornerindex
        print*,'cornervalue=',cornervalue
        call fatal_error('find_closest_grid_plane',&
            'A valid closest grid plane is not found!')
      endif
!
!  Due to roundoff errors the value of g_global might end up outside the
!  grid cell - in that case put it back in where it belongs
!
      if ((cornervalue(1,1)<=xxg(1).and.&
          cornervalue(1,2)>=xxg(1).or.nxgrid==1).and.&
          (cornervalue(2,1)<=xxg(2).and.&
          cornervalue(2,2)>=xxg(2).or.nygrid==1).and.&
          (cornervalue(3,1)<=xxg(3).and.&
          cornervalue(3,2)>=xxg(3).or.nzgrid==1)) then
! Everything okay
      else
        if (xxg(1)>cornervalue(1,2)) xxg(1)=cornervalue(1,2)
        if (xxg(1)<cornervalue(1,1)) xxg(1)=cornervalue(1,1)
        if (xxg(2)>cornervalue(2,2)) xxg(2)=cornervalue(2,2)
        if (xxg(2)<cornervalue(2,1)) xxg(2)=cornervalue(2,1)
        if (xxg(3)>cornervalue(3,2)) xxg(3)=cornervalue(3,2)
        if (xxg(3)<cornervalue(3,1)) xxg(3)=cornervalue(3,1)
      endif
!
!  Depending on the value of constdir the indeces of the corner points
!  specifying "gridplane" is given by lower_i, upper_i, lower_j ........
!
      if (constdir==1) then
        lower_i=cornerindex(constdir,topbot)
        upper_i=cornerindex(constdir,topbot)
        lower_j=cornerindex(vardir1,1)
        upper_j=cornerindex(vardir1,2)
        lower_k=cornerindex(3,1)
        upper_k=cornerindex(3,2)
      elseif (constdir==2) then
        lower_j=cornerindex(constdir,topbot)
        upper_j=cornerindex(constdir,topbot)
        lower_k=cornerindex(3,1)
        upper_k=cornerindex(3,2)
        lower_i=cornerindex(vardir1,1)
        upper_i=cornerindex(vardir1,2)
      endif
!
     inearg=(/lower_i,lower_j,lower_k/)
!
     end subroutine find_closest_grid_plane
!***********************************************************************
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
      real, dimension(3), optional      :: p_global
      real, dimension(3), optional      :: o_global
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
      subroutine freeze_solid_cells(df)
!
!  If we are in a solid cell (or in a cell where the value of the variables are
!  found from interpolation) set df=0 for all variables
!
!  19-nov-2008/nils: coded
!  23-sep-2015/chaoli: adapted to be available for stl geometry
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
!  temperature and the velocity components should be frozen? 
!  Modified as follows by maochaoli
!
          df(i,m,n,:)=0
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
      if (idiag_Nusselt /= 0) lpenc_requested(i_gTT)=.true.
      if (idiag_Nusselt /= 0) lpenc_requested(i_TT)=.true.
!     if (idiag_Nusselt /= 0) lpenc_requested(i_tcond)=.true.
!
      endsubroutine pencil_criteria_solid_cells
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
      endsubroutine solid_cells_clean_up
!************************************************************************
      subroutine output_solid_cells(f,df,p)
!   
      use Sub, only: dot
      use Viscosity, only: getnu
      use General, only: safe_character_append
      use General, only: linear_interpolate
      use EquationOfState, only: get_p0
!
      real, dimension (mx,my,mz,mfarray), intent(in):: f
      real, dimension (mx,my,mz,mvar), intent(in)   :: df
      type (pencil_case), intent(in)                :: p
!      
      integer :: iline, iobj
      integer :: ix0, iy0, iz0
      real    :: fpx, fpy, fpz
      real, dimension(3) :: nvec 
      real, dimension(3) :: fp_location
      real    :: fp_pressure
      real    :: force_x, force_y, force_z
      real    :: fp_stress(3,3)
!	  
      real, dimension(nx) :: nu, twonu
      real                :: fp_gradT_n
      real                :: Tobj
      real                :: deltaT
      real                :: loc_Nus, c_press
      real                :: loc_density
      character(len=60) :: Nus_string, temp_string, den_string
      character(len=120):: solid_cell_Nus, solid_cell_temp, solid_cell_den
      real              :: pp00
      real              :: drag_norm
      real              :: nusselt_norm
      integer, dimension(3) :: inear_error
      real, dimension(3) :: vel_error
      real, dimension(1) :: T_error
      integer :: lower_i, lower_j, lower_k, upper_i, upper_j, upper_k
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
      twonu=2.0*nu
!
!  Loop over all objects
!
!  Debug
!     open(444, file='data/fp_nearest', position='APPEND')
!     open(445, file='data/norm_param', position='APPEND')
      do iobj=1, nobjects
      do iline=1,n_lines(iobj)
        iy0=fpnearestgrid(iobj,iline,2)
        iz0=n
!
!  Test: Use this pencil for force calculation?
!
        if (iy0 == m .and. iz0 == n) then
          ix0=fpnearestgrid(iobj,iline,1)
          ! Test: ix0 in local domain?
          if (ix0 >= l1 .and. ix0 <= l2) then
!
!         if (it==1) write(444,*) x(ix0),y(iy0)
!
!  Acquire pressure and stress from grid point (ix0,iy0,iz0).
!  Shifting the location of the forcpoints in the theta direction
!  in order to avoid problems with autotesting
!
            fpx=centr_line(1,iline,iobj)
            fpy=centr_line(2,iline,iobj)
            fpz=z(n1)
!
            nvec(1)=norm_line(1,iline,iobj)
            nvec(2)=norm_line(2,iline,iobj)
            nvec(3)=0.0

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
              call get_p0(pp00)
              fp_stress(:,:)=twonu(ix0-nghost)*p%rho(ix0-nghost)*p%sij(ix0-nghost,:,:)
!
!  Force in x-,y-, and z-directions
!
              force_x = (-fp_pressure*nvec(1) &
                  +fp_stress(1,1)*nvec(1)+fp_stress(1,2)*nvec(2)+fp_stress(1,3)*nvec(3) &
                  ) * lineelement(iline,iobj)
!
              force_y = (-fp_pressure*nvec(2) &
                  +fp_stress(2,1)*nvec(1)+fp_stress(2,2)*nvec(2)+fp_stress(2,3)*nvec(3) &
                  ) * lineelement(iline,iobj)
!
              force_z = (-fp_pressure*nvec(3) &
                  +fp_stress(3,1)*nvec(1)+fp_stress(3,2)*nvec(2)+fp_stress(3,3)*nvec(3) &
                  ) * lineelement(iline,iobj)
!
              drag_norm=1.0/charac_len(iobj)
!             if (it==1) write(445,*) drag_norm
!
              c_dragx(iobj) = c_dragx(iobj) + force_x * drag_norm
              c_dragy(iobj) = c_dragy(iobj) + force_y * drag_norm
              c_dragz(iobj) = c_dragz(iobj) + force_z * drag_norm
!
              c_dragx_p(iobj)=c_dragx_p(iobj)-fp_pressure*nvec(1)*drag_norm*lineelement(iline,iobj)
              c_dragy_p(iobj)=c_dragy_p(iobj)-fp_pressure*nvec(2)*drag_norm*lineelement(iline,iobj)
              c_dragz_p(iobj)=c_dragz_p(iobj)-fp_pressure*nvec(3)*drag_norm*lineelement(iline,iobj)
!
              c_press=(fp_pressure-pp00)/(0.5*init_uu**2)
            endif
!
!  Local Nusselt number and average Nusselt number
!
            if (idiag_Nusselt /= 0) then
              nusselt_norm =1.0/sum_line(iobj)
!             if (it==1) write(445,*) drag_norm
              if (bc_thermo_type=='Dirichlet') then
                call dot(p%gTT(ix0-nghost,:),-nvec,fp_gradT_n)
                Tobj=solid_temp
              elseif (bc_thermo_type=='Neumann') then
                Tobj=p%TT(ix0-nghost)
                fp_gradT_n=solid_flux
              elseif (bc_thermo_type=='Robin') then
                call dot(p%gTT(ix0-nghost,:),-nvec,fp_gradT_n)
                Tobj=p%TT(ix0-nghost)
              endif
              if (.not. ltemperature_nolog) Tobj=exp(Tobj)
              deltaT=Tobj-T0
              loc_Nus=fp_gradT_n*charac_len(iobj)/deltaT
              Nusselt(iobj)=Nusselt(iobj)+loc_Nus * lineelement(iline,iobj)* nusselt_norm
            endif
!
!  Output local Nusselt number
!
            if (it==nt) then
              if (lLoc_Nusselt_output .and. ilnTT>0) then
!
                open(unit=83,file='data/local_Nusselt.dat',position='APPEND')
!
                write(83,102) fpx, fpy, fpz, loc_Nus
!
                close(83)
              endif
            endif ! it==nt
!
!  Output local temperature
            if (it==nt) then
              if(lLoc_Temp_output .and. ilnTT>0) then
                open(unit=84,file='data/local_Temperature.dat',position='APPEND')
!
                if (bc_thermo_type=='Neumann') then
		  call dot(p%gTT(ix0-nghost,:),-nvec,fp_gradT_n)
		endif
!
                write(84,100)  fpx, fpy, fpz, Tobj, fp_gradT_n
!
                close(84)
              endif
            endif ! it==nt
!
! Output Local density
            if(lLoc_density_output) then
              open(unit=85,file='data/local_density.dat',position='APPEND')
              write(solid_cell_den,98) it-1, t
!     
              write(den_string,96)  fpx, fpy, fpz, loc_density
! 
              call safe_character_append(solid_cell_den,den_string)
              write(85,*) trim(solid_cell_den)
              close(85)
            endif
!
! Output local pressure coefficient
            if (it==nt) then
              if(lLoc_c_press_output) then
                open(unit=88,file='data/local_c_press.dat',position='APPEND')
!
                write(88,102)  fpx, fpy, fpz, c_press
!
                close(88)
              endif
            endif ! it==nt
!
! Output error_norm
            if (it==nt) then
              if (lerror_norm) then
                open(unit=86,file='data/vel_error_norm.dat',position='APPEND')

                call find_near_indeces(lower_i,upper_i,lower_j,upper_j, &
                                       lower_k,upper_k,x,y,z,fp_location)
                inear_error=(/lower_i,lower_j,lower_k/)

                if(.not. linear_interpolate(f,iux,iuz,fp_location,vel_error,inear_error,.false.)) &
                  call fatal_error('linear_interpolate','')

                write(86,94) vel_error(1),vel_error(2)

                if (ilnTT>0) then
                  open(unit=87,file='data/tem_error_norm.dat',position='APPEND')

                  if(.not. linear_interpolate(f,ilnTT,ilnTT,fp_location,T_error,inear_error,.false.)) &
                    call fatal_error('linear_interpolate','')

                  write(87,*) T_error
                  close(87)

                endif

                close(86)

              endif 
            endif ! it==nt
!
102         format(3F10.5,1X,F10.5)
100         format(3F10.5,1X,2F10.5)
98          format(1I8,1F15.8)
96          format(3F15.8)
94          format(2F10.5)
!
          endif
        endif 
      enddo ! finalize loop over all lines
      enddo ! finalize loop over all objects
!     close(444)
!     close(445)
!
      call keep_compiler_quiet(df,f)
!   
      end subroutine output_solid_cells
!************************************************************************
      endmodule Solid_Cells
