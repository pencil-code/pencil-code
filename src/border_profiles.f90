! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lborder_profiles = .true.
!
!***************************************************************
module BorderProfiles
!
  use Cdata
  use Messages
!
  implicit none
!
  private
!
  include 'border_profiles.h'
!
!  border_prof_[x-z] could be of size n[x-z], but having the same
!  length as f (in the given dimension) gives somehow more natural code.
!
  real, dimension(mx) :: border_prof_x=1.0
  real, dimension(my) :: border_prof_y=1.0
  real, dimension(mz) :: border_prof_z=1.0
  real, dimension(mx, my, mz) :: border_prof_r=1.0
  real                :: tborder1=impossible
  real                :: fraction_tborder1=impossible
  real                :: fac_sqrt_gsum1=1.0
!
  real, dimension(:,:,:,:), allocatable :: f_init
!
  logical :: lborder_driving=.false.
  logical :: lborder_quenching=.false.
!
  contains
!***********************************************************************
    subroutine initialize_border_profiles
!
!  Position-dependent quenching factor that multiplies rhs of pde
!  by a factor that goes gradually to zero near the boundaries.
!  border_frac_[xyz] is a 2-D array, separately for all three directions.
!  border_frac_[xyz]=1 would affect everything between center and border.
!
!   9-nov-09/axel: set r_int_border and r_ext_border if still impossible
!  17-jun-10/wlad: moved r_int_border and r_ext_border to border drive
!                  because r_int and r_ext are not set until after all
!                  calls to initialize are done in Register.
!
      use HDF5_IO, only: output_profile
      use SharedVariables, only: get_shared_variable
      use Messages, only: fatal_error
!
      real, dimension(nx) :: xi
      real, dimension(ny) :: eta
      real, dimension(nz) :: zeta
      real :: border_width, lborder, uborder
      integer :: l, ierr
      real, pointer :: gsum
!
!  x-direction
!
      border_prof_x(l1:l2)=1
!
      if (border_frac_x(1)>0) then
        if (lperi(1)) call fatal_error('initialize_border_profiles', &
            'must have lperi(1)=F for border profile in x')
        border_width=border_frac_x(1)*Lxyz(1)/2
        lborder=xyz0(1)+border_width
        xi=1-max(lborder-x(l1:l2),0.0)/border_width
        border_prof_x(l1:l2)=min(border_prof_x(l1:l2),xi**2*(3-2*xi))
      endif
!
      if (border_frac_x(2)>0) then
        if (lperi(1)) call fatal_error('initialize_border_profiles', &
            'must have lperi(1)=F for border profile in x')
        border_width=border_frac_x(2)*Lxyz(1)/2
        uborder=xyz1(1)-border_width
        xi=1-max(x(l1:l2)-uborder,0.0)/border_width
        border_prof_x(l1:l2)=min(border_prof_x(l1:l2),xi**2*(3-2*xi))
      endif
!
!  y-direction
!
      border_prof_y(m1:m2)=1
!
      if (border_frac_y(1)>0) then
        if (lperi(2)) call fatal_error('initialize_border_profiles', &
            'must have lperi(2)=F for border profile in y')
        border_width=border_frac_y(1)*Lxyz(2)/2
        lborder=xyz0(2)+border_width
        eta=1-max(lborder-y(m1:m2),0.0)/border_width
        border_prof_y(m1:m2)=min(border_prof_y(m1:m2),eta**2*(3-2*eta))
      endif
!
      if (border_frac_y(2)>0) then
        if (lperi(2)) call fatal_error('initialize_border_profiles', &
            'must have lperi(2)=F for border profile in y')
        border_width=border_frac_y(2)*Lxyz(2)/2
        uborder=xyz1(2)-border_width
        eta=1-max(y(m1:m2)-uborder,0.0)/border_width
        border_prof_y(m1:m2)=min(border_prof_y(m1:m2),eta**2*(3-2*eta))
      endif
!
!  z-direction
!
      border_prof_z(n1:n2)=1
!
      if (border_frac_z(1)>0) then
        if (lperi(3)) call fatal_error('initialize_border_profiles', &
            'must have lperi(3)=F for border profile in z')
        border_width=border_frac_z(1)*Lxyz(3)/2
        lborder=xyz0(3)+border_width
        zeta=1-max(lborder-z(n1:n2),0.0)/border_width
        border_prof_z(n1:n2)=min(border_prof_z(n1:n2),zeta**2*(3-2*zeta))
      endif
!
      if (border_frac_z(2)>0) then
        if (lperi(3)) call fatal_error('initialize_border_profiles', &
            'must have lperi(3)=F for border profile in z')
        border_width=border_frac_z(2)*Lxyz(3)/2
        uborder=xyz1(3)-border_width
        zeta=1-max(z(n1:n2)-uborder,0.0)/border_width
        border_prof_z(n1:n2)=min(border_prof_z(n1:n2),zeta**2*(3-2*zeta))
      endif
!
!  r-direction to mask a spherical object at center like inner absorbing BC
!
      border_prof_r(l1:l2,m1:m2,n1:n2)=1
      if (border_frac_r(1)>0) then
        if (.not. lcartesian_coords) call fatal_error('initialize_border_profiles', &
            'must have lcartesian_coord=T for border profile in r')
        border_width=border_frac_r(1)
        lborder=4.0+border_width ! 3*rs
        do m=m1, m2
          do n=n1, n2
            border_prof_r(l1:l2,m,n)=0.5*(1+tanh((sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)- &
                                     lborder)/border_width))
          enddo
        enddo
      endif
!
      if (border_frac_r(2)>0) then
        if (lperi(3) .or. .not. lcartesian_coords) call fatal_error('initialize_border_profiles', &
            'must have lperi(3)=F for border profile in z')
        border_width=border_frac_r(2)*maxval(Lxyz)/2
        uborder=xyz1(3)-border_width
        border_width=border_frac_r(2)*maxval(Lxyz)/60
        do m=m1, m2
          do n=n1, n2
            border_prof_r(l1:l2,m,n)=0.5*(1-tanh(sqrt((x(l1:l2)**2+y(m)**2+z(n)**2)- &
                                     uborder)/border_width))
          enddo
        enddo
      endif
!
!  Write border profiles to file.
!
      if (lmonolithic_io) then
        if (lrun) then
          call output_profile('border', x, border_prof_x, 'x', lhas_ghost=.true.) ! , loverwrite=.true.)
          call output_profile('border', y, border_prof_y, 'y', lhas_ghost=.true.) ! , loverwrite=.true.)
          call output_profile('border', z, border_prof_z, 'z', lhas_ghost=.true.) ! , loverwrite=.true.)
        endif
      else
        open(1,file=trim(directory_dist)//'/border_prof_x.dat')
          do l=1,mx
            write(1,'(2f15.6)') x(l), border_prof_x(l)
          enddo
        close(1)
!
        open(1,file=trim(directory_dist)//'/border_prof_y.dat')
          do m=1,my
            write(1,'(2f15.6)') y(m), border_prof_y(m)
          enddo
        close(1)
!
        open(1,file=trim(directory_dist)//'/border_prof_z.dat')
          do n=1,mz
            write(1,'(2f15.6)') z(n), border_prof_z(n)
          enddo
        close(1)
      endif
!
!
!  Switch border quenching on if any border frac is non-zero
!
      if (any(border_frac_x/=0).or.&
          any(border_frac_y/=0).or.&
          any(border_frac_z/=0).or.&
          any(border_frac_r/=0))     lborder_quenching=.true.
!
!  Use a fixed timescale to drive the boundary, independently of radius.
!  Else use some fraction of the local orbital time.
!
      if (tborder/=0.) then
        tborder1=1./tborder
      else
        if (lgravr) then 
          call get_shared_variable('gsum',gsum,caller='initialize_border_profiles')
          !the inverse period is the inverse of 2pi/Omega =>  Omega/2pi = sqrt(GM/2pi) / r**1.5
          !gsum = GM; the 1/r**1.5 factor will be added in inverse_period
          if (gsum/=0) then
            fac_sqrt_gsum1 = sqrt(gsum)/(2*pi)
          else
            call fatal_error("initialize_border_profiles","gsum=0")
          endif
        else
          fac_sqrt_gsum1 = 1/(2*pi)
        endif
!
        fraction_tborder1=1./fraction_tborder
!
      endif
!
      if (lmeridional_border_drive.and..not.lspherical_coords) &
           call fatal_error("initialize_border_profiles", &
           "drive meridional borders for spherical coords only")
!
    endsubroutine initialize_border_profiles
!***********************************************************************
    subroutine request_border_driving(border_var)
!
!  Tell the BorderProfiles subroutine that we need border driving.
!  Used for requesting the right pencils. Also save the initial
!  conditions in case that is the border profile to be used.
!
!  25-aug-09/anders: coded
!  06-mar-11/wlad  : added IC functionality
!
      use IO, only: input_snap, input_snap_finalize
!
      character (len=labellen) :: border_var
      logical, save :: lread=.true.
      real :: t_old
!
      lborder_driving=.true.
!
!  Check if there is a variable requesting initial condition as border
!
      if (lread .and. (border_var=='initial-condition')) then

        if (.not.allocated(f_init)) allocate(f_init(mx,my,mz,mvar))
!
!  Read the data into an initial condition array f_init that will be saved
!
        t_old = t
        if (.not.allocated(f_init)) allocate(f_init(mx,my,mz,mvar))
        call input_snap ('VAR0', f_init, mvar, mode=0)
        if (t_old /= t) then
           if (lroot) then
              print*,"You have requested to read VAR0 to set the border profile."
              print*,"This action has modified the time stamp of the simulation."
              print*,"This is probably not what you want. Likely you are starting,"
              print*,"from a previous snapshot with less physics, and set ireset_start=0 in start.in."
              print*,"You probably want to reset it to the default value in run.in."
              print*,"Do it by setting ireset_start=2 in run_pars."
           endif
           call fatal_error("request_border_driving","")
        endif
        call input_snap_finalize
        lread=.false.
!
      endif
!
    endsubroutine request_border_driving
!***********************************************************************
    subroutine pencil_criteria_borderprofiles
!
!  All pencils that this module depends on are specified here.
!
!  25-dec-06/wolf: coded
!
      if (lborder_driving) then

        if (tborder==0.) lpenc_requested(i_uu)=.true.
        if (lcylinder_in_a_box.or.lcylindrical_coords) then
          lpenc_requested(i_rcyl_mn)=.true.
          if (tborder==0.) then
            lpenc_requested(i_phix)=.true.
            lpenc_requested(i_phiy)=.true.
          endif
        elseif (lsphere_in_a_box.or.lspherical_coords) then
          lpenc_requested(i_r_mn)=.true.
        else
          lpenc_requested(i_x_mn)=.true.
        endif
!
!  If r_int_border and/or r_ext_border are still set to impossible,
!  then put them equal to r_int and r_ext, respectively.
!  Can't be done in initialize_borderprofiles as any initialize routine can set r_[int|ext].
!
        if (r_int_border==impossible) r_int_border=r_int
        if (r_ext_border==impossible) r_ext_border=r_ext
!
      endif
!
    endsubroutine pencil_criteria_borderprofiles
!***********************************************************************
    subroutine calc_pencils_borderprofiles(f,p)
!
      use General, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_borderprofiles
!***********************************************************************
    subroutine set_border_initcond(f,ivar,fborder)
!
      use General, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray),intent(in) :: f
      real, dimension (nx), intent(out) :: fborder
      integer,intent(in) :: ivar
!
      if (lspherical_coords.or.lcylinder_in_a_box) then
        call set_border_xy(ivar,fborder)
      elseif (lcylindrical_coords) then
        call set_border_xz(ivar,fborder)
      else
        if (lroot) then
          print*,'The system has no obvious symmetry. It is    '
          print*,'better to stop and check how you want to save'
          print*,'the initial condition for border profiles    '
        endif
        call fatal_error('set_border_initcond','')
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine set_border_initcond
!***********************************************************************
    subroutine set_border_xy(ivar,fborder)
!
!  Save the initial condition for a quantity that is
!  symmetric in the z axis. That can be a vertically
!  symmetric box in cartesian coordinates or an
!  azimuthally symmetric box in spherical coordinates.
!
!  28-apr-09/wlad: coded
!
      real, dimension (nx,ny,mvar), save :: fsave_init
      real, dimension (nx), intent(out) :: fborder
      integer,intent(in) :: ivar
!
      if (lfirst .and. it==1) then
        fsave_init(:,m-m1+1,ivar)=f_init(l1:l2,m,npoint,ivar)
        if (headtt.and.ip <= 6) print*,'saving initial condition for ivar=',ivar
      endif
!
      fborder=fsave_init(:,m-m1+1,ivar)
!
    endsubroutine set_border_xy
!***********************************************************************
    subroutine set_border_xz(ivar,fborder)
!
!  Save the initial condition for a quantity that is
!  symmetric in the y axis. An azimuthally symmetric
!  box in cylindrical coordinates, for instance.
!
!  28-apr-09/wlad: coded
! 
      real, dimension (nx,nz,mvar), save :: fsave_init
      real, dimension (nx), intent(out) :: fborder
      integer,intent(in) :: ivar
!
      if (lfirst .and. it==1) then
        fsave_init(:,n-n1+1,ivar)=f_init(l1:l2,mpoint,n,ivar)
        if (headtt.and.ip <= 6) print*,'saving initial condition for ivar=',ivar
      endif
!
      fborder=fsave_init(:,n-n1+1,ivar)
!
    endsubroutine set_border_xz
!***********************************************************************
    subroutine border_driving(f,df,p,f_target,j)
!
!  Position-dependent driving term that attempts to drive pde
!  the variable toward some target solution on the boundary.
!
!  The driving is applied in the inner stripe between
!  r_int_border and r_int_border+2*w, and in the outer stripe between
!  r_ext_border-2*w and r_ext_border, as sketched below
!
!  Radial extent of the box:
!
!   -------------------------------------------
!  | border |        untouched        | border |
!   -------------------------------------------
! r_int_border    r_int_border        ...          r_ext_border    r_ext_border
!          +2*w                      -2*w
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension(nx) :: f_target, rborder_mn
      type (pencil_case) :: p
      real :: pborder,inverse_drive_time
      integer :: i,j
      logical :: lradial, lmeridional
!
      if (lcylinder_in_a_box.or.lcylindrical_coords) then
        rborder_mn = p%rcyl_mn
      elseif (lsphere_in_a_box.or.lspherical_coords) then
        rborder_mn = p%r_mn
      else
        rborder_mn = p%x_mn
      endif
!
!  Perform "border_driving" only if r < r_int_border or r > r_ext_border, but
!  take into acount that the profile further inside on both ends.
!  Note: instead of setting r_int_border=0 (to mask out the middle),
!  put it to a negative value instead, to avoid surprises at r=0.
!
      do i=1,nx
!
! conditions for driving of the radial border
!
        lradial=(rborder_mn(i)<=r_int_border+2*wborder_int).or. &  ! inner stripe
                (rborder_mn(i)>=r_ext_border-2*wborder_ext)        ! outer stripe
!
! conditions for driving of the meridional border
!
        lmeridional=lmeridional_border_drive.and.&
               ((y(m)<=theta_lower_border+2*wborder_theta_lower).or. & ! lower stripe
                (y(m)>=theta_upper_border-2*wborder_theta_upper))      ! upper stripe
!
        if (lradial.or.lmeridional) then
          call get_drive_time(inverse_drive_time,rborder_mn(i))
          call get_border(pborder,rborder_mn(i))
          df(i+l1-1,m,n,j) = df(i+l1-1,m,n,j) - (f(i+l1-1,m,n,j)-f_target(i))*pborder*inverse_drive_time
        endif

      enddo
!
    endsubroutine border_driving
!***********************************************************************
    subroutine get_border(pborder,rlim_mn)
!
! Apply a step function that smoothly goes from zero to one on both sides.
! In practice, means that the driving takes place
! from r_int_border to r_int_border+2*wborder_int, and
! from r_ext_border-2*wborder_ext to r_ext_border
!
! Regions away from these limits are unaffected, because we use SHIFT=+/-1.
!
! 28-Jul-06/wlad : coded
!
      use Sub, only: cubic_step
!
      real, intent(out) :: pborder
      real :: rlim_mn
!
! cint = 1-step_int , cext = step_ext
! pborder = cint+cext
!
      pborder = 1-cubic_step(rlim_mn,r_int_border,wborder_int,SHIFT= 1.) + &
                  cubic_step(rlim_mn,r_ext_border,wborder_ext,SHIFT=-1.)
!
      if (lmeridional_border_drive) pborder = pborder + &
           1-cubic_step(y(m),theta_lower_border,wborder_theta_lower,SHIFT= 1.) + &
             cubic_step(y(m),theta_upper_border,wborder_theta_upper,SHIFT=-1.)
!
    endsubroutine get_border
!***********************************************************************
    subroutine get_drive_time(inverse_drive_time,rborder_mn)
!
!  This is problem-dependent, since the driving should occur in the
!  typical time-scale of the problem. tborder can be specified as input.
!  Alternatively (tborder=0), the keplerian orbital time is used.
!
!  28-jul-06/wlad: coded
!  24-jun-09/axel: added tborder as input
!
      real, intent(out) :: inverse_drive_time
      real :: rborder_mn
      real :: inverse_period
!
!  calculate orbital time
!
      if (tborder==0.) then
!
!  If tborder is not hardcoded, use a fraction of the 
!  orbital period 2pi/Omega. This is of course specific
!  to Keplerian global disks. 
!
        inverse_period = rborder_mn**(-1.5) * fac_sqrt_gsum1
        inverse_drive_time = fraction_tborder1*inverse_period
!
!  specify tborder as input
!
      else
        inverse_drive_time=tborder1
      endif
!
    endsubroutine get_drive_time
!***********************************************************************
    subroutine border_quenching(f,df,dt_sub)
!
      use Sub, only: del6
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_sub
!
      real, dimension (nx) :: del6_fj,border_prof_pencil
      real :: border_diff=0.01
      integer :: j
!
!  Position-dependent quenching factor that multiplies rhs of pde
!  by a factor that goes gradually to zero near the boundaries.
!  border_frac_[xyz] is a 2-D array, separately for all three directions.
!  border_frac_[xyz]=1 would affect everything between center and border.
!
      if (lborder_quenching) then
        do j = 1, mvar
          do n = n1, n2
            do m = m1, m2
              border_prof_pencil=border_prof_x(l1:l2)*border_prof_y(m)*border_prof_z(n)* &
                                 border_prof_r(l1:l2,m,n)
!
              df(l1:l2,m,n,j) = df(l1:l2,m,n,j) * border_prof_pencil
!
              if (lborder_hyper_diff) then
                if (maxval(border_prof_pencil) < 1.) then
                  call del6(f,j,del6_fj,IGNOREDX=.true.)
                  df(l1:l2,m,n,j) = df(l1:l2,m,n,j)+border_diff*(1.-border_prof_pencil)*del6_fj/dt_sub
                endif
              endif
            enddo
          enddo
        enddo
      endif
!
    endsubroutine border_quenching
!***********************************************************************
   subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=10
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(tborder1,p_par(1))
    call copy_addr(fraction_tborder1,p_par(2))
    call copy_addr(fac_sqrt_gsum1,p_par(3))
   endsubroutine pushpars2c
!***********************************************************************

endmodule BorderProfiles
