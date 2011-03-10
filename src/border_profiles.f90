! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lborder_profiles = .true.
!
! PENCILS PROVIDED rborder_mn
!
!***************************************************************
module BorderProfiles
!
  use Cparam
  use Cdata
!
  implicit none
!
  private
!
  include 'border_profiles.h'
!
!  border_prof_[x-z] could be of size n[x-z], but having the same
!  length as f() (in the given dimension) gives somehow more natural code.
!
  real, dimension(mx) :: border_prof_x=1.0
  real, dimension(my) :: border_prof_y=1.0
  real, dimension(mz) :: border_prof_z=1.0
!
! WL: Ideally,this 4D array f_init should be allocatable, since it 
!     is only used in specific conditions in the code (only when the 
!     border profile chosen is 'initial-condition').
!
  real, dimension(mx,my,mz,mvar) :: f_init
!
  logical :: lborder_driving=.false.
  logical :: lborder_quenching=.false.
!
  contains
!***********************************************************************
    subroutine initialize_border_profiles()
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
      use Cdata
      use Messages
!
      real, dimension(nx) :: xi
      real, dimension(ny) :: eta
      real, dimension(nz) :: zeta
      real :: border_width, lborder, uborder
      integer :: l
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
!  Write border profiles to file.
!
      open(1,file=trim(directory_snap)//'/border_prof_x.dat')
        do l=1,mx
          write(1,'(2f15.6)') x(l), border_prof_x(l)
        enddo
      close(1)
!
      open(1,file=trim(directory_snap)//'/border_prof_y.dat')
        do m=1,my
          write(1,'(2f15.6)') y(m), border_prof_y(m)
        enddo
      close(1)
!
      open(1,file=trim(directory_snap)//'/border_prof_z.dat')
        do n=1,mz
          write(1,'(2f15.6)') z(n), border_prof_z(n)
        enddo
      close(1)
!
!  Switch border quenching on if any border frac is non-zero
!
      if (any(border_frac_x/=0).or.&
          any(border_frac_y/=0).or.&
          any(border_frac_z/=0)) &
        lborder_quenching=.true.
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
      character (len=labellen) :: border_var
      logical :: lread=.true.,exists
!
      lborder_driving=.true.
!
!  Check if there is a variable requesting initial condition as border
!
      if (lread.and.border_var=='initial-condition') then         
!
!  Check if VAR0 exists
!
        inquire(FILE=trim(directory_snap)//'/VAR0',EXIST=exists)
        if (.not.exists) then
          print*,'VAR0 for iproc=',iproc,'not found!'
          print*,'Check lwrite_ic in start.in. If it'
          print*,'is false, re-run start.csh with lwrite_ic=T' 
          STOP 1
        endif
!
!  Read the date into an initial condition array f_init that will be saved
!
        open(1,FILE=trim(directory_snap)//'/VAR0',FORM='unformatted')
        read(1) f_init
        close(1)
!
!  Switch lread to false so it does not read VAR0 multiple times
!
        lread=.false.
!
      endif
!
    endsubroutine request_border_driving
!***********************************************************************
    subroutine pencil_criteria_borderprofiles()
!
!  All pencils that this module depends on are specified here.
!
!  25-dec-06/wolf: coded
!
      use Cdata
!
      if (lborder_driving) then
        lpenc_requested(i_rborder_mn)=.true.
        if (lcylindrical_coords.or.lcylinder_in_a_box) then
          lpenc_requested(i_rcyl_mn)=.true.
          lpenc_requested(i_rcyl_mn1)=.true.
          lpenc_requested(i_phix)=.true.
          lpenc_requested(i_phiy)=.true.
        elseif (lspherical_coords.or.lsphere_in_a_box) then
          lpenc_requested(i_r_mn)=.true.
          lpenc_requested(i_r_mn1)=.true.
        else
          lpenc_requested(i_x_mn)=.true.
        endif
      endif
!
    endsubroutine pencil_criteria_borderprofiles
!***********************************************************************
    subroutine calc_pencils_borderprofiles(f,p)
!
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      if (lpencil(i_rborder_mn)) then 
        if (lcylinder_in_a_box.or.lcylindrical_coords) then
          p%rborder_mn = p%rcyl_mn
        elseif (lsphere_in_a_box.or.lspherical_coords) then
          p%rborder_mn = p%r_mn
        else
          p%rborder_mn = p%x_mn
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_borderprofiles
!***********************************************************************
    subroutine set_border_initcond(f,ivar,fborder)
!
      use Messages, only: fatal_error
      use Sub, only: keep_compiler_quiet
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
        print*,'The system has no obvious symmetry. It is    '
        print*,'better to stop and check how you want to save'
        print*,'the initial condition for border profiles    '
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
        if (headtt.and.ip <= 6) &
             print*,'saving initial condition for ivar=',ivar
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
        if (headtt.and.ip <= 6) &
             print*,'saving initial condition for ivar=',ivar
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
      real, dimension(nx) :: f_target
      type (pencil_case) :: p
      real :: pborder,inverse_drive_time
      integer :: i,j
      logical :: lfirstcall=.true.
!
!  if r_int_border and/or r_ext_border are still set to impossible,
!  then put them equal to r_int and r_ext, respectively.
!
      if (lfirstcall) then 
        if (r_int_border==impossible) r_int_border=r_int
        if (r_ext_border==impossible) r_ext_border=r_ext
        lfirstcall=.false.
      endif
!
!  Perform "border_driving" only if r < r_int_border or r > r_ext_border, but
!  take into acount that the profile further inside on both ends.
!  Note: instead of setting r_int_border=0 (to mask out the middle),
!  put it to a negative value instead, to avoid surprises at r=0.
!
      do i=1,nx
        if ( &
            !inner stripe
             (p%rborder_mn(i)<=r_int_border+2*wborder_int).or.&
            !outer stripe
             (p%rborder_mn(i)>=r_ext_border-2*wborder_ext)) then
!        
          call get_drive_time(p,inverse_drive_time,i)
          call get_border(p,pborder,i)
          df(i+l1-1,m,n,j) = df(i+l1-1,m,n,j) &
               - (f(i+l1-1,m,n,j) - f_target(i))*pborder*inverse_drive_time
        endif
        !else do nothing
      enddo
!
    endsubroutine border_driving
!***********************************************************************
    subroutine get_border(p,pborder,i)
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
      type (pencil_case) :: p
      real :: rlim_mn
      integer :: i
!
      if (lcylinder_in_a_box.or.lcylindrical_coords) then
         rlim_mn = p%rcyl_mn(i)
      elseif (lsphere_in_a_box.or.lspherical_coords) then
         rlim_mn = p%r_mn(i)
      else
         rlim_mn = p%x_mn(i)
      endif
!
! cint = 1-step_int , cext = step_ext
! pborder = cint+cext
!
      pborder = 1-cubic_step(rlim_mn,r_int_border,wborder_int,SHIFT=1.) + &
           cubic_step(rlim_mn,r_ext_border,wborder_ext,SHIFT=-1.)
!
    endsubroutine get_border
!***********************************************************************
    subroutine get_drive_time(p,inverse_drive_time,i)
!
!  This is problem-dependent, since the driving should occur in the
!  typical time-scale of the problem. tborder can be specified as input.
!  Alternatively (tborder=0), the keplerian orbital time is used.
!
!  28-jul-06/wlad: coded
!  24-jun-09/axel: added tborder as input
!
      real, intent(out) :: inverse_drive_time
      real :: uphi
      type (pencil_case) :: p
      integer :: i
!
!  calculate orbital time
!
      if (tborder==0.) then
        if (lcylinder_in_a_box.or.lcylindrical_coords) then
          uphi=p%uu(i,1)*p%phix(i)+p%uu(i,2)*p%phiy(i)
        elseif (lspherical_coords) then
          uphi=p%uu(i,3)
        endif
        inverse_drive_time = .5*pi_1*uphi/p%rborder_mn(i)
!
!  specify tborder as input
!
      else
        inverse_drive_time=1./tborder
      endif
!
    endsubroutine get_drive_time
!***********************************************************************
    subroutine border_quenching(df,j)
!
      real, dimension (mx,my,mz,mvar) :: df
      integer :: j
!
!  Position-dependent quenching factor that multiplies rhs of pde
!  by a factor that goes gradually to zero near the boundaries.
!  border_frac_[xyz] is a 2-D array, separately for all three directions.
!  border_frac_[xyz]=1 would affect everything between center and border.
!
      if (lborder_quenching) &
          df(l1:l2,m,n,j) = df(l1:l2,m,n,j) &
          *border_prof_x(l1:l2)*border_prof_y(m)*border_prof_z(n)
!
    endsubroutine border_quenching
!***********************************************************************
endmodule BorderProfiles
