! $Id: border_profiles.f90,v 1.8 2006-08-23 16:53:30 mee Exp $ 
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lborder_profiles = .true.
!
!***************************************************************

module BorderProfiles 

  use Cparam
  use Cdata

  implicit none

  private

  include 'border_profiles.h'
!:
!  border_prof_[x-z] could be of size n[x-z], but having the same
!  length as f() (in the given dimension) gives somehow more natural code.
!
  real, dimension(mx) :: border_prof_x=1.0
  real, dimension(my) :: border_prof_y=1.0
  real, dimension(mz) :: border_prof_z=1.0
!
  contains

!***********************************************************************
    subroutine initialize_border_profiles()
!
!  Position-dependent quenching factor that multiplies rhs of pde
!  by a factor that goes gradually to zero near the boundaries.
!  border_frac is a 3-D array, separately for all three directions.
!  border_frac=1 would affect everything between center and border.
!
      use Cdata

      real, dimension(nx) :: xi
      real, dimension(ny) :: eta
      real, dimension(nz) :: zeta
      real :: border_width,lborder,uborder
!
!  x-direction
!
      border_prof_x(l1:l2)=1

      if ((border_frac_x(1)>0) .and. (.not. lperi(1))) then
        border_width=border_frac_x(1)*Lxyz(1)/2
        lborder=xyz0(1)+border_width
        xi=1-max(lborder-x(l1:l2),0.0)/border_width
        border_prof_x(l1:l2)=min(border_prof_x(l1:l2),xi**2*(3-2*xi))
      endif

      if ((border_frac_x(2)>0) .and. (.not. lperi(1))) then
        border_width=border_frac_x(2)*Lxyz(1)/2
        uborder=xyz1(1)-border_width
        xi=1-max(x(l1:l2)-uborder,0.0)/border_width
        border_prof_x(l1:l2)=min(border_prof_x(l1:l2),xi**2*(3-2*xi))
      endif
!
!  y-direction
!
      border_prof_y(m1:m2)=1

      if ((border_frac_y(1)>0) .and. (.not. lperi(2))) then
        border_width=border_frac_y(1)*Lxyz(2)/2
        lborder=xyz0(2)+border_width
        eta=1-max(lborder-y(m1:m2),0.0)/border_width
        border_prof_y(m1:m2)=min(border_prof_y(m1:m2),eta**2*(3-2*eta))
      endif

      if ((border_frac_y(2)>0) .and. (.not. lperi(2))) then
        border_width=border_frac_y(2)*Lxyz(2)/2
        uborder=xyz1(2)-border_width
        eta=1-max(y(m1:m2)-uborder,0.0)/border_width
        border_prof_y(m1:m2)=min(border_prof_y(m1:m2),eta**2*(3-2*eta))
      endif
!
!  z-direction
!
      border_prof_z(n1:n2)=1

      if ((border_frac_z(1)>0) .and. (.not. lperi(3))) then
        border_width=border_frac_z(1)*Lxyz(3)/2
        lborder=xyz0(3)+border_width
        zeta=1-max(lborder-z(n1:n2),0.0)/border_width
        border_prof_z(n1:n2)=min(border_prof_z(n1:n2),zeta**2*(3-2*zeta))
      endif

      if ((border_frac_z(2)>0) .and. (.not. lperi(3))) then
        border_width=border_frac_z(2)*Lxyz(3)/2
        uborder=xyz1(3)-border_width
        zeta=1-max(z(n1:n2)-uborder,0.0)/border_width
        border_prof_z(n1:n2)=min(border_prof_z(n1:n2),zeta**2*(3-2*zeta))
      endif
!
    endsubroutine initialize_border_profiles
!***********************************************************************
    subroutine border_driving(f,df,f_target,j)
!                                                           
!  Position-dependent driving term that attempts to drive pde
!  the variable toward some target solution on the boundary.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension(nx) :: f_target,pborder,drive_time
      integer :: j
!
      call get_drive_time(drive_time)
      call get_border(pborder)
!     
      df(l1:l2,m,n,j) = df(l1:l2,m,n,j) &
           - (f(l1:l2,m,n,j) - f_target)*pborder/drive_time 
!
    endsubroutine border_driving
!***********************************************************************
    subroutine get_border(pborder)
!
! Apply a step function that smoothly goes from zero to one on both sides. 
! In practice, means that the driving takes place 
! from r_int to r_int+2*wborder_int, and
! from r_ext-2*wborder_ext to r_ext
!
! Regions away from these limits are unaffected. 
!
! 28-Jul-06/wlad : coded
!
      use Sub, only: cubic_step
!
      real, dimension(nx),intent(out) :: pborder
      real, dimension(nx) :: rlim_mn
      integer :: i
!
      if (lcylindrical) then
         rlim_mn = rcyl_mn
      else if (lspherical) then
         rlim_mn = r_mn
      else    
         rlim_mn = x_mn
      endif
!
      pborder = cubic_step(rlim_mn,r_int,wborder_int,SHIFT=1.)&
           *(1-cubic_step(rlim_mn,r_ext,wborder_ext,SHIFT=-1.))
!         
! Set pborder to zero at the limits of the border
!
      do i=1,nx
         if ((rlim_mn(i).lt.r_int).or.&
              ((rlim_mn(i).gt.r_int+2*wborder_int).and.(rlim_mn(i).lt.r_ext-2*wborder_ext)).or. &
               (rlim_mn(i).gt.r_ext)) &
               pborder(i)=0.
      enddo
!
    endsubroutine get_border
!***********************************************************************
    subroutine get_drive_time(drive_time)
!
! This is problem-dependent, since the driving should occur in the
! typical time-scale of the problem. Currently, only the keplerian
! orbital time is implemented.
! 
! 28-Jul-06/wlad : coded
!
      real, dimension(nx),intent(out) :: drive_time
!
      drive_time = 2*pi*(rcyl_mn)**(1.5)
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
!  border_frac is a 3-D array, separately for all three directions.
!  border_frac=1 would affect everything between center and border.
!
       df(l1:l2,m,n,j) = df(l1:l2,m,n,j) &
          *border_prof_x(l1:l2)*border_prof_y(m)*border_prof_z(n)
!
    endsubroutine border_quenching
!***********************************************************************
endmodule BorderProfiles
