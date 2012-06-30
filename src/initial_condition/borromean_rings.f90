!  Initial condition (density, magnetic field, velocity) 
!  for a particular configuration of magnetic tubes.
!
!  11-nov-10/simon:
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
! ampl = amplitude of the magnetic field
! width_ring = width of the flux tubes
! [xyz]scale = scale in each dimension

  real :: ampl=1.0, width_ring=0.6, minor_axis=1.0, major_axis=2.5
  real :: xscale = 1.0, yscale = 1.0, zscale = 1.0
  character (len=labellen) :: prof='constant'
!
  namelist /initial_condition_pars/ &
      ampl,width_ring,prof,minor_axis,major_axis,xscale,yscale,zscale
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  07-oct-09/wlad: coded
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho 
!  will take care of converting it to linear 
!  density if you use ldensity_nolog
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  11-nov-10/simon: coded
!
!  Magnetic flux ring which has the form of a n-foil knot.
!  NB: Curerently this works only for one CPU. For multi CPU usage
!  initialize on one CPU and do the remeshing.
!
!  Created 2010-11-11 by Simon Candelaresi (Iomsn)
!
      use Mpicomm, only: stop_it
      use Poisson
      use Sub
      
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ellipse_param, circle_param, circle_radius
      real :: delta_ellipse_param, delta_circle_param, delta_circle_radius
      real, dimension(3) :: ellipse_pos, circle_pos, tangent, normal
      integer :: domain_width, domain_depth, domain_height
      integer :: l, j, ju, ellipse_idx
      ! The next 2 variables are used for the uncurling.
      real, dimension (nx,ny,nz,3) :: jj, tmpJ  ! This is phi for poisson.f90
      
!
!  initialize the magnetic flux tube
!
      domain_width = l2-l1; domain_depth = m2-m1; domain_height = n2-n1
!
!  Calculate the minimum step size of the curve parameters 
!  to avoid discretisation issues, like mesh points without magnetic field
!
      delta_ellipse_param = 1./max(domain_width,domain_depth,domain_height) / major_axis
      delta_circle_param = delta_ellipse_param/(width_ring/2.)
      delta_circle_radius = delta_circle_param
!
      
!
!  loop which goes around the ellipses
!
      do ellipse_idx = 1,3
        ellipse_param = 0.
!
!  loop which moves along the ellipse
!
        do
          if (ellipse_param .gt. 2.*pi) exit              
          if (ellipse_idx .eq. 1) then
            ellipse_pos(1) = major_axis*sin(ellipse_param)
            ellipse_pos(2) = minor_axis*cos(ellipse_param)
            ellipse_pos(3) = 0
            tangent(1) = major_axis*cos(ellipse_param)
            tangent(2) = -major_axis*sin(ellipse_param)
            tangent(3) = 0
          endif
          if (ellipse_idx .eq. 2) then
            ellipse_pos(1) = minor_axis*cos(ellipse_param)
            ellipse_pos(2) = 0
            ellipse_pos(3) = major_axis*sin(ellipse_param)
            tangent(1) = -major_axis*sin(ellipse_param)
            tangent(2) = 0
            tangent(3) = major_axis*cos(ellipse_param)
          endif
          if (ellipse_idx .eq. 3) then
            ellipse_pos(1) = 0
            ellipse_pos(2) = major_axis*sin(ellipse_param)
            ellipse_pos(3) = minor_axis*cos(ellipse_param)
            tangent(1) = 0
            tangent(2) = major_axis*cos(ellipse_param)
            tangent(3) = -major_axis*sin(ellipse_param)
          endif
          tangent = tangent / sqrt(tangent(1)**2+tangent(2)**2+tangent(3)**2)
!
!  Find vector which is orthonormal to tangent vector.
!
          if (abs(tangent(1)) .le. 0.5) then
            normal(1) = tangent(1)**2 - 1.0
            normal(2) = tangent(2)*tangent(1)
            normal(3) = tangent(3)*tangent(1)
          elseif (abs(tangent(2)) .le. 0.5) then
            normal(1) = tangent(1)*tangent(2)
            normal(2) = tangent(2)**2 - 1.0
            normal(3) = tangent(3)*tangent(2)
          else
            normal(1) = tangent(1)*tangent(3)
            normal(2) = tangent(2)*tangent(3)
            normal(3) = tangent(3)**2 - 1.0
          endif
!
!  normalize the normal vector
!
          normal = normal / sqrt(normal(1)**2+normal(2)**2+normal(3)**2)

          circle_radius = 0.
!
!  loop which changes the circle's radius
!
          do
            if (circle_radius .gt. width_ring/2.) exit
            circle_param = 0.
!
!  loop which goes around the circle
!
            do
              if (circle_param .gt. 2.*pi) exit
              circle_pos(1) = ellipse_pos(1) + circle_radius * &
              ((tangent(1)*tangent(1)*(1-cos(circle_param))+cos(circle_param))*normal(1) + &
              (tangent(1)*tangent(2)*(1-cos(circle_param))-tangent(3)*sin(circle_param))*normal(2) + &
              (tangent(1)*tangent(3)*(1-cos(circle_param))+tangent(2)*sin(circle_param))*normal(3))
              circle_pos(2) = ellipse_pos(2) + circle_radius * &
              ((tangent(1)*tangent(2)*(1-cos(circle_param))+tangent(3)*sin(circle_param))*normal(1) + &
              (tangent(2)*tangent(2)*(1-cos(circle_param))+cos(circle_param))*normal(2) + &
              (tangent(2)*tangent(3)*(1-cos(circle_param))-tangent(1)*sin(circle_param))*normal(3))
              circle_pos(3) = ellipse_pos(3) + circle_radius * &
              ((tangent(1)*tangent(3)*(1-cos(circle_param))-tangent(2)*sin(circle_param))*normal(1) + &
              (tangent(2)*tangent(3)*(1-cos(circle_param))+tangent(1)*sin(circle_param))*normal(2) + &
              (tangent(3)*tangent(3)*(1-cos(circle_param))+cos(circle_param))*normal(3))
!
!  Find the corresponding mesh point to this position.
!
              l = nint(circle_pos(1)/(2.*pi)*domain_width*xscale + domain_width/2.0)
              m = nint(circle_pos(2)/(2.*pi)*domain_depth*yscale + domain_depth/2.0)
              n = nint(circle_pos(3)/(2.*pi)*domain_height*zscale + domain_height/2.0)
!
!  Write the magnetic field B.
!  Note that B is written in the f-array where A is stored. This is
!  corrected further in the code.
!
              f(l,m,n,iax:iaz) = tangent*ampl
              
              circle_param = circle_param + delta_circle_param
            enddo
            circle_radius = circle_radius + delta_circle_radius
          enddo
          ellipse_param = ellipse_param + delta_ellipse_param
        enddo
      enddo
      
!
!  Transform the magnetic field into a vector potential
!

!  Compute curl(B) = J for the Poisson solver

      do m=m1,m2
         do n=n1,n2
            call curl(f,iaa,jj(:,m-nghost,n-nghost,:))
         enddo
      enddo
      tmpJ = -jj
!  Use the Poisson solver to solve \nabla^2 A = -J for A
      do j=1,3
        call inverse_laplacian(f,tmpJ(:,:,:,j))
      enddo
      
!  Overwrite the f-array with the correct vector potential A
      do j=1,3
          ju=iaa-1+j
          f(l1:l2,m1:m2,n1:n2,ju) = tmpJ(:,:,:,j)
      enddo
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_pars(unit,iostat)
!
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=initial_condition_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=initial_condition_pars,ERR=99)
      endif
      print*, 'nfoil: read initial conditions'
!
99    return
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!     
      integer, intent(in) :: unit
!
      write(unit,NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
  endmodule InitialCondition
