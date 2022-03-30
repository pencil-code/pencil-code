!  Initial condition (density, magnetic field, velocity) 
!  for a particular configuration of a magnetic tube. 
!
!  17-jan-11/simon: created from nfoil.f90
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
! width_ring = width of the flux tube
! C = offset from the center of the knot, should be > 1
! D = extention coefficient for the z direction
! phase = relative phase of the z-variation
! twist = internal twist of the flux tube
! [xyz]scale = scale in each dimension
! [xyz]shift = shift in each dimension in a 2*pi box

  real :: ampl=1.0, width_ring=0.3, C=2.0, D=1.5, phase=4./3., twist=0
  real :: xscale=1.0, yscale=1.0, zscale=1.0
  real :: xshift=0.0, yshift=0.0, zshift=0.0
  character (len=labellen) :: prof='constant'
!
  namelist /initial_condition_pars/ &
      ampl,width_ring,prof,C,D,xscale,yscale,zscale,xshift,yshift,zshift,phase,twist
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
!  04-aug-10/simon: coded
!
!  Magnetic flux ring which has the form of a n-foil knot.
!
      use Mpicomm, only: stop_it
      use Poisson
      use Sub
      
      real, dimension (mx,my,mz,mfarray) :: f
      real :: knot_param, circle_param, circle_radius
      real :: delta_knot_param, delta_circle_param, delta_circle_radius
      real, dimension(3) :: knot_pos, circle_pos, tangent, normal, twist_vector
      integer :: domain_width, domain_depth, domain_height
      integer :: l, j, ju
      real :: dist
      ! The next 2 variables are used for the uncurling.
      real, dimension (nx,ny,nz,3) :: jj, tmpJ  ! This is phi for poisson.f90
      
!
!  initialize the magnetic flux tube
!
      domain_width = l2-l1; domain_depth = m2-m1; domain_height = n2-n1
!
!  Calculate the minimum step size of the curve parameters 
!  to avoid discretation issues, like mesh points without magnetic field
!
      delta_knot_param = .5/max(domain_width,domain_depth,domain_height)
      delta_circle_param = delta_knot_param/(width_ring/2.)
      delta_circle_radius = delta_circle_param
!
      knot_param = 0.
!
!  loop which moves along the n-foil knot
!
      do
        if (knot_param .gt. 2.*pi) exit
        knot_pos(1) = (C+sin(knot_param*4))*sin(knot_param*3)
        knot_pos(2) = (C+sin(knot_param*4))*cos(knot_param*3)
        knot_pos(3) = D*cos(8.*knot_param-phase*pi)
        tangent(1) = 4*cos(knot_param*4)*sin(knot_param*3)+&
               3*(C+sin(knot_param*4))*cos(knot_param*(4-1))
        tangent(2) = 4*cos(knot_param*4)*cos(knot_param*3)-&
               3*(C+sin(knot_param*4))*sin(knot_param*3)
        tangent(3) = -8*D*sin(8*(knot_param-pi/6))
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
!
!  Compute the position on the finite size flux tube.
!
            circle_pos(1) = knot_pos(1) + circle_radius * &
            ((tangent(1)*tangent(1)*(1-cos(circle_param))+cos(circle_param))*normal(1) + &
            (tangent(1)*tangent(2)*(1-cos(circle_param))-tangent(3)*sin(circle_param))*normal(2) + &
            (tangent(1)*tangent(3)*(1-cos(circle_param))+tangent(2)*sin(circle_param))*normal(3))
            circle_pos(2) = knot_pos(2) + circle_radius * &
            ((tangent(1)*tangent(2)*(1-cos(circle_param))+tangent(3)*sin(circle_param))*normal(1) + &
            (tangent(2)*tangent(2)*(1-cos(circle_param))+cos(circle_param))*normal(2) + &
            (tangent(2)*tangent(3)*(1-cos(circle_param))-tangent(1)*sin(circle_param))*normal(3))
            circle_pos(3) = knot_pos(3) + circle_radius * &
            ((tangent(1)*tangent(3)*(1-cos(circle_param))-tangent(2)*sin(circle_param))*normal(1) + &
            (tangent(2)*tangent(3)*(1-cos(circle_param))+tangent(1)*sin(circle_param))*normal(2) + &
            (tangent(3)*tangent(3)*(1-cos(circle_param))+cos(circle_param))*normal(3))
!
!  Compute the twist vector for the magnetic field.
            twist_vector(1) = tangent(2)*(circle_pos(3)-knot_pos(3)) - tangent(3)*(circle_pos(2)-knot_pos(2))
            twist_vector(2) = tangent(3)*(circle_pos(1)-knot_pos(1)) - tangent(1)*(circle_pos(3)-knot_pos(3))
            twist_vector(3) = tangent(1)*(circle_pos(3)-knot_pos(3)) - tangent(2)*(circle_pos(1)-knot_pos(1))
            twist_vector = twist_vector*circle_radius*twist
            if (sqrt((circle_pos(1)-knot_pos(1))**2 + (circle_pos(2)-knot_pos(2))**2 + &
            (circle_pos(3)-knot_pos(3))**2) > 0) then
                twist_vector = twist_vector/sqrt((circle_pos(1)-knot_pos(1))**2 + &
                (circle_pos(2)-knot_pos(2))**2 + &
                (circle_pos(3)-knot_pos(3))**2)
            endif
!
!
!  Scale and shift the circle position.
!
            circle_pos(1) = circle_pos(1)*xscale + xshift
            circle_pos(2) = circle_pos(2)*yscale + yshift
            circle_pos(3) = circle_pos(3)*zscale + zshift
!
!  Find the corresponding mesh point to this position.
!
            dist = x(l2+nghost) - x(1)
            do l = 1, l2+nghost
                if (abs(x(l) - circle_pos(1)) > dist) exit
                dist = abs(x(l) - circle_pos(1))
            end do
            dist = y(m2+nghost) - y(1)
            do m = 1, m2+nghost
                if (abs(y(m) - circle_pos(2)) > dist) exit
                dist = abs(y(m) - circle_pos(2))
            end do
            dist = z(n2+nghost) - z(1)
            do n = 1, n2+nghost
                if (abs(z(n) - circle_pos(3)) > dist) exit
                dist = abs(z(n) - circle_pos(3))
            end do
!
!  Write the magnetic field B.
!  Note that B is written in the f-array where A is stored. This is
!  corrected further in the code.
!
!              f(l,m,n,iax:iaz) = tangent*ampl
            if (l >= 1 .and. l <= l2+nghost .and. m >= 1 .and. m <= m2+nghost .and. n >= 1 .and. n <= n2+nghost) then
                if (prof == 'const') then
                    f(l,m,n,iax:iaz) = (tangent + twist_vector)*ampl
                else if (prof == 'smooth') then     
                    f(l,m,n,iax:iaz) = (tangent+twist_vector)*ampl*(1-(circle_radius/(width_ring/2.))**4)**4
                endif
            endif
            circle_param = circle_param + delta_circle_param
          enddo
          circle_radius = circle_radius + delta_circle_radius
        enddo
        knot_param = knot_param + delta_knot_param
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
        call inverse_laplacian(tmpJ(:,:,:,j))
      enddo
      
!  Overwrite the f-array with the correct vector potential A
      do j=1,3
          ju=iaa-1+j
          f(l1:l2,m1:m2,n1:n2,ju) = tmpJ(:,:,:,j)
      enddo
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
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
