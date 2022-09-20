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

  ! Note that spine_resolution refers to the number of spine points per ellipse. The total
  ! number of spine points will be three times that number.
  !integer :: spine_resolution = 200

  real :: ampl=1.0, width_ring_a=0.2, width_ring_b = 0.2,width_ring_c=0.2
  real :: radius=1
  real :: xscale = 1.0, yscale = 1.0, zscale = 1.0
  real :: twist_a = 0.0, twist_b = 0.0, twist_c = 0.0
  real :: xshift = 0, yshift = 0, zshift = 0.0
  real :: orientation_a=1, orientation_b=1, orientation_c=1
  character (len=labellen) :: prof='constant'
!
  namelist /initial_condition_pars/ &
      ampl,width_ring_a, width_ring_b,width_ring_c,&
      prof,xscale,yscale,zscale,&
      twist_a,twist_b,twist_c,xshift,yshift,zshift,radius,&
      orientation_a,orientation_b,orientation_c
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
      real, dimension (:), allocatable :: distance_to_spine
      real :: xspacing, yspacing, zspacing,delta_torus_param,torus_param
      real, dimension (3) :: origin, rho , B_0, c_0 , twist_vector, p
      real :: curvature_factor, rho_len , scaling , B_0_strength , twist , width_ring
      !real, dimension( 200 , 3):: spine, tangent, curvature, connection_vec
      real, dimension(: , :), allocatable :: spine, tangent, curvature, connection_vec
      integer :: l, j,  ju ,spine_index, min_index
      integer :: domain_width, domain_depth, domain_height
      ! The next variable is the number of spine points in the entire knot. spine_resolution is 
      !  the number of spine points per ellipse. So full_spine_resolution = 3*spine_resolution
      integer :: full_spine_resolution
      ! The next 2 variables are used for the uncurling.
      real, dimension (nx,ny,nz,3) :: jj, tmpJ  ! This is phi for poisson.f90
      integer :: spine_resolution = 200
      

      full_spine_resolution = 3*spine_resolution
      allocate(spine(full_spine_resolution,3), tangent(full_spine_resolution,3), &
               curvature(full_spine_resolution,3), connection_vec(full_spine_resolution,3), &
               distance_to_spine(full_spine_resolution))
      
!
!  initialize the magnetic flux tube
!
      domain_width = l2-l1; domain_depth = m2-m1; domain_height = n2-n1
      delta_torus_param = 2.*pi / real(spine_resolution)
!
      xspacing = 2.*pi / real(domain_width)
      yspacing = 2.*pi / real(domain_depth)
      zspacing = 2.*pi / real(domain_height)
      origin = (/ -real(domain_width - 1)*xspacing/2., &
                  -real(domain_depth - 1)*yspacing/2., -real(domain_height - 1)*zspacing/2. /)

!  We will now loop over ellipses a, b, c separately
!  Ellipse a has a major axis parallel to the x-axis and minor axis parallel to y
!  Ellipse b has a major axis parallel to the z-axis and minor axis parallel to x
!  Ellipse c has a major axis parallel to the y-axis and minor axis parallel to z
!  All ellipses are stored in the same array, with indeces as follows:
!   a: 1 to spine_resolution
!   b: spine_resolution+1 to 2*spine_resolution
!   c: 2*spine_resolution+1 to 3*spine_resolution

      do spine_index = 1, spine_resolution
            torus_param = real(spine_index)*delta_torus_param
    
    !  Loop over torus a
            spine(spine_index,1) = radius*cos(torus_param)
            spine(spine_index,2) = 0
            spine(spine_index,3) = radius*(sin(torus_param)-1.33)
    
            tangent(spine_index,1) = -radius*sin(torus_param)*orientation_a
            tangent(spine_index,2) = 0
            tangent(spine_index,3) = radius*cos(torus_param)*orientation_a
    
            curvature(spine_index,1) = -radius*cos(torus_param)
            curvature(spine_index,2) = 0
            curvature(spine_index,3) = -radius*sin(torus_param)
    
    !  Loop over torus b
            spine(spine_index + spine_resolution,1) = 0
            spine(spine_index + spine_resolution,2) = radius*cos(torus_param)
            spine(spine_index + spine_resolution,3) = radius*sin(torus_param)
    
            tangent(spine_index + spine_resolution,1) = 0
            tangent(spine_index + spine_resolution,2) = -radius*sin(torus_param)*orientation_b
            tangent(spine_index + spine_resolution,3) = radius*cos(torus_param)*orientation_b
    
            curvature(spine_index + spine_resolution,1) = 0
            curvature(spine_index + spine_resolution,2) = -radius*cos(torus_param)
            curvature(spine_index + spine_resolution,3) = -radius*sin(torus_param)
    
    !  Loop over torus c
            spine(spine_index + spine_resolution*2,1) = radius*cos(torus_param)
            spine(spine_index + spine_resolution*2,2) = 0
            spine(spine_index + spine_resolution*2,3) = radius*(sin(torus_param)+1.33)
    
            tangent(spine_index + spine_resolution*2,1) = -radius*sin(torus_param)*orientation_c
            tangent(spine_index + spine_resolution*2,2) = 0
            tangent(spine_index + spine_resolution*2,3) = radius*cos(torus_param)*orientation_c
    
            curvature(spine_index + spine_resolution*2,1) = -radius*cos(torus_param)
            curvature(spine_index + spine_resolution*2,2) = 0
            curvature(spine_index + spine_resolution*2,3) = -radius*sin(torus_param)
          enddo
      tangent = tangent / radius
! Variable names:
!   p = physical position of point in mesh
!   s = spine
!   t = tangent
!   c = curvature
!   rho = min(p-s); rho_len = vec_len(rho)

! B_0, c_0 are magnetic field (i.e. tangent) and curvature at closest pt on spine respectively
      do l=1,mx
        do m=1,my
          do n=1,mz
            ! p = (x,y,z)
            ! p is the physical position of the grid point
            p(1:3) = (/ x(l), y(m), z(n) /)

            ! rho_len = min(vector_len(p-s))
            ! distances from every point of the spine to p
            do spine_index=1,full_spine_resolution
              connection_vec(spine_index,:) =  p - spine(spine_index,:)
              distance_to_spine(spine_index) = sqrt( connection_vec(spine_index,1)**2. + &
              connection_vec(spine_index,2)**2. + connection_vec(spine_index,3)**2. )
            enddo

            ! let S_p be the spine point closest to p
            ! min_index is the locaton of S_p in the array "spine"
            min_index = minloc(distance_to_spine,dim=1)
            if (min_index .lt. spine_resolution) then
                  twist = twist_a
                  width_ring = width_ring_a
            else if (min_index .lt. 2*spine_resolution) then
                  twist = twist_b
                  width_ring = width_ring_b
            else
                  twist = twist_c
                  width_ring = width_ring_c
            endif
            ! rho is the vector from S_p to p, 
            rho(:) = connection_vec(min_index,:)
            rho_len = sqrt(rho(1)**2 + rho(2)**2 + rho(3)**2)

            ! B_0 and c_0 are tangent and curvature at corresp spine segment
            ! B_0_strength will be used to normalise B_0 and twist_vector below
            B_0(:) = tangent(min_index,:)
            B_0_strength = sqrt( B_0(1)**2 + B_0(2)**2 + B_0(3)**2)
            c_0(:) = curvature(min_index,:)

            ! Curvature based scaling by (1-dot(c,rho))
            ! i.e. by ratio of radius of curvature at p vs S_p
            curvature_factor = 1. - dot_product(c_0,rho)

            ! Twist vector perpendicular to the spine and rho
            ! twist_vector = vector_prod(B_0 , rho) * twist
            twist_vector(1) = B_0(2)*rho(3) - B_0(3)*rho(2)
            twist_vector(2) = B_0(3)*rho(1) - B_0(1)*rho(3)
            twist_vector(3) = B_0(1)*rho(2) - B_0(2)*rho(1)
            twist_vector = twist_vector*twist/rho_len

            ! Scaling function using tanh and dist from spine
            ! Note that width_ring = 2* radius_ring, so mult by 10 instead of 5
            scaling = tanh(-60.*rho_len/width_ring + 29.3)/2. + 0.5

            ! Calc magnetic field at p
            f(l,m,n,iax:iaz) = (B_0 + twist_vector)*scaling*ampl*curvature_factor/B_0_strength
          enddo
        enddo
      enddo

      deallocate(spine, tangent, curvature, connection_vec, distance_to_spine)  


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
