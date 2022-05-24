!  Initial condition (density, magnetic field, velocity) 
!  for a particular configuration of a magnetic tube. 
!
!  04-aug-10/simon: created from trefoil.f90
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
! [xyz]scale = scale in each dimension
! [xyz]shift = shift in each dimension in a 2*pi box
! twist = B_phi/B_0
    
  real :: ampl=1.0,width_ring=0.3, C=2.0, D=2.0
!   real :: xscale = 1.0, yscale = 1.0, zscale = 1.0
!   real :: xshift = 0.2, yshift = 0.5, zshift = 0.0
  real :: twist = 0
  integer :: spine_resolution = 200, n_foil = 3
  character (len=labellen) :: prof='constant'
!
  namelist /initial_condition_pars/ &
      ampl,width_ring,prof,C,D,n_foil,twist,spine_resolution
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
      use Mpicomm, only: stop_it
      use Poisson
      use Sub
          
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:), allocatable :: distance_to_spine
      real :: xspacing, yspacing, zspacing,delta_knot_param,dxyz,knot_param
      real, dimension (3) :: origin, rho , B_0, c_0 , twist_vector, p
      real :: curvature_factor, rho_len , scaling , B_0_strength
      !real, dimension( 200 , 3):: spine, tangent, curvature, connection_vec
      real, dimension(: , :), allocatable :: spine, tangent, curvature, connection_vec
      integer :: l, j, ju ,spine_index, min_index
      integer :: domain_width, domain_depth, domain_height
      ! The next 2 variables are used for the uncurling.
      real, dimension (nx,ny,nz,3) :: jj, tmpJ  ! This is phi for poisson.f90

      allocate(spine(spine_resolution,3), tangent(spine_resolution,3), &
               curvature(spine_resolution,3), connection_vec(spine_resolution,3), &
               distance_to_spine(spine_resolution))
!
!  initialize the magnetic flux tube
!
      delta_knot_param = 2.*pi / real(spine_resolution)
      domain_width = l2-l1; domain_depth = m2-m1; domain_height = (n2-n1)
      xspacing = 2.*pi / real(domain_width)
      yspacing = 2.*pi / real(domain_depth)
      zspacing = 2.*pi / real(domain_height)
      origin = (/ -real(domain_width - 1)*xspacing/2., &
                  -real(domain_depth - 1)*yspacing/2., -real(domain_height - 1)*zspacing/2. /)

!
!  loop which moves along the n-foil knot
!
!  Position along the central spine.
      do spine_index = 1, spine_resolution
        knot_param = real(spine_index)*delta_knot_param
        spine(spine_index,1) = (C+sin(knot_param*n_foil))*sin(knot_param*(n_foil-1.))
        spine(spine_index,2) = (C+sin(knot_param*n_foil))*cos(knot_param*(n_foil-1.))
        spine(spine_index,3) = D*cos(knot_param*n_foil)
!  Tangent to the trajectory and direction of the magnetic field on the spine
        tangent(spine_index,1) = n_foil*cos(knot_param*n_foil)*sin(knot_param*(n_foil-1)) + &
                                 (n_foil-1)*(C+sin(knot_param*n_foil))*cos(knot_param*(n_foil-1))
        tangent(spine_index,2) = n_foil*cos(knot_param*n_foil)*cos(knot_param*(n_foil-1)) - &
                                 (n_foil-1)*(C+sin(knot_param*n_foil))*sin(knot_param*(n_foil-1))
        tangent(spine_index,3) = -D*n_foil*sin(knot_param*n_foil)
    
!  Normal component using the length of the curve as parametrization.
        curvature(spine_index,1) = (n_foil*sin(knot_param*(n_foil - 1.0d0))*cos(knot_param*n_foil) + (C + &
                        sin(knot_param*n_foil))*(n_foil - 1)*cos(knot_param*(n_foil - &
                        1.0d0)))*(-D**2*n_foil**3*sin(knot_param*n_foil)*cos(knot_param* &
                        n_foil) - 1.0d0/2.0d0*(n_foil*sin(knot_param*(n_foil - 1.0d0))* &
                        cos(knot_param*n_foil) + (C + sin(knot_param*n_foil))*(n_foil - 1 &
                        )*cos(knot_param*(n_foil - 1.0d0)))*(-2*n_foil**2*sin(knot_param* &
                        n_foil)*sin(knot_param*(n_foil - 1.0d0)) + 4*n_foil*(n_foil - 1)* &
                        cos(knot_param*n_foil)*cos(knot_param*(n_foil - 1.0d0)) - 2*(C + &
                        sin(knot_param*n_foil))*(n_foil - 1)**2*sin(knot_param*(n_foil - &
                        1.0d0))) - 1.0d0/2.0d0*(n_foil*cos(knot_param*n_foil)*cos( &
                        knot_param*(n_foil - 1.0d0)) - (C + sin(knot_param*n_foil))*( &
                        n_foil - 1)*sin(knot_param*(n_foil - 1.0d0)))*(-2*n_foil**2*sin( &
                        knot_param*n_foil)*cos(knot_param*(n_foil - 1.0d0)) - 4*n_foil*( &
                        n_foil - 1)*sin(knot_param*(n_foil - 1.0d0))*cos(knot_param* &
                        n_foil) - 2*(C + sin(knot_param*n_foil))*(n_foil - 1)**2*cos( &
                        knot_param*(n_foil - 1.0d0))))/(D**2*n_foil**2*sin(knot_param* &
                        n_foil)**2 + (n_foil*sin(knot_param*(n_foil - 1.0d0))*cos( &
                        knot_param*n_foil) + (C + sin(knot_param*n_foil))*(n_foil - 1)* &
                        cos(knot_param*(n_foil - 1.0d0)))**2 + (n_foil*cos(knot_param* &
                        n_foil)*cos(knot_param*(n_foil - 1.0d0)) - (C + sin(knot_param* &
                        n_foil))*(n_foil - 1)*sin(knot_param*(n_foil - 1.0d0)))**2)**( &
                        3.0d0/2.0d0) + (-n_foil**2*sin(knot_param*n_foil)*sin(knot_param* &
                        (n_foil - 1.0d0)) + 2*n_foil*(n_foil - 1)*cos(knot_param*n_foil)* &
                        cos(knot_param*(n_foil - 1.0d0)) - (C + sin(knot_param*n_foil))*( &
                        n_foil - 1)**2*sin(knot_param*(n_foil - 1.0d0)))/sqrt(D**2*n_foil &
                        **2*sin(knot_param*n_foil)**2 + (n_foil*sin(knot_param*(n_foil - &
                        1.0d0))*cos(knot_param*n_foil) + (C + sin(knot_param*n_foil))*( &
                        n_foil - 1)*cos(knot_param*(n_foil - 1.0d0)))**2 + (n_foil*cos( &
                        knot_param*n_foil)*cos(knot_param*(n_foil - 1.0d0)) - (C + sin( &
                        knot_param*n_foil))*(n_foil - 1)*sin(knot_param*(n_foil - 1.0d0 &
                        )))**2)
        curvature(spine_index,2) = (n_foil*cos(knot_param*n_foil)*cos(knot_param*(n_foil - 1.0d0)) - (C + &
                        sin(knot_param*n_foil))*(n_foil - 1)*sin(knot_param*(n_foil - &
                        1.0d0)))*(-D**2*n_foil**3*sin(knot_param*n_foil)*cos(knot_param* &
                        n_foil) - 1.0d0/2.0d0*(n_foil*sin(knot_param*(n_foil - 1.0d0))* &
                        cos(knot_param*n_foil) + (C + sin(knot_param*n_foil))*(n_foil - 1 &
                        )*cos(knot_param*(n_foil - 1.0d0)))*(-2*n_foil**2*sin(knot_param* &
                        n_foil)*sin(knot_param*(n_foil - 1.0d0)) + 4*n_foil*(n_foil - 1)* &
                        cos(knot_param*n_foil)*cos(knot_param*(n_foil - 1.0d0)) - 2*(C + &
                        sin(knot_param*n_foil))*(n_foil - 1)**2*sin(knot_param*(n_foil - &
                        1.0d0))) - 1.0d0/2.0d0*(n_foil*cos(knot_param*n_foil)*cos( &
                        knot_param*(n_foil - 1.0d0)) - (C + sin(knot_param*n_foil))*( &
                        n_foil - 1)*sin(knot_param*(n_foil - 1.0d0)))*(-2*n_foil**2*sin( &
                        knot_param*n_foil)*cos(knot_param*(n_foil - 1.0d0)) - 4*n_foil*( &
                        n_foil - 1)*sin(knot_param*(n_foil - 1.0d0))*cos(knot_param* &
                        n_foil) - 2*(C + sin(knot_param*n_foil))*(n_foil - 1)**2*cos( &
                        knot_param*(n_foil - 1.0d0))))/(D**2*n_foil**2*sin(knot_param* &
                        n_foil)**2 + (n_foil*sin(knot_param*(n_foil - 1.0d0))*cos( &
                        knot_param*n_foil) + (C + sin(knot_param*n_foil))*(n_foil - 1)* &
                        cos(knot_param*(n_foil - 1.0d0)))**2 + (n_foil*cos(knot_param* &
                        n_foil)*cos(knot_param*(n_foil - 1.0d0)) - (C + sin(knot_param* &
                        n_foil))*(n_foil - 1)*sin(knot_param*(n_foil - 1.0d0)))**2)**( &
                        3.0d0/2.0d0) + (-n_foil**2*sin(knot_param*n_foil)*cos(knot_param* &
                        (n_foil - 1.0d0)) - 2*n_foil*(n_foil - 1)*sin(knot_param*(n_foil &
                        - 1.0d0))*cos(knot_param*n_foil) - (C + sin(knot_param*n_foil))*( &
                        n_foil - 1)**2*cos(knot_param*(n_foil - 1.0d0)))/sqrt(D**2*n_foil &
                        **2*sin(knot_param*n_foil)**2 + (n_foil*sin(knot_param*(n_foil - &
                        1.0d0))*cos(knot_param*n_foil) + (C + sin(knot_param*n_foil))*( &
                        n_foil - 1)*cos(knot_param*(n_foil - 1.0d0)))**2 + (n_foil*cos( &
                        knot_param*n_foil)*cos(knot_param*(n_foil - 1.0d0)) - (C + sin( &
                        knot_param*n_foil))*(n_foil - 1)*sin(knot_param*(n_foil - 1.0d0 &
                        )))**2)
        curvature(spine_index,3) = -D*n_foil**2*cos(knot_param*n_foil)/sqrt(D**2*n_foil**2*sin(knot_param* &
                        n_foil)**2 + (n_foil*sin(knot_param*(n_foil - 1.0d0))*cos( &
                        knot_param*n_foil) + (C + sin(knot_param*n_foil))*(n_foil - 1)* &
                        cos(knot_param*(n_foil - 1.0d0)))**2 + (n_foil*cos(knot_param* &
                        n_foil)*cos(knot_param*(n_foil - 1.0d0)) - (C + sin(knot_param* &
                        n_foil))*(n_foil - 1)*sin(knot_param*(n_foil - 1.0d0)))**2) - D* &
                        n_foil*(-D**2*n_foil**3*sin(knot_param*n_foil)*cos(knot_param* &
                        n_foil) - 1.0d0/2.0d0*(n_foil*sin(knot_param*(n_foil - 1.0d0))* &
                        cos(knot_param*n_foil) + (C + sin(knot_param*n_foil))*(n_foil - 1 &
                        )*cos(knot_param*(n_foil - 1.0d0)))*(-2*n_foil**2*sin(knot_param* &
                        n_foil)*sin(knot_param*(n_foil - 1.0d0)) + 4*n_foil*(n_foil - 1)* &
                        cos(knot_param*n_foil)*cos(knot_param*(n_foil - 1.0d0)) - 2*(C + &
                        sin(knot_param*n_foil))*(n_foil - 1)**2*sin(knot_param*(n_foil - &
                        1.0d0))) - 1.0d0/2.0d0*(n_foil*cos(knot_param*n_foil)*cos( &
                        knot_param*(n_foil - 1.0d0)) - (C + sin(knot_param*n_foil))*( &
                        n_foil - 1)*sin(knot_param*(n_foil - 1.0d0)))*(-2*n_foil**2*sin( &
                        knot_param*n_foil)*cos(knot_param*(n_foil - 1.0d0)) - 4*n_foil*( &
                        n_foil - 1)*sin(knot_param*(n_foil - 1.0d0))*cos(knot_param* &
                        n_foil) - 2*(C + sin(knot_param*n_foil))*(n_foil - 1)**2*cos( &
                        knot_param*(n_foil - 1.0d0))))*sin(knot_param*n_foil)/(D**2* &
                        n_foil**2*sin(knot_param*n_foil)**2 + (n_foil*sin(knot_param*( &
                        n_foil - 1.0d0))*cos(knot_param*n_foil) + (C + sin(knot_param* &
                        n_foil))*(n_foil - 1)*cos(knot_param*(n_foil - 1.0d0)))**2 + ( &
                        n_foil*cos(knot_param*n_foil)*cos(knot_param*(n_foil - 1.0d0)) - &
                        (C + sin(knot_param*n_foil))*(n_foil - 1)*sin(knot_param*(n_foil &
                        - 1.0d0)))**2)**(3.0d0/2.0d0)
      enddo

! Variable names:
!   p = physical position of point in mesh
!   s = spine
!   t = tangent
!   c = curvature
!   rho = min(p-s); rho_len = vec_len(rho)

! B_0, c_0 are magnetic field (i.e. tangent) and curvature at closest pt on spine respectively
!      do l=l1,l2
!        do m=m1,m2
!          do n=n1,n2
      do l=1,mx
        do m=1,my
          do n=1,mz
            ! p = (x,y,z)
            ! p is the physical position of the grid point
            p(1:3) = (/ x(l), y(m), z(n) /)

            ! rho_len = min(vector_len(p-s))
            ! distances from every point of the spine to p
            do spine_index=1,spine_resolution
              connection_vec(spine_index,:) =  p - spine(spine_index,:)
              distance_to_spine(spine_index) = sqrt( connection_vec(spine_index,1)**2. + &
              connection_vec(spine_index,2)**2. + connection_vec(spine_index,3)**2. )
            enddo

            ! let S_p be the spine point closest to p
            ! min_index is the locaton of S_p in the array "spine"
            min_index = minloc(distance_to_spine,dim=1)
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
            twist_vector = twist_vector*twist

            ! Scaling function using tanh and dist from spine
            ! Note that width_ring = 2* radius_ring, so mult by 10 instead of 5
            scaling = tanh(-20.*rho_len/width_ring + 9.3)/2. + 0.5

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
    
