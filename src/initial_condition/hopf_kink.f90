! $Id$
!
! Create an initial magnetic field of the Hopf kink type
! as describe in Finkelstein, Weil (1978) (DOI: https://doi.org/10.1007/BF00680372),
! section 2.2 C.
! Only works in Cartesian coordinates.
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
  ! ampl = amplitude of the field
  ! configuration = particular Hopf kink setup:
  !                 single = one single kink
  !                 pair = pair of kink and anti kink along the z-axis
  ! sign_pair = sign of the kinks of the pairs (1 for same, -1 for different)

  real :: ampl = 1.0
  character (len=labellen) :: configuration = 'single'
  integer :: sign_pair = -1

  namelist /initial_condition_pars/ &
      ampl, configuration, sign_pair
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize the density field.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrho
!********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize entropy.
!
      use Poisson
      use Sub

      real, dimension (mx,my,mz,mfarray) :: f      
      integer :: l, j, ju, s_idx, pair_idx, pair_final
      real :: A, r, theta, phi, rho, omega, B_r, B_theta, B_phi, integration ! auxiliary variables
      real :: z_shift
      logical :: pair
      real, dimension (400) :: s ! array for the integration
      real, dimension (nx,ny,nz,3) :: jj, tmpJ  ! This is phi for poisson.f90
!
! Prepare for a pair.
!
      pair = .false.
      z_shift = 0
      pair_final = 1
      if (configuration == 'pair') then
          pair = .true.
          z_shift = 2
          pair_final = 2
      endif

      do pair_idx=1,pair_final
        do m=m1,m2
          do n=n1,n2
            do l=l1,l2
                r = sqrt(x(l)**2 + y(m)**2 + (z(n)-z_shift)**2)
                rho = sqrt(x(l)**2 + y(m)**2)
                theta = atan2(rho, (z(n)-z_shift))
                phi = atan2(y(m), x(l))
                
                omega = 2*pi*sin(pi*r**3/2)
                
                do s_idx = 1, 400, 1
                    s(s_idx) = r/400*s_idx
                enddo
                integration = sum(sin(2*pi*sin(pi*s**3/2)/2)**2/s)*(s(2)-s(1))

                A = ampl*exp(-4*integration)

                ! Compute the magnetic field in spherical coordinates.
                B_r = A*cos(theta)
                B_theta = -A*sin(theta)*cos(omega)
                if (sign_pair == -1) then
                    B_phi = -A*sin(theta)*sin(omega)*sign(1.0, z_shift)
                else
                    B_phi = -A*sin(theta)*sin(omega)
                endif

                ! Transform the spherical field into Cartesian coordinates.
                if (pair_idx == 1) then
                  f(l,m,n,iax) = 0
                  f(l,m,n,iay) = 0
                  f(l,m,n,iaz) = 0
                endif
                f(l,m,n,iax) = f(l,m,n,iax) + B_r*sin(theta)*cos(phi) + B_theta*cos(theta)*cos(phi) - B_phi*sin(phi)
                f(l,m,n,iay) = f(l,m,n,iay) + B_r*sin(theta)*sin(phi) + B_theta*cos(theta)*sin(phi) + B_phi*cos(phi)
                f(l,m,n,iaz) = f(l,m,n,iaz) + B_r*cos(theta) - B_theta*sin(theta)
              enddo
            enddo
          enddo
          z_shift = -z_shift
        enddo
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
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_lncc(f)
!
!  Initialize the passive scalar.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lncc
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
!
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
