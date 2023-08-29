! $Id$
!
! Initial condition for the hopf kink magnetic field as
! discribed by Finkelstein and Weil (1978), eqs. in 2.2 C,
! DOI: 10.1007/BF00680372.

! Set the single magnetic field or array of fields.
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
  logical :: bz_correction = .false.

  namelist /initial_condition_pars/ &
      ampl, configuration, sign_pair, bz_correction
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
      ! Variables used for the div(B) = 0correction outside the sphere.
      integer :: z_0_idx
      real :: z_0, rho_0, psi_0, B_z_0, rho_section, psi_in, rho_integral, slope
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
        ! Add the field inside the sphere.
        do m=1,my
          do n=1,mz
            do l=1,mx
                r = sqrt(x(l)**2 + y(m)**2 + (z(n)-z_shift)**2)
                rho = sqrt(x(l)**2 + y(m)**2)
                theta = atan2(rho, (z(n)-z_shift))
                phi = atan2(y(m), x(l))
                
                omega = 2*pi*sin(pi*r**3/2)
                if (r > 1) then
                    omega = 2*pi
                endif
                
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
          
          ! Add the vertical field outside the spherefor div(B) = 0 condition.          
          if (bz_correction .eqv. .true.) then
            ! Perform some preliminarycalculations on the equator
            ! We use these numbers for thecorrections outside the sphere.
            z_0_idx = Int(mz/2)
            z_0 = z(z_0_idx)
            rho_0 = sqrt(1**2 - z_0**2)
            psi_0 = 0
            ! Compute the magnetic flux though the equator part.
            do m=m1,m2
                do l=l1,l2
                    if (x(l)**2 + y(m)**2 <= rho_0**2) then
                        psi_0 = psi_0 + f(l,m,z_0_idx,iaz)*dx*dy
                    endif
                enddo
            enddo
            
            ! Compute the magnetic field at the sphere's boundary.            
            do s_idx = 1, 400, 1
                s(s_idx) = 1.0/400*s_idx
            enddo
            integration = sum(sin(2*pi*sin(pi*s**3/2)/2)**2/s)*(s(2)-s(1))
            B_z_0 = ampl*exp(-4*integration)
            
            ! Peroform the corrections for everyhorizontal slice.
            do n=n1,n2
                ! Circular radius of the sphere on this xy-slice.
                if ((1**2 - z(n)**2) >= 0) then
                    rho_section = sqrt(1**2 - z(n)**2)
                else
                    rho_section = 0.0
                endif
                
                ! Compute the flux within the section.
                psi_in = 0
                do m=m1,m2
                    do l=l1,l2
                        if (x(l)**2 + y(m)**2 <= rho_section**2) then
                            psi_in = psi_in + f(l,m,n,iaz)*dx*dy
                        endif
                    enddo
                enddo
                
                ! Adjust the magnetic field outside section by linearly increasing it.
                rho_integral = 0
                do m=m1,m2
                    do l=l1,l2
                        rho_integral = rho_integral + sqrt(x(l)**2 + y(m)**2)*dx*dy
                    enddo
                enddo
                slope = (psi_0 - B_z_0*Lx*Ly + 2*pi*rho_section**2*B_z_0) / (rho_integral - 2.0/3*pi*rho_section**3)
                do m=m1,m2
                    do l=l1,l2
                        if (x(l)**2 + y(m)**2 > rho_section**2) then
                            f(l,m,n,iaz) = (sqrt(x(l)**2 + y(m)**2) - rho_section)*slope + B_z_0
                        endif
                    enddo
                enddo
            enddo
          endif
      enddo

    ! Removed the mean magnetic field in the z-direction for the uncurling.
    write(*,*) "add mean B_z in start.in under &magnetic_init_pars as b_ext = 0, 0, ", sum(f(:,:,:,iaz))/size(f(:,:,:,iaz))
    f(:,:,:,iaz) = f(:,:,:,iaz) - sum(f(:,:,:,iaz))/size(f(:,:,:,iaz))

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
