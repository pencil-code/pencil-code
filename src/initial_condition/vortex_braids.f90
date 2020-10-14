!  Initial condition (density, magnetic field, velocity)
!  for a particular configuration of a braided vortex field.
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
  use Mpicomm
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
! ampl = amplitude of the magnetic fie
! Omega = background rotation
! width_tube = width of the flux tube
! braid_margin = margin of outer most strands to the borders
! braid_shift_x = right shift of the braiding configuration in x-direction
! braid_shift_y = right shift of the braiding configuration in y-direction
! l_sigma = length of the twist region
! steepness = steepness of the braiding
! B_bkg = strength of the background field in z-direction
! word = sequence of the braid group
! prof = the amplitude profile across the tube
! rho0 = rho_0 from the density distribution to compensate the centrifugal force via a pressure gradient
! grad_p = boolean to activate pressure gradient to compensate for cetrifugal force in a fixed frame
!
! n_blobs = number of blobs for the blob configuration
! xc, yc, zc = position of the blobs
! blob_sgn = sign of the twist in the blob
! l_blob = length in z-direction of the blob
! blob_scale = scaling factor for the Gaussian
!
! inFile = Initial condition file from externally created field.
!
  real :: ampl = 1.0, B0 = 0.0, Omega_bkg = 0.0
  integer :: n_blobs = 0
  real, dimension (9) :: xc, yc, zc, kappa, l_blob, a_blob
  integer :: configuration = 0
  character (len=50) :: inFilePrefix='out'
  logical :: grad_p=.false.
  real :: rho0
!
  namelist /initial_condition_pars/ &
    ampl,Omega_bkg,n_blobs,xc,yc,zc,kappa,l_blob,a_blob,configuration,inFilePrefix,B0,rho0,grad_p
!
  contains
!***********************************************************************
  subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  01-july-11/simon: coded
!
!  Identify CVS/SVN version information.
!
    if (lroot) call svn_id( &
        "$Id: vortex_braid.f90,v 1.0 2017-05-31 16:43:18 iomsn Exp $")
!
  endsubroutine register_initial_condition
!***********************************************************************
  subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  31-may-17/iomsn: coded
!
    use Mpicomm, only: stop_it
    use Poisson
    use Sub
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
    real, dimension (mx,my,mz,3) :: uu
    real, dimension (9) :: rc, theta_c
    integer :: j, l
    character (len=10) :: intString
    character (len=50) :: inFile
!
!   clear the velocity field to zero
    f(:,:,:,iux:iuz) = 0.
!
!   Read the initial condition from an externally craeted file (Python routine), similar to import.f90.
    write(intString, "(I10.1)"), iproc
    write(inFile, "(A, A)") adjustr(trim(inFilePrefix)), adjustl(trim(intString))
    open(0, file = inFile, form = 'unformatted', action = 'read', access = 'stream')
    read(0) f(:,:,:,iux:iuz)
    write(*,*) 'file read'
!   Compute the velocity from the C field.
    do m=m1,m2
        do n=n1,n2
            call curl(f,iux,uu(l1:l2,m,n,:))
        enddo
    enddo

    f(:,:,:,iux:iuz) = ampl*uu

!   Add a background velocity field.
    do m=m1,m2
        do n=n1,n2
          f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + x(l1:l2)*Omega_bkg/2
        enddo
    enddo
!
  endsubroutine initial_condition_uu
!***********************************************************************
  subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho
!  will take care of converting it to linear
!  density if you use ldensity_nolog
!
!  07-sep-20/iomsn: coded
!
    use Mpicomm, only: stop_it
    use Poisson
    use Sub
    use EquationOfState, only: cs0
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    do m=m1,m2
        do n=n1,n2
            f(l1:l2,m,n,ilnrho) = log(rho0) + 0.5*Omega_bkg**2/cs0**2 * (x(l1:l2)**2 - x0**2)
        enddo
    enddo
!
    call keep_compiler_quiet(f)
!
  endsubroutine initial_condition_lnrho
!***********************************************************************
  subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  04-september-19/simon: coded
!
!  Homogeneous magnetic field compatible with periodic y-boundarier.
!
    real, dimension (mx,my,mz,mfarray) :: f
!        
    do m=m1,m2
        do n=n1,n2
            f(l1:l2,m,n,iaa+1) = B0*x(l1:l2)/2
        enddo
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
