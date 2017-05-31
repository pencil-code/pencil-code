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
! width_tube = width of the flux tube
! braid_margin = margin of outer most strands to the borders
! braid_shift_x = right shift of the braiding configuration in x-direction
! braid_shift_y = right shift of the braiding configuration in y-direction
! l_sigma = length of the twist region
! steepness = steepness of the braiding
! B_bkg = strength of the background field in z-direction
! word = sequence of the braid group
! prof = the amplitude profile across the tube
!
! n_blobs = number of blobs for the blob configuration
! xc, yc, zc = position of the blobs
! blob_sgn = sign of the twist in the blob
! l_blob = length in z-direction of the blob
! blob_scale = scaling factor for the Gaussian
!
  real :: ampl = 1.0
  integer :: n_blobs = 0
  real, dimension (9) :: xc, yc, zc, kappa, l_blob, a_blob
!
  namelist /initial_condition_pars/ &
    ampl,n_blobs,xc,yc,zc,kappa,l_blob,a_blob
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
        "$Id: braids.f90,v 1.9 2011-08-02 16:43:18 iomsn Exp $")
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
    real, dimension (9) :: rc, theta_c
    integer :: j, l
!
!   clear the velocity field to zero
    f(:,:,:,iux:iuz) = 0.
!
    if (n_blobs > 0) then
      do j=1,n_blobs
        ! Shift the origin for xc, yc and zc.
        xc(j) = xc(j) + xyz0(1) + (xyz1(1)-xyz0(1))/2
        ! Transform the center of the rotation.
        rc(j) = sqrt(xc(j)**2 + yc(j)**2)
        theta_c(j) = atan2(yc(j),xc(j))
        do l=1,mx
          do m=1,my
            do n=1,mz
              f(l,m,n,iuz) = f(l,m,n,iuz) + ampl*kappa(j) * exp(-1/a_blob(j)**2 * (x(l)**2 - &
                    2*rc(j)*x(l)*(cos(theta_c(j))*cos(y(m))+sin(theta_c(j))*sin(y(m))) + &
                    rc(j)**2) - ((z(n)-zc(j))/l_blob(j))**2)
            enddo
          enddo
        enddo
      enddo
    endif
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
