!  Initial condition (density, magnetic field, velocity) 
!  for a particular configuration of a twsited magnetic field tube(s).
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
  use Boundcond ! for the core boundary communication
!
  implicit none
!
!   include 'mpif.h'
  include '../initial_condition.h'
!
! B_bkg = strength of the background field in z-direction
! B_theta1,B_theta2 = strength of the azimuthal component of the field
! r1,r2 = radii of the twisted tubes
! x1, y1, x2, y2 = position of the tubes
! config = type of configuration (pairs or array)
! pert = magnitude of the random pertubation
!
  real :: B_bkg = 1.0, B_theta1 = 1., B_theta2 = 0.
  real :: r1 = 0., r2 = 0.
  real :: x1 = 0., y1 = 0., x2 = 1., y2 = 0.
  character (len=labellen) :: config = 'pair'
  real :: pert
!
  namelist /initial_condition_pars/ &
    B_bkg, B_theta1, B_theta2, r1, r2, x1, y1, x2, y2, config, pert
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
!  24-july-15/simon: coded
!
!  Twisted magnetic flux tube(s) starting from the lower xy-plane and
!  ending at the top plane.
!
!  Created 2015-07-24 by Simon Candelaresi (Iomsn)
!
    use Mpicomm, only: stop_it
    use Poisson
    use Sub
!    
    real, dimension (mx,my,mz,mfarray) :: f
    real :: r
    real, dimension (mz) :: p
!        
!   The next variables are used for the uncurling.
    integer :: l, j, ju, k
    real, dimension (nx,ny,nz,3) :: jj, tmpJ  ! This is phi for poisson.f90
!
!   clear the magnetic field to zero
    f(:,:,:,iax:iaz) = 0.
!
!   add the two first twisted tube
    do l=1,mx
      do m=1,my
        r = sqrt((x(l)-x1)**2+(y(m)-y1)**2)
        f(l,m,:,iax) = f(l,m,:,iax) - B_theta1*(-(tanh((r-r1)*10/r1)+1))*(y(m)-y1)/2/r1
        f(l,m,:,iay) = f(l,m,:,iay) + B_theta1*(-(tanh((r-r1)*10/r1)+1))*(x(l)-x1)/2/r1
        r = sqrt((x(l)-x2)**2+(y(m)-y2)**2)
        f(l,m,:,iax) = f(l,m,:,iax) - B_theta2*(-(tanh((r-r2)*10/r2)+1))*(y(m)-y2)/2/r2
        f(l,m,:,iay) = f(l,m,:,iay) + B_theta2*(-(tanh((r-r2)*10/r2)+1))*(x(l)-x2)/2/r2
!       add the pertubation
        call random_number(p)
        f(l,m,:,iax) = f(l,m,:,iax) + pert*p
        call random_number(p)
        f(l,m,:,iay) = f(l,m,:,iay) + pert*p
        call random_number(p)
        f(l,m,:,iaz) = f(l,m,:,iaz) + pert*p
      enddo
    enddo    
!
!   Transform the magnetic field into a vector potential
!
!   Compute curl(B) = J for the Poisson solver
    do m=m1,m2
      do n=n1,n2
        call curl(f,iaa,jj(:,m-nghost,n-nghost,:))
      enddo
    enddo
    tmpJ = -jj
!
!   Use the Poisson solver to solve \nabla^2 A = -J for A
    do j=1,3
      call inverse_laplacian(tmpJ(:,:,:,j))
    enddo
!
!   Overwrite the f-array with the correct vector potential A
    do j=1,3
      ju=iaa-1+j
      f(l1:l2,m1:m2,n1:n2,ju) = tmpJ(:,:,:,j)
    enddo
!
!   Add a background field to the braid
    do l=1,mx
      do m=1,my
        f(l,m,:,iax) = f(l,m,:,iax) - y(m)*B_bkg/2.
        f(l,m,:,iay) = f(l,m,:,iay) + x(l)*B_bkg/2.
      enddo
    enddo
!
  endsubroutine initial_condition_aa
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: get_unit
!
      integer, intent(out) :: iostat
      include "../parallel_unit.h"
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
