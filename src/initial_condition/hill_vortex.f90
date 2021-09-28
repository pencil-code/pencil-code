!  $Id$
!
!  Initial condition (density, magnetic field, velocity) 
!  for magnetohydrostatical equilibrium in a global accretion
!  disk with an imposed (cylindrically symmetric) sound speed 
!  profile in spherical coordinates. 
!
!  10-feb-11/axel: adapted from noinitial_condition.f90
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
  real :: a_hill, u_hill
!
  namelist /initial_condition_pars/ &
     a_hill, u_hill
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
!  10-sep-15/axel: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx) :: r, r2, r3, r5, pom2
      real :: a2_hill, a3_hill
      integer :: iy,iz
!
!  initialize Hill vortex at x=y=z=0
!
      a2_hill=a_hill**2
      a3_hill=a_hill*a2_hill
      do iy=m1,m2; do iz=n1,n2
        pom2=x**2+y(iy)**2
        r2=x**2+y(iy)**2+z(iz)**2
        r=sqrt(r2)
        r3=r2*r
        r5=r2*r3
!
        where (r<=a_hill)
          f(:,iy,iz,iux)=u_hill*(-1.5*x    *z(iz)/a2_hill)
          f(:,iy,iz,iuy)=u_hill*(-1.5*y(iy)*z(iz)/a2_hill)
          f(:,iy,iz,iuz)=u_hill*(-2.5+1.5*(pom2+r2)/a2_hill)
        elsewhere
          f(:,iy,iz,iux)=u_hill*(-1.5*x    *z(iz)*a3_hill/r5)
          f(:,iy,iz,iuy)=u_hill*(-1.5*y(iy)*z(iz)*a3_hill/r5)
          f(:,iy,iz,iuz)=u_hill*(-a3_hill/r3+1.5*pom2*a3_hill/r5)
        endwhere
        f(:,iy,iz,ilnrho)=f(:,iy,iz,ilnrho)-.5*( &
          f(:,iy,iz,iux)**2+f(:,iy,iz,iuy)**2+f(:,iy,iz,iuz)**2)
      enddo; enddo
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
      use EquationOfState, only: cs20
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: argum,term1,term2,press,del_lnrho
!
      if (lroot) print*,&
           'initial_condition_lnrho: hill_vortex'
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential. Constant plasma 
!  beta magnetic field. 
!
!  07-may-09/wlad: coded

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: argum,term1,term2,ax,az,ay,fp,fm,f0
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
