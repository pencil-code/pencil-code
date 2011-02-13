!  $Id$
!
!  Initial condition (density, magnetic field, velocity) 
!  for magnetohydrostatical equilibrium in a global accretion
!  disk with an imposed (cylindrically symmetric) sound speed 
!  profile in spherical coordinates. 
!
!  07-may-09/wlad: adapted from noinitial_condition.f90
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
  use Cdata
  use Cparam
  use Messages
  use Sub, only: erfunc,keep_compiler_quiet
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: b0,s0,width,p0
!
  namelist /initial_condition_pars/ &
      b0,s0,width,p0
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
           'initial_condition_lnrho: ring'
!
!  density for the magnetic flux flux ring
!
      argum=sqrt2*(x-s0)/width
      term1=s0*width*sqrtpi*sqrt2*erfunc(argum)
      term2=(2.*x**2-width**2)*exp(-argum**2)
      press=p0-(.5*b0/s0)**2*(term1+term2)
      del_lnrho=alog(press/cs20)
!
      do n=1,mz
        do m=1,my
          f(:,m,n,ilnrho)=f(:,m,n,ilnrho)+del_lnrho
        enddo
      enddo
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential. Constant plasma 
!  beta magnetic field. 
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: argum,term1,term2,az
!
!  vector potential for the magnetic flux flux ring
!
      argum=(x-s0)/width
      term1=s0*sqrtpi*erfunc(argum)
      term2=-width*exp(-argum**2)
      az=-(.5*b0/s0)*width*(term1+term2)
!
      do n=1,mz
        do m=1,my
          f(:,m,n,iaz)=f(:,m,n,iaz)+az
        enddo
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
