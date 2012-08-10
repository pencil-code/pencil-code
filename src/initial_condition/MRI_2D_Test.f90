!  $Id: fluxring_cylindrical.f90 19341 2012-08-01 12:11:20Z AxelBrandenburg $
!
!  Initial condition (density, magnetic field, velocity) 
!  for magnetohydrostatical equilibrium in a global accretion
!  disk with an imposed (cylindrically symmetric) sound speed 
!  profile in spherical coordinates. 
!
!  10-aug-12/axel: adapted from fluxring_cylindrical.f90
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
  use Sub, only: erfunc
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: pi2=2.*pi, pi6=6.*pi, B0=.2
  integer :: l
!
  namelist /initial_condition_pars/ B0
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
           "$Id: fluxring_cylindrical.f90 19341 2012-08-01 12:11:20Z AxelBrandenburg $")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  10-aug-12/axel: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  velocity for 2D MRI test
!
      do n=1,mz
      do l=1,mx
        f(l,:,n,iux)=cos(pi2*z(m))
        f(l,:,n,iuy)=cos(pi2*(x(l)+3*z(n)))
        f(l,:,n,iuz)=cos(pi2*x(l))
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
!  10-aug-12/axel: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: argum,term1,term2,press,del_lnrho
!
!  constant density, following convention of
!  http://www.astro.princeton.edu/~jstone/Athena/tests/orszag-tang/pagesource.html
!
      if (lroot) print*,'initial_condition_lnrho: 2D MRI test'
      f(:,:,:,ilnrho)=0.
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential. Constant plasma 
!  beta magnetic field. 
!
!  10-aug-12/axel: coded

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  vector potential for 2D MRI test
!
      do n=1,mz
      do l=1,mx
        f(l,:,n,iay)=B0*(sin(pi6*x(l))/pi6-sin(pi2*z(n))/pi2)
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
