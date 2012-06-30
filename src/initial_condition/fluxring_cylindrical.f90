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
  use Sub, only: erfunc
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: b0,s0,width,p0,eps=1.,mphi=1.,ampl=0.,om=1.,b1=0.,b2=0.,bz=0.,hel=1.,nohel=0.
  real :: omega_exponent=0.,ampl_diffrot=0.
  logical :: linitial_diffrot=.false.
!
  namelist /initial_condition_pars/ &
      b0,s0,width,p0,eps,mphi,ampl,om,b1,b2,linitial_diffrot,omega_exponent,&
     ampl_diffrot,hel,nohel
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
      real,dimension(mx) :: omega_diffrot
!
!
      if(linitial_diffrot) then
      do n=1,mz
        do m=1,my
          omega_diffrot = ampl_diffrot*x**(omega_exponent)
          f(:,m,n,iuy)=f(:,m,n,iuy)+ x*omega_diffrot
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
      del_lnrho=eps*alog(press/cs20)
!
      do n=1,mz
        do m=1,my
          f(:,m,n,ilnrho)=f(:,m,n,ilnrho)+del_lnrho
        enddo
      enddo
!
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential. Constant plasma 
!  beta magnetic field. 
!
!  07-may-09/wlad: coded
!  23-feb-12/fabio: added costant bz to the initial setup

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: argum,term1,term2,ax,az,ay
!
!  vector potential for the magnetic flux ring
!
      argum=(x-s0)/width
      term1=s0*sqrtpi*erfunc(argum)
      term2=-width*exp(-argum**2)
      az=-(.5*b0/s0)*width*(term1+term2)-b1*x-b2*log(x)
      ay=.5*bz*x
!
      do n=1,mz
        do m=1,my
!          f(:,m,n,iaz)=f(:,m,n,iaz)+az
          f(:,m,n,iaz)= az
          f(:,m,n,iay)=f(:,m,n,iay)+ay
        enddo
      enddo
!
!
!   perturbation for the initial field
!
print*,'ampl=',ampl

     do n=1,mz
      do m=1,my
        ax=ampl*x*(hel*cos(om*z(n))*sin(mphi*y(m))+nohel*sin(om*z(n))*sin(mphi*y(m)))
        az=ampl*x*cos(om*z(n))*cos(mphi*y(m))
        f(:,m,n,iax)=f(:,m,n,iax)+ax
        f(:,m,n,iaz)=f(:,m,n,iaz)+az
       enddo
     enddo


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
