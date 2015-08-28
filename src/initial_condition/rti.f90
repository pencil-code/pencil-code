! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
! developed by Mikhail Modestov from a simpler version of Michiel Lambrechts
module InitialCondition
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: amplnrho=0.0, rho_ratio=1.0, widthrho=1.0
  real :: kx_ro=1.0, ro_pert=0.0
  real :: ampluu=0.0, widthuu=1.0, kx_uy=0.0  
  character (len=labellen) :: initlnrho='nothing'
  character (len=labellen) :: inituu='nothing'
  logical :: lzerow=.false.
!
  namelist /initial_condition_pars/ initlnrho, &
      amplnrho,rho_ratio,widthrho,ro_pert, kx_ro, &
      inituu,ampluu,widthuu,kx_uy,lzerow
!
  contains
!***********************************************************************
    subroutine initial_condition_all(f,profiles)
!
!  Initializes all the f arrays in one call. This subroutine is called last.
!  It's tested in 2D in x-y plane with gravy
!
!  10-feb-15/MR: added optional parameter 'profiles' (intended to replace f)
!
      use EquationOfState, only: gamma1,cs20,rho0
      use Gravity, only: gravx, gravy, gravz
!      
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (:,:),              optional, intent(out)  :: profiles
!
      real :: pp0, rhoprof, pprof
      real :: w1ro, ro_a, ro_b, kxro
      real :: kxu, w1u
      real, dimension (mx) :: x_ro
      real, dimension (my) :: y_ro
      real, dimension (mz) :: z_ro
      integer :: l,n,m  !loop indices for x and z direction
! initial distribution for density
      if (lroot) print*, &
          'initial_condition_all: tanh density profile for RT'
      w1ro = 1.0/widthrho
      kxro = 2.0*pi*kx_ro/Lxyz(1)
      select case (initlnrho)
      case('tanhy')      ! usual setting of lnrho
        ro_a=0.0;  ro_b=amplnrho
      case('tanhy_cos')  ! usual setting + cos perturbations
        ro_a=0.0;  ro_b=amplnrho
        forall (l=l1:l2) x_ro(l)=ro_pert*cos(kxro*x(l))
      case('tanhy_sin')  ! usual setting + sin perturbations
        ro_a=0.0;  ro_b=amplnrho
        forall (l=l1:l2) x_ro(l)=ro_pert*sin(kxro*x(l))
      case('tanhy_1')    ! density varies from 1/rho_ratio to 1
        ro_a=0.5*(1.0+1.0/rho_ratio)  
        ro_b=0.5*(1.0-1.0/rho_ratio)
      case('tanhy_1_cos')  ! + cos perturbations
        ro_a=0.5*(1.0+1.0/rho_ratio)  
        ro_b=0.5*(1.0-1.0/rho_ratio)
        forall (l=l1:l2) x_ro(l)=ro_pert*cos(kxro*x(l))
      case('tanhy_1_sin')  ! + sin perturbations
        ro_a=0.5*(1.0+1.0/rho_ratio)  
        ro_b=0.5*(1.0-1.0/rho_ratio)
        forall (l=l1:l2) x_ro(l)=ro_pert*sin(kxro*x(l))
      case default  ! usual setting of lnrho
        ro_a=0.0;  ro_b=amplnrho
      end select
!  compute pressure and entropy distribution
      pp0 = gamma1*rho0*cs20
      do l=l1,l2;  do m=m1,m2;  do n=n1,n2
        rhoprof = ro_a + ro_b*tanh( ( y(m) + x_ro(l) )*w1ro )
! pressure is computing without density perturbations if any -> baroclinic <> 0
        if (abs(y(m)*w1ro).lt.50) then
          pprof = pp0+gravy*( ro_a*y(m) + ro_b*widthrho*log(cosh(y(m)*w1ro)) )
         else
          pprof = pp0+gravy*( ro_a*y(m) + ro_b*abs(y(m)) )
        endif
        f(l,m,n,ilnrho) = log(rhoprof)
        f(l,m,n,iss) = gamma1*log(pprof/pp0) - log(rhoprof)
      enddo;  enddo;  enddo
!  initial conditions/perturbation for velocity field      
      kxu = 2.0*pi*kx_uy/Lxyz(1)
      w1u = 1.0/widthuu
      select case (inituu)
      case ('exp')
        do l=l1,l2;  do m=m1,m2
          f(l,m,n1:n2,iuy) = ampluu*cos(kxu*x(l))*exp(-abs(y(m))*w1u)
        enddo;  enddo
      case ('sech')
        do l=l1,l2;  do m=m1,m2
          f(l,m,n1:n2,iuy) = ampluu*cos(kxu*x(l))/cosh((y(m))*w1u)
        enddo;  enddo
      case ('noy')
        do l=l1,l2;  do m=m1,m2
          f(l,m,n1:n2,iuy) = ampluu*cos(kxu*x(l))
          if (lzerow) f(l,m,n1:n2,iux) =-ampluu*kxu*y(m)*sin(kxu*x(l))
        enddo;  enddo
      case ('zero')
        f(l1:l2,m1:m2,n1:n2,iuy)=0.0
        f(l1:l2,m1:m2,n1:n2,iux)=0.0
        f(l1:l2,m1:m2,n1:n2,iuz)=0.0
      case default
        f(l1:l2,m1:m2,n1:n2,iuy)=0.0
        f(l1:l2,m1:m2,n1:n2,iux)=0.0
        f(l1:l2,m1:m2,n1:n2,iuz)=0.0
      end select
      if (.not.lzerow) f(l1:l2,m1:m2,n1:n2,iux)=0.0
!
     if (present(profiles)) then
       call fatal_error('initial_condition_all', &
                        'returning of profiles not implemented')
       call keep_compiler_quiet(profiles)
     endif
!
    endsubroutine initial_condition_all
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
