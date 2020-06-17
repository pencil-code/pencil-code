! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
!
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
  real :: eta0=0.0, k_eta, x0_drop, y0_drop
  real :: Omega_SB=1.0,gamma_parameter=1.0
!
  character (len=labellen), dimension(ninit) :: init_shallow_density='nothing'
  character (len=labellen), dimension(ninit) :: init_shallow_hydro='nothing'
!
  namelist /initial_condition_pars/  eta0, k_eta, x0_drop, y0_drop, Omega_SB, &
       init_shallow_density,init_shallow_hydro,gamma_parameter
!
  real :: planetary_radius=impossible ! derives from gamma_parameter and Omega
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
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
!  gamma_parameter = Omega/a**2 where a is the planetary radius
!
      planetary_radius = sqrt(Omega_SB/gamma_parameter)
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initial condition given by 
!
!     h = eta + Lb
!     rho is g*eta
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: eta,r2
      integer :: j
!
      do j=1,ninit
!
        select case (init_shallow_density(j))
!
        case('solid-body')   
!
          do n=n1,n2
            do m=m1,m2
!
! Equilibrium -> gh = gh0 + 3/2*Omega**2 * r**2 - Omega**2/(4*a**2) r**4
!
               r2 = x(l1:l2)**2 + y(m)**2
               eta = eta0 + Omega_SB**2*r2 * (1.5 - 0.25*gamma_parameter/Omega_SB * r2)
               f(l1:l2,m,n,ilnrho) = log(eta)
            enddo
          enddo
!
        case('gaussian-blob')   
!
!  with eta=a*exp(-k*(x-x0)**2) 
!
          do n=n1,n2
            do m=m1,m2
              eta = eta0 * exp(-k_eta * ( &  
                   (x(l1:l2)-x0_drop)**2 + (y(m)-y0_drop)**2 & 
                   ))
              f(l1:l2,m,n,ilnrho) = log(eta)
            enddo
          enddo
        endselect
!
      enddo
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initial condition given by 
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer :: j
!
      do j=1,ninit

        select case (init_shallow_hydro(j))
!
        case('solid-body')
          do n=n1,n2
            do m=m1,m2
              f(l1:l2,m,n,iux) = -Omega_SB * y(  m  )
              f(l1:l2,m,n,iuy) =  Omega_SB * x(l1:l2)
            enddo
          enddo
        endselect
      enddo
!
    endsubroutine initial_condition_uu
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
