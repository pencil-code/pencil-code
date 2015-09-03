! $Id$
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
  use General, only: keep_compiler_quiet
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: rhopload=1.
!
  namelist /initial_condition_pars/ &
    rhopload
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
          "$Id: noinitial_condition.f90 16806 2011-05-04 15:52:49Z dhruba.mitra $")
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
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f,profiles)
!
!  Initializes all the f arrays in one call. This subroutine is called last.
!
!  10-feb-15/MR: added optional parameter 'profiles' (intended to replace f)
!
      use EquationOfState, only: cs20, rho0 
      use Gravity, only: gravz
!      
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (:,:),              optional, intent(out)  :: profiles
!      
      real :: rhoprof
      integer :: l,n  !loop indices for x and z direction
!
!  Isothermal hydrotstatic equilibrium startification for rho0, 
!  in case of loading the box witth particles with therminal velocity
!  TODO replace rhopload (particle density) 
!  with variable from the particle module
!  TODO rho0 could be made independent variable
!  
          if (lroot) print*, &
              'initial_condition_all: isotherm hydrostat equi with part load'
          do n=n1,n2
              rhoprof= &
                  (rho0 + rhopload) * exp((gravz/cs20) * z(n)) - rhopload 
              if (ldensity_nolog) then
                f(l1:l2,m1:m2,n,irho) = rhoprof
              else
                f(l1:l2,m1:m2,n,ilnrho)= log(rhoprof) 
              endif
          enddo
!
      call keep_compiler_quiet(f)
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
