! $Id: waves_Bx.f90,v 1.0 2011-06-29 18:02:00 simon Exp $
!
! Initial condition for sin(kz) or cos(kz) magnetic field.
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
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../initial_condition.h'
!
  ! k = wave vector
  ! initaa = implitude of the field
  real :: k = 2.0, initaa = 1.0
  character (len=10) :: profile = 'cos'
  
  namelist /initial_condition_pars/ &
      k, initaa, profile
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  29-june-11/simon: coded
!
      if (lroot) call svn_id( &
           "$Id: bubbles_init.f90,v 1.0 2011-06-29 18:02:00 simon Exp $")
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
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
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
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  29-jun-11/simon: coded
!
!  Simple wave configuration in z-direction
!
!  Created 2011-06-29 by Simon Candelaresi (Iomsn)
!
      use Mpicomm, only: stop_it
      use Sub
      
      real, dimension (mx,my,mz,mfarray) :: f
!    
    f(:,:,:,iax) = 0.
    f(:,:,:,iaz) = 0.
    do n = n1, n2, 1        
        select case (profile)
          case ('sin')
            f(:,:,n,iay) = initaa * cos(k*z(n))/k
          case ('cos')
            f(:,:,n,iay) = -initaa * sin(k*z(n))/k
        endselect        
    enddo          
!
    endsubroutine initial_condition_aa
!***********************************************************************    
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!********************************************************************
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
