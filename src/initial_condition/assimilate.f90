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
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: dummy
!
  namelist /initial_condition_pars/ &
      dummy
!
contains
!***********************************************************************
  subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  04-sep-10/bing: coded
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
!  14-dec-10/bing: coded
!
    real, dimension (mx,my,mz,mfarray) :: f
!
    call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
  subroutine read_initial_condition_pars(unit,iostat)
!
!  04-sep-10/bing: coded
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
 99  return
!
  endsubroutine read_initial_condition_pars
!***********************************************************************
  subroutine write_initial_condition_pars(unit)
!
!  04-sep-10/bing: coded
!
    integer, intent(in) :: unit
!
    write(unit,NML=initial_condition_pars)
!
  endsubroutine write_initial_condition_pars
!***********************************************************************
  subroutine initial_condition_all(f,df)
!
!  Initialize logarithmic density.
!
!  04-sep-10/bing: coded
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
    real, dimension (mx,my,mz,mfarray), intent(in) :: df
!
print*,'f'
print*,f(l1:l2,m1:m2,n1:n2,:)
print*
print*,'df'
print*,df(l1:l2,m1:m2,n1:n2,:)
print*
    f=f+.0001*df
print*,'f'
print*,f(l1:l2,m1:m2,n1:n2,:)
print*
!
  endsubroutine initial_condition_all
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
