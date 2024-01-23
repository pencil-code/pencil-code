! $Id$
!
! CPARAM logical, parameter :: lpython = .true.
!
!***************************************************************
!
  module Python

    use Cparam, only: fnlen
    use Syscalls

    implicit none

    character(LEN=fnlen) :: pymodule='', pyfunction=''
    integer(KIND=ikind8) :: pModule, pFunction

    namelist /python_run_pars/ pymodule, pyfunction
!
    contains
!***************************************************************
    subroutine python_init

      call py_init

    endsubroutine python_init
!***************************************************************
    subroutine python_initialize

      call py_initialize(trim(pymodule)//char(0),trim(pyfunction)//char(0),pModule,pFunction)

    endsubroutine python_initialize
!***************************************************************
    subroutine python_call
     
      integer :: arg
      call py_call(pFunction,arg)

    endsubroutine python_call
!***************************************************************
    subroutine python_finalize

      call py_finalize(pModule,pFunction)

    endsubroutine python_finalize
!***************************************************************
  end module Python
