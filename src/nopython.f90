! $Id$
!
! CPARAM logical, parameter :: lpython = .false.
!
!***************************************************************
!
  module Python

    implicit none
!
    contains
!***************************************************************
    subroutine python_init
    endsubroutine python_init
!***************************************************************
    subroutine python_initialize
    endsubroutine python_initialize
!***************************************************************
    subroutine read_python_run_pars(iostat)
!
      use General, only: keep_compiler_quiet
      integer, intent(out) :: iostat
      
      call keep_compiler_quiet(iostat)

    endsubroutine read_python_run_pars
!***************************************************************
    subroutine write_python_run_pars(unit)
!
      integer, intent(in) :: unit

    endsubroutine write_python_run_pars
!***************************************************************
    subroutine python_call
    endsubroutine python_call
!***************************************************************
    subroutine python_finalize
    endsubroutine python_finalize
!***************************************************************
  end module Python
