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
    subroutine read_python_run_pars(iomsg)
!
      use General, only: keep_compiler_quiet
      character(LEN=*), intent(out) :: iomsg
      
      call keep_compiler_quiet(iomsg)

    endsubroutine read_python_run_pars
!***************************************************************
    subroutine write_python_run_pars(unit)
!
      use General, only: keep_compiler_quiet
      integer, intent(in) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine write_python_run_pars
!***************************************************************
    subroutine python_call
    endsubroutine python_call
!***************************************************************
    subroutine python_finalize
    endsubroutine python_finalize
!***************************************************************
  end module Python
