! $Id$
!
! CPARAM logical, parameter :: ltraining = .false.
!
!***************************************************************
!
  module Training

    use Cparam
    use General, only: keep_compiler_quiet

    implicit none

    contains
!***************************************************************
    subroutine initialize_training
 
    endsubroutine initialize_training
!***********************************************************************
    subroutine read_training_run_pars(iostat)
!
      integer, intent(out) :: iostat

      call keep_compiler_quiet(iostat)
!
    endsubroutine read_training_run_pars
!***************************************************************
    subroutine write_training_run_pars(unit)
!
      integer, intent(in) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine write_training_run_pars
!***************************************************************
    subroutine training_before_boundary(f)

      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)

    endsubroutine training_before_boundary
!***************************************************************
    subroutine finalize_training

    endsubroutine finalize_training
!***************************************************************
  end module Training
