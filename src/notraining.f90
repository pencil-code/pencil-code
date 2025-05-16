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
    include 'training.h'

    contains
!***************************************************************
    subroutine initialize_training
 
    endsubroutine initialize_training
!***********************************************************************
    subroutine register_training

    endsubroutine register_training
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
    subroutine training_after_boundary(f)

      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)

    endsubroutine training_after_boundary
!***************************************************************
    subroutine calc_diagnostics_training(f)

      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)

    endsubroutine calc_diagnostics_training
!***********************************************************************
    subroutine rprint_training(lreset)
!
      logical :: lreset
!
      call keep_compiler_quiet(lreset)

    endsubroutine rprint_training
!***************************************************************
    subroutine div_reynolds_stress(f,df)

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)

    endsubroutine div_reynolds_stress
!***************************************************************
    subroutine finalize_training

    endsubroutine finalize_training
!***********************************************************************
    subroutine get_slices_training(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices)
      
    endsubroutine get_slices_training
!***************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    endsubroutine pushpars2c
!***********************************************************************
  end module Training
