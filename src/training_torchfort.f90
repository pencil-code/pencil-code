! $Id$
!
! CPARAM logical, parameter :: ltraining = .true.
!
!***************************************************************
!
  module Training

    use Cparam, only: fnlen
    use Messages
    use Syscalls
    use Torchfort

    implicit none

    character(LEN=fnlen) :: model='mymodel', config="config_mlp_native.yaml"
    integer(KIND=ikind8) :: dummy

    real, allocatable, device :: input(:,:,:,:), label(:,:,:,:), output(:,:,:,:)
    
    namelist /training_run_pars/ config
!
    contains
!***************************************************************
    subroutine initialize_training
 
      use General, only: itoa

      integer :: istat

      istat = torchfort_create_model(model, config, 0)
      if (istat /= TORCHFORT_RESULT_SUCCESS) then
        call fatal_error("training_initialize","istat="//trim(itoa(istat)))
      else
        call information('training_initialize','TORCHFORT LOADED SUCCESFULLY')
      endif

      allocate(input(nx, ny, 1, 1))
      allocate(output(nx, ny, 1, 1))
      allocate(label(nx, ny, 1, 1))

    endsubroutine initialize_training
!***********************************************************************
    subroutine read_training_run_pars(iostat)
!
! 23-jan-24/MR: coded
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=training_run_pars, IOSTAT=iostat)

    endsubroutine read_training_run_pars
!***************************************************************
    subroutine write_training_run_pars(unit)
!
      integer, intent(in) :: unit

      write(unit, NML=training_run_pars)

    endsubroutine write_training_run_pars
!***************************************************************
    subroutine training_train(f)
     
      real, dimension (mx,my,mz,mfarray) :: f

      integer :: istat
      real :: loss_val

      ! Device to host
      input(:,:,1,1) = f(l1:l2,m1:m2,n1,iux)

      istat = torchfort_train(model, input, label, loss_val)
      if (istat /= TORCHFORT_RESULT_SUCCESS) then
        call fatal_error("training_train","istat="//trim(itoa(istat)))
      else
        if (lroot) print*, "training loss:", loss_val
      endif

    endsubroutine training_train
!***************************************************************
    subroutine finalize_training

    endsubroutine finalize_training
!***************************************************************
  end module Training
