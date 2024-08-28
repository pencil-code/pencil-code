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

    character(LEN=fnlen) :: model='mymodel', config_file="config_mlp_native.yaml"
    
    integer :: model_device=0
    integer :: it_train

    real, allocatable, device :: input(:,:,:,:), label(:,:,:,:), output(:,:,:,:)
    
    namelist /training_run_pars/ config_file, it_train
!
    contains
!***************************************************************
    subroutine initialize_training
 
      use General, only: itoa

      integer :: istat

      istat = torchfort_create_model(model, config_file, model_device)
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
    subroutine training_before_boundary(f)
     
      real, dimension (mx,my,mz,mfarray) :: f

      call train(f)

    endsubroutine training_before_boundary
!***************************************************************
    subroutine train(f)
     
      real, dimension (mx,my,mz,mfarray) :: f

      integer :: istat
      real :: loss_val

      if (mod(it,it_train)==0) then
        ! Device to host
        input(:,:,1,1) = f(l1:l2,m1:m2,n1,iux)

        istat = torchfort_train(model, input, label, loss_val)
        if (istat /= TORCHFORT_RESULT_SUCCESS) then
          call fatal_error("training_train","istat="//trim(itoa(istat)))
        else
          if (lroot) print*, "training loss = ", loss_val
        endif
      endif

    endsubroutine train
!***************************************************************
    subroutine finalize_training

    endsubroutine finalize_training
!***************************************************************
  end module Training
