! $Id$
!
! CPARAM logical, parameter :: ltraining = .true.
!
! MAUX CONTRIBUTION 6
!
!***************************************************************
!
  module Training

    use Cdata
    use General, only: itoa
    use Messages
    !use Torchfort

    implicit none

    integer :: model_device=0
    integer :: it_train=-1, it_train_chkpt=-1

    real, dimension(:,:,:,:,:), allocatable, device :: input, label, output   !, device :: input, label, output

    integer :: itau, itauxx, itauxy, itauxz, itauyy, itauyz, itauzz

    character(LEN=fnlen) :: model='mymodel', config_file="config_mlp_native.yaml"
    character(LEN=fnlen) :: model_output_dir='data/ml_models/', checkpoint_output_dir='data/ml_models/', &
                            model_file = "model.pt", chkpt_file=""

    logical :: luse_trained_tau

    integer :: idiag_tauerror=0        ! DIAG_DOC: $\sqrt{\left<(\sum_{i,j} u_i*u_j - tau_{ij})^2\right>}$

    namelist /training_run_pars/ config_file, model_file, it_train, it_train_chkpt, luse_trained_tau
!
    integer :: istat, train_step_ckpt, val_step_ckpt   !, TORCHFORT_RESULT_SUCCESS=0
    logical :: ltrained=.false.

    contains
!***************************************************************
    subroutine initialize_training

      use File_IO, only: file_exists
      use Mpicomm, only: mpibcast
      use Syscalls, only: system_cmd

      character(LEN=fnlen) :: modelfn

      modelfn=trim(model_output_dir)//trim(model_file)
      if (lroot) then
        if (.not.file_exists(model_output_dir)) then
          call system_cmd('mkdir '//trim(model_output_dir))
        else
          ltrained = file_exists(trim(modelfn))
        endif
      endif
      call mpibcast(ltrained)

      if (lmpicomm) then
        istat = torchfort_create_distributed_model(model, config_file, MPI_COMM_WORLD, iproc)  !multinode?
      else
        istat = torchfort_create_model(model, config_file, model_device)
      endif
      if (istat /= TORCHFORT_RESULT_SUCCESS) then
        call fatal_error("initialize_training","when creating model: istat="//trim(itoa(istat)))
      else
        call information('initialize_training','TORCHFORT LOADED SUCCESFULLY')
      endif

      if (ltrained) then
        istat = torchfort_load_model(model, modelfn)
        if (istat /= TORCHFORT_RESULT_SUCCESS) &
          call fatal_error("initialize_training","when loading model: istat="//trim(itoa(istat)))
      else
        if (file_exists(trim(model_output_dir)//trim(chkpt_file))) then

          istat = torchfort_load_checkpoint("mymodel", checkpoint_output_dir, train_step_ckpt, val_step_ckpt)
          if (istat /= TORCHFORT_RESULT_SUCCESS) &
            call fatal_error("initialize_training","when loading checkpoint: istat="//trim(itoa(istat)))

        endif
      endif

      luse_trained_tau = luse_trained_tau.and.ltrained

      allocate(input(nx, ny, nz, 3, 1))
      allocate(output(nx, ny, nz, 6, 1))
      allocate(label(nx, ny, nz, 3, 1))

    endsubroutine initialize_training
!***********************************************************************
    subroutine register_training
!
!  Register slots in f-array for the six independent components of the Reynolds stress tensor tau.
!
      use FArrayManager
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
           "$Id$")
!
      call farray_register_pde('tau',itau,vector=6)
!
!  Indices to access tau.
!
      itauxx=itau; itauxy=itau+1; itauxz=itau+2; itauyy=itau+3; itauyz=itau+4; itauzz=itau+5

    endsubroutine register_training
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

      if (ltrained) then
        call infer(f)
      else
        call train(f)
      endif

    endsubroutine training_before_boundary
!***************************************************************
    subroutine infer(f)
     
      real, dimension (mx,my,mz,mfarray) :: f

      ! Host to device
      input(:,:,:,:,1) = f(l1:l2,m1:m2,n1:n2,iux:iuz)

      istat = torchfort_inference(model, input(:,:,:,:,1), output(:,:,:,:,1))
      if (istat /= TORCHFORT_RESULT_SUCCESS) then
        call fatal_error("infer","istat="//trim(itoa(istat)))
      else
        ! Device to host
        f(l1:l2,m1:m2,n1:n2,itauxx:itauzz) = output(:,:,:,:,1)
      endif

    endsubroutine infer
!***************************************************************
    subroutine train(f)
     
      real, dimension (mx,my,mz,mfarray) :: f

      real :: loss_val

      if (mod(it,it_train)==0) then
        ! Host to device
        input(:,:,:,:,1) = f(l1:l2,m1:m2,n1:n2,iux:iuz)

        istat = torchfort_train(model, input, label, loss_val)
        if (istat /= TORCHFORT_RESULT_SUCCESS) then
          call fatal_error("train","istat="//trim(itoa(istat)))
        else
          if (lroot) print*, "training loss = ", loss_val
        endif

        if (mod(it,it_train_chkpt)==0) then
          istat = torchfort_save_checkpoint(model, checkpoint_output_dir)
          if (istat /= TORCHFORT_RESULT_SUCCESS) &
            call fatal_error("train","when saving checkpoint: istat="//trim(itoa(istat)))
        endif
      endif

    endsubroutine train
!***************************************************************
    subroutine div_reynolds_stress(f,df)

      use Sub, only: div_other

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df

      real, dimension(nx,3) :: divrey
      real, dimension(mx,my,mz,3) :: tmp

      if (luse_trained_tau) then 
        call div_other(f(:,:,:,itauxx:itauxz),divrey(nx,1))

        tmp=f(:,:,:,(/itauxy,itauyy,itauyz/))
        call div_other(tmp,divrey(nx,2))

        tmp=f(:,:,:,(/itauxz,itauyz,itauzz/))
        call div_other(tmp,divrey(nx,3))

        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - divrey
      endif

    endsubroutine div_reynolds_stress
!***************************************************************
    subroutine calc_diagnostics_training(f,p)

      use Diagnostics, only: sum_mn_name

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      integer :: i,j,jtau
      real, dimension(nx) :: error

      if (ldiagnos.and.ltrained) then
        if (idiag_tauerror>0) then

          jtau=0
          error=0.
          do i=1,3
            do j=i,3
              error=error+(p%uu(:,i)*p%uu(:,j)-f(l1:l2,m,n,itau+jtau))**2
              jtau=jtau+1
            enddo
          enddo
          call sum_mn_name(error,idiag_tauerror,lsqrt=.true.)

        endif 
      endif 

    endsubroutine calc_diagnostics_training
!***********************************************************************
    subroutine rprint_training(lreset)
!
!  reads and registers print parameters relevant for training
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset
!
      if (lreset) then
        idiag_tauerror=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_training: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'tauerror',idiag_tauerror)
      enddo

    endsubroutine rprint_training
!***************************************************************
    subroutine finalize_training

!  Saving trained model.

      if (ltrained) then
        istat = torchfort_save_model(model, trim(model_output_dir)//trim(model_file))
        if (istat /= TORCHFORT_RESULT_SUCCESS) &
          call fatal_error("finalize_training","when saving model: istat="//trim(itoa(istat)))
      endif
      deallocate(input,label,output)

    endsubroutine finalize_training
!***************************************************************
  end module Training
