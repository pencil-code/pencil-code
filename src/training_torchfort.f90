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
    use Cudafor
    use Torchfort
    use iso_c_binding

    implicit none

    integer :: model_device=0
    integer :: it_train=-1, it_train_chkpt=-1

    real, dimension(:,:,:,:,:), allocatable, device :: input, label, output

    integer :: itau, itauxx, itauxy, itauxz, itauyy, itauyz, itauzz

    character(LEN=fnlen) :: model_output_dir='training/', checkpoint_output_dir='training'
    character(LEN=fnlen) :: model='model', config_file="training/config_mlp_native.yaml", model_file

    logical :: luse_trained_tau
    real :: max_loss=1.e-4

    integer :: idiag_loss=0            ! DIAG_DOC: torchfort training loss
    integer :: idiag_tauerror=0        ! DIAG_DOC: $\sqrt{\left<(\sum_{i,j} u_i*u_j - tau_{ij})^2\right>}$

    namelist /training_run_pars/ config_file, model, it_train, it_train_chkpt, luse_trained_tau, max_loss
!
    integer :: istat, train_step_ckpt, val_step_ckpt
    logical :: ltrained=.false.

    contains
!***************************************************************
    subroutine initialize_training

      use File_IO, only: file_exists
      use Mpicomm, only: mpibcast, MPI_COMM_WORLD
      use Syscalls, only: system_cmd

      character(LEN=fnlen) :: modelfn

      if (.not.lhydro) call fatal_error('initialize_training','needs HYDRO module')
      istat = cudaSetDevice(iproc)
      if (istat /= CUDASUCCESS) call fatal_error('initialize_training','cudaSetDevice failed')
   
      model_file = trim(model)//'.pt'
      modelfn=trim(model_output_dir)//trim(model_file)

      if (lroot) then
        if (.not.file_exists(model_output_dir)) then
          call system_cmd('mkdir '//trim(model_output_dir))
        else
          ltrained = file_exists(trim(modelfn))
        endif
      endif
      call mpibcast(ltrained)
!
! TorchFort create model
!
      if (lmpicomm) then
        istat = torchfort_create_distributed_model(model, config_file, MPI_COMM_WORLD, iproc)
      else
        istat = torchfort_create_model(model, config_file, model_device)
      endif
      if (istat /= TORCHFORT_RESULT_SUCCESS) then
        call fatal_error("initialize_training","when creating model: istat="//trim(itoa(istat)))
      else
        call information('initialize_training','TORCHFORT LIB LOADED SUCCESFULLY')
      endif

      if (ltrained) then
        istat = torchfort_load_model(model, modelfn)
        if (istat /= TORCHFORT_RESULT_SUCCESS) then
          call fatal_error("initialize_training","when loading model: istat="//trim(itoa(istat)))
        else
          call information('initialize_training','TORCHFORT MODEL "'//trim(modelfn)//'" LOADED SUCCESFULLY')
        endif
      else
        if (file_exists(trim(model_output_dir)//trim(model)//'.ckpt')) then

          istat = torchfort_load_checkpoint(model, checkpoint_output_dir, train_step_ckpt, val_step_ckpt)
          if (istat /= TORCHFORT_RESULT_SUCCESS) then
            call fatal_error("initialize_training","when loading checkpoint: istat="//trim(itoa(istat)))
          else
            call information('initialize_training','TORCHFORT CHECKPOINT LOADED SUCCESFULLY')
          endif

        endif
      endif

      luse_trained_tau = luse_trained_tau.and.ltrained

      if (.not.lgpu) then
        allocate(input(mx, my, mz, 3, 1))
        allocate(output(mx, my, mz, 6, 1))
        allocate(label(mx, my, mz, 6, 1))
      endif

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
      call farray_register_auxiliary('tau',itau,vector=6)
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
    
      use Gpu, only: get_ptr_gpu

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:,:,:,:), pointer :: ptr_uu, ptr_tau

      ! Host to device
      if (.not.lgpu) then
        input(:,:,:,:,1) = f(:,:,:,iux:iuz)    ! host to device
        istat = torchfort_inference(model, input, output)
      else
        !call get_ptr_gpu(ptr_uu,iux,iuz)
        !call get_ptr_gpu(ptr_tau,tauxx,tauzz)
        istat = torchfort_inference(model, get_ptr_gpu(iux,iuz), get_ptr_gpu(itauxx,itauzz))
      endif

      if (istat /= TORCHFORT_RESULT_SUCCESS) then
        call fatal_error("infer","istat="//trim(itoa(istat)))
      elseif (.not.lgpu) then
        ! Device to host
        f(l1:l2,m1:m2,n1:n2,itauxx:itauzz) = output(:,:,:,:,1)
      endif

    endsubroutine infer
!***************************************************************
    subroutine train(f)
   
      use Diagnostics, only: save_name 
      use Gpu, only: get_ptr_gpu
      use Sub, only: smooth

      real, dimension (mx,my,mz,mfarray) :: f

      real :: loss_val

      if (mod(it,it_train)==0) then
!
!  Smooth velocity.
!
        if (lgpu) then
          !smoothing kernel for uu and tau
          istat = torchfort_inference(model, get_ptr_gpu(iux,iuz), get_ptr_gpu(itauxx,itauzz))
        else
          f(:,:,:,itau:itau+2) = f(:,:,:,iux:iuz)      ! use f(:,:,:,itau:itau+2) as work space
          call smooth(f,itau,itau+2)
          input(:,:,:,:,1) = f(:,:,:,itau:itau+2)      ! host to device
!
!  Calculate and smooth stress tensor.
!
          f(:,:,:,itauxx) = f(:,:,:,iux)**2
          f(:,:,:,itauyy) = f(:,:,:,iuy)**2
          f(:,:,:,itauzz) = f(:,:,:,iuz)**2
          f(:,:,:,itauxy) = f(:,:,:,iux)*f(:,:,:,iuy)
          f(:,:,:,itauyz) = f(:,:,:,iuy)*f(:,:,:,iuz)
          f(:,:,:,itauxz) = f(:,:,:,iux)*f(:,:,:,iuz)

          call smooth(f,itauxx,itauzz)
          label(:,:,:,:,1) = f(:,:,:,itauxx:itauzz)    ! host to device

          istat = torchfort_train(model, input, label, loss_val)
        endif

        if (istat /= TORCHFORT_RESULT_SUCCESS) then
          call fatal_error("train","istat="//trim(itoa(istat)))
        else
          if (ldiagnos) call save_name(loss_val, idiag_loss)
        endif

        if (loss_val <= max_loss) ltrained=.true.
        if (lroot.and.mod(it,it_train_chkpt)==0) then
          istat = torchfort_save_checkpoint(model, checkpoint_output_dir)
          if (istat /= TORCHFORT_RESULT_SUCCESS) &
            call fatal_error("train","when saving checkpoint: istat="//trim(itoa(istat)))
        endif
      endif

    endsubroutine train
!***************************************************************
    subroutine div_reynolds_stress(f,df)

      use Sub, only: div

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df

      real, dimension(nx,3) :: divrey

      if (luse_trained_tau) then 
        call div(f,itauxx,divrey(:,1))
        call div(f,0,divrey(:,2),inds=(/itauxy,itauyy,itauyz/))
        call div(f,0,divrey(:,3),inds=(/itauxz,itauyz,itauzz/))

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

      if (ltrained) then
        if (ldiagnos) then
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
      endif 
!
    endsubroutine calc_diagnostics_training
!***********************************************************************
    subroutine get_slices_training(f,slices)
!
!  Write slices for animation of predicted Reynolds stresses.
!
      use Slices_methods, only: assign_slices_scal

      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Velocity field.
!
        case ('tauxx'); call assign_slices_scal(slices,f,itauxx)
        case ('tauxy'); call assign_slices_scal(slices,f,itauxy)
        case ('tauxz'); call assign_slices_scal(slices,f,itauxz)
        case ('tauyy'); call assign_slices_scal(slices,f,itauyy)
        case ('tauyz'); call assign_slices_scal(slices,f,itauyz)
        case ('tauzz'); call assign_slices_scal(slices,f,itauzz)

      end select

    endsubroutine get_slices_training
!***********************************************************************
    subroutine rprint_training(lreset)
!
!  reads and registers print parameters relevant for training
!
      use Diagnostics, only: parse_name
!
      integer :: iname, inamev, idum
      logical :: lreset
!
      if (lreset) then
        idiag_tauerror=0; idiag_loss=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_training: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'loss',idiag_loss)
        call parse_name(iname,cname(iname),cform(iname),'tauerror',idiag_tauerror)
      enddo

      do inamev=1,nnamev
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauxx',idum)
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauxy',idum) 
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauxz',idum) 
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauyy',idum) 
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauyz',idum) 
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauzz',idum) 
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
      if (.not.lgpu) deallocate(input,label,output)

    endsubroutine finalize_training
!***************************************************************
    subroutine write_sample(sample, fname)
      
!      use HDF5

      real, dimension(:,:,:), intent(in) :: sample
      character(len=*), intent(in) :: fname

!      integer(HID_T) :: in_file_id
!      integer(HID_T) :: out_file_id
!      integer(HID_T) :: dset_id
!      integer(HID_T) :: dspace_id
    
      integer :: err
    
      !integer(HSIZE_T) :: dims(size(shape(sample)))
      integer :: dims(size(shape(sample)))
    
!      call h5open_f(err)
!      call h5fcreate_f (fname, H5F_ACC_TRUNC_F, out_file_id, err)
    
      dims = shape(sample)
!      call h5screate_simple_f(size(shape(sample)), dims, dspace_id, err)
!      call h5dcreate_f(out_file_id, "data", H5T_NATIVE_REAL, dspace_id, dset_id, err)
!      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sample, dims, err)
!      call h5dclose_f(dset_id, err)
!      call h5sclose_f(dspace_id, err)
!    
!      call h5fclose_f(out_file_id, err)
!      call h5close_f(err)
    
    endsubroutine write_sample
!***************************************************************
  end module Training
