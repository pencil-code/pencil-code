! $Id$
!
! CPARAM logical, parameter :: ltraining = .true.
!
! MAUX CONTRIBUTION 6
! COMMUNICATED AUXILIARIES 6
!
!***************************************************************
!
  module Training

    use Cdata
    use General, only: itoa
    use Messages
    use Cudafor, only: cudaSetDevice,CUDASUCCESS,cudaGetDeviceCount
    use Torchfort, only: torchfort_create_distributed_model, torchfort_create_model,&
                         torchfort_result_success,torchfort_load_model,torchfort_load_checkpoint,&
                         torchfort_save_model,torchfort_result_success,torchfort_save_checkpoint,&
                         torchfort_inference,torchfort_train
    !use iso_c_binding

    implicit none

    include 'training.h'

    integer :: model_device=0
    integer :: it_train=-1, it_train_chkpt=-1, it_train_start=1,it_train_end=-1

    !real(KIND=rkind4), dimension(:,:,:,:,:), allocatable, device :: input, label, output
    real, dimension(:,:,:,:,:), allocatable, device :: input, label, output
    real :: train_loss   !(KIND=rkind4) :: train_loss

    integer :: itau, itauxx, itauxy, itauxz, itauyy, itauyz, itauzz

    character(LEN=fnlen) :: model='model', config_file="config_mlp_native.yaml", model_file

    logical :: lroute_via_cpu=.false., lfortran_launched, luse_trained_tau, lwrite_sample=.false., lscale=.true.
    real :: max_loss=1.e-4, dt_train=1.e-10

    integer :: idiag_loss=0            ! DIAG_DOC: torchfort training loss
    integer :: idiag_tauerror=0        ! DIAG_DOC: $\sqrt{\left<(\sum_{i,j} u_i*u_j - tau_{ij})^2\right>}$

    namelist /training_run_pars/ config_file, model, it_train, it_train_start, it_train_chkpt, &
                                 luse_trained_tau, lscale, lwrite_sample, max_loss, lroute_via_cpu,&
                                 it_train_end, lrun_epoch, dt_train
!
    character(LEN=fnlen) :: model_output_dir, checkpoint_output_dir
    integer :: istat, train_step_ckpt, val_step_ckpt
    logical :: ltrained=.false., lckpt_written=.false.,lmodel_saved=.false., lrun_epoch=.false.
    real, dimension (mx,my,mz,3) :: uumean
    real :: tauerror, input_min, input_max, output_min, output_max
    real, dimension(mx, my, mz, 6) :: tau_pred
    real::  start_time, end_time

    contains
!***************************************************************
    subroutine initialize_training

      use File_IO, only: file_exists
      use Mpicomm, only: mpibcast, MPI_COMM_PENCIL
      use Syscalls, only: system_cmd

      character(LEN=fnlen) :: modelfn
      integer :: ndevs

      lfortran_launched = .not. lgpu .or. lroute_via_cpu
      if (lreloading) return

      if (.not.lhydro) call fatal_error('initialize_training','needs HYDRO module')

      istat = cudaGetDeviceCount(ndevs)
      if (istat /= CUDASUCCESS) call fatal_error('initialize_training','cudaGetDeviceCount failed')
      istat = cudaSetDevice(mod(iproc,ndevs))
      if (istat /= CUDASUCCESS) call fatal_error('initialize_training','cudaSetDevice failed')
  
      model_output_dir=trim(datadir)//'/training/' 
      checkpoint_output_dir=model_output_dir
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
      print*, 'CONFIG FILE=', trim(model_output_dir)//trim(config_file)
      if (lmpicomm) then
        istat = torchfort_create_distributed_model(trim(model), trim(model_output_dir)//trim(config_file), &
                                                   MPI_COMM_PENCIL, mod(iproc,ndevs))
      else
        istat = torchfort_create_model(trim(model), trim(model_output_dir)//trim(config_file), model_device)
      endif
      if (istat /= TORCHFORT_RESULT_SUCCESS) then
        call fatal_error("initialize_training","when creating model "//trim(model)//": istat="//trim(itoa(istat)))
      else
        call information('initialize_training','TORCHFORT LIB LOADED SUCCESFULLY')
      endif

!need this to be false for now but should be ltrained
      if (ltrained.and..not.lrun_epoch) then
        istat = torchfort_load_model(trim(model), trim(modelfn))
        if (istat /= TORCHFORT_RESULT_SUCCESS) then
          call fatal_error("initialize_training","when loading model: istat="//trim(itoa(istat)))
        else
          call information('initialize_training','TORCHFORT MODEL "'//trim(modelfn)//'" LOADED SUCCESFULLY')
        endif
      else
        if (file_exists(trim(checkpoint_output_dir)//'/'//trim(model)//'.pt').and.lroot) then
          print *, 'loaded checkpoint'
          ltrained=.false.
          istat = torchfort_load_checkpoint(trim(model), trim(checkpoint_output_dir), train_step_ckpt, val_step_ckpt)
          if (istat /= TORCHFORT_RESULT_SUCCESS) then
            call fatal_error("initialize_training","when loading checkpoint: istat="//trim(itoa(istat)))
          else
            call information('initialize_training','TORCHFORT CHECKPOINT LOADED SUCCESFULLY')
          endif

        endif
      endif

      ltrained = (ltrained .or. luse_trained_tau) .and. .not. lrun_epoch

      if (lrun .and. lfortran_launched) then
        allocate(input (mx, my, mz, 3, 1))
        allocate(output(mx, my, mz, 6, 1))
        allocate(label (mx, my, mz, 6, 1))
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
      call farray_register_auxiliary('tau',itau,vector=6,on_gpu=lgpu,communicated=.true.)
!
!  Indices to access tau.
!
      itauxx=itau; itauyy=itau+1; itauzz=itau+2; itauxy=itau+3; itauxz=itau+4; itauyz=itau+5

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
    subroutine training_after_boundary(f)
     
      use Sub, only: smooth

      real, dimension (mx,my,mz,mfarray) :: f

      if (ltrained) then

        call infer(f)

        ! added false since we dont need these
        !if ((ldiagnos.or.lvideo).and.lfirst) then
        if (.false.) then
          call calc_tau(f)
          if (lfortran_launched) then
!
! Copy data from device to host.
!
            f(:,:,:,itauxx:itauzz) = f(:,:,:,itauxx:itauzz) - output(:,:,:,:,1)
          endif
          tauerror = sum(f(l1:l2,m1:m2,n1:n2,itauxx:itauzz)**2)/nx
        else
          f(:,:,:,itauxx:itauzz) = output(:,:,:,:,1)
        endif
      else
        if (lfirst) call train(f)
      endif
!
! output for plotting
!
      ! added false since there is another way for writing samples
      !if (lvideo .or. lwrite_sample .and. mod(it, 50)==0) then
      if (.false.) then
!     
        call calc_tau(f)
         
        call infer(f)
        tau_pred = output(:,:,:,:,1)                 ! device to host
        if (lscale) call descale(tau_pred, output_min, output_max)

        if (lwrite_sample .and. mod(it, 50)==0) then
          call write_sample(f(:,:,:,itauxx), mx, my, mz, "target_"//trim(itoa(iproc))//".hdf5")
          call write_sample(tau_pred(:,:,:,1), mx, my, mz, "pred_"//trim(itoa(iproc))//".hdf5")
        endif

      endif

    endsubroutine training_after_boundary
!***************************************************************
    subroutine infer(f)
    
      use Gpu, only: get_ptr_gpu_training, infer_gpu
      use Sub, only: smooth
      use Mpicomm, only: mpiwtime

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:,:,:,:), pointer :: ptr_uu, ptr_tau

      if (lfortran_launched) then
!
! Smooth velocity -> "mean".
!
        uumean = f(:,:,:,iux:iuz)
        call smooth(uumean,1,3,lgauss=.true.)
        if (lscale) call scale(uumean, input_min, input_max)
!
! Copy data from host to device.
!
        input(:,:,:,:,1) = uumean
        istat = torchfort_inference(model, input, output)
      else
        !istat = torchfort_inference(model, get_ptr_gpu_training(iux,iuz), &
        !                                   get_ptr_gpu_training(itauxx,itauzz))

        start_time = mpiwtime()
        call infer_gpu(itsub)
        end_time = mpiwtime()
        inference_time = inference_time + end_time-start_time

      endif

      if (istat /= TORCHFORT_RESULT_SUCCESS) &
        call fatal_error("infer","istat="//trim(itoa(istat)))

    endsubroutine infer
!***************************************************************
    subroutine scale(f, minvalue, maxvalue)

      real, dimension (:,:,:,:) :: f
      real :: minvalue, maxvalue

      f = (f - minvalue)/(maxvalue - minvalue)

    endsubroutine
!***************************************************************
    subroutine descale(f, minvalue, maxvalue)

      real, dimension (:,:,:,:) :: f
      real :: minvalue, maxvalue

      f = f*(maxvalue - minvalue) + minvalue

    endsubroutine
!***************************************************************
    subroutine save_model
      if (iproc == root) print*,"Saving model to ",trim(model_output_dir)//trim(model_file)
       istat = torchfort_save_model(model, trim(model_output_dir)//trim(model_file))
       if (istat /= TORCHFORT_RESULT_SUCCESS) &
         call fatal_error("save_model","when saving model: istat="//trim(itoa(istat)))
       lmodel_saved = .true.
      if (iproc == root) print*,"Model saved to ",trim(model_output_dir)//trim(model_file)
    endsubroutine save_model
!***************************************************************
    subroutine train(f)
   
      use Gpu, only: get_ptr_gpu_training, train_gpu, infer_gpu
      use Mpicomm, only: mpiwtime
      real, save :: t_last_train = 0.0
  

      real, dimension (mx,my,mz,mfarray) :: f
    
      logical :: ldo_train = .false.

      if (it<it_train_start) return


      if((t-t_last_train) > dt_train) then
        t_last_train = t
        ldo_train = .true.
      endif


      if ((it_train /= -1).and.mod(it,it_train)==0) then
        ldo_train = .true.
      endif

      if(ldo_train) then
!
        if (.not. lfortran_launched) then
          !istat = torchfort_train(model, get_ptr_gpu_training(iux,iuz), &
          !                               get_ptr_gpu_training(itauxx,itauzz), train_loss)
          start_time = mpiwtime()
          call train_gpu(train_loss, itsub)
          end_time = mpiwtime()
          training_time = training_time + end_time-start_time
        else
          call calc_tau(f)
!
!  input scaling.
!
          if (lscale) then
            if (it == it_train_start) then
              input_min = minval(uumean)
              input_max = maxval(uumean)
            endif

            call scale(uumean, input_min, input_max)
!
! output scaling.
!
            if (it == it_train_start) then
              output_min = minval(f(:,:,:,itauxx:itauzz))
              output_max = maxval(f(:,:,:,itauxx:itauzz))
            endif
            call scale(f(:,:,:,itauxx:itauzz), output_min, output_max)
          endif

          ! print*, output_min, output_max, input_min, input_max
          input(:,:,:,:,1) = uumean                    ! host to device    !sngl(uumean)
          label(:,:,:,:,1) = f(:,:,:,itauxx:itauzz)    ! host to device

          istat = torchfort_train(model, input, label, train_loss)
!print*, 'TRAIN', it, train_loss

        endif

        if (istat /= TORCHFORT_RESULT_SUCCESS) call fatal_error("train","istat="//trim(itoa(istat)))

        if (train_loss <= max_loss) ltrained=.true.
        if ((it_train_end >= 0) .and. it >= it_train_end) ltrained=.true.
        if(ltrained) then
                call save_model
        else if (lroot.and.lfirst.and.mod(it,it_train_chkpt)==0) then
          istat = torchfort_save_checkpoint(trim(model), trim(checkpoint_output_dir))
          if (istat /= TORCHFORT_RESULT_SUCCESS) &
            call fatal_error("train","when saving checkpoint: istat="//trim(itoa(istat)))
          lckpt_written = .true.
          print*, 'it,it_train_chkpt=', it,it_train_chkpt, trim(model),istat, trim(checkpoint_output_dir), lckpt_written
        endif
      endif

    endsubroutine train
!***************************************************************
    subroutine calc_tau(f)

      use Sub, only: smooth

      real, dimension (mx,my,mz,mfarray) :: f
!
!  Smooth velocity.
!
      uumean = f(:,:,:,iux:iuz)
      call smooth(uumean,1,3,lgauss=.true.)
!
!  Calculate and smooth stress tensor.
!
      f(:,:,:,itauxx) = f(:,:,:,iux)**2
      f(:,:,:,itauyy) = f(:,:,:,iuy)**2
      f(:,:,:,itauzz) = f(:,:,:,iuz)**2
      f(:,:,:,itauxy) = f(:,:,:,iux)*f(:,:,:,iuy)
      f(:,:,:,itauyz) = f(:,:,:,iuy)*f(:,:,:,iuz)
      f(:,:,:,itauxz) = f(:,:,:,iux)*f(:,:,:,iuz)

      call smooth(f,itauxx,itauzz, lgauss=.true.)
!
!  Substract stresses from mean velocity.
!
      f(:,:,:,itauxx) = -uumean(:,:,:,1)**2 + f(:,:,:,itauxx)
      f(:,:,:,itauyy) = -uumean(:,:,:,2)**2 + f(:,:,:,itauyy)
      f(:,:,:,itauzz) = -uumean(:,:,:,3)**2 + f(:,:,:,itauzz)
      f(:,:,:,itauxy) = -uumean(:,:,:,1)*uumean(:,:,:,2) + f(:,:,:,itauxy)
      f(:,:,:,itauyz) = -uumean(:,:,:,2)*uumean(:,:,:,3) + f(:,:,:,itauyz)
      f(:,:,:,itauxz) = -uumean(:,:,:,1)*uumean(:,:,:,3) + f(:,:,:,itauxz)

    endsubroutine calc_tau
!***************************************************************
    subroutine div_reynolds_stress(f,df)

      use Sub, only: div

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df

      real, dimension(nx,3) :: divrey

      if (ltrained) then 
        call div(f,0,divrey(:,1),inds=(/itauxx,itauxy,itauxz/))
        call div(f,0,divrey(:,2),inds=(/itauxy,itauyy,itauyz/))
        call div(f,0,divrey(:,3),inds=(/itauxz,itauyz,itauzz/))

        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - divrey
      endif

    endsubroutine div_reynolds_stress
!***************************************************************
    subroutine calc_diagnostics_training(f)

      use Diagnostics, only: sum_mn_name, save_name

      real, dimension (mx,my,mz,mfarray) :: f

      integer :: i,j,jtau

      if (ldiagnos) then
        if (ltrained.and.lfirstpoint) then
          if (idiag_tauerror/=0) call sum_mn_name(spread(tauerror,1,nx),idiag_tauerror,lsqrt=.true.)
        else
          call save_name(real(train_loss), idiag_loss)
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
!  Original tau field.
!
        case ('tauxx'); call assign_slices_scal(slices,f,itauxx)
        case ('tauxy'); call assign_slices_scal(slices,f,itauxy)
        case ('tauxz'); call assign_slices_scal(slices,f,itauxz)
        case ('tauyy'); call assign_slices_scal(slices,f,itauyy)
        case ('tauyz'); call assign_slices_scal(slices,f,itauyz)
        case ('tauzz'); call assign_slices_scal(slices,f,itauzz)
!
!  Predicted tau field.
!
        case ('taupredxx'); call assign_slices_scal(slices,tau_pred,1)
        case ('taupredyy'); call assign_slices_scal(slices,tau_pred,2)
        case ('taupredzz'); call assign_slices_scal(slices,tau_pred,3)
        case ('taupredxy'); call assign_slices_scal(slices,tau_pred,4)
        case ('taupredxz'); call assign_slices_scal(slices,tau_pred,5)
        case ('taupredyz'); call assign_slices_scal(slices,tau_pred,6)

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
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'taupredxx',idum)
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'taupredxy',idum)
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'taupredxz',idum)
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'taupredyy',idum)
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'taupredyz',idum)
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'taupredzz',idum)
      enddo

    endsubroutine rprint_training
!***************************************************************
    subroutine finalize_training
!
!  Save trained model.
!
      if (.not.lstart) then
!print*, 'ltrained .or. .not. lckpt_written=', ltrained, lckpt_written
        !TP: now redundant (since model is saved immediately after ltrained becomes true in train) but for now lets keep it for safety
        if (ltrained .and. .not.lmodel_saved) then
                call save_model
        endif
        if (lfortran_launched) deallocate(input,label,output)
      endif

    endsubroutine finalize_training
!***************************************************************
    subroutine write_sample(sample, mx, my, mz, fname)

      use HDF5

      character(len=*) :: fname
      integer, intent(in) :: mx, my, mz
      real, intent(in) :: sample(mx, my, mz)
      integer(HID_T) :: in_file_id
      integer(HID_T) :: out_file_id
      integer(HID_T) :: dset_id
      integer(HID_T) :: dspace_id
      integer(HSIZE_T) :: dims(size(shape(sample)))
      integer :: err
    
      call h5open_f(err)
      call h5fcreate_f (fname, H5F_ACC_TRUNC_F, out_file_id, err)
    
      dims = shape(sample)
      call h5screate_simple_f(size(shape(sample)), dims, dspace_id, err)
      call h5dcreate_f(out_file_id, "data", H5T_NATIVE_REAL, dspace_id, dset_id, err)
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, sample, dims, err)
      call h5dclose_f(dset_id, err)
      call h5sclose_f(dspace_id, err)
    
      call h5fclose_f(out_file_id, err)
      call h5close_f(err)

    endsubroutine write_sample
!***************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=50
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(itauxx,p_par(1)) ! int
    call copy_addr(itauxy,p_par(2)) ! int
    call copy_addr(itauxz,p_par(3)) ! int
    call copy_addr(itauyy,p_par(4)) ! int
    call copy_addr(itauyz,p_par(5)) ! int
    call copy_addr(itauzz,p_par(6)) ! int
    call copy_addr(lscale,p_par(7)) ! bool
    call copy_addr(ltrained,p_par(8)) ! bool


    endsubroutine pushpars2c
!***********************************************************************
  endmodule Training
