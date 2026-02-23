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
    real :: t_train_start = 0.0, t_train_end = -1.0, t_train_chkpt=-1.0

    !real(KIND=rkind4), dimension(:,:,:,:,:), allocatable, device :: input, label, output
    real, dimension(:,:,:,:,:), allocatable, device :: input, label, output
    real :: train_loss   !(KIND=rkind4) :: train_loss

    integer :: itau_density, itau_densityx, itau_densityy, itau_densityz
    integer :: itau_hydro, itau_hydroxx, itau_hydroxy, itau_hydroxz, itau_hydroyy, itau_hydroyz, itau_hydrozz
    integer :: isgs_emf, isgs_emfx, isgs_emfy, isgs_emfz

    character(LEN=fnlen) :: model='model', config_file="config_mlp_native.yaml", model_file

    logical :: lroute_via_cpu=.false., lfortran_launched, luse_trained_tau, lwrite_sample=.false., lscale=.true.
    real :: max_loss=1.e-4, dt_train=1.e-10

    integer :: idiag_loss=0            ! DIAG_DOC: torchfort training loss
    integer :: idiag_tauerror=0        ! DIAG_DOC: $\sqrt{\left<(\sum_{i,j} u_i*u_j - tau_{ij})^2\right>}$

    namelist /training_run_pars/ config_file, model, it_train, it_train_start, it_train_chkpt, &
                                 luse_trained_tau, lscale, lwrite_sample, max_loss, lroute_via_cpu,&
                                 it_train_end, lrun_epoch, dt_train, t_train_start, t_train_end, t_train_chkpt,&
                                 ltrain_mag,ltrain_dens
!
    character(LEN=fnlen) :: model_output_dir, checkpoint_output_dir
    integer :: istat, train_step_ckpt, val_step_ckpt
    logical :: ltrained=.false., lckpt_written=.false.,lmodel_saved=.false., lrun_epoch=.false.
    real, dimension (mx,my,mz,3) :: uumean
    real :: tauerror, input_min, input_max, output_min, output_max
    real, dimension(mx, my, mz, 6) :: tau_pred
    real::  start_time, end_time
    real, save :: t_last_chkpt = 0.0
    integer :: input_channels  = 3
    integer :: output_channels = 6
    !TP: these are by default false now
    logical :: ltrain_mag  = .false.
    logical :: ltrain_dens = .false.

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
      f(:,:,:,itau_hydroxx:itau_hydroyz)   = 0.0
      if (ltrain_dens) then 
        f(:,:,:,itau_densityx:itau_densityz) = 0.0
        input_channels  = input_channels  + 1
        output_channels = output_channels + 3
      endif
      if (ltrain_mag) then
        f(:,:,:,isgs_emfx:isgs_emfz)       = 0.0
        input_channels  = input_channels  + 3
        output_channels = output_channels + 3
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

      ltrain_mag  = ltrain_mag  .and. lmagnetic
      ltrain_dens = ltrain_dens .and. ldensity
!
      call farray_register_auxiliary('tau_hydro',itau_hydro,vector=6,on_gpu=lgpu,communicated=.true.)
      if(ltrain_mag) call farray_register_auxiliary('sgs_emf',isgs_emf,vector=3,on_gpu=lgpu,communicated=.true.)
      if(ltrain_dens)  call farray_register_auxiliary('tau_density',itau_density,vector=3,on_gpu=lgpu,communicated=.true.)
!
!  Indices to access tau.
!
      itau_hydroxx=itau_hydro; itau_hydroyy=itau_hydro+1; itau_hydrozz=itau_hydro+2; itau_hydroxy=itau_hydro+3; itau_hydroxz=itau_hydro+4; itau_hydroyz=itau_hydro+5

      isgs_emfx=isgs_emf; isgs_emfy=isgs_emf+1; isgs_emfz=isgs_emf+2;

      itau_densityx=itau_density; itau_densityy=itau_density+1; itau_densityz=itau_density+2;

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
            f(:,:,:,itau_hydroxx:itau_hydrozz) = f(:,:,:,itau_hydroxx:itau_hydrozz) - output(:,:,:,:,1)
          endif
          tauerror = sum(f(l1:l2,m1:m2,n1:n2,itau_hydroxx:itau_hydrozz)**2)/nx
        else
          f(:,:,:,itau_hydroxx:itau_hydrozz) = output(:,:,:,:,1)
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
          call write_sample(f(:,:,:,itau_hydroxx), mx, my, mz, "target_"//trim(itoa(iproc))//".hdf5")
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
        !                                   get_ptr_gpu_training(itau_hydroxx,itau_hydrozz))

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
      !SG: should save the model instead of checkpoint in the end
       istat = torchfort_save_model(model, trim(model_output_dir)//trim(model_file))
       if (istat /= TORCHFORT_RESULT_SUCCESS) &
         call fatal_error("save_model","when saving model: istat="//trim(itoa(istat)))
       lmodel_saved = .true.
      if (iproc == root) print*,"Model saved to ",trim(model_output_dir)//trim(model_file)
    endsubroutine save_model
!***************************************************************
    subroutine save_chkpt
       if (lroot.and.lfirst.and.((mod(it,it_train_chkpt)==0).or.((t-t_last_chkpt) >= t_train_chkpt))) then
       print*, 'model: ', trim(model)
        istat = torchfort_save_checkpoint(trim(model), trim(checkpoint_output_dir))
        t_last_chkpt=t
        if (istat /= TORCHFORT_RESULT_SUCCESS) &
          call fatal_error("train","when saving checkpoint: istat="//trim(itoa(istat)))
        lckpt_written = .true.
        !print*, 'it, it_train_chkpt: , t:, t_train_chkpt: , t_last_chkpt: ', it, it_train_chkpt, t, t_train_chkpt, t_last_chkpt, trim(model), istat, trim(checkpoint_output_dir), lckpt_written
      endif


    endsubroutine save_chkpt
!***************************************************************
    subroutine train(f)
   
      use Gpu, only: get_ptr_gpu_training, train_gpu, infer_gpu
      use Mpicomm, only: mpiwtime
  

      real, dimension (mx,my,mz,mfarray) :: f
    
      real, save :: t_last_train = 0.0
      logical :: ldo_training_step

      if (it<it_train_start .or. t<t_train_start) return

      ldo_training_step = &
              ((it_train /= -1) .and. mod(it,it_train)==0) .or. &
              ((t-t_last_train) >= dt_train)
      if(.not. ldo_training_step) return
      
      t_last_train = t
      if (.not. lfortran_launched) then
        !istat = torchfort_train(model, get_ptr_gpu_training(iux,iuz), &
        !                               get_ptr_gpu_training(itau_hydroxx,itau_hydrozz), train_loss)
        start_time = mpiwtime()
        call train_gpu(train_loss, itsub, t)
        end_time = mpiwtime()
        training_time = training_time + end_time-start_time
      else
        call calc_tau(f)
!
!  inp scaling.
!
        if (lscale) then
          if (it == it_train_start) then
            input_min = minval(uumean)
            input_max = maxval(uumean)
          endif

          call scale(uumean, input_min, input_max)
!
! outp scaling.
!
          if (it == it_train_start) then
            output_min = minval(f(:,:,:,itau_hydroxx:itau_hydrozz))
            output_max = maxval(f(:,:,:,itau_hydroxx:itau_hydrozz))
          endif
          call scale(f(:,:,:,itau_hydroxx:itau_hydrozz), output_min, output_max)
        endif

        ! print*, output_min, output_max, input_min, input_max
        input(:,:,:,:,1) = uumean                    ! host to device    !sngl(uumean)
        label(:,:,:,:,1) = f(:,:,:,itau_hydroxx:itau_hydrozz)    ! host to device

        istat = torchfort_train(model, input, label, train_loss)
!print 'TRAIN', it, train_loss

      endif

      if (istat /= TORCHFORT_RESULT_SUCCESS) call fatal_error("train","istat="//trim(itoa(istat)))

      if (train_loss <= max_loss) ltrained=.true.
      if ((it_train_end >= 0) .and. it >= it_train_end) ltrained=.true.
      if ((t_train_end >= 0) .and. t >= t_train_end) ltrained=.true.
      if(ltrained) then
            call save_model
      else
            call save_chkpt
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
      f(:,:,:,itau_hydroxx) = f(:,:,:,iux)**2
      f(:,:,:,itau_hydroyy) = f(:,:,:,iuy)**2
      f(:,:,:,itau_hydrozz) = f(:,:,:,iuz)**2
      f(:,:,:,itau_hydroxy) = f(:,:,:,iux)*f(:,:,:,iuy)
      f(:,:,:,itau_hydroyz) = f(:,:,:,iuy)*f(:,:,:,iuz)
      f(:,:,:,itau_hydroxz) = f(:,:,:,iux)*f(:,:,:,iuz)

      call smooth(f,itau_hydroxx,itau_hydrozz, lgauss=.true.)
!
!  Substract stresses from mean velocity.
!
      f(:,:,:,itau_hydroxx) = -uumean(:,:,:,1)**2 + f(:,:,:,itau_hydroxx)
      f(:,:,:,itau_hydroyy) = -uumean(:,:,:,2)**2 + f(:,:,:,itau_hydroyy)
      f(:,:,:,itau_hydrozz) = -uumean(:,:,:,3)**2 + f(:,:,:,itau_hydrozz)
      f(:,:,:,itau_hydroxy) = -uumean(:,:,:,1)*uumean(:,:,:,2) + f(:,:,:,itau_hydroxy)
      f(:,:,:,itau_hydroyz) = -uumean(:,:,:,2)*uumean(:,:,:,3) + f(:,:,:,itau_hydroyz)
      f(:,:,:,itau_hydroxz) = -uumean(:,:,:,1)*uumean(:,:,:,3) + f(:,:,:,itau_hydroxz)

    endsubroutine calc_tau
!***************************************************************
    subroutine div_sgs_stresses(f,df)

      use Sub, only: div_tensor

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df

      real, dimension(nx,3) :: div_hydro_sgs
      real, dimension(nx)   :: div_dens_sgs

      if (ltrained) then 
        call div_tensor(f,div_hydro_sgs,itau_hydro)
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - div_hydro_sgs

        if(ltrain_mag) then
          df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + f(l1:l2,m,n,isgs_emf)
        endif

        if(ltrain_dens) then
          call div(f,itau_hydro,div_dens_sgs)
          df(l1:l2,m,n,ilnrho)  = df(l1:l2,m,n,ilnrho)  - div_dens_sgs
        endif
      endif

    endsubroutine div_sgs_stresses
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
        case ('tauxx'); call assign_slices_scal(slices,f,itau_hydroxx)
        case ('tauxy'); call assign_slices_scal(slices,f,itau_hydroxy)
        case ('tauxz'); call assign_slices_scal(slices,f,itau_hydroxz)
        case ('tauyy'); call assign_slices_scal(slices,f,itau_hydroyy)
        case ('tauyz'); call assign_slices_scal(slices,f,itau_hydroyz)
        case ('tauzz'); call assign_slices_scal(slices,f,itau_hydrozz)
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
        if (ltrained .or. .not.lmodel_saved) then
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

    call copy_addr(itau_hydroxx,p_par(1)) ! int
    call copy_addr(itau_hydroxy,p_par(2)) ! int
    call copy_addr(itau_hydroxz,p_par(3)) ! int
    call copy_addr(itau_hydroyy,p_par(4)) ! int
    call copy_addr(itau_hydroyz,p_par(5)) ! int
    call copy_addr(itau_hydrozz,p_par(6)) ! int
    call copy_addr(lscale,p_par(7)) ! bool
    call copy_addr(ltrained,p_par(8)) ! bool
    call copy_addr(itau_magxx,p_par(9)) ! int
    call copy_addr(itau_magxy,p_par(10)) ! int
    call copy_addr(itau_magxz,p_par(11)) ! int
    call copy_addr(itau_magyy,p_par(12)) ! int
    call copy_addr(itau_magyz,p_par(13)) ! int
    call copy_addr(itau_magzz,p_par(14)) ! int
    call copy_addr(itau_densityx,p_par(15)) ! int
    call copy_addr(itau_densityy,p_par(16)) ! int
    call copy_addr(itau_densityz,p_par(17)) ! int
    call copy_addr(input_channels,p_par(18)) ! int
    call copy_addr(output_channels,p_par(19)) ! int

    endsubroutine pushpars2c
!***********************************************************************
  endmodule Training
