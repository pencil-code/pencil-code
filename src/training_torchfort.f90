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
    use Cudafor, only: cudaSetDevice,CUDASUCCESS
    use Torchfort, only: torchfort_create_distributed_model, torchfort_create_model,&
                         torchfort_result_success,torchfort_load_model,torchfort_load_checkpoint,&
                         torchfort_save_model,torchfort_result_success,torchfort_save_checkpoint,&
                         torchfort_inference,torchfort_train
    !use iso_c_binding


    implicit none
    include 'training.h'

    include "training.h"

    integer :: model_device=0
    integer :: it_train=-1, it_train_chkpt=-1, it_train_start=1

    !real(KIND=rkind4), dimension(:,:,:,:,:), allocatable, device :: input, label, output
    real, dimension(:,:,:,:,:), allocatable, device :: input, label, output
    real :: train_loss   !(KIND=rkind4) :: train_loss

    integer :: itau, itauxx, itauxy, itauxz, itauyy, itauyz, itauzz

    character(LEN=fnlen) :: model='model', config_file="config_mlp_native.yaml", model_file

    logical :: lroute_via_cpu, lcpu_training=.true., luse_trained_tau, lwrite_sample=.false., lscale=.true.
    real :: max_loss=1.e-4

    integer :: idiag_loss=0            ! DIAG_DOC: torchfort training loss
    integer :: idiag_tauerror=0        ! DIAG_DOC: $\sqrt{\left<(\sum_{i,j} u_i*u_j - tau_{ij})^2\right>}$

    namelist /training_run_pars/ config_file, model, it_train, it_train_start, it_train_chkpt, &
                                 luse_trained_tau, lscale, lwrite_sample, max_loss, lroute_via_cpu
!
    character(LEN=fnlen) :: model_output_dir, checkpoint_output_dir
    integer :: istat, train_step_ckpt, val_step_ckpt
    logical :: ltrained=.false., lckpt_written=.false.
    real, dimension (mx,my,mz,3) :: uumean
    real :: tauerror, input_min, input_max, output_min, output_max

    contains
!***************************************************************
    subroutine initialize_training

      use File_IO, only: file_exists
      use Mpicomm, only: mpibcast, MPI_COMM_WORLD
      use Syscalls, only: system_cmd

      character(LEN=fnlen) :: modelfn

      lcpu_training = .not. lgpu .or. lroute_via_cpu
      if (lreloading) return

      if (.not.lhydro) call fatal_error('initialize_training','needs HYDRO module')
      istat = cudaSetDevice(iproc)
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
print*, 'ltrained, modelfn=', ltrained, modelfn
        endif
      endif
      call mpibcast(ltrained)
!
! TorchFort create model
!
      print*, 'CONFIG FILE=', trim(model_output_dir)//trim(config_file)
      if (lmpicomm) then
        istat = torchfort_create_distributed_model(trim(model), trim(model_output_dir)//trim(config_file), MPI_COMM_WORLD, iproc)
      else
        istat = torchfort_create_model(trim(model), trim(model_output_dir)//trim(config_file), model_device)
      endif
      if (istat /= TORCHFORT_RESULT_SUCCESS) then
        call fatal_error("initialize_training","when creating model "//trim(model)//": istat="//trim(itoa(istat)))
      else
        call information('initialize_training','TORCHFORT LIB LOADED SUCCESFULLY')
      endif

      if (ltrained) then
        istat = torchfort_load_model(trim(model), trim(modelfn))
        if (istat /= TORCHFORT_RESULT_SUCCESS) then
          call fatal_error("initialize_training","when loading model: istat="//trim(itoa(istat)))
        else
          call information('initialize_training','TORCHFORT MODEL "'//trim(modelfn)//'" LOADED SUCCESFULLY')
        endif
      else
        !if (file_exists(trim(checkpoint_output_dir)//'/'//trim(model)//'.pt')) then
        if (file_exists(trim(checkpoint_output_dir)//'/model.pt')) then

          !istat = torchfort_load_checkpoint(trim(model), trim(checkpoint_output_dir), train_step_ckpt, val_step_ckpt)
          istat = torchfort_load_checkpoint('model', trim(checkpoint_output_dir), train_step_ckpt, val_step_ckpt)
          if (istat /= TORCHFORT_RESULT_SUCCESS) then
            call fatal_error("initialize_training","when loading checkpoint: istat="//trim(itoa(istat)))
          else
            call information('initialize_training','TORCHFORT CHECKPOINT LOADED SUCCESFULLY')
          endif

        endif
      endif

      luse_trained_tau = luse_trained_tau.and.ltrained

      if (lrun .and. lcpu_training) then
        allocate(input(mx, my, mz, 3, 1))
        allocate(output(mx, my, mz, 6, 1))
        allocate(label(mx, my, mz, 6, 1))
print*, 'adresses:', loc(input), loc(output), loc(label)
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
          ! Device to host
        if ((ldiagnos.or.lvideo).and.lfirst) then
          call calc_tau(f)
          if (.not. lcpu_training) then
          else
            f(:,:,:,itauxx:itauzz) = f(:,:,:,itauxx:itauzz) - output(:,:,:,:,1)
          endif
          tauerror = sum(f(l1:l2,m1:m2,n1:n2,itauxx:itauzz)**2)/nx
        else
          f(:,:,:,itauxx:itauzz) = output(:,:,:,:,1)
        endif
      else
        if (lfirst) call train(f)
      endif

    endsubroutine training_before_boundary
!***************************************************************
    subroutine infer(f)
    
      use Gpu, only: get_ptr_gpu_for_training
      use Sub, only: smooth

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:,:,:,:), pointer :: ptr_uu, ptr_tau

      ! Host to device
      if (lcpu_training) then
        uumean = f(:,:,:,iux:iuz)
        call smooth(uumean,1,3,lgauss=.true.)
        input(:,:,:,:,1) = uumean    ! host to device
        istat = torchfort_inference(model, input, output)
      else
        !call get_ptr_gpu(ptr_uu,iux,iuz)
        !call get_ptr_gpu(ptr_tau,tauxx,tauzz)
        istat = torchfort_inference(model, get_ptr_gpu_for_training(iux,iuz), get_ptr_gpu_for_training(itauxx,itauzz))
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
    subroutine train(f)
   
      use Gpu, only: get_ptr_gpu_for_training

      real, dimension (mx,my,mz,mfarray) :: f

      if (it<it_train_start) return

      if (mod(it,it_train)==0) then
!
!  Smooth velocity.
!
        if (.not. lcpu_training) then
          !TODO: smoothing/scaling etc. for uu and tau
          istat = torchfort_train(model, get_ptr_gpu_for_training(iux,iuz), get_ptr_gpu_for_training(itauxx,itauzz), train_loss)
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
!print*, 'TRAIN', it, train_loss!   ltrained!, input_min, input_max
!
! output for plotting
!
          if (lwrite_sample .and. mod(it, 50)==0) then
            call write_sample(f(:,:,:,itauxx), mx, my, mz, "target_"//trim(itoa(iproc))//".hdf5")
            istat = torchfort_inference(model, input, output)
            f(:,:,:,itauxx:itauzz) = output(:,:,:,:,1)
            if (lscale) call descale(f(:,:,1:1,itauxx:itauzz), output_min, output_max)
            call write_sample(f(:,:,1:1,itauxx:itauzz), mx, my, mz, "pred_"//trim(itoa(iproc))//".hdf5")
          endif

        endif

        if (istat /= TORCHFORT_RESULT_SUCCESS) call fatal_error("train","istat="//trim(itoa(istat)))

        if (train_loss <= max_loss) ltrained=.true.
        if (lroot.and.lfirst.and.mod(it,it_train_chkpt)==0) then
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

      if (luse_trained_tau) then 
        call div(f,itauxx,divrey(:,1))
        call div(f,0,divrey(:,2),inds=(/itauxy,itauyy,itauyz/))
        call div(f,0,divrey(:,3),inds=(/itauxz,itauyz,itauzz/))

        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - divrey
      endif

    endsubroutine div_reynolds_stress
!***************************************************************
    subroutine calc_diagnostics_training(f,p)

      use Diagnostics, only: sum_mn_name, save_name

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      integer :: i,j,jtau

      if (ldiagnos) then
        if (ltrained.and.lfirstpoint) then
          call sum_mn_name(spread(tauerror,1,nx),idiag_tauerror,lsqrt=.true.)
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
!  Velocity field.
!
        case ('tauxxerr'); call assign_slices_scal(slices,f,itauxx)
        case ('tauxyerr'); call assign_slices_scal(slices,f,itauxy)
        case ('tauxzerr'); call assign_slices_scal(slices,f,itauxz)
        case ('tauyyerr'); call assign_slices_scal(slices,f,itauyy)
        case ('tauyzerr'); call assign_slices_scal(slices,f,itauyz)
        case ('tauzzerr'); call assign_slices_scal(slices,f,itauzz)

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
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauxxerr',idum)
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauxyerr',idum) 
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauxzerr',idum) 
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauyyerr',idum) 
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauyzerr',idum) 
        idum=0
        call parse_name(inamev,cnamev(inamev),cformv(inamev),'tauzzerr',idum) 
      enddo

    endsubroutine rprint_training
!***************************************************************
    subroutine finalize_training

!  Saving trained model.
      if (.not.lstart) then
print*, 'ltrained .or. .not. lckpt_written=', ltrained, lckpt_written
        if (ltrained .or. lckpt_written) then
          istat = torchfort_save_model(model, trim(model_output_dir)//trim(model_file))
          if (istat /= TORCHFORT_RESULT_SUCCESS) &
            call fatal_error("finalize_training","when saving model: istat="//trim(itoa(istat)))
        endif
        if (lcpu_training) deallocate(input,label,output)
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
  endmodule Training
