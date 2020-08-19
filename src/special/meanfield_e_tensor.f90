! $Id$
!
! This special module will contribute a mean electromotive
! force to the vector potential differential equation.
! 
! Calculation is done term by term, where each contribution
! is calculated using tensors supplied in a HDF5 file emftensors.h5.
!
! Simo Tuomisto, simo.tuomisto@aalto.fi
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED alpha_coefs(3,3); beta_coefs(3,3); gamma_coefs(3)
! PENCILS PROVIDED delta_coefs(3); kappa_coefs(3,3,3); utensor_coefs(3)
! PENCILS PROVIDED acoef_coefs(3,3); bcoef_coefs(3,3,3)
! PENCILS PROVIDED alpha_emf(3); beta_emf(3); gamma_emf(3)
! PENCILS PROVIDED delta_emf(3); kappa_emf(3); utensor_emf(3)
! PENCILS PROVIDED acoef_emf(3); bcoef_emf(3)
! PENCILS PROVIDED bij_symm(3,3)
! PENCILS PROVIDED emf(3)
!
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata
  use Diagnostics
  use General, only: keep_compiler_quiet,numeric_precision
  use Messages, only: svn_id, fatal_error, warning
  use Mpicomm, only: mpibarrier,MPI_COMM_WORLD,MPI_INFO_NULL,mpireduce_min, mpireduce_max

  use Sub, only: dot_mn, dot_mn_vm, curl_mn, cross_mn, vec_dot_3tensor,dot2_mn
  use HDF5
  use File_io, only: parallel_unit
!
  implicit none
!
  include '../special.h'

  real, dimension(nx) :: diffus_special, advec_special

  ! HDF debug parameters:

  integer :: hdferr
  logical :: hdf_exists

  ! Main HDF object ids and variables

  character (len=fnlen) :: hdf_emftensors_filename, emftensors_file='emftensors.h5'
  integer(HID_T) :: hdf_memtype, &
                    hdf_emftensors_file,  &
                    hdf_emftensors_plist, &
                    hdf_emftensors_group, &
                    hdf_grid_group

  ! Dataset HDF object ids

  integer, parameter :: ntensors = 8
  integer, parameter :: id_record_SPECIAL_ILOAD=1
  integer(HID_T),   dimension(ntensors) :: tensor_id_G, tensor_id_D, &
                                           tensor_id_S, tensor_id_memS
  integer(HSIZE_T), dimension(ntensors,10) :: tensor_dims, &
                                              tensor_offsets, tensor_counts, &
                                              tensor_memdims, &
                                              tensor_memoffsets, tensor_memcounts
  integer, parameter :: nscalars = 1
  integer(HID_T),   dimension(nscalars) :: scalar_id_D, scalar_id_S
  integer(HSIZE_T), dimension(nscalars) :: scalar_dims, &
                                           scalar_offsets, scalar_counts, &
                                           scalar_memdims, &
                                           scalar_memoffsets, scalar_memcounts
                                      
  ! Actual datasets

  real, dimension(:,:,:,:,:,:)  , allocatable :: alpha_data, beta_data, acoef_data
  real, dimension(:,:,:,:,:)    , allocatable :: gamma_data, delta_data, utensor_data
  real, dimension(:,:,:,:,:,:,:), allocatable :: kappa_data, bcoef_data
 
  real, dimension(:), allocatable :: tensor_times

  logical, dimension(3,3)::lalpha_sym=reshape((/.false.,.true. ,.false., &
                                                .true. ,.false.,.true. , &
                                                .false.,.true. ,.false. /), shape(lalpha_sym)), &
                           lbeta_sym =reshape((/.true. ,.false.,.true. , &
                                                .false.,.true. ,.false., &
                                                .true. ,.false.,.true.  /), shape(lbeta_sym))
  logical, dimension(3) :: lgamma_sym  =[.true. ,.false.,.true. ], &
                           ldelta_sym  =[.false.,.true. ,.false.], &
                           lutensor_sym=[.true. ,.false.,.true. ]
  !logical, dimension(3,3,3)::lkappa_sym=[[.false.,.true. ,.false.],[.true. ,.false.,.true. ],[.false.,.true.,.false.], &
  !                                       [.true. ,.false.,.true. ],[.false.,.true. ,.false.],[.true. ,.false.,.true.], &
  !                                       [.false.,.true. ,.false.],[.true. ,.false.,.true. ],[.false.,.true.,.false.]]
  logical, dimension(3,3,3)::lkappa_sym=reshape((/.false.,.true. ,.false.,.true. ,.false.,.true. ,.false.,.true.,.false., &
                                                  .true. ,.false.,.true. ,.false.,.true. ,.false.,.true. ,.false.,.true., &
                                                  .false.,.true. ,.false.,.true. ,.false.,.true. ,.false.,.true.,.false./), &
                                                shape(lkappa_sym))
  ! Dataset mappings

  integer,parameter :: alpha_id=1, beta_id=2,    &
                       gamma_id=3, delta_id=4,   &
                       kappa_id=5, utensor_id=6, & 
                       acoef_id=7, bcoef_id=8
  integer,parameter :: time_id=1

  integer, dimension(ntensors),parameter :: tensor_ndims = (/ 6, 6, 5, 5, 7, 5, 6, 7 /) 
  integer, dimension(nscalars),parameter :: scalar_ndims = (/ 1 /) 

  ! Dataset logical variables

  logical, dimension(3,3)  :: lalpha_arr, lbeta_arr, lacoef_arr
  logical, dimension(3)    :: lgamma_arr, ldelta_arr, lutensor_arr
  logical, dimension(3,3,3):: lkappa_arr, lbcoef_arr
  logical, dimension(6)    :: lalpha_c, lbeta_c, lacoef_c
  logical, dimension(3)    :: lgamma_c, ldelta_c, lutensor_c, lmeanuu_c
  logical, dimension(3,6)  :: lkappa_c, lbcoef_c
  logical :: lalpha=.true., lbeta=.true., lgamma=.true., ldelta=.true., lkappa=.true.
  logical :: lutensor=.true., lmeanuu=.true., lacoef=.false., lbcoef=.false., lusecoefs=.false.
  logical :: lread_datasets=.true., lread_time_series=.false., lloop=.false.
  real :: alpha_scale, beta_scale, gamma_scale, delta_scale, kappa_scale, utensor_scale, meanuu_scale, acoef_scale, bcoef_scale
  character (len=fnlen) :: defaultname
  character (len=fnlen) :: alpha_name, beta_name,    &
                           gamma_name, delta_name,   &
                           kappa_name, utensor_name, &
                           acoef_name, bcoef_name
  character (len=fnlen), dimension(ntensors) :: tensor_names
  character (len=fnlen), dimension(nscalars) :: scalar_names=(/'t'/)
  real, dimension(ntensors) :: tensor_scales, tensor_maxvals, tensor_minvals

  ! Diagnostic variables
  ! Max diagnostics

  integer :: idiag_alphaxmax=0
  integer :: idiag_alphaymax=0
  integer :: idiag_alphazmax=0
  integer :: idiag_betaxmax=0
  integer :: idiag_betaymax=0
  integer :: idiag_betazmax=0
  integer :: idiag_gammaxmax=0
  integer :: idiag_gammaymax=0
  integer :: idiag_gammazmax=0
  integer :: idiag_deltaxmax=0
  integer :: idiag_deltaymax=0
  integer :: idiag_deltazmax=0
  integer :: idiag_kappaxmax=0
  integer :: idiag_kappaymax=0
  integer :: idiag_kappazmax=0
  integer :: idiag_utensorxmax=0
  integer :: idiag_utensorymax=0
  integer :: idiag_utensorzmax=0
  integer :: idiag_acoefxmax=0
  integer :: idiag_acoefymax=0
  integer :: idiag_acoefzmax=0
  integer :: idiag_bcoefxmax=0
  integer :: idiag_bcoefymax=0
  integer :: idiag_bcoefzmax=0
  integer :: idiag_emfxmax=0
  integer :: idiag_emfymax=0
  integer :: idiag_emfzmax=0
  integer :: idiag_emfcoefxmax=0
  integer :: idiag_emfcoefymax=0
  integer :: idiag_emfcoefzmax=0
  integer :: idiag_emfdiffmax=0
  integer :: idiag_emfxdiffmax=0
  integer :: idiag_emfydiffmax=0
  integer :: idiag_emfzdiffmax=0
! RMS diagnostics
  integer :: idiag_alpharms=0
  integer :: idiag_betarms=0
  integer :: idiag_gammarms=0
  integer :: idiag_deltarms=0
  integer :: idiag_kapparms=0
  integer :: idiag_utensorrms=0
  integer :: idiag_acoefrms=0
  integer :: idiag_bcoefrms=0
  integer :: idiag_emfrms=0
  integer :: idiag_emfcoefrms=0
  integer :: idiag_emfdiffrms=0
! timestep diagnostics
  integer :: idiag_dtemf_ave=0
  integer :: idiag_dtemf_dif=0
!
!  2D diagnostics
!
  integer :: idiag_emfxmxy=0, idiag_emfymxy=0, idiag_emfzmxy=0
  integer :: idiag_emfcoefxmxy=0, idiag_emfcoefymxy=0, idiag_emfcoefzmxy=0
!
  ! Interpolation parameters
  character (len=fnlen) :: interpname
  integer :: dataload_len=1, tensor_times_len=-1, iload=0
  real, dimension(nx)     :: tmpline
  real, dimension(nx,3)   :: tmppencil,emftmp
  real, dimension(nx,3,3) :: tmptensor
!
! special symmetries
!
  logical :: lsymmetrize=.false.
  integer :: nsmooth_rbound=0, nsmooth_thbound=0
  integer :: field_symmetry=0
!
! regularize beta
!
  logical :: lregularize_beta=.false., &
             lreconstruct_tensors=.false., &
             lalt_decomp=.false.,&
             lremove_beta_negativ=.false.
  real :: rel_eta=.01    ! should be > 0

  ! Input dataset name
  ! Input namelist

  namelist /special_init_pars/ &
      emftensors_file, &
      lalpha,   lalpha_c,   alpha_name,   alpha_scale, &
      lbeta,    lbeta_c,    beta_name,    beta_scale, &
      lgamma,   lgamma_c,   gamma_name,   gamma_scale, &
      ldelta,   ldelta_c,   delta_name,   delta_scale, &
      lkappa,   lkappa_c,   kappa_name,   kappa_scale, &
      lutensor, lutensor_c, utensor_name, utensor_scale, &
      lmeanuu,  lmeanuu_c,                meanuu_scale, &
      lacoef,   lacoef_c,   acoef_name,   acoef_scale, &
      lbcoef,   lbcoef_c,   bcoef_name,   bcoef_scale, &
      interpname, defaultname, lusecoefs, lloop
  namelist /special_run_pars/ &
      emftensors_file, &
      lalpha,   lalpha_c,   alpha_name,   alpha_scale, &
      lbeta,    lbeta_c,    beta_name,    beta_scale, &
      lgamma,   lgamma_c,   gamma_name,   gamma_scale, &
      ldelta,   ldelta_c,   delta_name,   delta_scale, &
      lkappa,   lkappa_c,   kappa_name,   kappa_scale, &
      lutensor, lutensor_c, utensor_name, utensor_scale, &
      lmeanuu,  lmeanuu_c,                meanuu_scale, &
      lacoef,   lacoef_c,   acoef_name,   acoef_scale, &
      lbcoef,   lbcoef_c,   bcoef_name,   bcoef_scale, &
      interpname, defaultname, lusecoefs, lloop, lsymmetrize, field_symmetry, &
      nsmooth_rbound, nsmooth_thbound, lregularize_beta, lreconstruct_tensors, &
      lalt_decomp, lremove_beta_negativ, rel_eta

  interface loadDataset
    module procedure loadDataset_rank1
    module procedure loadDataset_rank2
    module procedure loadDataset_rank3
  end interface

  interface symmetrize
    module procedure symmetrize_3d
    module procedure symmetrize_4d
  end interface
 
  interface smooth_thbound
    module procedure smooth_thbound_4d
  end interface

  interface smooth_rbound
    module procedure smooth_rbound_4d
  end interface

  contains
!***********************************************************************
    subroutine register_special
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      integer :: i
!
      if (lroot) call svn_id( &
           "$Id$")
!
      if (lrun) then 

        call H5open_F(hdferr)                                              ! Initializes HDF5 library.

        if (numeric_precision() == 'S') then
          if (lroot) write(*,*) 'register_special: loading data as single precision'
          hdf_memtype = H5T_NATIVE_REAL
        else
          if (lroot) write(*,*) 'register_special: loading data as double precision'
          hdf_memtype = H5T_NATIVE_DOUBLE
        end if
          
        if (lroot) write (*,*) 'register_special: setting parallel HDF5 IO for data file reading'   !MR: Why parallel?
        call H5Pcreate_F(H5P_FILE_ACCESS_F, hdf_emftensors_plist, hdferr)   ! Creates porperty list for HDF5 file.

        if (lmpicomm) &     !MR: doesn't work for nompicomm
          call H5Pset_fapl_mpio_F(hdf_emftensors_plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

        hdf_emftensors_filename = trim(datadir_snap)//'/'//trim(emftensors_file)

        if (lroot) print *, 'register_special: opening datafile '//trim(hdf_emftensors_filename)// &
                            ' and loading relevant fields into memory...'

        inquire(file=hdf_emftensors_filename, exist=hdf_exists)

        if (.not. hdf_exists) then                                        ! If HDF5 file doesn't exist:
          call H5Pclose_F(hdf_emftensors_plist, hdferr)                   ! Terminates access to property list.
          call H5close_F(hdferr)                                          ! Frees resources used by library.
          call fatal_error('register_special','File '//trim(hdf_emftensors_filename)//' does not exist!')
        end if
!
! Opens HDF5 file for read access only, returns file identifier hdf_emftensors_file.
!
        call H5Fopen_F(hdf_emftensors_filename, H5F_ACC_RDONLY_F, hdf_emftensors_file, hdferr, access_prp = hdf_emftensors_plist)
!
! Checks whether in HDF5 file there is a link /emftensor/.
!
        call H5Lexists_F(hdf_emftensors_file,'/emftensor/', hdf_exists, hdferr)

        if (.not. hdf_exists) then                                         ! If link '/emftensor/' doesn't exist:
          call H5Fclose_F(hdf_emftensors_file, hdferr)                     ! Terminates access to HDF5 file.
          call H5Pclose_F(hdf_emftensors_plist, hdferr)
          call H5close_F(hdferr)
          call fatal_error('register_special','group /emftensor/ does not exist!')
        end if

        call H5Gopen_F(hdf_emftensors_file, 'emftensor', hdf_emftensors_group, hdferr) ! Opens group emftensor in HDF5 file.

      endif
!!      call farray_register_pde('special',ispecial)
!!      call farray_register_auxiliary('specaux',ispecaux)
!!      call farray_register_auxiliary('specaux',ispecaux,communicated=.true.)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j
!
      call keep_compiler_quiet(f)

      if (lrun) then 

        if (trim(defaultname) == 'time-series' .or. trim(defaultname) == 'time-crop') then
          lread_time_series=.true.
        else
        !  call fatal_error('initialize_special','Unknown dataset chosen!')
        end if

        if (.not.lreloading) then

          if (lread_time_series) then 

            call H5Lexists_F(hdf_emftensors_file,'/grid/', hdf_exists, hdferr)
            if (.not. hdf_exists) then
              call H5Fclose_F(hdf_emftensors_file, hdferr)
              call H5Pclose_F(hdf_emftensors_plist, hdferr)
              call H5close_F(hdferr)
              call fatal_error('initialize_special','group /grid/ does not exist!')
            end if

            call H5Gopen_F(hdf_emftensors_file, 'grid', hdf_grid_group, hdferr)
            if (hdferr /= 0) call fatal_error('initialize_special','error while opening /grid/')

            call openDataset_grid(time_id)
            allocate(tensor_times(scalar_dims(time_id))) 
          
            if (hdferr /= 0) call fatal_error('initialize special','cannot select grid/t')

            call H5Dread_F(scalar_id_D(time_id), hdf_memtype, tensor_times, &
                           [scalar_dims(time_id)], hdferr)
            if (hdferr /= 0) call fatal_error('initialize special','cannot read grid/t')

            t = tensor_times(1)
            if (lroot) then
              print*, 'min time array=',minval(tensor_times)
              print*, 'max time array=',maxval(tensor_times)
            endif
   
          end if
        endif

        ! Set dataset offsets
        tensor_offsets = 0
        do i=1,ntensors
          tensor_offsets(i,1:4) = [ 0, ipx*nx, ipy*ny , ipz*nz ]
        end do

        ! Set dataset counts
        tensor_counts(:,5:) = 1
        do i=1,ntensors
          tensor_counts(i,1:4) = [ dataload_len, nx, ny, nz ]
        end do
        tensor_dims=tensor_counts
        tensor_dims(:,5:) = 3

        ! Set dataset memory offsets
        tensor_memoffsets = 0
        ! Set dataset memory counts
        tensor_memcounts = 1
        do i=1,ntensors
          tensor_memcounts(i,1:4) = [ dataload_len, nx, ny, nz ]
        end do

        ! Set memory dimensions
        tensor_memdims = 3
        do i=1,ntensors
          tensor_memdims(i,1:4) = [ dataload_len, nx, ny, nz ]
        end do

        if (lalpha) then
          if (.not.allocated(alpha_data)) then
            allocate(alpha_data(dataload_len,nx,ny,nz,3,3))
            call openDataset('alpha',alpha_id)
          endif
          alpha_data = 0.
        elseif (allocated(alpha_data)) then
          deallocate(alpha_data)
          call closeDataset(alpha_id)
        endif

        if (lbeta) then
          if (.not.allocated(beta_data)) then
            allocate(beta_data(dataload_len,nx,ny,nz,3,3))
            call openDataset('beta',beta_id)
          endif
          beta_data = 0.
        elseif (allocated(beta_data)) then
          deallocate(beta_data)
          call closeDataset(beta_id)
        endif

        if (lgamma) then
          if (.not.allocated(gamma_data)) then
            allocate(gamma_data(dataload_len,nx,ny,nz,3))
            call openDataset('gamma',gamma_id)
          endif
          gamma_data = 0.
        elseif (allocated(gamma_data)) then
          deallocate(gamma_data)
          call closeDataset(gamma_id)
        endif

        if (ldelta) then
          if (.not.allocated(delta_data)) then
            allocate(delta_data(dataload_len,nx,ny,nz,3))
            call openDataset('delta',delta_id)
          endif
          delta_data = 0.
        elseif (allocated(delta_data)) then
          deallocate(delta_data)
          call closeDataset(delta_id)
        endif

        if (lkappa) then
          if (.not.allocated(kappa_data)) then
            allocate(kappa_data(dataload_len,nx,ny,nz,3,3,3))
            call openDataset('kappa',kappa_id)
          endif
          kappa_data = 0.
        elseif (allocated(kappa_data)) then
          deallocate(kappa_data)
          call closeDataset(kappa_id)
        endif

        if (lutensor .or. lmeanuu) then
          if (.not.allocated(utensor_data)) then 
            allocate(utensor_data(dataload_len,nx,ny,nz,3))
            call openDataset('utensor',utensor_id)
          endif
          utensor_data = 0.
        elseif (allocated(utensor_data)) then
          deallocate(utensor_data)
          call closeDataset(utensor_id)
        endif

        if (lacoef) then
          if (.not.allocated(acoef_data)) then
            allocate(acoef_data(dataload_len,nx,ny,nz,3,3))
            call openDataset('acoef', acoef_id)
          endif
          acoef_data = 0.
        elseif (allocated(acoef_data)) then
          deallocate(acoef_data)
          call closeDataset(acoef_id)
        endif

        if (lbcoef) then
          if (.not.allocated(bcoef_data)) then 
            allocate(bcoef_data(dataload_len,nx,ny,nz,3,3,3))
            call openDataset('bcoef', bcoef_id)
          endif
          bcoef_data = 0.
        elseif (allocated(bcoef_data)) then
          deallocate(bcoef_data)
          call closeDataset(bcoef_id)
        endif

        if (.not.lusecoefs) then
          idiag_emfcoefxmxy=0; idiag_emfcoefymxy=0; idiag_emfcoefzmxy=0
        endif
        
        if (lroot) then
          if (lalpha) then
            write (*,*) 'Alpha scale:   ', alpha_scale
            write (*,*) 'Alpha components used: '
            do i=1,3
              write (*,'(A3,3L3,A3)') '|', lalpha_arr(:,i), '|'
            end do
          endif
          if (lbeta) then
            write (*,*) 'Beta scale:   ', beta_scale
            write (*,*) 'Beta components used: '
            do i=1,3
              write (*,'(A3,3L3,A3)') '|', lbeta_arr(:,i), '|'
            end do
          endif
          if (lgamma) then
            write (*,*) 'Gamma scale:   ', gamma_scale
            write (*,*) 'Gamma components used: '
            write (*,'(A3,3L3,A3)') '|', lgamma_arr, '|'
          endif
          if (ldelta) then
            write (*,*) 'Delta scale:   ', delta_scale
            write (*,*) 'Delta components used: '
            write (*,'(A3,3L3,A3)') '|', ldelta_arr, '|'
          endif
          if (lkappa) then
            write (*,*) 'Kappa scale:   ', kappa_scale
            write (*,*) 'Kappa components used: '
            do j=1,3
              write(*,*) '|'
              do i=1,3
                write (*,'(A3,3L3,A3)') '||', lkappa_arr(:,j,i), '||'
              end do
            end do
            write(*,*) '|'
          endif
          if (lutensor) then
            write (*,*) 'U-tensor scale:  ', utensor_scale
            write (*,*) 'U-tensor components used: '
            write (*,'(A3,3L3,A3)') '|', lutensor_arr, '|'
          endif
          if (lmeanuu) then
            write (*,*) 'Mean U scale:  ', meanuu_scale
            write (*,*) 'Mean U components used: '
            write (*,'(A3,3L3,A3)') '|', lutensor_arr, '|'
          endif
          if (lacoef) then
            write (*,*) 'acoef scale:   ', acoef_scale
            write (*,*) 'acoef components used: '
            do i=1,3
              write (*,'(A3,3L3,A3)') '|', lacoef_arr(:,i), '|'
            end do
          endif
          if (lbcoef) then
            write (*,*) 'bcoef scale:   ', bcoef_scale
            write (*,*) 'bcoef components used: '
            do j=1,3
              write(*,*) '|'
              do i=1,3
                write (*,'(A3,3L3,A3)') '||', lbcoef_arr(:,j,i), '||'
              end do
            end do
          endif
        end if

        lread_datasets=.true.

      endif
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine smooth_thbound_4d(arr,nsmooth)
!
!  Simple smothing at theta boundaries: last nsmooth points follow linear trend,
!  defined by last unsmoothed point and zero in boundary point.
!
!  20-jan-2019/MR: coded
!
      real, dimension(:,:,:,:), intent(INOUT) :: arr
      integer, intent(IN) :: nsmooth

      integer :: len_theta,len_r,ir,is
      real, dimension(size(arr,1),size(arr,4)) :: del

      len_theta=size(arr,3); len_r=size(arr,2)

      do ir=1,len_r

        if (lfirst_proc_y) then
          del=arr(:,ir,nsmooth+1,:)/nsmooth
          do is=1,nsmooth
            arr(:,ir,nsmooth+1-is,:)=arr(:,ir,nsmooth+1,:)-is*del
          enddo
          arr(:,ir,1,:)=0.
        endif

        if (llast_proc_y) then
          del=arr(:,ir,len_theta-nsmooth,:)/nsmooth
          do is=1,nsmooth
            arr(:,ir,len_theta-nsmooth+is,:)=arr(:,ir,len_theta-nsmooth,:)-is*del
          enddo
          arr(:,ir,len_theta,:)=0.
        endif

      enddo

    endsubroutine smooth_thbound_4d
!***********************************************************************
    subroutine smooth_rbound_4d(arr,nsmooth)
!
!  Simple smothing at r boundaries: last nsmooth points follow linear trend,
!  defined by last unsmoothed point and zero in boundary point.
!
!  20-jan-2019/MR: coded
!
      real, dimension(:,:,:,:), intent(INOUT) :: arr
      integer, intent(IN) :: nsmooth

      integer :: len_theta,len_r,ith,is
      real, dimension(size(arr,1),size(arr,4)) :: del

      len_theta=size(arr,3); len_r=size(arr,2)

      do ith=1,len_theta

        if (lfirst_proc_x) then
          del=arr(:,nsmooth+1,ith,:)/nsmooth
          do is=1,nsmooth
            arr(:,nsmooth+1-is,ith,:)=arr(:,nsmooth+1,ith,:)-is*del
          enddo
!          arr(:,ir,1,:)=0.
        endif

        if (llast_proc_x) then
          del=arr(:,len_r-nsmooth,ith,:)/nsmooth
          do is=1,nsmooth
            arr(:,len_theta-nsmooth+is,ith,:)=arr(:,len_theta-nsmooth,ith,:)-is*del
          enddo
!          arr(:,ir,len_theta,:)=0.
        endif

      enddo

    endsubroutine smooth_rbound_4d
!***********************************************************************
    subroutine symmetrize_4d(arr,lsym)

      use Mpicomm, only: mpisendrecv_real,mpibarrier,MPI_ANY_TAG
      use General, only: find_proc

      real, dimension(:,:,:,:), intent(INOUT) :: arr
      logical                 , intent(IN)    :: lsym

      integer :: len_theta,len_theta_h,symthproc
      integer, dimension(4) :: sz
      logical :: lmiddle
      real, dimension(:,:,:,:), allocatable :: buffer

      len_theta=size(arr,3); len_theta_h=floor(len_theta/2.)
      lmiddle=mod(nprocy,2)/=0.and.ipy==floor(nprocy/2.)
      sz=(/size(arr,1),size(arr,2),len_theta,size(arr,4)/)

      if (lmiddle) then

        if (lsym) then
          arr(:,:,:len_theta_h,:) = 0.5*(arr(:,:,:len_theta_h,:)+arr(:,:,len_theta:len_theta_h:-1,:))
          arr(:,:,len_theta:len_theta_h:-1,:) = arr(:,:,:len_theta_h,:)
        else
          arr(:,:,:len_theta_h,:) = 0.5*(arr(:,:,:len_theta_h,:)-arr(:,:,len_theta:len_theta_h:-1,:))
          arr(:,:,len_theta:len_theta_h:-1,:) = -arr(:,:,:len_theta_h,:)
        endif

      else

        allocate(buffer(sz(1),sz(2),sz(3),sz(4)))
        symthproc=find_proc(ipx,nprocy-1-ipy,ipz)
        call mpisendrecv_real(arr,sz,symthproc,iproc,buffer,symthproc,symthproc)

        if (lsym) then
          arr = 0.5*(arr+buffer(:,:,len_theta:1:-1,:))
        else
          arr = 0.5*(arr-buffer(:,:,len_theta:1:-1,:))
        endif

      endif
      call mpibarrier

    endsubroutine symmetrize_4d
!***********************************************************************
    subroutine symmetrize_3d(arr,lsym)

      use Mpicomm, only: mpisendrecv_real,mpibarrier,MPI_ANY_TAG
      use General, only: find_proc

      real, dimension(:,:,:), intent(INOUT) :: arr
      logical               , intent(IN)    :: lsym

      integer :: len_theta,len_theta_h,symthproc
      integer, dimension(3) :: sz
      logical :: lmiddle
      real, dimension(:,:,:), allocatable :: buffer

      len_theta=size(arr,2); len_theta_h=floor(len_theta/2.)
      lmiddle=mod(nprocy,2)/=0.and.ipy==floor(nprocy/2.)
      sz=(/size(arr,1),len_theta,size(arr,3)/)

      if (lmiddle) then

        if (lsym) then
          arr(:,:len_theta_h,:) = 0.5*(arr(:,:len_theta_h,:)+arr(:,len_theta:len_theta_h:-1,:))
          arr(:,len_theta:len_theta_h:-1,:) = arr(:,:len_theta_h,:)
        else
          arr(:,:len_theta_h,:) = 0.5*(arr(:,:len_theta_h,:)-arr(:,len_theta:len_theta_h:-1,:))
          arr(:,len_theta:len_theta_h:-1,:) = -arr(:,:len_theta_h,:)
        endif

      else

        allocate(buffer(sz(1),sz(2),sz(3)))
        symthproc=find_proc(ipx,nprocy-1-ipy,ipz)
        call mpisendrecv_real(arr,sz,symthproc,iproc,buffer,symthproc,MPI_ANY_TAG)  ! symthproc

        if (lsym) then
          arr = 0.5*(arr+buffer(:,len_theta:1:-1,:))
        else
          arr = 0.5*(arr-buffer(:,len_theta:1:-1,:))
        endif

      endif
      call mpibarrier

    endsubroutine symmetrize_3d
!***********************************************************************
    subroutine finalize_special(f)
!
!    Called right before exiting.
!
!    14-aug-2011/Bourdin.KIS: coded
!
        real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
        call keep_compiler_quiet(f)

      if (lrun) then

        ! Deallocate data
        if (allocated(alpha_data)) deallocate(alpha_data)
        if (allocated(beta_data)) deallocate(beta_data)
        if (allocated(gamma_data)) deallocate(gamma_data)
        if (allocated(delta_data)) deallocate(delta_data)
        if (allocated(kappa_data)) deallocate(kappa_data)
        if (allocated(utensor_data)) deallocate(utensor_data)

        if (lalpha)   call closeDataset(alpha_id)
        if (lbeta)    call closeDataset(beta_id)
        if (lgamma)   call closeDataset(gamma_id)
        if (ldelta)   call closeDataset(delta_id)
        if (lkappa)   call closeDataset(kappa_id)
        if (lutensor .or. lmeanuu) call closeDataset(utensor_id)
        if (lacoef)   call closeDataset(acoef_id)
        if (lbcoef)   call closeDataset(bcoef_id)

        if (lread_time_series) then
          if (allocated(tensor_times)) deallocate(tensor_times)
          call closeDataset_grid(time_id)
          call H5Gclose_F(hdf_grid_group, hdferr)
        endif

        call H5Gclose_F(hdf_emftensors_group, hdferr)
        call H5Fclose_F(hdf_emftensors_file, hdferr)
        call H5Pclose_F(hdf_emftensors_plist, hdferr)
        call H5close_F(hdferr)
        call mpibarrier

      end if
  !
    endsubroutine finalize_special
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-jul-06/tony: coded
!  20-may-19/MR: added beta regularization, reconstruction of alpha and gamma from acoef and bcoef,
!                alternative decomposition of tensors
!
      use SharedVariables, only: get_shared_variable
      use Mpicomm, only: mpireduce_sum_int, mpiallreduce_sum_int, mpireduce_min
      use PolynomialRoots, only: cubicroots
      use Sub, only: dyadic2_other

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f
!
      real :: delt,minbeta,minbeta_,thmin,thmax,rmin,rmax,trmin,det,oldtrace,newtrace
      integer :: i,j,k,numzeros,numcomplex,numzeros_,mm,ll,iv,iv0,ik
      real, pointer :: eta
      integer, dimension(:,:,:,:), allocatable :: beta_mask
      real, dimension(4) :: polcoeffs
      complex, dimension(3) :: eigenvals
      real, dimension(3,3) :: ev, beta_sav

      call keep_compiler_quiet(f)

      if (lfirst) then
        if (lread_time_series) then
          if (t >= tensor_times(iload+1)) then
            lread_datasets=.true.
            iload = iload + 1  
          endif
        else 
          iload=1
        endif
        
        if (lread_datasets) then
          if (iload > tensor_times_len) then 
            if (lloop) then
              iload=1
              delt=2*tensor_times(tensor_times_len)-tensor_times(1)-tensor_times(tensor_times_len-1)
              tensor_times = tensor_times+delt
            else
              call fatal_error('special_before_boundary', 'no more data to load') 
            endif
          endif
!if (lread_time_series) &
!print *, 'loading: t,iload=',tensor_times(iload),iload 
          ! Load datasets
          if (lalpha) call loadDataset(alpha_data, lalpha_arr, alpha_id, iload-1,'Alpha')
          if (lbeta)  call loadDataset(beta_data,  lbeta_arr,  beta_id,  iload-1,'Beta')
          if (lgamma) call loadDataset(gamma_data, lgamma_arr, gamma_id, iload-1,'Gamma')
          if (ldelta) call loadDataset(delta_data, ldelta_arr, delta_id, iload-1,'Delta')
          if (lkappa) call loadDataset(kappa_data, lkappa_arr, kappa_id, iload-1,'Kappa')
          if (lutensor .or. lmeanuu) call loadDataset(utensor_data, lutensor_arr, utensor_id, iload-1,'Utensor')
          if (lacoef) call loadDataset(acoef_data, lacoef_arr, acoef_id, iload-1,'Acoef')
          if (lbcoef) call loadDataset(bcoef_data, lbcoef_arr, bcoef_id, iload-1,'Bcoef')
          lread_datasets=.false.

          if (lreconstruct_tensors.and.lacoef.and.lbcoef) then                 ! beta, delta, kappa supposed to be correct
            if (lroot) print*, 'Reconstruct alpha and gamma from raw tensors.'
            do mm=1,ny
              alpha_data(1,:,mm,1,1,1)=     acoef_data(1,:,mm,1,1,1)+bcoef_data(1,:,mm,1,1,2,2)/x(l1:l2)
              alpha_data(1,:,mm,1,1,2)=0.5*(acoef_data(1,:,mm,1,1,2)-bcoef_data(1,:,mm,1,1,1,2)/x(l1:l2)  &
                                           +acoef_data(1,:,mm,1,2,1)+bcoef_data(1,:,mm,1,2,2,2)/x(l1:l2))
              alpha_data(1,:,mm,1,2,2)=     acoef_data(1,:,mm,1,2,2)-bcoef_data(1,:,mm,1,2,1,2)/x(l1:l2)
              alpha_data(1,:,mm,1,1,3)=0.5*(acoef_data(1,:,mm,1,1,3)+acoef_data(1,:,mm,1,3,1) & 
                                           +bcoef_data(1,:,mm,1,3,2,2)/x(l1:l2))
              alpha_data(1,:,mm,1,2,3)=0.5*(acoef_data(1,:,mm,1,2,3)+acoef_data(1,:,mm,1,3,2) &
                                           -bcoef_data(1,:,mm,1,3,1,2)/x(l1:l2))
            
              gamma_data(1,:,mm,1,1)=0.5*(acoef_data(1,:,mm,1,3,2)-acoef_data(1,:,mm,1,2,3) &
                                         -bcoef_data(1,:,mm,1,3,1,2)/x(l1:l2))
              gamma_data(1,:,mm,1,2)=0.5*(acoef_data(1,:,mm,1,1,3)-acoef_data(1,:,mm,1,3,1) &
                                         -bcoef_data(1,:,mm,1,3,2,2)/x(l1:l2))
              gamma_data(1,:,mm,1,3)=0.5*(acoef_data(1,:,mm,1,2,1)-acoef_data(1,:,mm,1,1,2) &
                                         +bcoef_data(1,:,mm,1,1,1,2)/x(l1:l2)+bcoef_data(1,:,mm,1,2,2,2)/x(l1:l2))
            enddo
            
            alpha_data(:,:,:,:,2,1)=alpha_data(:,:,:,:,1,2)
            alpha_data(:,:,:,:,3,2)=alpha_data(:,:,:,:,2,3)
            alpha_data(:,:,:,:,3,1)=alpha_data(:,:,:,:,1,3)
            alpha_data=alpha_data*alpha_scale
            gamma_data=gamma_data*gamma_scale

          endif

          if (lalt_decomp.and.lacoef.and.lbcoef) then                    ! kappa supposed to be correct.
            if (lroot) print*, 'Construct alternative decomposition from raw tensors.'
            do mm=1,ny

              alpha_data(1,:,mm,1,1,1)=      acoef_data(1,:,mm,1,1,1)+bcoef_data(1,:,mm,1,1,2,2)/x(l1:l2)
              alpha_data(1,:,mm,1,1,2)=0.5*( acoef_data(1,:,mm,1,1,2)-bcoef_data(1,:,mm,1,1,1,2)/x(l1:l2)  &
                                            +acoef_data(1,:,mm,1,2,1)+bcoef_data(1,:,mm,1,2,2,2)/x(l1:l2))
              alpha_data(1,:,mm,1,2,2)=      acoef_data(1,:,mm,1,2,2)-bcoef_data(1,:,mm,1,2,1,2)/x(l1:l2)

              alpha_data(1,:,mm,1,1,3)=0.5*( acoef_data(1,:,mm,1,1,3)+acoef_data(1,:,mm,1,3,1) & 
                                            +(                  bcoef_data(1,:,mm,1,3,2,2) &
                                              +                 bcoef_data(1,:,mm,1,1,3,1) &
                                              +cotth(mm+nghost)*bcoef_data(1,:,mm,1,1,3,2))/x(l1:l2))

              alpha_data(1,:,mm,1,2,3)=0.5*( acoef_data(1,:,mm,1,2,3)+acoef_data(1,:,mm,1,3,2) &
                                            +(                  bcoef_data(1,:,mm,1,2,3,1) &  
                                              -                 bcoef_data(1,:,mm,1,3,1,2) &
                                              +cotth(mm+nghost)*bcoef_data(1,:,mm,1,2,3,2))/x(l1:l2))

              alpha_data(1,:,mm,1,3,3)= acoef_data(1,:,mm,1,3,3) + (                  bcoef_data(1,:,mm,1,3,3,1) &
                                                                    +cotth(mm+nghost)*bcoef_data(1,:,mm,1,3,3,2))/x(l1:l2)

              gamma_data(1,:,mm,1,1)=0.5*( acoef_data(1,:,mm,1,3,2)-acoef_data(1,:,mm,1,2,3) &
                                          -(                  bcoef_data(1,:,mm,1,2,3,1) &
                                            +                 bcoef_data(1,:,mm,1,3,1,2) &
                                            +cotth(mm+nghost)*bcoef_data(1,:,mm,1,2,3,2))/x(l1:l2))

              gamma_data(1,:,mm,1,2)=0.5*( acoef_data(1,:,mm,1,1,3)-acoef_data(1,:,mm,1,3,1) &
                                          +(                  bcoef_data(1,:,mm,1,1,3,1) &
                                            -                 bcoef_data(1,:,mm,1,3,2,2) &
                                            +cotth(mm+nghost)*bcoef_data(1,:,mm,1,1,3,2))/x(l1:l2))

              gamma_data(1,:,mm,1,3)=0.5*( acoef_data(1,:,mm,1,2,1)-acoef_data(1,:,mm,1,1,2) &
                                         +(bcoef_data(1,:,mm,1,1,1,2)+bcoef_data(1,:,mm,1,2,2,2))/x(l1:l2))

              delta_data(1,:,mm,1,1)=0.25*(bcoef_data(1,:,mm,1,2,1,2)-bcoef_data(1,:,mm,1,2,2,1)-2.*bcoef_data(1,:,mm,1,3,3,1))
                            
              delta_data(1,:,mm,1,2)=0.25*(bcoef_data(1,:,mm,1,1,2,1)-bcoef_data(1,:,mm,1,1,1,2)-2.*bcoef_data(1,:,mm,1,3,3,2))
                            
              delta_data(1,:,mm,1,3)=0.5*(bcoef_data(1,:,mm,1,1,3,1)+bcoef_data(1,:,mm,1,2,3,2))

              beta_data(1,:,mm,1,1,1)= bcoef_data(1,:,mm,1,1,3,2)
              beta_data(1,:,mm,1,2,2)=-bcoef_data(1,:,mm,1,2,3,1)
              beta_data(1,:,mm,1,3,3)=0.5*(bcoef_data(1,:,mm,1,3,2,1)-bcoef_data(1,:,mm,1,3,1,2))

              beta_data(1,:,mm,1,1,2)=0.5*(bcoef_data(1,:,mm,1,2,3,2)-bcoef_data(1,:,mm,1,1,3,1))
              beta_data(1,:,mm,1,1,3)=0.25*( 2.*bcoef_data(1,:,mm,1,3,3,2)-bcoef_data(1,:,mm,1,1,1,2)+bcoef_data(1,:,mm,1,1,2,1))
              beta_data(1,:,mm,1,2,3)=0.25*(-2.*bcoef_data(1,:,mm,1,3,3,1)-bcoef_data(1,:,mm,1,2,1,2)+bcoef_data(1,:,mm,1,2,2,1))
 
              do ik=1,3
                kappa_data(1,:,mm,1,ik,:,3)=0.; kappa_data(1,:,mm,1,ik,3,:)=0.
              enddo
            enddo
            
            alpha_data(:,:,:,:,2,1)=alpha_data(:,:,:,:,1,2)
            alpha_data(:,:,:,:,3,1)=alpha_data(:,:,:,:,1,3)
            alpha_data(:,:,:,:,3,2)=alpha_data(:,:,:,:,2,3)

            beta_data(:,:,:,:,2,1)=beta_data(:,:,:,:,1,2)
            beta_data(:,:,:,:,3,1)=beta_data(:,:,:,:,1,3)
            beta_data(:,:,:,:,3,2)=beta_data(:,:,:,:,2,3)

            alpha_data=alpha_data*alpha_scale
            gamma_data=gamma_data*gamma_scale
            delta_data=delta_data*delta_scale
            beta_data =beta_data*beta_scale

          endif

          if (nsmooth_rbound/=0) then
            do i=1,3
              if (lgamma) call smooth_rbound(gamma_data(:,:,:,:,i),nsmooth_rbound)
              if (ldelta) call smooth_rbound(delta_data(:,:,:,:,i),nsmooth_rbound)
              if (lutensor .or. lmeanuu) call smooth_rbound(utensor_data(:,:,:,:,i),nsmooth_rbound)
              do j=1,3
                if (lacoef) call smooth_rbound(acoef_data(:,:,:,:,i,j),nsmooth_rbound)
                if (lalpha) call smooth_rbound(alpha_data(:,:,:,:,i,j),nsmooth_rbound)
                if (lbeta) call smooth_rbound(beta_data(:,:,:,:,i,j),nsmooth_rbound)
                do k=1,3
                  if (lbcoef) call smooth_rbound(bcoef_data(:,:,:,:,i,j,k),nsmooth_rbound)
                  if (lkappa) call smooth_rbound(kappa_data(:,:,:,:,i,j,k),nsmooth_rbound)
                enddo
              enddo
            enddo
          endif

          if (nsmooth_thbound/=0) then
            do i=1,3
              if (lgamma) call smooth_thbound(gamma_data(:,:,:,:,i),nsmooth_thbound)
              if (ldelta) call smooth_thbound(delta_data(:,:,:,:,i),nsmooth_thbound)
              if (lutensor .or. lmeanuu) call smooth_thbound(utensor_data(:,:,:,:,i),nsmooth_thbound)
              do j=1,3
                if (lacoef) call smooth_thbound(acoef_data(:,:,:,:,i,j),nsmooth_thbound)
                if (lalpha) call smooth_thbound(alpha_data(:,:,:,:,i,j),nsmooth_thbound)
                if (lbeta) call smooth_thbound(beta_data(:,:,:,:,i,j),nsmooth_thbound)
                do k=1,3
                  if (lbcoef) call smooth_thbound(bcoef_data(:,:,:,:,i,j,k),nsmooth_thbound)
                  if (lkappa) call smooth_thbound(kappa_data(:,:,:,:,i,j,k),nsmooth_thbound)
                enddo
              enddo
            enddo
          endif



!if (lroot.and.lbeta) write(200,*) beta_data(1,:,:,1,:,:)

           if (lbeta.and.lremove_beta_negativ) then
           if (iload==1) call get_shared_variable('eta', eta)
             do i=1,3
               where(beta_data(:,:,:,:,i,i)<eta*rel_eta) beta_data(:,:,:,:,i,i)=eta*rel_eta
             enddo
           endif



!if (lroot.and.lbeta) write(300,*) beta_data(1,:,:,1,:,:)

          if (lbeta.and.lregularize_beta) then

            allocate(beta_mask(dataload_len,nx,ny,nz))
            if (iload==1) call get_shared_variable('eta', eta)
            do i=1,3
!
!  Look for vanishing diagonal elements of beta.
!
              beta_mask=0
              where(beta_data(:,:,:,:,i,i)+eta<=0.) beta_mask=1
              minbeta=minval(beta_data(:,:,:,:,i,i),beta_mask==1)
              !where(beta_mask==1) beta_data(:,:,:,:,i,i)=(-1.+rel_eta)*eta 

              numzeros=sum(beta_mask)
              call mpiallreduce_sum_int(numzeros,numzeros_)
              if (numzeros_>0) then
                call mpireduce_min(minbeta,minbeta_)
                if (lroot) then
                  print'(a,i1,a,i1,a,i15,a$)', 'beta(',i,',',i,')+eta<=0 at ', numzeros_, ' positions'
                  print'(a,i1,a,i1,a,e12.5)', ', min(beta(',i,',',i,')) = ', minbeta_
                endif
              endif
            enddo
            deallocate(beta_mask)
!
!  Look for locations of negative definite beta.
!
            numzeros=0; numcomplex=0; rmin=x(l2); rmax=x(l1); thmin=y(m2); thmax=y(m1); trmin=impossible
            do mm=1,ny; do ll=1,nx
              beta_sav=beta_data(1,ll,mm,1,:,:)
              polcoeffs=(/2.*beta_data(1,ll,mm,1,1,2)*beta_data(1,ll,mm,1,1,3)*beta_data(1,ll,mm,1,2,3) &
                           - beta_data(1,ll,mm,1,1,2)**2*beta_data(1,ll,mm,1,3,3) &
                           - beta_data(1,ll,mm,1,1,3)**2*beta_data(1,ll,mm,1,2,2) &
                           - beta_data(1,ll,mm,1,2,3)**2*beta_data(1,ll,mm,1,1,1) &
                           + beta_data(1,ll,mm,1,1,1)*beta_data(1,ll,mm,1,2,2)*beta_data(1,ll,mm,1,3,3), &
!
                             beta_data(1,ll,mm,1,1,2)**2+beta_data(1,ll,mm,1,1,3)**2+beta_data(1,ll,mm,1,2,3)**2 &
                           - beta_data(1,ll,mm,1,1,1)*beta_data(1,ll,mm,1,2,2) &
                           - beta_data(1,ll,mm,1,1,1)*beta_data(1,ll,mm,1,3,3) &
                           - beta_data(1,ll,mm,1,2,2)*beta_data(1,ll,mm,1,3,3), &
!
                             beta_data(1,ll,mm,1,1,1)+beta_data(1,ll,mm,1,2,2)+beta_data(1,ll,mm,1,3,3), &
                           -1. /)
              call cubicroots(polcoeffs, eigenvals)
              if (any(abs(imag(eigenvals))>0.e-9*abs(real(eigenvals)))) numcomplex=numcomplex+1

!              if (any(real(eigenvals)<-eta)) then
!              if (ll <= 3 .and. mm <= 84 .and. mm >= 82) print*, 'JOERN: eigenvalues-real:', real(eigenvals)
!              if (ll <= 3 .and. mm <= 84 .and. mm >= 82) print*, 'JOERN: eigenvalues-imag:', imag(eigenvals)
!
!  If there are negative eigenvalues of beta,
!
!write(100,*) iproc, ll,mm,real(eigenvals)  !, &
!sum(polcoeffs*(/1.d0,real(eigenvals(1)),real(eigenvals(1))**2,real(eigenvals(1))**3/)), &
!sum(polcoeffs*(/1.d0,real(eigenvals(2)),real(eigenvals(2))**2,real(eigenvals(2))**3/)), &
!sum(polcoeffs*(/1.d0,real(eigenvals(3)),real(eigenvals(3))**2,real(eigenvals(3))**3/))
                numzeros=numzeros+1
                trmin=min(trmin,sum(real(eigenvals)))
                rmin=min(rmin,x(ll+nghost)); rmax=max(rmax,x(ll+nghost))
                thmin=min(thmin,y(mm+nghost)); thmax=max(thmax,y(mm+nghost))
!
!  calculate eigenvectors (v1,v2,-1) for all three eigenvalues.
!
                ev(:,3)=-1.
                do iv=1,3
                  det = (beta_data(1,ll,mm,1,1,1)-real(eigenvals(iv)))*(beta_data(1,ll,mm,1,2,2)-real(eigenvals(1))) &
                        -beta_data(1,ll,mm,1,1,2)**2
                  ev(iv,1)= (beta_data(1,ll,mm,1,1,3)*(beta_data(1,ll,mm,1,2,2)-real(eigenvals(iv))) &
                            -beta_data(1,ll,mm,1,2,3)*beta_data(1,ll,mm,1,1,2))/det
                  ev(iv,2)=((beta_data(1,ll,mm,1,1,1)-real(eigenvals(iv)))*beta_data(1,ll,mm,1,2,3) &
                            -beta_data(1,ll,mm,1,1,2)*beta_data(1,ll,mm,1,1,3))/det
                  ev(iv,:)=ev(iv,:)/sqrt(sum(ev(iv,:)**2))    ! normalization
!
! Eigenvalues < -eta are set to (-1.+rel_eta)*eta with rel_eta>0.
!
! output here: eigenvalues JOERN
!
!                  if (real(eigenvals(iv))<-eta) eigenvals(iv)=cmplx((-1.+rel_eta)*eta,0.)
                enddo

                beta_data(1,ll,mm,1,:,:)=0.
!
!  Transform back with non-negative eigenvalues.
!
                oldtrace=0.
                do iv=1,3
                  beta_data(1,ll,mm,1,:,:)=beta_data(1,ll,mm,1,:,:)+real(eigenvals(iv))*dyadic2_other(ev(iv,:))
                  oldtrace=oldtrace+real(eigenvals(iv))
                enddo
                newtrace=beta_data(1,ll,mm,1,1,1)+beta_data(1,ll,mm,1,2,2)+beta_data(1,ll,mm,1,3,3)
                if (abs(oldtrace-newtrace)>1.e-9) print*, 'll,mm,oldtrace-newtrace=', ll,mm,abs(oldtrace-newtrace)
!              endif
!            if (maxval(abs((beta_sav-beta_data(1,ll,mm,1,:,:))/beta_sav)) > 1e-6) print*, 'JOERN', ll, mm, maxval(abs((beta_sav-beta_data(1,ll,mm,1,:,:))/beta_sav))
            enddo; enddo
            call mpiallreduce_sum_int(numzeros,numzeros_)
            if (numzeros_>0) then
              if (lroot) print'(a,i15,a)', 'beta is negative definite at', numzeros_,' positions.'
              print'(4(a,f6.3),a,i4)', &
                   'in r-theta region (',rmin,',',rmax,')x(',thmin,',',thmax,') of proc ',iproc,'.'
              call mpireduce_min(trmin,minbeta_)
              if (lroot) print'(a,e12.5)', 'minimal trace = ', minbeta_
            endif
            call mpireduce_sum_int(numcomplex,numzeros_)
            if (lroot.and.numzeros_>0) print'(a,i15,a)', 'beta has complex eigenvalues at', numzeros_,' positions.'
          endif

!if (lroot.and.lbeta) write(100,*) beta_data(1,:,:,1,:,:)

          if (lsymmetrize) then
            do i=1,3 
              if (lgamma) call symmetrize(gamma_data(:,:,:,:,i),lgamma_sym(i))
              if (ldelta) call symmetrize(delta_data(:,:,:,:,i),ldelta_sym(i))
              if (lutensor.or. lmeanuu) call symmetrize(utensor_data(:,:,:,:,i),lutensor_sym(i))
              do j=1,3
                if (lalpha) call symmetrize(alpha_data(:,:,:,:,i,j),lalpha_sym(i,j))
                if (lbeta) call symmetrize(beta_data(:,:,:,:,i,j),lbeta_sym(i,j))
                do k=1,3
                  if (lkappa) call symmetrize(kappa_data(:,:,:,:,i,j,k),lkappa_sym(i,j,k))
                enddo
              enddo
            enddo
          else
            field_symmetry=0
          endif

        end if
      end if

      if (field_symmetry==1) then
        call symmetrize_3d(f(:,:,:,iax),.false.)
        call symmetrize_3d(f(:,:,:,iay),.true.)
        call symmetrize_3d(f(:,:,:,iaz),.false.)
      elseif (field_symmetry==-1) then
        call symmetrize_3d(f(:,:,:,iax),.true.)
        call symmetrize_3d(f(:,:,:,iay),.false.)
        call symmetrize_3d(f(:,:,:,iaz),.true.)
      endif

    endsubroutine special_before_boundary
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
      lpenc_requested(i_bb)=.true.
      !lpenc_requested(i_bij)=.true.
      lpenc_requested(i_bijtilde)=.true.
      lpenc_requested(i_jj)=.true.
      !if (lbcoef.and.lusecoefs) lpenc_requested(i_bijtilde)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
      use General, only: notanumber

      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      integer :: i,j,k, ind(1)
!
      call keep_compiler_quiet(f)
!
! Calculate emf pencil
!
      p%emf = 0
!
      if (lacoef) then
        ! Calculate acoef B
        do j=1,3; do i=1,3
          if (lacoef_arr(i,j)) then
            p%acoef_coefs(:,i,j)=emf_interpolate(acoef_data(1:dataload_len,:,m-nghost,n-nghost,i,j))
          else
            p%acoef_coefs(:,i,j)=0
          end if
        end do; end do
        call dot_mn_vm(p%bb,p%acoef_coefs,p%acoef_emf)
      end if
!
      if (lbcoef) then
        ! Calculate bcoef (grad B)
        do k=1,3; do j=1,3; do i=1,3
          if (lbcoef_arr(i,j,k)) then
            p%bcoef_coefs(:,i,j,k)=emf_interpolate(bcoef_data(1:dataload_len,:,m-nghost,n-nghost,i,j,k))
          else
            p%bcoef_coefs(:,i,j,k)=0
          end if
        end do; end do; end do
        p%bcoef_emf = 0
! 
!  Use partial (non-covariant) derivatives of B in the form \partial B_{r,theta,phi}/\partial r, 
!  \partial B_{r,theta,phi}/(r \partial theta).
!
        do k=1,2; do j=1,3; do i=1,3
          if (lbcoef_arr(i,j,k)) &
            p%bcoef_emf(:,i)=p%bcoef_emf(:,i)+p%bcoef_coefs(:,i,j,k)*p%bijtilde(:,j,k)
        end do; end do; end do
      end if
!
      if (lalpha) then
        ! Calculate alpha B
        do j=1,3; do i=1,3
          if (lalpha_arr(i,j)) then
            p%alpha_coefs(:,i,j)=emf_interpolate(alpha_data(1:dataload_len,:,m-nghost,n-nghost,i,j))
          else
            p%alpha_coefs(:,i,j)=0
          end if
        end do; end do

        call dot_mn_vm(p%bb,p%alpha_coefs,p%alpha_emf)
        p%emf = p%emf + p%alpha_emf
      end if
!
      if (lbeta) then
        ! Calculate beta (curl B)
        do j=1,3; do i=1,3
          if (lbeta_arr(i,j)) then
            p%beta_coefs(:,i,j)=emf_interpolate(beta_data(1:dataload_len,:,m-nghost,n-nghost,i,j))
          else
            p%beta_coefs(:,i,j)=0
          end if
        end do; end do
        call dot_mn_vm(p%jj,p%beta_coefs,p%beta_emf)
        p%emf = p%emf - p%beta_emf
      end if
!
      if (lgamma) then
        ! Calculate gamma x B
        do i=1,3
          if (lgamma_arr(i)) then
            p%gamma_coefs(:,i)=emf_interpolate(gamma_data(1:dataload_len,:,m-nghost,n-nghost,i))
          else
            p%gamma_coefs(:,i)=0
          end if
        end do

        call cross_mn(p%gamma_coefs,p%bb,p%gamma_emf)
        p%emf = p%emf + p%gamma_emf
      end if
!
      if (ldelta) then
        ! Calculate delta x (curl B)
        do i=1,3
          if (ldelta_arr(i)) then
            p%delta_coefs(:,i)=emf_interpolate(delta_data(1:dataload_len,:,m-nghost,n-nghost,i))
          else
            p%delta_coefs(:,i)=0
          end if
        end do
        call cross_mn(p%delta_coefs,p%jj,p%delta_emf)
        p%emf = p%emf - p%delta_emf
      end if

      if (lkappa) then
        ! Calculate kappa (grad B)_symm
        do j=1,3; do i=1,3
          p%bij_symm(:,i,j)=0.5*(p%bijtilde(:,i,j)+p%bij_cov_corr(:,i,j) + &
                                 p%bijtilde(:,j,i)+p%bij_cov_corr(:,j,i))
        end do; end do

        do k=1,3; do j=1,3; do i=1,3
          if (lkappa_arr(i,j,k)) then
            p%kappa_coefs(:,i,j,k)=emf_interpolate(kappa_data(1:dataload_len,:,m-nghost,n-nghost,i,j,k))
          else
            p%kappa_coefs(:,i,j,k)=0
          endif
        end do; end do; end do

        p%kappa_emf = 0
        do k=1,3; do j=1,3; do i=1,3
          if (lkappa_arr(i,j,k)) &
            p%kappa_emf(:,i)=p%kappa_emf(:,i)+p%kappa_coefs(:,i,j,k)*p%bij_symm(:,j,k)
        end do; end do; end do
        p%emf = p%emf - p%kappa_emf
      end if
!
      if (lutensor .or. lmeanuu) then
        ! Calculate U x B
        do i=1,3
          if (lutensor_arr(i)) then
            p%utensor_coefs(:,i)=emf_interpolate(utensor_data(1:dataload_len,:,m-nghost,n-nghost,i))
          else
            p%utensor_coefs(:,i)=0
          end if
        end do
        call cross_mn(p%utensor_coefs,p%bb,p%utensor_emf)
        p%emf = p%emf + p%utensor_emf
      end if
!
if (.false.) then
if (maxval(abs(p%bcoef_coefs(:,2,2,1)-p%bcoef_coefs(:,2,1,2)-2.*(p%beta_coefs(:,2,3)-p%delta_coefs(:,1))))>1.e-9) print*, '2,3,m,n=', m,n
if (maxval(abs(p%bijtilde(:,2,1)-p%bijtilde(:,1,2)+p%bb(:,2)/x(l1:l2)-p%jj(:,3)))>1.e-10) print*, 'jphi,m=', m, &
maxval(abs(p%jj(:,3))), maxval(abs(p%bijtilde(:,2,1)-p%bijtilde(:,1,2)+p%bb(:,2)/x(l1:l2)))
endif
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  Several precalculated Pencils of information are passed for
!  efficiency.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df

      integer :: i, j
!
!  Identify module and boundary conditions.
!
      if (lroot.and.(headtt.or.ldebug)) print*,'dspecial_dt: SOLVE dspecial_dt'
!!      if (headtt) call identify_bcs('special',ispecial)
!
!
!! SAMPLE DIAGNOSTIC IMPLEMENTATION
!!
!!      if (ldiagnos) then
!!        if (idiag_SPECIAL_DIAGNOSTIC/=0) then
!!          call sum_mn_name(MATHEMATICAL EXPRESSION,idiag_SPECIAL_DIAGNOSTIC)
!!! see also integrate_mn_name
!!        endif
!!      endif
!      emftmp=0
!
      if (ldiagnos) then

        tmppencil = emftmp - p%emf
        !
        if (idiag_alphaxmax/=0) call max_mn_name(p%alpha_emf(:,1),idiag_alphaxmax)
        if (idiag_alphaymax/=0) call max_mn_name(p%alpha_emf(:,2),idiag_alphaymax)
        if (idiag_alphazmax/=0) call max_mn_name(p%alpha_emf(:,3),idiag_alphazmax)
        if (idiag_betaxmax/=0) call max_mn_name(p%beta_emf(:,1),idiag_betaxmax)
        if (idiag_betaymax/=0) call max_mn_name(p%beta_emf(:,2),idiag_betaymax)
        if (idiag_betazmax/=0) call max_mn_name(p%beta_emf(:,3),idiag_betazmax)
        if (idiag_gammaxmax/=0) call max_mn_name(p%gamma_emf(:,1),idiag_gammaxmax)
        if (idiag_gammaymax/=0) call max_mn_name(p%gamma_emf(:,2),idiag_gammaymax)
        if (idiag_gammazmax/=0) call max_mn_name(p%gamma_emf(:,3),idiag_gammazmax)
        if (idiag_deltaxmax/=0) call max_mn_name(p%delta_emf(:,1),idiag_deltaxmax)
        if (idiag_deltaymax/=0) call max_mn_name(p%delta_emf(:,2),idiag_deltaymax)
        if (idiag_deltazmax/=0) call max_mn_name(p%delta_emf(:,3),idiag_deltazmax)
        if (idiag_kappaxmax/=0) call max_mn_name(p%kappa_emf(:,1),idiag_kappaxmax)
        if (idiag_kappaymax/=0) call max_mn_name(p%kappa_emf(:,2),idiag_kappaymax)
        if (idiag_kappazmax/=0) call max_mn_name(p%kappa_emf(:,3),idiag_kappazmax)
        if (idiag_utensorxmax/=0) call max_mn_name(p%utensor_emf(:,1),idiag_utensorxmax)
        if (idiag_utensorymax/=0) call max_mn_name(p%utensor_emf(:,2),idiag_utensorymax)
        if (idiag_utensorzmax/=0) call max_mn_name(p%utensor_emf(:,3),idiag_utensorzmax)
        if (idiag_acoefxmax/=0) call max_mn_name(p%acoef_emf(:,1),idiag_acoefxmax)
        if (idiag_acoefymax/=0) call max_mn_name(p%acoef_emf(:,2),idiag_acoefymax)
        if (idiag_acoefzmax/=0) call max_mn_name(p%acoef_emf(:,3),idiag_acoefzmax)
        if (idiag_bcoefxmax/=0) call max_mn_name(p%bcoef_emf(:,1),idiag_bcoefxmax)
        if (idiag_bcoefymax/=0) call max_mn_name(p%bcoef_emf(:,2),idiag_bcoefymax)
        if (idiag_bcoefzmax/=0) call max_mn_name(p%bcoef_emf(:,3),idiag_bcoefzmax)
        if (idiag_emfxmax/=0) call max_mn_name(p%emf(:,1),idiag_emfxmax)
        if (idiag_emfymax/=0) call max_mn_name(p%emf(:,2),idiag_emfymax)
        if (idiag_emfzmax/=0) call max_mn_name(p%emf(:,3),idiag_emfzmax)
        if (idiag_emfcoefxmax/=0) call max_mn_name(emftmp(:,1),idiag_emfcoefxmax)
        if (idiag_emfcoefymax/=0) call max_mn_name(emftmp(:,2),idiag_emfcoefymax)
        if (idiag_emfcoefzmax/=0) call max_mn_name(emftmp(:,3),idiag_emfcoefzmax)
        if (idiag_emfdiffmax/=0) then
          call dot2_mn(tmppencil,tmpline)
          call max_mn_name(tmpline,idiag_emfdiffmax,lsqrt=.true.)
        end if
        if (idiag_emfxdiffmax/=0) &
          call max_mn_name(abs(tmppencil(:,1)),idiag_emfxdiffmax)
        if (idiag_emfydiffmax/=0) &
          call max_mn_name(abs(tmppencil(:,2)),idiag_emfydiffmax)
        if (idiag_emfzdiffmax/=0) &
          call max_mn_name(abs(tmppencil(:,3)),idiag_emfzdiffmax)
        if (idiag_alpharms/=0) then
          call dot2_mn(p%alpha_emf,tmpline)
          call sum_mn_name(tmpline,idiag_alpharms,lsqrt=.true.)
        end if
        if (idiag_betarms/=0) then
          call dot2_mn(p%beta_emf,tmpline)
          call sum_mn_name(tmpline,idiag_betarms,lsqrt=.true.)
        end if
        if (idiag_gammarms/=0) then
          call dot2_mn(p%gamma_emf,tmpline)
          call sum_mn_name(tmpline,idiag_gammarms,lsqrt=.true.)
        end if
        if (idiag_deltarms/=0) then
          call dot2_mn(p%delta_emf,tmpline)
          call sum_mn_name(tmpline,idiag_deltarms,lsqrt=.true.)
        end if
        if (idiag_kapparms/=0) then
          call dot2_mn(p%kappa_emf,tmpline)
          call sum_mn_name(tmpline,idiag_kapparms,lsqrt=.true.)
        end if
        if (idiag_utensorrms/=0) then
          call dot2_mn(p%utensor_emf,tmpline)
          call sum_mn_name(tmpline,idiag_utensorrms,lsqrt=.true.)
        end if
        if (idiag_acoefrms/=0) then
          call dot2_mn(p%acoef_emf,tmpline)
          call sum_mn_name(tmpline,idiag_acoefrms,lsqrt=.true.)
        end if
        if (idiag_bcoefrms/=0) then
          call dot2_mn(p%bcoef_emf,tmpline)
          call sum_mn_name(tmpline,idiag_bcoefrms,lsqrt=.true.)
        end if
        if (idiag_emfrms/=0) then
          call dot2_mn(p%emf,tmpline)
          call sum_mn_name(tmpline,idiag_emfrms,lsqrt=.true.)
        end if
        if (idiag_emfcoefrms/=0) then
          call dot2_mn(emftmp,tmpline)
          call sum_mn_name(tmpline,idiag_emfcoefrms,lsqrt=.true.)
        end if
        if (idiag_emfdiffrms/=0) then
          call dot2_mn(tmppencil,tmpline)
          call sum_mn_name(tmpline,idiag_emfdiffrms,lsqrt=.true.)
        end if
      end if 
!
      if (l2davgfirst) then
        call zsum_mn_name_xy(p%emf(:,1),idiag_emfxmxy)
        call zsum_mn_name_xy(p%emf(:,2),idiag_emfymxy)
        call zsum_mn_name_xy(p%emf(:,3),idiag_emfzmxy)
        call zsum_mn_name_xy(emftmp(:,1),idiag_emfcoefxmxy)
        call zsum_mn_name_xy(emftmp(:,2),idiag_emfcoefymxy)
        call zsum_mn_name_xy(emftmp(:,3),idiag_emfcoefzmxy)
      endif

      call keep_compiler_quiet(f,df)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
      call setParameterDefaults
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
      if (lroot) write (*,*) 'read_special_init_pars parameters read...'
      call parseParameters
      if (lroot) write (*,*) 'read_special_init_pars parameters parsed...'
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)

   endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
      
      call setParameterDefaults
      if (lroot) write (*,*) 'read_special_run_pars parameters read...'
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
      call parseParameters
      if (lroot) write (*,*) 'read_special_run_pars parameters parsed...'
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
!!      use FArrayManager, only: farray_index_append
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        ! Max diagnostics
        idiag_alphaxmax=0
        idiag_alphaymax=0
        idiag_alphazmax=0
        idiag_betaxmax=0
        idiag_betaymax=0
        idiag_betazmax=0
        idiag_gammaxmax=0
        idiag_gammaymax=0
        idiag_gammazmax=0
        idiag_deltaxmax=0
        idiag_deltaymax=0
        idiag_deltazmax=0
        idiag_kappaxmax=0
        idiag_kappaymax=0
        idiag_kappazmax=0
        idiag_utensorxmax=0
        idiag_utensorymax=0
        idiag_utensorzmax=0
        idiag_acoefxmax=0
        idiag_acoefymax=0
        idiag_acoefzmax=0
        idiag_bcoefxmax=0
        idiag_bcoefymax=0
        idiag_bcoefzmax=0
        idiag_emfxmax=0
        idiag_emfymax=0
        idiag_emfzmax=0
        idiag_emfcoefxmax=0
        idiag_emfcoefymax=0
        idiag_emfcoefzmax=0
        idiag_emfdiffmax=0
        idiag_emfxdiffmax=0
        idiag_emfydiffmax=0
        idiag_emfzdiffmax=0
!     RMS diagnostics
        idiag_alpharms=0
        idiag_betarms=0
        idiag_gammarms=0
        idiag_deltarms=0
        idiag_kapparms=0
        idiag_utensorrms=0
        idiag_acoefrms=0
        idiag_bcoefrms=0
        idiag_emfrms=0
        idiag_emfcoefrms=0
        idiag_emfdiffrms=0
!    timestep diagnostics
        idiag_dtemf_ave=0
        idiag_dtemf_dif=0
!    2D diagnostics
        idiag_emfxmxy=0; idiag_emfymxy=0; idiag_emfzmxy=0
        idiag_emfcoefxmxy=0; idiag_emfcoefymxy=0; idiag_emfcoefzmxy=0
      endif
!
      do iname=1,nname
        ! Maximum values of emf terms
        call parse_name(iname,cname(iname),cform(iname),'alphaxmax',idiag_alphaxmax)
        call parse_name(iname,cname(iname),cform(iname),'alphaymax',idiag_alphaymax)
        call parse_name(iname,cname(iname),cform(iname),'alphazmax',idiag_alphazmax)
        call parse_name(iname,cname(iname),cform(iname),'betaxmax',idiag_betaxmax)
        call parse_name(iname,cname(iname),cform(iname),'betaymax',idiag_betaymax)
        call parse_name(iname,cname(iname),cform(iname),'betazmax',idiag_betazmax)
        call parse_name(iname,cname(iname),cform(iname),'gammaxmax',idiag_gammaxmax)
        call parse_name(iname,cname(iname),cform(iname),'gammaymax',idiag_gammaymax)
        call parse_name(iname,cname(iname),cform(iname),'gammazmax',idiag_gammazmax)
        call parse_name(iname,cname(iname),cform(iname),'deltaxmax',idiag_deltaxmax)
        call parse_name(iname,cname(iname),cform(iname),'deltaymax',idiag_deltaymax)
        call parse_name(iname,cname(iname),cform(iname),'deltazmax',idiag_deltazmax)
        call parse_name(iname,cname(iname),cform(iname),'kappaxmax',idiag_kappaxmax)
        call parse_name(iname,cname(iname),cform(iname),'kappaymax',idiag_kappaymax)
        call parse_name(iname,cname(iname),cform(iname),'kappazmax',idiag_kappazmax)
        call parse_name(iname,cname(iname),cform(iname),'utensorxmax',idiag_utensorxmax)
        call parse_name(iname,cname(iname),cform(iname),'utensorymax',idiag_utensorymax)
        call parse_name(iname,cname(iname),cform(iname),'utensorzmax',idiag_utensorzmax)
        call parse_name(iname,cname(iname),cform(iname),'acoefxmax',idiag_acoefxmax)
        call parse_name(iname,cname(iname),cform(iname),'acoefymax',idiag_acoefymax)
        call parse_name(iname,cname(iname),cform(iname),'acoefzmax',idiag_acoefzmax)
        call parse_name(iname,cname(iname),cform(iname),'bcoefxmax',idiag_bcoefxmax)
        call parse_name(iname,cname(iname),cform(iname),'bcoefymax',idiag_bcoefymax)
        call parse_name(iname,cname(iname),cform(iname),'bcoefzmax',idiag_bcoefzmax)
        call parse_name(iname,cname(iname),cform(iname),'emfxmax',idiag_emfxmax)
        call parse_name(iname,cname(iname),cform(iname),'emfymax',idiag_emfymax)
        call parse_name(iname,cname(iname),cform(iname),'emfzmax',idiag_emfzmax)
        call parse_name(iname,cname(iname),cform(iname),'emfcoefxmax',idiag_emfcoefxmax)
        call parse_name(iname,cname(iname),cform(iname),'emfcoefymax',idiag_emfcoefymax)
        call parse_name(iname,cname(iname),cform(iname),'emfcoefzmax',idiag_emfcoefzmax)
        call parse_name(iname,cname(iname),cform(iname),'emfdiffmax',idiag_emfdiffmax)
        call parse_name(iname,cname(iname),cform(iname),'emfxdiffmax',idiag_emfxdiffmax)
        call parse_name(iname,cname(iname),cform(iname),'emfydiffmax',idiag_emfydiffmax)
        call parse_name(iname,cname(iname),cform(iname),'emfzdiffmax',idiag_emfzdiffmax)
        ! RMS values of emf terms
        call parse_name(iname,cname(iname),cform(iname),'alpharms',idiag_alpharms)
        call parse_name(iname,cname(iname),cform(iname),'betarms',idiag_betarms)
        call parse_name(iname,cname(iname),cform(iname),'gammarms',idiag_gammarms)
        call parse_name(iname,cname(iname),cform(iname),'deltarms',idiag_deltarms)
        call parse_name(iname,cname(iname),cform(iname),'kapparms',idiag_kapparms)
        call parse_name(iname,cname(iname),cform(iname),'utensorrms',idiag_utensorrms)
        call parse_name(iname,cname(iname),cform(iname),'acoefrms',idiag_acoefrms)
        call parse_name(iname,cname(iname),cform(iname),'bcoefrms',idiag_bcoefrms)
        call parse_name(iname,cname(iname),cform(iname),'emfrms',idiag_emfrms)
        call parse_name(iname,cname(iname),cform(iname),'emfcoefrms',idiag_emfcoefrms)
        call parse_name(iname,cname(iname),cform(iname),'emfdiffrms',idiag_emfdiffrms)
!    timestep diagnostics
        call parse_name(iname,cname(iname),cform(iname),'dtemf_ave',idiag_dtemf_ave)
        call parse_name(iname,cname(iname),cform(iname),'dtemf_dif',idiag_dtemf_dif)
      enddo

      do iname=1,nnamexy
        call parse_name(iname,cnamexy(iname),cformxy(iname),'EMFxmxy',idiag_emfxmxy)
        call parse_name(iname,cnamexy(iname),cformxy(iname),'EMFymxy',idiag_emfymxy)
        call parse_name(iname,cnamexy(iname),cformxy(iname),'EMFzmxy',idiag_emfzmxy)
        call parse_name(iname,cnamexy(iname),cformxy(iname),'EMFcoefxmxy',idiag_emfcoefxmxy)
        call parse_name(iname,cnamexy(iname),cformxy(iname),'EMFcoefymxy',idiag_emfcoefymxy)
        call parse_name(iname,cnamexy(iname),cformxy(iname),'EMFcoefzmxy',idiag_emfcoefzmxy)
      enddo
 
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  induction equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      real :: diffus_tmp
      type (pencil_case), intent(in) :: p
      integer :: i,j,k
!
      call keep_compiler_quiet(f)
!
! Overwrite with a and b coefs if needed
!
      if (lusecoefs) then
        emftmp=0
        if (lacoef) emftmp = emftmp + p%acoef_emf
        if (lbcoef) emftmp = emftmp - p%bcoef_emf
        if (lutensor .or. lmeanuu) emftmp = emftmp + p%utensor_emf
      else
        emftmp = p%emf
      end if

      df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+emftmp
!
      if (lfirst.and.ldt) then

        advec_special=0.0
        diffus_special=0.0
!
! Calculate advec_special
!
        advec_special=0.
        if (lalpha) then
          call dot_mn_vm(dline_1, abs(p%alpha_coefs), tmppencil)
          advec_special=advec_special+sum(tmppencil,2)
        end if
        if (lgamma) then
          call dot_mn(dline_1, abs(p%gamma_coefs), tmpline)
          advec_special=advec_special+tmpline
        end if
        if (lutensor .or. lmeanuu) then
          call dot_mn(dline_1, abs(p%utensor_coefs), tmpline)
          advec_special=advec_special+tmpline
        end if
!
        maxadvec=maxadvec+advec_special
!
        if (ldiagnos.and.idiag_dtemf_ave/=0) then
          call max_mn_name(advec_special/cdt,idiag_dtemf_ave,l_dt=.true.)
        endif
!
! Calculate diffus_special
!
        diffus_special=0.
        if (lbeta) then
          call dot_mn_vm(dline_1,abs(p%beta_coefs), tmppencil)
          call dot_mn(dline_1, tmppencil, tmpline)
          diffus_special=diffus_special+tmpline
        end if
!        
        if (ldelta) then
          call cross_mn(dline_1,abs(p%delta_coefs), tmppencil)
          call dot_mn(dline_1,tmppencil,tmpline)
          diffus_special=diffus_special+tmpline
        end if
!        
        if (lkappa) then
          call vec_dot_3tensor(dline_1, abs(p%kappa_coefs), tmptensor)
          call dot_mn_vm(dline_1,tmptensor, tmppencil)
          diffus_special=diffus_special+sum(tmppencil,2)
        end if

        maxdiffus=max(maxdiffus,diffus_special)
!
        if (ldiagnos.and.idiag_dtemf_dif/=0) then
          call max_mn_name(diffus_special/cdtv,idiag_dtemf_dif,l_dt=.true.)
        endif
!
      end if 
!
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine openDataset(datagroup,tensor_id)

      ! Open a dataset e.g. /emftensor/alpha/data and auxillary dataspaces

      character(len=*), intent(in)    :: datagroup     ! name of data group
      integer, intent(in)             :: tensor_id
!
      integer(HSIZE_T), dimension(10) :: dimsizes, maxdimsizes
      integer :: ndims
      integer(HSIZE_T) :: num
      logical :: hdf_exists
      character(len=fnlen)            :: dataset

      dataset = tensor_names(tensor_id)               ! name of dataset.
      ! Check that datagroup e.g. /emftensor/alpha exists
      call H5Lexists_F(hdf_emftensors_group, datagroup, hdf_exists, hdferr)
      if (.not. hdf_exists) then
        call fatal_error('openDataset','/emftensor/'//trim(datagroup)// &
                          ' does not exist')
      end if
      ! Open datagroup, returns identifier in tensor_id_G. 
      call H5Gopen_F(hdf_emftensors_group, datagroup, tensor_id_G(tensor_id),hdferr)
      if (hdferr /= 0) then
        call fatal_error('openDataset','Error opening /emftensor/'//trim(datagroup))
      end if
      ! Check that dataset e.g. /emftensor/alpha/mean exists
      call H5Lexists_F(tensor_id_G(tensor_id), dataset, hdf_exists, hdferr)
      if (.not. hdf_exists) then
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        call fatal_error('openDataset','/emftensor/'//trim(datagroup)// &
                          '/'//trim(dataset)//' does not exist')
      end if
      ! Open dataset <dataset> in group, returns identifier in tensor_id_D.
      call H5Dopen_F(tensor_id_G(tensor_id), dataset, tensor_id_D(tensor_id),hdferr)
      if (hdferr /= 0) then
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        call fatal_error('openDataset','Error opening /emftensor/'// &
                          trim(datagroup)//'/'//trim(dataset))
      end if
      ! Get dataspace of dataset (identifier of copy of dataspace in tensor_id_S).
      call H5Dget_space_F(tensor_id_D(tensor_id), tensor_id_S(tensor_id),hdferr)
      if (hdferr /= 0) then
        call H5Dclose_F(tensor_id_D(tensor_id), hdferr)
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        call fatal_error('openDataset','Error opening dataspace for/emftensor/'// &
                          trim(datagroup)//'/'//trim(dataset))
      end if
      ! Get dataspace dimensions in tensor_dims.
      ndims = tensor_ndims(tensor_id)
      call H5Sget_simple_extent_dims_F(tensor_id_S(tensor_id), &
                                       dimsizes(1:ndims), &
                                       maxdimsizes(1:ndims), &
                                       hdferr)                 !MR: hdferr/=0!
!print*, 'from H5Sget_simple_extent_dims_F, line 1330: hdferr=', hdferr
      call H5Sget_simple_extent_npoints_F(tensor_id_S(tensor_id),num,hdferr) ! This is to mask the error of the preceding call.
      tensor_dims(tensor_id,1:ndims)=dimsizes(1:ndims)
      if (tensor_times_len==-1) then
        tensor_times_len=dimsizes(1)
      elseif (tensor_times_len/=dimsizes(1)) then
        call fatal_error('openDataset','dataset emftensor/'//trim(datagroup)//'/'//trim(dataset)//' has deviating time extent')  
      endif
      if (hdferr /= 0) then
        call H5Sclose_F(tensor_id_S(tensor_id), hdferr)
        call H5Dclose_F(tensor_id_D(tensor_id), hdferr)
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        call fatal_error('openDataset','Cannot get dimensions extent '// &
                          'for /emftensor/'//trim(datagroup)//'/'//trim(dataset))
      end if

      ! Create a memory space mapping for input data (identifier in tensor_id_memS).
      call H5Screate_simple_F(ndims, tensor_memdims(tensor_id,1:ndims), &
                              tensor_id_memS(tensor_id), hdferr)
      if (hdferr /= 0) then
        call H5Sclose_F(tensor_id_S(tensor_id), hdferr)
        call H5Dclose_F(tensor_id_D(tensor_id), hdferr)
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        call fatal_error('openDataset','Error creating memory mapping '// &
                          'for /emftensor/'//trim(datagroup)//'/'//trim(dataset))
      end if

      if (lroot) write(*,*) 'Using dataset /emftensor/'//trim(datagroup)//'/'//trim(dataset)//'.'

    end subroutine openDataset
!***********************************************************************
    subroutine openDataset_grid(scalar_id)

      ! Open a dataset in group grid, e.g. /grid/t

      integer, intent(in)             :: scalar_id
!
      integer(HSIZE_T), dimension(10) :: maxdimsizes
      integer :: ndims
      logical :: hdf_exists
      character(len=fnlen)            :: dataset

      dataset = scalar_names(scalar_id)               ! name of dataset

      ! Check that dataset exists
      call H5Lexists_F(hdf_grid_group, dataset, hdf_exists, hdferr)
      if (.not. hdf_exists) then
        call H5Gclose_F(hdf_grid_group, hdferr)
        call fatal_error('openDataset','/grid/' &
                          //trim(dataset)//' does not exist')
      end if
      ! Open dataset dataset in group, returns identifier in scalar_id_D.
      call H5Dopen_F(hdf_grid_group, dataset, scalar_id_D(scalar_id),hdferr)
      if (hdferr /= 0) then
        call H5Gclose_F(hdf_grid_group, hdferr)
        call fatal_error('openDataset','Error opening /grid/' &
                          //trim(dataset))
      end if
      ! Get dataspace of dataset (identifier of copy of dataspace in scalar_id_S).
      call H5Dget_space_F(scalar_id_D(scalar_id), scalar_id_S(scalar_id),hdferr)
      if (hdferr /= 0) then
        call H5Dclose_F(scalar_id_D(scalar_id), hdferr)
        call H5Gclose_F(hdf_grid_group, hdferr)
        call fatal_error('openDataset','Error opening dataspace for/grid/' &
                          //trim(dataset))
      end if

      ! Get dataspace dimensions in scalar_dims.
      call H5Sget_simple_extent_npoints_F(scalar_id_S(scalar_id), &
                                          scalar_dims(scalar_id), &
                                          hdferr)
      if (hdferr /= 0) then
        call H5Sclose_F(scalar_id_S(scalar_id), hdferr)
        call H5Dclose_F(scalar_id_D(scalar_id), hdferr)
        call H5Gclose_F(hdf_grid_group, hdferr)
        call fatal_error('openDataset','Error in getting dimensions for grid/t')
      end if

    end subroutine openDataset_grid
!*********************************************************************** 
    subroutine closeDataset_grid(id)

      ! Close opened dataspaces, dataset and group

      integer :: id

      call H5Sclose_F(scalar_id_S(id), hdferr)
      call H5Dclose_F(scalar_id_D(id), hdferr)

    end subroutine closeDataset_grid
!*********************************************************************** 
    subroutine closeDataset(tensor_id)

      ! Close opened dataspaces, dataset and group

      integer :: tensor_id

      call H5Sclose_F(tensor_id_memS(tensor_id), hdferr)
      call H5Sclose_F(tensor_id_S(tensor_id), hdferr)
      call H5Dclose_F(tensor_id_D(tensor_id), hdferr)
      call H5Gclose_F(tensor_id_G(tensor_id), hdferr)

    end subroutine closeDataset
!*********************************************************************** 
    subroutine loadDataset_rank1(dataarray, datamask, tensor_id, loadstart,name)

      ! Load a chunk of data for a vector, beginning at loadstart
    
      use General, only: itoa

      real, dimension(:,:,:,:,:), intent(inout) :: dataarray
      logical, dimension(3), intent(in) :: datamask
      integer, intent(in) :: tensor_id
      integer, intent(in) :: loadstart
      character(*), optional :: name

      integer :: ndims,i
      integer(HSIZE_T) :: mask_i
      real :: globmin, globmax, sum, rms ! output for diagnostics

      ndims = tensor_ndims(tensor_id)
      tensor_offsets(tensor_id,1) = loadstart
      call H5Sselect_none_F(tensor_id_S(tensor_id), hdferr)     ! resets selection region.
      call H5Sselect_none_F(tensor_id_memS(tensor_id), hdferr)  ! perhaps dispensable when
                                                                ! H5S_SELECT_SET_F is used below.
      do mask_i=1,3
        ! Load only wanted datasets
        if (datamask(mask_i)) then
          ! Set the new offset for data reading
          tensor_offsets(tensor_id,ndims)    = mask_i-1
          tensor_memoffsets(tensor_id,ndims) = mask_i-1
          ! Select hyperslab for data.
!          print '(a,a,5(1x,i3))', 'tensor offset',name,tensor_offsets(tensor_id,1:ndims) 
!          print '(a,a,5(1x,i3))', 'tensor memoffset',name,tensor_memoffsets(tensor_id,1:ndims) 
!          print '(a,a,5(1x,i3))', 'tensor counts',name,tensor_counts(tensor_id,:ndims) 
!          print '(a,a,5(1x,i3))', 'tensor memcounts',name,tensor_memcounts(tensor_id,:ndims) 

          call H5Sselect_hyperslab_F(tensor_id_S(tensor_id), H5S_SELECT_OR_F, &
                                     tensor_offsets(tensor_id,1:ndims),       &
                                     tensor_counts(tensor_id,1:ndims),        &
                                     hdferr)
           if (hdferr /= 0) then
             call fatal_error('loadDataset_rank1','Error creating File mapping '// &
                           'for /grid/'//name)
           end if
          ! Select hyperslab for memory.
          call H5Sselect_hyperslab_F(tensor_id_memS(tensor_id), H5S_SELECT_OR_F, &
                                     tensor_memoffsets(tensor_id,1:ndims),       &
                                     tensor_memcounts(tensor_id,1:ndims),        &
                                     hdferr)
          if (hdferr /= 0) then
            call fatal_error('loadDataset_rank1','Error creating memory mapping '// &
                           'for /grid/'//name)
          end if
        end if
      end do
      ! Read data into memory.
      tensor_dims(tensor_id,ndims)=3
      call H5Dread_F(tensor_id_D(tensor_id), hdf_memtype, dataarray, &
                     tensor_dims(tensor_id,1:ndims), hdferr, &
                     tensor_id_memS(tensor_id), tensor_id_S(tensor_id))
      if (hdferr /= 0) &
        call fatal_error('loadDataset_rank1','Error creating reading dataset '// &
                         'for /grid/'//name)
      sum = 0.; rms = 0.
      do i=1,3
        tensor_maxvals(tensor_id) = maxval(dataarray(:,:,:,:,i))
        tensor_minvals(tensor_id) = minval(dataarray(:,:,:,:,i))
        if (present(name)) then
          call mpireduce_min(tensor_minvals(tensor_id),globmin)
          call mpireduce_max(tensor_maxvals(tensor_id),globmax)
          if (lroot) write (*,*) trim(name)//'(', trim(itoa(i)), ')  min/max:', globmin, globmax
        endif
      sum = sum + tensor_maxvals(tensor_id)*tensor_maxvals(tensor_id)
      enddo
      rms = sqrt(sum)
      if (present(name) .and. lroot) write (*,*) trim(name)//' (rms):', rms

      dataarray = tensor_scales(tensor_id) * dataarray

    end subroutine loadDataset_rank1
!*********************************************************************** 
    subroutine loadDataset_rank2(dataarray, datamask, tensor_id, loadstart,name)

      ! Load a chunk of data for a 2-rank tensor, beginning at loadstart

      use General, only: itoa

      real, dimension(:,:,:,:,:,:), intent(inout) :: dataarray
      logical, dimension(3,3), intent(in) :: datamask
      integer, intent(in) :: tensor_id
      integer, intent(in) :: loadstart
      character(*), optional :: name

      integer :: ndims, i, j
      integer(HSIZE_T) :: mask_i, mask_j
      real :: globmin, globmax, sum, rms ! output for diagnostics

      ndims = tensor_ndims(tensor_id)
      tensor_offsets(tensor_id,1) = loadstart
      call H5Sselect_none_F(tensor_id_S(tensor_id), hdferr)
      call H5Sselect_none_F(tensor_id_memS(tensor_id), hdferr)

      do mask_j=1,3; do mask_i=1,3

        ! Load only wanted datasets
        if (datamask(mask_i,mask_j)) then

          ! Set the new offset for data reading
          tensor_offsets(tensor_id,ndims-1)    = mask_i-1
          tensor_offsets(tensor_id,ndims)      = mask_j-1
          tensor_memoffsets(tensor_id,ndims-1) = mask_i-1
          tensor_memoffsets(tensor_id,ndims)   = mask_j-1

          ! Hyperslab for data
          !if (mask_i==1.and.mask_j==1) print '(a,i2,a,6(1x,i3))', 'iproc, tensor offset',iproc,name,tensor_offsets(tensor_id,1:ndims-2)
!          print '(a,i2,a,6(1x,i3))', 'iproc, tensor memoffset',iproc,name,tensor_memoffsets(tensor_id,1:ndims)
!          print '(a,a,6(1x,i3))', 'tensor counts',name,tensor_counts(tensor_id,:ndims)
!          print '(a,a,6(1x,i3))', 'tensor memcounts',name,tensor_memcounts(tensor_id,:ndims)

!print*, 'before H5Sselect_hyperslab_F for file'
          call H5Sselect_hyperslab_F(tensor_id_S(tensor_id), H5S_SELECT_OR_F, &
                                     tensor_offsets(tensor_id,1:ndims),       &
                                     tensor_counts(tensor_id,1:ndims),        &
                                     hdferr)
          if (hdferr /= 0) &
            call fatal_error('loadDataset_rank2','Error creating File mapping '// &
                             'for /grid/'//name)
          ! Hyperslab for memory
!print*, 'before H5Sselect_hyperslab_F for memory'
          call H5Sselect_hyperslab_F(tensor_id_memS(tensor_id), H5S_SELECT_OR_F, &
                                     tensor_memoffsets(tensor_id,1:ndims),        &
                                     tensor_memcounts(tensor_id,1:ndims),         &
                                     hdferr)
           if (hdferr /= 0) then
             call fatal_error('loadDataset_rank2','Error creating memory mapping '// &
                           'for /grid/'//name)
           end if
        end if
      end do; end do


      ! Read data into memory
      tensor_dims(tensor_id,ndims-1:ndims)=3
!print*, 'before H5Dread_F'
      call H5Dread_F(tensor_id_D(tensor_id), hdf_memtype, dataarray, &
                     tensor_dims(tensor_id,1:ndims), hdferr, &
                     tensor_id_memS(tensor_id), tensor_id_S(tensor_id))
      if (hdferr /= 0) &
        call fatal_error('loadDataset_rank1','Error reading dataset '// &
                         'for /grid/'//name)
      sum = 0.; rms = 0.
      do i=1,3 ; do j=1,3 
        tensor_maxvals(tensor_id) = maxval(dataarray(:,:,:,:,i,j))
        tensor_minvals(tensor_id) = minval(dataarray(:,:,:,:,i,j))
        if (present(name)) then
          call mpireduce_min(tensor_minvals(tensor_id),globmin)
          call mpireduce_max(tensor_maxvals(tensor_id),globmax)
          if (i == j .and. lroot) write (*,*) trim(name)// &
            '(', trim(itoa(i)), ',', trim(itoa(j)),') min/max: ', globmin, globmax
        endif
        sum = sum + tensor_maxvals(tensor_id)*tensor_maxvals(tensor_id)
      enddo ; enddo
      rms=sqrt(sum)
      if (present(name) .and. lroot) write (*,*) trim(name)//' (rms):', rms

      dataarray = tensor_scales(tensor_id) * dataarray

    end subroutine loadDataset_rank2
!*********************************************************************** 
    subroutine loadDataset_rank3(dataarray, datamask, tensor_id, loadstart,name)

      ! Load a chunk of data for a 3-rank tensor, beginning at loadstart

      real, dimension(:,:,:,:,:,:,:), intent(inout) :: dataarray
      logical, dimension(3,3,3), intent(in) :: datamask
      integer, intent(in) :: tensor_id
      integer, intent(in) :: loadstart
      character(*), optional :: name

      integer :: ndims
      integer(HSIZE_T) :: mask_i, mask_j, mask_k
      real :: globmin, globmax ! output for diagnostics

      ndims = tensor_ndims(tensor_id)
      tensor_offsets(tensor_id,1) = loadstart
      call H5Sselect_none_F(tensor_id_S(tensor_id), hdferr)
      call H5Sselect_none_F(tensor_id_memS(tensor_id), hdferr)

      do mask_k=1,3; do mask_j=1,3; do mask_i=1,3
        ! Load only wanted datasets
        if (datamask(mask_i,mask_j,mask_k)) then
          ! Set the new offset for data reading
          tensor_offsets(tensor_id,ndims-2)    = mask_i-1
          tensor_offsets(tensor_id,ndims-1)    = mask_j-1
          tensor_offsets(tensor_id,ndims)      = mask_k-1
          tensor_memoffsets(tensor_id,ndims-2) = mask_i-1
          tensor_memoffsets(tensor_id,ndims-1) = mask_j-1
          tensor_memoffsets(tensor_id,ndims)   = mask_k-1
          ! Hyperslab for data
!          print '(a,a,7(1x,i3))', 'tensor offset',name,tensor_offsets(tensor_id,1:ndims)
!          print '(a,a,7(1x,i3))', 'tensor memoffset',name,tensor_memoffsets(tensor_id,1:ndims)
!          print '(a,a,7(1x,i3))', 'tensor counts',name,tensor_counts(tensor_id,:ndims)
!          print '(a,a,7(1x,i3))', 'tensor memcounts',name,tensor_memcounts(tensor_id,:ndims)

          call H5Sselect_hyperslab_F(tensor_id_S(tensor_id), H5S_SELECT_OR_F, &
                                     tensor_offsets(tensor_id,1:ndims),       &
                                     tensor_counts(tensor_id,1:ndims),        &
                                     hdferr)
          ! Hyperslab for memory
          call H5Sselect_hyperslab_F(tensor_id_memS(tensor_id), H5S_SELECT_OR_F, &
                                     tensor_memoffsets(tensor_id,1:ndims),        &
                                     tensor_memcounts(tensor_id,1:ndims),         &
                                     hdferr)
        end if
      end do; end do; end do
      ! Read data into memory
      tensor_dims(tensor_id,ndims-1:ndims)=3
      call H5Dread_F(tensor_id_D(tensor_id), hdf_memtype, dataarray, &
                     tensor_dims(tensor_id,1:ndims), hdferr, &
                     tensor_id_memS(tensor_id), tensor_id_S(tensor_id))
      tensor_maxvals(tensor_id) = maxval(dataarray)
      tensor_minvals(tensor_id) = minval(dataarray)
      if (present(name)) then
        call mpireduce_min(tensor_minvals(tensor_id),globmin)
        call mpireduce_max(tensor_maxvals(tensor_id),globmax)
        if (lroot) write (*,*) trim(name)//' min/max: ', globmin, globmax
      endif
      dataarray = tensor_scales(tensor_id) * dataarray

    end subroutine loadDataset_rank3
!*********************************************************************** 
    function emf_interpolate(dataarray) result(interp_data)

      real, intent(in), dimension(dataload_len,nx) :: dataarray
      real, dimension(nx) :: interp_data

     ! interp_data=dataarray(:,1)
      interp_data=dataarray(1,:)

    end function emf_interpolate
!*********************************************************************** 
    subroutine setParameterDefaults

      ! alpha
      lalpha=.false.
      lalpha_c=.false.
      lalpha_arr=.false.
      alpha_scale=1.0
      alpha_name='data'
      ! beta
      lbeta=.false.
      lbeta_c=.false.
      lbeta_arr=.false.
      beta_scale=1.0
      beta_name='data'
      ! gamma
      lgamma=.false.
      lgamma_c=.false.
      lgamma_arr=.false.
      gamma_scale=1.0
      gamma_name='data'
      ! delta
      ldelta=.false.
      ldelta_c=.false.
      ldelta_arr=.false.
      delta_scale=1.0
      delta_name='data'
      ! kappa
      lkappa=.false.
      lkappa_c=.false.
      lkappa_arr=.false.
      kappa_scale=1.0
      kappa_name='data'
      ! utensor
      lutensor=.false.
      lutensor_c=.false.
      lutensor_arr=.false.
      utensor_scale=1.0
      utensor_name='data'
      !  meauu
      lmeanuu=.false.
      lmeanuu_c=.false.
      meanuu_scale=1.0
      ! acoef
      lacoef=.false.
      lacoef_c=.false.
      lacoef_arr=.false.
      acoef_scale=1.0
      acoef_name='data'
      ! bcoef
      lbcoef=.false.
      lbcoef_c=.false.
      lbcoef_arr=.false.
      bcoef_scale=1.0
      bcoef_name='data'
      ! other
      interpname  = 'none'
      defaultname = '' 
      tensor_maxvals=0.0
      tensor_minvals=0.0
      lusecoefs    = .false.
      lloop = .false.

    end subroutine setParameterDefaults
!***********************************************************************
    subroutine parseParameters
 
    integer :: i
!
! Load boolean array for alpha
!
      if (any(lalpha_c)) then
        lalpha = .true.
        lalpha_arr(1,1) = lalpha_c(1)
        lalpha_arr(2,1) = lalpha_c(2)
        lalpha_arr(1,2) = lalpha_c(2)
        lalpha_arr(3,1) = lalpha_c(3)
        lalpha_arr(1,3) = lalpha_c(3)
        lalpha_arr(2,2) = lalpha_c(4)
        lalpha_arr(2,3) = lalpha_c(5)
        lalpha_arr(3,2) = lalpha_c(5)
        lalpha_arr(3,3) = lalpha_c(6)
      else if (lalpha) then
        lalpha_arr = .true.
      end if
!
! Load boolean array for beta
!
      if (any(lbeta_c)) then
        lbeta = .true.
        lbeta_arr(1,1) = lbeta_c(1)
        lbeta_arr(2,1) = lbeta_c(2)
        lbeta_arr(1,2) = lbeta_c(2)
        lbeta_arr(3,1) = lbeta_c(3)
        lbeta_arr(1,3) = lbeta_c(3)
        lbeta_arr(2,2) = lbeta_c(4)
        lbeta_arr(2,3) = lbeta_c(5)
        lbeta_arr(3,2) = lbeta_c(5)
        lbeta_arr(3,3) = lbeta_c(6)
      else if (lbeta) then
        lbeta_arr = .true.
      end if
!
! Load boolean array for gamma
!
      if (any(lgamma_c)) then
        lgamma = .true.
        lgamma_arr  = lgamma_c
      else if (lgamma) then
        lgamma_arr  = .true.
      end if
!
! Load boolean array for delta
!
      if (any(ldelta_c)) then
        ldelta = .true.
        ldelta_arr  = ldelta_c
      else if (ldelta) then
        ldelta_arr  = .true.
      end if
!
! Load boolean array for kappa
!
      if (any(lkappa_c)) then
        lkappa = .true.
        do i=1,3
          lkappa_arr(i,1,1) = lkappa_c(i,1)
          lkappa_arr(i,2,1) = lkappa_c(i,2)
          lkappa_arr(i,1,2) = lkappa_c(i,2)
          lkappa_arr(i,3,1) = lkappa_c(i,3)
          lkappa_arr(i,1,3) = lkappa_c(i,3)
          lkappa_arr(i,2,2) = lkappa_c(i,4)
          lkappa_arr(i,2,3) = lkappa_c(i,5)
          lkappa_arr(i,3,2) = lkappa_c(i,5)
          lkappa_arr(i,3,3) = lkappa_c(i,6)
        enddo
      elseif (lkappa) then
        lkappa_arr = .true.
      end if
!
! Load boolean array for acoef
!
      if (any(lacoef_c)) then
        lacoef=.true.
        lacoef_arr(1,1) = lacoef_c(1)
        lacoef_arr(2,1) = lacoef_c(2)
        lacoef_arr(1,2) = lacoef_c(2)
        lacoef_arr(3,1) = lacoef_c(3)
        lacoef_arr(1,3) = lacoef_c(3)
        lacoef_arr(2,2) = lacoef_c(4)
        lacoef_arr(2,3) = lacoef_c(5)
        lacoef_arr(3,2) = lacoef_c(5)
        lacoef_arr(3,3) = lacoef_c(6)
        if (any([lalpha,lgamma])) then
          if (lroot) call warning('initialize_special', &
            'any lacoef_c=T overrides settings of lalpha and lgamma')     
          lalpha=.false.; lgamma=.false.
          lalpha_arr = .false.; lgamma_arr = .false.
        endif
      elseif (lacoef) then
        lacoef_arr = .true.
      end if
!
! Load boolean array for bcoef
!
      if (any(lbcoef_c)) then
        lbcoef = .true.
        do i=1,3
          lbcoef_arr(i,1,1) = lbcoef_c(i,1)
          lbcoef_arr(i,2,1) = lbcoef_c(i,2)
          lbcoef_arr(i,1,2) = lbcoef_c(i,2)
          lbcoef_arr(i,3,1) = lbcoef_c(i,3)
          lbcoef_arr(i,1,3) = lbcoef_c(i,3)
          lbcoef_arr(i,2,2) = lbcoef_c(i,4)
          lbcoef_arr(i,2,3) = lbcoef_c(i,5)
          lbcoef_arr(i,3,2) = lbcoef_c(i,5)
          lbcoef_arr(i,3,3) = lbcoef_c(i,6)
        enddo
        if (any([lbeta,ldelta,lkappa])) then
          if (lroot) call warning('initialize_special', &
            'any lbcoef_c=T overrides settings of lbeta,ldelta,lkappa')     
          lbeta=.false.; lbeta_arr = .false.
          ldelta=.false.; ldelta_arr = .false.
          lkappa=.false.; lkappa_arr = .false.
        endif
      elseif (lbcoef) then
        lbcoef_arr = .true.
      end if
!
! Load boolean array for utensor
!
      if (any(lutensor_c)) then
        lutensor = .true.
        lutensor_arr  = lutensor_c
      else if (lutensor) then
        lutensor_arr  = .true.
      end if

      if (any(lmeanuu_c)) then
        lmeanuu = .true.
        lutensor_arr  = lmeanuu_c
      else if (lmeanuu) then
        lutensor_arr  = .true.
      end if

      if (lmeanuu .and. (meanuu_scale /= 1.0)) then
        utensor_scale=meanuu_scale
      endif
!
! Store scales
!
      tensor_scales(alpha_id)   = alpha_scale
      tensor_scales(beta_id)    = beta_scale
      tensor_scales(gamma_id)   = gamma_scale
      tensor_scales(delta_id)   = delta_scale
      tensor_scales(kappa_id)   = kappa_scale
      tensor_scales(utensor_id) = utensor_scale
      tensor_scales(acoef_id)   = acoef_scale
      tensor_scales(bcoef_id)   = bcoef_scale
!
! Store names
!
      if (trim(defaultname) /= '') then
        alpha_name    = defaultname
        beta_name     = defaultname
        gamma_name    = defaultname
        delta_name    = defaultname
        kappa_name    = defaultname
        utensor_name  = defaultname
        acoef_name    = defaultname
        bcoef_name    = defaultname
      end if
      tensor_names(alpha_id)    = alpha_name
      tensor_names(beta_id)     = beta_name
      tensor_names(gamma_id)    = gamma_name
      tensor_names(delta_id)    = delta_name
      tensor_names(kappa_id)    = kappa_name
      tensor_names(utensor_id)  = utensor_name
      tensor_names(acoef_id)    = acoef_name
      tensor_names(bcoef_id)    = bcoef_name
            
    end subroutine parseParameters
!*****************************************************************************
    logical function output_persistent_special()
!
!  Writes out the time of the next SNI
!
!  13-Dec-2011/Bourdin.KIS: reworked
!  14-jul-2015/fred: removed obsolete Remnant persistant variable from current
!  write and added new cluster variables. All now consistent with any io
!
      use IO, only: write_persist
!
!      if (lcollective_IO) call fatal_error ('output_persistent_interstellar', &
!          "The interstellar persistent variables can't be written
!          collectively!")
!
      output_persistent_special = .true.
!
      if (write_persist ('SPECIAL_ILOAD', id_record_SPECIAL_ILOAD, iload)) return
!
      output_persistent_special = .false.
!
    endfunction output_persistent_special
!*****************************************************************************
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!*********************************************************************** 
endmodule Special
