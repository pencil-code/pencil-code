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
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
  use Mpicomm, only: mpibarrier
  use Sub, only: numeric_precision, dot_mn, dot_mn_vm, curl_mn, cross_mn, vec_dot_3tensor,dot2_mn
  use HDF5
  use File_io, only: parallel_unit
!
!
  implicit none
!
  include 'mpif.h'
  include '../special.h'
  ! HDF debug parameters:
  integer :: hdferr
  logical :: hdf_exists
  ! Main HDF object ids and variables
  character (len=fnlen) :: hdf_emftensors_filename
  integer(HID_T) :: hdf_memtype, &
                    hdf_emftensors_file,  &
                    hdf_emftensors_plist, &
                    hdf_emftensors_group
  ! Dataset HDF object ids
  integer, parameter :: ntensors = 8
  integer(HID_T),   dimension(ntensors) :: tensor_id_G, tensor_id_D, &
                                    tensor_id_S, tensor_id_memS
  integer(HSIZE_T), dimension(ntensors,10) :: tensor_dims, &
                                      tensor_offsets, tensor_counts, &
                                      tensor_memdims, &
                                      tensor_memoffsets, tensor_memcounts
                                      
  ! Actual datasets
  real, dimension(:,:,:,:,:,:)  , allocatable :: alpha_data, beta_data, &
                                                 acoef_data
  real, dimension(:,:,:,:,:)    , allocatable :: gamma_data, delta_data, &
                                                 utensor_data
  real, dimension(:,:,:,:,:,:,:), allocatable :: kappa_data, bcoef_data
  ! Dataset mappings
  integer,parameter :: alpha_id=1, beta_id=2,    &
                       gamma_id=3, delta_id=4,   &
                       kappa_id=5, utensor_id=6, & 
                       acoef_id=7, bcoef_id=8
  integer, dimension(ntensors),parameter :: tensor_ndims = (/ 6, 6, 5, 5, 7, 5, 6, 7 /) 
  ! Dataset logical variables
  logical, dimension(3,3)   :: lalpha_arr, lbeta_arr, lacoef_arr
  logical, dimension(3)     :: lgamma_arr, ldelta_arr, lutensor_arr
  logical, dimension(3,3,3) :: lkappa_arr, lbcoef_arr
  logical, dimension(6)     :: lalpha_c, lbeta_c, lacoef_c
  logical, dimension(3)     :: lgamma_c, ldelta_c, lutensor_c
  logical, dimension(3,3,3) :: lkappa_c, lbcoef_c
  logical :: lalpha, lbeta, lgamma, ldelta, lkappa, lutensor, lacoef, lbcoef, lusecoefs
  real :: alpha_scale, beta_scale, gamma_scale, delta_scale, kappa_scale, utensor_scale, acoef_scale, bcoef_scale
  character (len=fnlen) :: defaultname
  character (len=fnlen) :: alpha_name, beta_name,    &
                           gamma_name, delta_name,   &
                           kappa_name, utensor_name, &
                           acoef_name, bcoef_name
  character (len=fnlen), dimension(ntensors) :: tensor_names
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
  ! Interpolation parameters
  character (len=fnlen) :: interpname
  integer :: dataload_len
  real, dimension(nx)     :: tmpline
  real, dimension(nx,3)   :: tmppencil,emftmp
  real, dimension(nx,3,3) :: tmptensor
  ! Input dataset name
  ! Input namelist
  namelist /special_init_pars/ &
      lalpha,   lalpha_c,   alpha_name,   alpha_scale, &
      lbeta,    lbeta_c,    beta_name,    beta_scale, &
      lgamma,   lgamma_c,   gamma_name,   gamma_scale, &
      ldelta,   ldelta_c,   delta_name,   delta_scale, &
      lkappa,   lkappa_c,   kappa_name,   kappa_scale, &
      lutensor, lutensor_c, utensor_name, utensor_scale, &
      lacoef,   lacoef_c,   acoef_name,   acoef_scale, &
      lbcoef,   lbcoef_c,   bcoef_name,   bcoef_scale, &
      interpname, defaultname, lusecoefs
  namelist /special_run_pars/ &
      lalpha,   lalpha_c,   alpha_name,   alpha_scale, &
      lbeta,    lbeta_c,    beta_name,    beta_scale, &
      lgamma,   lgamma_c,   gamma_name,   gamma_scale, &
      ldelta,   ldelta_c,   delta_name,   delta_scale, &
      lkappa,   lkappa_c,   kappa_name,   kappa_scale, &
      lutensor, lutensor_c, utensor_name, utensor_scale, &
      lacoef,   lacoef_c,   acoef_name,   acoef_scale, &
      lbcoef,   lbcoef_c,   bcoef_name,   bcoef_scale, &
      interpname, defaultname, lusecoefs
  ! loadDataset interface
  interface loadDataset
    module procedure loadDataset_rank1
    module procedure loadDataset_rank2
    module procedure loadDataset_rank3
  end interface

  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
      integer :: i
!
      write(*,*) ' register_special running...'
      if (lroot) call svn_id( &
           "$Id$")
      if (trim(interpname) == 'none') then
        dataload_len = 1
      !else
      !  call fatal_error('register_special','Unknown interpolation chosen!')
      end if
      ! Set dataset offsets
      tensor_offsets = 0
      do i=1,ntensors
        tensor_offsets(i,1:3) = [ ipx*nx, ipy*ny , ipz*nz ]
      end do
      ! Set dataset counts
      tensor_counts = 1
      do i=1,ntensors
        tensor_counts(i,1:4) = [ nx, ny, nz, dataload_len ]
      end do
      ! Set dataset memory offsets
      tensor_memoffsets = 0
      ! Set dataset memory counts
      tensor_memcounts = 1
      do i=1,ntensors
        tensor_memcounts(i,1:4) = [ nx, ny, nz, dataload_len ]
      end do
      ! Set memory dimensions
      tensor_memdims = 3
      do i=1,ntensors
        tensor_memdims(i,1:4) = [ nx, ny, nz, dataload_len ]
      end do
!
!!      call farray_register_pde('special',ispecial)
!!      call farray_register_auxiliary('specaux',ispecaux)
!!      call farray_register_auxiliary('specaux',ispecaux,communicated=.true.)
!
    endsubroutine register_special
!!***********************************************************************
    subroutine register_particles_special(npvar)
!
!  Set up indices for particle variables in special modules.
!
!  4-jan-14/tony: coded
!
      integer :: npvar
!
      if (lroot) call svn_id( &
           "$Id$")
      call keep_compiler_quiet(npvar)
!
!
!!      iqp=npvar+1
!!      npvar=npvar+1
!
    endsubroutine register_particles_special
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

      write(*,*) 'initialize_special running...'

      if (lrun) then 

        hdf_emftensors_plist = -1
        hdf_emftensors_file  = -1

        call H5open_F(hdferr)

        hdf_memtype=-1
        if (numeric_precision() == 'S') then
          write(*,*) 'initialize special: loading data as single precision'
          hdf_memtype = H5T_NATIVE_REAL
        else
          write(*,*) 'initialize special: loading data as double precision'
          hdf_memtype = H5T_NATIVE_DOUBLE
        end if
        
        write (*,*) 'initialize special: setting parallel HDF5 IO for data file reading'
        call H5Pcreate_F(H5P_FILE_ACCESS_F, hdf_emftensors_plist, hdferr)
        call H5Pset_fapl_mpio_F(hdf_emftensors_plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

        print *, 'initialize special: opening emftensors.h5 and loading relevant fields of it into memory...'

        hdf_emftensors_filename = trim(datadir_snap)//'/emftensors.h5'

        inquire(file=hdf_emftensors_filename, exist=hdf_exists)

        if (.not. hdf_exists) then
          call H5Pclose_F(hdf_emftensors_plist, hdferr)
          call H5close_F(hdferr)
          call fatal_error('initialize_special','File '//trim(hdf_emftensors_filename)//' does not exist!')
        end if

        call H5Fopen_F(hdf_emftensors_filename, H5F_ACC_RDONLY_F, hdf_emftensors_file, hdferr, access_prp = hdf_emftensors_plist)

        call H5Lexists_F(hdf_emftensors_file,'/emftensor/', hdf_exists, hdferr)

        if (.not. hdf_exists) then
          call H5Fclose_F(hdf_emftensors_file, hdferr)
          call H5Pclose_F(hdf_emftensors_plist, hdferr)
          call H5close_F(hdferr)
          call fatal_error('initialize_special','group /emftensor/ does not exist!')
        end if

        call H5Gopen_F(hdf_emftensors_file, 'emftensor', hdf_emftensors_group, hdferr)

        ! Open datasets

        ! alpha
        if (lalpha) then
          call openDataset('alpha', alpha_id)
          write(*,*) 'initialize_special: Using dataset /emftensor/alpha/'//trim(alpha_name)//' for alpha'
        end if
        ! beta
        if (lbeta) then
          call openDataset('beta', beta_id)
          write(*,*) 'initialize_special: Using dataset /emftensor/beta/'//trim(beta_name)//' for beta'
        end if
        ! gamma
        if (lgamma) then
          call openDataset('gamma', gamma_id)
          write(*,*) 'initialize_special: Using dataset /emftensor/gamma/'//trim(gamma_name)//' for gamma'
        end if
        ! delta
        if (ldelta) then
          call openDataset('delta', delta_id)
          write(*,*) 'initialize_special: Using dataset /emftensor/delta/'//trim(delta_name)//' for delta'
        end if
        ! kappa
        if (lkappa) then
          call openDataset('kappa', kappa_id)
          write(*,*) 'initialize_special: Using dataset /emftensor/kappa/'//trim(kappa_name)//' for kappa'
        end if
        ! utensor
        if (lutensor) then
          call openDataset('utensor', utensor_id)
          write(*,*) 'initialize_special: Using dataset /emftensor/utensor/'//trim(utensor_name)//' for utensor'
        end if
        ! acoef
        if (lacoef) then
          call openDataset('acoef', acoef_id)
          write(*,*) 'initialize_special: Using dataset /emftensor/acoef/'//trim(acoef_name)//' for acoef'
        end if
        ! bcoef
        if (lbcoef) then
          call openDataset('bcoef', bcoef_id)
          write(*,*) 'initialize_special: Using dataset /emftensor/bcoef/'//trim(bcoef_name)//' for bcoef'
        end if
        
        
        
        ! Load initial dataset values

        ! Allocate data arrays
        allocate(alpha_data(nx,ny,nz,dataload_len,3,3))
        allocate(beta_data(nx,ny,nz,dataload_len,3,3))
        allocate(gamma_data(nx,ny,nz,dataload_len,3))
        allocate(delta_data(nx,ny,nz,dataload_len,3))
        allocate(kappa_data(nx,ny,nz,dataload_len,3,3,3))
        allocate(utensor_data(nx,ny,nz,dataload_len,3))
        allocate(acoef_data(nx,ny,nz,dataload_len,3,3))
        allocate(bcoef_data(nx,ny,nz,dataload_len,3,3,3))
        alpha_data    = 0
        beta_data     = 0
        gamma_data    = 0
        delta_data    = 0
        kappa_data    = 0
        utensor_data  = 0
        acoef_data    = 0
        bcoef_data    = 0
        ! Load datasets
        if (lalpha) then
          call loadDataset(alpha_data, lalpha_arr, alpha_id, 0)
          if (lroot) then
              write (*,*) 'Alpha scale:  ', alpha_scale
              write (*,*) 'Alpha maxval: ', maxval(alpha_data)
              write (*,*) 'Alpha components used: '
              do i=1,3
                  write (*,'(A3,3L3,A3)') '|', lalpha_arr(:,i), '|'
              end do
          end if
        end if
        if (lbeta) then
          call loadDataset(beta_data, lbeta_arr, beta_id, 0)
          if (lroot) then
              write (*,*) 'Beta scale:  ', beta_scale
              write (*,*) 'Beta maxval: ', maxval(beta_data)
              write (*,*) 'Beta components used: '
              do i=1,3
                  write (*,'(A3,3L3,A3)') '|', lbeta_arr(:,i), '|'
              end do
          end if
        end if
        if (lgamma) then
          call loadDataset(gamma_data, lgamma_arr, gamma_id, 0)
          if (lroot) then
              write (*,*) 'Gamma scale:  ', gamma_scale
              write (*,*) 'Gamma maxval: ', maxval(gamma_data)
              write (*,*) 'Gamma components used: '
              write (*,'(A3,3L3,A3)') '|', lgamma_arr, '|'
          end if
        end if
        if (ldelta) then
          call loadDataset(delta_data, ldelta_arr, delta_id, 0)
          if (lroot) then
              write (*,*) 'Delta scale:  ', delta_scale
              write (*,*) 'Delta maxval: ', maxval(delta_data)
              write (*,*) 'Delta components used: '
              write (*,'(A3,3L3,A3)') '|', ldelta_arr, '|'
          end if
        end if
        if (lkappa) then
          call loadDataset(kappa_data, lkappa_arr, kappa_id, 0)
          if (lroot) then
              write (*,*) 'Kappa scale:  ', kappa_scale
              write (*,*) 'Kappa maxval: ', maxval(kappa_data)
              write (*,*) 'Kappa components used: '
              do j=1,3
                write(*,*) '|'
                do i=1,3
                  write (*,'(A3,3L3,A3)') '||', lkappa_arr(:,j,i), '||'
                end do
              end do
              write(*,*) '|'
          end if
        end if
        if (lutensor) then
          call loadDataset(utensor_data, lutensor_arr, utensor_id, 0)
          if (lroot) then
              write (*,*) 'U-tensor scale:  ', utensor_scale
              write (*,*) 'U-tensor maxval: ', maxval(utensor_data)
              write (*,*) 'U-tensor components used: '
              write (*,'(A3,3L3,A3)') '|', lutensor_arr, '|'
          end if
        end if
        if (lacoef) then
          call loadDataset(acoef_data, lacoef_arr, acoef_id, 0)
          if (lroot) then
              write (*,*) 'acoef scale:  ', acoef_scale
              write (*,*) 'acoef maxval: ', maxval(acoef_data)
              write (*,*) 'acoef components used: '
              do i=1,3
                  write (*,'(A3,3L3,A3)') '|', lacoef_arr(:,i), '|'
              end do
          end if
        end if
        if (lbcoef) then
          call loadDataset(bcoef_data, lbcoef_arr, bcoef_id, 0)
          if (lroot) then
              write (*,*) 'bcoef scale:  ', bcoef_scale
              write (*,*) 'bcoef maxval: ', maxval(bcoef_data)
              write (*,*) 'bcoef components used: '
              do j=1,3
                write(*,*) '|'
                do i=1,3
                  write (*,'(A3,3L3,A3)') '||', lbcoef_arr(:,j,i), '||'
                end do
              end do
          end if
        end if
      end if
  !
    endsubroutine initialize_special
!***********************************************************************
    subroutine finalize_special(f)
!
!  Called right before exiting.
!
!  14-aug-2011/Bourdin.KIS: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
      if (lrun) then
        ! Deallocate data
        if (allocated(alpha_data)) then
          deallocate(alpha_data)
        end if
        if (allocated(beta_data)) then
          deallocate(beta_data)
        end if
        if (allocated(gamma_data)) then
          deallocate(gamma_data)
        end if
        if (allocated(delta_data)) then
          deallocate(delta_data)
        end if
        if (allocated(kappa_data)) then
          deallocate(kappa_data)
        end if
        if (allocated(utensor_data)) then
          deallocate(utensor_data)
        end if

        print *,'Closing emftensors.h5'
        
        if (lalpha) then
          call closeDataset(alpha_id)
        end if
        if (lbeta) then
          call closeDataset(beta_id)
        end if
        if (lgamma) then
          call closeDataset(gamma_id)
        end if
        if (ldelta) then
          call closeDataset(delta_id)
        end if
        if (lkappa) then
          call closeDataset(kappa_id)
        end if
        if (lutensor) then
          call closeDataset(utensor_id)
        end if
        if (lacoef) then
          call closeDataset(acoef_id)
        end if
        if (lbcoef) then
          call closeDataset(bcoef_id)
        end if

        call H5Gclose_F(hdf_emftensors_group, hdferr)
        call H5Fclose_F(hdf_emftensors_file, hdferr)
        call H5Pclose_F(hdf_emftensors_plist, hdferr)
        call H5close_F(hdferr)
        call mpibarrier()
      end if
  !
    endsubroutine finalize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!!
!!  SAMPLE IMPLEMENTATION
!!
!!      select case (initspecial)
!!        case ('nothing'); if (lroot) print*,'init_special: nothing'
!!        case ('zero', '0'); f(:,:,:,iSPECIAL_VARIABLE_INDEX) = 0.
!!        case default
!!          call fatal_error("init_special: No such value for initspecial:" &
!!              ,trim(initspecial))
!!      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_bij)=.true.
      lpenc_requested(i_jj)=.true.
      lpenc_requested(i_bij_symm)=.true.
      lpenc_requested(i_alpha_coefs)=.true.
      lpenc_requested(i_beta_coefs)=.true.
      lpenc_requested(i_gamma_coefs)=.true.
      lpenc_requested(i_delta_coefs)=.true.
      lpenc_requested(i_kappa_coefs)=.true.
      lpenc_requested(i_utensor_coefs)=.true.
      lpenc_requested(i_acoef_coefs)=.true.
      lpenc_requested(i_bcoef_coefs)=.true.
      lpenc_requested(i_alpha_emf)=.true.
      lpenc_requested(i_beta_emf)=.true.
      lpenc_requested(i_gamma_emf)=.true.
      lpenc_requested(i_delta_emf)=.true.
      lpenc_requested(i_kappa_emf)=.true.
      lpenc_requested(i_utensor_emf)=.true.
      lpenc_requested(i_acoef_emf)=.true.
      lpenc_requested(i_bcoef_emf)=.true.
      lpenc_requested(i_emf)=.true.


      write(*,*) 'pencil_criteria_special: Pencils requested'
!
!  18-07-06/tony: coded
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      integer i,j,k
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
      if (lalpha) then
        ! Calculate alpha B
        do j=1,3; do i=1,3
          if (lalpha_arr(i,j)) then
            p%alpha_coefs(1:nx,i,j)=emf_interpolate(alpha_data(1:nx,m-nghost,n-nghost,1:dataload_len,i,j))
          else
            p%alpha_coefs(1:nx,i,j)=0
          end if
        end do; end do
        call dot_mn_vm(p%bb,p%alpha_coefs,p%alpha_emf)
      end if
      if (lbeta) then
        ! Calculate beta (curl B)
        do j=1,3; do i=1,3
          if (lbeta_arr(i,j)) then
            p%beta_coefs(1:nx,i,j)=emf_interpolate(beta_data(1:nx,m-nghost,n-nghost,1:dataload_len,i,j))
          else
            p%beta_coefs(1:nx,i,j)=0
          end if
        end do; end do
        call dot_mn_vm(p%jj,p%beta_coefs,p%beta_emf)
      end if
      if (lgamma) then
        ! Calculate gamma x B
        do i=1,3
          if (lgamma_arr(i)) then
            p%gamma_coefs(1:nx,i)=emf_interpolate(gamma_data(1:nx,m-nghost,n-nghost,1:dataload_len,i))
          else
            p%gamma_coefs(1:nx,i)=0
          end if
        end do
        call cross_mn(p%gamma_coefs,p%bb,p%gamma_emf)
      end if
      if (ldelta) then
        ! Calculate delta x (curl B)
        do i=1,3
          if (ldelta_arr(i)) then
            p%delta_coefs(1:nx,i)=emf_interpolate(delta_data(1:nx,m-nghost,n-nghost,1:dataload_len,i))
          else
            p%delta_coefs(1:nx,i)=0
          end if
        end do
        call cross_mn(p%delta_coefs,p%jj,p%delta_emf)
      end if
      if (lkappa) then
        ! Calculate kappa (grad B)_symm
        do j=1,3; do i=1,3
          p%bij_symm(1:nx,i,j)=0.5*(p%bij(1:nx,i,j) + p%bij(1:nx,j,i))
        end do; end do
        do k=1,3; do j=1,3; do i=1,3
          p%kappa_coefs(1:nx,i,j,k)=emf_interpolate(kappa_data(1:nx,m-nghost,n-nghost,1:dataload_len,i,j,k))
        end do; end do; end do
        p%kappa_emf = 0
        do k=1,3; do j=1,3; do i=1,3
          p%kappa_emf(1:nx,i)=p%kappa_emf(1:nx,i)+p%kappa_coefs(1:nx,i,j,k)*p%bij_symm(1:nx,j,k)
        end do; end do; end do
      end if
      if (lutensor) then
        ! Calculate utensor x B
        do i=1,3
          if (lutensor_arr(i)) then
            p%utensor_coefs(1:nx,i)=emf_interpolate(utensor_data(1:nx,m-nghost,n-nghost,1:dataload_len,i))
          else
            p%utensor_coefs(1:nx,i)=0
          end if
        end do
        call cross_mn(p%utensor_coefs,p%bb,p%utensor_emf)
      end if
      if (lacoef) then
        ! Calculate acoef B
        do j=1,3; do i=1,3
          if (lacoef_arr(i,j)) then
            p%acoef_coefs(1:nx,i,j)=emf_interpolate(acoef_data(1:nx,m-nghost,n-nghost,1:dataload_len,i,j))
          else
            p%acoef_coefs(1:nx,i,j)=0
          end if
        end do; end do
        call dot_mn_vm(p%bb,p%acoef_coefs,p%acoef_emf)
      end if
      if (lbcoef) then
        ! Calculate bcoef (grad B)
        do k=1,3; do j=1,3; do i=1,3
          p%bcoef_coefs(1:nx,i,j,k)=emf_interpolate(bcoef_data(1:nx,m-nghost,n-nghost,1:dataload_len,i,j,k))
        end do; end do; end do
        p%bcoef_emf = 0
        do k=1,3; do j=1,3; do i=1,3
          p%bcoef_emf(1:nx,i)=p%bcoef_emf(1:nx,i)+p%bcoef_coefs(1:nx,i,j,k)*p%bij(1:nx,j,k)
        end do; end do; end do
      end if
!
! Calculate emf pencil
!
      p%emf = 0
      if (lalpha) then
        p%emf = p%emf + p%alpha_emf
      end if
      if (lbeta) then
        p%emf = p%emf - p%beta_emf
      end if
      if (lgamma) then
        p%emf = p%emf + p%gamma_emf
      end if
      if (ldelta) then
        p%emf = p%emf - p%delta_emf
      end if
      if (lkappa) then
        p%emf = p%emf - p%kappa_emf
      end if
      if (lutensor) then
        p%emf = p%emf + p%utensor_emf
      end if
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
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!!      if (headtt) call identify_bcs('special',ispecial)
!
!!
!! SAMPLE DIAGNOSTIC IMPLEMENTATION
!!
!!      if (ldiagnos) then
!!        if (idiag_SPECIAL_DIAGNOSTIC/=0) then
!!          call sum_mn_name(MATHEMATICAL EXPRESSION,idiag_SPECIAL_DIAGNOSTIC)
!!! see also integrate_mn_name
!!        endif
!!      endif
!      emftmp=0
      if (ldiagnos) then
        emftmp = p%acoef_emf + p%bcoef_emf + p%utensor_emf
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
    
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      integer, intent(out) :: iostat
      write (*,*) 'read_special_init_pars running...'
!
      iostat = 0
      call setParameterDefaults()
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
      write (*,*) 'read_special_init_pars parameters read...'
      call parseParameters()
      write (*,*) 'read_special_init_pars parameters parsed...'
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
      
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
      
      write (*,*) 'read_special_run_pars running...'
      call setParameterDefaults()
      write (*,*) 'read_special_run_pars parameters read...'
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
      call parseParameters()
      write (*,*) 'read_special_run_pars parameters parsed...'
!
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
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
        ! RMS diagnostics
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
      endif
!!
!!      do iname=1,nname
!!        call parse_name(iname,cname(iname),cform(iname),&
!!            'NAMEOFSPECIALDIAGNOSTIC',idiag_SPECIAL_DIAGNOSTIC)
!!      enddo
!!
!!!  write column where which magnetic variable is stored
!!      if (lwr) then
!!        write(3,*) 'idiag_SPECIAL_DIAGNOSTIC=',idiag_SPECIAL_DIAGNOSTIC
!!      endif
!!
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
      enddo
      
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  momentum equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  continuity equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_dustdensity(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  continuity equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_dustdensity
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  energy equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_energy
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
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
! Overwrite with a and b coefs if needed
!
      if (lusecoefs) then
        emftmp=0
        if (lacoef) then
          emftmp = emftmp + p%acoef_emf
        end if
        if (lbcoef) then
          emftmp = emftmp + p%bcoef_emf
        end if
        if (lutensor) then
          emftmp = emftmp + p%utensor_emf
        end if
      else
        emftmp = p%emf
      end if

      df(l1:l2,m,n,iax:iaz)=df(l1:l2,m,n,iax:iaz)+emftmp
!
      if (lfirst.and.ldt) then
!
! Calculate advec_special
!
        if (lalpha) then
          tmppencil=0
          call dot_mn_vm(dline_1, p%alpha_coefs, tmppencil)
          advec_special=advec_special+&
                            tmppencil(:,1)+& 
                            tmppencil(:,2)+& 
                            tmppencil(:,3)
        end if
        if (lgamma) then
          tmppencil=0
          call dot_mn(dline_1, p%gamma_coefs, tmpline)
          advec_special=advec_special+tmpline
        end if
        if (lutensor) then
          tmppencil=0
          call dot_mn(dline_1, p%utensor_coefs, tmpline)
          advec_special=advec_special+tmpline
        end if
!
! Calculate diffus_special
!
        if (lbeta) then
          tmppencil=0
          tmpline=0
          call dot_mn_vm(dline_1,p%beta_coefs, tmppencil)
          call dot_mn(dline_1, tmppencil, tmpline)
          diffus_special=diffus_special+tmpline
        end if
        
        if (ldelta) then
          tmppencil=0
          tmpline=0
          call cross_mn(dline_1,p%delta_coefs, tmppencil)
          call dot_mn(dline_1,tmppencil,tmpline)
          diffus_special=diffus_special+tmpline
        end if
        
        if (lkappa) then
          tmppencil=0
          tmpline=0
          tmptensor=0
          call vec_dot_3tensor(dline_1, p%kappa_coefs, tmptensor)
          call dot_mn_vm(dline_1,tmptensor, tmppencil)
          diffus_special=diffus_special+&
                            tmppencil(:,1)+& 
                            tmppencil(:,2)+& 
                            tmppencil(:,3)
        end if
      end if 
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine special_calc_pscalar(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  passive scalar equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  15-jun-09/anders: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ilncc) = df(l1:l2,m,n,ilncc) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_pscalar
!***********************************************************************
    subroutine special_calc_particles(f,fp,ineargrid)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (:,:), intent(in) :: fp
      integer, dimension(:,:) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine special_calc_particles
!***********************************************************************
    subroutine special_particles_bfre_bdary(f,fp,ineargrid)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (:,:), intent(in) :: fp
      integer, dimension(:,:) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine special_particles_bfre_bdary
!***********************************************************************
    subroutine special_calc_chemistry(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  induction equation.
!
!
!  15-sep-10/natalia: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_chemistry
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
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dt_)
!
    endsubroutine  special_after_timestep
!***********************************************************************
subroutine set_init_parameters(Ntot,dsize,init_distr,init_distr2)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(ndustspec) :: dsize,init_distr,init_distr2
      real :: Ntot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(dsize,init_distr,init_distr2)
      call keep_compiler_quiet(Ntot)
!
    endsubroutine  set_init_parameters
!***********************************************************************

    subroutine openDataset(datagroup,      &
                           tensor_id)
      ! Open a dataset e.g. /emftensor/alpha/data and auxillary dataspaces
      character(len=*), intent(in)    :: datagroup
      integer, intent(in)             :: tensor_id
!
      integer(HSIZE_T), dimension(10) :: maxdimsizes
      integer :: ndims
      logical :: hdf_exists
      character(len=fnlen)            :: dataname
      
      dataname = tensor_names(tensor_id)
      ! Check that datagroup e.g. /emftensor/alpha exists
      call H5Lexists_F(hdf_emftensors_group, datagroup, hdf_exists, hdferr)
      if (.not. hdf_exists) then
        call fatal_error('openDataset','/emftensor/'//trim(datagroup)// &
                          ' does not exist')
      end if
      ! Open datagroup
      call H5Gopen_F(hdf_emftensors_group, datagroup, tensor_id_G(tensor_id), hdferr)
      if (hdferr /= 0) then
        call fatal_error('openDataset','Error opening /emftensor/'//trim(datagroup))
      end if
      ! Check that dataset e.g. /emftensor/alpha/data exists
      call H5Lexists_F(tensor_id_G(tensor_id), dataname, hdf_exists, hdferr)
      if (.not. hdf_exists) then
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        call fatal_error('openDataset','/emftensor/'//trim(datagroup)// &
                          '/'//trim(dataname)//' does not exist')
      end if
      ! Open dataset
      call H5Dopen_F(tensor_id_G(tensor_id), dataname, tensor_id_D(tensor_id), hdferr)
      if (hdferr /= 0) then
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        call fatal_error('openDataset','Error opening /emftensor/'// &
                          trim(datagroup)//'/'//trim(dataname))
      end if
      ! Get dataspace
      call H5Dget_space_F(tensor_id_D(tensor_id), tensor_id_S(tensor_id), hdferr)
      if (hdferr /= 0) then
        call H5Dclose_F(tensor_id_D(tensor_id), hdferr)
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        call fatal_error('openDataset','Error opening dataspace for /emftensor/'// &
                          trim(datagroup)//'/'//trim(dataname))
      end if
      ! Get dataspace dimensions
      ndims = tensor_ndims(tensor_id)
      call H5Sget_simple_extent_dims_F(tensor_id_S(tensor_id), &
                                       tensor_dims(tensor_id,1:ndims), &
                                       maxdimsizes(1:ndims), &
                                       hdferr)
      ! Create a memory space mapping for input data
      call H5Screate_simple_F(ndims, tensor_memdims(tensor_id,1:ndims), &
                              tensor_id_memS(tensor_id), hdferr)
      if (hdferr /= 0) then
        call H5Sclose_F(tensor_id_S(tensor_id), hdferr)
        call H5Dclose_F(tensor_id_D(tensor_id), hdferr)
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        call fatal_error('openDataset','Error creating memory mapping & 
                          for /emftensor/'//trim(datagroup)//'/'//trim(dataname))
      end if

    end subroutine openDataset

    subroutine closeDataset(tensor_id)
      ! Close opened dataspaces, dataset and group
      implicit none
      integer :: tensor_id
      call H5Sclose_F(tensor_id_memS(tensor_id), hdferr)
      call H5Sclose_F(tensor_id_S(tensor_id), hdferr)
      call H5Dclose_F(tensor_id_D(tensor_id), hdferr)
      call H5Gclose_F(tensor_id_G(tensor_id), hdferr)

    end subroutine closeDataset

    subroutine loadDataset_rank1(dataarray, datamask, tensor_id, loadstart)
      ! Load a chunk of data for a vector, beginning at loadstart
      implicit none
      real, dimension(:,:,:,:,:), intent(inout) :: dataarray
      logical, dimension(3), intent(in) :: datamask
      integer, intent(in) :: tensor_id
      integer, intent(in) :: loadstart
      integer :: ndims
      integer(HSIZE_T) :: mask_i
      ndims = tensor_ndims(tensor_id)
      tensor_offsets(tensor_id,4) = loadstart
      call H5Sselect_none_F(tensor_id_S(tensor_id), hdferr)
      call H5Sselect_none_F(tensor_id_memS(tensor_id), hdferr)
      do mask_i=1,3
        ! Load only wanted datasets
        if (datamask(mask_i)) then
          ! Set the new offset for data reading
          tensor_offsets(tensor_id,ndims)    = mask_i-1
          tensor_memoffsets(tensor_id,ndims) = mask_i-1
          ! Hyperslab for data
          call H5Sselect_hyperslab_F(tensor_id_S(tensor_id), H5S_SELECT_OR_F, &
                                     tensor_offsets(tensor_id,1:ndims),       &
                                     tensor_counts(tensor_id,1:ndims),        &
                                     hdferr)
          ! Hyperslab for memory
          call H5Sselect_hyperslab_F(tensor_id_memS(tensor_id), H5S_SELECT_OR_F, &
                                     tensor_memoffsets(tensor_id,1:ndims),       &
                                     tensor_memcounts(tensor_id,1:ndims),        &
                                     hdferr)
        end if
      end do
      ! Read data into memory
      call H5Dread_F(tensor_id_D(tensor_id), hdf_memtype, dataarray, &
                     tensor_dims(tensor_id,1:ndims), hdferr, &
                     tensor_id_memS(tensor_id), tensor_id_S(tensor_id))
      dataarray = tensor_scales(tensor_id) * dataarray
      tensor_maxvals(tensor_id) = maxval(dataarray)
      tensor_minvals(tensor_id) = minval(dataarray)
    end subroutine loadDataset_rank1
    
    subroutine loadDataset_rank2(dataarray, datamask, tensor_id, loadstart)
      ! Load a chunk of data for a 2-rank tensor, beginning at loadstart
      implicit none
      real, dimension(:,:,:,:,:,:), intent(inout) :: dataarray
      logical, dimension(3,3), intent(in) :: datamask
      integer, intent(in) :: tensor_id
      integer, intent(in) :: loadstart
      integer :: ndims
      integer(HSIZE_T) :: mask_i, mask_j
      ndims = tensor_ndims(tensor_id)
      tensor_offsets(tensor_id,4) = loadstart
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
      end do; end do
      ! Read data into memory
      call H5Dread_F(tensor_id_D(tensor_id), hdf_memtype, dataarray, &
                     tensor_dims(tensor_id,1:ndims), hdferr, &
                     tensor_id_memS(tensor_id), tensor_id_S(tensor_id))
      dataarray = tensor_scales(tensor_id) * dataarray
      tensor_maxvals(tensor_id) = maxval(dataarray)
      tensor_minvals(tensor_id) = minval(dataarray)
    end subroutine loadDataset_rank2

    subroutine loadDataset_rank3(dataarray, datamask, tensor_id, loadstart)
      ! Load a chunk of data for a 3-rank tensor, beginning at loadstart
      implicit none
      real, dimension(:,:,:,:,:,:,:), intent(inout) :: dataarray
      logical, dimension(3,3,3), intent(in) :: datamask
      integer, intent(in) :: tensor_id
      integer, intent(in) :: loadstart
      integer :: ndims
      integer(HSIZE_T) :: mask_i, mask_j, mask_k
      ndims = tensor_ndims(tensor_id)
      tensor_offsets(tensor_id,4) = loadstart
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
      call H5Dread_F(tensor_id_D(tensor_id), hdf_memtype, dataarray, &
                     tensor_dims(tensor_id,1:ndims), hdferr, &
                     tensor_id_memS(tensor_id), tensor_id_S(tensor_id))
      dataarray = tensor_scales(tensor_id) * dataarray
      tensor_maxvals(tensor_id) = maxval(dataarray)
      tensor_minvals(tensor_id) = minval(dataarray)
    end subroutine loadDataset_rank3

    function emf_interpolate(dataarray) result(interp_data)
      real, intent(in), dimension(nx,dataload_len) :: dataarray
      real, dimension(nx) :: interp_data

      interp_data=dataarray(:,1)

    end function emf_interpolate

    subroutine setParameterDefaults
      implicit none
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
    end subroutine setParameterDefaults

    subroutine parseParameters
      implicit none
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
        lalpha     = .true.
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
        lbeta     = .true.
        lbeta_arr = .true.
      end if
!
! Load boolean array for gamma
!
      if (any(lgamma_c)) then
        lgamma_arr  = lgamma_c
      else if (lgamma) then
        lgamma      = .true.
        lgamma_arr  = .true.
      end if
!
! Load boolean array for delta
!
      if (any(ldelta_c)) then
        ldelta_arr  = ldelta_c
      else if (ldelta) then
        ldelta      = .true.
        ldelta_arr  = .true.
      end if
!
! Load boolean array for kappa
! TODO: implement kappa components
!
      if (lkappa) then
        lkappa_arr = .true.
      else
        lkappa_arr = .false.
      end if
!
! Load boolean array for acoef
! TODO: implement acoef components
!
      if (lacoef) then
        lacoef_arr = .true.
      else
        lacoef_arr = .false.
      end if
!
! Load boolean array for bcoef
! TODO: implement bcoef components
!
      if (lbcoef) then
        lbcoef_arr = .true.
      else
        lbcoef_arr = .false.
      end if
!
! Load boolean array for utensor
!
      if (any(lutensor_c)) then
        lutensor_arr  = lutensor_c
      else if (lutensor) then
        lutensor      = .true.
        lutensor_arr  = .true.
      end if
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
!****************************************************************************
  subroutine initialize_mult_special

  endsubroutine initialize_mult_special

!***********************************************************************
  subroutine finalize_mult_special


  endsubroutine finalize_mult_special
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!  Some precalculated pencils of data are passed in for efficiency,
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition), intent(in) :: bc
!
    endsubroutine special_boundconds
!*********************************************************************** 
  
endmodule Special
