! $Id$
!
!  This module provide a way for users to specify custom
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc.
!
!  The module provides a set of standard hooks into the Pencil-Code and
!  currently allows the following customizations:
!
!  Description                                     | Relevant function call
!  ---------------------------------------------------------------------------
!  Special variable registration                   | register_special
!    (pre parameter read)                          |
!  Special variable initialization                 | initialize_special
!    (post parameter read)                         |
!  Special variable finalization                   | finalize_special
!    (deallocation, etc.)                          |
!                                                  |
!  Special initial condition                       | init_special
!   this is called last so may be used to modify   |
!   the mvar variables declared by this module     |
!   or optionally modify any of the other f array  |
!   variables.  The latter, however, should be     |
!   avoided where ever possible.                   |
!                                                  |
!  Special term in the mass (density) equation     | special_calc_density
!  Special term in the momentum (hydro) equation   | special_calc_hydro
!  Special term in the energy equation             | special_calc_energy
!  Special term in the induction (magnetic)        | special_calc_magnetic
!     equation                                     |
!                                                  |
!  Special equation                                | dspecial_dt
!    NOT IMPLEMENTED FULLY YET - HOOKS NOT PLACED INTO THE PENCIL-CODE
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED alpha_emf(3,3); beta_emf(3,3); gamma_emf(3)
! PENCILS PROVIDED delta_emf(3); kappa_emf(3,3,3)
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!   lspecial = .true.
! to enable use of special hooks.
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code SVN repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module
! selections to say something like:
!
!   SPECIAL=special/geo_kws
!
! Where geo_kws it replaced by the filename of your new module
! upto and not including the .f90
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
  use Mpicomm, only: mpibarrier
  use Sub, only: numeric_precision
  use HDF5
!
  implicit none
!
  include 'mpif.h'
  include '../special.h'

  integer :: hdferr, loaderr
  integer(HID_T) :: hdf_memtype, &
                    hdf_emftensors_file,  &
                    hdf_emftensors_plist, &
                    hdf_emftensors_group
  !integer(HID_T), dimension(:), allocatable :: hdf_fieldgroups
  !integer(HID_T), dimension(:,:), allocatable :: hdf_spaces, phdf_chunkspaces, phdf_datasets, phdf_plist_ids
  integer(HSIZE_T),dimension(4) :: hdf_stride, hdf_count, hdf_block, hdf_offsets

  integer(HSIZE_T),dimension(:,:),allocatable :: phdf_dims, phdf_chunkdims

  !integer :: i_alpha, i_beta, i_gamma, i_delta, i_kappa
  
  !integer(HID_T) :: alpha_id_G, alpha_id_D, alpha_id_S, &
  !                  beta_id_G, beta_id_D, beta_id_S,    &
  !                  gamma_id_G, gamma_id_D, gamma_id_S, &
  !                  delta_id_G, delta_id_D, delta_id_S, &
  !                  kappa_id_G, kappa_id_D, kappa_id_S
  integer(HID_T),   dimension(5) :: tensor_id_G, tensor_id_D, &
                                    tensor_id_S, tensor_id_memS
  integer(HSIZE_T), dimension(5,7) :: tensor_dims
  integer, dimension(5),parameter :: tensor_ndims = (/ 6, 6, 5, 5, 7 /) 
  integer,parameter :: alpha_id=1, beta_id=2, &
                       gamma_id=3, delta_id=4, &
                       kappa_id=5
  !integer(HSIZE_T) :: alpha_dims(6), &
  !                    beta_dims(6),  &
  !                    gamma_dims(5), &
  !                    delta_dims(5), &
  !                    kappa_dims(7)
  
  integer(HSIZE_T) :: alpha_tlen

  character (len=10) :: interpname
  character (len=10),dimension(1),parameter :: interpnames = (/ 'none' /)

  character (len=10) :: dataname
                    
  logical :: hdf_exists

  integer :: dataload_len
  !real, dimension(3,3,  nx,ny,nz,dataload_len) :: alpha_data, beta_data
 ! real, dimension(3,    nx,ny,nz,dataload_len) :: gamma_data, delta_data
  !real, dimension(3,3,3,nx,ny,nz,dataload_len) :: kappa_data

  real, dimension(:,:,:,:,:,:)  , allocatable :: alpha_data, beta_data
  real, dimension(:,:,:,:,:)    , allocatable :: gamma_data, delta_data
  real, dimension(:,:,:,:,:,:,:), allocatable :: kappa_data
  integer :: t_index

  character (len=fnlen) :: hdf_emftensors_filename

  logical :: lalpha, lalpha_exists, lalpha_open, &
             lbeta,  lbeta_exists,  lbeta_open, &
             lgamma, lgamma_exists, lgamma_open, &
             ldelta, ldelta_exists, ldelta_open, &
             lkappa, lkappa_exists, lkappa_open
  
  namelist /special_init_pars/ &
       lalpha, lbeta, lgamma, lkappa, dataname, interpname
!
! Declare index of new variables in f array (if any).
!
!!   integer :: ispecial=0
!!   integer :: ispecaux=0
!
!! Diagnostic variables (needs to be consistent with reset list below).
!
!!   integer :: idiag_POSSIBLEDIAGNOSTIC=0
!
  interface loadDataset
    module procedure loadDataset_rank1
  end interface

  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
           "$Id$")
      interpname = 'none'
      lalpha=.false.
      lbeta =.false.
      lgamma=.false.
      lkappa=.false.
      dataname = 'data'

      if (trim(interpname) == 'none') then
        dataload_len = 100
      else
        call fatal_error('initialize_special','File '//trim(hdf_emftensors_filename)//' does not exist!')
      end if
!
!!      call farray_register_pde('special',ispecial)
!!      call farray_register_auxiliary('specaux',ispecaux)
!!      call farray_register_auxiliary('specaux',ispecaux,communicated=.true.)
!
    endsubroutine register_special
!***********************************************************************
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
!
      call keep_compiler_quiet(f)

      print *,'Initialize special <----------------'
  
      ! Allocate data arrays
      allocate(alpha_data(3,3,  nx,ny,nz,dataload_len))
      allocate( beta_data(3,3,  nx,ny,nz,dataload_len))
      allocate(gamma_data(3,    nx,ny,nz,dataload_len))
      allocate(delta_data(3,    nx,ny,nz,dataload_len))
      allocate(kappa_data(3,3,3,nx,ny,nz,dataload_len))

      hdf_emftensors_plist = -1
      hdf_emftensors_file  = -1

      call H5open_F(hdferr)

      hdf_memtype=-1
      if (numeric_precision() == 'S') then
        write(*,*) 'Register special: loading data as single precision'
        hdf_memtype = H5T_NATIVE_REAL
        !hdf_memtype = H5T_IEEE_F64LE
      else
        write(*,*) 'Register special: loading data as double precision'
        hdf_memtype = H5T_NATIVE_DOUBLE
        !hdf_memtype = H5T_IEEE_F64LE
      end if
      
      print *,'Setting parallel HDF5 IO for data file reading'
      call H5Pcreate_F(H5P_FILE_ACCESS_F, hdf_emftensors_plist, hdferr)
      call H5Pset_fapl_mpio_F(hdf_emftensors_plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdferr)

      print *, 'Opening emftensors.h5 and loading relevant fields of  it into memory...'

      hdf_emftensors_filename = trim(datadir_snap)//'/emftensors.h5'

      inquire(file=hdf_emftensors_filename, exist=hdf_exists)

      if (.not. hdf_exists) then
        call H5Pclose_F(hdf_emftensors_plist, hdferr)
        call H5close_F(hdferr)
        call fatal_error('initialize_special','File '//trim(hdf_emftensors_filename)//' does not exist!')
      end if

      call H5Fopen_F(hdf_emftensors_filename, H5F_ACC_RDONLY_F, hdf_emftensors_file, hdferr, access_prp = hdf_emftensors_plist)

      call H5Lexists_F(hdf_emftensors_file,'/zaver/', hdf_exists, hdferr)

      if (.not. hdf_exists) then
        call H5Fclose_F(hdf_emftensors_file, hdferr)
        call H5Pclose_F(hdf_emftensors_plist, hdferr)
        call H5close_F(hdferr)
        call fatal_error('initialize_special','group /zaver/ does not exist!')
      end if

      call H5Gopen_F(hdf_emftensors_file, 'zaver', hdf_emftensors_group, hdferr)

      ! Open datasets

      ! alpha
      lalpha_open = .false.
      call openDataset('alpha', alpha_id, &
                       lalpha_exists, & 
                       loaderr)
      if (loaderr == 0) then
        lalpha_open = .true.
        write(*,*) 'initialize_special: alpha found. Using dataset /zaver/alpha/'//trim(dataname)
      else if (lalpha) then
        call fatal_error('initialize_special','alpha not found, but required. Tried to use dataset /zaver/alpha/'//trim(dataname))
      end if

      ! beta
      lbeta_open = .false.
      call openDataset('beta', beta_id, &
                       lbeta_exists, & 
                       loaderr)
      if (loaderr == 0) then
        lbeta_open = .true.
        write(*,*) 'initialize_special: beta found. Using dataset /zaver/beta/'//trim(dataname)
      else if (lbeta) then
        call fatal_error('initialize_special','beta not found, but required. Tried to use dataset /zaver/beta/'//trim(dataname))
      end if

      ! gamma
      lgamma_open = .false.
      call openDataset('gamma', gamma_id, &
                       lgamma_exists, & 
                       loaderr)
      if (loaderr == 0) then
        lgamma_open = .true.
        write(*,*) 'initialize_special: gamma found. Using dataset /zaver/gamma/'//trim(dataname)
        write(*,*) 'emftensor: gamma found'
        write(*,*) '  using dataset /zaver/gamma/'//trim(dataname)
      else if (lgamma) then
        call fatal_error('initialize_special','gamma not found, but required. Tried to use dataset /zaver/gamma/'//trim(dataname))
      end if

      ! delta
      ldelta_open = .false.
      call openDataset('delta', delta_id, &
                       ldelta_exists, & 
                       loaderr)
      if (loaderr == 0) then
        ldelta_open = .true.
        write(*,*) 'initialize_special: delta found. Using dataset /zaver/delta/'//trim(dataname)
      else if (ldelta) then
        call fatal_error('initialize_special','delta not found, but required. Tried to use dataset /zaver/delta/'//trim(dataname))
      end if
      
      ! kappa
      lkappa_open = .false.
      call openDataset('kappa', kappa_id, &
                       lkappa_exists, & 
                       loaderr)
      if (loaderr == 0) then
        lkappa_open = .true.
        write(*,*) 'initialize_special: kappa found. Using dataset /zaver/kappa/'//trim(dataname)
      else if (lkappa) then
        call fatal_error('initialize_special','kappa not found, but required. Tried to use dataset /zaver/kappa/'//trim(dataname))
      end if
      
      ! Load initial dataset values

      call loadDataset(gamma_data, gamma_id, 0)
      write(*,*) sum(gamma_data), maxval(gamma_data), minval(gamma_data)
      write(*,*) sum(gamma_data(1,:,:,:,2)), maxval(gamma_data(1,:,:,:,2)), minval(gamma_data(1,:,:,:,2))
      write(*,*) gamma_data(1,10,10,1,10)
      write(*,*) gamma_data(1,11,10,1,10)
      write(*,*) gamma_data(1,10,11,1,10)
      write(*,*) gamma_data(1,10,10,1,11)
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

      ! Deallocate data
      deallocate(alpha_data)
      deallocate( beta_data)
      deallocate(gamma_data)
      deallocate(delta_data)
      deallocate(kappa_data)

      print *,'Closing emftensors.h5'
      
      if (lalpha_open) then
        call closeDataset(alpha_id)
      end if
      if (lbeta_open) then
        call closeDataset(beta_id)
      end if
      if (lgamma_open) then
        call closeDataset(gamma_id)
      end if
      if (ldelta_open) then
        call closeDataset(delta_id)
      end if
      if (lkappa_open) then
        call closeDataset(kappa_id)
      end if

      call H5Gclose_F(hdf_emftensors_group, hdferr)
      call H5Fclose_F(hdf_emftensors_file, hdferr)
      call H5Pclose_F(hdf_emftensors_plist, hdferr)
      call H5close_F(hdferr)
      call mpibarrier()
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
      lpenc_requested(i_uu)=.true.
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_jj)=.true.
      lpenc_requested(i_alpha_emf)=.true.
      lpenc_requested(i_beta_emf)=.true.
      lpenc_requested(i_gamma_emf)=.true.
      lpenc_requested(i_delta_emf)=.true.
      lpenc_requested(i_kappa_emf)=.true.

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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
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
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
      print *,'Read special init pars <--------'
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
      iostat = 0
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
      
      print *,'Read special run params <-----------'
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
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
!!      integer :: iname
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
!!        idiag_SPECIAL_DIAGNOSTIC=0
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
    subroutine calc_lspecial_pars(f)
!
!  Dummy routine.
!
!  15-jan-08/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lspecial_pars
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
      type (pencil_case), intent(in) :: p
      
      !print *,'Calc special magnetic <-----------'
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
    subroutine special_calc_particles_nbody(fsp)
!
!  Called before the loop, in case some massive particles value
!  is needed for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (:,:), intent(in) :: fsp
!
      call keep_compiler_quiet(fsp)
!
    endsubroutine special_calc_particles_nbody
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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc)
!
    endsubroutine special_boundconds
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
subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_after_boundary

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
                           tensor_id,     &
                           dataset_exists, &
                           loaderr)

      ! Load a certain timestep from input array
      character(len=*), intent(in)    :: datagroup
      integer, intent(in)             :: tensor_id
      logical, intent(inout)          :: dataset_exists
      integer, intent(out) :: loaderr
!
      integer(HSIZE_T), dimension(10) :: maxdimsizes
      integer :: ndims
      
      dataset_exists = .true.

      ! Check that datagroup e.g. /zaver/alpha exists
      call H5Lexists_F(hdf_emftensors_group, datagroup, dataset_exists, hdferr)
      if (hdferr /= 0) then
        loaderr = 1
        dataset_exists = .false.
        return
      end if
      ! Open datagroup
      call H5Gopen_F(hdf_emftensors_group, datagroup, tensor_id_G(tensor_id), hdferr)
      if (hdferr /= 0) then
        loaderr = 2
        return
      end if

      ! Check that dataset e.g. /zaver/alpha/data exists
      call H5Lexists_F(tensor_id_G(tensor_id), dataname, dataset_exists, hdferr)
      if (hdferr /= 0) then
        loaderr = 3
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        dataset_exists = .false.
        return
      end if
      ! Open dataset
      call H5Dopen_F(tensor_id_G(tensor_id), dataname, tensor_id_D(tensor_id), loaderr)
      if (hdferr /= 0) then
        loaderr = 4
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        return
      end if
      ! Get dataspace
      call H5Dget_space_F(tensor_id_D(tensor_id), tensor_id_S(tensor_id), hdferr)
      if (hdferr /= 0) then
        loaderr = 5
        call H5Dclose_F(tensor_id_D(tensor_id), hdferr)
        call H5Gclose_F(tensor_id_G(tensor_id), hdferr)
        return
      end if
      !call H5Sget_simple_extent_ndims_F(dataset_id_S, ndims, hdferr)
      ndims = tensor_ndims(tensor_id)
      call H5Sget_simple_extent_dims_F(tensor_id_S(tensor_id), &
                                       tensor_dims(tensor_id,1:ndims), &
                                       maxdimsizes, hdferr)
      print *,tensor_dims(tensor_id,1:tensor_ndims(tensor_id))
      print *,maxdimsizes(1:tensor_ndims(tensor_id))
      loaderr = 0

    end subroutine openDataset

    subroutine closeDataset(tensor_id)
      implicit none
      integer :: tensor_id
      ! Load a certain timestep from input array
      call H5Sclose_F(tensor_id_S(tensor_id), hdferr)
      call H5Dclose_F(tensor_id_D(tensor_id), hdferr)
      call H5Gclose_F(tensor_id_G(tensor_id), hdferr)

    end subroutine closeDataset

    subroutine loadDataset_rank1(dataarray, tensor_id, loadstart)
      implicit none
      real, dimension(:,:,:,:,:), intent(inout) :: dataarray
      integer, intent(in) :: tensor_id
      integer, intent(in) :: loadstart
      integer(HSIZE_T), dimension(5) :: dataarray_dims
      integer(HSIZE_T), dimension(5) :: hdf_offsets, hdf_stride, &
                                        hdf_block,   hdf_count, &
                                        hdf_memoffsets, hdf_memcount
      integer(HID_T)  :: hdf_memspace
!
!      call H5Dget_space_F(dataset_id_D, hdf_dataspace, hdferr)
      dataarray_dims = tensor_dims(tensor_id,1:5)
      hdf_offsets  = [ 0, ipx*nx, ipy*ny, ipz*nz, loadstart ]
      hdf_stride   = [ 1, 1, 1, 1, 1 ]
      hdf_block    = [ 1, 1, 1, 1, 1 ]
      hdf_count    = [ 3, nx, ny, nz, dataload_len ]
!
      call H5Sselect_hyperslab_F(tensor_id_S(tensor_id), H5S_SELECT_SET_F, &
                                 hdf_offsets, hdf_count, &
                                 hdferr, hdf_stride, hdf_block)

      hdf_memoffsets = [ 0 , 0 , 0 , 0 , 0 ]
      hdf_memcount  = [ 3, nx, ny , nz, dataload_len ]
      call H5Screate_simple_F(5, hdf_memcount, hdf_memspace, hdferr)
      call H5Sselect_hyperslab_F(hdf_memspace, H5S_SELECT_SET_F, &
                                 hdf_memoffsets, hdf_memcount, &
                                 hdferr)

      call H5Dread_F(tensor_id_D(tensor_id), hdf_memtype, dataarray, dataarray_dims, hdferr, &
                     hdf_memspace, tensor_id_S(tensor_id))
      call H5Sclose_F(hdf_memspace, hdferr)
      !call H5Sclose_F(hdf_dataspace, hdferr)

    end subroutine loadDataset_rank1

endmodule Special
