! $Id$
!
! MODULE_DOC: This module contains GPU related types and functions to be used with the ASTAROTH nucleus.
!
! CPARAM logical, parameter :: lgpu = .true.
!
!**************************************************************************
!
module GPU
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Mpicomm, only: stop_it
!$ use, intrinsic :: iso_c_binding

  implicit none

  external initialize_gpu_c
  external finalize_gpu_c
  external rhs_gpu_c
  external copy_farray_c
  external test_rhs_c

!$  interface
!$    subroutine random_initial_condition() bind(C)
!$    endsubroutine random_initial_condition
!$  end interface



  include 'gpu.h'

  integer(KIND=ikind8) :: pFarr_GPU_in, pFarr_GPU_out
  
contains

!***********************************************************************
    subroutine initialize_GPU
!
      character(LEN=512) :: str
!
      str=''
      if (lanelastic) str=trim(str)//', '//'anelastic'
      if (lboussinesq) str=trim(str)//', '//'boussinesq'
      if (ltemperature) str=trim(str)//', '//'temperature'
      if (lshock) str=trim(str)//', '//'shock'
      if (lheatflux) str=trim(str)//', '//'heatflux'
      if (lhyperresistivity_strict) str=trim(str)//', '//'hyperresi_strict'
      if (lhyperviscosity_strict) str=trim(str)//', '//'hypervisc_strict'
      if (lADI) str=trim(str)//', '//'implicit_physics'
      if (llorenz_gauge) str=trim(str)//', '//'lorenz_gauge'
      if (ldustvelocity) str=trim(str)//', '//'dustvelocity'
      if (ldustdensity) str=trim(str)//', '//'dustdensity'
      if (ltestscalar) str=trim(str)//', '//'testscalar'
      if (ltestfield) str=trim(str)//', '//'testfield'
      if (ltestflow) str=trim(str)//', '//'testflow'
      if (linterstellar) str=trim(str)//', '//'interstellar'
      if (lcosmicray) str=trim(str)//', '//'cosmicray'
      if (lcosmicrayflux) str=trim(str)//', '//'cosmicrayflux'
      if (lshear) str=trim(str)//', '//'shear'
      if (lpscalar) str=trim(str)//', '//'pscalar'
      if (lascalar) str=trim(str)//', '//'ascalar'
      if (lradiation) str=trim(str)//', '//'radiation'
      if (lchemistry) str=trim(str)//', '//'chemistry'
      if (lchiral) str=trim(str)//', '//'chiral'
      if (ldetonate) str=trim(str)//', '//'detonate'
      if (lneutralvelocity) str=trim(str)//', '//'neutralvelocity'
      if (lneutraldensity) str=trim(str)//', '//'neutraldensity'
      if (lopacity) str=trim(str)//', '//'opacity'
      if (lpolymer) str=trim(str)//', '//'polymer'
      if (lpointmasses) str=trim(str)//', '//'pointmasses'
      if (lpoisson) str=trim(str)//', '//'poisson'
      if (lselfgravity) str=trim(str)//', '//'selfgravity'
      if (lsolid_cells) str=trim(str)//', '//'solid_cells'
      if (lspecial) str=trim(str)//', '//'special'
      if (lpower_spectrum) str=trim(str)//', '//'power_spectrum'
      if (lparticles) str=trim(str)//', '//'particles'

      if (str/='') call stop_it('No GPU implementation for module(s) "'//trim(str(3:))//'"')
!
      call initialize_gpu_c(pFarr_GPU_in,pFarr_GPU_out)
!print'(a,1x,Z0,1x,Z0)', 'pFarr_GPU_in,pFarr_GPU_out=', pFarr_GPU_in,pFarr_GPU_out
    endsubroutine initialize_GPU
!**************************************************************************
    subroutine gpu_init
!
      call init_gpu_c
!
    endsubroutine gpu_init
!**************************************************************************
    subroutine register_GPU(f)
!
      real, dimension(:,:,:,:), intent(IN) :: f

      call register_gpu_c(f)
!
    endsubroutine register_GPU
!**************************************************************************
    subroutine finalize_GPU
!
      call finalize_gpu_c
!
    endsubroutine finalize_GPU
!**************************************************************************
    subroutine rhs_GPU(f,isubstep,early_finalize)
!
      use General, only: notanumber

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f
      integer,                            intent(IN)    :: isubstep
      logical,                            intent(IN)    :: early_finalize
!
      integer :: ll, mm, nn
      logical, save :: lvery_first=.true.

      call rhs_gpu_c(isubstep,lvery_first,early_finalize)
!
      lvery_first=.false.
    endsubroutine rhs_GPU
!**************************************************************************
    subroutine copy_farray_from_GPU(f)

      real, dimension (mx,my,mz,mfarray), intent(OUT) :: f

      call copy_farray_c(f(1,1,1,iux),f(1,1,1,iuy),f(1,1,1,iuz),f(1,1,1,ilnrho))

    endsubroutine copy_farray_from_GPU
!**************************************************************************
 subroutine test_rhs_gpu(f,df,p,mass_per_proc,early_finalize,cpu_version)
!  Used to test the CPU rhs vs the DSL code
!
!  13-nov-23/TP: Written
!
      use MPIcomm
      use Boundcond
!$    use ISO_fortran_env, only: stdout => output_unit
!$    use, intrinsic :: iso_c_binding
      real, dimension (mx,my,mz,mfarray) :: f,f_copy,f_copy_2
      real, dimension (mx,my,mz,mfarray) :: df,df_copy,ds
      type (pencil_case) :: p,p_copy
      real, dimension(1), intent(inout) :: mass_per_proc
      logical ,intent(in) :: early_finalize
      integer :: i,j,k,n
      logical :: passed
      real, parameter :: dt = 0.001
      real, parameter, dimension(3) :: alpha = (/0.0, -(5.0/9.0), -(153.0/128.0)/)
      real, parameter, dimension(3) :: beta = (/ 1. / 3., 15./ 16., 8. / 15. /)
      integer, parameter :: num_of_steps = 1
      interface
          subroutine cpu_version(f,df,p,mass_per_proc,early_finalize)
              import mx
              import my
              import mz
              import mfarray
              import pencil_case
              real, dimension (mx,my,mz,mfarray) :: f
              real, dimension (mx,my,mz,mfarray) :: df
              type (pencil_case) :: p
              real, dimension(1), intent(inout) :: mass_per_proc
              logical ,intent(in) :: early_finalize

              intent(inout) :: f
              intent(inout) :: p
              intent(out) :: df
          endsubroutine cpu_version
        endinterface
      !TP: uncomment if want to test from random initial condition
      ! call random_initial_condition()
      ! call copy_farray_from_GPU(f)
      df_copy = df
      p_copy = p
      f_copy = f
      f_copy_2 = f
      do n=1,num_of_steps
        print*,"GPU rhs test:\tcpu step: ",n
        ds = 0.0
        do i=1,3
          call boundconds_x(f_copy)
          call initiate_isendrcv_bdry(f_copy)
          call finalize_isendrcv_bdry(f_copy)
          call boundconds_y(f_copy)
          call boundconds_z(f_copy)
          df_copy = 0.0
          call cpu_version(f_copy,df_copy,p,mass_per_proc,early_finalize)
          ds = alpha(i)*ds + df_copy*dt
          f_copy = f_copy + beta(i)*ds
        enddo
      enddo

    call test_rhs_c(f_copy_2,f_copy);
    call die_gracefully
  end subroutine  test_rhs_gpu
!**************************************************************************
endmodule GPU
