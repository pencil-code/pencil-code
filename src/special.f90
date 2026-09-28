! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
!***************************************************************
!
  module Special

    use Cparam
    use Cdata, only: lroot, n_special_modules, special_modules, n_odevars
    use, intrinsic :: iso_c_binding, only: c_funptr, c_f_procpointer

    implicit none

    include 'special.h'
!
!  Interfaces of the hooks of the special modules; the check that the modules
!  implement them exactly is in check_special_hook_interfaces (nospecial.f90).
!
    include 'special_interfaces.inc'

    integer(KIND=ikind8), external :: dlopen_c, dlsym_c
    external dlclose_c
!
    integer, parameter :: I_REGISTER_SPECIAL=1,  &
                          I_REGISTER_PARTICLES_SPECIAL=2,  &
                          I_INITIALIZE_SPECIAL=3,  &
                          I_FINALIZE_SPECIAL=4,  &
                          I_READ_SPECIAL_INIT_PARS=5,  &
                          I_WRITE_SPECIAL_INIT_PARS=6,  &
                          I_READ_SPECIAL_RUN_PARS=7,  &
                          I_WRITE_SPECIAL_RUN_PARS=8,  &
                          I_RPRINT_SPECIAL=9,  &
                          I_GET_SLICES_SPECIAL=10,  &
                          I_INIT_SPECIAL=11,  &
                          I_DSPECIAL_DT=12,  &
                          I_DSPECIAL_DT_ODE=32,  &
                          I_CALC_PENCILS_SPECIAL=13,  &
                          I_PENCIL_CRITERIA_SPECIAL=14,  &
                          I_PENCIL_INTERDEP_SPECIAL=15,  &
                          I_SPECIAL_CALC_HYDRO=16, &
                          I_SPECIAL_CALC_DENSITY=17, &
                          I_SPECIAL_CALC_DUSTDENSITY=18, &
                          I_SPECIAL_CALC_ENERGY=19, &
                          I_SPECIAL_CALC_MAGNETIC=20, &
                          I_SPECIAL_CALC_PSCALAR=21, &
                          I_SPECIAL_CALC_PARTICLES=22, &
                          I_SPECIAL_CALC_CHEMISTRY=23, &
                          I_SPECIAL_BOUNDCONDS=24,  &
                          I_SPECIAL_BEFORE_BOUNDARY=25,  &
                          I_SPECIAL_PARTICLES_BFRE_BDARY=26,  &
                          I_SPECIAL_AFTER_BOUNDARY=27,  &
                          I_SPECIAL_AFTER_TIMESTEP=28,  &
                          I_SET_INIT_PARAMETERS=29, &
                          I_SPECIAL_CALC_SPECTRA=30, &
                          I_SPECIAL_CALC_SPECTRA_BYTE=31, &
                          I_INPUT_PERSIST_SPECIAL=33, &
                          I_INPUT_PERSIST_SPECIAL_ID=34, &
                          I_OUTPUT_PERSISTENT_SPECIAL=35, &
                          I_SPECIAL_PARTICLES_AFTER_DTSUB=36, &
                          I_CALC_DIAGNOSTICS_SPECIAL=37, &
                          I_CALC_ODE_DIAGNOSTICS_SPECIAL=38, &
                          I_PREP_RHS_SPECIAL=39, &
                          I_LOAD_VARIABLES_TO_GPU_SPECIAL=40, &
                          I_SPECIAL_BEFORE_BOUNDARY_DIAGNOSTICS=41
    
    integer, parameter :: n_subroutines=41
!
    character(LEN=256) :: special_modules_list = ''
    character(LEN=128), dimension(n_subroutines) :: special_subroutines=(/ &
                           'register_special                    ', &
                           'register_particles_special          ', &
                           'initialize_special                  ', &
                           'finalize_special                    ', &
                           'read_special_init_pars              ', &
                           'write_special_init_pars             ', &
                           'read_special_run_pars               ', &
                           'write_special_run_pars              ', &
                           'rprint_special                      ', &
                           'get_slices_special                  ', &
                           'init_special                        ', &
                           'dspecial_dt                         ', &
                           'calc_pencils_special                ', &
                           'pencil_criteria_special             ', &
                           'pencil_interdep_special             ', &
                           'special_calc_hydro                  ', &
                           'special_calc_density                ', &
                           'special_calc_dustdensity            ', &
                           'special_calc_energy                 ', &
                           'special_calc_magnetic               ', &
                           'special_calc_pscalar                ', &
                           'special_calc_particles              ', &
                           'special_calc_chemistry              ', &
                           'special_boundconds                  ', &
                           'special_before_boundary             ', &
                           'special_particles_bfre_bdary        ', &
                           'special_after_boundary              ', &
                           'special_after_timestep              ', &
                           'set_init_parameters                 ', &
                           'special_calc_spectra                ', &
                           'special_calc_spectra_byte           ', &
                           'dspecial_dt_ode                     ', &
                           'input_persist_special               ', &
                           'input_persist_special_id            ', &
                           'output_persistent_special           ', &
                           'special_particles_after_dtsub       ', &
                           'calc_diagnostics_special            ', &
                           'calc_ode_diagnostics_special        ', &
                           'prep_rhs_special                    ', &
                           'load_variables_to_gpu_special       ', &
                           'special_before_boundary_diagnostics '  &
                   /)

    integer(KIND=ikind8) :: libhandle
    integer(KIND=ikind8), dimension(n_special_modules,n_subroutines) :: special_sub_handles
!
!  The hooks of each special module as procedure pointers with explicit
!  interfaces, so that all calls below are checked by the compiler.
!
    type special_hooks
      procedure(iface_register_special),               pointer, nopass :: register_special
      procedure(iface_register_particles_special),     pointer, nopass :: register_particles_special
      procedure(iface_initialize_special),             pointer, nopass :: initialize_special
      procedure(iface_finalize_special),               pointer, nopass :: finalize_special
      procedure(iface_read_special_pars),              pointer, nopass :: read_special_init_pars
      procedure(iface_write_special_pars),             pointer, nopass :: write_special_init_pars
      procedure(iface_read_special_pars),              pointer, nopass :: read_special_run_pars
      procedure(iface_write_special_pars),             pointer, nopass :: write_special_run_pars
      procedure(iface_rprint_special),                 pointer, nopass :: rprint_special
      procedure(iface_get_slices_special),             pointer, nopass :: get_slices_special
      procedure(iface_init_special),                   pointer, nopass :: init_special
      procedure(iface_dspecial_dt),                    pointer, nopass :: dspecial_dt
      procedure(iface_special_noargs),                 pointer, nopass :: dspecial_dt_ode
      procedure(iface_calc_pencils_special),           pointer, nopass :: calc_pencils_special
      procedure(iface_special_noargs),                 pointer, nopass :: pencil_criteria_special
      procedure(iface_pencil_interdep_special),        pointer, nopass :: pencil_interdep_special
      procedure(iface_special_calc_rhs),               pointer, nopass :: special_calc_hydro
      procedure(iface_special_calc_rhs),               pointer, nopass :: special_calc_density
      procedure(iface_special_calc_rhs),               pointer, nopass :: special_calc_dustdensity
      procedure(iface_special_calc_rhs),               pointer, nopass :: special_calc_energy
      procedure(iface_special_calc_rhs),               pointer, nopass :: special_calc_magnetic
      procedure(iface_special_calc_rhs),               pointer, nopass :: special_calc_pscalar
      procedure(iface_special_calc_particles),         pointer, nopass :: special_calc_particles
      procedure(iface_special_calc_rhs),               pointer, nopass :: special_calc_chemistry
      procedure(iface_special_boundconds),             pointer, nopass :: special_boundconds
      procedure(iface_special_boundary),               pointer, nopass :: special_before_boundary
      procedure(iface_special_particles_bfre_bdary),   pointer, nopass :: special_particles_bfre_bdary
      procedure(iface_special_boundary),               pointer, nopass :: special_after_boundary
      procedure(iface_special_after_timestep),         pointer, nopass :: special_after_timestep
      procedure(iface_set_init_parameters),            pointer, nopass :: set_init_parameters
      procedure(iface_special_calc_spectra),           pointer, nopass :: special_calc_spectra
      procedure(iface_special_noargs),                 pointer, nopass :: input_persist_special
      procedure(iface_input_persist_special_id),       pointer, nopass :: input_persist_special_id
      procedure(iface_output_persistent_special),      pointer, nopass :: output_persistent_special
      procedure(iface_special_particles_after_dtsub),  pointer, nopass :: special_particles_after_dtsub
      procedure(iface_calc_diagnostics_special),       pointer, nopass :: calc_diagnostics_special
      procedure(iface_special_ode),                    pointer, nopass :: calc_ode_diagnostics_special
      procedure(iface_special_ode),                    pointer, nopass :: prep_rhs_special
      procedure(iface_special_noargs),                 pointer, nopass :: load_variables_to_gpu_special
      procedure(iface_special_boundary),               pointer, nopass :: special_before_boundary_diagnostics
    endtype special_hooks
    type(special_hooks), dimension(n_special_modules) :: hooks
    character(LEN=128) :: specific_subroutine

    contains
!****************************************************************************
  subroutine initialize_mult_special

    use Cdata, only: lroot, lreloading
    use General, only: parser, safe_string_replace
    use Messages, only: fatal_error
    use Mpicomm, only: mpibcast
    use Syscalls, only: extract_str, get_env_var

    integer, parameter :: RTLD_LAZY=0, RTLD_NOW=1

    character(LEN=128) :: line,parstr
    integer :: i,j,ipos,ind
    integer(KIND=ikind8) :: sub_handle

    if (lreloading) return

    call get_env_var("PC_MODULES_LIST", special_modules_list)
    ind = index(special_modules_list,'#')
    if (ind>0) special_modules_list(ind:)=''  ! remove trailing comment
    if (n_special_modules/=parser(trim(special_modules_list),special_modules,' ')) &
      call fatal_error('initialize_mult_special','number of names in $PC_MODULES_LIST /= n_special_modules')

    !Remove trailing newlines
    do i=1,n_special_modules
        if (len_trim(special_modules(i)) > 0 .and. &
        special_modules(i)(len_trim(special_modules(i)):len_trim(special_modules(i))) == CHAR(10)) then
                 special_modules(i) = special_modules(i)(:len_trim(special_modules(i))-1)
        end if
    enddo
!if (lroot) print*, 'special_modules_list=', trim(special_modules_list)//'<<<'

!
    libhandle=dlopen_c('src/special.so'//char(0),RTLD_NOW)
    if (libhandle==0) &
      call fatal_error('initialize_mult_special','library src/special.so could not be opened')

    if (lroot) then
      call extract_str("nm src/special.so|grep calc_pencils_special|grep "//trim(special_modules(1))// &
                       "|grep -i ' T '",parstr)
      ipos=index(trim(parstr),' ',back=.true.)
      line=parstr(ipos+1:)
    endif
    call mpibcast(line)

    do i=1,n_special_modules
      do j=1,n_subroutines
        specific_subroutine = trim(line)
        call safe_string_replace(specific_subroutine,trim(special_modules(1)),trim(special_modules(i)))
        call safe_string_replace(specific_subroutine,'calc_pencils_special',trim(special_subroutines(j)))
        if (index(specific_subroutine,'___')==1) specific_subroutine=specific_subroutine(2:)    ! MR: a hack needed on MacOS
        sub_handle=dlsym_c(libhandle,trim(specific_subroutine)//char(0))
!print*, 'sub_handle=', sub_handle
        if (sub_handle==0) &
          call fatal_error('initialize_mult_special','Error for symbol '// &
          trim(specific_subroutine)//' in module '//trim(special_modules(i))) 
        special_sub_handles(i,j) = sub_handle
      enddo
    enddo

!
!  Make the looked-up addresses callable through the typed procedure pointers.
!
    do i=1,n_special_modules
      call c_f_procpointer(hook(i,I_REGISTER_SPECIAL),             hooks(i)%register_special)
      call c_f_procpointer(hook(i,I_REGISTER_PARTICLES_SPECIAL),   hooks(i)%register_particles_special)
      call c_f_procpointer(hook(i,I_INITIALIZE_SPECIAL),           hooks(i)%initialize_special)
      call c_f_procpointer(hook(i,I_FINALIZE_SPECIAL),             hooks(i)%finalize_special)
      call c_f_procpointer(hook(i,I_READ_SPECIAL_INIT_PARS),       hooks(i)%read_special_init_pars)
      call c_f_procpointer(hook(i,I_WRITE_SPECIAL_INIT_PARS),      hooks(i)%write_special_init_pars)
      call c_f_procpointer(hook(i,I_READ_SPECIAL_RUN_PARS),        hooks(i)%read_special_run_pars)
      call c_f_procpointer(hook(i,I_WRITE_SPECIAL_RUN_PARS),       hooks(i)%write_special_run_pars)
      call c_f_procpointer(hook(i,I_RPRINT_SPECIAL),               hooks(i)%rprint_special)
      call c_f_procpointer(hook(i,I_GET_SLICES_SPECIAL),           hooks(i)%get_slices_special)
      call c_f_procpointer(hook(i,I_INIT_SPECIAL),                 hooks(i)%init_special)
      call c_f_procpointer(hook(i,I_DSPECIAL_DT),                  hooks(i)%dspecial_dt)
      call c_f_procpointer(hook(i,I_DSPECIAL_DT_ODE),              hooks(i)%dspecial_dt_ode)
      call c_f_procpointer(hook(i,I_CALC_PENCILS_SPECIAL),         hooks(i)%calc_pencils_special)
      call c_f_procpointer(hook(i,I_PENCIL_CRITERIA_SPECIAL),      hooks(i)%pencil_criteria_special)
      call c_f_procpointer(hook(i,I_PENCIL_INTERDEP_SPECIAL),      hooks(i)%pencil_interdep_special)
      call c_f_procpointer(hook(i,I_SPECIAL_CALC_HYDRO),           hooks(i)%special_calc_hydro)
      call c_f_procpointer(hook(i,I_SPECIAL_CALC_DENSITY),         hooks(i)%special_calc_density)
      call c_f_procpointer(hook(i,I_SPECIAL_CALC_DUSTDENSITY),     hooks(i)%special_calc_dustdensity)
      call c_f_procpointer(hook(i,I_SPECIAL_CALC_ENERGY),          hooks(i)%special_calc_energy)
      call c_f_procpointer(hook(i,I_SPECIAL_CALC_MAGNETIC),        hooks(i)%special_calc_magnetic)
      call c_f_procpointer(hook(i,I_SPECIAL_CALC_PSCALAR),         hooks(i)%special_calc_pscalar)
      call c_f_procpointer(hook(i,I_SPECIAL_CALC_PARTICLES),       hooks(i)%special_calc_particles)
      call c_f_procpointer(hook(i,I_SPECIAL_CALC_CHEMISTRY),       hooks(i)%special_calc_chemistry)
      call c_f_procpointer(hook(i,I_SPECIAL_BOUNDCONDS),           hooks(i)%special_boundconds)
      call c_f_procpointer(hook(i,I_SPECIAL_BEFORE_BOUNDARY),      hooks(i)%special_before_boundary)
      call c_f_procpointer(hook(i,I_SPECIAL_PARTICLES_BFRE_BDARY), hooks(i)%special_particles_bfre_bdary)
      call c_f_procpointer(hook(i,I_SPECIAL_AFTER_BOUNDARY),       hooks(i)%special_after_boundary)
      call c_f_procpointer(hook(i,I_SPECIAL_AFTER_TIMESTEP),       hooks(i)%special_after_timestep)
      call c_f_procpointer(hook(i,I_SET_INIT_PARAMETERS),          hooks(i)%set_init_parameters)
      call c_f_procpointer(hook(i,I_SPECIAL_CALC_SPECTRA),         hooks(i)%special_calc_spectra)
      call c_f_procpointer(hook(i,I_INPUT_PERSIST_SPECIAL),        hooks(i)%input_persist_special)
      call c_f_procpointer(hook(i,I_INPUT_PERSIST_SPECIAL_ID),     hooks(i)%input_persist_special_id)
      call c_f_procpointer(hook(i,I_OUTPUT_PERSISTENT_SPECIAL),    hooks(i)%output_persistent_special)
      call c_f_procpointer(hook(i,I_SPECIAL_PARTICLES_AFTER_DTSUB),hooks(i)%special_particles_after_dtsub)
      call c_f_procpointer(hook(i,I_CALC_DIAGNOSTICS_SPECIAL),     hooks(i)%calc_diagnostics_special)
      call c_f_procpointer(hook(i,I_CALC_ODE_DIAGNOSTICS_SPECIAL), hooks(i)%calc_ode_diagnostics_special)
      call c_f_procpointer(hook(i,I_PREP_RHS_SPECIAL),             hooks(i)%prep_rhs_special)
      call c_f_procpointer(hook(i,I_LOAD_VARIABLES_TO_GPU_SPECIAL),hooks(i)%load_variables_to_gpu_special)
      call c_f_procpointer(hook(i,I_SPECIAL_BEFORE_BOUNDARY_DIAGNOSTICS), &
                                                                   hooks(i)%special_before_boundary_diagnostics)
    enddo
  endsubroutine initialize_mult_special
!***********************************************************************
  function hook(imod,isub)
!
!  Address of subroutine isub of special module imod, as looked up with dlsym.
!
    use, intrinsic :: iso_c_binding, only: c_null_funptr
!
    integer, intent(in) :: imod, isub
    type(c_funptr) :: hook
!
    hook=transfer(special_sub_handles(imod,isub),c_null_funptr)
!
  endfunction hook
!***********************************************************************
    subroutine register_special
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%register_special()
      enddo
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      use Cdata, only: special_module_index
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      integer :: i
!
      do i=1,n_special_modules
        special_module_index=i
        call hooks(i)%initialize_special(f)
      enddo
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine finalize_special(f)
!
!  Called right before exiting.
!
!  14-aug-2011/Bourdin.KIS: coded
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%finalize_special(f)
      enddo

      call dlclose_c(libhandle)
!
    endsubroutine finalize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%init_special(f)
      enddo
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%pencil_criteria_special()
      enddo

    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%pencil_interdep_special(lpencil_in)
      enddo
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
      real, contiguous, dimension(:,:,:,:) :: f
      type(pencil_case) :: p
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%calc_pencils_special(f,p)
      enddo
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt_ode
!
!  calculate right hand side of ONE OR MORE extra coupled ODEs
!
!  07-sep-23/MR: coded
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%dspecial_dt_ode()
      enddo
!
    endsubroutine dspecial_dt_ode
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  06-oct-03/tony: coded
!
      use Cdata, only: lspecial_substepped, lsubstepping_in_time
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(mx,my,mz,mvar) :: df
      type(pencil_case) :: p
!
      integer :: i
!
      do i=1,n_special_modules
        if (lsubstepping_in_time .eqv. lspecial_substepped(i)) then
          call hooks(i)%dspecial_dt(f,df,p)
        endif
      enddo
!
    endsubroutine dspecial_dt
!****************************************************************************
    subroutine register_particles_special(npvar)
!
!  Set up indices for particle variables in special modules.
!
!  4-jan-14/tony: coded
!
      integer :: npvar
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%register_particles_special(npvar)
      enddo
!
    endsubroutine register_particles_special
!***********************************************************************
    subroutine read_special_init_pars(iomsg)
!
      use File_io, only: parallel_rewind
!
      character(len=iomsglen), intent(out) :: iomsg
!
      integer :: i
      character(len=iomsglen) :: msg
!
      iomsg=''
      do i=1,n_special_modules
        call hooks(i)%read_special_init_pars(msg)
        if (msg/='') then
          iomsg=trim(iomsg)//new_line('a')//trim(special_modules(i))//': '//trim(msg)
          call parallel_rewind
        endif
      enddo
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%write_special_init_pars(unit)
      enddo
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iomsg)
!
      use File_io, only: parallel_rewind
!
      character(len=iomsglen), intent(out) :: iomsg
!
      integer :: i
      character(len=iomsglen) :: msg
!
      iomsg=''
      do i=1,n_special_modules
        call hooks(i)%read_special_run_pars(msg)
        if (msg/='') then
          iomsg=trim(iomsg)//new_line('a')//trim(special_modules(i))//': '//trim(msg)
          call parallel_rewind
        endif
      enddo
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%write_special_run_pars(unit)
      enddo
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%rprint_special(lreset,lwrite)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
!
!  26-jun-06/tony: dummy
!
      real, contiguous, dimension(:,:,:,:) :: f
      type(slice_data) :: slices
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%get_slices_special(f,slices)
      enddo
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  momentum equation.
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_calc_hydro(f,df,p)
      enddo
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  continuity equation.
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_calc_density(f,df,p)
      enddo
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_dustdensity(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  dust continuity equation.
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_calc_dustdensity(f,df,p)
      enddo
!
    endsubroutine special_calc_dustdensity
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  energy equation.
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_calc_energy(f,df,p)
      enddo
!
    endsubroutine special_calc_energy
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  induction equation.
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_calc_magnetic(f,df,p)
      enddo
!
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine special_calc_pscalar(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  passive scalar equation.
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_calc_pscalar(f,df,p)
      enddo
!
    endsubroutine special_calc_pscalar
!***********************************************************************
    subroutine special_calc_chemistry(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  chemistry equation.
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_calc_chemistry(f,df,p)
      enddo
!
    endsubroutine special_calc_chemistry
!***********************************************************************
    subroutine special_calc_particles(f,df,fp,dfp,ineargrid)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(:,:) :: fp, dfp
      integer, dimension(:,:) :: ineargrid
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_calc_particles(f,df,fp,dfp,ineargrid)
      enddo
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
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(:,:) :: fp
      integer, dimension(:,:) :: ineargrid
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_particles_bfre_bdary(f,fp,ineargrid)
      enddo
!
    endsubroutine special_particles_bfre_bdary
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_before_boundary(f)
      enddo
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_after_boundary(f)
      enddo
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!  Some precalculated pencils of data are passed in for efficiency,
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, contiguous, dimension(:,:,:,:) :: f
      type(boundary_condition) :: bc
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_boundconds(f,bc)
      enddo
!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_,llast)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dt_
      logical, intent(in) :: llast
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_after_timestep(f,df,dt_,llast)
      enddo
!
    endsubroutine special_after_timestep
!***********************************************************************
    subroutine set_init_parameters(Ntot,dsize,init_distr,init_distr2)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      real :: Ntot
      real, dimension(ndustspec) :: dsize, init_distr, init_distr2
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%set_init_parameters(Ntot,dsize,init_distr,init_distr2)
      enddo
!
    endsubroutine set_init_parameters
!***********************************************************************
    subroutine special_calc_spectra(f,spec,spec_hel,spec_2d,spec_2d_hel,lfirstcall,kind)
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(:) :: spec, spec_hel
      real, dimension(:,:) :: spec_2d, spec_2d_hel
      logical :: lfirstcall
      character(len=3) :: kind
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_calc_spectra(f,spec,spec_hel,spec_2d,spec_2d_hel,lfirstcall,kind)
      enddo
!
    endsubroutine special_calc_spectra
!***********************************************************************
    subroutine special_calc_spectra_byte(f,spec,spec_hel,lfirstcall,kind,len)
!
!  Not dispatched: the special modules are called through special_calc_spectra
!  (the byte variant was needed only by the former C trampolines).
!
      use Quiet
!
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(:) :: spec, spec_hel
      logical :: lfirstcall
      character, dimension(3) :: kind
      integer :: len
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(spec)
      call keep_compiler_quiet(spec_hel)
      call keep_compiler_quiet(lfirstcall)
      call keep_compiler_quiet(len)
!
    endsubroutine special_calc_spectra_byte
!***********************************************************************
    subroutine input_persist_special_id(id,done)
!
      integer :: id
      logical :: done
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%input_persist_special_id(id,done)
      enddo
!
    endsubroutine input_persist_special_id
!***********************************************************************
    subroutine input_persist_special
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%input_persist_special()
      enddo
!
    endsubroutine input_persist_special
!***********************************************************************
    logical function output_persistent_special()
!
      integer :: i
!
      output_persistent_special=.false.
      do i=1,n_special_modules
        if (hooks(i)%output_persistent_special()) then
          output_persistent_special=.true.
          return
        endif
      enddo
!
    endfunction output_persistent_special
!***********************************************************************
    subroutine special_particles_after_dtsub(f,dtsub,fp,dfp,ineargrid)
!
!  Possibility to modify fp in the end of a sub-time-step.
!
!  28-aug-18/ccyang: coded
!
      real, contiguous, dimension(:,:,:,:) :: f
      real :: dtsub
      real, dimension(:,:) :: fp, dfp
      integer, dimension(:,:) :: ineargrid
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_particles_after_dtsub(f,dtsub,fp,dfp,ineargrid)
      enddo
!
    endsubroutine special_particles_after_dtsub
!***********************************************************************
    subroutine calc_diagnostics_special(f,p)
!
      real, contiguous, dimension(:,:,:,:) :: f
      type(pencil_case) :: p
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%calc_diagnostics_special(f,p)
      enddo
!
    endsubroutine calc_diagnostics_special
!***********************************************************************
    subroutine load_variables_to_gpu_special
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%load_variables_to_gpu_special()
      enddo
!
    endsubroutine load_variables_to_gpu_special
!***********************************************************************
    subroutine calc_ode_diagnostics_special(f_ode)
!
      real, dimension(n_odevars) :: f_ode
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%calc_ode_diagnostics_special(f_ode)
      enddo
!
    endsubroutine calc_ode_diagnostics_special
!***********************************************************************
    subroutine pushpars2c(p_par)
!
      use Messages, only: fatal_error
      use Quiet
!
      integer, parameter :: n_pars=0
      integer(KIND=ikind8), dimension(n_pars) :: p_par
!
      call fatal_error('pushpars2c_special','This function should not be called!')
      call keep_compiler_quiet(p_par)
!
    endsubroutine pushpars2c
!***********************************************************************
    subroutine prep_rhs_special(f_ode)
!
      real, dimension(n_odevars) :: f_ode
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%prep_rhs_special(f_ode)
      enddo
!
    endsubroutine prep_rhs_special
!***********************************************************************
    subroutine special_before_boundary_diagnostics(f)
!
      real, contiguous, dimension(:,:,:,:) :: f
!
      integer :: i
!
      do i=1,n_special_modules
        call hooks(i)%special_before_boundary_diagnostics(f)
      enddo
!
    endsubroutine special_before_boundary_diagnostics
!***********************************************************************
  endmodule Special
