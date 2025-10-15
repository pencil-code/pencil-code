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

    implicit none

    include 'special.h'

    external caller, caller0, caller1, caller2, caller3, caller4, caller5, caller5_str5, caller6, &
             func_int_caller0
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
                          I_LOAD_VARIABLES_TO_GPU_SPECIAL=40
    
    integer, parameter :: n_subroutines=40
    integer, parameter :: n_special_modules_max=2
!
    integer :: n_special_modules
    character(LEN=256) :: special_modules_list = ''
    character(LEN=30), dimension(n_subroutines) :: special_subroutines=(/ &
                           'register_special              ', &
                           'register_particles_special    ', &
                           'initialize_special            ', &
                           'finalize_special              ', &
                           'read_special_init_pars        ', &
                           'write_special_init_pars       ', &
                           'read_special_run_pars         ', &
                           'write_special_run_pars        ', &
                           'rprint_special                ', &
                           'get_slices_special            ', &
                           'init_special                  ', &
                           'dspecial_dt                   ', &
                           'calc_pencils_special          ', &
                           'pencil_criteria_special       ', &
                           'pencil_interdep_special       ', &
                           'special_calc_hydro            ', &
                           'special_calc_density          ', &
                           'special_calc_dustdensity      ', &
                           'special_calc_energy           ', &
                           'special_calc_magnetic         ', &
                           'special_calc_pscalar          ', &
                           'special_calc_particles        ', &
                           'special_calc_chemistry        ', &
                           'special_boundconds            ', &
                           'special_before_boundary       ', &
                           'special_particles_bfre_bdary  ', &
                           'special_after_boundary        ', &
                           'special_after_timestep        ', &
                           'set_init_parameters           ', &
                           'special_calc_spectra          ', &
                           'special_calc_spectra_byte     ', &
                           'dspecial_dt_ode               ', &
                           'input_persist_special         ', &
                           'input_persist_special_id      ', &
                           'output_persistent_special     ', &
                           'special_particles_after_dtsub ', &
                           'calc_diagnostics_special      ',  &
                           'calc_ode_diagnostics_special  ',  &
                           'prep_rhs_special              ',  &
                           'load_variables_to_gpu_special '  &
                   /)

    integer(KIND=ikind8) :: libhandle
    integer(KIND=ikind8), dimension(n_special_modules_max,n_subroutines) :: special_sub_handles
    character(LEN=80) :: specific_subroutine

    contains
!****************************************************************************
  subroutine initialize_mult_special

    use Cdata, only: lroot, lreloading, iproc
    use General, only: parser, safe_string_replace
    use Messages, only: fatal_error
    use Mpicomm, only: mpibcast
    use Syscalls, only: extract_str, get_env_var

    integer, parameter :: RTLD_LAZY=0, RTLD_NOW=1

    character(LEN=128) :: line,parstr
    integer :: i,j,ipos
    character(LEN=40), dimension(n_special_modules_max) :: special_modules
    character(LEN=8) :: mod_prefix, mod_infix, mod_suffix
    integer(KIND=ikind8) :: sub_handle

    if (lreloading) return

    call get_env_var("PC_MODULES_LIST", special_modules_list)
    n_special_modules=parser(trim(special_modules_list),special_modules,' ')
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
!if (lroot) print*, 'specific_subroutine=',specific_subroutine

        sub_handle=dlsym_c(libhandle,trim(specific_subroutine)//char(0))
        if (sub_handle==0) &
          call fatal_error('initialize_mult_special','Error for symbol '// &
          trim(specific_subroutine)//' in module '//trim(special_modules(i))) 
        special_sub_handles(i,j) = sub_handle
      enddo
    enddo

  endsubroutine initialize_mult_special
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
        call caller0(special_sub_handles(i,I_REGISTER_SPECIAL))
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
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: i
!
      do i=1,n_special_modules
        call caller1(special_sub_handles(i,I_INITIALIZE_SPECIAL),f)
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
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: i
!
      do i=1,n_special_modules
        call caller1(special_sub_handles(i,I_FINALIZE_SPECIAL),f)
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
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
      integer :: i
!
      do i=1,n_special_modules
        call caller1(special_sub_handles(i,I_INIT_SPECIAL),f)
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
        call caller0(special_sub_handles(i,I_PENCIL_CRITERIA_SPECIAL))
      enddo

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
      integer :: i
!
      do i=1,n_special_modules
        call caller1(special_sub_handles(i,I_PENCIL_INTERDEP_SPECIAL),lpencil_in)
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
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      integer :: i
!
      do i=1,n_special_modules
        !call caller(special_sub_handles(i,I_CALC_PENCILS_SPECIAL),2,f,p)
        call caller2(special_sub_handles(i,I_CALC_PENCILS_SPECIAL),f,p)
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
        call caller0(special_sub_handles(i,I_DSPECIAL_DT_ODE))
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      integer :: i
!
      call special_calc_3par(f,df,p,I_DSPECIAL_DT)

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
        call caller1(special_sub_handles(i,I_REGISTER_PARTICLES_SPECIAL),npvar)
      enddo
!
    endsubroutine register_particles_special
!*********************************************************************** 
    subroutine read_special_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      integer :: i
!
      do i=1,n_special_modules
        call caller1(special_sub_handles(i,I_READ_SPECIAL_INIT_PARS),iostat)
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
        call caller1(special_sub_handles(i,I_WRITE_SPECIAL_INIT_PARS),unit)
      enddo
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      integer :: i
!
      do i=1,n_special_modules
        call caller1(special_sub_handles(i,I_READ_SPECIAL_RUN_PARS),iostat)
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
        call caller1(special_sub_handles(i,I_WRITE_SPECIAL_RUN_PARS),unit)
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
      logical :: lreset,lwrite
      integer :: i
!
      do i=1,n_special_modules
        !call caller(special_sub_handles(i,I_RPRINT_SPECIAL),2,lreset,lwrite)
        call caller2(special_sub_handles(i,I_RPRINT_SPECIAL),lreset,lwrite)
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
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      integer :: i
!
      do i=1,n_special_modules
        !call caller(special_sub_handles(i,I_GET_SLICES_SPECIAL),2,f,slices)
        call caller2(special_sub_handles(i,I_GET_SLICES_SPECIAL),f,slices)
      enddo
!
    endsubroutine get_slices_special
!*********************************************************************** 
    subroutine special_calc_3par(f,df,p,modind)
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
      integer, intent(in) :: modind
!
      integer :: i
!
      do i=1,n_special_modules
        !call caller(special_sub_handles(i,modind),3,f,df,p)
        call caller3(special_sub_handles(i,modind),f,df,p)
      enddo
!
    endsubroutine special_calc_3par
!*********************************************************************** 
    subroutine special_calc_hydro(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  momentum equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

      call special_calc_3par(f,df,p,I_SPECIAL_CALC_HYDRO)

    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  continuity equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      call special_calc_3par(f,df,p,I_SPECIAL_CALC_DENSITY)
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
!
      call special_calc_3par(f,df,p,I_SPECIAL_CALC_DUSTDENSITY)
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
!
      call special_calc_3par(f,df,p,I_SPECIAL_CALC_ENERGY)
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
!
      call special_calc_3par(f,df,p,I_SPECIAL_CALC_MAGNETIC)
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
!
      call special_calc_3par(f,df,p,I_SPECIAL_CALC_PSCALAR)
!
    endsubroutine special_calc_pscalar
!*********************************************************************** 
    subroutine special_calc_particles(f,df,fp,dfp,ineargrid)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (:,:) :: fp
      real, dimension (:,:) :: dfp
      integer, dimension(:,:) :: ineargrid
!
      integer :: i
!
      do i=1,n_special_modules
        !call caller(special_sub_handles(i,I_SPECIAL_CALC_PARTICLES),5,f,df,fp,dfp,ineargrid)
        call caller5(special_sub_handles(i,I_SPECIAL_CALC_PARTICLES),f,df,fp,dfp,ineargrid)
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (:,:), intent(in) :: fp
      integer, dimension(:,:) :: ineargrid
!
      integer :: i
!
      do i=1,n_special_modules
        !call caller(special_sub_handles(i,I_SPECIAL_PARTICLES_BFRE_BDARY),3,f,fp,ineargrid)
        call caller3(special_sub_handles(i,I_SPECIAL_PARTICLES_BFRE_BDARY),f,fp,ineargrid)
      enddo
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
!
      call special_calc_3par(f,df,p,I_SPECIAL_CALC_CHEMISTRY)
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
      integer :: i
!
      do i=1,n_special_modules
        call caller1(special_sub_handles(i,I_SPECIAL_BEFORE_BOUNDARY),f)
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      integer :: i
!
      do i=1,n_special_modules
        call caller1(special_sub_handles(i,I_SPECIAL_AFTER_BOUNDARY),f)
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition), intent(in) :: bc
!
      integer :: i
!
      do i=1,n_special_modules
        !call caller(special_sub_handles(i,I_SPECIAL_BOUNDCONDS),2,f,bc)
        call caller2(special_sub_handles(i,I_SPECIAL_BOUNDCONDS),f,bc)
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
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_
      logical, intent(in) :: llast
!
      integer :: i
!
      do i=1,n_special_modules
        !call caller(special_sub_handles(i,I_SPECIAL_AFTER_TIMESTEP), &
        !            4,f,df,dt_,llast)
        call caller4(special_sub_handles(i,I_SPECIAL_AFTER_TIMESTEP), &
                    f,df,dt_,llast)
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
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(ndustspec) :: dsize,init_distr,init_distr2
      real :: Ntot
!
      integer :: i
!
      do i=1,n_special_modules
        !call caller(special_sub_handles(i,I_SET_INIT_PARAMETERS),4,Ntot,dsize,init_distr,init_distr2)
        call caller4(special_sub_handles(i,I_SET_INIT_PARAMETERS),Ntot,dsize,init_distr,init_distr2)
      enddo
!
    endsubroutine set_init_parameters
!***********************************************************************
    subroutine special_calc_spectra(f,spec,spec_hel,spec_2d,spec_2d_hel,lfirstcall,kind)

      use General, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:) :: spec,spec_hel
      real, dimension (:,:) :: spec_2d,spec_2d_hel
      logical :: lfirstcall
      character(LEN=3) :: kind
!
      integer :: i
      do i=1,n_special_modules
        call caller6(special_sub_handles(i,I_SPECIAL_CALC_SPECTRA_BYTE),f, &
                     spec,spec_hel,lfirstcall,kind,3)
!                    spec,spec_hel,spec_2d,spec_2d_hel,lfirstcall,kind,3)
!        call caller5_str5(special_sub_handles(i,I_SPECIAL_CALC_SPECTRA),f, &
!                         spec, spec_hel, lfirstcall, kind)
      enddo

      call keep_compiler_quiet(spec_2d,spec_2d_hel)

    endsubroutine special_calc_spectra
!***********************************************************************
    subroutine special_calc_spectra_byte(f,spec,spec_hel,lfirstcall,kind,len)

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:) :: spec,spec_hel
      logical :: lfirstcall
      integer :: len
      character, dimension(len) :: kind

    endsubroutine special_calc_spectra_byte
!*********************************************************************** 
    subroutine input_persist_special_id(id,done)
!
      integer :: id
      logical :: done

      integer :: i
      do i=1,n_special_modules
        call caller2(special_sub_handles(i,I_INPUT_PERSIST_SPECIAL_ID),id,done)
      enddo

    endsubroutine input_persist_special_id
!*********************************************************************** 
    subroutine input_persist_special
!
      integer :: i
      do i=1,n_special_modules
        call caller0(special_sub_handles(i,I_INPUT_PERSIST_SPECIAL))
      enddo

    endsubroutine input_persist_special
!*********************************************************************** 
    function output_persistent_special() result(ret)

      logical :: ret
      integer :: func_int_caller0
      integer :: i

      ret=.true.
      do i=1,n_special_modules
        ret=ret.and.(func_int_caller0(special_sub_handles(i,I_OUTPUT_PERSISTENT_SPECIAL))==1)
      enddo

    endfunction output_persistent_special        
!***********************************************************************
    subroutine special_particles_after_dtsub(f, dtsub, fp, dfp, ineargrid)
!
!  Possibility to modify fp in the end of a sub-time-step.
!
!  28-aug-18/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, intent(in) :: dtsub
      real, dimension(:,:), intent(in) :: fp, dfp
      integer, dimension(:,:), intent(in) :: ineargrid
!
      integer :: i
      do i=1,n_special_modules
        call caller5(special_sub_handles(i,I_SPECIAL_PARTICLES_AFTER_DTSUB),f,dtsub,fp,dfp,ineargrid)
      enddo
!
    endsubroutine special_particles_after_dtsub
!***********************************************************************
    subroutine calc_diagnostics_special(f,p)
      real, dimension(mx,my,mz,mfarray) :: f
      type(pencil_case) :: p
      
      integer :: i
      do i=1,n_special_modules
        call caller2(special_sub_handles(i,I_CALC_DIAGNOSTICS_SPECIAL),f,p)
      enddo
    endsubroutine calc_diagnostics_special
!***********************************************************************
    subroutine load_variables_to_gpu_special
      integer :: i
      do i=1,n_special_modules
        call caller0(special_sub_handles(i,I_LOAD_VARIABLES_TO_GPU_SPECIAL))
      enddo
    endsubroutine load_variables_to_gpu_special
!***********************************************************************
    subroutine calc_ode_diagnostics_special(f_ode)
      real, dimension(max_n_odevars), intent(IN) :: f_ode
      integer :: i
      do i=1,n_special_modules
        call caller1(special_sub_handles(i,I_CALC_ODE_DIAGNOSTICS_SPECIAL),f_ode)
      enddo
    endsubroutine calc_ode_diagnostics_special
!***********************************************************************
    subroutine pushpars2c(p_par)

      use Messages, only: not_implemented

      integer, parameter :: n_pars=0
      integer(KIND=ikind8), dimension(n_pars) :: p_par

      call not_implemented('special_pushpars2c','')

    endsubroutine pushpars2c
!*********************************************************************** 
    subroutine prep_rhs_special
      integer :: i
      do i=1,n_special_modules
        call caller0(special_sub_handles(i,I_PREP_RHS_SPECIAL))
      enddo
    endsubroutine prep_rhs_special
!***********************************************************************
  endmodule Special
