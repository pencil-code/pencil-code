!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!  special_interfaces.inc
!
!  Abstract interfaces of the hooks every Special module provides.
!
!  Included by the special module dispatcher (special.f90), which calls the
!  hooks of multiple special modules through procedure pointers of these
!  interfaces, and by check_special_hook_interfaces (nospecial.f90, copied into
!  every special module via its *_dummies.inc), which lets the compiler check
!  that the hook implementations have exactly these interfaces.
!
!  The f-array is passed as assumed-shape, contiguous array. Dummies without
!  INTENT here must be declared without INTENT in the implementations as well,
!  since the modules use them differently.
!
!  27-sep-26/TP: coded
!
  abstract interface
!
    subroutine iface_register_special()
    endsubroutine iface_register_special
!
    subroutine iface_register_particles_special(npvar)
      integer :: npvar
    endsubroutine iface_register_particles_special
!
    subroutine iface_initialize_special(f)
      real, contiguous, dimension(:,:,:,:) :: f
    endsubroutine iface_initialize_special
!
    subroutine iface_finalize_special(f)
      real, contiguous, dimension(:,:,:,:) :: f
    endsubroutine iface_finalize_special
!
    subroutine iface_read_special_pars(iomsg)
      character(len=*), intent(out) :: iomsg
    endsubroutine iface_read_special_pars
!
    subroutine iface_write_special_pars(unit)
      integer, intent(in) :: unit
    endsubroutine iface_write_special_pars
!
    subroutine iface_rprint_special(lreset,lwrite)
      logical :: lreset
      logical, optional :: lwrite
    endsubroutine iface_rprint_special
!
    subroutine iface_get_slices_special(f,slices)
      import :: slice_data
      real, contiguous, dimension(:,:,:,:) :: f
      type(slice_data) :: slices
    endsubroutine iface_get_slices_special
!
    subroutine iface_init_special(f)
      real, contiguous, dimension(:,:,:,:) :: f
    endsubroutine iface_init_special
!
    subroutine iface_dspecial_dt(f,df,p)
      import :: mx, my, mz, mvar, pencil_case
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(mx,my,mz,mvar) :: df
      type(pencil_case) :: p
    endsubroutine iface_dspecial_dt
!
    subroutine iface_calc_pencils_special(f,p)
      import :: pencil_case
      real, contiguous, dimension(:,:,:,:) :: f
      type(pencil_case) :: p
    endsubroutine iface_calc_pencils_special
!
    subroutine iface_pencil_interdep_special(lpencil_in)
      import :: npencils
      logical, dimension(npencils) :: lpencil_in
    endsubroutine iface_pencil_interdep_special
!
!  special_calc_hydro, _density, _dustdensity, _energy, _magnetic, _pscalar,
!  _chemistry (the pencil case is only read), and calc_diagnostics_special
!  (without df):
!
    subroutine iface_special_calc_rhs(f,df,p)
      import :: pencil_case
      real, contiguous, dimension(:,:,:,:) :: f
      real, contiguous, dimension(:,:,:,:) :: df
      type(pencil_case), intent(in) :: p
    endsubroutine iface_special_calc_rhs
!
    subroutine iface_calc_diagnostics_special(f,p)
      import :: pencil_case
      real, contiguous, dimension(:,:,:,:) :: f
      type(pencil_case) :: p
    endsubroutine iface_calc_diagnostics_special
!
    subroutine iface_special_calc_particles(f,df,fp,dfp,ineargrid)
      import :: mx, my, mz, mvar
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(:,:) :: fp, dfp
      integer, dimension(:,:) :: ineargrid
    endsubroutine iface_special_calc_particles
!
    subroutine iface_special_particles_bfre_bdary(f,fp,ineargrid)
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(:,:) :: fp
      integer, dimension(:,:) :: ineargrid
    endsubroutine iface_special_particles_bfre_bdary
!
    subroutine iface_special_particles_after_dtsub(f,dtsub,fp,dfp,ineargrid)
      real, contiguous, dimension(:,:,:,:) :: f
      real :: dtsub
      real, dimension(:,:) :: fp, dfp
      integer, dimension(:,:) :: ineargrid
    endsubroutine iface_special_particles_after_dtsub
!
!  special_before_boundary, special_before_boundary_diagnostics,
!  special_after_boundary:
!
    subroutine iface_special_boundary(f)
      real, contiguous, dimension(:,:,:,:) :: f
    endsubroutine iface_special_boundary
!
    subroutine iface_special_boundconds(f,bc)
      import :: boundary_condition
      real, contiguous, dimension(:,:,:,:) :: f
      type(boundary_condition) :: bc
    endsubroutine iface_special_boundconds
!
    subroutine iface_special_after_timestep(f,df,dt_,llast)
      import :: mx, my, mz, mvar
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dt_
      logical, intent(in) :: llast
    endsubroutine iface_special_after_timestep
!
    subroutine iface_set_init_parameters(Ntot,dsize,init_distr,init_distr2)
      import :: ndustspec
      real :: Ntot
      real, dimension(ndustspec) :: dsize, init_distr, init_distr2
    endsubroutine iface_set_init_parameters
!
    subroutine iface_special_calc_spectra(f,spectrum,spectrumhel,spectrum_2d,spectrum_2d_hel,lfirstcall,kind)
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(:) :: spectrum, spectrumhel
      real, dimension(:,:) :: spectrum_2d, spectrum_2d_hel
      logical :: lfirstcall
      character(len=3) :: kind
    endsubroutine iface_special_calc_spectra
!
    subroutine iface_special_calc_spectra_byte(f,spectrum,spectrumhel,lfirstcall,kind,len)
      real, contiguous, dimension(:,:,:,:) :: f
      real, dimension(:) :: spectrum, spectrumhel
      logical :: lfirstcall
      character, dimension(3) :: kind
      integer :: len
    endsubroutine iface_special_calc_spectra_byte
!
!  dspecial_dt_ode, pencil_criteria_special, input_persist_special,
!  load_variables_to_gpu_special:
!
    subroutine iface_special_noargs()
    endsubroutine iface_special_noargs
!
    subroutine iface_input_persist_special_id(id,done)
      integer :: id
      logical :: done
    endsubroutine iface_input_persist_special_id
!
    logical function iface_output_persistent_special()
    endfunction iface_output_persistent_special
!
!  calc_ode_diagnostics_special, prep_rhs_special:
!
    subroutine iface_special_ode(f_ode)
      import :: n_odevars
      real, dimension(n_odevars) :: f_ode
    endsubroutine iface_special_ode
!
  endinterface
