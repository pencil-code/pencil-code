! $Id$
!
!  subroutines in the chosen set of physics modules.
!
module Equ
!
  use Cdata
  use Messages
  use Boundcond
  use Mpicomm
  use Grid, only: calc_pencils_grid, get_grid_mn
!
  implicit none
!
  public :: pde, debug_imn_arrays, initialize_pencils
  public :: impose_floors_ceilings, finalize_diagnostics
  public :: read_diagnostics_accumulators, write_diagnostics
!
  private
!
!
  contains
!***********************************************************************
    include 'pencil_init.inc' ! defines subroutine initialize_pencils()
!***********************************************************************
    subroutine pde(f,df,p)
!
!  Call the different evolution equations.
!
!  10-sep-01/axel: coded
!  12-may-12/MR: call of density_before_boundary added for boussinesq;
!                moved call of timing after call of anelastic_after_mn
!  26-aug-13/MR: added call of diagnostic for imaginary parts
!   9-jun-15/MR: call of gravity_after_boundary added
!  24-sep-16/MR: added offset manipulation for second derivatives in complete one-sided fornulation.
!   5-jan-17/MR: removed mn-offset manipulation
!  14-feb-17/MR: adaptations for use of GPU kernels in calculating the rhss of the pde
!
      use Chiral
      use Chemistry
      use Density
      use Detonate, only: detonate_before_boundary
      use Diagnostics
      use Dustdensity
      use Energy
      use EquationOfState
      use Forcing, only: forcing_after_boundary
!                         
! To check ghost cell consistency, please uncomment the following line:
!     use Ghost_check, only: check_ghosts_consistency
      use GhostFold, only: fold_df, fold_df_3points
      use Gpu
      use Gravity
      use Hydro
      use Interstellar, only: interstellar_before_boundary
      use Magnetic
      use Magnetic_meanfield, only: meanfield_after_boundary
      use Hypervisc_strict, only: hyperviscosity_strict
      use Hyperresi_strict, only: hyperresistivity_strict
      use NeutralDensity, only: neutraldensity_after_boundary
      use NSCBC
      use Particles_main
      use Poisson
      use Pscalar
      use PointMasses
      use Polymer
      use Radiation
      use Selfgravity
      use Shear
      use Shock, only: shock_before_boundary, calc_shock_profile_simple
      use Solid_Cells, only: update_solid_cells, dsolid_dt_integrate
      use Special, only: special_before_boundary,special_after_boundary
      use Sub
      use Testfield
      use Testflow
      use Testscalar
      use Viscosity, only: viscosity_after_boundary
      use Grid, only: coarsegrid_interp
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(inout):: f       ! inout due to lshift_datacube_x,
                              ! density floor, or velocity ceiling
      intent(out)  :: df,p
!
      logical :: early_finalize
      real, dimension(1)  :: mass_per_proc
!
!  Print statements when they are first executed.
!
      headtt = headt .and. lfirst .and. lroot
!
      if (headtt.or.ldebug) print*,'pde: ENTER'
      if (headtt) call svn_id( &
           "$Id$")
!
!  Initialize counter for calculating and communicating print results.
!  Do diagnostics only in the first of the 3 (=itorder) substeps.
!
      ldiagnos   =lfirst.and.lout
      l1davgfirst=lfirst.and.l1davg
      l2davgfirst=lfirst.and.l2davg
!
!  Derived diagnostics switches.
!
      l1dphiavg=lcylinder_in_a_box.and.l1davgfirst
!
!  For chemistry with LSODE
!
      lchemonly=.false.
!
!  Record times for diagnostic and 2d average output.
!
      if (ldiagnos   ) tdiagnos  =t ! (diagnostics are for THIS time)
      if (l1davgfirst) t1ddiagnos=t ! (1-D averages are for THIS time)
      if (l2davgfirst)  then
        t2davgfirst=t ! (2-D averages are for THIS time)
!
!  [AB: Isn't it true that not all 2-D averages use rcyl_mn?
!  lwrite_phiaverages=T is required, and perhaps only that.]
!  [BD: add also the z_mn dependency]
!
        lpencil(i_rcyl_mn)=.true.
        lpencil(i_z_mn)=.true.
      endif
!
!  Shift entire data cube by one grid point at the beginning of each
!  time-step. Useful for smearing out possible x-dependent numerical
!  diffusion, e.g. in a linear shear flow.
!
      if (lfirst .and. lshift_datacube_x) then
        call boundconds_x(f)
        do  n=n1,n2; do m=m1,m2
          f(:,m,n,:)=cshift(f(:,m,n,:),1,1)
        enddo; enddo
      endif
!
!  Need to finalize communication early either for test purposes, or
!  when radiation transfer of global ionization is calculated.
!  This could in principle be avoided (but it not worth it now)
!
      early_finalize=test_nonblocking.or. &
                     leos_ionization.or.lradiation_ray.or. &
                     lhyperviscosity_strict.or.lhyperresistivity_strict.or. &
                     ltestscalar.or.ltestfield.or.ltestflow.or. &
                     lparticles_spin.or.lsolid_cells.or. &
                     lchemistry.or.lweno_transport .or. lbfield .or. & 
!                     lslope_limit_diff .or. lvisc_smag .or. &
                     lvisc_smag .or. &
                     lyinyang .or. lgpu .or. &   !!!
                     ncoarse>1 
!
!  Write crash snapshots to the hard disc if the time-step is very low.
!  The user must have set crash_file_dtmin_factor>0.0 in &run_pars for
!  this to be done.
!
      if (crash_file_dtmin_factor > 0.0) call output_crash_files(f)
!
!  For debugging purposes impose minimum or maximum value on certain variables.
!
      call impose_floors_ceilings(f)   !MR: too early, f modifications come below
!
!  Apply global boundary conditions to particle positions and communicate
!  migrating particles between the processors.
!
      if (lparticles) call particles_boundconds(f)
      if (lpointmasses) call boundconds_pointmasses
!
!  Calculate the potential of the self gravity. Must be done before
!  communication in order to be able to take the gradient of the potential
!  later.
!
      call calc_selfpotential(f)
!
!  Check for dust grain mass interval overflows
!  (should consider having possibility for all modules to fiddle with the
!   f array before boundary conditions are sent)
!
      if (.not. lchemistry) then
        if (ldustdensity) call null_dust_vars(f)
        if (ldustdensity .and. lmdvar .and. lfirst) call redist_mdbins(f)
      endif
!
!  Call "before_boundary" hooks (for f array precalculation)
!
      if (linterstellar) call interstellar_before_boundary(f)
      if (ldensity.or.lboussinesq) call density_before_boundary(f)
      if (lhydro.or.lhydro_kinematic) call hydro_before_boundary(f)
      if (lmagnetic)     call magnetic_before_boundary(f)
                         call energy_before_boundary(f)
      if (lshear)        call shear_before_boundary(f)
      if (lchiral)       call chiral_before_boundary(f)
      if (lspecial)      call special_before_boundary(f)
      if (ltestflow)     call testflow_before_boundary(f)
      if (ltestfield)    call testfield_before_boundary(f)
      if (lparticles)    call particles_before_boundary(f)
      if (lpscalar)      call pscalar_before_boundary(f)
      if (ldetonate)     call detonate_before_boundary(f)
      if (lchemistry)    call chemistry_before_boundary(f)
      if (lparticles.and.lspecial) call particles_special_bfre_bdary(f)
      if (lshock)        call shock_before_boundary(f)
!
!  Prepare x-ghost zones; required before f-array communication
!  AND shock calculation
!
      call boundconds_x(f)
!
!  Initiate (non-blocking) communication and do boundary conditions.
!  Required order:
!  1. x-boundaries (x-ghost zones will be communicated) - done above
!  2. communication
!  3. y- and z-boundaries
!
      if (nghost>0) then
        if (ldebug) print*,'pde: before initiate_isendrcv_bdry'
        call initiate_isendrcv_bdry(f)
        if (early_finalize) then
          call finalize_isendrcv_bdry(f)
          if (lcoarse) call coarsegrid_interp(f)   ! after boundconds_x???
          call boundconds_y(f)
          call boundconds_z(f)
        endif
      endif
!
! update solid cell "ghost points". This must be done in order to get the
! correct boundary layer close to the solid geometry, i.e. no-slip conditions.
!
      call update_solid_cells(f)
!
!  For sixth order momentum-conserving, symmetric hyperviscosity with positive
!  definite heating rate we need to precalculate the viscosity term. The
!  restivitity term for sixth order hyperresistivity with positive definite
!  heating rate must also be precalculated.
!
      if (lhyperviscosity_strict)   call hyperviscosity_strict(f)
      if (lhyperresistivity_strict) call hyperresistivity_strict(f)
!
!  Dynamically set the (hyper-)diffusion coefficients
!
      if (ldynamical_diffusion) call set_dyndiff_coeff(f)
!
!  Calculate the characteristic velocity
!  for slope limited diffusion
!
!      if (lslope_limit_diff.and.llast) then
!        f(2:mx-2,2:my-2,2:mz-2,iFF_char_c)=0.
!print*,'vor magnetic:', maxval(f(2:mx-2,2:my-2,2:mz-2,iFF_char_c))
!        call update_char_vel_energy(f)
!        call update_char_vel_magnetic(f)
!        call update_char_vel_hydro(f)
        !call update_char_vel_density(f)
        !f(2:mx-2,2:my-2,2:mz-2,iFF_char_c)=sqrt(f(2:mx-2,2:my-2,2:mz-2,iFF_char_c))
!  JW: for hydro it is done without sqrt
        !if (ldiagnos) print*, 'max(char_c)=', maxval(f(2:mx-2,2:my-2,2:mz-2,iFF_char_c))
!      endif
!
!  For calculating the pressure gradient directly from the pressure (which is
!  derived from the basic thermodynamical variables), we need to fill in the
!  pressure in the f array.
!
      call fill_farray_pressure(f)
!
!  Set inverse timestep to zero before entering loop over m and n.
!  If we want to have a logarithmic time advance, we want set this here
!  as the maximum. All other routines can then still make it shorter.
!
      if (lfirst.and.ldt) then
        if (dtmax/=0.0) then
          if (lfractional_tstep_advance) then
            dt1_max=1./(dt_incr*t)
          else
            dt1_max=1./dtmax
          endif
        else
          dt1_max=0.0
        endif
      endif
!
!  Calculate ionization degree (needed for thermodynamics)
!  Radiation transport along rays. If lsingle_ray, then this
!  is only used for visualization and only needed when lvideo
!  (but this is decided in radtransfer itself)
!
      if (leos_ionization.or.leos_temperature_ionization) call ioncalc(f)
      if (lradiation_ray) call radtransfer(f)
!
!  Calculate shock profile (simple).
!
      if (lshock) call calc_shock_profile_simple(f)
!
!  Call "after" hooks (for f array precalculation). This may imply
!  calculating averages (some of which may only be required for certain
!  settings in hydro of the testfield procedure (only when lsoca=.false.),
!  for example. They used to be or are still called hydro_after_boundary etc,
!  and will soon be renamed to hydro_after_boundary.
!
!  Important to note that the processor boundaries are not full updated 
!  at this point, even if the name 'after_boundary' suggesting this.
!  Use early_finalize in this case.
!  MR+joern+axel, 8.10.2015
!
      call timing('pde','before "after_boundary" calls')
!
      if (lhydro)                 call hydro_after_boundary(f)
      if (lviscosity)             call viscosity_after_boundary(f)
      if (lmagnetic)              call magnetic_after_boundary(f)
      if (ldustdensity)           call dustdensity_after_boundary(f)
      if (lenergy)                call energy_after_boundary(f)
      if (lgrav)                  call gravity_after_boundary(f)
      if (lforcing)               call forcing_after_boundary(f)
      if (lpolymer)               call calc_polymer_after_boundary(f)
      if (ltestscalar)            call testscalar_after_boundary(f)
      if (ltestfield)             call testfield_after_boundary(f)
!AB: quick fix
      !if (ltestfield)             call testfield_after_boundary(f,p)
      if (ldensity)               call density_after_boundary(f)
      if (lneutraldensity)        call neutraldensity_after_boundary(f)
      if (ltestflow)              call calc_ltestflow_nonlin_terms(f,df)  ! should not use df!
      if (lmagn_mf)               call meanfield_after_boundary(f)
      if (lspecial)               call special_after_boundary(f)
!
!  Calculate quantities for a chemical mixture. This is done after
!  communication has finalized since many of the arrays set up here
!  are not communicated, and in this subroutine also ghost zones are calculated.
!
      if (lchemistry .and. ldensity) call calc_for_chem_mixture(f)
!
      call timing('pde','after "after_boundary" calls')
!
      if (lgpu) then
        call rhs_gpu(f,itsub,early_finalize)
        if (ldiagnos.or.l1davgfirst.or.l1dphiavg.or.l2davgfirst) then
          call do_rest_diagnostics_tasks(f)
          call copy_farray_from_GPU(f)
          call calc_all_module_diagnostics(f,p)
        endif
      else
        call test_rhs(f,df,p,mass_per_proc,early_finalize,rhs_cpu,rhs_cpu)
        call rhs_cpu(f,df,p,mass_per_proc,early_finalize)
      endif
!
!  Doing df-related work which cannot be finished inside the main mn-loop.
!  (At the moment relevant for anelastic and Schur flows.)
!
      call density_after_mn(f, df, mass_per_proc)
!
      call timing('pde','after the end of the mn_loop')
!
!  Integrate diagnostics related to solid cells (e.g. drag and lift).
!
      if (lsolid_cells) call dsolid_dt_integrate
!
!  Calculate the gradient of the potential if there is room allocated in the
!  f-array.
!
      if (igpotselfx/=0) then
        call initiate_isendrcv_bdry(f,igpotselfx,igpotselfz)
        call finalize_isendrcv_bdry(f,igpotselfx,igpotselfz)
        call boundconds_x(f,igpotselfx,igpotselfz)
        call boundconds_y(f,igpotselfx,igpotselfz)
        call boundconds_z(f,igpotselfx,igpotselfz)
      endif
!
!  Change df and dfp according to the chosen particle modules.
!
      if (lparticles) then
        call particles_pde_blocks(f,df)
        call particles_pde(f,df)
      endif
!
      if (lpointmasses) call pointmasses_pde(f,df)
!
!  Electron inertia: our df(:,:,:,iax:iaz) so far is
!  (1 - l_e^2\Laplace) daa, thus to get the true daa, we need to invert
!  that operator.
!  [wd-aug-2007: This should be replaced by the more general stuff with the
!   Poisson solver (so l_e can be non-constant), so at some point, we can
!   remove/replace this]
!
!      if (lelectron_inertia .and. inertial_length/=0.) then
!        do iv = iax,iaz
!          call inverse_laplacian_semispectral(df(:,:,:,iv), H=linertial_2)
!        enddo
!        df(:,:,:,iax:iaz) = -df(:,:,:,iax:iaz) * linertial_2
!      endif
!
!  Take care of flux-limited diffusion
!  This is now commented out, because we always use radiation_ray instead.
!
!--   if (lradiation_fld) f(:,:,:,idd)=DFF_new
!
!  Fold df from first ghost zone into main df.
!
      if (lfold_df) then
        if (lhydro .and. (.not. lpscalar) .and. (.not. lchemistry)) then
          call fold_df(df,iux,iuz)
        else
          call fold_df(df,iux,mvar)
        endif
      endif
      if (lfold_df_3points) call fold_df_3points(df,iux,mvar)
!
!  -------------------------------------------------------------
!  NO CALLS MODIFYING DF BEYOND THIS POINT (APART FROM FREEZING)
!  -------------------------------------------------------------
!
!  Freezing must be done after the full (m,n) loop, as df may be modified
!  outside of the considered pencil.
!
      call freeze(df,p)

!  Boundary treatment of the df-array.
!
!  This is a way to impose (time-
!  dependent) boundary conditions by solving a so-called characteristic
!  form of the fluid equations on the boundaries, as opposed to setting
!  actual values of the variables in the f-array. The method is called
!  Navier-Stokes characteristic boundary conditions (NSCBC).
!
!  The treatment should be done after the y-z-loop, but before the Runge-
!  Kutta solver adds to the f-array.
!
      if (lnscbc) call nscbc_boundtreat(f,df)

      if (lfirst) then
!$      if(num_of_diag_iter_done==nyz .and. .not. lstarted_finalizing_diagnostics) then
          call finalize_diagnostics
!$      endif
      endif
!
!  Calculate rhoccm and cc2m (this requires that these are set in print.in).
!  Broadcast result to other processors. This is needed for calculating PDFs.
!
!      if (idiag_rhoccm/=0) then
!        if (iproc==0) rhoccm=fname(idiag_rhoccm)
!        call mpibcast_real(rhoccm)
!      endif
!
!      if (idiag_cc2m/=0) then
!        if (iproc==0) cc2m=fname(idiag_cc2m)
!        call mpibcast_real(cc2m)
!      endif
!
!      if (idiag_gcc2m/=0) then
!        if (iproc==0) gcc2m=fname(idiag_gcc2m)
!        call mpibcast_real(gcc2m)
!      endif
!
!  Reset lwrite_prof.
!
      lwrite_prof=.false.
!
    endsubroutine pde
!***********************************************************************
    subroutine pde_gpu(f,df,p)
!
!  Call the different evolution equations.
!
!  10-sep-01/axel: coded
!  12-may-12/MR: call of density_before_boundary added for boussinesq;
!                moved call of timing after call of anelastic_after_mn
!  26-aug-13/MR: added call of diagnostic for imaginary parts
!   9-jun-15/MR: call of gravity_after_boundary added
!  24-sep-16/MR: added offset manipulation for second derivatives in complete one-sided fornulation.
!   5-jan-17/MR: removed mn-offset manipulation
!  14-feb-17/MR: adaptations for use of GPU kernels in calculating the rhss of the pde
!
      use Chiral
      use Chemistry
      use Density
      use Detonate, only: detonate_before_boundary
      use Diagnostics
      use Dustdensity
      use Energy
      use EquationOfState
      use Forcing, only: forcing_after_boundary
!                         
! To check ghost cell consistency, please uncomment the following line:
!     use Ghost_check, only: check_ghosts_consistency
      use GhostFold, only: fold_df, fold_df_3points
      use Gpu
      use Gravity
      use Hydro
      use Interstellar, only: interstellar_before_boundary
      use Magnetic
      use Magnetic_meanfield, only: meanfield_after_boundary
      use Hypervisc_strict, only: hyperviscosity_strict
      use Hyperresi_strict, only: hyperresistivity_strict
      use NeutralDensity, only: neutraldensity_after_boundary
      use NSCBC
      use Particles_main
      use Poisson
      use Pscalar
      use PointMasses
      use Polymer
      use Radiation
      use Selfgravity
      use Shear
      use Shock, only: shock_before_boundary, calc_shock_profile_simple
      use Solid_Cells, only: update_solid_cells, dsolid_dt_integrate
      use Special, only: special_before_boundary,special_after_boundary
      use Sub
      use Testfield
      use Testflow
      use Testscalar
      use Viscosity, only: viscosity_after_boundary
      use Grid, only: coarsegrid_interp
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(inout):: f       ! inout due to lshift_datacube_x,
                              ! density floor, or velocity ceiling
      intent(out)  :: df,p
!
      logical :: early_finalize
      real, dimension(1)  :: mass_per_proc
!
!  Print statements when they are first executed.
!
      headtt = headt .and. lfirst .and. lroot
!
      if (headtt.or.ldebug) print*,'pde: ENTER'
      if (headtt) call svn_id( &
           "$Id$")
!
!  Initialize counter for calculating and communicating print results.
!  Do diagnostics only in the first of the 3 (=itorder) substeps.
!
      ldiagnos   =lfirst.and.lout
      l1davgfirst=lfirst.and.l1davg
      l2davgfirst=lfirst.and.l2davg
!
!  Derived diagnostics switches.
!
      l1dphiavg=lcylinder_in_a_box.and.l1davgfirst
!
!  For chemistry with LSODE
!
      lchemonly=.false.
!
!  Record times for diagnostic and 2d average output.
!
      if (ldiagnos   ) tdiagnos  =t ! (diagnostics are for THIS time)
      if (l1davgfirst) t1ddiagnos=t ! (1-D averages are for THIS time)
      if (l2davgfirst)  then
        t2davgfirst=t ! (2-D averages are for THIS time)
!
!  [AB: Isn't it true that not all 2-D averages use rcyl_mn?
!  lwrite_phiaverages=T is required, and perhaps only that.]
!  [BD: add also the z_mn dependency]
!
        lpencil(i_rcyl_mn)=.true.
        lpencil(i_z_mn)=.true.
      endif
!
!  Shift entire data cube by one grid point at the beginning of each
!  time-step. Useful for smearing out possible x-dependent numerical
!  diffusion, e.g. in a linear shear flow.
!
      if (lfirst .and. lshift_datacube_x) then
        call boundconds_x(f)
        do  n=n1,n2; do m=m1,m2
          f(:,m,n,:)=cshift(f(:,m,n,:),1,1)
        enddo; enddo
      endif
!
!  Need to finalize communication early either for test purposes, or
!  when radiation transfer of global ionization is calculated.
!  This could in principle be avoided (but it not worth it now)
!
      early_finalize=test_nonblocking.or. &
                     leos_ionization.or.lradiation_ray.or. &
                     lhyperviscosity_strict.or.lhyperresistivity_strict.or. &
                     ltestscalar.or.ltestfield.or.ltestflow.or. &
                     lparticles_spin.or.lsolid_cells.or. &
                     lchemistry.or.lweno_transport .or. lbfield .or. & 
!                     lslope_limit_diff .or. lvisc_smag .or. &
                     lvisc_smag .or. &
                     lyinyang .or. lgpu .or. &   !!!
                     ncoarse>1 
!
      call rhs_gpu(f,itsub,early_finalize)
      if (ldiagnos.or.l1davgfirst.or.l1dphiavg.or.l2davgfirst) then
        call do_rest_diagnostics_tasks(f)
        call copy_farray_from_GPU(f)
        call calc_all_module_diagnostics(f,p)
      endif

      if (lfirst) then
!$      if(num_of_diag_iter_done==nyz .and. .not. lstarted_finalizing_diagnostics) then
          call finalize_diagnostics
!$      endif
      endif
!
    endsubroutine pde_gpu
!***********************************************************************
    subroutine read_diagnostic_flags
    use Chemistry
    l1davgfirst = l1davgfirst_save 
    ldiagnos = ldiagnos_save 
    l1dphiavg = l1dphiavg_save 
    l2davgfirst = l2davgfirst_save 

    lout = lout_save 
    l1davg = l1davg_save
    l2davg = l2davg_save 
    lout_sound = lout_sound_save
    lvideo = lvideo_save
    lwrite_slices = lwrite_slices_save

    lchemistry_diag = lchemistry_diag_save
    it = it_save

    endsubroutine read_diagnostic_flags
!****************************************************************************
    subroutine read_diagnostics_accumulators

    use Chemistry
    use Diagnostics

    if (allocated(fname))      fname = p_fname
    if (allocated(fnamex))     fnamex = p_fnamex
    if (allocated(fnamey))     fnamey = p_fnamey 
    if (allocated(fnamez))     fnamez = p_fnamez
    if (allocated(fnamer))     fnamer = p_fnamer
    if (allocated(fnamexy))    fnamexy = p_fnamexy
    if (allocated(fnamexz))    fnamexz = p_fnamexz
    if (allocated(fnamerz))    fnamerz = p_fnamerz
    if (allocated(fname_keep)) fname_keep = p_fname_keep
    if (allocated(fname_sound))fname_sound = p_fname_sound
    if (allocated(ncountsz))   ncountsz = p_ncountsz
    call read_diagnostic_flags
    call diagnostics_read_diag_accum
    call chemistry_read_diag_accum


    endsubroutine read_diagnostics_accumulators
!***********************************************************************
   subroutine write_diagnostics(f)
    use Chemistry
    use Slices
    use Diagnostics
    real, dimension (mx,my,mz,mfarray) :: f

            call read_diagnostics_accumulators
            lstarted_writing_diagnostics = .true.
        !
        !  Print diagnostic averages to screen and file.
        !
        !task here
            if (lout) then
                call prints
                if (lchemistry_diag) call write_net_reaction
            endif
        !
            if (l1davg) call write_1daverages
            if (l2davg) call write_2daverages
        !
            if (lout_sound) then
                call write_sound(tsound)
                lout_sound = .false.
            endif
        !
        !  Write slices (for animation purposes).
        !
            if (lvideo .and. lwrite_slices) call wvid(f)
        !
            lwritten_diagnostics = .true.
    endsubroutine write_diagnostics
!***********************************************************************
  subroutine write_diagnostics_accumulators
!$  use OMP_lib
    integer :: imn
      if (ldiagnos .and. allocated(fname)) then
        p_fname = fname
      endif

      if (l1davgfirst) then
        if (allocated(fnamex)) p_fnamex =  fnamex
        if (allocated(fnamey)) p_fnamey =  fnamey
        if (allocated(fnamez)) p_fnamez =  fnamez
        if (allocated(fnamer)) p_fnamer =  fnamer
      endif

      if (l2davgfirst) then
        if (allocated(fnamexy)) p_fnamexy =  fnamexy
        if (allocated(fnamexz)) p_fnamexz =  fnamexz
        if (allocated(fnamerz)) p_fnamerz =  fnamerz
      endif

      if (allocated(fname_keep)) p_fname_keep =  fname_keep
      if (allocated(fname_sound)) p_fname_sound =  fname_sound
      if (allocated(ncountsz)) p_ncountsz =  ncountsz
    endsubroutine write_diagnostics_accumulators
!***********************************************************************
    subroutine do_rest_diagnostics_tasks(f)
    use Slices
    use Diagnostics
    use Chemistry
    real, dimension (mx,my,mz,mfarray) :: f

       !$omp taskwait
       !wait for diagnostics
        do while(num_of_diag_iter_done < nyz)
        enddo
        !if not started finalization do itself
        if(.not. lstarted_finalizing_diagnostics) call finalize_diagnostics
        !wait for finalization
        do while(.not. lfinalized_diagnostics)
        enddo
        !if not started writing diagnostics do itself
        if (.not. lstarted_writing_diagnostics) call write_diagnostics(f)
        !wait for writing
        do while(.not. lwritten_diagnostics)
        enddo
    endsubroutine do_rest_diagnostics_tasks
!***********************************************************************
    subroutine init_reduc_pointers
!
!  Initializes pointers used in diagnostics_reductions
!
!  30-mar-23/TP: Coded
!  
      use Diagnostics

      if (allocated(fname))      p_fname => fname
      if (allocated(fname_keep)) p_fname_keep => fname_keep
      if (allocated(fnamer))     p_fnamer => fnamer
      if (allocated(fname_sound))p_fname_sound => fname_sound
      if (allocated(fnamex))     p_fnamex => fnamex
      if (allocated(fnamey))     p_fnamey => fnamey
      if (allocated(fnamez))     p_fnamez => fnamez
      if (allocated(fnamexy))    p_fnamexy => fnamexy
      if (allocated(fnamexz))    p_fnamexz => fnamexz
      if (allocated(fnamerz))    p_fnamerz => fnamerz
      if (allocated(ncountsz))   p_ncountsz => ncountsz
      call diagnostics_init_reduc_pointers
 
    endsubroutine init_reduc_pointers
!***********************************************************************
   subroutine init_diagnostics_accumulators 
!    
!  Need to initialize accumulators since master thread does not take part in diagnostics
!
!  25-aug-23/TP: Coded
!
    use Chemistry
    if (allocated(fname))      p_fname = 0.
    if (allocated(fnamex))     p_fnamex = 0.
    if (allocated(fnamey))     p_fnamey = 0.
    if (allocated(fnamez))     p_fnamez = 0.
    if (allocated(fnamer))     p_fnamer = 0.
    if (allocated(fnamexy))    p_fnamexy = 0.
    if (allocated(fnamexz))    p_fnamexz = 0.
    if (allocated(fnamerz))    p_fnamerz = 0.
    if (allocated(fname_keep)) p_fname_keep = 0.
    if (allocated(fname_sound))p_fname_sound= 0.
    if (allocated(ncountsz))   p_ncountsz= 0

    num_of_diag_iter_done = 0
    lstarted_finalizing_diagnostics = .false.
    lfinalized_diagnostics = .false.
    lstarted_writing_diagnostics = .false.
    lwritten_diagnostics = .false.

    l1davgfirst_save = l1davgfirst
    ldiagnos_save = ldiagnos
    l1dphiavg_save = l1dphiavg
    l2davgfirst_save = l2davgfirst

    lout_save = lout
    l1davg_save = l1davg
    l2davg_save = l2davg
    lout_sound_save = lout_sound
    lvideo_save = lvideo
    lwrite_slices_save = lwrite_slices
    lchemistry_diag_save = lchemistry_diag
    it_save = it

    endsubroutine init_diagnostics_accumulators
!***********************************************************************
    subroutine diagnostics_reductions
!
!  Reduces accumulated diagnostic variables across threads. Only called if using OpenMP
!
!  30-mar-23/TP: Coded
!

!$  use OMP_lib
    use Diagnostics
    integer :: imn
      if (ldiagnos .and. allocated(fname)) then
        do imn=1,size(fname)
          if (any(inds_max_diags == imn) .and. fname(imn) /= 0.) then
            p_fname(imn) = max(p_fname(imn),fname(imn))
          else if (any(inds_sum_diags == imn)) then
            p_fname(imn) = p_fname(imn) + fname(imn)
          endif
        enddo
      endif

      if (l1davgfirst) then
        if (allocated(fnamex)) p_fnamex = p_fnamex + fnamex
        if (allocated(fnamey)) p_fnamey = p_fnamey + fnamey
        if (allocated(fnamez)) p_fnamez = p_fnamez + fnamez
        if (allocated(fnamer)) p_fnamer = p_fnamer + fnamer
      endif

      if (l2davgfirst) then
        if (allocated(fnamexy)) p_fnamexy = p_fnamexy + fnamexy
        if (allocated(fnamexz)) p_fnamexz = p_fnamexz + fnamexz
        if (allocated(fnamerz)) p_fnamerz = p_fnamerz + fnamerz
      endif

      if (allocated(fname_keep)) p_fname_keep = p_fname_keep + fname_keep
      if (allocated(fname_sound)) p_fname_sound = p_fname_sound + fname_sound
      if (allocated(ncountsz)) p_ncountsz = p_ncountsz + ncountsz

      call diagnostics_diag_reductions


    endsubroutine diagnostics_reductions
!***********************************************************************
    subroutine all_module_diags_slice(istart,iend,f,p)
!
!  Calculates module diagnostics (so far only density, energy, hydro, magnetic)
!
!  10-sep-2019/MR: coded
!
      use Ascalar, only: calc_diagnostics_ascalar
      use Chemistry, only: calc_diagnostics_chemistry
      use Chiral, only: calc_diagnostics_chiral
      use Cosmicray, only: calc_diagnostics_cosmicray
      use Density, only: calc_diagnostics_density
      use Dustdensity, only: calc_diagnostics_dustdensity
      use Dustvelocity, only: calc_diagnostics_dustvelocity
      use Energy, only: calc_diagnostics_energy
      use Forcing, only: calc_diagnostics_forcing
      use Gravity, only: calc_diagnostics_gravity
      use Heatflux, only: calc_diagnostics_heatflux
      use Hydro, only: calc_diagnostics_hydro
      use Interstellar, only: calc_diagnostics_interstellar
      use Lorenz_gauge, only: calc_diagnostics_lorenz_gauge
      use Magnetic, only: calc_diagnostics_magnetic
      use NeutralDensity, only: calc_diagnostics_neutraldens
      use NeutralVelocity, only: calc_diagnostics_neutralvel
      use Particles_main, only: particles_calc_pencil_diags
      use Pointmasses, only: calc_diagnostics_pointmasses
      use Pscalar, only: calc_diagnostics_pscalar
      use Polymer, only: calc_diagnostics_polymer
      !use Radiation, only: calc_diagnostics_radiation
      use Selfgravity, only: calc_diagnostics_selfgrav
      use Shear, only: calc_diagnostics_shear
      use Shock, only: calc_diagnostics_shock
      use Viscosity, only: calc_diagnostics_viscosity

      real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
      type (pencil_case)                ,intent(INOUT) :: p

      integer :: imn,istart,iend
!
!  This is the beginning of the famous mn-loop!
!  Here, m and n don't start with m1 and n1, as one would naively expect,
!  but they start at m1+nghost and n1+nghost, because those have not don't
!  rely on the correct calculation of any of the derivatives.
!  Once we reach m2-nghost and n2-nghost, this is the last moment
!  before we really need to make sure all the (concurrent) communication
!  has finished.
!
      lfirstpoint=.true.
      do imn=istart,iend

        n=nn(imn)
        m=mm(imn)
!
!  Skip points not belonging to coarse grid.
!
        lcoarse_mn=lcoarse.and.mexts(1)<=m.and.m<=mexts(2)
        if (lcoarse_mn) then
          lcoarse_mn=lcoarse_mn.and.ninds(0,m,n)>0
          if (ninds(0,m,n)<=0) cycle
        endif

        call calc_all_pencils(f,p)

        call calc_diagnostics_ascalar(p)
        call calc_diagnostics_chemistry(f,p)
        call calc_diagnostics_chiral(f,p)
        call calc_diagnostics_cosmicray(p)
        call calc_diagnostics_density(f,p)
        call calc_diagnostics_dustdensity(f,p)
        call calc_diagnostics_dustvelocity(p)
        call calc_diagnostics_energy(f,p)
        if (lforcing_cont) call calc_diagnostics_forcing(p)
        call calc_diagnostics_gravity(p)
        call calc_diagnostics_heatflux(p)
        call calc_diagnostics_hydro(f,p)
        call calc_diagnostics_interstellar(p)
        call calc_diagnostics_lorenz_gauge(f,p)
        call calc_diagnostics_magnetic(f,p)
        call calc_diagnostics_neutraldens(p)
        call calc_diagnostics_neutralvel(p)
        call particles_calc_pencil_diags(p)
        call calc_diagnostics_pointmasses(p)
        call calc_diagnostics_pscalar(p)
        call calc_diagnostics_polymer(p)
        !call calc_diagnostics_radiation(f,p)
        call calc_diagnostics_selfgrav(p)
        call calc_diagnostics_shear(p)
        call calc_diagnostics_shock(p)
        call calc_diagnostics_viscosity(p)

        lfirstpoint=.false.

      enddo
!$omp critical
!$   call prep_finalize_thread_diagnos
!$   if(omp_get_thread_num() /= 0) call diagnostics_reductions
      num_of_diag_iter_done = num_of_diag_iter_done + ((iend-istart)+1)
!$omp end critical

    endsubroutine all_module_diags_slice
!*****************************************************************************
    subroutine all_module_diags_threaded(f,p)
!
!
!
!$    use OMP_lib

      real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
      type (pencil_case)                ,intent(INOUT) :: p

      integer :: istart,iend,nper_thread,num_of_threads_to_use,i
!
        num_of_threads_to_use = 1
!TP: for now we are using all but one thread to do the diagnostics
!$      num_of_threads_to_use  = max(omp_get_num_procs()-1,1)
        call init_diagnostics_accumulators
        nper_thread = (nyz/num_of_threads_to_use)+1
        iend = 0
        do i=1,num_of_threads_to_use
          istart = iend+1
          iend = iend + nper_thread
!$omp task
          call all_module_diags_slice(istart,min(iend,nyz),f,p)
!$omp end task
        enddo
    endsubroutine all_module_diags_threaded
!*****************************************************************************
    subroutine calc_all_module_diagnostics(f,p)
!
!  Calculates module diagnostics (so far only density, energy, hydro, magnetic)
!
!  10-sep-2019/MR: coded
!
!$    use OMP_lib

      real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
      type (pencil_case)                ,intent(INOUT) :: p

      integer :: imn
!
!$  if(omp_in_parallel()) then
!$      call init_reduc_pointers
!$      call all_module_diags_threaded(f,p)
!$  else
        call all_module_diags_slice(1,nyz,f,p)
!$  endif
    endsubroutine calc_all_module_diagnostics
!****************************************************************************
    subroutine finalize_diagnostics
!
!  Finalizes all module diagnostics by MPI communication.
!
!  25-aug-23/TP: refactored from pde
!
      use Diagnostics
      use Hydro
      use Magnetic
      use Pscalar
!$    use OMP_lib

!  Set that the master thread knows we have started the finalization
!
!$    lstarted_finalizing_diagnostics = .true.
!
!  Restore options that were used when calc_all_module_diagnostics was called
!
!$   call read_diagnostics_accumulators
!
!TP: best would be to to pass p_fname vars as input but there are simply too
!many calls because of the calc_mfield etc.
!  0-D Diagnostics.
!
      if (lout) then
        call diagnostic(fname,nname)
        call diagnostic(fname_keep,nname,lcomplex=.true.)
      endif
!
!  1-D diagnostics.
!
      if (l1davg) then
        if (lwrite_xyaverages) call xyaverages_z(fnamez,ncountsz)
        if (lwrite_xzaverages) call xzaverages_y(fnamey)
        if (lwrite_yzaverages) call yzaverages_x(fnamex)
      endif
      if (lcylinder_in_a_box.and.l1davg) call phizaverages_r(fnamer)
!
!  2-D averages.
!
      if (l2davg) then
        if (lwrite_yaverages)   call yaverages_xz(fnamexz)
        if (lwrite_zaverages)   call zaverages_xy(fnamexy)
        if (lwrite_phiaverages) call phiaverages_rz(fnamerz)
      endif
!
!  Note: zaverages_xy are also needed if bmx and bmy are to be calculated
!  (of course, yaverages_xz does not need to be calculated for that).
!
      if (.not.l2davg.and.lout.and.ldiagnos_need_zaverages) then
        if (lwrite_zaverages) call zaverages_xy(fnamexy)
      endif
!
!  Calculate mean fields and diagnostics related to mean fields.
!
      if (lout) then
        if (lmagnetic) call calc_mfield
        if (lhydro)    call calc_mflow
        if (lpscalar)  call calc_mpscalar
      endif

!$    call write_diagnostics_accumulators
!$    lfinalized_diagnostics = .true.
 
    endsubroutine finalize_diagnostics
!****************************************************************************
    subroutine calc_all_pencils(f,p)

      use Diagnostics, only: calc_phiavg_profile

      use Ascalar
      use BorderProfiles, only: calc_pencils_borderprofiles
      use Chiral
      use Chemistry
      use Cosmicray
      use CosmicrayFlux
      use Density
      use Dustvelocity
      use Dustdensity
      use Energy
      use EquationOfState
      use Forcing, only: calc_pencils_forcing
      use Gravity
      use Heatflux
      use Hydro
      use Interstellar, only: calc_pencils_interstellar
      use Lorenz_gauge
      use Magnetic
      use NeutralDensity
      use NeutralVelocity
      use Particles_main
      use Pscalar
      use PointMasses
      use Polymer
      use Radiation
      use Selfgravity
      use Shear
      use Shock, only: calc_pencils_shock
      use Solid_Cells, only: update_solid_cells_pencil
      use Special, only: calc_pencils_special
      use Testfield
      use Testflow
      use Testscalar
      use Viscosity, only: calc_pencils_viscosity

      real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
      type (pencil_case)                ,intent(INOUT) :: p
!
!  The solid cells may have to be updated at the beginning of every
!  pencil calculation.
!
        call update_solid_cells_pencil(f)
!
!  Grid spacing. In case of equidistant grid and cartesian coordinates
!  this is calculated before the (m,n) loop.
!
        if (.not. lcartesian_coords .or. .not.all(lequidist)) call get_grid_mn
!
!  Calculate grid/geometry related pencils.
!
        call calc_pencils_grid(p)
!
!  Calculate profile for phi-averages if needed.
!
        if (((l2davgfirst.and.lwrite_phiaverages )  .or. &
             (l1dphiavg  .and.lwrite_phizaverages)) .and. &
            (lcylinder_in_a_box.or.lsphere_in_a_box)) call calc_phiavg_profile(p)
!
!  Calculate pencils for the pencil_case.
!  Note: some no-modules (e.g. nohydro) also calculate some pencils,
!  so it would be wrong to check for lhydro etc in such cases.
! DM : in the formulation of lambda effect due to Kitchanov and Olemski,
! DM : we need to have dsdr to calculate lambda. Hence, the way it is done now,
! DM : we need to have energy pencils calculated before viscosity pencils.
! DM : This is *bad* practice and must be corrected later.
!
! To check ghost cell consistency, please uncomment the following 2 lines:
!       if (.not. lpencil_check_at_work .and. necessary(imn)) &
!       call check_ghosts_consistency (f, 'before calc_pencils_*')
!
                              call calc_pencils_hydro(f,p)
                              call calc_pencils_density(f,p)
        if (lpscalar)         call calc_pencils_pscalar(f,p)
        if (lascalar)         call calc_pencils_ascalar(f,p)
                              call calc_pencils_eos(f,p)
        if (lshock)           call calc_pencils_shock(f,p)
        if (lchemistry)       call calc_pencils_chemistry(f,p)
                              call calc_pencils_energy(f,p)
        if (lviscosity)       call calc_pencils_viscosity(f,p)
        if (linterstellar)    call calc_pencils_interstellar(f,p)
        if (lforcing_cont)    call calc_pencils_forcing(f,p)
        if (llorenz_gauge)    call calc_pencils_lorenz_gauge(f,p)
        if (lmagnetic)        call calc_pencils_magnetic(f,p)
        if (lpolymer)         call calc_pencils_polymer(f,p)
        if (lgrav)            call calc_pencils_gravity(f,p)
        if (lselfgravity)     call calc_pencils_selfgravity(f,p)
        if (ldustvelocity)    call calc_pencils_dustvelocity(f,p)
        if (ldustdensity)     call calc_pencils_dustdensity(f,p)
        if (lneutralvelocity) call calc_pencils_neutralvelocity(f,p)
        if (lneutraldensity)  call calc_pencils_neutraldensity(f,p)
        if (lcosmicray)       call calc_pencils_cosmicray(f,p)
        if (lcosmicrayflux)   call calc_pencils_cosmicrayflux(f,p)
        if (lchiral)          call calc_pencils_chiral(f,p)
        if (lradiation)       call calc_pencils_radiation(f,p)
        if (lshear)           call calc_pencils_shear(f,p)
        if (lborder_profiles) call calc_pencils_borderprofiles(f,p)
        if (lpointmasses)     call calc_pencils_pointmasses(f,p)
        if (lparticles)       call particles_calc_pencils(f,p)
        if (lspecial)         call calc_pencils_special(f,p)
        if (lheatflux)        call calc_pencils_heatflux(f,p)

    endsubroutine calc_all_pencils
!***********************************************************************
    subroutine rhs_cpu(f,df,p,mass_per_proc,early_finalize)
!
!  Calculates rhss of the PDEs.
!
!  14-feb-17/MR: Carved out from pde.
!  21-feb-17/MR: Moved all module-specific estimators of the (inverse) possible timestep 
!                to the individual modules.
!
      use Ascalar
      use Chiral
      use Chemistry
      use Cosmicray
      use CosmicrayFlux
      use Density
      use Diagnostics
      use Dustvelocity
      use Dustdensity
      use Energy
      use EquationOfState
      use Forcing, only: calc_diagnostics_forcing
      use Gravity
      use Heatflux
      use Hydro
      use Lorenz_gauge
      use Magnetic
      use NeutralDensity
      use NeutralVelocity
      use Particles_main
      use Pscalar
      use PointMasses
      use Polymer
      use Radiation
      use Selfgravity
      use Shear
      use Solid_Cells, only: dsolid_dt
      use Special, only: dspecial_dt
      use Sub, only: sum_mn
      use Testfield
      use Testflow
      use Testscalar

      real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
      real, dimension (mx,my,mz,mvar)   ,intent(OUT  ) :: df
      type (pencil_case)                ,intent(INOUT) :: p
      real, dimension(1)                ,intent(INOUT) :: mass_per_proc
      logical                           ,intent(IN   ) :: early_finalize

      real, dimension (nx,3) :: df_iuu_pencil
      logical :: lcommunicate
!
      lfirstpoint=.true.
      lcommunicate=.not.early_finalize

      mn_loop: do imn=1,nyz

        n=nn(imn)
        m=mm(imn)

        !if (imn_array(m,n)==0) cycle
!
!  Skip points not belonging to coarse grid.
!
        lcoarse_mn=lcoarse.and.mexts(1)<=m.and.m<=mexts(2)
        if (lcoarse_mn) then
          lcoarse_mn=lcoarse_mn.and.ninds(0,m,n)>0
          if (ninds(0,m,n)<=0) cycle
        endif
!
!  Store the velocity part of df array in a temporary array
!  while solving the anelastic case.
!
        call timing('pde','before lanelastic',mnloop=.true.)

        if (lanelastic) then
          df_iuu_pencil = df(l1:l2,m,n,iux:iuz)
          df(l1:l2,m,n,iux:iuz)=0.0
        endif
!
!        if (loptimise_ders) der_call_count=0 !DERCOUNT
!
!  Make sure all ghost points are set.
!
        if (lcommunicate) then
          if (necessary(imn)) then
            call finalize_isendrcv_bdry(f)
            call boundconds_y(f)
            call boundconds_z(f)
            lcommunicate=.false.
          endif
        endif
        call timing('pde','finished boundconds_z',mnloop=.true.)
!
!  For each pencil, accumulate through the different modules
!  advec_XX and diffus_XX, which are essentially the inverse
!  advective and diffusive timestep for that module.
!  (note: advec2 is an inverse _squared_ timestep)
!  Note that advec_cs2 should always be initialized when leos.
!
        if (lfirst.and.ldt.and.(.not.ldt_paronly)) then
          advec_cs2=0.0
          maxadvec=0.
          if (lenergy.or.ldensity.or.lmagnetic.or.lradiation.or.lneutralvelocity.or.lcosmicray.or. &
              (ltestfield_z.and.iuutest>0)) &
            advec2=0.
          if (ldensity.or.lviscosity.or.lmagnetic.or.lenergy.or.ldustvelocity.or.ldustdensity) &
            advec2_hypermesh=0.0
          maxdiffus=0.
          maxdiffus2=0.
          maxdiffus3=0.
          maxsrc=0.
        endif

        call calc_all_pencils(f,p)
!
!  --------------------------------------------------------
!  NO CALLS MODIFYING PENCIL_CASE PENCILS BEYOND THIS POINT
!  --------------------------------------------------------
!
!  hydro, density, and entropy evolution
!  Note that pressure gradient is added in denergy_dt of noentropy to momentum,
!  even if lentropy=.false.
!
        call duu_dt(f,df,p)
        call dlnrho_dt(f,df,p)
        call denergy_dt(f,df,p)
!
!  Magnetic field evolution
!
        if (lmagnetic) call daa_dt(f,df,p)
!
!  Lorenz gauge evolution
        if (llorenz_gauge) call dlorenz_gauge_dt(f,df,p)
!
!  Polymer evolution
!
        if (lpolymer) call dpoly_dt(f,df,p)
!
!  Testscalar evolution
!
        if (ltestscalar) call dcctest_dt(f,df,p)
!
!  Testfield evolution
!
        if (ltestfield) call daatest_dt(f,df,p)
!
!  Testflow evolution
!
        if (ltestflow) call duutest_dt(f,df,p)
!
!  Passive scalar evolution
!
        if (lpscalar) call dlncc_dt(f,df,p)
!
!  Supersaturation evolution
        
        if (lascalar) call dacc_dt(f,df,p)
!
!  Dust evolution
!
        if (ldustvelocity) call duud_dt(f,df,p)
        if (ldustdensity) call dndmd_dt(f,df,p)
!
!  Neutral evolution
!
        if (lneutraldensity) call dlnrhon_dt(f,df,p)
        if (lneutralvelocity) call duun_dt(f,df,p)
!
!  Add gravity, if present
!
        if (lgrav) then
          if (lhydro.or.ldustvelocity.or.lneutralvelocity) &
               call duu_dt_grav(f,df,p)
        endif
!
!  Self-gravity
!
        if (lselfgravity) call duu_dt_selfgrav(f,df,p)
!
!  Cosmic ray energy density
!
        if (lcosmicray) call decr_dt(f,df,p)
!
!  Cosmic ray flux
!
        if (lcosmicrayflux) call dfcr_dt(f,df,p)
!
!  Chirality of left and right handed aminoacids
!
        if (lchiral) call dXY_chiral_dt(f,df,p)
!
!  Evolution of chemical species
!
        if (lchemistry) call dchemistry_dt(f,df,p)
!
!  Evolution of heatflux vector
!
        if (lheatflux) call dheatflux_dt(f,df,p)
!
!  Continuous forcing diagonstics.
!
        if (lforcing_cont) call calc_diagnostics_forcing(p)
!
!  Add and extra 'special' physics
!
        if (lspecial) call dspecial_dt(f,df,p)
!
!  Add radiative cooling and radiative pressure (for ray method)
!
        if (lradiation_ray.and.lenergy) call dradiation_dt(f,df,p)
!
!  Find diagnostics related to solid cells (e.g. drag and lift).
!  Integrating to the full result is done after loops over m and n.
!
        if (lsolid_cells) call dsolid_dt(f,df,p)
!
!  Add shear if present
!
        if (lshear) call shearing(f,df,p)
!
        if (lparticles) call particles_pde_pencil(f,df,p)
!
        if (lpointmasses) call pointmasses_pde_pencil(f,df,p)
!
!  Call diagnostics that involves the full right hand side
!  This must be done at the end of all calls that might modify df.
!
        if (ldiagnos) then
          if (lhydro) call df_diagnos_hydro(df,p)
          if (lmagnetic) call df_diagnos_magnetic(df,p)
        endif
!
!  General phiaverage quantities -- useful for debugging.
!  MR: Results are constant in time, so why here?
!
        if (l2davgfirst) then
          call phisum_mn_name_rz(p%rcyl_mn,idiag_rcylmphi)
          call phisum_mn_name_rz(p%phi_mn,idiag_phimphi)
          call phisum_mn_name_rz(p%z_mn,idiag_zmphi)
          call phisum_mn_name_rz(p%r_mn,idiag_rmphi)
        endif
!
!  Do the time integrations here, before the pencils are overwritten.
!
        if (ltime_integrals.and.llast) then
          if (lhydro)    call time_integrals_hydro(f,p)
          if (lmagnetic) call time_integrals_magnetic(f,p)
        endif

        call set_dt1_max(p)
!
!  Diagnostics showing how close to advective and diffusive time steps we are
!
        if (lfirst.and.ldt.and.(.not.ldt_paronly).and.ldiagnos) then
          if (idiag_dtv/=0) call max_mn_name(maxadvec/cdt,idiag_dtv,l_dt=.true.)
          if (idiag_dtdiffus/=0) call max_mn_name(maxdiffus/cdtv,idiag_dtdiffus,l_dt=.true.)
          if (idiag_dtdiffus2/=0) call max_mn_name(maxdiffus2/cdtv2,idiag_dtdiffus2,l_dt=.true.)
          if (idiag_dtdiffus3/=0) call max_mn_name(maxdiffus3/cdtv3,idiag_dtdiffus3,l_dt=.true.)
!
!  Regular and hyperdiffusive mesh Reynolds numbers
!
          if (idiag_Rmesh/=0) call max_mn_name(pi_1*maxadvec/(maxdiffus+tini),idiag_Rmesh)
          if (idiag_Rmesh3/=0) call max_mn_name(pi5_1*maxadvec/(maxdiffus3+tini),idiag_Rmesh3)
          call max_mn_name(maxadvec,idiag_maxadvec)
        endif
!
!  Display derivative info
!
!debug   if (loptimise_ders.and.lout) then                         !DERCOUNT
!debug     do iv=1,nvar                                            !DERCOUNT
!debug     do ider=1,8                                             !DERCOUNT
!debug     do j=1,3                                                !DERCOUNT
!debug     do k=1,3                                                !DERCOUNT
!debug       if (der_call_count(iv,ider,j,k) > 1) then             !DERCOUNT
!debug         print*,'DERCOUNT: '//varname(iv)//' derivative ', & !DERCOUNT
!debug                                                 ider,j,k, & !DERCOUNT
!debug                                               ' called ', & !DERCOUNT
!debug                              der_call_count(iv,ider,j,k), & !DERCOUNT
!debug                                                  'times!'   !DERCOUNT
!debug       endif                                                 !DERCOUNT
!debug     enddo                                                   !DERCOUNT
!debug     enddo                                                   !DERCOUNT
!debug     enddo                                                   !DERCOUNT
!debug     enddo                                                   !DERCOUNT
!debug     if (maxval(der_call_count)>1) call fatal_error( &        !DERCOUNT
!debug      'pde','ONE OR MORE DERIVATIVES HAS BEEN DOUBLE CALLED') !DERCOUNT
!debug   endif
!
!  Fill in the rhs of the poisson equation and restore the df(:,:,:,iuu) array
!  for anelastic case
!
        if (lanelastic) then
!          call calc_pencils_density(f,p)
          f(l1:l2,m,n,irhs)   = p%rho*df(l1:l2,m,n,iuu)
          f(l1:l2,m,n,irhs+1) = p%rho*df(l1:l2,m,n,iuu+1)
          f(l1:l2,m,n,irhs+2) = p%rho*df(l1:l2,m,n,iuu+2)
          df(l1:l2,m,n,iux:iuz) = df_iuu_pencil + df(l1:l2,m,n,iux:iuz)
          call sum_mn(p%rho,mass_per_proc(1))
        endif
        call timing('rhs_cpu','end of mn loop',mnloop=.true.)
!
!  End of loops over m and n.
!
        headtt=.false.
        lfirstpoint=.false.
!
      enddo mn_loop
!
!!$     call prep_finalize_thread_diagnos
!
      if (ltime_integrals.and.llast) then
        if (lhydro) call update_for_time_integrals_hydro
      endif
!
    endsubroutine rhs_cpu
!***********************************************************************
    subroutine debug_imn_arrays
!
!  For debug purposes: writes out the mm, nn, and necessary arrays.
!
!  23-nov-02/axel: coded
!
      open(1,file=trim(directory)//'/imn_arrays.dat')
      do imn=1,nyz
        if (necessary(imn)) write(1,'(a)') '----necessary=.true.----'
        write(1,'(4i6)') imn,mm(imn),nn(imn)
      enddo
      close(1)
!
    endsubroutine debug_imn_arrays
!***********************************************************************
    subroutine output_crash_files(f)
!
!  Write crash snapshots when time-step is low.
!
!  15-aug-2007/anders: coded
!
      use Snapshot
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      integer, save :: icrash=0
      character (len=10) :: filename
      character (len=1) :: icrash_string
!
      if ( (it>1) .and. lfirst .and. (dt<=crash_file_dtmin_factor*dtmin) ) then
        write(icrash_string, fmt='(i1)') icrash
        filename='crash'//icrash_string//'.dat'
        call wsnap(filename,f,mvar_io,ENUM=.false.)
        if (lroot) then
          print*, 'Time-step is very low - writing '//trim(filename)
          print*, '(it, itsub=', it, itsub, ')'
          print*, '(t, dt=', t, dt, ')'
        endif
!
!  Next crash index, cycling from 0-9 to avoid excessive writing of
!  snapshots to the hard disc.
!
        icrash=icrash+1
        icrash=mod(icrash,10)
      endif
!
    endsubroutine output_crash_files
!***********************************************************************
    subroutine set_dyndiff_coeff(f)
!
!  Set dynamical diffusion coefficients.
!
!  18-jul-14/ccyang: coded.
!  03-apr-16/ccyang: add switch ldyndiff_useumax
!
      use Density,   only: dynamical_diffusion
      use Energy,    only: dynamical_thermal_diffusion
      use Magnetic,  only: dynamical_resistivity
      use Sub,       only: find_max_fvec, find_rms_fvec
      use Viscosity, only: dynamical_viscosity
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      real :: uc
!
!  Find the characteristic speed, either the absolute maximum or the rms.
!
      if (ldyndiff_useumax) then
        uc = find_max_fvec(f, iuu)
      else
        uc = find_rms_fvec(f, iuu)
      endif
!
!  Ask each module to set the diffusion coefficients.
!
      if (ldensity)                      call dynamical_diffusion(uc)
      if (lmagnetic .and. .not. lbfield) call dynamical_resistivity(uc)
      if (lenergy)                       call dynamical_thermal_diffusion(uc)
      if (lviscosity)                    call dynamical_viscosity(uc)
!
    endsubroutine set_dyndiff_coeff
!***********************************************************************
    subroutine impose_floors_ceilings(f)
!
!  Impose floors or ceilings for implemented fields.
!
!  20-oct-14/ccyang: modularized from pde.
!
      use Cosmicray, only: impose_ecr_floor
      use Density, only: impose_density_floor,impose_density_ceiling
      use Dustdensity, only: impose_dustdensity_floor
      use Energy, only: impose_energy_floor
      use Hydro, only: impose_velocity_ceiling
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call impose_density_floor(f)
      call impose_density_ceiling(f)
      call impose_velocity_ceiling(f)
      call impose_energy_floor(f)
      call impose_dustdensity_floor(f)
      call impose_ecr_floor(f)
!
    endsubroutine impose_floors_ceilings
!***********************************************************************
    subroutine freeze(df,p)
!
      use Sub, only: quintic_step
      use Solid_Cells, only: freeze_solid_cells

      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case) :: p

      logical, dimension(npencils) :: lpenc_loc
      real, dimension(nx) :: pfreeze
      integer :: imn,iv
!
      !if (lgpu) then
      !  call freeze_gpu
      !  return
      !endif  

      lpenc_loc=.false. 
      if (lcylinder_in_a_box.or.lcylindrical_coords) then 
        lpenc_loc(i_rcyl_mn)=.true.
      else
        lpenc_loc(i_r_mn)=.true.
      endif
!
      headtt = headt .and. lfirst .and. lroot
!
!!$omp do private(p,pfreeze,iv)
!
      do imn=1,nyz

        n=nn(imn)
        m=mm(imn)
        !if (imn_array(m,n)==0) cycle
!
!  Recalculate grid/geometry related pencils. The r_mn and rcyl_mn are requested
!  in pencil_criteria_grid. Unfortunately we need to recalculate them here.
!
        if (any(lfreeze_varext).or.any(lfreeze_varint)) then

          call calc_pencils_grid(p,lpenc_loc)
!
!  Set df=0 for r_mn<rfreeze_int.
!
          if (any(lfreeze_varint)) then
            if (headtt) print*, 'pde: freezing variables for r < ', rfreeze_int, &
                                ' : ', lfreeze_varint
            if (lcylinder_in_a_box.or.lcylindrical_coords) then
              if (wfreeze_int==0.0) then
                where (p%rcyl_mn<=rfreeze_int)
                  pfreeze=0.0
                elsewhere  
                  pfreeze=1.0
                endwhere  
              else
                pfreeze=quintic_step(p%rcyl_mn,rfreeze_int,wfreeze_int,SHIFT=fshift_int)
              endif
            else
              if (wfreeze_int==0.0) then
                where (p%r_mn<=rfreeze_int)
                  pfreeze=0.0
                elsewhere 
                  pfreeze=1.0
                endwhere  
              else
                pfreeze=quintic_step(p%r_mn,rfreeze_int,wfreeze_int,SHIFT=fshift_int)
              endif
            endif
!
            do iv=1,nvar
              if (lfreeze_varint(iv)) df(l1:l2,m,n,iv)=pfreeze*df(l1:l2,m,n,iv)
            enddo
!
          endif
!
!  Set df=0 for r_mn>r_ext.
!
          if (any(lfreeze_varext)) then
            if (headtt) print*, 'pde: freezing variables for r > ', rfreeze_ext, &
                                ' : ', lfreeze_varext
            if (lcylinder_in_a_box.or.lcylindrical_coords) then
              if (wfreeze_ext==0.0) then
                where (p%rcyl_mn>=rfreeze_ext)
                  pfreeze=0.0
                elsewhere  
                  pfreeze=1.0
                endwhere  
              else
                pfreeze=1.0-quintic_step(p%rcyl_mn,rfreeze_ext,wfreeze_ext,SHIFT=fshift_ext)
              endif
            else
              if (wfreeze_ext==0.0) then
                where (p%r_mn>=rfreeze_ext)
                  pfreeze=0.0
                elsewhere 
                  pfreeze=1.0
                endwhere  
              else
                pfreeze=1.0-quintic_step(p%r_mn,rfreeze_ext,wfreeze_ext,SHIFT=fshift_ext)
              endif
            endif
!
            do iv=1,nvar
              if (lfreeze_varext(iv)) df(l1:l2,m,n,iv) = pfreeze*df(l1:l2,m,n,iv)
            enddo
          endif
        endif
!
!  Set df=0 inside square.
!
        if (any(lfreeze_varsquare)) then
          if (headtt) print*, 'pde: freezing variables inside square : ',lfreeze_varsquare
          pfreeze=1.0-quintic_step(x(l1:l2),xfreeze_square,wfreeze,SHIFT=-1.0)*&
                      quintic_step(spread(y(m),1,nx),yfreeze_square,-wfreeze,SHIFT=-1.0)
!
          do iv=1,nvar
            if (lfreeze_varsquare(iv)) df(l1:l2,m,n,iv) = pfreeze*df(l1:l2,m,n,iv)
          enddo
        endif
!
!  Freeze components of variables in boundary slice if specified by boundary
!  condition 'f'
!
!  Freezing boundary conditions in x.
!
        if (lfrozen_bcs_x) then ! are there any frozen vars at all?
!
!  Only need to do this on left/right-most
!  processor and in left/right--most pencils
!
          if (lfirst_proc_x) where (lfrozen_bot_var_x(1:nvar)) df(l1,m,n,1:nvar) = 0.
          if (llast_proc_x) where (lfrozen_top_var_x(1:nvar)) df(l2,m,n,1:nvar) = 0.
!
        endif
!
!  Freezing boundary conditions in y.
!
        if (lfrozen_bcs_y) then ! are there any frozen vars at all?
!
!  Only need to do this on bottom/top-most
!  processor and in bottom/top-most pencils.
!
          if (lfirst_proc_y .and. (m == m1)) then
            do iv=1,nvar
              if (lfrozen_bot_var_y(iv)) df(l1:l2,m,n,iv) = 0.
            enddo
          endif
          if (llast_proc_y .and. (m == m2)) then
            do iv=1,nvar
              if (lfrozen_top_var_y(iv)) df(l1:l2,m,n,iv) = 0.
            enddo
          endif
        endif
!
!  Freezing boundary conditions in z.
!
        if (lfrozen_bcs_z) then ! are there any frozen vars at all?
!
!  Only need to do this on bottom/top-most
!  processor and in bottom/top-most pencils.
!
          if (lfirst_proc_z .and. (n == n1)) then
            do iv=1,nvar
              if (lfrozen_bot_var_z(iv)) df(l1:l2,m,n,iv) = 0.
            enddo
          endif
          if (llast_proc_z .and. (n == n2)) then
            do iv=1,nvar
              if (lfrozen_top_var_z(iv)) df(l1:l2,m,n,iv) = 0.
            enddo
          endif
        endif
!
!  Set df=0 for all solid cells.
!
        call freeze_solid_cells(df)
!
        headtt=.false.
!
      enddo
!
!!$omp end do
!
    endsubroutine freeze
!***********************************************************************
    subroutine set_dt1_max(p)
!
!  In max_mn maximum values of u^2 (etc) are determined sucessively
!  va2 is set in magnetic (or nomagnetic)
!  In rms_mn sum of all u^2 (etc) is accumulated
!  Calculate maximum advection speed for timestep; needs to be done at
!  the first substep of each time step
!  Note that we are (currently) accumulating the maximum value,
!  not the maximum squared!
!
!  The dimension of the run ndim (=0, 1, 2, or 3) enters the viscous time step.
!  This has to do with the term on the diagonal, cdtv depends on order of scheme
!
      use General, only: notanumber

      type (pencil_case), intent(IN) :: p

      real, dimension(nx) :: dt1_max_loc, dt1_advec, dt1_diffus, dt1_src
      real :: dt1_preac

      if (lfirst.and.ldt.and.(.not.ldt_paronly)) then
!
!  sum or maximum of the advection terms?
!  (lmaxadvec_sum=.false. by default)
!
        advec2=advec2+advec_cs2
        if (lenergy.or.ldensity.or.lmagnetic.or.lradiation.or.lneutralvelocity.or.lcosmicray.or. &
            (ltestfield_z.and.iuutest>0)) maxadvec=maxadvec+sqrt(advec2)

        if (ldensity.or.lviscosity.or.lmagnetic.or.lenergy.or.ldustvelocity.or.ldustdensity) &
            maxadvec=maxadvec+sqrt(advec2_hypermesh)
!
!  Time step constraints from each module. (At the moment, magnetic and testfield use the same variable.)
!  cdt, cdtv, and cdtc are empirical non-dimensional coefficients.
!
!  Timestep constraint from advective terms.
!
        dt1_advec = maxadvec/cdt
!
!  Timestep constraint from diffusive terms.
!
        dt1_diffus = maxdiffus/cdtv + maxdiffus2/cdtv2 + maxdiffus3/cdtv3
!
!  Timestep constraint from source terms.
!
        dt1_src = maxsrc/cdtsrc
!
!  Timestep combination from advection, diffusion and "source". 
!
        dt1_max_loc = sqrt(dt1_advec**2 + dt1_diffus**2 + dt1_src**2)
!
!  time step constraint from the coagulation kernel
!
        if (ldustdensity) dt1_max_loc = max(dt1_max_loc,reac_dust/cdtc)
!
!  time step constraint from speed of chemical reactions
!
        if (lchemistry .and. .not.llsode) then
!           dt1_preac= reac_pchem/cdtc
          dt1_max_loc = max(dt1_max_loc,reac_chem/cdtc)
        endif
!
!  time step constraint from relaxation time of polymer
!
        if (lpolymer) dt1_max_loc = max(dt1_max_loc,1./(trelax_poly*cdt_poly))
!
!  Exclude the frozen zones from the time-step calculation.
!
        if (any(lfreeze_varint)) then
          if (lcylinder_in_a_box.or.lcylindrical_coords) then
            where (p%rcyl_mn<=rfreeze_int) 
              dt1_max_loc=0.; maxadvec=0.; maxdiffus=0.; maxdiffus2=0.; maxdiffus3=0.
            endwhere
          else
            where (p%r_mn<=rfreeze_int) 
              dt1_max_loc=0.; maxadvec=0.; maxdiffus=0.; maxdiffus2=0.; maxdiffus3=0.
            endwhere
          endif
        endif
!
        if (any(lfreeze_varext)) then
          if (lcylinder_in_a_box.or.lcylindrical_coords) then
            where (p%rcyl_mn>=rfreeze_ext)
              dt1_max_loc=0.; maxadvec=0.; maxdiffus=0.; maxdiffus2=0.; maxdiffus3=0.
            endwhere
          else
            where (p%r_mn>=rfreeze_ext)
              dt1_max_loc=0.; maxadvec=0.; maxdiffus=0.; maxdiffus2=0.; maxdiffus3=0.
            endwhere
          endif
        endif
!  MR: the next correct? freezes *outside* square
        if (any(lfreeze_varsquare).and.y(m)>yfreeze_square) then
          where (x(l1:l2)>xfreeze_square)
            dt1_max_loc=0.; maxadvec=0.; maxdiffus=0.; maxdiffus2=0.; maxdiffus3=0.
          endwhere
        endif

        dt1_max=max(dt1_max,dt1_max_loc)
!
!  Check for NaNs in the advection time-step.
!
        if (notanumber(maxadvec)) then
          print*, 'pde: maxadvec contains a NaN at iproc=', iproc_world
          if (lenergy) print*, 'advec_cs2  =',advec_cs2
          call fatal_error_local('set_dt1_max','')
        endif
      endif

    endsubroutine set_dt1_max
!***********************************************************************
    subroutine test_dt(f,df,p,rhs_1,rhs_2)

      use Mpicomm
      use Ascalar
      use Chiral
      use Chemistry
      use Cosmicray
      use CosmicrayFlux
      use Density
      use Diagnostics
      use Dustvelocity
      use Dustdensity
      use Energy
      use EquationOfState
      use Forcing, only: calc_diagnostics_forcing
      use Gravity
      use Heatflux
      use Hydro
      use Lorenz_gauge
      use Magnetic
      use NeutralDensity
      use NeutralVelocity
      use Particles_main
      use Pscalar
      use PointMasses
      use Polymer
      use Radiation
      use Selfgravity
      use Shear
      use Solid_Cells, only: dsolid_dt
      use Special, only: dspecial_dt
      use Sub, only: sum_mn
      use Testfield
      use Testflow
      use Testscalar

      real, dimension (mx,my,mz,mfarray) :: f,f_copy
      real, dimension (mx,my,mz,mfarray) :: df,df_copy
      integer :: i,j,k,n
      type (pencil_case) :: p,p_copy

      intent(inout) :: f
      intent(in) :: p
      intent(inout) :: df
      logical :: passed
      interface
          subroutine rhs_1(f,df,p)
              import mx
              import my
              import mz
              import mfarray
              import pencil_case
              real, dimension (mx,my,mz,mfarray) :: f
              real, dimension (mx,my,mz,mfarray) :: df
              type (pencil_case) :: p

              intent(inout) :: f
              intent(in) :: p
              intent(inout) :: df
          endsubroutine rhs_1 
        endinterface
        interface
          subroutine rhs_2(f,df,p)
              import mx
              import my
              import mz
              import mfarray
              import pencil_case
              real, dimension (mx,my,mz,mfarray) :: f
              real, dimension (mx,my,mz,mfarray) :: df
              type (pencil_case) :: p

              intent(inout) :: f
              intent(in) :: p
              intent(inout) :: df
          endsubroutine rhs_2 
      endinterface
      df_copy = df
      p_copy = p
      f_copy = f
      call rhs_1(f,df,p)
      call rhs_2(f_copy,df_copy,p_copy)
      passed = .true.
      do i=1,mx
        do j=1,my
          do k=1,mz
            do n=1,mfarray
              if(df_copy(i,j,k,n) /= df(i,j,k,n)) then
                print*,"Wrong at: ",i,j,k,n
                print*,"diff",df_copy(i,j,k,n) - df(i,j,k,n)
                passed = .false.
              endif
            enddo
          enddo
        enddo
      enddo
      if(passed) then
        print*,"passed test :)"
      else
        print*,"did not pass test :/"
      endif
      print*,iux,iuy,iuz,iss,ilnrho
      call die_gracefully
    endsubroutine test_dt
!***********************************************************************
subroutine test_rhs(f,df,p,mass_per_proc,early_finalize,rhs_1,rhs_2)

      use Mpicomm
      use Ascalar
      use Chiral
      use Chemistry
      use Cosmicray
      use CosmicrayFlux
      use Density
      use Diagnostics
      use Dustvelocity
      use Dustdensity
      use Energy
      use EquationOfState
      use Forcing, only: calc_diagnostics_forcing
      use Gravity
      use Heatflux
      use Hydro
      use Lorenz_gauge
      use Magnetic
      use NeutralDensity
      use NeutralVelocity
      use Particles_main
      use Pscalar
      use PointMasses
      use Polymer
      use Radiation
      use Selfgravity
      use Shear
      use Solid_Cells, only: dsolid_dt
      use Special, only: dspecial_dt
      use Sub, only: sum_mn
      use Testfield
      use Testflow
      use Testscalar

      real, dimension (mx,my,mz,mfarray) :: f,f_copy
      real, dimension (mx,my,mz,mfarray) :: df,df_copy
      type (pencil_case) :: p,p_copy
      real, dimension(1), intent(inout) :: mass_per_proc
      logical ,intent(in) :: early_finalize
      integer :: i,j,k,n
      logical :: passed
      interface
          subroutine rhs_1(f,df,p,mass_per_proc,early_finalize)
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
          endsubroutine rhs_1 
        endinterface
        interface
          subroutine rhs_2(f,df,p,mass_per_proc,early_finalize)
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
          endsubroutine rhs_2 
      endinterface
      df_copy = df
      p_copy = p
      f_copy = f
      call rhs_1(f,df,p,mass_per_proc,early_finalize)
      call rhs_2(f_copy,df_copy,p_copy,mass_per_proc,early_finalize)
      passed = .true.
      do i=1,mx
        do j=1,my
          do k=1,mz
            do n=1,mfarray
              if(df_copy(i,j,k,n) /= df(i,j,k,n)) then
                print*,"Wrong at: ",i,j,k,n
                print*,"diff",df_copy(i,j,k,n) - df(i,j,k,n)
                passed = .false.
              endif
            enddo
          enddo
        enddo
      enddo
      if(passed) then
        print*,"passed test :)"
      else
        print*,"did not pass test :/"
      endif
      print*,iux,iuy,iuz,iss,ilnrho
      call die_gracefully
    endsubroutine test_rhs
!***********************************************************************
subroutine rhs_cpu_test(f,df,p,mass_per_proc,early_finalize)
use Deriv
use General
use Energy
use Viscosity
use Hydro
use Density
use Gravity
use EquationOfState
real, dimension (mx,my,mz,mfarray),intent(inout) :: f
real, dimension (mx,my,mz,mfarray)   ,intent(out  ) :: df
type (pencil_case)              ,intent(inout) :: p
real, dimension(1)                ,intent(inout) :: mass_per_proc
logical                           ,intent(in   ) :: early_finalize
real, dimension (nxgrid/nprocx,3) :: df_iuu_pencil
real, dimension (nxgrid/nprocx)::tmp_53_54_55_82
real, dimension (nxgrid/nprocx)::dd_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_rho_53_54_55_82
real, dimension (nxgrid/nprocx,3) :: tmp3_53_54_55_82
real, dimension (nxgrid/nprocx,3,3) :: tmp33_53_54_55_82
real::cs201_53_54_55_82
real::outest_53_54_55_82
integer::i_53_54_55_82
integer::j_53_54_55_82
integer::ju_53_54_55_82
integer::jj_53_54_55_82
integer::kk_53_54_55_82
integer::jk_53_54_55_82
real, dimension (nxgrid/nprocx)::a_max_4_53_54_55_82
logical::fast_sqrt1_4_53_54_55_82
logical::precise_sqrt1_4_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_5_53_54_55_82
integer::i_5_53_54_55_82
integer::j_5_53_54_55_82
integer::k1_5_53_54_55_82
integer::i_8_53_54_55_82
integer::j_8_53_54_55_82
logical :: lshear_ros_8_53_54_55_82
logical :: loptest_return_value_7_8_53_54_55_82
integer::i_9_53_54_55_82
integer::j_9_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_10_53_54_55_82
integer::i_10_53_54_55_82
integer::j_10_53_54_55_82
integer::k1_10_53_54_55_82
integer::a1_12_53_54_55_82
integer::a2_12_53_54_55_82
logical :: loptest_return_value_11_12_53_54_55_82
real, dimension (nxgrid/nprocx)::a_max_13_53_54_55_82
logical::fast_sqrt1_13_53_54_55_82
logical::precise_sqrt1_13_53_54_55_82
integer :: i_15_53_54_55_82
logical :: loptest_return_value_14_15_53_54_55_82
real, dimension (nxgrid/nprocx)::a_max_17_53_54_55_82
logical::fast_sqrt1_17_53_54_55_82
logical::precise_sqrt1_17_53_54_55_82
real, dimension (nxgrid/nprocx,3)::ff_25_53_54_55_82
real, dimension (nxgrid/nprocx,3)::grad_f_tmp_25_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_25_53_54_55_82
integer::j_25_53_54_55_82
integer :: i_18_23_25_53_54_55_82
logical :: loptest_return_value_14_18_23_25_53_54_55_82
logical :: loptest_return_value_19_23_25_53_54_55_82
real, dimension(nxgrid/nprocx) :: del6f_upwind_22_23_25_53_54_55_82
integer :: msk_22_23_25_53_54_55_82
real, dimension(nxgrid/nprocx,3)::del6f_21_22_23_25_53_54_55_82
integer, dimension(nxgrid/nprocx)            :: indxs_21_22_23_25_53_54_55_82
integer::j_21_22_23_25_53_54_55_82
integer::msk_21_22_23_25_53_54_55_82
logical :: loptest_return_value_24_25_53_54_55_82
real, dimension (nxgrid/nprocx)::a_max_26_53_54_55_82
logical::fast_sqrt1_26_53_54_55_82
logical::precise_sqrt1_26_53_54_55_82
real, dimension (nxgrid/nprocx,3)::ff_27_53_54_55_82
real, dimension (nxgrid/nprocx,3)::grad_f_tmp_27_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_27_53_54_55_82
integer::j_27_53_54_55_82
integer :: i_18_23_27_53_54_55_82
logical :: loptest_return_value_14_18_23_27_53_54_55_82
logical :: loptest_return_value_19_23_27_53_54_55_82
real, dimension(nxgrid/nprocx) :: del6f_upwind_22_23_27_53_54_55_82
integer :: msk_22_23_27_53_54_55_82
real, dimension(nxgrid/nprocx,3)::del6f_21_22_23_27_53_54_55_82
integer, dimension(nxgrid/nprocx)            :: indxs_21_22_23_27_53_54_55_82
integer::j_21_22_23_27_53_54_55_82
integer::msk_21_22_23_27_53_54_55_82
logical :: loptest_return_value_24_27_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_29_53_54_55_82
integer::i_29_53_54_55_82
integer::k1_29_53_54_55_82
real, dimension (nxgrid/nprocx)::d4fdx_28_29_53_54_55_82
real, dimension (nxgrid/nprocx)::d4fdy_28_29_53_54_55_82
real, dimension (nxgrid/nprocx)::d4fdz_28_29_53_54_55_82
logical :: ignore_dx_28_29_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_33_53_54_55_82
integer::i_33_53_54_55_82
integer::k1_33_53_54_55_82
logical :: loptest_return_value_30_33_53_54_55_82
real, dimension(nxgrid/nprocx)::tmp_31_33_53_54_55_82
integer::i_31_33_53_54_55_82
integer::j_31_33_53_54_55_82
real, dimension (nxgrid/nprocx)::d6fdx_32_33_53_54_55_82
real, dimension (nxgrid/nprocx)::d6fdy_32_33_53_54_55_82
real, dimension (nxgrid/nprocx)::d6fdz_32_33_53_54_55_82
logical :: ignore_dx_32_33_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_34_53_54_55_82
integer::i_34_53_54_55_82
integer::k1_34_53_54_55_82
logical :: loptest_return_value_30_34_53_54_55_82
real, dimension(nxgrid/nprocx)::tmp_31_34_53_54_55_82
integer::i_31_34_53_54_55_82
integer::j_31_34_53_54_55_82
real, dimension (nxgrid/nprocx)::d6fdx_32_34_53_54_55_82
real, dimension (nxgrid/nprocx)::d6fdy_32_34_53_54_55_82
real, dimension (nxgrid/nprocx)::d6fdz_32_34_53_54_55_82
logical :: ignore_dx_32_34_53_54_55_82
real, dimension(nxgrid/nprocx) :: tmp_35_53_54_55_82
integer::ki_35_53_54_55_82
integer::kj_35_53_54_55_82
integer::i_35_53_54_55_82
integer::j_35_53_54_55_82
integer::k_35_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fjji_36_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fijj_36_53_54_55_82
real, dimension (nxgrid/nprocx,3) ::  fjik_36_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_36_53_54_55_82
integer::i_36_53_54_55_82
integer::j_36_53_54_55_82
integer::k1_36_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fjji_37_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fijj_37_53_54_55_82
real, dimension (nxgrid/nprocx,3) ::  fjik_37_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_37_53_54_55_82
integer::i_37_53_54_55_82
integer::j_37_53_54_55_82
integer::k1_37_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fjji_38_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fijj_38_53_54_55_82
real, dimension (nxgrid/nprocx,3) ::  fjik_38_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_38_53_54_55_82
integer::i_38_53_54_55_82
integer::j_38_53_54_55_82
integer::k1_38_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fjji_39_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fijj_39_53_54_55_82
real, dimension (nxgrid/nprocx,3) ::  fjik_39_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_39_53_54_55_82
integer::i_39_53_54_55_82
integer::j_39_53_54_55_82
integer::k1_39_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fjji_40_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fijj_40_53_54_55_82
real, dimension (nxgrid/nprocx,3) ::  fjik_40_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_40_53_54_55_82
integer::i_40_53_54_55_82
integer::j_40_53_54_55_82
integer::k1_40_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fjji_41_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fijj_41_53_54_55_82
real, dimension (nxgrid/nprocx,3) ::  fjik_41_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_41_53_54_55_82
integer::i_41_53_54_55_82
integer::j_41_53_54_55_82
integer::k1_41_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fjji_42_53_54_55_82
real, dimension (nxgrid/nprocx,3,3)::fijj_42_53_54_55_82
real, dimension (nxgrid/nprocx,3) ::  fjik_42_53_54_55_82
real, dimension (nxgrid/nprocx) :: tmp_42_53_54_55_82
integer::i_42_53_54_55_82
integer::j_42_53_54_55_82
integer::k1_42_53_54_55_82
real, dimension(nxgrid/nprocx,3) :: tmp_46_53_54_55_82
real, dimension(nxgrid/nprocx) :: tmp_47_53_54_55_82
integer::i_47_53_54_55_82
integer::j_47_53_54_55_82
integer::kincrement_47_53_54_55_82
integer :: i_51_52_53_54_55_82
logical :: loptest_return_value_14_51_52_53_54_55_82
integer :: i_65_66_67_82
integer :: i_18_57_65_66_67_82
logical :: loptest_return_value_14_18_57_65_66_67_82
logical :: loptest_return_value_19_57_65_66_67_82
real, dimension(nxgrid/nprocx) :: del6f_upwind_22_57_65_66_67_82
integer :: msk_22_57_65_66_67_82
real, dimension(nxgrid/nprocx,3)::del6f_21_22_57_65_66_67_82
integer, dimension(nxgrid/nprocx)            :: indxs_21_22_57_65_66_67_82
integer::j_21_22_57_65_66_67_82
integer::msk_21_22_57_65_66_67_82
real, dimension (nxgrid/nprocx)::a_max_58_65_66_67_82
logical::fast_sqrt1_58_65_66_67_82
logical::precise_sqrt1_58_65_66_67_82
real, dimension (nxgrid/nprocx)::d2fdx_59_65_66_67_82
real, dimension (nxgrid/nprocx)::d2fdy_59_65_66_67_82
real, dimension (nxgrid/nprocx)::d2fdz_59_65_66_67_82
real, dimension (nxgrid/nprocx)::tmp_59_65_66_67_82
real, dimension (nxgrid/nprocx) :: tmp_60_65_66_67_82
integer::i_60_65_66_67_82
integer::j_60_65_66_67_82
real, dimension (nxgrid/nprocx) :: tmp_62_65_66_67_82
integer::i_62_65_66_67_82
integer::j_62_65_66_67_82
logical :: loptest_return_value_61_62_65_66_67_82
real, dimension (nxgrid/nprocx) :: tmp_63_65_66_67_82
integer::i_63_65_66_67_82
integer::j_63_65_66_67_82
logical :: loptest_return_value_61_63_65_66_67_82
integer :: i_51_64_65_66_67_82
logical :: loptest_return_value_14_51_64_65_66_67_82
real, dimension(nxgrid/nprocx) :: tmp_72_73_82
integer::i_72_73_82
integer::j_72_73_82
real, dimension (nxgrid/nprocx) :: tmp_69_72_73_82
integer::i_69_72_73_82
integer::j_69_72_73_82
real, dimension (nxgrid/nprocx)::d2fdx_70_72_73_82
real, dimension (nxgrid/nprocx)::d2fdy_70_72_73_82
real, dimension (nxgrid/nprocx)::d2fdz_70_72_73_82
real, dimension (nxgrid/nprocx)::tmp_70_72_73_82
real, dimension (nxgrid/nprocx)::d6fdx_71_72_73_82
real, dimension (nxgrid/nprocx)::d6fdy_71_72_73_82
real, dimension (nxgrid/nprocx)::d6fdz_71_72_73_82
logical :: ignore_dx_71_72_73_82
integer :: j_79_82
integer :: i_18_74_79_82
logical :: loptest_return_value_14_18_74_79_82
logical :: loptest_return_value_19_74_79_82
real, dimension(nxgrid/nprocx) :: del6f_upwind_22_74_79_82
integer :: msk_22_74_79_82
real, dimension(nxgrid/nprocx,3)::del6f_21_22_74_79_82
integer, dimension(nxgrid/nprocx)            :: indxs_21_22_74_79_82
integer::j_21_22_74_79_82
integer::msk_21_22_74_79_82
integer :: i_18_75_79_82
logical :: loptest_return_value_14_18_75_79_82
logical :: loptest_return_value_19_75_79_82
real, dimension(nxgrid/nprocx) :: del6f_upwind_22_75_79_82
integer :: msk_22_75_79_82
real, dimension(nxgrid/nprocx,3)::del6f_21_22_75_79_82
integer, dimension(nxgrid/nprocx)            :: indxs_21_22_75_79_82
integer::j_21_22_75_79_82
integer::msk_21_22_75_79_82
real, dimension (nxgrid/nprocx) :: tmp_76_79_82
integer::i_76_79_82
integer::j_76_79_82
logical :: loptest_return_value_61_76_79_82
integer :: i_51_78_79_82
logical :: loptest_return_value_14_51_78_79_82
real, dimension (nxgrid/nprocx,3)::tmp_80_82
real, dimension (nxgrid/nprocx,3)::tmp2_80_82
real, dimension (nxgrid/nprocx,3)::gradnu_80_82
real, dimension (nxgrid/nprocx,3)::sgradnu_80_82
real, dimension (nxgrid/nprocx,3)::gradnu_shock_80_82
real, dimension (nxgrid/nprocx)::murho1_80_82
real, dimension (nxgrid/nprocx)::zetarho1_80_82
real, dimension (nxgrid/nprocx)::mutt_80_82
real, dimension (nxgrid/nprocx)::tmp3_80_82
real, dimension (nxgrid/nprocx)::tmp4_80_82
real, dimension (nxgrid/nprocx)::pnu_shock_80_82
real, dimension (nxgrid/nprocx)::lambda_phi_80_82
real, dimension (nxgrid/nprocx)::prof_80_82
real, dimension (nxgrid/nprocx)::prof2_80_82
real, dimension (nxgrid/nprocx)::derprof_80_82
real, dimension (nxgrid/nprocx)::derprof2_80_82
real, dimension (nxgrid/nprocx)::gradnu_effective_80_82
real, dimension (nxgrid/nprocx)::fac_80_82
real, dimension (nxgrid/nprocx)::advec_hypermesh_uu_80_82
real, dimension (nxgrid/nprocx,3)::deljskl2_80_82
real, dimension (nxgrid/nprocx,3)::fvisc_nnewton2_80_82
real, dimension (nxgrid/nprocx,3,3) :: d_sld_flux_80_82
integer::i_80_82
integer::j_80_82
integer::ju_80_82
integer::ii_80_82
integer::jj_80_82
integer::kk_80_82
integer::ll_80_82
logical::ldiffus_total_80_82
logical::ldiffus_total3_80_82
real, dimension (nxgrid/nprocx,3) :: uu1_89
real, dimension (nxgrid/nprocx)::tmp_89
real, dimension (nxgrid/nprocx)::ftot_89
real, dimension (nxgrid/nprocx)::ugu_schur_x_89
real, dimension (nxgrid/nprocx)::ugu_schur_y_89
real, dimension (nxgrid/nprocx)::ugu_schur_z_89
real, dimension (nxgrid/nprocx,3,3) :: puij_schur_89
integer::i_89
integer::j_89
integer::ju_89
real::c2_86_89
real::s2_86_89
integer :: i_88_89
real, dimension (nxgrid/nprocx)::reshock_87_88_89
real, dimension (nxgrid/nprocx)::fvisc2_87_88_89
real, dimension (nxgrid/nprocx)::uus_87_88_89
real, dimension (nxgrid/nprocx)::tmp_87_88_89
real, dimension (nxgrid/nprocx)::qfvisc_87_88_89
real, dimension (nxgrid/nprocx,3)::nud2uxb_87_88_89
real, dimension (nxgrid/nprocx,3)::fluxv_87_88_89
real, dimension (nxgrid/nprocx) :: fdiff_91
real, dimension (nxgrid/nprocx) :: tmp_91
real, dimension (nxgrid/nprocx,3) :: tmpv_91
real, dimension (nxgrid/nprocx)::density_rhs_91
real, dimension (nxgrid/nprocx)::advec_hypermesh_rho_91
integer :: j_91
logical :: ldt_up_91
real::ztop_99
real::xi_99
real::profile_cor_99
real, dimension(nxgrid/nprocx) :: tmp1_99
integer :: j_99
real, dimension (nxgrid/nprocx,3)::glnthcond_96_99
real, dimension (nxgrid/nprocx,3)::glhc_96_99
real, dimension (nxgrid/nprocx) :: chix_96_99
real, dimension (nxgrid/nprocx)::thdiff_96_99
real, dimension (nxgrid/nprocx)::g2_96_99
real, dimension (nxgrid/nprocx)::del2ss1_96_99
real, dimension (nxgrid/nprocx) :: glnrhoglnt_96_99
real, dimension (nxgrid/nprocx,3) :: gradchit_prof_96_99
real, dimension (nxgrid/nprocx,3,3) :: tmp_96_99
real::s2_96_99
real::c2_96_99
real::sc_96_99
integer :: j_96_99
real, dimension(nxgrid/nprocx)::r_mn_94_96_99
real, dimension(nxgrid/nprocx)::r_mn1_94_96_99
integer :: j_94_96_99
integer :: i_95_96_99
logical :: loptest_return_value_14_95_96_99
real, dimension (nxgrid/nprocx)::tmp_98_99
real, dimension (nxgrid/nprocx)::heat_98_99
real, dimension (nxgrid/nprocx)::tt_drive_98_99
real, dimension (nxgrid/nprocx)::prof_98_99
real :: profile_buffer_98_99
real::xi_98_99
real::rgas_98_99
real, dimension (nxgrid/nprocx)::prof_97_98_99
real::zbot_97_98_99
real::ztop_97_98_99
integer :: k_101
real, dimension(nxgrid/nprocx,3) :: gg_101
real, dimension(nxgrid/nprocx) :: refac_101
real, dimension(nxgrid/nprocx)::dt1_max_loc_102
real, dimension(nxgrid/nprocx)::dt1_advec_102
real, dimension(nxgrid/nprocx)::dt1_diffus_102
real, dimension(nxgrid/nprocx)::dt1_src_102
real :: dt1_preac_102
mn_loop: do imn=1,nygrid/nprocy*nzgrid/nprocz
n=nn(imn)
m=mm(imn)
if(lpencil(i_uu)) then
p%uu=f(1+3:l2,m,n,iux:iuz)
endif
if(lpencil(i_uij)) then
k1_5_53_54_55_82=iuu-1
do i_5_53_54_55_82=1,3
do j_5_53_54_55_82=1,3
call der(f,k1_5_53_54_55_82+i_5_53_54_55_82,tmp_5_53_54_55_82,j_5_53_54_55_82)
p%uij(:,i_5_53_54_55_82,j_5_53_54_55_82)=tmp_5_53_54_55_82
enddo
enddo
endif
if(lpencil(i_divu)) then
p%divu=p%uij(:,1,1)+p%uij(:,2,2)+p%uij(:,3,3)
endif
if(lpencil(i_sij)) then
do j_8_53_54_55_82=1,3
p%sij(:,j_8_53_54_55_82,j_8_53_54_55_82)=p%uij(:,j_8_53_54_55_82,j_8_53_54_55_82)-(1./3.)*p%divu
do i_8_53_54_55_82=j_8_53_54_55_82+1,3
p%sij(:,i_8_53_54_55_82,j_8_53_54_55_82)=.5*(p%uij(:,i_8_53_54_55_82,j_8_53_54_55_82)+p%uij(:,j_8_53_54_55_82,i_8_53_54_55_82))
p%sij(:,j_8_53_54_55_82,i_8_53_54_55_82)=p%sij(:,i_8_53_54_55_82,j_8_53_54_55_82)
enddo
enddo
endif
if(lpencil(i_sij2)) then
p%sij2 = p%sij(:,1,1)**2
do i_9_53_54_55_82 = 2, 3
p%sij2 = p%sij2 + p%sij(:,i_9_53_54_55_82,i_9_53_54_55_82)**2
do j_9_53_54_55_82 = 1, i_9_53_54_55_82-1
p%sij2 = p%sij2 + 2 * p%sij(:,i_9_53_54_55_82,j_9_53_54_55_82)**2
enddo
enddo
endif
if(lpencil(i_uij5)) then
k1_10_53_54_55_82=iuu-1
do i_10_53_54_55_82=1,3
do j_10_53_54_55_82=1,3
enddo
enddo
endif
if(lpencil(i_oo)) then
a1_12_53_54_55_82 = 1
a2_12_53_54_55_82 = nxgrid/nprocx
if(a2_12_53_54_55_82 == nxgrid/nprocx+2*3) then
a1_12_53_54_55_82 = 1+3
a2_12_53_54_55_82 = l2
endif
endif
if(lpencil(i_ou)) then
do i_15_53_54_55_82=2,3
enddo
endif
if(lpencil(i_ugu)) then
do j_25_53_54_55_82=1,3
grad_f_tmp_25_53_54_55_82 = p%uij(:,j_25_53_54_55_82,:)
tmp_25_53_54_55_82=p%uu(:,1)*grad_f_tmp_25_53_54_55_82(:,1)
do i_18_23_25_53_54_55_82=2,3
tmp_25_53_54_55_82=tmp_25_53_54_55_82+p%uu(:,i_18_23_25_53_54_55_82)*grad_f_tmp_25_53_54_55_82(:,i_18_23_25_53_54_55_82)
enddo
p%ugu(:,j_25_53_54_55_82)=tmp_25_53_54_55_82
enddo
endif
if(lpencil(i_ogu)) then
do j_27_53_54_55_82=1,3
do i_18_23_27_53_54_55_82=2,3
enddo
enddo
endif
if(lpencil(i_del4u)) then
k1_29_53_54_55_82=iuu-1
do i_29_53_54_55_82=1,3
enddo
endif
if(lpencil(i_del6u)) then
k1_33_53_54_55_82=iuu-1
do i_33_53_54_55_82=1,3
enddo
endif
if(lpencil(i_del6u_strict)) then
k1_34_53_54_55_82=iuu-1
do i_34_53_54_55_82=1,3
do i_31_34_53_54_55_82=1,3
do j_31_34_53_54_55_82=1,3
enddo
enddo
enddo
endif
if(lpencil(i_del4graddivu)) then
do i_35_53_54_55_82=1,3
ki_35_53_54_55_82 = iuu + (i_35_53_54_55_82-1)
do j_35_53_54_55_82=1,3
enddo
do j_35_53_54_55_82=1,3
if(j_35_53_54_55_82/=i_35_53_54_55_82) then
if((i_35_53_54_55_82==1).and.(j_35_53_54_55_82==2)) then
k_35_53_54_55_82=3
endif
if((i_35_53_54_55_82==1).and.(j_35_53_54_55_82==3)) then
k_35_53_54_55_82=2
endif
if((i_35_53_54_55_82==2).and.(j_35_53_54_55_82==1)) then
k_35_53_54_55_82=3
endif
if((i_35_53_54_55_82==2).and.(j_35_53_54_55_82==3)) then
k_35_53_54_55_82=1
endif
if((i_35_53_54_55_82==3).and.(j_35_53_54_55_82==1)) then
k_35_53_54_55_82=2
endif
if((i_35_53_54_55_82==3).and.(j_35_53_54_55_82==2)) then
k_35_53_54_55_82=1
endif
kj_35_53_54_55_82 = iuu+(j_35_53_54_55_82-1)
endif
enddo
enddo
endif
if(lpencil(i_der6u_res)) then
do j_53_54_55_82=1,3
ju_53_54_55_82=j_53_54_55_82+iuu-1
do i_53_54_55_82=1,3
enddo
enddo
endif
if(lpencil(i_del2u).and.lpencil(i_graddivu).and.lpencil(i_curlo)) then
k1_36_53_54_55_82=iuu-1
do i_36_53_54_55_82=1,3
do j_36_53_54_55_82=1,3
call der2 (f,k1_36_53_54_55_82+i_36_53_54_55_82,tmp_36_53_54_55_82,  j_36_53_54_55_82)
fijj_36_53_54_55_82(:,i_36_53_54_55_82,j_36_53_54_55_82)=tmp_36_53_54_55_82
call derij(f,k1_36_53_54_55_82+j_36_53_54_55_82,tmp_36_53_54_55_82,j_36_53_54_55_82,i_36_53_54_55_82)
fjji_36_53_54_55_82(:,i_36_53_54_55_82,j_36_53_54_55_82)=tmp_36_53_54_55_82
enddo
enddo
do i_36_53_54_55_82=1,3
p%del2u(:,i_36_53_54_55_82)=fijj_36_53_54_55_82(:,i_36_53_54_55_82,1)+fijj_36_53_54_55_82(:,i_36_53_54_55_82,2)+fijj_36_53_54_55_82(:,i_36_53_54_55_82,3)
enddo
do i_36_53_54_55_82=1,3
p%graddivu(:,i_36_53_54_55_82)=fjji_36_53_54_55_82(:,i_36_53_54_55_82,1)+fjji_36_53_54_55_82(:,i_36_53_54_55_82,2)+fjji_36_53_54_55_82(:,i_36_53_54_55_82,3)
enddo
else if(lpencil(i_del2u).and.lpencil(i_graddivu)) then
k1_37_53_54_55_82=iuu-1
do i_37_53_54_55_82=1,3
do j_37_53_54_55_82=1,3
call der2 (f,k1_37_53_54_55_82+i_37_53_54_55_82,tmp_37_53_54_55_82,  j_37_53_54_55_82)
fijj_37_53_54_55_82(:,i_37_53_54_55_82,j_37_53_54_55_82)=tmp_37_53_54_55_82
call derij(f,k1_37_53_54_55_82+j_37_53_54_55_82,tmp_37_53_54_55_82,j_37_53_54_55_82,i_37_53_54_55_82)
fjji_37_53_54_55_82(:,i_37_53_54_55_82,j_37_53_54_55_82)=tmp_37_53_54_55_82
enddo
enddo
do i_37_53_54_55_82=1,3
p%del2u(:,i_37_53_54_55_82)=fijj_37_53_54_55_82(:,i_37_53_54_55_82,1)+fijj_37_53_54_55_82(:,i_37_53_54_55_82,2)+fijj_37_53_54_55_82(:,i_37_53_54_55_82,3)
enddo
do i_37_53_54_55_82=1,3
p%graddivu(:,i_37_53_54_55_82)=fjji_37_53_54_55_82(:,i_37_53_54_55_82,1)+fjji_37_53_54_55_82(:,i_37_53_54_55_82,2)+fjji_37_53_54_55_82(:,i_37_53_54_55_82,3)
enddo
else if(lpencil(i_del2u).and.lpencil(i_curlo)) then
k1_38_53_54_55_82=iuu-1
do i_38_53_54_55_82=1,3
do j_38_53_54_55_82=1,3
call der2 (f,k1_38_53_54_55_82+i_38_53_54_55_82,tmp_38_53_54_55_82,  j_38_53_54_55_82)
fijj_38_53_54_55_82(:,i_38_53_54_55_82,j_38_53_54_55_82)=tmp_38_53_54_55_82
call derij(f,k1_38_53_54_55_82+j_38_53_54_55_82,tmp_38_53_54_55_82,j_38_53_54_55_82,i_38_53_54_55_82)
enddo
enddo
do i_38_53_54_55_82=1,3
p%del2u(:,i_38_53_54_55_82)=fijj_38_53_54_55_82(:,i_38_53_54_55_82,1)+fijj_38_53_54_55_82(:,i_38_53_54_55_82,2)+fijj_38_53_54_55_82(:,i_38_53_54_55_82,3)
enddo
do i_38_53_54_55_82=1,3
enddo
else if(lpencil(i_graddivu).and.lpencil(i_curlo)) then
k1_39_53_54_55_82=iuu-1
do i_39_53_54_55_82=1,3
do j_39_53_54_55_82=1,3
call der2 (f,k1_39_53_54_55_82+i_39_53_54_55_82,tmp_39_53_54_55_82,  j_39_53_54_55_82)
call derij(f,k1_39_53_54_55_82+j_39_53_54_55_82,tmp_39_53_54_55_82,j_39_53_54_55_82,i_39_53_54_55_82)
fjji_39_53_54_55_82(:,i_39_53_54_55_82,j_39_53_54_55_82)=tmp_39_53_54_55_82
enddo
enddo
do i_39_53_54_55_82=1,3
enddo
do i_39_53_54_55_82=1,3
p%graddivu(:,i_39_53_54_55_82)=fjji_39_53_54_55_82(:,i_39_53_54_55_82,1)+fjji_39_53_54_55_82(:,i_39_53_54_55_82,2)+fjji_39_53_54_55_82(:,i_39_53_54_55_82,3)
enddo
else if(lpencil(i_del2u)) then
k1_40_53_54_55_82=iuu-1
do i_40_53_54_55_82=1,3
do j_40_53_54_55_82=1,3
call der2 (f,k1_40_53_54_55_82+i_40_53_54_55_82,tmp_40_53_54_55_82,  j_40_53_54_55_82)
fijj_40_53_54_55_82(:,i_40_53_54_55_82,j_40_53_54_55_82)=tmp_40_53_54_55_82
call derij(f,k1_40_53_54_55_82+j_40_53_54_55_82,tmp_40_53_54_55_82,j_40_53_54_55_82,i_40_53_54_55_82)
enddo
enddo
do i_40_53_54_55_82=1,3
p%del2u(:,i_40_53_54_55_82)=fijj_40_53_54_55_82(:,i_40_53_54_55_82,1)+fijj_40_53_54_55_82(:,i_40_53_54_55_82,2)+fijj_40_53_54_55_82(:,i_40_53_54_55_82,3)
enddo
do i_40_53_54_55_82=1,3
enddo
else if(lpencil(i_graddivu)) then
k1_41_53_54_55_82=iuu-1
do i_41_53_54_55_82=1,3
do j_41_53_54_55_82=1,3
call der2 (f,k1_41_53_54_55_82+i_41_53_54_55_82,tmp_41_53_54_55_82,  j_41_53_54_55_82)
call derij(f,k1_41_53_54_55_82+j_41_53_54_55_82,tmp_41_53_54_55_82,j_41_53_54_55_82,i_41_53_54_55_82)
fjji_41_53_54_55_82(:,i_41_53_54_55_82,j_41_53_54_55_82)=tmp_41_53_54_55_82
enddo
enddo
do i_41_53_54_55_82=1,3
enddo
do i_41_53_54_55_82=1,3
p%graddivu(:,i_41_53_54_55_82)=fjji_41_53_54_55_82(:,i_41_53_54_55_82,1)+fjji_41_53_54_55_82(:,i_41_53_54_55_82,2)+fjji_41_53_54_55_82(:,i_41_53_54_55_82,3)
enddo
else if(lpencil(i_curlo)) then
k1_42_53_54_55_82=iuu-1
do i_42_53_54_55_82=1,3
do j_42_53_54_55_82=1,3
enddo
enddo
do i_42_53_54_55_82=1,3
enddo
do i_42_53_54_55_82=1,3
enddo
endif
if(lpencil(i_uijk)) then
do kincrement_47_53_54_55_82=0,2
do i_47_53_54_55_82=1,3
do j_47_53_54_55_82=i_47_53_54_55_82,3
enddo
enddo
enddo
endif
if(lpencil(i_grad5divu)) then
do i_53_54_55_82=1,3
do j_53_54_55_82=1,3
ju_53_54_55_82=iuu+j_53_54_55_82-1
enddo
enddo
endif
if(lpencil(i_uu_advec)) then
do j_53_54_55_82=1,3
do i_51_52_53_54_55_82=2,3
enddo
enddo
endif
p%lnrho=f(1+3:l2,m,n,ilnrho)
if(lpencil(i_rho1)) then
p%rho1=exp(-f(1+3:l2,m,n,ilnrho))
endif
if(lpencil(i_rho)) then
p%rho=1.0/p%rho1
endif
if(lpencil(i_glnrho).or.lpencil(i_grho)) then
call der(f,ilnrho,p%glnrho(:,1),1)
call der(f,ilnrho,p%glnrho(:,2),2)
call der(f,ilnrho,p%glnrho(:,3),3)
if(lpencil(i_grho)) then
do i_65_66_67_82=1,3
enddo
endif
endif
if(lpencil(i_uglnrho)) then
p%uglnrho=p%uu(:,1)*p%glnrho(:,1)
do i_18_57_65_66_67_82=2,3
p%uglnrho=p%uglnrho+p%uu(:,i_18_57_65_66_67_82)*p%glnrho(:,i_18_57_65_66_67_82)
enddo
msk_22_57_65_66_67_82=0
msk_21_22_57_65_66_67_82=0
msk_21_22_57_65_66_67_82=msk_22_57_65_66_67_82
do j_21_22_57_65_66_67_82=1,3
if(j_21_22_57_65_66_67_82==msk_21_22_57_65_66_67_82) then
del6f_21_22_57_65_66_67_82(:,j_21_22_57_65_66_67_82) = 0.
else
if(lequidist(j_21_22_57_65_66_67_82) .or. .false.) then
call der6(f,ilnrho,del6f_21_22_57_65_66_67_82(:,j_21_22_57_65_66_67_82),j_21_22_57_65_66_67_82,upwind=.true.)
else
where(p%uu(:,j_21_22_57_65_66_67_82)>=0)
indxs_21_22_57_65_66_67_82 = 7
elsewhere
indxs_21_22_57_65_66_67_82 = 8
endwhere
endif
del6f_21_22_57_65_66_67_82(:,j_21_22_57_65_66_67_82) = abs(p%uu(:,j_21_22_57_65_66_67_82))*del6f_21_22_57_65_66_67_82(:,j_21_22_57_65_66_67_82)
endif
enddo
del6f_upwind_22_57_65_66_67_82 = sum(del6f_21_22_57_65_66_67_82,2)
p%uglnrho = p%uglnrho-del6f_upwind_22_57_65_66_67_82
endif
if(lpencil(i_del2lnrho)) then
call der2(f,ilnrho,d2fdx_59_65_66_67_82,1)
call der2(f,ilnrho,d2fdy_59_65_66_67_82,2)
call der2(f,ilnrho,d2fdz_59_65_66_67_82,3)
p%del2lnrho=d2fdx_59_65_66_67_82+d2fdy_59_65_66_67_82+d2fdz_59_65_66_67_82
endif
if(lpencil(i_hlnrho)) then
do j_60_65_66_67_82=1,3
do i_60_65_66_67_82=j_60_65_66_67_82+1,3
enddo
enddo
endif
if(lpencil(i_sglnrho)) then
do i_62_65_66_67_82=1,3
j_62_65_66_67_82=1
tmp_62_65_66_67_82=p%sij(:,i_62_65_66_67_82,j_62_65_66_67_82)*p%glnrho(:,j_62_65_66_67_82)
do j_62_65_66_67_82=2,3
tmp_62_65_66_67_82=tmp_62_65_66_67_82+p%sij(:,i_62_65_66_67_82,j_62_65_66_67_82)*p%glnrho(:,j_62_65_66_67_82)
enddo
p%sglnrho(:,i_62_65_66_67_82)=tmp_62_65_66_67_82
enddo
endif
if(lpencil(i_uij5glnrho)) then
do i_63_65_66_67_82=1,3
j_63_65_66_67_82=1
do j_63_65_66_67_82=2,3
enddo
enddo
endif
if(lpencil(i_uuadvec_glnrho)) then
do i_51_64_65_66_67_82=2,3
enddo
endif
if(lpencil(i_cv)) then
p%cv=1/cv1
endif
if(lpencil(i_ss)) then
p%ss=f(1+3:l2,m,n,iss)
endif
if(lpencil(i_gss)) then
call der(f,iss,p%gss(:,1),1)
call der(f,iss,p%gss(:,2),2)
call der(f,iss,p%gss(:,3),3)
endif
if(lpencil(i_hss)) then
do j_69_72_73_82=1,3
do i_69_72_73_82=j_69_72_73_82+1,3
enddo
enddo
endif
if(lpencil(i_del2ss)) then
call der2(f,iss,d2fdx_70_72_73_82,1)
call der2(f,iss,d2fdy_70_72_73_82,2)
call der2(f,iss,d2fdz_70_72_73_82,3)
p%del2ss=d2fdx_70_72_73_82+d2fdy_70_72_73_82+d2fdz_70_72_73_82
endif
if(lpencil(i_cs2)) then
p%cs2=cs20*exp(cv1*p%ss+gamma_m1*(p%lnrho-lnrho0))
endif
if(lpencil(i_lntt)) then
p%lntt=lntt0+cv1*p%ss+gamma_m1*(p%lnrho-lnrho0)
endif
if(lpencil(i_tt)) then
p%tt=exp(p%lntt)
endif
if(lpencil(i_tt1)) then
p%tt1=exp(-p%lntt)
endif
if(lpencil(i_glntt)) then
p%glntt=gamma_m1*p%glnrho+cv1*p%gss
endif
if(lpencil(i_gtt)) then
do j_72_73_82=1,3
enddo
endif
if(lpencil(i_del2lntt)) then
p%del2lntt=gamma_m1*p%del2lnrho+cv1*p%del2ss
endif
if(lpencil(i_ugss)) then
p%ugss=p%uu(:,1)*p%gss(:,1)
do i_18_74_79_82=2,3
p%ugss=p%ugss+p%uu(:,i_18_74_79_82)*p%gss(:,i_18_74_79_82)
enddo
msk_22_74_79_82=0
msk_21_22_74_79_82=0
msk_21_22_74_79_82=msk_22_74_79_82
do j_21_22_74_79_82=1,3
if(j_21_22_74_79_82==msk_21_22_74_79_82) then
del6f_21_22_74_79_82(:,j_21_22_74_79_82) = 0.
else
if(lequidist(j_21_22_74_79_82) .or. .false.) then
call der6(f,iss,del6f_21_22_74_79_82(:,j_21_22_74_79_82),j_21_22_74_79_82,upwind=.true.)
else
where(p%uu(:,j_21_22_74_79_82)>=0)
indxs_21_22_74_79_82 = 7
elsewhere
indxs_21_22_74_79_82 = 8
endwhere
endif
del6f_21_22_74_79_82(:,j_21_22_74_79_82) = abs(p%uu(:,j_21_22_74_79_82))*del6f_21_22_74_79_82(:,j_21_22_74_79_82)
endif
enddo
del6f_upwind_22_74_79_82 = sum(del6f_21_22_74_79_82,2)
p%ugss = p%ugss-del6f_upwind_22_74_79_82
endif
if(lpencil(i_uglntt)) then
do i_18_75_79_82=2,3
enddo
msk_22_75_79_82=0
msk_21_22_75_79_82=0
msk_21_22_75_79_82=msk_22_75_79_82
do j_21_22_75_79_82=1,3
if(j_21_22_75_79_82==msk_21_22_75_79_82) then
else
if(lequidist(j_21_22_75_79_82) .or. .false.) then
else
where(p%uu(:,j_21_22_75_79_82)>=0)
indxs_21_22_75_79_82 = 7
elsewhere
indxs_21_22_75_79_82 = 8
endwhere
endif
endif
enddo
endif
if(lpencil(i_sglntt)) then
do i_76_79_82=1,3
j_76_79_82=1
do j_76_79_82=2,3
enddo
enddo
endif
if(lpencil(i_fpres)) then
do j_79_82=1,3
p%fpres(:,j_79_82)=-p%cs2*(p%glnrho(:,j_79_82) + p%glntt(:,j_79_82))*gamma1
enddo
endif
if(lpencil(i_uuadvec_gss)) then
do i_51_78_79_82=2,3
enddo
endif
p%fvisc=0.0
if(lpencil(i_visc_heat)) then
p%visc_heat=0.0
endif
fac_80_82=4.00000019e-03
do j_80_82=1,3
p%fvisc(:,j_80_82) = p%fvisc(:,j_80_82) + fac_80_82*(p%del2u(:,j_80_82) + 2.*p%sglnrho(:,j_80_82) + 1./3.*p%graddivu(:,j_80_82))
enddo
if(lpencil(i_visc_heat)) then
p%visc_heat=p%visc_heat+2*4.00000019e-03*p%sij2
endif
if(lpencil(i_gg)) then
p%gg(:,3) = gravz_zpencil(n)
endif
df(1+3:l2,m,n,iux:iuz)=df(1+3:l2,m,n,iux:iuz)-p%ugu
c2_86_89=2*0.100000001
df(1+3:l2,m,n,iux  )=df(1+3:l2,m,n,iux  )+c2_86_89*p%uu(:,2)
df(1+3:l2,m,n,iux+1)=df(1+3:l2,m,n,iux+1)-c2_86_89*p%uu(:,1)
df(1+3:l2,m,n,iux:iuz) = df(1+3:l2,m,n,iux:iuz) + p%fvisc
density_rhs_91= - p%divu
density_rhs_91 = density_rhs_91 - p%uglnrho
df(1+3:l2,m,n,ilnrho) = df(1+3:l2,m,n,ilnrho) + density_rhs_91
fdiff_91=0.0
tmp_91 = fdiff_91
forall(j_91 = iux:iuz) df(1+3:l2,m,n,j_91) = df(1+3:l2,m,n,j_91) - p%uu(:,j_91-iuu+1) * tmp_91
df(1+3:l2,m,n,iss) = df(1+3:l2,m,n,iss) - p%cv*tmp_91
df(1+3:l2,m,n,ilnrho) = df(1+3:l2,m,n,ilnrho) + fdiff_91
df(1+3:l2,m,n,iux:iuz) = df(1+3:l2,m,n,iux:iuz) + p%fpres
if(any(beta_glnrho_scaled/=0.0)) then
do j_99=1,3
df(1+3:l2,m,n,(iux-1)+j_99) = df(1+3:l2,m,n,(iux-1)+j_99) - p%cs2*beta_glnrho_scaled(j_99)
enddo
endif
df(1+3:l2,m,n,iss) = df(1+3:l2,m,n,iss) - p%ugss
do j_99=1,3
if(grads0_imposed(j_99)/=0.) then
df(1+3:l2,m,n,iss)=df(1+3:l2,m,n,iss)-grads0_imposed(j_99)*p%uu(:,j_99)
endif
enddo
df(1+3:l2,m,n,iss) = df(1+3:l2,m,n,iss) + p%tt1*p%visc_heat
if(hcond0/=0..or..false.) then
hcond=hcond_prof(n-3)
glhc_96_99(:,3)=dlnhcond_prof(n-3)
glhc_96_99(:,1:2)=0.
glnthcond_96_99 = p%glntt + glhc_96_99
g2_96_99=p%glntt(:,1)*glnthcond_96_99(:,1)
do i_95_96_99=2,3
g2_96_99=g2_96_99+p%glntt(:,i_95_96_99)*glnthcond_96_99(:,i_95_96_99)
enddo
thdiff_96_99 = p%rho1*hcond * (p%del2lntt + g2_96_99)
endif
if(hcond0/=0..or..false..or.0.00000000/=0.) then
df(1+3:l2,m,n,iss) = df(1+3:l2,m,n,iss) + thdiff_96_99
endif
heat_98_99=0.0
ztop_97_98_99=-0.680000007+2.00000000
prof_97_98_99 = spread(exp(-0.5*((ztop_97_98_99-z(n))/0.200000003)**2), 1, l2-4+1)
heat_98_99 = heat_98_99 - 15.0000000*prof_97_98_99*(p%cs2-cs2cool)/cs2cool
df(1+3:l2,m,n,iss) = df(1+3:l2,m,n,iss) + p%tt1*p%rho1*heat_98_99
gg_101(:,3)=p%gg(:,3)
df(1+3:l2,m,n,iuz)=df(1+3:l2,m,n,iuz)+gg_101(:,3)
enddo mn_loop
endsubroutine rhs_cpu_test
!***********************************************************************
endmodule Equ
