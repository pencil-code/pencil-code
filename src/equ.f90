! $Id$
!
!  subroutines in the chosen set of physics modules.
!
module Equ
!
  use Cdata
  use Messages
!
  implicit none
!
  private
!
  public :: pde, debug_imn_arrays, initialize_pencils
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
!
      use Boundcond
      use BorderProfiles, only: calc_pencils_borderprofiles
      use Chiral
      use Chemistry
      use Cosmicray
      use CosmicrayFlux
      use Density
      use Diagnostics
      use Dustvelocity
      use Dustdensity
      use Entropy
      use EquationOfState
      use Forcing, only: calc_pencils_forcing, calc_lforcing_cont_pars, &
                         forcing_continuous
! To check ghost cell consistency, please uncomment the following line:
!     use Ghost_check, only: check_ghosts_consistency
      use GhostFold, only: fold_df
      use Gravity
      use Grid, only: calc_pencils_grid
      use Hydro
      use Interstellar, only: interstellar_before_boundary
      use Lorenz_gauge
      use Magnetic
      use Hypervisc_strict, only: hyperviscosity_strict
      use Hyperresi_strict, only: hyperresistivity_strict
      use Mpicomm
      use NeutralDensity
      use NeutralVelocity
      use NSCBC
      use Particles_main
      use Poisson
      use Pscalar
      use Polymer
      use Radiation
      use Selfgravity
      use Shear
      use Shock, only: calc_pencils_shock, calc_shock_profile, &
                       calc_shock_profile_simple
      use Solid_Cells, only: update_solid_cells, freeze_solid_cells, &
          dsolid_dt,dsolid_dt_integrate
      use Special
      use Sub
      use Testfield
      use Testflow
      use Testscalar
      use Viscosity, only: calc_viscosity, calc_pencils_viscosity, dynamical_viscosity
!
      logical :: early_finalize
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: df_iuu_pencil
      type (pencil_case) :: p
      real, dimension (nx) :: maxadvec,maxdiffus,maxdiffus2,maxdiffus3
      real, dimension (nx) :: advec2,advec2_hypermesh
      real, dimension (nx) :: pfreeze,pfreeze_int,pfreeze_ext
      real, dimension(1) :: mass_per_proc
      integer :: iv
      integer :: ivar1,ivar2
      real :: umax = 0.
      intent(inout)  :: f       ! inout due to  lshift_datacube_x,
                                ! density floor, or velocity ceiling
      intent(out)    :: df, p
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
      if (ldiagnos   ) tdiagnos   =t ! (diagnostics are for THIS time)
      if (l1davgfirst) t1ddiagnos =t ! (1-D averages are for THIS time)
      if (l2davgfirst) t2davgfirst=t ! (2-D averages are for THIS time)
!
!  Grid spacing. For non equidistant grid or non-cartesian coordinates
!  the grid spacing is calculated in the (m,n) loop below.
!
      if (lcartesian_coords .and. all(lequidist)) then
        if (old_cdtv) then
          dxyz_2 = max(dx_1(l1:l2)**2,dy_1(m1)**2,dz_1(n1)**2)
        else
          dline_1(:,1)=dx_1(l1:l2)
          dline_1(:,2)=dy_1(m1)
          dline_1(:,3)=dz_1(n1)
          dxyz_2 = dline_1(:,1)**2+dline_1(:,2)**2+dline_1(:,3)**2
          dxyz_4 = dline_1(:,1)**4+dline_1(:,2)**4+dline_1(:,3)**4
          dxyz_6 = dline_1(:,1)**6+dline_1(:,2)**6+dline_1(:,3)**6
        endif
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
                     lparticles_prepencil_calc.or.lsolid_cells.or. &
                     lchemistry.or.lweno_transport
!
!  Write crash snapshots to the hard disc if the time-step is very low.
!  The user must have set crash_file_dtmin_factor>0.0 in &run_pars for
!  this to be done.
!
      if (crash_file_dtmin_factor > 0.0) call output_crash_files(f)
!
!  For debugging purposes impose minimum or maximum value on certain variables.
!
      call impose_density_floor(f)
      call impose_velocity_ceiling(f)
!
!  Apply global boundary conditions to particle positions and communicate
!  migrating particles between the processors.
!
      if (lparticles) call particles_boundconds(f)
!
!  Calculate the potential of the self gravity. Must be done before
!  communication in order to be able to take the gradient of the potential
!  later.
!
      call calc_selfpotential(f)
!
!  Remove mean emf in the radial direction if desired.
!  Useful as a simple way to remove the large scale
!  contribution from uphi x Bz from global disk simulations.
!  (This is only done if lremove_mean_emf=T,
!  to be set in magnetic_run_pars).
!
      if (lmagnetic) call remove_mean_emf(f,df)
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
      if (ldensity)      call density_before_boundary(f)
      if (lhydro)        call hydro_before_boundary(f)
      if (lshear)        call shear_before_boundary(f)
      if (lchiral)       call chiral_before_boundary(f)
      if (lspecial)      call special_before_boundary(f)
      if (ltestflow)     call testflow_before_boundary(f)
      if (lparticles)    call particles_before_boundary(f)
!
!  Fetch fp to the special module
!
      if (lparticles.and.lspecial) call particles_special
!
!  Initiate shock profile calculation and use asynchronous to handle
!  communication along processor/periodic boundaries.
!
      if (lshock) call calc_shock_profile(f)
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
      if (ldebug) print*,'pde: bef. initiate_isendrcv_bdry'
      call initiate_isendrcv_bdry(f)
      if (early_finalize) then
        call finalize_isendrcv_bdry(f)
        call boundconds_y(f)
        call boundconds_z(f)
      endif
!
! update solid cell "ghost points". This must be done in order to get the
! correct boundary layer close to the solid geometry, i.e. no-slip conditions.
!
      call update_solid_cells(f)
!
!  Give the particle modules a chance to do something special with a fully
!  communicated f array, for instance: the particles_spin module needs to
!  maintain the full vorticity field, including ghost zones, to be able to do
!  interpolation on the vorticity to subgrid particles positions.
!
      if (early_finalize.and.lparticles_prepencil_calc) then
        call particles_doprepencil_calc(f,ivar1,ivar2)
        if (ivar1>=0 .and. ivar2>=0) then
          call boundconds_x(f,ivar1,ivar2)
          call initiate_isendrcv_bdry(f,ivar1,ivar2)
          call finalize_isendrcv_bdry(f,ivar1,ivar2)
          call boundconds_y(f,ivar1,ivar2)
          call boundconds_z(f,ivar1,ivar2)
        endif
      endif
!
!  For sixth order momentum-conserving, symmetric hyperviscosity with positive
!  definite heating rate we need to precalculate the viscosity term. The
!  restivitity term for sixth order hyperresistivity with positive definite
!  heating rate must also be precalculated.
!
      if (lhyperviscosity_strict)   call hyperviscosity_strict(f)
      if (lhyperresistivity_strict) call hyperresistivity_strict(f)
!
!  Dynamical diffusion coefficients
!
      if (lfirst .and. ldynamical_diffusion) then
        call mpiallreduce_max(sqrt(maxval(f(l1:l2,m1:m2,n1:n2,iux)**2 &
                                        + f(l1:l2,m1:m2,n1:n2,iuy)**2 &
                                        + f(l1:l2,m1:m2,n1:n2,iuz)**2)), umax)
        if (lviscosity) call dynamical_viscosity(umax)
        if (ldensity) call dynamical_diffusion(umax)
        if (lmagnetic) call dynamical_resistivity(umax)
      endif
!
!  For calculating the pressure gradient directly from the pressure (which is
!  derived from the basic thermodynamical variables), we need to fill in the
!  pressure in the f array.
!
      call fill_farray_pressure(f)
!
!  Set inverse timestep to zero before entering loop over m and n.
!
      if (lfirst.and.ldt) then
        if (dtmax/=0.0) then
          dt1_max=1/dtmax
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
!  for example. The used to be or are still called calc_lhydro_pars etc,
!  and will soon be renamed to hydro_after_boundary.
!
      call timing('pde','before calc_lhydro_pars')
!DM I suggest the following lhydro_pars, lmagnetic_pars be renamed to
! hydro_after_boundary etc.
!AB: yes, we should rename these step by step
      if (lhydro)                 call calc_lhydro_pars(f)
      if (lmagnetic)              call calc_lmagnetic_pars(f)
      if (lentropy)               call calc_lentropy_pars(f)
      if (lforcing_cont)          call calc_lforcing_cont_pars(f)
      if (lpolymer)               call calc_polymer_after_boundary(f)
      if (ltestscalar)            call testscalar_after_boundary(f)
      if (ltestfield)             call testfield_after_boundary(f,p)
      if (ldensity)               call calc_ldensity_pars(f)
      if (ltestflow)              call calc_ltestflow_nonlin_terms(f,df)
      if (lspecial)               call calc_lspecial_pars(f)
!AB: could be renamed to special_after_boundary etc
      !if (lspecial)               call special_after_boundary(f)
!
!  Calculate quantities for a chemical mixture
!
      if (lchemistry .and. ldustdensity) then
        call chemspec_normalization(f)
!        call dustspec_normalization(f)
      endif
      if (lchemistry .and. ldensity) call calc_for_chem_mixture(f)
      call timing('pde','after calc_for_chem_mixture')
!
!------------------------------------------------------------------------------
!  Do loop over m and n.
!
      mn_loop: do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
        lfirstpoint=(imn==1)      ! true for very first m-n loop
        llastpoint=(imn==(ny*nz)) ! true for very last m-n loop
!
!  Store the velocity part of df array in a temporary array
!  while solving the anelastic case.
!
        call timing('pde','before lanelastic',mnloop=.true.)
        if (lanelastic) then
          df_iuu_pencil(1:nx,1:3) = df(l1:l2,m,n,iuu:iuu+2)
          df(l1:l2,m,n,iuu:iuu+2)=0.0
        endif
!
!        if (loptimise_ders) der_call_count=0 !DERCOUNT
!
!  Make sure all ghost points are set.
!
        if (.not.early_finalize.and.necessary(imn)) then
          call finalize_isendrcv_bdry(f)
          call boundconds_y(f)
          call boundconds_z(f)
        endif
        call timing('pde','finished boundconds_z',mnloop=.true.)
!
!  For each pencil, accumulate through the different modules
!  advec_XX and diffus_XX, which are essentially the inverse
!  advective and diffusive timestep for that module.
!  (note: advec_cs2 and advec_va2 are inverse _squared_ timesteps)
!
        if (lfirst.and.ldt.and.(.not.ldt_paronly)) then
          if (lhydro.or.lhydro_kinematic) then
            advec_uu=0.0; advec_hypermesh_uu=0.0
          endif
          if (ldensity) then
            diffus_diffrho=0.0; diffus_diffrho3=0.0; advec_hypermesh_rho=0.0
            if (leos) advec_cs2=0.0
          endif
          if (lenergy) then
            diffus_chi=0.0; diffus_chi3=0.0; advec_hypermesh_ss=0.0
          endif
          if (lmagnetic) then
            advec_va2=0.0; advec_hall=0.0; advec_hypermesh_aa=0.0
            diffus_eta=0.0; diffus_eta2=0.0; diffus_eta3=0.0
          endif
          if (ltestfield) then
            diffus_eta=0.0; diffus_eta3=0.0
          endif
          if (ltestscalar) then
            diffus_eta=0.0
          endif
          if (ldustvelocity) then
            advec_uud=0.0; diffus_nud=0.0; diffus_nud3=0.0
          endif
          if (lpscalar) then
            diffus_pscalar=0.0; diffus_pscalar3=0.0
          endif
          if (ldustdensity) then
            diffus_diffnd=0.0; diffus_diffnd3=0.0
          endif
          if (lviscosity) then
            diffus_nu=0.0; diffus_nu2=0.0; diffus_nu3=0.0
          endif
          if (lradiation) then
            advec_crad2=0.0
          endif
          if (lshear) then
            advec_shear=0.0
          endif
          if (lchemistry) then
            diffus_chem=0.0
          endif
          if (lchiral) then
            diffus_chiral=0.0
          endif
          if (lcosmicray) then
            diffus_cr=0.0
          endif
          if (lneutraldensity) then
            diffus_diffrhon=0.0; diffus_diffrhon3=0.0
          endif
          if (lneutralvelocity) then
            advec_uun=0.0; advec_csn2=0.0; diffus_nun=0.0; diffus_nun3=0.0
          endif
          if (lpolymer) then
            advec_poly=0.0; diffus_poly=0.0
          endif
        endif
!
!  Grid spacing. In case of equidistant grid and cartesian coordinates
!  this is calculated before the (m,n) loop.
!
        if (lspherical_coords.or.lcylindrical_coords.or. &
            .not.all(lequidist)) then
          if (old_cdtv) then
!
!  The following is only kept for backwards compatibility. Will be deleted in
!  the future.
!
            dxyz_2 = max(dx_1(l1:l2)**2,dy_1(m)**2,dz_1(n)**2)
          else
            if (lspherical_coords) then
              dline_1(:,1)=dx_1(l1:l2)
              dline_1(:,2)=r1_mn*dy_1(m)
              dline_1(:,3)=r1_mn*sin1th(m)*dz_1(n)
            else if (lcylindrical_coords) then
              dline_1(:,1)=dx_1(l1:l2)
              dline_1(:,2)=rcyl_mn1*dy_1(m)
              dline_1(:,3)=dz_1(n)
            else if (lcartesian_coords) then
              dline_1(:,1)=dx_1(l1:l2)
              dline_1(:,2)=dy_1(m)
              dline_1(:,3)=dz_1(n)
            endif
            dxmax_pencil = 0.
            if (nxgrid /= 1) dxmax_pencil=1./dx_1(l1:l2)
            if (nygrid /= 1) dxmax_pencil=max(1./dy_1(m),dxmax_pencil)
            if (nzgrid /= 1) dxmax_pencil=max(1./dz_1(n),dxmax_pencil)
!
            dxyz_2 = dline_1(:,1)**2+dline_1(:,2)**2+dline_1(:,3)**2
            dxyz_4 = dline_1(:,1)**4+dline_1(:,2)**4+dline_1(:,3)**4
            dxyz_6 = dline_1(:,1)**6+dline_1(:,2)**6+dline_1(:,3)**6
          endif
        endif
!
!  [AB: Isn't it true that not all 2-D averages use rcyl_mn?
!  lwrite_phiaverages=T is required, and perhaps only that.]
!  [BD: add also the z_mn dependency]
!
        if (l2davgfirst) then
          lpencil(i_rcyl_mn)=.true.
          lpencil(i_z_mn)=.true.
        endif
!
!  Calculate grid/geometry related pencils.
!
        call calc_pencils_grid(f,p)
!
!  Calculate profile for phi-averages if needed.
!
        if ((l2davgfirst.and.lwrite_phiaverages )  .or. &
            (l1dphiavg  .and.lwrite_phizaverages))  &
            call calc_phiavg_profile(p)
!
!  Calculate pencils for the pencil_case.
!  Note: some no-modules (e.g. nohydro) also calculate some pencils,
!  so it would be wrong to check for lhydro etc in such cases.
!
! To check ghost cell consistency, please uncomment the following 2 lines:
!       if (.not. lpencil_check_at_work .and. necessary(imn)) &
!       call check_ghosts_consistency (f, 'before calc_pencils_*')
                              call calc_pencils_hydro(f,p)
                              call calc_pencils_density(f,p)
        if (lpscalar)         call calc_pencils_pscalar(f,p)
                              call calc_pencils_eos(f,p)
        if (lshock)           call calc_pencils_shock(f,p)
        if (lchemistry)       call calc_pencils_chemistry(f,p)
        if (lviscosity)       call calc_pencils_viscosity(f,p)
        if (lforcing_cont)    call calc_pencils_forcing(f,p)
                              call calc_pencils_entropy(f,p)
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
        if (lspecial)         call calc_pencils_special(f,p)
        if (lborder_profiles) call calc_pencils_borderprofiles(f,p)
        if (lparticles)       call particles_calc_pencils(f,p)
!
!  --------------------------------------------------------
!  NO CALLS MODIFYING PENCIL_CASE PENCILS BEYOND THIS POINT
!  --------------------------------------------------------
!
!  hydro, density, and entropy evolution
!  Note that pressure gradient is added in dss_dt to momentum,
!  even if lentropy=.false.
!
        call duu_dt(f,df,p)
        call dlnrho_dt(f,df,p)
        call dss_dt(f,df,p)
!
!  Magnetic field evolution
!
        if (lmagnetic) call daa_dt(f,df,p)
!
!  Lorenz gauge evolution
!
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
          if (lhydro.or.ldustvelocity) call duu_dt_grav(f,df,p)
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
!  Evolution of radiative energy
!
        if (lradiation_fld) call de_dt(f,df,p,gamma)
!
!  Evolution of chemical species
!
        if (lchemistry) call dchemistry_dt(f,df,p)
!
!  Continuous forcing function (currently only for extra diagonstics)
!
        if (lforcing_cont) call forcing_continuous(df,p)
!
!  Add and extra 'special' physics
!
        if (lspecial) call dspecial_dt(f,df,p)
!
!  Add radiative cooling and radiative pressure (for ray method)
!
        if (lradiation_ray.and.(lentropy.or.ltemperature)) then
          call radiative_cooling(f,df,p)
          call radiative_pressure(f,df,p)
        endif
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
!  Call diagnostics that involves the full right hand side
!  This must be done at the end of all calls that might modify df.
!
        if (ldiagnos) then
          if (lmagnetic) call df_diagnos_magnetic(df,p)
        endif
!
!  General phiaverage quantities -- useful for debugging.
!
        if (l2davgfirst) then
          call phisum_mn_name_rz(p%rcyl_mn,idiag_rcylmphi)
          call phisum_mn_name_rz(p%phi_mn,idiag_phimphi)
          call phisum_mn_name_rz(p%z_mn,idiag_zmphi)
          call phisum_mn_name_rz(p%r_mn,idiag_rmphi)
        endif
!
!  Do the vorticity integration here, before the omega pencil is overwritten.
!
        if (ltime_integrals) then
          if (llast) then
            if (lhydro)    call time_integrals_hydro(f,p)
            if (lmagnetic) call time_integrals_magnetic(f,p)
          endif
        endif
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
        if (lfirst.and.ldt.and.(.not.ldt_paronly)) then
!
!  sum or maximum of the advection terms?
!  (lmaxadvec_sum=.false. by default)
!
          maxadvec=0.0
          maxdiffus=0.0
          maxdiffus2=0.0
          maxdiffus3=0.0
!--       if (lhydro) maxadvec=maxadvec+advec_uu
          if (lhydro.or.lhydro_kinematic) maxadvec=maxadvec+advec_uu
          if (lshear) maxadvec=maxadvec+advec_shear
          if (lneutralvelocity) maxadvec=maxadvec+advec_uun
          if (ldensity.or.lmagnetic.or.lradiation.or.lneutralvelocity) then
            advec2=0.0
            if (ldensity) advec2=advec2+advec_cs2
            if (lmagnetic) advec2=advec2+advec_va2
            if (lradiation) advec2=advec2+advec_crad2
            if (lneutralvelocity) advec2=advec2+advec_csn2
            if (lpolymer) advec2=advec2+advec_poly
            maxadvec=maxadvec+sqrt(advec2)
          endif
          if (ldensity.or.lhydro.or.lmagnetic.or.lenergy) then
            advec2_hypermesh=0.0
            if (ldensity)  advec2_hypermesh=advec2_hypermesh+advec_hypermesh_rho**2
            if (lhydro)    advec2_hypermesh=advec2_hypermesh+advec_hypermesh_uu**2
            if (lmagnetic) advec2_hypermesh=advec2_hypermesh+advec_hypermesh_aa**2
            if (lenergy)   advec2_hypermesh=advec2_hypermesh+advec_hypermesh_ss**2
            maxadvec=maxadvec+sqrt(advec2_hypermesh)
          endif
!
!  Time step constraints from each module.
!  (At the moment, magnetic and testfield use the same variable.)
!
          if (lviscosity)       maxdiffus=max(maxdiffus,diffus_nu)
          if (ldensity)         maxdiffus=max(maxdiffus,diffus_diffrho)
          if (lenergy)          maxdiffus=max(maxdiffus,diffus_chi)
          if (lmagnetic)        maxdiffus=max(maxdiffus,diffus_eta)
          if (ltestfield)       maxdiffus=max(maxdiffus,diffus_eta)
          if (ltestscalar)      maxdiffus=max(maxdiffus,diffus_eta)
          if (lpscalar)         maxdiffus=max(maxdiffus,diffus_pscalar)
          if (lcosmicray)       maxdiffus=max(maxdiffus,diffus_cr)
          if (ldustvelocity)    maxdiffus=max(maxdiffus,diffus_nud)
          if (ldustdensity)     maxdiffus=max(maxdiffus,diffus_diffnd)
          if (lchiral)          maxdiffus=max(maxdiffus,diffus_chiral)
          if (lchemistry)       maxdiffus=max(maxdiffus,diffus_chem)
          if (lneutralvelocity) maxdiffus=max(maxdiffus,diffus_nun)
          if (lneutraldensity)  maxdiffus=max(maxdiffus,diffus_diffrhon)
!
          if (lviscosity)       maxdiffus2=max(maxdiffus2,diffus_nu2)
          if (lmagnetic)        maxdiffus2=max(maxdiffus2,diffus_eta2)
!
          if (lviscosity)       maxdiffus3=max(maxdiffus3,diffus_nu3)
          if (ldensity)         maxdiffus3=max(maxdiffus3,diffus_diffrho3)
          if (lmagnetic)        maxdiffus3=max(maxdiffus3,diffus_eta3)
          if (ltestfield)       maxdiffus3=max(maxdiffus3,diffus_eta3)
          if (lenergy)          maxdiffus3=max(maxdiffus3,diffus_chi3)
          if (ldustvelocity)    maxdiffus3=max(maxdiffus3,diffus_nud3)
          if (ldustdensity)     maxdiffus3=max(maxdiffus3,diffus_diffnd3)
          if (lpscalar)         maxdiffus3=max(maxdiffus3,diffus_pscalar3)
          if (lneutralvelocity) maxdiffus3=max(maxdiffus3,diffus_nun3)
          if (lneutraldensity)  maxdiffus3=max(maxdiffus3,diffus_diffrhon3)
!
!  Exclude the frozen zones from the time-step calculation.
!
          if (any(lfreeze_varint)) then
            if (lcylinder_in_a_box.or.lcylindrical_coords) then
              where (p%rcyl_mn<=rfreeze_int)
                maxadvec=0.0
                maxdiffus=0.0
              endwhere
            else
              where (p%r_mn<=rfreeze_int)
                maxadvec=0.0
                maxdiffus=0.0
              endwhere
            endif
          endif
!
          if (any(lfreeze_varext)) then
            if (lcylinder_in_a_box.or.lcylindrical_coords) then
              where (p%rcyl_mn>=rfreeze_ext)
                maxadvec=0.0
                maxdiffus=0.0
              endwhere
            else
              where (p%r_mn>=rfreeze_ext)
                maxadvec=0.0
                maxdiffus=0.0
              endwhere
            endif
          endif
!
!  cdt, cdtv, and cdtc are empirical non-dimensional coefficients
!
          dt1_advec  = maxadvec/cdt
          dt1_diffus = maxdiffus/cdtv + maxdiffus2/cdtv2 + maxdiffus3/cdtv3
          dt1_max    = max(dt1_max,sqrt(dt1_advec**2+dt1_diffus**2))
!
!  time step constraint from speed of chemical reactions
!
          if (lchemistry .and. .not.llsode) then
            dt1_reac = reac_chem/cdtc
            dt1_max = max(dt1_max,dt1_reac)
          endif
!
!  time step constraint from relaxation time of polymer
!
          if (lpolymer) then
            dt1_poly_relax = trelax_poly/cdt_poly
            dt1_max = max(dt1_max,dt1_poly_relax)
          endif
!
!  Diagnostics showing how close to advective and diffusive time steps we are
!
          if (ldiagnos.and.idiag_dtv/=0) then
            call max_mn_name(maxadvec/cdt,idiag_dtv,l_dt=.true.)
          endif
          if (ldiagnos.and.idiag_dtdiffus/=0) then
            call max_mn_name(maxdiffus/cdtv,idiag_dtdiffus,l_dt=.true.)
          endif
!
!  Regular and hyperdiffusive mesh Reynolds numbers
!
          if (ldiagnos) then
            if (idiag_Rmesh/=0) &
                call max_mn_name(pi_1*maxadvec/(maxdiffus+tini),idiag_Rmesh)
            if (idiag_Rmesh3/=0) &
                call max_mn_name(pi5_1*maxadvec/(maxdiffus3+tini),idiag_Rmesh3)
            if (idiag_maxadvec/=0) &
                call max_mn_name(maxadvec,idiag_maxadvec)
          endif
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
          call calc_pencils_density(f,p)
          f(l1:l2,m,n,irhs) = p%rho*df(l1:l2,m,n,iuu)
          f(l1:l2,m,n,irhs+1) = p%rho*df(l1:l2,m,n,iuu+1)
          f(l1:l2,m,n,irhs+2) = p%rho*df(l1:l2,m,n,iuu+2)
          df(l1:l2,m,n,iuu:iuu+2) = df_iuu_pencil(1:nx,1:3) +&
                                    df(l1:l2,m,n,iuu:iuu+2)
          call sum_mn(p%rho,mass_per_proc(1))
        endif
        call timing('pde','end of mn loop',mnloop=.true.)
!
!  End of loops over m and n.
!
        headtt=.false.
      enddo mn_loop
      call timing('pde','at the end of the mn_loop')
!
!  Finish the job for the anelastic approximation
!
      if (lanelastic) call anelastic_after_mn(f, p, df, mass_per_proc)
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
      if (lparticles) call particles_pde_blocks(f,df)
      if (lparticles) call particles_pde(f,df)
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
!  Currently only needed for smoothed out particle drag force.
!
      if (lhydro .and. lfold_df) call fold_df(df,iux,iuz)
!
!  -------------------------------------------------------------
!  NO CALLS MODIFYING DF BEYOND THIS POINT (APART FROM FREEZING)
!  -------------------------------------------------------------
!
!  Freezing must be done after the full (m,n) loop, as df may be modified
!  outside of the considered pencil.
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
!
!  Recalculate grid/geometry related pencils. The r_mn and rcyl_mn are requested
!  in pencil_criteria_grid. Unfortunately we need to recalculate them here.
!
        if (any(lfreeze_varext).or.any(lfreeze_varint)) &
            call calc_pencils_grid(f,p)
!
!  Set df=0 for r_mn<r_int.
!
        if (any(lfreeze_varint)) then
          if (headtt) print*, 'pde: freezing variables for r < ', rfreeze_int, &
              ' : ', lfreeze_varint
          if (lcylinder_in_a_box.or.lcylindrical_coords) then
            if (wfreeze_int==0.0) then
              where (p%rcyl_mn<=rfreeze_int) pfreeze_int=0.0
              where (p%rcyl_mn> rfreeze_int) pfreeze_int=1.0
            else
              pfreeze_int=quintic_step(p%rcyl_mn,rfreeze_int,wfreeze_int, &
                  SHIFT=fshift_int)
            endif
          else
            if (wfreeze_int==0.0) then
              where (p%r_mn<=rfreeze_int) pfreeze_int=0.0
              where (p%r_mn> rfreeze_int) pfreeze_int=1.0
            else
              pfreeze_int=quintic_step(p%r_mn   ,rfreeze_int,wfreeze_int, &
                  SHIFT=fshift_int)
            endif
          endif
!
          do iv=1,nvar
            if (lfreeze_varint(iv)) &
                df(l1:l2,m,n,iv)=pfreeze_int*df(l1:l2,m,n,iv)
          enddo
!
        endif
!
!  Set df=0 for r_mn>r_ext.
!
        if (any(lfreeze_varext)) then
          if (headtt) print*, 'pde: freezing variables for r > ', rfreeze_ext, &
              ' : ', lfreeze_varext
          if (lcylinder_in_a_box) then
            if (wfreeze_ext==0.0) then
              where (p%rcyl_mn>=rfreeze_ext) pfreeze_ext=0.0
              where (p%rcyl_mn< rfreeze_ext) pfreeze_ext=1.0
            else
              pfreeze_ext=1.0-quintic_step(p%rcyl_mn,rfreeze_ext,wfreeze_ext, &
                SHIFT=fshift_ext)
            endif
          else
            if (wfreeze_ext==0.0) then
              where (p%r_mn>=rfreeze_ext) pfreeze_ext=0.0
              where (p%r_mn< rfreeze_ext) pfreeze_ext=1.0
            else
              pfreeze_ext=1.0-quintic_step(p%r_mn   ,rfreeze_ext,wfreeze_ext, &
                  SHIFT=fshift_ext)
            endif
          endif
!
          do iv=1,nvar
            if (lfreeze_varext(iv)) &
                df(l1:l2,m,n,iv) = pfreeze_ext*df(l1:l2,m,n,iv)
          enddo
        endif
!
!  Set df=0 inside square.
!
        if (any(lfreeze_varsquare)) then
          if (headtt) print*, 'pde: freezing variables inside square : ', &
              lfreeze_varsquare
          pfreeze=1.0-quintic_step(x(l1:l2),xfreeze_square,wfreeze,SHIFT=-1.0)*&
              quintic_step(spread(y(m),1,nx),yfreeze_square,-wfreeze,SHIFT=-1.0)
!
          do iv=1,nvar
            if (lfreeze_varsquare(iv)) &
                df(l1:l2,m,n,iv) = pfreeze*df(l1:l2,m,n,iv)
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
!  Only need to do this for nonperiodic x direction, on left/right-most
!  processor and in left/right--most pencils
!
          if (.not. lperi(1)) then
            if (lfirst_proc_x) then
              do iv=1,nvar
                if (lfrozen_bot_var_x(iv)) df(l1,m,n,iv) = 0.
              enddo
            endif
            if (llast_proc_x) then
              do iv=1,nvar
                if (lfrozen_top_var_x(iv)) df(l2,m,n,iv) = 0.
              enddo
            endif
          endif
!
        endif
!
!  Freezing boundary conditions in y.
!
        if (lfrozen_bcs_y) then ! are there any frozen vars at all?
!
!  Only need to do this for nonperiodic y direction, on bottom/top-most
!  processor and in bottom/top-most pencils.
!
          if (.not. lperi(2)) then
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
        endif
!
!  Freezing boundary conditions in z.
!
        if (lfrozen_bcs_z) then ! are there any frozen vars at all?
!
!  Only need to do this for nonperiodic z direction, on bottom/top-most
!  processor and in bottom/top-most pencils.
!
          if (.not. lperi(3)) then
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
        endif
!
!  Set df=0 for all solid cells.
!
        call freeze_solid_cells(df)
!
      enddo
!
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
!
!  Check for NaNs in the advection time-step.
!
      if (notanumber(dt1_advec)) then
        print*, 'pde: dt1_advec contains a NaN at iproc=', iproc
        if (lhydro)           print*, 'advec_uu   =',advec_uu
        if (lshear)           print*, 'advec_shear=',advec_shear
        if (lmagnetic)        print*, 'advec_hall =',advec_hall
        if (lneutralvelocity) print*, 'advec_uun  =',advec_uun
        if (lentropy)         print*, 'advec_cs2  =',advec_cs2
        if (lmagnetic)        print*, 'advec_va2  =',advec_va2
        if (lradiation)       print*, 'advec_crad2=',advec_crad2
        if (lneutralvelocity) print*, 'advec_csn2 =',advec_csn2
        if (lpolymer)         print*, 'advec_poly =',advec_poly
        call fatal_error_local('pde','')
      endif
!
!  Diagnostics.
!
      if (ldiagnos) call diagnostic(fname,nname)
!
!  1-D diagnostics.
!
      if (l1davgfirst) then
        if (lwrite_xyaverages) call xyaverages_z
        if (lwrite_xzaverages) call xzaverages_y
        if (lwrite_yzaverages) call yzaverages_x
      endif
      if (l1dphiavg) call phizaverages_r
!
!  2-D averages.
!
      if (l2davgfirst) then
        if (lwrite_yaverages)   call yaverages_xz
        if (lwrite_zaverages)   call zaverages_xy
        if (lwrite_phiaverages) call phiaverages_rz
      endif
!
!  Note: zaverages_xy are also needed if bmx and bmy are to be calculated
!  (of course, yaverages_xz does not need to be calculated for that).
!
      if (.not.l2davgfirst.and.ldiagnos.and.ldiagnos_need_zaverages) then
        if (lwrite_zaverages) call zaverages_xy
      endif
!
!  Calculate mean fields and diagnostics related to mean fields.
!
      if (ldiagnos) then
        if (lmagnetic) call calc_mfield
        if (lhydro)    call calc_mflow
        if (lpscalar)  call calc_mpscalar
      endif
!
!  Calculate rhoccm and cc2m (this requires that these are set in print.in).
!  Broadcast result to other processors. This is needed for calculating PDFs.
!
!      if (idiag_rhoccm/=0) then
!        if (iproc==0) rhoccm=fname(idiag_rhoccm)
!        call mpibcast_real(rhoccm,1)
!      endif
!
!      if (idiag_cc2m/=0) then
!        if (iproc==0) cc2m=fname(idiag_cc2m)
!        call mpibcast_real(cc2m,1)
!      endif
!
!      if (idiag_gcc2m/=0) then
!        if (iproc==0) gcc2m=fname(idiag_gcc2m)
!        call mpibcast_real(gcc2m,1)
!      endif
!
!  Reset lwrite_prof.
!
      lwrite_prof=.false.
!
    endsubroutine pde
!***********************************************************************
    subroutine debug_imn_arrays
!
!  For debug purposes: writes out the mm, nn, and necessary arrays.
!
!  23-nov-02/axel: coded
!
      open(1,file=trim(directory)//'/imn_arrays.dat')
      do imn=1,ny*nz
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
        call wsnap(trim(directory_snap)//'/'//filename,f,mvar_io,ENUM=.false.)
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
endmodule Equ
