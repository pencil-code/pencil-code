! $Id$
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
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
  interface add_interpolated
    module procedure add_interpolated_2D
    module procedure add_interpolated_3D
    module procedure add_interpolated_4D
  endinterface
!
  ! maximum number of granulation levels, technical maximum is 9
  integer, parameter :: max_gran_levels=3
  integer :: cool_RTV_cutoff=0
!
  real :: nc_tau=0., nc_alt=0., Kgpara=0., heat_cool=0., rho_diff_fix=0., cool_RTV=0., Kgpara2=0., nc_tau_rho=0., nc_alt_rho=0.
  real :: lntt0=0., wlntt=0., bmdi=0., hcond1=0., heatexp=0., heatamp=0., Ksat=0., Kc=0.
  real :: T_crit=0., deltaT_crit=0.
  real :: diffrho_hyper3=0., chi_hyper3=0., chi_hyper2=0., K_iso=0., b_tau=0., flux_tau=0.
  real :: Bavoid=0., vorticity_factor=5., tau_inv=1., Bz_flux=0., q0=1., qw=1., dq=0.1, dt_gran=0.
  logical :: lgranulation=.false., lgran_proc=.false., lgran_parallel=.false.
  logical :: luse_vel_field=.false., lquench=.false., lmassflux=.false.
  logical :: luse_mag_field=.false., luse_mag_vel_field=.false.
  logical :: luse_timedep_magnetogram=.false., lwrite_driver=.false.
  logical :: lnc_density_depend=.false., lnc_intrin_energy_depend=.false.
  logical :: lflux_emerg_bottom=.false.,lslope_limited_special=.false.,&
             lemerg_profx=.false.,lset_boundary_emf=.false.,lset_hotplate_lnTT=.false.,&
             lconv_vel_set_to_zero=.false.,ldispersion_acoustic=.false.
  integer :: irefz=n1, nglevel=max_gran_levels, cool_type=5
  real :: massflux=0., u_add
  real :: K_spitzer=0., hcond2=0., hcond3=0., init_time=0., init_time_hcond=0.
  real :: init_time_fade_start=0.0, init_time_hcond_fade_start=0.0
  real :: nc_z_max=0.0, nc_z_trans_width=0.0, nc_z_min=0.0, nc_z_min_trans_width=0.0
  real :: nc_lnrho_num_magn=0.0, nc_lnrho_trans_width=0.0
  real :: vel_time_offset=0.0, mag_time_offset=0.0
  real :: swamp_fade_start=0.0, swamp_fade_end=0.0
  real :: swamp_diffrho=0.0, swamp_chi=0.0, swamp_eta=0.0
  real :: lnrho_min=-max_real, lnrho_min_tau=1.0,uu_tau1_quench=0.0, lnTT_hotplate_tau=1.0
  real :: lnrho_max=max_real
  real, dimension(2) :: nwave=(/0.,0./),w_ff=(/0.,0./),z_ff=(/0.,0./)
  real, dimension(nx) :: glnTT_H, hlnTT_Bij, glnTT2, glnTT_abs, glnTT_abs_inv, glnTT_b
!
  integer :: nlf=4
!  real, dimension (mx) :: x12
!  real, dimension (my) :: y12
!  real, dimension (mz) :: z12
!
  real, dimension(nx,ny,2) :: A_init
!
  character(len=labellen), dimension(3) :: iheattype='nothing'
  real, dimension(2) :: heat_par_exp=(/0.,1./)
  real, dimension(2) :: heat_par_exp2=(/0.,1./)
  real, dimension(3) :: heat_par_gauss=(/0.,1.,0./)
!
  character(len=labellen) :: prof_type='nothing'
  real, dimension(mz) :: uu_init_z, lnrho_init_z, lnTT_init_z
  real, dimension(3) :: uu_emerg=0.0, bb_emerg=0.0, uu_drive=0.0
  real, dimension(:), allocatable :: deltaT_init_z, deltaE_init_z, E_init_z, deltarho_init_z, deltaH_part_init_z
  logical :: linit_uu=.false., linit_lnrho=.false., linit_lnTT=.false.
  logical :: lheatcond_cutoff=.false.
!
  ! file location settings
  character(len=*), parameter :: vel_times_dat = 'driver/vel_times.dat'
  character(len=*), parameter :: vel_field_dat = 'driver/vel_field.dat'
  character(len=*), parameter :: mag_times_dat = 'driver/mag_times.dat'
  character(len=*), parameter :: mag_field_dat = 'driver/mag_field.dat'
  character(len=*), parameter :: mag_vel_field_dat = 'driver/mag_vel_field.dat'
  character(len=labellen) :: flux_type='uniform'
!
  ! input parameters
  namelist /special_init_pars/ linit_uu,linit_lnrho,linit_lnTT,prof_type, &
           lslope_limited_special
!
  ! run parameters
  namelist /special_run_pars/ &
      nc_tau,nc_alt,Kgpara,heat_cool,rho_diff_fix,cool_RTV,lntt0,wlntt,bmdi,hcond1, &
      Kgpara2,K_spitzer,nc_tau_rho,nc_alt_rho,heatexp,heatamp,Ksat,Kc,diffrho_hyper3, &
      chi_hyper3,chi_hyper2,K_iso,lgranulation,lgran_parallel,irefz,tau_inv, &
      b_tau,flux_tau,Bavoid,nglevel,vorticity_factor,Bz_flux,init_time,init_time_hcond, &
      lquench,q0,qw,dq,massflux,luse_vel_field,luse_mag_vel_field,prof_type, &
      lmassflux,hcond2,hcond3,heat_par_gauss,heat_par_exp,heat_par_exp2, &
      iheattype,dt_gran,cool_type,luse_timedep_magnetogram,lwrite_driver, &
      nc_z_max,nc_z_min,nc_z_trans_width,nc_lnrho_num_magn,nc_lnrho_trans_width, nc_z_min_trans_width, &
      lnc_density_depend, lnc_intrin_energy_depend, &
      init_time_fade_start, init_time_hcond_fade_start, &
      swamp_fade_start, swamp_fade_end, swamp_diffrho, swamp_chi, swamp_eta, &
      vel_time_offset, mag_time_offset, lnrho_min, lnrho_min_tau, &
      cool_RTV_cutoff, T_crit, deltaT_crit, & 
      lflux_emerg_bottom, uu_emerg, bb_emerg, uu_drive,flux_type,lslope_limited_special, &
      lemerg_profx,lheatcond_cutoff,nwave,w_ff,z_ff,lset_boundary_emf,uu_tau1_quench, &
      lset_hotplate_lnTT,lnTT_hotplate_tau,lconv_vel_set_to_zero,ldispersion_acoustic, &
      lnrho_max
!
  integer :: ispecaux=0
  integer :: idiag_dtvel=0     ! DIAG_DOC: Velocity driver time step
  integer :: idiag_dtnewt=0    ! DIAG_DOC: Radiative cooling time step
  integer :: idiag_dtradloss=0 ! DIAG_DOC: Radiative losses time step
  integer :: idiag_dtchi2=0    ! DIAG_DOC: $\delta t / [c_{\delta t,{\rm v}}\,
                               ! DIAG_DOC:   \delta x^2/\chi_{\rm max}]$
                               ! DIAG_DOC:   \quad(time step relative to time
                               ! DIAG_DOC:   step based on heat conductivity;
                               ! DIAG_DOC:   see \S~\ref{time-step})
  integer :: idiag_dtspitzer=0 ! DIAG_DOC: Spitzer heat conduction time step
  integer :: idiag_mag_flux=0  ! DIAG_DOC: Total vertical magnetic flux at
                               ! bottom boundary: mag_flux=sum(|Bz(n1)|)*(dx*dy)
!
  integer :: ivid_logQ=0, ivid_rtv=0
!
  ! video slices
  real, target, dimension(:,:), allocatable :: rtv_xy, rtv_xy2, rtv_xy3, rtv_xy4
  real, target, dimension(:,:), allocatable :: rtv_xz, rtv_yz, rtv_xz2
  real, target, dimension(:,:,:,:,:), allocatable :: rtv_r
  real, target, dimension(:,:), allocatable :: logQ_xy, logQ_xy2, logQ_xy3, logQ_xy4
  real, target, dimension(:,:), allocatable :: logQ_xz, logQ_xz2, logQ_yz
  real, target, dimension(:,:,:,:,:), allocatable :: logQ_r
!
  ! Granule midpoint:
  type point
    ! X and Y position
    real :: pos_x, pos_y
    ! Current and maximum amplitude
    real :: amp, amp_max
    ! Time of amplitude maximum and total lifetime
    real :: t_amp_max, t_life
    type (point), pointer :: next
    type (point), pointer :: previous
  endtype point
!
  type (point), pointer :: first => null()
  type (point), pointer :: current => null()
!
  ! List of start pointers for each granulation level:
  type gran_list_start
    type (point), pointer :: first
  endtype gran_list_start
  type (gran_list_start), dimension(max_gran_levels) :: gran_list
!
  integer :: xrange, yrange, pow
  real :: ampl, dxdy2, ig, granr, pd, life_t, avoid
  real, dimension(:,:), allocatable :: w, vx, vy
  real, dimension(:,:), allocatable :: Ux, Uy
  real, dimension(:,:), allocatable :: Ux_ext, Uy_ext
  real, dimension(:,:), allocatable :: BB2
  real, dimension(:,:), allocatable :: BB2_local
  integer, dimension(:,:), allocatable :: avoid_gran
  real, save :: tsnap_uu=0., thresh
  integer, save :: isnap
  integer, save, dimension(mseed) :: points_rstate
  real, dimension(nx,ny), save :: Ux_local, Uy_local
!
  integer, save, dimension(mseed) :: nano_seed
  integer :: alloc_err
  real, dimension(nx) :: diffus_chi, diffus_chi3
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      use FArrayManager
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
! Called with lreloading indicating a RELOAD
!
!  06-oct-03/tony: coded
!
      use File_io, only: parallel_file_exists
      use Slices_methods, only: alloc_slice_buffers
!
      real, dimension(mx,my,mz,mfarray) :: f
!
! Consistency checks:
!
      ! Restoration half-time of initial magnetic field:
      if ((bmdi > 0.0) .and. (b_tau > 0.0) .and. (b_tau /= bmdi)) &
          call fatal_error ('solar_corona/mag_driver', &
              "Please use either 'b_tau' or 'bmdi', not both.")
      if ((b_tau == 0.0) .and. (bmdi == 0.0) .and. luse_mag_field) &
          call fatal_error ('solar_corona/mag_driver', &
              "With 'luse_mag_field', please set either 'b_tau' or 'bmdi'.")
      if (irefz > max (n1, nz)) &
          call fatal_error ('solar_corona', &
              "'irefz' must lie in the bottom layer of processors.")
      if (((nglevel < 1) .or. (nglevel > max_gran_levels))) &
          call fatal_error ('solar_corona/gran_driver', &
              "'nglevel' is invalid and/or larger than 'max_gran_levels'.")
      ! Using vorticity increase, the x- and y-directions must be equidistant.
      if (lgranulation .and. (vorticity_factor > 0.0) .and. ((dx /= dy) .or. any (.not. lequidist(1:2)))) &
          call fatal_error ('solar_corona/gran_driver', &
              "If 'vorticity_factor' is active, the grid must be equidistant in x and y.")
      ! For only one granulation level, no parallelization is required.
      if (lgran_parallel .and. (nglevel == 1)) &
          call fatal_error ('solar_corona/gran_driver', &
              "If 'nglevel' is 1, 'lgran_parallel' should be deactivated.")
      ! If not at least 3 procs above the ipz=0 plane are available,
      ! computing of granular velocities has to be done non-parallel.
      if (lgran_parallel .and. ((nprocz-1)*nprocxy < nglevel)) &
          call fatal_error ('solar_corona/gran_driver', &
              "There are not enough processors to activate 'lgran_parallel'.")
      if (lquench .and. (.not. lmagnetic)) &
          call fatal_error ('solar_corona', &
              "Quenching needs the magnetic module.")
      if ((Bavoid > 0.0) .and. (.not. lmagnetic)) &
          call fatal_error ('solar_corona', &
              "'Bavoid' needs the magnetic module.")
      if ((nc_tau > 0.0) .and. (nc_alt == 0.0) .and. (nc_lnrho_num_magn == 0.0)) &
          call fatal_error ('solar_corona', &
              "Please select decaying of Newton Cooling using 'nc_alt' or 'nc_lnrho_num_magn'.")
      ! Restoration half-time of initial total vertical flux:
      if ((flux_tau > 0.0) .and. (Bz_flux <= 0.0)) &
          call fatal_error ('solar_corona/mag_driver', &
              "Together with 'flux_tau', 'Bz_flux' needs to be set and positive.")
      ! Check if heat conduction terms are implemented:
      if ((K_spitzer /= 0.0) .and. (lentropy .or. (ltemperature .and. ltemperature_nolog))) &
          call fatal_error('solar_corona/calc_heatcond_tensor', &
              "Heat conduction 'K_spitzer' currently requirees logarithmic temperature.", .true.)
      if ((K_iso /= 0.0) .and. lentropy) &
          call fatal_error('solar_corona/calc_heatcond_grad', &
              "Heat conduction 'K_iso' is currently not implemented for entropy.", .true.)
      if (((nc_tau > 0.0) .or. (nc_alt /= 0.0) .or. (nc_lnrho_num_magn /= 0.0)) .and. lentropy) &
          call fatal_error('solar_corona/calc_heat_cool_newton', &
              "Newton cooling is currently not implemented for entropy (see 'nc_tau', 'nc_alt', 'nc_lnrho_num_magn').", .true.)
!
      if (lbfield) then
        if (luse_mag_field) call fatal_error("solar_corona", "luse_mag_field not implemented with bfield")
        if (lset_boundary_emf) call fatal_error("solar_corona", "lset_boundary_emf not implemented with bfield")
        if (lflux_emerg_bottom) call fatal_error("solar_corona", "lflux_emerg_bottom not implemented with bfield")
        if (swamp_eta > 0.0) call fatal_error("solar_corona", "not implemented with bfield. ")
      endif
!
      if ((.not. lreloading) .and. lrun) nano_seed = 0
!
      if (lpencil_check_at_work) return
!
      if (luse_timedep_magnetogram .and. .not. parallel_file_exists (mag_times_dat)) &
          call fatal_error ('solar_corona/mag_driver', &
              "Could not find file '"//trim (mag_times_dat)//"'.")
!
      if (lgranulation .and. lrun) then
        ! Define and initialize the processors that are computing the granulation
        lgran_proc = ((.not. lgran_parallel) .and. lroot) .or. &
            (lgran_parallel .and. (iproc >= nprocxy) .and. (iproc < nprocxy+nglevel))
        if (lroot .or. lgran_proc) call set_gran_params()
      endif
!
      ! Setup atmosphere stratification for later use
      call setup_profiles()
!
      if (lslope_limited_special) then
        if (lroot) print*,'initialize_special: Set up half grid x12, y12, z12'
        call generate_halfgrid(x12,y12,z12)
      endif
!
      if (ivid_rtv/=0) &
        call alloc_slice_buffers(rtv_xy,rtv_xz,rtv_yz,rtv_xy2,rtv_xy3,rtv_xy4,rtv_xz2,rtv_r)
      if (ivid_logQ/=0) &
        call alloc_slice_buffers(logQ_xy,logQ_xz,logQ_yz,logQ_xy2,logQ_xy3,logQ_xy4,logQ_xz2,logQ_r)

      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  Initialize special condition; called by start.f90.
!
!  27-aug-2010/Bourdin.KIS: coded
!
      use EquationOfState, only: lnrho0,gamma,gamma_m1,cs20,cs2top,cs2bot
      use Messages, only: warning
!
      real, dimension(mx,my,mz,mfarray), intent(out) :: f
!
      integer :: j
!
      if (linit_uu) then
        ! set initial vertical velocity profile values
        do j = 1, mz
          f(:,:,j,iuz) = uu_init_z(j)
        enddo
      endif
!
      if (linit_lnrho) then
        if (.not. ldensity) call fatal_error ('init_special', &
            "linit_lnrho=T needs the density module")
        ! set initial density profile values
        do j = 1, mz
          if (ldensity .and. ldensity_nolog) then
            f(:,:,j,irho) = exp (lnrho_init_z(j))
          elseif (ldensity) then
            f(:,:,j,ilnrho) = lnrho_init_z(j)
          endif
        enddo
      endif
!
      if (linit_lnTT) then
        if (pretend_lnTT) call fatal_error ('init_special', &
            "linit_lnTT=T not implemented for pretend_lnTT=T")
        ! set initial temperaure profile values
        do j = 1, mz
          if (ltemperature .and. ltemperature_nolog) then
            f(:,:,j,ilnTT) = exp (lnTT_init_z(j))
          elseif (ltemperature) then
            f(:,:,j,ilnTT) = lnTT_init_z(j)
          elseif (lentropy) then
            f(:,:,j,iss) = (alog(gamma_m1/cs20)+lnTT_init_z(j)- &
                gamma_m1*(f(l1,m1,j,ilnrho)-lnrho0))/gamma
          endif
        enddo
        ! set bottom and top boundary sound speed values
        cs2bot = gamma_m1*exp(lnTT_init_z(n1))
        cs2top = gamma_m1*exp(lnTT_init_z(n2))
      endif
!
    endsubroutine init_special
!***********************************************************************
    subroutine finalize_special(f)
!
!  Called right before exiting.
!
!  14-aug-2011/Bourdin.KIS: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      if (lgranulation .and. (lroot .or. lgran_proc) .and. lrun) then
        call free_points ()
        deallocate (Ux, Uy, w, vx, vy, avoid_gran)
      endif
      if (allocated (BB2)) deallocate (BB2)
      if (allocated (Ux_ext)) deallocate (Ux_ext)
      if (allocated (Uy_ext)) deallocate (Uy_ext)
!
      call special_before_boundary (f, .true.)
!
    endsubroutine finalize_special
!***********************************************************************
    subroutine setup_profiles()
!
!  Read and set vertical profiles for initial temperature and density.
!  Initial vertical velocity profile is given in [m/s] over z.
!  Initial density profile is given in ln(rho) [kg/m^3] over z.
!  Initial temperature profile is given in ln(T) [K] over z.
!  When using 'prof_ln*.dat' files, z is expected to be in SI units [m],
!  when using the 'stratification.dat' file, z is expected to be in [Mm].
!
!  25-aug-2010/Bourdin.KIS: coded
!
      use EquationOfState, only: lnrho0, rho0
      use Mpicomm, only: mpibcast_real
!
      logical :: lnewton_cooling=.false.
!
! Only read profiles if needed, e.g. for Newton cooling
!
      if (.not. (ltemperature .or. lentropy)) return
      lnewton_cooling = (nc_tau > 0.0) .or. (nc_tau_rho > 0.0)
      if (.not. (linit_uu .or. linit_lnrho .or. linit_lnTT .or. lnewton_cooling)) return
!
      ! default: read 'stratification.dat' with density and temperature
      if (prof_type == 'nothing') prof_type = 'lnrho_lnTT'
!
      ! check if density profile is read, when needed
      if ((nc_tau_rho > 0.0) .and. (index (prof_type, 'lnrho') < 1)) then
        call fatal_error ("setup_profiles", &
            "a density profile must be read to use density based newton cooling")
      endif
!
      ! read the profiles
      call read_profiles()
!
      if (linit_lnrho .or. lnewton_cooling) then
        ! some kind of density profile is actually in use,
        ! set lnrho0 and rho0 accordingly to the lower boundary value
        if (lroot) then
          if ((lnrho0 /= 0.0) .and. (abs(lnrho0 / lnrho_init_z(irefz) -1 ) > 1e-6)) then
            write (*,*) 'lnrho0 inconsistent: ', lnrho0, lnrho_init_z(irefz)
            call fatal_error ("setup_profiles", "conflicting manual lnrho0 setting", .true.)
          endif
          lnrho0 = lnrho_init_z(irefz)
          if ((rho0 /= 1.0) .and. (abs (rho0 / exp (lnrho0) - 1.0) > 1e-6)) then
            write (*,*) 'rho0 inconsistent: ', rho0, exp (lnrho0)
            call fatal_error ("setup_profiles", "conflicting manual rho0 setting", .true.)
          endif
          rho0 = exp (lnrho0)
        endif
        call mpibcast_real(lnrho0)
        call mpibcast_real (rho0)
      endif
!
    endsubroutine setup_profiles
!***********************************************************************
    subroutine read_profiles()
!
!  Read profiles for temperature, velocity, and/or density stratification.
!
!  21-oct-2010/Bourdin.KIS: coded
!
      use File_io, only: file_exists
      use Mpicomm, only: mpibcast_real
!
      integer :: i
      integer, parameter :: unit=12
      real :: var_lnrho, var_lnTT, var_z
      real, dimension(:), allocatable :: prof_lnrho, prof_lnTT, prof_z
      logical :: lread_prof_uu, lread_prof_lnrho, lread_prof_lnTT, lread_prof_deltaT, lread_prof_deltaE, lread_prof_deltarho, &
          lread_prof_E,lread_prof_deltaH_part
!
      ! file location settings
      character(len=*), parameter :: stratification_dat = 'stratification.dat'
      character(len=*), parameter :: lnrho_dat = 'driver/prof_lnrho.dat'
      character(len=*), parameter :: lnT_dat = 'driver/prof_lnT.dat'
      character(len=*), parameter :: uz_dat = 'driver/prof_uz.dat'
      character(len=*), parameter :: deltaT_dat = 'driver/prof_deltaT.dat'
      character(len=*), parameter :: deltaE_dat = 'driver/prof_deltaE.dat'
      character(len=*), parameter :: deltaH_part_dat = 'driver/prof_deltaH_part.dat'
      character(len=*), parameter :: E_dat = 'driver/prof_E.dat'
      character(len=*), parameter :: deltarho_dat = 'driver/prof_deltarho.dat'
!
! Check which stratification file should be used:
!
      lread_prof_uu = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_uu') > 0)
      lread_prof_lnrho = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_lnrho') > 0)
      lread_prof_lnTT = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_lnTT') > 0)
      lread_prof_deltaT = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_deltaT') > 0)
      lread_prof_deltaE = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_deltaE') > 0)
      lread_prof_deltaH_part = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_deltaH_part') > 0)
      lread_prof_E = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_E') > 0)
      lread_prof_deltarho = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_deltarho') > 0)
!
      if (prof_type == 'lnrho_lnTT') then
        allocate (prof_lnTT(nzgrid), prof_lnrho(nzgrid), prof_z(nzgrid), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('setup_profiles', &
            'Could not allocate memory for stratification variables', .true.)
!
        ! read stratification file only on the MPI root rank
        if (lroot) then
          if (.not. file_exists (stratification_dat)) call fatal_error ( &
              'setup_profiles', 'Stratification file not found', .true.)
          open (unit,file=stratification_dat)
          do i = 1, nzgrid
            read (unit,*) var_z, var_lnrho, var_lnTT
            prof_z(i) = var_z
            prof_lnTT(i) = var_lnTT
            prof_lnrho(i) = var_lnrho
          enddo
          close (unit)
        endif
!
        call mpibcast_real (prof_z, nzgrid)
        call mpibcast_real (prof_lnTT, nzgrid)
        call mpibcast_real (prof_lnrho, nzgrid)
!
        ! interpolate logarthmic data to Pencil grid profile
        call interpolate_profile (prof_lnTT, prof_z, nzgrid, lnTT_init_z)
        call interpolate_profile (prof_lnrho, prof_z, nzgrid, lnrho_init_z)
!
        deallocate (prof_lnTT, prof_lnrho, prof_z)
!
      elseif (lread_prof_uu .or. lread_prof_lnrho .or. lread_prof_lnTT .or. lread_prof_deltaT &
          .or. lread_prof_deltaE .or. lread_prof_deltarho .or. lread_prof_deltaH_part) then
!
        ! read vertical velocity profile for interpolation
        if (lread_prof_uu) &
            call read_profile (uz_dat, uu_init_z, real(unit_velocity), .false.)
!
        ! read logarithmic density profile for interpolation
        if (lread_prof_lnrho) &
            call read_profile (lnrho_dat, lnrho_init_z, real(unit_density), .true.)
!
        ! read logarithmic temperature profile
        if (lread_prof_lnTT) &
            call read_profile (lnT_dat, lnTT_init_z, real(unit_temperature), .true.)
!
        ! read temperature differences profile
        if (lread_prof_deltaT) then
          allocate (deltaT_init_z(mz))
          call read_profile (deltaT_dat, deltaT_init_z, real(unit_temperature/unit_time), .false.)
        endif
!
        ! read internal energy profile
        if (lread_prof_E) then
          allocate (E_init_z(mz))
          call read_profile (E_dat, E_init_z, real(unit_energy), .false.)
        endif
!
        ! read internal energy differences profile
        if (lread_prof_deltaE) then
          allocate (deltaE_init_z(mz))
          call read_profile (deltaE_dat, deltaE_init_z, real(unit_energy/unit_time), .false.)
        endif
!
        ! read heating per particle differences profile
        if (lread_prof_deltaH_part) then
          allocate (deltaH_part_init_z(mz))
          call read_profile (deltaH_part_dat, deltaH_part_init_z, real(unit_energy/unit_time), .false.)
        endif
        ! read density differences profile
        if (lread_prof_deltarho) then
          allocate (deltarho_init_z(mz))
          call read_profile (deltarho_dat, deltarho_init_z, real(unit_density/unit_time), .false.)
        endif
!
      elseif (index (prof_type, 'internal_') == 1) then
        call warning ('read_profiles', "using internal profile '"//trim(prof_type)//"'.")
      elseif (index (prof_type, 'initial_') == 1) then
        call fatal_error ('read_profiles', "prof_type='"//trim(prof_type)//"' is not yet implemented.")
      else
        call fatal_error ('read_profiles', "prof_type='"//trim(prof_type)//"' unknown.")
      endif
!
    endsubroutine read_profiles
!***********************************************************************
    subroutine read_profile(filename,profile,data_unit,llog)
!
!  Read vertical profile data.
!  Values are expected in SI units.
!
!  15-sept-2010/Bourdin.KIS: coded
!
      use File_io, only: parallel_file_exists, file_size
      use Mpicomm, only: mpibcast_int, mpibcast_real
!
      character(len=*), intent(in) :: filename
      real, dimension(mz), intent(out) :: profile
      real, intent(in) :: data_unit
      logical, intent(in) :: llog
!
      real, dimension(:), allocatable :: data, data_z
      integer :: n_data
!
      integer, parameter :: unit=12
      integer :: len_double
!
!
      inquire (IOLENGTH=len_double) 1.0d0
!
      if (.not. parallel_file_exists (filename)) &
          call fatal_error ('read_profile', "can't find "//filename)
      ! file access is only done on the MPI root rank
      if (lroot) then
        ! determine the number of data points in the profile
        n_data = (file_size (filename) - 2*2*4) / (len_double * 2)
      endif
      call mpibcast_int (n_data)
!
      ! allocate memory
      allocate (data(n_data), data_z(n_data), stat=alloc_err)
      if (alloc_err > 0) call fatal_error ('read_profile', &
          'Could not allocate memory for data and its z coordinate', .true.)
!
      if (lroot) then
        ! read profile
        open (unit, file=filename, form='unformatted', recl=len_double*n_data)
        read (unit) data
        read (unit) data_z
        close (unit)
!
        if (llog) then
          ! convert data from logarithmic SI to logarithmic Pencil units
          data = data - alog (data_unit)
        else
          ! convert data from SI to Pencil units
          data = data / data_unit
        endif
!
        ! convert z coordinates from SI to Pencil units
        data_z = data_z / unit_length
      endif
!
      ! broadcast profile
      call mpibcast_real (data, n_data)
      call mpibcast_real (data_z, n_data)
!
      ! interpolate logarthmic data to Pencil grid profile
      call interpolate_profile (data, data_z, n_data, profile)
!
      deallocate (data, data_z)
!
    endsubroutine read_profile
!***********************************************************************
    subroutine interpolate_profile(data,data_z,n_data,profile)
!
!  Interpolate profile data to Pencil grid.
!
!  15-sept-2010/Bourdin.KIS: coded
!
      use General, only: itoa
!
      integer :: n_data
      real, dimension(n_data), intent(in) :: data, data_z
      real, dimension(mz), intent(out) :: profile
!
      integer :: i, j, num_over, num_below
!
!
      ! linear interpolation of data
      num_below = 0
      num_over = 0
      do j = 1, mz
        if (z(j) < data_z(1) ) then
          ! extrapolate linarily below bottom
          num_below = num_below + 1
          profile(j) = data(1) + (data(2)-data(1))/(data_z(2)-data_z(1)) * (z(j)-data_z(1))
        elseif (z(j) >= data_z(n_data)) then
          ! extrapolate linarily over top
          num_over = num_over + 1
          profile(j) = data(n_data) + (data(n_data)-data(n_data-1))/(data_z(n_data)-data_z(n_data-1)) * (z(j)-data_z(n_data))
        else
          do i = 1, n_data-1
            if ((z(j) >= data_z(i)) .and. (z(j) < data_z(i+1))) then
              ! y = m*(x-x1) + y1
              profile(j) = (data(i+1)-data(i)) / (data_z(i+1)-data_z(i)) * (z(j)-data_z(i)) + data(i)
              exit
            endif
          enddo
        endif
      enddo
!
      if (lfirst_proc_xy .and. (num_below > 0)) then
        call warning ("interpolate_profile", &
            "extrapolated "//trim (itoa (num_below))//" grid points below bottom")
      endif
      if (lfirst_proc_xy .and. (num_over > 0)) then
        call warning ("interpolate_profile", &
            "extrapolated "//trim (itoa (num_over))//" grid points over top")
      endif
!
    endsubroutine interpolate_profile
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      integer :: i
      logical :: lprof_deltaH_part
!
      if (cool_RTV /= 0.0) then
        lpenc_requested(i_lnrho) = .true.
        lpenc_requested(i_lnTT) = .true.
        lpenc_requested(i_cp1) = .true.
      endif
!
      if (heat_cool /= 0.0) then
        lpenc_requested(i_lnrho) = .true.
        lpenc_requested(i_lnTT) = .true.
        lpenc_requested(i_cp1) = .true.
        lpenc_requested(i_rho1) = .true.
        lprof_deltaH_part = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_deltaH_part') > 0)
        if (lprof_deltaH_part) then
          lpenc_requested(i_TT1) = .true.
!         lpenc_requested(i_glnrho) = .true.
        endif
      endif
!
      if ((nc_tau > 0.0) .or. (nc_tau_rho > 0.0)) then
        lpenc_requested(i_lnrho) = .true.
        lpenc_requested(i_lnTT) = .true.
      endif
!
      if (swamp_diffrho > 0.0) then
        if (ldensity_nolog) then
          lpenc_requested(i_del2rho) = .true.
        else
          lpenc_requested(i_del2lnrho) = .true.
        endif
      endif
!
      if (swamp_chi > 0.0) then
        if (ltemperature) then
          if (ltemperature_nolog) then
            lpenc_requested(i_del2TT) = .true.
          else
            lpenc_requested(i_del2lnTT) = .true.
          endif
        elseif (lentropy) then
          lpenc_requested(i_del2ss) = .true.
        endif
      endif
!
      if (swamp_eta > 0.0) then
        lpenc_requested(i_del2a) = .true.
        lpenc_requested(i_diva) = .true.
      endif
!
      if (hcond1 /= 0.0) then
        lpenc_requested(i_b2) = .true.
        lpenc_requested(i_bij) = .true.
        lpenc_requested(i_bunit) = .true.
        lpenc_requested(i_lnTT) = .true.
        lpenc_requested(i_glnTT) = .true.
        lpenc_requested(i_hlnTT) = .true.
        lpenc_requested(i_del2lnTT) = .true.
        lpenc_requested(i_lnrho) = .true.
        lpenc_requested(i_glnrho) = .true.
      endif
!
      if (hcond2 /= 0.0) then
        lpenc_requested(i_b2) = .true.
        lpenc_requested(i_bij) = .true.
        lpenc_requested(i_bunit) = .true.
        lpenc_requested(i_glnTT) = .true.
        lpenc_requested(i_hlnTT) = .true.
        lpenc_requested(i_glnrho) = .true.
      endif
!
      if (hcond3 /= 0.0) then
        lpenc_requested(i_glnTT) = .true.
        lpenc_requested(i_hlnTT) = .true.
        lpenc_requested(i_del2lnTT) = .true.
        lpenc_requested(i_glnrho) = .true.
      endif
!
      if (K_iso /= 0.0) then
        lpenc_requested(i_glnrho) = .true.
        lpenc_requested(i_TT) = .true.
        lpenc_requested(i_lnTT) = .true.
        lpenc_requested(i_glnTT) = .true.
        lpenc_requested(i_hlnTT) = .true.
        lpenc_requested(i_del2lnTT) = .true.
      endif
!
      if (K_spitzer /= 0.0) then
        lpenc_requested(i_cp1) = .true.
        lpenc_requested(i_b2) = .true.
        lpenc_requested(i_bij) = .true.
        lpenc_requested(i_bunit) = .true.
        lpenc_requested(i_TT) = .true.
        lpenc_requested(i_glnTT) = .true.
        lpenc_requested(i_hlnTT) = .true.
        lpenc_requested(i_rho1) = .true.
        lpenc_requested(i_glnrho) = .true.
      endif
!
      if (Ksat /= 0.) then
        lpenc_requested(i_cp1) = .true.
        lpenc_requested(i_TT) = .true.
        lpenc_requested(i_glnTT) = .true.
        lpenc_requested(i_glnrho) = .true.
      endif
!
      if (Kc /= 0.) then
        lpenc_requested(i_cp1) = .true.
        lpenc_requested(i_lnrho) = .true.
        lpenc_requested(i_glnrho) = .true.
        lpenc_requested(i_glnTT) = .true.
      endif
!
      if (idiag_dtchi2 /= 0) then
        lpenc_diagnos(i_rho1) = .true.
        lpenc_diagnos(i_cv1) = .true.
        lpenc_diagnos(i_cs2) = .true.
      endif
!
      do i = 1,3
        select case(iheattype(i))
        case ('sven')
          lpenc_diagnos(i_cp1) = .true.
          lpenc_diagnos(i_TT1) = .true.
          lpenc_diagnos(i_rho1) = .true.
        case ('T<Tcrit')
          lpenc_requested(i_b2) = .true.
          lpenc_requested(i_TT) = .true.
          lpenc_requested(i_rho1) = .true.
          lpenc_requested(i_cp1) = .true.
        endselect
      enddo
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
      if (Kgpara /= 0.0) then
        call warning('calc_heatcond_grad', &
            'Please use K_spitzer instead of Kgpara')
        K_spitzer = Kgpara
      endif
!
      if (Kgpara2 /= 0.0) then
        if (K_iso /= 0.0) then
          call fatal_error('calc_heatcond_grad', &
              'Use only K_iso instead of Kgpara2')
        else
          call warning('calc_heatcond_grad', &
              'Please use K_iso instead of Kgpara2')
          K_iso = Kgpara2
        endif
      endif
!
      luse_mag_field = (b_tau > 0.0) .or. (bmdi > 0.0)
      ! Bz_flux is the sum of the absolute vertical flux provided in [T*m^2].
      ! After reading Bz_flux, convert SI units (from file) to Pencil units:
      Bz_flux = Bz_flux / (unit_magnetic * unit_length**2)
      ! In the 2D-case, units have to be [T*m] for backwards compatibility
      if ((nxgrid == 1) .or. (nygrid == 1)) Bz_flux = Bz_flux * unit_length
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
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Diagnostics, only: parse_name
      use FArrayManager, only: farray_index_append
!
      integer :: iname
      logical :: lreset, lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr = lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtvel = 0
        idiag_dtchi2 = 0
        idiag_dtnewt = 0
        idiag_dtradloss = 0
        idiag_dtspitzer = 0
        idiag_mag_flux = 0
        ivid_logQ=0
        ivid_rtv=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname = 1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtvel',idiag_dtvel)
        call parse_name(iname,cname(iname),cform(iname),'dtchi2',idiag_dtchi2)
        call parse_name(iname,cname(iname),cform(iname),'dtnewt',idiag_dtnewt)
        call parse_name(iname,cname(iname),cform(iname),'dtradloss',idiag_dtradloss)
        call parse_name(iname,cname(iname),cform(iname),'dtspitzer',idiag_dtspitzer)
        call parse_name(iname,cname(iname),cform(iname),'mag_flux',idiag_mag_flux)
      enddo
!
!  check for those quantities for which we want video slices
!
      do iname = 1,nnamev
        call parse_name(iname,cnamev(iname),cformv(iname),'logQ',ivid_logQ)
        call parse_name(iname,cnamev(iname),cformv(iname),'rtv',ivid_rtv)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        call farray_index_append('i_dtvel',idiag_dtvel)
        call farray_index_append('i_dtchi2',idiag_dtchi2)
        call farray_index_append('i_dtnewt',idiag_dtnewt)
        call farray_index_append('i_dtradloss',idiag_dtradloss)
        call farray_index_append('i_dtspitzer',idiag_dtspitzer)
        call farray_index_append('i_mag_flux',idiag_mag_flux)
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  26-jun-06/tony: dummy
! 
      use Slices_methods, only: assign_slices_scal

      real, dimension(mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
      case ('rtv'); call assign_slices_scal(slices,rtv_xy,rtv_xz,rtv_yz,rtv_xy2, &
                                            rtv_xy3,rtv_xy4,rtv_xz2,rtv_r)
!
      case ('logQ'); call assign_slices_scal(slices,logQ_xy,logQ_xz,logQ_yz,logQ_xy2, &
                                             logQ_xy3,logQ_xy4,logQ_xz2,logQ_r)
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      if (lgranulation .or. luse_vel_field) then
        ! Apply driving of velocity field.
        if (.not. lpencil_check_at_work) &
            call vel_driver (Ux_local, Uy_local, f, df)
      endif
!
      if (lmassflux) call force_solar_wind(df,p)
!
      if (swamp_eta > 0.0) call calc_swamp_eta(df,p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!  computes hyper diffusion for non equidistant grid
!  using the IGNOREDX keyword.
!  This includes an optional density floor;
!  should be moved to density.f90 at some point.
!
!  17-feb-10/bing: coded
!
      use Sub, only: del6
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension(nx) :: fdiff
!
      if (diffrho_hyper3 /= 0.0) then
        if (.not. ldensity_nolog) then
          call del6(f,ilnrho,fdiff,IGNOREDX=.true.)
        else
          call fatal_error('special_calc_density', &
              'not yet implemented for ldensity_nolog')
        endif
!
!        if (lfirst.and.ldt) diffus_diffrho3=diffus_diffrho3+diffrho_hyper3
!
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + diffrho_hyper3*fdiff
!
        if (headtt) print*,'special_calc_density: diffrho_hyper3=', &
            diffrho_hyper3
      endif
!
      if (swamp_diffrho > 0.0) call calc_swamp_density(df,p)
!
      if (lnrho_min > -max_real) then
        if (dt * lnrho_min_tau > 1.0) then
          if (lroot) print *, "ERROR: dt=", dt, " lnrho_min_tau=", lnrho_min_tau
          call fatal_error ('special_calc_density', "dt too large: dt * lnrho_min_tau > 1")
        endif
        fdiff = lnrho_min - p%lnrho
        where (fdiff < 0.0) fdiff = 0.0
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + lnrho_min_tau * fdiff
      endif
!
      if (lnrho_max < max_real) then
        if (dt * lnrho_min_tau > 1.0) then
          if (lroot) print *, "ERROR: dt=", dt, " lnrho_min_tau=", lnrho_min_tau
          call fatal_error ('special_calc_density', "dt too large: dt * lnrho_min_tau > 1")
        endif
        fdiff = lnrho_max - p%lnrho
        where (fdiff > 0.0) fdiff = 0.0
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + lnrho_min_tau * fdiff
      endif
!
      if (lslope_limited_special) then
        if (headtt) print*,'special_calc_density: call div_diff_flux'
        if (ibb==0)  &
          call fatal_error ('special_calc_density', "set lbb_as_aux=T in run.in")
        if (ldensity_nolog) then
          call div_diff_flux(f,irho,p,fdiff)
          if (lfirst_proc_x.and.lfrozen_bot_var_x(irho)) then
!
!  Density is frozen at boundary
!
            df(l1+1:l2,m,n,irho) = df(l1+1:l2,m,n,irho) + fdiff(2:nx)
          else
            df(l1:l2,m,n,irho) = df(l1:l2,m,n,irho) + fdiff
          endif
        else
!
! Store the diffusive flux in a special aux array to be added later to
! log-density in special_after_timestep
!
          call div_diff_flux(f,ilnrho,p,fdiff)
          if (lfirst_proc_x.and.lfrozen_bot_var_x(ilnrho)) then
            df(l1+1:l2,m,n,ilnrho) = df(l1+1:l2,m,n,ilnrho) + fdiff(2:nx)*p%rho1(2:nx)
          else
            df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + fdiff*p%rho1
          endif
        endif
      endif
!    
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy (or temperature) equation.
!
!  23-jun-08/bing: coded
!  17-feb-10/bing: added hyperdiffusion for non-equidistant grid
!
      use Sub, only: del6, del4, dot, dot2, multsv, multmv
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: hc, tmp, quenchfactor, b_abs_inv
      real, dimension(nx,3) :: hhh, tmpv
      integer :: i, j, k
!
      diffus_chi=0.; diffus_chi3=0.

      if (chi_hyper3 /= 0.0) then
        call del6(f,ilnTT,hc,IGNOREDX=.true.)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + chi_hyper3*hc
!
!  due to ignoredx chi_hyperx has [1/s]
!
        if (lfirst .and. ldt) diffus_chi3 = diffus_chi3 + chi_hyper3
      endif
!
      if (chi_hyper2 /= 0.0) then
        call del4(f,ilnTT,hc,IGNOREDX=.true.)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + chi_hyper2*hc
!
!  due to ignoredx chi_hyperx has [1/s]
!
        if (lfirst .and. ldt) diffus_chi3 = diffus_chi3 + chi_hyper2
      endif
!
      if ((K_spitzer /= 0.0) .or. (hcond1 /= 0.0) .or. (hcond2 /= 0.0)) then
        ! calculate inverse absolute value of bb
        b_abs_inv = 1./max(tini,sqrt(p%b2))
!
        ! calculate H_i = Sum_jk ( (delta_ik - 2*bunit_i*bunit_k)*bunit_j * dB_k/dj / |B| )
        do i = 1,3
          hhh(:,i) = 0.
          do j = 1,3
            tmp(:) = 0.
            do k = 1,3
              tmp(:) = tmp(:)-2.*p%bunit(:,k)*p%bij(:,k,j)
            enddo
            hhh(:,i) = hhh(:,i)+p%bunit(:,j)*(p%bij(:,i,j)+p%bunit(:,i)*tmp(:))
          enddo
        enddo
        call multsv(b_abs_inv,hhh,tmpv)
        ! calculate abs(h) limiting
        call dot2(tmpv,tmp,PRECISE_SQRT=.true.)
        ! limit the length of h
        quenchfactor = 1./max(1.,3.*tmp*dxmax)
        call multsv(quenchfactor,tmpv,hhh)
        call dot(hhh,p%glnTT,glnTT_H)
!
        ! dot Hessian matrix of lnTT with bi*bj
        call multmv(p%hlnTT,p%bunit,tmpv)
        call dot(tmpv,p%bunit,hlnTT_Bij)
!
        ! calculate Grad lnTT * bunit
        call dot(p%glnTT,p%bunit,glnTT_b)
      endif
!
      if ((K_spitzer /= 0.0) .or. (K_iso /= 0.0) .or. (hcond2 /= 0.0) .or. (hcond3 /= 0.0)) then
        call dot2(p%glnTT,glnTT2)
        if ((K_spitzer /= 0.0) .or. (K_iso /= 0.0)) then
          glnTT_abs = sqrt(glnTT2)
          glnTT_abs_inv = 1.0 / max(glnTT_abs, tini)
        endif
      endif
!
      if (K_spitzer /= 0.0) call calc_heatcond_tensor(df,p,K_spitzer,2.5)
      if (hcond1 /= 0.0) call calc_heatcond_constchi(df,p)
      if (hcond2 /= 0.0) call calc_heatcond_glnTT(df,p)
      if (hcond3 /= 0.0) call calc_heatcond_glnTT_iso(df,p)
      if (heat_cool /= 0.0) call calc_heat_cool(df,p)
      if (rho_diff_fix /= 0.0) call calc_rho_diff_fix(df,p)
      if (cool_RTV /= 0.0) call calc_heat_cool_RTV(df,p)
      if (max (nc_tau, nc_tau_rho) > 0.0) call calc_heat_cool_newton(df,p)
      if (K_iso /= 0.0) call calc_heatcond_grad(df,p)
      if (iheattype(1) /= 'nothing') call calc_artif_heating(df,p)
!
      if (swamp_chi > 0.0) call calc_swamp_temp(df,p)
!
      if (lfirst.and.ldt) then
        maxdiffus=max(maxdiffus,diffus_chi)
        maxdiffus3=max(maxdiffus3,diffus_chi3)
      endif

    endsubroutine special_calc_energy
!***********************************************************************
    subroutine special_before_boundary(f,lfinalize)
!
!   Possibility to modify the f array before the boundaries are
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      use Diagnostics, only: save_name
!
      real, dimension(mx,my,mz,mfarray) :: f
      logical, optional :: lfinalize
!
      real, save :: time_vel_l, time_vel_r
      real, save :: time_mag_l, time_mag_r, time_mag_vel_l, time_mag_vel_r
      real, dimension(:,:), allocatable, save :: vel_x_l, vel_y_l, vel_x_r, vel_y_r
      real, dimension(:,:), allocatable, save :: mag_x_l, mag_y_l, mag_x_r, mag_y_r
      real, dimension(:,:), allocatable, save :: gran_x, gran_y
      real, save :: Bz_total_flux=0.0
!
      if (present (lfinalize)) then
        if (lfinalize) then
          if (allocated (BB2_local)) deallocate (BB2_local)
          if (allocated (vel_x_l)) deallocate (vel_x_l, vel_y_l, vel_x_r, vel_y_r)
          if (allocated (mag_x_l)) deallocate (mag_x_l, mag_y_l, mag_x_r, mag_y_r)
          if (allocated (gran_x)) deallocate (gran_x, gran_y)
          call update_mag_field (0.0, "", "", A_init, time_mag_l, time_mag_r, .true.)
          return
        endif
      endif
!
      Ux_local = 0.0
      Uy_local = 0.0
!
      if (lpencil_check_at_work) return
!
      ! Prepare for footpoint quenching (lquench), flux conservation (Bz_flux),
      ! and the granulation computation (lgranulation and Bavoid).
      if ((lquench .and. lfirst_proc_z) .or. ((Bavoid > 0.) .and. lgran_proc) &
          .or. (lgranulation .and. (Bavoid > 0.) .and. lfirst_proc_z) &
          .or. (luse_mag_field .and. lfirst_proc_z)) then
        if (lfirst_proc_z .and. .not. allocated (BB2_local)) then
          allocate (BB2_local(nx,ny), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('special_before_boundary', &
              'Could not allocate memory for BB2', .true.)
        endif
        ! Compute magnetic energy BB2 for granulation driver.
        call set_BB2 (f, BB2_local, Bz_total_flux)
      endif
!
      ! External magnetic field driver
      if (luse_mag_field .and. lfirst_proc_z) then
        call update_mag_field (mag_time_offset, mag_times_dat, mag_field_dat, &
            A_init, time_mag_l, time_mag_r)
        call mag_driver (A_init, Bz_total_flux, f)
      endif
!
      ! External horizontal velocity driver
      if (luse_vel_field .and. lfirst_proc_z) then
        if (.not. allocated (vel_x_l)) then
          allocate (vel_x_l(nx,ny), vel_y_l(nx,ny), vel_x_r(nx,ny), vel_y_r(nx,ny), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('special_before_boundary', &
              'Could not allocate memory for velocity frames', .true.)
        endif
        call update_vel_field (vel_time_offset, vel_times_dat, vel_field_dat, &
            time_vel_l, time_vel_r, vel_x_l, vel_y_l, vel_x_r, vel_y_r)
        ! Quench the horizontal velocities
        if (lquench .and. lfirst_proc_z) &
            call vel_quench (BB2_local, Ux_local, Ux_local, f)
      endif
!
      ! External magnetic field horizontal LCT velocities (no quenching)
      if (luse_mag_vel_field .and. lfirst_proc_z) then
        if (.not. allocated (mag_x_l)) then
          allocate (mag_x_l(nx,ny), mag_y_l(nx,ny), mag_x_r(nx,ny), mag_y_r(nx,ny), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('special_before_boundary', &
              'Could not allocate memory for magnetic velocity frames', .true.)
        endif
        call update_vel_field (mag_time_offset, mag_times_dat, mag_vel_field_dat, &
            time_mag_vel_l, time_mag_vel_r, mag_x_l, mag_y_l, mag_x_r, mag_y_r)
      endif
!
      ! Compute photospheric granulation.
      if (lgranulation) then
        if (.not. allocated (gran_x)) then
          allocate (gran_x(nx,ny), gran_y(nx,ny), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('special_before_boundary', &
              'Could not allocate memory for gran_x/y', .true.)
        endif
        call gran_driver (f, gran_x, gran_y)
        Ux_local = Ux_local + gran_x
        Uy_local = Uy_local + gran_y
      endif
!
      if (lmassflux) call get_wind_speed_offset(f)
!
      if (ldiagnos .and. (idiag_mag_flux /= 0)) &
          call save_name (Bz_total_flux, idiag_mag_flux)
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
    subroutine special_after_timestep(f,df,dt_,llast)
!
!  Calculate an additional 'special' term on the right hand side of the
!  induction equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      use Deriv, only: der
      use Mpicomm, only: mpibcast_real,initiate_isendrcv_bdry, &
                         finalize_isendrcv_bdry
      use SharedVariables, only: get_shared_variable
      use Sub, only: cross,gij,curl_mn,step
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      logical, intent(in) :: llast
!
!      real, dimension(mx,my,mz):: rho_tmp
      real, intent(in) :: dt_
      real, dimension(nx,3) :: uu,bb,uxb,rn
      real, dimension(nx) :: va
      real, dimension(2) :: gn
      real :: uborder, border_width,nwave_disp
      real, pointer, dimension(:) :: B_ext
      real, pointer :: TTtop
      integer :: ig,ierr
!
! Flux emergence by driving an EMF at bottom boundary
! Presently only emf from constant field and constant velocity
!
      if (lset_boundary_emf) then
        call initiate_isendrcv_bdry(f,iax,iaz)
        call finalize_isendrcv_bdry(f,iax,iaz)
        do m=m1,m2
          call bc_emf_z(f,df,dt_,'top',iax)
          call bc_emf_z(f,df,dt_,'top',iay)
          call bc_emf_z(f,df,dt_,'top',iaz)
        enddo
      endif
!
! Allow putting a hot plate at z=uborder to heat the corona
!
      if (lset_hotplate_lnTT) then
          border_width=border_frac_z(2)*Lxyz(3)/2
          uborder=xyz1(3)-1.1*border_width
        call get_shared_variable('TTtop',TTtop,ierr)
        if (ierr/=0) call fatal_error('special_after_timestep:',&
             'failed to get TTtop from eos_temperature_ionization')
        do m=m1, m2; do n=n1, n2
          if (z(n) .gt. 0.9*uborder .and. z(n) .lt. 1.1*uborder) then
            f(l1:l2,m,n,ilnTT) = f(l1:l2,m,n,ilnTT) - &
                                 lnTT_hotplate_tau*(1.-&
                                 TTtop/exp(f(l1:l2,m,n,ilnTT)))* &
                                 step(z(n),uborder,0.2)*dt_
          endif     
        enddo; enddo
      endif
!
      if (lfirst_proc_z.and.lcartesian_coords) then
        if (lflux_emerg_bottom) then
          select case (flux_type)
          case ('uniform')
            do ig=0,nghost
              do m=m1, m2
                call get_shared_variable('B_ext', B_ext,ierr)
                if (ierr/=0) call fatal_error('special_after_timestep:',&
                'failed to get B_ext from magnetic')
                bb(:,1)=bb_emerg(1)+B_ext(1)
                bb(:,2)=bb_emerg(2)+B_ext(2)
                bb(:,3)=bb_emerg(3)+B_ext(3)
                if (iglobal_bx_ext/=0) bb(:,1)=bb(:,1)+f(l1:l2,m,n1-ig,iglobal_bx_ext)
                if (iglobal_by_ext/=0) bb(:,2)=bb(:,2)+f(l1:l2,m,n1-ig,iglobal_by_ext)
                if (iglobal_bz_ext/=0) bb(:,3)=bb(:,3)+f(l1:l2,m,n1-ig,iglobal_bz_ext)
                uu(:,1)=uu_emerg(1)*sin(nwave(1)*x(l1:l2)/Lxyz(1)-w_ff(1)*t)
                uu(:,2)=uu_emerg(2)*cos(nwave(1)*x(l1:l2)/Lxyz(1)-w_ff(1)*t)
                if (lemerg_profx) then
                  ! Hard coding the emerging velocity x-profile for testing
                  uu(:,3)=uu_emerg(3)*(1.0-0.29279746*abs((1+&
                          tanh((x(l1:l2)+4+0.5)/0.4))* &
                          (1-tanh((x(l1:l2)+4-0.5)/0.4))- &
                          (1-tanh((x(l1:l2)-4-0.5)/0.4))* &
                          (1+tanh((x(l1:l2)-4+0.5)/0.4))))
                else
                  uu(:,3)=uu_emerg(3)*cos(nwave(1)*x(l1:l2)/Lxyz(1)-w_ff(1)*t)
                endif
                f(l1:l2,m,n1-ig,iux) = uu(:,1)
                f(l1:l2,m,n1-ig,iuy) = uu(:,2)
                f(l1:l2,m,n1-ig,iuz) = uu(:,3)
                call cross(uu,bb,uxb)
                f(l1:l2,m,n1-ig,iax) = f(l1:l2,m,n1-ig,iax) + uxb(:,1)*dt_
                f(l1:l2,m,n1-ig,iay) = f(l1:l2,m,n1-ig,iay) + uxb(:,2)*dt_
                f(l1:l2,m,n1-ig,iaz) = f(l1:l2,m,n1-ig,iaz) + uxb(:,3)*dt_
              enddo
            enddo
          case ('noise')
            do ig=0,nghost
              do m=m1, m2
                uu(:,1)=uu_emerg(1)
                uu(:,2)=uu_emerg(2)
                if (lemerg_profx) then
!
! Hard coding the emerging velocity x-profile for testing
!
                  uu(:,3)=uu_emerg(3)*(1.0-0.29279746*abs((1+&
                          tanh((x(l1:l2)+4+0.5)/0.4))* &
                          (1-tanh((x(l1:l2)+4-0.5)/0.4))- &
                          (1-tanh((x(l1:l2)-4-0.5)/0.4))* &
                          (1+tanh((x(l1:l2)-4+0.5)/0.4))))
                else
                  uu(:,3)=uu_emerg(3)
                endif
                f(l1:l2,m,n1-ig,iuz) = uu(:,3)
!
                call random_number(gn)
                bb(:,1)=0.5*bb_emerg(1)*(1.+gn(1))
                bb(:,2)=bb_emerg(2)
                bb(:,3)=bb_emerg(3)*sin(pi*(2*x(l1:l2)/Lxyz(1)+gn(2)))
!                
                call cross(uu,bb,uxb)
                f(l1:l2,m,n1-ig,iax) = f(l1:l2,m,n1-ig,iax) + uxb(:,1)*dt_
                f(l1:l2,m,n1-ig,iay) = f(l1:l2,m,n1-ig,iay) + uxb(:,2)*dt_
                f(l1:l2,m,n1-ig,iaz) = f(l1:l2,m,n1-ig,iaz) + uxb(:,3)*dt_
              enddo
            enddo
          case default
            call fatal_error('special_after_timestep:','wrong flux_type')
          endselect
        endif
      endif
      if (any(uu_drive /= 0.0)) then
        do m=m1, m2
          do n=n1,n2
            if ((z(n) .lt. z_ff(2)) .and. (z(n) .gt. z_ff(1))) then
               uu(:,1)=uu_drive(1)*sin(nwave(2)*x(l1:l2)/Lxyz(1)+nwave(2)*z(n)/Lxyz(1)-w_ff(2)*t)
               uu(:,2)=uu_drive(2)*cos(nwave(2)*x(l1:l2)/Lxyz(1)+nwave(2)*z(n)/Lxyz(1)-w_ff(2)*t)
               if (Lxyz(2) /=  0.) then
                 if (ldispersion_acoustic) then
                   nwave_disp=w_ff(1)*nwave(2)/w_ff(2)
                   uu(:,3)=uu_drive(3)*cos(pi*nwave_disp*x(l1:l2)/Lxyz(1))*sin(w_ff(1)*t)*cos(pi*nwave_disp*y(m)/Lxyz(2))
                 else
                   uu(:,3)=uu_drive(3)*cos(pi*nwave(1)*x(l1:l2)/Lxyz(1))*sin(w_ff(1)*t)*cos(pi*nwave(1)*y(m)/Lxyz(2))
                 endif
               else
                 if (ldispersion_acoustic) then
                   nwave_disp=w_ff(1)*nwave(2)/w_ff(2)
                   uu(:,3)=uu_drive(3)*cos(pi*nwave_disp*x(l1:l2)/Lxyz(1))*sin(w_ff(1)*t)
!                   uu(:,3)=uu_drive(3)*cos(pi*nwave(1)*x(l1:l2)/Lxyz(1))* &
!                   (0.25*sin(.85714285714285714285*w_ff(1)*t)+0.5*sin(w_ff(1)*t)+0.25*sin(1.14285714285714285714*w_ff(1)*t))
                 else
                 uu(:,3)=uu_drive(3)*cos(pi*nwave(1)*x(l1:l2)/Lxyz(1))*sin(w_ff(1)*t)
                 endif
               endif
!               f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux)+uu(:,1)*w_ff(2)*dt_
!               f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy)+uu(:,2)*w_ff(2)*dt_
!               f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz)+uu(:,3)*w_ff(2)*dt_
               call random_number(rn)
               f(l1:l2,m,n,iux) = uu(:,1)+uu_emerg(1)*(2*rn(:,1)-1)
               f(l1:l2,m,n,iuy) = uu(:,2)+uu_emerg(2)*(2*rn(:,2)-1)
               f(l1:l2,m,n,iuz) = uu(:,3)+uu_emerg(3)*(2*rn(:,3)-1)
            else
              if (lconv_vel_set_to_zero .and. (z(n) .lt. z_ff(1))) then
                f(l1:l2,m,n,iux)=0.0
                f(l1:l2,m,n,iuy)=0.0
                f(l1:l2,m,n,iuz)=0.0 
              endif
            endif
          enddo
        enddo
     endif
!
!      if (.not.ldensity_nolog .and. lslope_limited_special) then
!        rho_tmp(l1:l2,m1:m2,n1:n2)=exp(f(l1:l2,m1:m2,n1:n2,ilnrho))+f(l1:l2,m1:m2,n1:n2,ispecaux)*dt_
!        f(l1:l2,m1:m2,n1:n2,ilnrho)=log(rho_tmp(l1:l2,m1:m2,n1:n2))
!      endif
!
    endsubroutine special_after_timestep
!***********************************************************************
    function get_swamp_fade_fact(height,deriv)
!
!   Get height-dependant smooth fading factor for swamp layer at top.
!
!   02-jun-11/Bourdin.KIS: coded
!
      real :: get_swamp_fade_fact
      real, intent(in) :: height
      real, intent(out), optional :: deriv
!
      real, save :: last_height = -1.0
      real, save :: last_fade_fact = 0.0
      real, save :: last_deriv = 0.0
      real :: tau, delta_inv
!
      if (last_height == height) then
        get_swamp_fade_fact = last_fade_fact
      else
        last_deriv = 0.0
        if (height <= swamp_fade_start) then
          get_swamp_fade_fact = 0.0
        elseif (height >= swamp_fade_end) then
          get_swamp_fade_fact = 1.0
        else
          ! tau is a normalized z, the transition interval is [-0.5, 0.5]:
          delta_inv = 1.0 / (swamp_fade_end-swamp_fade_start)
          tau = (height-swamp_fade_start) * delta_inv - 0.5
          if (tau <= -0.5) then
            get_swamp_fade_fact = 0.0
          elseif (tau >= 0.5) then
            get_swamp_fade_fact = 1.0
          else
            get_swamp_fade_fact = 0.5 + tau * (1.5 - 2.0*tau**2)
            last_deriv = (1.5 - 6.0*tau**2) * delta_inv
          endif
        endif
        last_height = height
        last_fade_fact = get_swamp_fade_fact
      endif
!
      if (present (deriv)) deriv = last_deriv
!
    endfunction get_swamp_fade_fact
!***********************************************************************
    function get_time_fade_fact()
!
!   Get time-dependant fading factor for all init_time based processes.
!
!   09-jun-11/Bourdin.KIS: coded
!
      use Sub, only: cubic_step
!
      real :: get_time_fade_fact
!
      real, save :: last_time = -max_real
      real, save :: last_fade_fact = 0.0
!
      if (last_time == t) then
        get_time_fade_fact = last_fade_fact
      elseif (t >= init_time_fade_start + init_time) then
        get_time_fade_fact = 1.0
        last_time = t
        last_fade_fact = get_time_fade_fact
      else
        get_time_fade_fact = cubic_step (real (t*unit_time) - init_time_fade_start, init_time, init_time)
        last_time = t
        last_fade_fact = get_time_fade_fact
      endif
!
    endfunction get_time_fade_fact
!***********************************************************************
    function get_hcond_fade_fact()
!
!   Get time-dependant heat conduction fading factor
!   for all init_time_hcond based processes.
!
!   09-jun-11/Bourdin.KIS: coded
!
      use Sub, only: cubic_step
!
      real :: get_hcond_fade_fact
!
      real, save :: last_time = -max_real
      real, save :: last_fade_fact = 0.0
!
      if (last_time == t) then
        get_hcond_fade_fact = last_fade_fact
      elseif (t >= init_time_hcond_fade_start + init_time_hcond) then
        get_hcond_fade_fact = 1.0
        last_time = t
        last_fade_fact = get_hcond_fade_fact
      else
        get_hcond_fade_fact = cubic_step (real (t*unit_time) - init_time_hcond_fade_start, init_time_hcond, init_time_hcond)
        last_time = t
        last_fade_fact = get_hcond_fade_fact
      endif
!
    endfunction get_hcond_fade_fact
!***********************************************************************
    subroutine calc_swamp_density(df,p)
!
!   Additional hight-dependant density diffusion (swamp layer at top).
!
!   02-jun-11/Bourdin.KIS: coded
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: fdiff
      real :: swamp_fade_fact
!
      swamp_fade_fact = get_swamp_fade_fact (z(n))
      if (swamp_fade_fact > 0.0) then
        if (ldensity_nolog) then
          fdiff = swamp_fade_fact * swamp_diffrho * p%del2rho
          df(l1:l2,m,n,irho) = df(l1:l2,m,n,irho) + fdiff
        else
          fdiff = swamp_fade_fact * swamp_diffrho * p%del2lnrho
          df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + fdiff
        endif
      endif
!
    endsubroutine calc_swamp_density
!***********************************************************************
    subroutine calc_swamp_temp(df,p)
!
!   Additional height-dependent temperature diffusion (swamp layer at top).
!
!   02-jun-11/Bourdin.KIS: coded
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: fdiff
      real :: swamp_fade_fact
!
      swamp_fade_fact = get_swamp_fade_fact (z(n))
      if (swamp_fade_fact > 0.0) then
        if (ltemperature) then
          if (ltemperature_nolog) then
            fdiff = swamp_fade_fact * swamp_chi * p%del2TT
            df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + fdiff
          else
            fdiff = swamp_fade_fact * swamp_chi * p%del2lnTT
            df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + fdiff
          endif
        elseif (lentropy) then
          fdiff = swamp_fade_fact * swamp_chi * p%del2ss
          df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) + fdiff
        endif
      endif
!
    endsubroutine calc_swamp_temp
!***********************************************************************
    subroutine calc_swamp_eta(df,p)
!
!   Additional hight-dependant magnetic diffusivity (swamp layer at top).
!
!   19-Dec-2011/Bourdin.KIS: coded
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real :: swamp_fade_fact, dfade_fact
!
      swamp_fade_fact = swamp_eta * get_swamp_fade_fact (z(n), dfade_fact)
      if (swamp_fade_fact > 0.0) then
        ! eta varies along z and requires the deta/dz term in the third component
        df(l1:l2,m,n,iax) = df(l1:l2,m,n,iax) + swamp_fade_fact * p%del2a(:,1)
        df(l1:l2,m,n,iay) = df(l1:l2,m,n,iay) + swamp_fade_fact * p%del2a(:,2)
        df(l1:l2,m,n,iaz) = df(l1:l2,m,n,iaz) + swamp_fade_fact * p%del2a(:,3) + swamp_eta * dfade_fact * p%diva
      endif
!
    endsubroutine calc_swamp_eta
!***********************************************************************
    subroutine update_vel_field (time_offset, times_dat, field_dat, time_l, time_r, Ux_l, Uy_l, Ux_r, Uy_r)
!
!  Check if an update of the velocity field is needed and load data from file.
!  An interpolated velocity field will be added to Ux and Uy.
!  The previous frame (l) and following frame (r) are updated.
!
!  18-jan-2011/Bourdin.KIS: coded
!
      real, intent(in) :: time_offset
      character(len=*), intent(in) :: times_dat, field_dat
      real, intent(inout) :: time_l, time_r
      real, dimension(nx,ny), intent(inout) :: Ux_l, Uy_l, Ux_r, Uy_r
!
      real :: time
      integer :: pos_l, pos_r
      logical, save :: lfirst_call=.true.
!
      time = t - time_offset
!
      if (lfirst_call) then
        ! Load previous (l) frame and store it in (r), will be shifted later
        call find_frame (time, times_dat, 'l', pos_l, time_l)
        if (pos_l == 0) then
          ! The simulation started before the first frame of the time series
          ! start from zero velocities
          Ux_r = 0.0
          Uy_r = 0.0
          time_l = -time_offset
        else
          call read_vel_field (pos_l, field_dat, Ux_r, Uy_r)
        endif
        ! Make sure that the following (r) frame will get loaded:
        time_r = time_l
        lfirst_call = .false.
      endif
!
      if (time >= time_r) then
        ! Shift data from following (r) to previous (l) frame
        Ux_l = Ux_r
        Uy_l = Uy_r
        time_l = time_r
        ! Read new following (r) frame
        call find_frame (time, times_dat, 'r', pos_r, time_r)
        call read_vel_field (pos_r, field_dat, Ux_r, Uy_r)
      endif
!
      ! Add interpolated values to global velocity field
      call add_interpolated (time, time_l, time_r, Ux_l, Ux_r, Ux_local)
      call add_interpolated (time, time_l, time_r, Uy_l, Uy_r, Uy_local)
!
    endsubroutine update_vel_field
!***********************************************************************
    subroutine update_mag_field (time_offset, times_dat, field_dat, A, time_l, time_r, lfinalize)
!
!  Check if an update of the vertical magnetic field is needed and load data.
!  An interpolated vector field will be added to Ax and Ay.
!  The previous frame (l) and following frame (r) are updated.
!
!  18-jan-2011/Bourdin.KIS: coded
!
      real, intent(in) :: time_offset
      character(len=*), intent(in) :: times_dat, field_dat
      real, dimension(nx,ny,2), intent(out) :: A
      real, intent(inout) :: time_l, time_r
      logical, optional :: lfinalize
!
      real, dimension(:,:,:,:), allocatable, save :: A_l, A_r
      logical, save :: lfirst_call=.true.
      real :: time
      integer :: pos_l, pos_r
!
      if (.not. allocated (A_l)) then
        allocate (A_l(nx,ny,1,2), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('update_mag_field', &
            'Could not allocate memory for A_l', .true.)
      endif
      if (.not. allocated (A_r) .and. luse_timedep_magnetogram) then
        allocate (A_r(nx,ny,1,2), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('update_mag_field', &
            'Could not allocate memory for A_r', .true.)
      endif
!
      if (present (lfinalize)) then
        if (lfinalize) then
          call read_mag_field (1, field_dat, A_l, .true.)
          if (allocated (A_l)) deallocate (A_l)
          if (allocated (A_r)) deallocate (A_r)
          return
        endif
      endif
!
      if (.not. luse_timedep_magnetogram) then
        ! Load first frame and use it as fixed magnetogram
        time_l = -time_offset
        time_r = huge (1.0)
        if (.not. lfirst_call) return
        call read_mag_field (1, field_dat, A_l)
        A = A_l(:,:,1,:)
        deallocate (A_l)
        lfirst_call = .false.
        return
      endif
!
      time = t - time_offset
!
      if (lfirst_call) then
        ! Load previous (l) frame and store it in (r), will be shifted later
        call find_frame (time, times_dat, 'l', pos_l, time_l)
        if (pos_l == 0) then
          ! The simulation started before the first frame of the time series
          ! start with a fixed magnetogram using the first frame
          pos_l = 1
          time_l = -time_offset
        endif
        call read_mag_field (pos_l, field_dat, A_r)
        ! Make sure that the following (r) frame will get loaded:
        time_r = time_l
        lfirst_call = .false.
      endif
!
      if (time >= time_r) then
        ! Shift data from following (r) to previous (l) frame
        A_l = A_r
        time_l = time_r
        ! Read new following (r) frame
        call find_frame (time, times_dat, 'r', pos_r, time_r)
        call read_mag_field (pos_r, field_dat, A_r)
      endif
!
      ! Store interpolated values to initial vector potential
      A = 0.0
      call add_interpolated (time, time_l, time_r, A_l(:,:,1,:), A_r(:,:,1,:), A)
!
    endsubroutine update_mag_field
!***********************************************************************
    subroutine read_vel_field (frame, filename, Ux, Uy)
!
!  Reads velocity field data from a time series for a given frame position
!  and distributes the results in the xy-plane.
!  The data is expected to be in SI units, not using F77 record markers.
!
!  07-jan-2011/Bourdin.KIS: coded
!
      use Mpicomm, only: distribute_xy
!
      integer, intent(in) :: frame
      character(len=*), intent(in) :: filename
      real, dimension(nx,ny), intent(out) :: Ux, Uy
!
      integer, parameter :: unit=12
      real, dimension(:,:), allocatable :: tmp_x, tmp_y
      integer :: rec_len
!
      if (lfirst_proc_xy) then
        allocate (tmp_x(nxgrid,nygrid), tmp_y(nxgrid,nygrid), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_vel_field', &
            'Could not allocate memory for tmp variables.', .true.)
!
        ! Read 2D vector field from file
        inquire (IOLENGTH=rec_len) 1.0d0
        rec_len = rec_len * nxgrid * nygrid
        open (unit, file=filename, form='unformatted', recl=rec_len, access='direct')
        read (unit, rec=2*(frame-1)+1) tmp_x
        read (unit, rec=2*(frame-1)+2) tmp_y
!
        call distribute_xy (Ux, tmp_x)
        call distribute_xy (Uy, tmp_y)
        deallocate (tmp_x, tmp_y)
      else
        call distribute_xy (Ux)
        call distribute_xy (Uy)
      endif
!
      ! velocity field is in [m/s] in the driver file, convert to PC units
      Ux = Ux / unit_velocity
      Uy = Uy / unit_velocity
!
    endsubroutine read_vel_field
!***********************************************************************
    subroutine read_mag_field (frame, filename, A, lfinalize)
!
!  Reads a Bz-magnetogram at the given position in a time series file.
!  The data is expected to be in Gauss units, not using F77 record markers.
!  After computing the vector field A, it will be distributed in the xy-plane.
!  This routine must be called from all processors in one ipz-plane.
!
!  'frame' specifies the frame number to be read.
!  'type' can be 'asc' or 'bin' to choose the file format.
!  'value' will be filled with the data devided by data_unit.
!
!  07-jan-2011/Bourdin.KIS: coded
!
      use Fourier, only: setup_extrapol_fact, field_extrapol_z_parallel
      use Mpicomm, only: mpisend_real, mpirecv_real
!
      integer, intent(in) :: frame
      character(len=*), intent(in) :: filename
      real, dimension(nx,ny,1,2), intent(out) :: A
      logical, optional :: lfinalize
!
      integer, parameter :: bnx=nxgrid, bny=ny/nprocx ! data in pencil shape
      integer, parameter :: enx=nygrid, eny=nx/nprocy ! transposed data in pencil shape
      integer, parameter :: unit=12, Bz_tag=356
      real, dimension(:,:,:), allocatable, save :: fact
      real, dimension(:,:), allocatable :: Bz
      integer :: py, partner, rec_len
!
      if (present (lfinalize)) then
        if (lfinalize) then
          if (allocated (fact)) deallocate (fact)
          return
        endif
      endif
!
      if (mod (nygrid, nprocxy) /= 0) call fatal_error ('read_mag_field', &
          'nygrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
      if (mod (nxgrid, nprocxy) /= 0) call fatal_error ('read_mag_field', &
          'nxgrid needs to be an integer multiple of nprocx*nprocy', lfirst_proc_xy)
!
      allocate (Bz(bnx,bny), stat=alloc_err)
      if (alloc_err > 0) call fatal_error ('read_mag_field', &
          'Could not allocate memory for Bz', .true.)
!
      if (lfirst_proc_xy) then
        ! Read Bz field from file and send to remote processors
        inquire (iolength=rec_len) 1.0d0
        rec_len = rec_len * bnx * bny
        open (unit, file=filename, form='unformatted', recl=rec_len, access='direct')
        do py = 1, nprocxy-1
          partner = py + ipz*nprocxy
          read (unit, rec=py+(frame-1)*nprocxy+1) Bz
          call mpisend_real (Bz, (/ bnx, bny /), partner, Bz_tag)
        enddo
        read (unit, rec=(frame-1)*nprocxy+1) Bz
        close (unit)
      else
        ! Receive Bz field
        call mpirecv_real (Bz, (/ bnx, bny /), ipz*nprocxy, Bz_tag)
      endif
!
      ! Convert Gauss (from file) to Tesla and then to PENCIL units
      Bz = Bz * 1.e-4 / unit_magnetic
!
      if (.not. allocated (fact)) then
        ! Setup exponential factor for bottom boundary
        allocate (fact(enx,eny,1), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('read_mag_field', &
            'Could not allocate memory for fact', .true.)
        call setup_extrapol_fact (z(n1:n1), z(n1), fact)
      endif
!
      call field_extrapol_z_parallel (Bz, A, fact)
!
      deallocate (Bz)
!
    endsubroutine read_mag_field
!***********************************************************************
    subroutine find_frame (time, filename, frame_type, frame_pos, frame_time)
!
!  Finds the position of the frame before/at (l) or after (r) the given time.
!  If a frame matches 'time', this frame is considered to be (l).
!  The time table is expected to be in SI units, not using F77 record markers.
!
!  'time' is the reference time that is being searched for.
!  'filename' is the time table file name.
!  'frame_type' indicates the desired frame type (l) or (r).
!  'frame_pos' is set to the position (record number) of the desired frame.
!  'frame_time' is set to the time of the corresponding frame.
!
!  07-jan-2011/Bourdin.KIS: coded
!
      use File_io, only: file_exists
      use Mpicomm, only: mpisend_int, mpirecv_int, mpisend_real, mpirecv_real
!
      real, intent(in) :: time
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: frame_type
      integer, intent(out) :: frame_pos
      real, intent(out) :: frame_time
!
      integer :: px, py, partner, rec_len, io_error
      real :: time_l, delta_t
      integer, parameter :: unit=12, tag_pos=314, tag_time=315
!
      if (lfirst_proc_xy) then
!
        if (.not. file_exists (filename)) then
          ! No time series => use only first frame, forever
          call warning ('solar_corona', '"'//trim (filename)//'" not found, using only first frame, forever!')
          frame_pos = 1
          frame_time = huge (0.0)
        else
          ! Read the time table from file
          inquire (iolength=rec_len) time
          open (unit, file=filename, form='unformatted', recl=rec_len, access='direct')
!
          io_error = 0
          delta_t = 0.0
          frame_pos = 1
          read (unit, rec=frame_pos) time_l
          time_l = time_l / unit_time
          if (time_l < 0.0) call fatal_error ('find_frame', trim (filename)//' first frame must be >= 0.', .true.)
          if (time < time_l) then
            ! 'time' is still before the first frame
            ! => set following (r) frame to point to the first frame
            io_error = -1
            frame_time = time_l
            time_l = 0.0
          endif
!
          do while (io_error == 0)
            frame_pos = frame_pos + 1
            read (unit, rec=frame_pos, iostat=io_error) frame_time
            if (io_error == 0) then
              frame_time = frame_time / unit_time + delta_t
              ! Test if correct time step has been reached
              if ((time >= time_l) .and. (time < frame_time)) exit
              ! If not, continue searching...
              time_l = frame_time
            else
              ! There was an error while reading, check why
              if (frame_pos <= 2) call fatal_error ('find_frame', &
                  trim (filename)//' contains less than two frames.', .true.)
              if (time_l <= 0.0) call fatal_error ('find_frame', &
                  trim (filename)//' last frame must have time > 0.', .true.)
              ! EOF reached => read from beginning
              delta_t = delta_t + time_l
              frame_pos = 0
              io_error = 0
            endif
          enddo
          close (unit)
!
          if (frame_type == 'l') then
            frame_pos = frame_pos - 1
            frame_time = time_l
          endif
        endif
!
        ! Distribute results in the xy-plane
        do px = 0, nprocx-1
          do py = 0, nprocy-1
            partner = px + py*nprocx + ipz*nprocxy
            if (partner == iproc) cycle
            call mpisend_int (frame_pos, partner, tag_pos)
            call mpisend_real (frame_time, partner, tag_time)
          enddo
        enddo
      else
        ! Receive results
        call mpirecv_int (frame_pos, ipz*nprocxy, tag_pos)
        call mpirecv_real (frame_time, ipz*nprocxy, tag_time)
      endif
!
    endsubroutine find_frame
!***********************************************************************
    subroutine vel_quench (BB2, Ux, Uy, f)
!
! Apply quenching to a given velocity field.
!
! 30-oct-2011/Bourdin.KIS: extracted from gran_driver
!
      use EquationOfState, only: gamma1, get_cp1, gamma_m1, lnrho0, cs20
      use Sub, only: cubic_step
!
      real, dimension(nx,ny), intent(in) :: BB2
      real, dimension(nx,ny), intent(inout) :: Ux, Uy
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      real, dimension(nx,ny) :: pp_tmp, beta, quench
      real :: cp1 = 1.0
      integer :: py
!
      ! Compute pressure for footpoint quenching.
      if (leos) call get_cp1(cp1)
      if (ltemperature) then
        if (ltemperature_nolog) then
          pp_tmp = gamma_m1*gamma1/cp1*f(l1:l2,m1:m2,irefz,iTT)
        else
          pp_tmp = gamma_m1*gamma1/cp1*exp(f(l1:l2,m1:m2,irefz,ilnTT))
        endif
      elseif (lentropy .and. pretend_lnTT) then
        pp_tmp = gamma_m1*gamma1/cp1*exp(f(l1:l2,m1:m2,irefz,ilnTT))
      elseif (lthermal_energy .or. lentropy) then
        call fatal_error('vel_quench', &
            'quenching not implemented for lthermal_energy or lentropy', .true.)
      else
        pp_tmp = gamma1*cs20
      endif
      if (ldensity) then
        if (ldensity_nolog) then
          pp_tmp = pp_tmp*f(l1:l2,m1:m2,irefz,irho)
        else
          pp_tmp = pp_tmp*exp(f(l1:l2,m1:m2,irefz,ilnrho))
        endif
      else
        pp_tmp = pp_tmp*exp(lnrho0)
      endif
!
      ! Compute plasma beta.
      beta = 2 * mu0 * pp_tmp / max (tini, BB2)
!
      ! Quench velocities to a fraction of the granule velocities.
      do py = 1, ny
        quench(:,py) = cubic_step (beta(:,py), q0, qw) * (1 - dq) + dq
      enddo
!
      Ux = Ux * quench
      Uy = Uy * quench
!
    endsubroutine vel_quench
!***********************************************************************
    subroutine vel_driver (Ux, Uy, f, df)
!
! Drive bottom boundary horizontal velocities with given velocity field.
!
! 22-jan-2011/Bourdin.KIS: coded
!
      use Diagnostics, only: max_mn_name
!
      real, dimension(nx,ny), intent(in) :: Ux, Uy
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
!
      real, dimension(nx) :: tmp
!
      if ((n == irefz) .and. lfirst_proc_z) then
        ! Push velocity field into the desired direction.
        df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) - tau_inv * (f(l1:l2,m,n,iux) - Ux(:,m-nghost))
        df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) - tau_inv * (f(l1:l2,m,n,iuy) - Uy(:,m-nghost))
      endif
!
      if (lfirst .and. ldt) then
        tmp(:) = tau_inv
        if (ldiagnos .and. (idiag_dtvel /= 0)) then
          itype_name(idiag_dtvel) = ilabel_max_dt
          call max_mn_name (tmp/cdts, idiag_dtvel, l_dt=.true.)
        endif
        dt1_max = max (dt1_max, tmp/cdts)
      endif
!
    endsubroutine vel_driver
!***********************************************************************
    subroutine mag_driver (A, Bz_total_flux, f)
!
! Drive bottom boundary magnetic field and apply flux conservation.
!
! 23-jan-2011/Bourdin.KIS: coded
!
      real, dimension(nx,ny,2), intent(in) :: A
      real, intent(in) :: Bz_total_flux
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      if (max (b_tau, bmdi, flux_tau) > 0.0) then
        if (dt * max (b_tau, bmdi, flux_tau) > 1.0) then
          if (lroot) print *, "ERROR: dt:", dt, " *", max (b_tau, bmdi, flux_tau), " > 1"
          call fatal_error ('solar_corona/mag_driver', &
              "dt too large: dt * max (b_tau, bmdi, flux_tau) > 1", lfirst_proc_xy)
        endif
      endif
!
      if ((Bz_flux > 0.0) .and. (Bz_total_flux > 0.0)) then
        ! Set sum(|Bz|) to the given sum of absolute vertical flux: Bz_flux
        if (flux_tau > 0.0) then
          ! Restore vertical flux with half-time flux_tau
          f(l1:l2,m1:m2,n1,iax:iaz) = f(l1:l2,m1:m2,n1,iax:iaz) * &
              (1.0 + dt*flux_tau * (Bz_flux / Bz_total_flux - 1.0))
        else
          ! Restore vertical flux immediately
          f(l1:l2,m1:m2,n1,iax:iaz) = f(l1:l2,m1:m2,n1,iax:iaz) * &
              Bz_flux / Bz_total_flux
        endif
      endif
!
      if (b_tau > 0.0) then
        ! Push vector potential back to initial setup with half-time b_tau
        f(l1:l2,m1:m2,n1,iax:iay) = f(l1:l2,m1:m2,n1,iax:iay) * (1.0 - dt*b_tau) + A * dt*b_tau
      elseif (bmdi > 0.0) then
        ! Push vector potential back to initial setup with half-time bmdi
        f(l1:l2,m1:m2,n1,iax:iay) = f(l1:l2,m1:m2,n1,iax:iay) * (1.0 - dt*bmdi) + A * dt*bmdi
      endif
!
    endsubroutine mag_driver
!***********************************************************************
    subroutine add_interpolated_2D (time, time_l, time_r, data_l, data_r, data)
!
!  Adds interpolated 2D data to a given field.
!
!  24-jan-2011/Bourdin.KIS: coded
!
      real, intent(in) :: time, time_l, time_r
      real, dimension(:,:), intent(in) :: data_l, data_r
      real, dimension(:,:), intent(inout) :: data
!
      real :: factor
!
      if (time <= time_l) then
        data = data + data_l
      elseif (time >= time_r) then
        data = data + data_r
      else
        ! Interpolate data
        factor = (time - time_l) / (time_r - time_l)
        data = data + data_l * (1.0 - factor) + data_r * factor
      endif
!
    endsubroutine add_interpolated_2D
!***********************************************************************
    subroutine add_interpolated_3D (time, time_l, time_r, data_l, data_r, data)
!
!  Adds interpolated 2D data to a given field.
!
!  24-jan-2011/Bourdin.KIS: coded
!
      real, intent(in) :: time, time_l, time_r
      real, dimension(:,:,:), intent(in) :: data_l, data_r
      real, dimension(:,:,:), intent(inout) :: data
!
      real :: factor
!
      if (time <= time_l) then
        data = data + data_l
      elseif (time >= time_r) then
        data = data + data_r
      else
        ! Interpolate data
        factor = (time - time_l) / (time_r - time_l)
        data = data + data_l * (1.0 - factor) + data_r * factor
      endif
!
    endsubroutine add_interpolated_3D
!***********************************************************************
    subroutine add_interpolated_4D (time, time_l, time_r, data_l, data_r, data)
!
!  Adds interpolated 4D data to a given field.
!
!  24-jan-2011/Bourdin.KIS: coded
!
      real, intent(in) :: time, time_l, time_r
      real, dimension(:,:,:,:), intent(in) :: data_l, data_r
      real, dimension(:,:,:,:), intent(inout) :: data
!
      real :: factor
!
      if (time <= time_l) then
        data = data + data_l
      elseif (time >= time_r) then
        data = data + data_r
      else
        ! Interpolate data
        factor = (time - time_l) / (time_r - time_l)
        data = data + data_l * (1.0 - factor) + data_r * factor
      endif
!
    endsubroutine add_interpolated_4D
!***********************************************************************
    subroutine calc_heat_cool_newton(df,p)
!
!  newton cooling dT/dt = -1/tau * (T-T0)
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: lnrho0
      use Sub, only: sine_step
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: newton, newton_rho, tmp_tau, lnTT_ref
!
      tmp_tau = 0.0
!
      ! Correction of density profile
      if (nc_tau_rho > 0.0) then
        ! Get reference density
        newton_rho = exp (lnrho_init_z(n) - p%lnrho) - 1.0
        ! nc_alt_rho is given in [Mm]
        tmp_tau = nc_tau_rho * exp (-nc_alt_rho * z(n)*unit_length*1e-6)
        ! Add correction term to density
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + newton_rho * tmp_tau
      endif
!
      ! Newton cooling of temperature profile
      if ((nc_tau /= 0.0) .and. (.not. lnc_intrin_energy_depend)) then
        if (lnc_density_depend) then
          ! Find reference temperature to actual density profile
          call find_ref_temp (p%lnrho, lnTT_ref)
          lnTT_ref = 0.5 * (lnTT_ref + lnTT_init_z(n))
        else
          ! Height dependant refernece temperaure profile
          lnTT_ref = lnTT_init_z(n)
        endif
        ! Calculate newton cooling factor to reference temperature
        newton = exp (lnTT_ref - p%lnTT) - 1.0
        if (nc_alt /= 0.0) then
          ! Calculate density-dependant inverse time scale and let cooling decay exponentially
          tmp_tau = nc_tau * exp (-nc_alt * (lnrho0 - p%lnrho))
        elseif (nc_lnrho_num_magn /= 0.0) then
          ! Calculate density-dependant inverse time scale with a smooth sine curve decay
          tmp_tau = nc_tau * sine_step (p%lnrho, lnrho0-nc_lnrho_num_magn, 0.5*nc_lnrho_trans_width, -1.0)
        endif
        ! Optional height dependant smooth cutoff
        ! Maximum height where cooling acts
        if (nc_z_max /= 0.0) &
            tmp_tau = tmp_tau * (1.0 - sine_step (z(n), nc_z_max, 0.5*nc_z_trans_width, -1.0))
        ! Minimum height where cooling acts
        if (nc_z_min /= 0.0) &
            tmp_tau = tmp_tau *  (sine_step (z(n), nc_z_min, 0.5*nc_z_min_trans_width, 1.0))
        ! Add newton cooling term to entropy
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + newton * tmp_tau
        endif
!
      ! Newton cooling that keeps the intrinsic energy constant
      if ((nc_tau /= 0.0) .and. lnc_intrin_energy_depend) then
        ! Calculate intrinsic-energy-dependant inverse time scale and let cooling decay exponentially
        tmp_tau = nc_tau * exp (-nc_alt * (lnrho0*lnTT_init_z(n1) - p%lnrho*p%lnTT))
        ! Optional height dependant smooth cutoff
        if (nc_z_max /= 0.0) &
            tmp_tau = tmp_tau * (1.0 - sine_step (z(n), nc_z_max, 0.5*nc_z_trans_width, -1.0))
        ! Calculate newton cooling factor to reference temperature
        newton = exp (lnrho_init_z(n)*lnTT_init_z(n)/p%lnrho - p%lnTT) - 1.0
        ! Calculate newton cooling factor to reference density
        newton_rho = exp (lnrho_init_z(n)*lnTT_init_z(n)/p%lnTT - p%lnrho) - 1.0
        ! Add newton cooling term to entropy and density
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + newton * tmp_tau
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + newton_rho * tmp_tau
      endif
!
      if (lfirst .and. ldt) then
        if (ldiagnos .and. (idiag_dtnewt /= 0)) then
          itype_name(idiag_dtnewt) = ilabel_max_dt
          call max_mn_name (tmp_tau/cdts, idiag_dtnewt, l_dt=.true.)
        endif
        dt1_max = max (dt1_max, tmp_tau/cdts)
      endif
!
    endsubroutine calc_heat_cool_newton
!***********************************************************************
    function interpol_tabulated (needle, haystack)
!
! Find the interpolated position of a given value in a tabulated values array.
! Bisection search algorithm with preset range guessing by previous value.
! Returns the interpolated position of the needle in the haystack.
! If needle is not inside the haystack, an extrapolated position is returned.
!
! 09-feb-2011/Bourdin.KIS: coded
!
      real :: interpol_tabulated
      real, intent(in) :: needle
      real, dimension(:), intent(in) :: haystack
!
      integer, save :: lower=1, upper=1
      integer :: mid, num, inc
!
      num = size (haystack, 1)
      if (num < 2) call fatal_error ('interpol_tabulated', "Too few tabulated values!", .true.)
      if (lower >= num) lower = num - 1
      if ((upper <= lower) .or. (upper > num)) upper = num
!
      if (haystack(lower) > haystack(upper)) then
!
!  Descending array:
!
        ! Search for lower limit, starting from last known position
        inc = 2
        do while ((lower > 1) .and. (needle > haystack(lower)))
          upper = lower
          lower = lower - inc
          if (lower < 1) lower = 1
          inc = inc * 2
        enddo
!
        ! Search for upper limit, starting from last known position
        inc = 2
        do while ((upper < num) .and. (needle < haystack(upper)))
          lower = upper
          upper = upper + inc
          if (upper > num) upper = num
          inc = inc * 2
        enddo
!
        if (needle < haystack(upper)) then
          ! Extrapolate needle value below range
          lower = num - 1
        elseif (needle > haystack(lower)) then
          ! Extrapolate needle value above range
          lower = 1
        else
          ! Interpolate needle value
          do while (lower+1 < upper)
            mid = lower + (upper - lower) / 2
            if (needle >= haystack(mid)) then
              upper = mid
            else
              lower = mid
            endif
          enddo
        endif
        upper = lower + 1
        interpol_tabulated = lower + (haystack(lower) - needle) / (haystack(lower) - haystack(upper))
!
      elseif (haystack(lower) < haystack(upper)) then
!
!  Ascending array:
!
        ! Search for lower limit, starting from last known position
        inc = 2
        do while ((lower > 1) .and. (needle < haystack(lower)))
          upper = lower
          lower = lower - inc
          if (lower < 1) lower = 1
          inc = inc * 2
        enddo
!
        ! Search for upper limit, starting from last known position
        inc = 2
        do while ((upper < num) .and. (needle > haystack(upper)))
          lower = upper
          upper = upper + inc
          if (upper > num) upper = num
          inc = inc * 2
        enddo
!
        if (needle > haystack(upper)) then
          ! Extrapolate needle value above range
          lower = num - 1
        elseif (needle < haystack(lower)) then
          ! Extrapolate needle value below range
          lower = 1
        else
          ! Interpolate needle value
          do while (lower+1 < upper)
            mid = lower + (upper - lower) / 2
            if (needle < haystack(mid)) then
              upper = mid
            else
              lower = mid
            endif
          enddo
        endif
        upper = lower + 1
        interpol_tabulated = lower + (needle - haystack(lower)) / (haystack(upper) - haystack(lower))
      else
        interpol_tabulated = -1.0
        call fatal_error ('interpol_tabulated', "Tabulated values are invalid!", .true.)
      endif
!
    endfunction interpol_tabulated
!***********************************************************************
    subroutine find_ref_temp (lnrho, lnTT_ref)
!
! Find reference temperature for a given density.
!
! 08-feb-2011/Bourdin.KIS: coded
!
      use Mpicomm, only: globalize_z, mpibcast_real
!
      real, dimension(nx), intent(in) :: lnrho
      real, dimension(nx), intent(out) :: lnTT_ref
!
      integer :: px, z_ref
      real :: pos, frac
      real, dimension(:), allocatable, save :: lnrho_global_z, lnTT_global_z
      logical, save :: lfirst_call=.true.
!
      if (lfirst_call) then
        allocate (lnrho_global_z(mzgrid), lnTT_global_z(mzgrid), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('find_ref_temp', &
            'Could not allocate memory for lnrho_/lnTT_global_z', .true.)
        call globalize_z (lnrho_init_z, lnrho_global_z)
        call globalize_z (lnTT_init_z, lnTT_global_z)
        call mpibcast_real (lnrho_global_z, mzgrid)
        call mpibcast_real (lnTT_global_z, mzgrid)
        lfirst_call = .false.
      endif
!
      do px = 1, nx
        pos = interpol_tabulated (lnrho(px), lnrho_global_z)
        z_ref = floor (pos)
        if (z_ref < 1) z_ref = 1
        if (z_ref >= mzgrid) z_ref = mzgrid - 1
        frac = pos - z_ref
        lnTT_ref(px) = lnTT_global_z(z_ref) * (1.0-frac) + lnTT_global_z(z_ref+1) * frac
      enddo
!
    endsubroutine find_ref_temp
!***********************************************************************
    subroutine calc_heatcond_tensor(df,p,Kpara,expo)
!
!    field-aligned heat conduction that is proportional to T^expo with b=B/|B|
!    L = Div (b K T^expo b*Grad(T))
!      = Div (b K T^(expo+1) b*Grad(lnT))
!    Spitzer-type coronal parameters: expo=2.5, K=9e-12 [kg*m/s^3/K^3.5]
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot, dot2, step
      use SharedVariables, only: get_shared_variable
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, intent(in) :: Kpara, expo
!
      real, dimension(nx,3) :: tmpv, gKp
      real, dimension(nx) :: cos_B_glnTT, gKp_b
      real, dimension(nx) :: chi_spitzer, chi_sat, chi_clight, fdiff
      real, pointer :: z_cutoff
      integer :: i, j
      integer :: ierr
!
      ! heatflux density vector: q = kappa * grad_T [W/m^2]
      ! thermal diffusivity: chi = gamma * kappa / (rho * cp)
      ! chi = gamma * chi_spitzer
      ! chi_spitzer = kappa / (rho * cp) ; kappa = K_spitzer * T^2.5
      chi_spitzer = Kpara * p%rho1 * p%TT**expo * p%cp1 * get_hcond_fade_fact()
!
      tmpv(:,:) = 0.
      do i = 1,3
        do j = 1,3
          tmpv(:,i) = tmpv(:,i)+p%glnTT(:,j)*p%hlnTT(:,j,i)
        enddo
      enddo
!
      gKp = (expo+1.) * p%glnTT
!
!  Limit heat condcution to a fraction of the available internal energy (Ksat)
!
      if (Ksat /= 0.) then
        chi_sat = sqrt(p%TT)*glnTT_abs_inv * p%cp1 * Ksat * 7.28e7*unit_temperature**1.5/unit_velocity**3
        where (chi_spitzer > chi_sat)
          gKp(:,1) = p%glnrho(:,1) + 1.5*p%glnTT(:,1) - tmpv(:,1)*glnTT_abs_inv
          gKp(:,2) = p%glnrho(:,2) + 1.5*p%glnTT(:,2) - tmpv(:,2)*glnTT_abs_inv
          gKp(:,3) = p%glnrho(:,3) + 1.5*p%glnTT(:,3) - tmpv(:,3)*glnTT_abs_inv
          chi_spitzer =  chi_sat
        endwhere
      endif
      if (lheatcond_cutoff) then 
        call get_shared_variable('z_cutoff',&
             z_cutoff,ierr)
        if (ierr/=0) call fatal_error('calc_heatcond_tensor:',&
             'failed to get z_cutoff from radiation_ray')
        chi_spitzer = chi_spitzer* &
          step(z(n),z_cutoff,0.2)
      endif
!
!  Limit heat conduction so that the diffusion speed
!  is smaller than a given fraction of the speed of light (Kc*c_light)
!
      if (Kc /= 0.) then
        chi_clight = Kc * c_light / max(dy_1(m),max(dz_1(n),dx_1(l1:l2)))
!
        where (chi_spitzer > chi_clight)
          chi_spitzer = chi_clight
          gKp(:,1) = p%glnrho(:,1)+p%glnTT(:,1)
          gKp(:,2) = p%glnrho(:,2)+p%glnTT(:,2)
          gKp(:,3) = p%glnrho(:,3)+p%glnTT(:,3)
        endwhere
      endif
!
      call dot(p%bunit,gKp,gKp_b)
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + &
          (glnTT_H + gKp_b*glnTT_b + hlnTT_Bij) * chi_spitzer*gamma
!
!  for timestep extension multiply with the
!  cosine between Grad(T) and bunit
!
      call dot(p%bunit,p%glnTT,cos_B_glnTT)
!
      where (glnTT_abs <= tini)
        cos_B_glnTT = 0.
      elsewhere
        cos_B_glnTT = cos_B_glnTT*glnTT_abs_inv
      endwhere
!
      if (lfirst .and. ldt) then
        fdiff = gamma*chi_spitzer * abs(cos_B_glnTT) * dxyz_2
        diffus_chi = diffus_chi+fdiff
        if (ldiagnos .and. (idiag_dtspitzer /= 0)) then
          call max_mn_name(fdiff/cdtv,idiag_dtspitzer,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_tensor
!***********************************************************************
    subroutine calc_heatcond_grad(df,p)
!
! isotropic heat conduction that is proportional to rho |Grad(T)|
! L = Div (K rho |Grad(T)| Grad(T))
! K = K_iso [m^5/s^3/K^2]
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx,3) :: tmpv
      real, dimension(nx) :: rhs, tmp, glnrho_glnTT, fdiff
      integer :: i, j
!
      call dot(p%glnrho,p%glnTT,glnrho_glnTT)
!
      tmpv(:,:) = 0.
      do i = 1,3
        do j = 1,3
          tmpv(:,i) = tmpv(:,i)+p%glnTT(:,j)*p%hlnTT(:,j,i)
        enddo
      enddo
      call dot(tmpv,p%glnTT,tmp)
!
      rhs = p%TT*(glnTT2*(p%del2lnTT+2.*glnTT2 + glnrho_glnTT)+tmp)*glnTT_abs_inv
!
      if (.not. ltemperature_nolog) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs * K_iso
      else
        df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + rhs * p%TT * K_iso
      endif
!
      if (lfirst .and. ldt) then
        fdiff = gamma*K_iso * p%TT * glnTT_abs * dxyz_2
        diffus_chi = diffus_chi+fdiff
        if (ldiagnos .and. (idiag_dtchi2 /= 0)) then
          call max_mn_name(fdiff/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_grad
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
! field-aligned heat conduction that is proportional to rho
! L = Div (b K rho b*Grad(T))
! K = chi [m^2/s] * cV [J/kg/K] = hcond1 [m^4/s^3/K]
! Define chi = K_0/rho
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: rhs, glnrho_b, fdiff
      real :: chi
!
      if (headtt) print*,'solar_corona/calc_heatcond_chiconst',hcond1
!
!  calculate Grad lnrho * bunit
!
      call dot(p%glnrho,p%bunit,glnrho_b)
!
!  calculate rhs
!
      chi = hcond1 * get_hcond_fade_fact()
      rhs = (glnTT_H + hlnTT_Bij + (glnrho_b + glnTT_b)*glnTT_b) * chi*gamma
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
!
      if (lfirst .and. ldt) then
        advec_cs2 = max(advec_cs2,chi*maxval(dxyz_2))
        fdiff = gamma*chi * dxyz_2
        diffus_chi = diffus_chi+fdiff
        if (ldiagnos .and. (idiag_dtchi2 /= 0)) then
          call max_mn_name(fdiff/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_glnTT(df,p)
!
! field-aligned heat conduction that is proportional to rho*Grad(lnT)^2
!  L = Div (b K rho Grad(lnT)^2 b*Grad(T))
!  => flux = K rho T Grad(lnT)^2
!    gflux = flux * (glnT + glnrho + Grad(Grad(lnT)^2))
!  => chi = Grad(lnT)^2
!  K = hcond2 [m^6/s^3/K]
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot, multsv
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx,3) :: tmpv, gflux
      real, dimension(nx) :: tmp, rhs, fdiff
      real :: chi
      integer :: i
!
      call multsv(glnTT2,p%glnTT + p%glnrho,gflux)
!
      do i = 1,3
        tmpv(:,i) = 2. * ( &
            p%glnTT(:,1)*p%hlnTT(:,i,1) + &
            p%glnTT(:,2)*p%hlnTT(:,i,2) + &
            p%glnTT(:,3)*p%hlnTT(:,i,3) )
      enddo
!
      call dot(gflux + tmpv,p%bunit,tmp)
!
      chi = hcond2 * get_hcond_fade_fact()
      rhs = (tmp*glnTT_b + glnTT2*(glnTT_H + hlnTT_Bij)) * chi*gamma
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
!
      if (lfirst .and. ldt) then
        fdiff = gamma*chi * glnTT2 * dxyz_2
        diffus_chi = diffus_chi+fdiff
        if (ldiagnos .and. (idiag_dtchi2 /= 0)) then
          call max_mn_name(fdiff/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_glnTT
!***********************************************************************
    subroutine calc_heatcond_glnTT_iso(df,p)
!
! isotropic heat conduction that is proportional to rho*Grad(lnT)^2
!  L = Div (K rho Grad(lnT)^2 Grad(T))
!  K = hcond3 [m^6/s^3/K]
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: gamma
      use Sub, only: dot
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx,3) :: tmpv
      real, dimension(nx) :: glnT_glnrho
      real, dimension(nx) :: rhs, tmp, fdiff
      real :: chi
      integer :: i
!
      call dot(p%glnTT,p%glnrho,glnT_glnrho)
!
      do i = 1,3
        tmpv(:,i) = &
            p%glnTT(:,1)*p%hlnTT(:,1,i) + &
            p%glnTT(:,2)*p%hlnTT(:,2,i) + &
            p%glnTT(:,3)*p%hlnTT(:,3,i)
      enddo
      call dot(p%glnTT,tmpv,tmp)
!
      chi = hcond3 * get_hcond_fade_fact()
      rhs = (2*tmp + glnTT2 * (glnTT2 + p%del2lnTT + glnT_glnrho)) * chi*gamma
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
!
      if (lfirst .and. ldt) then
        fdiff = gamma*chi * glnTT2*dxyz_2
        diffus_chi = diffus_chi+fdiff
        if (ldiagnos .and. (idiag_dtchi2 /= 0)) then
          call max_mn_name(fdiff/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_glnTT_iso
!***********************************************************************
    subroutine calc_heat_cool(df,p)
!
!  Apply external heating and cooling profile.
!
!  08-Mar-2017/PABourdin: coded
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      if (allocated (deltaT_init_z)) call calc_heat_cool_deltaT(df,p)
      if (allocated (deltaE_init_z)) call calc_heat_cool_deltaE(df,p)
      if (allocated (deltaH_part_init_z)) call calc_heat_cool_H_part(df,p)
      if (allocated (E_init_z)) call calc_heat_cool_E(df,p)
!
    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine calc_heat_cool_deltaT(df,p)
!
!  Apply external heating and cooling profile.
!
!  30-Sep-2016/PABourdin: coded
!
      use Diagnostics,     only: max_mn_name
      use Mpicomm,         only: stop_it
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: delta_lnTT, tmp
!
!     add to temperature equation
!
      if (ltemperature .and. ltemperature_nolog) then
        delta_lnTT = deltaT_init_z(n)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + delta_lnTT
      elseif (ltemperature) then
        delta_lnTT = -df(l1:l2,m,n,ilnTT)
        df(l1:l2,m,n,ilnTT) = alog (exp (df(l1:l2,m,n,ilnTT)) + deltaT_init_z(n))
        delta_lnTT = delta_lnTT + df(l1:l2,m,n,ilnTT)
      else
        if (lentropy) &
            call stop_it('solar_corona: calc_heat_cool_deltaT:lentropy=not implemented')
      endif
!
      if (lfirst .and. ldt) then
        tmp = delta_lnTT / cdts
        if (ldiagnos .and. (idiag_dtradloss /= 0)) then
          itype_name(idiag_dtradloss) = ilabel_max_dt
          call max_mn_name(tmp, idiag_dtradloss, l_dt=.true.)
        endif
        dt1_max = max(dt1_max, tmp)
      endif
!
    endsubroutine calc_heat_cool_deltaT
!***********************************************************************
    subroutine calc_heat_cool_deltaE(df,p)
!
!  Apply external heating and cooling profile.
!
!  08-Mar-2017/PABourdin: coded
!
      use Diagnostics,     only: max_mn_name
      use EquationOfState, only: gamma
      use Mpicomm,         only: stop_it
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: tmp
!
!     add to energy equation
!
!write (100+iproc,*) 'it    :', it
!write (100+iproc,*) 'pos_z :', n
!write (100+iproc,*) 'p%TT1 :', p%TT1
!write (100+iproc,*) 'p%rho1:', p%rho1
!write (100+iproc,*) 'deltaE_init_z(n):', deltaE_init_z(n)
      if (ltemperature .and. ltemperature_nolog) then
        tmp = p%rho1 * p%cp1 * gamma * deltaE_init_z(n)
!write (100+iproc,*) 'delta_T  :', tmp
! flush (100+iproc)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + tmp
        tmp = alog (1 + tmp)
      elseif (ltemperature) then
        tmp = p%TT1 * p%rho1 * p%cp1 * gamma * deltaE_init_z(n)
!write (100+iproc,*) 'delta_lnT:', tmp
! flush (100+iproc)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + tmp
      else
        if (lentropy) &
            call stop_it('solar_corona: calc_heat_cool_deltaE:lentropy=not implemented')
      endif
!
      if (lfirst .and. ldt) then
        tmp = tmp / cdts
        if (ldiagnos .and. (idiag_dtradloss /= 0)) then
          itype_name(idiag_dtradloss) = ilabel_max_dt
          call max_mn_name(tmp, idiag_dtradloss, l_dt=.true.)
        endif
        dt1_max = max(dt1_max, tmp)
      endif
!
    endsubroutine calc_heat_cool_deltaE
!***********************************************************************

subroutine calc_heat_cool_H_part(df,p)
!
!  Apply external heating and cooling profile for heating per particle.
!
      use Diagnostics,     only: max_mn_name
      use EquationOfState, only: gamma, getmu
      use Mpicomm,         only: stop_it
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: tmp,n_rho
      real :: mu 
!
!     add to energy equation
! Calculate number of particles
     !   n_rho = rho / (pc_get_parameter ('m_proton', label=quantity) * mu)
     call getmu(mu_tmp=mu) 
     n_rho= p%rho/(1.67d-27/unit_mass * mu)
      if (ltemperature .and. ltemperature_nolog) then
        tmp =heat_cool * p%rho1 * p%cp1 * gamma * deltaH_part_init_z(n) * n_rho 
        df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) + tmp
        tmp = alog (1 + tmp)
      elseif (ltemperature) then
        tmp =heat_cool *  p%TT1 * p%rho1 * p%cp1 * gamma * deltaH_part_init_z(n) * n_rho
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + tmp
      else
        if (lentropy) &
            call stop_it('solar_corona: calc_heat_cool_deltaH_part:lentropy=not implemented')
      endif
!
      if (lfirst .and. ldt) then
        tmp = tmp / cdts
        if (ldiagnos .and. (idiag_dtradloss /= 0)) then
          itype_name(idiag_dtradloss) = ilabel_max_dt
          call max_mn_name(tmp, idiag_dtradloss, l_dt=.true.)
        endif
        dt1_max = max(dt1_max, tmp)
      endif
!
    endsubroutine calc_heat_cool_H_part
!***********************************************************************
    subroutine calc_heat_cool_E(df,p)
!
!  Restore external thermal energy profile.
!
!  12-Mar-2017/PABourdin: coded
!
      use Diagnostics,     only: max_mn_name
      use EquationOfState, only: gamma
      use Mpicomm,         only: stop_it
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: tmp
!
!     change energy equation
!
! write (100+iproc+10*it,*) 'n:', n
! write (100+iproc+10*it,*) 'E:', E_init_z(n)
! write (100+iproc+10*it,*) 'P:', 1.0/(p%TT1 * p%rho1 * p%cp1 * gamma)
! write (100+iproc+10*it,*) 'ratio:', E_init_z(n) * (p%TT1 * p%rho1 * p%cp1 * gamma)
! write (100+iproc+10*it,*) 'TT1:', p%TT1
! write (100+iproc+10*it,*) 'rho1:', p%rho1
! write (100+iproc+10*it,*) 'unit_energy:', unit_energy   ! 1.0e12 ?
! write (100+iproc+10*it,*) 'cp1:', p%cp1   ! 1.0
! write (100+iproc+10*it,*) 'gamma:', gamma ! 5/3
      if (ltemperature .and. ltemperature_nolog) then
        tmp = heat_cool * E_init_z(n) * (p%TT1 * p%rho1 * p%cp1 * gamma)
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + tmp
        tmp = alog (1 + tmp)
      elseif (ltemperature) then
        tmp = alog (heat_cool * E_init_z(n) * (p%TT1 * p%rho1 * p%cp1 * gamma))
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + tmp
      else
        if (lentropy) &
            call stop_it('solar_corona: calc_heat_cool_E:lentropy=not implemented')
      endif
!
      if (lfirst .and. ldt) then
        tmp = tmp / cdts
        if (ldiagnos .and. (idiag_dtradloss /= 0)) then
          call max_mn_name(tmp,idiag_dtradloss,l_dt=.true.)
        endif
        dt1_max = max(dt1_max, tmp)
      endif
!
    endsubroutine calc_heat_cool_E
!***********************************************************************
    subroutine calc_rho_diff_fix(df,p)
!
!  Apply external density differences profile.
!
!  08-Oct-2016/PABourdin: coded
!
      use Diagnostics,     only: max_mn_name
      use Mpicomm,         only: stop_it
!
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: delta_lnrho
!
!     add to density
!
      if (ldensity .and. ldensity_nolog) then
        delta_lnrho = deltarho_init_z(n)
      elseif (ldensity) then
        delta_lnrho = alog (exp (df(l1:l2,m,n,ilnrho)) + deltarho_init_z(n)) - df(l1:l2,m,n,ilnrho)
      else
        call stop_it('solar_corona: calc_rho_diff_fix:need density to work on it')
      endif
      df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + delta_lnrho
!
    endsubroutine calc_rho_diff_fix
!***********************************************************************
    subroutine calc_heat_cool_RTV(df,p)
!
!  Electron Temperature should be used for the radiative loss
!  L = n_e * n_H * Q(T_e)
!
!  30-jan-08/bing: coded
!
      use EquationOfState, only: gamma
      use Diagnostics,     only: max_mn_name
      use Mpicomm,         only: stop_it
      use Sub,             only: cubic_step,step
      use SharedVariables, only: get_shared_variable
      use Slices_methods,  only: store_slices
!
      integer :: ierr
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: lnQ, rtv_cool, lnTT_SI, lnneni, delta_lnTT, tmp
      real :: unit_lnQ
      real, pointer :: z_cutoff
!
      unit_lnQ = 3*alog(real(unit_velocity))+ &
          5*alog(real(unit_length))+alog(real(unit_density))
      lnTT_SI = p%lnTT + alog(real(unit_temperature))
!
!     calculate ln(ne*ni) :
!          ln(ne*ni) = ln( 1.17*rho^2/(1.34*mp)^2)
!     lnneni = 2*p%lnrho + alog(1.17) - 2*alog(1.34)-2.*alog(real(m_p))
!
      lnneni = 2.*(p%lnrho+61.4412 +alog(real(unit_mass)))
!
      call get_lnQ(lnTT_SI, lnQ, delta_lnTT)
!
      rtv_cool = lnQ-unit_lnQ+lnneni-p%lnTT-p%lnrho
      rtv_cool = p%cv1*exp(rtv_cool)
!
      rtv_cool = rtv_cool*cool_RTV
!     for adjusting by setting cool_RTV in run.in
!
      select case (cool_RTV_cutoff)
      case(0)
        rtv_cool = rtv_cool*get_time_fade_fact() &
          *(1.-cubic_step(p%lnrho,-12.-alog(real(unit_density)),3.))
      case(1)
!
! Do nothing actually!
!
      case(2)
        if (lradiation) then
          call get_shared_variable('z_cutoff',&
             z_cutoff,ierr)
          if (ierr/=0) call fatal_error('calc_heat_cool_RTV:',&
             'failed to get z_cutoff from radiation_ray')
        rtv_cool = rtv_cool &
          *step(z(n),z_cutoff,0.2)
        else
        rtv_cool = rtv_cool &
          *step(z(n),1.2,0.2)
        endif
      case default
        call fatal_error('cool_RTV_cutoff:','wrong value')
      endselect
!
! slices
!
      if (ivid_rtv/=0) &
        call store_slices(rtv_cool,rtv_xy,rtv_xz,rtv_yz,rtv_xy2, &
                                   rtv_xy3,rtv_xy4,rtv_xz2,rtv_r)
!
!     add to temperature equation
!
      if (ltemperature) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)-rtv_cool
      else
        if (lentropy) &
            call stop_it('solar_corona: calc_heat_cool:lentropy=not implemented')
      endif
!
      if (lfirst .and. ldt) then
        tmp = max (rtv_cool/cdts, abs (rtv_cool/max (tini, delta_lnTT)))
        if (ldiagnos .and. idiag_dtradloss /= 0) then
          itype_name(idiag_dtradloss) = ilabel_max_dt
          call max_mn_name(tmp,idiag_dtradloss,l_dt=.true.)
        endif
        dt1_max = max(dt1_max,tmp)
      endif
!
      if (ivid_logQ/=0) &
        call store_slices(lnQ*0.43429448,logQ_xy,logQ_xz,logQ_yz,logQ_xy2, &
                                         logQ_xy3,logQ_xy4,logQ_xz2,logQ_r)
!
    endsubroutine calc_heat_cool_RTV
!***********************************************************************
    subroutine get_lnQ(lnTT,lnQ,delta_lnTT)
!
!  input: lnTT in SI units
!  output: lnP  [p]=W/s * m^3
!
      real, dimension(nx), intent(in) :: lnTT
      real, dimension(nx), intent(out) :: lnQ, delta_lnTT
!
      ! 37 points extracted from Cook et al. (1989)
      real, parameter, dimension(37) :: intlnT = (/ &
          8.74982, 8.86495, 8.98008, 9.09521, 9.21034, 9.44060, 9.67086, &
          9.90112, 10.1314, 10.2465, 10.3616, 10.5919, 10.8221, 11.0524, &
          11.2827, 11.5129, 11.7432, 11.9734, 12.2037, 12.4340, 12.6642, &
          12.8945, 13.1247, 13.3550, 13.5853, 13.8155, 14.0458, 14.2760, &
          14.5063, 14.6214, 14.7365, 14.8517, 14.9668, 15.1971, 15.4273, &
          15.6576, 69.0776 /)
      real, parameter, dimension(37) :: intlnQ = (/ &
          -93.9455, -91.1824, -88.5728, -86.1167, -83.8141, -81.6650, &
          -80.5905, -80.0532, -80.1837, -80.2067, -80.1837, -79.9765, &
          -79.6694, -79.2857, -79.0938, -79.1322, -79.4776, -79.4776, &
          -79.3471, -79.2934, -79.5159, -79.6618, -79.4776, -79.3778, &
          -79.4008, -79.5159, -79.7462, -80.1990, -80.9052, -81.3196, &
          -81.9874, -82.2023, -82.5093, -82.5477, -82.4172, -82.2637, &
          -0.66650 /)
!
      ! 16 points extracted from Cook et al. (1989)
      real, parameter, dimension(16) :: intlnT1 = (/ &
          8.98008, 9.44060, 9.90112, 10.3616, 10.8221, 11.2827, &
          11.5129, 11.8583, 12.4340, 12.8945, 13.3550, 13.8155, &
          14.2760, 14.9668, 15.8878, 18.4207 /)
      real, parameter, dimension(16) :: intlnQ1 = (/ &
          -83.9292, -81.2275, -80.0532, -80.1837, -79.6694, -79.0938, &
          -79.1322, -79.4776, -79.2934, -79.6618, -79.3778, -79.5159, &
          -80.1990, -82.5093, -82.1793, -78.6717 /)
!
      ! Four Gaussians plus one constant fitted to Cook et al. (1989)
      real, dimension(9) :: pars = (/ &
          2.12040e+00, 3.88284e-01, 2.02889e+00, 3.35665e-01, 6.34343e-01, &
          1.94052e-01, 2.54536e+00, 7.28306e-01, -2.40088e+01 /)
!
      real, dimension(nx) :: slope, ordinate
      real, dimension(nx) :: logT, logQ
      integer :: i, px, z_ref
      real :: pos, frac
!
!  select type for cooling function
!  1: 16-points logarithmic-piecewise-linear interpolation
!  2: 37-points logarithmic-piecewise-linear interpolation
!  3: four Gaussians fit
!  4: several fits
!  5: 37-points logarithmic-piecewise-linear interpolation with extrapolation
!
      select case (cool_type)
      case (1)
        lnQ = -max_real
        do i = 1, 15
          where ((lnTT >= intlnT1(i)) .and. (lnTT < intlnT1(i+1)))
            slope = (intlnQ1(i+1)-intlnQ1(i)) / (intlnT1(i+1)-intlnT1(i))
            ordinate = intlnQ1(i) - slope*intlnT1(i)
            lnQ = slope*lnTT + ordinate
          endwhere
        enddo
        delta_lnTT = max_real
!
      case (2)
        lnQ = -max_real
        do i = 1, 36
          where ((lnTT >= intlnT(i)) .and. (lnTT < intlnT(i+1)))
            slope = (intlnQ(i+1)-intlnQ(i)) / (intlnT(i+1)-intlnT(i))
            ordinate = intlnQ(i) - slope*intlnT(i)
            lnQ = slope*lnTT + ordinate
          endwhere
        enddo
        delta_lnTT = max_real
!
      case (3)
        call fatal_error('get_lnQ','this invokes exp() too often')
        logT = lnTT*alog10(exp(1.)) ! Pencil units
!
        lnQ = pars(1)*exp(-(logT-4.3)**2/pars(2)**2) &
            + pars(3)*exp(-(logT-4.9)**2/pars(4)**2) &
            + pars(5)*exp(-(logT-5.35)**2/pars(6)**2) &
            + pars(7)*exp(-(logT-5.85)**2/pars(8)**2) &
            + pars(9)
!
        lnQ = lnQ * (20.*(-tanh((logT-3.7)*10.))+21)
        lnQ = lnQ + (tanh((logT-6.9)*3.1)/2.+0.5)*3.
!
        lnQ = (lnQ+19.-32)*alog(10.)
        delta_lnTT = max_real
!
      case (4)
        logT = lnTT*alog10(exp(1.)) + alog(real(unit_temperature)) ! SI units
        where (logT <= 3.8)
          logQ = -max_real
        elsewhere ((logT > 3.8) .and. (logT <= 3.93))
          logQ = -7.155e1 + 9*logT
        elsewhere ((logT > 3.93) .and. (logT <= 4.55))
          logQ = +4.418916e+04 &
              -5.157164e+04 * logT &
              +2.397242e+04 * logT**2 &
              -5.553551e+03 * logT**3 &
              +6.413137e+02 * logT**4 &
              -2.953721e+01 * logT**5
        elsewhere ((logT > 4.55) .and. (logT <= 5.09))
          logQ = +8.536577e+02 &
              -5.697253e+02 * logT &
              +1.214799e+02 * logT**2 &
              -8.611106e+00 * logT**3
        elsewhere ((logT > 5.09) .and. (logT <= 5.63))
          logQ = +1.320434e+04 &
              -7.653183e+03 * logT &
              +1.096594e+03 * logT**2 &
              +1.241795e+02 * logT**3 &
              -4.224446e+01 * logT**4 &
              +2.717678e+00 * logT**5
        elsewhere ((logT > 5.63) .and. (logT <= 6.48))
          logQ = -2.191224e+04 &
              +1.976923e+04 * logT &
              -7.097135e+03 * logT**2 &
              +1.265907e+03 * logT**3 &
              -1.122293e+02 * logT**4 &
              +3.957364e+00 * logT**5
        elsewhere ((logT > 6.48) .and. (logT <= 6.62))
          logQ = +9.932921e+03 &
              -4.519940e+03 * logT &
              +6.830451e+02 * logT**2 &
              -3.440259e+01 * logT**3
        elsewhere (logT > 6.62)
          logQ = -3.991870e+01 + 6.169390e-01 * logT
        endwhere
        lnQ = (logQ+19.-32)*alog(10.)
        delta_lnTT = max_real
!
      case (5)
        do px = 1, nx
          pos = interpol_tabulated (lnTT(px), intlnT)
          z_ref = floor (pos)
          if (z_ref < 1) then
            lnQ(px) = -max_real
            delta_lnTT(px) = intlnT(2) - intlnT(1)
            cycle
          endif
          if (z_ref > 36) z_ref = 36
          frac = pos - z_ref
          lnQ(px) = intlnQ(z_ref) * (1.0-frac) + intlnQ(z_ref+1) * frac
          delta_lnTT(px) = intlnT(z_ref+1) - intlnT(z_ref)
        enddo
!
      case default
        call fatal_error('get_lnQ','wrong type')
      endselect
!
    endsubroutine get_lnQ
!***********************************************************************
    subroutine calc_artif_heating(df,p)
!
!  Subroutine to calculate intrisic heating.
!  Activated by setting iheattype = exp, exp2 and/or gauss
!  Maximum of 3 different possibible heating types
!  Also set the heating parameters for the chosen heating functions.
!
!  22-sept-10/Tijmen: coded
!
      use EquationOfState, only: gamma
      use Diagnostics, only: max_mn_name
      use General, only: random_number_wrapper,random_seed_wrapper, &
          normal_deviate
!
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx) :: heatinput, heat_flux
      real, dimension(nx) :: x_Mm, heat_nano, rhs
      real, dimension(nx) :: heat_event, heat_event1D
      integer, dimension(mseed) :: global_rstate
      real :: z_Mm, heat_unit
      real :: nano_sigma_t, nano_time=0., nano_start=0., nano_sigma_z
      real :: nano_flare_energy, nano_pos_x, nano_pos_z, nano_pos_y
      real :: nano_amplitude
      real, dimension(2) :: event_pos
      type (pencil_case) :: p
      integer :: i
!
      heat_unit = unit_density*unit_velocity**3/unit_length
      x_Mm = x(l1:l2)*unit_length*1e-6
      z_Mm = z(n)*unit_length*1e-6
!
      heatinput = 0.
      heat_flux = 0.
!
      do i = 1,3
        if (headtt) print*,'iheattype:',iheattype(i)
        select case(iheattype(i))
        case ('nothing')
!
        case ('exp')
          ! heat_par_exp(1) should be 530 W/m^3 (amplitude)
          ! heat_par_exp(2) should be 0.3 Mm (scale height)
!
          heatinput = heatinput + &
              heat_par_exp(1)*exp(-z_Mm/heat_par_exp(2))/heat_unit
!
          heat_flux = heat_flux +  heat_par_exp(1)*heat_par_exp(2)*1e6* &
              (1.-exp(-lz*unit_length*1e-6/heat_par_exp(2)))
!
          if (headtt) print*,'Flux of exp heating: ',heat_flux
!
        case ('exp2')
          ! A second exponential function
          ! For Sven et al. 2010 set:
          ! heat_par_exp= (1e3, 0.2 )
          ! heat_par_exp2= (1e-4, 10.)
!
          heatinput = heatinput + &
              heat_par_exp2(1)*exp(-z_Mm/heat_par_exp2(2))/heat_unit
!
          heat_flux = heat_flux + heat_par_exp2(1)*heat_par_exp2(2)*1e-6* &
              (1.-exp(-lz*unit_length*1e-6/heat_par_exp2(2)))
!
          if (headtt) print*,'Flux for exp2 heating: ', &
              heat_par_exp2(1)*heat_par_exp2(2)*1e-6* &
              (1.-exp(-lz*unit_length*1e-6/heat_par_exp2(2)))
!
        case ('gauss')
          ! heat_par_gauss(1) is Center (z in Mm)
          ! heat_par_gauss(2) is Width (sigma)
          ! heat_par_gauss(3) is the amplitude (Flux)
!
          heatinput = heatinput + &
              heat_par_gauss(3)*exp(-((z_Mm-heat_par_gauss(1))**2/ &
              (2*heat_par_gauss(2)**2)))/heat_unit
!
        case ('nanof')
          ! simulate nanoflare heating =)
          ! height dependend amplitude und duration
          ! idea: call random numbers, make condition when to flare, when
          ! flare get position
          ! then over time release a flare in shape of gaussian. Also
          ! distribution is gaussian around a point.
          ! gaussian distribution is gained via the dierfc function (double
          ! prec, inverser erf function)
          ! we draw for new flares as long as nano_time is /= 0. when done we
          ! reset nano_time to 0. again.
!
          ! Later we will implement more options for nanoflaring! =D
!
          ! SAVE GLOBAL SEED
          ! LOAD NANO_SEED
          call random_seed_wrapper(GET=global_rstate)
          call random_seed_wrapper(PUT=nano_seed)
!
          if (nano_start == 0.) then
            ! If 0 roll to see if a nanoflare occurs
            call random_number_wrapper(nano_start)
!
            if (nano_start > 0.95 ) then
              ! 5% chance for a nanoflare to occur, then get the location.
              call normal_deviate(nano_pos_z)
              nano_pos_z = nano_pos_z*lz
              call random_number_wrapper(nano_pos_y)
              call random_number_wrapper(nano_pos_x)
              nano_time = 60.
            else
              ! else no nanoflare, reset nano_start to 0. for the next
              ! timestep
              nano_start = 0.
            endif
          endif
!
          if (nano_start /= 0.) then
            ! if nano_start is not 0. then there is a nanoflare!
            ! 2nd assumption, nanoflare takes 60 seconds =)
            nano_flare_energy = 10.d17 ! joules
            nano_sigma_z = 0.5
            nano_sigma_t = 2.5
!
            nano_amplitude = nano_flare_energy/(pi/2*nano_sigma_t*nano_sigma_z*1.d6 )
!
            heat_nano = nano_amplitude*exp(-((nano_time-5.))**2/( 2*2.**2))* &
                exp(-((z_Mm-nano_pos_z)**2/ (2*0.5**2)))
            nano_time = nano_time-dt*unit_time
!
            heatinput = heatinput + heat_nano/heat_unit
!
            if (nano_time <= 0.) nano_start = 0.
          endif
!
          ! SAVE NANO_SEED
          ! RESTORE GLOBAL SEED
          call random_seed_wrapper(GET=nano_seed)
          call random_seed_wrapper(PUT=global_rstate)
!
! amp_t = exp(-((t-nano_time)**2/(2.*nano_dur**2)))
!-> no idea how to properly implement
! spread = exp(-((z_Mm-nano_pos)**2/(2.*nano_spread**2)))
!
! issue! How to put in timelike guassian for a random event starting
! at a random time?
!
        case ('event')
          ! one small point heating event (gaussian to prevent too strong gradients)
          ! one point in time, one point in space!
          if (t*unit_time > 150. .and. t*unit_time < 1000.) then
            event_pos(1) = 7.5
            event_pos(2) = 15.
            heat_event = 10.*exp(-((250.-t*unit_time))**2/(2*(20.*unit_time)**2))* &
                exp(-((x_Mm-event_pos(1))**2/ (2*0.2**2))) * &
                exp(-((z_Mm-event_pos(2))**2/ (2*0.2**2)))
            heatinput = heatinput + heat_event/heat_unit
          endif
!
        case ('event1D')
          ! one small point heating event (gaussian to prevent too strong gradients)
          ! one point in time, one point in space!
          if (t*unit_time > 300. .and. t*unit_time < 10000.) then
            if (t*unit_time > 300. .and. t*unit_time < 301.) &
                print*,'EVENTTTT!!!!!'
            event_pos(1) = 10.
            heat_event1D = 10.*exp(-((400.-t))**2/( 2*50.**2))* &
                exp(-((x_Mm-event_pos(1))**2/ (2*0.2**2)))
            heatinput = heatinput + heat_event1D/heat_unit
          endif
!
        case ('T<Tcrit')
!
          heatinput = heatinput + &
          0.5*(1.-tanh((p%TT-T_crit)/deltaT_crit))
!
        case default
          if (headtt) call fatal_error('calc_artif_heating', &
              'Please provide correct iheattype')
        endselect
      enddo
!
      if (headtt) print*,'Total flux for all types:',heat_flux
!
! Add to energy equation
!
      rhs = p%TT1*p%rho1*gamma*p%cp1 * heatinput * get_time_fade_fact()
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
!
      if (lfirst .and. ldt) then
        if (ldiagnos .and. idiag_dtnewt /= 0) then
          itype_name(idiag_dtnewt) = ilabel_max_dt
          call max_mn_name(rhs/cdts,idiag_dtnewt,l_dt=.true.)
        endif
        dt1_max = max(dt1_max,rhs/cdts)
      endif
!
    endsubroutine calc_artif_heating
!***********************************************************************
    subroutine set_gran_params()
!
      integer :: alloc_err_sum
!
! Every granule has 6 values associated with it: data(1-6).
! These contain,  x-position, y-position,
!    current amplitude, amplitude at t=t_0, t_0, and life_time.
!
! Gives intergranular area / (granular+intergranular area)
      ig = 0.3
!
! Gives average radius of granule + intergranular lane
! (no smaller than 6 grid points across)
!     here radius of granules is 0.8 Mm or bigger (3 times dx)
!
      if (unit_system == 'SI') then
        granr = max(real(0.8*1.e6/unit_length),3*dx,3*dy)
      elseif  (unit_system == 'cgs') then
        granr = max(real(0.8*1.e8/unit_length),3*dx,3*dy)
      endif
!
! Fractional difference in granule power
      pd = 0.15
!
! Gives exponential power of evolvement. Higher power faster growth/decay.
      pow = 2
!
! Fractional distance, where after emergence, no granule can emerge
! whithin this radius.(In order to not 'overproduce' granules).
! This limit is unset when granule reaches t_0.
      avoid = 0.8
!
! Lifetime of granule
! Now also resolution dependent(5 min for granular scale)
!
      life_t = (60.*5./unit_time)
      !*(granr/(0.8*1e8/u_l))**2
      !  removed since the life time was about 20 min !
!
      dxdy2 = dx**2+dy**2
!
! Typical central velocity of granule(1.5e5 cm/s=1.5km/s)
! Has to be multiplied by the smallest distance, since velocity=ampl/dist
! should now also be dependant on smallest granluar scale resolvable.
!
      if (unit_system == 'SI') then
        ampl = sqrt(dxdy2)/granr*0.28e4/unit_velocity
      elseif (unit_system == 'cgs') then
        ampl = sqrt(dxdy2)/granr*0.28e6/unit_velocity
      endif
!
! fraction of current amplitude to maximum amplitude to the beginning
! and when the granule disapears
      thresh = 0.78
!
      xrange = min(nint(1.5*granr*(1+ig)/dx),nint(nxgrid/2.0)-1)
      yrange = min(nint(1.5*granr*(1+ig)/dy),nint(nygrid/2.0)-1)
!
      if (lfirst_proc_xy) then
        print*,'| solar_corona: settings for granules'
        print*,'-----------------------------------'
        print*,'| radius [Mm]:',granr*unit_length*1e-6
        print*,'| lifetime [min]',life_t*unit_time/60.
        print*,'| update interval [s]',dt_gran*unit_time
        print*,'-----------------------------------'
      endif
!
! Don't reset if RELOAD is used
      if (.not. lreloading) then
!
        points_rstate(:) = 0
!
        isnap = ceiling (t/dsnap)
        tsnap_uu = (isnap+1) * dsnap
!
      endif
!
      alloc_err = 0
      if (.not. allocated (Ux)) allocate(Ux(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = abs(alloc_err)
      if (.not. allocated (Uy)) allocate(Uy(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (.not. allocated (Ux_ext)) allocate(Ux_ext(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (.not. allocated (Uy_ext)) allocate(Uy_ext(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (.not. allocated(w))  allocate (w(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (.not. allocated(vx))  allocate (vx(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (.not. allocated(vy))  allocate (vy(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (.not. allocated(avoid_gran)) &
          allocate(avoid_gran(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (alloc_err_sum > 0) call fatal_error ('set_gran_params', &
          'Could not allocate memory for the driver', .true.)
!
    endsubroutine set_gran_params
!***********************************************************************
    subroutine gran_driver(f, gran_x, gran_y)
!
! This routine replaces the external computing of granular velocity
! pattern initially written by B. Gudiksen (2004)
!
! It is invoked by setting lgranulation=T in run.in
! additional parameters are
! 'Bavoid': the magn. field strenght in Tesla at which no granule is allowed
! 'vorticity_factor': the strength by which the vorticity is enhanced
!
!  11-may-10/bing: coded
!
      use General, only: random_seed_wrapper
      use Mpicomm, only: mpisend_real, mpirecv_real, collect_xy, distribute_xy
      use Sub, only: cubic_step
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(nx,ny), intent(out) :: gran_x, gran_y
!
      real, dimension(:,:), allocatable :: buffer
      integer :: level, partner
      integer, dimension(mseed) :: global_rstate
      real, save :: next_time = 0.0
      integer, parameter :: tag_Ux=323, tag_Uy=324
!
      if (lwrite_driver .and. (nzgrid == 1)) then
        ! Stabilize 2D-runs
        f(:,:,:,iuz) = 0.0
        f(:,:,:,ilnTT) = lnTT_init_z(irefz)
        f(:,:,:,ilnrho) = lnrho_init_z(irefz)
      endif
!
! Update velocity field only every dt_gran after the first iteration
      if ((t < next_time) .and. .not.(lfirst .and. (it == 1))) return
      next_time = t + dt_gran
!
! Save global random number seed, will be restored after granulation
! is done
      call random_seed_wrapper(GET=global_rstate)
      call random_seed_wrapper(PUT=points_rstate)
!
! Collect external horizontal velocity field, if necessary.
      if (lfirst_proc_z .and. (luse_vel_field .or. luse_mag_vel_field)) then
        call collect_xy (Ux_local, Ux_ext)
        call collect_xy (Uy_local, Uy_ext)
        if (lgran_parallel) then
! In the parallel case, the root rank sends to the granulation processors.
          if (lroot) then
            do level = 1, nglevel
              partner = nprocxy + level - 1
              call mpisend_real (Ux_ext, (/nxgrid,nygrid/), partner, tag_Ux)
              call mpisend_real (Uy_ext, (/nxgrid,nygrid/), partner, tag_Uy)
            enddo
          elseif (lgran_proc) then
            call mpirecv_real (Ux_ext, (/nxgrid,nygrid/), 0, tag_Ux)
            call mpirecv_real (Uy_ext, (/nxgrid,nygrid/), 0, tag_Uy)
          endif
        endif
      endif
!
! Compute granular velocities.
      if (lgran_proc) call multi_gran_levels()
!
! In the parallel case, the root rank sums up the levels.
      if (lgran_parallel) then
        if (lroot) then
          allocate (buffer(nxgrid,nygrid), stat=alloc_err)
          if (alloc_err > 0) call fatal_error ('gran_driver', &
              'could not allocate memory for buffer', .true.)
          Ux = 0.0
          Uy = 0.0
          do level = 1, nglevel
            partner = nprocxy + level - 1
            call mpirecv_real (buffer, (/nxgrid,nygrid/), partner, tag_Ux)
            Ux = Ux + buffer
            call mpirecv_real (buffer, (/nxgrid,nygrid/), partner, tag_Uy)
            Uy = Uy + buffer
          enddo
          deallocate (buffer)
        elseif (lgran_proc) then
          call mpisend_real (Ux, (/nxgrid,nygrid/), 0, tag_Ux)
          call mpisend_real (Uy, (/nxgrid,nygrid/), 0, tag_Uy)
        endif
      endif
!
! Increase vorticity and normalize to given vrms.
      if (lroot) call enhance_vorticity()
!
! Distribute the results in the xy-plane.
      if (lfirst_proc_z) then
        call distribute_xy (gran_x, Ux)
        call distribute_xy (gran_y, Uy)
      endif
!
! Quench the horizontal velocities.
      if (lquench .and. lfirst_proc_z) &
          call vel_quench (BB2_local, gran_x, gran_y, f)
!
! Restore global seed and save seed list of the granulation
      call random_seed_wrapper(GET=points_rstate)
      call random_seed_wrapper(PUT=global_rstate)
!
      if (lwrite_driver .and. lroot) call write_driver (real(t), Ux, Uy)
!
    endsubroutine gran_driver
!***********************************************************************
    subroutine multi_gran_levels
!
      integer :: level
      real, parameter :: ldif=2.0
      real, dimension(max_gran_levels), save :: ampl_par, granr_par, life_t_par
      integer, dimension(max_gran_levels), save :: xrange_par, yrange_par
      logical, save :: first_call=.true.
!
      if (first_call) then
        do level = 1, nglevel
          nullify (gran_list(level)%first)
          granr_par(level) = granr * ldif**(level-1)
          ampl_par(level) = ampl / ldif**(level-1)
          life_t_par(level) = life_t * ldif**((level-1)**2)
          xrange_par(level) = min (nint (xrange * ldif**(level-1)), nint (nxgrid/2.-1))
          yrange_par(level) = min (nint (yrange * ldif**(level-1)), nint (nygrid/2.-1))
        enddo
        first_call = .false.
      endif
!
      ! Initialize velocity field
      Ux = 0.0
      Uy = 0.0
!
      ! Compute granulation levels
      do level = 1, nglevel
        ! Only let those processors work that are needed
        if (lgran_proc .and. (lroot .or. (iproc == (nprocxy-1+level)))) then
          ! Setup parameters
          ampl = ampl_par(level)
          granr = granr_par(level)
          life_t = life_t_par(level)
          xrange = xrange_par(level)
          yrange = yrange_par(level)
          ! Restore list of granules
          first => gran_list(level)%first
!
          ! Compute one granulation level
          call compute_gran_level (level)
!
          ! Store list of granules
          gran_list(level)%first => first
        endif
      enddo
!
    endsubroutine multi_gran_levels
!***********************************************************************
    subroutine compute_gran_level(level)
!
      use File_io, only: file_exists
!
      integer, intent(in) :: level
      logical :: lstop=.false.
!
      ! Initialize granule weight and velocities arrays
      w(:,:) = 0.0
      vx(:,:) = 0.0
      vy(:,:) = 0.0
      ! Initialize avoid_gran array to avoid granules at occupied places
      call fill_avoid_gran
!
      if (.not. associated(first)) then
        ! List is empty => try to read an old snapshot file
        call read_points (level)
        if (.not. associated(first)) then
          ! No snapshot available => initialize driver and draw granules
          call init_gran_driver
          call write_points (level, 0)
          call write_points (level)
        endif
      else
        ! List is filled => update amplitudes and remove dead granules
        call update_points
        ! Draw old granules to the granulation velocity field
        current => first
        do while (associated (current))
          call draw_update
          current => current%next
        enddo
        ! Fill free space with new granules and draw them
        do while (minval(avoid_gran) == 0)
          call add_point
          call find_free_place
          call draw_update
        enddo
      endif
!
      Ux = Ux + vx
      Uy = Uy + vy
!
      ! Evolve granule position by external velocity field, if necessary
      if (luse_vel_field .or. luse_mag_vel_field) call evolve_granules()
!
      ! Save regular granule snapshots
      if (t >= tsnap_uu) then
        call write_points (level, isnap)
        if (level == nglevel) then
          tsnap_uu = tsnap_uu + dsnap
          isnap  = isnap + 1
        endif
      endif
!
      ! On exit, save final granule snapshot
      if (itsub == 3) lstop = file_exists('STOP')
      if (lstop .or. (t >= tmax) .or. (it == nt) .or. (dt < dtmin) .or. &
          (mod(it,isave) == 0)) call write_points (level)
!
    endsubroutine compute_gran_level
!***********************************************************************
    subroutine write_driver (time, Ux, Uy)
!
! Write granulation driver velocity field to a file.
!
! 31-jan-2011/Bourdin.KIS: coded
!
      real, intent(in) :: time
      real, dimension(nxgrid,nygrid), intent(in) :: Ux, Uy
!
      character(len=*), parameter :: gran_times_dat = 'driver/gran_times.dat'
      character(len=*), parameter :: gran_field_dat = 'driver/gran_field.dat'
      integer, save :: gran_frame=0, unit=37
      integer :: rec_len
!
      inquire (iolength=rec_len) time
      open (unit, file=gran_times_dat, access="direct", recl=rec_len)
      write (unit, rec=gran_frame+1) time * unit_time
      close (unit)
!
      rec_len = rec_len * nxgrid*nygrid
      open (unit, file=gran_field_dat, access="direct", recl=rec_len)
      write (unit, rec=gran_frame*2+1) Ux * unit_velocity
      write (unit, rec=gran_frame*2+2) Uy * unit_velocity
      close (unit)
!
      gran_frame = gran_frame + 1
!
    endsubroutine write_driver
!***********************************************************************
    subroutine read_points(level)
!
! Read in points from file if existing.
!
! 12-aug-10/bing: coded
!
      integer, intent(in) :: level
!
      integer :: io_error, pos, len
      logical :: ex
      character(len=fnlen) :: filename
      integer, parameter :: unit=10
      real, dimension(6) :: buffer
!
      write (filename,'("driver/pts_",I1.1,".dat")') level
!
      inquire (file=filename, exist=ex)
!
      if (ex) then
        inquire (iolength=len) first%amp
        write (*,'(A,A)') 'Reading file: ', filename
        open (unit, file=filename, access="direct", recl=len*6)
!
        pos = 1
        read (unit, iostat=io_error, rec=pos) buffer
        do while (io_error == 0)
          call add_point
          current%pos_x = buffer(1)
          current%pos_y = buffer(2)
          current%amp = buffer(3)
          current%amp_max = buffer(4)
          current%t_amp_max = buffer(5)
          current%t_life = buffer(6)
          call draw_update
          pos = pos + 1
          read (unit, iostat=io_error, rec=pos) buffer
        enddo
!
        close (unit)
        write (*,'(A,I6,A,I1.1)') '=> ', pos-1, ' granules in level ', level
!
        if (level == nglevel) then
          write (filename,'("driver/seed_",I1.1,".dat")') level
          inquire (file=filename, exist=ex)
          if (.not. ex) call fatal_error ('read_points', &
              'cant find seed list for granules', .true.)
          open (unit, file=filename, access="direct", recl=len*mseed)
          read (unit, rec=1) points_rstate
          close (unit)
        endif
!
      endif
!
    endsubroutine read_points
!***********************************************************************
    subroutine write_points(level,issnap)
!
! Writes positions and amplitudes of granules to files. Also
! seed list are stored to be able to reproduce the run.
!
! 12-aug-10/bing: coded
!
      use General, only: itoa
!
      integer, intent(in) :: level
      integer, optional, intent(in) :: issnap
!
      integer :: pos, len
      character(len=fnlen) :: filename
      integer, parameter :: unit=10
      real, dimension(6) :: buffer
!
      inquire (iolength=len) first%amp
!
      if (present (issnap)) then
        write (filename,'("driver/pts_",I1.1,"_",A)') level, itoa (issnap)
      else
        write (filename,'("driver/pts_",I1.1,".dat")') level
      endif
!
      open (unit, file=filename, status="replace", access="direct", recl=len*6)
!
      pos = 1
      current => first
      do while (associated (current))
        buffer(1) = current%pos_x
        buffer(2) = current%pos_y
        buffer(3) = current%amp
        buffer(4) = current%amp_max
        buffer(5) = current%t_amp_max
        buffer(6) = current%t_life
        write (unit,rec=pos) buffer
        pos = pos + 1
        current => current%next
      enddo
!
      close (unit)
!
! Save seed list for each level. Is needed if levels are spread over 3 procs.
!
      if (present (issnap)) then
        write (filename,'("driver/seed_",I1.1,"_",A)') level, itoa (issnap)
      else
        write (filename,'("driver/seed_",I1.1,".dat")') level
      endif
!
      open (unit, file=filename, status="replace", access="direct", recl=len*mseed)
      write (unit,rec=1) points_rstate
      close (unit)
!
    endsubroutine write_points
!***********************************************************************
    subroutine add_point
!
! Add an entry to the list.
!
! 21-jan-2011/Bourdin.KIS: coded
!
      type (point), pointer, save :: new => null()
!
      allocate (new)
      nullify (new%next)
      nullify (new%previous)
!
      if (associated (first)) then
        ! Insert new entry before the first
        new%next => first
        first%previous => new
      endif
      first => new
      current => new
!
    endsubroutine add_point
!***********************************************************************
    subroutine del_point
!
! Remove an entry from the list.
!
! 21-jan-2011/Bourdin.KIS: coded
!
      type (point), pointer, save :: old => null()
!
      if (.not. associated (current)) return
!
      if (.not. associated (current%previous) .and. .not. associated (current%next)) then
        ! Current entry is the only entry
        deallocate (current)
        nullify (current)
        nullify (first)
      elseif (.not. associated (current%previous)) then
        ! Current entry is pointing to the first entry
        first => current%next
        nullify (first%previous)
        deallocate (current)
        current => first
      elseif (.not. associated (current%next)) then
        ! Current entry is pointing to the last entry
        old => current
        current => current%previous
        nullify (current%next)
        deallocate (old)
      else
        ! Current entry is between first and last entry
        current%next%previous => current%previous
        current%previous%next => current%next
        old => current
        current => current%next
        deallocate (old)
      endif
!
    endsubroutine del_point
!***********************************************************************
    subroutine init_gran_driver
!
! If no existing files are found initialize points.
! The lifetimes are randomly distribute around starting time.
!
! 12-aug-10/bing: coded
!
      use General, only: random_number_wrapper
!
      real :: rand
!
      do while (minval (avoid_gran) == 0)
!
        call add_point
        call find_free_place
!
! Set randomly some points t0 to the past so they already decay
!
        call random_number_wrapper(rand)
        current%t_amp_max = t+(rand*2-1)*current%t_life* &
            (-alog(thresh*ampl/current%amp_max))**(1./pow)
!
        current%amp = current%amp_max*exp(-((t-current%t_amp_max)/current%t_life)**pow)
!
! Update arrays with new data
!
        call draw_update
!
      enddo
!
    endsubroutine init_gran_driver
!***********************************************************************
    subroutine helmholtz(frx_r,fry_r)
!
! Extracts the rotational part of a 2d vector field
! to increase vorticity of the velocity field.
!
! 12-aug-10/bing: coded
!
      use Fourier, only: fourier_transform_other, kx_fft, ky_fft
!
      real, dimension(nxgrid,nygrid) :: kx, ky, k2, filter
      real, dimension(nxgrid,nygrid) :: fvx_r, fvy_r, fvx_i, fvy_i
      real, dimension(nxgrid,nygrid) :: frx_r, fry_r, frx_i, fry_i
      real, dimension(nxgrid,nygrid) :: fdx_r, fdy_r, fdx_i, fdy_i
      real :: k20
!
      fvx_r = vx
      fvx_i = 0.
      call fourier_transform_other(fvx_r,fvx_i)
!
      fvy_r = vy
      fvy_i = 0.
      call fourier_transform_other(fvy_r,fvy_i)
!
! Reference frequency is half the Nyquist frequency.
! This expression is only correct for grids with dx=dy.
      k20 = (kx_nyq/2.)**2.
!
      kx = spread(kx_fft,2,nygrid)
      ky = spread(ky_fft,1,nxgrid)
!
      k2 = max (kx**2 + ky**2, tini)
!
      frx_r = +ky*(ky*fvx_r - kx*fvy_r)/k2
      frx_i = +ky*(ky*fvx_i - kx*fvy_i)/k2
!
      fry_r = -kx*(ky*fvx_r - kx*fvy_r)/k2
      fry_i = -kx*(ky*fvx_i - kx*fvy_i)/k2
!
      fdx_r = +kx*(kx*fvx_r + ky*fvy_r)/k2
      fdx_i = +kx*(kx*fvx_i + ky*fvy_i)/k2
!
      fdy_r = +ky*(kx*fvx_r + ky*fvy_r)/k2
      fdy_i = +ky*(kx*fvx_i + ky*fvy_i)/k2
!
! Filter out large wave numbers.
      filter = exp(-(k2/k20)**2)
!
      frx_r = frx_r*filter
      frx_i = frx_i*filter
!
      fry_r = fry_r*filter
      fry_i = fry_i*filter
!
      fdx_r = fdx_r*filter
      fdx_i = fdx_i*filter
!
      fdy_r = fdy_r*filter
      fdy_i = fdy_i*filter
!
      call fourier_transform_other(fdx_r,fdx_i,linv=.true.)
      vx = fdx_r
      call fourier_transform_other(fdy_r,fdy_i,linv=.true.)
      vy = fdy_r
!
      call fourier_transform_other(frx_r,frx_i,linv=.true.)
      call fourier_transform_other(fry_r,fry_i,linv=.true.)
!
    endsubroutine helmholtz
!***********************************************************************
    subroutine draw_update
!
! Using a point from the list to update the velocity field.
!
! 12-aug-10/bing: coded
!
      real :: xdist, ydist, dist2, dist, wtmp, vv
      integer :: i, ii, j, jj
      real :: dist0, tmp
!
! Update weight and velocity for new granule
!
      do jj = int(current%pos_y)-yrange,int(current%pos_y)+yrange
        j = 1+mod(jj-1+nygrid,nygrid)
        do ii = int(current%pos_x)-xrange,int(current%pos_x)+xrange
          i = 1+mod(ii-1+nxgrid,nxgrid)
          xdist = dx*(ii-current%pos_x)
          ydist = dy*(jj-current%pos_y)
          dist2 = max(xdist**2+ydist**2,dxdy2)
          dist = sqrt(dist2)
!
          if (dist < avoid*granr .and. t < current%t_amp_max) avoid_gran(i,j) = 1
!
          wtmp = current%amp/dist
!
          dist0 = 0.53*granr
          tmp = (dist/dist0)**2
!
          vv = exp(1.)*current%amp*tmp*exp(-tmp)
!
          if (wtmp > w(i,j)*(1-ig)) then
            if (wtmp > w(i,j)*(1+ig)) then
              ! granular area
              vx(i,j) = vv*xdist/dist
              vy(i,j) = vv*ydist/dist
              w(i,j) = wtmp
            else
              ! intergranular area
              vx(i,j) = vx(i,j)+vv*xdist/dist
              vy(i,j) = vy(i,j)+vv*ydist/dist
              w(i,j) = max(w(i,j),wtmp)
            endif
          endif
          if (w(i,j) > ampl/(granr*(1+ig))) avoid_gran(i,j) = 1
        enddo
      enddo
!
    endsubroutine draw_update
!***********************************************************************
    subroutine find_free_place
!
! Find the position of a new granule midpoint.
!
! 12-aug-10/bing: coded
! 21-jan-2011/Bourdin.KIS: reduced runtime complexity from O(N^2) to O(N).
!
      use General, only: random_number_wrapper
!
      integer :: pos_x, pos_y, find_y, find_x, free_x, free_y
      integer :: count_x, count_y, num_free_y
      integer, dimension(nygrid) :: num_free_x
      real :: rand
!
      free_x = -1
      free_y = -1
!
      num_free_y = 0
      do pos_y = 1, nygrid
        num_free_x(pos_y) = nxgrid - sum (avoid_gran(:,pos_y))
        if (num_free_x(pos_y) > 0) num_free_y = num_free_y + 1
      enddo
      if (num_free_y == 0) return
!
! Find possible location for a new granule midpoint.
!
      call random_number_wrapper (rand)
      find_y = int (rand * num_free_y) + 1
!
      count_x = 0
      count_y = 0
      do pos_y = 1, nygrid
        if (num_free_x(pos_y) > 0) then
          count_y = count_y + 1
          if (count_y == find_y) then
!
            call random_number_wrapper (rand)
            find_x = int (rand * num_free_x(pos_y)) + 1
!
            do pos_x = 1, nxgrid
              if (avoid_gran(pos_x,pos_y) == 0) then
                count_x = count_x + 1
                if (count_x == find_x) then
                  free_x = pos_x
                  free_y = pos_y
                  exit
                endif
              endif
            enddo
!
            exit
          endif
        endif
      enddo
!
! Initialize granule properties
!
      current%pos_x = free_x
      current%pos_y = free_y
!
      call random_number_wrapper(rand)
      current%amp_max = ampl*(1+(2*rand-1)*pd)
!
      call random_number_wrapper(rand)
      current%t_life = life_t*(1+(2*rand-1)/10.)
!
      current%t_amp_max = t+current%t_life* &
          (-alog(thresh*ampl/current%amp_max))**(1./pow)
!
      current%amp = current%amp_max* &
          exp(-((t-current%t_amp_max)/current%t_life)**pow)
!
    endsubroutine find_free_place
!***********************************************************************
    subroutine update_points
!
! Update the amplitude/weight of a point.
!
! 12-aug-10/bing: coded
!
      current => first
      do while (associated (current))
!
! update amplitude
        current%amp = current%amp_max* &
            exp(-((t-current%t_amp_max)/current%t_life)**pow)
!
! remove point if amplitude is less than threshold
        if (current%amp/ampl < thresh) then
          call del_point
        else
          current => current%next
        endif
!
      enddo
!
    endsubroutine update_points
!***********************************************************************
    subroutine free_points
!
! Free list of points.
!
! 14-aug-2011/Bourdin.KIS: coded
!
      integer :: level
!
      do level = 1, nglevel
        ! Only let those processors work that are needed
        if (lgran_proc .and. (lroot .or. (iproc == (nprocxy-1+level)))) then
          ! Restore list of granules
          first => gran_list(level)%first
          ! Free list of granules
          current => first
          do while (associated (current))
            call del_point
          enddo
          nullify (gran_list(level)%first)
        endif
      enddo
!
    endsubroutine free_points
!***********************************************************************
    subroutine set_BB2(f,BB2_local,Bz_total_flux)
!
      use Mpicomm, only: collect_xy, sum_xy, mpisend_real, mpirecv_real
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,ny), intent(out) :: BB2_local
      real, intent(out) :: Bz_total_flux
!
      real, dimension(:,:), allocatable :: fac, bbx, bby, bbz
      integer :: partner, level
      integer, parameter :: BB2_tag=555
!
      if (.not. allocated (BB2) .and. (lroot .or. lgran_proc)) then
        allocate (BB2(nxgrid,nygrid), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('set_BB2', &
            'Could not allocate memory for BB2', .true.)
      endif
!
      if (lfirst_proc_z) then
!
        allocate (fac(nx,ny), bbx(nx,ny), bby(nx,ny), bbz(nx,ny), stat=alloc_err)
        if (alloc_err > 0) call fatal_error ('set_BB2', &
            'Could not allocate memory for fac and bbx/bby/bbz', .true.)
!
        if (lbfield) then
          bbx = f(l1:l2,m1:m2,irefz,ibx)
          bby = f(l1:l2,m1:m2,irefz,iby)
          bbz = f(l1:l2,m1:m2,irefz,ibz)
        else
!
! compute B = curl(A) for irefz layer
!
          bbx = 0
          bby = 0
          bbz = 0
!
          if (nxgrid /= 1) then
            fac = (1./60)*spread(dx_1(l1:l2),2,ny)
            bby = bby - fac * &
                ( 45.0*(f(l1+1:l2+1,m1:m2,irefz,iaz)-f(l1-1:l2-1,m1:m2,irefz,iaz)) &
                -  9.0*(f(l1+2:l2+2,m1:m2,irefz,iaz)-f(l1-2:l2-2,m1:m2,irefz,iaz)) &
                +      (f(l1+3:l2+3,m1:m2,irefz,iaz)-f(l1-3:l2-3,m1:m2,irefz,iaz)) )
            bbz = bbz + fac * &
                ( 45.0*(f(l1+1:l2+1,m1:m2,irefz,iay)-f(l1-1:l2-1,m1:m2,irefz,iay)) &
                -  9.0*(f(l1+2:l2+2,m1:m2,irefz,iay)-f(l1-2:l2-2,m1:m2,irefz,iay)) &
                +      (f(l1+3:l2+3,m1:m2,irefz,iay)-f(l1-3:l2-3,m1:m2,irefz,iay)) )
          else
            if (ip <= 5) print*, 'set_BB2: Degenerate case in x-direction'
          endif
!
          if (nygrid /= 1) then
            fac = (1./60)*spread(dy_1(m1:m2),1,nx)
            bbx = bbx + fac * &
                ( 45.0*(f(l1:l2,m1+1:m2+1,irefz,iaz)-f(l1:l2,m1-1:m2-1,irefz,iaz)) &
                -  9.0*(f(l1:l2,m1+2:m2+2,irefz,iaz)-f(l1:l2,m1-2:m2-2,irefz,iaz)) &
                +      (f(l1:l2,m1+3:m2+3,irefz,iaz)-f(l1:l2,m1-3:m2-3,irefz,iaz)) )
            bbz = bbz - fac * &
                ( 45.0*(f(l1:l2,m1+1:m2+1,irefz,iax)-f(l1:l2,m1-1:m2-1,irefz,iax)) &
                -  9.0*(f(l1:l2,m1+2:m2+2,irefz,iax)-f(l1:l2,m1-2:m2-2,irefz,iax)) &
                +      (f(l1:l2,m1+3:m2+3,irefz,iax)-f(l1:l2,m1-3:m2-3,irefz,iax)) )
          else
            if (ip <= 5) print*, 'set_BB2: Degenerate case in y-direction'
          endif
!
          if (nzgrid /= 1) then
            fac = (1./60)*spread(spread(dz_1(irefz),1,nx),2,ny)
            bbx = bbx - fac * &
                ( 45.0*(f(l1:l2,m1:m2,irefz+1,iay)-f(l1:l2,m1:m2,irefz-1,iay)) &
                -  9.0*(f(l1:l2,m1:m2,irefz+2,iay)-f(l1:l2,m1:m2,irefz-2,iay)) &
                +      (f(l1:l2,m1:m2,irefz+3,iay)-f(l1:l2,m1:m2,irefz-2,iay)) )
            bby = bby + fac * &
                ( 45.0*(f(l1:l2,m1:m2,irefz+1,iax)-f(l1:l2,m1:m2,irefz-1,iax)) &
                -  9.0*(f(l1:l2,m1:m2,irefz+2,iax)-f(l1:l2,m1:m2,irefz-2,iax)) &
                +      (f(l1:l2,m1:m2,irefz+3,iax)-f(l1:l2,m1:m2,irefz-3,iax)) )
          else
            if (ip <= 5) print*, 'set_BB2: Degenerate case in z-direction'
          endif
        endif
!
! Compute |B| and collect as global BB2 on root processor.
!
        BB2_local = bbx*bbx + bby*bby + bbz*bbz
        call collect_xy (BB2_local, BB2)
        call sum_xy (sum (abs (bbz)), Bz_total_flux)
        if (nxgrid > 1) Bz_total_flux = Bz_total_flux * dx
        if (nygrid > 1) Bz_total_flux = Bz_total_flux * dy
        deallocate (fac, bbx, bby, bbz)
      else
        Bz_total_flux = 0.0
      endif
!
      if (lgran_parallel) then
        ! Communicate BB2 to granulation computation processors.
        if (lroot) then
          do level = 1, nglevel
            partner = nprocxy + level - 1
            call mpisend_real (BB2, (/nxgrid,nygrid/), partner, BB2_tag)
          enddo
        elseif (lgran_proc) then
          call mpirecv_real (BB2, (/nxgrid,nygrid/), 0, BB2_tag)
        endif
      endif
!
    endsubroutine set_BB2
!***********************************************************************
    subroutine fill_avoid_gran
!
      integer :: i, j, itmp, jtmp
      integer :: il, ir, jl, jr
      integer :: ii, jj
      real :: BB2_limit
!
      avoid_gran = 0
      if (Bavoid <= 0.) return
!
      if (nxgrid == 1) then
        itmp = 0
      else
        itmp = nint(granr*(1-ig)/dx)
      endif
      if (nygrid == 1) then
        jtmp = 0
      else
        jtmp = nint(granr*(1-ig)/dy)
      endif
!
      BB2_limit = (Bavoid/unit_magnetic)**2
      do i = 1,nxgrid
        do j = 1,nygrid
          if (BB2(i,j) > BB2_limit) then
            il = max(1,i-itmp)
            ir = min(nxgrid,i+itmp)
            jl = max(1,j-jtmp)
            jr = min(nygrid,j+jtmp)
!
            do ii = il,ir
              do jj = jl,jr
                if ((ii-i)**2+(jj-j)**2 < itmp**2+jtmp**2) avoid_gran(ii,jj) = 1
              enddo
            enddo
          endif
        enddo
      enddo
!
    endsubroutine fill_avoid_gran
!***********************************************************************
    subroutine force_solar_wind(df,p)
!
      real, dimension(mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      if (n == n2 .and. llast_proc_z) &
          df(l1:l2,m,n2,iuz) = df(l1:l2,m,n2,iuz)-tau_inv*(p%uu(:,3)-u_add)
!
    endsubroutine force_solar_wind
!***********************************************************************
    subroutine get_wind_speed_offset(f)
!
!  Calculates u_0 so that rho*(u+u_0)=massflux.
!  Set 'win' for rho and
!  massflux can be set as fbcz1/2(rho) in run.in.
!
!  18-06-2008/bing: coded
!
      use Mpicomm, only: sum_xy
!
      real, dimension(mx,my,mz,mfarray) :: f
      real :: local_flux, local_mass
      real :: total_flux, total_mass
!
      local_flux = sum(exp(f(l1:l2,m1:m2,n2,ilnrho))*f(l1:l2,m1:m2,n2,iuz))
      local_mass = sum(exp(f(l1:l2,m1:m2,n2,ilnrho)))
!
      call sum_xy (local_flux, total_flux)
      call sum_xy (local_mass, total_mass)
!
!  Get u0 addition rho*(u+u0) = wind
!  rho*u + u0 *rho =wind
!  u0 = (wind-rho*u)/rho
!
      u_add = (massflux-total_flux) / total_mass
!
    endsubroutine get_wind_speed_offset
!***********************************************************************
    subroutine evolve_granules()
!
      integer :: xpos, ypos
!
      current => first
      do while (associated (current))
        xpos = int(current%pos_x)
        ypos = int(current%pos_y)
!
        current%pos_x =  current%pos_x + Ux_ext(xpos,ypos)*dt_gran
        current%pos_y =  current%pos_y + Uy_ext(xpos,ypos)*dt_gran
!
        if (current%pos_x < 0.5) current%pos_x = current%pos_x + nxgrid
        if (current%pos_y < 0.5) current%pos_y = current%pos_y + nygrid
!
        if (current%pos_x > nxgrid+0.5) current%pos_x = current%pos_x - nxgrid
        if (current%pos_y > nygrid+0.5) current%pos_y = current%pos_y - nygrid
!
        current => current%next
      enddo
!
    endsubroutine evolve_granules
!***********************************************************************
    subroutine enhance_vorticity()
!
      real, dimension(nxgrid,nygrid) :: wscr, wscr2
      real :: vrms, vtot
!
! Putting sum of velocities back into vx,vy
      vx = Ux
      vy = Uy
!
! Calculating and enhancing rotational part by factor 5
      if (vorticity_factor > 0.0) then
        call helmholtz(wscr,wscr2)
        vx = vx + vorticity_factor*wscr
        vy = vy + vorticity_factor*wscr2
      endif
!
! Normalize to given total rms-velocity
      vrms = sqrt(sum(vx**2+vy**2)/(nxgrid*nygrid)) + tini
!
      if (unit_system == 'SI') then
        vtot = 3.*1e3/unit_velocity
      elseif (unit_system == 'cgs') then
        vtot = 3.*1e5/unit_velocity
      else
        vtot = 0.
        call fatal_error('solar_corona','define a valid unit system')
      endif
!
! Reinserting rotationally enhanced velocity field
!
      Ux = vx * vtot/vrms
      Uy = vy * vtot/vrms
!
    endsubroutine enhance_vorticity
!***********************************************************************
    subroutine bc_emf_z(f,df,dt_,topbot,j)
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (mx,my,mz,mvar), intent(in) :: df
      real, dimension(mz) :: border_prof_z_uu=1.0
      real, dimension(nz) :: zeta
      real, intent(in) :: dt_
      real :: uborder, border_width
      integer, intent (in) :: j
      integer :: i
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        print*, "bc_emf_z: ", topbot, " should be 'top'"
      case ('top')               ! top boundary
        if (border_frac_z(2) /= 0) then 
          border_width=border_frac_z(2)*Lxyz(3)/2
          uborder=xyz1(3)-border_width
          zeta=1-max(z(n1:n2)-uborder,0.0)/border_width
          border_prof_z_uu(n1:n2)=min(border_prof_z_uu(n1:n2),zeta**2*(3-2*zeta))
          do n=n1, n2
            do i=0, 2
              f(l1:l2,m,n,iuu+i) = f(l1:l2,m,n,iuu+i) - &
              f(l1:l2,m,n,iuu+i)*(1.-border_prof_z_uu(n))*uu_tau1_quench*dt_
            enddo
          enddo
        endif
        if (llast_proc_z) then
          do i=1,nghost
            f(l1:l2,m,n2+i,j)=f(l1:l2,m,n2+i,j)+df(l1:l2,m,n2,j)*dt_
          enddo
        endif

      case default
        print*, "bc_emf_z: ", topbot, " should be 'top'"
      endselect
!
    endsubroutine bc_emf_z
!***********************************************************************

    subroutine generate_halfgrid(x12,y12,z12)
!
! x[l1:l2]+0.5/dx_1[l1:l2]-0.25*dx_tilde[l1:l2]/dx_1[l1:l2]^2
!
      real, dimension (mx), intent(out) :: x12
      real, dimension (my), intent(out) :: y12
      real, dimension (mz), intent(out) :: z12
!
      if (nxgrid == 1) then
        x12 = x
      else
        x12 = x + 0.5/dx_1 - 0.25*dx_tilde/dx_1**2
      endif
!
      if (nygrid == 1) then
        y12 = y
      else
        y12 = y + 0.5/dy_1 - 0.25*dy_tilde/dy_1**2
      endif
!
      if (nzgrid == 1) then
        z12 = z
      else
        z12 = z + 0.5/dz_1 - 0.25*dz_tilde/dz_1**2
      endif
!
    endsubroutine generate_halfgrid
!*******************************************************************************
    subroutine div_diff_flux(f,j,p,div_flux)
      intent(in) :: f,j
      intent(out) :: div_flux
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: flux_im12,flux_ip12
      real, dimension (nx,3) :: cmax_im12,cmax_ip12
      real, dimension (nx) :: div_flux
      real, dimension (nx) :: fim12_l,fim12_r,fip12_l,fip12_r
      real, dimension (nx) :: fim1,fip1,rfac,q1
      real :: alpha=1.25
      type (pencil_case), intent(in) :: p
      integer :: j,k
!
! First set the diffusive flux = cmax*(f_R-f_L) at half grid points
!
        fim1=0
        fip1=0
        fim12_l=0
        fim12_r=0
        fip12_l=0
        fip12_r=0
        cmax_ip12=0
        cmax_im12=0
!
!  Generate halfgrid points
!
      do k=1,3
!
        call slope_lim_lin_interpol(f,j,fim12_l,fim12_r,fip12_l,fip12_r,&
                                   fim1,fip1,k)      
        call characteristic_speed(f,p,cmax_im12,cmax_ip12,k)
        do ix=1,nx
          if (ldensity_nolog) then
            rfac(ix)=abs(fim12_r(ix)-fim12_l(ix))/(abs(f(ix+nghost,m,n,j)-&
                     fim1(ix))+tini)
          else
            rfac(ix)=abs(fim12_r(ix)-fim12_l(ix))/(abs(exp(f(ix+nghost,m,n,j))-&
                     fim1(ix))+tini)
          endif
          q1(ix)=(min1(1.0,alpha*rfac(ix)))**nlf
        enddo
        flux_im12(:,k)=0.5*cmax_im12(:,k)*q1*(fim12_r-fim12_l)
        do ix=1,nx
          if (ldensity_nolog) then
            rfac(ix)=abs(fip12_r(ix)-fip12_l(ix))/(abs(fip1(ix)-&
                     f(ix+nghost,m,n,j))+tini)
          else
            rfac(ix)=abs(fip12_r(ix)-fip12_l(ix))/(abs(fip1(ix)-&
                     exp(f(ix+nghost,m,n,j)))+tini)
          endif
          q1(ix)=(min1(1.0,alpha*rfac(ix)))**nlf
        enddo
        flux_ip12(:,k)=0.5*cmax_ip12(:,k)*q1*(fip12_r-fip12_l)
      enddo
!
! Now calculate the 2nd order divergence
!
      if (lspherical_coords) then
        div_flux=0.0
        div_flux=div_flux+(x12(l1:l2)**2*flux_ip12(:,1)-x12(l1-1:l2-1)**2*flux_im12(:,1))&
                 /(x(l1:l2)**2*(x12(l1:l2)-x12(l1-1:l2-1))) &
        +(sin(y(m)+0.5*dy)*flux_ip12(:,2)-&
                sin(y(m)-0.5*dy)*flux_im12(:,2))/(x(l1:l2)*sin(y(m))*dy) &
            +(flux_ip12(:,3)-flux_im12(:,3))/(x(l1:l2)*sin(y(m))*dz)
      else if (lcartesian_coords) then
        div_flux=0.0
        if (nygrid == 1) then
          div_flux=div_flux+(flux_ip12(:,1)-flux_im12(:,1))&
                   /(x12(l1:l2)-x12(l1-1:l2-1)) &
              +(flux_ip12(:,3)-flux_im12(:,3))/(z12(n)-z12(n-1))
        else
          div_flux=div_flux+(flux_ip12(:,1)-flux_im12(:,1))&
                   /(x12(l1:l2)-x12(l1-1:l2-1)) &
          +(flux_ip12(:,2)-flux_im12(:,2))/(y12(m)-y12(m-1)) &
            +(flux_ip12(:,3)-flux_im12(:,3))/(z12(n)-z12(n-1))
        endif
       else
         call fatal_error('solar_corona:div_diff_flux','Not coded for cylindrical')
       endif
!
    endsubroutine div_diff_flux
!*******************************************************************************
    subroutine slope_lim_lin_interpol(f,j,fim12_l,fim12_r,fip12_l,fip12_r,fim1,&
                                     fip1,k)
!
! Reconstruction of a scalar by slope limited linear interpolation
! Get values at half grid points l+1/2,m+1/2,n+1/2 depending on case(k)
!
      intent(in) :: f,k,j
      intent(out) :: fim12_l,fim12_r,fip12_l,fip12_r
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: fim12_l,fim12_r,fip12_l,fip12_r,fim1,fip1
      real, dimension (nx) :: delfy,delfz,delfyp1,delfym1,delfzp1,delfzm1
      real, dimension (0:nx+1) :: delfx
      integer :: j,k
      integer :: i,ix
      real :: tmp0,tmp1,tmp2,tmp3
! x-component
      select case (k)
        case(1)
! Special case of non uniform grid only in the radial direction
        if (ldensity_nolog) then
          do i=l1-1,l2+1
            ix=i-nghost
            tmp1=(f(i,m,n,j)-f(i-1,m,n,j))/(x(i)-x(i-1))
            tmp2=(f(i+1,m,n,j)-f(i,m,n,j))/(x(i+1)-x(i))
            delfx(ix) = minmod_alt(tmp1,tmp2)
          enddo
          do i=l1,l2
            ix=i-nghost
            fim12_l(ix) = f(i-1,m,n,j)+delfx(ix-1)*(x(i)-x(i-1))
            fim12_r(ix) = f(i,m,n,j)-delfx(ix)*(x(i)-x(i-1))
            fip12_l(ix) = f(i,m,n,j)+delfx(ix)*(x(i+1)-x(i))
            fip12_r(ix) = f(i+1,m,n,j)-delfx(ix+1)*(x(i+1)-x(i))
            fim1(ix) = f(i-1,m,n,j)
            fip1(ix) = f(i+1,m,n,j)
          enddo
        else
          do i=l1-1,l2+1
            ix=i-nghost
            tmp1=(exp(f(i,m,n,j))-exp(f(i-1,m,n,j)))/(x(i)-x(i-1))
            tmp2=(exp(f(i+1,m,n,j))-exp(f(i,m,n,j)))/(x(i+1)-x(i))
            delfx(ix) = minmod_alt(tmp1,tmp2)
          enddo
          do i=l1,l2
            ix=i-nghost
            fim12_l(ix) = exp(f(i-1,m,n,j))+delfx(ix-1)*(x(i)-x(i-1))
            fim12_r(ix) = exp(f(i,m,n,j))-delfx(ix)*(x(i)-x(i-1))
            fip12_l(ix) = exp(f(i,m,n,j))+delfx(ix)*(x(i+1)-x(i))
            fip12_r(ix) = exp(f(i+1,m,n,j))-delfx(ix+1)*(x(i+1)-x(i))
            fim1(ix) = exp(f(i-1,m,n,j))
            fip1(ix) = exp(f(i+1,m,n,j))
          enddo
        endif
! y-component
        case(2)
        if (ldensity_nolog) then
          do i=l1,l2
            ix=i-nghost
            tmp0=f(i,m-1,n,j)-f(i,m-2,n,j)
            tmp1=f(i,m,n,j)-f(i,m-1,n,j)
            delfym1(ix) = minmod_alt(tmp0,tmp1)
            tmp2=f(i,m+1,n,j)-f(i,m,n,j)
            delfy(ix) = minmod_alt(tmp1,tmp2)
            tmp3=f(i,m+2,n,j)-f(i,m+1,n,j)
            delfyp1(ix) = minmod_alt(tmp2,tmp3)
          enddo
          do i=l1,l2
            ix=i-nghost
            fim12_l(ix) = f(i,m-1,n,j)+delfym1(ix)
            fim12_r(ix) = f(i,m,n,j)-delfy(ix)
            fip12_l(ix) = f(i,m,n,j)+delfy(ix)
            fip12_r(ix) = f(i,m+1,n,j)-delfyp1(ix)
            fim1(ix) = f(i,m-1,n,j)
            fip1(ix) = f(i,m+1,n,j)
          enddo
        else
          do i=l1,l2
            ix=i-nghost
            tmp0=exp(f(i,m-1,n,j))-exp(f(i,m-2,n,j))
            tmp1=exp(f(i,m,n,j))-exp(f(i,m-1,n,j))
            delfym1(ix) = minmod_alt(tmp0,tmp1)
            tmp2=exp(f(i,m+1,n,j))-exp(f(i,m,n,j))
            delfy(ix) = minmod_alt(tmp1,tmp2)
            tmp3=exp(f(i,m+2,n,j))-exp(f(i,m+1,n,j))
            delfyp1(ix) = minmod_alt(tmp2,tmp3)
          enddo
          do i=l1,l2
            ix=i-nghost
            fim12_l(ix) = exp(f(i,m-1,n,j))+delfym1(ix)
            fim12_r(ix) = exp(f(i,m,n,j))-delfy(ix)
            fip12_l(ix) = exp(f(i,m,n,j))+delfy(ix)
            fip12_r(ix) = exp(f(i,m+1,n,j))-delfyp1(ix)
            fim1(ix) = exp(f(i,m-1,n,j))
            fip1(ix) = exp(f(i,m+1,n,j))
          enddo
        endif
! z-component
        case(3)
        if (ldensity_nolog) then
          do i=l1,l2
            ix=i-nghost
            tmp0=f(i,m,n-1,j)-f(i,m,n-2,j)
            tmp1=f(i,m,n,j)-f(i,m,n-1,j)
            delfzm1(ix) = minmod_alt(tmp0,tmp1)
            tmp2=f(i,m,n+1,j)-f(i,m,n,j)
            delfz(ix) = minmod_alt(tmp1,tmp2)
            tmp3=f(i,m,n+2,j)-f(i,m,n+1,j)
            delfzp1(ix) = minmod_alt(tmp2,tmp3)
          enddo
          do i=l1,l2
            ix=i-nghost
            fim12_l(ix) = f(i,m,n-1,j)+delfzm1(ix)
            fim12_r(ix) = f(i,m,n,j)-delfz(ix)
            fip12_l(ix) = f(i,m,n,j)+delfz(ix)
            fip12_r(ix) = f(i,m,n+1,j)-delfzp1(ix)
            fim1(ix) = f(i,m,n-1,j)
            fip1(ix) = f(i,m,n+1,j)
          enddo
        else
          do i=l1,l2
            ix=i-nghost
            tmp0=exp(f(i,m,n-1,j))-exp(f(i,m,n-2,j))
            tmp1=exp(f(i,m,n,j))-exp(f(i,m,n-1,j))
            delfzm1(ix) = minmod_alt(tmp0,tmp1)
            tmp2=exp(f(i,m,n+1,j))-exp(f(i,m,n,j))
            delfz(ix) = minmod_alt(tmp1,tmp2)
            tmp3=exp(f(i,m,n+2,j))-exp(f(i,m,n+1,j))
            delfzp1(ix) = minmod_alt(tmp2,tmp3)
          enddo
          do i=l1,l2
            ix=i-nghost
            fim12_l(ix) = exp(f(i,m,n-1,j))+delfzm1(ix)
            fim12_r(ix) = exp(f(i,m,n,j))-delfz(ix)
            fip12_l(ix) = exp(f(i,m,n,j))+delfz(ix)
            fip12_r(ix) = exp(f(i,m,n+1,j))-delfzp1(ix)
            fim1(ix) = exp(f(i,m,n-1,j))
            fip1(ix) = exp(f(i,m,n+1,j))
          enddo
        endif
        case default
          call fatal_error('twist_inject:slope_lim_lin_interpol','wrong component')
        endselect
!
    endsubroutine slope_lim_lin_interpol
!***********************************************************************
    subroutine characteristic_speed(f,p,cmax_im12,cmax_ip12,k)
      intent(in) :: f,k
      intent(out) :: cmax_im12,cmax_ip12
!
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: cmax_im12,cmax_ip12
      real, dimension (nx) :: b1_tmp,b2_tmp,b3_tmp,rho_tmp
      real, dimension (0:nx) :: b1_xtmp,b2_xtmp,b3_xtmp,rho_xtmp
      type (pencil_case), intent(in) :: p
      integer :: k
!
      select case (k)
        case(1)
          ! x-component
          b1_xtmp=0.5*(f(l1-1:l2,m,n,ibx)+f(l1:l2+1,m,n,ibx))
          b2_xtmp=0.5*(f(l1-1:l2,m,n,iby)+f(l1:l2+1,m,n,iby))
          b3_xtmp=0.5*(f(l1-1:l2,m,n,ibz)+f(l1:l2+1,m,n,ibz))
          if (ldensity_nolog) then
            rho_xtmp=0.5*(f(l1-1:l2,m,n,irho)+f(l1:l2+1,m,n,irho))
          else
            rho_xtmp=0.5*(exp(f(l1-1:l2,m,n,ilnrho))+exp(f(l1:l2+1,m,n,ilnrho)))
          endif
          cmax_im12(:,1)=sqrt(b1_xtmp(0:nx-1)**2+b2_xtmp(0:nx-1)**2+b3_xtmp(0:nx-1)**2)/sqrt(rho_xtmp(0:nx-1))+sqrt(p%cs2)
          cmax_ip12(:,1)=sqrt(b1_xtmp(1:nx)**2+b2_xtmp(1:nx)**2+b3_xtmp(1:nx)**2)/sqrt(rho_xtmp(1:nx))+sqrt(p%cs2)
        case(2)
          ! y-component
          b1_tmp=0.5*(f(l1:l2,m-1,n,ibx)+f(l1:l2,m,n,ibx))
          b2_tmp=0.5*(f(l1:l2,m-1,n,iby)+f(l1:l2,m,n,iby))
          b3_tmp=0.5*(f(l1:l2,m-1,n,ibz)+f(l1:l2,m,n,ibz))
          if (ldensity_nolog) then
            rho_tmp=0.5*(f(l1:l2,m-1,n,irho)+f(l1:l2,m,n,irho))
          else
            rho_tmp=0.5*(exp(f(l1:l2,m-1,n,ilnrho))+exp(f(l1:l2,m,n,ilnrho)))
          endif
          cmax_im12(:,2)=sqrt(b1_tmp**2+b2_tmp**2+b3_tmp**2)/sqrt(rho_tmp)+sqrt(p%cs2)
          b1_tmp=0.5*(f(l1:l2,m,n,ibx)+f(l1:l2,m+1,n,ibx))
          b2_tmp=0.5*(f(l1:l2,m,n,iby)+f(l1:l2,m+1,n,iby))
          b3_tmp=0.5*(f(l1:l2,m,n,ibz)+f(l1:l2,m+1,n,ibz))
          if (ldensity_nolog) then
            rho_tmp=0.5*(f(l1:l2,m,n,irho)+f(l1:l2,m+1,n,irho))
          else
            rho_tmp=0.5*(exp(f(l1:l2,m,n,ilnrho))+exp(f(l1:l2,m+1,n,ilnrho)))
          endif
          cmax_ip12(:,2)=sqrt(b1_tmp**2+b2_tmp**2+b3_tmp**2)/sqrt(rho_tmp)+sqrt(p%cs2)
        case(3)
          ! z-component
          b1_tmp=0.5*(f(l1:l2,m,n-1,ibx)+f(l1:l2,m,n,ibx))
          b2_tmp=0.5*(f(l1:l2,m,n-1,iby)+f(l1:l2,m,n,iby))
          b3_tmp=0.5*(f(l1:l2,m,n-1,ibz)+f(l1:l2,m,n,ibz))
          if (ldensity_nolog) then
            rho_tmp=0.5*(f(l1:l2,m,n-1,irho)+f(l1:l2,m,n,irho))
          else
            rho_tmp=0.5*(exp(f(l1:l2,m,n-1,ilnrho))+exp(f(l1:l2,m,n,ilnrho)))
          endif
          cmax_im12(:,3)=sqrt(b1_tmp**2+b2_tmp**2+b3_tmp**2)/sqrt(rho_tmp)+sqrt(p%cs2)
          b1_tmp=0.5*(f(l1:l2,m,n,ibx)+f(l1:l2,m,n+1,ibx))
          b2_tmp=0.5*(f(l1:l2,m,n,iby)+f(l1:l2,m,n+1,iby))
          b3_tmp=0.5*(f(l1:l2,m,n,ibz)+f(l1:l2,m,n+1,ibz))
          if (ldensity_nolog) then
            rho_tmp=0.5*(f(l1:l2,m,n,irho)+f(l1:l2,m,n+1,irho))
          else
            rho_tmp=0.5*(exp(f(l1:l2,m,n,ilnrho))+exp(f(l1:l2,m,n+1,ilnrho)))
          endif
          cmax_ip12(:,3)=sqrt(b1_tmp**2+b2_tmp**2+b3_tmp**2)/sqrt(rho_tmp)+sqrt(p%cs2)
        endselect
    endsubroutine characteristic_speed
!***********************************************************************
    real function minmod_alt(a,b)
!
!  minmod=max(0,min(abs(a),sign(1,a)*b))
!
      real :: a,b
!
      minmod_alt=sign(0.5,a)*max(0.0,min(abs(a),sign(1.0,a)*b))
!
    endfunction minmod_alt
!***********************************************************************
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
