! $Id$
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
!***************************************************************
!
module Special
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../special.h'
!
  integer, parameter :: max_gran_levels=3
!
  real :: tdown=0.,allp=0.,Kgpara=0.,cool_RTV=0.,Kgpara2=0.,tdownr=0.,allpr=0.
  real :: lntt0=0.,wlntt=0.,bmdi=0.,hcond1=0.,heatexp=0.,heatamp=0.,Ksat=0.
  real :: diffrho_hyper3=0.,chi_hyper3=0.,chi_hyper2=0.,K_iso=0.
  real :: Bavoid=0.,nvor=5.,tau_inv=1.,Bz_flux=0.,q0=1.,qw=1.,dq=0.1,dt_gran=0.
  logical :: lgranulation=.false.,lrotin=.true.,lgran_proc=.false.
  logical :: luse_ext_vel_field=.false.,lquench=.false.,lmassflux=.false.
  integer :: irefz=n1,nglevel=max_gran_levels,cool_type=2
  real :: massflux=0.,u_add,hcond2=0.,hcond3=0.,init_time=0.
!
  real, dimension (nx,ny) :: A_init_x, A_init_y
!
  character (len=labellen), dimension(3) :: iheattype='nothing'
  real, dimension(2) :: heat_par_exp=(/0.,1./)
  real, dimension(2) :: heat_par_exp2=(/0.,1./)
  real, dimension(3) :: heat_par_gauss=(/0.,1.,0./)
!
  character (len=labellen) :: prof_type='nothing'
  real, dimension (mz) :: uu_init_z, lnrho_init_z, lnTT_init_z
  logical :: linit_uu=.false., linit_lnrho=.false., linit_lnTT=.false.
!
! input parameters
  namelist /special_init_pars/ &
       linit_uu,linit_lnrho,linit_lnTT,prof_type
!
! run parameters
  namelist /special_run_pars/ &
       tdown,allp,Kgpara,cool_RTV,lntt0,wlntt,bmdi,hcond1,Kgpara2, &
       tdownr,allpr,heatexp,heatamp,Ksat,diffrho_hyper3, &
       chi_hyper3,chi_hyper2,K_iso,lgranulation,irefz, &
       Bavoid,nglevel,lrotin,nvor,tau_inv,Bz_flux,init_time, &
       lquench,q0,qw,dq,massflux,luse_ext_vel_field,prof_type, &
       lmassflux,hcond2,hcond3,heat_par_gauss,heat_par_exp,heat_par_exp2, &
       iheattype,dt_gran,cool_type
!
    integer :: idiag_dtnewt=0   ! DIAG_DOC: Radiative cooling time step
    integer :: idiag_dtchi2=0   ! DIAG_DOC: $\delta t / [c_{\delta t,{\rm v}}\,
                                ! DIAG_DOC:   \delta x^2/\chi_{\rm max}]$
                                ! DIAG_DOC:   \quad(time step relative to time
                                ! DIAG_DOC:   step based on heat conductivity;
                                ! DIAG_DOC:   see \S~\ref{time-step})
!
! video slices
    real, target, dimension (nx,ny) :: rtv_xy,rtv_xy2,rtv_xy3,rtv_xy4
    real, target, dimension (nx,nz) :: rtv_xz
    real, target, dimension (ny,nz) :: rtv_yz
    real, target, dimension (nx,ny) :: logQ_xy,logQ_xy2,logQ_xy3,logQ_xy4
    real, target, dimension (nx,nz) :: logQ_xz
    real, target, dimension (ny,nz) :: logQ_yz
!
!
    TYPE point
      real, dimension(2) :: pos
      real, dimension(4) :: data
      type(point),pointer :: next
      type(point),pointer :: previous
    end TYPE point
!
    type(point), pointer :: first => null()
    type(point), pointer :: current => null()
    type(point), pointer :: firstlev => null()
    type(point), pointer :: secondlev => null()
    type(point), pointer :: thirdlev => null()
!
    integer :: xrange,yrange,pow
    real :: ampl,dxdy2,ig,granr,pd,life_t,upd,avoid
    real, dimension(:,:), allocatable :: w,vx,vy
    real, dimension(:,:), allocatable :: Ux,Uy
    real, dimension(:,:), allocatable :: Ux_ext,Uy_ext
    real, dimension(:,:), allocatable :: BB2
    integer, dimension(:,:), allocatable :: avoidarr
    real, save :: tsnap_uu=0.,thresh
    integer, save :: isnap
    integer, save, dimension(mseed) :: points_rstate
    real, dimension(nx,ny), save :: ux_local,uy_local
    real, dimension(nx,ny), save :: ux_ext_local,uy_ext_local
    real :: Bzflux
    logical :: lgran_parallel=.false.
!
    integer, save, dimension(mseed) :: nano_seed
    integer :: alloc_err
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
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
! Called by start.f90 with lstarting=.true. or by
! run.f90 with lstarting=.false. and with lreloading indicating a RELOAD
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      if (lgranulation) then
!
! Consistency checks:
!
        if ((nglevel < 1) .or. (nglevel > max_gran_levels)) &
            call fatal_error ('initialize_special', 'nglevel is too large!')
        ! For only one granulation level, no parallelization is required.
        if (lgran_parallel .and. (nglevel == 1)) &
            call fatal_error ('initialize_special', &
                'if nglevel is 1, lgran_parallel should be set to false.')
        ! If not at least 3 procs above the ipz=0 plane are available,
        ! computing of granular velocities has to be done non-parallel.
        if (lgran_parallel .and. (nprocz-1)*nprocxy < 3) &
            call fatal_error ('initialize_special', &
                'you have not enough processors to activate lgran_parallel.')
!
! Define and initialize the processors that are computing the granulation:
!
        lgran_proc = (lroot .and. .not. lgran_parallel) .or. &
            (lgran_parallel .and. (iproc >= nprocxy) .and. (iproc < nprocxy+nglevel))
!
        if (lroot .or. lgran_proc) then
          call setdrparams()
        endif
      endif
!
      if ((.not. lreloading) .and. (.not. lstarting)) nano_seed = 0.
!
      call setup_magnetic()
      call setup_profiles()
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
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
      use Mpicomm, only: stop_it_if_any
      use Messages, only: warning
!
      real, dimension (mx,my,mz,mfarray), intent (out) :: f
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
        ! set initial density profile values
        do j = 1, mz
          f(:,:,j,ilnrho) = lnrho_init_z(j)
        enddo
      endif
!
      if (linit_lnTT) then
        if (pretend_lnTT) call stop_it_if_any (.true., &
            "init_special: linit_lnTT=T not implemented for pretend_lnTT=T")
        ! set initial temperaure profile values
        do j = 1, mz
          if (ltemperature) then
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
    subroutine setup_magnetic()
!
!  Compute and save initial magnetic vector potential A_init_x/_y.
!
!  25-mar-10/Bourdin.KIS: coded
!
      use Fourier, only: fourier_transform_other
      use Mpicomm, only: mpisend_real, mpirecv_real, stop_it_if_any
      use Syscalls, only: file_exists
!
      real, dimension(:,:), allocatable :: kx, ky, k2
      real, dimension(:,:), allocatable :: Bz0_i, Bz0_r
      real, dimension(:,:), allocatable :: Ax_i, Ay_i
      real, dimension(:,:), allocatable :: Ax_r, Ay_r
      real, dimension(:), allocatable :: kxp, kyp
!
      real :: dummy
      integer :: idx2,idy2,lend,ierr
      integer :: i,px,py
      integer, parameter :: unit=12,Ax_tag=366,Ay_tag=367
!
      ! file location settings
      character (len=*), parameter :: mag_field_txt = 'driver/mag_field.txt'
      character (len=*), parameter :: mag_field_dat = 'driver/mag_field.dat'
!
!
      ! don't touch anything during a RELOAD:
      if (lreloading) return
!
      inquire (IOLENGTH=lend) dummy
!
      if (.not. lfirst_proc_z .or. (bmdi == 0.0)) then
        A_init_x = 0.0
        A_init_y = 0.0
      else
        ! Magnetic field is set only in the bottom layer
        if (lroot) then
          allocate(kx(nxgrid,nygrid), ky(nxgrid,nygrid), k2(nxgrid,nygrid), kxp(nxgrid), kyp(nygrid), stat=alloc_err)
          if (alloc_err > 0) call stop_it_if_any (.true., 'setup_magnetic: '// &
              'Could not allocate memory for wave vector variables')
          allocate(Bz0_r(nxgrid,nygrid), Bz0_i(nxgrid,nygrid), stat=alloc_err)
          if (alloc_err > 0) call stop_it_if_any (.true., 'setup_magnetic: '// &
              'Could not allocate memory for vertical magnetic field variables')
          allocate(Ax_r(nxgrid,nygrid), Ay_r(nxgrid,nygrid), Ax_i(nxgrid,nygrid), Ay_i(nxgrid,nygrid), stat=alloc_err)
          if (alloc_err > 0) call stop_it_if_any (.true., 'setup_magnetic: '// &
              'Could not allocate memory for vector potential A variables')
          ! Auxiliary quantities:
          ! idx2 and idy2 are essentially =2, but this makes compilers
          ! complain if nygrid=1 (in which case this is highly unlikely to be
          ! correct anyway), so we try to do this better:
          idx2 = min(2,nxgrid)
          idy2 = min(2,nygrid)
!
          kxp=cshift((/(i-(nxgrid-1)/2,i=0,nxgrid-1)/),+(nxgrid-1)/2)*2*pi/Lx
          kyp=cshift((/(i-(nygrid-1)/2,i=0,nygrid-1)/),+(nygrid-1)/2)*2*pi/Ly
!
          kx=spread(kxp,2,nygrid)
          ky=spread(kyp,1,nxgrid)
!
          k2 = kx*kx + ky*ky
!
          ! Read in magnetogram
          if (file_exists(mag_field_txt)) then
            open (unit,file=mag_field_txt)
            read (unit,*,iostat=ierr) Bz0_r
            if (ierr /= 0) call stop_it_if_any(.true.,'setup_magnetic: '// &
                'Error reading magnetogram file: "'//trim(mag_field_txt)//'"')
            close (unit)
          elseif (file_exists(mag_field_dat)) then
            open (unit,file=mag_field_dat,form='unformatted',status='unknown', &
                recl=lend*nxgrid*nygrid,access='direct')
            read (unit,rec=1,iostat=ierr) Bz0_r
            if (ierr /= 0) call stop_it_if_any(.true.,'setup_magnetic: '// &
                'Error reading magnetogram file: "'//trim(mag_field_dat)//'"')
            close (unit)
          else
            call stop_it_if_any(.true., 'setup_magnetic: No magnetogram file found.')
          endif
!
          ! Gauss to Tesla and SI to PENCIL units
          Bz0_r = Bz0_r * 1e-4 / unit_magnetic
          Bz0_i = 0.
!
          ! Fourier Transform of Bz0:
          call fourier_transform_other(Bz0_r,Bz0_i)
!
          where (k2 /= 0)
            Ax_r = -Bz0_i*ky/k2*exp(-sqrt(k2)*z(n1) )
            Ax_i =  Bz0_r*ky/k2*exp(-sqrt(k2)*z(n1) )
!
            Ay_r =  Bz0_i*kx/k2*exp(-sqrt(k2)*z(n1) )
            Ay_i = -Bz0_r*kx/k2*exp(-sqrt(k2)*z(n1) )
          elsewhere
            Ax_r = -Bz0_i*ky/ky(1,idy2)*exp(-sqrt(k2)*z(n1) )
            Ax_i =  Bz0_r*ky/ky(1,idy2)*exp(-sqrt(k2)*z(n1) )
!
            Ay_r =  Bz0_i*kx/kx(idx2,1)*exp(-sqrt(k2)*z(n1) )
            Ay_i = -Bz0_r*kx/kx(idx2,1)*exp(-sqrt(k2)*z(n1) )
          endwhere
!
          deallocate(kx, ky, k2, kxp, kyp, Bz0_i, Bz0_r)
!
          call fourier_transform_other(Ax_r,Ax_i,linv=.true.)
          call fourier_transform_other(Ay_r,Ay_i,linv=.true.)
!
          ! Distribute initial A data
          do px=0, nprocx-1
            do py=0, nprocy-1
              if ((px /= 0) .or. (py /= 0)) then
                A_init_x = Ax_r(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
                A_init_y = Ay_r(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
                call mpisend_real (A_init_x, (/ nx, ny /), px+py*nprocx, Ax_tag)
                call mpisend_real (A_init_y, (/ nx, ny /), px+py*nprocx, Ay_tag)
              endif
            enddo
          enddo
          A_init_x = Ax_r(1:nx,1:ny)
          A_init_y = Ay_r(1:nx,1:ny)
!
          deallocate(Ax_r, Ay_r, Ax_i, Ay_i)
!
        else
          ! Receive initial A data
          call mpirecv_real (A_init_x, (/ nx, ny /), 0, Ax_tag)
          call mpirecv_real (A_init_y, (/ nx, ny /), 0, Ay_tag)
        endif
      endif
      ! globally catch eventual 'stop_it_if_any' calls from single MPI ranks
      call stop_it_if_any (.false.,'')
!
    endsubroutine setup_magnetic
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
!
!
      logical :: lnewton_cooling=.false.
!
! Only read profiles if needed, e.g. for Newton cooling
!
      if (.not. (ltemperature .or. lentropy)) return
      lnewton_cooling = (tdown/=0) .or. (tdownr/=0)
      if (.not. (linit_uu .or. linit_lnrho .or. linit_lnTT .or. lnewton_cooling)) return
!
      ! default: read 'stratification.dat' with density and temperature
      if (prof_type == 'nothing') prof_type = 'lnrho_lnTT'
!
      ! check if density profile is read, when needed
      if ((tdownr /= 0.0) .and. (index (prof_type, 'lnrho') < 1)) then
        call fatal_error ("setup_profiles", &
            "a density profile must be read to use density based newton cooling")
      endif
!
      ! on RELOAD we don't need to read the profiles again
      if (.not. lreloading) call read_profiles()
!
      ! check if any kind of density profile is in use...
      if (linit_lnrho .or. (tdownr /= 0.0)) then
        ! ...and set lnrho0 accordingly to the lower boundary value
        if ((lnrho0 /= 0.0) .and. (lnrho0 /= lnrho_init_z(n1))) then
          if (lroot) print *,'setup_profiles: WARNING: ', &
              'lnrho0 set to ', lnrho_init_z(n1), ' - was before ', lnrho0
          call warning ("setup_profiles", "overriding manual lnrho0 setting")
        endif
        lnrho0 = lnrho_init_z(n1)
        if ((rho0 /= 1.0) .and. (abs (rho0 / exp (lnrho0) - 1.0) > 1.e-6)) then
          if (lroot) print *,'setup_profiles: WARNING: ', &
              'rho0 set to ', exp (lnrho0), ' - was before ', rho0
          call warning ("setup_profiles", "overriding manual rho0 setting")
        endif
        rho0 = exp (lnrho0)
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
      use Mpicomm, only: mpibcast_real, stop_it_if_any
      use Syscalls, only: file_exists
!
      integer :: i, ierr
      integer, parameter :: unit=12, lnrho_tag=368, lnTT_tag=369
      real :: var0, var1, var2
      real, dimension (:), allocatable :: prof_lnrho, prof_lnTT
      logical :: lread_prof_uu, lread_prof_lnrho, lread_prof_lnTT
!
      ! file location settings
      character (len=*), parameter :: stratification_dat = 'stratification.dat'
      character (len=*), parameter :: lnrho_dat = 'driver/prof_lnrho.dat'
      character (len=*), parameter :: lnT_dat = 'driver/prof_lnT.dat'
      character (len=*), parameter :: uz_dat = 'driver/prof_uz.dat'
!
!
! Check which stratification file should be used:
!
      lread_prof_uu    = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_uu') > 0)
      lread_prof_lnrho = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_lnrho') > 0)
      lread_prof_lnTT  = (index (prof_type, 'prof_') == 1) .and. (index (prof_type, '_lnTT') > 0)
!
      if (prof_type=='lnrho_lnTT') then
        allocate (prof_lnTT(nzgrid), prof_lnrho(nzgrid), stat=ierr)
        if (ierr > 0) call stop_it_if_any (.true., 'setup_profiles: '// &
            'Could not allocate memory for stratification variables')
!
        ! read stratification file only on the MPI root rank
        if (lroot) then
          if (.not. file_exists (stratification_dat)) call stop_it_if_any ( &
              .true., 'setup_profiles: Stratification file not found')
          open (unit,file=stratification_dat)
          do i=1,nzgrid
            read (unit,*,iostat=ierr) var0,var1,var2
            if (ierr /= 0) call stop_it_if_any (.true., 'setup_profiles: '// &
                'Error reading stratification file: "'//trim(stratification_dat)//'"')
            prof_lnrho(i)=var1
            prof_lnTT(i) =var2
          enddo
          close(unit)
        endif
        call stop_it_if_any (.false.,'')
!
        call mpibcast_real (prof_lnTT,nzgrid)
        call mpibcast_real (prof_lnrho,nzgrid)
!
        lnTT_init_z(n1:n2) = prof_lnTT(ipz*nz+1:(ipz+1)*nz)
        lnrho_init_z(n1:n2) = prof_lnrho(ipz*nz+1:(ipz+1)*nz)
!
        deallocate (prof_lnTT, prof_lnrho)
!
      elseif (lread_prof_uu .or. lread_prof_lnrho .or. lread_prof_lnTT) then
!
        ! read vertical velocity profile for interpolation
        if (lread_prof_uu) &
            call read_profile (uz_dat, uu_init_z, unit_velocity, .false.)
!
        ! read logarithmic density profile for interpolation
        if (lread_prof_lnrho) &
            call read_profile (lnrho_dat, lnrho_init_z, unit_density, .true.)
!
        ! read logarithmic temperature profile
        if (lread_prof_lnTT) &
            call read_profile (lnT_dat, lnTT_init_z, unit_temperature, .true.)
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
      use Mpicomm, only: mpibcast_int, mpibcast_real, stop_it_if_any
      use Syscalls, only: file_exists, file_size
!
      character (len=*), intent (in) :: filename
      real, dimension (mz), intent (out) :: profile
      real, intent (in) :: data_unit
      logical, intent (in) :: llog
!
      real, dimension (:), allocatable :: data, data_z
      integer :: n_data
!
      integer, parameter :: unit=12
      integer :: lend, lend_b8, ierr
!
!
      inquire (IOLENGTH=lend) 1.0
      inquire (IOLENGTH=lend_b8) 1.0d0
!
      ! file access is only done on the MPI root rank
      if (lroot) then
        ! determine the number of data points in the profile
        if (.not. file_exists (filename)) &
            call stop_it_if_any (.true., "read_profile: can't find "//filename)
        n_data = (file_size (filename) - 2*2*4) / (lend*8/lend_b8 * 2)
      endif
      call stop_it_if_any (.false., '')
      call mpibcast_int (n_data, 1)
!
      ! allocate memory
      allocate (data(n_data), data_z(n_data), stat=ierr)
      if (ierr > 0) call stop_it_if_any (.true., 'read_profile: '// &
          'Could not allocate memory for data and its z coordinate')
!
      if (lroot) then
        ! read profile
        open (unit, file=filename, form='unformatted', recl=lend*n_data)
        read (unit, iostat=ierr) data
        read (unit, iostat=ierr) data_z
        if (ierr /= 0) call stop_it_if_any (.true., 'read_profile: '// &
            'Error reading profile data in "'//trim(filename)//'"')
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
      call stop_it_if_any (.false., '')
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
      real, dimension (:) :: data, data_z
      integer :: n_data
      real, dimension (mz) :: profile
!
      integer :: i, j, num_over, num_below
      character (len=12) :: num
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
      if (num_below > 0) then
        write (num, *) num_below
        call warning ("interpolate_profile", "extrapolated "//trim (adjustl (num))//" grid points below bottom")
      endif
      if (num_over > 0) then
        write (num, *) num_over
        call warning ("interpolate_profile", "extrapolated "//trim (adjustl (num))//" grid points over top")
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
!
      if (cool_RTV/=0) then
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_cp1)=.true.
      endif
!
      if ((tdown/=0.0) .or. (tdownr/=0.0)) then
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_lnTT)=.true.
      endif
!
      if (hcond1/=0) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
        lpenc_requested(i_lnrho)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (hcond2/=0) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (hcond3/=0) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_glnrho)=.true.
      endif
!
      if (K_iso/=0) then
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_lnTT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_del2lnTT)=.true.
      endif
!
      if (Kgpara/=0) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_bij)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_glnTT)=.true.
        lpenc_requested(i_hlnTT)=.true.
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_rho1)=.true.
      endif
!
      if (idiag_dtchi2/=0) then
        lpenc_diagnos(i_rho1)=.true.
        lpenc_diagnos(i_cv1) =.true.
        lpenc_diagnos(i_cs2)=.true.
      endif
!
      do i=1,3
        select case(iheattype(i))
        case ('sven')
          lpenc_diagnos(i_cp1)=.true.
          lpenc_diagnos(i_TT1)=.true.
          lpenc_diagnos(i_rho1)=.true.
        endselect
      enddo
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif
!
 99    return
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
      if (Kgpara2/=0) then
        if (K_iso/=0) then
          call fatal_error('calc_heatcond_grad', &
              'Use only K_iso instead of Kgpara2')
        else
          call warning('calc_heatcond_grad', &
              'Please use K_iso instead of Kgpara2')
          K_iso = Kgpara2
        endif
      endif
!
 99    return
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
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
        idiag_dtchi2=0.
        idiag_dtnewt=0.
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtchi2',idiag_dtchi2)
        call parse_name(iname,cname(iname),cform(iname),'dtnewt',idiag_dtnewt)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_dtchi2=',idiag_dtchi2
        write(3,*) 'i_dtnewt=',idiag_dtnewt
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
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
      case ('rtv')
        slices%yz =>rtv_yz
        slices%xz =>rtv_xz
        slices%xy =>rtv_xy
        slices%xy2=>rtv_xy2
        if (lwrite_slice_xy3) slices%xy3=>rtv_xy3
        if (lwrite_slice_xy4) slices%xy4=>rtv_xy4
        slices%ready=.true.
!
      case ('logQ')
        slices%yz =>logQ_yz
        slices%xz =>logQ_xz
        slices%xy =>logQ_xy
        slices%xy2=>logQ_xy2
        if (lwrite_slice_xy3) slices%xy3=>logQ_xy3
        if (lwrite_slice_xy4) slices%xy4=>logQ_xy4
        slices%ready=.true.
!
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
      use Diagnostics, only: max_mn_name
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension(nx) :: tmp
!
      if (lgranulation .and. (.not. lpencil_check_at_work)) then
        if (n.eq.n1.and.ipz.eq.0) then
          df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) - &
              tau_inv*(f(l1:l2,m,n,iux)-ux_local(:,m-nghost))
!
          df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) - &
              tau_inv*(f(l1:l2,m,n,iuy)-uy_local(:,m-nghost))
        endif
!
        tmp(:) = tau_inv
        if (lfirst.and.ldt) then
          if (ldiagnos.and.idiag_dtnewt/=0) then
            itype_name(idiag_dtnewt)=ilabel_max_dt
            call max_mn_name(tmp/cdts,idiag_dtnewt,l_dt=.true.)
          endif
          dt1_max=max(dt1_max,tmp/cdts)
        endif
      endif
!
      if (lmassflux) call force_solar_wind(df,p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!  computes hyper diffusion for non equidistant grid
!  using the IGNOREDX keyword.
!
!  17-feb-10/bing: coded
!
      use Deriv, only: der6
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: fdiff,tmp
!
      if (diffrho_hyper3/=0) then
        if (.not. ldensity_nolog) then
          call der6(f,ilnrho,fdiff,1,IGNOREDX=.true.)
          call der6(f,ilnrho,tmp,2,IGNOREDX=.true.)
          fdiff=fdiff + tmp
          call der6(f,ilnrho,tmp,3,IGNOREDX=.true.)
          fdiff=fdiff + tmp
          fdiff = diffrho_hyper3*fdiff
        else
          call fatal_error('special_calc_density', &
              'not yet implented for ldensity_nolog')
        endif
!
!        if (lfirst.and.ldt) diffus_diffrho3=diffus_diffrho3+diffrho_hyper3
!
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + fdiff
!
        if (headtt) print*,'special_calc_density: diffrho_hyper3=', &
            diffrho_hyper3
      endif
!
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy (or temperature) equation.
!
!  23-jun-08/bing: coded
!  17-feb-10/bing: added hyperdiffusion for non-equidistant grid
!
      use Deriv, only: der6,der4
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      real, dimension (nx) :: hc,tmp
!
      if (chi_hyper3/=0) then
        hc(:) = 0.
        call der6(f,ilnTT,tmp,1,IGNOREDX=.true.)
        hc = hc + tmp
        call der6(f,ilnTT,tmp,2,IGNOREDX=.true.)
        hc = hc + tmp
        call der6(f,ilnTT,tmp,3,IGNOREDX=.true.)
        hc = hc + tmp
        hc = chi_hyper3*hc
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + hc
!
!  due to ignoredx chi_hyperx has [1/s]
!
        if (lfirst.and.ldt) diffus_chi3=diffus_chi3  &
            + chi_hyper3
      endif
!
      if (chi_hyper2/=0) then
        hc(:) = 0.
        call der4(f,ilnTT,tmp,1,IGNOREDX=.true.)
        hc =  hc - chi_hyper2*tmp
        call der4(f,ilnTT,tmp,2,IGNOREDX=.true.)
        hc =  hc - chi_hyper2*tmp
        call der4(f,ilnTT,tmp,3,IGNOREDX=.true.)
        hc =  hc - chi_hyper2*tmp
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + hc
!
!  due to ignoredx chi_hyperx has [1/s]
!
        if (lfirst.and.ldt) diffus_chi3=diffus_chi3 &
            + chi_hyper2
      endif
!
      if (Kgpara/=0) call calc_heatcond_tensor(df,p,Kgpara,2.5)
      if (hcond1/=0) call calc_heatcond_constchi(df,p)
      if (hcond2/=0) call calc_heatcond_glnTT(df,p)
      if (hcond3/=0) call calc_heatcond_glnTT_iso(df,p)
      if (cool_RTV/=0) call calc_heat_cool_RTV(df,p)
      if (max (tdown, tdownr)/=0.0) call calc_heat_cool_newton(df,p)
      if (K_iso/=0) call calc_heatcond_grad(df,p)
      if (iheattype(1)/='nothing') call calc_artif_heating(df,p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mfarray) :: f
!
! Push vector potential back to initial vertical magnetic field values
      if (lfirst_proc_z .and. (bmdi /= 0.0)) then
        f(l1:l2,m1:m2,n1,iax)=f(l1:l2,m1:m2,n1,iax)*(1.-dt*bmdi) + &
            dt*bmdi * A_init_x
        f(l1:l2,m1:m2,n1,iay)=f(l1:l2,m1:m2,n1,iay)*(1.-dt*bmdi) + &
            dt*bmdi * A_init_y
!
        if (bmdi*dt > 1) call stop_it('special before boundary: bmdi*dt > 1 ')
      endif
!
! Read external velocity file. Has to be read before the granules are
! calculated.
      if (luse_ext_vel_field) call read_ext_vel_field()
!
! Compute photospheric granulation.
      if (lgranulation) then
        if (.not. lpencil_check_at_work) call uudriver(f)
      endif
!
      if (lmassflux) call get_wind_speed_offset(f)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine calc_heat_cool_newton(df,p)
!
!  newton cooling dT/dt = -1/tau * (T-T0)
!
      use Diagnostics, only: max_mn_name
      use EquationOfState, only: lnrho0
      use Sub, only: sine_step
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx) :: newton,newtonr,tmp_tau
      real, save :: lnrho_ref = -1.0
      integer :: pos, ref_pos
!
      if (headtt) &
          print *, 'special_calc_entropy: newton cooling active', tdown, tdownr
!
!
      tmp_tau = 0.0
!
      ! Correction of density profile
      if (tdownr /= 0.0) then
        ! Get reference density
        newtonr = exp (lnrho_init_z(n) - p%lnrho) - 1.0
        ! allpr is given in [Mm]
        tmp_tau = tdownr * exp (-allpr*unit_length*1e-6 * z(n))
        ! Add correction term to density
        df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + newtonr * tmp_tau
      endif
!
      ! Newton cooling of temperature profile
      if (tdown /= 0.0) then
        if (lnrho_ref == -1.0) then
          ! Get reference density
          ref_pos = 1
          do pos = 1, mz
            if (z(pos) <= 0.0) ref_pos = pos
          enddo
          lnrho_ref = lnrho_init_z(ref_pos)
          print *, 'calc_heat_cool_newton: reference density ', lnrho_ref
        endif
        ! Get reference temperature
        newton = exp (lnTT_init_z(n) - p%lnTT) - 1.0
!       tmp_tau = tdown * exp (-allp * (lnrho0 - p%lnrho))
        tmp_tau = tdown * sine_step (p%lnrho, lnrho_ref-allp, allp/2.0)
        ! Add newton cooling term to entropy
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + newton * tmp_tau
      endif
!
      if (lfirst.and.ldt) then
        if (ldiagnos.and.idiag_dtnewt/=0) then
          itype_name(idiag_dtnewt)=ilabel_max_dt
          call max_mn_name(tmp_tau/cdts,idiag_dtnewt,l_dt=.true.)
        endif
        dt1_max=max(dt1_max,tmp_tau/cdts)
      endif
!
    endsubroutine calc_heat_cool_newton
!***********************************************************************
    subroutine calc_heatcond_tensor(df,p,Kpara,expo)
!
!    anisotropic heat conduction with T^5/2
!    Div K T Grad ln T
!      =Grad(KT).Grad(lnT)+KT DivGrad(lnT)
!
      use Diagnostics,     only : max_mn_name
      use Sub,             only : dot2,dot,multsv,multmv,cubic_step
!--   use Io,              only : output_pencil
!AB: output_pencil is not currently used and breaks the auto-test
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: hhh,bunit,tmpv,gKp
      real, dimension (nx) :: tmpj,hhh2,quenchfactor
      real, dimension (nx) :: cosbgT,glnTT2,b2,bbb,b1,tmpk
      real, dimension (nx) :: chi_1,chi_2,rhs
      real :: Ksatb,Kpara,expo
      integer :: i,j,k
      type (pencil_case) :: p
!
!  calculate unit vector of bb
!
      call dot2(p%bb,bbb,PRECISE_SQRT=.true.)
      b1=1./max(tini,bbb)
      call multsv(b1,p%bb,bunit)
!
!  calculate H_i
!
      do i=1,3
        hhh(:,i)=0.
        do j=1,3
          tmpj(:)=0.
          do k=1,3
            tmpj(:)=tmpj(:)-2.*bunit(:,k)*p%bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+bunit(:,j)*(p%bij(:,i,j)+bunit(:,i)*tmpj(:))
        enddo
      enddo
      call multsv(b1,hhh,tmpv)
!
!  calculate abs(h) limiting
!
      call dot2(tmpv,hhh2,PRECISE_SQRT=.true.)
!
!  limit the length of h
!
      quenchfactor=1./max(1.,3.*hhh2*dxmax)
      call multsv(quenchfactor,tmpv,hhh)
!
      call dot(hhh,p%glnTT,rhs)
!
      chi_1 =  Kpara * p%rho1 * p%TT**expo* &
          cubic_step(real(t*unit_time),init_time,init_time)
!
      tmpv(:,:)=0.
      do i=1,3
        do j=1,3
          tmpv(:,i)=tmpv(:,i)+p%glnTT(:,j)*p%hlnTT(:,j,i)
        enddo
      enddo
!
      gKp = (expo+1.) * p%glnTT
!
      call dot2(p%glnTT,glnTT2)
!
      if (Ksat/=0.) then
        Ksatb = Ksat*7.28e7 /unit_velocity**3. * unit_temperature**1.5
!
        where (glnTT2 .le. tini)
          chi_2 =  0.
        elsewhere
          chi_2 =  Ksatb * sqrt(p%TT/max(tini,glnTT2))
        endwhere
!
        where (chi_1 .gt. chi_2)
          gKp(:,1)=p%glnrho(:,1) + 1.5*p%glnTT(:,1) - tmpv(:,1)/max(tini,glnTT2)
          gKp(:,2)=p%glnrho(:,2) + 1.5*p%glnTT(:,2) - tmpv(:,2)/max(tini,glnTT2)
          gKp(:,3)=p%glnrho(:,3) + 1.5*p%glnTT(:,3) - tmpv(:,3)/max(tini,glnTT2)
          chi_1 =  chi_2
        endwhere
      endif
!
      call dot(bunit,gKp,tmpj)
      call dot(bunit,p%glnTT,tmpk)
      rhs = rhs + tmpj*tmpk
!
      call multmv(p%hlnTT,bunit,tmpv)
      call dot(tmpv,bunit,tmpj)
      rhs = rhs + tmpj
!
      rhs = gamma*rhs*chi_1
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + rhs
!
!  for timestep extension multiply with the
!  cosine between grad T and bunit
!
      call dot(p%bb,p%glnTT,cosbgT)
      call dot2(p%bb,b2)
!
      where (glnTT2*b2.le.tini)
        cosbgT=0.
      elsewhere
        cosbgT=cosbgT/sqrt(glnTT2*b2)
      endwhere
!
      if (lfirst.and.ldt) then
        chi_1=abs(cosbgT)*chi_1
        diffus_chi=diffus_chi+gamma*chi_1*dxyz_2
        if (ldiagnos.and.idiag_dtchi2/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_tensor
!***********************************************************************
    subroutine calc_heatcond_grad(df,p)
!
! additional heat conduction where the heat flux is
! is proportional to \rho abs(gradT)
! L= Div(rho |Grad(T)| Grad(T))
!
      use Diagnostics,     only : max_mn_name
      use Sub,             only : dot,dot2
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: tmpv
      real, dimension (nx) :: tmpj,tmpi
      real, dimension (nx) :: rhs,g2,chix
      integer :: i,j
      type (pencil_case) :: p
!
      call dot2(p%glnTT,tmpi)
!
      tmpv(:,:)=0.
      do i=1,3
         do j=1,3
            tmpv(:,i)=tmpv(:,i)+p%glnTT(:,j)*p%hlnTT(:,j,i)
         enddo
      enddo
      call dot(tmpv,p%glnTT,tmpj)
!
      call dot(p%glnrho,p%glnTT,g2)
!
      rhs=p%TT*(tmpi*(p%del2lnTT+2.*tmpi + g2)+tmpj)/max(tini,sqrt(tmpi))
!
!      if (itsub .eq. 3 .and. ip .eq. 118) &
!          call output_pencil(trim(directory)//'/tensor3.dat',rhs,1)
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+ K_iso * rhs
!
      if (lfirst.and.ldt) then
        chix=K_iso*p%TT*sqrt(tmpi)
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        if (ldiagnos.and.idiag_dtchi2/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_grad
!***********************************************************************
    subroutine calc_heatcond_constchi(df,p)
!
! L = Div( K rho b*(b*Grad(T))
!
      use Diagnostics,     only : max_mn_name
      use Sub,             only : dot2,dot,multsv,multmv,cubic_step
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: bunit,hhh,tmpv
      real, dimension (nx) :: hhh2,quenchfactor
      real, dimension (nx) :: abs_b,b1
      real, dimension (nx) :: rhs,tmp,tmpi,tmpj,chix
      integer :: i,j,k
!
      intent(in) :: p
      intent(out) :: df
!
      if (headtt) print*,'special/calc_heatcond_chiconst',hcond1
!
!  Define chi= K_0/rho
!
!  calculate unit vector of bb
!
      call dot2(p%bb,abs_b,PRECISE_SQRT=.true.)
      b1=1./max(tini,abs_b)
      call multsv(b1,p%bb,bunit)
!
!  calculate first H_i
!
      do i=1,3
        hhh(:,i)=0.
        do j=1,3
          tmpj(:)=0.
          do k=1,3
            tmpj(:)=tmpj(:)-2.*bunit(:,k)*p%bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+bunit(:,j)*(p%bij(:,i,j)+bunit(:,i)*tmpj(:))
        enddo
      enddo
      call multsv(b1,hhh,tmpv)
!
!  calculate abs(h) for limiting H vector
!
      call dot2(tmpv,hhh2,PRECISE_SQRT=.true.)
!
!  limit the length of H
!
      quenchfactor=1./max(1.,3.*hhh2*dxmax)
      call multsv(quenchfactor,tmpv,hhh)
!
!  dot H with Grad lnTT
!
      call dot(hhh,p%glnTT,tmp)
!
!  dot Hessian matrix of lnTT with bi*bj, and add into tmp
!
      call multmv(p%hlnTT,bunit,tmpv)
      call dot(tmpv,bunit,tmpj)
      tmp = tmp+tmpj
!
!  calculate (Grad lnTT * bunit)^2 needed for lnecr form; also add into tmp
!
      call dot(p%glnTT,bunit,tmpi)
!
      call dot(p%glnrho,bunit,tmpj)
      tmp=tmp+(tmpj+tmpi)*tmpi
!
!  calculate rhs
!
      chix = hcond1*cubic_step(real(t*unit_time),init_time,init_time)
!
      rhs = gamma*chix*tmp
!
      if (.not.(ipz.eq.nprocz-1.and.n.ge.n2-3)) &
          df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+rhs
!
!      if (itsub .eq. 3 .and. ip .eq. 118) &
!          call output_pencil(trim(directory)//'/tensor2.dat',rhs,1)
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chix*dxyz_2
        advec_cs2=max(advec_cs2,maxval(chix*dxyz_2))
        if (ldiagnos.and.idiag_dtchi2/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_constchi
!***********************************************************************
    subroutine calc_heatcond_glnTT(df,p)
!
!  L = Div( Grad(lnT)^2  rho b*(b*Grad(T)))
!  => flux = T  Grad(lnT)^2  rho
!    gflux = flux * (glnT + glnrho +Grad( Grad(lnT)^2))
!  => chi =  Grad(lnT)^2
!
      use Diagnostics, only: max_mn_name
      use Sub, only: dot2,dot,multsv,multmv,cubic_step
      use EquationOfState, only: gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: glnT2
      real, dimension (nx) :: tmp,rhs,chi
      real, dimension (nx,3) :: bunit,hhh,tmpv,gflux
      real, dimension (nx) :: hhh2,quenchfactor
      real, dimension (nx) :: abs_b,b1
      real, dimension (nx) :: tmpj
      integer :: i,j,k
!
      intent(in) :: p
      intent(out) :: df
!
!  calculate unit vector of bb
!
      call dot2(p%bb,abs_b,PRECISE_SQRT=.true.)
      b1=1./max(tini,abs_b)
      call multsv(b1,p%bb,bunit)
!
!  calculate first H_i
!
      do i=1,3
        hhh(:,i)=0.
        do j=1,3
          tmpj(:)=0.
          do k=1,3
            tmpj(:)=tmpj(:)-2.*bunit(:,k)*p%bij(:,k,j)
          enddo
          hhh(:,i)=hhh(:,i)+bunit(:,j)*(p%bij(:,i,j)+bunit(:,i)*tmpj(:))
        enddo
      enddo
      call multsv(b1,hhh,tmpv)
!
!  calculate abs(h) for limiting H vector
!
      call dot2(tmpv,hhh2,PRECISE_SQRT=.true.)
!
!  limit the length of H
!
      quenchfactor=1./max(1.,3.*hhh2*dxmax)
      call multsv(quenchfactor,tmpv,hhh)
!
!  dot H with Grad lnTT
!
      call dot(hhh,p%glnTT,tmp)
!
!  dot Hessian matrix of lnTT with bi*bj, and add into tmp
!
      call multmv(p%hlnTT,bunit,tmpv)
      call dot(tmpv,bunit,tmpj)
      tmp = tmp+tmpj
!
      call dot2(p%glnTT,glnT2)
!
      tmpv = p%glnTT+p%glnrho
!
      call multsv(glnT2,tmpv,gflux)
!
      do i=1,3
        tmpv(:,i) = 2.*(&
            p%glnTT(:,1)*p%hlnTT(:,i,1) + &
            p%glnTT(:,2)*p%hlnTT(:,i,2) + &
            p%glnTT(:,3)*p%hlnTT(:,i,3) )
      enddo
!
      gflux  = gflux +tmpv
!
      call dot(gflux,bunit,rhs)
      call dot(p%glnTT,bunit,tmpj)
      rhs = rhs*tmpj
!
      chi = glnT2*hcond2
!
      rhs = hcond2*(rhs + glnT2*tmp)
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + &
          rhs*gamma*cubic_step(real(t*unit_time),init_time,init_time)
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chi*dxyz_2
        if (ldiagnos.and.idiag_dtchi2/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_glnTT
!***********************************************************************
    subroutine calc_heatcond_glnTT_iso(df,p)
!
!  L = Div( Grad(lnT)^2 Grad(T))
!
      use Diagnostics,     only : max_mn_name
      use Sub,             only : dot2,dot,multsv,multmv,cubic_step
      use EquationOfState, only : gamma
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: tmpv
      real, dimension (nx) :: glnT2,glnT_glnr
      real, dimension (nx) :: tmp,rhs,chi
      integer :: i
!
      intent(in) :: p
      intent(out) :: df
!
      call dot2(p%glnTT,glnT2)
      call dot(p%glnTT,p%glnrho,glnT_glnr)
!
      do i=1,3
        tmpv(:,i) = p%glnTT(:,1)*p%hlnTT(:,1,i) + &
                    p%glnTT(:,2)*p%hlnTT(:,2,i) + &
                    p%glnTT(:,3)*p%hlnTT(:,3,i)
      enddo
      call dot(p%glnTT,tmpv,tmp)
!
      chi = glnT2*hcond3
!
      rhs = 2*tmp+glnT2*(glnT2+p%del2lnTT+glnT_glnr)
!
      df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT) + rhs*gamma*hcond3
!
      if (lfirst.and.ldt) then
        diffus_chi=diffus_chi+gamma*chi*dxyz_2
        if (ldiagnos.and.idiag_dtchi2/=0) then
          call max_mn_name(diffus_chi/cdtv,idiag_dtchi2,l_dt=.true.)
        endif
      endif
!
    endsubroutine calc_heatcond_glnTT_iso
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
      use Sub,             only: cubic_step
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: lnQ,rtv_cool,lnTT_SI,lnneni
      real :: unit_lnQ
      type (pencil_case) :: p
!
      unit_lnQ=3*alog(real(unit_velocity))+&
          5*alog(real(unit_length))+alog(real(unit_density))
      lnTT_SI = p%lnTT + alog(real(unit_temperature))
!
!     calculate ln(ne*ni) :
!          ln(ne*ni) = ln( 1.17*rho^2/(1.34*mp)^2)
!     lnneni = 2*p%lnrho + alog(1.17) - 2*alog(1.34)-2.*alog(real(m_p))
!
      lnneni = 2.*(p%lnrho+61.4412 +alog(real(unit_mass)))
!
      lnQ   = get_lnQ(lnTT_SI)
!
      rtv_cool = lnQ-unit_lnQ+lnneni-p%lnTT-p%lnrho
      rtv_cool = gamma*p%cp1*exp(rtv_cool)
!
      rtv_cool = rtv_cool*cool_RTV *cubic_step(real(t*unit_time),init_time,init_time)
!     for adjusting by setting cool_RTV in run.in
!
      rtv_cool=rtv_cool &
          *(1.-cubic_step(p%lnrho,-12.-alog(real(unit_density)),3.))
!
! slices
      rtv_yz(m-m1+1,n-n1+1)=rtv_cool(ix_loc-l1+1)
      if (m==iy_loc)  rtv_xz(:,n-n1+1)= rtv_cool
      if (n==iz_loc)  rtv_xy(:,m-m1+1)= rtv_cool
      if (n==iz2_loc) rtv_xy2(:,m-m1+1)= rtv_cool
      if (n==iz3_loc) rtv_xy3(:,m-m1+1)= rtv_cool
      if (n==iz4_loc) rtv_xy4(:,m-m1+1)= rtv_cool
!
!     add to temperature equation
!
      if (ltemperature) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT)-rtv_cool
      else
        if (lentropy) &
            call stop_it('solar_corona: calc_heat_cool:lentropy=not implented')
      endif
!
      if (lfirst.and.ldt) then
        if (ldiagnos.and.idiag_dtnewt/=0) then
          itype_name(idiag_dtnewt)=ilabel_max_dt
          call max_mn_name(rtv_cool/cdts,idiag_dtnewt,l_dt=.true.)
        endif
        dt1_max=max(dt1_max,rtv_cool/cdts)
      endif
!
      logQ_yz(m-m1+1,n-n1+1)=lnQ(ix_loc-l1+1)*0.43429448
      if (m==iy_loc)  logQ_xz(:,n-n1+1)= lnQ*0.43429448
      if (n==iz_loc)  logQ_xy(:,m-m1+1)= lnQ*0.43429448
      if (n==iz2_loc) logQ_xy2(:,m-m1+1)= lnQ*0.43429448
      if (n==iz3_loc) logQ_xy3(:,m-m1+1)= lnQ*0.43429448
      if (n==iz4_loc) logQ_xy4(:,m-m1+1)= lnQ*0.43429448
!
    endsubroutine calc_heat_cool_RTV
!***********************************************************************
    function get_lnQ(lnTT)
!
!  input: lnTT in SI units
!  output: lnP  [p]=W/s * m^3
!
      real, parameter, dimension (37) :: intlnT = (/ &
          8.74982, 8.86495, 8.98008, 9.09521, 9.21034, 9.44060, 9.67086 &
          , 9.90112, 10.1314, 10.2465, 10.3616, 10.5919, 10.8221, 11.0524 &
          , 11.2827, 11.5129, 11.7432, 11.9734, 12.2037, 12.4340, 12.6642 &
          , 12.8945, 13.1247, 13.3550, 13.5853, 13.8155, 14.0458, 14.2760 &
          , 14.5063, 14.6214, 14.7365, 14.8517, 14.9668, 15.1971, 15.4273 &
          ,  15.6576,  69.0776 /)
      real, parameter, dimension (37) :: intlnQ = (/ &
          -93.9455, -91.1824, -88.5728, -86.1167, -83.8141, -81.6650 &
          , -80.5905, -80.0532, -80.1837, -80.2067, -80.1837, -79.9765 &
          , -79.6694, -79.2857, -79.0938, -79.1322, -79.4776, -79.4776 &
          , -79.3471, -79.2934, -79.5159, -79.6618, -79.4776, -79.3778 &
          , -79.4008, -79.5159, -79.7462, -80.1990, -80.9052, -81.3196 &
          , -81.9874, -82.2023, -82.5093, -82.5477, -82.4172, -82.2637 &
          , -0.66650 /)
!
      real, parameter, dimension (16) :: intlnT1 = (/ &
          8.98008,    9.44060,    9.90112,    10.3616,    10.8221,    11.2827 &
         ,11.5129,    11.8583,    12.4340,    12.8945,    13.3550,    13.8155 &
         ,14.2760,    14.9668,    15.8878,    18.4207 /)
      real, parameter, dimension (16) :: intlnQ1 = (/ &
          -83.9292,   -81.2275,   -80.0532,   -80.1837,   -79.6694,   -79.0938 &
         ,-79.1322,   -79.4776,   -79.2934,   -79.6618,   -79.3778,   -79.5159 &
         ,-80.1990,   -82.5093,   -82.1793,   -78.6717 /)
!
      real, dimension(9) :: pars=(/2.12040e+00,3.88284e-01,2.02889e+00, &
          3.35665e-01,6.34343e-01,1.94052e-01,2.54536e+00,7.28306e-01, &
          -2.40088e+01/)
!
      real, dimension (nx) :: lnTT,get_lnQ
      real, dimension (nx) :: slope,ordinate
      real, dimension (nx) :: logT_SI,logQ
      integer :: i
!
!  select type for cooling fxunction
!  1: 10 points interpolation
!  2: 37 points interpolation
!  3: four gaussian fit
!  4: several fits
!
      get_lnQ=-1000.
!
      select case(cool_type)
      case(1)
        do i=1,15
          where(lnTT .ge. intlnT1(i) .and. lnTT .lt. intlnT1(i+1))
            slope=(intlnQ1(i+1)-intlnQ1(i))/(intlnT1(i+1)-intlnT1(i))
            ordinate = intlnQ1(i) - slope*intlnT1(i)
            get_lnQ = slope*lnTT + ordinate
          endwhere
        enddo
!
      case(2)
        do i=1,36
          where(lnTT .ge. intlnT(i) .and. lnTT .lt. intlnT(i+1))
            slope=(intlnQ(i+1)-intlnQ(i))/(intlnT(i+1)-intlnT(i))
            ordinate = intlnQ(i) - slope*intlnT(i)
            get_lnQ = slope*lnTT + ordinate
          endwhere
        enddo
!
      case(3)
        call fatal_error('get_lnQ','this invokes epx() to often')
        lnTT  = lnTT*alog10(exp(1.))
        get_lnQ  = -1000.
        !
        get_lnQ = pars(1)*exp(-(lnTT-4.3)**2/pars(2)**2)  &
            +pars(3)*exp(-(lnTT-4.9)**2/pars(4)**2)  &
            +pars(5)*exp(-(lnTT-5.35)**2/pars(6)**2) &
            +pars(7)*exp(-(lnTT-5.85)**2/pars(8)**2) &
            +pars(9)
        !
        get_lnQ = get_lnQ * (20.*(-tanh((lnTT-3.7)*10.))+21)
        get_lnQ = get_lnQ+(tanh((lnTT-6.9)*3.1)/2.+0.5)*3.
        !
        get_lnQ = (get_lnQ +19.-32)*alog(10.)
!
      case(4)
        logT_SI = lnTT*0.43429448+alog(real(unit_temperature))
        where(logT_SI<=3.928)
          logQ =  -7.155e1 + 9*logT_SI
        elsewhere(logT_SI>3.93.and.logT_SI<=4.55)
          logQ = +4.418916e+04 &
              -5.157164e+04 * logT_SI &
              +2.397242e+04 * logT_SI**2 &
              -5.553551e+03 * logT_SI**3 &
              +6.413137e+02 * logT_SI**4 &
              -2.953721e+01 * logT_SI**5
        elsewhere(logT_SI>4.55.and.logT_SI<=5.09)
          logQ = +8.536577e+02 &
              -5.697253e+02 * logT_SI &
              +1.214799e+02 * logT_SI**2 &
              -8.611106e+00 * logT_SI**3
        elsewhere(logT_SI>5.09.and.logT_SI<=5.63)
          logQ = +1.320434e+04 &
              -7.653183e+03 * logT_SI &
              +1.096594e+03 * logT_SI**2 &
              +1.241795e+02 * logT_SI**3 &
              -4.224446e+01 * logT_SI**4 &
              +2.717678e+00 * logT_SI**5
        elsewhere(logT_SI>5.63.and.logT_SI<=6.48)
          logQ = -2.191224e+04 &
              +1.976923e+04 * logT_SI &
              -7.097135e+03 * logT_SI**2 &
              +1.265907e+03 * logT_SI**3 &
              -1.122293e+02 * logT_SI**4 &
              +3.957364e+00 * logT_SI**5
        elsewhere(logT_SI>6.48.and.logT_SI<=6.62)
          logQ = +9.932921e+03 &
              -4.519940e+03 * logT_SI &
              +6.830451e+02 * logT_SI**2 &
              -3.440259e+01 * logT_SI**3
        elsewhere(logT_SI>6.62)
          logQ = -3.991870e+01 + 6.169390e-01 * logT_SI
        endwhere
        get_lnQ = (logQ+19.-32)*2.30259
      case default
        call fatal_error('get_lnQ','wrong type')
      endselect
!
    endfunction get_lnQ
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
      use Sub, only: cubic_step
!
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: heatinput,heat_flux
      real, dimension (nx) :: x_Mm,heat_nano,rhs
      real, dimension (nx) :: heat_event,heat_event1D
      integer, dimension(mseed) :: global_rstate
      real :: z_Mm,heat_unit
      real :: nano_sigma_t,nano_time,nano_start,nano_sigma_z
      real :: nano_flare_energy,nano_pos_x,nano_pos_z,nano_pos_y
      real :: nano_amplitude
      real, dimension(2) :: event_pos
      type (pencil_case) :: p
      integer :: i
!
      heat_unit= unit_density*unit_velocity**3/unit_length
      x_Mm = x(l1:l2)*unit_length*1e-6
      z_Mm = z(n)*unit_length*1e-6
!
      heatinput = 0.
      heat_flux = 0.
!
      do i=1,3
        if (headtt) print*,'iheattype:',iheattype(i)
        select case(iheattype(i))
        case ('nothing')
          !
        case ('exp')
          ! heat_par_exp(1) should be 530 W/m^3 (amplitude)
          ! heat_par_exp(2) should be 0.3 Mm (scale height)
          !
          heatinput=heatinput + &
              heat_par_exp(1)*exp(-z_Mm/heat_par_exp(2))/heat_unit
!
          heat_flux=heat_flux +  heat_par_exp(1)*heat_par_exp(2)*1e6* &
              (1.-exp(-lz*unit_length*1e-6/heat_par_exp(2)))
!
          if (headtt) print*,'Flux of exp heating: ',heat_flux
!
        case ('exp2')
          ! A second exponential function
          ! For Sven et al. 2010 set:
          ! heat_par_exp= (1e3 , 0.2 )
          ! heat_par_exp2= (1e-4 , 10.)
          !
          heatinput=heatinput + &
              heat_par_exp2(1)*exp(-z_Mm/heat_par_exp2(2))/heat_unit
!
          heat_flux=heat_flux + heat_par_exp2(1)*heat_par_exp2(2)*1e-6* &
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
          heatinput=heatinput + &
              heat_par_gauss(3)*exp(-((z_Mm-heat_par_gauss(1))**2/ &
              (2*heat_par_gauss(2)**2)))/heat_unit
!
        case ('nanof')
          ! simulate nanoflare heating =)
          ! height dependend amplitude und duration
          ! idea: call random numbers , make condition when to flare, when
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
          if (nano_start .eq. 0.) then
            ! If 0 roll to see if a nanoflare occurs
            call random_number_wrapper(nano_start)
            !
            if (nano_start .gt. 0.95 ) then
              ! 5% chance for a nanoflare to occur, then get the location.
              call normal_deviate(nano_pos_z)
              nano_pos_z=nano_pos_z*lz
              call random_number_wrapper(nano_pos_y)
              call random_number_wrapper(nano_pos_x)
              nano_time=60.
            else
              ! else no nanoflare, reset nano_start to 0. for the next
              ! timestep
              nano_start=0.
            endif
          endif
          !
          if (nano_start .ne. 0.) then
            ! if nano_start is not 0. then there is a nanoflare!
            ! 2nd assumption , nanoflare takes 60 seconds =)
            nano_flare_energy=10.d17 !joules
            nano_sigma_z=0.5
            nano_sigma_t=2.5
!
            nano_amplitude=nano_flare_energy/(pi/2*nano_sigma_t*nano_sigma_z*1.d6 )
!
            heat_nano=nano_amplitude*exp(-((nano_time-5.))**2/( 2*2.**2))* &
                exp(-((z_Mm-nano_pos_z)**2/ (2*0.5**2)))
            nano_time=nano_time-dt*unit_time
!
            heatinput=heatinput + heat_nano/heat_unit
!
            if (nano_time .le. 0.) nano_start=0.
          end if
          !
          !SAVE NANO_SEED
          !RESTORE GLOBAL SEED
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
          if (t*unit_time .gt. 150. .AND. t*unit_time .lt. 1000.) then
            event_pos(1)=7.5
            event_pos(2)=15.
            heat_event=10.*exp(-((250.-t*unit_time))**2/(2*(20.*unit_time)**2))* &
                exp(-((x_Mm-event_pos(1))**2/ (2*0.2**2))) * &
                exp(-((z_Mm-event_pos(2))**2/ (2*0.2**2)))
            heatinput=heatinput + heat_event/heat_unit
          endif
!
        case ('event1D')
          ! one small point heating event (gaussian to prevent too strong gradients)
          ! one point in time, one point in space!
          if (t*unit_time .gt. 300. .AND. t*unit_time .lt. 10000.) then
            if (t*unit_time .gt. 300. .AND. t*unit_time .lt. 301.) &
                print*,'EVENTTTT!!!!!'
            event_pos(1)=10.
            heat_event1D=10.*exp(-((400.-t))**2/( 2*50.**2))* &
                exp(-((x_Mm-event_pos(1))**2/ (2*0.2**2)))
            heatinput=heatinput + heat_event1D/heat_unit
          endif
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
      rhs = p%TT1*p%rho1*gamma*p%cp1*heatinput* &
          cubic_step(real(t*unit_time),init_time,init_time)
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) + rhs
!
      if (lfirst.and.ldt) then
        if (ldiagnos.and.idiag_dtnewt/=0) then
          itype_name(idiag_dtnewt)=ilabel_max_dt
          call max_mn_name(rhs/cdts,idiag_dtnewt,l_dt=.true.)
        endif
        dt1_max=max(dt1_max,rhs/cdts)
      endif
!
    endsubroutine calc_artif_heating
!***********************************************************************
    subroutine setdrparams()
!
      Use Mpicomm, only: stop_it_if_any
!
      integer :: alloc_err_sum
!
! Every granule has 6 values associated with it: data(1-6).
! These contain,  x-position, y-position,
!    current amplitude, amplitude at t=t_0, t_0, and life_time.
!
! Gives intergranular area / (granular+intergranular area)
      ig=0.3
!
! Gives average radius of granule + intergranular lane
! (no smaller than 6 grid points across)
!     here radius of granules is 0.8 Mm or bigger (3 times dx)
!
      if (unit_system.eq.'SI') then
        granr=max(0.8*1.e6/unit_length,3*dx,3*dy)
      elseif  (unit_system.eq.'cgs') then
        granr=max(0.8*1.e8/unit_length,3*dx,3*dy)
      endif
!
! Fractional difference in granule power
      pd=0.15
!
! Gives exponential power of evolvement. Higher power faster growth/decay.
      pow=2
!
! Fractional distance, where after emergence, no granule can emerge
! whithin this radius.(In order to not 'overproduce' granules).
! This limit is unset when granule reaches t_0.
      avoid=0.8
!
! Lifetime of granule
! Now also resolution dependent(5 min for granular scale)
!
      life_t=(60.*5./unit_time)
      !*(granr/(0.8*1e8/u_l))**2
      !  removed since the life time was about 20 min !
!
      dxdy2=dx**2+dy**2
!
! Typical central velocity of granule(1.5e5 cm/s=1.5km/s)
! Has to be multiplied by the smallest distance, since velocity=ampl/dist
! should now also be dependant on smallest granluar scale resolvable.
!
      if (unit_system.eq.'SI') then
        ampl=sqrt(dxdy2)/granr*0.28e4/unit_velocity
      elseif (unit_system.eq.'cgs') then
        ampl=sqrt(dxdy2)/granr*0.28e6/unit_velocity
      endif
!
! fraction of current amplitude to maximum amplitude to the beginning
! and when the granule disapears
      thresh=0.78
!
      xrange=min(nint(1.5*granr*(1+ig)/dx),nint(nxgrid/2.0)-1)
      yrange=min(nint(1.5*granr*(1+ig)/dy),nint(nygrid/2.0)-1)
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
      if (.not.lreloading) then
!
        points_rstate(:)=0.
!
        isnap = nint (t/dsnap)
        tsnap_uu = (isnap+1) * dsnap
!
      endif
!
      alloc_err=0
      if (.not. allocated (Ux)) allocate(Ux(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = abs(alloc_err)
      if (.not. allocated (Uy)) allocate(Uy(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (.not. allocated(w))  allocate (w(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (.not. allocated(vx))  allocate (vx(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (.not. allocated(vy))  allocate (vy(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (.not. allocated(avoidarr)) &
          allocate(avoidarr(nxgrid,nygrid),stat=alloc_err)
      alloc_err_sum = alloc_err_sum + abs(alloc_err)
      if (alloc_err_sum > 0) call stop_it_if_any(.true., &
          'setdrparams: Could not allocate memory for the driver')
!
    endsubroutine setdrparams
!***********************************************************************
    subroutine uudriver(f)
!
! This routine replaces the external computing of granular velocity
! pattern initially written by B. Gudiksen (2004)
!
! It is invoked by setting lgranulation=T in run.in
! additional parameters are
!         Bavoid =0.01 : the magn. field strenght in Tesla at which
!                        no granule is allowed
!         nvor = 5.    : the strength by which the vorticity is
!                        enhanced
!
!  11-may-10/bing: coded
!
      use EquationOfState, only: gamma_inv,get_cp1,gamma_m1,lnrho0,cs20
      use General, only: random_seed_wrapper
      use Mpicomm, only: mpisend_real, mpirecv_real
      use Sub, only: cubic_step
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(:,:), allocatable :: uu_buffer
      integer :: i,j,ipt,main_proc
      real, dimension(nx,ny) :: pp_tmp,BB2_local,beta,quench
      real :: cp1=1.,dA
      integer, dimension(2) :: dims=(/nx,ny/)
      integer, dimension(mseed) :: global_rstate
      real, save :: next_time = 0.0
!
      if (.not.allocated(uu_buffer)) then
        allocate(uu_buffer(nxgrid,nygrid),stat=alloc_err)
        if (alloc_err>0) call fatal_error('uudriver', &
            'could not allocate memory for uu_buffer')
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
! Get magnetic field energy for footpoint quenching.
! The processors from the ipz=0 plane have to take part, too.
      if (lgran_proc .or. lfirst_proc_z) then
        if (lmagnetic) call set_B2(f,BB2_local)
!
! Set sum(abs(Bz)) to  a given flux.
        if (Bz_flux/=0) then
          if (nxgrid/=1.and.nygrid/=1) then
            dA=dx*dy*unit_length**2
          elseif (nygrid==1) then
            dA=dx*unit_length
          elseif (nxgrid==1) then
            dA=dy*unit_length
          endif
          f(l1:l2,m1:m2,n1,iax:iaz) = f(l1:l2,m1:m2,n1,iax:iaz) * &
              Bz_flux/(Bzflux*dA*unit_magnetic)
        endif
      endif
!
! Compute granular velocities. We use three levels.
!
!
! Either root processor or three procs with ipz>0 compute
! velocities for different levels in driver3().
      if (lgran_proc) call multi_drive3()
!
! In the parallel case, one proc has to sum up the levels.
      if (lgran_parallel .and. lgran_proc) then
        if ((iproc > nprocxy) .and. (iproc < nprocxy+nglevel)) then
          call mpisend_real(Ux,(/nxgrid,nygrid/),nprocxy,iproc)
          call mpisend_real(Uy,(/nxgrid,nygrid/),nprocxy,iproc)
        elseif (nglevel > 1) then
          do i=1, nglevel-1
            call mpirecv_real(uu_buffer,(/nxgrid,nygrid/),nprocxy+i,nprocxy+i)
            Ux = Ux + uu_buffer
            call mpirecv_real(uu_buffer,(/nxgrid,nygrid/),nprocxy+i,nprocxy+i)
            Uy = Uy + uu_buffer
          enddo
        endif
      endif
!
! Increase vorticity and normalize to given vrms.
      if (lgran_proc .and. lfirst_proc_xy) call enhance_vorticity()
!
! Distribute results, first select the proc which collected the data.
!
      if (lgran_parallel) then
        main_proc = nprocxy
      else
        main_proc = 0
      endif
!
! Then distribute.
      if (iproc==main_proc) then
        do i=0,nprocx-1
          do j=0,nprocy-1
            ipt = i+nprocx*j
            if (ipt == main_proc) then
              ux_local = Ux(i*nx+1:i*nx+nx,j*ny+1:j*ny+ny)
              uy_local = Uy(i*nx+1:i*nx+nx,j*ny+1:j*ny+ny)
            else
              call mpisend_real(Ux(i*nx+1:i*nx+nx,j*ny+1:j*ny+ny),dims,ipt,312+ipt)
              call mpisend_real(Uy(i*nx+1:i*nx+nx,j*ny+1:j*ny+ny),dims,ipt,313+ipt)
            endif
          enddo
        enddo
      elseif (lfirst_proc_z) then
        call mpirecv_real(ux_local,dims,main_proc,312+iproc)
        call mpirecv_real(uy_local,dims,main_proc,313+iproc)
      endif
!
! for footpoint quenching compute pressure
!
      if (leos) call get_cp1(cp1)
!
      if (lquench.and.lfirst_proc_z) then
        if (ltemperature.and..not.ltemperature_nolog) then
          if (ldensity_nolog) then
            call fatal_error('solar_corona', &
                'uudriver only implemented for ltemperature=true')
          else
            pp_tmp =gamma_m1*gamma_inv/cp1 * &
                exp(f(l1:l2,m1:m2,irefz,ilnrho)+f(l1:l2,m1:m2,irefz,ilnrho))
          endif
        else
          if (headt.and.lroot) call warning('solar_corona', &
              'uudriver only implemented for ltemperature=true')
          pp_tmp=gamma_inv*cs20*exp(lnrho0)
        endif
!
        beta =  pp_tmp/max(tini,BB2_local)*2*mu0
!
!  quench velocities to one percent of the granule velocities
        do i=1,ny
          quench(:,i) = cubic_step(beta(:,i),q0,qw)*(1.-dq)+dq
        enddo
!
        ux_local = ux_local * quench
        uy_local = uy_local * quench
      endif
!
      if (luse_ext_vel_field) then
        ux_local = ux_local + ux_ext_local
        uy_local = uy_local + uy_ext_local
      endif
!
      f(l1:l2,m1:m2,irefz,iuz) = 0.
!
! restore global seed and save seed list of the granulation
      call random_seed_wrapper(GET=points_rstate)
      call random_seed_wrapper(PUT=global_rstate)
!
      if (allocated(uu_buffer)) deallocate(uu_buffer)
!
    endsubroutine uudriver
!***********************************************************************
  subroutine multi_drive3
!
    real, parameter :: ldif=2.0
    real, dimension(max_gran_levels), save :: amplarr,granrarr,life_tarr
    integer, dimension(max_gran_levels), save :: xrangearr,yrangearr
    integer :: k
    logical, save :: first_call=.true.
!
    if (first_call) then
!
      granrarr(1)=granr
      granrarr(2)=granr*ldif
      granrarr(3)=granr*ldif*ldif
!
      amplarr(1)=ampl
      amplarr(2)=ampl/ldif
      amplarr(3)=ampl/(ldif*ldif)
!
      life_tarr(1)=life_t
      life_tarr(2)=ldif**2*life_t
      life_tarr(3)=ldif**4*life_t
!
      xrangearr(1)=xrange
      xrangearr(2)=min(nint(ldif*xrange),nint(nxgrid/2.-1.))
      xrangearr(3)=min(nint(ldif*ldif*xrange),nint(nxgrid/2.-1.))
      yrangearr(1)=yrange
      yrangearr(2)=min(nint(ldif*yrange),nint(nygrid/2-1.))
      yrangearr(3)=min(nint(ldif*ldif*yrange),nint(nygrid/2-1.))
!
      first_call = .false.
!
    endif
 !
    ! Initialize velocity field
    Ux = 0.0
    Uy = 0.0
!
    do k=1,nglevel
      select case (k)
      case (1)
        first => firstlev
      case (2)
        first => secondlev
      case (3)
        first => thirdlev
      end select
!
      ampl=amplarr(k)
      granr=granrarr(k)
      life_t=life_tarr(k)
      xrange=xrangearr(k)
      yrange=yrangearr(k)
!
      if ((iproc==nprocxy-1+k) .or. lroot) call drive3(k)
!
      select case (k)
      case (1)
        firstlev => first
      case (2)
        secondlev => first
      case (3)
        thirdlev => first
      end select
!
    enddo
!
    endsubroutine multi_drive3
!***********************************************************************
    subroutine drive3(level)
!
      use Syscalls, only: file_exists
!
      integer, intent(in) :: level
      logical :: lstop=.false.
!
      call resetarr
      if (Bavoid > 0.) call fill_B_avoidarr
!
      if (.not.associated(first)) then
        call rdpoints(level)
        if (.not.associated(first)) then
          call driveinit
          call wrpoints(level,0)
          call wrpoints(level)
        endif
      else
        call updatepoints
        current => first
        do while (associated (current))
          call drawupdate
          current => current%next
        enddo
        do while (minval(avoidarr) == 0)
          call add_point
          call make_newpoint
          call drawupdate
        enddo
      endif
!
      Ux = Ux + vx
      Uy = Uy + vy
!
! Move granule centers according to external velocity
! field. Needed to be done for each level
!
      if (luse_ext_vel_field) call evolve_granules()
!
! Save granules to file.
!
      if (t >= tsnap_uu) then
        call wrpoints(level,isnap)
        if (level == nglevel) then
          tsnap_uu = tsnap_uu + dsnap
          isnap  = isnap + 1
        endif
      endif
      if (itsub == 3) &
          lstop = file_exists('STOP')
      if (lstop .or. (t>=tmax) .or. (it == nt) .or. (dt < dtmin) .or. &
          (mod(it,isave) == 0)) call wrpoints(level)
!
    endsubroutine drive3
!***********************************************************************
    subroutine resetarr
!
! Reset arrays at the beginning of each call to the levels.
!
! 12-aug-10/bing: coded
!
      w(:,:) = 0.0
      vx(:,:) = 0.0
      vy(:,:) = 0.0
      avoidarr(:,:) = 0
!
    endsubroutine resetarr
!***********************************************************************
    subroutine rdpoints(level)
!
! Read in points from file if existing.
!
! 12-aug-10/bing: coded
!
      real, dimension(6) :: tmppoint
      integer :: iost,rn,iol
      integer, intent(in) :: level
      logical :: ex
      character(len=20) :: filename
!
      write (filename,'("driver/pts_",I1.1,".dat")') level
!
      inquire(file=filename,exist=ex)
!
      if (ex) then
        inquire(IOLENGTH=iol) dy
        print*,'Reading granule level ',level,':',filename
        open(10,file=filename,status="unknown",access="direct",recl=6*iol)
        iost=0
!
        rn=1
        read (10, iostat=iost, rec=rn) tmppoint
        do while (iost == 0)
          call add_point
          current%pos(:) =tmppoint(1:2)
          current%data(:)=tmppoint(3:6)
          call drawupdate
          rn=rn+1
          read (10, iostat=iost, rec=rn) tmppoint
        enddo
!
        close(10)
        print*,'read ',rn-1,' points'
!
        if (level==nglevel) then
          write (filename,'("driver/seed_",I1.1,".dat")') level
          inquire(file=filename,exist=ex)
          if (ex) then
            open(10,file=filename,status="unknown",access="direct",recl=mseed*iol)
            read(10,rec=1) points_rstate
            close(10)
          else
            call fatal_error('rdpoints','cant find seed list for granules')
          endif
        endif
!
      endif
!
    endsubroutine rdpoints
!***********************************************************************
    subroutine wrpoints(level,issnap)
!
! Writes positions and amplitudes of granules to files. Also
! seed list are stored to be able to reproduce the run.
!
! 12-aug-10/bing: coded
!
      integer :: rn,iol
      integer, optional, intent(in) :: issnap
      integer, intent(in) :: level
      real, dimension(6) :: posdata
      character(len=22) :: filename
!
      inquire(IOLENGTH=iol) dy
!
      if (present(issnap)) then
        write (filename,'("driver/pts_",I1.1,"_",I3.3,".dat")') level,issnap
      else
        write (filename,'("driver/pts_",I1.1,".dat")') level
      endif
!
      open(10,file=filename,status="replace",access="direct",recl=6*iol)
!
      rn=1
      current => first
      do while (associated (current))
        posdata(1:2)=current%pos
        posdata(3:6)=current%data
!
        write(10,rec=rn) posdata
        rn=rn+1
        current => current%next
      enddo
!
      close(10)
!
! Save seed list for each level. Is needed if levels are spread over 3 procs.
!
      if (present(issnap)) then
        write (filename,'("driver/seed_",I1.1,"_",I3.3,".dat")') level,issnap
      else
        write (filename,'("driver/seed_",I1.1,".dat")') level
      endif
      !
      open(10,file=filename,status="replace",access="direct",recl=mseed*iol)
      write(10,rec=1) points_rstate
      close(10)
!
    endsubroutine wrpoints
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
    subroutine driveinit
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
      do while (minval (avoidarr) == 0)
!
        call add_point
        call make_newpoint
!
! Set randomly some points t0 to the past so they already decay
!
        call random_number_wrapper(rand)
        current%data(3)=t+(rand*2-1)*current%data(4)* &
            (-alog(thresh*ampl/current%data(2)))**(1./pow)
!
        current%data(1)=current%data(2)*exp(-((t-current%data(3))/current%data(4))**pow)
!
! Update arrays with new data
!
        call drawupdate
!
      enddo
!
    endsubroutine driveinit
!***********************************************************************
    subroutine helmholtz(frx_r,fry_r)
!
! Extracts the rotational part of a 2d vector field
! to increase vorticity of the velocity field.
!
! 12-aug-10/bing: coded
!
      use Fourier, only: fourier_transform_other
!
      real, dimension(nxgrid,nygrid) :: kx,ky,k2,filter
      real, dimension(nxgrid,nygrid) :: fvx_r,fvy_r,fvx_i,fvy_i
      real, dimension(nxgrid,nygrid) :: frx_r,fry_r,frx_i,fry_i
      real, dimension(nxgrid,nygrid) :: fdx_r,fdy_r,fdx_i,fdy_i
      real :: k20
!
      fvx_r=vx
      fvx_i=0.
      call fourier_transform_other(fvx_r,fvx_i)
!
      fvy_r=vy
      fvy_i=0.
      call fourier_transform_other(fvy_r,fvy_i)
!
! Reference frequency is half the Nyquist frequency.
      k20 = (kx_ny/2.)**2.
!
      kx =spread(kx_fft,2,nygrid)
      ky =spread(ky_fft,1,nxgrid)
!
      k2 =kx**2 + ky**2 + tini
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
      vx=fdx_r
      call fourier_transform_other(fdy_r,fdy_i,linv=.true.)
      vy=fdy_r
!
      call fourier_transform_other(frx_r,frx_i,linv=.true.)
      call fourier_transform_other(fry_r,fry_i,linv=.true.)
!
    endsubroutine helmholtz
!***********************************************************************
    subroutine drawupdate
!
! Using a point from the list to update the velocity field.
!
! 12-aug-10/bing: coded
!
      real :: xdist,ydist,dist2,dist,wtmp,vv
      integer :: i,ii,j,jj
      real :: dist0,tmp
!
! Update weight and velocity for new granule
!
      do jj=int(current%pos(2))-yrange,int(current%pos(2))+yrange
        j = 1+mod(jj-1+nygrid,nygrid)
        do ii=int(current%pos(1))-xrange,int(current%pos(1))+xrange
          i = 1+mod(ii-1+nxgrid,nxgrid)
          xdist=dx*(ii-current%pos(1))
          ydist=dy*(jj-current%pos(2))
          dist2=max(xdist**2+ydist**2,dxdy2)
          dist=sqrt(dist2)
!
          if (dist.lt.avoid*granr.and.t.lt.current%data(3)) avoidarr(i,j)=1
!
          wtmp=current%data(1)/dist
!
          dist0 = 0.53*granr
          tmp = (dist/dist0)**2
!
          vv=exp(1.)*current%data(1)*tmp*exp(-tmp)
!
          if (wtmp.gt.w(i,j)*(1-ig)) then
            if (wtmp.gt.w(i,j)*(1+ig)) then
              ! granular area
              vx(i,j)=vv*xdist/dist
              vy(i,j)=vv*ydist/dist
              w(i,j) =wtmp
            else
              ! intergranular area
              vx(i,j)=vx(i,j)+vv*xdist/dist
              vy(i,j)=vy(i,j)+vv*ydist/dist
              w(i,j) =max(w(i,j),wtmp)
            end if
          endif
          if (w(i,j) .gt. ampl/(granr*(1+ig))) avoidarr(i,j)=1
        enddo
      enddo
!
    endsubroutine drawupdate
!***********************************************************************
    subroutine make_newpoint
!
! Find the position of a new point.
!
! 12-aug-10/bing: coded
!
      use General, only: random_number_wrapper
!
      integer :: kfind,count,ipos,jpos,i,j
      integer,dimension(nxgrid,nygrid) :: k
      real :: rand
!
      k(:,:)=0; ipos=0; jpos=0
!
      where (avoidarr.eq.0) k=1
!
! Choose and find location of one of them
!
      call random_number_wrapper(rand)
      kfind=int(rand*sum(k))+1
      count=0
      do i=1,nxgrid
        do j=1,nygrid
          if (k(i,j).eq.1) then
            count=count+1
            if (count.eq.kfind) then
              ipos=i
              jpos=j
            endif
          endif
        enddo
      enddo
!
! Create new data for new point
!
      current%pos(1)=ipos
      current%pos(2)=jpos
!
      call random_number_wrapper(rand)
      current%data(2)=ampl*(1+(2*rand-1)*pd)
!
      call random_number_wrapper(rand)
      current%data(4)=life_t*(1+(2*rand-1)/10.)
!
      current%data(3)=t+current%data(4)* &
          (-alog(thresh*ampl/current%data(2)))**(1./pow)
!
      current%data(1)=current%data(2)* &
          exp(-((t-current%data(3))/current%data(4))**pow)
!
    endsubroutine make_newpoint
!***********************************************************************
    subroutine updatepoints
!
! Update the amplitude/weight of a point.
!
! 12-aug-10/bing: coded
!
      use Sub, only: notanumber
!
      current => first
      do while (associated (current))
        if (notanumber(current%data)) &
            call fatal_error('update points','NaN found',.true.)
!
! update amplitude
        current%data(1)=current%data(2)* &
            exp(-((t-current%data(3))/current%data(4))**pow)
!
! remove point if amplitude is less than threshold
        if (current%data(1)/ampl.lt.thresh) then
          call del_point
        else
          current => current%next
        endif
!
      end do
!
    endsubroutine updatepoints
!***********************************************************************
    subroutine set_B2(f,BB2_local)
!
      use Mpicomm, only: mpisend_real, mpirecv_real, stop_it_if_any
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(nx,ny) :: bbx,bby,bbz
      real, dimension(nx,ny) :: fac,BB2_local,tmp
      integer :: px,py,partner,alloc_err
      integer, dimension(2) :: dims
      real :: temp
!
      intent(in) :: f
      intent(out) :: BB2_local
!
! compute B = curl(A) for irefz layer
!
      if (nygrid/=1) then
        fac=(1./60)*spread(dy_1(m1:m2),1,nx)
        bbx= fac*(+ 45.0*(f(l1:l2,m1+1:m2+1,irefz,iaz)-f(l1:l2,m1-1:m2-1,irefz,iaz)) &
            -  9.0*(f(l1:l2,m1+2:m2+2,irefz,iaz)-f(l1:l2,m1-2:m2-2,irefz,iaz)) &
            +      (f(l1:l2,m1+3:m2+3,irefz,iaz)-f(l1:l2,m1-3:m2-3,irefz,iaz)))
      else
        if (ip<=5) print*, 'set_B2: Degenerate case in y-direction'
      endif
      if (nzgrid/=1) then
        fac=(1./60)*spread(spread(dz_1(irefz),1,nx),2,ny)
        bbx= bbx -fac*(+ 45.0*(f(l1:l2,m1:m2,irefz+1,iay)-f(l1:l2,m1:m2,irefz-1,iay)) &
            -  9.0*(f(l1:l2,m1:m2,irefz+2,iay)-f(l1:l2,m1:m2,irefz-2,iay)) &
            +      (f(l1:l2,m1:m2,irefz+3,iay)-f(l1:l2,m1:m2,irefz-2,iay)))
      else
        if (ip<=5) print*, 'set_B2: Degenerate case in z-direction'
      endif
!
      if (nzgrid/=1) then
        fac=(1./60)*spread(spread(dz_1(irefz),1,nx),2,ny)
        bby= fac*(+ 45.0*(f(l1:l2,m1:m2,irefz+1,iax)-f(l1:l2,m1:m2,irefz-1,iax)) &
            -  9.0*(f(l1:l2,m1:m2,irefz+2,iax)-f(l1:l2,m1:m2,irefz-2,iax)) &
            +      (f(l1:l2,m1:m2,irefz+3,iax)-f(l1:l2,m1:m2,irefz-3,iax)))
      else
        if (ip<=5) print*, 'set_B2: Degenerate case in z-direction'
      endif
      if (nxgrid/=1) then
        fac=(1./60)*spread(dx_1(l1:l2),2,ny)
        bby=bby-fac*(+45.0*(f(l1+1:l2+1,m1:m2,irefz,iaz)-f(l1-1:l2-1,m1:m2,irefz,iaz)) &
            -  9.0*(f(l1+2:l2+2,m1:m2,irefz,iaz)-f(l1-2:l2-2,m1:m2,irefz,iaz)) &
            +      (f(l1+3:l2+3,m1:m2,irefz,iaz)-f(l1-3:l2-3,m1:m2,irefz,iaz)))
      else
        if (ip<=5) print*, 'set_B2: Degenerate case in x-direction'
      endif
      if (nxgrid/=1) then
        fac=(1./60)*spread(dx_1(l1:l2),2,ny)
        bbz= fac*(+ 45.0*(f(l1+1:l2+1,m1:m2,irefz,iay)-f(l1-1:l2-1,m1:m2,irefz,iay)) &
            -  9.0*(f(l1+2:l2+2,m1:m2,irefz,iay)-f(l1-2:l2-2,m1:m2,irefz,iay)) &
            +      (f(l1+3:l2+3,m1:m2,irefz,iay)-f(l1-3:l2-3,m1:m2,irefz,iay)))
      else
        if (ip<=5) print*, 'set_B2: Degenerate case in x-direction'
      endif
      if (nygrid/=1) then
        fac=(1./60)*spread(dy_1(m1:m2),1,nx)
        bbz=bbz-fac*(+45.0*(f(l1:l2,m1+1:m2+1,irefz,iax)-f(l1:l2,m1-1:m2-1,irefz,iax)) &
            -  9.0*(f(l1:l2,m1+2:m2+2,irefz,iax)-f(l1:l2,m1-2:m2-2,irefz,iax)) &
            +      (f(l1:l2,m1+3:m2+3,irefz,iax)-f(l1:l2,m1-3:m2-3,irefz,iax)))
      else
        if (ip<=5) print*, 'set_B2: Degenerate case in y-direction'
      endif
!
      BB2_local = bbx*bbx + bby*bby + bbz*bbz
      Bzflux = sum(abs(bbz))
!
! communicate to root processor
!
      if (lroot) then
        dims=(/nx,ny/)
        if (.not. allocated (BB2)) then
          allocate (BB2(nxgrid,nygrid), stat=alloc_err)
          if (alloc_err > 0) call stop_it_if_any (.true., 'set_B2: Could not allocate memory')
          BB2 = 0.0
        endif
        BB2(1:nx,1:ny) = BB2_local
        do px=0,nprocx-1
          do py=0,nprocy-1
            partner = px + py*nprocx
            if (partner == 0) cycle
            call mpirecv_real(tmp,dims,partner,555+partner)
            BB2(px*nx+1:px*nx+nx,py*ny+1:py*ny+ny) = tmp
            call mpirecv_real(temp,1,partner,556+partner)
            Bzflux = Bzflux+temp
          enddo
        enddo
      elseif (lfirst_proc_z) then
        call mpisend_real(BB2_local,dims,0,555+iproc)
        call mpisend_real(Bzflux,1,0,556+iproc)
      endif
!
      if (lgran_parallel) then
        dims=(/nxgrid,nygrid/)
        if (lroot) then
          call mpisend_real(BB2,dims,nprocxy,nprocxy)
          call mpisend_real(BB2,dims,nprocxy+1,nprocxy+1)
          call mpisend_real(BB2,dims,nprocxy+2,nprocxy+2)
        elseif (lgran_proc) then
          if (.not. allocated (BB2)) then
            allocate (BB2(nxgrid,nygrid), stat=alloc_err)
            if (alloc_err>0) call stop_it_if_any(.true., &
                'set_B2: Could not allocate memory')
          endif
          call mpirecv_real(BB2,dims,0,iproc)
        endif
      endif
!
    endsubroutine set_B2
!***********************************************************************
    subroutine fill_B_avoidarr
!
      integer :: i,j,itmp,jtmp
      integer :: il,ir,jl,jr
      integer :: ii,jj
!
      if (nxgrid==1) then
        itmp = 0
      else
        itmp = nint(granr*(1-ig)/dx)
      endif
      if (nygrid==1) then
        jtmp = 0
      else
        jtmp = nint(granr*(1-ig)/dy)
      endif
!
      do i=1,nxgrid
        do j=1,nygrid
          if (BB2(i,j) > (Bavoid/unit_magnetic)**2) then
            il=max(1,i-itmp); ir=min(nxgrid,i+itmp)
            jl=max(1,j-jtmp); jr=min(nygrid,j+jtmp)
!
            do ii=il,ir
              do jj=jl,jr
                if ((ii-i)**2+(jj-j)**2 < itmp**2+jtmp**2) avoidarr(ii,jj)=1
              enddo
            enddo
          endif
        enddo
      enddo
!
    endsubroutine fill_B_avoidarr
!***********************************************************************
    subroutine read_ext_vel_field()
!
      use Mpicomm, only: mpisend_real, mpirecv_real, stop_it_if_any
!
      real, dimension (:,:), save, allocatable :: uxl,uxr,uyl,uyr
      real, dimension (:,:), allocatable :: tmpl,tmpr
      integer, parameter :: tag_x=321,tag_y=322
      integer, parameter :: tag_tl=345,tag_tr=346,tag_dt=347
      integer :: lend=0,ierr,i,alloc_err,px,py
      real, save :: tl=0.,tr=0.,delta_t=0.
!
      character (len=*), parameter :: vel_times_dat = 'driver/vel_times.dat'
      character (len=*), parameter :: vel_field_dat = 'driver/vel_field.dat'
      integer :: unit=1
!
      if (.not. allocated (uxl)) then
        allocate (uxl(nx,ny), uxr(nx,ny), uyl(nx,ny), uyr(nx,ny), stat=alloc_err)
        if (alloc_err > 0) call stop_it_if_any (.true., 'read_ext_vel_field: '// &
            'Could not allocate memory for velocity field variables')
      endif
!
      allocate (tmpl(nxgrid,nygrid), tmpr(nxgrid,nygrid), stat=alloc_err)
      if (alloc_err > 0) call stop_it_if_any (.true., 'read_ext_vel_field: '// &
          'Could not allocate memory for tmp variables, please check')
!
!  Read the time table
!
      if ((t*unit_time<tl+delta_t) .or. (t*unit_time>=tr+delta_t)) then
!
        if (lroot) then
!
          if (.not. allocated (Ux_ext)) then
            allocate (Ux_ext(nxgrid,nygrid), Uy_ext(nxgrid,nygrid), stat=alloc_err)
            if (alloc_err > 0) call stop_it_if_any (.true., 'read_ext_vel_field: Could not allocate memory')
            Ux_ext = 0.0
            Uy_ext = 0.0
          endif
!
          inquire(IOLENGTH=lend) tl
          open (unit,file=vel_times_dat,form='unformatted',status='unknown',recl=lend,access='direct')
!
          ierr = 0
          i=0
          do while (ierr == 0)
            i=i+1
            read (unit,rec=i,iostat=ierr) tl
            read (unit,rec=i+1,iostat=ierr) tr
            if (ierr /= 0) then
              ! EOF is reached => read again
              i=1
              delta_t = t*unit_time
              read (unit,rec=i,iostat=ierr) tl
              read (unit,rec=i+1,iostat=ierr) tr
              ierr=-1
            else
              ! test if correct time step is reached
              if (t*unit_time>=tl+delta_t.and.t*unit_time<tr+delta_t) ierr=-1
            endif
          enddo
          close (unit)
!
          do px=0, nprocx-1
            do py=0, nprocy-1
              if ((px == 0) .and. (py == 0)) cycle
              call mpisend_real (tl, 1, px+py*nprocx, tag_tl)
              call mpisend_real (tr, 1, px+py*nprocx, tag_tr)
              call mpisend_real (delta_t, 1, px+py*nprocx, tag_dt)
            enddo
          enddo
!
! Read velocity field
!
          open (unit,file=vel_field_dat,form='unformatted',status='unknown',recl=lend*nxgrid*nygrid,access='direct')
!
          read (unit,rec=2*i-1) tmpl
          read (unit,rec=2*i+1) tmpr
          if (tr /= tl) then
            Ux_ext = (t*unit_time - (tl+delta_t)) * (tmpr - tmpl) / (tr - tl) + tmpl
          else
            Ux_ext = tmpr
          endif
          Ux_ext = Ux_ext/unit_velocity
!
          read (unit,rec=2*i)   tmpl
          read (unit,rec=2*i+2) tmpr
          if (tr /= tl) then
            Uy_ext = (t*unit_time - (tl+delta_t)) * (tmpr - tmpl) / (tr - tl) + tmpl
          else
            Uy_ext = tmpr
          endif
          Uy_ext = Uy_ext/unit_velocity
!
          do px=0, nprocx-1
            do py=0, nprocy-1
              if ((px == 0) .and. (py == 0)) cycle
              Ux_ext_local = Ux_ext(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
              Uy_ext_local = Uy_ext(px*nx+1:(px+1)*nx,py*ny+1:(py+1)*ny)
              call mpisend_real (Ux_ext_local, (/ nx, ny /), px+py*nprocx, tag_x)
              call mpisend_real (Uy_ext_local, (/ nx, ny /), px+py*nprocx, tag_y)
            enddo
          enddo
!
          Ux_ext_local = Ux_ext(1:nx,1:ny)
          Uy_ext_local = Uy_ext(1:nx,1:ny)
!
          close (unit)
        else
          if (lfirst_proc_z) then
            call mpirecv_real (tl, 1, 0, tag_tl)
            call mpirecv_real (tr, 1, 0, tag_tr)
            call mpirecv_real (delta_t, 1, 0, tag_dt)
            call mpirecv_real (Ux_ext_local, (/ nx, ny /), 0, tag_x)
            call mpirecv_real (Uy_ext_local, (/ nx, ny /), 0, tag_y)
          endif
        endif
!
      endif
!
      deallocate (tmpl, tmpr)
!
    endsubroutine read_ext_vel_field
!***********************************************************************
    subroutine force_solar_wind(df,p)
!
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      if (n==n2.and.llast_proc_z) &
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
      use Mpicomm, only: mpisend_real,mpirecv_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: i,j,ipt
      real :: local_flux,local_mass
      real :: total_flux,total_mass
      real :: get_lf,get_lm
!
      local_flux=sum(exp(f(l1:l2,m1:m2,n2,ilnrho))*f(l1:l2,m1:m2,n2,iuz))
      local_mass=sum(exp(f(l1:l2,m1:m2,n2,ilnrho)))
!
!  One  processor has to collect the data
!
      if (lfirst_proc_xy) then
        total_flux=local_flux
        total_mass=local_mass
        do i=0,nprocx-1
          do j=0,nprocy-1
            if ((i==0).and.(j==0)) cycle
            ipt = i+nprocx*j+ipz*nprocx*nprocy
            call mpirecv_real(get_lf,1,ipt,111+ipt)
            call mpirecv_real(get_lm,1,ipt,211+ipt)
            total_flux=total_flux+get_lf
            total_mass=total_mass+get_lm
          enddo
        enddo
!
!  Get u0 addition rho*(u+u0) = wind
!  rho*u + u0 *rho =wind
!  u0 = (wind-rho*u)/rho
!
        u_add = (massflux-total_flux) / total_mass
      else
        ! send to first processor at given height
        !
        call mpisend_real(local_flux,1,ipz*nprocx*nprocy,111+iproc)
        call mpisend_real(local_mass,1,ipz*nprocx*nprocy,211+iproc)
      endif
!
!  now distribute u_add
!
      if (lfirst_proc_xy) then
        do i=0,nprocx-1
          do j=0,nprocy-1
            if ((i==0).and.(j==0)) cycle
            ipt = i+nprocx*j+ipz*nprocx*nprocy
            call mpisend_real(u_add,1,ipt,311+ipt)
          enddo
        enddo
      else
        call mpirecv_real(u_add,1,ipz*nprocx*nprocy,311+iproc)
      endif
!
    endsubroutine get_wind_speed_offset
!***********************************************************************
    subroutine evolve_granules()
!
      integer :: xpos,ypos
!
      current => first
      do while (associated (current))
        xpos = int(current%pos(1))
        ypos = int(current%pos(2))
!
        current%pos(1) =  current%pos(1) + Ux_ext(xpos,ypos)*dt_gran
        current%pos(2) =  current%pos(2) + Uy_ext(xpos,ypos)*dt_gran
!
        if (current%pos(1)<0.5) current%pos(1) = current%pos(1) + nxgrid
        if (current%pos(2)<0.5) current%pos(2) = current%pos(2) + nygrid
!
        if (current%pos(1)>nxgrid+0.5) current%pos(1) = current%pos(1) - nxgrid
        if (current%pos(2)>nygrid+0.5) current%pos(2) = current%pos(2) - nygrid
!
        current => current%next
      enddo
!
    endsubroutine evolve_granules
!***********************************************************************
    subroutine enhance_vorticity()
!
      real,dimension(nxgrid,nygrid) :: wscr,wscr2
      real :: vrms,vtot
!
! Putting sum of velocities back into vx,vy
        vx=Ux
        vy=Uy
!
! Calculating and enhancing rotational part by factor 5
        if (lrotin) then
          call helmholtz(wscr,wscr2)
          !* war vorher 5 ; zum testen auf  50
          ! nvor is now keyword !!!
          vx=(vx+nvor*wscr )
          vy=(vy+nvor*wscr2)
        endif
!
! Normalize to given total rms-velocity
        vrms=sqrt(sum(vx**2+vy**2)/(nxgrid*nygrid))+tini
!
        if (unit_system.eq.'SI') then
          vtot=3.*1e3/unit_velocity
        elseif (unit_system.eq.'cgs') then
          vtot=3.*1e5/unit_velocity
        else
          vtot=0.
          call fatal_error('solar_corona','define a valid unit system')
        endif
!
! Reinserting rotationally enhanced velocity field
!
        Ux=vx*vtot/vrms
        Uy=vy*vtot/vrms
!
    endsubroutine enhance_vorticity
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
