! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
  character (len=labellen) :: stratitype='nothing'
  character (len=labellen) :: lnrho_init='nothing'
  character (len=labellen) :: lnTT_init='nothing'
  character (len=labellen) :: aa_init='nothing'
  logical :: set_lnTT_first =.true.
  real :: rho_init=0.,const_alfven=1.
  real :: T0=6000.,T1=1e6,z0_tanh=4e6,width_tanh=1e6
  character (len=labellen) :: direction='z'
  real, dimension(4) :: mpoly_special = (/1.3,1000.,-1.04,500./)
  real, dimension(3) :: zpoly = (/0.,3.,5./)
!
  namelist /initial_condition_pars/ &
      lnrho_init,lnTT_init,stratitype,rho_init,direction, &
      set_lnTT_first,T0,T1,z0_tanh,width_tanh,mpoly_special,zpoly, &
      const_alfven
!
  real :: gamma, gamma_m1, cp1

contains
!***********************************************************************
  subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  04-sep-10/bing: coded
!
    if (lroot) call svn_id( &
        "$Id$")
!
  endsubroutine register_initial_condition
!***********************************************************************
  subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  14-dec-10/bing: coded
!
    use EquationOfState, only: get_gamma_etc

    real, dimension (mx,my,mz,mfarray) :: f
    real :: cp
!
    if (iproc==0) then
      write(*,*) "-------------------------------------------------------------"
      write(*,*) "Parameters to be set in run.in:"
      write(*,*)
      write(*,'(A,ES10.2)') "Kpara=",2e-11 /unit_density/unit_velocity**3./ &
          unit_length*unit_temperature**3.5
      write(*,'(A,ES10.2)') "Kperp=",3.47e12/((unit_velocity**3*unit_magnetic**2 &
          *unit_length)/(unit_density * sqrt(unit_temperature)))
!
      write(*,'(A,ES10.2)') "eta <=",1e7/unit_velocity/unit_length
!
      write(*,'(A,ES10.2)') "hcond_grad=",1e9*dxmax**3*unit_temperature/unit_velocity**3
      write(*,'(A,ES10.2)') 'Assume Umax=150 km/s and R=1 => diff=',150e3/unit_velocity*dxmax
      write(*,'(A,ES10.2)') "nu_spitzer=",2.21e-17*unit_temperature**2.5/ &
          unit_density/unit_velocity/unit_length
      write(*,*) "-------------------------------------------------------------"
    endif
!
    call get_gamma_etc(gamma,cp)
    gamma_m1=gamma-1.
    cp1=1./cp

    call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
  subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density.
!
!  04-sep-10/bing: coded
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    if (lnrho_init=='hydrostatic'.and.lnTT_init=='hydrostatic') &
        call fatal_error('initial_condition_lnrho','doesnt work')
!
    if (stratitype=='hydrostatic'.and.direction=='x') &
        call hydrostatic_x(f)
!
    select case (lnTT_init)
    case ('nothing')
      ! do nothing
    case ('prof_lnTT')
      call setup_vert_profiles(f)
    case ('tanh')
      call setup_tanh(f)
    case ('piecewice_poly')
      call piecewice_poly(f)
    case default
      call fatal_error('initial_condition_lnrho', &
          'no such value for lnTT_init')
    endselect
!
    select case (lnrho_init)
    case ('nothing')
      ! do nothing
    case ('prof_lnrho')
      call setup_vert_profiles(f)
    case ('hydrostatic')
      call hydrostatic_lnTT(f)
    case ('alfven-const')
      ! do nothing, will be done after init_aa
    case default
      call fatal_error('initial_condition_lnrho', &
          'no such value for lnTT_init')
    endselect
!
    call write_stratification_dat(f)
!
  endsubroutine initial_condition_lnrho
!***********************************************************************
  subroutine initial_condition_aa(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      if (lnrho_init=='alfven-const') call const_alfven_speed(f)
!
  endsubroutine initial_condition_aa
!***********************************************************************
  subroutine setup_vert_profiles(f)
!
!  Read and set vertical profiles for initial temperature and density.
!  Initial temperature profile is given in ln(T) [K] over z [m]
!  Initial density profile is given in ln(rho) [kg/m^3] over z [m]
!
!  04-sep-10/bing: coded
!
      use EquationOfState, only: cs20
      use File_io, only: file_exists, file_size
      use Mpicomm, only: mpibcast_int, mpibcast_real, stop_it_if_any
      use Messages, only: warning
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: lend,lend_b8,ierr
      integer :: i,j
      integer, parameter :: unit=12
!
! file location settings
      character (len=*), parameter :: lnrho_dat = 'prof_lnrho.dat'
      character (len=*), parameter :: lnT_dat = 'prof_lnT.dat'
!
      integer :: prof_nz
      real, dimension (:), allocatable :: prof_z, prof_lnrho, prof_lnTT
      real, dimension (mz) :: profile_z
      real, dimension (mx) :: profile_x
      logical :: lread_lnrho=.false., lread_lnTT=.false.
!
      inquire(IOLENGTH=lend) 1.0
      inquire(IOLENGTH=lend_b8) 1.0d0
!
      lread_lnTT=(lnTT_init=='prof_lnTT')
      lread_lnrho=(lnrho_init=='prof_lnrho')
!
! read density profile for interpolation
      if (lread_lnrho) then
!
! file access is only done on the MPI root rank
        if (lroot) then
          if (.not. file_exists (lnrho_dat)) call stop_it_if_any ( &
              .true., 'setup_special: file not found: '//trim(lnrho_dat))
! find out, how many data points our profile file has
          prof_nz = (file_size (lnrho_dat) - 2*2*4) / (lend*8/lend_b8 * 2)
        endif
        call stop_it_if_any(.false.,'')
        call mpibcast_int (prof_nz)
!
        allocate (prof_z(prof_nz), prof_lnrho(prof_nz), stat=ierr)
!
        if (lroot) then
          open (unit,file=lnrho_dat,form='unformatted',status='unknown', &
              recl=lend*prof_nz)
          read (unit,iostat=ierr) prof_lnrho
          read (unit,iostat=ierr) prof_z
          if (ierr /= 0) call stop_it_if_any(.true.,'setup_special: '// &
              'Error reading stratification file: "'//trim(lnrho_dat)//'"')
          close (unit)
        endif
        call stop_it_if_any(.false.,'')
!
        call mpibcast_real (prof_lnrho,prof_nz)
        call mpibcast_real (prof_z,prof_nz)
!
! convert from logarithmic SI to Pencil units
        prof_lnrho = prof_lnrho - log(real(unit_density))
!
! convert z coordinates from [m] to Pencil units
        prof_z = prof_z / unit_length
!
! interpolate density profile to Pencil grid
        if (direction=='z') then
          do j = n1-nghost, n2+nghost
            if (z(j) < prof_z(1) ) then
              profile_z(j) = prof_lnrho(1)
            elseif (z(j) >= prof_z(prof_nz)) then
              profile_z(j) = prof_lnrho(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((z(j) >= prof_z(i)) .and. (z(j) < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  profile_z(j) = (prof_lnrho(i+1)-prof_lnrho(i)) / &
                      (prof_z(i+1)-prof_z(i)) * (z(j)-prof_z(i)) + prof_lnrho(i)
                  exit
                endif
              enddo
            endif
          enddo
        elseif (direction=='x') then
          do j = l1-nghost, l2+nghost
            if (x(j) < prof_z(1) ) then
              profile_x(j) = prof_lnrho(1)
            elseif (x(j) >= prof_z(prof_nz)) then
              profile_x(j) = prof_lnrho(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((x(j) >= prof_z(i)) .and. (x(j) < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  profile_x(j) = (prof_lnrho(i+1)-prof_lnrho(i)) / &
                      (prof_z(i+1)-prof_z(i)) * (x(j)-prof_z(i)) + prof_lnrho(i)
                  exit
                endif
              enddo
            endif
          enddo
        endif
!
        if (allocated (prof_z)) deallocate (prof_z)
        if (allocated (prof_lnrho)) deallocate (prof_lnrho)
!
! REMARK: f-array is filled with log value even for ldensity_nolog=T
!         because at the end of init_lnrho all values in the f-array are
!         converted by  f = exp(f)
!
        if (ldensity_nolog) then
          if (direction=='z') f(:,:,:,irho)=spread(spread(profile_z,1,mx),2,my)
          if (direction=='x') f(:,:,:,irho)=spread(spread(profile_x,2,my),3,mz)
        else
          if (direction=='z') f(:,:,:,ilnrho)=spread(spread(profile_z,1,mx),2,my)
          if (direction=='x') f(:,:,:,ilnrho)=spread(spread(profile_x,2,my),3,mz)
        endif
      endif
!
! read temperature profile for interpolation
      if (lread_lnTT) then
!
! file access is only done on the MPI root rank
        if (lroot) then
          if (.not. file_exists (lnT_dat)) call stop_it_if_any ( &
              .true., 'setup_special: file not found: '//trim(lnT_dat))
! find out, how many data points our profile file has
          prof_nz = (file_size (lnT_dat) - 2*2*4) / (lend*8/lend_b8 * 2)
        endif
        call stop_it_if_any(.false.,'')
        call mpibcast_int(prof_nz)
!
        allocate (prof_z(prof_nz), prof_lnTT(prof_nz), stat=ierr)
!
        if (lroot) then
          open (unit,file=lnT_dat,form='unformatted',status='unknown', &
              recl=lend*prof_nz)
          read (unit,iostat=ierr) prof_lnTT
          read (unit,iostat=ierr) prof_z
          if (ierr /= 0) call stop_it_if_any(.true.,'setup_special: '// &
              'Error reading stratification file: "'//trim(lnT_dat)//'"')
          close (unit)
        endif
        call stop_it_if_any(.false.,'')
!
        call mpibcast_real(prof_lnTT,prof_nz)
        call mpibcast_real(prof_z,prof_nz)
!
! convert from logarithmic SI to Pencil units
        prof_lnTT = prof_lnTT - log(real(unit_temperature))
!
! convert z coordinates from [m] to Pencil units
        prof_z = prof_z / unit_length
!
! interpolate temperature profile to Pencil grid
!
        if (direction=='z') then
          do j = n1-nghost, n2+nghost
            if (z(j) < prof_z(1) ) then
              profile_z(j) = prof_lnTT(1)
            elseif (z(j) >= prof_z(prof_nz)) then
              profile_z(j) = prof_lnTT(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((z(j) >= prof_z(i)) .and. (z(j) < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  profile_z(j) = (prof_lnTT(i+1)-prof_lnTT(i)) / &
                      (prof_z(i+1)-prof_z(i)) * (z(j)-prof_z(i)) + prof_lnTT(i)
                  exit
                endif
              enddo
            endif
          enddo
        else if (direction=='x') then
          do j = l1-nghost, l2+nghost
            if (x(j) < prof_z(1) ) then
              profile_x(j) = prof_lnTT(1)
            elseif (x(j) >= prof_z(prof_nz)) then
              profile_x(j) = prof_lnTT(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((x(j) >= prof_z(i)) .and. (x(j) < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  profile_x(j) = (prof_lnTT(i+1)-prof_lnTT(i)) / &
                      (prof_z(i+1)-prof_z(i)) * (x(j)-prof_z(i)) + prof_lnTT(i)
                  exit
                endif
              enddo
            endif
          enddo
        endif
!
        if (allocated (prof_z)) deallocate (prof_z)
        if (allocated (prof_lnTT)) deallocate (prof_lnTT)
!
        if (ltemperature) then
          if (ltemperature_nolog) then
            if (direction=='z') f(:,:,:,iTT)=spread(spread(exp(profile_z),1,mx),2,my)
            if (direction=='x') f(:,:,:,iTT)=spread(spread(exp(profile_x),2,my),3,mz)
          else
            if (direction=='z') f(:,:,:,ilnTT)=spread(spread(profile_z,1,mx),2,my)
            if (direction=='x') f(:,:,:,ilnTT)=spread(spread(profile_x,2,my),3,mz)
          endif
        else if (lthermal_energy) then
          if (direction=='z') f(:,:,:,ieth)=spread(spread(exp(profile_z),1,mx),2,my)
          if (direction=='x') f(:,:,:,ieth)=spread(spread(exp(profile_x),2,my),3,mz)
!
! No difference for ldensity_nolog because f-array values are still log-values
!
          if (ldensity_nolog) then
            f(:,:,:,ieth)=f(:,:,:,ieth)*exp(f(:,:,:,irho))  /(gamma*cp1)
          else
            f(:,:,:,ieth)=f(:,:,:,ieth)*exp(f(:,:,:,ilnrho))/(gamma*cp1)
          endif
        else if (lentropy) then
          if (direction=='z') then
            f(:,:,:,iss) = (log(gamma_m1/cs20/cp1)+spread(spread(profile_z,1,mx),2,my)- &
                gamma_m1*(f(:,:,:,ilnrho)-log(rho_init))) / cp1 /gamma
          endif
          if (direction=='x') then
            f(:,:,:,iss) = (log(gamma_m1/cs20/cp1)+spread(spread(profile_x,2,my),3,mz)- &
                gamma_m1*(f(:,:,:,ilnrho)-log(rho_init))) / cp1 /gamma
          endif
        else
          call fatal_error('setup_vert_profiles', &
              'Not implemented for current set of thermodynamic variables.')
        endif
      endif
!
    endsubroutine setup_vert_profiles
!***********************************************************************
    subroutine hydrostatic_z(f)
!
!  Intialize the density for given temperprofile in vertical
!  z direction by solving hydrostatic equilibrium.
!  dlnrho = - dlnTT + (cp-cv)/T g dz
!
!  The initial densitiy lnrho must be given in SI units.
!  Temperature given as function lnT(z) in SI units
!  [T] = K   &   [z] = Mm   & [rho] = kg/m^3
!
      use EquationOfState, only: cs2top,cs2bot
      use File_io, only: file_exists, file_size
      use Gravity, only: gravz
      use Mpicomm, only: mpibcast_real,mpibcast_int,stop_it_if_any
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: ztop,zbot
      integer :: prof_nz
      real, dimension(:), allocatable :: prof_lnTT,prof_z
      real :: tmp_lnrho,tmp_lnT,tmpdT,tmp_z,dz_step,lnrho_0
      integer :: i,lend,lend_b8,j,ierr,unit=1
!
! file location settings
      character (len=*), parameter :: lnT_dat = 'prof_lnT.dat'
!
      inquire(IOLENGTH=lend) 1.0
      inquire(IOLENGTH=lend_b8) 1.0d0
!
      if (lentropy.or.ltemperature_nolog.or.lthermal_energy.or.ldensity_nolog) &
          call fatal_error('hydrostatic','only implemented for ltemperature')
!
      lnrho_0 = log(rho_init)
!
! read temperature profile for interpolation
!
! file access is only done on the MPI root rank
      if (lroot) then
        if (.not. file_exists (lnT_dat)) call stop_it_if_any ( &
            .true., 'setup_special: file not found: '//trim(lnT_dat))
! find out, how many data points our profile file has
        prof_nz = (file_size (lnT_dat) - 2*2*4) / (lend*8/lend_b8 * 2)
      endif
!
      call stop_it_if_any(.false.,'')
      call mpibcast_int (prof_nz)
!
      allocate (prof_z(prof_nz), prof_lnTT(prof_nz), stat=ierr)
      if (ierr > 0) call stop_it_if_any (.true.,'setup_special: '// &
          'Could not allocate memory for z coordinate or lnTT profile')
!
      if (lroot) then
        open (unit,file=lnT_dat,form='unformatted',status='unknown', &
            recl=lend*prof_nz)
        read (unit,iostat=ierr) prof_lnTT
        read (unit,iostat=ierr) prof_z
        if (ierr /= 0) call stop_it_if_any(.true.,'setup_special: '// &
            'Error reading stratification file: "'//trim(lnT_dat)//'"')
        close (unit)
      endif
      call stop_it_if_any(.false.,'')
!
      call mpibcast_real (prof_lnTT,prof_nz)
      call mpibcast_real (prof_z,prof_nz)
!
      prof_z = prof_z*1.e6/unit_length
      prof_lnTT = prof_lnTT - log(real(unit_temperature))
!
! get step width
! should be smaler than grid width and
! data width
!
      dz_step = min((prof_z(2)-prof_z(1)),minval(1./dz_1))
      dz_step = dz_step/10.
!
      do j=n1,n2
        tmp_lnrho = lnrho_0
        tmp_lnT = prof_lnTT(1)
        tmp_z = prof_z(1)
!
        ztop = xyz0(3)+Lxyz(3)
        zbot = xyz0(3)
!
        do while (tmp_z <= ztop)
!
!  Set sound speed at the boundaries.
          if (abs(tmp_z-zbot) < dz_step) cs2bot = (gamma-1.)*exp(tmp_lnT)
          if (abs(tmp_z-ztop) < dz_step) cs2top = (gamma-1.)*exp(tmp_lnT)
!
          if (abs(tmp_z-z(j)) <= dz_step) then
            f(:,:,j,ilnrho) = tmp_lnrho
            f(:,:,j,ilnTT)  = tmp_lnT
          endif
!  new z coord
          tmp_z = tmp_z + dz_step
!  get T at new z
          do i=1,prof_nz-1
            if (tmp_z >= prof_z(i)  .and. tmp_z < prof_z(i+1) ) then
              tmpdT = (prof_lnTT(i+1)-prof_lnTT(i))/(prof_z(i+1)-prof_z(i)) * &
                  (tmp_z-prof_z(i)) + prof_lnTT(i) -tmp_lnT
              tmp_lnT = tmp_lnT + tmpdT
            elseif (tmp_z >= prof_z(prof_nz)) then
              tmpdT = prof_lnTT(prof_nz) - tmp_lnT
              tmp_lnT = tmp_lnT + tmpdT
            endif
          enddo
          tmp_lnrho=tmp_lnrho-tmpdT+gamma/(gamma-1.)*gravz*exp(-tmp_lnT)*dz_step
        enddo
      enddo
!
    endsubroutine hydrostatic_z
!***********************************************************************
  subroutine hydrostatic_lnTT(f)
!
! Solves the hydrostatic equilibrium using the
! temperature from the f-array.
! Integration is done using the trapezoidal rule.
!
    use Gravity, only: gravz
    use Mpicomm, only: mpisend_real,mpirecv_real
!
    real, dimension (mx,my,mz,mfarray) :: f
    real :: konst,lnrho_0,int
    real, dimension(mz) :: TT,lnTT
    integer :: i,j,k,ipt,ii
!
    konst = gamma*cp1/gamma_m1
!
    lnrho_0=log(rho_init)
    f(:,:,n1,ilnrho)= lnrho_0
!
    if (ltemperature_nolog) then
      TT = f(l1,m1,:,iTT)
      lnTT = log(f(l1,m1,:,iTT))
    else
      TT = exp(f(l1,m1,:,ilnTT))
      lnTT = f(l1,m1,:,ilnTT)
    endif
!
    do ii=0,nprocz-1
      if (ipz==ii) then
        do i=n1+1,n2+nghost
          int = 0.5 * (z(i)-z(i-1)) * (gravz*(1/TT(i-1)+1/TT(i)))
          f(:,:,i,ilnrho)=f(:,:,i-1,ilnrho)-lnTT(i)+ &
              lnTT(i-1)+konst*int
        enddo
        if (lfirst_proc_xy.and.ii<nprocz-1) then
          do j=0,nprocx-1
            do k=0,nprocy-1
              ipt=j + nprocx*k+nprocxy*(ii+1)
              call mpisend_real(f(l1,m1,n2+1,ilnrho),ipt,ipt)
            enddo
          enddo
        endif
      elseif (ipz==ii+1.and.ii<nprocz-1) then
        call mpirecv_real(lnrho_0,nprocxy*(ipz-1),iproc)
        f(:,:,n1,ilnrho) = lnrho_0
      endif
    enddo
!
!  Fill the lower most ghost cells. Can be overriden by standard
!  boundary conditions.
!
    if (ipz==0) then
      do i=n1-1,1,-1
        int = 0.5 * (z(i)-z(i+1)) * (gravz*(1/TT(i+1)+1/TT(i)))
        f(:,:,i,ilnrho)=f(:,:,i+1,ilnrho)-lnTT(i)+ &
            lnTT(i+1)+konst*int
      enddo
    endif
!
  endsubroutine hydrostatic_lnTT
!***********************************************************************
  subroutine hydrostatic_x(f)
!
!  Intialize the density for given temperprofile in vertical
!  z direction by solving hydrostatic equilibrium.
!  dlnrho = - dlnTT + (cp-cv)/T g dz
!
!  The initial densitiy lnrho0 must be given in SI units.
!  Temperature given as function lnT(z) in SI units
!  [T] = K   &   [z] = Mm   & [rho] = kg/m^3
!
    use File_io, only: file_exists, file_size
    use Gravity, only: get_xgravity
    use Mpicomm, only: mpibcast_real,mpibcast_int,stop_it_if_any
!
    real, dimension (mx,my,mz,mfarray) :: f
    integer :: prof_nx
    real, dimension(:), allocatable :: prof_lnTT,prof_x
    real :: tmp_lnrho,lnrho_0,ztmp,integrand
    integer :: i,lend,lend_b8,j,ierr,unit=1
    real, dimension(mx) :: xgrav,lnTT_loop
!
! file location settings
    character (len=*), parameter :: lnT_dat = 'prof_lnT.dat'
!
    call get_xgravity(xgrav)
!
    inquire(IOLENGTH=lend) 1.0
    inquire(IOLENGTH=lend_b8) 1.0d0
!
    if (lentropy.or.ltemperature_nolog) &
        call fatal_error('hydrostatic','only implemented for ltemperature')
!
    lnrho_0 = log(rho_init)
!
! read temperature profile for interpolation
!
! file access is only done on the MPI root rank
    if (lroot) then
      if (.not. file_exists (lnT_dat)) call stop_it_if_any ( &
          .true., 'setup_special: file not found: '//trim(lnT_dat))
! find out, how many data points our profile file has
      prof_nx = (file_size (lnT_dat) - 2*2*4) / (lend*8/lend_b8 * 2)
    endif
!
    call stop_it_if_any(.false.,'')
    call mpibcast_int (prof_nx)
!
    allocate (prof_x(prof_nx), prof_lnTT(prof_nx), stat=ierr)
    if (ierr > 0) call stop_it_if_any (.true.,'setup_special: '// &
        'Could not allocate memory for x coordinate or lnTT profile')
!
    if (lroot) then
      open (unit,file=lnT_dat,form='unformatted',status='unknown', &
          recl=lend*prof_nx)
      read (unit,iostat=ierr) prof_lnTT
      read (unit,iostat=ierr) prof_x
      if (ierr /= 0) call stop_it_if_any(.true.,'setup_special: '// &
          'Error reading stratification file: "'//trim(lnT_dat)//'"')
      close (unit)
    endif
    call stop_it_if_any(.false.,'')
!
    call mpibcast_real (prof_lnTT,prof_nx)
    call mpibcast_real (prof_x,prof_nx)
    !
    prof_x = prof_x/unit_length
    prof_lnTT = prof_lnTT - log(real(unit_temperature))
!
!  project T profile onto the loop
!
    do i=1,mx
      ztmp = Lxyz(1) / pi *sin( (x(i)-x(l1))/Lxyz(1)*pi )
!
      do j=1,prof_nx-1
        if ( ztmp .lt. prof_x(1) )  then
          lnTT_loop(i) = prof_lnTT(1)
        elseif( ztmp .ge. prof_x(j) .and.  ztmp .lt. prof_x(j+1)) then
          lnTT_loop(i) = lin_inpol(prof_lnTT(j+1),prof_lnTT(j), &
                                   prof_x(j+1),prof_x(j),ztmp)
        elseif  ( ztmp .ge. prof_x(prof_nx) )  then
          lnTT_loop(i) = prof_lnTT(prof_nx)
        endif
      enddo
!
    enddo
!
    tmp_lnrho = lnrho_0
!
    f(1,:,:,ilnrho) = tmp_lnrho
    f(1,:,:,ilnTT) = lnTT_loop(1)
!
    do j=2,mx
!
      integrand = 0.5*(x(j)-x(j-1))*(xgrav(j)*exp(-lnTT_loop(j)) + &
          xgrav(j-1)*exp(-lnTT_loop(j-1)) )
!
      tmp_lnrho = f(j-1,m1,n1,ilnrho)+ lnTT_loop(j-1)-lnTT_loop(j)+ &
          gamma/(gamma-1.)*cp1*integrand
!
      f(j,:,:,ilnrho) = tmp_lnrho
      f(j,:,:,ilnTT) = lnTT_loop(j)
!
    enddo
!
  endsubroutine hydrostatic_x
!***********************************************************************
  subroutine write_stratification_dat(f)
!
!  Writes the initial density temperature stratification into each
!  proc subfolder.
!
    real, dimension (mx,my,mz,mfarray), intent(in) :: f
    integer :: unit=12,lend
    real, dimension (mx) :: xwrite_density=0.,xwrite_energy=0.
    real, dimension (mz) :: zwrite_density=0.,zwrite_energy=0.
    real :: dummy=1.
!
    character (len=*), parameter :: filename='/strat.dat'
!
    inquire(IOLENGTH=lend) dummy
!
! Z - Direction
!
    if (direction=='z') then
      if (ldensity) zwrite_density=f(l1,m1,:,ilnrho)
      if (lentropy) zwrite_energy=f(l1,m1,:,iss)
      if (ltemperature) then
        if (ltemperature_nolog) then
          zwrite_energy=f(l1,m1,:,iTT)
        else
          zwrite_energy=f(l1,m1,:,ilnTT)
        endif
      endif
!
      open(unit,file=trim(directory_snap)//filename, &
          form='unformatted',status='unknown',recl=lend*mz)
      write(unit) z
      write(unit) zwrite_density
      write(unit) zwrite_energy
      close(unit)
!
! X - Direction
!
    elseif (direction=='x') then
      if (ldensity) xwrite_density=f(:,m1,n1,ilnrho)
      if (lentropy) xwrite_energy=f(:,m1,n1,iss)
      if (ltemperature) then
        if (ltemperature_nolog) then
          xwrite_energy=f(:,m1,n1,iTT)
        else
          xwrite_energy=f(:,m1,n1,ilnTT)
        endif
      endif
!
      open(unit,file=trim(directory_snap)//filename, &
          form='unformatted',status='unknown',recl=lend*mx)
      write(unit) x
      write(unit) xwrite_density
      write(unit) xwrite_energy
      close(unit)
    endif
!
  endsubroutine write_stratification_dat
!***********************************************************************
  function lin_inpol(f2,f1,x2,x1,yn)
!
    real :: f2,f1,x2,x1,yn
    real :: lin_inpol
!
    lin_inpol = (f2-f1)/(x2-x1)*(yn-x1) + f1
!
  endfunction lin_inpol
!***********************************************************************
  subroutine setup_tanh(f)
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f

    real, dimension (mz) :: TT,z_SI,TT_var
    integer :: i,j
!
    z_SI = z*unit_length
!
    TT = (T1-T0)*(0.5*tanh((z_SI-z0_tanh)/width_tanh)+0.5)+T0
!
    if (ltemperature) then
      if (ltemperature_nolog) then
        TT_var = TT / unit_temperature
      else
        TT_var = log(TT / unit_temperature)
      endif
    else
      TT_var = impossible
      call fatal_error('setup_tanh','only works for ltemperature=T')
    endif
!
    do i=1,mx
      do j=1,my
        f(i,j,:,ilnTT) = TT_var
      enddo
    enddo
!
  endsubroutine setup_tanh
!***********************************************************************
  subroutine piecewice_poly(f)
!
    use Gravity, only: gravz
!
    real, dimension(mx,my,mz,mfarray), intent(inout) :: f
    real :: Ttop,T2, T1, T0, temp
    real :: lnrhotop, lnrho2, lnrho1, lnrho0, ztop
    real :: lnrhobot,zbot,Tbot
    real, dimension(4) :: beta
    integer :: i
!
!  Top boundary values.
!
      ztop=xyz0(3)+Lxyz(3)
      zbot=xyz0(3)
!
!  Temperature gradients.
!
      beta = cp1*gravz/(mpoly_special+1.)*gamma/gamma_m1
!
!
      T0 = 6000./unit_temperature
      lnrho0 = log(3e-4/real(unit_density))
!
      Tbot = T0 - beta(1)*(zpoly(1)-zbot)
      T1   = T0 + beta(2)*(zpoly(2)-zpoly(1))
      T2   = T1 + beta(3)*(zpoly(3)-zpoly(2))
      Ttop = T2 + beta(4)*(ztop-zpoly(3))
!
!
      lnrhobot =  lnrho0+mpoly_special(1)*log(Tbot/T0)
      lnrho1   =  lnrho0+mpoly_special(2)*log(T1/T0)
      lnrho2   =  lnrho1+mpoly_special(3)*log(T2/T1)
      lnrhotop =  lnrho2+mpoly_special(4)*log(Ttop/T2)
!
      if (iproc==0) then
      print*,'########################################'
      print*,'Beta',beta
      print*,"TTbot",Tbot*unit_temperature
      print*,"TT1",T0*unit_temperature
      print*,"TT2",T1*unit_temperature
      print*,"TT3",T2*unit_temperature
      print*,"TTtop",Ttop*unit_temperature
      print*,"rhobot",exp(lnrhobot)*unit_density
      print*,"rho0",exp(lnrho0)*unit_density
      print*,"rho1",exp(lnrho1)*unit_density
      print*,"rho2",exp(lnrho2)*unit_density
      print*,"rhotop",exp(lnrhotop)*unit_density
      print*,'########################################'
      endif
      do  i=1,mz
!
        if (z(i) < zpoly(1)) then
          Temp = Tbot + beta(1)*(z(i)-zbot)
          f(:,:,i,ilnTT)=log(Temp)
          f(:,:,i,ilnrho)=lnrhobot+mpoly_special(1)*log(Temp/Tbot)
!
        elseif (z(i) >=zpoly(1) .and. z(i) < zpoly(2)) then
          Temp = T0 + beta(2)*(z(i)-zpoly(1))
          f(:,:,i,ilnTT)=log(Temp)
          f(:,:,i,ilnrho)=lnrho0+mpoly_special(2)*log(Temp/T0)
!
        elseif (z(i) >= zpoly(2) .and. z(i) <zpoly(3)) then
          Temp = T1 + beta(3)*(z(i)-zpoly(2))
          f(:,:,i,ilnTT)=log(Temp)
          f(:,:,i,ilnrho)=lnrho1+mpoly_special(3)*log(Temp/T1)
!
        elseif (z(i) >= zpoly(3)) then
          Temp = T2 + beta(4)*(z(i)-zpoly(3))
          f(:,:,i,ilnTT)=log(Temp)
          f(:,:,i,ilnrho)=lnrho2+mpoly_special(4)*log(Temp/T2)
!
        endif
!
      enddo
!
  endsubroutine piecewice_poly
!***********************************************************************
  subroutine const_alfven_speed(f)
!
! Toutine to set the alfven velocity constant everywhere in the box
! The alfven velocity is  v = sqrt(b^2/(mu0 rho))
!
    use Boundcond, only: boundconds_x, boundconds_y, boundconds_z 
    use Mpicomm, only: initiate_isendrcv_bdry, finalize_isendrcv_bdry
    use Sub, only: gij,curl_mn,dot2_mn
!
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    real, dimension(nx) :: b2
    real, dimension(nx,3) :: aa,bb
    real, dimension(nx,3,3) :: aij
!
! Make sure the internal boundaries are set properly 
! before computing the magnetic field.
!
    call boundconds_x(f,iax,iaz)
    call initiate_isendrcv_bdry(f,iax,iaz)
    call finalize_isendrcv_bdry(f,iax,iaz)
    call boundconds_y(f,iax,iaz)
    call boundconds_z(f,iax,iaz)
!
! Compute
!
    do m=m1,m2
      do n=n1,n2
        aa=f(l1:l2,m,n,iax:iaz)
        call gij(f,iaa,aij,1)
        call curl_mn(aij,bb,aa)
!
        call dot2_mn(bb,b2)
!
        f(l1:l2,m,n,ilnrho) = alog(b2/(mu0*const_alfven**2.))
!
      enddo
    enddo
!
  endsubroutine const_alfven_speed
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
