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
  use Messages
  use Sub, only: keep_compiler_quiet
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
  real :: rho_init=0.
  real :: T0=6000.,T1=1e6,z0_tanh=4e6,width_tanh=1e6
  character (len=labellen) :: direction='z'
  real, dimension(4) :: mpoly_special = (/1.3,1000.,-1.04,500./)
  real, dimension(3) :: zpoly = (/0.,3.,5./)
!
  namelist /initial_condition_pars/ &
      lnrho_init,lnTT_init,stratitype,rho_init,direction, &
      set_lnTT_first,T0,T1,z0_tanh,width_tanh,mpoly_special,zpoly
!
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
    real, dimension (mx,my,mz,mfarray) :: f
!
    if (iproc==0) then
      write(*,*) "------------------------------------------------------------------"
      write(*,*) "Parameters to be set in run.in:"
      write(*,*)
      write(*,*) "Kpara=",2e-11 /unit_density/unit_velocity**3./ &
          unit_length*unit_temperature**3.5
      write(*,*) "Kperp=",3.47e12 / ((unit_velocity**3*unit_magnetic**2*unit_length)/ &
          (unit_density * sqrt(unit_temperature) ))
      write(*,*) "hcond_grad=",unit_temperature
      write(*,*) "------------------------------------------------------------------"
    endif
!
    call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
  subroutine read_initial_condition_pars(unit,iostat)
!
!  04-sep-10/bing: coded
!
    integer, intent(in) :: unit
    integer, intent(inout), optional :: iostat
!
    if (present(iostat)) then
      read(unit,NML=initial_condition_pars,ERR=99, IOSTAT=iostat)
    else
      read(unit,NML=initial_condition_pars,ERR=99)
    endif
!
 99  return
!
  endsubroutine read_initial_condition_pars
!***********************************************************************
  subroutine write_initial_condition_pars(unit)
!
!  04-sep-10/bing: coded
!
    integer, intent(in) :: unit
!
    write(unit,NML=initial_condition_pars)
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
    case default
      call fatal_error('initial_condition_lnrho', &
          'no such value for lnTT_init')
    endselect
!
!    if (stratitype=='nothing') then
      !   call warning('initial_condition_lnrho','Nothing to do')
    ! elseif (stratitype=='hydrostatic') then
    !   if (direction=='z') then
    !     call hydrostatic_z(f)
    !   elseif (direction=='x') then
    !     call hydrostatic_x(f)
    !   endif
    ! else (stratitype=='hydrostatic_lnTT')

    ! else
    !   call setup_vert_profiles(f)
    ! endif
!
! save to file stratification.dat
!
    call write_stratification_dat(f)
!
  endsubroutine initial_condition_lnrho
!***********************************************************************
  subroutine setup_vert_profiles(f)
!
!  Read and set vertical profiles for initial temperature and density.
!  Initial temperature profile is given in ln(T) [K] over z [Mm]
!  Initial density profile is given in ln(rho) [kg/m^3] over z [Mm]
!
!  04-sep-10/bing: coded
!
      use Mpicomm, only: mpibcast_int, mpibcast_real, stop_it_if_any
      use Messages, only: warning
      use Syscalls, only: file_exists, file_size
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real :: dummy
      integer :: lend,ierr
      integer :: i,j
      integer, parameter :: unit=12
!
! file location settings
      character (len=*), parameter :: lnrho_dat = 'prof_lnrho.dat'
      character (len=*), parameter :: lnT_dat = 'prof_lnT.dat'
!
      integer :: prof_nz
      real, dimension (:), allocatable :: prof_z, prof_lnrho, prof_lnTT
      logical :: lread_lnrho=.false., lread_lnTT=.false.
!
      inquire(IOLENGTH=lend) dummy
!
      lread_lnTT=(lnTT_init=='prof_lnTT')
      lread_lnrho=(lnrho_init=='prof_lnrho')
!
! read temperature profile for interpolation
      if (lread_lnTT.and.ltemperature.and..not.ltemperature_nolog) then
!
! file access is only done on the MPI root rank
        if (lroot) then
          if (.not. file_exists (lnT_dat)) call stop_it_if_any ( &
              .true., 'setup_special: file not found: '//trim(lnT_dat))
! find out, how many data points our profile file has
          prof_nz = (file_size (lnT_dat) - 2*2*4) / (lend*4 * 2)
        endif
        call stop_it_if_any(.false.,'')
        call mpibcast_int(prof_nz,1)
!
        allocate (prof_z(prof_nz), prof_lnTT(prof_nz), stat=ierr)
!
        if (lroot) then
          open (unit,file=lnT_dat,form='unformatted',status='unknown',recl=lend*prof_nz)
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
        prof_lnTT = prof_lnTT - alog(real(unit_temperature))
!
! convert z coordinates from [Mm] to Pencil units
        prof_z = prof_z / unit_length
!
! interpolate temperature profile to Pencil grid
!
        if (direction=='z') then
          do j = n1-nghost, n2+nghost
            if (z(j) < prof_z(1) ) then
              !call warning("setup_special","extrapolated lnT below bottom of initial profile")
              f(:,:,j,ilnTT) = prof_lnTT(1)
            elseif (z(j) >= prof_z(prof_nz)) then
              !call warning("setup_special","extrapolated lnT over top of initial profile")
              f(:,:,j,ilnTT) = prof_lnTT(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((z(j) >= prof_z(i)) .and. (z(j) < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  f(:,:,j,ilnTT) = (prof_lnTT(i+1)-prof_lnTT(i)) / &
                      (prof_z(i+1)-prof_z(i)) * (z(j)-prof_z(i)) + prof_lnTT(i)
                  exit
                endif
              enddo
            endif
          enddo
        else if (direction=='x') then
          do j = l1-nghost, l2+nghost
            if (x(j) < prof_z(1) ) then
              !call warning("setup_special","extrapolated lnT below bottom of initial profile")
              f(j,:,:,ilnTT) = prof_lnTT(1)
            elseif (x(j) >= prof_z(prof_nz)) then
              !call warning("setup_special","extrapolated lnT over top of initial profile")
              f(j,:,:,ilnTT) = prof_lnTT(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((x(j) >= prof_z(i)) .and. (x(j) < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  f(j,:,:,ilnTT) = (prof_lnTT(i+1)-prof_lnTT(i)) / &
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
      endif
!
! read density profile for interpolation
      if (lread_lnrho) then
!
! file access is only done on the MPI root rank
        if (lroot) then
          if (.not. file_exists (lnrho_dat)) call stop_it_if_any ( &
              .true., 'setup_special: file not found: '//trim(lnrho_dat))
! find out, how many data points our profile file has
          prof_nz = (file_size (lnrho_dat) - 2*2*4) / (lend*4 * 2)
        endif
        call stop_it_if_any(.false.,'')
        call mpibcast_int (prof_nz,1)
!
        allocate (prof_z(prof_nz), prof_lnrho(prof_nz), stat=ierr)
!
        if (lroot) then
          open (unit,file=lnrho_dat,form='unformatted',status='unknown',recl=lend*prof_nz)
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
        prof_lnrho = prof_lnrho - alog(real(unit_density))
!
! convert z coordinates from [Mm] to Pencil units
        prof_z = prof_z / unit_length
!
! interpolate density profile to Pencil grid
        if (direction=='z') then
          do j = n1-nghost, n2+nghost
            if (z(j) < prof_z(1) ) then
              !call warning("setup_special","extrapolated density below bottom of initial profile")
              f(:,:,j,ilnrho) = prof_lnrho(1)
            elseif (z(j) >= prof_z(prof_nz)) then
              !call warning("setup_special","extrapolated density over top of initial profile")
              f(:,:,j,ilnrho) = prof_lnrho(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((z(j) >= prof_z(i)) .and. (z(j) < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  f(:,:,j,ilnrho) = (prof_lnrho(i+1)-prof_lnrho(i)) / &
                      (prof_z(i+1)-prof_z(i)) * (z(j)-prof_z(i)) + prof_lnrho(i)
                  exit
                endif
              enddo
            endif
          enddo
        elseif (direction=='x') then
          do j = l1-nghost, l2+nghost
            if (x(j) < prof_z(1) ) then
              !call warning("setup_special","extrapolated density below bottom of initial profile")
              f(j,:,:,ilnrho) = prof_lnrho(1)
            elseif (x(j) >= prof_z(prof_nz)) then
              !call warning("setup_special","extrapolated density over top of initial profile")
              f(j,:,:,ilnrho) = prof_lnrho(prof_nz)
            else
              do i = 1, prof_nz-1
                if ((x(j) >= prof_z(i)) .and. (x(j) < prof_z(i+1))) then
                  ! linear interpolation: y = m*(x-x1) + y1
                  f(j,:,:,ilnrho) = (prof_lnrho(i+1)-prof_lnrho(i)) / &
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
    use EquationOfState, only: gamma,cs2top,cs2bot
    use Gravity, only: gravz
    use Mpicomm, only: mpibcast_real,mpibcast_int,stop_it_if_any
    use Syscalls, only: file_exists, file_size
!
    real, dimension (mx,my,mz,mfarray) :: f
    real :: ztop,zbot
    integer :: prof_nz
    real, dimension(:), allocatable :: prof_lnTT,prof_z
    real :: tmp_lnrho,tmp_lnT,tmpdT,tmp_z,dz_step,lnrho_0
    integer :: i,lend,j,ierr,unit=1
!
! file location settings
    character (len=*), parameter :: lnT_dat = 'prof_lnT.dat'
!
    inquire(IOLENGTH=lend) lnrho_0

    if (lentropy.or.ltemperature_nolog) &
        call fatal_error('hydrostatic','only implemented for ltemperature')
!
    lnrho_0 = alog(rho_init)
!
! read temperature profile for interpolation
!
! file access is only done on the MPI root rank
    if (lroot) then
      if (.not. file_exists (lnT_dat)) call stop_it_if_any ( &
          .true., 'setup_special: file not found: '//trim(lnT_dat))
! find out, how many data points our profile file has
      prof_nz = (file_size (lnT_dat) - 2*2*4) / (lend*4 * 2)
    endif
!
    call stop_it_if_any(.false.,'')
    call mpibcast_int (prof_nz,1)
!
    allocate (prof_z(prof_nz), prof_lnTT(prof_nz), stat=ierr)
    if (ierr > 0) call stop_it_if_any (.true.,'setup_special: '// &
        'Could not allocate memory for z coordinate or lnTT profile')
!
    if (lroot) then
      open (unit,file=lnT_dat,form='unformatted',status='unknown',recl=lend*prof_nz)
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
    prof_lnTT = prof_lnTT - alog(real(unit_temperature))
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
            tmpdT = (prof_lnTT(i+1)-prof_lnTT(i))/(prof_z(i+1)-prof_z(i)) * (tmp_z-prof_z(i)) + prof_lnTT(i) -tmp_lnT
            tmp_lnT = tmp_lnT + tmpdT
          elseif (tmp_z >= prof_z(prof_nz)) then
            tmpdT = prof_lnTT(prof_nz) - tmp_lnT
            tmp_lnT = tmp_lnT + tmpdT
          endif
        enddo
        tmp_lnrho = tmp_lnrho - tmpdT + gamma/(gamma-1.)*gravz*exp(-tmp_lnT) * dz_step
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
    use EquationOfState, only: gamma,gamma_m1,get_cp1,lnrho0
    use Gravity, only: gravz
    use Mpicomm, only: mpisend_real,mpirecv_real
!
    real, dimension (mx,my,mz,mfarray) :: f
    real :: konst,cp1=1.,lnrho_0,int
    real, dimension(mz) :: TT,lnTT
    integer :: i,j,k,ipt
!
    if (leos) call get_cp1(cp1)
!
    konst = gamma*cp1/gamma_m1
!
    if (lroot) then
      lnrho_0=alog(rho_init)
    else
      lnrho_0=0.
    endif
    f(:,:,n1,ilnrho)= lnrho_0
!
    if (ltemperature_nolog) then
      TT = f(l1,m1,:,iTT)
      lnTT = alog(f(l1,m1,:,iTT))
    else
      TT = exp(f(l1,m1,:,ilnTT))
      lnTT = f(l1,m1,:,ilnTT)
    endif
!
    do i=n1+1,n2+nghost
      int = 0.5 * (z(i)-z(i-1)) * (gravz*(1/TT(i)+1/TT(i+1)))
      f(:,:,i,ilnrho)=f(:,:,i-1,ilnrho)-lnTT(i)+ &
          lnTT(i-1)+konst*int
    enddo
!
    do i=0,nprocz-2
      ipt = nprocxy*i
      if (iproc==ipt) then 
        do j=0,nprocx-1
          do k=0,nprocy-1
            ipt=j + nprocx*k+nprocxy*(ipz+1)
            call mpisend_real(f(l1,m1,n2+1,ilnrho),1,ipt,100*ipz+j*k)
          enddo
        enddo
      endif
    enddo
!
    if (ipz>0) then 
      call mpirecv_real(lnrho_0,1,nprocxy*(ipz-1),100*(ipz-1)+ipx*ipy)
      f(:,:,:,ilnrho) =  f(:,:,:,ilnrho)+lnrho_0
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
    use EquationOfState, only: gamma,cs2top,cs2bot,get_cp1
    use Gravity, only: get_xgravity
    use Mpicomm, only: mpibcast_real,mpibcast_int,stop_it_if_any
    use Syscalls, only: file_exists, file_size
!
    real, dimension (mx,my,mz,mfarray) :: f
    integer :: prof_nx
    real, dimension(:), allocatable :: prof_lnTT,prof_x
    real :: tmp_lnrho,lnrho_0,ztmp,integrand
    integer :: i,lend,j,ierr,unit=1
    real, dimension(mx) :: xgrav,lnTT_loop
    real :: cp1=1
!
! file location settings
    character (len=*), parameter :: lnT_dat = 'prof_lnT.dat'
!
    if (leos) call get_cp1(cp1)
    call get_xgravity(xgrav)
!
    inquire(IOLENGTH=lend) lnrho_0

    if (lentropy.or.ltemperature_nolog) &
        call fatal_error('hydrostatic','only implemented for ltemperature')
!
    lnrho_0 = alog(rho_init)
!
! read temperature profile for interpolation
!
! file access is only done on the MPI root rank
    if (lroot) then
      if (.not. file_exists (lnT_dat)) call stop_it_if_any ( &
          .true., 'setup_special: file not found: '//trim(lnT_dat))
! find out, how many data points our profile file has
      prof_nx = (file_size (lnT_dat) - 2*2*4) / (lend*4 * 2)
    endif
!
    call stop_it_if_any(.false.,'')
    call mpibcast_int (prof_nx,1)
!
    allocate (prof_x(prof_nx), prof_lnTT(prof_nx), stat=ierr)
    if (ierr > 0) call stop_it_if_any (.true.,'setup_special: '// &
        'Could not allocate memory for x coordinate or lnTT profile')
!
    if (lroot) then
      open (unit,file=lnT_dat,form='unformatted',status='unknown',recl=lend*prof_nx)
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
    prof_x = prof_x/unit_length*1e6
    prof_lnTT = prof_lnTT - alog(real(unit_temperature))
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
          lnTT_loop(i) = lin_inpol(prof_lnTT(j+1),prof_lnTT(j),prof_x(j+1),prof_x(j),ztmp)
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

    real :: f2,f1,x2,x1,yn
    real :: lin_inpol
!
    lin_inpol = (f2-f1)/(x2-x1)*(yn-x1) + f1
!
  endfunction lin_inpol
!***********************************************************************\
  subroutine setup_tanh(f)
!    
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
    real, dimension (mz) :: TT,lnTT,z_SI
    integer :: i,j
!
    z_SI = z*unit_length
!
    TT = (T1-T0)*(0.5*tanh((z_SI-z0_tanh)/width_tanh)+0.5)+T0
    TT = TT / unit_temperature
    
    lnTT = alog(TT)
!
    do i=1,mx
      do j=1,my
        f(i,j,:,ilnTT) = lnTT
      enddo
    enddo
!
  endsubroutine setup_tanh
!***********************************************************************
  subroutine piecewice_poly(f)
!    
    use EquationOfState, only: gamma, gamma_m1, get_cp1
    use Gravity, only: gravz
!
    real, dimension(mx,my,mz,mfarray), intent(inout) :: f
    real :: Ttop,T2, T1, T0, cp1=1, temp
    real :: lnrhotop, lnrho2, lnrho1, lnrho0, ztop
    real :: lnrhobot,zbot,Tbot
    real, dimension(4) :: beta
    integer :: i
!
    if (leos) call get_cp1(cp1)
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
      lnrho0 = alog(3e-4/real(unit_density))
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
